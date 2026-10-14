# -*- coding: utf-8 -*-
"""
backend_rep_period.py — Representative period identification.

Preserves the v1 algorithm from rep_period.py (780 lines).

Algorithm summary:
    1. Read wind time series CSV (date, speed, direction)
    2. Classify wind into 16 directional sectors × 12 Beaufort speed classes
    3. For user-specified search windows (e.g., 3-month, 6-month):
       - Slide window across the analysis period
       - For each position, compute:
         a. Speed-parametrised directional similarity (weighted correlation)
         b. Direction-parametrised speed similarity (weighted correlation)
         c. rep_score = sum of both (lower = more representative)
    4. Quality checks: R², Pearson correlation, MAE
    5. Sort by rep_score, output ranked candidates to CSV

v1 source: rep_period.py (780 lines)
v2 improvements:
    - classify_wind_data and classify_wind_data_all_params merged into one
    - Shadowed builtins fixed (v1 used `list` as variable name)
    - Uses logger instead of print()
    - Progress callback for GUI
    - Constants extracted (Beaufort bins, directional sectors)
    - Search windows are exact calendar offsets from each window start
      (so '1MS' is one month from any start date) and stop at the scanning
      end instead of running one window length past it

Usage:
    from backend_rep_period import identify_rep_period
"""

import re
import time
import logging
import threading
import warnings
import numpy as np
import pandas as pd
from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error
from sklearn import preprocessing
from scipy.stats import pearsonr, ConstantInputWarning
from pathlib import Path
from typing import Callable, List, Optional, Tuple

# Scanning windows that fall in data gaps produce constant histograms, which
# makes Pearson correlation undefined (returns nan). That is expected — the
# ranking falls back to the other QC metrics. Silence the flood of warnings.
warnings.filterwarnings("ignore", category=ConstantInputWarning)


class TaskCancelledError(Exception):
    """Raised when the backend observes a cancel_event set by the GUI."""
    pass


def _check_cancel(cancel_event: Optional[threading.Event]) -> None:
    if cancel_event is not None and cancel_event.is_set():
        raise TaskCancelledError("Cancelled by user")

# ─── CONSTANTS ──────────────────────────────────────────────────────────────
DIRECTIONAL_SECTOR_LABELS = np.array(
    'N NNE NE ENE E ESE SE SSE S SSW SW WSW W WNW NW NNW N'.split()
)
DIRECTIONAL_BINS = np.arange(11.25, 360, 22.5)

ALL_16_SECTORS = [
    'N', 'NNE', 'NE', 'ENE', 'E', 'ESE', 'SE', 'SSE',
    'S', 'SSW', 'SW', 'WSW', 'W', 'WNW', 'NW', 'NNW',
]
ALL_12_SPEED_CLASSES = list(range(12))

BEAUFORT_BINS = np.array([0.5, 1.5, 3.3, 5.5, 7.9, 10.7, 13.8, 17.1, 20.7, 24.4, 28.4])
BEAUFORT_CLASSES = np.array(list(range(12)))


# ─────────────────────────────────────────────────────────────────────────────
# MAIN FUNCTION
# ─────────────────────────────────────────────────────────────────────────────

def identify_rep_period(
    file_input: str,
    output_name: str,
    quad: List[str],
    spd: List[int],
    start_time_total: str,
    end_time_total: str,
    frequency: List[str],
    start_time: str,
    end_time: str,
    out_path: str,
    logger: Optional[logging.Logger] = None,
    progress_callback: Callable[[float], None] = lambda x: None,
    cancel_event: Optional[threading.Event] = None,
) -> dict:
    """
    Identify representative periods from wind time series.

    Args:
        file_input:       Path to wind CSV file (columns: date, speed, dir).
        output_name:      Name prefix for output files.
        quad:             List of directional sectors to analyse
                          e.g., ['SSE', 'S', 'SSW', 'SW', 'WSW', 'W']
        spd:              List of Beaufort speed classes to analyse
                          e.g., [3, 4, 5, 6, 7]
        start_time_total: Start of the total reference period 'YYYY-MM-DD HH:MM:SS'
        end_time_total:   End of the total reference period
        frequency:        List of search window sizes e.g., ['3MS', '6MS', '12MS']
        start_time:       Start of the scanning window
        end_time:         End of the scanning window
        out_path:         Output directory path.
        logger:           Logger instance.
        progress_callback: Progress callback (0.0–1.0).

    Returns:
        {'output_file': str, 'n_candidates': int}

    v1 source: rep_period.py (entire file).
    """
    if logger is None:
        logger = logging.getLogger('repd3d.rep_period')

    # ── Build the search windows first, so bad input fails before the slow
    #    scan starts and every window size is reported loudly ─────────────
    windows = []
    for time_freq in frequency:
        val_freq, _ = parse_window_freq(time_freq)      # ValueError on bad input
        first_date, second_date = time_window_list(start_time, end_time, time_freq)
        if not first_date:
            raise ValueError(
                f"Scanning period {start_time} to {end_time} is shorter than one "
                f"'{time_freq}' search window. Extend the scanning period or choose "
                f"a smaller window."
            )
        windows.append((time_freq, val_freq, first_date, second_date))
        logger.info(
            f"Search window {time_freq}: {len(first_date)} positions, "
            f"first {first_date[0]} to {second_date[0]}, "
            f"last {first_date[-1]} to {second_date[-1]}"
        )

    logger.info("Reading wind data...")

    # ── Read and clean wind data ─────────────────────────────────────────
    wind = pd.read_csv(file_input, delimiter=',', skiprows=1,
                       names=['date', 'speed', 'dir'])
    wind = wind[(wind['dir'] >= 0) & (wind['dir'] <= 360)]
    wind = wind[(wind['speed'] >= 0) & (wind['speed'] <= 100)]

    date = pd.to_datetime(wind['date'].astype(str), format='%Y%m%d%H')
    obs_data = pd.DataFrame({
        'datetime': date.values,
        'speed': wind['speed'].values,
        'dir': wind['dir'].values,
    }).set_index('datetime')

    logger.info(f"Wind data: {len(obs_data)} records, "
                f"{obs_data.index[0]} to {obs_data.index[-1]}")

    # ── Classify the total reference period ──────────────────────────────
    dir_speed_total, df_dir, df_speed = classify_wind_data(
        obs_data, start_time_total, end_time_total, quad, spd
    )
    dir_class_df = speed_class_of_direction(dir_speed_total, quad)
    speed_class_df = direction_sector_of_speed(dir_speed_total, spd, quad)

    progress_callback(0.05)

    # ── File naming ──────────────────────────────────────────────────────
    name_span = f'{start_time[:4]}-{int(end_time[:4]) + 1}'
    quad_range = f'{quad[0]}-{quad[-1]}'
    spd_range = f'{spd[0]}-{spd[-1]}'
    freq_label = f'{frequency[0]}-{frequency[-1]}'

    # ── Scan across all windows ──────────────────────────────────────────
    rows = []
    rows_qc = []
    total_freq = len(frequency)

    for freq_idx, (time_freq, val_freq, first_date, second_date) in enumerate(windows):
        _check_cancel(cancel_event)
        t_in = time.time()
        total_windows = len(first_date)
        logger.info(f"Scanning window {time_freq}: {total_windows} positions")

        for w_idx, (f_date, s_date) in enumerate(zip(first_date, second_date)):
            _check_cancel(cancel_event)
            # Classify this window
            ds_compare, df_dir_compare, df_speed_compare = classify_wind_data(
                obs_data, f_date, s_date, quad, spd
            )
            dc_compare = speed_class_of_direction(ds_compare, quad)
            sc_compare = direction_sector_of_speed(ds_compare, spd, quad)

            # Speed-parametrised directional similarity
            df_speed_corr = correlate(dir_class_df, dc_compare, quad)
            weight_total_d = df_dir.iloc[:, 1]
            weight_compare_d = df_dir_compare.iloc[:, 1]
            dir_weighted_avg = np.average(df_speed_corr.iloc[:, 1], weights=weight_total_d)
            euc_dist_d = np.sqrt(np.sum(np.square(weight_total_d - weight_compare_d)))
            speed_para = dir_weighted_avg * euc_dist_d

            # Direction-parametrised speed similarity
            df_dir_corr = correlate(speed_class_df, sc_compare, spd)
            weight_total_s = df_speed.iloc[:, 1]
            weight_compare_s = df_speed_compare.iloc[:, 1]
            speed_weighted_avg = np.average(df_dir_corr.iloc[:, 1], weights=weight_total_s)
            euc_dist_s = np.sqrt(np.sum(np.square(weight_total_s - weight_compare_s)))
            dir_para = speed_weighted_avg * euc_dist_s

            rows.append([f_date, s_date, val_freq, dir_para, speed_para])

            # Quality checks
            ds_qc, qc_dir, qc_speed = classify_wind_data(
                obs_data, start_time_total, end_time_total,
                ALL_16_SECTORS, ALL_12_SPEED_CLASSES,
            )
            ds_qc_c, qc_dir_c, qc_speed_c = classify_wind_data(
                obs_data, f_date, s_date,
                ALL_16_SECTORS, ALL_12_SPEED_CLASSES,
            )

            speed_pct = percent_sim(qc_speed, qc_speed_c)
            dir_pct = percent_sim(qc_dir, qc_dir_c)

            s_qc = qc_speed.iloc[:, 0]
            s_c = qc_speed_c.iloc[:, 0]
            r2_s = r2_score(s_qc, s_c)
            p_s = pearsonr(s_qc, s_c)

            d_qc = qc_dir.iloc[:, 0]
            d_c = qc_dir_c.iloc[:, 0]
            r2_d = r2_score(d_qc, d_c)
            p_d = pearsonr(d_qc, d_c)

            mae_dir = round(mean_absolute_error(d_c, d_qc), 3)
            mae_spd = round(mean_absolute_error(s_c, s_qc), 3)

            rows_qc.append([
                f_date, dir_pct, speed_pct,
                r2_d, p_d[0], p_d[1],
                r2_s, p_s[0], p_s[1],
                mae_dir, mae_spd,
            ])

            # Progress
            overall = (freq_idx + (w_idx + 1) / total_windows) / total_freq
            progress_callback(0.05 + 0.9 * overall)

        elapsed = (time.time() - t_in) / 60
        logger.info(f"  Completed {time_freq} in {elapsed:.2f} min")

    # ── Combine results ──────────────────────────────────────────────────
    df_results = pd.DataFrame(
        rows, columns=['start_point', 'end_point', 'period_freq', 'dir_para', 'speed_para']
    ).set_index('start_point')
    df_results['rep_score'] = df_results['dir_para'] + df_results['speed_para']

    df_qc = pd.DataFrame(
        rows_qc,
        columns=[
            'start_point', 'similarity_perc_direction_qc', 'similarity_perc_speed_qc',
            'dir_r2_qc', 'dir_corr_qc', 'p-val_dir_qc',
            'speed_r2_qc', 'speed_corr_qc', 'p-val_speed_qc',
            'mae_dir_qc', 'mae_speed_qc',
        ],
    ).set_index('start_point')

    rep_period = pd.concat([df_results, df_qc], axis=1)
    rep_period = rep_period.sort_values(by=[
        'rep_score', 'dir_r2_qc', 'similarity_perc_direction_qc',
        'speed_r2_qc', 'similarity_perc_speed_qc',
    ])

    # Write output
    out_file = f'{out_path}/rep_period_{output_name}_{name_span}_{quad_range}_{spd_range}_{freq_label}.txt'
    rep_period.to_csv(out_file)

    progress_callback(1.0)
    logger.info(f"Results written: {out_file} ({len(rep_period)} candidates)")

    return {
        'output_file': out_file,
        'n_candidates': len(rep_period),
    }


# ─────────────────────────────────────────────────────────────────────────────
# WIND CLASSIFICATION (merged from v1's two identical functions)
# ─────────────────────────────────────────────────────────────────────────────

def classify_wind_data(df_in, start, end, dir_sectors, speed_classes):
    """
    Classify wind data into directional sectors and speed classes.

    v1 source: rep_period.py lines 20-156 (classify_wind_data) and
    lines 158-288 (classify_wind_data_all_params — identical logic).
    v2 merges both into one function with sector/class parameters.

    Returns: (dir_speed_df, df_direction, df_speed)
    """
    req = df_in.loc[start:end]
    data_dir = req['dir']
    data_speed = req['speed']

    # Bin directions
    df_dir = pd.DataFrame(data_dir)
    df_dir['orientation'] = DIRECTIONAL_SECTOR_LABELS[np.digitize(df_dir['dir'], DIRECTIONAL_BINS)]

    # Bin speeds (Beaufort scale)
    df_spd = pd.DataFrame(data_speed)
    df_spd['class'] = BEAUFORT_CLASSES[np.digitize(df_spd['speed'], BEAUFORT_BINS)]

    # Filter to requested sectors and classes
    combined = pd.concat([df_dir, df_spd], axis=1)
    combined = combined[combined['orientation'].isin(dir_sectors)]
    combined = combined[combined['class'].isin(speed_classes)]

    # Direction frequencies
    dir_class = pd.DataFrame(combined.iloc[:, [0, 1]])
    df_direction = dir_class.groupby('orientation').count().reindex(dir_sectors).fillna(0)
    dir_arr = np.asarray(df_direction['dir'])
    df_direction['norm'] = np.transpose(preprocessing.normalize([dir_arr], norm='l2'))
    total_dir = df_direction['dir'].sum()
    df_direction['weightage'] = (df_direction['dir'] / total_dir * 100) if total_dir > 0 else 0
    df_direction = df_direction.iloc[:, [1, 2]]

    # Speed frequencies
    spd_class = pd.DataFrame(combined.iloc[:, [2, 3]])
    df_speed = spd_class.groupby('class').count()
    df_speed = df_speed.reset_index()
    existing = df_speed['class'].tolist()
    missing = [c for c in speed_classes if c not in existing]
    if missing:
        missing_df = pd.DataFrame({'class': missing, 'speed': [0] * len(missing)})
        df_speed = pd.concat([df_speed, missing_df]).sort_values('class').reset_index(drop=True)
    df_speed = df_speed.set_index('class')
    spd_arr = np.asarray(df_speed['speed'])
    df_speed['norm'] = np.transpose(preprocessing.normalize([spd_arr], norm='l2'))
    total_spd = df_speed['speed'].sum()
    df_speed['weightage'] = (df_speed['speed'] / total_spd * 100) if total_spd > 0 else 0
    df_speed = df_speed.iloc[:, [1, 2]]

    # Combined counts
    dir_speed = combined.groupby(['orientation', 'class']).count().reset_index()

    return dir_speed, df_direction, df_speed


# ─── Supporting functions preserved from v1 ─────────────────────────────────

def speed_class_of_direction(input_df, dir_sectors):
    """v1: rep_period.py lines 290-365."""
    df_list = []
    for q in dir_sectors:
        df = input_df[input_df['orientation'] == q].reset_index(drop=True)
        existing = df['class'].tolist()
        missing = [i for i in range(12) if i not in existing]
        if missing:
            rows = [[q, m, 0, 0] for m in missing]
            missing_df = pd.DataFrame(rows, columns=['orientation', 'class', 'dir', 'speed'])
            df = pd.concat([df, missing_df]).sort_values('class').reset_index(drop=True)
        total = df['dir'].sum()
        df['weightage'] = (df['dir'] / total * 100) if total > 0 else 0
        arr = np.asarray(df['dir'])
        norm = np.transpose(preprocessing.normalize([arr], norm='l2'))
        result = pd.DataFrame(norm, columns=['norm'])
        result['weight'] = df['weightage'].values
        result['orientation'] = df['orientation'].values
        df_list.append(result)
    return pd.concat(df_list)


def direction_sector_of_speed(input_df, speed_classes, dir_sectors):
    """v1: rep_period.py lines 367-407."""
    df_list = []
    for s in speed_classes:
        df = input_df[input_df['class'] == s].set_index('orientation').reindex(dir_sectors)
        df['dir'] = df['dir'].fillna(0)
        df['class'] = df['class'].fillna(s)
        total = df['dir'].sum()
        df['weightage'] = (df['dir'] / total * 100) if total > 0 else 0
        arr = np.asarray(df['dir'])
        norm = np.transpose(preprocessing.normalize([arr], norm='l2'))
        result = pd.DataFrame({'norm': norm.flatten()})
        result['weight'] = df['weightage'].values
        result['class'] = df['class'].values
        df_list.append(result)
    return pd.concat(df_list)


def correlate(input_df, compare_df, param_list):
    """v1: rep_period.py lines 409-465."""
    results = []
    for p in param_list:
        in_d = input_df[input_df.iloc[:, 2] == p]
        com_d = compare_df[compare_df.iloc[:, 2] == p]
        in_data = in_d.iloc[:, 0]
        com_data = com_d.iloc[:, 0]
        in_w = in_d.iloc[:, 1].fillna(0)
        com_w = com_d.iloc[:, 1].fillna(0)

        if com_data.sum() == 0:
            results.append([p, float('nan')])
        else:
            sq_diff = (in_data.values - com_data.values) ** 2
            wg_avg = np.average(sq_diff, weights=in_w.values) if in_w.sum() > 0 else 0
            euc = np.sqrt(np.sum(np.square(in_w.values - com_w.values)))
            results.append([p, wg_avg * euc])

    col_name = 'class' if isinstance(param_list[0], int) else 'orientation'
    df = pd.DataFrame(results, columns=[col_name, 'rep_metric'])
    max_val = df.iloc[:, 1].max()
    df = df.fillna((max_val + 1) if not np.isnan(max_val) else 1)
    return df


# Units accepted in a search-window string such as '1MS', '3MS', '2W', '45D'.
# Every month alias (MS, ME and the deprecated M) means the same thing here: a
# calendar-month shift from each window's own start date. The month-start /
# month-end anchoring these aliases carry in pandas is irrelevant for a
# sliding window and was the source of the v1 sizing bug described below.
_WINDOW_UNIT_TO_KWARG = {
    'MS': 'months', 'ME': 'months', 'M': 'months',
    'QS': 'quarters', 'QE': 'quarters', 'Q': 'quarters',
    'YS': 'years', 'YE': 'years', 'Y': 'years', 'AS': 'years', 'A': 'years',
    'W': 'weeks',
    'D': 'days',
    'H': 'hours',
}


def parse_window_freq(freq: str) -> Tuple[int, pd.DateOffset]:
    """
    Parse a search-window string into (count, calendar offset).

    '1MS' -> (1, DateOffset(months=1))     'MS'  -> (1, DateOffset(months=1))
    '3ME' -> (3, DateOffset(months=3))     '2W'  -> (2, DateOffset(weeks=2))
    '45D' -> (45, DateOffset(days=45))

    Raises ValueError with a readable message for anything else, so the GUI
    can reject bad input before a scan starts.
    """
    text = str(freq).strip()
    match = re.fullmatch(r'(\d*)\s*([A-Za-z]+)', text)
    unit = match.group(2).upper() if match else None
    if unit not in _WINDOW_UNIT_TO_KWARG:
        raise ValueError(
            f"Search window '{freq}' is not understood. Use a count followed by a "
            f"unit, e.g. '1MS' (one month), '3MS', '12MS', '2W' (weeks) or '45D' (days)."
        )
    count = int(match.group(1)) if match.group(1) else 1
    if count < 1:
        raise ValueError(f"Search window '{freq}' must be at least one unit long.")
    kwarg = _WINDOW_UNIT_TO_KWARG[unit]
    if kwarg == 'quarters':
        return count, pd.DateOffset(months=3 * count)
    return count, pd.DateOffset(**{kwarg: count})


def time_window_list(start, end, freq):
    """
    Build the sliding search windows for one window size.

    Window starts are every day from ``start`` to ``end``. Each window ends
    exactly one ``freq`` after its own start (one calendar month for '1MS',
    so 2015-01-15 -> 2015-02-15 and 2015-01-31 -> 2015-02-28). Windows that
    would end after ``end`` are dropped, so the last window finishes on or
    before the scanning end the user set.

    Returns two equal-length lists: window starts and window ends. Both are
    empty when the scanning period is shorter than one window.

    v1 (rep_period.py lines 467-493) took the window length from the second
    element of an anchored ``pd.date_range``. That silently added the rest
    of the start month to every window (a '1MS' scan starting on the 2nd of
    a month gave two-month windows, and '1ME' always gave two months), and
    the last windows ran one window length past the scanning end.
    """
    _, offset = parse_window_freq(freq)
    start_ts, end_ts = pd.Timestamp(start), pd.Timestamp(end)
    first, second = [], []
    for s in pd.date_range(start=start_ts, end=end_ts, freq='1D'):
        e = s + offset
        if e > end_ts:
            break                      # ends only grow from here on
        first.append(s)
        second.append(e)
    return first, second


def percent_sim(full_df, compare_df):
    """v1: rep_period.py lines 495-516."""
    abs_diff = abs(full_df.iloc[:, 1] - compare_df.iloc[:, 1])
    w = full_df.iloc[:, 1]
    err_pct = np.where(w > 0, (abs_diff * 100) / w, 0)
    sim_pct = np.nan_to_num(100 - err_pct)
    wavg = np.average(sim_pct, weights=w) if w.sum() > 0 else 0
    euc = np.sqrt(np.sum(np.square(full_df.iloc[:, 1] - compare_df.iloc[:, 1])))
    return wavg - euc


# ─────────────────────────────────────────────────────────────────────────────
# SELF-TEST
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    print("=== backend_rep_period.py self-test ===\n")
    all_pass = True

    def check(name, result, expected, tol=0.01):
        global all_pass
        if isinstance(result, float):
            ok = abs(result - expected) < tol
        else:
            ok = result == expected
        status = "PASS" if ok else "FAIL"
        if not ok:
            all_pass = False
            print(f"  {name}: [{status}]  got {result!r}")
        else:
            print(f"  {name}: [{status}]")

    # Generate synthetic wind data
    np.random.seed(42)
    dates = pd.date_range('2015-01-01', '2016-12-31', freq='h')
    speeds = np.random.rayleigh(5, len(dates))
    dirs = np.random.uniform(0, 360, len(dates))
    obs = pd.DataFrame({'speed': speeds, 'dir': dirs}, index=dates)

    print("classify_wind_data:")
    sectors = ['N', 'S', 'E', 'W']
    classes = [3, 4, 5]
    ds, df_dir, df_speed = classify_wind_data(obs, '2015-01-01', '2015-12-31', sectors, classes)
    check("dir df has 4 rows", len(df_dir), 4)
    check("speed df has 3 rows", len(df_speed), 3)
    check("weights sum ~100", abs(df_dir['weightage'].sum() - 100) < 1, True)

    print("\nspeed_class_of_direction:")
    scd = speed_class_of_direction(ds, sectors)
    check("scd has orientation col", 'orientation' in scd.columns, True)

    print("\ndirection_sector_of_speed:")
    dss = direction_sector_of_speed(ds, classes, sectors)
    check("dss has class col", 'class' in dss.columns, True)

    print("\nparse_window_freq:")
    check("'1MS' -> 1 month", parse_window_freq('1MS'), (1, pd.DateOffset(months=1)))
    check("'MS' defaults to 1 month", parse_window_freq('MS'), (1, pd.DateOffset(months=1)))
    check("'3ME' same as '3MS'", parse_window_freq('3ME'), parse_window_freq('3MS'))
    check("'12MS' -> 12 months", parse_window_freq('12MS'), (12, pd.DateOffset(months=12)))
    check("'2W' -> 2 weeks", parse_window_freq('2W'), (2, pd.DateOffset(weeks=2)))
    check("'45D' -> 45 days", parse_window_freq('45D'), (45, pd.DateOffset(days=45)))
    for bad in ('', '3', '1.5MS', '0MS', 'monthly'):
        try:
            parse_window_freq(bad)
            check(f"rejects {bad!r}", "no error", "ValueError")
        except ValueError:
            check(f"rejects {bad!r}", True, True)

    print("\ntime_window_list:")
    f, s = time_window_list('2015-01-01', '2016-01-01', '3MS')
    check("window list not empty", len(f) > 0, True)
    check("first and second same length", len(f), len(s))
    check("first 3MS window is 3 calendar months", s[0], pd.Timestamp('2015-04-01'))
    check("windows slide by one day", f[1] - f[0], pd.Timedelta(days=1))
    check("last window ends on the scan end", s[-1], pd.Timestamp('2016-01-01'))
    check("no window runs past the scan end", max(s) <= pd.Timestamp('2016-01-01'), True)
    f, s = time_window_list('2015-01-02', '2015-06-01', '1MS')
    check("1MS from the 2nd is one month, not two", s[0], pd.Timestamp('2015-02-02'))
    check("every 1MS window is one calendar month",
          all(e == b + pd.DateOffset(months=1) for b, e in zip(f, s)), True)
    check("1MS and 1ME give identical windows",
          time_window_list('2015-01-02', '2015-06-01', '1ME'), (f, s))
    f, s = time_window_list('2015-01-31', '2015-04-01', '1MS')
    check("month-end start clips to Feb 28", s[0], pd.Timestamp('2015-02-28'))
    check("too-short scan gives no windows",
          time_window_list('2015-01-01', '2015-01-20', '1MS'), ([], []))

    print("\npercent_sim (self-comparison):")
    sim = percent_sim(df_dir, df_dir)
    check("self-similarity high", sim > 50, True)

    print(f"\n{'=' * 50}")
    if all_pass:
        print("ALL TESTS PASSED")
    else:
        print("SOME TESTS FAILED")
    print(f"{'=' * 50}")
