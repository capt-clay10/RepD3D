# -*- coding: utf-8 -*-
"""
gui_windrose.py — Windrose comparison popup for the Rep Period module.

Three horizontal polar subplots:
    [ Reference period ] [ Prospect period ] [ Overlay ]

The overlay shows the prospect period as stacked bars with the reference
period overlaid as a dashed radial outline. Prev/Next + slider navigate the
top-K prospect periods.

Y-axis maximum and the speed-class colour legend are fixed across all
subplots so visual comparison is apples-to-apples.
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Tuple, Optional

import numpy as np
import pandas as pd
import customtkinter as ctk

import matplotlib
matplotlib.use('TkAgg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib import cm

from gui_theme import (
    COLORS, make_button, make_label, apply_window_icon,
)


# ─── Constants (mirrors backend_rep_period.py) ──────────────────────────────
SECTOR_LABELS = [
    'N', 'NNE', 'NE', 'ENE', 'E', 'ESE', 'SE', 'SSE',
    'S', 'SSW', 'SW', 'WSW', 'W', 'WNW', 'NW', 'NNW',
]
N_SECTORS = 16
DIR_BIN_EDGES = np.arange(11.25, 360, 22.5)  # 15 inner edges → 16 sectors (N wraps)
BEAUFORT_BIN_EDGES = np.array(
    [0.5, 1.5, 3.3, 5.5, 7.9, 10.7, 13.8, 17.1, 20.7, 24.4, 28.4]
)
N_SPEED_CLASSES = 12

SPEED_CMAP = cm.get_cmap('plasma', N_SPEED_CLASSES)
REF_OUTLINE_COLOR = COLORS.get('btn_warning', '#F67D65')

# Matplotlib colour hooks mapped to the app theme
FIG_BG = COLORS.get('bg', '#212729')
AXES_BG = COLORS.get('bg_secondary', '#2A2F31')
TEXT_COLOR = COLORS.get('text', '#FEF7F2')
TEXT_DIM = COLORS.get('text_dim', '#A0A4A7')


# ─── Binning helpers ────────────────────────────────────────────────────────

def _bin_windrose(obs: pd.DataFrame, start: str, end: str) -> np.ndarray:
    """Return a 16×12 count matrix (sector × Beaufort class)."""
    sl = obs.loc[start:end]
    out = np.zeros((N_SECTORS, N_SPEED_CLASSES))
    if len(sl) == 0:
        return out
    d = np.digitize(sl['dir'].values, DIR_BIN_EDGES) % N_SECTORS
    s = np.clip(np.digitize(sl['speed'].values, BEAUFORT_BIN_EDGES),
                0, N_SPEED_CLASSES - 1)
    for di, si in zip(d, s):
        out[di, si] += 1
    return out


def _counts_to_percent(counts: np.ndarray) -> np.ndarray:
    total = counts.sum()
    return counts / total * 100.0 if total > 0 else counts


def _load_obs(wind_csv_path: str) -> pd.DataFrame:
    """Load the wind CSV the same way backend_rep_period.py does."""
    wind = pd.read_csv(wind_csv_path, delimiter=',', skiprows=1,
                       names=['date', 'speed', 'dir'])
    wind = wind[(wind['dir'] >= 0) & (wind['dir'] <= 360)]
    wind = wind[(wind['speed'] >= 0) & (wind['speed'] <= 100)]
    date = pd.to_datetime(wind['date'].astype(str), format='%Y%m%d%H')
    return pd.DataFrame({
        'datetime': date.values,
        'speed': wind['speed'].values,
        'dir': wind['dir'].values,
    }).set_index('datetime')


# ─── Plotting ───────────────────────────────────────────────────────────────

def _style_polar_ax(ax, title: str):
    ax.set_facecolor(AXES_BG)
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.set_xticks(np.deg2rad(np.arange(0, 360, 45)))
    ax.set_xticklabels(['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW'],
                       color=TEXT_COLOR, fontsize=8)
    ax.tick_params(colors=TEXT_DIM, labelsize=7)
    for spine in ax.spines.values():
        spine.set_color(TEXT_DIM)
    ax.grid(color=TEXT_DIM, alpha=0.25, linewidth=0.5)
    ax.set_title(title, color=TEXT_COLOR, fontsize=10, pad=10)


def _draw_rose_bars(ax, pct: np.ndarray):
    """Draw stacked bars (sector × speed class) on a polar axes."""
    theta = np.deg2rad(np.arange(0, 360, 360 / N_SECTORS))
    width = np.deg2rad(360 / N_SECTORS) * 0.9
    bottom = np.zeros(N_SECTORS)
    for sc in range(N_SPEED_CLASSES):
        heights = pct[:, sc]
        if heights.sum() == 0:
            bottom += heights
            continue
        ax.bar(theta, heights, width=width, bottom=bottom,
               color=SPEED_CMAP(sc), edgecolor=AXES_BG, linewidth=0.4,
               zorder=2)
        bottom += heights


def _draw_rose_outline(ax, pct: np.ndarray, color: str):
    """Draw the per-sector total as a closed dashed polygon (reference contour)."""
    totals = pct.sum(axis=1)
    theta = np.deg2rad(np.arange(0, 360, 360 / N_SECTORS))
    t = np.append(theta, theta[0])
    r = np.append(totals, totals[0])
    ax.plot(t, r, color=color, linewidth=2.0, linestyle='--', zorder=3,
            label='Reference')


# ─── Popup window ───────────────────────────────────────────────────────────

class WindroseCompareWindow(ctk.CTkToplevel):
    """Popup that compares reference vs. prospect wind roses."""

    def __init__(
        self,
        parent,
        wind_csv_path: str,
        ref_start: str,
        ref_end: str,
        prospect_periods: List[Tuple[str, str, float]],
        title: str = "Wind Rose Comparison",
    ):
        super().__init__(parent)
        self.title(title)
        self.geometry("1280x620")
        self.configure(fg_color=FIG_BG)
        apply_window_icon(self)

        self._prospects = prospect_periods
        self._idx = 0

        # Load + pre-compute all data up front so navigation is instant
        obs = _load_obs(wind_csv_path)
        self._ref_pct = _counts_to_percent(_bin_windrose(obs, ref_start, ref_end))
        self._prospect_pcts = [
            _counts_to_percent(_bin_windrose(obs, s, e))
            for (s, e, _) in prospect_periods
        ]

        # Fixed y-axis max across every subplot (reference + all prospects)
        all_sector_totals = [self._ref_pct.sum(axis=1).max()] + [
            p.sum(axis=1).max() for p in self._prospect_pcts
        ]
        self._y_max = float(max(all_sector_totals)) * 1.1 or 10.0

        self._ref_label = f"Reference\n{ref_start}  →  {ref_end}"

        self._build_ui()
        self._render()

    def _build_ui(self):
        # Matplotlib figure — 1 row × 3 polar subplots + fixed colorbar on right
        self.fig = Figure(figsize=(12.5, 5.2), facecolor=FIG_BG)
        gs = self.fig.add_gridspec(
            1, 4, width_ratios=[1, 1, 1, 0.05], wspace=0.35,
        )
        self.ax_ref = self.fig.add_subplot(gs[0, 0], projection='polar')
        self.ax_pro = self.fig.add_subplot(gs[0, 1], projection='polar')
        self.ax_ovl = self.fig.add_subplot(gs[0, 2], projection='polar')
        self.cax = self.fig.add_subplot(gs[0, 3])

        # Fixed speed-class colour legend
        import matplotlib.colors as mcolors
        bounds = np.arange(N_SPEED_CLASSES + 1) - 0.5
        norm = mcolors.BoundaryNorm(bounds, SPEED_CMAP.N)
        cb = self.fig.colorbar(
            cm.ScalarMappable(norm=norm, cmap=SPEED_CMAP),
            cax=self.cax, ticks=np.arange(N_SPEED_CLASSES),
        )
        cb.set_label('Beaufort class', color=TEXT_COLOR, fontsize=9)
        cb.ax.yaxis.set_tick_params(color=TEXT_DIM, labelsize=8)
        for t in cb.ax.get_yticklabels():
            t.set_color(TEXT_COLOR)
        cb.outline.set_edgecolor(TEXT_DIM)

        self.canvas = FigureCanvasTkAgg(self.fig, master=self)
        self.canvas.get_tk_widget().pack(fill="both", expand=True,
                                          padx=10, pady=(10, 4))

        # Navigation row
        nav = ctk.CTkFrame(self, fg_color="transparent")
        nav.pack(fill="x", padx=10, pady=(4, 10))

        self.prev_btn = make_button(
            nav, "◀ Prev", command=self._on_prev,
            style="secondary", width=100,
        )
        self.prev_btn.pack(side="left", padx=(0, 10))

        slider_frame = ctk.CTkFrame(nav, fg_color="transparent")
        slider_frame.pack(side="left", fill="x", expand=True)
        self.slider_label = make_label(
            slider_frame, "", font_key="body_sm", color=COLORS["text"],
        )
        self.slider_label.pack(anchor="w")
        if len(self._prospects) > 1:
            self.slider = ctk.CTkSlider(
                slider_frame,
                from_=0, to=len(self._prospects) - 1,
                number_of_steps=len(self._prospects) - 1,
                command=self._on_slider,
            )
            self.slider.set(0)
            self.slider.pack(fill="x", pady=(2, 0))
        else:
            self.slider = None

        self.next_btn = make_button(
            nav, "Next ▶", command=self._on_next,
            style="secondary", width=100,
        )
        self.next_btn.pack(side="left", padx=(10, 0))

    def _on_prev(self):
        if self._idx > 0:
            self._idx -= 1
            if self.slider is not None:
                self.slider.set(self._idx)
            self._render()

    def _on_next(self):
        if self._idx < len(self._prospects) - 1:
            self._idx += 1
            if self.slider is not None:
                self.slider.set(self._idx)
            self._render()

    def _on_slider(self, value):
        new_idx = int(round(float(value)))
        if new_idx != self._idx:
            self._idx = new_idx
            self._render()

    def _render(self):
        for ax in (self.ax_ref, self.ax_pro, self.ax_ovl):
            ax.clear()

        # Reference
        _style_polar_ax(self.ax_ref, self._ref_label)
        self.ax_ref.set_ylim(0, self._y_max)
        _draw_rose_bars(self.ax_ref, self._ref_pct)

        # Prospect
        pstart, pend, score = self._prospects[self._idx]
        p_pct = self._prospect_pcts[self._idx]
        pro_title = (
            f"Prospect {self._idx + 1} / {len(self._prospects)}\n"
            f"{pstart}  →  {pend}\n"
            f"rep_score = {score:.4f}"
        )
        _style_polar_ax(self.ax_pro, pro_title)
        self.ax_pro.set_ylim(0, self._y_max)
        _draw_rose_bars(self.ax_pro, p_pct)

        # Overlay — prospect bars + reference outline contour
        _style_polar_ax(self.ax_ovl, "Overlay\n(prospect bars, reference = dashed outline)")
        self.ax_ovl.set_ylim(0, self._y_max)
        _draw_rose_bars(self.ax_ovl, p_pct)
        _draw_rose_outline(self.ax_ovl, self._ref_pct, REF_OUTLINE_COLOR)
        self.ax_ovl.legend(loc='upper right', bbox_to_anchor=(1.15, 1.1),
                            fontsize=8, facecolor=AXES_BG,
                            edgecolor=TEXT_DIM, labelcolor=TEXT_COLOR)

        self.slider_label.configure(
            text=f"Prospect {self._idx + 1} of {len(self._prospects)}"
        )
        self.canvas.draw_idle()


# ─── Helpers usable from gui_rep_period.py ──────────────────────────────────

def load_top_prospects(
    ranked_csv_path: str, k: int = 20,
) -> List[Tuple[str, str, float]]:
    """Parse the rep-period ranked CSV and return the top-K prospects."""
    df = pd.read_csv(ranked_csv_path, index_col=0)
    df = df.head(k)
    out = []
    for start_point, row in df.iterrows():
        out.append((
            str(start_point),
            str(row['end_point']),
            float(row.get('rep_score', 0.0)),
        ))
    return out


# ─── Self-test ──────────────────────────────────────────────────────────────
if __name__ == '__main__':
    # Minimal synthetic test — build a fake CSV and prospects, show window.
    import tempfile
    tmp = tempfile.mkdtemp()
    csv_path = Path(tmp) / 'wind.csv'
    rng = np.random.default_rng(0)
    n = 24 * 365 * 3
    dates = pd.date_range('2010-01-01', periods=n, freq='h')
    dir_vals = (rng.uniform(200, 280, n) + rng.normal(0, 30, n)) % 360
    spd_vals = np.clip(rng.gamma(3.5, 1.8, n), 0, 30)
    df = pd.DataFrame({
        'date': dates.strftime('%Y%m%d%H'),
        'speed': np.round(spd_vals, 2),
        'dir': np.round(dir_vals, 2),
    })
    with open(csv_path, 'w') as f:
        f.write('header line\n')
        df.to_csv(f, index=False, header=False)

    prospects = [
        ('2011-01-01 00:00:00', '2011-06-30 00:00:00', 0.12),
        ('2011-07-01 00:00:00', '2011-12-31 00:00:00', 0.18),
        ('2012-01-01 00:00:00', '2012-06-30 00:00:00', 0.22),
    ]

    from gui_theme import apply_theme
    apply_theme()
    root = ctk.CTk()
    root.withdraw()
    win = WindroseCompareWindow(
        root, str(csv_path),
        '2010-01-01 00:00:00', '2012-12-31 23:00:00',
        prospects,
    )
    win.protocol("WM_DELETE_WINDOW", root.destroy)
    root.mainloop()
