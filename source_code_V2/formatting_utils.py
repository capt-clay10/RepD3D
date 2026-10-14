# -*- coding: utf-8 -*-
"""
formatting_utils.py — Number formatting for Delft3D file output.

Consolidates identical formatting functions that were copy-pasted across
6 files in v1:
    - bct_generator.py (lines 59-87)
    - bct_year_overlap_file_generator.py (lines 70-98)
    - bcw_generator.py (lines 45-83, 231-236)
    - bcw_year_overlap_file_generator.py (same)
    - sea_level_change.py (lines 75-98)
    - mdw_writer.py (lines 20-35)

Every function preserves the EXACT output format of v1 for backward
compatibility — Delft3D is extremely sensitive to column alignment and
number formatting in its input files.

Usage:
    from formatting_utils import float_to_sci, list_to_sci, pad_positive
"""

import numpy as np
from typing import List, Union


# ─────────────────────────────────────────────────────────────────────────────
# SCIENTIFIC NOTATION — used by .bct, .bcc, .mdw, sea_level_change
# ─────────────────────────────────────────────────────────────────────────────

def float_to_sci(value: float, precision: int = 7, exp_digits: int = 3) -> str:
    """
    Convert a float to Delft3D scientific notation string.

    The format matches Delft3D-4 .bct / .bcc requirements:
        - Fixed number of digits after decimal in mantissa
        - Fixed number of digits in exponent (zero-padded)
        - Explicit sign on exponent

    Args:
        value:      The float to convert.
        precision:  Digits after decimal point in mantissa (default 7).
        exp_digits: Minimum digits in the exponent field (default 3).
                    The total exponent width is exp_digits + 1 (for the sign).

    Returns:
        Formatted string.

    Examples:
        >>> float_to_sci(20.0)
        '2.0000000e+001'
        >>> float_to_sci(-0.00345)
        '-3.4500000e-003'
        >>> float_to_sci(0.0)
        '0.0000000e+000'
        >>> float_to_sci(float('nan'))
        '0.0000000e+000'

    v1 source: bct_generator.py lines 59-69 (identical in 5 other files).
    """
    s = "%.*e" % (precision, value)
    # NaN safety — v1 replaces NaN with 0
    if s == 'nan':
        s = "%.*e" % (precision, 0)
    mantissa, exp = s.split('e')
    # exp_digits + 1 accounts for the sign character (+/-)
    return "%se%+0*d" % (mantissa, exp_digits + 1, int(exp))


def list_to_sci(
    values: Union[List[float], np.ndarray],
    precision: int = 7,
    exp_digits: int = 3,
) -> List[str]:
    """
    Convert a list/array of floats to scientific notation strings.

    Args:
        values:     Iterable of floats.
        precision:  Digits after decimal point in mantissa.
        exp_digits: Minimum digits in exponent.

    Returns:
        List of formatted strings.

    Example:
        >>> list_to_sci([0.0, 20.0, 40.0])
        ['0.0000000e+000', '2.0000000e+001', '4.0000000e+001']

    v1 source: bct_generator.py lines 71-78.
    """
    return [float_to_sci(v, precision, exp_digits) for v in values]


# ─────────────────────────────────────────────────────────────────────────────
# FIXED-POINT NOTATION — used by .bcw (SWAN wave boundary format)
# ─────────────────────────────────────────────────────────────────────────────

def float_to_fixed(value: float, decimal_digits: int = 4) -> str:
    """
    Convert a float to fixed-point string for SWAN format.

    Args:
        value:          The float to convert.
        decimal_digits: Number of digits after the decimal point.

    Returns:
        Fixed-point string.

    Examples:
        >>> float_to_fixed(1.5, 4)
        '1.5000'
        >>> float_to_fixed(300.0, 2)
        '300.00'
        >>> float_to_fixed(-0.123, 4)
        '-0.1230'

    v1 source: bcw_generator.py lines 231-236.
    """
    return format(value, f'.{decimal_digits}f')


def list_to_fixed(
    values: Union[List[float], np.ndarray],
    decimal_digits: int = 4,
) -> List[str]:
    """
    Convert a list/array of floats to fixed-point strings.

    Args:
        values:         Iterable of floats.
        decimal_digits: Number of digits after the decimal point.

    Returns:
        List of fixed-point strings.

    v1 source: bcw_generator.py lines 231-236 (called as convert_float_fstr).
    """
    return [float_to_fixed(v, decimal_digits) for v in values]


# ─────────────────────────────────────────────────────────────────────────────
# COLUMN PADDING — used by .bct and .bcw for Delft3D column alignment
# ─────────────────────────────────────────────────────────────────────────────

def pad_positive(value_str: str) -> str:
    """
    Add a leading blank space for positive values.

    Delft3D .bct format requires columns to be aligned. Negative values
    naturally have a '-' sign, so positive values need a leading space
    to maintain alignment.

    Args:
        value_str: Number string (e.g., '1.2345678e+000' or '-3.45e-003').

    Returns:
        Input string with leading space if positive, unchanged if negative.

    Examples:
        >>> pad_positive('1.2345678e+000')
        ' 1.2345678e+000'
        >>> pad_positive('-3.4500000e-003')
        '-3.4500000e-003'

    v1 source: bct_generator.py lines 80-87.
    """
    if not value_str.startswith('-'):
        return f' {value_str}'
    return value_str


def pad_to_integral_width(value_str: str, min_integral_width: int) -> str:
    """
    Pad a number string with leading spaces so the integer part
    (everything before the decimal point) occupies at least
    `min_integral_width` characters.

    Used in .bcw (SWAN) format for column alignment of time, wave height,
    period, direction, and directional spreading columns, where each column
    has a different expected integral width.

    Args:
        value_str:          Number string (e.g., '1.5000', '-12.34', '300.00').
        min_integral_width: Minimum character width of the integral part.

    Returns:
        Padded string.

    Examples:
        >>> pad_to_integral_width('1.5000', 4)
        '   1.5000'
        >>> pad_to_integral_width('-12.34', 4)
        ' -12.34'
        >>> pad_to_integral_width('300.00', 4)
        ' 300.00'
        >>> pad_to_integral_width('12345.00', 4)
        '12345.00'

    v1 source: bcw_generator.py lines 64-83.
    Note: v1 version had commented-out sign-handling code. This version
    purely pads based on integral width, matching the active v1 behaviour.
    """
    integral_len = len(value_str.split('.')[0])
    while integral_len < min_integral_width:
        value_str = f' {value_str}'
        integral_len = len(value_str.split('.')[0])
    return value_str


# ─────────────────────────────────────────────────────────────────────────────
# SELF-TEST
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    print("=== formatting_utils.py self-test ===\n")

    # ── Scientific notation ──────────────────────────────────────────────
    tests_sci = [
        (20.0,         '2.0000000e+001'),
        (0.0,          '0.0000000e+000'),
        (-0.00345,     '-3.4500000e-003'),
        (1234567.89,   '1.2345679e+006'),
        (float('nan'), '0.0000000e+000'),
        (-1.0,         '-1.0000000e+000'),
        (0.12345678,   '1.2345678e-001'),
    ]
    print("Scientific notation (precision=7, exp_digits=3):")
    all_pass = True
    for val, expected in tests_sci:
        result = float_to_sci(val)
        status = "PASS" if result == expected else "FAIL"
        if status == "FAIL":
            all_pass = False
        print(f"  float_to_sci({val:>15}) = '{result}'  expected '{expected}'  [{status}]")

    # List version
    result_list = list_to_sci([0.0, 20.0, 40.0])
    expected_list = ['0.0000000e+000', '2.0000000e+001', '4.0000000e+001']
    status = "PASS" if result_list == expected_list else "FAIL"
    if status == "FAIL":
        all_pass = False
    print(f"  list_to_sci([0, 20, 40]) = {result_list}  [{status}]")

    # ── Fixed-point notation ─────────────────────────────────────────────
    print("\nFixed-point notation:")
    tests_fixed = [
        ((1.5, 4),    '1.5000'),
        ((300.0, 2),  '300.00'),
        ((-0.123, 4), '-0.1230'),
        ((0.0, 2),    '0.00'),
    ]
    for (val, dec), expected in tests_fixed:
        result = float_to_fixed(val, dec)
        status = "PASS" if result == expected else "FAIL"
        if status == "FAIL":
            all_pass = False
        print(f"  float_to_fixed({val}, {dec}) = '{result}'  expected '{expected}'  [{status}]")

    # ── Positive padding ─────────────────────────────────────────────────
    print("\nPositive value padding:")
    tests_pad = [
        ('1.2345678e+000',  ' 1.2345678e+000'),
        ('-3.4500000e-003', '-3.4500000e-003'),
        ('0.0000000e+000',  ' 0.0000000e+000'),
    ]
    for val, expected in tests_pad:
        result = pad_positive(val)
        status = "PASS" if result == expected else "FAIL"
        if status == "FAIL":
            all_pass = False
        print(f"  pad_positive('{val}') = '{result}'  [{status}]")

    # ── Integral width padding ───────────────────────────────────────────
    print("\nIntegral width padding:")
    tests_intpad = [
        (('1.5000', 4),     '   1.5000'),
        (('-12.34', 4),     ' -12.34'),
        (('300.00', 4),     ' 300.00'),
        (('12345.00', 4),   '12345.00'),
        (('0.00', 2),       ' 0.00'),
    ]
    for (val, width), expected in tests_intpad:
        result = pad_to_integral_width(val, width)
        status = "PASS" if result == expected else "FAIL"
        if status == "FAIL":
            all_pass = False
        print(f"  pad_to_integral_width('{val}', {width}) = '{result}'  [{status}]")

    # ── Summary ──────────────────────────────────────────────────────────
    print(f"\n{'=' * 50}")
    if all_pass:
        print("ALL TESTS PASSED")
    else:
        print("SOME TESTS FAILED — check output above")
    print(f"{'=' * 50}")
