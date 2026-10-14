# -*- coding: utf-8 -*-
"""
gui_theme.py — Central theme definition for RepD3D v2.0.

All colours, fonts, and reusable widget factory functions live here.
Every GUI module imports from this file to ensure visual consistency
across the launcher and all sub-module windows.

Design specifications from handoff document Section 3:
    - Dark mode with custom colour palette
    - Sans-serif font family (Segoe UI / system sans-serif)
    - 4 button colour categories (primary, warning, secondary, neutral)
    - Modern, clean design — functionality over decoration

Usage:
    from gui_theme import (
        apply_theme, COLORS, FONTS,
        make_button, make_label, make_entry, make_optionmenu,
        make_textbox, make_progressbar, make_frame, make_scrollable_frame,
        make_checkbox, make_file_row,
        get_asset_path,
    )
"""

import customtkinter as ctk
from pathlib import Path
from typing import Optional, Callable
import sys

# ─────────────────────────────────────────────────────────────────────────────
# COLOUR PALETTE
# ─────────────────────────────────────────────────────────────────────────────

COLORS = {
    # ── Backgrounds ──────────────────────────────────────────────────────
    "bg":              "#1A1C1E",    # Main window background (deeper charcoal)
    "bg_secondary":    "#242729",    # Slightly lighter for nested frames
    "entry_bg":        "#2C3033",    # Entry / dropdown background
    "entry_border":    "#525A5E",    # Entry border (more visible)
    "log_bg":          "#141618",    # Log panel background (near-black)

    # ── Text ─────────────────────────────────────────────────────────────
    "text":            "#F5F0EB",    # Primary text on dark background
    "text_dim":        "#8A9299",    # Secondary / dimmed text
    "text_on_button":  "#0D0D0D",    # Text on coloured buttons

    # ── Buttons ──────────────────────────────────────────────────────────
    "btn_primary":     "#10D9A0",    # Primary / confirm (vivid mint)
    "btn_warning":     "#FF5A47",    # Warning / cancel (vivid coral)
    "btn_secondary":   "#E0A020",    # Secondary actions (bright amber)
    "btn_neutral":     "#C8C2BE",    # Neutral / info (warm grey)
    "btn_white":       "#FFFFFF",    # Launcher module buttons (white)
    "btn_settings":    "#5AADCF",    # Save/Load settings (saturated teal-blue)
    "btn_reset":       "#E03535",    # Reset/Stop (sharp red)

    # ── Accents ──────────────────────────────────────────────────────────
    "highlight":       "#10D9A0",    # Selection highlight
    "progress_trough": "#313538",    # Progress bar trough
    "separator":       "#3A4044",    # Visible but subdued divider lines
    "success":         "#10D9A0",    # Success messages
    "error":           "#FF5A47",    # Error messages
}

# ─────────────────────────────────────────────────────────────────────────────
# FONTS
# ─────────────────────────────────────────────────────────────────────────────
# Sans-serif family for a modern look.
# Segoe UI: native Windows system font, renders sharply with ClearType.
# macOS falls back to the system sans-serif transparently.

FONT_FAMILY = "Segoe UI"

FONTS = {
    "h1":        (FONT_FAMILY, 22, "bold"),     # Window titles
    "h2":        (FONT_FAMILY, 16, "bold"),     # Section headers
    "h3":        (FONT_FAMILY, 13, "bold"),     # Sub-section headers
    "body":      (FONT_FAMILY, 12),             # Default body text
    "body_sm":   (FONT_FAMILY, 11),             # Smaller body text
    "body_bold": (FONT_FAMILY, 12, "bold"),     # Emphasis
    "mono":      ("Consolas", 11),              # Log output (monospace)
    "mono_sm":   ("Consolas", 10),              # Smaller monospace
    "caption":   (FONT_FAMILY, 10),             # Captions / credits
}


# ─────────────────────────────────────────────────────────────────────────────
# THEME APPLICATION
# ─────────────────────────────────────────────────────────────────────────────

def apply_theme():
    """
    Apply the RepD3D dark theme globally.
    Call once at application startup before creating any windows.
    """
    ctk.set_appearance_mode("dark")
    ctk.set_default_color_theme("dark-blue")


def apply_window_icon(window) -> None:
    """
    Apply the RepD3D icon to a Tk/CTk window.

    Use the PNG as the portable runtime window icon across platforms.
    On Windows, also try the .ico so titlebar/taskbar integration remains
    compatible with native expectations.
    """
    try:
        from PIL import Image, ImageTk
        icon_path = get_asset_path("repd3d_logo.png")
        icon_image = ImageTk.PhotoImage(Image.open(icon_path))
        window.iconphoto(True, icon_image)
        window._repd3d_icon = icon_image
    except Exception:
        icon_image = None

    if sys.platform.startswith("win"):
        def _set_ico():
            try:
                window.iconbitmap(get_asset_path("repd3d.ico"))
            except Exception:
                pass

        _set_ico()
        # CustomTkinter's CTkToplevel schedules its own default icon 200 ms
        # after construction (ctk_toplevel.py, unconditional after(200, ...)),
        # which overwrites the icon set above. Re-apply once that has fired.
        try:
            window.after(400, _set_ico)
        except Exception:
            pass


# ─────────────────────────────────────────────────────────────────────────────
# COLOUR UTILITIES
# ─────────────────────────────────────────────────────────────────────────────

def _darken(hex_color: str, factor: float = 0.15) -> str:
    """Darken a hex colour by a factor (0.0–1.0)."""
    hex_color = hex_color.lstrip("#")
    r = int(int(hex_color[0:2], 16) * (1 - factor))
    g = int(int(hex_color[2:4], 16) * (1 - factor))
    b = int(int(hex_color[4:6], 16) * (1 - factor))
    return f"#{max(r,0):02x}{max(g,0):02x}{max(b,0):02x}"


def _lighten(hex_color: str, factor: float = 0.15) -> str:
    """Lighten a hex colour by a factor (0.0–1.0)."""
    hex_color = hex_color.lstrip("#")
    r = min(int(int(hex_color[0:2], 16) * (1 + factor)), 255)
    g = min(int(int(hex_color[2:4], 16) * (1 + factor)), 255)
    b = min(int(int(hex_color[4:6], 16) * (1 + factor)), 255)
    return f"#{r:02x}{g:02x}{b:02x}"


# ─────────────────────────────────────────────────────────────────────────────
# WIDGET FACTORIES
# ─────────────────────────────────────────────────────────────────────────────
# These ensure every widget across all windows looks identical.
# All factories return the widget AND grid/pack it — the caller only
# needs to position with .grid() or .pack() on the returned widget if
# they want to override.

def make_button(
    parent,
    text: str,
    command: Callable,
    style: str = "primary",
    width: int = 160,
    **kwargs,
) -> ctk.CTkButton:
    """
    Create a themed button.

    Args:
        parent:  Parent widget.
        text:    Button label.
        command: Callback function.
        style:   'primary' | 'warning' | 'secondary' | 'neutral' | 'settings' | 'reset'
        width:   Button width in pixels.

    Returns:
        CTkButton instance (not yet packed/gridded).
    """
    color_key = f"btn_{style}"
    fg_color = COLORS.get(color_key, COLORS["btn_primary"])
    text_color = COLORS["text"] if style == "reset" else COLORS["text_on_button"]
    return ctk.CTkButton(
        parent,
        text=text,
        command=command,
        width=width,
        fg_color=fg_color,
        text_color=text_color,
        hover_color=_darken(fg_color, 0.15),
        font=FONTS["body_bold"],
        corner_radius=8,
        **kwargs,
    )


def make_label(
    parent,
    text: str,
    font_key: str = "body",
    color: Optional[str] = None,
    **kwargs,
) -> ctk.CTkLabel:
    """
    Create a themed label.

    Args:
        parent:   Parent widget.
        text:     Label text.
        font_key: Key from FONTS dict (default 'body').
        color:    Override text colour (default: COLORS["text"]).
    """
    text_color = color or COLORS["text"]
    return ctk.CTkLabel(
        parent,
        text=text,
        text_color=text_color,
        font=FONTS.get(font_key, FONTS["body"]),
        **kwargs,
    )


def make_entry(
    parent,
    placeholder: str = "",
    width: int = 300,
    **kwargs,
) -> ctk.CTkEntry:
    """Create a themed text entry."""
    return ctk.CTkEntry(
        parent,
        placeholder_text=placeholder,
        width=width,
        fg_color=COLORS["entry_bg"],
        border_color=COLORS["entry_border"],
        text_color=COLORS["text"],
        font=FONTS["body"],
        corner_radius=6,
        **kwargs,
    )


def make_optionmenu(
    parent,
    values: list,
    command: Optional[Callable] = None,
    width: int = 200,
    **kwargs,
) -> ctk.CTkOptionMenu:
    """Create a themed dropdown menu."""
    return ctk.CTkOptionMenu(
        parent,
        values=values,
        command=command,
        width=width,
        fg_color=COLORS["entry_bg"],
        button_color=COLORS["btn_secondary"],
        button_hover_color=_darken(COLORS["btn_secondary"], 0.15),
        text_color=COLORS["text"],
        font=FONTS["body"],
        dropdown_fg_color=COLORS["bg_secondary"],
        dropdown_text_color=COLORS["text"],
        dropdown_hover_color=_darken(COLORS["btn_secondary"], 0.3),
        corner_radius=6,
        **kwargs,
    )


def make_textbox(
    parent,
    height: int = 300,
    **kwargs,
) -> ctk.CTkTextbox:
    """
    Create a themed read-only textbox for log output.
    State is set to 'disabled' — use gui_utils.append_log() to write.
    """
    return ctk.CTkTextbox(
        parent,
        height=height,
        fg_color=COLORS["log_bg"],
        text_color=COLORS["text"],
        font=FONTS["mono"],
        state="disabled",
        wrap="word",
        corner_radius=6,
        **kwargs,
    )


def make_progressbar(
    parent,
    width: int = 300,
    **kwargs,
) -> ctk.CTkProgressBar:
    """Create a themed progress bar. Call .set(0.0–1.0) to update."""
    bar = ctk.CTkProgressBar(
        parent,
        width=width,
        fg_color=COLORS["progress_trough"],
        progress_color=COLORS["btn_primary"],
        corner_radius=4,
        **kwargs,
    )
    bar.set(0)
    return bar


def make_frame(parent, **kwargs) -> ctk.CTkFrame:
    """Create a themed frame with default background."""
    return ctk.CTkFrame(
        parent,
        fg_color=COLORS["bg"],
        corner_radius=0,
        **kwargs,
    )


def make_scrollable_frame(
    parent,
    label_text: str = "",
    **kwargs,
) -> ctk.CTkScrollableFrame:
    """Create a themed scrollable frame."""
    return ctk.CTkScrollableFrame(
        parent,
        fg_color=COLORS["bg"],
        label_text=label_text,
        label_font=FONTS["h2"] if label_text else None,
        label_text_color=COLORS["text"] if label_text else None,
        corner_radius=0,
        **kwargs,
    )


def make_checkbox(
    parent,
    text: str,
    variable: Optional[ctk.Variable] = None,
    command: Optional[Callable] = None,
    **kwargs,
) -> ctk.CTkCheckBox:
    """Create a themed checkbox."""
    return ctk.CTkCheckBox(
        parent,
        text=text,
        variable=variable,
        command=command,
        text_color=COLORS["text"],
        font=FONTS["body"],
        fg_color=COLORS["btn_secondary"],
        hover_color=_darken(COLORS["btn_secondary"], 0.15),
        border_color=COLORS["entry_border"],
        checkmark_color=COLORS["text_on_button"],
        corner_radius=4,
        **kwargs,
    )


def make_radiobutton(
    parent,
    text: str,
    variable: ctk.Variable,
    value,
    command: Optional[Callable] = None,
    **kwargs,
) -> ctk.CTkRadioButton:
    """Create a themed radio button."""
    return ctk.CTkRadioButton(
        parent,
        text=text,
        variable=variable,
        value=value,
        command=command,
        text_color=COLORS["text"],
        font=FONTS["body"],
        fg_color=COLORS["btn_secondary"],
        hover_color=_darken(COLORS["btn_secondary"], 0.15),
        border_color=COLORS["entry_border"],
        **kwargs,
    )


def make_separator(parent, orientation: str = "horizontal") -> ctk.CTkFrame:
    """Create a thin separator line."""
    if orientation == "horizontal":
        sep = ctk.CTkFrame(
            parent,
            height=2,
            fg_color=COLORS["separator"],
            corner_radius=0,
            border_width=0,
        )
    else:
        sep = ctk.CTkFrame(
            parent,
            width=2,
            fg_color=COLORS["separator"],
            corner_radius=0,
            border_width=0,
        )
    sep.pack_propagate(False)
    return sep


# ─────────────────────────────────────────────────────────────────────────────
# COMPOSITE WIDGETS — common patterns across sub-module GUIs
# ─────────────────────────────────────────────────────────────────────────────

def make_step_header(
    parent,
    step_number: int,
    title: str,
) -> ctk.CTkLabel:
    """
    Create a step header label like "Step 1: Model Type".
    Used in the left panel workflow of every sub-module.
    """
    return make_label(
        parent,
        text=f"Step {step_number}:  {title}",
        font_key="h3",
    )


def make_file_row(
    parent,
    label_text: str,
    button_text: str = "Browse",
    command: Optional[Callable] = None,
) -> tuple:
    """
    Create a standard file-selection row: [Label] [Entry] [Browse Button]

    Returns:
        (frame, entry_widget, button_widget)
        Caller should .pack() or .grid() the frame.
    """
    row_frame = ctk.CTkFrame(parent, fg_color="transparent")

    label = make_label(row_frame, label_text, font_key="body_sm")
    label.pack(side="left", padx=(0, 8))

    entry = make_entry(row_frame, width=280, placeholder=f"Select {label_text.lower()}...")
    entry.pack(side="left", padx=(0, 8), fill="x", expand=True)

    btn = make_button(
        row_frame, button_text, command=command or (lambda: None),
        style="neutral", width=80,
    )
    btn.pack(side="left")

    return row_frame, entry, btn


def make_status_label(parent, initial_text: str = "Ready") -> ctk.CTkLabel:
    """Create a status label (e.g., for file selection feedback)."""
    return make_label(parent, initial_text, font_key="body_sm", color=COLORS["text_dim"])


# ─────────────────────────────────────────────────────────────────────────────
# ASSET PATH RESOLUTION
# ─────────────────────────────────────────────────────────────────────────────

def get_asset_path(filename: str) -> str:
    """
    Get the absolute path to an asset file (logo, images, etc.).

    Works both during development AND when frozen with auto-py-to-exe
    (PyInstaller sets sys._MEIPASS to the temp extraction directory).

    Args:
        filename: Asset filename (e.g., 'repd3d_logo.png').

    Returns:
        Absolute path string.

    v1 source: gui_generator.py lines 81-88 (resource_path method).
    """
    if getattr(sys, 'frozen', False):
        # Running as compiled executable
        base = Path(sys._MEIPASS)
    else:
        # Running in development
        base = Path(__file__).parent
    return str(base / filename)


# ─────────────────────────────────────────────────────────────────────────────
# SELF-TEST (visual — opens a window)
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    print("=== gui_theme.py visual self-test ===")
    print("Opening test window... Close it to finish.\n")

    apply_theme()

    root = ctk.CTk()
    root.title("gui_theme.py — Widget Gallery")
    root.geometry("700x800")
    root.configure(fg_color=COLORS["bg"])

    # Main scrollable area
    scroll = make_scrollable_frame(root, label_text="Widget Gallery")
    scroll.pack(fill="both", expand=True, padx=10, pady=10)

    # Step headers
    make_step_header(scroll, 1, "Buttons").pack(anchor="w", pady=(10, 5))
    btn_frame = ctk.CTkFrame(scroll, fg_color="transparent")
    btn_frame.pack(fill="x", pady=5)
    make_button(btn_frame, "Primary", lambda: print("primary"), "primary").pack(side="left", padx=5)
    make_button(btn_frame, "Warning", lambda: print("warning"), "warning").pack(side="left", padx=5)
    make_button(btn_frame, "Secondary", lambda: print("secondary"), "secondary").pack(side="left", padx=5)
    make_button(btn_frame, "Neutral", lambda: print("neutral"), "neutral").pack(side="left", padx=5)

    make_separator(scroll).pack(fill="x", pady=10)

    # Labels
    make_step_header(scroll, 2, "Labels & Text").pack(anchor="w", pady=(10, 5))
    make_label(scroll, "H1 Title", "h1").pack(anchor="w")
    make_label(scroll, "H2 Section Header", "h2").pack(anchor="w")
    make_label(scroll, "Body text (default)", "body").pack(anchor="w")
    make_label(scroll, "Dimmed caption text", "caption", color=COLORS["text_dim"]).pack(anchor="w")

    make_separator(scroll).pack(fill="x", pady=10)

    # Entries
    make_step_header(scroll, 3, "Input Widgets").pack(anchor="w", pady=(10, 5))
    make_entry(scroll, placeholder="Type something here...").pack(anchor="w", pady=3)
    make_optionmenu(scroll, ["EasyGSH", "TrilaWatt", "Copernicus", "Custom"]).pack(anchor="w", pady=3)

    # Checkboxes & radio buttons
    cb_var = ctk.BooleanVar()
    make_checkbox(scroll, "Enable feature X", variable=cb_var).pack(anchor="w", pady=3)
    rb_var = ctk.StringVar(value="A")
    make_radiobutton(scroll, "Option A", rb_var, "A").pack(anchor="w", pady=2)
    make_radiobutton(scroll, "Option B", rb_var, "B").pack(anchor="w", pady=2)

    make_separator(scroll).pack(fill="x", pady=10)

    # File row
    make_step_header(scroll, 4, "File Selection").pack(anchor="w", pady=(10, 5))
    file_row, file_entry, file_btn = make_file_row(scroll, ".mdf file")
    file_row.pack(fill="x", pady=3)

    make_separator(scroll).pack(fill="x", pady=10)

    # Log output
    make_step_header(scroll, 5, "Log & Progress").pack(anchor="w", pady=(10, 5))
    log = make_textbox(scroll, height=120)
    log.pack(fill="x", pady=5)
    # Write some sample text
    log.configure(state="normal")
    log.insert("end", "12:34:56  Loading EasyGSH data...\n")
    log.insert("end", "12:34:57  Extracting boundary conditions...\n")
    log.insert("end", "12:34:58  Processing 8 boundary points...\n")
    log.insert("end", "12:35:02  Writing output.bct\n")
    log.insert("end", "12:35:02  Complete!\n")
    log.configure(state="disabled")

    prog = make_progressbar(scroll, width=400)
    prog.pack(anchor="w", pady=5)
    prog.set(0.65)

    make_separator(scroll).pack(fill="x", pady=10)

    # Status label
    make_status_label(scroll, "✓ 3 files loaded successfully").pack(anchor="w", pady=5)

    # Colour palette display
    make_step_header(scroll, 6, "Colour Palette").pack(anchor="w", pady=(10, 5))
    palette_frame = ctk.CTkFrame(scroll, fg_color="transparent")
    palette_frame.pack(fill="x", pady=5)
    for name, color in COLORS.items():
        swatch = ctk.CTkFrame(palette_frame, width=50, height=25, fg_color=color, corner_radius=4)
        swatch.pack(side="left", padx=2, pady=2)

    print("All widgets created. Visual inspection required.")
    print(f"Asset path test: {get_asset_path('repd3d_logo.png')}")

    root.mainloop()
    print("Window closed. Self-test complete.")
