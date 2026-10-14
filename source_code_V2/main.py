# -*- coding: utf-8 -*-
"""
main.py — RepD3D v2.0 Launcher Window.

This is the entry point for the application. It displays the logo,
version info, and four buttons that open the sub-module windows.

Run this file to start RepD3D:
    python main.py

Each sub-module opens as a CTkToplevel window. The launcher remains
visible, and only one sub-module is opened at a time to keep long-
running tasks isolated.
"""

import customtkinter as ctk

from gui_theme import (
    apply_theme, COLORS, FONTS,
    make_button, make_label, make_separator,
    get_asset_path, apply_window_icon,
)

# Import sub-module windows
from gui_bc_extraction import BCExtractionWindow
from gui_wind_generator import WindGeneratorWindow
from gui_rep_period import RepPeriodWindow
from gui_sea_level import SeaLevelWindow
from gui_synthetic_bc import SyntheticBCWindow

# ─── Version ────────────────────────────────────────────────────────────────
VERSION = "2.0"


# ─────────────────────────────────────────────────────────────────────────────
# LAUNCHER WINDOW
# ─────────────────────────────────────────────────────────────────────────────

class RepD3DLauncher(ctk.CTk):
    """Main launcher window for RepD3D v2."""

    def __init__(self):
        super().__init__()

        self.title("RepD3D")
        self.geometry("420x720")
        self.resizable(False, False)
        self.configure(fg_color=COLORS["bg"])
        apply_window_icon(self)

        # Track open windows to avoid duplicates
        self._open_windows = {}
        self._module_buttons = {}
        self._closing = False
        self.protocol("WM_DELETE_WINDOW", self._on_app_close)

        self._build_ui()

    def _build_ui(self):
        # ── Logo ─────────────────────────────────────────────────────────
        # NOTE: We use ImageTk.PhotoImage instead of ctk.CTkImage because
        # CTkImage triggers a TclError ("image pyimageN doesn't exist") in
        # this customtkinter/Pillow version combo. The fixed 180x180 size
        # means we don't need CTkImage's HiDPI auto-scaling anyway.
        logo_path = get_asset_path("repd3d_logo.png")
        try:
            import warnings
            from PIL import Image, ImageTk
            pil_img = Image.open(logo_path).convert("RGBA").resize(
                (180, 180), Image.LANCZOS
            )
            logo_img = ImageTk.PhotoImage(pil_img, master=self)
            self._launcher_logo = logo_img
            with warnings.catch_warnings():
                warnings.filterwarnings(
                    "ignore",
                    message=r".*Given image is not CTkImage.*",
                    category=UserWarning,
                )
                logo_label = ctk.CTkLabel(self, image=logo_img, text="")
            logo_label.pack(pady=(25, 5))
        except Exception as e:
            # Logo file not found or PIL not available — show text instead
            print(f"[WARN] Logo failed to load from {logo_path!r}: {type(e).__name__}: {e}")
            make_label(self, "RepD3D", font_key="h1").pack(pady=(25, 5))

        # ── Title & version ──────────────────────────────────────────────
        make_label(self, f"RepD3D v{VERSION}", font_key="h2").pack(pady=(5, 2))
        make_label(
            self,
            "Representative Period Identifier\n& Boundary Condition Extractor for Delft3D",
            font_key="body_sm",
            color=COLORS["text_dim"],
            justify="center",
        ).pack(pady=(0, 20))

        make_separator(self).pack(fill="x", padx=30, pady=5)

        # ── Module buttons ───────────────────────────────────────────────
        btn_frame = ctk.CTkFrame(self, fg_color="transparent")
        btn_frame.pack(pady=10, padx=30, fill="x")

        modules = [
            ("Representative Period",         "white", RepPeriodWindow),
            ("Boundary Condition Extraction", "white", BCExtractionWindow),
            ("Meteo File Generator",          "white", WindGeneratorWindow),
            ("Sea Level Change",              "white", SeaLevelWindow),
            ("Generate Synthetic Files",      "white", SyntheticBCWindow),
        ]

        for text, style, window_class in modules:
            btn = make_button(
                btn_frame, text,
                command=lambda wc=window_class: self._open_module(wc),
                style=style, width=340,
            )
            btn.pack(pady=6)
            self._module_buttons[window_class.__name__] = btn

        # ── Credits ──────────────────────────────────────────────────────
        make_separator(self).pack(fill="x", padx=30, pady=(15, 5))

        credits_frame = ctk.CTkFrame(self, fg_color="transparent")
        credits_frame.pack(pady=(5, 15), padx=30)

        make_label(
            credits_frame,
            "Created by Clayton Soares\n"
            "Co-developer: Arne Knies\n"
            "github.com/capt-clay10/RepD3D",
            font_key="caption",
            color=COLORS["text_dim"],
            justify="center",
        ).pack()

        make_separator(credits_frame).pack(fill="x", pady=(10, 6))

        make_label(
            credits_frame,
            "If you use RepD3D, please cite:",
            font_key="caption",
            color=COLORS["text_dim"],
            justify="center",
        ).pack()

        make_label(
            credits_frame,
            "Soares, C. C., Knies, A., & Winter, C. (2025). RepD3D: A tool for "
            "representative period identification and associated boundary condition "
            "extraction. MethodsX, 14, 103109.",
            font_key="caption",
            color=COLORS["text_dim"],
            justify="center",
            wraplength=355,
        ).pack(pady=(4, 0))

    def _open_module(self, window_class):
        """
        Open a sub-module window. If one of the same type is already open,
        bring it to front instead of creating a duplicate.
        """
        if self._closing:
            return

        class_name = window_class.__name__

        # Only allow one sub-module window at a time to avoid overlapping
        # long-running tasks and to keep the workflow predictable.
        if self._open_windows and class_name not in self._open_windows:
            existing = next(iter(self._open_windows.values()))
            try:
                existing.focus()
                existing.lift()
            except Exception:
                pass
            return

        # Check if window already exists and is still open
        if class_name in self._open_windows:
            existing = self._open_windows[class_name]
            try:
                existing.focus()
                existing.lift()
                return
            except Exception:
                # Window was destroyed — remove reference
                del self._open_windows[class_name]

        # Create new window
        window = window_class(self)
        self._open_windows[class_name] = window
        self._set_module_button_lock(class_name)

        # Clean up reference when window is closed
        window.protocol(
            "WM_DELETE_WINDOW",
            lambda: self._on_module_close(class_name, window),
        )

    def _set_module_button_lock(self, active_class_name: str | None = None):
        """Disable launcher buttons while a sub-module window is open."""
        any_open = active_class_name is not None
        for class_name, button in self._module_buttons.items():
            try:
                button.configure(state="disabled" if any_open else "normal")
            except Exception:
                pass

    def _on_module_close(self, class_name: str, window):
        """Handle sub-module window close."""
        self._cancel_window_task(window)
        self._open_windows.pop(class_name, None)
        if not self._closing:
            self._set_module_button_lock(None)
        try:
            if window.winfo_exists():
                window.destroy()
        except Exception:
            pass

    def _cancel_window_task(self, window):
        """Request cooperative cancellation for a sub-window task, if running."""
        task = getattr(window, "_task", None)
        if task is not None and getattr(task, "is_running", False):
            try:
                task.cancel()
            except Exception:
                pass

    def _on_app_close(self):
        """
        Shut down child windows first, then destroy the launcher.

        This reduces the chance of Tk/CustomTkinter background callbacks
        firing against widgets that have already been torn down.
        """
        if self._closing:
            return

        self._closing = True

        try:
            self.withdraw()
        except Exception:
            pass

        for class_name, window in list(self._open_windows.items()):
            self._cancel_window_task(window)
            try:
                window.withdraw()
            except Exception:
                pass
            try:
                if window.winfo_exists():
                    window.destroy()
            except Exception:
                pass
            self._open_windows.pop(class_name, None)

        self._set_module_button_lock(None)

        try:
            self.after(25, self._finalize_close)
        except Exception:
            self._finalize_close()

    def _finalize_close(self):
        """Finish the shutdown after child windows have been closed."""
        try:
            self.quit()
        except Exception:
            pass
        try:
            if self.winfo_exists():
                self.destroy()
        except Exception:
            pass


# ─────────────────────────────────────────────────────────────────────────────
# ENTRY POINT
# ─────────────────────────────────────────────────────────────────────────────

def main():
    apply_theme()
    app = RepD3DLauncher()
    app.mainloop()


if __name__ == '__main__':
    main()
