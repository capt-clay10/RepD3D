# -*- coding: utf-8 -*-
"""
gui_rep_period.py — Representative Period Identifier GUI window.

Step-by-step workflow:
    1. Load wind CSV file
    2. Set total reference period
    3. Set scanning period (subset)
    4. Define search window sizes
    5. Select directional sectors
    6. Select speed classes
    7. Set output options
    8. Run identification

Can be run standalone: python gui_rep_period.py
"""

import customtkinter as ctk
from pathlib import Path
from typing import Optional

from gui_theme import (
    apply_theme, COLORS,
    make_button, make_label, make_entry, make_optionmenu,
    make_separator, make_step_header, make_file_row, make_status_label,
    make_checkbox, make_frame, get_asset_path,
)
from gui_utils import (
    create_sub_window_layout, BackgroundTask,
    browse_file, browse_directory,
    append_log, clear_log,
    save_json_settings, load_json_settings,
    save_settings_dialog, load_settings_dialog,
    show_success_in_log, show_error_dialog,
    set_widgets_state,
)


MODULE_ID = "rep_period"
SETTINGS_VERSION = 1


# All 16 compass sectors
ALL_SECTORS = [
    'N', 'NNE', 'NE', 'ENE', 'E', 'ESE', 'SE', 'SSE',
    'S', 'SSW', 'SW', 'WSW', 'W', 'WNW', 'NW', 'NNW',
]


class RepPeriodWindow(ctk.CTkToplevel):
    """Representative Period Identifier sub-module window."""

    def __init__(self, parent=None):
        super().__init__(parent)
        left, self.log_textbox, self.progress_bar, self.status_label, self.eta_label = \
            create_sub_window_layout(self, title="Representative Period Identifier")

        # Shrink log textbox to make room for the classes image below
        right_frame = self.log_textbox.master
        self.log_textbox.configure(height=300)
        right_frame.grid_rowconfigure(2, minsize=300)
        self._build_classes_image(right_frame)

        self._interactive_widgets = []
        self._sector_vars = {}
        self._speed_vars = {}
        self._task: Optional[BackgroundTask] = None
        self._build_workflow(left)

    def _build_classes_image(self, right_frame: ctk.CTkFrame) -> None:
        """Load and display classes.jpg at the bottom of the right panel."""
        try:
            from PIL import Image
            pil_img = Image.open(get_asset_path("classes.jpg"))
            max_w, max_h = 760, 320
            w, h = pil_img.size
            scale = min(max_w / w, max_h / h)
            display_size = (int(w * scale), int(h * scale))
            ctk_img = ctk.CTkImage(
                light_image=pil_img,
                dark_image=pil_img,
                size=display_size,
            )
            self._classes_image = ctk_img  # prevent garbage collection
            make_separator(right_frame).grid(
                row=6, column=0, sticky="ew", padx=14, pady=(4, 6),
            )
            ctk.CTkLabel(right_frame, image=ctk_img, text="").grid(
                row=7, column=0, sticky="", padx=14, pady=(0, 14),
            )
        except Exception:
            pass  # silently skip if PIL unavailable or file missing

    def _build_workflow(self, parent):
        pad = dict(anchor="w", padx=10)

        # ─── Step 1: Wind Data ───────────────────────────────────────────
        make_step_header(parent, 1, "Wind Data").pack(**pad, pady=(15, 5))
        make_label(
            parent,
            "CSV format: date (YYYYMMDDHH), speed (m/s), direction (°)",
            font_key="body_sm", color=COLORS["text_dim"],
        ).pack(**pad)
        row_wind, self.wind_entry, self.wind_btn = make_file_row(parent, "Wind CSV file")
        self.wind_btn.configure(command=lambda: browse_file(
            self.wind_entry, "Select wind CSV",
            [("CSV/Text files", "*.csv *.txt"), ("All files", "*.*")],
        ))
        row_wind.pack(fill="x", padx=10, pady=3)
        self._interactive_widgets.append(self.wind_btn)

        make_separator(parent).pack(fill="x", padx=10, pady=(12, 8))

        # ─── Step 2: Total Reference Period ──────────────────────────────
        make_step_header(parent, 2, "Total Reference Period").pack(**pad, pady=(0, 5))
        ref_frame = ctk.CTkFrame(parent, fg_color="transparent")
        ref_frame.pack(fill="x", padx=10, pady=3)
        make_label(ref_frame, "Start:", font_key="body_sm").pack(side="left", padx=(0, 5))
        self.ref_start = make_entry(ref_frame, placeholder="1995-01-01 00:00:00", width=180)
        self.ref_start.pack(side="left", padx=(0, 15))
        make_label(ref_frame, "End:", font_key="body_sm").pack(side="left", padx=(0, 5))
        self.ref_end = make_entry(ref_frame, placeholder="2020-12-31 00:00:00", width=180)
        self.ref_end.pack(side="left")
        self._interactive_widgets.extend([self.ref_start, self.ref_end])

        make_separator(parent).pack(fill="x", padx=10, pady=(12, 8))

        # ─── Step 3: Scanning Period ─────────────────────────────────────
        make_step_header(parent, 3, "Scanning Period").pack(**pad, pady=(0, 5))
        make_label(parent, "The algorithm scans within this subset for candidates",
                   font_key="body_sm", color=COLORS["text_dim"]).pack(**pad)
        scan_frame = ctk.CTkFrame(parent, fg_color="transparent")
        scan_frame.pack(fill="x", padx=10, pady=3)
        make_label(scan_frame, "Start:", font_key="body_sm").pack(side="left", padx=(0, 5))
        self.scan_start = make_entry(scan_frame, placeholder="2000-01-01 00:00:00", width=180)
        self.scan_start.pack(side="left", padx=(0, 15))
        make_label(scan_frame, "End:", font_key="body_sm").pack(side="left", padx=(0, 5))
        self.scan_end = make_entry(scan_frame, placeholder="2018-01-01 00:00:00", width=180)
        self.scan_end.pack(side="left")
        self._interactive_widgets.extend([self.scan_start, self.scan_end])

        make_separator(parent).pack(fill="x", padx=10, pady=(12, 8))

        # ─── Step 4: Search Windows ──────────────────────────────────────
        make_step_header(parent, 4, "Search Window Sizes").pack(**pad, pady=(0, 5))
        make_label(parent, "Comma-separated, e.g.: 1MS, 3MS, 6MS, 12MS",
                   font_key="body_sm", color=COLORS["text_dim"]).pack(**pad)
        make_label(parent, "MS = calendar months, W = weeks, D = days. Windows end by the scan end.",
                   font_key="body_sm", color=COLORS["text_dim"]).pack(**pad)
        self.freq_entry = make_entry(parent, placeholder="1MS, 3MS, 6MS", width=300)
        self.freq_entry.pack(**pad, pady=3)
        self._interactive_widgets.append(self.freq_entry)

        make_separator(parent).pack(fill="x", padx=10, pady=(12, 8))

        # ─── Step 5: Directional Sectors ─────────────────────────────────
        make_step_header(parent, 5, "Directional Sectors").pack(**pad, pady=(0, 5))

        # Select All / Deselect All buttons
        sec_btn_frame = ctk.CTkFrame(parent, fg_color="transparent")
        sec_btn_frame.pack(fill="x", padx=10, pady=2)
        make_button(sec_btn_frame, "Select All", self._select_all_sectors,
                    style="neutral", width=100).pack(side="left", padx=(0, 5))
        make_button(sec_btn_frame, "Deselect All", self._deselect_all_sectors,
                    style="neutral", width=100).pack(side="left")

        # Checkboxes in a grid (4 columns)
        sec_frame = ctk.CTkFrame(parent, fg_color="transparent")
        sec_frame.pack(fill="x", padx=10, pady=3)
        for i, sector in enumerate(ALL_SECTORS):
            var = ctk.BooleanVar(value=True)
            self._sector_vars[sector] = var
            cb = make_checkbox(sec_frame, sector, variable=var)
            cb.grid(row=i // 4, column=i % 4, sticky="w", padx=5, pady=1)

        make_separator(parent).pack(fill="x", padx=10, pady=(12, 8))

        # ─── Step 6: Speed Classes ───────────────────────────────────────
        make_step_header(parent, 6, "Speed Classes (Beaufort)").pack(**pad, pady=(0, 5))
        spd_frame = ctk.CTkFrame(parent, fg_color="transparent")
        spd_frame.pack(fill="x", padx=10, pady=3)
        for i in range(12):
            var = ctk.BooleanVar(value=(i <= 9))  # Default: classes 0-9
            self._speed_vars[i] = var
            cb = make_checkbox(spd_frame, str(i), variable=var)
            cb.grid(row=i // 6, column=i % 6, sticky="w", padx=5, pady=1)

        make_separator(parent).pack(fill="x", padx=10, pady=(12, 8))

        # ─── Step 7: Output ──────────────────────────────────────────────
        make_step_header(parent, 7, "Output").pack(**pad, pady=(0, 5))
        name_frame = ctk.CTkFrame(parent, fg_color="transparent")
        name_frame.pack(fill="x", padx=10, pady=3)
        make_label(name_frame, "Name:", font_key="body_sm").pack(side="left", padx=(0, 5))
        self.output_name_entry = make_entry(name_frame, placeholder="station_name", width=200)
        self.output_name_entry.pack(side="left")

        row_out, self.output_dir_entry, self.output_dir_btn = make_file_row(parent, "Output directory")
        self.output_dir_btn.configure(command=lambda: browse_directory(self.output_dir_entry))
        row_out.pack(fill="x", padx=10, pady=3)
        self._interactive_widgets.extend([self.output_name_entry, self.output_dir_btn])

        make_separator(parent).pack(fill="x", padx=10, pady=(12, 8))

        # ─── Step 8: Run ─────────────────────────────────────────────────
        make_step_header(parent, 8, "Run").pack(**pad, pady=(0, 5))
        self.run_btn = make_button(
            parent, "Find Representative Periods",
            command=self._on_run, style="primary", width=300,
        )
        self.run_btn.pack(**pad, pady=(5, 5))

        # ─── Settings / Reset row ────────────────────────────────────────
        btn_row = ctk.CTkFrame(parent, fg_color="transparent")
        btn_row.pack(fill="x", padx=10, pady=(5, 20))
        self.save_btn = make_button(
            btn_row, "Save Settings", self._on_save_settings,
            style="settings", width=140,
        )
        self.save_btn.pack(side="left", padx=(0, 5))
        self.load_btn = make_button(
            btn_row, "Load Settings", self._on_load_settings,
            style="settings", width=140,
        )
        self.load_btn.pack(side="left", padx=(0, 5))
        self.reset_btn = make_button(
            btn_row, "Reset / Stop", self._on_reset,
            style="reset", width=140,
        )
        self.reset_btn.pack(side="left")
        self._interactive_widgets.extend([self.save_btn, self.load_btn])

    # ─── Helpers ─────────────────────────────────────────────────────────

    def _select_all_sectors(self):
        for var in self._sector_vars.values():
            var.set(True)

    def _deselect_all_sectors(self):
        for var in self._sector_vars.values():
            var.set(False)

    def _get_selected_sectors(self):
        return [s for s, v in self._sector_vars.items() if v.get()]

    def _get_selected_speeds(self):
        return [i for i, v in self._speed_vars.items() if v.get()]

    def _on_run(self):
        wind_file = self.wind_entry.get().strip()
        if not wind_file:
            show_error_dialog(self, "Missing input", "Please select a wind CSV file.")
            return

        ref_start = self.ref_start.get().strip()
        ref_end = self.ref_end.get().strip()
        scan_start = self.scan_start.get().strip()
        scan_end = self.scan_end.get().strip()

        if not all([ref_start, ref_end, scan_start, scan_end]):
            show_error_dialog(self, "Missing input", "Please fill in all date fields.")
            return

        freq_raw = self.freq_entry.get().strip()
        if not freq_raw:
            show_error_dialog(self, "Missing input", "Please enter search window size(s).")
            return
        frequency = [f.strip() for f in freq_raw.split(',') if f.strip()]

        import backend_rep_period as brp
        try:
            for f in frequency:
                brp.parse_window_freq(f)
        except ValueError as exc:
            show_error_dialog(self, "Invalid search window", str(exc))
            return

        sectors = self._get_selected_sectors()
        speeds = self._get_selected_speeds()
        if not sectors:
            show_error_dialog(self, "Missing input", "Please select at least one directional sector.")
            return
        if not speeds:
            show_error_dialog(self, "Missing input", "Please select at least one speed class.")
            return

        output_name = self.output_name_entry.get().strip() or "rep_period"
        output_dir = self.output_dir_entry.get().strip() or str(Path(wind_file).parent)

        # Stash for the windrose popup triggered on completion
        self._last_run_inputs = {
            "wind_file": wind_file,
            "ref_start": ref_start,
            "ref_end": ref_end,
        }

        set_widgets_state(self._interactive_widgets + [self.run_btn], "disabled")
        self.status_label.configure(text="Scanning — this may take several minutes...")

        self._task = BackgroundTask(
            target=brp.identify_rep_period,
            kwargs=dict(
                file_input=wind_file,
                output_name=output_name,
                quad=sectors,
                spd=speeds,
                start_time_total=ref_start,
                end_time_total=ref_end,
                frequency=frequency,
                start_time=scan_start,
                end_time=scan_end,
                out_path=output_dir,
            ),
            log_textbox=self.log_textbox,
            progress_bar=self.progress_bar,
            eta_label=self.eta_label,
            on_complete=self._on_complete,
            on_error=self._on_error,
            on_cancel=self._on_cancel,
            root=self,
        )
        self._task.start()

    def _on_complete(self, result):
        n = result.get('n_candidates', '?')
        show_success_in_log(self.log_textbox, f"Complete — {n} candidates ranked")
        self.status_label.configure(text=f"✓ Complete — {n} candidates → {Path(result.get('output_file','')).name}")
        set_widgets_state(self._interactive_widgets + [self.run_btn], "normal")
        self._open_windrose_popup(result)

    def _open_windrose_popup(self, result):
        inputs = getattr(self, "_last_run_inputs", None)
        out_file = result.get("output_file")
        if not inputs or not out_file:
            return
        try:
            from gui_windrose import WindroseCompareWindow, load_top_prospects
            prospects = load_top_prospects(out_file, k=20)
            if not prospects:
                return
            WindroseCompareWindow(
                self,
                wind_csv_path=inputs["wind_file"],
                ref_start=inputs["ref_start"],
                ref_end=inputs["ref_end"],
                prospect_periods=prospects,
            )
        except Exception as exc:
            show_error_dialog(self, "Windrose Error",
                              f"Could not open windrose window:\n{exc}")

    def _on_error(self, exc):
        self.status_label.configure(text=f"✗ Error: {exc}")
        set_widgets_state(self._interactive_widgets + [self.run_btn], "normal")
        show_error_dialog(self, "Rep Period Error", str(exc))

    def _on_cancel(self):
        self.status_label.configure(text="Cancelled")
        set_widgets_state(self._interactive_widgets + [self.run_btn], "normal")

    # ─── Settings: Save / Load / Reset ───────────────────────────────────

    def _collect_settings(self) -> dict:
        return {
            "module": MODULE_ID,
            "version": SETTINGS_VERSION,
            "values": {
                "wind_file": self.wind_entry.get().strip(),
                "ref_start": self.ref_start.get().strip(),
                "ref_end": self.ref_end.get().strip(),
                "scan_start": self.scan_start.get().strip(),
                "scan_end": self.scan_end.get().strip(),
                "frequency": self.freq_entry.get().strip(),
                "sectors": [s for s, v in self._sector_vars.items() if v.get()],
                "speeds": [i for i, v in self._speed_vars.items() if v.get()],
                "output_name": self.output_name_entry.get().strip(),
                "output_dir": self.output_dir_entry.get().strip(),
            },
        }

    def _apply_settings(self, settings: dict) -> None:
        if settings.get("module") != MODULE_ID:
            raise ValueError(
                f"Settings file is for '{settings.get('module')}', "
                f"not '{MODULE_ID}'"
            )
        v = settings.get("values", {})

        def _set_entry(entry, value):
            entry.delete(0, "end")
            entry.insert(0, str(value))

        _set_entry(self.wind_entry, v.get("wind_file", ""))
        _set_entry(self.ref_start, v.get("ref_start", ""))
        _set_entry(self.ref_end, v.get("ref_end", ""))
        _set_entry(self.scan_start, v.get("scan_start", ""))
        _set_entry(self.scan_end, v.get("scan_end", ""))
        _set_entry(self.freq_entry, v.get("frequency", ""))
        _set_entry(self.output_name_entry, v.get("output_name", ""))
        _set_entry(self.output_dir_entry, v.get("output_dir", ""))

        saved_sectors = set(v.get("sectors", []))
        for sector, var in self._sector_vars.items():
            var.set(sector in saved_sectors)

        saved_speeds = set(int(s) for s in v.get("speeds", []))
        for speed, var in self._speed_vars.items():
            var.set(speed in saved_speeds)

    def _on_save_settings(self):
        path = save_settings_dialog(self, default_name=f"{MODULE_ID}_settings.json")
        if not path:
            return
        try:
            save_json_settings(path, self._collect_settings())
            self.status_label.configure(text=f"Settings saved → {Path(path).name}")
        except Exception as exc:
            show_error_dialog(self, "Save Failed", f"Could not save settings:\n{exc}")

    def _on_load_settings(self):
        path = load_settings_dialog(self)
        if not path:
            return
        try:
            data = load_json_settings(path)
            self._apply_settings(data)
            self.status_label.configure(text=f"Settings loaded ← {Path(path).name}")
        except Exception as exc:
            show_error_dialog(self, "Load Failed", f"Could not load settings:\n{exc}")

    def _on_reset(self):
        if self._task is not None and self._task.is_running:
            self._task.cancel()
            self.status_label.configure(text="Stopping — will reset when finished...")
            self.after(200, self._reset_when_idle)
            return
        self._reset_ui_to_defaults()

    def _reset_when_idle(self):
        if self._task is not None and self._task.is_running:
            self.after(200, self._reset_when_idle)
            return
        self._reset_ui_to_defaults()

    def _reset_ui_to_defaults(self):
        for entry in (
            self.wind_entry, self.ref_start, self.ref_end,
            self.scan_start, self.scan_end, self.freq_entry,
            self.output_name_entry, self.output_dir_entry,
        ):
            entry.delete(0, "end")

        # Default: all sectors selected
        for var in self._sector_vars.values():
            var.set(True)

        # Default: speed classes 0-9 selected
        for i, var in self._speed_vars.items():
            var.set(i <= 9)

        clear_log(self.log_textbox)
        if self.progress_bar:
            self.progress_bar.set(0)
        self.eta_label.configure(text="")
        self.status_label.configure(text="Ready")
        set_widgets_state(self._interactive_widgets + [self.run_btn], "normal")
        self._task = None


if __name__ == '__main__':
    apply_theme()
    root = ctk.CTk()
    root.withdraw()
    win = RepPeriodWindow(root)
    win.protocol("WM_DELETE_WINDOW", root.destroy)
    root.mainloop()
