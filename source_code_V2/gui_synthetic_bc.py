# -*- coding: utf-8 -*-
"""
gui_synthetic_bc.py - Synthetic Boundary Conditions GUI window.

Step-by-step workflow:
    0. Pick file type (water level / currents / waves / sediments / salinity)
    1. Choose source: model files (.mdf|.mdw + .bnd + .grd) or an existing
       file of the same type (non-wave only).
    2. Time window (start, end, timestep minutes).
    3. Number of intermediate control points (1-10) + per-point date entries.
    4. Value table (rows = control dates, cols = all endpoints x parameters).
    5. Interpolation mode (linear / constant).
    6. Output filename + wave-only .mdw companion option.
    7. Generate.

Can be run standalone: python gui_synthetic_bc.py
"""

import customtkinter as ctk
from pathlib import Path
from typing import Optional, Dict, List

from gui_theme import (
    apply_theme, COLORS,
    make_button, make_label, make_entry, make_optionmenu,
    make_separator, make_step_header, make_file_row,
    make_checkbox, make_radiobutton, make_scrollable_frame,
)
from gui_utils import (
    create_sub_window_layout, BackgroundTask,
    browse_file, append_log,
    clear_log,
    save_json_settings, load_json_settings,
    save_settings_dialog, load_settings_dialog,
    show_success_in_log, show_error_dialog,
    set_widgets_state,
)

import file_parsers as fp
from backend_synthetic_bc import (
    FILE_TYPES, COLUMN_SPECS, generate_synthetic_bc,
)


MODULE_ID = "synthetic_bc"
SETTINGS_VERSION = 1

FILE_TYPE_LABELS = {
    "Water level (.bct)":  "waterlevel_bct",
    "Currents (.bct)":     "currents_bct",
    "Waves (.bcw)":        "waves_bcw",
    "Sediments (.bcc)":    "sediments_bcc",
    "Salinity (.bcc)":     "salinity_bcc",
    "Wind point time series (.wnd)": "point_wnd",
}
POINT_WND_PSEUDO_BND = "wind_point"
FILE_TYPE_LABEL_LIST = list(FILE_TYPE_LABELS.keys())

SOURCE_MODES = ["Model files (.mdf/.mdw + .bnd + .grd)", "Existing file of same type"]
SRC_MODEL = SOURCE_MODES[0]
SRC_EXISTING = SOURCE_MODES[1]

MAX_MIDPOINTS = 10


class SyntheticBCWindow(ctk.CTkToplevel):
    """Synthetic Boundary Conditions sub-module window."""

    def __init__(self, parent=None):
        super().__init__(parent)
        left, self.log_textbox, self.progress_bar, self.status_label, self.eta_label = \
            create_sub_window_layout(self, title="Synthetic Boundary Conditions")

        self._interactive_widgets = []
        self._task: Optional[BackgroundTask] = None
        self._boundary_names: List[str] = []
        self._fraction_name_entries: List[ctk.CTkEntry] = []
        self._midpoint_date_entries: List[ctk.CTkEntry] = []
        self._value_entries: Dict[tuple, ctk.CTkEntry] = {}
        self._build_workflow(left)
        self._on_file_type_change(self.ftype_menu.get())

    # ─── Layout ─────────────────────────────────────────────────────────
    def _build_workflow(self, parent):
        pad = dict(anchor="w", padx=10)

        # Step 0 — File Type
        make_step_header(parent, 0, "File Type").pack(**pad, pady=(15, 5))
        self.ftype_menu = make_optionmenu(
            parent, FILE_TYPE_LABEL_LIST, command=self._on_file_type_change,
        )
        self.ftype_menu.pack(**pad, pady=3)
        self._interactive_widgets.append(self.ftype_menu)

        make_separator(parent).pack(fill="x", padx=10, pady=(12, 8))

        # Step 1 — Source
        make_step_header(parent, 1, "Source").pack(**pad, pady=(0, 5))
        self.source_var = ctk.StringVar(value=SRC_MODEL)
        self.source_radio_frame = ctk.CTkFrame(parent, fg_color="transparent")
        self.source_radio_frame.pack(fill="x", padx=10, pady=2)
        self.radio_model = make_radiobutton(
            self.source_radio_frame, SRC_MODEL, self.source_var, SRC_MODEL,
            command=self._on_source_change,
        )
        self.radio_model.pack(anchor="w", pady=2)
        self.radio_existing = make_radiobutton(
            self.source_radio_frame, SRC_EXISTING, self.source_var, SRC_EXISTING,
            command=self._on_source_change,
        )
        self.radio_existing.pack(anchor="w", pady=2)

        # Model-files rows
        self.model_frame = ctk.CTkFrame(parent, fg_color="transparent")
        self.model_frame.pack(fill="x", padx=10)
        row_mdf, self.mdf_entry, mdf_btn = make_file_row(
            self.model_frame, ".mdf (or .mdw for waves)",
        )
        mdf_btn.configure(command=lambda: browse_file(
            self.mdf_entry, "Select .mdf / .mdw file",
            [("Delft3D files", "*.mdf *.mdw"), ("All files", "*.*")],
        ))
        row_mdf.pack(fill="x", pady=2)

        row_bnd, self.bnd_entry, bnd_btn = make_file_row(self.model_frame, ".bnd")
        bnd_btn.configure(command=lambda: browse_file(
            self.bnd_entry, "Select .bnd file",
            [("Boundary files", "*.bnd"), ("All files", "*.*")],
        ))
        row_bnd.pack(fill="x", pady=2)

        row_grd, self.grd_entry, grd_btn = make_file_row(self.model_frame, ".grd")
        grd_btn.configure(command=lambda: browse_file(
            self.grd_entry, "Select .grd file",
            [("Grid files", "*.grd"), ("All files", "*.*")],
        ))
        row_grd.pack(fill="x", pady=2)

        # Point-wnd-only: model reference datetime (used to compute the
        # minutes-since-reference in the scientific-notation time column).
        self.wnd_ref_frame = ctk.CTkFrame(parent, fg_color="transparent")
        ref_row = ctk.CTkFrame(self.wnd_ref_frame, fg_color="transparent")
        ref_row.pack(fill="x", padx=10, pady=2)
        make_label(ref_row, "Model reference datetime:",
                   font_key="body_sm").pack(side="left", padx=(0, 5))
        self.wnd_ref_entry = make_entry(ref_row, placeholder="YYYY-MM-DD HH:MM:SS",
                                        width=220)
        self.wnd_ref_entry.pack(side="left")
        self._interactive_widgets.append(self.wnd_ref_entry)

        # Existing-file row
        self.existing_frame = ctk.CTkFrame(parent, fg_color="transparent")
        row_ex, self.existing_entry, ex_btn = make_file_row(
            self.existing_frame, "Existing .bct/.bcc",
        )
        ex_btn.configure(command=lambda: browse_file(
            self.existing_entry, "Select existing BC file",
            [("BC files", "*.bct *.bcc"), ("All files", "*.*")],
        ))
        row_ex.pack(fill="x", pady=2)

        # Load boundaries button
        self.load_row = ctk.CTkFrame(parent, fg_color="transparent")
        load_row = self.load_row
        load_row.pack(fill="x", padx=10, pady=(8, 2))
        self.load_btn = make_button(
            load_row, "Load boundaries", command=self._on_load_source,
            style="secondary", width=200,
        )
        self.load_btn.pack(side="left")
        self.bnd_info = make_label(load_row, "", font_key="body_sm",
                                   color=COLORS["text_dim"])
        self.bnd_info.pack(side="left", padx=(10, 0))
        self._interactive_widgets.extend([
            self.radio_model, self.radio_existing,
            self.mdf_entry, self.bnd_entry, self.grd_entry,
            self.existing_entry, self.load_btn,
        ])

        self._step1_end_sep = make_separator(parent)
        self._step1_end_sep.pack(fill="x", padx=10, pady=(12, 8))

        # Step 2 — Time window
        make_step_header(parent, 2, "Time Window").pack(**pad, pady=(0, 5))
        time_row = ctk.CTkFrame(parent, fg_color="transparent")
        time_row.pack(fill="x", padx=10, pady=2)
        make_label(time_row, "Start:", font_key="body_sm").pack(side="left", padx=(0, 5))
        self.start_entry = make_entry(time_row, placeholder="YYYY-MM-DD HH:MM:SS", width=180)
        self.start_entry.pack(side="left", padx=(0, 15))
        make_label(time_row, "End:", font_key="body_sm").pack(side="left", padx=(0, 5))
        self.end_entry = make_entry(time_row, placeholder="YYYY-MM-DD HH:MM:SS", width=180)
        self.end_entry.pack(side="left")

        ts_row = ctk.CTkFrame(parent, fg_color="transparent")
        ts_row.pack(fill="x", padx=10, pady=2)
        make_label(ts_row, "Timestep (min):", font_key="body_sm").pack(side="left", padx=(0, 5))
        self.ts_entry = make_entry(ts_row, placeholder="60", width=100)
        self.ts_entry.pack(side="left")
        self.ts_entry.insert(0, "60")

        ref_row = ctk.CTkFrame(parent, fg_color="transparent")
        ref_row.pack(fill="x", padx=10, pady=2)
        make_label(ref_row, "Reference date:", font_key="body_sm").pack(side="left", padx=(0, 5))
        self.ref_date_entry = make_entry(ref_row, placeholder="YYYYMMDD", width=120)
        self.ref_date_entry.pack(side="left", padx=(0, 10))
        make_label(ref_row, "(must match model Itdate)",
                   font_key="body_sm", color=COLORS["text_dim"]).pack(side="left")

        self._interactive_widgets.extend([self.start_entry, self.end_entry,
                                          self.ts_entry, self.ref_date_entry])

        make_separator(parent).pack(fill="x", padx=10, pady=(12, 8))

        # Step 3 — Mid points
        make_step_header(parent, 3, "Mid Points").pack(**pad, pady=(0, 5))
        mid_row = ctk.CTkFrame(parent, fg_color="transparent")
        mid_row.pack(fill="x", padx=10, pady=2)
        make_label(mid_row, "Number of mid points:", font_key="body_sm").pack(side="left", padx=(0, 5))
        self.nmid_menu = make_optionmenu(
            mid_row, [str(i) for i in range(0, MAX_MIDPOINTS + 1)],
            command=self._on_nmid_change, width=80,
        )
        self.nmid_menu.pack(side="left")
        self.nmid_menu.set("0")
        self.mid_frame = ctk.CTkFrame(parent, fg_color="transparent")
        self._interactive_widgets.append(self.nmid_menu)

        self._step3_end_sep = make_separator(parent)
        self._step3_end_sep.pack(fill="x", padx=10, pady=(12, 8))

        # Step 3b — Sediments fractions (only visible for sediments)
        self.fractions_frame = ctk.CTkFrame(parent, fg_color="transparent")
        make_separator(self.fractions_frame).pack(fill="x", padx=10, pady=(0, 8))
        make_step_header(self.fractions_frame, "3b", "Sediment Fractions").pack(**pad, pady=(0, 5))
        frac_count_row = ctk.CTkFrame(self.fractions_frame, fg_color="transparent")
        frac_count_row.pack(fill="x", padx=10, pady=2)
        make_label(frac_count_row, "Number of fractions (1-4):", font_key="body_sm").pack(side="left", padx=(0, 5))
        self.nfrac_menu = make_optionmenu(
            frac_count_row, ["1", "2", "3", "4"], command=self._on_nfrac_change, width=60,
        )
        self.nfrac_menu.pack(side="left")
        self.frac_entries_frame = ctk.CTkFrame(self.fractions_frame, fg_color="transparent")
        self.frac_entries_frame.pack(fill="x", padx=10, pady=2)
        self._interactive_widgets.append(self.nfrac_menu)
        # (packed/unpacked in _on_file_type_change)

        # Step 4 — Values
        make_step_header(parent, 4, "Values").pack(**pad, pady=(0, 5))
        make_label(parent, "Values at each control date per boundary endpoint:",
                   font_key="body_sm", color=COLORS["text_dim"]).pack(**pad)
        self.value_placeholder = make_label(
            parent,
            "Click 'Load boundaries' above to populate the value table.",
            font_key="body_sm", color=COLORS["text_dim"],
        )
        self.value_scroll = make_scrollable_frame(parent, label_text="")
        self.value_scroll.configure(height=260)

        self._step4_end_sep = make_separator(parent)
        self._step4_end_sep.pack(fill="x", padx=10, pady=(12, 8))

        # Step 5 — Interpolation
        make_step_header(parent, 5, "Interpolation").pack(**pad, pady=(0, 5))
        self.interp_var = ctk.StringVar(value="linear")
        interp_row = ctk.CTkFrame(parent, fg_color="transparent")
        interp_row.pack(fill="x", padx=10, pady=2)
        make_radiobutton(interp_row, "Linear", self.interp_var, "linear").pack(side="left", padx=(0, 15))
        make_radiobutton(interp_row, "Constant (step-hold)", self.interp_var, "constant").pack(side="left")

        make_separator(parent).pack(fill="x", padx=10, pady=(12, 8))

        # Step 6 — Output
        make_step_header(parent, 6, "Output").pack(**pad, pady=(0, 5))
        out_row = ctk.CTkFrame(parent, fg_color="transparent")
        out_row.pack(fill="x", padx=10, pady=2)
        make_label(out_row, "Filename:", font_key="body_sm").pack(side="left", padx=(0, 5))
        self.out_entry = make_entry(out_row, placeholder="synthetic_bc", width=220)
        self.out_entry.pack(side="left", padx=(0, 5))
        self.out_browse_btn = make_button(out_row, "Save As...",
                                          command=self._on_pick_output_dir,
                                          style="settings", width=100)
        self.out_browse_btn.pack(side="left")
        self.out_dir_entry = make_entry(parent, placeholder="Output directory", width=400)
        self.out_dir_entry.pack(**pad, pady=2)

        self.mdw_companion_var = ctk.BooleanVar(value=True)
        self.mdw_checkbox = make_checkbox(
            parent, "Also write companion .mdw (waves only)",
            variable=self.mdw_companion_var,
        )
        # packed/unpacked by _on_file_type_change
        self._interactive_widgets.extend([self.out_entry, self.out_dir_entry,
                                          self.out_browse_btn])

        make_separator(parent).pack(fill="x", padx=10, pady=(12, 8))

        # Step 7 — Generate
        make_step_header(parent, 7, "Generate").pack(**pad, pady=(0, 5))
        self.run_btn = make_button(
            parent, "Generate File", command=self._on_run,
            style="primary", width=280,
        )
        self.run_btn.pack(**pad, pady=(5, 5))

        btn_row = ctk.CTkFrame(parent, fg_color="transparent")
        btn_row.pack(fill="x", padx=10, pady=(5, 20))
        self.save_btn = make_button(btn_row, "Save Settings", self._on_save_settings,
                                    style="settings", width=140)
        self.save_btn.pack(side="left", padx=(0, 5))
        self.load_set_btn = make_button(btn_row, "Load Settings", self._on_load_settings,
                                        style="settings", width=140)
        self.load_set_btn.pack(side="left", padx=(0, 5))
        self.reset_btn = make_button(btn_row, "Reset / Stop", self._on_reset,
                                     style="reset", width=140)
        self.reset_btn.pack(side="left")
        self._interactive_widgets.extend([self.save_btn, self.load_set_btn])

        # Initialise dynamic regions
        self._on_nmid_change("0")
        self._on_nfrac_change("1")

    # ─── Dynamic UI handlers ────────────────────────────────────────────
    def _on_file_type_change(self, label: str):
        file_type = FILE_TYPE_LABELS[label]
        is_wave = (file_type == "waves_bcw")
        is_sediments = (file_type == "sediments_bcc")
        is_point_wnd = (file_type == "point_wnd")

        if is_wave:
            self.radio_existing.configure(state="disabled")
            self.source_var.set(SRC_MODEL)
        else:
            self.radio_existing.configure(state="normal")

        if is_wave:
            try:
                self.mdw_checkbox.pack(anchor="w", padx=10, pady=3)
            except Exception:
                pass
        else:
            try:
                self.mdw_checkbox.pack_forget()
            except Exception:
                pass

        if is_sediments:
            self.fractions_frame.pack(fill="x", padx=0, pady=(0, 10),
                                       before=self._step3_end_sep)
        else:
            self.fractions_frame.pack_forget()

        if is_point_wnd:
            # No boundary source files needed — single-point time series.
            # But user must supply the model reference datetime so the
            # sci-notation time column is relative to the correct zero.
            self.source_radio_frame.pack_forget()
            self.model_frame.pack_forget()
            self.existing_frame.pack_forget()
            self.load_row.pack_forget()
            self.wnd_ref_frame.pack(fill="x", padx=10, pady=2,
                                    before=self._step1_end_sep)
            self._boundary_names = [POINT_WND_PSEUDO_BND]
        else:
            self.wnd_ref_frame.pack_forget()
            self.source_radio_frame.pack(fill="x", padx=10, pady=2,
                                         before=self.model_frame)
            self.load_row.pack(fill="x", padx=10, pady=(8, 2))
            self._boundary_names = []
            self._on_source_change()

        self._rebuild_value_grid()

    def _on_source_change(self):
        if self.source_var.get() == SRC_MODEL:
            self.model_frame.pack(fill="x", padx=10)
            self.existing_frame.pack_forget()
        else:
            self.existing_frame.pack(fill="x", padx=10)
            self.model_frame.pack_forget()

    def _on_nmid_change(self, value: str):
        for w in self.mid_frame.winfo_children():
            w.destroy()
        self._midpoint_date_entries = []
        try:
            n = int(value)
        except ValueError:
            n = 0
        if n > 0:
            self.mid_frame.pack(fill="x", padx=10, pady=2,
                                before=self._step3_end_sep)
        else:
            self.mid_frame.pack_forget()
        for i in range(n):
            row = ctk.CTkFrame(self.mid_frame, fg_color="transparent")
            row.pack(fill="x", pady=2)
            make_label(row, f"Mid point {i + 1} date:",
                       font_key="body_sm").pack(side="left", padx=(0, 5))
            entry = make_entry(row, placeholder="YYYY-MM-DD HH:MM:SS", width=220)
            entry.pack(side="left")
            self._midpoint_date_entries.append(entry)
        self._rebuild_value_grid()

    def _on_nfrac_change(self, value: str):
        for w in self.frac_entries_frame.winfo_children():
            w.destroy()
        self._fraction_name_entries = []
        try:
            n = int(value)
        except ValueError:
            n = 1
        for i in range(n):
            row = ctk.CTkFrame(self.frac_entries_frame, fg_color="transparent")
            row.pack(fill="x", pady=2)
            make_label(row, f"Fraction {i + 1} name:",
                       font_key="body_sm").pack(side="left", padx=(0, 5))
            e = make_entry(row, placeholder=f"fraction_{i+1}", width=160)
            e.pack(side="left")
            self._fraction_name_entries.append(e)
        self._rebuild_value_grid()

    # ─── Boundary source parsing ────────────────────────────────────────
    def _on_load_source(self):
        """Parse source files to populate the boundary-name list."""
        try:
            file_type = FILE_TYPE_LABELS[self.ftype_menu.get()]
            if self.source_var.get() == SRC_EXISTING:
                path = self.existing_entry.get().strip()
                if not path:
                    raise ValueError("Select an existing file first.")
                if file_type.endswith("_bct"):
                    parsed = fp.parse_bct_sections(path)
                else:
                    parsed = fp.parse_bcc_sections(path)
                self._boundary_names = [s['boundary_name']
                                         for s in parsed['sections']]
                info = (f"{len(self._boundary_names)} section(s) · "
                        f"ref {parsed['reference_date'] or '?'}")
            else:
                bnd_path = self.bnd_entry.get().strip()
                if not bnd_path:
                    raise ValueError("Select a .bnd file first.")
                self._boundary_names = [b['bnd_name']
                                         for b in fp.parse_bnd(bnd_path)]
                info = f"{len(self._boundary_names)} boundaries"
            self.bnd_info.configure(text=info)
            self._rebuild_value_grid()
        except Exception as exc:
            self.bnd_info.configure(text="")
            show_error_dialog(self, "Load Failed", str(exc))

    # ─── Value grid ─────────────────────────────────────────────────────
    def _rebuild_value_grid(self):
        for w in self.value_scroll.winfo_children():
            w.destroy()
        self._value_entries = {}

        if not self._boundary_names:
            self.value_scroll.pack_forget()
            self.value_placeholder.pack(anchor="w", padx=10, pady=(0, 5),
                                        before=self._step4_end_sep)
            return
        self.value_placeholder.pack_forget()
        self.value_scroll.pack(fill="both", expand=False, padx=10, pady=5,
                               before=self._step4_end_sep)

        file_type = FILE_TYPE_LABELS[self.ftype_menu.get()]
        try:
            n_mid = int(self.nmid_menu.get())
        except ValueError:
            n_mid = 0

        # Sediments: repeat per fraction
        if file_type == "sediments_bcc":
            fraction_names = [e.get().strip() or f"fraction_{i + 1}"
                              for i, e in enumerate(self._fraction_name_entries)]
        else:
            fraction_names = [None]

        col_spec = COLUMN_SPECS[file_type]

        # Control date column headers
        date_labels = (["Start"]
                       + [f"Mid {i + 1}" for i in range(n_mid)]
                       + ["End"])
        n_dates = len(date_labels)

        # Vertical layout: each boundary is a card with parameters as rows
        # and control dates as columns.
        grid_row = 0
        for frac in fraction_names:
            for bnd in self._boundary_names:
                # Boundary header
                if bnd == POINT_WND_PSEUDO_BND:
                    header_text = "Single point"
                else:
                    header_text = f"[{frac}] {bnd}" if frac else bnd
                bnd_lbl = make_label(self.value_scroll, header_text,
                                     font_key="body_bold")
                bnd_lbl.grid(row=grid_row, column=0, padx=4, pady=(8, 2),
                             sticky="w")

                # Date column headers on the same row
                for c, dlabel in enumerate(date_labels):
                    dlbl = make_label(self.value_scroll, dlabel,
                                      font_key="body_sm", color=COLORS["text_dim"])
                    dlbl.grid(row=grid_row, column=c + 1, padx=2, pady=(8, 2))
                grid_row += 1

                # One row per parameter
                for key, unit in col_spec:
                    param_lbl = make_label(
                        self.value_scroll, f"{key} {unit}",
                        font_key="body_sm",
                    )
                    param_lbl.grid(row=grid_row, column=0, padx=4, pady=2,
                                   sticky="w")
                    for c in range(n_dates):
                        entry = make_entry(self.value_scroll, placeholder="",
                                           width=80)
                        entry.grid(row=grid_row, column=c + 1, padx=2, pady=2)
                        self._value_entries[(frac, bnd, key, c)] = entry
                    grid_row += 1

    # ─── Output directory ───────────────────────────────────────────────
    def _on_pick_output_dir(self):
        from gui_utils import browse_directory
        browse_directory(self.out_dir_entry)

    # ─── Collect values & run ───────────────────────────────────────────
    def _collect_control_dates(self) -> List[str]:
        start = self.start_entry.get().strip()
        end = self.end_entry.get().strip()
        mids = [e.get().strip() for e in self._midpoint_date_entries]
        return [start] + mids + [end]

    def _collect_value_matrix(self) -> Dict:
        """Walk the value grid and produce the shape expected by the backend."""
        file_type = FILE_TYPES_BY_KEY[FILE_TYPE_LABELS[self.ftype_menu.get()]]
        col_spec = COLUMN_SPECS[file_type]
        try:
            n_mid = int(self.nmid_menu.get())
        except ValueError:
            n_mid = 0
        n_rows = 2 + n_mid

        def _read(frac, bnd, key):
            vals = []
            for r in range(n_rows):
                raw = self._value_entries[(frac, bnd, key, r)].get().strip()
                if not raw:
                    raise ValueError(
                        f"Empty value at row {r + 1} for "
                        f"{bnd}/{key}{(' / ' + frac) if frac else ''}"
                    )
                try:
                    vals.append(float(raw))
                except ValueError:
                    raise ValueError(
                        f"Non-numeric value {raw!r} at row {r + 1} for "
                        f"{bnd}/{key}"
                    )
            return vals

        if file_type == "point_wnd":
            return (
                {key: _read(None, POINT_WND_PSEUDO_BND, key)
                 for key, _ in col_spec},
                None,
            )
        if file_type == "sediments_bcc":
            result = {}
            for frac_entry in self._fraction_name_entries:
                frac = frac_entry.get().strip() or None
            fraction_names = [e.get().strip() or f"fraction_{i + 1}"
                              for i, e in enumerate(self._fraction_name_entries)]
            for frac in fraction_names:
                result[frac] = {}
                for bnd in self._boundary_names:
                    result[frac][bnd] = {
                        key: _read(frac, bnd, key) for key, _ in col_spec
                    }
            return result, fraction_names
        else:
            result = {}
            for bnd in self._boundary_names:
                result[bnd] = {key: _read(None, bnd, key) for key, _ in col_spec}
            return result, None

    def _on_run(self):
        try:
            file_type = FILE_TYPE_LABELS[self.ftype_menu.get()]
            is_point_wnd = (file_type == "point_wnd")
            if not is_point_wnd and not self._boundary_names:
                raise ValueError("Load boundaries first (click 'Load boundaries').")
            source_mode = ("model_files" if self.source_var.get() == SRC_MODEL
                            else "existing_file")

            # --- Validate required fields ---
            missing = []
            if not self.start_entry.get().strip():
                missing.append("Start date")
            if not self.end_entry.get().strip():
                missing.append("End date")
            if not self.ts_entry.get().strip():
                missing.append("Timestep")
            if not self.ref_date_entry.get().strip():
                missing.append("Reference date (YYYYMMDD)")
            for i, e in enumerate(self._midpoint_date_entries):
                if not e.get().strip():
                    missing.append(f"Mid point {i + 1} date")
            if is_point_wnd:
                if not self.wnd_ref_entry.get().strip():
                    missing.append("Model reference datetime")
            elif source_mode == "model_files":
                if not self.mdf_entry.get().strip():
                    missing.append(".mdf / .mdw file path")
                if not self.bnd_entry.get().strip():
                    missing.append(".bnd file path")
                if not self.grd_entry.get().strip():
                    missing.append(".grd file path")
            else:
                if not self.existing_entry.get().strip():
                    missing.append("Existing BC file path")
            if not self.out_dir_entry.get().strip():
                missing.append("Output directory")
            if missing:
                raise ValueError(
                    "The following required fields are empty:\n\n• "
                    + "\n• ".join(missing)
                )

            ts_min = float(self.ts_entry.get().strip())
            control_dates = self._collect_control_dates()
            values, fraction_names = self._collect_value_matrix()

            out_name = self.out_entry.get().strip() or "synthetic_bc"
            out_dir = self.out_dir_entry.get().strip()
            ext = {
                "waterlevel_bct": ".bct", "currents_bct": ".bct",
                "waves_bcw": ".bcw",
                "sediments_bcc": ".bcc", "salinity_bcc": ".bcc",
                "point_wnd": ".wnd",
            }[file_type]
            if not out_name.endswith(ext):
                out_name += ext
            output_path = str(Path(out_dir) / out_name)

            ref_date_str = self.ref_date_entry.get().strip()

            kwargs = dict(
                file_type=file_type,
                source_mode=source_mode,
                start_date=self.start_entry.get().strip(),
                end_date=self.end_entry.get().strip(),
                timestep_min=ts_min,
                control_dates=control_dates,
                values=values,
                interpolation=self.interp_var.get(),
                output_path=output_path,
                fraction_names=fraction_names,
                write_companion_mdw=bool(self.mdw_companion_var.get()),
                user_reference_date=ref_date_str,
            )
            if is_point_wnd:
                kwargs["reference_datetime"] = self.wnd_ref_entry.get().strip()
            elif source_mode == "model_files":
                if file_type == "waves_bcw":
                    kwargs["mdw_path"] = self.mdf_entry.get().strip()
                else:
                    kwargs["mdf_path"] = self.mdf_entry.get().strip()
                kwargs["bnd_path"] = self.bnd_entry.get().strip()
                kwargs["grd_path"] = self.grd_entry.get().strip()
            else:
                kwargs["existing_file_path"] = self.existing_entry.get().strip()
        except Exception as exc:
            show_error_dialog(self, "Input Error", str(exc))
            return

        set_widgets_state(self._interactive_widgets + [self.run_btn], "disabled")
        self.status_label.configure(text=f"Generating {file_type} file...")

        self._task = BackgroundTask(
            target=generate_synthetic_bc,
            kwargs=kwargs,
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
        n = result.get("n_timesteps", "?")
        out_path = result.get("output_path", "")
        mdw = result.get("mdw_path")
        append_log(self.log_textbox, "")
        append_log(self.log_textbox, "─── File preview ───")
        for line in result.get("preview_lines", []):
            append_log(self.log_textbox, line)
        msg = f"Synthetic BC file complete — {n} timesteps · {Path(out_path).name}"
        if mdw:
            msg += f" · companion {Path(mdw).name}"
        show_success_in_log(self.log_textbox, msg)
        self.status_label.configure(text=f"Complete — {n} timesteps")
        set_widgets_state(self._interactive_widgets + [self.run_btn], "normal")

    def _on_error(self, exc):
        self.status_label.configure(text=f"Error: {exc}")
        set_widgets_state(self._interactive_widgets + [self.run_btn], "normal")
        show_error_dialog(self, "Generation Error", str(exc))

    def _on_cancel(self):
        self.status_label.configure(text="Cancelled")
        set_widgets_state(self._interactive_widgets + [self.run_btn], "normal")

    # ─── Settings ───────────────────────────────────────────────────────
    def _collect_settings(self) -> dict:
        return {
            "module": MODULE_ID,
            "version": SETTINGS_VERSION,
            "values": {
                "file_type": self.ftype_menu.get(),
                "source_mode": self.source_var.get(),
                "mdf_path": self.mdf_entry.get().strip(),
                "bnd_path": self.bnd_entry.get().strip(),
                "grd_path": self.grd_entry.get().strip(),
                "existing_path": self.existing_entry.get().strip(),
                "start": self.start_entry.get().strip(),
                "end": self.end_entry.get().strip(),
                "timestep": self.ts_entry.get().strip(),
                "reference_date": self.ref_date_entry.get().strip(),
                "nmid": self.nmid_menu.get(),
                "midpoints": [e.get().strip()
                              for e in self._midpoint_date_entries],
                "nfrac": self.nfrac_menu.get(),
                "fractions": [e.get().strip()
                              for e in self._fraction_name_entries],
                "interpolation": self.interp_var.get(),
                "output_name": self.out_entry.get().strip(),
                "output_dir": self.out_dir_entry.get().strip(),
                "write_companion_mdw": bool(self.mdw_companion_var.get()),
            },
        }

    def _apply_settings(self, settings: dict) -> None:
        if settings.get("module") != MODULE_ID:
            raise ValueError(
                f"Settings file is for '{settings.get('module')}', "
                f"not '{MODULE_ID}'"
            )
        values = settings.get("values", {})

        def _set_entry(entry, v):
            entry.delete(0, "end")
            entry.insert(0, str(v))

        ft = values.get("file_type", FILE_TYPE_LABEL_LIST[0])
        if ft in FILE_TYPE_LABEL_LIST:
            self.ftype_menu.set(ft)
            self._on_file_type_change(ft)

        sm = values.get("source_mode", SRC_MODEL)
        if sm in SOURCE_MODES:
            self.source_var.set(sm)
            self._on_source_change()

        _set_entry(self.mdf_entry, values.get("mdf_path", ""))
        _set_entry(self.bnd_entry, values.get("bnd_path", ""))
        _set_entry(self.grd_entry, values.get("grd_path", ""))
        _set_entry(self.existing_entry, values.get("existing_path", ""))
        _set_entry(self.start_entry, values.get("start", ""))
        _set_entry(self.end_entry, values.get("end", ""))
        _set_entry(self.ts_entry, values.get("timestep", "60"))
        _set_entry(self.ref_date_entry, values.get("reference_date", ""))

        nmid = values.get("nmid", "0")
        self.nmid_menu.set(nmid)
        self._on_nmid_change(nmid)
        for e, v in zip(self._midpoint_date_entries, values.get("midpoints", [])):
            _set_entry(e, v)

        nfrac = values.get("nfrac", "1")
        self.nfrac_menu.set(nfrac)
        self._on_nfrac_change(nfrac)
        for e, v in zip(self._fraction_name_entries, values.get("fractions", [])):
            _set_entry(e, v)

        interp = values.get("interpolation", "linear")
        if interp in ("linear", "constant"):
            self.interp_var.set(interp)

        _set_entry(self.out_entry, values.get("output_name", ""))
        _set_entry(self.out_dir_entry, values.get("output_dir", ""))
        self.mdw_companion_var.set(bool(values.get("write_companion_mdw", True)))

    def _on_save_settings(self):
        path = save_settings_dialog(self, default_name=f"{MODULE_ID}_settings.json")
        if not path:
            return
        try:
            save_json_settings(path, self._collect_settings())
            self.status_label.configure(text=f"Settings saved -> {Path(path).name}")
        except Exception as exc:
            show_error_dialog(self, "Save Failed", f"Could not save:\n{exc}")

    def _on_load_settings(self):
        path = load_settings_dialog(self)
        if not path:
            return
        try:
            data = load_json_settings(path)
            self._apply_settings(data)
            self.status_label.configure(text=f"Settings loaded <- {Path(path).name}")
        except Exception as exc:
            show_error_dialog(self, "Load Failed", f"Could not load:\n{exc}")

    def _on_reset(self):
        if self._task is not None and self._task.is_running:
            self._task.cancel()
            self.status_label.configure(text="Stopping...")
            self.after(200, self._reset_when_idle)
            return
        self._reset_ui_to_defaults()

    def _reset_when_idle(self):
        if self._task is not None and self._task.is_running:
            self.after(200, self._reset_when_idle)
            return
        self._reset_ui_to_defaults()

    def _reset_ui_to_defaults(self):
        self.ftype_menu.set(FILE_TYPE_LABEL_LIST[0])
        self._on_file_type_change(FILE_TYPE_LABEL_LIST[0])
        self.source_var.set(SRC_MODEL)
        self._on_source_change()

        for entry in (self.mdf_entry, self.bnd_entry, self.grd_entry,
                      self.existing_entry, self.start_entry, self.end_entry,
                      self.ref_date_entry, self.out_entry, self.out_dir_entry):
            entry.delete(0, "end")
        self.ts_entry.delete(0, "end")
        self.ts_entry.insert(0, "60")

        self.nmid_menu.set("0")
        self._on_nmid_change("0")
        self.nfrac_menu.set("1")
        self._on_nfrac_change("1")
        self.interp_var.set("linear")
        self.mdw_companion_var.set(True)
        self._boundary_names = []
        self.bnd_info.configure(text="")
        self._rebuild_value_grid()

        clear_log(self.log_textbox)
        if self.progress_bar:
            self.progress_bar.set(0)
        self.eta_label.configure(text="")
        self.status_label.configure(text="Ready")
        set_widgets_state(self._interactive_widgets + [self.run_btn], "normal")
        self._task = None


# Convenience map from internal file_type string back to itself (used by
# _collect_value_matrix for safety).
FILE_TYPES_BY_KEY = {k: k for k in FILE_TYPES}


if __name__ == "__main__":
    apply_theme()
    root = ctk.CTk()
    root.withdraw()
    win = SyntheticBCWindow(root)
    win.protocol("WM_DELETE_WINDOW", root.destroy)
    root.mainloop()
