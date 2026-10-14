# -*- coding: utf-8 -*-
"""
gui_bc_extraction.py — Boundary Condition Extraction GUI window.

Step-by-step workflow:
    1. Model type            (Delft3D-4 / Delft3D-FM)
    2. Extraction type       (Water Level / Currents / Waves / Sediments)
    3. Data source           (preset or Custom)
    4. Data files            (NetCDF — auto-inspected)
    5. User input files      (model definition + boundary, dynamic)
                             For Delft3D-4 + Waves: also a .mdw template
                             For Custom data source: a variable-mapping panel
                             with one dropdown per column of the target
                             output file (driven by bc_file_format_spec.md)
    6. Coordinate system
    7. Output timestep
    8. Output directory
    9. Run

Can be run standalone: python gui_bc_extraction.py
"""

import customtkinter as ctk
from pathlib import Path
from typing import Optional, List, Dict, Any

from gui_theme import (
    apply_theme, COLORS, FONTS,
    make_button, make_label, make_entry, make_optionmenu,
    make_separator, make_step_header, make_file_row, make_status_label,
    make_checkbox,
)
from gui_utils import (
    create_sub_window_layout, BackgroundTask,
    browse_file, browse_files, browse_directory,
    append_log, clear_log,
    save_json_settings, load_json_settings,
    save_settings_dialog, load_settings_dialog,
    show_success_in_log, show_error_dialog, show_warning_continue_dialog,
    set_widgets_state,
)
from presets import BC_PRESET_NAMES, get_bc_preset


MODULE_ID = "bc_extraction"
SETTINGS_VERSION = 4

EXTRACTION_TYPES = [
    "Water Level", "Currents", "Waves", "Sediments", "Salinity",
    "WL Gradient (Neumann)",
]
EXT_TYPE_MAP = {
    "Water Level": "water_level",
    "Currents":    "currents",
    "Waves":       "waves",
    "Sediments":   "sediments",
    "Salinity":    "salinity",
    "WL Gradient (Neumann)": "neumann",
}

WAVE_TIMING_SOURCES = ["Auto from .mdf", "Manual"]
WAVE_TIMING_SOURCE_AUTO = WAVE_TIMING_SOURCES[0]
WAVE_TIMING_SOURCE_MANUAL = WAVE_TIMING_SOURCES[1]

# Neumann probe-cloud methods. 'plane' fits a 3x3 cloud (9 probes) and needs
# >= 3 usable probes; 'centered' is a 2-probe finite difference sitting on the
# section line itself, which keeps coastal support points off the beach.
NEUMANN_METHODS = {
    "Plane fit (3x3 probes)": "plane",
    "Centered (2 probes)":    "centered",
}
NEUMANN_METHOD_DEFAULT = "Plane fit (3x3 probes)"
NEUMANN_METHOD_NOTES = {
    "plane": (
        "9 probes per support point, least-squares plane fit. Robust to parent-model "
        "noise and tolerates a couple of failed probes, but the 6 cross-shore probes "
        "can land on dry cells at coastal points."
    ),
    "centered": (
        "2 probes per support point, both on the section line (no cross-shore offset), "
        "so coastal points stay wet. Simple finite difference — more sensitive to noise, "
        "and needs both probes to pass QC."
    ),
}


def d3d4_wave_manual_timing_enabled(
    model_label: str,
    extract_label: str,
    timing_source_label: str,
) -> bool:
    return (
        model_label == "Delft3D-4"
        and extract_label == "Waves"
        and timing_source_label == WAVE_TIMING_SOURCE_MANUAL
    )

# Custom-mode column specs — keys must match BCPreset field names so the
# backend can read them straight from the preset object after the GUI fills
# them in. Source: bc_file_format_spec.md §1.4, §1.5, §2.3, §3.4.
COLUMN_SPECS: Dict[str, List[Dict[str, str]]] = {
    "water_level": [
        {"label": "Time coordinate",        "key": "time_coord",       "expected_unit": ""},
        {"label": "Latitude coordinate",    "key": "lat_coord",        "expected_unit": "degrees_north"},
        {"label": "Longitude coordinate",   "key": "lon_coord",        "expected_unit": "degrees_east"},
        {"label": "Water level",            "key": "water_level_var",  "expected_unit": "m"},
    ],
    "currents": [
        {"label": "Time coordinate",            "key": "time_coord",      "expected_unit": ""},
        {"label": "Latitude coordinate",        "key": "lat_coord",       "expected_unit": "degrees_north"},
        {"label": "Longitude coordinate",       "key": "lon_coord",       "expected_unit": "degrees_east"},
        {"label": "U velocity (east-west)",     "key": "current_u_var",   "expected_unit": "m/s"},
        {"label": "V velocity (north-south)",   "key": "current_v_var",   "expected_unit": "m/s"},
    ],
    "waves": [
        {"label": "Time coordinate",                  "key": "time_coord",          "expected_unit": ""},
        {"label": "Latitude coordinate",              "key": "lat_coord",           "expected_unit": "degrees_north"},
        {"label": "Longitude coordinate",             "key": "lon_coord",           "expected_unit": "degrees_east"},
        {"label": "Significant wave height (Hs/Hm0)", "key": "wave_hs_var",        "expected_unit": "m"},
        {"label": "Peak period (Tp)",                 "key": "wave_tp_var",        "expected_unit": "s"},
        # Direction handled separately (vector vs degrees radio + dropdowns)
        {"label": "Directional spreading",            "key": "wave_dir_spread_var", "expected_unit": "deg / cosine-power m"},
    ],
    "sediments": [
        {"label": "Time coordinate",          "key": "time_coord",  "expected_unit": ""},
        {"label": "Latitude coordinate",      "key": "lat_coord",   "expected_unit": "degrees_north"},
        {"label": "Longitude coordinate",     "key": "lon_coord",   "expected_unit": "degrees_east"},
        # Fraction rows are added dynamically based on the spinbox value.
    ],
    "salinity": [
        {"label": "Time coordinate",          "key": "time_coord",    "expected_unit": ""},
        {"label": "Latitude coordinate",      "key": "lat_coord",     "expected_unit": "degrees_north"},
        {"label": "Longitude coordinate",     "key": "lon_coord",     "expected_unit": "degrees_east"},
        {"label": "Salinity",                 "key": "salinity_var",  "expected_unit": "ppt"},
    ],
    # Neumann gradient generation samples parent-model water levels.
    "neumann": [
        {"label": "Time coordinate",        "key": "time_coord",       "expected_unit": ""},
        {"label": "Latitude coordinate",    "key": "lat_coord",        "expected_unit": "degrees_north"},
        {"label": "Longitude coordinate",   "key": "lon_coord",        "expected_unit": "degrees_east"},
        {"label": "Water level",            "key": "water_level_var",  "expected_unit": "m"},
    ],
}

WAVE_DIRECTION_SPECS: Dict[str, Dict[str, str]] = {
    "wave_dir_x_var": {
        "label": "Wave dir X-component",
        "expected_unit": "",
    },
    "wave_dir_y_var": {
        "label": "Wave dir Y-component",
        "expected_unit": "",
    },
    "wave_dir_var": {
        "label": "Wave direction (degN)",
        "expected_unit": "degrees",
    },
}


class BCExtractionWindow(ctk.CTkToplevel):
    """Boundary Condition Extraction sub-module window."""

    def __init__(self, parent=None):
        super().__init__(parent)

        left, self.log_textbox, self.progress_bar, self.status_label, self.eta_label = \
            create_sub_window_layout(self, title="Boundary Condition Extraction")

        # State ---------------------------------------------------------------
        self._data_files: List[str] = []
        self._source_groups: List[Dict[str, Any]] = []  # cross-border (files, preset) groups
        self._inspected_vars: List[str] = []      # NetCDF data variable names
        self._inspected_coords: List[str] = []    # NetCDF coordinate names
        self._inspected_var_meta: Dict[str, Dict[str, Any]] = {}
        self._inspected_coord_meta: Dict[str, Dict[str, Any]] = {}
        self._custom_widgets: Dict[str, ctk.CTkOptionMenu] = {}  # key -> dropdown
        self._custom_unit_entries: Dict[str, ctk.CTkEntry] = {}  # key -> unit entry
        self._fraction_rows: List[Dict[str, Any]] = []  # for sediment fractions
        self._interactive_widgets: List[Any] = []
        self._task: Optional[BackgroundTask] = None

        self._build_workflow(left)
        # Initialise dynamic sections
        self._refresh_user_input_files()
        self._refresh_custom_mapping()

    # ─── Build ──────────────────────────────────────────────────────────────

    def _build_workflow(self, parent):
        self._workflow_parent = parent
        pad = dict(anchor="w", padx=10)

        # ── Step 1: Model Type ──────────────────────────────────────────
        make_step_header(parent, 1, "Model Type").pack(**pad, pady=(15, 5))
        self.model_type_var = ctk.StringVar(value="Delft3D-4")
        self.model_type_menu = make_optionmenu(
            parent, ["Delft3D-4", "Delft3D-FM"],
            command=self._on_model_type_change,
        )
        self.model_type_menu.pack(**pad, pady=3)
        self._interactive_widgets.append(self.model_type_menu)

        make_separator(parent).pack(fill="x", padx=10, pady=(12, 8))

        # ── Step 2: Extraction Type ─────────────────────────────────────
        make_step_header(parent, 2, "Extraction Type").pack(**pad, pady=(0, 5))
        self.extract_type_menu = make_optionmenu(
            parent, EXTRACTION_TYPES,
            command=self._on_extract_type_change,
        )
        self.extract_type_menu.pack(**pad, pady=3)
        self._interactive_widgets.append(self.extract_type_menu)

        make_separator(parent).pack(fill="x", padx=10, pady=(12, 8))

        # ── Step 3: Data Source ─────────────────────────────────────────
        make_step_header(parent, 3, "Data Source").pack(**pad, pady=(0, 5))
        self.preset_menu = make_optionmenu(
            parent, BC_PRESET_NAMES,
            command=self._on_preset_change,
        )
        self.preset_menu.pack(**pad, pady=3)
        self._interactive_widgets.append(self.preset_menu)

        self.preset_info = make_label(
            parent, self._get_preset_info(BC_PRESET_NAMES[0]),
            font_key="body_sm", color=COLORS["text_dim"],
        )
        self.preset_info.pack(**pad, pady=(0, 5))

        make_separator(parent).pack(fill="x", padx=10, pady=(12, 8))

        # ── Step 4: Data Files ──────────────────────────────────────────
        make_step_header(parent, 4, "Data Files").pack(**pad, pady=(0, 5))
        row_nc, self.nc_entry, self.nc_btn = make_file_row(parent, "NetCDF file(s)")
        self.nc_btn.configure(command=self._browse_data_files)
        row_nc.pack(fill="x", padx=10, pady=3)
        self.nc_status = make_status_label(parent, "No files selected")
        self.nc_status.pack(**pad, pady=(0, 5))
        self._interactive_widgets.append(self.nc_btn)

        # Cross-border source groups: several (files + preset) regions served
        # per boundary point, in the order the groups are added.
        group_row = ctk.CTkFrame(parent, fg_color="transparent")
        self.add_group_btn = make_button(
            group_row, "+ Add as source group", self._on_add_source_group,
            style="secondary", width=170,
        )
        self.add_group_btn.pack(side="left", padx=(0, 6))
        self.clear_groups_btn = make_button(
            group_row, "Clear groups", self._on_clear_source_groups,
            style="neutral", width=110,
        )
        self.clear_groups_btn.pack(side="left")
        group_row.pack(**pad, pady=(0, 2))
        self.group_status = make_status_label(
            parent, "Cross-border mode off — single data source"
        )
        self.group_status.pack(**pad, pady=(0, 5))
        self._interactive_widgets.extend([self.add_group_btn, self.clear_groups_btn])

        make_separator(parent).pack(fill="x", padx=10, pady=(12, 8))

        # ── Step 5: User Input Files (DYNAMIC) ──────────────────────────
        make_step_header(parent, 5, "User Input Files").pack(**pad, pady=(0, 5))

        # Container for the model-definition file rows
        self.user_files_frame = ctk.CTkFrame(parent, fg_color="transparent")
        self.user_files_frame.pack(fill="x", padx=10)

        # D3D4 file rows (mdf, grd, bnd) — packed/unpacked together
        self.d3d4_frame = ctk.CTkFrame(self.user_files_frame, fg_color="transparent")
        self.row_mdf, self.mdf_entry, self.mdf_btn = make_file_row(self.d3d4_frame, ".mdf file")
        self.mdf_btn.configure(command=lambda: browse_file(
            self.mdf_entry, "Select .mdf file",
            [("MDF files", "*.mdf"), ("All files", "*.*")],
        ))
        self.row_mdf.pack(fill="x", pady=2)
        self.row_grd, self.grd_entry, self.grd_btn = make_file_row(self.d3d4_frame, ".grd file")
        self.grd_btn.configure(command=lambda: browse_file(
            self.grd_entry, "Select .grd file",
            [("Grid files", "*.grd"), ("All files", "*.*")],
        ))
        self.row_grd.pack(fill="x", pady=2)
        self.row_bnd, self.bnd_entry, self.bnd_btn = make_file_row(self.d3d4_frame, ".bnd file")
        self.bnd_btn.configure(command=lambda: browse_file(
            self.bnd_entry, "Select .bnd file",
            [("Boundary files", "*.bnd"), ("All files", "*.*")],
        ))
        self.row_bnd.pack(fill="x", pady=2)
        self._interactive_widgets.extend([self.mdf_btn, self.grd_btn, self.bnd_btn])

        # FM file rows (mdu, pli)
        self.fm_frame = ctk.CTkFrame(self.user_files_frame, fg_color="transparent")
        row_mdu, self.mdu_entry, self.mdu_btn = make_file_row(self.fm_frame, ".mdu file")
        self.mdu_btn.configure(command=lambda: browse_file(
            self.mdu_entry, "Select .mdu file",
            [("MDU files", "*.mdu"), ("All files", "*.*")],
        ))
        row_mdu.pack(fill="x", pady=2)
        row_pli, self.pli_entry, self.pli_btn = make_file_row(self.fm_frame, ".pli file")
        self.pli_btn.configure(command=lambda: browse_file(
            self.pli_entry, "Select .pli file",
            [("Polyline files", "*.pli"), ("All files", "*.*")],
        ))
        row_pli.pack(fill="x", pady=2)
        self._interactive_widgets.extend([self.mdu_btn, self.pli_btn])

        # Wave-only: .mdw template (D3D4 waves only)
        self.mdw_frame = ctk.CTkFrame(self.user_files_frame, fg_color="transparent")
        row_mdw, self.mdw_entry, self.mdw_btn = make_file_row(self.mdw_frame, ".mdw template")
        self.mdw_btn.configure(command=lambda: browse_file(
            self.mdw_entry, "Select .mdw template file",
            [("MDW files", "*.mdw"), ("All files", "*.*")],
        ))
        row_mdw.pack(fill="x", pady=2)
        self._interactive_widgets.append(self.mdw_btn)

        self.wave_timing_frame = ctk.CTkFrame(self.user_files_frame, fg_color="transparent")
        wave_timing_row = ctk.CTkFrame(self.wave_timing_frame, fg_color="transparent")
        wave_timing_row.pack(fill="x", pady=(4, 2))
        make_label(wave_timing_row, "Wave timing source", font_key="body_sm").pack(
            side="left", padx=(0, 8)
        )
        self.wave_timing_source_menu = make_optionmenu(
            wave_timing_row,
            WAVE_TIMING_SOURCES,
            command=self._on_wave_timing_source_change,
            width=180,
        )
        self.wave_timing_source_menu.pack(side="left")
        self.wave_timing_source_menu.set(WAVE_TIMING_SOURCE_AUTO)
        self.wave_timing_note = make_label(
            self.wave_timing_frame,
            "",
            font_key="caption",
            color=COLORS["text_dim"],
            justify="left",
            wraplength=560,
        )
        self.wave_timing_note.pack(anchor="w", pady=(0, 4))
        self.wave_manual_timing_frame = ctk.CTkFrame(
            self.wave_timing_frame,
            fg_color="transparent",
        )

        ref_row = ctk.CTkFrame(self.wave_manual_timing_frame, fg_color="transparent")
        ref_row.pack(fill="x", pady=2)
        make_label(ref_row, "Reference date", font_key="body_sm").pack(
            side="left", padx=(0, 8)
        )
        self.wave_reference_date_entry = make_entry(
            ref_row,
            placeholder="YYYYMMDD",
            width=180,
        )
        self.wave_reference_date_entry.pack(side="left")

        start_row = ctk.CTkFrame(self.wave_manual_timing_frame, fg_color="transparent")
        start_row.pack(fill="x", pady=2)
        make_label(start_row, "Simulation start", font_key="body_sm").pack(
            side="left", padx=(0, 8)
        )
        self.wave_start_datetime_entry = make_entry(
            start_row,
            placeholder="YYYY-MM-DD HH:MM",
            width=180,
        )
        self.wave_start_datetime_entry.pack(side="left")

        end_row = ctk.CTkFrame(self.wave_manual_timing_frame, fg_color="transparent")
        end_row.pack(fill="x", pady=2)
        make_label(end_row, "Simulation end", font_key="body_sm").pack(
            side="left", padx=(0, 8)
        )
        self.wave_end_datetime_entry = make_entry(
            end_row,
            placeholder="YYYY-MM-DD HH:MM",
            width=180,
        )
        self.wave_end_datetime_entry.pack(side="left")
        self._interactive_widgets.extend(
            [
                self.wave_timing_source_menu,
                self.wave_reference_date_entry,
                self.wave_start_datetime_entry,
                self.wave_end_datetime_entry,
            ]
        )

        # ── Neumann gradient options (shown for the Neumann extraction) ──
        self.neumann_frame = ctk.CTkFrame(self.user_files_frame, fg_color="transparent")
        neumann_method_row = ctk.CTkFrame(self.neumann_frame, fg_color="transparent")
        neumann_method_row.pack(fill="x", pady=(4, 2))
        make_label(neumann_method_row, "Gradient method", font_key="body_sm").pack(
            side="left", padx=(0, 8)
        )
        self.neumann_method_menu = make_optionmenu(
            neumann_method_row,
            list(NEUMANN_METHODS.keys()),
            command=self._on_neumann_method_change,
            width=220,
        )
        self.neumann_method_menu.pack(side="left")
        self.neumann_method_menu.set(NEUMANN_METHOD_DEFAULT)

        self.neumann_note = make_label(
            self.neumann_frame,
            "",
            font_key="caption",
            color=COLORS["text_dim"],
            justify="left",
            wraplength=560,
        )
        self.neumann_note.pack(anchor="w", pady=(0, 4))

        neumann_span_row = ctk.CTkFrame(self.neumann_frame, fg_color="transparent")
        neumann_span_row.pack(fill="x", pady=(0, 2))
        make_label(neumann_span_row, "Alongshore span", font_key="body_sm").pack(
            side="left", padx=(0, 6)
        )
        self.neumann_along_entry = make_entry(
            neumann_span_row, placeholder="auto", width=90
        )
        self.neumann_along_entry.pack(side="left", padx=(0, 4))
        make_label(neumann_span_row, "m", font_key="body_sm").pack(
            side="left", padx=(0, 14)
        )
        self.neumann_cross_label = make_label(
            neumann_span_row, "Cross-shore span", font_key="body_sm"
        )
        self.neumann_cross_label.pack(side="left", padx=(0, 6))
        self.neumann_cross_entry = make_entry(
            neumann_span_row, placeholder="auto", width=90
        )
        self.neumann_cross_entry.pack(side="left", padx=(0, 4))
        self.neumann_cross_unit = make_label(neumann_span_row, "m", font_key="body_sm")
        self.neumann_cross_unit.pack(side="left")

        make_label(
            self.neumann_frame,
            "Leave blank for automatic spans, scaled to the parent grid "
            "(4x spacing alongshore, 2x cross-shore).",
            font_key="caption",
            color=COLORS["text_dim"],
            justify="left",
            wraplength=560,
        ).pack(anchor="w", pady=(0, 4))

        self._interactive_widgets.extend(
            [
                self.neumann_method_menu,
                self.neumann_along_entry,
                self.neumann_cross_entry,
            ]
        )

        self.intertidal_var = ctk.BooleanVar(value=False)
        self.intertidal_checkbox = make_checkbox(
            self.user_files_frame,
            "Allow intertidal boundaries (wetting/drying zones)",
            variable=self.intertidal_var,
        )
        self.intertidal_checkbox.pack(anchor="w", pady=(6, 2))

        # ── Custom mapping panel (driven by extraction type + NetCDF inspection)
        self.custom_frame = ctk.CTkFrame(parent, fg_color="transparent")
        self.custom_frame.pack(fill="x", padx=10, pady=(8, 0))
        self.custom_header = make_label(
            self.custom_frame,
            "Custom variable mapping",
            font_key="body_sm",
            color=COLORS["text_dim"],
        )
        self.custom_header.pack(anchor="w", pady=(2, 2))
        self.custom_units_note = make_label(
            self.custom_frame,
            "Variable units are auto-filled from NetCDF metadata when available. "
            "Edit them if the dataset metadata is missing or ambiguous.",
            font_key="caption",
            color=COLORS["text_dim"],
            justify="left",
            wraplength=560,
        )
        self.custom_units_note.pack(anchor="w", pady=(0, 4))
        # Direction-format radio (waves only)
        self.wave_dir_format_var = ctk.StringVar(value="vector")
        self.wave_dir_format_frame = ctk.CTkFrame(self.custom_frame, fg_color="transparent")
        # Container that holds the dropdown rows themselves
        self.custom_rows_frame = ctk.CTkFrame(self.custom_frame, fg_color="transparent")
        self.custom_rows_frame.pack(fill="x")
        # Sediment fraction count spinner
        self.sed_count_frame = ctk.CTkFrame(self.custom_frame, fg_color="transparent")
        make_label(self.sed_count_frame, "Number of fractions:",
                   font_key="body_sm").pack(side="left", padx=(0, 6))
        self.sed_count_entry = make_entry(self.sed_count_frame, placeholder="1", width=60)
        self.sed_count_entry.insert(0, "1")
        self.sed_count_entry.pack(side="left")
        self.sed_count_apply = make_button(
            self.sed_count_frame, "Apply", command=self._refresh_custom_mapping,
            style="settings", width=80,
        )
        self.sed_count_apply.pack(side="left", padx=(6, 0))

        # Re-packing the custom frame later should always keep it inside
        # the Step 5 area, ahead of the following separator and Step 6.
        self.custom_panel_anchor = ctk.CTkFrame(parent, fg_color="transparent", height=1)
        self.custom_panel_anchor.pack(fill="x")

        make_separator(parent).pack(fill="x", padx=10, pady=(12, 8))

        # ── Step 6: Coordinate System ───────────────────────────────────
        make_step_header(parent, 6, "Coordinate System").pack(**pad, pady=(0, 5))
        make_label(parent, "UTM Zone (for Cartesian grids):", font_key="body_sm").pack(**pad)
        coord_frame = ctk.CTkFrame(parent, fg_color="transparent")
        coord_frame.pack(fill="x", padx=10, pady=3)
        self.utm_entry = make_entry(coord_frame, placeholder="32", width=80)
        self.utm_entry.pack(side="left", padx=(0, 10))
        self.utm_entry.insert(0, "32")
        self.hemisphere_menu = make_optionmenu(coord_frame, ["N", "S"], width=60)
        self.hemisphere_menu.pack(side="left")
        self._interactive_widgets.extend([self.utm_entry, self.hemisphere_menu])

        make_separator(parent).pack(fill="x", padx=10, pady=(12, 8))

        # ── Step 7: Timestep ────────────────────────────────────────────
        make_step_header(parent, 7, "Output Timestep").pack(**pad, pady=(0, 5))
        ts_frame = ctk.CTkFrame(parent, fg_color="transparent")
        ts_frame.pack(fill="x", padx=10, pady=3)
        self.timestep_entry = make_entry(ts_frame, placeholder="20", width=80)
        self.timestep_entry.pack(side="left", padx=(0, 5))
        self.timestep_entry.insert(0, "20")
        make_label(ts_frame, "minutes", font_key="body_sm").pack(side="left")
        self._interactive_widgets.append(self.timestep_entry)

        make_separator(parent).pack(fill="x", padx=10, pady=(12, 8))

        # ── Step 8: Output Directory ────────────────────────────────────
        make_step_header(parent, 8, "Output").pack(**pad, pady=(0, 5))
        row_out, self.output_entry, self.output_btn = make_file_row(parent, "Output directory")
        self.output_btn.configure(command=lambda: browse_directory(self.output_entry))
        row_out.pack(fill="x", padx=10, pady=3)
        self._interactive_widgets.append(self.output_btn)

        make_separator(parent).pack(fill="x", padx=10, pady=(12, 8))

        # ── Step 9: Run ─────────────────────────────────────────────────
        make_step_header(parent, 9, "Run").pack(**pad, pady=(0, 5))
        self.run_btn = make_button(
            parent, "Extract Boundary Conditions",
            command=self._on_run, style="primary", width=300,
        )
        self.run_btn.pack(**pad, pady=(5, 5))

        # ── Settings / Reset row ────────────────────────────────────────
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

    # ─── Dynamic refresh of Step 5 ──────────────────────────────────────────

    def _refresh_user_input_files(self, *_):
        """Pack only the file rows needed for the current model x extraction."""
        model = self.model_type_menu.get()
        ext = self.extract_type_menu.get()

        for f in (self.d3d4_frame, self.fm_frame, self.mdw_frame,
                  self.wave_timing_frame, self.neumann_frame):
            f.pack_forget()

        if model == "Delft3D-4":
            if ext == "Waves":
                self.wave_timing_frame.pack(fill="x", pady=(2, 2))
                self._refresh_wave_timing_fields()
            self.d3d4_frame.pack(fill="x", pady=2)
            if ext == "Waves":
                self.mdw_frame.pack(fill="x", pady=2)
            if ext == "WL Gradient (Neumann)":
                self.neumann_frame.pack(fill="x", pady=(2, 2))
                self._refresh_neumann_fields()
        else:
            self.fm_frame.pack(fill="x", pady=2)

        self.row_mdf.pack_forget()
        if model == "Delft3D-4" and not self._uses_manual_wave_timing():
            self.row_mdf.pack(fill="x", pady=2, before=self.row_grd)

    def _refresh_custom_mapping(self, *_):
        """
        Rebuild the variable-mapping dropdowns to match the current extraction
        type. Dropdowns are populated with all variables + coords found by the
        last NetCDF inspection. When the data source is not 'Custom', the
        whole panel is hidden (preset values are used directly by the backend).
        """
        # Clear existing rows + state
        for widget in self.custom_rows_frame.winfo_children():
            widget.destroy()
        for widget in self.wave_dir_format_frame.winfo_children():
            widget.destroy()
        self.wave_dir_format_frame.pack_forget()
        self.sed_count_frame.pack_forget()
        self._custom_widgets = {}
        self._custom_unit_entries = {}
        self._fraction_rows = []

        if self.preset_menu.get() != "Custom":
            self.custom_frame.pack_forget()
            return
        self.custom_frame.pack(
            fill="x",
            padx=10,
            pady=(8, 0),
            before=self.custom_panel_anchor,
        )

        ext_key = EXT_TYPE_MAP[self.extract_type_menu.get()]
        choices = self._dropdown_choices()

        # Static columns from the spec
        for col in COLUMN_SPECS[ext_key]:
            self._custom_widgets[col["key"]] = self._make_mapping_row(
                self.custom_rows_frame,
                col["key"],
                col["label"],
                choices,
                expected_unit=col.get("expected_unit", ""),
            )

        # Wave-specific direction-format radio + dependent dropdowns
        if ext_key == "waves":
            self.wave_dir_format_frame.pack(fill="x", pady=(4, 4))
            make_label(
                self.wave_dir_format_frame, "Direction format:",
                font_key="body_sm",
            ).pack(side="left", padx=(0, 6))
            ctk.CTkRadioButton(
                self.wave_dir_format_frame, text="Vector (x/y)",
                variable=self.wave_dir_format_var, value="vector",
                command=self._refresh_wave_dir_rows,
            ).pack(side="left", padx=4)
            ctk.CTkRadioButton(
                self.wave_dir_format_frame, text="Degrees",
                variable=self.wave_dir_format_var, value="degrees",
                command=self._refresh_wave_dir_rows,
            ).pack(side="left", padx=4)
            self._refresh_wave_dir_rows()

        # Sediment fraction rows
        if ext_key == "sediments":
            self.sed_count_frame.pack(fill="x", pady=(4, 4))
            try:
                n = max(1, int(self.sed_count_entry.get().strip()))
            except ValueError:
                n = 1
            for i in range(n):
                row = ctk.CTkFrame(self.custom_rows_frame, fg_color="transparent")
                row.pack(fill="x", pady=(4, 0))
                make_label(row, f"Fraction {i+1} name:", font_key="body_sm").pack(
                    side="left", padx=(0, 4))
                name_entry = make_entry(row, placeholder=f"frac_{i+1}", width=100)
                name_entry.pack(side="left", padx=(0, 8))
                make_label(row, "var:", font_key="body_sm").pack(side="left")
                conc_menu = make_optionmenu(row, choices or ["—"], width=180)
                conc_menu.pack(side="left", padx=(2, 8))
                make_label(row, "ratio:", font_key="body_sm").pack(side="left")
                ratio_entry = make_entry(row, placeholder="1.0", width=60)
                ratio_entry.insert(0, "1.0")
                ratio_entry.pack(side="left")
                self._fraction_rows.append({
                    "name_entry": name_entry,
                    "conc_menu": conc_menu,
                    "ratio_entry": ratio_entry,
                })

    def _refresh_wave_dir_rows(self):
        """Vector vs degrees toggles which direction dropdowns are shown."""
        # Remove any existing direction-related rows we previously added
        for key in ("wave_dir_var", "wave_dir_x_var", "wave_dir_y_var"):
            if key in self._custom_widgets:
                w = self._custom_widgets.pop(key)
                # Climb to the row frame and destroy it
                row_widget = w.master
                row_widget.destroy()

        choices = self._dropdown_choices()
        if self.wave_dir_format_var.get() == "vector":
            self._custom_widgets["wave_dir_x_var"] = self._make_mapping_row(
                self.custom_rows_frame, "Wave dir X-component", choices,
            )
            self._custom_widgets["wave_dir_y_var"] = self._make_mapping_row(
                self.custom_rows_frame, "Wave dir Y-component", choices,
            )
        else:
            self._custom_widgets["wave_dir_var"] = self._make_mapping_row(
                self.custom_rows_frame, "Wave direction (°N)", choices,
            )

    def _make_mapping_row(self, parent, label_text, choices):
        row = ctk.CTkFrame(parent, fg_color="transparent")
        row.pack(fill="x", pady=2)
        make_label(row, label_text, font_key="body_sm").pack(
            side="left", padx=(0, 6))
        menu = make_optionmenu(row, choices or ["—"], width=240)
        menu.pack(side="right")
        return menu

    def _dropdown_choices(self) -> List[str]:
        """Combined list of variable + coordinate names for dropdowns."""
        items = list(dict.fromkeys(self._inspected_coords + self._inspected_vars))
        return items if items else ["—"]

    # ─── Event handlers ─────────────────────────────────────────────────────

    def _choice_metadata(self, name: str) -> Dict[str, Any]:
        """Return cached metadata for an inspected variable or coordinate."""
        if not name or name in ("-", "â€”"):
            return {}
        if name in self._inspected_coord_meta:
            return self._inspected_coord_meta[name]
        return self._inspected_var_meta.get(name, {})

    def _choice_units(self, name: str) -> str:
        """Return the declared units for an inspected choice, if available."""
        info = self._choice_metadata(name)
        units = info.get("units", "")
        return str(units).strip() if units else ""

    def _set_units_entry_text(self, entry: ctk.CTkEntry, units: str) -> None:
        """Replace a units entry value without changing any other state."""
        entry.delete(0, "end")
        if units:
            entry.insert(0, units)

    def _on_custom_variable_change(self, key: str, value: str) -> None:
        """Auto-fill the units box when a custom mapping dropdown changes."""
        entry = self._custom_unit_entries.get(key)
        if entry is None:
            return
        self._set_units_entry_text(entry, self._choice_units(value))

    def _on_fraction_variable_change(self, row_frame, value: str) -> None:
        """Auto-fill sediment fraction units from the selected source variable."""
        for row in self._fraction_rows:
            if row.get("row_frame") is row_frame:
                self._set_units_entry_text(
                    row["units_entry"],
                    self._choice_units(value),
                )
                return

    def _refresh_custom_mapping(self, *_):
        """
        Rebuild the variable-mapping dropdowns to match the current extraction
        type and attach a source-units text field to each selectable variable.
        """
        for widget in self.custom_rows_frame.winfo_children():
            widget.destroy()
        for widget in self.wave_dir_format_frame.winfo_children():
            widget.destroy()
        self.wave_dir_format_frame.pack_forget()
        self.sed_count_frame.pack_forget()
        self._custom_widgets = {}
        self._custom_unit_entries = {}
        self._fraction_rows = []

        if self.preset_menu.get() != "Custom":
            self.custom_frame.pack_forget()
            return

        self.custom_frame.pack(
            fill="x",
            padx=10,
            pady=(8, 0),
            before=self.custom_panel_anchor,
        )

        ext_key = EXT_TYPE_MAP[self.extract_type_menu.get()]
        choices = self._dropdown_choices()

        for col in COLUMN_SPECS[ext_key]:
            self._custom_widgets[col["key"]] = self._make_mapping_row(
                self.custom_rows_frame,
                col["key"],
                col["label"],
                choices,
                expected_unit=col.get("expected_unit", ""),
            )

        if ext_key == "waves":
            self.wave_dir_format_frame.pack(fill="x", pady=(4, 4))
            make_label(
                self.wave_dir_format_frame,
                "Direction format:",
                font_key="body_sm",
            ).pack(side="left", padx=(0, 6))
            ctk.CTkRadioButton(
                self.wave_dir_format_frame,
                text="Vector (x/y)",
                variable=self.wave_dir_format_var,
                value="vector",
                command=self._refresh_wave_dir_rows,
            ).pack(side="left", padx=4)
            ctk.CTkRadioButton(
                self.wave_dir_format_frame,
                text="Degrees",
                variable=self.wave_dir_format_var,
                value="degrees",
                command=self._refresh_wave_dir_rows,
            ).pack(side="left", padx=4)
            self._refresh_wave_dir_rows()

        if ext_key == "sediments":
            self.sed_count_frame.pack(fill="x", pady=(4, 4))
            try:
                n = max(1, int(self.sed_count_entry.get().strip()))
            except ValueError:
                n = 1

            for i in range(n):
                row = ctk.CTkFrame(self.custom_rows_frame, fg_color="transparent")
                row.pack(fill="x", pady=(4, 0))
                make_label(row, f"Fraction {i+1} name:", font_key="body_sm").pack(
                    side="left", padx=(0, 4)
                )
                name_entry = make_entry(row, placeholder=f"frac_{i+1}", width=100)
                name_entry.pack(side="left", padx=(0, 8))
                make_label(row, "var:", font_key="body_sm").pack(side="left")
                conc_menu = make_optionmenu(
                    row,
                    choices or ["-"],
                    width=180,
                    command=lambda value, frac_row=row: self._on_fraction_variable_change(frac_row, value),
                )
                conc_menu.pack(side="left", padx=(2, 8))
                units_entry = make_entry(row, placeholder="kg/m3", width=110)
                units_entry.pack(side="left", padx=(0, 8))
                make_label(row, "ratio:", font_key="body_sm").pack(side="left")
                ratio_entry = make_entry(row, placeholder="1.0", width=60)
                ratio_entry.insert(0, "1.0")
                ratio_entry.pack(side="left")
                self._fraction_rows.append({
                    "row_frame": row,
                    "name_entry": name_entry,
                    "conc_menu": conc_menu,
                    "units_entry": units_entry,
                    "ratio_entry": ratio_entry,
                })

    def _refresh_wave_dir_rows(self):
        """Vector vs degrees toggles which direction dropdowns are shown."""
        for key in ("wave_dir_var", "wave_dir_x_var", "wave_dir_y_var"):
            if key in self._custom_widgets:
                widget = self._custom_widgets.pop(key)
                self._custom_unit_entries.pop(key, None)
                widget.master.destroy()

        choices = self._dropdown_choices()
        if self.wave_dir_format_var.get() == "vector":
            self._custom_widgets["wave_dir_x_var"] = self._make_mapping_row(
                self.custom_rows_frame,
                "wave_dir_x_var",
                WAVE_DIRECTION_SPECS["wave_dir_x_var"]["label"],
                choices,
                expected_unit=WAVE_DIRECTION_SPECS["wave_dir_x_var"]["expected_unit"],
            )
            self._custom_widgets["wave_dir_y_var"] = self._make_mapping_row(
                self.custom_rows_frame,
                "wave_dir_y_var",
                WAVE_DIRECTION_SPECS["wave_dir_y_var"]["label"],
                choices,
                expected_unit=WAVE_DIRECTION_SPECS["wave_dir_y_var"]["expected_unit"],
            )
        else:
            self._custom_widgets["wave_dir_var"] = self._make_mapping_row(
                self.custom_rows_frame,
                "wave_dir_var",
                WAVE_DIRECTION_SPECS["wave_dir_var"]["label"],
                choices,
                expected_unit=WAVE_DIRECTION_SPECS["wave_dir_var"]["expected_unit"],
            )

    def _make_mapping_row(self, parent, key, label_text, choices, expected_unit=""):
        """Build one custom mapping row with dropdown + auto-filled units."""
        row = ctk.CTkFrame(parent, fg_color="transparent")
        row.pack(fill="x", pady=2)
        row.grid_columnconfigure(1, weight=1)

        make_label(row, label_text, font_key="body_sm").grid(
            row=0, column=0, sticky="w", padx=(0, 8)
        )
        menu = make_optionmenu(
            row,
            choices or ["-"],
            width=240,
            command=lambda value, mapping_key=key: self._on_custom_variable_change(mapping_key, value),
        )
        menu.grid(row=0, column=1, sticky="ew", padx=(0, 8))

        units_entry = make_entry(
            row,
            placeholder=expected_unit or "units",
            width=130,
        )
        units_entry.grid(row=0, column=2, sticky="e")
        self._custom_unit_entries[key] = units_entry
        return menu

    def _dropdown_choices(self) -> List[str]:
        """Combined list of coordinate + variable names with an explicit blank option."""
        items = list(dict.fromkeys(self._inspected_coords + self._inspected_vars))
        return (["-"] + items) if items else ["-"]

    def _uses_manual_wave_timing(self) -> bool:
        return d3d4_wave_manual_timing_enabled(
            self.model_type_menu.get(),
            self.extract_type_menu.get(),
            self.wave_timing_source_menu.get(),
        )

    def _refresh_wave_timing_fields(self) -> None:
        self.wave_manual_timing_frame.pack_forget()
        if self._uses_manual_wave_timing():
            self.wave_timing_note.configure(
                text=(
                    "Manual mode makes .mdf optional for Delft3D-4 waves. "
                    "The .bcw reference date comes from the field below, while "
                    "data extraction still uses the simulation start and end window."
                )
            )
            self.wave_manual_timing_frame.pack(fill="x", pady=(2, 2))
        else:
            self.wave_timing_note.configure(
                text=(
                    "Auto mode reads reference date, simulation start, and end "
                    "from .mdf. The .bcw time column will start at that .mdf "
                    "offset, not always at 0.00."
                )
            )

    def _on_wave_timing_source_change(self, value):
        self._refresh_user_input_files()

    def _neumann_method(self) -> str:
        """Backend method key for the selected Neumann gradient method."""
        return NEUMANN_METHODS.get(self.neumann_method_menu.get(), "plane")

    def _refresh_neumann_fields(self) -> None:
        """Update the method note and hide cross-shore span for 'centered'."""
        method = self._neumann_method()
        self.neumann_note.configure(text=NEUMANN_METHOD_NOTES[method])
        # 'centered' places both probes on the section line, so a cross-shore
        # span is meaningless there (_probe_offsets ignores it).
        for widget in (self.neumann_cross_label, self.neumann_cross_entry,
                       self.neumann_cross_unit):
            if method == "centered":
                widget.pack_forget()
            else:
                widget.pack(side="left", padx=(0, 6) if widget is self.neumann_cross_label else (0, 4))

    def _on_neumann_method_change(self, value):
        self._refresh_neumann_fields()

    def _neumann_spans(self) -> Dict[str, Any]:
        """
        Probe spans in metres. Blank means 'let the backend scale them to the
        parent grid spacing', so the key is omitted rather than sent as None.
        """
        spans: Dict[str, Any] = {}
        fields = [("along_span_m", self.neumann_along_entry, "Alongshore span")]
        if self._neumann_method() != "centered":
            fields.append(("cross_span_m", self.neumann_cross_entry, "Cross-shore span"))
        for key, entry, label in fields:
            raw = entry.get().strip()
            if not raw:
                continue
            try:
                value = float(raw)
            except ValueError:
                raise ValueError(f"{label}: '{raw}' is not a number. Leave blank for auto.")
            if value <= 0.0:
                raise ValueError(f"{label} must be greater than 0 m (got {value:g}).")
            spans[key] = value
        return spans

    def _on_model_type_change(self, value):
        self._refresh_user_input_files()

    def _on_extract_type_change(self, value):
        self._refresh_user_input_files()
        self._refresh_custom_mapping()

    def _on_preset_change(self, value):
        self.preset_info.configure(text=self._get_preset_info(value))
        self._refresh_custom_mapping()

    def _get_preset_info(self, name):
        try:
            p = get_bc_preset(name)
            return f"  {p.description}  |  timestep: {p.default_timestep_minutes}min  |  offset: {p.time_offset_hours}h"
        except KeyError:
            return ""

    def _browse_data_files(self):
        files = browse_files(
            self.nc_entry, "Select NetCDF data file(s)",
            [("NetCDF files", "*.nc *.nc4 *.netcdf"), ("All files", "*.*")],
            self.nc_status,
        )
        if files:
            self._data_files = files
            self._inspect_first_file()

    # ─── Cross-border source groups ─────────────────────────────────────────

    def _on_add_source_group(self):
        """Snapshot the current file selection + preset as one source group."""
        if not self._data_files:
            show_error_dialog(
                self, "No data files",
                "Select the NetCDF file(s) for this region first, then add them as a group.",
            )
            return
        preset_name = self.preset_menu.get()
        group: Dict[str, Any] = {
            "data_files": list(self._data_files),
            "preset_name": preset_name,
            "label": f"{len(self._source_groups) + 1}:{Path(self._data_files[0]).stem}",
        }
        if preset_name == "Custom":
            # Snapshot this region's variable mapping now — each Custom group
            # keeps its own copy, applied to a cloned preset in the backend.
            ext_type = EXT_TYPE_MAP[self.extract_type_menu.get()]
            try:
                group["preset_overrides"] = self._collect_custom_overrides(ext_type)
            except ValueError as exc:
                show_error_dialog(
                    self, "Custom mapping incomplete",
                    f"Fill in the Custom variable mapping for this region before "
                    f"adding it as a group.\n\n{exc}",
                )
                return
        self._source_groups.append(group)
        self._data_files = []
        self.nc_entry.delete(0, "end")
        self.nc_status.configure(text="Group added — select the next region's file(s), or Run")
        self._update_group_status()

    def _collect_custom_overrides(self, ext_type: str) -> Dict[str, Any]:
        """
        Build a preset-override dict from the Custom mapping dropdowns
        (same validation as _apply_custom_mapping_to_preset, but without
        touching the shared Custom preset singleton).
        """
        overrides: Dict[str, Any] = {}

        def pick(key: str) -> str:
            menu = self._custom_widgets.get(key)
            if menu is None:
                return ""
            value = menu.get()
            return "" if value in ("", "-", "—", "â€”") else value

        for shared in ('time_coord', 'lat_coord', 'lon_coord'):
            value = pick(shared)
            if not value:
                raise ValueError(f"Custom mapping: '{shared}' is empty.")
            overrides[shared] = value

        if ext_type in ('water_level', 'neumann'):
            wl = pick('water_level_var')
            if not wl:
                raise ValueError("Custom mapping: 'Water level' is empty.")
            overrides['water_level_var'] = wl
        elif ext_type == 'currents':
            u = pick('current_u_var')
            v = pick('current_v_var')
            if not (u and v):
                raise ValueError("Custom mapping: U and V variables are required.")
            overrides['current_u_var'] = u
            overrides['current_v_var'] = v
        elif ext_type == 'waves':
            hs = pick('wave_hs_var')
            tp = pick('wave_tp_var')
            if not (hs and tp):
                raise ValueError("Custom mapping: Hs and Tp are required.")
            overrides['wave_hs_var'] = hs
            overrides['wave_tp_var'] = tp
            fmt = self.wave_dir_format_var.get()
            overrides['wave_dir_format'] = fmt
            if fmt == 'vector':
                dx = pick('wave_dir_x_var')
                dy = pick('wave_dir_y_var')
                if not (dx and dy):
                    raise ValueError("Custom mapping: wave direction X and Y are required.")
                overrides['wave_dir_x_var'] = dx
                overrides['wave_dir_y_var'] = dy
                overrides['wave_dir_var'] = None
            else:
                d = pick('wave_dir_var')
                if not d:
                    raise ValueError("Custom mapping: wave direction (degN) is required.")
                overrides['wave_dir_var'] = d
                overrides['wave_dir_x_var'] = None
                overrides['wave_dir_y_var'] = None
            spread = pick('wave_dir_spread_var')
            overrides['wave_dir_spread_var'] = spread or None
            overrides['wave_dir_spread_use_default'] = not spread
        elif ext_type == 'sediments':
            if not self._fraction_rows:
                raise ValueError("Custom mapping: at least one sediment fraction is required.")
            names, vars_, ratios = [], [], []
            for row in self._fraction_rows:
                name = row['name_entry'].get().strip() or f"frac_{len(names) + 1}"
                var = row['conc_menu'].get()
                if var in ("", "-", "—", "â€”"):
                    raise ValueError(f"Custom mapping: fraction '{name}' has no variable.")
                try:
                    ratio = float(row['ratio_entry'].get().strip() or "1.0")
                except ValueError as exc:
                    raise ValueError(f"Custom mapping: fraction '{name}' has invalid ratio.") from exc
                names.append(name)
                vars_.append(var)
                ratios.append(ratio)
            overrides['sediment_fraction_names'] = names
            overrides['sediment_conc_vars'] = vars_
            overrides['sediment_split_ratios'] = ratios
        elif ext_type == 'salinity':
            sal = pick('salinity_var')
            if not sal:
                raise ValueError("Custom mapping: 'Salinity' variable is empty.")
            overrides['salinity_var'] = sal

        return overrides

    def _on_clear_source_groups(self):
        self._source_groups = []
        self._update_group_status()

    def _update_group_status(self):
        if not self._source_groups:
            self.group_status.configure(text="Cross-border mode off — single data source")
            return
        parts = [
            f"{i + 1}) {g['preset_name']}: {Path(g['data_files'][0]).name}"
            + (f" (+{len(g['data_files']) - 1})" if len(g["data_files"]) > 1 else "")
            for i, g in enumerate(self._source_groups)
        ]
        self.group_status.configure(
            text="Cross-border groups (priority order): " + "   ".join(parts)
        )

    def _inspect_first_file(self):
        """Open the first selected NetCDF and cache its variables/coords + metadata."""
        try:
            from data_inspector import inspect_file, auto_detect_bc_preset
            info = inspect_file(self._data_files[0])
            self._inspected_var_meta = dict(info.get("variables", {}))
            self._inspected_coord_meta = dict(info.get("coords", {}))
            self._inspected_vars = list(self._inspected_var_meta.keys())
            self._inspected_coords = list(self._inspected_coord_meta.keys())
            n = len(self._inspected_vars)
            self.nc_status.configure(
                text=f"✓ {len(self._data_files)} file(s) — {n} variables found"
            )
            # Auto-pick a matching preset, but never override an existing custom mapping
            if self.preset_menu.get() != "Custom":
                detected = auto_detect_bc_preset(info)
                if detected and detected in BC_PRESET_NAMES:
                    self.preset_menu.set(detected)
                    self.preset_info.configure(text=self._get_preset_info(detected))
            self._refresh_custom_mapping()
        except Exception as exc:
            append_log(self.log_textbox, f"[inspect] {exc}\n")
            self.nc_status.configure(
                text=f"✓ {len(self._data_files)} file(s) (inspection failed)"
            )

    # ─── Run ────────────────────────────────────────────────────────────────

    def _inspect_first_file(self):
        """Open the first selected NetCDF and cache its variables/coords + metadata."""
        try:
            from data_inspector import inspect_file, auto_detect_bc_preset

            info = inspect_file(self._data_files[0])
            self._inspected_var_meta = dict(info.get("variables", {}))
            self._inspected_coord_meta = dict(info.get("coords", {}))
            self._inspected_vars = list(self._inspected_var_meta.keys())
            self._inspected_coords = list(self._inspected_coord_meta.keys())
            n = len(self._inspected_vars)
            self.nc_status.configure(
                text=f"âœ“ {len(self._data_files)} file(s) â€” {n} variables found"
            )
            if self.preset_menu.get() != "Custom":
                detected = auto_detect_bc_preset(info)
                if detected and detected in BC_PRESET_NAMES:
                    self.preset_menu.set(detected)
                    self.preset_info.configure(text=self._get_preset_info(detected))
            self._refresh_custom_mapping()
        except Exception as exc:
            self._inspected_vars = []
            self._inspected_coords = []
            self._inspected_var_meta = {}
            self._inspected_coord_meta = {}
            append_log(self.log_textbox, f"[inspect] {exc}\n")
            self.nc_status.configure(
                text=f"âœ“ {len(self._data_files)} file(s) (inspection failed)"
            )

    def _on_run(self):
        model_type = 'delft3d4' if self.model_type_menu.get() == "Delft3D-4" else 'delft3dfm'
        ext_label = self.extract_type_menu.get()
        ext_type = EXT_TYPE_MAP[ext_label]
        wave_timing_source = self.wave_timing_source_menu.get()
        manual_wave_timing = d3d4_wave_manual_timing_enabled(
            self.model_type_menu.get(),
            ext_label,
            wave_timing_source,
        )
        wave_reference_date = ""
        wave_start_datetime = ""
        wave_end_datetime = ""

        # Validate model files
        if model_type == 'delft3d4':
            model_def = self.mdf_entry.get().strip()
            grid = self.grd_entry.get().strip()
            boundary = self.bnd_entry.get().strip()
            if manual_wave_timing:
                if not all([grid, boundary]):
                    show_error_dialog(
                        self,
                        "Missing files",
                        "Manual Delft3D-4 wave timing requires .grd and .bnd files.",
                    )
                    return
                wave_reference_date = self.wave_reference_date_entry.get().strip()
                wave_start_datetime = self.wave_start_datetime_entry.get().strip()
                wave_end_datetime = self.wave_end_datetime_entry.get().strip()
                if not all([wave_reference_date, wave_start_datetime, wave_end_datetime]):
                    show_error_dialog(
                        self,
                        "Missing wave timing",
                        "Please fill in the wave reference date, simulation start, and simulation end.",
                    )
                    return
            elif not all([model_def, grid, boundary]):
                show_error_dialog(
                    self,
                    "Missing files",
                    "Please select .mdf, .grd, and .bnd files.",
                )
                return
            mdw_template = self.mdw_entry.get().strip() if ext_type == 'waves' else ""
            if ext_type == 'waves' and not mdw_template:
                show_error_dialog(
                    self, "Missing .mdw template",
                    "Wave extraction requires an existing .mdw file as a template.")
                return
        else:
            model_def = self.mdu_entry.get().strip()
            boundary = self.pli_entry.get().strip()
            grid = None
            mdw_template = ""
            if not all([model_def, boundary]):
                show_error_dialog(self, "Missing files",
                                  "Please select .mdu and .pli files.")
                return

        # Pending file selection while groups exist counts as the last group.
        if self._source_groups and self._data_files:
            self._on_add_source_group()
        source_groups = [dict(g) for g in self._source_groups]

        if not self._data_files and not source_groups:
            show_error_dialog(self, "Missing data",
                              "Please select NetCDF data file(s).")
            return

        # ── Pre-run data-source validation ──────────────────────────────
        preset_name = self.preset_menu.get()
        if source_groups:
            for group in source_groups:
                if group["preset_name"] == 'Custom':
                    continue  # mapping validated when the group was added
                if not self._validate_data_source(
                    group["preset_name"], ext_type, files=group["data_files"]
                ):
                    return
        else:
            if preset_name != 'Custom' and not self._validate_data_source(
                preset_name, ext_type
            ):
                return

            # Apply Custom mapping → preset object
            if preset_name == 'Custom':
                try:
                    self._apply_custom_mapping_to_preset(ext_type)
                except ValueError as exc:
                    show_error_dialog(self, "Custom mapping incomplete", str(exc))
                    return

        try:
            utm_zone = int(self.utm_entry.get().strip())
        except ValueError:
            utm_zone = 32
        try:
            timestep = float(self.timestep_entry.get().strip())
        except ValueError:
            timestep = 20.0

        output_dir = self.output_entry.get().strip() or '.'
        hemisphere = self.hemisphere_menu.get()

        # Disable controls
        set_widgets_state(self._interactive_widgets + [self.run_btn], "disabled")
        self.status_label.configure(text=f"Running {ext_label.lower()} extraction...")

        # Pick backend function
        import backend_bc_extractor as bce
        bce.set_intertidal_tolerance(self.intertidal_var.get())
        kwargs: Dict[str, Any] = dict(
            model_type=model_type,
            model_def_path=model_def,
            boundary_path=boundary,
            grid_path=grid,
            data_files=self._data_files or (
                source_groups[0]["data_files"] if source_groups else []
            ),
            preset_name=preset_name if preset_name != 'Custom' or not source_groups
            else source_groups[0]["preset_name"],
            timestep=timestep,
            utm_zone=utm_zone,
            hemisphere=hemisphere,
            output_dir=output_dir,
        )
        if source_groups:
            kwargs["source_groups"] = source_groups

        if ext_type == 'water_level':
            target_fn = bce.extract_water_level
        elif ext_type == 'currents':
            target_fn = bce.extract_currents
        elif ext_type == 'waves':
            target_fn = bce.extract_waves
            kwargs['mdw_template_path'] = mdw_template or None
            if model_type == 'delft3d4':
                kwargs['wave_timing_source'] = wave_timing_source
                if manual_wave_timing:
                    kwargs['wave_reference_date'] = wave_reference_date
                    kwargs['wave_start_datetime'] = wave_start_datetime
                    kwargs['wave_end_datetime'] = wave_end_datetime
        elif ext_type == 'sediments':
            target_fn = bce.extract_sediments
        elif ext_type == 'salinity':
            target_fn = bce.extract_salinity
        elif ext_type == 'neumann':
            if model_type != 'delft3d4':
                show_error_dialog(
                    self, "Delft3D-4 only",
                    "Neumann gradient generation currently supports Delft3D-4 "
                    "(.bnd type 'N' + .bct tables).",
                )
                set_widgets_state(self._interactive_widgets + [self.run_btn], "normal")
                return
            import backend_neumann
            target_fn = backend_neumann.generate_neumann_boundaries
            kwargs['method'] = self._neumann_method()
            try:
                spans = self._neumann_spans()
            except ValueError as exc:
                show_error_dialog(self, "Invalid probe span", str(exc))
                set_widgets_state(self._interactive_widgets + [self.run_btn], "normal")
                return
            kwargs.update(spans)
        else:
            show_error_dialog(self, "Internal error",
                              f"Unknown extraction type: {ext_type}")
            set_widgets_state(self._interactive_widgets + [self.run_btn], "normal")
            return

        self._task = BackgroundTask(
            target=target_fn,
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

    def _apply_custom_mapping_to_preset(self, ext_type: str) -> None:
        """
        Push values from the custom dropdowns into the singleton Custom preset
        so the backend reads them transparently. Raises ValueError if any
        required field is empty.
        """
        preset = get_bc_preset('Custom')

        def pick(key: str) -> str:
            menu = self._custom_widgets.get(key)
            if menu is None:
                return ""
            value = menu.get()
            return "" if value in ("", "—") else value

        # Always-shared fields
        for shared in ('time_coord', 'lat_coord', 'lon_coord'):
            v = pick(shared)
            if not v:
                raise ValueError(f"Custom mapping: '{shared}' is empty.")
            setattr(preset, shared, v)

        if ext_type == 'water_level':
            wl = pick('water_level_var')
            if not wl:
                raise ValueError("Custom mapping: 'Water level' is empty.")
            preset.water_level_var = wl
        elif ext_type == 'currents':
            u = pick('current_u_var')
            v = pick('current_v_var')
            if not (u and v):
                raise ValueError("Custom mapping: U and V variables are required.")
            preset.current_u_var = u
            preset.current_v_var = v
        elif ext_type == 'waves':
            hs = pick('wave_hs_var')
            tp = pick('wave_tp_var')
            if not (hs and tp):
                raise ValueError("Custom mapping: Hs and Tp are required.")
            preset.wave_hs_var = hs
            preset.wave_tp_var = tp
            fmt = self.wave_dir_format_var.get()
            preset.wave_dir_format = fmt
            if fmt == 'vector':
                dx = pick('wave_dir_x_var')
                dy = pick('wave_dir_y_var')
                if not (dx and dy):
                    raise ValueError("Custom mapping: wave direction X and Y are required.")
                preset.wave_dir_x_var = dx
                preset.wave_dir_y_var = dy
                preset.wave_dir_var = None
            else:
                d = pick('wave_dir_var')
                if not d:
                    raise ValueError("Custom mapping: wave direction (°N) is required.")
                preset.wave_dir_var = d
                preset.wave_dir_x_var = None
                preset.wave_dir_y_var = None
            spread = pick('wave_dir_spread_var')
            preset.wave_dir_spread_var = spread or None
            preset.wave_dir_spread_use_default = not spread
        elif ext_type == 'sediments':
            if not self._fraction_rows:
                raise ValueError("Custom mapping: at least one sediment fraction is required.")
            names, vars_, ratios = [], [], []
            for row in self._fraction_rows:
                name = row['name_entry'].get().strip() or f"frac_{len(names)+1}"
                var = row['conc_menu'].get()
                if var in ("", "—"):
                    raise ValueError(f"Custom mapping: fraction '{name}' has no variable.")
                try:
                    ratio = float(row['ratio_entry'].get().strip() or "1.0")
                except ValueError:
                    raise ValueError(f"Custom mapping: fraction '{name}' has invalid ratio.")
                names.append(name)
                vars_.append(var)
                ratios.append(ratio)
            preset.sediment_fraction_names = names
            preset.sediment_conc_vars = vars_
            preset.sediment_split_ratios = ratios
        elif ext_type == 'salinity':
            sal = pick('salinity_var')
            if not sal:
                raise ValueError("Custom mapping: 'Salinity' variable is empty.")
            preset.salinity_var = sal

    # ─── Task callbacks ─────────────────────────────────────────────────────

    def _apply_custom_mapping_to_preset(self, ext_type: str) -> None:
        """
        Push values from the custom dropdowns into the singleton Custom preset
        so the backend reads them transparently. Raises ValueError if any
        required field is empty.
        """
        preset = get_bc_preset('Custom')

        def pick(key: str) -> str:
            menu = self._custom_widgets.get(key)
            if menu is None:
                return ""
            value = menu.get()
            return "" if value in ("", "-", "â€”") else value

        for shared in ('time_coord', 'lat_coord', 'lon_coord'):
            value = pick(shared)
            if not value:
                raise ValueError(f"Custom mapping: '{shared}' is empty.")
            setattr(preset, shared, value)

        if ext_type in ('water_level', 'neumann'):
            wl = pick('water_level_var')
            if not wl:
                raise ValueError("Custom mapping: 'Water level' is empty.")
            preset.water_level_var = wl
        elif ext_type == 'currents':
            u = pick('current_u_var')
            v = pick('current_v_var')
            if not (u and v):
                raise ValueError("Custom mapping: U and V variables are required.")
            preset.current_u_var = u
            preset.current_v_var = v
        elif ext_type == 'waves':
            hs = pick('wave_hs_var')
            tp = pick('wave_tp_var')
            if not (hs and tp):
                raise ValueError("Custom mapping: Hs and Tp are required.")
            preset.wave_hs_var = hs
            preset.wave_tp_var = tp
            fmt = self.wave_dir_format_var.get()
            preset.wave_dir_format = fmt
            if fmt == 'vector':
                dx = pick('wave_dir_x_var')
                dy = pick('wave_dir_y_var')
                if not (dx and dy):
                    raise ValueError("Custom mapping: wave direction X and Y are required.")
                preset.wave_dir_x_var = dx
                preset.wave_dir_y_var = dy
                preset.wave_dir_var = None
            else:
                d = pick('wave_dir_var')
                if not d:
                    raise ValueError("Custom mapping: wave direction (degN) is required.")
                preset.wave_dir_var = d
                preset.wave_dir_x_var = None
                preset.wave_dir_y_var = None
            spread = pick('wave_dir_spread_var')
            preset.wave_dir_spread_var = spread or None
            preset.wave_dir_spread_use_default = not spread
        elif ext_type == 'sediments':
            if not self._fraction_rows:
                raise ValueError("Custom mapping: at least one sediment fraction is required.")
            names, vars_, ratios = [], [], []
            for row in self._fraction_rows:
                name = row['name_entry'].get().strip() or f"frac_{len(names)+1}"
                var = row['conc_menu'].get()
                if var in ("", "-", "â€”"):
                    raise ValueError(f"Custom mapping: fraction '{name}' has no variable.")
                try:
                    ratio = float(row['ratio_entry'].get().strip() or "1.0")
                except ValueError as exc:
                    raise ValueError(f"Custom mapping: fraction '{name}' has invalid ratio.") from exc
                names.append(name)
                vars_.append(var)
                ratios.append(ratio)
            preset.sediment_fraction_names = names
            preset.sediment_conc_vars = vars_
            preset.sediment_split_ratios = ratios
        elif ext_type == 'salinity':
            sal = pick('salinity_var')
            if not sal:
                raise ValueError("Custom mapping: 'Salinity' variable is empty.")
            preset.salinity_var = sal

    def _validate_data_source(
        self, preset_name: str, ext_type: str, files: Optional[List[str]] = None
    ) -> bool:
        """
        Inspect data files and warn if any don't match the chosen preset.
        Checks self._data_files unless an explicit file list (e.g. one
        cross-border source group) is given. Returns True if the user wants
        to continue (or no issue found), False if the user cancels.
        """
        from data_inspector import (
            inspect_file, auto_detect_bc_preset,
            validate_variables, get_required_vars_for_extraction,
        )
        has_issue = False

        for fpath in (files if files is not None else self._data_files):
            try:
                info = inspect_file(fpath)
            except Exception:
                continue

            detected = auto_detect_bc_preset(info)
            if detected and detected != preset_name:
                has_issue = True
                break

            try:
                required = get_required_vars_for_extraction(preset_name, ext_type)
                presence = validate_variables(info, required)
                if any(not ok for ok in presence.values()):
                    has_issue = True
                    break
            except Exception:
                pass

        if not has_issue:
            return True

        msg = (
            "The selected files don't appear to match the chosen data source.\n\n"
            "Continue anyway?"
        )
        return show_warning_continue_dialog(self, "Data Source Warning", msg)

    def _on_complete(self, result):
        files = result.get('output_files', [])
        repaired = result.get('repaired_points', [])
        rejected = result.get('rejected_boundaries', [])
        intertidal = result.get('intertidal_points', [])

        msg = (
            f"Extraction complete — {result.get('n_boundaries', '?')} boundaries, "
            f"{result.get('n_timesteps', '?')} timesteps"
        )
        if repaired:
            msg += f"\n\nQC repaired {len(repaired)} point(s):"
            for r in repaired:
                msg += f"\n  - {r}"
        if rejected:
            msg += f"\n\nQC rejected {len(rejected)} boundary section(s):"
            for r in rejected:
                msg += f"\n  - {r}"
        if intertidal:
            msg += f"\n\nIntertidal tolerance accepted {len(intertidal)} point(s):"
            for r in intertidal:
                msg += f"\n  - {r}"
        assignments = result.get('source_assignments', {})
        if assignments:
            counts: Dict[str, int] = {}
            for label in assignments.values():
                counts[label] = counts.get(label, 0) + 1
            msg += "\n\nCross-border source usage (see log for the per-point coverage report):"
            for label, count in counts.items():
                msg += f"\n  - {label}: {count} point(s)"

        show_success_in_log(self.log_textbox, msg)

        status = f"Complete — {', '.join(Path(f).name for f in files)}"
        if repaired or rejected:
            status += f" ({len(repaired)} repaired, {len(rejected)} rejected)"
        self.status_label.configure(text=status)
        set_widgets_state(self._interactive_widgets + [self.run_btn], "normal")

    def _on_error(self, exc):
        self.status_label.configure(text=f"✗ Error: {exc}")
        set_widgets_state(self._interactive_widgets + [self.run_btn], "normal")
        show_error_dialog(self, "Extraction Error", str(exc))

    def _on_cancel(self):
        self.status_label.configure(text="Cancelled")
        set_widgets_state(self._interactive_widgets + [self.run_btn], "normal")

    # ─── Settings save/load/reset ───────────────────────────────────────────

    def _collect_settings(self) -> dict:
        custom_values = {
            key: menu.get() for key, menu in self._custom_widgets.items()
        }
        fractions = [
            {
                "name": r["name_entry"].get().strip(),
                "var":  r["conc_menu"].get(),
                "ratio": r["ratio_entry"].get().strip(),
            }
            for r in self._fraction_rows
        ]
        return {
            "module": MODULE_ID,
            "version": SETTINGS_VERSION,
            "values": {
                "model_type":     self.model_type_menu.get(),
                "extract_type":   self.extract_type_menu.get(),
                "preset":         self.preset_menu.get(),
                "data_files":     list(self._data_files),
                "mdf_file":       self.mdf_entry.get().strip(),
                "grd_file":       self.grd_entry.get().strip(),
                "bnd_file":       self.bnd_entry.get().strip(),
                "mdu_file":       self.mdu_entry.get().strip(),
                "pli_file":       self.pli_entry.get().strip(),
                "mdw_file":       self.mdw_entry.get().strip(),
                "wave_timing_source": self.wave_timing_source_menu.get(),
                "wave_reference_date": self.wave_reference_date_entry.get().strip(),
                "wave_start_datetime": self.wave_start_datetime_entry.get().strip(),
                "wave_end_datetime": self.wave_end_datetime_entry.get().strip(),
                "utm_zone":       self.utm_entry.get().strip(),
                "hemisphere":     self.hemisphere_menu.get(),
                "timestep_minutes": self.timestep_entry.get().strip(),
                "output_dir":     self.output_entry.get().strip(),
                "custom_mapping": custom_values,
                "wave_dir_format": self.wave_dir_format_var.get(),
                "sediment_fractions": fractions,
                "intertidal_tolerance": bool(self.intertidal_var.get()),
                "neumann_method": self.neumann_method_menu.get(),
                "neumann_along_span": self.neumann_along_entry.get().strip(),
                "neumann_cross_span": self.neumann_cross_entry.get().strip(),
            },
        }

    def _apply_settings(self, settings: dict) -> None:
        if settings.get("module") != MODULE_ID:
            raise ValueError(
                f"Settings file is for '{settings.get('module')}', not '{MODULE_ID}'"
            )
        v = settings.get("values", {})

        def _set_entry(entry, value):
            entry.delete(0, "end")
            entry.insert(0, value)

        model_type = v.get("model_type", "Delft3D-4")
        self.model_type_menu.set(model_type)
        ext = v.get("extract_type", "Water Level")
        if ext in EXTRACTION_TYPES:
            self.extract_type_menu.set(ext)

        _set_entry(self.mdf_entry, v.get("mdf_file", ""))
        _set_entry(self.grd_entry, v.get("grd_file", ""))
        _set_entry(self.bnd_entry, v.get("bnd_file", ""))
        _set_entry(self.mdu_entry, v.get("mdu_file", ""))
        _set_entry(self.pli_entry, v.get("pli_file", ""))
        _set_entry(self.mdw_entry, v.get("mdw_file", ""))
        self.wave_timing_source_menu.set(
            v.get("wave_timing_source", WAVE_TIMING_SOURCE_AUTO)
        )
        _set_entry(self.wave_reference_date_entry, v.get("wave_reference_date", ""))
        _set_entry(self.wave_start_datetime_entry, v.get("wave_start_datetime", ""))
        _set_entry(self.wave_end_datetime_entry, v.get("wave_end_datetime", ""))

        preset = v.get("preset", BC_PRESET_NAMES[0])
        if preset in BC_PRESET_NAMES:
            self.preset_menu.set(preset)
            self.preset_info.configure(text=self._get_preset_info(preset))

        self._data_files = list(v.get("data_files", []))
        _set_entry(self.nc_entry, "")
        if self._data_files:
            self.nc_entry.insert(0, f"{len(self._data_files)} file(s) selected")
            self.nc_status.configure(text=f"✓ {len(self._data_files)} files loaded")
            try:
                self._inspect_first_file()
            except Exception:
                pass
        else:
            self.nc_status.configure(text="No files selected")

        _set_entry(self.utm_entry, str(v.get("utm_zone", "32")))
        hemi = v.get("hemisphere", "N")
        if hemi in ("N", "S"):
            self.hemisphere_menu.set(hemi)
        _set_entry(self.timestep_entry, str(v.get("timestep_minutes", "20")))
        _set_entry(self.output_entry, v.get("output_dir", ""))

        self.wave_dir_format_var.set(v.get("wave_dir_format", "vector"))
        self.intertidal_var.set(bool(v.get("intertidal_tolerance", False)))

        # Refresh dynamic sections, then re-apply saved dropdown selections
        self._refresh_user_input_files()
        self._refresh_custom_mapping()
        for key, value in v.get("custom_mapping", {}).items():
            menu = self._custom_widgets.get(key)
            if menu is not None and value:
                try:
                    menu.set(value)
                except Exception:
                    pass
        for i, frac in enumerate(v.get("sediment_fractions", [])):
            if i < len(self._fraction_rows):
                row = self._fraction_rows[i]
                row["name_entry"].delete(0, "end")
                row["name_entry"].insert(0, frac.get("name", ""))
                if frac.get("var"):
                    try:
                        row["conc_menu"].set(frac["var"])
                    except Exception:
                        pass
                row["ratio_entry"].delete(0, "end")
                row["ratio_entry"].insert(0, frac.get("ratio", "1.0"))

    def _collect_settings(self) -> dict:
        custom_values = {}
        for key, menu in self._custom_widgets.items():
            unit_entry = self._custom_unit_entries.get(key)
            custom_values[key] = {
                "variable": menu.get(),
                "units": unit_entry.get().strip() if unit_entry is not None else "",
            }

        fractions = [
            {
                "name": r["name_entry"].get().strip(),
                "var": r["conc_menu"].get(),
                "units": r["units_entry"].get().strip(),
                "ratio": r["ratio_entry"].get().strip(),
            }
            for r in self._fraction_rows
        ]

        return {
            "module": MODULE_ID,
            "version": SETTINGS_VERSION,
            "values": {
                "model_type": self.model_type_menu.get(),
                "extract_type": self.extract_type_menu.get(),
                "preset": self.preset_menu.get(),
                "data_files": list(self._data_files),
                "source_groups": [dict(g) for g in self._source_groups],
                "mdf_file": self.mdf_entry.get().strip(),
                "grd_file": self.grd_entry.get().strip(),
                "bnd_file": self.bnd_entry.get().strip(),
                "mdu_file": self.mdu_entry.get().strip(),
                "pli_file": self.pli_entry.get().strip(),
                "mdw_file": self.mdw_entry.get().strip(),
                "wave_timing_source": self.wave_timing_source_menu.get(),
                "wave_reference_date": self.wave_reference_date_entry.get().strip(),
                "wave_start_datetime": self.wave_start_datetime_entry.get().strip(),
                "wave_end_datetime": self.wave_end_datetime_entry.get().strip(),
                "utm_zone": self.utm_entry.get().strip(),
                "hemisphere": self.hemisphere_menu.get(),
                "timestep_minutes": self.timestep_entry.get().strip(),
                "output_dir": self.output_entry.get().strip(),
                "custom_mapping": custom_values,
                "wave_dir_format": self.wave_dir_format_var.get(),
                "sediment_fractions": fractions,
                "intertidal_tolerance": bool(self.intertidal_var.get()),
                "neumann_method": self.neumann_method_menu.get(),
                "neumann_along_span": self.neumann_along_entry.get().strip(),
                "neumann_cross_span": self.neumann_cross_entry.get().strip(),
            },
        }

    def _apply_settings(self, settings: dict) -> None:
        if settings.get("module") != MODULE_ID:
            raise ValueError(
                f"Settings file is for '{settings.get('module')}', not '{MODULE_ID}'"
            )

        values = settings.get("values", {})

        def _set_entry(entry, value):
            entry.delete(0, "end")
            entry.insert(0, value)

        model_type = values.get("model_type", "Delft3D-4")
        self.model_type_menu.set(model_type)

        extract_type = values.get("extract_type", "Water Level")
        if extract_type in EXTRACTION_TYPES:
            self.extract_type_menu.set(extract_type)

        _set_entry(self.mdf_entry, values.get("mdf_file", ""))
        _set_entry(self.grd_entry, values.get("grd_file", ""))
        _set_entry(self.bnd_entry, values.get("bnd_file", ""))
        _set_entry(self.mdu_entry, values.get("mdu_file", ""))
        _set_entry(self.pli_entry, values.get("pli_file", ""))
        _set_entry(self.mdw_entry, values.get("mdw_file", ""))
        self.wave_timing_source_menu.set(
            values.get("wave_timing_source", WAVE_TIMING_SOURCE_AUTO)
        )
        _set_entry(self.wave_reference_date_entry, values.get("wave_reference_date", ""))
        _set_entry(self.wave_start_datetime_entry, values.get("wave_start_datetime", ""))
        _set_entry(self.wave_end_datetime_entry, values.get("wave_end_datetime", ""))

        preset = values.get("preset", BC_PRESET_NAMES[0])
        if preset in BC_PRESET_NAMES:
            self.preset_menu.set(preset)
            self.preset_info.configure(text=self._get_preset_info(preset))

        self._data_files = list(values.get("data_files", []))
        self._source_groups = [
            dict(g) for g in values.get("source_groups", [])
            if g.get("data_files") and g.get("preset_name")
        ]
        self._update_group_status()
        _set_entry(self.nc_entry, "")
        if self._data_files:
            self.nc_entry.insert(0, f"{len(self._data_files)} file(s) selected")
            self.nc_status.configure(text=f"âœ“ {len(self._data_files)} files loaded")
            try:
                self._inspect_first_file()
            except Exception:
                pass
        else:
            self.nc_status.configure(text="No files selected")
            self._inspected_vars = []
            self._inspected_coords = []
            self._inspected_var_meta = {}
            self._inspected_coord_meta = {}

        _set_entry(self.utm_entry, str(values.get("utm_zone", "32")))
        hemisphere = values.get("hemisphere", "N")
        if hemisphere in ("N", "S"):
            self.hemisphere_menu.set(hemisphere)
        _set_entry(self.timestep_entry, str(values.get("timestep_minutes", "20")))
        _set_entry(self.output_entry, values.get("output_dir", ""))

        self.wave_dir_format_var.set(values.get("wave_dir_format", "vector"))
        self.intertidal_var.set(bool(values.get("intertidal_tolerance", False)))

        neumann_method = values.get("neumann_method", NEUMANN_METHOD_DEFAULT)
        if neumann_method in NEUMANN_METHODS:
            self.neumann_method_menu.set(neumann_method)
        _set_entry(self.neumann_along_entry, str(values.get("neumann_along_span", "")))
        _set_entry(self.neumann_cross_entry, str(values.get("neumann_cross_span", "")))

        self._refresh_user_input_files()
        self._refresh_custom_mapping()

        for key, raw_value in values.get("custom_mapping", {}).items():
            menu = self._custom_widgets.get(key)
            units_entry = self._custom_unit_entries.get(key)
            if menu is None:
                continue

            if isinstance(raw_value, dict):
                variable = raw_value.get("variable", "")
                saved_units = raw_value.get("units", "")
            else:
                variable = raw_value
                saved_units = ""

            if variable:
                try:
                    menu.set(variable)
                except Exception:
                    pass
                self._on_custom_variable_change(key, variable)

            if units_entry is not None and saved_units:
                self._set_units_entry_text(units_entry, str(saved_units))

        for i, frac in enumerate(values.get("sediment_fractions", [])):
            if i >= len(self._fraction_rows):
                break

            row = self._fraction_rows[i]
            row["name_entry"].delete(0, "end")
            row["name_entry"].insert(0, frac.get("name", ""))

            variable = frac.get("var", "")
            if variable:
                try:
                    row["conc_menu"].set(variable)
                except Exception:
                    pass
                self._on_fraction_variable_change(row["row_frame"], variable)

            saved_units = frac.get("units", "")
            if saved_units:
                self._set_units_entry_text(row["units_entry"], str(saved_units))

            row["ratio_entry"].delete(0, "end")
            row["ratio_entry"].insert(0, frac.get("ratio", "1.0"))

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
        self.model_type_menu.set("Delft3D-4")
        self.extract_type_menu.set("Water Level")
        self.preset_menu.set(BC_PRESET_NAMES[0])
        self.preset_info.configure(text=self._get_preset_info(BC_PRESET_NAMES[0]))

        for entry in (
            self.mdf_entry, self.grd_entry, self.bnd_entry,
            self.mdu_entry, self.pli_entry, self.mdw_entry,
            self.wave_reference_date_entry,
            self.wave_start_datetime_entry,
            self.wave_end_datetime_entry,
            self.nc_entry, self.output_entry,
        ):
            entry.delete(0, "end")

        self.utm_entry.delete(0, "end")
        self.utm_entry.insert(0, "32")
        self.hemisphere_menu.set("N")

        self.timestep_entry.delete(0, "end")
        self.timestep_entry.insert(0, "20")
        self.wave_timing_source_menu.set(WAVE_TIMING_SOURCE_AUTO)

        self._data_files = []
        self._source_groups = []
        self._update_group_status()
        self._inspected_vars = []
        self._inspected_coords = []
        self._inspected_var_meta = {}
        self._inspected_coord_meta = {}
        self.nc_status.configure(text="No files selected")

        self.intertidal_var.set(False)

        self.neumann_method_menu.set(NEUMANN_METHOD_DEFAULT)
        self.neumann_along_entry.delete(0, "end")
        self.neumann_cross_entry.delete(0, "end")

        self._refresh_user_input_files()
        self._refresh_custom_mapping()

        clear_log(self.log_textbox)
        if self.progress_bar:
            self.progress_bar.set(0)
        self.eta_label.configure(text="")
        self.status_label.configure(text="Ready")
        set_widgets_state(self._interactive_widgets + [self.run_btn], "normal")
        self._task = None

    def _inspect_first_file(self):
        """Open the first selected NetCDF and cache its variables/coords + metadata."""
        try:
            from data_inspector import inspect_file, auto_detect_bc_preset

            info = inspect_file(self._data_files[0])
            self._inspected_var_meta = dict(info.get("variables", {}))
            self._inspected_coord_meta = dict(info.get("coords", {}))
            self._inspected_vars = list(self._inspected_var_meta.keys())
            self._inspected_coords = list(self._inspected_coord_meta.keys())
            n = len(self._inspected_vars)
            self.nc_status.configure(
                text=f"Loaded {len(self._data_files)} file(s) - {n} variables found"
            )
            if self.preset_menu.get() != "Custom":
                detected = auto_detect_bc_preset(info)
                if detected and detected in BC_PRESET_NAMES:
                    self.preset_menu.set(detected)
                    self.preset_info.configure(text=self._get_preset_info(detected))
            self._refresh_custom_mapping()
        except Exception as exc:
            self._inspected_vars = []
            self._inspected_coords = []
            self._inspected_var_meta = {}
            self._inspected_coord_meta = {}
            append_log(self.log_textbox, f"[inspect] {exc}\n")
            self.nc_status.configure(
                text=f"Loaded {len(self._data_files)} file(s) (inspection failed)"
            )


# ─── Standalone entry point ──────────────────────────────────────────────────
if __name__ == '__main__':
    apply_theme()
    root = ctk.CTk()
    root.withdraw()
    win = BCExtractionWindow(root)
    win.protocol("WM_DELETE_WINDOW", root.destroy)
    root.mainloop()
