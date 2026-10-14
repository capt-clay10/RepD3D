# -*- coding: utf-8 -*-
"""
gui_wind_generator.py - Meteo File Generator GUI window.

Step-by-step workflow:
    0. Select extraction file type (Wind fields / Temperature)
    1. Select source (COSMO / ICON-EU / ICON-Global / ERA5 .nc)
    2. Load data files + auxiliary file
    3. Load target grid (.mat)
    4. Configure grid (CRS, UTM zone, dx/dy)
    5. Set reference time
    6. Set output options
    7. Run generation

Can be run standalone: python gui_wind_generator.py
"""

import customtkinter as ctk
from pathlib import Path
from typing import Any, Dict, List, Optional

from gui_theme import (
    apply_theme, COLORS,
    make_button, make_label, make_entry, make_optionmenu,
    make_separator, make_step_header, make_file_row,
    get_asset_path,
    make_checkbox,
)
from gui_utils import (
    create_sub_window_layout, BackgroundTask,
    browse_file, browse_directory,
    clear_log,
    save_json_settings, load_json_settings,
    save_settings_dialog, load_settings_dialog,
    show_success_in_log, show_error_dialog, show_warning_continue_dialog,
    set_widgets_state,
)
from presets import WIND_PRESET_NAMES, get_wind_preset


MODULE_ID = "wind_generator"
SETTINGS_VERSION = 2
MAX_VALIDATION_POINTS = 10
VALIDATION_POINT_COUNT_OPTIONS = [str(i) for i in range(MAX_VALIDATION_POINTS + 1)]

BUNDLED_AUX_FILES = {
    "COSMO REA6": "cosmo_latlon_rot.grb",
    "ICON-DREAM-EU": "ICON-DREAM-EU_grid.nc",
    "ICON-DREAM-Global": "ICON-DREAM-Global_grid.nc",
}

EXTRACTION_TYPES = ["Wind fields (u, v, p)", "Temperature"]
EXT_TYPE_WIND = "Wind fields (u, v, p)"
EXT_TYPE_TEMP = "Temperature"


class WindGeneratorWindow(ctk.CTkToplevel):
    """Meteo File Generator sub-module window."""

    def __init__(self, parent=None):
        super().__init__(parent)
        left, self.log_textbox, self.progress_bar, self.status_label, self.eta_label = \
            create_sub_window_layout(self, title="Meteo File Generator")

        self._interactive_widgets = []
        self._task: Optional[BackgroundTask] = None
        self._auto_aux_path: Optional[str] = None
        self._validation_rows: List[Dict[str, Any]] = []
        self._build_workflow(left)

    def _build_workflow(self, parent):
        pad = dict(anchor="w", padx=10)

        # Step 0: Extraction File Type
        make_step_header(parent, 0, "Extraction File Type").pack(**pad, pady=(15, 5))
        self.ext_type_menu = make_optionmenu(
            parent, EXTRACTION_TYPES, command=self._on_ext_type_change,
        )
        self.ext_type_menu.pack(**pad, pady=3)
        self.ext_type_info = make_label(
            parent,
            "  Wind fields: .amu, .amv, .amp (+ optional .xwind/.ywind)",
            font_key="body_sm", color=COLORS["text_dim"],
        )
        self.ext_type_info.pack(**pad)
        self._interactive_widgets.append(self.ext_type_menu)

        make_separator(parent).pack(fill="x", padx=10, pady=(12, 8))

        # Step 1: Source
        make_step_header(parent, 1, "Source").pack(**pad, pady=(0, 5))
        self.preset_menu = make_optionmenu(
            parent, WIND_PRESET_NAMES, command=self._on_preset_change,
        )
        self.preset_menu.pack(**pad, pady=3)
        self.preset_info = make_label(
            parent, "", font_key="body_sm", color=COLORS["text_dim"]
        )
        self.preset_info.pack(**pad)
        self._interactive_widgets.append(self.preset_menu)

        make_separator(parent).pack(fill="x", padx=10, pady=(12, 8))

        # Step 2: Data Files
        make_step_header(parent, 2, "Wind Data Files").pack(**pad, pady=(0, 5))
        row_data, self.data_entry, self.data_btn = make_file_row(parent, "Data directory")
        self.data_btn.configure(command=lambda: browse_directory(self.data_entry))
        row_data.pack(fill="x", padx=10, pady=2)

        self.aux_frame = ctk.CTkFrame(parent, fg_color="transparent")
        self.aux_frame.pack(fill="x", padx=10)
        self.row_aux, self.aux_entry, self.aux_btn = make_file_row(
            self.aux_frame, "Source grid / auxiliary file"
        )
        self.aux_btn.configure(command=lambda: browse_file(
            self.aux_entry, "Select source grid / auxiliary file",
            [("GRIB files", "*.grb"), ("NetCDF files", "*.nc"), ("All files", "*.*")],
        ))
        self.row_aux.pack(fill="x", pady=2)
        self.aux_hint = make_label(
            self.aux_frame, "", font_key="caption", color=COLORS["text_dim"]
        )
        self.aux_hint.pack(anchor="w", padx=2, pady=(0, 2))
        self._interactive_widgets.extend([self.data_btn, self.aux_btn])
        self._on_preset_change(self.preset_menu.get())

        make_separator(parent).pack(fill="x", padx=10, pady=(12, 8))

        # Step 3: Target Grid
        make_step_header(parent, 3, "Target Grid").pack(**pad, pady=(0, 5))
        row_grid, self.grid_entry, self.grid_btn = make_file_row(parent, ".mat grid file")
        self.grid_btn.configure(command=lambda: browse_file(
            self.grid_entry, "Select .mat grid file",
            [("MATLAB files", "*.mat"), ("All files", "*.*")],
        ))
        row_grid.pack(fill="x", padx=10, pady=3)
        self._interactive_widgets.append(self.grid_btn)

        make_separator(parent).pack(fill="x", padx=10, pady=(12, 8))

        # Step 4: Grid Configuration
        make_step_header(parent, 4, "Grid Configuration").pack(**pad, pady=(0, 5))
        config_frame = ctk.CTkFrame(parent, fg_color="transparent")
        config_frame.pack(fill="x", padx=10, pady=3)

        make_label(config_frame, "CRS:", font_key="body_sm").pack(side="left", padx=(0, 5))
        self.crs_menu = make_optionmenu(config_frame, ["Cartesian", "Spherical"], width=120)
        self.crs_menu.pack(side="left", padx=(0, 15))

        make_label(config_frame, "UTM Zone:", font_key="body_sm").pack(side="left", padx=(0, 5))
        self.utm_entry = make_entry(config_frame, placeholder="32", width=60)
        self.utm_entry.pack(side="left")
        self.utm_entry.insert(0, "32")

        spacing_frame = ctk.CTkFrame(parent, fg_color="transparent")
        spacing_frame.pack(fill="x", padx=10, pady=3)
        make_label(spacing_frame, "dx:", font_key="body_sm").pack(side="left", padx=(0, 5))
        self.dx_entry = make_entry(spacing_frame, placeholder="6000", width=80)
        self.dx_entry.pack(side="left", padx=(0, 15))
        self.dx_entry.insert(0, "6000")
        make_label(spacing_frame, "dy:", font_key="body_sm").pack(side="left", padx=(0, 5))
        self.dy_entry = make_entry(spacing_frame, placeholder="6000", width=80)
        self.dy_entry.pack(side="left")
        self.dy_entry.insert(0, "6000")
        self._interactive_widgets.extend([self.crs_menu, self.utm_entry, self.dx_entry, self.dy_entry])

        make_separator(parent).pack(fill="x", padx=10, pady=(12, 8))

        # Step 5: Reference Time
        make_step_header(parent, 5, "Reference Time").pack(**pad, pady=(0, 5))
        self.ref_time_entry = make_entry(parent, placeholder="YYYY-MM-DD HH:MM:SS", width=250)
        self.ref_time_entry.pack(**pad, pady=3)
        self._interactive_widgets.append(self.ref_time_entry)

        make_separator(parent).pack(fill="x", padx=10, pady=(12, 8))

        # Step 6: Output
        make_step_header(parent, 6, "Output").pack(**pad, pady=(0, 5))
        name_frame = ctk.CTkFrame(parent, fg_color="transparent")
        name_frame.pack(fill="x", padx=10, pady=3)
        make_label(name_frame, "Filename:", font_key="body_sm").pack(side="left", padx=(0, 5))
        self.name_entry = make_entry(name_frame, placeholder="wind_2015", width=200)
        self.name_entry.pack(side="left")

        row_outdir, self.outdir_entry, self.outdir_btn = make_file_row(parent, "Output directory")
        self.outdir_btn.configure(command=lambda: browse_directory(self.outdir_entry))
        row_outdir.pack(fill="x", padx=10, pady=3)
        self._interactive_widgets.append(self.outdir_btn)

        self.wnd_var = ctk.BooleanVar(value=True)
        self.wnd_checkbox = make_checkbox(
            parent, "Also generate .xwind/.ywind files", variable=self.wnd_var
        )
        self.wnd_checkbox.pack(**pad, pady=3)
        self._interactive_widgets.append(self.name_entry)

        make_separator(parent).pack(fill="x", padx=10, pady=(12, 8))

        # Step 7: Validation Points (optional)
        make_step_header(parent, 7, "Validation Points (optional)").pack(**pad, pady=(0, 5))
        make_label(
            parent, "Extract single-point time series at specified locations",
            font_key="body_sm", color=COLORS["text_dim"],
        ).pack(**pad)

        self._val_count_frame = ctk.CTkFrame(parent, fg_color="transparent")
        self._val_count_frame.pack(fill="x", padx=10, pady=3)
        val_count_frame = self._val_count_frame
        make_label(val_count_frame, "Validation points:", font_key="body_sm").pack(side="left", padx=(0, 5))
        self.val_count_menu = make_optionmenu(
            val_count_frame, VALIDATION_POINT_COUNT_OPTIONS,
            width=80, command=self._on_validation_count_change,
        )
        self.val_count_menu.pack(side="left")
        self.val_count_menu.set("0")
        self._interactive_widgets.append(self.val_count_menu)

        self.val_rows_frame = ctk.CTkFrame(parent, fg_color="transparent")
        self.val_rows_frame.pack(fill="x", padx=10, pady=3)
        for idx in range(MAX_VALIDATION_POINTS):
            row_frame = ctk.CTkFrame(self.val_rows_frame, fg_color="transparent")
            make_label(row_frame, f"Point {idx + 1}:", font_key="body_sm").pack(side="left", padx=(0, 6))
            name_entry = make_entry(row_frame, placeholder="Point name", width=140)
            name_entry.pack(side="left", padx=(0, 6))
            lat_entry = make_entry(row_frame, placeholder="Latitude", width=110)
            lat_entry.pack(side="left", padx=(0, 6))
            lon_entry = make_entry(row_frame, placeholder="Longitude", width=110)
            lon_entry.pack(side="left")
            self._validation_rows.append({
                "frame": row_frame,
                "name_entry": name_entry,
                "lat_entry": lat_entry,
                "lon_entry": lon_entry,
            })
            self._interactive_widgets.extend([name_entry, lat_entry, lon_entry])
        self._on_validation_count_change("0")

        make_separator(parent).pack(fill="x", padx=10, pady=(12, 8))

        # Step 8: Run
        make_step_header(parent, 8, "Run").pack(**pad, pady=(0, 5))
        self.run_btn = make_button(
            parent, "Generate Meteo Files",
            command=self._on_run, style="primary", width=280,
        )
        self.run_btn.pack(**pad, pady=(5, 5))

        # Settings / Reset row
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

    def _on_ext_type_change(self, value):
        """Show/hide wind-only controls based on extraction file type."""
        if value == EXT_TYPE_TEMP:
            self.ext_type_info.configure(
                text="  Temperature: .amt (air temperature in Celsius)"
            )
            try:
                self.wnd_checkbox.pack_forget()
            except Exception:
                pass
            try:
                self.run_btn.configure(text="Generate Temperature Field")
            except Exception:
                pass
        else:
            self.ext_type_info.configure(
                text="  Wind fields: .amu, .amv, .amp (+ optional .xwind/.ywind)"
            )
            try:
                self.wnd_checkbox.pack(anchor="w", padx=10, pady=3)
            except Exception:
                pass
            try:
                self.run_btn.configure(text="Generate Wind Fields")
            except Exception:
                pass

    def _resolve_bundled_aux_path(self, preset_name: str) -> Optional[str]:
        """Resolve a packaged auxiliary/grid file for presets that ship one."""
        filename = BUNDLED_AUX_FILES.get(preset_name)
        if not filename:
            return None
        path = Path(get_asset_path(filename))
        return str(path) if path.exists() else None

    def _set_aux_row_label(self, text: str) -> None:
        """Update the row label text created by make_file_row()."""
        try:
            label_widget = self.row_aux.winfo_children()[0]
            label_widget.configure(text=text)
        except Exception:
            pass

    def _on_preset_change(self, value):
        try:
            p = get_wind_preset(value)
            bundled_aux_path = self._resolve_bundled_aux_path(value)
            self._auto_aux_path = bundled_aux_path

            info_parts = [
                p.description,
                f"grid: {p.grid_type}",
                f"rotation: {'yes' if p.needs_rotation else 'no'}",
            ]
            if bundled_aux_path:
                info_parts.append(f"aux: packaged {Path(bundled_aux_path).name}")

            self.preset_info.configure(text="  " + "  |  ".join(info_parts))

            if not hasattr(self, "aux_frame"):
                return

            if bundled_aux_path:
                self.aux_entry.delete(0, "end")
                self.aux_hint.configure(
                    text=f"Using packaged file: {Path(bundled_aux_path).name}"
                )
                self.aux_frame.pack_forget()
            elif p.needs_aux_file:
                self._set_aux_row_label("Auxiliary file")
                self.aux_hint.configure(text="Select the required auxiliary/grid file.")
                self.aux_frame.pack(fill="x", padx=10)
            else:
                self.aux_entry.delete(0, "end")
                self.aux_hint.configure(text="")
                self.aux_frame.pack_forget()
        except KeyError:
            self._auto_aux_path = None

    def _on_validation_count_change(self, value):
        try:
            count = max(0, min(MAX_VALIDATION_POINTS, int(str(value).strip())))
        except ValueError:
            count = 0
        for idx, row in enumerate(self._validation_rows):
            if idx < count:
                row["frame"].pack(fill="x", pady=2)
            else:
                row["frame"].pack_forget()
        if count > 0:
            self.val_rows_frame.pack(fill="x", padx=10, pady=3,
                                     after=self._val_count_frame)
        else:
            self.val_rows_frame.pack_forget()

    def _get_validation_count(self) -> int:
        try:
            return max(0, min(MAX_VALIDATION_POINTS, int(self.val_count_menu.get().strip())))
        except ValueError:
            return 0

    def _collect_validation_points(self) -> Optional[List[Dict[str, Any]]]:
        count = self._get_validation_count()
        points: List[Dict[str, Any]] = []
        seen_names = set()

        for idx in range(count):
            row = self._validation_rows[idx]
            name = row["name_entry"].get().strip()
            lat_raw = row["lat_entry"].get().strip()
            lon_raw = row["lon_entry"].get().strip()

            if not name:
                show_error_dialog(
                    self, "Invalid validation point",
                    f"Validation point {idx + 1} is missing a point name.",
                )
                return None

            name_key = name.casefold()
            if name_key in seen_names:
                show_error_dialog(
                    self, "Invalid validation point",
                    f"Validation point name '{name}' is duplicated.",
                )
                return None
            seen_names.add(name_key)

            try:
                lat = float(lat_raw)
                lon = float(lon_raw)
            except ValueError:
                show_error_dialog(
                    self, "Invalid validation point",
                    f"Validation point '{name}' must have numeric latitude and longitude.",
                )
                return None

            if not -90.0 <= lat <= 90.0:
                show_error_dialog(
                    self, "Invalid validation point",
                    f"Validation point '{name}' has latitude {lat}, but latitude must be between -90 and 90.",
                )
                return None

            if not -180.0 <= lon <= 180.0:
                show_error_dialog(
                    self, "Invalid validation point",
                    f"Validation point '{name}' has longitude {lon}, but longitude must be between -180 and 180.",
                )
                return None

            points.append({"name": name, "lat": lat, "lon": lon})

        return points

    def _on_run(self):
        data_path = self.data_entry.get().strip()
        grid_path = self.grid_entry.get().strip()
        ref_time = self.ref_time_entry.get().strip()
        ext_type = self.ext_type_menu.get()
        is_temp = (ext_type == EXT_TYPE_TEMP)
        default_name = "temperature_output" if is_temp else "wind_output"
        output_name = self.name_entry.get().strip() or default_name
        preset_name = self.preset_menu.get()
        preset = get_wind_preset(preset_name)

        if not data_path:
            label = "temperature data" if is_temp else "wind data"
            show_error_dialog(self, "Missing input", f"Please select {label} directory.")
            return
        if not grid_path:
            show_error_dialog(self, "Missing input", "Please select target grid .mat file.")
            return
        if not ref_time:
            show_error_dialog(self, "Missing input", "Please enter reference time.")
            return

        if is_temp and not getattr(preset, "t_var", None):
            show_error_dialog(
                self,
                "Temperature not supported",
                f"Preset '{preset_name}' does not declare a temperature variable (t_var). "
                "Choose COSMO REA6, ICON-DREAM-EU, ICON-DREAM-Global, or ERA5 (.nc).",
            )
            return

        aux_path = self._auto_aux_path or (self.aux_entry.get().strip() or None)
        if preset.needs_aux_file and not aux_path:
            show_error_dialog(
                self,
                "Missing auxiliary file",
                f"{preset_name} requires an auxiliary/grid file, but the packaged file was not found.",
            )
            return
        crs_type = self.crs_menu.get()
        try:
            utm_zone = int(self.utm_entry.get().strip())
        except ValueError:
            utm_zone = 32
        try:
            dx = float(self.dx_entry.get().strip())
            dy = float(self.dy_entry.get().strip())
        except ValueError:
            dx, dy = 6000.0, 6000.0

        validation_points = None
        if not is_temp:
            validation_points = self._collect_validation_points()
            if validation_points is None:
                return

        output_dir = self.outdir_entry.get().strip() or str(Path(data_path).parent)

        # ── Pre-run data-source validation ──────────────────────────────
        if not self._validate_wind_source(data_path, preset_name, preset, is_temp):
            return

        set_widgets_state(self._interactive_widgets + [self.run_btn], "disabled")

        import backend_wind_processor as bwp

        if is_temp:
            self.status_label.configure(text="Generating temperature field...")
            target_fn = bwp.generate_temperature_field
            kwargs = dict(
                preset_name=preset_name,
                temperature_data_path=data_path,
                grid_mat_path=grid_path,
                reference_time=ref_time,
                output_dir=output_dir,
                output_name=output_name,
                aux_file_path=aux_path,
                crs_type=crs_type,
                utm_zone=utm_zone,
                grid_dx=dx,
                grid_dy=dy,
            )
        else:
            self.status_label.configure(text="Generating wind fields...")
            target_fn = bwp.generate_wind_fields
            kwargs = dict(
                preset_name=preset_name,
                wind_data_path=data_path,
                grid_mat_path=grid_path,
                reference_time=ref_time,
                output_dir=output_dir,
                output_name=output_name,
                aux_file_path=aux_path,
                crs_type=crs_type,
                utm_zone=utm_zone,
                grid_dx=dx,
                grid_dy=dy,
                write_wnd=self.wnd_var.get(),
                validation_points=validation_points,
            )

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

    def _validate_wind_source(
        self, data_path: str, preset_name: str, preset, is_temp: bool,
    ) -> bool:
        """
        Inspect a sample of actual data files that the backend will use and
        warn if they don't match the chosen wind preset.
        Returns True to proceed, False to cancel.
        """
        from data_inspector import inspect_file, auto_detect_wind_preset

        data_dir = Path(data_path)
        if getattr(preset, "name", preset_name) == 'ERA5':
            if data_dir.is_file():
                if data_dir.suffix.lower() != '.nc':
                    show_error_dialog(
                        self,
                        "ERA5 NetCDF required",
                        "Please select a pre-converted/downloaded ERA5 .nc file. "
                        "GRIB input requires ecCodes and is not used for the packaged app.",
                    )
                    return False
            elif data_dir.is_dir() and not any(p.suffix.lower() == '.nc' for p in data_dir.iterdir()):
                show_error_dialog(
                    self,
                    "ERA5 NetCDF required",
                    "The selected directory contains no ERA5 .nc files. "
                    "Use pre-converted/downloaded ERA5 NetCDF files.",
                )
                return False

        if data_dir.is_file():
            candidates = [data_dir]
        else:
            candidates = self._collect_validation_candidates(
                data_dir, preset_name, preset, is_temp,
            )
        if not candidates:
            return True

        # Only inspect the first file as a representative sample
        fpath = candidates[0]
        try:
            info = inspect_file(str(fpath))
        except Exception:
            return True

        detected = auto_detect_wind_preset(info)
        if detected and detected != preset_name:
            msg = (
                "The selected files don't appear to match the chosen data source.\n\n"
                "Continue anyway?"
            )
            return show_warning_continue_dialog(self, "Data Source Warning", msg)

        return True

    @staticmethod
    def _collect_validation_candidates(
        data_dir: Path, preset_name: str, preset, is_temp: bool,
    ) -> list:
        """
        Return only the data files the backend will actually use,
        based on the preset and extraction mode (wind vs temperature).
        """
        suffixes = ('.grb', '.grb2', '.grib', '.grib2', '.nc')

        if is_temp:
            # Temperature: look for T_2M (COSMO/ICON-EU) or T (ICON-Global) files
            if preset_name == 'COSMO REA6':
                return sorted(
                    p for p in data_dir.iterdir()
                    if p.name.startswith('T_2M') and p.suffix.lower() in suffixes
                )
            elif 'ICON' in preset_name:
                token = 'T_2M' if preset_name == 'ICON-DREAM-EU' else '_T_'
                return sorted(
                    p for p in data_dir.iterdir()
                    if token in p.name and p.suffix.lower() in suffixes
                )
            else:
                # ERA5: NetCDF files only for packaged builds.
                return sorted(
                    p for p in data_dir.iterdir()
                    if p.suffix.lower() == '.nc'
                )

        # Wind mode
        if preset_name == 'COSMO REA6':
            # COSMO: try UV/ subdirectory first, fall back to flat directory
            uv_dir = data_dir / 'UV'
            search_dir = uv_dir if uv_dir.is_dir() else data_dir
            return sorted(
                p for p in search_dir.iterdir()
                if p.name.startswith('U_') and p.suffix.lower() in suffixes
            )[:1]
        elif 'ICON' in preset_name:
            # ICON: U_10M, V_10M, PS files in flat directory
            return sorted(
                p for p in data_dir.iterdir()
                if 'U_10M' in p.name and p.suffix.lower() in suffixes
            )[:1]
        else:
            # ERA5: NetCDF files only for packaged builds.
            return sorted(
                p for p in data_dir.iterdir()
                if p.suffix.lower() == '.nc'
            )[:1]

    def _on_complete(self, result):
        n = result.get("n_timesteps", "?")
        is_temp = (self.ext_type_menu.get() == EXT_TYPE_TEMP)
        label = "Temperature field" if is_temp else "Wind fields"
        msg = f"{label} complete - {n} timesteps written"
        val_files = result.get("validation_files", [])
        if val_files:
            msg += f"\nValidation CSVs: {len(val_files)} file(s) in validation/ subfolder"
        show_success_in_log(self.log_textbox, msg)
        self.status_label.configure(text=f"Complete - {n} timesteps")
        set_widgets_state(self._interactive_widgets + [self.run_btn], "normal")

    def _on_error(self, exc):
        self.status_label.configure(text=f"Error: {exc}")
        set_widgets_state(self._interactive_widgets + [self.run_btn], "normal")
        show_error_dialog(self, "Meteo Generation Error", str(exc))

    def _on_cancel(self):
        self.status_label.configure(text="Cancelled")
        set_widgets_state(self._interactive_widgets + [self.run_btn], "normal")

    def _collect_settings(self) -> dict:
        return {
            "module": MODULE_ID,
            "version": SETTINGS_VERSION,
            "values": {
                "ext_type": self.ext_type_menu.get(),
                "preset": self.preset_menu.get(),
                "data_dir": self.data_entry.get().strip(),
                "aux_file": self.aux_entry.get().strip(),
                "grid_file": self.grid_entry.get().strip(),
                "crs_type": self.crs_menu.get(),
                "utm_zone": self.utm_entry.get().strip(),
                "grid_dx": self.dx_entry.get().strip(),
                "grid_dy": self.dy_entry.get().strip(),
                "reference_time": self.ref_time_entry.get().strip(),
                "output_name": self.name_entry.get().strip(),
                "output_dir": self.outdir_entry.get().strip(),
                "write_wnd": bool(self.wnd_var.get()),
                "validation_point_count": self._get_validation_count(),
                "validation_points": [
                    {
                        "name": row["name_entry"].get().strip(),
                        "lat": row["lat_entry"].get().strip(),
                        "lon": row["lon_entry"].get().strip(),
                    }
                    for row in self._validation_rows[:self._get_validation_count()]
                ],
            },
        }

    def _apply_settings(self, settings: dict) -> None:
        if settings.get("module") != MODULE_ID:
            raise ValueError(
                f"Settings file is for '{settings.get('module')}', "
                f"not '{MODULE_ID}'"
            )
        values = settings.get("values", {})

        def _set_entry(entry, value):
            entry.delete(0, "end")
            entry.insert(0, str(value))

        ext_type = values.get("ext_type", EXT_TYPE_WIND)
        if ext_type in EXTRACTION_TYPES:
            self.ext_type_menu.set(ext_type)
            self._on_ext_type_change(ext_type)

        preset = values.get("preset", WIND_PRESET_NAMES[0])
        if preset == "ERA5":
            preset = "ERA5 (.nc)"
        if preset in WIND_PRESET_NAMES:
            self.preset_menu.set(preset)
            self._on_preset_change(preset)

        _set_entry(self.data_entry, values.get("data_dir", ""))
        if get_wind_preset(preset).needs_aux_file and not self._auto_aux_path:
            _set_entry(self.aux_entry, values.get("aux_file", ""))
        else:
            self.aux_entry.delete(0, "end")
        _set_entry(self.grid_entry, values.get("grid_file", ""))

        crs_type = values.get("crs_type", "Cartesian")
        if crs_type in ("Cartesian", "Spherical"):
            self.crs_menu.set(crs_type)

        _set_entry(self.utm_entry, values.get("utm_zone", "32"))
        _set_entry(self.dx_entry, values.get("grid_dx", "6000"))
        _set_entry(self.dy_entry, values.get("grid_dy", "6000"))
        _set_entry(self.ref_time_entry, values.get("reference_time", ""))
        _set_entry(self.name_entry, values.get("output_name", ""))
        _set_entry(self.outdir_entry, values.get("output_dir", ""))
        self.wnd_var.set(bool(values.get("write_wnd", True)))

        for row in self._validation_rows:
            row["name_entry"].delete(0, "end")
            row["lat_entry"].delete(0, "end")
            row["lon_entry"].delete(0, "end")
        try:
            validation_count = max(
                0,
                min(MAX_VALIDATION_POINTS, int(values.get("validation_point_count", 0))),
            )
        except ValueError:
            validation_count = 0
        self.val_count_menu.set(str(validation_count))
        self._on_validation_count_change(str(validation_count))

        for idx, point in enumerate(values.get("validation_points", [])[:MAX_VALIDATION_POINTS]):
            row = self._validation_rows[idx]
            _set_entry(row["name_entry"], point.get("name", ""))
            _set_entry(row["lat_entry"], point.get("lat", ""))
            _set_entry(row["lon_entry"], point.get("lon", ""))

    def _on_save_settings(self):
        path = save_settings_dialog(self, default_name=f"{MODULE_ID}_settings.json")
        if not path:
            return
        try:
            save_json_settings(path, self._collect_settings())
            self.status_label.configure(text=f"Settings saved -> {Path(path).name}")
        except Exception as exc:
            show_error_dialog(self, "Save Failed", f"Could not save settings:\n{exc}")

    def _on_load_settings(self):
        path = load_settings_dialog(self)
        if not path:
            return
        try:
            data = load_json_settings(path)
            self._apply_settings(data)
            self.status_label.configure(text=f"Settings loaded <- {Path(path).name}")
        except Exception as exc:
            show_error_dialog(self, "Load Failed", f"Could not load settings:\n{exc}")

    def _on_reset(self):
        if self._task is not None and self._task.is_running:
            self._task.cancel()
            self.status_label.configure(text="Stopping - will reset when finished...")
            self.after(200, self._reset_when_idle)
            return
        self._reset_ui_to_defaults()

    def _reset_when_idle(self):
        if self._task is not None and self._task.is_running:
            self.after(200, self._reset_when_idle)
            return
        self._reset_ui_to_defaults()

    def _reset_ui_to_defaults(self):
        self.ext_type_menu.set(EXT_TYPE_WIND)
        self._on_ext_type_change(EXT_TYPE_WIND)
        self.preset_menu.set(WIND_PRESET_NAMES[0])
        self._on_preset_change(WIND_PRESET_NAMES[0])

        for entry in (
            self.data_entry, self.aux_entry, self.grid_entry,
            self.ref_time_entry, self.name_entry, self.outdir_entry,
        ):
            entry.delete(0, "end")

        self.crs_menu.set("Cartesian")

        for entry, default in (
            (self.utm_entry, "32"),
            (self.dx_entry, "6000"),
            (self.dy_entry, "6000"),
        ):
            entry.delete(0, "end")
            entry.insert(0, default)

        self.wnd_var.set(True)
        self.val_count_menu.set("0")
        self._on_validation_count_change("0")
        for row in self._validation_rows:
            row["name_entry"].delete(0, "end")
            row["lat_entry"].delete(0, "end")
            row["lon_entry"].delete(0, "end")

        clear_log(self.log_textbox)
        if self.progress_bar:
            self.progress_bar.set(0)
        self.eta_label.configure(text="")
        self.status_label.configure(text="Ready")
        set_widgets_state(self._interactive_widgets + [self.run_btn], "normal")
        self._task = None


if __name__ == "__main__":
    apply_theme()
    root = ctk.CTk()
    root.withdraw()
    win = WindGeneratorWindow(root)
    win.protocol("WM_DELETE_WINDOW", root.destroy)
    root.mainloop()
