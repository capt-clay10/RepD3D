# -*- coding: utf-8 -*-
"""
gui_utils.py — Thread-safe utilities for RepD3D v2 GUI.

This module solves the central GUI challenge: long-running backend
operations (BC extraction, wind processing) must not freeze the UI,
but CustomTkinter widgets can ONLY be modified from the main thread.

Architecture:
    1. BackgroundTask runs a backend function in a daemon thread
    2. The backend communicates via Python's logging module (thread-safe)
    3. A QueueHandler routes log records into a queue.Queue
    4. The main thread polls the queue via root.after(100, poll)
       and appends messages to the CTkTextbox
    5. Progress updates use a separate queue (also polled via .after)
    6. On completion/error, a callback fires ON THE MAIN THREAD

v1 had two problems:
    - ConsoleRedirector called widget.insert() from worker threads (unsafe)
    - stdout was globally redirected, breaking concurrent modules
    v2 fixes both by using logging → queue → after() polling.

Usage:
    from gui_utils import BackgroundTask, create_sub_window_layout, append_log

    # In a sub-module GUI:
    task = BackgroundTask(
        target=backend_bc_extractor.extract,
        kwargs={'model_def_path': '...', ...},
        log_textbox=self.log_textbox,
        progress_bar=self.progress_bar,
        eta_label=self.eta_label,
        on_complete=self._on_done,
        on_error=self._on_error,
        root=self,
    )
    task.start()
"""

import threading
import queue
import logging
import time
import json
from pathlib import Path
from typing import Callable, Optional, Any, Dict
from datetime import datetime

import customtkinter as ctk


# ─────────────────────────────────────────────────────────────────────────────
# CANCELLATION
# ─────────────────────────────────────────────────────────────────────────────

class TaskCancelledError(Exception):
    """Raised when a BackgroundTask is cancelled cooperatively via cancel_event."""
    pass


# ─────────────────────────────────────────────────────────────────────────────
# LOGGING INFRASTRUCTURE
# ─────────────────────────────────────────────────────────────────────────────

class QueueHandler(logging.Handler):
    """
    Logging handler that puts formatted log records into a queue.

    This is the bridge between the backend thread (which calls logger.info())
    and the GUI thread (which reads from the queue and updates the textbox).
    Thread-safe because queue.Queue is thread-safe.
    """

    def __init__(self, log_queue: queue.Queue):
        super().__init__()
        self.log_queue = log_queue

    def emit(self, record: logging.LogRecord):
        try:
            msg = self.format(record)
            self.log_queue.put(msg)
        except Exception:
            self.handleError(record)


def create_task_logger(name: str, log_queue: queue.Queue) -> logging.Logger:
    """
    Create a logger that routes to a queue (for GUI consumption).

    Args:
        name:      Logger name (should be unique per task).
        log_queue: Queue to receive formatted log strings.

    Returns:
        Configured logging.Logger instance.
    """
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    # Remove existing handlers to avoid duplicate messages
    logger.handlers.clear()

    handler = QueueHandler(log_queue)
    handler.setFormatter(logging.Formatter(
        "%(asctime)s  %(message)s",
        datefmt="%H:%M:%S",
    ))
    logger.addHandler(handler)

    # Don't propagate to root logger (avoids console spam)
    logger.propagate = False

    return logger


# ─────────────────────────────────────────────────────────────────────────────
# ETA FORMATTING
# ─────────────────────────────────────────────────────────────────────────────

def format_duration(seconds: float) -> str:
    """Format seconds as hours, minutes, and seconds."""
    if seconds is None or seconds < 0:
        return "0s"
    total_seconds = int(round(seconds))
    hours, remainder = divmod(total_seconds, 3600)
    minutes, secs = divmod(remainder, 60)
    if hours:
        return f"{hours}h {minutes:02d}m {secs:02d}s"
    if minutes:
        return f"{minutes}m {secs:02d}s"
    return f"{secs}s"


def format_eta(seconds: float) -> str:
    """
    Format a remaining-time estimate as a compact human string.

    Examples:
        5      → "ETA: 5s"
        18     → "ETA: 18s"
        135    → "ETA: 2m 15s"
        3725   → "ETA: 1h 02m"
    """
    if seconds is None or seconds < 0:
        return "ETA: calculating..."
    secs = int(round(seconds))
    if secs < 1:
        return "ETA: <1s"
    return f"ETA: {format_duration(secs)}"


def _widget_exists(widget) -> bool:
    """Safely check whether a Tk/CTk widget still exists."""
    try:
        return widget is not None and bool(widget.winfo_exists())
    except Exception:
        return False


# ─────────────────────────────────────────────────────────────────────────────
# JSON SETTINGS HELPERS
# ─────────────────────────────────────────────────────────────────────────────

def save_json_settings(file_path: str, data: Dict[str, Any]) -> None:
    """
    Write a settings dict to a JSON file.

    Raises:
        OSError:       If the file cannot be written.
        TypeError:     If the dict contains non-serialisable values.
    """
    path = Path(file_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, default=str)


def load_json_settings(file_path: str) -> Dict[str, Any]:
    """
    Read a settings dict from a JSON file.

    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError:        If the file contains invalid JSON or is not a dict.
    """
    path = Path(file_path)
    if not path.exists():
        raise FileNotFoundError(f"Settings file not found: {file_path}")
    with open(path, "r", encoding="utf-8") as f:
        try:
            data = json.load(f)
        except json.JSONDecodeError as exc:
            raise ValueError(f"Invalid JSON in settings file: {exc}") from exc
    if not isinstance(data, dict):
        raise ValueError("Settings file must contain a JSON object")
    return data


# ─────────────────────────────────────────────────────────────────────────────
# BACKGROUND TASK RUNNER
# ─────────────────────────────────────────────────────────────────────────────

class BackgroundTask:
    """
    Run a backend function in a background thread with thread-safe GUI updates.

    The backend function MUST accept these keyword arguments:
        logger: logging.Logger          — for log output
        progress_callback: Callable     — call with float 0.0–1.0

    Example:
        def my_backend(input_file, output_dir, logger, progress_callback):
            logger.info("Starting...")
            for i in range(100):
                # ... do work ...
                progress_callback(i / 100)
            logger.info("Done!")
            return {'output_file': 'result.bct'}

        task = BackgroundTask(
            target=my_backend,
            kwargs={'input_file': 'data.nc', 'output_dir': '/out'},
            log_textbox=self.log_textbox,
            progress_bar=self.progress_bar,
            on_complete=lambda result: print(f"Done: {result}"),
            on_error=lambda exc: print(f"Error: {exc}"),
            root=self,
        )
        task.start()
    """

    # Sentinel values for completion/error/cancellation signalling
    _COMPLETE = "__TASK_COMPLETE__"
    _ERROR = "__TASK_ERROR__"
    _CANCELLED = "__TASK_CANCELLED__"

    # Thresholds before a real ETA is shown (avoids noisy early estimates)
    _ETA_MIN_FRACTION = 0.02
    _ETA_MIN_ELAPSED_S = 1.0

    def __init__(
        self,
        target: Callable,
        kwargs: Dict[str, Any],
        log_textbox: ctk.CTkTextbox,
        progress_bar: Optional[ctk.CTkProgressBar] = None,
        eta_label: Optional[ctk.CTkLabel] = None,
        on_complete: Optional[Callable[[Any], None]] = None,
        on_error: Optional[Callable[[Exception], None]] = None,
        on_cancel: Optional[Callable[[], None]] = None,
        root: Optional[ctk.CTk] = None,
        poll_interval_ms: int = 100,
    ):
        """
        Args:
            target:           Backend function to run in background thread.
            kwargs:           Keyword arguments for target (excluding logger/progress_callback/cancel_event).
            log_textbox:      CTkTextbox widget for log output.
            progress_bar:     CTkProgressBar widget (optional).
            eta_label:        CTkLabel to display estimated remaining time (optional).
            on_complete:      Callback(result) called on main thread when done.
            on_error:         Callback(exception) called on main thread on error.
            on_cancel:        Callback() called on main thread when the task
                              finishes after cancel() was requested. If None,
                              on_error is NOT called for cancellation — the UI
                              is simply returned to idle cleanly.
            root:             CTk or CTkToplevel for .after() scheduling.
            poll_interval_ms: How often to poll queues (default 100ms).
        """
        self.target = target
        self.kwargs = kwargs
        self.log_textbox = log_textbox
        self.progress_bar = progress_bar
        self.eta_label = eta_label
        self.on_complete = on_complete
        self.on_error = on_error
        self.on_cancel = on_cancel
        self.root = root
        self.poll_interval_ms = poll_interval_ms

        self._log_queue: queue.Queue = queue.Queue()
        self._progress_queue: queue.Queue = queue.Queue()
        self._running = False
        self._thread: Optional[threading.Thread] = None
        self._start_time: Optional[float] = None
        self._last_fraction: float = 0.0

        # Cooperative cancellation — backend polls this event in long loops
        self.cancel_event = threading.Event()

        # Create a unique logger for this task
        task_id = f"repd3d.task.{id(self)}"
        self.logger = create_task_logger(task_id, self._log_queue)

    @property
    def is_running(self) -> bool:
        """Check if the task is currently running."""
        return self._running

    def cancel(self):
        """
        Request cooperative cancellation of the running task.

        Sets the cancel_event flag and immediately marks the task as not
        running so the GUI can reset without waiting for the backend thread
        to finish. The daemon thread will exit on its own once it observes
        the flag or completes the current blocking I/O call.
        """
        if not self._running:
            return
        if not self.cancel_event.is_set():
            self.cancel_event.set()
            self.logger.info("Cancellation requested...")
        self._running = False
        if self.on_cancel:
            self.on_cancel()

    def start(self):
        """Launch the backend function in a background thread."""
        if self._running:
            self.logger.warning("Task is already running")
            return

        self._running = True
        self._start_time = time.time()
        self._last_fraction = 0.0
        self.cancel_event.clear()

        # Reset progress bar
        if _widget_exists(self.progress_bar):
            self.progress_bar.set(0)

        # Seed ETA label
        if _widget_exists(self.eta_label):
            self.eta_label.configure(text="ETA: calculating...")

        # Clear log
        if _widget_exists(self.log_textbox):
            self.log_textbox.configure(state="normal")
            self.log_textbox.delete("1.0", "end")
            self.log_textbox.configure(state="disabled")

        self._thread = threading.Thread(
            target=self._run_target,
            daemon=True,
            name=f"RepD3D-{self.target.__name__}",
        )
        self.logger.info(f"Started {self.target.__name__} at {datetime.now():%H:%M:%S}")
        self._thread.start()
        self._poll()

    def _run_target(self):
        """Execute the backend function (runs in background thread)."""
        start_time = time.time()
        try:
            result = self.target(
                logger=self.logger,
                progress_callback=self._report_progress,
                cancel_event=self.cancel_event,
                **self.kwargs,
            )
            elapsed = time.time() - start_time
            self.logger.info(f"Completed in {format_duration(elapsed)}")
            self._log_queue.put((self._COMPLETE, result))

        except Exception as exc:
            elapsed = time.time() - start_time
            if self.cancel_event.is_set() or isinstance(exc, TaskCancelledError):
                self.logger.info(f"Cancelled after {format_duration(elapsed)}")
                self._log_queue.put((self._CANCELLED, exc))
            else:
                self.logger.error(f"Failed after {format_duration(elapsed)}: {exc}")
                self._log_queue.put((self._ERROR, exc))

    def _report_progress(self, fraction: float):
        """
        Called from the backend thread to update progress.
        Thread-safe: uses a queue.

        Args:
            fraction: Progress value in [0.0, 1.0].
        """
        self._progress_queue.put(max(0.0, min(1.0, fraction)))

    def _poll(self):
        """
        Poll log and progress queues from the main thread.
        Scheduled via root.after() — guaranteed main-thread execution.
        """
        # ── Drain log messages ───────────────────────────────────────────
        messages_processed = 0
        while not self._log_queue.empty() and messages_processed < 50:
            msg = self._log_queue.get_nowait()

            # Check for completion/error/cancellation sentinel
            if isinstance(msg, tuple) and len(msg) == 2:
                sentinel, payload = msg

                # If cancel() already fired on_cancel and reset _running,
                # silently discard any late sentinel from the dying thread.
                if self.cancel_event.is_set():
                    return

                if sentinel == self._COMPLETE:
                    self._running = False
                    if _widget_exists(self.progress_bar):
                        self.progress_bar.set(1.0)
                    if _widget_exists(self.eta_label):
                        self.eta_label.configure(text="")
                    if self.on_complete:
                        self.on_complete(payload)
                    return
                elif sentinel == self._ERROR:
                    self._running = False
                    if _widget_exists(self.eta_label):
                        self.eta_label.configure(text="")
                    if self.on_error:
                        self.on_error(payload)
                    return
                elif sentinel == self._CANCELLED:
                    self._running = False
                    if _widget_exists(self.progress_bar):
                        self.progress_bar.set(0)
                    if _widget_exists(self.eta_label):
                        self.eta_label.configure(text="")
                    if self.on_cancel:
                        self.on_cancel()
                    elif self.on_error:
                        # Fall back to error handler so UI is re-enabled
                        self.on_error(payload)
                    return

            # Regular log message — append to textbox
            if _widget_exists(self.log_textbox):
                append_log(self.log_textbox, str(msg))
            messages_processed += 1

        # ── Drain progress updates ───────────────────────────────────────
        latest_progress = None
        while not self._progress_queue.empty():
            latest_progress = self._progress_queue.get_nowait()

        if latest_progress is not None:
            self._last_fraction = latest_progress
            if _widget_exists(self.progress_bar):
                self.progress_bar.set(latest_progress)

        # ── Update ETA label ─────────────────────────────────────────────
        self._update_eta_label()

        # ── Schedule next poll ───────────────────────────────────────────
        if self._running and _widget_exists(self.root):
            try:
                self.root.after(self.poll_interval_ms, self._poll)
            except Exception:
                self._running = False

    def _update_eta_label(self):
        """
        Recompute ETA from elapsed time and last known progress fraction.
        Runs every poll tick so the displayed ETA stays fresh even when
        the backend doesn't emit new progress for a while.
        """
        if not _widget_exists(self.eta_label) or self._start_time is None:
            return

        elapsed = time.time() - self._start_time
        fraction = self._last_fraction

        if fraction < self._ETA_MIN_FRACTION or elapsed < self._ETA_MIN_ELAPSED_S:
            self.eta_label.configure(text="ETA: calculating...")
            return

        remaining = elapsed * (1.0 - fraction) / fraction
        self.eta_label.configure(text=format_eta(remaining))


# ─────────────────────────────────────────────────────────────────────────────
# LOG TEXTBOX HELPERS
# ─────────────────────────────────────────────────────────────────────────────

def append_log(textbox: ctk.CTkTextbox, text: str):
    """
    Append text to a CTkTextbox (main thread only).

    This is the ONLY safe way to write to the log textbox.
    It temporarily enables the widget, inserts text, scrolls to end,
    and re-disables it.

    Args:
        textbox: The CTkTextbox widget (must be in 'disabled' state).
        text:    Text to append (newline added automatically).
    """
    textbox.configure(state="normal")
    textbox.insert("end", text + "\n")
    textbox.see("end")
    textbox.configure(state="disabled")


def clear_log(textbox: ctk.CTkTextbox):
    """Clear all text from a CTkTextbox."""
    textbox.configure(state="normal")
    textbox.delete("1.0", "end")
    textbox.configure(state="disabled")


# ─────────────────────────────────────────────────────────────────────────────
# STANDARD SUB-MODULE WINDOW LAYOUT
# ─────────────────────────────────────────────────────────────────────────────

def create_sub_window_layout(
    window: ctk.CTkToplevel,
    title: str = "Sub-Module",
    min_width: int = 1100,
    min_height: int = 700,
) -> tuple:
    """
    Create the standard two-panel layout used by ALL sub-module windows.

    Layout:
    ┌───────────────────────────────┬─────────────────────────┐
    │                               │  Log Output  (h3 label) │
    │  Left Panel                   │  ┌─────────────────────┐│
    │  (CTkScrollableFrame)         │  │ CTkTextbox           ││
    │                               │  │ (read-only log)      ││
    │  Contains step-by-step        │  │                      ││
    │  workflow widgets             │  └─────────────────────┘│
    │                               │                         │
    │                               │  [====== Progress =====]│
    │                               │  Status: Ready          │
    └───────────────────────────────┴─────────────────────────┘

    Args:
        window: The CTkToplevel window to populate.
        title:  Window title and left panel header.
        min_width:  Minimum window width.
        min_height: Minimum window height.

    Returns:
        (left_frame, log_textbox, progress_bar, status_label, eta_label)

        left_frame:   CTkScrollableFrame — add your workflow widgets here
        log_textbox:  CTkTextbox — pass to BackgroundTask
        progress_bar: CTkProgressBar — pass to BackgroundTask
        status_label: CTkLabel — update with status messages
        eta_label:    CTkLabel — pass to BackgroundTask for live ETA display
    """
    # Import here to avoid circular imports
    from gui_theme import (
        COLORS, FONTS,
        apply_window_icon,
        make_textbox, make_progressbar, make_label, make_scrollable_frame,
    )

    window.title(f"RepD3D — {title}")
    window.configure(fg_color=COLORS["bg"])
    window.minsize(min_width, min_height)
    apply_window_icon(window)

    # ── Main container using grid ────────────────────────────────────────
    container = ctk.CTkFrame(window, fg_color=COLORS["bg"], corner_radius=0)
    container.pack(fill="both", expand=True, padx=10, pady=10)
    container.grid_columnconfigure(0, weight=3)
    container.grid_columnconfigure(1, weight=0, minsize=440)
    container.grid_rowconfigure(0, weight=1)

    # ── Left panel: scrollable workflow ──────────────────────────────────
    left_frame = make_scrollable_frame(
        container,
        label_text=title,
    )
    left_frame.grid(row=0, column=0, sticky="nsew", padx=(0, 5))

    # ── Right panel: log + progress ──────────────────────────────────────
    right_frame = ctk.CTkFrame(
        container, fg_color=COLORS["bg_secondary"], corner_radius=12, width=440
    )
    right_frame.grid(row=0, column=1, sticky="ns", padx=(5, 0))
    right_frame.grid_columnconfigure(0, weight=1, minsize=410)
    right_frame.grid_rowconfigure(2, minsize=360)

    # Log header
    make_label(right_frame, "Log Output", font_key="h3").grid(
        row=0, column=0, sticky="ew", padx=14, pady=(14, 2),
    )
    make_label(
        right_frame,
        "Live progress, timing, and processing details",
        font_key="caption",
        color=COLORS["text_dim"],
    ).grid(row=1, column=0, sticky="ew", padx=14, pady=(0, 8))

    # Log textbox
    log_textbox = make_textbox(right_frame, height=360, width=410)
    log_textbox.grid(row=2, column=0, sticky="ew", padx=14)

    # Progress bar
    progress_bar = make_progressbar(right_frame, width=300)
    progress_bar.grid(row=3, column=0, sticky="ew", padx=14, pady=(12, 6))

    # Status label
    status_label = make_label(
        right_frame, "Ready", font_key="body_sm",
        color=COLORS["text_dim"],
    )
    status_label.grid(row=4, column=0, sticky="ew", padx=14, pady=(0, 2))

    # ETA label (updated live by BackgroundTask while a task is running)
    eta_label = make_label(
        right_frame, "", font_key="body_sm",
        color=COLORS["text_dim"],
    )
    eta_label.grid(row=5, column=0, sticky="ew", padx=14, pady=(0, 14))

    return left_frame, log_textbox, progress_bar, status_label, eta_label


# ─────────────────────────────────────────────────────────────────────────────
# FILE DIALOG HELPERS
# ─────────────────────────────────────────────────────────────────────────────

def browse_file(
    entry: ctk.CTkEntry,
    title: str = "Select file",
    filetypes: Optional[list] = None,
    status_label: Optional[ctk.CTkLabel] = None,
    parent=None,
) -> Optional[str]:
    """
    Open a file dialog and insert the selected path into an entry widget.

    Args:
        entry:        CTkEntry to populate with the selected path.
        title:        Dialog title.
        filetypes:    List of (label, pattern) tuples.
                      e.g., [("MDF files", "*.mdf"), ("All files", "*.*")]
        status_label: Optional label to update with filename.

    Returns:
        Selected file path, or None if cancelled.
    """
    from tkinter import filedialog

    if filetypes is None:
        filetypes = [("All files", "*.*")]

    dialog_parent = parent or entry.winfo_toplevel()
    path = filedialog.askopenfilename(
        parent=dialog_parent,
        title=title,
        filetypes=filetypes,
    )
    if path:
        entry.delete(0, "end")
        entry.insert(0, path)
        if status_label:
            from pathlib import Path
            status_label.configure(text=f"✓ {Path(path).name}")
    return path or None


def browse_files(
    entry: ctk.CTkEntry,
    title: str = "Select files",
    filetypes: Optional[list] = None,
    status_label: Optional[ctk.CTkLabel] = None,
    parent=None,
) -> Optional[list]:
    """
    Open a multi-file dialog and insert the count into an entry widget.

    Returns:
        List of selected file paths, or None if cancelled.
    """
    from tkinter import filedialog

    if filetypes is None:
        filetypes = [("All files", "*.*")]

    dialog_parent = parent or entry.winfo_toplevel()
    paths = filedialog.askopenfilenames(
        parent=dialog_parent,
        title=title,
        filetypes=filetypes,
    )
    if paths:
        entry.delete(0, "end")
        entry.insert(0, f"{len(paths)} file(s) selected")
        if status_label:
            status_label.configure(text=f"✓ {len(paths)} files loaded")
        return list(paths)
    return None


def save_settings_dialog(
    parent,
    title: str = "Save Settings",
    default_name: str = "settings.json",
) -> Optional[str]:
    """Open a 'Save As' dialog for a JSON settings file. Returns path or None."""
    from tkinter import filedialog
    path = filedialog.asksaveasfilename(
        parent=parent,
        title=title,
        defaultextension=".json",
        filetypes=[("JSON files", "*.json"), ("All files", "*.*")],
        initialfile=default_name,
    )
    return path or None


def load_settings_dialog(
    parent,
    title: str = "Load Settings",
) -> Optional[str]:
    """Open an 'Open' dialog for a JSON settings file. Returns path or None."""
    from tkinter import filedialog
    path = filedialog.askopenfilename(
        parent=parent,
        title=title,
        filetypes=[("JSON files", "*.json"), ("All files", "*.*")],
    )
    return path or None


def browse_directory(
    entry: ctk.CTkEntry,
    title: str = "Select directory",
    status_label: Optional[ctk.CTkLabel] = None,
    parent=None,
) -> Optional[str]:
    """
    Open a directory dialog and insert the path into an entry widget.

    Returns:
        Selected directory path, or None if cancelled.
    """
    from tkinter import filedialog

    dialog_parent = parent or entry.winfo_toplevel()
    path = filedialog.askdirectory(parent=dialog_parent, title=title)
    if path:
        entry.delete(0, "end")
        entry.insert(0, path)
        if status_label:
            from pathlib import Path
            status_label.configure(text=f"✓ {Path(path).name}/")
    return path or None


# ─────────────────────────────────────────────────────────────────────────────
# UI STATE HELPERS
# ─────────────────────────────────────────────────────────────────────────────

def set_widgets_state(widgets: list, state: str = "normal"):
    """
    Enable or disable a list of widgets.

    Args:
        widgets: List of CTk widgets.
        state:   'normal' or 'disabled'.
    """
    for widget in widgets:
        try:
            widget.configure(state=state)
        except Exception:
            pass


def show_error_dialog(parent, title: str, message: str):
    """Show a simple error message dialog."""
    from tkinter import messagebox
    messagebox.showerror(title, message, parent=parent)


def show_info_dialog(parent, title: str, message: str):
    """Show a simple info message dialog."""
    from tkinter import messagebox
    messagebox.showinfo(title, message, parent=parent)


def show_warning_continue_dialog(parent, title: str, message: str) -> bool:
    """Show a warning with Continue / Cancel buttons. Returns True if user continues."""
    from tkinter import messagebox
    return messagebox.askokcancel(title, message, icon='warning', parent=parent)


def show_success_in_log(log_textbox: ctk.CTkTextbox, message: str):
    """Append a success message to the log."""
    append_log(log_textbox, f"{'─' * 40}")
    append_log(log_textbox, f"  ✓  {message}")
    append_log(log_textbox, f"{'─' * 40}")


# ─────────────────────────────────────────────────────────────────────────────
# SELF-TEST
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    import time as _time

    print("=== gui_utils.py self-test ===")
    print("This test creates a window and runs a mock background task.\n")

    # ── Non-visual unit tests first ──────────────────────────────────────
    all_pass = True

    def check(name, result, expected):
        global all_pass
        ok = result == expected
        status = "PASS" if ok else "FAIL"
        if not ok:
            all_pass = False
        print(f"  {name}: [{status}]")
        if not ok:
            print(f"    got: {result!r}, expected: {expected!r}")

    # QueueHandler test
    print("QueueHandler:")
    q = queue.Queue()
    logger = create_task_logger("test.unit", q)
    logger.info("Hello from test")
    check("message in queue", q.empty(), False)
    msg = q.get_nowait()
    check("message contains text", "Hello from test" in msg, True)
    check("message has timestamp", ":" in msg, True)

    # format_eta test
    print("\nformat_eta:")
    check("18s",      format_eta(18),   "ETA: 18s")
    check("2m 15s",   format_eta(135),  "ETA: 2m 15s")
    check("1h 02m",   format_eta(3725), "ETA: 1h 02m")
    check("negative", format_eta(-1),   "ETA: calculating...")

    # append_log (needs CTk context)
    print("\nVisual test:")

    from gui_theme import apply_theme, COLORS
    apply_theme()

    root = ctk.CTk()
    root.title("gui_utils.py — BackgroundTask Test")
    root.geometry("900x600")
    root.configure(fg_color=COLORS["bg"])

    left_frame, log_textbox, progress_bar, status_label, eta_label = \
        create_sub_window_layout(root, title="Test Module")

    # Add some dummy widgets to the left panel
    from gui_theme import make_step_header, make_button, make_label

    make_step_header(left_frame, 1, "Mock Workflow").pack(anchor="w", pady=10)
    make_label(left_frame, "This is a test of the BackgroundTask system.").pack(anchor="w", pady=5)

    # Mock backend function
    def mock_backend(duration=3, logger=None, progress_callback=None):
        logger.info("Starting mock extraction...")
        steps = 20
        for i in range(steps):
            _time.sleep(duration / steps)
            progress_callback((i + 1) / steps)
            logger.info(f"Processing step {i + 1}/{steps}")
        logger.info("Mock extraction complete!")
        return {"output_file": "test_output.bct", "n_boundaries": 8}

    def on_complete(result):
        show_success_in_log(log_textbox, f"Output: {result}")
        status_label.configure(text=f"✓ Complete — {result.get('n_boundaries', '?')} boundaries")
        run_btn.configure(state="normal")

    def on_error(exc):
        status_label.configure(text=f"✗ Error: {exc}")
        run_btn.configure(state="normal")

    def start_task():
        run_btn.configure(state="disabled")
        status_label.configure(text="Running...")
        task = BackgroundTask(
            target=mock_backend,
            kwargs={'duration': 3},
            log_textbox=log_textbox,
            progress_bar=progress_bar,
            eta_label=eta_label,
            on_complete=on_complete,
            on_error=on_error,
            root=root,
        )
        task.start()

    run_btn = make_button(left_frame, "Run Mock Task", start_task, "primary", width=200)
    run_btn.pack(anchor="w", pady=20)

    make_label(
        left_frame,
        "Click the button above to test threading.\n"
        "The progress bar and log should update smoothly\n"
        "without freezing the window.",
        font_key="body_sm",
    ).pack(anchor="w", pady=5)

    print("  Window opened. Click 'Run Mock Task' to test threading.")
    print("  Close the window when done.\n")

    root.mainloop()

    # ── Summary ──────────────────────────────────────────────────────────
    print(f"\n{'=' * 50}")
    if all_pass:
        print("ALL UNIT TESTS PASSED")
        print("Visual test requires manual verification.")
    else:
        print("SOME TESTS FAILED — check output above")
    print(f"{'=' * 50}")
