"""A small native Mac launcher for bes-search."""
import argparse
import csv
from datetime import datetime
import json
import os
from pathlib import Path
import queue
import shlex
import subprocess
import sys
import threading
import time
import traceback

from . import lab, lab_results, lab_storage, lab_workflow
from .lab_view import STATUS_LABELS, coverage_brief, coverage_facts, donor_scope, percent, reason_label, result_status

BG = "#ffffff"
CARD = "#ffffff"
HEADER_BG = "#f2f0ec"
TEXT = "#202c38"
MUTED = "#64717e"
ACCENT = "#2b4162"
ACCENT_DARK = "#172b43"
BORDER = "#e0e5e9"
SOFT = "#eef3f9"
WARNING = "#76591f"
MONO = "Menlo"
FONT = "Helvetica Neue"

ROW_CAP = 5000

STAGE_LABELS = [
    "Prepare inputs",
    "Check inputs and checksums",
    "Prepare reference chromosomes",
    "Search reference chromosomes (Cas-OFFinder)",
    "Reconstruct sample variant neighborhoods",
    "Compare loci and target sequences",
    "Write reports and check outputs",
]
VERIFY_STAGE = "Check expected example results"
VERIFY_MESSAGE = "Comparing scientific results with the accepted example"
MESSAGE_TO_STAGE = {
    "Checking input files and checksums": 1,
    "Preparing reference chromosomes": 2,
    "Searching reference chromosomes with Cas-OFFinder": 3,
    "Reconstructing sample variant neighborhoods": 4,
    "Comparing loci and target sequences": 5,
    "Writing reports and checking completed output files": 6,
}

TABLE_SOURCES = [
    ("Candidates", "candidates.tsv"),
    ("Sample differences", "locus-changes.tsv"),
    ("Excluded regions", "exclusions.tsv"),
]
STATUS_VALUES = ["retained", "gained", "lost", "unresolved", "unpaired_indel_path", "reference_only_run"]
COLUMN_SPECS = {
    "guide_id": ("Guide", 95), "location": ("Location · 0-based", 230),
    "strand": ("Strand", 80), "mismatches": ("Mismatches", 125),
    "target_sequence": ("Sequence · PAM", 300),
    "presence": ("Sequence in", 142),
    "reference_present": ("In reference", 104), "sample_copies": ("Donor copies", 110),
    "status": ("Locus status", 185), "sequence_changed": ("Sequence changed", 180),
    "reason": ("Reason", 260), "local_copy": ("Local copy", 100),
    "source_contig": ("Source chromosome", 150),
}
VIEW_COLUMNS = {
    "Candidates": ("guide_id", "location", "strand", "target_sequence", "mismatches"),
    "Sample differences": ("guide_id", "location", "strand", "status", "sequence_changed"),
    "Excluded regions": ("location", "reason"),
}
VIEW_NOTES = {
    "Candidates": "0-based coordinates, end excluded · target sequences include PAM",
    "Sample differences": "Changed sequence applies to retained loci. Unresolved loci are not losses.",
    "Excluded regions": "Excluded regions and phase notes. Records and reasons can overlap.",
}


def display_value(row, column):
    if column == "location":
        start, end = row.get("start0", row.get("start", "")), row.get("end0", row.get("end", ""))
        def coordinate(value):
            try:
                return f"{int(value):,}"
            except (ValueError, TypeError):
                return str(value) or "—"
        return f"{row.get('contig', '—')}:{coordinate(start)}–{coordinate(end)}"
    value = row.get(column, "")
    if column == "sequence_changed" and row.get("status") != "retained":
        return "—"
    if column in {"reference_present", "sequence_changed"}:
        return {"True": "Yes", "False": "No"}.get(str(value), value or "—")
    if column == "presence":
        try:
            copies = json.loads(row.get("sample_copies", "[]"))
        except (ValueError, TypeError):
            copies = []
        sources = (["Reference"] if str(row.get("reference_present")) == "True" else []) + (["sample"] if copies else [])
        return " + ".join(sources) or "—"
    if column == "target_sequence" and len(value) == 23:
        return value[:20] + " · " + value[20:]
    if column == "sample_copies":
        try:
            return ", ".join(str(item) for item in json.loads(value)) or "—"
        except (ValueError, TypeError):
            return value or "—"
    if column == "status":
        return STATUS_LABELS.get(value, value.replace("_", " "))
    if column == "reason":
        return reason_label(value)
    return value or "—"


def configure_tk():
    """Locate Tcl/Tk in the included relocatable Python runtime."""
    for name, directory in (("TCL_LIBRARY", "tcl9.0"), ("TK_LIBRARY", "tk9.0")):
        path = Path(sys.base_prefix) / "lib" / directory
        if not path.is_dir():
            raise ValueError("The bundled Tcl/Tk files are missing. Extract a fresh handoff and run Start.command.")
        os.environ[name] = str(path)


def apply_style(root):
    from tkinter import font, ttk
    style = ttk.Style(root)
    style.theme_use("clam")
    root.configure(background=BG)
    for name in ("TkDefaultFont", "TkTextFont", "TkHeadingFont", "TkMenuFont"):
        font.nametofont(name).configure(family=FONT, size=13)
    font.nametofont("TkFixedFont").configure(family=MONO, size=12)
    style.configure(".", background=BG, foreground=TEXT, font=(FONT, 13))
    for name, bg in (("TFrame", BG), ("Card.TFrame", CARD), ("Soft.TFrame", SOFT)):
        style.configure(name, background=bg, borderwidth=0, relief="flat")
    for name, bg, fg in (("TLabel", BG, TEXT), ("Card.TLabel", CARD, TEXT),
                         ("Muted.TLabel", BG, MUTED), ("CardMuted.TLabel", CARD, MUTED),
                         ("Scope.TLabel", SOFT, TEXT)):
        style.configure(name, background=bg, foreground=fg)
    style.configure("Header.TFrame", background=HEADER_BG)
    style.configure("Header.TLabel", background=HEADER_BG, foreground=MUTED)
    style.configure("Title.TLabel", background=HEADER_BG, foreground=ACCENT_DARK, font=(FONT, 16, "bold"))
    style.configure("Page.TLabel", foreground=ACCENT_DARK, font=(FONT, 23, "bold"))
    style.configure("H2.TLabel", font=(FONT, 16, "bold"))
    style.configure("SetupTitle.TLabel", foreground=ACCENT_DARK, font=(FONT, 19, "bold"))
    style.configure("SetupField.TLabel", foreground=ACCENT_DARK, font=(FONT, 14, "bold"))
    style.configure("SetupNote.TLabel", foreground=MUTED, font=(FONT, 12))
    style.configure("Setup.Quiet.TButton", padding=(0, 7), anchor="w")
    style.configure("Setup.Header.TButton", padding=(10, 6), font=(FONT, 12))
    style.configure("Info.Quiet.TButton", padding=(6, 4), font=(FONT, 12))
    style.configure("Metric.TLabel", background=CARD, foreground=ACCENT_DARK, font=(FONT, 20, "bold"))
    style.configure("Caption.TLabel", background=CARD, foreground=MUTED, font=(FONT, 11))
    style.configure("Status.TLabel", background=SOFT, foreground=ACCENT_DARK, padding=(10, 5))
    style.configure("Warning.TLabel", background="#f6f0e1", foreground=WARNING, padding=(10, 5))
    style.configure("TLabelframe", background=BG, borderwidth=0, relief="flat")
    style.configure("TLabelframe.Label", background=BG, foreground=TEXT, font=(FONT, 16, "bold"))
    style.configure("TButton", padding=(13, 8), background=CARD, foreground=TEXT,
                    bordercolor=BORDER, lightcolor=CARD, darkcolor=CARD, borderwidth=1, relief="flat")
    style.map("TButton", background=[("active", SOFT)], bordercolor=[("focus", ACCENT)],
              foreground=[("disabled", MUTED)])
    style.configure("File.TButton", background="#f3f6fa", foreground=ACCENT_DARK,
                    bordercolor="#b7c4d3", lightcolor="#f3f6fa", darkcolor="#f3f6fa", relief="raised")
    style.map("File.TButton", background=[("pressed", "#dae4ef"), ("active", "#e4ebf4")],
              bordercolor=[("focus", ACCENT), ("active", "#8598ac")])
    style.configure("Primary.TButton", background=ACCENT, foreground=CARD, bordercolor=ACCENT)
    style.map("Primary.TButton", background=[("disabled", "#dce2e9"), ("pressed", ACCENT_DARK), ("active", ACCENT_DARK)],
              foreground=[("disabled", "#64717e"), ("!disabled", CARD)],
              bordercolor=[("disabled", "#dce2e9"), ("!disabled", ACCENT)])
    style.configure("Quiet.TButton", background=BG, borderwidth=0, foreground=ACCENT)
    style.configure("Header.TButton", background=HEADER_BG, borderwidth=0, foreground=ACCENT)
    style.configure("TCheckbutton", background=BG, padding=(0, 4))
    style.configure("TRadiobutton", background=BG, padding=(0, 4))
    style.layout("Tab.TRadiobutton", [("Radiobutton.padding", {"sticky": "nswe", "children":
                  [("Radiobutton.focus", {"sticky": "nswe", "children": [("Radiobutton.label", {"sticky": "nswe"})]})]})])
    style.configure("Tab.TRadiobutton", padding=(16, 10), background=BG, foreground=MUTED)
    style.map("Tab.TRadiobutton", background=[("selected", SOFT), ("active", "#f5f6f7")],
              foreground=[("selected", ACCENT_DARK)], font=[("selected", (FONT, 13, "bold"))])
    style.layout("Mode.TRadiobutton", style.layout("Tab.TRadiobutton"))
    style.configure("Mode.TRadiobutton", padding=(12, 8), background="#f6f7f8",
                    foreground=MUTED, anchor="center", focuscolor=ACCENT, focusthickness=1)
    style.map("Mode.TRadiobutton", background=[("selected", "#e4ebf4"), ("active", "#eaf0f6")],
              foreground=[("selected", ACCENT_DARK), ("active", ACCENT)],
              font=[("selected", (FONT, 13, "bold"))])
    style.configure("TEntry", padding=(8, 7), fieldbackground=CARD, bordercolor=BORDER, lightcolor=CARD, darkcolor=CARD)
    style.map("TEntry", bordercolor=[("focus", ACCENT)])
    style.configure("TCombobox", padding=(8, 7), fieldbackground=CARD, background=CARD,
                    bordercolor=BORDER, lightcolor=CARD, darkcolor=CARD, arrowsize=12)
    style.map("TCombobox", fieldbackground=[("readonly", CARD)], foreground=[("readonly", TEXT)])
    style.configure("TMenubutton", padding=(10, 7), background=CARD, bordercolor=BORDER,
                    lightcolor=CARD, darkcolor=CARD, relief="flat")
    style.configure("Treeview", background=CARD, fieldbackground=CARD, foreground=TEXT,
                    rowheight=34, borderwidth=0, font=(FONT, 12))
    style.layout("Treeview", [("Treeview.treearea", {"sticky": "nswe"})])
    style.configure("Treeview.Heading", background="#f5f6f7", foreground=MUTED,
                    font=(FONT, 11, "bold"), padding=(10, 9), relief="flat", borderwidth=0)
    style.map("Treeview.Heading", background=[("active", "#e2eaf2")])
    style.map("Treeview", background=[("selected", SOFT)], foreground=[("selected", ACCENT_DARK)])
    style.configure("TScrollbar", background="#c5cfd9", troughcolor=BG, borderwidth=0,
                    arrowsize=12, lightcolor=BG, darkcolor=BG)
    for direction, stick in (("Vertical", "ns"), ("Horizontal", "ew")):
        style.layout(direction + ".TScrollbar", [(direction + ".Scrollbar.trough", {
            "sticky": stick, "children": [(direction + ".Scrollbar.thumb", {"sticky": "nswe"})]})])
    style.configure("TNotebook", background=BG, borderwidth=0, bordercolor=BG,
                    lightcolor=BG, darkcolor=BG)
    style.configure("TNotebook.Tab", padding=(16, 9), background=BG, borderwidth=0,
                    bordercolor=BG, lightcolor=BG, darkcolor=BG)
    style.map("TNotebook.Tab", background=[("selected", SOFT)])


def format_bytes(count):
    return f"{count / (1024 ** 3):.2f} GiB"


def format_duration(seconds):
    seconds = int(seconds)
    if seconds < 60:
        return f"{seconds} s"
    return f"{seconds // 60} min {seconds % 60} s"


def sort_key(value):
    try:
        return (0, float(value), "")
    except (TypeError, ValueError):
        return (1, 0.0, str(value).lower())


def row_count_text(shown, matched, total):
    if shown < matched:
        if matched == total:
            return f"{shown:,} of {total:,} records shown · sort shown rows only"
        return f"{shown:,} of {matched:,} matches shown · {total:,} total · sort shown rows only"
    if matched != total:
        return f"{matched:,} of {total:,} records match · all matches shown"
    return f"{total:,} {'record' if total == 1 else 'records'} · all shown"


def load_rows(path, query="", filters=None, cap=ROW_CAP, choices=None, guide_counts=None):
    """Stream a result TSV, returning (shown rows, matched count, total count)."""
    shown, matched, total = [], 0, 0
    with Path(path).open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            total += 1
            for column, values in (choices or {}).items():
                if row.get(column):
                    values.add(row[column])
            if guide_counts is not None and row.get("guide_id"):
                guide_counts[row["guide_id"]] = guide_counts.get(row["guide_id"], 0) + 1
            if not lab_results.row_matches(row, query, filters):
                continue
            matched += 1
            if len(shown) < cap:
                shown.append(row)
    return shown, matched, total


def distinct_values(path, column):
    values = set()
    with Path(path).open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            values.add(row.get(column, ""))
    values.discard("")
    return sorted(values)


class Disclosure:
    """A collapsible section with a keyboard-focusable toggle button."""

    def __init__(self, parent, title, ttk, open=False, managed=True):
        self.ttk = ttk
        self.title = title
        self.open = open
        self.section = ttk.Frame(parent)
        if managed:
            self.section.pack(fill="x", pady=(6, 0))
        self.button = ttk.Button(self.section, text=self._label(), command=self.toggle, takefocus=True, style="Quiet.TButton")
        self.button.pack(anchor="w", pady=(8, 2))
        self.body = ttk.Frame(self.section)
        if open:
            self.body.pack(fill="x", expand=True)

    def _label(self):
        return ("▾ " if self.open else "▸ ") + self.title

    def toggle(self):
        self.open = not self.open
        self.button.configure(text=self._label())
        if self.open:
            self.body.pack(fill="x", expand=True)
        else:
            self.body.pack_forget()


class Launcher:
    def __init__(self, root, bundle, *, preferences=None):
        import tkinter as tk
        from tkinter import ttk, filedialog, messagebox
        self.tk, self.ttk = tk, ttk
        self.dialog, self.messages = filedialog, messagebox
        self.root, self.bundle = root, Path(bundle)
        self.backend = self.bundle / "backend/cas-offinder"
        self.events, self.busy, self.result = queue.Queue(), False, None
        self.last_run = None
        self.pending_root = None
        self.screens = {}
        self.screen = None
        root.title("Base Edit Safety")
        root.geometry("1120x780")
        root.minsize(920, 640)

        self.flow = tk.StringVar(value="own")
        self.example_kind = tk.StringVar(value="small-reference")
        self.reference = tk.StringVar()
        self.vcf = tk.StringVar()
        self.sample = tk.StringVar()
        self.preferences = Path(preferences) if preferences is not None else lab_storage.preferences_path()
        self.preference_error = ""
        try:
            location = lab_storage.read_location(self.preferences)
        except ValueError as exc:
            location, self.preference_error = "", str(exc)
        self.remembered_output = location
        self.output = tk.StringVar(value=location)
        self.cache = tk.StringVar(value=str(Path.home() / "Documents/BES-public-cache"))
        self.diploid = tk.BooleanVar(value=False)
        self.table_source = tk.StringVar(value=TABLE_SOURCES[0][0])
        self.table_query = tk.StringVar()
        self.table_filter = tk.StringVar()
        self.sort_column, self.sort_reverse = None, False
        self.table_generation = 0
        self.shown_rows, self.row_records = [], {}
        self.table_cache = {}
        self.filter_display = tk.StringVar()
        self.guide_filter = tk.StringVar(value="All guides")
        self.guide_values, self.guide_counts = {"All guides": ""}, {}
        self.editor_label = tk.StringVar()
        self.cell_label = tk.StringVar()
        self.variant_source = tk.StringVar()
        self.export_busy = self.preparing = False
        self.applied_table = None
        # A monochrome open-folder icon keeps file actions distinct from text links.
        self.folder_icon = tk.BitmapImage(master=root, foreground=ACCENT, data="""
            #define folder_width 20
            #define folder_height 18
            static unsigned char folder_bits[] = {
                0x00,0x00,0x00, 0x00,0x00,0x00, 0x00,0x00,0x00, 0x00,0x00,0x00,
                0xfc,0x00,0x00, 0x04,0x01,0x00, 0x04,0xfe,0x01, 0x04,0x00,0x01,
                0xf4,0xff,0x07, 0x14,0x00,0x04, 0x0c,0x00,0x02, 0x0c,0x00,0x02,
                0x0c,0x00,0x01, 0x04,0x00,0x01, 0xfc,0xff,0x00, 0x00,0x00,0x00,
                0x00,0x00,0x00, 0x00,0x00,0x00};
        """)

        header = self.header = ttk.Frame(root, style="Header.TFrame", padding=(24, 5))
        header.pack(fill="x")
        self.brand = ttk.Label(header, text="Base Edit Safety", style="Title.TLabel")
        self.brand.configure(font=(FONT, 14, "bold"))
        self.brand.pack(side="left")
        self.quickstart_button = ttk.Button(header, text="Quickstart", style="Setup.Header.TButton", command=lambda: self.open_path(self.bundle / "QUICKSTART.html"))
        self.quickstart_button.pack(side="right")
        self.open_results_button = ttk.Button(header, text="Open results…", style="Setup.Header.TButton", command=self.open_results)
        self.open_results_button.pack(side="right", padx=(0, 4))
        self.container = ttk.Frame(root)
        self.build_progress()
        self.container.pack(fill="both", expand=True, padx=28, pady=(8, 0))

        self.build_inputs()
        self.show("inputs")
        root.protocol("WM_DELETE_WINDOW", self.close)
        root.after(200, self.poll)

    # ------------------------------------------------------------------ screens

    def show(self, name):
        if name == "results" and self.last_run is None:
            name = "inputs"
        self.screen = name
        if name != "run":
            self.hide_run_details()
        if name == "run":
            self.build_run()
        elif name == "results":
            self.build_results()
        for key, frame in self.screens.items():
            if key != name:
                frame.pack_forget()
        self.screens[name].pack(fill="both", expand=True)
        self.update_progress()
        self.layout_footer()
        if name == "inputs":
            self.root.minsize(720, 540)
        else:
            self.root.minsize(920, 640)

    def build_progress(self):
        """Show progress without adding navigation or keyboard stops."""
        self.progress_row = self.ttk.Frame(self.root, style="Header.TFrame", padding=(28, 12, 28, 16))
        self.progress_row.pack(side="bottom", fill="x")
        self.progress_row.columnconfigure((0, 2), weight=1, uniform="footer")
        steps = self.progress_steps = self.ttk.Frame(self.progress_row, style="Header.TFrame")
        steps.grid(row=1, column=1)
        self.progress_row.bind("<Configure>", lambda _: self.layout_footer())
        self.progress_widgets = []
        for i in range(3):
            if i:
                line = self.tk.Frame(steps, background=BORDER, width=42, height=1, takefocus=False)
                line.pack(side="left", padx=12)
            circle = self.tk.Canvas(steps, width=24, height=24, background=HEADER_BG,
                                    highlightthickness=0, takefocus=False)
            circle.pack(side="left")
            label = self.ttk.Label(steps, style="Header.TLabel", takefocus=False)
            label.pack(side="left", padx=(6, 0))
            self.progress_widgets.append((circle, label))

    def layout_footer(self):
        """Keep progress centered and the setup action at the window's right edge."""
        if not hasattr(self, "continue_button"):
            return
        setup = self.screen == "inputs"
        needed = self.progress_steps.winfo_reqwidth() + 2 * self.input_footer.winfo_reqwidth() + 88
        stacked = setup and self.progress_row.winfo_width() < needed
        layout = (setup, stacked)
        if getattr(self, "footer_layout", None) == layout:
            return
        self.footer_layout = layout
        self.progress_row.rowconfigure(1, minsize=self.continue_button.winfo_reqheight())
        self.progress_steps.grid_configure(row=2 if stacked else 1, column=0 if stacked else 1,
                                           columnspan=3 if stacked else 1, pady=(10 if stacked else 0, 0))
        if setup:
            self.input_footer.grid(row=1, column=2, sticky="e")
            self.input_notice.grid(row=0, column=0, columnspan=3, sticky="ew")
        else:
            self.input_footer.grid_remove()
            self.input_notice.grid_remove()

    def update_progress(self):
        current = ("inputs", "run", "results").index(self.screen)
        self.progress_states = ["complete" if i < current else "current" if i == current else "upcoming"
                                for i in range(3)]
        # A previous successful run must never complete the current Run step.
        if self.screen == "results" and self.saved_text != "Saved":
            self.progress_states[1] = "attention"
        for i, ((circle, label), state) in enumerate(zip(self.progress_widgets, self.progress_states)):
            circle.delete("all")
            color = WARNING if state == "attention" else MUTED if state == "upcoming" else ACCENT
            fill = ACCENT if state == "current" else SOFT if state == "complete" else HEADER_BG
            circle.create_oval(1, 1, 23, 23, fill=fill, outline=color, width=1)
            symbol = "✓" if state == "complete" else "!" if state == "attention" else str(i + 1)
            circle.create_text(12, 12, text=symbol, fill=CARD if state == "current" else color,
                               font=(FONT, 11, "bold" if state != "upcoming" else "normal"))
            title = ("Inputs", "Run", "Results")[i]
            label.configure(text=title + (" · needs attention" if state == "attention" else ""),
                            foreground=ACCENT_DARK if state in {"current", "complete"} else color,
                            font=(FONT, 12, "bold" if state == "current" else "normal"))

    def replace_screen(self, name, frame):
        old = self.screens.get(name)
        if old is not None:
            old.destroy()
        self.screens[name] = frame

    def scrollable(self, parent):
        canvas = self.tk.Canvas(parent, highlightthickness=0, background=BG)
        scrollbar = self.ttk.Scrollbar(parent, orient="vertical", command=canvas.yview)
        parent.rowconfigure(0, weight=1)
        parent.columnconfigure(0, weight=1)
        scrollbar.grid(row=0, column=1, sticky="ns")
        canvas.grid(row=0, column=0, sticky="nsew")
        canvas.configure(yscrollcommand=lambda first, last: self.set_scrollbar(scrollbar, first, last))
        pane = self.ttk.Frame(canvas, padding=4)
        window = canvas.create_window((0, 0), window=pane, anchor="nw")
        pane.bind("<Configure>", lambda event: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.bind("<Configure>", lambda event: canvas.itemconfigure(window, width=event.width))
        def scroll(event):
            if canvas.winfo_exists() and not isinstance(event.widget, (self.tk.Text, self.ttk.Combobox, self.ttk.Treeview)):
                canvas.yview_scroll(-1 if event.delta > 0 else 1, "units")
        parent.winfo_toplevel().bind("<MouseWheel>", scroll)
        return pane

    def set_scrollbar(self, scrollbar, first, last):
        scrollbar.set(first, last)
        if float(first) <= 0 and float(last) >= 1:
            scrollbar.grid_remove()
        else:
            scrollbar.grid()

    def caption(self, parent, text, card=False):
        label = self.ttk.Label(parent, text=text, style="CardMuted.TLabel" if card else "Muted.TLabel",
                               wraplength=780, justify="left")
        label.pack(fill="x", pady=(0, 4))
        label.bind("<Configure>", lambda event, widget=label: widget.configure(wraplength=max(240, event.width)))
        return label

    # ------------------------------------------------------------------ inputs

    def build_inputs(self):
        ttk = self.ttk
        frame = ttk.Frame(self.container)
        column = self.input_column = ttk.Frame(frame)
        column.pack(fill="both", expand=True)
        content = ttk.Frame(column)
        content.pack(fill="both", expand=True)
        canvas = self.input_canvas = self.tk.Canvas(content, highlightthickness=0, background=BG, takefocus=False)
        scroll = ttk.Scrollbar(content, command=canvas.yview)
        content.columnconfigure(0, weight=1)
        content.rowconfigure(0, weight=1)
        canvas.grid(row=0, column=0, sticky="nsew")
        scroll.grid(row=0, column=1, sticky="ns")
        canvas.configure(yscrollcommand=lambda a, b: self.set_scrollbar(scroll, a, b))
        pane = self.input_pane = ttk.Frame(canvas, padding=(0, 0, 12, 8))
        window = canvas.create_window((0, 0), window=pane, anchor="nw")
        canvas.bind("<Configure>", lambda event: canvas.itemconfigure(window, width=event.width))
        pane.bind("<Configure>", lambda _: canvas.configure(scrollregion=canvas.bbox("all")))

        top = ttk.Frame(pane)
        top.pack(fill="x", pady=(2, 8))
        title = ttk.Frame(top)
        title.pack(side="left")
        ttk.Label(title, text="Set up analysis", style="SetupTitle.TLabel").pack(side="left")
        self.settings_button = ttk.Button(title, text="ⓘ", style="Info.Quiet.TButton",
                                         takefocus=True, command=self.show_search_settings)
        self.settings_button.pack(side="left", padx=(6, 0))
        # Tk 9.1 can expose a separate screen-reader name. Older Tk keeps a text label.
        try:
            root = self.root
            root.tk.call("tk", "accessible", "set_acc_name", self.settings_button, "Search settings")
        except self.tk.TclError:
            self.settings_button.configure(text="ⓘ Search settings")
        self.settings_button.bind("<Return>", lambda _: self.settings_button.invoke())
        modes = self.mode_selector = self.tk.Frame(pane, background="#c3ccd6", padx=1, pady=1)
        modes.pack(anchor="w", pady=(0, 16))
        modes.columnconfigure((0, 2), weight=1, uniform="mode")
        self.mode_buttons, self.mode_indicators = {}, {}
        for segment, (title, value) in enumerate((("My analysis", "own"), ("Try an example", "example"))):
            button = ttk.Radiobutton(modes, text=title, variable=self.flow, value=value,
                                     width=17, takefocus=True, style="Mode.TRadiobutton", command=self.update_flow)
            button.grid(row=0, column=segment * 2, sticky="nsew")
            indicator = self.tk.Frame(modes, height=3, background="#f6f7f8", takefocus=False)
            indicator.grid(row=1, column=segment * 2, sticky="ew")
            self.mode_buttons[value], self.mode_indicators[value] = button, indicator
            button.bind("<Return>", lambda _, widget=button: widget.invoke())
            for key, direction in (("Left", -1), ("Right", 1)):
                button.bind(f"<{key}>", lambda _, widget=button, step=direction:
                            widget.tk.call("ttk::button::RadioTraverse", widget, step))
        self.tk.Frame(modes, background="#c3ccd6", width=1, takefocus=False).grid(row=0, column=1, rowspan=2, sticky="ns")

        self.example_frame = ttk.Frame(pane)
        for title, value, note in (("Reference only", "small-reference", "Artificial DNA · a few seconds · inspect and export candidates"),
                                   ("Reference + sample", "small", "Artificial DNA · a few seconds · checks the local backend"),
                                   ("Public HG003", "public", "23 published guides · GRCh38 · about 10–30 minutes, plus preparation")):
            choice = ttk.Frame(self.example_frame, padding=(0, 8))
            choice.pack(fill="x", pady=(0, 6))
            ttk.Radiobutton(choice, text=title, variable=self.example_kind, value=value,
                            command=self.update_flow).pack(anchor="w")
            self.caption(choice, note)
        self.public_note = ttk.Frame(self.example_frame)
        self.caption(self.public_note, "Donor model: 22 autosomes with local exclusions; X/Y remain unmodeled. Downloads 1.13 GB once.")
        self.cache_row = ttk.Frame(self.public_note)
        self.path_control(self.cache_row, "Public input cache", self.cache, "cache")
        self.caption(self.cache_row, "Downloaded files are verified before reuse.")
        self.cache_row.pack(fill="x")

        self.own_frame = ttk.Frame(pane)
        guides = ttk.Frame(self.own_frame)
        guides.pack(fill="x", pady=(0, 16))
        ttk.Label(guides, text="Guides", style="SetupField.TLabel").pack(anchor="w")
        ttk.Label(guides, text="1–32 spacers · one per line · 20 bases, without PAM", style="SetupNote.TLabel").pack(anchor="w", pady=(2, 6))
        guide_row = ttk.Frame(guides)
        guide_row.pack(fill="x")
        guide_box = ttk.Frame(guide_row)
        guide_box.pack(side="left", fill="x", expand=True)
        guide_box.columnconfigure(0, weight=1)
        self.guides = self.tk.Text(guide_box, height=5, width=25, wrap="word", font=(MONO, 12),
                                  background=CARD, foreground=TEXT, borderwidth=0, padx=10, pady=8,
                                  highlightthickness=1, highlightcolor=ACCENT, highlightbackground=BORDER)
        self.guides.grid(row=0, column=0, sticky="ew")
        guide_scroll = ttk.Scrollbar(guide_box, command=self.guides.yview)
        guide_scroll.grid(row=0, column=1, sticky="ns")
        self.guides.configure(yscrollcommand=lambda a, b: self.set_scrollbar(guide_scroll, a, b))
        self.guide_file_button = ttk.Button(guide_row, text="Open guide file…", style="File.TButton",
                                            image=self.folder_icon, compound="left", command=self.choose_guides)
        self.guide_file_button.pack(side="left", anchor="n", padx=(10, 0))
        self.guides.bind("<Tab>", lambda _: self.focus_after_guides())
        self.guides.bind("<Shift-Tab>", lambda _: self.focus_after_guides(backward=True))
        self.guide_error = self.input_error_label(guides)

        reference = ttk.Frame(self.own_frame)
        reference.pack(fill="x", pady=(0, 10))
        ttk.Label(reference, text="Reference genome", style="SetupField.TLabel").pack(anchor="w")
        self.reference_entry = self.path_control(reference, "GRCh38 reference FASTA", self.reference, "reference")
        prepared = ttk.Frame(reference)
        prepared.pack(fill="x", pady=(0, 2))
        ttk.Button(prepared, text="Use prepared public reference", style="Setup.Quiet.TButton",
                   command=self.use_public_reference).pack(side="left")
        ttk.Button(prepared, text="Prepare reference…", style="Setup.Quiet.TButton",
                   command=self.prepare_reference).pack(side="left", padx=(20, 0))
        self.reference_error = self.input_error_label(reference)

        self.sample_disclosure = Disclosure(self.own_frame, "Add sample variants — optional", ttk)
        self.sample_disclosure.button.configure(style="Setup.Quiet.TButton")
        sample_body = self.sample_disclosure.body
        self.vcf_entry = self.path_control(sample_body, "GRCh38 variant file (VCF)", self.vcf, "vcf")
        self.vcf_entry.bind("<FocusOut>", lambda _: self.load_typed_vcf())
        self.vcf_entry.bind("<Return>", lambda _: self.load_typed_vcf())
        self.donor_controls = ttk.Frame(sample_body)
        self.sample_row = ttk.Frame(self.donor_controls)
        self.sample_row.pack(fill="x", pady=(6, 0))
        selection = ttk.Frame(self.sample_row)
        selection.pack(side="left", fill="x", expand=True)
        ttk.Label(selection, text="Sample").pack(anchor="w", pady=(0, 4))
        self.sample_box = ttk.Combobox(selection, state="readonly", textvariable=self.sample, width=24)
        self.sample_box.pack(anchor="w")
        ttk.Button(self.sample_row, text="Clear VCF", style="Quiet.TButton", command=self.clear_vcf).pack(side="left", anchor="s")
        self.sample_note = self.caption(self.donor_controls, "Select a sample explicitly.")
        self.diploid_check = ttk.Checkbutton(self.donor_controls, variable=self.diploid,
            text="I confirm two copies of each autosome (chr1–22).")
        self.diploid_check.pack(anchor="w", pady=(2, 0))
        self.caption(self.donor_controls, "Required for this workflow. X/Y remain unresolved.")
        advanced = Disclosure(self.donor_controls, "Advanced genotype policy", ttk)
        self.caption(advanced.body, "Explicit PS phase sets; PASS and '.' filters; no GQ/DP cutoff. Missing calls remain unresolved. Other ploidy requires the advanced workflow.")
        self.variant_error = self.input_error_label(sample_body)

        self.additional_disclosure = Disclosure(self.own_frame, "Additional details — optional", ttk)
        self.additional_disclosure.button.configure(style="Setup.Quiet.TButton")
        details = self.additional_disclosure.body
        self.experiment_disclosure = Disclosure(details, "Experiment details", ttk)
        labels = self.experiment_disclosure.body
        self.caption(labels, "Labels are saved with your export; they do not change the search.")
        for title, variable in (("Editor", self.editor_label), ("Cell source", self.cell_label), ("Variant source", self.variant_source)):
            ttk.Label(labels, text=title).pack(anchor="w", pady=(6, 4))
            ttk.Entry(labels, textvariable=variable, width=34).pack(anchor="w")
        self.requirements_disclosure = Disclosure(details, "File requirements", ttk)
        self.caption(self.requirements_disclosure.body, "Guides: plain text or TSV with guide_id and spacer_sequence columns. Reference: full uncompressed GRCh38, chr1–22/X/Y; an index is created if needed. VCF: ##reference=GRCh38 and matching chromosome names. Normalization runs locally and preserves ALT indexes and phase.")

        self.input_footer = ttk.Frame(self.progress_row, style="Header.TFrame")
        self.continue_button = ttk.Button(self.input_footer, text="Review & run →", style="Primary.TButton",
                                         padding=(13, 6), command=self.continue_to_run)
        self.continue_button.pack(side="right")
        self.input_notice = ttk.Frame(self.progress_row, style="Header.TFrame")
        self.input_error = self.input_error_label(self.input_notice)
        self.input_error.configure(style="Header.TLabel")
        self.loaded_vcf = None
        self.vcf_read_error = ""
        self.vcf.trace_add("write", self.vcf_changed)
        self.bind_input_navigation(pane)
        self.replace_screen("inputs", frame)
        self.update_flow()

    def input_error_label(self, parent):
        label = self.ttk.Label(parent, foreground="#9a4038", justify="left", wraplength=700)
        label.bind("<Configure>", lambda event: label.configure(wraplength=max(160, event.width)))
        return label

    def set_input_error(self, label, message):
        label.configure(text=message)
        if message:
            label.pack(fill="x", pady=(4, 0))
        else:
            label.pack_forget()

    def show_search_settings(self):
        popup = self.tk.Toplevel(self.root)
        popup.title("Search settings")
        popup.resizable(False, False)
        pane = self.ttk.Frame(popup, padding=20)
        pane.pack(fill="both", expand=True)
        self.ttk.Label(pane, text="Search settings", style="H2.TLabel").pack(anchor="w", pady=(0, 10))
        self.ttk.Label(pane, text=self.search_scope(), wraplength=410, justify="left").pack(anchor="w")
        def close():
            popup.destroy()
            self.settings_button.focus_set()
        button = self.ttk.Button(pane, text="Close", command=close)
        button.pack(anchor="e", pady=(16, 0))
        popup.protocol("WM_DELETE_WINDOW", close)
        popup.bind("<Escape>", lambda _: close())
        popup.bind("<Return>", lambda _: close())
        popup.transient(self.root)
        popup.grab_set()
        button.focus_set()

    def focus_after_guides(self, backward=False):
        target = self.guides.tk_focusPrev() if backward else self.guide_file_button
        target.focus_set()
        return "break"

    def bind_input_navigation(self, widget):
        """Scroll the form under the pointer and reveal fields reached by Tab."""
        widget.bind("<FocusIn>", self.reveal_input, add="+")
        if not isinstance(widget, (self.tk.Text, self.ttk.Combobox)):
            widget.bind("<MouseWheel>", self.scroll_inputs)
        for child in widget.winfo_children():
            self.bind_input_navigation(child)
        self.input_canvas.bind("<MouseWheel>", self.scroll_inputs)

    def scroll_inputs(self, event):
        first, last = self.input_canvas.yview()
        if last - first < 1:
            self.input_canvas.yview_scroll(-1 if event.delta > 0 else 1, "units")
        return "break"

    def reveal_input(self, event):
        if self.screen != "inputs":
            return
        widget = event.widget
        canvas = self.input_canvas
        canvas.update_idletasks()
        top = widget.winfo_rooty() - self.input_pane.winfo_rooty()
        bottom = top + widget.winfo_height()
        visible = canvas.canvasy(0)
        height = canvas.winfo_height()
        total = max(1, self.input_pane.winfo_height())
        if top < visible:
            canvas.yview_moveto(max(0, top - 8) / total)
        elif bottom > visible + height:
            canvas.yview_moveto((bottom + 8 - height) / total)

    def update_flow(self):
        own = self.flow.get() == "own"
        for value, indicator in self.mode_indicators.items():
            indicator.configure(background=ACCENT if self.flow.get() == value else "#f6f7f8")
        public = not own and self.example_kind.get() == "public"
        self.own_frame.pack_forget()
        self.example_frame.pack_forget()
        (self.own_frame if own else self.example_frame).pack(fill="x", pady=4)
        if public:
            self.public_note.pack(fill="x")
        else:
            self.public_note.pack_forget()
        self.set_input_error(self.input_error, "")
        self.input_canvas.yview_moveto(0)

    def search_scope(self):
        limit = 2 if self.flow.get() == "example" and self.example_kind.get().startswith("small") else 6
        return f"20-base guides · NGG · ≤{limit} mismatches · both strands · no gaps"

    def vcf_changed(self, *_):
        if self.vcf.get().strip() == self.loaded_vcf:
            return
        self.loaded_vcf = None
        self.vcf_read_error = ""
        self.set_input_error(self.variant_error, "")
        self.sample.set("")
        self.sample_box.configure(values=[])
        self.diploid.set(False)
        self.sample_note.configure(text="Select a sample explicitly.")
        if self.vcf.get().strip():
            if not self.sample_disclosure.open:
                self.sample_disclosure.toggle()
            self.donor_controls.pack(fill="x", pady=(4, 0))
        else:
            self.donor_controls.pack_forget()

    def load_typed_vcf(self):
        path = self.vcf.get().strip()
        if path and path != self.loaded_vcf:
            self.on_vcf_chosen(path)

    def path_control(self, parent, label, variable, kind):
        ttk = self.ttk
        ttk.Label(parent, text=label).pack(anchor="w", pady=(6, 4))
        row = ttk.Frame(parent)
        row.pack(fill="x", pady=(0, 4))
        entry = ttk.Entry(row, textvariable=variable, width=18)
        entry.pack(side="left", fill="x", expand=True)
        ttk.Button(row, text="Browse…", style="File.TButton", image=self.folder_icon, compound="left",
                   command=lambda: self.choose_path(kind, variable)).pack(side="left", padx=(8, 0))
        return entry

    def choose_guides(self):
        path = self.dialog.askopenfilename(title="Choose guide text or TSV")
        if path:
            try:
                text = Path(path).read_text()
                lab.read_guides(text)
                self.guides.delete("1.0", "end")
                self.guides.insert("1.0", text)
                self.set_input_error(self.guide_error, "")
            except Exception as exc:
                self.set_input_error(self.guide_error, f"Cannot use {Path(path).name}: {exc}")
                self.guides.focus_set()

    def choose_path(self, kind, variable):
        path = (self.dialog.askdirectory(title="Choose a local folder", initialdir=variable.get()) if kind in {"cache", "output"}
                else self.dialog.askopenfilename(title="Choose " + kind))
        if path:
            variable.set(path)
            if kind == "vcf":
                self.on_vcf_chosen(path)

    def on_vcf_chosen(self, path):
        try:
            names = lab.samples(path)
        except Exception as exc:
            self.vcf_read_error = "Cannot read this VCF: " + str(exc)
            self.set_input_error(self.variant_error, self.vcf_read_error)
            return
        if path != self.loaded_vcf:
            self.sample.set("")
            self.diploid.set(False)
        self.loaded_vcf = path
        self.vcf_read_error = ""
        self.set_input_error(self.variant_error, "")
        self.donor_controls.pack(fill="x", pady=(4, 0))
        self.sample_box.configure(values=names)
        if not names:
            self.vcf_read_error = "This VCF has no sample genotype columns. Choose a VCF with sample genotypes."
            self.set_input_error(self.variant_error, self.vcf_read_error)
        else:
            self.sample_note.configure(text=f"{len(names)} sample(s) · select one explicitly.")

    def clear_vcf(self):
        self.vcf.set("")
        self.sample.set("")
        self.sample_box.configure(values=[])
        self.diploid.set(False)
        self.sample_note.configure(text="")

    def use_public_reference(self):
        path = Path(self.cache.get()).expanduser() / "out/primary.fa"
        if not path.is_file():
            self.prepare_reference()
        else:
            self.reference.set(str(path))

    def prepare_reference(self):
        if self.busy or self.preparing:
            return
        ttk = self.ttk
        popup = self.tk.Toplevel(self.root)
        popup.title("Prepare a reusable reference")
        popup.geometry("720x420")
        pane = ttk.Frame(popup, padding=24)
        pane.pack(fill="both", expand=True)
        ttk.Label(pane, text="Prepare once, reuse in your analyses", style="H2.TLabel").pack(anchor="w")
        self.caption(pane, "The existing public preparation downloads GRCh38 and the HG003 example inputs. "
                     "Only the reference is added to your analysis. No candidate search runs here.")
        self.path_control(pane, "Public input cache", self.cache, "cache")
        self.caption(pane, "1.13 GB download · 25 GiB free space required · cached inputs are verified")
        status = self.caption(pane, "Ready to prepare.")
        logs = Disclosure(pane, "Preparation log", ttk)
        log = self.detail_text(logs.body, "", height=5)
        button = ttk.Button(pane, text="Prepare & use reference", style="Primary.TButton")
        button.pack(anchor="e", pady=12)
        self.preparation_ui = (popup, status, log, button)
        def close():
            if not self.preparing:
                popup.destroy()
        popup.protocol("WM_DELETE_WINDOW", close)
        popup.bind("<Escape>", lambda _: close())
        def start():
            if self.preparing:
                return
            self.preparing = True
            button.configure(state="disabled")
            status.configure(text="Preparing public inputs…")
            cache = Path(self.cache.get()).expanduser().resolve()
            def worker():
                try:
                    request = lab.prepare_public(cache, self.backend, lambda line: self.events.put(("prepare_log", line)))
                    reference = str((cache / request["reference"]["path"]).resolve())
                    self.events.put(("prepare_done", reference))
                except Exception as exc:
                    self.events.put(("prepare_error", str(exc)))
            threading.Thread(target=worker, daemon=False).start()
        button.configure(command=start)
        popup.transient(self.root)
        popup.grab_set()
        popup.focus_set()

    def handle_preparation(self, kind, value):
        popup, status, log, button = self.preparation_ui
        if kind == "prepare_log":
            log.configure(state="normal")
            log.insert("end", value + "\n")
            log.see("end")
            log.configure(state="disabled")
            return
        self.preparing = False
        button.configure(state="normal")
        if kind == "prepare_done":
            self.reference.set(value)
            status.configure(text="Reference ready. It is selected for this and future analyses.")
            button.configure(text="Done", command=popup.destroy)
        else:
            status.configure(text=value)

    def open_results(self, folder=None):
        if self.busy or self.preparing or self.export_busy:
            return
        if folder is None:
            folder = self.dialog.askdirectory(title="Open a completed BES results folder")
        if not folder:
            return
        try:
            path, state, _, _ = lab_results.read_analysis(folder)
        except ValueError as exc:
            self.messages.showerror("Could not open results", str(exc))
            return
        self.last_run = (path, state)
        self.table_source.set("Candidates")
        self.table_query.set("")
        self.guide_filter.set("All guides")
        self.show("results")

    def collect_inputs(self):
        """Return (lab_mode, own kwargs) or None after showing what is missing."""
        self.input_error.configure(foreground="#9a4038")
        self.set_input_error(self.input_error, "")
        if self.flow.get() == "example":
            return self.example_kind.get(), {}
        guide_problem, reference_problem, variant_problems = "", "", []
        try:
            lab.read_guides(self.guides.get("1.0", "end"))
        except ValueError as exc:
            guide_problem = str(exc)
        reference = self.reference.get().strip()
        if not reference:
            reference_problem = "Choose a GRCh38 reference FASTA."
        else:
            path = Path(reference).expanduser()
            if not path.is_file():
                reference_problem = f"Reference file not found: {path}"
            elif path.suffix == ".gz":
                reference_problem = "Choose an uncompressed GRCh38 FASTA (.fa or .fasta)."
        vcf = self.vcf.get().strip()
        if vcf:
            self.load_typed_vcf()
            if not Path(vcf).expanduser().is_file():
                variant_problems.append(f"VCF file not found: {vcf}")
            elif self.vcf_read_error:
                variant_problems.append(self.vcf_read_error)
            else:
                if not self.sample.get():
                    variant_problems.append("Select a sample from the VCF.")
                if not self.diploid.get():
                    variant_problems.append("Confirm two copies of each autosome for this sample, or remove the VCF.")
        self.set_input_error(self.guide_error, guide_problem)
        self.set_input_error(self.reference_error, reference_problem)
        self.set_input_error(self.variant_error, " ".join(variant_problems))
        if variant_problems and not self.sample_disclosure.open:
            self.sample_disclosure.toggle()
        if guide_problem or reference_problem or variant_problems:
            self.set_input_error(self.input_error, "Check the highlighted inputs.")
            if guide_problem:
                target = self.guides
            elif reference_problem:
                target = self.reference_entry
            else:
                target = self.vcf_entry
            target.focus_set()
            return None
        return "own", dict(guide_text=self.guides.get("1.0", "end"), reference=reference,
                           vcf=vcf or None, sample=self.sample.get(), diploid=self.diploid.get())

    def continue_to_run(self):
        if self.busy or self.preparing:
            return
        collected = self.collect_inputs()
        if collected is None:
            return
        if not self.ensure_save_location():
            return
        self.lab_mode, self.own_kwargs = collected
        self.run_labels = {"editor": self.editor_label.get(), "cell_source": self.cell_label.get(),
                           "variant_source": self.variant_source.get()} if self.lab_mode == "own" else {}
        self.show("run")

    def remember_save_location(self, folder):
        try:
            folder = lab_storage.check_location(folder)
        except ValueError as exc:
            self.messages.showerror("Save location", str(exc))
            return False
        try:
            lab_storage.remember_location(folder, self.preferences)
        except OSError as exc:
            self.messages.showerror("Save location", f"Cannot remember this location. Check write access to {self.preferences.parent}: {exc}")
            return False
        self.output.set(str(folder))
        self.remembered_output = str(folder)
        self.preference_error = ""
        return True

    def choose_save_location(self):
        if self.busy:
            return False
        current = self.output.get()
        options = {"title": "Choose where BES saves analyses", "mustexist": True}
        if current and Path(current).is_dir():
            options["initialdir"] = current
        chosen = self.dialog.askdirectory(**options)
        return bool(chosen) and self.remember_save_location(chosen)

    def ensure_save_location(self):
        current = self.output.get()
        problem = self.preference_error
        if current:
            try:
                lab_storage.check_location(current)
            except ValueError as exc:
                problem = str(exc)
            else:
                return True if current == self.remembered_output else self.remember_save_location(current)
        if problem:
            self.messages.showerror("Save location needs attention", problem)
        if self.choose_save_location():
            self.set_input_error(self.input_error, "")
            return True
        self.set_input_error(self.input_error, problem or "Choose a save location to continue.")
        return False

    def change_save_location(self):
        if self.choose_save_location():
            self.run_preflight()

    # ------------------------------------------------------------------ run

    def build_run(self):
        ttk = self.ttk
        previous = getattr(self, "run_details_popup", None)
        if previous is not None and previous.winfo_exists():
            previous.destroy()
        frame = ttk.Frame(self.container)
        frame.rowconfigure(0, weight=1)
        frame.columnconfigure(0, weight=1)
        body = ttk.Frame(frame)
        body.grid(row=0, column=0, sticky="nsew")
        pane = self.scrollable(body)
        title = ttk.Frame(pane)
        title.pack(fill="x", pady=(0, 8))
        self.run_title = ttk.Label(title, text="Review & run", style="Page.TLabel")
        self.run_title.pack(side="left")
        self.run_info_button = self.info_button(title, self.show_run_details, label="Run details",
                                                state="normal", hint_to_right=True)
        self.run_info_button.pack(side="left", padx=(8, 0))
        self.run_summary = ttk.Label(pane, text=self.describe_run(), justify="left")
        self.run_summary.pack(fill="x", pady=(4, 8))
        self.run_summary.bind("<Configure>", lambda event: self.run_summary.configure(wraplength=max(240, event.width)))
        self.check_summary = ttk.Label(pane, text="Checking this Mac…", style="Caption.TLabel")
        self.check_summary.pack(anchor="w", pady=(8, 0))
        self.run_issues = ttk.Frame(pane)
        self.run_issues.pack(fill="x", pady=(4, 0))
        self.build_input_review(pane)
        self.stage_frame = ttk.LabelFrame(pane, text="Progress", padding=(0, 10))
        self.stage_list = ttk.Frame(self.stage_frame)
        self.stage_list.pack(fill="x")
        self.stage_rows = []
        self.stage_index = -1
        self.elapsed = ttk.Label(self.stage_frame, style="Muted.TLabel")
        self.elapsed.pack(anchor="w", pady=(6, 0))
        self.build_run_details_dialog()

        self.failure = ttk.LabelFrame(pane, text="Run did not finish", padding=12)
        self.failure_text = ttk.Label(self.failure, wraplength=780, justify="left")
        self.failure_text.pack(anchor="w")
        self.failure_partial = ttk.Label(self.failure, style="Muted.TLabel", wraplength=780, justify="left")
        self.failure_partial.pack(anchor="w", pady=(4, 0))
        failure_actions = ttk.Frame(self.failure)
        failure_actions.pack(fill="x", pady=(8, 0))
        ttk.Button(failure_actions, text="Back to inputs", style="Primary.TButton",
                   command=lambda: self.show("inputs")).pack(side="left")
        if self.last_run:
            ttk.Button(failure_actions, text="Return to previous results", command=lambda: self.show("results")).pack(side="left", padx=8)
        actions = ttk.Frame(frame)
        actions.grid(row=1, column=0, sticky="ew", pady=(8, 12))
        self.start = ttk.Button(actions, text="Start analysis", style="Primary.TButton", command=self.start_run)
        self.start.pack(side="right")
        self.back = ttk.Button(actions, text="Back to inputs", command=lambda: self.show("inputs"))
        self.back.pack(side="right", padx=8)

        self.replace_screen("run", frame)
        self.run_preflight()

    def review_guides(self):
        """Read only the small guide input; preparation stays in the runner."""
        if self.lab_mode == "own":
            return lab.read_guides(self.own_kwargs["guide_text"])
        if self.lab_mode.startswith("small"):
            path = lab.files() / "examples/small/guides.tsv"
        else:
            cache = Path(self.cache.get()).expanduser()
            request = cache / "request.json"
            path = (cache / json.loads(request.read_text())["guides"]["path"] if request.is_file()
                    else lab.files() / "examples/research/guides.tsv")
        return lab.read_guides(path.read_text())

    def build_input_review(self, pane):
        ttk = self.ttk
        review = self.input_review = ttk.Frame(pane)
        review.pack(fill="x", pady=(16, 0))
        ttk.Label(review, text="Guides to search", style="SetupField.TLabel").pack(anchor="w", pady=(0, 8))
        try:
            guides = self.review_guides()
            note = "Format checked · unique 20-base spacers · GC is the share of G/C bases, not a score."
        except (OSError, ValueError, KeyError, TypeError) as exc:
            guides = []
            note = f"Guide preview unavailable: {exc}. Inputs will be checked when the run starts."
        table = ttk.Frame(review)
        table.pack(fill="x")
        table.columnconfigure(0, weight=1)
        tree = self.review_guide_tree = ttk.Treeview(table, columns=("guide", "spacer", "gc"),
                                                     show="headings", height=min(4, max(1, len(guides))))
        for key, title, width, stretch in (("guide", "Guide", 180, True),
                                           ("spacer", "Spacer · 5′ → 3′", 320, True),
                                           ("gc", "GC", 90, False)):
            tree.heading(key, text=title, anchor="w")
            tree.column(key, width=width, minwidth=70 if key == "gc" else 140, stretch=stretch, anchor="w")
        tree.grid(row=0, column=0, sticky="ew")
        scrollbar = ttk.Scrollbar(table, command=tree.yview)
        scrollbar.grid(row=0, column=1, sticky="ns")
        tree.configure(yscrollcommand=lambda first, last: self.set_scrollbar(scrollbar, first, last))
        tree.tag_configure("alternate", background="#f8fafd")
        for index, guide in enumerate(guides):
            spacer = guide["spacer_sequence"]
            gc = (spacer.count("G") + spacer.count("C")) / len(spacer) * 100
            tree.insert("", "end", values=(guide["guide_id"], spacer, f"{gc:.0f}%"),
                        tags=("alternate",) if index % 2 else ())
        self.review_note = self.caption(review, note)
        self.review_note.configure(style="SetupNote.TLabel", foreground=MUTED if guides else WARNING)

        ttk.Label(review, text="What this run will do", style="SetupField.TLabel").pack(anchor="w", pady=(16, 8))
        small = self.lab_mode.startswith("small")
        sample = self.lab_mode in {"small", "public"} or bool(self.own_kwargs.get("vcf"))
        scope = "the artificial reference" if small else "the GRCh38 reference (chr1–22, X and Y)"
        plan = [("Search reference", f"Find potential off-target sites for each guide in {scope}.")]
        if sample:
            comparison = "Compare gained/lost loci and retained loci with changed sequences. Unresolved loci and unpaired paths stay distinct."
            if not small:
                comparison += " Sample model: two copies per autosome; X/Y copy counts are undeclared."
            plan.append(("Compare sample", comparison))
        else:
            plan.append(("Reference only", "No sample comparison is requested. You can inspect and export reference candidates."))
        plan.append(("Inspect & export", "Review site locations, strand, sequence/PAM and mismatches. Candidate counts and coverage are calculated during the run."))
        self.review_plan = ttk.Frame(review)
        self.review_plan.pack(fill="x")
        for title, description in plan:
            row = ttk.Frame(self.review_plan)
            row.pack(fill="x", pady=(0, 8))
            ttk.Label(row, text=title, width=18, style="SetupNote.TLabel").pack(side="left", anchor="n")
            detail = ttk.Label(row, text=description, justify="left", wraplength=640)
            detail.pack(side="left", fill="x", expand=True)
            detail.bind("<Configure>", lambda event, label=detail: label.configure(wraplength=max(240, event.width)))

    def build_run_details_dialog(self):
        ttk = self.ttk
        popup = self.run_details_popup = self.tk.Toplevel(self.root)
        popup.withdraw()
        popup.title("Run details")
        popup.geometry("880x620")
        popup.minsize(760, 500)
        popup.configure(background=BG)
        popup.protocol("WM_DELETE_WINDOW", self.hide_run_details)
        popup.bind("<Escape>", lambda _: self.hide_run_details())
        outer = ttk.Frame(popup, padding=20)
        outer.pack(fill="both", expand=True)
        ttk.Button(outer, text="Close", command=self.hide_run_details).pack(side="bottom", anchor="e", pady=(10, 0))
        tabs = self.run_info_tabs = ttk.Notebook(outer)
        tabs.pack(fill="both", expand=True)
        location_pane, checks_pane, log_pane = (ttk.Frame(tabs, padding=12) for _ in range(3))
        for panel, label in ((location_pane, "Save location"), (checks_pane, "Preflight checks"), (log_pane, "Technical log")):
            tabs.add(panel, text=label)
        location = ttk.Label(location_pane, textvariable=self.output, justify="left")
        location.pack(fill="x", pady=(0, 8))
        location.bind("<Configure>", lambda event: location.configure(wraplength=max(240, event.width)))
        self.caption(location_pane, "Runs are saved automatically in new dated folders. A location change applies to future runs.")
        self.change_output_button = ttk.Button(location_pane, text="Change save location", command=self.change_save_location)
        self.change_output_button.pack(anchor="w", pady=(8, 0))
        self.checks_frame = self.scrollable(checks_pane)
        self.check_labels = []
        self.details_disclosure = Disclosure(log_pane, "Error details", ttk, managed=False)
        self.details = self.detail_text(self.details_disclosure.body, "", height=8)
        self.log = self.detail_text(log_pane, "", height=12)

    def show_run_details(self):
        self.run_details_popup.transient(self.root)
        self.run_details_popup.deiconify()
        self.run_details_popup.lift()
        self.run_details_popup.focus_set()

    def hide_run_details(self):
        popup = getattr(self, "run_details_popup", None)
        if popup is not None and popup.winfo_exists():
            popup.withdraw()
            if self.screen == "run" and self.run_info_button.winfo_exists():
                self.run_info_button.focus_set()

    def describe_run(self):
        if self.lab_mode.startswith("small"):
            return "1 guide · artificial reference · " + ("sample SMALL" if self.lab_mode == "small" else "reference only")
        if self.lab_mode == "public":
            return "Public HG003 example — 23 published guides, all 24 chromosomes"
        parts = []
        try:
            parts.append(f"{len(lab.read_guides(self.own_kwargs['guide_text']))} guide(s)")
        except ValueError:
            parts.append("guides as entered")
        parts.append(Path(self.own_kwargs["reference"]).name)
        if self.own_kwargs.get("vcf"):
            parts.append(f"sample {self.own_kwargs['sample']}")
        else:
            parts.append("reference only, no sample VCF")
        return " — ".join(parts)

    def check_row(self, ok, label, detail):
        line = self.ttk.Frame(self.checks_frame)
        line.pack(fill="x", pady=1)
        self.ttk.Label(line, text="✓" if ok else "✗", width=3).pack(side="left", anchor="n")
        self.ttk.Label(line, text=label, width=28).pack(side="left", anchor="n")
        description = self.ttk.Label(line, text=detail, style="Muted.TLabel", wraplength=520, justify="left")
        description.pack(side="left", fill="x", expand=True)
        description.bind("<Configure>", lambda event: description.configure(wraplength=max(160, event.width)))
        if not ok:
            issue = self.ttk.Label(self.run_issues, text=f"{label}: {detail}", foreground=WARNING, justify="left")
            issue.pack(fill="x", pady=(4, 0))
            issue.bind("<Configure>", lambda event: issue.configure(wraplength=max(240, event.width)))
        self.check_labels.append(line)
        return ok

    def run_preflight(self):
        for child in self.checks_frame.winfo_children():
            child.destroy()
        for child in self.run_issues.winfo_children():
            child.destroy()
        self.check_labels = []
        mode = self.lab_mode
        backend_ok = self.backend.is_file() and os.access(self.backend, os.X_OK)
        self.check_row(backend_ok, "Cas-OFFinder backend",
                       "Found and executable." if backend_ok else
                       f"Missing or not executable at {self.backend}. Extract a fresh handoff and run Start.command.")
        try:
            lab_storage.check_location(self.output.get())
            info = lab.preflight(Path(self.output.get()).expanduser(), large=not mode.startswith("small"))
            self.check_row(True, "System resources",
                           f"macOS {info['os']} ({info['architecture']}); {format_bytes(info['available_memory_bytes'])} available RAM; "
                           f"{format_bytes(info['free_disk_bytes'])} free disk.")
            resources_ok = True
        except Exception as exc:
            self.check_row(False, "System resources", str(exc))
            resources_ok = False
        if mode.startswith("small"):
            present = all((lab.files() / "examples/small" / name).is_file()
                          for name in ("reference.fa", "guides.tsv", "sample.vcf"))
            self.check_row(present, "Bundled example files",
                           "Present." if present else "Missing — reinstall the package or handoff.")
            inputs_ok = present
        elif mode == "own":
            self.check_row(True, "Your inputs", "Guides, reference and optional VCF were checked on the previous screen.")
            inputs_ok = True
        else:
            self.check_row(True, "Public inputs",
                           "Missing inputs are downloaded and verified automatically. "
                           "Verified cached files are reused; changed files are never accepted.")
            inputs_ok = True
        ready = backend_ok and resources_ok and inputs_ok
        self.check_summary.configure(text="✓ Ready to run" if ready else "! Action needed before starting",
                                     foreground=ACCENT if ready else WARNING)
        self.check_summary.pack(anchor="w", pady=(8, 0), before=self.run_issues)
        self.start.configure(state="normal" if ready else "disabled")
        if not resources_ok:
            self.ttk.Button(self.run_issues, text="Change save location", command=self.change_save_location).pack(anchor="w", pady=(8, 0))
        return ready

    def stage_names(self):
        names = list(STAGE_LABELS)
        if self.lab_mode in {"small", "public"}:
            names.append(VERIFY_STAGE)
        return names

    def build_stages(self):
        self.check_summary.pack_forget()
        self.input_review.pack_forget()
        self.stage_frame.pack(fill="x", pady=(12, 0))
        for row in self.stage_rows:
            row["line"].destroy()
        self.stage_rows = []
        self.stage_index = -1
        for name in self.stage_names():
            line = self.ttk.Frame(self.stage_list)
            line.pack(fill="x", pady=5)
            marker = self.ttk.Label(line, text="○", width=3)
            marker.pack(side="left")
            label = self.ttk.Label(line, text=name)
            label.pack(side="left")
            timing = self.ttk.Label(line, style="Muted.TLabel")
            timing.pack(side="right")
            self.stage_rows.append({"line": line, "marker": marker, "label": label, "timing": timing,
                                    "state": "pending", "started": None})

    def mark_stage(self, index):
        now = time.monotonic()
        for i, row in enumerate(self.stage_rows):
            if i < index and row["state"] != "done":
                row["state"] = "done"
                row["marker"].configure(text="✓")
                if row["started"] is not None:
                    row["timing"].configure(text=format_duration(now - row["started"]))
            elif i == index and row["state"] == "pending":
                row["state"] = "active"
                row["started"] = now
                row["marker"].configure(text="…")
        self.stage_index = max(self.stage_index, index)

    def finish_stages(self):
        self.mark_stage(len(self.stage_rows))

    def append_log(self, message):
        self.log.configure(state="normal")
        self.log.insert("end", message + "\n")
        self.log.see("end")
        self.log.configure(state="disabled")

    def on_progress(self, message):
        self.events.put(("stage", message))

    def handle_stage(self, message):
        self.append_log(message)
        if message == VERIFY_MESSAGE and self.lab_mode in {"small", "public"}:
            self.mark_stage(len(STAGE_LABELS))
        elif message in MESSAGE_TO_STAGE:
            self.mark_stage(MESSAGE_TO_STAGE[message])
        elif self.stage_index < 0:
            self.mark_stage(0)

    def start_run(self):
        if self.busy or self.preparing or self.export_busy:
            return
        if not self.run_preflight():
            return
        self.failure.pack_forget()
        parent = Path(self.output.get()).expanduser()
        cache = Path(self.cache.get()).expanduser()
        root = parent / ("BES-" + datetime.now().strftime("%Y%m%d-%H%M%S-%f"))
        self.pending_root = root
        self.busy, self.started = True, time.monotonic()
        self.start.configure(state="disabled")
        self.back.configure(state="disabled")
        self.change_output_button.configure(state="disabled")
        self.run_title.configure(text="Analysis running")
        self.build_stages()
        self.mark_stage(0)
        mode = self.lab_mode
        kwargs = dict(self.own_kwargs) if mode == "own" else {}
        labels = dict(self.run_labels)
        def worker():
            try:
                original_vcf = lab.pin(kwargs["vcf"]) if kwargs.get("vcf") else None
                if mode == "small-reference":
                    value = lab_workflow.run_reference_example(root, self.backend, progress=self.on_progress)
                else:
                    value = lab.run_job(mode, root, self.backend, cache=cache, progress=self.on_progress, **kwargs)
                try:
                    lab_workflow.save_context(root, labels, original_vcf)
                except Exception as exc:
                    self.events.put(("notice", "The analysis completed, but optional context could not be saved: " + str(exc)))
                self.events.put(("done", value))
            except Exception as exc:
                self.events.put(("error", (str(exc), traceback.format_exc())))
        threading.Thread(target=worker, daemon=False).start()

    def poll(self):
        try:
            while True:
                kind, value = self.events.get_nowait()
                if kind == "notice":
                    self.messages.showinfo("Analysis context", value)
                    continue
                if kind.startswith("prepare_"):
                    self.handle_preparation(kind, value)
                    continue
                if kind.startswith("export_"):
                    self.handle_export(kind, value)
                    continue
                if kind == "table":
                    self.accept_table(value)
                    continue
                if kind == "stage":
                    if self.busy:
                        self.handle_stage(value)
                    continue
                self.busy = False
                self.start.configure(state="normal")
                self.back.configure(state="normal")
                self.change_output_button.configure(state="normal")
                if kind == "error":
                    message, details = value
                    self.show_failure(message, details)
                else:
                    self.result, state = value
                    self.last_run = (self.result, state)
                    self.finish_stages()
                    self.show("results")
        except queue.Empty:
            pass
        if self.busy:
            self.elapsed.configure(text=f"Elapsed  {format_duration(time.monotonic() - self.started)}   ·   Keep this window open")
        self.root.after(200, self.poll)

    def show_failure(self, message, details):
        self.hide_run_details()
        self.check_summary.pack_forget()
        self.run_title.configure(text="Run needs attention")
        self.failure_text.configure(text=message)
        partial = self.pending_root if self.pending_root is not None and self.pending_root.exists() else None
        self.failure_partial.configure(
            text="" if partial is None else
            f"The incomplete run folder is kept for inspection: {partial}. These files do not establish search completion.")
        self.details.configure(state="normal")
        self.details.delete("1.0", "end")
        self.details.insert("1.0", details)
        self.details.configure(state="disabled")
        self.details_disclosure.section.pack(side="bottom", fill="x", before=self.log.master)
        self.failure.pack(fill="x", pady=6, **({"before": self.stage_frame} if self.stage_frame.winfo_manager() else {}))

    # ------------------------------------------------------------------ results

    def info_button(self, parent, command, *, label="Selected row details", state="disabled", hint_to_right=False):
        if not hasattr(self, "row_info_icons"):
            data = """
                #define info_width 20
                #define info_height 20
                static unsigned char info_bits[] = {
                    0x00,0x00,0x00, 0x00,0x06,0x00, 0xc0,0x3f,0x00, 0x70,0xe0,0x00,
                    0x18,0x80,0x01, 0x08,0x06,0x01, 0x0c,0x06,0x03, 0x04,0x00,0x02,
                    0x04,0x00,0x02, 0x06,0x06,0x06, 0x06,0x06,0x06, 0x04,0x06,0x02,
                    0x04,0x06,0x02, 0x0c,0x06,0x03, 0x08,0x06,0x01, 0x18,0x80,0x01,
                    0x70,0xe0,0x00, 0xc0,0x3f,0x00, 0x00,0x06,0x00, 0x00,0x00,0x00};
            """
            self.row_info_icons = [self.tk.BitmapImage(master=self.root, data=data, foreground=color)
                                   for color in (ACCENT, "#b7c4d3")]
        button = self.ttk.Button(parent, text=label, compound="image",
                                 image=(self.row_info_icons[0], "disabled", self.row_info_icons[1]),
                                 style="Info.Quiet.TButton", takefocus=True, state=state, command=command)
        try:
            self.root.tk.call("tk", "accessible", "set_acc_name", button, label)
        except self.tk.TclError:
            pass  # Older Tk retains the button's text label alongside its image.
        hint = self.ttk.Label(parent, style="Caption.TLabel")
        def reveal(event=None):
            hint.configure(text="Select a row for details" if button.instate(["disabled"]) else label)
            hint.place(in_=button, relx=1 if hint_to_right else 0, rely=.5,
                       x=8 if hint_to_right else -8, anchor="w" if hint_to_right else "e")
            hint.lift()
        def dismiss(event=None):
            hint.place_forget()
        button.bind("<Enter>", reveal)
        button.bind("<FocusIn>", reveal)
        button.bind("<Leave>", lambda _: None if button.instate(["focus"]) else dismiss())
        button.bind("<FocusOut>", dismiss)
        button.bind("<Escape>", dismiss)
        button.bind("<Return>", lambda _: button.invoke())
        return button

    def build_results(self):
        ttk = self.ttk
        results, state = self.last_run
        self.results_dir = Path(results)
        def read(name, default):
            path = self.results_dir / name
            return json.loads(path.read_text()) if path.is_file() else default
        self.manifest = read("run-manifest.json", state)
        self.coverage = read("coverage.json", {})
        self.request = read("request.json", {})
        manifest, coverage, request = self.manifest, self.coverage, self.request
        self.guide_sequences = lab_results.guide_sequences(self.results_dir, request)
        self.guide_counts = {guide: 0 for guide in self.guide_sequences}
        self.guide_counts_ready = False
        self.guide_values = {"All guides": "", **{guide: guide for guide in self.guide_sequences}}
        guide_count = coverage.get("guide_count", len(self.guide_sequences))
        self.guide_filter.set(next(iter(self.guide_sequences)) if len(self.guide_sequences) == 1 else "All guides")
        self.table_source.set("Candidates")
        self.result_view = self.tk.StringVar(value="Overview" if guide_count > 1 else "Candidates")
        self.saved_text, self.saved_explanation = lab_results.saved_status(self.results_dir)
        frame = ttk.Frame(self.container)
        frame.rowconfigure(2, weight=1)
        frame.columnconfigure(0, weight=1)

        heading = ttk.Frame(frame)
        heading.grid(row=0, column=0, sticky="ew", pady=(4, 6))
        heading.columnconfigure(0, weight=1)
        sample = (request.get("sample") or {}).get("id", "Reference-only analysis")
        self.analysis_title = ttk.Label(heading, text=f"{sample} · {guide_count} {'guide' if guide_count == 1 else 'guides'} · "
                                        f"{request.get('assembly', 'Unknown assembly')}", style="H2.TLabel")
        self.analysis_title.grid(row=0, column=0, sticky="ew")
        self.analysis_title.bind("<Configure>", lambda event: self.analysis_title.configure(wraplength=max(240, event.width)))
        verification = manifest.get("example_verification", {})
        self.completion_text = "Example verification failed" if verification.get("status") == "mismatched" else result_status(manifest)

        tabs = ttk.Frame(frame)
        tabs.grid(row=1, column=0, sticky="ew", pady=(4, 10))
        self.result_tabs = {}
        for label in ("Overview", "Candidates", "Sample differences"):
            if label == "Sample differences" and not request.get("sample"):
                continue
            button = ttk.Radiobutton(tabs, text=label, variable=self.result_view, value=label,
                                      style="Tab.TRadiobutton", width={"Overview": 12, "Candidates": 14, "Sample differences": 22}[label],
                                      command=self.select_results_view)
            button.pack(side="left", padx=(0, 4))
            button.bind("<Return>", lambda _, widget=button: widget.invoke())
            self.result_tabs[label] = button
        content = ttk.Frame(frame)
        content.grid(row=2, column=0, sticky="nsew")
        self.overview_frame = ttk.Frame(content)
        self.build_guide_overview(self.overview_frame)
        self.candidate_frame = table_frame = ttk.Frame(content)
        self.results_filters = Disclosure(table_frame, "Search & filters", ttk)
        self.results_filters.section.pack_configure(pady=(0, 8))
        self.filters_button = self.results_filters.button
        self.filters_button.pack_configure(pady=(0, 2))
        self.filters_button.configure(command=self.toggle_result_filters)
        self.filters_button.bind("<Return>", lambda _: self.filters_button.invoke())
        self.filter_panel = self.results_filters.body
        controls = ttk.Frame(self.filter_panel)
        controls.pack(fill="x", pady=(8, 0))
        controls.columnconfigure(3, weight=1)
        ttk.Label(controls, text="Guide").grid(row=0, column=0, padx=(0, 6))
        self.guide_box = ttk.Combobox(controls, textvariable=self.guide_filter, state="readonly", width=19, values=list(self.guide_values))
        self.guide_box.grid(row=0, column=1, sticky="w", padx=(0, 16))
        self.guide_box.bind("<<ComboboxSelected>>", lambda _: self.refresh_table())
        ttk.Label(controls, text="Search").grid(row=0, column=2, padx=(0, 6))
        self.search = ttk.Entry(controls, textvariable=self.table_query, width=15)
        self.search.grid(row=0, column=3, sticky="ew")
        self.search.bind("<Return>", lambda _: self.refresh_table())
        self.search.bind("<Escape>", lambda _: self.toggle_result_filters())
        self.apply_button = ttk.Button(controls, text="Apply", command=self.refresh_table)
        self.apply_button.grid(row=0, column=4, padx=(6, 4))
        options = ttk.Frame(self.filter_panel)
        options.pack(fill="x", pady=(8, 0))
        ttk.Label(options, text="Show").pack(side="left", padx=(0, 8))
        self.filter_box = ttk.Combobox(options, textvariable=self.filter_display, state="readonly", width=30)
        self.filter_box.pack(side="left")
        self.filter_box.bind("<<ComboboxSelected>>", lambda _: self.apply_table_filter())
        self.filter_anchor = ttk.Frame(self.filter_panel)
        self.filter_anchor.pack(fill="x", pady=(4, 0))
        self.reset_filters_button = ttk.Button(self.filter_anchor, text="Reset", style="Quiet.TButton", padding=(8, 3), command=self.clear_table_search)
        self.reset_filters_button.pack(side="right")
        self.active_filters = ttk.Label(self.filter_anchor, style="Caption.TLabel", justify="left")
        self.active_filters.pack(side="left", fill="x", expand=True)
        self.active_filters.bind("<Configure>", lambda event: self.active_filters.configure(wraplength=max(180, event.width)))
        body = ttk.Frame(table_frame)
        self.table = ttk.Treeview(body, show="headings", selectmode="extended", height=4)
        yscroll = ttk.Scrollbar(body, orient="vertical", command=self.table.yview)
        xscroll = ttk.Scrollbar(body, orient="horizontal", command=self.table.xview)
        body.rowconfigure(0, weight=1)
        body.columnconfigure(0, weight=1)
        self.table.configure(yscrollcommand=lambda a, b: self.set_scrollbar(yscroll, a, b),
                             xscrollcommand=lambda a, b: self.set_scrollbar(xscroll, a, b))
        yscroll.grid(row=0, column=1, sticky="ns")
        xscroll.grid(row=1, column=0, sticky="ew")
        self.table.grid(row=0, column=0, sticky="nsew")
        self.table.tag_configure("stripe", background="#fafbfc")
        self.empty_table = ttk.Label(body, style="CardMuted.TLabel", justify="center")
        self.table.bind("<Double-1>", self.show_row_detail)
        self.table.bind("<Return>", self.show_row_detail)
        self.table.bind("<<TreeviewSelect>>", lambda _: self.row_button.configure(
            state="normal" if self.table.selection() else "disabled"))
        foot = ttk.Frame(table_frame)
        foot.pack(side="bottom", fill="x", pady=(6, 0))
        body.pack(fill="both", expand=True)
        self.row_button = self.info_button(foot, self.show_row_detail)
        self.row_button.pack(side="right")
        self.table_count = ttk.Label(foot, style="Caption.TLabel", justify="left")
        self.table_count.pack(side="left", fill="x", expand=True)
        self.table_count.bind("<Configure>", lambda event: self.table_count.configure(wraplength=max(180, event.width)))

        actions = ttk.Frame(frame)
        actions.grid(row=3, column=0, sticky="ew", pady=(8, 12))
        self.export_button = ttk.Button(actions, text="Export…", style="Primary.TButton", command=self.show_export)
        self.export_button.pack(side="right")
        self.more_button = ttk.Menubutton(actions, text="Learn more", takefocus=True)
        self.more_menu = self.tk.Menu(self.more_button, tearoff=False)
        self.more_menu.add_command(label="Coverage", command=lambda: self.show_result_details("Coverage"))
        self.more_menu.add_command(label="Run details", command=lambda: self.show_result_details("Run details"))
        self.more_menu.add_command(label="Full report", command=lambda: self.open_path(self.results_dir / "report.html"))
        self.more_menu.add_command(label="Open folder", command=lambda: self.open_path(self.results_dir))
        self.more_menu.add_separator()
        self.more_menu.add_command(label="New analysis", command=lambda: self.show("inputs"))
        if not request.get("sample"):
            self.more_menu.add_command(label="Add sample…", command=self.add_sample)
        self.more_button.configure(menu=self.more_menu)
        self.more_button.pack(side="left")
        self.replace_screen("results", frame)
        self.select_table()
        self.select_results_view()

    def build_guide_overview(self, pane):
        ttk = self.ttk
        body = ttk.Frame(pane)
        body.columnconfigure(0, weight=1)
        body.rowconfigure(0, weight=1)
        tree = self.guide_tree = ttk.Treeview(body, columns=("guide", "sequence", "count"), show="headings", selectmode="browse", height=4)
        for key, label, width in (("guide", "Guide", 140), ("sequence", "Input spacer", 310), ("count", "Candidate sequences", 180)):
            tree.heading(key, text=label, anchor="e" if key == "count" else "w")
            tree.column(key, width=width, minwidth=width, stretch=key == "sequence", anchor="e" if key == "count" else "w")
        scroll = ttk.Scrollbar(body, command=tree.yview)
        scroll.grid(row=0, column=1, sticky="ns")
        tree.configure(yscrollcommand=lambda a, b: self.set_scrollbar(scroll, a, b))
        tree.grid(row=0, column=0, sticky="nsew")
        self.overview_note = self.caption(pane, "Reading all candidate records…")
        self.overview_note.pack_forget()
        footer = ttk.Frame(pane)
        footer.pack(side="bottom", fill="x", pady=(8, 0))
        self.overview_note.pack(side="bottom", fill="x", pady=(8, 0))
        body.pack(fill="both", expand=True, pady=(4, 0))
        self.all_guides_button = ttk.Button(footer, text="All guides", command=self.inspect_all_guides)
        self.all_guides_button.pack(side="left")
        self.inspect_guide_button = ttk.Button(footer, text="Inspect candidates →", state="disabled", command=self.inspect_guide)
        self.inspect_guide_button.pack(side="right")
        tree.bind("<<TreeviewSelect>>", lambda _: self.inspect_guide_button.configure(state="normal" if tree.selection() else "disabled"))
        tree.bind("<Double-1>", self.inspect_guide)
        tree.bind("<Return>", self.inspect_guide)

    def fill_guide_overview(self, error=None):
        tree = self.guide_tree
        selected = tree.item(tree.selection()[0], "values")[0] if tree.selection() else None
        tree.delete(*tree.get_children())
        if error:
            self.overview_note.configure(text="Candidate counts unavailable. " + error)
            self.overview_note.pack(side="bottom", fill="x", pady=(8, 0), before=tree.master)
            self.inspect_guide_button.configure(state="disabled")
            return
        for guide, count in self.guide_counts.items():
            item = tree.insert("", "end", values=(guide, self.guide_sequences.get(guide, "Not available"), f"{count:,}"))
            if guide == selected:
                tree.selection_set(item)
                tree.focus(item)
        if not tree.selection() and tree.get_children():
            tree.selection_set(tree.get_children()[0])
            tree.focus(tree.get_children()[0])
        self.inspect_guide_button.configure(state="normal" if tree.selection() else "disabled")
        self.overview_note.pack_forget()
        if not self.guide_sequences:
            self.overview_note.configure(text="Input guide file unavailable; only guides with candidates are listed.")
            self.overview_note.pack(side="bottom", fill="x", pady=(8, 0), before=tree.master)

    def select_results_view(self):
        view = self.result_view.get()
        if view == "Overview":
            self.candidate_frame.pack_forget()
            self.overview_frame.pack(fill="both", expand=True)
            if not self.guide_counts_ready and self.table_source.get() != "Candidates":
                self.table_source.set("Candidates")
                self.select_table()
        else:
            self.overview_frame.pack_forget()
            self.candidate_frame.pack(fill="both", expand=True)
            if self.table_source.get() != view:
                self.table_source.set(view)
                self.select_table()

    def show_guides(self):
        self.result_view.set("Overview")
        self.select_results_view()
        self.guide_tree.focus_set()

    def inspect_guide(self, event=None):
        if not self.guide_tree.selection():
            return
        guide = self.guide_tree.item(self.guide_tree.selection()[0], "values")[0]
        self.guide_filter.set(next(key for key, value in self.guide_values.items() if value == guide))
        self.result_view.set("Candidates")
        self.table_source.set("Candidates")
        self.select_results_view()
        self.select_table()
        self.table.focus_set()

    def inspect_all_guides(self):
        self.guide_filter.set("All guides")
        self.result_view.set("Candidates")
        self.table_source.set("Candidates")
        self.select_results_view()
        self.select_table()

    def toggle_result_filters(self):
        self.results_filters.toggle()
        if not self.results_filters.open:
            self.filters_button.focus_set()
        else:
            self.search.focus_set()

    def add_sample(self):
        if self.request.get("input_kind") == "artificial":
            self.flow.set("example")
            self.example_kind.set("small")
        else:
            self.flow.set("own")
            if self.guide_sequences:
                self.guides.delete("1.0", "end")
                self.guides.insert("1.0", "guide_id\tspacer_sequence\n" + "\n".join(
                    f"{guide}\t{sequence}" for guide, sequence in self.guide_sequences.items()))
            self.reference.set(self.request.get("reference", {}).get("path", ""))
            self.clear_vcf()
            if not self.sample_disclosure.open:
                self.sample_disclosure.toggle()
        self.update_flow()
        self.show("inputs")
        self.input_error.configure(foreground=MUTED)
        self.set_input_error(self.input_error, "Adding a sample starts a new analysis. Your original results stay in their existing folder.")

    def show_export(self):
        if self.export_busy:
            return
        if getattr(self, "table_loading", False):
            self.messages.showinfo("Reading results", "Wait for the table to finish loading, then export.")
            return
        ttk = self.ttk
        popup = self.tk.Toplevel(self.root)
        popup.title("Export")
        popup.geometry("780x640")
        popup.minsize(700, 530)
        outer = ttk.Frame(popup, padding=24)
        outer.pack(fill="both", expand=True)
        footer = ttk.Frame(outer)
        footer.pack(side="bottom", fill="x", pady=(12, 0))
        body = ttk.Frame(outer)
        body.pack(fill="both", expand=True)
        pane = self.scrollable(body)
        ttk.Label(pane, text="Export results", style="H2.TLabel").pack(anchor="w")
        self.caption(pane, "Create a separate export. Your saved analysis stays in its original folder.")
        kind = self.tk.StringVar(value="candidates")
        formats = ttk.Frame(pane)
        formats.pack(fill="x", pady=(8, 8))
        candidate_options = ttk.Frame(pane)
        candidate_options.pack(fill="x")
        def change_kind():
            if kind.get() == "candidates":
                candidate_options.pack(fill="x", after=formats)
            else:
                candidate_options.pack_forget()
        for text, value in (("Candidate tables (TSV)", "candidates"), ("Report (HTML)", "report"), ("Analysis copy", "analysis")):
            ttk.Radiobutton(formats, text=text, value=value, variable=kind, command=change_kind).pack(side="left", padx=(0, 14))
        applied = self.applied_table
        candidate_view = bool(applied and applied["source"] == "Candidates" and self.result_view.get() == "Candidates")
        selected = [self.row_records[item]["configuration_id"] for item in self.table.selection()
                    if item in self.row_records] if candidate_view else []
        filtered = candidate_view and bool(applied["query"].strip() or applied["filters"])
        mode = self.tk.StringVar(value="filtered" if filtered else "all")
        choices = [(f"All candidates · {self.manifest.get('candidate_configurations', 0):,} sequences", "all", True),
                   (f"Current matches · {self.table_matches:,} sequences" if candidate_view else "Current matches · open Candidates to filter", "filtered", candidate_view),
                   (f"Selected rows · {len(selected):,} sequences", "selected", bool(selected))]
        for text, value, enabled in choices:
            ttk.Radiobutton(candidate_options, text=text, value=value, variable=mode, state="normal" if enabled else "disabled").pack(anchor="w", pady=(12, 0))
        criteria = ("Search: " + (applied["query"].strip() or "none") + " · " +
                    "; ".join(f"{key} = {value}" for key, value in applied["filters"].items())) if candidate_view else "All candidates from the complete result bundle."
        self.caption(candidate_options, criteria.rstrip(" ·"))
        self.caption(candidate_options, "All candidates and current matches include rows beyond the 5,000-row table display limit.")
        destination_parent = self.tk.StringVar(value=str(lab_results.analysis_folder(self.results_dir).parent))
        self.path_control(pane, "Save export in", destination_parent, "output")
        self.caption(pane, "Creates a new dated export folder. The complete analysis is preserved.")
        details = Disclosure(pane, "What is included", ttk)
        self.caption(details.body, "Candidates and their original locus comparisons; full-run coverage and exclusions; "
                     "source request and manifest; exact selection criteria; READ-ME.txt field dictionary. "
                     "The HTML report is copied as written. Analysis copy includes existing run files and prepared inputs; external reference files stay at their recorded paths.")
        status = ttk.Label(footer, wraplength=700, justify="left", style="Muted.TLabel")
        status.pack(fill="x", pady=(0, 12))
        status.bind("<Configure>", lambda event: status.configure(wraplength=max(240, event.width)))
        export_actions = ttk.Frame(footer)
        export_actions.pack(fill="x")
        action = ttk.Button(export_actions, text="Export", style="Primary.TButton")
        action.pack(side="right")
        self.export_ui = (popup, status, action)
        source = self.results_dir
        def close():
            if not self.export_busy:
                popup.destroy()
        cancel = self.export_cancel = ttk.Button(export_actions, text="Cancel", command=close)
        cancel.pack(side="right", padx=8)
        popup.protocol("WM_DELETE_WINDOW", close)
        popup.bind("<Escape>", lambda _: close())
        def start():
            if self.export_busy:
                return
            destination = Path(destination_parent.get()).expanduser() / ("BES-export-" + datetime.now().strftime("%Y%m%d-%H%M%S-%f"))
            selected_mode = mode.get()
            query = applied["query"] if candidate_view else ""
            filters = dict(applied["filters"]) if candidate_view else {}
            self.export_busy = True
            action.configure(state="disabled")
            cancel.configure(state="disabled")
            status.configure(text="Writing the export and checking source files…")
            selected_kind = kind.get()
            def worker():
                try:
                    if selected_kind == "candidates":
                        receipt = lab_results.export_candidates(source, destination, mode=selected_mode,
                            query=query, filters=filters, selected_ids=selected)
                    else:
                        receipt = lab_results.export_copy(source, destination, selected_kind)
                    self.events.put(("export_done", (destination, receipt)))
                except Exception as exc:
                    self.events.put(("export_error", str(exc)))
            threading.Thread(target=worker, daemon=False).start()
        action.configure(command=start)
        popup.transient(self.root)
        popup.grab_set()
        popup.focus_set()

    def handle_export(self, kind, value):
        popup, status, action = self.export_ui
        self.export_busy = False
        action.configure(state="normal")
        self.export_cancel.configure(state="normal", text="Close")
        if kind == "export_done":
            path, receipt = value
            self.last_export = path
            message = (f"Exported {receipt['candidate_rows']:,} candidate sequences." if "candidate_rows" in receipt
                       else "HTML report exported." if receipt["kind"] == "report" else "Analysis copied.")
            status.configure(text=f"{message}\n{path}")
            action.configure(text="Open export folder", command=lambda: self.open_path(path))
        else:
            status.configure(text="Export did not finish. " + value)

    def detail_text(self, parent, content, height=16):
        area = self.ttk.Frame(parent)
        area.pack(fill="both", expand=True, pady=(8, 0))
        text = self.tk.Text(area, wrap="word", height=height, font=(MONO, 11),
                            background=CARD, foreground=TEXT, borderwidth=0, padx=14, pady=12)
        scroll = self.ttk.Scrollbar(area, command=text.yview)
        text.configure(yscrollcommand=scroll.set)
        scroll.pack(side="right", fill="y")
        text.pack(fill="both", expand=True)
        text.insert("1.0", content)
        text.configure(state="disabled")
        return text

    def show_result_details(self, selected="Coverage"):
        ttk = self.ttk
        popup = self.tk.Toplevel(self.root)
        popup.title(selected)
        popup.geometry("980x700")
        popup.minsize(760, 560)
        popup.configure(background=BG)
        popup.bind("<Escape>", lambda _: popup.destroy())
        outer = ttk.Frame(popup, padding=20)
        outer.pack(fill="both", expand=True)
        ttk.Button(outer, text="Close", command=popup.destroy).pack(side="bottom", anchor="e", pady=(10, 0))
        if selected == "Run details":
            self.build_run_details(outer)
        else:
            notebook = ttk.Notebook(outer)
            notebook.pack(fill="both", expand=True)
            summary, exclusions = ttk.Frame(notebook), ttk.Frame(notebook, padding=12)
            notebook.add(summary, text="Summary")
            notebook.add(exclusions, text="Excluded regions")
            cover = self.scrollable(summary)
            self.caption(cover, coverage_brief(self.coverage, self.request))
            for label, value in coverage_facts(self.coverage, self.request):
                self.summary_row(cover, label, value)
            self.caption(cover, "Windows: 23 bases, both strands, per guide, before PAM or mismatch filtering.")
            self.summary_row(cover, "Reference callability", self.coverage.get("reference_callability", "not_assessed").replace("_", " "))
            if self.request.get("sample"):
                self.summary_row(cover, "Sample model basis", "Conditional on retained VCF calls; absent records do not establish callability or confident genotypes.")
                self.summary_row(cover, "Sample analysis status", self.coverage.get("sample_status", "Not recorded").replace("_", " "))
            ttk.Label(cover, text="By chromosome", style="SetupField.TLabel").pack(anchor="w", pady=(12, 6))
            body = ttk.Frame(cover)
            body.pack(fill="x")
            body.columnconfigure(0, weight=1)
            sample = bool(self.request.get("sample"))
            columns = ("contig", "bases", "windows", "copies") if sample else ("contig", "windows")
            tree = ttk.Treeview(body, columns=columns, show="headings", height=4)
            titles = ("Chromosome", "Modeled / reference bases", "Comparable / possible windows", "Sample copy windows") if sample else ("Chromosome", "Evaluable / possible reference windows")
            widths = (110, 250, 265, 175) if sample else (160, 540)
            for column, title, width in zip(columns, titles, widths):
                tree.heading(column, text=title)
                tree.column(column, width=width, minwidth=width)
            scroll = ttk.Scrollbar(body, command=tree.yview)
            horizontal = ttk.Scrollbar(body, orient="horizontal", command=tree.xview)
            tree.configure(yscrollcommand=lambda a, b: self.set_scrollbar(scroll, a, b),
                           xscrollcommand=lambda a, b: self.set_scrollbar(horizontal, a, b))
            scroll.grid(row=0, column=1, sticky="ns")
            horizontal.grid(row=1, column=0, sticky="ew")
            tree.grid(row=0, column=0, sticky="nsew")
            for row in self.coverage.get("chromosomes", []):
                pair = lambda a, b: f"{row.get(a, 0):,} / {row.get(b, 0):,}"
                values = (row["contig"], pair("modeled_reference_bases", "reference_bases"),
                          pair("comparable_oriented_windows_per_guide", "reference_oriented_windows_per_guide"),
                          f"{row.get('sample_evaluated_copy_oriented_windows_per_guide', 0):,}") if sample else (
                          row["contig"], pair("reference_evaluable_oriented_windows_per_guide", "reference_oriented_windows_per_guide")
                          if "reference_evaluable_oriented_windows_per_guide" in row else "Not recorded")
                tree.insert("", "end", values=values)
            if not self.coverage.get("chromosomes"):
                self.caption(cover, "Chromosome accounting is not recorded in this bundle.")
            reasons = self.coverage.get("exclusion_summary", {})
            if reasons:
                ttk.Label(cover, text="Exclusion records", style="SetupField.TLabel").pack(anchor="w", pady=(12, 4))
                self.caption(cover, "Records and reasons can overlap. Open Excluded regions for original row details.")
                for reason, count in reasons.items():
                    self.summary_row(cover, reason_label(reason), f"{count:,}")
            assumptions = Disclosure(cover, "Modeling details and complete accounting", ttk)
            if sample:
                self.caption(assumptions.body, "Local copy labels do not establish long-range phase. Missing calls and incompatible phase remain unresolved.")
            self.detail_text(assumptions.body, json.dumps(self.coverage, indent=2), height=10)
            self.build_exclusions(exclusions, popup)
            self.coverage_ui = {"popup": popup, "notebook": notebook, "summary": cover, "chromosomes": tree}
        popup.transient(self.root)
        popup.focus_set()

    def build_run_details(self, details):
        source = self.manifest.get("source", {})
        counts = self.manifest
        facts = [
            ("Analysis", lab_results.analysis_folder(self.results_dir).name),
            ("Analysis status", self.completion_text),
            ("Candidate sequences", f"{counts.get('candidate_configurations', 0):,} sequence configurations; multiple sequences can share a locus."),
            ("Reference loci", f"{counts.get('reference_loci', 0):,} reference candidate locations."),
            ("Sample paths", f"{counts.get('sample_loci', 0):,} sample candidate paths, including indel paths; a different counting unit." if self.request.get("sample") else "Not requested."),
            ("Locus comparisons", json.dumps(counts.get("locus_changes", {}), sort_keys=True)),
            ("Retained loci with changed sequences", f"{counts.get('changed_sequence_loci', 0):,}; not locus gains or losses."),
            ("Package version", source.get("package_version", "—")),
            ("Source SHA-256", source.get("source_sha256", "—")),
            ("Python", source.get("python", "—")), ("Run ID", self.manifest.get("run_id", "—")),
            ("Started UTC", self.manifest.get("started_utc", "—")),
            ("Finished UTC", self.manifest.get("finished_utc", "—")),
            ("Wall time", format_duration(self.manifest.get("wall_seconds", 0))),
            ("Peak parent memory", format_bytes(self.manifest.get("parent_peak_resident_bytes", 0))),
            ("Coverage basis", self.manifest.get("coverage_basis", "—")),
            ("Sample status", self.manifest.get("genotype_status", "—"))]
        self.caption(details, self.saved_text + ": " + self.saved_explanation)
        self.caption(details, "Search settings: 20-base guides · NGG · both strands · up to "
                     f"{self.request.get('max_mismatches', 6)} mismatches · no gaps")
        footer = self.ttk.Frame(details)
        footer.pack(side="bottom", fill="x", pady=(8, 0))
        command = self.reproduction_command(self.results_dir)
        self.ttk.Button(footer, text="Copy reproduction command", command=lambda: self.copy_command(command)).pack(anchor="w")
        self.caption(footer, "Replace NEW_OUTPUT_FOLDER with a new folder before running the command.")
        self.detail_text(details, "\n".join(f"{label}: {value}" for label, value in facts) +
                         "\n\nReading the tables\nGuide counts are sequence configurations. "
                         "0 means no candidate within the search scope; Coverage identifies unmodeled or unsearched regions.\n" +
                         "\n".join(f"{label}: {note}" for label, note in VIEW_NOTES.items()) +
                         "\n\nFull input identities and parameters\n" + json.dumps(self.request, indent=2) +
                         "\n\nRun metadata and checksums\n" + json.dumps(self.manifest, indent=2) +
                         "\n\nOptional experiment context\n" + json.dumps(lab_results.analysis_context(self.results_dir), indent=2))

    def build_exclusions(self, pane, popup):
        ttk = self.ttk
        self.caption(pane, "Excluded regions and phase notes · 0-based, end excluded · records may overlap")
        controls = ttk.Frame(pane)
        controls.pack(fill="x", pady=(4, 8))
        ttk.Label(controls, text="Search").pack(side="left", padx=(0, 8))
        query = self.tk.StringVar()
        entry = ttk.Entry(controls, textvariable=query)
        entry.pack(side="left", fill="x", expand=True)
        apply = ttk.Button(controls, text="Apply")
        apply.pack(side="left", padx=6)
        reset = ttk.Button(controls, text="Reset", style="Quiet.TButton")
        reset.pack(side="left")
        footer = ttk.Frame(pane)
        footer.pack(side="bottom", fill="x", pady=(8, 0))
        details = self.info_button(footer, None)
        details.pack(side="right")
        count = ttk.Label(footer, style="Caption.TLabel", justify="left")
        count.pack(side="left", fill="x", expand=True)
        count.bind("<Configure>", lambda event: count.configure(wraplength=max(180, event.width)))
        body = ttk.Frame(pane)
        body.pack(fill="both", expand=True)
        body.rowconfigure(0, weight=1)
        body.columnconfigure(0, weight=1)
        tree = ttk.Treeview(body, columns=VIEW_COLUMNS["Excluded regions"], show="headings", selectmode="browse", height=5)
        for column in VIEW_COLUMNS["Excluded regions"]:
            label, width = COLUMN_SPECS[column]
            tree.heading(column, text=label)
            tree.column(column, width=width, minwidth=width)
        scroll = ttk.Scrollbar(body, command=tree.yview)
        scroll.grid(row=0, column=1, sticky="ns")
        tree.configure(yscrollcommand=lambda a, b: self.set_scrollbar(scroll, a, b))
        tree.grid(row=0, column=0, sticky="nsew")
        empty = ttk.Label(body, style="Muted.TLabel", justify="center")
        records, pending = {}, queue.Queue()
        path = self.results_dir / "exclusions.tsv"
        generation = 0
        poll_id = None
        def inspect(event=None):
            if tree.selection() and tree.selection()[0] in records:
                self.show_row_detail(row=records[tree.selection()[0]], source="Excluded regions")
        def refresh():
            nonlocal generation
            generation += 1
            current, search = generation, query.get()
            tree.delete(*tree.get_children())
            records.clear()
            details.configure(state="disabled")
            count.configure(text="Reading exclusions…")
            empty.configure(text="Reading exclusions…")
            empty.place(relx=.5, rely=.5, anchor="center")
            def worker():
                try:
                    result = load_rows(path, search) if path.is_file() else ([], 0, 0)
                    pending.put((current, result, None))
                except Exception as exc:
                    pending.put((current, ([], 0, 0), str(exc)))
            threading.Thread(target=worker, daemon=True).start()
        def accept():
            nonlocal poll_id
            if not popup.winfo_exists():
                return
            try:
                while True:
                    current, (rows, matched, total), error = pending.get_nowait()
                    if current != generation:
                        continue
                    empty.place_forget()
                    for row in rows:
                        item = tree.insert("", "end", values=[display_value(row, c) for c in VIEW_COLUMNS["Excluded regions"]])
                        records[item] = row
                    count.configure(text=row_count_text(len(rows), matched, total).replace("sort shown rows only", "display limit 5,000"))
                    if not rows or error:
                        message = ("Could not read exclusions: " + error if error else "Exclusion table not recorded in this bundle." if not path.is_file()
                                   else "No matching exclusions." if total else "No exclusion records in this bundle.")
                        empty.configure(text=message)
                        empty.place(relx=.5, rely=.5, anchor="center")
            except queue.Empty:
                pass
            poll_id = popup.after(100, accept)
        def clear():
            query.set("")
            refresh()
        details.configure(command=inspect)
        apply.configure(command=refresh)
        reset.configure(command=clear)
        entry.bind("<Return>", lambda _: refresh())
        tree.bind("<Return>", inspect)
        tree.bind("<Double-1>", inspect)
        tree.bind("<<TreeviewSelect>>", lambda _: details.configure(state="normal" if tree.selection() else "disabled"))
        self.exclusion_ui = {"tree": tree, "query": query, "apply": apply, "reset": reset, "details": details, "count": count, "records": records}
        refresh()
        poll_id = popup.after(100, accept)
        popup.bind("<Destroy>", lambda event: popup.after_cancel(poll_id) if event.widget == popup and poll_id else None)

    def summary_row(self, parent, label, value, note=None):
        line = self.ttk.Frame(parent)
        line.pack(fill="x", pady=1)
        self.ttk.Label(line, text=label, width=24, wraplength=225,
                       justify="left", style="Muted.TLabel").pack(side="left", anchor="n", padx=(0, 12))
        text = value + (f"  —  {note}" if note else "")
        value = self.ttk.Label(line, text=text, wraplength=620, justify="left")
        value.pack(side="left", fill="x", expand=True)
        value.bind("<Configure>", lambda event: value.configure(wraplength=max(160, event.width)))

    def reproduction_command(self, results):
        return " ".join(["bes-search", "--config", shlex.quote(str(Path(results) / "request.json")),
                         "--output", "NEW_OUTPUT_FOLDER"])

    def copy_command(self, command):
        self.root.clipboard_clear()
        self.root.clipboard_append(command)
        self.messages.showinfo("Copied", "The reproduction command is on the clipboard. "
                               "Replace NEW_OUTPUT_FOLDER with a new folder before running it.")

    def table_file(self):
        return self.results_dir / dict(TABLE_SOURCES)[self.table_source.get()]

    def table_filters(self):
        label, value = self.table_source.get(), self.table_filter.get()
        filters = {}
        if label != "Excluded regions":
            filters["guide_id"] = self.guide_values.get(self.guide_filter.get(), "")
        if value:
            if label == "Sample differences":
                filters["difference" if value in {"differences", "changed"} else "status"] = value
            else:
                filters["contig"] = value
        return {key: value for key, value in filters.items() if value}

    def select_table(self):
        self.table_filter.set("differences" if self.table_source.get() == "Sample differences" else "")
        self.table_query.set("")
        self.sort_column, self.sort_reverse = None, False
        self.refresh_table()

    def apply_table_filter(self):
        self.table_filter.set(self.filter_values.get(self.filter_display.get(), ""))
        self.refresh_table()

    def refresh_table(self):
        path, label = self.table_file(), self.table_source.get()
        self.guide_box.configure(state="disabled" if label == "Excluded regions" else "readonly")
        column = "status" if label == "Sample differences" else "contig"
        default = "All loci" if column == "status" else "All chromosomes"
        self.set_filter_choices(self.table_cache.get(str(path), []), default)
        query, filters = self.table_query.get(), self.table_filters()
        self.table_generation += 1
        generation = self.table_generation
        self.table_loading = True
        self.export_button.configure(state="disabled")
        self.active_filters.configure(text="Applying selection…")
        self.reset_filters_button.configure(state="disabled")
        self.applied_table = None
        self.table.delete(*self.table.get_children())
        self.row_records = {}
        self.row_button.configure(state="disabled")
        self.empty_table.configure(text="Reading results…")
        self.empty_table.place(relx=.5, rely=.5, anchor="center")
        self.table_count.configure(text="Reading results…")
        def worker():
            try:
                choices, counts = {column: set()}, {}
                rows, matched, total = load_rows(path, query, filters, choices=choices,
                    guide_counts=counts if label == "Candidates" else None) if path.is_file() else ([], 0, 0)
                value = (generation, path, rows, matched, total, sorted(choices[column]), default, None, query, filters, counts)
            except Exception as exc:
                value = (generation, path, [], 0, 0, [], default, str(exc), query, filters, {})
            self.events.put(("table", value))
        threading.Thread(target=worker, daemon=True).start()

    def set_filter_choices(self, values, default):
        self.filter_values = {default: ""}
        if self.table_source.get() == "Sample differences":
            self.filter_values.update({"All differences": "differences", "Changed sequence": "changed"})
        self.filter_values.update({STATUS_LABELS.get(value, value): value for value in values})
        self.filter_box.configure(values=list(self.filter_values))
        self.filter_display.set(next((label for label, value in self.filter_values.items()
                                      if value == self.table_filter.get()), default))

    def accept_table(self, value):
        generation, path, rows, matched, total, choices, default, error, query, filters, counts = value
        if generation != self.table_generation or self.screen != "results":
            return
        self.table_loading = False
        self.export_button.configure(state="normal")
        self.table_cache[str(path)] = choices
        self.set_filter_choices(choices, default)
        self.shown_rows, self.table_matches, self.table_total = rows, matched, total
        if self.table_source.get() == "Candidates" and not error:
            self.guide_counts_ready = True
            self.guide_counts = {guide: counts.get(guide, 0) for guide in self.guide_sequences}
            self.guide_counts.update(counts)
            selected = filters.get("guide_id", "")
            self.guide_values = {"All guides": "", **{f"{guide} · {count:,}": guide for guide, count in self.guide_counts.items()}}
            self.guide_box.configure(values=list(self.guide_values))
            self.guide_filter.set(next((key for key, value in self.guide_values.items() if value == selected), "All guides"))
        if self.table_source.get() == "Candidates":
            self.fill_guide_overview(error)
        self.applied_table = {"source": self.table_source.get(), "query": query, "filters": filters} if not error else None
        criteria = []
        if filters.get("guide_id"):
            criteria.append("Guide: " + filters["guide_id"])
        if query.strip():
            criteria.append("Search: “" + query.strip() + "”")
        if self.table_filter.get():
            criteria.append(self.filter_display.get())
        self.results_filters.title = "Search & filters" + (f" · {len(criteria)} active" if criteria else "")
        self.filters_button.configure(text=self.results_filters._label())
        self.active_filters.configure(text=" · ".join(criteria) or "No filters")
        self.reset_filters_button.configure(state="normal" if criteria else "disabled")
        self.fill_table(rows, matched, total, path)
        if error:
            self.empty_table.configure(text="Could not read this table.\n" + error)
            self.empty_table.place(relx=.5, rely=.5, anchor="center")
            self.table_count.configure(text="Table unavailable; open the results folder to inspect the file.")

    def fill_table(self, rows, matched, total, path):
        self.table.delete(*self.table.get_children())
        columns = VIEW_COLUMNS[self.table_source.get()]
        self.table["columns"] = columns
        self.row_records = {}
        for column in columns:
            title, width = COLUMN_SPECS[column]
            arrow = (" ↓" if self.sort_reverse else " ↑") if column == self.sort_column else ""
            self.table.heading(column, text=title + arrow, anchor="w", command=lambda c=column: self.sort_rows(c))
            self.table.column(column, width=width, minwidth=width, stretch=column in {"location", "reason"})
        if self.sort_column:
            def key(row):
                if self.sort_column == "location":
                    return (row.get("contig", ""), sort_key(row.get("start0", row.get("start", ""))),
                            sort_key(row.get("end0", row.get("end", ""))))
                return sort_key(row.get(self.sort_column, ""))
            rows = sorted(rows, key=key, reverse=self.sort_reverse)
        for i, row in enumerate(rows):
            item = self.table.insert("", "end", values=[display_value(row, c) for c in columns],
                                     tags=("stripe",) if i % 2 else ())
            self.row_records[item] = row
        self.empty_table.place_forget()
        if not rows:
            if not path.is_file():
                message = "This table is not part of the run."
            elif total and not matched:
                message = "No matching rows\nClear the search or choose another filter."
            else:
                message = ("No candidates were nominated within this search scope.\nSee coverage for unresolved or unsearched regions."
                           if self.table_source.get() == "Candidates" else "No records in this view.")
            self.empty_table.configure(text=message)
            self.empty_table.place(relx=.5, rely=.5, anchor="center")
        self.table_count.configure(text=row_count_text(len(rows), matched, total))

    def sort_rows(self, column):
        if getattr(self, "table_loading", False):
            return
        if self.sort_column == column:
            self.sort_reverse = not self.sort_reverse
        else:
            self.sort_column, self.sort_reverse = column, False
        self.fill_table(self.shown_rows, self.table_matches, self.table_total, self.table_file())

    def clear_table_search(self):
        self.table_query.set("")
        self.table_filter.set("differences" if self.table_source.get() == "Sample differences" else "")
        self.guide_filter.set("All guides")
        self.refresh_table()

    def show_row_detail(self, event=None, *, row=None, source=None):
        if row is None:
            selected = self.table.selection()
            if not selected or selected[0] not in self.row_records:
                return
            row = self.row_records[selected[0]]
        source = source or self.table_source.get()
        previous_focus = self.root.focus_get()
        popup = self.tk.Toplevel(self.root)
        popup.title("Row details")
        popup.geometry("820x610")
        popup.minsize(620, 460)
        popup.configure(background=BG)
        pane = self.ttk.Frame(popup, padding=24)
        pane.pack(fill="both", expand=True)
        self.ttk.Label(pane, text=row.get("guide_id", "Excluded region"), style="H2.TLabel").pack(anchor="w")
        self.caption(pane, display_value(row, "location") + " · 0-based, end excluded")
        if row.get("target_sequence"):
            self.ttk.Label(pane, text=row["target_sequence"], font=(MONO, 18)).pack(anchor="w", pady=(12, 4))
            self.caption(pane, "Target sequence includes the PAM, in guide orientation.")
        extra = ("presence", "sample_copies") if "sample_copies" in row else ()
        for column in VIEW_COLUMNS[source] + extra:
            if column not in {"guide_id", "location", "target_sequence"}:
                self.summary_row(pane, COLUMN_SPECS[column][0], display_value(row, column))
        if "evidence_links" in row:
            try:
                links = json.loads(row["evidence_links"] or "[]")
            except ValueError:
                links = []
            self.summary_row(pane, "Evidence status", f"{len(links)} linked records" if links else "No linked evidence; untested in this bundle")
        self.caption(pane, "All fields · original export names and values")
        formatted = []
        for column, value in row.items():
            if value.startswith(("[", "{")):
                try:
                    value = json.dumps(json.loads(value), indent=2)
                except ValueError:
                    pass
            formatted.append(f"{column}: {value}")
        def close():
            popup.destroy()
            if previous_focus is not None and previous_focus.winfo_exists():
                previous_focus.focus_set()
        self.ttk.Button(pane, text="Close", command=close).pack(side="bottom", anchor="e", pady=(12, 0))
        self.detail_text(pane, "\n\n".join(formatted), height=10)
        popup.bind("<Escape>", lambda _: close())
        popup.protocol("WM_DELETE_WINDOW", close)
        popup.focus_set()

    # ------------------------------------------------------------------ shared

    def open_path(self, path):
        subprocess.run(["/usr/bin/open", str(path)], check=True)

    def close(self):
        if self.busy or self.export_busy or self.preparing:
            self.messages.showinfo("Operation in progress", "Keep this window open until the current operation finishes.")
        else:
            self.root.destroy()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, required=True)
    parser.add_argument("--run", choices=["small", "public"])
    parser.add_argument("--output", type=Path)
    parser.add_argument("--cache", type=Path)
    args = parser.parse_args(argv)
    os.umask(0o077)
    if args.run:
        if args.output is None or (args.run == "public" and args.cache is None):
            parser.error("--run requires --output; public also requires --cache")
        try:
            output, state = lab.run_job(args.run, args.output, args.bundle / "backend/cas-offinder", cache=args.cache)
            print(json.dumps({"status": state["status"], "output": str(output),
                              "example_verification": state.get("example_verification")}))
            if state.get("example_verification", {}).get("status") == "mismatched":
                return 2
            return 0 if state["status"] == "complete" else 3
        except Exception as exc:
            print("Action needed: " + str(exc), file=sys.stderr)
            return 2
    configure_tk()
    import tkinter as tk
    root = tk.Tk()
    apply_style(root)
    Launcher(root, args.bundle.resolve())
    root.mainloop()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
