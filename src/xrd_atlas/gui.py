from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

from .exporters import export_xrd_atlas_workbook
from .models import XrdAtlasExportPayload, XrdAtlasSettings
from .service import XrdAtlasService


@dataclass(frozen=True)
class SimpleGuiExportResult:
    output_path: Path
    total_peaks: int
    phase_rows: list[tuple[str, str, str, int, str]]


def _parse_float(value: str | float | int, field_name: str) -> float:
    try:
        number = float(str(value).strip())
    except ValueError as exc:
        raise ValueError(f"{field_name} must be a number.") from exc
    return number


def build_gui_settings(
    energy_keV: str | float,
    two_theta_min: str | float = 0.0,
    two_theta_max: str | float = 180.0,
) -> XrdAtlasSettings:
    energy = _parse_float(energy_keV, "X-ray energy keV")
    min_deg = _parse_float(two_theta_min, "2theta min")
    max_deg = _parse_float(two_theta_max, "2theta max")
    if energy <= 0:
        raise ValueError("X-ray energy keV must be greater than 0.")
    if min_deg < 0 or max_deg > 180 or min_deg >= max_deg:
        raise ValueError("2theta range must satisfy 0 <= min < max <= 180.")
    return XrdAtlasSettings(
        input_mode="energy",
        energy_keV=energy,
        two_theta_min_deg=min_deg,
        two_theta_max_deg=max_deg,
    )


def normalize_xlsx_output_path(output_path: str | Path) -> Path:
    path = Path(output_path).expanduser().resolve()
    if path.suffix.lower() != ".xlsx":
        path = path.with_suffix(".xlsx")
    return path


def run_simple_gui_export(
    cif_paths: Sequence[str | Path],
    output_path: str | Path,
    *,
    energy_keV: str | float,
    two_theta_min: str | float = 0.0,
    two_theta_max: str | float = 180.0,
) -> SimpleGuiExportResult:
    resolved_cifs = [Path(path).expanduser().resolve() for path in cif_paths]
    if not resolved_cifs:
        raise ValueError("Select at least one CIF file.")
    missing = [str(path) for path in resolved_cifs if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing CIF file(s): " + "; ".join(missing))

    settings = build_gui_settings(energy_keV, two_theta_min, two_theta_max)
    service = XrdAtlasService()
    phases = service.load_phases(resolved_cifs)
    service.simulate_phases(phases, settings)

    output = normalize_xlsx_output_path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    export_xrd_atlas_workbook(XrdAtlasExportPayload(phases, settings), output)

    phase_rows: list[tuple[str, str, str, int, str]] = []
    for phase in phases:
        peak_count = 0 if phase.result is None else len(phase.result.peaks)
        phase_rows.append(
            (
                phase.cif_path.name,
                phase.display_formula,
                phase.display_space_group,
                peak_count,
                phase.error or " | ".join(phase.warning_messages),
            )
        )
    return SimpleGuiExportResult(
        output_path=output,
        total_peaks=sum(row[3] for row in phase_rows),
        phase_rows=phase_rows,
    )


def _launch_tk_app() -> None:
    import threading
    import tkinter as tk
    from tkinter import filedialog, messagebox, ttk

    root = tk.Tk()
    root.title("XRD Atlas - CIF to Excel")
    root.geometry("860x560")
    root.minsize(760, 500)

    selected_paths: list[Path] = []
    energy_var = tk.StringVar(value="20.0")
    min_var = tk.StringVar(value="0")
    max_var = tk.StringVar(value="180")
    output_var = tk.StringVar(value=str(Path.home() / "Desktop" / "xrd_peak_reference.xlsx"))
    status_var = tk.StringVar(value="Add CIF files, set X-ray energy, then export.")

    root.columnconfigure(0, weight=1)
    root.rowconfigure(1, weight=1)

    header = ttk.Frame(root, padding=(14, 12, 14, 4))
    header.grid(row=0, column=0, sticky="ew")
    header.columnconfigure(0, weight=1)
    ttk.Label(header, text="XRD Atlas", font=("Segoe UI", 16, "bold")).grid(row=0, column=0, sticky="w")
    ttk.Label(header, text="Batch export theoretical powder XRD peak references from CIF files.").grid(
        row=1,
        column=0,
        sticky="w",
    )

    main = ttk.Frame(root, padding=14)
    main.grid(row=1, column=0, sticky="nsew")
    main.columnconfigure(0, weight=1)
    main.columnconfigure(1, weight=1)
    main.rowconfigure(1, weight=1)

    files_panel = ttk.LabelFrame(main, text="1. CIF input", padding=10)
    files_panel.grid(row=0, column=0, rowspan=2, sticky="nsew", padx=(0, 8))
    files_panel.columnconfigure(0, weight=1)
    files_panel.rowconfigure(1, weight=1)

    file_buttons = ttk.Frame(files_panel)
    file_buttons.grid(row=0, column=0, sticky="ew", pady=(0, 8))
    file_buttons.columnconfigure(2, weight=1)

    listbox = tk.Listbox(files_panel, height=14, selectmode=tk.EXTENDED)
    listbox.grid(row=1, column=0, sticky="nsew")
    scrollbar = ttk.Scrollbar(files_panel, orient="vertical", command=listbox.yview)
    scrollbar.grid(row=1, column=1, sticky="ns")
    listbox.configure(yscrollcommand=scrollbar.set)

    def refresh_list() -> None:
        listbox.delete(0, tk.END)
        for path in selected_paths:
            listbox.insert(tk.END, str(path))
        status_var.set(f"{len(selected_paths)} CIF file(s) selected.")

    def add_files() -> None:
        paths = filedialog.askopenfilenames(title="Select CIF files", filetypes=[("CIF files", "*.cif")])
        known = {path.resolve() for path in selected_paths}
        for item in paths:
            path = Path(item).resolve()
            if path not in known:
                selected_paths.append(path)
                known.add(path)
        refresh_list()

    def add_folder() -> None:
        folder = filedialog.askdirectory(title="Select folder containing CIF files")
        if not folder:
            return
        known = {path.resolve() for path in selected_paths}
        for path in sorted(Path(folder).rglob("*.cif"), key=lambda value: str(value).lower()):
            resolved = path.resolve()
            if resolved not in known:
                selected_paths.append(resolved)
                known.add(resolved)
        refresh_list()

    def remove_selected() -> None:
        selected = set(listbox.curselection())
        if not selected:
            return
        selected_paths[:] = [path for index, path in enumerate(selected_paths) if index not in selected]
        refresh_list()

    ttk.Button(file_buttons, text="Add files", command=add_files).grid(row=0, column=0, padx=(0, 6))
    ttk.Button(file_buttons, text="Add folder", command=add_folder).grid(row=0, column=1, padx=(0, 6))
    ttk.Button(file_buttons, text="Remove", command=remove_selected).grid(row=0, column=2, sticky="w")

    settings_panel = ttk.LabelFrame(main, text="2. X-ray and export settings", padding=10)
    settings_panel.grid(row=0, column=1, sticky="nsew", padx=(8, 0))
    settings_panel.columnconfigure(1, weight=1)

    ttk.Label(settings_panel, text="X-ray energy (keV)").grid(row=0, column=0, sticky="w", pady=4)
    ttk.Entry(settings_panel, textvariable=energy_var).grid(row=0, column=1, sticky="ew", pady=4)
    ttk.Label(settings_panel, text="2theta min (deg)").grid(row=1, column=0, sticky="w", pady=4)
    ttk.Entry(settings_panel, textvariable=min_var).grid(row=1, column=1, sticky="ew", pady=4)
    ttk.Label(settings_panel, text="2theta max (deg)").grid(row=2, column=0, sticky="w", pady=4)
    ttk.Entry(settings_panel, textvariable=max_var).grid(row=2, column=1, sticky="ew", pady=4)
    ttk.Label(settings_panel, text="Output Excel").grid(row=3, column=0, sticky="w", pady=4)

    output_frame = ttk.Frame(settings_panel)
    output_frame.grid(row=3, column=1, sticky="ew", pady=4)
    output_frame.columnconfigure(0, weight=1)
    ttk.Entry(output_frame, textvariable=output_var).grid(row=0, column=0, sticky="ew", padx=(0, 6))

    def choose_output() -> None:
        path = filedialog.asksaveasfilename(
            title="Save Excel workbook",
            defaultextension=".xlsx",
            filetypes=[("Excel workbook", "*.xlsx")],
            initialfile=Path(output_var.get()).name or "xrd_peak_reference.xlsx",
        )
        if path:
            output_var.set(path)

    ttk.Button(output_frame, text="Browse", command=choose_output).grid(row=0, column=1)

    preview_panel = ttk.LabelFrame(main, text="3. Export summary", padding=10)
    preview_panel.grid(row=1, column=1, sticky="nsew", padx=(8, 0), pady=(10, 0))
    preview_panel.columnconfigure(0, weight=1)
    preview_panel.rowconfigure(0, weight=1)

    columns = ("cif", "formula", "space_group", "peaks", "warning")
    tree = ttk.Treeview(preview_panel, columns=columns, show="headings", height=9)
    for column, label, width in (
        ("cif", "CIF", 170),
        ("formula", "Formula", 90),
        ("space_group", "Space group", 90),
        ("peaks", "Peaks", 60),
        ("warning", "Warnings / errors", 220),
    ):
        tree.heading(column, text=label)
        tree.column(column, width=width, anchor="w")
    tree.grid(row=0, column=0, sticky="nsew")
    tree_scroll = ttk.Scrollbar(preview_panel, orient="vertical", command=tree.yview)
    tree_scroll.grid(row=0, column=1, sticky="ns")
    tree.configure(yscrollcommand=tree_scroll.set)

    footer = ttk.Frame(root, padding=(14, 4, 14, 12))
    footer.grid(row=2, column=0, sticky="ew")
    footer.columnconfigure(0, weight=1)
    ttk.Label(footer, textvariable=status_var).grid(row=0, column=0, sticky="w")

    export_button = ttk.Button(footer, text="Export Excel")
    export_button.grid(row=0, column=1, sticky="e")

    def set_busy(is_busy: bool) -> None:
        export_button.configure(state=tk.DISABLED if is_busy else tk.NORMAL)

    def export_now() -> None:
        set_busy(True)
        status_var.set("Running XRD simulation and exporting workbook...")
        tree.delete(*tree.get_children())

        def worker() -> None:
            try:
                result = run_simple_gui_export(
                    selected_paths,
                    output_var.get(),
                    energy_keV=energy_var.get(),
                    two_theta_min=min_var.get(),
                    two_theta_max=max_var.get(),
                )
            except Exception as exc:
                root.after(0, lambda: (set_busy(False), status_var.set("Export failed."), messagebox.showerror("Export failed", str(exc))))
                return

            def finish() -> None:
                for row in result.phase_rows:
                    tree.insert("", tk.END, values=row)
                set_busy(False)
                status_var.set(f"Exported {result.total_peaks} peaks to {result.output_path}")
                messagebox.showinfo("Export complete", f"Exported {result.total_peaks} peaks.\n\n{result.output_path}")

            root.after(0, finish)

        threading.Thread(target=worker, daemon=True).start()

    export_button.configure(command=export_now)
    root.mainloop()


def main() -> int:
    _launch_tk_app()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
