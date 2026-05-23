from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

from .constants import DEFAULT_XRD_SOURCE, XRD_SOURCE_PRESETS
from .exporters import export_peak_reference_csv, export_xrd_atlas_json, export_xrd_atlas_workbook
from .models import XrdAtlasExportPayload, XrdAtlasSettings
from .service import XrdAtlasService


def collect_cif_paths(inputs: Sequence[str | Path], *, recursive: bool = True) -> list[Path]:
    paths: list[Path] = []
    seen: set[Path] = set()
    for item in inputs:
        path = Path(item).expanduser().resolve()
        candidates: list[Path]
        if path.is_dir():
            iterator = path.rglob("*.cif") if recursive else path.glob("*.cif")
            candidates = sorted(iterator, key=lambda value: str(value).lower())
        else:
            candidates = [path]
        for candidate in candidates:
            resolved = candidate.expanduser().resolve()
            if resolved.suffix.lower() == ".cif" and resolved not in seen:
                paths.append(resolved)
                seen.add(resolved)
    return paths


def batch_export_peak_reference(
    inputs: Sequence[str | Path],
    output_path: str | Path,
    settings: XrdAtlasSettings | None = None,
    *,
    recursive: bool = True,
) -> list:
    cif_paths = collect_cif_paths(inputs, recursive=recursive)
    if not cif_paths:
        raise FileNotFoundError("未找到 CIF 文件。")

    service = XrdAtlasService()
    resolved_settings = settings or XrdAtlasSettings()
    phases = service.load_phases(cif_paths)
    service.simulate_phases(phases, resolved_settings)
    payload = XrdAtlasExportPayload(phases, resolved_settings)

    output = Path(output_path).expanduser().resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.suffix.lower() == ".xlsx":
        export_xrd_atlas_workbook(payload, output)
        export_xrd_atlas_json(payload, output.with_suffix(".json"))
    else:
        export_peak_reference_csv(payload, output)
    return phases


def _settings_from_args(args: argparse.Namespace) -> XrdAtlasSettings:
    input_mode = "source"
    if args.energy_keV is not None:
        input_mode = "energy"
    elif args.wavelength_A is not None:
        input_mode = "wavelength"

    return XrdAtlasSettings(
        input_mode=input_mode,  # type: ignore[arg-type]
        source_preset=args.source,
        wavelength_A=args.wavelength_A if args.wavelength_A is not None else XRD_SOURCE_PRESETS[DEFAULT_XRD_SOURCE],
        energy_keV=args.energy_keV if args.energy_keV is not None else XrdAtlasSettings().energy_keV,
        two_theta_min_deg=args.two_theta_min,
        two_theta_max_deg=args.two_theta_max,
    )


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Batch-export phase/hkl/d/2theta peak reference tables from CIF files.",
    )
    parser.add_argument("inputs", nargs="+", help="CIF files or directories containing CIF files.")
    parser.add_argument("-o", "--output", default="xrd_peak_reference.xlsx", help="Output .xlsx or .csv path.")
    parser.add_argument("--no-recursive", action="store_true", help="Do not recurse into input directories.")
    parser.add_argument("--source", default=DEFAULT_XRD_SOURCE, choices=sorted(XRD_SOURCE_PRESETS), help="X-ray source preset.")
    parser.add_argument("--energy-keV", type=float, default=None, help="Use a custom X-ray energy in keV.")
    parser.add_argument("--wavelength-A", type=float, default=None, help="Use a custom wavelength in Angstrom.")
    parser.add_argument("--two-theta-min", type=float, default=0.0, help="Minimum 2theta in degrees.")
    parser.add_argument("--two-theta-max", type=float, default=180.0, help="Maximum 2theta in degrees.")
    args = parser.parse_args(argv)

    settings = _settings_from_args(args)
    phases = batch_export_peak_reference(args.inputs, args.output, settings, recursive=not args.no_recursive)
    peak_count = sum(0 if phase.result is None else len(phase.result.peaks) for phase in phases)
    failed = [phase.cif_path.name for phase in phases if phase.error]
    print(f"Exported {peak_count} peaks from {len(phases)} CIF files to {Path(args.output).resolve()}")
    if failed:
        print("Failed CIF files: " + ", ".join(failed))
    return 0 if peak_count else 1


if __name__ == "__main__":
    raise SystemExit(main())
