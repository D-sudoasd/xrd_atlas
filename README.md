# XRD Atlas

XRD Atlas is a lightweight GUI and CLI tool for converting CIF crystal
structures into theoretical powder XRD peak reference tables.

It is designed for materials researchers who need a practical way to batch
export phase, hkl, d-spacing, 2theta, q, g, relative intensity, and warnings
into Excel for follow-up work in Excel, Origin, Python, or lab notebooks.

## Features

- Load one CIF file, many CIF files, or a folder of CIF files.
- Export one combined Excel workbook with:
  - `Summary`
  - `Combined Peaks`
  - one sheet per phase
- Use the desktop GUI for non-programming workflows.
- Use the CLI for reproducible batch processing.
- Choose common X-ray sources or enter custom energy / wavelength.
- Keep going when one CIF fails; errors and warnings are written to `Summary`.
- Handles common multi-block CIF files by selecting the structural block.

## Install

Python 3.11 or newer is required.

```powershell
cd C:\Users\AORUS\Desktop\xrd_atlas
py -3.11 -m pip install -e .[dev]
```

## GUI Quick Start

Start the desktop app:

```powershell
xrd-atlas-gui
```

Or run it from the project folder:

```powershell
py -3.11 -m xrd_atlas.gui
```

Basic workflow:

1. Click `Add files` or `Add folder`.
2. Enter X-ray energy in keV.
3. Set the 2theta range.
4. Choose the output `.xlsx` path.
5. Click `Export Excel`.

## CLI Examples

Export all CIF files in a folder:

```powershell
xrd-atlas "C:\path\to\cif_folder" -o result.xlsx
```

Export several CIF files:

```powershell
xrd-atlas phase1.cif phase2.cif phase3.cif -o result.xlsx
```

Use a custom X-ray energy:

```powershell
xrd-atlas "C:\path\to\cif_folder" -o result.xlsx --energy-keV 20
```

Use a custom wavelength:

```powershell
xrd-atlas "C:\path\to\cif_folder" -o result.xlsx --wavelength-A 1.5406
```

Limit the 2theta range:

```powershell
xrd-atlas "C:\path\to\cif_folder" -o result.xlsx --source "Cu Ka" --two-theta-min 20 --two-theta-max 100
```

Export CSV instead of Excel:

```powershell
xrd-atlas "C:\path\to\cif_folder" -o result.csv
```

## Output Columns

The peak tables include:

- `phase_name`
- `cif_name`
- `formula`
- `space_group`
- `hkl`
- `d_A`
- `two_theta_current_deg`
- `relative_intensity`
- `multiplicity`
- `family_label`
- `h`, `k`, `l`
- `g_1_over_A`
- `q_1_over_A`
- `theta_deg`
- `two_theta_cu_ka_deg`
- `warnings`

## Scientific Scope

XRD Atlas exports theoretical powder XRD peak references from CIF structures.

It is not:

- an experimental pattern fitting program
- a phase identification database
- a Rietveld refinement tool
- a replacement for instrument calibration

For phase and peak-position comparison, prioritize `phase_name`, `hkl`, `d_A`,
and `two_theta_current_deg`. Treat relative intensity as a theoretical
reference, not as a refined experimental quantity.

## Tests

```powershell
cd C:\Users\AORUS\Desktop\xrd_atlas
py -3.11 -m pytest -q
```

## License

MIT License. See [LICENSE](LICENSE).
