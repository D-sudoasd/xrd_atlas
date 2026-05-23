from __future__ import annotations

import csv
import os
import subprocess
import sys
import textwrap
import tomllib
import xml.etree.ElementTree as ET
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pytest

from xrd_atlas.batch import batch_export_peak_reference
from xrd_atlas.exporters import combined_peak_rows, export_peak_reference_csv, export_xrd_atlas_workbook
from xrd_atlas.models import XrdAtlasExportPayload, XrdAtlasSettings
from xrd_atlas.service import XrdAtlasService


ROOT = Path(__file__).resolve().parents[1]
EXAMPLES_CIF_DIR = ROOT / "examples" / "cif"
TI_BETA_CIF = EXAMPLES_CIF_DIR / "ti_beta_bcc_im3m.cif"
TI_NB_HCP_CIF = EXAMPLES_CIF_DIR / "ti_nb_hcp_p63mmc.cif"

NI_OCCUPANCY_CIF = """
#======================================================================
# CRYSTAL DATA
#----------------------------------------------------------------------
data_VESTA_phase_1

_chemical_name_common                  'Ni                                   '
_cell_length_a                         3.598240
_cell_length_b                         3.598240
_cell_length_c                         3.598240
_cell_angle_alpha                      90.000000
_cell_angle_beta                       90.000000
_cell_angle_gamma                      90.000000
_cell_volume                           46.587601
_space_group_name_H-M_alt              'P 2 3'
_space_group_IT_number                 195

loop_
_space_group_symop_operation_xyz
   'x, y, z'
   '-x, -y, z'
   '-x, y, -z'
   'x, -y, -z'
   'z, x, y'
   'z, -x, -y'
   '-z, -x, y'
   '-z, x, -y'
   'y, z, x'
   '-y, z, -x'
   'y, -z, -x'
   '-y, -z, x'

loop_
   _atom_site_label
   _atom_site_occupancy
   _atom_site_fract_x
   _atom_site_fract_y
   _atom_site_fract_z
   _atom_site_adp_type
   _atom_site_U_iso_or_equiv
   _atom_site_type_symbol
   Ni1        1.0     0.000000     0.000000     -0.000000    Uiso  ? Ni
   Ni2        1.0     0.500000     0.500000     0.000000    Uiso  ? Ni
   Ni3        1.0     -0.000000     0.500000     0.500000    Uiso  ? Ni
   Ni4        1.0     0.500000     -0.000000     0.500000    Uiso  ? Ni
"""

MULTI_BLOCK_STANDARDIZED_CIF = """
data_sm_global
_audit_creation_method 'metadata only'

data_probe-standardized_unitcell
_cell_length_a 3.0
_cell_length_b 3.0
_cell_length_c 3.0
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_symmetry_space_group_name_H-M 'P 1'
_symmetry_Int_Tables_number 1
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
Si1 Si 0 0 0 1

data_probe-published_cell
_cell_length_a 4.0
_cell_length_b 4.0
_cell_length_c 4.0
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_symmetry_space_group_name_H-M 'P 2 3'
_symmetry_Int_Tables_number 195
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
Si1 Si 0 0 0 1

data_probe-niggli_reduced_cell
_cell_length_a 2.0
_cell_length_b 2.0
_cell_length_c 2.0
_cell_angle_alpha 60
_cell_angle_beta 60
_cell_angle_gamma 60
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
? ? ? ? ? ?
"""


def _write_ni_occupancy_cif(path: Path) -> Path:
    path.write_text(textwrap.dedent(NI_OCCUPANCY_CIF).strip() + "\n", encoding="utf-8")
    return path


def _write_multi_block_standardized_cif(path: Path) -> Path:
    path.write_text(textwrap.dedent(MULTI_BLOCK_STANDARDIZED_CIF).strip() + "\n", encoding="utf-8")
    return path


def _nb_hea_cif_paths() -> list[Path]:
    root = Path.home() / "Desktop" / "Nb_HEA_peak_separation"
    names = ["AlNi.cif", "Cr2Nb.cif", "FeCr.cif", "Ni.cif", "Ni3Al.cif"]
    found = {path.name: path for path in root.rglob("*.cif")} if root.exists() else {}
    paths = [found.get(name) for name in names]
    if any(path is None for path in paths):
        pytest.skip("Nb_HEA_peak_separation CIF fixtures are not available on this machine.")
    return [path for path in paths if path is not None]


def _zr_hydride_cif_paths() -> list[Path]:
    root = Path.home() / "Desktop" / "ZrNb_SXRD_deformation" / "Cif"
    names = [
        "\N{GREEK SMALL LETTER ALPHA}-Zr.cif",
        "\N{GREEK SMALL LETTER BETA}-Zr.cif",
        "\N{GREEK SMALL LETTER GAMMA}-ZrH.cif",
        "\N{GREEK SMALL LETTER DELTA}-ZrH1.66.cif",
    ]
    paths = [root / name for name in names]
    if any(not path.exists() for path in paths):
        pytest.skip("Zr/ZrH CIF fixtures are not available on this machine.")
    return paths


def _subprocess_env() -> dict[str, str]:
    env = os.environ.copy()
    pythonpath = str(ROOT / "src")
    if env.get("PYTHONPATH"):
        pythonpath = pythonpath + os.pathsep + env["PYTHONPATH"]
    env["PYTHONPATH"] = pythonpath
    return env


def _worksheet_rows(workbook_path: Path, sheet_index: int) -> list[list[str]]:
    namespace = {"main": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    with ZipFile(workbook_path) as archive:
        xml = archive.read(f"xl/worksheets/sheet{sheet_index}.xml")
    root = ET.fromstring(xml)
    rows: list[list[str]] = []
    for row in root.findall(".//main:row", namespace):
        values: list[str] = []
        for cell in row.findall("main:c", namespace):
            inline = cell.find("main:is/main:t", namespace)
            numeric = cell.find("main:v", namespace)
            values.append("" if inline is None and numeric is None else ((inline.text or "") if inline is not None else numeric.text or ""))
        rows.append(values)
    return rows


def test_project_is_cli_only_without_gui_dependencies() -> None:
    project = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    dependencies = "\n".join(project["project"]["dependencies"])
    scripts = project["project"]["scripts"]
    requirements = (ROOT / "requirements.txt").read_text(encoding="utf-8")

    assert "PySide6" not in dependencies
    assert "matplotlib" not in dependencies
    assert "PySide6" not in requirements
    assert "matplotlib" not in requirements
    assert scripts["xrd-atlas"] == "xrd_atlas.batch:main"
    assert scripts["xrd-atlas-peaks"] == "xrd_atlas.batch:main"
    assert scripts["xrd-atlas-gui"] == "xrd_atlas.gui:main"


def test_xrd_atlas_single_cif_peak_table_and_energy_shift() -> None:
    service = XrdAtlasService()
    phase = service.load_phase(TI_BETA_CIF)
    settings = XrdAtlasSettings(input_mode="energy", energy_keV=8.0478, wavelength_A=1.5406)
    service.simulate_phase(phase, settings)
    assert phase.result is not None
    rows = combined_peak_rows([phase])
    assert rows
    first = rows[0]
    for key in ("d_A", "two_theta_current_deg", "q_1_over_A", "g_1_over_A", "multiplicity"):
        assert key in first

    high_energy_phase = service.load_phase(TI_BETA_CIF)
    service.simulate_phase(
        high_energy_phase,
        XrdAtlasSettings(input_mode="energy", energy_keV=20.0, wavelength_A=1.5406),
    )
    assert high_energy_phase.result is not None
    assert high_energy_phase.result.peaks[0].two_theta_deg < phase.result.peaks[0].two_theta_deg
    assert np.isclose(high_energy_phase.result.peaks[0].d_spacing_A, phase.result.peaks[0].d_spacing_A)


def test_xrd_atlas_loads_occupancy_conflict_cif_with_warning(tmp_path: Path) -> None:
    cif_path = _write_ni_occupancy_cif(tmp_path / "Ni.cif")
    service = XrdAtlasService()

    phase = service.load_phase(cif_path)
    service.simulate_phase(phase, XrdAtlasSettings())

    assert phase.error is None
    assert phase.result is not None
    assert len(phase.result.peaks) == 8
    assert any("occupancy" in warning.lower() for warning in phase.warning_messages)


def test_xrd_atlas_loads_standardized_unitcell_from_multi_block_cif(tmp_path: Path) -> None:
    cif_path = _write_multi_block_standardized_cif(tmp_path / "multi_block.cif")
    service = XrdAtlasService()

    phase = service.load_phase(cif_path)

    assert phase.error is None
    assert phase.crystal is not None
    assert phase.display_space_group == "P 1"
    assert phase.crystal.cell_parameters[:3] == (3.0, 3.0, 3.0)


def test_xrd_atlas_batch_loads_real_cifs_and_keeps_occupancy_warning() -> None:
    service = XrdAtlasService()
    phases = service.load_phases(_nb_hea_cif_paths())
    service.simulate_phases(phases, XrdAtlasSettings())

    by_name = {phase.cif_path.name: phase for phase in phases}
    assert list(by_name) == ["AlNi.cif", "Cr2Nb.cif", "FeCr.cif", "Ni.cif", "Ni3Al.cif"]
    expected_peak_counts = {
        "AlNi.cif": 12,
        "Cr2Nb.cif": 22,
        "FeCr.cif": 7,
        "Ni.cif": 8,
        "Ni3Al.cif": 86,
    }
    for name, expected_count in expected_peak_counts.items():
        phase = by_name[name]
        assert phase.error is None
        assert phase.result is not None
        assert len(phase.result.peaks) == expected_count
    assert any("occupancy" in warning.lower() for warning in by_name["Ni.cif"].warning_messages)


def test_xrd_atlas_exports_simple_phase_hkl_d_two_theta_csv(tmp_path: Path) -> None:
    service = XrdAtlasService()
    phases = service.load_phases(_nb_hea_cif_paths())
    settings = XrdAtlasSettings()
    service.simulate_phases(phases, settings)
    output = tmp_path / "phase_hkl_d_2theta.csv"

    export_peak_reference_csv(XrdAtlasExportPayload(phases, settings), output)

    with output.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert rows
    assert {"phase_name", "cif_name", "formula", "space_group", "hkl", "d_A", "two_theta_deg"}.issubset(rows[0])
    assert {row["cif_name"] for row in rows} == {"AlNi.cif", "Cr2Nb.cif", "FeCr.cif", "Ni.cif", "Ni3Al.cif"}
    assert any(row["cif_name"] == "Ni.cif" and row["hkl"] == "(1 1 1)" for row in rows)
    assert all(float(row["d_A"]) > 0 and float(row["two_theta_deg"]) > 0 for row in rows)


def test_xrd_atlas_batch_export_defaults_to_excel_from_many_cifs(tmp_path: Path) -> None:
    output = tmp_path / "many_cif_reference.xlsx"

    phases = batch_export_peak_reference(_nb_hea_cif_paths(), output)

    assert len(phases) == 5
    assert output.exists()
    headers = _worksheet_rows(output, 2)[0]
    expected_front = [
        "phase_name",
        "cif_name",
        "formula",
        "space_group",
        "hkl",
        "d_A",
        "two_theta_current_deg",
        "relative_intensity",
        "multiplicity",
        "warnings",
    ]
    assert headers[: len(expected_front)] == expected_front
    data_rows = _worksheet_rows(output, 2)[1:]
    assert len(data_rows) == sum(len(phase.result.peaks) for phase in phases if phase.result is not None)
    assert all(float(row[5]) > 0 and float(row[6]) > 0 for row in data_rows)


def test_xrd_atlas_batch_export_zr_hydride_multi_block_cifs(tmp_path: Path) -> None:
    output = tmp_path / "zr_hydride_reference.xlsx"

    phases = batch_export_peak_reference(_zr_hydride_cif_paths(), output)

    assert len(phases) == 4
    assert output.exists()
    assert all(phase.error is None for phase in phases)
    assert all(phase.result is not None and len(phase.result.peaks) > 0 for phase in phases)
    headers = _worksheet_rows(output, 2)[0]
    for header in ("hkl", "d_A", "two_theta_current_deg", "relative_intensity"):
        assert header in headers
    data_rows = _worksheet_rows(output, 2)[1:]
    assert len(data_rows) == sum(len(phase.result.peaks) for phase in phases if phase.result is not None)
    assert {row[1] for row in data_rows} == {path.name for path in _zr_hydride_cif_paths()}
    assert all(float(row[5]) > 0 and float(row[6]) > 0 for row in data_rows)


def test_xrd_atlas_batch_export_can_still_write_csv(tmp_path: Path) -> None:
    output = tmp_path / "many_cif_reference.csv"

    phases = batch_export_peak_reference(_nb_hea_cif_paths(), output)

    with output.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert len(phases) == 5
    assert output.exists()
    assert len(rows) == sum(len(phase.result.peaks) for phase in phases if phase.result is not None)
    assert {row["phase_name"] for row in rows} >= {"AlNi", "Cr2Nb", "FeCr", "Ni", "Ni3Al"}


def test_xrd_atlas_batch_load_isolates_unrecoverable_bad_cif(tmp_path: Path) -> None:
    bad_cif = tmp_path / "bad.cif"
    bad_cif.write_text("data_bad\n_cell_length_a 3\n", encoding="utf-8")
    service = XrdAtlasService()

    phases = service.load_phases([TI_BETA_CIF, bad_cif])
    service.simulate_phases(phases, XrdAtlasSettings())
    output = tmp_path / "with_bad_cif.xlsx"
    export_xrd_atlas_workbook(XrdAtlasExportPayload(phases, XrdAtlasSettings()), output)

    assert len(phases) == 2
    assert phases[0].result is not None
    assert phases[0].error is None
    assert phases[1].result is None
    assert phases[1].error
    assert not phases[1].enabled
    summary_text = "\n".join("|".join(row) for row in _worksheet_rows(output, 1))
    assert "bad.cif" in summary_text
    assert "error" in summary_text.lower() or "invalid" in summary_text.lower()


def test_xrd_atlas_multi_phase_workbook_export(tmp_path: Path) -> None:
    service = XrdAtlasService()
    phases = [service.load_phase(TI_BETA_CIF), service.load_phase(TI_NB_HCP_CIF)]
    settings = XrdAtlasSettings()
    service.simulate_phases(phases, settings)
    output = tmp_path / "xrd_atlas_export.xlsx"
    export_xrd_atlas_workbook(XrdAtlasExportPayload(phases, settings), output)
    with ZipFile(output) as archive:
        workbook = archive.read("xl/workbook.xml").decode("utf-8")
        combined = archive.read("xl/worksheets/sheet2.xml").decode("utf-8")
    assert "Summary" in workbook
    assert "Combined Peaks" in workbook
    assert "ti_beta_bcc_im3m" in workbook
    assert "ti_nb_hcp_p63mmc" in workbook
    assert combined.count("<row ") > 2


def test_batch_module_help_and_package_main_cli(tmp_path: Path) -> None:
    help_result = subprocess.run(
        [sys.executable, "-m", "xrd_atlas.batch", "--help"],
        cwd=ROOT,
        env=_subprocess_env(),
        text=True,
        capture_output=True,
        timeout=10,
        check=False,
    )
    assert help_result.returncode == 0
    assert "Batch-export" in help_result.stdout

    output = tmp_path / "package_main.xlsx"
    run_result = subprocess.run(
        [sys.executable, "-m", "xrd_atlas", str(TI_BETA_CIF), "-o", str(output)],
        cwd=ROOT,
        env=_subprocess_env(),
        text=True,
        capture_output=True,
        timeout=15,
        check=False,
    )
    assert run_result.returncode == 0, run_result.stderr
    assert output.exists()
    assert "Exported" in run_result.stdout


def test_batch_cli_default_output_is_excel(tmp_path: Path) -> None:
    result = subprocess.run(
        [sys.executable, "-m", "xrd_atlas.batch", str(TI_BETA_CIF)],
        cwd=tmp_path,
        env=_subprocess_env(),
        text=True,
        capture_output=True,
        timeout=15,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert (tmp_path / "xrd_peak_reference.xlsx").exists()
    assert not (tmp_path / "xrd_peak_reference.csv").exists()


def test_simple_gui_builds_energy_settings_from_user_inputs() -> None:
    from xrd_atlas.gui import build_gui_settings

    settings = build_gui_settings("20.0", "5", "120")

    assert settings.input_mode == "energy"
    assert settings.energy_keV == 20.0
    assert settings.two_theta_min_deg == 5.0
    assert settings.two_theta_max_deg == 120.0


def test_simple_gui_export_writes_xlsx_without_json_sidecar(tmp_path: Path) -> None:
    from xrd_atlas.gui import run_simple_gui_export

    output_without_suffix = tmp_path / "gui_reference"

    result = run_simple_gui_export([TI_BETA_CIF], output_without_suffix, energy_keV="20.0")

    assert result.output_path == output_without_suffix.with_suffix(".xlsx")
    assert result.output_path.exists()
    assert not result.output_path.with_suffix(".json").exists()
    assert result.total_peaks > 0
    assert result.phase_rows[0][0] == "ti_beta_bcc_im3m.cif"
    assert result.phase_rows[0][1] == "Ti"
    assert result.phase_rows[0][3] == result.total_peaks
    assert result.phase_rows[0][4] == ""
