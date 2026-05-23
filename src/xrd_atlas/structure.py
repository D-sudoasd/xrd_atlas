from __future__ import annotations

from pathlib import Path
import warnings

import gemmi
import numpy as np
import spglib
from pymatgen.core import Structure as PymatgenStructure
from pymatgen.io.cif import CifParser

from .models import CrystalModel, StructureValidationReport
from .utils import file_sha256


OCCUPANCY_FALLBACK_WARNING = (
    "CIF occupancy 严格解析失败，已使用 occupancy_tolerance=4.0 兼容导入并由 pymatgen "
    "归一化占位；相对强度仅作理论近似参考。"
)


def _clean_cif_scalar(value: str | None) -> str | None:
    if value is None:
        return None
    text = str(value).strip().strip("'\"")
    return text or None


def _read_formula(block: gemmi.cif.Block, structure: PymatgenStructure) -> str:
    formula = _clean_cif_scalar(block.find_value("_chemical_formula_sum")) or _clean_cif_scalar(
        block.find_value("_chemical_name_systematic")
    )
    if formula:
        return formula
    return structure.composition.reduced_formula


def _has_known_cif_value(value: str | None) -> bool:
    return _clean_cif_scalar(value) not in {None, "?", "."}


def _block_has_cell_parameters(block: gemmi.cif.Block) -> bool:
    return all(
        _has_known_cif_value(block.find_value(tag))
        for tag in (
            "_cell_length_a",
            "_cell_length_b",
            "_cell_length_c",
            "_cell_angle_alpha",
            "_cell_angle_beta",
            "_cell_angle_gamma",
        )
    )


def _valid_atom_site_count(block: gemmi.cif.Block) -> int:
    labels = block.find_loop("_atom_site_label")
    types = block.find_loop("_atom_site_type_symbol")
    x_values = block.find_loop("_atom_site_fract_x")
    y_values = block.find_loop("_atom_site_fract_y")
    z_values = block.find_loop("_atom_site_fract_z")
    row_count = min(len(x_values), len(y_values), len(z_values))
    valid_count = 0
    for index in range(row_count):
        has_site_id = index < len(labels) and _has_known_cif_value(labels[index])
        has_type = index < len(types) and _has_known_cif_value(types[index])
        has_position = all(
            _has_known_cif_value(values[index])
            for values in (x_values, y_values, z_values)
        )
        if has_position and (has_site_id or has_type):
            valid_count += 1
    return valid_count


def _is_structure_block(block: gemmi.cif.Block) -> bool:
    return _block_has_cell_parameters(block) and _valid_atom_site_count(block) > 0


def _select_structure_block(document: gemmi.cif.Document) -> gemmi.cif.Block:
    candidates = [block for block in document if _is_structure_block(block)]
    if not candidates:
        raise ValueError("CIF file contains no data block with valid cell parameters and atom coordinates.")

    for preferred_name in ("standardized_unitcell", "published_cell"):
        for block in candidates:
            if preferred_name in block.name.lower():
                return block
    return candidates[0]


def _spglib_dataset(structure: PymatgenStructure) -> tuple[int | None, str | None, list[str]]:
    warnings: list[str] = []
    lattice = np.asarray(structure.lattice.matrix, dtype=float)
    positions = np.asarray(structure.frac_coords, dtype=float)
    numbers: list[int] = []
    for site in structure:
        # spglib needs a single integer species per site. For disordered sites,
        # use the most occupied species only for symmetry detection; XRD itself
        # still uses pymatgen's disordered structure.
        species = site.species
        dominant = max(species.items(), key=lambda item: float(item[1]))[0]
        numbers.append(int(dominant.Z))
    dataset = spglib.get_symmetry_dataset((lattice, positions, np.asarray(numbers, dtype=int)), symprec=1e-3)
    if dataset is None:
        warnings.append("spglib 未能识别该结构的空间群；XRD 仍按 CIF 中的结构计算。")
        return None, None, warnings
    return int(dataset.number), str(dataset.international), warnings


def _has_partial_occupancy(structure: PymatgenStructure) -> bool:
    return any(any(abs(float(occ) - 1.0) > 1e-8 for occ in site.species.values()) for site in structure)


def _is_occupancy_parse_failure(exc: Exception, warning_messages: list[str]) -> bool:
    text = "\n".join([str(exc), *warning_messages]).lower()
    return "occupancy" in text or "occupancies" in text


def _load_pymatgen_structure(path: Path) -> tuple[PymatgenStructure, list[str]]:
    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            return PymatgenStructure.from_file(str(path)), []
    except Exception as exc:
        strict_warnings = [str(item.message) for item in caught]
        if not _is_occupancy_parse_failure(exc, strict_warnings):
            raise

    with warnings.catch_warnings(record=True):
        warnings.simplefilter("ignore")
        parser = CifParser(str(path), occupancy_tolerance=4.0)
        structures = parser.parse_structures(primitive=False)
    if not structures:
        raise ValueError(f"CIF 文件未解析出有效结构：{path}")
    return structures[0], [OCCUPANCY_FALLBACK_WARNING]


def load_crystal_model(cif_path: str | Path) -> CrystalModel:
    path = Path(cif_path).expanduser().resolve()
    if not path.exists():
        raise FileNotFoundError(f"CIF 文件不存在：{path}")

    document = gemmi.cif.read_file(str(path))
    block = _select_structure_block(document)
    structure, parser_warnings = _load_pymatgen_structure(path)
    detected_number, detected_symbol, warnings = _spglib_dataset(structure)
    cell = structure.lattice
    report = StructureValidationReport(
        warnings=[*parser_warnings, *warnings],
        space_group_detected=detected_symbol,
        space_group_from_cif=_clean_cif_scalar(block.find_value("_symmetry_space_group_name_H-M"))
        or _clean_cif_scalar(block.find_value("_space_group_name_H-M_alt")),
        occupancy_summary="存在部分占位，按平均结构计算。" if _has_partial_occupancy(structure) else "所有位点按全占位处理。",
    )
    return CrystalModel(
        cif_path=path,
        cif_hash=file_sha256(path),
        formula=_read_formula(block, structure),
        space_group_number=detected_number,
        space_group_symbol=report.space_group_from_cif,
        detected_space_group_number=detected_number,
        detected_space_group_symbol=detected_symbol,
        cell_parameters=(float(cell.a), float(cell.b), float(cell.c), float(cell.alpha), float(cell.beta), float(cell.gamma)),
        pymatgen_structure=structure,
        validation_report=report,
        has_partial_occupancy=_has_partial_occupancy(structure),
    )
