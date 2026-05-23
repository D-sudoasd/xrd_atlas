from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Literal
from uuid import uuid4

import numpy as np


XrdAxisMode = Literal["two_theta", "d_spacing", "q", "g"]
XrayInputMode = Literal["source", "wavelength", "energy"]


@dataclass
class StructureValidationReport:
    warnings: list[str] = field(default_factory=list)
    space_group_detected: str | None = None
    space_group_from_cif: str | None = None
    occupancy_summary: str = ""


@dataclass
class CrystalModel:
    cif_path: Path
    cif_hash: str
    formula: str
    space_group_number: int | None
    space_group_symbol: str | None
    detected_space_group_number: int | None
    detected_space_group_symbol: str | None
    cell_parameters: tuple[float, float, float, float, float, float]
    pymatgen_structure: Any
    validation_report: StructureValidationReport
    has_partial_occupancy: bool = False


@dataclass(frozen=True)
class XRDRequest:
    cif_path: Path
    two_theta_min_deg: float
    two_theta_max_deg: float
    step_deg: float
    input_mode: XrayInputMode = "source"
    source_preset: str = "Cu Ka"
    wavelength_A: float | None = None
    energy_keV: float | None = None
    profile_model: str = "pseudo_voigt"
    fwhm_deg: float = 0.15
    show_hkl_labels: bool = True
    show_sticks: bool = True

    def cache_key(self) -> tuple[Any, ...]:
        return (
            str(self.cif_path.resolve()),
            round(self.two_theta_min_deg, 6),
            round(self.two_theta_max_deg, 6),
            round(self.step_deg, 6),
            self.input_mode,
            self.source_preset,
            None if self.wavelength_A is None else round(self.wavelength_A, 6),
            None if self.energy_keV is None else round(self.energy_keV, 6),
            self.profile_model,
            round(self.fwhm_deg, 6),
            bool(self.show_hkl_labels),
            bool(self.show_sticks),
        )


@dataclass(frozen=True)
class XRDPeakRecord:
    hkl: tuple[int, int, int]
    d_spacing_A: float
    theta_deg: float
    two_theta_deg: float
    q_invA: float
    g_invA: float
    intensity: float
    normalized_intensity: float
    multiplicity: int
    label: str
    family_label: str


@dataclass
class XRDResult:
    two_theta_grid: np.ndarray
    intensity_profile: np.ndarray
    stick_positions_deg: np.ndarray
    stick_intensities: np.ndarray
    peaks: list[XRDPeakRecord]
    metadata: dict[str, Any]
    warnings: list[str]
    summary_sections: dict[str, list[str]] = field(default_factory=dict)
    quality_report: dict[str, Any] | None = None


@dataclass
class XrdAtlasSettings:
    input_mode: XrayInputMode = "source"
    source_preset: str = "Cu Ka"
    wavelength_A: float = 1.5406
    energy_keV: float = 8.047786474944927
    two_theta_min_deg: float = 0.0
    two_theta_max_deg: float = 180.0
    step_deg: float = 0.02
    fwhm_deg: float = 0.15
    profile_model: str = "pseudo_voigt"
    show_labels: bool = False
    show_sticks: bool = True
    axis_mode: XrdAxisMode = "two_theta"


@dataclass
class XrdPhase:
    cif_path: Path
    phase_name: str
    phase_id: str = field(default_factory=lambda: uuid4().hex)
    enabled: bool = True
    crystal: CrystalModel | None = None
    result: XRDResult | None = None
    error: str | None = None

    @property
    def display_space_group(self) -> str:
        if self.crystal is None:
            return "-"
        return self.crystal.space_group_symbol or self.crystal.detected_space_group_symbol or "-"

    @property
    def display_formula(self) -> str:
        return "-" if self.crystal is None else self.crystal.formula

    @property
    def warning_messages(self) -> list[str]:
        warnings: list[str] = []
        if self.crystal is not None:
            warnings.extend(self.crystal.validation_report.warnings)
        if self.result is not None:
            warnings.extend(self.result.warnings)
        deduplicated: list[str] = []
        seen: set[str] = set()
        for warning in warnings:
            if warning not in seen:
                deduplicated.append(warning)
                seen.add(warning)
        return deduplicated


@dataclass
class ExperimentalPattern:
    path: Path
    label: str
    x_values: np.ndarray
    intensity: np.ndarray
    axis_mode: XrdAxisMode


@dataclass(frozen=True)
class XrdAtlasPeakRow:
    phase_name: str
    cif_name: str
    h: int
    k: int
    l: int
    family_label: str
    d_A: float
    g_1_over_A: float
    q_1_over_A: float
    theta_deg: float
    two_theta_current_deg: float
    two_theta_cu_ka_deg: float | None
    relative_intensity: float
    multiplicity: int


@dataclass
class XrdAtlasExportPayload:
    phases: list[XrdPhase]
    settings: XrdAtlasSettings
    experimental_patterns: list[ExperimentalPattern] = field(default_factory=list)
