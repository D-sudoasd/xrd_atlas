from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .constants import DEFAULT_XRD_WAVELENGTH_A, X_RAY_ENERGY_WAVELENGTH_KEV_A
from .models import XRDPeakRecord, XRDRequest, XrdAtlasPeakRow, XrdAtlasSettings, XrdPhase
from .structure import load_crystal_model
from .xrd import XRDService


CU_KA_WAVELENGTH_A = DEFAULT_XRD_WAVELENGTH_A


@dataclass
class XrdAtlasService:
    xrd_service: XRDService

    def __init__(self) -> None:
        self.xrd_service = XRDService()

    def load_phase(self, cif_path: str | Path) -> XrdPhase:
        path = Path(cif_path).expanduser().resolve()
        crystal = load_crystal_model(path)
        return XrdPhase(
            cif_path=path,
            phase_name=path.stem,
            crystal=crystal,
        )

    def load_phases(self, cif_paths: list[str | Path]) -> list[XrdPhase]:
        phases: list[XrdPhase] = []
        for cif_path in cif_paths:
            path = Path(cif_path).expanduser().resolve()
            try:
                phases.append(self.load_phase(path))
            except Exception as exc:
                phases.append(
                    XrdPhase(
                        cif_path=path,
                        phase_name=path.stem,
                        enabled=False,
                        error=str(exc),
                    )
                )
        return phases

    def build_request(self, phase: XrdPhase, settings: XrdAtlasSettings) -> XRDRequest:
        return XRDRequest(
            cif_path=phase.cif_path,
            input_mode=settings.input_mode,
            source_preset=settings.source_preset,
            wavelength_A=settings.wavelength_A,
            energy_keV=settings.energy_keV,
            two_theta_min_deg=settings.two_theta_min_deg,
            two_theta_max_deg=settings.two_theta_max_deg,
            step_deg=settings.step_deg,
            profile_model=settings.profile_model,
            fwhm_deg=settings.fwhm_deg,
            show_hkl_labels=settings.show_labels,
            show_sticks=settings.show_sticks,
        )

    def simulate_phase(self, phase: XrdPhase, settings: XrdAtlasSettings) -> XrdPhase:
        request = self.build_request(phase, settings)
        crystal = load_crystal_model(request.cif_path)
        result = self.xrd_service.simulate(crystal, request)
        phase.crystal = crystal
        phase.result = result
        phase.error = None
        return phase

    def simulate_phases(self, phases: list[XrdPhase], settings: XrdAtlasSettings) -> list[XrdPhase]:
        for phase in phases:
            if not phase.enabled or phase.crystal is None:
                continue
            try:
                self.simulate_phase(phase, settings)
            except Exception as exc:  # pragma: no cover - surfaced in CLI summaries and tests inspect phase.error
                phase.result = None
                phase.error = str(exc)
        return phases


def two_theta_for_wavelength(d_spacing_A: float, wavelength_A: float) -> float | None:
    if d_spacing_A <= 0 or wavelength_A <= 0:
        return None
    argument = wavelength_A / (2.0 * d_spacing_A)
    if argument < -1.0 or argument > 1.0:
        return None
    return float(np.rad2deg(2.0 * np.arcsin(argument)))


def peak_to_atlas_row(phase: XrdPhase, peak: XRDPeakRecord) -> XrdAtlasPeakRow:
    return XrdAtlasPeakRow(
        phase_name=phase.phase_name,
        cif_name=phase.cif_path.name,
        h=peak.hkl[0],
        k=peak.hkl[1],
        l=peak.hkl[2],
        family_label=peak.family_label,
        d_A=peak.d_spacing_A,
        g_1_over_A=peak.g_invA,
        q_1_over_A=peak.q_invA,
        theta_deg=peak.theta_deg,
        two_theta_current_deg=peak.two_theta_deg,
        two_theta_cu_ka_deg=two_theta_for_wavelength(peak.d_spacing_A, CU_KA_WAVELENGTH_A),
        relative_intensity=peak.normalized_intensity,
        multiplicity=peak.multiplicity,
    )


def phase_peak_rows(phase: XrdPhase) -> list[XrdAtlasPeakRow]:
    if phase.result is None:
        return []
    return [peak_to_atlas_row(phase, peak) for peak in phase.result.peaks]
