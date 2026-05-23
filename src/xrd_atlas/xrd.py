from __future__ import annotations

from pathlib import Path

import numpy as np
from pymatgen.analysis.diffraction.xrd import XRDCalculator

from .constants import XRD_SOURCE_PRESETS, X_RAY_ENERGY_WAVELENGTH_KEV_A
from .models import CrystalModel, XRDPeakRecord, XRDRequest, XRDResult
from .utils import now_iso, package_versions


def _gaussian_profile(grid: np.ndarray, center: float, fwhm: float) -> np.ndarray:
    sigma = max(fwhm, 1e-6) / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    return np.exp(-0.5 * ((grid - center) / sigma) ** 2)


def _lorentzian_profile(grid: np.ndarray, center: float, fwhm: float) -> np.ndarray:
    gamma = max(fwhm, 1e-6) / 2.0
    return gamma**2 / ((grid - center) ** 2 + gamma**2)


def _peak_profile(grid: np.ndarray, center: float, fwhm: float, model: str) -> np.ndarray:
    if model == "gaussian":
        return _gaussian_profile(grid, center, fwhm)
    if model == "lorentzian":
        return _lorentzian_profile(grid, center, fwhm)
    return 0.5 * _gaussian_profile(grid, center, fwhm) + 0.5 * _lorentzian_profile(grid, center, fwhm)


def build_quality_report(result: XRDResult) -> dict[str, object]:
    warnings = list(result.warnings)
    if result.metadata.get("step_deg", 0.02) > 0.05:
        warnings.append("步长偏大，峰位标定建议使用 step <= 0.02 deg。")
    return {
        "geometry_reliability": "high",
        "intensity_reliability": "medium",
        "model_identity": "theoretical_powder_xrd_profile",
        "convergence_status": "conditionally_stable" if warnings else "stable",
        "warnings": warnings,
        "recommended_next_actions": ["用于峰位标定时优先使用 d、q 和 2θ 表格，不要把相对强度当作精修结果。"],
    }


class XRDService:
    def __init__(self) -> None:
        self._result_cache: dict[tuple[object, ...], XRDResult] = {}

    @staticmethod
    def _build_result_cache_key(crystal: CrystalModel, request: XRDRequest) -> tuple[object, ...]:
        return (crystal.cif_hash, *request.cache_key())

    def resolve_wavelength(self, request: XRDRequest) -> tuple[float, float | None, str]:
        input_mode = request.input_mode or "source"
        if input_mode == "energy":
            if request.energy_keV is None or request.energy_keV <= 0:
                raise ValueError("XRD 能量输入必须为正数，单位 keV。")
            energy = float(request.energy_keV)
            return X_RAY_ENERGY_WAVELENGTH_KEV_A / energy, energy, "energy_keV"

        if input_mode == "wavelength":
            if request.wavelength_A is None or request.wavelength_A <= 0:
                raise ValueError("XRD 波长输入必须为正数，单位 Å。")
            return float(request.wavelength_A), None, "wavelength_A"

        value = XRD_SOURCE_PRESETS.get(request.source_preset)
        if value is None:
            if request.wavelength_A is None or request.wavelength_A <= 0:
                raise ValueError("选择自定义 XRD 光源时必须提供正的波长，单位 Å。")
            return float(request.wavelength_A), None, "source_preset_custom_wavelength"

        wavelength = float(value)
        return wavelength, X_RAY_ENERGY_WAVELENGTH_KEV_A / wavelength, f"source_preset:{request.source_preset}"

    def simulate(self, crystal: CrystalModel, request: XRDRequest) -> XRDResult:
        wavelength, energy_keV, wavelength_source = self.resolve_wavelength(request)
        resolved_request = XRDRequest(
            cif_path=Path(request.cif_path).expanduser().resolve(),
            two_theta_min_deg=request.two_theta_min_deg,
            two_theta_max_deg=request.two_theta_max_deg,
            step_deg=request.step_deg,
            input_mode=request.input_mode or "source",
            source_preset=request.source_preset,
            wavelength_A=wavelength,
            energy_keV=energy_keV if request.input_mode != "wavelength" else request.energy_keV,
            profile_model=request.profile_model,
            fwhm_deg=request.fwhm_deg,
            show_hkl_labels=request.show_hkl_labels,
            show_sticks=request.show_sticks,
        )
        cache_key = self._build_result_cache_key(crystal, resolved_request)
        cached = self._result_cache.get(cache_key)
        if cached is not None:
            return cached

        calculator = XRDCalculator(wavelength=wavelength, symprec=0.0)
        pattern = calculator.get_pattern(
            crystal.pymatgen_structure,
            scaled=True,
            two_theta_range=(resolved_request.two_theta_min_deg, resolved_request.two_theta_max_deg),
        )

        grid = np.arange(
            resolved_request.two_theta_min_deg,
            resolved_request.two_theta_max_deg + resolved_request.step_deg * 0.5,
            resolved_request.step_deg,
            dtype=float,
        )
        profile = np.zeros_like(grid, dtype=float)
        peaks: list[XRDPeakRecord] = []
        pattern_x = np.asarray(pattern.x, dtype=float)
        pattern_y = np.asarray(pattern.y, dtype=float)
        d_hkls = np.asarray(pattern.d_hkls, dtype=float)
        max_intensity = float(np.max(pattern_y)) if pattern_y.size else 1.0

        for two_theta, intensity, families, d_spacing in zip(pattern_x, pattern_y, pattern.hkls, d_hkls, strict=True):
            if not families:
                continue
            representative = max(families, key=lambda item: int(item.get("multiplicity", 1)))
            hkl = tuple(int(value) for value in representative["hkl"])
            multiplicity = sum(int(item.get("multiplicity", 1)) for item in families)
            label = " / ".join(f"({item['hkl'][0]} {item['hkl'][1]} {item['hkl'][2]})" for item in families[:3])
            family_label = "{" + " / ".join(f"{item['hkl'][0]} {item['hkl'][1]} {item['hkl'][2]}" for item in families[:6]) + "}"
            theta = float(two_theta) / 2.0
            g_invA = 0.0 if float(d_spacing) == 0 else 1.0 / float(d_spacing)
            q_invA = 4.0 * np.pi * np.sin(np.deg2rad(theta)) / wavelength
            peaks.append(
                XRDPeakRecord(
                    hkl=hkl,
                    d_spacing_A=float(d_spacing),
                    theta_deg=theta,
                    two_theta_deg=float(two_theta),
                    q_invA=float(q_invA),
                    g_invA=float(g_invA),
                    intensity=float(intensity),
                    normalized_intensity=0.0 if max_intensity == 0 else float(intensity / max_intensity),
                    multiplicity=multiplicity,
                    label=label,
                    family_label=family_label,
                )
            )
            profile += float(intensity) * _peak_profile(grid, float(two_theta), resolved_request.fwhm_deg, resolved_request.profile_model)

        if np.max(profile) > 0:
            profile = profile / np.max(profile) * 100.0

        warnings = list(crystal.validation_report.warnings)
        if crystal.has_partial_occupancy:
            warnings.append("当前 XRD 为平均结构近似；部分占位按 CIF 平均占位参与结构因子计算。")
        if not peaks:
            warnings.append("给定 2θ 范围内没有理论衍射峰。")

        result = XRDResult(
            two_theta_grid=grid,
            intensity_profile=profile,
            stick_positions_deg=pattern_x,
            stick_intensities=pattern_y if pattern_y.size else np.array([], dtype=float),
            peaks=peaks,
            metadata={
                "timestamp_utc": now_iso(),
                "cif_hash": crystal.cif_hash,
                "structure_version_token": crystal.cif_hash,
                "cif_path": str(crystal.cif_path),
                "formula": crystal.formula,
                "space_group_from_cif": crystal.space_group_symbol,
                "space_group_detected": crystal.detected_space_group_symbol,
                "xray_input_mode": resolved_request.input_mode,
                "energy_keV": energy_keV,
                "wavelength_source": wavelength_source,
                "source_preset": resolved_request.source_preset,
                "wavelength_A": wavelength,
                "q_definition": "4*pi*sin(theta)/lambda",
                "two_theta_range_deg": [resolved_request.two_theta_min_deg, resolved_request.two_theta_max_deg],
                "step_deg": resolved_request.step_deg,
                "profile_model": resolved_request.profile_model,
                "fwhm_deg": resolved_request.fwhm_deg,
                "software_versions": package_versions(),
            },
            warnings=warnings,
        )
        result.quality_report = build_quality_report(result)
        result.summary_sections = {
            "实验事实": ["输入为 CIF 结构，输出为理论粉末 XRD。"],
            "计算结果": [f"峰数量：{len(result.peaks)}", f"波长：{wavelength:.5f} Å"],
            "模型限制": ["不包含实验仪器函数拟合、物相检索或 Rietveld 精修。"],
            "可信度提示": result.quality_report["warnings"] or ["峰位几何用于标定较可靠；相对强度仅作理论参考。"],
        }
        self._result_cache[cache_key] = result
        return result
