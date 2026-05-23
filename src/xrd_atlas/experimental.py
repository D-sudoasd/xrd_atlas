from __future__ import annotations

from pathlib import Path

import numpy as np

from .models import ExperimentalPattern, XrdAxisMode


def _read_numeric_table(path: Path) -> np.ndarray:
    attempts = [
        {"delimiter": ","},
        {"delimiter": None},
        {"delimiter": "\t"},
        {"delimiter": ";"},
    ]
    for kwargs in attempts:
        data = np.genfromtxt(
            path,
            comments="#",
            invalid_raise=False,
            dtype=float,
            **kwargs,
        )
        data = np.asarray(data, dtype=float)
        if data.ndim == 1 and data.size >= 2:
            data = data.reshape(-1, data.size)
        if data.ndim == 2 and data.shape[1] >= 2:
            data = data[:, :2]
            data = data[np.all(np.isfinite(data), axis=1)]
            if len(data):
                return data
    raise ValueError("实验谱线必须包含至少两列可解析的数值数据。")


def load_experimental_pattern(path: str | Path, axis_mode: XrdAxisMode) -> ExperimentalPattern:
    resolved = Path(path).expanduser().resolve()
    data = _read_numeric_table(resolved)
    y = np.asarray(data[:, 1], dtype=float)
    if np.max(np.abs(y)) > 0:
        y = y / np.max(np.abs(y)) * 100.0
    return ExperimentalPattern(
        path=resolved,
        label=resolved.stem,
        x_values=np.asarray(data[:, 0], dtype=float),
        intensity=y,
        axis_mode=axis_mode,
    )
