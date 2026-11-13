from __future__ import annotations

import logging
import platform
from pathlib import Path

import numpy as np

from replication.io import check_output, output_directory, write_json

log = logging.getLogger(__name__)

PROJECT_ROOT = Path(__file__).resolve().parents[1]
AGGREGATED_DIR = PROJECT_ROOT / "data/processed/antspy/aggregated"

GROUPS = ("CN", "MCI", "Dementia")

TR_IN_SECONDS = 3.0
BANDPASS_LOW_HZ = 0.01
BANDPASS_HIGH_HZ = 0.1

ORTHONORMAL_TOL = 1e-10
IN_BAND_KEEP_TOL = 1e-12
RESIDUAL_IN_BAND_MAX = 1e-6
NEAR_NULL_RATIO = 1e-2
DEGENERATE_REL_TOL = 1e-9


def check_matches_config() -> None:
    from . import settings as pp_config

    mismatches = [
        f"{name}: this script {here} vs config {there}"
        for name, here, there in (
            ("TR", TR_IN_SECONDS, pp_config.TR_SECONDS),
            ("bandpass low", BANDPASS_LOW_HZ, pp_config.BANDPASS_LOW),
            ("bandpass high", BANDPASS_HIGH_HZ, pp_config.BANDPASS_HIGH),
        )
        if not np.isclose(here, there)
    ]
    if mismatches:
        raise RuntimeError(
            "adni.preprocess.config 와 어긋난다:\n  " + "\n  ".join(mismatches)
        )


def passband_basis(n_time: int, tr: float, lo: float, hi: float) -> np.ndarray:
    if lo <= 0.0:
        raise ValueError(
            f"고역 하한이 0 이하다 ({lo}). DC 를 포함한 대역은 아직 지원하지 "
            "않는다 — verify 조건 2 의 평균 처리가 사영과 어긋난다"
        )
    freqs = np.fft.rfftfreq(n_time, d=tr)
    t = np.arange(n_time) * tr
    cols = []
    for k, f in enumerate(freqs):
        if not (lo <= f <= hi):
            continue
        omega = 2.0 * np.pi * f
        cols.append(np.cos(omega * t))
        real_only = k == 0 or (n_time % 2 == 0 and k == n_time // 2)
        if not real_only:
            cols.append(np.sin(omega * t))
    if not cols:
        raise ValueError(
            f"통과대역 [{lo}, {hi}] Hz 에 드는 성분이 없다 (n_time={n_time}, tr={tr})"
        )
    basis, _ = np.linalg.qr(np.column_stack(cols))
    return basis


def time_axis_gram(Y: np.ndarray) -> np.ndarray:
    n, p_roi, _ = Y.shape
    Yc = Y - Y.mean(axis=0, keepdims=True)
    return np.einsum("ijk,ijl->kl", Yc, Yc) / (n * p_roi)


def in_band_mean_square(Y: np.ndarray, tr: float, lo: float, hi: float) -> float:
    n_time = Y.shape[2]
    freqs = np.fft.rfftfreq(n_time, d=tr)
    spec = np.abs(np.fft.rfft(Y - Y.mean(axis=2, keepdims=True), axis=2)) ** 2
    spec = spec.mean(axis=(0, 1))
    weight = np.full(freqs.shape, 2.0)
    weight[0] = 1.0
    if n_time % 2 == 0:
        weight[-1] = 1.0
    mask = (freqs >= lo) & (freqs <= hi)
    return float((weight[mask] * spec[mask]).sum() / n_time**2)


def spectrum_summary(Y: np.ndarray) -> dict:
    ev = np.clip(np.linalg.eigvalsh(time_axis_gram(Y))[::-1], 0.0, None)
    lam_max, lam_min = float(ev[0]), float(ev[-1])
    scale = float((Y**2).mean())
    degenerate = not np.isfinite(lam_max) or lam_max <= DEGENERATE_REL_TOL * scale
    return {
        "p_second_mode": int(Y.shape[2]),
        "degenerate": bool(degenerate),
        "cond": float(lam_max / lam_min) if lam_min > 0 else float("inf"),
        "lambda_max": lam_max,
        "lambda_min": lam_min,
        "pr_effective_rank": (
            float(ev.sum() ** 2 / (ev**2).sum()) if not degenerate else None
        ),
        "n_near_null": (
            int((ev < NEAR_NULL_RATIO * lam_max).sum()) if not degenerate else None
        ),
    }


def verify(
    basis: np.ndarray,
    Y: np.ndarray,
    Y_proj: np.ndarray,
    tr: float = TR_IN_SECONDS,
    lo: float = BANDPASS_LOW_HZ,
    hi: float = BANDPASS_HIGH_HZ,
) -> list[str]:
    problems = []
    if spectrum_summary(Y)["degenerate"] or spectrum_summary(Y_proj)["degenerate"]:
        problems.append(
            "스펙트럼이 퇴화다 (Gram 최대 고유값이 자료 규모 대비 0). "
            "피험자 간 정보가 없다"
        )
        return problems

    dof = basis.shape[1]
    gram_err = np.abs(basis.T @ basis - np.eye(dof)).max()
    if not (gram_err <= ORTHONORMAL_TOL):
        problems.append(f"기저가 정규직교가 아니다 (최대 이탈 {gram_err:.3e})")

    before = in_band_mean_square(Y, tr, lo, hi)
    if not (before > 0):
        problems.append("입력의 대역내 파워가 0 이다 — 퇴화 입력")
        return problems

    proj_ms = float((Y_proj**2).mean() * dof / Y.shape[2])
    keep = proj_ms / before
    if not (abs(keep - 1.0) <= IN_BAND_KEEP_TOL):
        problems.append(
            f"대역내 파워 보존비가 {keep:.6g} 로 1 에서 "
            f"{abs(keep - 1.0):.3g} 벗어났다 (허용 {IN_BAND_KEEP_TOL})"
        )

    residual = Y - np.einsum("ijc,tc->ijt", Y_proj, basis)
    resid_in_band = in_band_mean_square(residual, tr, lo, hi) / before
    if not (resid_in_band <= RESIDUAL_IN_BAND_MAX):
        problems.append(
            f"버린 잔차에 대역내 파워가 {resid_in_band:.3e} 남았다 "
            f"(허용 {RESIDUAL_IN_BAND_MAX}). 기저가 통과대역을 못 덮었다"
        )

    after = spectrum_summary(Y_proj)
    if after["n_near_null"] > 0:
        problems.append(
            f"사영 뒤에도 빈 방향이 {after['n_near_null']} 개다. "
            "사영의 존재 이유가 그것을 없애는 것이다"
        )
    return problems


def process_group(group: str, in_dir: Path, out_dir: Path) -> dict:
    src = in_dir / f"Y_{group}_combat.npy"
    if not src.exists():
        raise FileNotFoundError(f"입력 텐서가 없다: {src}")
    Y = np.load(src).astype(np.float64)
    if Y.ndim != 3:
        raise ValueError(f"(n, p_roi, T) 3 차원이어야 한다. 받은 것: {Y.shape}")

    basis = passband_basis(Y.shape[2], TR_IN_SECONDS, BANDPASS_LOW_HZ, BANDPASS_HIGH_HZ)
    Y_proj = np.einsum("ijt,tc->ijc", Y, basis)

    before = spectrum_summary(Y)
    after = spectrum_summary(Y_proj)
    problems = verify(basis, Y, Y_proj)

    log.info(
        "%s: p %d -> %d | 빈 방향 %d -> %d | cond %.3g -> %.3g | PR %.1f -> %.1f",
        group,
        before["p_second_mode"],
        after["p_second_mode"],
        before["n_near_null"],
        after["n_near_null"],
        before["cond"],
        after["cond"],
        before["pr_effective_rank"],
        after["pr_effective_rank"],
    )

    dst = out_dir / f"Y_{group}_combat.npy"
    saved = not problems
    if not saved:
        log.error("%s: 검증 실패라 저장하지 않는다 — %s", group, "; ".join(problems))
    else:
        out_dir.mkdir(parents=True, exist_ok=True)
        np.save(dst, Y_proj)

    return {
        "group": group,
        "source": src.name,
        "output": dst.name if saved else None,
        "n_subjects": int(Y.shape[0]),
        "p_roi": int(Y.shape[1]),
        "n_components": int(basis.shape[1]),
        "before": before,
        "after": after,
        "problems": problems,
    }


def run(input_dir=None, output_dir=None, groups=GROUPS, option="C", out_suffix="pb"):
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    check_matches_config()
    if output_dir is not None:
        check_output(output_dir)
    if (
        not groups
        or len(set(groups)) != len(groups)
        or not set(groups).issubset(GROUPS)
    ):
        raise ValueError("groups must be distinct ADNI group names")

    in_dir = (
        Path(input_dir) if input_dir is not None else AGGREGATED_DIR / option / "combat"
    )
    groups = tuple(groups)

    n_times = {}
    for group in groups:
        src = in_dir / f"Y_{group}_combat.npy"
        if not src.exists():
            raise FileNotFoundError(f"입력 텐서가 없다: {src}")
        n_times[group] = int(np.load(src, mmap_mode="r").shape[2])
    if len(set(n_times.values())) != 1:
        raise ValueError(
            f"군마다 시간 길이가 다르다: {n_times}. 한 디렉토리에 폭이 다른 "
            "배열이 섞이면 하류가 어느 쪽을 믿을지 알 수 없다"
        )
    n_comp = passband_basis(
        next(iter(n_times.values())),
        TR_IN_SECONDS,
        BANDPASS_LOW_HZ,
        BANDPASS_HIGH_HZ,
    ).shape[1]
    out_dir = (
        Path(output_dir)
        if output_dir is not None
        else AGGREGATED_DIR / option / f"combat_{out_suffix}{n_comp}"
    )

    manifest = {
        "transform": "passband_projection",
        "tr_in_seconds": TR_IN_SECONDS,
        "bandpass_hz": [BANDPASS_LOW_HZ, BANDPASS_HIGH_HZ],
        "n_components": n_comp,
        "near_null_ratio": NEAR_NULL_RATIO,
        "option": option,
        "package_versions": {
            "numpy": np.__version__,
            "python": platform.python_version(),
        },
        "groups": [],
    }
    with output_directory(out_dir) as staging:
        for group in groups:
            record = process_group(group, in_dir, staging)
            if record["problems"]:
                raise ValueError(
                    f"projection validation failed for {group}: {record['problems']}"
                )
            for stage in ("before", "after"):
                if not np.isfinite(record[stage]["cond"]):
                    record[stage]["cond"] = None
            manifest["groups"].append(record)
        write_json(staging / "resample_manifest.json", manifest)
    return out_dir
