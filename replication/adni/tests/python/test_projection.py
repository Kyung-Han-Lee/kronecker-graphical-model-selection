from __future__ import annotations


import numpy as np
import pytest

from replication.adni.python import projection as ap

BAND = (ap.BANDPASS_LOW_HZ, ap.BANDPASS_HIGH_HZ)
TR = ap.TR_IN_SECONDS
REAL_P_ROI = 116


def band_limited(n=12, p_roi=8, n_time=187, tr=TR, band=BAND, seed=0):
    lo, hi = band
    rng = np.random.default_rng(seed)
    freqs = np.fft.rfftfreq(n_time, d=tr)
    keep = (freqs >= lo) & (freqs <= hi)
    spec = rng.standard_normal((n, p_roi, len(freqs))) + 1j * rng.standard_normal(
        (n, p_roi, len(freqs))
    )
    spec[..., ~keep] = 0
    return np.fft.irfft(spec, n=n_time, axis=2)


def project(Y, basis):
    return np.einsum("ijt,tc->ijc", Y, basis)


def test_basis_matches_rfft_band_mask():
    n_time, tr = 187, TR
    basis = ap.passband_basis(n_time, tr, *BAND)
    rng = np.random.default_rng(1)
    x = rng.standard_normal(n_time)
    freqs = np.fft.rfftfreq(n_time, d=tr)
    spec = np.fft.rfft(x)
    spec[~((freqs >= BAND[0]) & (freqs <= BAND[1]))] = 0
    want = np.fft.irfft(spec, n=n_time)
    got = basis @ (basis.T @ x)
    assert np.abs(got - want).max() < 1e-10


@pytest.mark.parametrize(
    "n_time,tr,band,expect_cols",
    [
        (8, 1.0, (0.1, 0.5), 7),
        (64, 2.0, (0.05, 0.25), 51),
        (187, TR, BAND, 102),
    ],
)
def test_basis_column_count(n_time, tr, band, expect_cols):
    basis = ap.passband_basis(n_time, tr, *band)
    assert basis.shape == (n_time, expect_cols)
    assert basis.shape[1] <= n_time


def test_basis_is_orthonormal():
    basis = ap.passband_basis(187, TR, *BAND)
    dof = basis.shape[1]
    assert np.abs(basis.T @ basis - np.eye(dof)).max() < ap.ORTHONORMAL_TOL


def test_basis_projection_is_idempotent_and_symmetric():
    basis = ap.passband_basis(187, TR, *BAND)
    P = basis @ basis.T
    assert np.abs(P @ P - P).max() < 1e-10
    assert np.abs(P - P.T).max() < 1e-12


def test_basis_guard_drops_the_null_sine_column():
    n_time, tr, band = 64, 2.0, (0.05, 0.25)
    guarded = ap.passband_basis(n_time, tr, *band)
    freqs = np.fft.rfftfreq(n_time, d=tr)
    t = np.arange(n_time) * tr
    unguarded = []
    for f in freqs:
        if not (band[0] <= f <= band[1]):
            continue
        unguarded += [np.cos(2 * np.pi * f * t), np.sin(2 * np.pi * f * t)]
    raw = np.column_stack(unguarded)
    assert raw.shape[1] == guarded.shape[1] + 1
    assert np.linalg.matrix_rank(raw) == guarded.shape[1]


def test_basis_rejects_dc_band():
    with pytest.raises(ValueError, match="고역 하한"):
        ap.passband_basis(187, TR, 0.0, 0.1)


def test_basis_rejects_empty_band():
    with pytest.raises(ValueError, match="드는 성분이 없다"):
        ap.passband_basis(187, TR, 0.40, 0.45)


def test_verify_passes_on_clean_projection():
    Y = band_limited(n=20, p_roi=30)
    basis = ap.passband_basis(Y.shape[2], TR, *BAND)
    assert ap.verify(basis, Y, project(Y, basis)) == []


def test_verify_pins_in_band_power_to_unity():
    Y = band_limited(n=20, p_roi=30)
    basis = ap.passband_basis(Y.shape[2], TR, *BAND)
    Yp = project(Y, basis)
    keep = float((Yp**2).mean() * basis.shape[1] / Y.shape[2]) / ap.in_band_mean_square(
        Y, TR, *BAND
    )
    assert abs(keep - 1.0) < ap.IN_BAND_KEEP_TOL


def test_verify_rejects_non_orthonormal_basis():
    Y = band_limited(n=20, p_roi=30)
    basis = ap.passband_basis(Y.shape[2], TR, *BAND)
    bad = basis * 1.05
    problems = ap.verify(bad, Y, project(Y, bad))
    assert problems and "정규직교" in problems[0]


def test_verify_rejects_basis_missing_an_in_band_bin():
    Y = band_limited(n=20, p_roi=30)
    full = ap.passband_basis(Y.shape[2], TR, *BAND)
    short = full[:, :-2]
    problems = ap.verify(short, Y, project(Y, short))
    assert problems and any("잔차" in p for p in problems)


def test_verify_rejects_zero_tensor():
    Y = np.zeros((20, REAL_P_ROI, 187))
    basis = ap.passband_basis(187, TR, *BAND)
    problems = ap.verify(basis, Y, project(Y, basis))
    assert problems and "퇴화" in problems[0]


def test_verify_rejects_near_degenerate_at_real_size():
    rng = np.random.default_rng(3)
    base = band_limited(n=1, p_roi=REAL_P_ROI)
    Y = np.repeat(base, 20, axis=0) + 1e-12 * rng.standard_normal((20, REAL_P_ROI, 187))
    basis = ap.passband_basis(187, TR, *BAND)
    Yp = project(Y, basis)
    assert ap.spectrum_summary(Yp)["degenerate"] is True
    problems = ap.verify(basis, Y, Yp)
    assert problems and "퇴화" in problems[0]


def test_verify_honours_explicit_band_arguments():
    band = (0.02, 0.07)
    Y = band_limited(n=20, p_roi=30, band=band)
    basis = ap.passband_basis(Y.shape[2], TR, *band)
    assert ap.verify(basis, Y, project(Y, basis), TR, *band) == []


def test_config_values_match_preprocessing_package():
    ap.check_matches_config()


def test_component_count_agrees_with_probe_arm():
    n_comp = ap.passband_basis(187, TR, *BAND).shape[1]
    assert n_comp == 102, "R 쪽 ARM_SPECS$projected$p_time 을 같이 고쳐야 한다"


def test_spectrum_summary_reports_zero_near_null_after_projection():
    Y = band_limited(n=20, p_roi=30)
    basis = ap.passband_basis(Y.shape[2], TR, *BAND)
    assert ap.spectrum_summary(Y)["n_near_null"] > 0
    assert ap.spectrum_summary(project(Y, basis))["n_near_null"] == 0


def test_occupied_output_is_rejected_before_loading_inputs(tmp_path):
    output = tmp_path / "occupied"
    output.mkdir()
    with pytest.raises(FileExistsError):
        ap.run(tmp_path / "missing", output, groups=("CN",))


def test_failed_projection_does_not_publish_partial_output(tmp_path):
    source = tmp_path / "source"
    source.mkdir()
    np.save(source / "Y_CN_combat.npy", np.zeros((20, REAL_P_ROI, 187)))
    output = tmp_path / "projected"
    with pytest.raises(ValueError, match="projection validation"):
        ap.run(source, output, groups=("CN",))
    assert not output.exists()


def test_valid_projection_writes_the_tensor_and_manifest(tmp_path):
    import json

    source = tmp_path / "source"
    source.mkdir()
    data = band_limited(n=20, p_roi=30)
    np.save(source / "Y_CN_combat.npy", data)
    output = tmp_path / "projected"
    ap.run(source, output, groups=("CN",))
    expected = project(data, ap.passband_basis(187, TR, *BAND))
    np.testing.assert_array_equal(np.load(output / "Y_CN_combat.npy"), expected)
    manifest = json.loads((output / "resample_manifest.json").read_text())
    assert manifest["groups"][0]["problems"] == []
