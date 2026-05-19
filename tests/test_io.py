"""Round-trip tests for the IO layer (HDF5 + manifest)."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from thznerve.io.hdf5 import (
    AxialProfile,
    AxialSlice,
    Field3D,
    read_result,
    write_result,
)
from thznerve.io.manifest import MANIFEST_COLUMNS, append, query
from thznerve.io.provenance import config_hash


def _sample_arrays():
    z = np.linspace(0.0, 240.0, 100)
    e_profile = np.cos(0.1 * z)
    x = np.linspace(-20.0, 20.0, 50)
    e_slice = np.outer(np.ones_like(x), e_profile)
    return z, e_profile, x, e_slice


def test_hdf5_roundtrip_2d_only(tmp_path: Path):
    z, e_profile, x, e_slice = _sample_arrays()
    scalars = {"peak_e_at_node": 1.234, "abs_power_node": 0.456}
    metadata = {
        "comsol_version": "6.3",
        "git_sha": "abc123",
        "config_hash": "deadbeef0000",
        "timestamp_iso": "2026-05-19T10:00:00",
        "params": {"frequency_THz": 0.6, "node_sigma_S_per_m": 0.0},
    }

    path = tmp_path / "result.h5"
    write_result(
        path,
        scalars=scalars,
        axial_profile=(z, e_profile),
        axial_slice=(x, z, e_slice),
        metadata=metadata,
    )

    result = read_result(path)
    assert result.scalars == pytest.approx(scalars)
    np.testing.assert_allclose(result.axial_profile.z, z)
    np.testing.assert_allclose(result.axial_profile.e_mag, e_profile)
    np.testing.assert_allclose(result.axial_slice.x, x)
    np.testing.assert_allclose(result.axial_slice.z, z)
    np.testing.assert_allclose(result.axial_slice.e_mag, e_slice)
    assert result.field_3d is None
    assert result.metadata["comsol_version"] == "6.3"
    assert result.metadata["params"] == metadata["params"]


def test_hdf5_roundtrip_with_3d(tmp_path: Path):
    z, e_profile, x, e_slice = _sample_arrays()
    y = np.linspace(-20.0, 20.0, 20)
    field3d = np.random.default_rng(0).random((len(x), len(y), len(z)))

    path = tmp_path / "result_3d.h5"
    write_result(
        path,
        scalars={"peak": 1.0},
        axial_profile=AxialProfile(z=z, e_mag=e_profile),
        axial_slice=AxialSlice(x=x, z=z, e_mag=e_slice),
        metadata={"git_sha": "abc"},
        field_3d=Field3D(x=x, y=y, z=z, e_mag=field3d),
    )

    result = read_result(path)
    assert result.field_3d is not None
    np.testing.assert_allclose(result.field_3d.e_mag, field3d)


def test_manifest_append_and_query(tmp_path: Path):
    path = tmp_path / "manifest.csv"

    base = {
        "sweep_name": "freq_sweep_001",
        "params_json": json.dumps({"f": 0.6}),
        "output_path": "results/freq_sweep_001/001.h5",
        "comsol_version": "6.3",
        "git_sha": "abc123",
        "config_hash": "h1",
        "timestamp": "2026-05-19T10:00:00",
        "status": "ok",
    }
    append(path, {**base, "sim_id": "001"})
    append(path, {**base, "sim_id": "002", "status": "failed"})
    append(path, {**base, "sim_id": "003"})

    df_all = query(path)
    assert len(df_all) == 3
    assert list(df_all.columns) == MANIFEST_COLUMNS

    df_ok = query(path, status="ok")
    assert len(df_ok) == 2
    assert set(df_ok["sim_id"]) == {"001", "003"}


def test_manifest_rejects_unknown_columns(tmp_path: Path):
    path = tmp_path / "manifest.csv"
    bad = {col: "x" for col in MANIFEST_COLUMNS}
    bad["sim_id"] = "001"
    bad["surprise"] = "extra"
    with pytest.raises(ValueError, match="unknown columns"):
        append(path, bad)


def test_config_hash_is_order_independent():
    a = {"name": "x", "model": {"r": 5, "n": 40}}
    b = {"model": {"n": 40, "r": 5}, "name": "x"}
    assert config_hash(a) == config_hash(b)
