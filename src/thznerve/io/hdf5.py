"""Per-simulation HDF5 storage (ADR-0004).

One HDF5 file per simulation, with:
    /scalars                  group, scalar metrics stored as attrs
    /axial_profile/z          1D dataset
    /axial_profile/E_mag      1D dataset
    /axial_slice/x            1D dataset
    /axial_slice/z            1D dataset
    /axial_slice/E_mag        2D dataset, shape (len(x), len(z))
    /field_3d/{x,y,z,E_mag}   optional, opt-in via sweep config
    /metadata                 group, provenance stored as attrs
                              (params is a JSON string attr)
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import h5py
import numpy as np


@dataclass
class AxialProfile:
    z: np.ndarray
    e_mag: np.ndarray


@dataclass
class AxialSlice:
    x: np.ndarray
    z: np.ndarray
    e_mag: np.ndarray  # shape (len(x), len(z))


@dataclass
class Field3D:
    x: np.ndarray
    y: np.ndarray
    z: np.ndarray
    e_mag: np.ndarray  # shape (len(x), len(y), len(z))


@dataclass
class SimResult:
    scalars: dict[str, float]
    axial_profile: AxialProfile
    axial_slice: AxialSlice
    metadata: dict[str, Any]
    field_3d: Field3D | None = None


def write_result(
    path: Path | str,
    *,
    scalars: dict[str, float],
    axial_profile: tuple[np.ndarray, np.ndarray] | AxialProfile,
    axial_slice: tuple[np.ndarray, np.ndarray, np.ndarray] | AxialSlice,
    metadata: dict[str, Any],
    field_3d: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray] | Field3D | None = None,
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    if isinstance(axial_profile, tuple):
        axial_profile = AxialProfile(*axial_profile)
    if isinstance(axial_slice, tuple):
        axial_slice = AxialSlice(*axial_slice)
    if field_3d is not None and isinstance(field_3d, tuple):
        field_3d = Field3D(*field_3d)

    with h5py.File(path, "w") as f:
        g_scalars = f.create_group("scalars")
        for k, v in scalars.items():
            g_scalars.attrs[k] = float(v)

        g_prof = f.create_group("axial_profile")
        g_prof.create_dataset("z", data=np.asarray(axial_profile.z))
        g_prof.create_dataset("E_mag", data=np.asarray(axial_profile.e_mag))

        g_slice = f.create_group("axial_slice")
        g_slice.create_dataset("x", data=np.asarray(axial_slice.x))
        g_slice.create_dataset("z", data=np.asarray(axial_slice.z))
        g_slice.create_dataset("E_mag", data=np.asarray(axial_slice.e_mag))

        if field_3d is not None:
            g_3d = f.create_group("field_3d")
            g_3d.create_dataset("x", data=np.asarray(field_3d.x))
            g_3d.create_dataset("y", data=np.asarray(field_3d.y))
            g_3d.create_dataset("z", data=np.asarray(field_3d.z))
            g_3d.create_dataset("E_mag", data=np.asarray(field_3d.e_mag))

        g_meta = f.create_group("metadata")
        for k, v in metadata.items():
            # dicts/lists get JSON-encoded; primitives go in directly
            if isinstance(v, (dict, list)):
                g_meta.attrs[k] = json.dumps(v, sort_keys=True)
            else:
                g_meta.attrs[k] = v


def read_result(path: Path | str) -> SimResult:
    path = Path(path)

    with h5py.File(path, "r") as f:
        scalars = {k: float(v) for k, v in f["scalars"].attrs.items()}

        axial_profile = AxialProfile(
            z=f["axial_profile/z"][...],
            e_mag=f["axial_profile/E_mag"][...],
        )
        axial_slice = AxialSlice(
            x=f["axial_slice/x"][...],
            z=f["axial_slice/z"][...],
            e_mag=f["axial_slice/E_mag"][...],
        )

        field_3d: Field3D | None = None
        if "field_3d" in f:
            field_3d = Field3D(
                x=f["field_3d/x"][...],
                y=f["field_3d/y"][...],
                z=f["field_3d/z"][...],
                e_mag=f["field_3d/E_mag"][...],
            )

        metadata: dict[str, Any] = {}
        for k, v in f["metadata"].attrs.items():
            # Attempt to JSON-decode dict/list-shaped strings, else passthrough
            if isinstance(v, (bytes, np.bytes_)):
                v = v.decode()
            if isinstance(v, str) and v and v[0] in "{[":
                try:
                    metadata[k] = json.loads(v)
                    continue
                except json.JSONDecodeError:
                    pass
            metadata[k] = v

    return SimResult(
        scalars=scalars,
        axial_profile=axial_profile,
        axial_slice=axial_slice,
        metadata=metadata,
        field_3d=field_3d,
    )
