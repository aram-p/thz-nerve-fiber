"""Independent quasi-static check of the node field (no COMSOL).

Axisymmetric finite-volume solver for div(eps grad phi) = 0 with complex eps,
uniform applied axial field E0 (phi = -E0 z on the outer boundary). Geometry:
axon (r < a), myelin tube (a < r < b) along the whole domain except a node gap
|z| < Ln/2 filled with node material, outer medium for r > b.

Returns the node field factor F = <|E|^2>_node / |E0|^2 and absorbed powers.
Valid when every size is << the wavelength in water (node, axon, myelin): the
node absorbed power is then P = 0.5 w eps0 eps''_node |E0|^2 V_node F.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


def graded_edges(breaks, h_fine, h_max, growth=1.15):
    """1D edges with spacing h_fine at every break point, growing geometrically to h_max."""
    breaks = np.unique(np.asarray(breaks, dtype=float))
    edges = [breaks[0]]
    for lo, hi in zip(breaks[:-1], breaks[1:]):
        seg = [lo]
        # grow from both ends towards the middle
        left, right = [lo], [hi]
        h = h_fine
        while left[-1] + h < (lo + hi) / 2:
            left.append(left[-1] + h)
            h = min(h * growth, h_max)
        h = h_fine
        while right[-1] - h > (lo + hi) / 2:
            right.append(right[-1] - h)
            h = min(h * growth, h_max)
        seg = left + right[::-1]
        edges.extend(seg[1:])
    return np.unique(np.asarray(edges))


@dataclass
class QSResult:
    F_node: float
    F_myelin: float
    F_axon_at_node: float
    r: np.ndarray
    z: np.ndarray
    E2: np.ndarray


def solve_node_qs(a, b, Ln, eps_axon, eps_myelin, eps_node, eps_out,
                  R=None, Z=None, h_fine=None, h_max=None, field="axial"):
    """Solve the quasi-static node problem. Lengths in any consistent unit.

    field="axial" applies E0 along the fibre axis (E parallel to fibre, m = 0).
    field="transverse" applies E0 along x (E perpendicular to fibre): potential
    psi(r, z) cos(phi) with psi = -E0 r far away, and an extra eps m^2/r^2 term (m = 1).
    F is then averaged over phi: <|E|^2> = (|d_r psi|^2 + |psi/r|^2 + |d_z psi|^2) / 2.
    """
    if field not in ("axial", "transverse"):
        raise ValueError(field)
    m = 0 if field == "axial" else 1
    t = b - a
    R = R or 12 * b
    Z = Z or max(12 * b, 20 * Ln)
    h_fine = h_fine or min(Ln, t, a) / 20
    h_max = h_max or b / 2
    re = graded_edges([0, a, b, R], h_fine, h_max)
    ze = graded_edges([-Z, -Ln / 2, 0, Ln / 2, Z], h_fine, h_max)
    rc, zc = 0.5 * (re[1:] + re[:-1]), 0.5 * (ze[1:] + ze[:-1])
    dr, dz = np.diff(re), np.diff(ze)
    nr, nz = len(rc), len(zc)
    RC, ZC = np.meshgrid(rc, zc, indexing="ij")
    eps = np.full((nr, nz), eps_out, dtype=complex)
    eps[RC < a] = eps_axon
    shell = (RC > a) & (RC < b)
    eps[shell] = eps_myelin
    eps[shell & (np.abs(ZC) < Ln / 2)] = eps_node

    idx = np.arange(nr * nz).reshape(nr, nz)
    rows, cols, vals = [], [], []
    rhs = np.zeros(nr * nz, dtype=complex)
    diag = np.zeros(nr * nz, dtype=complex)

    def phi_bc(r, z):
        return -z if m == 0 else -r * np.ones_like(z)  # E0 = 1 along +z or +x

    # r-faces between i and i+1
    for i in range(nr - 1):
        rf = re[i + 1]
        area = 2 * np.pi * rf * dz
        d = rc[i + 1] - rc[i]
        ef = 2 * eps[i] * eps[i + 1] / (eps[i] + eps[i + 1])
        g = ef * area / d
        p, q = idx[i], idx[i + 1]
        rows += [p, q]; cols += [q, p]; vals += [-g, -g]
        np.add.at(diag, p, g); np.add.at(diag, q, g)
    # outer r boundary (Dirichlet)
    i = nr - 1
    g = eps[i] * 2 * np.pi * re[-1] * dz / (re[-1] - rc[i])
    np.add.at(diag, idx[i], g)
    rhs[idx[i]] += g * phi_bc(re[-1], zc)
    # z-faces
    ring = np.pi * (re[1:] ** 2 - re[:-1] ** 2)
    if m:
        np.add.at(diag, idx, eps * m ** 2 / RC ** 2 * np.outer(ring, dz))
    for j in range(nz - 1):
        d = zc[j + 1] - zc[j]
        ef = 2 * eps[:, j] * eps[:, j + 1] / (eps[:, j] + eps[:, j + 1])
        g = ef * ring / d
        p, q = idx[:, j], idx[:, j + 1]
        rows += [p, q]; cols += [q, p]; vals += [-g, -g]
        np.add.at(diag, p, g); np.add.at(diag, q, g)
    for j, zb in ((0, ze[0]), (nz - 1, ze[-1])):
        g = eps[:, j] * ring / abs(zb - zc[j])
        np.add.at(diag, idx[:, j], g)
        rhs[idx[:, j]] += g * phi_bc(rc, zb)
    rows = np.concatenate([np.atleast_1d(x) for x in rows] + [np.arange(nr * nz)])
    cols = np.concatenate([np.atleast_1d(x) for x in cols] + [np.arange(nr * nz)])
    vals = np.concatenate([np.atleast_1d(x) for x in vals] + [diag])
    A = sp.csr_matrix((vals, (rows, cols)), shape=(nr * nz, nr * nz))
    phi = spla.spsolve(A.tocsc(), rhs).reshape(nr, nz)

    # cell-centred gradient (one-sided at edges)
    Er = -np.gradient(phi, rc, axis=0)
    Ez = -np.gradient(phi, zc, axis=1)
    E2 = np.abs(Er) ** 2 + np.abs(Ez) ** 2
    if m:
        E2 = 0.5 * (E2 + np.abs(phi / RC) ** 2)
    vol = np.outer(ring, dz)
    node = shell & (np.abs(ZC) < Ln / 2)
    my = shell & (np.abs(ZC) > Ln / 2 + 5 * t) & (np.abs(ZC) < Z / 2)
    ax = (RC < a) & (np.abs(ZC) < Ln / 2)
    F = lambda m: float(np.sum(E2[m] * vol[m]) / np.sum(vol[m]))
    return QSResult(F(node), F(my), F(ax), rc, zc, E2)
