"""Quasi-static checks of the single-fibre node law that a referee will ask for.

layered_cylinder_factor: <|E|^2>/|E0|^2 in one layer of an infinite N-layer cylinder in a host,
    E perpendicular to the axis (exact 2D series, potential (A r + B/r) cos(theta) per layer).
    Used for the axon membrane (5 nm lipid bilayer) at the node.
periodic_shell_factor: the same mean in the myelin shell when the fibre sits in a square
    lattice at area fraction p (dense white matter), normalised to the macroscopic field.
    Finite-volume Laplace solver on the periodic cell, potential phi = -x + psi(x, y).
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


def layered_cylinder_factor(radii, eps_layers, eps_host, layer):
    """radii: increasing outer radii r_1..r_N of the N layers (layer 0 is the core).
    Returns <|E|^2>/|E0|^2 averaged over the area of layer `layer` (0-based)."""
    n = len(radii)
    eps = list(eps_layers) + [eps_host]
    # unknowns: A_0, (A_i, B_i) for i = 1..N-1, B_N (host, A_N = -1 for E0 = 1 along x)
    nu = 1 + 2 * (n - 1) + 1
    M = np.zeros((2 * n, nu), complex); rhs = np.zeros(2 * n, complex)

    def idx(i):            # column of (A_i, B_i); None where the coefficient is fixed
        if i == 0:
            return 0, None
        if i == n:
            return None, nu - 1
        return 1 + 2 * (i - 1), 2 + 2 * (i - 1)

    for k, r in enumerate(radii):          # interface between layer k and k+1
        for side, sgn in ((k, 1), (k + 1, -1)):
            ia, ib = idx(side)
            e = eps[side]
            # continuity of phi: A r + B / r ; of eps dphi/dr: eps (A - B / r^2)
            if ia is not None:
                M[2 * k, ia] += sgn * r; M[2 * k + 1, ia] += sgn * e
            if ib is not None:
                M[2 * k, ib] += sgn / r; M[2 * k + 1, ib] -= sgn * e / r ** 2
        if k + 1 == n:                     # host: A_N = -1 moves to the right-hand side
            rhs[2 * k] -= -1 * (-1) * r; rhs[2 * k + 1] -= -1 * (-1) * eps[n]
    x = np.linalg.solve(M, rhs)
    ia, ib = idx(layer)
    A = x[ia] if ia is not None else -1.0
    B = x[ib] if ib is not None else 0.0
    r0 = 0.0 if layer == 0 else radii[layer - 1]
    r = np.linspace(max(r0, 1e-12), radii[layer], 2001)
    q = (abs(A) ** 2 + abs(B) ** 2 / r ** 4) * r
    return float(np.trapezoid(q, r) / np.trapezoid(r, r))


def periodic_shell_factor(a, b, eps_ax, eps_my, eps_w, fraction, N=400):
    """<|E|^2> in the myelin shell / |<E>|^2 for a square lattice of coated cylinders, E along x."""
    L = np.sqrt(np.pi * b ** 2 / fraction)
    h = L / N
    c = (np.arange(N) + 0.5) * h - L / 2
    X, Y = np.meshgrid(c, c, indexing="ij")
    R = np.hypot(X, Y)
    eps = np.where(R < a, eps_ax, np.where(R < b, eps_my, eps_w)).astype(complex)
    # face permittivities (harmonic mean), periodic neighbours
    ex = 2 * eps * np.roll(eps, -1, 0) / (eps + np.roll(eps, -1, 0))     # face i+1/2, j
    ey = 2 * eps * np.roll(eps, -1, 1) / (eps + np.roll(eps, -1, 1))
    idx = np.arange(N * N).reshape(N, N)
    ip, jp = np.roll(idx, -1, 0), np.roll(idx, -1, 1)
    rows, cols, vals = [], [], []
    for nb, e in ((ip, ex), (jp, ey)):
        for (p, q) in ((idx, nb), (nb, idx)):
            rows += [p.ravel(), p.ravel()]; cols += [p.ravel(), q.ravel()]
            vals += [e.ravel(), -e.ravel()]
    A = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(N * N,) * 2)
    # phi = -x + psi: flux of -x across x-faces is -eps_face * h / h * ... -> source term
    flux_x = ex * (-1.0) * h          # eps * d(-x)/dx * face length (h) ... per unit depth
    src = flux_x - np.roll(flux_x, 1, 0)          # -div of the applied-field flux (A psi = -div eps grad psi)
    A = A.tolil(); A[0, :] = 0; A[0, 0] = 1; A = A.tocsc()      # fix psi gauge
    b_ = src.ravel().astype(complex); b_[0] = 0
    psi = spla.spsolve(A, b_).reshape(N, N)
    # cell-centred field: E = -grad(phi) = (1 - dpsi/dx, -dpsi/dy) by central differences
    Ex = 1 - (np.roll(psi, -1, 0) - np.roll(psi, 1, 0)) / (2 * h)
    Ey = -(np.roll(psi, -1, 1) - np.roll(psi, 1, 1)) / (2 * h)
    E2 = abs(Ex) ** 2 + abs(Ey) ** 2
    shell = (R >= a) & (R < b)
    return float(E2[shell].mean() / abs(Ex.mean()) ** 2)


def periodic_shell_factor_extrap(a, b, eps_ax, eps_my, eps_w, fraction, N=400):
    """Staircase error is first order in h: Richardson 2 F(N) - F(N/2). Converged to ~0.2 %
    (N = 200/400 and 400/600 extrapolations agree to 0.05-0.2 % at p = 0.05 and 0.6)."""
    return 2 * periodic_shell_factor(a, b, eps_ax, eps_my, eps_w, fraction, N) \
        - periodic_shell_factor(a, b, eps_ax, eps_my, eps_w, fraction, N // 2)
