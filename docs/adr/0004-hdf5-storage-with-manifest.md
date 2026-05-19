# HDF5 per-simulation storage with a central manifest; 2D-slice always, 3D opt-in

Each simulation writes one HDF5 file containing scalar metrics, a 1D axial profile (`|E|(z)`), and a 2D axial slice (`|E|(x, z)` at y=0). Full 3D fields are opt-in via a sweep-config flag and used sparingly for hero figures. Every simulation also appends one row to a top-level `results/manifest.csv` with parameters, output path, timestamp, COMSOL version, git SHA, and YAML config hash. We chose this over CSV-only or always-full-3D because it covers ~90% of expected thesis figures (spectra, axial profiles, cross-section heatmaps) at modest storage cost (~1 MB/sim), while keeping the door open for richer figures when needed.

## Considered Options

- **CSV everywhere.** Rejected: pandas-friendly, but painful for sweeps with 3+ varying dimensions and for 2D slice arrays. HDF5 + xarray is the standard in computational physics for a reason.
- **Always save full 3D fields.** Rejected: ~50 MB/sim × hundreds of sims = bloat that mostly never gets plotted. Better to default to slice + 1D and opt in to 3D for the handful of "hero" simulations.
- **Scalars only, regenerate fields on demand.** Rejected: COMSOL solves are minutes to hours; we don't want to re-solve to make a plot, and disk is much cheaper than compute.
- **SQLite as the manifest.** Rejected for now: CSV is grep-able, diff-able, openable in any spreadsheet. Move to SQLite if the manifest exceeds a few thousand rows or query patterns get complex.

## Consequences

- Every figure script reads from `manifest.csv` (filter by params) and loads the matching HDF5 — a uniform contract.
- HDF5 files carry their own provenance (COMSOL version, git SHA, config hash) in metadata, so a single HDF5 is self-describing even if separated from the manifest.
- "3D opt-in" decisions are made at sweep-config time, not retroactively — if you later want 3D for an old sweep, you re-run that sweep with the flag.
