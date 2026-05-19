# Python script as source of truth for the COMSOL model

The COMSOL model (geometry, materials, mesh, studies) is defined in Python code using the `mph` library and rebuilt from scratch on every run. The `.mph` file is a disposable build artifact, gitignored. We chose this over keeping a hand-built `.mph` as the source of truth so that geometry and physics variations can be parameterised and version-controlled like ordinary code, and so the entire thesis pipeline — model build, sweep, export, plot — lives in one diffable language.

## Considered Options

- **`.mph` as source of truth, scripts only override parameters.** Rejected: locks us out of geometry sweeps (varying node length, fiber radius, number of nodes) without manually re-editing the `.mph` in the GUI, which is exactly what this project is trying to avoid. Also makes the model un-diffable in git.
- **Java API directly (no Python wrapper).** Rejected: the `mph` Python package exposes the same Java API with much less ceremony, and Python is already the language of the analysis and plotting layers — keeping a single language end-to-end is worth more than the (small) loss of fidelity to COMSOL's native Java surface.

## Consequences

- Onboarding cost: porting the existing physics model into a Python builder is ~one week of upfront work before any sweep can run.
- Debugging: when a Python-built model misbehaves, we may occasionally need to open the `.mph` build artifact in the COMSOL GUI to inspect geometry/mesh visually. RDP into the Windows box is the escape hatch.
