# Sweeps are defined as YAML configs, not Python scripts

Each experiment is described by a YAML file in `experiments/` declaring what parameters to vary, what to hold fixed, and where to write results. A single `run.py` dispatcher loads any YAML, validates it (Pydantic), expands the parameter grid, and drives the sweep. We chose this over imperative Python sweep scripts so that experiments are flat, readable, comparable side-by-side, and easy to discuss with a thesis advisor without reading Python.

## Considered Options

- **Imperative Python — one `.py` per sweep.** Rejected despite some real advantages (natural NumPy expressions for non-uniform grids, no schema layer to maintain). The readability and discussability of YAML for advisor reviews and for future-me reading old experiments outweighs the loss of expressive power; a small sweep expression vocabulary (`linspace`, `logspace`, `list`, nested keys) covers what NumPy gave us.

## Consequences

- We need to maintain a YAML schema (Pydantic models) and a small expression vocabulary; this is a non-trivial up-front investment but stabilises quickly.
- Sweeps that genuinely need computed parameters (e.g., a frequency grid whose spacing depends on another variable) become awkward. The escape hatch is to add a new expression type to the schema rather than dropping back to Python.
- The YAML file is the experimental record alongside the git SHA — what was *intended* to run, separable from what was actually computed (in the HDF5 metadata).
