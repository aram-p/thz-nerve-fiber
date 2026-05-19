"""Pydantic schema for YAML sweep configs.

Wire format (per ADR-0003):

    name: freq_sweep_001
    description: ...
    model:
      axon_radius_um: 5
      myelin_radius_um: 7
      node_length_um: 40
      internode_length_um: 100
      external_half_width_um: 20
      node_sigma_S_per_m: 0
    sweep:
      frequency_THz:
        linspace: {start: 0.1, stop: 2.0, n: 26}
      # optional zipped group:
      # zip:
      #   frequency_THz: {list: [0.6, 1.2]}
      #   node_sigma_S_per_m: {list: [0, 0.5]}
    output:
      dir: results/freq_sweep_001/
      save_3d_field: false
"""

from typing import Union

from pydantic import BaseModel, ConfigDict, Field, model_validator


class Linspace(BaseModel):
    start: float
    stop: float
    n: int = Field(ge=1)


class Logspace(BaseModel):
    """Numpy convention: start/stop are EXPONENTS, output is base ** linspace(start, stop, n)."""

    start: float
    stop: float
    n: int = Field(ge=1)
    base: float = 10.0


class SweepExpr(BaseModel):
    """A single-parameter sweep expression. Exactly one of linspace / logspace / list."""

    model_config = ConfigDict(extra="forbid")

    linspace: Linspace | None = None
    logspace: Logspace | None = None
    list_: list[float] | None = Field(default=None, alias="list")

    @model_validator(mode="after")
    def _exactly_one(self) -> "SweepExpr":
        provided = sum(x is not None for x in (self.linspace, self.logspace, self.list_))
        if provided != 1:
            raise ValueError("SweepExpr must have exactly one of: linspace, logspace, list")
        return self


class ZipGroup(BaseModel):
    """A `zip:` block — its inner parameters vary together (same length)."""

    model_config = ConfigDict(extra="forbid")

    zip: dict[str, SweepExpr]


SweepValue = Union[SweepExpr, ZipGroup]


class ModelParams(BaseModel):
    """Geometry + per-node-conductivity baseline (fixed for the sweep unless swept)."""

    model_config = ConfigDict(extra="forbid", protected_namespaces=())

    axon_radius_um: float
    myelin_radius_um: float
    node_length_um: float
    internode_length_um: float
    external_half_width_um: float
    node_sigma_S_per_m: float


class OutputSpec(BaseModel):
    model_config = ConfigDict(extra="forbid")

    dir: str
    save_3d_field: bool = False


class SweepConfig(BaseModel):
    """Top-level YAML sweep config."""

    model_config = ConfigDict(extra="forbid", protected_namespaces=())

    name: str
    description: str = ""
    model: ModelParams
    sweep: dict[str, SweepValue]
    output: OutputSpec
