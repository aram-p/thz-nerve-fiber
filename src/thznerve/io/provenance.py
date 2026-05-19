"""Provenance helpers — record where a simulation came from.

`get_comsol_version()` accepts an existing `mph` client to avoid the
~25s cold-start cost; pass `client=mph.start()` if you already have one.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path
from typing import Any


def get_comsol_version(client: Any | None = None) -> str:
    """COMSOL version string. Starts a new mph client if none supplied."""

    if client is None:
        import mph

        client = mph.start()
    return str(client.version)


def get_git_sha(repo_root: Path | str | None = None) -> str:
    """Current HEAD SHA of the repo. Defaults to the package's own repo."""

    if repo_root is None:
        # repo_root = <pkg>/src/thznerve/io/provenance.py -> repo root is 4 parents up
        repo_root = Path(__file__).resolve().parents[3]
    result = subprocess.run(
        ["git", "-C", str(repo_root), "rev-parse", "HEAD"],
        capture_output=True,
        text=True,
        check=True,
    )
    return result.stdout.strip()


def config_hash(yaml_dict: dict[str, Any]) -> str:
    """Stable 12-char hash of a sweep config dict. Canonical JSON, sha256, hex prefix."""

    canonical = json.dumps(yaml_dict, sort_keys=True, separators=(",", ":"), default=str)
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()[:12]
