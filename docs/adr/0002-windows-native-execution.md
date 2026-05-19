# Everything runs on the Windows machine

Code, simulations, data, and plots all live on the Windows box where COMSOL is licensed and installed. The Mac is a thin editing client that SSHs into Windows; Claude Code runs natively on Windows. We chose this over Mac-driven workflows (COMSOL Server, network shares, cloud-synced results) because the simulation, analysis, and figure-generation toolchain stays co-located on a single machine with no client/server marshalling, license server complications, or sync conflicts.

## Considered Options

- **Mac-driven via COMSOL Server.** Rejected: requires a floating network license, the `mph` package's server-mode story is less mature than its standalone-mode story, and LAN latency makes interactive debugging painful.
- **Mac for code + cloud sync (Dropbox/OneDrive) for results.** Rejected: adds a sync layer that can silently fail or conflict, and Jupyter/plotting still wants to be near the data.
- **Mac for code + git/git-LFS for results.** Rejected: a thesis-scale sweep produces hundreds of MB to multi-GB of HDF5 per run; committing that to git history is wrong shape, and git-LFS adds cost/quota for public repos.

## Consequences

- The Mac cannot make progress on simulations or analysis when the Windows box is offline; all serious work routes through Windows.
- Interactive plotting (Jupyter, matplotlib windows) happens either over an SSH-tunnelled Jupyter port or over RDP — neither is as snappy as native Mac.
- Backups become a Windows-side concern (the manifest + HDF5 files are the canonical thesis data, not anything on the Mac).
