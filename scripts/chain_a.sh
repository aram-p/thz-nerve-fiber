#!/bin/sh
# Wait for Sim 31 L2 perp, stop that run before L3 (would exceed RAM), then run the queue.
cd /c/Users/aramp/dev/thz-nerve-fiber
until grep -q "conv_L2 perp" results/sim30/sim31_mesh.log; do sleep 5; done
pid=$(powershell -NoProfile -Command "(Get-CimInstance Win32_Process -Filter \"Name='python.exe'\" | Where-Object { \$_.CommandLine -like '*sim31_convergence.py mesh*' }).ProcessId")
for p in $pid; do taskkill //F //PID $p >/dev/null 2>&1; done
sleep 3
uv run python -u experiments/prod/sim31_convergence.py short > results/sim30/sim31_short.log 2>&1
uv run python -u experiments/prod/sim31_convergence.py bc > results/sim30/sim31_bc.log 2>&1
uv run python -u experiments/prod/sim32_scaling.py all > results/sim30/sim32.log 2>&1
echo CHAIN_A_DONE
