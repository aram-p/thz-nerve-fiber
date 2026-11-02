"""What the tissue itself must change for the two-term waveform fits (Sims 42, 43, 46).

The fits need the unbiased tissue transfer T raised to q: an extra loss (q - 1)|T|^2_dB together
with an extra delay (q - 1)(-arg T). Sim 43 fits the traced tissue phase (numpy FFT, the same
exp(+jwt) sign as the TMM: a delay makes arg T negative) over 0.2-1.4 THz as
    arg T = -2 pi f TAU - OFFSET,  TAU = 0.611 ps, OFFSET = 0.22 rad,
so the delay phase -arg T = 2 pi f TAU + OFFSET. Round 12 flipped the offset's sign; this module
keeps it in one place.
"""
import numpy as np

TAU_S = 0.611e-12
OFFSET_RAD = 0.22
LOSS_DB = {0.6e12: 12.1, 1.0e12: 14.9, 1.4e12: 19.5}   # unbiased |T|^2 below the empty cuvette, Sim 43 (c)


def delay_phase(f_hz, offset=True):
    """-arg T of the unbiased slice (rad) at f_hz."""
    return 2 * np.pi * np.asarray(f_hz, dtype=float) * TAU_S + (OFFSET_RAD if offset else 0.0)


def required_change(f_hz, q, loss_db=None, offset=True):
    """(extra loss dB, extra delay rad) the slice must add for T -> T^q."""
    loss = np.array([LOSS_DB[float(f)] for f in np.atleast_1d(f_hz)]) if loss_db is None else np.asarray(loss_db)
    return (q - 1) * loss, (q - 1) * delay_phase(f_hz, offset)

# Shift s of the tissue pulse with the second pulse free (Sims 43, 47), ps: best fit and the one-sided
# 95 % lower bound from a moving-block bootstrap of the correlated residual (Sim 47). q - 1 ~ s / TAU_S.
SHIFT_PS = {"2 min": (0.043, 0.0275), "4 min": (0.053, 0.0345), "5 min": (0.066, 0.0399)}


def min_required_delay(f_hz, scans=("2 min", "4 min", "5 min")):
    """Smallest extra delay (rad) of the tissue pulse that the traces allow (95 %, one-sided) at f_hz, over the scans."""
    return 2 * np.pi * np.asarray(f_hz, dtype=float) * min(SHIFT_PS[s][1] for s in scans) * 1e-12
