import numpy as np
from thznerve.physics.materials import WATER_KS1996, water_liebe
f = np.array([0.6, 1.0, 1.5, 2.0, 2.5, 3.0]) * 1e12
F = f / 1e12
def L(F):  # Liebe 1991 broadband Lorentzians, exp(+jwt): A/(f0^2 - f^2 + j f g)
    return 25.03 / (5.11**2 - F**2 + 1j * F * 4.46) + 282.4 / (18.2**2 - F**2 + 1j * F * 15.4)
ks = WATER_KS1996(f); li = water_liebe(25)(f)
v1 = li + L(F) - L(0)
for name, e in (("Liebe25 DD", li), ("Liebe25 DD + Lorentz - static", v1)):
    print(name, "re/KS", np.round(e.real / ks.real, 3), "im/KS", np.round(e.imag / ks.imag, 3))
print("KS eps", np.round(ks, 2))
