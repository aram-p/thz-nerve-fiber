"""Analytic physics: materials, slab TMM, multilayer-cylinder Mie series, closed-form node model.

Convention everywhere: time dependence exp(+j w t) (COMSOL / engineering),
so a lossy relative permittivity is eps = eps' - j eps'' with eps'' > 0, and a
wave travelling in +x is exp(-j k x) with Im(k) < 0.
"""
