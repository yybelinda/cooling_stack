# -*- coding: utf-8 -*-
"""
=================================================================
 gas_props.py — embedded property library (standalone operation)
=================================================================
Position:
  Minimal property set required to run cooling_stack on its own.
  All formulas are identical to the gas-properties module of the
  parent 1D combustor suite (checked formula by formula).

Formulation:
  cp    : truncated NASA polynomial air (valid 250-1500 K)
  mu    : Sutherland (mu0=1.716e-5, T0=273 K, S=111 K)
  k     : k = cp*mu/Pr, Pr = 0.71
  rho   : ideal gas P/(R*T)
  metal k: Haynes International published datasheets
=================================================================
"""

import numpy as np


def air_props(T, P=101325.0):
    """Air properties.

    Input: T [K], P [Pa]
    Output: dict(cp, mu, k, Pr, rho, R)
    """
    T = max(min(T, 1500.0), 200.0)
    cp = (1031.5 - 0.28217 * T + 7.8044e-4 * T**2
          - 6.9777e-7 * T**3 + 2.2737e-10 * T**4)
    R = 287.0
    rho = P / (R * T)
    mu = 1.716e-5 * (273.0 + 111.0) / (T + 111.0) * (T / 273.0)**1.5
    Pr = 0.71
    k = cp * mu / Pr
    return {'cp': cp, 'mu': mu, 'k': k, 'Pr': Pr, 'rho': rho, 'R': R}


# =================================================================
# Metal thermal conductivity (Haynes published datasheets).
# If the parent suite's table is updated, keep the copies in sync.
# =================================================================

MATERIAL_K = {
    'Hastelloy X': {   # GH3536 / UNS N06002
        'T_k':  [20,  100, 300, 600, 700, 800, 900],        # [degC]
        'k':    [9.7, 11.1, 14.7, 20.6, 22.8, 25.0, 27.4],  # [W/m/K]
        'source': 'Haynes International datasheet',
    },
    'Haynes 188': {    # GH5188 / UNS R30188
        'T_k':  [20,  100, 300, 600, 700, 800, 900, 1000],
        'k':    [10.6, 12.4, 16.9, 22.8, 24.9, 27.1, 29.1, 31.2],
        'source': 'Haynes International datasheet',
    },
}


def metal_k(material, T_K):
    """Metal conductivity lookup (linear interpolation, degC domain).
    Clamped to table ends outside the data range."""
    mat = MATERIAL_K[material]
    T_c = T_K - 273.15
    return float(np.interp(T_c, mat['T_k'], mat['k']))


# Creep allowables [K] (30 MPa / 1000 h basis, back-calculated from
# creep-rupture data): Hastelloy X 1198 K / Haynes 188 1225 K.
# Used by diagnose() for the wall-temperature margin check. At a
# different stress level, re-derive from the material's
# creep-rupture curves.
MATERIAL_TALLOW = {
    'Hastelloy X': 1198.0,
    'Haynes 188':  1225.0,
}
