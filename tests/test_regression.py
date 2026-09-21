# -*- coding: utf-8 -*-
"""Regression tests for cooling_stack.

Numerical anchors were verified on 2026-09-21 with the default
reference-engine design point (cooling_params.py). Tolerances are
deliberately tight (+/- 1 K on wall temperatures) so that any formula
or parameter drift fails loudly instead of silently shifting results.

Run:  pytest tests/ -v
"""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

from cooling_params import (build_bc, scheme_A0, scheme_A1, scheme_B,
                            scheme_C, L_TOTAL)
from cooling_modules import BacksideConvection, FilmRow
from cooling_solver import solve_scheme, diagnose

BC = build_bc()

# T_metal anchors [K], single-node lumped basis
ANCHORS = {'A0': 1158.7, 'A1': 760.9, 'B': 662.1, 'C': 544.5}
RESULTS = {
    'A0': solve_scheme(scheme_A0(), BC, verbose=False),
    'A1': solve_scheme(scheme_A1(), BC, verbose=False),
    'B':  solve_scheme(scheme_B(),  BC, verbose=False),
    'C':  solve_scheme(scheme_C(),  BC, verbose=False),
}


@pytest.mark.parametrize('name', list(ANCHORS.keys()))
def test_wall_temperature_anchor(name):
    """T_metal must match the recorded anchor within +/- 1 K."""
    assert abs(RESULTS[name]['T_metal'] - ANCHORS[name]) <= 1.0


@pytest.mark.parametrize('name', list(ANCHORS.keys()))
def test_energy_residual(name):
    """Energy balance residual must stay below 2%."""
    assert RESULTS[name]['energy_residual_pct'] < 2.0


def test_A1_pressure_drop():
    """Scheme A1 dP is dominated by the film orifice: 0.31 kPa."""
    assert abs(RESULTS['A1']['dP_total'] / 1e3 - 0.31) < 0.02


def test_B_saturation_advice():
    """Scheme B is coolant-capacity saturated; the advice layer must
    say so (coolant rise >= 0.9 x (T_metal - T2) and NTU >> 1)."""
    res = RESULTS['B']
    assert res['T_c_rise'] >= 0.9 * (res['T_metal'] - BC['T2'])
    advice = diagnose(res, BC, verbose=False)
    assert any('saturation' in a for a in advice)


def test_A0_low_margin_advice():
    """A0 has a 39 K margin to the 1198 K creep allowable -> the
    advice layer must warn and (no film active) suggest a film row."""
    res = RESULTS['A0']
    assert res['eta_overall'] == 0.0
    assert 0 < (1198.0 - res['T_metal']) < 100
    advice = diagnose(res, BC, verbose=False)
    assert any('margin' in a for a in advice)
    assert any('FilmRow' in a for a in advice)


def test_switch_disabled_modules():
    """span=0 modules must be filtered out and reported as disabled."""
    res = RESULTS['A0']      # film row declared with span=0
    assert 'film row 1' in res['disabled']
    assert res['eta_overall'] == 0.0


def test_all_cold_disabled_raises():
    """A scheme with every cold-side module off must fail loudly."""
    bad = [FilmRow('film row 1', 0.0244, 0.030, d=2e-3, n_rows=12,
                   n_per_row=20)]
    with pytest.raises(ValueError):
        solve_scheme(bad, BC, verbose=False)


def test_zone_partition_warning():
    """Overlapping spans must be detected: >5% partition error."""
    overlap = [BacksideConvection('smooth backside', 0.0, L_TOTAL),
               BacksideConvection('smooth backside (aft)', 0.04, L_TOTAL)]
    res = solve_scheme(overlap, BC, verbose=False)
    assert res['partition_err_pct'] > 5.0
