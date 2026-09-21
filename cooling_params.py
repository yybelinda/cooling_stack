# -*- coding: utf-8 -*-
"""
=================================================================
 cooling_params.py — general input (the single file to edit)
=================================================================
Position:
  - Standalone use: edit this file, then run `python cooling_main.py`
  - Library use inside a larger suite:
        from cooling_stack.cooling_solver import solve_scheme
        from cooling_stack.cooling_params import build_bc, scheme_B
        res = solve_scheme(scheme_B(), build_bc())
  - Default values = reference-engine cruise design point (the same
    published 1D design case, see README "Provenance").

Ribbed-channel policy:
  Annulus ribs are rarely machined at this engine scale, so the default
  schemes keep the ribbed channel disabled (span=0). Scheme A1 is a
  correlation cross-check configuration that enables it deliberately.
=================================================================
"""

# =========================================================
# 1. Global boundary conditions (reference engine design point,
#    hot-station basis)
# =========================================================
P2        = 3.338e5      # [Pa]  compressor exit total pressure (3.338 bar)
T2        = 494.1        # [K]   compressor exit total temp (= coolant supply)
T_GAS     = 1334.0       # [K]   hot-station gas temperature (suite-consistent)
M_GAS     = 0.13422      # [kg/s] hot-station gas mass flow
M_COOL    = 0.01796      # [kg/s] total cooling-air flow (=12% of air, 17.96 g/s)
D_LINER   = 0.0637       # [m]   liner diameter (63.7 mm)
D_CASING  = 0.0850       # [m]   casing inner diameter (85.0 mm)
L_TOTAL   = 0.172        # [m]   liner axial length (172 mm)
MATERIAL  = 'Hastelloy X'
Z_EVAL    = 0.0344       # [m]   wall-temperature evaluation station (s=34.4 mm)
DP_BUDGET = 0.04 * P2    # [Pa]  available cooling-air pressure budget (4%)

# ── Module master switch ─────────────────────────────────
# A cooling module with span=0 (or span < MIN_SPAN) does NOT exist:
# that axial region is not covered and the module is switched off.
# The solver filters cold-side modules and prints a [SWITCH] note;
# film rows self-check their span and report eta=0 when disabled.
#   e.g. FilmRow('film row 1', 0.0244)        -> span=0  -> off
#        FilmRow('film row 1', 0.0244, 0.030) -> span=5.6mm > MIN_SPAN -> on
MIN_SPAN  = 1.0e-3       # [m] modules shorter than 1 mm count as disabled

# Wall layers (hot face -> cold face order, pure thermal resistances):
#   (name, thickness [m], type, parameter)  type: 'tbc'=k given / 'metal'=lookup
WALL_LAYERS = [
    ('TBC (YSZ)',    0.2e-3,  'tbc',   1.5),      # k=1.5 W/m/K
    ('Metal (HX)',   1.5e-3,  'metal', MATERIAL),
]

# Axial zones (station positions [m] from dome) — used by the layout figure
ZONES = [
    ('Dome',    0.000),
    ('Primary', 0.0344),
    ('Dilution', 0.138),
    ('Exit',    0.172),
]


# =========================================================
# 2. Scheme definitions (module parameters are edited here;
#    the module library lives in cooling_modules.py)
#    Rule: list order = coolant flow direction; film rows go at
#    the end of the path; span length 0 = module not used
#    (small engines: annulus ribs rarely machined, thin dP budget
#    rarely affords film — modules stay in the library, set span>0
#    to re-enable).
# =========================================================

def scheme_A0():
    """A0 — baseline liner configuration (smooth backside, no ribs,
    no film). The film row is listed but span=0 keeps it disabled;
    to enable it, give the second axial position:
    FilmRow('film row 1', 0.0244, 0.030)  (span=5.6 mm)."""
    from cooling_modules import BacksideConvection, FilmRow
    return [
        BacksideConvection('smooth backside', 0.0, L_TOTAL),
        FilmRow('film row 1', 0.0244, d=2e-3, n_rows=12, n_per_row=20),
    ]

def scheme_A1():
    """A1 — correlation cross-check configuration: full-length ribbed
    channel + film row ENABLED. This scheme reproduces the anchor of an
    independent multi-station reference calculation (T_metal = 750 K)
    within the stated tolerance; see tests/test_regression.py.
    IMPORTANT: the rib enhancement multiplies the smooth-channel
    coefficient (single-channel basis) — RibbedChannel must REPLACE
    BacksideConvection over the same span. Stacking both over the full
    length double-counts wetted area (measured +8% over-absorption
    in the first release)."""
    from cooling_modules import RibbedChannel, FilmRow
    return [
        RibbedChannel('ribbed channel', 0.0, L_TOTAL,
                      rib_height=0.5e-3, rib_pitch=10e-3, rib_angle=90.0),
        FilmRow('film row 1', 0.0244, 0.030,
                d=2e-3, n_rows=12, n_per_row=20),
    ]

def scheme_B():
    """B — two-layer wall with impingement (no film): dome-region
    impingement + smooth backside over the aft section. Coolant is
    exhausted overboard after the aft section (no film on the gas side)."""
    from cooling_modules import BacksideConvection, ImpingementJet, FilmRow
    return [
        ImpingementJet('impingement (dome)', 0.0, 0.040, dP_j=5e3, n_j=15,
                       H=12.8e-3),
        BacksideConvection('smooth backside (aft)', 0.040, L_TOTAL),
        FilmRow('film row 1', 0.0244, d=2e-3, n_rows=12, n_per_row=20),
    ]

def scheme_C():
    """C — impingement + two staggered film rows: row 1 (z=14.4) +
    row 2 (z=24.4). Demonstrates multi-row film superposition — with the
    'no-film' baseline accepted, this scheme is the B + dual-film
    sensitivity case: how much wall temperature is bought back and at
    what dP cost."""
    from cooling_modules import BacksideConvection, ImpingementJet, FilmRow
    return [
        ImpingementJet('impingement (dome)', 0.0, 0.040, dP_j=5e3, n_j=15,
                       H=12.8e-3),
        BacksideConvection('smooth backside (aft)', 0.040, L_TOTAL),
        FilmRow('film row 1', 0.0144, 0.020, d=2e-3, n_rows=12, n_per_row=20),
        FilmRow('film row 2', 0.0244, 0.030, d=2e-3, n_rows=12, n_per_row=20),
    ]

SCHEMES = [
    ('A0_baseline_smooth_nofilm', scheme_A0),
    ('A1_ribbed_film_crosscheck', scheme_A1),
    ('B_impingement_nofilm',      scheme_B),
    ('C_impingement_dual_film',   scheme_C),
]


# =========================================================
# 3. Assemble the global boundary dict (solver input)
# =========================================================

def build_bc():
    import numpy as np
    A_annulus = np.pi * (D_CASING**2 - D_LINER**2) / 4.0
    return {
        'P2': P2, 'T2': T2, 'T_gas': T_GAS, 'm_gas': M_GAS,
        'm_cool': M_COOL, 'D_liner': D_LINER, 'D_casing': D_CASING,
        'L_total': L_TOTAL, 'A_wall': np.pi * D_LINER * L_TOTAL,
        'A_annulus': A_annulus, 'D_hyd': D_CASING - D_LINER,
        'wall_layers': WALL_LAYERS, 'z_eval': Z_EVAL,
        'material': MATERIAL, 'dp_budget': DP_BUDGET, 'zones': ZONES,
        'min_span': MIN_SPAN,
    }
