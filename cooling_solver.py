# -*- coding: utf-8 -*-
"""
=================================================================
 cooling_solver.py — combinator and solver
=================================================================
Responsibilities (contains NO cooling formula — all live in
cooling_modules):
  1. Serializes the scheme (module list) along the coolant path
  2. Solves the energy-balance equation for the single unknown
     Tw_cold (cold-side wall temperature) by bisection
  3. Outputs: T_ad (adiabatic wall temperature) / layer-by-layer wall
     temperatures / per-module h*A*q*dP / energy residual

Solver format (why not a fixed-point averaging scheme):
  Fixed-point schemes with frozen coolant temperature assume a
  LINEAR cold side. With serial coolant heating the cold-side demand
  gain w.r.t. Tw_cold can exceed 1 at high NTU (impingement
  NTU = h*A/(m*cp) ~ 3.8), the fixed-point map stops contracting and
  diverges (first release: schemes B/C ran away to 18,000 K with
  >1500% energy residual). This solver instead bisects
  F(Tw_cold) = cold-side demand - hot-side supply: demand grows with
  Tw, supply decreases with Tw, F necessarily crosses zero — robust
  for any module combination. Scheme A: both methods agree within
  2 K (the averaging scheme carries an intrinsic ~3% energy residual).

Physics conventions (honest accounting):
  - Single-node design-point lumped model (evaluation station z_eval
    defaults to the hot station); station-by-station axial extension
    is a listed limitation (see README)
  - Coolant heats up along the path (a frozen T_c = T2 convention
    ignores this; this model includes it — more physical; cross-
    checks against the frozen-Tc reference run ~10-14 K higher wall
    temperature, a convention difference, not an error)
  - Second-law clamp: any cold-side module outlet temperature <=
    current wall temperature (heat-exchanger effectiveness <= 1)
  - Multi-row superposition: eta_overall = 1 - prod(1 - eta_i)
=================================================================
"""

import numpy as np
from gas_props import air_props, metal_k, MATERIAL_TALLOW


# =================================================================
# Automatic diagnosis and advice layer
# Design decision: the tool reports suggestions based on results —
# it never modifies parameters automatically. Parameter changes
# remain a human decision (edit cooling_params.py yourself).
# All thresholds are explicit below; changing a threshold = editing
# this function.
# =================================================================

def diagnose(res, bc, verbose=True):
    """Rule-based diagnosis of a solve result; returns the advice list
    and stores it in res['advice'].

    Checks and thresholds:
      1. Wall-temperature margin: T_metal vs creep allowable
         (30 MPa / 1000 h, gas_props) — <0 over allowable / <100 K warning
      2. Cold-side saturation: coolant temperature rise >= 0.9*(T_metal-T2)
         and (T_metal-T2) > 50 K -> lengthening/strengthening cold-side
         modules is ineffective; bottleneck = coolant heat capacity
      3. dP budget: >100% over budget / >80% near limit
      4. No film and margin <150 K -> suggest trying a FilmRow
         (see README Example 2 reference values)
      5. Film supply temperature > 0.8*T_gas -> film approaching
         ineffectiveness; cut upstream heat absorption
    """
    adv = []
    T_metal, T2 = res['T_metal'], bc['T2']
    T_allow = MATERIAL_TALLOW.get(bc.get('material', ''), 1200.0)
    margin = T_allow - T_metal
    if margin < 0:
        adv.append(f"Wall temperature OVER allowable: T_metal={T_metal:.0f}K > "
                   f"creep allowable {T_allow:.0f}K (30MPa/1000h) — "
                   "cooling must be strengthened or T_gas reduced")
    elif margin < 100:
        adv.append(f"Wall temperature margin only {margin:.0f}K (<100K warning "
                   "line) — recheck the actual stress level (allowable is "
                   "for 30MPa/1000h; re-derive at your stress)")
    if (T_metal - T2) > 50 and res['T_c_rise'] >= 0.9 * (T_metal - T2):
        adv.append("Cold-side saturation (coolant temperature rise ~= "
                   "T_metal - T2): lengthening/strengthening cold-side modules "
                   "is ineffective; bottleneck = coolant heat capacity. "
                   "Effective levers: more coolant / a film row / lower T_gas")
    ratio = res['dP_total'] / bc['dp_budget']
    if ratio > 1.0:
        adv.append(f"dP over budget ({ratio*100:.0f}%) — increase hole "
                   "diameter / reduce hole count / lower impingement dP_j")
    elif ratio > 0.8:
        adv.append(f"dP has consumed {ratio*100:.0f}% of the budget, "
                   "approaching the limit")
    if res['eta_overall'] == 0 and 0 <= margin < 150:
        adv.append("No film active (eta=0) and wall margin <150K: if the dP "
                   "budget allows, try enabling a FilmRow (set span>0) — "
                   "README Example 2 reference: -390K for 0.31kPa")
    if res['eta_overall'] > 0 and res['T_c_film'] > 0.8 * bc['T_gas']:
        adv.append(f"Film supply temperature {res['T_c_film']:.0f}K exceeds "
                   "0.8*T_gas — too much upstream heat absorption, film "
                   "approaching ineffectiveness; add coolant flow or shorten "
                   "the impingement/backside absorption length")
    res['advice'] = adv
    if verbose:
        for a in adv:
            print(f"    [ADVICE] {a}")
        if not adv:
            print("    [ADVICE] All indicators within thresholds; no changes suggested")
    return adv


def _wall_resistance(wall_layers, bc):
    """Wall thermal resistances in series: [(name, thickness m, type, param)]
    -> R_total [m2K/W]. Types: 'tbc' (k given) / 'metal' (k from material
    table, temperature-dependent -> updated by the outer loop)."""
    layers, R = [], 0.0
    for (nm, t, typ, prm) in wall_layers:
        if typ == 'tbc':
            R += t / prm
            layers.append({'name': nm, 't': t, 'R': t / prm, 'k': prm})
        else:  # metal
            k_m = metal_k(prm, 750.0)   # initial guess, updated with wall T
            R += t / k_m
            layers.append({'name': nm, 't': t, 'R': t / k_m, 'k': k_m,
                           'material': prm})
    return R, layers


def solve_scheme(scheme, bc, verbose=False):
    """Solve one cooling scheme (single-node design point).

    Parameters:
      scheme : [CoolingModule...] ordered along coolant flow (cold/film
               mixed; film rows conventionally at the path end; cold
               order = flow direction)
      bc     : dict of global boundaries (see cooling_params.build_bc)
    Returns:
      dict: T_ad, T_surf_hot, T_surf_cold, T_metal, q_flux, eta_overall,
            T_c_film, path (per-module results), film_res, dP_total,
            energy residual, wall_layers, convergence info
    """
    # --- MIN_SPAN master switch ---
    # Modules with span=0 or span<MIN_SPAN do not exist: switched off.
    # Cold side filtered here (film rows self-check span).
    min_span = bc.get('min_span', 1e-3)
    cold_all = [m for m in scheme if m.kind == 'cold']
    films = [m for m in scheme if m.kind == 'film']
    cold = [m for m in cold_all if m.span >= min_span]
    disabled = [m.name for m in cold_all if m.span < min_span] + \
               [f.name for f in films if f.span < min_span]
    if not cold:
        raise ValueError(
            'No enabled cold-side module in the scheme (span >= '
            f'{min_span*1e3:.1f}mm) — wall heat has nowhere to go, '
            'physically invalid. Enable at least one '
            'BacksideConvection / impingement / ribbed channel.')
    if disabled:
        print(f"  [SWITCH] Modules with span<{min_span*1e3:.1f}mm treated "
              f"as disabled: {', '.join(disabled)}")

    # --- Derived boundary quantities ---
    T2, P2 = bc['T2'], bc['P2']
    T_gas, m_gas = bc['T_gas'], bc['m_gas']
    m_cool = bc['m_cool']
    A_wall = bc['A_wall']
    A_ref = np.pi * bc['D_liner']**2 / 4.0
    prop_g = air_props(T_gas, bc['P2'] * 0.99)   # gas ~ air properties (same basis)
    V_g = m_gas / (prop_g['rho'] * A_ref)
    bc = dict(bc, A_ref=A_ref, rho_g=prop_g['rho'], V_g=V_g)

    R_wall, wall_layers = _wall_resistance(bc['wall_layers'], bc)

    # Hot-side convection (same convention as the parent suite:
    # Tw=0.6Tg -> temperature-ratio exponent = 1)
    Re_g = prop_g['rho'] * V_g * bc['D_liner'] / prop_g['mu']
    Nu_g = 0.021 * Re_g**0.8 * prop_g['Pr']**0.43
    h_g = Nu_g * prop_g['k'] / bc['D_liner']

    def _path_walk(Tw_cold):
        """Walk the coolant along the serial path -> (film supply
        temperature, per-module results, total heat absorbed)"""
        T_c, path, q_sum = T2, [], 0.0
        for m in cold:
            st = {'T_c': T_c, 'P_c': P2, 'm_c': m_cool,
                  'T_wall_cold': Tw_cold}
            res = m.solve(st, bc)
            if res['T_out'] > Tw_cold:      # second-law clamp: outlet <= wall T
                res['T_out'] = Tw_cold
                res['q'] = m_cool * res['cp'] * (Tw_cold - T_c)
                res['clamped'] = True
            path.append((m, st, res))
            q_sum += res['q']
            T_c = res['T_out']
        return T_c, path, q_sum

    def _films(T_c_film):
        """Film rows: multi-row eta superposition -> (eta_overall, results)"""
        etas, film_res = [], {}
        for f in films:
            st = {'T_c': T_c_film, 'P_c': P2, 'm_film': m_cool}
            fr = f.solve(st, bc)
            film_res[f.name] = fr
            if fr.get('active'):
                etas.append(fr['eta'])
        eta_ov = 1.0 - np.prod([1.0 - e for e in etas]) if etas else 0.0
        return eta_ov, film_res

    def _hot_flux(T_ad, Tw_cold, R_wall):
        """Hot-side supply flux [W/m2]: solve q=h_g(T_ad-Tw_hot)+rad(Tw_hot)
        =(Tw_hot-Tw_cold)/R_wall for Tw_hot (fixed point, small R converges
        fast)"""
        Tw_hot = Tw_cold + 5.0
        q = 0.0
        for _ in range(50):
            K = 0.5 * (0.1 + 0.9 * T_gas / 1000.0)
            eps_g = 1.0 - np.exp(-K * (bc['P2'] / 101325.0)
                                 * (0.6 * bc['D_liner']))
            q_rad = (5.67e-8 * eps_g * (T_gas**4 - Tw_hot**4)
                     / (1.0 / eps_g + 1.0 / 0.8 - 1.0))
            q = h_g * (T_ad - Tw_hot) + q_rad
            Tw_hot_new = Tw_cold + q * R_wall
            if abs(Tw_hot_new - Tw_hot) < 0.01:
                Tw_hot = Tw_hot_new
                break
            Tw_hot = 0.5 * (Tw_hot + Tw_hot_new)   # relaxation 0.5
        return q, Tw_hot, q_rad

    def _F(Tw_cold):
        """Energy balance function: F = cold-side demand - hot-side supply.
        Demand grows with Tw, supply decreases with Tw -> monotone zero
        crossing, bisection robust. Returns (F, all intermediates)."""
        T_c_film, path, q_demand = _path_walk(Tw_cold)
        eta_ov, film_res = _films(T_c_film)
        T_ad = T_gas - eta_ov * (T_gas - T_c_film)
        q_hot, Tw_hot, q_rad = _hot_flux(T_ad, Tw_cold, R_wall)
        return q_demand - q_hot * A_wall, dict(
            T_c_film=T_c_film, path=path, q_demand=q_demand,
            eta_ov=eta_ov, film_res=film_res, T_ad=T_ad,
            q_hot=q_hot, Tw_hot=Tw_hot, q_rad=q_rad)

    # =========================================================
    # Outer loop: update metal k with wall temperature (3 rounds suffice,
    # k varies slowly)
    # Inner loop: bisection on Tw_cold
    # =========================================================
    lo, hi = T2 - 5.0, T_gas            # F(lo)<0 (demand~0), F(hi)>0 (demand max)
    art = None
    n_iter = 0
    for _ in range(3):
        lo, hi = T2 - 5.0, T_gas
        for it in range(60):
            mid = 0.5 * (lo + hi)
            F, art = _F(mid)
            n_iter += 1
            if F > 0:
                hi = mid
            else:
                lo = mid
            if hi - lo < 0.05:
                break
        # refresh metal thermal resistance with this round's wall temperature
        T_metal_mid = 0.5 * (art['Tw_hot'] + mid)
        for lay in wall_layers:
            if 'material' in lay:
                k_new = metal_k(lay['material'], T_metal_mid)
                R_wall = R_wall - lay['R'] + lay['t'] / k_new
                lay['R'] = lay['t'] / k_new
                lay['k'] = k_new

    Tw_cold = 0.5 * (lo + hi)
    F_final, art = _F(Tw_cold)
    T_ad = art['T_ad']
    Tw_hot = art['Tw_hot']
    q_in = art['q_hot']
    q_rad = art['q_rad']
    T_c_film = art['T_c_film']
    eta_overall = art['eta_ov']
    path, film_res = art['path'], art['film_res']
    q_demand = art['q_demand']

    # ---- Summary ----
    hA_sum = sum(r['h'] * r['A_eff'] for _, _, r in path)
    h_c_eff = hA_sum / A_wall
    # Area-partition check: cold-side module spans should be mutually
    # exclusive and jointly cover the full wall. Overlap -> double-counted
    # heat absorption (falsely low wall T); gap -> uncooled wall (falsely
    # high). Hard warning beyond 5% (never silent).
    A_cold = sum(r['A_eff'] for _, _, r in path)
    partition_err = (A_cold - A_wall) / A_wall * 100.0
    if abs(partition_err) > 5.0:
        print(f"  [WARNING] Cold-side axial coverage deviation {partition_err:+.1f}% "
              f"(sum A={A_cold*1e3:.1f} vs A_wall={A_wall*1e3:.1f} m2*1e-3) — "
              "check module spans for overlap or gaps")
    dP_total = sum(r['dP'] for _, _, r in path) + \
               sum(fr.get('dP', 0.0) for fr in film_res.values())
    T_metal_hot = Tw_hot - q_in * R_wall
    T_metal = 0.5 * (T_metal_hot + Tw_cold)
    q_supply = q_in * A_wall
    e_res = abs(q_demand - q_supply) / max(q_supply, 1.0) * 100.0

    if verbose:
        print(f"  Converged: bisection, {n_iter} evals  T_ad={T_ad:.1f}K  "
              f"Tw_hot={Tw_hot:.1f}K  Tw_cold={Tw_cold:.1f}K  "
              f"T_metal={T_metal:.1f}K")
        print(f"  eta_overall={eta_overall:.3f}  T_c_film={T_c_film:.1f}K "
              f"(rise {T_c_film - T2:.1f}K)  q={q_in/1e3:.2f}kW/m2  "
              f"dP_total={dP_total/1e3:.2f}kPa")
        for m, st, r in path:
            print(f"    [cold] {m.name:22s} h={r['h']:8.1f}  "
                  f"A={r['A_eff']*1e3:6.1f}m2*1e-3  q={r['q']:7.1f}W  "
                  f"T_out={r['T_out']:6.1f}K  dP={r['dP']/1e3:.2f}kPa"
                  + ('  [outlet clamped to wall T: effectiveness->1]' if r.get('clamped') else '')
                  + (f"  ({r.get('regime','')})" if r.get('regime') else ''))
        for nm, fr in film_res.items():
            if fr.get('active'):
                print(f"    [film] {nm:22s} x/s={fr['x_s']:.1f}  "
                      f"M_annulus={fr['M_annulus']:.3f}  M_jet={fr['M_jet']:.3f}  "
                      f"eta={fr['eta']:.3f}  dP={fr['dP']/1e3:.2f}kPa")
            else:
                print(f"    [film] {nm:22s} {fr.get('note','')}")
        print(f"    Energy check: cold-side absorption sum q={q_demand:.1f}W vs "
              f"hot-side supply q*A_wall={q_supply:.1f}W (residual {e_res:.1f}%)")

    return {
        'T_ad': T_ad, 'T_surf_hot': Tw_hot, 'T_surf_cold': Tw_cold,
        'T_metal': T_metal, 'T_metal_hot': T_metal_hot,
        'q_flux': q_in, 'eta_overall': eta_overall,
        'T_c_film': T_c_film, 'T_c_rise': T_c_film - T2,
        'path': path, 'film_res': film_res, 'h_g': h_g,
        'h_c_eff': h_c_eff, 'R_wall': R_wall, 'wall_layers': wall_layers,
        'dP_total': dP_total, 'n_iter': n_iter,
        'q_rad': q_rad, 'Re_g': Re_g, 'energy_residual_pct': e_res,
        'converged': (hi - lo) < 0.05, 'disabled': disabled,
        'partition_err_pct': partition_err,
    }
