# -*- coding: utf-8 -*-
"""
=================================================================
 cooling_main.py — main entry (read schemes -> solve -> figures)
=================================================================
Three ways to use:
  1. Standalone run:  python cooling_main.py
     (runs the default reference-engine schemes from cooling_params.py,
      prints all results and generates 3 kinds of figures)
  2. Library use:
         from cooling_solver import solve_scheme
         from cooling_params import build_bc, scheme_B
         res = solve_scheme(scheme_B(), build_bc())
  3. Custom scheme: add a scheme_X() to cooling_params.py (module list
     ordered along coolant flow) and register it in SCHEMES
     (composition rules at the top of cooling_modules.py)

Figures (all generated from actual computation):
  fig1_layout_*.png   axial layout: colored blocks above/below the wall
                      = axial coverage of each cooling method; film rows
                      span both sides of the wall
  fig2_coolant_path_*.png  along the coolant path: per-module temperature
                      steps + dP accumulation
  fig3_scheme_comparison.png  per-scheme T_metal / T_ad / dP / coolant
                      temperature rise, four indicators
=================================================================
"""

import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# Make this folder importable regardless of the working directory;
# figures are saved to the folder containing this script (= current
# folder). If that folder is not writable, fall back to the current
# working directory and print the location clearly.
_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

from cooling_solver import solve_scheme, diagnose
from cooling_params import build_bc, SCHEMES, ZONES, L_TOTAL, D_LINER

OUTDIR = _HERE if os.access(_HERE, os.W_OK) else os.getcwd()

plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['axes.unicode_minus'] = False

# Module colors (shared by layout and path figures)
COLOR = {
    'smooth backside':      '#93c5fd',
    'smooth backside (aft)': '#93c5fd',
    'ribbed channel':       '#0d9488',
    'impingement (dome)':   '#3b82f6',
    'film row 1':           '#ef4444',
    'film row 2':           '#f97316',
    'TBC (YSZ)':            '#a855f7',
    'Metal (HX)':           '#9ca3af',
}


# =========================================================
# Main flow: solve scheme by scheme
# =========================================================

def run_all():
    bc = build_bc()
    results = {}
    print('=' * 68)
    print('  cooling_stack V1 — serial-path cooling scheme solver '
          '(reference-engine design point defaults)')
    print('=' * 68)
    print(f"Boundary: P2={bc['P2']/1e3:.1f}kPa  T2={bc['T2']:.1f}K  "
          f"T_gas={bc['T_gas']:.0f}K  m_gas={bc['m_gas']:.4f}kg/s  "
          f"m_cool={bc['m_cool']*1e3:.2f}g/s")
    print(f"Geometry: D_liner={bc['D_liner']*1e3:.1f}mm  "
          f"L={bc['L_total']*1e3:.0f}mm"
          f"  wall: TBC 0.2mm + {bc['material']} 1.5mm  "
          f"dP budget={bc['dp_budget']/1e3:.1f}kPa")
    for name, fn in SCHEMES:
        print(f"\n--- Scheme {name} ---")
        scheme = fn()
        for m in scheme:
            m.describe()
        res = solve_scheme(scheme, bc, verbose=True)
        res['scheme'] = scheme
        diagnose(res, bc, verbose=True)      # automatic diagnosis / advice layer
        results[name] = res
    return bc, results


# =========================================================
# Figure 1: axial layout (colored blocks + film spanning both sides)
# =========================================================

def fig1_layout(bc, results, name):
    res = results[name]
    scheme = res['scheme']
    fig, ax = plt.subplots(figsize=(12, 5))
    z_mm = lambda z: z * 1e3

    # --- centerline: liner wall (thickness exaggerated x8, noted) ---
    # physical stack order (review fix): cold face (y=0) -> metal -> TBC
    # -> hot face (top). Drawing in list order from y=0 put the TBC on
    # the cold side and misplaced the gas zone.
    wall_y = 0.0
    draw_seq = ([lay for lay in bc['wall_layers'] if lay[2] != 'tbc'] +
                [lay for lay in bc['wall_layers'] if lay[2] == 'tbc'])
    for (nm, t, typ, prm) in draw_seq:
        h_draw = max(t * 1e3 * 8.0, 0.8)   # x8 exaggeration; min visible height
        ax.add_patch(Rectangle((0, wall_y), z_mm(L_TOTAL), h_draw,
                               facecolor=COLOR.get(nm, '#9ca3af'),
                               edgecolor='#334155', lw=0.8))
        ax.annotate(f'{nm} {t*1e3:.1f}mm', xy=(z_mm(L_TOTAL)*0.72,
                    wall_y + h_draw/2), fontsize=8, va='center')
        wall_y += h_draw
    wall_y_top = wall_y                 # hot face = top of the stack
    ax.annotate('liner wall (thickness x8, not to scale)', xy=(z_mm(L_TOTAL)*0.13,
                -2.55), fontsize=7.5, color='#64748b', ha='center')

    # --- hot side (above): gas background + film rows spanning the wall ---
    ax.add_patch(Rectangle((0, wall_y_top), z_mm(L_TOTAL), 2.6,
                           facecolor='#fee2e2', alpha=0.5, zorder=0))
    ax.annotate('gas side', xy=(z_mm(L_TOTAL)*0.05, wall_y_top + 2.2),
                fontsize=10, color='#b91c1c')
    min_span = bc.get('min_span', 1e-3)
    for m in scheme:
        if m.kind == 'film' and m.span >= min_span:   # disabled modules not drawn
            c = COLOR.get(m.name, '#ef4444')
            zc = z_mm(m.z0)
            # spanning block: cold side (below) -> wall -> gas side (above)
            for (y0, y1) in [(wall_y_top, wall_y_top + 2.0),
                             (-1.8, 0.0)]:
                ax.add_patch(Rectangle((zc - 3, y0), 6, y1 - y0,
                                       facecolor=c, alpha=0.85,
                                       edgecolor='#7f1d1d', lw=1.0,
                                       zorder=5))
            ax.annotate('', xy=(zc, wall_y_top + 2.3),
                        xytext=(zc, -2.1),
                        arrowprops=dict(arrowstyle='->', color=c, lw=1.6))
            ax.annotate(f"{m.name}\nz={zc:.0f}mm", xy=(zc, wall_y_top + 2.45),
                        fontsize=8, ha='center', color='#7f1d1d')

    # --- cold side (below): colored blocks per cold-side module ---
    path_by_mod = {id(mm): r for mm, _, r in res['path']}
    for m in scheme:
        if m.kind != 'cold' or m.span < min_span:
            continue
        c = COLOR.get(m.name, '#93c5fd')
        ax.add_patch(Rectangle((z_mm(m.z0), -1.8), z_mm(m.span), 1.8,
                               facecolor=c, alpha=0.85, edgecolor='#1e3a5f',
                               lw=1.0, zorder=4))
        label = f'{m.name}\n{z_mm(m.z0):.0f}-{z_mm(m.z1):.0f}mm ' \
                f'({m.span*1e3:.0f}mm)'
        r = path_by_mod.get(id(m))
        if r is not None:
            label += f"\nh={r['h']:.0f} W/m2K"
        ax.annotate(label, xy=(z_mm((m.z0 + m.z1)/2), -0.9),
                    fontsize=8, ha='center', va='center', zorder=6)
    if not any(m.kind == 'cold' and m.span >= min_span for m in scheme):
        ax.annotate('(no cold-side module — physically invalid)',
                    xy=(z_mm(L_TOTAL)/2, -0.9),
                    fontsize=9, ha='center', color='red')

    # --- axial zone lines ---
    for zname, z in ZONES:
        ax.axvline(z_mm(z), color='#94a3b8', ls='--', lw=0.8, alpha=0.7)
        ax.annotate(zname, xy=(z_mm(z), wall_y_top + 2.9), fontsize=8,
                    ha='center', color='#64748b')

    ax.set_xlim(-4, z_mm(L_TOTAL) + 4)
    ax.set_ylim(-4.6, wall_y_top + 3.4)
    ax.set_xlabel('axial distance from dome s [mm]')
    ax.set_yticks([])
    ax.set_title(f'Scheme layout: {name}\n'
                 f'above wall = gas side (film rows span both sides)   '
                 f'below wall = serial cold-side modules')
    ax.spines[['left', 'right', 'top']].set_visible(False)
    fig.tight_layout()
    fp = os.path.join(OUTDIR, f'fig1_layout_{name}.png')
    fig.savefig(fp, dpi=150)
    plt.close(fig)
    print(f'  [fig1] {fp}')
    return fp


# =========================================================
# Figure 2: coolant path waterfall (temperature steps + dP accumulation)
# =========================================================

def fig2_path(bc, res, name):
    steps = [('T2 inlet', bc['T2'], 0.0)]
    for m, st, r in res['path']:
        steps.append((m.name, r['T_out'], r['dP']))
    for nm, fr in res['film_res'].items():
        if fr.get('active'):
            steps.append((nm + ' (orifice)', steps[-1][1], fr.get('dP', 0.0)))
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.6))
    xs = np.arange(len(steps))
    # temperature steps
    ax1.step(xs, [s[1] for s in steps], where='post', color='#1f77b4',
             lw=2, marker='o')
    for i, s in enumerate(steps):
        ax1.annotate(f'{s[1]:.0f}K', xy=(i, s[1]), xytext=(0, 6),
                     textcoords='offset points', ha='center', fontsize=8)
    ax1.set_xticks(xs)
    ax1.set_xticklabels([s[0] for s in steps], fontsize=8,
                        rotation=18, ha='right')
    ax1.set_ylabel('coolant temperature [K]')
    ax1.set_title(f'Coolant temperature along the path (serial, feeds film)\n{name}')
    ax1.grid(alpha=0.3)
    # dP accumulation
    cum = np.cumsum([s[2] for s in steps]) / 1e3
    ax2.bar(xs, [s[2]/1e3 for s in steps], color='#f59e0b', alpha=0.8,
            label='per-module dP')
    ax2.plot(xs, cum, 'o-', color='#b91c1c', lw=2, label='cumulative dP')
    ax2.axhline(bc['dp_budget']/1e3, color='r', ls=':', lw=1.5)
    ax2.annotate(f'dP budget {bc["dp_budget"]/1e3:.1f}kPa',
                 xy=(0.02, bc['dp_budget']/1e3),
                 xycoords=('axes fraction', 'data'), va='bottom',
                 fontsize=9, color='r')
    ax2.set_xticks(xs)
    ax2.set_xticklabels([s[0] for s in steps], fontsize=8,
                        rotation=18, ha='right')
    ax2.set_ylabel('dP [kPa]')
    ax2.set_title('Per-module pressure drop and accumulation (serial-path cost)')
    ax2.legend(fontsize=9)
    ax2.grid(alpha=0.3)
    fig.tight_layout()
    fp = os.path.join(OUTDIR, f'fig2_coolant_path_{name}.png')
    fig.savefig(fp, dpi=150)
    plt.close(fig)
    print(f'  [fig2] {fp}')


# =========================================================
# Figure 3: scheme comparison (four indicators)
# =========================================================

def fig3_compare(bc, results):
    names = list(results.keys())
    short = [n.split('_')[0] for n in names]
    T_metal = [results[n]['T_metal'] for n in names]
    T_ad = [results[n]['T_ad'] for n in names]
    dP = [results[n]['dP_total']/1e3 for n in names]
    Tc_rise = [results[n]['T_c_rise'] for n in names]
    fig, axes = plt.subplots(1, 4, figsize=(14, 4))
    for ax, vals, ttl, fmt in zip(
            axes,
            [T_metal, T_ad, dP, Tc_rise],
            ['T_metal [K]\n(wall temperature, lower is better)',
             'T_ad [K]\n(adiabatic wall temperature)',
             'dP_total [kPa]\n(lower is better)',
             'coolant temperature rise [K]\n(film supply penalty)'],
            ['{:.0f}', '{:.0f}', '{:.2f}', '{:.1f}']):
        bars = ax.bar(short, vals, color=['#0d9488', '#94a3b8',
                                          '#3b82f6', '#f59e0b'][:len(names)],
                      alpha=0.85)
        for b, v in zip(bars, vals):
            ax.annotate(fmt.format(v), xy=(b.get_x() + b.get_width()/2, v),
                        xytext=(0, 4), textcoords='offset points',
                        ha='center', fontsize=9)
        if 'T_metal' in ttl:
            ax.axhline(1198, color='r', ls=':', lw=1.5)
            ax.annotate('T_allow=1198K', xy=(0.55, 1198),
                        xycoords=('axes fraction', 'data'),
                        va='bottom', fontsize=8, color='r')
        ax.set_title(ttl, fontsize=10)
        ax.grid(alpha=0.3, axis='y')
    fig.suptitle('Scheme comparison (reference-engine design point, '
                 'single-node lumped basis)', fontsize=12)
    fig.tight_layout()
    fp = os.path.join(OUTDIR, 'fig3_scheme_comparison.png')
    fig.savefig(fp, dpi=150)
    plt.close(fig)
    print(f'  [fig3] {fp}')


# =========================================================
# Self-test: A1 cross-check against the multi-station reference
# =========================================================

def selftest(results):
    print('\n' + '=' * 68)
    print('  Self-test: A1 (correlation cross-check) vs independent '
          'multi-station reference')
    print('=' * 68)
    a1 = next(v for k, v in results.items() if k.startswith('A1'))
    diff = a1['T_metal'] - 750.0
    ok = abs(diff) <= 25.0
    print(f"  Multi-station reference (station-by-station, row-by-row): "
          f"T_metal_max=750K (frozen T_c=T2, ribbed)")
    print(f"  A1 (single-node lumped): T_metal={a1['T_metal']:.1f}K  "
          f"difference={diff:+.1f}K")
    print(f"  Sources of the difference (convention, not error):")
    print(f"    1) the reference freezes T_c=T2; A1 includes the coolant "
          f"temperature rise (+{a1['T_c_rise']:.1f}K) -> wall T runs higher")
    print(f"    2) single-node lumping vs station-by-station (same hot-station "
          f"h_g, area lumping approximation)")
    print(f"  Criterion: |difference| <= 25K (coolant-rise dominated) -> "
          f"{'PASS' if ok else 'FAIL'}")
    return ok


if __name__ == '__main__':
    bc, results = run_all()
    print('\n' + '=' * 68)
    print('  Generating figures')
    print('=' * 68)
    fig1_layout(bc, results, list(results.keys())[1])   # A1 layout (ribbed + film)
    fig1_layout(bc, results, list(results.keys())[2])   # B layout (most complete blocks)
    for name in results:
        fig2_path(bc, results[name], name)
    fig3_compare(bc, results)
    selftest(results)
    print('\n--- cooling_stack V1 run complete ---')
