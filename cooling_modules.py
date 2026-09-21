# -*- coding: utf-8 -*-
"""
=================================================================
 cooling_modules.py — cooling-module library
=================================================================
Architecture:
  - Serial-path topology: coolant flows through the scheme's module
    list in order, absorbing heat and accumulating dP segment by
    segment; the path-end state feeds the film module (T_c_film).
  - Two module kinds:
      kind='cold'  cold-side module (backside convection / ribbed
                   channel / impingement) — changes h and absorbs heat
      kind='film'  hot-side film row — changes the driving temperature
                   T_ad
    Wall layers (TBC / metal) are pure thermal resistances and live in
    cooling_solver, not here.
  - Each module: own parameters + axial span (z0,z1) + solve() returning
    h / A_eff / dP / q / outlet state / formula provenance.
    span length 0 = module not used.
  - Modules are zero-coupled: add, remove, re-parameterize and reorder
    freely; the solver contains no cooling formula at all.

Formula provenance policy: every correlation carries source / expected
error / bias direction, printed in verbose output.
=================================================================
"""

import numpy as np
from gas_props import air_props

# Generic orifice discharge coefficient (sharp-edged, typical value,
# error +/-10%; low Cd -> high dP, i.e. conservative side)
CD_DEFAULT = 0.8


# =================================================================
# Base class
# =================================================================

class CoolingModule:
    """Cooling-module base class. span=(z0,z1) [m] axial coverage;
    zero length = not used."""
    kind = 'base'

    def __init__(self, name, z0=0.0, z1=0.0, **params):
        self.name = name
        self.z0, self.z1 = z0, z1
        self.params = params

    @property
    def span(self):
        return self.z1 - self.z0

    def describe(self):
        print(f"  [{self.kind:5s}] {self.name}  span={self.span*1e3:.1f}mm")
        for k, v in self.params.items():
            if v is None:      # derived quantity (back-calculated), skip
                continue
            print(f"          {k} = {v}")


# =================================================================
# Cold-side module 1: smooth backside convection (annulus passage)
# =================================================================

class BacksideConvection(CoolingModule):
    """Coolant flows through the annular gap between liner outer wall
    and casing inner wall; smooth-wall convection.

    Correlation: Dittus-Boelter  Nu = 0.023*Re^0.8*Pr^0.4  (Re>2300)
        source: Incropera, Fundamentals of Heat and Mass Transfer,
                turbulent smooth pipe
        error: +/-25%
        bias: non-circular annulus cross-section and heat-flux
              direction corrections not included -> h optimistic
              (wall temperature biased high)
    Friction dP: Blasius  f=0.316*Re^-0.25, dP=f*(L/D)*rho*V^2/2
        source: same, error +/-15% (smooth pipe)
    """

    kind = 'cold'

    def __init__(self, name='smooth backside', z0=0.0, z1=0.0,
                 D_hyd=None, A_pass=None, **kw):
        super().__init__(name, z0, z1,
                         D_hyd=D_hyd, A_pass=A_pass, **kw)

    def solve(self, st, bc):
        """st: dict(T_c, P_c, m_c);  bc: global boundaries. Returns dict."""
        T_c, P_c, m_c = st['T_c'], st['P_c'], st['m_c']
        prop = air_props(T_c, P_c)
        D_h = self.params['D_hyd'] or bc['D_hyd']
        A_p = self.params['A_pass'] or bc['A_annulus']
        V = m_c / (prop['rho'] * A_p)
        Re = prop['rho'] * V * D_h / prop['mu']
        if Re < 2300:
            # laminar: Nu=3.66 constant wall temperature (Incropera);
            # flagged only — designs should avoid this regime
            Nu = 3.66
            regime = 'laminar (Re<2300, D-B invalid)'
        else:
            Nu = 0.023 * Re**0.8 * prop['Pr']**0.4
            regime = 'turbulent'
        h = Nu * prop['k'] / D_h
        A_wall = np.pi * bc['D_liner'] * self.span
        # Heat-flux split: serial absorption, coolant temperature rise
        # self-limiting (denominator holds hA/2mcp correction)
        hA = h * A_wall
        q = hA * (st['T_wall_cold'] - T_c) / (1.0 + hA / (2.0 * m_c * prop['cp']))
        T_out = T_c + q / (m_c * prop['cp'])
        f = 0.316 * Re**-0.25 if Re > 2300 else 64.0 / max(Re, 1.0)
        dP = f * (self.span / D_h) * prop['rho'] * V**2 / 2.0
        return {'h': h, 'A_eff': A_wall, 'q': q, 'dP': dP,
                'T_out': T_out, 'cp': prop['cp'], 'Re': Re, 'Nu': Nu, 'regime': regime,
                'formula': f'D-B {regime}: Nu={Nu:.1f}'}


# =================================================================
# Cold-side module 2: ribbed channel (ribbed annulus)
# =================================================================

class RibbedChannel(BacksideConvection):
    """Smooth backside + rib enhancement on the inner wall.

    Enhancement factor (Webb-type):
        enh = (1 + 2.64*e^0.36*(p/e)^-0.46*Pr^0.24) * angle_factor
        angle_factor = 0.7 + 0.3*sin(angle)
        source: Webb 1971-type correlation
        error: +/-30%
        bias: extrapolated beyond rib-parameter range, enh optimistic
    Design note: annulus ribs are rarely machined at this engine scale,
    so default schemes keep span=0 (module retained for sensitivity
    studies); the cross-check scheme A1 enables it deliberately.
    """

    kind = 'cold'

    def __init__(self, name='ribbed channel', z0=0.0, z1=0.0,
                 rib_height=0.5e-3, rib_pitch=10e-3, rib_angle=90.0, **kw):
        super().__init__(name, z0, z1, rib_height=rib_height,
                         rib_pitch=rib_pitch, rib_angle=rib_angle, **kw)

    def solve(self, st, bc):
        res = super().solve(st, bc)
        e = self.params['rib_height']
        p = self.params['rib_pitch']
        ang = self.params['rib_angle']
        prop = air_props(st['T_c'], st['P_c'])
        enh = (1.0 + 2.64 * (e / bc['D_hyd'])**0.36
               * (p / e)**-0.46 * prop['Pr']**0.24)
        enh *= 0.7 + 0.3 * np.sin(ang * np.pi / 180.0)
        res['h'] *= enh
        res['q'] *= enh
        res['T_out'] = st['T_c'] + res['q'] / (st['m_c'] * prop['cp'])
        res['enhancement'] = enh
        res['formula'] += f'  x Webb enhancement {enh:.3f}'
        return res


# =================================================================
# Cold-side module 3: impingement cooling
# =================================================================

class ImpingementJet(CoolingModule):
    """Coolant impinges normally on the liner outer wall through one or
    more rows of impingement holes.

    Correlation: mean Nu = 0.72*Re_j^0.75*(H/d)^-0.1*Pr^(1/3)
        source: Lefebvre, Gas Turbine Combustion, Ch.5 (single-row
        mean correlation, Behbahani & Goldstein type)
        validity: H/d ~ 5-15, Re_j ~ 1e3-5e4
        error: +/-25%
        bias: single-row flat-plate basis; array + crossflow
              corrections not included -> h optimistic; small engines
              run at H/d below the validity range (extrapolated)
    Impingement-hole pressure drop: dP_j = rho/2*(V_j/Cd)^2, Cd=0.8
    Hole diameter back-calculated from the coolant budget:
        d_j = sqrt( 4*m_c/(pi*n_j*rho*V_j) )
        (flow and diameter are not independent — given V_j or dP_j,
        the other is back-calculated)

    Parameters:
        dP_j   impingement-hole allowable pressure drop [Pa] (default driver)
        n_j    total number of impingement holes
        H      jet-to-target-plate distance [m]
        V_j    (optional) jet velocity given directly, dP back-calculated
    """

    kind = 'cold'

    def __init__(self, name='impingement', z0=0.0, z1=0.0,
                 dP_j=5e3, n_j=15, H=12.8e-3, V_j=None, **kw):
        super().__init__(name, z0, z1, dP_j=dP_j, n_j=n_j, H=H, V_j=V_j, **kw)

    def solve(self, st, bc):
        T_c, P_c, m_c = st['T_c'], st['P_c'], st['m_c']
        prop = air_props(T_c, P_c)
        Cd = CD_DEFAULT
        if self.params['V_j'] is not None:
            V_j = self.params['V_j']
            dP_j = prop['rho'] / 2.0 * (V_j / Cd)**2
        else:
            V_j = Cd * np.sqrt(2.0 * self.params['dP_j'] / prop['rho'])
            dP_j = self.params['dP_j']
        n_j = self.params['n_j']
        A_j_total = m_c / (prop['rho'] * V_j)          # total jet area
        d_j = np.sqrt(4.0 * A_j_total / (np.pi * n_j))  # single-hole diameter
        H = self.params['H']
        Re_j = prop['rho'] * V_j * d_j / prop['mu']
        H_d = H / d_j
        Nu = 0.72 * Re_j**0.75 * H_d**-0.1 * prop['Pr']**(1.0/3.0)
        h = Nu * prop['k'] / d_j
        A_wall = np.pi * bc['D_liner'] * self.span
        hA = h * A_wall
        q = hA * (st['T_wall_cold'] - T_c) / (1.0 + hA / (2.0 * m_c * prop['cp']))
        T_out = T_c + q / (m_c * prop['cp'])
        return {'h': h, 'A_eff': A_wall, 'q': q, 'dP': dP_j,
                'T_out': T_out, 'cp': prop['cp'], 'Re': Re_j, 'Nu': Nu,
                'd_j': d_j, 'V_j': V_j, 'H_d': H_d,
                'regime': f'impingement Re_j={Re_j:.0f}, H/d={H_d:.1f}',
                'formula': f'impingement: Nu={Nu:.0f} (0.72Re^0.75(H/d)^-0.1Pr^1/3)'}


# =================================================================
# Hot-side module: film-cooling row
# =================================================================

class FilmRow(CoolingModule):
    """One film-cooling row, located at the end of the coolant path:
    coolant heated by the upstream cold-side modules leaves through
    this row into the gas path and forms the film.

    Film effectiveness (Lefebvre simplified form):
        M < 0.5:  eta = exp(-0.05*x/s)
        M >= 0.5: eta = 1/(1 + 0.1*(x/s)^0.8*M^0.5)
        source: Lefebvre 2010 simplified correlation; error +/-20-30%
        bias: the two branches are discontinuous at M=0.5
              (0.607 -> 0.691), a known feature; x/s uses an
              equivalent-slot-height basis (s_ref default 1 mm)
    Blowing ratio M, two parallel bases:
        A annulus basis:  M = (m_c/m_g)*(A_ref/A_annulus)
        B orifice basis:  M = rho_c*V_j/(rho_g*V_g),
                          V_j = m_film/(rho_c*Cd*A_holes)
          (physical orifice jet; the primary basis is to be fixed once
          the 3D geometry is frozen)
    Row pressure drop: dP = rho_c/2*(V_j/Cd)^2 (from hole geometry and
        flow, for budget checking)

    Multi-row superposition (done in solver): eta_overall = 1 - prod(1-eta_i)
        (source: Lefebvre multi-row staggered superposition, error +/-25%)
    """

    kind = 'film'

    def __init__(self, name='film row', z0=0.0, z1=None,
                 d=2e-3, n_rows=12, n_per_row=20, Cd=CD_DEFAULT,
                 s_ref=1e-3, **kw):
        # span semantics unified with the other modules (same treatment
        # as the ribbed channel):
        #   span=0 -> row disabled (thin dP budgets rarely afford film;
        #             module retained, set span>0 to re-enable)
        #   span>0 -> enabled; row axial position = z0 (upstream hole
        #             edge), x = z_eval - z0
        if z1 is None:
            z1 = z0
        super().__init__(name, z0, z1, d=d, n_rows=n_rows,
                         n_per_row=n_per_row, Cd=Cd, s_ref=s_ref, **kw)

    def solve(self, st, bc):
        """st requires m_film (coolant flow through this row);
        bc requires z_eval. Returns single-row eta, dual-basis M, row dP."""
        if self.span < bc.get('min_span', 1e-3):
            return {'eta': 0.0, 'active': False, 'x': 0.0,
                    'note': f"span={self.span*1e3:.1f}mm<MIN_SPAN, disabled "
                            "(film rows are sparingly used on small "
                            "engines; module retained for reference)"}
        z_eval = bc['z_eval']
        x = z_eval - self.z0          # downstream distance row -> eval point
        if x <= 0:
            return {'eta': 0.0, 'active': False, 'x': x,
                    'note': 'film row is downstream of the evaluation '
                            'station; no film contribution at this station'}
        s = self.params['s_ref']
        d = self.params['d']
        n_tot = self.params['n_rows'] * self.params['n_per_row']
        Cd = self.params['Cd']
        m_film = st['m_film']
        prop_c = air_props(st['T_c'], st['P_c'])
        A_holes = n_tot * np.pi * d**2 / 4.0
        V_jet = m_film / (prop_c['rho'] * Cd * A_holes)
        dP = prop_c['rho'] / 2.0 * (V_jet / Cd)**2
        # M, two bases
        M_annulus = (m_film / bc['m_gas']) * (bc['A_ref'] / bc['A_annulus'])
        M_jet = prop_c['rho'] * V_jet / (bc['rho_g'] * bc['V_g'])
        out = {'active': True, 'x': x, 's': s, 'x_s': x / s,
               'M_annulus': M_annulus, 'M_jet': M_jet,
               'V_jet': V_jet, 'dP': dP, 'A_holes': A_holes,
               'd_eq': float(np.sqrt(4*A_holes/(np.pi*n_tot))),
               'n_tot': n_tot}
        for tag, M in (('annulus', M_annulus), ('jet', M_jet)):
            if M < 0.5:
                eta = np.exp(-0.05 * x / s)
            else:
                eta = 1.0 / (1.0 + 0.1 * (x / s)**0.8 * M**0.5)
            out[f'eta_{tag}'] = float(eta)
        out['eta'] = out['eta_annulus']   # primary = annulus basis; jet basis parallel
        return out
