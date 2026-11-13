"""
Simplified physics scheme: all-or-nothing saturation adjustment,
K-gradient turbulence, no radiation, only vapour and cloud liquid water.
"""

import numpy

from . import wasp
from .physics_base import PhysicsBase


class SimplePhysics(PhysicsBase):
    """
    Simplified physics scheme.

    Features:
      - all-or-nothing saturation adjustment (liquid only)
      - K-gradient turbulence with constant K
      - mass-flux shallow convection (prescribed half-sine profile)
      - WASP or prescribed surface fluxes
      - only vapour (v) and cloud liquid water (c)
    """

    def __init__(self, dt, method, name, tag, case, grid, prognostic_variables,
                 attrs=None,
                 mf_max=0.05, mf_dtheta=.5, mf_rt_factor=1.2,
                 do_adjustment=True, do_shallow=True, do_turb=True):
        super().__init__(dt, method, name, tag, case, grid, prognostic_variables)
        self.mf_max = mf_max
        self.mf_dtheta = mf_dtheta
        self.mf_rt_factor = mf_rt_factor
        self.do_adjustment = do_adjustment
        self.do_shallow = do_shallow
        self.do_turb = do_turb

    def setup(self, init_state, duration):
        pass

    # ------------------------------------------------------------------
    #  step helpers
    # ------------------------------------------------------------------
    def _k_diff(self, theta, rv, u, v, ctx, sfth, sfrv, sfu, sfv, timestep):
        """
        K-gradient diffusion tendencies for all variables in one pass.
        Uses an implicit tridiagonal solver following PHYEX's TRIDIAG_THERMO.

        Parameters
        ----------
        theta : ndarray (nz,) — potential temperature
        rv : ndarray (nz,) — water vapour mixing ratio
        u, v : ndarray (nz,) — wind components
        ctx : dict — context from _preparation
        sfth, sfrv, sfu, sfv : float — surface fluxes
        timestep : float — time step (s)

        Returns
        -------
        tuple of ndarray — (tend_theta, tend_rv, tend_u, tend_v)
        """
        nz = ctx['nkt']
        dz_mass = ctx['dz_mass']
        z_flux = ctx['z_flux']
        dzz = ctx['dzz']

        # Height above ground for flux levels
        if self.grid.ascending:
            z_ground = z_flux[0]
        else:
            z_ground = z_flux[-1]
        height = numpy.abs(z_flux - z_ground)

        # Blackadar length
        lasymp = 200.
        l = self.cst.Karman * height[1:-1] / (1. + self.cst.Karman * height[1:-1] / lasymp)

        # Shear
        du_dz = ((u[1:] - u[:-1]) / dz_mass) ** 2 + \
                ((v[1:] - v[:-1]) / dz_mass) ** 2
        du_dz = numpy.sqrt(du_dz)

        # Virtual potential temperature (buoyancy)
        theta_v = theta * (1.0 + (self.cst.Rv / self.cst.Rd - 1.0) * rv)
        theta_v_flux = 0.5 * (theta_v[1:] + theta_v[:-1])

        # Brunt-Väisälä frequency squared at flux points
        N2 = (self.cst.g / theta_v_flux) * (theta_v[1:] - theta_v[:-1]) / dz_mass

        # Gradient Richardson number (only stable stratification reduces K)
        Ri = numpy.maximum(N2, 0.0) / (du_dz**2 + 1e-10)

        # Stability function (Louis 1979, simplified)
        f_Ri = 1.0 / (1.0 + 10.0 * Ri)**2

        K = l**2 * numpy.abs(du_dz) * f_Ri
        K = numpy.maximum(K, 1e-4)
        K = numpy.minimum(K, 200)

        rhodJ = ctx['rhod'] * dzz

        # MZM_RHODJ at flux points (0.5*(rhodJ[k] + rhodJ[k-1]))
        MZM_RHODJ = numpy.zeros(nz + 1)
        MZM_RHODJ[1:-1] = 0.5 * (rhodJ[1:] + rhodJ[:-1])
        MZM_RHODJ[0] = rhodJ[0]
        MZM_RHODJ[-1] = rhodJ[-1]

        # PDZZ at flux points: distance between adjacent mass levels
        PDZZ = numpy.ones(nz + 1)
        PDZZ[1:-1] = dz_mass
        PDZZ[0] = dz_mass[0]
        PDZZ[-1] = dz_mass[-1]

        # PDFDDTDZ = K at flux points
        PDFDDTDZ = numpy.zeros(nz + 1)
        PDFDDTDZ[1:-1] = -K
        PDFDDTDZ[0] = -K[0]
        PDFDDTDZ[-1] = -K[-1]

        ZRHODJ_DF = MZM_RHODJ * PDFDDTDZ / (PDZZ ** 2)

        PIMPL = 1.0

        def _compute_fluxes(field, surface_flux):
            PF = numpy.zeros(nz + 1)
            PF[0] = surface_flux
            PF[-1] = 0.
            PF[1:-1] = - K * (field[1:] - field[:-1]) / dz_mass
            return PF

        def _tridiag_solver_thermo(PVARM, PF):
            A = numpy.zeros(nz)
            B = numpy.zeros(nz)
            C = numpy.zeros(nz)
            Y = numpy.zeros(nz)

            # Bottom level
            Y[0] = rhodJ[0] * PVARM[0] / timestep \
                - MZM_RHODJ[1] * PF[1] / PDZZ[1] \
                + MZM_RHODJ[0] * PF[0] / PDZZ[0] \
                + ZRHODJ_DF[1] * PIMPL * PVARM[1] \
                - ZRHODJ_DF[1] * PIMPL * PVARM[0]

            B[0] = rhodJ[0] / timestep - ZRHODJ_DF[1] * PIMPL
            C[0] = ZRHODJ_DF[1] * PIMPL

            # Interior levels
            for k in range(1, nz - 1):
                Y[k] = rhodJ[k] * PVARM[k] / timestep \
                    - MZM_RHODJ[k + 1] * PF[k + 1] / PDZZ[k + 1] \
                    + MZM_RHODJ[k] * PF[k] / PDZZ[k] \
                    + ZRHODJ_DF[k + 1] * PIMPL * PVARM[k + 1] \
                    - ZRHODJ_DF[k + 1] * PIMPL * PVARM[k] \
                    - ZRHODJ_DF[k] * PIMPL * PVARM[k] \
                    + ZRHODJ_DF[k] * PIMPL * PVARM[k - 1]
                A[k] = ZRHODJ_DF[k] * PIMPL
                B[k] = rhodJ[k] / timestep \
                    - ZRHODJ_DF[k + 1] * PIMPL - ZRHODJ_DF[k] * PIMPL
                C[k] = ZRHODJ_DF[k + 1] * PIMPL

            # Top level
            Y[-1] = rhodJ[-1] * PVARM[-1] / timestep \
                - MZM_RHODJ[-1] * PF[-1] / PDZZ[-1] \
                + MZM_RHODJ[-2] * PF[-2] / PDZZ[-2] \
                - ZRHODJ_DF[-2] * PIMPL * PVARM[-1] \
                + ZRHODJ_DF[-2] * PIMPL * PVARM[-2]
            A[-1] = ZRHODJ_DF[-2] * PIMPL
            B[-1] = rhodJ[-1] / timestep - ZRHODJ_DF[-2] * PIMPL

            # Thomas algorithm
            gam = numpy.zeros(nz)
            bet = B[0]
            PVARP = numpy.zeros(nz)
            PVARP[0] = Y[0] / bet

            for k in range(1, nz):
                gam[k] = C[k - 1] / bet
                bet = B[k] - A[k] * gam[k]
                PVARP[k] = (Y[k] - A[k] * PVARP[k - 1]) / bet

            for k in range(nz - 2, -1, -1):
                PVARP[k] = PVARP[k] - gam[k + 1] * PVARP[k + 1]

            return (PVARP - PVARM) / timestep

        PF_theta = _compute_fluxes(theta, sfth)
        PF_rv = _compute_fluxes(rv, sfrv)
        PF_u = _compute_fluxes(u, sfu)
        PF_v = _compute_fluxes(v, sfv)

        return (
            _tridiag_solver_thermo(theta, PF_theta),
            _tridiag_solver_thermo(rv, PF_rv),
            _tridiag_solver_thermo(u, PF_u),
            _tridiag_solver_thermo(v, PF_v),
        )

    def _preparation(self, state):
        """Compute grid, density and reference fields, return a context dict."""
        ctx = {}

        ctx['pressure'] = self.grid.get_pressure(
            'MASS', state, self.prognostic_variables)
        ctx['z_mass'] = self.grid.get_altitude(
            'MASS', state, self.prognostic_variables)
        ctx['z_flux'] = self.grid.get_altitude(
            'FLUX', state, self.prognostic_variables)

        dzz = numpy.diff(ctx['z_flux'])
        if self.grid.ascending:
            if dzz[-1] == numpy.inf:
                dzz[-1] = dzz[-2]
        else:
            dzz = -dzz
            if dzz[0] == numpy.inf:
                dzz[0] = dzz[1]
        ctx['dzz'] = dzz
        dz_mass = numpy.diff(ctx['z_mass'])
        if not self.grid.ascending:
            dz_mass = -dz_mass
        ctx['dz_mass'] = dz_mass
        ctx['nkt'] = len(ctx['pressure'])

        ctx['exner'] = (ctx['pressure'] / 1.E5) ** (self.cst.Rd / self.cst.Cpd)

        if 'T' in self.prognostic_variables:
            ctx['theta'] = state['T'] / ctx['exner']
            ctx['temperature'] = state['T']
        else:
            ctx['theta'] = state['Theta']
            ctx['temperature'] = ctx['exner'] * ctx['theta']

        if 'qv' in self.prognostic_variables:
            qt = state['qv'] + state['qc']
            ctx['qdm'] = 1. - qt
            gas_constant = (self.cst.Rd
                            + state['qv'] * (self.cst.Rv - self.cst.Rd)
                            - state['qc'] * self.cst.Rd)
            ctx['rho'] = ctx['pressure'] / (gas_constant * ctx['temperature'])
            ctx['rhod'] = ctx['rho'] * ctx['qdm']
        else:
            rt = state['rv'] + state['rc']
            ctx['qdm'] = 1. / (1. + rt)
            gas_constant = (self.cst.Rd + state['rv'] * self.cst.Rv) / (1. + rt)
            ctx['rho'] = ctx['pressure'] / (gas_constant * ctx['temperature'])
            ctx['rhod'] = ctx['pressure'] / (
                (self.cst.Rd + state['rv'] * self.cst.Rv)
                * ctx['theta'] * ctx['exner'])

        return ctx

    def _adjust(self, theta, rv, rc, ctx):
        """
        Adjustment
        """
        T_cur = ctx['exner'] * theta
        logT = numpy.log(T_cur.clip(min=1.))
        es = numpy.exp(self.cst.alpw - self.cst.betaw / T_cur
                       - self.cst.gamw * logT)
        es = numpy.minimum(es, 0.99 * ctx['pressure'])

        eps = self.cst.Rd / self.cst.Rv
        rsat = eps * es / (ctx['pressure'] - (1. - eps) * es)

        drc = rv - rsat
        drc = numpy.clip(drc, -rc, rv)
        return drc

    def _surface(self, state, ctx):
        """Compute surface fluxes, return (sfu, sfv, sfth, sfrv)."""
        klevgrd = 0 if self.grid.ascending else -1

        if self.case.surface_forcing_temp in ('none', 'ts') or \
           self.case.surface_forcing_moisture in ('none', 'beta', 'mrsos') or \
           self.case.surface_forcing_wind in ('none',):
            need_scheme = True
        else:
            need_scheme = False

        def ustar2fluxes(ustar, u, v):
            alpha = numpy.arctan2(v, u)
            return -ustar ** 2 * numpy.cos(alpha), \
                   -ustar ** 2 * numpy.sin(alpha)

        if 'Ts' in state:
            surface_temperature = state['Ts']
        else:
            if self.grid.ascending:
                surface_temperature = ctx['temperature'][0]
            else:
                surface_temperature = ctx['temperature'][-1]
        if surface_temperature > self.cst.Tt:
            latent_heat = (self.cst.LvTt
                           - (self.cst.Cpv - self.cst.Cl)
                           * (surface_temperature - self.cst.Tt))
        else:
            latent_heat = (self.cst.LsTt
                           - (self.cst.Cpv - self.cst.Ci)
                           * (surface_temperature - self.cst.Tt))

        sshf_scheme = None
        swf_scheme = None
        sfu_scheme = None
        sfv_scheme = None
        if not need_scheme:
            pass
        else:
            if self.case.surface_forcing_temp == 'none':
                from . import Phyex1DError
                raise Phyex1DError(
                    'surface_forcing_temp == none not implemented')
            if self.case.surface_forcing_moisture in ('beta', 'mrsos'):
                from . import Phyex1DError
                raise Phyex1DError(
                    'surface_forcing_moisture == beta or mrsos '
                    'not implemented')
            if self.case.surface_forcing_wind != 'none':
                from . import Phyex1DError
                raise Phyex1DError(
                    'surface_forcing_wind != none not implemented')
            windgrd = numpy.sqrt(state['u'][klevgrd] ** 2
                                 + state['v'][klevgrd] ** 2)
            exner_surf = (state['Ps'] / 1.E5) ** (self.cst.Rd / self.cst.Cpd)
            if 'qv' in self.prognostic_variables:
                qv_grd = state['qv'][klevgrd]
            else:
                qv_grd = state['rv'][klevgrd] * ctx['qdm'][klevgrd]
            precip = 0.
            wave_height = 0.
            wave_peak_period = 0.
            result = wasp.WASP_FLUX(
                ctx['temperature'][klevgrd], qv_grd,
                ctx['exner'][klevgrd], ctx['rho'][klevgrd],
                windgrd, ctx['z_mass'][klevgrd], ctx['z_mass'][klevgrd],
                state['Ts' if self.case.surface_type == 'ocean' else 'Tskin'],
                exner_surf, state['Ps'], precip,
                wave_height, wave_peak_period)
            sshf_scheme, swf_scheme, ustar = result
            sfu_scheme, sfv_scheme = ustar2fluxes(
                ustar, state['u'][klevgrd], state['v'][klevgrd])

        if self.case.surface_forcing_temp == 'none':
            if sshf_scheme is None:
                from . import Phyex1DError
                raise Phyex1DError("Surface scheme didn't compute sshf")
        elif self.case.surface_forcing_temp == 'kinematic':
            raise NotImplementedError('surface_forcing_temp=kinematic')
        elif self.case.surface_forcing_temp == 'surface_flux':
            sshf = state['sshf']
        elif self.case.surface_forcing_temp == 'ts':
            sshf = sshf_scheme
        else:
            from . import Phyex1DError
            raise Phyex1DError('Wrong surface_forcing_temp option')

        if self.case.surface_forcing_moisture == 'none':
            if swf_scheme is None:
                from . import Phyex1DError
                raise Phyex1DError("Surface scheme didn't compute swf")
            swf = swf_scheme
        elif self.case.surface_forcing_moisture == 'kinematic':
            raise NotImplementedError('surface_forcing_moisture=kinematic')
        elif self.case.surface_forcing_moisture == 'surface_flux':
            swf = state['slhf'] / latent_heat
        elif self.case.surface_forcing_moisture == 'beta':
            raise NotImplementedError('surface_forcing_moisture=beta')
        elif self.case.surface_forcing_moisture == 'mrsos':
            raise NotImplementedError('surface_forcing_moisture=mrsos')
        else:
            from . import Phyex1DError
            raise Phyex1DError('Wrong surface_forcing_moisture option')

        if self.case.surface_forcing_wind == 'none':
            sfu = sfu_scheme
            sfv = sfv_scheme
        elif self.case.surface_forcing_wind == 'z0':
            ustar = (numpy.sqrt(state['u'][klevgrd] ** 2
                                + state['v'][klevgrd] ** 2)
                     * self.cst.Karman
                     / numpy.log(ctx['z_mass'][klevgrd] / state['z0']))
            sfu, sfv = ustar2fluxes(
                ustar, state['u'][klevgrd], state['v'][klevgrd])
        elif self.case.surface_forcing_wind == 'ustar':
            sfu, sfv = ustar2fluxes(
                state['ustar'], state['u'][klevgrd], state['v'][klevgrd])
        else:
            from . import Phyex1DError
            raise Phyex1DError('Wrong surface_forcing_wind option')

        sfrv = swf / ctx['rhod'][0 if self.grid.ascending else -1]
        sfth = sshf / (self.cst.Cpd * ctx['rho'][klevgrd])

        state['sfrv'] = sfrv
        state['swf'] = swf
        state['slhf'] = swf * latent_heat
        state['sfth'] = sfth
        state['sshf'] = sshf

        return sfu, sfv, sfth, sfrv

    def _massflux(self, state, ctx, timestep):
        """
        Mass-flux shallow convection scheme with implicit tridiagonal solver
        following PHYEX's TRIDIAG_MASSFLUX.

        Updraft properties equal near-surface values (no perturbation).
        Cloud top is the first flux level above ground where environmental
        theta_l equals the updraft theta_l.  The mass-flux profile is a
        half-sine (0 at surface, 0 at cloud top, max in the middle).

        Parameters
        ----------
        timestep : float — time step (s)

        Returns
        -------
        tend_theta_l : ndarray (nkt,) — theta_l tendency (K/s)
        tend_rt      : ndarray (nkt,) — total-water tendency (kg/kg/s)
        tend_u       : ndarray (nkt,) — u-wind tendency (m/s²)
        tend_v       : ndarray (nkt,) — v-wind tendency (m/s²)
        """
        nkt = ctx['nkt']

        # Height above ground for flux levels
        if self.grid.ascending:
            z_ground = ctx['z_flux'][0]
            klev = 0
        else:
            z_ground = ctx['z_flux'][-1]
            klev = -1
        height = numpy.abs(ctx['z_flux'] - z_ground)

        # Current theta, rv, rc from state
        if 'qv' in self.prognostic_variables:
            rv = state['qv'] / ctx['qdm']
            rc = state['qc'] / ctx['qdm']
        else:
            rv = state['rv']
            rc = state['rc']

        # Liquid water potential temperature
        lv = (self.cst.LvTt
              + (self.cst.Cpv - self.cst.Cl) * (ctx['temperature'] - self.cst.Tt))
        cp = ctx['qdm'] * (self.cst.Cpd + self.cst.Cpv * rv + self.cst.Cl * rc)
        theta_l = ctx['theta'] - lv / (cp * ctx['exner']) * rc

        # Total water
        rt = rv + rc

        # Updraft properties (near-surface values with excess/factor)
        theta_l_up = theta_l[klev] + self.mf_dtheta
        rt_up = rt[klev] * self.mf_rt_factor
        u_up = state['u'][klev]
        v_up = state['v'][klev]

        # Interpolate environmental profiles to flux levels
        def to_flux(var):
            vf = numpy.zeros(nkt + 1)
            vf[1:-1] = 0.5 * (var[:-1] + var[1:])
            vf[0] = var[klev]
            vf[-1] = var[-1 if self.grid.ascending else 0]
            return vf

        theta_l_env = to_flux(theta_l)
        rt_env = to_flux(rt)
        u_env = to_flux(state['u'])
        v_env = to_flux(state['v'])

        # PBL top: first flux level above ground where theta_l_env >= theta_l_up
        if self.grid.ascending:
            ktop = nkt
            for k in range(1, nkt + 1):
                if theta_l_env[k] >= theta_l_up:
                    ktop = k
                    break
        else:
            ktop = 0
            for k in range(nkt - 1, -1, -1):
                if theta_l_env[k] >= theta_l_up:
                    ktop = k
                    break
        pbl_height = height[ktop]

        # Mass flux profile at interfaces (0 at surface, 0 at cloud top)
        mf = numpy.zeros(nkt + 1)
        mask = (height > 0.) & (height < pbl_height)
        mf[mask] = self.mf_max * numpy.sin(
            numpy.pi * height[mask] / pbl_height)

        # Explicit fluxes at flux points: F = mf * (updraft - env)
        PF_theta_l = mf * (theta_l_up - theta_l_env)
        PF_rt = mf * (rt_up - rt_env)
        PF_u = mf * (u_up - u_env)
        PF_v = mf * (v_up - v_env)

        # PDFDT = dF/d(env) at flux points = -mf
        PDFDT = -mf

        # Grid metrics for tridiagonal solver
        dzz = ctx['dzz']
        rhodJ = ctx['rhod'] * dzz

        MZM_RHODJ = numpy.zeros(nkt + 1)
        MZM_RHODJ[1:-1] = 0.5 * (rhodJ[1:] + rhodJ[:-1])
        MZM_RHODJ[0] = rhodJ[klev]
        MZM_RHODJ[-1] = rhodJ[-1 if self.grid.ascending else 0]

        PDZZ = numpy.ones(nkt + 1)
        PDZZ[1:-1] = ctx['dz_mass']
        PDZZ[0] = ctx['dz_mass'][0]
        PDZZ[-1] = ctx['dz_mass'][-1]

        ZRHODJ_DFDT_O_DZ = MZM_RHODJ * PDFDT / PDZZ

        PIMPL = 1.0

        def _tridiag_solver(PVARM, PF):
            A = numpy.zeros(nkt)
            B = numpy.zeros(nkt)
            C = numpy.zeros(nkt)
            Y = numpy.zeros(nkt)

            # Bottom level
            Y[0] = rhodJ[0] * PVARM[0] / timestep \
                - MZM_RHODJ[1] * PF[1] / PDZZ[1] \
                + MZM_RHODJ[0] * PF[0] / PDZZ[0] \
                + ZRHODJ_DFDT_O_DZ[1] * 0.5 * PIMPL * PVARM[1] \
                + ZRHODJ_DFDT_O_DZ[1] * 0.5 * PIMPL * PVARM[0]

            B[0] = rhodJ[0] / timestep \
                + ZRHODJ_DFDT_O_DZ[1] * 0.5 * PIMPL
            C[0] = ZRHODJ_DFDT_O_DZ[1] * 0.5 * PIMPL

            # Interior levels
            for k in range(1, nkt - 1):
                Y[k] = rhodJ[k] * PVARM[k] / timestep \
                    - MZM_RHODJ[k + 1] * PF[k + 1] / PDZZ[k + 1] \
                    + MZM_RHODJ[k] * PF[k] / PDZZ[k] \
                    + ZRHODJ_DFDT_O_DZ[k + 1] * 0.5 * PIMPL * PVARM[k + 1] \
                    + ZRHODJ_DFDT_O_DZ[k + 1] * 0.5 * PIMPL * PVARM[k] \
                    - ZRHODJ_DFDT_O_DZ[k] * 0.5 * PIMPL * PVARM[k] \
                    - ZRHODJ_DFDT_O_DZ[k] * 0.5 * PIMPL * PVARM[k - 1]
                A[k] = -ZRHODJ_DFDT_O_DZ[k] * 0.5 * PIMPL
                B[k] = rhodJ[k] / timestep \
                    + ZRHODJ_DFDT_O_DZ[k + 1] * 0.5 * PIMPL \
                    - ZRHODJ_DFDT_O_DZ[k] * 0.5 * PIMPL
                C[k] = ZRHODJ_DFDT_O_DZ[k + 1] * 0.5 * PIMPL

            # Top level
            Y[-1] = rhodJ[-1] * PVARM[-1] / timestep
            B[-1] = rhodJ[-1] / timestep

            # Thomas algorithm
            gam = numpy.zeros(nkt)
            bet = B[0]
            PVARP = numpy.zeros(nkt)
            PVARP[0] = Y[0] / bet

            for k in range(1, nkt):
                gam[k] = C[k - 1] / bet
                bet = B[k] - A[k] * gam[k]
                PVARP[k] = (Y[k] - A[k] * PVARP[k - 1]) / bet

            for k in range(nkt - 2, -1, -1):
                PVARP[k] = PVARP[k] - gam[k + 1] * PVARP[k + 1]

            return (PVARP - PVARM) / timestep

        return (
            _tridiag_solver(theta_l, PF_theta_l),
            _tridiag_solver(rt, PF_rt),
            _tridiag_solver(state['u'], PF_u),
            _tridiag_solver(state['v'], PF_v),
        )

    # ------------------------------------------------------------------
    #  step
    # ------------------------------------------------------------------
    def step(self, state, timestep, timestep_number):
        state['tke'] = numpy.maximum(1e-6, state['tke'])

        ctx = self._preparation(state)

        ###############################################################
        #  1. SATURATION ADJUSTMENT  (all-or-nothing, liquid only)
        ###############################################################

        if self.do_adjustment:
            if 'qv' in self.prognostic_variables:
                rv = state['qv'] / ctx['qdm']
                rc = state['qc'] / ctx['qdm']
            else:
                rv = state['rv'].copy()
                rc = state['rc'].copy()
            theta = ctx['theta'].copy()

            drc = self._adjust(theta, rv, rc, ctx)
            rv = rv - drc
            rc = rc + drc

            lv = (self.cst.LvTt
                  + (self.cst.Cpv - self.cst.Cl) * (ctx['exner'] * theta - self.cst.Tt))
            cp = ctx['qdm'] * (self.cst.Cpd + self.cst.Cpv * rv
                               + self.cst.Cl * rc)

            theta = ctx['theta'] + lv / (cp * ctx['exner']) * rc

            if 'qv' in self.prognostic_variables:
                state['qv'] = rv * ctx['qdm']
                state['qc'] = rc * ctx['qdm']
                state['T'] = theta * ctx['exner']
            else:
                state['rv'] = rv
                state['rc'] = rc
                state['Theta'] = theta

            # Recompute derived fields after adjustment
            ctx = self._preparation(state)

        ###############################################################
        #  2. SURFACE
        ###############################################################

        sfu, sfv, sfth, sfrv = self._surface(state, ctx)

        ###############################################################
        #  3. SHALLOW CONVECTION  (mass flux)
        ###############################################################
        if self.do_shallow:
            tend_theta_l, tend_rt, tend_u, tend_v = self._massflux(state, ctx, timestep)

            # Apply tendencies
            if 'Theta' in self.prognostic_variables:
                state['Theta'] += tend_theta_l * timestep
            else:
                state['T'] += tend_theta_l * ctx['exner'] * timestep

            if 'qv' in self.prognostic_variables:
                state['qv'] += tend_rt * ctx['qdm'] * timestep
                state['qv'] = numpy.maximum(state['qv'], 0.)
            else:
                state['rv'] += tend_rt * timestep
                state['rv'] = numpy.maximum(state['rv'], 0.)

            state['u'] += tend_u * timestep
            state['v'] += tend_v * timestep

        ###############################################################
        #  4. K-GRADIENT TURBULENCE  (K = (du/dz)^2)
        ###############################################################
        if self.do_turb:
            # Prepare variables in theta/r space
            if 'Theta' in self.prognostic_variables:
                theta = state['Theta']
            else:
                theta = state['T'] / ctx['exner']

            if 'rv' in self.prognostic_variables:
                rv = state['rv']
            elif 'qv' in self.prognostic_variables:
                rv = state['qv'] / ctx['qdm']

            tend_theta, tend_rv, tend_u, tend_v = self._k_diff(
                theta, rv, state['u'], state['v'], ctx,
                sfth, sfrv, sfu, sfv, timestep)

            if 'Theta' in self.prognostic_variables:
                state['Theta'] += tend_theta * timestep
            else:
                state['T'] += tend_theta * ctx['exner'] * timestep

            if 'rv' in self.prognostic_variables:
                state['rv'] += tend_rv * timestep
            elif 'qv' in self.prognostic_variables:
                state['qv'] += tend_rv * ctx['qdm'] * timestep

            state['u'] += tend_u * timestep
            state['v'] += tend_v * timestep

        return state


class SimplePhysicsTQ(SimplePhysics):
    """
    Simplified physics with T and q as prognostic variables
    """

    def __init__(self, dt, method, name, tag, case, grid,
                 attrs=None,
                 mf_max=0.05, mf_dtheta=.5, mf_rt_factor=1.2,
                 do_adjustment=True, do_shallow=True, do_turb=True):
        super().__init__(dt, method, name, tag, case, grid,
                         ['T', 'qv', 'qc', 'u', 'v', 'w', 'tke'],
                         attrs,
                         mf_max=mf_max, mf_dtheta=mf_dtheta, mf_rt_factor=mf_rt_factor,
                         do_adjustment=do_adjustment, do_shallow=do_shallow, do_turb=do_turb)


class SimplePhysicsThetaR(SimplePhysics):
    """
    Simplified physics with Theta and r as prognostic variables
    """

    def __init__(self, dt, method, name, tag, case, grid,
                 attrs=None,
                 mf_max=0.05, mf_dtheta=.5, mf_rt_factor=1.2,
                 do_adjustment=True, do_shallow=True, do_turb=True):
        super().__init__(dt, method, name, tag, case, grid,
                         ['Theta', 'rv', 'rc', 'u', 'v', 'w', 'tke'],
                         attrs,
                         mf_max=mf_max, mf_dtheta=mf_dtheta, mf_rt_factor=mf_rt_factor,
                         do_adjustment=do_adjustment, do_shallow=do_shallow, do_turb=do_turb)
