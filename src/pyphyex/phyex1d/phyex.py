"""
Implementation of PHYEX physics schemes
"""
# pylint: disable=no-member

import datetime
import importlib.util
import logging
import os
import sys
import tempfile

import ephem
import f90nml
import numpy
from ctypesForFortran import MISSING

from . import Phyex1DError
from . import wasp
from .physics_base import PhysicsBase


class PHYEX(PhysicsBase):
    """
    PHYEX-based physics scheme

    This class handles the initialisation / setup of PHYEX and ECRAD.
    """

    def __init__(self, dt, method, name, tag, case, grid, prognostic_variables,
                 pyphyex=None, pyecrad=None, namel=None, dx=0, dy=0, attrs=None):
        super().__init__(dt, method, name, tag, case, grid, prognostic_variables)
        self.pyphyex = pyphyex
        self._pyphyex = None
        self.pyecrad = pyecrad
        self._pyecrad = None
        self.namel = namel
        self.full_phyex_namel = None
        self.dx = dx
        self.dy = dy
        if self.namel is not None and not isinstance(namel, dict):
            if self.namel.endswith('.namel') and os.path.exists(self.namel):
                # namelist is a filename provided by the user
                namelist_file = self.namel
            else:
                # look for a namelis with the name provided in the user's directory,
                # then in the package directory
                user_dir = os.path.join(os.environ['HOME'], '.phyex1d', 'namelists')
                package_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'namelists')
                namelist_file = None
                for namel_dir in (user_dir, package_dir):
                    namelist_file = os.path.join(namel_dir, self.namel + '.namel')
                    if os.path.exists(namelist_file):
                        break
            if namelist_file is None:
                raise Phyex1DError('Namelist not found')

            with open(namelist_file, 'r', encoding='UTF-8') as f:
                self.namel = f90nml.read(f)

    def setup(self, init_state, duration):
        """Initialise PHYEX and ECRAD"""
        try:
            if self.pyphyex is None:
                # pyphyex.py should be found by python
                import pyphyex  # pylint: disable=import-outside-toplevel
            else:
                # We import pyphyex.py by file path
                spec = importlib.util.spec_from_file_location('pyphyex', self.pyphyex)
                pyphyex = importlib.util.module_from_spec(spec)
                sys.modules['pyphyex'] = pyphyex
                spec.loader.exec_module(pyphyex)
        except ModuleNotFoundError:
            logging.error("The module pyphyex has not been found. " +
                          "This module can be built with the PHYEX package " +
                          "(available on github). Details on the compilation " +
                          "process can be found in the PHYEX documentation.")
            raise
        self._pyphyex = pyphyex

        def set_default(nml):
            """
            Set default values for schemes
            """
            if 'PHYEX' not in nml:
                nml['PHYEX'] = {}
            nml['PHYEX']['CMICRO'] = nml['PHYEX'].get('CMICRO', 'ICE3')
            nml['PHYEX']['CSCONV'] = nml['PHYEX'].get('CSCONV', 'EDKF')
            nml['PHYEX']['CTURB'] = nml['PHYEX'].get('CTURB', 'TKEL')
            nml['PHYEX']['RAD'] = nml['PHYEX'].get('RAD', 'ECRAD')
            nml['PHYEX']['SURFACE'] = nml['PHYEX'].get('SURFACE', 'WASP')

        with tempfile.NamedTemporaryFile() as namel:
            nml = f90nml.read(namel.name)
            if self.namel is not None:
                for k, v in self.namel.items():
                    nml[k] = v
            set_default(nml)
            nml.write(namel.name, force=True)
            # First call to really initialize
            try:
                zdzmin = self.grid.get_altitude('FLUX', init_state, self.prognostic_variables).min()
            except (KeyError, Phyex1DError):
                logging.warning('Unable to compute true minimal thickness, use 20m instead')
                zdzmin = 20.
            numnml = 20
            # First call to really initialize
            self._pyphyex.PYINI_PHYEX(
                'PHX1D', 33, namel.name, False, numnml, 0, 1, self._dt, zdzmin,
                nml['PHYEX']['CMICRO'], nml['PHYEX']['CSCONV'],
                nml['PHYEX']['CTURB'],
                LDCHANGEMODEL=True, LDDEFAULTVAL=True, LDREADNAM=True,
                LDCHECK=True, KPRINT=0, LDINIT=True)
            os.remove(f'fort.{numnml}')
            numnml = 21
            # Second call to only write the namelist
            self._pyphyex.PYINI_PHYEX(
                'PHX1D', 33, namel.name, False, numnml, 0, 1, self._dt, zdzmin,
                nml['PHYEX']['CMICRO'], nml['PHYEX']['CSCONV'],
                nml['PHYEX']['CTURB'],
                LDCHANGEMODEL=False, LDDEFAULTVAL=False, LDREADNAM=False,
                LDCHECK=False, KPRINT=1, LDINIT=False)
            self.full_phyex_namel = f90nml.read(f'fort.{numnml}')
            set_default(self.full_phyex_namel)
            os.remove(f'fort.{numnml}')

        # ECRAD SETUP
        if nml['PHYEX']['RAD'] == 'ECRAD' and self.case.radiation == 'on':
            try:
                if self.pyecrad is None:
                    # pyecrad.py should be found by python
                    import pyecrad  # pylint: disable=import-outside-toplevel
                else:
                    # We import pyecrad.py by file path
                    spec = importlib.util.spec_from_file_location('pyecrad', self.pyecrad)
                    pyecrad = importlib.util.module_from_spec(spec)
                    sys.modules['pyecrad'] = pyecrad
                    spec.loader.exec_module(pyecrad)
            except ModuleNotFoundError:
                logging.error("The module pyecrad has not been found.")
                raise

            with tempfile.NamedTemporaryFile() as namel:
                nml = f90nml.read(namel.name)
                if self.namel is not None:
                    for k, v in self.namel.items():
                        nml[k] = v
                set_default(nml)
                nml.write(namel.name, force=True)
                self._pyecrad = pyecrad.Ecrad(namel.name)

    def finalize(self):
        """Close PHYEX resources"""
        super().finalize()
        if self._pyphyex is not None:
            self._pyphyex.close()


class PhysicsArome(PHYEX):
    """
    Implementation of a PHYEX timestep using the AROME way
    This class is not intended to be used directly
    """

    def step(self, state, timestep, timestep_number):
        """Advance the state by one PHYEX timestep"""

        # Save initial state and control
        state0 = {k: v.copy() for (k, v) in state.items()}
        state['tke'] = numpy.maximum(self.full_phyex_namel['NAM_TURBn']['XTKEMIN'],
                                     state['tke'])

        # Preparation: grid and other dimensions
        pressure = self.grid.get_pressure('MASS', state, self.prognostic_variables)
        z_mass = self.grid.get_altitude('MASS', state, self.prognostic_variables)
        z_flux = self.grid.get_altitude('FLUX', state, self.prognostic_variables)
        if self.grid.ascending:
            dzz = numpy.diff(z_flux)
            if dzz[-1] == numpy.inf:
                dzz[-1] = dzz[-2]
        else:
            dzz = -numpy.diff(z_flux)
            if dzz[0] == numpy.inf:
                dzz[0] = dzz[1]
        nijt = 1
        nkt = len(pressure)
        if self.full_phyex_namel['PHYEX']['CMICRO'] == 'ICE3':
            krr = 6
        elif self.full_phyex_namel['PHYEX']['CMICRO'] == 'ICE4':
            krr = 7
        else:
            raise Phyex1DError('Unknown microphysics scheme')
        ksv = 0

        # Preparation: conversions
        exner = (pressure / 1.E5) ** (self.cst.Rd / self.cst.Cpd)
        if 'T' in self.prognostic_variables:
            theta = state['T'] / exner
            temperature = state['T']
        else:
            theta = state['Theta']
            temperature = exner * theta
        if 'qv' in self.prognostic_variables:
            qdm = 1.
            for var in ('qv', 'qc', 'qr', 'qi', 'qs', 'qg', 'qh'):
                qdm -= state[var]
            gas_constant = self.cst.Rd + state['qv'] * (self.cst.Rv - self.cst.Rd)
            for var in ('qc', 'qr', 'qi', 'qs', 'qg', 'qh'):
                gas_constant += - state[var] * self.cst.Rd
            rho = pressure / (gas_constant * temperature)
            rhodref = rho * qdm
            rv = state['qv'] / qdm
            rc = state['qc'] / qdm
            rr = state['qr'] / qdm
            ri = state['qi'] / qdm
            rs = state['qs'] / qdm
            rg = state['qg'] / qdm
            rh = state['qh'] / qdm
        else:
            qdm = 1.
            for var in ('rv', 'rc', 'rr', 'ri', 'rs', 'rg', 'rh'):
                qdm += state[var]
            qdm = 1. / qdm
            rv = state['rv']
            rc = state['rc']
            rr = state['rr']
            ri = state['ri']
            rs = state['rs']
            rg = state['rg']
            rh = state['rh']
            div = 1 + state['rv']
            for var in ('rc', 'rr', 'ri', 'rs', 'rg', 'rh'):
                if var in self.prognostic_variables:
                    div += state[var]
            gas_constant = (self.cst.Rd + state['rv'] * self.cst.Rv) / div
            rho = pressure / (gas_constant * temperature)
            rhodref = pressure / ((self.cst.Rd + rv * self.cst.Rv) * theta * exner)

        # Preparation: derived fields
        rhodj = dzz * rho

        # Preparation: tendencies
        us = state['u'] / timestep
        vs = state['v'] / timestep
        ws = state['w'] / timestep
        tkes = state['tke'] / timestep
        rvs = rv / timestep
        rcs = rc / timestep
        rrs = rr / timestep
        ris = ri / timestep
        rss = rs / timestep
        rgs = rg / timestep
        rhs = rh / timestep
        thetas = theta / timestep
        if 'qv' in self.prognostic_variables:
            dqv = numpy.zeros((nkt, ))
            dqc = numpy.zeros((nkt, ))
            dqr = numpy.zeros((nkt, ))
            dqi = numpy.zeros((nkt, ))
            dqs = numpy.zeros((nkt, ))
            dqg = numpy.zeros((nkt, ))
            dqh = numpy.zeros((nkt, ))
        if 'T' in self.prognostic_variables:
            dtemperature = numpy.zeros((nkt, ))

        # Preparation: misc fields
        sigqsat = numpy.ones((nijt, )) * self.full_phyex_namel['NAM_NEBn']['VSIGQSAT']
        mfconv = numpy.zeros((nkt, ))
        sv = numpy.zeros((ksv, nkt))
        svs = sv / timestep

        ##################################################################
        ##################################################################
        #                       ADJUSTMENT
        ##################################################################
        ##################################################################
        if self.full_phyex_namel['PHYEX']['CMICRO'] == 'NONE':
            pass
        elif self.full_phyex_namel['PHYEX']['CMICRO'] == 'ICE3':
            if 'qv' in self.prognostic_variables:
                rvsin = rvs.copy()
                rcsin = rcs.copy()
                risin = ris.copy()
            if 'T' in self.prognostic_variables:
                thsin = thetas.copy()
            x = numpy.newaxis
            result = self._pyphyex.PYICE_ADJUST(
                nijt, nkt, 1 if self.grid.ascending else -1, 0,
                krr, 'DEPO', timestep,
                sigqsat, rhodj[:, x],
                exner[:, x], rhodref[:, x], state['sigs'][:, x], False,
                mfconv[:, x], pressure[:, x], z_mass[:, x], exner[:, x],
                state['CF_MF'][:, x], state['rc_MF'][:, x],
                state['ri_MF'][:, x], state['WEIGHT_MF_CLOUD'][:, x], rv[:, x],
                rc[:, x], rvs[:, x], rcs[:, x],
                theta[:, x], thetas[:, x],
                True, rr[:, x], ri[:, x],
                ris[:, x], rs[:, x], rg[:, x], 0,
                PRH=rh if krr == 7 else MISSING,
                PHLC_HRC_MF=state['HLC_HRC_MF'][:, x],
                PHLC_HCF_MF=state['HLC_HCF_MF'][:, x],
                PHLI_HRI_MF=state['HLI_HRI_MF'][:, x],
                PHLI_HCF_MF=state['HLI_HCF_MF'][:, x])
            result = [array[:, 0] for array in result]
            (_, _, _, _, _, rvs, rcs, thetas, src, state['CF'], ris, _, _, _, _,
             hlc_hrc, hlc_hcf, hli_hri, hli_hcf) = result
            theta = thetas * timestep
            temperature = theta * exner
            rv = rvs * timestep
            rc = rcs * timestep
            ri = ris * timestep

            if 'qv' in self.prognostic_variables:
                qv = rv * qdm
                qc = rc * qdm
                qr = rr * qdm
                qi = ri * qdm
                qs = rs * qdm
                qg = rg * qdm
                qh = rh * qdm
                gas_constant = self.cst.Rd + state['qv'] * (self.cst.Rv - self.cst.Rd)
                for var in ('qc', 'qr', 'qi', 'qs', 'qg', 'qh'):
                    if var in self.prognostic_variables:
                        gas_constant += - state[var] * self.cst.Rd
                rho = pressure / (gas_constant * temperature)
                rhodref = rho * qdm
                dqv += (rvs - rvsin) * qdm
                dqc += (rcs - rcsin) * qdm
                dqi += (ris - risin) * qdm
            else:
                div = 1. + state['rv']
                for var in ('rc', 'rr', 'ri', 'rs', 'rg', 'rh'):
                    if var in self.prognostic_variables:
                        div += state[var]
                gas_constant = (self.cst.Rd + state['rv'] * self.cst.Rv) / div
                rho = pressure / (gas_constant * temperature)
                rhodref = pressure / ((self.cst.Rd + rv * self.cst.Rv) * theta * exner)

            if 'T' in self.prognostic_variables:
                dtemperature += (thetas - thsin) * exner

            # Update state to be able to recompute pressure and altitude
            if 'qv' in self.prognostic_variables:
                state['qv'] = state0['qv'] + dqv * timestep
                state['qc'] = state0['qc'] + dqc * timestep
                state['qr'] = state0['qr'] + dqr * timestep
                state['qi'] = state0['qi'] + dqi * timestep
                state['qs'] = state0['qs'] + dqs * timestep
                state['qg'] = state0['qg'] + dqg * timestep
                state['qh'] = state0['qh'] + dqh * timestep
            else:
                state['rv'] = rvs * timestep
                state['rc'] = rcs * timestep
                state['rr'] = rrs * timestep
                state['ri'] = ris * timestep
                state['rs'] = rss * timestep
                state['rg'] = rgs * timestep
                state['rh'] = rhs * timestep
            if 'T' in self.prognostic_variables:
                state['T'] = state0['T'] + dtemperature * timestep
            else:
                state['Theta'] = thetas * timestep
            pressure = self.grid.get_pressure('MASS', state, self.prognostic_variables)
            z_mass = self.grid.get_altitude('MASS', state, self.prognostic_variables)
            z_flux = self.grid.get_altitude('FLUX', state, self.prognostic_variables)
            if self.grid.ascending:
                dzz = numpy.diff(z_flux)
                if dzz[-1] == numpy.inf:
                    dzz[-1] = dzz[-2]
            else:
                dzz = -numpy.diff(z_flux)
                if dzz[0] == numpy.inf:
                    dzz[0] = dzz[1]
        else:
            raise Phyex1DError('Wrong CMICRO scheme choice for adjustment')

        ##################################################################
        ##################################################################
        #                       RADIATION
        ##################################################################
        ##################################################################
        if self.case.radiation == 'on':
            if self.full_phyex_namel['PHYEX']['RAD'] == 'NONE':
                pass
            elif self.full_phyex_namel['PHYEX']['RAD'] == 'ECRAD':
                order = -1 if self.grid.ascending else 1
                pressure_flux = self.grid.get_pressure('FLUX', state, self.prognostic_variables)
                temperature_flux = numpy.ndarray(pressure_flux.shape)
                temperature_flux[1:-1] = .5 * (temperature[1:] + temperature[:-1])
                tskin = state['Ts' if self.case.surface_type == 'ocean' else 'Tskin']
                if self.grid.ascending:
                    temperature_flux[0] = tskin
                    temperature_flux[-1] = temperature[-1]
                else:
                    temperature_flux[-1] = tskin
                    temperature_flux[0] = temperature[0]

                x = numpy.newaxis
                solar_irradiance = state.get('SolarIrradiance', 1366.)
                spectral_solar_cycle_multiplier = 0.
                if 'sza' in state:
                    cos_solar_zenith_angle = numpy.ones((nijt, )) * state['sza']
                else:
                    obs = ephem.Observer()
                    obs.lat = str(state['lat'])
                    obs.long = str(state['lon'])
                    obs.date = datetime.datetime.fromisoformat(self.case.start_date) + \
                        datetime.timedelta(seconds=timestep * timestep_number)
                    sun = ephem.Sun(obs)
                    sun.compute(obs)
                    cos_solar_zenith_angle = numpy.ones((nijt, )) * \
                        numpy.cos(numpy.pi / 2. - sun.alt)
                fractional_std = numpy.ones((nkt, )) * 1.  # could be computed by the cloud scheme
                q_liquid = rc * qdm
                # Simple distinction between land (10um) and ocean (13um) by Zhang and Rossow
                re_liquid = (10. if self.case.surface_type == 'ocean' else 13.)
                re_liquid = numpy.ones((nkt, )) * re_liquid * 1.E-6  # convert in m
                q_ice = ri * qdm
                # Liou and Ou
                temperature_c = numpy.minimum(temperature - self.cst.Tt, -0.1)
                re_ice = 326.3 + temperature_c * \
                    (12.42 + temperature_c * (0.197 + temperature_c * 0.0012))
                re_ice = numpy.minimum(numpy.maximum(re_ice, 40), 130) * 1.E-6
                overlap_param = numpy.ones((nkt - 1, ))  # 1 for max
                tskin = numpy.ones((nijt, )) * tskin
                sw_albedo = numpy.ones((1, )) * state.get('albedo', 1.)  # albedo diffuse
                sw_albedo_direct = sw_albedo.copy()  # albedo direct
                lw_emissivity = numpy.ones((1, )) * state.get('emissivity', 1.)  # emissivity
                q = rv * qdm

                result = self._pyecrad.exec(
                    pressure_flux[::order, x], temperature_flux[::order, x],
                    solar_irradiance, spectral_solar_cycle_multiplier, cos_solar_zenith_angle,
                    state['CF'][::order, x], fractional_std[::order, x], q_liquid[::order, x],
                    re_liquid[::order, x], q_ice[::order, x], re_ice[::order, x],
                    overlap_param[::order, x],
                    tskin, sw_albedo[:, x], lw_emissivity[:, x],
                    sw_albedo_direct=sw_albedo_direct[:, x],
                    q_unit=self._pyecrad.IMassMixingRatio, q=q[::order, x])
                result = {k: v[::order, 0] if len(v.shape) == 2 else v
                          for (k, v) in result.items()}
                for var in ('lw_up', 'lw_dn', 'sw_up', 'sw_dn'):
                    state[var] = result[var]
                cp = qdm * (self.cst.Cpd +
                            self.cst.Cpv * rv +
                            self.cst.Cl * (rc + rr) +
                            self.cst.Ci * (ri + rs + rg + rh))
                dt_rad_lw = -numpy.diff((state['lw_dn'] - state['lw_up'])[::order])[::order] / \
                    (rho * cp * dzz)
                dt_rad_sw = -numpy.diff((state['sw_dn'] - state['sw_up'])[::order])[::order] / \
                    (rho * cp * dzz)
                state['dt_rad_lw'] = dt_rad_lw
                state['dt_rad_sw'] = dt_rad_sw
                if 'T' in self.prognostic_variables:
                    dtemperature += dt_rad_lw + dt_rad_sw
                else:
                    thetas += (dt_rad_lw + dt_rad_sw) * exner
            else:
                raise Phyex1DError('Wrong RAD scheme choice')

        ##################################################################
        ##################################################################
        #                       SURFACE
        ##################################################################
        ##################################################################
        klevgrd = 0 if self.grid.ascending else -1
        if self.case.surface_forcing_temp in ('none', 'ts') or \
           self.case.surface_forcing_moisture in ('none', 'beta', 'mrsos') or \
           self.case.surface_forcing_wind in ('none', ):
            need_scheme = True
        else:
            need_scheme = False

        def ustar2fluxes(ustar, u, v):
            """
            Convert ustar into momentum fluxes
            :param ustar: ustar value
            :param u, v: wind components
            :return: w'u', w'v'
            """
            alpha = numpy.arctan2(v, u)
            return - ustar**2 * numpy.cos(alpha), - ustar**2 * numpy.sin(alpha)

        if 'Ts' in state:
            surface_temperature = state['Ts']
        else:
            if self.grid.ascending:
                surface_temperature = temperature[0]
            else:
                surface_temperature = temperature[-1]
        if surface_temperature > self.cst.Tt:
            latent_heat = self.cst.LvTt - (self.cst.Cpv - self.cst.Cl) * \
                (surface_temperature - self.cst.Tt)
        else:
            latent_heat = self.cst.LsTt - (self.cst.Cpv - self.cst.Ci) * \
                (surface_temperature - self.cst.Tt)

        # Surface schemes output must be: sshf (W m-2) and swf (kg m-2 s-1)
        sshf_scheme = None
        swf_scheme = None
        sfu_scheme = None
        sfv_scheme = None
        if self.full_phyex_namel['PHYEX']['SURFACE'] == 'NONE':
            sshf_scheme = 0.
            swf_scheme = 0.
            sfu_scheme = 0.
            sfv_scheme = 0.
        elif not need_scheme:
            pass
        elif self.full_phyex_namel['PHYEX']['SURFACE'] == 'WASP':
            if self.case.surface_forcing_temp == 'none':
                raise Phyex1DError('surface_forcing_temp == none not implemented by WASP')
            if self.case.surface_forcing_moisture in ('beta', 'mrsos'):
                raise Phyex1DError(
                    'surface_forcing_moisture == beta or mrsos not implemented by WASP')
            if self.case.surface_forcing_wind != 'none':
                raise Phyex1DError('surface_forcing_wind != none not implemented by WASP')
            windgrd = numpy.sqrt(state['u'][klevgrd]**2 + state['v'][klevgrd]**2)
            exner_surf = (state['Ps'] / 1.E5) ** (self.cst.Rd / self.cst.Cpd)
            if 'qv' in self.prognostic_variables:
                qv_grd = state['rv'][klevgrd] * qdm[klevgrd]
            else:
                qv_grd = state['qv'][klevgrd]
            # TODO, filled precip in precipitation rate (kg/s/m2)
            precip = 0.
            wave_height = 0.
            wave_peak_period = 0.
            result = wasp.WASP_FLUX(
                state['T'][klevgrd], qv_grd, exner[klevgrd], rho[klevgrd], windgrd,
                z_mass[klevgrd], z_mass[klevgrd],
                state['Ts' if self.case.surface_type == 'ocean' else 'Tskin'],
                exner_surf, state['Ps'], precip, wave_height, wave_peak_period)
            sshf_scheme, swf_scheme, ustar = result
            sfu_scheme, sfv_scheme = ustar2fluxes(ustar, state['u'][klevgrd],
                                                   state['v'][klevgrd])
        else:
            raise Phyex1DError('Wrong SURFACE scheme choice')

        if self.case.surface_forcing_temp == 'none':
            if sshf_scheme is None:
                raise Phyex1DError("Surface scheme didn't compute sshf")
        elif self.case.surface_forcing_temp == 'kinematic':
            raise NotImplementedError('surface_forcing_temp=kinematic')
        elif self.case.surface_forcing_temp == 'surface_flux':
            sshf = state['sshf']
        elif self.case.surface_forcing_temp == 'ts':
            sshf = sshf_scheme
        else:
            raise Phyex1DError('Wrong surface_forcing_temp option')

        if self.case.surface_forcing_moisture == 'none':
            if swf_scheme is None:
                raise Phyex1DError("Surface scheme didn't compute swf")
            swf = swf_scheme
        elif self.case.surface_forcing_moisture == 'kinematic':
            raise NotImplementedError('surface_forcing_moisture=kinematic')
        elif self.case.surface_forcing_moisture == 'surface_flux':
            # W m-2 --> kg m-2 s-1
            swf = state['slhf'] / latent_heat
        elif self.case.surface_forcing_moisture == 'beta':
            raise NotImplementedError('surface_forcing_moisture=beta')
        elif self.case.surface_forcing_moisture == 'mrsos':
            raise NotImplementedError('surface_forcing_moisture=mrsos')
        else:
            raise Phyex1DError('Wrong surface_forcing_moisture option')

        if self.case.surface_forcing_wind == 'none':
            sfu = sfu_scheme
            sfv = sfv_scheme
        elif self.case.surface_forcing_wind == 'z0':
            ustar = numpy.sqrt(state['u'][klevgrd]**2 + state['v'][klevgrd]**2) * \
                self.cst.Karman / numpy.log(z_mass[klevgrd] / state['z0'])
            sfu, sfv = ustar2fluxes(ustar, state['u'][klevgrd], state['v'][klevgrd])
        elif self.case.surface_forcing_wind == 'ustar':
            sfu, sfv = ustar2fluxes(state['ustar'], state['u'][klevgrd], state['v'][klevgrd])
        else:
            raise Phyex1DError('Wrong surface_forcing_wind option')
        sfsv = numpy.zeros((ksv, ))

        # The turbulence scheme waits for w'r' and w'th'
        # for the water flux, we want w'r' (division by rhod) and not w'q' (division by rho)
        sfrv = swf / rhodref[0 if self.grid.ascending else -1]  # kg m-2 s-1 --> kg/kg m s-1
        cp = qdm * (self.cst.Cpd + self.cst.Cpv * rv + self.cst.Cl * (rc + rr) +
                    self.cst.Ci * (ri + rs + rg + rh))
        sfth = sshf / (self.cst.Cpd * rho[0 if self.grid.ascending else -1])  # W m-2 --> K m s-1

        state['sfrv'] = sfrv  # kg/kg m s-1 (w'r' mixing ratio)
        state['swf'] = swf  # kg m-2 s-1
        state['slhf'] = swf * latent_heat  # W m-2 using Lv or Ls depending on temperature °C sign
        state['sfth'] = sfth  # K m s-1 (w'th' theta)
        state['sshf'] = sshf  # W m-2

        ##################################################################
        ##################################################################
        #                       SHALLOW CONVECTION
        ##################################################################
        ##################################################################
        if self.full_phyex_namel['PHYEX']['CSCONV'] == 'NONE':
            flxzthvmf = flxzumf = flxzvmf = numpy.zeros((nkt, ))
        elif self.full_phyex_namel['PHYEX']['CSCONV'] == 'EDKF':
            # Updraft properties
            pthl_up = numpy.zeros((nkt, ))
            prt_up = numpy.zeros((nkt, ))
            prv_up = numpy.zeros((nkt, ))
            prc_up = numpy.zeros((nkt, ))
            pri_up = numpy.zeros((nkt, ))
            pu_up = numpy.zeros((nkt, ))
            pv_up = numpy.zeros((nkt, ))
            ptke_up = numpy.zeros((nkt, ))
            pthv_up = numpy.zeros((nkt, ))
            pw_up = numpy.zeros((nkt, ))
            pfrac_up = numpy.zeros((nkt, ))
            pemf = numpy.zeros((nkt, ))

            # Other arrays
            rm = numpy.array([rv, rc, rr, ri, rs, rg, rh])
            if self.grid.ascending:
                z_flux_trunc = z_flux[:-1]
            else:
                z_flux_trunc = z_flux[1:]

            x = numpy.newaxis
            result = self._pyphyex.PYSHALLOW_MF(
                nijt, nkt, 1 if self.grid.ascending else -1,
                0, krr, 2, 3, 0, False, 0, 0,
                timestep, dzz[:, x], z_flux_trunc[:, x], rhodj[:, x],
                rhodref[:, x], pressure[:, x], exner[:, x],
                numpy.ones((nijt, )) * sfth, numpy.ones((nijt, )) * sfrv, theta[:, x],
                rm[:krr, :, x], state['u'][:, x], state['v'][:, x], state['tke'][:, x],
                sv[:, :, x], pthl_up[:, x], prt_up[:, x], prv_up[:, x], prc_up[:, x],
                pri_up[:, x], pu_up[:, x], pv_up[:, x], ptke_up[:, x], pthv_up[:, x],
                pw_up[:, x], pfrac_up[:, x], pemf[:, x],
                self.dx, self.dy, PRSVS=svs[:, :, x], KBUDGETS=0)
            result = [array[..., 0] for array in result]
            (du_mf, dv_mf, dtke_mf, dthl_mf, drt_mf, dsv_mf, sigs_mf,
             state['rc_MF'], state['ri_MF'], state['CF_MF'], state['HLC_HRC_MF'],
             state['HLC_HCF_MF'], state['HLI_HRI_MF'], state['HLI_HCF_MF'],
             state['WEIGHT_MF_CLOUD'],
             flxzthvmf, flxzthmf, flxzrmf, flxzumf, flxzvmf, flxztkemf) = result[:21]
            us += du_mf
            vs += dv_mf
            tkes += dtke_mf
            svs += dsv_mf
            if 'Theta' in self.prognostic_variables:
                thetas += dthl_mf
            else:
                dtemperature += dthl_mf * exner
            if 'rv' in self.prognostic_variables:
                rvs += drt_mf
            else:
                dqv += drt_mf * qdm
        else:
            raise Phyex1DError('Wrong CSCONV scheme choice')

        ##################################################################
        ##################################################################
        #                       TURBULENCE
        ##################################################################
        ##################################################################
        if self.full_phyex_namel['PHYEX']['CTURB'] == 'NONE':
            pass
        elif self.full_phyex_namel['PHYEX']['CTURB'] == 'TKEL':
            hlbcx = hlbcy = numpy.array(['cycl', 'cycl'], dtype=('S', 4))
            kgradientsleo, kgradientsgog, khalo, ksplit = 0, 0, 1, 1
            ocloudmodiflm = False
            ksv_lgbeg, ksv_lgend = 0, 0
            ksv_lima_nr, ksv_lima_ns, ksv_lima_ng, ksv_lima_nh = 0, 0, 0, 0
            o2d, onomixlg, oflat, ocouples = False, False, False, False
            oblowsnow, oibm, oflyer, ocompute_src = False, False, True, True
            rsnow = 1.
            oocean, odeepoc, odiag_in_run = False, False, False
            hturblen_cl, helec = 'DELT', 'NONE'
            dxx = dyy = dzx = dzy = numpy.ndarray((nkt, ))
            dircosxw = dircosyw = dircoszw = cosslope = 1.
            sinslope = 0.
            hgradleo = numpy.ndarray((kgradientsleo, nkt))
            hgradgog = numpy.ndarray((kgradientsgog, nkt))
            lengthm, lengthh = numpy.zeros((nkt, )), numpy.zeros((nkt, ))
            mfmoist = numpy.zeros((nkt, ))
            cei = numpy.zeros((nkt, ))
            cei_min, cei_max = 0.001E-06, 0.01E-06
            coef_ampl_sat = 5.
            kbudgets = 12
            bl_depth, sbl_depth = 0., 0.
            prm = numpy.array([rv, rc, rr, ri, rs, rg, rh])
            prs = numpy.array([rvs, rcs, rrs, ris, rss, rgs, rhs])
            thvref = theta * (1 + rv * self.cst.Rv / self.cst.Rd) / (1 + prm.sum(axis=0))
            nvext_turb = 1
            krrl = 2
            krri = 3 if krr == 6 else 4
            psea_ucu = psea_vcu = 0.

            def x(array, reverse=False):
                """Extend vertical dimension and add horizontal one"""
                if not reverse:
                    try:
                        len(array)
                        if len(array.shape) == 1:
                            new = numpy.ndarray((array.shape[0] + 2, ))
                            new[1:-1] = array
                            new[0] = new[1]
                            new[-1] = new[-2]
                            return new[:, numpy.newaxis]
                        if len(array.shape) == 2:
                            new = numpy.ndarray((array.shape[0], array.shape[1] + 2))
                            new[:, 1:-1] = array
                            new[:, 0] = new[:, 1]
                            new[:, -1] = new[:, -2]
                            return new[:, :, numpy.newaxis]
                    except TypeError:
                        return numpy.ones((nijt, )) * array
                    raise ValueError(type(array), str(array))
                else:
                    if len(array.shape) == 3:
                        return array[:, 1:-1, 0]
                    if len(array.shape) == 2:
                        return array[1:-1, 0]

            us = us * rhodj
            vs = vs * rhodj
            ws = ws * rhodj
            thetas = thetas * rhodj
            tkes = tkes * rhodj
            prs = prs * rhodj[numpy.newaxis, :]
            svs = svs * rhodj
            if self.grid.ascending:
                z_flux_trunc = z_flux[:-1]
            else:
                z_flux_trunc = z_flux[1:]

            result = self._pyphyex.PYTURB(
                nijt, nkt + 2 * nvext_turb, 1 if self.grid.ascending else -1, nvext_turb,
                krr, krrl, krri, hlbcx, hlbcy, kgradientsleo, kgradientsgog,
                khalo, ksplit, ocloudmodiflm, ksv, ksv_lgbeg, ksv_lgend,
                ksv_lima_nr, ksv_lima_ns, ksv_lima_ng, ksv_lima_nh,
                o2d, onomixlg, oflat, ocouples, oblowsnow, oibm, oflyer,
                ocompute_src, rsnow, oocean, odeepoc, odiag_in_run,
                hturblen_cl, self.full_phyex_namel['PHYEX']['CMICRO'], helec, timestep, 999,
                x(dxx), x(dyy), x(dzz), x(dzx), x(dzy), x(z_flux_trunc),
                x(dircosxw), x(dircosyw), x(dircoszw), x(cosslope), x(sinslope),
                x(rhodj), x(thvref), x(hgradleo), x(hgradgog), x(state['Zs']),
                x(sfth), x(sfrv), sfsv[:, numpy.newaxis], x(sfu), x(sfv),
                x(psea_ucu), x(psea_vcu),
                x(pressure), x(state['u']), x(state['v']), x(state['w']),
                x(state['tke']), x(sv), x(src), x(lengthm), x(lengthh), x(mfmoist),
                x(bl_depth), x(sbl_depth), x(cei), cei_min, cei_max, coef_ampl_sat,
                x(theta), x(prm[:krr, ...]), x(us), x(vs), x(ws), x(thetas),
                x(prs[:krr, ...]), x(svs), x(tkes),
                x(flxzthvmf), x(flxzumf), x(flxzvmf),
                kbudgets, missingOUT=['PEDR', 'PLEM', 'PDPMF', 'PTPMF',
                'PTR', 'PDISS', 'PIBM_XMUT', 'PCURRENT_TKE_DISS'])

            result = [x(array, reverse=True) for array in result]
            (_, _, _, _, us, vs, ws, thetas, prs, svs, tkes, sigs_turb, _,
             _, _, _, _, _, _, _, _, dthetal_turb, drt_turb, _) = result

            us = us / rhodj
            vs = vs / rhodj
            ws = ws / rhodj
            thetas = thetas / rhodj
            tkes = tkes / rhodj
            dthetal_turb = dthetal_turb / rhodj
            drt_turb = drt_turb / rhodj
            prs = prs / rhodj[numpy.newaxis, :]
            svs = svs / rhodj

            if 'T' in self.prognostic_variables:
                dtemperature += dthetal_turb * exner
            if 'qv' in self.prognostic_variables:
                print('We must take into account tendencies on rc and ri. Here and in arome')
                dqv += drt_turb * qdm
            else:
                rvs = prs[0]
                rcs = prs[1]
                rrs = prs[2]
                ris = prs[3]
                rss = prs[4]
                rgs = prs[5]
                if krr == 7:
                    rhs = prs[6]
            state['sigs'] = numpy.sqrt(sigs_mf**2 + sigs_turb**2)
        else:
            raise Phyex1DError('Wrong CTURB scheme choice')

        ##################################################################
        ##################################################################
        #                       MICROPHYSICS
        ##################################################################
        ##################################################################
        if self.full_phyex_namel['PHYEX']['CMICRO'] == 'NONE':
            pass
        elif self.full_phyex_namel['PHYEX']['CMICRO'] == 'ICE3':
            if 'qv' in self.prognostic_variables:
                rvsin = rvs.copy()
                rcsin = rcs.copy()
                rrsin = rrs.copy()
                risin = ris.copy()
                rssin = rss.copy()
                rgsin = rgs.copy()
                rhsin = rhs.copy()
            if 'T' in self.prognostic_variables:
                thsin = thetas.copy()

            x = numpy.newaxis

            if krr == 7:
                prm = numpy.array([rv[:, x], rc[:, x], rr[:, x], ri[:, x],
                                   rs[:, x], rg[:, x], rh[:, x]])
                prs = numpy.array([rvs[:, x], rcs[:, x], rrs[:, x], ris[:, x],
                                   rss[:, x], rgs[:, x], rhs[:, x]])
            else:
                prm = numpy.array([rv[:, x], rc[:, x], rr[:, x], ri[:, x],
                                   rs[:, x], rg[:, x]])
                prs = numpy.array([rvs[:, x], rcs[:, x], rrs[:, x], ris[:, x],
                                   rss[:, x], rgs[:, x]])

            dummy = ris.copy() * 0.

            result = self._pyphyex.PYRAIN_ICE(
                nijt, nkt, 1 if self.grid.ascending else -1, 0, False, False,
                MISSING, timestep, krr, exner[:, x], dzz[:, x], rhodj[:, x],
                rhodref[:, x], exner[:, x], pressure[:, x], state['ni'][:, x],
                state['CF'][:, x],
                dummy[:, x], dummy[:, x], dummy[:, x], dummy[:, x],
                hlc_hrc[:, x], hlc_hcf[:, x], hli_hri[:, x], hli_hcf[:, x],
                theta[:, x], prm, thetas[:, x], prs, state['sigs'][:, x],
                0, MISSING, MISSING, MISSING, MISSING, MISSING, MISSING,
                MISSING, MISSING, MISSING, MISSING, MISSING, MISSING, MISSING, MISSING,
                MISSING, MISSING, MISSING, MISSING, MISSING,
                MISSING, MISSING)

            (state['ni'], hlc_hrc, hlc_hcf, hli_hri, hli_hcf,
             thetas, prs, _, _, _, _, _, _, _, _, _, ) = result
            state['ni'] = state['ni'][:, 0]
            hlc_hrc = hlc_hrc[:, 0]
            hlc_hcf = hlc_hcf[:, 0]
            hli_hri = hli_hri[:, 0]
            hli_hcf = hli_hcf[:, 0]
            thetas = thetas[:, 0]

            rvs = prs[0, :, 0]
            rcs = prs[1, :, 0]
            rrs = prs[2, :, 0]
            ris = prs[3, :, 0]
            rss = prs[4, :, 0]
            rgs = prs[5, :, 0]
            if krr == 7:
                rhs = prs[6, :, 0] * timestep

            if 'qv' in self.prognostic_variables:
                dqv += (prs[0, :, 0] - rvsin) * qdm
                dqc += (prs[1, :, 0] - rcsin) * qdm
                dqr += (prs[2, :, 0] - rrsin) * qdm
                dqi += (prs[3, :, 0] - risin) * qdm
                dqs += (prs[4, :, 0] - rssin) * qdm
                dqg += (prs[5, :, 0] - rgsin) * qdm
                if krr == 7:
                    dqh += (prs[6, :, 0] - rhsin) * qdm

            if 'T' in self.prognostic_variables:
                dtemperature += (thetas - thsin) * exner
        else:
            raise Phyex1DError('Wrong CMICRO scheme choice for microphysics')

        ##################################################################
        ##################################################################
        #                       NEW VALUES
        ##################################################################
        ##################################################################
        if 'qv' in self.prognostic_variables:
            state['qv'] = state0['qv'] + dqv * timestep
            state['qc'] = state0['qc'] + dqc * timestep
            state['qr'] = state0['qr'] + dqr * timestep
            state['qi'] = state0['qi'] + dqi * timestep
            state['qs'] = state0['qs'] + dqs * timestep
            state['qg'] = state0['qg'] + dqg * timestep
            state['qh'] = state0['qh'] + dqh * timestep
        else:
            state['rv'] = rvs * timestep
            state['rc'] = rcs * timestep
            state['rr'] = rrs * timestep
            state['ri'] = ris * timestep
            state['rs'] = rss * timestep
            state['rg'] = rgs * timestep
            state['rh'] = rhs * timestep
        if 'T' in self.prognostic_variables:
            state['T'] = state0['T'] + dtemperature * timestep
        else:
            state['Theta'] = thetas * timestep
        state['u'] = us * timestep
        state['v'] = vs * timestep
        state['tke'] = tkes * timestep
        return state


class PhysicsAromeTQ(PhysicsArome):
    """
    AROME physics with T and q as prognostic variables
    """

    def __init__(self, dt, method, name, tag, case, grid, pyphyex=None, pyecrad=None,
                 namel=None, dx=0, dy=0, attrs=None):
        """
        :param dt: timestep (s)
        :param method: 'step-by-step' or 'one-step'
        :param name: full name of the execution
        :param tag: tag to identify the execution
        :param case: Case instance describing the case
        :param grid: Grid object instance
        :param phyex: path to the pyphyex.py file generated by the PHYEX package
        :param pyecrad: path to the pyecrad.py file
        :param namel: namelist to use: as a file name or a dictionnary
        :param dx, dy: mesh size
        :param attrs: dictionnary holding caes attributes to override
        """
        super().__init__(dt, method, name, tag, case, grid,
                         ['T', 'qv', 'qc', 'qi', 'qr', 'qs', 'qg', 'qh', 'u', 'v', 'w', 'tke'],
                         pyphyex, pyecrad, namel, dx, dy, attrs)


class PhysicsAromeThetaR(PhysicsArome):
    """
    AROME physics with Theta and r as prognostic variables
    """

    def __init__(self, dt, method, name, tag, case, grid, pyphyex=None, pyecrad=None,
                 namel=None, dx=0, dy=0, attrs=None):
        """
        :param dt: timestep (s)
        :param method: 'step-by-step' or 'one-step'
        :param name: full name of the execution
        :param tag: tag to identify the execution
        :param case: Case instance describing the case
        :param grid: Grid object instance
        :param phyex: path to the pyphyex.py file generated by the PHYEX package
        :param pyecrad: path to the pyecrad.py file
        :param namel: namelist to use: as a file name or a dictionnary
        :param dx, dy: mesh size
        :param attrs: dictionnary holding caes attributes to override
        """
        super().__init__(dt, method, name, tag, case, grid,
                         ['Theta', 'rv', 'rc', 'ri', 'rr', 'rs', 'rg', 'rh', 'u', 'v', 'w', 'tke'],
                         pyphyex, pyecrad, namel, dx, dy, attrs)


class PhysicsForcingTQ(PHYEX):
    """
    Only the forcing with T and q as prognostic variables
    """

    def __init__(self, dt, method, name, tag, case, grid, pyphyex=None, pyecrad=None,
                 namel=None, dx=0, dy=0, attrs=None):
        """
        :param dt: timestep (s)
        :param method: 'step-by-step' or 'one-step'
        :param name: full name of the execution
        :param tag: tag to identify the execution
        :param case: Case instance describing the case
        :param grid: Grid object instance
        :param phyex: path to the pyphyex.py file generated by the PHYEX package
        :param pyecrad: path to the pyecrad.py file
        :param namel: namelist to use: as a file name or a dictionnary
        :param dx, dy: mesh size
        :param attrs: dictionnary holding caes attributes to override
        """
        super().__init__(dt, method, name, tag, case, grid,
                         ['T', 'qv', 'qc', 'qi', 'qr', 'qs', 'qg', 'qh', 'u', 'v', 'w', 'tke'],
                         pyphyex, pyecrad, namel, dx, dy, attrs)


class PhysicsForcingThetaR(PHYEX):
    """
    Only the forcing with Theta and r as prognostic variables
    """

    def __init__(self, dt, method, name, tag, case, grid, pyphyex=None, pyecrad=None,
                 namel=None, dx=0, dy=0, attrs=None):
        """
        :param dt: timestep (s)
        :param method: 'step-by-step' or 'one-step'
        :param name: full name of the execution
        :param tag: tag to identify the execution
        :param case: Case instance describing the case
        :param grid: Grid object instance
        :param phyex: path to the pyphyex.py file generated by the PHYEX package
        :param pyecrad: path to the pyecrad.py file
        :param namel: namelist to use: as a file name or a dictionnary
        :param dx, dy: mesh size
        :param attrs: dictionnary holding caes attributes to override
        """
        super().__init__(dt, method, name, tag, case, grid,
                         ['Theta', 'rv', 'rc', 'ri', 'rr', 'rs', 'rg', 'rh', 'u', 'v', 'w', 'tke'],
                         pyphyex, pyecrad, namel, dx, dy, attrs)
