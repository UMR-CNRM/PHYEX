"""
Minimal base class for physics schemes
"""

import logging

import numpy
from pppy import PPPY
from scipy.interpolate import RegularGridInterpolator

from . import Cst
from . import Phyex1DError


class PhysicsBase(PPPY):
    """
    Base class for 1D model physics schemes

    This class is not intended to be used directly.
    """

    def __init__(self, dt, method, name, tag, case, grid, prognostic_variables):
        super().__init__(dt, method, name, tag, inputfile=case.description,
                         grid=grid.filename)
        self.case = case
        self.grid = grid
        self.prognostic_variables = prognostic_variables
        self.cst = Cst()

    def setup(self, init_state, duration):
        """Set up the scheme. Override in subclasses."""

    def finalize(self):
        """Tear down the scheme. Override in subclasses."""

    def execute(self, previous_state, timestep, timestep_number):
        """
        Advance the state by one timestep

        Parameters
        ----------
        previous_state : dict
            Dictionary holding all the variables.
        timestep : float
            Timestep (s)
        timestep_number : int
            Current time step number

        Returns
        -------
        dict
            State after time integration
        """
        new_state = self.step(previous_state, timestep, timestep_number)
        new_state = self.forcing(new_state, timestep, timestep_number)
        self.add_vertical_coordinate(new_state)
        self.add_conversions(new_state)
        return new_state

    def step(self, state, timestep, timestep_number):
        """
        Physics step. Override in subclasses.

        Returns the state unchanged by default.
        """
        return state

    def add_vertical_coordinate(self, state):
        """
        Add vertical coordinate to the state

        Parameters
        ----------
        state : dict
            Current state
        """
        state['Z_mass'] = self.grid.get_altitude('MASS', state, self.prognostic_variables)
        state['P'] = self.grid.get_pressure('MASS', state, self.prognostic_variables)
        state['Z_flux'] = self.grid.get_altitude('FLUX', state, self.prognostic_variables)
        top = -1 if self.grid.ascending else 0
        if numpy.all(state['Z_flux'][top] == numpy.inf):
            state['Z_flux'][top] = state['Z_flux'][-2 if self.grid.ascending else 1]

    def add_conversions(self, state):
        """
        Add derived variables to the state

        Parameters
        ----------
        state : dict
            Current state
        """
        pressure_mass = self.grid.get_pressure('MASS', state, self.prognostic_variables)
        exner = (pressure_mass / 1.E5) ** (self.cst.Rd / self.cst.Cpd)
        if 'T' in self.prognostic_variables:
            state['Theta'] = state['T'] / exner
        elif 'Theta' in self.prognostic_variables:
            state['T'] = state['Theta'] * exner

        if 'qv' in self.prognostic_variables:
            qt = state['qv'].copy()
            for var in ('qc', 'qi', 'qr', 'qs', 'qg', 'qh'):
                qt += state[var]
            state['qt'] = qt
            for v in ('v', 'c', 'i', 'r', 's', 'g', 'h'):
                state['r' + v] = state['q' + v] / (1. - qt)
        elif 'rv' in self.prognostic_variables:
            rt = state['rv'].copy()
            for var in ('rc', 'ri', 'rr', 'rs', 'rg', 'rh'):
                rt += state[var]
            state['rt'] = rt
            for v in ('v', 'c', 'i', 'r', 's', 'g', 'h'):
                state['q' + v] = state['r' + v] / (1. + rt)

    def build_init_state(self, init_state):
        """
        Build the initial state

        Parameters
        ----------
        init_state : dict
            Initial state dictionary to populate

        Returns
        -------
        dict
            Populated initial state
        """
        init_state = super().build_init_state(init_state)

        with self.case:
            # Surface fields
            for var in ('Zs', 'Ps', 'sshf', 'slhf', 'z0', 'z0h', 'z0q', 'Tskin', 'Ts', 'ustar',
                        'SolarIrradiance', 'albedo', 'emissivity', 'sza', 'lat', 'lon'):
                if self.case.contains(var):
                    init_state[var] = self.case.get(var, 0)
                else:
                    logging.warning('%s not found in the case description!', var)

            # We perform the interpolation on the natural vertical coordinate
            if self.grid.kind in ('H', 'hybridH'):
                input_coord = self.case.get('zh', 0) + self.case.get('orog', 0)
                output_coord = self.grid.get_altitude('MASS', init_state, self.prognostic_variables)
            elif self.grid.kind in ('P', 'hybridP'):
                input_coord = self.case.get('pa', 0)
                output_coord = self.grid.get_pressure('MASS', init_state, self.prognostic_variables)
            else:
                raise Phyex1DError('Wrong grid kind')

            # Interpolation of initial profiles
            for var in self.prognostic_variables:
                if var in ('qc', 'qi', 'qr', 'qs', 'qg', 'qh',
                           'rc', 'ri', 'rr', 'rs', 'rg', 'rh',
                           'w') and not self.case.contains(var):
                    init_state[var] = numpy.zeros(output_coord.shape)
                    logging.warning('%s not found in the case description!', var)
                else:
                    interp = RegularGridInterpolator((input_coord, ), self.case.get(var, 0),
                                                     bounds_error=False, fill_value=None)
                    init_state[var] = interp(output_coord)
            init_state['tke'] = numpy.maximum(self.full_phyex_namel['NAM_TURBn']['XTKEMIN'],
                                              init_state['tke'])

        # Vertical coordinates
        self.add_vertical_coordinate(init_state)

        # Conversions
        self.add_conversions(init_state)

        # Arrays needed to keep memory from one timestep to the other
        for var in ('ni', 'sigs', 'WEIGHT_MF_CLOUD', 'CF_MF', 'rc_MF', 'ri_MF',
                    'HLC_HRC_MF', 'HLC_HCF_MF', 'HLI_HRI_MF', 'HLI_HCF_MF'):
            init_state[var] = numpy.zeros(init_state['P'].shape)

        # Initial value
        for var in ('CF', 'dt_rad_lw', 'dt_rad_sw'):
            init_state[var] = numpy.zeros(init_state['P'].shape)
        for var in ('lw_up', 'lw_dn', 'sw_up', 'sw_dn'):
            init_state[var] = numpy.zeros((init_state['P'].shape[0] + 1, ))
        for var in ('sfrv', 'swf', 'sfth', 'sshf', 'slhf'):
            if var not in init_state:
                init_state[var] = numpy.zeros((1, ))

        return init_state

    def forcing(self, state, timestep, timestep_number):
        """
        Apply large-scale forcings to the state

        Parameters
        ----------
        state : dict
            Dictionary holding all the variables.
        timestep : float
            Timestep (s)
        timestep_number : int
            Current time step number

        Returns
        -------
        dict
            Updated state after forcing
        """
        time = timestep * timestep_number
        # All forcing variable are taken on mass levels
        if self.grid.kind in ('H', 'hybridH'):
            output_coord = self.grid.get_altitude('MASS', state, self.prognostic_variables)
        elif self.grid.kind in ('P', 'hybridP'):
            output_coord = self.grid.get_pressure('MASS', state, self.prognostic_variables)
        else:
            raise Phyex1DError('Wrong grid kind')

        with self.case:

            # Temperature advection and nudging, radiative tendency
            if 'Theta' in self.prognostic_variables:
                if self.case.adv_theta == 1:
                    interp = self.case.get_interpolator('tntheta_adv', self.grid)
                    state['Theta'] += interp(numpy.column_stack(([time] * len(output_coord),
                                                                 output_coord))) * timestep
                if self.case.nudging_theta == 1:
                    raise NotImplementedError('Nudging theta')
                if self.case.radiation == 'tend':
                    interp = self.case.get_interpolator('tntheta_rad', self.grid)
                    state['Theta'] += interp(numpy.column_stack(([time] * len(output_coord),
                                                                 output_coord))) * timestep
            elif 'T' in self.prognostic_variables:
                if self.case.adv_ta == 1:
                    interp = self.case.get_interpolator('tnta_adv', self.grid)
                    state['T'] += interp(numpy.column_stack(([time] * len(output_coord),
                                                             output_coord))) * timestep
                if self.case.nudging_ta == 1:
                    raise NotImplementedError('Nudging ta')
                if self.case.radiation == 'tend':
                    interp = self.case.get_interpolator('tnta_rad', self.grid)
                    state['T'] += interp(numpy.column_stack(([time] * len(output_coord),
                                                             output_coord))) * timestep

            # Hydrometeors advection and nudging
            if 'rv' in self.prognostic_variables:
                if self.case.adv_rv == 1:
                    interp = self.case.get_interpolator('tnrv_adv', self.grid)
                    state['rv'] += interp(numpy.column_stack(([time] * len(output_coord),
                                                              output_coord))) * timestep
                if self.case.nudging_rv == 1:
                    raise NotImplementedError('Nudging rv')
            elif 'qv' in self.prognostic_variables:
                if self.case.adv_qv == 1:
                    interp = self.case.get_interpolator('tnqv_adv', self.grid)
                    state['qv'] += interp(numpy.column_stack(([time] * len(output_coord),
                                                              output_coord))) * timestep
                if self.case.nudging_qv == 1:
                    raise NotImplementedError('Nudging qv')

            # Wind advection and nudging
            if self.case.adv_ua == 1:
                interp = self.case.get_interpolator('tnua_adv', self.grid)
                state['u'] += interp(numpy.column_stack(([time] * len(output_coord),
                                                         output_coord))) * timestep
            if self.case.adv_va == 1:
                interp = self.case.get_interpolator('tnva_adv', self.grid)
                state['v'] += interp(numpy.column_stack(([time] * len(output_coord),
                                                         output_coord))) * timestep
            if self.case.nudging_ua == 1:
                raise NotImplementedError('Nudging u')
            if self.case.nudging_va == 1:
                raise NotImplementedError('Nudging v')
            if self.case.forc_wa == 1:
                interp = self.case.get_interpolator('wa', self.grid)
                w = interp(numpy.column_stack(([time] * len(output_coord),
                                               output_coord)))
                z_mass = self.grid.get_altitude('MASS', state, self.prognostic_variables)
                for var in [var for var in self.prognostic_variables if var != 'w']:
                    grad = numpy.gradient(state[var], z_mass)
                    state[var] += - grad * w * timestep
                    if var in ('qv', 'qc', 'qr', 'qi', 'qs', 'qg', 'qh',
                               'rv', 'rc', 'rr', 'ri', 'rs', 'rg', 'rh',
                               'tke'):
                        state[var] = numpy.maximum(state[var], 0.)
            if self.case.forc_geo == 1:
                interp_ug = self.case.get_interpolator('ug', self.grid)
                interp_vg = self.case.get_interpolator('vg', self.grid)
                f = 2. * self.cst.omega * numpy.sin(numpy.radians(state['lat']))
                state['u'] += +f * (state['v'] - interp_vg(
                    numpy.column_stack(([time] * len(output_coord), output_coord)))) * timestep
                state['v'] += -f * (state['u'] - interp_ug(
                    numpy.column_stack(([time] * len(output_coord), output_coord)))) * timestep

            # Surface
            for surf_var in ('Ts', 'Tskin', 'Ps'):
                if surf_var in state:
                    interp = self.case.get_interpolator(
                        {'Ts': 'ts_forc', 'Tskin': 'tskin', 'Ps': 'ps_forc'}[surf_var], self.grid)
                    state[surf_var] = interp([time])[0]
            if self.case.surface_forcing_temp == 'surface_flux':
                interp = self.case.get_interpolator('hfss', self.grid)
                state['sshf'] = interp([time])[0]
            if self.case.surface_forcing_moisture == 'surface_flux':
                interp = self.case.get_interpolator('hfls', self.grid)
                state['slhf'] = interp([time])[0]
            for var in ('z0', 'z0h', 'z0q'):
                if var in state:
                    interp = self.case.get_interpolator(var, self.grid)
                    state[var] = interp([time])[0]

            # Other fields
            interp = self.case.get_interpolator('lat', self.grid)
            state['lat'] = interp([time])[0]
            interp = self.case.get_interpolator('lon', self.grid)
            state['lon'] = interp([time])[0]

        return state
