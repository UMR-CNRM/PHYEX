"""
The grib module defines the Grib class that handles the vertical grids:
- reading
- conversion between half/full and Height/Pressure
"""

import os
import glob
import json
import numpy
from . import Phyex1DError
from . import Cst

FLUX = ['FLUX', 'HALF', 'INTERFACE']
MASS = ['MASS', 'FULL']


class Grid():
    """
    The Grid class handles vertical grids

    The definition of the vertical grids are stored in files with the '.grid' extension.
    Theses files contains a dictionary as a json string with the folowing content:
    - kind: 'hybridP' for an hybrid-pressure coordinate (pressure on level i is P(i)=a(i)+b(i)*Ps)
            'hybridH' for an hybrid-height coordinate (altitude on level i is Z(i)=a(i)+b(i)*Zs)
            'P' for a pressure coordinate
            'H' for a height (above sea level) coordinate
    - position: 'FLUX' if the grid refers to flux levels
                'MASS' otherwise
    - mean: 'linear' or 'quadratic'
            way to compute flux levels from mass levels if grid is defined on mass levels
            or the reverse if grid is defined on flux levels
    - description: a single list of values if kind is 'P' or 'H'
                   a list containing a list for the a(i) coefficients and a list for the b(i)
                   coefficients if kind is 'hybridP' or 'hybridH'.
    """

    def __init__(self, grid):
        """
        Instanciates a Grid object from a file containing a grid description or from an
        already known grid.

        :param grid: a filename (if ending with the '.grid' extension or a grid name
        """
        if grid.endswith('.grid') and os.path.exists(grid):
            # grid is a filename provided by the user
            grid_descr = self._read_gridfile(grid)
            self.filename = grid
        else:
            # look for a grid with the name provided in the user's directory,
            # then in the package directory
            user_dir = os.path.join(os.environ['HOME'], '.phyex1d', 'grids')
            package_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'grids')
            for grid_dir in (user_dir, package_dir):
                gridfile = os.path.join(grid_dir, grid + '.grid')
                if os.path.exists(gridfile):
                    grid_descr = self._read_gridfile(gridfile)
                    self.filename = gridfile
                    break
        if grid_descr is None:
            raise Phyex1DError('Grid not found')
        self.position = grid_descr['position']
        self.mean = grid_descr['mean']
        self.kind = grid_descr['kind']
        self.description = grid_descr['description']

        assert self.position in FLUX + MASS, 'grid must be defined on flux or mass points'
        assert self.mean in ('linear', 'quadratic'), \
            "grid mean attribute must be 'quadratic' or 'linear'"
        assert self.kind in ('H', 'P', 'hybridH', 'hybridP'), "grid kind unknown"
        if self.kind in ('H', 'P'):
            assert isinstance(self.description, list), 'grid description must be a list'
            self.len = len(self.description)
            self.description = numpy.array(self.description)
        else:
            assert isinstance(self.description, list) and len(self.description) == 2, \
                "grid description must be a two-element list"
            assert isinstance(self.description[0], list) and \
                   isinstance(self.description[1], list), \
                   "Both elements of grid description must be lists"
            self.len = len(self.description[0])
            assert self.len == len(self.description[1]), \
                "Both elements of grid description must have the same length"
        assert self.len > 0, 'grid must contain at least one level'

        self.ascending = None

        self.cst = Cst()

    def _read_gridfile(self, gridfile):
        """
        Reads a grid
        """
        with open(gridfile, 'r', encoding='UTF-8') as f:
            grid = json.load(f)
        return grid

    def _mean(self, values, method):
        """
        Computes the mean
        :param method: must be 'linear' or 'quadartic'
        """
        if method == 'linear':
            values = numpy.array([.5 * (values[i] + values[i + 1])
                                  for i in range(self.len - 1)])
        elif method == 'quadratic':
            if values[0] == 0.:
                update = (0, values[1] / 2.)  # Linear when quadratic is not possible
            elif values[-1] == 0.:
                update = (-1, values[-2] / 2.)  # Linear when quadratic is not possible
            else:
                update = None
            values = numpy.array([numpy.sqrt(values[i] * values[i + 1])
                                  for i in range(self.len - 1)])
            if update is not None:
                values[update[0]] = update[1]
        else:
            raise Phyex1DError("method must be 'linear' or 'quadartic'")
        return values

    def get_pressure(self, position, state, pronostic_variables):
        """
        Returns the level pressures
        :param position: must be among grid.MASS + .grid.FLUX
        :param state: dictionary containing the current state
        :param pronostic_variables: list of pronostic variables
        """
        if self.kind == 'hybridP' and position == self.position:
            if 'Ps' not in state:
                raise Phyex1DError("Surface pressure ('Ps') is needed to compute the pressures")
            pressure = numpy.array([self.description[0][i] + self.description[1][i] * state['Ps']
                                    for i in range(self.len)])
            if self.ascending is None:
                self.ascending = numpy.all(pressure[1:] < pressure[:-1])
        elif self.kind == 'P' and position == self.position:
            pressure = numpy.array(self.description)
            if self.ascending is None:
                self.ascending = numpy.all(pressure[1:] < pressure[:-1])
        elif self.kind in ('P', 'hybridP'):
            pressure = self._mean(self.get_pressure(self.position, state, pronostic_variables),
                                  self.mean)
        elif self.kind in ('H', 'hybridH'):
            # We assume hydrostatic equilibrium
            # leading to dln(P) / dz = -g / (R * T) if T is the prognostic variable
            # or d[(P / P0) ** (Rd / Cpd)] / dz = -g / (R * Theta) * (Rd / Cpd) with Theta
            if 'qv' in pronostic_variables:
                gas_constant = self.cst.Rd + state['qv'] * (self.cst.Rv - self.cst.Rd)
                for var in ('qc', 'qr', 'qi', 'qs', 'qg', 'qh'):
                    if var in pronostic_variables:
                        gas_constant += - state[var] * self.cst.Rd
            elif 'rv' in pronostic_variables:
                div = 1 + state['rv']
                for var in ('rc', 'rr', 'ri', 'rs', 'rg', 'rh'):
                    if var in pronostic_variables:
                        div += state[var]
                gas_constant = (self.cst.Rd + state['rv'] * self.cst.Rv) / div
            z_flux = self.get_altitude('FLUX', state, pronostic_variables)
            dx = - self.cst.g / gas_constant * numpy.diff(z_flux)
            if 'T' in pronostic_variables:
                dx = dx / state['T']
            elif 'Theta' in pronostic_variables:
                dx = dx / state['Theta'] * self.cst.Rd / self.cst.Cpd
            else:
                raise Phyex1DError("T or Theta must be prognostic")
            x_flux = numpy.zeros(z_flux.shape)
            if self.ascending:
                x_flux[1:] = numpy.cumsum(dx)
            else:
                x_flux[:-1] = numpy.cumsum(-dx[::-1])[::-1]
            if 'T' in pronostic_variables:
                x_flux += numpy.log(state['Ps'])
                pressure_flux = numpy.exp(x_flux)
            else:
                x_flux += (state['Ps'] / 1.E5) ** (self.cst.Rd / self.cst.Cpd)
                pressure_flux = x_flux ** (self.cst.Cpd / self.cst.Rd) * 1.E5
            if position == 'FLUX':
                pressure = pressure_flux
            else:
                # We add the contribution between flux and mass levels assuming constant profiles
                z_mass = self.get_altitude('MASS', state, pronostic_variables)
                dx = - self.cst.g / gas_constant
                if self.ascending:
                    dx = dx * (z_mass - z_flux[:-1])
                else:
                    dx = dx * (z_mass - z_flux[1:])
                if self.ascending:
                    temp = pressure_flux[:-1]
                else:
                    temp = pressure_flux[1:]
                if 'T' in pronostic_variables:
                    dx = dx / state['T']
                    pressure = numpy.exp(numpy.log(temp) + dx)
                elif 'Theta' in pronostic_variables:
                    dx = dx / state['Theta'] * self.cst.Rd / self.cst.Cpd
                    pressure = ((temp / 1.E5) ** (self.cst.Rd / self.cst.Cpd) + dx
                               ) ** (self.cst.Cpd / self.cst.Rd) * 1.E5
                else:
                    raise Phyex1DError("T or Theta must be prognostic")
        else:
            raise Phyex1DError("Wrong grid's kind")

        return pressure

    def get_altitude(self, position, state, pronostic_variables):
        """
        Returns the level altitudes
        :param position: must be among grid.MASS + .grid.FLUX
        :param state: dictionary containing the current state
        :param pronostic_variables: list of pronostic variables
        """
        if self.kind == 'hybridH' and position == self.position:
            if 'Zs' not in state:
                raise Phyex1DError("Surface altitude ('Zs') is needed to compute the altitudes")
            altitude = numpy.array([self.description[0][i] + self.description[1][i] * state['Zs']
                                    for i in range(self.len)])
            if self.ascending is None:
                self.ascending = numpy.all(altitude[1:] > altitude[:-1])
        elif self.kind == 'H' and position == self.position:
            altitude = numpy.array(self.description)
            if self.ascending is None:
                self.ascending = numpy.all(altitude[1:] > altitude[:-1])
        elif self.kind in ('H', 'hybridH'):
            altitude = self._mean(self.get_altitude(self.position, state, pronostic_variables),
                                  self.mean)
        elif self.kind in ('P', 'hybridP'):
            # We assume hydrostatic equilibrium
            # leading to dln(P) / dz = -g / (R*T)
            pressure_mass = None
            if 'T' in pronostic_variables:
                temperature = state['T']
            elif 'Theta' in pronostic_variables:
                pressure_mass = self.get_pressure('MASS', state, pronostic_variables)
                exner = (pressure_mass / 1.E5) ** (self.cst.Rd / self.cst.Cpd)
                temperature = exner * state['Theta']
            else:
                raise Phyex1DError("T or Theta must be prognostic")
            if 'qv' in pronostic_variables:
                gas_constant = self.cst.Rd + state['qv'] * (self.cst.Rv - self.cst.Rd)
                for var in ('qc', 'qr', 'qi', 'qs', 'qg', 'qh'):
                    if var in pronostic_variables:
                        gas_constant += - state[var] * self.cst.Rd
            elif 'rv' in pronostic_variables:
                div = 1 + state['rv']
                for var in ('rc', 'rr', 'ri', 'rs', 'rg', 'rh'):
                    if var in pronostic_variables:
                        div += state[var]
                gas_constant = (self.cst.Rd + state['rv'] * self.cst.Rv) / div
            pressure_flux = self.get_pressure('FLUX', state, pronostic_variables)
            dz = -gas_constant * temperature / self.cst.g * numpy.diff(numpy.log(pressure_flux))
            z_flux = numpy.zeros(pressure_flux.shape)
            if self.ascending:
                z_flux[1:] = numpy.cumsum(dz)
            else:
                z_flux[:-1] = numpy.cumsum(-dz[::-1])[::-1]
            z_flux += state['Zs']
            if position == 'FLUX':
                altitude = z_flux
            else:
                # We add the contribution between flux and mass levels assuming constant profiles
                if pressure_mass is None:
                    pressure_mass = self.get_pressure('MASS', state, pronostic_variables)
                if self.ascending:
                    dlnp = numpy.log(pressure_mass) - numpy.log(pressure_flux[:-1])
                else:
                    dlnp = numpy.log(pressure_mass) - numpy.log(pressure_flux[1:])
                dz = -gas_constant * temperature / self.cst.g * dlnp
                if self.ascending:
                    altitude = z_flux[:-1] + dz
                else:
                    altitude = z_flux[1:] + dz
        else:
            raise Phyex1DError("Wrong grid's kind")
        return altitude
