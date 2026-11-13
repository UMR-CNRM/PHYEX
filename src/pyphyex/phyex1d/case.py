"""
Case classes representing a simulation case
"""

import netCDF4
from scipy.interpolate import RegularGridInterpolator

from . import Phyex1DError


class Case():
    """
    Base class for case implementations

    All subclasses must provide:
      - description (str property)
      - duration (float property, seconds)
      - contains(varname) -> bool
      - get(varname, slice_idx=None) -> numpy.ndarray
      - get_interpolator(varname, grid) -> RegularGridInterpolator
      - context manager support (__enter__ / __exit__)
    """

    def contains(self, varname):
        """Check if a variable exists in the case data"""
        raise NotImplementedError

    def get(self, varname, slice_idx=None):
        """Read a variable from the case data"""
        raise NotImplementedError

    @property
    def duration(self):
        """Simulation duration in seconds"""
        raise NotImplementedError

    @property
    def description(self):
        """Case description string"""
        raise NotImplementedError

    def get_interpolator(self, varname, grid):
        """Get or build an interpolator for a forcing variable"""
        raise NotImplementedError

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    def close(self):
        """Close any open resources. Override in subclasses."""


class CaseCommonFormat(Case):
    """
    Case backed by a netCDF driver file in the DEPHY common format
    """

    nc_names = {'Theta': 'theta',
                'T': 'ta',
                'Zs': 'orog',
                'Ps': 'ps',
                'qv': 'qv',
                'qc': 'ql',
                'qi': 'qi',
                'qr': 'qr',
                'qs': 'qs',
                'qg': 'qg',
                'qh': 'qh',
                'rv': 'rv',
                'rc': 'rl',
                'ri': 'ri',
                'rr': 'rr',
                'rs': 'rs',
                'rg': 'rg',
                'rh': 'rh',
                'u': 'ua',
                'v': 'va',
                'w': 'wa',
                'sshf': 'hfss',
                'slhf': 'hfls',
                'z0': 'z0',
                'z0h': 'z0h',
                'z0q': 'z0q',
                'tke': 'tke',
                'Tskin': 'tskin',
                'ustar': 'ustar',
                'Ts': 'ts',
                'SolarIrradiance': 'i0',
                'albedo': 'alb',
                'emissivity': 'emis',
                'sza': 'sza',
                'o3': 'o3',
                'lat': 'lat',
                'lon': 'lon',
                'u_geo': 'ug',
                'v_geo': 'vg',
                'ni': 'ni',
                'zh': 'zh',
                'orog': 'orog',
                'pa': 'pa',
                'zh_forc': 'zh_forc',
                'pa_forc': 'pa_forc',
                'time': 'time',
                }

    def __init__(self, inputfile, attrs=None):
        """
        :param inputfile: path to a netCDF driver file
        :param attrs: dictionary of attributes to override
        """
        self._inputfile = inputfile
        self._nc = None
        self._interpolators = {}
        self.adv_ua = 0
        self.adv_va = 0
        with netCDF4.Dataset(inputfile, 'r') as nc:
            for k in nc.ncattrs():
                setattr(self, k, getattr(nc, k))
        if attrs is not None:
            for k, v in attrs.items():
                setattr(self, k, v)

    @property
    def is_open(self):
        """Whether the netCDF file is currently open"""
        return self._nc is not None

    def open(self):
        """
        Open the netCDF file if not already open

        Returns
        -------
        netCDF4.Dataset
            The opened dataset
        """
        if self._nc is None:
            self._nc = netCDF4.Dataset(self._inputfile, 'r')
        return self._nc

    def close(self):
        """Close the netCDF file if open"""
        if self._nc is not None:
            self._nc.close()
            self._nc = None

    def contains(self, varname):
        """
        Check if a variable exists in the netCDF file

        Parameters
        ----------
        varname : str
            Variable name (translated via nc_names if present)

        Returns
        -------
        bool
            True if the variable exists
        """
        self.open()
        ncf_name = self.nc_names.get(varname, varname)
        return ncf_name in self._nc.variables

    def get(self, varname, slice_idx=None):
        """
        Read a variable from the netCDF file

        Parameters
        ----------
        varname : str
            Variable name (translated via nc_names if present)
        slice_idx : int or None, optional
            Slice index to read. None returns the full array.

        Returns
        -------
        numpy.ndarray
            The data array
        """
        self.open()
        ncf_name = self.nc_names.get(varname, varname)
        data = self._nc[ncf_name]
        return data[...] if slice_idx is None else data[slice_idx]

    @property
    def duration(self):
        """
        Simulation duration

        Returns
        -------
        float
            Duration in seconds (time[-1] - time[0])
        """
        times = self.get('time')
        self.close()
        return times[-1] - times[0]

    def get_interpolator(self, varname, grid):
        """
        Get from cache or build an interpolator for a forcing variable

        Parameters
        ----------
        varname : str
            Variable name
        grid : Grid
            Grid instance for vertical coordinate

        Returns
        -------
        RegularGridInterpolator
            The interpolator
        """
        if varname in self._interpolators:
            return self._interpolators[varname]
        input_times = self.get('time')
        if len(self.get(varname).shape) == 2:
            if grid.kind in ('H', 'hybridH'):
                input_coord = self.get('zh_forc', 0)
            elif grid.kind in ('P', 'hybridP'):
                input_coord = self.get('pa_forc', 0)
            else:
                raise Phyex1DError('Wrong grid kind')
            interp = RegularGridInterpolator((input_times, input_coord),
                                             self.get(varname),
                                             bounds_error=False, fill_value=None)
        else:
            interp = RegularGridInterpolator((input_times, ),
                                             self.get(varname),
                                             bounds_error=False, fill_value=None)
        self._interpolators[varname] = interp
        return interp

    @property
    def description(self):
        """
        Case description

        Returns
        -------
        str
            Case description string (file path)
        """
        return self._inputfile


class CaseXarray(Case):
    """
    Case backed by an in-memory xarray Dataset

    Parameters
    ----------
    dataset : xarray.Dataset
        Dataset containing all case variables. Variable names must match the
        phyex1d internal names (e.g. 'T', 'qv', 'u', 'zh', 'Ps', 'lat', ...).
    attrs : dict, optional
        Additional attributes to set on the case (overrides dataset attrs).
    """

    def __init__(self, dataset, attrs=None):
        self._dataset = dataset
        self._interpolators = {}
        self.adv_ua = 0
        self.adv_va = 0
        for k, v in dataset.attrs.items():
            setattr(self, k, v)
        if attrs is not None:
            for k, v in attrs.items():
                setattr(self, k, v)

    def contains(self, varname):
        """
        Check if a variable exists in the dataset

        Parameters
        ----------
        varname : str
            Variable name

        Returns
        -------
        bool
            True if the variable exists
        """
        return varname in self._dataset

    def get(self, varname, slice_idx=None):
        """
        Read a variable from the dataset

        Parameters
        ----------
        varname : str
            Variable name
        slice_idx : int or None, optional
            Slice index to read. None returns the full array.

        Returns
        -------
        numpy.ndarray
            The data array
        """
        data = self._dataset[varname].values
        return data if slice_idx is None else data[slice_idx]

    @property
    def duration(self):
        """
        Simulation duration

        Returns
        -------
        float
            Duration in seconds (time[-1] - time[0])
        """
        times = self._dataset['time'].values
        return float(times[-1] - times[0])

    @property
    def description(self):
        """
        Case description

        Returns
        -------
        str
            Case description string
        """
        return getattr(self, '_description', 'CaseXarray_' + str(id(self._dataset)))

    @description.setter
    def description(self, value):
        """Set a custom description string"""
        self._description = value

    def get_interpolator(self, varname, grid):
        """
        Get from cache or build an interpolator for a forcing variable

        Parameters
        ----------
        varname : str
            Variable name
        grid : Grid
            Grid instance for vertical coordinate

        Returns
        -------
        RegularGridInterpolator
            The interpolator
        """
        if varname in self._interpolators:
            return self._interpolators[varname]
        input_times = self.get('time')
        ndim = len(self.get(varname).shape)
        if ndim == 2:
            if grid.kind in ('H', 'hybridH'):
                input_coord = self.get('zh_forc', 0)
            elif grid.kind in ('P', 'hybridP'):
                input_coord = self.get('pa_forc', 0)
            else:
                raise Phyex1DError('Wrong grid kind')
            interp = RegularGridInterpolator((input_times, input_coord),
                                             self.get(varname),
                                             bounds_error=False, fill_value=None)
        elif ndim == 1:
            interp = RegularGridInterpolator((input_times, ),
                                             self.get(varname),
                                             bounds_error=False, fill_value=None)
        else:
            raise Phyex1DError(f'Unexpected variable dimensions: {ndim}')
        self._interpolators[varname] = interp
        return interp
