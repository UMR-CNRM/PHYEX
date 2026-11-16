"""
Case classes representing a simulation case
"""

import argparse
import glob
import os
import re
import subprocess
import sys

import netCDF4
from scipy.interpolate import RegularGridInterpolator

from . import Phyex1DError
from . import user_dir

# DEPHY-SCM repository (cases in the common format)
DEPHY_SCM_URL = 'https://github.com/GdR-DEPHY/DEPHY-SCM.git'
# Driver files of the cases: <CASE>/<SUBCASE>/<CASE>_<SUBCASE>_SCM_driver.nc
_DRIVER_RE = re.compile(r'^(?P<case>[^/]+)/(?P<subcase>[^/]+)/[^/]+_SCM_driver\.nc$')


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


####################################################################################################
# DEPHY-SCM cases
####################################################################################################

def _git(args, directory=None):
    """Executes a git command (in directory), returns its standard output"""
    cmd = ['git'] + ([] if directory is None else ['-C', directory]) + args
    try:
        return subprocess.run(cmd, check=True, stdout=subprocess.PIPE, text=True).stdout
    except (OSError, subprocess.CalledProcessError) as err:
        raise Phyex1DError(f"git command failed: {' '.join(cmd)}") from err


def cases_directory():
    """Directory of the DEPHY-SCM repository: $PHYEXCONF/phyex1d/DEPHY-SCM (PHYEXCONF
    defaults to ~/.phyex)"""
    return user_dir('DEPHY-SCM')


def install_cases(directory=None):
    """
    Installs (clone) or updates (fetch and reset to the remote branch) the DEPHY-SCM
    repository containing the cases in the common format (the repository is big, about
    1 GB).

    :param directory: directory of the repository (cases_directory() by default)
    :return: directory of the repository
    """
    directory = cases_directory() if directory is None else directory
    if os.path.islink(directory):
        raise Phyex1DError(f"{directory} is a link to {os.readlink(directory)}: the "
                           "repository it points to must be updated by its owner (or the "
                           "link must be removed before downloading)")
    if os.path.isdir(os.path.join(directory, '.git')):
        # Update robust to a rewritten history (forced push): the local branch is reset
        # to the remote one (local modifications of the repository are lost)
        _git(['fetch', 'origin'], directory)
        _git(['reset', '--hard', '@{upstream}'], directory)
    else:
        os.makedirs(os.path.dirname(os.path.abspath(directory)), exist_ok=True)
        _git(['clone', DEPHY_SCM_URL, directory])
    return directory


def link_cases(path, directory=None):
    """
    Uses an existing DEPHY-SCM repository (or a copy of it): a symbolic link to it is
    created instead of downloading the repository.

    :param path: path to the existing DEPHY-SCM repository
    :param directory: path of the link (cases_directory() by default); an existing link
                      is replaced, an existing directory (downloaded repository) is kept
                      and an error is raised
    :return: path of the link
    """
    directory = cases_directory() if directory is None else directory
    path = os.path.abspath(path)
    if not os.path.isdir(path):
        raise Phyex1DError(f"{path} is not a directory")
    if not glob.glob(os.path.join(path, '*', '*', '*_SCM_driver.nc')):
        raise Phyex1DError(f"{path} does not contain DEPHY-SCM cases "
                           "(<CASE>/<SUBCASE>/*_SCM_driver.nc files)")
    if os.path.islink(directory):
        os.remove(directory)
    elif os.path.exists(directory):
        raise Phyex1DError(f"{directory} already exists and is not a link (downloaded "
                           "repository?), remove it before creating the link")
    os.makedirs(os.path.dirname(os.path.abspath(directory)), exist_ok=True)
    os.symlink(path, directory)
    return directory


def list_cases(directory=None):
    """
    Lists the cases available in the DEPHY-SCM repository.

    :param directory: directory of the repository (cases_directory() by default)
    :return: list of tuples (case, subcase, path of the driver file)
    """
    directory = cases_directory() if directory is None else directory
    if not os.path.isdir(directory):
        raise Phyex1DError(f"The DEPHY-SCM repository is not installed in {directory}, "
                           "use phyex1d-cases_install")
    result = []
    for path in glob.glob(os.path.join(directory, '*', '*', '*_SCM_driver.nc')):
        match = _DRIVER_RE.match(os.path.relpath(path, directory))
        if match:
            result.append((match.group('case'), match.group('subcase'), path))
    return sorted(result)


def main_cases_install():
    """Entry point of the phyex1d-cases_install command"""
    parser = argparse.ArgumentParser(
        description='Installs the DEPHY-SCM cases (in the common format) in '
                    f'{cases_directory()} (the PHYEXCONF environment variable can be used to '
                    'change the root directory, ~/.phyex by default)')
    subparsers = parser.add_subparsers(dest='command', required=True)
    subparsers.add_parser('download',
                          help='download (clone) or update the DEPHY-SCM repository (big: '
                               'about 1 GB)')
    link = subparsers.add_parser('link', help='use an existing DEPHY-SCM repository (a '
                                              'symbolic link to it is created)')
    link.add_argument('path', help='path to the existing DEPHY-SCM repository')
    args = parser.parse_args()
    try:
        if args.command == 'download':
            directory = install_cases()
            print(f'DEPHY-SCM cases installed in {directory}')
        else:
            directory = link_cases(args.path)
            print(f'{directory} is a link to {os.path.abspath(args.path)}')
    except Phyex1DError as err:
        sys.exit(f'Error: {err}')


def main_cases_list():
    """Entry point of the phyex1d-cases_list command"""
    parser = argparse.ArgumentParser(
        description='Lists the DEPHY-SCM cases installed in '
                    f'{cases_directory()} (see phyex1d-cases_install) with the path of their '
                    'driver file (the case name or the path can be given to phyex1d-exp_run)')
    parser.add_argument('-n', '--names', action='store_true',
                        help='list only the case names (CASE/SUBCASE), without the paths')
    args = parser.parse_args()
    try:
        cases = list_cases()
    except Phyex1DError as err:
        sys.exit(f'Error: {err}')
    names = [case + '/' + subcase for case, subcase, _ in cases]
    width = max((len(name) for name in names), default=0)
    for name, (_, _, path) in zip(names, cases):
        print(name if args.names else f'{name:{width}s}  {path}')
