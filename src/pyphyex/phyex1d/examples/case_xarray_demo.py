"""
Demonstration of CaseXarray with idealised profiles

This example shows how to create a CaseXarray from scratch using xarray,
run a simple experiment, and produce plots. It requires the PHYEX package
to be compiled and importable (see the main README for installation).

Usage:
    source PHYEX/tools/env.sh
    python examples/case_xarray_demo.py
"""

import logging
import os
import sys  # noqa: F401

import numpy as np
import xarray as xr
from matplotlib import pyplot as plt

from phyex1d import CaseXarray, ExperimentRunner

logging.basicConfig(level=logging.INFO)


def make_idealised_dataset():
    """
    Create an xarray Dataset with idealised atmospheric profiles
    on a 90-level height grid up to 20 km.
    """
    nz = 90
    z = np.linspace(0, 20000, nz)  # height above ground (m)
    time = np.array([0., 3600.])    # start and end times (s)

    # Lapse-rate temperature profile
    T = np.broadcast_to(290. - z * 0.0065, (2, nz))
    # Exponentially decreasing moisture
    qv = np.broadcast_to(0.01 * np.exp(-z / 3000.), (2, nz))
    # Zero wind
    u = np.zeros((2, nz))
    v = np.zeros((2, nz))
    w = np.zeros((2, nz))
    tke = np.zeros((2, nz))

    ds = xr.Dataset(
        {
            'Zs': ('time', np.array([0., 0.])),
            'zh': (('time', 'z'), np.broadcast_to(z, (2, nz))),
            'orog': ('time', np.array([0., 0.])),
            'Ps': ('time', np.array([101325., 101325.])),
            'Ts': ('time', np.array([290., 290.])),
            'T': (('time', 'z'), T),
            'qv': (('time', 'z'), qv),
            'u': (('time', 'z'), u),
            'v': (('time', 'z'), v),
            'w': (('time', 'z'), w),
            'tke': (('time', 'z'), tke),
            'lat': ('time', np.array([0., 0.])),
            'lon': ('time', np.array([0., 0.])),
            'ts_forc': ('time', np.array([290., 290.])),
            'ps_forc': ('time', np.array([101325., 101325.])),
            'zh_forc': (('time', 'z'), np.broadcast_to(z, (2, nz))),
        },
        coords={'time': time, 'z': z},
        attrs={
            'adv_theta': 0,
            'adv_ta': 0,
            'adv_rv': 0,
            'adv_qv': 0,
            'adv_ua': 0,
            'adv_va': 0,
            'nudging_theta': 0,
            'nudging_ta': 0,
            'nudging_rv': 0,
            'nudging_qv': 0,
            'nudging_ua': 0,
            'nudging_va': 0,
            'forc_wa': 0,
            'forc_geo': 0,
            'radiation': 'none',
            'surface_forcing_temp': 'ts',
            'surface_forcing_moisture': 'none',
            'surface_forcing_wind': 'none',
            'surface_type': 'ocean',
            'start_date': '2000-01-01T00:00:00',
        },
    )
    return ds


def make_plots(runner, output_dir):
    """Create plots similar to the easy_install script."""
    plots = [
        ('T.png', ['T'], {'y_var': 'P'}),
        ('qv.png', ['qv'], {'y_var': 'P'}),
        ('qc.png', ['qc'], {'y_var': 'P', 'vmin': 0., 'vmax': 5.e-5}),
        ('Theta.png', ['Theta'], {'y_var': 'P'}),
    ]
    for filename, var_names, options in plots:
        fig, ax = plt.subplots()
        kwargs_contourf = {}
        vmin = options.pop('vmin', None)
        vmax = options.pop('vmax', None)
        if vmin is not None:
            kwargs_contourf['vmin'] = vmin
        if vmax is not None:
            kwargs_contourf['vmax'] = vmax
        if kwargs_contourf:
            options['kwargs_contourf'] = kwargs_contourf
        runner.plot_evol(ax, var_names, **options)
        if options.get('y_var', 'X') == 'P':
            ax.invert_yaxis()
        filepath = os.path.join(output_dir, filename)
        fig.savefig(filepath)
        plt.close(fig)
        print(f'  Saved {filepath}')


def main():
    ds = make_idealised_dataset()
    case = CaseXarray(ds)

    print(f'Created CaseXarray (duration: {case.duration} s)')
    print(f'Variables: {list(ds.data_vars)}')

    output_dir = './output_casexarray_demo'
    os.makedirs(output_dir, exist_ok=True)

    runner = ExperimentRunner(
        case,
        [
            {
                'dt': 60,
                'grid': 'L90mesonh',
                'name': 'idealised_demo',
                'class': 'PhysicsAromeTQ',
            },
        ],
        output_dir=output_dir,
        comp_name='demo',
    )

    print('\nRunning experiment...')
    runner.run()
    print('Done!')

    print('\nCreating plots...')
    make_plots(runner, output_dir)
    print('\nAll plots saved in', output_dir)


if __name__ == '__main__':
    main()
