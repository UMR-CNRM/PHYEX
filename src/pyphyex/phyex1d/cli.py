"""
Command line clients: phyex1d-exp_run (execution of experiments), phyex1d-exp_plot (plot of
the results) and phyex1d-exp_list (list of the existing results)
"""

import argparse
import os
import sys

from matplotlib import pyplot as plt
from .experiment_runner import (ExperimentRunner, ExistingResults, PHYSICS_CLASSES,
                                physics_options, list_experiments, DEFAULT_COMP_NAME)
from .case import list_cases
from . import Phyex1DError


def convert_option(cls_name, key, value, default):
    """
    Conversion of a command line value, according to the type of the default value
    of the option in the physics class
    :param cls_name: name of the physics class (for error messages)
    :param key: option name
    :param value: string to convert
    :param default: default value of the option
    :return: converted value
    """
    if isinstance(default, bool):
        if value not in ('True', 'False'):
            raise Phyex1DError(f"Option {key} of {cls_name} must be True or False")
        return value == 'True'
    if isinstance(default, (int, float)):
        try:
            return type(default)(value)
        except ValueError as err:
            raise Phyex1DError(f"Option {key} of {cls_name} must be a "
                               f"{type(default).__name__}") from err
    return value


def physics_options_help():
    """Description of the options of the different physics classes"""
    lines = []
    for cls_name, cls in PHYSICS_CLASSES.items():
        options = ', '.join(f'{key}={default!r}'
                            for key, default in physics_options(cls).items())
        lines.append(f'- {cls_name}: {options}')
    return '\n'.join(lines)


EXPERIMENTS_HELP = """
The case is given by the path to a netCDF file in the common format or by the name of
an installed DEPHY-SCM case (CASE/SUBCASE, e.g. ARMCU/REF, see phyex1d-cases_install and
phyex1d-cases_list).

Each --exp option is followed by a list of key=value strings (without --exp, one
experiment is defined with the default values).
The keys common to all the physics classes are:
- grid: Grid name or file name containing a grid description
        (with '.grid' extension).
- dt: Timestep (s)
- name: experiment name
- class: name of the class that implements a timestep among:
         - PhysicsAromeTQ (default), PhysicsAromeThetaR (PHYEX physics)
         - SimplePhysicsTQ, SimplePhysicsThetaR (simplified physics, without PHYEX)
- attr: override an attribute definig the case, eg.: attr=radiation:on
        (only some of the attributes are implemented)

The other keys are options of the physics class (with their default values):
""" + physics_options_help() + """

with, for the PHYEX physics:
- dx, dy: mesh size (for scale aware parametrizations)
- pyecrad: path to the pyecrad/__init__.py file to use. If not provided, the package must be found
           in the PYTHONPATH in order to use ecrad.
- namel: namelist name or file name (with a .namel extension)

and, for the simplified physics:
- mf_max, mf_dtheta, mf_rt_factor: parameters of the prescribed mass flux
- do_adjustment, do_shallow, do_turb: True or False to activate or not the
  saturation adjustment, the shallow convection and the turbulence

Results are stored in OUTPUT_DIR/<comparison name>/<experiment name>.hdf5 (OUTPUT_DIR is
the current directory by default; quotes are removed, spaces and slashes are replaced by
underscores).
"""

RUN_EXAMPLES = """Examples:
- without comparison name (--name) and without experiment name (name key):
    phyex1d-exp_run ARMCU/REF -o out --exp dt=60
  results in out/""" + DEFAULT_COMP_NAME + """/Experiment.hdf5
- with both names, and a case file:
    phyex1d-exp_run case.nc -o out --name armcu --exp dt=60 name=simple class=SimplePhysicsTQ
  results in out/armcu/simple.hdf5
An experiment is not run again if its result file already exists (use --force).
"""

PLOT_HELP = """
The results are read in OUTPUT_DIR/<comparison name> (see phyex1d-exp_run); all the
experiments having results are plotted, unless some are selected with --exp (their list is
given by phyex1d-exp_list).

Each --plot option is followed by the plot file name, then by options (key=value) and by
the list of variables to plot. Options:
- y_var: variable used for the y-axis (e.g. y_var=P)
- title, xlabel, ylabel, legend, only_param_tag: strings
- mfactor, linewidth: floats
- enable_contourf, switch_cls, colorbar, clabel: True or False
- ymin, ymax: limits of the y-axis
- vmin, vmax: limits of the plotted values
Example:
    phyex1d-exp_plot -o out --exp simple --plot rc.png y_var=P rc
"""


def input_file(case):
    """
    Path of the netCDF file describing the case
    :param case: path to a netCDF file or name of an installed DEPHY-SCM case
                 (CASE/SUBCASE)
    :return: path to the netCDF file
    """
    if os.path.isfile(case):
        return case
    for case_name, subcase_name, path in list_cases():
        if case.strip('/') == case_name + '/' + subcase_name:
            return path
    raise Phyex1DError(f"{case} is neither a file nor an installed DEPHY-SCM case "
                       "(see phyex1d-cases_list)")


def experiments_parser(description, epilog):
    """Parser with the arguments describing the experiments (common to run and plot)"""
    parser = argparse.ArgumentParser(description=description, epilog=epilog,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('CASE',
                        help='Netcdf file describing the case or name of an installed '
                             'DEPHY-SCM case (CASE/SUBCASE)')
    parser.add_argument('--exp', nargs='*', action='append', dest='experiments',
                        help="List of key=value strings defining an experiment. " +
                             "The list of keys is provided at the end of the help message.")
    parser.add_argument('-o', dest='output_dir', type=str, default='.',
                        help="Output directory (default: current directory); it must exist")
    parser.add_argument('--name', default=DEFAULT_COMP_NAME, type=str,
                        help='Name of the comparison (default: %(default)s); results are '
                             'stored in OUTPUT_DIR/NAME/<experiment name>.hdf5')
    return parser


def experiments_from_args(args):
    """List of experiment descriptions (for ExperimentRunner) from the --exp options"""
    experiments = []
    for expe in args.experiments or [[]]:
        values = dict(arg.split('=', 1) for arg in expe if not arg.startswith('attr='))
        cls_name = values.get('class', 'PhysicsAromeTQ')
        if cls_name not in PHYSICS_CLASSES:
            raise Phyex1DError(f"Unknown physics class {cls_name}, allowed classes are: " +
                               ', '.join(PHYSICS_CLASSES))
        options = physics_options(PHYSICS_CLASSES[cls_name])
        expe_descr = {}
        attrs = {}
        for arg in expe:
            key, val = arg.split('=', 1)
            if key == 'attr':
                attr, attrval = val.split(':', 1)
                if attr in ('radiation', ):
                    attrs[attr] = str(attrval)
                else:
                    raise Phyex1DError(f"This experiment attribute is not implemented: {attr}")
            elif key == 'dt':
                expe_descr[key] = float(val)
            elif key in ('name', 'grid', 'class'):
                expe_descr[key] = val
            elif key in options:
                expe_descr[key] = convert_option(cls_name, key, val, options[key])
            else:
                raise Phyex1DError(f"Option {key} does not exist for {cls_name}, allowed "
                                   "options are: dt, name, grid, class, attr, " +
                                   ', '.join(options))
        expe_descr['attrs'] = attrs
        experiments.append(expe_descr)
    return experiments


def runner_from_args(args):
    """ExperimentRunner corresponding to the command line arguments"""
    return ExperimentRunner(input_file(args.CASE), experiments_from_args(args),
                            args.output_dir, args.name)


def plot(ph1d, plot_args):
    """
    Plots the evolution of variables
    :param ph1d: ExperimentRunner instance
    :param plot_args: arguments of a --plot option (file name, options and variables)
    """
    fig, ax = plt.subplots()
    filename = plot_args[0]
    var_names = []
    options = {}
    ylim = [None, None]  # min max for y-axis
    vlim = [None, None]  # min and max value
    for arg in plot_args[1:]:
        if '=' in arg:
            key, value = arg.split('=', 1)
            if key in ('y_var', 'only_param_tag', 'title', 'xlabel', 'ylabel',
                       'legend'):
                # Strings
                options[key] = value
            elif key in ('mfactor', 'linewidth'):
                # Floats
                options[key] = float(value)
            elif key in ('enable_contourf', 'switch_cls', 'colorbar', 'clabel'):
                # Booleans
                options[key] = value == 'True'
            elif key in ('ymin', 'ymax'):
                if key == 'ymin':
                    ylim[0] = float(value)
                else:
                    ylim[1] = float(value)
            elif key in ('vmin', 'vmax'):
                if key == 'vmin':
                    vlim[0] = float(value)
                else:
                    vlim[1] = float(value)
            else:
                raise Phyex1DError(f"This plot option does not exist: {key}")
        else:
            var_names.append(arg)
    kwargs_contourf = {}
    if vlim[0] is not None:
        kwargs_contourf['vmin'] = vlim[0]
    if vlim[1] is not None:
        kwargs_contourf['vmax'] = vlim[1]
    if len(kwargs_contourf) != 0:
        options['kwargs_contourf'] = kwargs_contourf
    ph1d.plot_evol(ax, var_names, **options)
    if options.get('y_var', 'X') == 'P':
        ax.invert_yaxis()
    ax.set_ylim(ylim)
    fig.savefig(filename)
    plt.close(fig)


def main_run():
    """Entry point of the phyex1d-exp_run command: execution of experiments"""
    parser = experiments_parser('Execution of 1D experiments using PHYEX',
                                EXPERIMENTS_HELP + '\n' + RUN_EXAMPLES)
    parser.add_argument('--force', default=False, action='store_true',
                        help="Force run to be executed even if output already exists")
    args = parser.parse_args()
    try:
        ph1d = runner_from_args(args)
    except Phyex1DError as err:
        sys.exit(f'Error: {err}')
    ph1d.run(force=args.force)


def main_plot():
    """Entry point of the phyex1d-exp_plot command: plot of the results of experiments"""
    parser = argparse.ArgumentParser(
        description='Plot of the results of 1D experiments (run with phyex1d-exp_run)',
        epilog=PLOT_HELP, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-o', dest='output_dir', type=str, default='.',
                        help="Output directory (default: current directory)")
    parser.add_argument('--name', default=DEFAULT_COMP_NAME, type=str,
                        help='Name of the comparison (default: %(default)s)')
    parser.add_argument('--exp', action='append', dest='experiments', metavar='EXP_NAME',
                        help='name of an experiment to plot (option can be repeated), all '
                             'the experiments by default')
    parser.add_argument('--plot', nargs='+', action='append', dest='plots', required=True,
                        help="plot file name followed by options and then the list of "
                             "variables to plot (one plot by --plot usage). " +
                             "example: --plot theta.png y_var=P Theta")
    args = parser.parse_args()
    try:
        results = ExistingResults(args.output_dir, args.name, args.experiments)
        for plot_args in args.plots:
            plot(results, plot_args)
    except Phyex1DError as err:
        sys.exit(f'Error: {err}')


def main_list():
    """Entry point of the phyex1d-exp_list command: list of the existing results"""
    parser = argparse.ArgumentParser(
        description='List of the experiments (name and physics class) having results in '
                    'OUTPUT_DIR/<comparison name> (see phyex1d-exp_run); they can be plotted '
                    'with phyex1d-exp_plot')
    parser.add_argument('-o', dest='output_dir', type=str, default='.',
                        help="Output directory (default: current directory)")
    parser.add_argument('--name', default=DEFAULT_COMP_NAME, type=str,
                        help='Name of the comparison (default: %(default)s)')
    args = parser.parse_args()
    try:
        experiments = list_experiments(args.output_dir, args.name)
    except Phyex1DError as err:
        sys.exit(f'Error: {err}')
    if experiments:
        width_name = max(len(name) for name, _, _ in experiments)
        width_cls = max(len(cls) for _, cls, _ in experiments)
        for name, cls, filename in experiments:
            print(f'{name:{width_name}s}  {cls:{width_cls}s}  {filename}')
