"""
ExperimentRunner is the main object of the phyex1d package; it deals with the execution
"""

import glob
import inspect
import os

import h5py
from pppy import PPPY, PPPYComp

from .phyex import PhysicsAromeTQ, PhysicsAromeThetaR
from .simple_physics import SimplePhysicsTQ, SimplePhysicsThetaR
from .case import Case, CaseCommonFormat
from .grid import Grid
from . import Phyex1DError

# Physics classes that can be designated by their name
PHYSICS_CLASSES = {'PhysicsAromeTQ': PhysicsAromeTQ,
                   'PhysicsAromeThetaR': PhysicsAromeThetaR,
                   'SimplePhysicsTQ': SimplePhysicsTQ,
                   'SimplePhysicsThetaR': SimplePhysicsThetaR}

# Default comparison name (the results are stored in a subdirectory of the output
# directory named after the comparison)
DEFAULT_COMP_NAME = 'phyex1dcomp'

# Arguments of the physics classes that are set by ExperimentRunner
RUNNER_ARGUMENTS = ('dt', 'method', 'name', 'tag', 'case', 'grid', 'attrs')


def string2filename(string):
    """
    Tag (usable in file names) corresponding to a comparison or an experiment name:
    quotes are removed, spaces and slashes are replaced by underscores
    """
    return string.replace('"', '').replace("'", "").replace(' ', '_').replace('/', '_')


def physics_options(cls):
    """
    Options of a physics class: the arguments of its constructor having a default
    value, except the ones set by ExperimentRunner
    :param cls: physics class
    :return: dictionary {option name: default value}
    """
    return {name: param.default
            for name, param in inspect.signature(cls.__init__).parameters.items()
            if name not in RUNNER_ARGUMENTS and param.default is not inspect.Parameter.empty}


class ExperimentRunner(PPPYComp):
    """
    Deals with the execution of the 1D model
    """
    def __init__(self, inputfile_or_case, experiments, output_dir, comp_name=DEFAULT_COMP_NAME):
        """
        :param inputfile_or_case: path to a netCDF driver file or a Case instance
        :param experiments: list of dictionaries describing the experiments
                            allowed keys are:
                              - class: name of the class (key of PHYSICS_CLASSES,
                                       PhysicsAromeTQ by default) or class implementing
                                       the physics
                              - grid: Grid name or file name containing a grid description
                                      (with '.grid' extension).
                              - dt: Timestep (s)
                              - name: experiment name
                              - attrs: dictionary of case attributes to override
                            the other keys are options given to the class (see
                            physics_options)
        :param comp_name: comparison name (DEFAULT_COMP_NAME if empty)
        The results are stored in <output_dir>/<comparison tag>/<experiment tag>.hdf5, the
        tags being the comparison and experiment names with quotes removed and spaces and
        slashes replaced by underscores.
        """
        if not comp_name:
            comp_name = DEFAULT_COMP_NAME

        if isinstance(inputfile_or_case, Case):
            case = inputfile_or_case
        else:
            case = CaseCommonFormat(inputfile_or_case)

        self.case = case
        schemes = []
        for exp in experiments:
            cls = exp.get('class', 'PhysicsAromeTQ')
            if isinstance(cls, str):
                cls = PHYSICS_CLASSES[cls]
            # The experiment keys are given to the class
            kw = dict(exp)
            kw.pop('class', None)
            kw['dt'] = float(kw.get('dt', 1.))  # pppy needs a float
            kw['method'] = 'step-by-step'
            kw['name'] = kw.get('name', 'Experiment')
            kw['tag'] = string2filename(kw['name'])
            kw['case'] = case
            kw['grid'] = Grid(kw.get('grid', 'L90arome'))
            kw['attrs'] = kw.get('attrs', {})
            schemes.append(cls(**kw))
        super().__init__(schemes=schemes,
                         output_dir=output_dir,
                         duration=self.case.duration,
                         init_state={},
                         name=comp_name,
                         tag=string2filename(comp_name))


def list_experiments(output_dir='.', comp_name=DEFAULT_COMP_NAME):
    """
    Existing results of a comparison
    :param output_dir: output directory
    :param comp_name: comparison name
    :return: list of tuples (experiment name, physics class, result file); the name and
             class are '?' for results produced by older versions
    """
    directory = os.path.join(output_dir, string2filename(comp_name or DEFAULT_COMP_NAME))
    if not os.path.isdir(directory):
        raise Phyex1DError(f"No results for the comparison {comp_name} ({directory} does "
                           "not exist)")
    result = []
    for filename in sorted(glob.glob(os.path.join(directory, '*.hdf5'))):
        with h5py.File(filename, 'r') as hdf5:
            result.append((str(hdf5.attrs.get('experiment_name', '?')),
                           str(hdf5.attrs.get('physics_class', '?')), filename))
    return result


class _ExistingResult(PPPY):
    """Experiment known only by its result file (it cannot be run)"""

    def __init__(self, name, tag, dt):
        super().__init__(dt, 'step-by-step', name, tag)

    def step(self, state, timestep, timestep_number):
        raise Phyex1DError(f"Experiment {self.name} cannot be run, only its results can be used")


class ExistingResults(PPPYComp):
    """
    Comparison built from the existing result files (for plotting): neither the case nor
    the description of the experiments are needed
    """
    def __init__(self, output_dir='.', comp_name=DEFAULT_COMP_NAME, experiments=None):
        """
        :param output_dir: output directory
        :param comp_name: comparison name
        :param experiments: list of the names of the experiments to use (all the experiments
                            having results by default); for results produced by older
                            versions (without stored name), the name is the file name
                            (without .hdf5)
        """
        schemes = {}
        for name, _, filename in list_experiments(output_dir, comp_name):
            tag = os.path.splitext(os.path.basename(filename))[0]
            if name == '?':
                name = tag
            with h5py.File(filename, 'r') as hdf5:
                dt = float(hdf5.attrs['dt'])
            schemes[name] = _ExistingResult(name, tag, dt)
        if experiments:
            unknown = [name for name in experiments if name not in schemes]
            if unknown:
                raise Phyex1DError(f"No results for the experiment(s) {', '.join(unknown)}, "
                                   f"available experiments are: {', '.join(schemes)}")
            schemes = {name: schemes[name] for name in experiments}
        if not schemes:
            raise Phyex1DError(f"No results for the comparison {comp_name}")
        comp_name = comp_name or DEFAULT_COMP_NAME
        # The duration is only used to run the experiments
        super().__init__(schemes=list(schemes.values()), output_dir=output_dir,
                         duration=1., init_state={}, name=comp_name,
                         tag=string2filename(comp_name))
