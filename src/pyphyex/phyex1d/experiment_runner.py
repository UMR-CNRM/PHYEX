"""
ExperimentRunner is the main object of the phyex1d package; it deals with the execution
"""

import inspect

from pppy import PPPYComp

from .phyex import PhysicsAromeTQ, PhysicsAromeThetaR
from .simple_physics import SimplePhysicsTQ, SimplePhysicsThetaR
from .case import Case, CaseCommonFormat
from .grid import Grid

# Physics classes that can be designated by their name
PHYSICS_CLASSES = {'PhysicsAromeTQ': PhysicsAromeTQ,
                   'PhysicsAromeThetaR': PhysicsAromeThetaR,
                   'SimplePhysicsTQ': SimplePhysicsTQ,
                   'SimplePhysicsThetaR': SimplePhysicsThetaR}

# Arguments of the physics classes that are set by ExperimentRunner
RUNNER_ARGUMENTS = ('dt', 'method', 'name', 'tag', 'case', 'grid', 'attrs')


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
    def __init__(self, inputfile_or_case, experiments, output_dir, comp_name=''):
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
        :param comp_name: comparison name
        """
        def string2filename(s):
            return s.replace('"', '').replace("'", "").replace(' ', '_')

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
