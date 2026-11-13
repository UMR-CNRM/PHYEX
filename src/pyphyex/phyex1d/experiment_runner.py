"""
ExperimentRunner is the main object of the phyex1d package; it deals with the execution
"""

from pppy import PPPYComp

from .phyex import PhysicsAromeTQ, PhysicsAromeThetaR
from .simple_physics import SimplePhysicsTQ, SimplePhysicsThetaR
from .case import Case, CaseCommonFormat
from .grid import Grid


class ExperimentRunner(PPPYComp):
    """
    Deals with the execution of the 1D model
    """
    def __init__(self, inputfile_or_case, experiments, output_dir, comp_name=''):
        """
        :param inputfile_or_case: path to a netCDF driver file or a Case instance
        :param experiments: list of dictionaries describing the experiments
                            allowed keys are:
                              - class: name of the class (PhysicsAromeTQ by default) or
                                       class implementing the physics
                              - grid: Grid name or file name containing a grid description
                                      (with '.grid' extension).
                              - dt: Timestep (s)
                              - name: experiment name
                              - attrs: dictionary of case attributes to override
                            the other keys are given to the class
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
            class_spec = exp.get('class', 'PhysicsAromeTQ')
            if isinstance(class_spec, str):
                cls = {
                    'PhysicsAromeTQ': PhysicsAromeTQ,
                    'PhysicsAromeThetaR': PhysicsAromeThetaR,
                    'SimplePhysicsTQ': SimplePhysicsTQ,
                    'SimplePhysicsThetaR': SimplePhysicsThetaR,
                }[class_spec]
            else:
                cls = class_spec
            # The experiment keys are given to the class (values can be strings when
            # given on the command line)
            kw = dict(exp)
            kw.pop('class', None)
            kw['dt'] = float(kw.get('dt', 1.))
            kw['method'] = 'step-by-step'
            kw['name'] = kw.get('name', 'Experiment')
            kw['tag'] = string2filename(kw['name'])
            kw['case'] = case
            kw['grid'] = Grid(kw.get('grid', 'L90arome'))
            kw['attrs'] = kw.get('attrs', {})
            for k in ('dx', 'dy', 'mf_max', 'mf_dtheta', 'mf_rt_factor'):
                if k in kw:
                    kw[k] = float(kw[k])
            for k in kw:
                if k.startswith('do_') and not isinstance(kw[k], bool):
                    if kw[k] == 'True':
                        kw[k] = True
                    elif kw[k] == 'False':
                        kw[k] = False
                    else:
                        raise ValueError(f"Option {k} must be 'True' or 'False'")
            schemes.append(cls(**kw))
        super().__init__(schemes=schemes,
                         output_dir=output_dir,
                         duration=self.case.duration,
                         init_state={},
                         name=comp_name,
                         tag=string2filename(comp_name))
