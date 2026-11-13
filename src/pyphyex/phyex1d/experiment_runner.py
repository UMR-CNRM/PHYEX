"""
ExperimentRunner is the main object of the phyex1d package; it deals with the execution
"""

from pppy import PPPYComp

from .phyex import PhysicsAromeTQ, PhysicsAromeThetaR
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
                              - grid: Grid name or file name containing a grid description
                                      (with '.grid' extension).
                              - dt: Timestep (s)
                              - name: experiment name
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
            cls = {'PhysicsAromeTQ': PhysicsAromeTQ,
                   'PhysicsAromeThetaR': PhysicsAromeThetaR,
                  }[exp['class'] if 'class' in exp else 'PhysicsAromeTQ']
            schemes.append(cls(dt=float(exp.get('dt', 1.)),
                               method='step-by-step',
                               name=exp.get('name', 'Experiment'),
                               tag=string2filename(exp.get('name', 'Experiment')),
                               case=case,
                               grid=Grid(exp.get('grid', 'L90arome')),
                               pyphyex=exp.get('pyphyex', None),
                               pyecrad=exp.get('pyecrad', None),
                               namel=exp.get('namel', 'default'),
                               dx=exp.get('dx', 0),
                               dy=exp.get('dy', 0),
                               attrs=exp.get('attrs', {})))
        super().__init__(schemes=schemes,
                         output_dir=output_dir,
                         duration=self.case.duration,
                         init_state={},
                         name=comp_name,
                         tag=string2filename(comp_name))
