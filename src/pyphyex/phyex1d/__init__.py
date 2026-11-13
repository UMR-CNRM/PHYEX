"""
1D model using PHYEX
"""

import numpy

__version__ = "0.1.0"


class Phyex1DError(Exception):
    """
    phyex1d error
    """

class Cst():  # pylint: disable=too-few-public-methods,too-many-instance-attributes
    """Empty class to store constants"""
    def __init__(self):
        # pylint: disable=invalid-name
        self.pi = 2. * numpy.arcsin(1.)
        self.Karman = 0.4
        self.Planck = 6.6260755E-34
        self.boltz = 1.380658E-23
        self.avogadro= 6.0221367e+23

        self.P0 = 100000.
        self.g = 9.80665
        self.day = 86400.
        self.siyear = 365.25 * self.day * 2. * self.pi / 6.283076
        self.siday = self.day / (1. + self.day / self.siyear)
        self.omega = 2. * self.pi / self.siday

        self.md= 28.9644e-3
        self.mv= 18.0153E-3
        self.Rd = self.avogadro * self.boltz / self.md
        self.Rv  = self.avogadro * self.boltz / self.mv
        self.Cpd = 7.* self.Rd / 2.
        self.Cpv = 4.* self.Rv
        self.Cl = 4.218E+3
        self.Ci = 2.106E+3
        self.Tt = 273.16
        self.LvTt = 2.5008E+6
        self.LsTt = 2.8345E+6
        self.Lv_0 = self.LvTt - (self.Cpv - self.Cl) * self.Tt
        self.Ls_0 = self.LsTt - (self.Cpv - self.Ci) * self.Tt
        self.EsTt = 611.14
        self.gamw = (self.Cl - self.Cpv) / self.Rv
        self.betaw = (self.LvTt / self.Rv) + (self.gamw * self.Tt)
        self.alpw = numpy.log(self.EsTt) + (self.betaw / self.Tt) + (self.gamw * numpy.log(self.Tt))
        self.gami = (self.Ci - self.Cpv) / self.Rv
        self.betai = (self.LsTt / self.Rv) + (self.gami * self.Tt)
        self.alpi = numpy.log(self.EsTt) + (self.betai / self.Tt) + (self.gami * numpy.log(self.Tt))


from .case import Case, CaseCommonFormat, CaseXarray
from .grid import Grid
from .experiment_runner import ExperimentRunner
from .physics_base import PhysicsBase
from .phyex import (PHYEX, PhysicsArome, PhysicsAromeTQ, PhysicsAromeThetaR,
                    PhysicsForcingTQ, PhysicsForcingThetaR)
