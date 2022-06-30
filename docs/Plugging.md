This document must still be written

Some elements to include (or not):

The PHYEX parameterizations can be called from different models (eg. Meso-NH, AROME) and from
a driver (which will be included in this repository).
Moreover, PHYEX parameterizations call externals subroutines, which are dependencies.

This document aims at listing the interfaces to call the PHYEX parameterizations and the
interfaces of the external modules called by PHYEX.

PHYEX interfaces:
- lima_adjust
- ice_adjust
- shallow_mf
- turb
- lima, lima_warm, lima_cold et lima_mixed
- rain_ice, rain_ice_old
- ini_...

Dependencies:
- mode_budget
- mode_msg, modd_io
- modd_precision
- yomhook, parkind1

Specificities:
- in AROME, BUDGET_SORE_INIT does nothing: it is impossible to compute a tendencie
  from the difference of a temporary variable.
  Invalid:
    budget_store_init(tempo_s)
    modification of tempo_s
    budget_store_end(tempo_s)
  Valid:
    budget_store_init(pronostic_s) #useless but valid
    modification of pronostic_s
    budget_store_end(pronostic_s)
  Also valid:
    budget_store_add(delta tempo_s)

