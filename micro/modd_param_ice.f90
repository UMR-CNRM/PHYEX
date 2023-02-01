!MNH_LIC Copyright 1995-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ######spl
      MODULE MODD_PARAM_ICE
!     #####################
!
!!****  *MODD_PARAM_ICE* - declaration of the control parameters for the
!!                           mixed phase cloud parameterization
!!
!!    PURPOSE
!!    -------
!!      The purpose of this declarative module is to define the set of space
!!    and time control parameters for the microphysics.
!!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (module MODD_PARAM_ICE)
!!
!!    AUTHOR
!!    ------
!!     J.-P. Pinty   *Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original      14/12/95
!!      Jan 2015 S. Riette: new ICE3/ICE4 parameters
!!      01/10/16 (C.Lac)  Add droplet deposition for fog
!!
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE
!
TYPE PARAM_ICE_t
LOGICAL :: LWARM       ! When .TRUE. activates the formation of rain by
                             ! the warm microphysical processes
LOGICAL :: LSEDIC      ! TRUE to enable the droplet sedimentation
LOGICAL :: LDEPOSC     ! TRUE to enable cloud droplet deposition
REAL    :: XVDEPOSC    ! Droplet deposition velocity
!
CHARACTER(LEN=4) :: CPRISTINE_ICE ! Pristine ice type PLAT, COLU or BURO
CHARACTER(LEN=4) :: CSEDIM        ! Sedimentation calculation mode      
!
LOGICAL :: LRED       ! To use modified ICE3/ICE4 to reduce time step dependency
LOGICAL :: LFEEDBACKT ! When .TRUE. feed back on temperature is taken into account
LOGICAL :: LEVLIMIT   ! When .TRUE. water vapour pressure is limited by saturation
LOGICAL :: LNULLWETG  ! When .TRUE. graupel wet growth is activated with null rate (to allow water shedding)
LOGICAL :: LWETGPOST  ! When .TRUE. graupel wet growth is activated with positive temperature (to allow water shedding)
LOGICAL :: LNULLWETH  ! Same as LNULLWETG but for hail
LOGICAL :: LWETHPOST  ! Same as LWETGPOST but for hail
CHARACTER(LEN=4) :: CSNOWRIMING ! OLD or M90 for Murakami 1990 formulation
REAL :: XFRACM90 ! Fraction used for the Murakami 1990 formulation
INTEGER :: NMAXITER ! Maximum number of iterations for mixing ratio or time splitting
REAL :: XMRSTEP ! maximum mixing ratio step for mixing ratio splitting
LOGICAL :: LCONVHG ! TRUE to allow the conversion from hail to graupel
LOGICAL :: LCRFLIMIT !True to limit rain contact freezing to possible heat exchange
!
REAL :: XTSTEP_TS ! Approximative time step for time-splitting (0 for no time-splitting)
!
CHARACTER(LEN=80) :: CSUBG_RC_RR_ACCR ! subgrid rc-rr accretion
CHARACTER(LEN=80) :: CSUBG_RR_EVAP ! subgrid rr evaporation
CHARACTER(LEN=80) :: CSUBG_PR_PDF ! pdf for subgrid precipitation
!
LOGICAL :: LADJ_BEFORE ! must we perform an adjustment before rain_ice call
LOGICAL :: LADJ_AFTER ! must we perform an adjustment after rain_ice call
CHARACTER(LEN=1) :: CFRAC_ICE_ADJUST ! ice fraction for adjustments
CHARACTER(LEN=1) :: CFRAC_ICE_SHALLOW_MF ! ice fraction for shallow_mf
LOGICAL :: LSEDIM_AFTER ! sedimentation done before (.FALSE.) or after (.TRUE.) microphysics
!
REAL :: XSPLIT_MAXCFL ! Maximum CFL number allowed for SPLIT scheme
LOGICAL :: LSNOW_T         ! Snow parameterization from Wurtz (2021)
!
LOGICAL :: LPACK_INTERP !To pack arrays before computing the different interpolations (kernels and other)
LOGICAL :: LPACK_MICRO !To pack arrays before computing the process tendencies
END TYPE PARAM_ICE_t
!
TYPE(PARAM_ICE_t), SAVE, TARGET :: PARAM_ICE
!
LOGICAL, POINTER :: LWARM => NULL(), &
                    LSEDIC => NULL(), &
                    LDEPOSC => NULL(), &
                    LRED => NULL(), &
                    LFEEDBACKT => NULL(), &
                    LEVLIMIT => NULL(), &
                    LNULLWETG => NULL(), &
                    LWETGPOST => NULL(), &
                    LNULLWETH => NULL(), &
                    LWETHPOST => NULL(), &
                    LCONVHG => NULL(), &
                    LCRFLIMIT => NULL(), &
                    LADJ_BEFORE => NULL(), &
                    LADJ_AFTER => NULL(), &
                    LSEDIM_AFTER => NULL(),&
                    LSNOW_T => NULL(),&
                    LPACK_INTERP => NULL(),&
                    LPACK_MICRO => NULL()

REAL, POINTER :: XVDEPOSC => NULL(), &
                 XFRACM90 => NULL(), &
                 XMRSTEP => NULL(), &
                 XTSTEP_TS => NULL(), &
                 XSPLIT_MAXCFL => NULL()

INTEGER, POINTER :: NMAXITER => NULL()

CHARACTER(LEN=1), POINTER :: CFRAC_ICE_ADJUST => NULL()
CHARACTER(LEN=1), POINTER :: CFRAC_ICE_SHALLOW_MF => NULL()
CHARACTER(LEN=4), POINTER :: CPRISTINE_ICE => NULL()
CHARACTER(LEN=4), POINTER :: CSEDIM => NULL()
CHARACTER(LEN=4), POINTER :: CSNOWRIMING => NULL()
CHARACTER(LEN=80),POINTER :: CSUBG_RC_RR_ACCR => NULL()
CHARACTER(LEN=80),POINTER :: CSUBG_RR_EVAP => NULL()
CHARACTER(LEN=80),POINTER :: CSUBG_PR_PDF => NULL()
!
!-------------------------------------------------------------------------------
!
CONTAINS
SUBROUTINE PARAM_ICE_ASSOCIATE()
  IMPLICIT NONE
  LWARM => PARAM_ICE%LWARM
  LSEDIC => PARAM_ICE%LSEDIC
  LDEPOSC => PARAM_ICE%LDEPOSC
  LRED => PARAM_ICE%LRED
  LFEEDBACKT => PARAM_ICE%LFEEDBACKT
  LEVLIMIT => PARAM_ICE%LEVLIMIT
  LNULLWETG => PARAM_ICE%LNULLWETG
  LWETGPOST => PARAM_ICE%LWETGPOST
  LNULLWETH => PARAM_ICE%LNULLWETH
  LWETHPOST => PARAM_ICE%LWETHPOST
  LCONVHG => PARAM_ICE%LCONVHG
  LCRFLIMIT => PARAM_ICE%LCRFLIMIT
  LADJ_BEFORE => PARAM_ICE%LADJ_BEFORE
  LADJ_AFTER => PARAM_ICE%LADJ_AFTER
  LSEDIM_AFTER => PARAM_ICE%LSEDIM_AFTER
  LSNOW_T => PARAM_ICE%LSNOW_T
  LPACK_INTERP => PARAM_ICE%LPACK_INTERP
  LPACK_MICRO => PARAM_ICE%LPACK_MICRO
  !
  XVDEPOSC => PARAM_ICE%XVDEPOSC
  XFRACM90 => PARAM_ICE%XFRACM90
  XMRSTEP => PARAM_ICE%XMRSTEP
  XTSTEP_TS => PARAM_ICE%XTSTEP_TS
  XSPLIT_MAXCFL => PARAM_ICE%XSPLIT_MAXCFL
  !
  NMAXITER => PARAM_ICE%NMAXITER
  !
  CFRAC_ICE_ADJUST => PARAM_ICE%CFRAC_ICE_ADJUST
  CFRAC_ICE_SHALLOW_MF => PARAM_ICE%CFRAC_ICE_SHALLOW_MF
  CPRISTINE_ICE => PARAM_ICE%CPRISTINE_ICE
  CSEDIM => PARAM_ICE%CSEDIM
  CSNOWRIMING => PARAM_ICE%CSNOWRIMING
  CSUBG_RC_RR_ACCR => PARAM_ICE%CSUBG_RC_RR_ACCR
  CSUBG_RR_EVAP => PARAM_ICE%CSUBG_RR_EVAP
  CSUBG_PR_PDF => PARAM_ICE%CSUBG_PR_PDF
END SUBROUTINE PARAM_ICE_ASSOCIATE
END MODULE MODD_PARAM_ICE
