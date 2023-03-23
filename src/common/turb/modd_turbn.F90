!MNH_LIC Copyright 1995-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ##################
      MODULE MODD_TURB_n
!     ##################
!> @file
!!****  *MODD_TURB$n* - declaration of turbulence scheme free parameters
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to declare the
!     variables that may be set by namelist for the turbulence scheme
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (module MODD_PARAMn)
!!          
!!    AUTHOR
!!    ------
!!	    J. Cuxart and J. Stein       * I.N.M. and Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    January 9, 1995                   
!!      J.Cuxart    February 15, 1995 add the switches for diagnostic storages
!!      J.M. Carriere May  15, 1995 add the subgrid condensation
!!      M. Tomasini Jul  05, 2001 add the subgrid autoconversion
!!      P. Bechtold Feb 11, 2002    add switch for Sigma_s computation
!!      P. Jabouille Apr 4, 2002    add switch for Sigma_s convection
!!      V. Masson    Nov 13 2002    add switch for SBL lengths
!!                   May   2006    Remove KEPS
!!      C.Lac        Nov 2014      add terms of TKE production for LES diag
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!      D. Ricard     May 2021      add the switches for Leonard terms
!!    JL Redelsperger  03/2021   Add O-A flux for auto-coupled LES case
!!
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_PARAMETERS, ONLY: JPMODELMAX
IMPLICIT NONE

TYPE TURB_t
! 
! 
  REAL               :: XIMPL       !< implicitness degree for the vertical terms of the turbulence scheme
  REAL               :: XTKEMIN     !< mimimum value for the TKE                                  
  REAL               :: XCED        !< Constant for dissipation of Tke                            
  REAL               :: XCTP        !< Constant for temperature and vapor pressure-correlations
  REAL               :: XCSHF       !< constant for the sensible heat flux 
  REAL               :: XCHF        !< constant for the humidity flux 
  REAL               :: XCTV        !< constant for the temperature variance
  REAL               :: XCHV        !< constant for the humidity variance
  REAL               :: XCHT1       !< first ct. for the humidity-temperature correlation
  REAL               :: XCHT2       !< second ct. for the humidity-temperature correlation
  REAL               :: XCPR1       !< first ct. for the turbulent Prandtl numbers
  REAL               :: XCADAP      !< Coefficient for ADAPtative mixing length
  CHARACTER (LEN=4)  :: CTURBLEN    !< type of length used for the closure:
                                    !! 'BL89' Bougeault and Lacarrere scheme;
                                    !! 'DELT' length = ( volum) ** 1/3
  CHARACTER (LEN=4)  :: CTURBDIM    !< dimensionality of the turbulence scheme:
                                    !! '1DIM' for purely vertical computations;
                                    !! '3DIM' for computations in the 3 directions
  LOGICAL            :: LTURB_FLX   !< logical switch for the storage of all the turbulent fluxes
  LOGICAL            :: LTURB_DIAG  !< logical switch for the storage of some turbulence related diagnostics
  LOGICAL            :: LSIG_CONV   !< Switch for computing Sigma_s due to convection
!
  LOGICAL            :: LHARAT      !< if true RACMO turbulence is used
  LOGICAL            :: LRMC01      !< Switch for computing separate mixing and dissipative length in the SBL
                                    !! according to Redelsperger, Mahe & Carlotti 2001
  CHARACTER(LEN=4)   :: CTOM        !< type of Third Order Moments:
                                    !! 'NONE' none;
                                    !! 'TM06' Tomas Masson 2006

!  REAL, DIMENSION(:,:), POINTER :: XBL_DEPTH=>NULL() ! BL depth for TOMS computations
!  REAL, DIMENSION(:,:), POINTER :: XSBL_DEPTH=>NULL()! SurfaceBL depth for RMC01 computations
!  REAL, DIMENSION(:,:,:), POINTER :: XWTHVMF=>NULL()! Mass Flux vert. transport of buoyancy
  REAL, DIMENSION(:,:,:), POINTER :: XDYP=>NULL()     !< Dynamical production of Kinetic energy
  REAL, DIMENSION(:,:,:), POINTER :: XTHP=>NULL()     !< Thermal production of Kinetic energy
  REAL, DIMENSION(:,:,:), POINTER :: XTR=>NULL()      !< Transport production of Kinetic energy
  REAL, DIMENSION(:,:,:), POINTER :: XDISS=>NULL()    !< Dissipation of Kinetic energy
  REAL, DIMENSION(:,:,:), POINTER :: XLEM=>NULL()     !< Mixing length
  REAL, DIMENSION(:,:,:), POINTER :: XSSUFL_C=>NULL() !< O-A interface flux for u
  REAL, DIMENSION(:,:,:), POINTER :: XSSVFL_C=>NULL() !< O-A interface flux for v
  REAL, DIMENSION(:,:,:), POINTER :: XSSTFL_C=>NULL() !< O-A interface flux for theta
  REAL, DIMENSION(:,:,:), POINTER :: XSSRFL_C=>NULL() !< O-A interface flux for vapor
  LOGICAL            :: LLEONARD      !< logical switch for the computation of the Leornard Terms
  REAL               :: XCOEFHGRADTHL !< coeff applied to thl contribution
  REAL               :: XCOEFHGRADRM  !< coeff applied to mixing ratio contribution
  REAL               :: XALTHGRAD  !< altitude from which to apply the Leonard terms
  REAL               :: XCLDTHOLD  !< cloud threshold to apply the Leonard terms:
                                   !!  negative value to apply everywhere;
                                   !!  0.000001 applied only inside the clouds ri+rc > 10**-6 kg/kg
  REAL               :: XLINI      !< initial value for BL mixing length
!  
END TYPE TURB_t

TYPE(TURB_t), DIMENSION(JPMODELMAX), TARGET, SAVE :: TURB_MODEL
TYPE(TURB_t), POINTER, SAVE :: TURBN => NULL()
!
REAL, POINTER :: XIMPL=>NULL()
REAL, POINTER :: XTKEMIN=>NULL()
REAL, POINTER :: XCED=>NULL()
REAL, POINTER :: XCTP=>NULL()
REAL, POINTER :: XCSHF=>NULL()
REAL, POINTER :: XCHF=>NULL()
REAL, POINTER :: XCTV=>NULL()
REAL, POINTER :: XCHV=>NULL()
REAL, POINTER :: XCHT1=>NULL()
REAL, POINTER :: XCHT2=>NULL()
REAL, POINTER :: XCPR1=>NULL()
REAL, POINTER :: XCADAP=>NULL()
CHARACTER (LEN=4), POINTER :: CTURBLEN=>NULL()
CHARACTER (LEN=4), POINTER :: CTURBDIM=>NULL()
LOGICAL, POINTER :: LTURB_FLX=>NULL()
LOGICAL, POINTER :: LTURB_DIAG=>NULL()
LOGICAL, POINTER :: LSIG_CONV=>NULL()
LOGICAL, POINTER :: LRMC01=>NULL()
LOGICAL, POINTER :: LHARAT=>NULL()
CHARACTER(LEN=4),POINTER :: CTOM=>NULL()
REAL, DIMENSION(:,:), POINTER :: XBL_DEPTH=>NULL()
REAL, DIMENSION(:,:), POINTER :: XSBL_DEPTH=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XWTHVMF=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XDYP=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XTHP=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XTR=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XDISS=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLEM=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XSSUFL_C=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XSSVFL_C=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XSSTFL_C=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XSSRFL_C=>NULL()
LOGICAL, POINTER :: LLEONARD=>NULL()
REAL, POINTER :: XCOEFHGRADTHL=>NULL()
REAL, POINTER :: XCOEFHGRADRM=>NULL()
REAL, POINTER :: XALTHGRAD=>NULL()
REAL, POINTER :: XCLDTHOLD=>NULL()
REAL, POINTER :: XLINI=>NULL()
!
NAMELIST/NAM_TURBn/XIMPL,CTURBLEN,CTURBDIM,LTURB_FLX,LTURB_DIAG,  &
                   LSIG_CONV,LRMC01,CTOM,&
                   XTKEMIN,XCED,XCTP,XCADAP,&
                   LLEONARD,XCOEFHGRADTHL, XCOEFHGRADRM, &
                   XALTHGRAD, XCLDTHOLD, XLINI,  LHARAT
!
!-------------------------------------------------------------------------------
!
CONTAINS

SUBROUTINE TURB_GOTO_MODEL(KFROM, KTO)
!! This subroutine associate all the pointers to the right component of
!! the right strucuture. A value can be accessed through the structure TURBN
!! or through the strucuture TURB_MODEL(KTO) or directly through these pointers.
IMPLICIT NONE
INTEGER, INTENT(IN) :: KFROM, KTO
!
IF(.NOT. ASSOCIATED(TURBN, TURB_MODEL(KTO))) THEN
!
TURBN => TURB_MODEL(KTO)
!
! Save current state for allocated arrays
!
IF(KFROM>0 .AND. KFROM<=JPMODELMAX) THEN
  !TURB_MODEL(KFROM)%XBL_DEPTH=>XBL_DEPTH !Done in FIELDLIST_GOTO_MODEL
  !TURB_MODEL(KFROM)%XSBL_DEPTH=>XSBL_DEPTH !Done in FIELDLIST_GOTO_MODEL
  !TURB_MODEL(KFROM)%XWTHVMF=>XWTHVMF !Done in FIELDLIST_GOTO_MODEL
  TURB_MODEL(KFROM)%XDYP=>XDYP 
  TURB_MODEL(KFROM)%XTHP=>XTHP 
  TURB_MODEL(KFROM)%XTR=>XTR 
  TURB_MODEL(KFROM)%XDISS=>XDISS
  TURB_MODEL(KFROM)%XLEM=>XLEM
  TURB_MODEL(KFROM)%XSSUFL_C=>XSSUFL_C
  TURB_MODEL(KFROM)%XSSVFL_C=>XSSVFL_C
  TURB_MODEL(KFROM)%XSSTFL_C=>XSSTFL_C
  TURB_MODEL(KFROM)%XSSRFL_C=>XSSRFL_C
ENDIF
!
! Current model is set to model KTO
XIMPL=>TURB_MODEL(KTO)%XIMPL
XTKEMIN=>TURB_MODEL(KTO)%XTKEMIN
XCED=>TURB_MODEL(KTO)%XCED
XCTP=>TURB_MODEL(KTO)%XCTP
XCSHF=>TURB_MODEL(KTO)%XCSHF
XCHF=>TURB_MODEL(KTO)%XCHF
XCTV=>TURB_MODEL(KTO)%XCTV
XCHV=>TURB_MODEL(KTO)%XCHV
XCHT1=>TURB_MODEL(KTO)%XCHT1
XCHT2=>TURB_MODEL(KTO)%XCHT2
XCPR1=>TURB_MODEL(KTO)%XCPR1
XCADAP=>TURB_MODEL(KTO)%XCADAP
CTURBLEN=>TURB_MODEL(KTO)%CTURBLEN
CTURBDIM=>TURB_MODEL(KTO)%CTURBDIM
LTURB_FLX=>TURB_MODEL(KTO)%LTURB_FLX
LHARAT=>TURB_MODEL(KTO)%LHARAT
LTURB_DIAG=>TURB_MODEL(KTO)%LTURB_DIAG
LSIG_CONV=>TURB_MODEL(KTO)%LSIG_CONV
LRMC01=>TURB_MODEL(KTO)%LRMC01
CTOM=>TURB_MODEL(KTO)%CTOM
!XBL_DEPTH=>TURB_MODEL(KTO)%XBL_DEPTH !Done in FIELDLIST_GOTO_MODEL
!XSBL_DEPTH=>TURB_MODEL(KTO)%XSBL_DEPTH !Done in FIELDLIST_GOTO_MODEL
!XWTHVMF=>TURB_MODEL(KTO)%XWTHVMF !Done in FIELDLIST_GOTO_MODEL
XDYP=>TURB_MODEL(KTO)%XDYP 
XTHP=>TURB_MODEL(KTO)%XTHP 
XTR=>TURB_MODEL(KTO)%XTR  
XDISS=>TURB_MODEL(KTO)%XDISS
XLEM=>TURB_MODEL(KTO)%XLEM
XSSUFL_C=>TURB_MODEL(KTO)%XSSUFL_C
XSSVFL_C=>TURB_MODEL(KTO)%XSSVFL_C
XSSTFL_C=>TURB_MODEL(KTO)%XSSTFL_C
XSSRFL_C=>TURB_MODEL(KTO)%XSSRFL_C
LLEONARD=>TURB_MODEL(KTO)%LLEONARD
XCOEFHGRADTHL=>TURB_MODEL(KTO)%XCOEFHGRADTHL
XCOEFHGRADRM=>TURB_MODEL(KTO)%XCOEFHGRADRM
XALTHGRAD=>TURB_MODEL(KTO)%XALTHGRAD
XCLDTHOLD=>TURB_MODEL(KTO)%XCLDTHOLD
XLINI=>TURB_MODEL(KTO)%XLINI
!
ENDIF
!
END SUBROUTINE TURB_GOTO_MODEL

SUBROUTINE TURBN_INIT(HPROGRAM, KUNITNML, LDNEEDNAM, KLUOUT, &
                     &LDDEFAULTVAL, LDREADNAM, LDCHECK, KPRINT)
!!*** *TURBN_INIT* - Code needed to initialize the MODD_TURB_n module
!!
!!*   PURPOSE
!!    -------
!!    Sets the default values, reads the namelist, performs the checks and prints
!!
!!*   METHOD
!!    ------
!!    0. Declarations
!!       1. Declaration of arguments
!!       2. Declaration of local variables
!!    1. Default values
!!    2. Namelist
!!    3. Checks
!!    4. Prints
!!
!!    AUTHOR
!!    ------
!!    S. Riette
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    Feb 2023
!-------------------------------------------------------------------------------
!
!*      0. DECLARATIONS
!       ---------------
!
USE MODE_POSNAM_PHY, ONLY: POSNAM_PHY
USE MODE_MSG, ONLY: PRINT_MSG, NVERB_FATAL
USE MODE_CHECK_NAM_VAL, ONLY: CHECK_NAM_VAL_CHAR, CHECK_NAM_VAL_REAL, CHECK_NAM_VAL_INT
USE MODD_PARAMETERS, ONLY: XUNDEF
!
IMPLICIT NONE
!
!* 0.1. Declaration of arguments
!       ------------------------
!
CHARACTER(LEN=6),  INTENT(IN) :: HPROGRAM     !< Name of the calling program
INTEGER,           INTENT(IN) :: KUNITNML     !< Logical unit to access the namelist
LOGICAL,           INTENT(IN) :: LDNEEDNAM    !< True to abort if namelist is absent
INTEGER,           INTENT(IN) :: KLUOUT       !< Logical unit for outputs
LOGICAL, OPTIONAL, INTENT(IN) :: LDDEFAULTVAL !< Must we initialize variables with default values (defaults to .TRUE.)
LOGICAL, OPTIONAL, INTENT(IN) :: LDREADNAM    !< Must we read the namelist (defaults to .TRUE.)
LOGICAL, OPTIONAL, INTENT(IN) :: LDCHECK      !< Must we perform some checks on values (defaults to .TRUE.)
INTEGER, OPTIONAL, INTENT(IN) :: KPRINT       !< Print level (defaults to 0): 0 for no print, 1 to safely print namelist,
                                              !! 2 to print informative messages
!
!* 0.2 Declaration of local variables
!      ------------------------------
!
LOGICAL :: LLDEFAULTVAL, LLREADNAM, LLCHECK, LLFOUND
INTEGER :: IPRINT

LLDEFAULTVAL=.TRUE.
LLREADNAM=.TRUE.
LLCHECK=.TRUE.
IPRINT=0
IF(PRESENT(LDDEFAULTVAL)) LLDEFAULTVAL=LDDEFAULTVAL
IF(PRESENT(LDREADNAM   )) LLREADNAM   =LDREADNAM
IF(PRESENT(LDCHECK     )) LLCHECK     =LDCHECK
IF(PRESENT(KPRINT      )) IPRINT      =KPRINT
!
!*      1. DEFAULT VALUES
!       -----------------
!
IF(LLDEFAULTVAL) THEN
  XIMPL     = 1.
  XTKEMIN   = 0.01
  XCED      = XUNDEF
  XCTP      = XUNDEF
  XCADAP    = 0.5
  CTURBLEN  = 'BL89'
  CTURBDIM  = '1DIM'
  LTURB_FLX =.FALSE.
  LTURB_DIAG=.FALSE.
  LSIG_CONV =.FALSE.
  LRMC01    =.FALSE.
  CTOM      ='NONE'
  LLEONARD =.FALSE.
  XCOEFHGRADTHL = 1.0
  XCOEFHGRADRM = 1.0
  XALTHGRAD = 2000.0
  XCLDTHOLD = -1.0
  XLINI=0.1 !old value: 10.
  LHARAT=.FALSE.

  IF(HPROGRAM=='AROME') THEN
    XTKEMIN=1.E-6
    XLINI=0.
  ENDIF
ENDIF
!
!*      2. NAMELIST
!       -----------
!
IF(LLREADNAM) THEN
  CALL POSNAM_PHY(KUNITNML, 'NAM_TURBN', LDNEEDNAM, LLFOUND, KLUOUT) 
  IF(LLFOUND) READ(UNIT=KUNITNML, NML=NAM_TURBn) 
ENDIF
!
!*      3. CHECKS
!       ---------
!
IF(LLCHECK) THEN
  CALL CHECK_NAM_VAL_CHAR(KLUOUT, 'CTURBDIM', CTURBDIM, '1DIM', '3DIM')
  CALL CHECK_NAM_VAL_CHAR(KLUOUT, 'CTURBLEN', CTURBLEN, 'DELT', 'BL89', 'RM17', 'DEAR', 'BLKR', 'ADAP')
  CALL CHECK_NAM_VAL_CHAR(KLUOUT, 'CTOM', CTOM, 'NONE', 'TM06')
ENDIF
!
!*      3. PRINTS
!       ---------
!
IF(IPRINT>=1) THEN
  WRITE(UNIT=KLUOUT,NML=NAM_TURBn) 
ENDIF
!
END SUBROUTINE TURBN_INIT
!
END MODULE MODD_TURB_n
