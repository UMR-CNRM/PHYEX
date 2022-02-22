!MNH_LIC Copyright 1995-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ##################
      MODULE MODD_TURB_n
!     ##################
!
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
  REAL               :: XIMPL     ! implicitness degree for the vertical terms of 
                                     ! the turbulence scheme
  REAL               :: XKEMIN      ! mimimum value for the TKE                                  
  REAL               :: XCEDIS      ! Constant for dissipation of Tke                            
  REAL               :: XCADAP      ! Coefficient for ADAPtative mixing length
  CHARACTER (LEN=4)  :: CTURBLEN  ! type of length used for the closure
                                     ! 'BL89' Bougeault and Lacarrere scheme
                                     ! 'DELT' length = ( volum) ** 1/3
  CHARACTER (LEN=4)  :: CTURBDIM  ! dimensionality of the turbulence scheme
                                     ! '1DIM' for purely vertical computations
                                     ! '3DIM' for computations in the 3 
                                     ! directions
  LOGICAL            :: LTURB_FLX ! logical switch for the storage of all  
                                     ! the turbulent fluxes
  LOGICAL            :: LTURB_DIAG! logical switch for the storage of some 
                                     ! turbulence related diagnostics
  LOGICAL            :: LSUBG_COND! Switch for subgrid condensation 
  LOGICAL            :: LSIGMAS   ! Switch for using Sigma_s from turbulence scheme
  LOGICAL            :: LSIG_CONV ! Switch for computing Sigma_s due to convection
!
  LOGICAL            :: LRMC01    ! Switch for computing separate mixing
!                                    ! and dissipative length in the SBL
!                                    ! according to Redelsperger, Mahe &
!                                    ! Carlotti 2001
  CHARACTER(LEN=4)   :: CTOM      ! type of Third Order Moments
                                  ! 'NONE' none
                                  ! 'TM06' Tomas Masson 2006
  CHARACTER(LEN=4)   :: CSUBG_AUCV ! type of subgrid rc->rr autoconv. method
  CHARACTER(LEN=80)  :: CSUBG_AUCV_RI ! type of subgrid ri->rs autoconv. method
  CHARACTER(LEN=80)  :: CCONDENS ! subrgrid condensation PDF
  CHARACTER(LEN=4)   :: CLAMBDA3 ! lambda3 choice for subgrid cloud scheme
  CHARACTER(LEN=80)  :: CSUBG_MF_PDF ! PDF to use for MF cloud autoconversions

!  REAL, DIMENSION(:,:), POINTER :: XBL_DEPTH=>NULL() ! BL depth for TOMS computations
!  REAL, DIMENSION(:,:), POINTER :: XSBL_DEPTH=>NULL()! SurfaceBL depth for RMC01 computations
!  REAL, DIMENSION(:,:,:), POINTER :: XWTHVMF=>NULL()! Mass Flux vert. transport of buoyancy
  REAL               :: VSIGQSAT  ! coeff applied to qsat variance contribution
  REAL, DIMENSION(:,:,:), POINTER :: XDYP=>NULL()    ! Dynamical production of Kinetic energy
  REAL, DIMENSION(:,:,:), POINTER :: XTHP=>NULL()    ! Thermal production of Kinetic energy
  REAL, DIMENSION(:,:,:), POINTER :: XTR=>NULL()    ! Transport production of Kinetic energy
  REAL, DIMENSION(:,:,:), POINTER :: XDISS=>NULL()    ! Dissipation of Kinetic energy
  REAL, DIMENSION(:,:,:), POINTER :: XLEM=>NULL()    ! Mixing length
  REAL, DIMENSION(:,:,:), POINTER :: XSSUFL_C=>NULL() ! O-A interface flux for u
  REAL, DIMENSION(:,:,:), POINTER :: XSSVFL_C=>NULL() ! O-A interface flux for v
  REAL, DIMENSION(:,:,:), POINTER :: XSSTFL_C=>NULL() ! O-A interface flux for theta
  REAL, DIMENSION(:,:,:), POINTER :: XSSRFL_C=>NULL() ! O-A interface flux for vapor
  LOGICAL            :: LHGRAD ! logical switch for the computation of the Leornard Terms
  REAL               :: XCOEFHGRADTHL  ! coeff applied to thl contribution
  REAL               :: XCOEFHGRADRM  ! coeff applied to mixing ratio contribution
  REAL               :: XALTHGRAD  ! altitude from which to apply the Leonard terms
  REAL               :: XCLDTHOLD  ! cloud threshold to apply the Leonard terms
                                   ! negative value : applied everywhere
                                   ! 0.000001 applied only inside the clouds ri+rc > 10**-6 kg/kg
!
END TYPE TURB_t

TYPE(TURB_t), DIMENSION(JPMODELMAX), TARGET, SAVE :: TURB_MODEL

REAL, POINTER :: XIMPL=>NULL()
REAL, POINTER :: XKEMIN=>NULL()
REAL, POINTER :: XCEDIS=>NULL()
REAL, POINTER :: XCADAP=>NULL()
CHARACTER (LEN=4), POINTER :: CTURBLEN=>NULL()
CHARACTER (LEN=4), POINTER :: CTURBDIM=>NULL()
LOGICAL, POINTER :: LTURB_FLX=>NULL()
LOGICAL, POINTER :: LTURB_DIAG=>NULL()
LOGICAL, POINTER :: LSUBG_COND=>NULL()
LOGICAL, POINTER :: LSIGMAS=>NULL()
LOGICAL, POINTER :: LSIG_CONV=>NULL()
LOGICAL, POINTER :: LRMC01=>NULL()
CHARACTER(LEN=4),POINTER :: CTOM=>NULL()
CHARACTER(LEN=4),POINTER :: CSUBG_AUCV=>NULL()
CHARACTER(LEN=80),POINTER :: CSUBG_AUCV_RI=>NULL()
CHARACTER(LEN=80),POINTER :: CCONDENS=>NULL()
CHARACTER(LEN=4),POINTER :: CLAMBDA3=>NULL()
CHARACTER(LEN=80),POINTER :: CSUBG_MF_PDF=>NULL()
REAL, DIMENSION(:,:), POINTER :: XBL_DEPTH=>NULL()
REAL, DIMENSION(:,:), POINTER :: XSBL_DEPTH=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XWTHVMF=>NULL()
REAL, POINTER :: VSIGQSAT=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XDYP=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XTHP=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XTR=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XDISS=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLEM=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XSSUFL_C=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XSSVFL_C=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XSSTFL_C=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XSSRFL_C=>NULL()
LOGICAL, POINTER :: LHGRAD=>NULL()
REAL, POINTER :: XCOEFHGRADTHL=>NULL()
REAL, POINTER :: XCOEFHGRADRM=>NULL()
REAL, POINTER :: XALTHGRAD=>NULL()
REAL, POINTER :: XCLDTHOLD=>NULL()

CONTAINS

SUBROUTINE TURB_GOTO_MODEL(KFROM, KTO)
INTEGER, INTENT(IN) :: KFROM, KTO
!
! Save current state for allocated arrays
!
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
!
! Current model is set to model KTO
XIMPL=>TURB_MODEL(KTO)%XIMPL
XKEMIN=>TURB_MODEL(KTO)%XKEMIN
XCEDIS=>TURB_MODEL(KTO)%XCEDIS
XCADAP=>TURB_MODEL(KTO)%XCADAP
CTURBLEN=>TURB_MODEL(KTO)%CTURBLEN
CTURBDIM=>TURB_MODEL(KTO)%CTURBDIM
LTURB_FLX=>TURB_MODEL(KTO)%LTURB_FLX
LTURB_DIAG=>TURB_MODEL(KTO)%LTURB_DIAG
LSUBG_COND=>TURB_MODEL(KTO)%LSUBG_COND
LSIGMAS=>TURB_MODEL(KTO)%LSIGMAS
LSIG_CONV=>TURB_MODEL(KTO)%LSIG_CONV
LRMC01=>TURB_MODEL(KTO)%LRMC01
CTOM=>TURB_MODEL(KTO)%CTOM
CSUBG_AUCV=>TURB_MODEL(KTO)%CSUBG_AUCV
CSUBG_AUCV_RI=>TURB_MODEL(KTO)%CSUBG_AUCV_RI
CCONDENS=>TURB_MODEL(KTO)%CCONDENS
CLAMBDA3=>TURB_MODEL(KTO)%CLAMBDA3
CSUBG_MF_PDF=>TURB_MODEL(KTO)%CSUBG_MF_PDF
!XBL_DEPTH=>TURB_MODEL(KTO)%XBL_DEPTH !Done in FIELDLIST_GOTO_MODEL
!XSBL_DEPTH=>TURB_MODEL(KTO)%XSBL_DEPTH !Done in FIELDLIST_GOTO_MODEL
!XWTHVMF=>TURB_MODEL(KTO)%XWTHVMF !Done in FIELDLIST_GOTO_MODEL
VSIGQSAT=>TURB_MODEL(KTO)%VSIGQSAT
XDYP=>TURB_MODEL(KTO)%XDYP 
XTHP=>TURB_MODEL(KTO)%XTHP 
XTR=>TURB_MODEL(KTO)%XTR  
XDISS=>TURB_MODEL(KTO)%XDISS
XLEM=>TURB_MODEL(KTO)%XLEM
XSSUFL_C=>TURB_MODEL(KTO)%XSSUFL_C
XSSVFL_C=>TURB_MODEL(KTO)%XSSVFL_C
XSSTFL_C=>TURB_MODEL(KTO)%XSSTFL_C
XSSRFL_C=>TURB_MODEL(KTO)%XSSRFL_C
LHGRAD=>TURB_MODEL(KTO)%LHGRAD
XCOEFHGRADTHL=>TURB_MODEL(KTO)%XCOEFHGRADTHL
XCOEFHGRADRM=>TURB_MODEL(KTO)%XCOEFHGRADRM
XALTHGRAD=>TURB_MODEL(KTO)%XALTHGRAD
XCLDTHOLD=>TURB_MODEL(KTO)%XCLDTHOLD

END SUBROUTINE TURB_GOTO_MODEL

END MODULE MODD_TURB_n
