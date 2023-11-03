!MNH_LIC Copyright 2009-2023 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ######################
      MODULE MODI_INI_ELEC_n
!     ######################
!
INTERFACE
      SUBROUTINE INI_ELEC_n (KLUOUT, HELEC, HCLOUD, TPINIFILE, &
                             PTSTEP, PZZ,                      &
                             PDXX, PDYY, PDZZ, PDZX, PDZY      )
!
USE MODD_IO,  ONLY : TFILEDATA
!
INTEGER,           INTENT(IN) :: KLUOUT   ! Logical unit number for prints
CHARACTER (LEN=4), INTENT(IN) :: HELEC    ! atmospheric electricity scheme
CHARACTER (LEN=4), INTENT(IN) :: HCLOUD   ! microphysics scheme
TYPE(TFILEDATA),   INTENT(IN) :: TPINIFILE! Initial file
REAL,              INTENT(IN) :: PTSTEP   ! Time STEP
!
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PZZ     ! height z
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PDXX    ! metric coefficient dxx
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PDYY    ! metric coefficient dyy
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PDZZ    ! metric coefficient dzz
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PDZX    ! metric coefficient dzx
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PDZY    ! metric coefficient dzy
!
END SUBROUTINE INI_ELEC_n
END INTERFACE
END MODULE MODI_INI_ELEC_n
!
!     #########################################################
      SUBROUTINE INI_ELEC_n(KLUOUT, HELEC, HCLOUD, TPINIFILE, &
                            PTSTEP, PZZ,                      &
                            PDXX, PDYY, PDZZ, PDZX, PDZY      )
!     #########################################################
!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to initialize the variables
!     of the atmospheric electricity scheme
!
!!    METHOD
!!    ------
!!      The initialization of the scheme is performed as follows :
!!   
!!    EXTERNAL
!!    --------
!!      CLEANLIST_ll : deaalocate a list
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!
!!    REFERENCE
!!    ---------
!!      
!!    AUTHOR
!!    ------
!!      C. Barthe     * Laboratoire de l'Atmosphère et des Cyclones *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original     09/11/09
!!      M. Chong     13/05/11  Add computation of specific parameters for solving
!!                             the electric field equation (elements of tri-diag
!!                             matrix) 
!!      J.-P. Pinty  13/04/12  Add elec_trid to initialise the tridiagonal syst.
!!      J.-P. Pinty  01/07/12  Add a non-homogeneous Neuman fair-weather 
!!                             boundary condition at the top
!!      J.-P. Pinty  15/11/13  Initialize the flash maps
!!                    10/2016 (C.Lac) Add droplet deposition
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!  P. Wautelet 14/02/2019: remove CLUOUT/CLUOUT0 and associated variables
!  P. Wautelet 10/04/2019: replace ABORT and STOP calls by Print_msg
!!      C. Barthe    04/02/22  Remove call to ini_rain_ice_elec
!!                             Initialization of cloud microphysics and cloud electricity
!!                             is now done separately
!!      C. Barthe    07/07/23  New data structures for some variables
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_ARGSLIST_ll, ONLY: LIST_ll
USE MODD_CLOUDPAR_n, ONLY: NSPLITR
USE MODD_CONF,       ONLY: CEQNSYS, CCONF, CPROGRAM
USE MODD_CONF_n,     ONLY: NRR
USE MODD_CST
USE MODD_DIM_n,      ONLY: NIMAX_ll, NJMAX_ll
USE MODD_DYN
USE MODD_DYN_n,      ONLY: XRHOM, XTRIGSX, XTRIGSY, XAF, XCF, XBFY, XBFB, XDXHATM, &
                           XDYHATM, NIFAXX, NIFAXY, XBF_SXP2_YP1_Z
USE MODD_ELEC_DESCR
USE MODD_ELEC_PARAM
USE MODD_ELEC_FLASH
USE MODD_ELEC_n,     ONLY: XRHOM_E, XAF_E, XCF_E, XBFY_E, XBFB_E, XBF_SXP2_YP1_Z_E
USE MODD_GET_n,      ONLY: CGETSVT, CGETINPRC, CGETINPRR, CGETINPRS, CGETINPRG, CGETINPRH, &
                           CGETCLOUD
USE MODD_GRID_n,     ONLY: XMAP, XDXHAT, XDYHAT
USE MODD_IO,         ONLY: TFILEDATA
USE MODD_LBC_n,      ONLY: CLBCX, CLBCY
USE MODD_LUNIT_n,    ONLY: TLUOUT
USE MODD_PARAMETERS, ONLY: JPVEXT, JPHEXT
USE MODD_PARAM_ICE_n,ONLY : LDEPOSC, LRED
USE MODD_PRECIP_n,   ONLY : XINPRR, XACPRR, XINPRS, XACPRS, XINPRG, XACPRG, &
                            XINPRH, XACPRH, XINPRC, XACPRC, XINPRR3D, XEVAP3D,&
                            XINDEP,XACDEP
USE MODD_REF
USE MODD_REF_n,      ONLY: XRHODJ, XTHVREF
USE MODD_TIME
!
USE MODE_ll
use mode_msg
!
USE MODI_ELEC_TRIDZ
USE MODI_INI_CLOUD
USE MODI_INI_FIELD_ELEC
USE MODI_INI_FLASH_GEOM_ELEC
USE MODI_INI_PARAM_ELEC
USE MODI_INI_RAIN_ICE_ELEC
USE MODI_READ_PRECIP_FIELD
!
!
IMPLICIT NONE
!
!*       0.1   declarations of dummy arguments
!
INTEGER,           INTENT(IN) :: KLUOUT   ! Logical unit number for prints
CHARACTER (LEN=4), INTENT(IN) :: HELEC    ! atmospheric electricity scheme
CHARACTER (LEN=4), INTENT(IN) :: HCLOUD   ! microphysics scheme
TYPE(TFILEDATA),   INTENT(IN) :: TPINIFILE! Initial file
REAL,              INTENT(IN) :: PTSTEP   ! Time STEP
!
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PZZ     ! height z
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PDXX    ! metric coefficient dxx
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PDYY    ! metric coefficient dyy
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PDZZ    ! metric coefficient dzz
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PDZX    ! metric coefficient dzx
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PDZY    ! metric coefficient dzy
!
!*       0.2   declarations of local variables
!
INTEGER :: ILUOUT  ! Logical unit number of output-listing
!
INTEGER :: IIU     ! Upper dimension in x direction (local)
INTEGER :: IJU     ! Upper dimension in y direction (local)
INTEGER :: IKU     ! Upper dimension in z direction
INTEGER :: IKB, IKE
INTEGER :: JK      ! Loop vertical index
INTEGER :: IINFO_ll ! Return code of // routines
INTEGER :: IINTVL   ! Number of intervals to integrate the kernels
REAL    :: ZFDINFTY ! Factor used to define the "infinite" diameter
REAL    :: ZRHO00   ! Surface reference air density
REAL    :: ZDZMIN
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZDZ    ! mesh size
CHARACTER (LEN=3) :: YEQNSYS
!
!
!-------------------------------------------------------------------------------
!
!*       1.    PROLOGUE
!              --------
!
ILUOUT = TLUOUT%NLU
!
CALL GET_DIM_EXT_ll('B',IIU,IJU)
IKU = SIZE(PZZ,3)
IKB = 1 + JPVEXT
IKE = SIZE(PZZ,3) - JPVEXT
!
IF (.NOT.ASSOCIATED(XFCI)) CALL ELEC_PARAM_ASSOCIATE()
IF (.NOT.ASSOCIATED(XFC))  CALL ELEC_DESCR_ASSOCIATE()
!
!-------------------------------------------------------------------------------
!++cb++ 26/04/2023 this part is needed to run the old version of the electrical scheme
! --> use of rain_ice_elec
! --> should be removed when the new scheme is fully validated
!
!*       2.    INITIALIZE THE PARAMETERS FOR THE MICROPHYSICS AND THE ELECTRICITY
!*             IN THE "OLD" ELECTRICAL SCHEME
!              ------------------------------------------------------------------
!
!*       2.1   Allocate module modd_precip_n
!
IF (HELEC == 'ELE3' .AND. (HCLOUD(1:3) == 'ICE' .AND. .NOT.(LRED))) THEN
  ALLOCATE( XINPRR(IIU,IJU) ) ; XINPRR(:,:) = 0.0
  ALLOCATE( XINPRR3D(IIU,IJU,IKU) ) ; XACPRR(:,:) = 0.0
  ALLOCATE( XEVAP3D(IIU,IJU,IKU) ) ; XINPRR3D(:,:,:) = 0.0
  ALLOCATE( XACPRR(IIU,IJU) ) ; XEVAP3D(:,:,:) = 0.0
  ALLOCATE( XINPRC(IIU,IJU) ) ; XINPRC(:,:) = 0.0
  ALLOCATE( XACPRC(IIU,IJU) ) ; XACPRC(:,:) = 0.0
  ALLOCATE( XINPRS(IIU,IJU) ) ; XINPRS(:,:) = 0.0
  ALLOCATE( XACPRS(IIU,IJU) ) ; XACPRS(:,:) = 0.0
  ALLOCATE( XINPRG(IIU,IJU) ) ; XINPRG(:,:) = 0.0
  ALLOCATE( XACPRG(IIU,IJU) ) ; XACPRG(:,:) = 0.0
!
  IF (HCLOUD == 'ICE4') THEN
    ALLOCATE( XINPRH(IIU,IJU) ) ; XINPRH(:,:) = 0.0
    ALLOCATE( XACPRH(IIU,IJU) ) ; XACPRH(:,:) = 0.0
  ELSE
    ALLOCATE( XINPRH(0,0) )
    ALLOCATE( XACPRH(0,0) )
  END IF
!
  IF ( LDEPOSC) THEN
    ALLOCATE(XINDEP(IIU,IJU)) ; XINDEP(:,:) = 0.0
    ALLOCATE(XACDEP(IIU,IJU)) ; XACDEP(:,:) = 0.0
  ELSE
    ALLOCATE(XINDEP(0,0))
    ALLOCATE(XACDEP(0,0))
  END IF
!
  IF(SIZE(XINPRR) == 0) RETURN
!
!*       2.2   Initialize modd_precip_n variables
!
  CALL READ_PRECIP_FIELD (TPINIFILE, CPROGRAM, CCONF,                           &
                          CGETINPRC,CGETINPRR,CGETINPRS,CGETINPRG,CGETINPRH,    &
                          XINPRC,XACPRC,XINDEP,XACDEP,XINPRR,XINPRR3D,XEVAP3D,  &
                          XACPRR, XINPRS, XACPRS, XINPRG, XACPRG, XINPRH, XACPRH)
!
!*       2.3   Initialize the parameters for the microphysics and the electrification
!
! Compute the minimun vertical mesh size
  ALLOCATE( ZDZ(IIU,IJU,IKU) )
  ZDZ(:,:,:) = 0.
!
  DO JK = IKB, IKE
    ZDZ(:,:,JK) = PZZ(:,:,JK+1) - PZZ(:,:,JK)
  END DO
  ZDZMIN = MIN_ll (ZDZ,IINFO_ll,1,1,IKB,NIMAX_ll+2*JPHEXT,NJMAX_ll+2*JPHEXT,IKE )
!
  DEALLOCATE(ZDZ)
!
! initialize the parameters for the mixed-phase microphysics and the electrification
  CALL INI_RAIN_ICE_ELEC (KLUOUT, PTSTEP, ZDZMIN, NSPLITR, HCLOUD, &
                          IINTVL, ZFDINFTY)
END IF
!
!--cb--
!-------------------------------------------------------------------------------
!
!*       2.    INITIALIZE THE PARAMETERS FOR THE ELECTRICITY
!              ---------------------------------------------
!
IF (HELEC(1:3) == 'ELE') THEN
!
!*       2.1   Initialize the electrical parameters for cloud electrification
!
  ZRHO00 = XP00 / (XRD * XTHVREFZ(IKB))
!
  CALL INI_PARAM_ELEC (TPINIFILE, CGETSVT, HCLOUD, HELEC, &
                       ZRHO00, NRR, IIU, IJU, IKU)
!
!
!*       2.2   Initialize the parameters for the electric field
!
  IF (LINDUCTIVE .OR. ((.NOT. LOCG) .AND. LELEC_FIELD)) THEN
    CALL INI_FIELD_ELEC (PDXX, PDYY, PDZZ, PDZX, PDZY, PZZ)
  END IF
!
!
!*       2.3   Initialize the parameters for the lightning flashes
!
  IF (.NOT. LOCG) THEN
    IF (LFLASH_GEOM) THEN
      CALL INI_FLASH_GEOM_ELEC (HCLOUD)
    ELSE
      call Print_msg( NVERB_FATAL, 'GEN', 'INI_ELEC_n', 'INI_LIGHTNING_ELEC not yet developed' )
    END IF
  END IF
!
ELSE IF (HELEC /= 'NONE') THEN
  call Print_msg( NVERB_FATAL, 'GEN', 'INI_ELEC_n', 'not yet developed for CELEC='//trim(HELEC) )
END IF
!
!
!*       2.4   Initialize the parameters for the resolution of the electric field
!
YEQNSYS = CEQNSYS
CEQNSYS = 'LHE'
! Force any CEQNSYS (DUR, MAE, LHE) to LHE to obtain a unique set of coefficients
!    for the flat laplacian operator and Return to the original CEQNSYS
!
ALLOCATE (XRHOM_E(SIZE(XRHOM)))
ALLOCATE (XAF_E(SIZE(XAF)))
ALLOCATE (XCF_E(SIZE(XCF)))
ALLOCATE (XBFY_E(SIZE(XBFY,1),SIZE(XBFY,2),SIZE(XBFY,3)))
ALLOCATE (XBFB_E(SIZE(XBFB,1),SIZE(XBFB,2),SIZE(XBFB,3)))
ALLOCATE (XBF_SXP2_YP1_Z_E(SIZE(XBF_SXP2_YP1_Z,1),SIZE(XBF_SXP2_YP1_Z,2),&
                           SIZE(XBF_SXP2_YP1_Z,3)))
!
CALL ELEC_TRIDZ (CLBCX, CLBCY,                                           &
                 XMAP, XDXHAT, XDYHAT, XDXHATM, XDYHATM, XRHOM_E, XAF_E, & 
                 XCF_E, XTRIGSX, XTRIGSY, NIFAXX, NIFAXY,                &
                 XRHODJ, XTHVREF, PZZ, XBFY_E, XEPOTFW_TOP,              &
                 XBFB_E, XBF_SXP2_YP1_Z_E)
!
CEQNSYS = YEQNSYS
!
!
!*       2.5   Initialize the flash maps
!
ALLOCATE( NMAP_TRIG_IC(IIU,IJU) ); NMAP_TRIG_IC(:,:) = 0
ALLOCATE( NMAP_IMPACT_CG(IIU,IJU) ); NMAP_IMPACT_CG(:,:) = 0
ALLOCATE( NMAP_2DAREA_IC(IIU,IJU) ); NMAP_2DAREA_IC(:,:) = 0
ALLOCATE( NMAP_2DAREA_CG(IIU,IJU) ); NMAP_2DAREA_CG(:,:) = 0
ALLOCATE( NMAP_3DIC(IIU,IJU,IKU) ); NMAP_3DIC(:,:,:) = 0
ALLOCATE( NMAP_3DCG(IIU,IJU,IKU) ); NMAP_3DCG(:,:,:) = 0
!
!-------------------------------------------------------------------------------
! 
END SUBROUTINE INI_ELEC_n
