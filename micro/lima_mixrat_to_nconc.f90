!MNH_LIC Copyright 2016-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ################################
      MODULE MODI_LIMA_MIXRAT_TO_NCONC
!     ################################
INTERFACE
SUBROUTINE LIMA_MIXRAT_TO_NCONC(PPABST, PTHT, PRVT, PSVT)
!
REAL, DIMENSION(:,:,:),   INTENT(IN)      :: PPABST ! Absolute pressure
REAL, DIMENSION(:,:,:),   INTENT(IN)      :: PTHT   ! Potential temperature
REAL, DIMENSION(:,:,:),   INTENT(IN)      :: PRVT   ! Water Vapor mix. ratio
REAL, DIMENSION(:,:,:,:), INTENT(INOUT)   :: PSVT   ! Mixing ratios IN, conc. OUT
!
END SUBROUTINE LIMA_MIXRAT_TO_NCONC
END INTERFACE
END MODULE MODI_LIMA_MIXRAT_TO_NCONC
!
!     ########################################################
      SUBROUTINE LIMA_MIXRAT_TO_NCONC(PPABST, PTHT, PRVT, PSVT)
!     ########################################################
!
!
!!****  *LIMA_MIXRAT_TO_NCONC* - converts CAMS aerosol mixing ratios into
!!                                 number concentrations 
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    23/01/16 (J.-P. Pinty) 
!-------------------------------------------------------------------------------
!
!*      0. DECLARATIONS
!
USE MODD_CST,        ONLY : XP00, XMD, XMV, XRD, XCPD, XTT, XPI, XRHOLW, &
                            XALPW, XBETAW, XGAMW, XALPI, XBETAI, XGAMI
USE MODD_NSV,        ONLY : NSV_LIMA_CCN_FREE, NSV_LIMA_IFN_FREE
USE MODD_PARAM_LIMA, ONLY : NMOD_CCN, NMOD_IFN,                   &
                            XR_MEAN_CCN, XLOGSIG_CCN, XRHO_CCN,   &
                            XMDIAM_IFN, XSIGMA_IFN, XRHO_IFN,     & 
                            NSPECIE, XFRAC,                       &
                            CCCN_MODES, CIFN_SPECIES
!
IMPLICIT NONE
!
!* 0.1. Declaration of arguments
!       ------------------------
!
REAL, DIMENSION(:,:,:),   INTENT(IN)      :: PPABST ! Absolute pressure
REAL, DIMENSION(:,:,:),   INTENT(IN)      :: PTHT   ! Potential temperature
REAL, DIMENSION(:,:,:),   INTENT(IN)      :: PRVT   ! Water Vapor mix. ratio
REAL, DIMENSION(:,:,:,:), INTENT(INOUT)   :: PSVT   ! Mixing ratios IN, conc. OUT
!
!* 0.2 Declaration of local variables
!      ------------------------------
!
REAL,DIMENSION(SIZE(PTHT,1),SIZE(PTHT,2),SIZE(PTHT,3)) :: ZT    ! Temperature
REAL,DIMENSION(SIZE(PTHT,1),SIZE(PTHT,2),SIZE(PTHT,3)) :: ZREHU ! Relat. Humid.
REAL,DIMENSION(SIZE(PTHT,1),SIZE(PTHT,2),SIZE(PTHT,3)) :: ZGROWTH_FACT
REAL,DIMENSION(SIZE(PTHT,1),SIZE(PTHT,2),SIZE(PTHT,3)) :: ZRHO_CCN_WET
REAL,DIMENSION(SIZE(PTHT,1),SIZE(PTHT,2),SIZE(PTHT,3)) :: ZWORK
!
INTEGER       :: JLOC, JCCN, JIFN, JSPECIE
REAL          :: ZFACT_CCN, ZFACT_IFN
!
!----------------------------------------------------------------------
!
! Temperature to compute the relative humidity
!
ZT(:,:,:) = PTHT(:,:,:)*(PPABST(:,:,:)/XP00)**(XRD/XCPD)
ZWORK(:,:,:) = PRVT(:,:,:)*PPABST(:,:,:)/((XMV/XMD)+PRVT(:,:,:))
                   ! water vapor partial pressure
ZREHU(:,:,:) = ZWORK(:,:,:)/EXP( XALPW-XBETAW/ZT(:,:,:)-XGAMW*ALOG(ZT(:,:,:)) )
                   ! saturation over water
WHERE ( ZT(:,:,:)<XTT )
  ZREHU(:,:,:) = ZWORK(:,:,:)/EXP(XALPI-XBETAI/ZT(:,:,:)-XGAMI*ALOG(ZT(:,:,:)))
                   ! saturation over ice
END WHERE
ZREHU(:,:,:) = MIN( 0.99, MAX( 0.01,ZREHU(:,:,:) ) )
!
! All size distribution parameters are XLOGSIG_CCN and XR_MEAN_CCN (radii)
! Treatment of the soluble aerosols (CCN)
!
! All CAMS aerosol mr are given for dry particles, except for sea-salt (given at Hu=80%) 
!
!

!IF( NAERO_TYPE=="CCN" ) THEN
!
! sea-salt, sulfate, hydrophilic (GADS data)
!
!  NMOD_CCN=3
  IF (.NOT.(ALLOCATED(XR_MEAN_CCN))) ALLOCATE(XR_MEAN_CCN(NMOD_CCN))
  IF (.NOT.(ALLOCATED(XLOGSIG_CCN))) ALLOCATE(XLOGSIG_CCN(NMOD_CCN))
  IF (.NOT.(ALLOCATED(XRHO_CCN)))    ALLOCATE(XRHO_CCN(NMOD_CCN))
  IF( CCCN_MODES=='CAMS_ACC') THEN
      XR_MEAN_CCN(:) = (/ 0.2E-6   , 0.5E-6    , 0.4E-6 /)
      XLOGSIG_CCN(:) = (/ 0.693    , 0.476     , 0.788  /)
      XRHO_CCN(:)    = (/ 2200.    , 1700.     , 1800.  /)
  END IF
!
  IF( CCCN_MODES=='CAMS_AIT') THEN
      XR_MEAN_CCN(:) = (/ 0.2E-6   , 0.05E-6   , 0.02E-6 /)
      XLOGSIG_CCN(:) = (/ 0.693    , 0.693     , 0.788   /)
      XRHO_CCN(:)    = (/ 2200.    , 1700.     , 1800.   /)
  END IF
!
DO JCCN = 1,NMOD_CCN
!
  JLOC = NSV_LIMA_CCN_FREE + JCCN-1 ! CCN free then CCN acti
!
  ZFACT_CCN = ( (0.75/XPI)*EXP(-4.5*(XLOGSIG_CCN(JCCN))**2) )/XR_MEAN_CCN(JCCN)**3
!
! JCCN=1 is for Sea Salt
! JCCN=2 is for Sulphate
! JCCN=3 is for Hydrophilic OC and BC (sulphate coating)
!
  IF( JCCN==1 ) THEN ! Sea salt : convert mass at Hu=80% to dry mass
     PSVT(:,:,:,JLOC) = PSVT(:,:,:,JLOC) / 4.302
  END IF
!
! compute the CCN number concentration
!
! Pourquoi 0.5* ?
! PSVT(:,:,:,JLOC) =0.5* ZFACT_CCN*(PSVT(:,:,:,JLOC)/XRHO_CCN(JCCN)) ! Result 
  PSVT(:,:,:,JLOC) = ZFACT_CCN*(PSVT(:,:,:,JLOC)/XRHO_CCN(JCCN)) ! Result 
                                                       ! is in #/Kg of dry air
END DO
!
! All size distribution parameters are XSIGMA_IFN and XMDIAM_IFN (diameters)
! Treatment of the insoluble aerosols (IFN)
!
!ELSE IF( NAERO_TYPE=="IFN" ) THEN
!
! dust, hydrophobic BIO+ORGA (GADS data)
!
!  NMOD_IFN=2
  NSPECIE=4
  IF (.NOT.(ALLOCATED(XMDIAM_IFN))) ALLOCATE(XMDIAM_IFN(NSPECIE))
  IF (.NOT.(ALLOCATED(XSIGMA_IFN))) ALLOCATE(XSIGMA_IFN(NSPECIE))
  IF (.NOT.(ALLOCATED(XRHO_IFN)))   ALLOCATE(XRHO_IFN(NSPECIE))
  IF( CIFN_SPECIES=='CAMS_ACC') THEN
      XMDIAM_IFN = (/0.8E-6, 3.0E-6, 0.04E-6, 0.8E-6 /)
      XSIGMA_IFN = (/2.0,    2.15,   2.0,     2.2    /)
      XRHO_IFN   = (/2600.,  2600.,  1000.,   2000.  /)
  END IF
  IF( CIFN_SPECIES=='CAMS_AIT') THEN
      XMDIAM_IFN = (/0.8E-6, 3.0E-6, 0.04E-6, 0.04E-6/)
      XSIGMA_IFN = (/2.0,    2.15,   2.0,     2.2 /)
      XRHO_IFN   = (/2600.,  2600.,  1000.,   1800./)
  END IF
  IF (.NOT.(ALLOCATED(XFRAC))) ALLOCATE(XFRAC(NSPECIE,NMOD_IFN))
      XFRAC(1,1)=1.0
      XFRAC(2,1)=0.0
      XFRAC(3,1)=0.0
      XFRAC(4,1)=0.0
      XFRAC(1,2)=0.0
      XFRAC(2,2)=0.0
      XFRAC(3,2)=0.0
      XFRAC(4,2)=1.0
!
DO JIFN = 1,NMOD_IFN
!
! compute the number concentration assuming no deposition of water
! IFN are considered as insoluble dry aerosols
!
  ZFACT_IFN = 0.0
  DO JSPECIE = 1,NSPECIE ! Conversion factor is weighted by XFRAC
    ZFACT_IFN = ZFACT_IFN + XFRAC(JSPECIE,JIFN)*                       &
                ( (6/XPI)*EXP(-(9.0/2.0)*LOG(XSIGMA_IFN(JSPECIE))**2) ) / &
                ( XRHO_IFN(JSPECIE)*XMDIAM_IFN(JSPECIE)**3 ) 
  END DO
  JLOC = NSV_LIMA_IFN_FREE + JIFN-1 ! IFN free then IFN nucl
! Pourquoi 0.5* ?
!  PSVT(:,:,:,JLOC) = 0.5* ZFACT_IFN*PSVT(:,:,:,JLOC) ! Result is in #/Kg of dry air
  PSVT(:,:,:,JLOC) = ZFACT_IFN*PSVT(:,:,:,JLOC) ! Result is in #/Kg of dry air
END DO
!
END SUBROUTINE LIMA_MIXRAT_TO_NCONC
