!MNH_LIC Copyright 2002-2023 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!      #######################
       MODULE  MODD_ELEC_PARAM
!      #######################
!
!!****  *MODD_ELEC_PARAM* - declaration of some electrical factors
!!			    extensively used in the electrical scheme.
!!
!!      PURPOSE
!!      -------
!	  The purpose of this declarative module is to declare some precomputed
!	electrical parameters directly used in routines related to cloud electricity
!!
!!**    IMPLICIT ARGUMENTS
!!      ------------------
!!	  None
!!
!!      REFERENCE
!!      ---------
!!
!!      AUTHOR
!!      ------
!!        Gilles Molinie    * Laboratoire d'Aerologie *
!!        
!!
!!      MODIFICATIONS
!!      -------------
!!        Original      14/11/02
!!        C. Barthe     31/01/2022   add XFQUPDNCI 
!!        C. Barthe     07/06/2022   add parameters for charge sedimentation in LIMA
!!        C. Barthe     28/03/2023   add parameters for sedimentation of cloud droplets charge
!!        C. Barthe     05/07/2023   new data structures for PHYEX - for sedimentation in ICE3
!!
!-------------------------------------------------------------------------------
!
!*	0.	DECLARATIONS
!		------------
!
USE MODD_PARAMETERS, ONLY: JPMODELMAX
!
IMPLICIT NONE
!
SAVE
!
REAL :: XCOEF_RQ_V, XCOEF_RQ_C, &   ! Constants for proportionality
        XCOEF_RQ_R, XCOEF_RQ_I, &   ! between mass transfer and
        XCOEF_RQ_S, XCOEF_RQ_G, &   ! charge transfer
        XCOEF_RQ_H
!
REAL :: XQHON    ! Constant for spontaneous freezing of droplets if T<-35°
!
REAL, DIMENSION(:), ALLOCATABLE :: XFQSED  ! Constant for sedimentation of 
                                           ! electric charge in LIMA
REAL, DIMENSION(:), ALLOCATABLE :: XDQ     ! Exponent for sedimentation of
                                           ! electric charge in LIMA
REAL :: XFQUPDNCI ! constant used to update e_i for sedimentation where
                  ! N_i follows McFarquhar and Heysmfield (1997)
!
REAL :: XQSRIMCG, XEXQSRIMCG       ! Constant for riming of cloud droplets
                                   ! on snow 
REAL, DIMENSION(:), ALLOCATABLE :: XGAMINC_RIM3  
!
REAL :: XQRCFRIG, XEXQRCFRIG       ! Constant for contact freezing between
                                   ! raindrops and pristine ice
REAL :: XFQRACCS                   ! Constant in RACCS
!
REAL :: XFQIAGGSBH,                &         ! Constant for IAGGS charging
        XFQIAGGSBG, XEXFQIAGGSBG,  &         ! process for HELFA, GARDI,
        XFQIAGGSBS,                &         ! SAUND and TAKAH
        XFQIAGGSBT1, XFQIAGGSBT2, XFQIAGGSBT3 
!
REAL :: XLBQRACCS1, XLBQRACCS2, XLBQRACCS3    ! Integral of normalization 
REAL :: XLBQSACCRG1, XLBQSACCRG2, XLBQSACCRG3 ! in accretion of raindrops
                                              ! on snow process 
!                         
REAL, DIMENSION(:,:), ALLOCATABLE  &               ! Normalized kernel for
     :: XKER_Q_RACCS, XKER_Q_RACCSS, XKER_Q_SACCRG ! RACCS, RACCSS, SACCRG
!
REAL :: XFQSDRYG, XFQSDRYGB, XFQRDRYG        ! Constant in SDRYG and RDRYG
!
! charge separation
!
REAL, DIMENSION(:,:), ALLOCATABLE :: XKER_Q_LIMSG
REAL, DIMENSION(:,:), ALLOCATABLE :: XKER_Q_SDRYGB, XKER_Q_SDRYGB1, XKER_Q_SDRYGB2 
!
! Helsdon-Farley
!
REAL :: XHIDRYG     ! Constant charge separated 
REAL :: XHSDRYG     ! Constant charge separated 
REAL :: XLBQSDRYGB4H, XLBQSDRYGB5H, XLBQSDRYGB6H    ! Constants in QIDRYGB
REAL :: XFQSDRYGBH                                  !
!
! Gardiner
!
REAL :: XLWCC  ! LWC critic in Gardiner NI charging
REAL :: XFQIDRYGBG, XLBQIDRYGBG                     ! Constants in QIDRYGB
REAL :: XFQSDRYGBG                                  ! Constants in QSDRYGB
REAL :: XLBQSDRYGB4G, XLBQSDRYGB5G, XLBQSDRYGB6G    ! 
!
! Saunders
!
REAL :: XIMP, XINP, XIKP, &   ! Parameters m, n and k
        XIMN, XINN, XIKN, &   ! for the NI processes
        XSMP, XSNP, XSKP, &   ! following
        XSMN, XSNN, XSKN      ! Saunders et al. (1991)
REAL :: XFQIAGGSP, XFQIAGGSN,         & ! Auxiliary parameters
        XFQIDRYGBSP, XFQIDRYGBSN,     & ! containing MOMG function
        XLBQSDRYGB1SP, XLBQSDRYGB1SN, &
        XLBQSDRYGB2SP, XLBQSDRYGB2SN, &
        XLBQSDRYGB3SP, XLBQSDRYGB3SN, &
        XAIGAMMABI
REAL :: XIKP_TAK, XIKN_TAK, XSKP_TAK, XSKN_TAK ! Using Takahashi charge
REAL :: XFQIAGGSP_TAK, XFQIAGGSN_TAK, XFQIDRYGBSP_TAK, XFQIDRYGBSN_TAK
REAL :: XVSCOEF, XVGCOEF
REAL :: XFQIDRYGBS,   XLBQIDRYGBS    ! Constants in QIDRYGB
REAL :: XFQSDRYGBS                   ! Constants in QSDRYGB
REAL :: XLBQSDRYGB1S, XLBQSDRYGB2S   !
!
! Takahashi
!
INTEGER :: NIND_TEMP ! number of indexes for temperature
INTEGER :: NIND_LWC  ! number of indexes for liquid water content
REAL, DIMENSION(:,:), ALLOCATABLE :: XMANSELL ! F(LWC, T) for Takahashi(1978) /Mansell
REAL, DIMENSION(:,:), ALLOCATABLE :: XSAUNDER ! F(LWC, T) for SAUN1/SAUN2, BSMP1/BSMP2
REAL, DIMENSION(:,:), ALLOCATABLE :: XTAKA_TM ! F(LWC, T) for Takahashi/Tsenova and Mitzeva
!
REAL :: XFQIDRYGBT1,  XFQIDRYGBT2,  XFQIDRYGBT3, &  ! IDRYGB
        XFQSDRYGBT1,  XFQSDRYGBT2,  XFQSDRYGBT3, &  ! SDRYGB
        XFQSDRYGBT4,  XFQSDRYGBT5,  XFQSDRYGBT6, &  ! SDRYGB
        XFQSDRYGBT7,  XFQSDRYGBT8,  XFQSDRYGBT9, &  ! SDRYGB
        XFQSDRYGBT10, XFQSDRYGBT11, XFQSDRYGBT12    ! SDRYGB
!
REAL :: XLBQRDRYG1, XLBQRDRYG2, XLBQRDRYG3  ! Integral of normalization in
REAL :: XLBQSDRYG1, XLBQSDRYG2, XLBQSDRYG3  ! the accretion of graupel on
                                                  ! raindrop and snow process
!
REAL, DIMENSION(:,:), ALLOCATABLE &  
           :: XKER_Q_SDRYG, XKER_Q_RDRYG ! Normalized kernel for SDRYG and RDRYG
!
REAL :: XQREVAV1, XQREVAV2                  ! Raindrops evaporation
!
! Add variables to limit the exchanged charge
!
REAL :: XAUX_LIM
REAL :: XAUX_LIM1, XAUX_LIM2, XAUX_LIM3
!
!
! Inductive charging process
!
REAL :: XCOLCG_IND     ! collision effiency
REAL :: XEBOUND        ! rebound efficiency
REAL :: XALPHA_IND     ! fraction of droplets with grazing trajectories
REAL :: XCOS_THETA     ! average cosine of the angle of rebounding collision
REAL :: XIND1, XIND2, XIND3
!
! lightning
!
REAL :: XFQLIGHTC, XFQLIGHTR, XFQLIGHTI, &
        XFQLIGHTS, XFQLIGHTG, XFQLIGHTH     ! Constant for charge redistribution
REAL :: XEXQLIGHTR, XEXQLIGHTI, &
        XEXQLIGHTS, XEXQLIGHTG, XEXQLIGHTH  ! Exponent for charge redistribution
!
! The following variables must be declared with a derived type to match with PHYEX requirements
TYPE ELEC_PARAM_t
  REAL :: XFCI    ! Constant for sedimentation of the mixing ratio of ice
                  ! which the computation is modified in regard of rain_ice.f90
  !
  REAL :: XFQSEDC, XEXQSEDC,    &    ! Constant for sedimentation of cloud droplets
          XFQSEDR, XEXQSEDR,    &    !  rain
          XFQSEDI, XEXQSEDI,    &    !  ice
          XFQSEDS, XEXQSEDS,    &    !  snow
          XFQSEDG, XEXQSEDG,    &    !  graupel
          XFQSEDH, XEXQSEDH          !  hail
  !
  REAL :: XEGMIN, XEGMAX, XESMIN, XESMAX, &  ! Max and min values for
          XEIMIN, XEIMAX, XECMIN, XECMAX, &  ! e_x in q=e_x D^f_x
          XERMIN, XERMAX, XEHMIN, XEHMAX
  !
  REAL :: XFQUPDC, XFQUPDR, XFQUPDI,&         ! Update Q=f(D)
          XEXFQUPDI, XFQUPDS, XFQUPDG, XFQUPDH
END TYPE ELEC_PARAM_t
!
TYPE(ELEC_PARAM_t), SAVE, TARGET :: ELEC_PARAM
!
REAL, POINTER :: XFCI => NULL(),     &
                 XFQSEDC => NULL(),  &
                 XEXQSEDC => NULL(), &
                 XFQSEDR => NULL(),  &
                 XEXQSEDR => NULL(), &
                 XFQSEDI => NULL(),  &
                 XEXQSEDI => NULL(), &
                 XFQSEDS => NULL(),  &
                 XEXQSEDS => NULL(), &
                 XFQSEDG => NULL(),  &
                 XEXQSEDG => NULL(), &
                 XFQSEDH => NULL(),  &
                 XEXQSEDH => NULL(), &
                 XEGMIN => NULL(),   &
                 XEGMAX => NULL(),   &
                 XESMIN => NULL(),   &
                 XESMAX => NULL(),   &
                 XEIMIN => NULL(),   &
                 XEIMAX => NULL(),   &
                 XECMIN => NULL(),   &
                 XECMAX => NULL(),   &
                 XERMIN => NULL(),   &
                 XERMAX => NULL(),   &
                 XEHMIN => NULL(),   &
                 XEHMAX => NULL(),   &
                 XFQUPDC => NULL(),  &
                 XFQUPDR => NULL(),  &
                 XFQUPDI => NULL(),  &
                 XEXFQUPDI => NULL(),&
                 XFQUPDS => NULL(),  &
                 XFQUPDG => NULL(),  &
                 XFQUPDH => NULL()
!
CONTAINS
!
SUBROUTINE ELEC_PARAM_ASSOCIATE()
  IMPLICIT NONE
  !
  XFCI => ELEC_PARAM%XFCI
  XFQSEDC => ELEC_PARAM%XFQSEDC
  XEXQSEDC => ELEC_PARAM%XEXQSEDC
  XFQSEDR => ELEC_PARAM%XFQSEDR 
  XEXQSEDR => ELEC_PARAM%XEXQSEDR
  XFQSEDI => ELEC_PARAM%XFQSEDI
  XEXQSEDI => ELEC_PARAM%XEXQSEDI
  XFQSEDS => ELEC_PARAM%XFQSEDS
  XEXQSEDS => ELEC_PARAM%XEXQSEDS
  XFQSEDG => ELEC_PARAM%XFQSEDG
  XEXQSEDG => ELEC_PARAM%XEXQSEDG
  XFQSEDH => ELEC_PARAM%XFQSEDH
  XEXQSEDH => ELEC_PARAM%XEXQSEDH
  XEGMIN => ELEC_PARAM%XEGMIN
  XEGMAX => ELEC_PARAM%XEGMAX
  XESMIN => ELEC_PARAM%XESMIN
  XESMAX => ELEC_PARAM%XESMAX
  XEIMIN => ELEC_PARAM%XEIMIN
  XEIMAX => ELEC_PARAM%XEIMAX
  XECMIN => ELEC_PARAM%XECMIN
  XECMAX => ELEC_PARAM%XECMAX
  XERMIN => ELEC_PARAM%XERMIN
  XERMAX => ELEC_PARAM%XERMAX
  XEHMIN => ELEC_PARAM%XEHMIN
  XEHMAX => ELEC_PARAM%XEHMAX
  XFQUPDC => ELEC_PARAM%XFQUPDC
  XFQUPDR => ELEC_PARAM%XFQUPDR
  XFQUPDI => ELEC_PARAM%XFQUPDI
  XEXFQUPDI => ELEC_PARAM%XEXFQUPDI
  XFQUPDS => ELEC_PARAM%XFQUPDS
  XFQUPDG => ELEC_PARAM%XFQUPDG
  XFQUPDH => ELEC_PARAM%XFQUPDH
END SUBROUTINE ELEC_PARAM_ASSOCIATE
!
END MODULE MODD_ELEC_PARAM
