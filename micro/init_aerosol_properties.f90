!MNH_LIC Copyright 2013-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!      ####################
       MODULE MODI_INIT_AEROSOL_PROPERTIES
INTERFACE
   SUBROUTINE INIT_AEROSOL_PROPERTIES
   END SUBROUTINE INIT_AEROSOL_PROPERTIES
END INTERFACE
END MODULE MODI_INIT_AEROSOL_PROPERTIES
!      ####################
!
!     #############################################################
      SUBROUTINE INIT_AEROSOL_PROPERTIES
!     #############################################################

!!
!!
!!      PURPOSE
!!      -------
!!      
!!      Define the aerosol properties
!!
!!
!!    AUTHOR
!!    ------
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!      S.    Berthet    * Laboratoire d'Aerologie*
!!      B.    Vié        * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original             ??/??/13 
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!  Philippe Wautelet: 22/01/2019: bugs correction: incorrect writes + unauthorized goto
!  P. Wautelet 10/04/2019: replace ABORT and STOP calls by Print_msg
!  P. Wautelet 30/03/2021: move NINDICE_CCN_IMM and NIMM initializations from init_aerosol_properties to ini_nsv
!  B. Vié         06/2021: kappa-kohler CCN activation parameters
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!USE MODD_LUNIT,           ONLY : TLUOUT0
USE MODD_PARAM_LIMA,      ONLY : NMOD_CCN, HINI_CCN, HTYPE_CCN,        &
                                 XR_MEAN_CCN, XLOGSIG_CCN, XRHO_CCN,                 &
                                 XKHEN_MULTI, XMUHEN_MULTI, XBETAHEN_MULTI,          &
                                 XLIMIT_FACTOR, CCCN_MODES, LSCAV,                    &
                                 XACTEMP_CCN, XFSOLUB_CCN,                           &
                                 NMOD_IFN, NSPECIE, CIFN_SPECIES,       &
                                 XMDIAM_IFN, XSIGMA_IFN, XRHO_IFN, XFRAC, XFRAC_REF, &
                                 CINT_MIXING, NPHILLIPS,           &
                                 NIMM, NMOD_IMM, NINDICE_CCN_IMM
!
use mode_msg
!
USE MODI_GAMMA
USE MODE_LIMA_INIT_CCN_ACTIVATION_SPECTRUM, ONLY: LIMA_INIT_CCN_ACTIVATION_SPECTRUM
!
IMPLICIT NONE
!
REAL             ::  XKHEN0     
REAL             ::  XLOGSIG0   
REAL             ::  XALPHA1     
REAL             ::  XMUHEN0    
REAL             ::  XALPHA2   
REAL             ::  XBETAHEN0  
REAL             ::  XR_MEAN0  
REAL             ::  XALPHA3    
REAL             ::  XALPHA4     
REAL             ::  XALPHA5      
REAL             ::  XACTEMP0    
REAL             ::  XALPHA6   
!
REAL, DIMENSION(6) :: XKHEN_TMP   = (/1.56, 1.56, 1.56, 1.56, 1.56, 1.56 /)
REAL, DIMENSION(6) :: XMUHEN_TMP  = (/0.80, 0.80, 0.80, 0.80, 0.80, 0.80 /)
REAL, DIMENSION(6) :: XBETAHEN_TMP= (/136., 136., 136., 136., 136., 136. /)
!
REAL, DIMENSION(3) :: RCCN
REAL, DIMENSION(3) :: LOGSIGCCN
REAL, DIMENSION(3) :: RHOCCN
!
INTEGER            :: I,J,JMOD
!
!INTEGER  :: ILUOUT0 ! Logical unit number for output-listing
!INTEGER  :: IRESP   ! Return code of FM-routines
!
REAL :: X1, X2, X3, X4, X5
! REAL, DIMENSION(7) :: diameters=(/ 0.01E-6, 0.05E-6, 0.1E-6, 0.2E-6, 0.5E-6, 1.E-6, 2.E-6 /)
! REAL, DIMENSION(3) :: sigma=(/ 2., 2.5, 3. /)
! CHARACTER(LEN=7), DIMENSION(3) :: types=(/ 'NH42SO4', 'NaCl   ', '       ' /)
!REAL, DIMENSION(1) :: diameters=(/ 0.25E-6 /)
!CHARACTER(LEN=7), DIMENSION(1) :: types=(/ '       ' /)
INTEGER :: II, IJ, IK
!
!-------------------------------------------------------------------------------
!
!ILUOUT0 = TLUOUT0%NLU
!
!!!!!!!!!!!!!!!!
! CCN properties
!!!!!!!!!!!!!!!!
!
IF ( NMOD_CCN .GE. 1 ) THEN
!
   IF (.NOT.(ALLOCATED(XR_MEAN_CCN))) ALLOCATE(XR_MEAN_CCN(NMOD_CCN))
   IF (.NOT.(ALLOCATED(XLOGSIG_CCN))) ALLOCATE(XLOGSIG_CCN(NMOD_CCN))
   IF (.NOT.(ALLOCATED(XRHO_CCN)))    ALLOCATE(XRHO_CCN(NMOD_CCN))
!
   SELECT CASE (CCCN_MODES)
   CASE ('JUNGFRAU')
      RCCN(:)   = (/ 0.02E-6 , 0.058E-6 , 0.763E-6 /)
      LOGSIGCCN(:) = (/ 0.28 , 0.57 , 0.34 /)
      RHOCCN(:) = (/ 1500. , 1500. , 1500. /)
   CASE ('COPT')
      RCCN(:)   = (/ 0.125E-6 , 0.4E-6 , 1.0E-6 /)
      LOGSIGCCN(:) = (/ 0.69 , 0.41 , 0.47 /)
      RHOCCN(:) = (/ 1000. , 1000. , 1000. /)
   CASE ('CAMS')
      RCCN(:)   = (/ 0.4E-6 , 0.25E-6 , 0.1E-6 /)
      LOGSIGCCN(:) = (/ 0.64 , 0.47 , 0.47 /)
      RHOCCN(:) = (/ 2160. , 2000. , 1750. /)
   CASE ('CAMS_JPP')
! sea-salt, sulfate, hydrophilic (GADS data)
      RCCN(:)      = (/ 0.209E-6 , 0.0695E-6 , 0.0212E-6 /)
      LOGSIGCCN(:) = (/ 0.708    , 0.708     , 0.806     /)
      RHOCCN(:)    = (/ 2200.    , 1700.     , 1800.     /)
   CASE ('CAMS_ACC')
! sea-salt, sulfate, hydrophilic (GADS data)
      RCCN(:) = (/ 0.2E-6   , 0.5E-6    , 0.4E-6 /)
      LOGSIGCCN(:) = (/ 0.693    , 0.476     , 0.788  /)
      RHOCCN(:)    = (/ 2200.    , 1700.     , 1800.  /)
   CASE ('CAMS_AIT')
! sea-salt, sulfate, hydrophilic (GADS data)
      RCCN(:) = (/ 0.2E-6   , 0.05E-6   , 0.02E-6 /)
      LOGSIGCCN(:) = (/ 0.693    , 0.693     , 0.788   /)
      RHOCCN(:)    = (/ 2200.    , 1700.     , 1800.   /)
   CASE ('SIRTA')
      RCCN(:)   = (/ 0.153E-6 , 0.058E-6 , 0.763E-6 /)
      LOGSIGCCN(:) = (/ 0.846 , 0.57 , 0.34 /)
      RHOCCN(:) = (/ 1500. , 1500. , 1500. /)
   CASE ('CPS00')
      RCCN(:)   = (/ 0.0218E-6 , 0.058E-6 , 0.763E-6 /)
      LOGSIGCCN(:) = (/ 1.16 , 0.57 , 0.34 /)
      RHOCCN(:) = (/ 1500. , 1500. , 1500. /)
   CASE ('MOCAGE') ! ordre : sulfates, sels marins, BC+O
      RCCN(:)   = (/ 0.01E-6 , 0.05E-6 , 0.008E-6 /)
      LOGSIGCCN(:) = (/ 0.788 , 0.993 , 0.916 /)
      RHOCCN(:) = (/ 1000. , 2200. , 1000. /)
   CASE ('FREETROP') ! d'après Jaenicke 1993, aerosols troposphere libre, masse volumique typique
      RCCN(:)   = (/ 0.0035E-6 , 0.125E-6 , 0.26E-6 /)
      LOGSIGCCN(:) = (/ 0.645 , 0.253 , 0.425 /)
      RHOCCN(:) = (/ 1000. , 1000. , 1000. /)
   CASE DEFAULT
      call Print_msg(NVERB_FATAL,'GEN','INIT_AEROSOL_PROPERTIES','CCN_MODES must be JUNGFRAU, COPT, CAMS, CAMS_JPP,'// &
                                                                 'CAMS_ACC, CAMS_AIT, SIRTA, CPS00, MOCAGE or FREETROP')
   ENDSELECT
!
  DO I=1, MIN(NMOD_CCN,3)
    XR_MEAN_CCN(I) = RCCN(I)
    XLOGSIG_CCN(I) = LOGSIGCCN(I)
    XRHO_CCN(I)    = RHOCCN(I)
  END DO
!
  IF (NMOD_CCN .EQ. 4) THEN
! default values as coarse sea salt mode
    XR_MEAN_CCN(4) = 1.75E-6
    XLOGSIG_CCN(4) = 0.708
    XRHO_CCN(4)    = 2200.   
  END IF
!
!
!		Compute CCN spectra parameters from CCN characteristics
!
!* INPUT : XBETAHEN_TEST is in 'percent' and XBETAHEN_MULTI in 'no units', 
!  XK... and XMU... are invariant
!
  IF (.NOT.(ALLOCATED(XKHEN_MULTI)))    ALLOCATE(XKHEN_MULTI(NMOD_CCN))  
  IF (.NOT.(ALLOCATED(XMUHEN_MULTI)))   ALLOCATE(XMUHEN_MULTI(NMOD_CCN))
  IF (.NOT.(ALLOCATED(XBETAHEN_MULTI))) ALLOCATE(XBETAHEN_MULTI(NMOD_CCN))
  IF (.NOT.(ALLOCATED(XLIMIT_FACTOR)))  ALLOCATE(XLIMIT_FACTOR(NMOD_CCN))
!
  IF (HINI_CCN == 'CCN') THEN
!!$    IF (LSCAV) THEN
!!$! Attention !
!!$      WRITE(UNIT=ILUOUT0,FMT='("You are using a numerical initialization &
!!$        &not depending on the aerosol properties, however you need it for &
!!$        &scavenging.                                                      &
!!$        &With LSCAV = true, HINI_CCN should be set to AER for consistency")')
!!$    END IF
! Numerical initialization without dependence on AP physical properties
    DO JMOD = 1, NMOD_CCN
      XKHEN_MULTI(JMOD)    = XKHEN_TMP(JMOD)  
      XMUHEN_MULTI(JMOD)   = XMUHEN_TMP(JMOD)
      XBETAHEN_MULTI(JMOD) = XBETAHEN_TMP(JMOD)*(100.)**2                    
! no units relative to smax 
      XLIMIT_FACTOR(JMOD)  = ( GAMMA_X0D(0.5*XKHEN_MULTI(JMOD)+1.)&
           *GAMMA_X0D(XMUHEN_MULTI(JMOD)-0.5*XKHEN_MULTI(JMOD)) ) &
           /( XBETAHEN_MULTI(JMOD)**(0.5*XKHEN_MULTI(JMOD))       &
           *GAMMA_X0D(XMUHEN_MULTI(JMOD)) ) ! N/C
    END DO
  ELSE IF (HINI_CCN == 'AER') THEN
! 
! Initialisation depending on aerosol physical properties
!
! First, computing k, mu, beta, and XLIMIT_FACTOR as in CPS2000 (eqs 9a-9c)
!
! XLIMIT_FACTOR replaces C, because C depends on the CCN number concentration
! which is therefore determined at each grid point and time step as 
! Nccn / XLIMIT_FACTOR
!
    DO JMOD = 1, NMOD_CCN 
!
!!$       SELECT CASE (HTYPE_CCN(JMOD))
!!$       CASE ('M') ! CCN marins
!!$          XKHEN0     = 3.251
!!$          XLOGSIG0   = 0.4835
!!$          XALPHA1    = -1.297
!!$          XMUHEN0    = 2.589
!!$          XALPHA2    = -1.511
!!$          XBETAHEN0  = 621.689
!!$          XR_MEAN0   = 0.133E-6
!!$          XALPHA3    = 3.002
!!$          XALPHA4    = 1.081
!!$          XALPHA5    = 1.0
!!$          XACTEMP0   = 290.16
!!$          XALPHA6    = 2.995
!!$       CASE ('C') ! CCN continentaux
!!$          XKHEN0     = 1.403
!!$          XLOGSIG0   = 1.16
!!$          XALPHA1    = -1.172
!!$          XMUHEN0    = 0.834
!!$          XALPHA2    = -1.350
!!$          XBETAHEN0  = 25.499
!!$          XR_MEAN0   = 0.0218E-6
!!$          XALPHA3    = 3.057
!!$          XALPHA4    = 4.092
!!$          XALPHA5    = 1.011
!!$          XACTEMP0   = 290.16
!!$          XALPHA6    = 3.076
!!$       CASE DEFAULT
!!$          call Print_msg(NVERB_FATAL,'GEN','INIT_AEROSOL_PROPERTIES','HTYPE_CNN(JMOD)=C or M must be specified'// &
!!$                                                                     ' in EXSEG1.nam for each CCN mode')
!!$       ENDSELECT
!!$!
!!$      XKHEN_MULTI(JMOD)   =   XKHEN0*(XLOGSIG_CCN(JMOD)/XLOGSIG0)**XALPHA1
!!$      XMUHEN_MULTI(JMOD)  =  XMUHEN0*(XLOGSIG_CCN(JMOD)/XLOGSIG0)**XALPHA2
!!$      XBETAHEN_MULTI(JMOD)=XBETAHEN0*(XR_MEAN_CCN(JMOD)/XR_MEAN0)**XALPHA3 &
!!$           * EXP( XALPHA4*((XLOGSIG_CCN(JMOD)/XLOGSIG0)-1.) )              &
!!$           * XFSOLUB_CCN**XALPHA5                                          &
!!$           * (XACTEMP_CCN/XACTEMP0)**XALPHA6
!!$      XLIMIT_FACTOR(JMOD)  = ( GAMMA_X0D(0.5*XKHEN_MULTI(JMOD)+1.) &
!!$           *GAMMA_X0D(XMUHEN_MULTI(JMOD)-0.5*XKHEN_MULTI(JMOD)) )  &
!!$           /( XBETAHEN_MULTI(JMOD)**(0.5*XKHEN_MULTI(JMOD))        &
!!$           *GAMMA_X0D(XMUHEN_MULTI(JMOD)) )
!!$
!!$
       CALL LIMA_INIT_CCN_ACTIVATION_SPECTRUM (HTYPE_CCN(JMOD),XR_MEAN_CCN(JMOD)*2.,EXP(XLOGSIG_CCN(JMOD)),X1,X2,X3,X4,X5)
       !
       ! LIMA_INIT_CCN_ACTIVATION_SPECTRUM returns X1=C/Nccn (instead of XLIMIT_FACTOR), X2=k, X3=mu, X4=beta, X5=kappa
       ! So XLIMIT_FACTOR = 1/X1
       ! Nc = Nccn/XLIMIT_FACTOR * S^k *F() = Nccn * X1 * S^k *F()
       !
       XLIMIT_FACTOR(JMOD) = 1./X1
       XKHEN_MULTI(JMOD)   = X2
       XMUHEN_MULTI(JMOD)  = X3
       XBETAHEN_MULTI(JMOD)= X4
    ENDDO
!
! These parameters are correct for a nucleation spectra 
! Nccn(Smax) = C Smax^k F(mu,k/2,1+k/2,-beta Smax^2)
! with Smax expressed in % (Smax=1 for a supersaturation of 1%).
!
! All the computations in LIMA are done for an adimensional Smax (Smax=0.01 for 
! a 1% supersaturation). So beta and C (XLIMIT_FACTOR) are changed :
! new_beta = beta * 100^2
! new_C = C * 100^k (ie XLIMIT_FACTOR = XLIMIT_FACTOR / 100^k)
!
    XBETAHEN_MULTI(:) = XBETAHEN_MULTI(:) * 10000
    XLIMIT_FACTOR(:)  = XLIMIT_FACTOR(:) / (100**XKHEN_MULTI(:))
  END IF
END IF ! NMOD_CCN > 0
!
!!!!!!!!!!!!!!!!
! IFN properties
!!!!!!!!!!!!!!!!
!
IF ( NMOD_IFN .GE. 1 ) THEN
  SELECT CASE (CIFN_SPECIES)
   CASE ('MOCAGE')
         NSPECIE = 4
         IF (.NOT.(ALLOCATED(XMDIAM_IFN))) ALLOCATE(XMDIAM_IFN(NSPECIE))
         IF (.NOT.(ALLOCATED(XSIGMA_IFN))) ALLOCATE(XSIGMA_IFN(NSPECIE))
         IF (.NOT.(ALLOCATED(XRHO_IFN)))   ALLOCATE(XRHO_IFN(NSPECIE))
         XMDIAM_IFN = (/ 0.05E-6 , 3.E-6 , 0.016E-6 , 0.016E-6 /)
         XSIGMA_IFN = (/ 2.4 , 1.6 , 2.5 , 2.5 /)
         XRHO_IFN   = (/ 2650. , 2650. , 1000. , 1000. /)
   CASE ('CAMS_JPP')
! sea-salt, sulfate, hydrophilic (GADS data)
! 2 species, dust-metallic and hydrophobic (as BC)
! (Phillips et al. 2013 and GADS data)
      NSPECIE = 4 ! DM1, DM2, BC, BIO+(O)
      IF (.NOT.(ALLOCATED(XMDIAM_IFN))) ALLOCATE(XMDIAM_IFN(NSPECIE))
      IF (.NOT.(ALLOCATED(XSIGMA_IFN))) ALLOCATE(XSIGMA_IFN(NSPECIE))
      IF (.NOT.(ALLOCATED(XRHO_IFN)))   ALLOCATE(XRHO_IFN(NSPECIE))
      XMDIAM_IFN = (/0.8E-6, 3.0E-6, 0.025E-6, 0.2E-6/)
      XSIGMA_IFN = (/2.0, 2.15, 2.0, 1.6 /)
      XRHO_IFN   = (/2600., 2600., 1000., 1500./) 
   CASE ('CAMS_ACC')
! sea-salt, sulfate, hydrophilic (GADS data)
! 2 species, dust-metallic and hydrophobic (as BC)
! (Phillips et al. 2013 and GADS data)
      NSPECIE = 4 ! DM1, DM2, BC, BIO+(O)
      IF (.NOT.(ALLOCATED(XMDIAM_IFN))) ALLOCATE(XMDIAM_IFN(NSPECIE))
      IF (.NOT.(ALLOCATED(XSIGMA_IFN))) ALLOCATE(XSIGMA_IFN(NSPECIE))
      IF (.NOT.(ALLOCATED(XRHO_IFN)))   ALLOCATE(XRHO_IFN(NSPECIE))
      XMDIAM_IFN = (/0.8E-6, 3.0E-6, 0.04E-6, 0.8E-6 /)
      XSIGMA_IFN = (/2.0,    2.15,   2.0,     2.2    /)
      XRHO_IFN   = (/2600.,  2600.,  1000.,   2000.  /)
   CASE ('CAMS_AIT')
! sea-salt, sulfate, hydrophilic (GADS data)
! 2 species, dust-metallic and hydrophobic (as BC)
! (Phillips et al. 2013 and GADS data)
      NSPECIE = 4 ! DM1, DM2, BC, BIO+(O)
      IF (.NOT.(ALLOCATED(XMDIAM_IFN))) ALLOCATE(XMDIAM_IFN(NSPECIE))
      IF (.NOT.(ALLOCATED(XSIGMA_IFN))) ALLOCATE(XSIGMA_IFN(NSPECIE))
      IF (.NOT.(ALLOCATED(XRHO_IFN)))   ALLOCATE(XRHO_IFN(NSPECIE))
      XMDIAM_IFN = (/0.8E-6, 3.0E-6, 0.04E-6, 0.04E-6/)
      XSIGMA_IFN = (/2.0,    2.15,   2.0,     2.2 /)
      XRHO_IFN   = (/2600.,  2600.,  1000.,   1800./)
   CASE DEFAULT
      IF (NPHILLIPS == 8) THEN
! 4 species, according to Phillips et al. 2008
         NSPECIE = 4
         IF (.NOT.(ALLOCATED(XMDIAM_IFN))) ALLOCATE(XMDIAM_IFN(NSPECIE))
         IF (.NOT.(ALLOCATED(XSIGMA_IFN))) ALLOCATE(XSIGMA_IFN(NSPECIE))
         IF (.NOT.(ALLOCATED(XRHO_IFN)))   ALLOCATE(XRHO_IFN(NSPECIE))
         XMDIAM_IFN = (/0.8E-6, 3.0E-6, 0.2E-6, 0.2E-6/)
         XSIGMA_IFN = (/1.9, 1.6, 1.6, 1.6 /)
         XRHO_IFN   = (/2300., 2300., 1860., 1500./)
      ELSE IF (NPHILLIPS == 13) THEN
! 4 species, according to Phillips et al. 2013
         NSPECIE = 4
         IF (.NOT.(ALLOCATED(XMDIAM_IFN))) ALLOCATE(XMDIAM_IFN(NSPECIE))
         IF (.NOT.(ALLOCATED(XSIGMA_IFN))) ALLOCATE(XSIGMA_IFN(NSPECIE))
         IF (.NOT.(ALLOCATED(XRHO_IFN)))   ALLOCATE(XRHO_IFN(NSPECIE))
         XMDIAM_IFN = (/0.8E-6, 3.0E-6, 90.E-9, 0.163E-6/)
         XSIGMA_IFN = (/1.9, 1.6, 1.6, 2.54 /)
         XRHO_IFN   = (/2300., 2300., 1860., 1000./)
      END IF
   ENDSELECT
!
! internal mixing
!
   IF (.NOT.(ALLOCATED(XFRAC))) ALLOCATE(XFRAC(NSPECIE,NMOD_IFN))
   XFRAC(:,:)=0. 
   SELECT CASE (CINT_MIXING)
   CASE ('DM1')
      XFRAC(1,:)=1.
   CASE ('DM2')
      XFRAC(2,:)=1.
   CASE ('BC')
      XFRAC(3,:)=1.
   CASE ('O')
      XFRAC(4,:)=1.
   CASE ('CAMS')
      XFRAC(1,1)=0.99
      XFRAC(2,1)=0.01
      XFRAC(3,1)=0.
      XFRAC(4,1)=0.
      XFRAC(1,2)=0.
      XFRAC(2,2)=0.
      XFRAC(3,2)=0.5
      XFRAC(4,2)=0.5
   CASE ('CAMS_JPP')
      XFRAC(1,1)=1.0
      XFRAC(2,1)=0.0
      XFRAC(3,1)=0.0
      XFRAC(4,1)=0.0
      XFRAC(1,2)=0.0
      XFRAC(2,2)=0.0
      XFRAC(3,2)=0.5
      XFRAC(4,2)=0.5
   CASE ('CAMS_ACC')
      XFRAC(1,1)=1.0
      XFRAC(2,1)=0.0
      XFRAC(3,1)=0.0
      XFRAC(4,1)=0.0
      XFRAC(1,2)=0.0
      XFRAC(2,2)=0.0
      XFRAC(3,2)=0.0
      XFRAC(4,2)=1.0
   CASE ('CAMS_AIT')
      XFRAC(1,1)=1.0
      XFRAC(2,1)=0.0
      XFRAC(3,1)=0.0
      XFRAC(4,1)=0.0
      XFRAC(1,2)=0.0
      XFRAC(2,2)=0.0
      XFRAC(3,2)=0.0
      XFRAC(4,2)=1.0
   CASE ('MOCAGE')
      XFRAC(1,1)=1.
      XFRAC(2,1)=0.
      XFRAC(3,1)=0.
      XFRAC(4,1)=0.
      XFRAC(1,2)=0.
      XFRAC(2,2)=0.
      XFRAC(3,2)=0.7
      XFRAC(4,2)=0.3
   CASE DEFAULT
      XFRAC(1,:)=0.6
      XFRAC(2,:)=0.009
      XFRAC(3,:)=0.33
      XFRAC(4,:)=0.06
   ENDSELECT
!
! Phillips 08 alpha (table 1) 
   IF (.NOT.(ALLOCATED(XFRAC_REF))) ALLOCATE(XFRAC_REF(4))
      IF (NPHILLIPS == 13) THEN
         XFRAC_REF(1)=0.66
         XFRAC_REF(2)=0.66
         XFRAC_REF(3)=0.31
         XFRAC_REF(4)=0.03
      ELSE IF (NPHILLIPS == 8) THEN
         XFRAC_REF(1)=0.66
         XFRAC_REF(2)=0.66
         XFRAC_REF(3)=0.28
         XFRAC_REF(4)=0.06
      END IF
!
! Immersion modes
!
   IF (.NOT.(ALLOCATED(NIMM))) ALLOCATE(NIMM(NMOD_CCN))
   NIMM(:)=0
   IF (ALLOCATED(NINDICE_CCN_IMM)) DEALLOCATE(NINDICE_CCN_IMM)
   ALLOCATE(NINDICE_CCN_IMM(MAX(1,NMOD_IMM)))
   IF (NMOD_IMM .GE. 1) THEN
      DO J = 0, NMOD_IMM-1
         NIMM(NMOD_CCN-J)=1
         NINDICE_CCN_IMM(NMOD_IMM-J) = NMOD_CCN-J
      END DO
!   ELSE IF (NMOD_IMM == 0) THEN ! PNIS existe mais vaut 0, pour l'appel à resolved_cloud
!      NMOD_IMM = 1
!      NINDICE_CCN_IMM(1) = 0
   END IF
!
END IF ! NMOD_IFN > 0
! 
END SUBROUTINE INIT_AEROSOL_PROPERTIES
