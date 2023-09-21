!MNH_LIC Copyright 2007-2020 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!      #################
       MODULE MODI_LIDAR 
!      #################
!
INTERFACE
      SUBROUTINE LIDAR(HCLOUD,HVIEW,PALT,PWVL,PZZ,PRHO,PT,PCLDFR,PRT, &
                       PLIDAROUT,PLIPAROUT,PCT,PDSTC,PDSTD,PDSTS)
!
CHARACTER(LEN=*),         INTENT(IN) :: HCLOUD  ! Name of the cloud scheme
CHARACTER(LEN=*),         INTENT(IN) :: HVIEW   ! Upward or Downward integration
REAL,                     INTENT(IN) :: PALT    ! Altitude of the lidar source
REAL,                     INTENT(IN) :: PWVL    ! Wavelength of the lidar source
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PZZ     ! Altitude
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PRHO    ! Air density
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PT      ! Air temperature
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PCLDFR  ! Cloud fraction
REAL, DIMENSION(:,:,:,:), INTENT(IN) :: PRT     ! Moist variables at t
REAL, DIMENSION(:,:,:),  INTENT(OUT) :: PLIDAROUT ! Lidar output
REAL, DIMENSION(:,:,:),  INTENT(OUT) :: PLIPAROUT ! Lidar output (particle only)

REAL, DIMENSION(:,:,:,:), OPTIONAL, INTENT(IN) :: PCT ! Concentration 
                                                      ! (C2R2 and C1R3) 
REAL, DIMENSION(:,:,:,:), OPTIONAL, INTENT(IN) :: PDSTC ! Dust Concentration 
REAL, DIMENSION(:,:,:,:), OPTIONAL, INTENT(IN) :: PDSTD ! Dust Diameter
REAL, DIMENSION(:,:,:,:), OPTIONAL, INTENT(IN) :: PDSTS ! Dust Sigma
!

!
END SUBROUTINE LIDAR
!
END INTERFACE
!
END MODULE MODI_LIDAR
!     #########################################################
      SUBROUTINE LIDAR(HCLOUD,HVIEW,PALT,PWVL,PZZ,PRHO,PT,PCLDFR,PRT, &
                       PLIDAROUT,PLIPAROUT,PCT,PDSTC,PDSTD,PDSTS)
!     #########################################################
!
!!****  *LIDAR * - computes pertinent lidar parameters
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the normalized backscattered
!!    signal of an upward or downward looking lidar in an atmosperic column
!!    containing air molecules, aerosols, cloud particles and hydrometeors.
!!
!!**  METHOD
!!    ------
!!      The reflectivities are computed using the n(D) * D**6 formula. The 
!!    equivalent reflectiviy is the sum of the reflectivity produced by the
!!    the raindrops and the equivalent reflectivities of the ice crystals.
!!    The latter are computed using the melted diameter. The Doppler 
!!    reflectivity is the 'fall speed'-moment of individual particle
!!    reflectivity. Ice crystal are assumed to have no preferred orientation.
!!    the Z_VV formula is taken from Brandes et al. (MWR, 1995).
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST
!!        XPI                  !
!!        XRHOLW               ! Liquid water density
!!      Module MODD_RAIN_ICE_DESCR
!!      Module MODD_RAIN_ICE_PARAM
!!
!!    REFERENCE
!!    ---------
!!    Chaboureau et al. 2011: Long-range transport of Saharan dust and its
!!    radiative impact on precipitation forecast over western Europe: a case
!!    study during COPS. Quart. J. Roy. Meteor. Soc., 137, 236-251
!!
!!    AUTHOR
!!    ------
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original      04/10/07
!!      JP Chaboureau 12/02/10 change dust refraction index
!!                             add inputs (lidar charact. and cloud fraction)
!!   J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1
!!   B.VIE  2016 : LIMA
!  P. Wautelet 18/03/2020: remove ICE2 option
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAMETERS
USE MODD_CST
USE MODD_RAIN_C2R2_DESCR, ONLY : XLBEXC, XLBEXR, &
                                 XRTMIN, XCTMIN
USE MODD_PARAM_C2R2,      ONLY : YALPHAC=>XALPHAC,YNUC=>XNUC, &
                                 YALPHAR=>XALPHAR,YNUR=>XNUR
USE MODD_PARAM_ICE_n,        ONLY: WSNOW_T=>LSNOW_T
USE MODD_RAIN_ICE_DESCR_n,  ONLY : XCCR, WLBEXR=>XLBEXR, XLBR, &
                                 XCCS, XCXS,   XLBEXS, XLBS, WNS=>XNS, WBS=>XBS, &
                                 XCCG, XCXG,   XLBEXG, XLBG, &
                                 XCCH, XCXH,   XLBEXH, XLBH, &
                                 WRTMIN=>XRTMIN,   &
                                 WLBDAS_MAX=>XLBDAS_MAX,WLBDAS_MIN=>XLBDAS_MIN,WTRANS_MP_GAMMAS=>XTRANS_MP_GAMMAS
USE MODD_ICE_C1R3_DESCR,  ONLY : XLBEXI,                      &
                                 YRTMIN=>XRTMIN, YCTMIN=>XCTMIN
!
USE MODD_PARAM_LIMA,      ONLY : URTMIN=>XRTMIN, UCTMIN=>XCTMIN, &
                                 UALPHAC=>XALPHAC,UNUC=>XNUC, &
                                 UALPHAR=>XALPHAR,UNUR=>XNUR, &
                                 UALPHAI=>XALPHAI,UNUI=>XNUI, &
                                 USNOW_T=>LSNOW_T
USE MODD_PARAM_LIMA_COLD, ONLY : UCCS=>XCCS, UCXS=>XCXS, ULBEXS=>XLBEXS, & 
                                                ULBS=>XLBS, UNS=>XNS, UBS=>XBS,   &
                                 ULBDAS_MAX=>XLBDAS_MAX,ULBDAS_MIN=>XLBDAS_MIN,UTRANS_MP_GAMMAS=>XTRANS_MP_GAMMAS
USE MODD_PARAM_LIMA_MIXED,ONLY : UCCG=>XCCG, UCXG=>XCXG, ULBEXG=>XLBEXG, &
                                                         ULBG=>XLBG

use mode_tools_ll,        only: GET_INDICE_ll

USE MODI_BHMIE_WATER    ! Gamma or mono dispersed size distributions
USE MODI_BHMIE_AEROSOLS ! Lognormal or mono dispersed size distributions
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
CHARACTER(LEN=*),         INTENT(IN) :: HCLOUD  ! Name of the cloud scheme
CHARACTER(LEN=*),         INTENT(IN) :: HVIEW   ! Upward or Downward integration
REAL,                     INTENT(IN) :: PALT    ! Altitude of the lidar source
REAL,                     INTENT(IN) :: PWVL    ! Wavelength of the lidar source
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PZZ     ! Altitude
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PRHO    ! Air density
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PT      ! Air temperature
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PCLDFR  ! Cloud fraction
REAL, DIMENSION(:,:,:,:), INTENT(IN) :: PRT     ! Moist variables at t
REAL, DIMENSION(:,:,:),  INTENT(OUT) :: PLIDAROUT ! Lidar output
REAL, DIMENSION(:,:,:),  INTENT(OUT) :: PLIPAROUT ! Lidar output (particle only)

REAL, DIMENSION(:,:,:,:), OPTIONAL, INTENT(IN) :: PCT  ! Concentration 
                                                       ! (C2R2 and C1R3) 
REAL, DIMENSION(:,:,:,:), OPTIONAL, INTENT(IN) :: PDSTC ! Dust Concentration 
REAL, DIMENSION(:,:,:,:), OPTIONAL, INTENT(IN) :: PDSTD ! Dust Diameter
REAL, DIMENSION(:,:,:,:), OPTIONAL, INTENT(IN) :: PDSTS ! Dust Sigma
!
!*       0.2   Declarations of local variables :
!
INTEGER :: JI, JJ, JK
INTEGER :: IIB, IIE, IJB, IJE, IKB, IKE
INTEGER :: IALT
INTEGER, DIMENSION(3) :: IKMIN
!
REAL, PARAMETER    :: ZCLDFRMIN = 1.0E-03                  ! Cloud fraction minimum
!
!
COMPLEX, PARAMETER :: ZZREFIND_WAT = (1.337E+00,1.818E-09) ! Refraction Index
                                                           ! of pure water
COMPLEX, PARAMETER :: ZZREFIND_ICE = (1.312E+00,2.614E-09) ! Refraction Index
                                                           ! of pure ice
!COMPLEX, PARAMETER :: ZZREFIND_DUST= (1.530E+00,8.000E-03) ! Refraction Index
!                                                           ! of mineral dust
! West, R. A., L. R. Doose, A. M. Eibl, M. G. Tomasko, and M. I. Mishchenko
! (1997), Laboratory measurements of mineral dust scattering phase function 
! and linear polarization, J. Geophys. Res., 102(D14), 16,871-16,882.
COMPLEX :: ZZREFIND_DUST

! Tulet, P., M. Mallet, V. Pont, J. Pelon, and A. Boone (2008), The 7-13
! March 2006 dust storm over West Africa: Generation, transport, and vertical
! stratification, J. Geophys. Res., 113, D00C08, doi:10.1029/2008JD009871.
!! Ri = 1.448-0.00292i for wavelengths between 0.185 and 0.69um.
!! Ri = 1.44023-0.00116i for wavelengths between 0.69 and 1.19um.
!! Ri = 1.41163-0.00106i for wavelengths between 1.19 and 4.0um.
COMPLEX, PARAMETER :: ZZREFIND_DSTL= (1.448,2.92E-03)
COMPLEX, PARAMETER :: ZZREFIND_DSTM= (1.44023,1.16E-03)
COMPLEX, PARAMETER :: ZZREFIND_DSTH= (1.41163,1.06E-03)

!
! COMPLEX, PARAMETER :: ZZREFIND_WAT = (1.321E+00,1.280E-06) ! Refraction Index
!                                                            ! of pure water
! COMPLEX, PARAMETER :: ZZREFIND_ICE = (1.300E+00,1.898E-06) ! Refraction Index
!                                                            ! of pure ice
!
COMPLEX, PARAMETER :: ZZREFIND_COAT= (1.337E+00,1.818E-09) ! Refraction Index
                                                           ! of coating material
COMPLEX, PARAMETER :: ZZREFIND_BC  = (1.870E+00,0.569E+00) ! Refraction Index
                                                           ! of black carbone
REAL               :: ZCXR=-1.0                  ! for rain N ~ 1/N_0 
                                                 ! (in Kessler parameterization)
!
REAL :: ZCMOL
REAL :: ZWAVE_LENGTH
! BETA: backscattering coefficient
! ALPHA:    extinction coefficient
REAL,    DIMENSION(SIZE(PRHO,1),SIZE(PRHO,2),SIZE(PRHO,3))  :: ZBETA_MOL
REAL,    DIMENSION(SIZE(PRHO,1),SIZE(PRHO,2),SIZE(PRHO,3))  :: ZALPH_MOL
REAL,    DIMENSION(SIZE(PRHO,1),SIZE(PRHO,2),SIZE(PRHO,3))  :: ZBETA_PAR
REAL,    DIMENSION(SIZE(PRHO,1),SIZE(PRHO,2),SIZE(PRHO,3))  :: ZALPH_PAR
REAL, ALLOCATABLE, DIMENSION(:,:) :: ZOPTD_TOT ! Optical depths
REAL, ALLOCATABLE, DIMENSION(:,:) :: ZOPTD_MOL  
REAL, ALLOCATABLE, DIMENSION(:,:) :: ZOPTD_PAR 
! 
CHARACTER (LEN=5) :: YDSD
INTEGER :: IRADIUS, IANGLE
REAL    :: ZRADIUS, ZCONC, ZLWC, ZIWC
REAL    :: ZREFF_FACT
REAL    :: ZEXT_COEF, ZBAK_COEF
REAL    :: ZLBDAR, ZLBDAS, ZLBDAG, ZLBDAH, ZTCEL
REAL    :: ZFRACVOL_CORE, ZDMODAL, ZSIG
REAL    :: ZFRACVOL_BC
!
REAL    :: ZETACLD, ZETAAER       ! Multiple diffusion paramter for cloud and dust
!
REAL, DIMENSION(5) :: ZPOLC, ZPOLR, ZPOLI ! BackScat. Coefficients
! 
REAL, DIMENSION(10) :: ZRTMIN, ZCTMIN
REAL                :: ZLBEXR
!
INTEGER :: JL
REAL :: ZALPHAC, ZNUC, ZALPHAR, ZNUR, ZALPHAI, ZNUI
REAL :: ZCCS, ZCXS, ZLBEXS, ZLBS, ZNS
REAL :: ZCCG, ZCXG, ZLBEXG, ZLBG
! 
! -----------------------------------------------------------------------------
!
!*       1.     COMPUTE THE LOOP BOUNDS
!   	        -----------------------
!
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
IKB=1+JPVEXT
IKE=SIZE(PRHO,3) - JPVEXT
!
ZWAVE_LENGTH = PWVL
ZCMOL=5.45E-32*((ZWAVE_LENGTH)/0.55E-6)**(-4.09)
IF (ZWAVE_LENGTH<0.69E-6) THEN
  PRINT *,'Tulet et al. refractive index - low wavelength'
  ZZREFIND_DUST = ZZREFIND_DSTL
ELSEIF (ZWAVE_LENGTH<1.00E-6) THEN
  PRINT *,'Tulet et al. refractive index - medium wavelength'
  ZZREFIND_DUST = ZZREFIND_DSTM
ELSE
  PRINT *,'Tulet et al. refractive index - high wavelength'
  ZZREFIND_DUST = ZZREFIND_DSTH
END IF
ZPOLC = (/ 2.6980E-8,-3.7701E-6, 1.6594E-4,-0.0024, 0.0626 /)
ZPOLR(:) = ZPOLC(:) 
ZPOLI = (/-1.0176E-8, 1.7615E-6,-1.0480E-4, 0.0019, 0.0460 /)
!
! Multiple diffusion parameter
ZETAAER=1.0
ZETACLD=1.0
! a multiple scattering correction for lidar in space; Platt 73
IF (HVIEW=='NADIR'.AND.PALT==0.) ZETACLD=0.5
PRINT *,'Multiple diffusion parameter for aerosol ',ZETAAER
PRINT *,'Multiple diffusion parameter for cloud   ',ZETACLD
!
!
!*       1.     MORE INITIALIZATION
!               -------------------
!
SELECT CASE ( HCLOUD )
  CASE('KESS')
    ZRTMIN(1) = 1.0E-20
    ZRTMIN(2) = 1.0E-20
    ZRTMIN(3) = 1.0E-20
    ZLBEXR = 1.0/(-1.0-3.0)
  CASE('ICE3','ICE4')
    ZRTMIN(1:SIZE(WRTMIN)) = WRTMIN(1:SIZE(WRTMIN))
    ZLBEXR = WLBEXR
    ZCCS    = XCCS
    ZCXS    = XCXS
    ZLBEXS  = XLBEXS
    ZLBS    = XLBS
    ZNS     = WNS
    ZCCG    = XCCG
    ZCXG    = XCXG
    ZLBEXG  = XLBEXG
    ZLBG    = XLBG
  CASE('C2R2')
    ZRTMIN(1:SIZE(XRTMIN)) = XRTMIN(1:SIZE(XRTMIN))
    ZCTMIN(1:SIZE(XCTMIN)) = XCTMIN(1:SIZE(XCTMIN))
    ZLBEXR  = XLBEXR
    ZALPHAC = YALPHAC
    ZNUR    = YNUR
    ZALPHAR = YALPHAR
    ZNUC    = YNUC
  CASE('C3R5')
    ZRTMIN(1:SIZE(YRTMIN)) = YRTMIN(1:SIZE(YRTMIN))
    ZCTMIN(1:SIZE(YCTMIN)) = YCTMIN(1:SIZE(YCTMIN))
    ZALPHAC = YALPHAC
    ZNUR    = YNUR
    ZALPHAR = YALPHAR
    ZNUC    = YNUC
    ZALPHAI = ZALPHAC
    ZNUI    = ZNUC
    ZCCS    = XCCS
    ZCXS    = XCXS
    ZLBEXS  = XLBEXS
    ZLBS    = XLBS
    ZCCG    = XCCG
    ZCXG    = XCXG
    ZLBEXG  = XLBEXG
    ZLBG    = XLBG
  CASE('LIMA')
      ZRTMIN(1:SIZE(URTMIN)) = URTMIN(1:SIZE(URTMIN))
      ZCTMIN(1:SIZE(UCTMIN)) = UCTMIN(1:SIZE(UCTMIN))
      ZALPHAC = UALPHAC
      ZNUR    = UNUR
      ZALPHAR = UALPHAR
      ZNUC    = UNUC
      ZALPHAI = UALPHAI
      ZNUI    = UNUI
      ZCCS    = UCCS
      ZCXS    = UCXS 
      ZLBEXS  = ULBEXS 
      ZLBS    = ULBS
      ZNS     = UNS
      ZCCG    = UCCG
      ZCXG    = UCXG 
      ZLBEXG  = ULBEXG
      ZLBG    = ULBG
END SELECT
!
! -----------------------------------------------------------------------------
!
!*       2.    INITIALIZES THE MEAN-LAYER VARIABLES
!              ------------------------------------
!
!
! MOLECULAR CONTRIBUTION
!
ZBETA_MOL(:,:,:) = ( PRHO(:,:,:)*XAVOGADRO/XMD )*ZCMOL         
ZALPH_MOL(:,:,:) = ZBETA_MOL(:,:,:)*(8.0*XPI/3.0)
!
! PARTICULAR CONTRIBUTION
!
ZBETA_PAR(:,:,:) = 0.
ZALPH_PAR(:,:,:) = 0.
!
! AEROSOL CONTRIBUTION     ! call bhmie_aerosols
!
IF (PRESENT(PDSTC)) THEN
  DO JL = 1, SIZE(PDSTD,4)
    DO JK = IKB, IKE
      DO JJ = IJB, IJE
        DO JI = IIB, IIE
          IF ( PDSTD(JI,JJ,JK,JL)>0.1 ) THEN
 !
 ! Desert dust particles
 !
            YDSD    = 'MONOD'
            ZCONC   = PDSTC(JI,JJ,JK,JL)
            ZFRACVOL_CORE = 1.0
            ZRADIUS = PDSTD(JI,JJ,JK,JL)*1.0E-6
            IF( ZRADIUS .GE. 1.0E-3 ) ZRADIUS = ZRADIUS * 1.0E-6
            CALL BHMIE_AEROSOLS( ZWAVE_LENGTH, ZZREFIND_DUST, ZZREFIND_DUST,  &
                                 YDSD, ZCONC, ZFRACVOL_CORE, ZEXT_COEF,       &
                                 ZBAK_COEF, PRADIUS=ZRADIUS )
            ZALPH_PAR(JI,JJ,JK) = ZALPH_PAR(JI,JJ,JK) + ZETAAER * ZEXT_COEF
            ZBETA_PAR(JI,JJ,JK) = ZBETA_PAR(JI,JJ,JK) + ZBAK_COEF
          END IF
        END DO
      END DO
    END DO
  END DO
END IF
!
!
! HYDROMETEOR CONTRIBUTION ! call bhmie_water
!
! LIQUID WATER
!
! Some Prefactors:  Assume Martin et al. (1994, JAS) for Reff
!
ZREFF_FACT = 1.0E-3*(3.E3/(4.0*XPI*0.67E-3))**0.33 ! Continental N=500
ZREFF_FACT = 1.0E-3*(3.E3/(4.0*XPI*0.80E-3))**0.33 ! Maritime    N=150
!
SELECT CASE ( HCLOUD )
  CASE('KESS','ICE3','ICE4')
    DO JK = IKB, IKE
      DO JJ = IJB, IJE
        DO JI = IIB, IIE
          IF ( PRT(JI,JJ,JK,2)>ZRTMIN(2) .AND. PCLDFR(JI,JJ,JK)>ZCLDFRMIN) THEN
!
! Cloud droplets
!
            YDSD = 'MONOD'
            ZCONC   = 200.E6 ! Continental case
            ZLWC    = PRHO(JI,JJ,JK)*PRT(JI,JJ,JK,2) / PCLDFR(JI,JJ,JK)
            ZRADIUS = MIN( 16.0E-6,MAX( 4.0E-6,ZREFF_FACT*(ZLWC/ZCONC)**0.33 ) )
            IANGLE  = 11
            CALL BHMIE_WATER( ZWAVE_LENGTH, ZZREFIND_WAT, YDSD, ZCONC,      &
                              IANGLE, ZEXT_COEF, ZBAK_COEF, PRADIUS=ZRADIUS )
            ZALPH_PAR(JI,JJ,JK) = ZALPH_PAR(JI,JJ,JK) + ZETACLD * ZEXT_COEF &
                                                        * PCLDFR(JI,JJ,JK)
            ZBETA_PAR(JI,JJ,JK) = ZBETA_PAR(JI,JJ,JK) + ZBAK_COEF           &
                                                        * PCLDFR(JI,JJ,JK)
          END IF
        END DO
      END DO
    END DO
    DO JK = IKB, IKE
      DO JJ = IJB, IJE
        DO JI = IIB, IIE
          IF ( PRT(JI,JJ,JK,3)>ZRTMIN(3) ) THEN
!
! Rain drops
!
            YDSD = 'MONOD'
            ZLWC    = PRHO(JI,JJ,JK)*PRT(JI,JJ,JK,3)
            ZLBDAR  = XLBR*(ZLWC)**ZLBEXR
            ZCONC   = XCCR*(ZLBDAR)**ZCXR
            ZRADIUS = 0.5*(3.0/ZLBDAR) ! Assume Marshall-Palmer law for Reff
            IANGLE  = 11
            CALL BHMIE_WATER( ZWAVE_LENGTH, ZZREFIND_WAT, YDSD, ZCONC,      &
                              IANGLE, ZEXT_COEF, ZBAK_COEF, PRADIUS=ZRADIUS )
            ZALPH_PAR(JI,JJ,JK) = ZALPH_PAR(JI,JJ,JK) + ZETACLD * ZEXT_COEF
            ZBETA_PAR(JI,JJ,JK) = ZBETA_PAR(JI,JJ,JK) + ZBAK_COEF
          END IF
        END DO
      END DO
    END DO
  CASE ('C2R2','C3R5','LIMA')
    DO JK = IKB, IKE
      DO JJ = IJB, IJE
        DO JI = IIB, IIE
          IF (PRT(JI,JJ,JK,2)>ZRTMIN(2) .AND. PCT(JI,JJ,JK,2)>ZCTMIN(2)) THEN
!
! Cloud droplets
!
            YDSD = 'GAMMA'
            ZCONC   = PCT(JI,JJ,JK,2)
            IRADIUS  = 20
            ZLWC    = PRHO(JI,JJ,JK)*PRT(JI,JJ,JK,2)
            IANGLE  = 11
            CALL BHMIE_WATER( ZWAVE_LENGTH, ZZREFIND_WAT, YDSD, ZCONC,       &
                              IANGLE, ZEXT_COEF, ZBAK_COEF, KRADIUS=IRADIUS, & 
                              PALPHA=ZALPHAC, PNU=ZNUC, PLWC=ZLWC            )
            ZALPH_PAR(JI,JJ,JK) = ZALPH_PAR(JI,JJ,JK) + ZETACLD * ZEXT_COEF
            ZBETA_PAR(JI,JJ,JK) = ZBETA_PAR(JI,JJ,JK) + ZBAK_COEF
          END IF
        END DO
      END DO
    END DO
    DO JK = IKB, IKE
      DO JJ = IJB, IJE
        DO JI = IIB, IIE
          IF (PRT(JI,JJ,JK,3)>ZRTMIN(3) .AND. PCT(JI,JJ,JK,3)>ZCTMIN(3)) THEN
!
! Rain drops
!
            YDSD = 'GAMMA'
            ZCONC   = PCT(JI,JJ,JK,3)
            IRADIUS  = 20
            ZLWC    = PRHO(JI,JJ,JK)*PRT(JI,JJ,JK,3)
            IANGLE  = 11
            CALL BHMIE_WATER( ZWAVE_LENGTH, ZZREFIND_WAT, YDSD, ZCONC,       &
                              IANGLE, ZEXT_COEF, ZBAK_COEF, KRADIUS=IRADIUS, & 
                              PALPHA=ZALPHAR, PNU=ZNUR, PLWC=ZLWC            )
            ZALPH_PAR(JI,JJ,JK) = ZALPH_PAR(JI,JJ,JK) + ZETACLD * ZEXT_COEF
            ZBETA_PAR(JI,JJ,JK) = ZBETA_PAR(JI,JJ,JK) + ZBAK_COEF
          END IF
        END DO
      END DO
    END DO
END SELECT
!
! SOLID ICE
!
SELECT CASE ( HCLOUD )
  CASE('ICE3','ICE4')
    DO JK = IKB, IKE
      DO JJ = IJB, IJE
        DO JI = IIB, IIE
          IF ( PRT(JI,JJ,JK,4)>ZRTMIN(4) .AND. PCLDFR(JI,JJ,JK)>ZCLDFRMIN) THEN
!
! Pristine crystals
!
            YDSD = 'MONOD'
            ZCONC   = 10.E3 ! Continental case
            ZIWC    = PRHO(JI,JJ,JK)*PRT(JI,JJ,JK,4) / PCLDFR(JI,JJ,JK)
            ZTCEL   = 10.0-0.0065*PZZ(JI,JJ,JK)  ! A rough estimate
            ZRADIUS = MIN( 350.0E-6,MAX( 45.0E-6,0.5E-6*(1.2351+0.0105*ZTCEL)* &
                                (5.8966*(ZIWC*1.0E3)**0.2214 +                 &
                                (0.7957*(ZIWC*1.0E3)**0.2535)*(ZTCEL+190.0)) ) )
            IANGLE  = 11
            CALL BHMIE_WATER( ZWAVE_LENGTH, ZZREFIND_ICE, YDSD, ZCONC,      &
                              IANGLE, ZEXT_COEF, ZBAK_COEF, PRADIUS=ZRADIUS )
            ZALPH_PAR(JI,JJ,JK) = ZALPH_PAR(JI,JJ,JK) + ZETACLD * ZEXT_COEF &
                                                        *PCLDFR(JI,JJ,JK)
            ZBETA_PAR(JI,JJ,JK) = ZBETA_PAR(JI,JJ,JK) + ZBAK_COEF           &
                                                        *PCLDFR(JI,JJ,JK)
          END IF
        END DO
      END DO
    END DO
  CASE ('C3R5','LIMA')
    DO JK = IKB, IKE
      DO JJ = IJB, IJE
        DO JI = IIB, IIE
          IF (PRT(JI,JJ,JK,4)>ZRTMIN(4) .AND. PCT(JI,JJ,JK,4)>ZCTMIN(4)) THEN
!
! Pristine crystals
!
            YDSD = 'GAMMA'
            ZCONC   = PCT(JI,JJ,JK,4)
            IRADIUS  = 20
            ZIWC    = PRHO(JI,JJ,JK)*PRT(JI,JJ,JK,4)
            IANGLE  = 11
            CALL BHMIE_WATER( ZWAVE_LENGTH, ZZREFIND_ICE, YDSD, ZCONC,       &
                              IANGLE, ZEXT_COEF, ZBAK_COEF, KRADIUS=IRADIUS, &
                              PALPHA=ZALPHAI, PNU=ZNUI, PLWC=ZIWC            )
            ZALPH_PAR(JI,JJ,JK) = ZALPH_PAR(JI,JJ,JK) + ZETACLD * ZEXT_COEF
            ZBETA_PAR(JI,JJ,JK) = ZBETA_PAR(JI,JJ,JK) + ZBAK_COEF
          END IF
        END DO
      END DO
    END DO
END SELECT
SELECT CASE ( HCLOUD )
  CASE('ICE3','ICE4','C3R5','LIMA')
    DO JK = IKB, IKE
      DO JJ = IJB, IJE
        DO JI = IIB, IIE
          IF ( PRT(JI,JJ,JK,5)>ZRTMIN(5) ) THEN
!
! Snow flakes
!
            YDSD = 'MONOD'
            ZIWC    = PRHO(JI,JJ,JK)*PRT(JI,JJ,JK,5)
            IF (HCLOUD=='LIMA' .AND. USNOW_T) THEN
               IF (PT(JI,JJ,JK)>263.15) THEN
                  ZLBDAS = MAX(MIN(ULBDAS_MAX, 10**(14.554-0.0423*PT(JI,JJ,JK))),ULBDAS_MIN)*UTRANS_MP_GAMMAS
               ELSE
                  ZLBDAS = MAX(MIN(ULBDAS_MAX, 10**(6.226-0.0106*PT(JI,JJ,JK))),ULBDAS_MIN)*UTRANS_MP_GAMMAS
               END IF
               ZCONC=ZNS*ZIWC*ZLBDAS**UBS
            ELSE IF (HCLOUD=='ICE3' .AND. WSNOW_T) THEN
               IF (PT(JI,JJ,JK)>263.15) THEN
                  ZLBDAS = MAX(MIN(WLBDAS_MAX, 10**(14.554-0.0423*PT(JI,JJ,JK))),WLBDAS_MIN)*WTRANS_MP_GAMMAS
               ELSE
                  ZLBDAS = MAX(MIN(WLBDAS_MAX, 10**(6.226-0.0106*PT(JI,JJ,JK))),WLBDAS_MIN)*WTRANS_MP_GAMMAS
               END IF
               ZCONC=ZNS*ZIWC*ZLBDAS**WBS
            ELSE
               ZLBDAS  = ZLBS*(ZIWC)**ZLBEXS
               ZCONC   = ZCCS*(ZLBDAS)**ZCXS
            END IF
            IF (ZLBDAS .GT. 0) THEN
              ZRADIUS = 0.5*(3.0/ZLBDAS) ! Assume Marshall-Palmer law for Reff
              IANGLE  = 11
              CALL BHMIE_WATER( ZWAVE_LENGTH, ZZREFIND_ICE, YDSD, ZCONC,      &
                                IANGLE, ZEXT_COEF, ZBAK_COEF, PRADIUS=ZRADIUS )
              ZALPH_PAR(JI,JJ,JK) = ZALPH_PAR(JI,JJ,JK) + ZETACLD * ZEXT_COEF
              ZBETA_PAR(JI,JJ,JK) = ZBETA_PAR(JI,JJ,JK) + ZBAK_COEF
            END IF
          END IF
        END DO
      END DO
    END DO
END SELECT
SELECT CASE ( HCLOUD )
  CASE('ICE3','ICE4','C3R5','LIMA')
    DO JK = IKB, IKE
      DO JJ = IJB, IJE
        DO JI = IIB, IIE
          IF ( PRT(JI,JJ,JK,6)>ZRTMIN(6) ) THEN
!
! Graupel particles
!
            YDSD = 'MONOD'
            ZIWC    = PRHO(JI,JJ,JK)*PRT(JI,JJ,JK,6)
            ZLBDAG  = ZLBG*(ZIWC)**ZLBEXG
            ZCONC   = ZCCG*(ZLBDAG)**ZCXG
            ZRADIUS = 0.5*(3.0/ZLBDAG) ! Assume Marshall-Palmer law for Reff
            IANGLE  = 11
            CALL BHMIE_WATER( ZWAVE_LENGTH, ZZREFIND_ICE, YDSD, ZCONC,      &
                              IANGLE, ZEXT_COEF, ZBAK_COEF, PRADIUS=ZRADIUS )
            ZALPH_PAR(JI,JJ,JK) = ZALPH_PAR(JI,JJ,JK) + ZETACLD * ZEXT_COEF
            ZBETA_PAR(JI,JJ,JK) = ZBETA_PAR(JI,JJ,JK) + ZBAK_COEF
          END IF
        END DO
      END DO
    END DO
END SELECT
SELECT CASE ( HCLOUD )
  CASE('ICE4')
    DO JK = IKB, IKE
      DO JJ = IJB, IJE
        DO JI = IIB, IIE
          IF ( PRT(JI,JJ,JK,7)>ZRTMIN(7) ) THEN
!
! Hailstones
!
            YDSD = 'MONOD'
            ZIWC    = PRHO(JI,JJ,JK)*PRT(JI,JJ,JK,7)
            ZLBDAH  = XLBH*(ZIWC)**XLBEXH
            ZCONC   = XCCH*(ZLBDAH)**XCXH
            ZRADIUS = 0.5*(3.0/ZLBDAH) ! Assume Marshall-Palmer law for Reff
            IANGLE  = 11
            CALL BHMIE_WATER( ZWAVE_LENGTH, ZZREFIND_ICE, YDSD, ZCONC,      &
                              IANGLE, ZEXT_COEF, ZBAK_COEF, PRADIUS=ZRADIUS )
            ZALPH_PAR(JI,JJ,JK) = ZALPH_PAR(JI,JJ,JK) + ZETACLD * ZEXT_COEF
            ZBETA_PAR(JI,JJ,JK) = ZBETA_PAR(JI,JJ,JK) + ZBAK_COEF
          END IF
        END DO
      END DO
    END DO
END SELECT
!
! -----------------------------------------------------------------------------
!
!*       3.    PERFORMS THE BOTTOM-UP OR TOP-DOWN VERTICAL INTEGRATION
!              -------------------------------------------------------
!
!
ALLOCATE(ZOPTD_TOT(SIZE(PRHO,1),SIZE(PRHO,2)))
ALLOCATE(ZOPTD_MOL(SIZE(PRHO,1),SIZE(PRHO,2)))
ALLOCATE(ZOPTD_PAR(SIZE(PRHO,1),SIZE(PRHO,2)))
ZOPTD_TOT(:,:) = 0.
ZOPTD_MOL(:,:) = 0.
ZOPTD_PAR(:,:) = 0.
!
IF( HVIEW=='ZENIT' ) THEN
  IALT=IKB
  IF (PALT/=0.) THEN
    IKMIN=MINLOC(ABS(PZZ(:,:,:)-PALT))
    IALT=MIN(MAX(IKB,IKMIN(3)),IKE)
  ENDIF
  DO JK=IALT,IKE
!
! molecular optical depth
!
    ZOPTD_MOL(:,:) = ZOPTD_MOL(:,:) &
                   + ZALPH_MOL(:,:,JK)*(PZZ(:,:,JK)-PZZ(:,:,JK-1))
!
! Particular optical depth
!
    ZOPTD_PAR(:,:) = ZOPTD_PAR(:,:) &
                   + ZALPH_PAR(:,:,JK)*(PZZ(:,:,JK)-PZZ(:,:,JK-1))
!
! Total optical depth
!
    ZOPTD_TOT(:,:) = ZOPTD_MOL(:,:) + ZOPTD_PAR(:,:)
!
! Normalized Lidar profile
!
    PLIDAROUT(:,:,JK) = ( ZBETA_MOL(:,:,JK)+ZBETA_PAR(:,:,JK) ) &
                        * EXP( -2.0*ZOPTD_TOT(:,:) )
!
! Normalized Lidar particle profile
!
    PLIPAROUT(:,:,JK) = ZBETA_PAR(:,:,JK) * EXP( -2.0*ZOPTD_PAR(:,:) )
  END DO
ELSE IF( HVIEW=='NADIR' ) THEN
  IALT=IKE
  IF (PALT/=0.) THEN
    IKMIN=MINLOC(ABS(PZZ(:,:,:)-PALT))
    IALT=MIN(MAX(IKB,IKMIN(3)),IKE)
  ENDIF
  DO JK=IALT,IKB,-1
!
! molecular optical depth
!
    ZOPTD_MOL(:,:) = ZOPTD_MOL(:,:) &
                   + ZALPH_MOL(:,:,JK)*(PZZ(:,:,JK)-PZZ(:,:,JK-1))
!
! Particular optical depth
!
    ZOPTD_PAR(:,:) = ZOPTD_PAR(:,:) &
                   + ZALPH_PAR(:,:,JK)*(PZZ(:,:,JK)-PZZ(:,:,JK-1))
!
! Total optical depth
!
    ZOPTD_TOT(:,:) = ZOPTD_MOL(:,:) + ZOPTD_PAR(:,:)
!
! Normalized Lidar profile
!
    PLIDAROUT(:,:,JK) = ( ZBETA_MOL(:,:,JK)+ZBETA_PAR(:,:,JK) ) &
                        * EXP( -2.0*ZOPTD_TOT(:,:) )
!
! Normalized Lidar particle profile
!
    PLIPAROUT(:,:,JK) = ZBETA_PAR(:,:,JK) * EXP( -2.0*ZOPTD_PAR(:,:) )
  END DO
ENDIF
!
DEALLOCATE(ZOPTD_TOT,ZOPTD_MOL,ZOPTD_PAR)
!
!------------------------------------------------------------------------------
!
END SUBROUTINE LIDAR
