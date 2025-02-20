!MNH_LIC Copyright 2013-2024 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!      ###############################
       MODULE MODE_INI_LIMA_COLD_MIXED
!      ###############################
!
IMPLICIT NONE
!
CONTAINS
!
!     ###############################################
      SUBROUTINE INI_LIMA_COLD_MIXED (PTSTEP, PDZMIN)
!     ###############################################
!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to initialize the constants used in the 
!!    microphysical scheme LIMA for the cold and mixed phase variables
!!    and processes. 
!!
!!    AUTHOR
!!    ------
!!      J.-M. Cohard     * Laboratoire d'Aerologie*
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!      S.    Berthet    * Laboratoire d'Aerologie*
!!      B.    Vié        * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original             ??/??/13 
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!  P. Wautelet 10/04/2019: replace ABORT and STOP calls by Print_msg
!  P. Wautelet 26/04/2019: replace non-standard FLOAT function by REAL function
!  C. Barthe   14/03/2022: add CIBU and RDSF
!  J. Wurtz       03/2022: new snow characteristics
!  M. Taufour     07/2022: add concentration for snow, graupel, hail
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
!USE MODD_LUNIT, ONLY: TLUOUT0
USE MODD_PARAMETERS
USE MODD_PARAM_LIMA
USE MODD_PARAM_LIMA_WARM
USE MODD_PARAM_LIMA_COLD
USE MODD_PARAM_LIMA_MIXED
!
USE MODE_MSG
USE MODD_PRECISION, ONLY: MNHREAL64
!
USE MODE_LIMA_FUNCTIONS, ONLY: MOMG, GAUHER
USE MODI_GAMMA
USE MODI_GAMMA_INC
USE MODE_LIMA_RRCOLSS, ONLY: LIMA_RRCOLSS
USE MODE_LIMA_RZCOLX, ONLY: LIMA_RZCOLX
USE MODE_LIMA_RSCOLRG, ONLY: LIMA_RSCOLRG
USE MODE_NRCOLSS, ONLY: NRCOLSS
USE MODE_NZCOLX, ONLY: NZCOLX 
USE MODE_NSCOLRG, ONLY: NSCOLRG
USE MODE_LIMA_READ_XKER_RACCS, ONLY: LIMA_READ_XKER_RACCS
USE MODE_LIMA_READ_XKER_SDRYG, ONLY: LIMA_READ_XKER_SDRYG
USE MODE_LIMA_READ_XKER_RDRYG, ONLY: LIMA_READ_XKER_RDRYG
USE MODE_LIMA_READ_XKER_SWETH, ONLY: LIMA_READ_XKER_SWETH
USE MODE_LIMA_READ_XKER_GWETH, ONLY: LIMA_READ_XKER_GWETH
USE YOMHOOK, ONLY:LHOOK, DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
REAL,                    INTENT(IN) :: PTSTEP    ! Effective Time step 
REAL,                    INTENT(IN) :: PDZMIN    ! minimun vertical mesh size
!
!*       0.2   Declarations of local variables :
!
CHARACTER(LEN=13) :: YVAL     ! String for error message
INTEGER :: IKB                ! Coordinates of the first  physical 
                              ! points along z
INTEGER :: I1, ISH, ISH2      ! Internal loop indexes
!
REAL, DIMENSION(8)  :: ZGAMI  ! parameters involving various moments
REAL, DIMENSION(2)  :: ZGAMS  ! of the generalized gamma law
REAL, DIMENSION(2,4) :: ZGAMI_SHAPE  ! replace ZGAMI(6) and ZGAMI(7)   !++cb--
!
REAL :: ZRHO00                ! Surface reference air density
REAL :: ZRATE                 ! Geometrical growth of Lbda in the tabulated
                              ! functions and kernels
REAL :: ZBOUND                ! XDCSLIM*Lbda_s: upper bound for the partial
                              ! integration of the riming rate of the aggregates
REAL :: ZZEGS, ZZEGR, ZZEHS, ZZEHG! Bulk collection efficiencies
!
INTEGER :: IND                ! Number of interval to integrate the kernels
REAL :: ZZESR, ZZESS, ZZEII   ! Mean efficiency of rain-aggregate collection, aggregate-aggregate collection
REAL :: ZZFDINFTY             ! Factor used to define the "infinite" diameter
!
!
!INTEGER  :: ILUOUT0 ! Logical unit number for output-listing
!LOGICAL  :: GFLAG   ! Logical flag for printing the constatnts on the output
!                    ! listing
REAL     :: ZCONC_MAX ! Maximal concentration for snow
REAL     :: ZFACT_NUCL! Amplification factor for the minimal ice concentration
!  
INTEGER  :: IKND
INTEGER  :: IACCLBDAS,IACCLBDAR,IDRYLBDAG,IDRYLBDAS,IDRYLBDAR
REAL     :: ZALPHAR,ZALPHAS,ZALPHAG,ZALPHAH
REAL     :: ZNUR,ZNUS,ZNUG,ZNUH
REAL     :: ZBR,ZBS,ZBG
REAL     :: ZCR,ZCS,ZCG,ZCH
REAL     :: ZDR,ZDS,ZFVELOS,ZDG,ZDH
REAL     :: ZESR,ZEGS,ZEGR,ZEHS,ZEHG
REAL     :: ZFDINFTY
REAL     :: ZACCLBDAS_MAX,ZACCLBDAR_MAX,ZACCLBDAS_MIN,ZACCLBDAR_MIN
REAL     :: ZDRYLBDAG_MAX,ZDRYLBDAS_MAX,ZDRYLBDAG_MIN,ZDRYLBDAS_MIN
REAL     :: ZDRYLBDAR_MAX,ZDRYLBDAR_MIN
REAL     :: ZWETLBDAS_MAX,ZWETLBDAG_MAX,ZWETLBDAS_MIN,ZWETLBDAG_MIN
REAL     :: ZWETLBDAH_MAX,ZWETLBDAH_MIN
INTEGER  :: IWETLBDAS,IWETLBDAG,IWETLBDAH
!
REAL     :: ZFAC_ZRNIC ! Zrnic factor used to decrease Long Kernels
!
REAL :: ZBOUND_CIBU_SMIN    ! XDCSLIM*Lbda_s : lower bound used in the tabulated function
REAL :: ZBOUND_CIBU_GMIN    ! XDCGLIM*Lbda_g : lower bound used in the tabulated function
REAL :: ZRATE_S             ! Geometrical growth of Lbda_s in the tabulated function
REAL :: ZRATE_G             ! Geometrical growth of Lbda_g in the tabulated function
!
REAL :: ZBOUND_RDSF_RMIN    ! XDCRLIM*Lbda_r : lower bound used in the tabulated function
REAL :: ZRATE_R             ! Geometrical growth of Lbda_r in the tabulated function
REAL :: ZKHI_LWM            ! Coefficient of Lawson et al. (2015)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
REAL :: ZRHOIW ! ice density
!
!-------------------------------------------------------------------------------
!
!
!ILUOUT0 = TLUOUT0%NLU
IF (LHOOK) CALL DR_HOOK('INI_LIMA_COLD_MIXED', 0, ZHOOK_HANDLE)
CALL PARAM_LIMA_COLD_ASSOCIATE()
CALL PARAM_LIMA_MIXED_ASSOCIATE()
!
!
!*       1.     CHARACTERISTICS OF THE SPECIES
!               ------------------------------
!
!CALL RAIN_ICE_PARAM_ASSOCIATE()
!
!*       1.2    Ice crystal characteristics
!
SELECT CASE (CPRISTINE_ICE_LIMA)
  CASE('PLAT')
    XAI = 0.82      ! Plates
    XBI = 2.5       ! Plates 
    XC_I = 747.     ! Plates
    XDI = 1.0       ! Plates
    XGAMMAI = 0.096
    XDELTAI = 1.83    
    XC1I = 1./XPI   ! Plates
  CASE('COLU')
    XAI = 2.14E-3   ! Columns
    XBI = 1.7       ! Columns
    XC_I = 1.96E5   ! Columns
    XDI = 1.585     ! Columns
    XGAMMAI = 0.659
    XDELTAI = 2.0
    XC1I = 0.8      ! Columns 
  CASE('BURO')
    XAI = 44.0      ! Bullet rosettes
    XBI = 3.0       ! Bullet rosettes
    XC_I = 4.E5     ! Bullet rosettes
    XDI = 1.663     ! Bullet rosettes
    XGAMMAI = 0.062
    XDELTAI = 1.81    
    XC1I = 0.5      ! Bullet rosettes
  CASE('POIR')      ! polycristaux/irreguliers
    XAI = 10.0
    XBI = 2.7 
    XC_I = 7000.0
    XDI = 1.235
    XC1I = 0.566
  CASE('YPLA')
    XAI = 0.744      ! Plates_from Yang et al (2013)
    XBI = 2.47       ! Plates_from Yang et al (2013)
    XC_I = 48.9     ! Plates_from Yang et al (2013)
    XDI = 0.66       ! Plates_from Yang et al (2013)
    XGAMMAI = 0.106
    XDELTAI = 1.84    
    XC1I = 1./XPI   ! Plates_from Yang et al (2013)
  CASE('YCOL')
    XAI = 6.72   ! Columns_from Yang et al (2013)
    XBI = 2.6       ! Columns_from Yang et al (2013)
    XC_I = 99.3   ! Columns_from Yang et al (2013)
    XDI = 0.62     ! Columns_from Yang et al (2013)
    XGAMMAI = 0.154
    XDELTAI = 1.84    
    XC1I = 0.8      ! Columns_from Yang et al (2013)
  CASE('YBUR')
    XAI = 1.059      ! Bullet rosettes_from Yang et al (2013)
    XBI = 2.57       ! Bullet rosettes_from Yang et al (2013)
    XC_I = 99.1     ! Bullet rosettes_from Yang et al (2013)
    XDI = 0.7     ! Bullet rosettes_from Yang et al (2013)
    XGAMMAI = 0.057
    XDELTAI = 1.8    
    XC1I = 0.5      ! Bullet rosettes_from Yang et al (2013)
  CASE('YDRO')
    XAI = 346.906      ! Droxtals_from Yang et al (2013)
    XBI = 3       ! Droxtals_from Yang et al (2013)
    XC_I = 695.83     ! Droxtals_from Yang et al (2013)
    XDI = 0.61     ! Droxtals_from Yang et al (2013)
    XGAMMAI = 0.673
    XDELTAI = 2.0    
    XC1I = 0.5      ! Droxtals_from Yang et al (2013)
  CASE('YHCO')
    XAI = 5.6   ! Hollow_Columns_from Yang et al (2013)
    XBI = 2.6       ! Hollow_Columns_from Yang et al (2013)
    XC_I = 97.7   ! Hollow_Columns_from Yang et al (2013)
    XDI = 0.63     ! Hollow_Columns_from Yang et al (2013)
    XGAMMAI = 0.154
    XDELTAI = 1.84    
    XC1I = 0.8      ! Hollow_Columns_from Yang et al (2013)
  CASE('YHBU')
    XAI = 1.033      ! Hollow_Bullet rosettes_from Yang et al (2013)
    XBI = 2.59       ! Hollow_Bullet rosettes_from Yang et al (2013)
    XC_I = 112.5     ! Hollow_Bullet rosettes_from Yang et al (2013)
    XDI = 0.72     ! Hollow_Bullet rosettes_from Yang et al (2013)
    XGAMMAI = 0.057
    XDELTAI = 1.8    
    XC1I = 0.5      ! Hollow_Bullet rosettes_from Yang et al (2013)    
END SELECT
!
IF (LCRYSTAL_SHAPE) THEN
  CALL PARAM_LIMA_COLD_ALLOCATE_1D('XAI_SHAPE',  NNB_CRYSTAL_SHAPE)
  CALL PARAM_LIMA_COLD_ALLOCATE_1D('XBI_SHAPE',  NNB_CRYSTAL_SHAPE)
  CALL PARAM_LIMA_COLD_ALLOCATE_1D('XC_I_SHAPE', NNB_CRYSTAL_SHAPE)
  CALL PARAM_LIMA_COLD_ALLOCATE_1D('XDI_SHAPE',  NNB_CRYSTAL_SHAPE)
  CALL PARAM_LIMA_COLD_ALLOCATE_1D('XC1I_SHAPE', NNB_CRYSTAL_SHAPE)
  CALL PARAM_LIMA_COLD_ALLOCATE_1D('XGAMMAI_SHAPE', NNB_CRYSTAL_SHAPE)
  CALL PARAM_LIMA_COLD_ALLOCATE_1D('XDELTAI_SHAPE', NNB_CRYSTAL_SHAPE)
  !
  ZRHOIW = 916.0 ! kg/m-3 : ice density
! In order : plates, columns, irregulars, dendrites
! a, b, c et d : Woods et al 2007
! Irrreguliers : bullets de Woods et al 2007 (compris entre plates et columns)
! capacitance : choisir un rapport d'aspect pour plates et columns issu de Westbrook et al
!  Note that XCCI=N_i (a locally predicted value) and XCXI=0.0, implicitly
!

  DO ISH = 1, NNB_CRYSTAL_SHAPE
  SELECT CASE (HTYPE_CRYSTAL_SHAPE(ISH))
    CASE('CPLA')
      XAI_SHAPE(ISH)  = 0.82    ! relation masse maxdim
      XBI_SHAPE(ISH)  = 2.5     ! relation masse maxdim
      XC_I_SHAPE(ISH) = 747.0   ! vitesse de chute maxdim
      XDI_SHAPE(ISH)  = 1.0     ! vitesse de chute maxdim
      XC1I_SHAPE(ISH) = 0.339
      XGAMMAI_SHAPE(ISH) = 0.106 ! value from Yang et al. (2013)
      XDELTAI_SHAPE(ISH) = 1.84  ! value from Yang et al. (2013)
    CASE('CCOL')
      XAI_SHAPE(ISH)  = 110.798 ! Mitchell 1(d < 100 um) 1996
      XBI_SHAPE(ISH)  = 2.91    ! Mitchell 1(d < 100 um) 1996
      XC_I_SHAPE(ISH) = 1.96E5
      XDI_SHAPE(ISH)  = 1.585
      XC1I_SHAPE(ISH) = 0.24
      XGAMMAI_SHAPE(ISH) = 0.154 ! value from Yang et al. (2013)
      XDELTAI_SHAPE(ISH) = 1.84  ! value from Yang et al. (2013)
    CASE('POIR')
      XAI_SHAPE(ISH)  = 10.0    ! moyenne entre plates et colonnes (jpp)
      XBI_SHAPE(ISH)  = 2.7     ! relation masse diamètre des columns
      XC_I_SHAPE(ISH) = 7000.0
      XDI_SHAPE(ISH)  = 1.235
      XC1I_SHAPE(ISH) = 0.566
      XGAMMAI_SHAPE(ISH) = 0.057 ! value from Yang et al. (2013)
      XDELTAI_SHAPE(ISH) = 1.8  ! value from Yang et al. (2013)
    CASE('YPLA')
      XAI_SHAPE(ISH)  = 0.744      ! Plates_from Yang et al (2013)
      XBI_SHAPE(ISH)  = 2.47 
      XC_I_SHAPE(ISH) = 48.9    
      XDI_SHAPE(ISH)  = 0.66     
      XGAMMAI_SHAPE(ISH) = 0.106
      XDELTAI_SHAPE(ISH) = 1.84    
      XC1I_SHAPE(ISH) = 1./XPI   ! values from LIMA
    CASE('YCOL')
      XAI_SHAPE(ISH)  = 6.72   ! Columns_from Yang et al (2013)
      XBI_SHAPE(ISH)  = 2.6  
      XC_I_SHAPE(ISH) = 99.3
      XDI_SHAPE(ISH)  = 0.62 
      XGAMMAI_SHAPE(ISH) = 0.154
      XDELTAI_SHAPE(ISH) = 1.84    
      XC1I_SHAPE(ISH) = 0.8      ! values from LIMA
    CASE('YBUR')
      XAI_SHAPE(ISH)  = 1.059      ! Bullet rosettes_from Yang et al (2013)
      XBI_SHAPE(ISH)  = 2.57 
      XC_I_SHAPE(ISH) = 99.1
      XDI_SHAPE(ISH)  = 0.7  
      XGAMMAI_SHAPE(ISH) = 0.057
      XDELTAI_SHAPE(ISH) = 1.8    
      XC1I_SHAPE(ISH) = 0.5      !  values from LIMA      
    CASE('YDRO')
      XAI_SHAPE(ISH)  = 346.906      ! Droxtals_from Yang et al (2013)
      XBI_SHAPE(ISH)  = 3      
      XC_I_SHAPE(ISH) = 695.83 
      XDI_SHAPE(ISH)  = 0.61    
      XGAMMAI_SHAPE(ISH) = 0.673
      XDELTAI_SHAPE(ISH) = 2.0    
      XC1I_SHAPE(ISH) = 0.5      !  values from LIMA
    CASE('YHCO')
      XAI_SHAPE(ISH)  = 5.6   ! Hollow_Columns_from Yang et al (2013)
      XBI_SHAPE(ISH)  = 2.6   
      XC_I_SHAPE(ISH) = 97.7 
      XDI_SHAPE(ISH)  = 0.63  
      XGAMMAI_SHAPE(ISH) = 0.154
      XDELTAI_SHAPE(ISH) = 1.84    
      XC1I_SHAPE(ISH) = 0.8      ! values from LIMA
    CASE('YHBU')
      XAI_SHAPE(ISH)  = 1.033      ! Hollow_Bullet rosettes_from Yang et al (2013)
      XBI_SHAPE(ISH)  = 2.59  
      XC_I_SHAPE(ISH) = 112.5  
      XDI_SHAPE(ISH)  = 0.72   
      XGAMMAI_SHAPE(ISH) = 0.057
      XDELTAI_SHAPE(ISH) = 1.8    
      XC1I_SHAPE(ISH) = 0.5      ! values from LIMA            
  END SELECT
  ENDDO
!
! capacitances calculees a partir de Westbrook et al. (2007) pour des
! rapports d'aspect fixes : 
! A(plate) = 0.1 ; A(columns) = 5. ; A(irr) = 1. ; A'(dendrites) = 0.2
!
END IF
!
!  Note that XCCI=N_i (a locally predicted value) and XCXI=0.0, implicitly
!
XF0I = 1.00
! Correction BVIE XF2I from Pruppacher 1997 eq 13-88
!XF2I = 0.103
XF2I = 0.14
XF0IS = 0.86
XF1IS = 0.28
!
!*       1.3    Snowflakes/aggregates characteristics
!
XAS = 0.02
XBS = 1.9

IF (LSNOW_T) THEN
!Cas Gamma generalisee
XCS = 11.52
XDS = 0.39
XFVELOS =0.097
!Cas MP
!XCS = 13.2
!XDS = 0.423       
!XFVELOS = 25.14
ELSE
XCS = 5.
XDS = 0.27
XFVELOS = 0.
END IF

XCCS = 5.0
XCXS = 1.0

XF0S = 0.86
XF1S = 0.28
!
XC1S = 1./XPI
!
!*       1.4    Graupel characteristics
!
XAG = 19.6  ! Lump graupel case
XBG = 2.8   ! Lump graupel case
XCG = 122.  ! Lump graupel case
XDG = 0.66  ! Lump graupel case
!
XCCG = 5.E5
XCXG = -0.5
! XCCG = 4.E4 ! Test of Ziegler (1988)
! XCXG = -1.0 ! Test of Ziegler (1988)
!
XF0G = 0.86
XF1G = 0.28
!
XC1G = 1./2.
!
!*       2.5    Hailstone characteristics
!
!
XAH = 470.
XBH = 3.0
XCH = 201.
XDH = 0.64
!
!XCCH = 5.E-4
!XCXH = 2.0
!!!!!!!!!!!!
XCCH = 4.E4 ! Test of Ziegler (1988)
XCXH = -1.0 ! Test of Ziegler (1988)
!!!    XCCH = 5.E5 ! Graupel_like
!!!    XCXH = -0.5 ! Graupel_like
!!!!!!!!!!!!
!
XF0H = 0.86
XF1H = 0.28
!
XC1H = 1./2.
!  
!-------------------------------------------------------------------------------
!
!
!*       2.     DIMENSIONAL DISTRIBUTIONS OF THE SPECIES
!               ----------------------------------------
!
!
!*       2.1    Ice, snow, graupel and hail distribution
!
!
XALPHAI = 3.0  ! Gamma law for the ice crystal volume
XNUI    = 3.0  ! Gamma law with little dispersion
!
IF (LSNOW_T) THEN
!Cas GAMMAGEN
!   XALPHAS = .214   ! Generalized gamma law
!   XNUS    = 43.7   ! Generalized gamma law
   XALPHAS = 1.   ! Generalized gamma law
   XNUS    = 1.3   ! Generalized gamma law
   XTRANS_MP_GAMMAS = SQRT( ( GAMMA(DBLE(XNUS + 2._MNHREAL64/XALPHAS))*GAMMA(DBLE(XNUS + 4._MNHREAL64/XALPHAS)) ) / &
                            ( 8._MNHREAL64* GAMMA(DBLE(XNUS + 1._MNHREAL64/XALPHAS))*GAMMA(DBLE(XNUS + 3._MNHREAL64/XALPHAS)) ) )
ELSE IF (NMOM_S.EQ.2) THEN
   XALPHAS = 1.0  ! Gamma law 
   XNUS    = 2.0  !
   XTRANS_MP_GAMMAS = SQRT( ( GAMMA(DBLE(XNUS + 2._MNHREAL64/XALPHAS))*GAMMA(DBLE(XNUS + 4._MNHREAL64/XALPHAS)) ) / &
                            ( 8._MNHREAL64* GAMMA(DBLE(XNUS + 1._MNHREAL64/XALPHAS))*GAMMA(DBLE(XNUS + 3._MNHREAL64/XALPHAS)) ) )
ELSE
   XALPHAS = 1.0  ! Exponential law
   XNUS    = 1.0  ! Exponential law
   XTRANS_MP_GAMMAS = 1.
END IF
!
IF (NMOM_G.GE.2) then
   XALPHAG = 1.0  ! 
   XNUG    = 2.0  !
ELSE
   XALPHAG = 1.0  ! Exponential law
   XNUG    = 1.0  ! Exponential law
END IF
!
IF (NMOM_H.GE.2) then
   XALPHAH = 1.0  ! Gamma law
   XNUH    = 5.0  ! Gamma law with little dispersion
ELSE
   XALPHAH = 1.0  ! Gamma law
   XNUH    = 8.0  ! Gamma law with little dispersion
END IF
!
!*       2.2    Constants for shape parameter
!
XLBEXI = 1.0/XBI
XLBI   = XAI*MOMG(XALPHAI,XNUI,XBI)
!
IF (LCRYSTAL_SHAPE) THEN
  CALL PARAM_LIMA_COLD_ALLOCATE_1D('XLBEXI_SHAPE', NNB_CRYSTAL_SHAPE)
  CALL PARAM_LIMA_COLD_ALLOCATE_1D('XLBI_SHAPE', NNB_CRYSTAL_SHAPE)
  DO ISH = 1, NNB_CRYSTAL_SHAPE
    XLBEXI_SHAPE(ISH) = 1.0 / XBI_SHAPE(ISH)
    XLBI_SHAPE(ISH)   = XAI_SHAPE(ISH) * MOMG(XALPHAI,XNUI,XBI_SHAPE(ISH))
  END DO
END IF
!
XNS = 1.0/(XAS*MOMG(XALPHAS,XNUS,XBS))
IF (NMOM_S.EQ.1) THEN
   XLBEXS = 1.0/(XCXS-XBS)
   XLBS   = ( XAS*XCCS*MOMG(XALPHAS,XNUS,XBS) )**(-XLBEXS)
ELSE
   XLBEXS = 1./XBS
   XLBS   = XAS * MOMG(XALPHAS,XNUS,XBS)
END IF
!
XNG = 1.0/(XAG*MOMG(XALPHAG,XNUG,XBG))
IF (NMOM_G.EQ.1) THEN
   XLBEXG = 1.0/(XCXG-XBG)
   XLBG   = ( XAG*XCCG*MOMG(XALPHAG,XNUG,XBG))**(-XLBEXG)
ELSE
   XLBEXG = 1./XBG
   XLBG   = XAG * MOMG(XALPHAG,XNUG,XBG)
END IF
!
IF (NMOM_H.EQ.1) THEN
   XLBEXH = 1.0/(XCXH-XBH)
   XLBH   = ( XAH*XCCH*MOMG(XALPHAH,XNUH,XBH) )**(-XLBEXH)
ELSE
   XLBEXH = 1./XBH
   XLBH   = XAH * MOMG(XALPHAH,XNUH,XBH)
END IF
!
!!$GFLAG = .TRUE.
!!$IF (GFLAG) THEN
!!$  WRITE(UNIT=ILUOUT0,FMT='("      Shape Parameters")')
!!$  WRITE(UNIT=ILUOUT0,FMT='(" XLBEXI =",E13.6," XLBI =",E13.6)') XLBEXI,XLBI
!!$  WRITE(UNIT=ILUOUT0,FMT='(" XLBEXS =",E13.6," XLBS =",E13.6)') XLBEXS,XLBS
!!$  WRITE(UNIT=ILUOUT0,FMT='(" XLBEXG =",E13.6," XLBG =",E13.6)') XLBEXG,XLBG
!!$  WRITE(UNIT=ILUOUT0,FMT='(" XLBEXH =",E13.6," XLBH =",E13.6)') XLBEXH,XLBH
!!$END IF
!
XLBDAS_MAX = 1.E7 ! (eq to r~1E-7kg/kg) (for non MP PSD, use conversion XTRANS_MP_GAMMAS)
XLBDAS_MIN = 1.   ! (eq to r~0.18kg/kg) (for non MP PSD, use conversion XTRANS_MP_GAMMAS)
XLBDAG_MAX = 100000.0
!
ZCONC_MAX  = 1.E6 ! Maximal concentration for falling particules set to 1 per cc
!XLBDAS_MAX = ( ZCONC_MAX/XCCS )**(1./XCXS) 
!XLBDAG_MAX = ( ZCONC_MAX/XCCG )**(1./XCXG) 
!XLBDAH_MAX = ( ZCONC_MAX/XCCH )**(1./XCXH) 
!  
!    constants for ecRad effective radius
ZRHOIW = 0.917
XREFFI = (3*XAI/(2*ZRHOIW*10.**3*XGAMMAI)*MOMG(XALPHAI,XNUI,XBI)/MOMG(XALPHAI,XNUI,XDELTAI))*1E6
!
!-------------------------------------------------------------------------------
!
!
!*       3.     CONSTANTS FOR THE SEDIMENTATION
!               -------------------------------
!
!
!*       3.1    Exponent of the fall-speed air density correction
!
IKB = 1 + JPVEXT
! Correction
! ZRHO00 = XP00/(XRD*XTHVREFZ(IKB))
ZRHO00 = 1.2041 ! at P=1013.25hPa and T=20°C
!
!*       3.2    Constants for sedimentation
!
!! XEXRSEDI = (XBI+XDI)/XBI
!! XEXCSEDI = 1.0-XEXRSEDI
!! XFSEDI   = (4.*XPI*900.)**(-XEXCSEDI) *                         &
!!            XC_I*XAI*MOMG(XALPHAI,XNUI,XBI+XDI) *                &
!!            ((XAI*MOMG(XALPHAI,XNUI,XBI)))**(-XEXRSEDI) *        &
!!            (ZRHO00)**XCEXVT
!! !
!! !  Computations made for Columns
!! !
!! XEXRSEDI = 1.9324
!! XEXCSEDI =-0.9324
!! XFSEDI   = 3.89745E11*MOMG(XALPHAI,XNUI,3.285)*                          &
!!                       MOMG(XALPHAI,XNUI,1.7)**(-XEXRSEDI)*(ZRHO00)**XCEXVT
!! XEXCSEDI =-0.9324*3.0
!! WRITE (ILUOUT0,FMT=*)' PRISTINE ICE SEDIMENTATION for columns XFSEDI=',XFSEDI
!
!
XFSEDRI  = XC_I*GAMMA_X0D(XNUI+(XDI+XBI)/XALPHAI)/GAMMA_X0D(XNUI+XBI/XALPHAI)*     &
            (ZRHO00)**XCEXVT
XFSEDCI  = XC_I*GAMMA_X0D(XNUI+XDI/XALPHAI)/GAMMA_X0D(XNUI)*     &
            (ZRHO00)**XCEXVT
!
IF (LCRYSTAL_SHAPE) THEN
  CALL PARAM_LIMA_COLD_ALLOCATE_1D('XFSEDCI_SHAPE', NNB_CRYSTAL_SHAPE)
  CALL PARAM_LIMA_COLD_ALLOCATE_1D('XFSEDRI_SHAPE', NNB_CRYSTAL_SHAPE)
  CALL PARAM_LIMA_COLD_ALLOCATE_1D('XFSEDRI_TOT_SHAPE', NNB_CRYSTAL_SHAPE)
!
  DO ISH = 1, NNB_CRYSTAL_SHAPE
    XFSEDCI_SHAPE(ISH) = XC_I_SHAPE(ISH) * GAMMA_X0D(XNUI+XDI_SHAPE(ISH)/XALPHAI) / &
                         GAMMA_X0D(XNUI) * (ZRHO00)**XCEXVT
    XFSEDRI_SHAPE(ISH) = XC_I_SHAPE(ISH) * GAMMA_X0D(XNUI+(XDI_SHAPE(ISH)+XBI_SHAPE(ISH))/XALPHAI) / &
                         GAMMA_X0D(XNUI+XBI_SHAPE(ISH)/XALPHAI) * (ZRHO00)**XCEXVT
! marine : correction de l'air ?
    XFSEDRI_TOT_SHAPE(ISH) = XAI_SHAPE(ISH) * GAMMA_X0D(XNUI+XBI_SHAPE(ISH)/XALPHAI) / &
                                              GAMMA_X0D(XNUI)
!    XFSEDRI_TOT_SHAPE(ISH) = XAI_SHAPE(ISH) * GAMMA_X0D(XNUI+XBI_SHAPE(ISH)/XALPHAI) / &
!                                              GAMMA_X0D(XNUI) * (ZRHO00)**XCEXVT
  END DO
END IF
!
IF (LSNOW_T) THEN
!HOUZE/HAIC
   !XEXSEDS = -XDS !(2*XBS+XDS)
   !XFSEDS  = XCS*MOMG(XALPHAS,XNUS,XBS+XDS)/(MOMG(XALPHAS,XNUS,XBS))    &
   !            *(ZRHO00)**XCEXVT
!LH_EXTENDED
   XEXSEDS = -XDS-XBS
   XFSEDS  = XCS*MOMG(XALPHAS,XNUS,XBS+XDS)/(MOMG(XALPHAS,XNUS,XBS))    &
            *(ZRHO00)**XCEXVT
ELSE
   XEXSEDS = (XBS+XDS-XCXS)/(XBS-XCXS)
   XFSEDS  = XCS*XAS*XCCS*MOMG(XALPHAS,XNUS,XBS+XDS)*                         &
            (XAS*XCCS*MOMG(XALPHAS,XNUS,XBS))**(-XEXSEDS)*(ZRHO00)**XCEXVT
END IF
!
XEXSEDG = (XBG+XDG-XCXG)/(XBG-XCXG)
XFSEDG  = XCG*XAG*XCCG*MOMG(XALPHAG,XNUG,XBG+XDG)*                         &
            (XAG*XCCG*MOMG(XALPHAG,XNUG,XBG))**(-XEXSEDG)*(ZRHO00)**XCEXVT
!
XEXSEDH = (XBH+XDH-XCXH)/(XBH-XCXH)
XFSEDH  = XCH*XAH*XCCH*MOMG(XALPHAH,XNUH,XBH+XDH)*                         &
            (XAH*XCCH*MOMG(XALPHAH,XNUH,XBH))**(-XEXSEDH)*(ZRHO00)**XCEXVT
IF (NMOM_S.GE.2) THEN
  XFSEDRS  = XCS*GAMMA_X0D(XNUS+(XDS+XBS)/XALPHAS)/GAMMA_X0D(XNUS+XBS/XALPHAS)*  &   
            (ZRHO00)**XCEXVT                                                       
  XFSEDCS  = XCS*GAMMA_X0D(XNUS+XDS/XALPHAS)/GAMMA_X0D(XNUS)*                     & 
            (ZRHO00)**XCEXVT
END IF
IF (NMOM_G.GE.2) THEN
  XFSEDRG  = XCG*GAMMA_X0D(XNUG+(XDG+XBG)/XALPHAG)/GAMMA_X0D(XNUG+XBG/XALPHAG)*  &  
               (ZRHO00)**XCEXVT                                                      
  XFSEDCG  = XCG*GAMMA_X0D(XNUG+XDG/XALPHAG)/GAMMA_X0D(XNUG)*                     & 
               (ZRHO00)**XCEXVT 
END IF
IF (NMOM_H.GE.2) THEN
  XFSEDRH  = XCH*GAMMA_X0D(XNUH+(XDH+XBH)/XALPHAH)/GAMMA_X0D(XNUH+XBH/XALPHAH)*  &  
            (ZRHO00)**XCEXVT                                                        
  XFSEDCH  = XCH*GAMMA_X0D(XNUH+XDH/XALPHAH)/GAMMA_X0D(XNUH)*                     & 
            (ZRHO00)**XCEXVT
END IF
!
!
!
XLB(4)    = XLBI
XLBEX(4)  = XLBEXI
XD(4)     = XDI
XFSEDR(4) = XFSEDRI
XFSEDC(4) = XFSEDCI
!
XLB(5)    = XLBS
XLBEX(5)  = XLBEXS
XD(5)     = XDS
XFSEDR(5) = XCS*GAMMA_X0D(XNUS+(XDS+XBS)/XALPHAS)/GAMMA_X0D(XNUS+XBS/XALPHAS)*     &
            (ZRHO00)**XCEXVT
XFSEDC(5) = XCS*GAMMA_X0D(XNUS+XDS/XALPHAS)/GAMMA_X0D(XNUS)*                     & 
            (ZRHO00)**XCEXVT
!
XLB(6)    = XLBG
XLBEX(6)  = XLBEXG
XD(6)     = XDG
XFSEDR(6) = XCG*GAMMA_X0D(XNUG+(XDG+XBG)/XALPHAG)/GAMMA_X0D(XNUG+XBG/XALPHAG)*     &
            (ZRHO00)**XCEXVT
XFSEDC(6) = XCG*GAMMA_X0D(XNUG+XDG/XALPHAG)/GAMMA_X0D(XNUG)*                     & 
            (ZRHO00)**XCEXVT
!
XLB(7)    = XLBH
XLBEX(7)  = XLBEXH
XD(7)     = XDH
XFSEDR(7) = XCH*GAMMA_X0D(XNUH+(XDH+XBH)/XALPHAH)/GAMMA_X0D(XNUH+XBH/XALPHAH)*     &
            (ZRHO00)**XCEXVT
XFSEDC(7) = XCH*GAMMA_X0D(XNUH+XDH/XALPHAH)/GAMMA_X0D(XNUH)*                     & 
            (ZRHO00)**XCEXVT
!
!-------------------------------------------------------------------------------
!
!
!*       4.     CONSTANTS FOR HETEROGENEOUS NUCLEATION
!               --------------------------------------
!
!
!                 ***************
!*          4.1   LIMA_NUCLEATION
!                 ***************   
!*               4.1.1  Constants for the computation of the number concentration
!                       of active IN 
!
XRHO_CFDC  = 0.76 
!
XGAMMA     = 2.
!
IF (NPHILLIPS == 13) THEN
   XAREA1(1)  = 2.0E-6    !DM1
   XAREA1(2)  = XAREA1(1) !DM2
   XAREA1(3)  = 1.0E-7    !BC
   XAREA1(4)  = 8.9E-7    !BIO
ELSE IF (NPHILLIPS == 8) THEN
   XAREA1(1)  = 2.0E-6    !DM1
   XAREA1(2)  = XAREA1(1) !DM2
   XAREA1(3)  = 2.7E-7    !BC
   XAREA1(4)  = 9.1E-7    !BIO
ELSE
  CALL PRINT_MSG( NVERB_FATAL, 'GEN', 'INI_LIMA_COLD_MIXED', 'NPHILLIPS SHOULD BE EQUAL TO 8 OR 13' )
END IF
!
!*               4.1.2  Constants for the computation of H_X (the fraction-redu-
!                       cing IN activity at low S_i and warm T) for X={DM1,DM2,BC,BIO} 
!                       
!
IF (NPHILLIPS == 13) THEN
   XDT0(1)   =   5.  !DM1
   XDT0(2)   =   5.  !DM2
   XDT0(3)   =  10.  !BC
   XDT0(4)   =   5.  !BIOO
!
   XT0(1) = -40. +273.15 !DM1
   XT0(2) = XT0(1)       !DM2
   XT0(3) = -50. +273.15 !BC
   XT0(4) = -20. +273.15 !BIO
!
   XSW0   = 0.97
!
   XDSI0(1)  = 0.1 !DM1
   XDSI0(2)  = 0.1 !DM2
   XDSI0(3)  = 0.1 !BC
   XDSI0(4)  = 0.2 !BIO
!
   XH(1) = 0.15          !DM1
   XH(2) = 0.15          !DM2
   XH(3) = 0.            !BC
   XH(4) = 0.            !O
!
   XTX1(1) = -30. +273.15 !DM1
   XTX1(2) = XTX1(1)      !DM2
   XTX1(3) = -25. +273.15 !BC
   XTX1(4) =  -5. +273.15 !BIO
!
   XTX2(1) = -10. +273.15 !DM1
   XTX2(2) = XTX2(1)      !DM2
   XTX2(3) = -15. +273.15 !BC
   XTX2(4) =  -2. +273.15 !BIO
ELSE IF (NPHILLIPS == 8) THEN
   XDT0(1)   =   5.  !DM1
   XDT0(2)   =   5.  !DM2
   XDT0(3)   =   5.  !BC
   XDT0(4)   =   5.  !O
!
   XT0(1) = -40. +273.15 !DM1
   XT0(2) = XT0(1)       !DM2
   XT0(3) = -50. +273.15 !BC
   XT0(4) = -50. +273.15 !BIO
!
   XSW0   = 0.97
!
   XDSI0(1)  = 0.1 !DM1
   XDSI0(2)  = 0.1 !DM2
   XDSI0(3)  = 0.1 !BC
   XDSI0(4)  = 0.1 !BIO
!
   XH(1) = 0.15          !DM1
   XH(2) = 0.15          !DM2
   XH(3) = 0.            !BC
   XH(4) = 0.            !O
!
   XTX1(1) =  -5. +273.15 !DM1
   XTX1(2) = XTX1(1)      !DM2
   XTX1(3) =  -5. +273.15 !BC
   XTX1(4) =  -5. +273.15 !BIO
!
   XTX2(1) =  -2. +273.15 !DM1
   XTX2(2) = XTX2(1)      !DM2
   XTX2(3) =  -2. +273.15 !BC
   XTX2(4) =  -2. +273.15 !BIO
END IF
!
!*               4.1.3  Constants for the computation of the Gauss Hermitte 
!                       quadrature method used for the integration of the total
!                       crystal number over T>-35°C
!
NDIAM = 70
!
CALL PARAM_LIMA_ALLOCATE('XABSCISS', NDIAM)
CALL PARAM_LIMA_ALLOCATE('XWEIGHT',  NDIAM)
!
CALL GAUHER(XABSCISS, XWEIGHT, NDIAM)
!
!                 ***************** 
!*          4.2   MEYERS NUCLEATION
!                 ***************** 
!
ZFACT_NUCL =  1.0    ! Plates, Columns and Bullet rosettes
!
!*               5.2.1  Constants for nucleation from ice nuclei
!
XNUC_DEP  = XFACTNUC_DEP*1000.*ZFACT_NUCL
XEXSI_DEP = 12.96E-2
XEX_DEP   = -0.639
!
XNUC_CON   = XFACTNUC_CON*1000.*ZFACT_NUCL
XEXTT_CON  = -0.262
XEX_CON    = -2.8
!
XMNU0 = 6.88E-13
!
!!$IF (LMEYERS) THEN
!!$  WRITE(UNIT=ILUOUT0,FMT='("      Heterogeneous nucleation")')
!!$  WRITE(UNIT=ILUOUT0,FMT='(" XNUC_DEP=",E13.6," XEXSI=",E13.6," XEX=",E13.6)') &
!!$                                                      XNUC_DEP,XEXSI_DEP,XEX_DEP
!!$  WRITE(UNIT=ILUOUT0,FMT='(" XNUC_CON=",E13.6," XEXTT=",E13.6," XEX=",E13.6)') &
!!$                                                      XNUC_CON,XEXTT_CON,XEX_CON
!!$  WRITE(UNIT=ILUOUT0,FMT='(" mass of embryo XMNU0=",E13.6)') XMNU0
!!$END IF
!
!                 ***************** 
!*          4.3   NUCLEATION for NMOM_I=1
!                 ***************** 
!
XNU10 = 50.*ZFACT_NUCL
XALPHA1 = 4.5
XBETA1 = 0.6
!
XNU20 = 1000.*ZFACT_NUCL
XALPHA2 = 12.96
XBETA2 = 0.639
!
!XMNU0 = 6.88E-13
!-------------------------------------------------------------------------------
!
!
!*       5.     CONSTANTS FOR THE SLOW COLD PROCESSES
!               -------------------------------------
!
!
!*       5.1.2  Constants for homogeneous nucleation from haze particules
!
XRHOI_HONH = 925.0
XCEXP_DIFVAP_HONH = 1.94
XCOEF_DIFVAP_HONH = (2.0*XPI)*0.211E-4*XP00/XTT**XCEXP_DIFVAP_HONH
XCRITSAT1_HONH = 2.583
XCRITSAT2_HONH = 207.83
XTMIN_HONH = 180.0
XTMAX_HONH = 240.0
XDLNJODT1_HONH = 4.37
XDLNJODT2_HONH = 0.03
XC1_HONH = 100.0
XC2_HONH = 22.6
XC3_HONH = 0.1
XRCOEF_HONH = (XPI/6.0)*XRHOI_HONH
!
!
!*       5.1.3  Constants for homogeneous nucleation from cloud droplets
!
XTEXP1_HONC = -606.3952*LOG(10.0)
XTEXP2_HONC =  -52.6611*LOG(10.0)
XTEXP3_HONC =   -1.7439*LOG(10.0)
XTEXP4_HONC =   -0.0265*LOG(10.0)
XTEXP5_HONC = -1.536E-4*LOG(10.0)
IF (XALPHAC == 3.0) THEN
  XC_HONC   = XPI/6.0
  XR_HONC   = XPI/6.0
ELSE
  WRITE ( YVAL, '( E13.6 )' ) XALPHAC
  CALL PRINT_MSG( NVERB_FATAL, 'GEN', 'INI_LIMA_COLD_MIXED', 'HOMOGENEOUS NUCLEATION: XALPHAC='//TRIM(YVAL)// &
                  '/= 3. NO ALGORITHM DEVELOPED FOR THIS CASE' )
END IF
!
!!$GFLAG = .TRUE.
!!$IF (GFLAG) THEN
!!$  WRITE(UNIT=ILUOUT0,FMT='("      Homogeneous nucleation")')
!!$  WRITE(UNIT=ILUOUT0,FMT='(" XTEXP1_HONC=",E13.6)') XTEXP1_HONC
!!$  WRITE(UNIT=ILUOUT0,FMT='(" XTEXP2_HONC=",E13.6)') XTEXP2_HONC
!!$  WRITE(UNIT=ILUOUT0,FMT='(" XTEXP3_HONC=",E13.6)') XTEXP3_HONC
!!$  WRITE(UNIT=ILUOUT0,FMT='(" XTEXP4_HONC=",E13.6)') XTEXP4_HONC
!!$  WRITE(UNIT=ILUOUT0,FMT='(" XTEXP5_HONC=",E13.6)') XTEXP5_HONC
!!$  WRITE(UNIT=ILUOUT0,FMT='("XC_HONC=",E13.6," XR_HONC=",E13.6)') XC_HONC,XR_HONC
!!$END IF
!
!
!*       5.2    Constants for vapor deposition on ice
!
XSCFAC = (0.63**(1./3.))*SQRT((ZRHO00)**XCEXVT) ! One assumes Sc=0.63
!
X0DEPI = (4.0*XPI)*XC1I*XF0I*MOMG(XALPHAI,XNUI,1.)
X2DEPI = (4.0*XPI)*XC1I*XF2I*XC_I*MOMG(XALPHAI,XNUI,XDI+2.0)
!
! Harrington parameterization for ice to snow conversion
IF (LCRYSTAL_SHAPE) THEN
  CALL PARAM_LIMA_COLD_ALLOCATE_1D('X0DEPI_SHAPE', NNB_CRYSTAL_SHAPE)
  CALL PARAM_LIMA_COLD_ALLOCATE_1D('X2DEPI_SHAPE', NNB_CRYSTAL_SHAPE)
  DO ISH = 1, NNB_CRYSTAL_SHAPE
    X0DEPI_SHAPE(ISH) = (4.0 * XPI) * XC1I_SHAPE(ISH) * XF0I * MOMG(XALPHAI,XNUI,1.)
    X2DEPI_SHAPE(ISH) = (4.0 * XPI) * XC1I_SHAPE(ISH) * XF2I * XC_I_SHAPE(ISH) * &
                        MOMG(XALPHAI,XNUI,XDI_SHAPE(ISH)+2.0)
  END DO
END IF
!
! Harrington parameterization for ice to snow conversion
!
XDICNVS_LIM = 125.E-6  ! size in microns
XLBDAICNVS_LIM = (50.0**(1.0/(XALPHAI)))/XDICNVS_LIM  ! ZLBDAI Limitation
XC0DEPIS = ((4.0*XPI)/(XAI*XBI))*XC1I*XF0IS*                        &
           (XALPHAI/GAMMA_X0D(XNUI))*XDICNVS_LIM**(1.0-XBI)
XC1DEPIS = ((4.0*XPI)/(XAI*XBI))*XC1I*XF1IS*SQRT(XC_I)*             &
           (XALPHAI/GAMMA_X0D(XNUI))*XDICNVS_LIM**(1.0-XBI+(XDI+1.0)/2.0)
XR0DEPIS = XC0DEPIS *(XAI*XDICNVS_LIM**XBI)
XR1DEPIS = XC1DEPIS *(XAI*XDICNVS_LIM**XBI)
!
IF (LCRYSTAL_SHAPE) THEN
  CALL PARAM_LIMA_COLD_ALLOCATE_1D('XC0DEPIS_SHAPE', NNB_CRYSTAL_SHAPE)
  CALL PARAM_LIMA_COLD_ALLOCATE_1D('XC1DEPIS_SHAPE', NNB_CRYSTAL_SHAPE)
  CALL PARAM_LIMA_COLD_ALLOCATE_1D('XR0DEPIS_SHAPE', NNB_CRYSTAL_SHAPE)
  CALL PARAM_LIMA_COLD_ALLOCATE_1D('XR1DEPIS_SHAPE', NNB_CRYSTAL_SHAPE)
  DO ISH = 1, NNB_CRYSTAL_SHAPE
    XC0DEPIS_SHAPE(ISH) = ((4.0 * XPI) / (XAI_SHAPE(ISH) * XBI_SHAPE(ISH))) *     &
                          XC1I_SHAPE(ISH) * XF0IS * (XALPHAI / GAMMA_X0D(XNUI)) * &
                          XDICNVS_LIM**(1.0-XBI_SHAPE(ISH))
    XC1DEPIS_SHAPE(ISH) = ((4.0 * XPI) / (XAI_SHAPE(ISH) * XBI_SHAPE(ISH))) *  &
                          XC1I_SHAPE(ISH) * XF1IS * SQRT(XC_I_SHAPE(ISH)) *    &
                          (XALPHAI / GAMMA_X0D(XNUI)) *                        &
                          XDICNVS_LIM**(1.0-XBI_SHAPE(ISH)+(XDI_SHAPE(ISH)+1.0)/2.0)
    XR0DEPIS_SHAPE(ISH) = XC0DEPIS_SHAPE(ISH) * &
                         (XAI_SHAPE(ISH) * XDICNVS_LIM**XBI_SHAPE(ISH))
    XR1DEPIS_SHAPE(ISH) = XC1DEPIS_SHAPE(ISH) * &
                         (XAI_SHAPE(ISH) * XDICNVS_LIM**XBI_SHAPE(ISH))
  END DO
END IF
!
! Harrington parameterization for snow to ice conversion
!
XLBDASCNVI_MAX = 6000.*XTRANS_MP_GAMMAS ! lbdas max after Field (1999)
!
XDSCNVI_LIM    = 125.E-6  ! size in microns
XLBDASCNVI_LIM = (50.0**(1.0/(XALPHAS)))/XDSCNVI_LIM  ! ZLBDAS Limitation
XC0DEPSI = ((4.0*XPI)/(XAS*XBS))*XC1S*XF0IS*                        &
           (XALPHAS/GAMMA_X0D(XNUS))*XDSCNVI_LIM**(1.0-XBS)
XC1DEPSI = ((4.0*XPI)/(XAS*XBS))*XC1S*XF1IS*SQRT(XCS)*             &
           (XALPHAS/GAMMA_X0D(XNUS))*XDSCNVI_LIM**(1.0-XBS+(XDS+1.0)/2.0)
XR0DEPSI = XC0DEPSI *(XAS*XDSCNVI_LIM**XBS)
XR1DEPSI = XC1DEPSI *(XAS*XDSCNVI_LIM**XBS)
!
! Vapor deposition on snow and graupel and hail
!
X0DEPS = (4.0*XPI)*XC1S*XF0S*MOMG(XALPHAS,XNUS,1.)
X1DEPS = (4.0*XPI)*XC1S*XF1S*SQRT(XCS)*MOMG(XALPHAS,XNUS,0.5*XDS+1.5)
XEX0DEPS = -1.0
XEX1DEPS = -0.5*(XDS+3.0)
!
X0DEPG = (4.0*XPI)*XC1G*XF0G*MOMG(XALPHAG,XNUG,1.)  
X1DEPG = (4.0*XPI)*XC1G*XF1G*SQRT(XCG)*MOMG(XALPHAG,XNUG,0.5*XDG+1.5) 
XEX0DEPG = -1.0                                                       
XEX1DEPG = -0.5*(XDG+3.0)
!
X0DEPH = (4.0*XPI)*XC1H*XF0H*MOMG(XALPHAH,XNUH,1.)   
X1DEPH = (4.0*XPI)*XC1H*XF1H*SQRT(XCH)*MOMG(XALPHAH,XNUH,0.5*XDH+1.5) 
XEX0DEPH = -1.0                                                        
XEX1DEPH = -0.5*(XDH+3.0)
!
!-------------------------------------------------------------------------------
!
!
!*       6.     CONSTANTS FOR THE COALESCENCE PROCESSES
!               ---------------------------------------
!
!
!*       6.0    Precalculation of the gamma function momentum
!
ZGAMI(1) = GAMMA_X0D(XNUI)
ZGAMI(2) = MOMG(XALPHAI,XNUI,3.)
ZGAMI(3) = MOMG(XALPHAI,XNUI,6.)
ZGAMI(4) = ZGAMI(3)-ZGAMI(2)**2  ! useful for Sig_I
ZGAMI(5) = MOMG(XALPHAI,XNUI,9.)
ZGAMI(6) = MOMG(XALPHAI,XNUI,3.+XBI)
ZGAMI(7) = MOMG(XALPHAI,XNUI,XBI)
ZGAMI(8) = MOMG(XALPHAI,XNUI,3.)/MOMG(XALPHAI,XNUI,2.)
!
IF (LCRYSTAL_SHAPE) THEN
  DO ISH = 1, NNB_CRYSTAL_SHAPE
    ZGAMI_SHAPE(1,ISH) = MOMG(XALPHAI,XNUI,3.+XBI_SHAPE(ISH))
    ZGAMI_SHAPE(2,ISH) = MOMG(XALPHAI,XNUI,XBI_SHAPE(ISH))
  END DO
END IF
!
ZGAMS(1) = GAMMA_X0D(XNUS)
ZGAMS(2) = MOMG(XALPHAS,XNUS,3.)
!
!
!*       6.1    Csts for the coalescence processes
!
ZFAC_ZRNIC = 0.1
XKER_ZRNIC_A1 = 2.59E15*ZFAC_ZRNIC**2! From Long  a1=9.44E9 cm-3 
                                     ! so XKERA1= 9.44E9*1E6*(PI/6)**2
XKER_ZRNIC_A2 = 3.03E3*ZFAC_ZRNIC    ! From Long  a2=5.78E3      
                                     ! so XKERA2= 5.78E3*    (PI/6)
!
!
!*       6.2    Csts for the pristine ice selfcollection process
!
XSELFI = XKER_ZRNIC_A1*ZGAMI(3)
XCOLEXII = 0.025   !  Temperature factor of the I+I collection efficiency
!
!
!*       6.3    Constants for pristine ice autoconversion
!
XTEXAUTI = 0.025   !  Temperature factor of the I+I collection efficiency
!
!!$GFLAG = .TRUE.
!!$IF (GFLAG) THEN
!!$  WRITE(UNIT=ILUOUT0,FMT='("      pristine ice autoconversion")')
!!$  WRITE(UNIT=ILUOUT0,FMT='(" Temp. factor    XTEXAUTI=",E13.6)') XTEXAUTI
!!$END IF
!
XAUTO3 = 6.25E18*(ZGAMI(2))**(1./3.)*SQRT(ZGAMI(4))
XAUTO4 = 0.5E6*(ZGAMI(4))**(1./6.)
XLAUTS = 2.7E-2
XLAUTS_THRESHOLD  = 0.4
XITAUTS= 0.27 ! (Notice that T2 of BR74 is uncorrect and that 0.27=1./3.7
XITAUTS_THRESHOLD = 7.5
!
!
!*       6.4    Constants for snow aggregation
!
XCOLEXIS = 0.05    ! Temperature factor of the I+S collection efficiency
XFIAGGS  = (XPI/4.0)*0.25*XCS*(ZRHO00**XCEXVT)*MOMG(XALPHAS,XNUS,XDS+2.0)
XEXIAGGS = -XDS - 2.0
XAGGS_CLARGE1 = XKER_ZRNIC_A2*ZGAMI(2)
XAGGS_CLARGE2 = XKER_ZRNIC_A2*ZGAMS(2)
XAGGS_RLARGE1 = XKER_ZRNIC_A2*ZGAMI(6)*XAI
XAGGS_RLARGE2 = XKER_ZRNIC_A2*ZGAMI(7)*ZGAMS(2)*XAI
!
IF (LCRYSTAL_SHAPE) THEN
  CALL PARAM_LIMA_COLD_ALLOCATE_1D('XAGGS_RLARGE1_SHAPE', NNB_CRYSTAL_SHAPE)
  CALL PARAM_LIMA_COLD_ALLOCATE_1D('XAGGS_RLARGE2_SHAPE', NNB_CRYSTAL_SHAPE)
  DO ISH = 1, NNB_CRYSTAL_SHAPE
    XAGGS_RLARGE1_SHAPE(ISH) = XKER_ZRNIC_A2 * ZGAMI_SHAPE(1,ISH) * XAI_SHAPE(ISH)
    XAGGS_RLARGE2_SHAPE(ISH) = XKER_ZRNIC_A2 * ZGAMI_SHAPE(2,ISH) * ZGAMS(2) * XAI_SHAPE(ISH)
  END DO
END IF
!
!!$GFLAG = .TRUE.
!!$IF (GFLAG) THEN
!!$  WRITE(UNIT=ILUOUT0,FMT='("      snow aggregation")')
!!$  WRITE(UNIT=ILUOUT0,FMT='(" Temp. factor     XCOLEXIS=",E13.6)') XCOLEXIS
!!$END IF
!  
!-------------------------------------------------------------------------------
!
!
!*       7.     CONSTANTS FOR THE FAST COLD PROCESSES FOR THE AGGREGATES
!               --------------------------------------------------------
!
!
!*       7.1    Constants for the riming of the aggregates
!
XDCSLIM  = 0.007 ! D_cs^lim = 7 mm as suggested by Farley et al. (1989)
XCOLCS   = 1.0
XEXCRIMSS= -XDS-2.0
XCRIMSS  = (XPI/4.0)*XCOLCS*XCS*(ZRHO00**XCEXVT)*MOMG(XALPHAS,XNUS,XDS+2.0)

XEXCRIMSG= XEXCRIMSS
XCRIMSG  = XCRIMSS
XSRIMCG  = XAS*MOMG(XALPHAS,XNUS,XBS)
XEXSRIMCG= -XBS
!
!Pour Murakami 1990
XSRIMCG2 = XAG*MOMG(XALPHAS,XNUS,XBG)
XSRIMCG3 = 0.1
XEXSRIMCG2=XBG
!!$GFLAG = .TRUE.
!!$IF (GFLAG) THEN
!!$  WRITE(UNIT=ILUOUT0,FMT='("      riming of the aggregates")')
!!$  WRITE(UNIT=ILUOUT0,FMT='(" D_cs^lim (Farley et al.) XDCSLIM=",E13.6)') XDCSLIM
!!$  WRITE(UNIT=ILUOUT0,FMT='(" Coll. efficiency          XCOLCS=",E13.6)') XCOLCS
!!$END IF
!!$!
PARAM_LIMA_MIXED%NGAMINC = 80
PARAM_LIMA_MIXED%XGAMINC_BOUND_MIN = (1000.*XTRANS_MP_GAMMAS*XDCSLIM)**XALPHAS !1.0E-1 ! Minimal value of (Lbda * D_cs^lim)**alpha
PARAM_LIMA_MIXED%XGAMINC_BOUND_MAX = (50000.*XTRANS_MP_GAMMAS*XDCSLIM)**XALPHAS !1.0E7 ! Maximal value of (Lbda * D_cs^lim)**alpha
ZRATE = EXP(LOG(PARAM_LIMA_MIXED%XGAMINC_BOUND_MAX/PARAM_LIMA_MIXED%XGAMINC_BOUND_MIN)/FLOAT(PARAM_LIMA_MIXED%NGAMINC-1))
!
CALL PARAM_LIMA_MIXED_ALLOCATE('XGAMINC_RIM1', NGAMINC)
CALL PARAM_LIMA_MIXED_ALLOCATE('XGAMINC_RIM2', NGAMINC)
CALL PARAM_LIMA_MIXED_ALLOCATE('XGAMINC_RIM4', NGAMINC)
!
DO I1=1,NGAMINC
  ZBOUND = PARAM_LIMA_MIXED%XGAMINC_BOUND_MIN*ZRATE**(I1-1)
  XGAMINC_RIM1(I1) = GAMMA_INC(XNUS+(2.0+XDS)/XALPHAS,ZBOUND)
  XGAMINC_RIM2(I1) = GAMMA_INC(XNUS+XBS/XALPHAS      ,ZBOUND)
  XGAMINC_RIM4(I1) = GAMMA_INC(XNUS+XBG/XALPHAS      ,ZBOUND) ! Pour Murakami 1990
END DO
!
XRIMINTP1 = XALPHAS / LOG(ZRATE)
XRIMINTP2 = 1.0 + XRIMINTP1*LOG( XDCSLIM/(XGAMINC_BOUND_MIN)**(1.0/XALPHAS) )
!
!*       7.1.1  Defining the constants for the Hallett-Mossop 
!               secondary ice nucleation process
!
XHMTMIN = XTT - 8.0
XHMTMAX = XTT - 3.0
XHM1 = 9.3E-3            ! Obsolete parameterization
XHM2 = 1.5E-3/LOG(10.0)  ! from Ferrier (1995)
XHM_YIELD = 5.E-3 ! A splinter is produced after the riming of 200 droplets
XHM_COLLCS= 1.0   ! Collision efficiency snow/droplet (with Dc>25 microns)
XHM_FACTS = XHM_YIELD*(XHM_COLLCS/XCOLCS)
!
! Notice: One magnitude of lambda discretized over 10 points for the droplets 
!
PARAM_LIMA_MIXED%XGAMINC_HMC_BOUND_MIN = 1.0E-3 ! Min value of (Lbda * (12,25) microns)**alpha
PARAM_LIMA_MIXED%XGAMINC_HMC_BOUND_MAX = 1.0E5  ! Max value of (Lbda * (12,25) microns)**alpha
ZRATE = EXP(LOG(PARAM_LIMA_MIXED%XGAMINC_HMC_BOUND_MAX/PARAM_LIMA_MIXED%XGAMINC_HMC_BOUND_MIN)/REAL(PARAM_LIMA_MIXED%NGAMINC-1))
!
CALL PARAM_LIMA_MIXED_ALLOCATE('XGAMINC_HMC', NGAMINC)
!
DO I1=1,NGAMINC
  ZBOUND = PARAM_LIMA_MIXED%XGAMINC_HMC_BOUND_MIN*ZRATE**(I1-1)
  XGAMINC_HMC(I1) = GAMMA_INC(XNUC,ZBOUND)
END DO
!
XHMSINTP1 = XALPHAC / LOG(ZRATE)
XHMSINTP2 = 1.0 + XHMSINTP1*LOG( 12.E-6/(XGAMINC_HMC_BOUND_MIN)**(1.0/XALPHAC) )
XHMLINTP1 = XALPHAC / LOG(ZRATE)
XHMLINTP2 = 1.0 + XHMLINTP1*LOG( 25.E-6/(XGAMINC_HMC_BOUND_MIN)**(1.0/XALPHAC) )
!
!
!*       7.2    Constants for the accretion of raindrops onto aggregates
!
XFRACCSS  = XPI/4.0*XAR*(ZRHO00**XCEXVT)
XFNRACCSS = XPI/4.0*(ZRHO00**XCEXVT)
!
XLBRACCS1   =    MOMG(XALPHAS,XNUS,2.)*MOMG(XALPHAR,XNUR,3.)
XLBRACCS2   = 2.*MOMG(XALPHAS,XNUS,1.)*MOMG(XALPHAR,XNUR,4.)
XLBRACCS3   =                          MOMG(XALPHAR,XNUR,5.)
XLBNRACCS1  =    MOMG(XALPHAS,XNUS,2.)                             
XLBNRACCS2  = 2.*MOMG(XALPHAS,XNUS,1.)*MOMG(XALPHAR,XNUR,1.)       
XLBNRACCS3  =                          MOMG(XALPHAR,XNUR,2.)
!
XFSACCRG  = (XPI/4.0)*XAS*(ZRHO00**XCEXVT)
XFNSACCRG = (XPI/4.0)*(ZRHO00**XCEXVT)
!
XLBSACCR1   =    MOMG(XALPHAR,XNUR,2.)*MOMG(XALPHAS,XNUS,XBS)
XLBSACCR2   = 2.*MOMG(XALPHAR,XNUR,1.)*MOMG(XALPHAS,XNUS,XBS+1.)
XLBSACCR3   =                          MOMG(XALPHAS,XNUS,XBS+2.)
XLBNSACCR1  =    MOMG(XALPHAR,XNUR,2.)                             
XLBNSACCR2  = 2.*MOMG(XALPHAR,XNUR,1.)*MOMG(XALPHAS,XNUS,1.)      
XLBNSACCR3  =                          MOMG(XALPHAS,XNUS,2.)
!
!*       7.2.1  Defining the ranges for the computation of the kernels
!
! Notice: One magnitude of lambda discretized over 10 points for rain
! Notice: One magnitude of lambda discretized over 10 points for snow
!
PARAM_LIMA_MIXED%NACCLBDAS = 40
PARAM_LIMA_MIXED%XACCLBDAS_MIN = 5.0E1*XTRANS_MP_GAMMAS !5.0E1*XTRANS_MP_GAMMAS ! Minimal value of Lbda_s to tabulate XKER_RACCS
PARAM_LIMA_MIXED%XACCLBDAS_MAX = 5.0E5*XTRANS_MP_GAMMAS !5.0E5*XTRANS_MP_GAMMAS ! Maximal value of Lbda_s to tabulate XKER_RACCS
ZRATE = LOG(XACCLBDAS_MAX/XACCLBDAS_MIN)/FLOAT(NACCLBDAS-1)
PARAM_LIMA_MIXED%XACCINTP1S = 1.0 / ZRATE
PARAM_LIMA_MIXED%XACCINTP2S = 1.0 - LOG( XACCLBDAS_MIN ) / ZRATE
PARAM_LIMA_MIXED%NACCLBDAR = 40
PARAM_LIMA_MIXED%XACCLBDAR_MIN = 1.0E3 ! Minimal value of Lbda_r to tabulate XKER_RACCS
PARAM_LIMA_MIXED%XACCLBDAR_MAX = 1.0E7 ! Maximal value of Lbda_r to tabulate XKER_RACCS
ZRATE = LOG(PARAM_LIMA_MIXED%XACCLBDAR_MAX/PARAM_LIMA_MIXED%XACCLBDAR_MIN)/REAL(NACCLBDAR-1)
XACCINTP1R = 1.0 / ZRATE
XACCINTP2R = 1.0 - LOG( PARAM_LIMA_MIXED%XACCLBDAR_MIN ) / ZRATE
!
!*       7.2.2  Computations of the tabulated normalized kernels
!
IND      = 50    ! Interval number, collection efficiency and infinite diameter
ZZESR     = 1.0   ! factor used to integrate the dimensional distributions when
ZZFDINFTY = 20.0  ! computing the kernels XKER_RACCSS, XKER_RACCS and XKER_SACCRG
!
CALL PARAM_LIMA_MIXED_ALLOCATE('XKER_RACCSS', NACCLBDAS,NACCLBDAR)
CALL PARAM_LIMA_MIXED_ALLOCATE('XKER_RACCS', NACCLBDAS,NACCLBDAR)
CALL PARAM_LIMA_MIXED_ALLOCATE('XKER_SACCRG', NACCLBDAR,NACCLBDAS)
CALL PARAM_LIMA_MIXED_ALLOCATE('XKER_N_RACCSS', NACCLBDAS,NACCLBDAR)
CALL PARAM_LIMA_MIXED_ALLOCATE('XKER_N_RACCS', NACCLBDAS,NACCLBDAR)
CALL PARAM_LIMA_MIXED_ALLOCATE('XKER_N_SACCRG', NACCLBDAR,NACCLBDAS)
IF (NMOM_R.GE.2 .AND. .FALSE.) THEN 
CALL NRCOLSS ( NACCLBDAS,NACCLBDAR,IND, XALPHAS, XNUS, XALPHAR, XNUR,     & 
               ZZESR, XCS, XDS, XFVELOS, XCR, XDR,                        & 
               XACCLBDAS_MAX, XACCLBDAR_MAX, XACCLBDAS_MIN, XACCLBDAR_MIN, & 
               ZZFDINFTY, XKER_N_RACCSS, XAG, XBS, XAS                      )
ELSE
   XKER_N_RACCSS(:,:)=0.
END IF
IF (NMOM_R.GE.2) THEN 
CALL NZCOLX  ( NACCLBDAS,NACCLBDAR, IND, XALPHAS, XNUS, XALPHAR, XNUR,      & 
               ZZESR, XCS, XDS, XFVELOS, XCR, XDR, 0.,                      & ! 
               XACCLBDAS_MAX, XACCLBDAR_MAX, XACCLBDAS_MIN, XACCLBDAR_MIN, & 
               ZZFDINFTY, XKER_N_RACCS                                      )
ELSE
   XKER_N_RACCS(:,:)=0.
END IF
IF (NMOM_S.GE.2 .OR. NMOM_G.GE.2) THEN 
CALL NSCOLRG ( NACCLBDAR,NACCLBDAS, IND, XALPHAS, XNUS, XALPHAR, XNUR,     & 
               ZZESR, XCS, XDS, XFVELOS, XCR, XDR,                          & 
               XACCLBDAS_MAX, XACCLBDAR_MAX, XACCLBDAS_MIN, XACCLBDAR_MIN, & 
               ZZFDINFTY, XKER_N_SACCRG,XAG, XBS, XAS                       )
ELSE
   XKER_N_SACCRG(:,:)=0.
END IF
!!$WRITE(UNIT=ILUOUT0,FMT='("!")')
!!$WRITE(UNIT=ILUOUT0,FMT='("IF( PRESENT(PKER_N_RACCSS) ) THEN")')          
!!$DO I1 = 1 , NACCLBDAS                                                    
!!$  DO J2 = 1 , NACCLBDAR                                                  
!!$  WRITE(UNIT=ILUOUT0,FMT='("  PKER_N_RACCSS(",I3,",",I3,") = ",E13.6)') &
!!$                    I1,J2,XKER_N_RACCSS(I1,J2)                          
!!$  END DO                                                                
!!$END DO                                                                   
!!$WRITE(UNIT=ILUOUT0,FMT='("!")')                                           
!!$WRITE(UNIT=ILUOUT0,FMT='("!")')
!!$WRITE(UNIT=ILUOUT0,FMT='("IF( PRESENT(PKER_N_RACCS) ) THEN")')          
!!$DO I1 = 1 , NACCLBDAS                                                    
!!$  DO J2 = 1 , NACCLBDAR                                                  
!!$  WRITE(UNIT=ILUOUT0,FMT='("  PKER_N_RACCS(",I3,",",I3,") = ",E13.6)') &
!!$                    I1,J2,XKER_N_RACCS(I1,J2)                          
!!$  END DO                                                                
!!$END DO                                                                   
!!$WRITE(UNIT=ILUOUT0,FMT='("!")') 
!!$WRITE(UNIT=ILUOUT0,FMT='("!")')
!!$WRITE(UNIT=ILUOUT0,FMT='("IF( PRESENT(PKER_N_SACCRG) ) THEN")')          
!!$DO I1 = 1 , NACCLBDAR                                                    
!!$  DO J2 = 1 , NACCLBDAS                                                  
!!$  WRITE(UNIT=ILUOUT0,FMT='("  PKER_N_SACCRG(",I3,",",I3,") = ",E13.6)') &
!!$                    I1,J2,XKER_N_SACCRG(I1,J2)                          
!!$  END DO                                                                
!!$END DO                                                                   
!!$WRITE(UNIT=ILUOUT0,FMT='("!")') 

!
CALL LIMA_READ_XKER_RACCS (IACCLBDAS,IACCLBDAR,IKND,                                        &
                           ZALPHAS,ZNUS,ZALPHAR,ZNUR,ZESR,ZBS,ZBR,ZCS,ZDS,ZFVELOS,ZCR,ZDR, &
                           ZACCLBDAS_MAX,ZACCLBDAR_MAX,ZACCLBDAS_MIN,ZACCLBDAR_MIN,        &
                           ZFDINFTY                                                        )
IF( (IACCLBDAS/=NACCLBDAS) .OR. (IACCLBDAR/=NACCLBDAR) .OR. (IKND/=IND) .OR. &
    (ZALPHAS/=XALPHAS) .OR. (ZNUS/=XNUS)                               .OR. &
    (ZALPHAR/=XALPHAR) .OR. (ZNUR/=XNUR)                               .OR. &
    (ZESR/=ZZESR) .OR. (ZBS/=XBS) .OR. (ZBR/=XBR)                       .OR. &
    (ZCS/=XCS) .OR. (ZDS/=XDS) .OR. (ZFVELOS/=XFVELOS) .OR. (ZCR/=XCR) .OR. (ZDR/=XDR)         .OR. &
    (ABS(ZACCLBDAS_MAX-XACCLBDAS_MAX).GE.(0.001*ZACCLBDAS_MAX)) .OR. (ZACCLBDAR_MAX/=XACCLBDAR_MAX) .OR. &
    (ABS(ZACCLBDAS_MIN-XACCLBDAS_MIN).GE.(0.001*ZACCLBDAS_MIN)) .OR. (ZACCLBDAR_MIN/=XACCLBDAR_MIN) .OR. &
    (ZFDINFTY/=ZZFDINFTY)                                               ) THEN
  CALL LIMA_RRCOLSS ( IACCLBDAS, IACCLBDAR, IND, XALPHAS, XNUS, XALPHAR, XNUR, &
                 ZZESR, XBR, XCS, XDS, XFVELOS, XCR, XDR,                     &
                 XACCLBDAS_MAX, XACCLBDAR_MAX, XACCLBDAS_MIN, XACCLBDAR_MIN, &
                 ZZFDINFTY, XKER_RACCSS, XAG, XBS, XAS                        )
  CALL LIMA_RZCOLX  ( IACCLBDAS, IACCLBDAR, IND, XALPHAS, XNUS, XALPHAR, XNUR, &
                 ZZESR, XBR, XCS, XDS, XFVELOS, XCR, XDR, 0.,                 &
                 XACCLBDAS_MAX, XACCLBDAR_MAX, XACCLBDAS_MIN, XACCLBDAR_MIN, &
                 ZZFDINFTY, XKER_RACCS                                        )
  CALL LIMA_RSCOLRG ( IACCLBDAS, IACCLBDAR, IND, XALPHAS, XNUS, XALPHAR, XNUR, &
                 ZZESR, XBS, XCS, XDS, XFVELOS, XCR, XDR,                     &
                 XACCLBDAS_MAX, XACCLBDAR_MAX, XACCLBDAS_MIN, XACCLBDAR_MIN, &
                 ZZFDINFTY, XKER_SACCRG,XAG, XBS, XAS                         )
!!$  WRITE(UNIT=ILUOUT0,FMT='("*****************************************")')
!!$  WRITE(UNIT=ILUOUT0,FMT='("**** UPDATE NEW SET OF RACSS KERNELS ****")')
!!$  WRITE(UNIT=ILUOUT0,FMT='("**** UPDATE NEW SET OF RACS  KERNELS ****")')
!!$  WRITE(UNIT=ILUOUT0,FMT='("**** UPDATE NEW SET OF SACRG KERNELS ****")')
!!$  WRITE(UNIT=ILUOUT0,FMT='("*****************************************")')
!!$  WRITE(UNIT=ILUOUT0,FMT='("!")')
!!$  WRITE(UNIT=ILUOUT0,FMT='("IKND=",I3)') IND
!!$  WRITE(UNIT=ILUOUT0,FMT='("IACCLBDAS=",I3)') NACCLBDAS
!!$  WRITE(UNIT=ILUOUT0,FMT='("IACCLBDAR=",I3)') NACCLBDAR
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZALPHAS=",E13.6)') XALPHAS
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZNUS=",E13.6)') XNUS
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZALPHAR=",E13.6)') XALPHAR
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZNUR=",E13.6)') XNUR
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZESR=",E13.6)') ZZESR
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZBS=",E13.6)') XBS
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZBR=",E13.6)') XBR
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZCS=",E13.6)') XCS
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZDS=",E13.6)') XDS
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZFVELOS=",E13.6)') XFVELOS
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZCR=",E13.6)') XCR
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZDR=",E13.6)') XDR
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZACCLBDAS_MAX=",E13.6)') &
!!$                                                    XACCLBDAS_MAX
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZACCLBDAR_MAX=",E13.6)') &
!!$                                                    XACCLBDAR_MAX
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZACCLBDAS_MIN=",E13.6)') &
!!$                                                    XACCLBDAS_MIN
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZACCLBDAR_MIN=",E13.6)') &
!!$                                                    XACCLBDAR_MIN
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZFDINFTY=",E13.6)') ZZFDINFTY
!!$  WRITE(UNIT=ILUOUT0,FMT='("!")')
!!$  WRITE(UNIT=ILUOUT0,FMT='("IF( PRESENT(PKER_RACCSS) ) THEN")')
!!$  DO I1 = 1 , NACCLBDAS
!!$    DO J2 = 1 , NACCLBDAR
!!$    WRITE(UNIT=ILUOUT0,FMT='("  PKER_RACCSS(",I3,",",I3,") = ",E13.6)') &
!!$                        I1,J2,XKER_RACCSS(I1,J2)
!!$    END DO
!!$  END DO
!!$  WRITE(UNIT=ILUOUT0,FMT='("END IF")')
!!$  WRITE(UNIT=ILUOUT0,FMT='("!")')
!!$  WRITE(UNIT=ILUOUT0,FMT='("IF( PRESENT(PKER_RACCS ) ) THEN")')
!!$  DO I1 = 1 , NACCLBDAS
!!$    DO J2 = 1 , NACCLBDAR
!!$    WRITE(UNIT=ILUOUT0,FMT='("  PKER_RACCS (",I3,",",I3,") = ",E13.6)') &
!!$                        I1,J2,XKER_RACCS (I1,J2)
!!$    END DO
!!$  END DO
!!$  WRITE(UNIT=ILUOUT0,FMT='("END IF")')
!!$  WRITE(UNIT=ILUOUT0,FMT='("!")')
!!$  WRITE(UNIT=ILUOUT0,FMT='("IF( PRESENT(PKER_SACCRG) ) THEN")')
!!$  DO I1 = 1 , NACCLBDAR
!!$    DO J2 = 1 , NACCLBDAS
!!$    WRITE(UNIT=ILUOUT0,FMT='("  PKER_SACCRG(",I3,",",I3,") = ",E13.6)') &
!!$                        I1,J2,XKER_SACCRG(I1,J2)
!!$    END DO
!!$  END DO
!!$  WRITE(UNIT=ILUOUT0,FMT='("END IF")')
  ELSE
  CALL LIMA_READ_XKER_RACCS (IACCLBDAS,IACCLBDAR,IKND,                                        &
                             ZALPHAS,ZNUS,ZALPHAR,ZNUR,ZESR,ZBS,ZBR,ZCS,ZDS,ZFVELOS,ZCR,ZDR, &
                             ZACCLBDAS_MAX,ZACCLBDAR_MAX,ZACCLBDAS_MIN,ZACCLBDAR_MIN,        &
                             ZFDINFTY,XKER_RACCSS,XKER_RACCS,XKER_SACCRG                     )
!!$  WRITE(UNIT=ILUOUT0,FMT='(" Read XKER_RACCSS")')
!!$  WRITE(UNIT=ILUOUT0,FMT='(" Read XKER_RACCS ")')
!!$  WRITE(UNIT=ILUOUT0,FMT='(" Read XKER_SACCRG")')
END IF
!
!*       7.2N  Computations of the tabulated normalized kernels Snow Self Collection !! 
!
!
!if (NMOM_S.GE.2) then
  XCOLSS   = 0.1 ! Collection efficiency of S+S
  XCOLEXSS = 0.1  ! Temperature factor of the S+S collection efficiency
  XFNSSCS = (XPI/4.0)*XCOLSS*(ZRHO00**XCEXVT)/2.0       
!
  XLBNSSCS1   = 2.0*MOMG(XALPHAS,XNUS,2.)                        
  XLBNSSCS2   = 2.0*MOMG(XALPHAS,XNUS,1.)**2  
!
!*       7.2N.1  Defining the ranges for the computation of the kernels
!
! Notice: One magnitude of lambda discretized over 10 points for snow
!
PARAM_LIMA_COLD%NSCLBDAS = 80
PARAM_LIMA_MIXED%XSCLBDAS_MIN = 1.0E0*XTRANS_MP_GAMMAS  ! Minimal value of Lbda_s to tabulate XKER_RSCS
PARAM_LIMA_MIXED%XSCLBDAS_MAX = 5.0E10*XTRANS_MP_GAMMAS ! Maximal value of Lbda_s to tabulate XKER_RSCS
ZRATE = LOG(PARAM_LIMA_MIXED%XSCLBDAS_MAX/PARAM_LIMA_MIXED%XSCLBDAS_MIN)/FLOAT(PARAM_LIMA_COLD%NSCLBDAS-1)
XSCINTP1S = 1.0 / ZRATE
XSCINTP2S = 1.0 - LOG( XSCLBDAS_MIN ) / ZRATE
!
!
  IND      = 50    ! Interval number, collection efficiency and infinite diameter
  ZZESS     = 1.0   ! factor used to integrate the dimensional distributions when
  ZZFDINFTY = 20.0  ! computing the kernels XKER_SSCSS
!
  CALL PARAM_LIMA_COLD_ALLOCATE('XKER_N_SSCS', NSCLBDAS,NSCLBDAS)
!
IF (NMOM_S.GE.2) THEN 
  CALL NZCOLX  ( NSCLBDAS,NSCLBDAS,IND, XALPHAS, XNUS, XALPHAS, XNUS,      & 
               ZZESS, XCS, XDS, XFVELOS, XCS, XDS, XFVELOS,                   & 
               XSCLBDAS_MAX, XSCLBDAS_MAX, XSCLBDAS_MIN, XSCLBDAS_MIN, & 
               ZZFDINFTY, XKER_N_SSCS                                      ) 
ELSE
   XKER_N_SSCS(:,:)=0.
END IF
!
!!$  WRITE(UNIT=ILUOUT0,FMT='("*****************************************")')
!!$  WRITE(UNIT=ILUOUT0,FMT='("**** UPDATE NEW SET OF SSCS  KERNELS ***")')
!!$  WRITE(UNIT=ILUOUT0,FMT='("*****************************************")')
!!$  WRITE(UNIT=ILUOUT0,FMT='("!")')
!!$  WRITE(UNIT=ILUOUT0,FMT='("IKND=",I3)') IND
!!$  WRITE(UNIT=ILUOUT0,FMT='("KSCLBDAS=",I3)') NSCLBDAS
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZALPHAS=",E13.6)') XALPHAS
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZNUS=",E13.6)') XNUS
!!$  WRITE(UNIT=ILUOUT0,FMT='("PESS=",E13.6)') ZZESS
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZBS=",E13.6)') XBS
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZCS=",E13.6)') XCS
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZDS=",E13.6)') XDS
!!$  WRITE(UNIT=ILUOUT0,FMT='("PSCLBDAS_MAX=",E13.6)') &
!!$                                                  XSCLBDAS_MAX
!!$  WRITE(UNIT=ILUOUT0,FMT='("PSCLBDAS_MIN=",E13.6)') &
!!$                                                  XSCLBDAS_MIN
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZFDINFTY=",E13.6)') ZZFDINFTY
!!$  WRITE(UNIT=ILUOUT0,FMT='("!")')
!!$  WRITE(UNIT=ILUOUT0,FMT='("!")')
!!$  WRITE(UNIT=ILUOUT0,FMT='("IF( PRESENT(PKER_N_SSCS ) ) THEN")')           
!!$  DO I1 = 1 , NSCLBDAS                                                     
!!$    WRITE(UNIT=ILUOUT0,FMT='("  PKER_N_SSCS (",I3,",",I3,") = ",E13.6)') & 
!!$                      I1,I1,XKER_N_SSCS (I1,I1)                          
!!$  END DO                                                                    
!!$  WRITE(UNIT=ILUOUT0,FMT='("!")')  
!
!*       7.2N2    Constants for the 'spontaneous' break-up
!
  XACCS1=MOMG(XALPHAS,XNUS,XBS)
  XSPONBUDS1 = 2.5E-3
  XSPONBUDS2 = 3.5E-3
  XSPONBUDS3 = 4.5E-3
  XSPONCOEFS2 = ((XSPONBUDS3/XSPONBUDS2)**3 - 1.0)/(XSPONBUDS3-XSPONBUDS1)**2
!
!end if
!
!*       7.2TER  Computations of the tabulated normalized kernels ice Self Collection !! 
!
!
IF (NMOM_S.GE.2) THEN
XCOLII   = 0.005 ! Collection efficiency of I+I
XCOLEXII = 0.1  ! Temperature factor of the I+I collection efficiency
XFNISCI = (XPI/4.0)*(ZRHO00**XCEXVT)/2.0       
!
XLBNISCI1   = 2.0*MOMG(XALPHAI,XNUI,2.)                        
XLBNISCI2   = 2.0*MOMG(XALPHAI,XNUI,1.)**2 
!
!
!*       7.2TER.1  Defining the ranges for the computation of the kernels
!
!
PARAM_LIMA_COLD%NSCLBDAI = 80
PARAM_LIMA_COLD%XSCLBDAI_MIN = 1.0E0
PARAM_LIMA_COLD%XSCLBDAI_MAX = 5.0E10
ZRATE = LOG(PARAM_LIMA_COLD%XSCLBDAI_MAX/PARAM_LIMA_COLD%XSCLBDAI_MIN)/FLOAT(PARAM_LIMA_COLD%NSCLBDAI-1)
XSCINTP1I = 1.0 / ZRATE
XSCINTP2I = 1.0 - LOG( XSCLBDAI_MIN ) / ZRATE
!
!
IND      = 50    ! Interval number, collection efficiency and infinite diameter
ZZEII     = 1.0   ! factor used to integrate the dimensional distributions when
ZZFDINFTY = 20.0  ! computing the kernels XKER_SSCSS
!
CALL PARAM_LIMA_COLD_ALLOCATE('XKER_N_ISCI', NSCLBDAI,NSCLBDAI)
!
CALL NZCOLX  ( NSCLBDAI, NSCLBDAI, IND, XALPHAI, XNUI, XALPHAI, XNUI,&
             ZZEII, XC_I, XDI, 0., XC_I, XDI, 0.,                    &
             XSCLBDAI_MAX, XSCLBDAI_MAX, XSCLBDAI_MIN, XSCLBDAI_MIN, &
             ZZFDINFTY, XKER_N_ISCI                                  )
!! 
IF(LCRYSTAL_SHAPE) THEN
    CALL PARAM_LIMA_COLD_ALLOCATE_1D('XLBISCS1',NNB_CRYSTAL_SHAPE)
    CALL PARAM_LIMA_COLD_ALLOCATE_1D('XLBISCS2',NNB_CRYSTAL_SHAPE)
    CALL PARAM_LIMA_COLD_ALLOCATE_1D('XLBISCS3',NNB_CRYSTAL_SHAPE)
    CALL PARAM_LIMA_COLD_ALLOCATE_1D('XFISCS',NNB_CRYSTAL_SHAPE)
   
    CALL PARAM_LIMA_COLD_ALLOCATE_3D('XKER_R_ISCI_SHAPE',NNB_CRYSTAL_SHAPE, NSCLBDAI,NSCLBDAI)
    CALL PARAM_LIMA_COLD_ALLOCATE_4D('XKER_N_ISCI_SHAPE',NNB_CRYSTAL_SHAPE,NNB_CRYSTAL_SHAPE, NSCLBDAI,NSCLBDAI)    
    DO ISH = 1, NNB_CRYSTAL_SHAPE
        XFISCS(ISH) = XAI_SHAPE(3)*(XPI/4.0)*(ZRHO00**XCEXVT)/(MOMG(XALPHAI,XNUI,XBI_SHAPE(ISH)))

        XLBISCS1(ISH) =    MOMG(XALPHAI,XNUI,2.)*MOMG(XALPHAI,XNUI,XBI_SHAPE(ISH))
        XLBISCS2(ISH) = 2.*MOMG(XALPHAI,XNUI,1.)*MOMG(XALPHAI,XNUI,XBI_SHAPE(ISH)+1.)
        XLBISCS3(ISH) =                          MOMG(XALPHAI,XNUI,XBI_SHAPE(ISH)+2.)   

        CALL LIMA_RZCOLX  ( NSCLBDAI, NSCLBDAI, IND, XALPHAI, XNUI, XALPHAI, XNUI,                       &
               ZZEII, XBI_SHAPE(3), XC_I_SHAPE(ISH), XDI_SHAPE(ISH), 0., XC_I_SHAPE(3), XDI_SHAPE(3), 0., &
               XSCLBDAI_MAX, XSCLBDAI_MAX, XSCLBDAI_MIN, XSCLBDAI_MIN,                                 &
               ZZFDINFTY, XKER_R_ISCI_SHAPE(ISH,:,:)                                                    )
        DO ISH2 = 1, NNB_CRYSTAL_SHAPE
            CALL NZCOLX  ( NSCLBDAI, NSCLBDAI, IND, XALPHAI, XNUI, XALPHAI, XNUI,                 &
               ZZEII, XC_I_SHAPE(ISH), XDI_SHAPE(ISH), 0., XC_I_SHAPE(ISH2), XDI_SHAPE(ISH2), 0., &
               XSCLBDAI_MAX, XSCLBDAI_MAX, XSCLBDAI_MIN, XSCLBDAI_MIN,                         &
               ZZFDINFTY, XKER_N_ISCI_SHAPE(ISH,ISH2,:,:)                                       )
        ENDDO
    ENDDO
ENDIF
ENDIF
!
!*       7.3    Constant for the conversion-melting rate
!
XFSCVMG = 2.0
!
!!$GFLAG = .TRUE.
!!$IF (GFLAG) THEN
!!$  WRITE(UNIT=ILUOUT0,FMT='("      conversion-melting of the aggregates")')
!!$  WRITE(UNIT=ILUOUT0,FMT='(" Conv. factor XFSCVMG=",E13.6)') XFSCVMG
!!$END IF
!
!
!*       7.4 Constants for Ice-Ice collision process (CIBU)
!
XDCSLIM_CIBU_MIN = 2.0E-4 ! D_cs lim min
XDCSLIM_CIBU_MAX = 1.0E-3 ! D_cs lim max
XDCGLIM_CIBU_MIN = 2.0E-3 ! D_cg lim min
!
!!$GFLAG = .TRUE.
!!$IF (GFLAG) THEN
!!$  WRITE(UNIT=ILUOUT0,FMT='(" Ice-ice collision process")')
!!$  WRITE(UNIT=ILUOUT0,FMT='(" D_cs^lim min-max =",E13.6)') XDCSLIM_CIBU_MIN,XDCSLIM_CIBU_MAX
!!$  WRITE(UNIT=ILUOUT0,FMT='(" D_cg^lim min     =",E13.6)') XDCGLIM_CIBU_MIN
!!$END IF
!
PARAM_LIMA_MIXED%NGAMINC = 80
!
!Note : Boundaries are rounded at 5.0 or 1.0 (down for Bound_min and up for Bound_max)
PARAM_LIMA_MIXED%XGAMINC_BOUND_CIBU_SMIN = 1.0E-5 * XTRANS_MP_GAMMAS**XALPHAS ! Minimal value of (Lbda_s * D_cs^lim)**alpha) 0.2 mm
PARAM_LIMA_MIXED%XGAMINC_BOUND_CIBU_SMAX = 5.0E+2 * XTRANS_MP_GAMMAS**XALPHAS ! Maximal value of (Lbda_s * D_cs^lim)**alpha) 1 mm
ZRATE_S = EXP(LOG(PARAM_LIMA_MIXED%XGAMINC_BOUND_CIBU_SMAX/PARAM_LIMA_MIXED%XGAMINC_BOUND_CIBU_SMIN)/&
              FLOAT(PARAM_LIMA_MIXED%NGAMINC-1))
!
PARAM_LIMA_MIXED%XGAMINC_BOUND_CIBU_GMIN = 1.0E-1 ! Minimal value of (Lbda_g * D_cg^lim)**alpha) 2 mm
PARAM_LIMA_MIXED%XGAMINC_BOUND_CIBU_GMAX = 5.0E+1 ! Maximal value of (Lbda_g * D_cg^lim)**alpha) 2 mm
ZRATE_G = EXP(LOG(PARAM_LIMA_MIXED%XGAMINC_BOUND_CIBU_GMAX/PARAM_LIMA_MIXED%XGAMINC_BOUND_CIBU_GMIN)/&
              FLOAT(PARAM_LIMA_MIXED%NGAMINC-1))
!
CALL PARAM_LIMA_MIXED_ALLOCATE('XGAMINC_CIBU_S', 4,PARAM_LIMA_MIXED%NGAMINC)
CALL PARAM_LIMA_MIXED_ALLOCATE('XGAMINC_CIBU_G', 2,PARAM_LIMA_MIXED%NGAMINC)
!
DO I1 = 1, NGAMINC
  ZBOUND_CIBU_SMIN = PARAM_LIMA_MIXED%XGAMINC_BOUND_CIBU_SMIN * ZRATE_S**(I1-1)
  ZBOUND_CIBU_GMIN = PARAM_LIMA_MIXED%XGAMINC_BOUND_CIBU_GMIN * ZRATE_G**(I1-1)
!
! For ZNI_CIBU
  XGAMINC_CIBU_S(1,I1) = GAMMA_INC(XNUS,ZBOUND_CIBU_SMIN)
  XGAMINC_CIBU_S(2,I1) = GAMMA_INC(XNUS+(XDS/XALPHAS),ZBOUND_CIBU_SMIN)
!
  XGAMINC_CIBU_G(1,I1) = GAMMA_INC(XNUG+((2.0+XDG)/XALPHAG),ZBOUND_CIBU_GMIN)
  XGAMINC_CIBU_G(2,I1) = GAMMA_INC(XNUG+(2.0/XALPHAG),ZBOUND_CIBU_GMIN)
!
! For ZRI_CIBU
  XGAMINC_CIBU_S(3,I1) = GAMMA_INC(XNUS+(XBS/XALPHAS),ZBOUND_CIBU_SMIN)
  XGAMINC_CIBU_S(4,I1) = GAMMA_INC(XNUS+((XBS+XDS)/XALPHAS),ZBOUND_CIBU_SMIN)
END DO
!
XCIBUINTP_S = XALPHAS / LOG(ZRATE_S)
XCIBUINTP1_S = 1.0 + XCIBUINTP_S * LOG(XDCSLIM_CIBU_MIN/(XGAMINC_BOUND_CIBU_SMIN)**(1.0/XALPHAS))
XCIBUINTP2_S = 1.0 + XCIBUINTP_S * LOG(XDCSLIM_CIBU_MAX/(XGAMINC_BOUND_CIBU_SMIN)**(1.0/XALPHAS))
!
XCIBUINTP_G = XALPHAG / LOG(ZRATE_G)
XCIBUINTP1_G = 1.0 + XCIBUINTP_G * LOG(XDCGLIM_CIBU_MIN/(XGAMINC_BOUND_CIBU_GMIN)**(1.0/XALPHAG))
!
! For ZNI_CIBU
XFACTOR_CIBU_NI = (XPI / 4.0) * (ZRHO00**XCEXVT)
XMOMGG_CIBU_1 = MOMG(XALPHAG,XNUG,2.0+XDG)
XMOMGG_CIBU_2 = MOMG(XALPHAG,XNUG,2.0)
XMOMGS_CIBU_1 = MOMG(XALPHAS,XNUS,XDS)
!
! For ZRI_CIBU
XFACTOR_CIBU_RI = XAS * (XPI / 4.0) * (ZRHO00**XCEXVT)
XMOMGS_CIBU_2 = MOMG(XALPHAS,XNUS,XBS)
XMOMGS_CIBU_3 = MOMG(XALPHAS,XNUS,XBS+XDS)
!
!
!*       7.5 Constants for raindrop shattering by freezing process (RDSF)
!
XDCRLIM_RDSF_MIN = 0.1E-3 ! D_cr lim min
!
!!$GFLAG = .TRUE.
!!$IF (GFLAG) THEN
!!$  WRITE(UNIT=ILUOUT0,FMT='(" Ice-rain collision process")')
!!$  WRITE(UNIT=ILUOUT0,FMT='(" D_cr^lim min     =",E13.6)') XDCRLIM_RDSF_MIN
!!$END IF
!
PARAM_LIMA_MIXED%NGAMINC = 80
!
PARAM_LIMA_MIXED%XGAMINC_BOUND_RDSF_RMIN = 1.0E-5    ! Minimal value of (Lbda_r * D_cr^lim)**alpha) 0.1 mm
PARAM_LIMA_MIXED%XGAMINC_BOUND_RDSF_RMAX = 5.0E-3    ! Maximal value of (Lbda_r * D_cr^lim)**alpha) 1 mm
ZRATE_R = EXP(LOG(PARAM_LIMA_MIXED%XGAMINC_BOUND_RDSF_RMAX/PARAM_LIMA_MIXED%XGAMINC_BOUND_RDSF_RMIN)/&
              FLOAT(PARAM_LIMA_MIXED%NGAMINC-1))
!
CALL PARAM_LIMA_MIXED_ALLOCATE('XGAMINC_RDSF_R', NGAMINC)
!
DO I1 = 1, NGAMINC
  ZBOUND_RDSF_RMIN = PARAM_LIMA_MIXED%XGAMINC_BOUND_RDSF_RMIN * ZRATE_R**(I1-1)
!
! For ZNI_RDSF
  XGAMINC_RDSF_R(I1) = GAMMA_INC(XNUR+((6.0+XDR)/XALPHAR),ZBOUND_RDSF_RMIN)
END DO
!
XRDSFINTP_R = XALPHAR / LOG(ZRATE_R)
XRDSFINTP1_R = 1.0 + XRDSFINTP_R * LOG( XDCRLIM_RDSF_MIN/(XGAMINC_BOUND_RDSF_RMIN)**(1.0/XALPHAR) )
!
! For ZNI_RDSF
ZKHI_LWM = 2.5E13 ! Coeff. in Lawson-Woods-Morrison for the number of splinters
                  ! N_DF = XKHI_LWM * D_R**4
XFACTOR_RDSF_NI = ZKHI_LWM * (XPI / 4.0) * XCR * (ZRHO00**XCEXVT)
XMOMGR_RDSF = MOMG(XALPHAR,XNUR,6.0+XDR)
!
XTM_PSH = 258.     ! Mean Temperature of shattering probability normal distribution
XSIG_PSH = 5.      ! SIGMA of shattering probability normal distribution
!
!-------------------------------------------------------------------------------
!
!
!*       8.     CONSTANTS FOR THE FAST COLD PROCESSES FOR THE GRAUPELN
!               --------------------------------------------------------
!
!
!*       8.1    Constants for the rain contact freezing
!
XCOLIR    = 1.0
!
! values of these coeficients differ from the single-momemt rain_ice case
!
XEXRCFRI  = -XDR-5.0
XRCFRI    = ((XPI**2)/24.0)*XRHOLW*XCOLIR*XCR*(ZRHO00**XCEXVT)     &
                                                     *MOMG(XALPHAR,XNUR,XDR+5.0)
XEXICFRR  = -XDR-2.0
XICFRR    = (XPI/4.0)*XCOLIR*XCR*(ZRHO00**XCEXVT)                  &
                                                     *MOMG(XALPHAR,XNUR,XDR+2.0)
!
!!$GFLAG = .TRUE.
!!$IF (GFLAG) THEN
!!$  WRITE(UNIT=ILUOUT0,FMT='("      rain contact freezing")')
!!$  WRITE(UNIT=ILUOUT0,FMT='(" Coll. efficiency          XCOLIR=",E13.6)') XCOLIR
!!$END IF
!
!
!*       8.2    Constants for the dry growth of the graupeln
!
! Min val. of R for dry and wet for D = 125 mum
XMINDG  =  XAG*XCCG*XMVDMIN_G**(XBG-XCXG)*MOMG(XALPHAG,XNUG,XBG)**(XCXG/XBG)
! Max val of lambda for D = 125 mum
XMAXLBDG = MOMG(XALPHAG,XNUG,XBG)**(1./XBG)/XMVDMIN_G
!
!*       8.2.1  Constants for the cloud droplet collection by the graupeln
!               and for the Hallett-Mossop process
!
XCOLCG  = 0.6  !  Estimated from Cober and List (1993)
XFCDRYG = (XPI/4.0)*XCOLCG*XCG*(ZRHO00**XCEXVT)*MOMG(XALPHAG,XNUG,XDG+2.0)
!
XHM_COLLCG= 0.9   ! Collision efficiency graupel/droplet (with Dc>25 microns)
XHM_FACTG = XHM_YIELD*(XHM_COLLCG/XCOLCG)
!
!*       8.2.2  Constants for the cloud ice collection by the graupeln
!
XCOLIG    = 0.25 ! Collection efficiency of I+G
XCOLEXIG  = 0.05 ! Temperature factor of the I+G collection efficiency
XCOLIG   = 0.01 ! Collection efficiency of I+G
XCOLEXIG = 0.1  ! Temperature factor of the I+G collection efficiency
!!$WRITE (ILUOUT0, FMT=*) ' NEW Constants for the cloud ice collection by the graupeln'
!!$WRITE (ILUOUT0, FMT=*) ' XCOLIG, XCOLEXIG  = ',XCOLIG,XCOLEXIG
XFIDRYG = (XPI/4.0)*XCOLIG*XCG*(ZRHO00**XCEXVT)*MOMG(XALPHAG,XNUG,XDG+2.0)
!
!!$GFLAG = .TRUE.
!!$IF (GFLAG) THEN
!!$  WRITE(UNIT=ILUOUT0,FMT='("      cloud ice collection by the graupeln")')
!!$  WRITE(UNIT=ILUOUT0,FMT='(" Coll. efficiency XCOLIG=",E13.6)') XCOLIG
!!$  WRITE(UNIT=ILUOUT0,FMT='(" Temp. factor     XCOLEXIG=",E13.6)') XCOLEXIG
!!$END IF
!
!*       8.2.3  Constants for the aggregate collection by the graupeln
!
XCOLSG    = 0.25 ! Collection efficiency of S+G
XCOLEXSG  = 0.05 ! Temperature factor of the S+G collection efficiency
XCOLSG   = 0.01 ! Collection efficiency of S+G
XCOLEXSG = 0.1  ! Temperature factor of the S+G collection efficiency
!!$WRITE (ILUOUT0, FMT=*) ' NEW Constants for the aggregate collection by the graupeln'
!!$WRITE (ILUOUT0, FMT=*) ' XCOLSG, XCOLEXSG  = ',XCOLSG,XCOLEXSG
XFSDRYG = XNS*(XPI/4.0)*XCOLSG*XAS*(ZRHO00**XCEXVT)
XFNSDRYG= (XPI/4.0)*XCOLSG*(ZRHO00**XCEXVT)
!
XLBNSDRYG1  =    MOMG(XALPHAG,XNUG,2.)                        
XLBNSDRYG2  = 2.*MOMG(XALPHAG,XNUG,1.)*MOMG(XALPHAS,XNUS,1.) 
XLBNSDRYG3  =                          MOMG(XALPHAS,XNUS,2.)
XLBSDRYG1   =    MOMG(XALPHAG,XNUG,2.)*MOMG(XALPHAS,XNUS,XBS)
XLBSDRYG2   = 2.*MOMG(XALPHAG,XNUG,1.)*MOMG(XALPHAS,XNUS,XBS+1.)
XLBSDRYG3   =                          MOMG(XALPHAS,XNUS,XBS+2.)
!
!!$GFLAG = .TRUE.
!!$IF (GFLAG) THEN
!!$  WRITE(UNIT=ILUOUT0,FMT='("      aggregate collection by the graupeln")')
!!$  WRITE(UNIT=ILUOUT0,FMT='(" Coll. efficiency XCOLSG=",E13.6)') XCOLSG
!!$  WRITE(UNIT=ILUOUT0,FMT='(" Temp. factor     XCOLEXSG=",E13.6)') XCOLEXSG
!!$END IF
!
!*       8.2.4  Constants for the raindrop collection by the graupeln
!
XFRDRYG = XNR*(XPI/4.0)*XAR*(ZRHO00**XCEXVT)
XFNRDRYG= (XPI/4.0)*(ZRHO00**XCEXVT)
!
XLBRDRYG1   =    MOMG(XALPHAG,XNUG,2.)*MOMG(XALPHAR,XNUR,XBR)
XLBRDRYG2   = 2.*MOMG(XALPHAG,XNUG,1.)*MOMG(XALPHAR,XNUR,XBR+1.)
XLBRDRYG3   =                          MOMG(XALPHAR,XNUR,XBR+2.)
XLBNRDRYG1  =    MOMG(XALPHAG,XNUG,2.)                      
XLBNRDRYG2  = 2.*MOMG(XALPHAG,XNUG,1.)*MOMG(XALPHAR,XNUR,1.) 
XLBNRDRYG3  =                          MOMG(XALPHAR,XNUR,2.)
!
! Notice: One magnitude of lambda discretized over 10 points
!
PARAM_LIMA_MIXED%NDRYLBDAR = 40
PARAM_LIMA_MIXED%XDRYLBDAR_MIN = 1.0E3 ! Minimal value of Lbda_r to tabulate XKER_RDRYG
PARAM_LIMA_MIXED%XDRYLBDAR_MAX = 1.0E7 ! Maximal value of Lbda_r to tabulate XKER_RDRYG
ZRATE = LOG(PARAM_LIMA_MIXED%XDRYLBDAR_MAX/PARAM_LIMA_MIXED%XDRYLBDAR_MIN)/REAL(PARAM_LIMA_MIXED%NDRYLBDAR-1)
XDRYINTP1R = 1.0 / ZRATE
XDRYINTP2R = 1.0 - LOG( XDRYLBDAR_MIN ) / ZRATE
PARAM_LIMA_MIXED%NDRYLBDAS = 80
PARAM_LIMA_MIXED%XDRYLBDAS_MIN = 5.0E1*XTRANS_MP_GAMMAS ! Minimal value of Lbda_s to tabulate XKER_SDRYG
PARAM_LIMA_MIXED%XDRYLBDAS_MAX = 5.0E8*XTRANS_MP_GAMMAS ! Maximal value of Lbda_s to tabulate XKER_SDRYG
ZRATE = LOG(PARAM_LIMA_MIXED%XDRYLBDAS_MAX/PARAM_LIMA_MIXED%XDRYLBDAS_MIN)/REAL(PARAM_LIMA_MIXED%NDRYLBDAS-1)
XDRYINTP1S = 1.0 / ZRATE
XDRYINTP2S = 1.0 - LOG( XDRYLBDAS_MIN ) / ZRATE
PARAM_LIMA_MIXED%NDRYLBDAG = 40
PARAM_LIMA_MIXED%XDRYLBDAG_MIN = 1.0E3 ! Min value of Lbda_g to tabulate XKER_SDRYG,XKER_RDRYG
PARAM_LIMA_MIXED%XDRYLBDAG_MAX = 1.0E7 ! Max value of Lbda_g to tabulate XKER_SDRYG,XKER_RDRYG
ZRATE = LOG(PARAM_LIMA_MIXED%XDRYLBDAG_MAX/PARAM_LIMA_MIXED%XDRYLBDAG_MIN)/REAL(PARAM_LIMA_MIXED%NDRYLBDAG-1)
XDRYINTP1G = 1.0 / ZRATE
XDRYINTP2G = 1.0 - LOG( PARAM_LIMA_MIXED%XDRYLBDAG_MIN ) / ZRATE
!
!*       8.2.5  Computations of the tabulated normalized kernels
!
IND      = 50    ! Interval number, collection efficiency and infinite diameter
ZZEGS     = 1.0   ! factor used to integrate the dimensional distributions when
ZZFDINFTY = 20.0  ! computing the kernels XKER_SDRYG
!
CALL PARAM_LIMA_MIXED_ALLOCATE('XKER_SDRYG', NDRYLBDAG,NDRYLBDAS)
!if (NMOM_S.GE.2) then
  CALL PARAM_LIMA_MIXED_ALLOCATE('XKER_N_SDRYG', NDRYLBDAG,NDRYLBDAS)
IF (NMOM_S.GE.2 ) THEN
  CALL NZCOLX  ( NDRYLBDAG,NDRYLBDAS, IND, XALPHAG, XNUG, XALPHAS, XNUS,    & 
              ZZEGS, XCG, XDG, 0., XCS, XDS, XFVELOS,                        & 
              XDRYLBDAG_MAX, XDRYLBDAS_MAX, XDRYLBDAG_MIN, XDRYLBDAS_MIN, & 
              ZZFDINFTY, XKER_N_SDRYG                                      ) 
ELSE
   XKER_N_SDRYG(:,:)=0.
END IF
!!$  WRITE(UNIT=ILUOUT0,FMT='("!")')
!!$  WRITE(UNIT=ILUOUT0,FMT='("IF( PRESENT(PKER_N_SDRYG) ) THEN")')          
!!$  DO I1 = 1 , NDRYLBDAG                                                    
!!$    DO J2 = 1 , NDRYLBDAS                                                  
!!$    WRITE(UNIT=ILUOUT0,FMT='("  PKER_N_SDRYG(",I3,",",I3,") = ",E13.6)') &
!!$                      I1,J2,XKER_N_SDRYG(I1,J2)                          
!!$    END DO                                                                
!!$  END DO                                                                   
!!$  WRITE(UNIT=ILUOUT0,FMT='("!")') 
!end if
!
CALL LIMA_READ_XKER_SDRYG (IDRYLBDAG,IDRYLBDAS,IKND,                                    &
                           ZALPHAG,ZNUG,ZALPHAS,ZNUS,ZEGS,ZBS,ZCG,ZDG,ZCS,ZDS,ZFVELOS, &
                           ZDRYLBDAG_MAX,ZDRYLBDAS_MAX,ZDRYLBDAG_MIN,ZDRYLBDAS_MIN,    &
                           ZFDINFTY                                                    )
IF( (IDRYLBDAG/=NDRYLBDAG) .OR. (IDRYLBDAS/=NDRYLBDAS) .OR. (IKND/=IND) .OR. &
    (ZALPHAG/=XALPHAG) .OR. (ZNUG/=XNUG)                               .OR. &
    (ZALPHAS/=XALPHAS) .OR. (ZNUS/=XNUS)                               .OR. &
    (ZEGS/=ZZEGS) .OR. (ZBS/=XBS)                                       .OR. &
    (ZCG/=XCG) .OR. (ZDG/=XDG) .OR. (ZCS/=XCS) .OR. (ZDS/=XDS) .OR. (ZFVELOS/=XFVELOS) .OR. &
    (ZDRYLBDAG_MAX/=XDRYLBDAG_MAX) .OR. (ABS(ZDRYLBDAS_MAX-XDRYLBDAS_MAX).GE.(0.001*ZDRYLBDAS_MAX)) .OR. &
    (ZDRYLBDAG_MIN/=XDRYLBDAG_MIN) .OR. (ABS(ZDRYLBDAS_MIN-XDRYLBDAS_MIN).GE.(0.001*ZDRYLBDAS_MIN)) .OR. &
    (ZFDINFTY/=ZZFDINFTY)                                               ) THEN
  CALL LIMA_RZCOLX ( IDRYLBDAG, IDRYLBDAS, IND, XALPHAG, XNUG, XALPHAS, XNUS, &
                ZZEGS, XBS, XCG, XDG, 0., XCS, XDS, XFVELOS,                 &
                XDRYLBDAG_MAX, XDRYLBDAS_MAX, XDRYLBDAG_MIN, XDRYLBDAS_MIN, &
                ZZFDINFTY, XKER_SDRYG                                        )
!!$  WRITE(UNIT=ILUOUT0,FMT='("*****************************************")')
!!$  WRITE(UNIT=ILUOUT0,FMT='("**** UPDATE NEW SET OF SDRYG KERNELS ****")')
!!$  WRITE(UNIT=ILUOUT0,FMT='("*****************************************")')
!!$  WRITE(UNIT=ILUOUT0,FMT='("!")')
!!$  WRITE(UNIT=ILUOUT0,FMT='("IKND=",I3)') IND
!!$  WRITE(UNIT=ILUOUT0,FMT='("IDRYLBDAG=",I3)') NDRYLBDAG
!!$  WRITE(UNIT=ILUOUT0,FMT='("IDRYLBDAS=",I3)') NDRYLBDAS
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZALPHAG=",E13.6)') XALPHAG
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZNUG=",E13.6)') XNUG
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZALPHAS=",E13.6)') XALPHAS
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZNUS=",E13.6)') XNUS
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZEGS=",E13.6)') ZZEGS
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZBS=",E13.6)') XBS
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZCG=",E13.6)') XCG
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZDG=",E13.6)') XDG
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZCS=",E13.6)') XCS
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZDS=",E13.6)') XDS
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZFVELOS=",E13.6)') XFVELOS
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZDRYLBDAG_MAX=",E13.6)') &
!!$                                                    XDRYLBDAG_MAX
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZDRYLBDAS_MAX=",E13.6)') &
!!$                                                    XDRYLBDAS_MAX
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZDRYLBDAG_MIN=",E13.6)') &
!!$                                                    XDRYLBDAG_MIN
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZDRYLBDAS_MIN=",E13.6)') &
!!$                                                    XDRYLBDAS_MIN
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZFDINFTY=",E13.6)') ZZFDINFTY
!!$  WRITE(UNIT=ILUOUT0,FMT='("!")')
!!$  WRITE(UNIT=ILUOUT0,FMT='("IF( PRESENT(PKER_SDRYG) ) THEN")')
!!$  DO I1 = 1 , NDRYLBDAG
!!$    DO J2 = 1 , NDRYLBDAS
!!$    WRITE(UNIT=ILUOUT0,FMT='("PKER_SDRYG(",I3,",",I3,") = ",E13.6)') &
!!$                        I1,J2,XKER_SDRYG(I1,J2)
!!$    END DO
!!$  END DO
!!$  WRITE(UNIT=ILUOUT0,FMT='("END IF")')
  ELSE
  CALL LIMA_READ_XKER_SDRYG (IDRYLBDAG,IDRYLBDAS,IKND,                                    &
                             ZALPHAG,ZNUG,ZALPHAS,ZNUS,ZEGS,ZBS,ZCG,ZDG,ZCS,ZDS,ZFVELOS, &
                             ZDRYLBDAG_MAX,ZDRYLBDAS_MAX,ZDRYLBDAG_MIN,ZDRYLBDAS_MIN,    &
                             ZFDINFTY,XKER_SDRYG                                         )
!!$  WRITE(UNIT=ILUOUT0,FMT='(" Read XKER_SDRYG")')
END IF
!
!
IND      = 50    ! Number of interval used to integrate the dimensional
ZZEGR     = 1.0   ! distributions when computing the kernel XKER_RDRYG
ZZFDINFTY = 20.0
!
CALL PARAM_LIMA_MIXED_ALLOCATE('XKER_RDRYG', NDRYLBDAG,NDRYLBDAR)
!if ( NMOM_R.GE.2 ) then
  CALL PARAM_LIMA_MIXED_ALLOCATE('XKER_N_RDRYG', NDRYLBDAG,NDRYLBDAR)
IF (NMOM_R.GE.2) THEN
  CALL NZCOLX  ( NDRYLBDAG,NDRYLBDAR, IND, XALPHAS, XNUS, XALPHAR, XNUR,    & 
               ZZEGR, XCG, XDG, 0., XCR, XDR, 0.,                            & 
               XDRYLBDAG_MAX, XDRYLBDAR_MAX, XDRYLBDAG_MIN, XDRYLBDAR_MIN, &
               ZZFDINFTY, XKER_N_RDRYG                                      )
ELSE
   XKER_N_RDRYG(:,:)=0.
END IF
!!$  WRITE(UNIT=ILUOUT0,FMT='("!")')
!!$  WRITE(UNIT=ILUOUT0,FMT='("IF( PRESENT(PKER_N_RDRYG) ) THEN")')          
!!$  DO I1 = 1 , NDRYLBDAG                                                    
!!$    DO J2 = 1 , NDRYLBDAR                                                  
!!$    WRITE(UNIT=ILUOUT0,FMT='("  PKER_N_RDRYG(",I3,",",I3,") = ",E13.6)') &
!!$                      I1,J2,XKER_N_RDRYG(I1,J2)                          
!!$    END DO                                                                
!!$  END DO                                                                   
!!$  WRITE(UNIT=ILUOUT0,FMT='("!")') 
!end if
!
CALL LIMA_READ_XKER_RDRYG (IDRYLBDAG,IDRYLBDAR,IKND,                                 &
                           ZALPHAG,ZNUG,ZALPHAR,ZNUR,ZEGR,ZBR,ZCG,ZDG,ZCR,ZDR,      &
                           ZDRYLBDAG_MAX,ZDRYLBDAR_MAX,ZDRYLBDAG_MIN,ZDRYLBDAR_MIN, &
                           ZFDINFTY                                                 )
IF( (IDRYLBDAG/=NDRYLBDAG) .OR. (IDRYLBDAR/=NDRYLBDAR) .OR. (IKND/=IND) .OR. &
    (ZALPHAG/=XALPHAG) .OR. (ZNUG/=XNUG)                               .OR. &
    (ZALPHAR/=XALPHAR) .OR. (ZNUR/=XNUR)                               .OR. &
    (ZEGR/=ZZEGR) .OR. (ZBR/=XBR)                                       .OR. &
    (ZCG/=XCG) .OR. (ZDG/=XDG) .OR. (ZCR/=XCR) .OR. (ZDR/=XDR)         .OR. &
    (ZDRYLBDAG_MAX/=XDRYLBDAG_MAX) .OR. (ZDRYLBDAR_MAX/=XDRYLBDAR_MAX) .OR. &
    (ZDRYLBDAG_MIN/=XDRYLBDAG_MIN) .OR. (ZDRYLBDAR_MIN/=XDRYLBDAR_MIN) .OR. &
    (ZFDINFTY/=ZZFDINFTY)                                               ) THEN
  CALL LIMA_RZCOLX ( IDRYLBDAG, IDRYLBDAR, IND, XALPHAG, XNUG, XALPHAR, XNUR, &
                ZZEGR, XBR, XCG, XDG, 0., XCR, XDR, 0.,                      &
                XDRYLBDAG_MAX, XDRYLBDAR_MAX, XDRYLBDAG_MIN, XDRYLBDAR_MIN, &
                ZZFDINFTY, XKER_RDRYG                                        )
!!$  WRITE(UNIT=ILUOUT0,FMT='("*****************************************")')
!!$  WRITE(UNIT=ILUOUT0,FMT='("**** UPDATE NEW SET OF RDRYG KERNELS ****")')
!!$  WRITE(UNIT=ILUOUT0,FMT='("*****************************************")')
!!$  WRITE(UNIT=ILUOUT0,FMT='("!")')
!!$  WRITE(UNIT=ILUOUT0,FMT='("IKND=",I3)') IND
!!$  WRITE(UNIT=ILUOUT0,FMT='("IDRYLBDAG=",I3)') NDRYLBDAG
!!$  WRITE(UNIT=ILUOUT0,FMT='("IDRYLBDAR=",I3)') NDRYLBDAR
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZALPHAG=",E13.6)') XALPHAG
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZNUG=",E13.6)') XNUG
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZALPHAR=",E13.6)') XALPHAR
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZNUR=",E13.6)') XNUR
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZEGR=",E13.6)') ZZEGR
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZBR=",E13.6)') XBR
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZCG=",E13.6)') XCG
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZDG=",E13.6)') XDG
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZCR=",E13.6)') XCR
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZDR=",E13.6)') XDR
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZDRYLBDAG_MAX=",E13.6)') &
!!$                                                    XDRYLBDAG_MAX
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZDRYLBDAR_MAX=",E13.6)') &
!!$                                                    XDRYLBDAR_MAX
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZDRYLBDAG_MIN=",E13.6)') &
!!$                                                    XDRYLBDAG_MIN
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZDRYLBDAR_MIN=",E13.6)') &
!!$                                                    XDRYLBDAR_MIN
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZFDINFTY=",E13.6)') ZZFDINFTY
!!$  WRITE(UNIT=ILUOUT0,FMT='("!")')
!!$  WRITE(UNIT=ILUOUT0,FMT='("IF( PRESENT(PKER_RDRYG) ) THEN")')
!!$  DO I1 = 1 , NDRYLBDAG
!!$    DO J2 = 1 , NDRYLBDAR
!!$    WRITE(UNIT=ILUOUT0,FMT='("PKER_RDRYG(",I3,",",I3,") = ",E13.6)') &
!!$                        I1,J2,XKER_RDRYG(I1,J2)
!!$    END DO
!!$  END DO
!!$  WRITE(UNIT=ILUOUT0,FMT='("END IF")')
  ELSE
  CALL LIMA_READ_XKER_RDRYG (IDRYLBDAG,IDRYLBDAR,IKND,                                 &
                             ZALPHAG,ZNUG,ZALPHAR,ZNUR,ZEGR,ZBR,ZCG,ZDG,ZCR,ZDR,      &
                             ZDRYLBDAG_MAX,ZDRYLBDAR_MAX,ZDRYLBDAG_MIN,ZDRYLBDAR_MIN, &
                             ZFDINFTY,XKER_RDRYG                                      )
!!$  WRITE(UNIT=ILUOUT0,FMT='(" Read XKER_RDRYG")')
END IF
!  
!-------------------------------------------------------------------------------
!
!*       9.     CONSTANTS FOR THE FAST COLD PROCESSES FOR THE HAILSTONES
!               --------------------------------------------------------
!
!*       9.2    Constants for the wet growth of the hailstones
!
!
!*       9.2.1  Constant for the cloud droplet and cloud ice collection
!               by the hailstones
!
XFWETH = (XPI/4.0)*XCH*(ZRHO00**XCEXVT)*MOMG(XALPHAH,XNUH,XDH+2.0)
!
!*       9.2.2  Constants for the aggregate collection by the hailstones
!
XFSWETH = XNS*(XPI/4.0)*XAS*(ZRHO00**XCEXVT)
XFNSWETH= (XPI/4.0)*(ZRHO00**XCEXVT)
!
XLBSWETH1   =    MOMG(XALPHAH,XNUH,2.)*MOMG(XALPHAS,XNUS,XBS)
XLBSWETH2   = 2.*MOMG(XALPHAH,XNUH,1.)*MOMG(XALPHAS,XNUS,XBS+1.)
XLBSWETH3   =                          MOMG(XALPHAS,XNUS,XBS+2.)
XLBNSWETH1  =    MOMG(XALPHAH,XNUH,2.)                       
XLBNSWETH2  = 2.*MOMG(XALPHAH,XNUH,1.)*MOMG(XALPHAS,XNUS,1.) 
XLBNSWETH3  =                          MOMG(XALPHAS,XNUS,2.)
!
!*       9.2.3  Constants for the graupel collection by the hailstones
!
XFGWETH = XNG*(XPI/4.0)*XAG*(ZRHO00**XCEXVT)
XFNGWETH= (XPI/4.0)*(ZRHO00**XCEXVT)
!
XLBGWETH1   =    MOMG(XALPHAH,XNUH,2.)*MOMG(XALPHAG,XNUG,XBG)
XLBGWETH2   = 2.*MOMG(XALPHAH,XNUH,1.)*MOMG(XALPHAG,XNUG,XBG+1.)
XLBGWETH3   =                          MOMG(XALPHAG,XNUG,XBG+2.)
XLBNGWETH1  =    MOMG(XALPHAH,XNUH,2.)                      
XLBNGWETH2  = 2.*MOMG(XALPHAH,XNUH,1.)*MOMG(XALPHAG,XNUG,1.)
XLBNGWETH3  =                          MOMG(XALPHAG,XNUG,2.)
!
! Notice: One magnitude of lambda discretized over 10 points
!
PARAM_LIMA_MIXED%NWETLBDAS = 80
PARAM_LIMA_MIXED%XWETLBDAS_MIN = 5.0E1*XTRANS_MP_GAMMAS ! Minimal value of Lbda_s to tabulate XKER_SWETH
PARAM_LIMA_MIXED%XWETLBDAS_MAX = 5.0E8*XTRANS_MP_GAMMAS ! Maximal value of Lbda_s to tabulate XKER_SWETH
ZRATE = LOG(PARAM_LIMA_MIXED%XWETLBDAS_MAX/PARAM_LIMA_MIXED%XWETLBDAS_MIN)/REAL(PARAM_LIMA_MIXED%NWETLBDAS-1)
XWETINTP1S = 1.0 / ZRATE
XWETINTP2S = 1.0 - LOG( PARAM_LIMA_MIXED%XWETLBDAS_MIN ) / ZRATE
PARAM_LIMA_MIXED%NWETLBDAG = 40
PARAM_LIMA_MIXED%XWETLBDAG_MIN = 1.0E3 ! Min value of Lbda_g to tabulate XKER_GWETH
PARAM_LIMA_MIXED%XWETLBDAG_MAX = 1.0E7 ! Max value of Lbda_g to tabulate XKER_GWETH
ZRATE = LOG(PARAM_LIMA_MIXED%XWETLBDAG_MAX/PARAM_LIMA_MIXED%XWETLBDAG_MIN)/REAL(PARAM_LIMA_MIXED%NWETLBDAG-1)
XWETINTP1G = 1.0 / ZRATE
XWETINTP2G = 1.0 - LOG( XWETLBDAG_MIN ) / ZRATE
PARAM_LIMA_MIXED%NWETLBDAH = 40
PARAM_LIMA_MIXED%XWETLBDAH_MIN = 1.0E3 ! Min value of Lbda_h to tabulate XKER_SWETH,XKER_GWETH
PARAM_LIMA_MIXED%XWETLBDAH_MAX = 1.0E7 ! Max value of Lbda_h to tabulate XKER_SWETH,XKER_GWETH
ZRATE = LOG(PARAM_LIMA_MIXED%XWETLBDAH_MAX/PARAM_LIMA_MIXED%XWETLBDAH_MIN)/REAL(PARAM_LIMA_MIXED%NWETLBDAH-1)
XWETINTP1H = 1.0 / ZRATE
XWETINTP2H = 1.0 - LOG( PARAM_LIMA_MIXED%XWETLBDAH_MIN ) / ZRATE
!
!*       9.2.4  Computations of the tabulated normalized kernels
!
IND      = 50    ! Interval number, collection efficiency and infinite diameter
ZZEHS     = 1.0   ! factor used to integrate the dimensional distributions when
ZZFDINFTY = 20.0  ! computing the kernels XKER_SWETH
!
!if ( NMOM_S.GE.2 ) then
  IF( .NOT.ASSOCIATED(XKER_N_SWETH) ) CALL PARAM_LIMA_MIXED_ALLOCATE('XKER_N_SWETH', NWETLBDAH,NWETLBDAS)
IF (NMOM_S.GE.2) THEN
  CALL NZCOLX ( NWETLBDAH,NWETLBDAS,IND, XALPHAH, XNUH, XALPHAS, XNUS,     &  
              ZZEHS, XCH, XDH, 0., XCS, XDS, XFVELOS,                       &  ! 
              XWETLBDAH_MAX, XWETLBDAS_MAX, XWETLBDAH_MIN, XWETLBDAS_MIN, &  ! 
              ZZFDINFTY, XKER_N_SWETH                                        )
ELSE
   XKER_N_SWETH(:,:)=0.
END IF
!!$  WRITE(UNIT=ILUOUT0,FMT='("!")')
!!$  WRITE(UNIT=ILUOUT0,FMT='("IF( PRESENT(PKER_N_SWETH) ) THEN")')          
!!$  DO I1 = 1 , NWETLBDAH                                                    
!!$    DO J2 = 1 , NWETLBDAS                                                  
!!$    WRITE(UNIT=ILUOUT0,FMT='("  PKER_N_SWETH(",I3,",",I3,") = ",E13.6)') &
!!$                      I1,J2,XKER_N_SWETH(I1,J2)                          
!!$    END DO                                                                
!!$  END DO                                                                   
!!$  WRITE(UNIT=ILUOUT0,FMT='("!")') 
!end  if
IF( .NOT.ASSOCIATED(XKER_SWETH) ) CALL PARAM_LIMA_MIXED_ALLOCATE('XKER_SWETH', NWETLBDAH,NWETLBDAS)
!
CALL LIMA_READ_XKER_SWETH (IWETLBDAH,IWETLBDAS,IKND,                                    &
                           ZALPHAH,ZNUH,ZALPHAS,ZNUS,ZEHS,ZBS,ZCH,ZDH,ZCS,ZDS,ZFVELOS, &
                           ZWETLBDAH_MAX,ZWETLBDAS_MAX,ZWETLBDAH_MIN,ZWETLBDAS_MIN,    &
                           ZFDINFTY                                                    )
IF( (IWETLBDAH/=NWETLBDAH) .OR. (IWETLBDAS/=NWETLBDAS) .OR. (IKND/=IND) .OR. &
    (ZALPHAH/=XALPHAH) .OR. (ZNUH/=XNUH)                               .OR. &
    (ZALPHAS/=XALPHAS) .OR. (ZNUS/=XNUS)                               .OR. &
    (ZEHS/=ZZEHS) .OR. (ZBS/=XBS)                                       .OR. &
    (ZCH/=XCH) .OR. (ZDH/=XDH) .OR. (ZCS/=XCS) .OR. (ZDS/=XDS) .OR. (ZFVELOS/=XFVELOS) .OR. &
    (ZWETLBDAH_MAX/=XWETLBDAH_MAX) .OR. (ABS(ZWETLBDAS_MAX-XWETLBDAS_MAX).GE.(0.001*ZWETLBDAS_MAX)) .OR. &
    (ZWETLBDAH_MIN/=XWETLBDAH_MIN) .OR. (ABS(ZWETLBDAS_MIN-XWETLBDAS_MIN).GE.(0.001*ZWETLBDAS_MIN)) .OR. &
    (ZFDINFTY/=ZZFDINFTY)                                               ) THEN
  CALL LIMA_RZCOLX ( IWETLBDAH, IWETLBDAS, IND, XALPHAH, XNUH, XALPHAS, XNUS, &
                ZZEHS, XBS, XCH, XDH, 0., XCS, XDS, XFVELOS,                 &
                XWETLBDAH_MAX, XWETLBDAS_MAX, XWETLBDAH_MIN, XWETLBDAS_MIN, &
                ZZFDINFTY, XKER_SWETH                                        )
!!$  WRITE(UNIT=ILUOUT0,FMT='("*****************************************")')
!!$  WRITE(UNIT=ILUOUT0,FMT='("**** UPDATE NEW SET OF SWETH KERNELS ****")')
!!$  WRITE(UNIT=ILUOUT0,FMT='("*****************************************")')
!!$  WRITE(UNIT=ILUOUT0,FMT='("!")')
!!$  WRITE(UNIT=ILUOUT0,FMT='("IKND=",I3)') IND
!!$  WRITE(UNIT=ILUOUT0,FMT='("IWETLBDAH=",I3)') NWETLBDAH
!!$  WRITE(UNIT=ILUOUT0,FMT='("IWETLBDAS=",I3)') NWETLBDAS
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZALPHAH=",E13.6)') XALPHAH
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZNUH=",E13.6)') XNUH
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZALPHAS=",E13.6)') XALPHAS
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZNUS=",E13.6)') XNUS
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZEHS=",E13.6)') ZZEHS
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZBS=",E13.6)') XBS
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZCH=",E13.6)') XCH
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZDH=",E13.6)') XDH
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZCS=",E13.6)') XCS
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZDS=",E13.6)') XDS
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZFVELOS=",E13.6)') XFVELOS
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZWETLBDAH_MAX=",E13.6)') &
!!$                                                    XWETLBDAH_MAX
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZWETLBDAS_MAX=",E13.6)') &
!!$                                                    XWETLBDAS_MAX
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZWETLBDAH_MIN=",E13.6)') &
!!$                                                    XWETLBDAH_MIN
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZWETLBDAS_MIN=",E13.6)') &
!!$                                                    XWETLBDAS_MIN
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZFDINFTY=",E13.6)') ZZFDINFTY
!!$  WRITE(UNIT=ILUOUT0,FMT='("!")')
!!$  WRITE(UNIT=ILUOUT0,FMT='("IF( PRESENT(PKER_SWETH) ) THEN")')
!!$  DO I1 = 1 , NWETLBDAH
!!$    DO J2 = 1 , NWETLBDAS
!!$    WRITE(UNIT=ILUOUT0,FMT='("PKER_SWETH(",I3,",",I3,") = ",E13.6)') &
!!$                        I1,J2,XKER_SWETH(I1,J2)
!!$    END DO
!!$  END DO
!!$  WRITE(UNIT=ILUOUT0,FMT='("END IF")')
  ELSE
  CALL LIMA_READ_XKER_SWETH (IWETLBDAH,IWETLBDAS,IKND,                                    &
                             ZALPHAH,ZNUH,ZALPHAS,ZNUS,ZEHS,ZBS,ZCH,ZDH,ZCS,ZDS,ZFVELOS, &
                             ZWETLBDAH_MAX,ZWETLBDAS_MAX,ZWETLBDAH_MIN,ZWETLBDAS_MIN,    &
                             ZFDINFTY,XKER_SWETH                                         )
!!$  WRITE(UNIT=ILUOUT0,FMT='(" Read XKER_SWETH")')
END IF
!
!
IND      = 50    ! Number of interval used to integrate the dimensional
ZZEHG     = 1.0   ! distributions when computing the kernel XKER_GWETH
ZZFDINFTY = 20.0
!
!if ( NMOM_G.GE.2 ) then
  IF( .NOT.ASSOCIATED(XKER_N_GWETH) ) CALL PARAM_LIMA_MIXED_ALLOCATE('XKER_N_GWETH', NWETLBDAH,NWETLBDAG)
IF (NMOM_G.GE.2) THEN
  CALL NZCOLX ( NWETLBDAH,NWETLBDAG,IND, XALPHAH, XNUH, XALPHAG, XNUG,   & 
              ZZEHG, XCH, XDH, 0., XCG, XDG, 0.,                            & 
              XWETLBDAH_MAX, XWETLBDAG_MAX, XWETLBDAH_MIN, XWETLBDAG_MIN, & 
              ZZFDINFTY, XKER_N_GWETH                                      )
ELSE
   XKER_N_GWETH(:,:)=0.
END IF
!!$  WRITE(UNIT=ILUOUT0,FMT='("!")')
!!$  WRITE(UNIT=ILUOUT0,FMT='("IF( PRESENT(PKER_N_GWETH) ) THEN")')          
!!$  DO I1 = 1 , NWETLBDAH                                                    
!!$    DO J2 = 1 , NWETLBDAG                                                  
!!$    WRITE(UNIT=ILUOUT0,FMT='("  PKER_N_GWETH(",I3,",",I3,") = ",E13.6)') &
!!$                      I1,J2,XKER_N_GWETH(I1,J2)                          
!!$    END DO                                                                
!!$  END DO                                                                   
!!$  WRITE(UNIT=ILUOUT0,FMT='("!")') 
!end if
IF( .NOT.ASSOCIATED(XKER_GWETH) ) CALL PARAM_LIMA_MIXED_ALLOCATE('XKER_GWETH', NWETLBDAH,NWETLBDAG)
!
CALL LIMA_READ_XKER_GWETH (IWETLBDAH,IWETLBDAG,IKND,                                 &
                           ZALPHAH,ZNUH,ZALPHAG,ZNUG,ZEHG,ZBG,ZCH,ZDH,ZCG,ZDG,      &
                           ZWETLBDAH_MAX,ZWETLBDAG_MAX,ZWETLBDAH_MIN,ZWETLBDAG_MIN, &
                           ZFDINFTY                                                 )
IF( (IWETLBDAH/=NWETLBDAH) .OR. (IWETLBDAG/=NWETLBDAG) .OR. (IKND/=IND) .OR. &
    (ZALPHAH/=XALPHAH) .OR. (ZNUH/=XNUH)                               .OR. &
    (ZALPHAG/=XALPHAG) .OR. (ZNUG/=XNUG)                               .OR. &
    (ZEHG/=ZZEHG) .OR. (ZBG/=XBG)                                       .OR. &
    (ZCH/=XCH) .OR. (ZDH/=XDH) .OR. (ZCG/=XCG) .OR. (ZDG/=XDG)         .OR. &
    (ZWETLBDAH_MAX/=XWETLBDAH_MAX) .OR. (ZWETLBDAG_MAX/=XWETLBDAG_MAX) .OR. &
    (ZWETLBDAH_MIN/=XWETLBDAH_MIN) .OR. (ZWETLBDAG_MIN/=XWETLBDAG_MIN) .OR. &
    (ZFDINFTY/=ZZFDINFTY)                                               ) THEN
  CALL LIMA_RZCOLX ( IWETLBDAH, IWETLBDAG, IND, XALPHAH, XNUH, XALPHAG, XNUG, &
                ZZEHG, XBG, XCH, XDH, 0., XCG, XDG, 0.,                      &
                XWETLBDAH_MAX, XWETLBDAG_MAX, XWETLBDAH_MIN, XWETLBDAG_MIN, &
                ZZFDINFTY, XKER_GWETH                                        )
!!$  WRITE(UNIT=ILUOUT0,FMT='("*****************************************")')
!!$  WRITE(UNIT=ILUOUT0,FMT='("**** UPDATE NEW SET OF GWETH KERNELS ****")')
!!$  WRITE(UNIT=ILUOUT0,FMT='("*****************************************")')
!!$  WRITE(UNIT=ILUOUT0,FMT='("!")')
!!$  WRITE(UNIT=ILUOUT0,FMT='("IKND=",I3)') IND
!!$  WRITE(UNIT=ILUOUT0,FMT='("IWETLBDAH=",I3)') NWETLBDAH
!!$  WRITE(UNIT=ILUOUT0,FMT='("IWETLBDAG=",I3)') NWETLBDAG
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZALPHAH=",E13.6)') XALPHAH
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZNUH=",E13.6)') XNUH
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZALPHAG=",E13.6)') XALPHAG
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZNUG=",E13.6)') XNUG
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZEHG=",E13.6)') ZZEHG
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZBG=",E13.6)') XBG
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZCH=",E13.6)') XCH
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZDH=",E13.6)') XDH
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZCG=",E13.6)') XCG
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZDG=",E13.6)') XDG
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZWETLBDAH_MAX=",E13.6)') &
!!$                                                    XWETLBDAH_MAX
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZWETLBDAG_MAX=",E13.6)') &
!!$                                                    XWETLBDAG_MAX
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZWETLBDAH_MIN=",E13.6)') &
!!$                                                    XWETLBDAH_MIN
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZWETLBDAG_MIN=",E13.6)') &
!!$                                                    XWETLBDAG_MIN
!!$  WRITE(UNIT=ILUOUT0,FMT='("ZFDINFTY=",E13.6)') ZZFDINFTY
!!$  WRITE(UNIT=ILUOUT0,FMT='("!")')
!!$  WRITE(UNIT=ILUOUT0,FMT='("IF( PRESENT(PKER_GWETH) ) THEN")')
!!$  DO I1 = 1 , NWETLBDAH
!!$    DO J2 = 1 , NWETLBDAG
!!$    WRITE(UNIT=ILUOUT0,FMT='("PKER_GWETH(",I3,",",I3,") = ",E13.6)') &
!!$                        I1,J2,XKER_GWETH(I1,J2)
!!$    END DO
!!$  END DO
!!$  WRITE(UNIT=ILUOUT0,FMT='("END IF")')
  ELSE
  CALL LIMA_READ_XKER_GWETH (IWETLBDAH,IWETLBDAG,IKND,                                 &
                             ZALPHAH,ZNUH,ZALPHAG,ZNUG,ZEHG,ZBG,ZCH,ZDH,ZCG,ZDG,      &
                             ZWETLBDAH_MAX,ZWETLBDAG_MAX,ZWETLBDAH_MIN,ZWETLBDAG_MIN, &
                             ZFDINFTY,XKER_GWETH                                      )
!!$  WRITE(UNIT=ILUOUT0,FMT='(" Read XKER_GWETH")')
END IF
!
!
!
!-------------------------------------------------------------------------------
!
!*      10.     SET-UP RADIATIVE PARAMETERS
!               ---------------------------
!
!
! R_eff_i = XFREFFI * (rho*r_i/N_i)**(1/3)
!
XFREFFI = 0.5 * ZGAMI(8) * (1.0/XLBI)**XLBEXI
!  
!-------------------------------------------------------------------------------
IF (LCRYSTAL_SHAPE) THEN
  CALL PARAM_LIMA_COLD_ALLOCATE_1D('XFREFFI_SHAPE', NNB_CRYSTAL_SHAPE)
  DO ISH = 1, NNB_CRYSTAL_SHAPE
    XFREFFI_SHAPE(ISH) = 0.5 * ZGAMI(8) * (1.0/XLBI_SHAPE(ISH))**XLBEXI_SHAPE(ISH)
  END DO
END IF
!
!
!*      11.     SOME PRINTS FOR CONTROL
!               -----------------------
!
!
!!$GFLAG = .TRUE.
!!$IF (GFLAG) THEN
!!$  WRITE(UNIT=ILUOUT0,FMT='(" Summary of the ice particule characteristics")')
!!$  WRITE(UNIT=ILUOUT0,FMT='("      PRISTINE ICE")')
!!$  WRITE(UNIT=ILUOUT0,FMT='("                   masse: A=",E13.6," B=",E13.6)') &
!!$                                                      XAI,XBI
!!$  WRITE(UNIT=ILUOUT0,FMT='("                 vitesse: C=",E13.6," D=",E13.6)') &
!!$                                                      XC_I,XDI
!!$  WRITE(UNIT=ILUOUT0,FMT='("            distribution:AL=",E13.6,"NU=",E13.6)') &
!!$                                                      XALPHAI,XNUI
!!$  WRITE(UNIT=ILUOUT0,FMT='("              SNOW")')
!!$  WRITE(UNIT=ILUOUT0,FMT='("                   masse: A=",E13.6," B=",E13.6)') &
!!$                                                      XAS,XBS
!!$  WRITE(UNIT=ILUOUT0,FMT='("                 vitesse: C=",E13.6," D=",E13.6)') &
!!$                                                      XCS,XDS
!!$  WRITE(UNIT=ILUOUT0,FMT='("           concentration:CC=",E13.6," x=",E13.6)') &
!!$                                                      XCCS,XCXS
!!$  WRITE(UNIT=ILUOUT0,FMT='("            distribution:AL=",E13.6,"NU=",E13.6)') &
!!$                                                      XALPHAS,XNUS
!!$  WRITE(UNIT=ILUOUT0,FMT='("            GRAUPEL")')
!!$  WRITE(UNIT=ILUOUT0,FMT='("                   masse: A=",E13.6," B=",E13.6)') &
!!$                                                      XAG,XBG
!!$  WRITE(UNIT=ILUOUT0,FMT='("                 vitesse: C=",E13.6," D=",E13.6)') &
!!$                                                      XCG,XDG
!!$  WRITE(UNIT=ILUOUT0,FMT='("           concentration:CC=",E13.6," x=",E13.6)') &
!!$                                                      XCCG,XCXG
!!$  WRITE(UNIT=ILUOUT0,FMT='("            distribution:AL=",E13.6,"NU=",E13.6)') &
!!$                                                      XALPHAG,XNUG
!!$END IF
!
!------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('INI_LIMA_COLD_MIXED', 1, ZHOOK_HANDLE)
END SUBROUTINE INI_LIMA_COLD_MIXED
!
END MODULE MODE_INI_LIMA_COLD_MIXED
