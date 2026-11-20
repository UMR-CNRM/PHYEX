!@no_insert_drhook
!     ######spl
      MODULE MODD_TIWMX

!$ACDC methods

!     ###############
!
!!****  *MODD_TIWMX* - 
!!
!!    PURPOSE
!!    -------
!       Lookup tables (saturation pressures and their derivatives, ice crystal and water
!     droplet growth functions...) as a function of temperature, and their initialisation
!     (INI_TIWMX). The functions reading the tables are in tiwmx.func.h.
!
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (ha ha)
!!          
!!    AUTHOR
!!    ------
!!      K. I. Ivarsson   *SMHI*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    20/11/14  
!!      S. Riette   Nov. 2026: INI_TIWMX (formerly in mode_ini_tiwmx.F90) and the analytic formulas
!!                             (formerly in mode_tiwmx_fun.F90) moved here
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE

REAL, PARAMETER :: XNDEGR = 100.0
INTEGER, PARAMETER :: NSTART = 13200 ! A too small value may result into a FPE in single precision mode. REK.
INTEGER, PARAMETER :: NSTOP = 37316

TYPE TIWMX_t
! Saturation tables and derivatives
REAL ::  ESTABW(NSTART:NSTOP)
REAL :: DESTABW(NSTART:NSTOP)
REAL ::  ESTABI(NSTART:NSTOP)
REAL :: DESTABI(NSTART:NSTOP)

! Ice crystal- or water droplet tables
REAL ::   A2TAB(NSTART:NSTOP)
REAL ::  BB3TAB(NSTART:NSTOP)
REAL ::  AM3TAB(NSTART:NSTOP)
REAL ::  AF3TAB(NSTART:NSTOP)
REAL ::  A2WTAB(NSTART:NSTOP)
REAL :: BB3WTAB(NSTART:NSTOP)
REAL :: REDINTAB(NSTART:NSTOP)
END TYPE TIWMX_t

TYPE(TIWMX_t), SAVE, TARGET :: TIWMX
!
CONTAINS
!
SUBROUTINE INI_TIWMX
! Fills the tables with the analytic formulas (internal functions below)
  USE YOMHOOK, ONLY : LHOOK, DR_HOOK, JPHOOK
  IMPLICIT NONE

  INTEGER :: JK
  REAL :: ZTEMP

  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

  IF (LHOOK) CALL DR_HOOK('INI_TIWMX',0,ZHOOK_HANDLE)

  DO JK=NSTART,NSTOP
     ZTEMP = REAL(JK)/XNDEGR

     TIWMX%ESTABW(JK) = ESATW(ZTEMP)
     TIWMX%DESTABW(JK) = DESDTW(ZTEMP)

     TIWMX%ESTABI(JK) = ESATI(ZTEMP)
     TIWMX%DESTABI(JK) = DESDTI(ZTEMP)

!    Functions for ice crystals or water droplets:
     TIWMX%A2TAB(JK) = AA2(ZTEMP)
     TIWMX%A2WTAB(JK) = AA2W(ZTEMP)
     TIWMX%AM3TAB(JK) = AM3(ZTEMP)
     TIWMX%AF3TAB(JK) = AF3(ZTEMP)
     TIWMX%BB3TAB(JK) = BB3(ZTEMP)
     TIWMX%BB3WTAB(JK) = BB3W(ZTEMP)
     TIWMX%REDINTAB(JK) = REDIN(ZTEMP)

     IF (ZTEMP.GT.273.15) THEN
!       Avoid ice calculation above freezing level:
        TIWMX%ESTABI(JK) = ESATW(ZTEMP)
        TIWMX%DESTABI(JK) = DESDTW(ZTEMP)

!       Functions for ice crystal growth:  
        TIWMX%A2TAB(JK) = AA2(ZTEMP)
        TIWMX%A2WTAB(JK) = AA2W(ZTEMP)
        TIWMX%AM3TAB(JK) = 0.
        TIWMX%AF3TAB(JK) = 0.
        TIWMX%BB3TAB(JK) = BB3(ZTEMP)
        TIWMX%BB3WTAB(JK) = BB3W(ZTEMP)
     ENDIF
  END DO

  IF (LHOOK) CALL DR_HOOK('INI_TIWMX',1,ZHOOK_HANDLE)

!
CONTAINS
!
! Analytic formulas, used only to fill the tables (the functions used during the simulation,
! which read the tables, are in tiwmx.func.h)
!
!
  ELEMENTAL REAL FUNCTION ESATW(TT)
    USE MODD_CST, ONLY : XALPW, XBETAW, XGAMW
    REAL,INTENT(IN) :: TT
    ESATW = EXP( XALPW - XBETAW/TT - XGAMW*LOG(TT) )
  END FUNCTION ESATW
!      
!     pure saturation pressure over ice for tt <0 C,
!     esatw otherwise.
!
  ELEMENTAL REAL FUNCTION ESATI(TT)
    USE MODD_CST, ONLY : XALPW, XBETAW, XGAMW, XALPI,XBETAI,XGAMI, XTT
    REAL,INTENT(IN) :: TT
    ESATI = ( 0.5 + SIGN(0.5,TT-XTT) )*EXP( XALPW - XBETAW/TT - XGAMW*LOG(TT) ) - &
         &   ( SIGN(0.5,TT-XTT) - 0.5)*EXP( XALPI - XBETAI/TT - XGAMI*LOG(TT) )
  END FUNCTION ESATI
!
!     pure saturation pressure over water
  ELEMENTAL REAL FUNCTION DESDTW(TT)
    USE MODD_CST, ONLY : XBETAW, XGAMW
    REAL,INTENT(IN) :: TT
    DESDTW = ESATW(TT)*(XBETAW/TT - XGAMW)/TT
  END FUNCTION DESDTW

  ELEMENTAL REAL FUNCTION DESDTI(TT)
    USE MODD_CST, ONLY : XBETAI, XGAMI, XTT
    REAL,INTENT(IN) :: TT
    DESDTI = ( 0.5 + SIGN(0.5,TT-XTT) )*DESDTW(TT) - &
         & ( SIGN(0.5,TT-XTT) - 0.5)*ESATI(TT)*(XBETAI/TT - XGAMI)/TT
  END FUNCTION DESDTI

!     Ice crystal function:
  ELEMENTAL REAL FUNCTION AA2(TT)
    USE MODD_CST, ONLY : XLSTT, XCPV, XCI, XTT, XRV
    REAL,INTENT(IN) :: TT
    AA2 =  ( XLSTT + (XCPV-XCI)*(TT-XTT) )**2 / &
     & (2.38E-2 + 0.0071E-2 *(TT - XTT))/(TT**2*XRV)
  END FUNCTION AA2

!     Water droplet function:
  ELEMENTAL REAL FUNCTION AA2W(TT)
    USE MODD_CST, ONLY : XLVTT, XCPV, XCL, XTT, XRV
    REAL,INTENT(IN) :: TT
    AA2W = ( (XLVTT+ (XCPV-XCL)*(TT-XTT))**2)/ &
     &  (2.38E-2 + 0.0071E-2 *(TT - XTT))/(TT**2*XRV)
  END FUNCTION AA2W

!     Ice crystal function:
  ELEMENTAL REAL FUNCTION BB3(TT)
    USE MODD_CST, ONLY : XTT, XRV, XP00
    REAL,INTENT(IN) :: TT
    BB3 = XRV/(0.211E-4 * (TT/XTT)**1.94 * XP00)*TT/ESATI(TT)
  END FUNCTION BB3

!     Water droplet function:
  ELEMENTAL REAL FUNCTION BB3W(TT)
    USE MODD_CST, ONLY : XTT, XRV, XP00
    REAL,INTENT(IN) :: TT
    BB3W =  XRV/(0.211E-4 * (TT/XTT)**1.94 * XP00)*TT/ESATW(TT)
  END FUNCTION BB3W

! Meyers IN concentration function:
  ELEMENTAL REAL FUNCTION AM3(TT)
    REAL,INTENT(IN) :: TT
    AM3 = 1000.*EXP(12.96*(ESATW(TT)/ESATI(TT) -1.) -0.639)
  END FUNCTION AM3

! Fletchers IN concentration function:
  ELEMENTAL REAL FUNCTION AF3(TT)
    USE MODD_CST, ONLY : XTT
    REAL,INTENT(IN) :: TT
    AF3 = 0.01*EXP(0.6*(XTT-TT))
  END FUNCTION AF3

! Function for IN concentration reduction between 0 and -20 C:
  ELEMENTAL REAL FUNCTION REDIN(TT)
    USE MODD_CST, ONLY : XTT
    REAL,INTENT(IN) :: TT
    REAL ZZT
    ZZT = MAX(0., MIN(1.,(XTT - TT)/20.))
    REDIN = 1.- (1.-ZZT)/(ZZT**3 + (1.-ZZT)**3)**.333    
  END FUNCTION REDIN
END SUBROUTINE INI_TIWMX
!
END MODULE MODD_TIWMX
