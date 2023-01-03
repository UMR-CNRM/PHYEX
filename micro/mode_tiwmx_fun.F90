!@no_insert_drhook
!     ######spl
      MODULE MODE_TIWMX_FUN
!     ###############
!
!!****  *MODD_TIWMX_FUN* - 
!!
!!    PURPOSE
!!    -------
!       The purpose of this  ...
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
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_CST, ONLY : XALPW,XBETAW,XGAMW,XTT,XALPI,XBETAI,XGAMI,XLSTT,XRV,XLVTT,&
     &XLSTT,XP00,XCPV,XCI,XCL

IMPLICIT NONE

CONTAINS
!
  REAL FUNCTION ESATW(TT)
    REAL,INTENT(IN) :: TT
    ESATW = EXP( XALPW - XBETAW/TT - XGAMW*ALOG(TT) )
  END FUNCTION ESATW
!      
!     pure saturation pressure over ice for tt <0 C,
!     esatw otherwise.
!
  REAL FUNCTION ESATI(TT)
    REAL,INTENT(IN) :: TT
    ESATI = ( 0.5 + SIGN(0.5,TT-XTT) )*EXP( XALPW - XBETAW/TT - XGAMW*ALOG(TT) ) - &
         &   ( SIGN(0.5,TT-XTT) - 0.5)*EXP( XALPI - XBETAI/TT - XGAMI*ALOG(TT) )
  END FUNCTION ESATI
!
!     pure saturation pressure over water
  REAL FUNCTION DESDTW(TT)
    REAL,INTENT(IN) :: TT
    DESDTW = ESATW(TT)*(XBETAW/TT - XGAMW)/TT
  END FUNCTION DESDTW

  REAL FUNCTION DESDTI(TT)
    REAL,INTENT(IN) :: TT
    DESDTI = ( 0.5 + SIGN(0.5,TT-XTT) )*DESDTW(TT) - &
         & ( SIGN(0.5,TT-XTT) - 0.5)*ESATI(TT)*(XBETAI/TT - XGAMI)/TT
  END FUNCTION DESDTI

!     Ice crystal function:
  REAL FUNCTION AA2(TT)
    REAL,INTENT(IN) :: TT
    AA2 =  ( XLSTT + (XCPV-XCI)*(TT-XTT) )**2 / &
     & (2.38E-2 + 0.0071E-2 *(TT - XTT))/(TT**2*XRV)
  END FUNCTION AA2

!     Water droplet function:
  REAL FUNCTION AA2W(TT)
    REAL,INTENT(IN) :: TT
    AA2W = ( (XLVTT+ (XCPV-XCL)*(TT-XTT))**2)/ &
     &  (2.38E-2 + 0.0071E-2 *(TT - XTT))/(TT**2*XRV)
  END FUNCTION AA2W

!     Ice crystal function:
  REAL FUNCTION BB3(TT)
    REAL,INTENT(IN) :: TT
    BB3 = XRV/(0.211E-4 * (TT/XTT)**1.94 * XP00)*TT/ESATI(TT)
  END FUNCTION BB3

!     Water droplet function:
  REAL FUNCTION BB3W(TT)
    REAL,INTENT(IN) :: TT
    BB3W =  XRV/(0.211E-4 * (TT/XTT)**1.94 * XP00)*TT/ESATW(TT)
  END FUNCTION BB3W

! Meyers IN concentration function:
  REAL FUNCTION AM3(TT)
    REAL,INTENT(IN) :: TT
    AM3 = 1000.*EXP(12.96*(ESATW(TT)/ESATI(TT) -1.) -0.639)
  END FUNCTION AM3

! Fletchers IN concentration function:
  REAL FUNCTION AF3(TT)
    REAL,INTENT(IN) :: TT
    AF3 = 0.01*EXP(0.6*(XTT-TT))
  END FUNCTION AF3

! Function for IN concentration reduction between 0 and -20 C:
  REAL FUNCTION REDIN(TT)
    REAL,INTENT(IN) :: TT
    REAL ZZT
    ZZT = MAX(0., MIN(1.,(XTT - TT)/20.))
    REDIN = 1.- (1.-ZZT)/(ZZT**3 + (1.-ZZT)**3)**.333    
  END FUNCTION REDIN

END MODULE MODE_TIWMX_FUN
