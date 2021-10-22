!     ######spl
     SUBROUTINE CH_AER_EQM_CORMASS(PSVT) 
     USE PARKIND1, ONLY : JPRB
     USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!!   ############################################################
!!
!!    PURPOSE
!!    -------
!!    Realise la conservation de la masse 
!!    Filtre les valeurs des moments 0 et 6 inferieures aux valeurs 
!!    minimales (Rg et SIG) introduites en nameliste 
!!
!!    REFERENCE
!!    ---------
!!    none
!!
!!    AUTHOR
!!    ------
!!    Pierre TULET (GMEI)
!!
!!    MODIFICATIONS
!!    -------------
!!    none
!!
!!    EXTERNAL
!!    --------
!!    None
!!
!!
IMPLICIT NONE
REAL,       DIMENSION(:,:,:,:),  INTENT(INOUT)   :: PSVT
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
!*      0.1    declarations of arguments
!
!
!*      0.2    declarations local variables
!
!-------------------------------------------------------------------------------

!
!*       1.     TRANSFER FROM GAS TO AEROSOL MODULE
!               ------------------------------------
!
  REAL(KIND=JPRB) :: ZHOOK_HANDLE
  IF (LHOOK) CALL DR_HOOK('CH_AER_EQM_CORMASS',0,ZHOOK_HANDLE)
  PSVT(:,:,:,:) = MAX(PSVT(:,:,:,:),1.E-80)

!
IF (LHOOK) CALL DR_HOOK('CH_AER_EQM_CORMASS',1,ZHOOK_HANDLE)
END SUBROUTINE CH_AER_EQM_CORMASS
