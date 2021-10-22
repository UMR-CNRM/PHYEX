      MODULE MODI_INV_LEVELS
!     ########################
INTERFACE
      SUBROUTINE  INV_LEVELS(KIDIA,KFDIA,KLEV,KINCR,PARRAYALD,  &
                                  PARRAYMNH)
      USE PARKIND1, ONLY : JPRB, JPIM
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!     ##########################################################################
!
!!****  * -  transforms  from MNH (KLEV+2) arrays to
!!      Aladin (KLEV) arrays and reversely
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to reshuffle Aladin (MNH) arrays to
!!      MNH (Aladin) ones
!!
!!**  METHOD
!!    ------
!!    To go from MNH to ALADIN (KINCR=1): 
!!    Suppress first and last level. levels
!!    To go from ALADIN to MNH (KINCR=-1):
!!    Add a vertical level at top and bottom and then reverse orders
!!
!!
!!
!!    EXTERNAL
!!    --------
!!      None
!!     
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      None     
!! 
!!    REFERENCE 
!!    ---------
!!
!!      Documentation "New data flux for diagnostics in Arome"
!!
!!    AUTHOR
!!    ------
!!    O.Riviere
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    18/06/08

!!                                  
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!

!
INTEGER(KIND=JPIM),                  INTENT(IN)   :: KINCR  ! =1 if ALD->MNH and -1 if MNH->ALD 
INTEGER(KIND=JPIM),                  INTENT(IN)   :: KIDIA,KFDIA
INTEGER(KIND=JPIM),                  INTENT(IN)   :: KLEV     !Number of vertical levels 
REAL(KIND=JPRB), DIMENSION(KFDIA,KLEV),   INTENT(INOUT) :: PARRAYALD  ! Aladin type arry
REAL(KIND=JPRB), DIMENSION(KFDIA,1,KLEV+2), INTENT(INOUT) :: PARRAYMNH  ! MNH type array

END SUBROUTINE INV_LEVELS
END INTERFACE
END MODULE MODI_INV_LEVELS


