!     ######spl
      SUBROUTINE  INV_LEVELS(KIDIA,KFDIA,KLEV,KINCR,PARRAYALD,PARRAYMNH)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!     ##########################################################################
!
!!****  * -  transforms  from MNH (KLEV+2) arrays to
!!           Aladin (KLEV) arrays and reversely
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
!!      "derived from mpa/micro/externals/invert_vlev.f90"
!!
!!    AUTHOR
!!    ------
!!    Y.Seity
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    18/09/12
!!    M. Mokhtari & A. Ambar 09/2016: adptation for dust in both version aladin and arome                         
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!

!
INTEGER,                  INTENT(IN)   :: KINCR  ! =1 if ALD->MNH and -1 if MNH->ALD
INTEGER,                  INTENT(IN)   :: KIDIA,KFDIA
INTEGER,                  INTENT(IN)   :: KLEV     !Number of vertical levels
REAL, DIMENSION(KFDIA,KLEV),   INTENT(INOUT) :: PARRAYALD  ! Aladin type arry
REAL, DIMENSION(KFDIA,1,KLEV+2), INTENT(INOUT) :: PARRAYMNH  ! MNH type array

!*      0.2   Local variables:
INTEGER:: JLEV
INTEGER:: JLON
!---------------------------------------------------------------

!write(*,*) 'debut inv lev'

REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('INV_LEVELS',0,ZHOOK_HANDLE)
IF (KINCR==1) THEN

       DO JLEV = 2 , KLEV+1
           DO JLON = KIDIA,KFDIA
              PARRAYMNH(JLON,1,JLEV)= PARRAYALD(JLON,KLEV+2-JLEV)
           ENDDO
       ENDDO
       DO JLON = KIDIA,KFDIA
           PARRAYMNH(JLON,1,1)= PARRAYMNH(JLON,1,2)
           PARRAYMNH(JLON,1,KLEV+2)= PARRAYMNH(JLON,1,KLEV+1)
       ENDDO



ElSE IF  (KINCR==-1) THEN
       DO JLEV = 1 ,KLEV
            DO JLON = KIDIA,KFDIA
              PARRAYALD(JLON,JLEV)=PARRAYMNH(JLON,1,KLEV+2-JLEV)
            ENDDO
       ENDDO

ELSE
  WRITE(*,*) 'VALUE FOR KINCR NE TO 1 OR -1 !!!'
  STOP
ENDIF



IF (LHOOK) CALL DR_HOOK('INV_LEVELS',1,ZHOOK_HANDLE)
END SUBROUTINE INV_LEVELS
