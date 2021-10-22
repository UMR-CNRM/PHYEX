!     ######spl
       SUBROUTINE CH_JVALUES_n(KVECNPT, KIINDEX, KJINDEX, KMODELLEVEL, KRATES, PJVALUES, PRATES)
       USE PARKIND1, ONLY : JPRB
       USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!      ###############################################################
! 
!!
!!*** *CH_JVALUES* extracts photolysis rates
!!
!!    PURPOSE
!!    -------
!!      Extract photolysis rates from MODD_CH_JVALUES at a given model level
!!
!!**  METHOD
!!    ------
!!      After the photolysis rates have been set with CH_UPDATE_JVALUES_n,
!!    the different photolysis rates may be extracted on the model levels.
!       
!!    REFERENCE
!!    ---------
!!    MesoNH documentation
!!
!!    AUTHOR
!!    ------
!!    Karsten Suhre (LA)
!!    
!!    MODIFICATIONS
!!    -------------
!!    Original 05/03/97
!!    C. Mari  20/03/01 3D version of J values + vectorisation
!!    Modification   01/12/03  (Gazen)   change Chemical scheme interface
!!
!!------------------------------------------------------------------------------
!!
!!    EXTERNAL
!!    --------
!!    none
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
USE MODD_CH_INIT_JVALUES, ONLY: JPJVMAX
!!
!!    EXPLICIT ARGUMENTS
!!    ------------------
IMPLICIT NONE
INTEGER, DIMENSION(:),   INTENT(IN)  :: KMODELLEVEL
INTEGER,                 INTENT(IN)  :: KVECNPT 
INTEGER, DIMENSION(:),   INTENT(IN)  :: KIINDEX,KJINDEX 
                                      ! current model grid point
INTEGER,                 INTENT(IN)  :: KRATES          
                                      ! dimension of PRATES
REAL, DIMENSION(KVECNPT,KRATES), INTENT(OUT) :: PRATES          
                                      ! photolysis rates 
REAL,DIMENSION(:,:,:,:), INTENT(IN) :: PJVALUES    ! Tuv coefficients
!!
!!    LOCAL VARIABLES
!!    ---------------
INTEGER :: JI,JN
!!
!------------------------------------------------------------------------------
!!
!!    EXECUTABLE STATEMENTS
!!    ---------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('CH_JVALUES_N',0,ZHOOK_HANDLE)
DO JI = 1, KRATES
 DO JN = 1, KVECNPT
   PRATES(JN,JI) = PJVALUES(KIINDEX(JN),KJINDEX(JN),KMODELLEVEL(JN),JI)
 ENDDO
ENDDO
!
IF (LHOOK) CALL DR_HOOK('CH_JVALUES_N',1,ZHOOK_HANDLE)
END SUBROUTINE CH_JVALUES_n
