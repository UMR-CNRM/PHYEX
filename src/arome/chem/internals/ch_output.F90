!     ######spl
      SUBROUTINE CH_OUTPUT(PCONC, PAERO, PMI, TPM, KMI, KVECNPT)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!!    ##############################################
!!
!!***  *CH_OUTPUT*
!!
!!    PURPOSE
!!    -------
!!    save variables to disk for later visualisation
!!
!!    METHOD
!!    ------
!!    write all variables specified in CH_INIT_OUTPUT in the format
!!    CRESULTFORMAT to unit NRESULTIO, which should have been opened
!!    previously by CH_INIT_OUTPUT
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!    K. Suhre
!!
!!    MODIFICATIONS
!!    -------------
!!    Original 02/03/95
!!    27/07/96 (K. Suhre) restructured
!!    01/12/03 (D. Gazen)   change Chemical scheme interface
!!
!!    EXTERNAL
!!    --------
USE MODI_CH_GET_RATES
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
USE MODD_CH_MODEL0D,  ONLY: NRESULTIO, XTSIMUL, CRESULTFORMAT, NVERB
USE MODD_CH_AEROSOL,  ONLY: LORILAM
USE MODD_CH_M9,       ONLY: NEQ, NREAC, NMETEOVARS, METEOTRANSTYPE  
!!

!!
!------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!        -----------------
IMPLICIT NONE
!
!*       0.1  declaration of arguments
!
INTEGER, INTENT(IN)              :: KVECNPT
REAL, DIMENSION(:,:), INTENT(IN) :: PCONC ! the chem species
REAL, DIMENSION(:,:), INTENT(IN) :: PAERO ! the chem species
REAL, DIMENSION(:,:), INTENT(IN) :: PMI   ! the molecular mass
TYPE(METEOTRANSTYPE), DIMENSION(KVECNPT), INTENT(IN) :: TPM   ! the meteo variables
INTEGER, INTENT(IN)              :: KMI
!
!*       0.2  declaration of local variables
!!
INTEGER                :: JI     ! loop control
!REAL, DIMENSION(KVECNPT,NREAC) :: ZRATE  ! for retrieval of reaction rates and names
REAL, DIMENSION(:,:), ALLOCATABLE :: ZRATE  ! for retrieval of reaction rates and names
INTEGER :: NAERO, NMI
!
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('CH_OUTPUT',0,ZHOOK_HANDLE)
ALLOCATE(ZRATE(KVECNPT,NREAC))
!
!*    EXECUTABLE STATEMENTS
!     ---------------------
!
IF (LORILAM) THEN
NAERO=SIZE(PAERO,2)
NMI=SIZE(PMI,2)
ELSE
NAERO=0
NMI=0
END IF
IF (NVERB >= 5) &
   PRINT *, "CH_OUTPUT: saving results to disk at XTSIMUL = ", XTSIMUL
!
! get reaction rates
CALL CH_GET_RATES(ZRATE,KMI,KVECNPT,NREAC)
!
! write all variables to file
WRITE(NRESULTIO,CRESULTFORMAT) XTSIMUL,        &
       (PCONC(1,JI), JI = 1, NEQ),             &
       (PAERO(1,JI), JI = 1, NAERO),           &
       (PMI(1,JI), JI = 1, NMI),               &
       (ZRATE(1,JI), JI = 1, NREAC),           &
       (TPM(1)%XMETEOVAR(JI), JI = 1, NMETEOVARS)
!
DEALLOCATE(ZRATE)
!
IF (LHOOK) CALL DR_HOOK('CH_OUTPUT',1,ZHOOK_HANDLE)
END SUBROUTINE CH_OUTPUT
