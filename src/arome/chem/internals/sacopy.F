!     ######spl
      SUBROUTINE SACOPY (NROW, NCOL, A, NROWA, B, NROWB)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!     ##################################################
      REAL A, B
      INTEGER NROW, NCOL, NROWA, NROWB
      DIMENSION A(NROWA,NCOL), B(NROWB,NCOL)
!-----------------------------------------------------------------------
! Call sequence input -- NROW, NCOL, A, NROWA, NROWB
! Call sequence output -- B
! COMMON block variables accessed -- None
!
! Subroutines called by SACOPY.. CH_SCOPY
! Function routines called by SACOPY.. None
!-----------------------------------------------------------------------
! This routine copies one rectangular array, A, to another, B,
! where A and B may have different row dimensions, NROWA and NROWB.
! The data copied consists of NROW rows and NCOL columns.
!-----------------------------------------------------------------------
      INTEGER IC
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('SACOPY',0,ZHOOK_HANDLE)
      DO 20 IC = 1,NCOL
        CALL CH_SCOPY (NROW, A(1,IC), 1, B(1,IC), 1)
 20     CONTINUE
!
      IF (LHOOK) CALL DR_HOOK('SACOPY',1,ZHOOK_HANDLE)
      RETURN
      END SUBROUTINE SACOPY
