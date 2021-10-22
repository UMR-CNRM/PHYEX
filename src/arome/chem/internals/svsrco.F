!     ######spl
      SUBROUTINE SVSRCO (RSAV, ISAV, JOB)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!     ###################################
      REAL RSAV
      INTEGER ISAV, JOB
      DIMENSION RSAV(*), ISAV(*)
!-----------------------------------------------------------------------
! Call sequence input -- RSAV, ISAV, JOB
! Call sequence output -- RSAV, ISAV
! COMMON block variables accessed -- All of /SVOD01/ and /SVOD02/
!
! Subroutines/functions called by SVSRCO.. None
!-----------------------------------------------------------------------
! This routine saves or restores (depending on JOB) the contents of the
! COMMON blocks SVOD01 and SVOD02, which are used internally by SVODE.
!
! RSAV = real array of length 49 or more.
! ISAV = integer array of length 41 or more.
! JOB  = flag indicating to save or restore the COMMON blocks..
!        JOB  = 1 if COMMON is to be saved (written to RSAV/ISAV).
!        JOB  = 2 if COMMON is to be restored (read from RSAV/ISAV).
!        A call with JOB = 2 presumes a prior call with JOB = 1.
!-----------------------------------------------------------------------
      REAL RVOD1, RVOD2
      INTEGER IVOD1, IVOD2
      INTEGER I, LENIV1, LENIV2, LENRV1, LENRV2
!-----------------------------------------------------------------------
! The following Fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to this integrator.
!-----------------------------------------------------------------------
      SAVE LENRV1, LENIV1, LENRV2, LENIV2
!
      COMMON /SVOD01/ RVOD1(48), IVOD1(33)
      COMMON /SVOD02/ RVOD2(1), IVOD2(8)
      DATA LENRV1/48/, LENIV1/33/, LENRV2/1/, LENIV2/8/
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('SVSRCO',0,ZHOOK_HANDLE)
      IF (JOB .EQ. 2) GO TO 100
      DO 10 I = 1,LENRV1
 10     RSAV(I) = RVOD1(I)
      DO 15 I = 1,LENRV2
 15     RSAV(LENRV1+I) = RVOD2(I)
!
      DO 20 I = 1,LENIV1
 20     ISAV(I) = IVOD1(I)
      DO 25 I = 1,LENIV2
 25     ISAV(LENIV1+I) = IVOD2(I)
!
      IF (LHOOK) CALL DR_HOOK('SVSRCO',1,ZHOOK_HANDLE)
      RETURN
!
 100  CONTINUE
      DO 110 I = 1,LENRV1
 110     RVOD1(I) = RSAV(I)
      DO 115 I = 1,LENRV2
 115     RVOD2(I) = RSAV(LENRV1+I)
!
      DO 120 I = 1,LENIV1
 120     IVOD1(I) = ISAV(I)
      DO 125 I = 1,LENIV2
 125     IVOD2(I) = ISAV(LENIV1+I)
!
      IF (LHOOK) CALL DR_HOOK('SVSRCO',1,ZHOOK_HANDLE)
      RETURN
      END SUBROUTINE SVSRCO
