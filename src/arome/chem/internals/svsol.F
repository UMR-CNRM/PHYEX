!     ######spl
      SUBROUTINE SVSOL (WM, IWM, X, IERSL)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!     ####################################
      REAL WM, X
      INTEGER IWM, IERSL
      DIMENSION WM(*), IWM(*), X(*)
!-----------------------------------------------------------------------
! Call sequence input -- WM, IWM, X
! Call sequence output -- X, IERSL
! COMMON block variables accessed..
!     /SVOD01/ -- H, RL1, MITER, N
!
! Subroutines called by SVSOL.. SGESL, SGBSL
! Function routines called by SVSOL.. None
!-----------------------------------------------------------------------
! This routine manages the solution of the linear system arising from
! a chord iteration.  It is called if MITER .ne. 0.
! If MITER is 1 or 2, it calls SGESL to accomplish this.
! If MITER = 3 it updates the coefficient H*RL1 in the diagonal
! matrix, and then computes the solution.
! If MITER is 4 or 5, it calls SGBSL.
! Communication with SVSOL uses the following variables..
! WM    = Real work space containing the inverse diagonal matrix if
!         MITER = 3 and the LU decomposition of the matrix otherwise.
!         Storage of matrix elements starts at WM(3).
!         WM also contains the following matrix-related data..
!         WM(1) = SQRT(UROUND) (not used here),
!         WM(2) = HRL1, the previous value of H*RL1, used if MITER = 3.
! IWM   = Integer work space containing pivot information, starting at
!         IWM(31), if MITER is 1, 2, 4, or 5.  IWM also contains band
!         parameters ML = IWM(1) and MU = IWM(2) if MITER is 4 or 5.
! X     = The right-hand side vector on input, and the solution vector
!         on output, of length N.
! IERSL = Output flag.  IERSL = 0 if no trouble occurred.
!         IERSL = 1 if a singular matrix arose with MITER = 3.
!-----------------------------------------------------------------------
!
! Type declarations for labeled COMMON block SVOD01 --------------------
!
      REAL ACNRM, CCMXJ, CONP, CRATE, DRC, EL,                          &
     &     ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,               &
     &     RC, RL1, TAU, TQ, TN, UROUND
      INTEGER ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,          &
     &        L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,              &
     &        LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,         &
     &        N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,            &
     &        NSLP, NYH
!
! Type declarations for local variables --------------------------------
!
      INTEGER I, MEBAND, ML, MU
      REAL DI, HRL1, ONE, PHRL1, R, ZERO
!-----------------------------------------------------------------------
! The following Fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to this integrator.
!-----------------------------------------------------------------------
      SAVE ONE, ZERO
!
      COMMON /SVOD01/ ACNRM, CCMXJ, CONP, CRATE, DRC, EL(13),           &
     &                ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,    &
     &                RC, RL1, TAU(13), TQ(5), TN, UROUND,              &
     &                ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,  &
     &                L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,      &
     &                LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP, &
     &                N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,    &
     &                NSLP, NYH
!
      DATA ONE /1.0E0/, ZERO /0.0E0/
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('SVSOL',0,ZHOOK_HANDLE)
      IERSL = 0
      GO TO (100, 100, 300, 400, 400), MITER
 100  CALL SGESL (WM(3), N, N, IWM(31), X, 0)
      IF (LHOOK) CALL DR_HOOK('SVSOL',1,ZHOOK_HANDLE)
      RETURN
!
 300  PHRL1 = WM(2)
      HRL1 = H*RL1
      WM(2) = HRL1
      IF (HRL1 .EQ. PHRL1) GO TO 330
      R = HRL1/PHRL1
      DO 320 I = 1,N
        DI = ONE - R*(ONE - ONE/WM(I+2))
        IF (ABS(DI) .EQ. ZERO) GO TO 390
 320    WM(I+2) = ONE/DI
!
 330  DO 340 I = 1,N
 340    X(I) = WM(I+2)*X(I)
      IF (LHOOK) CALL DR_HOOK('SVSOL',1,ZHOOK_HANDLE)
      RETURN
 390  IERSL = 1
      IF (LHOOK) CALL DR_HOOK('SVSOL',1,ZHOOK_HANDLE)
      RETURN
!
 400  ML = IWM(1)
      MU = IWM(2)
      MEBAND = 2*ML + MU + 1
      CALL SGBSL (WM(3), MEBAND, N, ML, MU, IWM(31), X, 0)
      IF (LHOOK) CALL DR_HOOK('SVSOL',1,ZHOOK_HANDLE)
      RETURN
      END SUBROUTINE SVSOL
