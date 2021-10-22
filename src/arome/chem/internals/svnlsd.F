!     ######spl
      SUBROUTINE SVNLSD (Y, YH, LDYH, VSAV, SAVF, EWT, ACOR, IWM, WM,   &
     &                 F, JAC, PDUM, NFLAG, RPAR, IPAR, KMI, KINDEX)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!     ###############################################################
      EXTERNAL F, JAC, PDUM
      REAL Y, YH, VSAV, SAVF, EWT, ACOR, WM, RPAR
      INTEGER LDYH, IWM, NFLAG, IPAR
      DIMENSION Y(*), YH(LDYH,*), VSAV(*), SAVF(*), EWT(*), ACOR(*),    &
     &          IWM(*), WM(*), RPAR(*), IPAR(*)
      INTEGER KMI,KINDEX
!-----------------------------------------------------------------------
! Call sequence input -- Y, YH, LDYH, SAVF, EWT, ACOR, IWM, WM,
!                        F, JAC, NFLAG, RPAR, IPAR
! Call sequence output -- YH, ACOR, WM, IWM, NFLAG
! COMMON block variables accessed..
!     /SVOD01/ ACNRM, CRATE, DRC, H, RC, RL1, TQ(5), TN, ICF,
!                JCUR, METH, MITER, N, NSLP
!     /SVOD02/ HU, NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST
!
! Subroutines called by SVNLSD.. F, SAXPY, CH_SCOPY, SSCAL, SVJAC, SVSOL
! Function routines called by SVNLSD.. SVNORM
!-----------------------------------------------------------------------
! Subroutine SVNLSD is a nonlinear system solver, which uses functional
! iteration or a chord (modified Newton) method.  For the chord method
! direct linear algebraic system solvers are used.  Subroutine SVNLSD
! then handles the corrector phase of this integration package.
!
! Communication with SVNLSD is done with the following variables. (For
! more details, please see the comments in the driver subroutine.)
!
! Y          = The dependent variable, a vector of length N, input.
! YH         = The Nordsieck (Taylor) array, LDYH by LMAX, input
!              and output.  On input, it contains predicted values.
! LDYH       = A constant .ge. N, the first dimension of YH, input.
! VSAV       = Unused work array.
! SAVF       = A work array of length N.
! EWT        = An error weight vector of length N, input.
! ACOR       = A work array of length N, used for the accumulated
!              corrections to the predicted y vector.
! WM,IWM     = Real and integer work arrays associated with matrix
!              operations in chord iteration (MITER .ne. 0).
! F          = Dummy name for user supplied routine for f.
! JAC        = Dummy name for user supplied Jacobian routine.
! PDUM       = Unused dummy subroutine name.  Included for uniformity
!              over collection of integrators.
! NFLAG      = Input/output flag, with values and meanings as follows..
!              INPUT
!                  0 first call for this time step.
!                 -1 convergence failure in previous call to SVNLSD.
!                 -2 error test failure in SVSTEP.
!              OUTPUT
!                  0 successful completion of nonlinear solver.
!                 -1 convergence failure or singular matrix.
!                 -2 unrecoverable error in matrix preprocessing
!                    (cannot occur here).
!                 -3 unrecoverable error in solution (cannot occur
!                    here).
! RPAR, IPAR = Dummy names for user's real and integer work arrays.
!
! IPUP       = Own variable flag with values and meanings as follows..
!              0,            do not update the Newton matrix.
!              MITER .ne. 0, update Newton matrix, because it is the
!                            initial step, order was changed, the error
!                            test failed, or an update is indicated by
!                            the scalar RC or step counter NST.
!
! For more details, see comments in driver subroutine.
!-----------------------------------------------------------------------
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
! Type declarations for labeled COMMON block SVOD02 --------------------
!
      REAL HU
      INTEGER NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST
!
! Type declarations for local variables --------------------------------
!
      REAL CCMAX, CRDOWN, CSCALE, DCON, DEL, DELP, ONE,                 &
     &     RDIV, TWO, ZERO
      INTEGER I, IERPJ, IERSL, M, MAXCOR, MSBP
!
! Type declaration for function subroutines called ---------------------
!
      REAL SVNORM
!-----------------------------------------------------------------------
! The following Fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to this integrator.
!-----------------------------------------------------------------------
      SAVE CCMAX, CRDOWN, MAXCOR, MSBP, RDIV, ONE, TWO, ZERO
!
      COMMON /SVOD01/ ACNRM, CCMXJ, CONP, CRATE, DRC, EL(13),           &
     &                ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,    &
     &                RC, RL1, TAU(13), TQ(5), TN, UROUND,              &
     &                ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,  &
     &                L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,      &
     &                LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP, &
     &                N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,    &
     &                NSLP, NYH
      COMMON /SVOD02/ HU, NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST
!
      DATA CCMAX /0.3E0/, CRDOWN /0.3E0/, MAXCOR /3/, MSBP /20/,        &
     &     RDIV  /2.0E0/
      DATA ONE /1.0E0/, TWO /2.0E0/, ZERO /0.0E0/
!-----------------------------------------------------------------------
! On the first step, on a change of method order, or after a
! nonlinear convergence failure with NFLAG = -2, set IPUP = MITER
! to force a Jacobian update when MITER .ne. 0.
!-----------------------------------------------------------------------
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('SVNLSD',0,ZHOOK_HANDLE)
      IF (JSTART .EQ. 0) NSLP = 0
      IF (NFLAG .EQ. 0) ICF = 0
      IF (NFLAG .EQ. -2) IPUP = MITER
      IF ( (JSTART .EQ. 0) .OR. (JSTART .EQ. -1) ) IPUP = MITER
! If this is functional iteration, set CRATE .eq. 1 and drop to 220
      IF (MITER .EQ. 0) THEN
        CRATE = ONE
        GO TO 220
      ENDIF
!-----------------------------------------------------------------------
! RC is the ratio of new to old values of the coefficient H/EL(2)=h/l1.
! When RC differs from 1 by more than CCMAX, IPUP is set to MITER
! to force SVJAC to be called, if a Jacobian is involved.
! In any case, SVJAC is called at least every MSBP steps.
!-----------------------------------------------------------------------
      DRC = ABS(RC-ONE)
      IF (DRC .GT. CCMAX .OR. NST .GE. NSLP+MSBP) IPUP = MITER
!-----------------------------------------------------------------------
! Up to MAXCOR corrector iterations are taken.  A convergence test is
! made on the r.m.s. norm of each correction, weighted by the error
! weight vector EWT.  The sum of the corrections is accumulated in the
! vector ACOR(i).  The YH array is not altered in the corrector loop.
!-----------------------------------------------------------------------
 220  M = 0
      DELP = ZERO
      CALL CH_SCOPY (N, YH(1,1), 1, Y, 1 )
!
!*UPG*MNH
!
      CALL F (N, TN, Y, SAVF, RPAR, IPAR, KMI, KINDEX)
!
!*UPG*MNH
!
      NFE = NFE + 1
      IF (IPUP .LE. 0) GO TO 250
!-----------------------------------------------------------------------
! If indicated, the matrix P = I - h*rl1*J is reevaluated and
! preprocessed before starting the corrector iteration.  IPUP is set
! to 0 as an indicator that this has been done.
!-----------------------------------------------------------------------
      CALL SVJAC (Y, YH, LDYH, EWT, ACOR, SAVF, WM, IWM, F, JAC, IERPJ, &
     &           RPAR, IPAR, KMI, KINDEX)
      IPUP = 0
      RC = ONE
      DRC = ZERO
      CRATE = ONE
      NSLP = NST
! If matrix is singular, take error return to force cut in step size. --
      IF (IERPJ .NE. 0) GO TO 430
 250  DO 260 I = 1,N
 260    ACOR(I) = ZERO
! This is a looping point for the corrector iteration. -----------------
 270  IF (MITER .NE. 0) GO TO 350
!-----------------------------------------------------------------------
! In the case of functional iteration, update Y directly from
! the result of the last function evaluation.
!-----------------------------------------------------------------------
      DO 280 I = 1,N
 280    SAVF(I) = RL1*(H*SAVF(I) - YH(I,2))
      DO 290 I = 1,N
 290    Y(I) = SAVF(I) - ACOR(I)
      DEL = SVNORM (N, Y, EWT)
      DO 300 I = 1,N
 300    Y(I) = YH(I,1) + SAVF(I)
      CALL CH_SCOPY (N, SAVF, 1, ACOR, 1)
      GO TO 400
!-----------------------------------------------------------------------
! In the case of the chord method, compute the corrector error,
! and solve the linear system with that as right-hand side and
! P as coefficient matrix.  The correction is scaled by the factor
! 2/(1+RC) to account for changes in h*rl1 since the last SVJAC call.
!-----------------------------------------------------------------------
 350  DO 360 I = 1,N
 360    Y(I) = (RL1*H)*SAVF(I) - (RL1*YH(I,2) + ACOR(I))
      CALL SVSOL (WM, IWM, Y, IERSL)
      NNI = NNI + 1
      IF (IERSL .GT. 0) GO TO 410
      IF (METH .EQ. 2 .AND. RC .NE. ONE) THEN
        CSCALE = TWO/(ONE + RC)
        CALL SSCAL (N, CSCALE, Y, 1)
      ENDIF
      DEL = SVNORM (N, Y, EWT)
      CALL SAXPY (N, ONE, Y, 1, ACOR, 1)
      DO 380 I = 1,N
 380    Y(I) = YH(I,1) + ACOR(I)
!-----------------------------------------------------------------------
! Test for convergence.  If M .gt. 0, an estimate of the convergence
! rate constant is stored in CRATE, and this is used in the test.
!-----------------------------------------------------------------------
 400  IF (M .NE. 0) CRATE = MAX(CRDOWN*CRATE,DEL/DELP)
      DCON = DEL*MIN(ONE,CRATE)/TQ(4)
      IF (DCON .LE. ONE) GO TO 450
      M = M + 1
      IF (M .EQ. MAXCOR) GO TO 410
      IF (M .GE. 2 .AND. DEL .GT. RDIV*DELP) GO TO 410
      DELP = DEL
!
!*UPG*MNH
!
      CALL F (N, TN, Y, SAVF, RPAR, IPAR, KMI, KINDEX)
!
!*UPG*MNH
!
      NFE = NFE + 1
      GO TO 270
!
 410  IF (MITER .EQ. 0 .OR. JCUR .EQ. 1) GO TO 430
      ICF = 1
      IPUP = MITER
      GO TO 220
!
 430  CONTINUE
      NFLAG = -1
      ICF = 2
      IPUP = MITER
      IF (LHOOK) CALL DR_HOOK('SVNLSD',1,ZHOOK_HANDLE)
      RETURN
!
! Return for successful step. ------------------------------------------
 450  NFLAG = 0
      JCUR = 0
      ICF = 0
      IF (M .EQ. 0) ACNRM = DEL
      IF (M .GT. 0) ACNRM = SVNORM (N, ACOR, EWT)
      IF (LHOOK) CALL DR_HOOK('SVNLSD',1,ZHOOK_HANDLE)
      RETURN
      END SUBROUTINE SVNLSD
