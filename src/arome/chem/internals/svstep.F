!     ######spl
      SUBROUTINE SVSTEP (Y, YH, LDYH, YH1, EWT, SAVF, VSAV, ACOR,       &
     &         WM, IWM, F, JAC, PSOL, VNLS, RPAR, IPAR, KMI, KINDEX)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!     ###########################################################
!
      EXTERNAL F, JAC, PSOL, VNLS
      REAL Y, YH, YH1, EWT, SAVF, VSAV, ACOR, WM, RPAR
      INTEGER LDYH, IWM, IPAR
      DIMENSION Y(*), YH(LDYH,*), YH1(*), EWT(*), SAVF(*), VSAV(*),     &
     &   ACOR(*), WM(*), IWM(*), RPAR(*), IPAR(*)
      INTEGER KMI, KINDEX
!-----------------------------------------------------------------------
! Call sequence input -- Y, YH, LDYH, YH1, EWT, SAVF, VSAV,
!                        ACOR, WM, IWM, F, JAC, PSOL, VNLS, RPAR, IPAR
! Call sequence output -- YH, ACOR, WM, IWM
! COMMON block variables accessed..
!     /SVOD01/  ACNRM, EL(13), H, HMIN, HMXI, HNEW, HSCAL, RC, TAU(13),
!               TQ(5), TN, JCUR, JSTART, KFLAG, KUTH,
!               L, LMAX, MAXORD, MITER, N, NEWQ, NQ, NQWAIT
!     /SVOD02/  HU, NCFN, NETF, NFE, NQU, NST
!
! Subroutines called by SVSTEP.. F, SAXPY, CH_SCOPY, SSCAL,
!                               SVJUST, VNLS, SVSET
! Function routines called by SVSTEP.. SVNORM
!-----------------------------------------------------------------------
! SVSTEP performs one step of the integration of an initial value
! problem for a system of ordinary differential equations.
! SVSTEP calls subroutine VNLS for the solution of the nonlinear system
! arising in the time step.  Thus it is independent of the problem
! Jacobian structure and the type of nonlinear system solution method.
! SVSTEP returns a completion flag KFLAG (in COMMON).
! A return with KFLAG = -1 or -2 means either ABS(H) = HMIN or 10
! consecutive failures occurred.  On a return with KFLAG negative,
! the values of TN and the YH array are as of the beginning of the last
! step, and H is the last step size attempted.
!
! Communication with SVSTEP is done with the following variables..
!
! Y      = An array of length N used for the dependent variable vector.
! YH     = An LDYH by LMAX array containing the dependent variables
!          and their approximate scaled derivatives, where
!          LMAX = MAXORD + 1.  YH(i,j+1) contains the approximate
!          j-th derivative of y(i), scaled by H**j/factorial(j)
!          (j = 0,1,...,NQ).  On entry for the first step, the first
!          two columns of YH must be set from the initial values.
! LDYH   = A constant integer .ge. N, the first dimension of YH.
!          N is the number of ODEs in the system.
! YH1    = A one-dimensional array occupying the same space as YH.
! EWT    = An array of length N containing multiplicative weights
!          for local error measurements.  Local errors in y(i) are
!          compared to 1.0/EWT(i) in various error tests.
! SAVF   = An array of working storage, of length N.
!          also used for input of YH(*,MAXORD+2) when JSTART = -1
!          and MAXORD .lt. the current order NQ.
! VSAV   = A work array of length N passed to subroutine VNLS.
! ACOR   = A work array of length N, used for the accumulated
!          corrections.  On a successful return, ACOR(i) contains
!          the estimated one-step local error in y(i).
! WM,IWM = Real and integer work arrays associated with matrix
!          operations in VNLS.
! F      = Dummy name for the user supplied subroutine for f.
! JAC    = Dummy name for the user supplied Jacobian subroutine.
! PSOL   = Dummy name for the subroutine passed to VNLS, for
!          possible use there.
! VNLS   = Dummy name for the nonlinear system solving subroutine,
!          whose real name is dependent on the method used.
! RPAR, IPAR = Dummy names for user's real and integer work arrays.
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
! Type declarations for labeled COMMON block SVOD02 --------------------
!
      REAL HU
      INTEGER NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST
!
! Type declarations for local variables --------------------------------
!
      REAL ADDON, BIAS1,BIAS2,BIAS3, CNQUOT, DDN, DSM, DUP,             &
     &     ETACF, ETAMIN, ETAMX1, ETAMX2, ETAMX3, ETAMXF,               &
     &     ETAQ, ETAQM1, ETAQP1, FLOTL, ONE, ONEPSM,                    &
     &     R, THRESH, TOLD, ZERO
      INTEGER I, I1, I2, IBACK, J, JB, KFC, KFH, MXNCF, NCF, NFLAG
!
! Type declaration for function subroutines called ---------------------
!
      REAL SVNORM
!-----------------------------------------------------------------------
! The following Fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to this integrator.
!-----------------------------------------------------------------------
      SAVE ADDON, BIAS1, BIAS2, BIAS3,                                  &
     &     ETACF, ETAMIN, ETAMX1, ETAMX2, ETAMX3, ETAMXF,               &
     &     KFC, KFH, MXNCF, ONEPSM, THRESH, ONE, ZERO
!-----------------------------------------------------------------------
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
      DATA KFC/-3/, KFH/-7/, MXNCF/10/
      DATA ADDON  /1.0E-6/,    BIAS1  /6.0E0/,     BIAS2  /6.0E0/,      &
     &     BIAS3  /10.0E0/,    ETACF  /0.25E0/,    ETAMIN /0.1E0/,      &
     &     ETAMXF /0.2E0/,     ETAMX1 /1.0E4/,     ETAMX2 /10.0E0/,     &
     &     ETAMX3 /10.0E0/,    ONEPSM /1.00001E0/, THRESH /1.5E0/
      DATA ONE/1.0E0/, ZERO/0.0E0/
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('SVSTEP',0,ZHOOK_HANDLE)
      KFLAG = 0
      TOLD = TN
      NCF = 0
      JCUR = 0
      NFLAG = 0
      IF (JSTART .GT. 0) GO TO 20
      IF (JSTART .EQ. -1) GO TO 100
!-----------------------------------------------------------------------
! On the first call, the order is set to 1, and other variables are
! initialized.  ETAMAX is the maximum ratio by which H can be increased
! in a single step.  It is normally 1.5, but is larger during the
! first 10 steps to compensate for the small initial H.  If a failure
! occurs (in corrector convergence or error test), ETAMAX is set to 1
! for the next increase.
!-----------------------------------------------------------------------
      LMAX = MAXORD + 1
      NQ = 1
      L = 2
      NQNYH = NQ*LDYH
      TAU(1) = H
      PRL1 = ONE
      RC = ZERO
      ETAMAX = ETAMX1
      NQWAIT = 2
      HSCAL = H
      GO TO 200
!-----------------------------------------------------------------------
! Take preliminary actions on a normal continuation step (JSTART.GT.0).
! If the driver changed H, then ETA must be reset and NEWH set to 1.
! If a change of order was dictated on the previous step, then
! it is done here and appropriate adjustments in the history are made.
! On an order decrease, the history array is adjusted by SVJUST.
! On an order increase, the history array is augmented by a column.
! On a change of step size H, the history array YH is rescaled.
!-----------------------------------------------------------------------
 20   CONTINUE
      IF (KUTH .EQ. 1) THEN
        ETA = MIN(ETA,H/HSCAL)
        NEWH = 1
        ENDIF
 50   IF (NEWH .EQ. 0) GO TO 200
      IF (NEWQ .EQ. NQ) GO TO 150
      IF (NEWQ .LT. NQ) THEN
        CALL SVJUST (YH, LDYH, -1)
        NQ = NEWQ
        L = NQ + 1
        NQWAIT = L
        GO TO 150
        ENDIF
      IF (NEWQ .GT. NQ) THEN
        CALL SVJUST (YH, LDYH, 1)
        NQ = NEWQ
        L = NQ + 1
        NQWAIT = L
        GO TO 150
      ENDIF
!-----------------------------------------------------------------------
! The following block handles preliminaries needed when JSTART = -1.
! If N was reduced, zero out part of YH to avoid undefined references.
! If MAXORD was reduced to a value less than the tentative order NEWQ,
! then NQ is set to MAXORD, and a new H ratio ETA is chosen.
! Otherwise, we take the same preliminary actions as for JSTART .gt. 0.
! In any case, NQWAIT is reset to L = NQ + 1 to prevent further
! changes in order for that many steps.
! The new H ratio ETA is limited by the input H if KUTH = 1,
! by HMIN if KUTH = 0, and by HMXI in any case.
! Finally, the history array YH is rescaled.
!-----------------------------------------------------------------------
 100  CONTINUE
      LMAX = MAXORD + 1
      IF (N .EQ. LDYH) GO TO 120
      I1 = 1 + (NEWQ + 1)*LDYH
      I2 = (MAXORD + 1)*LDYH
      IF (I1 .GT. I2) GO TO 120
      DO 110 I = I1, I2
 110    YH1(I) = ZERO
 120  IF (NEWQ .LE. MAXORD) GO TO 140
      FLOTL = REAL(LMAX)
      IF (MAXORD .LT. NQ-1) THEN
        DDN = SVNORM (N, SAVF, EWT)/TQ(1)
        ETA = ONE/((BIAS1*DDN)**(ONE/FLOTL) + ADDON)
        ENDIF
      IF (MAXORD .EQ. NQ .AND. NEWQ .EQ. NQ+1) ETA = ETAQ
      IF (MAXORD .EQ. NQ-1 .AND. NEWQ .EQ. NQ+1) THEN
        ETA = ETAQM1
        CALL SVJUST (YH, LDYH, -1)
        ENDIF
      IF (MAXORD .EQ. NQ-1 .AND. NEWQ .EQ. NQ) THEN
        DDN = SVNORM (N, SAVF, EWT)/TQ(1)
        ETA = ONE/((BIAS1*DDN)**(ONE/FLOTL) + ADDON)
        CALL SVJUST (YH, LDYH, -1)
        ENDIF
      ETA = MIN(ETA,ONE)
      NQ = MAXORD
      L = LMAX
 140  IF (KUTH .EQ. 1) ETA = MIN(ETA,ABS(H/HSCAL))
      IF (KUTH .EQ. 0) ETA = MAX(ETA,HMIN/ABS(HSCAL))
      ETA = ETA/MAX(ONE,ABS(HSCAL)*HMXI*ETA)
      NEWH = 1
      NQWAIT = L
      IF (NEWQ .LE. MAXORD) GO TO 50
! Rescale the history array for a change in H by a factor of ETA. ------
 150  R = ONE
      DO 180 J = 2, L
        R = R*ETA
        CALL SSCAL (N, R, YH(1,J), 1 )
 180    CONTINUE
      H = HSCAL*ETA
      HSCAL = H
      RC = RC*ETA
      NQNYH = NQ*LDYH
!-----------------------------------------------------------------------
! This section computes the predicted values by effectively
! multiplying the YH array by the Pascal triangle matrix.
! SVSET is called to calculate all integration coefficients.
! RC is the ratio of new to old values of the coefficient H/EL(2)=h/l1.
!-----------------------------------------------------------------------
 200  TN = TN + H
      I1 = NQNYH + 1
      DO 220 JB = 1, NQ
        I1 = I1 - LDYH
        DO 210 I = I1, NQNYH
 210      YH1(I) = YH1(I) + YH1(I+LDYH)
 220  CONTINUE
      CALL SVSET
      RL1 = ONE/EL(2)
      RC = RC*(RL1/PRL1)
      PRL1 = RL1
!
! Call the nonlinear system solver. ------------------------------------
!
      CALL VNLS (Y, YH, LDYH, VSAV, SAVF, EWT, ACOR, IWM, WM,           &
     &           F, JAC, PSOL, NFLAG, RPAR, IPAR, KMI, KINDEX)
!
      IF (NFLAG .EQ. 0) GO TO 450
!-----------------------------------------------------------------------
! The VNLS routine failed to achieve convergence (NFLAG .NE. 0).
! The YH array is retracted to its values before prediction.
! The step size H is reduced and the step is retried, if possible.
! Otherwise, an error exit is taken.
!-----------------------------------------------------------------------
        NCF = NCF + 1
        NCFN = NCFN + 1
        ETAMAX = ONE
        TN = TOLD
        I1 = NQNYH + 1
        DO 430 JB = 1, NQ
          I1 = I1 - LDYH
          DO 420 I = I1, NQNYH
 420        YH1(I) = YH1(I) - YH1(I+LDYH)
 430      CONTINUE
        IF (NFLAG .LT. -1) GO TO 680
        IF (ABS(H) .LE. HMIN*ONEPSM) GO TO 670
        IF (NCF .EQ. MXNCF) GO TO 670
        ETA = ETACF
        ETA = MAX(ETA,HMIN/ABS(H))
        NFLAG = -1
        GO TO 150
!-----------------------------------------------------------------------
! The corrector has converged (NFLAG = 0).  The local error test is
! made and control passes to statement 500 if it fails.
!-----------------------------------------------------------------------
 450  CONTINUE
      DSM = ACNRM/TQ(2)
      IF (DSM .GT. ONE) GO TO 500
!-----------------------------------------------------------------------
! After a successful step, update the YH and TAU arrays and decrement
! NQWAIT.  If NQWAIT is then 1 and NQ .lt. MAXORD, then ACOR is saved
! for use in a possible order increase on the next step.
! If ETAMAX = 1 (a failure occurred this step), keep NQWAIT .ge. 2.
!-----------------------------------------------------------------------
      KFLAG = 0
      NST = NST + 1
      HU = H
      NQU = NQ
      DO 470 IBACK = 1, NQ
        I = L - IBACK
 470    TAU(I+1) = TAU(I)
      TAU(1) = H
      DO 480 J = 1, L
        CALL SAXPY (N, EL(J), ACOR, 1, YH(1,J), 1 )
 480    CONTINUE
      NQWAIT = NQWAIT - 1
      IF ((L .EQ. LMAX) .OR. (NQWAIT .NE. 1)) GO TO 490
      CALL CH_SCOPY (N, ACOR, 1, YH(1,LMAX), 1 )
      CONP = TQ(5)
 490  IF (ETAMAX .NE. ONE) GO TO 560
      IF (NQWAIT .LT. 2) NQWAIT = 2
      NEWQ = NQ
      NEWH = 0
      ETA = ONE
      HNEW = H
      GO TO 690
!-----------------------------------------------------------------------
! The error test failed.  KFLAG keeps track of multiple failures.
! Restore TN and the YH array to their previous values, and prepare
! to try the step again.  Compute the optimum step size for the
! same order.  After repeated failures, H is forced to decrease
! more rapidly.
!-----------------------------------------------------------------------
 500  KFLAG = KFLAG - 1
      NETF = NETF + 1
      NFLAG = -2
      TN = TOLD
      I1 = NQNYH + 1
      DO 520 JB = 1, NQ
        I1 = I1 - LDYH
        DO 510 I = I1, NQNYH
 510      YH1(I) = YH1(I) - YH1(I+LDYH)
 520  CONTINUE
      IF (ABS(H) .LE. HMIN*ONEPSM) GO TO 660
      ETAMAX = ONE
      IF (KFLAG .LE. KFC) GO TO 530
! Compute ratio of new H to current H at the current order. ------------
      FLOTL = REAL(L)
      ETA = ONE/((BIAS2*DSM)**(ONE/FLOTL) + ADDON)
      ETA = MAX(ETA,HMIN/ABS(H),ETAMIN)
      IF ((KFLAG .LE. -2) .AND. (ETA .GT. ETAMXF)) ETA = ETAMXF
      GO TO 150
!-----------------------------------------------------------------------
! Control reaches this section if 3 or more consecutive failures
! have occurred.  It is assumed that the elements of the YH array
! have accumulated errors of the wrong order.  The order is reduced
! by one, if possible.  Then H is reduced by a factor of 0.1 and
! the step is retried.  After a total of 7 consecutive failures,
! an exit is taken with KFLAG = -1.
!-----------------------------------------------------------------------
 530  IF (KFLAG .EQ. KFH) GO TO 660
      IF (NQ .EQ. 1) GO TO 540
      ETA = MAX(ETAMIN,HMIN/ABS(H))
      CALL SVJUST (YH, LDYH, -1)
      L = NQ
      NQ = NQ - 1
      NQWAIT = L
      GO TO 150
 540  ETA = MAX(ETAMIN,HMIN/ABS(H))
      H = H*ETA
      HSCAL = H
      TAU(1) = H
!
!*UPG*MNH
!
      CALL F (N, TN, Y, SAVF, RPAR, IPAR, KMI, KINDEX)
!
!*UPG*MNH
!
      NFE = NFE + 1
      DO 550 I = 1, N
 550    YH(I,2) = H*SAVF(I)
      NQWAIT = 10
      GO TO 200
!-----------------------------------------------------------------------
! If NQWAIT = 0, an increase or decrease in order by one is considered.
! Factors ETAQ, ETAQM1, ETAQP1 are computed by which H could
! be multiplied at order q, q-1, or q+1, respectively.
! The largest of these is determined, and the new order and
! step size set accordingly.
! A change of H or NQ is made only if H increases by at least a
! factor of THRESH.  If an order change is considered and rejected,
! then NQWAIT is set to 2 (reconsider it after 2 steps).
!-----------------------------------------------------------------------
! Compute ratio of new H to current H at the current order. ------------
 560  FLOTL = REAL(L)
      ETAQ = ONE/((BIAS2*DSM)**(ONE/FLOTL) + ADDON)
      IF (NQWAIT .NE. 0) GO TO 600
      NQWAIT = 2
      ETAQM1 = ZERO
      IF (NQ .EQ. 1) GO TO 570
! Compute ratio of new H to current H at the current order less one. ---
      DDN = SVNORM (N, YH(1,L), EWT)/TQ(1)
      ETAQM1 = ONE/((BIAS1*DDN)**(ONE/(FLOTL - ONE)) + ADDON)
 570  ETAQP1 = ZERO
      IF (L .EQ. LMAX) GO TO 580
! Compute ratio of new H to current H at current order plus one. -------
      CNQUOT = (TQ(5)/CONP)*(H/TAU(2))**L
      DO 575 I = 1, N
 575    SAVF(I) = ACOR(I) - CNQUOT*YH(I,LMAX)
      DUP = SVNORM (N, SAVF, EWT)/TQ(3)
      ETAQP1 = ONE/((BIAS3*DUP)**(ONE/(FLOTL + ONE)) + ADDON)
 580  IF (ETAQ .GE. ETAQP1) GO TO 590
      IF (ETAQP1 .GT. ETAQM1) GO TO 620
      GO TO 610
 590  IF (ETAQ .LT. ETAQM1) GO TO 610
 600  ETA = ETAQ
      NEWQ = NQ
      GO TO 630
 610  ETA = ETAQM1
      NEWQ = NQ - 1
      GO TO 630
 620  ETA = ETAQP1
      NEWQ = NQ + 1
      CALL CH_SCOPY (N, ACOR, 1, YH(1,LMAX), 1)
! Test tentative new H against THRESH, ETAMAX, and HMXI, then exit. ----
 630  IF (ETA .LT. THRESH .OR. ETAMAX .EQ. ONE) GO TO 640
      ETA = MIN(ETA,ETAMAX)
      ETA = ETA/MAX(ONE,ABS(H)*HMXI*ETA)
      NEWH = 1
      HNEW = H*ETA
      GO TO 690
 640  NEWQ = NQ
      NEWH = 0
      ETA = ONE
      HNEW = H
      GO TO 690
!-----------------------------------------------------------------------
! All returns are made through this section.
! On a successful return, ETAMAX is reset and ACOR is scaled.
!-----------------------------------------------------------------------
 660  KFLAG = -1
      GO TO 720
 670  KFLAG = -2
      GO TO 720
 680  IF (NFLAG .EQ. -2) KFLAG = -3
      IF (NFLAG .EQ. -3) KFLAG = -4
      GO TO 720
 690  ETAMAX = ETAMX3
      IF (NST .LE. 10) ETAMAX = ETAMX2
 700  R = ONE/TQ(2)
      CALL SSCAL (N, R, ACOR, 1)
 720  JSTART = 1
      IF (LHOOK) CALL DR_HOOK('SVSTEP',1,ZHOOK_HANDLE)
      RETURN
      END SUBROUTINE SVSTEP
