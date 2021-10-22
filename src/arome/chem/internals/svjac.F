!     ######spl
      SUBROUTINE SVJAC (Y, YH, LDYH, EWT, FTEM, SAVF, WM, IWM, F, JAC,  &
     &                 IERPJ, RPAR, IPAR, KMI, KINDEX)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!     ################################################################
      EXTERNAL F, JAC
      REAL Y, YH, EWT, FTEM, SAVF, WM, RPAR
      INTEGER LDYH, IWM, IERPJ, IPAR
      DIMENSION Y(*), YH(LDYH,*), EWT(*), FTEM(*), SAVF(*),             &
     &   WM(*), IWM(*), RPAR(*), IPAR(*)
      INTEGER KMI, KINDEX
!-----------------------------------------------------------------------
! Call sequence input -- Y, YH, LDYH, EWT, FTEM, SAVF, WM, IWM,
!                        F, JAC, RPAR, IPAR
! Call sequence output -- WM, IWM, IERPJ
! COMMON block variables accessed..
!     /SVOD01/  CCMXJ, DRC, H, RL1, TN, UROUND, ICF, JCUR, LOCJS,
!               MSBJ, NSLJ
!     /SVOD02/  NFE, NST, NJE, NLU
!
! Subroutines called by SVJAC.. F, JAC, SACOPY, CH_SCOPY, SGBFA, SGEFA,
!                              SSCAL
! Function routines called by SVJAC.. SVNORM
!-----------------------------------------------------------------------
! SVJAC is called by SVSTEP to compute and process the matrix
! P = I - h*rl1*J , where J is an approximation to the Jacobian.
! Here J is computed by the user-supplied routine JAC if
! MITER = 1 or 4, or by finite differencing if MITER = 2, 3, or 5.
! If MITER = 3, a diagonal approximation to J is used.
! If JSV = -1, J is computed from scratch in all cases.
! If JSV = 1 and MITER = 1, 2, 4, or 5, and if the saved value of J is
! considered acceptable, then P is constructed from the saved J.
! J is stored in wm and replaced by P.  If MITER .ne. 3, P is then
! subjected to LU decomposition in preparation for later solution
! of linear systems with P as coefficient matrix. This is done
! by SGEFA if MITER = 1 or 2, and by SGBFA if MITER = 4 or 5.
!
! Communication with SVJAC is done with the following variables.  (For
! more details, please see the comments in the driver subroutine.)
! Y          = Vector containing predicted values on entry.
! YH         = The Nordsieck array, an LDYH by LMAX array, input.
! LDYH       = A constant .ge. N, the first dimension of YH, input.
! EWT        = An error weight vector of length N.
! SAVF       = Array containing f evaluated at predicted y, input.
! WM         = Real work space for matrices.  In the output, it containS
!              the inverse diagonal matrix if MITER = 3 and the LU
!              decomposition of P if MITER is 1, 2 , 4, or 5.
!              Storage of matrix elements starts at WM(3).
!              Storage of the saved Jacobian starts at WM(LOCJS).
!              WM also contains the following matrix-related data..
!              WM(1) = SQRT(UROUND), used in numerical Jacobian step.
!              WM(2) = H*RL1, saved for later use if MITER = 3.
! IWM        = Integer work space containing pivot information,
!              starting at IWM(31), if MITER is 1, 2, 4, or 5.
!              IWM also contains band parameters ML = IWM(1) and
!              MU = IWM(2) if MITER is 4 or 5.
! F          = Dummy name for the user supplied subroutine for f.
! JAC        = Dummy name for the user supplied Jacobian subroutine.
! RPAR, IPAR = Dummy names for user's real and integer work arrays.
! RL1        = 1/EL(2) (input).
! IERPJ      = Output error flag,  = 0 if no trouble, 1 if the P
!              matrix is found to be singular.
! JCUR       = Output flag to indicate whether the Jacobian matrix
!              (or approximation) is now current.
!              JCUR = 0 means J is not current.
!              JCUR = 1 means J is current.
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
      REAL CON, DI, FAC, HRL1, ONE, PT1, R, R0, SRUR, THOU,             &
     &     YI, YJ, YJJ, ZERO
      INTEGER I, I1, I2, IER, II, J, J1, JJ, JOK, LENP, MBA, MBAND,     &
     &        MEB1, MEBAND, ML, ML3, MU, NP1
!
! Type declaration for function subroutines called ---------------------
!
      REAL SVNORM
!-----------------------------------------------------------------------
! The following Fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to this subroutine.
!-----------------------------------------------------------------------
      SAVE ONE, PT1, THOU, ZERO
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
      DATA ONE /1.0E0/, THOU /1000.0E0/, ZERO /0.0E0/, PT1 /0.1E0/
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('SVJAC',0,ZHOOK_HANDLE)
      IERPJ = 0
      HRL1 = H*RL1
! See whether J should be evaluated (JOK = -1) or not (JOK = 1). -------
      JOK = JSV
      IF (JSV .EQ. 1) THEN
        IF (NST .EQ. 0 .OR. NST .GT. NSLJ+MSBJ) JOK = -1
        IF (ICF .EQ. 1 .AND. DRC .LT. CCMXJ) JOK = -1
        IF (ICF .EQ. 2) JOK = -1
      ENDIF
! End of setting JOK. --------------------------------------------------
!
      IF (JOK .EQ. -1 .AND. MITER .EQ. 1) THEN
! If JOK = -1 and MITER = 1, call JAC to evaluate Jacobian. ------------
      NJE = NJE + 1
      NSLJ = NST
      JCUR = 1
      LENP = N*N
      DO 110 I = 1,LENP
 110    WM(I+2) = ZERO
!
!*UPG*MNH
!
      CALL JAC (N, TN, Y, 0, 0, WM(3), N, RPAR, IPAR, KMI, KINDEX)
!
!*UPG*MNH
!
      IF (JSV .EQ. 1) CALL CH_SCOPY (LENP, WM(3), 1, WM(LOCJS), 1)
      ENDIF
!
      IF (JOK .EQ. -1 .AND. MITER .EQ. 2) THEN
! If MITER = 2, make N calls to F to approximate the Jacobian. ---------
      NJE = NJE + 1
      NSLJ = NST
      JCUR = 1
      FAC = SVNORM (N, SAVF, EWT)
      R0 = THOU*ABS(H)*UROUND*REAL(N)*FAC
      IF (R0 .EQ. ZERO) R0 = ONE
      SRUR = WM(1)
      J1 = 2
      DO 230 J = 1,N
        YJ = Y(J)
        R = MAX(SRUR*ABS(YJ),R0/EWT(J))
        Y(J) = Y(J) + R
        FAC = ONE/R
!
!*UPG*MNH
!
        CALL F (N, TN, Y, FTEM, RPAR, IPAR, KMI, KINDEX)
!
!*UPG*MNH
!
        DO 220 I = 1,N
 220      WM(I+J1) = (FTEM(I) - SAVF(I))*FAC
        Y(J) = YJ
        J1 = J1 + N
 230    CONTINUE
      NFE = NFE + N
      LENP = N*N
      IF (JSV .EQ. 1) CALL CH_SCOPY (LENP, WM(3), 1, WM(LOCJS), 1)
      ENDIF
!
      IF (JOK .EQ. 1 .AND. (MITER .EQ. 1 .OR. MITER .EQ. 2)) THEN
      JCUR = 0
      LENP = N*N
      CALL CH_SCOPY (LENP, WM(LOCJS), 1, WM(3), 1)
      ENDIF
!
      IF (MITER .EQ. 1 .OR. MITER .EQ. 2) THEN
! Multiply Jacobian by scalar, add identity, and do LU decomposition. --
      CON = -HRL1
      CALL SSCAL (LENP, CON, WM(3), 1)
      J = 3
      NP1 = N + 1
      DO 250 I = 1,N
        WM(J) = WM(J) + ONE
 250    J = J + NP1
      NLU = NLU + 1
      CALL SGEFA (WM(3), N, N, IWM(31), IER)
      IF (IER .NE. 0) IERPJ = 1
      IF (LHOOK) CALL DR_HOOK('SVJAC',1,ZHOOK_HANDLE)
      RETURN
      ENDIF
! End of code block for MITER = 1 or 2. --------------------------------
!
      IF (MITER .EQ. 3) THEN
! If MITER = 3, construct a diagonal approximation to J and P. ---------
      NJE = NJE + 1
      JCUR = 1
      WM(2) = HRL1
      R = RL1*PT1
      DO 310 I = 1,N
 310    Y(I) = Y(I) + R*(H*SAVF(I) - YH(I,2))
!
!*UPG*MNH
!
      CALL F (N, TN, Y, WM(3), RPAR, IPAR, KMI, KINDEX)
!
!*UPG*MNH
!
      NFE = NFE + 1
      DO 320 I = 1,N
        R0 = H*SAVF(I) - YH(I,2)
        DI = PT1*R0 - H*(WM(I+2) - SAVF(I))
        WM(I+2) = ONE
        IF (ABS(R0) .LT. UROUND/EWT(I)) GO TO 320
        IF (ABS(DI) .EQ. ZERO) GO TO 330
        WM(I+2) = PT1*R0/DI
 320    CONTINUE
      IF (LHOOK) CALL DR_HOOK('SVJAC',1,ZHOOK_HANDLE)
      RETURN
 330  IERPJ = 1
      IF (LHOOK) CALL DR_HOOK('SVJAC',1,ZHOOK_HANDLE)
      RETURN
      ENDIF
! End of code block for MITER = 3. -------------------------------------
!
! Set constants for MITER = 4 or 5. ------------------------------------
      ML = IWM(1)
      MU = IWM(2)
      ML3 = ML + 3
      MBAND = ML + MU + 1
      MEBAND = MBAND + ML
      LENP = MEBAND*N
!
      IF (JOK .EQ. -1 .AND. MITER .EQ. 4) THEN
! If JOK = -1 and MITER = 4, call JAC to evaluate Jacobian. ------------
      NJE = NJE + 1
      NSLJ = NST
      JCUR = 1
      DO 410 I = 1,LENP
 410    WM(I+2) = ZERO
!
!*UPG*MNH
!
      CALL JAC (N, TN, Y, ML, MU, WM(ML3), MEBAND, RPAR, IPAR,          &
     &          KMI, KINDEX)
!
!*UPG*MNH
!
      IF (JSV .EQ. 1)                                                   &
     &   CALL SACOPY (MBAND, N, WM(ML3), MEBAND, WM(LOCJS), MBAND)
      ENDIF
!
      IF (JOK .EQ. -1 .AND. MITER .EQ. 5) THEN
! If MITER = 5, make N calls to F to approximate the Jacobian. ---------
      NJE = NJE + 1
      NSLJ = NST
      JCUR = 1
      MBA = MIN(MBAND,N)
      MEB1 = MEBAND - 1
      SRUR = WM(1)
      FAC = SVNORM (N, SAVF, EWT)
      R0 = THOU*ABS(H)*UROUND*REAL(N)*FAC
      IF (R0 .EQ. ZERO) R0 = ONE
      DO 560 J = 1,MBA
        DO 530 I = J,N,MBAND
          YI = Y(I)
          R = MAX(SRUR*ABS(YI),R0/EWT(I))
 530      Y(I) = Y(I) + R
!
!*UPG*MNH
!
        CALL F (N, TN, Y, FTEM, RPAR, IPAR, KMI, KINDEX)
!
!*UPG*MNH
!
        DO 550 JJ = J,N,MBAND
          Y(JJ) = YH(JJ,1)
          YJJ = Y(JJ)
          R = MAX(SRUR*ABS(YJJ),R0/EWT(JJ))
          FAC = ONE/R
          I1 = MAX(JJ-MU,1)
          I2 = MIN(JJ+ML,N)
          II = JJ*MEB1 - ML + 2
          DO 540 I = I1,I2
 540        WM(II+I) = (FTEM(I) - SAVF(I))*FAC
 550      CONTINUE
 560    CONTINUE
      NFE = NFE + MBA
      IF (JSV .EQ. 1)                                                   &
     &   CALL SACOPY (MBAND, N, WM(ML3), MEBAND, WM(LOCJS), MBAND)
      ENDIF
!
      IF (JOK .EQ. 1) THEN
      JCUR = 0
      CALL SACOPY (MBAND, N, WM(LOCJS), MBAND, WM(ML3), MEBAND)
      ENDIF
!
! Multiply Jacobian by scalar, add identity, and do LU decomposition.
      CON = -HRL1
      CALL SSCAL (LENP, CON, WM(3), 1 )
      II = MBAND + 2
      DO 580 I = 1,N
        WM(II) = WM(II) + ONE
 580    II = II + MEBAND
      NLU = NLU + 1
      CALL SGBFA (WM(3), MEBAND, N, ML, MU, IWM(31), IER)
      IF (IER .NE. 0) IERPJ = 1
      IF (LHOOK) CALL DR_HOOK('SVJAC',1,ZHOOK_HANDLE)
      RETURN
! End of code block for MITER = 4 or 5. --------------------------------
!
      END SUBROUTINE SVJAC
