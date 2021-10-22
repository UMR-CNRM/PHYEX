!     ######spl
      SUBROUTINE SVINDY (T, K, YH, LDYH, DKY, IFLAG)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!     ##############################################
      REAL T, YH, DKY
      INTEGER K, LDYH, IFLAG
      DIMENSION YH(LDYH,*), DKY(*)
!-----------------------------------------------------------------------
! Call sequence input -- T, K, YH, LDYH
! Call sequence output -- DKY, IFLAG
! COMMON block variables accessed..
!     /SVOD01/ --  H, TN, UROUND, L, N, NQ
!     /SVOD02/ --  HU
!
! Subroutines called by SVINDY.. SSCAL, XERRWV
! Function routines called by SVINDY.. None
!-----------------------------------------------------------------------
! SVINDY computes interpolated values of the K-th derivative of the
! dependent variable vector y, and stores it in DKY.  This routine
! is called within the package with K = 0 and T = TOUT, but may
! also be called by the user for any K up to the current order.
! (See detailed instructions in the usage documentation.)
!-----------------------------------------------------------------------
! The computed values in DKY are gotten by interpolation using the
! Nordsieck history array YH.  This array corresponds uniquely to a
! vector-valued polynomial of degree NQCUR or less, and DKY is set
! to the K-th derivative of this polynomial at T.
! The formula for DKY is..
!              q
!  DKY(i)  =  sum  c(j,K) * (T - TN)**(j-K) * H**(-j) * YH(i,j+1)
!             j=K
! where  c(j,K) = j*(j-1)*...*(j-K+1), q = NQCUR, TN = TCUR, H = HCUR.
! The quantities  NQ = NQCUR, L = NQ+1, N, TN, and H are
! communicated by COMMON.  The above sum is done in reverse order.
! IFLAG is returned negative if either K or T is out of bounds.
!
! Discussion above and comments in driver explain all variables.
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
      REAL C, HUN, R, S, TFUZZ, TN1, TP, ZERO
      INTEGER I, IC, J, JB, JB2, JJ, JJ1, JP1
      CHARACTER*80 MSG
!-----------------------------------------------------------------------
! The following Fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to this integrator.
!-----------------------------------------------------------------------
      SAVE HUN, ZERO
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
      DATA HUN /100.0E0/, ZERO /0.0E0/
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('SVINDY',0,ZHOOK_HANDLE)
      IFLAG = 0
      IF (K .LT. 0 .OR. K .GT. NQ) GO TO 80
      TFUZZ = HUN*UROUND*(TN + HU)
      TP = TN - HU - TFUZZ
      TN1 = TN + TFUZZ
      IF ((T-TP)*(T-TN1) .GT. ZERO) GO TO 90
!
      S = (T - TN)/H
      IC = 1
      IF (K .EQ. 0) GO TO 15
      JJ1 = L - K
      DO 10 JJ = JJ1, NQ
 10     IC = IC*JJ
 15   C = REAL(IC)
      DO 20 I = 1, N
 20     DKY(I) = C*YH(I,L)
      IF (K .EQ. NQ) GO TO 55
      JB2 = NQ - K
      DO 50 JB = 1, JB2
        J = NQ - JB
        JP1 = J + 1
        IC = 1
        IF (K .EQ. 0) GO TO 35
        JJ1 = JP1 - K
        DO 30 JJ = JJ1, J
 30       IC = IC*JJ
 35     C = REAL(IC)
        DO 40 I = 1, N
 40       DKY(I) = C*YH(I,JP1) + S*DKY(I)
 50     CONTINUE
      IF (K .EQ. 0)  THEN
        IF (LHOOK) CALL DR_HOOK('SVINDY',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
 55   R = H**(-K)
      CALL SSCAL (N, R, DKY, 1)
      IF (LHOOK) CALL DR_HOOK('SVINDY',1,ZHOOK_HANDLE)
      RETURN
!
 80   MSG = 'SVINDY-- K (=I1) illegal      '
      CALL XERRWV (MSG, 30, 51, 1, 1, K, 0, 0, ZERO, ZERO)
      IFLAG = -1
      IF (LHOOK) CALL DR_HOOK('SVINDY',1,ZHOOK_HANDLE)
      RETURN
 90   MSG = 'SVINDY-- T (=R1) illegal      '
      CALL XERRWV (MSG, 30, 52, 1, 0, 0, 0, 1, T, ZERO)
      MSG='      T not in interval TCUR - HU (= R1) to TCUR (=R2)      '
      CALL XERRWV (MSG, 60, 52, 1, 0, 0, 0, 2, TP, TN)
      IFLAG = -2
      IF (LHOOK) CALL DR_HOOK('SVINDY',1,ZHOOK_HANDLE)
      RETURN
      END SUBROUTINE SVINDY
