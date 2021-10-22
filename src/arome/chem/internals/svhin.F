!     ######spl
      SUBROUTINE SVHIN (N, T0, Y0, YDOT, F, RPAR, IPAR, TOUT, UROUND,   &
     &   EWT, ITOL, ATOL, Y, TEMP, H0, NITER, IER, KMI, KINDEX)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!     ###############################################################
      EXTERNAL F
      REAL T0, Y0, YDOT, RPAR, TOUT, UROUND, EWT, ATOL, Y,              &
     &   TEMP, H0
      INTEGER N, IPAR, ITOL, NITER, IER
      DIMENSION Y0(*), YDOT(*), EWT(*), ATOL(*), Y(*),                  &
     &   TEMP(*), RPAR(*), IPAR(*)
      INTEGER KMI, KINDEX
!-----------------------------------------------------------------------
! Call sequence input -- N, T0, Y0, YDOT, F, RPAR, IPAR, TOUT, UROUND,
!                        EWT, ITOL, ATOL, Y, TEMP
! Call sequence output -- H0, NITER, IER
! COMMON block variables accessed -- None
!
! Subroutines called by SVHIN.. F
! Function routines called by SVHIN.. SVNORM
!-----------------------------------------------------------------------
! This routine computes the step size, H0, to be attempted on the
! first step, when the user has not supplied a value for this.
!
! First we check that TOUT - T0 differs significantly from zero.  Then
! an iteration is done to approximate the initial second derivative
! and this is used to define h from w.r.m.s.norm(h**2 * yddot / 2) = 1.
! A bias factor of 1/2 is applied to the resulting h.
! The sign of H0 is inferred from the initial values of TOUT and T0.
!
! Communication with SVHIN is done with the following variables..
!
! N      = Size of ODE system, input.
! T0     = Initial value of independent variable, input.
! Y0     = Vector of initial conditions, input.
! YDOT   = Vector of initial first derivatives, input.
! F      = Name of subroutine for right-hand side f(t,y), input.
! RPAR, IPAR = Dummy names for user's real and integer work arrays.
! TOUT   = First output value of independent variable
! UROUND = Machine unit roundoff
! EWT, ITOL, ATOL = Error weights and tolerance parameters
!                   as described in the driver routine, input.
! Y, TEMP = Work arrays of length N.
! H0     = Step size to be attempted, output.
! NITER  = Number of iterations (and of f evaluations) to compute H0,
!          output.
! IER    = The error flag, returned with the value
!          IER = 0  if no trouble occurred, or
!          IER = -1 if TOUT and T0 are considered too close to proceed.
!-----------------------------------------------------------------------
!
! Type declarations for local variables --------------------------------
!
      REAL AFI, ATOLI, DELYI, HALF, HG, HLB, HNEW, HRAT,                &
     &     HUB, HUN, PT1, T1, TDIST, TROUND, TWO, YDDNRM
      INTEGER I, ITER
!
! Type declaration for function subroutines called ---------------------
!
      REAL SVNORM
!-----------------------------------------------------------------------
! The following Fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to this integrator.
!-----------------------------------------------------------------------
      SAVE HALF, HUN, PT1, TWO
      DATA HALF /0.5E0/, HUN /100.0E0/, PT1 /0.1E0/, TWO /2.0E0/
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('SVHIN',0,ZHOOK_HANDLE)
      NITER = 0
      TDIST = ABS(TOUT - T0)
      TROUND = UROUND*MAX(ABS(T0),ABS(TOUT))
      IF (TDIST .LT. TWO*TROUND) GO TO 100
!
! Set a lower bound on h based on the roundoff level in T0 and TOUT. ---
      HLB = HUN*TROUND
! Set an upper bound on h based on TOUT-T0 and the initial Y and YDOT. -
      HUB = PT1*TDIST
      ATOLI = ATOL(1)
      DO 10 I = 1, N
        IF (ITOL .EQ. 2 .OR. ITOL .EQ. 4) ATOLI = ATOL(I)
        DELYI = PT1*ABS(Y0(I)) + ATOLI
        AFI = ABS(YDOT(I))
        IF (AFI*HUB .GT. DELYI) HUB = DELYI/AFI
 10     CONTINUE
!
! Set initial guess for h as geometric mean of upper and lower bounds. -
      ITER = 0
      HG = SQRT(HLB*HUB)
! If the bounds have crossed, exit with the mean value. ----------------
      IF (HUB .LT. HLB) THEN
        H0 = HG
        GO TO 90
      ENDIF
!
! Looping point for iteration. -----------------------------------------
 50   CONTINUE
! Estimate the second derivative as a difference quotient in f. --------
      T1 = T0 + HG
      DO 60 I = 1, N
 60     Y(I) = Y0(I) + HG*YDOT(I)
!
!*UPG*MNH
!
      CALL F (N, T1, Y, TEMP, RPAR, IPAR, KMI, KINDEX)
!
!*UPG*MNH
!
      DO 70 I = 1, N
 70     TEMP(I) = (TEMP(I) - YDOT(I))/HG
      YDDNRM = SVNORM (N, TEMP, EWT)
! Get the corresponding new value of h. --------------------------------
      IF (YDDNRM*HUB*HUB .GT. TWO) THEN
        HNEW = SQRT(TWO/YDDNRM)
      ELSE
        HNEW = SQRT(HG*HUB)
      ENDIF
      ITER = ITER + 1
!-----------------------------------------------------------------------
! Test the stopping conditions.
! Stop if the new and previous h values differ by a factor of .lt. 2.
! Stop if four iterations have been done.  Also, stop with previous h
! if HNEW/HG .gt. 2 after first iteration, as this probably means that
! the second derivative value is bad because of cancellation error.
!-----------------------------------------------------------------------
      IF (ITER .GE. 4) GO TO 80
      HRAT = HNEW/HG
      IF ( (HRAT .GT. HALF) .AND. (HRAT .LT. TWO) ) GO TO 80
      IF ( (ITER .GE. 2) .AND. (HNEW .GT. TWO*HG) ) THEN
        HNEW = HG
        GO TO 80
      ENDIF
      HG = HNEW
      GO TO 50
!
! Iteration done.  Apply bounds, bias factor, and sign.  Then exit. ----
 80   H0 = HNEW*HALF
      IF (H0 .LT. HLB) H0 = HLB
      IF (H0 .GT. HUB) H0 = HUB
 90   H0 = SIGN(H0, TOUT - T0)
      NITER = ITER
      IER = 0
      IF (LHOOK) CALL DR_HOOK('SVHIN',1,ZHOOK_HANDLE)
      RETURN
! Error return for TOUT - T0 too small. --------------------------------
 100  IER = -1
      IF (LHOOK) CALL DR_HOOK('SVHIN',1,ZHOOK_HANDLE)
      RETURN
      END SUBROUTINE SVHIN
