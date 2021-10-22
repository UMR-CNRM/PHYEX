!     ######spl
      SUBROUTINE SVODE (F, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK,    &
     &            ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF,        &
     &            RPAR, IPAR, KMI, KINDEX)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!
!
      EXTERNAL F, JAC
      !
!
!*UPG*MNH
!
      INTEGER KMI, KINDEX
!
!*UPG*MNH
!
      REAL Y, T, TOUT, RTOL, ATOL, RWORK, RPAR
      INTEGER NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, IWORK, LIW,          &
     &        MF, IPAR
      DIMENSION Y(*), RTOL(*), ATOL(*), RWORK(LRW), IWORK(LIW),         &
     &          RPAR(*), IPAR(*)
!-----------------------------------------------------------------------
! SVODE.. Variable-coefficient Ordinary Differential Equation solver,
! with fixed-leading coefficient implementation.
! This version is in single precision.
!
! SVODE solves the initial value problem for stiff or nonstiff
! systems of first order ODEs,
!     dy/dt = f(t,y) ,  or, in component form,
!     dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(NEQ)) (i = 1,...,NEQ).
! SVODE is a package based on the EPISODE and EPISODEB packages, and
! on the ODEPACK user interface standard, with minor modifications.
!-----------------------------------------------------------------------
! Revision History (YYMMDD)
!   890615  Date Written
!   890922  Added interrupt/restart ability, minor changes throughout.
!   910228  Minor revisions in line format,  prologue, etc.
!   920227  Modifications by D. Pang:
!           (1) Applied subgennam to get generic intrinsic names.
!           (2) Changed intrinsic names to generic in comments.
!           (3) Added *DECK lines before each routine.
!   920721  Names of routines and labeled Common blocks changed, so as
!           to be unique in combined single/double precision code (ACH).
!   920722  Minor revisions to prologue (ACH).
!-----------------------------------------------------------------------
! References..
!
! 1. P. N. Brown, G. D. Byrne, and A. C. Hindmarsh, "VODE: A Variable
!    Coefficient ODE Solver," SIAM J. Sci. Stat. Comput., 10 (1989),
!    pp. 1038-1051.  Also, LLNL Report UCRL-98412, June 1988.
! 2. G. D. Byrne and A. C. Hindmarsh, "A Polyalgorithm for the
!    Numerical Solution of Ordinary Differential Equations,"
!    ACM Trans. Math. Software, 1 (1975), pp. 71-96.
! 3. A. C. Hindmarsh and G. D. Byrne, "EPISODE: An Effective Package
!    for the Integration of Systems of Ordinary Differential
!    Equations," LLNL Report UCID-30112, Rev. 1, April 1977.
! 4. G. D. Byrne and A. C. Hindmarsh, "EPISODEB: An Experimental
!    Package for the Integration of Systems of Ordinary Differential
!    Equations with Banded Jacobians," LLNL Report UCID-30132, April
!    1976.
! 5. A. C. Hindmarsh, "ODEPACK, a Systematized Collection of ODE
!    Solvers," in Scientific Computing, R. S. Stepleman et al., eds.,
!    North-Holland, Amsterdam, 1983, pp. 55-64.
! 6. K. R. Jackson and R. Sacks-Davis, "An Alternative Implementation
!    of Variable Step-Size Multistep Formulas for Stiff ODEs," ACM
!    Trans. Math. Software, 6 (1980), pp. 295-318.
!-----------------------------------------------------------------------
! Authors..
!
!               Peter N. Brown and Alan C. Hindmarsh
!               Computing and Mathematics Research Division, L-316
!               Lawrence Livermore National Laboratory
!               Livermore, CA 94550
! and
!               George D. Byrne
!               Exxon Research and Engineering Co.
!               Clinton Township
!               Route 22 East
!               Annandale, NJ 08801
!-----------------------------------------------------------------------
! Summary of usage.
!
! Communication between the user and the SVODE package, for normal
! situations, is summarized here.  This summary describes only a subset
! of the full set of options available.  See the full description for
! details, including optional communication, nonstandard options,
! and instructions for special situations.  See also the example
! problem (with program and output) following this summary.
!
! A. First provide a subroutine of the form..
!
!           SUBROUTINE F (NEQ, T, Y, YDOT, RPAR, IPAR)
!           REAL T, Y, YDOT, RPAR
!           DIMENSION Y(NEQ), YDOT(NEQ)
!
! which supplies the vector function f by loading YDOT(i) with f(i).
!
! B. Next determine (or guess) whether or not the problem is stiff.
! Stiffness occurs when the Jacobian matrix df/dy has an eigenvalue
! whose real part is negative and large in magnitude, compared to the
! reciprocal of the t span of interest.  If the problem is nonstiff,
! use a method flag MF = 10.  If it is stiff, there are four standard
! choices for MF (21, 22, 24, 25), and SVODE requires the Jacobian
! matrix in some form.  In these cases (MF .gt. 0), SVODE will use a
! saved copy of the Jacobian matrix.  If this is undesirable because of
! storage limitations, set MF to the corresponding negative value
! (-21, -22, -24, -25).  (See full description of MF below.)
! The Jacobian matrix is regarded either as full (MF = 21 or 22),
! or banded (MF = 24 or 25).  In the banded case, SVODE requires two
! half-bandwidth parameters ML and MU.  These are, respectively, the
! widths of the lower and upper parts of the band, excluding the main
! diagonal.  Thus the band consists of the locations (i,j) with
! i-ML .le. j .le. i+MU, and the full bandwidth is ML+MU+1.
!
! C. If the problem is stiff, you are encouraged to supply the Jacobian
! directly (MF = 21 or 24), but if this is not feasible, SVODE will
! compute it internally by difference quotients (MF = 22 or 25).
! If you are supplying the Jacobian, provide a subroutine of the form..
!
!           SUBROUTINE JAC (NEQ, T, Y, ML, MU, PD, NROWPD, RPAR, IPAR)
!           REAL T, Y, PD, RPAR
!           DIMENSION Y(NEQ), PD(NROWPD,NEQ)
!
! which supplies df/dy by loading PD as follows..
!     For a full Jacobian (MF = 21), load PD(i,j) with df(i)/dy(j),
! the partial derivative of f(i) with respect to y(j).  (Ignore the
! ML and MU arguments in this case.)
!     For a banded Jacobian (MF = 24), load PD(i-j+MU+1,j) with
! df(i)/dy(j), i.e. load the diagonal lines of df/dy into the rows of
! PD from the top down.
!     In either case, only nonzero elements need be loaded.
!
! D. Write a main program which calls subroutine SVODE once for
! each point at which answers are desired.  This should also provide
! for possible use of logical unit 6 for output of error messages
! by SVODE.  On the first call to SVODE, supply arguments as follows..
! F      = Name of subroutine for right-hand side vector f.
!          This name must be declared external in calling program.
! NEQ    = Number of first order ODE-s.
! Y      = Array of initial values, of length NEQ.
! T      = The initial value of the independent variable.
! TOUT   = First point where output is desired (.ne. T).
! ITOL   = 1 or 2 according as ATOL (below) is a scalar or array.
! RTOL   = Relative tolerance parameter (scalar).
! ATOL   = Absolute tolerance parameter (scalar or array).
!          The estimated local error in Y(i) will be controlled so as
!          to be roughly less (in magnitude) than
!             EWT(i) = RTOL*abs(Y(i)) + ATOL     if ITOL = 1, or
!             EWT(i) = RTOL*abs(Y(i)) + ATOL(i)  if ITOL = 2.
!          Thus the local error test passes if, in each component,
!          either the absolute error is less than ATOL (or ATOL(i)),
!          or the relative error is less than RTOL.
!          Use RTOL = 0.0 for pure absolute error control, and
!          use ATOL = 0.0 (or ATOL(i) = 0.0) for pure relative error
!          control.  Caution.. Actual (global) errors may exceed these
!          local tolerances, so choose them conservatively.
! ITASK  = 1 for normal computation of output values of Y at t = TOUT.
! ISTATE = Integer flag (input and output).  Set ISTATE = 1.
! IOPT   = 0 to indicate no optional input used.
! RWORK  = Real work array of length at least..
!             20 + 16*NEQ                      for MF = 10,
!             22 +  9*NEQ + 2*NEQ**2           for MF = 21 or 22,
!             22 + 11*NEQ + (3*ML + 2*MU)*NEQ  for MF = 24 or 25.
! LRW    = Declared length of RWORK (in user's DIMENSION statement).
! IWORK  = Integer work array of length at least..
!             30        for MF = 10,
!             30 + NEQ  for MF = 21, 22, 24, or 25.
!          If MF = 24 or 25, input in IWORK(1),IWORK(2) the lower
!          and upper half-bandwidths ML,MU.
! LIW    = Declared length of IWORK (in user's DIMENSION).
! JAC    = Name of subroutine for Jacobian matrix (MF = 21 or 24).
!          If used, this name must be declared external in calling
!          program.  If not used, pass a dummy name.
! MF     = Method flag.  Standard values are..
!          10 for nonstiff (Adams) method, no Jacobian used.
!          21 for stiff (BDF) method, user-supplied full Jacobian.
!          22 for stiff method, internally generated full Jacobian.
!          24 for stiff method, user-supplied banded Jacobian.
!          25 for stiff method, internally generated banded Jacobian.
! RPAR,IPAR = user-defined real and integer arrays passed to F and JAC.
! Note that the main program must declare arrays Y, RWORK, IWORK,
! and possibly ATOL, RPAR, and IPAR.
!
! E. The output from the first call (or any call) is..
!      Y = Array of computed values of y(t) vector.
!      T = Corresponding value of independent variable (normally TOUT).
! ISTATE = 2  if SVODE was successful, negative otherwise.
!          -1 means excess work done on this call. (Perhaps wrong MF.)
!          -2 means excess accuracy requested. (Tolerances too small.)
!          -3 means illegal input detected. (See printed message.)
!          -4 means repeated error test failures. (Check all input.)
!          -5 means repeated convergence failures. (Perhaps bad
!             Jacobian supplied or wrong choice of MF or tolerances.)
!          -6 means error weight became zero during problem. (Solution
!             component i vanished, and ATOL or ATOL(i) = 0.)
!
! F. To continue the integration after a successful return, simply
! reset TOUT and call SVODE again.  No other parameters need be reset.
!
!-----------------------------------------------------------------------
! EXAMPLE PROBLEM
!
! The following is a simple example problem, with the coding
! needed for its solution by SVODE.  The problem is from chemical
! kinetics, and consists of the following three rate equations..
!     dy1/dt = -.04*y1 + 1.e4*y2*y3
!     dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e7*y2**2
!     dy3/dt = 3.e7*y2**2
! on the interval from t = 0.0 to t = 4.e10, with initial conditions
! y1 = 1.0, y2 = y3 = 0.  The problem is stiff.
!
! The following coding solves this problem with SVODE, using MF = 21
! and printing results at t = .4, 4., ..., 4.e10.  It uses
! ITOL = 2 and ATOL much smaller for y2 than y1 or y3 because
! y2 has much smaller values.
! At the end of the run, statistical quantities of interest are
! printed. (See optional output in the full description below.)
! To generate Fortran source code, replace C in column 1 with a blank
! in the coding below.
!
!     EXTERNAL FEX, JEX
!     REAL ATOL, RPAR, RTOL, RWORK, T, TOUT, Y
!     DIMENSION Y(3), ATOL(3), RWORK(67), IWORK(33)
!     NEQ = 3
!     Y(1) = 1.0E0
!     Y(2) = 0.0E0
!     Y(3) = 0.0E0
!     T = 0.0E0
!     TOUT = 0.4E0
!     ITOL = 2
!     RTOL = 1.E-4
!     ATOL(1) = 1.E-8
!     ATOL(2) = 1.E-14
!     ATOL(3) = 1.E-6
!     ITASK = 1
!     ISTATE = 1
!     IOPT = 0
!     LRW = 67
!     LIW = 33
!     MF = 21
!     DO 40 IOUT = 1,12
!       CALL SVODE(FEX,NEQ,Y,T,TOUT,ITOL,RTOL,ATOL,ITASK,ISTATE,
!    1            IOPT,RWORK,LRW,IWORK,LIW,JEX,MF,RPAR,IPAR)
!       WRITE(6,20)T,Y(1),Y(2),Y(3)
! 20    FORMAT(' At t =',E12.4,'   y =',3E14.6)
!       IF (ISTATE .LT. 0) GO TO 80
! 40    TOUT = TOUT*10.
!     WRITE(6,60) IWORK(11),IWORK(12),IWORK(13),IWORK(19),
!    1            IWORK(20),IWORK(21),IWORK(22)
! 60  FORMAT(/' No. steps =',I4,'   No. f-s =',I4,
!    1       '   No. J-s =',I4,'   No. LU-s =',I4/
!    2       '  No. nonlinear iterations =',I4/
!    3       '  No. nonlinear convergence failures =',I4/
!    4       '  No. error test failures =',I4/)
!     STOP
! 80  WRITE(6,90)ISTATE
! 90  FORMAT(///' Error halt.. ISTATE =',I3)
!     STOP
!     END
!
!     SUBROUTINE FEX (NEQ, T, Y, YDOT, RPAR, IPAR)
!     REAL RPAR, T, Y, YDOT
!     DIMENSION Y(NEQ), YDOT(NEQ)
!     YDOT(1) = -.04E0*Y(1) + 1.E4*Y(2)*Y(3)
!     YDOT(3) = 3.E7*Y(2)*Y(2)
!     YDOT(2) = -YDOT(1) - YDOT(3)
!     RETURN
!     END
!
!     SUBROUTINE JEX (NEQ, T, Y, ML, MU, PD, NRPD, RPAR, IPAR)
!     REAL PD, RPAR, T, Y
!     DIMENSION Y(NEQ), PD(NRPD,NEQ)
!     PD(1,1) = -.04E0
!     PD(1,2) = 1.E4*Y(3)
!     PD(1,3) = 1.E4*Y(2)
!     PD(2,1) = .04E0
!     PD(2,3) = -PD(1,3)
!     PD(3,2) = 6.E7*Y(2)
!     PD(2,2) = -PD(1,2) - PD(3,2)
!     RETURN
!     END
!
! The following output was obtained from the above program on a
! Cray-1 computer with the CFT compiler.
!
! At t =  4.0000e-01   y =  9.851680e-01  3.386314e-05  1.479817e-02
! At t =  4.0000e+00   y =  9.055255e-01  2.240539e-05  9.445214e-02
! At t =  4.0000e+01   y =  7.158108e-01  9.184883e-06  2.841800e-01
! At t =  4.0000e+02   y =  4.505032e-01  3.222940e-06  5.494936e-01
! At t =  4.0000e+03   y =  1.832053e-01  8.942690e-07  8.167938e-01
! At t =  4.0000e+04   y =  3.898560e-02  1.621875e-07  9.610142e-01
! At t =  4.0000e+05   y =  4.935882e-03  1.984013e-08  9.950641e-01
! At t =  4.0000e+06   y =  5.166183e-04  2.067528e-09  9.994834e-01
! At t =  4.0000e+07   y =  5.201214e-05  2.080593e-10  9.999480e-01
! At t =  4.0000e+08   y =  5.213149e-06  2.085271e-11  9.999948e-01
! At t =  4.0000e+09   y =  5.183495e-07  2.073399e-12  9.999995e-01
! At t =  4.0000e+10   y =  5.450996e-08  2.180399e-13  9.999999e-01
!
! No. steps = 595   No. f-s = 832   No. J-s =  13   No. LU-s = 112
!  No. nonlinear iterations = 831
!  No. nonlinear convergence failures =   0
!  No. error test failures =  22
!-----------------------------------------------------------------------
! Full description of user interface to SVODE.
!
! The user interface to SVODE consists of the following parts.
!
! i.   The call sequence to subroutine SVODE, which is a driver
!      routine for the solver.  This includes descriptions of both
!      the call sequence arguments and of user-supplied routines.
!      Following these descriptions is
!        * a description of optional input available through the
!          call sequence,
!        * a description of optional output (in the work arrays), and
!        * instructions for interrupting and restarting a solution.
!
! ii.  Descriptions of other routines in the SVODE package that may be
!      (optionally) called by the user.  These provide the ability to
!      alter error message handling, save and restore the internal
!      COMMON, and obtain specified derivatives of the solution y(t).
!
! iii. Descriptions of COMMON blocks to be declared in overlay
!      or similar environments.
!
! iv.  Description of two routines in the SVODE package, either of
!      which the user may replace with his own version, if desired.
!      these relate to the measurement of errors.
!
!-----------------------------------------------------------------------
! Part i.  Call Sequence.
!
! The call sequence parameters used for input only are
!     F, NEQ, TOUT, ITOL, RTOL, ATOL, ITASK, IOPT, LRW, LIW, JAC, MF,
! and those used for both input and output are
!     Y, T, ISTATE.
! The work arrays RWORK and IWORK are also used for conditional and
! optional input and optional output.  (The term output here refers
! to the return from subroutine SVODE to the user's calling program.)
!
! The legality of input parameters will be thoroughly checked on the
! initial call for the problem, but not checked thereafter unless a
! change in input parameters is flagged by ISTATE = 3 in the input.
!
! The descriptions of the call arguments are as follows.
!
! F      = The name of the user-supplied subroutine defining the
!          ODE system.  The system must be put in the first-order
!          form dy/dt = f(t,y), where f is a vector-valued function
!          of the scalar t and the vector y.  Subroutine F is to
!          compute the function f.  It is to have the form
!               SUBROUTINE F (NEQ, T, Y, YDOT, RPAR, IPAR)
!               REAL T, Y, YDOT, RPAR
!               DIMENSION Y(NEQ), YDOT(NEQ)
!          where NEQ, T, and Y are input, and the array YDOT = f(t,y)
!          is output.  Y and YDOT are arrays of length NEQ.
!          (In the DIMENSION statement above, NEQ  can be replaced by
!          *  to make  Y  and  YDOT  assumed size arrays.)
!          Subroutine F should not alter Y(1),...,Y(NEQ).
!          F must be declared EXTERNAL in the calling program.
!
!          Subroutine F may access user-defined real and integer
!          work arrays RPAR and IPAR, which are to be dimensioned
!          in the main program.
!
!          If quantities computed in the F routine are needed
!          externally to SVODE, an extra call to F should be made
!          for this purpose, for consistent and accurate results.
!          If only the derivative dy/dt is needed, use SVINDY instead.
!
! NEQ    = The size of the ODE system (number of first order
!          ordinary differential equations).  Used only for input.
!          NEQ may not be increased during the problem, but
!          can be decreased (with ISTATE = 3 in the input).
!
! Y      = A real array for the vector of dependent variables, of
!          length NEQ or more.  Used for both input and output on the
!          first call (ISTATE = 1), and only for output on other calls.
!          On the first call, Y must contain the vector of initial
!          values.  In the output, Y contains the computed solution
!          evaluated at T.  If desired, the Y array may be used
!          for other purposes between calls to the solver.
!
!          This array is passed as the Y argument in all calls to
!          F and JAC.
!
! T      = The independent variable.  In the input, T is used only on
!          the first call, as the initial point of the integration.
!          In the output, after each call, T is the value at which a
!          computed solution Y is evaluated (usually the same as TOUT).
!          On an error return, T is the farthest point reached.
!
! TOUT   = The next value of t at which a computed solution is desired.
!          Used only for input.
!
!          When starting the problem (ISTATE = 1), TOUT may be equal
!          to T for one call, then should .ne. T for the next call.
!          For the initial T, an input value of TOUT .ne. T is used
!          in order to determine the direction of the integration
!          (i.e. the algebraic sign of the step sizes) and the rough
!          scale of the problem.  Integration in either direction
!          (forward or backward in t) is permitted.
!
!          If ITASK = 2 or 5 (one-step modes), TOUT is ignored after
!          the first call (i.e. the first call with TOUT .ne. T).
!          Otherwise, TOUT is required on every call.
!
!          If ITASK = 1, 3, or 4, the values of TOUT need not be
!          monotone, but a value of TOUT which backs up is limited
!          to the current internal t interval, whose endpoints are
!          TCUR - HU and TCUR.  (See optional output, below, for
!          TCUR and HU.)
!
! ITOL   = An indicator for the type of error control.  See
!          description below under ATOL.  Used only for input.
!
! RTOL   = A relative error tolerance parameter, either a scalar or
!          an array of length NEQ.  See description below under ATOL.
!          Input only.
!
! ATOL   = An absolute error tolerance parameter, either a scalar or
!          an array of length NEQ.  Input only.
!
!          The input parameters ITOL, RTOL, and ATOL determine
!          the error control performed by the solver.  The solver will
!          control the vector e = (e(i)) of estimated local errors
!          in Y, according to an inequality of the form
!                      rms-norm of ( e(i)/EWT(i) )   .le.   1,
!          where       EWT(i) = RTOL(i)*abs(Y(i)) + ATOL(i),
!          and the rms-norm (root-mean-square norm) here is
!          rms-norm(v) = sqrt(sum v(i)**2 / NEQ).  Here EWT = (EWT(i))
!          is a vector of weights which must always be positive, and
!          the values of RTOL and ATOL should all be non-negative.
!          The following table gives the types (scalar/array) of
!          RTOL and ATOL, and the corresponding form of EWT(i).
!
!             ITOL    RTOL       ATOL          EWT(i)
!              1     scalar     scalar     RTOL*ABS(Y(i)) + ATOL
!              2     scalar     array      RTOL*ABS(Y(i)) + ATOL(i)
!              3     array      scalar     RTOL(i)*ABS(Y(i)) + ATOL
!              4     array      array      RTOL(i)*ABS(Y(i)) + ATOL(i)
!
!          When either of these parameters is a scalar, it need not
!          be dimensioned in the user's calling program.
!
!          If none of the above choices (with ITOL, RTOL, and ATOL
!          fixed throughout the problem) is suitable, more general
!          error controls can be obtained by substituting
!          user-supplied routines for the setting of EWT and/or for
!          the norm calculation.  See Part iv below.
!
!          If global errors are to be estimated by making a repeated
!          run on the same problem with smaller tolerances, then all
!          components of RTOL and ATOL (i.e. of EWT) should be scaled
!          down uniformly.
!
! ITASK  = An index specifying the task to be performed.
!          Input only.  ITASK has the following values and meanings.
!          1  means normal computation of output values of y(t) at
!             t = TOUT (by overshooting and interpolating).
!          2  means take one step only and return.
!          3  means stop at the first internal mesh point at or
!             beyond t = TOUT and return.
!          4  means normal computation of output values of y(t) at
!             t = TOUT but without overshooting t = TCRIT.
!             TCRIT must be input as RWORK(1).  TCRIT may be equal to
!             or beyond TOUT, but not behind it in the direction of
!             integration.  This option is useful if the problem
!             has a singularity at or beyond t = TCRIT.
!          5  means take one step, without passing TCRIT, and return.
!             TCRIT must be input as RWORK(1).
!
!          Note..  If ITASK = 4 or 5 and the solver reaches TCRIT
!          (within roundoff), it will return T = TCRIT (exactly) to
!          indicate this (unless ITASK = 4 and TOUT comes before TCRIT,
!          in which case answers at T = TOUT are returned first).
!
! ISTATE = an index used for input and output to specify the
!          the state of the calculation.
!
!          In the input, the values of ISTATE are as follows.
!          1  means this is the first call for the problem
!             (initializations will be done).  See note below.
!          2  means this is not the first call, and the calculation
!             is to continue normally, with no change in any input
!             parameters except possibly TOUT and ITASK.
!             (If ITOL, RTOL, and/or ATOL are changed between calls
!             with ISTATE = 2, the new values will be used but not
!             tested for legality.)
!          3  means this is not the first call, and the
!             calculation is to continue normally, but with
!             a change in input parameters other than
!             TOUT and ITASK.  Changes are allowed in
!             NEQ, ITOL, RTOL, ATOL, IOPT, LRW, LIW, MF, ML, MU,
!             and any of the optional input except H0.
!             (See IWORK description for ML and MU.)
!          Note..  A preliminary call with TOUT = T is not counted
!          as a first call here, as no initialization or checking of
!          input is done.  (Such a call is sometimes useful to include
!          the initial conditions in the output.)
!          Thus the first call for which TOUT .ne. T requires
!          ISTATE = 1 in the input.
!
!          In the output, ISTATE has the following values and meanings.
!           1  means nothing was done, as TOUT was equal to T with
!              ISTATE = 1 in the input.
!           2  means the integration was performed successfully.
!          -1  means an excessive amount of work (more than MXSTEP
!              steps) was done on this call, before completing the
!              requested task, but the integration was otherwise
!              successful as far as T.  (MXSTEP is an optional input
!              and is normally 500.)  To continue, the user may
!              simply reset ISTATE to a value .gt. 1 and call again.
!              (The excess work step counter will be reset to 0.)
!              In addition, the user may increase MXSTEP to avoid
!              this error return.  (See optional input below.)
!          -2  means too much accuracy was requested for the precision
!              of the machine being used.  This was detected before
!              completing the requested task, but the integration
!              was successful as far as T.  To continue, the tolerance
!              parameters must be reset, and ISTATE must be set
!              to 3.  The optional output TOLSF may be used for this
!              purpose.  (Note.. If this condition is detected before
!              taking any steps, then an illegal input return
!              (ISTATE = -3) occurs instead.)
!          -3  means illegal input was detected, before taking any
!              integration steps.  See written message for details.
!              Note..  If the solver detects an infinite loop of calls
!              to the solver with illegal input, it will cause
!              the run to stop.
!          -4  means there were repeated error test failures on
!              one attempted step, before completing the requested
!              task, but the integration was successful as far as T.
!              The problem may have a singularity, or the input
!              may be inappropriate.
!          -5  means there were repeated convergence test failures on
!              one attempted step, before completing the requested
!              task, but the integration was successful as far as T.
!              This may be caused by an inaccurate Jacobian matrix,
!              if one is being used.
!          -6  means EWT(i) became zero for some i during the
!              integration.  Pure relative error control (ATOL(i)=0.0)
!              was requested on a variable which has now vanished.
!              The integration was successful as far as T.
!
!          Note..  Since the normal output value of ISTATE is 2,
!          it does not need to be reset for normal continuation.
!          Also, since a negative input value of ISTATE will be
!          regarded as illegal, a negative output value requires the
!          user to change it, and possibly other input, before
!          calling the solver again.
!
! IOPT   = An integer flag to specify whether or not any optional
!          input is being used on this call.  Input only.
!          The optional input is listed separately below.
!          IOPT = 0 means no optional input is being used.
!                   Default values will be used in all cases.
!          IOPT = 1 means optional input is being used.
!
! RWORK  = A real working array (single precision).
!          The length of RWORK must be at least
!             20 + NYH*(MAXORD + 1) + 3*NEQ + LWM    where
!          NYH    = the initial value of NEQ,
!          MAXORD = 12 (if METH = 1) or 5 (if METH = 2) (unless a
!                   smaller value is given as an optional input),
!          LWM = length of work space for matrix-related data..
!          LWM = 0             if MITER = 0,
!          LWM = 2*NEQ**2 + 2  if MITER = 1 or 2, and MF.gt.0,
!          LWM = NEQ**2 + 2    if MITER = 1 or 2, and MF.lt.0,
!          LWM = NEQ + 2       if MITER = 3,
!          LWM = (3*ML+2*MU+2)*NEQ + 2 if MITER = 4 or 5, and MF.gt.0,
!          LWM = (2*ML+MU+1)*NEQ + 2   if MITER = 4 or 5, and MF.lt.0.
!          (See the MF description for METH and MITER.)
!          Thus if MAXORD has its default value and NEQ is constant,
!          this length is..
!             20 + 16*NEQ                    for MF = 10,
!             22 + 16*NEQ + 2*NEQ**2         for MF = 11 or 12,
!             22 + 16*NEQ + NEQ**2           for MF = -11 or -12,
!             22 + 17*NEQ                    for MF = 13,
!             22 + 18*NEQ + (3*ML+2*MU)*NEQ  for MF = 14 or 15,
!             22 + 17*NEQ + (2*ML+MU)*NEQ    for MF = -14 or -15,
!             20 +  9*NEQ                    for MF = 20,
!             22 +  9*NEQ + 2*NEQ**2         for MF = 21 or 22,
!             22 +  9*NEQ + NEQ**2           for MF = -21 or -22,
!             22 + 10*NEQ                    for MF = 23,
!             22 + 11*NEQ + (3*ML+2*MU)*NEQ  for MF = 24 or 25.
!             22 + 10*NEQ + (2*ML+MU)*NEQ    for MF = -24 or -25.
!          The first 20 words of RWORK are reserved for conditional
!          and optional output.
!
!          The following word in RWORK is a conditional input..
!            RWORK(1) = TCRIT = critical value of t which the solver
!                       is not to overshoot.  Required if ITASK is
!                       4 or 5, and ignored otherwise.  (See ITASK.)
!
! LRW    = The length of the array RWORK, as declared by the user.
!          (This will be checked by the solver.)
!
! IWORK  = An integer work array.  The length of IWORK must be at least
!             30        if MITER = 0 or 3 (MF = 10, 13, 20, 23), or
!             30 + NEQ  otherwise (abs(MF) = 11,12,14,15,21,22,24,25).
!          The first 30 words of IWORK are reserved for conditional and
!          optional input and optional output.
!
!          The following 2 words in IWORK are conditional input..
!            IWORK(1) = ML     These are the lower and upper
!            IWORK(2) = MU     half-bandwidths, respectively, of the
!                       banded Jacobian, excluding the main diagonal.
!                       The band is defined by the matrix locations
!                       (i,j) with i-ML .le. j .le. i+MU.  ML and MU
!                       must satisfy  0 .le.  ML,MU  .le. NEQ-1.
!                       These are required if MITER is 4 or 5, and
!                       ignored otherwise.  ML and MU may in fact be
!                       the band parameters for a matrix to which
!                       df/dy is only approximately equal.
!
! LIW    = the length of the array IWORK, as declared by the user.
!          (This will be checked by the solver.)
!
! Note..  The work arrays must not be altered between calls to SVODE
! for the same problem, except possibly for the conditional and
! optional input, and except for the last 3*NEQ words of RWORK.
! The latter space is used for internal scratch space, and so is
! available for use by the user outside SVODE between calls, if
! desired (but not for use by F or JAC).
!
! JAC    = The name of the user-supplied routine (MITER = 1 or 4) to
!          compute the Jacobian matrix, df/dy, as a function of
!          the scalar t and the vector y.  It is to have the form
!               SUBROUTINE JAC (NEQ, T, Y, ML, MU, PD, NROWPD,
!                               RPAR, IPAR)
!               REAL T, Y, PD, RPAR
!               DIMENSION Y(NEQ), PD(NROWPD, NEQ)
!          where NEQ, T, Y, ML, MU, and NROWPD are input and the array
!          PD is to be loaded with partial derivatives (elements of the
!          Jacobian matrix) in the output.  PD must be given a first
!          dimension of NROWPD.  T and Y have the same meaning as in
!          Subroutine F.  (In the DIMENSION statement above, NEQ can
!          be replaced by  *  to make Y and PD assumed size arrays.)
!               In the full matrix case (MITER = 1), ML and MU are
!          ignored, and the Jacobian is to be loaded into PD in
!          columnwise manner, with df(i)/dy(j) loaded into PD(i,j).
!               In the band matrix case (MITER = 4), the elements
!          within the band are to be loaded into PD in columnwise
!          manner, with diagonal lines of df/dy loaded into the rows
!          of PD. Thus df(i)/dy(j) is to be loaded into PD(i-j+MU+1,j).
!          ML and MU are the half-bandwidth parameters. (See IWORK).
!          The locations in PD in the two triangular areas which
!          correspond to nonexistent matrix elements can be ignored
!          or loaded arbitrarily, as they are overwritten by SVODE.
!               JAC need not provide df/dy exactly.  A crude
!          approximation (possibly with a smaller bandwidth) will do.
!               In either case, PD is preset to zero by the solver,
!          so that only the nonzero elements need be loaded by JAC.
!          Each call to JAC is preceded by a call to F with the same
!          arguments NEQ, T, and Y.  Thus to gain some efficiency,
!          intermediate quantities shared by both calculations may be
!          saved in a user COMMON block by F and not recomputed by JAC,
!          if desired.  Also, JAC may alter the Y array, if desired.
!          JAC must be declared external in the calling program.
!               Subroutine JAC may access user-defined real and integer
!          work arrays, RPAR and IPAR, whose dimensions are set by the
!          user in the main program.
!
! MF     = The method flag.  Used only for input.  The legal values of
!          MF are 10, 11, 12, 13, 14, 15, 20, 21, 22, 23, 24, 25,
!          -11, -12, -14, -15, -21, -22, -24, -25.
!          MF is a signed two-digit integer, MF = JSV*(10*METH + MITER).
!          JSV = SIGN(MF) indicates the Jacobian-saving strategy..
!            JSV =  1 means a copy of the Jacobian is saved for reuse
!                     in the corrector iteration algorithm.
!            JSV = -1 means a copy of the Jacobian is not saved
!                     (valid only for MITER = 1, 2, 4, or 5).
!          METH indicates the basic linear multistep method..
!            METH = 1 means the implicit Adams method.
!            METH = 2 means the method based on backward
!                     differentiation formulas (BDF-s).
!          MITER indicates the corrector iteration method..
!            MITER = 0 means functional iteration (no Jacobian matrix
!                      is involved).
!            MITER = 1 means chord iteration with a user-supplied
!                      full (NEQ by NEQ) Jacobian.
!            MITER = 2 means chord iteration with an internally
!                      generated (difference quotient) full Jacobian
!                      (using NEQ extra calls to F per df/dy value).
!            MITER = 3 means chord iteration with an internally
!                      generated diagonal Jacobian approximation
!                      (using 1 extra call to F per df/dy evaluation).
!            MITER = 4 means chord iteration with a user-supplied
!                      banded Jacobian.
!            MITER = 5 means chord iteration with an internally
!                      generated banded Jacobian (using ML+MU+1 extra
!                      calls to F per df/dy evaluation).
!          If MITER = 1 or 4, the user must supply a subroutine JAC
!          (the name is arbitrary) as described above under JAC.
!          For other values of MITER, a dummy argument can be used.
!
! RPAR     User-specified array used to communicate real parameters
!          to user-supplied subroutines.  If RPAR is a vector, then
!          it must be dimensioned in the user's main program.  If it
!          is unused or it is a scalar, then it need not be
!          dimensioned.
!
! IPAR     User-specified array used to communicate integer parameter
!          to user-supplied subroutines.  The comments on dimensioning
!          RPAR apply to IPAR.
!-----------------------------------------------------------------------
! Optional Input.
!
! The following is a list of the optional input provided for in the
! call sequence.  (See also Part ii.)  For each such input variable,
! this table lists its name as used in this documentation, its
! location in the call sequence, its meaning, and the default value.
! The use of any of this input requires IOPT = 1, and in that
! case all of this input is examined.  A value of zero for any
! of these optional input variables will cause the default value to be
! used.  Thus to use a subset of the optional input, simply preload
! locations 5 to 10 in RWORK and IWORK to 0.0 and 0 respectively, and
! then set those of interest to nonzero values.
!
! NAME    LOCATION      MEANING AND DEFAULT VALUE
!
! H0      RWORK(5)  The step size to be attempted on the first step.
!                   The default value is determined by the solver.
!
! HMAX    RWORK(6)  The maximum absolute step size allowed.
!                   The default value is infinite.
!
! HMIN    RWORK(7)  The minimum absolute step size allowed.
!                   The default value is 0.  (This lower bound is not
!                   enforced on the final step before reaching TCRIT
!                   when ITASK = 4 or 5.)
!
! MAXORD  IWORK(5)  The maximum order to be allowed.  The default
!                   value is 12 if METH = 1, and 5 if METH = 2.
!                   If MAXORD exceeds the default value, it will
!                   be reduced to the default value.
!                   If MAXORD is changed during the problem, it may
!                   cause the current order to be reduced.
!
! MXSTEP  IWORK(6)  Maximum number of (internally defined) steps
!                   allowed during one call to the solver.
!                   The default value is 500.
!
! MXHNIL  IWORK(7)  Maximum number of messages printed (per problem)
!                   warning that T + H = T on a step (H = step size).
!                   This must be positive to result in a non-default
!                   value.  The default value is 10.
!
!-----------------------------------------------------------------------
! Optional Output.
!
! As optional additional output from SVODE, the variables listed
! below are quantities related to the performance of SVODE
! which are available to the user.  These are communicated by way of
! the work arrays, but also have internal mnemonic names as shown.
! Except where stated otherwise, all of this output is defined
! on any successful return from SVODE, and on any return with
! ISTATE = -1, -2, -4, -5, or -6.  On an illegal input return
! (ISTATE = -3), they will be unchanged from their existing values
! (if any), except possibly for TOLSF, LENRW, and LENIW.
! On any error return, output relevant to the error will be defined,
! as noted below.
!
! NAME    LOCATION      MEANING
!
! HU      RWORK(11) The step size in t last used (successfully).
!
! HCUR    RWORK(12) The step size to be attempted on the next step.
!
! TCUR    RWORK(13) The current value of the independent variable
!                   which the solver has actually reached, i.e. the
!                   current internal mesh point in t.  In the output,
!                   TCUR will always be at least as far from the
!                   initial value of t as the current argument T,
!                   but may be farther (if interpolation was done).
!
! TOLSF   RWORK(14) A tolerance scale factor, greater than 1.0,
!                   computed when a request for too much accuracy was
!                   detected (ISTATE = -3 if detected at the start of
!                   the problem, ISTATE = -2 otherwise).  If ITOL is
!                   left unaltered but RTOL and ATOL are uniformly
!                   scaled up by a factor of TOLSF for the next call,
!                   then the solver is deemed likely to succeed.
!                   (The user may also ignore TOLSF and alter the
!                   tolerance parameters in any other way appropriate.)
!
! NST     IWORK(11) The number of steps taken for the problem so far.
!
! NFE     IWORK(12) The number of f evaluations for the problem so far.
!
! NJE     IWORK(13) The number of Jacobian evaluations so far.
!
! NQU     IWORK(14) The method order last used (successfully).
!
! NQCUR   IWORK(15) The order to be attempted on the next step.
!
! IMXER   IWORK(16) The index of the component of largest magnitude in
!                   the weighted local error vector ( e(i)/EWT(i) ),
!                   on an error return with ISTATE = -4 or -5.
!
! LENRW   IWORK(17) The length of RWORK actually required.
!                   This is defined on normal returns and on an illegal
!                   input return for insufficient storage.
!
! LENIW   IWORK(18) The length of IWORK actually required.
!                   This is defined on normal returns and on an illegal
!                   input return for insufficient storage.
!
! NLU     IWORK(19) The number of matrix LU decompositions so far.
!
! NNI     IWORK(20) The number of nonlinear (Newton) iterations so far.
!
! NCFN    IWORK(21) The number of convergence failures of the nonlinear
!                   solver so far.
!
! NETF    IWORK(22) The number of error test failures of the integrator
!                   so far.
!
! The following two arrays are segments of the RWORK array which
! may also be of interest to the user as optional output.
! For each array, the table below gives its internal name,
! its base address in RWORK, and its description.
!
! NAME    BASE ADDRESS      DESCRIPTION
!
! YH      21             The Nordsieck history array, of size NYH by
!                        (NQCUR + 1), where NYH is the initial value
!                        of NEQ.  For j = 0,1,...,NQCUR, column j+1
!                        of YH contains HCUR**j/factorial(j) times
!                        the j-th derivative of the interpolating
!                        polynomial currently representing the
!                        solution, evaluated at t = TCUR.
!
! ACOR     LENRW-NEQ+1   Array of size NEQ used for the accumulated
!                        corrections on each step, scaled in the output
!                        to represent the estimated local error in Y
!                        on the last step.  This is the vector e in
!                        the description of the error control.  It is
!                        defined only on a successful return from SVODE.
!
!-----------------------------------------------------------------------
! Interrupting and Restarting
!
! If the integration of a given problem by SVODE is to be
! interrrupted and then later continued, such as when restarting
! an interrupted run or alternating between two or more ODE problems,
! the user should save, following the return from the last SVODE call
! prior to the interruption, the contents of the call sequence
! variables and internal COMMON blocks, and later restore these
! values before the next SVODE call for that problem.  To save
! and restore the COMMON blocks, use subroutine SVSRCO, as
! described below in part ii.
!
! In addition, if non-default values for either LUN or MFLAG are
! desired, an extra call to XSETUN and/or XSETF should be made just
! before continuing the integration.  See Part ii below for details.
!
!-----------------------------------------------------------------------
! Part ii.  Other Routines Callable.
!
! The following are optional calls which the user may make to
! gain additional capabilities in conjunction with SVODE.
! (The routines XSETUN and XSETF are designed to conform to the
! SLATEC error handling package.)
!
!     FORM OF CALL                  FUNCTION
!  CALL XSETUN(LUN)           Set the logical unit number, LUN, for
!                             output of messages from SVODE, if
!                             the default is not desired.
!                             The default value of LUN is 6.
!
!  CALL XSETF(MFLAG)          Set a flag to control the printing of
!                             messages by SVODE.
!                             MFLAG = 0 means do not print. (Danger..
!                             This risks losing valuable information.)
!                             MFLAG = 1 means print (the default).
!
!                             Either of the above calls may be made at
!                             any time and will take effect immediately.
!
!  CALL SVSRCO(RSAV,ISAV,JOB) Saves and restores the contents of
!                             the internal COMMON blocks used by
!                             SVODE. (See Part iii below.)
!                             RSAV must be a real array of length 49
!                             or more, and ISAV must be an integer
!                             array of length 40 or more.
!                             JOB=1 means save COMMON into RSAV/ISAV.
!                             JOB=2 means restore COMMON from RSAV/ISAV.
!                                SVSRCO is useful if one is
!                             interrupting a run and restarting
!                             later, or alternating between two or
!                             more problems solved with SVODE.
!
!  CALL SVINDY(,,,,,)         Provide derivatives of y, of various
!        (See below.)         orders, at a specified point T, if
!                             desired.  It may be called only after
!                             a successful return from SVODE.
!
! The detailed instructions for using SVINDY are as follows.
! The form of the call is..
!
!  CALL SVINDY (T, K, RWORK(21), NYH, DKY, IFLAG)
!
! The input parameters are..
!
! T         = Value of independent variable where answers are desired
!             (normally the same as the T last returned by SVODE).
!             For valid results, T must lie between TCUR - HU and TCUR.
!             (See optional output for TCUR and HU.)
! K         = Integer order of the derivative desired.  K must satisfy
!             0 .le. K .le. NQCUR, where NQCUR is the current order
!             (see optional output).  The capability corresponding
!             to K = 0, i.e. computing y(T), is already provided
!             by SVODE directly.  Since NQCUR .ge. 1, the first
!             derivative dy/dt is always available with SVINDY.
! RWORK(21) = The base address of the history array YH.
! NYH       = Column length of YH, equal to the initial value of NEQ.
!
! The output parameters are..
!
! DKY       = A real array of length NEQ containing the computed value
!             of the K-th derivative of y(t).
! IFLAG     = Integer flag, returned as 0 if K and T were legal,
!             -1 if K was illegal, and -2 if T was illegal.
!             On an error return, a message is also written.
!-----------------------------------------------------------------------
! Part iii.  COMMON Blocks.
! If SVODE is to be used in an overlay situation, the user
! must declare, in the primary overlay, the variables in..
!   (1) the call sequence to SVODE,
!   (2) the two internal COMMON blocks
!         /SVOD01/  of length  81  (48 single precision words
!                         followed by 33 integer words),
!         /SVOD02/  of length  9  (1 single precision word
!                         followed by 8 integer words),
!
! If SVODE is used on a system in which the contents of internal
! COMMON blocks are not preserved between calls, the user should
! declare the above two COMMON blocks in his main program to insure
! that their contents are preserved.
!
!-----------------------------------------------------------------------
! Part iv.  Optionally Replaceable Solver Routines.
!
! Below are descriptions of two routines in the SVODE package which
! relate to the measurement of errors.  Either routine can be
! replaced by a user-supplied version, if desired.  However, since such
! a replacement may have a major impact on performance, it should be
! done only when absolutely necessary, and only with great caution.
! (Note.. The means by which the package version of a routine is
! superseded by the user's version may be system-dependent.)
!
! (a) SEWSET.
! The following subroutine is called just before each internal
! integration step, and sets the array of error weights, EWT, as
! described under ITOL/RTOL/ATOL above..
!     SUBROUTINE SEWSET (NEQ, ITOL, RTOL, ATOL, YCUR, EWT)
! where NEQ, ITOL, RTOL, and ATOL are as in the SVODE call sequence,
! YCUR contains the current dependent variable vector, and
! EWT is the array of weights set by SEWSET.
!
! If the user supplies this subroutine, it must return in EWT(i)
! (i = 1,...,NEQ) a positive quantity suitable for comparison with
! errors in Y(i).  The EWT array returned by SEWSET is passed to the
! SVNORM routine (See below.), and also used by SVODE in the computation
! of the optional output IMXER, the diagonal Jacobian approximation,
! and the increments for difference quotient Jacobians.
!
! In the user-supplied version of SEWSET, it may be desirable to use
! the current values of derivatives of y.  Derivatives up to order NQ
! are available from the history array YH, described above under
! Optional Output.  In SEWSET, YH is identical to the YCUR array,
! extended to NQ + 1 columns with a column length of NYH and scale
! factors of h**j/factorial(j).  On the first call for the problem,
! given by NST = 0, NQ is 1 and H is temporarily set to 1.0.
! NYH is the initial value of NEQ.  The quantities NQ, H, and NST
! can be obtained by including in SEWSET the statements..
!     REAL RVOD, H, HU
!     COMMON /SVOD01/ RVOD(48), IVOD(33)
!     COMMON /SVOD02/ HU, NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST
!     NQ = IVOD(28)
!     H = RVOD(21)
! Thus, for example, the current value of dy/dt can be obtained as
! YCUR(NYH+i)/H  (i=1,...,NEQ)  (and the division by H is
! unnecessary when NST = 0).
!
! (b) SVNORM.
! The following is a real function routine which computes the weighted
! root-mean-square norm of a vector v..
!     D = SVNORM (N, V, W)
! where..
!   N = the length of the vector,
!   V = real array of length N containing the vector,
!   W = real array of length N containing weights,
!   D = sqrt( (1/N) * sum(V(i)*W(i))**2 ).
! SVNORM is called with N = NEQ and with W(i) = 1.0/EWT(i), where
! EWT is as set by subroutine SEWSET.
!
! If the user supplies this function, it should return a non-negative
! value of SVNORM suitable for use in the error control in SVODE.
! None of the arguments should be altered by SVNORM.
! For example, a user-supplied SVNORM routine might..
!   -substitute a max-norm of (V(i)*W(i)) for the rms-norm, or
!   -ignore some components of V in the norm, with the effect of
!    suppressing the error control on those components of Y.
!-----------------------------------------------------------------------
! Other Routines in the SVODE Package.
!
! In addition to subroutine SVODE, the SVODE package includes the
! following subroutines and function routines..
!  SVHIN     computes an approximate step size for the initial step.
!  SVINDY    computes an interpolated value of the y vector at t = TOUT.
!  SVSTEP    is the core integrator, which does one step of the
!            integration and the associated error control.
!  SVSET     sets all method coefficients and test constants.
!  SVNLSD    solves the underlying nonlinear system -- the corrector.
!  SVJAC     computes and preprocesses the Jacobian matrix J = df/dy
!            and the Newton iteration matrix P = I - (h/l1)*J.
!  SVSOL     manages solution of linear system in chord iteration.
!  SVJUST    adjusts the history array on a change of order.
!  SEWSET    sets the error weight vector EWT before each step.
!  SVNORM    computes the weighted r.m.s. norm of a vector.
!  SVSRCO    is a user-callable routines to save and restore
!            the contents of the internal COMMON blocks.
!  SACOPY    is a routine to copy one two-dimensional array to another.
!  SGEFA and SGESL   are routines from LINPACK for solving full
!            systems of linear algebraic equations.
!  SGBFA and SGBSL   are routines from LINPACK for solving banded
!            linear systems.
!  SAXPY, SSCAL, and CH_SCOPY are basic linear algebra modules (BLAS).
!  R1MACH    sets the unit roundoff of the machine.
!  XERRWV, XSETUN, XSETF, LUNSAV, and MFLGSV handle the printing of all
!            error messages and warnings.  XERRWV is machine-dependent.
! Note..  SVNORM, R1MACH, LUNSAV, and MFLGSV are function routines.
! All the others are subroutines.
!
! The intrinsic and external routines used by the SVODE package are..
! ABS, MAX, MIN, REAL, SIGN, SQRT, and WRITE.
!
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
      EXTERNAL SVNLSD
      LOGICAL IHIT
      REAL ATOLI, BIG, EWTI, FOUR, H0, HMAX, HMX, HUN, ONE,             &
     &   PT2, RH, RTOLI, SIZE, TCRIT, TNEXT, TOLSF, TP, TWO, ZERO
      INTEGER I, IER, IFLAG, IMXER, JCO, KGO, LENIW, LENJ, LENP, LENRW, &
     &   LENWM, LF0, MBAND, ML, MORD, MU, MXHNL0, MXSTP0, NITER, NSLAST
      CHARACTER*80 MSG
!
! Type declaration for function subroutines called ---------------------
!
      REAL R1MACH, SVNORM
!
      DIMENSION MORD(2)
!-----------------------------------------------------------------------
! The following Fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to SVODE.
!-----------------------------------------------------------------------
      SAVE MORD, MXHNL0, MXSTP0
      SAVE ZERO, ONE, TWO, FOUR, PT2, HUN
!-----------------------------------------------------------------------
! The following internal COMMON blocks contain variables which are
! communicated between subroutines in the SVODE package, or which are
! to be saved between calls to SVODE.
! In each block, real variables precede integers.
! The block /SVOD01/ appears in subroutines SVODE, SVINDY, SVSTEP,
! SVSET, SVNLSD, SVJAC, SVSOL, SVJUST and SVSRCO.
! The block /SVOD02/ appears in subroutines SVODE, SVINDY, SVSTEP,
! SVNLSD, SVJAC, and SVSRCO.
!
! The variables stored in the internal COMMON blocks are as follows..
!
! ACNRM  = Weighted r.m.s. norm of accumulated correction vectors.
! CCMXJ  = Threshhold on DRC for updating the Jacobian. (See DRC.)
! CONP   = The saved value of TQ(5).
! CRATE  = Estimated corrector convergence rate constant.
! DRC    = Relative change in H*RL1 since last SVJAC call.
! EL     = Real array of integration coefficients.  See SVSET.
! ETA    = Saved tentative ratio of new to old H.
! ETAMAX = Saved maximum value of ETA to be allowed.
! H      = The step size.
! HMIN   = The minimum absolute value of the step size H to be used.
! HMXI   = Inverse of the maximum absolute value of H to be used.
!          HMXI = 0.0 is allowed and corresponds to an infinite HMAX.
! HNEW   = The step size to be attempted on the next step.
! HSCAL  = Stepsize in scaling of YH array.
! PRL1   = The saved value of RL1.
! RC     = Ratio of current H*RL1 to value on last SVJAC call.
! RL1    = The reciprocal of the coefficient EL(1).
! TAU    = Real vector of past NQ step sizes, length 13.
! TQ     = A real vector of length 5 in which SVSET stores constants
!          used for the convergence test, the error test, and the
!          selection of H at a new order.
! TN     = The independent variable, updated on each step taken.
! UROUND = The machine unit roundoff.  The smallest positive real number
!          such that  1.0 + UROUND .ne. 1.0
! ICF    = Integer flag for convergence failure in SVNLSD..
!            0 means no failures.
!            1 means convergence failure with out of date Jacobian
!                   (recoverable error).
!            2 means convergence failure with current Jacobian or
!                   singular matrix (unrecoverable error).
! INIT   = Saved integer flag indicating whether initialization of the
!          problem has been done (INIT = 1) or not.
! IPUP   = Saved flag to signal updating of Newton matrix.
! JCUR   = Output flag from SVJAC showing Jacobian status..
!            JCUR = 0 means J is not current.
!            JCUR = 1 means J is current.
! JSTART = Integer flag used as input to SVSTEP..
!            0  means perform the first step.
!            1  means take a new step continuing from the last.
!            -1 means take the next step with a new value of MAXORD,
!                  HMIN, HMXI, N, METH, MITER, and/or matrix parameters.
!          On return, SVSTEP sets JSTART = 1.
! JSV    = Integer flag for Jacobian saving, = sign(MF).
! KFLAG  = A completion code from SVSTEP with the following meanings..
!               0      the step was succesful.
!              -1      the requested error could not be achieved.
!              -2      corrector convergence could not be achieved.
!              -3, -4  fatal error in VNLS (can not occur here).
! KUTH   = Input flag to SVSTEP showing whether H was reduced by the
!          driver.  KUTH = 1 if H was reduced, = 0 otherwise.
! L      = Integer variable, NQ + 1, current order plus one.
! LMAX   = MAXORD + 1 (used for dimensioning).
! LOCJS  = A pointer to the saved Jacobian, whose storage starts at
!          WM(LOCJS), if JSV = 1.
! LYH, LEWT, LACOR, LSAVF, LWM, LIWM = Saved integer pointers
!          to segments of RWORK and IWORK.
! MAXORD = The maximum order of integration method to be allowed.
! METH/MITER = The method flags.  See MF.
! MSBJ   = The maximum number of steps between J evaluations, = 50.
! MXHNIL = Saved value of optional input MXHNIL.
! MXSTEP = Saved value of optional input MXSTEP.
! N      = The number of first-order ODEs, = NEQ.
! NEWH   = Saved integer to flag change of H.
! NEWQ   = The method order to be used on the next step.
! NHNIL  = Saved counter for occurrences of T + H = T.
! NQ     = Integer variable, the current integration method order.
! NQNYH  = Saved value of NQ*NYH.
! NQWAIT = A counter controlling the frequency of order changes.
!          An order change is about to be considered if NQWAIT = 1.
! NSLJ   = The number of steps taken as of the last Jacobian update.
! NSLP   = Saved value of NST as of last Newton matrix update.
! NYH    = Saved value of the initial value of NEQ.
! HU     = The step size in t last used.
! NCFN   = Number of nonlinear convergence failures so far.
! NETF   = The number of error test failures of the integrator so far.
! NFE    = The number of f evaluations for the problem so far.
! NJE    = The number of Jacobian evaluations so far.
! NLU    = The number of matrix LU decompositions so far.
! NNI    = Number of nonlinear iterations so far.
! NQU    = The method order last used.
! NST    = The number of steps taken for the problem so far.
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
      DATA  MORD(1) /12/, MORD(2) /5/, MXSTP0 /500/, MXHNL0 /10/
      DATA ZERO /0.0E0/, ONE /1.0E0/, TWO /2.0E0/, FOUR /4.0E0/,        &
     &     PT2 /0.2E0/, HUN /100.0E0/
!-----------------------------------------------------------------------
! Block A.
! This code block is executed on every call.
! It tests ISTATE and ITASK for legality and branches appropriately.
! If ISTATE .gt. 1 but the flag INIT shows that initialization has
! not yet been done, an error return occurs.
! If ISTATE = 1 and TOUT = T, return immediately.
!-----------------------------------------------------------------------
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('SVODE',0,ZHOOK_HANDLE)
      IF (ISTATE .LT. 1 .OR. ISTATE .GT. 3) GO TO 601
      IF (ITASK .LT. 1 .OR. ITASK .GT. 5) GO TO 602
      IF (ISTATE .EQ. 1) GO TO 10
      IF (INIT .NE. 1) GO TO 603
      IF (ISTATE .EQ. 2) GO TO 200
      GO TO 20
 10   INIT = 0
      IF (TOUT .EQ. T)  THEN
        IF (LHOOK) CALL DR_HOOK('SVODE',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
!-----------------------------------------------------------------------
! Block B.
! The next code block is executed for the initial call (ISTATE = 1),
! or for a continuation call with parameter changes (ISTATE = 3).
! It contains checking of all input and various initializations.
!
! First check legality of the non-optional input NEQ, ITOL, IOPT,
! MF, ML, and MU.
!-----------------------------------------------------------------------
 20   IF (NEQ .LE. 0) GO TO 604
      IF (ISTATE .EQ. 1) GO TO 25
      IF (NEQ .GT. N) GO TO 605
 25   N = NEQ
      IF (ITOL .LT. 1 .OR. ITOL .GT. 4) GO TO 606
      IF (IOPT .LT. 0 .OR. IOPT .GT. 1) GO TO 607
      JSV = SIGN(1,MF)
      MF = ABS(MF)
      METH = MF/10
      MITER = MF - 10*METH
      IF (METH .LT. 1 .OR. METH .GT. 2) GO TO 608
      IF (MITER .LT. 0 .OR. MITER .GT. 5) GO TO 608
      IF (MITER .LE. 3) GO TO 30
      ML = IWORK(1)
      MU = IWORK(2)
      IF (ML .LT. 0 .OR. ML .GE. N) GO TO 609
      IF (MU .LT. 0 .OR. MU .GE. N) GO TO 610
 30   CONTINUE
! Next process and check the optional input. ---------------------------
      IF (IOPT .EQ. 1) GO TO 40
      MAXORD = MORD(METH)
      MXSTEP = MXSTP0
      MXHNIL = MXHNL0
      IF (ISTATE .EQ. 1) H0 = ZERO
      HMXI = ZERO
      HMIN = ZERO
      GO TO 60
 40   MAXORD = IWORK(5)
      IF (MAXORD .LT. 0) GO TO 611
      IF (MAXORD .EQ. 0) MAXORD = 100
      MAXORD = MIN(MAXORD,MORD(METH))
      MXSTEP = IWORK(6)
      IF (MXSTEP .LT. 0) GO TO 612
      IF (MXSTEP .EQ. 0) MXSTEP = MXSTP0
      MXHNIL = IWORK(7)
      IF (MXHNIL .LT. 0) GO TO 613
      IF (MXHNIL .EQ. 0) MXHNIL = MXHNL0
      IF (ISTATE .NE. 1) GO TO 50
      H0 = RWORK(5)
      IF ((TOUT - T)*H0 .LT. ZERO) GO TO 614
 50   HMAX = RWORK(6)
      IF (HMAX .LT. ZERO) GO TO 615
      HMXI = ZERO
      IF (HMAX .GT. ZERO) HMXI = ONE/HMAX
      HMIN = RWORK(7)
      IF (HMIN .LT. ZERO) GO TO 616
!-----------------------------------------------------------------------
! Set work array pointers and check lengths LRW and LIW.
! Pointers to segments of RWORK and IWORK are named by prefixing L to
! the name of the segment.  E.g., the segment YH starts at RWORK(LYH).
! Segments of RWORK (in order) are denoted  YH, WM, EWT, SAVF, ACOR.
! Within WM, LOCJS is the location of the saved Jacobian (JSV .gt. 0).
!-----------------------------------------------------------------------
 60   LYH = 21
      IF (ISTATE .EQ. 1) NYH = N
      LWM = LYH + (MAXORD + 1)*NYH
      JCO = MAX(0,JSV)
      IF (MITER .EQ. 0) LENWM = 0
      IF (MITER .EQ. 1 .OR. MITER .EQ. 2) THEN
        LENWM = 2 + (1 + JCO)*N*N
        LOCJS = N*N + 3
      ENDIF
      IF (MITER .EQ. 3) LENWM = 2 + N
      IF (MITER .EQ. 4 .OR. MITER .EQ. 5) THEN
        MBAND = ML + MU + 1
        LENP = (MBAND + ML)*N
        LENJ = MBAND*N
        LENWM = 2 + LENP + JCO*LENJ
        LOCJS = LENP + 3
        ENDIF
      LEWT = LWM + LENWM
      LSAVF = LEWT + N
      LACOR = LSAVF + N
      LENRW = LACOR + N - 1
      IWORK(17) = LENRW
      LIWM = 1
      LENIW = 30 + N
      IF (MITER .EQ. 0 .OR. MITER .EQ. 3) LENIW = 30
      IWORK(18) = LENIW
      IF (LENRW .GT. LRW) GO TO 617
      IF (LENIW .GT. LIW) GO TO 618
! Check RTOL and ATOL for legality. ------------------------------------
      RTOLI = RTOL(1)
      ATOLI = ATOL(1)
      DO 70 I = 1,N
        IF (ITOL .GE. 3) RTOLI = RTOL(I)
        IF (ITOL .EQ. 2 .OR. ITOL .EQ. 4) ATOLI = ATOL(I)
        IF (RTOLI .LT. ZERO) GO TO 619
        IF (ATOLI .LT. ZERO) GO TO 620
 70     CONTINUE
      IF (ISTATE .EQ. 1) GO TO 100
! If ISTATE = 3, set flag to signal parameter changes to SVSTEP. -------
      JSTART = -1
      IF (NQ .LE. MAXORD) GO TO 90
! MAXORD was reduced below NQ.  Copy YH(*,MAXORD+2) into SAVF. ---------
      CALL CH_SCOPY (N, RWORK(LWM), 1, RWORK(LSAVF), 1)
! Reload WM(1) = RWORK(LWM), since LWM may have changed. ---------------
 90   IF (MITER .GT. 0) RWORK(LWM) = SQRT(UROUND)
!-----------------------------------------------------------------------
! Block C.
! The next block is for the initial call only (ISTATE = 1).
! It contains all remaining initializations, the initial call to F,
! and the calculation of the initial step size.
! The error weights in EWT are inverted after being loaded.
!-----------------------------------------------------------------------
 100  UROUND = R1MACH(4)
      TN = T
      IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GO TO 110
      TCRIT = RWORK(1)
      IF ((TCRIT - TOUT)*(TOUT - T) .LT. ZERO) GO TO 625
      IF (H0 .NE. ZERO .AND. (T + H0 - TCRIT)*H0 .GT. ZERO)             &
     &   H0 = TCRIT - T
 110  JSTART = 0
      IF (MITER .GT. 0) RWORK(LWM) = SQRT(UROUND)
      CCMXJ = PT2
      MSBJ = 50
      NHNIL = 0
      NST = 0
      NJE = 0
      NNI = 0
      NCFN = 0
      NETF = 0
      NLU = 0
      NSLJ = 0
      NSLAST = 0
      HU = ZERO
      NQU = 0
! Initial call to F.  (LF0 points to YH(*,2).) -------------------------
      LF0 = LYH + NYH
!
!*UPG*MNH
!
      CALL F (N, T, Y, RWORK(LF0), RPAR, IPAR, KMI, KINDEX)
!
!*UPG*MNH
!
      NFE = 1
! Load the initial value vector in YH. ---------------------------------
      CALL CH_SCOPY (N, Y, 1, RWORK(LYH), 1)
! Load and invert the EWT array.  (H is temporarily set to 1.0.) -------
      NQ = 1
      H = ONE
      CALL SEWSET (N, ITOL, RTOL, ATOL, RWORK(LYH), RWORK(LEWT))
      DO 120 I = 1,N
        IF (RWORK(I+LEWT-1) .LE. ZERO) GO TO 621
 120    RWORK(I+LEWT-1) = ONE/RWORK(I+LEWT-1)
      IF (H0 .NE. ZERO) GO TO 180
! Call SVHIN to set initial step size H0 to be attempted. --------------
      CALL SVHIN (N, T, RWORK(LYH), RWORK(LF0), F, RPAR, IPAR, TOUT,    &
     &   UROUND, RWORK(LEWT), ITOL, ATOL, Y, RWORK(LACOR), H0,          &
     &   NITER, IER, KMI, KINDEX)
      NFE = NFE + NITER
      IF (IER .NE. 0) GO TO 622
! Adjust H0 if necessary to meet HMAX bound. ---------------------------
 180  RH = ABS(H0)*HMXI
      IF (RH .GT. ONE) H0 = H0/RH
! Load H with H0 and scale YH(*,2) by H0. ------------------------------
      H = H0
      CALL SSCAL (N, H0, RWORK(LF0), 1)
      GO TO 270
!-----------------------------------------------------------------------
! Block D.
! The next code block is for continuation calls only (ISTATE = 2 or 3)
! and is to check stop conditions before taking a step.
!-----------------------------------------------------------------------
 200  NSLAST = NST
      KUTH = 0
      GO TO (210, 250, 220, 230, 240), ITASK
 210  IF ((TN - TOUT)*H .LT. ZERO) GO TO 250
      CALL SVINDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      IF (IFLAG .NE. 0) GO TO 627
      T = TOUT
      GO TO 420
 220  TP = TN - HU*(ONE + HUN*UROUND)
      IF ((TP - TOUT)*H .GT. ZERO) GO TO 623
      IF ((TN - TOUT)*H .LT. ZERO) GO TO 250
      GO TO 400
 230  TCRIT = RWORK(1)
      IF ((TN - TCRIT)*H .GT. ZERO) GO TO 624
      IF ((TCRIT - TOUT)*H .LT. ZERO) GO TO 625
      IF ((TN - TOUT)*H .LT. ZERO) GO TO 245
      CALL SVINDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      IF (IFLAG .NE. 0) GO TO 627
      T = TOUT
      GO TO 420
 240  TCRIT = RWORK(1)
      IF ((TN - TCRIT)*H .GT. ZERO) GO TO 624
 245  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. HUN*UROUND*HMX
      IF (IHIT) GO TO 400
      TNEXT = TN + HNEW*(ONE + FOUR*UROUND)
      IF ((TNEXT - TCRIT)*H .LE. ZERO) GO TO 250
      H = (TCRIT - TN)*(ONE - FOUR*UROUND)
      KUTH = 1
!-----------------------------------------------------------------------
! Block E.
! The next block is normally executed for all calls and contains
! the call to the one-step core integrator SVSTEP.
!
! This is a looping point for the integration steps.
!
! First check for too many steps being taken, update EWT (if not at
! start of problem), check for too much accuracy being requested, and
! check for H below the roundoff level in T.
!-----------------------------------------------------------------------
 250  CONTINUE
      IF ((NST-NSLAST) .GE. MXSTEP) GO TO 500
      CALL SEWSET (N, ITOL, RTOL, ATOL, RWORK(LYH), RWORK(LEWT))
      DO 260 I = 1,N
        IF (RWORK(I+LEWT-1) .LE. ZERO) GO TO 510
 260    RWORK(I+LEWT-1) = ONE/RWORK(I+LEWT-1)
 270  TOLSF = UROUND*SVNORM (N, RWORK(LYH), RWORK(LEWT))
      IF (TOLSF .LE. ONE) GO TO 280
      TOLSF = TOLSF*TWO
      IF (NST .EQ. 0) GO TO 626
      GO TO 520
 280  IF ((TN + H) .NE. TN) GO TO 290
!KS:  strange things happen on HP f90 (this error message is
!KS:  often printed but results are the same than on Linux)
!KS:  let's jump over these prints
      GOTO 290
      NHNIL = NHNIL + 1
      IF (NHNIL .GT. MXHNIL) GO TO 290
      MSG = 'SVODE--  Warning..internal T (=R1) and H (=R2) are'
      CALL XERRWV (MSG, 50, 101, 1, 0, 0, 0, 0, ZERO, ZERO)
      MSG='      such that in the machine, T + H = T on the next step  '
      CALL XERRWV (MSG, 60, 101, 1, 0, 0, 0, 0, ZERO, ZERO)
      MSG = '      (H = step size). solver will continue anyway'
      CALL XERRWV (MSG, 50, 101, 1, 0, 0, 0, 2, TN, H)
      IF (NHNIL .LT. MXHNIL) GO TO 290
      MSG = 'SVODE--  Above warning has been issued I1 times.  '
      CALL XERRWV (MSG, 50, 102, 1, 0, 0, 0, 0, ZERO, ZERO)
      MSG = '      it will not be issued again for this problem'
      CALL XERRWV (MSG, 50, 102, 1, 1, MXHNIL, 0, 0, ZERO, ZERO)
 290  CONTINUE
!-----------------------------------------------------------------------
! CALL SVSTEP (Y, YH, NYH, YH, EWT, SAVF, VSAV, ACOR,
!              WM, IWM, F, JAC, F, SVNLSD, RPAR, IPAR)
!-----------------------------------------------------------------------
      CALL SVSTEP (Y, RWORK(LYH), NYH, RWORK(LYH), RWORK(LEWT),         &
     &   RWORK(LSAVF), Y, RWORK(LACOR), RWORK(LWM), IWORK(LIWM),        &
     &   F, JAC, F, SVNLSD, RPAR, IPAR, KMI, KINDEX)
      KGO = 1 - KFLAG
! Branch on KFLAG.  Note..In this version, KFLAG can not be set to -3.
!  KFLAG .eq. 0,   -1,  -2
      GO TO (300, 530, 540), KGO
!-----------------------------------------------------------------------
! Block F.
! The following block handles the case of a successful return from the
! core integrator (KFLAG = 0).  Test for stop conditions.
!-----------------------------------------------------------------------
 300  INIT = 1
      KUTH = 0
      GO TO (310, 400, 330, 340, 350), ITASK
! ITASK = 1.  If TOUT has been reached, interpolate. -------------------
 310  IF ((TN - TOUT)*H .LT. ZERO) GO TO 250
      CALL SVINDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      T = TOUT
      GO TO 420
! ITASK = 3.  Jump to exit if TOUT was reached. ------------------------
 330  IF ((TN - TOUT)*H .GE. ZERO) GO TO 400
      GO TO 250
! ITASK = 4.  See if TOUT or TCRIT was reached.  Adjust H if necessary.
 340  IF ((TN - TOUT)*H .LT. ZERO) GO TO 345
      CALL SVINDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      T = TOUT
      GO TO 420
 345  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. HUN*UROUND*HMX
      IF (IHIT) GO TO 400
      TNEXT = TN + HNEW*(ONE + FOUR*UROUND)
      IF ((TNEXT - TCRIT)*H .LE. ZERO) GO TO 250
      H = (TCRIT - TN)*(ONE - FOUR*UROUND)
      KUTH = 1
      GO TO 250
! ITASK = 5.  See if TCRIT was reached and jump to exit. ---------------
 350  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. HUN*UROUND*HMX
!-----------------------------------------------------------------------
! Block G.
! The following block handles all successful returns from SVODE.
! If ITASK .ne. 1, Y is loaded from YH and T is set accordingly.
! ISTATE is set to 2, and the optional output is loaded into the work
! arrays before returning.
!-----------------------------------------------------------------------
 400  CONTINUE
      CALL CH_SCOPY (N, RWORK(LYH), 1, Y, 1)
      T = TN
      IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GO TO 420
      IF (IHIT) T = TCRIT
 420  ISTATE = 2
      RWORK(11) = HU
      RWORK(12) = HNEW
      RWORK(13) = TN
      IWORK(11) = NST
      IWORK(12) = NFE
      IWORK(13) = NJE
      IWORK(14) = NQU
      IWORK(15) = NEWQ
      IWORK(19) = NLU
      IWORK(20) = NNI
      IWORK(21) = NCFN
      IWORK(22) = NETF
      IF (LHOOK) CALL DR_HOOK('SVODE',1,ZHOOK_HANDLE)
      RETURN
!-----------------------------------------------------------------------
! Block H.
! The following block handles all unsuccessful returns other than
! those for illegal input.  First the error message routine is called.
! if there was an error test or convergence test failure, IMXER is set.
! Then Y is loaded from YH, T is set to TN, and the illegal input
! The optional output is loaded into the work arrays before returning.
!-----------------------------------------------------------------------
! The maximum number of steps was taken before reaching TOUT. ----------
 500  MSG = 'SVODE--  At current T (=R1), MXSTEP (=I1) steps   '
      CALL XERRWV (MSG, 50, 201, 1, 0, 0, 0, 0, ZERO, ZERO)
      MSG = '      taken on this call before reaching TOUT     '
      CALL XERRWV (MSG, 50, 201, 1, 1, MXSTEP, 0, 1, TN, ZERO)
      ISTATE = -1
      GO TO 580
! EWT(i) .le. 0.0 for some i (not at start of problem). ----------------
 510  EWTI = RWORK(LEWT+I-1)
      MSG = 'SVODE--  At T (=R1), EWT(I1) has become R2 .le. 0.'
      CALL XERRWV (MSG, 50, 202, 1, 1, I, 0, 2, TN, EWTI)
      ISTATE = -6
      GO TO 580
! Too much accuracy requested for machine precision. -------------------
 520  MSG = 'SVODE--  At T (=R1), too much accuracy requested  '
      CALL XERRWV (MSG, 50, 203, 1, 0, 0, 0, 0, ZERO, ZERO)
      MSG = '      for precision of machine..  see TOLSF (=R2) '
      CALL XERRWV (MSG, 50, 203, 1, 0, 0, 0, 2, TN, TOLSF)
      RWORK(14) = TOLSF
      ISTATE = -2
      GO TO 580
! KFLAG = -1.  Error test failed repeatedly or with ABS(H) = HMIN. -----
 530  MSG = 'SVODE--  At T(=R1) and step size H(=R2), the error'
      CALL XERRWV (MSG, 50, 204, 1, 0, 0, 0, 0, ZERO, ZERO)
      MSG = '      test failed repeatedly or with abs(H) = HMIN'
      CALL XERRWV (MSG, 50, 204, 1, 0, 0, 0, 2, TN, H)
      ISTATE = -4
      GO TO 560
! KFLAG = -2.  Convergence failed repeatedly or with abs(H) = HMIN. ----
 540  MSG = 'SVODE--  At T (=R1) and step size H (=R2), the    '
      CALL XERRWV (MSG, 50, 205, 1, 0, 0, 0, 0, ZERO, ZERO)
      MSG = '      corrector convergence failed repeatedly     '
      CALL XERRWV (MSG, 50, 205, 1, 0, 0, 0, 0, ZERO, ZERO)
      MSG = '      or with abs(H) = HMIN   '
      CALL XERRWV (MSG, 30, 205, 1, 0, 0, 0, 2, TN, H)
      ISTATE = -5
! Compute IMXER if relevant. -------------------------------------------
 560  BIG = ZERO
      IMXER = 1
      DO 570 I = 1,N
        SIZE = ABS(RWORK(I+LACOR-1)*RWORK(I+LEWT-1))
        IF (BIG .GE. SIZE) GO TO 570
        BIG = SIZE
        IMXER = I
 570    CONTINUE
      IWORK(16) = IMXER
! Set Y vector, T, and optional output. --------------------------------
 580  CONTINUE
      CALL CH_SCOPY (N, RWORK(LYH), 1, Y, 1)
      T = TN
      RWORK(11) = HU
      RWORK(12) = H
      RWORK(13) = TN
      IWORK(11) = NST
      IWORK(12) = NFE
      IWORK(13) = NJE
      IWORK(14) = NQU
      IWORK(15) = NQ
      IWORK(19) = NLU
      IWORK(20) = NNI
      IWORK(21) = NCFN
      IWORK(22) = NETF
      IF (LHOOK) CALL DR_HOOK('SVODE',1,ZHOOK_HANDLE)
      RETURN
!-----------------------------------------------------------------------
! Block I.
! The following block handles all error returns due to illegal input
! (ISTATE = -3), as detected before calling the core integrator.
! First the error message routine is called.   If the illegal input
! is a negative ISTATE, the run is aborted (apparent infinite loop).
!-----------------------------------------------------------------------
 601  MSG = 'SVODE--  ISTATE (=I1) illegal '
      CALL XERRWV (MSG, 30, 1, 1, 1, ISTATE, 0, 0, ZERO, ZERO)
      IF (ISTATE .LT. 0) GO TO 800
      GO TO 700
 602  MSG = 'SVODE--  ITASK (=I1) illegal  '
      CALL XERRWV (MSG, 30, 2, 1, 1, ITASK, 0, 0, ZERO, ZERO)
      GO TO 700
 603  MSG='SVODE--  ISTATE (=I1) .gt. 1 but SVODE not initialized      '
      CALL XERRWV (MSG, 60, 3, 1, 1, ISTATE, 0, 0, ZERO, ZERO)
      GO TO 700
 604  MSG = 'SVODE--  NEQ (=I1) .lt. 1     '
      CALL XERRWV (MSG, 30, 4, 1, 1, NEQ, 0, 0, ZERO, ZERO)
      GO TO 700
 605  MSG = 'SVODE--  ISTATE = 3 and NEQ increased (I1 to I2)  '
      CALL XERRWV (MSG, 50, 5, 1, 2, N, NEQ, 0, ZERO, ZERO)
      GO TO 700
 606  MSG = 'SVODE--  ITOL (=I1) illegal   '
      CALL XERRWV (MSG, 30, 6, 1, 1, ITOL, 0, 0, ZERO, ZERO)
      GO TO 700
 607  MSG = 'SVODE--  IOPT (=I1) illegal   '
      CALL XERRWV (MSG, 30, 7, 1, 1, IOPT, 0, 0, ZERO, ZERO)
      GO TO 700
 608  MSG = 'SVODE--  MF (=I1) illegal     '
      CALL XERRWV (MSG, 30, 8, 1, 1, MF, 0, 0, ZERO, ZERO)
      GO TO 700
 609  MSG = 'SVODE--  ML (=I1) illegal.. .lt.0 or .ge.NEQ (=I2)'
      CALL XERRWV (MSG, 50, 9, 1, 2, ML, NEQ, 0, ZERO, ZERO)
      GO TO 700
 610  MSG = 'SVODE--  MU (=I1) illegal.. .lt.0 or .ge.NEQ (=I2)'
      CALL XERRWV (MSG, 50, 10, 1, 2, MU, NEQ, 0, ZERO, ZERO)
      GO TO 700
 611  MSG = 'SVODE--  MAXORD (=I1) .lt. 0  '
      CALL XERRWV (MSG, 30, 11, 1, 1, MAXORD, 0, 0, ZERO, ZERO)
      GO TO 700
 612  MSG = 'SVODE--  MXSTEP (=I1) .lt. 0  '
      CALL XERRWV (MSG, 30, 12, 1, 1, MXSTEP, 0, 0, ZERO, ZERO)
      GO TO 700
 613  MSG = 'SVODE--  MXHNIL (=I1) .lt. 0  '
      CALL XERRWV (MSG, 30, 13, 1, 1, MXHNIL, 0, 0, ZERO, ZERO)
      GO TO 700
 614  MSG = 'SVODE--  TOUT (=R1) behind T (=R2)      '
      CALL XERRWV (MSG, 40, 14, 1, 0, 0, 0, 2, TOUT, T)
      MSG = '      integration direction is given by H0 (=R1)  '
      CALL XERRWV (MSG, 50, 14, 1, 0, 0, 0, 1, H0, ZERO)
      GO TO 700
 615  MSG = 'SVODE--  HMAX (=R1) .lt. 0.0  '
      CALL XERRWV (MSG, 30, 15, 1, 0, 0, 0, 1, HMAX, ZERO)
      GO TO 700
 616  MSG = 'SVODE--  HMIN (=R1) .lt. 0.0  '
      CALL XERRWV (MSG, 30, 16, 1, 0, 0, 0, 1, HMIN, ZERO)
      GO TO 700
 617  CONTINUE
      MSG='SVODE--  RWORK length needed, LENRW (=I1), exceeds LRW (=I2)'
      CALL XERRWV (MSG, 60, 17, 1, 2, LENRW, LRW, 0, ZERO, ZERO)
      GO TO 700
 618  CONTINUE
      MSG='SVODE--  IWORK length needed, LENIW (=I1), exceeds LIW (=I2)'
      CALL XERRWV (MSG, 60, 18, 1, 2, LENIW, LIW, 0, ZERO, ZERO)
      GO TO 700
 619  MSG = 'SVODE--  RTOL(I1) is R1 .lt. 0.0        '
      CALL XERRWV (MSG, 40, 19, 1, 1, I, 0, 1, RTOLI, ZERO)
      GO TO 700
 620  MSG = 'SVODE--  ATOL(I1) is R1 .lt. 0.0        '
      CALL XERRWV (MSG, 40, 20, 1, 1, I, 0, 1, ATOLI, ZERO)
      GO TO 700
 621  EWTI = RWORK(LEWT+I-1)
      MSG = 'SVODE--  EWT(I1) is R1 .le. 0.0         '
      CALL XERRWV (MSG, 40, 21, 1, 1, I, 0, 1, EWTI, ZERO)
      GO TO 700
 622  CONTINUE
      MSG='SVODE--  TOUT (=R1) too close to T(=R2) to start integration'
      CALL XERRWV (MSG, 60, 22, 1, 0, 0, 0, 2, TOUT, T)
      GO TO 700
 623  CONTINUE
      MSG='SVODE--  ITASK = I1 and TOUT (=R1) behind TCUR - HU (= R2)  '
      CALL XERRWV (MSG, 60, 23, 1, 1, ITASK, 0, 2, TOUT, TP)
      GO TO 700
 624  CONTINUE
      MSG='SVODE--  ITASK = 4 or 5 and TCRIT (=R1) behind TCUR (=R2)   '
      CALL XERRWV (MSG, 60, 24, 1, 0, 0, 0, 2, TCRIT, TN)
      GO TO 700
 625  CONTINUE
      MSG='SVODE--  ITASK = 4 or 5 and TCRIT (=R1) behind TOUT (=R2)   '
      CALL XERRWV (MSG, 60, 25, 1, 0, 0, 0, 2, TCRIT, TOUT)
      GO TO 700
 626  MSG = 'SVODE--  At start of problem, too much accuracy   '
      CALL XERRWV (MSG, 50, 26, 1, 0, 0, 0, 0, ZERO, ZERO)
      MSG='      requested for precision of machine..  see TOLSF (=R1) '
      CALL XERRWV (MSG, 60, 26, 1, 0, 0, 0, 1, TOLSF, ZERO)
      RWORK(14) = TOLSF
      GO TO 700
 627  MSG='SVODE--  Trouble from SVINDY.  ITASK = I1, TOUT = R1.       '
      CALL XERRWV (MSG, 60, 27, 1, 1, ITASK, 0, 1, TOUT, ZERO)
!
 700  CONTINUE
      ISTATE = -3
      IF (LHOOK) CALL DR_HOOK('SVODE',1,ZHOOK_HANDLE)
      RETURN
!
 800  MSG = 'SVODE--  Run aborted.. apparent infinite loop     '
      CALL XERRWV (MSG, 50, 303, 2, 0, 0, 0, 0, ZERO, ZERO)
      IF (LHOOK) CALL DR_HOOK('SVODE',1,ZHOOK_HANDLE)
      RETURN
      END SUBROUTINE SVODE
