!     ######spl
SUBROUTINE CH_CRANCK(PTSIMUL, PDTACT, PCONC, PNEWCONC, KEQ, KVECNPT, KMI, &
                       PALPHA)
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!!    #########################################################################
!!
!!*** *CH_CRANCK*
!!
!!    PURPOSE
!!    -------
!     solve one time-step of the chemical differential equation d/dt y = f(y)
!!
!!**  METHOD
!!    ------
!!    Cranck-Nicholson method (CRANCK):
!!    y^n+1 = y^n + dt*( alpha*f^n + (1-alpha)f^n-1 )
!!    for alpha = 1   the method reduces to EULER EXPLICIT
!!    for alpha = 0   the method reduces to EULER IMPLICIT
!!    for alpha = 1/2 the method reduces to EULER SEMI-IMPLICIT
!!    The implicit equation is solved by Newton-Raphson iteration:
!!    z^m+1 = z^m - (DF)(z^m)^-1 * F(z^m)
!!    where F(z) = z - y^n - dt*( alpha*f^n + (1-alpha)f(z) )
!!    the iteration process should yield:  z^m --> y^n+1 for increasing m
!!
!!    REFERENCE
!!    ---------
!!    J. Stoer: Einf\"uhrung in die Numerische Mathematik I & II,
!!    Heidelberger Taschenb\"ucher, Springer Verlag, Berlin, 1983 & 1978.
!!
!!    AUTHOR
!!    ------
!!    K. Suhre   *Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!    Original 12/06/95
!!    31/07/96 (K. Suhre) restructured
!!    19/04/02 add PALPHA argument
!!    01/12/03  (Gazen)   change Chemical scheme interface
!!    EXTERNAL
!!    --------
USE MODI_CH_FCN
USE MODI_CH_JAC
USE MODI_CH_GAUSS
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!        -----------------
IMPLICIT NONE
!
!*       0.1  declaration of arguments
!
REAL,    INTENT(IN) :: PTSIMUL  ! time of simulation
REAL,    INTENT(IN) :: PDTACT   ! actual time-step
INTEGER, INTENT(IN) :: KEQ      ! dimension of the problem to solve
INTEGER, INTENT(IN) :: KVECNPT
REAL, INTENT(IN),  DIMENSION(KVECNPT,KEQ) :: PCONC    ! concentration vector at PTSIMUL
REAL, INTENT(OUT), DIMENSION(KVECNPT,KEQ) :: PNEWCONC ! solution at PTSIMUL + PDTACT
INTEGER, INTENT(IN) :: KMI      ! model number
REAL, INTENT(IN) :: PALPHA
!
!*       0.2  declaration of local variables
!
INTEGER :: IFAIL, JI, JJ
INTEGER :: IITERCOUNT      ! counter for Newton-Raphson iteration
INTEGER :: IMAXITER = 20   ! maximal # of iterations before faillure
REAL, DIMENSION(KVECNPT) :: ZERR            ! error in iteration
REAL    :: ZMAXERR = 1E-6  ! maximal relative error for solution of implicit eqn
REAL    :: ZNORM           ! vector norm for calculation of relative error
REAL, DIMENSION(KVECNPT,KEQ) :: ZYN      ! stores y^n + dt**alpha*f^n
REAL, DIMENSION(KVECNPT,KEQ) :: ZF       ! f(y)
REAL, DIMENSION(KEQ) :: ZFTRAPEZ ! F(z^m), then z^m+1 - z^m
REAL, DIMENSION(KVECNPT,KEQ) :: ZFTRAPEZVECT ! F(z^m), then z^m+1 - z^m
REAL, DIMENSION(KEQ,KEQ) :: ZB,ZC   ! working matrice for the iteration
REAL, DIMENSION(KVECNPT,KEQ,KEQ) :: ZBVECT 

!------------------------------------------------------------------------------
!
!*       1.   CALCULATE FIRST GUESS FOR ITERATION (ZYN)
!        ----------------------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('CH_CRANCK',0,ZHOOK_HANDLE)
CALL CH_FCN(PTSIMUL,PCONC,ZF,KMI,KVECNPT,KEQ)
ZYN(:,:) = PCONC(:,:) + PALPHA*PDTACT*ZF(:,:)
PNEWCONC(:,:) = PCONC(:,:)
!
!*       2.   NEWTON RAPHSON ITERATION
!        -----------------------------
!
ZERR(:) = 2.*ZMAXERR
IITERCOUNT = 0
newton: DO WHILE (MAXVAL(ZERR).GT.ZMAXERR)
!
  IITERCOUNT = IITERCOUNT + 1
  IF (IITERCOUNT.GT.IMAXITER) THEN
    STOP "CH_CRANCK ERROR: no convergence of Newton-Raphson iteration obtained"
  ENDIF
!
!*       2.1  calculate derivative F for next iteration
!
  CALL CH_FCN(PTSIMUL+PDTACT,PNEWCONC,ZF,KMI,KVECNPT,KEQ)
  
  ZFTRAPEZVECT(:,:) = PNEWCONC(:,:) - ZYN(:,:)  - (1.0-PALPHA)*PDTACT*ZF(:,:)
!
!*       2.2  calculate the Jacobien B
!
  CALL CH_JAC(PTSIMUL+PDTACT,PNEWCONC,ZBVECT,KMI,KVECNPT,KEQ)
!
!*       2.3  modify JAC after Cranck-Nicholson method
!
  ZBVECT(:,:,:) = -(1.0-PALPHA)*PDTACT*ZBVECT(:,:,:)
  DO JI = 1, KEQ
    ZBVECT(:,JI,JI) = 1.0 + ZBVECT(:,JI,JI)
  ENDDO
!
!%kloos%
! Modification of matrix inversion (25/04/97)
!
!*       2.4 calculate LU factorization for ZB (result is put in ZB) 
!
  DO JI=1,KVECNPT
  
  ZB(:,:)=ZBVECT(JI,:,:)
  ZFTRAPEZ(:)=ZFTRAPEZVECT(JI,:)
  
  IFAIL = 1
  CALL CH_GAUSS(ZB,ZC,KEQ,IFAIL)
  IF (IFAIL.NE.0) THEN
    STOP 'CH_CRANCK ERROR: matrix cannot be inverted by CH_GAUSS'
  ENDIF
!
!*       2.5  calculate dY = ZB F (result is put in ZFTRAPEZ)
!
  ZFTRAPEZ(:)=MATMUL(ZC(:,:),ZFTRAPEZ(:))
!
!*       2.6  calculate Y (n+1) 
!
  ZERR(JI) = 0.0
  ZNORM = 0.0
  DO JJ=1,KEQ
    ZERR(JI) = ZERR(JI) + ABS(ZFTRAPEZ(JJ))
    PNEWCONC(JI,JJ) = PNEWCONC(JI,JJ) - ZFTRAPEZ(JJ)
    ZNORM = ZNORM + 0.5*ABS(PCONC(JI,JJ)+PNEWCONC(JI,JJ))
  ENDDO
!%
  ZERR(JI) = ZERR(JI) / ZNORM
  
  ENDDO
!
END DO newton
!
IF (LHOOK) CALL DR_HOOK('CH_CRANCK',1,ZHOOK_HANDLE)
END SUBROUTINE CH_CRANCK
