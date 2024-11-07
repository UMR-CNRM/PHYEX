!MNH_LIC Copyright 2016-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
! Modifications:
!  P. Wautelet 22/01/2019: replace double precision declarations by real(kind(0.0d0)) (to allow compilation by NAG compiler)
!  P. Wautelet 19/04/2019: use modd_precision kinds
!  P. Wautelet 26/04/2019: replace non-standard FLOAT function by REAL function
!  P. Wautelet 28/05/2019: move COUNTJV function to tools.f90
!-----------------------------------------------------------------
MODULE MODE_LIMA_FUNCTIONS
  IMPLICIT NONE
CONTAINS
!
!------------------------------------------------------------------------------
!
  FUNCTION MOMG (PALPHA,PNU,PP) RESULT (PMOMG)
! Pth moment order of the generalized gamma law
    USE MODI_GAMMA
    IMPLICIT NONE
    REAL, INTENT(IN)     :: PALPHA ! first shape parameter of the dimensionnal distribution
    REAL, INTENT(IN)     :: PNU    ! second shape parameter of the dimensionnal distribution
    REAL, INTENT(IN)     :: PP     ! order of the moment
    REAL     :: PMOMG  ! result: moment of order ZP
    PMOMG = GAMMA_X0D(PNU+PP/PALPHA)/GAMMA_X0D(PNU)
  END FUNCTION MOMG
!
!------------------------------------------------------------------------------
!
  FUNCTION RECT(PA,PB,PX,PX1,PX2)  RESULT(PRECT)
!     PRECT takes the value PA if PX1<=PX<PX2, and PB outside the [PX1;PX2[ interval
    IMPLICIT NONE
    REAL,                        INTENT(IN) :: PA
    REAL,                        INTENT(IN) :: PB
    REAL, DIMENSION(:),          INTENT(IN) :: PX
    REAL,                        INTENT(IN) :: PX1
    REAL,                        INTENT(IN) :: PX2
    REAL, DIMENSION(SIZE(PX,1))             :: PRECT
    PRECT(:) = PB
    WHERE (PX(:).GE.PX1 .AND. PX(:).LT.PX2)
       PRECT(:) = PA
    END WHERE
    RETURN
  END FUNCTION RECT
!
!-------------------------------------------------------------------------------
!
  FUNCTION DELTA(PA,PB,PX,PX1,PX2)  RESULT(PDELTA)
!     PDELTA takes the value PA if PX<PX1, and PB if PX>=PX2
!     PDELTA is a cubic interpolation between PA and PB for PX between PX1 and PX2 
    IMPLICIT NONE
    REAL,                        INTENT(IN) :: PA
    REAL,                        INTENT(IN) :: PB
    REAL, DIMENSION(:),          INTENT(IN) :: PX
    REAL,                        INTENT(IN) :: PX1
    REAL,                        INTENT(IN) :: PX2
    REAL, DIMENSION(SIZE(PX,1))             :: PDELTA
    REAL                                    :: ZA
    ZA = 6.0*(PA-PB)/(PX2-PX1)**3
    WHERE     (PX(:).LT.PX1)
       PDELTA(:) = PA
    ELSEWHERE (PX(:).GE.PX2)
       PDELTA(:) = PB
    ELSEWHERE
       PDELTA(:) =   PA + ZA*PX1**2*(PX1/6.0 - 0.5*PX2)              &
            + ZA*PX1*PX2*                          (PX(:))    &
            - (0.5*ZA*(PX1+PX2))*                  (PX(:)**2) & 
            + (ZA/3.0)*                            (PX(:)**3)
    END WHERE
    RETURN
!
  END FUNCTION DELTA
!
!-------------------------------------------------------------------------------
!
  FUNCTION DELTA_VEC(PA,PB,PX,PX1,PX2)  RESULT(PDELTA_VEC)
!    Same as DELTA for vectorized PX1 and PX2 arguments
    IMPLICIT NONE
    REAL,                        INTENT(IN) :: PA
    REAL,                        INTENT(IN) :: PB
    REAL, DIMENSION(:),          INTENT(IN) :: PX
    REAL, DIMENSION(:),          INTENT(IN) :: PX1
    REAL, DIMENSION(:),          INTENT(IN) :: PX2
    REAL, DIMENSION(SIZE(PX,1))             :: PDELTA_VEC
    REAL, DIMENSION(SIZE(PX,1))             :: ZA
    ZA(:) = 0.0
    WHERE     (PX(:)<=PX1(:))
       PDELTA_VEC(:) = PA
    ELSEWHERE (PX(:)>=PX2(:))
       PDELTA_VEC(:) = PB
    ELSEWHERE
       ZA(:)         = 6.0*(PA-PB)/(PX2(:)-PX1(:))**3
       PDELTA_VEC(:) =   PA + ZA(:)*PX1(:)**2*(PX1(:)/6.0 - 0.5*PX2(:))           &
            + ZA(:)*PX1(:)*PX2(:)*                          (PX(:))    &
            - (0.5*ZA(:)*(PX1(:)+PX2(:)))*                  (PX(:)**2) & 
            + (ZA(:)/3.0)*                                  (PX(:)**3)
    END WHERE
    RETURN
  END FUNCTION DELTA_VEC
!
!-------------------------------------------------------------------------------
!
SUBROUTINE GAULAG(X,W,N,ALF)
  USE MODD_PRECISION, only: MNHREAL64
  INTEGER, INTENT(IN) :: N
  INTEGER MAXIT
  REAL, INTENT(IN)  :: ALF
  REAL, INTENT(OUT) :: W(N),X(N)
  REAL(KIND=MNHREAL64) :: EPS
  PARAMETER (EPS=3.D-14,MAXIT=10)
  INTEGER I,ITS,J
  REAL AI
  REAL(KIND=MNHREAL64) :: P1,P2,P3,PP,Z,Z1
!
  REAL SUMW
!
  DO 13 I=1,N
     IF(I.EQ.1)then
        Z=(1.+ALF)*(3.+.92*ALF)/(1.+2.4*N+1.8*ALF)
     ELSE IF(I.EQ.2)then
        Z=Z+(15.+6.25*ALF)/(1.+.9*ALF+2.5*N)
     ELSE
        AI=I-2
        Z=Z+((1.+2.55*AI)/(1.9*AI)+1.26*AI*ALF/(1.+3.5*AI))* &
             (Z-X(I-2))/(1.+.3*ALF)
     ENDIF
     DO 12 ITS=1,MAXIT
        P1=1.D0
        P2=0.D0
        DO 11 J=1,N
           P3=P2
           P2=P1
           P1=((2*J-1+ALF-Z)*P2-(J-1+ALF)*P3)/J
11      CONTINUE
        PP=(N*P1-(N+ALF)*P2)/Z
        Z1=Z
        Z=Z1-P1/PP
        IF(ABS(Z-Z1).LE.EPS)GOTO 1
12   CONTINUE
1    X(I)=Z
     W(I)=-EXP(GAMMLN(ALF+N)-GAMMLN(REAL(N)))/(PP*N*P2)
13 CONTINUE
! NORMALISATION
  SUMW = 0.0
  DO 14 I=1,N
     SUMW = SUMW + W(I)
14 CONTINUE
  DO 15 I=1,N
     W(I) = W(I)/SUMW
15 CONTINUE
!
  RETURN
END SUBROUTINE GAULAG
!
!------------------------------------------------------------------------------
!
SUBROUTINE GAUHER(X,W,N)
  USE MODD_PRECISION, only: MNHREAL64
  INTEGER, INTENT(IN) :: N
  INTEGER MAXIT
  REAL, INTENT(OUT) ::  W(N),X(N)
  REAL(KIND=MNHREAL64) :: EPS,PIM4
  PARAMETER (EPS=3.D-14,PIM4=.7511255444649425D0,MAXIT=10)
  INTEGER I,ITS,J,M
  REAL(KIND=MNHREAL64) :: P1,P2,P3,PP,Z,Z1
!
  REAL SUMW
!
  M=(N+1)/2
  DO 13 I=1,M
     IF(I.EQ.1)then
        Z=SQRT(REAL(2*N+1))-1.85575*(2*N+1)**(-.16667)
     ELSE IF(I.EQ.2)then
        Z=Z-1.14*N**.426/Z
     ELSE IF (I.EQ.3)then
        Z=1.86*Z-.86*X(1)
     ELSE IF (I.EQ.4)then
        Z=1.91*Z-.91*X(2)
     ELSE
        Z=2.*Z-X(I-2)
     ENDIF
     DO 12 ITS=1,MAXIT
        P1=PIM4
        P2=0.D0
        DO 11 J=1,N
           P3=P2
           P2=P1
           P1=Z*SQRT(2.D0/J)*P2-SQRT(DBLE(J-1)/DBLE(J))*P3
11      CONTINUE
        PP=SQRT(2.D0*N)*P2
        Z1=Z
        Z=Z1-P1/PP
        IF(ABS(Z-Z1).LE.EPS)GOTO 1
12   CONTINUE
1    X(I)=Z
     X(N+1-I)=-Z
     PP=PP/PIM4 ! NORMALIZATION
     W(I)=2.0/(PP*PP)
     W(N+1-I)=W(I)
13 CONTINUE
! NORMALISATION
  SUMW = 0.0
  DO 14 I=1,N
     SUMW = SUMW + W(I)
14 CONTINUE
  DO 15 I=1,N
     W(I) = W(I)/SUMW
15 CONTINUE
!
  RETURN
END SUBROUTINE GAUHER
!
!------------------------------------------------------------------------------
!
FUNCTION ARTH(FIRST,INCREMENT,N)
  REAL,INTENT(IN) :: FIRST,INCREMENT
  INTEGER,INTENT(IN) :: N
  REAL,DIMENSION(N) :: ARTH
  INTEGER :: K
  DO K=1,N
     ARTH(K)=FIRST+INCREMENT*(K-1)
  END DO
END FUNCTION ARTH
!
!------------------------------------------------------------------------------
!
FUNCTION GAMMLN(XX)
  IMPLICIT NONE
  REAL, INTENT(IN) :: XX
  REAL :: GAMMLN
  REAL :: TMP,X
  REAL :: STP = 2.5066282746310005
  REAL, DIMENSION(6) :: COEF = (/76.18009172947146,& 
       -86.50532032941677,24.01409824083091,& 
       -1.231739572450155,0.1208650973866179E-2,&
       -0.5395239384953E-5/)
  X=XX
  TMP=X+5.5
  TMP=(X+0.5)*LOG(TMP)-TMP
  GAMMLN=TMP+LOG(STP*(1.000000000190015+&
       SUM(COEF(:)/ARTH(X+1.,1.,SIZE(COEF))))/X)
END FUNCTION GAMMLN
!
!------------------------------------------------------------------------------
!
END MODULE MODE_LIMA_FUNCTIONS
