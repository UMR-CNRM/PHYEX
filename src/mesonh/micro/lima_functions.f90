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
!#################################
        MODULE MODI_LIMA_FUNCTIONS
!#################################
!
INTERFACE
!
FUNCTION MOMG (PALPHA,PNU,PP) RESULT (PMOMG)
  REAL,                        INTENT(IN)           :: PALPHA 
  REAL,                        INTENT(IN)           :: PNU    
  REAL,                        INTENT(IN)           :: PP     
  REAL                                              :: PMOMG  
END FUNCTION MOMG
!
FUNCTION RECT(PA,PB,PX,PX1,PX2)  RESULT(PRECT)
  REAL,                        INTENT(IN)           :: PA
  REAL,                        INTENT(IN)           :: PB
  REAL, DIMENSION(:),          INTENT(IN)           :: PX
  REAL,                        INTENT(IN)           :: PX1
  REAL,                        INTENT(IN)           :: PX2
  REAL, DIMENSION(SIZE(PX,1))                       :: PRECT
END FUNCTION RECT
!
FUNCTION DELTA(PA,PB,PX,PX1,PX2)  RESULT(PDELTA)
  REAL,                        INTENT(IN)           :: PA
  REAL,                        INTENT(IN)           :: PB
  REAL, DIMENSION(:),          INTENT(IN)           :: PX
  REAL,                        INTENT(IN)           :: PX1
  REAL,                        INTENT(IN)           :: PX2
  REAL, DIMENSION(SIZE(PX,1))                       :: PDELTA
END FUNCTION DELTA
!
FUNCTION DELTA_VEC(PA,PB,PX,PX1,PX2)  RESULT(PDELTA_VEC)
  REAL,                        INTENT(IN)           :: PA
  REAL,                        INTENT(IN)           :: PB
  REAL, DIMENSION(:),          INTENT(IN)           :: PX
  REAL, DIMENSION(:),          INTENT(IN)           :: PX1
  REAL, DIMENSION(:),          INTENT(IN)           :: PX2
  REAL, DIMENSION(SIZE(PX,1))                       :: PDELTA_VEC
END FUNCTION DELTA_VEC
!
SUBROUTINE GAULAG(x,w,n,alf)
  INTEGER,                     INTENT(IN)           :: n
  REAL,                        INTENT(IN)           :: alf
  REAL, DIMENSION(n),          INTENT(INOUT)        :: w, x
END SUBROUTINE GAULAG
!
SUBROUTINE GAUHER(x,w,n)
  INTEGER,                     INTENT(IN)           :: n
  REAL, DIMENSION(n),          INTENT(INOUT)        :: w, x
END SUBROUTINE GAUHER
!
END INTERFACE
!
END MODULE MODI_LIMA_FUNCTIONS
!
!------------------------------------------------------------------------------
!
!###########################################
FUNCTION MOMG (PALPHA,PNU,PP) RESULT (PMOMG)
!###########################################
!
! auxiliary routine used to compute the Pth moment order of the generalized
! gamma law
!
  USE MODI_GAMMA
!
  IMPLICIT NONE
!
  REAL     :: PALPHA ! first shape parameter of the dimensionnal distribution
  REAL     :: PNU    ! second shape parameter of the dimensionnal distribution
  REAL     :: PP     ! order of the moment
  REAL     :: PMOMG  ! result: moment of order ZP
!
  PMOMG = GAMMA_X0D(PNU+PP/PALPHA)/GAMMA_X0D(PNU)
!
END FUNCTION MOMG
!
!------------------------------------------------------------------------------
!
!#############################################
FUNCTION RECT(PA,PB,PX,PX1,PX2)  RESULT(PRECT)
!#############################################
!
!     PRECT takes the value PA if PX1<=PX<PX2, and PB outside the [PX1;PX2[ interval
!
  IMPLICIT NONE
!
  REAL,                        INTENT(IN)           :: PA
  REAL,                        INTENT(IN)           :: PB
  REAL, DIMENSION(:),          INTENT(IN)           :: PX
  REAL,                        INTENT(IN)           :: PX1
  REAL,                        INTENT(IN)           :: PX2
  REAL, DIMENSION(SIZE(PX,1))                       :: PRECT
!
  PRECT(:) = PB
  WHERE (PX(:).GE.PX1 .AND. PX(:).LT.PX2)
     PRECT(:) = PA
  END WHERE
  RETURN
!
END FUNCTION RECT
!
!-------------------------------------------------------------------------------
!
!###############################################
FUNCTION DELTA(PA,PB,PX,PX1,PX2)  RESULT(PDELTA)
!###############################################
!
!     PDELTA takes the value PA if PX<PX1, and PB if PX>=PX2
!     PDELTA is a cubic interpolation between PA and PB for PX between PX1 and PX2 
!
  IMPLICIT NONE
!
  REAL,                        INTENT(IN)           :: PA
  REAL,                        INTENT(IN)           :: PB
  REAL, DIMENSION(:),          INTENT(IN)           :: PX
  REAL,                        INTENT(IN)           :: PX1
  REAL,                        INTENT(IN)           :: PX2
  REAL, DIMENSION(SIZE(PX,1))                       :: PDELTA
!
!*       local variable
!
  REAL                                              :: ZA
!
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
!#######################################################
FUNCTION DELTA_VEC(PA,PB,PX,PX1,PX2)  RESULT(PDELTA_VEC)
!#######################################################
!
!    Same as DELTA for vectorized PX1 and PX2 arguments
!
  IMPLICIT NONE
!
  REAL,                        INTENT(IN)           :: PA
  REAL,                        INTENT(IN)           :: PB
  REAL, DIMENSION(:),          INTENT(IN)           :: PX
  REAL, DIMENSION(:),          INTENT(IN)           :: PX1
  REAL, DIMENSION(:),          INTENT(IN)           :: PX2
  REAL, DIMENSION(SIZE(PX,1))                       :: PDELTA_VEC
!
!*       local variable
!
  REAL, DIMENSION(SIZE(PX,1))                       :: ZA
!
  ZA(:) = 0.0
  wHERE     (PX(:)<=PX1(:))
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
!
END FUNCTION DELTA_VEC
!
!-------------------------------------------------------------------------------
!
!###########################
SUBROUTINE gaulag(x,w,n,alf)
!###########################
  use modd_precision, only: MNHREAL64

  INTEGER n,MAXIT
  REAL alf,w(n),x(n)
  REAL(kind=MNHREAL64) :: EPS
  PARAMETER (EPS=3.D-14,MAXIT=10)
  INTEGER i,its,j
  REAL ai
  REAL(kind=MNHREAL64) :: p1,p2,p3,pp,z,z1
!
  REAL SUMW
!
  do 13 i=1,n
     if(i.eq.1)then
        z=(1.+alf)*(3.+.92*alf)/(1.+2.4*n+1.8*alf)
     else if(i.eq.2)then
        z=z+(15.+6.25*alf)/(1.+.9*alf+2.5*n)
     else
        ai=i-2
        z=z+((1.+2.55*ai)/(1.9*ai)+1.26*ai*alf/(1.+3.5*ai))* &
             (z-x(i-2))/(1.+.3*alf)
     endif
     do 12 its=1,MAXIT
        p1=1.d0
        p2=0.d0
        do 11 j=1,n
           p3=p2
           p2=p1
           p1=((2*j-1+alf-z)*p2-(j-1+alf)*p3)/j
11      continue
        pp=(n*p1-(n+alf)*p2)/z
        z1=z
        z=z1-p1/pp
        if(abs(z-z1).le.EPS)goto 1
12   continue
1    x(i)=z
     w(i)=-exp(gammln(alf+n)-gammln(real(n)))/(pp*n*p2)
13 continue
!
! NORMALISATION
!
  SUMW = 0.0
  DO 14 I=1,N
     SUMW = SUMW + W(I)
14 CONTINUE
  DO 15 I=1,N
     W(I) = W(I)/SUMW
15 CONTINUE
!
  return
END SUBROUTINE gaulag
!
!------------------------------------------------------------------------------
!
!##########################################
SUBROUTINE gauher(x,w,n)
!##########################################
  use modd_precision, only: MNHREAL64

  INTEGER n,MAXIT
  REAL w(n),x(n)
  REAL(kind=MNHREAL64) :: EPS,PIM4
  PARAMETER (EPS=3.D-14,PIM4=.7511255444649425D0,MAXIT=10)
  INTEGER i,its,j,m
  REAL(kind=MNHREAL64) :: p1,p2,p3,pp,z,z1
!
  REAL SUMW
!
  m=(n+1)/2
  do 13 i=1,m
     if(i.eq.1)then
        z=sqrt(real(2*n+1))-1.85575*(2*n+1)**(-.16667)
     else if(i.eq.2)then
        z=z-1.14*n**.426/z
     else if (i.eq.3)then
        z=1.86*z-.86*x(1)
     else if (i.eq.4)then
        z=1.91*z-.91*x(2)
     else
        z=2.*z-x(i-2)
     endif
     do 12 its=1,MAXIT
        p1=PIM4
        p2=0.d0
        do 11 j=1,n
           p3=p2
           p2=p1
           p1=z*sqrt(2.d0/j)*p2-sqrt(dble(j-1)/dble(j))*p3
11      continue
        pp=sqrt(2.d0*n)*p2
        z1=z
        z=z1-p1/pp
        if(abs(z-z1).le.EPS)goto 1
12   continue
1    x(i)=z
     x(n+1-i)=-z
     pp=pp/PIM4 ! NORMALIZATION
     w(i)=2.0/(pp*pp)
     w(n+1-i)=w(i)
13 continue
!
! NORMALISATION
!
  SUMW = 0.0
  DO 14 I=1,N
     SUMW = SUMW + W(I)
14 CONTINUE
  DO 15 I=1,N
     W(I) = W(I)/SUMW
15 CONTINUE
!
  return
END SUBROUTINE gauher
!
!------------------------------------------------------------------------------
