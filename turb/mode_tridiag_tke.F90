!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
MODULE MODE_TRIDIAG_TKE
IMPLICIT NONE
CONTAINS       
SUBROUTINE TRIDIAG_TKE(D,PVARM,PA,PTSTEP,PEXPL,PIMPL, &
                                  PRHODJ,PSOURCE,PDIAG,PVARP )
       USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
!      ########################################################
!
!
!!****   *TRIDIAG_TKE* - routine to solve a time implicit scheme
!!
!!
!!     PURPOSE
!!     -------
!        The purpose of this routine is to give a field PVARP at t+1, by 
!      solving an implicit tridiagonal system obtained by the 
!      discretization of the vertical turbulent diffusion. It should be noted 
!      that the degree of implicitness can be varied (PIMPL parameter) and the
!      sources of evolution other than the turbulent diffusion can be taken
!      into account through the PSOURCE field. PVARP is localized at a mass 
!      point.
!
!!**   METHOD
!!     ------
!!        First, the Right Hand Side of the implicit equation is computed.
!!     It is build as follows:
!!        ZY = PVARM + PTSTEP*PSOURCE + DIFF_EXPLI
!!     where PVARM is the variable at t-dt, PSOURCE the supplementary sources of
!!     PVAR ( and not PVAR * PRHODJ !!) and  DIFF_EXPLI is the explicit part
!!     of the vertical turbulent diffusion. This operator is spatially 
!!     discretized as the implicit one, thus:
!!        DIFF_EXPLI(k) = - PEXPL / PRHODJ(k) * 
!!                       ( PA(k+1) * (PVARM(k+1) - PVARM(k)  )
!!                        -PA(k)   * (PVARM(k)   - PVARM(k-1)) )
!!     For the first level, only the upper part is considered, the lower one 
!!     is replaced by the turbulent surface flux (taken into account in the 
!!     PSOURCE(ikb) term).
!!        DIFF_EXPLI(ikb) = - PEXPL / PRHODJ(ikb) * 
!!                       ( PA(ikb+1) * (PVARM(ikb+1) - PVARM(ikb))  )
!!     For the last level, only the lower part is considered, the upper one 
!!     is replaced by the turbulent flux which is taken equal to 0 
!!     (taken into account in the PSOURCE(ike) term).
!!
!!        DIFF_EXPLI(ike) = + PEXPL / PRHODJ(ike) * 
!!                       ( PA(ike) * (PVARM(ike) - PVARM(ike-1))  )
!!                      
!!        Then, the classical tridiagonal algorithm is used to invert the 
!!     implicit operator. Its matrix is given by:
!!
!!     ( b(ikb)   c(ikb)      0        0        0         0        0        0  )
!!     (   0      a(ikb+1) b(ikb+1) c(ikb+1)    0  ...    0        0        0  ) 
!!     (   0         0     a(ikb+2) b(ikb+2) c(ikb+2).    0        0        0  ) 
!!      .......................................................................
!!     (   0   ...   0     a(k)     b(k)     c(k)         0   ...  0        0  ) 
!!      .......................................................................
!!     (   0         0        0        0        0 ...a(ike-1) b(ike-1) c(ike-1))
!!     (   0         0        0        0        0 ...     0   a(ike)   b(ike)  )
!!
!!     ikb and ike represent the first and the last inner mass levels of the
!!     model. The coefficients are:
!!         
!!          a(k) = PIMPL * PA(k)/PRHODJ(k)   
!!          b(k) = 1 - PIMPL * PA(k)/PRHODJ(k) - PIMPL * PA(k+1)/PRHODJ(k)
!!          c(k) = PIMPL * PA(k+1)/PRHODJ(k)
!!
!!          for all k /= ikb or ike
!!
!!          b(ikb) = 1 - PIMPL * PA(ikb+1)/PRHODJ(ikb)
!!          c(ikb) = PIMPL * PA(ikb+1)/PRHODJ(ikb)
!!             (discretization of the upper part of the implicit operator)
!!          b(ike) = 1 - PIMPL * PA(ike)/PRHODJ(ike)
!!          a(ike) = PIMPL * PA(ike)/PRHODJ(ike)
!!             (discretization of the lower part of the implicit operator)
!!       Finally, the marginal points are prescribed.
!!
!!       All these computations are purely vertical and vectorizations are 
!!     easely achieved by processing all the verticals in parallel.
!!
!!     EXTERNAL
!!     --------
!!
!!       NONE
!!
!!     IMPLICIT ARGUMENTS
!!     ------------------
!!       MODD_PARAMETERS
!!            JPVEXT_TURB: number of vertical external points
!!
!!     REFERENCE
!!     ---------
!!       Book 1 of Meso-NH documentation (chapter Turbulence)
!!       Press et al: Numerical recipes (1986) Cambridge Univ. Press
!!
!!     AUTHOR
!!     ------
!!       Joan Cuxart       * INM and Meteo-France *   
!! 
!!     MODIFICATIONS
!!     -------------
!!       Original        August 29, 1994
!!       Modification :  January 29, 1995 Algorithm written with two 
!!                                        local variables less
!!      (Cuxart, Stein)  August 21, 1995  Bug correction for PRHODJ
!!      (Stein)         November 16, 1995 new version
!!      (Stein)         February 28, 1995 no inversion in the explicit case
!! ---------------------------------------------------------------------
!
!*       0. DECLARATIONS
!
USE MODD_PARAMETERS
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
!
IMPLICIT NONE
!
!
!*       0.1 declarations of arguments
!
TYPE(DIMPHYEX_t),                 INTENT(IN)   :: D
REAL, DIMENSION(D%NIJT,D%NKT),    INTENT(IN)  :: PVARM       ! variable at t-1  
REAL, DIMENSION(D%NIJT,D%NKT),    INTENT(IN)  :: PA          ! upper diag. elements
REAL,                             INTENT(IN)  :: PTSTEP      ! Double time step
REAL,                             INTENT(IN)  :: PEXPL,PIMPL ! weights of the temporal scheme
REAL, DIMENSION(D%NIJT,D%NKT),    INTENT(IN)  :: PRHODJ      ! (dry rho)*J
REAL, DIMENSION(D%NIJT,D%NKT),    INTENT(IN)  :: PSOURCE     ! source term of PVAR    
REAL, DIMENSION(D%NIJT,D%NKT),    INTENT(IN)  :: PDIAG       ! diagonal term linked to
                                                      ! the implicit dissipation
!
REAL, DIMENSION(D%NIJT,D%NKT),    INTENT(OUT) :: PVARP       ! variable at t+1        
!
!*       0.2 declarations of local variables
!
REAL, DIMENSION(D%NIJT,D%NKT)  :: ZY ,ZGAM 
                                         ! RHS of the equation, 3D work array
REAL, DIMENSION(D%NIJT)        :: ZBET
                                         ! 2D work array
INTEGER             :: JIJ,JK     ! loop counter
INTEGER             :: IKB,IKE      ! inner vertical limits
INTEGER             :: IKT,IKA,IKU  ! array size in k direction
INTEGER             :: IKTB,IKTE    ! start, end of k loops in physical domain
INTEGER             :: IIJB, IIJE   ! start, end of ij loops in physical domain
INTEGER             :: IKL
!
! ---------------------------------------------------------------------------
!                                              
!*      1.  COMPUTE THE RIGHT HAND SIDE
!           ---------------------------
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('TRIDIAG_TKE',0,ZHOOK_HANDLE)
!


IKT=D%NKT
IKTB=D%NKTB
IKTE=D%NKTE
IKB=D%NKB
IKE=D%NKE
IKA=D%NKA
IKU=D%NKU
IKL=D%NKL
IIJB=D%NIJB
IIJE=D%NIJE
!

DO JIJ=IIJB, IIJE
  ZY(JIJ, IKB) = PVARM(JIJ, IKB)  + PTSTEP*PSOURCE(JIJ, IKB) -   &
    PEXPL / PRHODJ(JIJ, IKB) * PA(JIJ, IKB+IKL) *  & 
    (PVARM(JIJ, IKB+IKL) - PVARM(JIJ, IKB))
END DO

!


DO JK=IKTB+1,IKTE-1
  DO JIJ=IIJB, IIJE
    ZY(JIJ, JK)= PVARM(JIJ, JK)  + PTSTEP*PSOURCE(JIJ, JK) -               &
        PEXPL / PRHODJ(JIJ, JK) *                                           &
                               ( PVARM(JIJ, JK-IKL)*PA(JIJ, JK)                &
                                -PVARM(JIJ, JK)*(PA(JIJ, JK)+PA(JIJ, JK+IKL))   &
                                +PVARM(JIJ, JK+IKL)*PA(JIJ, JK+IKL)              &
                               ) 
  END DO
END DO

! 

DO JIJ=IIJB, IIJE
  ZY(JIJ, IKE)= PVARM(JIJ, IKE) + PTSTEP*PSOURCE(JIJ, IKE) +               &
    PEXPL / PRHODJ(JIJ, IKE) * PA(JIJ, IKE) * (PVARM(JIJ, IKE)-PVARM(JIJ, IKE-IKL))
END DO

!
!
!*       2.  INVERSION OF THE TRIDIAGONAL SYSTEM
!            -----------------------------------
!
IF ( PIMPL > 1.E-10 ) THEN

  !
  !  going up
  !
  DO JIJ=IIJB, IIJE
    ZBET(JIJ) = 1. + PIMPL * (PDIAG(JIJ, IKB)-PA(JIJ, IKB+IKL) / PRHODJ(JIJ, IKB))
                                                      ! bet = b(ikb)
    PVARP(JIJ, IKB) = ZY(JIJ, IKB) / ZBET(JIJ)                
  END DO
  !



  DO JK = IKB+IKL,IKE-IKL,IKL
    DO JIJ=IIJB, IIJE
      ZGAM(JIJ, JK) = PIMPL * PA(JIJ, JK) / PRHODJ(JIJ, JK-IKL) / ZBET(JIJ)  
                                                      ! gam(k) = c(k-1) / bet
      ZBET(JIJ)    = 1. + PIMPL * ( PDIAG(JIJ, JK) -                     &
                                   (  PA(JIJ, JK) * (1. + ZGAM(JIJ, JK))  &
                                    + PA(JIJ, JK+IKL)                      &
                                   ) / PRHODJ(JIJ, JK)                   &
                                  )                   ! bet = b(k) - a(k)* gam(k)  
      PVARP(JIJ, JK)= ( ZY(JIJ, JK) - PIMPL * PA(JIJ, JK) / PRHODJ(JIJ, JK) &
                      * PVARP(JIJ, JK-IKL)                                 &
                     ) / ZBET(JIJ)
                                          ! res(k) = (y(k) -a(k)*res(k-1))/ bet 
    END DO
  END DO


  DO JIJ=IIJB, IIJE
    ! special treatment for the last level
    ZGAM(JIJ, IKE) = PIMPL * PA(JIJ, IKE) / PRHODJ(JIJ, IKE-IKL) / ZBET(JIJ) 
                                                      ! gam(k) = c(k-1) / bet
    ZBET(JIJ)    = 1. + PIMPL * ( PDIAG(JIJ, IKE) -                   &
           (  PA(JIJ, IKE) * (1. + ZGAM(JIJ, IKE)) ) / PRHODJ(JIJ, IKE) &
                                )  
                                                      ! bet = b(k) - a(k)* gam(k)  
    PVARP(JIJ, IKE)= ( ZY(JIJ, IKE) - PIMPL * PA(JIJ, IKE) / PRHODJ(JIJ, IKE) &
                                  * PVARP(JIJ, IKE-IKL)                      &
                   ) / ZBET(JIJ)
                                         ! res(k) = (y(k) -a(k)*res(k-1))/ bet 
  END DO

  !
  !  going down
  !


  DO JK = IKE-IKL,IKB,-1*IKL
  DO JIJ=IIJB, IIJE
      PVARP(JIJ, JK) = PVARP(JIJ, JK) - ZGAM(JIJ, JK+IKL) * PVARP(JIJ, JK+IKL) 
  END DO
  END DO

!
ELSE
!

  DO JK=IKTB,IKTE
    DO JIJ=IIJB, IIJE
      PVARP(JIJ, JK) = ZY(JIJ, JK)
    END DO
  END DO

!
END IF 
!
!
!*       3.  FILL THE UPPER AND LOWER EXTERNAL VALUES
!            ----------------------------------------
!

DO JIJ=IIJB, IIJE
  PVARP(JIJ, IKA)=PVARP(JIJ, IKB)
  PVARP(JIJ, IKU)=PVARP(JIJ, IKE)
END DO

!

!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('TRIDIAG_TKE',1,ZHOOK_HANDLE)
END SUBROUTINE TRIDIAG_TKE
END MODULE MODE_TRIDIAG_TKE
