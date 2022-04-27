!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
MODULE MODE_TRIDIAG_WIND
IMPLICIT NONE
CONTAINS       
SUBROUTINE TRIDIAG_WIND(D,PVARM,PA,PCOEFS,PTSTEP,PEXPL,PIMPL, &
                                             PRHODJA,PSOURCE,PVARP )
       USE PARKIND1, ONLY : JPRB
       USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!      #############################################################
!
!
!!****   *TRIDIAG_WIND* - routine to solve a time implicit scheme
!!
!!
!!     PURPOSE
!!     -------
!        The purpose of this routine is to give a field PVARP at t+1, by 
!      solving an implicit tridiagonal system obtained by the 
!      discretization of the vertical turbulent diffusion. It should be noted 
!      that the degree of implicitness can be varied (PIMPL parameter) and the
!      sources of evolution other than the turbulent diffusion can be taken
!      into account through the PSOURCE field. PVARP is localized at a wind 
!      point either U or V, PRHODJA is averaged to be localized at the same 
!      point. The surface flux is also implicitly computed. 
!
!!**   METHOD
!!     ------
!!        First, the Right Hand Side of the implicit equation is computed.
!!     It is build as follows:
!!        ZY = PVARM + PTSTEP*PSOURCE + DIFF_EXPLI
!!     where PVARM is the variable at t-dt, PSOURCE the supplementary sources of
!!     PVAR ( and not PVAR * PRHODJA !!) and  DIFF_EXPLI is the explicit part
!!     of the vertical turbulent diffusion. This operator is spatially 
!!     discretized as the implicit one, thus:
!!        DIFF_EXPLI(k) = - PEXPL / PRHODJA(k) * 
!!                       ( PA(k+1) * (PVARM(k+1) - PVARM(k)  )
!!                        -PA(k)   * (PVARM(k)   - PVARM(k-1)) )
!!     For the first level, only the upper part is considered, the lower one 
!!     is replaced by the turbulent surface flux (taken into account in the 
!!     PSOURCE(ikb) term).
!!        DIFF_EXPLI(ikb) = - PEXPL / PRHODJA(ikb) * 
!!                       ( PA(ikb+1) * (PVARM(ikb+1) - PVARM(ikb))  )
!!     For the last level, only the lower part is considered, the upper one 
!!     is replaced by the turbulent flux which is taken equal to 0 
!!     (taken into account in the PSOURCE(ike) term).
!!
!!        DIFF_EXPLI(ike) = + PEXPL / PRHODJA(ike) * 
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
!!          a(k) = PIMPL * PA(k)/PRHODJA(k)   
!!          b(k) = 1 - PIMPL * PA(k)/PRHODJA(k) - PIMPL * PA(k+1)/PRHODJA(k)
!!          c(k) = PIMPL * PA(k+1)/PRHODJA(k)
!!
!!          for all k /= ikb or ike
!!
!!          b(ikb) = 1 - PIMPL * PA(ikb+1)/PRHODJA(ikb) - PIMPL * PCOEFS
!!          c(ikb) = PIMPL * PA(ikb+1)/PRHODJA(ikb)
!!             (discretization of the upper part of the implicit operator)
!!          b(ike) = 1 - PIMPL * PA(ike)/PRHODJA(ike)
!!          a(ike) = PIMPL * PA(ike)/PRHODJA(ike)
!!             (discretization of the lower part of the implicit operator)
!!
!!          The surface flux is given by:
!!             <w'u'> = <w'u'>EXPL + PIMPL * PCOEFS * PVARP    
!!          The explicit part is taken into account in PSOURCE(ikb) and the 
!!          implicit one is present in the LHS of the equation in b(ikb)
!!       
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
!!       Joel Stein       * Meteo-France *   
!! 
!!     MODIFICATIONS
!!     -------------
!!       Original         November 16, 1995
!!      (Stein)           February 28, 1995 no inversion in the explicit case
!!      (Seity)           February 2012 add possibility to run with reversed
!!                            vertical levels
!! ---------------------------------------------------------------------
!
!*       0. DECLARATIONS
!
USE MODD_PARAMETERS
USE MODD_DIMPHYEX,   ONLY: DIMPHYEX_t
!
IMPLICIT NONE
!
!
!*       0.1 declarations of arguments
!
TYPE(DIMPHYEX_t),     INTENT(IN)   :: D
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),    INTENT(IN)  :: PVARM       ! variable at t-1  
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),    INTENT(IN)  :: PA          ! upper diag. elements
REAL, DIMENSION(D%NIT,D%NJT),      INTENT(IN)  :: PCOEFS      ! implicit coeff for the
                                                      ! surface flux
REAL,                      INTENT(IN)  :: PTSTEP      ! Double time step
REAL,                      INTENT(IN)  :: PEXPL,PIMPL ! weights of the temporal scheme
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),    INTENT(IN)  :: PRHODJA     ! (dry rho)*J averaged 
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),    INTENT(IN)  :: PSOURCE     ! source term of PVAR    
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),    INTENT(OUT) :: PVARP       ! variable at t+1        
!
!*       0.2 declarations of local variables
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT)  :: ZY ,ZGAM 
                                         ! RHS of the equation, 3D work array
REAL, DIMENSION(D%NIT,D%NJT)                :: ZBET
                                         ! 2D work array
INTEGER             :: JI,JJ,JK     ! loop counter
INTEGER             :: IKB,IKE      ! inner vertical limits
INTEGER             :: IKT          ! array size in k direction
INTEGER             :: IKTB,IKTE    ! start, end of k loops in physical domain 
INTEGER             :: IIE,IIB,IJE,IJB
!
! ---------------------------------------------------------------------------
!                                              
!*      1.  COMPUTE THE RIGHT HAND SIDE
!           ---------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('TRIDIAG_WIND',0,ZHOOK_HANDLE)
!
IKT=D%NKT
IKTB=D%NKTB
IKTE=D%NKTE
IKB=D%NKB
IKE=D%NKE
IIE=D%NIEC
IIB=D%NIBC
IJE=D%NJEC
IJB=D%NJBC
!
!$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE)
ZY(IIB:IIE,IJB:IJE,IKB) = PVARM(IIB:IIE,IJB:IJE,IKB)  + PTSTEP*PSOURCE(IIB:IIE,IJB:IJE,IKB) -   &
  PEXPL / PRHODJA(IIB:IIE,IJB:IJE,IKB) * PA(IIB:IIE,IJB:IJE,IKB+D%NKL) * (PVARM(IIB:IIE,IJB:IJE,IKB+D%NKL) - PVARM(IIB:IIE,IJB:IJE,IKB))
!$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE)
!
DO JK=IKTB+1,IKTE-1
  !$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE)
  ZY(IIB:IIE,IJB:IJE,JK)= PVARM(IIB:IIE,IJB:IJE,JK)  + PTSTEP*PSOURCE(IIB:IIE,IJB:IJE,JK) -               &
      PEXPL / PRHODJA(IIB:IIE,IJB:IJE,JK) *                                          &
                             ( PVARM(IIB:IIE,IJB:IJE,JK-D%NKL)*PA(IIB:IIE,IJB:IJE,JK)                &
                              -PVARM(IIB:IIE,IJB:IJE,JK)*(PA(IIB:IIE,IJB:IJE,JK)+PA(IIB:IIE,IJB:IJE,JK+D%NKL))   &
                              +PVARM(IIB:IIE,IJB:IJE,JK+D%NKL)*PA(IIB:IIE,IJB:IJE,JK+D%NKL)              &
                             ) 
  !$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE)
END DO
! 
!$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE)
ZY(IIB:IIE,IJB:IJE,IKE)= PVARM(IIB:IIE,IJB:IJE,IKE) + PTSTEP*PSOURCE(IIB:IIE,IJB:IJE,IKE) +               &
  PEXPL / PRHODJA(IIB:IIE,IJB:IJE,IKE) * PA(IIB:IIE,IJB:IJE,IKE) * (PVARM(IIB:IIE,IJB:IJE,IKE)-PVARM(IIB:IIE,IJB:IJE,IKE-D%NKL))
!$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE)
!
!
!*       2.  INVERSION OF THE TRIDIAGONAL SYSTEM
!            -----------------------------------
!
IF ( PIMPL > 1.E-10 ) THEN
!
  !
  !  going up
  !
  !$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE)
  ZBET(IIB:IIE,IJB:IJE) = 1. - PIMPL * (  PA(IIB:IIE,IJB:IJE,IKB+D%NKL) / PRHODJA(IIB:IIE,IJB:IJE,IKB) &  
                            + PCOEFS(IIB:IIE,IJB:IJE) *  PTSTEP        )   ! bet = b(ikb)
  PVARP(IIB:IIE,IJB:IJE,IKB) = ZY(IIB:IIE,IJB:IJE,IKB) / ZBET(IIB:IIE,IJB:IJE)
  !$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE)               
  !
  DO JK = IKB+D%NKL,IKE-D%NKL,D%NKL
    !$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE)
    ZGAM(IIB:IIE,IJB:IJE,JK) = PIMPL * PA(IIB:IIE,IJB:IJE,JK) / PRHODJA(IIB:IIE,IJB:IJE,JK-D%NKL) / ZBET(IIB:IIE,IJB:IJE)  
                                                    ! gam(k) = c(k-1) / bet
    ZBET(IIB:IIE,IJB:IJE)    = 1. - PIMPL * (  PA(IIB:IIE,IJB:IJE,JK) * (1. + ZGAM(IIB:IIE,IJB:IJE,JK))  &
                                 + PA(IIB:IIE,IJB:IJE,JK+D%NKL)                      &
                                ) / PRHODJA(IIB:IIE,IJB:IJE,JK)  
                                                    ! bet = b(k) - a(k)* gam(k)  
    PVARP(IIB:IIE,IJB:IJE,JK)= ( ZY(IIB:IIE,IJB:IJE,JK) - PIMPL * PA(IIB:IIE,IJB:IJE,JK) / PRHODJA(IIB:IIE,IJB:IJE,JK) &
                    * PVARP(IIB:IIE,IJB:IJE,JK-D%NKL)                                 &
                   ) / ZBET(IIB:IIE,IJB:IJE)
                                        ! res(k) = (y(k) -a(k)*res(k-1))/ bet 
    !$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE)
  END DO
  !$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE)
  ! special treatment for the last level
  ZGAM(IIB:IIE,IJB:IJE,IKE) = PIMPL * PA(IIB:IIE,IJB:IJE,IKE) / PRHODJA(IIB:IIE,IJB:IJE,IKE-D%NKL) / ZBET(IIB:IIE,IJB:IJE) 
                                                    ! gam(k) = c(k-1) / bet
  ZBET(IIB:IIE,IJB:IJE)    = 1. - PIMPL * (  PA(IIB:IIE,IJB:IJE,IKE) * (1. + ZGAM(IIB:IIE,IJB:IJE,IKE))  &
                              ) / PRHODJA(IIB:IIE,IJB:IJE,IKE)  
                                                    ! bet = b(k) - a(k)* gam(k)  
  PVARP(IIB:IIE,IJB:IJE,IKE)= ( ZY(IIB:IIE,IJB:IJE,IKE) - PIMPL * PA(IIB:IIE,IJB:IJE,IKE) / PRHODJA(IIB:IIE,IJB:IJE,IKE) &
                                 * PVARP(IIB:IIE,IJB:IJE,IKE-D%NKL)                      &
                  ) / ZBET(IIB:IIE,IJB:IJE)
                                        ! res(k) = (y(k) -a(k)*res(k-1))/ bet 
  !
  !  going down
  !
  !$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE)
  DO JK = IKE-D%NKL,IKB,-1*D%NKL
  !$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE)
    PVARP(IIB:IIE,IJB:IJE,JK) = PVARP(IIB:IIE,IJB:IJE,JK) - ZGAM(IIB:IIE,IJB:IJE,JK+D%NKL) * PVARP(IIB:IIE,IJB:IJE,JK+D%NKL) 
  !$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE)
  END DO
!
ELSE
! 
  DO JK=IKTB,IKTE
    !$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE)
    PVARP(IIB:IIE,IJB:IJE,JK) = ZY(IIB:IIE,IJB:IJE,JK)
    !$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE)
  END DO
!
END IF 
!
!
!*       3.  FILL THE UPPER AND LOWER EXTERNAL VALUES
!            ----------------------------------------
!
!$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE)
PVARP(IIB:IIE,IJB:IJE,D%NKA)=PVARP(IIB:IIE,IJB:IJE,IKB)
PVARP(IIB:IIE,IJB:IJE,D%NKU)=PVARP(IIB:IIE,IJB:IJE,IKE)
!$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE)
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('TRIDIAG_WIND',1,ZHOOK_HANDLE)
END SUBROUTINE TRIDIAG_WIND
END MODULE MODE_TRIDIAG_WIND
