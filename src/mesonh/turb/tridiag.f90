!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 turb 2006/06/06 09:55:03
!-----------------------------------------------------------------
!     ###################
      MODULE MODI_TRIDIAG
!     ###################
INTERFACE
!
       SUBROUTINE TRIDIAG(KKA,KKU,KKL,PVARM,PA,PTSTEP,PEXPL,PIMPL, &
                                  PRHODJ,PSOURCE,PVARP )
!
INTEGER,                INTENT(IN)   :: KKA           !near ground array index  
INTEGER,                INTENT(IN)   :: KKU           !uppest atmosphere array index
INTEGER,                INTENT(IN)   :: KKL           !vert. levels type 1=MNH -1=AR
REAL, DIMENSION(:,:,:),    INTENT(IN)  :: PVARM       ! variable at t-1  
REAL, DIMENSION(:,:,:),    INTENT(IN)  :: PA          ! upper diag. elements
REAL,                      INTENT(IN)  :: PTSTEP      ! Double time step
REAL,                      INTENT(IN)  :: PEXPL,PIMPL ! weights of the temporal scheme
REAL, DIMENSION(:,:,:),    INTENT(IN)  :: PRHODJ      ! (dry rho)*J
REAL, DIMENSION(:,:,:),    INTENT(IN)  :: PSOURCE     ! source term of PVAR    
!
REAL, DIMENSION(:,:,:),    INTENT(OUT) :: PVARP       ! variable at t+1        
!
END SUBROUTINE TRIDIAG
!
END INTERFACE
!
END MODULE MODI_TRIDIAG 
!
!
!
!      #################################################
       SUBROUTINE TRIDIAG(KKA,KKU,KKL,PVARM,PA,PTSTEP,PEXPL,PIMPL, &
                                  PRHODJ,PSOURCE,PVARP )
!      #################################################
!
!
!!****   *TRIDIAG* - routine to solve a time implicit scheme
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
!!      (Seity)         February 2012 add possibility to run with reversed 
!!                            vertical levels
!! ---------------------------------------------------------------------
!
!*       0. DECLARATIONS
!
USE MODD_PARAMETERS
!
IMPLICIT NONE
!
!
!*       0.1 declarations of arguments
!
INTEGER,                   INTENT(IN)  :: KKA        !near ground array index  
INTEGER,                   INTENT(IN)  :: KKU        !uppest atmosphere array index
INTEGER,                   INTENT(IN)  :: KKL        !vert. levels type 1=MNH -1=ARO
REAL, DIMENSION(:,:,:),    INTENT(IN)  :: PVARM       ! variable at t-1  
REAL, DIMENSION(:,:,:),    INTENT(IN)  :: PA          ! upper diag. elements
REAL,                      INTENT(IN)  :: PTSTEP      ! Double time step
REAL,                      INTENT(IN)  :: PEXPL,PIMPL ! weights of the temporal scheme
REAL, DIMENSION(:,:,:),    INTENT(IN)  :: PRHODJ      ! (dry rho)*J
REAL, DIMENSION(:,:,:),    INTENT(IN)  :: PSOURCE     ! source term of PVAR    
!
REAL, DIMENSION(:,:,:),    INTENT(OUT) :: PVARP       ! variable at t+1        
!
!*       0.2 declarations of local variables
!
REAL, DIMENSION(SIZE(PVARM,1),SIZE(PVARM,2),SIZE(PVARM,3))  :: ZY ,ZGAM 
                                         ! RHS of the equation, 3D work array
REAL, DIMENSION(SIZE(PVARM,1),SIZE(PVARM,2))                :: ZBET
                                         ! 2D work array
INTEGER                              :: JK            ! loop counter
INTEGER                              :: IKB,IKE       ! inner vertical limits
INTEGER                              :: IKT           ! array size in k direction
INTEGER                              :: IKTB,IKTE     ! start, end of k loops in physical domain 

!
! ---------------------------------------------------------------------------
!                                              
!*      1.  COMPUTE THE RIGHT HAND SIDE
!           ---------------------------
!
IKTB=1+JPVEXT_TURB
IKT=SIZE(PVARM,3)
IKTE=IKT-JPVEXT_TURB 
IKB=KKA+JPVEXT_TURB*KKL
IKE=KKU-JPVEXT_TURB*KKL
!
!
ZY(:,:,IKB) = PVARM(:,:,IKB)  + PTSTEP*PSOURCE(:,:,IKB) -   &
  PEXPL / PRHODJ(:,:,IKB) * PA(:,:,IKB+KKL) * (PVARM(:,:,IKB+KKL) - PVARM(:,:,IKB))
!
DO JK=IKTB+1,IKTE-1
  ZY(:,:,JK)= PVARM(:,:,JK)  + PTSTEP*PSOURCE(:,:,JK) -               &
      PEXPL / PRHODJ(:,:,JK) *                                           &
                             ( PVARM(:,:,JK-KKL)*PA(:,:,JK)                &
                              -PVARM(:,:,JK)*(PA(:,:,JK)+PA(:,:,JK+KKL))   &
                              +PVARM(:,:,JK+KKL)*PA(:,:,JK+KKL)              &
                             ) 
END DO
! 
ZY(:,:,IKE)= PVARM(:,:,IKE) + PTSTEP*PSOURCE(:,:,IKE) +               &
  PEXPL / PRHODJ(:,:,IKE) * PA(:,:,IKE) * (PVARM(:,:,IKE)-PVARM(:,:,IKE-KKL))
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
  ZBET(:,:) = 1. - PIMPL * PA(:,:,IKB+KKL) / PRHODJ(:,:,IKB)  ! bet = b(ikb)
  PVARP(:,:,IKB) = ZY(:,:,IKB) / ZBET(:,:)                
  !
  DO JK = IKB+KKL,IKE-KKL,KKL
      ZGAM(:,:,JK) = PIMPL * PA(:,:,JK) / PRHODJ(:,:,JK-KKL) / ZBET(:,:)  
                                                    ! gam(k) = c(k-1) / bet
    ZBET(:,:)    = 1. - PIMPL * (  PA(:,:,JK) * (1. + ZGAM(:,:,JK))  &
                                 + PA(:,:,JK+KKL)                      &
                                ) / PRHODJ(:,:,JK)  
                                                    ! bet = b(k) - a(k)* gam(k)  
    PVARP(:,:,JK)= ( ZY(:,:,JK) - PIMPL * PA(:,:,JK) / PRHODJ(:,:,JK) &
                    * PVARP(:,:,JK-KKL)                                 &
                   ) / ZBET(:,:)
                                        ! res(k) = (y(k) -a(k)*res(k-1))/ bet 
  END DO 
  ! special treatment for the last level
  ZGAM(:,:,IKE) = PIMPL * PA(:,:,IKE) / PRHODJ(:,:,IKE-KKL) / ZBET(:,:) 
                                                    ! gam(k) = c(k-1) / bet
  ZBET(:,:)    = 1. - PIMPL * (  PA(:,:,IKE) * (1. + ZGAM(:,:,IKE))  &
                              ) / PRHODJ(:,:,IKE)  
                                                    ! bet = b(k) - a(k)* gam(k)  
  PVARP(:,:,IKE)= ( ZY(:,:,IKE) - PIMPL * PA(:,:,IKE) / PRHODJ(:,:,IKE) &
                                * PVARP(:,:,IKE-KKL)                      &
                 ) / ZBET(:,:)
                                       ! res(k) = (y(k) -a(k)*res(k-1))/ bet 
  !
  !  going down
  !
  DO JK = IKE-KKL,IKB,-1*KKL
    PVARP(:,:,JK) = PVARP(:,:,JK) - ZGAM(:,:,JK+KKL) * PVARP(:,:,JK+KKL) 
  END DO
!
ELSE
! 
  PVARP(:,:,IKTB:IKTE) = ZY(:,:,IKTB:IKTE)
!
END IF 
!
!
!*       3.  FILL THE UPPER AND LOWER EXTERNAL VALUES
!            ----------------------------------------
!
PVARP(:,:,KKA)=PVARP(:,:,IKB)
PVARP(:,:,KKU)=PVARP(:,:,IKE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE TRIDIAG
