!MNH_LIC Copyright 2003-2020 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
MODULE MODE_TRIDIAG_THERMO
IMPLICIT NONE
CONTAINS       
SUBROUTINE TRIDIAG_THERMO(D,PVARM,PF,PDFDDTDZ,PTSTEP,PIMPL,  &
                                 PDZZ,PRHODJ,PVARP             )
!      #################################################
!
!
!!****   *TRIDIAG_THERMO* - routine to solve a time implicit scheme
!!
!!
!!     PURPOSE
!!     -------
!        The purpose of this routine is to give a field PVARP at t+1, by 
!      solving an implicit TRIDIAGonal system obtained by the 
!      discretization of the vertical turbulent diffusion. It should be noted 
!      that the degree of implicitness can be varied (PIMPL parameter) and that
!      the function of F(dT/dz) must have been linearized.
!      PVARP is localized at a mass point.
!
!!**   METHOD
!!     ------
!!
!!        [T(+) - T(-)]/2Dt = -d{ F + dF/d(dT/dz) * [impl*dT/dz(+) + expl* dT/dz(-)] }/dz
!!
!!     It is discretized as follows:
!!
!!    PRHODJ(k)*PVARP(k)/PTSTEP
!!              = 
!!    PRHODJ(k)*PVARM(k)/PTSTEP 
!!  - (PRHODJ(k+1)+PRHODJ(k)  )/2. * PF(k+1)/PDZZ(k+1)
!!  + (PRHODJ(k)  +PRHODJ(k-1))/2. * PF(k)  /PDZZ(k)
!!  + (PRHODJ(k+1)+PRHODJ(k)  )/2. * ZEXPL* PDFDDTDZ(k+1) * PVARM(k+1)/PDZZ(k+1)**2
!!  - (PRHODJ(k+1)+PRHODJ(k)  )/2. * PIMPL* PDFDDTDZ(k+1) * PVARP(k+1)/PDZZ(k+1)**2
!!  - (PRHODJ(k+1)+PRHODJ(k)  )/2. * ZEXPL* PDFDDTDZ(k+1) * PVARM(k)  /PDZZ(k+1)**2
!!  + (PRHODJ(k+1)+PRHODJ(k)  )/2. * PIMPL* PDFDDTDZ(k+1) * PVARP(k)  /PDZZ(k+1)**2
!!  - (PRHODJ(k)  +PRHODJ(k-1))/2. * ZEXPL* PDFDDTDZ(k)   * PVARM(k)  /PDZZ(k)**2
!!  + (PRHODJ(k)  +PRHODJ(k-1))/2. * PIMPL* PDFDDTDZ(k)   * PVARP(k)  /PDZZ(k)**2
!!  + (PRHODJ(k)  +PRHODJ(k-1))/2. * ZEXPL* PDFDDTDZ(k)   * PVARM(k-1)/PDZZ(k)**2
!!  - (PRHODJ(k)  +PRHODJ(k-1))/2. * PIMPL* PDFDDTDZ(k)   * PVARP(k-1)/PDZZ(k)**2
!!
!!
!!    The system to solve is:
!!
!!      A*PVARP(k-1) + B*PVARP(k) + C*PVARP(k+1) = Y(k)
!!
!!
!!    The RHS of the linear system in PVARP writes:
!!
!! y(k)    = PRHODJ(k)*PVARM(k)/PTSTEP
!!  - (PRHODJ(k+1)+PRHODJ(k)  )/2. * PF(k+1)/PDZZ(k+1)
!!  + (PRHODJ(k)  +PRHODJ(k-1))/2. * PF(k)  /PDZZ(k)
!!  + (PRHODJ(k+1)+PRHODJ(k)  )/2. * ZEXPL* PDFDDTDZ(k+1) * PVARM(k+1)/PDZZ(k+1)**2
!!  - (PRHODJ(k+1)+PRHODJ(k)  )/2. * ZEXPL* PDFDDTDZ(k+1) * PVARM(k)  /PDZZ(k+1)**2
!!  - (PRHODJ(k)  +PRHODJ(k-1))/2. * ZEXPL* PDFDDTDZ(k)   * PVARM(k)  /PDZZ(k)**2
!!  + (PRHODJ(k)  +PRHODJ(k-1))/2. * ZEXPL* PDFDDTDZ(k)   * PVARM(k-1)/PDZZ(k)**2
!!
!!                      
!!        Then, the classical TRIDIAGonal algorithm is used to invert the 
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
!! a(k) = + (PRHODJ(k)  +PRHODJ(k-1))/2. * PIMPL* PDFDDTDZ(k)  /PDZZ(k)**2
!! b(k) =    PRHODJ(k) / PTSTEP
!!        - (PRHODJ(k+1)+PRHODJ(k)  )/2. * PIMPL* PDFDDTDZ(k+1)/PDZZ(k+1)**2
!!        - (PRHODJ(k)  +PRHODJ(k-1))/2. * PIMPL* PDFDDTDZ(k)  /PDZZ(k)**2
!! c(k) = + (PRHODJ(k+1)+PRHODJ(k)  )/2. * PIMPL* PDFDDTDZ(k+1)/PDZZ(k+1)**2
!!
!!          for all k /= ikb or ike
!!
!!
!! b(ikb) =  PRHODJ(ikb) / PTSTEP
!!          -(PRHODJ(ikb+1)+PRHODJ(ikb))/2.*PIMPL*PDFDDTDZ(ikb+1)/PDZZ(ikb+1)**2
!! c(ikb) = +(PRHODJ(ikb+1)+PRHODJ(ikb))/2.*PIMPL*PDFDDTDZ(ikb+1)/PDZZ(ikb+1)**2
!!
!! b(ike) =  PRHODJ(ike) / PTSTEP
!!          -(PRHODJ(ike)+PRHODJ(ike-1))/2.*PIMPL*PDFDDTDZ(ike)/PDZZ(ike)**2
!! a(ike) = +(PRHODJ(ike)+PRHODJ(ike-1))/2.*PIMPL*PDFDDTDZ(ike)/PDZZ(ike)**2
!!
!!
!!     EXTERNAL
!!     --------
!!
!!       NONE
!!
!!     IMPLICIT ARGUMENTS
!!     ------------------
!!
!!     REFERENCE
!!     ---------
!!       Press et al: Numerical recipes (1986) Cambridge Univ. Press
!!
!!     AUTHOR
!!     ------
!!       V. Masson         * Meteo-France *   
!! 
!!     MODIFICATIONS
!!     -------------
!!       Original        04/2003 (from tridiag.f90)
!! ---------------------------------------------------------------------
!
!*       0. DECLARATIONS
!
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
USE MODD_DIMPHYEX, ONLY : DIMPHYEX_t
USE MODD_PARAMETERS, ONLY : JPVEXT_TURB
!
USE MODI_SHUMAN, ONLY : MZM
!
IMPLICIT NONE
!
!
!*       0.1 declarations of arguments
!
TYPE(DIMPHYEX_t),     INTENT(IN)   :: D
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PVARM   ! variable at t-1      at mass point
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PF      ! flux in dT/dt=-dF/dz at flux point
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PDFDDTDZ! dF/d(dT/dz)          at flux point
REAL,                   INTENT(IN) :: PTSTEP  ! Double time step
REAL,                   INTENT(IN) :: PIMPL   ! implicit weight
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PDZZ    ! Dz                   at flux point
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PRHODJ  ! (dry rho)*J          at mass point
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT):: PVARP   ! variable at t+1      at mass point
!
!
!*       0.2 declarations of local variables
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT)  :: ZRHODJ_DFDDTDZ_O_DZ2
REAL, DIMENSION(D%NIT,D%NJT,D%NKT)  :: ZMZM_RHODJ
REAL, DIMENSION(D%NIT,D%NJT,D%NKT)  :: ZA, ZB, ZC
REAL, DIMENSION(D%NIT,D%NJT,D%NKT)  :: ZY ,ZGAM 
                                         ! RHS of the equation, 3D work array
REAL, DIMENSION(D%NIT,D%NJT)                :: ZBET
                                         ! 2D work array
INTEGER             :: JI,JJ,JK            ! loop counter
INTEGER             :: IKB,IKE       ! inner vertical limits
INTEGER             :: IKT          ! array size in k direction
INTEGER             :: IKTB,IKTE    ! start, end of k loops in physical domain 
!
! ---------------------------------------------------------------------------
!                                              
!*      1.  Preliminaries
!           -------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('TRIDIAG_THERMO',0,ZHOOK_HANDLE)
IKT=D%NKT  
IKTB=D%NKTB          
IKTE=D%NKTE
IKB=D%NKB
IKE=D%NKE

!
ZMZM_RHODJ  = MZM(PRHODJ,D%NKA,D%NKU,D%NKL)
!$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
ZRHODJ_DFDDTDZ_O_DZ2(:,:,:) = ZMZM_RHODJ(:,:,:)*PDFDDTDZ(:,:,:)/PDZZ(:,:,:)**2
!$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
!
ZA=0.
ZB=0.
ZC=0.
ZY=0.
!
!
!*      2.  COMPUTE THE RIGHT HAND SIDE
!           ---------------------------
!
!$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT)
ZY(:,:,IKB) = PRHODJ(:,:,IKB)*PVARM(:,:,IKB)/PTSTEP                  &
    - ZMZM_RHODJ(:,:,IKB+D%NKL) * PF(:,:,IKB+D%NKL)/PDZZ(:,:,IKB+D%NKL)    &
    + ZMZM_RHODJ(:,:,IKB  ) * PF(:,:,IKB  )/PDZZ(:,:,IKB  )          &
    + ZRHODJ_DFDDTDZ_O_DZ2(:,:,IKB+D%NKL) * PIMPL * PVARM(:,:,IKB+D%NKL) &
    - ZRHODJ_DFDDTDZ_O_DZ2(:,:,IKB+D%NKL) * PIMPL * PVARM(:,:,IKB  )
!$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT)
!
!$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT)
DO JK=IKTB+1,IKTE-1
  ZY(:,:,JK) = PRHODJ(:,:,JK)*PVARM(:,:,JK)/PTSTEP                 &
    - ZMZM_RHODJ(:,:,JK+D%NKL) * PF(:,:,JK+D%NKL)/PDZZ(:,:,JK+D%NKL)     &
    + ZMZM_RHODJ(:,:,JK  ) * PF(:,:,JK  )/PDZZ(:,:,JK  )           &
    + ZRHODJ_DFDDTDZ_O_DZ2(:,:,JK+D%NKL) * PIMPL * PVARM(:,:,JK+D%NKL) &
    - ZRHODJ_DFDDTDZ_O_DZ2(:,:,JK+D%NKL) * PIMPL * PVARM(:,:,JK  )   &
    - ZRHODJ_DFDDTDZ_O_DZ2(:,:,JK    ) * PIMPL * PVARM(:,:,JK  )   &
    + ZRHODJ_DFDDTDZ_O_DZ2(:,:,JK    ) * PIMPL * PVARM(:,:,JK-D%NKL)
END DO
!$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT)
! 
!$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT)
ZY(:,:,IKE) = PRHODJ(:,:,IKE)*PVARM(:,:,IKE)/PTSTEP               &
    - ZMZM_RHODJ(:,:,IKE+D%NKL) * PF(:,:,IKE+D%NKL)/PDZZ(:,:,IKE+D%NKL) &
    + ZMZM_RHODJ(:,:,IKE  ) * PF(:,:,IKE  )/PDZZ(:,:,IKE  )       &
    - ZRHODJ_DFDDTDZ_O_DZ2(:,:,IKE ) * PIMPL * PVARM(:,:,IKE  )   &
    + ZRHODJ_DFDDTDZ_O_DZ2(:,:,IKE ) * PIMPL * PVARM(:,:,IKE-D%NKL)
!$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT)
!
!
!*       3.  INVERSION OF THE TRIDIAGONAL SYSTEM
!            -----------------------------------
!
IF ( PIMPL > 1.E-10 ) THEN
!
!*       3.1 arrays A, B, C
!            --------------
!
  !$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT)
  ZB(:,:,IKB) =   PRHODJ(:,:,IKB)/PTSTEP                   &
                - ZRHODJ_DFDDTDZ_O_DZ2(:,:,IKB+D%NKL) * PIMPL
  ZC(:,:,IKB) =   ZRHODJ_DFDDTDZ_O_DZ2(:,:,IKB+D%NKL) * PIMPL
!
  DO JK=IKTB+1,IKTE-1
    ZA(:,:,JK) =   ZRHODJ_DFDDTDZ_O_DZ2(:,:,JK) * PIMPL
    ZB(:,:,JK) =   PRHODJ(:,:,JK)/PTSTEP                        &
                            - ZRHODJ_DFDDTDZ_O_DZ2(:,:,JK+D%NKL) * PIMPL &
                            - ZRHODJ_DFDDTDZ_O_DZ2(:,:,JK) * PIMPL
    ZC(:,:,JK) =   ZRHODJ_DFDDTDZ_O_DZ2(:,:,JK+D%NKL) * PIMPL
  END DO
!
  ZA(:,:,IKE) =   ZRHODJ_DFDDTDZ_O_DZ2(:,:,IKE  ) * PIMPL
  ZB(:,:,IKE) =   PRHODJ(:,:,IKE)/PTSTEP                   &
                - ZRHODJ_DFDDTDZ_O_DZ2(:,:,IKE  ) * PIMPL
!
!*       3.2 going up
!            --------
!
  ZBET(:,:) = ZB(:,:,IKB)  ! bet = b(ikb)
  PVARP(:,:,IKB) = ZY(:,:,IKB) / ZBET(:,:)

  !
  DO JK = IKB+D%NKL,IKE-D%NKL,D%NKL
    ZGAM(:,:,JK) = ZC(:,:,JK-D%NKL) / ZBET(:,:)  
                                                    ! gam(k) = c(k-1) / bet
    ZBET(:,:)    = ZB(:,:,JK) - ZA(:,:,JK) * ZGAM(:,:,JK)
                                                    ! bet = b(k) - a(k)* gam(k)  
    PVARP(:,:,JK)= ( ZY(:,:,JK) - ZA(:,:,JK) * PVARP(:,:,JK-D%NKL) ) / ZBET(:,:)
                                        ! res(k) = (y(k) -a(k)*res(k-1))/ bet 
  END DO 
  ! special treatment for the last level
  ZGAM(:,:,IKE) = ZC(:,:,IKE-D%NKL) / ZBET(:,:) 
                                                    ! gam(k) = c(k-1) / bet
  ZBET(:,:)     = ZB(:,:,IKE) - ZA(:,:,IKE) * ZGAM(:,:,IKE)
                                                    ! bet = b(k) - a(k)* gam(k)  
  PVARP(:,:,IKE)= ( ZY(:,:,IKE) - ZA(:,:,IKE) * PVARP(:,:,IKE-D%NKL) ) / ZBET(:,:)
                                       ! res(k) = (y(k) -a(k)*res(k-1))/ bet 
!
!*       3.3 going down
!            ----------
!
  DO JK = IKE-D%NKL,IKB,-1*D%NKL
    PVARP(:,:,JK) = PVARP(:,:,JK) - ZGAM(:,:,JK+D%NKL) * PVARP(:,:,JK+D%NKL)
  END DO
!
  !$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT)
ELSE
! 
  !$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT)
  DO JK=IKTB,IKTE
    PVARP(:,:,JK) = ZY(:,:,JK) * PTSTEP / PRHODJ(:,:,JK)
  END DO
  !$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT)
!
END IF 
!
!
!*       4.  FILL THE UPPER AND LOWER EXTERNAL VALUES
!            ----------------------------------------
!
!$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT)
PVARP(:,:,D%NKA)=PVARP(:,:,IKB)
PVARP(:,:,D%NKU)=PVARP(:,:,IKE)
!$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT)
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('TRIDIAG_THERMO',1,ZHOOK_HANDLE)
END SUBROUTINE TRIDIAG_THERMO
END MODULE MODE_TRIDIAG_THERMO
