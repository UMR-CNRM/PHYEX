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
USE SHUMAN_PHY, ONLY: MZM_PHY
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
INTEGER             :: IKB,IKE,IIB,IIE,IJB,IJE ! inner limits
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
IIE=D%NIEC
IIB=D%NIBC
IJE=D%NJEC
IJB=D%NJBC

!
CALL MZM_PHY(D,PRHODJ,ZMZM_RHODJ)
!$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)
ZRHODJ_DFDDTDZ_O_DZ2(IIB:IIE,IJB:IJE,1:D%NKT) = ZMZM_RHODJ(IIB:IIE,IJB:IJE,1:D%NKT)*PDFDDTDZ(IIB:IIE,IJB:IJE,1:D%NKT) &
                                                /PDZZ(IIB:IIE,IJB:IJE,1:D%NKT)**2
!$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)
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
!$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE)
ZY(IIB:IIE,IJB:IJE,IKB) = PRHODJ(IIB:IIE,IJB:IJE,IKB)*PVARM(IIB:IIE,IJB:IJE,IKB)/PTSTEP                  &
    - ZMZM_RHODJ(IIB:IIE,IJB:IJE,IKB+D%NKL) * PF(IIB:IIE,IJB:IJE,IKB+D%NKL)/PDZZ(IIB:IIE,IJB:IJE,IKB+D%NKL)    &
    + ZMZM_RHODJ(IIB:IIE,IJB:IJE,IKB  ) * PF(IIB:IIE,IJB:IJE,IKB  )/PDZZ(IIB:IIE,IJB:IJE,IKB  )          &
    + ZRHODJ_DFDDTDZ_O_DZ2(IIB:IIE,IJB:IJE,IKB+D%NKL) * PIMPL * PVARM(IIB:IIE,IJB:IJE,IKB+D%NKL) &
    - ZRHODJ_DFDDTDZ_O_DZ2(IIB:IIE,IJB:IJE,IKB+D%NKL) * PIMPL * PVARM(IIB:IIE,IJB:IJE,IKB  )
!$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE)
!
!$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE)
DO JK=IKTB+1,IKTE-1
  ZY(IIB:IIE,IJB:IJE,JK) = PRHODJ(IIB:IIE,IJB:IJE,JK)*PVARM(IIB:IIE,IJB:IJE,JK)/PTSTEP                 &
    - ZMZM_RHODJ(IIB:IIE,IJB:IJE,JK+D%NKL) * PF(IIB:IIE,IJB:IJE,JK+D%NKL)/PDZZ(IIB:IIE,IJB:IJE,JK+D%NKL)     &
    + ZMZM_RHODJ(IIB:IIE,IJB:IJE,JK  ) * PF(IIB:IIE,IJB:IJE,JK  )/PDZZ(IIB:IIE,IJB:IJE,JK  )           &
    + ZRHODJ_DFDDTDZ_O_DZ2(IIB:IIE,IJB:IJE,JK+D%NKL) * PIMPL * PVARM(IIB:IIE,IJB:IJE,JK+D%NKL) &
    - ZRHODJ_DFDDTDZ_O_DZ2(IIB:IIE,IJB:IJE,JK+D%NKL) * PIMPL * PVARM(IIB:IIE,IJB:IJE,JK  )   &
    - ZRHODJ_DFDDTDZ_O_DZ2(IIB:IIE,IJB:IJE,JK    ) * PIMPL * PVARM(IIB:IIE,IJB:IJE,JK  )   &
    + ZRHODJ_DFDDTDZ_O_DZ2(IIB:IIE,IJB:IJE,JK    ) * PIMPL * PVARM(IIB:IIE,IJB:IJE,JK-D%NKL)
END DO
!$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE)
! 
!$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE)
ZY(IIB:IIE,IJB:IJE,IKE) = PRHODJ(IIB:IIE,IJB:IJE,IKE)*PVARM(IIB:IIE,IJB:IJE,IKE)/PTSTEP               &
    - ZMZM_RHODJ(IIB:IIE,IJB:IJE,IKE+D%NKL) * PF(IIB:IIE,IJB:IJE,IKE+D%NKL)/PDZZ(IIB:IIE,IJB:IJE,IKE+D%NKL) &
    + ZMZM_RHODJ(IIB:IIE,IJB:IJE,IKE  ) * PF(IIB:IIE,IJB:IJE,IKE  )/PDZZ(IIB:IIE,IJB:IJE,IKE  )       &
    - ZRHODJ_DFDDTDZ_O_DZ2(IIB:IIE,IJB:IJE,IKE ) * PIMPL * PVARM(IIB:IIE,IJB:IJE,IKE  )   &
    + ZRHODJ_DFDDTDZ_O_DZ2(IIB:IIE,IJB:IJE,IKE ) * PIMPL * PVARM(IIB:IIE,IJB:IJE,IKE-D%NKL)
!$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE)
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
  !$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE)
  ZB(IIB:IIE,IJB:IJE,IKB) =   PRHODJ(IIB:IIE,IJB:IJE,IKB)/PTSTEP                   &
                - ZRHODJ_DFDDTDZ_O_DZ2(IIB:IIE,IJB:IJE,IKB+D%NKL) * PIMPL
  ZC(IIB:IIE,IJB:IJE,IKB) =   ZRHODJ_DFDDTDZ_O_DZ2(IIB:IIE,IJB:IJE,IKB+D%NKL) * PIMPL
  !$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE)
!
  DO JK=IKTB+1,IKTE-1
    !$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE)
    ZA(IIB:IIE,IJB:IJE,JK) =   ZRHODJ_DFDDTDZ_O_DZ2(IIB:IIE,IJB:IJE,JK) * PIMPL
    ZB(IIB:IIE,IJB:IJE,JK) =   PRHODJ(IIB:IIE,IJB:IJE,JK)/PTSTEP                        &
                            - ZRHODJ_DFDDTDZ_O_DZ2(IIB:IIE,IJB:IJE,JK+D%NKL) * PIMPL &
                            - ZRHODJ_DFDDTDZ_O_DZ2(IIB:IIE,IJB:IJE,JK) * PIMPL
    ZC(IIB:IIE,IJB:IJE,JK) =   ZRHODJ_DFDDTDZ_O_DZ2(IIB:IIE,IJB:IJE,JK+D%NKL) * PIMPL
    !$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE)
  END DO
!
  !$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE)
  ZA(IIB:IIE,IJB:IJE,IKE) =   ZRHODJ_DFDDTDZ_O_DZ2(IIB:IIE,IJB:IJE,IKE  ) * PIMPL
  ZB(IIB:IIE,IJB:IJE,IKE) =   PRHODJ(IIB:IIE,IJB:IJE,IKE)/PTSTEP                   &
                - ZRHODJ_DFDDTDZ_O_DZ2(IIB:IIE,IJB:IJE,IKE  ) * PIMPL
!
!*       3.2 going up
!            --------
!
  ZBET(IIB:IIE,IJB:IJE) = ZB(IIB:IIE,IJB:IJE,IKB)  ! bet = b(ikb)
  PVARP(IIB:IIE,IJB:IJE,IKB) = ZY(IIB:IIE,IJB:IJE,IKB) / ZBET(IIB:IIE,IJB:IJE)
  !$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE)

  !
  DO JK = IKB+D%NKL,IKE-D%NKL,D%NKL
    !$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE)
    ZGAM(IIB:IIE,IJB:IJE,JK) = ZC(IIB:IIE,IJB:IJE,JK-D%NKL) / ZBET(IIB:IIE,IJB:IJE)  
                                                    ! gam(k) = c(k-1) / bet
    ZBET(IIB:IIE,IJB:IJE)    = ZB(IIB:IIE,IJB:IJE,JK) - ZA(IIB:IIE,IJB:IJE,JK) * ZGAM(IIB:IIE,IJB:IJE,JK)
                                                    ! bet = b(k) - a(k)* gam(k)  
    PVARP(IIB:IIE,IJB:IJE,JK)= ( ZY(IIB:IIE,IJB:IJE,JK) - ZA(IIB:IIE,IJB:IJE,JK) * PVARP(IIB:IIE,IJB:IJE,JK-D%NKL) ) &
                               / ZBET(IIB:IIE,IJB:IJE)
                                        ! res(k) = (y(k) -a(k)*res(k-1))/ bet 
    !$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE)
  END DO 
  ! special treatment for the last level
  !$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE)
  ZGAM(IIB:IIE,IJB:IJE,IKE) = ZC(IIB:IIE,IJB:IJE,IKE-D%NKL) / ZBET(IIB:IIE,IJB:IJE) 
                                                    ! gam(k) = c(k-1) / bet
  ZBET(IIB:IIE,IJB:IJE)     = ZB(IIB:IIE,IJB:IJE,IKE) - ZA(IIB:IIE,IJB:IJE,IKE) * ZGAM(IIB:IIE,IJB:IJE,IKE)
                                                    ! bet = b(k) - a(k)* gam(k)  
  PVARP(IIB:IIE,IJB:IJE,IKE)= ( ZY(IIB:IIE,IJB:IJE,IKE) - ZA(IIB:IIE,IJB:IJE,IKE) * PVARP(IIB:IIE,IJB:IJE,IKE-D%NKL) ) &
                              / ZBET(IIB:IIE,IJB:IJE)
                                       ! res(k) = (y(k) -a(k)*res(k-1))/ bet 
  !$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE)
!
!*       3.3 going down
!            ----------
!
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
    PVARP(IIB:IIE,IJB:IJE,JK) = ZY(IIB:IIE,IJB:IJE,JK) * PTSTEP / PRHODJ(IIB:IIE,IJB:IJE,JK)
    !$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE)
  END DO
!
END IF 
!
!
!*       4.  FILL THE UPPER AND LOWER EXTERNAL VALUES
!            ----------------------------------------
!
!$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE)
PVARP(IIB:IIE,IJB:IJE,D%NKA)=PVARP(IIB:IIE,IJB:IJE,IKB)
PVARP(IIB:IIE,IJB:IJE,D%NKU)=PVARP(IIB:IIE,IJB:IJE,IKE)
!$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE)
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('TRIDIAG_THERMO',1,ZHOOK_HANDLE)
END SUBROUTINE TRIDIAG_THERMO
END MODULE MODE_TRIDIAG_THERMO
