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
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PVARM   ! variable at t-1      at mass point
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PF      ! flux in dT/dt=-dF/dz at flux point
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PDFDDTDZ! dF/d(dT/dz)          at flux point
REAL,                   INTENT(IN) :: PTSTEP  ! Double time step
REAL,                   INTENT(IN) :: PIMPL   ! implicit weight
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PDZZ    ! Dz                   at flux point
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PRHODJ  ! (dry rho)*J          at mass point
!
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT):: PVARP   ! variable at t+1      at mass point
!
!
!*       0.2 declarations of local variables
!
REAL, DIMENSION(D%NIJT,D%NKT)  :: ZRHODJ_DFDDTDZ_O_DZ2
REAL, DIMENSION(D%NIJT,D%NKT)  :: ZMZM_RHODJ
REAL, DIMENSION(D%NIJT,D%NKT)  :: ZA, ZB, ZC
REAL, DIMENSION(D%NIJT,D%NKT)  :: ZY ,ZGAM 
                                         ! RHS of the equation, 3D work array
REAL, DIMENSION(D%NIJT)                :: ZBET
                                         ! 2D work array
INTEGER             :: JIJ,JK            ! loop counter
INTEGER             :: IKB,IKE ! inner limits
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
CALL MZM_PHY(D,PRHODJ,ZMZM_RHODJ)
!$mnh_expand_array(JIJ=1:D%NIJT,JK=1:D%NKT)
ZRHODJ_DFDDTDZ_O_DZ2(1:D%NIJT,1:D%NKT) = ZMZM_RHODJ(1:D%NIJT,1:D%NKT)*PDFDDTDZ(1:D%NIJT,1:D%NKT) &
                                                /PDZZ(1:D%NIJT,1:D%NKT)**2
!$mnh_end_expand_array(JIJ=1:D%NIJT,JK=1:D%NKT)
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
!$mnh_expand_array(JIJ=1:D%NIJT)
ZY(1:D%NIJT,IKB) = PRHODJ(1:D%NIJT,IKB)*PVARM(1:D%NIJT,IKB)/PTSTEP                  &
    - ZMZM_RHODJ(1:D%NIJT,IKB+D%NKL) * PF(1:D%NIJT,IKB+D%NKL)/PDZZ(1:D%NIJT,IKB+D%NKL)    &
    + ZMZM_RHODJ(1:D%NIJT,IKB  ) * PF(1:D%NIJT,IKB  )/PDZZ(1:D%NIJT,IKB  )          &
    + ZRHODJ_DFDDTDZ_O_DZ2(1:D%NIJT,IKB+D%NKL) * PIMPL * PVARM(1:D%NIJT,IKB+D%NKL) &
    - ZRHODJ_DFDDTDZ_O_DZ2(1:D%NIJT,IKB+D%NKL) * PIMPL * PVARM(1:D%NIJT,IKB  )
!$mnh_end_expand_array(JIJ=1:D%NIJT)
!
DO JK=IKTB+1,IKTE-1
  !$mnh_expand_array(JIJ=1:D%NIJT)
  ZY(1:D%NIJT,JK) = PRHODJ(1:D%NIJT,JK)*PVARM(1:D%NIJT,JK)/PTSTEP                 &
    - ZMZM_RHODJ(1:D%NIJT,JK+D%NKL) * PF(1:D%NIJT,JK+D%NKL)/PDZZ(1:D%NIJT,JK+D%NKL)     &
    + ZMZM_RHODJ(1:D%NIJT,JK  ) * PF(1:D%NIJT,JK  )/PDZZ(1:D%NIJT,JK  )           &
    + ZRHODJ_DFDDTDZ_O_DZ2(1:D%NIJT,JK+D%NKL) * PIMPL * PVARM(1:D%NIJT,JK+D%NKL) &
    - ZRHODJ_DFDDTDZ_O_DZ2(1:D%NIJT,JK+D%NKL) * PIMPL * PVARM(1:D%NIJT,JK  )   &
    - ZRHODJ_DFDDTDZ_O_DZ2(1:D%NIJT,JK    ) * PIMPL * PVARM(1:D%NIJT,JK  )   &
    + ZRHODJ_DFDDTDZ_O_DZ2(1:D%NIJT,JK    ) * PIMPL * PVARM(1:D%NIJT,JK-D%NKL)
  !$mnh_end_expand_array(JIJ=1:D%NIJT)
END DO
! 
!$mnh_expand_array(JIJ=1:D%NIJT)
ZY(1:D%NIJT,IKE) = PRHODJ(1:D%NIJT,IKE)*PVARM(1:D%NIJT,IKE)/PTSTEP               &
    - ZMZM_RHODJ(1:D%NIJT,IKE+D%NKL) * PF(1:D%NIJT,IKE+D%NKL)/PDZZ(1:D%NIJT,IKE+D%NKL) &
    + ZMZM_RHODJ(1:D%NIJT,IKE  ) * PF(1:D%NIJT,IKE  )/PDZZ(1:D%NIJT,IKE  )       &
    - ZRHODJ_DFDDTDZ_O_DZ2(1:D%NIJT,IKE ) * PIMPL * PVARM(1:D%NIJT,IKE  )   &
    + ZRHODJ_DFDDTDZ_O_DZ2(1:D%NIJT,IKE ) * PIMPL * PVARM(1:D%NIJT,IKE-D%NKL)
!$mnh_end_expand_array(JIJ=1:D%NIJT)
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
  !$mnh_expand_array(JIJ=1:D%NIJT)
  ZB(1:D%NIJT,IKB) =   PRHODJ(1:D%NIJT,IKB)/PTSTEP                   &
                - ZRHODJ_DFDDTDZ_O_DZ2(1:D%NIJT,IKB+D%NKL) * PIMPL
  ZC(1:D%NIJT,IKB) =   ZRHODJ_DFDDTDZ_O_DZ2(1:D%NIJT,IKB+D%NKL) * PIMPL
  !$mnh_end_expand_array(JIJ=1:D%NIJT)
!
  DO JK=IKTB+1,IKTE-1
    !$mnh_expand_array(JIJ=1:D%NIJT)
    ZA(1:D%NIJT,JK) =   ZRHODJ_DFDDTDZ_O_DZ2(1:D%NIJT,JK) * PIMPL
    ZB(1:D%NIJT,JK) =   PRHODJ(1:D%NIJT,JK)/PTSTEP                        &
                            - ZRHODJ_DFDDTDZ_O_DZ2(1:D%NIJT,JK+D%NKL) * PIMPL &
                            - ZRHODJ_DFDDTDZ_O_DZ2(1:D%NIJT,JK) * PIMPL
    ZC(1:D%NIJT,JK) =   ZRHODJ_DFDDTDZ_O_DZ2(1:D%NIJT,JK+D%NKL) * PIMPL
    !$mnh_end_expand_array(JIJ=1:D%NIJT)
  END DO
!
  !$mnh_expand_array(JIJ=1:D%NIJT)
  ZA(1:D%NIJT,IKE) =   ZRHODJ_DFDDTDZ_O_DZ2(1:D%NIJT,IKE  ) * PIMPL
  ZB(1:D%NIJT,IKE) =   PRHODJ(1:D%NIJT,IKE)/PTSTEP                   &
                - ZRHODJ_DFDDTDZ_O_DZ2(1:D%NIJT,IKE  ) * PIMPL
!
!*       3.2 going up
!            --------
!
  ZBET(1:D%NIJT) = ZB(1:D%NIJT,IKB)  ! bet = b(ikb)
  PVARP(1:D%NIJT,IKB) = ZY(1:D%NIJT,IKB) / ZBET(1:D%NIJT)
  !$mnh_end_expand_array(JIJ=1:D%NIJT)

  !
  DO JK = IKB+D%NKL,IKE-D%NKL,D%NKL
    !$mnh_expand_array(JIJ=1:D%NIJT)
    ZGAM(1:D%NIJT,JK) = ZC(1:D%NIJT,JK-D%NKL) / ZBET(1:D%NIJT)  
                                                    ! gam(k) = c(k-1) / bet
    ZBET(1:D%NIJT)    = ZB(1:D%NIJT,JK) - ZA(1:D%NIJT,JK) * ZGAM(1:D%NIJT,JK)
                                                    ! bet = b(k) - a(k)* gam(k)  
    PVARP(1:D%NIJT,JK)= ( ZY(1:D%NIJT,JK) - ZA(1:D%NIJT,JK) * PVARP(1:D%NIJT,JK-D%NKL) ) &
                               / ZBET(1:D%NIJT)
                                        ! res(k) = (y(k) -a(k)*res(k-1))/ bet 
    !$mnh_end_expand_array(JIJ=1:D%NIJT)
  END DO 
  ! special treatment for the last level
  !$mnh_expand_array(JIJ=1:D%NIJT)
  ZGAM(1:D%NIJT,IKE) = ZC(1:D%NIJT,IKE-D%NKL) / ZBET(1:D%NIJT) 
                                                    ! gam(k) = c(k-1) / bet
  ZBET(1:D%NIJT)     = ZB(1:D%NIJT,IKE) - ZA(1:D%NIJT,IKE) * ZGAM(1:D%NIJT,IKE)
                                                    ! bet = b(k) - a(k)* gam(k)  
  PVARP(1:D%NIJT,IKE)= ( ZY(1:D%NIJT,IKE) - ZA(1:D%NIJT,IKE) * PVARP(1:D%NIJT,IKE-D%NKL) ) &
                              / ZBET(1:D%NIJT)
                                       ! res(k) = (y(k) -a(k)*res(k-1))/ bet 
  !$mnh_end_expand_array(JIJ=1:D%NIJT)
!
!*       3.3 going down
!            ----------
!
  DO JK = IKE-D%NKL,IKB,-1*D%NKL
    !$mnh_expand_array(JIJ=1:D%NIJT)
    PVARP(1:D%NIJT,JK) = PVARP(1:D%NIJT,JK) - ZGAM(1:D%NIJT,JK+D%NKL) * PVARP(1:D%NIJT,JK+D%NKL)
    !$mnh_end_expand_array(JIJ=1:D%NIJT)
  END DO
!
ELSE
! 
  DO JK=IKTB,IKTE
    !$mnh_expand_array(JIJ=1:D%NIJT)
    PVARP(1:D%NIJT,JK) = ZY(1:D%NIJT,JK) * PTSTEP / PRHODJ(1:D%NIJT,JK)
    !$mnh_end_expand_array(JIJ=1:D%NIJT)
  END DO
!
END IF 
!
!
!*       4.  FILL THE UPPER AND LOWER EXTERNAL VALUES
!            ----------------------------------------
!
!$mnh_expand_array(JIJ=1:D%NIJT)
PVARP(1:D%NIJT,D%NKA)=PVARP(1:D%NIJT,IKB)
PVARP(1:D%NIJT,D%NKU)=PVARP(1:D%NIJT,IKE)
!$mnh_end_expand_array(JIJ=1:D%NIJT)
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('TRIDIAG_THERMO',1,ZHOOK_HANDLE)
END SUBROUTINE TRIDIAG_THERMO
END MODULE MODE_TRIDIAG_THERMO
