!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
MODULE MODE_TRIDIAG_MASSFLUX
IMPLICIT NONE
CONTAINS
SUBROUTINE TRIDIAG_MASSFLUX(D,PVARM,PF,PDFDT,PTSTEP,PIMPL,  &
                                 PDZZ,PRHODJ,PVARP             )

       USE PARKIND1, ONLY : JPRB
       USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!      #################################################
!
!
!!****   *TRIDIAG_MASSFLUX* - routine to solve a time implicit scheme
!!
!!
!!     PURPOSE
!!     -------
!        The purpose of this routine is to give a field PVARP at t+1, by 
!      solving an implicit TRIDIAGonal system obtained by the 
!      discretization of the vertical turbulent diffusion. It should be noted 
!      that the degree of implicitness can be varied (PIMPL parameter) and that
!      the function of F(T) must have been linearized.
!      PVARP is localized at a mass point.
!
!!**   METHOD
!!     ------
!!
!!        [T(+) - T(-)]/2Dt = -d{ F + dF/dT *impl*[T(+) + T(-)] }/dz
!!
!!     It is discretized as follows:
!!
!!    PRHODJ(k)*PVARP(k)/PTSTEP
!!              = 
!!    PRHODJ(k)*PVARM(k)/PTSTEP 
!!  - (PRHODJ(k+1)+PRHODJ(k)  )/2. * PF(k+1)/PDZZ(k+1)
!!  + (PRHODJ(k)  +PRHODJ(k-1))/2. * PF(k)  /PDZZ(k)
!!  + (PRHODJ(k+1)+PRHODJ(k)  )/2. * 0.5*PIMPL* PDFDT(k+1) * PVARM(k+1)/PDZZ(k+1)
!!  - (PRHODJ(k+1)+PRHODJ(k)  )/2. * 0.5*PIMPL* PDFDT(k+1) * PVARP(k+1)/PDZZ(k+1)
!!  + (PRHODJ(k+1)+PRHODJ(k)  )/2. * 0.5*PIMPL* PDFDT(k+1) * PVARM(k)  /PDZZ(k+1)
!!  - (PRHODJ(k+1)+PRHODJ(k)  )/2. * 0.5*PIMPL* PDFDT(k+1) * PVARP(k)  /PDZZ(k+1)
!!  - (PRHODJ(k)  +PRHODJ(k-1))/2. * 0.5*PIMPL* PDFDT(k)   * PVARM(k)  /PDZZ(k)
!!  + (PRHODJ(k)  +PRHODJ(k-1))/2. * 0.5*PIMPL* PDFDT(k)   * PVARP(k)  /PDZZ(k)
!!  - (PRHODJ(k)  +PRHODJ(k-1))/2. * 0.5*PIMPL* PDFDT(k)   * PVARM(k-1)/PDZZ(k)
!!  + (PRHODJ(k)  +PRHODJ(k-1))/2. * 0.5*PIMPL* PDFDT(k)   * PVARP(k-1)/PDZZ(k)
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
!!  + (PRHODJ(k+1)+PRHODJ(k)  )/2. * 0.5*PIMPL* PDFDT(k+1) * PVARM(k+1)/PDZZ(k+1)
!!  + (PRHODJ(k+1)+PRHODJ(k)  )/2. * 0.5*PIMPL* PDFDT(k+1) * PVARM(k)  /PDZZ(k+1)
!!  - (PRHODJ(k)  +PRHODJ(k-1))/2. * 0.5*PIMPL* PDFDT(k)   * PVARM(k)  /PDZZ(k)
!!  - (PRHODJ(k)  +PRHODJ(k-1))/2. * 0.5*PIMPL* PDFDT(k)   * PVARM(k-1)/PDZZ(k)
!!
!!                      
!!        Then, the classical TRIDIAGonal algorithm is used to invert the 
!!     implicit operator. Its matrix is given by:
!!
!!     ( b(KKB)   c(KKB)      0        0        0         0        0        0  )
!!     ( a(KKB+1) b(KKB+1) c(KKB+1)    0  ...    0        0        0        0  ) 
!!     (   0      a(KKB+2) b(KKB+2) c(KKB+2).    0        0        0        0  ) 
!!      .......................................................................
!!     (   0   ...   0     a(k)     b(k)     c(k)         0   ...  0        0  ) 
!!      .......................................................................
!!     (   0         0        0        0        0 ...a(KKE-1) b(KKE-1) c(KKE-1))
!!     (   0         0        0        0        0 ...     0   a(KKE)   b(KKE)  )
!!
!!     KKB and KKE represent the first and the last inner mass levels of the
!!     model. The coefficients are:
!!         
!! a(k) = - (PRHODJ(k)  +PRHODJ(k-1))/2. * 0.5*PIMPL* PDFDT(k)  /PDZZ(k)
!! b(k) =    PRHODJ(k) / PTSTEP
!!        + (PRHODJ(k+1)+PRHODJ(k)  )/2. * 0.5*PIMPL* PDFDT(k+1)/PDZZ(k+1)
!!        - (PRHODJ(k)  +PRHODJ(k-1))/2. * 0.5*PIMPL* PDFDT(k)  /PDZZ(k)
!! c(k) = + (PRHODJ(k+1)+PRHODJ(k)  )/2. * 0.5*PIMPL* PDFDT(k+1)/PDZZ(k+1)
!!
!!          for all k /= KKB or KKE
!!
!!
!! b(KKB) =  PRHODJ(KKB) / PTSTEP
!!          +(PRHODJ(KKB+1)+PRHODJ(KKB))/2.*0.5*PIMPL*PDFDT(KKB+1)/PDZZ(KKB+1)
!! c(KKB) = +(PRHODJ(KKB+1)+PRHODJ(KKB))/2.*0.5*PIMPL*PDFDT(KKB+1)/PDZZ(KKB+1)
!!
!! b(KKE) =  PRHODJ(KKE) / PTSTEP
!!          -(PRHODJ(KKE)+PRHODJ(KKE-1))/2.*0.5*PIMPL*PDFDT(KKE)/PDZZ(KKE)
!! a(KKE) = -(PRHODJ(KKE)+PRHODJ(KKE-1))/2.*0.5*PIMPL*PDFDT(KKE)/PDZZ(KKE)
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
!!       V. Masson and S. Malardel         * Meteo-France *   
!! 
!!     MODIFICATIONS
!!     -------------
!!       Original        07/2006
!!       V.Masson : Optimization
!!       S. Riette Jan 2012: support for both order of vertical levels
!! ---------------------------------------------------------------------
!
!*       0. DECLARATIONS
!
USE MODD_DIMPHYEX,        ONLY: DIMPHYEX_t
!
USE MODI_SHUMAN_MF, ONLY: MZM_MF
!
IMPLICIT NONE
!
!
!*       0.1 declarations of arguments
!
TYPE(DIMPHYEX_t),       INTENT(IN)   :: D
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PVARM   ! variable at t-1      at mass point
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PF      ! flux in dT/dt=-dF/dz at flux point
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PDFDT   ! dF/dT                at flux point
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
REAL, DIMENSION(D%NIJT,D%NKT)  :: ZRHODJ_DFDT_O_DZ
REAL, DIMENSION(D%NIJT,D%NKT)  :: ZMZM_RHODJ
REAL, DIMENSION(D%NIJT,D%NKT)  :: ZA, ZB, ZC
REAL, DIMENSION(D%NIJT,D%NKT)  :: ZY ,ZGAM 
                                         ! RHS of the equation, 3D work array
REAL, DIMENSION(D%NIJT)                :: ZBET
                                         ! 2D work array
INTEGER                              :: JK, JI            ! loop counter
!
! ---------------------------------------------------------------------------
!                                              
!*      1.  Preliminaries
!           -------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('TRIDIAG_MASSFLUX',0,ZHOOK_HANDLE)
CALL MZM_MF(D, PRHODJ, ZMZM_RHODJ)
DO JK=1,D%NKT 
  DO JI=D%NIJB,D%NIJE 
    ZRHODJ_DFDT_O_DZ(JI,JK) = ZMZM_RHODJ(JI,JK)*PDFDT(JI,JK)/PDZZ(JI,JK)
  ENDDO
ENDDO
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
DO JI=D%NIJB,D%NIJE 
  ZY(JI,D%NKB) = PRHODJ(JI,D%NKB)*PVARM(JI,D%NKB)/PTSTEP             &
  - ZMZM_RHODJ(JI,D%NKB+D%NKL) * PF(JI,D%NKB+D%NKL)/PDZZ(JI,D%NKB+D%NKL)     &
  + ZMZM_RHODJ(JI,D%NKB) * PF(JI,D%NKB)/PDZZ(JI,D%NKB)     &
  + ZRHODJ_DFDT_O_DZ(JI,D%NKB+D%NKL) * 0.5*PIMPL * PVARM(JI,D%NKB+D%NKL)    &
  + ZRHODJ_DFDT_O_DZ(JI,D%NKB+D%NKL) * 0.5*PIMPL * PVARM(JI,D%NKB)
ENDDO
!
DO JK=1+D%NKTB,D%NKTE-1
  DO JI=D%NIJB,D%NIJE 
    ZY(JI,JK) = PRHODJ(JI,JK)*PVARM(JI,JK)/PTSTEP          &
    - ZMZM_RHODJ(JI,JK+D%NKL) * PF(JI,JK+D%NKL)/PDZZ(JI,JK+D%NKL)    &
    + ZMZM_RHODJ(JI,JK) * PF(JI,JK)/PDZZ(JI,JK)    &
    + ZRHODJ_DFDT_O_DZ(JI,JK+D%NKL) * 0.5*PIMPL * PVARM(JI,JK+D%NKL)  &
    + ZRHODJ_DFDT_O_DZ(JI,JK+D%NKL) * 0.5*PIMPL * PVARM(JI,JK)  &
    - ZRHODJ_DFDT_O_DZ(JI,JK) * 0.5*PIMPL * PVARM(JI,JK)  &
    - ZRHODJ_DFDT_O_DZ(JI,JK) * 0.5*PIMPL * PVARM(JI,JK-D%NKL)
  ENDDO
END DO
! 
IF (D%NKE==D%NKU) THEN
  DO JI=D%NIJB,D%NIJE 
    ZY(JI,D%NKE) = PRHODJ(JI,D%NKE)*PVARM(JI,D%NKE)/PTSTEP
  ENDDO
ELSE
  DO JI=D%NIJB,D%NIJE 
    ZY(JI,D%NKE) = PRHODJ(JI,D%NKE)*PVARM(JI,D%NKE)/PTSTEP &
    - ZMZM_RHODJ(JI,D%NKE+D%NKL) * PF(JI,D%NKE+D%NKL)/PDZZ(JI,D%NKE+D%NKL) &
    + ZMZM_RHODJ(JI,D%NKE) * PF(JI,D%NKE)/PDZZ(JI,D%NKE) &
    - ZRHODJ_DFDT_O_DZ(JI,D%NKE) * 0.5*PIMPL * PVARM(JI,D%NKE) &
    - ZRHODJ_DFDT_O_DZ(JI,D%NKE) * 0.5*PIMPL * PVARM(JI,D%NKE-D%NKL)
  ENDDO
ENDIF
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
  DO JI=D%NIJB,D%NIJE 
    ZB(JI,D%NKB) =   PRHODJ(JI,D%NKB)/PTSTEP                   &
    + ZRHODJ_DFDT_O_DZ(JI,D%NKB+D%NKL) * 0.5*PIMPL
    ZC(JI,D%NKB) =   ZRHODJ_DFDT_O_DZ(JI,D%NKB+D%NKL) * 0.5*PIMPL
  ENDDO

  DO JK=1+D%NKTB,D%NKTE-1
    DO JI=D%NIJB,D%NIJE 
      ZA(JI,JK) = - ZRHODJ_DFDT_O_DZ(JI,JK) * 0.5*PIMPL
      ZB(JI,JK) =   PRHODJ(JI,JK)/PTSTEP                   &
      + ZRHODJ_DFDT_O_DZ(JI,JK+D%NKL) * 0.5*PIMPL &
      - ZRHODJ_DFDT_O_DZ(JI,JK) * 0.5*PIMPL
      ZC(JI,JK) =   ZRHODJ_DFDT_O_DZ(JI,JK+D%NKL) * 0.5*PIMPL
    ENDDO
  END DO
  DO JI=D%NIJB,D%NIJE 
    ZA(JI,D%NKE) = - ZRHODJ_DFDT_O_DZ(JI,D%NKE) * 0.5*PIMPL
    ZB(JI,D%NKE) =   PRHODJ(JI,D%NKE)/PTSTEP                   &
    - ZRHODJ_DFDT_O_DZ(JI,D%NKE) * 0.5*PIMPL
  ENDDO
!
!*       3.2 going up
!            --------
!
  DO JI=D%NIJB,D%NIJE 
    ZBET(JI) = ZB(JI,D%NKB)  ! bet = b(D%NKB)
    PVARP(JI,D%NKB) = ZY(JI,D%NKB) / ZBET(JI)
  ENDDO

  !
  DO JK = D%NKB+D%NKL,D%NKE-D%NKL,D%NKL
    DO JI=D%NIJB,D%NIJE 
      ZGAM(JI,JK) = ZC(JI,JK-D%NKL) / ZBET(JI)
                                                    ! gam(k) = c(k-1) / bet
      ZBET(JI)    = ZB(JI,JK) - ZA(JI,JK) * ZGAM(JI,JK)
                                                    ! bet = b(k) - a(k)* gam(k)  
      PVARP(JI,JK)= ( ZY(JI,JK) - ZA(JI,JK) * PVARP(JI,JK-D%NKL) ) / ZBET(JI)
                                        ! res(k) = (y(k) -a(k)*res(k-1))/ bet 
    ENDDO
  END DO 
  DO JI=D%NIJB,D%NIJE 
  ! special treatment for the last level
    ZGAM(JI,D%NKE) = ZC(JI,D%NKE-D%NKL) / ZBET(JI)
                                                    ! gam(k) = c(k-1) / bet
    ZBET(JI)     = ZB(JI,D%NKE) - ZA(JI,D%NKE) * ZGAM(JI,D%NKE)
                                                    ! bet = b(k) - a(k)* gam(k)  
    PVARP(JI,D%NKE)= ( ZY(JI,D%NKE) - ZA(JI,D%NKE) * PVARP(JI,D%NKE-D%NKL) ) / &
    &ZBET(JI)
                                       ! res(k) = (y(k) -a(k)*res(k-1))/ bet 
  ENDDO
!
!*       3.3 going down
!            ----------
!
  DO JK = D%NKE-D%NKL,D%NKB,-D%NKL
    DO JI=D%NIJB,D%NIJE 
      PVARP(JI,JK) = PVARP(JI,JK) - ZGAM(JI,JK+D%NKL) * PVARP(JI,JK+D%NKL)
    ENDDO
  END DO
!
!
ELSE
  !!! EXPLICIT FORMULATION
  !
  DO JK=D%NKTB,D%NKTE
    DO JI=D%NIJB,D%NIJE 
      PVARP(JI,JK) = ZY(JI,JK) * PTSTEP / PRHODJ(JI,JK)
    ENDDO
  ENDDO
  !
END IF 
!
!
!*       4.  FILL THE UPPER AND LOWER EXTERNAL VALUES
!            ----------------------------------------
!
DO JI=D%NIJB,D%NIJE 
  PVARP(JI,D%NKA)=PVARP(JI,D%NKB)
  PVARP(JI,D%NKU)=PVARP(JI,D%NKE)
ENDDO
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('TRIDIAG_MASSFLUX',1,ZHOOK_HANDLE)
END SUBROUTINE TRIDIAG_MASSFLUX
END MODULE MODE_TRIDIAG_MASSFLUX
