!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     ###################
      MODULE MODI_TRIDIAG_MASSFLUX
!     ###################
INTERFACE
!
       SUBROUTINE TRIDIAG_MASSFLUX(KKA,KKB,KKE,KKU,KKL,PVARM,PF,PDFDT,PTSTEP,PIMPL,  &
                                 PDZZ,PRHODJ,PVARP             )
!
INTEGER,                INTENT(IN)   :: KKA          ! near ground array index
INTEGER,                INTENT(IN)   :: KKB          ! near ground physical index
INTEGER,                INTENT(IN)   :: KKE          ! uppest atmosphere physical index
INTEGER,                INTENT(IN)   :: KKU          ! uppest atmosphere array index
INTEGER,                INTENT(IN)   :: KKL          ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(:,:), INTENT(IN) :: PVARM   ! variable at t-1      at mass point
REAL, DIMENSION(:,:), INTENT(IN) :: PF      ! flux in dT/dt=-dF/dz at flux point
REAL, DIMENSION(:,:), INTENT(IN) :: PDFDT   ! dF/dT                at flux point
REAL,                   INTENT(IN) :: PTSTEP  ! Double time step
REAL,                   INTENT(IN) :: PIMPL   ! implicit weight
REAL, DIMENSION(:,:), INTENT(IN) :: PDZZ    ! Dz                   at flux point
REAL, DIMENSION(:,:), INTENT(IN) :: PRHODJ  ! (dry rho)*J          at mass point
!
REAL, DIMENSION(:,:), INTENT(OUT):: PVARP   ! variable at t+1      at mass point
!
END SUBROUTINE TRIDIAG_MASSFLUX
!
END INTERFACE
!
END MODULE MODI_TRIDIAG_MASSFLUX 


!      #################################################
       SUBROUTINE TRIDIAG_MASSFLUX(KKA,KKB,KKE,KKU,KKL,PVARM,PF,PDFDT,PTSTEP,PIMPL,  &
                                 PDZZ,PRHODJ,PVARP             )
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
USE MODD_PARAMETERS, ONLY: JPVEXT
USE MODI_SHUMAN_MF
!
IMPLICIT NONE
!
!
!*       0.1 declarations of arguments
!
INTEGER,                INTENT(IN)   :: KKA          ! near ground array index
INTEGER,                INTENT(IN)   :: KKB          ! near ground physical index
INTEGER,                INTENT(IN)   :: KKE          ! uppest atmosphere physical index
INTEGER,                INTENT(IN)   :: KKU          ! uppest atmosphere array index
INTEGER,                INTENT(IN)   :: KKL          ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(:,:), INTENT(IN) :: PVARM   ! variable at t-1      at mass point
REAL, DIMENSION(:,:), INTENT(IN) :: PF      ! flux in dT/dt=-dF/dz at flux point
REAL, DIMENSION(:,:), INTENT(IN) :: PDFDT   ! dF/dT                at flux point
REAL,                   INTENT(IN) :: PTSTEP  ! Double time step
REAL,                   INTENT(IN) :: PIMPL   ! implicit weight
REAL, DIMENSION(:,:), INTENT(IN) :: PDZZ    ! Dz                   at flux point
REAL, DIMENSION(:,:), INTENT(IN) :: PRHODJ  ! (dry rho)*J          at mass point
!
REAL, DIMENSION(:,:), INTENT(OUT):: PVARP   ! variable at t+1      at mass point
!
!
!*       0.2 declarations of local variables
!
REAL, DIMENSION(SIZE(PVARM,1),SIZE(PVARM,2))  :: ZRHODJ_DFDT_O_DZ
REAL, DIMENSION(SIZE(PVARM,1),SIZE(PVARM,2))  :: ZMZM_RHODJ
REAL, DIMENSION(SIZE(PVARM,1),SIZE(PVARM,2))  :: ZA, ZB, ZC
REAL, DIMENSION(SIZE(PVARM,1),SIZE(PVARM,2))  :: ZY ,ZGAM 
                                         ! RHS of the equation, 3D work array
REAL, DIMENSION(SIZE(PVARM,1))                :: ZBET
                                         ! 2D work array
INTEGER                              :: JK            ! loop counter
!
! ---------------------------------------------------------------------------
!                                              
!*      1.  Preliminaries
!           -------------
!
ZMZM_RHODJ  = MZM_MF(KKA,KKU,KKL,PRHODJ)
ZRHODJ_DFDT_O_DZ = ZMZM_RHODJ*PDFDT/PDZZ
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
ZY(:,KKB) = PRHODJ(:,KKB)*PVARM(:,KKB)/PTSTEP             &
    - ZMZM_RHODJ(:,KKB+KKL) * PF(:,KKB+KKL)/PDZZ(:,KKB+KKL)     &
    + ZMZM_RHODJ(:,KKB  ) * PF(:,KKB  )/PDZZ(:,KKB  )     &
    + ZRHODJ_DFDT_O_DZ(:,KKB+KKL) * 0.5*PIMPL * PVARM(:,KKB+KKL)    &
    + ZRHODJ_DFDT_O_DZ(:,KKB+KKL) * 0.5*PIMPL * PVARM(:,KKB  )
!
DO JK=2+JPVEXT,SIZE(ZY,2)-JPVEXT-1
  ZY(:,JK) = PRHODJ(:,JK)*PVARM(:,JK)/PTSTEP          &
    - ZMZM_RHODJ(:,JK+KKL) * PF(:,JK+KKL)/PDZZ(:,JK+KKL)    &
    + ZMZM_RHODJ(:,JK  ) * PF(:,JK  )/PDZZ(:,JK  )    &
    + ZRHODJ_DFDT_O_DZ(:,JK+KKL) * 0.5*PIMPL * PVARM(:,JK+KKL)  &
    + ZRHODJ_DFDT_O_DZ(:,JK+KKL) * 0.5*PIMPL * PVARM(:,JK  )  &
    - ZRHODJ_DFDT_O_DZ(:,JK  ) * 0.5*PIMPL * PVARM(:,JK  )  &
    - ZRHODJ_DFDT_O_DZ(:,JK  ) * 0.5*PIMPL * PVARM(:,JK-KKL)
END DO
! 
IF (JPVEXT==0) THEN
  ZY(:,KKE) = PRHODJ(:,KKE)*PVARM(:,KKE)/PTSTEP 
ELSE
  ZY(:,KKE) = PRHODJ(:,KKE)*PVARM(:,KKE)/PTSTEP &
   - ZMZM_RHODJ(:,KKE+KKL) * PF(:,KKE+KKL)/PDZZ(:,KKE+KKL) &
   + ZMZM_RHODJ(:,KKE  ) * PF(:,KKE  )/PDZZ(:,KKE  ) &
   - ZRHODJ_DFDT_O_DZ(:,KKE ) * 0.5*PIMPL * PVARM(:,KKE  ) &
   - ZRHODJ_DFDT_O_DZ(:,KKE ) * 0.5*PIMPL * PVARM(:,KKE-KKL)
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
  ZB(:,KKB) =   PRHODJ(:,KKB)/PTSTEP                   &
                + ZRHODJ_DFDT_O_DZ(:,KKB+KKL) * 0.5*PIMPL
  ZC(:,KKB) =   ZRHODJ_DFDT_O_DZ(:,KKB+KKL) * 0.5*PIMPL

  DO JK=2+JPVEXT,SIZE(ZY,2)-JPVEXT-1
    ZA(:,JK) = - ZRHODJ_DFDT_O_DZ(:,JK  ) * 0.5*PIMPL
    ZB(:,JK) =   PRHODJ(:,JK)/PTSTEP                   &
                 + ZRHODJ_DFDT_O_DZ(:,JK+KKL) * 0.5*PIMPL &
                 - ZRHODJ_DFDT_O_DZ(:,JK  ) * 0.5*PIMPL
    ZC(:,JK) =   ZRHODJ_DFDT_O_DZ(:,JK+KKL) * 0.5*PIMPL
  END DO

  ZA(:,KKE) = - ZRHODJ_DFDT_O_DZ(:,KKE  ) * 0.5*PIMPL
  ZB(:,KKE) =   PRHODJ(:,KKE)/PTSTEP                   &
                - ZRHODJ_DFDT_O_DZ(:,KKE  ) * 0.5*PIMPL
!
!*       3.2 going up
!            --------
!
  ZBET(:) = ZB(:,KKB)  ! bet = b(KKB)
  PVARP(:,KKB) = ZY(:,KKB) / ZBET(:)

  !
  DO JK = KKB+KKL,KKE-KKL,KKL
    ZGAM(:,JK) = ZC(:,JK-KKL) / ZBET(:)
                                                    ! gam(k) = c(k-1) / bet
    ZBET(:)    = ZB(:,JK) - ZA(:,JK) * ZGAM(:,JK)
                                                    ! bet = b(k) - a(k)* gam(k)  
    PVARP(:,JK)= ( ZY(:,JK) - ZA(:,JK) * PVARP(:,JK-KKL) ) / ZBET(:)
                                        ! res(k) = (y(k) -a(k)*res(k-1))/ bet 
  END DO 
  ! special treatment for the last level
  ZGAM(:,KKE) = ZC(:,KKE-KKL) / ZBET(:)
                                                    ! gam(k) = c(k-1) / bet
  ZBET(:)     = ZB(:,KKE) - ZA(:,KKE) * ZGAM(:,KKE)
                                                    ! bet = b(k) - a(k)* gam(k)  
  PVARP(:,KKE)= ( ZY(:,KKE) - ZA(:,KKE) * PVARP(:,KKE-KKL) ) / ZBET(:)
                                       ! res(k) = (y(k) -a(k)*res(k-1))/ bet 
!
!*       3.3 going down
!            ----------
!
  DO JK = KKE-KKL,KKB,-KKL
    PVARP(:,JK) = PVARP(:,JK) - ZGAM(:,JK+KKL) * PVARP(:,JK+KKL)
  END DO
!
!
ELSE
  !!! EXPLICIT FORMULATION
  !
  DO JK=1+JPVEXT,SIZE(PVARP,2)-JPVEXT
    PVARP(:,JK) = ZY(:,JK) * PTSTEP / PRHODJ(:,JK)
  ENDDO
  !
END IF 
!
!
!*       4.  FILL THE UPPER AND LOWER EXTERNAL VALUES
!            ----------------------------------------
!
PVARP(:,KKA)=PVARP(:,KKB)
PVARP(:,KKU)=PVARP(:,KKE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE TRIDIAG_MASSFLUX
