!MNH_LIC Copyright 2003-2020 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ###################
      MODULE MODI_TRIDIAG_THERMO
!     ###################
INTERFACE
!
       SUBROUTINE TRIDIAG_THERMO(KKA,KKU,KKL,PVARM,PF,PDFDDTDZ,PTSTEP,PIMPL,  &
                                 PDZZ,PRHODJ,PVARP             )
!
INTEGER,                INTENT(IN)   :: KKA           !near ground array index  
INTEGER,                INTENT(IN)   :: KKU           !uppest atmosphere array index
INTEGER,                INTENT(IN)   :: KKL           !vert. levels type 1=MNH -1=AR
REAL, DIMENSION(:,:,:), INTENT(IN) :: PVARM   ! variable at t-1      at mass point
REAL, DIMENSION(:,:,:), INTENT(IN) :: PF      ! flux in dT/dt=-dF/dz at flux point
REAL, DIMENSION(:,:,:), INTENT(IN) :: PDFDDTDZ! dF/d(dT/dz)          at flux point
REAL,                   INTENT(IN) :: PTSTEP  ! Double time step
REAL,                   INTENT(IN) :: PIMPL   ! implicit weight
REAL, DIMENSION(:,:,:), INTENT(IN) :: PDZZ    ! Dz                   at flux point
REAL, DIMENSION(:,:,:), INTENT(IN) :: PRHODJ  ! (dry rho)*J          at mass point
!
REAL, DIMENSION(:,:,:), INTENT(OUT):: PVARP   ! variable at t+1      at mass point
!
END SUBROUTINE TRIDIAG_THERMO
!
END INTERFACE
!
END MODULE MODI_TRIDIAG_THERMO 
!
!
!

!      #################################################
       SUBROUTINE TRIDIAG_THERMO(KKA,KKU,KKL,PVARM,PF,PDFDDTDZ,PTSTEP,PIMPL,  &
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
USE MODD_PARAMETERS, ONLY : JPVEXT_TURB
!
USE MODI_SHUMAN
!
IMPLICIT NONE
!
!
!*       0.1 declarations of arguments
!
INTEGER,              INTENT(IN)   :: KKA     !near ground array index  
INTEGER,              INTENT(IN)   :: KKU     !uppest atmosphere array index
INTEGER,              INTENT(IN)   :: KKL     !vert. levels type 1=MNH -1=ARO
REAL, DIMENSION(:,:,:), INTENT(IN) :: PVARM   ! variable at t-1      at mass point
REAL, DIMENSION(:,:,:), INTENT(IN) :: PF      ! flux in dT/dt=-dF/dz at flux point
REAL, DIMENSION(:,:,:), INTENT(IN) :: PDFDDTDZ! dF/d(dT/dz)          at flux point
REAL,                   INTENT(IN) :: PTSTEP  ! Double time step
REAL,                   INTENT(IN) :: PIMPL   ! implicit weight
REAL, DIMENSION(:,:,:), INTENT(IN) :: PDZZ    ! Dz                   at flux point
REAL, DIMENSION(:,:,:), INTENT(IN) :: PRHODJ  ! (dry rho)*J          at mass point
!
REAL, DIMENSION(:,:,:), INTENT(OUT):: PVARP   ! variable at t+1      at mass point
!
!
!*       0.2 declarations of local variables
!
REAL, DIMENSION(SIZE(PVARM,1),SIZE(PVARM,2),SIZE(PVARM,3))  :: ZRHODJ_DFDDTDZ_O_DZ2
REAL, DIMENSION(SIZE(PVARM,1),SIZE(PVARM,2),SIZE(PVARM,3))  :: ZMZM_RHODJ
REAL, DIMENSION(SIZE(PVARM,1),SIZE(PVARM,2),SIZE(PVARM,3))  :: ZA, ZB, ZC
REAL, DIMENSION(SIZE(PVARM,1),SIZE(PVARM,2),SIZE(PVARM,3))  :: ZY ,ZGAM 
                                         ! RHS of the equation, 3D work array
REAL, DIMENSION(SIZE(PVARM,1),SIZE(PVARM,2))                :: ZBET
                                         ! 2D work array
INTEGER             :: JK            ! loop counter
INTEGER             :: IKB,IKE       ! inner vertical limits
INTEGER             :: IKT          ! array size in k direction
INTEGER             :: IKTB,IKTE    ! start, end of k loops in physical domain 
!
! ---------------------------------------------------------------------------
!                                              
!*      1.  Preliminaries
!           -------------
!
IKT=SIZE(PVARM,3)          
IKTB=1+JPVEXT_TURB              
IKTE=IKT-JPVEXT_TURB
IKB=KKA+JPVEXT_TURB*KKL
IKE=KKU-JPVEXT_TURB*KKL

!
ZMZM_RHODJ  = MZM(PRHODJ)
ZRHODJ_DFDDTDZ_O_DZ2 = ZMZM_RHODJ*PDFDDTDZ/PDZZ**2
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
ZY(:,:,IKB) = PRHODJ(:,:,IKB)*PVARM(:,:,IKB)/PTSTEP                  &
    - ZMZM_RHODJ(:,:,IKB+KKL) * PF(:,:,IKB+KKL)/PDZZ(:,:,IKB+KKL)    &
    + ZMZM_RHODJ(:,:,IKB  ) * PF(:,:,IKB  )/PDZZ(:,:,IKB  )          &
    + ZRHODJ_DFDDTDZ_O_DZ2(:,:,IKB+KKL) * PIMPL * PVARM(:,:,IKB+KKL) &
    - ZRHODJ_DFDDTDZ_O_DZ2(:,:,IKB+KKL) * PIMPL * PVARM(:,:,IKB  )
!
  ZY(:,:,IKTB+1:IKTE-1) = PRHODJ(:,:,IKTB+1:IKTE-1)*PVARM(:,:,IKTB+1:IKTE-1)/PTSTEP                 &
    - ZMZM_RHODJ(:,:,IKTB+1+KKL:IKTE-1+KKL) * PF(:,:,IKTB+1+KKL:IKTE-1+KKL)/PDZZ(:,:,IKTB+1+KKL:IKTE-1+KKL)     &
    + ZMZM_RHODJ(:,:,IKTB+1:IKTE-1  ) * PF(:,:,IKTB+1:IKTE-1  )/PDZZ(:,:,IKTB+1:IKTE-1  )           &
    + ZRHODJ_DFDDTDZ_O_DZ2(:,:,IKTB+1+KKL:IKTE-1+KKL) * PIMPL * PVARM(:,:,IKTB+1+KKL:IKTE-1+KKL) &
    - ZRHODJ_DFDDTDZ_O_DZ2(:,:,IKTB+1+KKL:IKTE-1+KKL) * PIMPL * PVARM(:,:,IKTB+1:IKTE-1  )   &
    - ZRHODJ_DFDDTDZ_O_DZ2(:,:,IKTB+1:IKTE-1    ) * PIMPL * PVARM(:,:,IKTB+1:IKTE-1  )   &
    + ZRHODJ_DFDDTDZ_O_DZ2(:,:,IKTB+1:IKTE-1    ) * PIMPL * PVARM(:,:,IKTB+1-KKL:IKTE-1-KKL)
! 
ZY(:,:,IKE) = PRHODJ(:,:,IKE)*PVARM(:,:,IKE)/PTSTEP               &
    - ZMZM_RHODJ(:,:,IKE+KKL) * PF(:,:,IKE+KKL)/PDZZ(:,:,IKE+KKL) &
    + ZMZM_RHODJ(:,:,IKE  ) * PF(:,:,IKE  )/PDZZ(:,:,IKE  )       &
    - ZRHODJ_DFDDTDZ_O_DZ2(:,:,IKE ) * PIMPL * PVARM(:,:,IKE  )   &
    + ZRHODJ_DFDDTDZ_O_DZ2(:,:,IKE ) * PIMPL * PVARM(:,:,IKE-KKL)
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
  ZB(:,:,IKB) =   PRHODJ(:,:,IKB)/PTSTEP                   &
                - ZRHODJ_DFDDTDZ_O_DZ2(:,:,IKB+KKL) * PIMPL
  ZC(:,:,IKB) =   ZRHODJ_DFDDTDZ_O_DZ2(:,:,IKB+KKL) * PIMPL
!
  ZA(:,:,IKTB+1:IKTE-1) =   ZRHODJ_DFDDTDZ_O_DZ2(:,:,IKTB+1:IKTE-1) * PIMPL
  ZB(:,:,IKTB+1:IKTE-1) =   PRHODJ(:,:,IKTB+1:IKTE-1)/PTSTEP                        &
                          - ZRHODJ_DFDDTDZ_O_DZ2(:,:,IKTB+1+KKL:IKTE-1+KKL) * PIMPL &
                          - ZRHODJ_DFDDTDZ_O_DZ2(:,:,IKTB+1:IKTE-1) * PIMPL
  ZC(:,:,IKTB+1:IKTE-1) =   ZRHODJ_DFDDTDZ_O_DZ2(:,:,IKTB+1+KKL:IKTE-1+KKL) * PIMPL
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
  DO JK = IKB+KKL,IKE-KKL,KKL
    ZGAM(:,:,JK) = ZC(:,:,JK-KKL) / ZBET(:,:)  
                                                    ! gam(k) = c(k-1) / bet
    ZBET(:,:)    = ZB(:,:,JK) - ZA(:,:,JK) * ZGAM(:,:,JK)
                                                    ! bet = b(k) - a(k)* gam(k)  
    PVARP(:,:,JK)= ( ZY(:,:,JK) - ZA(:,:,JK) * PVARP(:,:,JK-KKL) ) / ZBET(:,:)
                                        ! res(k) = (y(k) -a(k)*res(k-1))/ bet 
  END DO 
  ! special treatment for the last level
  ZGAM(:,:,IKE) = ZC(:,:,IKE-KKL) / ZBET(:,:) 
                                                    ! gam(k) = c(k-1) / bet
  ZBET(:,:)     = ZB(:,:,IKE) - ZA(:,:,IKE) * ZGAM(:,:,IKE)
                                                    ! bet = b(k) - a(k)* gam(k)  
  PVARP(:,:,IKE)= ( ZY(:,:,IKE) - ZA(:,:,IKE) * PVARP(:,:,IKE-KKL) ) / ZBET(:,:)
                                       ! res(k) = (y(k) -a(k)*res(k-1))/ bet 
!
!*       3.3 going down
!            ----------
!
  DO JK = IKE-KKL,IKB,-1*KKL
    PVARP(:,:,JK) = PVARP(:,:,JK) - ZGAM(:,:,JK+KKL) * PVARP(:,:,JK+KKL)
  END DO
!
ELSE
! 
  PVARP(:,:,IKTB:IKTE) = ZY(:,:,IKTB:IKTE) * PTSTEP / PRHODJ(:,:,IKTB:IKTE)
!
END IF 
!
!
!*       4.  FILL THE UPPER AND LOWER EXTERNAL VALUES
!            ----------------------------------------
!
PVARP(:,:,KKA)=PVARP(:,:,IKB)
PVARP(:,:,KKU)=PVARP(:,:,IKE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE TRIDIAG_THERMO
