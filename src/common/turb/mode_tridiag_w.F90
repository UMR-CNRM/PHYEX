!MNH_LIC Copyright 2011-2020 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ###################
      MODULE MODE_TRIDIAG_W
!     ###################
IMPLICIT NONE
CONTAINS
!
!
!

!      #################################################
       SUBROUTINE TRIDIAG_W(PVARM,PF,PDFDDWDZ,PTSTEP, &
                                 PMZF_DZZ,PRHODJ,PVARP)
!      #################################################
!
!
!!****   *TRIDIAG_W* - routine to solve a time implicit scheme
!!
!!
!!     PURPOSE
!!     -------
!        The purpose of this routine is to give a field PVARP at t+1, by 
!      solving an implicit TRIDIAGonal system obtained by the 
!      discretization of the vertical turbulent diffusion. 
!      The function of F(dT/dz) must have been linearized.
!      PVARP is localized at a flux point.
!
!!**   METHOD
!!     ------
!!
!!        [W(+) - W(-)]/Dt = -d{ F + dF/d(dW/dz) * [dW/dz(+)-dW/dz(-)] }/dz
!!
!!     It is discretized as follows:
!!
!!    (PRHODJ(k)+PRHODJ(k-1))/2.*PVARP(k)/PTSTEP
!!              = 
!!    (PRHODJ(k)+PRHODJ(k-1))/2.*PVARM(k)/PTSTEP 
!!  - PRHODJ(k)   * PF(k  )/PMZF_PDZZ(k  )
!!  + PRHODJ(k-1) * PF(k-1)/PMZF_PDZZ(k-1)
!!  - PRHODJ(k)   * PDFDDWDZ(k)   * PVARP(k+1)/PMZF_DZZ(k)**2
!!  + PRHODJ(k)   * PDFDDWDZ(k)   * PVARP(k)  /PMZF_DZZ(k)**2
!!  + PRHODJ(k-1) * PDFDDWDZ(k-1) * PVARP(k)  /PMZF_DZZ(k-1)**2
!!  - PRHODJ(k-1) * PDFDDWDZ(k-1) * PVARP(k-1)/PMZF_DZZ(k-1)**2
!!  + PRHODJ(k)   * PDFDDWDZ(k)   * PVARM(k+1)/PMZF_DZZ(k)**2
!!  - PRHODJ(k)   * PDFDDWDZ(k)   * PVARM(k)  /PMZF_DZZ(k)**2
!!  - PRHODJ(k-1) * PDFDDWDZ(k-1) * PVARM(k)  /PMZF_DZZ(k-1)**2
!!  + PRHODJ(k-1) * PDFDDWDZ(k-1) * PVARM(k-1)/PMZF_DZZ(k-1)**2
!!
!!
!!    The system to solve is:
!!
!!      A*PVARP(k-1) + B*PVARP(k) + C*PVARP(k+1) = Y(k)
!!
!!
!!    The RHS of the linear system in PVARP writes:
!!
!! y(k)    = (PRHODJ(k)+PRHODJ(k-1))/2.*PVARM(k)/PTSTEP
!!        - PRHODJ(k)   * PF(k  )/PMZF_PDZZ(k  )
!!        + PRHODJ(k-1) * PF(k-1)/PMZF_PDZZ(k-1)
!!        + PRHODJ(k)   * PDFDDWDZ(k)   * PVARM(k+1)/PMZF_DZZ(k)**2
!!        - PRHODJ(k)   * PDFDDWDZ(k)   * PVARM(k)  /PMZF_DZZ(k)**2
!!        - PRHODJ(k-1) * PDFDDWDZ(k-1) * PVARM(k)  /PMZF_DZZ(k-1)**2
!!        + PRHODJ(k-1) * PDFDDWDZ(k-1) * PVARM(k-1)/PMZF_DZZ(k-1)**2
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
!! a(k) = + PRHODJ(k-1) * PDFDDWDZ(k-1)/PMZF_DZZ(k-1)**2
!! b(k) =   (PRHODJ(k)+PRHODJ(k-1))/2. / PTSTEP
!!        - PRHODJ(k)   * PDFDDWDZ(k)  /PMZF_DZZ(k)**2
!!        - PRHODJ(k-1) * PDFDDWDZ(k-1)/PMZF_DZZ(k-1)**2
!! c(k) = + PRHODJ(k)   * PDFDDWDZ(k)/PMZF_DZZ(k)**2
!!
!!          for all k /= ikb or ike
!!
!! with the boundary conditions:
!!      PVARP(ikb-1) = PVARP(ikb)
!!      PVARP(ike+1) = 0.
!!
!! This induces:
!!
!! b(ikb) =   (PRHODJ(ikb)+PRHODJ(ikb-1))/2. / PTSTEP
!!          - PRHODJ(ikb)   * PDFDDWDZ(ikb)  /PMZF_DZZ(ikb)**2
!! c(ikb) = + PRHODJ(ikb)   * PDFDDWDZ(ikb)/PMZF_DZZ(ikb)**2
!!
!! b(ike) =  (PRHODJ(ike)+PRHODJ(ike-1))/2. / PTSTEP
!!        - PRHODJ(ike-1) * PDFDDWDZ(ike-1)/PMZF_DZZ(ike-1)**2
!!        - PRHODJ(ike  ) * PDFDDWDZ(ike  )/PMZF_DZZ(ike  )**2
!! a(ike) = + PRHODJ(ike-1) * PDFDDWDZ(ike-1)/PMZF_DZZ(ike-1)**2
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
!!       Original        04/2011 (from tridiag_thermo.f90)
!!                       03/2014 modification of upper boundary condition
!! ---------------------------------------------------------------------
!
!*       0. DECLARATIONS
!
USE MODD_PARAMETERS, ONLY : JPVEXT
!
USE MODI_SHUMAN
!
IMPLICIT NONE
!
!
!*       0.1 declarations of arguments
!
REAL, DIMENSION(:,:,:), INTENT(IN) :: PVARM   ! variable at t-1      at flux point
REAL, DIMENSION(:,:,:), INTENT(IN) :: PF      ! flux in dT/dt=-dF/dz at mass point
REAL, DIMENSION(:,:,:), INTENT(IN) :: PDFDDWDZ! dF/d(dW/dz)          at mass point
REAL,                   INTENT(IN) :: PTSTEP  ! Double time step
REAL, DIMENSION(:,:,:), INTENT(IN) :: PMZF_DZZ! Dz                   at mass point
REAL, DIMENSION(:,:,:), INTENT(IN) :: PRHODJ  ! (dry rho)*J          at mass point
!
REAL, DIMENSION(:,:,:), INTENT(OUT):: PVARP   ! variable at t+1      at flux point
!
!
!*       0.2 declarations of local variables
!
REAL, DIMENSION(SIZE(PVARM,1),SIZE(PVARM,2),SIZE(PVARM,3))  :: ZRHODJ_DFDDWDZ_O_DZ2
REAL, DIMENSION(SIZE(PVARM,1),SIZE(PVARM,2),SIZE(PVARM,3))  :: ZMZM_RHODJ
REAL, DIMENSION(SIZE(PVARM,1),SIZE(PVARM,2),SIZE(PVARM,3))  :: ZA, ZB, ZC
REAL, DIMENSION(SIZE(PVARM,1),SIZE(PVARM,2),SIZE(PVARM,3))  :: ZY ,ZGAM 
                                         ! RHS of the equation, 3D work array
REAL, DIMENSION(SIZE(PVARM,1),SIZE(PVARM,2))                :: ZBET
                                         ! 2D work array
INTEGER                              :: JK            ! loop counter
INTEGER                              :: IKB,IKE       ! inner vertical limits
!
! ---------------------------------------------------------------------------
!                                              
!*      1.  Preliminaries
!           -------------
!
IKB=1+JPVEXT
IKE=SIZE(PVARM,3)-JPVEXT 
!
ZMZM_RHODJ  = MZM(PRHODJ)
ZRHODJ_DFDDWDZ_O_DZ2 = PRHODJ*PDFDDWDZ/PMZF_DZZ**2
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
!! y(k)    = (PRHODJ(k)+PRHODJ(k-1))/2.*PVARM(k)/PTSTEP
!!        - PRHODJ(k)   * PF(k  )/PMZF_PDZZ(k  )
!!        + PRHODJ(k-1) * PF(k-1)/PMZF_PDZZ(k-1)
!!        + PRHODJ(k)   * PDFDDWDZ(k)   * PVARM(k+1)/PMZF_DZZ(k)**2
!!        - PRHODJ(k)   * PDFDDWDZ(k)   * PVARM(k)  /PMZF_DZZ(k)**2
!!        - PRHODJ(k-1) * PDFDDWDZ(k-1) * PVARM(k)  /PMZF_DZZ(k-1)**2
!!        + PRHODJ(k-1) * PDFDDWDZ(k-1) * PVARM(k-1)/PMZF_DZZ(k-1)**2
!
ZY(:,:,IKB) = ZMZM_RHODJ(:,:,IKB)*PVARM(:,:,IKB)/PTSTEP              &
    - PRHODJ(:,:,IKB  ) * PF(:,:,IKB  )/PMZF_DZZ(:,:,IKB  )           &
    + PRHODJ(:,:,IKB-1) * PF(:,:,IKB-1)/PMZF_DZZ(:,:,IKB-1)           &
    + ZRHODJ_DFDDWDZ_O_DZ2(:,:,IKB) * PVARM(:,:,IKB+1)&
    - ZRHODJ_DFDDWDZ_O_DZ2(:,:,IKB) * PVARM(:,:,IKB  )
!
  ZY(:,:,IKB+1:IKE-1) = ZMZM_RHODJ(:,:,IKB+1:IKE-1)*PVARM(:,:,IKB+1:IKE-1)/PTSTEP               &
    - PRHODJ(:,:,IKB+1:IKE-1  ) * PF(:,:,IKB+1:IKE-1  )/PMZF_DZZ(:,:,IKB+1:IKE-1  )              &
    + PRHODJ(:,:,IKB:IKE-2) * PF(:,:,IKB:IKE-2)/PMZF_DZZ(:,:,IKB:IKE-2)              &
    + ZRHODJ_DFDDWDZ_O_DZ2(:,:,IKB+1:IKE-1  ) * PVARM(:,:,IKB+2:IKE)  &
    - ZRHODJ_DFDDWDZ_O_DZ2(:,:,IKB+1:IKE-1  ) * PVARM(:,:,IKB+1:IKE-1  )  &
    - ZRHODJ_DFDDWDZ_O_DZ2(:,:,IKB:IKE-2) * PVARM(:,:,IKB+1:IKE-1  )  &
    + ZRHODJ_DFDDWDZ_O_DZ2(:,:,IKB:IKE-2) * PVARM(:,:,IKB:IKE-2)
! 
ZY(:,:,IKE) = ZMZM_RHODJ(:,:,IKE)*PVARM(:,:,IKE)/PTSTEP              &
    - PRHODJ(:,:,IKE  ) * PF(:,:,IKE  )/PMZF_DZZ(:,:,IKE  )           &
    + PRHODJ(:,:,IKE-1) * PF(:,:,IKE-1)/PMZF_DZZ(:,:,IKE-1)           &
    - ZRHODJ_DFDDWDZ_O_DZ2(:,:,IKE  ) * PVARM(:,:,IKE )  &
    - ZRHODJ_DFDDWDZ_O_DZ2(:,:,IKE-1) * PVARM(:,:,IKE  ) &
    + ZRHODJ_DFDDWDZ_O_DZ2(:,:,IKE-1) * PVARM(:,:,IKE-1)
!
!
!*       3.  INVERSION OF THE TRIDIAGONAL SYSTEM
!            -----------------------------------
!
!
!*       3.1 arrays A, B, C
!            --------------
!
!! a(k) = + PRHODJ(k-1) * PDFDDWDZ(k-1)/PMZF_DZZ(k-1)**2
!! b(k) =   (PRHODJ(k)+PRHODJ(k-1))/2. / PTSTEP
!!        - PRHODJ(k)   * PDFDDWDZ(k)  /PMZF_DZZ(k)**2
!!        - PRHODJ(k-1) * PDFDDWDZ(k-1)/PMZF_DZZ(k-1)**2
!! c(k) = + PRHODJ(k)   * PDFDDWDZ(k)/PMZF_DZZ(k)**2
!
  ZB(:,:,IKB) =   ZMZM_RHODJ(:,:,IKB)/PTSTEP      &
                - ZRHODJ_DFDDWDZ_O_DZ2(:,:,IKB)
  ZC(:,:,IKB) =   ZRHODJ_DFDDWDZ_O_DZ2(:,:,IKB)

    ZA(:,:,IKB+1:IKE-1) =   ZRHODJ_DFDDWDZ_O_DZ2(:,:,IKB:IKE-2)
    ZB(:,:,IKB+1:IKE-1) =   ZMZM_RHODJ(:,:,IKB+1:IKE-1)/PTSTEP      &
                 - ZRHODJ_DFDDWDZ_O_DZ2(:,:,IKB+1:IKE-1  ) &
                 - ZRHODJ_DFDDWDZ_O_DZ2(:,:,IKB:IKE-2)
    ZC(:,:,IKB+1:IKE-1) =   ZRHODJ_DFDDWDZ_O_DZ2(:,:,IKB+1:IKE-1  )

  ZA(:,:,IKE) =   ZRHODJ_DFDDWDZ_O_DZ2(:,:,IKE-1)
  ZB(:,:,IKE) =   ZMZM_RHODJ(:,:,IKE)/PTSTEP      &
                - ZRHODJ_DFDDWDZ_O_DZ2(:,:,IKE  ) &
                - ZRHODJ_DFDDWDZ_O_DZ2(:,:,IKE-1)
!
!*       3.2 going up
!            --------
!
  ZBET(:,:) = ZB(:,:,IKB)  ! bet = b(ikb)
  PVARP(:,:,IKB) = ZY(:,:,IKB) / ZBET(:,:)

  !
  DO JK = IKB+1,IKE-1
    ZGAM(:,:,JK) = ZC(:,:,JK-1) / ZBET(:,:)  
                                                    ! gam(k) = c(k-1) / bet
    ZBET(:,:)    = ZB(:,:,JK) - ZA(:,:,JK) * ZGAM(:,:,JK)
                                                    ! bet = b(k) - a(k)* gam(k)  
    PVARP(:,:,JK)= ( ZY(:,:,JK) - ZA(:,:,JK) * PVARP(:,:,JK-1) ) / ZBET(:,:)
                                        ! res(k) = (y(k) -a(k)*res(k-1))/ bet 
  END DO 
  ! special treatment for the last level
  ZGAM(:,:,IKE) = ZC(:,:,IKE-1) / ZBET(:,:) 
                                                    ! gam(k) = c(k-1) / bet
  ZBET(:,:)     = ZB(:,:,IKE) - ZA(:,:,IKE) * ZGAM(:,:,IKE)
                                                    ! bet = b(k) - a(k)* gam(k)  
  PVARP(:,:,IKE)= ( ZY(:,:,IKE) - ZA(:,:,IKE) * PVARP(:,:,IKE-1) ) / ZBET(:,:)
                                       ! res(k) = (y(k) -a(k)*res(k-1))/ bet 
!
!*       3.3 going down
!            ----------
!
  DO JK = IKE-1,IKB,-1
    PVARP(:,:,JK) = PVARP(:,:,JK) - ZGAM(:,:,JK+1) * PVARP(:,:,JK+1)
  END DO
!
!
!*       4.  FILL THE UPPER AND LOWER EXTERNAL VALUES
!            ----------------------------------------
!
PVARP(:,:,IKB-1)=PVARP(:,:,IKB)
PVARP(:,:,IKE+1)=0.
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE TRIDIAG_W
!
END MODULE MODE_TRIDIAG_W
