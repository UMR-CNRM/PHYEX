!MNH_LIC Copyright 2011-2025 CNRS, Meteo-France and Universite Paul Sabatier
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
       SUBROUTINE TRIDIAG_W(D,PVARM,PF,PDFDDWDZ,PTSTEP, &
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
!!       M.Moge          04/2016 Use openACC directives to port the TURB part of Meso-NH on GPU
!! ---------------------------------------------------------------------
!
!*       0. DECLARATIONS
!
USE MODD_DIMPHYEX,   ONLY : DIMPHYEX_t
USE MODD_PARAMETERS, ONLY : JPVEXT
USE MODD_DIMPHYEX,   ONLY: DIMPHYEX_t
USE MODE_SHUMAN_PHY, ONLY: MZM_PHY

#ifdef MNH_COMPILER_CCE
!$mnh_undef(LOOP)
!$mnh_undef(OPENACC)
#endif
!
IMPLICIT NONE
!
!
!*       0.1 declarations of arguments
!
TYPE(DIMPHYEX_t),         INTENT(IN)    :: D
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PVARM   ! variable at t-1      at flux point
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PF      ! flux in dT/dt=-dF/dz at mass point
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PDFDDWDZ! dF/d(dW/dz)          at mass point
REAL,                   INTENT(IN) :: PTSTEP  ! Double time step
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PMZF_DZZ! Dz                   at mass point
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PRHODJ  ! (dry rho)*J          at mass point
!
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT):: PVARP   ! variable at t+1      at flux point
!
!*       0.2 declarations of local variables
!
REAL, DIMENSION(D%NIJT,D%NKT) :: ZRHODJ_DFDDWDZ_O_DZ2
REAL, DIMENSION(D%NIJT,D%NKT) :: ZMZM_RHODJ
REAL, DIMENSION(D%NIJT,D%NKT) :: ZA, ZB, ZC
REAL, DIMENSION(D%NIJT,D%NKT) :: ZY ,ZGAM ! RHS of the equation, 3D work array
REAL, DIMENSION(D%NIJT)       :: ZBET     ! 2D work array
!
INTEGER                              :: JK,JIJ        ! loop counter
INTEGER                              :: IKB,IKE       ! inner vertical limits
!
INTEGER  :: IIJT, IKT
! ---------------------------------------------------------------------------




!
!*      1.  Preliminaries
!           -------------
!
IKB=1+JPVEXT
IKE=SIZE(PVARM,2)-JPVEXT 
IIJT=D%NIJT
IKT=D%NKT
!
CALL MZM_PHY(D, PRHODJ, ZMZM_RHODJ)


DO JK=1, IKT
     DO JIJ=1, IIJT
       ZRHODJ_DFDDWDZ_O_DZ2(JIJ, JK) = PRHODJ(JIJ, JK)*PDFDDWDZ(JIJ, JK)/PMZF_DZZ(JIJ, JK)**2
  END DO
   END DO

!

ZA=0.
ZB=0.
ZC=0.
ZY=0.

!
! acc wait
!
!
!*      2.  COMPUTE THE RIGHT HAND SIDE
!           ---------------------------
!
!! y(k)    = (PRHODJ(k)+PRHODJ(k-1))/2.*PVARM(k)/PTSTEP
!!        - PRHODJ(k)   * PF(k  )/PMZF_PDZZ(k  )
!!        + PRHODJ(k-1) * PF(k-1)/PMZF_PDZZ(k-1)
!!#ifndef MNH_BITREP
!!        + PRHODJ(k)   * PDFDDWDZ(k)   * PVARM(k+1)/PMZF_DZZ(k)**2
!!        - PRHODJ(k)   * PDFDDWDZ(k)   * PVARM(k)  /PMZF_DZZ(k)**2
!!        - PRHODJ(k-1) * PDFDDWDZ(k-1) * PVARM(k)  /PMZF_DZZ(k-1)**2
!!        + PRHODJ(k-1) * PDFDDWDZ(k-1) * PVARM(k-1)/PMZF_DZZ(k-1)**2
!!#else
!!        + PRHODJ(k)   * PDFDDWDZ(k)   * PVARM(k+1)/BR_P2(PMZF_DZZ(k))
!!        - PRHODJ(k)   * PDFDDWDZ(k)   * PVARM(k)  /BR_P2(PMZF_DZZ(k))
!!        - PRHODJ(k-1) * PDFDDWDZ(k-1) * PVARM(k)  /BR_P2(PMZF_DZZ(k-1))
!!        + PRHODJ(k-1) * PDFDDWDZ(k-1) * PVARM(k-1)/BR_P2(PMZF_DZZ(k-1))
!!#endif
!

!$mnh_do_concurrent(JIJ=1:IIJT)
DO JIJ=1,IIJT
   ZY(JIJ,IKB) = ZMZM_RHODJ(JIJ,IKB)*PVARM(JIJ,IKB)/PTSTEP              &
        - PRHODJ(JIJ,IKB  ) * PF(JIJ,IKB  )/PMZF_DZZ(JIJ,IKB  )           &
        + PRHODJ(JIJ,IKB-1) * PF(JIJ,IKB-1)/PMZF_DZZ(JIJ,IKB-1)           &
        + ZRHODJ_DFDDWDZ_O_DZ2(JIJ,IKB) * PVARM(JIJ,IKB+1)&
        - ZRHODJ_DFDDWDZ_O_DZ2(JIJ,IKB) * PVARM(JIJ,IKB  )
END DO
!$mnh_end_do()

!

!$mnh_do_concurrent(JIJ=1:IIJT,JK=IKB+1:IKE-1)
DO JK=IKB+1,IKE-1
  DO JIJ=1,IIJT
    ZY(JIJ,JK) = ZMZM_RHODJ(JIJ,JK)*PVARM(JIJ,JK)/PTSTEP               &
       - PRHODJ(JIJ,JK  ) * PF(JIJ,JK  )/PMZF_DZZ(JIJ,JK  )          &
       + PRHODJ(JIJ,JK-1) * PF(JIJ,JK-1)/PMZF_DZZ(JIJ,JK-1)              &
       + ZRHODJ_DFDDWDZ_O_DZ2(JIJ,JK  ) * PVARM(JIJ,JK+1)  &
       - ZRHODJ_DFDDWDZ_O_DZ2(JIJ,JK  ) * PVARM(JIJ,JK  )  &
       - ZRHODJ_DFDDWDZ_O_DZ2(JIJ,JK-1) * PVARM(JIJ,JK  )  &
       + ZRHODJ_DFDDWDZ_O_DZ2(JIJ,JK-1) * PVARM(JIJ,JK-1)
  END DO
END DO
!$mnh_end_do()  

! 

!$mnh_do_concurrent(JIJ=1:IIJT)
DO JIJ=1,IIJT
   ZY(JIJ,IKE) = ZMZM_RHODJ(JIJ,IKE)*PVARM(JIJ,IKE)/PTSTEP              &
        - PRHODJ(JIJ,IKE  ) * PF(JIJ,IKE  )/PMZF_DZZ(JIJ,IKE  )           &
        + PRHODJ(JIJ,IKE-1) * PF(JIJ,IKE-1)/PMZF_DZZ(JIJ,IKE-1)           &
        - ZRHODJ_DFDDWDZ_O_DZ2(JIJ,IKE  ) * PVARM(JIJ,IKE )  &
        - ZRHODJ_DFDDWDZ_O_DZ2(JIJ,IKE-1) * PVARM(JIJ,IKE  ) &
        + ZRHODJ_DFDDWDZ_O_DZ2(JIJ,IKE-1) * PVARM(JIJ,IKE-1)
END DO
!$mnh_end_do()

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

!$mnh_do_concurrent(JIJ=1:IIJT)
DO JIJ=1,IIJT
  ZB(JIJ,IKB) =   ZMZM_RHODJ(JIJ,IKB)/PTSTEP      &
       - ZRHODJ_DFDDWDZ_O_DZ2(JIJ,IKB)
END DO
!$mnh_end_do()  


!$mnh_do_concurrent(JIJ=1:IIJT)
DO JIJ=1,IIJT
   ZC(JIJ,IKB) =   ZRHODJ_DFDDWDZ_O_DZ2(JIJ,IKB)
END DO
!$mnh_end_do()   



!$mnh_do_concurrent(JIJ=1:IIJT,JK=IKB+1:IKE-1)
DO JK=IKB+1,IKE-1
  DO JIJ=1,IIJT
    ZA(JIJ,JK) =   ZRHODJ_DFDDWDZ_O_DZ2(JIJ,JK-1)
    ZB(JIJ,JK) =   ZMZM_RHODJ(JIJ,JK)/PTSTEP      &
                 - ZRHODJ_DFDDWDZ_O_DZ2(JIJ,JK  ) &
                 - ZRHODJ_DFDDWDZ_O_DZ2(JIJ,JK-1)
    ZC(JIJ,JK) =   ZRHODJ_DFDDWDZ_O_DZ2(JIJ,JK  )
  END DO
END DO
!$mnh_end_do()    



!$mnh_do_concurrent(JIJ=1:IIJT)
DO JIJ=1,IIJT
   ZA(JIJ,IKE) =   ZRHODJ_DFDDWDZ_O_DZ2(JIJ,IKE-1)
END DO
!$mnh_end_do()   


!$mnh_do_concurrent(JIJ=1:IIJT)
DO JIJ=1,IIJT
  ZB(JIJ,IKE) =   ZMZM_RHODJ(JIJ,IKE)/PTSTEP      &
                - ZRHODJ_DFDDWDZ_O_DZ2(JIJ,IKE  ) &
                - ZRHODJ_DFDDWDZ_O_DZ2(JIJ,IKE-1)
END DO
!$mnh_end_do()  

!
!
! acc wait
!
!*       3.2 going up
!            --------
!

!$mnh_do_concurrent(JIJ=1:IIJT)
DO JIJ=1,IIJT
  ZBET(JIJ) = ZB(JIJ,IKB)  ! bet = b(ikb)
  PVARP(JIJ,IKB) = ZY(JIJ,IKB) / ZBET(JIJ)
END DO
!$mnh_end_do()

!


DO JK = IKB+1,IKE-1
#ifdef MNH_COMPILER_CCE
   
#endif
   !$mnh_do_concurrent(JIJ=1:IIJT)
  DO JIJ=1,IIJT
      ZGAM(JIJ,JK) = ZC(JIJ,JK-1) / ZBET(JIJ)  
      ! gam(k) = c(k-1) / bet
      ZBET(JIJ)    = ZB(JIJ,JK) - ZA(JIJ,JK) * ZGAM(JIJ,JK)
      ! bet = b(k) - a(k)* gam(k)  
      PVARP(JIJ,JK)= ( ZY(JIJ,JK) - ZA(JIJ,JK) * PVARP(JIJ,JK-1) ) / ZBET(JIJ)
      ! res(k) = (y(k) -a(k)*res(k-1))/ bet
  END DO
   !$mnh_end_do()
END DO

! special treatment for the last level

!$mnh_do_concurrent(JIJ=1:IIJT)
DO JIJ=1,IIJT
   ZGAM(JIJ,IKE) = ZC(JIJ,IKE-1) / ZBET(JIJ) 
   ! gam(k) = c(k-1) / bet
   ZBET(JIJ)     = ZB(JIJ,IKE) - ZA(JIJ,IKE) * ZGAM(JIJ,IKE)
   ! bet = b(k) - a(k)* gam(k)  
   PVARP(JIJ,IKE)= ( ZY(JIJ,IKE) - ZA(JIJ,IKE) * PVARP(JIJ,IKE-1) ) / ZBET(JIJ)
   ! res(k) = (y(k) -a(k)*res(k-1))/ bet
END DO
!$mnh_end_do()  

!
!*       3.3 going down
!            ----------
!


DO JK = IKE-1,IKB,-1
#ifdef MNH_COMPILER_CCE
   
#endif
   !$mnh_do_concurrent(JIJ=1:IIJT)
  DO JIJ=1,IIJT
      PVARP(JIJ,JK) = PVARP(JIJ,JK) - ZGAM(JIJ,JK+1) * PVARP(JIJ,JK+1)
  END DO
   !$mnh_end_do()   
END DO

!
!
!*       4.  FILL THE UPPER AND LOWER EXTERNAL VALUES
!            ----------------------------------------
!

!$mnh_do_concurrent(JIJ=1:IIJT)
DO JIJ=1,IIJT
   PVARP(JIJ,IKB-1)=PVARP(JIJ,IKB)
   PVARP(JIJ,IKE+1)=0.
END DO
!$mnh_end_do()



!-------------------------------------------------------------------------------
!
END SUBROUTINE TRIDIAG_W
!
END MODULE MODE_TRIDIAG_W
