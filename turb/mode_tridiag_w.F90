!MNH_LIC Copyright 2011-2022 CNRS, Meteo-France and Universite Paul Sabatier
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
USE MODD_PARAMETERS, ONLY : JPVEXT
USE MODD_DIMPHYEX,   ONLY: DIMPHYEX_t
USE MODE_SHUMAN_PHY, ONLY: MZM_PHY

#ifdef MNH_COMPILER_CCE
!$mnh_undef(LOOP)
!$mnh_undef(OPENACC)
#endif
!
#ifdef MNH_OPENACC
USE MODE_MNH_ZWORK,      ONLY: MNH_MEM_GET, MNH_MEM_POSITION_PIN, MNH_MEM_RELEASE
#endif
!
IMPLICIT NONE
!
!
!*       0.1 declarations of arguments
!
TYPE(DIMPHYEX_t),         INTENT(IN)    :: D
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PVARM   ! variable at t-1      at flux point
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PF      ! flux in dT/dt=-dF/dz at mass point
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PDFDDWDZ! dF/d(dW/dz)          at mass point
REAL,                   INTENT(IN) :: PTSTEP  ! Double time step
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PMZF_DZZ! Dz                   at mass point
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PRHODJ  ! (dry rho)*J          at mass point
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT):: PVARP   ! variable at t+1      at flux point
!
!*       0.2 declarations of local variables
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: ZRHODJ_DFDDWDZ_O_DZ2
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: ZMZM_RHODJ
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: ZA, ZB, ZC
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: ZY ,ZGAM ! RHS of the equation, 3D work array
REAL, DIMENSION(D%NIT,D%NJT)       :: ZBET     ! 2D work array
!
INTEGER                              :: JK            ! loop counter
INTEGER                              :: IKB,IKE       ! inner vertical limits
!
INTEGER  :: JIU,JJU,JKU
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZMZM3D_WORK1
INTEGER  :: JI,JJ
! ---------------------------------------------------------------------------



JIU =  size( pvarm, 1 )
JJU =  size( pvarm, 2 )
JKU =  size( pvarm, 3 )



!
!*      1.  Preliminaries
!           -------------
!
IKB=1+JPVEXT
IKE=SIZE(PVARM,3)-JPVEXT 
!
CALL MZM_PHY(D, PRHODJ(:, :, :), ZMZM3D_WORK1)

DO JK=1, JKU
  DO JJ=1, JJU
    DO JI=1, JIU
      ZMZM_RHODJ(JI, JJ, JK)  = ZMZM3D_WORK1(JI, JJ, JK)
    END DO
  END DO
END DO

!

DO JK=1, JKU
     DO JJ=1, JJU
       DO JI=1, JIU
         ZRHODJ_DFDDWDZ_O_DZ2(JI, JJ, JK) = PRHODJ(JI, JJ, JK)*PDFDDWDZ(JI, JJ, JK)/PMZF_DZZ(JI, JJ, JK)**2
    END DO
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

!$mnh_do_concurrent(JI=1:JIU,JJ=1:JJU)
   ZY(JI,JJ,IKB) = ZMZM_RHODJ(JI,JJ,IKB)*PVARM(JI,JJ,IKB)/PTSTEP              &
        - PRHODJ(JI,JJ,IKB  ) * PF(JI,JJ,IKB  )/PMZF_DZZ(JI,JJ,IKB  )           &
        + PRHODJ(JI,JJ,IKB-1) * PF(JI,JJ,IKB-1)/PMZF_DZZ(JI,JJ,IKB-1)           &
        + ZRHODJ_DFDDWDZ_O_DZ2(JI,JJ,IKB) * PVARM(JI,JJ,IKB+1)&
        - ZRHODJ_DFDDWDZ_O_DZ2(JI,JJ,IKB) * PVARM(JI,JJ,IKB  )
!$mnh_end_do()

!

!$mnh_do_concurrent(JI=1:JIU,JJ=1:JJU,JK=IKB+1:IKE-1)
  ZY(JI,JJ,JK) = ZMZM_RHODJ(JI,JJ,JK)*PVARM(JI,JJ,JK)/PTSTEP               &
       - PRHODJ(JI,JJ,JK  ) * PF(JI,JJ,JK  )/PMZF_DZZ(JI,JJ,JK  )          &
       + PRHODJ(JI,JJ,JK-1) * PF(JI,JJ,JK-1)/PMZF_DZZ(JI,JJ,JK-1)              &
       + ZRHODJ_DFDDWDZ_O_DZ2(JI,JJ,JK  ) * PVARM(JI,JJ,JK+1)  &
       - ZRHODJ_DFDDWDZ_O_DZ2(JI,JJ,JK  ) * PVARM(JI,JJ,JK  )  &
       - ZRHODJ_DFDDWDZ_O_DZ2(JI,JJ,JK-1) * PVARM(JI,JJ,JK  )  &
       + ZRHODJ_DFDDWDZ_O_DZ2(JI,JJ,JK-1) * PVARM(JI,JJ,JK-1)
!$mnh_end_do()  

! 

!$mnh_do_concurrent(JI=1:JIU,JJ=1:JJU)
   ZY(JI,JJ,IKE) = ZMZM_RHODJ(JI,JJ,IKE)*PVARM(JI,JJ,IKE)/PTSTEP              &
        - PRHODJ(JI,JJ,IKE  ) * PF(JI,JJ,IKE  )/PMZF_DZZ(JI,JJ,IKE  )           &
        + PRHODJ(JI,JJ,IKE-1) * PF(JI,JJ,IKE-1)/PMZF_DZZ(JI,JJ,IKE-1)           &
        - ZRHODJ_DFDDWDZ_O_DZ2(JI,JJ,IKE  ) * PVARM(JI,JJ,IKE )  &
        - ZRHODJ_DFDDWDZ_O_DZ2(JI,JJ,IKE-1) * PVARM(JI,JJ,IKE  ) &
        + ZRHODJ_DFDDWDZ_O_DZ2(JI,JJ,IKE-1) * PVARM(JI,JJ,IKE-1)
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

!$mnh_do_concurrent(JI=1:JIU,JJ=1:JJU)
  ZB(JI,JJ,IKB) =   ZMZM_RHODJ(JI,JJ,IKB)/PTSTEP      &
       - ZRHODJ_DFDDWDZ_O_DZ2(JI,JJ,IKB)
!$mnh_end_do()  


!$mnh_do_concurrent(JI=1:JIU,JJ=1:JJU)
   ZC(JI,JJ,IKB) =   ZRHODJ_DFDDWDZ_O_DZ2(JI,JJ,IKB)
!$mnh_end_do()   



!$mnh_do_concurrent(JI=1:JIU,JJ=1:JJU,JK=IKB+1:IKE-1)
    ZA(JI,JJ,JK) =   ZRHODJ_DFDDWDZ_O_DZ2(JI,JJ,JK-1)
    ZB(JI,JJ,JK) =   ZMZM_RHODJ(JI,JJ,JK)/PTSTEP      &
                 - ZRHODJ_DFDDWDZ_O_DZ2(JI,JJ,JK  ) &
                 - ZRHODJ_DFDDWDZ_O_DZ2(JI,JJ,JK-1)
    ZC(JI,JJ,JK) =   ZRHODJ_DFDDWDZ_O_DZ2(JI,JJ,JK  )
!$mnh_end_do()    



!$mnh_do_concurrent(JI=1:JIU,JJ=1:JJU)
   ZA(JI,JJ,IKE) =   ZRHODJ_DFDDWDZ_O_DZ2(JI,JJ,IKE-1)
!$mnh_end_do()   


!$mnh_do_concurrent(JI=1:JIU,JJ=1:JJU)
  ZB(JI,JJ,IKE) =   ZMZM_RHODJ(JI,JJ,IKE)/PTSTEP      &
                - ZRHODJ_DFDDWDZ_O_DZ2(JI,JJ,IKE  ) &
                - ZRHODJ_DFDDWDZ_O_DZ2(JI,JJ,IKE-1)
!$mnh_end_do()  

!
!
! acc wait
!
!*       3.2 going up
!            --------
!

!$mnh_do_concurrent(JI=1:JIU,JJ=1:JJU)
  ZBET(JI,JJ) = ZB(JI,JJ,IKB)  ! bet = b(ikb)
  PVARP(JI,JJ,IKB) = ZY(JI,JJ,IKB) / ZBET(JI,JJ)
!$mnh_end_do()

!


DO JK = IKB+1,IKE-1
#ifdef MNH_COMPILER_CCE
   
#endif
   !$mnh_do_concurrent(JI=1:JIU,JJ=1:JJU)
      ZGAM(JI,JJ,JK) = ZC(JI,JJ,JK-1) / ZBET(JI,JJ)  
      ! gam(k) = c(k-1) / bet
      ZBET(JI,JJ)    = ZB(JI,JJ,JK) - ZA(JI,JJ,JK) * ZGAM(JI,JJ,JK)
      ! bet = b(k) - a(k)* gam(k)  
      PVARP(JI,JJ,JK)= ( ZY(JI,JJ,JK) - ZA(JI,JJ,JK) * PVARP(JI,JJ,JK-1) ) / ZBET(JI,JJ)
      ! res(k) = (y(k) -a(k)*res(k-1))/ bet
   !$mnh_end_do()
END DO

! special treatment for the last level

!$mnh_do_concurrent(JI=1:JIU,JJ=1:JJU)
   ZGAM(JI,JJ,IKE) = ZC(JI,JJ,IKE-1) / ZBET(JI,JJ) 
   ! gam(k) = c(k-1) / bet
   ZBET(JI,JJ)     = ZB(JI,JJ,IKE) - ZA(JI,JJ,IKE) * ZGAM(JI,JJ,IKE)
   ! bet = b(k) - a(k)* gam(k)  
   PVARP(JI,JJ,IKE)= ( ZY(JI,JJ,IKE) - ZA(JI,JJ,IKE) * PVARP(JI,JJ,IKE-1) ) / ZBET(JI,JJ)
   ! res(k) = (y(k) -a(k)*res(k-1))/ bet
!$mnh_end_do()  

!
!*       3.3 going down
!            ----------
!


DO JK = IKE-1,IKB,-1
#ifdef MNH_COMPILER_CCE
   
#endif
   !$mnh_do_concurrent(JI=1:JIU,JJ=1:JJU)
      PVARP(JI,JJ,JK) = PVARP(JI,JJ,JK) - ZGAM(JI,JJ,JK+1) * PVARP(JI,JJ,JK+1)
   !$mnh_end_do()   
END DO

!
!
!*       4.  FILL THE UPPER AND LOWER EXTERNAL VALUES
!            ----------------------------------------
!

!$mnh_do_concurrent(JI=1:JIU,JJ=1:JJU)
   PVARP(JI,JJ,IKB-1)=PVARP(JI,JJ,IKB)
   PVARP(JI,JJ,IKE+1)=0.
!$mnh_end_do()



!-------------------------------------------------------------------------------
!
END SUBROUTINE TRIDIAG_W
!
END MODULE MODE_TRIDIAG_W
