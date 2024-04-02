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
!!       M.Moge          04/2016 Use openACC directives to port the TURB part of Meso-NH on GPU
!! ---------------------------------------------------------------------
!
!*       0. DECLARATIONS
!
USE MODD_PARAMETERS, ONLY : JPVEXT


#ifndef MNH_OPENACC
USE MODI_SHUMAN
#else
USE MODI_SHUMAN_DEVICE
#endif
#if defined(MNH_BITREP) || defined(MNH_BITREP_OMP)
USE MODI_BITREP
#endif
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
REAL, DIMENSION(:,:,:), INTENT(IN) :: PVARM   ! variable at t-1      at flux point
REAL, DIMENSION(:,:,:), INTENT(IN) :: PF      ! flux in dT/dt=-dF/dz at mass point
REAL, DIMENSION(:,:,:), INTENT(IN) :: PDFDDWDZ! dF/d(dW/dz)          at mass point
REAL,                   INTENT(IN) :: PTSTEP  ! Double time step
REAL, DIMENSION(:,:,:), INTENT(IN) :: PMZF_DZZ! Dz                   at mass point
REAL, DIMENSION(:,:,:), INTENT(IN) :: PRHODJ  ! (dry rho)*J          at mass point
!
REAL, DIMENSION(:,:,:), INTENT(OUT):: PVARP   ! variable at t+1      at flux point
!
!*       0.2 declarations of local variables
!
REAL, DIMENSION(:,:,:), pointer , contiguous :: ZRHODJ_DFDDWDZ_O_DZ2
REAL, DIMENSION(:,:,:), pointer , contiguous :: ZMZM_RHODJ
REAL, DIMENSION(:,:,:), pointer , contiguous :: ZA, ZB, ZC
REAL, DIMENSION(:,:,:), pointer , contiguous :: ZY ,ZGAM ! RHS of the equation, 3D work array
REAL, DIMENSION(:,:),   pointer , contiguous :: ZBET     ! 2D work array
!
INTEGER                              :: JK            ! loop counter
INTEGER                              :: IKB,IKE       ! inner vertical limits
!
INTEGER  :: JIU,JJU,JKU
INTEGER  :: JI,JJ
! ---------------------------------------------------------------------------

!$acc data present_crm( PVARM, PF, PDFDDWDZ, PMZF_DZZ, PRHODJ, PVARP )

JIU =  size( pvarm, 1 )
JJU =  size( pvarm, 2 )
JKU =  size( pvarm, 3 )

#ifndef MNH_OPENACC
allocate( zrhodj_dfddwdz_o_dz2(JIU,JJU,JKU ) )
allocate( zmzm_rhodj          (JIU,JJU,JKU ) )
allocate( za                  (JIU,JJU,JKU ) )
allocate( zb                  (JIU,JJU,JKU ) )
allocate( zc                  (JIU,JJU,JKU ) )
allocate( zy                  (JIU,JJU,JKU ) )
allocate( zgam                (JIU,JJU,JKU ) )
allocate( zbet                (JIU,JJU ) )
#else
!Pin positions in the pools of MNH memory
CALL MNH_MEM_POSITION_PIN( 'TRIDIAG_W' )

CALL MNH_MEM_GET( zrhodj_dfddwdz_o_dz2,JIU,JJU,JKU )
CALL MNH_MEM_GET( zmzm_rhodj          ,JIU,JJU,JKU )
CALL MNH_MEM_GET( za                  ,JIU,JJU,JKU )
CALL MNH_MEM_GET( zb                  ,JIU,JJU,JKU )
CALL MNH_MEM_GET( zc                  ,JIU,JJU,JKU )
CALL MNH_MEM_GET( zy                  ,JIU,JJU,JKU )
CALL MNH_MEM_GET( zgam                ,JIU,JJU,JKU )
CALL MNH_MEM_GET( zbet                ,JIU,JJU )
#endif 

!$acc data present( ZRHODJ_DFDDWDZ_O_DZ2, ZMZM_RHODJ, ZA, ZB, ZC, ZY, ZGAM, ZBET )

!
!*      1.  Preliminaries
!           -------------
!
IKB=1+JPVEXT
IKE=SIZE(PVARM,3)-JPVEXT 
!
#ifndef MNH_OPENACC
ZMZM_RHODJ  = MZM(PRHODJ)
#else
CALL MZM_DEVICE(PRHODJ,ZMZM_RHODJ)
#endif
!$acc kernels ! async 
#if !defined(MNH_BITREP) && !defined(MNH_BITREP_OMP)
ZRHODJ_DFDDWDZ_O_DZ2 = PRHODJ*PDFDDWDZ/PMZF_DZZ**2
#else
!$mnh_expand_array(JI=1:JIU,JJ=1:JJU,JK=1:JKU)
   ZRHODJ_DFDDWDZ_O_DZ2(:,:,:) = PRHODJ(:,:,:)*PDFDDWDZ(:,:,:)/BR_P2(PMZF_DZZ(:,:,:))
!$mnh_end_expand_array()
#endif
!$acc end kernels
!
!$acc kernels ! async
ZA=0.
ZB=0.
ZC=0.
ZY=0.
!$acc end kernels
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
!$acc kernels ! async
!$mnh_do_concurrent(JI=1:JIU,JJ=1:JJU)
   ZY(JI,JJ,IKB) = ZMZM_RHODJ(JI,JJ,IKB)*PVARM(JI,JJ,IKB)/PTSTEP              &
        - PRHODJ(JI,JJ,IKB  ) * PF(JI,JJ,IKB  )/PMZF_DZZ(JI,JJ,IKB  )           &
        + PRHODJ(JI,JJ,IKB-1) * PF(JI,JJ,IKB-1)/PMZF_DZZ(JI,JJ,IKB-1)           &
        + ZRHODJ_DFDDWDZ_O_DZ2(JI,JJ,IKB) * PVARM(JI,JJ,IKB+1)&
        - ZRHODJ_DFDDWDZ_O_DZ2(JI,JJ,IKB) * PVARM(JI,JJ,IKB  )
!$mnh_end_do()
!$acc end kernels
!
!$acc kernels ! async
!$mnh_do_concurrent(JI=1:JIU,JJ=1:JJU,JK=IKB+1:IKE-1)
  ZY(JI,JJ,JK) = ZMZM_RHODJ(JI,JJ,JK)*PVARM(JI,JJ,JK)/PTSTEP               &
       - PRHODJ(JI,JJ,JK  ) * PF(JI,JJ,JK  )/PMZF_DZZ(JI,JJ,JK  )          &
       + PRHODJ(JI,JJ,JK-1) * PF(JI,JJ,JK-1)/PMZF_DZZ(JI,JJ,JK-1)              &
       + ZRHODJ_DFDDWDZ_O_DZ2(JI,JJ,JK  ) * PVARM(JI,JJ,JK+1)  &
       - ZRHODJ_DFDDWDZ_O_DZ2(JI,JJ,JK  ) * PVARM(JI,JJ,JK  )  &
       - ZRHODJ_DFDDWDZ_O_DZ2(JI,JJ,JK-1) * PVARM(JI,JJ,JK  )  &
       + ZRHODJ_DFDDWDZ_O_DZ2(JI,JJ,JK-1) * PVARM(JI,JJ,JK-1)
!$mnh_end_do()  
!$acc end kernels
! 
!$acc kernels ! async
!$mnh_do_concurrent(JI=1:JIU,JJ=1:JJU)
   ZY(JI,JJ,IKE) = ZMZM_RHODJ(JI,JJ,IKE)*PVARM(JI,JJ,IKE)/PTSTEP              &
        - PRHODJ(JI,JJ,IKE  ) * PF(JI,JJ,IKE  )/PMZF_DZZ(JI,JJ,IKE  )           &
        + PRHODJ(JI,JJ,IKE-1) * PF(JI,JJ,IKE-1)/PMZF_DZZ(JI,JJ,IKE-1)           &
        - ZRHODJ_DFDDWDZ_O_DZ2(JI,JJ,IKE  ) * PVARM(JI,JJ,IKE )  &
        - ZRHODJ_DFDDWDZ_O_DZ2(JI,JJ,IKE-1) * PVARM(JI,JJ,IKE  ) &
        + ZRHODJ_DFDDWDZ_O_DZ2(JI,JJ,IKE-1) * PVARM(JI,JJ,IKE-1)
!$mnh_end_do()
!$acc end kernels
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
!$acc kernels ! async
!$mnh_do_concurrent(JI=1:JIU,JJ=1:JJU)
  ZB(JI,JJ,IKB) =   ZMZM_RHODJ(JI,JJ,IKB)/PTSTEP      &
       - ZRHODJ_DFDDWDZ_O_DZ2(JI,JJ,IKB)
!$mnh_end_do()  
!$acc end kernels
!$acc kernels ! async
!$mnh_do_concurrent(JI=1:JIU,JJ=1:JJU)
   ZC(JI,JJ,IKB) =   ZRHODJ_DFDDWDZ_O_DZ2(JI,JJ,IKB)
!$mnh_end_do()   
!$acc end kernels

!$acc kernels ! async
!$mnh_do_concurrent(JI=1:JIU,JJ=1:JJU,JK=IKB+1:IKE-1)
    ZA(JI,JJ,JK) =   ZRHODJ_DFDDWDZ_O_DZ2(JI,JJ,JK-1)
    ZB(JI,JJ,JK) =   ZMZM_RHODJ(JI,JJ,JK)/PTSTEP      &
                 - ZRHODJ_DFDDWDZ_O_DZ2(JI,JJ,JK  ) &
                 - ZRHODJ_DFDDWDZ_O_DZ2(JI,JJ,JK-1)
    ZC(JI,JJ,JK) =   ZRHODJ_DFDDWDZ_O_DZ2(JI,JJ,JK  )
!$mnh_end_do()    
!$acc end kernels

!$acc kernels ! async
!$mnh_do_concurrent(JI=1:JIU,JJ=1:JJU)
   ZA(JI,JJ,IKE) =   ZRHODJ_DFDDWDZ_O_DZ2(JI,JJ,IKE-1)
!$mnh_end_do()   
!$acc end kernels
!$acc kernels ! async
!$mnh_do_concurrent(JI=1:JIU,JJ=1:JJU)
  ZB(JI,JJ,IKE) =   ZMZM_RHODJ(JI,JJ,IKE)/PTSTEP      &
                - ZRHODJ_DFDDWDZ_O_DZ2(JI,JJ,IKE  ) &
                - ZRHODJ_DFDDWDZ_O_DZ2(JI,JJ,IKE-1)
!$mnh_end_do()  
!$acc end kernels
!
!
! acc wait
!
!*       3.2 going up
!            --------
!
!$acc kernels
!$mnh_do_concurrent(JI=1:JIU,JJ=1:JJU)
  ZBET(JI,JJ) = ZB(JI,JJ,IKB)  ! bet = b(ikb)
  PVARP(JI,JJ,IKB) = ZY(JI,JJ,IKB) / ZBET(JI,JJ)
!$mnh_end_do()
!$acc end kernels
!
!$acc parallel 
!$acc loop seq
DO JK = IKB+1,IKE-1
#ifdef MNH_COMPILER_CCE
   !$acc loop independent
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
!$acc end parallel
! special treatment for the last level
!$acc kernels
!$mnh_do_concurrent(JI=1:JIU,JJ=1:JJU)
   ZGAM(JI,JJ,IKE) = ZC(JI,JJ,IKE-1) / ZBET(JI,JJ) 
   ! gam(k) = c(k-1) / bet
   ZBET(JI,JJ)     = ZB(JI,JJ,IKE) - ZA(JI,JJ,IKE) * ZGAM(JI,JJ,IKE)
   ! bet = b(k) - a(k)* gam(k)  
   PVARP(JI,JJ,IKE)= ( ZY(JI,JJ,IKE) - ZA(JI,JJ,IKE) * PVARP(JI,JJ,IKE-1) ) / ZBET(JI,JJ)
   ! res(k) = (y(k) -a(k)*res(k-1))/ bet
!$mnh_end_do()  
!$acc end kernels
!
!*       3.3 going down
!            ----------
!
!$acc parallel   
!$acc loop seq
DO JK = IKE-1,IKB,-1
#ifdef MNH_COMPILER_CCE
   !$acc loop independent
#endif
   !$mnh_do_concurrent(JI=1:JIU,JJ=1:JJU)
      PVARP(JI,JJ,JK) = PVARP(JI,JJ,JK) - ZGAM(JI,JJ,JK+1) * PVARP(JI,JJ,JK+1)
   !$mnh_end_do()   
END DO
!$acc end parallel
!
!
!*       4.  FILL THE UPPER AND LOWER EXTERNAL VALUES
!            ----------------------------------------
!
!$acc kernels
!$mnh_do_concurrent(JI=1:JIU,JJ=1:JJU)
   PVARP(JI,JJ,IKB-1)=PVARP(JI,JJ,IKB)
   PVARP(JI,JJ,IKE+1)=0.
!$mnh_end_do()
!$acc end kernels

!$acc end data

#ifndef MNH_OPENACC
deallocate (ZRHODJ_DFDDWDZ_O_DZ2, ZMZM_RHODJ, ZA, ZB, ZC, ZY, ZGAM, ZBET)
#else
!Release all memory allocated with MNH_MEM_GET calls since last call to MNH_MEM_POSITION_PIN
CALL MNH_MEM_RELEASE( 'TRIDIAG_W' )
#endif


!$acc end data

!-------------------------------------------------------------------------------
!
END SUBROUTINE TRIDIAG_W
!
END MODULE MODE_TRIDIAG_W
