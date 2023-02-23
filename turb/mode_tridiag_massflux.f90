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
INTEGER                              :: JK, JIJ            ! loop counter
INTEGER :: IIJB,IIJE ! physical horizontal domain indices
INTEGER :: IKTB,IKTE
INTEGER :: IKT,IKB,IKA,IKU,IKE
INTEGER :: IKL
!
! ---------------------------------------------------------------------------
!                                              
!*      1.  Preliminaries
!           -------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('TRIDIAG_MASSFLUX',0,ZHOOK_HANDLE)
!
IIJE=D%NIJE
IIJB=D%NIJB
IKT=D%NKT
IKB=D%NKB
IKL=D%NKL
IKA=D%NKA
IKU=D%NKU
IKE=D%NKE
IKTB=D%NKTB
IKTE=D%NKTE
!
CALL MZM_MF(D, PRHODJ, ZMZM_RHODJ)
!$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
ZRHODJ_DFDT_O_DZ(:,:) = ZMZM_RHODJ(:,:)*PDFDT(:,:)/PDZZ(:,:)
!$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
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
!$mnh_expand_array(JIJ=IIJB:IIJE)
ZY(:,IKB) = PRHODJ(:,IKB)*PVARM(:,IKB)/PTSTEP             &
    - ZMZM_RHODJ(:,IKB+IKL) * PF(:,IKB+IKL)/PDZZ(:,IKB+IKL)     &
    + ZMZM_RHODJ(:,IKB  ) * PF(:,IKB  )/PDZZ(:,IKB  )     &
    + ZRHODJ_DFDT_O_DZ(:,IKB+IKL) * 0.5*PIMPL * PVARM(:,IKB+IKL)    &
    + ZRHODJ_DFDT_O_DZ(:,IKB+IKL) * 0.5*PIMPL * PVARM(:,IKB  )
!$mnh_end_expand_array(JIJ=IIJB:IIJE)
!
DO JK=1+IKTB,IKTE-1
  !$mnh_expand_array(JIJ=IIJB:IIJE)
  ZY(:,JK) = PRHODJ(:,JK)*PVARM(:,JK)/PTSTEP          &
    - ZMZM_RHODJ(:,JK+IKL) * PF(:,JK+IKL)/PDZZ(:,JK+IKL)    &
    + ZMZM_RHODJ(:,JK  ) * PF(:,JK  )/PDZZ(:,JK  )    &
    + ZRHODJ_DFDT_O_DZ(:,JK+IKL) * 0.5*PIMPL * PVARM(:,JK+IKL)  &
    + ZRHODJ_DFDT_O_DZ(:,JK+IKL) * 0.5*PIMPL * PVARM(:,JK  )  &
    - ZRHODJ_DFDT_O_DZ(:,JK  ) * 0.5*PIMPL * PVARM(:,JK  )  &
    - ZRHODJ_DFDT_O_DZ(:,JK  ) * 0.5*PIMPL * PVARM(:,JK-IKL)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE)
END DO
! 
IF (IKE==IKU) THEN
  !$mnh_expand_array(JIJ=IIJB:IIJE)
  ZY(:,IKE) = PRHODJ(:,IKE)*PVARM(:,IKE)/PTSTEP
  !$mnh_end_expand_array(JIJ=IIJB:IIJE)
ELSE
  !$mnh_expand_array(JIJ=IIJB:IIJE)
  ZY(:,IKE) = PRHODJ(:,IKE)*PVARM(:,IKE)/PTSTEP &
   - ZMZM_RHODJ(:,IKE+IKL) * PF(:,IKE+IKL)/PDZZ(:,IKE+IKL) &
   + ZMZM_RHODJ(:,IKE  ) * PF(:,IKE  )/PDZZ(:,IKE  ) &
   - ZRHODJ_DFDT_O_DZ(:,IKE ) * 0.5*PIMPL * PVARM(:,IKE  ) &
   - ZRHODJ_DFDT_O_DZ(:,IKE ) * 0.5*PIMPL * PVARM(:,IKE-IKL)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE)
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
  !$mnh_expand_array(JIJ=IIJB:IIJE)
  ZB(:,IKB) =   PRHODJ(:,IKB)/PTSTEP                   &
                + ZRHODJ_DFDT_O_DZ(:,IKB+IKL) * 0.5*PIMPL
  ZC(:,IKB) =   ZRHODJ_DFDT_O_DZ(:,IKB+IKL) * 0.5*PIMPL
  !$mnh_end_expand_array(JIJ=IIJB:IIJE)

  DO JK=1+IKTB,IKTE-1
    !$mnh_expand_array(JIJ=IIJB:IIJE)
    ZA(:,JK) = - ZRHODJ_DFDT_O_DZ(:,JK  ) * 0.5*PIMPL
    ZB(:,JK) =   PRHODJ(:,JK)/PTSTEP                   &
                 + ZRHODJ_DFDT_O_DZ(:,JK+IKL) * 0.5*PIMPL &
                 - ZRHODJ_DFDT_O_DZ(:,JK  ) * 0.5*PIMPL
    ZC(:,JK) =   ZRHODJ_DFDT_O_DZ(:,JK+IKL) * 0.5*PIMPL
    !$mnh_end_expand_array(JIJ=IIJB:IIJE)
  END DO
  !$mnh_expand_array(JIJ=IIJB:IIJE)
  ZA(:,IKE) = - ZRHODJ_DFDT_O_DZ(:,IKE  ) * 0.5*PIMPL
  ZB(:,IKE) =   PRHODJ(:,IKE)/PTSTEP                   &
                - ZRHODJ_DFDT_O_DZ(:,IKE  ) * 0.5*PIMPL
  !$mnh_end_expand_array(JIJ=IIJB:IIJE)
!
!*       3.2 going up
!            --------
!
  !$mnh_expand_array(JIJ=IIJB:IIJE)
  ZBET(:) = ZB(:,IKB)  ! bet = b(IKB)
  PVARP(:,IKB) = ZY(:,IKB) / ZBET(:)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE)

  !
  DO JK = IKB+IKL,IKE-IKL,IKL
    !$mnh_expand_array(JIJ=IIJB:IIJE)
    ZGAM(:,JK) = ZC(:,JK-IKL) / ZBET(:)
                                                    ! gam(k) = c(k-1) / bet
    ZBET(:)    = ZB(:,JK) - ZA(:,JK) * ZGAM(:,JK)
                                                    ! bet = b(k) - a(k)* gam(k)  
    PVARP(:,JK)= ( ZY(:,JK) - ZA(:,JK) * PVARP(:,JK-IKL) ) / ZBET(:)
                                        ! res(k) = (y(k) -a(k)*res(k-1))/ bet 
    !$mnh_end_expand_array(JIJ=IIJB:IIJE)
  END DO 
  !$mnh_expand_array(JIJ=IIJB:IIJE)
  ! special treatment for the last level
  ZGAM(:,IKE) = ZC(:,IKE-IKL) / ZBET(:)
                                                    ! gam(k) = c(k-1) / bet
  ZBET(:)     = ZB(:,IKE) - ZA(:,IKE) * ZGAM(:,IKE)
                                                    ! bet = b(k) - a(k)* gam(k)  
  PVARP(:,IKE)= ( ZY(:,IKE) - ZA(:,IKE) * PVARP(:,IKE-IKL) ) / &
                              &ZBET(:)
                                       ! res(k) = (y(k) -a(k)*res(k-1))/ bet 
  !$mnh_end_expand_array(JIJ=IIJB:IIJE)
!
!*       3.3 going down
!            ----------
!
  DO JK = IKE-IKL,IKB,-IKL
    !$mnh_expand_array(JIJ=IIJB:IIJE)
    PVARP(:,JK) = PVARP(:,JK) - ZGAM(:,JK+IKL) * PVARP(:,JK+IKL)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE)
  END DO
!
!
ELSE
  !!! EXPLICIT FORMULATION
  !
  DO JK=IKTB,IKTE
    !$mnh_expand_array(JIJ=IIJB:IIJE)
    PVARP(:,JK) = ZY(:,JK) * PTSTEP / PRHODJ(:,JK)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE)
  ENDDO
  !
END IF 
!
!
!*       4.  FILL THE UPPER AND LOWER EXTERNAL VALUES
!            ----------------------------------------
!
!$mnh_expand_array(JIJ=IIJB:IIJE)
PVARP(:,IKA)=PVARP(:,IKB)
PVARP(:,IKU)=PVARP(:,IKE)
!$mnh_end_expand_array(JIJ=IIJB:IIJE)
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('TRIDIAG_MASSFLUX',1,ZHOOK_HANDLE)
END SUBROUTINE TRIDIAG_MASSFLUX
END MODULE MODE_TRIDIAG_MASSFLUX
