!MNH_LIC Copyright 2019-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!
MODULE MODE_IBM_MIXINGLENGTH
IMPLICIT NONE
CONTAINS
SUBROUTINE IBM_MIXINGLENGTH(D,PLM,PLEPS,PMU,PHI,PTKE)
  !     ###################################################
  !
  !****  *IBM_MIXINGLENGTH* - Alteration of the mixing lenght (IBM)
  !
  !    PURPOSE
  !    -------
  !     The limitation is corrected for the immersed bonudary method:
  !        => using the level set phi 
  !        => LM < k(-phi)
  !
  !    METHOD
  !    ------
  !
  !    INDEX
  !    -----
  !
  !    IMPLICIT ARGUMENTS
  !    ------------------
  !
  !    REFERENCE
  !    ---------
  !
  !    AUTHOR
  !    ------
  !	
  !      Franck Auguste       * CERFACS(AE) *
  !
  !    MODIFICATIONS
  !    -------------
  !      Original    01/01/2019
  !
  !-------------------------------------------------------------------------------
  !
  !**** 0.    DECLARATIONS
  !     ------------------
  !
  ! module
  USE MODE_ll
  USE MODD_DIMPHYEX,   ONLY: DIMPHYEX_t
  !
  ! declaration
  USE MODD_FIELD_n
  USE MODD_PARAMETERS
  USE MODD_IBM_PARAM_n
  USE MODD_REF_n, ONLY: XRHODJ,XRHODREF
  USE MODD_CTURB
  USE MODD_CST
  USE MODD_GRID_n, ONLY: XZZ
  !
  ! interface
  !
  IMPLICIT NONE
  !
  !------------------------------------------------------------------------------
  !
  !       0.1   Declaration of arguments
  TYPE(DIMPHYEX_t),       INTENT(IN)    :: D
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(INOUT) :: PLM
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(INOUT) :: PLEPS
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT)   :: PMU
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)    :: PHI
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)    :: PTKE
  !
  !------------------------------------------------------------------------------
  !
  !       0.2   Declaration of local variables
  REAL,   DIMENSION(D%NIT,D%NJT,D%NKT) :: ZALPHA,ZBETA
  REAL,   DIMENSION(D%NIT,D%NJT,D%NKT) :: ZLM,ZMU,ZLN
  INTEGER                                                :: IKU,IKB,IKE,IIB,IIE,IJB,IJE
  REAL                                                   :: ZKARMAN
  !
  !-------------------------------------------------------------------------------
  !
  IKU = D%NKT
  IKE = D%NKE
  IKB = D%NKB
  IIB = D%NIB
  IIE = D%NIE
  IJB = D%NJB
  IJE = D%NJE
  !
  !  Turbulent velocity
  !
  ZMU(:,:,:)= 2.*XIBM_CNU**0.25*(PTKE(:,:,:))**(1./2.)  !2 correspond to KTKE
  ZKARMAN = XKARMAN*XCED/XIBM_CNU**0.75
  !
  ! Mesh scale
  !
  ZLN(:,:,:) = (XRHODJ(:,:,:)/XRHODREF(:,:,:))**(1./3.)
  ZLM(:,:,:) = PLM(:,:,:)
  !
  ! limit domain 
  !
  ZBETA(:,:,:)= XZZ (:,:,:)
  ZBETA(:,:,IKB:IKE) = 0.5*(XZZ(:,:,IKB+1:IKE+1)+XZZ(:,:,IKB:IKE))
  ZBETA(:,:,IKB-1) = -0.5*(XZZ(:,:,IKB+1)+XZZ(:,:,IKB))
  ZBETA(:,:,IKE+1) = ZBETA(:,:,IKE)
  ZLM(:,:,:) = MIN(ZLM(:,:,:),+ZKARMAN*ZBETA(:,:,:))
  !
  ! limit immersed wall
  !
  ZLM(:,:,:) = MIN(ZLM(:,:,:),-ZKARMAN*PHI(:,:,:))
  !
  ! limit physical scale
  ZALPHA(:,:,:) =  MIN(9.8*XIBM_RUG,0.5*ZKARMAN*ZLN(:,:,:))
  ZLM(:,:,:) = MAX(ZALPHA(:,:,:),ZLM(:,:,:))
  !
  !  Boundary condition
  ZMU(:,:,IKB-1)=  ZMU(:,:,IKB)
  ZLM(:,:,IKB-1)=  ZLM(:,:,IKB)
  ZMU(:,:,IKE+1)=  ZMU(:,:,IKE)
  ZLM(:,:,IKE+1)=  ZLM(:,:,IKE)
  IF (LEAST_ll()) THEN
     ZMU(IIE+1,:,:)=  ZMU(IIE,:,:)
     ZLM(IIE+1,:,:)=  ZLM(IIE,:,:)
  ENDIF
  IF (LWEST_ll())  THEN
     ZMU(IIB-1,:,:)=  ZMU(IIB,:,:)
     ZLM(IIB-1,:,:)=  ZLM(IIB,:,:)
  ENDIF
  IF (LNORTH_ll()) THEN
     ZMU(:,IJE+1,:)=  ZMU(:,IJE,:)
     ZLM(:,IJE+1,:)=  ZLM(:,IJE,:)
  ENDIF
  IF (LSOUTH_ll()) THEN
     ZMU(:,IJB-1,:)=  ZMU(:,IJB,:)
     ZLM(:,IJB-1,:)=  ZLM(:,IJB,:)
  ENDIF
  !
  !Communication
  PLM(:,:,:) = ZLM(:,:,:)
  PLEPS(:,:,:) = PLM(:,:,:) 
  PMU(:,:,:) = ZMU(:,:,:)
  !
END SUBROUTINE IBM_MIXINGLENGTH
END MODULE MODE_IBM_MIXINGLENGTH
