!MNH_LIC Copyright 2019-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!
!       #######################
MODULE MODI_IBM_FORCING
  !       #######################
  !
  INTERFACE
     !
     SUBROUTINE IBM_FORCING(PRUS,PRVS,PRWS,PTHS,PRRS,PSVS,PTKS)
       !
       REAL, DIMENSION(:,:,:)  ,INTENT(INOUT)           :: PRUS,PRVS,PRWS
       REAL, DIMENSION(:,:,:)  ,INTENT(INOUT)           :: PTHS
       REAL, DIMENSION(:,:,:,:),INTENT(INOUT), OPTIONAL :: PRRS
       REAL, DIMENSION(:,:,:,:),INTENT(INOUT), OPTIONAL :: PSVS
       REAL, DIMENSION(:,:,:)  ,INTENT(INOUT), OPTIONAL :: PTKS
       !
     END SUBROUTINE IBM_FORCING
     !
  END INTERFACE
  !
END MODULE MODI_IBM_FORCING
!
!       ##########################################################
SUBROUTINE IBM_FORCING(PRUS,PRVS,PRWS,PTHS,PRRS,PSVS,PTKS)
  !       ##########################################################
  !
  !!****     *IBM_FORCING*  - routine to force all desired fields 
  !!
  !!      PURPOSE
  !!      -------
  !         The purpose of this routine is to compute variables in the virtual
  !         embedded solid region in regard of variables computed in the real
  !         fluid region
  !
  !!      METHOD
  !!      ------
  !!
  !!      EXTERNAL
  !!      --------
  !!        NONE
  !!
  !!      IMPLICIT ARGUMENTS
  !!      ------------------
  !!
  !!      REFERENCE
  !!      ---------
  !!
  !!      AUTHOR
  !!      ------
  !!        Franck Auguste       * CERFACS(AE) *
  !!
  !!      MODIFICATIONS
  !!      -------------
  !!        Original          01/01/2019
  !!
  !-----------------------------------------------------------------------------
  !       
  !**** 0. DECLARATIONS
  !     ---------------
  !
  ! module
  USE MODE_POS
  USE MODE_ll
  USE MODE_IO
  USE MODD_ARGSLIST_ll, ONLY : LIST_ll
  !
  ! declaration
  USE MODD_CST
  USE MODD_FIELD_n
  USE MODD_REF
  USE MODD_REF_n, ONLY: XRHODJ,XRHODREF,XTHVREF,XEXNREF,XRVREF
  USE MODD_PARAMETERS, ONLY: JPVEXT,JPHEXT
  USE MODD_IBM_PARAM_n
  USE MODD_LBC_n
  USE MODD_CONF
  USE MODD_CONF_n
  USE MODD_NSV
  USE MODD_TURB_n, ONLY: XTKEMIN
  USE MODD_PARAM_n
  USE MODD_DYN_n, ONLY: XTSTEP
  !
  ! interface
  USE MODI_IBM_AFFECTV
  USE MODI_IBM_AFFECTP
  USE MODI_SHUMAN
  !
  IMPLICIT NONE
  !
  !-----------------------------------------------------------------------------
  !
  !       0.1  declarations of arguments
  !
  REAL, DIMENSION(:,:,:)  ,INTENT(INOUT) :: PRUS,PRVS,PRWS
  REAL, DIMENSION(:,:,:)  ,INTENT(INOUT) :: PTHS
  REAL, DIMENSION(:,:,:,:),INTENT(INOUT), OPTIONAL :: PRRS
  REAL, DIMENSION(:,:,:,:),INTENT(INOUT), OPTIONAL :: PSVS
  REAL, DIMENSION(:,:,:)  ,INTENT(INOUT), OPTIONAL :: PTKS
  !
  !-----------------------------------------------------------------------------
  !
  !       0.2  declaration of local variables
  REAL, DIMENSION(:,:,:)  , ALLOCATABLE :: ZTMP,ZXMU,ZDIV,ZTKE
  REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZTMU,ZTRY
  INTEGER                               :: IIU,IJU,IKU,IKB,IKE
  INTEGER                               :: JRR,JSV
  TYPE(LIST_ll), POINTER                :: TZFIELDS_ll   
  INTEGER                               :: IINFO_ll
  !
  !-----------------------------------------------------------------------------
  !
  !**** 0. ALLOCATIONS
  !     --------------
  !
  IIU = SIZE(PRUS,1)
  IJU = SIZE(PRVS,2)
  IKU = SIZE(PRWS,3)
  !
  ALLOCATE(ZTMU(IIU,IJU,IKU,3),ZTMP(IIU,IJU,IKU),ZTRY(IIU,IJU,IKU,3), &
       ZXMU(IIU,IJU,IKU),ZDIV(IIU,IJU,IKU),ZTKE(IIU,IJU,IKU))
  !
  ZTMU=0.
  ZXMU=0.
  ZDIV=0.
  ZTMP=0.
  ZTRY=0.
  !
  IKB =   1 + JPVEXT
  IKE = IKU - JPVEXT
  !
  !-----------------------------------------------------------------------------
  !       
  !**** 1. PRELIMINARIES
  !     ----------------
  IF (NSV>=1) THEN
     !
     DO JSV=1,NSV
        WHERE (XIBM_LS(:,:,:,1).GT.XIBM_EPSI) PSVS(:,:,:,JSV) = XIBM_EPSI**1.5
     ENDDO
     !
  ENDIF
  !
  WHERE (XIBM_LS(:,:,:,1).GT.XIBM_EPSI) PTHS(:,:,:) = XTHVREF(:,:,:)
  !
  IF (NRR>=1) THEN
     WHERE (XIBM_LS(:,:,:,1).GT.XIBM_EPSI)  
        PRRS(:,:,:,1) = XRVREF(:,:,:)
        PTHS(:,:,:) = XTHVREF(:,:,:)/(1.+XRD/XRV*XRVREF(:,:,:))
     ENDWHERE
  ENDIF
  IF (NRR>=2) THEN
     DO JRR=2,NRR
        WHERE (XIBM_LS(:,:,:,1).GT.XIBM_EPSI) PRRS(:,:,:,JRR) = XIBM_EPSI
     ENDDO
  ENDIF
  !
  WHERE (XIBM_LS(:,:,:,2).GT.XIBM_EPSI) PRUS(:,:,:) = XIBM_EPSI
  WHERE (XIBM_LS(:,:,:,3).GT.XIBM_EPSI) PRVS(:,:,:) = XIBM_EPSI
  WHERE (XIBM_LS(:,:,:,4).GT.XIBM_EPSI) PRWS(:,:,:) = XIBM_EPSI
  IF (CTURB/='NONE') WHERE (XIBM_LS(:,:,:,1).GT.XIBM_EPSI) PTKS(:,:,:) = XTKEMIN
  !
  !**** 2. EXECUTIONS
  !     -------------
  !
  ! ======================
  ! === SCALAR FORCING ===
  ! ======================
  !
  IF (CTURB/='NONE') THEN
     ZTMP(:,:,:) = PTKS(:,:,:)
     ZTMP(:,:,IKB-1)=ZTMP(:,:,IKB)
     ZTMP(:,:,IKE+1)=ZTMP(:,:,IKE)
     ZXMU(:,:,:) = XIBM_XMUT(:,:,:)
     ZDIV(:,:,:) = XIBM_CURV(:,:,:)
     CALL IBM_AFFECTP(ZTMP,NIBM_LAYER_E,XIBM_RADIUS_E,XIBM_POWERS_E,&
          CIBM_MODE_INTE1_E,CIBM_MODE_INTE3_E,&
          CIBM_TYPE_BOUND_E,CIBM_MODE_BOUND_E,&
          CIBM_FORC_BOUND_E,XIBM_FORC_BOUND_E,ZXMU,ZDIV)
     ZTMP(:,:,IKB-1)=ZTMP(:,:,IKB)
     ZTMP(:,:,IKE+1)=XTKEMIN
     PTKS(:,:,:)=MAX(XTKEMIN,ZTMP(:,:,:))
  ENDIF
  !
  ZTMP(:,:,:) = PTHS(:,:,:)
  ZTMP(:,:,IKB-1)=ZTMP(:,:,IKB)
  ZTMP(:,:,IKE+1)=ZTMP(:,:,IKE)
  CALL IBM_AFFECTP(ZTMP,NIBM_LAYER_T,XIBM_RADIUS_T,XIBM_POWERS_T,&
       CIBM_MODE_INTE1_T,CIBM_MODE_INTE3_T,&
       CIBM_TYPE_BOUND_T,CIBM_MODE_BOUND_T,&
       CIBM_FORC_BOUND_T,XIBM_FORC_BOUND_T,ZXMU,ZDIV)
  ZTMP(:,:,:) =  ZTMP(:,:,:)
  ZTMP(:,:,IKB-1)=ZTMP(:,:,IKB)
  ZTMP(:,:,IKE+1)=ZTMP(:,:,IKE)
  PTHS(:,:,:) = MAX(ZTMP(:,:,:),250.)
  !
  IF (NRR>=1) THEN
     DO JRR=1,NRR
        ZTMP(:,:,:) = PRRS(:,:,:,JRR)
        ZTMP(:,:,IKB-1)=ZTMP(:,:,IKB)
        ZTMP(:,:,IKE+1)=ZTMP(:,:,IKE)
        CALL IBM_AFFECTP(ZTMP,NIBM_LAYER_R,XIBM_RADIUS_R,XIBM_POWERS_R,&
             CIBM_MODE_INTE1_R,CIBM_MODE_INTE3_R,&
             CIBM_TYPE_BOUND_R,CIBM_MODE_BOUND_R,&
             CIBM_FORC_BOUND_R,XIBM_FORC_BOUND_R,ZXMU,ZDIV)
        ZTMP(:,:,IKB-1)=ZTMP(:,:,IKB)
        ZTMP(:,:,IKE+1)=ZTMP(:,:,IKE)
        PRRS(:,:,:,JRR) = ZTMP(:,:,:)
     ENDDO
  ENDIF
  !
  IF (NSV>=1) THEN 
     DO JSV=1,NSV
        ZTMP(:,:,:) = PSVS(:,:,:,JSV)
        ZTMP(:,:,IKB-1)=ZTMP(:,:,IKB)
        ZTMP(:,:,IKE+1)=ZTMP(:,:,IKE)
        CALL IBM_AFFECTP(ZTMP,NIBM_LAYER_S,XIBM_RADIUS_S,XIBM_POWERS_S,&
             CIBM_MODE_INTE1_S,CIBM_MODE_INTE3_S,&
             CIBM_TYPE_BOUND_S,CIBM_MODE_BOUND_S,&
             CIBM_FORC_BOUND_S,XIBM_FORC_BOUND_S,ZXMU,ZDIV) 
        ZTMP(:,:,:) = MAX(XIBM_EPSI**1.5,ZTMP(:,:,:))
        ZTMP(:,:,IKB-1)=ZTMP(:,:,IKB)
        ZTMP(:,:,IKE+1)=ZTMP(:,:,IKE)
        PSVS(:,:,:,JSV) = ZTMP(:,:,:)
     ENDDO
  ENDIF
  !
  !=======================
  ! === VECTOR FORCING ===
  ! ======================
  !
  PRUS(:,:,IKB-1)=PRUS(:,:,IKB)
  PRUS(:,:,IKE+1)=PRUS(:,:,IKE)
  PRVS(:,:,IKB-1)=PRVS(:,:,IKB)
  PRVS(:,:,IKE+1)=PRVS(:,:,IKE)
  PRWS(:,:,IKB-1)=0.
  PRWS(:,:,IKE+1)=0.
  !
  ZTMU(:,:,:,1) = PRUS(:,:,:)
  ZTMU(:,:,:,2) = PRVS(:,:,:)
  ZTMU(:,:,:,3) = PRWS(:,:,:)
  !
  ZTMP(:,:,:) = PRUS(:,:,:)
  ZXMU(:,:,:) = MXM(XIBM_XMUT(:,:,:))
  ZDIV(:,:,:) = MXM(XIBM_CURV(:,:,:))
  CALL IBM_AFFECTV(ZTMP,ZTMU,ZTRY,'U',NIBM_LAYER_V,CIBM_MODE_INTE3_V,&
       CIBM_FORC_BOUNR_V,XIBM_RADIUS_V,XIBM_POWERS_V,&
       CIBM_MODE_INTE1NV,CIBM_TYPE_BOUNN_V,CIBM_MODE_BOUNN_V,CIBM_FORC_BOUNN_V ,XIBM_FORC_BOUNN_V,&
       CIBM_MODE_INTE1TV,CIBM_TYPE_BOUNT_V,CIBM_MODE_BOUNT_V,CIBM_FORC_BOUNT_V ,XIBM_FORC_BOUNT_V,&
       CIBM_MODE_INTE1CV,CIBM_TYPE_BOUNC_V,CIBM_MODE_BOUNC_V,CIBM_FORC_BOUNC_V ,XIBM_FORC_BOUNC_V,ZXMU,ZDIV)
  PRUS(:,:,:) = ZTMP(:,:,:)
  ZTMP(:,:,:) = PRVS(:,:,:)
  ZXMU(:,:,:) = MYM(XIBM_XMUT(:,:,:))
  ZDIV(:,:,:) = MYM(XIBM_CURV(:,:,:))
  CALL IBM_AFFECTV(ZTMP,ZTMU,ZTRY,'V',NIBM_LAYER_V,CIBM_MODE_INTE3_V,&
       CIBM_FORC_BOUNR_V,XIBM_RADIUS_V,XIBM_POWERS_V,&
       CIBM_MODE_INTE1NV,CIBM_TYPE_BOUNN_V,CIBM_MODE_BOUNN_V,CIBM_FORC_BOUNN_V ,XIBM_FORC_BOUNN_V,&
       CIBM_MODE_INTE1TV,CIBM_TYPE_BOUNT_V,CIBM_MODE_BOUNT_V,CIBM_FORC_BOUNT_V ,XIBM_FORC_BOUNT_V,&
       CIBM_MODE_INTE1CV,CIBM_TYPE_BOUNC_V,CIBM_MODE_BOUNC_V,CIBM_FORC_BOUNC_V ,XIBM_FORC_BOUNC_V,ZXMU,ZDIV)
  PRVS(:,:,:) = ZTMP(:,:,:)
  ZTMP(:,:,:) = PRWS(:,:,:)
  ZXMU(:,:,:) = MZM(XIBM_XMUT(:,:,:))
  ZDIV(:,:,:) = MZM(XIBM_CURV(:,:,:))
  CALL IBM_AFFECTV(ZTMP,ZTMU,ZTRY,'W',NIBM_LAYER_V,CIBM_MODE_INTE3_V,&
       CIBM_FORC_BOUNR_V,XIBM_RADIUS_V,XIBM_POWERS_V,&
       CIBM_MODE_INTE1NV,CIBM_TYPE_BOUNN_V,CIBM_MODE_BOUNN_V,CIBM_FORC_BOUNN_V ,XIBM_FORC_BOUNN_V,&
       CIBM_MODE_INTE1TV,CIBM_TYPE_BOUNT_V,CIBM_MODE_BOUNT_V,CIBM_FORC_BOUNT_V ,XIBM_FORC_BOUNT_V,&
       CIBM_MODE_INTE1CV,CIBM_TYPE_BOUNC_V,CIBM_MODE_BOUNC_V,CIBM_FORC_BOUNC_V ,XIBM_FORC_BOUNC_V,ZXMU,ZDIV)
  PRWS(:,:,:) = ZTMP(:,:,:)
  PRUS(:,:,IKB-1)=PRUS(:,:,IKB)
  PRUS(:,:,IKE+1)=PRUS(:,:,IKE)
  PRVS(:,:,IKB-1)=PRVS(:,:,IKB)
  PRVS(:,:,IKE+1)=PRVS(:,:,IKE)
  PRWS(:,:,IKB-1)=0.
  PRWS(:,:,IKB)  =0.
  PRWS(:,:,IKE+1)=0.
  !
  !**** 3. COMMUNICATIONS 
  !     -----------------
  !
  IF (.NOT. LIBM_TROUBLE) THEN
     !
     NULLIFY(TZFIELDS_ll)
     CALL ADD3DFIELD_ll(TZFIELDS_ll,PTHS(:,:,:),'IBM_FORCING::PTHS')
     IF (CTURB/='NONE') CALL ADD3DFIELD_ll(TZFIELDS_ll,PTKS(:,:,:),'IBM_FORCING::PTKS')
     CALL ADD3DFIELD_ll(TZFIELDS_ll,PRUS(:,:,:),'IBM_FORCING::PRUS')
     CALL ADD3DFIELD_ll(TZFIELDS_ll,PRVS(:,:,:),'IBM_FORCING::PRVS')
     CALL ADD3DFIELD_ll(TZFIELDS_ll,PRWS(:,:,:),'IBM_FORCING::PRWS')
     IF (NRR>=1) THEN
        DO JRR=1,NRR
           CALL ADD3DFIELD_ll(TZFIELDS_ll,PRRS(:,:,:,JRR),'IBM_FORCING::PRRS')
        ENDDO
     ENDIF
     IF (NSV>=1) THEN
        DO JSV=1,NSV
           CALL ADD3DFIELD_ll(TZFIELDS_ll,PSVS(:,:,:,JSV),'IBM_FORCING::PSVS')
        ENDDO
     ENDIF
     !
     CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
     CALL CLEANLIST_ll(TZFIELDS_ll)
     !
  ENDIF
  !
  !**** 4. DEALLOCATIONS
  !     ----------------
  !
  DEALLOCATE(ZTMP,ZTMU,ZTRY,ZXMU,ZDIV,ZTKE)
  !
  RETURN
  ! 
END SUBROUTINE IBM_FORCING
