!MNH_LIC Copyright 2019-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!
!       ##########################
MODULE MODI_IBM_FORCING_TR
  !       ##########################
  !
  INTERFACE
     !
     SUBROUTINE IBM_FORCING_TR(PRUS,PRVS,PRWS,PTHS,PRRS,PSVS,PTKS)
       !
       REAL, DIMENSION(:,:,:)   ,INTENT(INOUT)          :: PRUS,PRVS,PRWS
       REAL, DIMENSION(:,:,:)   ,INTENT(INOUT)          :: PTHS
       REAL, DIMENSION(:,:,:,:) ,INTENT(INOUT),OPTIONAL :: PRRS
       REAL, DIMENSION(:,:,:,:) ,INTENT(INOUT),OPTIONAL :: PSVS
       REAL, DIMENSION(:,:,:)   ,INTENT(INOUT),OPTIONAL :: PTKS
       !
     END SUBROUTINE IBM_FORCING_TR
     !
  END INTERFACE
  !
END MODULE MODI_IBM_FORCING_TR
!
!
!       #############################################################
SUBROUTINE IBM_FORCING_TR(PRUS,PRVS,PRWS,PTHS,PRRS,PSVS,PTKS)
  !       #############################################################
  !
  !!****     *IBM_FORCING_TR*  - routine to force all desired fields 
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
  !------------------------------------------------------------------------------
  !       
  !**** 0. DECLARATIONS
  !     ---------------
  !
  ! module
  USE MODE_POS
  USE MODE_ll
  USE MODE_IO
  USE MODD_ARGSLIST_ll, ONLY: LIST_ll
  !
  ! declaration
  USE MODD_CST, ONLY: XRD,XRV
  USE MODD_REF_n, ONLY: XRHODJ,XRHODREF,XTHVREF,XRVREF
  USE MODD_PARAMETERS, ONLY: JPVEXT,JPHEXT
  USE MODD_IBM_PARAM_n
  USE MODD_LBC_n
  USE MODD_CONF
  USE MODD_CONF_n
  USE MODD_NSV
  USE MODD_TURB_n, ONLY: XTKEMIN
  USE MODD_PARAM_n
  !
  ! interface
  !
  IMPLICIT NONE
  !
  !-----------------------------------------------------------------------------
  !
  !       0.1  declarations of arguments
  !
  REAL, DIMENSION(:,:,:)  ,INTENT(INOUT)          :: PRUS,PRVS,PRWS
  REAL, DIMENSION(:,:,:)  ,INTENT(INOUT)          :: PTHS
  REAL, DIMENSION(:,:,:,:),INTENT(INOUT),OPTIONAL :: PRRS
  REAL, DIMENSION(:,:,:,:),INTENT(INOUT),OPTIONAL :: PSVS
  REAL, DIMENSION(:,:,:)  ,INTENT(INOUT),OPTIONAL :: PTKS
  !
  !-----------------------------------------------------------------------------
  !
  !       0.2  declaration of local variables
  INTEGER                         :: JI,JJ,JK,JI2,JJ2,JK2,IIU,IJU,IKU,JL
  INTEGER                         :: JIM1,JJM1,JKM1,JIP1,JJP1,JKP1
  INTEGER                         :: IIE,IIB,IJE,IJB,IKB,IKE
  REAL                            :: ZSUM1,ZSUM2,ZSUM4
  REAL, DIMENSION(:), ALLOCATABLE :: ZSUM3,ZSUM5
  TYPE(LIST_ll), POINTER          :: TZFIELDS_ll   
  INTEGER                         :: IINFO_ll
  !
  !-----------------------------------------------------------------------------
  !
  !**** 0. ALLOCATIONS
  !     --------------
  CALL GET_INDICE_ll(IIB,IJB,IIE,IJE)
  IIU = SIZE(PRUS,1)
  IJU = SIZE(PRUS,2)
  IKU = SIZE(PRUS,3)
  IKB =   1 + JPVEXT
  IKE = SIZE(PRUS,3) - JPVEXT
  !
  !-----------------------------------------------------------------------------
  ! 
  ! Problems in GCT ? => imposition of the adjacent value
  DO JI=IIB,IIE
     DO JJ=IJB,IJE
        DO JK=IKB,IKE
           !
           IF (XIBM_SUTR(JI,JJ,JK,1).LT.0.5) THEN
              !
              JIM1 = JI-1
              JJM1 = JJ-1
              JKM1 = JK-1
              JIP1 = JI+1
              JJP1 = JJ+1
              JKP1 = JK+1
              ZSUM1 = 0.
              ZSUM2 = 0.
              IF (NSV>=1) ALLOCATE(ZSUM3(NSV))
              ZSUM3 = 0.
              ZSUM4 = 0.
              IF (NRR>=1) ALLOCATE(ZSUM5(NRR))
              ZSUM5 = 0.
              !
              DO JI2=JIM1,JIP1
                 DO JJ2=JJM1,JJP1
                    DO JK2=JKM1,JKP1
                       !
                       ZSUM1 = ZSUM1 + (XIBM_SUTR(JI2,JJ2,JK2,1))
                       ZSUM2 = ZSUM2 + (XIBM_SUTR(JI2,JJ2,JK2,1))*PTHS(JI2,JJ2,JK2)
                       IF (NRR>=1) THEN
                          DO JL = 1,NRR
                             ZSUM5(JL) = ZSUM5(JL) + (XIBM_SUTR(JI2,JJ2,JK2,1))*PRRS(JI2,JJ2,JK2,JL)
                          ENDDO
                       ENDIF
                       IF (NSV>=1) THEN
                          DO JL = 1,NSV
                             ZSUM3(JL) = ZSUM3(JL) + (XIBM_SUTR(JI2,JJ2,JK2,1))*PSVS(JI2,JJ2,JK2,JL)
                          ENDDO
                       ENDIF
                       IF (CTURB/='NONE') ZSUM4 = ZSUM4 + (XIBM_SUTR(JI2,JJ2,JK2,1))*PTKS(JI2,JJ2,JK2)
                       !
                    ENDDO
                 ENDDO
              ENDDO
              !
              PTHS(JI,JJ,JK) = XTHVREF(JI,JJ,JK)
              IF (NRR>=1) PTHS(JI,JJ,JK) = XTHVREF(JI,JJ,JK)/(1.+XRD/XRV*XRVREF(JI,JJ,JK))
              IF (ZSUM1.GT.XIBM_EPSI)     PTHS(JI,JJ,JK) =  ZSUM2/ZSUM1
              IF (NRR>=1) THEN
                 PRRS(JI,JJ,JK,1) = XRVREF(JI,JJ,JK)
                 IF (ZSUM1.GT.XIBM_EPSI) PRRS(JI,JJ,JK,1) = ZSUM5(1)/ZSUM1
                 IF (NRR>=2) THEN
                    DO JL = 2,NRR
                       PRRS(JI,JJ,JK,JL) = 0.
                       IF (ZSUM1.GT.XIBM_EPSI) PRRS(JI,JJ,JK,JL) = ZSUM5(JL)/ZSUM1
                    ENDDO
                 ENDIF
              ENDIF
              !
              IF (NSV>=1) THEN
                 DO JL = 1,NSV
                    PSVS(JI,JJ,JK,JL) = 0.
                    IF (ZSUM1.GT.XIBM_EPSI) PSVS(JI,JJ,JK,JL) = ZSUM3(JL)/ZSUM1
                 ENDDO
              ENDIF
              !
              IF (CTURB/='NONE') PTKS(JI,JJ,JK) = XTKEMIN
              IF (ZSUM1.GT.XIBM_EPSI.AND.(CTURB/='NONE'))   PTKS(JI,JJ,JK) = ZSUM4/ZSUM1
              IF (NSV>=1) DEALLOCATE(ZSUM3)
              IF (NRR>=1) DEALLOCATE(ZSUM5)
              !
           ENDIF
           !
           IF (XIBM_SUTR(JI,JJ,JK,2).LT.0.5) THEN
              !
              JIM1 = JI-1
              JJM1 = JJ-1
              JKM1 = JK-1
              JIP1 = JI+1
              JJP1 = JJ+1
              JKP1 = JK+1
              ZSUM1 = 0.
              ZSUM2 = 0.
              !
              DO JI2=JIM1,JIP1
                 DO JJ2=JJM1,JJP1
                    DO JK2=JKM1,JKP1
                       ZSUM1 = ZSUM1 + (XIBM_SUTR(JI2,JJ2,JK2,2))
                       ZSUM2 = ZSUM2 + (XIBM_SUTR(JI2,JJ2,JK2,2))*PRUS(JI2,JJ2,JK2)
                    ENDDO
                 ENDDO
              ENDDO
              !
              PRUS(JI,JJ,JK) = 0.
              IF (ZSUM1.GT.XIBM_EPSI) PRUS(JI,JJ,JK) =  ZSUM2/ZSUM1
              !
           ENDIF
           !
           IF (XIBM_SUTR(JI,JJ,JK,3).LT.0.5) THEN
              !
              JIM1 = JI-1
              JJM1 = JJ-1
              JKM1 = JK-1
              JIP1 = JI+1
              JJP1 = JJ+1
              JKP1 = JK+1
              ZSUM1 = 0.
              ZSUM2 = 0.
              !
              DO JI2=JIM1,JIP1
                 DO JJ2=JJM1,JJP1
                    DO JK2=JKM1,JKP1
                       ZSUM1 = ZSUM1 + (XIBM_SUTR(JI2,JJ2,JK2,3))
                       ZSUM2 = ZSUM2 + (XIBM_SUTR(JI2,JJ2,JK2,3))*PRVS(JI2,JJ2,JK2)
                    ENDDO
                 ENDDO
              ENDDO
              !
              PRVS(JI,JJ,JK) = 0.
              IF (ZSUM1.GT.XIBM_EPSI)     PRVS(JI,JJ,JK) =  ZSUM2/ZSUM1
              !
           ENDIF
           !
           IF (XIBM_SUTR(JI,JJ,JK,4).LT.0.5) THEN
              !
              JIM1 = JI-1
              JJM1 = JJ-1
              JKM1 = JK-1
              JIP1 = JI+1
              JJP1 = JJ+1
              JKP1 = JK+1
              ZSUM1 = 0.
              ZSUM2 = 0.
              !
              DO JI2=JIM1,JIP1
                 DO JJ2=JJM1,JJP1
                    DO JK2=JKM1,JKP1
                       ZSUM1 = ZSUM1 + (XIBM_SUTR(JI2,JJ2,JK2,4))
                       ZSUM2 = ZSUM2 + (XIBM_SUTR(JI2,JJ2,JK2,4))*PRWS(JI2,JJ2,JK2)
                    ENDDO
                 ENDDO
              ENDDO
              !
              PRWS(JI,JJ,JK) = 0.
              IF (ZSUM1.GT.XIBM_EPSI)     PRWS(JI,JJ,JK) =  ZSUM2/ZSUM1
              !
           ENDIF
        ENDDO
     ENDDO
  ENDDO
  !
  PTHS(:,:,IKB-1)=PTHS(:,:,IKB)
  PTHS(:,:,IKE+1)=PTHS(:,:,IKE)
  IF (CTURB/='NONE') PTKS(:,:,IKB-1)=PTKS(:,:,IKB)
  IF (CTURB/='NONE') PTKS(:,:,IKE+1)=PTKS(:,:,IKE)
  IF (NSV>=1) PSVS(:,:,IKB-1,:)=PSVS(:,:,IKB,:)
  IF (NSV>=1) PSVS(:,:,IKE+1,:)=PSVS(:,:,IKE,:)
  IF (NRR>=1) PRRS(:,:,IKB-1,:)=PRRS(:,:,IKB,:)
  IF (NRR>=1) PRRS(:,:,IKE+1,:)=PRRS(:,:,IKE,:)
  PRUS(:,:,IKB-1)=PRUS(:,:,IKB)
  PRUS(:,:,IKE+1)=PRUS(:,:,IKE)
  PRVS(:,:,IKB-1)=PRVS(:,:,IKB)
  PRVS(:,:,IKE+1)=PRVS(:,:,IKE)
  PRWS(:,:,IKB-1)=0.
  PRWS(:,:,IKB)  =0.
  PRWS(:,:,IKE+1)=0.
  !
  NULLIFY(TZFIELDS_ll)
  CALL ADD3DFIELD_ll(TZFIELDS_ll,PTHS(:,:,:),'IBM_FORCING_TR::PTHS')
  IF (CTURB/='NONE') CALL ADD3DFIELD_ll(TZFIELDS_ll,PTKS(:,:,:),'IBM_FORCING_TR::PTKS')
  CALL ADD3DFIELD_ll(TZFIELDS_ll,PRUS(:,:,:),'IBM_FORCING_TR::PRUS')
  CALL ADD3DFIELD_ll(TZFIELDS_ll,PRVS(:,:,:),'IBM_FORCING_TR::PRVS')
  CALL ADD3DFIELD_ll(TZFIELDS_ll,PRWS(:,:,:),'IBM_FORCING_TR::PRWS')
  IF (NSV>=1) THEN
     DO JL=1,NSV
        CALL ADD3DFIELD_ll(TZFIELDS_ll,PSVS(:,:,:,JL),'IBM_FORCING_TR::PSVS')
     ENDDO
  ENDIF
  IF (NRR>=1) THEN
     DO JL=1,NRR
        CALL ADD3DFIELD_ll(TZFIELDS_ll,PRRS(:,:,:,JL),'IBM_FORCING_TR::PRRS')
     ENDDO
  ENDIF
  CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
  CALL CLEANLIST_ll(TZFIELDS_ll)
  !
  ! Problems on corners ? => imposition of the adjacent value
  !
  DO JI=IIB,IIE
     DO JJ=IJB,IJE
        DO JK=IKB,IKE
           !
           IF (XIBM_LS(JI,JJ,JK,2).GT.XIBM_EPSI) THEN
              !
              ZSUM1 = (XIBM_CURV(JI,JJ,JK)+XIBM_CURV(JI-1,JJ,JK))/2.
              ZSUM1 = ABS(ZSUM1)
              ZSUM1 = MIN(1.,ZSUM1)
              !
              JIM1 = JI-1
              JJM1 = JJ-1
              JKM1 = JK-1
              JIP1 = JI+1
              JJP1 = JJ+1
              JKP1 = JK+1
              ZSUM2 = 0.
              !
              DO JI2=JIM1,JIP1
                 DO JJ2=JJM1,JJP1
                    DO JK2=JKM1,JKP1
                       ZSUM2 = ZSUM2 + PRUS(JI2,JJ2,JK2)
                    ENDDO
                 ENDDO
              ENDDO
              !
              PRUS(JI,JJ,JK) = (1.-ZSUM1)*PRUS(JI,JJ,JK)+ZSUM1*ZSUM2/27.
              !
           ENDIF
           !
           IF (XIBM_LS(JI,JJ,JK,3).GT.XIBM_EPSI) THEN
              !
              ZSUM1 = (XIBM_CURV(JI,JJ,JK)+XIBM_CURV(JI,JJ-1,JK))/2.
              ZSUM1 = ABS(ZSUM1)
              ZSUM1 = MIN(1.,ZSUM1)
              !
              JIM1 = JI-1
              JJM1 = JJ-1
              JKM1 = JK-1
              JIP1 = JI+1
              JJP1 = JJ+1
              JKP1 = JK+1
              ZSUM2 = 0.
              !
              DO JI2=JIM1,JIP1
                 DO JJ2=JJM1,JJP1
                    DO JK2=JKM1,JKP1
                       ZSUM2 = ZSUM2 + PRVS(JI2,JJ2,JK2)
                    ENDDO
                 ENDDO
              ENDDO
              !
              PRVS(JI,JJ,JK) =   (1.-ZSUM1)*PRVS(JI,JJ,JK)+ZSUM1*ZSUM2/27.
              !
           ENDIF
           !
           IF (XIBM_LS(JI,JJ,JK,4).GT.XIBM_EPSI) THEN
              !
              ZSUM1 = (XIBM_CURV(JI,JJ,JK)+XIBM_CURV(JI,JJ,JK-1))/2.
              ZSUM1 = ABS(ZSUM1)
              ZSUM1 = MIN(1.,ZSUM1)
              !
              JIM1 = JI-1
              JJM1 = JJ-1
              JKM1 = JK-1
              JIP1 = JI+1
              JJP1 = JJ+1
              JKP1 = JK+1
              ZSUM2 = 0.
              !
              DO JI2=JIM1,JIP1
                 DO JJ2=JJM1,JJP1
                    DO JK2=JKM1,JKP1
                       ZSUM2 = ZSUM2 + PRWS(JI2,JJ2,JK2)
                    ENDDO
                 ENDDO
              ENDDO
              !
              PRWS(JI,JJ,JK) = (1.-ZSUM1)*PRWS(JI,JJ,JK)+ZSUM1*ZSUM2/27.
              !
           ENDIF
        ENDDO
     ENDDO
  ENDDO
  !
  PRUS(:,:,IKB-1)=PRUS(:,:,IKB)
  PRUS(:,:,IKE+1)=PRUS(:,:,IKE)
  PRVS(:,:,IKB-1)=PRVS(:,:,IKB)
  PRVS(:,:,IKE+1)=PRVS(:,:,IKE)
  PRWS(:,:,IKB-1)=0.
  PRWS(:,:,IKB)  =0.
  PRWS(:,:,IKE+1)=0.
  !
  NULLIFY(TZFIELDS_ll)
  CALL ADD3DFIELD_ll(TZFIELDS_ll,PRUS(:,:,:),'IBM_FORCING_TR::PRUS')
  CALL ADD3DFIELD_ll(TZFIELDS_ll,PRVS(:,:,:),'IBM_FORCING_TR::PRVS')
  CALL ADD3DFIELD_ll(TZFIELDS_ll,PRWS(:,:,:),'IBM_FORCING_TR::PRWS')
  CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
  CALL CLEANLIST_ll(TZFIELDS_ll)
  !
  RETURN
  !
END SUBROUTINE IBM_FORCING_TR
