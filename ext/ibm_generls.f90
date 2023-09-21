!MNH_LIC Copyright 2021-2022 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!
!    #######################
MODULE MODI_IBM_GENERLS  
  !    ####################### 
  !
  INTERFACE
     !
     SUBROUTINE IBM_GENERLS(PIBM_FACES,PNORM_FACES,PV1,PV2,PV3,PX_MIN,PY_MIN,PX_MAX,PY_MAX,PPHI)
       !
       REAL, DIMENSION(:,:,:)   ,INTENT(IN)    :: PIBM_FACES  
       REAL, DIMENSION(:,:)     ,INTENT(IN)    :: PNORM_FACES,PV1,PV2,PV3
       REAL, DIMENSION(:,:,:,:) ,INTENT(INOUT) :: PPHI
       REAL                     ,INTENT(IN)    :: PX_MIN,PY_MIN,PX_MAX,PY_MAX
       !
     END SUBROUTINE IBM_GENERLS
     !
  END INTERFACE
  !
END MODULE MODI_IBM_GENERLS
!
!     #####################################
SUBROUTINE IBM_GENERLS(PIBM_FACES,PNORM_FACES,PV1,PV2,PV3,PX_MIN,PY_MIN,PX_MAX,PY_MAX,PPHI)
  !     #####################################
  !
  !
  !****  IBM_GENERLS computes the Level Set function for any surface       
  !                
  !    PURPOSE
  !    -------
  !****  The purpose of this routine is to estimate the level set
  !      containing XYZ minimalisation interface locations

  !    METHOD
  !    ------
  !****  Iterative system and minimization of the interface distance
  !
  !    EXTERNAL
  !    --------
  !      SUBROUTINE ?
  !
  !    IMPLICIT ARGUMENTS
  !    ------------------
  !       MODD_?   
  !
  !    REFERENCE
  !    ---------
  !    The method is based on '3D Distance from a Point to a Triangle'
  !    a technical report from Mark W. Jones, University of Wales Swansea
  !
  !    AUTHORS
  !    ------
  !      Tim Nagel, Valéry Masson & Robert Schoetter 
  !
  !    MODIFICATIONS
  !    -------------
  !      Original         01/06/2021
  !
  !------------------------------------------------------------------------------
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
  USE MODD_IBM_PARAM_n
  USE MODD_IBM_LSF
  USE MODD_DIM_n, ONLY: NIMAX,NJMAX,NKMAX
  USE MODD_PARAMETERS, ONLY: JPVEXT,JPHEXT,XUNDEF
  USE MODD_METRICS_n, ONLY: XDXX,XDYY,XDZZ
  USE MODD_VAR_ll, ONLY: IP
  USE MODD_CST, ONLY: XMNH_EPSILON 
  !
  ! interface
  USE MODI_SHUMAN
  USE MODI_IBM_INTERPOS
  USE MODI_IBM_DETECT
  !
  IMPLICIT NONE
  !
  !       0.1  declarations of arguments
  !                                      
  REAL, DIMENSION(:,:,:)   ,INTENT(IN)    :: PIBM_FACES     !faces coordinates
  REAL, DIMENSION(:,:)     ,INTENT(IN)    :: PNORM_FACES    !normal
  REAL, DIMENSION(:,:)     ,INTENT(IN)    :: PV1,PV2,PV3
  REAL, DIMENSION(:,:,:,:) ,INTENT(INOUT) :: PPHI           ! LS functions
  REAL                     ,INTENT(IN)    :: PX_MIN,PY_MIN,PX_MAX,PY_MAX           
  !
  !------------------------------------------------------------------------------
  !
  !       0.2  declaration of local variables
  !
  INTEGER :: JI,JJ,JK,JN,JM,JI2,JJ2,JK2                                       ! loop index
  INTEGER :: JI_MIN,JI_MAX,JJ_MIN,JJ_MAX,JK_MIN,JK_MAX,IIU,IJU,IKU            ! loop boundaries
  REAL                                :: Z_DIST_TEST1,Z_DIST_TEST2            ! saving distances  
  REAL                                :: Z_DIST_TEST3,Z_DIST_TEST4,ZDIST_REF0
  INTEGER                             :: INUMB_FACES                          ! number of faces 
  REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZXHATM,ZYHATM,ZZHATM,ZDP0PP0PAST 
  CHARACTER(LEN=1)                    :: YPOS
  REAL, DIMENSION(3)                  :: ZP1P0,ZP1P2,ZP0PP0,ZP1PP0,ZP2PP0,ZP3PP0,ZPP0P1,ZPP0P2,ZPP0P3
  REAL, DIMENSION(3)                  :: ZPP0PPP0,ZPPP0P1,ZPPP0P2,ZP2P1,ZP2P0,ZP2P3,ZP3P2,ZP3P1
  REAL, DIMENSION(3)                  :: ZPP0,ZFT1,ZFT2,ZFT3,ZFT1B,ZFT2B,ZFT3B,ZR,ZPPP0,ZP3P0,ZP0P1
  REAL, DIMENSION(3)                  :: ZPPP0P3,ZP1P3,ZPCP0,ZR0
  REAL, DIMENSION(:), ALLOCATABLE     :: ZSTEMP,ZRDIR,ZVECTDISTPLUS,ZVECTDISTMOINS,ZVECTDIST!,ZFACE
  REAL, DIMENSION(:,:), ALLOCATABLE   :: ZC
  REAL                                :: ZF1,ZF2,ZF3,ZF1B,ZF2B,ZF3B,ZDPP0PPP0
  REAL                                :: ZT,ZSIGN,ZS,ZDIST,ZDP0PP0,ZNNORM,ZRN,ZPHI_OLD
  TYPE(LIST_ll), POINTER              :: TZFIELDS_ll   ! list of fields to exchange
  INTEGER                             :: IINFO_ll,IMI  ! return code of parallel routine
  INTEGER                             :: IIE,IIB,IJB,IJE,IKE,IKB,ZBPLUS
  LOGICAL                             :: GABOVE_ROOF,LFACE,LDZ
  LOGICAL, DIMENSION(:), ALLOCATABLE  :: ZFACE
  INTEGER                             :: ZCOUNT,ZIDX,ZII,ZCHANGE,ZCHANGE1
  REAL                                :: ZDIFF,ZMIN_DIFF,ZDX
  !
  !------------------------------------------------------------------------------
  !   
  !       0.3 allocation
  !
  NULLIFY(TZFIELDS_ll)
  IIU = SIZE(PPHI,1)
  IJU = SIZE(PPHI,2)
  IKU = SIZE(PPHI,3)
  IIB=1+JPHEXT
  IIE=IIU-JPHEXT
  IJB=1+JPHEXT
  IJE=IJU-JPHEXT
  IKB=1+JPVEXT
  IKE=IKU-JPVEXT
  !
  JK_MIN = 1 + JPVEXT
  JK_MAX = IKU - JPVEXT
  !
  CALL GET_INDICE_ll (JI_MIN,JJ_MIN,JI_MAX,JJ_MAX)
  !
  ALLOCATE(ZXHATM(IIU,IJU,IKU))
  ALLOCATE(ZYHATM(IIU,IJU,IKU))
  ALLOCATE(ZZHATM(IIU,IJU,IKU))
  !
  !-------------------------------------------------------------------------------
  !
  !**** 1. PRELIMINARIES
  !     ----------------
  !
  INUMB_FACES = SIZE(PIBM_FACES,1)
  ALLOCATE(ZC(INUMB_FACES,3))
  ALLOCATE(ZSTEMP(1))
  ALLOCATE(ZRDIR(1))  
  PPHI = -XUNDEF
  ALLOCATE(ZDP0PP0PAST(IIU,IJU,IKU))
  ZDP0PP0PAST = 0.
  ALLOCATE(ZVECTDIST(10000))
  ALLOCATE(ZVECTDISTPLUS(10000))
  ALLOCATE(ZVECTDISTMOINS(10000))
  ALLOCATE(ZFACE(10000))
  ZFACE=.FALSE.
  !      
  !-------------------------------------------------------------------------------
  !
  !**** 2. EXECUTIONS
  !     -------------
  !
  JM=1
  YPOS = 'P'
  !
  CALL IBM_INTERPOS(ZXHATM,ZYHATM,ZZHATM,YPOS)
  ZDX = ZXHATM(JI_MIN+1,JJ_MIN,JK_MIN)-ZXHATM(JI_MIN,JJ_MIN,JK_MIN)
  !
  DO JK = JK_MIN,JK_MAX
     DO JJ = JJ_MIN,JJ_MAX
        DO JI = JI_MIN,JI_MAX
           ZCOUNT = 1
           ZVECTDIST = -999.
           DO JN = 1,INUMB_FACES
              LFACE=.FALSE.
              !***Calcul of the face center
              ZC(JN,1)=(PIBM_FACES(JN,1,1)+PIBM_FACES(JN,2,1)+PIBM_FACES(JN,3,1))/3.
              ZC(JN,2)=(PIBM_FACES(JN,1,2)+PIBM_FACES(JN,2,2)+PIBM_FACES(JN,3,2))/3.
              ZC(JN,3)=(PIBM_FACES(JN,1,3)+PIBM_FACES(JN,2,3)+PIBM_FACES(JN,3,3))/3.
              !***Norm normalization
              ZNNORM = SQRT(PNORM_FACES(JN,1)**2+PNORM_FACES(JN,2)**2+PNORM_FACES(JN,3)**2)
              !***Vector between the face center and the current grid point
              ZPCP0(1) = ZXHATM(JI,JJ,JK)-ZC(JN,1)
              ZPCP0(2) = ZYHATM(JI,JJ,JK)-ZC(JN,2)
              ZPCP0(3) = ZZHATM(JI,JJ,JK)-ZC(JN,3)
              ZSIGN = ZPCP0(1)*PNORM_FACES(JN,1)+ &
                      ZPCP0(2)*PNORM_FACES(JN,2)+ &
                      ZPCP0(3)*PNORM_FACES(JN,3)
              !***Various vectors
              ZP1P0(1) = ZXHATM(JI,JJ,JK)-PIBM_FACES(JN,1,1)
              ZP1P0(2) = ZYHATM(JI,JJ,JK)-PIBM_FACES(JN,1,2)
              ZP1P0(3) = ZZHATM(JI,JJ,JK)-PIBM_FACES(JN,1,3)
              ZP3P0(1) = ZXHATM(JI,JJ,JK)-PIBM_FACES(JN,3,1)
              ZP3P0(2) = ZYHATM(JI,JJ,JK)-PIBM_FACES(JN,3,2)
              ZP3P0(3) = ZZHATM(JI,JJ,JK)-PIBM_FACES(JN,3,3)
              ZP0P1(1) = PIBM_FACES(JN,1,1)-ZXHATM(JI,JJ,JK)
              ZP0P1(2) = PIBM_FACES(JN,1,2)-ZYHATM(JI,JJ,JK)
              ZP0P1(3) = PIBM_FACES(JN,1,3)-ZZHATM(JI,JJ,JK)
              ZP2P0(1) = ZXHATM(JI,JJ,JK)-PIBM_FACES(JN,2,1)
              ZP2P0(2) = ZYHATM(JI,JJ,JK)-PIBM_FACES(JN,2,2)
              ZP2P0(3) = ZZHATM(JI,JJ,JK)-PIBM_FACES(JN,2,3)
              !***Equation (3) of Jones (1995)
              IF(ZP1P0(1)==0.AND.ZP1P0(2)==0.AND.ZP1P0(3)==0) THEN
                 WRITE(*,*) 'ZP1P0(1,2,3)',ZP1P0(1),ZP1P0(2),ZP1P0(3)
                 ZDP0PP0 = 0.
              ELSE
                 ZDP0PP0 = SQRT(ZP0P1(1)**2+ZP0P1(2)**2+ZP0P1(3)**2)* &
                        ((ZP1P0(1)*PNORM_FACES(JN,1)+ZP1P0(2)*PNORM_FACES(JN,2)+&
                         ZP1P0(3)*PNORM_FACES(JN,3))/( &
                        SQRT((ZP1P0(1))**2+(ZP1P0(2))**2+(ZP1P0(3))**2)*ZNNORM))
              END IF
              !***Equation (4) of Jones (1995)
              ZP0PP0(1) = -ZDP0PP0*(PNORM_FACES(JN,1)/ZNNORM)
              ZP0PP0(2) = -ZDP0PP0*(PNORM_FACES(JN,2)/ZNNORM)
              ZP0PP0(3) = -ZDP0PP0*(PNORM_FACES(JN,3)/ZNNORM)
              !***Equation (5) of Jones (1995)
              ZPP0(1) = ZXHATM(JI,JJ,JK)+ZP0PP0(1)
              ZPP0(2) = ZYHATM(JI,JJ,JK)+ZP0PP0(2)
              ZPP0(3) = ZZHATM(JI,JJ,JK)+ZP0PP0(3)
              !
              ZP1PP0(1)=ZPP0(1)-PIBM_FACES(JN,1,1)
              ZP1PP0(2)=ZPP0(2)-PIBM_FACES(JN,1,2)
              ZP1PP0(3)=ZPP0(3)-PIBM_FACES(JN,1,3)
              !
              ZP2PP0(1)=ZPP0(1)-PIBM_FACES(JN,2,1)
              ZP2PP0(2)=ZPP0(2)-PIBM_FACES(JN,2,2)
              ZP2PP0(3)=ZPP0(3)-PIBM_FACES(JN,2,3)
              !
              ZP3PP0(1)=ZPP0(1)-PIBM_FACES(JN,3,1)
              ZP3PP0(2)=ZPP0(2)-PIBM_FACES(JN,3,2)
              ZP3PP0(3)=ZPP0(3)-PIBM_FACES(JN,3,3)
              !
              ZPP0P1(1)=PIBM_FACES(JN,1,1)-ZPP0(1)
              ZPP0P1(2)=PIBM_FACES(JN,1,2)-ZPP0(2)
              ZPP0P1(3)=PIBM_FACES(JN,1,3)-ZPP0(3)
              !
              ZPP0P2(1)=PIBM_FACES(JN,2,1)-ZPP0(1)
              ZPP0P2(2)=PIBM_FACES(JN,2,2)-ZPP0(2)
              ZPP0P2(3)=PIBM_FACES(JN,2,3)-ZPP0(3)
              !
              ZPP0P3(1)=PIBM_FACES(JN,3,1)-ZPP0(1)
              ZPP0P3(2)=PIBM_FACES(JN,3,2)-ZPP0(2)
              ZPP0P3(3)=PIBM_FACES(JN,3,3)-ZPP0(3)
              !
              !***Calculation of f1,f2,f3 (Jones (1995))
              ZFT1= CROSSPRODUCT(PV1(JN,:),ZP1PP0)
              ZFT2= CROSSPRODUCT(PV2(JN,:),ZP2PP0)
              ZFT3= CROSSPRODUCT(PV3(JN,:),ZP3PP0)

              ZF1 =ZFT1(1)*PNORM_FACES(JN,1)+ &
                   ZFT1(2)*PNORM_FACES(JN,2)+ &
                   ZFT1(3)*PNORM_FACES(JN,3)

              ZF2 =ZFT2(1)*PNORM_FACES(JN,1)+ &
                   ZFT2(2)*PNORM_FACES(JN,2)+ &
                   ZFT2(3)*PNORM_FACES(JN,3)

              ZF3 =ZFT3(1)*PNORM_FACES(JN,1)+ &
                   ZFT3(2)*PNORM_FACES(JN,2)+ &
                   ZFT3(3)*PNORM_FACES(JN,3)
              !***Point anticlockwise of V1 and clockwise of V2
              IF (ZF1.GE.0.AND.ZF2.LE.0) THEN
                 ZFT1B = CROSSPRODUCT(ZPP0P1,ZPP0P2)
                 ZF1B = ZFT1B(1)*PNORM_FACES(JN,1)+ &
                      ZFT1B(2)*PNORM_FACES(JN,2)+ &
                      ZFT1B(3)*PNORM_FACES(JN,3)
                 IF (ZF1B<0) THEN
                    ZP1P2(:) = PIBM_FACES(JN,2,:)-PIBM_FACES(JN,1,:)
                    ZR = CROSSPRODUCT(CROSSPRODUCT(ZPP0P2,ZPP0P1),ZP1P2)
                    ZRN = SQRT(ZR(1)**2+ZR(2)**2+ZR(3)**2)    
                    !***Eq. (10) of Jones(1995)
                    ZDPP0PPP0 = SQRT(ZPP0P1(1)**2+ZPP0P1(2)**2+ZPP0P1(3)**2)* &
                                ((ZPP0P1(1)*ZR(1)+ZPP0P1(2)*ZR(2)+ZPP0P1(3)*ZR(3))/( &
                                SQRT(ZPP0P1(1)**2+ZPP0P1(2)**2+ZPP0P1(3)**2)*ZRN))! &
                    ZPP0PPP0 = ZDPP0PPP0*(ZR/ZRN)
                    ZPPP0 = ZPP0+ZPP0PPP0
                    ZPPP0P1 = PIBM_FACES(JN,1,:)-ZPPP0
                    ZP2P1 = PIBM_FACES(JN,1,:)-PIBM_FACES(JN,2,:)
                    ZRDIR = SIGN(1.,SCALPRODUCT(ZPPP0P1,ZP2P1))
                    ZT = SQRT(ZPPP0P1(1)**2+ZPPP0P1(2)**2+ZPPP0P1(3)**2)/ &
                         SQRT(ZP2P1(1)**2+ZP2P1(2)**2+ZP2P1(3)**2)*ZRDIR(1)
                    IF (ZT.GE.0.AND.ZT.LE.1) THEN
                       ZDIST =SQRT(ZDPP0PPP0**2+ZDP0PP0**2)
                    ELSEIF (ZT<0.) THEN
                       ZDIST = SQRT(ZP1P0(1)**2+ZP1P0(2)**2+ZP1P0(3)**2)
                    ELSEIF (ZT>1.) THEN
                       ZDIST = SQRT(ZP2P0(1)**2+ZP2P0(2)**2+ZP2P0(3)**2)
                    ELSE
                       call Print_msg( NVERB_FATAL, 'GEN', 'IBM_PREP_LS', 'Error in ZT calculation' )
                    ENDIF
                 ELSE
                    ZDIST = ZDP0PP0 
                    LFACE = .TRUE.                    
                 ENDIF
              !***Point anticlockwise of V2 and clockwise of V3
              ELSEIF (ZF2.GE.0.AND.ZF3.LE.0) THEN
                 ZFT2B = CROSSPRODUCT(ZPP0P2,ZPP0P3)
                 ZF2B = ZFT2B(1)*PNORM_FACES(JN,1)+ &
                      ZFT2B(2)*PNORM_FACES(JN,2)+ &
                      ZFT2B(3)*PNORM_FACES(JN,3)
                 IF (ZF2B<0) THEN
                    ZP2P3(:) = PIBM_FACES(JN,3,:)-PIBM_FACES(JN,2,:)
                    ZR = CROSSPRODUCT(CROSSPRODUCT(ZPP0P3,ZPP0P2),ZP2P3)
                    ZRN = SQRT(ZR(1)**2+ZR(2)**2+ZR(3)**2)
                    ZDPP0PPP0 = SQRT(ZPP0P2(1)**2+ZPP0P2(2)**2+ZPP0P2(3)**2)* &
                         ((ZPP0P2(1)*ZR(1)+ZPP0P2(2)*ZR(2)+ZPP0P2(3)*ZR(3))/( &
                         SQRT(ZPP0P2(1)**2+ZPP0P2(2)**2+ZPP0P2(3)**2)*ZRN))! &
                    ZPP0PPP0 = ZDPP0PPP0*(ZR/ZRN)
                    ZPPP0 = ZPP0+ZPP0PPP0
                    ZPPP0P2 = PIBM_FACES(JN,2,:)-ZPPP0
                    ZP3P2 = PIBM_FACES(JN,2,:)-PIBM_FACES(JN,3,:)
                    ZRDIR = SIGN(1.,SCALPRODUCT(ZPPP0P2,ZP3P2))
                    ZT = SQRT(ZPPP0P2(1)**2+ZPPP0P2(2)**2+ZPPP0P2(3)**2)/ &
                         SQRT(ZP3P2(1)**2+ZP3P2(2)**2+ZP3P2(3)**2)*ZRDIR(1)
                    IF (ZT.GE.0.AND.ZT.LE.1) THEN
                       ZDIST = SQRT(ZDPP0PPP0**2+ZDP0PP0**2)
                    ELSEIF (ZT<0.) THEN
                       ZDIST = SQRT(ZP2P0(1)**2+ZP2P0(2)**2+ZP2P0(3)**2)
                    ELSEIF (ZT>1.) THEN
                       ZDIST = SQRT(ZP3P0(1)**2+ZP3P0(2)**2+ZP3P0(3)**2)
                    ELSE
                       call Print_msg( NVERB_FATAL, 'GEN', 'IBM_PREP_LS', 'Error in ZT calculation' )
                    ENDIF
                 ELSE
                    ZDIST = ZDP0PP0
                    LFACE = .TRUE.                    
                 ENDIF
              !***Point anticlockwise of V3 and clockwise of V1
              ELSEIF (ZF3.GE.0.AND.ZF1.LE.0) THEN
                 ZFT3B = CROSSPRODUCT(ZPP0P3,ZPP0P1)
                 ZF3B = ZFT3B(1)*PNORM_FACES(JN,1)+ &
                      ZFT3B(2)*PNORM_FACES(JN,2)+ &
                      ZFT3B(3)*PNORM_FACES(JN,3)
                 IF (ZF3B<0) THEN
                    ZP3P1(:) = PIBM_FACES(JN,1,:)-PIBM_FACES(JN,3,:)
                    ZR = CROSSPRODUCT(CROSSPRODUCT(ZPP0P1,ZPP0P3),ZP3P1)
                    ZRN = SQRT(ZR(1)**2+ZR(2)**2+ZR(3)**2)
                    ZDPP0PPP0 = SQRT(ZPP0P3(1)**2+ZPP0P3(2)**2+ZPP0P3(3)**2)* &
                         ((ZPP0P3(1)*ZR(1)+ZPP0P3(2)*ZR(2)+ZPP0P3(3)*ZR(3))/( &
                         SQRT((ZPP0P3(1))**2+(ZPP0P3(2))**2+(ZPP0P3(3))**2)*ZRN))! &
                    ZPP0PPP0 = ZDPP0PPP0*(ZR/ZRN)
                    ZPPP0 = ZPP0+ZPP0PPP0
                    ZPPP0P3 = PIBM_FACES(JN,3,:)-ZPPP0
                    ZP1P3 = PIBM_FACES(JN,3,:)-PIBM_FACES(JN,1,:)
                    ZRDIR = SIGN(1.,SCALPRODUCT(ZPPP0P3,ZP1P3))
                    ZT = SQRT(ZPPP0P3(1)**2+ZPPP0P3(2)**2+ZPPP0P3(3)**2)/ &
                         SQRT(ZP1P3(1)**2+ZP1P3(2)**2+ZP1P3(3)**2)*ZRDIR(1)
                    IF (ZT.GE.0.AND.ZT.LE.1) THEN
                       ZDIST = SQRT(ZDPP0PPP0**2+ZDP0PP0**2)
                    ELSEIF (ZT<0.) THEN
                       ZDIST = SQRT(ZP3P0(1)**2+ZP3P0(2)**2+ZP3P0(3)**2)
                    ELSEIF (ZT>1.) THEN
                       ZDIST = SQRT(ZP1P0(1)**2+ZP1P0(2)**2+ZP1P0(3)**2)
                    ELSE
                       call Print_msg( NVERB_FATAL, 'GEN', 'IBM_PREP_LS', 'Error in ZT calculation' )
                    ENDIF
                 ELSE
                    ZDIST = ZDP0PP0
                    LFACE = .TRUE.
                 ENDIF
              ELSE
                 call Print_msg( NVERB_FATAL, 'GEN', 'IBM_PREP_LS', 'Error in ZF instruction' )
              ENDIF
              ZDIST = SIGN(ZDIST,-ZSIGN)
              ZDIST = ANINT(ZDIST*10.E5) / 10.E5
              PPHI(JI,JJ,JK,JM) = ANINT(PPHI(JI,JJ,JK,JM)*10.E5) / 10.E5
              IF (ABS(ZDIST).LE.ABS(PPHI(JI,JJ,JK,JM))) THEN
                 ZPHI_OLD = PPHI(JI,JJ,JK,JM)
                 IF (ABS(ZDIST)==ABS(PPHI(JI,JJ,JK,JM))) THEN
                    IF (ABS(ZDP0PP0).GT.ABS(ZDP0PP0PAST(JI,JJ,JK))) THEN
                       PPHI(JI,JJ,JK,JM) = ZDIST
                       ZDP0PP0PAST(JI,JJ,JK) = ZDP0PP0
                    ENDIF
                 ELSE
                    PPHI(JI,JJ,JK,JM) = ZDIST
                 ENDIF
                 IF (ABS(ZDIST).LT.ABS(ZPHI_OLD)) THEN        
                    ZDP0PP0PAST(JI,JJ,JK) = ZDP0PP0
                 ENDIF
              ENDIF
              IF (ABS(PPHI(JI,JJ,JK,JM)).GT.(SQRT(3.)*4.)) THEN
                 PPHI(JI,JJ,JK,JM) = -999.
              ENDIF
              IF (ABS(ZDIST).LT.(SQRT(3.)*4.)) THEN
                 ZVECTDIST(ZCOUNT)=ZDIST
                 ZFACE(ZCOUNT)=LFACE
                 ZCOUNT = ZCOUNT +1
              ENDIF
           ENDDO
           ZVECTDISTPLUS=ZVECTDIST
           ZVECTDISTMOINS=ZVECTDIST
           WHERE (ZVECTDIST.GT.0)
                 ZVECTDISTMOINS=-999. 
           ENDWHERE
           WHERE (ZVECTDIST.LT.0)
                 ZVECTDISTPLUS=999.
           ENDWHERE
           IF (ANY(ZVECTDIST.GT.0.).AND.(ABS(ABS(MINVAL(ZVECTDISTPLUS))-ABS(MAXVAL(ZVECTDISTMOINS))).LT.10.E-6)) THEN
              ZMIN_DIFF = 1.
              ZIDX = 0
              DO ZII = 1, SIZE(ZVECTDIST)
                 ZDIFF = ABS(ZVECTDIST(ZII)-MINVAL(ZVECTDISTPLUS))
                 IF ( ZDIFF < ZMIN_DIFF) THEN
                    ZIDX = ZII
                    ZMIN_DIFF = ZDIFF
                 ENDIF
               ENDDO
               IF (ZFACE(ZIDX)) THEN
                  PPHI(JI,JJ,JK,JM) = MINVAL(ZVECTDISTPLUS)
               ENDIF   
           ENDIF
        ENDDO
     ENDDO
  ENDDO

DO JJ=JJ_MIN,JJ_MAX
DO JI=JI_MIN,JI_MAX
GABOVE_ROOF=.FALSE.
DO JK=IKB, IKE
  ! check if point is flagged as not calculated
  IF (PPHI(JI,JJ,JK,JM)==-999.) THEN
     ! check if point is already above a point that encountered a point near the
     ! surface (that can be outside or inside a building)
     ! check if that point was inside (if outside, the value of the levelset
     ! stays at -999.)
     IF (GABOVE_ROOF .AND. PPHI(JI,JJ,JK-1,JM) > XIBM_EPSI) THEN
       PPHI(JI,JJ,JK,JM) = 999.
       CYCLE
     END IF
    ! check if the point of the column have not encoutered a near-building
    ! surface point with a physical value of the level set
    IF (.NOT. GABOVE_ROOF) THEN
      ! if the point above has a physical value for the level set, then the
      ! status inside (999) or outside (-999) is given to all points below,
      ! depending if this point above (that needs not to be the point at the top
      ! of the model!) is inside or outside
      ! checks if the point above has a physical value for the levelset                   
      IF (JK<IKE .AND. ABS (PPHI(JI,JJ,JK+1,JM)) < 900.) THEN
         ! if the point above is inside, all points below are set inside
         IF (PPHI(JI,JJ,JK+1,JM)>XIBM_EPSI) PPHI(JI,JJ,IKB:JK,JM) = 999.
         ! indicate for further processing of points above the current point
         ! that we have encountered a physical value of the level set, near the
         ! surface building
         GABOVE_ROOF = .TRUE.
      END IF
      CYCLE
    ENDIF
  END IF
  ! if we have never encoutered a roof or point near a building form above,
  ! then, we are outside, and nothing is changed (value -999 kept)
  END DO
  PPHI(JI,JJ,IKB-1,JM) = PPHI(JI,JJ,IKB,JM)
  PPHI(JI,JJ,IKE+1,JM) = PPHI(JI,JJ,IKE,JM)
END DO
END DO


JN=1
PPHI(:,:,IKB-1,JN)=2*PPHI(:,:,IKB,JN)-PPHI(:,:,IKB+1,JN)
PPHI(:,:,IKE+1,JN)=2*PPHI(:,:,IKE,JN)-PPHI(:,:,IKE-1,JN)
PPHI(IIB-1,:,:,JN) = PPHI(    IIB  ,:,:,JN)
PPHI(IIE+1,:,:,JN) = PPHI(    IIE  ,:,:,JN)
PPHI(:,IJB-1,:,JN) = PPHI(:,    IJB  ,:,JN)
PPHI(:,IJE+1,:,JN) = PPHI(:,    IJE  ,:,JN)

PPHI(:,:,:,2)=MXM(PPHI(:,:,:,1))
PPHI(:,:,:,3)=MYM(PPHI(:,:,:,1))
PPHI(:,:,:,4)=MZM(PPHI(:,:,:,1))

NULLIFY(TZFIELDS_ll)
DO JN=2,4
  PPHI(:,:,IKB-1,JN)=2*PPHI(:,:,IKB,JN)-PPHI(:,:,IKB+1,JN)
  PPHI(:,:,IKE+1,JN)=2*PPHI(:,:,IKE,JN)-PPHI(:,:,IKE-1,JN)
    PPHI(IIB-1,:,:,JN) = PPHI(    IIB  ,:,:,JN)
    PPHI(IIE+1,:,:,JN) = PPHI(    IIE  ,:,:,JN)
    PPHI(:,IJB-1,:,JN) = PPHI(:,    IJB  ,:,JN)
    PPHI(:,IJE+1,:,JN) = PPHI(:,    IJE  ,:,JN)
ENDDO

PPHI(:,:,:,5)=MYM(PPHI(:,:,:,2))
PPHI(:,:,:,6)=MXM(PPHI(:,:,:,4))
PPHI(:,:,:,7)=MYM(PPHI(:,:,:,4))
NULLIFY(TZFIELDS_ll)
DO JN=5,7
  PPHI(:,:,IKB-1,JN)=2*PPHI(:,:,IKB,JN)-PPHI(:,:,IKB+1,JN)
  PPHI(:,:,IKE+1,JN)=2*PPHI(:,:,IKE,JN)-PPHI(:,:,IKE-1,JN)
    PPHI(IIB-1,:,:,JN) = PPHI(    IIB  ,:,:,JN)
    PPHI(IIE+1,:,:,JN) = PPHI(    IIE  ,:,:,JN)
    PPHI(:,IJB-1,:,JN) = PPHI(:,    IJB  ,:,JN)
    PPHI(:,IJE+1,:,JN) = PPHI(:,    IJE  ,:,JN)
ENDDO
WHERE (ABS(PPHI(:,:,:,:)).LT.XIBM_EPSI) PPHI(:,:,:,:)=2.*XIBM_EPSI


  !COMPLETE PPHI ON THE HALO OF EACH SUBDOMAINS
  DO JN=1,7
     CALL ADD3DFIELD_ll(TZFIELDS_ll,PPHI(:,:,:,JN),'IBM_GENERLS::PPHI')
  ENDDO
  CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
  CALL CLEANLIST_ll(TZFIELDS_ll)

  !
  !-------------------------------------------------------------------------------
  !
  !**** X. DEALLOCATIONS/CLOSES
  !     -----------------------
  !  
  !DEALLOCATE(ZDP0PP0,ZDIST,ZC,ZSTEMP)
  DEALLOCATE(ZC,ZSTEMP)
  DEALLOCATE(ZXHATM,ZYHATM,ZZHATM)              
  !
  RETURN
  !
CONTAINS
  !
  FUNCTION CROSSPRODUCT(PA,PB) RESULT(CROSS)
    !
    REAL, DIMENSION(3)             :: CROSS
    REAL, DIMENSION(3), INTENT(IN) :: PA, PB
    CROSS(1) = PA(2) * PB(3) - PA(3) * PB(2)
    CROSS(2) = PA(3) * PB(1) - PA(1) * PB(3)
    CROSS(3) = PA(1) * PB(2) - PA(2) * PB(1)
  END FUNCTION CROSSPRODUCT

  FUNCTION SCALPRODUCT(PA,PB) RESULT(SCAL)
    !
    REAL                           :: SCAL
    REAL, DIMENSION(3), INTENT(IN) :: PA, PB
    SCAL = PA(1)*PB(1)+PA(2)*PB(2)+PA(3)*PB(3)
  END FUNCTION SCALPRODUCT

END SUBROUTINE IBM_GENERLS
