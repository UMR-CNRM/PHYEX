!MNH_LIC Copyright 2018-2020 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_LIMA_SEDIMENTATION
  IMPLICIT NONE
CONTAINS
!     ######################################################################
  SUBROUTINE LIMA_SEDIMENTATION (LIMAP, LIMAW, LIMAC, LIMAM, D, CST, &
                                 HPHASE, KMOMENTS, KID, KISHAPE, KSPLITG, PTSTEP, OELEC, ELECP, ELECD, &
                                 PDZZ, PRHODREF, PTHVREFZIKB, &
                                 PPABST, PT, PRT_SUM, PCPT, PRS, PCS, PINPR, PFPR, &
                                 PEFIELDW, PQS, PLBDAI_SHAPE)
!     ######################################################################
!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the sedimentation of any hydrometeor,
!!    also accounting for the transport of heat
!!
!!    METHOD
!!    ------
!!      The sedimentation rates are computed with a time spliting technique: 
!!    an upstream scheme, written as a difference of non-advective fluxes. 
!!    This source term is added to the next coming time step (split-implicit 
!!    process).
!!
!!    AUTHOR
!!    ------
!!      J.-M. Cohard     * Laboratoire d'Aerologie*
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!      S.    Berthet    * Laboratoire d'Aerologie*
!!      B.    Vié        * CNRM *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original             15/03/2018
!!
!!      B.Vie  02/2019  Desactivate (comment) the heat transport by droplets
!  P. Wautelet 26/04/2019: replace non-standard FLOAT function by REAL function
!  P. Wautelet 28/05/2019: move COUNTJV function to tools.f90
!  B. Vie         03/2020: disable temperature change of droplets by air temperature
!  J. Wurtz       03/2022: new snow characteristics
!  C. Barthe   03/06/2022: add sedimentation for electric charges
!  C. Barthe   02/06/2023: add the Beard effect (electric field)
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_DIMPHYEX,         ONLY: DIMPHYEX_T
USE MODD_CST,              ONLY: CST_T
USE MODD_RAIN_ICE_DESCR_N, ONLY: RAIN_ICE_DESCR_T
USE MODD_ELEC_DESCR,       ONLY: ELEC_DESCR_t
USE MODD_ELEC_PARAM,       ONLY: ELEC_PARAM_t

USE MODE_TOOLS,            only: COUNTJV

USE MODE_ELEC_COMPUTE_EX,   ONLY: ELEC_COMPUTE_EX
USE MODE_ELEC_BEARD_EFFECT, ONLY: ELEC_BEARD_EFFECT
USE MODD_PARAM_LIMA_MIXED, ONLY:PARAM_LIMA_MIXED_T
USE MODD_PARAM_LIMA_COLD, ONLY:PARAM_LIMA_COLD_T
USE MODD_PARAM_LIMA_WARM, ONLY:PARAM_LIMA_WARM_T
USE MODD_PARAM_LIMA, ONLY:PARAM_LIMA_T
USE YOMHOOK, ONLY:LHOOK, DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(DIMPHYEX_T),              INTENT(IN)    :: D
TYPE(CST_T),                   INTENT(IN)    :: CST
CHARACTER(1),                  INTENT(IN)    :: HPHASE    ! Liquid or solid hydrometeors
INTEGER,                       INTENT(IN)    :: KMOMENTS  ! Number of moments 
INTEGER,                       INTENT(IN)    :: KID       ! Hydrometeor ID
INTEGER,                       INTENT(IN)    :: KISHAPE   ! Ice shape ID if LCRYSTAL_SHAPE (0 otherwise)
INTEGER,                       INTENT(IN)    :: KSPLITG   !  
REAL,                          INTENT(IN)    :: PTSTEP    ! Time step  
LOGICAL,                       INTENT(IN)    :: OELEC     ! if true, cloud electrification is activated
TYPE(ELEC_PARAM_t),            INTENT(IN)    :: ELECP   ! electrical parameters
TYPE(ELEC_DESCR_t),            INTENT(IN)    :: ELECD   ! electrical descriptive csts
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PDZZ      ! Height (z)
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PRHODREF  ! Reference density
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PPABST    ! abs. pressure at time t
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(INOUT) :: PT        ! Temperature
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(INOUT) :: PRT_SUM   ! total water mixing ratio
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(INOUT) :: PCPT      ! Cp
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(INOUT) :: PRS       ! m.r. source
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(INOUT) :: PCS       ! C. source
REAL, DIMENSION(D%NIJT),       INTENT(OUT)   :: PINPR     ! Instant precip rate
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)   :: PFPR      ! Precip. fluxes in altitude
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN),    OPTIONAL :: PEFIELDW  ! Vertical component of the electric field
REAL,                          INTENT(IN)    :: PTHVREFZIKB ! Reference thv at IKB for electricity
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(INOUT), OPTIONAL :: PQS ! Elec. charge density source
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN), OPTIONAL :: PLBDAI_SHAPE ! lambda for one ice crystal shape
!
!*       0.2   Declarations of local variables :
!
INTEGER :: IK, IL, IN                     ! Loop index
INTEGER :: ISEDIM                         ! Case number of sedimentation
!
LOGICAL, DIMENSION(D%NIJT, D%NKT) :: GSEDIM      ! Test where to compute the SED processes
REAL,    DIMENSION(D%NIJT, D%NKT) :: ZW,       & ! Work array
                                     ZWDT        ! Temperature change
REAL,    DIMENSION(D%NIJT,0:D%NKT+1) &
                           :: ZWSEDR,   & ! Sedimentation of MMR
                              ZWSEDC      ! Sedimentation of number conc.
!
REAL, DIMENSION(:), ALLOCATABLE         &
                           :: ZRS,      & ! m.r. source
                              ZCS,      & ! conc. source
                              ZRHODREF, & ! RHO Dry REFerence
                              ZPABST,   & ! Pressure
                              ZT,       & ! Temperature
                              ZZW,      & ! Work array
                              ZZX,      & ! Work array
                              ZZY,      & ! Work array
                              ZLBDA,    & ! Slope parameter
                              ZCC         ! Cunningham corrective term for droplets fall speed
!
INTEGER , DIMENSION(D%NIJT*D%NKT) :: I1,I2,I3 ! Indexes for PACK replacement
!
REAL    :: ZTSPLITG                       ! Small time step for rain sedimentation
REAL    :: ZC                             ! Cpl or Cpi
INTEGER :: IMOMENTS
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
! Variables for cloud electricity
REAL :: ZCX, ZXX  ! C and x parameters for N-lambda relationship
REAL, DIMENSION(:),     ALLOCATABLE :: ZQS, &    ! Electric charge density source
                                       ZZQ, &    ! Work array
                                       ZES       ! e in q-D relationship
REAL, DIMENSION(MERGE(D%NIJT, 0, OELEC), &
               &MERGE(D%NKT, 0, OELEC)) :: ZWSEDQ, &   ! Sedimentation of electric charge density
                                         & ZLBDA3
TYPE(PARAM_LIMA_MIXED_T),INTENT(IN)::LIMAM
TYPE(PARAM_LIMA_COLD_T),INTENT(IN)::LIMAC
TYPE(PARAM_LIMA_WARM_T),INTENT(IN)::LIMAW
TYPE(PARAM_LIMA_T),INTENT(IN)::LIMAP
REAL, DIMENSION(D%NIJT, D%NKT):: ZBEARDCOEFF ! effect of electrical forces on terminal fall speed
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('LIMA_SEDIMENTATION', 0, ZHOOK_HANDLE)
IMOMENTS=KMOMENTS
!
! Time splitting
!
ZTSPLITG= PTSTEP / REAL(LIMAP%NSPLITSED(KID))
!
ZWDT(:,:)=0.
PINPR(:) = 0.
PFPR(:,:) = 0.
ZWSEDR(:,:) = 0.
ZWSEDC(:,:) = 0.
!
PRS(:,:) = PRS(:,:) * PTSTEP
IF (KMOMENTS==2) PCS(:,:) = PCS(:,:) * PTSTEP
IF (OELEC)       PQS(:,:) = PQS(:,:) * PTSTEP
DO IK = D%NKTB , D%NKTE
   ZW(:,IK)=ZTSPLITG/PDZZ(:,IK)
END DO
!
IF (HPHASE=='L') ZC=CST%XCL
IF (HPHASE=='I') ZC=CST%XCI
!
! When pristine ice is 1-moment, nb concentration is parameterized following 
! McFarquhar and Heymsfield (1997) for columns as in ICE3
IF (KID==4 .AND. IMOMENTS==1) THEN
   IMOMENTS=2
   WHERE(PRS(:,:)>0) PCS(:,:)=1/(4*CST%XPI*900.) * PRS(:,:) * &
        MAX(0.05E6,-0.15319E6-0.021454E6*ALOG(PRHODREF(:,:)*PRS(:,:)))**3
END IF
!
! ################################
! Compute the sedimentation fluxes
! ################################
!
DO IN = 1 ,  LIMAP%NSPLITSED(KID)
  ! Computation only where enough ice, snow, graupel or hail
   GSEDIM(:,:) = .FALSE.
   GSEDIM(D%NIJB:D%NIJE,D%NKTB:D%NKTE) = PRS(D%NIJB:D%NIJE,D%NKTB:D%NKTE)>LIMAP%XRTMIN(KID)
   IF (IMOMENTS==2)  GSEDIM(:,:) = GSEDIM(:,:) .AND. PCS(:,:)>LIMAP%XCTMIN(KID)
   ISEDIM = COUNTJV( GSEDIM(:,:),I1(:),I3(:))
!
   IF( ISEDIM >= 1 ) THEN
!
      ALLOCATE(ZRHODREF(ISEDIM))
      ALLOCATE(ZPABST(ISEDIM))
      ALLOCATE(ZT(ISEDIM))
      ALLOCATE(ZRS(ISEDIM))
      ALLOCATE(ZCS(ISEDIM))
      ALLOCATE(ZLBDA(ISEDIM)) ; ZLBDA(:) = 1.E10
      ALLOCATE(ZCC(ISEDIM))   ; ZCC(:) = 1.0
      ALLOCATE(ZZW(ISEDIM))   ; ZZW(:) = 0.0
      ALLOCATE(ZZX(ISEDIM))   ; ZZX(:) = 0.0
      ALLOCATE(ZZY(ISEDIM))   ; ZZY(:) = 0.0
      !
      IF (OELEC) THEN
        ZWSEDQ(:,:) = 0.
        ALLOCATE(ZES(ISEDIM)) ; ZES(:) = 0.0
        ALLOCATE(ZQS(ISEDIM)) ; ZQS(:) = 0.0
        ALLOCATE(ZZQ(ISEDIM)) ; ZZQ(:) = 0.0
      END IF      
!
      DO IL = 1,ISEDIM
         ZRHODREF(IL) = PRHODREF(I1(IL),I3(IL))
         ZPABST(IL) = PPABST(I1(IL),I3(IL))
         ZT(IL) = PT(I1(IL),I3(IL))
         ZRS(IL) = PRS(I1(IL),I3(IL))
         IF (IMOMENTS==2) ZCS(IL) = PCS(I1(IL),I3(IL))
         IF (OELEC)       ZQS(IL) = PQS(I1(IL),I3(IL))
         IF (LIMAP%LCRYSTAL_SHAPE .AND. IMOMENTS == 2 .AND. KID == 4) &
           ZLBDA(IL) = PLBDAI_SHAPE(I1(IL),I3(IL))
      END DO
!
! Compute lambda
      IF (KID == 5 .AND. LIMAP%NMOM_S.EQ.1 .AND. LIMAP%LSNOW_T) THEN
         ZLBDA(:) = 1.E10
         WHERE(ZT(:)>263.15 .AND. ZRS(:)>LIMAP%XRTMIN(5))
            ZLBDA(:) = MAX(MIN(LIMAC%XLBDAS_MAX, 10**(14.554-0.0423*ZT(:))),LIMAC%XLBDAS_MIN)
         END WHERE
         WHERE(ZT(:)<=263.15 .AND. ZRS(:)>LIMAP%XRTMIN(5))
            ZLBDA(:) = MAX(MIN(LIMAC%XLBDAS_MAX, 10**(6.226-0.0106*ZT(:))),LIMAC%XLBDAS_MIN)
         END WHERE
         ZLBDA(:) = ZLBDA(:)*LIMAC%XTRANS_MP_GAMMAS
         ZZW(:) = LIMAP%XFSEDR(KID) * ZRHODREF(:)**(1.-LIMAP%XCEXVT)*ZRS(:)* &
              (1 + (LIMAC%XFVELOS/ZLBDA(:))**LIMAP%XALPHAS)**(-LIMAP%XNUS-(LIMAP%XD(KID)+LIMAC%XBS)/LIMAP%XALPHAS) &
              * ZLBDA(:)**(-LIMAP%XD(KID))
      ELSE IF (KID == 4 .AND. IMOMENTS==2 .AND. LIMAP%LCRYSTAL_SHAPE) THEN
         ZZY(:) = ZRHODREF(:)**(-LIMAP%XCEXVT) * ZLBDA(:)**(-LIMAC%XDI_SHAPE(KISHAPE))
         ZZW(:) = LIMAC%XFSEDRI_SHAPE(KISHAPE) * ZRS(:) * ZZY(:) * ZRHODREF(:)
         ZZX(:) = LIMAC%XFSEDCI_SHAPE(KISHAPE) * ZCS(:) * ZZY(:) * ZRHODREF(:)
      ELSE
         IF (IMOMENTS==1) ZLBDA(:) = LIMAP%XLB(KID) * ( ZRHODREF(:) * ZRS(:) )**LIMAP%XLBEX(KID)
         IF (IMOMENTS==2) ZLBDA(:) = ( LIMAP%XLB(KID)*ZCS(:) / ZRS(:) )**LIMAP%XLBEX(KID)
         ZZY(:) = ZRHODREF(:)**(-LIMAP%XCEXVT) * ZLBDA(:)**(-LIMAP%XD(KID))
         IF (LIMAP%LSNOW_T .AND. KID==5) &
              ZZY(:) = ZZY(:) * (1 + (LIMAC%XFVELOS/ZLBDA(:))**LIMAP%XALPHAS)**(-LIMAP%XNUS-(LIMAP%XD(KID)+LIMAC%XBS)/LIMAP%XALPHAS)
         ZZW(:) = LIMAP%XFSEDR(KID) * ZRS(:) * ZZY(:) * ZRHODREF(:)
      END IF ! Wurtz
!
      IF (KMOMENTS==2 .AND. .NOT.(LIMAP%LCRYSTAL_SHAPE .AND. KID==4)) &
         ZZX(:) = LIMAP%XFSEDC(KID) * ZCS(:) * ZZY(:) * ZRHODREF(:)

      IF (KID==2) THEN
         ! mean cloud droplet diameter
         ZCC(:) = LIMAW%XGCC/ZLBDA(:)
         ! correction factor for cloud droplet terminal fall speed
         ZCC(:) = 1.+1.26*6.6E-8*(101325./ZPABST(:))*(ZT(:)/293.15)/ZCC(:)
         ZZW(:) = ZCC(:) * ZZW(:)
         ZZX(:) = ZCC(:) * ZZX(:)
      END IF
!
! If the electrical scheme is activated, the electric field can impact the sedimentation
      ZBEARDCOEFF(:,:) = 1.0
      IF (OELEC .AND. ELECD%LSEDIM_BEARD) THEN
        ZLBDA3(:,:) = UNPACK( ZLBDA(:),MASK=GSEDIM(:,:),FIELD=0.0 )
        CALL ELEC_BEARD_EFFECT(D, CST, 'LIMA', KID, GSEDIM, PT, PRHODREF, PTHVREFZIKB, &
                               PRS, PQS, PEFIELDW, ZLBDA3, ZBEARDCOEFF, &
                               LIMAP=LIMAP, LIMAPC=LIMAC, LIMAPW=LIMAW, LIMAPM=LIMAM)
      END IF
!
      ZWSEDR(:,1:D%NKT) = UNPACK( ZZW(:),MASK=GSEDIM(:,:),FIELD=0.0 )
      ZWSEDR(:,D%NKTB:D%NKTE) = MIN( ZWSEDR(:,D%NKTB:D%NKTE) * ZBEARDCOEFF(:,D%NKTB:D%NKTE), &
                                       PRS(:,D%NKTB:D%NKTE) * PRHODREF(:,D%NKTB:D%NKTE) /      &
                                                                ZW(:,D%NKTB:D%NKTE) )
      IF (KMOMENTS==2) THEN
         ZWSEDC(:,1:D%NKT) = UNPACK( ZZX(:),MASK=GSEDIM(:,:),FIELD=0.0 )
         ZWSEDC(:,D%NKTB:D%NKTE) = MIN( ZWSEDC(:,D%NKTB:D%NKTE) * ZBEARDCOEFF(:,D%NKTB:D%NKTE), &
                                          PCS(:,D%NKTB:D%NKTE) * PRHODREF(:,D%NKTB:D%NKTE) /      &
                                                                   ZW(:,D%NKTB:D%NKTE) )
      END IF
!
! Sedimentation of electric charges
      IF (OELEC) THEN
        ! compute e of the q-D relationship
        IF (IMOMENTS == 2) THEN  ! 2-moment species
          CALL ELEC_COMPUTE_EX (KID, IMOMENTS, ISEDIM, 'LIMA', PTSTEP, ZRHODREF, LIMAP%XRTMIN(KID), &
                                ZRS, ZQS, ZES, PLBDX=ZLBDA, PCX=ZCS)
        ELSE                     ! 1-moment species
          CALL ELEC_COMPUTE_EX (KID, IMOMENTS, ISEDIM, 'LIMA', PTSTEP, ZRHODREF, LIMAP%XRTMIN(KID), &
                                ZRS, ZQS, ZES, PLBDX=ZLBDA)
        END IF
        !
        ! number concentration for 1-moment species
        ! for precipitating hydrometeors, N=C\lambda^x, except for snow if lsnow_t=t
        IF (IMOMENTS == 1) THEN
          IF (KID == 5) THEN
            ZCX = LIMAC%XCCS
            ZXX = LIMAC%XCXS
          ELSE IF (KID == 6) THEN
            ZCX = LIMAM%XCCG
            ZXX = LIMAM%XCXG
          ELSE IF (KID == 7) THEN
            ZCX = LIMAM%XCCH
            ZXX = LIMAM%XCXH
          END IF
          ZCS(:) = ZCX * ZLBDA(:)**ZXX
        END IF
        !
        ZZQ(:) = ZRHODREF(:)**(1.-LIMAP%XCEXVT) * ZES(:) * ZCS(:) * ELECP%XFQSED(KID) * ZLBDA(:)**(-ELECP%XDQ(KID))
        !
        ! correction for cloud droplet terminal fall speed
        IF (KID == 2)  ZZQ(:) = ZZQ(:) * ZCC(:)
        !
        ZWSEDQ(:,1:D%NKT) = UNPACK( ZZQ(:),MASK=GSEDIM(:,:),FIELD=0.0 )
        ZWSEDQ(:,1:D%NKT) = ZWSEDQ(:,1:D%NKT) * ZBEARDCOEFF(:,1:D%NKT)
        ZWSEDQ(:,D%NKTB:D%NKTE) = SIGN(MIN(ABS(ZWSEDQ(:,D%NKTB:D%NKTE)),                                                 &
                                             ABS(PQS(:,D%NKTB:D%NKTE)*PRHODREF(:,D%NKTB:D%NKTE)/ZW(:,D%NKTB:D%NKTE))), &
                                         ZWSEDQ(:,D%NKTB:D%NKTE))
      END IF      
      
      DO IK = D%NKTB , D%NKTE
         PRS(:,IK) = PRS(:,IK) + ZW(:,IK)*    &
              (ZWSEDR(:,IK+D%NKL)-ZWSEDR(:,IK))/PRHODREF(:,IK)
         PFPR(:,IK) = ZWSEDR(:,IK)
         IF (KMOMENTS==2) PCS(:,IK) = PCS(:,IK) + ZW(:,IK)*    &
              (ZWSEDC(:,IK+D%NKL)-ZWSEDC(:,IK))/PRHODREF(:,IK)
         ! Heat transport
         !PRT_SUM(:,IK-D%NKL) = PRT_SUM(:,IK-D%NKL) + ZW(:,IK-D%NKL)*ZWSEDR(:,IK)/PRHODREF(:,IK-D%NKL)
         !PRT_SUM(:,IK) = PRT_SUM(:,IK) - ZW(:,IK)*ZWSEDR(:,IK)/PRHODREF(:,IK)
         !PCPT(:,IK-D%NKL) = PCPT(:,IK-D%NKL) + ZC * (ZW(:,IK-D%NKL)*ZWSEDR(:,IK)/PRHODREF(:,IK-D%NKL))
         !PCPT(:,IK) = PCPT(:,IK) - ZC * (ZW(:,IK)*ZWSEDR(:,IK)/PRHODREF(:,IK))
         !ZWDT(:,IK) =(PRHODREF(:,IK+D%NKL)*(1.+PRT_SUM(:,IK))*PCPT(:,IK)*PT(:,IK) + &
         !     ZW(:,IK)*ZWSEDR(:,IK+1)*ZC*PT(:,IK+D%NKL)) / &
         !     (PRHODREF(:,IK+D%NKL)*(1.+PRT_SUM(:,IK))*PCPT(:,IK) + ZW(:,IK)*ZWSEDR(:,IK+D%NKL)*ZC)
         !ZWDT(:,IK) = ZWDT(:,IK) - PT(:,IK)
         IF (OELEC) PQS(:,IK) = PQS(:,IK) + ZW(:,IK) *    &
                                 (ZWSEDQ(:,IK+D%NKL) - ZWSEDQ(:,IK)) / PRHODREF(:,IK)
         
      END DO
      !
      DEALLOCATE(ZRHODREF)
      DEALLOCATE(ZPABST)
      DEALLOCATE(ZT)
      DEALLOCATE(ZRS)
      DEALLOCATE(ZCS)
      DEALLOCATE(ZCC)
      DEALLOCATE(ZLBDA)
      DEALLOCATE(ZZW)
      DEALLOCATE(ZZX)
      DEALLOCATE(ZZY)
      IF (ALLOCATED(ZQS))    DEALLOCATE(ZQS)
      IF (ALLOCATED(ZZQ))    DEALLOCATE(ZZQ)
      IF (ALLOCATED(ZES))    DEALLOCATE(ZES)      
      !      
      PINPR(:) = PINPR(:) + ZWSEDR(:,D%NKB)/CST%XRHOLW/LIMAP%NSPLITSED(KID)                          ! in m/s
      !PT(:,:) = PT(:,:) + ZWDT(:,:)
      
   END IF
END DO
!
PRS(:,:) = PRS(:,:) / PTSTEP
IF (KMOMENTS==2) PCS(:,:) = PCS(:,:) / PTSTEP
IF (OELEC) PQS(:,:) = PQS(:,:) / PTSTEP
!
IF (LHOOK) CALL DR_HOOK('LIMA_SEDIMENTATION', 1, ZHOOK_HANDLE)
END SUBROUTINE LIMA_SEDIMENTATION
END MODULE MODE_LIMA_SEDIMENTATION
