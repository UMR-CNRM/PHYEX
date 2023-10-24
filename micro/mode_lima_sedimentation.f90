!MNH_LIC Copyright 2018-2020 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_LIMA_SEDIMENTATION
  IMPLICIT NONE
CONTAINS
!     ######################################################################
  SUBROUTINE LIMA_SEDIMENTATION (D, CST, &
                                 HPHASE, KMOMENTS, KID, KSPLITG, PTSTEP, OELEC, &
                                 PDZZ, PRHODREF, &
                                 PPABST, PT, PRT_SUM, PCPT, PRS, PCS, PINPR, PFPR, &
                                 PEFIELDW, PQS)
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
USE MODD_DIMPHYEX,         ONLY: DIMPHYEX_t
USE MODD_CST,              ONLY: CST_t
USE MODD_ELEC_DESCR,       ONLY: LSEDIM_BEARD
USE MODD_ELEC_PARAM,       ONLY: XFQSED, XDQ
USE MODD_PARAM_LIMA,       ONLY: XCEXVT, XRTMIN, XCTMIN, NSPLITSED,           &
                                 XLB, XLBEX, XD, XFSEDR, XFSEDC,              &
                                 XALPHAC, XNUC, XALPHAS, XNUS, LSNOW_T,       &
                                 NMOM_S
USE MODD_PARAM_LIMA_COLD,  ONLY: XLBEXI, XLBI, XDI, XLBDAS_MAX, XBS, XEXSEDS, &
                                 XLBDAS_MIN, XTRANS_MP_GAMMAS, XFVELOS,       &
                                 XCCS, XCXS
USE MODD_PARAM_LIMA_MIXED, ONLY: XCCG, XCXG, XCCH, XCXH

use mode_tools,            only: Countjv

USE MODI_GAMMA,            ONLY: GAMMA_X0D
USE MODI_ELEC_COMPUTE_EX
USE MODE_ELEC_BEARD_EFFECT, ONLY: ELEC_BEARD_EFFECT
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(DIMPHYEX_t),         INTENT(IN)    :: D
TYPE(CST_t),              INTENT(IN)    :: CST
CHARACTER(1),             INTENT(IN)    :: HPHASE    ! Liquid or solid hydrometeors
INTEGER,                  INTENT(IN)    :: KMOMENTS  ! Number of moments 
INTEGER,                  INTENT(IN)    :: KID       ! Hydrometeor ID
INTEGER,                  INTENT(IN)    :: KSPLITG   !  
REAL,                     INTENT(IN)    :: PTSTEP    ! Time step  
LOGICAL,                  INTENT(IN)    :: OELEC     ! if true, cloud electrification is activated
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PDZZ      ! Height (z)
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF  ! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABST    ! abs. pressure at time t
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PT        ! Temperature
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRT_SUM   ! total water mixing ratio
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCPT      ! Cp
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRS       ! m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCS       ! C. source
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PINPR     ! Instant precip rate
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PFPR      ! Precip. fluxes in altitude
REAL, DIMENSION(:,:,:),   INTENT(IN),    OPTIONAL :: PEFIELDW  ! Vertical component of the electric field
REAL, DIMENSION(:,:,:),   INTENT(INOUT), OPTIONAL :: PQS ! Elec. charge density source
!
!*       0.2   Declarations of local variables :
!
INTEGER :: JK, JL, JN                     ! Loop index
INTEGER :: ISEDIM                         ! Case number of sedimentation
!
LOGICAL, DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3)) &
                           :: GSEDIM      ! Test where to compute the SED processes
REAL,    DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3)) &
                           :: ZW,       & ! Work array
                              ZWDT        ! Temperature change
REAL,    DIMENSION(D%NIT,D%NJT,0:D%NKT+1) &
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
INTEGER , DIMENSION(SIZE(PRHODREF)) :: I1,I2,I3 ! Indexes for PACK replacement
!
REAL    :: ZTSPLITG                       ! Small time step for rain sedimentation
REAL    :: ZC                             ! Cpl or Cpi
INTEGER :: IMOMENTS
!
! Variables for cloud electricity
REAL :: ZCX, ZXX  ! C and x parameters for N-lambda relationship
REAL, DIMENSION(:),     ALLOCATABLE :: ZQS, &    ! Electric charge density source
                                       ZZQ, &    ! Work array
                                       ZES       ! e in q-D relationship
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZWSEDQ    ! Sedimentation of electric charge density
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZLBDA3
REAL, DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3)) :: ZBEARDCOEFF ! effect of
                                                 ! electrical forces on terminal fall speed
!
!-------------------------------------------------------------------------------
!
IMOMENTS=KMOMENTS
!
! Time splitting
!
ZTSPLITG= PTSTEP / REAL(NSPLITSED(KID))
!
ZWDT=0.
PINPR(:,:) = 0.
ZWSEDR(:,:,:) = 0.
ZWSEDC(:,:,:) = 0.
!
PRS(:,:,:) = PRS(:,:,:) * PTSTEP
IF (KMOMENTS==2) PCS(:,:,:) = PCS(:,:,:) * PTSTEP
IF (OELEC)       PQS(:,:,:) = PQS(:,:,:) * PTSTEP
DO JK = D%NKTB , D%NKTE
   ZW(:,:,JK)=ZTSPLITG/PDZZ(:,:,JK)
END DO
!
IF (HPHASE=='L') ZC=CST%XCL
IF (HPHASE=='I') ZC=CST%XCI
!
! When pristine ice is 1-moment, nb concentration is parameterized following 
! McFarquhar and Heymsfield (1997) for columns as in ICE3
IF (KID==4 .AND. IMOMENTS==1) THEN
   IMOMENTS=2
   WHERE(PRS(:,:,:)>0) PCS(:,:,:)=1/(4*CST%XPI*900.) * PRS(:,:,:) * &
        MAX(0.05E6,-0.15319E6-0.021454E6*ALOG(PRHODREF(:,:,:)*PRS(:,:,:)))**3
END IF
!
! ################################
! Compute the sedimentation fluxes
! ################################
!
DO JN = 1 ,  NSPLITSED(KID)
  ! Computation only where enough ice, snow, graupel or hail
   GSEDIM(:,:,:) = .FALSE.
   GSEDIM(D%NIB:D%NIE,D%NJB:D%NJE,D%NKTB:D%NKTE) = PRS(D%NIB:D%NIE,D%NJB:D%NJE,D%NKTB:D%NKTE)>XRTMIN(KID)
   IF (IMOMENTS==2)  GSEDIM(:,:,:) = GSEDIM(:,:,:) .AND. PCS(:,:,:)>XCTMIN(KID)
   ISEDIM = COUNTJV( GSEDIM(:,:,:),I1(:),I2(:),I3(:))
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
        ALLOCATE(ZWSEDQ(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3))) ; ZWSEDQ(:,:,:) = 0.
        ALLOCATE(ZES(ISEDIM)) ; ZES(:) = 0.0
        ALLOCATE(ZQS(ISEDIM)) ; ZQS(:) = 0.0
        ALLOCATE(ZZQ(ISEDIM)) ; ZZQ(:) = 0.0
      END IF      
!
      DO JL = 1,ISEDIM
         ZRHODREF(JL) = PRHODREF(I1(JL),I2(JL),I3(JL))
         ZPABST(JL) = PPABST(I1(JL),I2(JL),I3(JL))
         ZT(JL) = PT(I1(JL),I2(JL),I3(JL))
         ZRS(JL) = PRS(I1(JL),I2(JL),I3(JL))
         IF (IMOMENTS==2) ZCS(JL) = PCS(I1(JL),I2(JL),I3(JL))
         IF (OELEC)       ZQS(JL) = PQS(I1(JL),I2(JL),I3(JL))
      END DO
!
! Compute lambda
      IF (KID == 5 .AND. NMOM_S.EQ.1 .AND. LSNOW_T) THEN
         ZLBDA(:) = 1.E10
         WHERE(ZT(:)>263.15 .AND. ZRS(:)>XRTMIN(5))
            ZLBDA(:) = MAX(MIN(XLBDAS_MAX, 10**(14.554-0.0423*ZT(:))),XLBDAS_MIN)
         END WHERE
         WHERE(ZT(:)<=263.15 .AND. ZRS(:)>XRTMIN(5))
            ZLBDA(:) = MAX(MIN(XLBDAS_MAX, 10**(6.226-0.0106*ZT(:))),XLBDAS_MIN)
         END WHERE
         ZLBDA(:) = ZLBDA(:)*XTRANS_MP_GAMMAS
         ZZW(:) = XFSEDR(KID) * ZRHODREF(:)**(1.-XCEXVT)*ZRS(:)* &
              (1 + (XFVELOS/ZLBDA(:))**XALPHAS)**(-XNUS-(XD(KID)+XBS)/XALPHAS) * ZLBDA(:)**(-XD(KID))
      ELSE
         IF (IMOMENTS==1) ZLBDA(:) = XLB(KID) * ( ZRHODREF(:) * ZRS(:) )**XLBEX(KID)
         IF (IMOMENTS==2) ZLBDA(:) = ( XLB(KID)*ZCS(:) / ZRS(:) )**XLBEX(KID)
         ZZY(:) = ZRHODREF(:)**(-XCEXVT) * ZLBDA(:)**(-XD(KID))
         IF (LSNOW_T .AND. KID==5) &
              ZZY(:) = ZZY(:) * (1 + (XFVELOS/ZLBDA(:))**XALPHAS)**(-XNUS-(XD(KID)+XBS)/XALPHAS)
         ZZW(:) = XFSEDR(KID) * ZRS(:) * ZZY(:) * ZRHODREF(:)
      END IF ! Wurtz
!
      IF (KMOMENTS==2) ZZX(:) = XFSEDC(KID) * ZCS(:) * ZZY(:) * ZRHODREF(:)

      IF (KID==2) THEN
         ! mean cloud droplet diameter
         ZCC(:) = 0.5*GAMMA_X0D(XNUC+1./XALPHAC)/(GAMMA_X0D(XNUC)*ZLBDA(:))
         ! correction factor for cloud droplet terminal fall speed
         ZCC(:) = 1.+1.26*6.6E-8*(101325./ZPABST(:))*(ZT(:)/293.15)/ZCC(:)
         ZZW(:) = ZCC(:) * ZZW(:)
         ZZX(:) = ZCC(:) * ZZX(:)
      END IF
!
! If the electrical scheme is activated, the electric field can impact the sedimentation
      ZBEARDCOEFF(:,:,:) = 1.0
      IF (OELEC .AND. LSEDIM_BEARD) THEN
        ALLOCATE(ZLBDA3(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3)))
        ZLBDA3(:,:,:) = UNPACK( ZLBDA(:),MASK=GSEDIM(:,:,:),FIELD=0.0 )
        CALL ELEC_BEARD_EFFECT(D, KID, GSEDIM, PT, PRHODREF, &
                               PRS, PQS, PEFIELDW, ZLBDA3, ZBEARDCOEFF)
        DEALLOCATE(ZLBDA3)
      END IF
!
      ZWSEDR(:,:,1:D%NKT) = UNPACK( ZZW(:),MASK=GSEDIM(:,:,:),FIELD=0.0 )
      ZWSEDR(:,:,D%NKTB:D%NKTE) = MIN( ZWSEDR(:,:,D%NKTB:D%NKTE) * ZBEARDCOEFF(:,:,D%NKTB:D%NKTE), &
                                       PRS(:,:,D%NKTB:D%NKTE) * PRHODREF(:,:,D%NKTB:D%NKTE) /      &
                                                                ZW(:,:,D%NKTB:D%NKTE) )
      IF (KMOMENTS==2) THEN
         ZWSEDC(:,:,1:D%NKT) = UNPACK( ZZX(:),MASK=GSEDIM(:,:,:),FIELD=0.0 )
         ZWSEDC(:,:,D%NKTB:D%NKTE) = MIN( ZWSEDC(:,:,D%NKTB:D%NKTE) * ZBEARDCOEFF(:,:,D%NKTB:D%NKTE), &
                                          PCS(:,:,D%NKTB:D%NKTE) * PRHODREF(:,:,D%NKTB:D%NKTE) /      &
                                                                   ZW(:,:,D%NKTB:D%NKTE) )
      END IF
!
! Sedimentation of electric charges
      IF (OELEC) THEN
        ! compute e of the q-D relationship
        IF (IMOMENTS == 2) THEN  ! 2-moment species
          CALL ELEC_COMPUTE_EX (KID, IMOMENTS, ISEDIM, PTSTEP, ZRHODREF, XRTMIN(KID), &
                                ZRS, ZQS, ZES, PLBDX=ZLBDA, PCX=ZCS)
        ELSE                     ! 1-moment species
          CALL ELEC_COMPUTE_EX (KID, IMOMENTS, ISEDIM, PTSTEP, ZRHODREF, XRTMIN(KID), &
                                ZRS, ZQS, ZES, PLBDX=ZLBDA)
        END IF
        !
        ! number concentration for 1-moment species
        ! for precipitating hydrometeors, N=C\lambda^x, except for snow if lsnow_t=t
        IF (IMOMENTS == 1) THEN
          IF (KID == 5) THEN
            ZCX = XCCS
            ZXX = XCXS
          ELSE IF (KID == 6) THEN
            ZCX = XCCG
            ZXX = XCXG
          ELSE IF (KID == 7) THEN
            ZCX = XCCH
            ZXX = XCXH
          END IF
          ZCS(:) = ZCX * ZLBDA(:)**ZXX
        END IF
        !
        ZZQ(:) = ZRHODREF(:)**(1.-XCEXVT) * ZES(:) * ZCS(:) * XFQSED(KID) * ZLBDA(:)**(-XDQ(KID))
        !
        ! correction for cloud droplet terminal fall speed
        IF (KID == 2)  ZZQ(:) = ZZQ(:) * ZCC(:)
        !
        ZWSEDQ(:,:,1:D%NKT) = UNPACK( ZZQ(:),MASK=GSEDIM(:,:,:),FIELD=0.0 )
        ZWSEDQ(:,:,1:D%NKT) = ZWSEDQ(:,:,1:D%NKT) * ZBEARDCOEFF(:,:,1:D%NKT)
        ZWSEDQ(:,:,D%NKTB:D%NKTE) = SIGN(MIN(ABS(ZWSEDQ(:,:,D%NKTB:D%NKTE)),                                                 &
                                             ABS(PQS(:,:,D%NKTB:D%NKTE)*PRHODREF(:,:,D%NKTB:D%NKTE)/ZW(:,:,D%NKTB:D%NKTE))), &
                                         ZWSEDQ(:,:,D%NKTB:D%NKTE))
      END IF      
      
      DO JK = D%NKTB , D%NKTE
         PRS(:,:,JK) = PRS(:,:,JK) + ZW(:,:,JK)*    &
              (ZWSEDR(:,:,JK+D%NKL)-ZWSEDR(:,:,JK))/PRHODREF(:,:,JK)
         PFPR(:,:,JK) = ZWSEDR(:,:,JK)
         IF (KMOMENTS==2) PCS(:,:,JK) = PCS(:,:,JK) + ZW(:,:,JK)*    &
              (ZWSEDC(:,:,JK+D%NKL)-ZWSEDC(:,:,JK))/PRHODREF(:,:,JK)
         ! Heat transport
         !PRT_SUM(:,:,JK-D%NKL) = PRT_SUM(:,:,JK-D%NKL) + ZW(:,:,JK-D%NKL)*ZWSEDR(:,:,JK)/PRHODREF(:,:,JK-D%NKL)
         !PRT_SUM(:,:,JK) = PRT_SUM(:,:,JK) - ZW(:,:,JK)*ZWSEDR(:,:,JK)/PRHODREF(:,:,JK)
         !PCPT(:,:,JK-D%NKL) = PCPT(:,:,JK-D%NKL) + ZC * (ZW(:,:,JK-D%NKL)*ZWSEDR(:,:,JK)/PRHODREF(:,:,JK-D%NKL))
         !PCPT(:,:,JK) = PCPT(:,:,JK) - ZC * (ZW(:,:,JK)*ZWSEDR(:,:,JK)/PRHODREF(:,:,JK))
         !ZWDT(:,:,JK) =(PRHODREF(:,:,JK+D%NKL)*(1.+PRT_SUM(:,:,JK))*PCPT(:,:,JK)*PT(:,:,JK) + &
         !     ZW(:,:,JK)*ZWSEDR(:,:,JK+1)*ZC*PT(:,:,JK+D%NKL)) / &
         !     (PRHODREF(:,:,JK+D%NKL)*(1.+PRT_SUM(:,:,JK))*PCPT(:,:,JK) + ZW(:,:,JK)*ZWSEDR(:,:,JK+D%NKL)*ZC)
         !ZWDT(:,:,JK) = ZWDT(:,:,JK) - PT(:,:,JK)
         IF (OELEC) PQS(:,:,JK) = PQS(:,:,JK) + ZW(:,:,JK) *    &
                                 (ZWSEDQ(:,:,JK+D%NKL) - ZWSEDQ(:,:,JK)) / PRHODREF(:,:,JK)
         
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
      IF (ALLOCATED(ZWSEDQ)) DEALLOCATE(ZWSEDQ)
      IF (ALLOCATED(ZQS))    DEALLOCATE(ZQS)
      IF (ALLOCATED(ZZQ))    DEALLOCATE(ZZQ)
      IF (ALLOCATED(ZES))    DEALLOCATE(ZES)      
      !      
      PINPR(:,:) = PINPR(:,:) + ZWSEDR(:,:,D%NKB)/CST%XRHOLW/NSPLITSED(KID)                          ! in m/s
      !PT(:,:,:) = PT(:,:,:) + ZWDT(:,:,:)
      
   END IF
END DO
!
PRS(:,:,:) = PRS(:,:,:) / PTSTEP
IF (KMOMENTS==2) PCS(:,:,:) = PCS(:,:,:) / PTSTEP
IF (OELEC) PQS(:,:,:) = PQS(:,:,:) / PTSTEP
!
END SUBROUTINE LIMA_SEDIMENTATION
END MODULE MODE_LIMA_SEDIMENTATION
