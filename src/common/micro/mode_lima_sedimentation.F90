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
!  C. Barthe   26/01/2024: add several ice crystal shapes
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_DIMPHYEX,         ONLY: DIMPHYEX_T
USE MODD_CST,              ONLY: CST_T
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
INTEGER :: IK, IL, IN, JK, JIJ            ! Loop index
INTEGER :: ISEDIM                         ! Case number of sedimentation
!
LOGICAL :: GSEDIM      ! Test where to compute the SED processes
REAL,    DIMENSION(D%NIJT, D%NKT) :: ZW,       & ! Work array
                                     ZWDT        ! Temperature change
REAL,    DIMENSION(D%NIJT,0:D%NKT+1) &
                           :: ZWSEDR,   & ! Sedimentation of MMR
                              ZWSEDC      ! Sedimentation of number conc.
!
REAL                       :: ZZW,      & ! Work array
                              ZZX,      & ! Work array
                              ZZY,      & ! Work array
                              ZLBDA,    & ! Slope parameter
                              ZCC         ! Cunningham corrective term for droplets fall speed
!
INTEGER , DIMENSION(D%NIJT*D%NKT) :: I1,I3 ! Indexes for PACK replacement
!
REAL    :: ZC                             ! Cpl or Cpi
INTEGER :: IMOMENTS
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
! Variables for cloud electricity
REAL                                :: ZZQ, &    ! Work array
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
REAL, DIMENSION(D%NIJT)       :: ZMAX_TSTEP1D ! Maximum CFL in column
REAL, DIMENSION(D%NIJT,D%NKT) :: ZMAX_TSTEP2D ! Maximum CFL in column
REAL, DIMENSION(D%NIJT)       :: ZREMAINT   ! Remaining time until the timestep end
LOGICAL :: ZANYREMAINT
REAL                            :: ZINVTSTEP
REAL                            :: ZMRCHANGE
REAL                            :: ZCCHANGE
REAL                            :: ZQCHANGE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('LIMA_SEDIMENTATION', 0, ZHOOK_HANDLE)
!
ZINVTSTEP=1./PTSTEP
IMOMENTS=KMOMENTS
PINPR(:) = 0.
PRS(D%NIJB:D%NIJE,:) = PRS(D%NIJB:D%NIJE,:) * PTSTEP
IF (KMOMENTS==2) PCS(D%NIJB:D%NIJE,:) = PCS(D%NIJB:D%NIJE,:) * PTSTEP
IF (OELEC)       PQS(D%NIJB:D%NIJE,:) = PQS(D%NIJB:D%NIJE,:) * PTSTEP
DO IK = D%NKTB , D%NKTE
   ZW(D%NIJB:D%NIJE,IK)=1./(PDZZ(D%NIJB:D%NIJE,IK)*PRHODREF(D%NIJB:D%NIJE,IK))
END DO
IF (HPHASE=='L') ZC=CST%XCL
IF (HPHASE=='I') ZC=CST%XCI
!
! Time splitting
!
ZREMAINT(:) = 0.
ZREMAINT(:) = PTSTEP
ZANYREMAINT = .TRUE.
DO WHILE (ZANYREMAINT)
!
   ZWDT(:,:)=0.
   PFPR(:,:) = 0.
   ZWSEDR(:,:) = 0.
   ZWSEDC(:,:) = 0.
   ZWSEDQ(:,:) = 0.
   ZLBDA3(:,:) = 0.
   !
      
   DO JK = D%NKTB,D%NKTE
      DO JIJ = D%NIJB,D%NIJE
!
! When pristine ice is 1-moment, nb concentration is parameterized following 
! McFarquhar and Heymsfield (1997) for columns as in ICE3
         IF (KID==4 .AND. IMOMENTS==1) THEN
            IMOMENTS=2
            IF (PRS(JIJ,JK)>LIMAP%XRTMIN(KID) .AND. ZREMAINT(JIJ)>0.) PCS(JIJ,JK)=1/(4*CST%XPI*900.) * PRS(JIJ,JK) * &
                 MAX(0.05E6,-0.15319E6-0.021454E6*LOG(PRHODREF(JIJ,JK)*PRS(JIJ,JK)))**3
         END IF
!
! ################################
! Compute the sedimentation fluxes
! ################################
!
         GSEDIM = PRS(JIJ,JK)>LIMAP%XRTMIN(KID) .AND. ZREMAINT(JIJ)>0.
         IF (IMOMENTS==2)  GSEDIM = GSEDIM .AND. PCS(JIJ,JK)>LIMAP%XCTMIN(KID)
!
         IF( GSEDIM ) THEN
!
            ZLBDA = 1.E10
            ZCC = 1.0
            ZZW = 0.0
            ZZX = 0.0
            ZZY = 0.0
            IF (OELEC) THEN
               ZES = 0.0
               ZZQ = 0.0
            END IF
! Compute lambda
            IF (KID == 5 .AND. LIMAP%NMOM_S.EQ.1 .AND. LIMAP%LSNOW_T) THEN
               ZLBDA= 1.E10
               IF (PT(JIJ,JK)>263.15 .AND. PRS(JIJ,JK)>LIMAP%XRTMIN(5)) THEN
                  ZLBDA = MAX(MIN(LIMAC%XLBDAS_MAX, 10**(14.554-0.0423*PT(JIJ,JK))),LIMAC%XLBDAS_MIN)
               ELSE IF(PT(JIJ,JK)<=263.15 .AND. PRS(JIJ,JK)>LIMAP%XRTMIN(5)) THEN
                  ZLBDA = MAX(MIN(LIMAC%XLBDAS_MAX, 10**(6.226-0.0106*PT(JIJ,JK))),LIMAC%XLBDAS_MIN)
               END IF
               ZLBDA = ZLBDA*LIMAC%XTRANS_MP_GAMMAS
               ZWSEDR(JIJ,JK) = LIMAP%XFSEDR(KID) * PRHODREF(JIJ,JK)**(1.-LIMAP%XCEXVT)*PRS(JIJ,JK)* &
                    (1 + (LIMAC%XFVELOS/ZLBDA)**LIMAP%XALPHAS)**(-LIMAP%XNUS-(LIMAP%XD(KID)+LIMAC%XBS)/LIMAP%XALPHAS) &
                    * ZLBDA**(-LIMAP%XD(KID))
            ELSE IF (KID == 4 .AND. IMOMENTS==2 .AND. LIMAP%LCRYSTAL_SHAPE) THEN
               ZZY = PRHODREF(JIJ,JK)**(-LIMAP%XCEXVT) * ZLBDA**(-LIMAC%XDI_SHAPE(KISHAPE))
               ZWSEDR(JIJ,JK) = LIMAC%XFSEDRI_SHAPE(KISHAPE) * PRS(JIJ,JK) * ZZY * PRHODREF(JIJ,JK)
               ZWSEDC(JIJ,JK) = LIMAC%XFSEDCI_SHAPE(KISHAPE) * PCS(JIJ,JK) * ZZY * PRHODREF(JIJ,JK)
               ZLBDA = PLBDAI_SHAPE(JIJ,JK)
            ELSE
               IF (IMOMENTS==1) ZLBDA = LIMAP%XLB(KID) * ( PRHODREF(JIJ,JK) * PRS(JIJ,JK) )**LIMAP%XLBEX(KID)
               IF (IMOMENTS==2) ZLBDA = ( LIMAP%XLB(KID)*PCS(JIJ,JK) / PRS(JIJ,JK) )**LIMAP%XLBEX(KID)
               ZZY = PRHODREF(JIJ,JK)**(-LIMAP%XCEXVT) * ZLBDA**(-LIMAP%XD(KID))
               IF (LIMAP%LSNOW_T .AND. KID==5) &
                    ZZY = ZZY * (1 + (LIMAC%XFVELOS/ZLBDA)**LIMAP%XALPHAS)**(-LIMAP%XNUS-(LIMAP%XD(KID)+LIMAC%XBS)/LIMAP%XALPHAS)
               ZWSEDR(JIJ,JK) = LIMAP%XFSEDR(KID) * PRS(JIJ,JK) * ZZY * PRHODREF(JIJ,JK)
            END IF ! Wurtz
            !
            IF (KMOMENTS==2 .AND. .NOT.(LIMAP%LCRYSTAL_SHAPE .AND. KID==4)) &
                 ZWSEDC(JIJ,JK) = LIMAP%XFSEDC(KID) * PCS(JIJ,JK) * ZZY * PRHODREF(JIJ,JK)
            
            IF (KID==2) THEN
               ! mean cloud droplet diameter
               ZCC = LIMAW%XGCC/ZLBDA
               ! correction factor for cloud droplet terminal fall speed
               ZCC = 1.+1.26*6.6E-8*(101325./PPABST(JIJ,JK))*(PT(JIJ,JK)/293.15)/ZCC
               ZWSEDR(JIJ,JK) = ZCC * ZWSEDR(JIJ,JK)
               ZWSEDC(JIJ,JK) = ZCC * ZWSEDC(JIJ,JK)
            END IF
            !
            !
!            IF (OELEC) THEN
!               ! compute e of the q-D relationship
!               IF (IMOMENTS == 2) THEN  ! 2-moment species
!                  CALL ELEC_COMPUTE_EX (KID, IMOMENTS, ISEDIM, 'LIMA', PTSTEP, ZRHODREF, LIMAP%XRTMIN(KID), &
!                       PRS(JIJ,JK), PQS(JIJ,JK), ZES, PLBDX=ZLBDA, PCX=PCS(JIJ,JK))
!               ELSE                     ! 1-moment species
!                  CALL ELEC_COMPUTE_EX (KID, IMOMENTS, ISEDIM, 'LIMA', PTSTEP, ZRHODREF, LIMAP%XRTMIN(KID), &
!                       PRS(JIJ,JK), PQS(JIJ,JK), ZES, PLBDX=ZLBDA)
!               END IF
!               !
!               IF (ELECD%LSEDIM_BEARD) ZLBDA3(JIJ,JK) = ZLBDA
!               !
!               IF (IMOMENTS == 1) THEN
!                  IF (KID == 5) THEN
!                     ZWSEDQ(JIJ,JK) = PRHODREF(JIJ,JK)**(1.-LIMAP%XCEXVT) * ZES * LIMAC%XCCS*ZLBDA**LIMAC%XCXS * ELECP%XFQSED(KID) * ZLBDA**(-ELECP%XDQ(KID))
!                  ELSE IF (KID == 6) THEN
!                     ZWSEDQ(JIJ,JK) = PRHODREF(JIJ,JK)**(1.-LIMAP%XCEXVT) * ZES * LIMAC%XCCG*ZLBDA**LIMAC%XCXG * ELECP%XFQSED(KID) * ZLBDA**(-ELECP%XDQ(KID))
!                  ELSE IF (KID == 7) THEN
!                     ZWSEDQ(JIJ,JK) = PRHODREF(JIJ,JK)**(1.-LIMAP%XCEXVT) * ZES * LIMAC%XCCH*ZLBDA**LIMAC%XCXH * ELECP%XFQSED(KID) * ZLBDA**(-ELECP%XDQ(KID))
!                  ELSE
!                     ZWSEDQ(JIJ,JK) = 0.
!                  END IF
!               ELSE
!                  ZWSEDQ(JIJ,JK) = PRHODREF(JIJ,JK)**(1.-LIMAP%XCEXVT) * ZES * PCS(JIJ,JK) * ELECP%XFQSED(KID) * ZLBDA**(-ELECP%XDQ(KID))
!               END IF
!               IF (KID == 2) ZWSEDQ(JIJ,JK) = ZWSEDQ(JIJ,JK) * ZCC
!            END IF
            !
         END IF
      END DO
   END DO
!
! If the electrical scheme is activated, the electric field can impact the sedimentation
!   ZBEARDCOEFF = 1.0
!   IF (OELEC .AND. ELECD%LSEDIM_BEARD) THEN
!      CALL ELEC_BEARD_EFFECT(D, CST, 'LIMA', KID, GSEDIM, PT, PRHODREF, PTHVREFZIKB, &
!           PRS, PQS, PEFIELDW, ZLBDA3, ZBEARDCOEFF, &
!           LIMAP=LIMAP, LIMAPC=LIMAC, LIMAPW=LIMAW, LIMAPM=LIMAM)
!      ZWSEDR(:,:)=ZWSEDR(:,:)*ZBEARDCOEFF(:,:)
!      ZWSEDC(:,:)=ZWSEDC(:,:)*ZBEARDCOEFF(:,:)
!      ZWSEDQ(:,:)=ZWSEDQ(:,:)*ZBEARDCOEFF(:,:)
!   END IF
   !
   ZMAX_TSTEP1D(:) = ZREMAINT(:)
   DO JK = D%NKTB,D%NKTE
      DO JIJ = D%NIJB,D%NIJE
         ZMAX_TSTEP2D(JIJ,JK) = ZREMAINT(JIJ)
      END DO
   END DO
   !
   DO JK = D%NKTB,D%NKTE
      DO JIJ = D%NIJB,D%NIJE
         IF(PRS(JIJ,JK)>LIMAP%XRTMIN(KID) .AND. ZWSEDR(JIJ, JK)>1.E-20 .AND. ZREMAINT(JIJ)>0.) THEN
            ZMAX_TSTEP2D(JIJ,JK) = 0.8 * PRHODREF(JIJ, JK) * PRS(JIJ, JK) * PDZZ(JIJ, JK) / ZWSEDR(JIJ, JK)
         ENDIF
      ENDDO
   ENDDO
   !
   DO JK = D%NKTB,D%NKTE
      DO JIJ = D%NIJB,D%NIJE
         IF(PRS(JIJ,JK)>LIMAP%XRTMIN(KID) .AND. ZWSEDR(JIJ, JK)>1.E-20 .AND. ZREMAINT(JIJ)>0.) THEN
            ZMAX_TSTEP1D(JIJ) = MIN(ZMAX_TSTEP1D(JIJ), ZMAX_TSTEP2D(JIJ,JK))
         ENDIF
      ENDDO
   ENDDO
   !
   DO JIJ = D%NIJB, D%NIJE
      ZREMAINT(JIJ) = ZREMAINT(JIJ) - ZMAX_TSTEP1D(JIJ)
      PINPR(JIJ) = PINPR(JIJ) + ZWSEDR(JIJ,D%NKB) / CST%XRHOLW * (ZMAX_TSTEP1D(JIJ) * ZINVTSTEP)
   ENDDO
   !
   DO JK = D%NKTB, D%NKTE
      DO JIJ = D%NIJB, D%NIJE
         ZMRCHANGE = ZMAX_TSTEP1D(JIJ) * ZW(JIJ,JK)*(ZWSEDR(JIJ,JK+D%NKL)-ZWSEDR(JIJ,JK))
         PRS(JIJ,JK) = PRS(JIJ,JK) + ZMRCHANGE
         PFPR(JIJ,JK) = PFPR(JIJ,JK) + ZWSEDR(JIJ,JK) * (ZMAX_TSTEP1D(JIJ) * ZINVTSTEP)
         IF (KMOMENTS==2) THEN
            ZCCHANGE = ZMAX_TSTEP1D(JIJ) * ZW(JIJ,JK)*(ZWSEDC(JIJ,JK+D%NKL)-ZWSEDC(JIJ,JK))
            PCS(JIJ,JK) = PCS(JIJ,JK) + ZCCHANGE
         END IF
         IF (OELEC) THEN
            ZQCHANGE = ZMAX_TSTEP1D(JIJ) * ZW(JIJ,JK) * (ZWSEDQ(JIJ,JK+D%NKL) - ZWSEDQ(JIJ,JK))
            PQS(JIJ,JK) = PQS(JIJ,JK) + ZQCHANGE
         ENDIF
      ENDDO
   ENDDO
   !   
   ZANYREMAINT = .FALSE.
   DO JIJ=D%NIJB,D%NIJE
      IF(ZREMAINT(JIJ)>0.) THEN
         ZANYREMAINT = .TRUE.
      END IF
   END DO
   !
END DO
!
PRS(D%NIJB:D%NIJE,:) = PRS(D%NIJB:D%NIJE,:) / PTSTEP
IF (KMOMENTS==2) PCS(D%NIJB:D%NIJE,:) = PCS(D%NIJB:D%NIJE,:) / PTSTEP
IF (OELEC) PQS(D%NIJB:D%NIJE,:) = PQS(D%NIJB:D%NIJE,:) / PTSTEP
!
IF (LHOOK) CALL DR_HOOK('LIMA_SEDIMENTATION', 1, ZHOOK_HANDLE)
END SUBROUTINE LIMA_SEDIMENTATION
END MODULE MODE_LIMA_SEDIMENTATION
