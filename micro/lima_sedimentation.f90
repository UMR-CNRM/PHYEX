!MNH_LIC Copyright 2018-2020 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!      ###################################
       MODULE MODI_LIMA_SEDIMENTATION
!      ###################################
!
INTERFACE
      SUBROUTINE LIMA_SEDIMENTATION (KIB, KIE, KIT, KJB, KJE, KJT, KKB, KKE, KKTB, KKTE, KKT, KKL, &
                                     HPHASE, KMOMENTS, KID, KSPLITG, PTSTEP, PDZZ, PRHODREF,       &
                                     PPABST, PT, PRT_SUM, PCPT, PRS, PCS, PINPR )
!
INTEGER,                  INTENT(IN)    :: KIB, KIE, KIT, KJB, KJE, KJT, KKB, KKE, KKTB, KKTE, KKT, KKL
CHARACTER(1),             INTENT(IN)    :: HPHASE    ! Liquid or solid hydrometeors
INTEGER,                  INTENT(IN)    :: KMOMENTS  ! Number of moments 
INTEGER,                  INTENT(IN)    :: KID       ! Hydrometeor ID
INTEGER,                  INTENT(IN)    :: KSPLITG   !  
REAL,                     INTENT(IN)    :: PTSTEP    ! Time step          
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PDZZ      ! Height (z)
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF  ! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABST    ! abs. pressure at time t
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PT        ! Temperature
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRT_SUM   ! total water mixing ratio
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCPT      ! Cp
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRS       ! m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCS       ! C. source
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PINPR     ! Instant precip rate
!
END SUBROUTINE LIMA_SEDIMENTATION
END INTERFACE
END MODULE MODI_LIMA_SEDIMENTATION
!
!
!     ######################################################################
      SUBROUTINE LIMA_SEDIMENTATION (KIB, KIE, KIT, KJB, KJE, KJT, KKB, KKE, KKTB, KKTE, KKT, KKL, &
                                     HPHASE, KMOMENTS, KID, KSPLITG, PTSTEP, PDZZ, PRHODREF,       &
                                     PPABST, PT, PRT_SUM, PCPT, PRS, PCS, PINPR )
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
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST,              ONLY: XRHOLW, XCL, XCI, XPI
USE MODD_PARAMETERS,       ONLY: JPHEXT, JPVEXT
USE MODD_PARAM_LIMA,       ONLY: XCEXVT, XRTMIN, XCTMIN, NSPLITSED,           &
                                 XLB, XLBEX, XD, XFSEDR, XFSEDC,              &
                                 XALPHAC, XNUC, XALPHAS, XNUS, LSNOW_T,       &
                                 NMOM_S
USE MODD_PARAM_LIMA_COLD,  ONLY: XLBEXI, XLBI, XDI, XLBDAS_MAX, XBS, XEXSEDS, &
                                 XLBDAS_MIN, XTRANS_MP_GAMMAS, XFVELOS

use mode_tools,            only: Countjv

USE MODI_GAMMA,            ONLY: GAMMA_X0D
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER,                  INTENT(IN)    :: KIB, KIE, KIT, KJB, KJE, KJT, KKB, KKE, KKTB, KKTE, KKT, KKL
CHARACTER(1),             INTENT(IN)    :: HPHASE    ! Liquid or solid hydrometeors
INTEGER,                  INTENT(IN)    :: KMOMENTS  ! Number of moments 
INTEGER,                  INTENT(IN)    :: KID       ! Hydrometeor ID
INTEGER,                  INTENT(IN)    :: KSPLITG   !  
REAL,                     INTENT(IN)    :: PTSTEP    ! Time step          
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PDZZ      ! Height (z)
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF  ! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABST    ! abs. pressure at time t
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PT        ! Temperature
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRT_SUM   ! total water mixing ratio
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCPT      ! Cp
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRS       ! m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCS       ! C. source
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PINPR     ! Instant precip rate
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
                              ZWSEDR,   & ! Sedimentation of MMR
                              ZWSEDC,   & ! Sedimentation of number conc.
                              ZWDT        ! Temperature change
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
INTEGER :: ZMOMENTS
!
!
!-------------------------------------------------------------------------------
!
ZMOMENTS=KMOMENTS
!
! Time splitting
!
ZTSPLITG= PTSTEP / REAL(NSPLITSED(KID))
!
ZWDT=0.
PINPR(:,:) = 0.
!
PRS(:,:,:) = PRS(:,:,:) * PTSTEP
IF (KMOMENTS==2) PCS(:,:,:) = PCS(:,:,:) * PTSTEP
DO JK = KKTB , KKTE
   ZW(:,:,JK)=ZTSPLITG/PDZZ(:,:,JK)
END DO
!
IF (HPHASE=='L') ZC=XCL
IF (HPHASE=='I') ZC=XCI
!
IF (KID==4 .AND. ZMOMENTS==1) THEN
   ZMOMENTS=2
   WHERE(PRS(:,:,:)>0) PCS(:,:,:)=1/(4*XPI*900.) * PRS(:,:,:) * &
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
   GSEDIM(KIB:KIE,KJB:KJE,KKTB:KKTE) = PRS(KIB:KIE,KJB:KJE,KKTB:KKTE)>XRTMIN(KID)
   IF (ZMOMENTS==2)  GSEDIM(:,:,:) = GSEDIM(:,:,:) .AND. PCS(:,:,:)>XCTMIN(KID)
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
      DO JL = 1,ISEDIM
         ZRHODREF(JL) = PRHODREF(I1(JL),I2(JL),I3(JL))
         ZPABST(JL) = PPABST(I1(JL),I2(JL),I3(JL))
         ZT(JL) = PT(I1(JL),I2(JL),I3(JL))
         ZRS(JL) = PRS(I1(JL),I2(JL),I3(JL))
         IF (ZMOMENTS==2) ZCS(JL) = PCS(I1(JL),I2(JL),I3(JL))
      END DO
!
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
         IF (ZMOMENTS==1) ZLBDA(:) = XLB(KID) * ( ZRHODREF(:) * ZRS(:) )**XLBEX(KID)
         IF (ZMOMENTS==2) ZLBDA(:) = ( XLB(KID)*ZCS(:) / ZRS(:) )**XLBEX(KID)
         ZZY(:) = ZRHODREF(:)**(-XCEXVT) * ZLBDA(:)**(-XD(KID))
         IF (LSNOW_T .AND. KID==5) &
              ZZY(:) = ZZY(:) * (1 + (XFVELOS/ZLBDA(:))**XALPHAS)**(-XNUS-(XD(KID)+XBS)/XALPHAS)
         ZZW(:) = XFSEDR(KID) * ZRS(:) * ZZY(:) * ZRHODREF(:)
      END IF ! Wurtz
!
      IF (KMOMENTS==2) ZZX(:) = XFSEDC(KID) * ZCS(:) * ZZY(:) * ZRHODREF(:)

      IF (KID==2) THEN
         ZCC(:) = 0.5*GAMMA_X0D(XNUC+1./XALPHAC)/(GAMMA_X0D(XNUC)*ZLBDA(:))
         ZCC(:) = 1.+1.26*6.6E-8*(101325./ZPABST(:))*(ZT(:)/293.15)/ZCC(:)
         ZZW(:) = ZCC(:) * ZZW(:)
         ZZX(:) = ZCC(:) * ZZX(:)
      END IF

      ZWSEDR(:,:,:) = UNPACK( ZZW(:),MASK=GSEDIM(:,:,:),FIELD=0.0 )
      ZWSEDR(:,:,KKTB:KKTE) = MIN( ZWSEDR(:,:,KKTB:KKTE), PRS(:,:,KKTB:KKTE) * PRHODREF(:,:,KKTB:KKTE) / ZW(:,:,KKTB:KKTE) )
      IF (KMOMENTS==2) THEN
         ZWSEDC(:,:,:) = UNPACK( ZZX(:),MASK=GSEDIM(:,:,:),FIELD=0.0 )
         ZWSEDC(:,:,KKTB:KKTE) = MIN( ZWSEDC(:,:,KKTB:KKTE), PCS(:,:,KKTB:KKTE) * PRHODREF(:,:,KKTB:KKTE) / ZW(:,:,KKTB:KKTE) )
      END IF
      
      DO JK = KKTB , KKTE
         PRS(:,:,JK) = PRS(:,:,JK) + ZW(:,:,JK)*    &
              (ZWSEDR(:,:,JK+KKL)-ZWSEDR(:,:,JK))/PRHODREF(:,:,JK)
         IF (KMOMENTS==2) PCS(:,:,JK) = PCS(:,:,JK) + ZW(:,:,JK)*    &
              (ZWSEDC(:,:,JK+KKL)-ZWSEDC(:,:,JK))/PRHODREF(:,:,JK)
         ! Heat transport
         !PRT_SUM(:,:,JK-KKL) = PRT_SUM(:,:,JK-KKL) + ZW(:,:,JK-KKL)*ZWSEDR(:,:,JK)/PRHODREF(:,:,JK-KKL)
         !PRT_SUM(:,:,JK) = PRT_SUM(:,:,JK) - ZW(:,:,JK)*ZWSEDR(:,:,JK)/PRHODREF(:,:,JK)
         !PCPT(:,:,JK-KKL) = PCPT(:,:,JK-KKL) + ZC * (ZW(:,:,JK-KKL)*ZWSEDR(:,:,JK)/PRHODREF(:,:,JK-KKL))
         !PCPT(:,:,JK) = PCPT(:,:,JK) - ZC * (ZW(:,:,JK)*ZWSEDR(:,:,JK)/PRHODREF(:,:,JK))
         !ZWDT(:,:,JK) =(PRHODREF(:,:,JK+KKL)*(1.+PRT_SUM(:,:,JK))*PCPT(:,:,JK)*PT(:,:,JK) + &
         !     ZW(:,:,JK)*ZWSEDR(:,:,JK+1)*ZC*PT(:,:,JK+KKL)) / &
         !     (PRHODREF(:,:,JK+KKL)*(1.+PRT_SUM(:,:,JK))*PCPT(:,:,JK) + ZW(:,:,JK)*ZWSEDR(:,:,JK+KKL)*ZC)
         !ZWDT(:,:,JK) = ZWDT(:,:,JK) - PT(:,:,JK)
      END DO
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
      !      
      PINPR(:,:) = PINPR(:,:) + ZWSEDR(:,:,KKB)/XRHOLW/NSPLITSED(KID)                          ! in m/s
      !PT(:,:,:) = PT(:,:,:) + ZWDT(:,:,:)
      
   END IF
END DO
!
PRS(:,:,:) = PRS(:,:,:) / PTSTEP
IF (KMOMENTS==2) PCS(:,:,:) = PCS(:,:,:) / PTSTEP
!
END SUBROUTINE LIMA_SEDIMENTATION
!
!-------------------------------------------------------------------------------
