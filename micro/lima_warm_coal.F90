!      ##########################
       MODULE MODI_LIMA_WARM_COAL
!      ##########################
!
INTERFACE
      SUBROUTINE LIMA_WARM_COAL (PTSTEP, KMI, HFMFILE, HLUOUT, OCLOSE_OUT,   &
                                 PRHODREF, ZWLBDC3, ZWLBDC, ZWLBDR3, ZWLBDR, &
                                 PRCT, PRRT, PCCT, PCRT,                     &
                                 PRCS, PRRS, PCCS, PCRS,                     &
                                 PRHODJ, &
                            YDDDH, YDLDDH, YDMDDH                              )
USE DDH_MIX, ONLY  : TYP_DDH
USE YOMLDDH, ONLY  : TLDDH
USE YOMMDDH, ONLY  : TMDDH
!
!
REAL,                     INTENT(IN)    :: PTSTEP     ! Double Time step
                                                      ! (single if cold start)
INTEGER,                  INTENT(IN)    :: KMI        ! Model index 
CHARACTER(LEN=*),         INTENT(IN)    :: HFMFILE    ! Name of the output FM-file
CHARACTER(LEN=*),         INTENT(IN)    :: HLUOUT     ! Output-listing name for
                                                      ! model n
LOGICAL,                  INTENT(IN)    :: OCLOSE_OUT ! Conditional closure of 
                                                      ! the tput FM fileoutp
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF   ! Reference density
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: ZWLBDC3    ! Lambda(cloud) **3
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: ZWLBDC     ! Lambda(cloud)
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: ZWLBDR3    ! Lambda(rain) **3
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: ZWLBDR     ! Lambda(rain)
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRCT       ! Cloud water m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRRT       ! Rain water m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCCT       ! Cloud water C. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCRT       ! Rain water C. at t
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRCS       ! Cloud water m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRRS       ! Rain water m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCCS       ! Cloud water C. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCRS       ! Rain water C. source
!
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PRHODJ
!
TYPE(TYP_DDH), INTENT(INOUT) :: YDDDH
TYPE(TLDDH), INTENT(IN) :: YDLDDH
TYPE(TMDDH), INTENT(IN) :: YDMDDH
!
      END SUBROUTINE LIMA_WARM_COAL
END INTERFACE
END MODULE MODI_LIMA_WARM_COAL
!     #############################################################################
      SUBROUTINE LIMA_WARM_COAL (PTSTEP, KMI, HFMFILE, HLUOUT, OCLOSE_OUT,   &
                                 PRHODREF, ZWLBDC3, ZWLBDC, ZWLBDR3, ZWLBDR, &
                                 PRCT, PRRT, PCCT, PCRT,                     &
                                 PRCS, PRRS, PCCS, PCRS,                     &
                                 PRHODJ, &
                            YDDDH, YDLDDH, YDMDDH                              )
!     #############################################################################
!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the microphysical sources:
!!    nucleation, sedimentation, autoconversion, accretion, self-collection 
!!    and vaporisation which are parameterized according to Cohard and Pinty 
!!    QJRMS, 2000
!!
!!
!!**  METHOD
!!    ------
!!      Assuming a generalized gamma distribution law for the cloud droplets 
!!    and the raindrops, the zeroth and third order moments tendencies 
!!    are evaluated for all the coalescence terms by integrating the
!!    Stochastic Collection Equation. As autoconversion is a process that 
!!    cannot be resolved analytically, the Berry-Reinhardt parameterisation 
!!    is employed with modifications to initiate the raindrop spectrum mode.
!!     
!!    Computation steps :
!!      1- Check where computations are necessary, pack variables
!!      2- Self collection of cloud droplets
!!      3- Autoconversion of cloud droplets (Berry-Reinhardt parameterization)
!!      4- Accretion sources
!!      5- Self collection - Coalescence/Break-up
!!      6- Unpack variables, clean
!!
!!
!!    REFERENCE
!!    ---------
!!
!!      Cohard, J.-M. and J.-P. Pinty, 2000: A comprehensive two-moment warm 
!!      microphysical bulk scheme. 
!!        Part I: Description and tests
!!        Part II: 2D experiments with a non-hydrostatic model
!!      Accepted for publication in Quart. J. Roy. Meteor. Soc. 
!!
!!    AUTHOR
!!    ------
!!      J.-M. Cohard     * Laboratoire d'Aerologie*
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!      S.    Berthet    * Laboratoire d'Aerologie*
!!      B.    Vié        * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original             ??/??/13 
!!      C. Barthe  * LACy *  jan. 2014   add budgets
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAMETERS,     ONLY : JPHEXT, JPVEXT
USE MODD_PARAM_LIMA
USE MODD_PARAM_LIMA_WARM
!
USE MODD_NSV, ONLY : NSV_LIMA_NC, NSV_LIMA_NR
USE MODD_BUDGET
USE MODE_BUDGET, ONLY: BUDGET_DDH
!
USE MODI_LIMA_FUNCTIONS, ONLY : COUNTJV
!
USE DDH_MIX, ONLY  : TYP_DDH
USE YOMLDDH, ONLY  : TLDDH
USE YOMMDDH, ONLY  : TMDDH
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
REAL,                     INTENT(IN)    :: PTSTEP     ! Double Time step
                                                      ! (single if cold start)
INTEGER,                  INTENT(IN)    :: KMI        ! Model index 
CHARACTER(LEN=*),         INTENT(IN)    :: HFMFILE    ! Name of the output FM-file
CHARACTER(LEN=*),         INTENT(IN)    :: HLUOUT     ! Output-listing name for
                                                      ! model n
LOGICAL,                  INTENT(IN)    :: OCLOSE_OUT ! Conditional closure of 
                                                      ! the tput FM fileoutp
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF   ! Reference density
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: ZWLBDC3    ! Lambda(cloud) **3
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: ZWLBDC     ! Lambda(cloud)
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: ZWLBDR3    ! Lambda(rain) **3
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: ZWLBDR     ! Lambda(rain)
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRCT       ! Cloud water m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRRT       ! Rain water m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCCT       ! Cloud water C. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCRT       ! Rain water C. at t
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRCS       ! Cloud water m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRRS       ! Rain water m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCCS       ! Cloud water C. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCRS       ! Rain water C. source
!
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PRHODJ
!
TYPE(TYP_DDH), INTENT(INOUT) :: YDDDH
TYPE(TLDDH), INTENT(IN) :: YDLDDH
TYPE(TMDDH), INTENT(IN) :: YDMDDH
!
!*       0.1   Declarations of local variables :
!
! Packing variables
LOGICAL, DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3)) :: GMICRO 
INTEGER :: IMICRO
INTEGER , DIMENSION(SIZE(GMICRO))   :: I1,I2,I3 ! Used to replace the COUNT
INTEGER                             :: JL       ! and PACK intrinsics 
!
! Packed micophysical variables
REAL, DIMENSION(:)  , ALLOCATABLE :: ZRCT    ! Cloud water m.r. at t 
REAL, DIMENSION(:)  , ALLOCATABLE :: ZRRT    ! Rain water m.r. at t 
REAL, DIMENSION(:)  , ALLOCATABLE :: ZCCT    ! cloud conc. at t
REAL, DIMENSION(:)  , ALLOCATABLE :: ZCRT    ! rain conc. at t
!
REAL, DIMENSION(:)  , ALLOCATABLE :: ZRCS    ! Cloud water m.r. source
REAL, DIMENSION(:)  , ALLOCATABLE :: ZRRS    ! Rain water m.r. source
REAL, DIMENSION(:)  , ALLOCATABLE :: ZCCS    ! cloud conc. source
REAL, DIMENSION(:)  , ALLOCATABLE :: ZCRS    ! rain conc. source
!
! Other packed variables
REAL, DIMENSION(:)  , ALLOCATABLE  :: ZRHODREF ! RHO Dry REFerence
REAL, DIMENSION(:)  , ALLOCATABLE  :: ZLBDC3
REAL, DIMENSION(:)  , ALLOCATABLE  :: ZLBDC
REAL, DIMENSION(:)  , ALLOCATABLE  :: ZLBDR3
REAL, DIMENSION(:)  , ALLOCATABLE  :: ZLBDR
!
! Work arrays
REAL,    DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3)) :: ZW
!
REAL, DIMENSION(:), ALLOCATABLE    :: ZZW1, ZZW2, ZZW3, ZZW4, ZSCBU
LOGICAL, DIMENSION(:), ALLOCATABLE :: GSELF,               &
                                      GACCR,               &
                                      GSCBU,               &
                                      GENABLE_ACCR_SCBU
! 
!
INTEGER :: ISELF, IACCR, ISCBU
INTEGER :: IIB, IIE, IJB, IJE, IKB, IKE        ! Physical domain
!
!-------------------------------------------------------------------------------
!
!
!*       1.     PREPARE COMPUTATIONS - PACK
!   	        ---------------------------
!
!
IIB=1+JPHEXT
IIE=SIZE(PRHODREF,1) - JPHEXT
IJB=1+JPHEXT
IJE=SIZE(PRHODREF,2) - JPHEXT
IKB=1+JPVEXT
IKE=SIZE(PRHODREF,3) - JPVEXT
!
GMICRO(:,:,:) = .FALSE.
GMICRO(IIB:IIE,IJB:IJE,IKB:IKE) =                  &
     PRCT(IIB:IIE,IJB:IJE,IKB:IKE)>XRTMIN(2) .OR.  &
     PRRT(IIB:IIE,IJB:IJE,IKB:IKE)>XRTMIN(3)
!
IMICRO = COUNTJV( GMICRO(:,:,:),I1(:),I2(:),I3(:))
!
IF( IMICRO >= 1 ) THEN
   ALLOCATE(ZRCT(IMICRO))
   ALLOCATE(ZRRT(IMICRO))
   ALLOCATE(ZCCT(IMICRO))
   ALLOCATE(ZCRT(IMICRO))
!
   ALLOCATE(ZRCS(IMICRO))
   ALLOCATE(ZRRS(IMICRO))
   ALLOCATE(ZCCS(IMICRO))
   ALLOCATE(ZCRS(IMICRO))
!
   ALLOCATE(ZLBDC(IMICRO)) 
   ALLOCATE(ZLBDC3(IMICRO))
   ALLOCATE(ZLBDR(IMICRO)) 
   ALLOCATE(ZLBDR3(IMICRO))
! 
   ALLOCATE(ZRHODREF(IMICRO))
   DO JL=1,IMICRO
      ZCCT(JL) = PCCT(I1(JL),I2(JL),I3(JL))
      ZRCT(JL) = PRCT(I1(JL),I2(JL),I3(JL))
      ZRRT(JL) = PRRT(I1(JL),I2(JL),I3(JL))
      ZCRT(JL) = PCRT(I1(JL),I2(JL),I3(JL))
      ZCCS(JL) = PCCS(I1(JL),I2(JL),I3(JL))
      ZRCS(JL) = PRCS(I1(JL),I2(JL),I3(JL))
      ZRRS(JL) = PRRS(I1(JL),I2(JL),I3(JL))
      ZCRS(JL) = PCRS(I1(JL),I2(JL),I3(JL))
      ZLBDR(JL) = ZWLBDR(I1(JL),I2(JL),I3(JL))
      ZLBDR3(JL) = ZWLBDR3(I1(JL),I2(JL),I3(JL))
      ZLBDC(JL) = ZWLBDC(I1(JL),I2(JL),I3(JL))
      ZLBDC3(JL) = ZWLBDC3(I1(JL),I2(JL),I3(JL))
      ZRHODREF(JL) = PRHODREF(I1(JL),I2(JL),I3(JL))
   END DO
! 
   ALLOCATE(GSELF(IMICRO))
   ALLOCATE(GACCR(IMICRO))
   ALLOCATE(GSCBU(IMICRO))
   ALLOCATE(ZZW1(IMICRO))
   ALLOCATE(ZZW2(IMICRO))
   ALLOCATE(ZZW3(IMICRO))
!
!
!-------------------------------------------------------------------------------
!
!
!*       2. Self-collection of cloud droplets    
!   	 ------------------------------------
!
!
   GSELF(:) = ZCCT(:)>XCTMIN(2)
   ISELF = COUNT(GSELF(:))
   IF( ISELF>0 ) THEN
      ZZW1(:) = XSELFC*(ZCCT(:)/ZLBDC3(:))**2 * ZRHODREF(:) ! analytical integration
      WHERE( GSELF(:) )
         ZCCS(:) = ZCCS(:) - MIN( ZCCS(:),ZZW1(:) )
      END WHERE
   END IF
!
!
  ZW(:,:,:) = PCCS(:,:,:)
  IF (LBUDGET_SV) CALL BUDGET_DDH (                                 &
                   UNPACK(ZCCS(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:))&
                   &*PRHODJ(:,:,:),12+NSV_LIMA_NC,'SELF_BU_RSV',YDDDH, YDLDDH, YDMDDH) 
!
!
!-------------------------------------------------------------------------------
!
!
!*       3. Autoconversion of cloud droplets (Berry-Reinhardt parameterization)
!   	 ----------------------------------------------------------------------
!
!
IF (LRAIN_LIMA) THEN
!
   ZZW2(:) = 0.0
   ZZW1(:) = 0.0
   WHERE( ZRCT(:)>XRTMIN(2) )
      ZZW2(:) = MAX( 0.0,XLAUTR*ZRHODREF(:)*ZRCT(:)*             &
                           (XAUTO1/ZLBDC(:)**4-XLAUTR_THRESHOLD) ) ! L 
!
      ZZW3(:) = MIN( ZRCS(:), MAX( 0.0,XITAUTR*ZZW2(:)*ZRCT(:)*  &
                           (XAUTO2/ZLBDC(:)-XITAUTR_THRESHOLD) ) ) ! L/tau
!
      ZRCS(:) = ZRCS(:) - ZZW3(:)
      ZRRS(:) = ZRRS(:) + ZZW3(:)
!
      ZZW1(:) = MIN( MIN( 1.2E4,(XACCR4/ZLBDC(:)-XACCR5)/XACCR3),   &
                           ZLBDR(:)/XACCR1 ) ! D**-1 threshold diameter for 
                                             ! switching the autoconversion regimes
                                             ! min (80 microns, D_h, D_r)
      ZZW3(:) = ZZW3(:) * MAX( 0.0,ZZW1(:) )**3 / XAC 
      ZCRS(:) = ZCRS(:) + ZZW3(:)
   END WHERE
!
!
   ZW(:,:,:) = PRCS(:,:,:)
   IF (LBUDGET_RC) CALL BUDGET_DDH (                                  &
               UNPACK(ZRCS(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:)) &
                            *PRHODJ(:,:,:),7 ,'AUTO_BU_RRC',YDDDH, YDLDDH, YDMDDH)

   ZW(:,:,:) = PRRS(:,:,:)
   IF (LBUDGET_RR) CALL BUDGET_DDH (                                  &
               UNPACK(ZRRS(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:)) &
                            *PRHODJ(:,:,:),8 ,'AUTO_BU_RRR',YDDDH, YDLDDH, YDMDDH)
   ZW(:,:,:) = PCRS(:,:,:)
   IF (LBUDGET_SV) THEN
      ZW(:,:,:) = PCRS(:,:,:)
      CALL BUDGET_DDH (UNPACK(ZCRS(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:)) &
               *PRHODJ(:,:,:),12+NSV_LIMA_NR,'AUTO_BU_RSV',YDDDH, YDLDDH, YDMDDH)
      ZW(:,:,:) = PCCS(:,:,:)
      CALL BUDGET_DDH (UNPACK(ZCCS(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:)) &
               *PRHODJ(:,:,:),12+NSV_LIMA_NC,'AUTO_BU_RSV',YDDDH, YDLDDH, YDMDDH)
   END IF
!
!
!-------------------------------------------------------------------------------
!
!
!*       4. Accretion sources
!   	 --------------------
!
!
   GACCR(:) = ZRRT(:)>XRTMIN(3) .AND. ZCRT(:)>XCTMIN(3) 
   IACCR = COUNT(GACCR(:))
   IF( IACCR>0 ) THEN
      ALLOCATE(ZZW4(IMICRO)); ZZW4(:) = XACCR1/ZLBDR(:)
      ALLOCATE(GENABLE_ACCR_SCBU(IMICRO))
      GENABLE_ACCR_SCBU(:) = ZRRT(:)>1.2*ZZW2(:)/ZRHODREF(:) .OR.           &
                       ZZW4(:)>=MAX( XACCR2,XACCR3/(XACCR4/ZLBDC(:)-XACCR5) )
      GACCR(:) = GACCR(:) .AND. ZRCT(:)>XRTMIN(2) .AND. GENABLE_ACCR_SCBU(:)
   END IF
!
   IACCR = COUNT(GACCR(:))
   IF( IACCR>0 ) THEN
      WHERE( GACCR(:).AND.(ZZW4(:)>1.E-4) ) ! Accretion for D>100 10-6 m
         ZZW3(:) = ZLBDC3(:) / ZLBDR3(:)
         ZZW1(:) = ( ZCCT(:)*ZCRT(:) / ZLBDC3(:) )*ZRHODREF(:)
         ZZW2(:) = MIN( ZZW1(:)*(XACCR_CLARGE1+XACCR_CLARGE2*ZZW3(:)),ZCCS(:) )
         ZCCS(:) = ZCCS(:) - ZZW2(:)
!
         ZZW1(:) = ( ZZW1(:) / ZLBDC3(:) )*ZRHODREF(:)
         ZZW2(:) = MIN( ZZW1(:)*(XACCR_RLARGE1+XACCR_RLARGE2*ZZW3(:)),ZRCS(:) )
         ZRCS(:) = ZRCS(:) - ZZW2(:)
         ZRRS(:) = ZRRS(:) + ZZW2(:)
      END WHERE
      WHERE( GACCR(:).AND.(ZZW4(:)<=1.E-4) ) ! Accretion for D<100 10-6 m
         ZZW3(:) = ZLBDC3(:) / ZLBDR3(:)
         ZZW1(:) = ( ZCCT(:)*ZCRT(:) / ZLBDC3(:)**2 )*ZRHODREF(:)
         ZZW3(:) = ZZW3(:)**2
         ZZW2(:) = MIN( ZZW1(:)*(XACCR_CSMALL1+XACCR_CSMALL2*ZZW3(:)),ZCCS(:) )
         ZCCS(:) = ZCCS(:) - ZZW2(:)
!
         ZZW1(:) = ZZW1(:) / ZLBDC3(:)
         ZZW2(:) = MIN( ZZW1(:)*(XACCR_RSMALL1+XACCR_RSMALL2*ZZW3(:))             &
                                                          *ZRHODREF(:),ZRCS(:) )
         ZRCS(:) = ZRCS(:) - ZZW2(:)
         ZRRS(:) = ZRRS(:) + ZZW2(:)
      END WHERE
   END IF
!
!
   ZW(:,:,:) = PRCS(:,:,:)
   IF (LBUDGET_RC) CALL BUDGET_DDH (                                  &
               UNPACK(ZRCS(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:)) &
                              *PRHODJ(:,:,:),7 ,'ACCR_BU_RRC',YDDDH, YDLDDH, YDMDDH)
   ZW(:,:,:) = PRRS(:,:,:)
   IF (LBUDGET_RR) CALL BUDGET_DDH (                                  &
               UNPACK(ZRRS(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:)) &
                              *PRHODJ(:,:,:),8 ,'ACCR_BU_RRR',YDDDH, YDLDDH, YDMDDH)
   ZW(:,:,:) = PCCS(:,:,:)
   IF (LBUDGET_SV) CALL BUDGET_DDH (                                  &
               UNPACK(ZCCS(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:)) &
                  *PRHODJ(:,:,:),12+NSV_LIMA_NC,'ACCR_BU_RSV',YDDDH, YDLDDH, YDMDDH)
!
!
!-------------------------------------------------------------------------------
!
!
!*       5. Self collection - Coalescence/Break-up
!   	 -----------------------------------------
!
!
   IF( IACCR>0 ) THEN
      GSCBU(:) = ZCRT(:)>XCTMIN(3) .AND. GENABLE_ACCR_SCBU(:)
      ISCBU = COUNT(GSCBU(:))
   ELSE
      ISCBU = 0.0
   END IF
   IF( ISCBU>0 ) THEN
!
!*       5.1  efficiencies
!
      IF (.NOT.ALLOCATED(ZZW4)) ALLOCATE(ZZW4(IMICRO))
      ZZW4(:)  = XACCR1 / ZLBDR(:)                ! Mean diameter
      ALLOCATE(ZSCBU(IMICRO))
      ZSCBU(:) = 1.0
      WHERE (ZZW4(:)>=XSCBU_EFF1 .AND. GSCBU(:))   ZSCBU(:) = &  ! Coalescence
                            EXP(XSCBUEXP1*(ZZW4(:)-XSCBU_EFF1))  ! efficiency
      WHERE (ZZW4(:)>=XSCBU_EFF2) ZSCBU(:) = 0.0  ! Break-up
!
!*       5.2  integration
!
      ZZW1(:) = 0.0
      ZZW2(:) = 0.0
      ZZW3(:) = 0.0
      ZZW4(:) = XACCR1 / ZLBDR(:)                 ! Mean volume drop diameter
      WHERE (GSCBU(:).AND.(ZZW4(:)>1.E-4))              ! analytical integration
         ZZW1(:) = XSCBU2 * ZCRT(:)**2 / ZLBDR3(:)   ! D>100 10-6 m
         ZZW3(:) = ZZW1(:)*ZSCBU(:)
      END WHERE
      WHERE (GSCBU(:).AND.(ZZW4(:)<=1.E-4))
         ZZW2(:) = XSCBU3 *(ZCRT(:) / ZLBDR3(:))**2  ! D<100 10-6 m
         ZZW3(:) = ZZW2(:)
      END WHERE
      ZCRS(:) = ZCRS(:) - MIN( ZCRS(:),ZZW3(:) * ZRHODREF(:) )
      DEALLOCATE(ZSCBU)
   END IF
!
!
   ZW(:,:,:) = PCRS(:,:,:)
   IF (LBUDGET_SV) CALL BUDGET_DDH (                                  &
               UNPACK(ZCRS(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:)) &
                  *PRHODJ(:,:,:),12+NSV_LIMA_NR,'SCBU_BU_RSV',YDDDH, YDLDDH, YDMDDH)
!
END IF ! LRAIN_LIMA
!
!
!-------------------------------------------------------------------------------
!
!
!*       6. Unpack and clean
!   	 -------------------
!
!
   ZW(:,:,:) = PRCS(:,:,:)
   PRCS(:,:,:) = UNPACK( ZRCS(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
   ZW(:,:,:) = PRRS(:,:,:)
   PRRS(:,:,:) = UNPACK( ZRRS(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
   ZW(:,:,:) = PCCS(:,:,:)
   PCCS(:,:,:) = UNPACK( ZCCS(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
   ZW(:,:,:) = PCRS(:,:,:)
   PCRS(:,:,:) = UNPACK( ZCRS(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
!
   DEALLOCATE(ZRCT)
   DEALLOCATE(ZRRT)
   DEALLOCATE(ZCCT)
   DEALLOCATE(ZCRT)
   DEALLOCATE(ZRCS)
   DEALLOCATE(ZRRS)
   DEALLOCATE(ZCRS)
   DEALLOCATE(ZCCS)
   DEALLOCATE(ZRHODREF) 
   DEALLOCATE(GSELF)
   DEALLOCATE(GACCR)
   DEALLOCATE(GSCBU)
   IF( ALLOCATED(GENABLE_ACCR_SCBU) ) DEALLOCATE(GENABLE_ACCR_SCBU)
   DEALLOCATE(ZZW1)
   DEALLOCATE(ZZW2)
   DEALLOCATE(ZZW3)
   IF( ALLOCATED(ZZW4) ) DEALLOCATE(ZZW4)
   DEALLOCATE(ZLBDR3)
   DEALLOCATE(ZLBDC3)
   DEALLOCATE(ZLBDR)
   DEALLOCATE(ZLBDC)
!
!
!-------------------------------------------------------------------------------
!
ELSE
!*       7. Budgets are forwarded
!        ------------------------
!
!
   IF (LBUDGET_SV) CALL BUDGET_DDH (PCCS(:,:,:)*PRHODJ(:,:,:),12+NSV_LIMA_NC,'SELF_BU_RSV',YDDDH, YDLDDH, YDMDDH)
!
   IF (LBUDGET_RC) CALL BUDGET_DDH (PRCS(:,:,:)*PRHODJ(:,:,:),7 ,'AUTO_BU_RRC',YDDDH, YDLDDH, YDMDDH)
   IF (LBUDGET_RR) CALL BUDGET_DDH (PRRS(:,:,:)*PRHODJ(:,:,:),8 ,'AUTO_BU_RRR',YDDDH, YDLDDH, YDMDDH)
   IF (LBUDGET_SV) CALL BUDGET_DDH (PCRS(:,:,:)*PRHODJ(:,:,:),12+NSV_LIMA_NR,'AUTO_BU_RSV',YDDDH, YDLDDH, YDMDDH)
   IF (LBUDGET_SV) CALL BUDGET_DDH (PCCS(:,:,:)*PRHODJ(:,:,:),12+NSV_LIMA_NC,'AUTO_BU_RSV',YDDDH, YDLDDH, YDMDDH)
!
   IF (LBUDGET_RC) CALL BUDGET_DDH (PRCS(:,:,:)*PRHODJ(:,:,:),7 ,'ACCR_BU_RRC',YDDDH, YDLDDH, YDMDDH)
   IF (LBUDGET_RR) CALL BUDGET_DDH (PRRS(:,:,:)*PRHODJ(:,:,:),8 ,'ACCR_BU_RRR',YDDDH, YDLDDH, YDMDDH)
   IF (LBUDGET_SV) CALL BUDGET_DDH (PCCS(:,:,:)*PRHODJ(:,:,:),12+NSV_LIMA_NC,'ACCR_BU_RSV',YDDDH, YDLDDH, YDMDDH)
!
   IF (LBUDGET_SV) CALL BUDGET_DDH (PCRS(:,:,:)*PRHODJ(:,:,:),12+NSV_LIMA_NR,'SCBU_BU_RSV',YDDDH, YDLDDH, YDMDDH)

END IF ! IMICRO
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_WARM_COAL
