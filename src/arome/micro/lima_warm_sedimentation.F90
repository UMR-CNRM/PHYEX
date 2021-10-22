!      ###################################
       MODULE MODI_LIMA_WARM_SEDIMENTATION
!      ###################################
!
INTERFACE
      SUBROUTINE LIMA_WARM_SEDIMENTATION (OSEDC, KSPLITR, PTSTEP, KMI,  &
                                          HFMFILE, HLUOUT, OCLOSE_OUT,  &
                                          PZZ, PRHODREF, PPABST, ZT,    &
                                          ZWLBDC,                       &
                                          PRCT, PRRT, PCCT, PCRT,       &
                                          PRCS, PRRS, PCCS, PCRS,       &
                                          PINPRC, PINPRR, PINPRR3D )
!
LOGICAL,                  INTENT(IN)    :: OSEDC      ! switch to activate the 
                                                      ! cloud droplet sedimentation
INTEGER,                  INTENT(IN)    :: KSPLITR    ! Number of small time step 
                                                      ! for sedimendation
REAL,                     INTENT(IN)    :: PTSTEP     ! Double Time step
                                                      ! (single if cold start)
INTEGER,                  INTENT(IN)    :: KMI        ! Model index 
CHARACTER(LEN=*),         INTENT(IN)    :: HFMFILE    ! Name of the output FM-file
CHARACTER(LEN=*),         INTENT(IN)    :: HLUOUT     ! Output-listing name for
                                                      ! model n
LOGICAL,                  INTENT(IN)    :: OCLOSE_OUT ! Conditional closure of 
                                                      ! the tput FM fileoutp
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PZZ        ! Height (z)
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF   ! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABST     ! abs. pressure at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: ZT         ! Temperature
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: ZWLBDC     ! libre parcours moyen
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
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PINPRC     ! Cloud instant precip
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PINPRR     ! Rain instant precip
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PINPRR3D   ! Rain inst precip 3D
!
END SUBROUTINE LIMA_WARM_SEDIMENTATION
END INTERFACE
END MODULE MODI_LIMA_WARM_SEDIMENTATION
!     #####################################################################
      SUBROUTINE LIMA_WARM_SEDIMENTATION (OSEDC, KSPLITR, PTSTEP, KMI,  &
                                          HFMFILE, HLUOUT, OCLOSE_OUT,  &
                                          PZZ, PRHODREF, PPABST, ZT,    &
                                          ZWLBDC,                       &
                                          PRCT, PRRT, PCCT, PCRT,       &
                                          PRCS, PRRS, PCCS, PCRS,       &
                                          PINPRC, PINPRR, PINPRR3D )
!     #####################################################################
!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the sedimentation
!!    of cloud droplets and rain drops
!!
!!
!!**  METHOD
!!    ------
!!      The sedimentation rates are computed with a time spliting technique: 
!!    an upstream scheme, written as a difference of non-advective fluxes. 
!!    This source term is added to the next coming time step (split-implicit 
!!    process).
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
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAMETERS,      ONLY : JPHEXT, JPVEXT
USE MODD_CST,             ONLY : XRHOLW
USE MODD_PARAM_LIMA,      ONLY : XRTMIN, XCTMIN, XALPHAC, XNUC, XCEXVT
USE MODD_PARAM_LIMA_WARM, ONLY : XLBC, XLBEXC, XLBR, XLBEXR,        &
                                 XFSEDRC, XFSEDCC, XFSEDRR, XFSEDCR,&
                                 XDC, XDR
USE MODI_LIMA_FUNCTIONS,  ONLY : COUNTJV
USE MODI_GAMMA,           ONLY : GAMMA_X0D
!
USE YOMLUN   , ONLY : NULOUT
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
LOGICAL,                  INTENT(IN)    :: OSEDC      ! switch to activate the 
                                                      ! cloud droplet sedimentation
INTEGER,                  INTENT(IN)    :: KSPLITR    ! Number of small time step 
                                                      ! for sedimendation
REAL,                     INTENT(IN)    :: PTSTEP     ! Double Time step
                                                      ! (single if cold start)
INTEGER,                  INTENT(IN)    :: KMI        ! Model index 
CHARACTER(LEN=*),         INTENT(IN)    :: HFMFILE    ! Name of the output FM-file
CHARACTER(LEN=*),         INTENT(IN)    :: HLUOUT     ! Output-listing name for
                                                      ! model n
LOGICAL,                  INTENT(IN)    :: OCLOSE_OUT ! Conditional closure of 
                                                      ! the tput FM fileoutp
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PZZ        ! Height (z)
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF   ! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABST     ! abs. pressure at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: ZT         ! Temperature
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: ZWLBDC     ! libre parcours moyen
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
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PINPRC     ! Cloud instant precip
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PINPRR     ! Rain instant precip
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PINPRR3D   ! Rain inst precip 3D
!
!
!*       0.2   Declarations of local variables :
!
! Packing variables
LOGICAL, DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3)) :: GSEDIM 
INTEGER :: ISEDIM
INTEGER , DIMENSION(SIZE(GSEDIM)) :: I1,I2,I3 ! Used to replace the COUNT
INTEGER                           :: JL       ! and PACK intrinsics 
!
! Packed micophysical variables
REAL, DIMENSION(:)  , ALLOCATABLE :: ZRCT     ! Cloud water m.r. at t 
REAL, DIMENSION(:)  , ALLOCATABLE :: ZRRT     ! Rain water m.r. at t 
REAL, DIMENSION(:)  , ALLOCATABLE :: ZCCT     ! cloud conc. at t
REAL, DIMENSION(:)  , ALLOCATABLE :: ZCRT     ! rain conc. at t
!
REAL, DIMENSION(:)  , ALLOCATABLE :: ZRCS     ! Cloud water m.r. source
REAL, DIMENSION(:)  , ALLOCATABLE :: ZRRS     ! Rain water m.r. source
REAL, DIMENSION(:)  , ALLOCATABLE :: ZCCS     ! cloud conc. source
REAL, DIMENSION(:)  , ALLOCATABLE :: ZCRS     ! rain conc. source
!
! Other packed variables
REAL, DIMENSION(:)  , ALLOCATABLE :: ZRHODREF ! RHO Dry REFerence
REAL, DIMENSION(:)  , ALLOCATABLE :: ZLBDC    
REAL, DIMENSION(:)  , ALLOCATABLE :: ZLBDR    
!
! Work arrays
REAL,    DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3))   &
                                  :: ZW,             &
                                     ZWLBDA,         &  ! Mean free path
                                     ZRAY,           &  ! Mean volumic radius
                                     ZCC                ! Terminal vertical velocity
REAL,    DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3)+1)   &
                                  :: ZWSEDR, ZWSEDC     ! Sedim. fluxes
!
REAL, DIMENSION(:), ALLOCATABLE   :: ZZW1, ZZW2, ZZW3, &
                                     ZTCC,             &
                                     ZRTMIN, ZCTMIN
!
!
INTEGER :: JK            ! Vertical loop index for the rain sedimentation 
INTEGER :: JN            ! Temporal loop index for the rain sedimentation
INTEGER :: IIB, IIE, IJB, IJE, IKB, IKE           ! Physical domain
REAL    :: ZTSPLITR      ! Small time step for rain sedimentation
!
INTEGER :: IKMAX 
!
!!!!!!!!!!! Entiers pour niveaux inversés dans AROME !!!!!!!!!!!!!!!!!!!!!!!!!!!
INTEGER :: IBOTTOM, INVLVL
!
!-------------------------------------------------------------------------------
!
!        0. Prepare computations
!        -----------------------
!
!
ALLOCATE(ZRTMIN(SIZE(XCTMIN)))
ALLOCATE(ZCTMIN(SIZE(XCTMIN)))
ZRTMIN(:) = XRTMIN(:) / PTSTEP
ZCTMIN(:) = XCTMIN(:) / PTSTEP
!
IIB=1+JPHEXT
IIE=SIZE(PZZ,1) - JPHEXT
IJB=1+JPHEXT
IJE=SIZE(PZZ,2) - JPHEXT
IKB=1+JPVEXT
IKE=SIZE(PZZ,3) - JPVEXT
!
!!!!!!!!!!! Entiers pour niveaux inversés dans AROME !!!!!!!!!!!!!!!!!!!!!!!!!!!
IBOTTOM=IKE
INVLVL=-1
!
ZWSEDR(:,:,:)=0.
ZWSEDC(:,:,:)=0.
IKMAX=SIZE(PRHODREF,3)
!
ZTSPLITR= PTSTEP / FLOAT(KSPLITR)
!
PINPRC(:,:) = 0.
PINPRR(:,:) = 0.
PINPRR3D(:,:,:) = 0.
!
IF (OSEDC) THEN
   ZWLBDA(:,:,:) = 0.
   ZRAY(:,:,:)   = 0.
   ZCC(:,:,:) = 1. 
   DO JK=IKB,IKE
      ZWLBDA(:,:,JK) = 6.6E-8*(101325./PPABST(:,:,JK))*(ZT(:,:,JK)/293.15)
   END DO
   WHERE (PRCT(:,:,:)>XRTMIN(2) .AND. PCCT(:,:,:)>XCTMIN(2))
      ZRAY(:,:,:) = 0.5*GAMMA_X0D(XNUC+1./XALPHAC)/(GAMMA_X0D(XNUC)*ZWLBDC(:,:,:))
      ! ZCC : Corrective Cunningham term for the terminal velocity
      ZCC(:,:,:)=1.+1.26*ZWLBDA(:,:,:)/ZRAY(:,:,:)
   END WHERE
END IF
!
!-------------------------------------------------------------------------------
!
!
!        1. Computations only where necessary
!        ------------------------------------
!
!
DO JN = 1 , KSPLITR
   GSEDIM(:,:,:) = .FALSE.
   GSEDIM(IIB:IIE,IJB:IJE,IKB:IKE) = PRRS(IIB:IIE,IJB:IJE,IKB:IKE)>ZRTMIN(3) &
                               .AND. PCRS(IIB:IIE,IJB:IJE,IKB:IKE)>ZCTMIN(3)
   IF( OSEDC ) THEN
      GSEDIM(IIB:IIE,IJB:IJE,IKB:IKE) = GSEDIM(IIB:IIE,IJB:IJE,IKB:IKE) .OR.    &
                                       (PRCS(IIB:IIE,IJB:IJE,IKB:IKE)>ZRTMIN(2) &
                                  .AND. PCCS(IIB:IIE,IJB:IJE,IKB:IKE)>ZCTMIN(2) )
   END IF
!
   ISEDIM = COUNTJV( GSEDIM(:,:,:),I1(:),I2(:),I3(:))
   IF( ISEDIM >= 1 ) THEN
!
      IF( JN==1 ) THEN
         IF( OSEDC ) THEN
            PRCS(:,:,:) = PRCS(:,:,:) * PTSTEP
            PCCS(:,:,:) = PCCS(:,:,:) * PTSTEP
         END IF
         PRRS(:,:,:) = PRRS(:,:,:) * PTSTEP
         PCRS(:,:,:) = PCRS(:,:,:) * PTSTEP
         DO JK = IKB , IKE
!Dans AROME, PZZ = épaisseur de la couche
!            ZW(:,:,JK)=ZTSPLITR/(PZZ(:,:,JK+1)-PZZ(:,:,JK))
            ZW(:,:,JK)=ZTSPLITR/(PZZ(:,:,JK))
         END DO
      END IF
!
      ALLOCATE(ZRHODREF(ISEDIM))
      DO JL = 1,ISEDIM
         ZRHODREF(JL) = PRHODREF(I1(JL),I2(JL),I3(JL))
      END DO
!
      ALLOCATE(ZZW1(ISEDIM)) 
      ALLOCATE(ZZW2(ISEDIM)) 
      ALLOCATE(ZZW3(ISEDIM)) 
!
!
!-------------------------------------------------------------------------------
!
!
!        2. Cloud droplets sedimentation
!        -------------------------------
!
!
      IF( OSEDC .AND. MAXVAL(PRCS(:,:,:))>ZRTMIN(2) ) THEN
         ZZW1(:) = 0.0
         ZZW2(:) = 0.0
         ZZW3(:) = 0.0
         ALLOCATE(ZRCS(ISEDIM))
         ALLOCATE(ZCCS(ISEDIM))
         ALLOCATE(ZRCT(ISEDIM))
         ALLOCATE(ZCCT(ISEDIM))
         ALLOCATE(ZTCC(ISEDIM))
         ALLOCATE(ZLBDC(ISEDIM))
         DO JL = 1,ISEDIM
            ZRCS(JL) = PRCS(I1(JL),I2(JL),I3(JL))
            ZCCS(JL) = PCCS(I1(JL),I2(JL),I3(JL))
            ZRCT(JL) = PRCT(I1(JL),I2(JL),I3(JL))
            ZCCT(JL) = PCCT(I1(JL),I2(JL),I3(JL))
            ZTCC(JL) = ZCC (I1(JL),I2(JL),I3(JL))
         END DO
         ZLBDC(:) = 1.E15
         WHERE (ZRCT(:)>XRTMIN(2) .AND. ZCCT(:)>XCTMIN(2))
            ZLBDC(:) = ( XLBC*ZCCT(:) / ZRCT(:) )**XLBEXC
         END WHERE
         WHERE( ZRCS(:)>ZRTMIN(2) )
            ZZW3(:) = ZRHODREF(:)**(-XCEXVT) * ZLBDC(:)**(-XDC)
            ZZW1(:) = ZTCC(:) * XFSEDRC * ZRCS(:) * ZZW3(:) * ZRHODREF(:)
            ZZW2(:) = ZTCC(:) * XFSEDCC * ZCCS(:) * ZZW3(:) * ZRHODREF(:) 
         END WHERE
         ZWSEDR(:,:,1:IKMAX) = UNPACK( ZZW1(:),MASK=GSEDIM(:,:,:),FIELD=0.0 )
         ZWSEDC(:,:,1:IKMAX) = UNPACK( ZZW2(:),MASK=GSEDIM(:,:,:),FIELD=0.0 )
         DO JK = IKB+1 , IKE
            PRCS(:,:,JK) = PRCS(:,:,JK) + ZW(:,:,JK)*                       &
                 (ZWSEDR(:,:,JK+1*INVLVL)-ZWSEDR(:,:,JK))/PRHODREF(:,:,JK)
            PCCS(:,:,JK) = PCCS(:,:,JK) + ZW(:,:,JK)*                       &
                 (ZWSEDC(:,:,JK+1*INVLVL)-ZWSEDC(:,:,JK))/PRHODREF(:,:,JK)
         END DO
            PRCS(:,:,1) = PRCS(:,:,1) + ZW(:,:,1)*                       &
                 (0.-ZWSEDR(:,:,1))/PRHODREF(:,:,1)
            PCCS(:,:,1) = PCCS(:,:,1) + ZW(:,:,1)*                       &
                 (0.-ZWSEDC(:,:,1))/PRHODREF(:,:,1)
         DEALLOCATE(ZRCS)
         DEALLOCATE(ZCCS)
         DEALLOCATE(ZRCT)
         DEALLOCATE(ZCCT)
         DEALLOCATE(ZTCC)
         DEALLOCATE(ZLBDC)
!
         PINPRC(:,:) = PINPRC(:,:) + ZWSEDR(:,:,IBOTTOM)/XRHOLW/KSPLITR                        ! in m/s
      ELSE
         ZWSEDR(:,:,IBOTTOM) = 0.0
      END IF ! OSEDC
!
!
!-------------------------------------------------------------------------------
!
!
!        2. Rain drops sedimentation
!        ---------------------------
!
!
      IF( MAXVAL(PRRS(:,:,:))>ZRTMIN(3) ) THEN
         ZZW1(:) = 0.0
         ZZW2(:) = 0.0
         ZZW3(:) = 0.0
         ALLOCATE(ZRRS(ISEDIM)) 
         ALLOCATE(ZCRS(ISEDIM))
         ALLOCATE(ZRRT(ISEDIM)) 
         ALLOCATE(ZCRT(ISEDIM))
         ALLOCATE(ZLBDR(ISEDIM))
         DO JL = 1,ISEDIM
            ZRRS(JL) = PRRS(I1(JL),I2(JL),I3(JL))
            ZCRS(JL) = PCRS(I1(JL),I2(JL),I3(JL))
            ZRRT(JL) = PRRT(I1(JL),I2(JL),I3(JL))
            ZCRT(JL) = PCRT(I1(JL),I2(JL),I3(JL))
         END DO
         ZLBDR(:) = 1.E10
         WHERE (ZRRT(:)>XRTMIN(3) .AND. ZCRT(:)>XCTMIN(3))
            ZLBDR(:) = ( XLBR*ZCRT(:) / ZRRT(:) )**XLBEXR
         END WHERE
         WHERE( ZRRS(:)>ZRTMIN(3) )
            ZZW3(:) = ZRHODREF(:)**(-XCEXVT) * (ZLBDR(:)**(-XDR))
            ZZW1(:) = XFSEDRR * ZRRS(:) * ZZW3(:) * ZRHODREF(:)
            ZZW2(:) = XFSEDCR * ZCRS(:) * ZZW3(:) * ZRHODREF(:) 
         END WHERE
         ZWSEDR(:,:,1:IKMAX) = UNPACK( ZZW1(:),MASK=GSEDIM(:,:,:),FIELD=0.0 )
         ZWSEDC(:,:,1:IKMAX) = UNPACK( ZZW2(:),MASK=GSEDIM(:,:,:),FIELD=0.0 )
         DO JK = IKB+1 , IKE
            PRRS(:,:,JK) = PRRS(:,:,JK) + ZW(:,:,JK)*                      &
                 (ZWSEDR(:,:,JK+1*INVLVL)-ZWSEDR(:,:,JK))/PRHODREF(:,:,JK)
            PCRS(:,:,JK) = PCRS(:,:,JK) + ZW(:,:,JK)*                      &
                 (ZWSEDC(:,:,JK+1*INVLVL)-ZWSEDC(:,:,JK))/PRHODREF(:,:,JK)
         END DO
            PRRS(:,:,1) = PRRS(:,:,1) + ZW(:,:,1)*                      &
                 (0.-ZWSEDR(:,:,1))/PRHODREF(:,:,1)
            PCRS(:,:,1) = PCRS(:,:,1) + ZW(:,:,1)*                      &
                 (0.-ZWSEDC(:,:,1))/PRHODREF(:,:,1)
         DEALLOCATE(ZRRS)
         DEALLOCATE(ZCRS)
         DEALLOCATE(ZRRT)
         DEALLOCATE(ZCRT)
         DEALLOCATE(ZLBDR)
      ELSE
         ZWSEDR(:,:,IBOTTOM) = 0.0
      END IF ! max PRRS > ZRTMIN(3)
!    
      PINPRR(:,:) = PINPRR(:,:) + ZWSEDR(:,:,IBOTTOM)/XRHOLW/KSPLITR              ! in m/s
      PINPRR3D(:,:,:) = PINPRR3D(:,:,:) + ZWSEDR(:,:,1:IKMAX)/XRHOLW/KSPLITR        ! in m/s
!
      DEALLOCATE(ZRHODREF)
      DEALLOCATE(ZZW1)
      DEALLOCATE(ZZW2)
      DEALLOCATE(ZZW3)
      IF( JN==KSPLITR ) THEN
         IF( OSEDC ) THEN
            PRCS(:,:,:) = PRCS(:,:,:) / PTSTEP
            PCCS(:,:,:) = PCCS(:,:,:) / PTSTEP
         END IF
         PRRS(:,:,:) = PRRS(:,:,:) / PTSTEP
         PCRS(:,:,:) = PCRS(:,:,:) / PTSTEP
      END IF
   END IF ! ISEDIM
END DO ! KSPLITR
!
!++cb++
DEALLOCATE(ZRTMIN)
DEALLOCATE(ZCTMIN)
!--cb--

!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_WARM_SEDIMENTATION
