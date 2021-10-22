
SUBROUTINE ARO_WETDEP(KLON,KLEV,NSV,KRR,PTSTEP   &
       ,PSVT        & !I [moments/molec_{air}] Transported moments of dust
       ,PZZ         & !I [m] height of layers
       ,PPABST      & ! Pressure
       ,PTHT        & ! Potential temperature
       ,PRHODREF    & !I [kg/m3] density of air
       ,KTCOUNT     & ! number of time step
       ,PRT         & ! moist field
       ,PEVAP       & ! evaporation profile
       ,KSPLITR     & ! rain sedimentation time splitting
       )
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!
!*** *ARO_WETDEP*
!
!    PURPOSE
!    -------

!     Interface routine for aerosol scavenging
  
!     AUTHOR.
!     -------
!      P. Tulet
!      10-03-08

!    MODIFICATIONS
!    -------------
!    01-02-2011 M. Mokhtari : Adaptation of ARO_RAINAERO under ARO_WETDEP for Aladin
!    M. Mokhtari & A. Ambar 09/2016: inversion levels, call inv_levels.F90
!
!    EXTERNAL
!     -------



USE MODI_AER_WET_DEP
USE MODD_DUST
USE MODD_CSTS_DUST
USE MODD_CST
USE MODD_NSV, ONLY : NSV_DSTBEG, NSV_DSTEND, &
                     NSV_DSTDEPBEG, NSV_DSTDEPEND
USE MODE_DUST_PSD
USE MODI_INV_LEVELS

 IMPLICIT NONE
 
!INPUT
INTEGER,   INTENT(IN)   :: KLON     ! NPROMA under CPG
INTEGER,   INTENT(IN)   :: KLEV     ! Number of vertical levels
INTEGER,   INTENT(IN)   :: NSV      ! Number of passive scalar
INTEGER,   INTENT(IN)   :: KRR      ! Number of moist variables
REAL,      INTENT(IN)   :: PTSTEP   ! Time step

REAL, DIMENSION(KLON,1,KLEV,NSV),INTENT(INOUT) :: PSVT 
REAL, DIMENSION(KLON,1,KLEV,KRR), INTENT(IN)   :: PRT   ! Moist variables at time t
             !I [moments/molec_{air}] transported moments of dust
REAL, DIMENSION(KLON,1,KLEV),INTENT(IN)        :: PZZ        !I [m] height of layers
REAL, DIMENSION(KLON,1,KLEV),INTENT(IN)        :: PRHODREF   !I [kg/m3] density of air
REAL, DIMENSION(KLON,1,KLEV),INTENT(IN)        :: PTHT       !I [K] potentiel temperature
REAL, DIMENSION(KLON,1,KLEV),INTENT(IN)        :: PPABST     !I [Pa] pressure
REAL, DIMENSION(KLON,1,KLEV),INTENT(IN)        :: PEVAP      !I Evaporation
INTEGER, INTENT(IN)        :: KTCOUNT    ! Number of time step
INTEGER, INTENT(IN)        :: KSPLITR    ! Rain sedimentation time spliting
   
INTEGER :: JLEV,JK, JL, JK1, IMODEIDX
INTEGER :: II, IKE, JSV, JRR
REAL    :: ZSIGMIN
!  local variable
REAL, DIMENSION(:),       ALLOCATABLE  :: ZMASSMIN, ZINIRADIUS
REAL,DIMENSION(KLON,1,KLEV+2)         :: ZZZ        ! Local value of PZZ
REAL,DIMENSION(KLON,1,KLEV+2)         :: ZPABST     ! Local value of ZPABST
REAL,DIMENSION(KLON,1,KLEV+2)         :: ZTHT       ! Local value of ZTHT
REAL,DIMENSION(KLON,1,KLEV+2)         :: ZRHODREF   ! Local value of ZRHODREF
REAL,DIMENSION(KLON,1,KLEV+2)         :: ZEVAP      ! Evaporation
REAL,DIMENSION(KLON,1,KLEV+2,KRR)     :: ZRT        ! Moist variables at time t
REAL,DIMENSION(KLON,1,KLEV+2,NSV)     :: ZSVT       ! Local value of ZSVT
!
REAL,DIMENSION(KLON,1,KLEV)         :: ZZPABST     ! Local value of ZPABST
REAL,DIMENSION(KLON,1,KLEV)         :: ZZTHT       ! Local value of ZTHT
REAL,DIMENSION(KLON,1,KLEV)         :: ZZRHODREF   ! Local value of ZRHODREF
REAL,DIMENSION(KLON,1,KLEV)         :: ZZEVAP      !I Evaporation
REAL,DIMENSION(KLON,1,KLEV,KRR)     :: ZZRT        ! Moist variables at time t
REAL,DIMENSION(KLON,1,KLEV,NSV)     :: ZZSVT       ! Local value of ZSVT
!
REAL, DIMENSION(KLON,1,KLEV+2,NMODE_DST)   :: ZSIGDST, ZRGDST, ZNDST,ZVMASSMIN
REAL, DIMENSION(KLON,1,KLEV+2,NMODE_DST*3) :: ZSVDST
!--------------------------------------------------------------------

! 1) Initialization of aerosols acqueous variables
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('ARO_WETDEP',0,ZHOOK_HANDLE)

ALLOCATE(ZMASSMIN(NMODE_DST))
ALLOCATE(ZINIRADIUS(NMODE_DST))

!
!*       1.     PRELIMINARY COMPUTATIONS
!*   INVERS LEVEL FOR AROME
  !initialisation de ZZZ
DO JL = 1,KLON
   DO JK = 2 , KLEV+1
      ZZZ(JL,1,JK)=PZZ(JL,1,KLEV+2-JK)
   ENDDO
   ZZZ(JL,1,1) = 2*ZZZ(JL,1,2)-ZZZ(JL,1,3)
   ZZZ(JL,1,KLEV+2) = 2*ZZZ(JL,1,KLEV+1)-ZZZ(JL,1,KLEV)
ENDDO
!
ZZRHODREF(:,:,:) = PRHODREF(:,:,:)
ZZPABST(:,:,:)   = PPABST(:,:,:)
ZZTHT(:,:,:)     = PTHT(:,:,:)
ZZEVAP(:,:,:)    = PEVAP(:,:,:)
ZZRT(:,:,:,:)    = PRT(:,:,:,:)
ZZSVT(:,:,:,:)   = PSVT(:,:,:,:)
ZSVDST(:,:,:,:) = 0.
!
CALL INV_LEVELS(1,KLON,KLEV,1,ZZRHODREF,ZRHODREF)
CALL INV_LEVELS(1,KLON,KLEV,1,ZZPABST,ZPABST)
CALL INV_LEVELS(1,KLON,KLEV,1,ZZTHT,ZTHT)
CALL INV_LEVELS(1,KLON,KLEV,1,ZZEVAP,ZEVAP)
!
DO JSV=1,NSV
  CALL INV_LEVELS(1,KLON,KLEV,1,ZZSVT(:,:,:,JSV),ZSVT(:,:,:,JSV))
ENDDO
!
DO JRR=1,KRR
  CALL INV_LEVELS(1,KLON,KLEV,1,ZZRT(:,:,:,JRR),ZRT(:,:,:,JRR))
ENDDO

IF (KTCOUNT == 1) THEN
  ZSVT(:,:,:,NSV_DSTDEPBEG:NSV_DSTDEPEND) =  0.
END IF

!     3.1 Minimum mass to transfer between dry mass or in-cloud droplets

DO JSV=1,NMODE_DST
  IMODEIDX = JPDUSTORDER(JSV)
   IF (CRGUNITD=="MASS") THEN
    ZINIRADIUS(JSV) = XINIRADIUS(IMODEIDX) * EXP(-3.*(LOG(XINISIG(IMODEIDX)))**2)
   ELSE
    ZINIRADIUS(JSV) = XINIRADIUS(IMODEIDX)
   END IF
   IF (LVARSIG) THEN
    ZSIGMIN = XSIGMIN
   ELSE
    ZSIGMIN = XINISIG(IMODEIDX)
   ENDIF
   ZMASSMIN(JSV) = XN0MIN(IMODEIDX) * (ZINIRADIUS(JSV)**3)*EXP(4.5 * LOG(ZSIGMIN)**2)
! volume/um3 =>  #/molec_{air}
   ZVMASSMIN(:,:,:,JSV)=  ZMASSMIN(JSV) * XMD * XPI * 4./3. * XDENSITY_DUST  / &
                 (XMOLARWEIGHT_DUST*XM3TOUM3*ZRHODREF(:,:,:))
ENDDO
  
!
!     3.3 Compute and store Standard deviation and mean radius 
!         from moments
  CALL PPP2DUST(ZSVT(:,:,:,NSV_DSTBEG:NSV_DSTEND), &
                ZRHODREF(:,:,:),                   &
                ZSIGDST(:,:,:,:),                  &
                ZRGDST(:,:,:,:),                   &
                ZNDST(:,:,:,:))

!     3.4 Compute acquous aerosol mass vector from moment scalar vector
!
  DO JSV= 1, NMODE_DST
   IF (LVARSIG) THEN
    ZSVDST(:,:,:,JSV) = ZSVT(:,:,:,NSV_DSTBEG+1+(JSV-1)*3)
   ELSE IF (LRGFIX_DST) THEN
    ZSVDST(:,:,:,JSV) = ZSVT(:,:,:,NSV_DSTBEG+JSV-1)
   ELSE
    ZSVDST(:,:,:,JSV) = ZSVT(:,:,:,NSV_DSTBEG+1+(JSV-1)*2)
   ENDIF
  ENDDO    
  DO JSV=1,2*NMODE_DST
    ZSVDST(:,:,:,NMODE_DST+JSV) = ZSVT(:,:,:,NSV_DSTDEPBEG-1+JSV)
  ENDDO

! One moment cloud scheme
  CALL AER_WET_DEP           (KSPLITR,          &
                              PTSTEP,           &
                              ZZZ(:,:,:),       &
                              ZRHODREF(:,:,:),  &
                              ZRT(:,:,:,2),     &
                              ZRT(:,:,:,3),     &
                              ZRT(:,:,:,2),     &
                              ZRT(:,:,:,3),     &
                              ZSVDST(:,:,:,:),  &
                              ZTHT(:,:,:),      &
                              ZPABST(:,:,:),    &
                              ZRGDST(:,:,:,:),  &
                              ZEVAP(:,:,:),     &
                              NMODE_DST,        &
                              XDENSITY_DUST,    &
                              ZVMASSMIN         )

!      Compute return to moment vector

  DO JSV=1,NMODE_DST
   IF (LVARSIG) THEN
    ZSVT(:,:,:,NSV_DSTBEG+1+(JSV-1)*3) = ZSVDST(:,:,:,JSV)
   ELSE IF (LRGFIX_DST) THEN
    ZSVT(:,:,:,NSV_DSTBEG+JSV-1) = ZSVDST(:,:,:,JSV)
   ELSE
    ZSVT(:,:,:,NSV_DSTBEG+1+(JSV-1)*2) = ZSVDST(:,:,:,JSV)
   ENDIF
  ENDDO    
                             
!      Return to lognormal distribution (compute M0 and M6 using RG, SIG and
!     new mass from M3)

  CALL DUST2PPP(ZSVT(:,:,:,NSV_DSTBEG:NSV_DSTEND), &
                ZRHODREF(:,:,:),                   &
                ZSIGDST(:,:,:,:),                  &
                ZRGDST(:,:,:,:))

  DO JSV=1,2*NMODE_DST
     ZSVT(:,:,:,NSV_DSTDEPBEG-1+JSV) =  ZSVDST(:,:,:,NMODE_DST+JSV) 
  ENDDO

! Return to aladin/arome
DO JSV=1,NSV
  CALL INV_LEVELS(1,KLON,KLEV,-1,PSVT(:,:,:,JSV),ZSVT(:,:,:,JSV))
ENDDO

  DEALLOCATE(ZINIRADIUS)
  DEALLOCATE(ZMASSMIN)

IF (LHOOK) CALL DR_HOOK('ARO_WETDEP',1,ZHOOK_HANDLE)
END SUBROUTINE ARO_WETDEP

