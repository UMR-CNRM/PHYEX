!     ######spl
      SUBROUTINE AER_WET_DEP          (KSPLITR, PTSTEP, PZZ,      &
                            PRHODREF, PRCT, PRRT,                 &
                            PRCS, PRRS,  PSVT, PTHT,              &
                            PPABST, PRGAER, PEVAP3D, KMODE,       &
                            PDENSITY_AER, PMASSMIN, PSEA, PTOWN,  &
                            PCCT, PCRT )
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!     ###############################################################
!
!!****  * -  compute the explicit microphysical processes involved in the
!!***   * -  wet deposition of aerosols species in mixed clouds
!!
!!    PURPOSE
!!    -------
!!
!!  The purpose of this subroutine is to calculate the mass transfer
!!  of aerosol species between cloud hydrometeors.
!! 
!!
!!
!!**  METHOD
!!    ------
!!      Aerosols mass are dissolved into the cloud water and rain
!!      drops, it is subject to transfer through the microphysical processes
!!      that affect the parent hydrometeor [Rutledge et al., 1986].
!!      Aerosol mass transfer has been computed using scavenging coefficient 
!!      and brownian nucleation scavenging coefficient (Seinfeld and Pandis,
!!      1998; Tost et al, 2006).
!!
!!      The sedimentation rate is computed with a time spliting technique and
!!      an upstream scheme, written as a difference of non-advective fluxes.
!!
!!     KMODE: Number of aerosol modes (lognormal, bin..)
!!     PSVT : 1 => KMODE          : dry aerosol mass
!!     PSVT : KMODE+1 => 2*KMODE  : aerosol mass in cloud
!!     PSVT : 2*KMODE+1 => 3*KMODE: aerosol mass in rain

!!
!!    EXTERNAL
!!    --------
!!      None
!!     
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST     
!!          XP00               ! Reference pressure
!!          XRD,XRV            ! Gaz  constant for dry air, vapor
!!          XMD,XMV            ! Molecular weight for dry air, vapor
!!          XCPD               ! Cpd (dry air)
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      P. Tulet & K. Crahan-Kaku      * CNRM * 
!!
!!      Based on rain_ice.f90 and ch_wet_dep_kmt_warm.f90
!!      from C. Mari & J.P. Pinty  * LA*
!!
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    09/05/07
!!      01-02-2011 M. Mokhtari: Adaptation of aer_wet_dep_kmt_warm.f90 under aer_wet_dep for Aladin  
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
USE MODD_WET_DEP_PARAM
USE MODD_WET_DEP_DESCR
USE MODI_AER_VELGRAV
USE MODI_AER_EFFIC_DEP
USE MODI_GAMMA
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER,                  INTENT(IN)    :: KSPLITR ! Number of small time step 
                                      ! integration for  rain sedimendation
REAL,                     INTENT(IN)    :: PTSTEP  ! Time step          
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PZZ     ! Height (z)
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF! Reference [kg/m3] air density
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRCT    ! Cloud water m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRRT    ! Rain water m.r. at t 
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PSVT    ! Tracer m.r. at t
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRCS    ! Cloud water m.r. from source term
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRRS    ! Rain water m.r. from source term 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PEVAP3D ! Instantaneous 3D Rain Evaporation flux (KG/KG/S) 
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT     ! Potential temp
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABST   ! [Pa] pressure
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRGAER   ! Aerosols radius (um)
INTEGER,                  INTENT(IN)    :: KMODE    ! Nb aerosols mode
REAL,                     INTENT(IN)    :: PDENSITY_AER ! Begin Index for aerosol in cloud
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PMASSMIN     ! Aerosol mass minimum value 
REAL, DIMENSION(:,:),OPTIONAL, INTENT(IN)   :: PSEA  ! Sea mask
REAL, DIMENSION(:,:),OPTIONAL, INTENT(IN)   :: PTOWN ! Town mask
REAL, DIMENSION(:,:,:),OPTIONAL, INTENT(IN) :: PCCT   ! Cloud water concentration
REAL, DIMENSION(:,:,:),OPTIONAL, INTENT(IN) :: PCRT   ! Rain water concentration

!
!*       0.2   Declarations of local variables :
!
INTEGER :: JK            ! Vertical loop index for the rain sedimentation 
INTEGER :: JN            ! Temporal loop index for the rain sedimentation
INTEGER :: JJ, JI        ! Loop index for the interpolation
!
REAL    :: ZTSPLITR      ! Small time step for rain sedimentation
!
REAL,  DIMENSION(:,:), ALLOCATABLE :: ZEFC  !efficiency factor [unitless]
!
!Declaration of Dust Variables
!
INTEGER :: ICLOUD, IRAIN
! Case number of sedimentation, T>0 (for HEN)
                         ! and r_x>0 locations
LOGICAL, DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3)) &
                    :: GRAIN, GCLOUD ! Test where to compute all processes
 ! Test where to compute the SED/EVAP processes
REAL,    DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3))   &
                                  :: ZW, ZZW1, ZZW2, ZZW4 ! work array
REAL,    DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3))   &
                                  :: ZWEVAP ! sedimentation fluxes
REAL,    DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3))   &
                                  :: ZWSED ! sedimentation fluxes
REAL,    DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3))   &
                                  :: ZWSEDW1,  ZWSEDW2 ! sedimentation speed
REAL,    DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3)) :: ZLBDAR
! Slope parameter of the raindrop  distribution                                  
REAL,    DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3))   &
                                  :: ZZRCT, ZZEVAP, ZMASK
REAL,    DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3))   &
                                  :: ZRAY,   & ! Mean radius
                                     ZNRT,   & ! Number of rain droplets
                                     ZLBC ,  & ! XLBC weighted by sea fraction
                                     ZFSEDC
REAL,    DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2)) :: ZCONC_TMP
REAL,    DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3)) :: ZCONC
REAL,    DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3)) :: ZRRS
!
REAL, DIMENSION(:,:), ALLOCATABLE :: ZSVT  ! Tracer m.r. concentration
!
REAL, DIMENSION(:,:), ALLOCATABLE :: ZVGG, ZDPG  !aerosol velocity [m/s], diffusivity [m2/s]
REAL, DIMENSION(:,:), ALLOCATABLE :: ZRG    !Dust R[�m]
REAL, DIMENSION(:,:), ALLOCATABLE :: ZCOR   !Cunningham correction factor [unitless]
REAL, DIMENSION(:,:), ALLOCATABLE :: ZMASSMIN ! Aerosol mass minimum value
!
REAL, DIMENSION(:), ALLOCATABLE &
               :: ZRHODREF,   &  ! RHO Dry REFerence
                  ZTHT,       &  ! Potential temp
                  ZPABST,     &  ! Pressure [Pa]
                  ZZW,        &  ! Work array
                  ZTEMP,      &  ! Air Temp [K]
                  ZRC,        &  ! Cloud radius [m]
                  ZRCT,       &  ! Cloud water
                  ZRR,        &  ! Rain radius [m]
                  ZNT,        &  ! Rain droplets number 
                  ZRRT,       &  ! Rain water
                  ZMU,ZMUW,   &  ! viscosity aerosol, water [Pa s] 
                  ZFLUX,      &  ! Effective precipitation flux (kg.m-2.s-1)
                  ZCONC1D,    &  ! Weighted droplets concentration
                  ZWLBDC,     &  ! Slope parameter of the droplet  distribution
                  ZGAMMA         ! scavenging coefficient
REAL, DIMENSION(:),   ALLOCATABLE :: ZW1 ! Work arrays 
REAL, DIMENSION(SIZE(XRTMIN))     :: ZRTMIN

INTEGER                           :: JL  ! and PACK intrinsics
!
INTEGER                           :: JKAQ, JSV
!
REAL  :: A0, A1, A2, A3             ! Constants for computing viscocity
INTEGER :: IKE, IJE, IIE, IKB

!  
!-------------------------------------------------------------------------------
!
!*       0.   Initialize work array
!             ---------------------
!
! Compute Effective cloud radius 
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('AER_WET_DEP',0,ZHOOK_HANDLE)

ZRAY(:,:,:)   = 0.
ZLBC(:,:,:)   = XLBC(1)
ZFSEDC(:,:,:) = XFSEDC(1)
ZCONC_TMP(:,:)= XCONC_LAND
ZCONC(:,:,:)  = XCONC_LAND
!
IF (PRESENT(PCCT)) THEN  ! case KHKO, C2R2, C3R5 (two moments schemes)
 ZRAY(:,:,:)  = 3.* PRCT(:,:,:) / (4.*XPI*XRHOLW*PCCT(:,:,:))
 ZRAY(:,:,:)  = ZRAY(:,:,:)**(1./3.) ! Cloud mean radius in m
ELSE IF (PRESENT(PSEA)) THEN ! Case ICE3, RAVE, KESS, ..
ZCONC_TMP(:,:)=PSEA(:,:)*XCONC_SEA+(1.-PSEA(:,:))*XCONC_LAND
DO JK=1,SIZE(PRHODREF,3)
  ZLBC(:,:,JK)   = PSEA(:,:)*XLBC(2)+(1.-PSEA(:,:))*XLBC(1)
  ZFSEDC(:,:,JK) = (PSEA(:,:)*XFSEDC(2)+(1.-PSEA(:,:))*XFSEDC(1))
  ZFSEDC(:,:,JK) = MAX(MIN(XFSEDC(1),XFSEDC(2)),ZFSEDC(:,:,JK))
  ZCONC(:,:,JK)  = (1.-PTOWN(:,:))*ZCONC_TMP(:,:)+PTOWN(:,:)*XCONC_URBAN
  ZRAY(:,:,JK)   = 0.5*((1.-PSEA(:,:))*GAMMA(XNUC+1.0/XALPHAC)/(GAMMA(XNUC)) + &
  PSEA(:,:)*GAMMA(XNUC2+1.0/XALPHAC2)/(GAMMA(XNUC2)))
END DO
ZRAY(:,:,:)      = MAX(1.,ZRAY(:,:,:))
ZLBC(:,:,:)      = MAX(MIN(XLBC(1),XLBC(2)),ZLBC(:,:,:))
ELSE
ZRAY(:,:,:) = 30. ! default value for cloud radius
END IF
!
!ZNRT(:,:,:)  = 0.
!IF (PRESENT(PCRT)) THEN ! case KHKO, C2R2, C3R5
! Transfert Number of rain droplets
  ZNRT(:,:,:)  = PRRT(:,:,:)      
!END IF

!-------------------------------------------------------------------------------
!
!*       1.     COMPUTE THE AEROSOL/CLOUD-RAIN MASS TRANSFER
!              ----------------------------------------------
CALL AER_WET_MASS_TRANSFER_DEP
!-------------------------------------------------------------------------------
!
!*       2.     COMPUTE THE SEDIMENTATION (RS) SOURCE
!              -------------------------------------
!
CALL AER_WET_SEDIMENT_DEP
!
!-------------------------------------------------------------------------------
!!
!!*      3.     COMPUTES THE SLOW WARM PROCESS SOURCES 
!!              --------------------------------------
!!
CALL AER_WET_WARM_DEP

!
!-------------------------------------------------------------------------------
!!*      4.     COMPUTES EVAPORATION PROCESS
!!              ----------------------------
!!
CALL AER_WET_EVAP_DEP   
!

!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('AER_WET_DEP',1,ZHOOK_HANDLE)
!
CONTAINS
!
!
!-------------------------------------------------------------------------------
!
SUBROUTINE AER_WET_MASS_TRANSFER_DEP
!
!*      0. DECLARATIONS
!          ------------
!
IMPLICIT NONE
!
!*       0.2  declaration of local variables
!
!
INTEGER , DIMENSION(SIZE(GCLOUD)) :: I1C,I2C,I3C! Used to replace the COUNT
INTEGER , DIMENSION(SIZE(GRAIN))  :: I1R,I2R,I3R ! Used to replace the COUNT
INTEGER                           :: JL       ! and PACK intrinsics
INTEGER                           :: JKAQ     ! counter for chemistry
!  
!
!  1 Mass transfer Aerosol to cloud (Tost et al., 2006)
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('AER_WET_DEP:AER_WET_MASS_TRANSFER_DEP',0,ZHOOK_HANDLE)
!
GCLOUD(:,:,:) = .FALSE.
GCLOUD(:,:,:) =  PRCS(:,:,:)>XRTMIN(2)
ICLOUD = COUNTJV( GCLOUD(:,:,:),I1C(:),I2C(:),I3C(:))
IF( ICLOUD >= 1 ) THEN   
ALLOCATE(ZSVT(ICLOUD,KMODE*3)) 
ALLOCATE(ZRHODREF(ICLOUD))
ALLOCATE(ZTHT(ICLOUD))
ALLOCATE(ZRC(ICLOUD))
ALLOCATE(ZPABST(ICLOUD))
ALLOCATE(ZRG(ICLOUD,KMODE))    
ALLOCATE(ZTEMP(ICLOUD))
ALLOCATE(ZMU(ICLOUD))
ALLOCATE(ZRCT(ICLOUD))
ALLOCATE(ZVGG(ICLOUD,KMODE))
ALLOCATE(ZDPG(ICLOUD,KMODE))
ALLOCATE(ZGAMMA(ICLOUD))
ALLOCATE(ZW1(ICLOUD))    
ALLOCATE(ZCOR(ICLOUD,KMODE))
ALLOCATE(ZMASSMIN(ICLOUD,KMODE))
ALLOCATE(ZWLBDC(ICLOUD))
ALLOCATE(ZCONC1D(ICLOUD))
ZSVT(:,:) = 0.
DO JL=1,ICLOUD
   DO JKAQ = 1, KMODE
     ZRG(JL,JKAQ) = PRGAER(I1C(JL),I2C(JL),I3C(JL),JKAQ)  
     ZSVT(JL,JKAQ) = PSVT(I1C(JL),I2C(JL),I3C(JL),JKAQ) 
     ZSVT(JL,KMODE+JKAQ) = PSVT(I1C(JL),I2C(JL),I3C(JL),KMODE+JKAQ)
     ZMASSMIN(JL,JKAQ) =  PMASSMIN(I1C(JL),I2C(JL),I3C(JL),JKAQ)
     !
  ZTHT(JL) = PTHT(I1C(JL),I2C(JL),I3C(JL))
  ZRC(JL)  = ZRAY(I1C(JL),I2C(JL),I3C(JL))
  ZPABST(JL) = PPABST(I1C(JL),I2C(JL),I3C(JL))
  ZRCT(JL) = PRCS(I1C(JL),I2C(JL),I3C(JL)) 
  ZRHODREF(JL) =  PRHODREF(I1C(JL),I2C(JL),I3C(JL))   
  ZWLBDC(JL)     = ZLBC(I1C(JL),I2C(JL),I3C(JL))
  ZCONC1D(JL)    = ZCONC(I1C(JL),I2C(JL),I3C(JL))
  ENDDO
END DO

IF (ANY(ZWLBDC(:)/=0.)) THEN ! case one moments
 ! On calcule Rc a partir de M(3) car c'est le seul moment indt de alpha et nu
 ! Rho_air * Rc / (Pi/6 * Rho_eau * Nc) =  M(3) = 1/ (Lambda**3 * rapport des
 ! gamma)
 ZWLBDC(:) = ZWLBDC(:) * ZCONC1D(:) / (ZRHODREF(:) * ZRCT(:))
 ZWLBDC(:) = ZWLBDC(:)**XLBEXC
 ZRC(:) = ZRC(:) / ZWLBDC(:)
END IF

!   
! initialize temperature
 ZTEMP(:)=ZTHT(:)*(ZPABST(:)/XP00)**(XRD/XCPD)

! compute diffusion and gravitation velocity

 CALL AER_VELGRAV(ZRG(:,:),  ZPABST,              &
                  KMODE, ZMU(:), ZVGG(:,:),       &
                  ZDPG(:,:),ZTEMP(:),ZCOR(:,:),   &
                  PDENSITY_AER)
DO JL = 1, ICLOUD
DO JKAQ = 1, KMODE
! Browninan nucleation scavenging (Pruppacher and Klett, 2000, p723) 
  ZGAMMA(JL)  = 1.35 * ZRCT(JL)*ZRHODREF(JL)*1.E-3 * ZDPG(JL,JKAQ) /&
               (ZRC(JL)*ZRC(JL))
  ZW1(JL) = ZSVT(JL,JKAQ) * EXP(-ZGAMMA(JL) * PTSTEP)
  ZW1(JL) = MAX(ZW1(JL), ZMASSMIN(JL,JKAQ))
  ZW1(JL) = MIN(ZW1(JL),ZSVT(JL,JKAQ))
! Aerosol mass in cloud
  ZSVT(JL,KMODE+JKAQ) = ZSVT(JL,KMODE+JKAQ) + ZSVT(JL,JKAQ) - ZW1(JL)
! New aerosol mass 
  ZSVT(JL,JKAQ)  = ZW1(JL)
ENDDO
ENDDO
DO JKAQ = 1, KMODE
! Return in 3D
  PSVT(:,:,:,JKAQ) = &
      UNPACK(ZSVT(:,JKAQ),MASK=GCLOUD(:,:,:),FIELD=PSVT(:,:,:,JKAQ))
  PSVT(:,:,:,KMODE+JKAQ) = &
     UNPACK(ZSVT(:,KMODE+JKAQ),MASK=GCLOUD(:,:,:),FIELD=PSVT(:,:,:,KMODE+JKAQ))
ENDDO

DEALLOCATE(ZSVT) 
DEALLOCATE(ZRHODREF)
DEALLOCATE(ZTHT)
DEALLOCATE(ZRC)
DEALLOCATE(ZPABST)
DEALLOCATE(ZRG)    
DEALLOCATE(ZTEMP)
DEALLOCATE(ZMU)
DEALLOCATE(ZRCT)
DEALLOCATE(ZVGG)
DEALLOCATE(ZDPG)
DEALLOCATE(ZGAMMA)
DEALLOCATE(ZW1)    
DEALLOCATE(ZCOR)
DEALLOCATE(ZMASSMIN)
DEALLOCATE(ZWLBDC)
DEALLOCATE(ZCONC1D)
END IF
!
!  2 Mass transfer Aerosol to Rain (Seinfeld and Pandis, 1998, Tost et al., 2006)
!
GRAIN(:,:,:) = .FALSE.
GRAIN(:,:,:) =  PRRS(:,:,:)>XRTMIN(3)
IRAIN = COUNTJV( GRAIN(:,:,:),I1R(:),I2R(:),I3R(:))
IF( IRAIN >= 1 ) THEN   
! 
 ALLOCATE(ZRRT(IRAIN))
 ALLOCATE(ZSVT(IRAIN,3*KMODE)) 
 ALLOCATE(ZRHODREF(IRAIN))
 ALLOCATE(ZTHT(IRAIN))
 ALLOCATE(ZRR(IRAIN))
 ALLOCATE(ZNT(IRAIN))
 ALLOCATE(ZPABST(IRAIN))
 ALLOCATE(ZRG(IRAIN,KMODE))    
 ALLOCATE(ZCOR(IRAIN,KMODE))
 ALLOCATE(ZTEMP(IRAIN))
 ALLOCATE(ZMU(IRAIN))
 ALLOCATE(ZVGG(IRAIN,KMODE))
 ALLOCATE(ZDPG(IRAIN,KMODE))
 ALLOCATE(ZMUW(IRAIN))
 ALLOCATE(ZEFC(IRAIN,KMODE))
 ALLOCATE(ZW1(IRAIN))
 ALLOCATE(ZFLUX(IRAIN))
 ALLOCATE(ZGAMMA(IRAIN))
 ALLOCATE(ZMASSMIN(IRAIN,KMODE))

 ZSVT(:,:) = 0.
 DO JL=1,IRAIN    
   DO JKAQ = 1, KMODE
    ZRG(JL,JKAQ)   = PRGAER(I1R(JL),I2R(JL),I3R(JL),JKAQ )  
    ZSVT(JL,JKAQ)  = PSVT(I1R(JL),I2R(JL),I3R(JL),JKAQ)
    ZSVT(JL,KMODE*2+JKAQ) = PSVT(I1R(JL),I2R(JL),I3R(JL),KMODE*2+JKAQ) 
    ZMASSMIN(JL,JKAQ) =  PMASSMIN(I1R(JL),I2R(JL),I3R(JL),JKAQ)
     !
   ZTHT(JL) = PTHT(I1R(JL),I2R(JL),I3R(JL))
   ZPABST(JL) = PPABST(I1R(JL),I2R(JL),I3R(JL))
   ZRRT(JL) = PRRS(I1R(JL),I2R(JL),I3R(JL)) 
   ZRHODREF(JL) =  PRHODREF(I1R(JL),I2R(JL),I3R(JL))   
   ZNT(JL) =  ZNRT(I1R(JL),I2R(JL),I3R(JL))   
   ENDDO
 ENDDO

!
CALL AER_WET_EFFIC_DEP

! Compute scavenging coefficient
ZFLUX(:) = 0.
ZRRT(:) = MAX(ZRRT(:), 0.)
! Effective precipitation flux (kg.m-2.s-1)
ZFLUX(:) = XFSEDR  * ZRRT(:)**(XEXSEDR )   &
           * ZRHODREF(:)**(XEXSEDR-XCEXVT)  
ZFLUX(:) =  MAX(ZFLUX(:), 0.)

!IF (ALL(ZNT(:) == 0.)) THEN ! case one moments
!Number concentration NT=No/lbda   p. 415 Jacobson
!4/3 *pi *r�*NT*rho_eau(kg/m3) =rho(lwc)=rho(air)* qc(kg/kg)
ZNT (:) = XCCR/(XLBR*( ZRHODREF(:)* ZRRT(:) )**XLBEXR)
!END IF

ZRR(:) =  (ZRRT(:)*ZRHODREF(:)/(XRHOLW*ZNT(:)*4./3.*XPI))**(1./3.)

DO JL=1,IRAIN    
DO JKAQ = 1, KMODE
  ! Tost et al, 2006
  ZGAMMA(JL)  = 0.75 * ZEFC(JL,JKAQ) * ZFLUX(JL) / (ZRR(JL)*1E3)
  ZW1(JL) = ZSVT(JL,JKAQ) * EXP(-ZGAMMA(JL) * PTSTEP)
  ZW1(JL) = MAX(ZW1(JL), ZMASSMIN(JL,JKAQ))
  ZW1(JL) = MIN(ZW1(JL),ZSVT(JL,JKAQ))
  ! Aerosol mass in rain
  ZSVT(JL,KMODE*2+JKAQ) = ZSVT(JL,KMODE*2+JKAQ) + ZSVT(JL,JKAQ) - ZW1(JL)
  ! New aerosol mass 
  ZSVT(JL,JKAQ) =  ZW1(JL)
ENDDO
ENDDO

DO JKAQ = 1, KMODE
  ! Return to 3D
  PSVT(:,:,:,JKAQ) = &
    UNPACK(ZSVT(:,JKAQ),MASK=GRAIN(:,:,:),FIELD=PSVT(:,:,:,JKAQ))
  PSVT(:,:,:,KMODE*2+JKAQ) = &
    UNPACK(ZSVT(:,KMODE*2+JKAQ),MASK=GRAIN(:,:,:),FIELD=PSVT(:,:,:,KMODE*2+JKAQ))
ENDDO

 DEALLOCATE(ZRRT)
 DEALLOCATE(ZSVT) 
 DEALLOCATE(ZRHODREF)
 DEALLOCATE(ZTHT)
 DEALLOCATE(ZRR)
 DEALLOCATE(ZNT)
 DEALLOCATE(ZPABST)
 DEALLOCATE(ZRG)    
 DEALLOCATE(ZCOR)
 DEALLOCATE(ZTEMP)
 DEALLOCATE(ZMU)
 DEALLOCATE(ZVGG)
 DEALLOCATE(ZDPG)
 DEALLOCATE(ZMUW)
 DEALLOCATE(ZEFC)
 DEALLOCATE(ZW1)
 DEALLOCATE(ZFLUX)
 DEALLOCATE(ZGAMMA)
 DEALLOCATE(ZMASSMIN)
 END IF 
IF (LHOOK) CALL DR_HOOK('AER_WET_DEP:AER_WET_MASS_TRANSFER_DEP',1,ZHOOK_HANDLE)
END SUBROUTINE AER_WET_MASS_TRANSFER_DEP
!
SUBROUTINE AER_WET_SEDIMENT_DEP
!
!*         Sedimentation of aerosol in rain droplets
!
!*      0. DECLARATIONS
!          ------------
!
IMPLICIT NONE
!
!*        declaration of local variables
!
!
REAL                              :: P1,P2,ZQP,ZH
INTEGER                           :: JL       ! and PACK intrinsics
INTEGER                           :: JKAQ     ! counter for acquous aerosols
!  
!-------------------------------------------------------------------------------
!
!*         Time splitting initialization
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('AER_WET_DEP:AER_WET_SEDIMENT_DEP',0,ZHOOK_HANDLE)
ZTSPLITR = PTSTEP / FLOAT(KSPLITR)
!
ZWSEDW1(:,:,:)=0.
ZWSEDW2(:,:,:)=0.
ZWSED(:,:,:) = 0.
ZW(:,:,:)=0.
ZRRS(:,:,:) = MAX(PTSTEP*PRRS(:,:,:), 1E-80)
IKE = SIZE(PRRS,3)-1
IJE = SIZE(PRRS,2)
IIE = SIZE(PRRS,1)
IKB = 2
!
DO JK = 1 , SIZE(PZZ,3)-1
  ZW(:,:,JK) =ZTSPLITR/(( PZZ(:,:,JK+1)-PZZ(:,:,JK) ))
END DO
!WHERE (ZRRS(:,:,:)<=XRTMIN(3))    
!  ZW(:,:,:)=0.
!END WHERE
! 
! Flux mass aerosol in rain droplets = 
! Flux mass rain water * Mass aerosol in rain / Mass rain water
!print*,'XEXSEDR=',XEXSEDR,'*****','XCEXVT=',XCEXVT,'***','KSPLITR=',KSPLITR

DO JKAQ = 1,KMODE
   DO JK = IKE,IKB, -1
     DO JJ = 1, IJE
       DO JI = 1, IIE
         !estimation of q' taking into account incomming ZWSED
           ZQP=ZWSED(JI,JJ,JK+1)*ZW(JI,JJ,JK)
         !calculation of w
         IF ( PRRS(JI,JJ,JK) > ZRTMIN(3) ) THEN
           ZWSEDW1 (JI,JJ,JK)= XFSEDR *PRRS(JI,JJ,JK)**(XEXSEDR-1)* &
           PRHODREF(JI,JJ,JK)**(XEXSEDR-XCEXVT-1)
         ENDIF
         IF ( ZQP > ZRTMIN(3) ) THEN
           ZWSEDW2 (JI,JJ,JK)= XFSEDR *(ZQP)**(XEXSEDR-1)* &
           PRHODREF(JI,JJ,JK)**(XEXSEDR-XCEXVT-1)
         ENDIF
       ENDDO
       DO JI = 1, IIE
         ZH=PZZ(JI,JJ,JK+1) - PZZ(JI,JJ,JK)
         P1 = MIN(1., ZWSEDW1(JI,JJ,JK) * PTSTEP / ZH )
         IF (ZWSEDW2(JI,JJ,JK) /= 0.) THEN
           P2 = MAX(0.,1 -  ZH &
         & / (PTSTEP*ZWSEDW2(JI,JJ,JK)) )
         ELSE
           P2 = 0.
         ENDIF
         ZWSED (JI,JJ,JK)=P1*PRHODREF(JI,JJ,JK)*&
         &ZH*PRRS(JI,JJ,JK)&
         &/PTSTEP+ P2 * ZWSED (JI,JJ,JK+1)
!
         ZWSED (JI,JJ,JK)=ZWSED (JI,JJ,JK)*PSVT(JI,JJ,JK,KMODE*2+JKAQ)
!
       ENDDO
     ENDDO
   ENDDO
ENDDO
!
DO JKAQ = 1,KMODE
   DO JK = IKB, IKE
    PSVT(:,:,JK,KMODE*2+JKAQ)=  PSVT(:,:,JK,KMODE*2+JKAQ) + &
                                ZW(:,:,JK)*(ZWSED(:,:,JK+1)-ZWSED(:,:,JK))
   ! Aerosol mass in rain droplets need to be positive
    PSVT(:,:,JK,KMODE*2+JKAQ)= MAX(PSVT(:,:,JK,KMODE*2+JKAQ), 0.)
   END DO
END DO
  
  IF (LHOOK) CALL DR_HOOK('AER_WET_DEP:AER_WET_SEDIMENT_DEP',1,ZHOOK_HANDLE)
!  
  END SUBROUTINE AER_WET_SEDIMENT_DEP
!
!-------------------------------------------------------------------------------
!
  SUBROUTINE AER_WET_WARM_DEP
!
!*      0. DECLARATIONS
!
IMPLICIT NONE
!-------------------------------------------------------------------------------
!
!*       1.    compute the autoconversion of r_c for r_r production: RCAUTR
!    
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('AER_WET_DEP:AER_WET_WARM_DEP',0,ZHOOK_HANDLE)
ZZW4(:,:,:)=0.0
! to be sure no division by zero in case of ZZRCT = 0.
ZZRCT(:,:,:) = PRCT(:,:,:)
ZZRCT(:,:,:) = MAX(ZZRCT(:,:,:), XRTMIN(2)/2.)

WHERE( (ZZRCT(:,:,:)>XRTMIN(2)) .AND. (PRCS(:,:,:)>0.0 ) ) 
      ZZW4(:,:,:) = MIN( PRCS(:,:,:),XTIMAUTC* &
               MAX((ZZRCT(:,:,:)-XCRIAUTC/ PRHODREF(:,:,:)),0.0))
END WHERE 

DO JKAQ = 1,KMODE
    ZZW2(:,:,:) =0.0   
    ZZW2(:,:,:)=ZZW4(:,:,:) * PSVT(:,:,:,KMODE+JKAQ)/ZZRCT(:,:,:) * PTSTEP
    ZZW2(:,:,:) = MAX(MIN( ZZW2(:,:,:), PSVT(:,:,:,KMODE+JKAQ)),0.0)

! For rain - Increase the aerosol conc in rain
    PSVT(:,:,:,KMODE*2+JKAQ) = &
                   PSVT(:,:,:,KMODE*2+JKAQ) +  ZZW2(:,:,:)
! For Cloud Decrease the aerosol conc in cloud
    PSVT(:,:,:,KMODE+JKAQ) = &
                   PSVT(:,:,:,KMODE+JKAQ) - ZZW2(:,:,:)
ENDDO
! 
!*       2.    compute the accretion of r_c for r_r production: RCACCR
!
ZZW4(:,:,:)=0.0
ZLBDAR(:,:,:)=0.0
WHERE ( (ZZRCT(:,:,:)>XRTMIN(2)) .AND. (PRRT(:,:,:)>XRTMIN(3))    &
                           .AND. (PRCS(:,:,:)> 0.0 ) )      
     ZLBDAR(:,:,:)  = XLBR*( PRHODREF(:,:,:)* PRRT(:,:,:) )**XLBEXR
     ZZW4(:,:,:) = MIN( PRCS(:,:,:),XFCACCR * ZZRCT(:,:,:) &
                                    * ZLBDAR(:,:,:)**XEXCACCR      &
                                    * PRHODREF(:,:,:)**(-XCEXVT) )
END WHERE                           
!
DO JKAQ = 1,KMODE
  ZZW2(:,:,:)=0.0
  ZZW2(:,:,:)=ZZW4(:,:,:) * PSVT(:,:,:,KMODE+JKAQ)/ZZRCT(:,:,:) * PTSTEP
  ZZW2(:,:,:) = MAX(MIN(ZZW2(:,:,:),PSVT(:,:,:,KMODE+JKAQ)), 0.0)

!*       3.    compute the new acquous aerosol mass
!
! For rain - Increase the aerosol conc in rain
   PSVT(:,:,:,KMODE*2+JKAQ) = PSVT(:,:,:,KMODE*2+JKAQ) + ZZW2(:,:,:)
! For Cloud Decrease the aerosol conc in cloud
   PSVT(:,:,:,KMODE+JKAQ)   = PSVT(:,:,:,KMODE+JKAQ) - ZZW2(:,:,:)
ENDDO
 IF (LHOOK) CALL DR_HOOK('AER_WET_DEP:AER_WET_WARM_DEP',1,ZHOOK_HANDLE)                
 END SUBROUTINE AER_WET_WARM_DEP
!---------------------------------------------------------------------------------------
  SUBROUTINE AER_WET_EVAP_DEP
!
!*            COMPUTES THE EVAPORATION OF CLOUD-RAIN FOR THE 
!*             RE-RELEASE OF AER INTO THE ENVIRONMENT
!                --------------------------------------
!
  
!
!*      0. DECLARATIONS
!          ------------
!
IMPLICIT NONE
!*         declaration of local variables
!
INTEGER    :: JKAQ     ! counter for aerosols

!*       1.    compute the evaporation of r_r: RREVAV 
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('AER_WET_DEP:AER_WET_EVAP_DEP',0,ZHOOK_HANDLE)         
!When partial reevaporation of precip takes place, the fraction of 
!tracer precipitating form above is reevaporated is equal to
!half of the evaporation rate of water
!
! Rain water evaporated during PTSTEP in kg/kg
ZZEVAP(:,:,:) = PEVAP3D(:,:,:) * PTSTEP  
! Fraction of rain water evaporated 
! at this stage (bulk), we consider that the flux of evaporated aerosol
! is a ratio of the evaporated rain water. 
! It will interested to calculate with a two moment scheme (C2R2 or C3R5)
! the complete evaporation of rain droplet to use it for the compuation
! of the evaporated aerosol flux.
ZWEVAP(:,:,:)=0.0
WHERE( PRRT(:,:,:) .GT. XRTMIN(3) )
    ZWEVAP(:,:,:) = ZZEVAP(:,:,:)/(PRRT(:,:,:))  
END WHERE
ZWEVAP(:,:,:)=MIN(ZWEVAP(:,:,:),1.0)
ZWEVAP(:,:,:)=MAX(ZWEVAP(:,:,:),0.0)

!*       2.    compute the mask of r_c evaporation : all cloud is evaporated
!              no partial cloud evaporation at this stage
ZMASK(:,:,:) = 0.
WHERE( PRCS(:,:,:) .LT. XRTMIN(2) )
   ZMASK(:,:,:) = 1.
END WHERE
!
!
DO JKAQ = 1,KMODE       
   ZZW1(:,:,:) = ZMASK(:,:,:)*PSVT(:,:,:,KMODE+JKAQ)

   ZZW2(:,:,:) = ZWEVAP(:,:,:)*PSVT(:,:,:,KMODE*2+JKAQ)

!     3.   New dry aerosol mass
!
   PSVT(:,:,:,JKAQ) =  PSVT(:,:,:,JKAQ) + ZZW2(:,:,:) + ZZW1(:,:,:)
!     4.  New cloud aerosol mass
!
   PSVT(:,:,:,KMODE+JKAQ) = PSVT(:,:,:,KMODE+JKAQ) - ZZW1(:,:,:)
       
!     5.  New rain aerosol mass
!     
   PSVT(:,:,:,KMODE*2+JKAQ) = PSVT(:,:,:,KMODE*2+JKAQ) - ZZW2(:,:,:)
END DO      
!   
  IF (LHOOK) CALL DR_HOOK('AER_WET_DEP:AER_WET_EVAP_DEP',1,ZHOOK_HANDLE)
  END SUBROUTINE AER_WET_EVAP_DEP
!---------------------------------------------------------------------------------------
  SUBROUTINE AER_WET_EFFIC_DEP
!
!*            COMPUTES THE EFFICIENCY FACTOR 
!             ------------------------------
!
!
!*      0. DECLARATIONS
!          ------------
!
IMPLICIT NONE
!
!
!*       1.     COMPUTES THE EFFICIENCY FACTOR 
!                 --------------------------------------
!
!*       1.1   compute gravitational velocities
!   
!initialize
 REAL(KIND=JPRB) :: ZHOOK_HANDLE
 IF (LHOOK) CALL DR_HOOK('AER_WET_DEP:AER_WET_EFFIC_DEP',0,ZHOOK_HANDLE)

 ZTEMP(:)=ZTHT(:)*(ZPABST(:)/XP00)**(XRD/XCPD)
 ZTEMP(:)=MAX(ZTEMP(:),1.e-12)

 CALL AER_VELGRAV(ZRG(:,:), ZPABST, KMODE,  &
                  ZMU(:), ZVGG(:,:),        &
                  ZDPG(:,:),ZTEMP(:),       &
                  ZCOR(:,:), PDENSITY_AER)

!  Above gives mu (ZMU), v(aerosol)(PVGG, m/s), diffusion (ZDPG, m2/s)
!
!*      1.2   Compute Water Viscocity in kg/m/s Prup. & Klett, p.95
!
!
 A0=1.76
 A1=-5.5721e-2
 A2=-1.3943e-3
 A3=-4.3015e-5
 ZMUW(:)=A0*EXP(A1*(ZTEMP(:)-273.15) &
        +A2*(ZTEMP(:)-273.15) + A3*(ZTEMP(:)-273.15))*1.e-3
      
 A1=-3.5254e-2
 A2=4.7163e-4
 A3=-6.0667e-6

 WHERE  (ZTEMP(:)>273.15) 
         ZMUW(:)=A0*EXP(A1*(ZTEMP(:)-273.15) &
        +A2*(ZTEMP(:)-273.15) + A3*(ZTEMP(:)-273.15))*1.e-3
 END WHERE
 ZMUW(:)=MAX(ZMUW(:),1.e-12)

!
!*       1.3   compute efficiency factor
!
! This gives aerosol collection efficiency by calculating Reynolds number
! schmidt number, stokes number, etc
 CALL AER_EFFIC_DEP(ZRG(:,:), ZVGG(:,:),   & !aerosol radius/velocity
                ZRHODREF(:),               & !Air density
                ZMUW(:), ZMU(:),           & !mu water/air
                ZDPG(:,:), ZEFC(:,:),      & !diffusivity, efficiency
                ZRRT(:), KMODE,            & !Rain water, nb aerosols modes
                ZTEMP(:),ZCOR(:,:),        & ! Temperature, Cunnimgham coeff
                PDENSITY_AER  )              ! aerosol density
!
IF (LHOOK) CALL DR_HOOK('AER_WET_DEP:AER_WET_EFFIC_DEP',1,ZHOOK_HANDLE)
END SUBROUTINE AER_WET_EFFIC_DEP

!
!-------------------------------------------------------------------------------
!
!
  FUNCTION COUNTJV(LTAB,I1,I2,I3) RESULT(IC)
!
!*      0. DECLARATIONS
!          ------------
!
IMPLICIT NONE
!
!*        declaration of local variables
!
!
LOGICAL, DIMENSION(:,:,:) :: LTAB ! Mask
INTEGER, DIMENSION(:) :: I1,I2,I3 ! Used to replace the COUNT and PACK
INTEGER :: JI,JJ,JK,IC
!  
!-------------------------------------------------------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('AER_WET_DEP:COUNTJV',0,ZHOOK_HANDLE)
I1(:)=0 
I2(:)=0 
I3(:)=0                
!ktk
IC = 0
DO JK = 1,SIZE(LTAB,3)
  DO JJ = 1,SIZE(LTAB,2)
    DO JI = 1,SIZE(LTAB,1)
      IF( LTAB(JI,JJ,JK) ) THEN
        IC = IC +1
        I1(IC) = JI
        I2(IC) = JJ
        I3(IC) = JK
      END IF
    END DO
  END DO
END DO
!
IF (LHOOK) CALL DR_HOOK('AER_WET_DEP:COUNTJV',1,ZHOOK_HANDLE)
END FUNCTION COUNTJV
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE AER_WET_DEP
