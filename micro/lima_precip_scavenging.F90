!MNH_LIC Copyright 2013-2024 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!########################################################################
   SUBROUTINE LIMA_PRECIP_SCAVENGING (TNSV, D, CST, BUCONF, TBUDGETS, KBUDGETS, &
                                      HCLOUD, HDCONF, KLUOUT, KTCOUNT, PTSTEP,    &
                                      PRRT, PRHODREF, PRHODJ, PZZ,        &
                                      PPABST, PTHT, PSVT, PRSVS, PINPAP )
!########################################################################x
!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the total number
!!      below-cloud scavenging rate. 
!!
!!
!!**  METHOD
!!    ------
!!      We assume a generalized gamma distribution law for the raindrop.
!!      The aerosols particles distribution follows a log-normal law.
!!      First, we have to compute the Collision Efficiency, which takes 
!!      account of the three most important wet removal mechanism : 
!!      Brownian diffusion, interception and inertial impaction. 
!!      It is a function of several number (like Reynolds, Schmidt 
!!      or Stokes number for instance). Consequently,
!!      we need first to calculate these numbers. 
!!
!!      Then the scavenging coefficient is deduced from the integration 
!!      of the droplet size distribution, the falling velocity of 
!!      raindrop and aerosol, their diameter, and the collision 
!!      (or collection) efficiency, over the spectrum of droplet
!!      diameters.
!!
!!      The total scavenging rate of aerosol is computed from the 
!!      integration, over all the spectrum of particles aerosols 
!!      diameters, of the scavenging coefficient.
!! 
!!
!!    EXTERNAL
!!    --------
!!      Subroutine SCAV_MASS_SEDIMENTATION
!!
!!      Function COLL_EFFIC    : computes the collision efficiency
!!
!!      Function CONTJV   |
!!      Function GAUHER   |
!!      Function GAULAG   |-> in lima_functions.f90
!!      Function GAMMLN   |
!!     
!!
!!   REFERENCES
!!   ----------
!!   Seinfeld and Pandis
!!   Andronache
!!   
!!   AUTHOR
!!   ------
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!      S.    Berthet    * Laboratoire d'Aerologie*
!!      B.    Vié        * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original             ??/??/13 
!!
!  P. Wautelet 28/05/2018: corrected truncated integer division (3/2 -> 1.5)
!  P. Wautelet 26/04/2019: replace non-standard FLOAT function by REAL function
!  P. Wautelet 28/05/2019: move COUNTJV function to tools.f90
!  P. Wautelet    03/2020: use the new data structures and subroutines for budgets
!  P. Wautelet 03/06/2020: bugfix: correct array starts for PSVT and PRSVS
!  P. Wautelet 11/02/2021: bugfix: ZRTMIN was of wrong size (replaced by a scalar)
!  P. Wautelet 11/02/2021: budgets: add missing term SCAV for NSV_LIMA_SCAVMASS budget
!-------------------------------------------------------------------------------
!
!*                  0.DECLARATIONS          
!                   --------------
!
USE MODD_DIMPHYEX,        ONLY: DIMPHYEX_T
USE MODD_BUDGET,          ONLY: TBUDGETDATA_PTR, TBUDGETCONF_T, NBUDGET_SV1
USE MODD_CST,             ONLY: CST_T
USE MODD_NSV,             ONLY: NSV_T
USE MODD_PARAM_LIMA,      ONLY: NMOD_IFN, NSPECIE, XFRAC,                         &
                                XMDIAM_IFN, XSIGMA_IFN, XRHO_IFN,                 &
                                NMOD_CCN, XR_MEAN_CCN, XLOGSIG_CCN, XRHO_CCN,     &
                                XALPHAR, XNUR,                                    &
                                LAERO_MASS, NDIAMR, NDIAMP, XT0SCAV, XTREF,       &
                                XMUA0, XT_SUTH_A, XMFPA0, XVISCW, XRHO00,         &
                                XRTMIN, XCTMIN
USE MODD_PARAM_LIMA_WARM, ONLY: XCR, XDR

USE MODE_TOOLS,           only: COUNTJV

USE MODI_GAMMA
USE MODE_LIMA_FUNCTIONS, ONLY: GAUHER, GAULAG
USE YOMHOOK, ONLY:LHOOK, DR_HOOK, JPHOOK

IMPLICIT NONE
!
!*                 0.1 declarations of dummy arguments :
!
TYPE(NSV_T),              INTENT(IN)    :: TNSV
TYPE(DIMPHYEX_T),         INTENT(IN)    :: D
TYPE(CST_T),              INTENT(IN)    :: CST
TYPE(TBUDGETCONF_T),      INTENT(IN)    :: BUCONF
TYPE(TBUDGETDATA_PTR), DIMENSION(KBUDGETS), INTENT(INOUT) :: TBUDGETS
INTEGER,                  INTENT(IN)    :: KBUDGETS
!
CHARACTER(LEN=4),       INTENT(IN)    :: HCLOUD   ! cloud paramerization
CHARACTER(LEN=5),       INTENT(IN)    :: HDCONF   ! CCONF from MODD_CONF
INTEGER,                INTENT(IN)    :: KLUOUT   ! unit for output listing
INTEGER,                INTENT(IN)    :: KTCOUNT  ! iteration count
REAL,                   INTENT(IN)    :: PTSTEP   ! Double timestep except 
                                                  ! for the first time step
!
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PRRT     ! Rain mixing ratio at t
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PRHODREF ! Air Density [kg/m**3]
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PRHODJ   ! Dry Density [kg]
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PZZ      ! Altitude
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PPABST   ! Absolute pressure at t
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PTHT     ! Theta at time t 
!
REAL, DIMENSION(D%NIJT,D%NKT,TNSV%NSV_LIMA), INTENT(IN)    :: PSVT   ! Particle Concentration [/m**3]
REAL, DIMENSION(D%NIJT,D%NKT,TNSV%NSV_LIMA), INTENT(INOUT) :: PRSVS  ! Total Number Scavenging Rate
!
REAL, DIMENSION(D%NIJT),   INTENT(INOUT) :: PINPAP
!
!*       0.2   Declarations of local variables :
!
INTEGER :: IIB           !  Define the domain where is 
INTEGER :: IIE           !  the microphysical sources have to be computed
INTEGER :: IJB           ! 
INTEGER :: IJE           !
INTEGER :: IKB           ! 
INTEGER :: IKE           !
!
INTEGER :: ISV               ! CCN or IFN mode 
INTEGER :: J1, J2, IMOD
!
LOGICAL, DIMENSION(D%NIJT,D%NKT) &
                                 :: GRAIN,  &! Test where rain is present
                                    GSCAV    ! Test where rain is present
INTEGER , DIMENSION(SIZE(GSCAV)) :: I1,I2,I3 ! Used to replace the COUNT
INTEGER                          :: IL       ! and PACK intrinsics
INTEGER                          :: ISCAV
!
REAL                     :: ZDENS_RATIO, & !density ratio 
                            ZNUM,        & !PNU-1.               
                            ZSHAPE_FACTOR
!
REAL,    DIMENSION(D%NIJT,D%NKT)  :: ZCRT   ! cloud droplet conc.
!
REAL, DIMENSION(:), ALLOCATABLE :: ZLAMBDAR,      &  !slope parameter of the 
                                                     ! generalized Gamma 
                                                     !distribution law for the 
                                                     !raindrop
                                   ZVISC_RATIO,   &  !viscosity ratio 
                                   ZMFPA,         &  !Mean Free Path
                                   ZRHODREF,      &  !Air Density [kg/m**3]
                                   ZVISCA,        &  !Viscosity of Air [kg/(m*s)]
                                   ZT,            &  !Absolute Temperature
                                   ZPABST,        &
                                   ZRRT,          &
                                   ZCONCP,        &
                                   ZCONCR,        &
                                   ZTOT_SCAV_RATE1D,&
                                   ZTOT_MASS_RATE1D,&
                                   ZMEAN_SCAV_COEF1D
!
REAL, DIMENSION(:,:), ALLOCATABLE :: &
                      ZVOLDR,        &  !Mean volumic Raindrop diameter [m]
                      ZBC_SCAV_COEF2D, &
                      ZCUNSLIP,      &  !CUnningham SLIP correction factor 
                      ZST_STAR,      &  !critical Stokes number for impaction
                      ZSC,           &  !aerosol particle Schmidt number
                      ZRE,           &  !raindrop Reynolds number (for radius)
                      ZFVELR,        &  !Falling VELocity of the Raindrop
                      ZRELT,         &  !RELaxation Time of the particle [s]
                      ZDIFF             !Particle Diffusivity
!
REAL, DIMENSION(NDIAMP) :: ZVOLDP,          & !Mean volumic diameter [m]
                           ZABSCISSP,       & !Aerosol Abscisses
                           ZWEIGHTP           !Aerosol Weights
REAL, DIMENSION(NDIAMR) :: ZABSCISSR,       & !Raindrop Abscisses
                           ZWEIGHTR           !Raindrop Weights
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZCOL_EF,     &! Collision efficiency
                                       ZSIZE_RATIO, &! Size Ratio
                                       ZST           ! Stokes number
!
!
REAL, DIMENSION(D%NIJT,D%NKT) &
                                    :: ZMEAN_SCAV_COEF, & !Mean Scavenging 
                                                          ! Coefficient
                                       ZTOT_SCAV_RATE,  & !Total Number 
                                                          ! Scavenging Rate
                                       ZTOT_MASS_RATE     !Total Mass
                                                          ! Scavenging Rate
REAL, DIMENSION(D%NIJT,D%NKT,NDIAMP) &
                                    ::ZBC_SCAV_COEF  !Scavenging Coefficient
REAL, DIMENSION(:), ALLOCATABLE :: ZKNUDSEN ! Knuudsen number
!
! Opt. BVIE
REAL, DIMENSION(D%NIJT,D%NKT)        &
                   :: ZT_3D, ZCONCR_3D, ZVISCA_3D, ZMFPA_3D,               &
                      ZVISC_RATIO_3D, ZLAMBDAR_3D, ZFACTOR_3D
REAL, DIMENSION(D%NIJT,D%NKT,NDIAMP) &
                   :: ZVOLDR_3D, ZVOLDR_3D_INV, ZVOLDR_3D_POW,             &
                      ZFVELR_3D, ZRE_3D, ZRE_3D_SQRT, ZST_STAR_3D
REAL, DIMENSION(:), ALLOCATABLE   :: ZFACTOR 
REAL, DIMENSION(:,:), ALLOCATABLE ::     &
                      ZRE_SQRT,          &  ! SQRT of raindrop Reynolds number
                      ZRE_INV,           &  ! INV of raindrop Reynolds number
                      ZSC_INV,           &  ! INV of aerosol particle Schmidt number
                      ZSC_SQRT,          &  ! SQRT of aerosol particle Schmidt number
                      ZSC_3SQRT,         &  ! aerosol particle Schmidt number**(1./3.)
                      ZVOLDR_POW,        &  ! Mean volumic Raindrop diameter [m] **(2+ZDR)
                      ZVOLDR_INV            ! INV of Mean volumic Raindrop diameter [m]
REAL               :: ZDENS_RATIO_SQRT 
INTEGER :: ISV_VAR, INM, IJM
INTEGER :: IDX
REAL :: ZMDIAMP 
REAL :: ZSIGMAP  
REAL :: ZRHOP   
REAL :: ZFRACP
!
INTEGER :: ISV_LIMA_NR
INTEGER :: ISV_LIMA_SCAVMASS
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
!------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('LIMA_PRECIP_SCAVENGING', 0, ZHOOK_HANDLE)
ISV_LIMA_NR       = TNSV%NSV_LIMA_NR       - TNSV%NSV_LIMA_BEG + 1
ISV_LIMA_SCAVMASS = TNSV%NSV_LIMA_SCAVMASS - TNSV%NSV_LIMA_BEG + 1

IF ( BUCONF%LBUDGET_SV ) then
  DO IL = 1, NMOD_CCN
    IDX = TNSV%NSV_LIMA_CCN_FREE - 1 + IL
    CALL TBUDGETS(NBUDGET_SV1 - 1 + IDX)%PTR%INIT_PHY(D, 'SCAV', PRSVS( :, :, IDX) )
  END DO
  DO IL = 1, NMOD_IFN
    IDX = TNSV%NSV_LIMA_IFN_FREE - 1 + IL
    CALL TBUDGETS(NBUDGET_SV1 - 1 + IDX)%PTR%INIT_PHY(D, 'SCAV', PRSVS( :, :, IDX) )
  END DO
  IF ( LAERO_MASS ) then
     CALL TBUDGETS(NBUDGET_SV1 - 1 + TNSV%NSV_LIMA_SCAVMASS)%PTR%INIT_PHY(D, 'SCAV', &
          PRSVS( :, :, TNSV%NSV_LIMA_SCAVMASS) )
  END IF
END IF
!
!*       1.     PRELIMINARY COMPUTATIONS
!               ------------------------
!
!
! ZCRT
ZCRT(:,:)=PSVT(:,:,ISV_LIMA_NR)
!
! Rain mask 
GRAIN(:,:) = .FALSE.
GRAIN(D%NIJB:D%NIJE,D%NKB:D%NKE) = (PRRT(D%NIJB:D%NIJE,D%NKB:D%NKE)>XRTMIN(3) &
                            .AND. ZCRT(D%NIJB:D%NIJE,D%NKB:D%NKE)>XCTMIN(3) )
!
! Initialize the total mass scavenging rate if LAERO_MASS=T
IF (LAERO_MASS) ZTOT_MASS_RATE(:,:)  = 0.
!
! Quadrature method: compute absissae and weights
CALL GAUHER(ZABSCISSP,ZWEIGHTP,NDIAMP)
ZNUM = XNUR-1.0E0
CALL GAULAG(ZABSCISSR,ZWEIGHTR,NDIAMR,ZNUM)
!
!
!------------------------------------------------------------------------------
!
!
!*       2.     NUMERICAL OPTIMIZATION
!               ----------------------
!
!
! Optimization : compute in advance parameters depending on rain particles and
! environment conditions only, to avoid multiple identical computations in loops
!
!
ZSHAPE_FACTOR = GAMMA_X0D(XNUR+3./XALPHAR)/GAMMA_X0D(XNUR) 
!
WHERE ( GRAIN(:,:) )
   !
   ZT_3D(:,:)      = PTHT(:,:) * ( PPABST(:,:)/CST%XP00 )**(CST%XRD/CST%XCPD)
   ZCONCR_3D(:,:)  = ZCRT(:,:) * PRHODREF(:,:)                    ![/m3]
   ! Sutherland law for viscosity of air
   ZVISCA_3D(:,:)  = XMUA0*(ZT_3D(:,:)/XTREF)**1.5*(XTREF+XT_SUTH_A) &
                                                        /(XT_SUTH_A+ZT_3D(:,:))
   ! Air mean free path
   ZMFPA_3D(:,:)   = XMFPA0*(CST%XP00*ZT_3D(:,:))/(PPABST(:,:)*XT0SCAV)           
   ! Viscosity ratio
   ZVISC_RATIO_3D(:,:) = ZVISCA_3D(:,:)/XVISCW   !!!!! inversé par rapport à orig. !
   ! Rain drops parameters
   ZLAMBDAR_3D(:,:) = ( ((CST%XPI/6.)*ZSHAPE_FACTOR*CST%XRHOLW*ZCONCR_3D(:,:))  &
                          /(PRHODREF(:,:)*PRRT(:,:)) )**(1./3.)         ![/m]
   ZFACTOR_3D(:,:) = CST%XPI*0.25*ZCONCR_3D(:,:)*XCR*(XRHO00/PRHODREF(:,:))**(0.4)
   !
END WHERE
!
DO J2=1,NDIAMR
   WHERE ( GRAIN(:,:) )
      ! exchange of variables: [m]
      ZVOLDR_3D(:,:,J2) = ZABSCISSR(J2)**(1./XALPHAR)/ZLAMBDAR_3D(:,:)
      ZVOLDR_3D_INV(:,:,J2) = 1./ZVOLDR_3D(:,:,J2)
      ZVOLDR_3D_POW(:,:,J2) = ZVOLDR_3D(:,:,J2)**(2.+XDR)
      ! Raindrop Falling VELocity [m/s]
      ZFVELR_3D(:,:,J2) = XCR*(ZVOLDR_3D(:,:,J2)**XDR)*(XRHO00/PRHODREF(:,:))**(0.4)
      ! Reynolds number
      ZRE_3D(:,:,J2) = ZVOLDR_3D(:,:,J2)*ZFVELR_3D(:,:,J2)          &
                                            *PRHODREF(:,:)/(2.0*ZVISCA_3D(:,:))   
      ZRE_3D_SQRT(:,:,J2) = SQRT( ZRE_3D(:,:,J2) )   
      ! Critical Stokes number
      ZST_STAR_3D(:,:,J2) = (1.2+(LOG(1.+ZRE_3D(:,:,J2)))/12.)        &
                                           /(1.+LOG(1.+ZRE_3D(:,:,J2)))
   END WHERE
END DO
!
!
!------------------------------------------------------------------------------
!
!
!*       3.     AEROSOL SCAVENGING
!               ------------------
!
!
! Iteration over the aerosol type and mode
!
DO ISV = 1, NMOD_CCN+NMOD_IFN
!
   IF (ISV .LE. NMOD_CCN) THEN
      IMOD = ISV
      ISV_VAR = TNSV%NSV_LIMA_CCN_FREE - TNSV%NSV_LIMA_BEG + IMOD ! Variable number in PSVT
      INM = 1                               ! Number of species (for IFN int. mixing)
   ELSE
      IMOD = ISV - NMOD_CCN
      ISV_VAR = TNSV%NSV_LIMA_IFN_FREE - TNSV%NSV_LIMA_BEG + IMOD
      INM = NSPECIE
   END IF
!
   ZBC_SCAV_COEF(:,:,:) = 0. 
   ZMEAN_SCAV_COEF(:,:) = 0.
   ZTOT_SCAV_RATE(:,:)  = 0.
!
   GSCAV(:,:) = .FALSE.
   GSCAV(D%NIJB:D%NIJE,D%NKB:D%NKE) =GRAIN(D%NIJB:D%NIJE,D%NKB:D%NKE) .AND.        &
        (PSVT(D%NIJB:D%NIJE,D%NKB:D%NKE,ISV_VAR)>1.0E-2)
   ISCAV = COUNTJV(GSCAV(:,:),I1(:),I3(:))
!
   IF( ISCAV>=1 ) THEN
      ALLOCATE(ZVISC_RATIO(ISCAV))
      ALLOCATE(ZRHODREF(ISCAV))
      ALLOCATE(ZVISCA(ISCAV))
      ALLOCATE(ZT(ISCAV))
      ALLOCATE(ZRRT(ISCAV))
      ALLOCATE(ZCONCR(ISCAV))
      ALLOCATE(ZLAMBDAR(ISCAV))  
      ALLOCATE(ZCONCP(ISCAV))
      ALLOCATE(ZMFPA(ISCAV))
      ALLOCATE(ZTOT_SCAV_RATE1D(ISCAV))
      ALLOCATE(ZTOT_MASS_RATE1D(ISCAV))
      ALLOCATE(ZMEAN_SCAV_COEF1D(ISCAV))
      ALLOCATE(ZPABST(ISCAV))
      ALLOCATE(ZKNUDSEN(ISCAV))
      ALLOCATE(ZFACTOR(ISCAV))
!
      ALLOCATE(ZCUNSLIP(ISCAV,NDIAMP))
      ALLOCATE(ZBC_SCAV_COEF2D(ISCAV,NDIAMP))
      ALLOCATE(ZSC(ISCAV,NDIAMP))
      ALLOCATE(ZSC_INV(ISCAV,NDIAMP))
      ALLOCATE(ZSC_SQRT(ISCAV,NDIAMP))
      ALLOCATE(ZSC_3SQRT(ISCAV,NDIAMP))
      ALLOCATE(ZRELT(ISCAV,NDIAMP))
      ALLOCATE(ZDIFF(ISCAV,NDIAMP))
      ALLOCATE(ZVOLDR(ISCAV,NDIAMR))
      ALLOCATE(ZVOLDR_POW(ISCAV,NDIAMR))
      ALLOCATE(ZVOLDR_INV(ISCAV,NDIAMR))
      ALLOCATE(ZRE(ISCAV,NDIAMR))
      ALLOCATE(ZRE_INV(ISCAV,NDIAMR))
      ALLOCATE(ZRE_SQRT(ISCAV,NDIAMR))
      ALLOCATE(ZST_STAR(ISCAV,NDIAMR))
      ALLOCATE(ZFVELR(ISCAV,NDIAMR)) 
      ALLOCATE(ZST(ISCAV,NDIAMP,NDIAMR))
      ALLOCATE(ZCOL_EF(ISCAV,NDIAMP,NDIAMR))
      ALLOCATE(ZSIZE_RATIO(ISCAV,NDIAMP,NDIAMR))
!
      ZMEAN_SCAV_COEF1D(:)=0.
      ZTOT_SCAV_RATE1D(:) =0.
      ZTOT_MASS_RATE1D(:) =0.
      DO IL=1,ISCAV
         ZRHODREF(IL) =  PRHODREF(I1(IL),I3(IL))
         ZT(IL)       =     ZT_3D(I1(IL),I3(IL))
         ZRRT(IL)     =      PRRT(I1(IL),I3(IL))
         ZPABST(IL)   =    PPABST(I1(IL),I3(IL))
         ZCONCP(IL)   =      PSVT(I1(IL),I3(IL),ISV_VAR)*ZRHODREF(IL)![/m3]
         ZCONCR(IL)   = ZCONCR_3D(I1(IL),I3(IL))                       ![/m3]
         ZVISCA(IL)   = ZVISCA_3D(I1(IL),I3(IL))
         ZMFPA(IL)    =  ZMFPA_3D(I1(IL),I3(IL))
         ZVISC_RATIO(IL) = ZVISC_RATIO_3D(I1(IL),I3(IL))
         ZLAMBDAR(IL) = ZLAMBDAR_3D(I1(IL),I3(IL))
         ZFACTOR(IL)   = ZFACTOR_3D(I1(IL),I3(IL))
         ZVOLDR(IL,:) = ZVOLDR_3D(I1(IL),I3(IL),:)
         ZVOLDR_POW(IL,:) = ZVOLDR_3D_POW(I1(IL),I3(IL),:)
         ZVOLDR_INV(IL,:) = ZVOLDR_3D_INV(I1(IL),I3(IL),:)
         ZFVELR(IL,:) = ZFVELR_3D(I1(IL),I3(IL),:)
         ZRE(IL,:)    = ZRE_3D(I1(IL),I3(IL),:)
         ZRE_SQRT(IL,:) = ZRE_3D_SQRT(I1(IL),I3(IL),:)
         ZST_STAR(IL,:) = ZST_STAR_3D(I1(IL),I3(IL),:)
      ENDDO
      ZRE_INV(:,:) = 1./ZRE(:,:)

      IF (ANY(ZCONCR .EQ. 0.)) PRINT *, 'VALEUR NULLE DANS ZLAMBDAR !' 
      IF (ANY(ZLAMBDAR .EQ. 0.)) PRINT *, 'VALEUR NULLE DANS ZLAMBDAR !' 
!
!------------------------------------------------------------------------------------
!
! Loop over the different species (for IFN int. mixing)
!
      DO IJM = 1, INM  ! species (DM1,DM2,BC,O) for IFN
         IF ( ISV .LE. NMOD_CCN ) THEN          ! CCN case
            ZRHOP   = XRHO_CCN(IMOD)
            ZMDIAMP = 2*XR_MEAN_CCN(IMOD)
            ZSIGMAP = EXP(XLOGSIG_CCN(IMOD))
            ZFRACP  = 1.0
         ELSE                                   ! IFN case
            ZRHOP   = XRHO_IFN(IJM)
            ZMDIAMP = XMDIAM_IFN(IJM)
            ZSIGMAP = XSIGMA_IFN(IJM)
            ZFRACP  = XFRAC(IJM,IMOD)
         END IF
      !-----------------------------------------------------------------------------
      ! Loop over the aerosols particles diameters (log normal distribution law) :
      ! 
         DO J1=1,NDIAMP                        
            ! exchange of variables: [m]
            ZVOLDP(J1) = ZMDIAMP * EXP(ZABSCISSP(J1)*SQRT(2.)*LOG(ZSIGMAP))
            ! Cunningham slip correction factor (1+alpha*Knudsen) 
            ZKNUDSEN(:) = MIN( 20.,ZVOLDP(J1)/ZMFPA(:) )
            ZCUNSLIP(:,J1) = 1.0+2.0/ZKNUDSEN(:)*(1.257+0.4*EXP(-0.55*ZKNUDSEN(:)))
            ! Diffusion coefficient
            ZDIFF(:,J1) = CST%XBOLTZ*ZT(:)*ZCUNSLIP(:,J1)/(3.*CST%XPI*ZVISCA(:)*ZVOLDP(J1))
            ! Schmidt number
            ZSC(:,J1)       = ZVISCA(:)/(ZRHODREF(:)*ZDIFF(:,J1))  
            ZSC_INV(:,J1)   = 1./ZSC(:,J1)
            ZSC_SQRT(:,J1)  = SQRT( ZSC(:,J1) ) 
            ZSC_3SQRT(:,J1) = ZSC(:,J1)**(1./3.)  
            ! Characteristic Time Required for reaching terminal velocity 
            ZRELT(:,J1) = (ZVOLDP(J1)**2)*ZCUNSLIP(:,J1)*ZRHOP/(18.*ZVISCA(:))
            ! Density number
            ZDENS_RATIO = ZRHOP/CST%XRHOLW
            ZDENS_RATIO_SQRT = SQRT(ZDENS_RATIO)
            ! Initialisation
            ZBC_SCAV_COEF2D(:,J1)=0.
         !-------------------------------------------------------------------------
         ! Loop over the drops diameters (generalized Gamma distribution) :
         !
            DO J2=1,NDIAMR
               ! Stokes number
               ZST(:,J1,J2) = 2.*ZRELT(:,J1)*(ZFVELR(:,J2)-ZRELT(1,J1)*CST%XG)          &
                                            *ZVOLDR_INV(:,J2)
               ! Size Ratio
               ZSIZE_RATIO(:,J1,J2) = ZVOLDP(J1)*ZVOLDR_INV(:,J2)
               ! Collision Efficiency
               ZCOL_EF(:,J1,J2) = COLL_EFFI(ISCAV, NDIAMP, NDIAMR, &
                    ZRE, ZRE_INV, ZRE_SQRT, ZSC, ZSC_INV,   &
                    ZSC_SQRT, ZSC_3SQRT, ZST, ZST_STAR,          &
                    ZSIZE_RATIO, ZVISC_RATIO, ZDENS_RATIO_SQRT) 
               ! Below-Cloud Scavenging Coefficient for a fixed ZVOLDP: [/s]
               ZBC_SCAV_COEF2D(:,J1) = ZBC_SCAV_COEF2D(:,J1) +                          &
                                     ZCOL_EF(:,J1,J2) * ZWEIGHTR(J2) * ZFACTOR(:) * ZVOLDR_POW(:,J2)
            END DO 
         ! End of the loop over the drops diameters
         !--------------------------------------------------------------------------

            ! Total NUMBER Scavenging Rate of aerosol [m**-3.s**-1]
            ZTOT_SCAV_RATE1D(:) = ZTOT_SCAV_RATE1D(:) -                                 &
                                ZWEIGHTP(J1)*ZFRACP*ZCONCP(:)*ZBC_SCAV_COEF2D(:,J1)
            ! Total MASS Scavenging Rate of aerosol [kg.m**-3.s**-1]
            ZTOT_MASS_RATE1D(:) = ZTOT_MASS_RATE1D(:) +                                 &
                                ZWEIGHTP(J1)*ZFRACP*ZCONCP(:)*ZBC_SCAV_COEF2D(:,J1)   &
                                *CST%XPI/6.*ZRHOP*(ZVOLDP(J1)**3)  
         END DO
      ! End of the loop over the drops diameters
      !--------------------------------------------------------------------------

         ! Total NUMBER Scavenging Rate of aerosol [m**-3.s**-1]
         ZTOT_SCAV_RATE(:,:)=UNPACK(ZTOT_SCAV_RATE1D(:),MASK=GSCAV(:,:),FIELD=0.0)
         ! Free particles (CCN or IFN) [/s]:
         PRSVS(:,:,ISV_VAR) = MAX(PRSVS(:,:,ISV_VAR)+ZTOT_SCAV_RATE(:,:)  &
                                         * PRHODJ(:,:)/PRHODREF(:,:) , 0.0 )
         ! Total MASS Scavenging Rate of aerosol which REACH THE FLOOR because of 
         ! rain sedimentation [kg.m**-3.s**-1]
         IF (LAERO_MASS)THEN
            ZTOT_MASS_RATE(:,:) = ZTOT_MASS_RATE(:,:) +                         &
                 UNPACK(ZTOT_MASS_RATE1D(:), MASK=GSCAV(:,:), FIELD=0.0)
            CALL SCAV_MASS_SEDIMENTATION( D, HCLOUD, HDCONF, PTSTEP, KTCOUNT, PZZ, PRHODJ,     &
                                      PRHODREF, PRRT, PSVT(:,:,ISV_LIMA_SCAVMASS),&
                                      PRSVS(:,:,ISV_LIMA_SCAVMASS), PINPAP        )
            PRSVS(:,:,ISV_LIMA_SCAVMASS)=PRSVS(:,:,ISV_LIMA_SCAVMASS) +         &
                             ZTOT_MASS_RATE(:,:)*PRHODJ(:,:)/PRHODREF(:,:)
         END IF
      ENDDO
! End of the loop over the aerosol species
!--------------------------------------------------------------------------
!
!
!
      DEALLOCATE(ZFACTOR)
      DEALLOCATE(ZSC_INV)
      DEALLOCATE(ZSC_SQRT)
      DEALLOCATE(ZSC_3SQRT)
      DEALLOCATE(ZRE_INV)
      DEALLOCATE(ZRE_SQRT)
      DEALLOCATE(ZVOLDR_POW)
      DEALLOCATE(ZVOLDR_INV)
!
      DEALLOCATE(ZFVELR)
      DEALLOCATE(ZRE)
      DEALLOCATE(ZST_STAR)
      DEALLOCATE(ZST)
      DEALLOCATE(ZSIZE_RATIO)
      DEALLOCATE(ZCOL_EF)
      DEALLOCATE(ZVOLDR)
      DEALLOCATE(ZDIFF)
      DEALLOCATE(ZRELT)
      DEALLOCATE(ZSC)
      DEALLOCATE(ZCUNSLIP)
      DEALLOCATE(ZBC_SCAV_COEF2D)
!
      DEALLOCATE(ZTOT_SCAV_RATE1D)
      DEALLOCATE(ZTOT_MASS_RATE1D)
      DEALLOCATE(ZMEAN_SCAV_COEF1D)
!
      DEALLOCATE(ZRRT)
      DEALLOCATE(ZCONCR)
      DEALLOCATE(ZLAMBDAR)
      DEALLOCATE(ZCONCP)
      DEALLOCATE(ZVISC_RATIO)
      DEALLOCATE(ZRHODREF)
      DEALLOCATE(ZVISCA)
      DEALLOCATE(ZPABST)
      DEALLOCATE(ZKNUDSEN)
      DEALLOCATE(ZT)
      DEALLOCATE(ZMFPA)
   ENDIF
ENDDO
!
IF ( BUCONF%LBUDGET_SV ) then
  DO IL = 1, NMOD_CCN
    IDX = TNSV%NSV_LIMA_CCN_FREE - 1 + IL
    CALL TBUDGETS(NBUDGET_SV1 - 1 + IDX)%PTR%END_PHY(D, 'SCAV', PRSVS(:, :, IDX) )
  END DO
  DO IL = 1, NMOD_IFN
    IDX = TNSV%NSV_LIMA_IFN_FREE - 1 + IL
    CALL TBUDGETS(NBUDGET_SV1 - 1 + IDX)%PTR%END_PHY(D, 'SCAV', PRSVS(:, :, IDX) )
  END DO
  IF ( LAERO_MASS ) then
     CALL TBUDGETS(NBUDGET_SV1 - 1 + TNSV%NSV_LIMA_SCAVMASS)%PTR%END_PHY(D, &
          'SCAV', PRSVS(:, :, TNSV%NSV_LIMA_SCAVMASS) )
  END IF
END IF
!------------------------------------------------------------------------------
!
!
!*       3.     SUBROUTINE AND FUNCTION
!               -----------------------
!
!
IF (LHOOK) CALL DR_HOOK('LIMA_PRECIP_SCAVENGING', 1, ZHOOK_HANDLE)
CONTAINS
!
!------------------------------------------------------------------------------
!     ##########################################################################
      SUBROUTINE SCAV_MASS_SEDIMENTATION( D, HCLOUD, HDCONF, PTSTEP, KTCOUNT, PZZ, PRHODJ,&
                                PRHODREF, PRAIN, PSVT_MASS, PRSVS_MASS, PINPAP )
!     ##########################################################################
!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the total mass of aerosol 
!!    scavenged by precipitations
!!
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!      None
!!     
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_PARAMETERS
!!          JPHEXT       : Horizontal external points number
!!          JPVEXT       : Vertical external points number
!!
!!    REFERENCE
!!    ---------
!!      Book1 of the documentation ( routine CH_AQUEOUS_SEDIMENTATION )
!!
!!    AUTHOR
!!    ------
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    22/07/07
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_DIMPHYEX,        ONLY: DIMPHYEX_T
USE MODD_PARAMETERS
!
USE MODD_PARAM_LIMA,      ONLY : XCEXVT, XRTMIN
USE MODD_PARAM_LIMA_WARM, ONLY : XBR, XDR, XFSEDRR
USE YOMHOOK, ONLY:LHOOK, DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
TYPE(DIMPHYEX_T),         INTENT(IN)    :: D
CHARACTER (LEN=4),        INTENT(IN)    :: HCLOUD  ! Cloud parameterization
CHARACTER(LEN=5),         INTENT(IN)    :: HDCONF
REAL,                     INTENT(IN)    :: PTSTEP  ! Time step  
INTEGER,                  INTENT(IN)    :: KTCOUNT ! Current time step number
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PZZ     ! Height (z)
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PRHODJ  ! Dry Density [kg]
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PRAIN   ! Rain water m.r. source
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PSVT_MASS  ! Precip. aerosols at t
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PRSVS_MASS ! Precip. aerosols source
!
REAL, DIMENSION(D%NIJT),     INTENT(INOUT) :: PINPAP
!
!*       0.2   Declarations of local variables :
!
INTEGER :: IK, IN                         ! Loop indexes 
INTEGER :: IIB, IIE, IJB, IJE, IKB, IKE   ! Physical domain
!
REAL    :: ZTSPLITR      ! Small time step for rain sedimentation
REAL    :: ZTSTEP        ! Large time step for rain sedimentation
!
!
LOGICAL, DIMENSION(D%NIJT,D%NKT) &
                                :: GSEDIM   ! where to compute the SED processes
INTEGER :: ISEDIM 
INTEGER , DIMENSION(SIZE(GSEDIM)) :: I1,I2,I3 ! Used to replace the COUNT
INTEGER                           :: IL       ! and PACK intrinsics
!
!
REAL,    DIMENSION(D%NIJT,D%NKT)   &
                                :: ZW,    & ! work array
                                   ZWSED, & ! sedimentation fluxes
                                   ZZS      ! Rain water m.r. source
!
REAL, DIMENSION(:), ALLOCATABLE :: ZRRS,     & ! Rain water m.r. source
                                   ZRHODREF, & ! RHO Dry REFerence
                                   ZZW         ! Work array
!
REAL                            :: ZRTMIN3
!
!
REAL                            :: ZVTRMAX, ZDZMIN, ZT
REAL,    SAVE                   :: ZEXSEDR
LOGICAL, SAVE                   :: GSFIRSTCALL = .TRUE.
INTEGER, SAVE                   :: ISPLITR
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!*       1.     COMPUTE THE LOOP BOUNDS
!               -----------------------
!
IF (LHOOK) CALL DR_HOOK('SCAV_MASS_SEDIMENTATION', 0, ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
!*       2.     COMPUTE THE SEDIMENTATION (RS) SOURCE
!            -------------------------------------
!
!*       2.1    splitting factor for high Courant number C=v_fall*(del_Z/del_T)
!  
FIRSTCALL : IF (GSFIRSTCALL) THEN
   GSFIRSTCALL = .FALSE.
   ZVTRMAX = 10.                          
   ZDZMIN = MINVAL(PZZ(D%NIJB:D%NIJE,D%NKB+1:D%NKE+1)-PZZ(D%NIJB:D%NIJE,D%NKB:D%NKE))
   ISPLITR = 1
   SPLIT : DO
      ZT = 2.* PTSTEP / REAL(ISPLITR)
      IF ( ZT * ZVTRMAX / ZDZMIN .LT. 1.) EXIT SPLIT
      ISPLITR = ISPLITR + 1
   END DO SPLIT
!
   ZEXSEDR = (XBR+XDR+1.0)/(XBR+1.0) 
!
END IF FIRSTCALL
!
!*       2.2    time splitting loop initialization        
!
IF( (KTCOUNT==1) .AND. (HDCONF=='START') ) THEN
  ZTSPLITR = PTSTEP / REAL(ISPLITR)       ! Small time step
  ZTSTEP   = PTSTEP                        ! Large time step
  ELSE
  ZTSPLITR= 2. * PTSTEP / REAL(ISPLITR)
  ZTSTEP  = 2. * PTSTEP
END IF
!
!*       2.3    compute the fluxes
! 
!  optimization by looking for locations where
!  the precipitating fields are larger than a minimal value only !!!
!
ZRTMIN3 = XRTMIN(3) / ZTSTEP
ZZS(:,:) = PRAIN(:,:)
DO IN = 1 , ISPLITR
   GSEDIM(:,:) = .FALSE.
   GSEDIM(D%NIJB:D%NIJE,D%NKB:D%NKE) = ZZS(D%NIJB:D%NIJE,D%NKB:D%NKE) > ZRTMIN3
! 
   ISEDIM = COUNTJV( GSEDIM(:,:),I1(:),I3(:))
   IF( ISEDIM >= 1 ) THEN
      IF( IN==1 ) THEN
         ZZS(:,:) = ZZS(:,:) * ZTSTEP
         DO IK = IKB , IKE-1
            ZW(:,IK) =ZTSPLITR*2./(PRHODREF(:,IK)*(PZZ(:,IK+2)-PZZ(:,IK)))
         END DO
         ZW(:,IKE)  =ZTSPLITR/(PRHODREF(:,IKE)*(PZZ(:,IKE+1)-PZZ(:,IKE)))
      END IF
      ALLOCATE(ZRRS(ISEDIM)) 
      ALLOCATE(ZRHODREF(ISEDIM))
      DO IL=1,ISEDIM
         ZRRS(IL) = ZZS(I1(IL),I3(IL))
         ZRHODREF(IL) =  PRHODREF(I1(IL),I3(IL))
      ENDDO
      ALLOCATE(ZZW(ISEDIM)) ; ZZW(:) = 0.0
!
!*       2.2.1   for rain
!
      ZZW(:) = XFSEDRR * ZRRS(:)**(ZEXSEDR) * ZRHODREF(:)**(ZEXSEDR-XCEXVT)
      ZWSED(:,:) = UNPACK( ZZW(:),MASK=GSEDIM(:,:),FIELD=0.0 )
      DO IK = IKB , IKE
         ZZS(:,IK) = ZZS(:,IK) + ZW(:,IK)*(ZWSED(:,IK+1)-ZWSED(:,IK))
      END DO
      IF( IN==1 ) THEN
         PINPAP(:) = ZWSED(:,IKB)*                                            &
                             ( PSVT_MASS(:,IKB)/MAX(ZRTMIN3,PRRT(:,IKB)) )
      END IF
      DEALLOCATE(ZRHODREF)
      DEALLOCATE(ZRRS)
      DEALLOCATE(ZZW)
      IF( IN==ISPLITR ) THEN
         GSEDIM(:,:) = .FALSE.
         GSEDIM(D%NIJB:D%NIJE,D%NKB:D%NKE) = ZZS(D%NIJB:D%NIJE,D%NKB:D%NKE) > ZRTMIN3
         ZWSED(:,:) = 0.0
         WHERE( GSEDIM(:,:) ) 
            ZWSED(:,:) = 1.0/ZTSTEP - PRAIN(:,:)/ZZS(:,:)
         END WHERE
      END IF
   END IF
END DO
!
! Apply the rain sedimentation rate to the WR_xxx aqueous species
!
PRSVS_MASS(:,:) = PRSVS_MASS(:,:) + ZWSED(:,:)*PSVT_MASS(:,:)
!
IF (LHOOK) CALL DR_HOOK('SCAV_MASS_SEDIMENTATION', 1, ZHOOK_HANDLE)
END SUBROUTINE SCAV_MASS_SEDIMENTATION
!
!------------------------------------------------------------------------------
!
!###################################################################
FUNCTION COLL_EFFI (KSCAV,KDIAMP,KDIAMR, &
                    PRE, PRE_INV, PRE_SQRT, PSC, PSC_INV, PSC_SQRT, &
                    PSC_3SQRT, PST, PST_STAR, PSIZE_RATIO,          &
                    PVISC_RATIO, PDENS_RATIO_SQRT) RESULT(PCOL_EF)
!###################################################################
!
!Compute the Raindrop-Aerosol Collision Efficiency
!
!*             0. DECLARATIONS  
!              ---------------
!
  IMPLICIT NONE
!
  INTEGER :: I
  !
  INTEGER,              INTENT(IN)    :: KSCAV
  INTEGER,              INTENT(IN)    :: KDIAMP
  INTEGER,              INTENT(IN)    :: KDIAMR
  REAL, DIMENSION(KSCAV,KDIAMR), INTENT(IN)    :: PRE         
  REAL, DIMENSION(KSCAV,KDIAMR), INTENT(IN)    :: PRE_INV         
  REAL, DIMENSION(KSCAV,KDIAMR), INTENT(IN)    :: PRE_SQRT         
  REAL, DIMENSION(KSCAV,KDIAMP), INTENT(IN)    :: PSC     
  REAL, DIMENSION(KSCAV,KDIAMP), INTENT(IN)    :: PSC_INV    
  REAL, DIMENSION(KSCAV,KDIAMP), INTENT(IN)    :: PSC_SQRT     
  REAL, DIMENSION(KSCAV,KDIAMP), INTENT(IN)    :: PSC_3SQRT     
  REAL, DIMENSION(KSCAV,KDIAMR), INTENT(IN)    :: PST_STAR   
!
  REAL, DIMENSION(KSCAV,KDIAMP,KDIAMR), INTENT(IN)  :: PST         
  REAL, DIMENSION(KSCAV,KDIAMP,KDIAMR), INTENT(IN)  :: PSIZE_RATIO 
!
  REAL, DIMENSION(KSCAV), INTENT(IN)    :: PVISC_RATIO 
  REAL,               INTENT(IN)    :: PDENS_RATIO_SQRT 
!
  REAL, DIMENSION(KSCAV)      :: PCOL_EF   !result : collision efficiency
!
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
  IF (LHOOK) CALL DR_HOOK('COLL_EFFI', 0, ZHOOK_HANDLE)
  PCOL_EF(:) = (4.*PSC_INV(:,J1)*PRE_INV(:,J2)*(1.+0.4*PRE_SQRT(:,J2)              &
                  *PSC_3SQRT(:,J1)+0.16*PRE_SQRT(:,J2)*PSC_SQRT(:,J1)))      &
              +(4.*PSIZE_RATIO(:,J1,J2)*(PVISC_RATIO(:)                    & 
              +(1.+2.*PRE_SQRT(:,J2))*PSIZE_RATIO(:,J1,J2)))  
  DO I=1,KSCAV
    IF (PST(I,J1,J2)>PST_STAR(I,J2)) THEN            
            PCOL_EF(I) = PCOL_EF(I)                                           &
                    +(PDENS_RATIO_SQRT*((PST(I,J1,J2)-PST_STAR(I,J2))        &
                    /(PST(I,J1,J2)-PST_STAR(I,J2)+2./3.))**(3./2.))
    ENDIF
  ENDDO
  IF (LHOOK) CALL DR_HOOK('COLL_EFFI', 1, ZHOOK_HANDLE)
  END FUNCTION COLL_EFFI
!
!------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_PRECIP_SCAVENGING
