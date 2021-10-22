SUBROUTINE ARO_MNHDUST(KKL,KLON,KLEV,KSV,PTSTEP &
       ,PSVTIN      & !I [moments/molec_{air}] Transported moments of dust
       ,PZZ         & !I [m] height of layers
       ,PDZZ        & !I [m] height of layers
       ,PPABST      & !I Pressure
       ,PTHT        & !I Potential temperature
       ,PRHODREF    & !I [kg/m3] density of air
       ,KSWB        & !I [nbr] number of shortwave bands
       ,KTCOUNT     & !I number of time step
       ,PSVT        & !O [moments/molec_{air}] Transported moments of dust
       ,PPIZA_WVL   & !IO [-] single scattering albedo of dust layer for all SW wavelength
       ,PCGA_WVL    & !IO [-] assymetry factor for dust layer for all SW wavelength
       ,PTAUREL_WVL & !IO [-] opt.depth/opt.depth(550) for dust layer for all SWwvl
       ,PAER        & !IO [-]  ext coeff at 550 for dust layer
       ,NDIAG       & !I [-] nb of diagnostics
       ,PPEZDIAG    & !IO [-] diag Nb/m3,ug/m3,rg(nb;um),rg(m;um),SSA,assym,AOD/550,mode & wvl
        )
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK


!*** *ARO_MNHDUST*
!
!    PURPOSE
!    -------

!     Interface routine for initialisation of dust optical properties
!     before radiation scheme call

!     AUTHOR.
!     -------
!      Y. Seity (CNRM/GMAP)
!      10-10-05

!    MODIFICATIONS
!    -------------
!    P. Tulet 10/02/06
!    P. Tulet add scavenging 10/02/08
!    M. Mokhtari & A. Ambar 09/2016: inversion levels, call inv_levels.F90
!
!    EXTERNAL
!     -------



USE MODE_DUSTOPT, ONLY: DUSTOPT_GET
USE MODE_DUST_PSD      ! include PPP2DUST
USE MODI_INIT_DUST
USE MODI_SEDIM_DUST
USE MODD_DUST          ! LDUST=.FALSE. ; NMODE_DST=3
USE MODD_NSV, ONLY : NSV_DSTBEG, NSV_DSTEND
USE MODI_INV_LEVELS
 IMPLICIT NONE

!INPUT
INTEGER, INTENT(IN)   :: KKL      ! vertical levels ordering 1: MNH -1: ARPEGE
INTEGER, INTENT(IN)   :: KLON     ! NPROMA under CPG
INTEGER, INTENT(IN)   :: KLEV     ! Number of vertical levels
INTEGER, INTENT(IN)   :: KSV      ! Number of passive scalar
INTEGER, INTENT(IN)   :: KSWB     ! Number of shortwave wavelengths
INTEGER, INTENT(IN)   :: KTCOUNT  ! Number of time step
REAL,    INTENT(IN)   :: PTSTEP   ! Time step in s
REAL, DIMENSION(KLON,1,KLEV,KSV),INTENT(IN) :: PSVTIN !I [moments/molec_{air}] transported moments of dust
INTEGER, INTENT(IN)   :: NDIAG   ! nb of diagnostics
REAL, DIMENSION(KLON,1,KLEV),INTENT(IN)        :: PZZ   !I [m] height of layers
REAL, DIMENSION(KLON,1,KLEV),INTENT(IN)        :: PDZZ  !I [m] layers thikness
REAL, DIMENSION(KLON,1,KLEV),INTENT(IN)        :: PRHODREF   !I [kg/m3] density of air
REAL, DIMENSION(KLON,1,KLEV),INTENT(IN)        :: PTHT       !I [K] potentiel temperature
REAL, DIMENSION(KLON,1,KLEV),INTENT(IN)        :: PPABST     !I [Pa] pressure

!OUTPUT
REAL, DIMENSION(KLON,1,KLEV,KSV),INTENT(OUT) :: PSVT !O [moments/molec_{air}] transported moments of dust
REAL, DIMENSION(KLON,KLEV,KSWB),INTENT(INOUT)    :: PPIZA_WVL !O [-] SSA of dust layer for all SW wvl
REAL, DIMENSION(KLON,KLEV,KSWB),INTENT(INOUT)    :: PCGA_WVL  !O [-] assym factor for dust layer for all SW
REAL, DIMENSION(KLON,KLEV,KSWB),INTENT(INOUT)    :: PTAUREL_WVL !O [-] AOD/AOD(550) for dust layer for all SW
REAL, DIMENSION(KLON,KLEV,NDIAG),INTENT(INOUT)   :: PPEZDIAG !O [-] Diagnostics table
REAL, DIMENSION(KLON,KLEV),INTENT(INOUT)         :: PAER  !O [-] AOD/AOD(550) for dust layer for all SW wvl
!
!*       0.2   Declarations of local variables :
INTEGER :: JLEV,JK,JK1,IMOD,INSWB,WVL_IDX
INTEGER :: IKE, JSV, JL
REAL    :: ZSIGMIN
REAL, DIMENSION(KLON,1,KLEV+2,KSWB) :: ZPIZA_DST_TMP, ZCGA_DST_TMP, ZTAUREL_DST_TMP, &
         &                               ZPIZAZ,        ZCGAZ,        ZTAUAZ
REAL, DIMENSION(KLON,1,KLEV+2,NMODE_DST) :: ZSIGDST,       ZRGDST,       ZNDST,           &
         &                                              ZRGMDST,      ZMDST
REAL, DIMENSION(KLON,1,KLEV+2) :: ZAER

REAL, DIMENSION(KLON,KLEV,KSWB)       :: ZPIZAZ_WVL
REAL, DIMENSION(KLON,KLEV,KSWB)       :: ZCGAZ_WVL
REAL, DIMENSION(KLON,KLEV,KSWB)       :: ZTAUAZ_WVL
REAL, DIMENSION(KLON,KLEV,NMODE_DST)  :: ZRGDST_MOD
REAL, DIMENSION(KLON,KLEV,NMODE_DST)  :: ZRGMDST_MOD
REAL, DIMENSION(KLON,KLEV,NMODE_DST)  :: ZNDST_MOD
REAL, DIMENSION(KLON,KLEV,NMODE_DST)  :: ZMDST_MOD
!
REAL,DIMENSION(KLON,1,KLEV+2)         :: ZZZ        ! Local value of PZZ
REAL,DIMENSION(KLON,1,KLEV+2)         :: ZDZZ       ! Local value of PDZZ
REAL,DIMENSION(KLON,1,KLEV+2)         :: ZPABST     ! Local value of ZPABST
REAL,DIMENSION(KLON,1,KLEV+2)         :: ZTHT       ! Local value of ZTHT
REAL,DIMENSION(KLON,1,KLEV+2)         :: ZRHODREF   ! Local value of ZRHODREF
REAL,DIMENSION(KLON,1,KLEV+2,KSV)         :: ZSVT       ! Local value of ZSVT
!
REAL,DIMENSION(KLON,1,KLEV)         :: ZZPABST     ! Local value of ZPABST
REAL,DIMENSION(KLON,1,KLEV)         :: ZZTHT       ! Local value of ZTHT
REAL,DIMENSION(KLON,1,KLEV)         :: ZZRHODREF   ! Local value of ZRHODREF
REAL,DIMENSION(KLON,1,KLEV,KSV)     :: ZZSVT       ! Local value of ZSVT
!--------------------------------------------------------------------


! 1) Initialization of dust

REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('ARO_MNHDUST',0,ZHOOK_HANDLE)
!
!*       1.     PRELIMINARY COMPUTATIONS
  !initialisation de ZZZ
DO JL = 1,KLON
   DO JK = 2 , KLEV+1
      ZZZ(JL,1,JK)=PZZ(JL,1,KLEV+2-JK)
   ENDDO
   ZZZ(JL,1,1) = 2*ZZZ(JL,1,2)-ZZZ(JL,1,3)
   ZZZ(JL,1,KLEV+2) = 2*ZZZ(JL,1,KLEV+1)-ZZZ(JL,1,KLEV)
ENDDO
  !initialisation de ZDZZ
DO JL = 1,KLON
   DO JK = 1 , KLEV+1
      ZDZZ(JL,1,JK)=ZZZ(JL,1,JK+1)-ZZZ(JL,1,JK)
   ENDDO
   ZDZZ(JL,1,KLEV+2)=ZZZ(JL,1,KLEV+2)-ZZZ(JL,1,KLEV+1)
ENDDO
!
ZZRHODREF(:,:,:) = PRHODREF(:,:,:)
ZZPABST(:,:,:)   = PPABST(:,:,:)
ZZTHT(:,:,:)     = PTHT(:,:,:)
ZZSVT(:,:,:,:)   = PSVTIN(:,:,:,:)
!
CALL INV_LEVELS(1,KLON,KLEV,1,ZZRHODREF,ZRHODREF)
CALL INV_LEVELS(1,KLON,KLEV,1,ZZPABST,ZPABST)
CALL INV_LEVELS(1,KLON,KLEV,1,ZZTHT,ZTHT)
!initialisation de ZSVT
DO JSV=1,KSV
CALL INV_LEVELS(1,KLON,KLEV,1,ZZSVT(:,:,:,JSV),ZSVT(:,:,:,JSV))
ENDDO
!
IF ((KTCOUNT == 1).AND.(LDSTINIT)) THEN
  CALL INIT_DUST(ZSVT(:,:,:,NSV_DSTBEG:NSV_DSTEND), ZRHODREF)
END IF

!Get dust optical properties from look up tables

   CALL DUSTOPT_GET(                             &
        ZSVT(:,:,:,NSV_DSTBEG:NSV_DSTEND)        &  !I [ppp] Dust scalar concentration
        ,ZZZ(:,:,:)                              &  !I [m] height of layers
        ,ZRHODREF(:,:,:)                         &  !I [kg/m3] density of air
        ,ZPIZA_DST_TMP                           &  !O [-] single scattering albedo of dust
        ,ZCGA_DST_TMP                            &  !O [-] assymetry factor for dust
        ,ZTAUREL_DST_TMP                         &  !O [-] opt.depth(wvl=lambda)/opt.depth(wvl=550nm)
        ,ZAER(:,:,:)                             &  !O [-] optical depth of dust at wvl=550nm
        ,KSWB                                    &  !I |nbr] number of shortwave bands
        )

! Compute SSA, AOD, assymetry factor for clear sky (dust aerosols)

ZTAUAZ(:,:,:,:)=0.
ZPIZAZ(:,:,:,:)=0.
ZCGAZ(:,:,:,:)=0.
DO WVL_IDX=1,KSWB
! Ponderation of aerosol optical depht in case of explicit factor
! ti
            ZTAUAZ(:,:,:,WVL_IDX) = ZTAUAZ(:,:,:,WVL_IDX) + &
                        ZAER(:,:,:)                       * &
                        ZTAUREL_DST_TMP(:,:,:,WVL_IDX)
! wi*ti
            ZPIZAZ(:,:,:,WVL_IDX) = ZPIZAZ(:,:,:,WVL_IDX) + &
                        ZAER(:,:,:)                       * &
                        ZTAUREL_DST_TMP(:,:,:,WVL_IDX)    * &
                        ZPIZA_DST_TMP(:,:,:,WVL_IDX)
! gi*wi*ti
        ZCGAZ(:,:,:,WVL_IDX) = ZCGAZ(:,:,:,WVL_IDX)              + &
                               ZAER(:,:,:)                       * &
                               ZTAUREL_DST_TMP(:,:,:,WVL_IDX)    * &
                               ZPIZA_DST_TMP(:,:,:,WVL_IDX)      * &
                               ZCGA_DST_TMP(:,:,:,WVL_IDX)
ENDDO
! Ponderation of assymetry factor
ZCGAZ(:,:,:,:) = ZCGAZ(:,:,:,:) / ZPIZAZ(:,:,:,:)
! Ponderation of SSA
ZPIZAZ(:,:,:,:) = ZPIZAZ(:,:,:,:) / ZTAUAZ(:,:,:,:)

! Compute and store Standard deviation, median radius, concentration of dust from moments

    CALL PPP2DUST (ZSVT(:,:,:,NSV_DSTBEG:NSV_DSTEND), &
                   ZRHODREF(:,:,:),                   &
                   ZSIGDST(:,:,:,:),                  &
                   ZRGDST(:,:,:,:),                   &
                   ZNDST(:,:,:,:))
ZRGMDST(:,:,:,:)=0.
ZMDST(:,:,:,:)=0.
! Compute mass median radius in um in all dust mode
ZRGMDST(:,:,:,:) = ZRGDST(:,:,:,:) / (exp(-3.*(log(ZSIGDST(:,:,:,:)))**2))
! Compute integrated mass concentration ug/m3 in all dust mode
ZMDST(:,:,:,:) =  ZNDST(:,:,:,:)*4./3.*3.14*2500.*1e9  *&   !kg-->ug
                  (ZRGDST(:,:,:,:)**3)*1.d-18          *&   !um-->m
                  exp(4.5*log(ZSIGDST(:,:,:,:))*log(ZSIGDST(:,:,:,:)))

!Transform from vector of type #lon #lat #lev #wvl or #mod
!to vectors of type #points, #inverted levs, #wavelengths or #modes

JLEV             = KLEV
ZTAUAZ_WVL(:,:,:)= 0.
ZPIZAZ_WVL(:,:,:)= 0.
ZCGAZ_WVL(:,:,:) = 0.
ZRGDST_MOD(:,:,:)= 0.
ZRGMDST_MOD(:,:,:)=0.
ZNDST_MOD(:,:,:) = 0.
ZMDST_MOD(:,:,:) = 0.


DO JK=1,JLEV
  JK1=JLEV+2-JK
  PPIZA_WVL(:,JK,:)   = ZPIZA_DST_TMP(:,1,JK1,:)
  PCGA_WVL(:,JK,:)    = ZCGA_DST_TMP(:,1,JK1,:)
  PTAUREL_WVL(:,JK,:) = ZTAUREL_DST_TMP(:,1,JK1,:)
  PAER(:,JK)          = ZAER(:,1,JK1)
  ZPIZAZ_WVL(:,JK,:)  = ZPIZAZ(:,1,JK1,:)
  ZCGAZ_WVL(:,JK,:)   = ZCGAZ(:,1,JK1,:)
  ZTAUAZ_WVL(:,JK,:)  = ZTAUAZ(:,1,JK1,:)
  ZRGDST_MOD(:,JK,:)  = ZRGDST(:,1,JK1,:)
  ZRGMDST_MOD(:,JK,:) = ZRGMDST(:,1,JK1,:)
  ZNDST_MOD(:,JK,:)   = ZNDST(:,1,JK1,:)
  ZMDST_MOD(:,JK,:)   = ZMDST(:,1,JK1,:)
ENDDO

! somme on every mode (different dust diameter)

IF (SIZE(PPEZDIAG, 3) .GE. 4) THEN
PPEZDIAG(:,:,1:4) = 0.
DO IMOD = 1,NMODE_DST       ! NMODE_DST!=3  ! PPEZDIAG<25
  PPEZDIAG(:,:,1) = PPEZDIAG(:,:,1) + ZNDST_MOD(:,:,IMOD)    ! Nb/m3
  PPEZDIAG(:,:,2) = PPEZDIAG(:,:,2) + ZMDST_MOD(:,:,IMOD)    ! Mass(ug)/m3
  PPEZDIAG(:,:,3) = PPEZDIAG(:,:,3) + ZRGDST_MOD(:,:,IMOD)   ! RG nb (um)
  PPEZDIAG(:,:,4) = PPEZDIAG(:,:,4) + ZRGMDST_MOD(:,:,IMOD)  ! RG m (um)
ENDDO
END IF
IF (SIZE(PPEZDIAG, 3) .GE. 7) THEN
! somme on every level (integreted value)
PPEZDIAG(:,:,5:7) = 0.
DO JK=1,JLEV
! to have lwv=550nm : NSWB = 3     !NSWB!=6  ! PPEZDIAG<25
  PPEZDIAG(:,1,5) = PPEZDIAG(:,1,5) + ZPIZAZ_WVL(:,JK,3) ! SSA(550)
  PPEZDIAG(:,1,6) = PPEZDIAG(:,1,6) + ZCGAZ_WVL(:,JK,3)  ! assymetry(550)
  PPEZDIAG(:,1,7) = PPEZDIAG(:,1,7) + ZTAUAZ_WVL(:,JK,3) ! AOD/AOD(550)
ENDDO
DO JK=1,JLEV
  PPEZDIAG(:,JK,5) = PPEZDIAG(:,1,5)
  PPEZDIAG(:,JK,6) = PPEZDIAG(:,1,6)
  PPEZDIAG(:,JK,7) = PPEZDIAG(:,1,7)
ENDDO
END IF

IF (LSEDIMDUST) THEN
   IKE = KLEV+2
   CALL SEDIM_DUST(ZTHT(:,:,1:IKE), PTSTEP,  &
                  ZRHODREF(:,:,1:IKE),       &
                  ZPABST(:,:,1:IKE),         &
                  ZDZZ(:,:,1:IKE),          &
                  ZSVT(:,:,1:IKE,NSV_DSTBEG:NSV_DSTEND))
ENDIF
!return to aladin/arome ZSVT
DO JSV=1,KSV
CALL INV_LEVELS(1,KLON,KLEV,-1,PSVT(:,:,:,JSV),ZSVT(:,:,:,JSV))
ENDDO
!
IF (LHOOK) CALL DR_HOOK('ARO_MNHDUST',1,ZHOOK_HANDLE)
END SUBROUTINE ARO_MNHDUST

