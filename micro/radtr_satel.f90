!MNH_LIC Copyright 2000-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!    #######################
     MODULE MODI_RADTR_SATEL 
!    #######################
INTERFACE
!
     SUBROUTINE RADTR_SATEL(KYEARF, KMONTHF, KDAYF, PSECF,         &
                KDLON, KFLEV, KSTATM, KRAD_COLNBR, PEMIS, PCCO2,   &
                PTSRAD, PSTATM, PTHT, PRT, PPABST, PZZ,            &
                PSIGS, PMFCONV, PCLDFR, OUSERI, OSIGMAS,           &
                OSUBG_COND, ORAD_SUBG_COND, PIRBT, PWVBT, KGEO,PSIGQSAT )
!
INTEGER, INTENT(IN) :: KYEARF  ! year of Final date
INTEGER, INTENT(IN) :: KMONTHF ! month of Final date
INTEGER, INTENT(IN) :: KDAYF   ! day of Final date
REAL,    INTENT(IN) :: PSECF   ! number of seconds since date at 00 UTC 
!
INTEGER, INTENT(IN)   :: KDLON !number of columns where the
                               !radiation calculations are performed
INTEGER, INTENT(IN)   :: KFLEV !number of vertical levels where the
                               !radiation calculations are performed
INTEGER, INTENT(IN)   :: KSTATM  !index of the standard atmosphere level
                                 !just above the model top
INTEGER, INTENT(IN)   :: KRAD_COLNBR !factor by which the memory is split
!
REAL, DIMENSION(:,:),     INTENT(IN) :: PEMIS  !Surface IR EMISsivity
REAL,                     INTENT(IN) :: PCCO2  !CO2 content
REAL, DIMENSION(:,:),     INTENT(IN) :: PTSRAD !RADiative Surface Temperature
REAL, DIMENSION(:,:),     INTENT(IN) :: PSTATM !selected standard atmosphere
!
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PTHT   !THeta at t
REAL, DIMENSION(:,:,:,:), INTENT(IN) :: PRT    !moist variables at t
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PPABST !pressure at t
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PZZ    !Model level heights
!
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PSIGS  ! Sigma_s from turbulence scheme
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PMFCONV! convective mass flux (kg /s m^2)
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PCLDFR  ! cloud fraction
!
LOGICAL, INTENT(IN)                  :: OUSERI ! logical switch to compute both
                                               ! liquid and solid condensate (OUSERI=.TRUE.)
                                               ! or only liquid condensate (OUSERI=.FALSE.)
LOGICAL, INTENT(IN)                  :: OSIGMAS! use present global Sigma_s values
                                               ! or that from turbulence scheme
LOGICAL, INTENT(IN)                  :: OSUBG_COND ! Switch for Subgrid Condensation
                                                   ! (prognotic mode)
LOGICAL, INTENT(IN)                  :: ORAD_SUBG_COND ! Switch for Subgrid Condensation
                                                       ! (diagnostic mode)
!
REAL, DIMENSION(:,:),     INTENT(OUT):: PIRBT  !IR Brightness Temp. (K)
REAL, DIMENSION(:,:),     INTENT(OUT):: PWVBT  !WV Brightness Temp. (K)
!
INTEGER, INTENT(IN)   :: KGEO   !SATELLITE INDEX
REAL, INTENT(IN)                            :: PSIGQSAT ! use an extra "qsat" variance contribution (OSIGMAS case)
!
END SUBROUTINE RADTR_SATEL
END INTERFACE
END MODULE MODI_RADTR_SATEL
!    #####################################################################
     SUBROUTINE RADTR_SATEL(KYEARF, KMONTHF, KDAYF, PSECF,         &
                KDLON, KFLEV, KSTATM, KRAD_COLNBR, PEMIS, PCCO2,   &
                PTSRAD, PSTATM, PTHT, PRT, PPABST, PZZ,            &
                PSIGS, PMFCONV, PCLDFR, OUSERI, OSIGMAS,           &
                OSUBG_COND, ORAD_SUBG_COND, PIRBT, PWVBT, KGEO,PSIGQSAT)
!    #####################################################################
!
!!****  *RADTR_SATEL* - 
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!    Chaboureau, J.-P., J.-P. Cammas, P. Mascart, J.-P. Pinty, C. Claud, R. Roca,
!!       and J.-J. Morcrette, 2000: Evaluation of a cloud system life-cycle simulated
!!       by Meso-NH during FASTEX using METEOSAT radiances and TOVS-3I cloud retrievals.
!!       Q. J. R. Meteorol. Soc., 126, 1735-1750.
!!    Chaboureau, J.-P. and P. Bechtold, 2002: A simple cloud parameterization from
!!       cloud resolving model data: Theory and application. J. Atmos. Sci., 59, 2362-2372.
!!
!!    AUTHOR
!!    ------
!!      J.-P. Chaboureau       *L.A.*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    29/03/00
!!      J.-P. Chaboureau 15/04/03  add call to the subgrid condensation scheme
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!!      G.Delautier 04/2016 : BUG JPHEXT
!!      S. Riette 11/2016 : Condensation interface changed
!  P. Wautelet 26/04/2019: replace non-standard FLOAT function by REAL function
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
USE MODD_PARAMETERS
USE MODD_GRID_n
USE MODD_RAIN_ICE_PARAM,   ONLY: RAIN_ICE_PARAM
USE MODD_NEB,              ONLY: NEB
USE MODD_TURB_n,           ONLY: TURBN
USE MODD_DIMPHYEX,         ONLY: DIMPHYEX_t
!
USE MODD_RAD_TRANSF
USE MODE_ll
USE MODE_FILL_DIMPHYEX,   ONLY: FILL_DIMPHYEX
!
USE MODI_INIT_NBMOD
USE MODI_DETER_ANGLE
USE MODI_MAKE_RADSAT
!
USE MODI_CONDENSATION
!
IMPLICIT NONE
!
!*       0.1   DECLARATIONS OF DUMMY ARGUMENTS :
!
INTEGER, INTENT(IN) :: KYEARF  ! year of Final date
INTEGER, INTENT(IN) :: KMONTHF ! month of Final date
INTEGER, INTENT(IN) :: KDAYF   ! day of Final date
REAL,    INTENT(IN) :: PSECF   ! number of seconds since date at 00 UTC 
!
INTEGER, INTENT(IN)   :: KDLON   !number of columns where the
                                 ! radiation calculations are performed
INTEGER, INTENT(IN)   :: KFLEV   !number of vertical levels where the
                                 ! radiation calculations are performed
INTEGER, INTENT(IN)   :: KSTATM  !index of the standard atmosphere level
                                 !just above the model top
INTEGER, INTENT(IN)   :: KRAD_COLNBR !factor by which the memory is split
!
REAL, DIMENSION(:,:),     INTENT(IN) :: PEMIS  !Surface IR EMISsivity
REAL,                     INTENT(IN) :: PCCO2  !CO2 content
REAL, DIMENSION(:,:),     INTENT(IN) :: PTSRAD !RADiative Surface Temperature
REAL, DIMENSION(:,:),     INTENT(IN) :: PSTATM !selected standard atmosphere
!
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PTHT   !THeta at t
REAL, DIMENSION(:,:,:,:), INTENT(IN) :: PRT    !moist variables at t
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PPABST !pressure at t
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PZZ    !Model level heights
!
!
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PSIGS  ! Sigma_s from turbulence scheme
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PMFCONV! convective mass flux (kg /s m^2)
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PCLDFR  ! cloud fraction
!
LOGICAL, INTENT(IN)                  :: OUSERI ! logical switch to compute both
                                               ! liquid and solid condensate (OUSERI=.TRUE.)
                                               ! or only liquid condensate (OUSERI=.FALSE.)
LOGICAL, INTENT(IN)                  :: OSIGMAS! use present global Sigma_s values
                                               ! or that from turbulence scheme
LOGICAL, INTENT(IN)                  :: OSUBG_COND ! Switch for Subgrid Condensation
                                                   ! (prognotic mode)
LOGICAL, INTENT(IN)                  :: ORAD_SUBG_COND ! Switch for Subgrid Condensation
                                                       ! (diagnostic mode)
!
REAL, DIMENSION(:,:),     INTENT(OUT):: PIRBT  !IR Brightness Temp. (K)
REAL, DIMENSION(:,:),     INTENT(OUT):: PWVBT  !WV Brightness Temp. (K)
!
INTEGER, INTENT(IN)   :: KGEO   !SATELLITE INDEX
REAL, INTENT(IN)                            :: PSIGQSAT ! use an extra "qsat" variance contribution (OSIGMAS case)
!
!*       0.2   DECLARATIONS OF LOCAL VARIABLES
!
LOGICAL :: GPTDEP, GPVOIGT
!
!   reference state
!from inprof
INTEGER :: IGL, ICABS, ING1, IUABS, IINIS, IENDS, ICONF, ICLOUD, IOVLP
INTEGER :: IH2O, ICO2, IO3, ICNT, IN2O, ICH4, ICO, IC11, IC12, ICFC
!
LOGICAL, DIMENSION(KDLON)       :: GDOIT_2D ! .TRUE. for the larger scale
LOGICAL, DIMENSION(KDLON,KFLEV) :: GDOIT  ! .TRUE. for all the levels of the 
                                          !        larger scale columns
!
INTEGER :: JI,JJ,JK,JK1,JK2,JKRAD ! loop indexes
!
INTEGER :: IIB,IIE        ! I index value of the first/last inner mass point
INTEGER :: IJB,IJE        ! J index value of the first/last inner mass point
INTEGER :: IKB,IKE        ! K index value of the first/last inner mass point
INTEGER :: IIU          ! array size for the first  index
INTEGER :: IJU          ! array size for the second index
INTEGER :: IKU          ! array size for the third  index
INTEGER :: IIJ          ! reformatted array index
INTEGER :: IKSTAE       ! level number of the STAndard atmosphere array
INTEGER :: IKUP         ! vertical level above which STAndard atmosphere data
INTEGER :: IDOIT_COL    ! number of larger scale columns
INTEGER :: IDOIT        ! number of levels corresponding of the larger scale 
                        ! columns are filled in
INTEGER :: IDIM         ! effective number of columns for which the radiation
                        !                code is run
INTEGER, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,3)) :: IKKOZ ! indice array used to
               ! vertically interpolate the ozone content on the model grid
!
REAL, DIMENSION(:,:), ALLOCATABLE :: ZTAVE  ! mean-layer temperature
REAL, DIMENSION(:,:), ALLOCATABLE :: ZQVAVE ! mean-layer specific humidity
REAL, DIMENSION(:,:), ALLOCATABLE :: ZO3AVE ! mean-layer ozone content
REAL, DIMENSION(:,:), ALLOCATABLE :: ZPRES_HL ! half-level pressure
REAL, DIMENSION(:,:), ALLOCATABLE :: ZT_HL    ! half-level temperature
REAL, DIMENSION(:,:), ALLOCATABLE :: ZCLDLD ! Downward cloud emissivity
REAL, DIMENSION(:,:), ALLOCATABLE :: ZCLDLU ! Upward cloud emissivity
REAL, DIMENSION(:),   ALLOCATABLE :: ZVIEW  ! cosecant of viewing angle
REAL, DIMENSION(:),   ALLOCATABLE :: ZREMIS ! Reformatted PEMIS array
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZEXNT ! Exner function
REAL, DIMENSION(SIZE(PSTATM,1)) :: ZSTAZZ,ZSTAOZ ! STAndard atmosphere height
                                                 !    and OZone content
REAL :: ZOZ ! variable used to interpolate the ozone profile
!
REAL, DIMENSION(:),   ALLOCATABLE   ::  ZDT0     ! surface discontinuity
REAL, DIMENSION(:,:), ALLOCATABLE   ::  ZRADBT
REAL, DIMENSION(:,:), ALLOCATABLE   ::  ZRADBC
REAL, DIMENSION(:,:), ALLOCATABLE   ::  ZRADFT
REAL, DIMENSION(:),   ALLOCATABLE   ::  ZULAT
REAL, DIMENSION(:),   ALLOCATABLE   ::  ZULON
!
REAL, DIMENSION(:,:), ALLOCATABLE   ::  ZZRADFT
!
REAL, DIMENSION(:),   ALLOCATABLE :: ZWORK1, ZWORK3
!
!  split arrays used to split the memory required by the ECMWF_radiation 
!  subroutine, the fields have the same meaning as their complete counterpart
REAL, DIMENSION(:),     ALLOCATABLE :: ZREMIS_SPLIT
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZO3AVE_SPLIT
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZT_HL_SPLIT
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZPRES_HL_SPLIT
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZQVAVE_SPLIT
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZTAVE_SPLIT
REAL, DIMENSION(:,:), ALLOCATABLE :: ZCLDLD_SPLIT
REAL, DIMENSION(:,:), ALLOCATABLE :: ZCLDLU_SPLIT
REAL, DIMENSION(:),     ALLOCATABLE :: ZVIEW_SPLIT
REAL, DIMENSION(:),   ALLOCATABLE   ::  ZDT0_SPLIT
REAL, DIMENSION(:,:), ALLOCATABLE   ::  ZRADBT_SPLIT
REAL, DIMENSION(:,:), ALLOCATABLE   ::  ZRADBC_SPLIT
!
INTEGER  :: JI_SPLIT          ! loop on the split array
INTEGER  :: INUM_CALL         ! number of CALL of the radiation scheme
INTEGER  :: IDIM_EFF          ! effective number of air-columns to compute
INTEGER  :: IDIM_RESIDUE      ! number of remaining air-columns to compute
INTEGER  :: IBEG, IEND        ! auxiliary indices
!
! Other arrays for emissivity
REAL :: ZFLWP, ZFIWP, ZANGCOR, ZRADLP, ZMULTS, ZTMP, ZKI
!
! Other arrays for condensation
REAL, DIMENSION(:,:,:), ALLOCATABLE  :: ZTEMP  ! Temperature
REAL, DIMENSION(:,:,:), ALLOCATABLE  :: ZSIGRC ! s r_c / sig_s^2
REAL, DIMENSION(:,:,:), ALLOCATABLE  :: ZNCLD  ! grid scale cloud fraction
REAL, DIMENSION(:,:,:), ALLOCATABLE  :: ZRC_IN, ZRC_OUT ! grid scale r_c mixing ratio (kg/kg)
REAL, DIMENSION(:,:,:), ALLOCATABLE  :: ZRI_IN, ZRI_OUT ! grid scale r_i (kg/kg)
REAL, DIMENSION(:,:,:), ALLOCATABLE  :: ZRV_IN, ZRV_OUT ! grid scale r_v (kg/kg)
REAL, DIMENSION(:,:,:), ALLOCATABLE  :: ZRHO
REAL, DIMENSION(SIZE(PPABST,1),SIZE(PPABST,2)) :: ZSIGQSAT2D, ZDUM
TYPE(DIMPHYEX_t)    :: D
!----------------------------------------------------------------------------
!
!*       1.    INITIALIZATION OF CONSTANTS FOR TRANSFERT CODE
!              ----------------------------------------------
!
CALL INIT_NBMOD(KFLEV, IGL, ICABS, ING1, IUABS, IINIS, IENDS, &
       IH2O, ICO2, IO3, ICNT, IN2O, ICH4, ICO, IC11, IC12, ICFC, &
       ICONF, ICLOUD, IOVLP, GPVOIGT, GPTDEP)
X1CO2 = PCCO2 / 44.0 * XMD
!
!----------------------------------------------------------------------------
!
!*       2.    COMPUTE DIMENSIONS OF ARRAYS AND OTHER INDICES
!              ----------------------------------------------
!
IIU = SIZE(PTHT,1)
IJU = SIZE(PTHT,2)
IKU = SIZE(PTHT,3)
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
IKB = 1 + JPVEXT
IKE = IKU - JPVEXT
CALL FILL_DIMPHYEX(D, IIU, IJU, IKU)
!
IKSTAE = SIZE(PSTATM,1)
IKUP   = IKE-JPVEXT+1
!
!----------------------------------------------------------------------------
!
!*       3.    INITIALIZES THE MEAN-LAYER VARIABLES
!              ------------------------------------
!
ALLOCATE(ZEXNT(SIZE(PTHT,1),SIZE(PTHT,2),SIZE(PTHT,3)))
ZEXNT(:,:,:)= ( PPABST(:,:,:)/XP00 ) ** (XRD/XCPD)
! 
ALLOCATE(ZTAVE(KDLON,KFLEV))
ALLOCATE(ZQVAVE(KDLON,KFLEV))
!
ZQVAVE(:,:) = 0.0
!
DO JK=IKB,IKE
  JKRAD = JK-JPVEXT
  DO JJ=IJB,IJE
    DO JI=IIB,IIE
      IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
      ZTAVE(IIJ,JKRAD)  = PTHT(JI,JJ,JK)*ZEXNT(JI,JJ,JK)
    END DO
  END DO
END DO
!
!  Check if the humidity mixing ratio is available
!
IF( SIZE(PRT(:,:,:,:),4) >= 1 ) THEN
  DO JK=IKB,IKE
    JKRAD = JK-JPVEXT
    DO JJ=IJB,IJE
      DO JI=IIB,IIE
        IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
        ZQVAVE(IIJ,JKRAD) = PRT(JI,JJ,JK,1)
      END DO
    END DO
  END DO
END IF
!
!  Standard atmosphere extension
!
DO JK=IKUP,KFLEV
  JK1 = (KSTATM-1)+(JK-IKUP)
  JK2 = JK1+1
  ZTAVE(:,JK)  = 0.5*( PSTATM(JK1,3)+PSTATM(JK2,3) )
  ZQVAVE(:,JK) = 0.5*( PSTATM(JK1,5)/PSTATM(JK1,4)+   &
                       PSTATM(JK2,5)/PSTATM(JK2,4)    )
END DO
!
!----------------------------------------------------------------------------
!
!*       4.    INITIALIZES THE HALF-LEVEL VARIABLES
!  	           ------------------------------------
!
ALLOCATE(ZPRES_HL(KDLON,KFLEV+1))
ALLOCATE(ZT_HL(KDLON,KFLEV+1))
!
DO JK=IKB,IKE+1
  JKRAD = JK-JPVEXT
  DO JJ=IJB,IJE
    DO JI=IIB,IIE
      IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
      ZPRES_HL(IIJ,JKRAD) = XP00 * &
                           (0.5*(ZEXNT(JI,JJ,JK)+ZEXNT(JI,JJ,JK-1)))**(XCPD/XRD)
    END DO
  END DO
END DO
!
!  Standard atmosphere extension
! begining at ikup+1 level allows to use a model domain higher than 50km
!
DO JK=IKUP+1,KFLEV+1
  JK1 = (KSTATM-1)+(JK-IKUP)
  ZPRES_HL(:,JK) = PSTATM(JK1,2)*100.0
END DO
!
!  Surface temperature at the first level
DO JJ=IJB,IJE
  DO JI=IIB,IIE
    IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
    ZT_HL(IIJ,1) = PTSRAD(JI,JJ)
  END DO
END DO
!
!  Temperature at half levels
ZT_HL(:,2:IKE-JPVEXT) = 0.5*(ZTAVE(:,1:IKE-JPVEXT-1)+ZTAVE(:,2:IKE-JPVEXT))
!
DO JJ=IJB,IJE
  DO JI=IIB,IIE
    IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
    ZT_HL(IIJ,IKE-JPVEXT+1)  =  0.5*PTHT(JI,JJ,IKE  )*ZEXNT(JI,JJ,IKE  ) &
                              + 0.5*PTHT(JI,JJ,IKE+1)*ZEXNT(JI,JJ,IKE+1)
  END DO
END DO
!
!  Standard atmosphere extension
! begining at ikup+1 level allows to use a model domain higher than 50km
!
DO JK=IKUP+1,KFLEV+1
  JK1 = (KSTATM-1)+(JK-IKUP)
  ZT_HL(:,JK) = PSTATM(JK1,3)
END DO
!
!----------------------------------------------------------------------------
!
!*       5.    INITIALIZES THE OZONE PROFILES from the standard atmosphere
!	           ------------------------------
!
ALLOCATE(ZO3AVE(KDLON,KFLEV))
!
ZSTAOZ(:) = PSTATM(:,6)/PSTATM(:,4)
ZSTAZZ(:) = 1000.0*PSTATM(:,1)
!
DO JJ = IJB,IJE
  DO JK2 = IKB,IKE
    JKRAD = JK2-JPVEXT
    IKKOZ(:,JK2) = IKB-1
    DO JK1 = 1,IKSTAE
      DO JI = IIB,IIE
        IKKOZ(JI,JK2)=IKKOZ(JI,JK2) + NINT(0.5 + SIGN(0.5,    &
                                          -ZSTAZZ(JK1)+0.5*(PZZ(JI,JJ,JK2)+PZZ(JI,JJ,JK2+1)) ))
      END DO
    END DO
    DO JI = IIB,IIE
      ZOZ=(0.5*(PZZ(JI,JJ,JK2)+PZZ(JI,JJ,JK2+1))- ZSTAZZ(IKKOZ(JI,JK2))) &
           /( ZSTAZZ(IKKOZ(JI,JK2)+1)           - ZSTAZZ(IKKOZ(JI,JK2)))
      IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
      ZO3AVE(IIJ,JKRAD) =( (1.- ZOZ) * ZSTAOZ(IKKOZ(JI,JK2))    &
                         +  ZOZ     * ZSTAOZ(IKKOZ(JI,JK2)+1))
    END DO
  END DO
END DO
!
DO JK=IKUP,KFLEV
  JK1 = (KSTATM)+(JK-IKUP)
  ZO3AVE(:,JK) = ZSTAOZ(JK1)
END DO
!
!----------------------------------------------------------------------------
!
!*       6.    CALLS THE E.C.M.W.F. RADIATION CODE
!	           -----------------------------------
!
!*       6.1   INITIALIZES 2D AND SURFACE FIELDS
!
ALLOCATE(ZREMIS(KDLON))
DO JJ=IJB,IJE
  DO JI=IIB,IIE
    IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
    ZREMIS(IIJ)   = PEMIS(JI,JJ)
  END DO
END DO
!
! initializes surface discontinuity field
ALLOCATE(ZDT0(KDLON))
DO JJ=IJB,IJE
  DO JI=IIB,IIE
    IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
    ZDT0(IIJ) = PTSRAD(JI,JJ) - PTHT(JI,JJ,1)*ZEXNT(JI,JJ,1)
  END DO
END DO
!
ALLOCATE(ZULAT(KDLON))
ALLOCATE(ZULON(KDLON))
DO JJ=IJB,IJE
  DO JI=IIB,IIE
    IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
    ZULON(IIJ) = XLON(JI,JJ)
    ZULAT(IIJ) = XLAT(JI,JJ)
  END DO
END DO
ALLOCATE(ZVIEW(KDLON))
CALL DETER_ANGLE(KGEO, KDLON, ZULAT, ZULON, ZVIEW)
DEALLOCATE(ZULAT)
DEALLOCATE(ZULON)
!
!
ALLOCATE(ZCLDLD(KDLON,KFLEV))
ALLOCATE(ZCLDLU(KDLON,KFLEV))
ZCLDLD = 0.
ZCLDLU = 0.
!
IF( SIZE(PRT(:,:,:,:),4) >= 2 ) THEN
  ALLOCATE(ZNCLD(IIU,IJU,IKU))
  ALLOCATE(ZRC_IN(IIU,IJU,IKU))
  ALLOCATE(ZRC_OUT(IIU,IJU,IKU))
  ZRC_IN=PRT(:,:,:,2)
  ALLOCATE(ZRI_IN(IIU,IJU,IKU))
  ALLOCATE(ZRI_OUT(IIU,IJU,IKU))
  ZRI_IN=0.
  IF( OUSERI ) ZRI_IN=PRT(:,:,:,4)
  IF ( .NOT. OSUBG_COND .AND. ORAD_SUBG_COND) THEN
    PRINT*,' THE SUBGRID CONDENSATION SCHEME IN DIAGNOSTIC MODE IS ACTIVATED'
    ALLOCATE(ZTEMP(IIU,IJU,IKU))
    ZTEMP=PTHT*ZEXNT
    ALLOCATE(ZSIGRC(IIU,IJU,IKU))
    ALLOCATE(ZRV_IN(IIU,IJU,IKU))

    ZRV_IN=PRT(:,:,:,1)
    ALLOCATE(ZRHO(IIU,IJU,IKU))
    ZRHO=1. !unused
    ZSIGQSAT2D(:,:)=PSIGQSAT
    !CALL CONDENSATION( IIU, IJU, IKU, IIB, IIE, IJB, IJE, IKB, IKE, 1, &
    !     'T', 'CB02', 'CB',&
    !     PPABST, PZZ, ZRHO, ZTEMP, ZRV_IN, ZRV_OUT, ZRC_IN, ZRC_OUT, ZRI_IN, ZRI_OUT, &
    !     PRT(:,:,:,2), PRT(:,:,:,5), PRT(:,:,:,6), PSIGS, PMFCONV, ZNCLD, &
    !     ZSIGRC, OUSERI, OSIGMAS, .FALSE., .FALSE., &
    !     ZDUM, ZDUM, ZDUM, ZDUM, ZDUM, ZSIGQSAT2D )
    CALL CONDENSATION(D, CST, RAIN_ICE_PARAM, NEB, TURBN, &                                                                         
                     &'T', 'CB02', 'CB',                                                  &                                         
                     &PPABST, PZZ, ZRHO, ZTEMP, ZRV_IN, ZRV_OUT, ZRC_IN, ZRC_OUT, ZRI_IN, ZRI_OUT,    &                                             
                     &PRT(:,:,:,2), PRT(:,:,:,5), PRT(:,:,:,6), PSIGS, .FALSE., PMFCONV, ZNCLD, ZSIGRC, .FALSE.,                 &     
                     &OSIGMAS, .FALSE., .FALSE.,                                                        &                           
                     &ZDUM, ZDUM, ZDUM, ZDUM, ZDUM, ZSIGQSAT2D)
    DEALLOCATE(ZTEMP,ZSIGRC)
    DEALLOCATE(ZRV_OUT)
  ELSE
    ZNCLD=PCLDFR
  END IF
  DO JK=IKB,IKE-1
    JKRAD = JK-JPVEXT
    DO JJ=IJB,IJE
      DO JI=IIB,IIE
        IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
        IF ( ZVIEW(IIJ) /= XUNDEF .AND. &
             (ZRC_OUT(JI,JJ,JK)  > 0. .OR. ZRI_OUT(JI,JJ,JK) > 0. ) ) THEN
          ZFLWP = ZRC_OUT(JI,JJ,JK) / XG /MAX(1.E-10,ZNCLD(JI,JJ,JK)) &
               * (PPABST(JI,JJ,JK)-PPABST(JI,JJ,JK+1))
          ZFIWP = ZRI_OUT(JI,JJ,JK) / XG /MAX(1.E-10,ZNCLD(JI,JJ,JK)) &
               * (PPABST(JI,JJ,JK)-PPABST(JI,JJ,JK+1))
          ZANGCOR = ZVIEW(IIJ) / 1.66
          !!!Parametrization following Ou and Chou, 1995 (Atmos. Res.)
        ZTMP = ZTAVE(IIJ,JKRAD)-XTT !ZTMP in Celsius degree
        ZRADLP = 326.3+12.42*ZTMP+0.197*(ZTMP**2)+0.0012*(ZTMP**3)
        ZRADLP = MIN(140., MAX(20., ZRADLP))
!!! Parametrization following Ebert and Curry, 1992 (JGR-d)
        ZKI = 0.3 + 1290. / ZRADLP
         ZCLDLD(IIJ,JKRAD) = ZNCLD(JI,JJ,JK)*(1.-EXP &
                             ( -158.*ZFLWP *ZANGCOR-ZKI*ZFIWP*ZVIEW(IIJ)))
         ZCLDLU(IIJ,JKRAD) = ZNCLD(JI,JJ,JK)*(1.-EXP &
                             ( -130.*ZFLWP *ZANGCOR-ZKI*ZFIWP*ZVIEW(IIJ)))
          END IF
        END DO
      END DO
    END DO
  DEALLOCATE(ZNCLD,ZRC_OUT,ZRI_OUT)
END IF
!
DEALLOCATE(ZEXNT)
!
GDOIT_2D(:) = .FALSE.
!
! Flags the columns for which the computations have to be performed
!
DO JJ=IJB,IJE
  DO JI=IIB,IIE
    IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
    IF (ZVIEW(IIJ) /= XUNDEF) GDOIT_2D(IIJ) = .TRUE.
  END DO
END DO
IDOIT_COL = COUNT( GDOIT_2D(:) )  ! number of larger scale columns
!
GDOIT(:,:) = SPREAD( GDOIT_2D(:),DIM=2,NCOPIES=KFLEV )
IDOIT = IDOIT_COL*KFLEV
ALLOCATE(ZWORK1(IDOIT))
!
! temperature profiles
ZWORK1(:) = PACK( ZTAVE(:,:),MASK=GDOIT(:,:) )
DEALLOCATE(ZTAVE)
ALLOCATE(ZTAVE(IDOIT_COL,KFLEV))
ZTAVE(:,:) = RESHAPE( ZWORK1(:),(/IDOIT_COL,KFLEV/) )
!
! vapor mixing ratio profiles
ZWORK1(:) = PACK( ZQVAVE(:,:),MASK=GDOIT(:,:) )
DEALLOCATE(ZQVAVE)
ALLOCATE(ZQVAVE(IDOIT_COL,KFLEV))
ZQVAVE(:,:) = RESHAPE( ZWORK1(:),(/IDOIT_COL,KFLEV/) )
!
! cloud emissivities
ZWORK1(:) = PACK( ZCLDLD(:,:),MASK=GDOIT(:,:) )
DEALLOCATE(ZCLDLD)
ALLOCATE(ZCLDLD(IDOIT_COL,KFLEV))
ZCLDLD(:,:) = RESHAPE( ZWORK1(:),(/IDOIT_COL,KFLEV/) )
!
ZWORK1(:) = PACK( ZCLDLU(:,:),MASK=GDOIT(:,:) )
DEALLOCATE(ZCLDLU)
ALLOCATE(ZCLDLU(IDOIT_COL,KFLEV))
ZCLDLU(:,:) = RESHAPE( ZWORK1(:),(/IDOIT_COL,KFLEV/) )
!
! ozone content profiles
ZWORK1(:) = PACK( ZO3AVE(:,:),MASK=GDOIT(:,:) )
DEALLOCATE(ZO3AVE)
ALLOCATE(ZO3AVE(IDOIT_COL,KFLEV))
ZO3AVE(:,:) = RESHAPE( ZWORK1(:),(/IDOIT_COL,KFLEV/) )
!
! half-level variables
ZWORK1(:) = PACK(  ZPRES_HL(:,1:KFLEV),MASK=GDOIT(:,:) )
DEALLOCATE(ZPRES_HL)
ALLOCATE(ZPRES_HL(IDOIT_COL,KFLEV+1))
ZPRES_HL(:,1:KFLEV) = RESHAPE( ZWORK1(:),(/IDOIT_COL,KFLEV/) )
ZPRES_HL(:,KFLEV+1) = PSTATM(IKSTAE,2)*100.0
!
ZWORK1(:) = PACK(  ZT_HL(:,1:KFLEV),MASK=GDOIT(:,:) )
DEALLOCATE(ZT_HL)
ALLOCATE(ZT_HL(IDOIT_COL,KFLEV+1))
ZT_HL(:,1:KFLEV) = RESHAPE( ZWORK1(:),(/IDOIT_COL,KFLEV/) )
ZT_HL(:,KFLEV+1) = PSTATM(IKSTAE,3)
!
! surface fields
ALLOCATE(ZWORK3(IDOIT_COL))
ZWORK3(:) = PACK( ZVIEW(:),MASK=GDOIT_2D(:) )
DEALLOCATE(ZVIEW)
ALLOCATE(ZVIEW(IDOIT_COL))
ZVIEW(:) = ZWORK3(:)
!
ZWORK3(:) = PACK( ZREMIS(:),MASK=GDOIT_2D(:) )
DEALLOCATE(ZREMIS)
ALLOCATE(ZREMIS(IDOIT_COL))
ZREMIS(:) = ZWORK3(:)
!
ZWORK3(:) = PACK( ZDT0(:),MASK=GDOIT_2D(:) )
DEALLOCATE(ZDT0)
ALLOCATE(ZDT0(IDOIT_COL))
ZDT0(:) = ZWORK3(:)
!
DEALLOCATE(ZWORK1)
DEALLOCATE(ZWORK3)
!
! radiation fields
ALLOCATE(ZRADBC(IDOIT_COL,JPWVINT))
ALLOCATE(ZRADBT(IDOIT_COL,JPWVINT))
!
IDIM = IDOIT_COL
PRINT *,'KGEO =',KGEO,' IDIM =',IDIM
!
!*       6.2   CALLS THE ECMWF_RADIATION ROUTINES
!
!  ***********************************************************
!  *CAUTION: Routine  nbmvec  is written in FORTRAN 77*
!  ***********************************************************
!
!  mixing ratio -> specific humidity conversion
ZQVAVE(:,:) = ZQVAVE(:,:) / (1.+ZQVAVE(:,:))
!
IF( IDIM <= KRAD_COLNBR ) THEN 
   !
   !  there is less than KRAD_COLNBR verticals to be considered therefore
   ! no split of the arrays is performed
   !
   CALL NBMVEC( 1, IDIM, IDIM, KFLEV, IGL, ICABS, ING1, IUABS, &
        IH2O, ICO2, IO3, ICNT, IN2O, ICH4, ICO, IC11, IC12, ICFC, &
        IINIS, IENDS, ICONF, ICLOUD, IOVLP, GPVOIGT, GPTDEP, &
        ZTAVE, ZQVAVE, ZO3AVE, ZPRES_HL, ZT_HL, &
        ZVIEW, ZCLDLD, ZCLDLU, ZDT0, ZREMIS, ZRADBC, ZRADBT)
ELSE
   !
   ! the splitting of the arrays will be performed
   !
   INUM_CALL = CEILING( REAL( IDIM ) / REAL( KRAD_COLNBR ) )
   IDIM_RESIDUE = IDIM
   DO JI_SPLIT = 1 , INUM_CALL
     IDIM_EFF = MIN( IDIM_RESIDUE,KRAD_COLNBR )
     !
     IF( JI_SPLIT == 1 .OR. JI_SPLIT == INUM_CALL ) THEN
       ALLOCATE(  ZREMIS_SPLIT(IDIM_EFF))
       ALLOCATE(  ZO3AVE_SPLIT(IDIM_EFF,KFLEV))
       ALLOCATE(  ZT_HL_SPLIT(IDIM_EFF,KFLEV+1))
       ALLOCATE(  ZPRES_HL_SPLIT(IDIM_EFF,KFLEV+1))
       ALLOCATE(  ZQVAVE_SPLIT(IDIM_EFF,KFLEV))
       ALLOCATE(  ZTAVE_SPLIT(IDIM_EFF,KFLEV))
       ALLOCATE(  ZCLDLU_SPLIT(IDIM_EFF,KFLEV))
       ALLOCATE(  ZCLDLD_SPLIT(IDIM_EFF,KFLEV))
       ALLOCATE(  ZVIEW_SPLIT(IDIM_EFF))
       ALLOCATE(  ZDT0_SPLIT(IDIM_EFF))
       ALLOCATE(  ZRADBT_SPLIT(IDIM_EFF,JPWVINT))
       ALLOCATE(  ZRADBC_SPLIT(IDIM_EFF,JPWVINT))
     END IF
     !
     ! fill the split arrays with their values
     ! taken from the full arrays 
     !
     IBEG = IDIM-IDIM_RESIDUE+1
     IEND = IBEG+IDIM_EFF-1
     ZREMIS_SPLIT(:)   = ZREMIS( IBEG:IEND )
     ZO3AVE_SPLIT(:,:) = ZO3AVE( IBEG:IEND ,:)
     ZT_HL_SPLIT(:,:)    = ZT_HL( IBEG:IEND ,:)
     ZPRES_HL_SPLIT(:,:) = ZPRES_HL( IBEG:IEND ,:)
     ZQVAVE_SPLIT(:,:) = ZQVAVE( IBEG:IEND ,:)
     ZTAVE_SPLIT(:,:)  = ZTAVE ( IBEG:IEND ,:)
     ZCLDLU_SPLIT(:,:)  = ZCLDLU ( IBEG:IEND ,:)
     ZCLDLD_SPLIT(:,:)  = ZCLDLD ( IBEG:IEND ,:)
     ZVIEW_SPLIT(:)    = ZVIEW ( IBEG:IEND )
     ZDT0_SPLIT(:)    = ZDT0 ( IBEG:IEND )
     !  
     ! call ECMWF_radiation with the split arrays
     !
     CALL NBMVEC( 1, IDIM_EFF, IDIM_EFF, KFLEV, IGL, ICABS, ING1, IUABS,&
          IH2O, ICO2, IO3, ICNT, IN2O, ICH4, ICO, IC11, IC12, ICFC, &
          IINIS, IENDS, ICONF, ICLOUD, IOVLP, GPVOIGT, GPTDEP, &
          ZTAVE_SPLIT, ZQVAVE_SPLIT, ZO3AVE_SPLIT, &
          ZPRES_HL_SPLIT, ZT_HL_SPLIT, &
          ZVIEW_SPLIT, ZCLDLD_SPLIT, ZCLDLU_SPLIT, ZDT0_SPLIT, &
          ZREMIS_SPLIT, ZRADBC_SPLIT, ZRADBT_SPLIT)
     !
     ! fill the full output arrays with the split arrays
     !
     ZRADBT( IBEG:IEND ,:)  = ZRADBT_SPLIT(:,:)  
     ZRADBC( IBEG:IEND ,:)  = ZRADBC_SPLIT(:,:)  
     !
     IDIM_RESIDUE = IDIM_RESIDUE - IDIM_EFF
     !
     ! desallocation of the split arrays
     !
     IF( JI_SPLIT >= INUM_CALL-1 ) THEN
       DEALLOCATE(ZREMIS_SPLIT)
       DEALLOCATE(ZO3AVE_SPLIT)
       DEALLOCATE(ZT_HL_SPLIT)
       DEALLOCATE(ZPRES_HL_SPLIT)
       DEALLOCATE(ZQVAVE_SPLIT)
       DEALLOCATE(ZTAVE_SPLIT)
       DEALLOCATE(ZCLDLU_SPLIT)
       DEALLOCATE(ZCLDLD_SPLIT)
       DEALLOCATE(ZVIEW_SPLIT)
       DEALLOCATE(ZDT0_SPLIT)
       DEALLOCATE(ZRADBT_SPLIT)
       DEALLOCATE(ZRADBC_SPLIT)
     END IF
   END DO
END IF
!
DEALLOCATE(ZTAVE,ZQVAVE,ZO3AVE)
DEALLOCATE(ZPRES_HL,ZT_HL)
DEALLOCATE(ZREMIS)
DEALLOCATE(ZDT0)
DEALLOCATE(ZCLDLD,ZCLDLU)
DEALLOCATE(ZVIEW)
!
ZRADBT = ZRADBT / XPI
ALLOCATE(ZRADFT(IDIM,JPCAN))
CALL MAKE_RADSAT(KYEARF, KMONTHF, KDAYF, PSECF, &
                 KGEO, IDIM, ZRADBT, ZRADFT)
DEALLOCATE(ZRADBT)
DEALLOCATE(ZRADBC)
!
ALLOCATE(ZWORK1(IDIM*JPCAN))
ZWORK1(:) = PACK( ZRADFT(:,:),MASK=.TRUE. )
ALLOCATE(ZZRADFT(KDLON,JPCAN))
ZZRADFT(:,:) = UNPACK( ZWORK1(:),MASK=GDOIT(:,1:JPCAN),FIELD=XUNDEF )
DEALLOCATE(ZRADFT)
DEALLOCATE(ZWORK1)
!
PIRBT = XUNDEF
PWVBT = XUNDEF
DO JJ=IJB,IJE
  DO JI=IIB,IIE
    IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
    PIRBT(JI,JJ) = ZZRADFT(IIJ,1)
    PWVBT(JI,JJ) = ZZRADFT(IIJ,2)
  END DO
END DO
DEALLOCATE(ZZRADFT)
!
END SUBROUTINE RADTR_SATEL
