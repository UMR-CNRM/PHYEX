!MNH_LIC Copyright 2013-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!      ######################
       MODULE MODI_LIMA_MIXED
!      ######################
!
IMPLICIT NONE
INTERFACE
      SUBROUTINE LIMA_MIXED (OSEDI, OHHONI, KSPLITG, PTSTEP, KMI, &
                             KRR, PZZ, PRHODJ,                    &
                             PRHODREF, PEXNREF, PPABST, PW_NU,    &
                             PTHT, PRT, PSVT,                     &
                             PTHS, PRS, PSVS)
!
USE MODD_NSV,   only: NSV_LIMA_BEG
IMPLICIT NONE
!
LOGICAL,                  INTENT(IN)    :: OSEDI   ! switch to activate the
                                                   ! cloud ice sedimentation
LOGICAL,                  INTENT(IN)    :: OHHONI  ! enable haze freezing
INTEGER,                  INTENT(IN)    :: KSPLITG ! Number of small time step 
                                      ! integration for  ice sedimendation
REAL,                     INTENT(IN)    :: PTSTEP  ! Time step          
INTEGER,                  INTENT(IN)    :: KMI     ! Model index 
!
INTEGER,                  INTENT(IN)    :: KRR     ! Number of moist variables
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PZZ     ! Height (z)
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PEXNREF ! Reference Exner function
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABST  ! abs. pressure at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PW_NU   ! updraft velocity used for
                                                   ! the nucleation param.
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT    ! Theta at time t
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRT     ! m.r. at t 
REAL, DIMENSION(:,:,:,NSV_LIMA_BEG:), INTENT(IN)    :: PSVT ! Concentrations at time t
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PTHS    ! Theta source
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRS     ! m.r. source
REAL, DIMENSION(:,:,:,NSV_LIMA_BEG:), INTENT(INOUT) :: PSVS ! Concentration sources
!
END SUBROUTINE LIMA_MIXED
END INTERFACE
END MODULE MODI_LIMA_MIXED
!
!     #######################################################################
      SUBROUTINE LIMA_MIXED (OSEDI, OHHONI, KSPLITG, PTSTEP, KMI, &
                             KRR, PZZ, PRHODJ,                    &
                             PRHODREF, PEXNREF, PPABST, PW_NU,    &
                             PTHT, PRT, PSVT,                     &
                             PTHS, PRS, PSVS                      )
!     #######################################################################
!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the mixed-phase 
!!    microphysical processes
!!
!!
!!**  METHOD
!!    ------
!!
!!    REFERENCE
!!    ---------
!!
!!      Most of the parameterizations come from the ICE3 scheme, described in
!!    the MESO-NH scientific documentation.
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
!  P. Wautelet 05/2016-04/2018: new data structures and calls for I/O
!  P. Wautelet    03/2020: use the new data structures and subroutines for budgets (no more call to budget in this subroutine)
!  P. Wautelet 28/05/2020: bugfix: correct array start for PSVT and PSVS
!  P. Wautelet 02/02/2021: budgets: add missing source terms for SV budgets in LIMA
!  J. Wurtz       03/2022: new snow characteristics
!  M. Taufour     07/2022: add concentration for snow, graupel, hail
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
use modd_budget,           only: lbu_enable, nbumod
USE MODD_CST,              ONLY: XP00, XRD, XRV, XMV, XMD, XCPD, XCPV,       &
                                 XCL, XCI, XTT, XLSTT, XLVTT,                &
                                 XALPI, XBETAI, XGAMI
USE MODD_NSV
USE MODD_PARAMETERS,       ONLY: JPHEXT, JPVEXT
USE MODD_PARAM_LIMA,       ONLY: NMOD_IFN, XRTMIN, XCTMIN,     &
                                 NMOD_CCN, NMOD_IMM, LSNOW_T,    &
                                 NMOM_C, NMOM_R, NMOM_I, NMOM_S, NMOM_G, NMOM_H
USE MODD_PARAM_LIMA_WARM,  ONLY: XLBC, XLBEXC, XLBR, XLBEXR
USE MODD_PARAM_LIMA_COLD,  ONLY: XLBI, XLBEXI, XLBS, XLBEXS, XSCFAC, &
                                 XLBDAS_MAX, XLBDAS_MIN, XTRANS_MP_GAMMAS, &
                                 XBS, XCCS, XCXS, XNS
USE MODD_PARAM_LIMA_MIXED, ONLY: XLBG, XLBEXG, XCCG, XCXG, XLBH, XLBEXH, XCCH, XCXH

use mode_tools,            only: Countjv

USE MODI_LIMA_MIXED_FAST_PROCESSES
USE MODI_LIMA_MIXED_SLOW_PROCESSES

IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
LOGICAL,                  INTENT(IN)    :: OSEDI   ! switch to activate the 
                                                   ! cloud ice sedimentation
LOGICAL,                  INTENT(IN)    :: OHHONI  ! enable haze freezing
INTEGER,                  INTENT(IN)    :: KSPLITG ! Number of small time step 
                                      ! integration for  ice sedimendation
REAL,                     INTENT(IN)    :: PTSTEP  ! Time step          
INTEGER,                  INTENT(IN)    :: KMI     ! Model index 
!
INTEGER,                  INTENT(IN)    :: KRR     ! Number of moist variables
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PZZ     ! Height (z)
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PEXNREF ! Reference Exner function
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABST  ! abs. pressure at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PW_NU   ! updraft velocity used for
                                                   ! the nucleation param.
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT    ! Theta at time t
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRT     ! m.r. at t 
REAL, DIMENSION(:,:,:,NSV_LIMA_BEG:), INTENT(IN)    :: PSVT ! Concentrations at time t
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PTHS    ! Theta source
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRS     ! m.r. source
REAL, DIMENSION(:,:,:,NSV_LIMA_BEG:), INTENT(INOUT) :: PSVS ! Concentration sources
!
!*       0.2   Declarations of local variables :
!
!3D microphysical variables
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3))  &
                                    :: PRVT,    & ! Water vapor m.r. at t 
                                       PRCT,    & ! Cloud water m.r. at t 
                                       PRRT,    & ! Rain water m.r. at t 
                                       PRIT,    & ! Cloud ice m.r. at t 
                                       PRST,    & ! Snow/aggregate m.r. at t 
                                       PRGT,    & ! Graupel m.r. at t 
                                       PRHT,    & ! Hail m.r. at t 
                                       !
                                       PRVS,    & ! Water vapor m.r. source
                                       PRCS,    & ! Cloud water m.r. source
                                       PRRS,    & ! Rain water m.r. source
                                       PRIS,    & ! Pristine ice m.r. source
                                       PRSS,    & ! Snow/aggregate m.r. source
                                       PRGS,    & ! Graupel m.r. source
                                       PRHS,    & ! Hail m.r. source
                                       !
                                       PCCT,    & ! Cloud water C. at t
                                       PCRT,    & ! Rain water C. at t
                                       PCIT,    & ! Ice crystal C. at t
                                       PCST,    & ! Snow/aggregates C. at t   
                                       PCGT,    & ! Graupel C. at t           
                                       PCHT,    & ! Hail C. at t 
                                       !
                                       PCCS,    & ! Cloud water C. source
                                       PCRS,    & ! Rain water C. source
                                       PCIS,    & ! Ice crystal C. source
                                       PCSS,    & ! Snow/aggregates C. source
                                       PCGS,    & ! Graupel C. source
                                       PCHS       ! Hail C. source
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: PNFS     ! CCN C. available source
                                                  !used as Free ice nuclei for
                                                  !HOMOGENEOUS nucleation of haze
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: PNAS     ! Cloud  C. nuclei C. source
                                                  !used as Free ice nuclei for
                                                  !IMMERSION freezing
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: PIFS     ! Free ice nuclei C. source 
                                                  !for DEPOSITION and CONTACT
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: PINS     ! Activated ice nuclei C. source
                                                  !for DEPOSITION and CONTACT
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: PNIS     ! Activated ice nuclei C. source
                                                  !for IMMERSION
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: PNHS     ! Hom. freezing of CCN
!
! Replace PACK
LOGICAL, DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3)) :: GMICRO
INTEGER :: IMICRO
INTEGER , DIMENSION(SIZE(GMICRO)) :: I1,I2,I3 ! Used to replace the COUNT
INTEGER                           :: JL       ! and PACK intrinsics
!
! Packed microphysical variables
REAL, DIMENSION(:), ALLOCATABLE   :: ZRVT    ! Water vapor m.r. at t
REAL, DIMENSION(:), ALLOCATABLE   :: ZRCT    ! Cloud water m.r. at t
REAL, DIMENSION(:), ALLOCATABLE   :: ZRRT    ! Rain water m.r. at t
REAL, DIMENSION(:), ALLOCATABLE   :: ZRIT    ! Pristine ice m.r. at t
REAL, DIMENSION(:), ALLOCATABLE   :: ZRST    ! Snow/aggregate m.r. at t
REAL, DIMENSION(:), ALLOCATABLE   :: ZRGT    ! Graupel m.r. at t
REAL, DIMENSION(:), ALLOCATABLE   :: ZRHT    ! Hail m.r. at t
REAL, DIMENSION(:), ALLOCATABLE   :: ZCCT    ! Cloud water conc. at t
REAL, DIMENSION(:), ALLOCATABLE   :: ZCRT    ! Rain water conc. at t
REAL, DIMENSION(:), ALLOCATABLE   :: ZCIT    ! Pristine ice conc. at t
REAL, DIMENSION(:), ALLOCATABLE   :: ZCST    ! Snow/aggregates conc. at t 
REAL, DIMENSION(:), ALLOCATABLE   :: ZCGT    ! Graupel conc. at t
REAL, DIMENSION(:), ALLOCATABLE   :: ZCHT    ! Hail conc. at t
!
REAL, DIMENSION(:), ALLOCATABLE   :: ZRVS    ! Water vapor m.r. source
REAL, DIMENSION(:), ALLOCATABLE   :: ZRCS    ! Cloud water m.r. source
REAL, DIMENSION(:), ALLOCATABLE   :: ZRRS    ! Rain water m.r. source
REAL, DIMENSION(:), ALLOCATABLE   :: ZRIS    ! Pristine ice m.r. source
REAL, DIMENSION(:), ALLOCATABLE   :: ZRSS    ! Snow/aggregate m.r. source
REAL, DIMENSION(:), ALLOCATABLE   :: ZRGS    ! Graupel m.r. source
REAL, DIMENSION(:), ALLOCATABLE   :: ZRHS    ! Hail m.r. source
!
REAL, DIMENSION(:), ALLOCATABLE   :: ZTHS    ! Theta source
!
REAL, DIMENSION(:),   ALLOCATABLE :: ZCCS    ! Cloud water conc. source
REAL, DIMENSION(:),   ALLOCATABLE :: ZCRS    ! Rain water conc. source
REAL, DIMENSION(:),   ALLOCATABLE :: ZCIS    ! Pristine ice conc. source
REAL, DIMENSION(:),   ALLOCATABLE :: ZCSS    ! Snow/aggregates conc. at t
REAL, DIMENSION(:),   ALLOCATABLE :: ZCGS    ! Graupel conc. at t
REAL, DIMENSION(:),   ALLOCATABLE :: ZCHS    ! Hail conc. at t
!
REAL, DIMENSION(:,:), ALLOCATABLE :: ZIFS    ! Free Ice nuclei conc. source
REAL, DIMENSION(:,:), ALLOCATABLE :: ZINS    ! Nucleated Ice nuclei conc. source
!
! Other packed variables
REAL, DIMENSION(:), ALLOCATABLE &
                   :: ZRHODREF, & ! RHO Dry REFerence
                      ZRHODJ,   & ! RHO times Jacobian
                      ZZT,      & ! Temperature
                      ZPRES,    & ! Pressure
                      ZEXNREF,  & ! EXNer Pressure REFerence
                      ZZW,      & ! Work array
                      ZLSFACT,  & ! L_s/(Pi_ref*C_ph)
                      ZLVFACT,  & ! L_v/(Pi_ref*C_ph)
                      ZSSI,     & ! Supersaturation over ice
                      ZLBDAC,   & ! Slope parameter of the cloud droplet distr.
                      ZLBDAR,   & ! Slope parameter of the raindrop  distr.
                      ZLBDAI,   & ! Slope parameter of the ice crystal distr.
                      ZLBDAS,   & ! Slope parameter of the aggregate distr.
                      ZLBDAG,   & ! Slope parameter of the graupel   distr.
                      ZLBDAH,   & ! Slope parameter of the hail   distr.
                      ZAI,      & ! Thermodynamical function
                      ZCJ,      & ! used to compute the ventilation coefficient
                      ZKA,      & ! Thermal conductivity of the air
                      ZDV         ! Diffusivity of water vapor in the air
!
! 3D Temperature
REAL,    DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3)) :: ZT, ZW
!
!
INTEGER :: IIB, IIE, IJB, IJE, IKB, IKE        ! Physical domain
INTEGER :: JMOD_IFN                            ! Loop index 
!
!-------------------------------------------------------------------------------
!
!
!*       0.     3D MICROPHYSCAL VARIABLES
!	        -------------------------
!
!
! Prepare 3D water mixing ratios
PRVT(:,:,:) = PRT(:,:,:,1)
PRVS(:,:,:) = PRS(:,:,:,1)
!
PRCT(:,:,:) = 0.
PRCS(:,:,:) = 0.
PRRT(:,:,:) = 0.
PRRS(:,:,:) = 0.
PRIT(:,:,:) = 0.
PRIS(:,:,:) = 0.
PRST(:,:,:) = 0.
PRSS(:,:,:) = 0.
PRGT(:,:,:) = 0.
PRGS(:,:,:) = 0.
PRHT(:,:,:) = 0.
PRHS(:,:,:) = 0.
!
IF ( KRR .GE. 2 ) PRCT(:,:,:) = PRT(:,:,:,2)
IF ( KRR .GE. 2 ) PRCS(:,:,:) = PRS(:,:,:,2)
IF ( KRR .GE. 3 ) PRRT(:,:,:) = PRT(:,:,:,3)
IF ( KRR .GE. 3 ) PRRS(:,:,:) = PRS(:,:,:,3)
IF ( KRR .GE. 4 ) PRIT(:,:,:) = PRT(:,:,:,4)
IF ( KRR .GE. 4 ) PRIS(:,:,:) = PRS(:,:,:,4)
IF ( KRR .GE. 5 ) PRST(:,:,:) = PRT(:,:,:,5) 
IF ( KRR .GE. 5 ) PRSS(:,:,:) = PRS(:,:,:,5)
IF ( KRR .GE. 6 ) PRGT(:,:,:) = PRT(:,:,:,6)
IF ( KRR .GE. 6 ) PRGS(:,:,:) = PRS(:,:,:,6)
IF ( KRR .GE. 7 ) PRHT(:,:,:) = PRT(:,:,:,7)
IF ( KRR .GE. 7 ) PRHS(:,:,:) = PRS(:,:,:,7)
!
! Prepare 3D number concentrations
PCCT(:,:,:) = 0.
PCRT(:,:,:) = 0.
PCIT(:,:,:) = 0.
PCST(:,:,:) = 0. 
PCGT(:,:,:) = 0. 
PCHT(:,:,:) = 0.
PCCS(:,:,:) = 0.
PCRS(:,:,:) = 0.
PCIS(:,:,:) = 0.
PCSS(:,:,:) = 0. 
PCGS(:,:,:) = 0. 
PCHS(:,:,:) = 0. 
!
IF ( NMOM_C.GE.2 ) PCCT(:,:,:) = PSVT(:,:,:,NSV_LIMA_NC) 
IF ( NMOM_R.GE.2 ) PCRT(:,:,:) = PSVT(:,:,:,NSV_LIMA_NR)
IF ( NMOM_I.GE.2 ) PCIT(:,:,:) = PSVT(:,:,:,NSV_LIMA_NI)
IF ( NMOM_S.GE.2 ) PCST(:,:,:) = PSVT(:,:,:,NSV_LIMA_NS)
IF ( NMOM_G.GE.2 ) PCGT(:,:,:) = PSVT(:,:,:,NSV_LIMA_NG)
IF ( NMOM_H.GE.2 ) PCHT(:,:,:) = PSVT(:,:,:,NSV_LIMA_NH)
!
IF ( NMOM_C.GE.2 ) PCCS(:,:,:) = PSVS(:,:,:,NSV_LIMA_NC)
IF ( NMOM_R.GE.2 ) PCRS(:,:,:) = PSVS(:,:,:,NSV_LIMA_NR)
IF ( NMOM_I.GE.2 ) PCIS(:,:,:) = PSVS(:,:,:,NSV_LIMA_NI)
IF ( NMOM_S.GE.2 ) PCSS(:,:,:) = PSVS(:,:,:,NSV_LIMA_NS)
IF ( NMOM_G.GE.2 ) PCGS(:,:,:) = PSVS(:,:,:,NSV_LIMA_NG)
IF ( NMOM_H.GE.2 ) PCHS(:,:,:) = PSVS(:,:,:,NSV_LIMA_NH)
!
IF ( NMOD_CCN .GE. 1 ) THEN
   ALLOCATE( PNFS(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3),NMOD_CCN) )
   ALLOCATE( PNAS(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3),NMOD_CCN) )
   PNFS(:,:,:,:) = PSVS(:,:,:,NSV_LIMA_CCN_FREE:NSV_LIMA_CCN_FREE+NMOD_CCN-1)
   PNAS(:,:,:,:) = PSVS(:,:,:,NSV_LIMA_CCN_ACTI:NSV_LIMA_CCN_ACTI+NMOD_CCN-1)
ELSE
   ALLOCATE( PNFS(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3),1) )
   ALLOCATE( PNAS(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3),1) )
   PNFS(:,:,:,:) = 0.
   PNAS(:,:,:,:) = 0.
END IF
!
IF ( NMOD_IFN .GE. 1 ) THEN
   ALLOCATE( PIFS(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3),NMOD_IFN) )
   ALLOCATE( PINS(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3),NMOD_IFN) )
   PIFS(:,:,:,:) = PSVS(:,:,:,NSV_LIMA_IFN_FREE:NSV_LIMA_IFN_FREE+NMOD_IFN-1)
   PINS(:,:,:,:) = PSVS(:,:,:,NSV_LIMA_IFN_NUCL:NSV_LIMA_IFN_NUCL+NMOD_IFN-1)
ELSE
   ALLOCATE( PIFS(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3),1) )
   ALLOCATE( PINS(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3),1) )
   PIFS(:,:,:,:) = 0.
   PINS(:,:,:,:) = 0.
END IF
!
IF ( NMOD_IMM .GE. 1 ) THEN
   ALLOCATE( PNIS(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3),NMOD_IMM) )
   PNIS(:,:,:,:) = PSVS(:,:,:,NSV_LIMA_IMM_NUCL:NSV_LIMA_IMM_NUCL+NMOD_IMM-1)
ELSE
   ALLOCATE( PNIS(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3),1) )
   PNIS(:,:,:,:) = 0.0
END IF
!
IF ( OHHONI ) THEN
   ALLOCATE( PNHS(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3)) )
   PNHS(:,:,:) = PSVS(:,:,:,NSV_LIMA_HOM_HAZE)
ELSE
   ALLOCATE( PNHS(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3)) )
   PNHS(:,:,:) = 0.0
END IF
!
!-------------------------------------------------------------------------------
!
!
!*       1.     Pack variables, computations only where necessary
!	        -------------------------------------------------
!
! Physical domain
!
IIB=1+JPHEXT
IIE=SIZE(PZZ,1) - JPHEXT
IJB=1+JPHEXT
IJE=SIZE(PZZ,2) - JPHEXT
IKB=1+JPVEXT
IKE=SIZE(PZZ,3) - JPVEXT
!
! Temperature
ZT(:,:,:)  = PTHT(:,:,:) * ( PPABST(:,:,:)/XP00 ) ** (XRD/XCPD)
!
! Looking for regions where computations are necessary
GMICRO(:,:,:) = .FALSE.
GMICRO(IIB:IIE,IJB:IJE,IKB:IKE) = PRCT(IIB:IIE,IJB:IJE,IKB:IKE)>XRTMIN(2) .OR. &
                                  PRRT(IIB:IIE,IJB:IJE,IKB:IKE)>XRTMIN(3) .OR. &
                                  PRIT(IIB:IIE,IJB:IJE,IKB:IKE)>XRTMIN(4) .OR. &
                                  PRST(IIB:IIE,IJB:IJE,IKB:IKE)>XRTMIN(5) .OR. &
                                  PRGT(IIB:IIE,IJB:IJE,IKB:IKE)>XRTMIN(6) .OR. &
                                  PRHT(IIB:IIE,IJB:IJE,IKB:IKE)>XRTMIN(7)
!
IMICRO = COUNTJV( GMICRO(:,:,:),I1(:),I2(:),I3(:))
!
IF( IMICRO >= 0 ) THEN
!
   ALLOCATE(ZRVT(IMICRO)) 
   ALLOCATE(ZRCT(IMICRO))
   ALLOCATE(ZRRT(IMICRO))  
   ALLOCATE(ZRIT(IMICRO)) 
   ALLOCATE(ZRST(IMICRO)) 
   ALLOCATE(ZRGT(IMICRO))  
   ALLOCATE(ZRHT(IMICRO))  
   !
   ALLOCATE(ZCCT(IMICRO)) 
   ALLOCATE(ZCRT(IMICRO)) 
   ALLOCATE(ZCIT(IMICRO))
   ALLOCATE(ZCST(IMICRO)) 
   ALLOCATE(ZCGT(IMICRO))
   ALLOCATE(ZCHT(IMICRO))
   !
   ALLOCATE(ZRVS(IMICRO))  
   ALLOCATE(ZRCS(IMICRO)) 
   ALLOCATE(ZRRS(IMICRO)) 
   ALLOCATE(ZRIS(IMICRO))
   ALLOCATE(ZRSS(IMICRO))
   ALLOCATE(ZRGS(IMICRO)) 
   ALLOCATE(ZRHS(IMICRO)) 
   ALLOCATE(ZTHS(IMICRO))
   !
   ALLOCATE(ZCCS(IMICRO)) 
   ALLOCATE(ZCRS(IMICRO)) 
   ALLOCATE(ZCIS(IMICRO))
   ALLOCATE(ZCSS(IMICRO)) 
   ALLOCATE(ZCGS(IMICRO)) 
   ALLOCATE(ZCHS(IMICRO))
   ALLOCATE(ZIFS(IMICRO,NMOD_IFN))
   ALLOCATE(ZINS(IMICRO,NMOD_IFN))
   !
   ALLOCATE(ZRHODREF(IMICRO)) 
   ALLOCATE(ZZT(IMICRO)) 
   ALLOCATE(ZPRES(IMICRO)) 
   ALLOCATE(ZEXNREF(IMICRO))
   DO JL=1,IMICRO   
      ZRVT(JL) = PRVT(I1(JL),I2(JL),I3(JL))
      ZRCT(JL) = PRCT(I1(JL),I2(JL),I3(JL))
      ZRRT(JL) = PRRT(I1(JL),I2(JL),I3(JL))
      ZRIT(JL) = PRIT(I1(JL),I2(JL),I3(JL))
      ZRST(JL) = PRST(I1(JL),I2(JL),I3(JL))
      ZRGT(JL) = PRGT(I1(JL),I2(JL),I3(JL))
      ZRHT(JL) = PRHT(I1(JL),I2(JL),I3(JL))
      !
      ZCCT(JL) = PCCT(I1(JL),I2(JL),I3(JL))
      ZCRT(JL) = PCRT(I1(JL),I2(JL),I3(JL))
      ZCIT(JL) = PCIT(I1(JL),I2(JL),I3(JL))
      ZCST(JL) = PCST(I1(JL),I2(JL),I3(JL))
      ZCGT(JL) = PCGT(I1(JL),I2(JL),I3(JL)) 
      ZCHT(JL) = PCHT(I1(JL),I2(JL),I3(JL))
      !
      ZRVS(JL) = PRVS(I1(JL),I2(JL),I3(JL))
      ZRCS(JL) = PRCS(I1(JL),I2(JL),I3(JL))
      ZRRS(JL) = PRRS(I1(JL),I2(JL),I3(JL))
      ZRIS(JL) = PRIS(I1(JL),I2(JL),I3(JL))
      ZRSS(JL) = PRSS(I1(JL),I2(JL),I3(JL))
      ZRGS(JL) = PRGS(I1(JL),I2(JL),I3(JL))
      ZRHS(JL) = PRHS(I1(JL),I2(JL),I3(JL))
      ZTHS(JL) = PTHS(I1(JL),I2(JL),I3(JL))
      !
      ZCCS(JL) = PCCS(I1(JL),I2(JL),I3(JL))
      ZCRS(JL) = PCRS(I1(JL),I2(JL),I3(JL))
      ZCIS(JL) = PCIS(I1(JL),I2(JL),I3(JL))
      ZCSS(JL) = PCSS(I1(JL),I2(JL),I3(JL)) 
      ZCGS(JL) = PCGS(I1(JL),I2(JL),I3(JL)) 
      ZCHS(JL) = PCHS(I1(JL),I2(JL),I3(JL))       
      DO JMOD_IFN = 1, NMOD_IFN
         ZIFS(JL,JMOD_IFN) = PIFS(I1(JL),I2(JL),I3(JL),JMOD_IFN)
         ZINS(JL,JMOD_IFN) = PINS(I1(JL),I2(JL),I3(JL),JMOD_IFN)
      ENDDO
      !
      ZRHODREF(JL) = PRHODREF(I1(JL),I2(JL),I3(JL))
      ZZT(JL)      = ZT(I1(JL),I2(JL),I3(JL))
      ZPRES(JL)    = PPABST(I1(JL),I2(JL),I3(JL))
      ZEXNREF(JL)  = PEXNREF(I1(JL),I2(JL),I3(JL))
   ENDDO
   IF (NBUMOD==KMI .AND. LBU_ENABLE) THEN
      ALLOCATE(ZRHODJ(IMICRO))
      ZRHODJ(:) = PACK( PRHODJ(:,:,:),MASK=GMICRO(:,:,:) )
   END IF
!
! Atmospheric parameters 
!
   ALLOCATE(ZZW(IMICRO))
   ALLOCATE(ZLSFACT(IMICRO))
   ALLOCATE(ZLVFACT(IMICRO))
   ALLOCATE(ZSSI(IMICRO))
   ALLOCATE(ZAI(IMICRO))
   ALLOCATE(ZCJ(IMICRO))
   ALLOCATE(ZKA(IMICRO))
   ALLOCATE(ZDV(IMICRO))
!
   ZZW(:)  = ZEXNREF(:)*( XCPD+XCPV*ZRVT(:)+XCL*(ZRCT(:)+ZRRT(:)) &
        +XCI*(ZRIT(:)+ZRST(:)+ZRGT(:)) )
!
   ZLSFACT(:) = (XLSTT+(XCPV-XCI)*(ZZT(:)-XTT))/ZZW(:) ! L_s/(Pi_ref*C_ph)
   ZLVFACT(:) = (XLVTT+(XCPV-XCL)*(ZZT(:)-XTT))/ZZW(:) ! L_v/(Pi_ref*C_ph)
!
   ZZW(:) = EXP( XALPI - XBETAI/ZZT(:) - XGAMI*ALOG(ZZT(:) ) )
   ZSSI(:) = ZRVT(:)*( ZPRES(:)-ZZW(:) ) / ( (XMV/XMD) * ZZW(:) ) - 1.0
                                                       ! Supersaturation over ice
!
   ZKA(:) = 2.38E-2 + 0.0071E-2 * ( ZZT(:) - XTT )          ! k_a
   ZDV(:) = 0.211E-4 * (ZZT(:)/XTT)**1.94 * (XP00/ZPRES(:)) ! D_v
!
! Thermodynamical function ZAI = A_i(T,P)
   ZAI(:) = ( XLSTT + (XCPV-XCI)*(ZZT(:)-XTT) )**2 / (ZKA(:)*XRV*ZZT(:)**2) &
                                         + ( XRV*ZZT(:) ) / (ZDV(:)*ZZW(:))
! ZCJ = c^prime_j (in the ventilation factor)
   ZCJ(:) = XSCFAC * ZRHODREF(:)**0.3 / SQRT( 1.718E-5+0.0049E-5*(ZZT(:)-XTT) )
!
!
! Particle distribution parameters
!
   ALLOCATE(ZLBDAC(IMICRO)) 
   ALLOCATE(ZLBDAR(IMICRO))
   ALLOCATE(ZLBDAI(IMICRO)) 
   ALLOCATE(ZLBDAS(IMICRO))
   ALLOCATE(ZLBDAG(IMICRO))
   ALLOCATE(ZLBDAH(IMICRO))
   ZLBDAC(:)  = 1.E10
   WHERE (ZRCT(:)>XRTMIN(2) .AND. ZCCT(:)>XCTMIN(2))
      ZLBDAC(:) = ( XLBC*ZCCT(:) / ZRCT(:) )**XLBEXC
   END WHERE
   ZLBDAR(:)  = 1.E10
   WHERE (ZRRT(:)>XRTMIN(3) .AND. ZCRT(:)>XCTMIN(3))
      ZLBDAR(:) = ( XLBR*ZCRT(:) / ZRRT(:) )**XLBEXR
   END WHERE
   ZLBDAI(:)  = 1.E10
   WHERE (ZRIT(:)>XRTMIN(4) .AND. ZCIT(:)>XCTMIN(4))
      ZLBDAI(:) = ( XLBI*ZCIT(:) / ZRIT(:) )**XLBEXI
   END WHERE
   IF (LSNOW_T .AND. NMOM_S.EQ.1) THEN
      WHERE(ZZT(:)>263.15 .AND. ZRST(:)>XRTMIN(5)) 
         ZLBDAS(:) = MAX(MIN(XLBDAS_MAX, 10**(14.554-0.0423*ZZT(:))),XLBDAS_MIN)
      END WHERE
      WHERE(ZZT(:)<=263.15 .AND. ZRST(:)>XRTMIN(5)) 
         ZLBDAS(:) = MAX(MIN(XLBDAS_MAX, 10**(6.226-0.0106*ZZT(:))),XLBDAS_MIN)
      END WHERE
      ZLBDAS(:) = ZLBDAS(:) * XTRANS_MP_GAMMAS
      ZCST(:) = (XNS*ZRST(:)*ZLBDAS(:)**XBS)
      ZCSS(:) = (XNS*ZRSS(:)*ZLBDAS(:)**XBS)
   ELSE IF (NMOM_S.GE.2) THEN
	WHERE (ZRST(:)>XRTMIN(5) .AND. ZCST(:)>XCTMIN(5))
                ZLBDAS(:) = ( XLBS*ZCST(:) / ZRST(:) )**XLBEXS 
        END WHERE
   ELSE
      ZLBDAS(:)  = 1.E10
      WHERE (ZRSS(:)>XRTMIN(5)/PTSTEP)
         ZLBDAS(:) = MAX(MIN(XLBDAS_MAX,XLBS*( ZRHODREF(:)*ZRSS(:)*PTSTEP )**XLBEXS),XLBDAS_MIN)
      END WHERE
      ZCSS(:) = XCCS*ZLBDAS(:)**XCXS / ZRHODREF(:) / PTSTEP
      ZLBDAS(:)  = 1.E10
      WHERE (ZRST(:)>XRTMIN(5) )
         ZLBDAS(:) = MAX(MIN(XLBDAS_MAX,XLBS*( ZRHODREF(:)*ZRST(:) )**XLBEXS),XLBDAS_MIN)
      END WHERE
      ZCST(:) = XCCS*ZLBDAS(:)**XCXS / ZRHODREF(:)
   END IF
   ZLBDAG(:)  = 1.E10
   if (NMOM_G.GE.2) then
       WHERE (ZRGT(:)>XRTMIN(6) .AND. ZCGT(:)>XCTMIN(6) )  
           ZLBDAG(:) =  ( XLBG*ZCGT(:) / ZRGT(:) )**XLBEXG 
       END WHERE
   else
       WHERE (ZRGT(:)>XRTMIN(6) )
           ZLBDAG(:) = XLBG*( ZRHODREF(:)*ZRGT(:) )**XLBEXG
       END WHERE
      ZCGT(:) = XCCG*ZLBDAG(:)**XCXG / ZRHODREF(:)
      ZCGS(:) = XCCG*ZLBDAG(:)**XCXG / ZRHODREF(:)
   end if
   ZLBDAH(:)  = 1.E10
   if (NMOM_H.GE.2) then
       WHERE (ZRHT(:)>XRTMIN(7) .AND.  ZCHT(:)>XCTMIN(7))  
           ZLBDAH(:) =  ( XLBH*ZCHT(:) / ZRHT(:) )**XLBEXH 
       END WHERE
   else
       WHERE (ZRHT(:)>XRTMIN(7) )
           ZLBDAH(:) = XLBH*( ZRHODREF(:)*ZRHT(:) )**XLBEXH
       END WHERE
      ZCHT(:) = XCCH*ZLBDAH(:)**XCXH / ZRHODREF(:)
      ZCHS(:) = XCCH*ZLBDAH(:)**XCXH / ZRHODREF(:)
   end if
! 
!-------------------------------------------------------------------------------
!
!
!*       2.     Compute the slow processes involving cloud water and graupel
!	        ------------------------------------------------------------
!
   CALL LIMA_MIXED_SLOW_PROCESSES(ZRHODREF, ZZT, ZSSI, PTSTEP,        &
                                  ZLSFACT, ZLVFACT, ZAI, ZCJ,         &
                                  ZRGT, ZRHT, ZCIT, ZCGT, ZCHT,       &
                                  ZRVS, ZRCS, ZRIS, ZRGS, ZRHS, ZTHS, &
                                  ZCCS, ZCIS, ZCGS, ZIFS, ZINS,       &
                                  ZLBDAI, ZLBDAG, ZLBDAH,             &
                                  ZRHODJ, GMICRO, PRHODJ, KMI,        &
                                  PTHS, PRVS, PRCS, PRIS, PRGS, PRHS, &
                                  PCCS, PCIS                          )
! 
!-------------------------------------------------------------------------------
!
!
!        3.     Compute the fast RS and RG processes
!   	        ------------------------------------
!
IF (NMOM_S.GE.1) THEN
   CALL LIMA_MIXED_FAST_PROCESSES(ZRHODREF, ZZT, ZPRES, PTSTEP,                   &
                                  ZLSFACT, ZLVFACT, ZKA, ZDV, ZCJ,                &
                                  ZRVT, ZRCT, ZRRT, ZRIT, ZRST, ZRGT,             &
                                  ZRHT, ZCCT, ZCRT, ZCIT, ZCST, ZCGT, ZCHT,       &
                                  ZRCS, ZRRS, ZRIS, ZRSS, ZRGS, ZRHS,             &
                                  ZTHS, ZCCS, ZCRS, ZCIS, ZCSS, ZCGS, ZCHS,       &
                                  ZLBDAC, ZLBDAR, ZLBDAI, ZLBDAS, ZLBDAG, ZLBDAH, &
                                  ZRHODJ, GMICRO, PRHODJ, KMI, PTHS,              &
                                  PRCS, PRRS, PRIS, PRSS, PRGS, PRHS,             &
                                  PCCS, PCRS, PCIS, PCSS, PCGS, PCHS              )
END IF
!
!-------------------------------------------------------------------------------
!
!
!
!        4.     Unpack variables
!   	        ----------------
!
!
   ZW(:,:,:) = PRVS(:,:,:)
   PRVS(:,:,:) = UNPACK( ZRVS(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
   ZW(:,:,:) = PRCS(:,:,:)
   PRCS(:,:,:) = UNPACK( ZRCS(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
   ZW(:,:,:) = PRRS(:,:,:)
   PRRS(:,:,:) = UNPACK( ZRRS(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
   ZW(:,:,:) = PRIS(:,:,:)
   PRIS(:,:,:) = UNPACK( ZRIS(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
   ZW(:,:,:) = PRSS(:,:,:)
   PRSS(:,:,:) = UNPACK( ZRSS(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
   ZW(:,:,:) = PRGS(:,:,:)
   PRGS(:,:,:) = UNPACK( ZRGS(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
   ZW(:,:,:) = PRHS(:,:,:)
   PRHS(:,:,:) = UNPACK( ZRHS(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
!
   ZW(:,:,:) = PTHS(:,:,:)
   PTHS(:,:,:) = UNPACK( ZTHS(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
!
   ZW(:,:,:) = PCCS(:,:,:)
   PCCS(:,:,:) = UNPACK( ZCCS(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
   ZW(:,:,:) = PCRS(:,:,:)
   PCRS(:,:,:) = UNPACK( ZCRS(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
   ZW(:,:,:) = PCIS(:,:,:)
   PCIS(:,:,:) = UNPACK( ZCIS(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
   if (NMOM_S.GE.2) then
       ZW(:,:,:) = PCSS(:,:,:)
       PCSS(:,:,:) = UNPACK( ZCSS(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) ) 
   end if
   if (NMOM_G.GE.2) then
       ZW(:,:,:) = PCGS(:,:,:)
       PCGS(:,:,:) = UNPACK( ZCGS(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
   end if
   if (NMOM_H.GE.2) then
       ZW(:,:,:) = PCHS(:,:,:)      
       PCHS(:,:,:) = UNPACK( ZCHS(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) ) 
   end if
!
   DO JMOD_IFN = 1, NMOD_IFN
      ZW(:,:,:) = PIFS(:,:,:,JMOD_IFN)
      PIFS(:,:,:,JMOD_IFN) = UNPACK( ZIFS(:,JMOD_IFN),MASK=GMICRO(:,:,:),         &
                                                    FIELD=ZW(:,:,:) )
      ZW(:,:,:) = PINS(:,:,:,JMOD_IFN)
      PINS(:,:,:,JMOD_IFN) = UNPACK( ZINS(:,JMOD_IFN),MASK=GMICRO(:,:,:),         &
                                                    FIELD=ZW(:,:,:) )
   ENDDO
!
   DEALLOCATE(ZRVT) 
   DEALLOCATE(ZRCT)
   DEALLOCATE(ZRRT)  
   DEALLOCATE(ZRIT) 
   DEALLOCATE(ZRST) 
   DEALLOCATE(ZRGT)
   DEALLOCATE(ZRHT)
!  
   DEALLOCATE(ZCCT) 
   DEALLOCATE(ZCRT) 
   DEALLOCATE(ZCIT)
   DEALLOCATE(ZCST) 
   DEALLOCATE(ZCGT)
   DEALLOCATE(ZCHT) 
! 
   DEALLOCATE(ZRVS)  
   DEALLOCATE(ZRCS) 
   DEALLOCATE(ZRRS) 
   DEALLOCATE(ZRIS)
   DEALLOCATE(ZRSS)
   DEALLOCATE(ZRGS) 
   DEALLOCATE(ZRHS) 
   DEALLOCATE(ZTHS)
!
   DEALLOCATE(ZCCS) 
   DEALLOCATE(ZCRS) 
   DEALLOCATE(ZCIS)  
   DEALLOCATE(ZCSS)  
   DEALLOCATE(ZCGS)     
   DEALLOCATE(ZCHS)
   DEALLOCATE(ZIFS)
   DEALLOCATE(ZINS)
!
   DEALLOCATE(ZRHODREF) 
   DEALLOCATE(ZZT) 
   DEALLOCATE(ZPRES) 
   DEALLOCATE(ZEXNREF)
!
   DEALLOCATE(ZZW)
   DEALLOCATE(ZLSFACT)
   DEALLOCATE(ZLVFACT)
   DEALLOCATE(ZSSI)
   DEALLOCATE(ZAI)
   DEALLOCATE(ZCJ)
   DEALLOCATE(ZKA)
   DEALLOCATE(ZDV)
!
   DEALLOCATE(ZLBDAC) 
   DEALLOCATE(ZLBDAR)
   DEALLOCATE(ZLBDAI) 
   DEALLOCATE(ZLBDAS)
   DEALLOCATE(ZLBDAG)
   DEALLOCATE(ZLBDAH)
!
   IF (NBUMOD==KMI .AND. LBU_ENABLE) DEALLOCATE(ZRHODJ)
!
END IF ! IMICRO >= 1
!
!------------------------------------------------------------------------------
!
!
!*       5.    REPORT 3D MICROPHYSICAL VARIABLES IN PRS AND PSVS
!              -------------------------------------------------
!
PRS(:,:,:,1) = PRVS(:,:,:)
IF ( KRR .GE. 2 ) PRS(:,:,:,2) = PRCS(:,:,:)
IF ( KRR .GE. 3 ) PRS(:,:,:,3) = PRRS(:,:,:)
IF ( KRR .GE. 4 ) PRS(:,:,:,4) = PRIS(:,:,:)
IF ( KRR .GE. 5 ) PRS(:,:,:,5) = PRSS(:,:,:)
IF ( KRR .GE. 6 ) PRS(:,:,:,6) = PRGS(:,:,:)
IF ( KRR .GE. 7 ) PRS(:,:,:,7) = PRHS(:,:,:)
!
! Prepare 3D number concentrations
!
IF ( NMOM_C.GE.2 ) PSVS(:,:,:,NSV_LIMA_NC) = PCCS(:,:,:)
IF ( NMOM_R.GE.2 ) PSVS(:,:,:,NSV_LIMA_NR) = PCRS(:,:,:)
IF ( NMOM_I.GE.2 ) PSVS(:,:,:,NSV_LIMA_NI) = PCIS(:,:,:)
IF ( NMOM_S.GE.2 ) PSVS(:,:,:,NSV_LIMA_NS) = PCSS(:,:,:)
IF ( NMOM_G.GE.2 ) PSVS(:,:,:,NSV_LIMA_NG) = PCGS(:,:,:)
IF ( NMOM_H.GE.2 ) PSVS(:,:,:,NSV_LIMA_NH) = PCHS(:,:,:)
!
IF ( NMOD_CCN .GE. 1 ) THEN
   PSVS(:,:,:,NSV_LIMA_CCN_FREE:NSV_LIMA_CCN_FREE+NMOD_CCN-1) = PNFS(:,:,:,:)
   PSVS(:,:,:,NSV_LIMA_CCN_ACTI:NSV_LIMA_CCN_ACTI+NMOD_CCN-1) = PNAS(:,:,:,:)
END IF
!
IF ( NMOD_IFN .GE. 1 ) THEN
   PSVS(:,:,:,NSV_LIMA_IFN_FREE:NSV_LIMA_IFN_FREE+NMOD_IFN-1) = PIFS(:,:,:,:)
   PSVS(:,:,:,NSV_LIMA_IFN_NUCL:NSV_LIMA_IFN_NUCL+NMOD_IFN-1) = PINS(:,:,:,:)
END IF
!
IF ( NMOD_IMM .GE. 1 ) THEN
   PSVS(:,:,:,NSV_LIMA_IMM_NUCL:NSV_LIMA_IMM_NUCL+NMOD_IMM-1) = PNIS(:,:,:,:)
END IF
!
!++cb++
IF (ALLOCATED(PNFS)) DEALLOCATE(PNFS)
IF (ALLOCATED(PNAS)) DEALLOCATE(PNAS)
IF (ALLOCATED(PIFS)) DEALLOCATE(PIFS)
IF (ALLOCATED(PINS)) DEALLOCATE(PINS)
IF (ALLOCATED(PNIS)) DEALLOCATE(PNIS)
IF (ALLOCATED(PNHS)) DEALLOCATE(PNHS)
!--cb--
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_MIXED
