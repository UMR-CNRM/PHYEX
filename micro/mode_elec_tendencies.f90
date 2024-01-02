!MNH_LIC Copyright 2022-2023 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_ELEC_TENDENCIES
!
IMPLICIT NONE
CONTAINS
!
!     #########################################################################################
      SUBROUTINE ELEC_TENDENCIES (D, CST, ICED, ICEP, ELECD, ELECP,                           &
                                  KRR, KMICRO, PTSTEP, ODMICRO,                               &
                                  BUCONF, TBUDGETS, KBUDGETS,                                 &
                                  HCLOUD, PTHVREFZIKB,                                        &
                                  PRHODREF, PRHODJ, PZT, PCIT,                                &
                                  PRVT, PRCT, PRRT, PRIT, PRST, PRGT,                         &
                                  PQPIT, PQCT, PQRT, PQIT, PQST, PQGT, PQNIT,                 &
                                  PQPIS, PQCS, PQRS, PQIS, PQSS, PQGS, PQNIS,                 &
                                  PRVHENI, PRRHONG, PRIMLTC,                                  &
                                  PRCHONI, PRVDEPS, PRIAGGS, PRIAUTS, PRVDEPG,                &
                                  PRCAUTR, PRCACCR, PRREVAV,                                  &
                                  PRCRIMSS, PRCRIMSG, PRSRIMCG, PRRACCSS, PRRACCSG, PRSACCRG, &
                                  PRSMLTG, PRICFRRG, PRRCFRIG,                                &
                                  PRCWETG, PRIWETG, PRRWETG, PRSWETG,                         &
                                  PRCDRYG, PRIDRYG, PRRDRYG, PRSDRYG,                         &
                                  PRGMLTR, PRCBERI,                                           & 
                                  PRCMLTSR, PRICFRR,                                          & !- opt. param. for ICE3 
                                  PCCT, PCRT, PCST, PCGT,                                     & !-- optional
                                  PRVHENC, PRCHINC, PRVHONH,                                  & !| parameters 
                                  PRRCVRC, PRICNVI, PRVDEPI, PRSHMSI, PRGHMGI,                & !|    for
                                  PRICIBU, PRIRDSF,                                           & !|    LIMA
                                  PRCCORR2, PRRCORR2, PRICORR2,                               & !--
                                  PRWETGH, PRCWETH, PRIWETH, PRSWETH, PRGWETH, PRRWETH,       & !--  optional
                                  PRCDRYH, PRIDRYH, PRSDRYH, PRRDRYH, PRGDRYH,                & !|  parameters
                                  PRHMLTR, PRDRYHG,                                           & !|     for
                                  PRHT, PRHS, PCHT, PQHT, PQHS)                                 !--    hail
!     ##########################################################################################
!
!!****  * - compute the explicit cloud electrification sources
!!
!!    This routine is adapted from rain_ice_elec.f90.
!!    To avoid duplicated routines, the cloud electrification routine is now CALLed 
!!    at the end of the microphysics scheme but needs the microphysical tendencies as arguments.
!!    The sedimentation source for electric charges is treated separately.
!!
!!    AUTHOR
!!    ------
!!      C. Barthe    * LAERO *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    February 2022
!!
!!      Modifications
!! C. Barthe   12/04/2022   include electrification from LIMA
!! C. Barthe   22/03/2023   5-6: take into account news from LIMA (Ns, Ng, Nh, CIBU and RDSF) and PHYEX
!! C. Barthe   13/07/2023   5-6: Ns, Ng and Nh can be pronostic variables (LIMA2)
!! C. Barthe   22/11/2023   initialize Nx to 0 when 1-moment
!!
!------------------------------------------------------------------
!
!*      0.    DECLARATIONS
!             ------------
!
USE MODD_BUDGET,     ONLY:  NBUDGET_TH, NBUDGET_RV, NBUDGET_RC, NBUDGET_RR, NBUDGET_RI, &
                            NBUDGET_RS, NBUDGET_RG, NBUDGET_RH, NBUDGET_SV1,            &
                            TBUDGETDATA, TBUDGETCONF_t
!
USE MODD_CST,              ONLY: CST_t
USE MODD_DIMPHYEX,         ONLY: DIMPHYEX_t
USE MODD_RAIN_ICE_DESCR_n, ONLY: RAIN_ICE_DESCR_t
USE MODD_RAIN_ICE_PARAM_n, ONLY: RAIN_ICE_PARAM_t
USE MODD_NSV, ONLY: NSV_ELECBEG, NSV_ELECEND
USE MODD_ELEC_DESCR      
USE MODD_ELEC_PARAM      
USE MODD_ELEC_n
USE MODD_PARAM_LIMA,       ONLY: XALPHAI_L=>XALPHAI, XNUI_L=>XNUI,   &
                                 XCEXVT_L=>XCEXVT, XRTMIN_L=>XRTMIN, &
                                 LCIBU, LRDSF,                       &
                                 NMOM_C, NMOM_R, NMOM_I, NMOM_S, NMOM_G, NMOM_H
USE MODD_PARAM_LIMA_COLD,  ONLY: XAI_L=>XAI, XBI_L=>XBI,   &
                                 XDS_L=>XDS, XCXS_L=>XCXS, &
                                 XCOLEXIS_L=>XCOLEXIS
USE MODD_PARAM_LIMA_MIXED, ONLY: XDG_L=>XDG, XCXG_L=>XCXG,                           &
                                 XCOLIG_L=>XCOLIG, XCOLEXIG_L=>XCOLEXIG,             &
                                 XCOLSG_L=>XCOLSG, XCOLEXSG_L=>XCOLEXSG,             &
                                 NGAMINC_L=>NGAMINC,                                 &
                                 NACCLBDAR_L=>NACCLBDAR, NACCLBDAS_L=>NACCLBDAS,     &
                                 XACCINTP1S_L=>XACCINTP1S, XACCINTP2S_L=>XACCINTP2S, &
                                 XACCINTP1R_L=>XACCINTP1R, XACCINTP2R_L=>XACCINTP2R, &
                                 NDRYLBDAR_L=>NDRYLBDAR, NDRYLBDAS_L=>NDRYLBDAS,     &
                                 NDRYLBDAG_L=>NDRYLBDAG,                             &
                                 XDRYINTP1R_L=>XDRYINTP1R, XDRYINTP2R_L=>XDRYINTP2R, &
                                 XDRYINTP1S_L=>XDRYINTP1S, XDRYINTP2S_L=>XDRYINTP2S, &
                                 XDRYINTP1G_L=>XDRYINTP1G, XDRYINTP2G_L=>XDRYINTP2G, &
                                 XRIMINTP1_L=>XRIMINTP1,   XRIMINTP2_L=>XRIMINTP2
!
!#ifdef MNH_PGI
!USE MODE_PACK_PGI
!#endif
use mode_tools,           only: Countjv
USE MODE_BUDGET_PHY,      ONLY: BUDGET_STORE_ADD_PHY, BUDGET_STORE_INIT_PHY, BUDGET_STORE_END_PHY
!
USE MODE_COMPUTE_LAMBDA, ONLY: COMPUTE_LAMBDA
USE MODE_ELEC_COMPUTE_EX,ONLY: ELEC_COMPUTE_EX
USE MODI_MOMG
!
IMPLICIT NONE
!
!
!*      0.1   Declaration of dummy arguments
!
TYPE(DIMPHYEX_t),                      INTENT(IN)     :: D
TYPE(CST_t),                           INTENT(IN)     :: CST
TYPE(TBUDGETCONF_t),                   INTENT(IN)     :: BUCONF        ! budget structure
TYPE(RAIN_ICE_DESCR_t),                INTENT(IN)     :: ICED
TYPE(RAIN_ICE_PARAM_t),                INTENT(IN)     :: ICEP
TYPE(ELEC_PARAM_t),                    INTENT(IN)     :: ELECP   ! electrical parameters
TYPE(ELEC_DESCR_t),                    INTENT(IN)     :: ELECD   ! electrical descriptive csts
TYPE(TBUDGETDATA), DIMENSION(KBUDGETS),INTENT(INOUT)  :: TBUDGETS
INTEGER,                               INTENT(IN)     :: KBUDGETS
!
INTEGER,                              INTENT(IN)    :: KMICRO
REAL,                                 INTENT(IN)    :: PTSTEP  ! Double Time step
                                                               ! (single if cold start)
INTEGER,                              INTENT(IN)    :: KRR     ! Number of moist variable
!
LOGICAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   :: ODMICRO ! mask to limit computation
CHARACTER (LEN=4),      INTENT(IN)   ::  HCLOUD       ! Kind of microphysical scheme
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    :: PZT     ! Temperature (K)
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    :: PCIT    ! Pristine ice n.c. at t
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    :: PRVT    ! Water vapor m.r. at t
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    :: PRCT    ! Cloud water m.r. at t
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    :: PRRT    ! Rain water m.r. at t
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    :: PRIT    ! Pristine ice m.r. at t
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    :: PRST    ! Snow/aggregate m.r. at t
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    :: PRGT    ! Graupel/hail m.r. at t
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    :: PQPIT   ! Positive ion (Nb/kg) at t
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    :: PQNIT   ! Negative ion (Nb/kg) at t
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    :: PQCT    ! Cloud water charge at t
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    :: PQRT    ! Raindrops charge at t
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    :: PQIT    ! Pristine ice charge at t
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    :: PQST    ! Snow/aggregates charge at t
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    :: PQGT    ! Graupel charge at t
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(INOUT) :: PQPIS   ! Positive ion (Nb/kg) source
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(INOUT) :: PQNIS   ! Negative ion (Nb/kg) source
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(INOUT) :: PQCS    ! Cloud water charge source
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(INOUT) :: PQRS    ! Raindrops charge source
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(INOUT) :: PQIS    ! Pristine ice charge source
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(INOUT) :: PQSS    ! Snow/aggregates charge source
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(INOUT) :: PQGS    ! Graupel charge source
!
! microphysics rates common to ICE3 and LIMA
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    :: PRVHENI, &  ! heterogeneous nucleation mixing ratio change (HIND for LIMA)
                                                       PRCHONI, &  ! Homogeneous nucleation
                                                       PRRHONG, &  ! Spontaneous freezing mixing ratio change
                                                       PRVDEPS, &  ! Deposition on r_s,
                                                       PRIAGGS, &  ! Aggregation on r_s
                                                       PRIAUTS, &  ! Autoconversion of r_i for r_s production (CNVS for LIMA)
                                                       PRVDEPG, &  ! Deposition on r_g
                                                       PRCAUTR, &  ! Autoconversion of r_c for r_r production
                                                       PRCACCR, &  ! Accretion of r_c for r_r production
                                                       PRREVAV, &  ! Evaporation of r_r
                                                       PRIMLTC, &  ! Cloud ice melting mixing ratio change
                                                       PRCBERI, &  ! Bergeron-Findeisen effect
                                                       PRSMLTG, &  ! Conversion-Melting of the aggregates
                                                       PRRACCSS, PRRACCSG, PRSACCRG, & ! Rain accretion onto the aggregates
                                                       PRCRIMSS, PRCRIMSG, PRSRIMCG, & ! Cloud droplet riming of the aggregates
                                                       PRICFRRG, PRRCFRIG,           & ! Rain contact freezing
                                                       PRCWETG, PRIWETG, PRRWETG, PRSWETG, &  ! Graupel wet growth
                                                       PRCDRYG, PRIDRYG, PRRDRYG, PRSDRYG, &  ! Graupel dry growth
                                                       PRGMLTR     ! Melting of the graupel
! microphysics rates specific to ICE3 (knmoments==1)
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), OPTIONAL, INTENT(IN)    :: PRCMLTSR,&  ! Cld droplet collection onto aggregates by pos. temp.
                                                               PRICFRR     ! Rain contact freezing (part of ice crystals converted to rain)
! microphysics rates specific to LIMA (knmoments==2)
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), OPTIONAL, INTENT(IN)    :: PRVHENC, &  ! Cld droplet formation
                                                               PRCHINC, &  ! Heterogeneous nucleation of coated IFN
                                                               PRVHONH, &  ! Nucleation of haze
                                                               PRRCVRC, &  ! Conversion of small drops into droplets
                                                               PRICNVI, &  ! Conversion snow --> ice
                                                               PRVDEPI, &  ! Deposition on r_i
                                                               PRSHMSI, PRGHMGI, & ! Hallett Mossop for snow and graupel
                                                               PRICIBU, &  ! Collisional ice breakup
                                                               PRIRDSF, &  ! Raindrop shattering by freezing
                                                               PRCCORR2, PRRCORR2, PRICORR2 ! Correction inside LIMA splitting
! microphysics rates related to hail (krr == 7, lhail = .t.)
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), OPTIONAL, INTENT(IN)    :: PRWETGH, &  ! Conversion of graupel into hail
                                                               PRCWETH, PRIWETH, PRSWETH, PRGWETH, PRRWETH, & ! Dry growth of hail
                                                               PRCDRYH, PRIDRYH, PRSDRYH, PRRDRYH, PRGDRYH, & ! Wet growth of hail
                                                               PRHMLTR, &                                     ! Melting of hail  
                                                               PRDRYHG     ! Conversion of hail into graupel
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), OPTIONAL, INTENT(IN)    :: PCCT   ! Cloud droplets conc. at t
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), OPTIONAL, INTENT(IN)    :: PCRT   ! Raindrops conc. at t
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), OPTIONAL, INTENT(IN)    :: PCST   ! Snow conc. at t
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), OPTIONAL, INTENT(IN)    :: PCGT   ! Graupel conc. at t
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), OPTIONAL, INTENT(IN)    :: PCHT   ! Hail conc. at t
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), OPTIONAL, INTENT(IN)    :: PRHT   ! Hail m.r. at t
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), OPTIONAL, INTENT(INOUT) :: PRHS   ! Hail m.r. source
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), OPTIONAL, INTENT(IN)    :: PQHT   ! Hail charge at t
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), OPTIONAL, INTENT(INOUT) :: PQHS   ! Hail charge source
REAL, INTENT(IN)                :: PTHVREFZIKB ! Reference thv at IKB for electricity
!
!
!*      0.2   Declaration of local variables    
!
INTEGER :: II, JJ, JL     ! Loop indexes
INTEGER :: IIB, IIE, &    ! Define the domain
           IJB, IJE, &    ! where the microphysical sources
           IKB, IKE       ! must be computed
INTEGER                    :: IMICRO   ! nb of pts where r_x > 0
INTEGER, DIMENSION(KMICRO) :: I1
INTEGER, DIMENSION(KMICRO) :: II1, II2, II3
!
LOGICAL, DIMENSION(KMICRO) :: GMASK    !          Mask 
!REAL,    DIMENSION(KMICRO) :: ZMASK    !       to reduce
INTEGER                    :: IGMASK   ! the computation domain
!
REAL, DIMENSION(KMICRO) :: ZRHODREF  ! Reference density
REAL, DIMENSION(KMICRO) :: ZRHODJ    ! RHO times Jacobian
REAL, DIMENSION(KMICRO) :: ZZT     ! Temperature
!
REAL, DIMENSION(KMICRO) :: ZRVT    ! Water vapor m.r. at t 
REAL, DIMENSION(KMICRO) :: ZRCT    ! Cloud water m.r. at t 
REAL, DIMENSION(KMICRO) :: ZRRT    ! Rain water m.r. at t 
REAL, DIMENSION(KMICRO) :: ZRIT    ! Pristine ice m.r. at t
REAL, DIMENSION(KMICRO) :: ZRST    ! Snow/aggregate m.r. at t
REAL, DIMENSION(KMICRO) :: ZRGT    ! Graupel m.r. at t
REAL, DIMENSION(KMICRO) :: ZRHT    ! Hail m.r. at t
REAL, DIMENSION(KMICRO) :: ZCCT    ! Cloud water conc. at t
REAL, DIMENSION(KMICRO) :: ZCRT    ! Raindrops conc. at t
REAL, DIMENSION(KMICRO) :: ZCIT    ! Pristine ice conc. at t
REAL, DIMENSION(KMICRO) :: ZCST    ! Snow/aggregate conc. at t
REAL, DIMENSION(KMICRO) :: ZCGT    ! Graupel conc. at t
REAL, DIMENSION(KMICRO) :: ZCHT    ! Hail conc. at t
!
REAL, DIMENSION(KMICRO) :: ZQPIT   ! Positive ion (/kg) at t
REAL, DIMENSION(KMICRO) :: ZQNIT   ! Negative ion (/kg) at t
REAL, DIMENSION(KMICRO) :: ZQCT    ! Cloud water charge at t
REAL, DIMENSION(KMICRO) :: ZQRT    ! Raindrops charge at t
REAL, DIMENSION(KMICRO) :: ZQIT    ! Pristine ice charge at t
REAL, DIMENSION(KMICRO) :: ZQST    ! Snow/aggregate charge at t
REAL, DIMENSION(KMICRO) :: ZQGT    ! Graupel charge at t
REAL, DIMENSION(KMICRO) :: ZQHT    ! Hail charge at t
!
REAL, DIMENSION(KMICRO) :: ZQPIS   ! Positive ion (/kg) source
REAL, DIMENSION(KMICRO) :: ZQNIS   ! Negative ion (/kg) source
REAL, DIMENSION(KMICRO) :: ZQCS    ! Cloud water charge source
REAL, DIMENSION(KMICRO) :: ZQRS    ! Raindrops charge source
REAL, DIMENSION(KMICRO) :: ZQIS    ! Pristine ice charge source
REAL, DIMENSION(KMICRO) :: ZQSS    ! Snow/aggregate charge source
REAL, DIMENSION(KMICRO) :: ZQGS    ! Graupel charge source
REAL, DIMENSION(KMICRO) :: ZQHS    ! Hail charge source
!
REAL, DIMENSION(KMICRO) :: ZLBDAC  ! Slope parameter of the droplets distribution
REAL, DIMENSION(KMICRO) :: ZLBDAR  ! Slope parameter of the raindrop distribution
REAL, DIMENSION(KMICRO) :: ZLBDAI  ! Slope parameter of the pristine ice distribution
REAL, DIMENSION(KMICRO) :: ZLBDAS  ! Slope parameter of the aggregate distribution
REAL, DIMENSION(KMICRO) :: ZLBDAG  ! Slope parameter of the graupel distribution
REAL, DIMENSION(KMICRO) :: ZLBDAH  ! Slope parameter of the hail distribution
!
REAL, DIMENSION(KMICRO) :: ZECT    !
REAL, DIMENSION(KMICRO) :: ZERT    !     e_x coef
REAL, DIMENSION(KMICRO) :: ZEIT    !      in the 
REAL, DIMENSION(KMICRO) :: ZEST    ! q_x - D_x relation
REAL, DIMENSION(KMICRO) :: ZEGT    !
REAL, DIMENSION(KMICRO) :: ZEHT    !
!
LOGICAL, DIMENSION(KMICRO,4) :: GELEC  ! Mask for non-inductive charging
!
REAL, DIMENSION(:), ALLOCATABLE :: ZDQ, ZDQ_IS, ZDQ_IG, ZDQ_SG
!
! Non-inductive charging process following Gardiner et al. (1995)
REAL, DIMENSION(:), ALLOCATABLE :: ZDELTALWC ! Gap between LWC and a critical LWC
REAL, DIMENSION(:), ALLOCATABLE :: ZFT       ! Fct depending on temperature
!
! Non-inductive charging process following Saunders et al. (1991) / EW
REAL, DIMENSION(:), ALLOCATABLE :: ZEW      ! Effective liquid water content
REAL, DIMENSION(:), ALLOCATABLE :: ZSAUNSK  ! constant B 
REAL, DIMENSION(:), ALLOCATABLE :: ZSAUNIM  !  d_i exponent
REAL, DIMENSION(:), ALLOCATABLE :: ZSAUNIN  !  v_g/s-v_i
REAL, DIMENSION(:), ALLOCATABLE :: ZSAUNSM  !  d_s exponent
REAL, DIMENSION(:), ALLOCATABLE :: ZSAUNSN  !  v_g-v_s
REAL, DIMENSION(:), ALLOCATABLE :: ZFQIAGGS, ZFQIDRYGBS
REAL, DIMENSION(:), ALLOCATABLE :: ZLBQSDRYGB1S, ZLBQSDRYGB2S, ZLBQSDRYGB3S
!
! Non-inductive charging process following Saunders and Peck (1998) / RAR
REAL, DIMENSION(:), ALLOCATABLE :: ZRAR        ! Rime accretion rate
REAL, DIMENSION(:), ALLOCATABLE :: ZSAUNIM_IS  !  d_i exponent
REAL, DIMENSION(:), ALLOCATABLE :: ZSAUNIN_IS  !  v_g/s-v_i
REAL, DIMENSION(:), ALLOCATABLE :: ZSAUNIM_IG  !  d_i exponent
REAL, DIMENSION(:), ALLOCATABLE :: ZSAUNIN_IG  !  v_g/s-v_i
REAL, DIMENSION(:), ALLOCATABLE :: ZSAUNSK_SG  ! constant B 
REAL, DIMENSION(:), ALLOCATABLE :: ZSAUNSM_SG  !  d_s exponent
REAL, DIMENSION(:), ALLOCATABLE :: ZSAUNSN_SG  !  v_g-v_s
!
! Inductive charging process (Ziegler et al., 1991)
REAL, DIMENSION(:), ALLOCATABLE :: ZEFIELDW  ! Vertical component of the electric field
!
REAL, DIMENSION(KMICRO) :: ZLIMIT   ! Used to limit the charge separated during NI process 
REAL, DIMENSION(KMICRO) :: ZQCOLIS  ! Collection efficiency between ice and snow
REAL, DIMENSION(KMICRO) :: ZQCOLIG  ! Collection efficiency between ice and graupeln
REAL, DIMENSION(KMICRO) :: ZQCOLSG  ! Collection efficiency between snow and graupeln
!
REAL                    :: ZRHO00, ZCOR00   ! Surface reference air density
REAL, DIMENSION(KMICRO) :: ZRHOCOR  ! Density correction for fallspeed
!
INTEGER, DIMENSION(:), ALLOCATABLE :: IVEC1, IVEC2                   ! Vectors of indices for interpolation
REAL,    DIMENSION(:), ALLOCATABLE :: ZVEC1, ZVEC2, ZVEC3            ! Work vectors for interpolation
REAL,    DIMENSION(:), ALLOCATABLE :: ZVECQ1, ZVECQ2, ZVECQ3, ZVECQ4 ! Work vectors for interpolation
!
REAL,    DIMENSION(KMICRO)   :: ZWQ, ZWQ_NI                    !  Work arrays
REAL,    DIMENSION(KMICRO)   :: ZWQ1, ZWQ2, ZWQ3, ZWQ4         !      for 
REAL,    DIMENSION(KMICRO,9) :: ZWQ5                           ! charge transfer
!
! variables used to select between common parameters between ICEx and LIMA
INTEGER :: IMOM_C, IMOM_R, IMOM_I, IMOM_S, IMOM_G, IMOM_H ! number of moments for each hydrometeor
INTEGER :: IGAMINC,              &
           IACCLBDAR, IACCLBDAS, &
           IDRYLBDAR, IDRYLBDAS, IDRYLBDAG
!
REAL    :: ZCEXVT,                                               &
           ZALPHAI, ZNUI, ZAI, ZBI, ZDS, ZDG, ZCXS, ZCXG,        &
           ZCOLIS, ZCOLEXIS, ZCOLIG, ZCOLEXIG, ZCOLSG, ZCOLEXSG, &
           ZACCINTP1S, ZACCINTP2S, ZACCINTP1R, ZACCINTP2R,       &
           ZDRYINTP1R, ZDRYINTP2R, ZDRYINTP1S, ZDRYINTP2S,       &
           ZDRYINTP1G, ZDRYINTP2G,                               &
           ZRIMINTP1,  ZRIMINTP2
REAL, DIMENSION(:), ALLOCATABLE :: ZRTMIN
!
! microphysical tendencies have to be transformed in 1D arrays
REAL, DIMENSION(KMICRO) :: ZRVHENI, ZRCHONI, ZRRHONG, ZRVDEPS, ZRIAGGS,      &
                           ZRIAUTS, ZRVDEPG, ZRCAUTR, ZRCACCR, ZRREVAV,      &
                           ZRIMLTC, ZRCBERI, ZRSMLTG, ZRRACCSS, ZRRACCSG,    &
                           ZRSACCRG, ZRCRIMSS, ZRCRIMSG, ZRSRIMCG, ZRICFRRG, &
                           ZRRCFRIG, ZRCWETG, ZRIWETG, ZRRWETG, ZRSWETG,     &
                           ZRCDRYG, ZRIDRYG, ZRRDRYG, ZRSDRYG, ZRGMLTR
! optional microphysical tendencies
REAL, DIMENSION(:), ALLOCATABLE :: ZRCMLTSR, ZRICFRR, ZRVHENC, ZRCHINC, ZRVHONH,     &
                                   ZRRCVRC, ZRICNVI, ZRVDEPI, ZRSHMSI, ZRGHMGI,      &
                                   ZRICIBU, ZRIRDSF, ZRCCORR2, ZRRCORR2, ZRICORR2,   &
                                   ZRWETGH, ZRCWETH, ZRIWETH, ZRSWETH, ZRGWETH,      &
                                   ZRRWETH, ZRCDRYH, ZRIDRYH, ZRSDRYH, ZRRDRYH,      &
                                   ZRGDRYH, ZRHMLTR, ZRDRYHG 
!
!------------------------------------------------------------------
ASSOCIATE(XCEXVT_I=>ICED%XCEXVT, XRTMIN_I=>ICED%XRTMIN,                       &
                                XALPHAI_I=>ICED%XALPHAI, XNUI_I=>ICED%XNUI, XAI_I=>ICED%XAI, XBI_I=>ICED%XBI, &
                                XDS_I=>ICED%XDS, XDG_I=>ICED%XDG,                             &
                                XCXS_I=>ICED%XCXS, XCXG_I=>ICED%XCXG,                         &
                                XCOLIS_I=>ICEP%XCOLIS, XCOLEXIS_I=>ICEP%XCOLEXIS,             &
                                XCOLIG_I=>ICEP%XCOLIG, XCOLEXIG_I=>ICEP%XCOLEXIG,             &
                                XCOLSG_I=>ICEP%XCOLSG, XCOLEXSG_I=>ICEP%XCOLEXSG,             &
                                NGAMINC_I=>ICEP%NGAMINC,                                 &                                
                                NACCLBDAR_I=>ICEP%NACCLBDAR, NACCLBDAS_I=>ICEP%NACCLBDAS,     &
                                XACCINTP1S_I=>ICEP%XACCINTP1S, XACCINTP2S_I=>ICEP%XACCINTP2S, &
                                XACCINTP1R_I=>ICEP%XACCINTP1R, XACCINTP2R_I=>ICEP%XACCINTP2R, &
                                NDRYLBDAR_I=>ICEP%NDRYLBDAR, NDRYLBDAS_I=>ICEP%NDRYLBDAS,     &
                                NDRYLBDAG_I=>ICEP%NDRYLBDAG,                             &
                                XDRYINTP1R_I=>ICEP%XDRYINTP1R, XDRYINTP2R_I=>ICEP%XDRYINTP2R, &
                                XDRYINTP1S_I=>ICEP%XDRYINTP1S, XDRYINTP2S_I=>ICEP%XDRYINTP2S, &
                                XDRYINTP1G_I=>ICEP%XDRYINTP1G, XDRYINTP2G_I=>ICEP%XDRYINTP2G, &
                                XRIMINTP1_I=>ICEP%XRIMINTP1,   XRIMINTP2_I=>ICEP%XRIMINTP2    )
!
!*      1.    INITIALIZATIONS
!             ---------------
!
!*      1.1   compute the loop bounds
!
IIB = D%NIB
IIE = D%NIE
IJB = D%NJB
IJE = D%NJE
IKB = D%NKB
IKE = D%NKE
!
!
!*      1.2   select parameters between ICEx and LIMA
!
IF (HCLOUD(1:3) == 'ICE') THEN
  ZCEXVT = XCEXVT_I
  IMOM_C = 1
  IMOM_R = 1
  IMOM_I = 2 ! Ni is diagnostic and always available
  IMOM_S = 1
  IMOM_G = 1
  IF (KRR == 7) THEN
    IMOM_H = 1
  ELSE
    IMOM_H = 0
  END IF
ELSE IF (HCLOUD == 'LIMA') THEN
  ZCEXVT = XCEXVT_L
  IMOM_C = NMOM_C
  IMOM_R = NMOM_R
  IMOM_I = 2 ! Ni is diagnostic and always available
  IMOM_S = NMOM_S
  IMOM_G = NMOM_G
  IMOM_H = NMOM_H
END IF
!
ZRHO00 = CST%XP00 / (CST%XRD * PTHVREFZIKB)
ZCOR00 = ZRHO00**ZCEXVT
!
IF (LINDUCTIVE) ALLOCATE (ZEFIELDW(KMICRO))
!
!
!*      1.3   packing
!
!  optimization by looking for locations where
!  the microphysical fields are larger than a minimal value only !!!
!
IF (KMICRO >= 0) THEN
  IMICRO = COUNTJV(ODMICRO(:,:,:), II1(:), II2(:), II3(:))
  !
  ! some microphysical tendencies are optional: the corresponding 1D arrays must be allocated
  IF (HCLOUD(1:3) == 'ICE') THEN  ! ICE3 scheme
    ALLOCATE(ZRCMLTSR(IMICRO))
    ALLOCATE(ZRICFRR(IMICRO))
  END IF
  IF (HCLOUD == 'LIMA') THEN  ! LIMA scheme
    ALLOCATE(ZRVHENC(IMICRO))
    ALLOCATE(ZRCHINC(IMICRO))
    ALLOCATE(ZRVHONH(IMICRO))
    ALLOCATE(ZRRCVRC(IMICRO))
    ALLOCATE(ZRICNVI(IMICRO))
    ALLOCATE(ZRVDEPI(IMICRO))
    ALLOCATE(ZRSHMSI(IMICRO))
    ALLOCATE(ZRGHMGI(IMICRO))
    ALLOCATE(ZRICIBU(IMICRO))
    ALLOCATE(ZRIRDSF(IMICRO))
    ALLOCATE(ZRCCORR2(IMICRO))
    ALLOCATE(ZRRCORR2(IMICRO))
    ALLOCATE(ZRICORR2(IMICRO))
  END IF
  IF (KRR == 7) THEN ! hail activated
    ALLOCATE(ZRWETGH(IMICRO))
    ALLOCATE(ZRCWETH(IMICRO))
    ALLOCATE(ZRIWETH(IMICRO))
    ALLOCATE(ZRSWETH(IMICRO))
    ALLOCATE(ZRGWETH(IMICRO))
    ALLOCATE(ZRRWETH(IMICRO))
    ALLOCATE(ZRCDRYH(IMICRO))
    ALLOCATE(ZRRDRYH(IMICRO))
    ALLOCATE(ZRIDRYH(IMICRO))
    ALLOCATE(ZRSDRYH(IMICRO))
    ALLOCATE(ZRGDRYH(IMICRO))
    ALLOCATE(ZRHMLTR(IMICRO))
    ALLOCATE(ZRDRYHG(IMICRO))
  END IF
  !
  DO JL = 1, IMICRO
    ZZT(JL)      = PZT(II1(JL),II2(JL),II3(JL))
    ZRHODREF(JL) = PRHODREF(II1(JL),II2(JL),II3(JL))
    ZRHOCOR(JL)  = (ZRHO00 / ZRHODREF(JL))**ZCEXVT
    ZRHODJ(JL)   = PRHODJ(II1(JL),II2(JL),II3(JL))
    !
    ZCIT(JL)     = PCIT(II1(JL),II2(JL),II3(JL))
    IF (IMOM_C == 2) ZCCT(JL) = PCCT(II1(JL),II2(JL),II3(JL))
    IF (IMOM_R == 2) ZCRT(JL) = PCRT(II1(JL),II2(JL),II3(JL))
    IF (IMOM_S == 2) ZCST(JL) = PCST(II1(JL),II2(JL),II3(JL))
    IF (IMOM_G == 2) ZCGT(JL) = PCGT(II1(JL),II2(JL),II3(JL))
    IF (IMOM_H == 2) ZCHT(JL) = PCHT(II1(JL),II2(JL),II3(JL))
    IF (IMOM_C == 1) ZCCT(JL) = 0.
    IF (IMOM_R == 1) ZCRT(JL) = 0.
    IF (IMOM_S == 1) ZCST(JL) = 0.
    IF (IMOM_G == 1) ZCGT(JL) = 0.
    IF (IMOM_H == 1) ZCHT(JL) = 0.
    !
    ZRVT(JL)  = PRVT(II1(JL),II2(JL),II3(JL))
    ZRCT(JL)  = PRCT(II1(JL),II2(JL),II3(JL))
    ZRRT(JL)  = PRRT(II1(JL),II2(JL),II3(JL))
    ZRIT(JL)  = PRIT(II1(JL),II2(JL),II3(JL))
    ZRST(JL)  = PRST(II1(JL),II2(JL),II3(JL))
    ZRGT(JL)  = PRGT(II1(JL),II2(JL),II3(JL))
    IF (KRR == 7) ZRHT(JL) = PRHT(II1(JL),II2(JL),II3(JL))
    !
    ZQPIT(JL) = PQPIT(II1(JL),II2(JL),II3(JL))
    ZQNIT(JL) = PQNIT(II1(JL),II2(JL),II3(JL))
    ZQCT(JL)  = PQCT(II1(JL),II2(JL),II3(JL))
    ZQRT(JL)  = PQRT(II1(JL),II2(JL),II3(JL))
    ZQIT(JL)  = PQIT(II1(JL),II2(JL),II3(JL))
    ZQST(JL)  = PQST(II1(JL),II2(JL),II3(JL))
    ZQGT(JL)  = PQGT(II1(JL),II2(JL),II3(JL))
    IF (KRR == 7) ZQHT(JL) = PQHT(II1(JL),II2(JL),II3(JL))
    !
    ZQPIS(JL) = PQPIS(II1(JL), II2(JL), II3(JL))
    ZQNIS(JL) = PQNIS(II1(JL), II2(JL), II3(JL))
    ZQCS(JL)  = PQCS(II1(JL), II2(JL), II3(JL))
    ZQRS(JL)  = PQRS(II1(JL), II2(JL), II3(JL))
    ZQIS(JL)  = PQIS(II1(JL), II2(JL), II3(JL))
    ZQSS(JL)  = PQSS(II1(JL), II2(JL), II3(JL))
    ZQGS(JL)  = PQGS(II1(JL), II2(JL), II3(JL))
    IF (KRR == 7) ZQHS(JL) = PQHS(II1(JL), II2(JL), II3(JL))
    !
    IF (LINDUCTIVE) ZEFIELDW(JL) = XEFIELDW(II1(JL), II2(JL), II3(JL))
    !
    ! microphysical tendencies
    ZRVHENI(JL) = PRVHENI(II1(JL), II2(JL), II3(JL))
    ZRRHONG(JL) = PRRHONG(II1(JL), II2(JL), II3(JL))
    ZRIMLTC(JL) = PRIMLTC(II1(JL), II2(JL), II3(JL))
    ZRCHONI(JL) = PRCHONI(II1(JL), II2(JL), II3(JL))
    ZRVDEPS(JL) = PRVDEPS(II1(JL), II2(JL), II3(JL))
    ZRIAGGS(JL) = PRIAGGS(II1(JL), II2(JL), II3(JL))
    ZRIAUTS(JL) = PRIAUTS(II1(JL), II2(JL), II3(JL))
    ZRVDEPG(JL) = PRVDEPG(II1(JL), II2(JL), II3(JL))
    ZRCAUTR(JL) = PRCAUTR(II1(JL), II2(JL), II3(JL))
    ZRCACCR(JL) = PRCACCR(II1(JL), II2(JL), II3(JL))
    ZRREVAV(JL) = PRREVAV(II1(JL), II2(JL), II3(JL))
    ZRCRIMSS(JL) = PRCRIMSS(II1(JL), II2(JL), II3(JL))
    ZRCRIMSG(JL) = PRCRIMSG(II1(JL), II2(JL), II3(JL))
    ZRSRIMCG(JL) = PRSRIMCG(II1(JL), II2(JL), II3(JL))
    ZRRACCSS(JL) = PRRACCSS(II1(JL), II2(JL), II3(JL))
    ZRRACCSG(JL) = PRRACCSG(II1(JL), II2(JL), II3(JL))
    ZRSACCRG(JL) = PRSACCRG(II1(JL), II2(JL), II3(JL))
    ZRSMLTG(JL) = PRSMLTG(II1(JL), II2(JL), II3(JL))
    ZRICFRRG(JL) = PRICFRRG(II1(JL), II2(JL), II3(JL))
    ZRRCFRIG(JL) = PRRCFRIG(II1(JL), II2(JL), II3(JL))
    ZRCWETG(JL) = PRCWETG(II1(JL), II2(JL), II3(JL))
    ZRIWETG(JL) = PRIWETG(II1(JL), II2(JL), II3(JL))
    ZRRWETG(JL) = PRRWETG(II1(JL), II2(JL), II3(JL))
    ZRSWETG(JL) = PRSWETG(II1(JL), II2(JL), II3(JL))
    ZRCDRYG(JL) = PRCDRYG(II1(JL), II2(JL), II3(JL))
    ZRIDRYG(JL) = PRIDRYG(II1(JL), II2(JL), II3(JL))
    ZRRDRYG(JL) = PRRDRYG(II1(JL), II2(JL), II3(JL))
    ZRSDRYG(JL) = PRSDRYG(II1(JL), II2(JL), II3(JL))
    ZRGMLTR(JL) = PRGMLTR(II1(JL), II2(JL), II3(JL))
    ZRCBERI(JL) = PRCBERI(II1(JL), II2(JL), II3(JL))
    IF (HCLOUD(1:3) == 'ICE') THEN
      ZRCMLTSR(JL) = PRCMLTSR(II1(JL), II2(JL), II3(JL))
      ZRICFRR(JL)  = PRICFRR(II1(JL), II2(JL), II3(JL))
    END IF
    IF (HCLOUD == 'LIMA') THEN
      ZCST(JL)    = PCST(II1(JL), II2(JL), II3(JL))
      ZCGT(JL)    = PCGT(II1(JL), II2(JL), II3(JL))
      ZRVHENC(JL) = PRVHENC(II1(JL), II2(JL), II3(JL))
      ZRCHINC(JL) = PRCHINC(II1(JL), II2(JL), II3(JL))
      ZRVHONH(JL) = PRVHONH(II1(JL), II2(JL), II3(JL))
      ZRRCVRC(JL) = PRRCVRC(II1(JL), II2(JL), II3(JL))
      ZRICNVI(JL) = PRICNVI(II1(JL), II2(JL), II3(JL))
      ZRVDEPI(JL) = PRVDEPI(II1(JL), II2(JL), II3(JL))
      ZRSHMSI(JL) = PRSHMSI(II1(JL), II2(JL), II3(JL))
      ZRGHMGI(JL) = PRGHMGI(II1(JL), II2(JL), II3(JL))
      ZRICIBU(JL) = PRICIBU(II1(JL), II2(JL), II3(JL))
      ZRIRDSF(JL) = PRIRDSF(II1(JL), II2(JL), II3(JL))
      ZRCCORR2(JL) = PRCCORR2(II1(JL), II2(JL), II3(JL))
      ZRRCORR2(JL) = PRRCORR2(II1(JL), II2(JL), II3(JL))
      ZRICORR2(JL) = PRICORR2(II1(JL), II2(JL), II3(JL))
    END IF
    IF (KRR == 7) THEN
      ZCHT(JL)    = PCHT(II1(JL), II2(JL), II3(JL))
      ZRWETGH(JL) = PRWETGH(II1(JL), II2(JL), II3(JL))
      ZRCWETH(JL) = PRCWETH(II1(JL), II2(JL), II3(JL))
      ZRIWETH(JL) = PRIWETH(II1(JL), II2(JL), II3(JL))
      ZRSWETH(JL) = PRSWETH(II1(JL), II2(JL), II3(JL))
      ZRGWETH(JL) = PRGWETH(II1(JL), II2(JL), II3(JL))
      ZRRWETH(JL) = PRRWETH(II1(JL), II2(JL), II3(JL))
      ZRCDRYH(JL) = PRCDRYH(II1(JL), II2(JL), II3(JL))
      ZRRDRYH(JL) = PRRDRYH(II1(JL), II2(JL), II3(JL))
      ZRIDRYH(JL) = PRIDRYH(II1(JL), II2(JL), II3(JL))
      ZRSDRYH(JL) = PRSDRYH(II1(JL), II2(JL), II3(JL))
      ZRGDRYH(JL) = PRGDRYH(II1(JL), II2(JL), II3(JL))
      ZRHMLTR(JL) = PRHMLTR(II1(JL), II2(JL), II3(JL))
      ZRDRYHG(JL) = PRDRYHG(II1(JL), II2(JL), II3(JL))
    END IF
  END DO
  !
  ZRHOCOR(:) = (ZRHO00 / ZRHODREF(:))**ZCEXVT
!
!
!*      1.4   allocations for the non-inductive parameterizations
!
  IF (CNI_CHARGING == 'GARDI') THEN
    ALLOCATE( ZDELTALWC(KMICRO) )
    ALLOCATE( ZFT(KMICRO) )
  END IF
!
  IF (CNI_CHARGING == 'SAUN1' .OR. CNI_CHARGING == 'SAUN2' .OR. &
      CNI_CHARGING == 'TAKAH' .OR.                              &
      CNI_CHARGING == 'BSMP1' .OR. CNI_CHARGING == 'BSMP2' .OR. &
      CNI_CHARGING == 'TEEWC' .OR. CNI_CHARGING == 'TERAR') THEN
    ALLOCATE( ZEW(KMICRO) )
  END IF
!
  IF (CNI_CHARGING == 'SAUN1' .OR. CNI_CHARGING == 'SAUN2' .OR. &
      CNI_CHARGING == 'TAKAH' .OR. CNI_CHARGING == 'TEEWC') THEN
    ALLOCATE( ZDQ(KMICRO) )
  END IF
!
  IF (CNI_CHARGING == 'SAUN1' .OR. CNI_CHARGING == 'SAUN2' .OR. &
      CNI_CHARGING == 'TEEWC' )  THEN
    ALLOCATE( ZSAUNSK(KMICRO) )
    ALLOCATE( ZSAUNIM(KMICRO) )
    ALLOCATE( ZSAUNIN(KMICRO) )
    ALLOCATE( ZSAUNSM(KMICRO) )
    ALLOCATE( ZSAUNSN(KMICRO) )
  END IF
!
  IF (CNI_CHARGING == 'SAUN1' .OR. CNI_CHARGING == 'SAUN2' .OR. &
      CNI_CHARGING == 'SAP98' .OR.                              &
      CNI_CHARGING == 'BSMP1' .OR. CNI_CHARGING == 'BSMP2' .OR. &
      CNI_CHARGING == 'TEEWC' .OR. CNI_CHARGING == 'TERAR') THEN
    ALLOCATE( ZFQIAGGS(KMICRO) )
    ALLOCATE( ZFQIDRYGBS(KMICRO) )
    ALLOCATE( ZLBQSDRYGB1S(KMICRO) )
    ALLOCATE( ZLBQSDRYGB2S(KMICRO) )
    ALLOCATE( ZLBQSDRYGB3S(KMICRO) )
  END IF
!
  IF (CNI_CHARGING == 'SAP98' .OR. CNI_CHARGING == 'TERAR' .OR. &
      CNI_CHARGING == 'BSMP1' .OR. CNI_CHARGING == 'BSMP2') THEN
    ALLOCATE( ZRAR(KMICRO) )
    ALLOCATE( ZDQ_IS(KMICRO) )
    ALLOCATE( ZDQ_IG(KMICRO) )
    ALLOCATE( ZDQ_SG(KMICRO) )
    ALLOCATE( ZSAUNIM_IS(KMICRO) )
    ALLOCATE( ZSAUNIN_IS(KMICRO) )
    ALLOCATE( ZSAUNIM_IG(KMICRO) )
    ALLOCATE( ZSAUNIN_IG(KMICRO) )
    ALLOCATE( ZSAUNSK_SG(KMICRO) )
    ALLOCATE( ZSAUNSM_SG(KMICRO) )
    ALLOCATE( ZSAUNSN_SG(KMICRO) )
  END IF
!
!
!*      1.5   select parameters between ICEx and LIMA
!
  ALLOCATE(ZRTMIN(KRR))
  IF (HCLOUD(1:3) == 'ICE') THEN
! in ini_rain_ice, xrtmin is initialized with dimension 6 (hail not activated) or 7 (hail activated)
    ZRTMIN(1:KRR) = XRTMIN_I(1:KRR)
    !
    ZALPHAI = XALPHAI_I
    ZNUI = XNUI_I
    ZAI = XAI_I
    ZBI = XBI_I
    ZDS = XDS_I
    ZDG = XDG_I
    ZCXS = XCXS_I
    ZCXG = XCXG_I
    !
    ZCOLIS   = XCOLIS_I
    ZCOLEXIS = XCOLEXIS_I
    ZCOLIG   = XCOLIG_I
    ZCOLEXIG = XCOLEXIG_I
    ZCOLSG   = XCOLSG_I
    ZCOLEXSG = XCOLEXSG_I
    !
    IGAMINC = NGAMINC_I
    !
    IACCLBDAR = NACCLBDAR_I
    IACCLBDAS = NACCLBDAS_I
    ZACCINTP1S = XACCINTP1S_I
    ZACCINTP2S = XACCINTP2S_I
    ZACCINTP1R = XACCINTP1R_I
    ZACCINTP2R = XACCINTP2R_I
    ! 
    IDRYLBDAR = NDRYLBDAR_I
    IDRYLBDAS = NDRYLBDAS_I
    IDRYLBDAG = NDRYLBDAG_I
    ZDRYINTP1R = XDRYINTP1R_I
    ZDRYINTP2R = XDRYINTP2R_I
    ZDRYINTP1S = XDRYINTP1S_I
    ZDRYINTP2S = XDRYINTP2S_I
    ZDRYINTP1G = XDRYINTP1G_I
    ZDRYINTP2G = XDRYINTP2G_I
    !
    ZRIMINTP1 = XRIMINTP1_I
    ZRIMINTP2 = XRIMINTP2_I
    !
  ELSE IF (HCLOUD == 'LIMA') THEN
! in ini_lima, xrtmin is initialized with dimension 7
    ZRTMIN(1:KRR) = XRTMIN_L(1:KRR)
    !
    ZALPHAI = XALPHAI_L
    ZNUI = XNUI_L
    ZAI = XAI_L
    ZBI = XBI_L
    ZDS = XDS_L
    ZDG = XDG_L
    ZCXS = XCXS_L
    ZCXG = XCXG_L
    !
    ZCOLIS   = 0.25  ! variable not defined in LIMA, the value of ICEx is used 
    ZCOLEXIS = XCOLEXIS_L
    ZCOLIG   = XCOLIG_L
    ZCOLEXIG = XCOLEXIG_L
    ZCOLSG   = XCOLSG_L
    ZCOLEXSG = XCOLEXSG_L
    !
    IGAMINC = NGAMINC_L
    !
    IACCLBDAR = NACCLBDAR_L
    IACCLBDAS = NACCLBDAS_L
    ZACCINTP1S = XACCINTP1S_L
    ZACCINTP2S = XACCINTP2S_L
    ZACCINTP1R = XACCINTP1R_L
    ZACCINTP2R = XACCINTP2R_L
    ! 
    IDRYLBDAR = NDRYLBDAR_L
    IDRYLBDAS = NDRYLBDAS_L
    IDRYLBDAG = NDRYLBDAG_L
    ZDRYINTP1R = XDRYINTP1R_L
    ZDRYINTP2R = XDRYINTP2R_L
    ZDRYINTP1S = XDRYINTP1S_L
    ZDRYINTP2S = XDRYINTP2S_L
    ZDRYINTP1G = XDRYINTP1G_L
    ZDRYINTP2G = XDRYINTP2G_L
    !
    ZRIMINTP1 = XRIMINTP1_L
    ZRIMINTP2 = XRIMINTP2_L
  END IF
!
!
!*      1.6   update the slope parameter of the distribution
!*            and compute N_x if necessary
!
  IF (HCLOUD(1:3) == 'ICE') ZCCT(:) = 0.
  CALL COMPUTE_LAMBDA(2, IMOM_C, KMICRO, HCLOUD, ZRHODREF, ZRTMIN(2), ZRCT, ZCCT, ZLBDAC)
  CALL COMPUTE_LAMBDA(3, IMOM_R, KMICRO, HCLOUD, ZRHODREF, ZRTMIN(3), ZRRT, ZCRT, ZLBDAR)
  CALL COMPUTE_LAMBDA(4, IMOM_I, KMICRO, HCLOUD, ZRHODREF, ZRTMIN(4), ZRIT, ZCIT, ZLBDAI)
  CALL COMPUTE_LAMBDA(5, IMOM_S, KMICRO, HCLOUD, ZRHODREF, ZRTMIN(5), ZRST, ZCST, ZLBDAS)
  CALL COMPUTE_LAMBDA(6, IMOM_G, KMICRO, HCLOUD, ZRHODREF, ZRTMIN(6), ZRGT, ZCGT, ZLBDAG)
  IF (KRR == 7) CALL COMPUTE_LAMBDA(7, IMOM_H, KMICRO, HCLOUD, ZRHODREF, ZRTMIN(7), ZRHT, ZCHT, ZLBDAH)
!
!
!*      1.7   update the parameter e in the charge-diameter relationship
!
! Compute e_x at time t
  IF (HCLOUD == 'LIMA') THEN
    CALL ELEC_COMPUTE_EX(2, IMOM_C, KMICRO, HCLOUD, 1., ZRHODREF, ZRTMIN(2), ZRCT, ZQCT, ZECT, PLBDX=ZLBDAC, PCX=ZCCT)
    CALL ELEC_COMPUTE_EX(3, IMOM_R, KMICRO, HCLOUD, 1., ZRHODREF, ZRTMIN(3), ZRRT, ZQRT, ZERT, PLBDX=ZLBDAR, PCX=ZCRT)
    CALL ELEC_COMPUTE_EX(4, IMOM_I, KMICRO, HCLOUD, 1., ZRHODREF, ZRTMIN(4), ZRIT, ZQIT, ZEIT, PLBDX=ZLBDAI, PCX=ZCIT)
    CALL ELEC_COMPUTE_EX(5, IMOM_S, KMICRO, HCLOUD, 1., ZRHODREF, ZRTMIN(5), ZRST, ZQST, ZEST, PLBDX=ZLBDAS, PCX=ZCST)
    CALL ELEC_COMPUTE_EX(6, IMOM_G, KMICRO, HCLOUD, 1., ZRHODREF, ZRTMIN(6), ZRGT, ZQGT, ZEGT, PLBDX=ZLBDAG, PCX=ZCGT)
    IF (KRR == 7) CALL ELEC_COMPUTE_EX(7, IMOM_H, KMICRO, HCLOUD, 1., ZRHODREF, ZRTMIN(7), ZRHT, ZQHT, ZEHT, PLBDX=ZLBDAH, PCX=ZCHT)
  ELSE IF (HCLOUD(1:3) == 'ICE') THEN
    CALL ELEC_COMPUTE_EX(2, 1, KMICRO, HCLOUD, 1., ZRHODREF, ZRTMIN(2), ZRCT, ZQCT, ZECT)
    CALL ELEC_COMPUTE_EX(3, 1, KMICRO, HCLOUD, 1., ZRHODREF, ZRTMIN(3), ZRRT, ZQRT, ZERT, PLBDX=ZLBDAR)
    CALL ELEC_COMPUTE_EX(4, 1, KMICRO, HCLOUD, 1., ZRHODREF, ZRTMIN(4), ZRIT, ZQIT, ZEIT, PCX=ZCIT)
    CALL ELEC_COMPUTE_EX(5, 1, KMICRO, HCLOUD, 1., ZRHODREF, ZRTMIN(5), ZRST, ZQST, ZEST, PLBDX=ZLBDAS)
    CALL ELEC_COMPUTE_EX(6, 1, KMICRO, HCLOUD, 1., ZRHODREF, ZRTMIN(6), ZRGT, ZQGT, ZEGT, PLBDX=ZLBDAG)
    IF (KRR == 7) CALL ELEC_COMPUTE_EX(7, 1, KMICRO, HCLOUD, 1., ZRHODREF, ZRTMIN(7), ZRHT, ZQHT, ZEHT, PLBDX=ZLBDAH)
  END IF
!
!
!*      1.8   initialization for the non-inductive charging process
!
  SELECT CASE (CNI_CHARGING)
    ! Initialization for the parameterization of Gardiner et al. (1995)
    CASE ('GARDI')
      CALL ELEC_INIT_NOIND_GARDI()
      ! Save the effective water content
      DO JL = 1, KMICRO
        XEW(II1(JL),II2(JL),II3(JL)) = ZDELTALWC(JL) ! 
      END DO      
    !
    ! Initialization for the parameterizations of Saunders et al. (1991)
    ! with and without anomalies, and Tsenova and Mitzeva (2009)
    CASE ('SAUN1', 'SAUN2', 'TEEWC')
      CALL ELEC_INIT_NOIND_EWC()
      ! Save the effective water content
      DO JL = 1, KMICRO
        XEW(II1(JL),II2(JL),II3(JL)) = ZEW(JL) ! g/m3
      END DO
    !
    ! Initialization for the parameterizations of Saunders and Peck (1998), 
    ! Brooks et al. (1997) and Tsenova and Mitzeva (2011)
    CASE ('SAP98', 'BSMP1', 'BSMP2', 'TERAR')
      CALL ELEC_INIT_NOIND_RAR()
      ! Save the rime accretion rate (not recorded properly: 3 different RAR are computed !!!)
      DO JL = 1, KMICRO
        XEW(II1(JL),II2(JL),II3(JL)) = ZRAR(JL) ! g/m3
      END DO
    !
    ! Initialization for the parameterization of Takahashi (1978)
    CASE ('TAKAH')
      CALL ELEC_INIT_NOIND_TAKAH()
      ! Save the effective water content
      DO JL = 1, KMICRO
        XEW(II1(JL),II2(JL),II3(JL)) = ZEW(JL) ! g/m3
      END DO
  END SELECT
!
!
!------------------------------------------------------------------
!
!*      2.    COMPUTE THE SLOW COLD PROCESS SOURCES
!             -------------------------------------
!
!*      2.1   heterogeneous nucleation
!
! --> rien n'est fait pour l'elec pour le moment
! ICE3/4 : rvheni/rvhind
! LIMA : rvhenc, rchinc, rvhonh
!
!
!*      2.2   spontaneous freezing (rhong)
!
  if ( BUCONF%LBUDGET_SV ) then
    CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 2 ), 'SFR', &
                            Unpack( zqrs(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
    CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 5 ), 'SFR', &
                            Unpack( zqgs(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
  end if
!  
  ZWQ(:) = 0.
  WHERE (ZRRHONG(:) > 0. .AND. &
         ZRRT(:) > XRTMIN_ELEC(3) .AND. ABS(ZQRT(:)) > XQTMIN(3))
    ZWQ(:) = ZQRS(:)
    !
    ZQGS(:) = ZQGS(:) + ZQRS(:)
    ZQRS(:) = 0.
  END WHERE
!
  if ( BUCONF%LBUDGET_SV ) then
    CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 2 ), 'SFR', &
                           Unpack( zqrs(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
    CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 5 ), 'SFR', &
                           Unpack( zqgs(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
  end if
!
!
!*      2.3   cloud ice melting (rimltc)
!
  if ( BUCONF%LBUDGET_SV ) then
    CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 1 ), 'IMLT', &
                            Unpack( zqcs(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
    CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 3 ), 'IMLT', &
                            Unpack( zqis(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
  end if
!
  WHERE (ZRIMLTC(:) > 0.)
    ZQCS(:) = ZQCS(:) + ZQIS(:)
    ZQIS(:) = 0.
  END WHERE
!
  if ( BUCONF%LBUDGET_SV ) then
    CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 1 ), 'IMLT', &
                           Unpack( zqcs(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
    CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 3 ), 'IMLT', &
                           Unpack( zqis(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
  end if
!
!
!*      2.4   riming-conversion of the large sized aggregates into graupel ???
! ancienne param => on calcule plutot cette tendance un peu plus loin ?
!
!
!*      2.5   homogeneous nucleation (rchoni)
!
! CB : traitement different entre ice3 et lima --> a modifier eventuellement
!
  ZWQ(:) = 0.
  WHERE (ZRCHONI(:) > 0. .AND.           &
         ZRCT(:) > XRTMIN_ELEC(2) .AND.  &
         ABS(ZQCT(:)) > XQTMIN(2) .AND. ABS(ZECT(:)) > ELECP%XECMIN)
    ZWQ(:) = XQHON * ZECT(:) * ZRCHONI(:)
    ZWQ(:) = SIGN( MIN( ABS(ZQCT(:)/PTSTEP),ABS(ZWQ(:)) ),ZQCS(:) )
    !
    ZQIS(:) = ZQIS(:) + ZWQ(:)
    ZQCS(:) = ZQCS(:) - ZWQ(:)
  END WHERE
!
  if ( BUCONF%LBUDGET_SV ) then
    CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 1 ), 'HON', &
                           Unpack( -zwq(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
    CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 3 ), 'HON', &
                           Unpack(  zwq(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
  end if
!
!
!*      2.6   deposition on snow/aggregates (rvdeps)
!
  if ( BUCONF%LBUDGET_SV ) then
    CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG ), 'DEPS', &
                            Unpack( zqpis(:) * zrhodj(:), mask = odmicro(:, :, :), field =  0. ) )
    CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECEND ), 'DEPS', &
                            Unpack( zqnis(:) * zrhodj(:), mask = odmicro(:, :, :), field =  0. ) )
  end if
!
  ZWQ(:) = 0.
  !
  ! Only the sublimation of snow/aggregates is considered (negative part of PRVDEPS)
  WHERE (ZRVDEPS(:) < 0. .AND. & 
         ZRST(:) > XRTMIN_ELEC(5) .AND. ABS(ZQST(:)) > XQTMIN(5))
    ZWQ(:) = XCOEF_RQ_S * ZQST(:) * ZRVDEPS(:) / ZRST(:)
    ZWQ(:) = SIGN( MIN( ABS(ZQST(:)/PTSTEP),ABS(ZWQ(:)) ),ZQSS(:) )
    !
    ZQSS(:)  = ZQSS(:)  - ZWQ(:)
    ZQPIS(:) = ZQPIS(:) + MAX( 0.0,ZWQ(:)/ELECD%XECHARGE )
    ZQNIS(:) = ZQNIS(:) - MIN( 0.0,ZWQ(:)/ELECD%XECHARGE )
  END WHERE
!
  if ( BUCONF%LBUDGET_SV ) then
    CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG ), 'DEPS', &
                           Unpack( zqpis(:) * zrhodj(:), mask = odmicro(:, :, :), field =  0. ) )
    CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECEND ), 'DEPS', &
                           Unpack( zqnis(:) * zrhodj(:), mask = odmicro(:, :, :), field =  0. ) )
    CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 4 ), 'DEPS', &
                           Unpack( -zwq(:) * zrhodj(:),  mask = odmicro(:, :, :), field =  0. ) )
  end if
!
!
!*      2.7   aggregation on snow/aggregates (riaggs)
!
  CALL COMPUTE_CHARGE_TRANSFER (ZRIAGGS, ZRIT, ZQIT, PTSTEP, &
                                XRTMIN_ELEC(4), XQTMIN(4), XCOEF_RQ_I, &
                                ZWQ, ZQIS, ZQSS)
!
  if ( BUCONF%LBUDGET_SV ) then
    CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 3 ), 'AGGS', &
                           Unpack( -zwq(:), mask = odmicro(:, :, :), field =  0. ) )
    CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 4 ), 'AGGS', &
                           Unpack(  zwq(:), mask = odmicro(:, :, :), field =  0. ) )
  end if
!
!
!*      2.8   non-inductive charging during ice - snow collisions
!
  CALL ELEC_IAGGS_B()
!
  if ( BUCONF%LBUDGET_SV ) then
    CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 3 ), 'NIIS', &
                           Unpack( -zwq_ni(:), mask = odmicro(:, :, :), field =  0. ) )
    CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 4 ), 'NIIS', &
                           Unpack(  zwq_ni(:), mask = odmicro(:, :, :), field =  0. ) )
  end if
!
! Save the NI charging rate
  DO JL = 1, KMICRO
    XNI_IAGGS(II1(JL),II2(JL),II3(JL)) = ZWQ_NI(JL) * ZRHODREF(JL) ! C/m3/s
  END DO
!
!
!*      2.9   autoconversion of r_i for r_s production (riauts/ricnvs)
!
  CALL COMPUTE_CHARGE_TRANSFER (ZRIAUTS, ZRIT, ZQIT, PTSTEP, &
                                XRTMIN_ELEC(4), XQTMIN(4), XCOEF_RQ_I, &
                                ZWQ, ZQIS, ZQSS)
!
  if ( BUCONF%LBUDGET_SV ) then
    CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 3 ), 'AUTS', &
                           Unpack( -zwq(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
    CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 4 ), 'AUTS', &
                           Unpack(  zwq(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
  end if
!
!
!*      2.10  snow --> ice conversion (rscnvi)
!
  IF (HCLOUD == 'LIMA') THEN
    CALL COMPUTE_CHARGE_TRANSFER (ZRICNVI, ZRST, ZQST, PTSTEP, &
                                  XRTMIN_ELEC(5), XQTMIN(5), XCOEF_RQ_S, &
                                  ZWQ, ZQSS, ZQIS)          
!
    if ( BUCONF%LBUDGET_SV ) then
      CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 3 ), 'CNVI', &
                             Unpack(  zwq(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
      CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 4 ), 'CNVI', &
                             Unpack( -zwq(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
    end if
!
!*      2.11  water vapor deposition on ice crystals (rvdepi)
!
    if ( BUCONF%LBUDGET_SV ) then
      CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG ), 'SUBI', &
                              Unpack( zqpis(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
      CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECEND ), 'SUBI', &
                              Unpack( zqnis(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
    end if
    !
    ZWQ(:) = 0.
    !
    ! Only the sublimation of ice crystals is considered (negative part of PRVDEPI)
    WHERE (ZRVDEPI(:) < 0. .AND. &
           ZRIT(:) > XRTMIN_ELEC(4) .AND. ABS(ZQIT(:)) > XQTMIN(4))
      ZWQ(:) = XCOEF_RQ_I * ZQIT(:) * ZRVDEPI(:) / ZRIT(:)
      ZWQ(:) = SIGN( MIN( ABS(ZQIT(:)/PTSTEP),ABS(ZWQ(:)) ),ZQIS(:) )
      !
      ZQIS(:)  = ZQIS(:)  - ZWQ(:)
      ZQPIS(:) = ZQPIS(:) + MAX( 0.0,ZWQ(:)/ELECD%XECHARGE )
      ZQNIS(:) = ZQNIS(:) - MIN( 0.0,ZWQ(:)/ELECD%XECHARGE )
    END WHERE
!
    if ( BUCONF%LBUDGET_SV ) then
      CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG ), 'SUBI', &
                             Unpack( zqpis(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
      CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECEND ), 'SUBI', &
                             Unpack( zqnis(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )           
      CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 3 ), 'SUBI', &
                             Unpack( -zwq(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
    end if    
  END IF
!
!
!*      2.12  water vapor deposition on graupel (rvdepg)
!
  if ( BUCONF%LBUDGET_SV ) then
    CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG ), 'DEPG', &
                            Unpack( zqpis(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
    CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECEND ), 'DEPG', &
                            Unpack( zqnis(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
  end if
!
  ZWQ(:) = 0.
  !
  ! Only the sublimation of graupel is considered (negative part of PRVDEPG)
  WHERE (ZRVDEPG(:) < 0. .AND. &
         ZRGT(:) > XRTMIN_ELEC(6) .AND. ABS(ZQGT(:)) > XQTMIN(6))
    ZWQ(:) = XCOEF_RQ_G * ZQGT(:) * ZRVDEPG(:) / ZRGT(:)
    ZWQ(:) = SIGN( MIN( ABS(ZQGT(:)/PTSTEP),ABS(ZWQ(:)) ),ZQGS(:) )
    !
    ZQGS(:)  = ZQGS(:)  - ZWQ(:) 
    ZQPIS(:) = ZQPIS(:) + MAX( 0.0,ZWQ(:)/ELECD%XECHARGE )
    ZQNIS(:) = ZQNIS(:) - MIN( 0.0,ZWQ(:)/ELECD%XECHARGE )
  END WHERE
!
  if ( BUCONF%LBUDGET_SV ) then
    CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG ), 'DEPG', &
                           Unpack( zqpis(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
    CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECEND ), 'DEPG', &
                           Unpack( zqnis(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
    CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 5 ), 'DEPG', &
                           Unpack( -zwq(:) * zrhodj(:),  mask = odmicro(:, :, :), field = 0. ) )
  end if
!
!
!------------------------------------------------------------------
!
!*      3.    COMPUTE THE WARM PROCESS SOURCES
!             --------------------------------
!
!*      3.1   autoconversion of r_c for r_r production (rcautr)
!
  CALL COMPUTE_CHARGE_TRANSFER (ZRCAUTR, ZRCT, ZQCT, PTSTEP,           &
                                XRTMIN_ELEC(2), XQTMIN(2), XCOEF_RQ_C, &
                                ZWQ, ZQCS, ZQRS)
!
  if ( BUCONF%LBUDGET_SV ) then
    CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 1 ), 'AUTO', &
                           Unpack( -zwq(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
    CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 2 ), 'AUTO', &
                           Unpack(  zwq(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
  end if
!
!
!*      3.2   accretion of r_c for r_r production (rcaccr)
!
  CALL COMPUTE_CHARGE_TRANSFER (ZRCACCR, ZRCT, ZQCT, PTSTEP,           &
                                XRTMIN_ELEC(2), XQTMIN(2), XCOEF_RQ_C, &
                                ZWQ, ZQCS, ZQRS)
!
  if ( BUCONF%LBUDGET_SV ) then
    CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 1 ), 'ACCR', &
                           Unpack( -zwq(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
    CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 2 ), 'ACCR', &
                           Unpack(  zwq(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
  end if
!
!
!*      3.3   evaporation of raindrops (rrevav)
!
  if ( BUCONF%LBUDGET_SV ) then
    CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG ), 'REVA', &
                            Unpack( zqpis(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
    CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECEND ), 'REVA', &
                            Unpack( zqnis(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
  end if
!
  ZWQ(:) = 0.
  WHERE (ZRREVAV(:) > 0. .AND. &
         ZRRT(:) > XRTMIN_ELEC(3) .AND. ABS(ZQRT(:)) > XQTMIN(3))
    ZWQ(:) = XCOEF_RQ_R * ZQRT(:) * ZRREVAV(:) / ZRRT(:)
    ZWQ(:) = SIGN( MIN( ABS(ZQRT(:)/PTSTEP),ABS(ZWQ(:)) ),ZQRS(:) )
    !
    ZQRS(:)  = ZQRS(:)  - ZWQ(:)
    ZQPIS(:) = ZQPIS(:) + MAX( 0.0,ZWQ(:)/ELECD%XECHARGE )
    ZQNIS(:) = ZQNIS(:) - MIN( 0.0,ZWQ(:)/ELECD%XECHARGE )
  END WHERE
!
  if ( BUCONF%LBUDGET_SV ) then
    CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG ), 'REVA', &
                           Unpack( zqpis(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
    CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECEND ), 'REVA', &
                           Unpack( zqnis(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
    CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 2 ), 'REVA', &
                           Unpack( -zwq(:) * zrhodj(:),  mask = odmicro(:, :, :), field = 0. ) )
  end if
!
!
!*      3.4   conversion of drops to droplets (rrcvrc)
!
  IF (HCLOUD == 'LIMA') THEN
    if ( BUCONF%LBUDGET_SV ) then
      CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 1 ), 'R2C1', &
                              Unpack( zqcs(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
      CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 2 ), 'R2C1', &
                              Unpack( zqrs(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
    end if
!
    CALL COMPUTE_CHARGE_TRANSFER (ZRRCVRC, ZRRT, ZQRT, PTSTEP,           &
                                  XRTMIN_ELEC(3), XQTMIN(3), XCOEF_RQ_R, &
                                  ZWQ, ZQRS, ZQCS)
!
    if ( BUCONF%LBUDGET_SV ) then
      CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 1 ), 'R2C1', &
                             Unpack( zqcs(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
      CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 2 ), 'R2C1', &
                             Unpack( zqrs(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
    end if
  END IF
!
!------------------------------------------------------------------
!
!*      4.    COMPUTE THE FAST COLD PROCESS SOURCES FOR r_s
!             ---------------------------------------------
!
!*      4.1   cloud droplet riming of the aggregates
!
  if ( BUCONF%LBUDGET_SV ) then
    CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 1 ), 'RIM', &
                           Unpack( zqcs(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
    CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 4 ), 'RIM', &
                           Unpack( zqss(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
    CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 5 ), 'RIM', &
                           Unpack( zqgs(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
  end if
!
!*      4.1.1 riming of the small sized aggregates (rcrimss)
!
  CALL COMPUTE_CHARGE_TRANSFER (ZRCRIMSS, ZRCT, ZQCT, PTSTEP,          &
                                XRTMIN_ELEC(2), XQTMIN(2), XCOEF_RQ_C, &
                                ZWQ, ZQCS, ZQSS)
!
!
!*      4.1.2 riming conversion of the large sized aggregates into graupel (rcrimsg)
!
  ZWQ(:) = 0.
  WHERE (ZRCRIMSG(:) > 0. .AND.                                        &
         ZRCT(:) > XRTMIN_ELEC(2) .AND. ZRST(:) > XRTMIN_ELEC(5) .AND. &
         ABS(ZQCT(:)) > XQTMIN(2))
    ZWQ(:) = XCOEF_RQ_C * ZQCT(:) * ZRCRIMSG(:) / ZRCT(:)      ! QCRIMSG
    ZWQ(:) = SIGN( MIN( ABS(ZQCT(:)/PTSTEP),ABS(ZWQ(:)) ),ZQCS(:) )
    !
    ZQGS(:) = ZQGS(:) + ZWQ(:)
    ZQCS(:) = ZQCS(:) - ZWQ(:)
  END WHERE
!
!
!*      4.1.3 riming conversion of the large sized aggregates into graupel (rsrimcg)
!
  GMASK(:) = .FALSE.
  IGMASK = 0
  DO JJ = 1, SIZE(GMASK)
    IF (ZRSRIMCG(JJ) > 0. .AND. ZZT(JJ) < CST%XTT .AND. &
        ZRCT(JJ) > XRTMIN_ELEC(2) .AND. ZRST(JJ) > XRTMIN_ELEC(5) .AND. &
        ZLBDAS(JJ) > 0.) THEN  !++cb-- 12/07/23 condition ajoutee pour eviter log(0)
      IGMASK = IGMASK + 1
      I1(IGMASK) = JJ
      GMASK(JJ) = .TRUE.
    ELSE
      GMASK(JJ) = .FALSE.
    END IF
  END DO
  !
  ALLOCATE(ZVEC1(IGMASK))
  ALLOCATE(ZVEC2(IGMASK))
  ALLOCATE(IVEC2(IGMASK))
  !
  ! select the ZLBDAS
  DO JJ = 1, IGMASK
    ZVEC1(JJ) = ZLBDAS(I1(JJ))
  END DO
  ! find the next lower indice for the ZLBDAS in the geometrical set of Lbda_s 
  ! used to tabulate some moments of the incomplete gamma function
  ZVEC2(1:IGMASK) = MAX( 1.00001, MIN( REAL(IGAMINC)-0.00001,           &
                        ZRIMINTP1 * LOG( ZVEC1(1:IGMASK) ) + ZRIMINTP2 ) )
  IVEC2(1:IGMASK) = INT( ZVEC2(1:IGMASK) )
  ZVEC2(1:IGMASK) = ZVEC2(1:IGMASK) - REAL( IVEC2(1:IGMASK) )
  !
  ! perform the linear interpolation of the normalized "XFS"-moment of
  ! the incomplete gamma function
  ZVEC1(1:IGMASK) =  XGAMINC_RIM3( IVEC2(1:IGMASK)+1 ) *  ZVEC2(1:IGMASK)      &
                   - XGAMINC_RIM3( IVEC2(1:IGMASK)   ) * (ZVEC2(1:IGMASK) - 1.0)
  !
  ZWQ(:) = 0.
  DO JJ = 1, IGMASK
    ZWQ(I1(JJ)) = ZVEC1(JJ)
  END DO
  !
  DEALLOCATE(ZVEC1)
  DEALLOCATE(ZVEC2)
  DEALLOCATE(IVEC2)
  !
  ! riming-conversion of the large sized aggregates into graupeln (rsrimcg)
  WHERE (ZRSRIMCG(:) > 0. .AND.                                        &
         ZRCT(:) > XRTMIN_ELEC(2) .AND. ZRST(:) > XRTMIN_ELEC(5) .AND. &
         ABS(ZQCT(:)) > XQTMIN(2) .AND. ABS(ZEST(:)) > XESMIN)
    ZWQ(:) = XQSRIMCG * ZEST(:) * ZCST(:) *          &    ! QSRIMCG
             ZLBDAS(:)**XEXQSRIMCG * (1. - ZWQ(:)) / &
             (PTSTEP * ZRHODREF(:))
    ZWQ(:) = SIGN( MIN( ABS(ZQST(:)/PTSTEP),ABS(ZWQ(:)) ),ZQSS(:) )
    !
    ZQGS(:) = ZQGS(:) + ZWQ(:)
    ZQSS(:) = ZQSS(:) - ZWQ(:)
  END WHERE
!
  if ( BUCONF%LBUDGET_SV ) then
    CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 1 ), 'RIM', &
                           Unpack( zqcs(:) * zrhodj(:),  mask = odmicro(:, :, :), field = 0. ) )
    CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 4 ), 'RIM', &
                           Unpack( zqss(:) * zrhodj(:),  mask = odmicro(:, :, :), field = 0. ) )
    CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 5 ), 'RIM', &
                           Unpack( zqgs(:) * zrhodj(:),  mask = odmicro(:, :, :), field = 0. ) )
  end if
!
!
!*      4.2   Hallett-Mossop ice multiplication process due to snow riming (rhmsi)
!
  IF (HCLOUD == 'LIMA') THEN
    CALL COMPUTE_CHARGE_TRANSFER (ZRSHMSI, ZRST, ZQST, PTSTEP,          &
                                  XRTMIN_ELEC(4), XQTMIN(4), XCOEF_RQ_S, &
                                  ZWQ, ZQSS, ZQIS)       
!
    if ( BUCONF%LBUDGET_SV ) then
      CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 4 ), 'HMS', &
                             Unpack( -zwq(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
      CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 3 ), 'HMS', &
                             Unpack(  zwq(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
    end if
  END IF
!
!
!*      4.3   Raindrop accretion onto the aggregates
!
  IGMASK = 0
  DO JJ = 1, SIZE(GMASK)
    IF (ZRRT(JJ) > ZRTMIN(3) .AND. ZLBDAR(JJ) > 0. .AND. &
        ZRST(JJ) > ZRTMIN(5) .AND. ZLBDAS(JJ) > 0.) THEN
      IGMASK = IGMASK + 1
      I1(IGMASK) = JJ
      GMASK(JJ) = .TRUE.
    ELSE
      GMASK(JJ) = .FALSE.
    END IF
  END DO
  !
  IF (IGMASK > 0) THEN
    ALLOCATE(ZVEC1(IGMASK))
    ALLOCATE(ZVEC2(IGMASK))
    ALLOCATE(IVEC1(IGMASK))
    ALLOCATE(IVEC2(IGMASK))
    ALLOCATE(ZVECQ1(IGMASK))
    ALLOCATE(ZVECQ2(IGMASK))
    ALLOCATE(ZVECQ3(IGMASK))
    !
    ! select the (ZLBDAS,ZLBDAR) couplet
    DO JJ = 1, IGMASK
      ZVEC1(JJ) = ZLBDAS(I1(JJ))
      ZVEC2(JJ) = ZLBDAR(I1(JJ))
    END DO
    !
    ! find the next lower indice for the ZLBDAS and for the ZLBDAR in the geometrical 
    ! set of (Lbda_s,Lbda_r) couplet use to tabulate the kernels
    ZVEC1(1:IGMASK) = MAX( 1.00001, MIN( REAL(IACCLBDAS)-0.00001,           &
                          ZACCINTP1S * LOG( ZVEC1(1:IGMASK) ) + ZACCINTP2S ) )
    IVEC1(1:IGMASK) = INT( ZVEC1(1:IGMASK) )
    ZVEC1(1:IGMASK) = ZVEC1(1:IGMASK) - REAL( IVEC1(1:IGMASK) )
    !
    ZVEC2(1:IGMASK) = MAX( 1.00001, MIN( REAL(IACCLBDAR)-0.00001,           &
                          ZACCINTP1R * LOG( ZVEC2(1:IGMASK) ) + ZACCINTP2R ) )
    IVEC2(1:IGMASK) = INT( ZVEC2(1:IGMASK) )
    ZVEC2(1:IGMASK) = ZVEC2(1:IGMASK) - REAL( IVEC2(1:IGMASK) )
    !
    ! perform the bilinear interpolation of the normalized kernels
    ZVECQ1(:) = BI_LIN_INTP_V(XKER_Q_RACCSS, IVEC1, IVEC2, ZVEC1, ZVEC2, IGMASK)
    ZVECQ2(:) = BI_LIN_INTP_V(XKER_Q_RACCS,  IVEC1, IVEC2, ZVEC1, ZVEC2, IGMASK)
    ZVECQ3(:) = BI_LIN_INTP_V(XKER_Q_SACCRG, IVEC1, IVEC2, ZVEC1, ZVEC2, IGMASK)
    ZWQ1(:) = 0.
    ZWQ2(:) = 0.
    ZWQ3(:) = 0.
    DO JJ = 1, IGMASK
      ZWQ1(I1(JJ)) = ZVECQ1(JJ)
      ZWQ2(I1(JJ)) = ZVECQ2(JJ)
      ZWQ3(I1(JJ)) = ZVECQ3(JJ)
    END DO
!
    DEALLOCATE(ZVEC1)
    DEALLOCATE(ZVEC2)
    DEALLOCATE(IVEC1)
    DEALLOCATE(IVEC2)
    DEALLOCATE(ZVECQ1)
    DEALLOCATE(ZVECQ2)
    DEALLOCATE(ZVECQ3)
!
!
!*      4.3.1 raindrop accretion onto the small sized aggregates (rraccss)
!
    ZWQ4(:)   = 0.
    ZWQ5(:,:) = 0.
    WHERE (ZRRACCSS(:) > 0. .AND. &
           ZRRT(:) > XRTMIN_ELEC(3) .AND. ZRST(:) > XRTMIN_ELEC(5) .AND. &
           ZCRT(:) > 0.             .AND. ZCST(:) > 0.             .AND. &
           ZLBDAR(:) > 0.           .AND. ZLBDAS(:) > 0.           .AND. &
           ABS(ZERT(:)) > ELECP%XERMIN)  ! and zzt(:) < xtt ?
      ZWQ4(:) = XFQRACCS * ZERT(:) * ZRHOCOR(:) / (ZCOR00 * ZRHODREF(:)) * &
                ZCRT(:) * ZCST(:)                                        * &
               (XLBQRACCS1 * ZLBDAR(:)**(-2.0 - XFR)                     + &
                XLBQRACCS2 * ZLBDAR(:)**(-1.0 - XFR) * ZLBDAS(:)**(-1.0) + &
                XLBQRACCS3 * ZLBDAR(:)**(-XFR)       * ZLBDAS(:)**(-2.0))
      ZWQ5(:,1) = ZWQ4(:) * ZWQ1(:)                       ! QRACCSS
      ZWQ5(:,1) = SIGN( MIN( ABS(ZQRT(:)/PTSTEP),ABS(ZWQ5(:,1)) ),ZQRS(:) )
      !
      ZQRS(:) = ZQRS(:) - ZWQ5(:,1)
      ZQSS(:) = ZQSS(:) + ZWQ5(:,1)
    END WHERE
!
!
!*      4.3.2 raindrop accretion-conversion of the large sized aggregates into graupel 
!*            (rsaccrg & rraccsg)
!
    ZWQ5(:,2) = ZWQ2(:) * ZWQ4(:) ! QRACCS
    WHERE (ZRRACCSG(:) > 0. .AND. &
           ZRRT(:) > XRTMIN_ELEC(3) .AND. ZRST(:) > XRTMIN_ELEC(5) .AND. &
           ZLBDAR(:) > 0.           .AND. ZLBDAS(:) > 0.)
      ZWQ5(:,3) = ZWQ5(:,2) - ZWQ5(:,1)         ! QRACCSG
      ZWQ5(:,3) = SIGN( MIN( ABS(ZQRT(:)/PTSTEP),ABS(ZWQ5(:,3)) ),ZQRS(:) )
      !
      ZQRS(:) = ZQRS(:) - ZWQ5(:,3)
      ZQGS(:) = ZQGS(:) + ZWQ5(:,3)
    END WHERE
!
    WHERE (ZRSACCRG(:) > 0. .AND.                                        &
           ZRRT(:) > XRTMIN_ELEC(3) .AND. ZRST(:) > XRTMIN_ELEC(5) .AND. &
           ZCRT(:) > 0.             .AND. ZCST(:) > 0.             .AND. &
           ZLBDAR(:) > 0.           .AND. ZLBDAS(:) > 0.           .AND. &
           ABS(ZEST) > XESMIN) 
      ZWQ5(:,4) = ZWQ3(:) * XFQRACCS * ZEST(:) *                               &
                  ZRHOCOR(:) / (ZCOR00 * ZRHODREF(:)) *                        &
                  ZCRT(:) * ZCST(:) *                                          &
                 (XLBQSACCRG1 * ZLBDAS(:)**(-2.0 - XFS) +                      &
                  XLBQSACCRG2 * ZLBDAS(:)**(-1.0 - XFS) * ZLBDAR(:)**(-1.0) +  &
                  XLBQSACCRG3 * ZLBDAS(:)**(-XFS)       * ZLBDAR(:)**(-2.0)) ! QSACCR
      ZWQ5(:,4) = SIGN( MIN( ABS(ZQST(:)/PTSTEP),ABS(ZWQ5(:,4)) ),ZQSS(:) )
      !
      ZQSS(:) = ZQSS(:) - ZWQ5(:,4)
      ZQGS(:) = ZQGS(:) + ZWQ5(:,4)
    END WHERE
    !
    if ( BUCONF%LBUDGET_SV ) then
      CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 2 ), 'ACC', &
                             Unpack( (-zwq5(:,1) - zwq5(:,3)) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
      CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 4 ), 'ACC', &
                             Unpack( ( zwq5(:,1) - zwq5(:,4)) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
      CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 5 ), 'ACC', &
                             Unpack( ( zwq5(:,3) + zwq5(:,4)) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
    end if    
    !
  END IF  ! end if igmask>0
!
!
!*      4.4   conversion-melting of the aggregates (rsmltg)
!
  CALL COMPUTE_CHARGE_TRANSFER (ZRSMLTG, ZRST, ZQST, PTSTEP,           &
                                XRTMIN_ELEC(5), XQTMIN(5), XCOEF_RQ_S, &
                                ZWQ, ZQSS, ZQGS)
!
  if ( BUCONF%LBUDGET_SV ) then
    CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 4 ), 'CMEL', &
                           Unpack( -zwq(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
    CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 5 ), 'CMEL', &
                           Unpack(  zwq(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
  end if
!
!
!*      4.5   cloud droplet collection onto aggregates by positive temperature (rcmltsr)
!
  IF (HCLOUD(1:3) == 'ICE') THEN
    CALL COMPUTE_CHARGE_TRANSFER (ZRCMLTSR, ZRCT, ZQCT, PTSTEP,           &
                                  XRTMIN_ELEC(2), XQTMIN(2), XCOEF_RQ_C, &
                                  ZWQ, ZQCS, ZQRS)
!
    if ( BUCONF%LBUDGET_SV ) then
      CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 1 ), 'CMEL', &
                             Unpack( -zwq(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
      CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 2 ), 'CMEL', &
                             Unpack(  zwq(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
    end if
  END IF
!
!
!------------------------------------------------------------------
!
!*      5.    COMPUTE THE FAST COLD PROCESS SOURCES FOR r_g
!             ---------------------------------------------
!
!*      5.1   rain contact freezing (ricfrrg, rrcfrig, ricfrr)
!
  if ( BUCONF%LBUDGET_SV ) then
    CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 2 ), 'CFRZ', &
                           Unpack( zqrs(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
    CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 3 ), 'CFRZ', &
                           Unpack( zqis(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
    CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 5 ), 'CFRZ', &
                           Unpack( zqgs(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
  end if
  !
  ZWQ(:) = 0.
  WHERE (ZRRCFRIG(:) > 0. .AND.                                        &
         ZRIT(:) > XRTMIN_ELEC(4) .AND. ZRRT(:) > XRTMIN_ELEC(3) .AND. &
         ZCRT(:) > 0.             .AND.                                &
         ABS(ZERT(:)) > ELECP%XERMIN    .AND. ABS(ZQRT(:)) > XQTMIN(3))
    ZWQ(:) = XQRCFRIG * ZLBDAR(:)**XEXQRCFRIG * ZCIT(:) * ZCRT(:) * &
             ZERT(:) * ZRHOCOR(:) / (ZCOR00 * ZRHODREF(:))  ! QRCFRIG
    ZWQ(:) = SIGN( MIN( ABS(ZQRT(:)/PTSTEP),ABS(ZWQ(:)) ),ZQRS(:) )
    !
    ZQGS(:) = ZQGS(:) + ZWQ(:)
    ZQRS(:) = ZQRS(:) - ZWQ(:)
  END WHERE
  !
  ZWQ(:) = 0.
  WHERE (ZRICFRRG(:) > 0. .AND.                                        &
         ZRIT(:) > XRTMIN_ELEC(4) .AND. ZRRT(:) > XRTMIN_ELEC(3) .AND. &
         ABS(ZQIT(:)) > XQTMIN(4))
    ZWQ(:) = XCOEF_RQ_I * ZQIT(:) * ZRICFRRG(:) / ZRIT(:)     ! QICFRRG
    ZWQ(:) = SIGN( MIN( ABS(ZQIT(:)/PTSTEP),ABS(ZWQ(:)) ),ZQIS(:) )
    !
    ZQGS(:) = ZQGS(:) + ZWQ(:)
    ZQIS(:) = ZQIS(:) - ZWQ(:)
  ENDWHERE
!
!++CB-- 16/06/2022 il manque le traitement de qricfrr
!
  if ( BUCONF%LBUDGET_SV ) then
    CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 2 ), 'CFRZ', &
                           Unpack( zqrs(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
    CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 3 ), 'CFRZ', &
                           Unpack( zqis(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
    CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 5 ), 'CFRZ', &
                           Unpack( zqgs(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
  end if
!
!
!*      5.2   graupel dry growth (qcdryg, qrdryg, qidryg & qsdryg)
!
  ZWQ5(:,:) = 0.
!
!*      5.2.1 compute qcdryg
!
  WHERE (ZRCDRYG(:) > 0. .AND.                                         &
         ZRCT(:) > XRTMIN_ELEC(2) .AND. ZRGT(:) > XRTMIN_ELEC(6) .AND. &
         ABS(ZQCT(:)) > XQTMIN(2)) 
    ZWQ5(:,1) = XCOEF_RQ_C * ZQCT(:) * ZRCDRYG(:) / ZRCT(:)
    ZWQ5(:,1) = SIGN( MIN( ABS(ZQCT(:)/PTSTEP),ABS(ZWQ5(:,1)) ),ZQCS(:) )
    !
    ZQCS(:) = ZQCS(:) - ZWQ5(:,1)
    ZQGS(:) = ZQGS(:) + ZWQ5(:,1)
  ENDWHERE
!
!
!*      5.2.2 compute qidryg = qidryg_coal + qidryg_boun
!
  WHERE (ZRIDRYG(:) > 0. .AND.                                         &
         ZRIT(:) > XRTMIN_ELEC(4) .AND. ZRGT(:) > XRTMIN_ELEC(6) .AND. &
         ABS(ZQIT(:)) > XQTMIN(4))
    ZWQ5(:,2) = XCOEF_RQ_I * ZQIT(:) * ZRIDRYG(:) / ZRIT(:)         ! QIDRYG_coal
    ZWQ5(:,2) = SIGN( MIN( ABS(ZQIT(:)/PTSTEP),ABS(ZWQ5(:,2)) ),ZQIS(:) )
    !
    ZQIS(:) = ZQIS(:) - ZWQ5(:,2)
    ZQGS(:) = ZQGS(:) + ZWQ5(:,2)
  END WHERE
!
!
!*      5.2.3 compute non-inductive charging durig ice - graupel collisions
!
  ! charge separation during collision between ice and graupel
  CALL ELEC_IDRYG_B()  ! QIDRYG_boun
  !
  if ( BUCONF%LBUDGET_SV ) then
    CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 3 ), 'NIIG', &
                            Unpack( -zwq_ni(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
    CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 5 ), 'NIIG', &
                            Unpack(  zwq_ni(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
  end if
  !
  ! Save the NI charging rate
  DO JL = 1, KMICRO
    XNI_IDRYG(II1(JL),II2(JL),II3(JL)) = ZWQ_NI(JL) * ZRHODREF(JL) ! C/m3/s
  END DO
!
!
!*      5.2.4 compute qsdryg
!
  IGMASK = 0
  DO JJ = 1, SIZE(GMASK)
    IF (ZRST(JJ) > ZRTMIN(5) .AND. ZLBDAS(JJ) > 0. .AND. &
        ZRGT(JJ) > ZRTMIN(6) .AND. ZLBDAG(JJ) > 0.) THEN
      IGMASK = IGMASK + 1
      I1(IGMASK) = JJ
      GMASK(JJ) = .TRUE.
    ELSE
      GMASK(JJ) = .FALSE.
    END IF
  END DO
  !
  IF (IGMASK > 0) THEN
  !
    ALLOCATE(ZVEC1(IGMASK))
    ALLOCATE(ZVEC2(IGMASK))
    ALLOCATE(IVEC1(IGMASK))
    ALLOCATE(IVEC2(IGMASK))
    ALLOCATE(ZVECQ1(IGMASK))
    ALLOCATE(ZVECQ2(IGMASK))
    ALLOCATE(ZVECQ3(IGMASK))
    ALLOCATE(ZVECQ4(IGMASK))
  !
  ! select the (ZLBDAG,ZLBDAS) couplet
    DO JJ = 1, IGMASK
      ZVEC1(JJ) = ZLBDAG(I1(JJ))
      ZVEC2(JJ) = ZLBDAS(I1(JJ))
    END DO
  !
  ! find the next lower indice for the ZLBDAG and for the ZLBDAS in the geometrical set 
  ! of (Lbda_g,Lbda_s) couplet use to tabulate the SDRYG-kernel
    ZVEC1(1:IGMASK) = MAX(1.00001, MIN(REAL(IDRYLBDAG)-0.00001,  &
                          ZDRYINTP1G*LOG(ZVEC1(1:IGMASK))+ZDRYINTP2G))
    IVEC1(1:IGMASK) = INT(ZVEC1(1:IGMASK) )
    ZVEC1(1:IGMASK) = ZVEC1(1:IGMASK) - REAL(IVEC1(1:IGMASK))
    !
    ZVEC2(1:IGMASK) = MAX(1.00001, MIN( REAL(IDRYLBDAS)-0.00001, &
                          ZDRYINTP1S*LOG(ZVEC2(1:IGMASK))+ZDRYINTP2S))
    IVEC2(1:IGMASK) = INT(ZVEC2(1:IGMASK))
    ZVEC2(1:IGMASK) = ZVEC2(1:IGMASK) - REAL(IVEC2(1:IGMASK))
    !
    ! perform the bilinear interpolation of the normalized QSDRYG-kernels
    ! normalized Q-SDRYG-kernel 
    ZVECQ1(:) = BI_LIN_INTP_V(XKER_Q_SDRYG, IVEC1, IVEC2, ZVEC1, ZVEC2, IGMASK)
    ZWQ5(:,3) = 0.  ! normalement pas utile
    DO JJ = 1, IGMASK
      ZWQ5(I1(JJ),3) = ZVECQ1(JJ)
    END DO    
    !
    ! normalized Q-???-kernel
    IF (CNI_CHARGING == 'TAKAH' .OR. CNI_CHARGING == 'SAUN1' .OR. &
        CNI_CHARGING == 'SAUN2' .OR. CNI_CHARGING == 'SAP98' .OR. &
        CNI_CHARGING == 'GARDI' .OR.                              &
        CNI_CHARGING == 'BSMP1' .OR. CNI_CHARGING == 'BSMP2' .OR. &
        CNI_CHARGING == 'TEEWC' .OR. CNI_CHARGING == 'TERAR') THEN
      ZVECQ2(:)  = BI_LIN_INTP_V(XKER_Q_LIMSG, IVEC1, IVEC2, ZVEC1, ZVEC2, IGMASK) 
      ZWQ5(:,4) = 0.  ! normalement pas utile
      DO JJ = 1, IGMASK
        ZWQ5(I1(JJ),4) = ZVECQ2(JJ)
      END DO
    END IF
    !
    ! normalized Q-SDRYG-bouncing kernel
    IF (CNI_CHARGING == 'TAKAH' .OR. CNI_CHARGING == 'HELFA' .OR. &
        CNI_CHARGING == 'GARDI') THEN
      ZVECQ3(:)  = BI_LIN_INTP_V(XKER_Q_SDRYGB,IVEC1,IVEC2,ZVEC1,ZVEC2,IGMASK)
      ZWQ5(:,5) = 0.  ! normalement pas utile
      DO JJ = 1, IGMASK
        ZWQ5(I1(JJ),5) = ZVECQ3(JJ)
      END DO
    ELSE
      ZVECQ3(:) = BI_LIN_INTP_V(XKER_Q_SDRYGB1,IVEC1,IVEC2,ZVEC1,ZVEC2,IGMASK)
      ZVECQ4(:) = BI_LIN_INTP_V(XKER_Q_SDRYGB2,IVEC1,IVEC2,ZVEC1,ZVEC2,IGMASK)
      ZWQ5(:,6:7) = 0.  ! normalement pas utile
      DO JJ = 1, IGMASK
        ZWQ5(I1(JJ),6) = ZVECQ3(JJ) ! Dvqsgmn if charge>0
        ZWQ5(I1(JJ),7) = ZVECQ4(JJ) ! Dvqsgmn if charge<0
      END DO
    ENDIF
    !
    DEALLOCATE(ZVEC1)
    DEALLOCATE(ZVEC2)
    DEALLOCATE(IVEC1)
    DEALLOCATE(IVEC2)
    DEALLOCATE(ZVECQ1)
    DEALLOCATE(ZVECQ2)
    DEALLOCATE(ZVECQ3)
    DEALLOCATE(ZVECQ4)
!
!++CB-- CALCULER E_SG ICI POUR EVITER DES CALCULS REDONDANTS
    !
    ! compute QSDRYG_coal
    WHERE (ZRSDRYG(:) > 0 .AND.                                & !GDRY(:) .AND. &
           ZRST(:) > XRTMIN_ELEC(5) .AND.                      &
           ZLBDAS(:) > 0.           .AND. ZLBDAG(:) > 0. .AND. &
           ABS(ZQST(:)) > XQTMIN(5) .AND. ABS(ZEST(:)) > XESMIN)
      ZWQ5(:,3) = ZWQ5(:,3) * XFQSDRYG *                                   &
                  ZCOLSG * EXP(ZCOLEXSG * (ZZT(:) - CST%XTT)) *                &
                  ZEST(:) * ZRHOCOR(:) / (ZCOR00 * ZRHODREF(:)) *          & 
                  ZCGT(:) * ZCST(:) *                                      &
                 (XLBQSDRYG1 * ZLBDAS(:)**(-2.0-XFS) +                     &
                  XLBQSDRYG2 * ZLBDAS(:)**(-1.0-XFS) * ZLBDAG(:)**(-1.0) + &
                  XLBQSDRYG3 * ZLBDAS(:)**(-XFS)     * ZLBDAG(:)**(-2.0)) ! QSDRYG_coal
      ZWQ5(:,3) = SIGN( MIN( ABS(ZQST(:)/PTSTEP),ABS(ZWQ5(:,3)) ),ZQSS(:) )
      !
      ZQSS(:) = ZQSS(:) - ZWQ5(:,3)
      ZQGS(:) = ZQGS(:) + ZWQ5(:,3)
    ELSEWHERE
      ZWQ5(:,3) = 0.
    END WHERE
!
!
!*      5.2.5 compute non-inductive charging during snow - graupel collisions
!
    ! compute QSDRYG_boun
    CALL ELEC_SDRYG_B()
    !
    if ( BUCONF%LBUDGET_SV ) then
      CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 4 ), 'NISG', &
                              Unpack( -zwq_ni(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
      CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 5 ), 'NISG', &
                              Unpack(  zwq_ni(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
    end if
    !
    ! Save the NI charging rate
    DO JL = 1, KMICRO
      XNI_SDRYG(II1(JL),II2(JL),II3(JL)) = ZWQ_NI(JL) * ZRHODREF(JL) ! C/m3/s
    END DO
  END IF  ! end if igmask>0
!
!
!*      5.2.6 compute qrdryg
!
  IGMASK = 0
  GMASK(:) = .FALSE.
  DO JJ = 1, SIZE(GMASK)
    IF (ZRRT(JJ) > ZRTMIN(3) .AND. ZLBDAR(JJ) > 0. .AND. &
        ZRGT(JJ) > ZRTMIN(6) .AND. ZLBDAG(JJ) > 0.) THEN
      IGMASK = IGMASK + 1
      I1(IGMASK) = JJ
      GMASK(JJ) = .TRUE.
    ELSE
      GMASK(JJ) = .FALSE.
    END IF
  END DO
  !
  IF (IGMASK > 0) THEN
    !
    ALLOCATE(ZVEC1(IGMASK))
    ALLOCATE(ZVEC2(IGMASK))
    ALLOCATE(IVEC1(IGMASK))
    ALLOCATE(IVEC2(IGMASK))
    ALLOCATE(ZVECQ1(IGMASK))
    !
    ! select the (ZLBDAG,ZLBDAR) couplet
    DO JJ = 1, IGMASK
      ZVEC1(JJ) = ZLBDAG(I1(JJ))
      ZVEC2(JJ) = ZLBDAR(I1(JJ))
    END DO
    !
    ! find the next lower indice for the ZLBDAG and for the ZLBDAR in the geometrical set 
    ! of (Lbda_g,Lbda_r) couplet use to tabulate the QDRYG-kernel
    ZVEC1(1:IGMASK) = MAX(1.00001, MIN( REAL(IDRYLBDAG)-0.00001,           &
                         ZDRYINTP1G*LOG(ZVEC1(1:IGMASK))+ZDRYINTP2G))
    IVEC1(1:IGMASK) = INT(ZVEC1(1:IGMASK))
    ZVEC1(1:IGMASK) = ZVEC1(1:IGMASK) - REAL(IVEC1(1:IGMASK))
    !
    ZVEC2(1:IGMASK) = MAX(1.00001, MIN( REAL(IDRYLBDAR)-0.00001,           &
                          ZDRYINTP1R*LOG(ZVEC2(1:IGMASK))+ZDRYINTP2R))
    IVEC2(1:IGMASK) = INT(ZVEC2(1:IGMASK))
    ZVEC2(1:IGMASK) = ZVEC2(1:IGMASK) - REAL(IVEC2(1:IGMASK))
    !
    ! perform the bilinear interpolation of the normalized RDRYG-kernel
    ZVECQ1(:) = BI_LIN_INTP_V(XKER_Q_RDRYG, IVEC1, IVEC2, ZVEC1, ZVEC2, IGMASK)
    ZWQ5(:,4) = 0.
    DO JJ = 1, IGMASK
      ZWQ5(I1(JJ),4) = ZVECQ1(JJ)
    END DO
    !
    DEALLOCATE(ZVEC1)
    DEALLOCATE(ZVEC2)
    DEALLOCATE(IVEC1)
    DEALLOCATE(IVEC2)
    DEALLOCATE(ZVECQ1)  
    !
    ! compute QRDRYG
    WHERE (ZRRDRYG(:) > 0. .AND.                                         &
           ZRRT(:) > XRTMIN_ELEC(3) .AND. ZRGT(:) > XRTMIN_ELEC(6) .AND. & 
           ZCRT(:) > 0.             .AND. ZCGT(:) > 0.             .AND. & 
           ZLBDAR(:) > 0.           .AND. ZLBDAG(:) > 0.           .AND. &
           ABS(ZERT(:)) > ELECP%XERMIN .AND. ABS(ZQRT(:)) > XQTMIN(3))
      ZWQ5(:,4) = ZWQ5(:,4) * XFQRDRYG *                                       &
                  ZRHODREF(:)**(-ZCEXVT) *                                     &
                  ZERT(:) * ZCGT(:) * ZCRT(:) *                                &
                 (XLBQRDRYG1 * ZLBDAR(:)**(-2.0 - XFR) +                       &
                  XLBQRDRYG2 * ZLBDAR(:)**(-1.0 - XFR) * ZLBDAG(:)**(-1.0) +   &
                  XLBQRDRYG3 * ZLBDAR(:)**(-XFR)       * ZLBDAG(:)**(-2.0))
      ZWQ5(:,4) = SIGN( MIN( ABS(ZQRT(:)/PTSTEP),ABS(ZWQ5(:,4)) ),ZQRS(:) )
      !
      ZQRS(:) = ZQRS(:) - ZWQ5(:,4)
      ZQGS(:) = ZQGS(:) + ZWQ5(:,4)
    ELSEWHERE
      ZWQ5(:,4) = 0.
    ENDWHERE
!    ZRDRYG(:) = ZWQ5(:,1) + ZWQ5(:,2) + ZWQ5(:,3) + ZWQ5(:,4)
!
    if ( BUCONF%LBUDGET_SV ) then
      CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 1 ), 'DRYG', &
                              Unpack( -zwq5(:,1) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
      CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 2 ), 'DRYG', &
                              Unpack( -zwq5(:,4) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
      CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 3 ), 'DRYG', &
                              Unpack( -zwq5(:,2) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
      CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 4 ), 'DRYG', &
                              Unpack( -zwq5(:,3) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
      CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 5 ), 'DRYG', &
                              Unpack( (zwq5(:,1) + zwq5(:,2) + zwq5(:,3) + zwq5(:,4)) * zrhodj(:), &
                              mask = odmicro(:, :, :), field = 0. ) )
    end if
!
  END IF  ! end if igmask>0
!
!
!*      5.3   Hallett-Mossop ice multiplication process due to graupel riming (rhmgi)
!
  IF (HCLOUD == 'LIMA') THEN
    CALL COMPUTE_CHARGE_TRANSFER (ZRGHMGI, ZRGT, ZQGT, PTSTEP,          &
                                  XRTMIN_ELEC(6), XQTMIN(6), XCOEF_RQ_G, &
                                  ZWQ, ZQGS, ZQIS)          
!
    if ( BUCONF%LBUDGET_SV ) then
      CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 5 ), 'HMG', &
                             Unpack( -zwq(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
      CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 3 ), 'HMG', &
                             Unpack(  zwq(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
    end if
  END IF
!
!
!*      5.4   graupel wet growth (rcwetg, rrwetg, riwetg & rswetg)
!
!*      5.4.1 compute qcwetg
!
  ZWQ5(:,5) = 0.
  WHERE (ZRCWETG(:) > 0. .AND. ZRCT(:) > XRTMIN_ELEC(2) .AND. ABS(ZQCT(:)) > XQTMIN(2) .AND. &
         ZRGT(:) > XRTMIN_ELEC(6))
    ZWQ5(:,5) = XCOEF_RQ_C * ZRCWETG(:) * ZQCT(:) / ZRCT(:)
    ZWQ5(:,5) = SIGN( MIN( ABS(ZQCT(:)/PTSTEP),ABS(ZWQ5(:,5)) ),ZQCS(:) )
  END WHERE
!
!
!*      5.4.1 compute qiwetg
!
  ZWQ5(:,6) = 0.
  WHERE (ZRIWETG(:) > 0. .AND. ZRIT(:) > XRTMIN_ELEC(4) .AND. ABS(ZQIT(:)) > XQTMIN(4) .AND. &
         ZRGT(:) > XRTMIN_ELEC(6))
    ZWQ5(:,6) = XCOEF_RQ_I * ZRIWETG(:) * ZQIT(:) / ZRIT(:)
    ZWQ5(:,6) = SIGN( MIN( ABS(ZQIT(:)/PTSTEP),ABS(ZWQ5(:,6)) ),ZQIS(:) )
  END WHERE
!
!
!*      5.4.2 compute qswetg
!
  ZWQ5(:,7) = 0.
  WHERE (ZRSWETG(:) > 0. .AND. ZRST(:) > XRTMIN_ELEC(5) .AND. ABS(ZQST(:)) > XQTMIN(5) .AND. &
         ZRGT(:) > XRTMIN_ELEC(6))
    ZWQ5(:,7) = XCOEF_RQ_S * ZRSWETG(:) * ZQST(:) / ZRST(:)
    ZWQ5(:,7) = SIGN( MIN( ABS(ZQST(:)/PTSTEP),ABS(ZWQ5(:,7)) ),ZQSS(:) )
  END WHERE
!
!
!*      5.4.3 compute qrwetg
!
  ZWQ5(:,8) = 0.
  WHERE (ZRRWETG(:) > 0. .AND. ZRRT(:) > XRTMIN_ELEC(3) .AND. ABS(ZQRT(:)) > XQTMIN(3) .AND. &
         ZRGT(:) > XRTMIN_ELEC(6))
    ZWQ5(:,8) = XCOEF_RQ_R * ZQRT(:) * ZRRWETG(:) / ZRRT(:)
    ZWQ5(:,8) = SIGN( MIN( ABS(ZQRT(:)/PTSTEP),ABS(ZWQ5(:,8)) ),ZQRS(:) )
  ENDWHERE
!
!
!*      5.4.4 conversion of graupel into hail (rwetgh)
!
  if ( BUCONF%LBUDGET_SV ) then
    CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 5 ), 'WETG', &
                            Unpack( zqgs(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
  end if
!
  IF (KRR == 7) THEN
    ZWQ5(:,9) = 0.
    WHERE (ZRWETGH(:) > 0. .AND. ZRGT(:) > XRTMIN_ELEC(6) .AND. ABS(ZQGT(:)) > XQTMIN(6))
      ZWQ5(:,9) = XCOEF_RQ_G * ZQGT(:) * ZRWETGH(:) / ZRGT(:)
      ZWQ5(:,9) = SIGN( MIN( ABS(ZQGT(:)/PTSTEP),ABS(ZWQ5(:,9)) ),ZQGS(:) )
    END WHERE
    !
    WHERE (ZRCWETG(:) > 0. .OR. ZRRWETG(:) > 0. .OR. ZRIWETG(:) > 0. .OR. &
           ZRSWETG(:) > 0. .OR. ZRWETGH(:) > 0.)
      ZQCS(:) = ZQCS(:) - ZWQ5(:,5)
      ZQRS(:) = ZQRS(:) - ZWQ5(:,8)
      ZQIS(:) = ZQIS(:) - ZWQ5(:,6)
      ZQSS(:) = ZQSS(:) - ZWQ5(:,7)
      ZQGS(:) = ZQGS(:) + ZWQ5(:,5) + ZWQ5(:,8) + ZWQ5(:,6) + ZWQ5(:,7) - ZWQ5(:,9)
      ZQHS(:) = ZQHS(:) + ZWQ5(:,9)
    END WHERE
  ELSE IF (KRR == 6) THEN
    WHERE (ZRCWETG(:) > 0. .OR. ZRRWETG(:) > 0. .OR. ZRIWETG(:) > 0. .OR. &
           ZRSWETG(:) > 0.)
      ZQCS(:) = ZQCS(:) - ZWQ5(:,5)
      ZQRS(:) = ZQRS(:) - ZWQ5(:,8)
      ZQIS(:) = ZQIS(:) - ZWQ5(:,6)
      ZQSS(:) = ZQSS(:) - ZWQ5(:,7)
      ZQGS(:) = ZQGS(:) + ZWQ5(:,5) + ZWQ5(:,8) + ZWQ5(:,6) + ZWQ5(:,7)
    END WHERE
  END IF
!
  if ( BUCONF%LBUDGET_SV ) then
    CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 1 ), 'WETG', &
                           Unpack( -zwq5(:,5) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
    CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 2 ), 'WETG', &
                           Unpack( -zwq5(:,8) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
    CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 3 ), 'WETG', &
                           Unpack( -zwq5(:,6) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
    CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 4 ), 'WETG', &
                           Unpack( -zwq5(:,7) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
    CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 5 ), 'WETG', &
                           Unpack(  zqgs(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
    if ( krr == 7 ) &
      CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 6 ), 'WETG', &
                           Unpack( zwq5(:,9) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
  end if
!
!
!*      5.5   compute charge separation by the inductive mechanism
!
! Computation of the charge transfer rate during inductive mechanism
! Only the bouncing droplet-graupel collision when the graupel is in the dry 
! growth mode is considered
! The electric field is limited to 100 kV/m
!
  IF (LINDUCTIVE) THEN
    ZWQ(:) = 0.
    GMASK(:) = ZRCDRYG(:) > 0. 
    IGMASK   = COUNT(GMASK(:))
    !
    IF (IGMASK > 0) THEN
      ZWQ(:) = 0.
      !
      WHERE (GMASK(:) .AND.                                      &
             ZEFIELDW(:) /= 0. .AND. ABS(ZEGT(:)) > ELECP%XEGMIN .AND. &
             ZLBDAG(:) > 0. .AND. ZCGT(:) > 0. .AND.             &
             ZRGT(:) > XRTMIN_ELEC(6) .AND. ZRCT(:) > XRTMIN_ELEC(2))
        ZWQ(:) = XIND1 * ZCGT(:) * ZRHOCOR(:) *                             &
                (XIND2 * SIGN(MIN(100.E3, ABS(ZEFIELDW(:))), ZEFIELDW(:)) * &
                 ZLBDAG(:) **(-2.-ZDG) -                                    &
                 XIND3 * ZEGT(:) * ZLBDAG(:)**(-XFG-ZDG))
        ZWQ(:) = ZWQ(:) / ZRHODREF(:)
        !
        ZQGS(:) = ZQGS(:) + ZWQ(:)
        ZQCS(:) = ZQCS(:) - ZWQ(:)
      END WHERE
    !
      if ( BUCONF%LBUDGET_SV ) then
        CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 1 ), 'INCG', &
                               Unpack( -zwq(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
        CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 5 ), 'INCG', &
                               Unpack(  zwq(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
      end if    
      !
      ! Save the inductive charging rate
      DO JL = 1, KMICRO
        XIND_RATE(II1(JL),II2(JL),II3(JL)) = ZWQ(JL) * ZRHODREF(JL) ! C/m3/s
      END DO
    END IF
    !
    ! Save the inductive charging rate
    DO JL = 1, KMICRO
      XIND_RATE(II1(JL),II2(JL),II3(JL)) = ZWQ(JL) * ZRHODREF(JL) ! C/m3/s
    END DO
  END IF
!
!
!*      5.6   melting of the graupel (rgmltr)
!
  CALL COMPUTE_CHARGE_TRANSFER (ZRGMLTR, ZRGT, ZQGT, PTSTEP, &
                                XRTMIN_ELEC(6), XQTMIN(6), XCOEF_RQ_G, &
                                ZWQ, ZQGS, ZQRS)
!
  if ( BUCONF%LBUDGET_SV ) then
    CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 2 ), 'GMLT', &
                           Unpack(  zwq(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
    CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 5 ), 'GMLT', &
                           Unpack( -zwq(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
  end if
!
!
!------------------------------------------------------------------
!
!*      6.    COMPUTE THE OPTIONAL SECONDARY ICE PRODUCTION
!             ---------------------------------------------
!
! dans un premier temps, on considere que la charge echangee est proportionnelle
! a la masse echangee
!
!*      6.1   collisional ice breakup (cibu)
!
  IF (HCLOUD == 'LIMA' .AND. LCIBU) &
    CALL COMPUTE_CHARGE_TRANSFER (ZRICIBU, ZRST, ZQST, PTSTEP,           &
                                  XRTMIN_ELEC(5), XQTMIN(5), XCOEF_RQ_S, &
                                  ZWQ, ZQSS, ZQIS)
!
!*      6.2   raindrop shattering freezing (rdsf)
!
  IF (HCLOUD == 'LIMA' .AND. LRDSF) &
    CALL COMPUTE_CHARGE_TRANSFER (ZRIRDSF, ZRRT, ZQRT, PTSTEP,           &
                                  XRTMIN_ELEC(3), XQTMIN(3), XCOEF_RQ_R, &
                                  ZWQ, ZQRS, ZQIS)
!
!
!------------------------------------------------------------------
!
!*      7.    COMPUTE THE FAST COLD PROCESS SOURCES FOR r_h
!             ---------------------------------------------
!
  IF (KRR == 7) THEN
!
!*      7.1   wet growth of hail (qcweth, qrweth, qiweth, qsweth, qgweth)
!
    ZWQ5(:,:) = 0.
    !
    WHERE (ZRCWETH(:) > 0. .AND. ZRCT(:) > XRTMIN_ELEC(2))
      ZWQ5(:,1) = XCOEF_RQ_C * ZQCT(:) * ZRCWETH(:) / ZRCT(:)        ! QCWETH
      ZWQ5(:,1) = SIGN( MIN( ABS(ZQCT(:)/PTSTEP),ABS(ZWQ5(:,1)) ),ZQCS(:) )
      !
      ZQCS(:) = ZQCS(:) - ZWQ5(:,1)
      ZQHS(:) = ZQHS(:) + ZWQ5(:,1)
    END WHERE
    !
    WHERE (ZRIWETH(:) > 0. .AND. ZRIT(:) > XRTMIN_ELEC(4))
      ZWQ5(:,2) = XCOEF_RQ_I * ZQIT(:) * ZRIWETH(:) / ZRIT(:)        ! QIWETH
      ZWQ5(:,2) = SIGN( MIN( ABS(ZQIT(:)/PTSTEP),ABS(ZWQ5(:,2)) ),ZQIS(:) )
      !
      ZQIS(:) = ZQIS(:) - ZWQ5(:,2)
      ZQHS(:) = ZQHS(:) + ZWQ5(:,2)
    END WHERE
    !
    WHERE (ZRSWETH(:) > 0. .AND. ZRST(:) > XRTMIN_ELEC(5))
      ZWQ5(:,3) = XCOEF_RQ_S * ZQST(:) * ZRSWETH(:) / ZRST(:)        ! QSWETH
      ZWQ5(:,3) = SIGN( MIN( ABS(ZQST(:)/PTSTEP),ABS(ZWQ5(:,3)) ),ZQSS(:) )
      !
      ZQSS(:) = ZQSS(:) - ZWQ5(:,3)
      ZQHS(:) = ZQHS(:) + ZWQ5(:,3)
    END WHERE
    !
    WHERE (ZRGWETH(:) > 0. .AND. ZRGT(:) > XRTMIN_ELEC(6))
      ZWQ5(:,5) = XCOEF_RQ_G * ZQGT(:) * ZRGWETH(:) / ZRGT(:)        ! QGWETH
      ZWQ5(:,5) = SIGN( MIN( ABS(ZQGT(:)/PTSTEP),ABS(ZWQ5(:,5)) ),ZQGS(:) )
      !
      ZQGS(:) = ZQGS(:) - ZWQ5(:,5)
      ZQHS(:) = ZQHS(:) + ZWQ5(:,5)
    END WHERE
    !
    WHERE (ZRRWETH(:) > 0. .AND. ZRRT(:) > XRTMIN_ELEC(3))
      ZWQ5(:,4) = XCOEF_RQ_R * ZQRT(:) * ZRRWETH(:) / ZRRT(:)        ! QRWETH
      ZWQ5(:,4) = SIGN( MIN( ABS(ZQRT(:)/PTSTEP),ABS(ZWQ5(:,4)) ),ZQRS(:) )
      !
      ZQRS(:) = ZQRS(:) - ZWQ5(:,4)
      ZQHS(:) = ZQHS(:) + ZWQ5(:,4)
    END WHERE
!
    if ( BUCONF%LBUDGET_SV ) then
      CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 1 ), 'WETH', &
                             Unpack( -zwq5(:, 1) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
      CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 2 ), 'WETH', &
                             Unpack( -zwq5(:, 4) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
      CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 3 ), 'WETH', &
                             Unpack( -zwq5(:, 2) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
      CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 4 ), 'WETH', &
                             Unpack( -zwq5(:, 3) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
      CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 5 ), 'WETH', &
                             Unpack( -zwq5(:, 5) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
      CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 6 ), 'WETH',                      &
                             Unpack( ( zwq5(:, 1) + zwq5(:, 2) + zwq5(:, 3) + zwq5(:, 4) + zwq5(:, 5) ) &
                                       * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
    end if
!
!
!*      7.2   dry growth of hail (qcdryh, qrdryh, qidryh, qsdryh, qgdryh)
!
    ZWQ5(:,:) = 0.
    !
    WHERE (ZRCDRYH(:) > 0. .AND. ZRCT(:) > XRTMIN_ELEC(2))
      ZWQ5(:,1) = XCOEF_RQ_C * ZQCT(:) * ZRCDRYH(:) / ZRCT(:)        ! QCDRYH
      ZWQ5(:,1) = SIGN( MIN( ABS(ZQCT(:)/PTSTEP),ABS(ZWQ5(:,1)) ),ZQCS(:) )
      !
      ZQCS(:) = ZQCS(:) - ZWQ5(:,1)
      ZQHS(:) = ZQHS(:) + ZWQ5(:,1)
    END WHERE
    !
    WHERE (ZRIDRYH(:) > 0. .AND. ZRIT(:) > XRTMIN_ELEC(4))
      ZWQ5(:,2) = XCOEF_RQ_I * ZQIT(:) * ZRIDRYH(:) / ZRIT(:)        ! QIDRYH
      ZWQ5(:,2) = SIGN( MIN( ABS(ZQIT(:)/PTSTEP),ABS(ZWQ5(:,2)) ),ZQIS(:) )
      !
      ZQIS(:) = ZQIS(:) - ZWQ5(:,2)
      ZQHS(:) = ZQHS(:) + ZWQ5(:,2)
    END WHERE
    !
    WHERE (ZRSDRYH(:) > 0. .AND. ZRST(:) > XRTMIN_ELEC(5))
      ZWQ5(:,3) = XCOEF_RQ_S * ZQST(:) * ZRSDRYH(:) / ZRST(:)        ! QSDRYH
      ZWQ5(:,3) = SIGN( MIN( ABS(ZQST(:)/PTSTEP),ABS(ZWQ5(:,3)) ),ZQSS(:) )
      !
      ZQSS(:) = ZQSS(:) - ZWQ5(:,3)
      ZQHS(:) = ZQHS(:) + ZWQ5(:,3)
    END WHERE
    !
    WHERE (ZRGDRYH(:) > 0. .AND. ZRGT(:) > XRTMIN_ELEC(6))
      ZWQ5(:,5) = XCOEF_RQ_G * ZQGT(:) * ZRGDRYH(:) / ZRGT(:)        ! QGDRYH
      ZWQ5(:,5) = SIGN( MIN( ABS(ZQGT(:)/PTSTEP),ABS(ZWQ5(:,5)) ),ZQGS(:) )
      !
      ZQGS(:) = ZQGS(:) - ZWQ5(:,5)
      ZQHS(:) = ZQHS(:) + ZWQ5(:,5)
    END WHERE
    !
    WHERE (ZRRDRYH(:) > 0. .AND. ZRRT(:) > XRTMIN_ELEC(3))
      ZWQ5(:,4) = XCOEF_RQ_R * ZQRT(:) * ZRRDRYH(:) / ZRRT(:)        ! QRDRYH
      ZWQ5(:,4) = SIGN( MIN( ABS(ZQRT(:)/PTSTEP),ABS(ZWQ5(:,4)) ),ZQRS(:) )
      !
      ZQRS(:) = ZQRS(:) - ZWQ5(:,4)
      ZQHS(:) = ZQHS(:) + ZWQ5(:,4)
    END WHERE
!
!
!*      7.3   conversion of hail into graupel (qdryhg)
!
    CALL COMPUTE_CHARGE_TRANSFER (ZRDRYHG, ZRHT, ZQHT, PTSTEP, &
                                  XRTMIN_ELEC(7), XQTMIN(7), XCOEF_RQ_H, &
                                  ZWQ, ZQHS, ZQGS)
!
    if ( BUCONF%LBUDGET_SV ) then
      CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 1 ), 'DRYH', &
                             Unpack( -zwq5(:, 1) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
      CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 2 ), 'DRYH', &
                             Unpack( -zwq5(:, 4) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
      CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 3 ), 'DRYH', &
                             Unpack( -zwq5(:, 2) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
      CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 4 ), 'DRYH', &
                             Unpack( -zwq5(:, 3) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
      CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 5 ), 'DRYH', &
                             Unpack( (-zwq5(:, 5) - zwq(:)) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
      CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 6 ), 'DRYH',                      &
                             Unpack( ( zwq5(:, 1) + zwq5(:, 2) + zwq5(:, 3) + zwq5(:, 4) + zwq5(:, 5) + zwq(:) ) &
                                       * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
    end if
!
!
!*      7.4   melting of hail (qhmltr)
!
    CALL COMPUTE_CHARGE_TRANSFER (ZRHMLTR, ZRHT, ZQHT, PTSTEP, &
                                  XRTMIN_ELEC(7), XQTMIN(7), XCOEF_RQ_H, &
                                  ZWQ, ZQHS, ZQRS)
!
    if ( BUCONF%LBUDGET_SV ) then
      CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 2 ), 'HMLT', &
                             Unpack(  zwq(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
      CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 6 ), 'HMLT', &
                             Unpack( -zwq(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
    end if
!
  END IF  ! end if krr==7
!
!
!------------------------------------------------------------------
!
!*      8.    COMPUTE THE FAST COLD PROCESS SOURCES FOR r_i
!             ---------------------------------------------
!
!*      8.1   Bergeron-Findeisen effect (qcberi)
!
  CALL COMPUTE_CHARGE_TRANSFER (ZRCBERI, ZRCT, ZQCT, PTSTEP, &
                                XRTMIN_ELEC(2), XQTMIN(2), XCOEF_RQ_C, &
                                ZWQ, ZQCS, ZQIS)
!
  if ( BUCONF%LBUDGET_SV ) then
    CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 1 ), 'BERFI', &
                           Unpack( -zwq(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
    CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 3 ), 'BERFI', &
                           Unpack(  zwq(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
  end if
!
!
!------------------------------------------------------------------
!
!*      9.    COMPUTE THE CHARGE TRANSFER ASSOCIATED WITH THE CORRECTION TERM
!             ---------------------------------------------------------------
!
  IF (HCLOUD == 'LIMA') THEN
!
    if ( BUCONF%LBUDGET_SV ) then
      CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG ), 'CORR2', &
                              Unpack( zqpis(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
      CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECEND ), 'CORR2', &
                              Unpack( zqnis(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
    end if
!
    ZWQ1(:) = 0.
    WHERE (ZRCCORR2(:) .NE. 0. .AND. ZRCT(:) .GT. XRTMIN_ELEC(2))
      ZWQ1(:) = XCOEF_RQ_C * ZQCT(:) * ZRCCORR2(:) / ZRCT(:)
      ZWQ(:) = SIGN( MIN( ABS(ZQCT(:)/PTSTEP),ABS(ZWQ1(:)) ),ZQCS(:) )
      !
      ZQCS(:)  = ZQCS(:)  - ZWQ1(:)
      ZQPIS(:) = ZQPIS(:) + MAX( 0.0,ZWQ1(:)/ELECD%XECHARGE )
      ZQNIS(:) = ZQNIS(:) - MIN( 0.0,ZWQ1(:)/ELECD%XECHARGE )      
    END WHERE
    !
    !
    ZWQ2(:) = 0.
    WHERE (ZRRCORR2(:) .NE. 0. .AND. ZRRT(:) .GT. XRTMIN_ELEC(3))
      ZWQ2(:) = XCOEF_RQ_R * ZQRT(:) * ZRRCORR2(:) / ZRRT(:)
      ZWQ2(:) = SIGN( MIN( ABS(ZQRT(:)/PTSTEP),ABS(ZWQ2(:)) ),ZQRS(:) )
      !
      ZQRS(:)  = ZQRS(:)  - ZWQ2(:)
      ZQPIS(:) = ZQPIS(:) + MAX( 0.0,ZWQ2(:)/ELECD%XECHARGE )
      ZQNIS(:) = ZQNIS(:) - MIN( 0.0,ZWQ2(:)/ELECD%XECHARGE )      
    END WHERE
    !
    ZWQ3(:) = 0.
    WHERE (ZRICORR2(:) .NE. 0. .AND. ZRIT(:) .GT. XRTMIN_ELEC(4))
      ZWQ3(:) = XCOEF_RQ_I * ZQIT(:) * ZRICORR2(:) / ZRIT(:)
      ZWQ3(:) = SIGN( MIN( ABS(ZQIT(:)/PTSTEP),ABS(ZWQ3(:)) ),ZQIS(:) )
      !
      ZQIS(:)  = ZQIS(:)  - ZWQ3(:)
      ZQPIS(:) = ZQPIS(:) + MAX( 0.0,ZWQ3(:)/ELECD%XECHARGE )
      ZQNIS(:) = ZQNIS(:) - MIN( 0.0,ZWQ3(:)/ELECD%XECHARGE )      
    END WHERE
!
    if ( BUCONF%LBUDGET_SV ) then
      CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG ), 'CORR2', &
                             Unpack( zqpis(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
      CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECEND ), 'CORR2', &
                             Unpack( zqnis(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
      CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 1 ), 'CORR2', &
                             Unpack( zwq1(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
      CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 2 ), 'CORR2', &
                             Unpack( zwq2(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )                     
      CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + NSV_ELECBEG + 3 ), 'CORR2', &
                             Unpack( zwq3(:) * zrhodj(:), mask = odmicro(:, :, :), field = 0. ) )
    end if    
  END IF
!
!
!------------------------------------------------------------------
!
!*      X.    COMPUTE THE SEDIMENTATION SOURCE FOR Q_x
!             ----------------------------------------
!
! The sedimentation for electric charges is computed directly 
! in the microphysics scheme
!
!
!------------------------------------------------------------------
!
!*      10.   UPDATE VOLUMETRIC CHARGE CONCENTRATIONS
!             ---------------------------------------
!
  DO JL = 1, KMICRO
    PQPIS(II1(JL),II2(JL),II3(JL)) = ZQPIS(JL)
    PQNIS(II1(JL),II2(JL),II3(JL)) = ZQNIS(JL)
    PQCS (II1(JL),II2(JL),II3(JL)) = ZQCS(JL)
    PQRS (II1(JL),II2(JL),II3(JL)) = ZQRS(JL)
    PQIS (II1(JL),II2(JL),II3(JL)) = ZQIS(JL)
    PQSS (II1(JL),II2(JL),II3(JL)) = ZQSS(JL)
    PQGS (II1(JL),II2(JL),II3(JL)) = ZQGS(JL)
  END DO
  IF ( KRR == 7 ) THEN
    DO JL = 1, KMICRO
      PQHS(II1(JL),II2(JL),II3(JL)) = ZQHS(JL)
    END DO
  END IF
END IF  ! end if kmicro>0
!
!
!------------------------------------------------------------------
!
!*      11.   DEALLOCATE
!             ----------
!
IF (ALLOCATED( ZDELTALWC )) DEALLOCATE( ZDELTALWC )   
IF (ALLOCATED( ZFT ))       DEALLOCATE( ZFT )
!  
IF (ALLOCATED( ZEW ))          DEALLOCATE( ZEW )
IF (ALLOCATED( ZSAUNSK ))      DEALLOCATE( ZSAUNSK )
IF (ALLOCATED( ZSAUNIM ))      DEALLOCATE( ZSAUNIM )
IF (ALLOCATED( ZSAUNIN ))      DEALLOCATE( ZSAUNIN )
IF (ALLOCATED( ZSAUNSM ))      DEALLOCATE( ZSAUNSM )
IF (ALLOCATED( ZSAUNSN ))      DEALLOCATE( ZSAUNSN )
IF (ALLOCATED( ZFQIAGGS ))     DEALLOCATE( ZFQIAGGS )
IF (ALLOCATED( ZFQIDRYGBS ))   DEALLOCATE( ZFQIDRYGBS )
IF (ALLOCATED( ZLBQSDRYGB1S )) DEALLOCATE( ZLBQSDRYGB1S )
IF (ALLOCATED( ZLBQSDRYGB2S )) DEALLOCATE( ZLBQSDRYGB2S )
IF (ALLOCATED( ZLBQSDRYGB3S )) DEALLOCATE( ZLBQSDRYGB3S )
!
IF (ALLOCATED( ZDQ ))        DEALLOCATE( ZDQ )
IF (ALLOCATED( ZRAR ))       DEALLOCATE( ZRAR )
IF (ALLOCATED( ZDQ_IS ))     DEALLOCATE( ZDQ_IS )
IF (ALLOCATED( ZSAUNIM_IS )) DEALLOCATE( ZSAUNIM_IS )
IF (ALLOCATED( ZSAUNIN_IS )) DEALLOCATE( ZSAUNIN_IS )
IF (ALLOCATED( ZDQ_IG ))     DEALLOCATE( ZDQ_IG )
IF (ALLOCATED( ZSAUNIM_IG )) DEALLOCATE( ZSAUNIM_IG )
IF (ALLOCATED( ZSAUNIN_IG )) DEALLOCATE( ZSAUNIN_IG )
IF (ALLOCATED( ZDQ_SG ))     DEALLOCATE( ZDQ_SG )
IF (ALLOCATED( ZSAUNSK_SG )) DEALLOCATE( ZSAUNSK_SG )
IF (ALLOCATED( ZSAUNSM_SG )) DEALLOCATE( ZSAUNSM_SG )
IF (ALLOCATED( ZSAUNSN_SG )) DEALLOCATE( ZSAUNSN_SG )
!
IF (ALLOCATED( ZEFIELDW ))   DEALLOCATE( ZEFIELDW )
!
IF (ALLOCATED(ZRCMLTSR)) DEALLOCATE(ZRCMLTSR)
IF (ALLOCATED(ZRICFRR))  DEALLOCATE(ZRICFRR)
IF (ALLOCATED(ZRVHENC))  DEALLOCATE(ZRVHENC)
IF (ALLOCATED(ZRCHINC))  DEALLOCATE(ZRCHINC)
IF (ALLOCATED(ZRVHONH))  DEALLOCATE(ZRVHONH)
IF (ALLOCATED(ZRRCVRC))  DEALLOCATE(ZRRCVRC)
IF (ALLOCATED(ZRICNVI))  DEALLOCATE(ZRICNVI)
IF (ALLOCATED(ZRVDEPI))  DEALLOCATE(ZRVDEPI)
IF (ALLOCATED(ZRSHMSI))  DEALLOCATE(ZRSHMSI)
IF (ALLOCATED(ZRGHMGI))  DEALLOCATE(ZRGHMGI)
IF (ALLOCATED(ZRICIBU))  DEALLOCATE(ZRICIBU)
IF (ALLOCATED(ZRIRDSF))  DEALLOCATE(ZRIRDSF)
IF (ALLOCATED(ZRCCORR2)) DEALLOCATE(ZRCCORR2)
IF (ALLOCATED(ZRRCORR2)) DEALLOCATE(ZRRCORR2)
IF (ALLOCATED(ZRICORR2)) DEALLOCATE(ZRICORR2)
IF (ALLOCATED(ZRWETGH))  DEALLOCATE(ZRWETGH)
IF (ALLOCATED(ZRCWETH))  DEALLOCATE(ZRCWETH)
IF (ALLOCATED(ZRIWETH))  DEALLOCATE(ZRIWETH)
IF (ALLOCATED(ZRSWETH))  DEALLOCATE(ZRSWETH)
IF (ALLOCATED(ZRGWETH))  DEALLOCATE(ZRGWETH)
IF (ALLOCATED(ZRRWETH))  DEALLOCATE(ZRRWETH)
IF (ALLOCATED(ZRCDRYH))  DEALLOCATE(ZRCDRYH)
IF (ALLOCATED(ZRRDRYH))  DEALLOCATE(ZRRDRYH)
IF (ALLOCATED(ZRIDRYH))  DEALLOCATE(ZRIDRYH)
IF (ALLOCATED(ZRSDRYH))  DEALLOCATE(ZRSDRYH)
IF (ALLOCATED(ZRGDRYH))  DEALLOCATE(ZRGDRYH)
IF (ALLOCATED(ZRHMLTR))  DEALLOCATE(ZRHMLTR)
IF (ALLOCATED(ZRDRYHG))  DEALLOCATE(ZRDRYHG)
!
!------------------------------------------------------------------
END ASSOCIATE
!
CONTAINS
!
! - routines to initialize the non-inductive charging
! - routines to compute the non-inductive charging
! - various useful routines
!
!------------------------------------------------------------------
!
! ##################################
  SUBROUTINE ELEC_INIT_NOIND_GARDI()
! ##################################
!
!
! Purpose : initialization for the non-inductive charging process
!           following Gardiner et al. (1985)
!
!
!*	0.	DECLARATIONS
!		------------
!
IMPLICIT NONE
!
!
!*     1.      COMPUTE f(DeltaT) AND (LWC - LWC_crit)
!              --------------------------------------
!
GELEC(:,:) = .FALSE.
!
ZDELTALWC(:) = 0.
ZFT(:) = 0.
!
GELEC(:,3) = ZZT(:) > (CST%XTT - 40.) .AND. ZZT(:) < CST%XTT
GELEC(:,1) = GELEC(:,3) .AND.                                              &
             ZRIT(:) > XRTMIN_ELEC(4) .AND. ZRST(:) > XRTMIN_ELEC(5) .AND. &
             ZCIT(:) > 0.0 .AND. ZLBDAS(:) > 0. .AND. ZLBDAS(:) < XLBDAS_MAXE
GELEC(:,2) = GELEC(:,3) .AND.                                              &
             ZRIT(:) > XRTMIN_ELEC(4) .AND. ZRGT(:) > XRTMIN_ELEC(6) .AND. &
             ZCIT(:) > 0.0 .AND. ZLBDAG(:) > 0. .AND. ZLBDAG(:) < XLBDAG_MAXE
GELEC(:,3) = GELEC(:,3) .AND.                                              &
             ZRST(:) > XRTMIN_ELEC(5) .AND. ZRGT(:) > XRTMIN_ELEC(6) .AND. &
             ZLBDAS(:) > 0.           .AND. ZLBDAS(:) < XLBDAS_MAXE  .AND. &
             ZLBDAG(:) > 0.           .AND. ZLBDAG(:) < XLBDAG_MAXE
GELEC(:,4) = GELEC(:,1) .OR. GELEC(:,2) .OR. GELEC(:,3)
!
WHERE (GELEC(:,4))
  ! f(DeltaT)
  ZFT(:) = - 1.7E-5 * ((-21 / (XQTC - CST%XTT)) * (ZZT(:) - CST%XTT))**3   &
           - 0.003  * ((-21 / (XQTC - CST%XTT)) * (ZZT(:) - CST%XTT))**2   &
           - 0.05   * ((-21 / (XQTC - CST%XTT)) * (ZZT(:) - CST%XTT))      &
           + 0.13
  !
  ! LWC - LWC_crit
  ZDELTALWC(:) = (ZRCT(:) * ZRHODREF(:) * 1.E3) - XLWCC   ! (g m^-3)
ENDWHERE
!
END SUBROUTINE ELEC_INIT_NOIND_GARDI
!
!-----------------------------------------------------------------
!
! ################################
  SUBROUTINE ELEC_INIT_NOIND_EWC()
! ################################
!
!
! Purpose : initialization for the non-inductive charging process
!           following Saunders et al. (1991)
!
!
!*	0.	DECLARATIONS
!		------------
!
IMPLICIT NONE
!
!
!*      1.      PARAMETERS FOR POSITIVE NI CHARGING
!               -----------------------------------
!
GELEC(:,:) = .FALSE.
ZDQ(:) = 0.
ZEW(:) = 0.
!
! positive case is the default value
IF (CNI_CHARGING == 'SAUN1' .OR. CNI_CHARGING == 'SAUN2') THEN
  ZFQIAGGS(:)   = XFQIAGGSP
  ZFQIDRYGBS(:) = XFQIDRYGBSP
ELSE IF (CNI_CHARGING == 'TEEWC') THEN
  ZFQIAGGS(:)   = XFQIAGGSP_TAK
  ZFQIDRYGBS(:) = XFQIDRYGBSP_TAK
END IF
ZLBQSDRYGB1S(:) = XLBQSDRYGB1SP
ZLBQSDRYGB2S(:) = XLBQSDRYGB2SP
ZLBQSDRYGB3S(:) = XLBQSDRYGB3SP
ZSAUNIM(:) = XIMP     !3.76
ZSAUNIN(:) = XINP     !2.5
ZSAUNSK(:) = XSKP     !52.8
ZSAUNSM(:) = XSMP     !0.44
ZSAUNSN(:) = XSNP     !2.5
!
!
!*      2.      PARAMETERS FOR NEGATIVE NI CHARGING
!               -----------------------------------
!
! Mansell et al. (2005, JGR): droplet collection efficiency of the graupel ~ 0.6-1.0
WHERE (ZLBDAG(:) > 0. .AND. ZRCT(:) > 0.)
  ZEW(:) = 0.8 * ZRCT(:) * ZRHODREF(:) * 1.E3   ! (g m^-3)
END WHERE
!
GELEC(:,3) = ZZT(:) > (CST%XTT - 40.) .AND. ZZT(:) <= CST%XTT .AND. &
             ZEW(:) >= 0.01 .AND. ZEW(:) <= 10.
GELEC(:,1) = GELEC(:,3) .AND.                                              &
             ZRIT(:) > XRTMIN_ELEC(4) .AND. ZRST(:) > XRTMIN_ELEC(5) .AND. &
             ZCIT(:) > 0.0 .AND. ZLBDAS(:) > 0. .AND. ZLBDAS(:) < XLBDAS_MAXE
GELEC(:,2) = GELEC(:,3) .AND.                                              &
             ZRIT(:) > XRTMIN_ELEC(4) .AND. ZRGT(:) > XRTMIN_ELEC(6) .AND. &
             ZCIT(:) > 0.0 .AND. ZLBDAG(:) > 0. .AND. ZLBDAG(:) < XLBDAG_MAXE
GELEC(:,3) = GELEC(:,3) .AND.                                              &
             ZRST(:) > XRTMIN_ELEC(5) .AND. ZRGT(:) > XRTMIN_ELEC(6) .AND. &
             ZLBDAS(:) > 0.           .AND. ZLBDAS(:) < XLBDAS_MAXE  .AND. &
             ZLBDAG(:) > 0.           .AND. ZLBDAG(:) < XLBDAG_MAXE
!
GELEC(:,4) = GELEC(:,1) .OR. GELEC(:,2) .OR. GELEC(:,3)
!
IF (COUNT(GELEC(:,4)) > 0) THEN
  IF (CNI_CHARGING == 'SAUN1' .OR. CNI_CHARGING == 'SAUN2') THEN
    CALL INTERP_DQ_TABLE (KMICRO, NIND_TEMP, NIND_LWC, GELEC(:,4), &
                          ZEW, ZZT, XSAUNDER, ZDQ)
    !
    WHERE (ZDQ(:) < 0.)
      ZFQIAGGS(:)   = XFQIAGGSN
      ZFQIDRYGBS(:) = XFQIDRYGBSN
    END WHERE
  ELSE IF (CNI_CHARGING == 'TEEWC') THEN
    CALL INTERP_DQ_TABLE (KMICRO, NIND_TEMP, NIND_LWC, GELEC(:,4), &
                          ZEW, ZZT, XTAKA_TM, ZDQ)
    !
    WHERE (ZDQ(:) < 0.)
      ZFQIAGGS(:)   = XFQIAGGSN_TAK
      ZFQIDRYGBS(:) = XFQIDRYGBSN_TAK
    END WHERE
  END IF
!
! value of the parameters for the negative case
  WHERE (ZDQ(:) < 0.)
    ZLBQSDRYGB1S(:) = XLBQSDRYGB1SN
    ZLBQSDRYGB2S(:) = XLBQSDRYGB2SN
    ZLBQSDRYGB3S(:) = XLBQSDRYGB3SN
    ZSAUNIM(:) = XIMN      !2.54
    ZSAUNIN(:) = XINN      !2.8
    ZSAUNSK(:) = XSKN      !24.
    ZSAUNSM(:) = XSMN      !0.5
    ZSAUNSN(:) = XSNN      !2.8
  ENDWHERE
ENDIF
!
END SUBROUTINE ELEC_INIT_NOIND_EWC
!
!------------------------------------------------------------------
!
! ################################
  SUBROUTINE ELEC_INIT_NOIND_RAR()
! ################################
!
!
! Purpose : initialization for the non-inductive charging process
!           following Saunders and Peck (1998)
!
!
!*	0.	DECLARATIONS
!		------------
!
IMPLICIT NONE
!
! local variables
REAL, DIMENSION(KMICRO) :: ZRAR_CRIT     ! critical rime accretion rate
REAL, DIMENSION(KMICRO) :: ZVGMEAN, &    ! mean velocity of graupel 
                           ZVSMEAN       ! mean velocity of snow
!
!
!*      1.      COMPUTE THE CRITICAL RIME ACCRETION RATE
!               ----------------------------------------
!
ZRAR_CRIT(:) = 0.
!
IF (CNI_CHARGING == 'SAP98') THEN
!
  WHERE (ZZT(:) <= CST%XTT .AND. ZZT(:) >= (CST%XTT - 23.7)) ! Original from SAP98
    ZRAR_CRIT(:) = 1.0 + 7.93E-2 * (ZZT(:) - CST%XTT) +    &
                         4.48E-2 * (ZZT(:) - CST%XTT)**2 + &
                         7.48E-3 * (ZZT(:) - CST%XTT)**3 + &
                         5.47E-4 * (ZZT(:) - CST%XTT)**4 + &
                         1.67E-5 * (ZZT(:) - CST%XTT)**5 + &
                         1.76E-7 * (ZZT(:) - CST%XTT)**6
  END WHERE
  !
  WHERE (ZZT(:) < (CST%XTT - 23.7) .AND. ZZT(:) > (CST%XTT - 40.)) ! Added by Mansell
    ZRAR_CRIT(:) = 3.4 * (1.0 - (ABS(ZZT(:) - CST%XTT + 23.7) / & ! et al. (2005)
                   (-23.7 + 40.))**3.)
  END WHERE
  !
  GELEC(:,3) = ZZT(:) >= (CST%XTT - 40.) .AND. ZZT(:) <= CST%XTT
!
ELSE IF (CNI_CHARGING == 'BSMP1' .OR. CNI_CHARGING == 'BSMP2') THEN
!
  WHERE (ZZT(:) > (CST%XTT - 10.7))
    ZRAR_CRIT(:) = 0.66
  END WHERE
  WHERE (ZZT(:) <= (CST%XTT - 10.7) .AND. ZZT(:) >= (CST%XTT - 23.7))
    ZRAR_CRIT(:) = -1.47 - 0.2 * (ZZT(:) - CST%XTT)
  END WHERE
  WHERE (ZZT(:) < (CST%XTT - 23.7) .AND. ZZT(:) > (CST%XTT - 40.))
    ZRAR_CRIT(:) = 3.3
  END WHERE
  !
  GELEC(:,3) = ZZT(:) > (CST%XTT - 40.) .AND. ZZT(:) <= CST%XTT .AND. &
               ZEW(:) >= 0.01 .AND. ZEW(:) <= 10.
!
ELSE IF (CNI_CHARGING == 'TERAR') THEN
!
  GELEC(:,3) = ZZT(:) >= (CST%XTT - 40.) .AND. ZZT(:) <= CST%XTT
END IF
!
GELEC(:,1) = GELEC(:,3) .AND.                                              &
             ZRIT(:) > XRTMIN_ELEC(4) .AND. ZRST(:) > XRTMIN_ELEC(5) .AND. &
             ZCIT(:) > 0.0 .AND. ZLBDAS(:) > 0. .AND. ZLBDAS(:) < XLBDAS_MAXE
GELEC(:,2) = GELEC(:,3) .AND.                                              &
             ZRIT(:) > XRTMIN_ELEC(4) .AND. ZRGT(:) > XRTMIN_ELEC(6) .AND. &
             ZCIT(:) > 0.0 .AND. ZLBDAG(:) > 0. .AND. ZLBDAG(:) < XLBDAG_MAXE
GELEC(:,3) = GELEC(:,3) .AND.                                              &
             ZRST(:) > XRTMIN_ELEC(5) .AND. ZRGT(:) > XRTMIN_ELEC(6) .AND. &
             ZLBDAS(:) > 0.           .AND. ZLBDAS(:) < XLBDAS_MAXE  .AND. &
             ZLBDAG(:) > 0.           .AND. ZLBDAG(:) < XLBDAG_MAXE
!
!
!*      2.      INITIALIZATION FOR ICE CRYSTAL - GRAUPEL COLLISIONS
!               ---------------------------------------------------
!
ZDQ_IG(:) = 0.
!
! positive case is the default value
ZSAUNIM_IG(:) = XIMP
ZSAUNIN_IG(:) = XINP
!
! Compute the Rime Accretion Rate 
ZRAR(:) = 0.
ZVGMEAN(:) = 0.
WHERE (ZLBDAG(:) > 0. .AND. ZRCT(:) > 0.)
  ZVGMEAN(:) = XVGCOEF * ZRHOCOR(:) * ZLBDAG(:)**(-ZDG)
  ZRAR(:)    = 0.8 * ZRHODREF(:) * ZRCT(:) * ZVGMEAN(:) * 1.E3
END WHERE
!
IF (CNI_CHARGING == 'TERAR') THEN
  GELEC(:,2) = GELEC(:,2) .AND. ZRAR(:) > 0.01 .AND. ZRAR(:) <= 80.
ELSE    
  GELEC(:,2) = GELEC(:,2) .AND. ZRAR(:) > 0.1 
END IF
GELEC(:,4) = GELEC(:,2)
!
IF (COUNT(GELEC(:,4)) .GT. 0) THEN
! compute the coefficients for I-G collisions
  IF (CNI_CHARGING == 'SAP98') THEN
    CALL ELEC_INI_NI_SAP98 (KMICRO, GELEC(:,4), ZRAR, ZRAR_CRIT, ZDQ_IG)
    !
  ELSE IF (CNI_CHARGING == 'BSMP1' .OR. CNI_CHARGING == 'BSMP2') THEN
    ZRAR(:) = ZRAR(:) / 3
    CALL INTERP_DQ_TABLE (KMICRO, NIND_TEMP, NIND_LWC, GELEC(:,4), &
                          ZRAR, ZZT, XSAUNDER, ZDQ_IG)
    !
  ELSE IF (CNI_CHARGING == 'TERAR') THEN
    ZRAR(:) = ZRAR(:) / 8.
    CALL INTERP_DQ_TABLE (KMICRO, NIND_TEMP, NIND_LWC, GELEC(:,4), &
                          ZRAR, ZZT, XTAKA_TM, ZDQ_IG)
    !
  END IF
    !
  WHERE (ZDQ_IG(:) < 0.)
    ZSAUNIM_IG(:) = XIMN
    ZSAUNIN_IG(:) = XINN
  ENDWHERE
ENDIF
!
!
!*      3.      INITIALIZATION FOR ICE CRYSTAL - SNOW COLLISIONS
!               ------------------------------------------------
!
ZDQ_IS(:)  = 0.
!
! positive case is the default value
ZSAUNIM_IS(:) = XIMP
ZSAUNIN_IS(:) = XINP
!
! Compute the Rime Accretion Rate 
ZRAR(:) = 0.
ZVSMEAN(:) = 0.
WHERE (ZLBDAS(:) > 0. .AND. ZRCT(:) > 0.)
  ZVSMEAN(:) = XVSCOEF * ZRHOCOR(:) * ZLBDAS(:)**(-ZDS)
  ZRAR(:)    = 0.8 * ZRHODREF(:) * ZRCT(:) * ZVSMEAN(:) * 1.E3
END WHERE
!
IF (CNI_CHARGING == 'TERAR') THEN
  GELEC(:,1) = GELEC(:,1) .AND. ZRAR(:) > 0.01 .AND. ZRAR(:) <= 80.
ELSE    
  GELEC(:,1) = GELEC(:,1) .AND. ZRAR(:) > 0.1 
END IF
GELEC(:,4) = GELEC(:,1)
!
IF (COUNT(GELEC(:,4)) .GT. 0) THEN
! compute the coefficients for I-S collisions
  IF (CNI_CHARGING == 'SAP98') THEN
    CALL ELEC_INI_NI_SAP98 (KMICRO, GELEC(:,4), ZRAR, ZRAR_CRIT, ZDQ_IS)
    !
  ELSE IF (CNI_CHARGING == 'BSMP1' .OR. CNI_CHARGING == 'BSMP2') THEN
    ZRAR(:) = ZRAR(:) / 3
    CALL INTERP_DQ_TABLE (KMICRO, NIND_TEMP, NIND_LWC, GELEC(:,4), &
                          ZRAR, ZZT, XSAUNDER, ZDQ_IS)
    !
  ELSE IF (CNI_CHARGING == 'TERAR') THEN
    ZRAR(:) = ZRAR(:) / 8.
    CALL INTERP_DQ_TABLE (KMICRO, NIND_TEMP, NIND_LWC, GELEC(:,4), &
                          ZRAR, ZZT, XTAKA_TM, ZDQ_IS)
    !
  END IF
  !
  WHERE (ZDQ_IS(:) < 0.)
    ZSAUNIM_IS(:) = XIMN
    ZSAUNIN_IS(:) = XINN
  ENDWHERE
ENDIF
!
!
!*      4.      INITIALIZATION FOR GRAUPEL - SNOW COLLISIONS
!               --------------------------------------------
!
ZDQ_SG(:)  = 0.
!
! positive case is the default value
ZSAUNSK_SG(:) = XSKP
ZSAUNSM_SG(:) = XSMP
ZSAUNSN_SG(:) = XSNP
!
! Compute the Rime Accretion Rate 
ZRAR(:) = 0.
WHERE (ZVSMEAN(:) > 0. .AND. ZVGMEAN(:) > 0.)
  ZRAR(:) = 0.8 * ZRHODREF(:) * ZRCT(:) * ABS(ZVGMEAN(:) - ZVSMEAN(:)) * 1.E3
END WHERE
!
IF (CNI_CHARGING == 'TERAR') THEN
  GELEC(:,3) = GELEC(:,3) .AND. ZRAR(:) > 0.01 .AND. ZRAR(:) <= 80.
ELSE    
  GELEC(:,3) = GELEC(:,3) .AND. ZRAR(:) > 0.1 
END IF
GELEC(:,4) = GELEC(:,3)
!
IF( COUNT(GELEC(:,4)) .GT. 0) THEN
! compute the coefficients for S-G collisions
  IF (CNI_CHARGING == 'SAP98') THEN
    CALL ELEC_INI_NI_SAP98 (KMICRO, GELEC(:,4), ZRAR, ZRAR_CRIT, ZDQ_SG)
    !
  ELSE IF (CNI_CHARGING == 'BSMP1' .OR. CNI_CHARGING == 'BSMP2') THEN
    ZRAR(:) = ZRAR(:) / 3
    CALL INTERP_DQ_TABLE (KMICRO, NIND_TEMP, NIND_LWC, GELEC(:,4), &
                          ZRAR, ZZT, XSAUNDER, ZDQ_SG)
    !
  ELSE IF (CNI_CHARGING == 'TERAR') THEN
    ZRAR(:) = ZRAR(:) / 8.
    CALL INTERP_DQ_TABLE (KMICRO, NIND_TEMP, NIND_LWC, GELEC(:,4), &
                          ZRAR, ZZT, XTAKA_TM, ZDQ_SG)
    !
  END IF
  !
  WHERE (ZDQ_SG(:) < 0.)
    ZSAUNSK_SG(:) = XSKN
    ZSAUNSM_SG(:) = XSMN
    ZSAUNSN_SG(:) = XSNN
  ENDWHERE
ENDIF
!
END SUBROUTINE ELEC_INIT_NOIND_RAR
!
!------------------------------------------------------------------
!
! ##################################
  SUBROUTINE ELEC_INIT_NOIND_TAKAH()
! ################################## 
!
! Purpose : initialization for the non-inductive charging process
!           following Takahashi (1978)
!
!
!*	0.	DECLARATIONS
!		------------
!
IMPLICIT NONE
!
!
!*      1.      COMPUTE f(T, LWC)
!               -----------------
!
ZDQ(:) = 0.
!
ZEW(:) = ZRCT(:) * ZRHODREF(:) * 1.E3      ! (g m^-3)
!
GELEC(:,3) = ZZT(:) > (CST%XTT - 40.) .AND. ZZT(:) <= CST%XTT .AND. &
             ZEW(:) >= 0.01 .AND. ZEW(:) <= 10.
GELEC(:,1) = GELEC(:,3) .AND.                                              &
             ZRIT(:) > XRTMIN_ELEC(4) .AND. ZRST(:) > XRTMIN_ELEC(5) .AND. &
             ZCIT(:) > 0.0 .AND. ZLBDAS(:) > 0. .AND. ZLBDAS(:) < XLBDAS_MAXE
GELEC(:,2) = GELEC(:,3) .AND.                                              &
             ZRIT(:) > XRTMIN_ELEC(4) .AND. ZRGT(:) > XRTMIN_ELEC(6) .AND. &
             ZCIT(:) > 0.0 .AND. ZLBDAG(:) > 0. .AND. ZLBDAG(:) < XLBDAG_MAXE
GELEC(:,3) = GELEC(:,3) .AND.                                              &
             ZRST(:) > XRTMIN_ELEC(5) .AND. ZRGT(:) > XRTMIN_ELEC(6) .AND. &
             ZLBDAS(:) > 0.           .AND. ZLBDAS(:) < XLBDAS_MAXE  .AND. &
             ZLBDAG(:) > 0.           .AND. ZLBDAG(:) < XLBDAG_MAXE
GELEC(:,4) = GELEC(:,1) .OR. GELEC(:,2) .OR. GELEC(:,3)
!
IF (COUNT(GELEC(:,4)) > 0) THEN
  CALL INTERP_DQ_TABLE (KMICRO, NIND_TEMP, NIND_LWC, GELEC(:,4), &
                        ZEW, ZZT, XMANSELL, ZDQ)        
ENDIF
!
END SUBROUTINE ELEC_INIT_NOIND_TAKAH
!
!------------------------------------------------------------------
!
! ##################################
  SUBROUTINE ELEC_INIT_NOIND_TEEWC()
! ##################################
!
!
! Purpose : initialization for the non-inductive charging process
!           following Tsenova and Mitzeva (2009)
!
!
!*	0.	DECLARATIONS
!		------------
!
IMPLICIT NONE
!
!
!*      1.      PARAMETERS FOR POSITIVE NI CHARGING
!               -----------------------------------
!
GELEC(:,:) = .FALSE.
ZDQ(:) = 0.
!
! positive case is the default value
ZFQIAGGS(:)   = XFQIAGGSP_TAK
ZFQIDRYGBS(:) = XFQIDRYGBSP_TAK
ZLBQSDRYGB1S(:) = XLBQSDRYGB1SP
ZLBQSDRYGB2S(:) = XLBQSDRYGB2SP
ZLBQSDRYGB3S(:) = XLBQSDRYGB3SP
ZSAUNIM(:) = XIMP     !3.76
ZSAUNIN(:) = XINP     !2.5
ZSAUNSK(:) = XSKP_TAK !6.5
ZSAUNSM(:) = XSMP     !0.44
ZSAUNSN(:) = XSNP     !2.5
!
!
!*      2.      PARAMETERS FOR NEGATIVE NI CHARGING
!               -----------------------------------   
!
! Compute the effective water content
ZEW(:) = 0.
!
WHERE (ZLBDAG(:) > 0. .AND. ZRCT(:) > 0.)
  ZEW(:) = 0.8 * ZRHODREF(:) * ZRCT(:) * 1.E3 
END WHERE
!
GELEC(:,3) = ZZT(:) >= (CST%XTT - 40.) .AND. ZZT(:) <= CST%XTT .AND.  &
             ZEW(:) >= 0.01 .AND. ZEW(:) <= 10.
GELEC(:,1) = GELEC(:,3) .AND.                                              &
             ZRIT(:) > XRTMIN_ELEC(4) .AND. ZRST(:) > XRTMIN_ELEC(5) .AND. &
             ZCIT(:) > 0.0 .AND. ZLBDAS(:) > 0. .AND. ZLBDAS(:) < XLBDAS_MAXE
GELEC(:,2) = GELEC(:,3) .AND.                                              &
             ZRIT(:) > XRTMIN_ELEC(4) .AND. ZRGT(:) > XRTMIN_ELEC(6) .AND. &
             ZCIT(:) > 0.0 .AND. ZLBDAG(:) > 0. .AND. ZLBDAG(:) < XLBDAG_MAXE
GELEC(:,3) = GELEC(:,3) .AND.                                              &
             ZRST(:) > XRTMIN_ELEC(5) .AND. ZRGT(:) > XRTMIN_ELEC(6) .AND. &
             ZLBDAS(:) > 0.           .AND. ZLBDAS(:) < XLBDAS_MAXE  .AND. &
             ZLBDAG(:) > 0.           .AND. ZLBDAG(:) < XLBDAG_MAXE
!
GELEC(:,4) = GELEC(:,1) .OR. GELEC(:,2) .OR. GELEC(:,3)
!
IF (COUNT(GELEC(:,4)) > 0) THEN
  CALL INTERP_DQ_TABLE (KMICRO, NIND_TEMP, NIND_LWC, GELEC(:,4), &
                        ZEW, ZZT, XTAKA_TM, ZDQ)        
!
  WHERE (ZDQ(:) < 0.)
    ZFQIAGGS(:)   = XFQIAGGSN_TAK
    ZFQIDRYGBS(:) = XFQIDRYGBSN_TAK
    ZLBQSDRYGB1S(:) = XLBQSDRYGB1SN
    ZLBQSDRYGB2S(:) = XLBQSDRYGB2SN
    ZLBQSDRYGB3S(:) = XLBQSDRYGB3SN
    ZSAUNIM(:) = XIMN      !2.54
    ZSAUNIN(:) = XINN      !2.8
    ZSAUNSK(:) = XSKN_TAK  !2.0
    ZSAUNSM(:) = XSMN      !0.5
    ZSAUNSN(:) = XSNN      !2.8
 ENDWHERE
ENDIF
!
END SUBROUTINE ELEC_INIT_NOIND_TEEWC
!
!------------------------------------------------------------------
!
! #################################################################
  SUBROUTINE ELEC_INI_NI_SAP98(KMICRO, OMASK, PRAR, PRAR_CRIT, PDQ)
! #################################################################
!
!
! Purpose : compute dQ(RAR,T) in the parameterization of Saunders and Peck (1998)
!
!
!*	0.	DECLARATIONS
!		------------
!
IMPLICIT NONE
!
INTEGER,                    INTENT(IN)    :: KMICRO
LOGICAL, DIMENSION(KMICRO), INTENT(IN)    :: OMASK
REAL,    DIMENSION(KMICRO), INTENT(IN)    :: PRAR       ! Rime accretion rate 
REAL,    DIMENSION(KMICRO), INTENT(IN)    :: PRAR_CRIT  ! Critical rime accretion rate
REAL,    DIMENSION(KMICRO), INTENT(INOUT) :: PDQ        ! interpolated dQ
!
!
!*      1.      COMPUTE dQ(RAR, T)
!		------------------
!
PDQ(:) = 0.
!
! positive region : Mansell et al., 2005
WHERE (OMASK(:) .AND. PRAR(:) > PRAR_CRIT(:))
  PDQ(:) = MAX(0., 6.74 * (PRAR(:) - PRAR_CRIT(:)) * 1.E-15)
ENDWHERE
!
! negative region : Mansell et al. 2005
WHERE (OMASK(:) .AND. PRAR(:) < PRAR_CRIT(:))
  PDQ(:) = MIN(0., 3.9 * (PRAR_CRIT(:) - 0.1) *                   &
                  (4.0 * ((PRAR(:) - (PRAR_CRIT(:) + 0.1) / 2.) / &
                         (PRAR_CRIT(:) - 0.1))**2 - 1.) * 1.E-15)
ENDWHERE
!
END SUBROUTINE ELEC_INI_NI_SAP98
!
!------------------------------------------------------------------
!
! #################################################################
  SUBROUTINE INTERP_DQ_TABLE (KMICRO, KIND_TEMP, KIND_LWC, OMASK, &
                              PLIQ, PTEMP, PTABLE, PDQ)
! ################################################################# 
!
!
! Purpose : interpolate dQ from a lookup table at each gridpoint
!
!
!*	0.	DECLARATIONS
!		------------
!
IMPLICIT NONE
!
INTEGER,                    INTENT(IN)    :: KMICRO
INTEGER,                    INTENT(IN)    :: KIND_TEMP, KIND_LWC
LOGICAL, DIMENSION(KMICRO), INTENT(IN)    :: OMASK
REAL,    DIMENSION(KMICRO), INTENT(IN)    :: PLIQ   ! effective water content or rime accretion rate
REAL,    DIMENSION(KMICRO), INTENT(IN)    :: PTEMP  ! temperature
REAL,    DIMENSION(KIND_LWC+1,KIND_TEMP+1), INTENT(IN)    :: PTABLE  ! lookup table
REAL,    DIMENSION(KMICRO), INTENT(INOUT) :: PDQ    ! interpolated dQ
!
! declaration of local variables
INTEGER                 :: IGAUX
REAL,    DIMENSION(:), ALLOCATABLE :: ZDQ_INTERP
REAL,    DIMENSION(:), ALLOCATABLE :: ZVECT1, ZVECT2
INTEGER, DIMENSION(:), ALLOCATABLE :: IVECT1, IVECT2
!
!
!*	1.      FIND THE INDEXES FOR RAR/EW AND T
!		---------------------------------
!
PDQ(:) = 0.
!
IGAUX = 0
DO II = 1, KMICRO
  IF (OMASK(II)) THEN
    IGAUX = IGAUX + 1
    I1(IGAUX) = II
  END IF 
END DO
!
IF (IGAUX > 0) THEN
  ALLOCATE(ZDQ_INTERP(IGAUX))
  ALLOCATE(ZVECT1(IGAUX))
  ALLOCATE(ZVECT2(IGAUX))
  ALLOCATE(IVECT1(IGAUX))
  ALLOCATE(IVECT2(IGAUX))
  ZDQ_INTERP(:) = 0.
  IVECT1(:) = 0
  IVECT2(:) = 0
!
  DO II = 1, IGAUX
    ZVECT1(II) = PTEMP(I1(II))
    ZVECT2(II) = PLIQ(I1(II))
    ZDQ_INTERP(II) = PDQ(I1(II))
  END DO
!
! Temperature index (0C --> -40C)
  ZVECT1(1:IGAUX) = MAX( 1.00001, MIN( REAL(KIND_TEMP)-0.00001, &
                                    (ZVECT1(1:IGAUX) - CST%XTT - 1.)/(-1.) ) )
  IVECT1(1:IGAUX) = INT( ZVECT1(1:IGAUX) )
  ZVECT1(1:IGAUX) = ZVECT1(1:IGAUX) - REAL(IVECT1(1:IGAUX))
!
! LWC index (0.01 g.m^-3 --> 10 g.m^-3)
  WHERE (ZVECT2(:) >= 0.01 .AND. ZVECT2(:) < 0.1)
    ZVECT2(:) = MAX( 1.00001, MIN( REAL(10)-0.00001, &
                                  ZVECT2(:) * 100. ))
    IVECT2(:) = INT(ZVECT2(:))
    ZVECT2(:) = ZVECT2(:) - REAL(IVECT2(:))
  ENDWHERE
!
  WHERE (ZVECT2(:) >= 0.1 .AND. ZVECT2(:) < 1. .AND. IVECT2(:) == 0)
    ZVECT2(:) = MAX( 10.00001, MIN( REAL(19)-0.00001, &
                                   ZVECT2(:) * 10. + 9. ) )
    IVECT2(:) = INT(ZVECT2(:))
    ZVECT2(:) = ZVECT2(:) - REAL(IVECT2(:))
  ENDWHERE
!
  WHERE ((ZVECT2(:) >= 1.) .AND. ZVECT2(:) <= 10.)
    ZVECT2(:) = MAX( 19.00001, MIN( REAL(KIND_LWC)-0.00001, &
                                   ZVECT2(:) + 18. ) )
    IVECT2(:) = INT(ZVECT2(:))
    ZVECT2(:) = ZVECT2(:) - REAL(IVECT2(:))
  ENDWHERE
!
!
!*	2.      INTERPOLATE dQ(RAR or EW,T)
!		---------------------------
!
  ZDQ_INTERP(:) = BI_LIN_INTP_V( PTABLE, IVECT2, IVECT1, ZVECT2, ZVECT1, &
                                 IGAUX )
!
  DO II = 1, IGAUX
    PDQ(I1(II)) = ZDQ_INTERP(II)
  END DO
END IF
!
DEALLOCATE(ZDQ_INTERP)
DEALLOCATE(ZVECT1)
DEALLOCATE(ZVECT2)
DEALLOCATE(IVECT1)
DEALLOCATE(IVECT2)
!
END SUBROUTINE INTERP_DQ_TABLE
!
!------------------------------------------------------------------
!
! #########################
  SUBROUTINE ELEC_IAGGS_B()
! #########################
!
!
! Purpose : compute charge separation process during the collision
!           between ice and snow
!
!
!*      0.     DECLARATIONS
!              ------------
!
IMPLICIT NONE
!
!
!*      1.     COMPUTE THE COLLISION EFFICIENCY
!              --------------------------------
!
ZQCOLIS(:) = ZCOLIS * EXP(ZCOLEXIS * (ZZT(:) - CST%XTT))
!
ZWQ_NI(:) = 0.
ZLIMIT(:) = 0.
!
!*      2.     COMPUTE THE RATE OF SEPARATED CHARGE
!              ------------------------------------
!
!*      2.1    Charging process following Helsdon and Farley (1987)
!
IF (CNI_CHARGING == 'HELFA') THEN
  !
  WHERE (ZCIT(:) > 0.0 .AND. &
         ZRIT(:) > XRTMIN_ELEC(4) .AND. ZRST(:) > XRTMIN_ELEC(5))
    ZWQ_NI(:) = XFQIAGGSBH * ZRIAGGS(:) * ZCIT(:) / ZRIT(:)
    ZWQ_NI(:) = ZWQ_NI(:) * (1. - ZQCOLIS(:)) / ZQCOLIS(:) 
!
! Temperature dependance of the charge transferred
    ZWQ_NI(:) = ZWQ_NI(:) * (ZZT(:) - XQTC) / ABS(ZZT(:) - XQTC)
    ZWQ_NI(:) = ZWQ_NI(:) / ZRHODREF(:)
  END WHERE
!
ELSE
!
!
!*      2.2    Charging process following Gardiner et al. (1985)
!
  IF (CNI_CHARGING == 'GARDI') THEN
    WHERE (GELEC(:,1) .AND. ZDELTALWC(:) > 0.)
      ZWQ_NI(:) = XFQIAGGSBG  * (1 - ZQCOLIS(:)) *         &
                  ZRHODREF(:)**(-4. * ZCEXVT + 4. / ZBI) * &
                  ZCIT(:)**(1 - 4. / ZBI) *                &
                  ZDELTALWC(:) * ZFT(:) *                  &
                  ZCST(:) * ZLBDAS(:)**(-2. - 4. * ZDS) *  &
                  (ZAI * MOMG(ZALPHAI, ZNUI, ZBI) /        &
                  ZRIT(:))**(-4 / ZBI)
    ENDWHERE
  END IF
!
!
!*      2.3    Charging process based on EW: SAUN1/SAUN2, TEEWC
!*             following Saunders et al. (1991), Takahashi via Tsenova and Mitzeva (2009)
!
  IF (CNI_CHARGING == 'SAUN1' .OR. CNI_CHARGING == 'SAUN2' .OR. &
      CNI_CHARGING == 'TEEWC') THEN
    WHERE (GELEC(:,1) .AND. ZDQ(:) /= 0.)
      ZWQ_NI(:) = XFQIAGGSBS * (1 - ZQCOLIS(:)) *                       &
                  ZRHOCOR(:)**(1 + ZSAUNIN(:)) *                        &
                  ZFQIAGGS(:) * ZDQ(:) *                                &
                  ZCIT(:)**(1 - ZSAUNIM(:) / ZBI) *                     &
                  ZCST(:) * ZLBDAS(:)**(-2.- ZDS * (1. + ZSAUNIN(:))) * &
                 (ZRHODREF(:) * ZRIT(:) / XAIGAMMABI)**(ZSAUNIM(:) / ZBI)
    ENDWHERE
  END IF
!
!
!*      2.4    Charging process based on RAR (=EW*V): SAP98, BSMP1/BSMP2, TERAR
!*             following Saunders and Peck (1998) or Brooks et al., 1997 (with/out anomalies) 
!*             or Takahashi via Tsenova and Mitzeva (2011)
!
  IF (CNI_CHARGING == 'SAP98' .OR. CNI_CHARGING == 'BSMP1' .OR.   &
      CNI_CHARGING == 'BSMP2' .OR. CNI_CHARGING == 'TERAR') THEN
    IF (CNI_CHARGING /= 'TERAR') THEN
      ZFQIAGGS(:) = XFQIAGGSP
      WHERE (ZDQ_IS(:) < 0.)
        ZFQIAGGS(:) = XFQIAGGSN
      ENDWHERE
    ELSE
      ZFQIAGGS(:) = XFQIAGGSP_TAK
      WHERE (ZDQ_IS(:) < 0.)
        ZFQIAGGS(:) = XFQIAGGSN_TAK
      ENDWHERE
    ENDIF
!
    WHERE (GELEC(:,1) .AND. ZDQ_IS(:) /= 0.)
      ZWQ_NI(:) = XFQIAGGSBS * (1 - ZQCOLIS(:)) *                          &
                  ZRHOCOR(:)**(1 + ZSAUNIN_IS(:)) *                        &
                  ZFQIAGGS(:) * ZDQ_IS(:) *                                &
                  ZCIT(:)**(1 - ZSAUNIM_IS(:) / ZBI) *                     &
                  ZCST(:) * ZLBDAS(:)**(-2.- ZDS * (1. + ZSAUNIN_IS(:))) * &
                 (ZRHODREF(:) * ZRIT(:) / XAIGAMMABI)**(ZSAUNIM_IS(:) / ZBI)
    ENDWHERE
  END IF
!
!
!*      2.5    Charging process following Takahashi (1978)
!
  IF (CNI_CHARGING == 'TAKAH') THEN
    WHERE (GELEC(:,1) .AND. ZDQ(:) /= 0.)
      ZWQ_NI(:) = XFQIAGGSBT1 * (1.0 - ZQCOLIS(:)) * ZRHOCOR(:) *         &
                  ZCIT(:) * ZCST(:) * ZDQ(:) *                            &
                  MIN( XFQIAGGSBT2 / (ZLBDAS(:)**(2. + ZDS)) ,            &
                       XFQIAGGSBT3 * ZRHOCOR(:) * ZRHODREF(:)**(2./ZBI) * &
                       ZRIT(:)**(2. / ZBI) /                              &
                      (ZCIT(:)**(2. / ZBI) * ZLBDAS(:)**(2. + 2. * ZDS)))
    ENDWHERE
  END IF
!
!
!*      3.     LIMITATION OF THE SEPARATED CHARGE
!              ----------------------------------
!
! Dq is limited to XLIM_NI_IS
  WHERE (ZWQ_NI(:) .NE. 0.)
    ZLIMIT(:) = XLIM_NI_IS * ZRIAGGS(:) * ZCIT(:) * &
                (1 - ZQCOLIS(:)) / (ZRIT(:) * ZQCOLIS(:))
    ZWQ_NI(:) = SIGN( MIN( ABS(ZLIMIT(:)), ABS(ZWQ_NI(:) ) ), ZWQ_NI(:) )   
    ZWQ_NI(:) = ZWQ_NI(:) / ZRHODREF(:)
  ENDWHERE
!
! For temperatures lower than -30C --> linear interpolation
  WHERE (ZWQ_NI(:) /= 0. .AND. ZZT(:) < (CST%XTT-30.) .AND. ZZT(:) >= (CST%XTT-40.))
    ZWQ_NI(:) = ZWQ_NI(:) * (ZZT(:) - CST%XTT + 40.) / 10.
  ENDWHERE
!
END IF
!
ZQSS(:) = ZQSS(:) + ZWQ_NI(:)
ZQIS(:) = ZQIS(:) - ZWQ_NI(:)
!
END SUBROUTINE ELEC_IAGGS_B
!
!------------------------------------------------------------------
!
! #########################
  SUBROUTINE ELEC_IDRYG_B()
! #########################
!
!
!   Purpose : compute charge separation process during the dry collision
!             between ice and graupeln
!
!
!*      0.     DECLARATIONS
!              ------------
!
IMPLICIT NONE
!
!
!*      1.     COMPUTE THE COLLISION EFFICIENCY
!              --------------------------------
!
ZQCOLIG(:) = ZCOLIG * EXP(ZCOLEXIG * (ZZT(:) - CST%XTT))
!
ZWQ_NI(:) = 0.
ZLIMIT(:) = 0.
!
!*      2.     COMPUTE THE RATE OF SEPARATED CHARGE
!              ------------------------------------
!
!*      2.1    Charging process following Helsdon and Farley (1987)
!
IF (CNI_CHARGING == 'HELFA') THEN
  !
  WHERE (ZRIT(:) > XRTMIN_ELEC(4) .AND. ZCIT(:) > 0. .AND. &
         ZRGT(:) > XRTMIN_ELEC(6))
    ZWQ_NI(:) = XHIDRYG * ZRIDRYG(:) * ZCIT(:) / ZRIT(:)    
    ZWQ_NI(:) = ZWQ_NI(:) * (1. - ZQCOLIG(:)) / ZQCOLIG(:)        ! QIDRYG_boun
!
! Temperature dependance of the charge transfered
    ZWQ_NI(:) = ZWQ_NI(:) * (ZZT(:) - XQTC) / ABS(ZZT(:) - XQTC)
    ZWQ_NI(:) = ZWQ_NI(:) / ZRHODREF(:)
  END WHERE
!
ELSE
!
!
!*      2.2    Charging process following Gardiner et al. (1985)
!
  IF (CNI_CHARGING == 'GARDI') THEN
    WHERE (GELEC(:,2) .AND. ZDELTALWC(:) > 0.) 
      ZWQ_NI(:) = XFQIDRYGBG * XLBQIDRYGBG * (1 - ZQCOLIG) * &
                  ZRHODREF(:)**(-4. * ZCEXVT + 4. / ZBI) *   &
                  ZCIT(:)**(1 - 4. / ZBI) *                  &
                  ZDELTALWC(:) * ZFT(:) *                    &
                  ZCGT(:) * ZLBDAG(:)**(-2. - 4. * ZDG) *    &
                  (ZAI * MOMG(ZALPHAI, ZNUI, ZBI) /          &
                  ZRIT(:))**(-4 / ZBI)
    ENDWHERE
  END IF
!
!
!*      2.3    Charging process based on EW: SAUN1/SAUN2, TEEWC following
!*             Saunders et al. (1991), Takahashi via Tsenova and Mitzeva(2009)
!
  IF (CNI_CHARGING == 'SAUN1' .OR. CNI_CHARGING == 'SAUN2' .OR. &
      CNI_CHARGING == 'TEEWC') THEN
    WHERE (GELEC(:,2) .AND. ZDQ(:) /= 0.)
      ZWQ_NI(:) = XFQIDRYGBS * (1. - ZQCOLIG(:)) *                       &
                  ZRHOCOR(:)**(1. + ZSAUNIN(:)) *                        &
                  ZFQIDRYGBS(:) * ZDQ(:) *                               &
                  ZCIT(:)**(1. - ZSAUNIM(:) / ZBI) *                     &
                  ZCGT(:) * ZLBDAG(:)**(-2. - ZDG * (1. + ZSAUNIN(:))) * &
                 (ZRHODREF(:) * ZRIT(:) / XAIGAMMABI)**(ZSAUNIM(:) / ZBI)
    ENDWHERE
  END IF
!
!
!*      2.4    Charging process based on RAR (=EW*V): SAP98, BSMP1/BSMP2, TERAR
!*             following Saunders and Peck (1998) or Brooks et al., 1997 (with/out anomalies)
!*             or Takahashi via Tsenova and Mitzeva (2011)
!
  IF (CNI_CHARGING == 'SAP98' .OR. CNI_CHARGING == 'BSMP1' .OR.   &
      CNI_CHARGING == 'BSMP2' .OR. CNI_CHARGING == 'TERAR') THEN
    IF (CNI_CHARGING /= 'TERAR') THEN
       ZFQIDRYGBS(:) = XFQIDRYGBSP
       WHERE (ZDQ_IG(:) < 0.)
         ZFQIDRYGBS(:) = XFQIDRYGBSN
       ENDWHERE
    ELSE
       ZFQIDRYGBS(:) = XFQIDRYGBSP_TAK
       WHERE (ZDQ_IG(:) <0.)
         ZFQIDRYGBS(:) = XFQIDRYGBSN_TAK
       ENDWHERE
    END IF
!
    WHERE (GELEC(:,2) .AND. ZDQ_IG(:) /= 0.)
      ZWQ_NI(:) = XFQIDRYGBS * (1. - ZQCOLIG(:)) *                          &
                  ZRHOCOR(:)**(1 + ZSAUNIN_IG(:)) *                         &
                  ZFQIDRYGBS(:) * ZDQ_IG(:) *                               &
                  ZCIT(:)**(1 - ZSAUNIM_IG(:) / ZBI) *                      &
                  ZCGT(:) * ZLBDAG(:)**(-2. - ZDG * (1. + ZSAUNIN_IG(:))) * &
                 (ZRHODREF(:) * ZRIT(:) / XAIGAMMABI)**(ZSAUNIM_IG(:) / ZBI)
    ENDWHERE
  END IF
!
!
!*      2.5     Charging process following Takahashi (1978)
!
  IF (CNI_CHARGING == 'TAKAH') THEN
    WHERE (GELEC(:,2) .AND. ZDQ(:) /= 0.)
      ZWQ_NI(:) = XFQIDRYGBT1 * (1. - ZQCOLIG(:)) * ZRHOCOR(:) *          &
                  ZCIT(:) * ZCGT(:) * ZDQ(:) *                            &
                  MIN( XFQIDRYGBT2 / (ZLBDAG(:)**(2. + ZDG)),             &
                       XFQIDRYGBT3 * ZRHOCOR(:) * ZRHODREF(:)**(2./ZBI) * &
                       ZRIT(:)**(2. / ZBI) / (ZCIT(:)**(2. / ZBI) *       &
                       ZLBDAG(:)**(2. + 2. * ZDG)) )
    ENDWHERE
  END IF
!
!
!*      3.     LIMITATION OF THE SEPARATED CHARGE
!              ----------------------------------
!
! Dq is limited to XLIM_NI_IG
  WHERE (ZWQ_NI(:) .NE. 0. .AND. ZRIT(:) > 0.)
    ZLIMIT(:) = XLIM_NI_IG * ZRIDRYG(:) * ZCIT(:) * (1 - ZQCOLIG(:)) / &
                (ZRIT(:) * ZQCOLIG(:))
    ZWQ_NI(:) = SIGN( MIN( ABS(ZLIMIT(:)), ABS(ZWQ_NI(:) ) ), ZWQ_NI(:) )   
    ZWQ_NI(:) = ZWQ_NI(:) / ZRHODREF(:)
  ENDWHERE
!
! For temperatures lower than -30C --> linear interpolation
  WHERE (ZWQ_NI(:) /= 0. .AND. ZZT(:) < (CST%XTT-30.) .AND. ZZT(:) >= (CST%XTT-40.))
    ZWQ_NI(:) = ZWQ_NI(:) * (ZZT(:) - CST%XTT + 40.) / 10.
  ENDWHERE
!
END IF
!
WHERE (ZRIDRYG(:) > 0.)
  ZQGS(:) = ZQGS(:) + ZWQ_NI(:)
  ZQIS(:) = ZQIS(:) - ZWQ_NI(:)
END WHERE
!
END SUBROUTINE ELEC_IDRYG_B
!
!------------------------------------------------------------------
!
! #########################
  SUBROUTINE ELEC_SDRYG_B()
! #########################
!
!
! Purpose : compute the charge separation during the dry collision
!           between snow and graupel
!
!
!*	0.	DECLARATIONS
!               ------------
!
IMPLICIT NONE
!
!
!*      1.      COMPUTE THE COLLECTION EFFICIENCY
!               ---------------------------------
!
ZQCOLSG(:) = ZCOLSG * EXP (ZCOLEXSG * (ZZT(:) - CST%XTT))
!
ZWQ_NI(:) = 0.
ZLIMIT(:) = 0.
!
!*      2.      COMPUTE THE RATE OF SEPARATED CHARGE
!               ------------------------------------
!
!*      2.1     Charge separation following Helsdon and Farley (1987)
! 
IF (CNI_CHARGING == 'HELFA') THEN
!
  WHERE (ZRGT(:) > XRTMIN_ELEC(6) .AND. ZLBDAG(:) > 0. .AND. &
         ZRST(:) > XRTMIN_ELEC(5) .AND. ZLBDAS(:) > 0.)
    ZWQ_NI(:) = ZWQ5(:,5) * XFQSDRYGBH * ZRHODREF(:)**(-ZCEXVT) *    &
               (1. - ZQCOLSG(:)) *                                   &
                ZCST(:) * ZCGT(:) *                                  &
               (XLBQSDRYGB4H * ZLBDAS(:)**(-2.) +                    &
                XLBQSDRYGB5H * ZLBDAS(:)**(-1.) * ZLBDAG(:)**(-1.) + &
                XLBQSDRYGB6H * ZLBDAG(:)**(-2.))
!
! Temperature dependance of the charge transfered
    ZWQ_NI(:) = ZWQ_NI(:) * (ZZT(:) - XQTC) / ABS(ZZT(:) - XQTC)
    ZWQ_NI(:) = ZWQ_NI(:) / ZRHODREF(:)
  ENDWHERE
!
ELSE
!
!
!*      2.2     Charge separation following Gardiner et al. (1985)
!
  IF (CNI_CHARGING == 'GARDI') THEN 
    WHERE (GELEC(:,3) .AND. ZDELTALWC(:) > 0.) 
      ZWQ_NI(:) = XFQSDRYGBG * (1. - ZQCOLSG(:)) *                     &
                  ZRHODREF(:)**(-4. * ZCEXVT) *                        &
                  ZFT(:) * ZDELTALWC(:) *                              &
                  ZCST(:) * ZCGT(:) *                                  &
                 (XLBQSDRYGB4G * ZLBDAS(:)**(-4.) * ZLBDAG(:)**(-2.) + &
                  XLBQSDRYGB5G * ZLBDAS(:)**(-5.) * ZLBDAG(:)**(-1.) + &
                  XLBQSDRYGB6G * ZLBDAS(:)**(-6.)) *                   &
                  ZWQ5(:,5)
    ENDWHERE
  END IF
!
!
!*      2.3     Charging process based on EW: SAUN1/SAUN2, TEEWC following
!*              Saunders et al. (1991), Takahashi via Tsenova and Mitzeva(2009)
!
  IF (CNI_CHARGING == 'SAUN1' .OR. CNI_CHARGING == 'SAUN2' .OR. &
      CNI_CHARGING == 'TEEWC') THEN
    WHERE (GELEC(:,3) .AND. ZDQ(:) /= 0.) 
!    ZWQ_NI(:) = ZWQ5(:,6) If graupel gains positive charge ZDQ(:) > 0.
!    ZWQ_NI(:) = ZWQ5(:,7) If graupel gains negative charge ZDQ(:) < 0.
      ZWQ_NI(:) = ZWQ5(:,6) * (0.5 + SIGN(0.5,ZDQ(:))) + &
                  ZWQ5(:,7) * (0.5 - SIGN(0.5,ZDQ(:)))
!
      ZWQ_NI(:) = ZWQ_NI(:) * XFQSDRYGBS * (1. - ZQCOLSG(:)) *                  &
                  ZRHOCOR(:)**(1. + ZSAUNSN(:)) * ZSAUNSK(:) * ZDQ(:) *         &
                  ZCST(:) * ZCGT(:) *                                           &
                ( ZLBQSDRYGB1S(:) / (ZLBDAS(:)**ZSAUNSM(:) * ZLBDAG(:)**2) +    &
                  ZLBQSDRYGB2S(:) / (ZLBDAS(:)**( 1.+ZSAUNSM(:)) * ZLBDAG(:)) + &
                  ZLBQSDRYGB3S(:) /  ZLBDAS(:)**(2.+ZSAUNSM(:)) )
    ENDWHERE
  END IF
!
!
!*      2.4     Charging process based on RAR (=EW*V): SAP98, BSMP1/BSMP2, TERAR
!*              following Saunders and Peck (1998) or Brooks et al., 1997 (with/out anomalies)
!*              or Takahashi via Tsenova and Mitzeva (2011)
!
  IF (CNI_CHARGING == 'SAP98' .OR. CNI_CHARGING == 'BSMP1' .OR.   &
      CNI_CHARGING == 'BSMP2' .OR. CNI_CHARGING == 'TERAR') THEN
    ZLBQSDRYGB1S(:) = XLBQSDRYGB1SP
    ZLBQSDRYGB2S(:) = XLBQSDRYGB2SP
    ZLBQSDRYGB3S(:) = XLBQSDRYGB3SP
    WHERE (ZDQ_SG(:) < 0.)
      ZLBQSDRYGB1S(:) = XLBQSDRYGB1SN
      ZLBQSDRYGB2S(:) = XLBQSDRYGB2SN
      ZLBQSDRYGB3S(:) = XLBQSDRYGB3SN
    ENDWHERE
!
    WHERE (GELEC(:,3) .AND. ZDQ_SG(:) /= 0.)
      ZWQ_NI(:) = ZWQ5(:,6) * (0.5+SIGN(0.5,ZDQ_SG(:))) + &
                  ZWQ5(:,7) * (0.5-SIGN(0.5,ZDQ_SG(:)))
!
      ZWQ_NI(:) = ZWQ_NI(:) * XFQSDRYGBS * (1. - ZQCOLSG(:)) *                       &
                  ZRHOCOR(:)**(1. + ZSAUNSN_SG(:)) * ZSAUNSK_SG(:) * ZDQ_SG(:) *     &
                  ZCST(:) * ZCGT(:) *                                                &
                 (ZLBQSDRYGB1S(:) / (ZLBDAS(:)**ZSAUNSM_SG(:)      * ZLBDAG(:)**2) + &
                  ZLBQSDRYGB2S(:) / (ZLBDAS(:)**(1.+ZSAUNSM_SG(:)) * ZLBDAG(:)) +    &
                  ZLBQSDRYGB3S(:) /  ZLBDAS(:)**(2.+ZSAUNSM_SG(:)) )
    ENDWHERE
  END IF
!
!
!*      2.5     Charging process following Takahashi (1978)
!
  IF (CNI_CHARGING == 'TAKAH') THEN 
    WHERE (GELEC(:,3) .AND. ZDQ(:) /= 0.)
      ZWQ_NI(:) = XFQSDRYGBT1 * (1. - ZQCOLSG(:)) * ZRHOCOR(:) *            &
                  ZCGT(:) * ZCST(:) * ZDQ(:) *                              &
                  MIN(10. * (                                               &
                   ABS(XFQSDRYGBT2 / (ZLBDAG(:)**ZDG * ZLBDAS(:)**2.) -     & 
                       XFQSDRYGBT3 / (ZLBDAS(:)**(2. + ZDS))) +             &
                   ABS(XFQSDRYGBT4 / (ZLBDAG(:)**(2. + ZDG)) -              &
                       XFQSDRYGBT5 / (ZLBDAS(:)**ZDS * ZLBDAG(:)**2.)) +    &
                   ABS(XFQSDRYGBT6 / (ZLBDAG(:)**(1. + ZDG) * ZLBDAS(:)) -  &
                       XFQSDRYGBT7 / (ZLBDAS(:)**(1. + ZDS) * ZLBDAG(:)))), &
                   XFQSDRYGBT8 * ZRHOCOR(:) * ZWQ5(:,5) *                   &
                  (XFQSDRYGBT9  / (ZLBDAS(:)**2. * ZLBDAG(:)**2.) +         &
                   XFQSDRYGBT10 / (ZLBDAS(:)**4.) +                         &
                   XFQSDRYGBT11 / (ZLBDAS(:)**3. * ZLBDAG(:))))
    ENDWHERE
  END IF
!
!
!*      3.     LIMITATION OF THE SEPARATED CHARGE
!              ----------------------------------
!
! Dq is limited to XLIM_NI_SG
  WHERE (ZWQ_NI(:) .NE. 0.)
    ZLIMIT(:) = XLIM_NI_SG * ZWQ5(:,4) * XAUX_LIM *  &
                ZRHOCOR(:) * (1. - ZQCOLSG(:)) *     &
                ZCST(:) * ZCGT(:) *                  &
              ( XAUX_LIM1 / ZLBDAS(:)**2           + &
                XAUX_LIM2 /(ZLBDAS(:) * ZLBDAG(:)) + &
                XAUX_LIM3 / ZLBDAG(:)**2 )
!          
    ZWQ_NI(:) = SIGN( MIN( ABS(ZLIMIT(:)), ABS(ZWQ_NI(:) ) ), ZWQ_NI(:) )   
    ZWQ_NI(:) = ZWQ_NI(:) / ZRHODREF(:)
  ENDWHERE
!
! For temperatures lower than -30C --> linear interpolation
  WHERE (ZWQ_NI(:) /= 0. .AND. ZZT(:) < (CST%XTT-30.) .AND. ZZT(:) >= (CST%XTT-40.))
    ZWQ_NI(:) = ZWQ_NI(:) * (ZZT(:) - CST%XTT + 40.) / 10.
  ENDWHERE
!
END IF
!
WHERE (ZRSDRYG(:) > 0.)
  ZQGS(:) = ZQGS(:) + ZWQ_NI(:)
  ZQSS(:) = ZQSS(:) - ZWQ_NI(:)
END WHERE
!
END SUBROUTINE ELEC_SDRYG_B
!
!------------------------------------------------------------------
!
! #########################################################################
  SUBROUTINE COMPUTE_CHARGE_TRANSFER (PR_RATE, PRXT, PQXT, PTSTEP,        &
                                      PRX_THRESH, PQX_THRESH, PCOEF_RQ_X, &
                                      PQ_RATE, PQXS, PQYS                 )
! #########################################################################
!
!  Purpose : compute the charge transfer rate in proportion of the mass transfer rate
!  x --> y
!  q_rate_xy = r_rate_xy * coef_rq_x * qx_t / rx_t
!
!
!*      0.      DECLARATIONS
!               ------------
!
IMPLICIT NONE
!
!*      0.1     Declaration of dummy arguments
!
REAL, INTENT(IN),    DIMENSION(:)  :: PR_RATE     ! Mass exchange rate from x to y
REAL, INTENT(IN),    DIMENSION(:)  :: PRXT        ! Mixing ratio of x at t
REAL, INTENT(IN),    DIMENSION(:)  :: PQXT        ! Electric charge of x at t
REAL, INTENT(IN)                   :: PRX_THRESH  ! Threshold on mixing ratio
REAL, INTENT(IN)                   :: PQX_THRESH  ! Threshold on electric charge
REAL, INTENT(IN)                   :: PCOEF_RQ_X  ! Coefficient for charge exchange
REAL, INTENT(IN)                   :: PTSTEP      ! Time step
REAL, INTENT(INOUT), DIMENSION(:)  :: PQ_RATE     ! Charge exchange rate from x to y
REAL, INTENT(INOUT), DIMENSION(:)  :: PQXS        ! Electric charge of x - source term
REAL, INTENT(INOUT), DIMENSION(:)  :: PQYS        ! Electric charge of y - source term
!
!
!*      0.2     Declaration of local variables
!
!
PQ_RATE(:) = 0.
!
WHERE (PR_RATE(:) > 0. .AND. &
       PRXT(:) > PRX_THRESH .AND. ABS(PQXT(:)) > PQX_THRESH)
! Compute the charge exchanged during the mass tranfer from species x to y
  PQ_RATE(:) = PCOEF_RQ_X * PR_RATE(:) * PQXT(:) / PRXT(:)
! Limit the charge exchanged to the charge available on x at t
  PQ_RATE(:) = SIGN( MIN( ABS(PQXT(:)/PTSTEP),ABS(PQ_RATE(:)) ), PQXT(:)/PTSTEP )
  !
! Update the source terms of x and y
  PQXS(:) = PQXS(:) - PQ_RATE(:)
  PQYS(:) = PQYS(:) + PQ_RATE(:)
END WHERE
!
END SUBROUTINE COMPUTE_CHARGE_TRANSFER
!
!------------------------------------------------------------------
!
! ###########################################################
  FUNCTION BI_LIN_INTP_V(ZT, KI, KJ, PDX, PDY, KN)  RESULT(PY)
! ###########################################################          
!
! Purpose : 
!
!                 |                   |
! ZT(KI(1),KJ(2))-|-------------------|-ZT(KI(2),KJ(2))
!                 |                   | 
!                 |                   |
!              x2-|-------|y(x1,x2)   |
!                 |       |           |
!              PDY|       |           |
!                 |       |           |
!                 |       |           |
!ZT( KI(1),KJ(1))-|-------------------|-ZT(KI(2),KJ(1))
!                 |  PDX  |x1         |
!                 |                   |
!
!*	0.	DECLARATIONS
!          	------------
!
IMPLICIT NONE
!
!*	0.1	Declaration of dummy arguments
!
INTEGER, INTENT(IN)                 :: KN        ! Size of the result vector
INTEGER, INTENT(IN), DIMENSION(KN)  :: KI        ! Tabulated  coordinate
INTEGER, INTENT(IN), DIMENSION(KN)  :: KJ        ! Tabulated  coordinate
REAL,    INTENT(IN), DIMENSION(:,:) :: ZT        ! Tabulated data
REAL,    INTENT(IN), DIMENSION(KN)  :: PDX, PDY  ! 
!
REAL,                DIMENSION(KN)  :: PY         ! Interpolated value
!
!*	0.2	Declaration of local variables
!
INTEGER :: JJ        ! Loop index
!
!
!*	1.	INTERPOLATION
!		-------------
!  
DO JJ = 1, KN
  PY(JJ) = (1.0 - PDX(JJ)) * (1.0 - PDY(JJ)) * ZT(KI(JJ),  KJ(JJ))   + &
            PDX(JJ)        * (1.0 - PDY(JJ)) * ZT(KI(JJ)+1,KJ(JJ))   + &
            PDX(JJ)        * PDY(JJ)         * ZT(KI(JJ)+1,KJ(JJ)+1) + &
           (1.0 - PDX(JJ)) * PDY(JJ)         * ZT(KI(JJ)  ,KJ(JJ)+1)
ENDDO
!
END FUNCTION BI_LIN_INTP_V
!
!------------------------------------------------------------------
!
END SUBROUTINE ELEC_TENDENCIES
END MODULE MODE_ELEC_TENDENCIES
