!MNH_LIC Copyright 2018-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!      ########################################
       MODULE MODI_LIMA_PHILLIPS_IFN_NUCLEATION
!      ########################################
!
INTERFACE
   SUBROUTINE LIMA_PHILLIPS_IFN_NUCLEATION (PTSTEP,                                   &
                                            PRHODREF, PEXNREF, PPABST,                &
                                            PTHT, PRVT, PRCT, PRRT, PRIT, PRST, PRGT, &
                                            PCCT, PCIT, PNAT, PIFT, PINT, PNIT,       &
                                            P_TH_HIND, P_RI_HIND, P_CI_HIND,          &
                                            P_TH_HINC, P_RC_HINC, P_CC_HINC,          &
                                            PICEFR                                    )
!
REAL,                     INTENT(IN)    :: PTSTEP 
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PEXNREF ! Reference Exner function
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABST  ! abs. pressure at time t
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PTHT    ! Theta at time t
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRVT    ! Water vapor m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRCT    ! Cloud water m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRRT    ! Rain water m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRIT    ! Cloud ice m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRST    ! Snow/aggregate m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRGT    ! Graupel m.r. at t
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCCT    ! Cloud water conc. at t 
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCIT    ! Cloud water conc. at t 
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PNAT    ! CCN conc. used for immersion nucl.
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PIFT    ! Free IFN conc.
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PINT    ! Nucleated IFN conc.
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PNIT    ! Nucleated (by immersion) CCN conc.
!
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: P_TH_HIND
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: P_RI_HIND
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: P_CI_HIND
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: P_TH_HINC
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: P_RC_HINC
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: P_CC_HINC
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PICEFR
!
END SUBROUTINE LIMA_PHILLIPS_IFN_NUCLEATION
END INTERFACE
END MODULE MODI_LIMA_PHILLIPS_IFN_NUCLEATION
!
!     #################################################################################
   SUBROUTINE LIMA_PHILLIPS_IFN_NUCLEATION (PTSTEP,                                   &
                                            PRHODREF, PEXNREF, PPABST,                &
                                            PTHT, PRVT, PRCT, PRRT, PRIT, PRST, PRGT, &
                                            PCCT, PCIT, PNAT, PIFT, PINT, PNIT,       &
                                            P_TH_HIND, P_RI_HIND, P_CI_HIND,          &
                                            P_TH_HINC, P_RC_HINC, P_CC_HINC,          &
                                            PICEFR                                    )
!     #################################################################################
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the heterogeneous nucleation
!!    following Phillips (2008) for the time-split version of LIMA
!!
!!
!!**  METHOD
!!    ------
!!      The parameterization of Phillips (2008) is based on observed nucleation
!!    in the CFDC for a range of T and Si values. Phillips therefore defines a 
!!    reference activity spectrum, that is, for given T and Si values, the 
!!    reference concentration of primary ice crystals.
!!      
!!      The activation of IFN is closely related to their total surface. Thus, 
!!    the activable fraction of each IFN specie is determined by an integration
!!    over the particle size distributions.
!!
!!    Subroutine organisation :
!!
!!      1- Preliminary computations
!!      2- Check where computations are necessary, and pack variables
!!      3- Compute the saturation over water and ice
!!      4- Compute the reference activity spectrum
!!             -> CALL LIMA_PHILLIPS_REF_SPECTRUM
!!         Integrate over the size distributions to compute the IFN activable fraction
!!             -> CALL LIMA_PHILLIPS_INTEG
!!      5- Heterogeneous nucleation of insoluble IFN
!!      6- Heterogeneous nucleation of coated IFN
!!      7- Unpack variables & deallocations
!! 
!!
!!    REFERENCE
!!    ---------
!!
!!      Phillips et al., 2008: An empirical parameterization of heterogeneous
!!        ice nucleation for multiple chemical species of aerosols, J. Atmos. Sci. 
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
!!      Original             15/03/2018
!  P. Wautelet 28/05/2019: move COUNTJV function to tools.f90
!  P. Wautelet 27/02/2020: bugfix: P_TH_HIND was not accumulated (will affect budgets) + add P_TH_HINC dummy argument
!                          + change intent of *_HIND and *_HINC dummy arguments (INOUT->OUT)
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST,             ONLY : XP00, XRD, XMV, XMD, XCPD, XCPV, XCL, XCI,        &
                                 XTT, XLSTT, XLVTT, XALPI, XBETAI, XGAMI,          &
                                 XALPW, XBETAW, XGAMW, XPI
USE MODD_NSV, ONLY : NSV_LIMA_NC, NSV_LIMA_NI, NSV_LIMA_IFN_FREE
USE MODD_PARAMETERS,      ONLY : JPHEXT, JPVEXT
USE MODD_PARAM_LIMA,      ONLY : NMOD_IFN, NSPECIE, XFRAC,                         &
                                 NMOD_CCN, NMOD_IMM, NIND_SPECIE, NINDICE_CCN_IMM,  &
                                 XDSI0, XRTMIN, XCTMIN, NPHILLIPS
USE MODD_PARAM_LIMA_COLD, ONLY : XMNU0

use mode_tools,           only: Countjv

USE MODI_LIMA_PHILLIPS_INTEG
USE MODI_LIMA_PHILLIPS_REF_SPECTRUM

IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
REAL,                     INTENT(IN)    :: PTSTEP 
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PEXNREF ! Reference Exner function
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABST  ! abs. pressure at time t
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PTHT    ! Theta at time t
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRVT    ! Water vapor m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRCT    ! Cloud water m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRRT    ! Rain water m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRIT    ! Cloud ice m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRST    ! Snow/aggregate m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRGT    ! Graupel m.r. at t
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCCT    ! Cloud water conc. at t 
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCIT    ! Cloud water conc. at t 
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PNAT    ! CCN conc. used for immersion nucl.
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PIFT    ! Free IFN conc.
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PINT    ! Nucleated IFN conc.
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PNIT    ! Nucleated (by immersion) CCN conc.
!
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: P_TH_HIND
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: P_RI_HIND
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: P_CI_HIND
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: P_TH_HINC
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: P_RC_HINC
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: P_CC_HINC
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PICEFR
!
!
!*       0.2   Declarations of local variables :
!
!
INTEGER :: IIB, IIE, IJB, IJE, IKB, IKE               ! Physical domain
INTEGER :: JL, JMOD_CCN, JMOD_IFN, JSPECIE, JMOD_IMM  ! Loop index
INTEGER :: INEGT  ! Case number of sedimentation, nucleation,
!
LOGICAL, DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3)) &
			  :: GNEGT  ! Test where to compute the nucleation
!
INTEGER, DIMENSION(SIZE(PRHODREF))  :: I1,I2,I3 ! Indexes for PACK replacement
!
REAL, DIMENSION(:),   ALLOCATABLE :: ZRVT    ! Water vapor m.r. at t
REAL, DIMENSION(:),   ALLOCATABLE :: ZRCT    ! Cloud water m.r. at t
REAL, DIMENSION(:),   ALLOCATABLE :: ZRRT    ! Rain water m.r. at t
REAL, DIMENSION(:),   ALLOCATABLE :: ZRIT    ! Pristine ice m.r. at t
REAL, DIMENSION(:),   ALLOCATABLE :: ZRST    ! Snow/aggregate m.r. at t
REAL, DIMENSION(:),   ALLOCATABLE :: ZRGT    ! Graupel/hail m.r. at t
!
REAL, DIMENSION(:),   ALLOCATABLE :: ZCCT    ! Cloud water conc. at t
!
REAL, DIMENSION(:,:), ALLOCATABLE :: ZNAT    ! Cloud Cond. nuclei conc. source
REAL, DIMENSION(:,:), ALLOCATABLE :: ZIFT    ! Free Ice nuclei conc. source
REAL, DIMENSION(:,:), ALLOCATABLE :: ZINT    ! Nucleated Ice nuclei conc. source
                                             !by Deposition/Contact
REAL, DIMENSION(:,:), ALLOCATABLE :: ZNIT    ! Nucleated Ice nuclei conc. source 
                                             !by Immersion
!
REAL, DIMENSION(:), ALLOCATABLE &
                           :: ZRHODREF, & ! RHO Dry REFerence
                              ZRHODJ,   & ! RHO times Jacobian
                              ZZT,      & ! Temperature
                              ZPRES,    & ! Pressure
                              ZEXNREF,  & ! EXNer Pressure REFerence
                              ZZW,      & ! Work array
                              ZZX,      & ! Work array
                              ZZY,      & ! Work array
                              ZLSFACT,  & ! L_s/(Pi_ref*C_ph)
                              ZLVFACT,  & ! L_v/(Pi_ref*C_ph)
                              ZLBDAC,   & ! Slope parameter of the cloud droplet distr.
                              ZSI,      &
                              ZSW,      &
                              ZSI_W
!
REAL,    DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3))   &
                                  :: ZW, ZT ! work arrays
!
REAL,    DIMENSION(:,:), ALLOCATABLE :: ZSI0, &    ! Si threshold in H_X for X={DM,BC,O}
                                        Z_FRAC_ACT ! Activable frac. of each AP species
REAL,    DIMENSION(:),   ALLOCATABLE :: ZTCELSIUS, ZZT_SI0_BC
!
!-------------------------------------------------------------------------------
!
!
!*       1.     PRELIMINARY COMPUTATIONS
!	        ------------------------
!
P_TH_HIND(:,:,:) = 0.
P_RI_HIND(:,:,:) = 0.
P_CI_HIND(:,:,:) = 0.
P_TH_HINC(:,:,:) = 0.
P_RC_HINC(:,:,:) = 0.
P_CC_HINC(:,:,:) = 0.
!
! Physical domain
!
IIB=1+JPHEXT
IIE=SIZE(PTHT,1) - JPHEXT
IJB=1+JPHEXT
IJE=SIZE(PTHT,2) - JPHEXT
IKB=1+JPVEXT
IKE=SIZE(PTHT,3) - JPVEXT
!
! Temperature
!
ZT(:,:,:)  = PTHT(:,:,:) * ( PPABST(:,:,:)/XP00 ) ** (XRD/XCPD)
!
! Saturation over ice
!
ZW(:,:,:) = EXP( XALPI - XBETAI/ZT(:,:,:) - XGAMI*ALOG(ZT(:,:,:) ) )
ZW(:,:,:) = PRVT(:,:,:)*( PPABST(:,:,:)-ZW(:,:,:) ) / ( (XMV/XMD) * ZW(:,:,:) )
!
!
!-------------------------------------------------------------------------------
!
!
!*       2.     COMPUTATIONS ONLY WHERE NECESSARY : PACK
!	        ----------------------------------------
!
!
GNEGT(:,:,:) = .FALSE.
GNEGT(IIB:IIE,IJB:IJE,IKB:IKE) = ZT(IIB:IIE,IJB:IJE,IKB:IKE)<XTT-2.0 .AND. &
                                 ZW(IIB:IIE,IJB:IJE,IKB:IKE)>0.95 
!
INEGT = COUNTJV( GNEGT(:,:,:),I1(:),I2(:),I3(:))
!
IF (INEGT > 0) THEN
!
   ALLOCATE(ZRVT(INEGT)) 
   ALLOCATE(ZRCT(INEGT)) 
   ALLOCATE(ZRRT(INEGT)) 
   ALLOCATE(ZRIT(INEGT)) 
   ALLOCATE(ZRST(INEGT)) 
   ALLOCATE(ZRGT(INEGT)) 
!
   ALLOCATE(ZCCT(INEGT)) 
!
   ALLOCATE(ZNAT(INEGT,NMOD_CCN))
   ALLOCATE(ZIFT(INEGT,NMOD_IFN))
   ALLOCATE(ZINT(INEGT,NMOD_IFN))
   ALLOCATE(ZNIT(INEGT,NMOD_IMM))
!
   ALLOCATE(ZRHODREF(INEGT)) 
   ALLOCATE(ZZT(INEGT)) 
   ALLOCATE(ZPRES(INEGT)) 
   ALLOCATE(ZEXNREF(INEGT))
!
   DO JL=1,INEGT
      ZRVT(JL) = PRVT(I1(JL),I2(JL),I3(JL))
      ZRCT(JL) = PRCT(I1(JL),I2(JL),I3(JL))
      ZRRT(JL) = PRRT(I1(JL),I2(JL),I3(JL))
      ZRIT(JL) = PRIT(I1(JL),I2(JL),I3(JL))
      ZRST(JL) = PRST(I1(JL),I2(JL),I3(JL))
      ZRGT(JL) = PRGT(I1(JL),I2(JL),I3(JL))
!
      ZCCT(JL) = PCCT(I1(JL),I2(JL),I3(JL))
!
      DO JMOD_CCN = 1, NMOD_CCN
         ZNAT(JL,JMOD_CCN) = PNAT(I1(JL),I2(JL),I3(JL),JMOD_CCN)
      ENDDO
      DO JMOD_IFN = 1, NMOD_IFN
         ZIFT(JL,JMOD_IFN) = PIFT(I1(JL),I2(JL),I3(JL),JMOD_IFN)
         ZINT(JL,JMOD_IFN) = PINT(I1(JL),I2(JL),I3(JL),JMOD_IFN)
      ENDDO
      DO JMOD_IMM = 1, NMOD_IMM
         ZNIT(JL,JMOD_IMM) = PNIT(I1(JL),I2(JL),I3(JL),JMOD_IMM)
      ENDDO
      ZRHODREF(JL) = PRHODREF(I1(JL),I2(JL),I3(JL))
      ZZT(JL)      = ZT(I1(JL),I2(JL),I3(JL))
      ZPRES(JL)    = PPABST(I1(JL),I2(JL),I3(JL))
      ZEXNREF(JL)  = PEXNREF(I1(JL),I2(JL),I3(JL))
   ENDDO
!
! PACK : done
! Prepare computations
!
   ALLOCATE( ZLSFACT    (INEGT) )
   ALLOCATE( ZLVFACT    (INEGT) )
   ALLOCATE( ZSI        (INEGT) )
   ALLOCATE( ZTCELSIUS  (INEGT) )
   ALLOCATE( ZZT_SI0_BC (INEGT) )
   ALLOCATE( ZLBDAC     (INEGT) )
   ALLOCATE( ZSI0       (INEGT,NSPECIE) )
   ALLOCATE( Z_FRAC_ACT (INEGT,NSPECIE) ) ; Z_FRAC_ACT(:,:) = 0.0
   ALLOCATE( ZSW        (INEGT) )
   ALLOCATE( ZSI_W      (INEGT) )
!
   ALLOCATE( ZZW (INEGT) ) ; ZZW(:) = 0.0
   ALLOCATE( ZZX (INEGT) ) ; ZZX(:) = 0.0
   ALLOCATE( ZZY (INEGT) ) ; ZZY(:) = 0.0
!
!
!-------------------------------------------------------------------------------
!
!
!*       3.     COMPUTE THE SATURATION OVER WATER AND ICE
!	        -----------------------------------------
!
!
   ZTCELSIUS(:) = ZZT(:)-XTT                                    ! T [°C]
   ZZW(:)  = ZEXNREF(:)*( XCPD+XCPV*ZRVT(:)+XCL*(ZRCT(:)+ZRRT(:)) &
        +XCI*(ZRIT(:)+ZRST(:)+ZRGT(:)) )
   ZLSFACT(:) = (XLSTT+(XCPV-XCI)*ZTCELSIUS(:))/ZZW(:)          ! L_s/(Pi_ref*C_ph)
   ZLVFACT(:) = (XLVTT+(XCPV-XCL)*ZTCELSIUS(:))/ZZW(:)          ! L_v/(Pi_ref*C_ph)
!
   ZZW(:)  = EXP( XALPI - XBETAI/ZZT(:) - XGAMI*ALOG(ZZT(:) ) ) ! es_i
   ZSI(:)  = ZRVT(:)*(ZPRES(:)-ZZW(:))/((XMV/XMD)*ZZW(:))       ! Saturation over ice
!
   ZZY(:)  = EXP( XALPW - XBETAW/ZZT(:) - XGAMW*ALOG(ZZT(:) ) ) ! es_w
   ZSW(:)  = ZRVT(:)*(ZPRES(:)-ZZY(:))/((XMV/XMD)*ZZY(:))       ! Saturation over water
!
   ZSI_W(:)= ZZY(:)/ZZW(:)     ! Saturation over ice at water saturation: es_w/es_i
!
! Saturation parameters for H_X, with X={Dust/Metallic (2 modes), Black Carbon, Organic}
!
   ZSI0(:,1) = 1.0 + 10.0**( -1.0261 + 3.1656E-3* ZTCELSIUS(:)     &
                                     + 5.3938E-4*(ZTCELSIUS(:)**2) &
                                     + 8.2584E-6*(ZTCELSIUS(:)**3) ) ! with T [°C]
   ZSI0(:,2) = ZSI0(:,1) ! DM2 = DM1
   ZSI0(:,3) = 0.0       ! BC
   ZZT_SI0_BC(:) = MAX( 198.0, MIN( 239.0,ZZT(:) ) )
   ZSI0(:,3) = (-3.118E-5*ZZT_SI0_BC(:)+1.085E-2)*ZZT_SI0_BC(:)+0.5652 - XDSI0(3)
   IF (NPHILLIPS == 8) THEN
      ZSI0(:,4) = ZSI0(:,3) ! O = BC
   ELSE IF (NPHILLIPS == 13) THEN
      ZSI0(:,4) = 1.15      ! BIO
   END IF
!
!
!-------------------------------------------------------------------------------
!
!
!*       4.     COMPUTE THE ACTIVABLE FRACTION OF EACH IFN SPECIE
!	        -------------------------------------------------
!
!
! Computation of the reference activity spectrum ( ZZY = N_{IN,1,*} )
!
   CALL LIMA_PHILLIPS_REF_SPECTRUM(ZZT, ZSI, ZSI_W, ZZY)
!
! For each aerosol species (DM1, DM2, BC, O), compute the fraction that may be activated
! Z_FRAC_ACT(INEGT,NSPECIE) = fraction of each species that may be activated
!
   CALL LIMA_PHILLIPS_INTEG(ZZT, ZSI, ZSI0, ZSW, ZZY, Z_FRAC_ACT)
!
!
!-------------------------------------------------------------------------------
!
!
!*       5.     COMPUTE THE HETEROGENEOUS NUCLEATION OF INSOLUBLE IFN
!	        -----------------------------------------------------
!
!
!
   DO JMOD_IFN = 1,NMOD_IFN    ! IFN modes
      ZZX(:)=0.
      DO JSPECIE = 1, NSPECIE  ! Each IFN mode is mixed with DM1, DM2, BC, O
         ZZX(:)=ZZX(:)+XFRAC(JSPECIE,JMOD_IFN)*(ZIFT(:,JMOD_IFN)+ZINT(:,JMOD_IFN))* &
                                               Z_FRAC_ACT(:,JSPECIE)
      END DO
! Now : ZZX(:) = number conc. of activable AP.
! Activated AP at this time step = activable AP - already activated AP 
      ZZX(:) = MIN( ZIFT(:,JMOD_IFN), MAX( (ZZX(:)-ZINT(:,JMOD_IFN)),0.0 ))
      ZZW(:) = MIN( XMNU0*ZZX(:), ZRVT(:) )
! Now : ZZX(:) = number conc. of AP activated at this time step (#/kg) from IFN mode JMOD_IFN
! Now : ZZW(:) = mmr of ice nucleated at this time step (kg/kg) from IFN mode JMOD_IFN
!
! Update the concentrations and MMR
!
      ZW(:,:,:) = UNPACK( ZZX(:), MASK=GNEGT(:,:,:), FIELD=0. )
      PIFT(:,:,:,JMOD_IFN) = PIFT(:,:,:,JMOD_IFN) - ZW(:,:,:)
      PINT(:,:,:,JMOD_IFN) = PINT(:,:,:,JMOD_IFN) + ZW(:,:,:)
!
      P_CI_HIND(:,:,:) = P_CI_HIND(:,:,:) + ZW(:,:,:)
      PCIT(:,:,:) = PCIT(:,:,:) + ZW(:,:,:)
!
      ZW(:,:,:) = UNPACK( ZZW(:), MASK=GNEGT(:,:,:), FIELD=0. )
      P_RI_HIND(:,:,:) = P_RI_HIND(:,:,:) + ZW(:,:,:)
      PRVT(:,:,:) = PRVT(:,:,:) - ZW(:,:,:)
      PRIT(:,:,:) = PRIT(:,:,:) + ZW(:,:,:)
!
      ZW(:,:,:) = UNPACK( ZZW(:)*ZLSFACT(:), MASK=GNEGT(:,:,:), FIELD=0. )
      P_TH_HIND(:,:,:) = P_TH_HIND(:,:,:) + ZW(:,:,:)
      PTHT(:,:,:) = PTHT(:,:,:) + ZW(:,:,:)
   END DO
!
!
!-------------------------------------------------------------------------------
!
!
!*       6.     COMPUTE THE HETEROGENEOUS NUCLEATION OF COATED IFN
!	        --------------------------------------------------
!
!
! Heterogeneous nucleation by immersion of the activated CCN
! Currently, we represent coated IFN as a pure aerosol type (NIND_SPECIE)
!
!
   DO JMOD_IMM = 1,NMOD_IMM  ! Coated IFN modes
      JMOD_CCN = NINDICE_CCN_IMM(JMOD_IMM) ! Corresponding CCN mode
      IF (JMOD_CCN .GT. 0) THEN
!
! OLD LIMA : Compute the appropriate mean diameter and sigma      
!      XMDIAM_IMM = MIN( XMDIAM_IFN(NIND_SPECIE) , XR_MEAN_CCN(JMOD_CCN)*2. )
!      XSIGMA_IMM = MIN( XSIGMA_IFN(JSPECIE) , EXP(XLOGSIG_CCN(JMOD_CCN)) )
!
         ZZW(:) = MIN( ZCCT(:) , ZNAT(:,JMOD_CCN) )
         ZZX(:)=  ( ZZW(:)+ZNIT(:,JMOD_IMM) ) * Z_FRAC_ACT(:,NIND_SPECIE)
! Now : ZZX(:) = number of activable AP.
! Activated AP at this time step = activable AP - already activated AP 
         ZZX(:) = MIN( ZZW(:), MAX( (ZZX(:)-ZNIT(:,JMOD_IMM)),0.0 ) )
         ZZY(:) = MIN( XMNU0*ZZX(:) , ZRVT(:), ZRCT(:) )
!
! Update the concentrations and MMR
!
         ZW(:,:,:) = UNPACK( ZZX(:), MASK=GNEGT(:,:,:), FIELD=0. )
         PNIT(:,:,:,JMOD_IMM) = PNIT(:,:,:,JMOD_IMM) + ZW(:,:,:)
         PNAT(:,:,:,JMOD_CCN) = PNAT(:,:,:,JMOD_CCN) - ZW(:,:,:)
!
         P_CC_HINC(:,:,:) = P_CC_HINC(:,:,:) - ZW(:,:,:) 
         PCCT(:,:,:) = PCCT(:,:,:) - ZW(:,:,:)
         PCIT(:,:,:) = PCIT(:,:,:) + ZW(:,:,:)
!
         ZW(:,:,:) = UNPACK( ZZY(:), MASK=GNEGT(:,:,:), FIELD=0. )
         P_RC_HINC(:,:,:) = P_RC_HINC(:,:,:) - ZW(:,:,:)
         PRCT(:,:,:) = PRCT(:,:,:) - ZW(:,:,:)
         PRIT(:,:,:) = PRIT(:,:,:) + ZW(:,:,:)
!
         ZW(:,:,:) = UNPACK( ZZY(:)*ZLSFACT(:), MASK=GNEGT(:,:,:), FIELD=0. )
         P_TH_HINC(:,:,:) = P_TH_HINC(:,:,:) + ZW(:,:,:)
         PTHT(:,:,:) = PTHT(:,:,:) + ZW(:,:,:)
      END IF
   END DO
!
!-------------------------------------------------------------------------------
!
!
!*       7.     CLEAN
!	        -----
!
!
   DEALLOCATE(ZRVT) 
   DEALLOCATE(ZRCT) 
   DEALLOCATE(ZRRT) 
   DEALLOCATE(ZRIT) 
   DEALLOCATE(ZRST) 
   DEALLOCATE(ZRGT) 
   DEALLOCATE(ZCCT) 
   DEALLOCATE(ZNAT)
   DEALLOCATE(ZIFT)
   DEALLOCATE(ZINT)
   DEALLOCATE(ZNIT)
   DEALLOCATE(ZRHODREF) 
   DEALLOCATE(ZZT) 
   DEALLOCATE(ZPRES) 
   DEALLOCATE(ZEXNREF)
   DEALLOCATE(ZLSFACT)
   DEALLOCATE(ZLVFACT)
   DEALLOCATE(ZSI)
   DEALLOCATE(ZTCELSIUS)
   DEALLOCATE(ZZT_SI0_BC)
   DEALLOCATE(ZLBDAC)
   DEALLOCATE(ZSI0)
   DEALLOCATE(Z_FRAC_ACT)
   DEALLOCATE(ZSW)
   DEALLOCATE(ZZW)
   DEALLOCATE(ZZX)
   DEALLOCATE(ZZY)
   DEALLOCATE(ZSI_W)
!
END IF ! INEGT > 0
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_PHILLIPS_IFN_NUCLEATION
