!MNH_LIC Copyright 2018-2024 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_LIMA_PHILLIPS_IFN_NUCLEATION
  IMPLICIT NONE
CONTAINS
!     #################################################################################
  SUBROUTINE LIMA_PHILLIPS_IFN_NUCLEATION (LIMAP, LIMAC, D, CST, PTSTEP,             &
                                           PRHODREF, PEXNREF, PPABST,                &
                                           PTHT, PRVT, PRCT, PRRT, PRIT, PRST, PRGT, &
                                           PCCT, PCIT, PCIT_SHAPE, PNAT, PIFT, PINT, PNIT, &
                                           P_TH_HIND, P_RI_HIND, P_CI_HIND, P_SHCI_HIND, &
                                           P_TH_HINC, P_RC_HINC, P_CC_HINC, P_SHCI_HINC, &
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
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_T
USE MODD_CST,            ONLY: CST_T

USE MODE_TOOLS,           only: COUNTJV

USE MODE_LIMA_PHILLIPS_INTEG, ONLY : LIMA_PHILLIPS_INTEG
USE MODE_LIMA_PHILLIPS_REF_SPECTRUM, ONLY : LIMA_PHILLIPS_REF_SPECTRUM
USE MODD_PARAM_LIMA_COLD, ONLY:PARAM_LIMA_COLD_T
USE MODD_PARAM_LIMA, ONLY:PARAM_LIMA_T
USE YOMHOOK, ONLY:LHOOK, DR_HOOK, JPHOOK

IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(PARAM_LIMA_COLD_T),INTENT(IN)::LIMAC
TYPE(PARAM_LIMA_T),INTENT(IN)::LIMAP
TYPE(DIMPHYEX_T),       INTENT(IN)    :: D
TYPE(CST_T),            INTENT(IN)    :: CST
REAL,                   INTENT(IN)    :: PTSTEP 
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PEXNREF ! Reference Exner function
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PPABST  ! abs. pressure at time t
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PTHT    ! Theta at time t
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PRVT    ! Water vapor m.r. at t 
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PRCT    ! Cloud water m.r. at t 
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PRRT    ! Rain water m.r. at t 
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PRIT    ! Cloud ice m.r. at t 
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PRST    ! Snow/aggregate m.r. at t 
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PRGT    ! Graupel m.r. at t
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PCCT    ! Cloud water conc. at t 
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PCIT    ! Cloud water conc. at t 
REAL, DIMENSION(D%NIJT,D%NKT,LIMAP%NNB_CRYSTAL_SHAPE),   INTENT(INOUT) :: PCIT_SHAPE ! Cloud water conc. at t 
REAL, DIMENSION(D%NIJT,D%NKT,LIMAP%NMOD_CCN), INTENT(INOUT) :: PNAT    ! CCN conc. used for immersion nucl.
REAL, DIMENSION(D%NIJT,D%NKT,LIMAP%NMOD_IFN), INTENT(INOUT) :: PIFT    ! Free IFN conc.
REAL, DIMENSION(D%NIJT,D%NKT,LIMAP%NMOD_IFN), INTENT(INOUT) :: PINT    ! Nucleated IFN conc.
REAL, DIMENSION(D%NIJT,D%NKT,LIMAP%NMOD_IMM), INTENT(INOUT) :: PNIT    ! Nucleated (by immersion) CCN conc.
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(OUT)   :: P_TH_HIND
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(OUT)   :: P_RI_HIND
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(OUT)   :: P_CI_HIND
REAL, DIMENSION(D%NIJT,D%NKT,LIMAP%NNB_CRYSTAL_SHAPE),   INTENT(OUT)   :: P_SHCI_HIND
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(OUT)   :: P_TH_HINC
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(OUT)   :: P_RC_HINC
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(OUT)   :: P_CC_HINC
REAL, DIMENSION(D%NIJT,D%NKT,LIMAP%NNB_CRYSTAL_SHAPE),   INTENT(OUT)   :: P_SHCI_HINC
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PICEFR
!
!
!*       0.2   Declarations of local variables :
!
!
INTEGER :: IIJB, IIJE, IKB, IKE               ! Physical domain
INTEGER :: IL, IMOD_CCN, IMOD_IFN, ISPECIE, IMOD_IMM, ISH  ! Loop index
INTEGER :: INEGT  ! Case number of sedimentation, nucleation,
!
LOGICAL, DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2)) &
                       :: GNEGT  ! Test where to compute the nucleation
!
INTEGER, DIMENSION(SIZE(PRHODREF))  :: I1,I3 ! Indexes for PACK replacement
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
REAL,    DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2))   &
                                  :: ZW, ZT ! work arrays
!
REAL,    DIMENSION(:,:), ALLOCATABLE :: ZSI0, &    ! Si threshold in H_X for X={DM,BC,O}
                                        Z_FRAC_ACT ! Activable frac. of each AP species
REAL,    DIMENSION(:),   ALLOCATABLE :: ZTCELSIUS, ZZT_SI0_BC
REAL, DIMENSION(:,:), ALLOCATABLE :: ZTC3D, & ! arrays of temperature and Si
                                     ZSI3D    ! --> used if lcrystal_shape=t
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
!
!*       1.     PRELIMINARY COMPUTATIONS
!               ------------------------
!
IF (LHOOK) CALL DR_HOOK('LIMA_PHILLIPS_IFN_NUCLEATION', 0, ZHOOK_HANDLE)
P_TH_HIND(:,:) = 0.
P_RI_HIND(:,:) = 0.
P_CI_HIND(:,:) = 0.
P_TH_HINC(:,:) = 0.
P_RC_HINC(:,:) = 0.
P_CC_HINC(:,:) = 0.
!
! Temperature
!
ZT(:,:)  = PTHT(:,:) * ( PPABST(:,:)/CST%XP00 ) ** (CST%XRD/CST%XCPD)
!
! Saturation over ice
!
ZW(:,:) = EXP( CST%XALPI - CST%XBETAI/ZT(:,:) - CST%XGAMI*ALOG(ZT(:,:) ) )
ZW(:,:) = PRVT(:,:)*( PPABST(:,:)-ZW(:,:) ) / ( (CST%XMV/CST%XMD) * ZW(:,:) )
!
!
!-------------------------------------------------------------------------------
!
!
!*       2.     COMPUTATIONS ONLY WHERE NECESSARY : PACK
!               ----------------------------------------
!
!
GNEGT(:,:) = .FALSE.
GNEGT(D%NIJB:D%NIJE,D%NKB:D%NKE) = ZT(D%NIJB:D%NIJE,D%NKB:D%NKE)<CST%XTT-2.0 .AND. &
                           ZW(D%NIJB:D%NIJE,D%NKB:D%NKE)>0.95 
!
INEGT = COUNTJV( GNEGT(:,:),I1(:),I3(:))
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
   ALLOCATE(ZNAT(INEGT,LIMAP%NMOD_CCN))
   ALLOCATE(ZIFT(INEGT,LIMAP%NMOD_IFN))
   ALLOCATE(ZINT(INEGT,LIMAP%NMOD_IFN))
   ALLOCATE(ZNIT(INEGT,LIMAP%NMOD_IMM))
!
   ALLOCATE(ZRHODREF(INEGT)) 
   ALLOCATE(ZZT(INEGT)) 
   ALLOCATE(ZPRES(INEGT)) 
   ALLOCATE(ZEXNREF(INEGT))
!
   DO IL=1,INEGT
      ZRVT(IL) = PRVT(I1(IL),I3(IL))
      ZRCT(IL) = PRCT(I1(IL),I3(IL))
      ZRRT(IL) = PRRT(I1(IL),I3(IL))
      ZRIT(IL) = PRIT(I1(IL),I3(IL))
      ZRST(IL) = PRST(I1(IL),I3(IL))
      ZRGT(IL) = PRGT(I1(IL),I3(IL))
!
      ZCCT(IL) = PCCT(I1(IL),I3(IL))
!
      DO IMOD_CCN = 1, LIMAP%NMOD_CCN
         ZNAT(IL,IMOD_CCN) = PNAT(I1(IL),I3(IL),IMOD_CCN)
      ENDDO
      DO IMOD_IFN = 1, LIMAP%NMOD_IFN
         ZIFT(IL,IMOD_IFN) = PIFT(I1(IL),I3(IL),IMOD_IFN)
         ZINT(IL,IMOD_IFN) = PINT(I1(IL),I3(IL),IMOD_IFN)
      ENDDO
      DO IMOD_IMM = 1, LIMAP%NMOD_IMM
         ZNIT(IL,IMOD_IMM) = PNIT(I1(IL),I3(IL),IMOD_IMM)
      ENDDO
      ZRHODREF(IL) = PRHODREF(I1(IL),I3(IL))
      ZZT(IL)      = ZT(I1(IL),I3(IL))
      ZPRES(IL)    = PPABST(I1(IL),I3(IL))
      ZEXNREF(IL)  = PEXNREF(I1(IL),I3(IL))
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
   ALLOCATE( ZSI0       (INEGT,LIMAP%NSPECIE) )
   ALLOCATE( Z_FRAC_ACT (INEGT,LIMAP%NSPECIE) )
   ALLOCATE( ZSW        (INEGT) )
   ALLOCATE( ZSI_W      (INEGT) )
!
   ALLOCATE( ZZW (INEGT) ) ; ZZW(:) = 0.0
   ALLOCATE( ZZX (INEGT) ) ; ZZX(:) = 0.0
   ALLOCATE( ZZY (INEGT) )
!
!
!-------------------------------------------------------------------------------
!
!
!*       3.     COMPUTE THE SATURATION OVER WATER AND ICE
!               -----------------------------------------
!
!
   ZTCELSIUS(:) = ZZT(:)-CST%XTT                                    ! T [°C]
   ZZW(:)  = ZEXNREF(:)*( CST%XCPD+CST%XCPV*ZRVT(:)+CST%XCL*(ZRCT(:)+ZRRT(:)) &
        +CST%XCI*(ZRIT(:)+ZRST(:)+ZRGT(:)) )
   ZLSFACT(:) = (CST%XLSTT+(CST%XCPV-CST%XCI)*ZTCELSIUS(:))/ZZW(:)          ! L_s/(Pi_ref*C_ph)
   ZLVFACT(:) = (CST%XLVTT+(CST%XCPV-CST%XCL)*ZTCELSIUS(:))/ZZW(:)          ! L_v/(Pi_ref*C_ph)
!
   ZZW(:)  = EXP( CST%XALPI - CST%XBETAI/ZZT(:) - CST%XGAMI*ALOG(ZZT(:) ) ) ! es_i
   ZSI(:)  = ZRVT(:)*(ZPRES(:)-ZZW(:))/((CST%XMV/CST%XMD)*ZZW(:))       ! Saturation over ice
!
   ZZY(:)  = EXP( CST%XALPW - CST%XBETAW/ZZT(:) - CST%XGAMW*ALOG(ZZT(:) ) ) ! es_w
   ZSW(:)  = ZRVT(:)*(ZPRES(:)-ZZY(:))/((CST%XMV/CST%XMD)*ZZY(:))       ! Saturation over water
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
   ZSI0(:,3) = (-3.118E-5*ZZT_SI0_BC(:)+1.085E-2)*ZZT_SI0_BC(:)+0.5652 - LIMAP%XDSI0(3)
   IF (LIMAP%NPHILLIPS == 8) THEN
      ZSI0(:,4) = ZSI0(:,3) ! O = BC
   ELSE IF (LIMAP%NPHILLIPS == 13) THEN
      ZSI0(:,4) = 1.15      ! BIO
   END IF
!
! if lcrystal_shape=t, 3D array of temperature and supersaturation are needed
   IF (LIMAP%LCRYSTAL_SHAPE) THEN
     ALLOCATE(ZTC3D(SIZE(PRHODREF,1),SIZE(PRHODREF,2)))
     ALLOCATE(ZSI3D(SIZE(PRHODREF,1),SIZE(PRHODREF,2)))
     ZTC3D(:,:) = ZT(:,:) - CST%XTT
     ZSI3D(:,:) = UNPACK( ZSI(:), MASK=GNEGT(:,:), FIELD=0. )
   END IF
!
!-------------------------------------------------------------------------------
!
!
!*       4.     COMPUTE THE ACTIVABLE FRACTION OF EACH IFN SPECIE
!               -------------------------------------------------
!
!
! Computation of the reference activity spectrum ( ZZY = N_{IN,1,*} )
!
   CALL LIMA_PHILLIPS_REF_SPECTRUM(LIMAP, CST, INEGT, ZZT, ZSI, ZSI_W, ZZY)
!
! For each aerosol species (DM1, DM2, BC, O), compute the fraction that may be activated
! Z_FRAC_ACT(INEGT,LIMAP%NSPECIE) = fraction of each species that may be activated
!
   CALL LIMA_PHILLIPS_INTEG(LIMAP, CST, INEGT, ZZT, ZSI, ZSI0, ZSW, ZZY, Z_FRAC_ACT)
!
!
!-------------------------------------------------------------------------------
!
!
!*       5.     COMPUTE THE HETEROGENEOUS NUCLEATION OF INSOLUBLE IFN
!               -----------------------------------------------------
!
!
!
   DO IMOD_IFN = 1,LIMAP%NMOD_IFN    ! IFN modes
      ZZX(:)=0.
      DO ISPECIE = 1, LIMAP%NSPECIE  ! Each IFN mode is mixed with DM1, DM2, BC, O
         ZZX(:)=ZZX(:)+LIMAP%XFRAC(ISPECIE,IMOD_IFN)*(ZIFT(:,IMOD_IFN)+ZINT(:,IMOD_IFN))* &
                                               Z_FRAC_ACT(:,ISPECIE)
      END DO
! Now : ZZX(:) = number conc. of activable AP.
! Activated AP at this time step = activable AP - already activated AP 
      ZZX(:) = MIN( ZIFT(:,IMOD_IFN), MAX( (ZZX(:)-ZINT(:,IMOD_IFN)),0.0 ))
      ZZW(:) = MIN( LIMAC%XMNU0*ZZX(:), ZRVT(:) )
! Now : ZZX(:) = number conc. of AP activated at this time step (#/kg) from IFN mode IMOD_IFN
! Now : ZZW(:) = mmr of ice nucleated at this time step (kg/kg) from IFN mode IMOD_IFN
!
! Update the concentrations and MMR
!
      ZW(:,:) = UNPACK( ZZX(:), MASK=GNEGT(:,:), FIELD=0. )
      PIFT(:,:,IMOD_IFN) = PIFT(:,:,IMOD_IFN) - ZW(:,:)
      PINT(:,:,IMOD_IFN) = PINT(:,:,IMOD_IFN) + ZW(:,:)
!
      P_CI_HIND(:,:) = P_CI_HIND(:,:) + ZW(:,:)
      IF (.NOT. LIMAP%LCRYSTAL_SHAPE) THEN
        PCIT(:,:) = PCIT(:,:) + ZW(:,:)
      ELSE
!NOTE : p_shci_hinX est utile uniquement pour les bilans --> peut-etre mettre une condition pour leur calcul ?
        ! different crystal habits are generated depending on the temperature
!++cb++ 18/04/24 la formation des cristaux produit uniquement des formes primaires
! on se base sur la figure 5 de Bailey et Hallett (2009)
! Plates: -1<T<-3, -9<T<-20, -20<T<-40, -40<T si SSI < 0.05
        WHERE (((ZTC3D(:,:) .GT. -3.0)  .AND.  (ZTC3D(:,:) .LT. 0.0))   .OR. &
               ((ZTC3D(:,:) .LE. -9.0)  .AND.  (ZTC3D(:,:) .GT. -40.0)) .OR. &
               ((ZTC3D(:,:) .LE. -40.0) .AND. ((ZSI3D(:,:)-1.) .LT. 0.05)) )
          P_SHCI_HIND(:,:,1) = P_SHCI_HIND(:,:,1) + ZW(:,:)
          PCIT_SHAPE(:,:,1)  = PCIT_SHAPE(:,:,1)  + ZW(:,:)
        END WHERE
!
! Columns: -3<T<-9, -40<T<-70 si SSI > 0.05
        WHERE (((ZTC3D(:,:) .LE. -3.0)   .AND.  (ZTC3D(:,:) .GT. -9.0)) .OR. &
               ((ZTC3D(:,:) .LE. -40.0)  .AND. ((ZSI3D(:,:)-1.) .GE. 0.05)))
          P_SHCI_HIND(:,:,2) = P_SHCI_HIND(:,:,2) + ZW(:,:)
          PCIT_SHAPE(:,:,2)  = PCIT_SHAPE(:,:,2)  + ZW(:,:)
        END WHERE
        !
        PCIT(:,:) = SUM(PCIT_SHAPE, DIM=3)
      END IF
!
      ZW(:,:) = UNPACK( ZZW(:), MASK=GNEGT(:,:), FIELD=0. )
      P_RI_HIND(:,:) = P_RI_HIND(:,:) + ZW(:,:)
      PRVT(:,:) = PRVT(:,:) - ZW(:,:)
      PRIT(:,:) = PRIT(:,:) + ZW(:,:)
!
      ZW(:,:) = UNPACK( ZZW(:)*ZLSFACT(:), MASK=GNEGT(:,:), FIELD=0. )
      P_TH_HIND(:,:) = P_TH_HIND(:,:) + ZW(:,:)
      PTHT(:,:) = PTHT(:,:) + ZW(:,:)
   END DO
!
!
!-------------------------------------------------------------------------------
!
!
!*       6.     COMPUTE THE HETEROGENEOUS NUCLEATION OF COATED IFN
!               --------------------------------------------------
!
!
! Heterogeneous nucleation by immersion of the activated CCN
! Currently, we represent coated IFN as a pure aerosol type (LIMAP%NIND_SPECIE)
!
!
   DO IMOD_IMM = 1,LIMAP%NMOD_IMM  ! Coated IFN modes
      IMOD_CCN = LIMAP%NINDICE_CCN_IMM(IMOD_IMM) ! Corresponding CCN mode
      IF (IMOD_CCN .GT. 0) THEN
!
! OLD LIMA : Compute the appropriate mean diameter and sigma      
!      XMDIAM_IMM = MIN( LIMAP%XMDIAM_IFN(LIMAP%NIND_SPECIE) , LIMAP%XR_MEAN_CCN(IMOD_CCN)*2. )
!      XSIGMA_IMM = MIN( LIMAP%XSIGMA_IFN(ISPECIE) , EXP(LIMAP%XLOGSIG_CCN(IMOD_CCN)) )
!
         ZZW(:) = MIN( ZCCT(:) , ZNAT(:,IMOD_CCN) )
         ZZX(:)=  ( ZZW(:)+ZNIT(:,IMOD_IMM) ) * Z_FRAC_ACT(:,LIMAP%NIND_SPECIE)
! Now : ZZX(:) = number of activable AP.
! Activated AP at this time step = activable AP - already activated AP 
         ZZX(:) = MIN( ZZW(:), MAX( (ZZX(:)-ZNIT(:,IMOD_IMM)),0.0 ) )
         ZZY(:) = MIN( LIMAC%XMNU0*ZZX(:) , ZRVT(:), ZRCT(:) )
!
! Update the concentrations and MMR
!
         ZW(:,:) = UNPACK( ZZX(:), MASK=GNEGT(:,:), FIELD=0. )
         PNIT(:,:,IMOD_IMM) = PNIT(:,:,IMOD_IMM) + ZW(:,:)
         PNAT(:,:,IMOD_CCN) = PNAT(:,:,IMOD_CCN) - ZW(:,:)
!
         P_CC_HINC(:,:) = P_CC_HINC(:,:) - ZW(:,:) 
         PCCT(:,:) = PCCT(:,:) - ZW(:,:)
         IF (.NOT. LIMAP%LCRYSTAL_SHAPE) THEN
           PCIT(:,:) = PCIT(:,:) + ZW(:,:)
         ELSE
           ! different crystal habits are generated depending on the temperature

! Plates: -1<T<-3, -9<T<-20, -20<T<-40, -40<T si SSI < 0.05
!           WHERE (((ZTC3D(:,:) .GT. -3.0) .AND.  (ZTC3D(:,:) .LT. 0.0))   .OR. &
!                 ((ZTC3D(:,:) .LE. -9.0)  .AND.  (ZTC3D(:,:) .GT. -40.0)) .OR. &
!                 ((ZTC3D(:,:) .LE. -40.0) .AND. ((ZSI3D(:,:)-1.) .LT. 0.05)) )
!             P_SHCI_HINC(:,:,1) = P_SHCI_HINC(:,:,1) + ZW(:,:)
!             PCIT_SHAPE(:,:,1)  = PCIT_SHAPE(:,:,1)  + ZW(:,:)
!           END WHERE
!
! Columns: -3<T<-9, -40<T<-70 si SSI > 0.05
!           WHERE (((ZTC3D(:,:) .LE. -3.0)   .AND.  (ZTC3D(:,:) .GT. -9.0)) .OR. &
!                  ((ZTC3D(:,:) .LE. -40.0)  .AND. ((ZSI3D(:,:)-1.) .GE. 0.05)))
!             P_SHCI_HINC(:,:,2) = P_SHCI_HINC(:,:,2) + ZW(:,:)
!             PCIT_SHAPE(:,:,2)  = PCIT_SHAPE(:,:,2)  + ZW(:,:)
!           END WHERE
! Droxtals
           P_SHCI_HINC(:,:,4) = P_SHCI_HINC(:,:,4) + ZW(:,:)
           PCIT_SHAPE(:,:,4)  = PCIT_SHAPE(:,:,4)  + ZW(:,:)
           !
           PCIT(:,:) = SUM(PCIT_SHAPE, DIM=3)
         END IF
!
         ZW(:,:) = UNPACK( ZZY(:), MASK=GNEGT(:,:), FIELD=0. )
         P_RC_HINC(:,:) = P_RC_HINC(:,:) - ZW(:,:)
         PRCT(:,:) = PRCT(:,:) - ZW(:,:)
         PRIT(:,:) = PRIT(:,:) + ZW(:,:)
!
         ZW(:,:) = UNPACK( ZZY(:)*ZLSFACT(:), MASK=GNEGT(:,:), FIELD=0. )
         P_TH_HINC(:,:) = P_TH_HINC(:,:) + ZW(:,:)
         PTHT(:,:) = PTHT(:,:) + ZW(:,:)
      END IF
   END DO
!
!-------------------------------------------------------------------------------
!
!
!*       7.     CLEAN
!               -----
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
   IF (ALLOCATED(ZTC3D)) DEALLOCATE(ZTC3D)
   IF (ALLOCATED(ZSI3D)) DEALLOCATE(ZSI3D)
!
END IF ! INEGT > 0
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('LIMA_PHILLIPS_IFN_NUCLEATION', 1, ZHOOK_HANDLE)
END SUBROUTINE LIMA_PHILLIPS_IFN_NUCLEATION
END MODULE MODE_LIMA_PHILLIPS_IFN_NUCLEATION
