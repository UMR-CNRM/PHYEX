!      #########################
       MODULE MODI_LIMA_PHILLIPS
!      #########################
!
INTERFACE
      SUBROUTINE LIMA_PHILLIPS (OHHONI, PTSTEP, KMI, HFMFILE, HLUOUT, OCLOSE_OUT, &
                                PZZ, PRHODJ, PRHODREF, PEXNREF, PPABST,           &
                                PTHT, PRVT, PRCT, PRRT, PRIT, PRST, PRGT,         &
                                PTHS, PRVS, PRCS, PRIS,                           &
                                PCIT, PCCS, PCIS,                                 &
                                PNAS, PIFS, PINS, PNIS, &
                            YDDDH, YDLDDH, YDMDDH   )
!
USE DDH_MIX, ONLY  : TYP_DDH
USE YOMLDDH, ONLY  : TLDDH
USE YOMMDDH, ONLY  : TMDDH
!
LOGICAL,                  INTENT(IN)    :: OHHONI  ! enable haze freezing
REAL,                     INTENT(IN)    :: PTSTEP  ! Time step          
INTEGER,                  INTENT(IN)    :: KMI     ! Model index 
CHARACTER(LEN=*),         INTENT(IN)    :: HFMFILE ! Name of the output FM-file
CHARACTER(LEN=*),         INTENT(IN)    :: HLUOUT  ! Output-listing name for
                                                   ! model n
LOGICAL,                  INTENT(IN)    :: OCLOSE_OUT ! Conditional closure of 
                                                   ! the tput FM fileoutp
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PZZ     ! Height (z)
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PEXNREF ! Reference Exner function
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABST  ! abs. pressure at time t
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT    ! Theta at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRVT    ! Water vapor m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRCT    ! Cloud water m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRRT    ! Rain water m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRIT    ! Cloud ice m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRST    ! Snow/aggregate m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRGT    ! Graupel m.r. at t 
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PTHS    ! Theta source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRVS    ! Water vapor m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRCS    ! Cloud water m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRIS    ! Pristine ice m.r. source
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCIT    ! Ice crystal C. at t
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCCS    ! Cloud water C. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCIS    ! Ice crystal C. source
!
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PNAS    ! Cloud  C. nuclei C. source
                                                   !used as Free ice nuclei for
                                                   !IMMERSION freezing
!
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PIFS    ! Free ice nuclei C. source 
                                                   !for DEPOSITION and CONTACT
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PINS    ! Activated ice nuclei C. source
                                                   !for DEPOSITION and CONTACT
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PNIS    ! Activated ice nuclei C. source
                                                   !for IMMERSION
!
TYPE(TYP_DDH), INTENT(INOUT) :: YDDDH
TYPE(TLDDH), INTENT(IN) :: YDLDDH
TYPE(TMDDH), INTENT(IN) :: YDMDDH
!
END SUBROUTINE LIMA_PHILLIPS
END INTERFACE
END MODULE MODI_LIMA_PHILLIPS
!
!     ######################################################################
      SUBROUTINE LIMA_PHILLIPS (OHHONI, PTSTEP, KMI, HFMFILE, HLUOUT, OCLOSE_OUT, &
                                PZZ, PRHODJ, PRHODREF, PEXNREF, PPABST,           &
                                PTHT, PRVT, PRCT, PRRT, PRIT, PRST, PRGT,         &
                                PTHS, PRVS, PRCS, PRIS,                           &
                                PCIT, PCCS, PCIS,                                 &
                                PNAS, PIFS, PINS, PNIS, &
                            YDDDH, YDLDDH, YDMDDH   )
!     ######################################################################
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the heterogeneous nucleation
!!    following Phillips (2008).
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
!!      Original             ??/??/13 
!!      C. Barthe  * LACy *  jan. 2014   add budgets
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAMETERS,      ONLY : JPHEXT, JPVEXT
USE MODD_CST,             ONLY : XP00, XRD, XMV, XMD, XCPD, XCPV, XCL, XCI,        &
                                 XTT, XLSTT, XLVTT, XALPI, XBETAI, XGAMI,          &
                                 XALPW, XBETAW, XGAMW, XPI
USE MODD_PARAM_LIMA,      ONLY : NMOD_IFN, NSPECIE, XFRAC,                         &
                                 NMOD_CCN, NMOD_IMM, NIND_SPECIE, NINDICE_CCN_IMM,  & 
                                 XDSI0, XRTMIN, XCTMIN, NPHILLIPS
USE MODD_PARAM_LIMA_COLD, ONLY : XMNU0
!
USE MODI_LIMA_FUNCTIONS,  ONLY : COUNTJV
USE MODI_LIMA_PHILLIPS_REF_SPECTRUM
USE MODI_LIMA_PHILLIPS_INTEG
!
USE MODD_BUDGET
USE MODE_BUDGET, ONLY: BUDGET_DDH
USE MODD_NSV, ONLY : NSV_LIMA_NC, NSV_LIMA_NI, NSV_LIMA_IFN_FREE
!
USE DDH_MIX, ONLY  : TYP_DDH
USE YOMLDDH, ONLY  : TLDDH
USE YOMMDDH, ONLY  : TMDDH
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
LOGICAL,                  INTENT(IN)    :: OHHONI  ! enable haze freezing
REAL,                     INTENT(IN)    :: PTSTEP  ! Time step          
INTEGER,                  INTENT(IN)    :: KMI     ! Model index 
CHARACTER(LEN=*),         INTENT(IN)    :: HFMFILE ! Name of the output FM-file
CHARACTER(LEN=*),         INTENT(IN)    :: HLUOUT  ! Output-listing name for
                                                   ! model n
LOGICAL,                  INTENT(IN)    :: OCLOSE_OUT ! Conditional closure of 
                                                   ! the tput FM fileoutp
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PZZ     ! Height (z)
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PEXNREF ! Reference Exner function
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABST  ! abs. pressure at time t
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT    ! Theta at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRVT    ! Water vapor m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRCT    ! Cloud water m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRRT    ! Rain water m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRIT    ! Cloud ice m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRST    ! Snow/aggregate m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRGT    ! Graupel m.r. at t 
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PTHS    ! Theta source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRVS    ! Water vapor m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRCS    ! Cloud water m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRIS    ! Pristine ice m.r. source
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCIT    ! Ice crystal C. at t
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCCS    ! Cloud water C. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCIS    ! Ice crystal C. source
!
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PNAS    ! Cloud  C. nuclei C. source
                                                   !used as Free ice nuclei for
                                                   !IMMERSION freezing
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PIFS    ! Free ice nuclei C. source 
                                                   !for DEPOSITION and CONTACT
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PINS    ! Activated ice nuclei C. source
                                                   !for DEPOSITION and CONTACT
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PNIS    ! Activated ice nuclei C. source
                                                   !for IMMERSION
!
TYPE(TYP_DDH), INTENT(INOUT) :: YDDDH
TYPE(TLDDH), INTENT(IN) :: YDLDDH
TYPE(TMDDH), INTENT(IN) :: YDMDDH
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
REAL, DIMENSION(:),   ALLOCATABLE :: ZCIT    ! Pristine ice conc. at t
!
REAL, DIMENSION(:),   ALLOCATABLE :: ZRVS    ! Water vapor m.r. source
REAL, DIMENSION(:),   ALLOCATABLE :: ZRCS    ! Cloud water m.r. source
REAL, DIMENSION(:),   ALLOCATABLE :: ZRIS    ! Pristine ice m.r. source
REAL, DIMENSION(:),   ALLOCATABLE :: ZCCS    ! Cloud water conc. source
REAL, DIMENSION(:),   ALLOCATABLE :: ZCIS    ! Pristine ice conc. source
!
REAL, DIMENSION(:),   ALLOCATABLE :: ZTHS    ! Theta source
!
REAL, DIMENSION(:,:), ALLOCATABLE :: ZNAS    ! Cloud Cond. nuclei conc. source
REAL, DIMENSION(:,:), ALLOCATABLE :: ZIFS    ! Free Ice nuclei conc. source
REAL, DIMENSION(:,:), ALLOCATABLE :: ZINS    ! Nucleated Ice nuclei conc. source
                                             !by Deposition/Contact
REAL, DIMENSION(:,:), ALLOCATABLE :: ZNIS    ! Nucleated Ice nuclei conc. source 
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
REAL,    DIMENSION(:),   ALLOCATABLE :: ZRTMIN, ZCTMIN
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
! Physical limitations
!
ALLOCATE(ZRTMIN(SIZE(XRTMIN)))
ALLOCATE(ZCTMIN(SIZE(XCTMIN)))
ZRTMIN(:) = XRTMIN(:) / PTSTEP
ZCTMIN(:) = XCTMIN(:) / PTSTEP
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
ALLOCATE(ZCIT(INEGT))
!
ALLOCATE(ZRVS(INEGT)) 
ALLOCATE(ZRCS(INEGT))
ALLOCATE(ZRIS(INEGT))
!
ALLOCATE(ZTHS(INEGT))
!
ALLOCATE(ZCCS(INEGT))
ALLOCATE(ZCIS(INEGT))
!
ALLOCATE(ZNAS(INEGT,NMOD_CCN))
ALLOCATE(ZIFS(INEGT,NMOD_IFN))
ALLOCATE(ZINS(INEGT,NMOD_IFN))
ALLOCATE(ZNIS(INEGT,NMOD_IMM))
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
   ZCIT(JL) = PCIT(I1(JL),I2(JL),I3(JL))
!
   ZRVS(JL) = PRVS(I1(JL),I2(JL),I3(JL))
   ZRCS(JL) = PRCS(I1(JL),I2(JL),I3(JL))
   ZRIS(JL) = PRIS(I1(JL),I2(JL),I3(JL))
!
   ZTHS(JL) = PTHS(I1(JL),I2(JL),I3(JL))
!
   ZCCS(JL) = PCCS(I1(JL),I2(JL),I3(JL))
   ZCIS(JL) = PCIS(I1(JL),I2(JL),I3(JL))
!
   DO JMOD_CCN = 1, NMOD_CCN
      ZNAS(JL,JMOD_CCN) = PNAS(I1(JL),I2(JL),I3(JL),JMOD_CCN)
   ENDDO
   DO JMOD_IFN = 1, NMOD_IFN
      ZIFS(JL,JMOD_IFN) = PIFS(I1(JL),I2(JL),I3(JL),JMOD_IFN)
      ZINS(JL,JMOD_IFN) = PINS(I1(JL),I2(JL),I3(JL),JMOD_IFN)
   ENDDO
   DO JMOD_IMM = 1, NMOD_IMM
      ZNIS(JL,JMOD_IMM) = PNIS(I1(JL),I2(JL),I3(JL),JMOD_IMM)
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
      ZZX(:)=ZZX(:)+XFRAC(JSPECIE,JMOD_IFN)*(ZIFS(:,JMOD_IFN)+ZINS(:,JMOD_IFN))* &
                                            Z_FRAC_ACT(:,JSPECIE)
   END DO                   
! Now : ZZX(:) = number of activable AP.
! Activated AP at this time step = activable AP - already activated AP 
   ZZX(:) = MIN( ZIFS(:,JMOD_IFN), MAX( (ZZX(:)-ZINS(:,JMOD_IFN)),0.0 ))
! Correction BVIE division by PTSTEP ?
!   ZZW(:) = MIN( XMNU0*ZZX(:) / PTSTEP , ZRVS(:)     )
   ZZW(:) = MIN( XMNU0*ZZX(:), ZRVS(:)     )
!
! Update the concentrations and MMR
!   
   ZIFS(:,JMOD_IFN)     = ZIFS(:,JMOD_IFN) - ZZX(:)
   ZW(:,:,:)            = PIFS(:,:,:,JMOD_IFN)
   PIFS(:,:,:,JMOD_IFN) = UNPACK( ZIFS(:,JMOD_IFN), MASK=GNEGT(:,:,:),        &
                                                               FIELD=ZW(:,:,:) )
!
   ZINS(:,JMOD_IFN)     = ZINS(:,JMOD_IFN) + ZZX(:)
   ZW(:,:,:)            = PINS(:,:,:,JMOD_IFN)
   PINS(:,:,:,JMOD_IFN) = UNPACK( ZINS(:,JMOD_IFN), MASK=GNEGT(:,:,:),        &
                                                               FIELD=ZW(:,:,:) )
!
   ZRVS(:) = ZRVS(:) - ZZW(:)
   ZRIS(:) = ZRIS(:) + ZZW(:)
   ZTHS(:) = ZTHS(:) + ZZW(:)*ZLSFACT(:)    !-ZLVFACT(:)) ! f(L_s*(RVHNDI))
   ZCIS(:) = ZCIS(:) + ZZX(:)
END DO
!
!
! Budget storage
IF (NBUMOD==KMI .AND. LBU_ENABLE) THEN
  IF (LBUDGET_TH) CALL BUDGET_DDH (                                              &
                  UNPACK(ZTHS(:),MASK=GNEGT(:,:,:),FIELD=PTHS)*PRHODJ(:,:,:),&
                                                              4,'HIND_BU_RTH',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RV) CALL BUDGET_DDH (                                              &
                  UNPACK(ZRVS(:),MASK=GNEGT(:,:,:),FIELD=PRVS)*PRHODJ(:,:,:),&
                                                              6,'HIND_BU_RRV',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RI) CALL BUDGET_DDH (                                              &
                  UNPACK(ZRIS(:),MASK=GNEGT(:,:,:),FIELD=PRIS)*PRHODJ(:,:,:),&
                                                              9,'HIND_BU_RRI',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_SV) THEN
    CALL BUDGET_DDH ( UNPACK(ZCIS(:),MASK=GNEGT(:,:,:),FIELD=PCIS)*PRHODJ(:,:,:),&
                                                       12+NSV_LIMA_NI,'HIND_BU_RSV',YDDDH, YDLDDH, YDMDDH)
    IF (NMOD_IFN.GE.1) THEN
       DO JL=1, NMOD_IFN
          CALL BUDGET_DDH ( PIFS(:,:,:,JL)*PRHODJ(:,:,:),12+NSV_LIMA_IFN_FREE+JL-1,'HIND_BU_RSV',YDDDH, YDLDDH, YDMDDH) 
       END DO
    END IF
  END IF
END IF
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
      ZZW(:) = MIN( ZCCS(:) , ZNAS(:,JMOD_CCN) )
      ZZX(:)=  ( ZZW(:)+ZNIS(:,JMOD_IMM) ) * Z_FRAC_ACT(:,NIND_SPECIE)
! Now : ZZX(:) = number of activable AP.
! Activated AP at this time step = activable AP - already activated AP 
      ZZX(:) = MIN( ZZW(:), MAX( (ZZX(:)-ZNIS(:,JMOD_IMM)),0.0 ) )
! Correction BVIE division by PTSTEP ?
!      ZZY(:) = MIN( XMNU0*ZZX(:) / PTSTEP , ZRVS(:)     )
      ZZY(:) = MIN( XMNU0*ZZX(:) , ZRVS(:)     )
!
! Update the concentrations and MMR
!   
      ZNAS(:,JMOD_CCN)     = ZNAS(:,JMOD_CCN) - ZZX(:)
      ZW(:,:,:)            = PNAS(:,:,:,JMOD_CCN)
      PNAS(:,:,:,JMOD_CCN) = UNPACK(ZNAS(:,JMOD_CCN),MASK=GNEGT(:,:,:), &
                                                     FIELD=ZW(:,:,:))
      ZNIS(:,JMOD_IMM)     = ZNIS(:,JMOD_IMM) + ZZX(:)
      ZW(:,:,:)            = PNIS(:,:,:,JMOD_IMM)
      PNIS(:,:,:,JMOD_IMM) = UNPACK(ZNIS(:,JMOD_IMM),MASK=GNEGT(:,:,:), &
                                                     FIELD=ZW(:,:,:))
!
      ZRCS(:) = ZRCS(:) - ZZY(:)
      ZRIS(:) = ZRIS(:) + ZZY(:)
      ZTHS(:) = ZTHS(:) + ZZY(:)*ZLSFACT(:) !-ZLVFACT(:)) ! f(L_s*(RVHNCI))
      ZCCS(:) = ZCCS(:) - ZZX(:)
      ZCIS(:) = ZCIS(:) + ZZX(:)
   END IF
END DO
!
! Budget storage
IF (NBUMOD==KMI .AND. LBU_ENABLE) THEN
  IF (LBUDGET_TH) CALL BUDGET_DDH (                                                 &
                  UNPACK(ZTHS(:),MASK=GNEGT(:,:,:),FIELD=PTHS)*PRHODJ(:,:,:),&
                                                              4,'HINC_BU_RTH',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RC) CALL BUDGET_DDH (                                                 &
                  UNPACK(ZRCS(:),MASK=GNEGT(:,:,:),FIELD=PRCS)*PRHODJ(:,:,:),&
                                                              7,'HINC_BU_RRC',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RI) CALL BUDGET_DDH (                                                 &
                  UNPACK(ZRIS(:),MASK=GNEGT(:,:,:),FIELD=PRIS)*PRHODJ(:,:,:),&
                                                              9,'HINC_BU_RRI',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_SV) THEN
    CALL BUDGET_DDH ( UNPACK(ZCCS(:),MASK=GNEGT(:,:,:),FIELD=PCCS)*PRHODJ(:,:,:),&
                                                       12+NSV_LIMA_NC,'HINC_BU_RSV',YDDDH, YDLDDH, YDMDDH)
    CALL BUDGET_DDH ( UNPACK(ZCIS(:),MASK=GNEGT(:,:,:),FIELD=PCIS)*PRHODJ(:,:,:),&
                                                       12+NSV_LIMA_NI,'HINC_BU_RSV',YDDDH, YDLDDH, YDMDDH)
  END IF
END IF
!
!
!-------------------------------------------------------------------------------
!
!
!*       7.     UNPACK VARIABLES AND CLEAN
!	        --------------------------
!
!
! End of the heterogeneous nucleation following Phillips 08
! Unpack variables, deallocate...
!
!
ZW(:,:,:)   = PRVS(:,:,:)
PRVS(:,:,:) = UNPACK( ZRVS(:),MASK=GNEGT(:,:,:),FIELD=ZW(:,:,:) )
ZW(:,:,:)   = PRCS(:,:,:)
PRCS(:,:,:) = UNPACK( ZRCS(:),MASK=GNEGT(:,:,:),FIELD=ZW(:,:,:) )
ZW(:,:,:)   = PRIS(:,:,:)
PRIS(:,:,:) = UNPACK( ZRIS(:),MASK=GNEGT(:,:,:),FIELD=ZW(:,:,:) )
ZW(:,:,:)   = PTHS(:,:,:)
PTHS(:,:,:) = UNPACK( ZTHS(:),MASK=GNEGT(:,:,:),FIELD=ZW(:,:,:) )
ZW(:,:,:)   = PCCS(:,:,:)
PCCS(:,:,:) = UNPACK( ZCCS(:),MASK=GNEGT(:,:,:),FIELD=ZW(:,:,:) )
ZW(:,:,:)   = PCIS(:,:,:)
PCIS(:,:,:) = UNPACK( ZCIS(:),MASK=GNEGT(:,:,:),FIELD=ZW(:,:,:) )
!
DEALLOCATE(ZRTMIN)
DEALLOCATE(ZCTMIN)
DEALLOCATE(ZRVT) 
DEALLOCATE(ZRCT) 
DEALLOCATE(ZRRT) 
DEALLOCATE(ZRIT) 
DEALLOCATE(ZRST) 
DEALLOCATE(ZRGT) 
DEALLOCATE(ZCIT)
DEALLOCATE(ZRVS) 
DEALLOCATE(ZRCS)
DEALLOCATE(ZRIS)
DEALLOCATE(ZTHS)
DEALLOCATE(ZCCS)
DEALLOCATE(ZCIS)
DEALLOCATE(ZNAS)
DEALLOCATE(ZIFS)
DEALLOCATE(ZINS)
DEALLOCATE(ZNIS)
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
!
ELSE
!
! Advance the budget calls
!
  IF (NBUMOD==KMI .AND. LBU_ENABLE) THEN
    IF (LBUDGET_TH) CALL BUDGET_DDH (PTHS(:,:,:)*PRHODJ(:,:,:),4,'HIND_BU_RTH',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RV) CALL BUDGET_DDH (PRVS(:,:,:)*PRHODJ(:,:,:),6,'HIND_BU_RRV',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RI) CALL BUDGET_DDH (PRIS(:,:,:)*PRHODJ(:,:,:),9,'HIND_BU_RRI',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_SV) THEN
       !print*, 'LBUDGET_SV dans lima_phillips = ', LBUDGET_SV
       CALL BUDGET_DDH (PCIS(:,:,:)*PRHODJ(:,:,:),12+NSV_LIMA_NI,'HIND_BU_RSV',YDDDH, YDLDDH, YDMDDH)
       IF (NMOD_IFN.GE.1) THEN
         DO JL=1, NMOD_IFN
            CALL BUDGET_DDH ( PIFS(:,:,:,JL)*PRHODJ(:,:,:),12+NSV_LIMA_IFN_FREE+JL-1,'HIND_BU_RSV',YDDDH, YDLDDH, YDMDDH) 
         END DO
       END IF
    END IF

    IF (LBUDGET_TH) CALL BUDGET_DDH (PTHS(:,:,:)*PRHODJ(:,:,:),4,'HINC_BU_RTH',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RC) CALL BUDGET_DDH (PRCS(:,:,:)*PRHODJ(:,:,:),7,'HINC_BU_RRC',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RI) CALL BUDGET_DDH (PRIS(:,:,:)*PRHODJ(:,:,:),9,'HINC_BU_RRI',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_SV) THEN
       !print*, 'LBUDGET_SV dans lima_phillips = ', LBUDGET_SV
       CALL BUDGET_DDH (PCCS(:,:,:)*PRHODJ(:,:,:),12+NSV_LIMA_NC,'HINC_BU_RSV',YDDDH, YDLDDH, YDMDDH)
       CALL BUDGET_DDH (PCIS(:,:,:)*PRHODJ(:,:,:),12+NSV_LIMA_NI,'HINC_BU_RSV',YDDDH, YDLDDH, YDMDDH)
    END IF
  END IF
!
!
END IF ! INEGT > 0
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_PHILLIPS
