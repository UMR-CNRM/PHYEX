!      #######################
       MODULE MODI_LIMA_MEYERS
!      #######################
!
INTERFACE
      SUBROUTINE LIMA_MEYERS   (OHHONI, PTSTEP, KMI, HFMFILE, HLUOUT, OCLOSE_OUT, &
                                PZZ, PRHODJ, PRHODREF, PEXNREF, PPABST,           &
                                PTHT, PRVT, PRCT, PRRT, PRIT, PRST, PRGT, PCCT,   &
                                PTHS, PRVS, PRCS, PRIS,                           &
                                PCCS, PCIS, PINS, &
                            YDDDH, YDLDDH, YDMDDH )
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
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCCT    ! Cloud water C. at t
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PTHS    ! Theta source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRVS    ! Water vapor m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRCS    ! Cloud water m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRIS    ! Pristine ice m.r. source
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCCS    ! Cloud water C. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCIS    ! Ice crystal C. source
!
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PINS    ! Activated ice nuclei C. source
                                                   !for DEPOSITION and CONTACT
                                                   !for IMMERSION
TYPE(TYP_DDH), INTENT(INOUT) :: YDDDH
TYPE(TLDDH), INTENT(IN) :: YDLDDH
TYPE(TMDDH), INTENT(IN) :: YDMDDH
!
END SUBROUTINE LIMA_MEYERS
END INTERFACE
END MODULE MODI_LIMA_MEYERS
!
!     ######################################################################
      SUBROUTINE LIMA_MEYERS   (OHHONI, PTSTEP, KMI, HFMFILE, HLUOUT, OCLOSE_OUT, &
                                PZZ, PRHODJ, PRHODREF, PEXNREF, PPABST,           &
                                PTHT, PRVT, PRCT, PRRT, PRIT, PRST, PRGT, PCCT,   &
                                PTHS, PRVS, PRCS, PRIS,                           &
                                PCCS, PCIS, PINS, &
                            YDDDH, YDLDDH, YDMDDH )
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
USE MODD_PARAMETERS
USE MODD_CST
USE MODD_PARAM_LIMA
USE MODD_PARAM_LIMA_COLD
USE MODD_BUDGET
USE MODE_BUDGET, ONLY: BUDGET_DDH
USE MODD_NSV, ONLY : NSV_LIMA_NC, NSV_LIMA_NI
!
USE DDH_MIX, ONLY  : TYP_DDH
USE YOMLDDH, ONLY  : TLDDH
USE YOMMDDH, ONLY  : TMDDH
!
USE MODI_LIMA_FUNCTIONS,  ONLY : COUNTJV
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
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCCT    ! Cloud water C. at t
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PTHS    ! Theta source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRVS    ! Water vapor m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRCS    ! Cloud water m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRIS    ! Pristine ice m.r. source
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCCS    ! Cloud water C. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCIS    ! Ice crystal C. source
!
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PINS    ! Activated ice nuclei C. source
                                                   !for DEPOSITION and CONTACT
TYPE(TYP_DDH), INTENT(INOUT) :: YDDDH
TYPE(TLDDH), INTENT(IN) :: YDLDDH
TYPE(TMDDH), INTENT(IN) :: YDMDDH
!
!
!*       0.2   Declarations of local variables :
!
!
INTEGER :: IIB, IIE, IJB, IJE, IKB, IKE               ! Physical domain
INTEGER :: JL     ! Loop index
INTEGER :: INEGT  ! Case number of nucleation
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
REAL, DIMENSION(:),   ALLOCATABLE :: ZRVS    ! Water vapor m.r. source
REAL, DIMENSION(:),   ALLOCATABLE :: ZRCS    ! Cloud water m.r. source
REAL, DIMENSION(:),   ALLOCATABLE :: ZRIS    ! Pristine ice m.r. source
REAL, DIMENSION(:),   ALLOCATABLE :: ZCCS    ! Cloud water conc. source
REAL, DIMENSION(:),   ALLOCATABLE :: ZCIS    ! Pristine ice conc. source
!
REAL, DIMENSION(:),   ALLOCATABLE :: ZTHS    ! Theta source
!
REAL, DIMENSION(:,:), ALLOCATABLE :: ZINS    ! Nucleated Ice nuclei conc. source
                                             ! by Deposition/Contact
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
                              ZSSI
!
REAL,    DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3))   &
                                  :: ZW, ZT ! work arrays
!
REAL,    DIMENSION(:),   ALLOCATABLE :: ZTCELSIUS
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
!  optimization by looking for locations where
!  the temperature is negative only !!!
!
GNEGT(:,:,:) = .FALSE.
GNEGT(IIB:IIE,IJB:IJE,IKB:IKE) = ZT(IIB:IIE,IJB:IJE,IKB:IKE)<XTT .AND. &
                                 ZW(IIB:IIE,IJB:IJE,IKB:IKE)>0.8
INEGT = COUNTJV( GNEGT(:,:,:),I1(:),I2(:),I3(:))
IF( INEGT >= 1 ) THEN
  ALLOCATE(ZRVT(INEGT)) 
  ALLOCATE(ZRCT(INEGT)) 
  ALLOCATE(ZRRT(INEGT)) 
  ALLOCATE(ZRIT(INEGT)) 
  ALLOCATE(ZRST(INEGT)) 
  ALLOCATE(ZRGT(INEGT)) 
!
  ALLOCATE(ZCCT(INEGT))
!
  ALLOCATE(ZRVS(INEGT)) 
  ALLOCATE(ZRCS(INEGT))
  ALLOCATE(ZRIS(INEGT))
!
  ALLOCATE(ZTHS(INEGT))
!
  ALLOCATE(ZCCS(INEGT))
  ALLOCATE(ZINS(INEGT,1))
  ALLOCATE(ZCIS(INEGT))
!
  ALLOCATE(ZRHODREF(INEGT)) 
  ALLOCATE(ZZT(INEGT)) 
  ALLOCATE(ZPRES(INEGT)) 
  ALLOCATE(ZEXNREF(INEGT))
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
    ZRVS(JL) = PRVS(I1(JL),I2(JL),I3(JL))
    ZRCS(JL) = PRCS(I1(JL),I2(JL),I3(JL))
    ZRIS(JL) = PRIS(I1(JL),I2(JL),I3(JL))
!
    ZTHS(JL) = PTHS(I1(JL),I2(JL),I3(JL))
!
    ZCCS(JL) = PCCS(I1(JL),I2(JL),I3(JL))
    ZCIS(JL) = PCIS(I1(JL),I2(JL),I3(JL))
!
    ZRHODREF(JL) = PRHODREF(I1(JL),I2(JL),I3(JL))
    ZZT(JL)      = ZT(I1(JL),I2(JL),I3(JL))
    ZPRES(JL)    = PPABST(I1(JL),I2(JL),I3(JL))
    ZEXNREF(JL)  = PEXNREF(I1(JL),I2(JL),I3(JL))
  ENDDO
  ALLOCATE(ZZW(INEGT))
  ALLOCATE(ZZX(INEGT))
  ALLOCATE(ZZY(INEGT))
  ALLOCATE(ZLSFACT(INEGT))
  ALLOCATE(ZLVFACT(INEGT))
  ALLOCATE(ZSSI(INEGT))
  ALLOCATE(ZTCELSIUS(INEGT))
!
  ZZW(:)  = ZEXNREF(:)*( XCPD+XCPV*ZRVT(:)+XCL*(ZRCT(:)+ZRRT(:)) &
                                  +XCI*(ZRIT(:)+ZRST(:)+ZRGT(:)) )
  ZTCELSIUS(:) = MAX( ZZT(:)-XTT,-50.0 )
  ZLSFACT(:) = (XLSTT+(XCPV-XCI)*ZTCELSIUS(:))/ZZW(:) ! L_s/(Pi_ref*C_ph)
  ZLVFACT(:) = (XLVTT+(XCPV-XCL)*ZTCELSIUS(:))/ZZW(:) ! L_v/(Pi_ref*C_ph)
!
  ZZW(:) = EXP( XALPI - XBETAI/ZZT(:) - XGAMI*ALOG(ZZT(:)) ) ! es_i
  ZSSI(:) = ZRVT(:)*(ZPRES(:)-ZZW(:))/((XMV/XMD)*ZZW(:)) - 1.0
                                                    ! Supersaturation over ice
!
!*            compute the heterogeneous nucleation by deposition: RVHNDI
!
  DO JL=1,INEGT
    ZINS(JL,1) = PINS(I1(JL),I2(JL),I3(JL),1)
  END DO
  ZZW(:) = 0.0
  ZZX(:) = 0.0
  ZZY(:) = 0.0
!
  WHERE( ZZT(:)<XTT-5.0 .AND. ZSSI(:)>0.0 )
    ZZY(:) = XNUC_DEP*EXP( XEXSI_DEP*100.*MIN(1.,ZSSI(:))+XEX_DEP)/(PTSTEP*ZRHODREF(:))
    ZZX(:) = MAX( ZZY(:)-ZINS(:,1) , 0.0 )
    ZZW(:) = MIN( XMNU0*ZZX(:) , ZRVS(:) )
  END WHERE
!
  ZINS(:,1)     = ZINS(:,1) + ZZX(:)
  ZW(:,:,:)     = PINS(:,:,:,1)
  PINS(:,:,:,1) = UNPACK( ZINS(:,1), MASK=GNEGT(:,:,:), FIELD=ZW(:,:,:)  )
!
  ZRVS(:) = ZRVS(:) - ZZW(:)
  ZRIS(:) = ZRIS(:) + ZZW(:)
  ZTHS(:) = ZTHS(:) + ZZW(:) * (ZLSFACT(:)-ZLVFACT(:)) ! f(L_s*(RVHNDI))
  ZCIS(:) = ZCIS(:) + ZZX(:)
!
!
! Budget storage
  IF (NBUMOD==KMI .AND. LBU_ENABLE) THEN
    IF (LBUDGET_TH) CALL BUDGET_DDH (                                                 &
                    UNPACK(ZTHS(:),MASK=GNEGT(:,:,:),FIELD=PTHS)*PRHODJ(:,:,:),&
                                                                4,'HIND_BU_RTH',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RV) CALL BUDGET_DDH (                                                 &
                    UNPACK(ZRVS(:),MASK=GNEGT(:,:,:),FIELD=PRVS)*PRHODJ(:,:,:),&
                                                                6,'HIND_BU_RRV',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RI) CALL BUDGET_DDH (                                                 &
                    UNPACK(ZRIS(:),MASK=GNEGT(:,:,:),FIELD=PRIS)*PRHODJ(:,:,:),&
                                                                9,'HIND_BU_RRI',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_SV) THEN
      CALL BUDGET_DDH ( UNPACK(ZCIS(:),MASK=GNEGT(:,:,:),FIELD=PCIS)*PRHODJ(:,:,:),&
                                                   12+NSV_LIMA_NI,'HIND_BU_RSV',YDDDH, YDLDDH, YDMDDH)
    END IF
  END IF
!
!*            compute the heterogeneous nucleation by contact: RVHNCI
!
  DO JL=1,INEGT
    ZINS(JL,1) = PINS(I1(JL),I2(JL),I3(JL),1)
  END DO
  ZZW(:) = 0.0
  ZZX(:) = 0.0
  ZZY(:) = 0.0
!
  WHERE( ZZT(:)<XTT-2.0 .AND. ZCCT(:)>XCTMIN(2) .AND. ZRCT(:)>XRTMIN(2) )
    ZZY(:) = MIN( XNUC_CON * EXP( XEXTT_CON*ZTCELSIUS(:)+XEX_CON )             &
                                               /(PTSTEP*ZRHODREF(:)) , ZCCS(:) )
    ZZX(:) = MAX( ZZY(:)-ZINS(:,1),0.0 )
    ZZW(:) = MIN( (ZRCT(:)/ZCCT(:))*ZZX(:),ZRCS(:) )
  END WHERE
!
  ZINS(:,1)     = ZINS(:,1) + ZZX(:)
  ZW(:,:,:)     = PINS(:,:,:,1)
  PINS(:,:,:,1) = UNPACK( ZINS(:,1), MASK=GNEGT(:,:,:), FIELD=ZW(:,:,:)  )
!
  ZRCS(:) = ZRCS(:) - ZZW(:)
  ZRIS(:) = ZRIS(:) + ZZW(:)
  ZTHS(:) = ZTHS(:) + ZZW(:)*(ZLSFACT(:)-ZLVFACT(:)) ! f(L_s*(RVHNCI))
  ZCCS(:) = ZCCS(:) - ZZX(:)
  ZCIS(:) = ZCIS(:) + ZZX(:)
!
!*            unpack variables
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
! Budget storage
  IF (NBUMOD==KMI .AND. LBU_ENABLE) THEN
    IF (LBUDGET_TH) CALL BUDGET_DDH (PTHS(:,:,:)*PRHODJ(:,:,:), 4,'HINC_BU_RTH',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RC) CALL BUDGET_DDH (PRCS(:,:,:)*PRHODJ(:,:,:), 7,'HINC_BU_RRC',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RI) CALL BUDGET_DDH (PRIS(:,:,:)*PRHODJ(:,:,:), 9,'HINC_BU_RRI',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_SV) THEN
      CALL BUDGET_DDH ( PCCS(:,:,:)*PRHODJ(:,:,:),12+NSV_LIMA_NC,'HINC_BU_RSV',YDDDH, YDLDDH, YDMDDH)
      CALL BUDGET_DDH ( PCIS(:,:,:)*PRHODJ(:,:,:),12+NSV_LIMA_NI,'HINC_BU_RSV',YDDDH, YDLDDH, YDMDDH)
    END IF
  END IF

!
  DEALLOCATE(ZRVT) 
  DEALLOCATE(ZRCT) 
  DEALLOCATE(ZRRT) 
  DEALLOCATE(ZRIT) 
  DEALLOCATE(ZRST) 
  DEALLOCATE(ZRGT) 
!
  DEALLOCATE(ZCCT)
!
  DEALLOCATE(ZRVS)
  DEALLOCATE(ZRCS)
  DEALLOCATE(ZRIS)
!
  DEALLOCATE(ZTHS)
!
  DEALLOCATE(ZCCS)
  DEALLOCATE(ZINS)
  DEALLOCATE(ZCIS)
!
  DEALLOCATE(ZRHODREF) 
  DEALLOCATE(ZZT) 
  DEALLOCATE(ZTCELSIUS)
  DEALLOCATE(ZPRES) 
  DEALLOCATE(ZEXNREF)
  DEALLOCATE(ZSSI)
  DEALLOCATE(ZZW)
  DEALLOCATE(ZZX)
  DEALLOCATE(ZZY)
  DEALLOCATE(ZLSFACT)
  DEALLOCATE(ZLVFACT)
!
ELSE
!
! Advance the budget calls
!
  IF (NBUMOD==KMI .AND. LBU_ENABLE) THEN
    IF (LBUDGET_TH)  CALL BUDGET_DDH (PTHS(:,:,:)*PRHODJ(:,:,:),4,'HIND_BU_RTH',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RV)  CALL BUDGET_DDH (PRVS(:,:,:)*PRHODJ(:,:,:),6,'HIND_BU_RRV',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RI)  CALL BUDGET_DDH (PRIS(:,:,:)*PRHODJ(:,:,:),9,'HIND_BU_RRI',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_SV)  CALL BUDGET_DDH (PCIS(:,:,:)*PRHODJ(:,:,:),12+NSV_LIMA_NI,'HIND_BU_RSV',YDDDH, YDLDDH, YDMDDH)

    IF (LBUDGET_TH)  CALL BUDGET_DDH (PTHS(:,:,:)*PRHODJ(:,:,:),4,'HINC_BU_RTH',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RC)  CALL BUDGET_DDH (PRCS(:,:,:)*PRHODJ(:,:,:),7,'HINC_BU_RRC',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RI)  CALL BUDGET_DDH (PRIS(:,:,:)*PRHODJ(:,:,:),9,'HINC_BU_RRI',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_SV) THEN
      CALL BUDGET_DDH (PCCS(:,:,:)*PRHODJ(:,:,:),12+NSV_LIMA_NC,'HINC_BU_RSV',YDDDH, YDLDDH, YDMDDH)
      CALL BUDGET_DDH (PCIS(:,:,:)*PRHODJ(:,:,:),12+NSV_LIMA_NI,'HINC_BU_RSV',YDDDH, YDLDDH, YDMDDH)
    END IF
  END IF
!
END IF




!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_MEYERS
