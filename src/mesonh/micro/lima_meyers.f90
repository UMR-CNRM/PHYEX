!MNH_LIC Copyright 2013-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!      #######################
       MODULE MODI_LIMA_MEYERS
!      #######################
!
IMPLICIT NONE
INTERFACE
      SUBROUTINE LIMA_MEYERS   (OHHONI, PTSTEP, KMI,                            &
                                PZZ, PRHODJ, PRHODREF, PEXNREF, PPABST,         &
                                PTHT, PRVT, PRCT, PRRT, PRIT, PRST, PRGT, PCCT, &
                                PTHS, PRVS, PRCS, PRIS,                         &
                                PCCS, PCIS, PINS )
IMPLICIT NONE
!
LOGICAL,                  INTENT(IN)    :: OHHONI  ! enable haze freezing
REAL,                     INTENT(IN)    :: PTSTEP  ! Time step          
INTEGER,                  INTENT(IN)    :: KMI     ! Model index 
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
!
END SUBROUTINE LIMA_MEYERS
END INTERFACE
END MODULE MODI_LIMA_MEYERS
!
!     ###########################################################################
      SUBROUTINE LIMA_MEYERS   (OHHONI, PTSTEP, KMI,                            &
                                PZZ, PRHODJ, PRHODREF, PEXNREF, PPABST,         &
                                PTHT, PRVT, PRCT, PRRT, PRIT, PRST, PRGT, PCCT, &
                                PTHS, PRVS, PRCS, PRIS,                         &
                                PCCS, PCIS, PINS )
!     ###########################################################################
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
!  P. Wautelet 05/2016-04/2018: new data structures and calls for I/O
!  P. Wautelet 28/05/2019: move COUNTJV function to tools.f90
!  P. Wautelet    02/2020: use the new data structures and subroutines for budgets
!  P. Wautelet 02/02/2021: budgets: add missing source terms for SV budgets in LIMA
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
use modd_budget,          only: lbu_enable, nbumod,                                          &
                                lbudget_th, lbudget_rv, lbudget_rc, lbudget_ri, lbudget_sv,  &
                                NBUDGET_TH, NBUDGET_RV, NBUDGET_RC, NBUDGET_RI, NBUDGET_SV1, &
                                tbudgets
USE MODD_CST
USE MODD_NSV,             ONLY: NSV_LIMA_NC, NSV_LIMA_NI, NSV_LIMA_IFN_NUCL
USE MODD_PARAMETERS
USE MODD_PARAM_LIMA
USE MODD_PARAM_LIMA_COLD

use mode_budget,          only: Budget_store_init, Budget_store_end
use mode_tools,           only: Countjv

IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
LOGICAL,                  INTENT(IN)    :: OHHONI  ! enable haze freezing
REAL,                     INTENT(IN)    :: PTSTEP  ! Time step          
INTEGER,                  INTENT(IN)    :: KMI     ! Model index 
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
  if ( nbumod == kmi .and. lbu_enable ) then
    if ( lbudget_th ) call Budget_store_init( tbudgets(NBUDGET_TH),                    'HIND', pths(:, :, :) * prhodj(:, :, :) )
    if ( lbudget_rv ) call Budget_store_init( tbudgets(NBUDGET_RV),                    'HIND', prvs(:, :, :) * prhodj(:, :, :) )
    if ( lbudget_ri ) call Budget_store_init( tbudgets(NBUDGET_RI),                    'HIND', pris(:, :, :) * prhodj(:, :, :) )
    if ( lbudget_sv ) call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ni), 'HIND', pcis(:, :, :) * prhodj(:, :, :) )
    if ( lbudget_sv .and. nmod_ifn > 0) &
      call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ifn_nucl), 'HIND', pins(:, :, :, 1) * prhodj(:, :, :) )
  end if

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
!
  ZRVS(:) = ZRVS(:) - ZZW(:)
  ZRIS(:) = ZRIS(:) + ZZW(:)
  ZTHS(:) = ZTHS(:) + ZZW(:) * (ZLSFACT(:)-ZLVFACT(:)) ! f(L_s*(RVHNDI))
  ZCIS(:) = ZCIS(:) + ZZX(:)
!
!
! Budget storage
  if ( nbumod == kmi .and. lbu_enable ) then
    if ( lbudget_th ) call Budget_store_end( tbudgets(NBUDGET_TH), 'HIND', &
                                             Unpack ( zths(:), mask = gnegt(:, :, :), field = pths(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rv ) call Budget_store_end( tbudgets(NBUDGET_RV), 'HIND', &
                                             Unpack ( zrvs(:), mask = gnegt(:, :, :), field = prvs(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_ri ) call Budget_store_end( tbudgets(NBUDGET_RI), 'HIND', &
                                             Unpack ( zris(:), mask = gnegt(:, :, :), field = pris(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_sv ) call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ni), 'HIND', &
                                             Unpack ( zcis(:), mask = gnegt(:, :, :), field = pcis(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_sv .and. nmod_ifn > 0 ) call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ifn_nucl), 'HIND', &
                                       Unpack ( zins(:, 1), mask = gnegt(:, :, :), field = pins(:, :, :, 1) ) * prhodj(:, :, :) )
  end if
!
!*            compute the heterogeneous nucleation by contact: RVHNCI
!
  if ( nbumod == kmi .and. lbu_enable ) then
    if ( lbudget_th ) call Budget_store_init( tbudgets(NBUDGET_TH), 'HINC', &
                                             Unpack ( zths(:), mask = gnegt(:, :, :), field = pths(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rc ) call Budget_store_init( tbudgets(NBUDGET_RC), 'HINC', prcs(:, :, :) * prhodj(:, :, :) )
    if ( lbudget_ri ) call Budget_store_init( tbudgets(NBUDGET_RI), 'HINC', &
                                             Unpack ( zris(:), mask = gnegt(:, :, :), field = pris(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_sv ) then
      call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nc), 'HINC', pccs(:, :, :) * prhodj(:, :, :) )
      call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ni), 'HINC', &
                                       Unpack ( zcis(:), mask = gnegt(:, :, :), field = pcis(:, :, :) ) * prhodj(:, :, :) )
      if ( nmod_ifn > 0 ) &
        call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ifn_nucl), 'HINC', &
                                 Unpack ( zins(:, 1), mask = gnegt(:, :, :), field = pins(:, :, :, 1) ) * prhodj(:, :, :) )
    end if
  end if

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
  if ( nbumod == kmi .and. lbu_enable ) then
    if ( lbudget_th ) call Budget_store_end( tbudgets(NBUDGET_TH), 'HINC', pths(:, :, :) * prhodj(:, :, :) )
    if ( lbudget_rc ) call Budget_store_end( tbudgets(NBUDGET_RC), 'HINC', prcs(:, :, :) * prhodj(:, :, :) )
    if ( lbudget_ri ) call Budget_store_end( tbudgets(NBUDGET_RI), 'HINC', pris(:, :, :) * prhodj(:, :, :) )
    if ( lbudget_sv ) then
      call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nc), 'HINC', pccs(:, :, :) * prhodj(:, :, :) )
      call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ni), 'HINC', pcis(:, :, :) * prhodj(:, :, :) )
      if ( nmod_ifn > 0 ) &
        call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ifn_nucl), 'HINC', pins(:, :, :, 1) * prhodj(:, :, :) )
    end if
  end if

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
END IF
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_MEYERS
