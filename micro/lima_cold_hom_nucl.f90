!MNH_LIC Copyright 2013-2020 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!      ######################
       MODULE MODI_LIMA_COLD_HOM_NUCL
!      ######################
!
INTERFACE
      SUBROUTINE LIMA_COLD_HOM_NUCL (OHHONI, PTSTEP, KMI,             &
                           PZZ, PRHODJ,                               &
                           PRHODREF, PEXNREF, PPABST, PW_NU,          &
                           PTHT, PRVT, PRCT, PRRT, PRIT, PRST, PRGT,  &
                           PTHS, PRVS, PRCS, PRRS, PRIS, PRGS,        &
                           PCCT,                                      &
                           PCCS, PCRS, PNFS, PCIS, PNHS               )
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
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PW_NU   ! updraft velocity used for
                                                   ! the nucleation param.
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
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRRS    ! Rain water m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRIS    ! Pristine ice m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRGS    ! Graupel/hail m.r. source
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCCT    ! Cloud water C. at t
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCCS    ! Cloud water C. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCRS    ! Rain water C. source
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PNFS    ! CCN C. available source
                                                   !used as Free ice nuclei for
                                                   !HOMOGENEOUS nucleation of haze
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCIS    ! Ice crystal C. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PNHS    ! haze homogeneous freezing
!
END SUBROUTINE LIMA_COLD_HOM_NUCL
END INTERFACE
END MODULE MODI_LIMA_COLD_HOM_NUCL
!
!     ######################################################################
      SUBROUTINE LIMA_COLD_HOM_NUCL (OHHONI, PTSTEP, KMI,             &
                           PZZ, PRHODJ,                               &
                           PRHODREF, PEXNREF, PPABST, PW_NU,          &
                           PTHT, PRVT, PRCT, PRRT, PRIT, PRST, PRGT,  &
                           PTHS, PRVS, PRCS, PRRS, PRIS, PRGS,        &
                           PCCT,                                      &
                           PCCS, PCRS, PNFS, PCIS, PNHS               )
!     ######################################################################
!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the cold-phase homogeneous
!!    freezing of CCN, droplets and drops (T<-35°C)
!!
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
!!      C. Barthe  * LACy*   jan. 2014  add budgets
!!      B.Vie 10/2016 Bug zero division
!  P. Wautelet 05/2016-04/2018: new data structures and calls for I/O
!  M. Leriche     05/2019: suppress unused actived coated IN (immersion) source 
!  P. Wautelet 28/05/2019: move COUNTJV function to tools.f90
!  P. Wautelet    03/2020: use the new data structures and subroutines for budgets
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
use modd_budget,          only: lbu_enable, nbumod,                                                                  &
                                lbudget_th, lbudget_rv, lbudget_rc, lbudget_rr, lbudget_ri, lbudget_rg, lbudget_sv,  &
                                NBUDGET_TH, NBUDGET_RV, NBUDGET_RC, NBUDGET_RR, NBUDGET_RI, NBUDGET_RG, NBUDGET_SV1, &
                                tbudgets
USE MODD_CST,             ONLY: XP00, XRD, XRV, XMV, XMD, XCPD, XCPV, XCL, XCI,   &
                                XTT, XLSTT, XLVTT, XALPI, XBETAI, XGAMI,          &
                                XG
USE MODD_NSV
USE MODD_PARAMETERS,      ONLY: JPHEXT, JPVEXT
USE MODD_PARAM_LIMA,      ONLY: NMOD_CCN, NMOD_IMM, XRTMIN, XCTMIN, XNUC, NMOM_C, NMOM_R
USE MODD_PARAM_LIMA_COLD, ONLY: XRCOEF_HONH, XCEXP_DIFVAP_HONH, XCOEF_DIFVAP_HONH,&
                                XCRITSAT1_HONH, XCRITSAT2_HONH, XTMAX_HONH,       &
                                XTMIN_HONH, XC1_HONH, XC2_HONH, XC3_HONH,         &
                                XDLNJODT1_HONH, XDLNJODT2_HONH, XRHOI_HONH,       &
                                XC_HONC, XTEXP1_HONC, XTEXP2_HONC, XTEXP3_HONC,   &
                                XTEXP4_HONC, XTEXP5_HONC
USE MODD_PARAM_LIMA_WARM, ONLY: XLBC

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
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PW_NU   ! updraft velocity used for
                                                   ! the nucleation param.
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
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRRS    ! Rain water m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRIS    ! Pristine ice m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRGS    ! Graupel/hail m.r. source
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCCT    ! Cloud water C. at t
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCCS    ! Cloud water C. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCRS    ! Rain water C. source
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PNFS    ! CCN C. available source
                                                   !used as Free ice nuclei for
                                                   !HOMOGENEOUS nucleation of haze
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCIS    ! Ice crystal C. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PNHS    ! haze homogeneous freezing
!
!*       0.2   Declarations of local variables :
!
REAL, DIMENSION(:), ALLOCATABLE :: ZRVT    ! Water vapor m.r. at t
REAL, DIMENSION(:), ALLOCATABLE :: ZRCT    ! Cloud water m.r. at t
REAL, DIMENSION(:), ALLOCATABLE :: ZRRT    ! Rain water m.r. at t
REAL, DIMENSION(:), ALLOCATABLE :: ZRIT    ! Pristine ice m.r. at t
REAL, DIMENSION(:), ALLOCATABLE :: ZRST    ! Snow/aggregate m.r. at t
REAL, DIMENSION(:), ALLOCATABLE :: ZRGT    ! Graupel/hail m.r. at t
!
REAL, DIMENSION(:), ALLOCATABLE :: ZCCT    ! Cloud water conc. at t
!
REAL, DIMENSION(:), ALLOCATABLE :: ZTHS    ! Theta source
!
REAL, DIMENSION(:), ALLOCATABLE :: ZRVS    ! Water vapor m.r. source
REAL, DIMENSION(:), ALLOCATABLE :: ZRCS    ! Cloud water m.r. source
REAL, DIMENSION(:), ALLOCATABLE :: ZRRS    ! Rain water m.r. source
REAL, DIMENSION(:), ALLOCATABLE :: ZRIS    ! Pristine ice m.r. source
REAL, DIMENSION(:), ALLOCATABLE :: ZRGS    ! Graupel/hail m.r. source
!
REAL, DIMENSION(:),   ALLOCATABLE :: ZCCS    ! Cloud water conc. source
REAL, DIMENSION(:),   ALLOCATABLE :: ZCRS    ! Rain water conc. source
REAL, DIMENSION(:,:), ALLOCATABLE :: ZNFS    ! available nucleus conc. source
REAL, DIMENSION(:),   ALLOCATABLE :: ZCIS    ! Pristine ice conc. source
REAL, DIMENSION(:),   ALLOCATABLE :: ZZNHS   ! Nucleated Ice nuclei conc. source
                                             !by Homogeneous freezing
!
REAL, DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3))   &
                                  :: ZNHS  ! Nucleated Ice nuclei conc. source
                                           ! by Homogeneous freezing of haze
REAL, DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3))   &
                                  :: ZW, ZT ! work arrays
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
                              ZSI,      & ! Saturation over ice
                              ZTCELSIUS,&
                              ZLS,      &
                              ZPSI1,    &
                              ZPSI2,    &
                              ZTAU,     &
                              ZBFACT,   &
                              ZW_NU,    &
                              ZFREECCN, &
                              ZCCNFROZEN
!
INTEGER :: IIB, IIE, IJB, IJE, IKB, IKE   ! Physical domain
INTEGER :: JL, JMOD_CCN, JMOD_IMM         ! Loop index
!
INTEGER :: INEGT                          ! Case number of hom. nucleation
integer :: idx
LOGICAL, DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3)) &
			  :: GNEGT        ! Test where to compute the hom. nucleation
INTEGER , DIMENSION(SIZE(GNEGT)) :: I1,I2,I3 ! Used to replace the COUNT
!
REAL    :: ZEPS                           ! molar mass ratio
!
!-------------------------------------------------------------------------------
!
!
!*       1.     PRELIMINARY COMPUTATIONS
!	        ------------------------
!
!
! Physical domain
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
IF( OHHONI ) THEN
   ZNHS(:,:,:) = PNHS(:,:,:)
ELSE
   ZNHS(:,:,:) = 0.0 
END IF
!
! Computations only where the temperature is below -35°C
! PACK variables
!
GNEGT(:,:,:) = .FALSE.
GNEGT(IIB:IIE,IJB:IJE,IKB:IKE) = ZT(IIB:IIE,IJB:IJE,IKB:IKE)<XTT-35.0
INEGT = COUNTJV( GNEGT(:,:,:),I1(:),I2(:),I3(:))
!
IF (INEGT.GT.0) THEN

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
   ALLOCATE(ZRRS(INEGT))
   ALLOCATE(ZRIS(INEGT))
   ALLOCATE(ZRGS(INEGT))
   !
   ALLOCATE(ZTHS(INEGT))
   !
   ALLOCATE(ZCCS(INEGT))
   ALLOCATE(ZCRS(INEGT))
   ALLOCATE(ZCIS(INEGT))
   !
   ALLOCATE(ZNFS(INEGT,NMOD_CCN))
   ALLOCATE(ZZNHS(INEGT))
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
      ZRVS(JL) = PRVS(I1(JL),I2(JL),I3(JL))
      ZRCS(JL) = PRCS(I1(JL),I2(JL),I3(JL))
      ZRRS(JL) = PRRS(I1(JL),I2(JL),I3(JL))
      ZRIS(JL) = PRIS(I1(JL),I2(JL),I3(JL))
      ZRGS(JL) = PRGS(I1(JL),I2(JL),I3(JL))
      !
      ZTHS(JL) = PTHS(I1(JL),I2(JL),I3(JL))
      !
      ZCCS(JL) = PCCS(I1(JL),I2(JL),I3(JL))
      ZCRS(JL) = PCRS(I1(JL),I2(JL),I3(JL))
      ZCIS(JL) = PCIS(I1(JL),I2(JL),I3(JL))
      !
      DO JMOD_CCN = 1, NMOD_CCN
         ZNFS(JL,JMOD_CCN) = PNFS(I1(JL),I2(JL),I3(JL),JMOD_CCN)
      ENDDO
      ZZNHS(JL) = ZNHS(I1(JL),I2(JL),I3(JL))
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
   ALLOCATE( ZLBDAC     (INEGT) )
!
   ALLOCATE( ZZW (INEGT) ) ; ZZW(:) = 0.0
   ALLOCATE( ZZX (INEGT) ) ; ZZX(:) = 0.0
   ALLOCATE( ZZY (INEGT) ) ; ZZY(:) = 0.0
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
!
!-------------------------------------------------------------------------------
!
!
!*       2.     Haze homogeneous freezing
!	        ------------------------
!
  if ( nbumod == kmi .and. lbu_enable .and. ohhoni .and. nmod_ccn > 0 ) then
    if ( lbudget_th ) call Budget_store_init( tbudgets(NBUDGET_TH), 'HONH', pths(:, :, :) * prhodj(:, :, :) )
    if ( lbudget_rv ) call Budget_store_init( tbudgets(NBUDGET_RV), 'HONH', prvs(:, :, :) * prhodj(:, :, :) )
    if ( lbudget_ri ) call Budget_store_init( tbudgets(NBUDGET_RI), 'HONH', pris(:, :, :) * prhodj(:, :, :) )
    if ( lbudget_sv ) then
      call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ni), 'HONH', pcis(:, :, :) * prhodj(:, :, :) )
      if ( nmod_ccn > 0 ) then
        do jl = 1, nmod_ccn
          idx = NBUDGET_SV1 - 1 + nsv_lima_ccn_free - 1 + jl
          call Budget_store_init( tbudgets(idx), 'HONH', pnfs(:, :, :, jl) * prhodj(:, :, :) )
        end do
        call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_hom_haze), 'HONH', znhs(:, :, :) * prhodj(:, :, :) )
      end if
    end if
  end if
!
!  Compute the haze homogeneous nucleation source: RHHONI
!
   IF( OHHONI .AND. NMOD_CCN.GT.0 ) THEN

! Sum of the available CCN
      ALLOCATE( ZFREECCN(INEGT) )
      ALLOCATE( ZCCNFROZEN(INEGT) )
      ZFREECCN(:)=0.
      ZCCNFROZEN(:)=0.
      DO JMOD_CCN = 1, NMOD_CCN
         ZFREECCN(:) = ZFREECCN(:) + ZNFS(:,JMOD_CCN)
      END DO
!
      ALLOCATE(ZW_NU(INEGT))
      DO JL=1,INEGT
         ZW_NU(JL) = PW_NU(I1(JL),I2(JL),I3(JL))
      END DO
!
      ZZW(:)  = 0.0
      ZZX(:)  = 0.0
      ZEPS    = XMV / XMD
      ZZY(:)  = XCRITSAT1_HONH -                              &  ! Critical Sat.
              (MIN( XTMAX_HONH,MAX( XTMIN_HONH,ZZT(:) ) )/XCRITSAT2_HONH)
!
      ALLOCATE(ZLS(INEGT))
      ALLOCATE(ZPSI1(INEGT))
      ALLOCATE(ZPSI2(INEGT))
      ALLOCATE(ZTAU(INEGT))
      ALLOCATE(ZBFACT(INEGT))
!
      WHERE( (ZZT(:)<XTT-35.0) .AND. (ZSI(:)>ZZY(:)) )
            ZLS(:)   = XLSTT+(XCPV-XCI)*ZTCELSIUS(:)          ! Ls
!
            ZPSI1(:) = ZZY(:) * (XG/(XRD*ZZT(:)))*(ZEPS*ZLS(:)/(XCPD*ZZT(:))-1.)
!                                                         ! Psi1 (a1*Scr in KL01)
! BV correction PSI2 enlever 1/ZEPS ?
!            ZPSI2(:) = ZSI(:) * (1.0/ZEPS+1.0/ZRVT(:)) +                           &
            ZPSI2(:) = ZSI(:) * (1.0/ZRVT(:)) +                           &
                 ZZY(:) * ((ZLS(:)/ZZT(:))**2)/(XCPD*XRV) 
!                                                         ! Psi2 (a2+a3*Scr in KL01)
            ZTAU(:) = 1.0 / ( MAX( XC1_HONH,XC1_HONH*(XC2_HONH-XC3_HONH*ZZT(:)) ) *&
                 ABS( (XDLNJODT1_HONH - XDLNJODT2_HONH*ZZT(:))       *             &
                 ((ZPRES(:)/XP00)**(XRD/XCPD))*ZTHS(:) ) )
!
            ZBFACT(:) = (XRHOI_HONH/ZRHODREF(:)) * (ZSI(:)/(ZZY(:)-1.0))           &
! BV correction ZBFACT enlever 1/ZEPS ?
!                 * (1.0/ZRVT(:)+1.0/ZEPS)                                          &
                 * (1.0/ZRVT(:))                                          &
                 / (XCOEF_DIFVAP_HONH*(ZZT(:)**XCEXP_DIFVAP_HONH /ZPRES(:)))
!
! BV correction ZZX rho_i{-1} ?
!            ZZX(:) = MAX( MIN( XRHOI_HONH*ZBFACT(:)**1.5 * (ZPSI1(:)/ZPSI2(:))     &
            ZZX(:) = MAX( MIN( (1/XRHOI_HONH)*ZBFACT(:)**1.5 * (ZPSI1(:)/ZPSI2(:))     &
! BV correction ZZX PTSTEP wrong place ?
!                 * (ZW_NU(:)/SQRT(ZTAU(:))), ZNFS(:,JMOD_CCN) )/PTSTEP , 0.)
                 * (ZW_NU(:)/SQRT(ZTAU(:)))/PTSTEP , ZFREECCN(:) ) , 0.)
!
            ZZW(:) = MIN( XRCOEF_HONH*ZZX(:)*(ZTAU(:)/ZBFACT(:))**1.5 , ZRVS(:) )
      END WHERE
!
! Apply the changes to ZNFS,
    DO JMOD_CCN = 1, NMOD_CCN
        WHERE(ZFREECCN(:)>1.)
            ZCCNFROZEN(:) = ZZX(:) * ZNFS(:,JMOD_CCN)/ZFREECCN(:)
            ZNFS(:,JMOD_CCN) = ZNFS(:,JMOD_CCN) - ZCCNFROZEN(:)
        END WHERE
        ZW(:,:,:) = PNFS(:,:,:,JMOD_CCN)
        PNFS(:,:,:,JMOD_CCN)=UNPACK( ZNFS(:,JMOD_CCN), MASK=GNEGT(:,:,:),FIELD=ZW(:,:,:))
    END DO
      ZZNHS(:)    = ZZNHS(:) + ZZX(:)
      ZNHS(:,:,:) = UNPACK( ZZNHS(:), MASK=GNEGT(:,:,:),FIELD=0.0)
      PNHS(:,:,:) = ZNHS(:,:,:)
!
      DEALLOCATE(ZFREECCN)
      DEALLOCATE(ZCCNFROZEN)
      DEALLOCATE(ZLS)
      DEALLOCATE(ZPSI1)
      DEALLOCATE(ZPSI2)
      DEALLOCATE(ZTAU)
      DEALLOCATE(ZBFACT)
      DEALLOCATE(ZW_NU)
!
      ZRVS(:) = ZRVS(:) - ZZW(:)
      ZRIS(:) = ZRIS(:) + ZZW(:)
      ZTHS(:) = ZTHS(:) + ZZW(:) * (ZLSFACT(:)-ZLVFACT(:)) ! f(L_s*(RHHONI))
      ZCIS(:) = ZCIS(:) + ZZX(:)
!
   END IF ! OHHONI
!
! Budget storage
  if ( nbumod == kmi .and. lbu_enable .and. ohhoni .and. nmod_ccn > 0 ) then
    if ( lbudget_th ) call Budget_store_end( tbudgets(NBUDGET_TH), 'HONH', &
                                             Unpack( zths(:), mask = gnegt(:, :, :), field = pths(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rv ) call Budget_store_end( tbudgets(NBUDGET_RV), 'HONH', &
                                             Unpack( zrvs(:), mask = gnegt(:, :, :), field = prvs(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_ri ) call Budget_store_end( tbudgets(NBUDGET_RI), 'HONH', &
                                             Unpack( zris(:), mask = gnegt(:, :, :), field = pris(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_sv ) then
      call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ni), 'HONH', &
                             Unpack( zcis(:), mask = gnegt(:, :, :), field = pcis(:, :, :) ) * prhodj(:, :, :) )
      if ( nmod_ccn > 0 ) then
        do jl = 1, nmod_ccn
          idx = NBUDGET_SV1 - 1 + nsv_lima_ccn_free - 1 + jl
          call Budget_store_end( tbudgets(idx), 'HONH', &
                                 Unpack( znfs(:, jl), mask = gnegt(:, :, :), field = pnfs(:, :, :, jl) ) * prhodj(:, :, :) )
        end do
        call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_hom_haze), 'HONH', &
                               Unpack( zznhs(:), mask = gnegt(:, :, :), field = znhs(:, :, :) ) * prhodj(:, :, :) )
      end if
    end if
  end if
!
!-------------------------------------------------------------------------------
!
!
!*       3.     Cloud droplets homogeneous freezing
!	        -----------------------------------
!
!
!  Compute the droplet homogeneous nucleation source: RCHONI
!                 -> Pruppacher(1995)
!
IF (NMOM_C.GE.2) THEN
  if ( nbumod == kmi .and. lbu_enable ) then
    if ( lbudget_th ) call Budget_store_init( tbudgets(NBUDGET_TH), 'HONC', &
                                             Unpack( zths(:), mask = gnegt(:, :, :), field = pths(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rc ) call Budget_store_init( tbudgets(NBUDGET_RC), 'HONC', prcs(:, :, :) * prhodj(:, :, :) )
    if ( lbudget_ri ) call Budget_store_init( tbudgets(NBUDGET_RI), 'HONC', &
                                             Unpack( zris(:), mask = gnegt(:, :, :), field = pris(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_sv ) then
      call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nc), 'HONC', pccs(:, :, :) * prhodj(:, :, :) )
      call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ni), 'HONC', &
                             Unpack( zcis(:), mask = gnegt(:, :, :), field = pcis(:, :, :) ) * prhodj(:, :, :) )
    end if
  end if

   ZZW(:) = 0.0
   ZZX(:) = 0.0
   WHERE( (ZZT(:)<XTT-35.0) .AND. (ZCCT(:)>XCTMIN(2)) .AND. (ZRCT(:)>XRTMIN(2)) )
      ZLBDAC(:) = XLBC*ZCCT(:) / (ZRCT(:)) ! Lambda_c**3
      ZZX(:) = 1.0 / ( 1.0 + (XC_HONC/ZLBDAC(:))*PTSTEP*                      &
           EXP( XTEXP1_HONC + ZTCELSIUS(:)*(                        &
           XTEXP2_HONC + ZTCELSIUS(:)*(                        &
           XTEXP3_HONC + ZTCELSIUS(:)*(                        &
           XTEXP4_HONC + ZTCELSIUS(:)*XTEXP5_HONC))) ) )**XNUC
      ZZW(:) = ZCCS(:) * (1.0 - ZZX(:))                                  ! CCHONI
!                                                                       
      ZCCS(:) = ZCCS(:) - ZZW(:)
      ZCIS(:) = ZCIS(:) + ZZW(:)
!
      ZZW(:) = ZRCS(:) * (1.0 - ZZX(:))                                  ! RCHONI
!                                                                       
      ZRCS(:) = ZRCS(:) - ZZW(:)
      ZRIS(:) = ZRIS(:) + ZZW(:)
      ZTHS(:) = ZTHS(:) + ZZW(:) * (ZLSFACT(:)-ZLVFACT(:)) ! f(L_f*(RCHONI))
   END WHERE
!
! Budget storage
  if ( nbumod == kmi .and. lbu_enable ) then
    if ( lbudget_th ) call Budget_store_end( tbudgets(NBUDGET_TH), 'HONC', &
                                             Unpack( zths(:), mask = gnegt(:, :, :), field = pths(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rc ) call Budget_store_end( tbudgets(NBUDGET_RC), 'HONC', &
                                             Unpack( zrcs(:), mask = gnegt(:, :, :), field = prcs(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_ri ) call Budget_store_end( tbudgets(NBUDGET_RI), 'HONC', &
                                             Unpack( zris(:), mask = gnegt(:, :, :), field = pris(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_sv ) then
      call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nc), 'HONC', &
                             Unpack( zccs(:), mask = gnegt(:, :, :), field = pccs(:, :, :) ) * prhodj(:, :, :) )
      call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ni), 'HONC', &
                             Unpack( zcis(:), mask = gnegt(:, :, :), field = pcis(:, :, :) ) * prhodj(:, :, :) )
    end if
  end if
END IF
!
!
!-------------------------------------------------------------------------------
!
!
!*       4.     Rain drops homogeneous freezing
!	        -------------------------------
!
!
!  Compute the drop homogeneous nucleation source: RRHONG
!
IF (NMOM_R.GE.2) THEN
  if ( nbumod == kmi .and. lbu_enable ) then
    if ( lbudget_th ) call Budget_store_init( tbudgets(NBUDGET_TH), 'HONR', &
                                             Unpack( zths(:), mask = gnegt(:, :, :), field = pths(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rr ) call Budget_store_init( tbudgets(NBUDGET_RR), 'HONR', prrs(:, :, :) * prhodj(:, :, :) )
    if ( lbudget_rg ) call Budget_store_init( tbudgets(NBUDGET_RG), 'HONR', prgs(:, :, :) * prhodj(:, :, :) )
    if ( lbudget_sv ) then
      call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nr), 'HONR', pcrs(:, :, :) * prhodj(:, :, :) )
    end if
  end if

   ZZW(:) = 0.0
   WHERE( (ZZT(:)<XTT-35.0) .AND. (ZRRS(:)>XRTMIN(3)/PTSTEP) )
      ZZW(:)  = ZRRS(:) ! Instantaneous freezing of the raindrops
      ZRRS(:) = ZRRS(:) - ZZW(:)
      ZRGS(:) = ZRGS(:) + ZZW(:)
      ZTHS(:) = ZTHS(:) + ZZW(:)*(ZLSFACT(:)-ZLVFACT(:)) ! f(L_f*(RRHONG))
!
      ZCRS(:) = 0.0     ! No more raindrops when T<-35 C
   ENDWHERE
!
! Budget storage
  if ( nbumod == kmi .and. lbu_enable ) then
    if ( lbudget_th ) call Budget_store_end( tbudgets(NBUDGET_TH), 'HONR', &
                                             Unpack( zths(:), mask = gnegt(:, :, :), field = pths(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rr ) call Budget_store_end( tbudgets(NBUDGET_RR), 'HONR', &
                                             Unpack( zrrs(:), mask = gnegt(:, :, :), field = prrs(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rg ) call Budget_store_end( tbudgets(NBUDGET_RG), 'HONR', &
                                             Unpack( zrgs(:), mask = gnegt(:, :, :), field = prgs(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_sv ) then
      call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nr), 'HONR', &
                             Unpack( zcrs(:), mask = gnegt(:, :, :), field = pcrs(:, :, :) ) * prhodj(:, :, :) )
    end if
  end if
END IF
!
!
!-------------------------------------------------------------------------------
!
!
!*       4.     Unpack variables, clean
!	        -----------------------
!
!
! End of homogeneous nucleation processes
!
   ZW(:,:,:)   = PRVS(:,:,:)
   PRVS(:,:,:) = UNPACK( ZRVS(:),MASK=GNEGT(:,:,:),FIELD=ZW(:,:,:) )
   ZW(:,:,:)   = PRCS(:,:,:)
   PRCS(:,:,:) = UNPACK( ZRCS(:),MASK=GNEGT(:,:,:),FIELD=ZW(:,:,:) )
   ZW(:,:,:)   = PRRS(:,:,:)
   PRRS(:,:,:) = UNPACK( ZRRS(:),MASK=GNEGT(:,:,:),FIELD=ZW(:,:,:) )
   ZW(:,:,:)   = PRIS(:,:,:)
   PRIS(:,:,:) = UNPACK( ZRIS(:),MASK=GNEGT(:,:,:),FIELD=ZW(:,:,:) )
   ZW(:,:,:)   = PRGS(:,:,:)
   PRGS(:,:,:) = UNPACK( ZRGS(:),MASK=GNEGT(:,:,:),FIELD=ZW(:,:,:) )
   ZW(:,:,:)   = PTHS(:,:,:)
   PTHS(:,:,:) = UNPACK( ZTHS(:),MASK=GNEGT(:,:,:),FIELD=ZW(:,:,:) )
   ZW(:,:,:)   = PCCS(:,:,:)
   PCCS(:,:,:) = UNPACK( ZCCS(:),MASK=GNEGT(:,:,:),FIELD=ZW(:,:,:) )
   ZW(:,:,:)   = PCRS(:,:,:)
   PCRS(:,:,:) = UNPACK( ZCRS(:),MASK=GNEGT(:,:,:),FIELD=ZW(:,:,:) )
   ZW(:,:,:)   = PCIS(:,:,:)
   PCIS(:,:,:) = UNPACK( ZCIS(:),MASK=GNEGT(:,:,:),FIELD=ZW(:,:,:) )
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
   DEALLOCATE(ZRRS)
   DEALLOCATE(ZRIS)
   DEALLOCATE(ZRGS)
!
   DEALLOCATE(ZTHS)
!
   DEALLOCATE(ZCCS)
   DEALLOCATE(ZCRS)
   DEALLOCATE(ZCIS)
!
   DEALLOCATE(ZNFS)
   DEALLOCATE(ZZNHS)
!
   DEALLOCATE(ZRHODREF) 
   DEALLOCATE(ZZT) 
   DEALLOCATE(ZPRES) 
   DEALLOCATE(ZEXNREF)
!
   DEALLOCATE(ZLSFACT)
   DEALLOCATE(ZLVFACT)
   DEALLOCATE(ZSI)
   DEALLOCATE(ZTCELSIUS)
   DEALLOCATE(ZLBDAC)
!
   DEALLOCATE(ZZW) 
   DEALLOCATE(ZZX)
   DEALLOCATE(ZZY)
!
END IF ! INEGT>0
!
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_COLD_HOM_NUCL
