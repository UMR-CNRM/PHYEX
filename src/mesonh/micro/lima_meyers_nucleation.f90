!MNH_LIC Copyright 2018-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!      ##################################
       MODULE MODI_LIMA_MEYERS_NUCLEATION
!      ##################################
!
INTERFACE
   SUBROUTINE LIMA_MEYERS_NUCLEATION (PTSTEP,                                     &
                                      PRHODREF, PEXNREF, PPABST,                  &
                                      PTHT, PRVT, PRCT, PRRT, PRIT, PRST, PRGT,   &
                                      PCCT, PCIT, PINT,                           &
                                      P_TH_HIND, P_RI_HIND, P_CI_HIND,            &
                                      P_TH_HINC, P_RC_HINC, P_CC_HINC,            &
                                      PICEFR                                      )
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
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCCT    ! Cloud water C. at t
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCIT    ! Ice crystal C. source
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PINT    ! Activated ice nuclei C.
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
END SUBROUTINE LIMA_MEYERS_NUCLEATION
END INTERFACE
END MODULE MODI_LIMA_MEYERS_NUCLEATION
!
!     #############################################################################
   SUBROUTINE LIMA_MEYERS_NUCLEATION (PTSTEP,                                     &
                                      PRHODREF, PEXNREF, PPABST,                  &
                                      PTHT, PRVT, PRCT, PRRT, PRIT, PRST, PRGT,   &
                                      PCCT, PCIT, PINT,                           &
                                      P_TH_HIND, P_RI_HIND, P_CI_HIND,            &
                                      P_TH_HINC, P_RC_HINC, P_CC_HINC,            &
                                      PICEFR                                      )
!     #############################################################################
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the heterogeneous nucleation
!!    following Meyers (1992).
!!
!!
!!    AUTHOR
!!    ------
!!      J.-M. Cohard     * Laboratoire d'Aerologie*
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!      S.    Berthet    * Laboratoire d'Aerologie*
!!      B.    Vié        * CNRM *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original             15/03/2018
!  P. Wautelet 28/05/2019: move COUNTJV function to tools.f90
!  P. Wautelet 27/02/2020: add P_TH_HINC dummy argument + change intent of *_HIND and *_HINC dummy arguments (INOUT->OUT)
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
USE MODD_NSV,             ONLY: NSV_LIMA_NC, NSV_LIMA_NI
USE MODD_PARAMETERS
USE MODD_PARAM_LIMA
USE MODD_PARAM_LIMA_COLD

use mode_tools,           only: Countjv

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
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCCT    ! Cloud water C. at t
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCIT    ! Ice crystal C. source
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PINT    ! Activated ice nuclei C.
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
REAL, DIMENSION(:),   ALLOCATABLE :: ZCIT    ! Pristine ice conc. source
!
REAL, DIMENSION(:),   ALLOCATABLE :: ZTHT    ! Theta source
!
REAL, DIMENSION(:,:), ALLOCATABLE :: ZINT    ! Nucleated Ice nuclei conc. source
                                             ! by Deposition/Contact
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
  ALLOCATE(ZTHT(INEGT))
!
  ALLOCATE(ZCCT(INEGT))
  ALLOCATE(ZINT(INEGT,1))
  ALLOCATE(ZCIT(INEGT))
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
    ZTHT(JL) = PTHT(I1(JL),I2(JL),I3(JL))
!
    ZCCT(JL) = PCCT(I1(JL),I2(JL),I3(JL))
    ZCIT(JL) = PCIT(I1(JL),I2(JL),I3(JL))
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
!---------------------------------------------------------------------------
!
!*            compute the heterogeneous nucleation by deposition: RVHNDI
!
  DO JL=1,INEGT
    ZINT(JL,1) = PINT(I1(JL),I2(JL),I3(JL),1)
  END DO
  ZZW(:) = 0.0
  ZZX(:) = 0.0
  ZZY(:) = 0.0
!
  WHERE( ZZT(:)<XTT-5.0 .AND. ZSSI(:)>0.0 )
    ZZY(:) = XNUC_DEP*EXP( XEXSI_DEP*100.*MIN(1.,ZSSI(:))+XEX_DEP)/ZRHODREF(:)
    ZZX(:) = MAX( ZZY(:)-ZINT(:,1) , 0.0 ) ! number of ice crystals formed at this time step #/kg
    ZZW(:) = MIN( XMNU0*ZZX(:) , ZRVT(:) ) ! mass of ice formed at this time step (kg/kg)
  END WHERE
  !
  P_CI_HIND(:,:,:) = UNPACK( ZZX(:), MASK=GNEGT(:,:,:), FIELD=0. )
  P_RI_HIND(:,:,:) = UNPACK( ZZW(:), MASK=GNEGT(:,:,:), FIELD=0. )
  P_TH_HIND(:,:,:) = UNPACK( ZZW(:)*(ZLSFACT(:)-ZLVFACT(:)), MASK=GNEGT(:,:,:), FIELD=0. )
  PTHT(:,:,:) = PTHT(:,:,:) + P_TH_HIND(:,:,:)
  PRVT(:,:,:) = PRVT(:,:,:) - P_RI_HIND(:,:,:)
  PRIT(:,:,:) = PRIT(:,:,:) + P_RI_HIND(:,:,:)
  PCIT(:,:,:) = PCIT(:,:,:) + P_CI_HIND(:,:,:)
  PINT(:,:,:,1) = PINT(:,:,:,1) + P_CI_HIND(:,:,:)
!
!---------------------------------------------------------------------------
!
!*            compute the heterogeneous nucleation by contact: RVHNCI
!
!
  DO JL=1,INEGT
    ZINT(JL,1) = PINT(I1(JL),I2(JL),I3(JL),1)
  END DO
  ZZW(:) = 0.0
  ZZX(:) = 0.0
  ZZY(:) = 0.0
!
  WHERE( ZZT(:)<XTT-2.0 .AND. ZCCT(:)>XCTMIN(2) .AND. ZRCT(:)>XRTMIN(2) )
    ZZY(:) = MIN( XNUC_CON * EXP( XEXTT_CON*ZTCELSIUS(:)+XEX_CON )             &
                                               /ZRHODREF(:) , ZCCT(:) )
    ZZX(:) = MAX( ZZY(:)-ZINT(:,1),0.0 )
    ZZW(:) = MIN( (ZRCT(:)/ZCCT(:))*ZZX(:),ZRCT(:) )
  END WHERE
!
  P_RC_HINC(:,:,:) = - UNPACK( ZZW(:), MASK=GNEGT(:,:,:), FIELD=0. )
  P_CC_HINC(:,:,:) = - UNPACK( ZZX(:), MASK=GNEGT(:,:,:), FIELD=0. )
  P_TH_HINC(:,:,:) =   UNPACK( ZZW(:)*(ZLSFACT(:)-ZLVFACT(:)), MASK=GNEGT(:,:,:), FIELD=0. )
  PTHT(:,:,:) = PTHT(:,:,:) + P_TH_HINC(:,:,:)
  PRCT(:,:,:) = PRCT(:,:,:) + P_RC_HINC(:,:,:)
  PRIT(:,:,:) = PRIT(:,:,:) - P_RC_HINC(:,:,:)
  PCCT(:,:,:) = PCCT(:,:,:) + P_CC_HINC(:,:,:)
  PCIT(:,:,:) = PCIT(:,:,:) - P_CC_HINC(:,:,:)
  PINT(:,:,:,1) = PINT(:,:,:,1) - P_CC_HINC(:,:,:)
!
  DEALLOCATE(ZRVT) 
  DEALLOCATE(ZRCT) 
  DEALLOCATE(ZRRT) 
  DEALLOCATE(ZRIT) 
  DEALLOCATE(ZRST) 
  DEALLOCATE(ZRGT) 
!
  DEALLOCATE(ZTHT)
!
  DEALLOCATE(ZCCT)
  DEALLOCATE(ZINT)
  DEALLOCATE(ZCIT)
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
END SUBROUTINE LIMA_MEYERS_NUCLEATION
