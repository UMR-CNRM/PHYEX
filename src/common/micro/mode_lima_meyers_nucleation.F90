!MNH_LIC Copyright 2018-2024 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_LIMA_MEYERS_NUCLEATION
  IMPLICIT NONE
CONTAINS
!     #############################################################################
  SUBROUTINE LIMA_MEYERS_NUCLEATION (LIMAP, LIMAC, D, CST, PTSTEP,               &
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
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_T
USE MODD_CST,            ONLY: CST_T

USE MODE_TOOLS,           only: COUNTJV
USE MODD_PARAM_LIMA,      ONLY:PARAM_LIMA_T
USE MODD_PARAM_LIMA_COLD, ONLY:PARAM_LIMA_COLD_T
USE YOMHOOK, ONLY:LHOOK, DR_HOOK, JPHOOK

IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(PARAM_LIMA_T),INTENT(IN)::LIMAP
TYPE(PARAM_LIMA_COLD_T),INTENT(IN)::LIMAC
TYPE(DIMPHYEX_T),         INTENT(IN)    :: D
TYPE(CST_T),              INTENT(IN)    :: CST
REAL,                     INTENT(IN)    :: PTSTEP
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
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PCCT    ! Cloud water C. at t
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PCIT    ! Ice crystal C. source
REAL, DIMENSION(D%NIJT,D%NKT,LIMAP%NMOD_IFN), INTENT(INOUT) :: PINT    ! Activated ice nuclei C.
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(OUT)   :: P_TH_HIND
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(OUT)   :: P_RI_HIND
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(OUT)   :: P_CI_HIND
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(OUT)   :: P_TH_HINC
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(OUT)   :: P_RC_HINC
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(OUT)   :: P_CC_HINC
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PICEFR
!
!
!*       0.2   Declarations of local variables :
!
!
INTEGER :: IIJB, IIJE, IKB, IKE               ! Physical domain
INTEGER :: IL     ! Loop index
INTEGER :: INEGT  ! Case number of nucleation
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
REAL,    DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2))   &
                                  :: ZW, ZT ! work arrays
!
REAL,    DIMENSION(:),   ALLOCATABLE :: ZTCELSIUS
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
!
!*       1.     PRELIMINARY COMPUTATIONS
!               ------------------------
!
IF (LHOOK) CALL DR_HOOK('LIMA_MEYERS_NUCLEATION', 0, ZHOOK_HANDLE)
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
!  optimization by looking for locations where
!  the temperature is negative only !!!
!
GNEGT(:,:) = .FALSE.
GNEGT(D%NIJB:D%NIJE,D%NKB:D%NKE) = ZT(D%NIJB:D%NIJE,D%NKB:D%NKE)<CST%XTT .AND. &
                           ZW(D%NIJB:D%NIJE,D%NKB:D%NKE)>0.8 
INEGT = COUNTJV( GNEGT(:,:),I1(:),I3(:))
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
    ZTHT(IL) = PTHT(I1(IL),I3(IL))
!
    ZCCT(IL) = PCCT(I1(IL),I3(IL))
    ZCIT(IL) = PCIT(I1(IL),I3(IL))
!
    ZRHODREF(IL) = PRHODREF(I1(IL),I3(IL))
    ZZT(IL)      = ZT(I1(IL),I3(IL))
    ZPRES(IL)    = PPABST(I1(IL),I3(IL))
    ZEXNREF(IL)  = PEXNREF(I1(IL),I3(IL))
  ENDDO
  ALLOCATE(ZZW(INEGT))
  ALLOCATE(ZZX(INEGT))
  ALLOCATE(ZZY(INEGT))
  ALLOCATE(ZLSFACT(INEGT))
  ALLOCATE(ZLVFACT(INEGT))
  ALLOCATE(ZSSI(INEGT))
  ALLOCATE(ZTCELSIUS(INEGT))
!
  ZZW(:)  = ZEXNREF(:)*( CST%XCPD+CST%XCPV*ZRVT(:)+CST%XCL*(ZRCT(:)+ZRRT(:)) &
                                  +CST%XCI*(ZRIT(:)+ZRST(:)+ZRGT(:)) )
  ZTCELSIUS(:) = MAX( ZZT(:)-CST%XTT,-50.0 )
  ZLSFACT(:) = (CST%XLSTT+(CST%XCPV-CST%XCI)*ZTCELSIUS(:))/ZZW(:) ! L_s/(Pi_ref*C_ph)
  ZLVFACT(:) = (CST%XLVTT+(CST%XCPV-CST%XCL)*ZTCELSIUS(:))/ZZW(:) ! L_v/(Pi_ref*C_ph)
!
  ZZW(:) = EXP( CST%XALPI - CST%XBETAI/ZZT(:) - CST%XGAMI*ALOG(ZZT(:)) ) ! es_i
  ZSSI(:) = ZRVT(:)*(ZPRES(:)-ZZW(:))/((CST%XMV/CST%XMD)*ZZW(:)) - 1.0
                                                    ! Supersaturation over ice
!
!---------------------------------------------------------------------------
!
!*            compute the heterogeneous nucleation by deposition: RVHNDI
!
  DO IL=1,INEGT
    ZINT(IL,1) = PINT(I1(IL),I3(IL),1)
  END DO
  ZZW(:) = 0.0
  ZZX(:) = 0.0
  ZZY(:) = 0.0
!
  WHERE( ZZT(:)<CST%XTT-5.0 .AND. ZSSI(:)>0.0 )
    ZZY(:) = LIMAC%XNUC_DEP*EXP( LIMAC%XEXSI_DEP*100.*MIN(1.,ZSSI(:))+LIMAC%XEX_DEP)/ZRHODREF(:)
    ZZX(:) = MAX( ZZY(:)-ZINT(:,1) , 0.0 ) ! number of ice crystals formed at this time step #/kg
    ZZW(:) = MIN( LIMAC%XMNU0*ZZX(:) , ZRVT(:) ) ! mass of ice formed at this time step (kg/kg)
  END WHERE
  !
  P_CI_HIND(:,:) = UNPACK( ZZX(:), MASK=GNEGT(:,:), FIELD=0. )
  P_RI_HIND(:,:) = UNPACK( ZZW(:), MASK=GNEGT(:,:), FIELD=0. )
  P_TH_HIND(:,:) = UNPACK( ZZW(:)*(ZLSFACT(:)-ZLVFACT(:)), MASK=GNEGT(:,:), FIELD=0. )
  PTHT(:,:) = PTHT(:,:) + P_TH_HIND(:,:)
  PRVT(:,:) = PRVT(:,:) - P_RI_HIND(:,:)
  PRIT(:,:) = PRIT(:,:) + P_RI_HIND(:,:)
  PCIT(:,:) = PCIT(:,:) + P_CI_HIND(:,:)
  PINT(:,:,1) = PINT(:,:,1) + P_CI_HIND(:,:)
!
!---------------------------------------------------------------------------
!
!*            compute the heterogeneous nucleation by contact: RVHNCI
!
!
  DO IL=1,INEGT
    ZINT(IL,1) = PINT(I1(IL),I3(IL),1)
  END DO
  ZZW(:) = 0.0
  ZZX(:) = 0.0
  ZZY(:) = 0.0
!
  WHERE( ZZT(:)<CST%XTT-2.0 .AND. ZCCT(:)>LIMAP%XCTMIN(2) .AND. ZRCT(:)>LIMAP%XRTMIN(2) )
    ZZY(:) = MIN( LIMAC%XNUC_CON * EXP( LIMAC%XEXTT_CON*ZTCELSIUS(:)+LIMAC%XEX_CON )             &
                                               /ZRHODREF(:) , ZCCT(:) )
    ZZX(:) = MAX( ZZY(:)-ZINT(:,1),0.0 )
    ZZW(:) = MIN( (ZRCT(:)/ZCCT(:))*ZZX(:),ZRCT(:) )
  END WHERE
!
  P_RC_HINC(:,:) = - UNPACK( ZZW(:), MASK=GNEGT(:,:), FIELD=0. )
  P_CC_HINC(:,:) = - UNPACK( ZZX(:), MASK=GNEGT(:,:), FIELD=0. )
  P_TH_HINC(:,:) =   UNPACK( ZZW(:)*(ZLSFACT(:)-ZLVFACT(:)), MASK=GNEGT(:,:), FIELD=0. )
  PTHT(:,:) = PTHT(:,:) + P_TH_HINC(:,:)
  PRCT(:,:) = PRCT(:,:) + P_RC_HINC(:,:)
  PRIT(:,:) = PRIT(:,:) - P_RC_HINC(:,:)
  PCCT(:,:) = PCCT(:,:) + P_CC_HINC(:,:)
  PCIT(:,:) = PCIT(:,:) - P_CC_HINC(:,:)
  PINT(:,:,1) = PINT(:,:,1) - P_CC_HINC(:,:)
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
IF (LHOOK) CALL DR_HOOK('LIMA_MEYERS_NUCLEATION', 1, ZHOOK_HANDLE)
END SUBROUTINE LIMA_MEYERS_NUCLEATION
END MODULE MODE_LIMA_MEYERS_NUCLEATION
