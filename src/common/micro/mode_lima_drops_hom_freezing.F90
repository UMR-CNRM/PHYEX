!MNH_LIC Copyright 2013-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
MODULE MODE_LIMA_DROPS_HOM_FREEZING
  IMPLICIT NONE
CONTAINS
!     ###############################################################################
  SUBROUTINE LIMA_DROPS_HOM_FREEZING (CST, LIMAP, KSIZE, PTSTEP, ODCOMPUTE,     &
                                      PEXNREF, PPABST,                          &
                                      PTHT, PRVT, PRCT, PRRT, PRIT, PRST, PRGT, &
                                      PCRT,                                     &
                                      P_TH_HONR, P_RR_HONR, P_CR_HONR,          &
                                      PB_TH, PB_RR, PB_CR, PB_RG, PB_CG         )
!     ###############################################################################
!
!!    PURPOSE
!!    -------
!!      Homogeneous freezing of rain drops below -35°C
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
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAM_LIMA, ONLY:PARAM_LIMA_T
USE MODD_CST, ONLY:CST_T
USE YOMHOOK, ONLY:LHOOK, DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER,               INTENT(IN)    :: KSIZE
REAL,                  INTENT(IN)    :: PTSTEP 
LOGICAL, DIMENSION(KSIZE), INTENT(IN)    :: ODCOMPUTE
!
REAL, DIMENSION(KSIZE),    INTENT(IN)    :: PEXNREF ! Reference Exner function
REAL, DIMENSION(KSIZE),    INTENT(IN)    :: PPABST  ! abs. pressure at time t
!
REAL, DIMENSION(KSIZE),    INTENT(IN)    :: PTHT    ! Theta at time t
REAL, DIMENSION(KSIZE),    INTENT(IN)    :: PRVT    ! Water vapor m.r. at t 
REAL, DIMENSION(KSIZE),    INTENT(IN)    :: PRCT    ! Cloud water m.r. at t 
REAL, DIMENSION(KSIZE),    INTENT(IN)    :: PRRT    ! Rain water m.r. at t 
REAL, DIMENSION(KSIZE),    INTENT(IN)    :: PRIT    ! Cloud ice m.r. at t 
REAL, DIMENSION(KSIZE),    INTENT(IN)    :: PRST    ! Snow/aggregate m.r. at t 
REAL, DIMENSION(KSIZE),    INTENT(IN)    :: PRGT    ! Graupel m.r. at t 
!
REAL, DIMENSION(KSIZE),    INTENT(IN)    :: PCRT    ! Rain water C. at t
!
REAL, DIMENSION(KSIZE),    INTENT(OUT) :: P_TH_HONR
REAL, DIMENSION(KSIZE),    INTENT(OUT) :: P_RR_HONR
REAL, DIMENSION(KSIZE),    INTENT(OUT) :: P_CR_HONR
REAL, DIMENSION(KSIZE),    INTENT(INOUT) :: PB_TH
REAL, DIMENSION(KSIZE),    INTENT(INOUT) :: PB_RR
REAL, DIMENSION(KSIZE),    INTENT(INOUT) :: PB_CR
REAL, DIMENSION(KSIZE),    INTENT(INOUT) :: PB_RG
REAL, DIMENSION(KSIZE),    INTENT(INOUT) :: PB_CG
!
!*       0.2   Declarations of local variables :
!
TYPE(PARAM_LIMA_T),INTENT(IN)::LIMAP
TYPE(CST_T),INTENT(IN)::CST
REAL, DIMENSION(SIZE(PTHT)) ::  &
     ZW,       &
     ZT,       &
     ZLSFACT,  &
     ZLVFACT,  &
     ZTCELSIUS
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('LIMA_DROPS_HOM_FREEZING', 0, ZHOOK_HANDLE)
P_TH_HONR(:) = 0.
P_RR_HONR(:) = 0.
P_CR_HONR(:) = 0.
!
! Temperature
ZT(:)  = PTHT(:) * ( PPABST(:)/CST%XP00 ) ** (CST%XRD/CST%XCPD)
ZTCELSIUS(:) = ZT(:)-CST%XTT                                    ! T [°C]
!
ZW(:) = PEXNREF(:)*( CST%XCPD+CST%XCPV*PRVT(:)+CST%XCL*(PRCT(:)+PRRT(:)) &
     +CST%XCI*(PRIT(:)+PRST(:)+PRGT(:)) )
ZLSFACT(:) = (CST%XLSTT+(CST%XCPV-CST%XCI)*ZTCELSIUS(:))/ZW(:)          ! L_s/(Pi_ref*C_ph)
ZLVFACT(:) = (CST%XLVTT+(CST%XCPV-CST%XCL)*ZTCELSIUS(:))/ZW(:)          ! L_v/(Pi_ref*C_ph)
!
ZW(:) = 0.0
!
WHERE( (ZT(:)<CST%XTT-35.0) .AND. (PRRT(:)>LIMAP%XRTMIN(3)) .AND. ODCOMPUTE(:) )
   P_TH_HONR(:) = PRRT(:)*(ZLSFACT(:)-ZLVFACT(:))
   P_RR_HONR(:) = - PRRT(:)
   P_CR_HONR(:) = - PCRT(:)
   PB_TH(:) = PB_TH(:) + P_TH_HONR(:)
   PB_RR(:) = PB_RR(:) - PRRT(:)
   PB_CR(:) = PB_CR(:) - PCRT(:)
   PB_RG(:) = PB_RG(:) + PRRT(:)
   PB_CG(:) = PB_CG(:) + PCRT(:)
ENDWHERE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('LIMA_DROPS_HOM_FREEZING', 1, ZHOOK_HANDLE)
END SUBROUTINE LIMA_DROPS_HOM_FREEZING
END MODULE MODE_LIMA_DROPS_HOM_FREEZING
