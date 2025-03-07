!MNH_LIC Copyright 2018-2024 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-------------------------------------------------------------------------------
MODULE MODE_LIMA_GRAUPEL_DEPOSITION
  IMPLICIT NONE
CONTAINS
!     ###########################################################################
  SUBROUTINE LIMA_GRAUPEL_DEPOSITION (LIMAP, LIMAM, KSIZE, ODCOMPUTE, PRHODREF,                 &
                                      PRGT, PCGT, PSSI, PLBDG, PAI, PCJ, PLSFACT, &
                                      P_TH_DEPG, P_RG_DEPG                        )
!     ###########################################################################
!
!!    PURPOSE
!!    -------
!!      Deposition of water vapour on graupel
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
!       M. Taufour              07/2022 add concentration for snow, graupel, hail        
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAM_LIMA_MIXED, ONLY:PARAM_LIMA_MIXED_T
USE MODD_PARAM_LIMA, ONLY:PARAM_LIMA_T
USE YOMHOOK, ONLY:LHOOK, DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER, INTENT(IN) :: KSIZE
LOGICAL, DIMENSION(KSIZE),INTENT(IN)    :: ODCOMPUTE
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRHODREF ! 
!
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRGT     ! graupel mr
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PCGT     ! graupel conc
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PSSI     ! 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLBDG    ! 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PAI      ! 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PCJ      ! 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLSFACT  ! 
!
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_TH_DEPG
TYPE(PARAM_LIMA_MIXED_T),INTENT(IN)::LIMAM
TYPE(PARAM_LIMA_T),INTENT(IN)::LIMAP
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_RG_DEPG
!
!
!*       0.2   Declarations of local variables :
!
REAL,    DIMENSION(SIZE(PRGT))   :: ZSIGMOIDE
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!
!*       1.     Deposition of vapour on graupel
!               -------------------------------
!
IF (LHOOK) CALL DR_HOOK('LIMA_GRAUPEL_DEPOSITION', 0, ZHOOK_HANDLE)
P_TH_DEPG(:) = 0.0
P_RG_DEPG(:) = 0.0
WHERE ( PRGT(:)>LIMAP%XRTMIN(6) .AND. PCGT(:)>LIMAP%XCTMIN(6) .AND. ODCOMPUTE(:) )
   P_RG_DEPG(:) = PSSI(:) / PAI(:) * PCGT(:) *                      &
                ( LIMAM%X0DEPG*PLBDG(:)**LIMAM%XEX0DEPG + LIMAM%X1DEPG*PCJ(:)*PLBDG(:)**LIMAM%XEX1DEPG )
   P_TH_DEPG(:) = P_RG_DEPG(:)*PLSFACT(:)
END WHERE
!
IF (LIMAP%LSIGMOIDE_G) THEN
     ZSIGMOIDE(:)      =  1/(1 +  exp(-LIMAP%XSIGMOIDE_G*(PRGT(:)-LIMAM%XMINDG/PRHODREF(:))))
     P_TH_DEPG(:) = P_TH_DEPG(:) * ZSIGMOIDE(:)
     P_RG_DEPG(:) = P_RG_DEPG(:) * ZSIGMOIDE(:)
END IF
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('LIMA_GRAUPEL_DEPOSITION', 1, ZHOOK_HANDLE)
END SUBROUTINE LIMA_GRAUPEL_DEPOSITION
END MODULE MODE_LIMA_GRAUPEL_DEPOSITION
