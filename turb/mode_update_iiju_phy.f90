!MNH_LIC Copyright 1994-2022 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
MODULE MODE_UPDATE_IIJU_PHY
IMPLICIT NONE
CONTAINS
SUBROUTINE UPDATE_IIJU_PHY(D,PVAR)
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
!     ##############################################################
!
!!****  *MODE_UPDATE_IIJU_PHY* -
!!
!!    PURPOSE
!!    -------
!!    This routine update IIU-1 and IJU-1 values to (IIU,IJU) values in the PHYEX
!!    package where all arrays have a single dimension for the horizontal coordinates
!!    i.e. ni*nj
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!
!!      Q. Rodier  - Meteo-France -
!!
!!    MODIFICATIONS
!!    -------------
!!     Original     18/08/22
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
!
IMPLICIT NONE
!
!*       0.1   Declaration of arguments
TYPE(DIMPHYEX_t),         INTENT(IN)   :: D
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(INOUT)    :: PVAR  ! working variable
!
INTEGER :: IIE,IIB,IJE,IJB,IIU,IJU,IKT
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('UPDATE_IIJU_PHY',0,ZHOOK_HANDLE)
IIE=D%NIEC
IIB=D%NIBC
IJE=D%NJEC
IJB=D%NJBC
IIU=D%NIT
IJU=D%NJT
IKT=D%NKT
!
PVAR(IIU,IJB:IJE,:) = PVAR(IIU-1,IJB:IJE,:)
PVAR(IIB:IIE,IJU,:) = PVAR(IIB:IIE,IJU-1,:)
!
IF (LHOOK) CALL DR_HOOK('UPDATE_IIJU_PHY',1,ZHOOK_HANDLE)
END SUBROUTINE UPDATE_IIJU_PHY
END MODULE MODE_UPDATE_IIJU_PHY
