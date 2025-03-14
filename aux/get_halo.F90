!MNH_LIC Copyright 1994-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
! Modifications:
!  P. Wautelet 20/05/2019: add name argument to ADDnFIELD_ll + new ADD4DFIELD_ll subroutine
!-----------------------------------------------------------------
!     ####################
      MODULE MODI_GET_HALO
!     ####################
!
IMPLICIT NONE
INTERFACE
!
SUBROUTINE GET_HALO_PHY(D,PSRC)
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
IMPLICIT NONE
!
TYPE(DIMPHYEX_t),        INTENT(IN)   :: D
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)  :: PSRC    ! variable at t
!
END SUBROUTINE GET_HALO_PHY
!
SUBROUTINE GET_HALO(PSRC)
IMPLICIT NONE
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PSRC    ! variable at t
!
END SUBROUTINE GET_HALO
!
END INTERFACE
!
END MODULE MODI_GET_HALO         
!
!-------------------------------------------------------------------------------
!     #########################
      SUBROUTINE GET_HALO(PSRC)
!     #########################
!
IMPLICIT NONE
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PSRC    ! variable at t
!
END SUBROUTINE GET_HALO
!-----------------------------------------------------------------------
!     #########################
      SUBROUTINE GET_HALO_PHY(D,PSRC)
!     #########################
!
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
IMPLICIT NONE
!
TYPE(DIMPHYEX_t),        INTENT(IN)   :: D
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)  :: PSRC    ! variable at t
!
END SUBROUTINE GET_HALO_PHY
