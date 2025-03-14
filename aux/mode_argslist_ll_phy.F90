!MNH_LIC Copyright 2023-2023 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!
!!    #######################
MODULE MODE_ARGSLIST_ll_PHY
!
 USE MODE_ll
 USE MODD_ARGSLIST_ll, ONLY : LIST_ll
 USE MODD_DIMPHYEX,       ONLY: DIMPHYEX_t
!
 CONTAINS
!
 SUBROUTINE ADD3DFIELD_ll_PHY(D, TPLIST, PFIELD, HNAME)
!!    ###############################################
!
!!****  *ADD3DFIELD_ll_PHY* -
!
!!    Purpose
!!    -------
!     This routine is used as an interface to ADD3DFIELD_ll for
!     unpacking horizontal dimensions
!
!!    Reference
!!    ---------
!
!     see PHYEX documentation
!
!!    Implicit Arguments
!!    ------------------
!
!     Module MODD_ARGSLIST :
!         LIST_ll : list of fields
!         DIMPHYEX_t: PHYEX dimensions
!
!!    Author
!!    ------
!
!     Q.Rodier
!
!!    Modifications
!!    -------------
!     Original    August, 3, 2023
!
!-------------------------------------------------------------------------------
!
  IMPLICIT NONE
!
  TYPE(DIMPHYEX_t),        INTENT(IN)   :: D
  TYPE(LIST_ll), POINTER         :: TPLIST   ! list of fields
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), TARGET :: PFIELD   ! field which is unpaked here
!                                              of fields
  CHARACTER(LEN=*), INTENT(IN) :: HNAME ! Name of the field to be added

  CALL ADD3DFIELD_ll(TPLIST, PFIELD, HNAME)
              
 END SUBROUTINE ADD3DFIELD_ll_PHY
END MODULE MODE_ARGSLIST_ll_PHY
