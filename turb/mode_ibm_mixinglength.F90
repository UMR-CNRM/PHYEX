!MNH_LIC Copyright 2019-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!
MODULE MODE_IBM_MIXINGLENGTH
IMPLICIT NONE
CONTAINS
SUBROUTINE IBM_MIXINGLENGTH(D,PLM,PLEPS,PMU,PHI,PTKE)
  !     ###################################################
  !
  !****  *IBM_MIXINGLENGTH* - Alteration of the mixing lenght (IBM)
  !
  !    PURPOSE
  !    -------
  !     The limitation is corrected for the immersed bonudary method:
  !        => using the level set phi 
  !        => LM < k(-phi)
  !
  !    METHOD
  !    ------
  !
  !    INDEX
  !    -----
  !
  !    IMPLICIT ARGUMENTS
  !    ------------------
  !
  !    REFERENCE
  !    ---------
  !
  !    AUTHOR
  !    ------
  !	
  !      Franck Auguste       * CERFACS(AE) *
  !
  !    MODIFICATIONS
  !    -------------
  !      Original    01/01/2019
  !
  !-------------------------------------------------------------------------------
  !
  !**** 0.    DECLARATIONS
  !     ------------------
  !
  ! module
  USE MODD_DIMPHYEX,   ONLY: DIMPHYEX_t
  !
  IMPLICIT NONE
  !
  !------------------------------------------------------------------------------
  !
  !       0.1   Declaration of arguments
  TYPE(DIMPHYEX_t),       INTENT(IN)    :: D
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(INOUT) :: PLM
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(INOUT) :: PLEPS
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT)   :: PMU
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)    :: PHI
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)    :: PTKE
  !
  ! IBM CAN NOT BE USED WITH AROME
  !
END SUBROUTINE IBM_MIXINGLENGTH
END MODULE MODE_IBM_MIXINGLENGTH
