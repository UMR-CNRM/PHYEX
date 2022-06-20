!MNH_LIC Copyright 2019-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!
!     #######################
MODULE MODD_IBM_PARAM_n
  !   #######################
  IMPLICIT NONE
  LOGICAL          :: LIBM,LIBM_TROUBLE
  REAL, DIMENSION(:,:,:,:)       , POINTER :: XIBM_LS=>NULL()      ! LSF for MNH
  REAL, DIMENSION(:,:,:)          , POINTER :: XIBM_XMUT=>NULL()
END MODULE MODD_IBM_PARAM_n
!
