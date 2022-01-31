!MNH_LIC Copyright 2019-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!
!     ############################
MODULE MODI_IBM_MIXINGLENGTH
  !     ############################
  !
  INTERFACE 
     !
     SUBROUTINE IBM_MIXINGLENGTH(PLM,PLEPS,PMU,PHI,PTKE)
       !
       REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PLM
       REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PLEPS
       REAL, DIMENSION(:,:,:), INTENT(OUT)   :: PMU
       REAL, DIMENSION(:,:,:), INTENT(IN)    :: PHI
       REAL, DIMENSION(:,:,:), INTENT(IN)    :: PTKE
       !
     END SUBROUTINE IBM_MIXINGLENGTH
     !
  END INTERFACE
END MODULE MODI_IBM_MIXINGLENGTH
  !
  SUBROUTINE IBM_MIXINGLENGTH(PLM,PLEPS,PMU,PHI,PTKE)
       !
       REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PLM
       REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PLEPS
       REAL, DIMENSION(:,:,:), INTENT(OUT)   :: PMU
       REAL, DIMENSION(:,:,:), INTENT(IN)    :: PHI
       REAL, DIMENSION(:,:,:), INTENT(IN)    :: PTKE
       !
     END SUBROUTINE IBM_MIXINGLENGTH
     !
