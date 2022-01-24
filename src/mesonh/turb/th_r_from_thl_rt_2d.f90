!MNH_LIC Copyright 2006-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ######spl
      MODULE MODI_TH_R_FROM_THL_RT_2D
      INTERFACE
      SUBROUTINE TH_R_FROM_THL_RT_2D(HFRAC_ICE,PFRAC_ICE,PP,             &
                                  PTHL, PRT, PTH, PRV, PRL, PRI,         &
                                  PRSATW, PRSATI, PRR, PRS, PRG, PRH     )
CHARACTER(len=1),     INTENT(IN) :: HFRAC_ICE
REAL, DIMENSION(:,:), INTENT(INOUT) :: PFRAC_ICE
REAL, DIMENSION(:,:), INTENT(IN) :: PP          ! Pressure
REAL, DIMENSION(:,:), INTENT(IN) :: PTHL    ! thetal to transform into th
REAL, DIMENSION(:,:),INTENT(IN)  :: PRT    ! Total mixing ratios to transform into rv,rc and ri
REAL, DIMENSION(:),OPTIONAL,INTENT(IN) :: PRR, PRS, PRG, PRH
REAL, DIMENSION(:,:), INTENT(OUT):: PTH    ! th
REAL, DIMENSION(:,:), INTENT(OUT):: PRV    ! vapor mixing ratio
REAL, DIMENSION(:,:), INTENT(INOUT):: PRL    ! vapor mixing ratio
REAL, DIMENSION(:,:), INTENT(INOUT):: PRI    ! vapor mixing ratio
REAL, DIMENSION(:,:), INTENT(OUT)  :: PRSATW ! estimated mixing ration at saturation over water
REAL, DIMENSION(:,:), INTENT(OUT)  :: PRSATI ! estimated mixing ration at saturation over ice

      END SUBROUTINE TH_R_FROM_THL_RT_2D

      END INTERFACE

      END MODULE MODI_TH_R_FROM_THL_RT_2D

!     ######spl
      SUBROUTINE TH_R_FROM_THL_RT_2D(HFRAC_ICE,PFRAC_ICE,PP,             &
                                  PTHL, PRT, PTH, PRV, PRL, PRI,         &
                                  PRSATW, PRSATI, PRR, PRS, PRG, PRH     )
!     #################################################################
!
!
!!****  *TH_R_FROM_THL_RT_2D* - computes the non-conservative variables
!!                          from conservative variables
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!    
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
!!      Julien PERGAUD      * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original         13/03/06
!!      Sébastien Riette April 2011: code moved in th_r_from_thl_rt_3D
!!
!! --------------------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
!
!
USE MODI_TH_R_FROM_THL_RT_3D

IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
CHARACTER(len=1),     INTENT(IN) :: HFRAC_ICE
REAL, DIMENSION(:,:), INTENT(INOUT) :: PFRAC_ICE
REAL, DIMENSION(:,:), INTENT(IN) :: PP     ! Pressure
REAL, DIMENSION(:,:), INTENT(IN) :: PTHL   ! Liquid pot. temp.
REAL, DIMENSION(:,:), INTENT(IN) :: PRT    ! Total mixing ratios
REAL, DIMENSION(:,:),OPTIONAL,INTENT(IN) :: PRR, PRS, PRG, PRH
REAL, DIMENSION(:,:), INTENT(OUT):: PTH    ! Potential temp.
REAL, DIMENSION(:,:), INTENT(OUT):: PRV    ! vapor mixing ratio
REAL, DIMENSION(:,:), INTENT(INOUT):: PRL  ! cloud mixing ratio
REAL, DIMENSION(:,:), INTENT(INOUT):: PRI  ! ice   mixing ratio
REAL, DIMENSION(:,:), INTENT(OUT)  :: PRSATW ! estimated mixing ration at saturation over water
REAL, DIMENSION(:,:), INTENT(OUT)  :: PRSATI ! estimated mixing ration at saturation over ice

!
!-------------------------------------------------------------------------------
!
!       0.2  declaration of local variables
!
!----------------------------------------------------------------------------
!
REAL, DIMENSION(SIZE(PP,1),SIZE(PP,2)) :: ZRR, ZRS, ZRG, ZRH
INTEGER :: JK
!----------------------------------------------------------------------------
!
!*      1 Initialisation
!         --------------
!
ZRR(:,:)=0.
ZRS(:,:)=0.
ZRG(:,:)=0.
ZRH(:,:)=0.
IF(PRESENT(PRR)) ZRR(:,:)=PRR(:,:)
IF(PRESENT(PRS)) ZRS(:,:)=PRS(:,:)
IF(PRESENT(PRG)) ZRG(:,:)=PRG(:,:)
IF(PRESENT(PRH)) ZRH(:,:)=PRH(:,:)
!
!
!       2 Call of 1d version
!         ------------------
DO JK=1, SIZE(PTHL,2)
  CALL TH_R_FROM_THL_RT_1D(HFRAC_ICE,PFRAC_ICE(:,JK),PP(:,JK),             &
                                PTHL(:,JK), PRT(:,JK), PTH(:,JK),       &
                                PRV(:,JK), PRL(:,JK), PRI(:,JK),        &
                                PRSATW(:,JK), PRSATI(:,JK),                &
                                ZRR(:,JK), ZRS(:,JK), ZRG(:,JK), ZRH(:,JK))
ENDDO


END SUBROUTINE TH_R_FROM_THL_RT_2D
