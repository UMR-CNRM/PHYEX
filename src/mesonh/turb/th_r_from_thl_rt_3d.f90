!MNH_LIC Copyright 2011-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ######spl
      MODULE MODI_TH_R_FROM_THL_RT_3D
!     ###############################
INTERFACE
!
      SUBROUTINE TH_R_FROM_THL_RT_3D(HFRAC_ICE,PFRAC_ICE,PP,             &
                                  PTHL, PRT, PTH, PRV, PRL, PRI, &
                                  PRSATW, PRSATI, PRR, PRS, PRG, PRH       )

CHARACTER(len=1),       INTENT(IN) :: HFRAC_ICE
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PFRAC_ICE
REAL, DIMENSION(:,:,:), INTENT(IN) :: PP          ! Pressure
REAL, DIMENSION(:,:,:), INTENT(IN) :: PTHL    ! thetal to transform into th
REAL, DIMENSION(:,:,:),INTENT(IN)  :: PRT    ! Total mixing ratios to transform into rv,rc and ri
REAL, DIMENSION(:,:,:),OPTIONAL,INTENT(IN) :: PRR, PRS, PRG, PRH
REAL, DIMENSION(:,:,:), INTENT(OUT):: PTH    ! th
REAL, DIMENSION(:,:,:), INTENT(OUT):: PRV    ! vapor mixing ratio
REAL, DIMENSION(:,:,:), INTENT(INOUT):: PRL    ! vapor mixing ratio
REAL, DIMENSION(:,:,:), INTENT(INOUT):: PRI    ! vapor mixing ratio
REAL, DIMENSION(:,:,:), INTENT(OUT)  :: PRSATW ! estimated mixing ration at saturation over water
REAL, DIMENSION(:,:,:), INTENT(OUT)  :: PRSATI ! estimated mixing ration at saturation over ice
!
END SUBROUTINE TH_R_FROM_THL_RT_3D
!
END INTERFACE
!
END MODULE MODI_TH_R_FROM_THL_RT_3D
!     ######spl
      SUBROUTINE TH_R_FROM_THL_RT_3D(HFRAC_ICE,PFRAC_ICE,PP,             &
                                  PTHL, PRT, PTH, PRV, PRL, PRI, &
                                  PRSATW, PRSATI, PRR, PRS, PRG, PRH       )
!     #################################################################
!
!
!!****  *TH_R_FROM_THL_RT_3D* - computes the non-conservative variables
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
!!      Original         15/03/2011
!!      S. Riette April 2011 : code moved in th_r_from_thl_rt_1d
!!
!! --------------------------------------------------------------------------
!
!*      0. DECLARATIONS
!          ------------
!
USE MODI_TH_R_FROM_THL_RT_1D
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
CHARACTER(len=1),       INTENT(IN) :: HFRAC_ICE
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PFRAC_ICE
REAL, DIMENSION(:,:,:), INTENT(IN) :: PP          ! Pressure
REAL, DIMENSION(:,:,:), INTENT(IN) :: PTHL    ! thetal to transform into th
REAL, DIMENSION(:,:,:),INTENT(IN)  :: PRT    ! Total mixing ratios to transform into rv,rc and ri
REAL, DIMENSION(:,:,:),OPTIONAL,INTENT(IN) :: PRR, PRS, PRG, PRH
REAL, DIMENSION(:,:,:), INTENT(OUT):: PTH    ! th
REAL, DIMENSION(:,:,:), INTENT(OUT):: PRV    ! vapor mixing ratio
REAL, DIMENSION(:,:,:), INTENT(INOUT):: PRL    ! vapor mixing ratio
REAL, DIMENSION(:,:,:), INTENT(INOUT):: PRI    ! vapor mixing ratio
REAL, DIMENSION(:,:,:), INTENT(OUT)  :: PRSATW ! estimated mixing ration at saturation over water
REAL, DIMENSION(:,:,:), INTENT(OUT)  :: PRSATI ! estimated mixing ration at saturation over ice
!
!-------------------------------------------------------------------------------
!
!       0.2  declaration of local variables
REAL, DIMENSION(SIZE(PTHL,1),SIZE(PTHL,2),SIZE(PTHL,3)) :: ZRR, ZRS, ZRG, ZRH
INTEGER :: JJ, JK
!----------------------------------------------------------------------------
!
!*      1 Initialisation
!         --------------
!
ZRR(:,:,:)=0.
ZRS(:,:,:)=0.
ZRG(:,:,:)=0.
ZRH(:,:,:)=0.
IF(PRESENT(PRR)) ZRR(:,:,:)=PRR(:,:,:)
IF(PRESENT(PRS)) ZRS(:,:,:)=PRS(:,:,:)
IF(PRESENT(PRG)) ZRG(:,:,:)=PRG(:,:,:)
IF(PRESENT(PRH)) ZRH(:,:,:)=PRH(:,:,:)
!
!
!       2 Call of 1d version
!         ------------------
DO JK=1, SIZE(PTHL,3)
  DO JJ=1, SIZE(PTHL,2)
    CALL TH_R_FROM_THL_RT_1D(HFRAC_ICE,PFRAC_ICE(:,JJ,JK),PP(:,JJ,JK),             &
                                  PTHL(:,JJ,JK), PRT(:,JJ,JK), PTH(:,JJ,JK),       &
                                  PRV(:,JJ,JK), PRL(:,JJ,JK), PRI(:,JJ,JK),        &
                                  PRSATW(:,JJ,JK), PRSATI(:,JJ,JK),                &
                                  ZRR(:,JJ,JK), ZRS(:,JJ,JK), ZRG(:,JJ,JK), ZRH(:,JJ,JK))
  ENDDO
ENDDO

END SUBROUTINE TH_R_FROM_THL_RT_3D
