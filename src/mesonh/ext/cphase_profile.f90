!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!#########################
MODULE MODI_CPHASE_PROFILE
!#########################
!
INTERFACE
!
      SUBROUTINE CPHASE_PROFILE (PZHAT,PCPHASE,PCPHASE_PBL,PCPHASE_PROFILE,PTKEM)
! 
REAL, DIMENSION(:)              ,  INTENT(IN)  :: PZHAT       ! height level without orography
REAL                            ,  INTENT(IN)  :: PCPHASE     ! prescribed phase velocity
REAL                            ,  INTENT(IN)  :: PCPHASE_PBL ! prescribed phase velocity
REAL, DIMENSION(:,:)            ,  INTENT(OUT) :: PCPHASE_PROFILE ! profile of Cphase speed
REAL, DIMENSION(:,:),OPTIONAL   ,  INTENT(IN)  :: PTKEM ! TKE at t-dt
!
END SUBROUTINE CPHASE_PROFILE
!
END INTERFACE
!
END MODULE MODI_CPHASE_PROFILE
!
!     ##########################################################################
      SUBROUTINE CPHASE_PROFILE (PZHAT,PCPHASE,PCPHASE_PBL,PCPHASE_PROFILE,PTKEM)
!     ##########################################################################
!
!!****  *CPHASE_PROFILE* - defines a non-constant vertical profile for Cphase
!!                         velocity
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
!!   
!!        
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!	V. Masson & C. Lac     * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    08/2010
!!      Escobar     9/11/2010 : array bound problem if NO Turb =>  PTKEM optional 
!!      C.Lac       06/2013   : correction and introduction of PCPHASE_PBL
!!      
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CTURB
USE MODD_PARAMETERS
!
IMPLICIT NONE
!
!
!*       0.1   declarations of arguments
! 
! 
! 
REAL, DIMENSION(:)              ,  INTENT(IN)  :: PZHAT       ! height level without orography
REAL                            ,  INTENT(IN)  :: PCPHASE     ! prescribed phase velocity
REAL                            ,  INTENT(IN)  :: PCPHASE_PBL ! prescribed phase velocity
REAL, DIMENSION(:,:)            ,  INTENT(OUT) :: PCPHASE_PROFILE ! profile of Cphase speed
REAL, DIMENSION(:,:),OPTIONAL   ,  INTENT(IN)  :: PTKEM ! TKE at t-dt
!
!*       0.2   declarations of local variables
!
INTEGER             :: IKB       ! indice K Beginning in z direction 
INTEGER             :: IKE       ! indice K End       in z direction 
! 
REAL, DIMENSION(SIZE(PCPHASE_PROFILE,1)) :: ZTKE, ZTKEMIN           
INTEGER                                  :: JL,JK,JKTKE
!
!-------------------------------------------------------------------------------
!
!*       1.    PROLOGUE
!              --------
!
!*       1.1  Compute dimensions of arrays and other indices
! 
IKB = 1 + JPVEXT
IKE = SIZE(PCPHASE_PROFILE,2) - JPVEXT
!
!
!*       1.2  Initializations
!
!
PCPHASE_PROFILE   = 0.0
ZTKEMIN           = PZHAT(IKE)
ZTKE              = PZHAT(IKE-1)
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
    IF (PRESENT(PTKEM)) THEN
!
     DO JL = 1,SIZE(PCPHASE_PROFILE,1)
      JKTKE=IKE-1
      DO JK = IKB, IKE-1 
       IF (PTKEM(JL,JK) < 5.*XTKEMIN ) THEN
           ZTKE (JL) = PZHAT (JK)
           JKTKE = JK
           EXIT
        END IF
      END DO
      DO JK = JKTKE+1,IKE
       IF (PTKEM(JL,JK) == XTKEMIN ) THEN
           ZTKEMIN (JL) = PZHAT (JK)
           EXIT
        END IF
      END DO
     END DO
!
    ELSE 
      ZTKE (:) = 1000.
      ZTKEMIN (:) = 2000.
    END IF
!
    DO JL = 1,SIZE(PCPHASE_PROFILE,1)
     DO JK = IKB, IKE
      IF (PZHAT(JK) > ZTKEMIN (JL) ) THEN
        PCPHASE_PROFILE(JL,JK) = PCPHASE
      ELSE IF  (PZHAT(JK) < ZTKE (JL) ) THEN
        PCPHASE_PROFILE(JL,JK) = PCPHASE_PBL
      ELSE 
        PCPHASE_PROFILE(JL,JK) = 1./(ZTKEMIN (JL) - ZTKE (JL)) * &
         ((PZHAT(JK) - ZTKE(JL)) * PCPHASE + (ZTKEMIN (JL) - PZHAT(JK)) * PCPHASE_PBL )
      END IF
     END DO
    END DO
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE CPHASE_PROFILE
