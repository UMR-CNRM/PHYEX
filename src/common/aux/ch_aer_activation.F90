!ORILAM_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!ORILAM_LIC This is part of the ORILAM software governed by the CeCILL-C licence
!ORILAM_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!ORILAM_LIC for details.
!-----------------------------------------------------------------
!!   ##############################
     MODULE MODI_CH_AER_ACTIVATION
!!   ##############################
!!
INTERFACE
!
SUBROUTINE CH_AER_ACTIVATION(PSVT,PTEMP,  PWT, PTRAD, &
                             PRHODREF,PPABST, PNCN, PMCN,&
                             PSOLORG, PMI, PSMAX)

IMPLICIT NONE

REAL,  DIMENSION(:,:),  INTENT(INOUT) :: PSVT  ! Aerosol concentration
REAL,  DIMENSION(:),    INTENT(IN)    :: PTEMP ! Air temperature (K)
REAL,  DIMENSION(:),    INTENT(IN)    :: PRHODREF ! Air density (kg/m3)
REAL,  DIMENSION(:),    INTENT(IN)    :: PWT   ! Activation vertical velocity (m/s)
REAL,  DIMENSION(:),    INTENT(IN)    :: PTRAD ! Activation cooling radiative tendency (K/s)
REAL,  DIMENSION(:),    INTENT(IN)    :: PPABST! Air pressure (Pa)
REAL,  DIMENSION(:,:),  INTENT(IN)    :: PSOLORG ! ratio of SOA in acqueous phase
REAL,  DIMENSION(:,:),  INTENT(IN)    :: PMI   ! Molecular weight (g/mol)
REAL,  DIMENSION(:),    INTENT(OUT)   :: PNCN  ! Number of activated aerosol (#/m3) 
REAL,  DIMENSION(:),    INTENT(OUT)   :: PMCN  ! Mass of activated aerosol (ug/m3)
REAL,  DIMENSION(:),    INTENT(INOUT) :: PSMAX ! Maximum supersaturation


END SUBROUTINE CH_AER_ACTIVATION
!!
END INTERFACE
!!
END MODULE MODI_CH_AER_ACTIVATION
!!
!!   #######################################
SUBROUTINE CH_AER_ACTIVATION(PSVT,PTEMP,  PWT, PTRAD, &
                             PRHODREF,PPABST, PNCN, PMCN,&
                             PSOLORG, PMI, PSMAX)
!!   #######################################
!!
!!   PURPOSE
!!   -------
!!   Input of Abdul-Razzak activation scheme
!!   Here we compute the size distribution of aerosols together with the 
!!   dissociative ions and the soluble fraction need in the parameterization.
!!   Dynamical variables are also passed by argument.
!!   REFERENCE
!!   ---------
!!   none
!!
!!   AUTHOR
!!    ------
!!    Pierre TULET (GMEI)  
!!
!!   MODIFICATIONS
!!    -------------
!!   Original
!!
!!   Empty routines for PHYEX-offline (Q. Rodier 12/2024)
!!
!!   IMPLICIT ARGUMENTS
!
!!
IMPLICIT NONE

!! Arguments variables

REAL,  DIMENSION(:,:),  INTENT(INOUT) :: PSVT  ! Aerosol concentration
REAL,  DIMENSION(:),    INTENT(IN)    :: PTEMP ! Air temperature (K)
REAL,  DIMENSION(:),    INTENT(IN)    :: PRHODREF ! Air density (kg/m3)
REAL,  DIMENSION(:),    INTENT(IN)    :: PWT   ! Activation vertical velocity (m/s)
REAL,  DIMENSION(:),    INTENT(IN)    :: PTRAD ! Activation cooling radiative tendency (K/s)
REAL,  DIMENSION(:),    INTENT(IN)    :: PPABST! Air pressure (Pa)
REAL,  DIMENSION(:,:),  INTENT(IN)    :: PSOLORG ! ratio of SOA in acqueous phase
REAL,  DIMENSION(:,:),  INTENT(IN)    :: PMI   ! Molecular weight (g/mol)
REAL,  DIMENSION(:),    INTENT(OUT)   :: PNCN  ! Number of activated aerosol (#/m3) 
REAL,  DIMENSION(:),    INTENT(OUT)   :: PMCN  ! Mass of activated aerosol (ug/m3)
REAL,  DIMENSION(:),    INTENT(INOUT) :: PSMAX ! Maximum supersaturation
!
END SUBROUTINE CH_AER_ACTIVATION
