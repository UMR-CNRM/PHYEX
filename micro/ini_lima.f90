!MNH_LIC Copyright 2013-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!      ####################
       MODULE MODI_INI_LIMA
!      ####################
!
INTERFACE
      SUBROUTINE INI_LIMA (PTSTEP, PDZMIN, KSPLITR, KSPLITG)
!
INTEGER,                 INTENT(OUT):: KSPLITR   ! Number of small time step
                                                 ! integration for  rain
                                                 ! sedimendation
INTEGER,                 INTENT(OUT):: KSPLITG   ! Number of small time step
                                                 ! integration for graupel
                                                 ! sedimendation
REAL,                    INTENT(IN) :: PTSTEP    ! Effective Time step 
REAL,                    INTENT(IN) :: PDZMIN    ! minimun vertical mesh size
!
END SUBROUTINE INI_LIMA
!
END INTERFACE
!
END MODULE MODI_INI_LIMA
!     ######################################################
      SUBROUTINE INI_LIMA (PTSTEP, PDZMIN, KSPLITR, KSPLITG)
!     ######################################################
!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to initialize the constants used in the 
!!    microphysical scheme LIMA. 
!!
!!    AUTHOR
!!    ------
!!      J.-M. Cohard     * Laboratoire d'Aerologie*
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!      S.    Berthet    * Laboratoire d'Aerologie*
!!      B.    Vié        * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original             ??/??/13 
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!  P. Wautelet 26/04/2019: replace non-standard FLOAT function by REAL function
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
USE MODD_REF
USE MODD_PARAM_LIMA
USE MODD_PARAMETERS
!USE MODD_LUNIT, ONLY : TLUOUT0
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
INTEGER,                 INTENT(OUT):: KSPLITR   ! Number of small time step
                                                 ! integration for  rain
                                                 ! sedimendation
INTEGER,                 INTENT(OUT):: KSPLITG   ! Number of small time step
                                                 ! integration for graupel or hail
                                                 ! sedimendation
REAL,                    INTENT(IN) :: PTSTEP    ! Effective Time step 
REAL,                    INTENT(IN) :: PDZMIN    ! minimun vertical mesh size
!
!*       0.2   Declarations of local variables :
!
REAL     :: ZT      ! Work variable
REAL, DIMENSION(7)  :: ZVTRMAX
!
INTEGER  :: JI
INTEGER  :: ILUOUT0 ! Logical unit number for output-listing
INTEGER  :: IRESP   ! Return code of FM-routines
!  
!-------------------------------------------------------------------------------
!
!
!*       1.     INIT OUTPUT LISTING, COMPUTE KSPLITR AND KSPLITG
!   	        ------------------------------------------------
!
!
! Init output listing
!ILUOUT0 = TLUOUT0%NLU
!
!
ZVTRMAX(2) = 0.3         ! Maximum cloud droplet fall speed
ZVTRMAX(3) = 15.         ! Maximum rain drop fall speed
ZVTRMAX(4) = 1.5         ! Maximum ice crystal fall speed
ZVTRMAX(5) = 3.0         ! Maximum snow fall speed
ZVTRMAX(6) = 15.         ! Maximum graupel fall speed
ZVTRMAX(7) = 30.         ! Maximum hail fall speed
!
! NSPLITSED
!
DO JI=2,7
   NSPLITSED(JI) = 1
   SPLIT : DO
      ZT = PTSTEP / REAL(NSPLITSED(JI))
      IF ( ZT * ZVTRMAX(JI) / PDZMIN < 1.0) EXIT SPLIT
      NSPLITSED(JI) = NSPLITSED(JI) + 1
   END DO SPLIT
END DO
!
! KSPLITR
!
KSPLITR = 1
SPLITR : DO
   ZT = PTSTEP / REAL(KSPLITR)
   IF ( ZT * ZVTRMAX(7) / PDZMIN < 1.0) EXIT SPLITR
   KSPLITR = KSPLITR + 1
END DO SPLITR
!
!
! KSPLITG
!
KSPLITG = 1
SPLITG : DO
   ZT = 2.* PTSTEP / REAL(KSPLITG)
   IF ( ZT * ZVTRMAX(7) / PDZMIN .LT. 1.) EXIT SPLITG
   KSPLITG = KSPLITG + 1
END DO SPLITG
!
!
!
IF (ALLOCATED(XRTMIN)) RETURN    ! In case of nesting microphysics, constants of
                                 ! MODD_RAIN_C2R2_PARAM are computed only once.
!
!
! Set bounds for mixing ratios and concentrations
ALLOCATE( XRTMIN(7) )
XRTMIN(1) = 1.0E-10   ! rv
XRTMIN(2) = 1.0E-10   ! rc
XRTMIN(3) = 1.0E-10   ! rr
XRTMIN(4) = 1.0E-10   ! ri
XRTMIN(5) = 1.0E-10   ! rs
XRTMIN(6) = 1.0E-10   ! rg
XRTMIN(7) = 1.0E-10   ! rh
ALLOCATE( XCTMIN(7) )
XCTMIN(1) = 1.0       ! Not used
XCTMIN(2) = 1.0E-3    ! Nc
XCTMIN(3) = 1.0E-3    ! Nr
XCTMIN(4) = 1.0E-3    ! Ni
XCTMIN(5) = 1.0E-3    ! Ns
XCTMIN(6) = 1.0E-3    ! Ng
XCTMIN(7) = 1.0E-3    ! Nh
!
!
! Air density fall speed correction
XCEXVT = 0.4
!
!------------------------------------------------------------------------------
!
!
!
!*       2.     DEFINE SPECIES CHARACTERISTICS AND PROCESSES CONSTANTS
!   	        ------------------------------------------------------
!
!
CALL INI_LIMA_WARM(PTSTEP, PDZMIN)
!
CALL INI_LIMA_COLD_MIXED(PTSTEP, PDZMIN)
!
!------------------------------------------------------------------------------
!
END SUBROUTINE INI_LIMA
