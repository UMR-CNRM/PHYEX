!      ####################
       MODULE MODI_INI_LIMA
!      ####################
!
INTERFACE
      SUBROUTINE INI_LIMA (KULOUT, PTSTEP, PDZMIN, KSPLITR, KSPLITG)
!
INTEGER,                 INTENT(IN) :: KULOUT    ! output logical unit number
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
      SUBROUTINE INI_LIMA (KULOUT, PTSTEP, PDZMIN, KSPLITR, KSPLITG)
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
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
USE MODD_REF
USE MODD_PARAM_LIMA
USE MODD_PARAMETERS
USE MODD_LUNIT
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
INTEGER,                 INTENT(IN) :: KULOUT    ! output logical unit number
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
REAL     :: ZVTRMAX
!
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
CALL FMLOOK_ll(CLUOUT0,CLUOUT0,ILUOUT0,IRESP)
!
!
! KSPLITR
ZVTRMAX = 30.         ! Maximum rain drop fall speed
!
KSPLITR = 1
SPLITR : DO
   ZT = PTSTEP / FLOAT(KSPLITR)
   IF ( ZT * ZVTRMAX / PDZMIN < 1.0) EXIT SPLITR
   KSPLITR = KSPLITR + 1
END DO SPLITR
!
!
! KSPLITG
ZVTRMAX = 30.                          
IF( LHAIL_LIMA ) THEN
   ZVTRMAX = 60. ! Hail case 
END IF
!
KSPLITG = 1
SPLITG : DO
   ZT = 2.* PTSTEP / FLOAT(KSPLITG)
   IF ( ZT * ZVTRMAX / PDZMIN .LT. 1.) EXIT SPLITG
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
XRTMIN(1) = 1.0E-20   ! rv
XRTMIN(2) = 1.0E-20   ! rc
!XRTMIN(3) = 1.0E-20   ! rr
XRTMIN(3) = 1.0E-17   ! rr
XRTMIN(4) = 1.0E-20   ! ri
XRTMIN(5) = 1.0E-15   ! rs
XRTMIN(6) = 1.0E-15   ! rg
XRTMIN(7) = 1.0E-15   ! rh
ALLOCATE( XCTMIN(7) )
XCTMIN(1) = 1.0       ! Not used
XCTMIN(2) = 1.0E+4    ! Nc
!XCTMIN(3) = 1.0E+1    ! Nr
XCTMIN(3) = 1.0E-3    ! Nr
XCTMIN(4) = 1.0E-3    ! Ni
XCTMIN(5) = 1.0E-3    ! Not used
XCTMIN(6) = 1.0E-3    ! Not used
XCTMIN(7) = 1.0E-3    ! Not used
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
CALL INI_LIMA_WARM(KULOUT, PTSTEP, PDZMIN)
!
CALL INI_LIMA_COLD_MIXED(KULOUT, PTSTEP, PDZMIN)
!
!------------------------------------------------------------------------------
!
END SUBROUTINE INI_LIMA
