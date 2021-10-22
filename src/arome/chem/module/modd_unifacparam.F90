!     ######spl
MODULE MODD_UNIFACPARAM

  USE MODD_GLO
!
  implicit none
!
!*************************************************************************
!
!Include file : unifacparma.h
!
!Purpose: Unifac parameters (replace file read)
!
!Include dependencies:  Included in Unidriver.c
!
!Notes: 10 Type AQ SOA + H2O at 25 C
!       (Data of Lyman, Reehl, Rosenblatt, 1990)
!
!       NMOL and NFUNC need to match DIMMOL and DIMFUN in unidriver.c
!
!       Parameters inputted as in the original fortran input file, therefore
!       a transpose is needed in the unidriver C program for matrices A and NU
!       in order to pass them properly into the Fortran unifac routine.
!
!revision History:  1. Developed by Betty Pun, AER, December, 1999
!                      under CARB funding
!
!**************************************************************************
!
!***************************************************************************
!****************NUMBERS NEEDED FOR GRIFFINS AQUOUS PHASE******************
!**************************************************************************
!
integer, parameter ::  NMOL_AQ = NBSP + 1  ! Number of aquous molecules is number of species+water
integer, parameter ::  NFUNC_AQ = 17       ! Number of functional groups
REAL, dimension(NFUNC_AQ) :: RG_AQ         !Group volume parameters
REAL, DIMENSION(NFUNC_AQ) :: QG_AQ         ! Group surface parameters
INTEGER, DIMENSION(NMOL_AQ, NFUNC_AQ)  ::  NU_AQ  !Number of groups in each molecule
REAL, DIMENSION(NFUNC_AQ,NFUNC_AQ)  :: A_AQ   !Interaction parameters between different groups
!
!Parameters needed in unifac code which are precalculated in unifac_ini
!**********************************************************************
REAL, DIMENSION(NMOL_AQ,NFUNC_AQ)   :: THTAGP_AQ    !surface area ratio of groups in molecules
REAL, DIMENSION(NMOL_AQ)         :: Q_AQ            !surface parameters for molecules
REAL, DIMENSION(NMOL_AQ)         :: R_AQ            !volume parameter for molecules
REAL, DIMENSION(NMOL_AQ)         :: L_AQ            !unifac parameter for molecules

!***************************************************************************
!****************NUMBERS NEEDED FOR GRIFFINS ORGANIC  PHASE******************
!**************************************************************************

integer, parameter     ::  NMOL_ORG = NBSP + NBSPOA  !Total number of moles in organic phase
integer, parameter     ::  NFUNC_ORG = 16            !Number of functional groups
REAL, DIMENSION(NFUNC_ORG) :: RG_ORG                 !Group volume parameters
REAL, DIMENSION(NFUNC_ORG) :: QG_ORG                 !Group surface parameters
INTEGER, DIMENSION(NMOL_ORG, NFUNC_ORG) :: NU_ORG    !Number of groups in each molecule
REAL, DIMENSION(NFUNC_ORG,NFUNC_ORG) ::  A_ORG       !Interaction parameters between different groups

!Parameters needed in unifac code which are precalculated in unifac_ini
!**********************************************************************
REAL, DIMENSION(NMOL_ORG,NFUNC_ORG)   :: THTAGP_ORG    !surface area ratio of groups in molecules
REAL, DIMENSION(NMOL_ORG)         :: Q_ORG             !surface parameters for molecules
REAL, DIMENSION(NMOL_ORG)         :: R_ORG             !volume parameter for molecules
REAL, DIMENSION(NMOL_ORG)         :: L_ORG             !unifac parameter for molecules

!***************************************************************************
!****************NUMBERS NEEDED FOR PUN'S AQUOUS (A) PHASE******************
!**************************************************************************

integer, parameter ::  NMOL_A = NBSPA +1           !Number of molecules (the "+1" is for water)
integer, parameter ::  NFUNC_A = 11                !Number of functional groups
REAL, dimension(NFUNC_A) :: RG_A                   !Group volume parameters
REAL, DIMENSION(NFUNC_A) :: QG_A                   !Group surface parameters
INTEGER, DIMENSION(NMOL_A, NFUNC_A)  ::  NU_A      !Number of groups in each molecule
REAL, DIMENSION(NFUNC_A,NFUNC_A)  :: A_A           !Interaction parameters between different groups

!Parameters needed in unifac code which are precalculated in unifac_ini
!**********************************************************************
REAL, DIMENSION(NMOL_A,NFUNC_A)   :: THTAGP_A     !surface area ratio of groups in molecules
REAL, DIMENSION(NMOL_A)         :: Q_A            !surface parameters for molecules
REAL, DIMENSION(NMOL_A)         :: R_A            !volume parameter for molecules
REAL, DIMENSION(NMOL_A)         :: L_A            !unifac parameter for molecules

!***************************************************************************
!****************NUMBERS NEEDED FOR PUN'S ORGANIC (B) PHASE******************
!**************************************************************************
integer, parameter     ::  NMOL_B = 10              !Number of molecules
integer, parameter     ::  NFUNC_B = 16             !Number of functional groups
REAL, DIMENSION(NFUNC_B) :: RG_B                    !Group volume parameters
REAL, DIMENSION(NFUNC_B) :: QG_B                    !Group surface parameters
INTEGER, DIMENSION(NMOL_B, NFUNC_B) :: NU_B         !Number of groups in each molecule
REAL, DIMENSION(NFUNC_B,NFUNC_B) ::  A_B            !Interaction parameters between different groups

!Parameters needed in unifac code which are precalculated in unifac_ini
!**********************************************************************
REAL, DIMENSION(NMOL_B,NFUNC_B)   :: THTAGP_B    !surface area ratio of groups in molecules
REAL, DIMENSION(NMOL_B)         :: Q_B         !surface parameters for molecules
REAL, DIMENSION(NMOL_B)         :: R_B         !volume parameter for molecules
REAL, DIMENSION(NMOL_B)         :: L_B         !unifac parameter for molecules

END MODULE MODD_UNIFACPARAM
