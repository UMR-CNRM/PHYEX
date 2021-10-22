!     ######spl
module modd_bunifacparam

  use modd_glo

  implicit none

  !************************************************************************
  !Purpose: Unifac parameters Type B (replace file read)

  !Notes: 5 primary compounds (EPRI 99): C23H47COOH, C8H17CH=CHC7H14COOH,
  !     4-(2-propio)-syringone, C29H60, 2-carboxybenzoic acid
  !     + 5 Type B compounds at 25C
  !     (Data of Lyman, Reehl, Rosenblatt, 1990)
  !
  !     NMOL and NFUNC need to match DIMMOL and DIMFUN (used in unidriver.c)
  !
  !     Parameters inputted as in the original fortran input file, therefore
  !     a transpose is needed in the unidriver C program for matrices A and NU
  !     in order to pass them properly into the Fortran unifac routine.
  !
  !revision History:  1. Developed by Betty Pun, AER, December, 1999
  !                    under CARB funding
  !                 2. Modified by Betty Pun, AER, December, 1999
  !                    under CARB funding for Type B module
  !************************************************************************


! no. of molecules
integer, parameter     ::  NMOL_ORG = NBSP + NBSPOA

! no. of functional groups
integer, parameter     ::  NFUNC_ORG = 16

!/* Z = 10 is a fixed parameter in Unifac */
real, parameter        :: Z = 10.0

!/* group volume parameters */
!/* dimension of RG is the same as NFUNC */
REAL, DIMENSION(NFUNC_ORG) :: RG_ORG

!/* group surface area parameters */
!/* dimension of QG is the same as NFUNC */
REAL, DIMENSION(NFUNC_ORG) :: QG_ORG

!/* no. of groups in each molecule*/
INTEGER, DIMENSION(NMOL_ORG, NFUNC_ORG) :: NU_ORG

!Interaction parameters between different groups
REAL, DIMENSION(NFUNC_ORG,NFUNC_ORG) ::  A_ORG

!Parameters needed in unifac code which are precalculated in unifac_ini
!**********************************************************************
REAL, DIMENSION(NMOL_ORG,NFUNC_ORG)   :: THTAGP_ORG    !surface area ratio of groups in molecules
REAL, DIMENSION(NMOL_ORG)         :: Q_ORG         !surface parameters for molecules
REAL, DIMENSION(NMOL_ORG)         :: R_ORG         !volume parameter for molecules
REAL, DIMENSION(NMOL_ORG)         :: L_ORG         !unifac parameter for molecules

end module modd_bunifacparam
