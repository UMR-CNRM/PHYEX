!     ######spl
module modd_aunifacparam

  use modd_glo

  implicit none

!/*************************************************************************

!Include file : unifacparma.h

!Purpose: Unifac parameters (replace file read)

!Include dependencies:  Included in Unidriver.c

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
!**************************************************************************/

!/* no. of molecules */
integer, parameter ::  NMOL_AQ = NBSP + 1

!/* no. of functional groups */
integer, parameter ::  NFUNC_AQ = 17

!/* Z = 10 is a fixed parameter in Unifac */
real, parameter    ::  Z = 10.0

!/* group volume parameters */
!/* dimension of RG is the same as NFUNC */
REAL, dimension(NFUNC_AQ) :: RG_AQ

!/* group surface area parameters */
!/* dimension of QG is the same as NFUNC */
REAL, DIMENSION(NFUNC_AQ) :: QG_AQ

!/* no. of groups in each molecule*/
INTEGER, DIMENSION(NMOL_AQ, NFUNC_AQ)  ::  NU_AQ

!Interaction parameters between different groups
REAL, DIMENSION(NFUNC_AQ,NFUNC_AQ)  :: A_AQ

!Parameters needed in unifac code which are precalculated in unifac_ini
!**********************************************************************
REAL, DIMENSION(NMOL_AQ,NFUNC_AQ)   :: THTAGP_AQ    !surface area ratio of groups in molecules
REAL, DIMENSION(NMOL_AQ)         :: Q_AQ            !surface parameters for molecules
REAL, DIMENSION(NMOL_AQ)         :: R_AQ            !volume parameter for molecules
REAL, DIMENSION(NMOL_AQ)         :: L_AQ            !unifac parameter for molecules

end module modd_aunifacparam
