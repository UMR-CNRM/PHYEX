!     ######spl
      MODULE MODD_BUDGET
!     ##################
!
!!****  *MODD_BUDGET* - declaration of budget variables
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to specify  the budget
!     variables.
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      MODD_PARAMETERS: JPBUMAX, JPBUPROCMAX
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (module MODD_BUDGET)
!!
!!    AUTHOR
!!    ------
!!      P. Hereil   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original        23/02/95
!!      J.-P. Lafore    10/02/98    adding of rhodj declaration for budget
!!      V. Ducrocq      4/06/99     //
!!      J.-P. Pinty     25/09/00    additional budget terms for C2R2 scheme
!!      D. Gazen        22/01/01    add NCHEMSV
!!      V. Masson       06/11/02    new flags for budget calls and time counters
!!      V. Masson       27/11/02    add 2way nesting effect
!!      P. Jabouille    07/07/04    add budget terms for microphysics
!!      C. Barthe       19/11/09    add budget terms for electricity
!!      S. Riette       July 2022   simplification for PHYEX
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE

SAVE
!
INTEGER, PARAMETER:: NBUDGET_RHO = 0  ! Reference number for budget of RhoJ
INTEGER, PARAMETER:: NBUDGET_U   = 1  ! Reference number for budget of RhoJu  and/or LES budgets with u
INTEGER, PARAMETER:: NBUDGET_V   = 2  ! Reference number for budget of RhoJv  and/or LES budgets with u
INTEGER, PARAMETER:: NBUDGET_W   = 3  ! Reference number for budget of RhoJw  and/or LES budgets with u
INTEGER, PARAMETER:: NBUDGET_TH  = 4  ! Reference number for budget of RhoJTh and/or LES budgets with th
INTEGER, PARAMETER:: NBUDGET_TKE = 5  ! Reference number for budget of RhoJTke and/or LES budgets with Tke
INTEGER, PARAMETER:: NBUDGET_RV  = 6  ! Reference number for budget of RhoJrv and/or LES budgets with rv
INTEGER, PARAMETER:: NBUDGET_RC  = 7  ! Reference number for budget of RhoJrc and/or LES budgets with rc
INTEGER, PARAMETER:: NBUDGET_RR  = 8  ! Reference number for budget of RhoJrr and/or LES budgets with rr
INTEGER, PARAMETER:: NBUDGET_RI  = 9  ! Reference number for budget of RhoJri and/or LES budgets with ri
INTEGER, PARAMETER:: NBUDGET_RS  = 10 ! Reference number for budget of RhoJrs and/or LES budgets with rs
INTEGER, PARAMETER:: NBUDGET_RG  = 11 ! Reference number for budget of RhoJrg and/or LES budgets with rg
INTEGER, PARAMETER:: NBUDGET_RH  = 12 ! Reference number for budget of RhoJrh and/or LES budgets with rh
INTEGER, PARAMETER:: NBUDGET_SV1 = 13 ! Reference number for 1st budget of RhoJsv and/or LES budgets with sv
!
TYPE TBUDGETDATA
  INTEGER :: NBUDGET
ENDTYPE
!
!                       General variables
LOGICAL :: LBU_ENABLE=.FALSE.
!
INTEGER :: NBUMOD=0                    ! model in which budget is calculated
!
LOGICAL :: LBUDGET_U=.FALSE.  ! flag to compute budget of RhoJu  and/or LES budgets with u
LOGICAL :: LBUDGET_V=.FALSE.  ! flag to compute budget of RhoJv  and/or LES budgets with u
LOGICAL :: LBUDGET_W=.FALSE.  ! flag to compute budget of RhoJw  and/or LES budgets with u
LOGICAL :: LBUDGET_TH=.FALSE. ! flag to compute budget of RhoJTh and/or LES budgets with th
LOGICAL :: LBUDGET_TKE=.FALSE.! flag to compute budget of RhoJTke and/or LES budgets with Tke
LOGICAL :: LBUDGET_RV=.FALSE. ! flag to compute budget of RhoJrv and/or LES budgets with rv
LOGICAL :: LBUDGET_RC=.FALSE. ! flag to compute budget of RhoJrc and/or LES budgets with rc
LOGICAL :: LBUDGET_RR=.FALSE. ! flag to compute budget of RhoJrr and/or LES budgets with rr
LOGICAL :: LBUDGET_RI=.FALSE. ! flag to compute budget of RhoJri and/or LES budgets with ri
LOGICAL :: LBUDGET_RS=.FALSE. ! flag to compute budget of RhoJrs and/or LES budgets with rs
LOGICAL :: LBUDGET_RG=.FALSE. ! flag to compute budget of RhoJrg and/or LES budgets with rg
LOGICAL :: LBUDGET_RH=.FALSE. ! flag to compute budget of RhoJrh and/or LES budgets with rh
LOGICAL :: LBUDGET_SV=.FALSE. ! flag to compute budget of RhoJsv and/or LES budgets with sv
!
END MODULE MODD_BUDGET
