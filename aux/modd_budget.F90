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
TYPE, ABSTRACT :: TBUDGETDATA
  INTEGER :: NBUDGET
  CONTAINS
  PROCEDURE (TBUDGETDATA_STORE_INIT),     DEFERRED :: INIT
  PROCEDURE (TBUDGETDATA_STORE_INIT_PHY), DEFERRED :: INIT_PHY
  PROCEDURE (TBUDGETDATA_STORE_END),      DEFERRED :: END
  PROCEDURE (TBUDGETDATA_STORE_END_PHY),  DEFERRED :: END_PHY
  PROCEDURE (TBUDGETDATA_STORE_ADD),      DEFERRED :: ADD
  PROCEDURE (TBUDGETDATA_STORE_ADD_PHY),  DEFERRED :: ADD_PHY
ENDTYPE TBUDGETDATA

ABSTRACT INTERFACE
SUBROUTINE TBUDGETDATA_STORE_INIT(SELF, HSOURCE, PVARS)
  IMPORT TBUDGETDATA
  CLASS(TBUDGETDATA),      INTENT(INOUT) :: SELF ! Budget datastructure
  CHARACTER(LEN=*),       INTENT(IN)    :: HSOURCE  ! Name of the source term
  REAL, DIMENSION(:,:), INTENT(IN)    :: PVARS    ! Current value to be stored
END SUBROUTINE TBUDGETDATA_STORE_INIT
!
SUBROUTINE TBUDGETDATA_STORE_INIT_PHY(SELF, D, HSOURCE, PVARS)
  USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
  IMPORT TBUDGETDATA
  CLASS(TBUDGETDATA),      INTENT(INOUT) :: SELF ! Budget datastructure
  TYPE(DIMPHYEX_t),       INTENT(IN)    :: D
  CHARACTER(LEN=*),       INTENT(IN)    :: HSOURCE  ! Name of the source term
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PVARS    ! Current value to be stored
END SUBROUTINE TBUDGETDATA_STORE_INIT_PHY
!
SUBROUTINE TBUDGETDATA_STORE_END(SELF, HSOURCE, PVARS)
  IMPORT TBUDGETDATA
  CLASS(TBUDGETDATA),      INTENT(INOUT) :: SELF ! Budget datastructure
  CHARACTER(LEN=*),       INTENT(IN)    :: HSOURCE  ! Name of the source term
  REAL, DIMENSION(:,:), INTENT(IN)    :: PVARS    ! Current value to be stored
END SUBROUTINE TBUDGETDATA_STORE_END
!
SUBROUTINE TBUDGETDATA_STORE_END_PHY(SELF, D, HSOURCE, PVARS)
  USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
  IMPORT TBUDGETDATA
  CLASS(TBUDGETDATA),      INTENT(INOUT) :: SELF ! Budget datastructure
  TYPE(DIMPHYEX_t),       INTENT(IN)    :: D
  CHARACTER(LEN=*),       INTENT(IN)    :: HSOURCE  ! Name of the source term
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PVARS    ! Current value to be stored
END SUBROUTINE TBUDGETDATA_STORE_END_PHY
!
SUBROUTINE TBUDGETDATA_STORE_ADD_PHY(SELF, D, HSOURCE, PVARS)
  USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
  IMPORT TBUDGETDATA
  CLASS(TBUDGETDATA),      INTENT(INOUT) :: SELF ! Budget datastructure
  TYPE(DIMPHYEX_t),       INTENT(IN)    :: D
  CHARACTER(LEN=*),       INTENT(IN)    :: HSOURCE  ! Name of the source term
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PVARS    ! Current value to be stored
END SUBROUTINE TBUDGETDATA_STORE_ADD_PHY
!
SUBROUTINE TBUDGETDATA_STORE_ADD(SELF, HSOURCE, PVARS)
  IMPORT TBUDGETDATA
  CLASS(TBUDGETDATA),      INTENT(INOUT) :: SELF ! Budget datastructure
  CHARACTER(LEN=*),       INTENT(IN)    :: HSOURCE  ! Name of the source term
  REAL, DIMENSION(:,:), INTENT(IN)    :: PVARS    ! Current value to be stored
END SUBROUTINE TBUDGETDATA_STORE_ADD

END INTERFACE

TYPE TBUDGETDATA_PTR
  CLASS (TBUDGETDATA), POINTER :: PTR => NULL()
END TYPE TBUDGETDATA_PTR
!
TYPE TBUDGETCONF_t
  LOGICAL :: LBU_ENABLE=.FALSE.
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
END TYPE TBUDGETCONF_t
!
TYPE(TBUDGETCONF_t), TARGET :: TBUCONF
!
!                       General variables
LOGICAL, POINTER :: LBU_ENABLE=>NULL()
!
INTEGER :: NBUMOD=0                    ! model in which budget is calculated
!
LOGICAL, POINTER :: LBUDGET_U=>NULL()    ! flag to compute budget of RhoJu  and/or LES budgets with u
LOGICAL, POINTER :: LBUDGET_V=>NULL()    ! flag to compute budget of RhoJv  and/or LES budgets with u
LOGICAL, POINTER :: LBUDGET_W=>NULL()    ! flag to compute budget of RhoJw  and/or LES budgets with u
LOGICAL, POINTER :: LBUDGET_TH=>NULL()   ! flag to compute budget of RhoJTh and/or LES budgets with th
LOGICAL, POINTER :: LBUDGET_TKE=>NULL()  ! flag to compute budget of RhoJTke and/or LES budgets with Tke
LOGICAL, POINTER :: LBUDGET_RV=>NULL()   ! flag to compute budget of RhoJrv and/or LES budgets with rv
LOGICAL, POINTER :: LBUDGET_RC=>NULL()   ! flag to compute budget of RhoJrc and/or LES budgets with rc
LOGICAL, POINTER :: LBUDGET_RR=>NULL()   ! flag to compute budget of RhoJrr and/or LES budgets with rr
LOGICAL, POINTER :: LBUDGET_RI=>NULL()   ! flag to compute budget of RhoJri and/or LES budgets with ri
LOGICAL, POINTER :: LBUDGET_RS=>NULL()   ! flag to compute budget of RhoJrs and/or LES budgets with rs
LOGICAL, POINTER :: LBUDGET_RG=>NULL()   ! flag to compute budget of RhoJrg and/or LES budgets with rg
LOGICAL, POINTER :: LBUDGET_RH=>NULL()   ! flag to compute budget of RhoJrh and/or LES budgets with rh
LOGICAL, POINTER :: LBUDGET_SV=>NULL()   ! flag to compute budget of RhoJsv and/or LES budgets with sv

CONTAINS

SUBROUTINE TBUCONF_ASSOCIATE()
  IMPLICIT NONE
  LBU_ENABLE=>TBUCONF%LBU_ENABLE

  LBUDGET_U=>TBUCONF%LBUDGET_U
  LBUDGET_V=>TBUCONF%LBUDGET_V
  LBUDGET_W=>TBUCONF%LBUDGET_W
  LBUDGET_TH=>TBUCONF%LBUDGET_TH
  LBUDGET_TKE=>TBUCONF%LBUDGET_TKE
  LBUDGET_RV=>TBUCONF%LBUDGET_RV
  LBUDGET_RC=>TBUCONF%LBUDGET_RC
  LBUDGET_RR=>TBUCONF%LBUDGET_RR
  LBUDGET_RI=>TBUCONF%LBUDGET_RI
  LBUDGET_RS=>TBUCONF%LBUDGET_RS
  LBUDGET_RG=>TBUCONF%LBUDGET_RG
  LBUDGET_RH=>TBUCONF%LBUDGET_RH
  LBUDGET_SV=>TBUCONF%LBUDGET_SV
END SUBROUTINE TBUCONF_ASSOCIATE
!
END MODULE MODD_BUDGET
