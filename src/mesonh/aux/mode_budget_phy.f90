!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
! Modifications
!  P. Wautelet 28/01/2020: new SUBROUTINEs: Budget_store_init, Budget_store_end and Budget_source_id_find in new module mode_budget
!  P. Wautelet 17/08/2020: treat LES budgets correctly
!  P. Wautelet 05/03/2021: measure cpu_time for budgets
!-----------------------------------------------------------------

!#################
MODULE MODE_BUDGET_PHY
!#################

USE MODD_BUDGET, ONLY: TBUDGETDATA
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t

IMPLICIT NONE

PRIVATE

PUBLIC :: Budget_store_init_phy
PUBLIC :: Budget_store_end_phy
PUBLIC :: Budget_store_add_phy

CONTAINS

SUBROUTINE Budget_store_init_phy(D, tpbudget, hsource, pvars)
  USE MODE_BUDGET, ONLY: BUDGET_STORE_INIT
  TYPE(DIMPHYEX_t),       INTENT(IN)    :: D
  type(tbudgetdata),      intent(inout) :: tpbudget ! Budget datastructure
  character(len=*),       intent(in)    :: hsource  ! Name of the source term
  real, dimension(D%NIT,D%NJT,D%NKT), intent(in)    :: pvars    ! Current value to be stored
!
  CALL Budget_store_init(tpbudget, hsource, pvars)
!
END SUBROUTINE Budget_store_init_phy
!
SUBROUTINE Budget_store_end_phy(D, tpbudget, hsource, pvars)
  USE MODE_BUDGET, ONLY: BUDGET_STORE_END
  TYPE(DIMPHYEX_t),       INTENT(IN)    :: D
  type(tbudgetdata),      intent(inout) :: tpbudget ! Budget datastructure
  character(len=*),       intent(in)    :: hsource  ! Name of the source term
  real, dimension(D%NIT,D%NJT,D%NKT), intent(in)    :: pvars    ! Current value to be stored
!
  CALL Budget_store_end(tpbudget, hsource, pvars)
!
END SUBROUTINE Budget_store_end_phy
!
SUBROUTINE Budget_store_add_phy(D, tpbudget, hsource, pvars)
  USE MODE_BUDGET,   ONLY: BUDGET_STORE_ADD
  TYPE(DIMPHYEX_t),       INTENT(IN)    :: D
  type(tbudgetdata),      intent(inout) :: tpbudget ! Budget datastructure
  character(len=*),       intent(in)    :: hsource  ! Name of the source term
  real, dimension(D%NIT,D%NJT,D%NKT), intent(in)    :: pvars    ! Current value to be stored
!
  CALL Budget_store_add(tpbudget, hsource, pvars)
!
END SUBROUTINE Budget_store_add_phy
!
END MODULE MODE_BUDGET_PHY
