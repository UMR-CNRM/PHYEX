MODULE SHUMAN_PHY
IMPLICIT NONE
CONTAINS
!     ###############################
      SUBROUTINE MZM_PHY(D,PA,PMZM)
!     ###############################
!
!!****  *MZM* -  Shuman operator : mean operator in z direction for a 
!!                                 mass variable 
!!
!!    PURPOSE
!!    -------
!       The purpose of this function  is to compute a mean
!     along the z direction (K index) for a field PA localized at a mass
!     point. The result is localized at a z-flux point (w point).
!
!!**  METHOD
!!    ------ 
!!        The result PMZM(:,:,k) is defined by 0.5*(PA(:,:,k)+PA(:,:,k-1))
!!        At k=1, PMZM(:,:,1) is defined by -999.
!!    
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      NONE
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (SHUMAN operators)
!!      Technical specifications Report of The Meso-NH (chapters 3)  
!!
!!
!!    AUTHOR
!!    ------
!!	V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/07/94 
!!                   optimisation                 20/08/00 J. Escobar
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
IMPLICIT NONE
!
!*       0.1   Declarations of argument and result
!              ------------------------------------
!
TYPE(DIMPHYEX_t),       INTENT(IN)  :: D
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)  :: PA     ! variable at mass localization
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT) :: PMZM   ! result at flux localization 

!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: JK             ! Loop index in z direction
INTEGER :: IKU            ! upper bound in z direction of PA
!           
INTEGER :: IIU,IJU
INTEGER :: JIJ,JI,JJ
INTEGER :: JIJK,JIJKOR,JIJKEND
!
!
!-------------------------------------------------------------------------------
!
!*       1.    DEFINITION OF MZM
!              ------------------
!
IIU = SIZE(PA,1)
IJU = SIZE(PA,2)
IKU = SIZE(PA,3)
!
JIJKOR  = 1 + IIU*IJU
JIJKEND = IIU*IJU*IKU
!
!CDIR NODEP
!OCL NOVREC
DO JIJK=JIJKOR , JIJKEND
   PMZM(JIJK,1,1) = 0.5*( PA(JIJK,1,1)+PA(JIJK-IIU*IJU,1,1) )
END DO
!
!CDIR NODEP
!OCL NOVREC
DO JIJ=1,IIU*IJU
   PMZM(JIJ,1,1)    = -999.
END DO
!-------------------------------------------------------------------------------
!
END SUBROUTINE MZM_PHY
END MODULE SHUMAN_PHY
