!      #################################
       MODULE MODI_LIMA_DROPS_SELF_COLLECTION
!      #################################
!
INTERFACE
   SUBROUTINE LIMA_DROPS_SELF_COLLECTION (HFMFILE, OCLOSE_OUT, LDCOMPUTE, &
                                             PRHODREF,                       &
                                             PCRT, PLBDR, PLBDR3,            &
                                             P_CR_SCBU,                      &
                                             PA_CR                           )
!
CHARACTER(LEN=*),         INTENT(IN)    :: HFMFILE 
LOGICAL,                  INTENT(IN)    :: OCLOSE_OUT 
LOGICAL, DIMENSION(:),INTENT(IN)    :: LDCOMPUTE
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRHODREF ! Reference Exner function
!
REAL, DIMENSION(:),   INTENT(IN)    :: PCRT    ! Cloud water C. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDR   ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDR3  ! 
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CR_SCBU
!
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_CR
!
END SUBROUTINE LIMA_DROPS_SELF_COLLECTION
END INTERFACE
END MODULE MODI_LIMA_DROPS_SELF_COLLECTION
!
!     ######################################################################
      SUBROUTINE LIMA_DROPS_SELF_COLLECTION (HFMFILE, OCLOSE_OUT, LDCOMPUTE, &
                                             PRHODREF,                       &
                                             PCRT, PLBDR, PLBDR3,            &
                                             P_CR_SCBU,                      &
                                             PA_CR                           )
!     ######################################################################
!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the cold-phase homogeneous
!!    freezing of CCN, droplets and drops (T<-35°C)
!!
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
!!      C. Barthe  * LACy*   jan. 2014  add budgets
!!      B.Vie 10/2016 Bug zero division
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAM_LIMA,      ONLY : XCTMIN
USE MODD_PARAM_LIMA_WARM, ONLY : XACCR1, XSCBUEXP1, XSCBU_EFF1, XSCBU_EFF2, &
                                 XSCBU2, XSCBU3
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
CHARACTER(LEN=*),         INTENT(IN)    :: HFMFILE 
LOGICAL,                  INTENT(IN)    :: OCLOSE_OUT 
LOGICAL, DIMENSION(:),INTENT(IN)    :: LDCOMPUTE
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRHODREF ! Reference Exner function
!
REAL, DIMENSION(:),   INTENT(IN)    :: PCRT    ! Cloud water C. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDR   ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDR3  ! 
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CR_SCBU
!
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_CR
!
!*       0.2   Declarations of local variables :
!
REAL, DIMENSION(SIZE(PCRT)) :: &
                                           ZW1, & ! work arrays
                                           ZW2, &
                                           ZW3, &
                                           ZW4, &
                                           ZSCBU
!
!-------------------------------------------------------------------------------
!
!
!*       1.     PRELIMINARY COMPUTATIONS
!	        ------------------------
!
!
P_CR_SCBU(:)=0.
!
WHERE( PCRT(:)>XCTMIN(3) .AND. LDCOMPUTE(:) )
   ZW4(:)  = XACCR1 / PLBDR(:)                ! Mean diameter
END WHERE
ZSCBU(:)=0.
WHERE (ZW4(:)>=XSCBU_EFF1 .AND. PCRT(:)>XCTMIN(3) .AND. LDCOMPUTE(:)) &
     ZSCBU(:) = EXP(XSCBUEXP1*(ZW4(:)-XSCBU_EFF1))  ! coalescence efficiency
WHERE (ZW4(:)>=XSCBU_EFF2 .AND. LDCOMPUTE(:)) ZSCBU(:) = 0.0  ! Break-up
!
ZW1(:) = 0.0
ZW2(:) = 0.0
ZW3(:) = 0.0
!
WHERE (PCRT(:)>XCTMIN(3) .AND. (ZW4(:)>1.E-4) .AND. LDCOMPUTE(:))  ! analytical integration
   ZW1(:) = XSCBU2 * PCRT(:)**2 / PLBDR3(:)   ! D>100 10-6 m
   ZW3(:) = ZW1(:)*ZSCBU(:)
END WHERE
!
WHERE (PCRT(:)>XCTMIN(3) .AND. (ZW4(:)<=1.E-4) .AND. LDCOMPUTE(:))
   ZW2(:) = XSCBU3 *(PCRT(:) / PLBDR3(:))**2  ! D<100 10-6 m
   ZW3(:) = ZW2(:)
END WHERE
!
P_CR_SCBU(:) = - ZW3(:) * PRHODREF(:)
!
PA_CR(:) = PA_CR(:) + P_CR_SCBU(:)
!
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_DROPS_SELF_COLLECTION
