!      #################################
       MODULE MODI_LIMA_CONVERSION_MELTING_SNOW
!      #################################
!
INTERFACE
   SUBROUTINE LIMA_CONVERSION_MELTING_SNOW (HFMFILE, OCLOSE_OUT, LDCOMPUTE,     &
                                            PRHODREF, PPRES, PT, PKA, PDV, PCJ, &
                                            PRVT, PRST, PLBDS,                  &
                                            P_RS_CMEL,                          &
                                            PA_RS, PA_RG                        )
!
CHARACTER(LEN=*),         INTENT(IN)    :: HFMFILE 
LOGICAL,                  INTENT(IN)    :: OCLOSE_OUT 
LOGICAL, DIMENSION(:),INTENT(IN)    :: LDCOMPUTE
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRHODREF ! Reference Exner function
REAL, DIMENSION(:),   INTENT(IN)    :: PPRES    !
REAL, DIMENSION(:),   INTENT(IN)    :: PT       !
REAL, DIMENSION(:),   INTENT(IN)    :: PKA      !
REAL, DIMENSION(:),   INTENT(IN)    :: PDV      !
REAL, DIMENSION(:),   INTENT(IN)    :: PCJ      !
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRVT    ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PRST    ! Cloud water C. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDS   ! 
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RS_CMEL
!
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_RS
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_RG
!
END SUBROUTINE LIMA_CONVERSION_MELTING_SNOW
END INTERFACE
END MODULE MODI_LIMA_CONVERSION_MELTING_SNOW
!
!     ######################################################################
      SUBROUTINE LIMA_CONVERSION_MELTING_SNOW (HFMFILE, OCLOSE_OUT, LDCOMPUTE,     &
                                               PRHODREF, PPRES, PT, PKA, PDV, PCJ, &
                                               PRVT, PRST, PLBDS,                  &
                                               P_RS_CMEL,                          &
                                               PA_RS, PA_RG                        )
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
USE MODD_CST,              ONLY : XTT, XMV, XMD, XLVTT, XCPV, XCL, XESTT, XRV
USE MODD_PARAM_LIMA,       ONLY : XRTMIN
USE MODD_PARAM_LIMA_MIXED, ONLY : XFSCVMG
USE MODD_PARAM_LIMA_COLD,  ONLY : X0DEPS, XEX0DEPS, X1DEPS, XEX1DEPS
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
REAL, DIMENSION(:),   INTENT(IN)    :: PPRES    !
REAL, DIMENSION(:),   INTENT(IN)    :: PT       !
REAL, DIMENSION(:),   INTENT(IN)    :: PKA      !
REAL, DIMENSION(:),   INTENT(IN)    :: PDV      !
REAL, DIMENSION(:),   INTENT(IN)    :: PCJ      !
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRVT    ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PRST    ! Cloud water C. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDS   ! 
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RS_CMEL
!
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_RS
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_RG
!
!*       0.2   Declarations of local variables :
!
REAL, DIMENSION(SIZE(PRST)) :: ZW ! work arrays
!
!-------------------------------------------------------------------------------
!
!
!*       1.     PRELIMINARY COMPUTATIONS
!	        ------------------------
!
!
P_RS_CMEL(:)=0.
!
ZW(:) = 0.0
WHERE( (PRST(:)>XRTMIN(5)) .AND. (PT(:)>XTT) .AND. LDCOMPUTE(:) )
   ZW(:) = PRVT(:)*PPRES(:)/((XMV/XMD)+PRVT(:)) ! Vapor pressure
   ZW(:) = PKA(:)*(XTT-PT(:)) +                                 &
              ( PDV(:)*(XLVTT + ( XCPV - XCL ) * ( PT(:) - XTT )) &
                          *(XESTT-ZW(:))/(XRV*PT(:))             )
!
! compute RSMLT
!
   ZW(:)  = XFSCVMG*MAX( 0.0,( -ZW(:) *             &
                          ( X0DEPS*           PLBDS(:)**XEX0DEPS +     &
                            X1DEPS*PCJ(:)*PLBDS(:)**XEX1DEPS ) ))!-    &
!!! BVIE
!!! ZZW1(1) et ZZW1(4) sont nuls là où PT>XTT !!!!!!!!!!!!!!
!!! On ne tient pas compte de la collection de pluie et gouttelettes par la neige si T>0 !!!! 
!                                    ( ZZW1(:,1)+ZZW1(:,4) ) *       &
!                             ( ZRHODREF(:)*XCL*(XTT-ZZT(:))) ) /    &
!                                            ( ZRHODREF(:)*XLMTT ) )
!
! note that RSCVMG = RSMLT*XFSCVMG but no heat is exchanged (at the rate RSMLT)
! because the graupeln produced by this process are still icy!!!
!
   P_RS_CMEL(:) = - ZW(:)
!
END WHERE
!
PA_RS(:) = PA_RS(:) + P_RS_CMEL(:)
PA_RG(:) = PA_RG(:) - P_RS_CMEL(:)
!
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_CONVERSION_MELTING_SNOW
