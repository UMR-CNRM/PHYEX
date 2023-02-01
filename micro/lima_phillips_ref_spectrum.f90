!      ######################################
       MODULE MODI_LIMA_PHILLIPS_REF_SPECTRUM
!      ######################################
!
INTERFACE
      SUBROUTINE LIMA_PHILLIPS_REF_SPECTRUM (ZZT, ZSI, ZSI_W, ZZY)
!
REAL, DIMENSION(:), INTENT(IN)    :: ZZT    ! Temperature
REAL, DIMENSION(:), INTENT(IN)    :: ZSI    ! Saturation over ice
REAL, DIMENSION(:), INTENT(IN)    :: ZSI_W  ! Saturation over ice at water sat.
REAL, DIMENSION(:), INTENT(INOUT) :: ZZY    ! Reference activity spectrum
!
END SUBROUTINE LIMA_PHILLIPS_REF_SPECTRUM
END INTERFACE
END MODULE MODI_LIMA_PHILLIPS_REF_SPECTRUM
!
!     ######################################################################
      SUBROUTINE LIMA_PHILLIPS_REF_SPECTRUM (ZZT, ZSI, ZSI_W, ZZY)
!     ######################################################################
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the reference activation spectrum
!!    described by Phillips (2008)
!!
!!
!!    REFERENCE
!!    ---------
!!
!!      Phillips et al., 2008: An empirical parameterization of heterogeneous
!!        ice nucleation for multiple chemical species of aerosols, J. Atmos. Sci. 
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
USE MODD_CST,             ONLY : XTT
USE MODD_PARAM_LIMA,      ONLY : XGAMMA, XRHO_CFDC
USE MODI_LIMA_FUNCTIONS,  ONLY : RECT, DELTA
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
REAL, DIMENSION(:), INTENT(IN)    :: ZZT    ! Temperature
REAL, DIMENSION(:), INTENT(IN)    :: ZSI    ! Saturation over ice
REAL, DIMENSION(:), INTENT(IN)    :: ZSI_W  ! Saturation over ice at water sat.
REAL, DIMENSION(:), INTENT(INOUT) :: ZZY    ! Reference activity spectrum
!
!*       0.2   Declarations of local variables :
!
REAL, DIMENSION(:), ALLOCATABLE   :: ZMAX, &
                                     ZMOY, &
                                     ZZY1, &
                                     ZZY2, &
                                     Z1,   &
                                     Z2,   &
                                     ZSI2
!
REAL                              :: XPSI
!
!-------------------------------------------------------------------------------
!
ALLOCATE(ZMAX(SIZE(ZZT))) ; ZMAX(:)= 0.0
ALLOCATE(ZMOY(SIZE(ZZT))) ; ZMOY(:)= 0.0
ALLOCATE(ZZY1(SIZE(ZZT))) ; ZZY1(:)= 0.0
ALLOCATE(ZZY2(SIZE(ZZT))) ; ZZY2(:)= 0.0
ALLOCATE(Z1(SIZE(ZZT)))   ; Z1(:)  = 0.0
ALLOCATE(Z2(SIZE(ZZT)))   ; Z2(:)  = 0.0
ALLOCATE(ZSI2(SIZE(ZZT))) ; ZSI2(:)= 0.0
!
ZZY(:) = 0.0   
!
XPSI   = 0.058707*XGAMMA/XRHO_CFDC
!
ZSI2(:)=min(ZSI(:),ZSI_W(:))
!
WHERE( ZSI(:)>1.0 )
!
!* T <= -35 C 
!
   ZZY(:)  =1000.*XGAMMA/XRHO_CFDC                  &
        * ( EXP(12.96*(MIN(ZSI2(:),7.)-1.1)) )**0.3          &
        * RECT(1.,0.,ZZT(:),(XTT-80.),(XTT-35.))
!
!* -35 C < T <= -25 C (in Appendix A) 
!
   ZZY1(:) =1000.*XGAMMA/XRHO_CFDC                  &
        * ( EXP(12.96*(MIN(ZSI2(:),7.)-1.1)) )**0.3
   ZZY2(:) =1000.*XPSI                              &
        *   EXP(12.96*(MIN(ZSI2(:),7.)-1.0)-0.639)
!
!* -35 C < T <= -30 C
!
   ZMAX(:) =1000.*XGAMMA/XRHO_CFDC                  &
        * ( EXP(12.96*(ZSI_W(:)-1.1)) )**0.3        &
        * RECT(1.,0.,ZZT(:),(XTT-35.),(XTT-30.))
!
!* -30 C < T <= -25 C
!
   ZMAX(:) = ZMAX(:) +1000.*XPSI                    &
        * EXP( 12.96*(ZSI_W(:)-1.0)-0.639 )         &
        * RECT(1.,0.,ZZT(:),(XTT-30.),(XTT-25.))
   Z1(:)   = MIN(ZZY1(:), ZMAX(:)) 
   Z2(:)   = MIN(ZZY2(:), ZMAX(:)) 
!
!* T > -25 C 
!
   ZZY(:)  = ZZY(:) + 1000.*XPSI                    &
        * EXP( 12.96*(MIN(ZSI2(:),7.)-1.0)-0.639 )           &
        * RECT(1.,0.,ZZT(:),(XTT-25.),(XTT-2.))
END WHERE
!
WHERE (Z2(:)>0.0 .AND. Z1(:)>0.0)
   ZMOY(:) = Z2(:)*(Z1(:)/Z2(:))**DELTA(1.,0.,ZZT(:),(XTT-35.),(XTT-25.))
   ZZY(:)  = ZZY(:) + MIN(ZMOY(:),ZMAX(:))  ! N_{IN,1,*}
END WHERE
!
!++cb++
DEALLOCATE(ZMAX)
DEALLOCATE(ZMOY)
DEALLOCATE(ZZY1)
DEALLOCATE(ZZY2)
DEALLOCATE(Z1)
DEALLOCATE(Z2)
!--cb--
!
END SUBROUTINE LIMA_PHILLIPS_REF_SPECTRUM
