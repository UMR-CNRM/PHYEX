MODULE MODE_LIMA_PHILLIPS_REF_SPECTRUM
  IMPLICIT NONE
CONTAINS
!     ######################################################################
  SUBROUTINE LIMA_PHILLIPS_REF_SPECTRUM (LIMAP, CST, ISIZE, PZT, PSI, PSI_W, PZY)
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
USE MODD_CST,            ONLY: CST_t
USE MODD_PARAM_LIMA,      ONLY : XGAMMA, XRHO_CFDC
USE MODE_LIMA_FUNCTIONS,  ONLY : RECT, DELTA
USE MODD_PARAM_LIMA, ONLY:PARAM_LIMA_t
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(CST_t),            INTENT(IN)    :: CST
INTEGER,                INTENT(IN)    :: ISIZE
REAL, DIMENSION(ISIZE), INTENT(IN)    :: PZT    ! Temperature
REAL, DIMENSION(ISIZE), INTENT(IN)    :: PSI    ! Saturation over ice
REAL, DIMENSION(ISIZE), INTENT(IN)    :: PSI_W  ! Saturation over ice at water sat.
REAL, DIMENSION(ISIZE), INTENT(OUT)   :: PZY    ! Reference activity spectrum
!
!*       0.2   Declarations of local variables :
!
REAL, DIMENSION(ISIZE)   :: ZMAX, &
                            ZMOY, &
                            ZZY1, &
                            ZZY2, &
                            Z1,   &
                            Z2,   &
                            ZSI2
!
TYPE(PARAM_LIMA_t),INTENT(IN)::LIMAP
REAL                     :: ZPSI
!
!-------------------------------------------------------------------------------
!
ZMAX(:)= 0.0
ZMOY(:)= 0.0
ZZY1(:)= 0.0
ZZY2(:)= 0.0
Z1(:)  = 0.0
Z2(:)  = 0.0
ZSI2(:)= 0.0
!
PZY(:) = 0.0   
!
ZPSI   = 0.058707*XGAMMA/XRHO_CFDC
!
ZSI2(:)=min(PSI(:),PSI_W(:))
!
WHERE( ZSI2(:)>1.0 )
!
!* T <= -35 C 
!
   PZY(:)  =1000.*XGAMMA/XRHO_CFDC                  &
        * ( EXP(12.96*(MIN(ZSI2(:),7.)-1.1)) )**0.3          &
        * RECT(1.,0.,PZT(:),(CST%XTT-80.),(CST%XTT-35.))
!
!* -35 C < T <= -25 C (in Appendix A) 
!
   ZZY1(:) =1000.*XGAMMA/XRHO_CFDC                  &
        * ( EXP(12.96*(MIN(ZSI2(:),7.)-1.1)) )**0.3
   ZZY2(:) =1000.*ZPSI                              &
        *   EXP(12.96*(MIN(ZSI2(:),7.)-1.0)-0.639)
!
!* -35 C < T <= -30 C
!
   ZMAX(:) =1000.*XGAMMA/XRHO_CFDC                  &
        * ( EXP(12.96*(PSI_W(:)-1.1)) )**0.3        &
        * RECT(1.,0.,PZT(:),(CST%XTT-35.),(CST%XTT-30.))
!
!* -30 C < T <= -25 C
!
   ZMAX(:) = ZMAX(:) +1000.*ZPSI                    &
        * EXP( 12.96*(PSI_W(:)-1.0)-0.639 )         &
        * RECT(1.,0.,PZT(:),(CST%XTT-30.),(CST%XTT-25.))
   Z1(:)   = MIN(ZZY1(:), ZMAX(:)) 
   Z2(:)   = MIN(ZZY2(:), ZMAX(:)) 
!
!* T > -25 C 
!
   PZY(:)  = PZY(:) + 1000.*ZPSI                    &
        * EXP( 12.96*(MIN(ZSI2(:),7.)-1.0)-0.639 )           &
        * RECT(1.,0.,PZT(:),(CST%XTT-25.),(CST%XTT-2.))
END WHERE
!
WHERE (Z2(:)>0.0 .AND. Z1(:)>0.0)
   ZMOY(:) = Z2(:)*(Z1(:)/Z2(:))**DELTA(1.,0.,PZT(:),(CST%XTT-35.),(CST%XTT-25.))
   PZY(:)  = PZY(:) + MIN(ZMOY(:),ZMAX(:))  ! N_{IN,1,*}
END WHERE
!
END SUBROUTINE LIMA_PHILLIPS_REF_SPECTRUM
END MODULE MODE_LIMA_PHILLIPS_REF_SPECTRUM
