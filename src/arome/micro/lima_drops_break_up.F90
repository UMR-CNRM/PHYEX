!      ###############################
       MODULE MODI_LIMA_DROPS_BREAK_UP
!      ###############################
!
INTERFACE
   SUBROUTINE LIMA_DROPS_BREAK_UP (HFMFILE, OCLOSE_OUT, LDCOMPUTE,    &
                                   PCRT, PRRT,                        &
                                   P_CR_BRKU,                         &
                                   PB_CR                              )

!
CHARACTER(LEN=*),      INTENT(IN)    :: HFMFILE    
LOGICAL,               INTENT(IN)    :: OCLOSE_OUT 
LOGICAL, DIMENSION(:), INTENT(IN)    :: LDCOMPUTE  
!
REAL, DIMENSION(:),    INTENT(IN)    :: PCRT             !
REAL, DIMENSION(:),    INTENT(IN)    :: PRRT             !
!
REAL, DIMENSION(:),    INTENT(INOUT) :: P_CR_BRKU        ! Concentration change (#/kg)
REAL, DIMENSION(:),    INTENT(INOUT) :: PB_CR            ! Cumulated concentration change (#/kg)
!
END SUBROUTINE LIMA_DROPS_BREAK_UP
END INTERFACE
END MODULE MODI_LIMA_DROPS_BREAK_UP
!     #####################################################################
   SUBROUTINE LIMA_DROPS_BREAK_UP (HFMFILE, OCLOSE_OUT, LDCOMPUTE,    &
                                   PCRT, PRRT,                        &
                                   P_CR_BRKU,                         &
                                   PB_CR                              )

!     #####################################################################
!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the warm microphysical 
!!    sources: nucleation, sedimentation, autoconversion, accretion,  
!!    self-collection and vaporisation which are parameterized according  
!!    to Cohard and Pinty, QJRMS, 2000
!!
!!
!!**  METHOD
!!    ------
!!      The activation of CCN is checked for quasi-saturated air parcels 
!!    to update the cloud droplet number concentration. Then assuming a 
!!    generalized gamma distribution law for the cloud droplets and the 
!!    raindrops, the zeroth and third order moments tendencies are evaluated
!!    for all the coalescence terms by integrating the Stochastic Collection 
!!    Equation. As autoconversion is a process that cannot be resolved 
!!    analytically, the Berry-Reinhardt parameterisation is employed with
!!    modifications to initiate the raindrop spectrum mode. The integration
!!    of the raindrop evaporation below clouds is straightforward.
!!
!!      The sedimentation rates are computed with a time spliting technique: 
!!    an upstream scheme, written as a difference of non-advective fluxes. 
!!    This source term is added to the next coming time step (split-implicit 
!!    process).
!!
!!    REFERENCE
!!    ---------
!!
!!      Cohard, J.-M. and J.-P. Pinty, 2000: A comprehensive two-moment warm 
!!      microphysical bulk scheme. 
!!        Part I: Description and tests
!!        Part II: 2D experiments with a non-hydrostatic model
!!      Accepted for publication in Quart. J. Roy. Meteor. Soc. 
!!
!!    AUTHOR
!!    ------
!!      J.-M. Cohard     * Laboratoire d'Aerologie*
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!      S.    Berthet    * Laboratoire d'Aerologie*
!!      B.    Vié        * Laboratoire d'Aerologie*
!!
!!
!!    MODIFICATIONS
!!    -------------
!!      Original             ??/??/13 
!!      C. Barthe  * LACy *  jan. 2014   add budgets
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAM_LIMA,      ONLY : XCTMIN, XRTMIN
USE MODD_PARAM_LIMA_WARM, ONLY : XACCR1, XLBEXR, XLBR, XSPONBUD1, XSPONBUD3, XSPONCOEF2
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
CHARACTER(LEN=*),      INTENT(IN)    :: HFMFILE    
LOGICAL,               INTENT(IN)    :: OCLOSE_OUT 
LOGICAL, DIMENSION(:), INTENT(IN)    :: LDCOMPUTE  
!
REAL, DIMENSION(:),    INTENT(IN)    :: PCRT             !
REAL, DIMENSION(:),    INTENT(IN)    :: PRRT             !
!
REAL, DIMENSION(:),    INTENT(INOUT) :: P_CR_BRKU        ! Concentration change (#/kg)
REAL, DIMENSION(:),    INTENT(INOUT) :: PB_CR            ! Cumulated concentration change (#/kg)
!
!*       0.2   Declarations of local variables :
!
REAL,    DIMENSION(SIZE(PCRT)) :: ZWLBDR,ZWLBDR3
INTEGER :: JL
!
!-------------------------------------------------------------------------------
!
!              SPONTANEOUS BREAK-UP (NUMERICAL FILTER)
!              --------------------
!
P_CR_BRKU(:)=0.
!
ZWLBDR3(:) = 1.E30
ZWLBDR(:) = 1.E10
WHERE ( PRRT(:)>XRTMIN(3) .AND. PCRT(:)>XCTMIN(3) .AND. LDCOMPUTE(:) )
   ZWLBDR3(:) = XLBR * PCRT(:) / PRRT(:)
   ZWLBDR(:)  = ZWLBDR3(:)**XLBEXR
END WHERE
WHERE (ZWLBDR(:)<(XACCR1/XSPONBUD1) .AND. LDCOMPUTE(:))
   P_CR_BRKU(:) = PCRT(:)*( MAX((1.+XSPONCOEF2*(XACCR1/ZWLBDR(:)-XSPONBUD1)**2),&
                                                     (XACCR1/ZWLBDR(:)/XSPONBUD3)**3) -1. )
END WHERE
!
PB_CR(:) = PB_CR(:) + P_CR_BRKU(:)
!-------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_DROPS_BREAK_UP
