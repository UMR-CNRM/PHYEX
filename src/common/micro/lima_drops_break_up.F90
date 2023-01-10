!MNH_LIC Copyright 2013-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!      ###############################
       MODULE MODI_LIMA_DROPS_BREAK_UP
!      ###############################
!
INTERFACE
   SUBROUTINE LIMA_DROPS_BREAK_UP (LDCOMPUTE,  &
                                   PCRT, PRRT, &
                                   P_CR_BRKU,  &
                                   PB_CR       )

!
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
!
!
!     ##########################################
   SUBROUTINE LIMA_DROPS_BREAK_UP (LDCOMPUTE,  &
                                   PCRT, PRRT, &
                                   P_CR_BRKU,  &
                                   PB_CR       )

!     ##########################################
!
!!
!!    PURPOSE
!!    -------
!!      Numerical filter to prevent drops from growing too much
!!
!!    AUTHOR
!!    ------
!!      J.-M. Cohard     * Laboratoire d'Aerologie*
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!      S.    Berthet    * Laboratoire d'Aerologie*
!!      B.    Vié        * CNRM *
!!
!!
!!    MODIFICATIONS
!!    -------------
!!      Original             15/03/2018 
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
!              ---------------------------------------
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
