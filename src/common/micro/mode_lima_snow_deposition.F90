!MNH_LIC Copyright 2018-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_LIMA_SNOW_DEPOSITION
  IMPLICIT NONE
CONTAINS
!     ##########################################################################
  SUBROUTINE LIMA_SNOW_DEPOSITION (LDCOMPUTE,                                &
                                   PRHODREF,  PSSI, PAI, PCJ, PLSFACT,       &
                                   PRST, PCST, PLBDS,                        &
                                   P_RI_CNVI, P_CI_CNVI,                     &
                                   P_TH_DEPS, P_RS_DEPS                      )
!     ##########################################################################
!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the microphysical sources
!!    for slow cold processes :
!!      - conversion of snow to ice
!!      - deposition of vapor on snow
!!      - conversion of ice to snow (Harrington 1995)
!!
!!
!!    AUTHOR
!!    ------
!!      J.-M. Cohard     * Laboratoire d'Aerologie*
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!      S.    Berthet    * Laboratoire d'Aerologie*
!!      B.    Vié        * CNRM *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original             15/03/2018
!!
!  J. Wurtz       03/2022: new snow characteristics
!  B. Vie         03/2022: Add option for 1-moment pristine ice
!  M. Taufour     07/2022: add snow concentration
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAM_LIMA,      ONLY : XRTMIN, XCTMIN, XALPHAS, XNUS, NMOM_I
USE MODD_PARAM_LIMA_COLD, ONLY : XDSCNVI_LIM, XLBDASCNVI_MAX,     &
                                 XC0DEPSI, XC1DEPSI, XR0DEPSI, XR1DEPSI,      &
                                 X1DEPS, X0DEPS, XEX1DEPS, XEX0DEPS,  &
                                 XFVELOS

!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
LOGICAL, DIMENSION(:),INTENT(IN)    :: LDCOMPUTE
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(:),   INTENT(IN)    :: PSSI  ! abs. pressure at time t
REAL, DIMENSION(:),   INTENT(IN)    :: PAI  ! abs. pressure at time t
REAL, DIMENSION(:),   INTENT(IN)    :: PCJ  ! abs. pressure at time t
REAL, DIMENSION(:),   INTENT(IN)    :: PLSFACT  ! abs. pressure at time t
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRST    ! Snow/aggregate m.r. at t 
REAL, DIMENSION(:),   INTENT(IN)    :: PCST    ! Snow/aggregate m.r. at t 
!
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDS    ! Graupel m.r. at t 
!
REAL, DIMENSION(:),   INTENT(OUT)   :: P_RI_CNVI
REAL, DIMENSION(:),   INTENT(OUT)   :: P_CI_CNVI
REAL, DIMENSION(:),   INTENT(OUT)   :: P_TH_DEPS
REAL, DIMENSION(:),   INTENT(OUT)   :: P_RS_DEPS
!
!*       0.2   Declarations of local variables :
!
LOGICAL, DIMENSION(SIZE(PRHODREF)) :: GMICRO ! Computations only where necessary
REAL,    DIMENSION(SIZE(PRHODREF)) :: ZZW, ZZW2, ZZX ! Work array
!
!-------------------------------------------------------------------------------
!
P_RI_CNVI(:) = 0.
P_CI_CNVI(:) = 0.
P_TH_DEPS(:) = 0.
P_RS_DEPS(:) = 0.
!
! Looking for regions where computations are necessary
GMICRO(:) = LDCOMPUTE(:) .AND. PRST(:)>XRTMIN(5)
!
IF (NMOM_I.EQ.1) THEN
   WHERE( GMICRO )
!
! Deposition of water vapor on r_s: RVDEPS
!
      ZZW(:) = 0.0
      WHERE ( PRST(:)>XRTMIN(5) )
         ZZW(:) = PCST(:) * PSSI(:) / PAI(:) * &
              ( X0DEPS*PLBDS(:)**XEX0DEPS +             &
                X1DEPS*PLBDS(:)**XEX1DEPS *PCJ(:) *     &
                     (1+0.5*(XFVELOS/PLBDS(:))**XALPHAS)**(-XNUS+XEX1DEPS/XALPHAS) )
         ZZW(:) =    ZZW(:)*(0.5+SIGN(0.5,ZZW(:))) - ABS(ZZW(:))*(0.5-SIGN(0.5,ZZW(:)))
      END WHERE
      P_RS_DEPS(:) = ZZW(:)
   END WHERE
ELSE
   WHERE( GMICRO )
!
!*       2.1    Conversion of snow to r_i: RSCNVI
!        ----------------------------------------
!
!
      ZZW2(:) = 0.0
      ZZW(:) = 0.0
      WHERE ( PLBDS(:)<XLBDASCNVI_MAX .AND. PRST(:)>XRTMIN(5) .AND. PCST(:)>XCTMIN(5) &
                                      .AND. PSSI(:)<0.0                               )
         ZZW(:) = (PLBDS(:)*XDSCNVI_LIM)**(XALPHAS)
         ZZX(:) = ( -PSSI(:)/PAI(:) ) * PCST(:) * (ZZW(:)**XNUS) * EXP(-ZZW(:))
!
         ZZW(:) = ( XR0DEPSI+XR1DEPSI*PCJ(:) )*ZZX(:)
!
         ZZW2(:)= ( XC0DEPSI+XC1DEPSI*PCJ(:) )*ZZX(:)
      END WHERE
!
      P_RI_CNVI(:) = ZZW(:)
      P_CI_CNVI(:) = ZZW2(:)
!
!
!*       2.2    Deposition of water vapor on r_s: RVDEPS
!        -----------------------------------------------
!
!
      ZZW(:) = 0.0
      WHERE ( PRST(:)>XRTMIN(5) .AND. PCST(:)>XCTMIN(5) )
         ZZW(:) = ( PCST(:)*PSSI(:)/PAI(:) ) *     &
              ( X0DEPS*PLBDS(:)**XEX0DEPS +        &
              ( X1DEPS*PCJ(:)*PLBDS(:)**XEX1DEPS * &
                   (1+0.5*(XFVELOS/PLBDS(:))**XALPHAS)**(-XNUS+XEX1DEPS/XALPHAS)) )
         ZZW(:) =    ZZW(:)*(0.5+SIGN(0.5,ZZW(:))) - ABS(ZZW(:))*(0.5-SIGN(0.5,ZZW(:)))
      END WHERE
!
      P_RS_DEPS(:) = ZZW(:)
! 
   END WHERE
END IF
!
END SUBROUTINE LIMA_SNOW_DEPOSITION
END MODULE MODE_LIMA_SNOW_DEPOSITION
