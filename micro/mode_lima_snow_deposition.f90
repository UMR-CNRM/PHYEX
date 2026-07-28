!MNH_LIC Copyright 2018-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_LIMA_SNOW_DEPOSITION
  IMPLICIT NONE
CONTAINS
!     ##########################################################################
  SUBROUTINE LIMA_SNOW_DEPOSITION (LIMAP, LIMAC, KSIZE, ODCOMPUTE,                         &
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
USE MODD_PARAM_LIMA_COLD, ONLY:PARAM_LIMA_COLD_T
USE MODD_PARAM_LIMA, ONLY:PARAM_LIMA_T
USE YOMHOOK, ONLY:LHOOK, DR_HOOK, JPHOOK

!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER, INTENT(IN) :: KSIZE
LOGICAL, DIMENSION(KSIZE),INTENT(IN)    :: ODCOMPUTE
!
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PSSI  ! abs. pressure at time t
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PAI  ! abs. pressure at time t
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PCJ  ! abs. pressure at time t
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLSFACT  ! abs. pressure at time t
!
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRST    ! Snow/aggregate m.r. at t 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PCST    ! Snow/aggregate m.r. at t 
!
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLBDS    ! Graupel m.r. at t 
!
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_RI_CNVI
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_CI_CNVI
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_TH_DEPS
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_RS_DEPS
!
!*       0.2   Declarations of local variables :
!
LOGICAL, DIMENSION(SIZE(PRHODREF)) :: GMICRO ! Computations only where necessary
TYPE(PARAM_LIMA_COLD_T),INTENT(IN)::LIMAC
TYPE(PARAM_LIMA_T),INTENT(IN)::LIMAP
REAL,    DIMENSION(SIZE(PRHODREF)) :: ZZW, ZZW2, ZZX ! Work array
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('LIMA_SNOW_DEPOSITION', 0, ZHOOK_HANDLE)
P_RI_CNVI(:) = 0.
P_CI_CNVI(:) = 0.
P_TH_DEPS(:) = 0.
P_RS_DEPS(:) = 0.
!
! Looking for regions where computations are necessary
GMICRO(:) = ODCOMPUTE(:) .AND. PRST(:)>LIMAP%XRTMIN(5)
IF (LIMAP%NMOM_I.GE.2) GMICRO(:) = GMICRO(:) .AND. PCST(:)>LIMAP%XCTMIN(5)
!
IF (LIMAP%NMOM_I.EQ.1) THEN
!
! Deposition of water vapor on r_s: RVDEPS
!
   ZZW(:) = 0.0
   WHERE( GMICRO )
      ZZW(:) = PCST(:) * PSSI(:) / PAI(:) * &
           ( LIMAC%X0DEPS*PLBDS(:)**LIMAC%XEX0DEPS +             &
           LIMAC%X1DEPS*PLBDS(:)**LIMAC%XEX1DEPS *PCJ(:) *     &
           (1+0.5*(LIMAC%XFVELOS/PLBDS(:))**LIMAP%XALPHAS)**(-LIMAP%XNUS+LIMAC%XEX1DEPS/LIMAP%XALPHAS) )
      ZZW(:) =    ZZW(:)*(0.5+SIGN(0.5,ZZW(:))) - ABS(ZZW(:))*(0.5-SIGN(0.5,ZZW(:)))
      P_RS_DEPS(:) = ZZW(:)
   END WHERE
ELSE
!
!*       2.1    Conversion of snow to r_i: RSCNVI
!        ----------------------------------------
!
   ZZW2(:) = 0.0
   ZZW(:) = 0.0
   WHERE ( GMICRO .AND. PLBDS(:)<LIMAC%XLBDASCNVI_MAX .AND. PSSI(:)<0.0 )
      ZZW(:) = (PLBDS(:)*LIMAC%XDSCNVI_LIM)**(LIMAP%XALPHAS)
      ZZX(:) = ( -PSSI(:)/PAI(:) ) * PCST(:) * (ZZW(:)**LIMAP%XNUS) * EXP(-ZZW(:))
!
      ZZW(:) = ( LIMAC%XR0DEPSI+LIMAC%XR1DEPSI*PCJ(:) )*ZZX(:)
!
      ZZW2(:)= ( LIMAC%XC0DEPSI+LIMAC%XC1DEPSI*PCJ(:) )*ZZX(:)
      P_RI_CNVI(:) = ZZW(:)
      P_CI_CNVI(:) = ZZW2(:)
   END WHERE
!
!*       2.2    Deposition of water vapor on r_s: RVDEPS
!        -----------------------------------------------
!
      ZZW(:) = 0.0
      WHERE ( GMICRO )
         ZZW(:) = ( PCST(:)*PSSI(:)/PAI(:) ) *     &
              ( LIMAC%X0DEPS*PLBDS(:)**LIMAC%XEX0DEPS +        &
              ( LIMAC%X1DEPS*PCJ(:)*PLBDS(:)**LIMAC%XEX1DEPS * &
                   (1+0.5*(LIMAC%XFVELOS/PLBDS(:))**LIMAP%XALPHAS)**(-LIMAP%XNUS+LIMAC%XEX1DEPS/LIMAP%XALPHAS)) )
         ZZW(:) =    ZZW(:)*(0.5+SIGN(0.5,ZZW(:))) - ABS(ZZW(:))*(0.5-SIGN(0.5,ZZW(:)))
         P_RS_DEPS(:) = ZZW(:)
      END WHERE
! 
END IF
!
IF (LHOOK) CALL DR_HOOK('LIMA_SNOW_DEPOSITION', 1, ZHOOK_HANDLE)
END SUBROUTINE LIMA_SNOW_DEPOSITION
END MODULE MODE_LIMA_SNOW_DEPOSITION
