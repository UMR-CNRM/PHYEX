!MNH_LIC Copyright 2018-2024 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-------------------------------------------------------------------------------
MODULE MODE_LIMA_CONVERSION_MELTING_SNOW
  IMPLICIT NONE
CONTAINS
!     ##############################################################################
  SUBROUTINE LIMA_CONVERSION_MELTING_SNOW (CST, LIMAP, LIMAC, LIMAM, KSIZE, ODCOMPUTE,&
                                           PRHODREF, PPRES, PT, PKA, PDV, PCJ, &
                                           PRVT, PRST, PCST, PLBDS,            &
                                           P_RS_CMEL, P_CS_CMEL                )
!     ##############################################################################
!
!!    PURPOSE
!!    -------
!!      Compute the conversion-melting of snow into graupel
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
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAM_LIMA_MIXED, ONLY:PARAM_LIMA_MIXED_T
USE MODD_PARAM_LIMA_COLD, ONLY:PARAM_LIMA_COLD_T
USE MODD_PARAM_LIMA, ONLY:PARAM_LIMA_T
USE MODD_CST, ONLY:CST_T
USE YOMHOOK, ONLY:LHOOK, DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(CST_T),INTENT(IN)::CST
TYPE(PARAM_LIMA_MIXED_T),INTENT(IN)::LIMAM
TYPE(PARAM_LIMA_COLD_T),INTENT(IN)::LIMAC
TYPE(PARAM_LIMA_T),INTENT(IN)::LIMAP
INTEGER, INTENT(IN) :: KSIZE
LOGICAL, DIMENSION(KSIZE),INTENT(IN)    :: ODCOMPUTE
!
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRHODREF ! Reference Exner function
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PPRES    !
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PT       !
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PKA      !
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PDV      !
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PCJ      !
!
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRVT    ! 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRST    ! Snow mr at t
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PCST    ! Snow C. at t
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLBDS   ! 
!
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_RS_CMEL
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_CS_CMEL
!
!*       0.2   Declarations of local variables :
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
REAL, DIMENSION(SIZE(PRST)) :: ZW ! work arrays
!
!-------------------------------------------------------------------------------
!
!
!*       1.     Conversion-melting of snow
!               --------------------------
!
!
IF (LHOOK) CALL DR_HOOK('LIMA_CONVERSION_MELTING_SNOW', 0, ZHOOK_HANDLE)
P_RS_CMEL(:)=0.
P_CS_CMEL(:)=0.
!
ZW(:) = 0.0
WHERE( PRST(:)>LIMAP%XRTMIN(5) .AND. PCST(:)>LIMAP%XCTMIN(5) .AND. PT(:)>CST%XTT .AND. ODCOMPUTE(:) )
   ZW(:) = PRVT(:)*PPRES(:)/((CST%XMV/CST%XMD)+PRVT(:)) ! Vapor pressure
   ZW(:) = PKA(:)*(CST%XTT-PT(:)) +                                 &
              ( PDV(:)*(CST%XLVTT + ( CST%XCPV - CST%XCL ) * ( PT(:) - CST%XTT )) &
                          *(CST%XESTT-ZW(:))/(CST%XRV*PT(:))             )
!
! compute RSMLT
!
   ZW(:)  = LIMAM%XFSCVMG*MAX( 0.0,( -ZW(:) * PCST(:) *                        &
                               ( LIMAC%X0DEPS*PLBDS(:)**LIMAC%XEX0DEPS +             &
                                 LIMAC%X1DEPS*PCJ(:)*PLBDS(:)**LIMAC%XEX1DEPS *      &
                                   (1+0.5*(LIMAC%XFVELOS/PLBDS(:))**LIMAP%XALPHAS)**(-LIMAP%XNUS+LIMAC%XEX1DEPS/LIMAP%XALPHAS)) ))
! On ne tient pas compte de la collection de pluie et gouttelettes par la neige si T>0 !!!! 
! Note that no heat is exchanged because the graupeln produced are still icy!!!
   P_RS_CMEL(:) = - ZW(:)
   P_CS_CMEL(:) = - ZW(:) * PCST(:) / PRST(:)
!
END WHERE
!
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('LIMA_CONVERSION_MELTING_SNOW', 1, ZHOOK_HANDLE)
END SUBROUTINE LIMA_CONVERSION_MELTING_SNOW
END MODULE MODE_LIMA_CONVERSION_MELTING_SNOW
