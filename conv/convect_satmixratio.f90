!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 conv 2006/05/18 13:07:25
!-----------------------------------------------------------------
!     #################
      MODULE MODI_CONVECT_SATMIXRATIO
!     #################
!
INTERFACE
!
       SUBROUTINE CONVECT_SATMIXRATIO( KLON,                          &
                                      PPRES, PT, PEW, PLV, PLS, PCPH )
!
INTEGER,                INTENT(IN) :: KLON    ! horizontal loop index
REAL, DIMENSION(KLON),  INTENT(IN) :: PPRES   ! pressure
REAL, DIMENSION(KLON),  INTENT(IN) :: PT      ! temperature   
!
REAL, DIMENSION(KLON),  INTENT(OUT):: PEW     ! vapor saturation mixing ratio
REAL, DIMENSION(KLON),  INTENT(OUT):: PLV     ! latent heat L_v    
REAL, DIMENSION(KLON),  INTENT(OUT):: PLS     ! latent heat L_s  
REAL, DIMENSION(KLON),  INTENT(OUT):: PCPH    ! specific heat C_ph  
!
END SUBROUTINE CONVECT_SATMIXRATIO
!
END INTERFACE
!
END MODULE MODI_CONVECT_SATMIXRATIO
!     ######spl
      SUBROUTINE CONVECT_SATMIXRATIO( KLON,                          &
                                      PPRES, PT, PEW, PLV, PLS, PCPH )      
!     ################################################################
!
!!**** Compute vapor saturation mixing ratio over liquid water
!!
!!
!!    PDRPOSE
!!    -------
!!     The purpose of this routine is to determine saturation mixing ratio
!!     and to return values for L_v L_s and C_ph
!!
!!
!!**  METHOD
!!    ------
!!     
!!
!!    EXTERNAL
!!    --------
!!     None
!!     
!!
!!    IMPLICIT ARGUMENTS    
!!    ------------------
!!      Module MODD_CST
!!          XALPW, XBETAW, XGAMW ! constants for water saturation pressure
!!          XRD, XRV             ! gaz  constants for dry air and water vapor
!!          XCPD, XCPV           ! specific heat for dry air and water vapor
!!          XCL, XCI             ! specific heat for liquid water and ice
!!          XTT                  ! triple point temperature
!!          XLVTT, XLSTT         ! vaporization, sublimation heat constant
!!
!!
!!    REFERENCE
!!    ---------
!!
!!      Book1,2 of documentation ( routine CONVECT_SATMIXRATIO)
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    07/11/95 
!!   Last modified  04/10/97
!------------------------- ------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
!
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
INTEGER,                INTENT(IN) :: KLON    ! horizontal loop index
REAL, DIMENSION(KLON),  INTENT(IN) :: PPRES   ! pressure
REAL, DIMENSION(KLON),  INTENT(IN) :: PT      ! temperature   
!
REAL, DIMENSION(KLON),  INTENT(OUT):: PEW     ! vapor saturation mixing ratio
REAL, DIMENSION(KLON),  INTENT(OUT):: PLV     ! latent heat L_v    
REAL, DIMENSION(KLON),  INTENT(OUT):: PLS     ! latent heat L_s  
REAL, DIMENSION(KLON),  INTENT(OUT):: PCPH    ! specific heat C_ph   
!
!*       0.2   Declarations of local variables :
!
REAL, DIMENSION(KLON)              :: ZT      ! temperature   
REAL    :: ZEPS           ! R_d / R_v
!
!
!-------------------------------------------------------------------------------
!
    ZEPS      = XRD / XRV
!
    ZT(:)     = MIN( 400., MAX( PT(:), 10. ) ) ! overflow bound
    PEW(:)    = EXP( XALPW - XBETAW / ZT(:) - XGAMW * ALOG( ZT(:) ) )
    PEW(:)    = ZEPS * PEW(:) / ( PPRES(:) - PEW(:) )
!
    PLV(:)    = XLVTT + ( XCPV - XCL ) * ( ZT(:) - XTT ) ! compute L_v
    PLS(:)    = XLSTT + ( XCPV - XCI ) * ( ZT(:) - XTT ) ! compute L_i
!    
    PCPH(:)   = XCPD + XCPV * PEW(:)                     ! compute C_ph 
!
END SUBROUTINE CONVECT_SATMIXRATIO
