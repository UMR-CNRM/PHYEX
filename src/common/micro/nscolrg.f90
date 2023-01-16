!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 init 2006/11/23 10:43:02
!-----------------------------------------------------------------
!     ###################
      MODULE MODI_NSCOLRG
!     ###################
!
INTERFACE
!
      SUBROUTINE NSCOLRG( KND, PALPHAS, PZNUS, PALPHAR, PNUR,                &
                         PESR, PFALLS, PEXFALLS, PFALLEXPS, PFALLR, PEXFALLR,&
                         PLBDASMAX, PLBDARMAX, PLBDASMIN, PLBDARMIN,         &
                         PDINFTY, PNSCOLRG,PAG, PBS, PAS                     )
!
INTEGER, INTENT(IN) :: KND    ! Number of discrete size intervals in DS and DR  
!
REAL, INTENT(IN) :: PALPHAS   ! First shape parameter of the aggregates 
                              ! size distribution (generalized gamma law)
REAL, INTENT(IN) :: PZNUS     ! Second shape parameter of the aggregates
                              ! size distribution (generalized gamma law)
REAL, INTENT(IN) :: PALPHAR   ! First shape parameter of the rain  
                              ! size distribution (generalized gamma law)
REAL, INTENT(IN) :: PNUR      ! Second shape parameter of the rain 
                              ! size distribution (generalized gamma law)
REAL, INTENT(IN) :: PESR      ! Efficiency of the aggregates collecting rain 
REAL, INTENT(IN) :: PFALLS    ! Fall speed constant of the aggregates
REAL, INTENT(IN) :: PEXFALLS  ! Fall speed exponent of the aggregates
REAL, INTENT(IN) :: PFALLEXPS  ! Fall speed exponential constant of the aggregates
REAL, INTENT(IN) :: PFALLR    ! Fall speed constant of rain 
REAL, INTENT(IN) :: PEXFALLR  ! Fall speed exponent of rain 
REAL, INTENT(IN) :: PLBDASMAX ! Maximun slope of size distribution of the aggregates
REAL, INTENT(IN) :: PLBDARMAX ! Maximun slope of size distribution of rain 
REAL, INTENT(IN) :: PLBDASMIN ! Minimun slope of size distribution of the aggregates
REAL, INTENT(IN) :: PLBDARMIN ! Minimun slope of size distribution of rain 
REAL, INTENT(IN) :: PDINFTY   ! Factor to define the largest diameter up to
                              ! which the diameter integration is performed
REAL, INTENT(IN) :: PAG, PBS, PAS
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PNSCOLRG! Scaled fall speed difference in
                                               ! the mass collection kernel as a
                                               ! function of LAMBDAX and LAMBDAZ
!
      END SUBROUTINE NSCOLRG
!
END INTERFACE
!
      END MODULE MODI_NSCOLRG
!     ########################################################################
      SUBROUTINE NSCOLRG( KND, PALPHAS, PZNUS, PALPHAR, PNUR,                &
                         PESR, PFALLS, PEXFALLS, PFALLEXPS, PFALLR, PEXFALLR,&
                         PLBDASMAX, PLBDARMAX, PLBDASMIN, PLBDARMIN,         &
                         PDINFTY, PNSCOLRG,PAG, PBS, PAS                     )
!     ########################################################################
!
!
!
!!****  * -  Build up a look-up table containing the scaled fall speed
!!           difference between size distributed particles of the aggregates and Z
!!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to integrate numerically the scaled fall
!!      speed difference between aggregates and rain  for use in collection
!!      kernels. A first integral of the form
!!
!!       infty  Dz_max
!!           / /
!!           |{|                                                 }
!!           |{| E_xz (Dx+Dz)^2 |cxDx^dx-czDz^dz| n(Dz) dDz} n(Dx) dDx
!!           |{|                                                 }
!!           / /
!!          0 Dz_min
!!
!!      is evaluated and normalised by a second integral of the form
!!
!!              infty
!!             / /
!!             |{|                          }
!!             |{| (Dx+Dz)^2 n(Dz) dDz} n(Dx) dDx
!!             |{|                          }
!!             / /
!!              0
!!
!!      The result is stored in a two-dimensional array.
!! 
!!**  METHOD
!!    ------
!!      The free parameters of the size distribution function of the aggregates
!!      and Z (slope parameter LAMBDA) are discretized with a geometrical rate 
!!      in a specific range
!!            LAMBDA = exp( (Log(LAMBDA_max) - Log(LAMBDA_min))/N_interval )
!!      The two above integrals are performed using the trapezoidal scheme.
!!
!!    EXTERNAL
!!    --------
!!      MODI_GENERAL_GAMMA: Generalized gamma distribution law 
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      MODD_CST           : XPI,XRHOLW
!!      MODD_RAIN_ICE_DESCR: XAS,XAS,XBS
!!
!!    REFERENCE
!!    ---------
!!      B.S. Ferrier , 1994 : A Double-Moment Multiple-Phase Four-Class
!!                            Bulk Ice Scheme,JAS,51,249-280.
!!
!!    AUTHOR
!!    ------
!!      J.-P. Pinty     * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS 
!!    -------------
!!      Original    8/11/95
!!      M. Taufour    03/2022  Adapted from rscolrg for concentration
!!      J. Wurtz      03/2022  New snow characteristics        
!!
!-------------------------------------------------------------------------------
!
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODI_GENERAL_GAMMA
!
USE MODD_CST
USE MODD_RAIN_ICE_DESCR
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments 
!              ------------------------------- 
!
!
INTEGER, INTENT(IN) :: KND    ! Number of discrete size intervals in DS and DR  
!
REAL, INTENT(IN) :: PALPHAS   ! First shape parameter of the aggregates 
                              ! size distribution (generalized gamma law)
REAL, INTENT(IN) :: PZNUS     ! Second shape parameter of the aggregates
                              ! size distribution (generalized gamma law)
REAL, INTENT(IN) :: PALPHAR   ! First shape parameter of the rain  
                              ! size distribution (generalized gamma law)
REAL, INTENT(IN) :: PNUR      ! Second shape parameter of the rain 
                              ! size distribution (generalized gamma law)
REAL, INTENT(IN) :: PESR      ! Efficiency of the aggregates collecting rain 
REAL, INTENT(IN) :: PFALLS    ! Fall speed constant of the aggregates
REAL, INTENT(IN) :: PEXFALLS  ! Fall speed exponent of the aggregates
REAL, INTENT(IN) :: PFALLEXPS ! Fall speed exponential constant of the aggregates
REAL, INTENT(IN) :: PFALLR    ! Fall speed constant of rain 
REAL, INTENT(IN) :: PEXFALLR  ! Fall speed exponent of rain 
REAL, INTENT(IN) :: PLBDASMAX ! Maximun slope of size distribution of the aggregates
REAL, INTENT(IN) :: PLBDARMAX ! Maximun slope of size distribution of rain 
REAL, INTENT(IN) :: PLBDASMIN ! Minimun slope of size distribution of the aggregates
REAL, INTENT(IN) :: PLBDARMIN ! Minimun slope of size distribution of rain 
REAL, INTENT(IN) :: PDINFTY   ! Factor to define the largest diameter up to
                              ! which the diameter integration is performed
REAL, INTENT(IN) :: PAG, PBS, PAS
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PNSCOLRG! Scaled fall speed difference in
                                               ! the mass collection kernel as a
                                               ! function of LAMBDAX and LAMBDAZ
!
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
!
INTEGER :: JLBDAS  ! Slope index of the size distribution of the aggregates
INTEGER :: JLBDAR  ! Slope index of the size distribution of rain 
INTEGER :: JDS     ! Diameter index of a particle of the aggregates
INTEGER :: JDR     ! Diameter index of a particle of rain 
!
INTEGER :: INR     ! Number of diameter step for the partial integration
!
REAL :: ZLBDAS  ! Current slope parameter LAMBDA of the aggregates
REAL :: ZLBDAR  ! Current slope parameter LAMBDA of rain 
REAL :: ZDLBDAS ! Growth rate of the slope parameter LAMBDA of the aggregates
REAL :: ZDLBDAR ! Growth rate of the slope parameter LAMBDA of rain 
REAL :: ZDDS    ! Integration step of the diameter of the aggregates
REAL :: ZDDSCALR! Integration step of the diameter of rain  (scaling integral)
REAL :: ZDDCOLLR! Integration step of the diameter of rain  (fallspe integral)
REAL :: ZDS     ! Current diameter of the particle aggregates
REAL :: ZDR     ! Current diameter of the raindrops 
REAL :: ZDRMIN  ! Minimal diameter of the raindrops where the integration starts
REAL :: ZDRMAX  ! Maximal diameter of the raindrops where the integration ends
REAL :: ZCOLLR  ! Single integral of the mass weighted fall speed difference 
                ! over the spectrum of rain 
REAL :: ZCOLLDRMIN ! Minimum ending point for the partial integral
REAL :: ZCOLLSR ! Double integral of the mass weighted fall speed difference
                ! over the spectra of the aggregates and rain 
REAL :: ZSCALR  ! Single integral of the scaling factor over 
                ! the spectrum of rain 
REAL :: ZSCALSR ! Double integral of the scaling factor over
                ! the spectra of the aggregates and rain 
REAL :: ZFUNC   ! Ancillary function
REAL :: ZCST1
!
!
!-------------------------------------------------------------------------------
!
!
!*       1       COMPUTE THE SCALED VELOCITY DIFFERENCE IN THE MASS
!*                               COLLECTION KERNEL,
!                -------------------------------------------------
!
!
!*       1.0     Initialization
!
PNSCOLRG(:,:) = 0.0
ZCST1  = (3.0/XPI)/XRHOLW
!
!*       1.1     Compute the growth rate of the slope factors LAMBDA
!
ZDLBDAR = EXP( LOG(PLBDARMAX/PLBDARMIN)/REAL(SIZE(PNSCOLRG(:,:),1)-1) )
ZDLBDAS = EXP( LOG(PLBDASMAX/PLBDASMIN)/REAL(SIZE(PNSCOLRG(:,:),2)-1) )
!
!*       1.2     Scan the slope factors LAMBDAX and LAMBDAZ
!
DO JLBDAR = 1,SIZE(PNSCOLRG(:,:),1)
  ZLBDAR = PLBDARMIN * ZDLBDAR ** (JLBDAR-1) 
  ZDRMAX = PDINFTY / ZLBDAR
!
!*       1.3     Compute the diameter steps
!
  ZDDSCALR = PDINFTY / (REAL(KND) * ZLBDAR)
  DO JLBDAS = 1,SIZE(PNSCOLRG(:,:),2)
    ZLBDAS = PLBDASMIN * ZDLBDAS ** (JLBDAS-1)
!
!*       1.4     Initialize the collection integrals
!
    ZSCALSR = 0.0
    ZCOLLSR = 0.0
!
!*       1.5     Compute the diameter steps
!
    ZDDS     = PDINFTY / (REAL(KND) * ZLBDAS)
!
!*       1.6     Scan over the diameters DS and DR
!
    DO JDS = 1,KND-1
      ZDS = ZDDS * REAL(JDS)
      ZSCALR = 0.0
      ZCOLLR = 0.0
      DO JDR = 1,KND-1
        ZDR = ZDDSCALR * REAL(JDR)
!
!*       1.7     Compute the normalization factor by integration over the
!                dimensional spectrum of rain   
!
        ZSCALR = ZSCALR + (ZDS+ZDR)**2 * GENERAL_GAMMA(PALPHAR,PNUR,ZLBDAR,ZDR)
      END DO
!
!*       1.8     Compute the scaled fall speed difference by partial
!                integration over the dimensional spectrum of rain   
!
      ZFUNC = PAG - PAS*ZDS**(PBS-3.0) ! approximate limit is Ds=240 microns
      IF( ZFUNC>0.0 ) THEN
        ZDRMIN = ZDS*( ZCST1*ZFUNC )**0.3333333
        ELSE
        ZDRMIN = 0.0
      END IF
      IF( ZDS>1.0E-4 ) THEN            ! allow computation if Ds>100 microns 
            ! corresponding to a maximal density of the aggregates of XRHOLW
        IF( (ZDRMAX-ZDRMIN) >= 0.5*ZDDSCALR ) THEN
          INR = CEILING( (ZDRMAX-ZDRMIN)/ZDDSCALR )
          ZDDCOLLR = (ZDRMAX-ZDRMIN) / REAL(INR)
          DO JDR = 1,INR-1
            ZDR = ZDDCOLLR * REAL(JDR) + ZDRMIN
            ZCOLLR = ZCOLLR + (ZDS+ZDR)**2                                     &
                       * GENERAL_GAMMA(PALPHAR,PNUR,ZLBDAR,ZDR)                &
                         * PESR * ABS(PFALLS*ZDS**PEXFALLS*EXP(-(ZDS*PFALLEXPS)**PALPHAS)-PFALLR*ZDR**PEXFALLR)
          END DO
          IF( ZDRMIN>0.0 ) THEN
            ZCOLLDRMIN = (ZDS+ZDRMIN)**2                                       &
                      * GENERAL_GAMMA(PALPHAR,PNUR,ZLBDAR,ZDRMIN)              &
                      * PESR * ABS(PFALLS*ZDS**PEXFALLS*EXP(-(ZDS*PFALLEXPS)**PALPHAS)-PFALLR*ZDRMIN**PEXFALLR)
            ELSE
            ZCOLLDRMIN = 0.0
          END IF 
          ZCOLLR = (ZCOLLR + 0.5*ZCOLLDRMIN)*(ZDDCOLLR/ZDDSCALR)
!
!*       1.9     Compute the normalization factor by integration over the
!                dimensional spectrum of the aggregates
!
          ZFUNC   = GENERAL_GAMMA(PALPHAS,PZNUS,ZLBDAS,ZDS)  ! MTaufour : !*(ZDS**PEXMASSS)
          ZSCALSR = ZSCALSR + ZSCALR * ZFUNC
!
!*       1.10    Compute the scaled fall speed difference by integration over
!                the dimensional spectrum of the aggregates
!
          ZCOLLSR = ZCOLLSR + ZCOLLR * ZFUNC
!
! Otherwise ZDRMIN>ZDRMAX so PRRCOLSS(JLBDAS,JLBDAR) = 0.0        !
!
        END IF
!
! Otherwise ZDRMAX = 0.0 so the density of the graupel cannot be reached
!                    and so PRRCOLSS(JLBDAS,JLBDAR) = 0.0        !
!
      END IF
    END DO
!
!*       1.10    Scale the fall speed difference
!
    IF( ZSCALSR>0.0 ) PNSCOLRG(JLBDAR,JLBDAS) = ZCOLLSR / ZSCALSR
  END DO
END DO
!
END SUBROUTINE NSCOLRG
