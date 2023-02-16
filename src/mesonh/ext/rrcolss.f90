!MNH_LIC Copyright 1995-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ###################
      MODULE MODI_RRCOLSS
!     ###################
!
INTERFACE
!
      SUBROUTINE RRCOLSS( KND, PALPHAS, PNUS, PALPHAR, PNUR,                 &
                         PESR, PEXMASSR, PFALLS, PEXFALLS, PFALLEXPS, PFALLR, PEXFALLR, &
                         PLBDASMAX, PLBDARMAX, PLBDASMIN, PLBDARMIN,         &
                         PDINFTY, PRRCOLSS, PAG, PBS, PAS                    )
!
INTEGER, INTENT(IN) :: KND    ! Number of discrete size intervals in DS and DR  
!
REAL, INTENT(IN) :: PALPHAS   ! First shape parameter of the aggregates 
                              ! size distribution (generalized gamma law)
REAL, INTENT(IN) :: PNUS      ! Second shape parameter of the aggregates
                              ! size distribution (generalized gamma law)
REAL, INTENT(IN) :: PALPHAR   ! First shape parameter of the rain  
                              ! size distribution (generalized gamma law)
REAL, INTENT(IN) :: PNUR      ! Second shape parameter of the rain 
                              ! size distribution (generalized gamma law)
REAL, INTENT(IN) :: PESR      ! Efficiency of aggregates collecting rain 
REAL, INTENT(IN) :: PEXMASSR  ! Mass exponent of rain 
REAL, INTENT(IN) :: PFALLS    ! Fall speed constant of aggregates
REAL, INTENT(IN) :: PEXFALLS  ! Fall speed exponent of aggregates
REAL, INTENT(IN) :: PFALLEXPS ! Fall speed exponential of aggregates (Thompson 2008)
REAL, INTENT(IN) :: PFALLR    ! Fall speed constant of rain 
REAL, INTENT(IN) :: PEXFALLR  ! Fall speed exponent of rain 
REAL, INTENT(IN) :: PLBDASMAX ! Maximun slope of size distribution of aggregates
REAL, INTENT(IN) :: PLBDARMAX ! Maximun slope of size distribution of rain 
REAL, INTENT(IN) :: PLBDASMIN ! Minimun slope of size distribution of aggregates
REAL, INTENT(IN) :: PLBDARMIN ! Minimun slope of size distribution of rain 
REAL, INTENT(IN) :: PDINFTY   ! Factor to define the largest diameter up to
                              ! which the diameter integration is performed
REAL, INTENT(IN) :: PAG, PBS, PAS
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PRRCOLSS! Scaled fall speed difference in
                                               ! the mass collection kernel as a
                                               ! function of LAMBDAX and LAMBDAZ
!
      END SUBROUTINE RRCOLSS
!
END INTERFACE
!
      END MODULE MODI_RRCOLSS
!     ########################################################################
      SUBROUTINE RRCOLSS( KND, PALPHAS, PNUS, PALPHAR, PNUR,                 &
                         PESR, PEXMASSR, PFALLS, PEXFALLS, PFALLEXPS, PFALLR, PEXFALLR, &
                         PLBDASMAX, PLBDARMAX, PLBDASMIN, PLBDARMIN,         &
                         PDINFTY, PRRCOLSS, PAG, PBS, PAS                    )
!     ########################################################################
!
!
!
!!****  * -  Build up a look-up table containing the scaled fall speed
!!           difference between size distributed particles of aggregates and Z
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
!!           |{| E_xz (Dx+Dz)^2 |cxDx^dx-czDz^dz| Dz^bz n(Dz) dDz} n(Dx) dDx
!!           |{|                                                 }
!!           / /
!!          0 Dz_min
!!
!!      is evaluated and normalised by a second integral of the form
!!
!!              infty
!!             / /
!!             |{|                          }
!!             |{| (Dx+Dz)^2 Dz^bz n(Dz) dDz} n(Dx) dDx
!!             |{|                          }
!!             / /
!!              0
!!
!!      The result is stored in a two-dimensional array.
!! 
!!**  METHOD
!!    ------
!!      The free parameters of the size distribution function of aggregates and Z
!!      (slope parameter LAMBDA) are discretized with a geometrical rate in a
!!      specific range
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
!!
!  P. Wautelet 26/04/2019: replace non-standard FLOAT function by REAL function
!  J. Wurtz       03/2022: new snow characteristics
!
!-------------------------------------------------------------------------------
!
!
!*       0.    DECLARATIONS
!              ------------
!
!
USE MODI_GENERAL_GAMMA
!
USE MODD_CST
USE MODD_RAIN_ICE_DESCR
!
IMPLICIT NONE
!
!
!*       0.1   Declarations of dummy arguments 
!              ------------------------------- 
!
!
INTEGER, INTENT(IN) :: KND    ! Number of discrete size intervals in DS and DR  
!
REAL, INTENT(IN) :: PALPHAS   ! First shape parameter of the aggregates 
                              ! size distribution (generalized gamma law)
REAL, INTENT(IN) :: PNUS      ! Second shape parameter of the aggregates
                              ! size distribution (generalized gamma law)
REAL, INTENT(IN) :: PALPHAR   ! First shape parameter of the rain  
                              ! size distribution (generalized gamma law)
REAL, INTENT(IN) :: PNUR      ! Second shape parameter of the rain 
                              ! size distribution (generalized gamma law)
REAL, INTENT(IN) :: PESR      ! Efficiency of aggregates collecting rain 
REAL, INTENT(IN) :: PEXMASSR  ! Mass exponent of rain 
REAL, INTENT(IN) :: PFALLS    ! Fall speed constant of aggregates
REAL, INTENT(IN) :: PEXFALLS  ! Fall speed exponent of aggregates
REAL, INTENT(IN) :: PFALLEXPS ! Fall speed exponential of aggregates  (Thompson 2008)
REAL, INTENT(IN) :: PFALLR    ! Fall speed constant of rain 
REAL, INTENT(IN) :: PEXFALLR  ! Fall speed exponent of rain 
REAL, INTENT(IN) :: PLBDASMAX ! Maximun slope of size distribution of aggregates
REAL, INTENT(IN) :: PLBDARMAX ! Maximun slope of size distribution of rain 
REAL, INTENT(IN) :: PLBDASMIN ! Minimun slope of size distribution of aggregates
REAL, INTENT(IN) :: PLBDARMIN ! Minimun slope of size distribution of rain 
REAL, INTENT(IN) :: PDINFTY   ! Factor to define the largest diameter up to
                              ! which the diameter integration is performed
REAL, INTENT(IN) :: PAG, PBS, PAS
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PRRCOLSS! Scaled fall speed difference in
                                               ! the mass collection kernel as a
                                               ! function of LAMBDAX and LAMBDAZ
!
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
!
INTEGER :: JLBDAS  ! Slope index of the size distribution of aggregates
INTEGER :: JLBDAR  ! Slope index of the size distribution of rain 
INTEGER :: JDS     ! Diameter index of a particle of aggregates
INTEGER :: JDR     ! Diameter index of a particle of rain 
!
INTEGER :: INR     ! Number of diameter step for the partial integration
!
!
REAL :: ZLBDAS  ! Current slope parameter LAMBDA of aggregates
REAL :: ZLBDAR  ! Current slope parameter LAMBDA of rain 
REAL :: ZDLBDAS ! Growth rate of the slope parameter LAMBDA of aggregates
REAL :: ZDLBDAR ! Growth rate of the slope parameter LAMBDA of rain 
REAL :: ZDDS    ! Integration step of the diameter of aggregates
REAL :: ZDDSCALR! Integration step of the diameter of rain  (scaling integral)
REAL :: ZDDCOLLR! Integration step of the diameter of rain  (fallspe integral)
REAL :: ZDS     ! Current diameter of the particle aggregates
REAL :: ZDR     ! Current diameter of the rain 
REAL :: ZDRMAX  ! Maximal diameter of the raindrops where the integration ends 
REAL :: ZCOLLR  ! Single integral of the mass weighted fall speed difference 
                ! over the spectrum of rain 
REAL :: ZCOLLDRMAX ! Maximum ending point for the partial integral
REAL :: ZCOLLSR ! Double integral of the mass weighted fall speed difference
                ! over the spectra of aggregates and rain 
REAL :: ZSCALR  ! Single integral of the scaling factor over 
                ! the spectrum of rain 
REAL :: ZSCALSR ! Double integral of the scaling factor over
                ! the spectra of aggregates and rain 
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
!
!*       1.0     Initialization
!
PRRCOLSS(:,:) = 0.0
ZCST1 = (3.0/XPI)/XRHOLW
!
!*       1.1     Compute the growth rate of the slope factors LAMBDA
!
ZDLBDAS = EXP( LOG(PLBDASMAX/PLBDASMIN)/REAL(SIZE(PRRCOLSS(:,:),1)-1) )
ZDLBDAR = EXP( LOG(PLBDARMAX/PLBDARMIN)/REAL(SIZE(PRRCOLSS(:,:),2)-1) )
!
!*       1.2     Scan the slope factors LAMBDAX and LAMBDAZ
!
DO JLBDAS = 1,SIZE(PRRCOLSS(:,:),1)
  ZLBDAS = PLBDASMIN * ZDLBDAS ** (JLBDAS-1) 
!
!*       1.3     Compute the diameter steps
!
  ZDDS   = PDINFTY / (REAL(KND) * ZLBDAS)
  DO JLBDAR = 1,SIZE(PRRCOLSS(:,:),2)
    ZLBDAR = PLBDARMIN * ZDLBDAR ** (JLBDAR-1)
!
!*       1.4     Initialize the collection integrals
!
    ZSCALSR = 0.0
    ZCOLLSR = 0.0
!
!*       1.5     Compute the diameter steps
!
    ZDDSCALR = PDINFTY / (REAL(KND) * ZLBDAR)
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
        ZSCALR = ZSCALR + (ZDS+ZDR)**2 * ZDR**PEXMASSR                         &
                                      * GENERAL_GAMMA(PALPHAR,PNUR,ZLBDAR,ZDR)
      END DO
!
!*       1.8     Compute the scaled fall speed difference by partial
!                integration over the dimensional spectrum of rain   
!
      ZFUNC = PAG - PAS*ZDS**(PBS-3.0) ! approximate limit is Ds=240 microns
      IF( ZFUNC>0.0 ) THEN
        ZDRMAX = ZDS*( ZCST1*ZFUNC )**0.3333333
        ELSE
        ZDRMAX = PDINFTY / ZLBDAR
      END IF
      IF( ZDS>1.0E-4 ) THEN            ! allow computation if Ds>100 microns
            ! corresponding to a maximal density of the aggregates of XRHOLW
        IF( ZDRMAX >= 0.5*ZDDSCALR ) THEN
          INR = CEILING( ZDRMAX/ZDDSCALR )
          ZDDCOLLR = ZDRMAX / REAL(INR)
          IF (INR>=KND ) THEN
            INR = KND
            ZDDCOLLR = ZDDSCALR
          END IF
          DO JDR = 1,INR-1
            ZDR = ZDDCOLLR * REAL(JDR)
            ZCOLLR = ZCOLLR + (ZDS+ZDR)**2 * ZDR**PEXMASSR                     &
                       * PESR * ABS(PFALLS*ZDS**PEXFALLS * EXP(-(PFALLEXPS*ZDS)**PALPHAS)-PFALLR*ZDR**PEXFALLR) &
                                      * GENERAL_GAMMA(PALPHAR,PNUR,ZLBDAR,ZDR)
          END DO
          ZCOLLDRMAX = (ZDS+ZDRMAX)**2 * ZDRMAX**PEXMASSR                      &
                    * PESR * ABS(PFALLS*ZDS**PEXFALLS* EXP(-(PFALLEXPS*ZDS)**PALPHAS)-PFALLR*ZDRMAX**PEXFALLR) &
                                   * GENERAL_GAMMA(PALPHAR,PNUR,ZLBDAR,ZDRMAX)
          ZCOLLR = (ZCOLLR + 0.5*ZCOLLDRMAX)*(ZDDCOLLR/ZDDSCALR)
!
!*       1.9     Compute the normalization factor by integration over the
!                dimensional spectrum of aggregates
!
          ZFUNC   = GENERAL_GAMMA(PALPHAS,PNUS,ZLBDAS,ZDS)
          ZSCALSR = ZSCALSR + ZSCALR * ZFUNC
!
!*       1.10    Compute the scaled fall speed difference by integration over
!                the dimensional spectrum of aggregates
!
          ZCOLLSR = ZCOLLSR + ZCOLLR * ZFUNC
        END IF
!
! Otherwise ZDRMAX = 0.0 so the density of the graupel cannot be reached 
!                    and so PRRCOLSS(JLBDAS,JLBDAR) = 0.0        !
!
      END IF
    END DO
!
!*       1.11    Scale the fall speed difference
!
    IF( ZSCALSR>0.0 ) PRRCOLSS(JLBDAS,JLBDAR) = ZCOLLSR / ZSCALSR
  END DO
END DO
!
END SUBROUTINE RRCOLSS
