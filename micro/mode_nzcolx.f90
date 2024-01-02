!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_NZCOLX
  IMPLICIT NONE
CONTAINS
!     ################################################################
  SUBROUTINE NZCOLX( KND, PALPHAX, PNUX, PALPHAZ, PNUZ,          &
                     PEXZ, PFALLX, PEXFALLX, PFALLEXPX,          &
                     PFALLZ, PEXFALLZ, PFALLEXPZ,                &
                     PLBDAXMAX, PLBDAZMAX, PLBDAXMIN, PLBDAZMIN, &
                     PDINFTY, PNZCOLX                            )
!     ################################################################
!
!
!
!!****  * -  Build up a look-up table containing the scaled fall speed
!!           difference between size distributed particles of specy X and Z
!!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to integrate numerically the scaled fall
!!      speed difference between specy X and specy Z for use in collection
!!      kernels. A first integral of the form
!!
!!        infty
!!       / /
!!       |{|                                                 }
!!       |{| E_xz (Dx+Dz)^2 |cxDx^dx-czDz^dz| g(Dz) dDz} g(Dx) dDx
!!       |{|                                                 }
!!       / /
!!        0
!!
!!      is evaluated and normalised by a second integral of the form
!!
!!        infty
!!       / /
!!       |{|                          }
!!       |{| (Dx+Dz)^2 g(Dz) dDz} g(Dx) dDx
!!       |{|                          }
!!       / /
!!        0
!!
!!      where E_xz is a collection efficiency, g(D) is the generalized Gamma
!!      distribution law. The 'infty' diameter is defined according to the
!!      current value of the Lbda that is D_x=PDINFTY/Lbda_x or
!!      D_z=PINFTY/Lbda_z. 
!!      The result is stored in a two-dimensional array.
!! 
!!**  METHOD
!!    ------
!!      The free parameters of the size distribution function of specy X and Z
!!      (slope parameter LAMBDA) are discretized with a geometrical rate in a
!!      specific range
!!            LAMBDA = exp( (Log(LAMBDA_max) - Log(LAMBDA_min))/N_interval )
!!      The two above integrals are performed using the trapezoidal scheme and
!!      the [0,infty] interval is discretized over KND values of D_x or D_z.
!!
!!    EXTERNAL
!!    --------
!!      MODI_GENERAL_GAMMA: Generalized gamma distribution law 
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      None
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
!!      M. Taufour     03/2022: adapted from rzcolx for concentration
!!      J. Wurtz       03/2022: new snow characteristics
!!
!-------------------------------------------------------------------------------
!
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODI_GENERAL_GAMMA
!
IMPLICIT NONE
!
!
!*       0.1   Declarations of dummy arguments 
!              ------------------------------- 
!
!
INTEGER, INTENT(IN) :: KND    ! Number of discrete size intervals in DX and DZ  
!
!
REAL, INTENT(IN) :: PALPHAX   ! First shape parameter of the specy X 
			      ! size distribution (generalized gamma law)
REAL, INTENT(IN) :: PNUX      ! Second shape parameter of the specy X
			      ! size distribution (generalized gamma law)
REAL, INTENT(IN) :: PALPHAZ   ! First shape parameter of the specy Z 
			      ! size distribution (generalized gamma law)
REAL, INTENT(IN) :: PNUZ      ! Second shape parameter of the specy Z
			      ! size distribution (generalized gamma law)
REAL, INTENT(IN) :: PEXZ      ! Efficiency of specy X collecting specy Z
REAL, INTENT(IN) :: PFALLX    ! Fall speed constant of specy X
REAL, INTENT(IN) :: PEXFALLX  ! Fall speed exponent of specy X
REAL, INTENT(IN) :: PFALLEXPX ! Fall speed exponential constant of specy X
REAL, INTENT(IN) :: PFALLZ    ! Fall speed constant of specy Z
REAL, INTENT(IN) :: PEXFALLZ  ! Fall speed exponent of specy Z
REAL, INTENT(IN) :: PFALLEXPZ ! Fall speed exponential constant of specy Z
REAL, INTENT(IN) :: PLBDAXMAX ! Maximun slope of size distribution of specy X
REAL, INTENT(IN) :: PLBDAZMAX ! Maximun slope of size distribution of specy Z
REAL, INTENT(IN) :: PLBDAXMIN ! Minimun slope of size distribution of specy X
REAL, INTENT(IN) :: PLBDAZMIN ! Minimun slope of size distribution of specy Z
REAL, INTENT(IN) :: PDINFTY   ! Factor to define the largest diameter up to
			      ! which the diameter integration is performed
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PNZCOLX ! Scaled fall speed difference in
				               ! the mass collection kernel as a
					       ! function of LAMBDAX and LAMBDAZ
!
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
!
INTEGER :: JLBDAX  ! Slope index of the size distribution of specy X
INTEGER :: JLBDAZ  ! Slope index of the size distribution of specy Z
INTEGER :: JDX     ! Diameter index of a particle of specy X
INTEGER :: JDZ     ! Diameter index of a particle of specy Z
!
!
REAL    :: ZLBDAX  ! Current slope parameter LAMBDA of specy X
REAL    :: ZLBDAZ  ! Current slope parameter LAMBDA of specy Z
REAL    :: ZDLBDAX ! Growth rate of the slope parameter LAMBDA of specy X
REAL    :: ZDLBDAZ ! Growth rate of the slope parameter LAMBDA of specy Z
REAL    :: ZDDX    ! Integration step of the diameter of specy X
REAL    :: ZDDZ    ! Integration step of the diameter of specy Z
REAL    :: ZDX     ! Current diameter of the particle specy X
REAL    :: ZDZ     ! Current diameter of the particle specy Z
REAL    :: ZCOLLZ  ! Single integral of the mass weighted fall speed difference 
		   ! over the spectrum of specy Z
REAL    :: ZCOLLXZ ! Double integral of the mass weighted fall speed difference
		   ! over the spectra of specy X and specy Z
REAL    :: ZSCALZ  ! Single integral of the scaling factor over 
		   ! the spectrum of specy Z
REAL    :: ZSCALXZ ! Double integral of the scaling factor over
		   ! the spectra of specy X and specy Z
REAL    :: ZFUNC   ! Ancillary function
!
!
!-------------------------------------------------------------------------------
!
!
!*       1       COMPUTE THE SCALED VELOCITZ DIFFERENCE IN THE MASS
!*                               COLLECTION KERNEL,
!                -------------------------------------------------
!
!
!
!*       1.1     Compute the growth rate of the slope factors LAMBDA
!
ZDLBDAX = EXP( LOG(PLBDAXMAX/PLBDAXMIN)/REAL(SIZE(PNZCOLX(:,:),1)-1) )
ZDLBDAZ = EXP( LOG(PLBDAZMAX/PLBDAZMIN)/REAL(SIZE(PNZCOLX(:,:),2)-1) )
!
!*       1.2     Scan the slope factors LAMBDAX and LAMBDAZ
!
DO JLBDAX = 1,SIZE(PNZCOLX(:,:),1)
  ZLBDAX = PLBDAXMIN * ZDLBDAX ** (JLBDAX-1) 
  DO JLBDAZ = 1,SIZE(PNZCOLX(:,:),2)
    ZLBDAZ = PLBDAZMIN * ZDLBDAZ ** (JLBDAZ-1)
!
!*       1.3     Initialize the collection integrals
!
    ZSCALXZ = 0.0
    ZCOLLXZ = 0.0
!
!*       1.4     Compute the diameter steps
!
    ZDDX = PDINFTY / (REAL(KND) * ZLBDAX)
    ZDDZ = PDINFTY / (REAL(KND) * ZLBDAZ)
!
!*       1.5     Scan over the diameters DX and DZ
!
    DO JDX = 1,KND-1
      ZDX = ZDDX * REAL(JDX)
!
      ZSCALZ = 0.0
      ZCOLLZ = 0.0
      DO JDZ = 1,KND-1
        ZDZ = ZDDZ * REAL(JDZ)
!
!*       1.6     Compute the normalization factor by integration over the
!                dimensional spectrum of specy Z  
!
	ZFUNC  = (ZDX+ZDZ)**2 * GENERAL_GAMMA(PALPHAZ,PNUZ,ZLBDAZ,ZDZ)
	ZSCALZ = ZSCALZ + ZFUNC
!
!*       1.7     Compute the scaled fall speed difference by integration over
!                the dimensional spectrum of specy Z
!
         ZCOLLZ = ZCOLLZ + ZFUNC * PEXZ * ABS( PFALLX*ZDX**PEXFALLX * EXP(-(ZDX*PFALLEXPX)**PALPHAX) &
                                             - PFALLZ*ZDZ**PEXFALLZ * EXP(-(ZDZ*PFALLEXPZ)**PALPHAZ))
      END DO
!
!*       1.8     Compute the normalization factor by integration over the
!                dimensional spectrum of specy X
!
      ZFUNC   = GENERAL_GAMMA(PALPHAX,PNUX,ZLBDAX,ZDX)
      ZSCALXZ = ZSCALXZ + ZSCALZ * ZFUNC
!
!*       1.9     Compute the scaled fall speed difference by integration over
!                the dimensional spectrum of specy X
!
      ZCOLLXZ = ZCOLLXZ + ZCOLLZ * ZFUNC
    END DO
!
!*       1.10    Scale the fall speed difference
!
    PNZCOLX(JLBDAX,JLBDAZ) = ZCOLLXZ / ZSCALXZ
  END DO
END DO
!
END SUBROUTINE NZCOLX
END MODULE MODE_NZCOLX
