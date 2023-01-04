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
      MODULE MODI_CONVECT_TSTEP_PREF
!     #################
!
INTERFACE
!
       SUBROUTINE CONVECT_TSTEP_PREF( KLON, KLEV,                           &
                                     PU, PV, PPRES, PZ, PDXDY, KLCL, KCTL, &
                                     PTIMEA, PPREF )
!
INTEGER, INTENT(IN)                    :: KLON   ! horizontal dimension
INTEGER, INTENT(IN)                    :: KLEV   ! vertical dimension
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PPRES  ! pressure (Pa) 
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PU     ! grid scale horiz. wind u 
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PV     ! grid scale horiz. wind v
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PZ     ! height of model layer (m) 
REAL, DIMENSION(KLON),      INTENT(IN) :: PDXDY  ! grid area (m^2)
INTEGER, DIMENSION(KLON),   INTENT(IN) :: KLCL   ! lifting condensation level index
INTEGER, DIMENSION(KLON),   INTENT(IN) :: KCTL   ! cloud top level index
!
REAL, DIMENSION(KLON),      INTENT(OUT):: PTIMEA ! advective time period
REAL, DIMENSION(KLON),      INTENT(OUT):: PPREF  ! precipitation efficiency 
!
END SUBROUTINE CONVECT_TSTEP_PREF
!
END INTERFACE
!
END MODULE MODI_CONVECT_TSTEP_PREF
!     ######################################################################
      SUBROUTINE CONVECT_TSTEP_PREF( KLON, KLEV,                           &
                                     PU, PV, PPRES, PZ, PDXDY, KLCL, KCTL, &
                                     PTIMEA, PPREF )
!     ######################################################################
!
!!**** Routine to compute convective advection time step and precipitation 
!!     efficiency 
!!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to determine the convective
!!      advection time step PTIMEC as a function of the mean ambient 
!!      wind as well as the precipitation efficiency as a function
!!      of wind shear and cloud base height.
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
!!
!!     Module MODD_CONVPAREXT
!!          JCVEXB, JCVEXT     ! extra levels on the vertical boundaries
!!
!!    REFERENCE
!!    ---------
!!
!!      Book1,2 of documentation 
!!      Fritsch and Chappell, 1980, J. Atmos. Sci.
!!      Kain and Fritsch, 1993, Meteor. Monographs, Vol.
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    07/11/95 
!!   Last modified  04/10/97
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CONVPAREXT
!
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER, INTENT(IN)                    :: KLON   ! horizontal dimension
INTEGER, INTENT(IN)                    :: KLEV   ! vertical dimension
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PPRES  ! pressure (Pa) 
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PU     ! grid scale horiz. wind u 
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PV     ! grid scale horiz. wind v
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PZ     ! height of model layer (m) 
REAL, DIMENSION(KLON),      INTENT(IN) :: PDXDY  ! grid area (m^2)
INTEGER, DIMENSION(KLON),   INTENT(IN) :: KLCL   ! lifting condensation level index
INTEGER, DIMENSION(KLON),   INTENT(IN) :: KCTL   ! cloud top level index
!
REAL, DIMENSION(KLON),      INTENT(OUT):: PTIMEA ! advective time period
REAL, DIMENSION(KLON),      INTENT(OUT):: PPREF  ! precipitation efficiency 
!
!
!*       0.2   Declarations of local variables KLON
!
INTEGER :: IIE, IKB, IKE                      ! horizontal + vertical loop bounds
INTEGER :: JI                                 ! horizontal loop index
INTEGER :: JK, JKLC, JKP5, JKCT               ! vertical loop index
!
INTEGER, DIMENSION(KLON)  :: IP500       ! index of 500 hPa levels
REAL, DIMENSION(KLON)     :: ZCBH        ! cloud base height 
REAL, DIMENSION(KLON)     :: ZWORK1, ZWORK2, ZWORK3  ! work arrays
!
!
!-------------------------------------------------------------------------------
!
!        0.3   Set loop bounds
!              ---------------
!
IIE = KLON
IKB = 1 + JCVEXB 
IKE = KLEV - JCVEXT 
!
!
!*       1.     Determine vertical index for 500 hPa levels 
!               ------------------------------------------
!
!
IP500(:) = IKB
DO JK = IKB, IKE
    WHERE ( PPRES(:,JK) >= 500.E2 ) IP500(:) = JK
END DO
!
!
!*       2.     Compute convective time step 
!               ----------------------------
!
	    ! compute wind speed at LCL, 500 hPa, CTL

DO JI = 1, IIE
   JKLC = KLCL(JI)
   JKP5 = IP500(JI)
   JKCT = KCTL(JI)
   ZWORK1(JI) = SQRT( PU(JI,JKLC) * PU(JI,JKLC) +           &
		      PV(JI,JKLC) * PV(JI,JKLC)  ) 
   ZWORK2(JI) = SQRT( PU(JI,JKP5) * PU(JI,JKP5) +           &
		      PV(JI,JKP5) * PV(JI,JKP5)  ) 
   ZWORK3(JI) = SQRT( PU(JI,JKCT) * PU(JI,JKCT) +           &
		      PV(JI,JKCT) * PV(JI,JKCT)  ) 
END DO
!
ZWORK2(:) = MAX( 0.1, 0.5 * ( ZWORK1(:) + ZWORK2(:) ) )
!
PTIMEA(:) = SQRT( PDXDY(:) ) / ZWORK2(:) 
!
!
!*       3.     Compute precipitation efficiency 
!               -----------------------------------
!
!*       3.1    Precipitation efficiency as a function of wind shear
!               ----------------------------------------------------
!
ZWORK2(:) = SIGN( 1., ZWORK3(:) - ZWORK1(:) )
DO JI = 1, IIE
    JKLC = KLCL(JI)
    JKCT = KCTL(JI)
    ZWORK1(JI) = ( PU(JI,JKCT) - PU(JI,JKLC) )  *          &
                 ( PU(JI,JKCT) - PU(JI,JKLC) )  +          &
                 ( PV(JI,JKCT) - PV(JI,JKLC) )  *          &
                 ( PV(JI,JKCT) - PV(JI,JKLC) )  
    ZWORK1(JI) = 1.E3 * ZWORK2(JI) * SQRT( ZWORK1(JI) ) /  &
	         MAX( 1.E-2, PZ(JI,JKCT) - PZ(JI,JKLC) )
END DO
!
PPREF(:)  = 1.591 + ZWORK1(:) * ( -.639 + ZWORK1(:) * (        &
				9.53E-2 - ZWORK1(:) * 4.96E-3 ) ) 
PPREF(:)  = MAX( .4, MIN( PPREF(:), .9 ) )
!
!*       3.2    Precipitation efficiency as a function of cloud base height 
!               ----------------------------------------------------------
!
DO JI = 1, IIE
   JKLC = KLCL(JI)
   ZCBH(JI)   = MAX( 3., ( PZ(JI,JKLC) - PZ(JI,IKB) ) * 3.281E-3 ) 
END DO
ZWORK1(:) = .9673 + ZCBH(:) * ( -.7003 + ZCBH(:) * ( .1622 + &
	      ZCBH(:) *  ( -1.2570E-2 + ZCBH(:) * ( 4.2772E-4 -  &
              ZCBH(:) * 5.44E-6 ) ) ) )
ZWORK1(:) = MAX( .4, MIN( .9, 1./ ( 1. + ZWORK1(:) ) ) )
!
!*       3.3    Mean precipitation efficiency is used to compute rainfall 
!               ----------------------------------------------------------
!
PPREF(:) = 0.5 * ( PPREF(:) + ZWORK1(:) )
!
!
END SUBROUTINE CONVECT_TSTEP_PREF
