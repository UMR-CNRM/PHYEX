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
      MODULE MODI_CONVECT_MIXING_FUNCT
!     #################
!
INTERFACE
!
       SUBROUTINE CONVECT_MIXING_FUNCT( KLON,                &
                                       PMIXC, KMF, PER, PDR )
!
INTEGER,               INTENT(IN) :: KLON   ! horizontal dimension
INTEGER,               INTENT(IN) :: KMF    ! switch for dist. function
REAL, DIMENSION(KLON), INTENT(IN) :: PMIXC  ! critical mixed fraction
!
REAL, DIMENSION(KLON), INTENT(OUT):: PER    ! normalized entrainment rate
REAL, DIMENSION(KLON), INTENT(OUT):: PDR    ! normalized detrainment rate
!
END SUBROUTINE CONVECT_MIXING_FUNCT
!
END INTERFACE
!
END MODULE MODI_CONVECT_MIXING_FUNCT
!     ######spl
      SUBROUTINE CONVECT_MIXING_FUNCT( KLON,                &
                                       PMIXC, KMF, PER, PDR ) 
!     #######################################################
!
!!**** Determine the area under the distribution function
!!     KMF = 1 : gaussian  KMF = 2 : triangular distribution function
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to determine the entrainment and
!!      detrainment rate by evaluating the are under the distribution 
!!      function. The integration interval is limited by the critical
!!      mixed fraction PMIXC
!!   
!!
!!
!!**  METHOD
!!    ------
!!      Use handbook of mathemat. functions by Abramowitz and Stegun, 1968
!!      
!!     
!!
!!    EXTERNAL
!!    --------
!!      None
!!     
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      None
!!
!!
!!    REFERENCE
!!    ---------
!!
!!      Book2 of documentation ( routine MIXING_FUNCT)
!!      Abramovitz and Stegun (1968), handbook of math. functions 
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
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER,               INTENT(IN) :: KLON   ! horizontal dimension
INTEGER,               INTENT(IN) :: KMF    ! switch for dist. function
REAL, DIMENSION(KLON), INTENT(IN) :: PMIXC  ! critical mixed fraction
!
REAL, DIMENSION(KLON), INTENT(OUT):: PER    ! normalized entrainment rate
REAL, DIMENSION(KLON), INTENT(OUT):: PDR    ! normalized detrainment rate
!
!*       0.2   Declarations of local variables :
!
REAL    :: ZSIGMA = 0.166666667                   ! standard deviation 
REAL    :: ZFE    = 4.931813949                   ! integral normalization 
REAL    :: ZSQRTP = 2.506628,  ZP  = 0.33267      ! constants
REAL    :: ZA1    = 0.4361836, ZA2 =-0.1201676    ! constants
REAL    :: ZA3    = 0.9372980, ZT1 = 0.500498     ! constants
REAL    :: ZE45   = 0.01111                       ! constant
!
REAL, DIMENSION(KLON) :: ZX, ZY, ZW1, ZW2         ! work variables
REAL    :: ZW11
!
!
!-------------------------------------------------------------------------------
!
!       1.     Use gaussian function for KMF=1
!              -------------------------------
!
IF( KMF == 1 ) THEN 
    ! ZX(:)  = ( PMIXC(:) - 0.5 ) / ZSIGMA
      ZX(:)  = 6. * PMIXC(:) - 3.
      ZW1(:) = 1. / ( 1.+ ZP * ABS ( ZX(:) ) )
      ZY(:)  = EXP( -0.5 * ZX(:) * ZX(:) )
      ZW2(:) = ZA1 * ZW1(:) + ZA2 * ZW1(:) * ZW1(:) +                   &
		 ZA3 * ZW1(:) * ZW1(:) * ZW1(:)
      ZW11   = ZA1 * ZT1 + ZA2 * ZT1 * ZT1 + ZA3 * ZT1 * ZT1 * ZT1
ENDIF 
!
WHERE ( KMF == 1 .AND. ZX(:) >= 0. )
	PER(:) = ZSIGMA * ( 0.5 * ( ZSQRTP - ZE45 * ZW11                 &
		 - ZY(:) * ZW2(:) ) + ZSIGMA * ( ZE45 - ZY(:) ) )        &
		 - 0.5 * ZE45 * PMIXC(:) * PMIXC(:)
	PDR(:) = ZSIGMA*( 0.5 * ( ZY(:) * ZW2(:) - ZE45 * ZW11   )       &
		 + ZSIGMA * ( ZE45 - ZY(:) ) )                           &
		 - ZE45 * ( 0.5 + 0.5 * PMIXC(:) * PMIXC(:) - PMIXC(:) )
END WHERE
WHERE ( KMF == 1 .AND. ZX(:) < 0. ) 
	PER(:) = ZSIGMA*( 0.5 * ( ZY(:) * ZW2(:) - ZE45 * ZW11   )       &
		 + ZSIGMA * ( ZE45 - ZY(:) ) )                           &
		 - 0.5 * ZE45 * PMIXC(:) * PMIXC(:)
	PDR(:) = ZSIGMA * ( 0.5 * ( ZSQRTP - ZE45 * ZW11 - ZY(:)         &
		 * ZW2(:) ) + ZSIGMA * ( ZE45 - ZY(:) ) )                &
		 - ZE45 * ( 0.5 + 0.5 * PMIXC(:) * PMIXC(:) - PMIXC(:) )
END WHERE
!
      PER(:) = PER(:) * ZFE
      PDR(:) = PDR(:) * ZFE
!
!
!       2.     Use triangular function KMF=2
!              -------------------------------
!
!     not yet released
!
!
END SUBROUTINE CONVECT_MIXING_FUNCT
