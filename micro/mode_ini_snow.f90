!     ######spl
MODULE MODE_INI_SNOW
IMPLICIT NONE
CONTAINS
      SUBROUTINE INI_SNOW ( KLUOUT )
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
!     ###########################################################
!
!!****  *INI_SNOW * - re-initialize the constants based on snow-size distubutio
!!                        cold microphysical schemes.
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to reinitialize the constants for snow used to
!!    resolve the mixed phase microphysical scheme.
!!    EXTERNAL
!!    --------
!!      GAMMA    :  gamma function
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST
!!        XPI                  !
!!        XP00                 ! Reference pressure
!!        XRD                  ! Gaz constant for dry air
!!        XRHOLW               ! Liquid water density
!!      Module MODD_REF
!!        XTHVREFZ             ! Reference virtual pot.temp. without orography
!!      Module MODD_PARAMETERS
!!        JPVEXT               !
!!      Module MODD_RAIN_ICE_DESCR
!!      Module MODD_RAIN_ICE_PARAM
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation ( routine INI_RAIN_ICE )
!!
!!    ORIGINAL AUTHOR (from ini_rain_ice)
!!    --------------
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      2018-02
!!      Karl-Ivar Ivarsson
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
USE MODD_PARAMETERS
USE MODD_PARAM_ICE_n
USE MODD_RAIN_ICE_DESCR_n
USE MODD_RAIN_ICE_PARAM_n
!
USE MODI_GAMMA
USE MODI_GAMMA_INC
USE MODE_RRCOLSS
USE MODE_RZCOLX
USE MODE_RSCOLRG
USE MODE_READ_XKER_RACCS
USE MODE_READ_XKER_SDRYG
USE MODE_READ_XKER_RDRYG
USE MODE_READ_XKER_SWETH
USE MODE_READ_XKER_GWETH


IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
INTEGER,                 INTENT(IN) :: KLUOUT   ! Logical unit number for prints
!*       0.2   Declarations of local variables :
!
REAL :: ZRHO00                ! Surface reference air density

REAL :: ZCONC_MAX ! Maximal concentration for snow


REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('INI_RAIN_ICE',0,ZHOOK_HANDLE)


XCCS = XFRMIN(16)
XCXS = XFRMIN(17)
ZRHO00 = XP00/(XRD*300.0)
!     recalculate ini_rain_ice stuff:

!     3.4    Constants for shape parameter
XLBEXS = 1.0/(XCXS-XBS)
XLBS   = ( XAS*XCCS*MOMG(XALPHAS,XNUS,XBS) )**(-XLBEXS)
ZCONC_MAX  = 1.E6 ! Maximal concentration for falling particules set to 1 per cc
IF(XCCS>0. .AND. XCXS>0. )XLBDAS_MAX = ( ZCONC_MAX/XCCS )**(1./XCXS)

!     4.2    Constants for sedimentation
XEXSEDS = (XBS+XDS-XCXS)/(XBS-XCXS)

XFSEDS  = XCS*XAS*XCCS*MOMG(XALPHAS,XNUS,XBS+XDS)*                         &
     (XAS*XCCS*MOMG(XALPHAS,XNUS,XBS))**(-XEXSEDS)*(ZRHO00)**XCEXVT

!     5.2    Constants for vapor deposition on ice
X0DEPS = (4.0*XPI)*XCCS*XC1S*XF0S*MOMG(XALPHAS,XNUS,1.)
X1DEPS = (4.0*XPI)*XCCS*XC1S*XF1S*SQRT(XCS)*MOMG(XALPHAS,XNUS,0.5*XDS+1.5)
XEX0DEPS = XCXS-1.0
XEX1DEPS = XCXS-0.5*(XDS+3.0)

!     5.4    Constants for snow aggregation
XFIAGGS  = (XPI/4.0)*XCOLIS*XCCS*XCS*(ZRHO00**XCEXVT)*MOMG(XALPHAS,XNUS,XDS+2.0)
XEXIAGGS = XCXS-XDS-2.0

!     7.1    Constants for the riming of the aggregates
XEXCRIMSS= XCXS-XDS-2.0
XCRIMSS  = (XPI/4.0)*XCOLCS*XCCS*XCS*(ZRHO00**XCEXVT)*MOMG(XALPHAS,XNUS,XDS+2.0)
XEXCRIMSG= XEXCRIMSS
XCRIMSG  = XCRIMSS
XSRIMCG  = XCCS*XAS*MOMG(XALPHAS,XNUS,XBS)
XEXSRIMCG= XCXS-XBS

!     7.2    Constants for the accretion of raindrops onto aggregates

XFRACCSS = ((XPI**2)/24.0)*XCCS*XCCR*XRHOLW*(ZRHO00**XCEXVT)

XFSACCRG = (XPI/4.0)*XAS*XCCS*XCCR*(ZRHO00**XCEXVT)

!     8.2.3  Constants for the aggregate collection by the graupeln
XFSDRYG = (XPI/4.0)*XCOLSG*XCCG*XCCS*XAS*(ZRHO00**XCEXVT)

!     9.2.2  Constants for the aggregate collection by the hailstones
XFSWETH = (XPI/4.0)*XCCH*XCCS*XAS*(ZRHO00**XCEXVT)

WRITE(UNIT=KLUOUT,FMT='("  updated snow concentration:C=",E13.6," x=",E13.6)') &
                                                      XCCS,XCXS

IF (LHOOK) CALL DR_HOOK('INI_SNOW',1,ZHOOK_HANDLE)

CONTAINS
!
!------------------------------------------------------------------------------
!
  FUNCTION MOMG(PALPHA,PNU,PP) RESULT (PMOMG)
!
! auxiliary routine used to compute the Pth moment order of the generalized
! gamma law
!
  USE MODI_GAMMA
!
  IMPLICIT NONE
!
  REAL, INTENT(IN)     :: PALPHA ! first shape parameter of the dimensionnal distribution
  REAL, INTENT(IN)     :: PNU    ! second shape parameter of the dimensionnal distribution
  REAL, INTENT(IN)     :: PP     ! order of the moment
  REAL     :: PMOMG  ! result: moment of order ZP
!
!------------------------------------------------------------------------------
!
!
  PMOMG = GAMMA(PNU+PP/PALPHA)/GAMMA(PNU)
!
  END FUNCTION MOMG
!
!-------------------------------------------------------------------------------
!
!
END SUBROUTINE INI_SNOW
END MODULE MODE_INI_SNOW
