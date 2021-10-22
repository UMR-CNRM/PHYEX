!=======================================================================
!
! *** ISORROPIA CODE
! *** INCLUDE FILE 'ISRPIA.INC'
! *** THIS FILE CONTAINS THE DECLARATIONS OF THE GLOBAL CONSTANTS
!     AND VARIABLES.
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!=======================================================================
!
!      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
      INTEGER :: NCOMP,NIONS,NGASAQ,NSLDS,NPAIR,NZSR,NERRMX
      PARAMETER (NCOMP=5,NIONS=7,NGASAQ=3,NSLDS=9,NPAIR=13,NZSR=100,    &
     &           NERRMX=25)
!
! *** INPUT VARIABLES **************************************************
!
      INTEGER :: IPROB,METSTBL
      REAL(KIND=8) :: W,WAER,TEMP,RH
      COMMON /INPT/ W(NCOMP), WAER(NCOMP), TEMP, RH, IPROB, METSTBL
!
! *** WATER ACTIVITIES OF PURE SALT SOLUTIONS **************************
!
      REAL(KIND=8) :: AWAS,AWSS,AWAC,AWSC,AWAN,AWSN,AWSB,AWAB,AWSA,AWLC
      COMMON /ZSR / AWAS(NZSR), AWSS(NZSR), AWAC(NZSR), AWSC(NZSR),     &
     &              AWAN(NZSR), AWSN(NZSR), AWSB(NZSR), AWAB(NZSR),     &
     &              AWSA(NZSR), AWLC(NZSR)
!
! *** DELIQUESCENCE RELATIVE HUMIDITIES ********************************
!
      INTEGER :: WFTYP
      REAL(KIND=8) :: DRH2SO4,DRNH42S4,DRNAHSO4,DRNACL,DRNANO3
      REAL(KIND=8) :: DRNA2SO4,DRNH4HS4,DRLC,DRNH4NO3,DRNH4CL
      REAL(KIND=8) :: DRMLCAB,DRMLCAS,DRMASAN,DRMG1,DRMG2
      REAL(KIND=8) :: DRMG3,DRMH1,DRMH2,DRMI1,DRMI2
      REAL(KIND=8) :: DRMI3,DRMQ1,DRMR1,DRMR2,DRMR3
      REAL(KIND=8) :: DRMR4,DRMR5,DRMR6,DRMR7,DRMR8
      REAL(KIND=8) :: DRMR9,DRMR10,DRMR11,DRMR12,DRMR13
      COMMON /DRH / DRH2SO4,  DRNH42S4, DRNAHSO4, DRNACL,   DRNANO3,    &
     &              DRNA2SO4, DRNH4HS4, DRLC,     DRNH4NO3, DRNH4CL
      COMMON /MDRH/ DRMLCAB,  DRMLCAS,  DRMASAN,  DRMG1,    DRMG2,      &
     &              DRMG3,    DRMH1,    DRMH2,    DRMI1,    DRMI2,      &
     &              DRMI3,    DRMQ1,    DRMR1,    DRMR2,    DRMR3,      &
     &              DRMR4,    DRMR5,    DRMR6,    DRMR7,    DRMR8,      &
     &              DRMR9,    DRMR10,   DRMR11,   DRMR12,   DRMR13,     &
     &              WFTYP
!
! *** VARIABLES FOR LIQUID AEROSOL PHASE *******************************
!
      REAL(KIND=8) :: MOLAL, MOLALR, M0
      REAL :: IONIC
      REAL(KIND=8) :: GAMA,ZZ,Z,GAMOU,GAMIN,GASAQ,EPSACT
      REAL(KIND=8) :: COH,CHNO3,CHCL,WATER
      INTEGER :: IACALC
      LOGICAL :: CALAOU, CALAIN, FRST, DRYF
      COMMON /IONS/ MOLAL(NIONS), MOLALR(NPAIR), GAMA(NPAIR), ZZ(NPAIR),&
     &              Z(NIONS),     GAMOU(NPAIR),  GAMIN(NPAIR),M0(NPAIR),&
     &              GASAQ(NGASAQ),                                      &
     &              EPSACT,       COH,           CHNO3,       CHCL,     &
     &              WATER,        IONIC,         IACALC,                &
     &              FRST,         CALAIN,        CALAOU,      DRYF
!
! *** VARIABLES FOR SOLID AEROSOL PHASE ********************************
!
      REAL(KIND=8) :: CH2SO4,CNH42S4,CNH4HS4,CNACL,CNA2SO4
      REAL(KIND=8) :: CNANO3,CNH4NO3,CNH4CL,CNAHSO4,CLC
      COMMON /SALT/ CH2SO4,  CNH42S4, CNH4HS4, CNACL,   CNA2SO4,        &
     &              CNANO3,  CNH4NO3, CNH4CL,  CNAHSO4, CLC
!
! *** VARIABLES FOR GAS PHASE ******************************************
!
      REAL(KIND=8) :: GNH3,GHNO3,GHCL
      COMMON /GAS / GNH3, GHNO3, GHCL
!
! *** EQUILIBRIUM CONSTANTS ********************************************
!
      REAL(KIND=8) :: XK1,XK2,XK3,XK4,XK5,XK6,XK7,XK8,XK9,XK10,XK42
      REAL(KIND=8) :: XK11,XK12,XK13,XK14,XKW, XK21,XK22,XK31,XK32,XK41
      COMMON /EQUK/ XK1, XK2, XK3, XK4, XK5, XK6, XK7, XK8, XK9, XK10,  &
     &              XK11,XK12,XK13,XK14,XKW, XK21,XK22,XK31,XK32,XK41,  &
     &              XK42
!
! *** MOLECULAR WEIGHTS ************************************************
!
      REAL(KIND=8) :: IMW
      REAL(KIND=8) :: R,WMW,SMW
      COMMON /OTHR/ R, IMW(NIONS), WMW(NCOMP), SMW(NPAIR)
!
! *** SOLUTION/INFO VARIABLES ******************************************
!
      CHARACTER SCASE*15
      REAL(KIND=8) :: SULRATW,SULRAT,SODRAT
      COMMON /CASE/ SULRATW, SULRAT, SODRAT, SCASE
!
      REAL(KIND=8) :: EPS
      INTEGER :: MAXIT,NSWEEP,NDIV,ICLACT
      COMMON /SOLN/ EPS, MAXIT, NSWEEP, NDIV, ICLACT
!
! *** ERROR SYSTEM *****************************************************
!
      CHARACTER ERRMSG*40
      INTEGER   ERRSTK, NOFER
      LOGICAL   STKOFL
      COMMON /EROR/ STKOFL, NOFER, ERRSTK(NERRMX), ERRMSG(NERRMX)
!
! *** GENERIC VARIABLES ************************************************
!
      CHARACTER VERSION*14
      REAL(KIND=8) :: GREAT,TINY,TINY2,ZERO,ONE
      COMMON /CGEN/ GREAT, TINY, TINY2, ZERO, ONE, VERSION
!
! *** END OF INCLUDE FILE **********************************************
!
