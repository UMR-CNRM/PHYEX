!     ######spl
SUBROUTINE AROINI_CSTMNH(KULOUT,OWTOTL,OCARTESIAN)
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK

!**** INI_CSTMNH *  - Initiallize MesoNH Physics constantes
!**   Interface.
!     ----------
!        *CALL* *INI_CSTMNH(KULOUT)

!        Explicit arguments :
!        --------------------
!        KULOUT : Logical unit for the output
!        OWTWTL : T if SL2TL scheme is used in AROME
!        OCARTESIAN : T for academic cases

!        Implicit arguments :
!        --------------------
!        

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        Documentation AROME 

!     Author.
!     -------
!        Y. Seity 

!     Modifications.
!     --------------
!        Original : 03-12-12
!     ------------------------------------------------------------------


USE MODD_PARAMETERS 
USE MODD_CST
USE MODD_CONF
USE MODD_LUNIT
USE MODI_INI_CST

IMPLICIT NONE

!   ---------------------------------------------------------------
!     DUMMY INTEGER SCALARS
INTEGER, INTENT(IN) :: KULOUT
LOGICAL, INTENT(IN) :: OWTOTL
LOGICAL, INTENT(IN) :: OCARTESIAN
!   ---------------------------------------------------------------
!*       1.    Set default values.
!              -------------------


!        1.1 Set implicit default values for MODD_PARAMETERS
!       les variables sont initialisées dans le module lui même 
!        1.2 Set implicit default values for MODD_CST
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('AROINI_CSTMNH',0,ZHOOK_HANDLE)
CALL INI_CST

!        1.3 Set implicit default values for MODD_CONF
IF (OWTOTL) THEN
  CCONF='RESTA'       !for 2TLSL scheme
ELSE
  CCONF='START'
ENDIF

  LTHINSHELL  = .TRUE.
  L2D        = .FALSE.
  L1D        = .FALSE.
  LFLAT      = OCARTESIAN
  NMODEL     = 1
  CEQNSYS    = 'DUR'
  CEXP       = 'AROME'
  CSEG       = 'SEG01'
  LLG        = .FALSE.
  LINIT_LG   = .FALSE.
  LNOMIXLG   = .FALSE.
  LCARTESIAN = OCARTESIAN
  CPROGRAM   = 'AROME '

!*      2.    Print final values.
!              -------------------
WRITE(UNIT=KULOUT,FMT='('' COMMON MODD_CONF Meso_NH '')')
WRITE(UNIT=KULOUT,FMT='('' LCARTESIAN = '',L2,'' CPROGRAM = '',A6,&
     &'' CCONF = '',A5,/, '' LFLAT = '',L2,'' L1D = '',L2, '' L2D = '',L2)')&
     &LCARTESIAN,CPROGRAM,CCONF,LFLAT,L1D,L2D

WRITE(UNIT=KULOUT,FMT='('' COMMON MODD_CST Meso_NH '')')
WRITE(UNIT=KULOUT,FMT='('' XPI = '',E10.4,'' XKARMAN = '',E10.4,&
     &'' XLIGHTSPEED = '',E10.4 &
     &,/,'' XPLANCK = '',E10.4 &
     &,'' XBOLTZ = '',E10.4,'' XAVOGADRO = '',E10.4,/&
     &,'' XDAY = '',E10.4 &
     &,'' XSIYEA = '',E10.4,'' XSIDAY = '',E10.4,/&
     &,'' XOMEGA = '',E10.4 &
     &,'' XRADIUS = '',E10.4,'' XG = '',E10.4,/&
     &,'' XP00 = '',E10.4 &
     &,'' XT00 = '',E10.4,'' XSTEFAN = '',E10.4,/&
     &,'' XIO = '',E10.4 &
     &,'' XMD = '',E10.4,'' XMV = '',E10.4,/&
     &,'' XRD = '',E10.4 &
     &,'' XRV = '',E10.4,'' XCPD = '',E10.4,/&
     &,'' XCPV = '',E10.4 &
     &,'' XRHOLW = '',E10.4,'' XRHOLI = '',E10.4,/&
     &,'' XCL = '',E10.4 &
     &,'' XCI = '',E10.4,'' XTT = '',E10.4,/&
     &,'' XLVTT = '',E10.4 &
     &,'' XLSTT = '',E10.4,'' XLMTT = '',E10.4,/&
     &,'' XESTT = '',E10.4 &
     &,'' XGAMW = '',E10.4,'' XBETAW = '',E10.4,/&
     &,'' XALPW = '',E10.4 &
     &,'' XGAMI = '',E10.4,'' XBETAI = '',E10.4,/&
     &,'' XALPI = '',E10.4)')&
 &XPI,XKARMAN,XLIGHTSPEED,XPLANCK,XBOLTZ,XAVOGADRO,XDAY,XSIYEA,&
 &XSIDAY,XOMEGA,XRADIUS,XG,XP00,XTH00,XSTEFAN,XI0,XMD,XMV,XRD,&
 &XRV,XCPD,XCPV,XRHOLW,XRHOLI,XCL,XCI,XTT,XLVTT,XLSTT,XLMTT,XESTT,&
 &XGAMW,XBETAW,XALPW,XGAMI,XBETAI,XALPI 

!   ---------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('AROINI_CSTMNH',1,ZHOOK_HANDLE)
RETURN
END SUBROUTINE AROINI_CSTMNH 
