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
USE MODD_CONF

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

!   ---------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('AROINI_CSTMNH',1,ZHOOK_HANDLE)
RETURN
END SUBROUTINE AROINI_CSTMNH 
