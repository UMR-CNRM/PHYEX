!     ######spl
SUBROUTINE AROINI_CONF(KULOUT,OWTOTL,OCARTESIAN)
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK

!**** INI_CSTMNH *  - Initiallize MesoNH Physics configuration module
!**   Interface.
!     ----------
!        *CALL* *INI_CONF(KULOUT)

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
!       the variables are initialised in the module itself
!        1.2 Set implicit default values for MODD_CST
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('AROINI_CONF',0,ZHOOK_HANDLE)

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

!   ---------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('AROINI_CONF',1,ZHOOK_HANDLE)
RETURN
END SUBROUTINE AROINI_CONF 
