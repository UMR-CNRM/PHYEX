!     ######spl
SUBROUTINE AROINI_MICRO(KULOUT,PTSTEP,LDWARM,CMICRO,KSPLITR,CCSEDIM,LDCRIAUTI,&
              PCRIAUTI,PT0CRIAUTI,PCRIAUTC,PTSTEP_TS, CCSNOWRIMING, PMRSTEP, KMAXITER, &
              LDFEEDBACKT, LDEVLIMIT, LDNULLWETG, LDWETGPOST, LDNULLWETH, LDWETHPOST, &
              PFRACM90, LDCONVHG, CCSUBG_RC_RR_ACCR, CCSUBG_RR_EVAP, CCSUBG_PR_PDF, &
              LDCRFLIMIT, CCFRAC_ICE_ADJUST, PSPLIT_MAXCFL,&
              CCFRAC_ICE_SHALLOW_MF, LDSEDIM_AFTER,LDDEPOSC,PVDEPOSC, PFRMIN,&
              LDDEPSG,PRDEPSRED,PRDEPGRED)

USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!**** *INI_MICRO*   - Initialize common meso_NH MODD_ used in microphysics for AROME

!     Purpose.
!     --------
!           Initialize 
!           MODD_RAIN_ICE_DESCR, MODD_RAIN_ICE_PARAM and MODD_PARAM_ICE  
!           parameters used in AROME microphysics 

!**   Interface.
!     ----------
!        *CALL* *INI_MICRO (KULOUT,KSTEP,KSPLITR)

!        Explicit arguments :
!        --------------------
!        KULOUT : Logical unit for the output
!        PTSTEP  : Time step
!        KSPLITR : Number of small time step interation for rain sedimentation 
!        LDWARM : value assigned to LWARM       

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
!        05-08-25 Kovacic  added LDWARM
!        Jan 2015 S. Riette: LFEEDBACKT, LEVLIMIT, LNULLWETG, LWETGPOST, CSNOWRIMING,
!                            XFRACM90, NMRSITER, XMRSTEP, LSIMULSG, XTSTEP_TS
!                            LNULLWETH, LWETHPOST added
!        Oct 2016 S. Riette: LDCRFLIMIT, CCFRAC_ICE_ADJUST
!                            and CCFRAC_ICE_SHALLOW_MF added
!        Dec 2020 Y. Seity : Add Fog deposition term
!        Jan 2020 C.Wittmann: Add LDDEPSG,PRDEPSRED,PRDEPGRED
!     ------------------------------------------------------------------

USE MODD_RAIN_ICE_DESCR
USE MODD_RAIN_ICE_PARAM
USE MODD_PARAM_ICE
USE MODD_PARAM_C1R3

USE MODI_INI_RAIN_ICE
USE MODI_INI_TIWMX

IMPLICIT NONE
! -----------------------------------------------------------------------
!     DUMMY INTEGER SCALARS
INTEGER, INTENT (IN) :: KULOUT
REAL, INTENT (IN) :: PTSTEP
LOGICAL, INTENT (IN) :: LDWARM
CHARACTER(4), INTENT (IN) :: CMICRO 
CHARACTER(4), INTENT (IN) :: CCSEDIM
INTEGER, INTENT (OUT) :: KSPLITR
LOGICAL, INTENT (IN) :: LDCRIAUTI
REAL, INTENT (IN) :: PCRIAUTI
REAL, INTENT (IN) :: PT0CRIAUTI
REAL, INTENT (IN) :: PCRIAUTC
REAL, INTENT (IN) :: PTSTEP_TS
CHARACTER(4), INTENT (IN) :: CCSNOWRIMING
REAL, INTENT (IN) :: PMRSTEP
INTEGER, INTENT (IN) :: KMAXITER
LOGICAL, INTENT (IN) :: LDFEEDBACKT
LOGICAL, INTENT (IN) :: LDEVLIMIT
LOGICAL, INTENT (IN) :: LDNULLWETG
LOGICAL, INTENT (IN) :: LDWETGPOST
LOGICAL, INTENT (IN) :: LDNULLWETH
LOGICAL, INTENT (IN) :: LDWETHPOST
REAL, INTENT (IN) :: PFRACM90
LOGICAL, INTENT (IN) :: LDCONVHG
CHARACTER(LEN=80), INTENT(IN) :: CCSUBG_RC_RR_ACCR
CHARACTER(LEN=80), INTENT(IN) :: CCSUBG_RR_EVAP
CHARACTER(LEN=80), INTENT(IN) :: CCSUBG_PR_PDF
LOGICAL, INTENT (IN) :: LDCRFLIMIT
CHARACTER(LEN=1), INTENT(IN) :: CCFRAC_ICE_ADJUST
REAL, INTENT (IN) :: PSPLIT_MAXCFL
CHARACTER(LEN=1), INTENT(IN) :: CCFRAC_ICE_SHALLOW_MF
LOGICAL, INTENT (IN) :: LDSEDIM_AFTER
LOGICAL, INTENT (IN) :: LDDEPOSC
REAL, INTENT(IN):: PVDEPOSC
REAL, OPTIONAL, INTENT (IN) :: PFRMIN(40)
LOGICAL, INTENT (IN) :: LDDEPSG 
REAL, INTENT (IN) :: PRDEPSRED, PRDEPGRED

!-----------------------------------------------------------------------
!    LOCAL VARIABLES
REAL :: ZCRI0, ZTCRI0
REAL(KIND=JPRB) :: ZHOOK_HANDLE
! -----------------------------------------------------------------------
!        1. Set implicit default values for MODD_PARAM_ICE

IF (LHOOK) CALL DR_HOOK('AROINI_MICRO',0,ZHOOK_HANDLE)
CALL PARAM_ICE_ASSOCIATE()
!
LWARM=LDWARM
CPRISTINE_ICE='PLAT'
CPRISTINE_ICE_C1R3='PLAT'
CHEVRIMED_ICE_C1R3='GRAU'
CSEDIM=CCSEDIM
CSUBG_RC_RR_ACCR=CCSUBG_RC_RR_ACCR
CSUBG_RR_EVAP=CCSUBG_RR_EVAP
CSUBG_PR_PDF=CCSUBG_PR_PDF
LFEEDBACKT=LDFEEDBACKT ! When .TRUE. feed back on temperature is taken into account
LEVLIMIT=LDEVLIMIT   ! When .TRUE. water vapour pressure is limited by saturation
LNULLWETG=LDNULLWETG  ! When .TRUE. graupel wet growth is activated with null rate (to allow water shedding)
LWETGPOST=LDWETGPOST  ! When .TRUE. graupel wet growth is activated with positive temperature (to allow water shedding)
LNULLWETH=LDNULLWETH  ! Same as LNULLWETG but for hail
LWETHPOST=LDWETHPOST  ! Same as LWETGPOST but for hail
CSNOWRIMING=CCSNOWRIMING ! OLD or M90 for Murakami 1990 formulation
XFRACM90=PFRACM90 ! Fraction used for the Murakami 1990 formulation
NMAXITER=KMAXITER ! Maximum number of iterations for mixing ratio  or time splitting
XMRSTEP=PMRSTEP ! maximum mixing ratio step for mixing ratio splitting
LCONVHG=LDCONVHG ! TRUE to allow the conversion from hail to graupel
LCRFLIMIT=LDCRFLIMIT !True to limit rain contact freezing to possible heat exchange
CFRAC_ICE_ADJUST=CCFRAC_ICE_ADJUST !Choice of solid/liquid partition in adjustements
CFRAC_ICE_SHALLOW_MF=CCFRAC_ICE_SHALLOW_MF !Choice of solid/liquid partition in shallow_mf
XSPLIT_MAXCFL=PSPLIT_MAXCFL
LSEDIM_AFTER=LDSEDIM_AFTER ! sedimentation done before or after microphysics
!
XTSTEP_TS=PTSTEP_TS ! Approximative time step for time-splitting (0 for no time-splitting)
!
!        2. Set implicit default values for MODD_RAIN_ICE_DESCR 
!                     et MODD_RAIN_ICE_PARAM

CALL INI_RAIN_ICE (KULOUT, PTSTEP, 20.,KSPLITR,CMICRO)
CALL INI_TIWMX

IF(PRESENT(PFRMIN))THEN
   XFRMIN = PFRMIN
   WRITE(UNIT=KULOUT,FMT='('' UPDATED VALUES OF XFRMIN FROM NAMPARAR :'')')
   WRITE(UNIT=KULOUT,FMT='('' XFRMIN = '',40E10.3)') XFRMIN
   IF(XFRMIN(16) > 0.) THEN
      CALL INI_SNOW(KULOUT) ! Recalculate snow parameters :  XCCS = XFRMIN(16),XCXS = XFRMIN(17)
   ENDIF
ENDIF

!update values from namparar
LDEPOSC=LDDEPOSC
XVDEPOSC=PVDEPOSC
IF (LDCRIAUTI) THEN

  XCRIAUTI=PCRIAUTI
  XCRIAUTC=PCRIAUTC
  XT0CRIAUTI=PT0CRIAUTI
  !second point to determine 10**(aT+b) law
  ZTCRI0=-40.0
  ZCRI0=1.25E-6
  
  XBCRIAUTI=-( LOG10(XCRIAUTI) - LOG10(ZCRI0)*PT0CRIAUTI/ZTCRI0 )&
                   *ZTCRI0/(XT0CRIAUTI-ZTCRI0)
  XACRIAUTI=(LOG10(ZCRI0)-XBCRIAUTI)/ZTCRI0
  
  !        3. Write NSPLITR,updated CRIAUTI  
      
  WRITE(UNIT=KULOUT,FMT='('' NSPLITR = '',I8.4)')KSPLITR
  WRITE(UNIT=KULOUT,FMT='('' UPDATED VALUES FROM NAMPARAR :'')')
  WRITE(UNIT=KULOUT,FMT='('' LCRIAUTI = '',L5)')LDCRIAUTI
  WRITE(UNIT=KULOUT,FMT='('' XCRIAUTI = '',E13.6)')XCRIAUTI
  WRITE(UNIT=KULOUT,FMT='('' XACRIAUTI = '',E13.6)')XACRIAUTI
  WRITE(UNIT=KULOUT,FMT='('' XBCRIAUTI = '',E13.6)')XBCRIAUTI
  WRITE(UNIT=KULOUT,FMT='('' XT0CRIAUTI = '',E13.6)')XT0CRIAUTI
  WRITE(UNIT=KULOUT,FMT='('' XCRIAUTC = '',E13.6)')XCRIAUTC
  WRITE(UNIT=KULOUT,FMT='('' XVDEPOSC = '',E13.6)')XVDEPOSC
  WRITE(UNIT=KULOUT,FMT='('' LDEPOSC = '',L5)')LDEPOSC
ENDIF

XRDEPSRED=PRDEPSRED
XRDEPGRED=PRDEPGRED

! -----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('AROINI_MICRO',1,ZHOOK_HANDLE)
RETURN
END SUBROUTINE AROINI_MICRO
