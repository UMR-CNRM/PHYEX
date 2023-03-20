!     ######spl
SUBROUTINE AROINI_MICRO_LIMA(KULOUT,KULNAM,PTSTEP,CMICRO,KSPLITR,KSPLITG)

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
!        B. Vie 

!     Modifications.
!     --------------
!        Original : 17-10-09
!     ------------------------------------------------------------------

USE MODD_NSV
USE MODD_LIMA_PRECIP_SCAVENGING_n
USE MODD_PARAM_LIMA_COLD
USE MODD_PARAM_LIMA
USE MODD_PARAM_LIMA_MIXED
USE MODD_PARAM_LIMA_WARM
 
USE MODI_INI_LIMA
USE MODI_INIT_AEROSOL_PROPERTIES

USE MODD_LUNIT, ONLY : ILUOUT

IMPLICIT NONE
! -----------------------------------------------------------------------
!     DUMMY INTEGER SCALARS
INTEGER, INTENT (IN) :: KULOUT
INTEGER, INTENT (IN) :: KULNAM
REAL, INTENT (IN) :: PTSTEP
CHARACTER(4), INTENT (IN) :: CMICRO 
INTEGER, INTENT (OUT) :: KSPLITR
INTEGER, INTENT (OUT) :: KSPLITG
!-----------------------------------------------------------------------
!    LOCAL VARIABLES
REAL :: ZCRI0, ZTCRI0   
INTEGER :: ISV 
REAL(KIND=JPRB) :: ZHOOK_HANDLE
! -----------------------------------------------------------------------
!
#include "namlima.nam.h"
#include "posnam.intfb.h"
!
! -----------------------------------------------------------------------

ILUOUT = KULOUT


! -----------------------------------------------------------------------
! lecture Valeurs par défaut pour les paramètres de la namelist LIMA
!
LPTSPLIT           = .FALSE.
LFEEDBACKT         = .TRUE.
NMAXITER           = 5
XMRSTEP            = 0.
XTSTEP_TS          = 0.
!
NMOM_C             = 2
NMOM_R             = 2
NMOM_I             = 2
NMOM_S             = 1
NMOM_G             = 1
NMOM_H             = 0
!
LNUCL              = .TRUE.
LSEDI              = .TRUE.
LSNOW_T            = .FALSE.
LHHONI             = .FALSE.
LMEYERS            = .FALSE.
LCIBU              = .FALSE.
LRDSF              = .FALSE.
LMURAKAMI          = .FALSE.
NMOD_IFN           = 1
XIFN_CONC(1)       = 1000
LIFN_HOM           = .TRUE.
CIFN_SPECIES       = 'PHILLIPS'
CINT_MIXING        = ''
NMOD_IMM           = 0
NIND_SPECIE        = 1
CPRISTINE_ICE_LIMA = 'PLAT'
CHEVRIMED_ICE_LIMA = 'GRAU'
XALPHAI            = 0.
XNUI               = 0.
XALPHAS            = 0.
XNUS               = 0.
XALPHAG            = 0.
XNUG               = 0.
XFACTNUC_DEP       = 1.
XFACTNUC_CON       = 1.
NPHILLIPS          = 8
!
LACTI              = .TRUE.
LSEDC              = .TRUE.
LDEPOC             = .TRUE.
LACTIT             = .FALSE.
LACTTKE            = .TRUE.
LADJ               = .TRUE.
LSPRO              = .FALSE.
LKHKO              = .FALSE.
LKESSLERAC         = .FALSE.
NMOD_CCN           = 1
XCCN_CONC(1)       = 350.
LCCN_HOM           = .TRUE.
CCCN_MODES         = ''
HINI_CCN           = 'XXX'
HTYPE_CCN          = 'X'
XALPHAC            = 3.
XNUC               = 1.
XALPHAR            = 1.
XNUR               = 2.
XFSOLUB_CCN        = 1.
XACTEMP_CCN        = 280.
XAERDIFF           = 0.
XAERHEIGHT         = 2000.
LSCAV              = .FALSE.
LAERO_MASS         = .FALSE.
! -----------------------------------------------------------------------
! lecture de la namelist LIMA
  CALL POSNAM(KULNAM,'NAMLIMA')
  READ(KULNAM,NAMLIMA)
! -----------------------------------------------------------------------
! initialisation des NSV
  ISV = 1

  NSV_LIMA_BEG = ISV
! Nc
  IF (NMOM_C.GE.2) THEN
     NSV_LIMA_NC = ISV
     ISV = ISV+1
  END IF
! Nr
  IF (NMOM_R.GE.2) THEN
     NSV_LIMA_NR = ISV
     ISV = ISV+1
  END IF
! CCN
  IF (NMOD_CCN .GT. 0) THEN
     NSV_LIMA_CCN_FREE = ISV
     ISV = ISV + NMOD_CCN
     NSV_LIMA_CCN_ACTI = ISV
     ISV = ISV + NMOD_CCN
  END IF
! Scavenging
  IF (LSCAV .AND. LAERO_MASS) THEN
     NSV_LIMA_SCAVMASS = ISV
     ISV = ISV+1
  END IF ! LSCAV
! 
! Ni
  IF (NMOM_I.GE.2) THEN
     NSV_LIMA_NI = ISV
     ISV = ISV+1
  END IF ! LCOLD_LIMA
! Ns
  IF (NMOM_S.GE.2) THEN
     NSV_LIMA_NS = ISV
     ISV = ISV+1
  END IF ! LCOLD_LIMA
! Ng
  IF (NMOM_G.GE.2) THEN
     NSV_LIMA_NG = ISV
     ISV = ISV+1
  END IF ! LCOLD_LIMA
! Nh
  IF (NMOM_H.GE.2) THEN
     NSV_LIMA_NH = ISV
     ISV = ISV+1
  END IF ! LCOLD_LIMA
! IFN
  IF (NMOD_IFN .GT. 0) THEN
     NSV_LIMA_IFN_FREE = ISV
     ISV = ISV + NMOD_IFN
     NSV_LIMA_IFN_NUCL = ISV
     ISV = ISV + NMOD_IFN
  END IF
! IMM
  IF (NMOD_IMM .GT. 0) THEN
     NSV_LIMA_IMM_NUCL = ISV
     ISV = ISV + MAX(1,NMOD_IMM)
  END IF
! Homogeneous freezing of CCN
  IF (NMOM_I.GE.1 .AND. LHHONI) THEN
     NSV_LIMA_HOM_HAZE = ISV
     ISV = ISV + 1
  END IF
! End and total variables
  ISV = ISV - 1
  NSV_LIMA_END = ISV
  NSV_LIMA = NSV_LIMA_END - NSV_LIMA_BEG + 1

NSV=NSV_LIMA

! -----------------------------------------------------------------------
! initialisation de LIMA
CALL INIT_AEROSOL_PROPERTIES
! PDZMIN = 20 comme dans l'appel à INI_RAIN_ICE !
CALL INI_LIMA(PTSTEP, 20., KSPLITR, KSPLITG)

!!$WRITE(UNIT=KULOUT,FMT='(''LIMA SCHEME TUNING  VARIABLES :'')')
!!$WRITE(UNIT=KULOUT,FMT='('' LCOLD_LIMA = '',L5)') LCOLD_LIMA
!!$WRITE(UNIT=KULOUT,FMT='('' LNUCL_LIMA = '',L5)') LNUCL_LIMA
!!$WRITE(UNIT=KULOUT,FMT='('' LSEDI_LIMA = '',L5)') LSEDI_LIMA
!!$WRITE(UNIT=KULOUT,FMT='('' LSNOW_LIMA = '',L5)') LSNOW_LIMA
!!$WRITE(UNIT=KULOUT,FMT='('' LHAIL_LIMA = '',L5)') LHAIL_LIMA
!!$WRITE(UNIT=KULOUT,FMT='('' LHHONI_LIMA = '',L5)') LHHONI_LIMA
!!$WRITE(UNIT=KULOUT,FMT='('' LMEYERS_LIMA = '',L5)') LMEYERS_LIMA
!!$WRITE(UNIT=KULOUT,FMT='('' LIFN_HOM = '',L5)') LIFN_HOM
!!$WRITE(UNIT=KULOUT,FMT='('' LWARM_LIMA = '',L5)') LWARM_LIMA
!!$WRITE(UNIT=KULOUT,FMT='('' LACTI_LIMA = '',L5)') LACTI_LIMA
!!$WRITE(UNIT=KULOUT,FMT='('' LRAIN_LIMA = '',L5)') LRAIN_LIMA
!!$WRITE(UNIT=KULOUT,FMT='('' LSEDC_LIMA = '',L5)') LSEDC_LIMA
!!$WRITE(UNIT=KULOUT,FMT='('' LACTIT_LIMA = '',L5)') LACTIT_LIMA
!!$WRITE(UNIT=KULOUT,FMT='('' LCCN_HOM = '',L5)') LCCN_HOM
!!$WRITE(UNIT=KULOUT,FMT='('' LSCAV = '',L5)') LSCAV
!!$WRITE(UNIT=KULOUT,FMT='('' LAERO_MASS = '',L5)') LAERO_MASS
!!$WRITE(UNIT=KULOUT,FMT='('' CIFN_SPECIES = '',A8,''CINT_MIXING = '',A8)')&
!!$&CIFN_SPECIES,CINT_MIXING
!!$WRITE(UNIT=KULOUT,FMT='('' CPRISTINE_ICE_LIMA = '',A4,''CHEVRIMED_ICE_LIMA = '',A4)')&
!!$&CPRISTINE_ICE_LIMA, CHEVRIMED_ICE_LIMA
!!$WRITE(UNIT=KULOUT,FMT='('' CCCN_MODES = '',A8)')CCCN_MODES
!!$WRITE(UNIT=KULOUT,FMT='('' HINI_CCN = '',A3,''HTYPE_CCN = '',A1)')&
!!$&HINI_CCN,HTYPE_CCN
!!$WRITE(UNIT=KULOUT,FMT='('' NMOD_IFN = '',I5)') NMOD_IFN
!!$WRITE(UNIT=KULOUT,FMT='('' NMOD_IMM = '',I5)') NMOD_IMM
!!$WRITE(UNIT=KULOUT,FMT='('' NIND_SPECIE = '',I5)') NIND_SPECIE
!!$WRITE(UNIT=KULOUT,FMT='('' NPHILLIPS = '',I5)') NPHILLIPS
!!$WRITE(UNIT=KULOUT,FMT='('' NMOD_CCN = '',I5)') NMOD_CCN
!!$WRITE(UNIT=KULOUT,FMT='('' XIFN_CONC = '',f6.2)') XIFN_CONC
!!$WRITE(UNIT=KULOUT,FMT='('' XALPHAI = '',f6.2)') XALPHAI
!!$WRITE(UNIT=KULOUT,FMT='('' XNUI = '',f6.2)') XNUI
!!$WRITE(UNIT=KULOUT,FMT='('' XALPHAS = '',f6.2)') XALPHAS
!!$WRITE(UNIT=KULOUT,FMT='('' XNUS = '',f6.2)') XNUS
!!$WRITE(UNIT=KULOUT,FMT='('' XALPHAG = '',f6.2)') XALPHAG
!!$WRITE(UNIT=KULOUT,FMT='('' XNUG = '',f6.2)') XNUG
!!$WRITE(UNIT=KULOUT,FMT='('' XCCN_CONC = '',f6.2)') XCCN_CONC
!!$WRITE(UNIT=KULOUT,FMT='('' XALPHAC = '',f6.2)') XALPHAC
!!$WRITE(UNIT=KULOUT,FMT='('' XNUC = '',f6.2)') XNUC
!!$WRITE(UNIT=KULOUT,FMT='('' XALPHAR = '',f6.2)') XALPHAR
!!$WRITE(UNIT=KULOUT,FMT='('' XNUR = '',f6.2)') XNUR
!!$WRITE(UNIT=KULOUT,FMT='('' XFSOLUB_CCN = '',f6.2)') XFSOLUB_CCN
!!$WRITE(UNIT=KULOUT,FMT='('' XACTEMP_CCN = '',f6.2)') XACTEMP_CCN
!!$WRITE(UNIT=KULOUT,FMT='('' XAERDIFF = '',f6.2)') XAERDIFF
!!$WRITE(UNIT=KULOUT,FMT='('' XAERHEIGHT = '',f6.2)') XAERHEIGHT



RETURN
END SUBROUTINE AROINI_MICRO_LIMA
