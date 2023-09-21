!MNH_LIC Copyright 2001-2023 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ###################
      MODULE MODI_INI_NSV
!     ###################
INTERFACE 
!
  SUBROUTINE INI_NSV(KMI)
  INTEGER, INTENT(IN)            :: KMI ! model index
  END SUBROUTINE INI_NSV
!
END INTERFACE
!
END MODULE MODI_INI_NSV
!
!
!     ###########################
      SUBROUTINE INI_NSV(KMI)
!     ###########################
!
!!****   *INI_NSV* - compute NSV_* values and indices for model KMI
!!
!!    PURPOSE
!!    -------
!     
!
!     
!!**  METHOD
!!    ------
!!
!!    This routine is called from any routine which stores values in 
!!    the first model module (for example READ_EXSEG).
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_NSV     : contains NSV_A array variable
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!      D. Gazen              * LA *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original   01/02/01
!!      Modification   29/11/02  (Pinty)  add SV for C3R5 and ELEC
!!      Modification   01/2004   (Masson) add scalar names
!!      Modification   03/2006   (O.Geoffroy) add KHKO scheme
!!      Modification   04/2007   (Leriche) add SV for aqueous chemistry
!!      M. Chong       26/01/10   Add Small ions
!!      Modification   07/2010   (Leriche) add SV for ice chemistry
!!      X.Pialat & J.Escobar 11/2012 remove deprecated line NSV_A(KMI) = ISV
!!      Modification   15/02/12  (Pialat/Tulet) Add SV for ForeFire scalars
!!                     03/2013   (C.Lac) add supersaturation as 
!!                               the 4th C2R2 scalar variable
!!       J.escobar     04/08/2015 suit Pb with writ_lfin JSA increment , modif in ini_nsv to have good order initialization
!!      Modification    01/2016  (JP Pinty) Add LIMA and LUSECHEM condition
!!      Modification    07/2017  (V. Vionnet) Add blowing snow condition
!  P. Wautelet 09/03/2021: move some chemistry initializations to ini_nsv
!  P. Wautelet 10/03/2021: move scalar variable name initializations to ini_nsv
!  P. Wautelet 30/03/2021: move NINDICE_CCN_IMM and NIMM initializations from init_aerosol_properties to ini_nsv
!  B. Vie         06/2021: add prognostic supersaturation for LIMA
!  P. Wautelet 26/11/2021: initialize TSVLIST_A
!  A. Costes      12/2021: smoke tracer for fire model
!  P. Wautelet 14/01/2022: add CSV_CHEM_LIST(_A) to store the list of all chemical variables
!                          + NSV_CHEM_LIST(_A) the size of the list
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_BLOWSNOW,        ONLY: CSNOWNAMES, LBLOWSNOW, NBLOWSNOW3D, YPSNOW_INI
USE MODD_CH_AEROSOL
! USE MODD_CH_AEROSOL,      ONLY: CAERONAMES, CDEAERNAMES, JPMODE, LAERINIT, LDEPOS_AER, LORILAM, &
!                                 LVARSIGI, LVARSIGJ, NCARB, NM6_AER, NSOA, NSP
USE MODD_CH_M9_n,         ONLY: CICNAMES, CNAMES, NEQ, NEQAQ
USE MODD_CH_MNHC_n,       ONLY: LCH_PH, LUSECHEM, LUSECHAQ, LUSECHIC, CCH_SCHEME, LCH_CONV_LINOX
USE MODD_CONDSAMP,        ONLY: LCONDSAMP, NCONDSAMP
USE MODD_CONF,            ONLY: LLG, CPROGRAM, NVERB
USE MODD_CST,             ONLY: XMNH_TINY
USE MODD_DIAG_FLAG,       ONLY: LCHEMDIAG, LCHAQDIAG
USE MODD_DUST,            ONLY: CDEDSTNAMES, CDUSTNAMES, JPDUSTORDER, LDEPOS_DST, LDSTINIT, LDSTPRES, LDUST, &
                                LRGFIX_DST, LVARSIG, NMODE_DST, YPDEDST_INI, YPDUST_INI
USE MODD_DYN_n,           ONLY: LHORELAX_SV,LHORELAX_SVC2R2,LHORELAX_SVC1R3,   &
                                LHORELAX_SVFIRE, LHORELAX_SVLIMA,              &
                                LHORELAX_SVELEC,LHORELAX_SVCHEM,LHORELAX_SVLG, &
                                LHORELAX_SVDST,LHORELAX_SVAER, LHORELAX_SVSLT, &
                                LHORELAX_SVPP,LHORELAX_SVCS, LHORELAX_SVCHIC,  &
                                LHORELAX_SVSNW
#ifdef MNH_FOREFIRE
USE MODD_DYN_n,           ONLY: LHORELAX_SVFF
#endif
USE MODD_ELEC_DESCR,      ONLY: LLNOX_EXPLICIT
USE MODD_ELEC_DESCR,      ONLY: CELECNAMES
USE MODD_FIELD,           ONLY: TFIELDMETADATA, TYPEREAL
USE MODD_FIRE_n
#ifdef MNH_FOREFIRE
USE MODD_FOREFIRE
#endif
USE MODD_ICE_C1R3_DESCR,  ONLY: C1R3NAMES
USE MODD_LG,              ONLY: CLGNAMES, XLG1MIN, XLG2MIN, XLG3MIN
USE MODD_LUNIT_n,         ONLY: TLUOUT
USE MODD_NSV
USE MODD_PARAM_C2R2,      ONLY: LSUPSAT
USE MODD_PARAMETERS,      ONLY: NCOMMENTLGTMAX, NLONGNAMELGTMAX, NUNITLGTMAX
USE MODD_PARAM_LIMA,      ONLY: NINDICE_CCN_IMM, NIMM, NMOD_CCN, NMOD_IFN, NMOD_IMM, PARAM_LIMA_ALLOCATE, PARAM_LIMA_DEALLOCATE
USE MODD_PARAM_LIMA_COLD, ONLY: CLIMA_COLD_NAMES
USE MODD_PARAM_LIMA_WARM, ONLY: CAERO_MASS, CLIMA_WARM_NAMES
USE MODD_PARAM_n,         ONLY: CCLOUD, CELEC
USE MODD_PASPOL,          ONLY: LPASPOL, NRELEASE
USE MODD_PREP_REAL,       ONLY: XT_LS
USE MODD_RAIN_C2R2_DESCR, ONLY: C2R2NAMES
USE MODD_SALT,            ONLY: CSALTNAMES, CDESLTNAMES, JPSALTORDER, &
                                LRGFIX_SLT, LSALT, LSLTINIT, LSLTPRES, LDEPOS_SLT, LVARSIG_SLT, NMODE_SLT, YPDESLT_INI, YPSALT_INI

USE MODE_MSG
USE MODE_LIMA_UPDATE_NSV, ONLY: LIMA_UPDATE_NSV

USE MODI_CH_AER_INIT_SOA,  ONLY: CH_AER_INIT_SOA
USE MODI_CH_INIT_SCHEME_n, ONLY: CH_INIT_SCHEME_n
USE MODI_UPDATE_NSV,       ONLY: UPDATE_NSV
!
IMPLICIT NONE 
!
!-------------------------------------------------------------------------------
!
!*       0.1   Declarations of arguments
!
INTEGER, INTENT(IN)             :: KMI ! model index
!
!*       0.2   Declarations of local variables
!
CHARACTER(LEN=2) :: YNUM2
CHARACTER(LEN=3) :: YNUM3
CHARACTER(LEN=NCOMMENTLGTMAX) :: YCOMMENT
CHARACTER(LEN=NUNITLGTMAX)    :: YUNITS
CHARACTER(LEN=NLONGNAMELGTMAX), DIMENSION(:), ALLOCATABLE :: YAEROLONGNAMES
CHARACTER(LEN=NLONGNAMELGTMAX), DIMENSION(:), ALLOCATABLE :: YDUSTLONGNAMES
CHARACTER(LEN=NLONGNAMELGTMAX), DIMENSION(:), ALLOCATABLE :: YSALTLONGNAMES
INTEGER :: ILUOUT
INTEGER :: ICHIDX ! Index for position in CSV_CHEM_LIST_A array
INTEGER :: ISV ! total number of scalar variables
INTEGER :: IMODEIDX
INTEGER :: JAER
INTEGER :: JI, JJ, JSV
INTEGER :: JMODE, JMOM, JSV_NAME
INTEGER :: INMOMENTS_DST, INMOMENTS_SLT !Number of moments for dust or salt
!
!-------------------------------------------------------------------------------
!

!Associate the pointers
CALL NSV_ASSOCIATE
!
LINI_NSV(KMI) = .TRUE.

ILUOUT = TLUOUT%NLU

ICHIDX = 0
NSV_CHEM_LIST_A(KMI) = 0
!
! Users scalar variables are first considered
!
NSV_USER_A(KMI) = NSV_USER
ISV = NSV_USER
!
! scalar variables used in microphysical schemes C2R2,KHKO and C3R5
!
IF (CCLOUD == 'C2R2' .OR. CCLOUD == 'C3R5' .OR. CCLOUD == 'KHKO' ) THEN
  IF ((CCLOUD == 'C2R2' .AND. LSUPSAT) .OR. (CCLOUD == 'KHKO'.AND. LSUPSAT)) THEN
   ! 4th scalar field = supersaturation
    NSV_C2R2_A(KMI)    = 4
  ELSE
    NSV_C2R2_A(KMI)    = 3
  END IF
  NSV_C2R2BEG_A(KMI) = ISV+1
  NSV_C2R2END_A(KMI) = ISV+NSV_C2R2_A(KMI)
  ISV                = NSV_C2R2END_A(KMI)
  IF (CCLOUD == 'C3R5') THEN  ! the SVs for C2R2 and C1R3 must be contiguous
    NSV_C1R3_A(KMI)    = 2
    NSV_C1R3BEG_A(KMI) = ISV+1
    NSV_C1R3END_A(KMI) = ISV+NSV_C1R3_A(KMI)
    ISV                = NSV_C1R3END_A(KMI)
  ELSE
    NSV_C1R3_A(KMI)    = 0
  ! force First index to be superior to last index
  ! in order to create a null section
    NSV_C1R3BEG_A(KMI) = 1
    NSV_C1R3END_A(KMI) = 0
  END IF
ELSE
  NSV_C2R2_A(KMI)    = 0
  NSV_C1R3_A(KMI)    = 0
! force First index to be superior to last index
! in order to create a null section
  NSV_C2R2BEG_A(KMI) = 1
  NSV_C2R2END_A(KMI) = 0
  NSV_C1R3BEG_A(KMI) = 1
  NSV_C1R3END_A(KMI) = 0
END IF
!
! scalar variables used in the LIMA microphysical scheme
!
CALL LIMA_UPDATE_NSV(LDINIT=.TRUE., KMI=KMI, KSV=ISV, CDCLOUD=CCLOUD, LDUPDATE=.FALSE.)
IF (CCLOUD == 'LIMA' ) THEN

  IF ( NMOD_IFN > 0 ) THEN
    IF ( .NOT. ASSOCIATED( NIMM ) ) CALL PARAM_LIMA_ALLOCATE('NIMM', NMOD_CCN)
    NIMM(:) = 0
    IF ( ASSOCIATED( NINDICE_CCN_IMM ) ) CALL PARAM_LIMA_DEALLOCATE('NINDICE_CCN_IMM')
    CALL PARAM_LIMA_ALLOCATE('NINDICE_CCN_IMM', MAX( 1, NMOD_IMM ))
    IF (NMOD_IMM > 0 ) THEN
      DO JI = 0, NMOD_IMM - 1
        NIMM(NMOD_CCN - JI) = 1
        NINDICE_CCN_IMM(NMOD_IMM - JI) = NMOD_CCN - JI
      END DO
!     ELSE IF (NMOD_IMM == 0) THEN ! PNIS exists but is 0 for the call to resolved_cloud
!       NMOD_IMM = 1
!       NINDICE_CCN_IMM(1) = 0
    END IF
  END IF
END IF ! CCLOUD = LIMA
!
!
!  Add one scalar for negative ion
!   First variable: positive ion (NSV_ELECBEG_A index number)
!   Last  --------: negative ion (NSV_ELECEND_A index number)
! Correspondence for ICE3:
! Relative index    1       2        3       4      5      6       7
! Charge for     ion+     cloud    rain     ice   snow  graupel  ion-
!
! Correspondence for ICE4:
! Relative index    1       2        3       4      5      6       7       8
! Charge for     ion+     cloud    rain     ice   snow  graupel   hail   ion-
!
IF (CELEC /= 'NONE') THEN
  IF (CCLOUD == 'ICE3') THEN
    NSV_ELEC_A(KMI)   = 7 
    NSV_ELECBEG_A(KMI)= ISV+1
    NSV_ELECEND_A(KMI)= ISV+NSV_ELEC_A(KMI)
    ISV               = NSV_ELECEND_A(KMI)
    CELECNAMES(7) = CELECNAMES(8) 
  ELSE IF (CCLOUD == 'ICE4') THEN
    NSV_ELEC_A(KMI)   = 8 
    NSV_ELECBEG_A(KMI)= ISV+1
    NSV_ELECEND_A(KMI)= ISV+NSV_ELEC_A(KMI)
    ISV               = NSV_ELECEND_A(KMI)
  END IF
ELSE
  NSV_ELEC_A(KMI)    = 0
! force First index to be superior to last index
! in order to create a null section
  NSV_ELECBEG_A(KMI) = 1
  NSV_ELECEND_A(KMI) = 0
END IF
!
! scalar variables used as lagragian variables
!
IF (LLG) THEN
  NSV_LG_A(KMI)     = 3
  NSV_LGBEG_A(KMI)  = ISV+1
  NSV_LGEND_A(KMI)  = ISV+NSV_LG_A(KMI)
  ISV               = NSV_LGEND_A(KMI)
ELSE
  NSV_LG_A(KMI)     = 0
! force First index to be superior to last index
! in order to create a null section
  NSV_LGBEG_A(KMI)  = 1
  NSV_LGEND_A(KMI)  = 0
END IF
!
! scalar variables used as LiNOX passive tracer
!
! In case without chemistry
IF (LPASPOL) THEN
  NSV_PP_A(KMI)   = NRELEASE
  NSV_PPBEG_A(KMI)= ISV+1
  NSV_PPEND_A(KMI)= ISV+NSV_PP_A(KMI)
  ISV               = NSV_PPEND_A(KMI)
ELSE
  NSV_PP_A(KMI)   = 0
! force First index to be superior to last index
! in order to create a null section
  NSV_PPBEG_A(KMI)= 1
  NSV_PPEND_A(KMI)= 0
END IF
!
#ifdef MNH_FOREFIRE
! ForeFire tracers
IF (LFOREFIRE .AND. NFFSCALARS .GT. 0) THEN
  NSV_FF_A(KMI)    = NFFSCALARS
  NSV_FFBEG_A(KMI) = ISV+1
  NSV_FFEND_A(KMI) = ISV+NSV_FF_A(KMI)
  ISV              = NSV_FFEND_A(KMI)
ELSE
  NSV_FF_A(KMI)   = 0
! force First index to be superior to last index
! in order to create a null section
  NSV_FFBEG_A(KMI)= 1
  NSV_FFEND_A(KMI)= 0
END IF
#endif
! Blaze tracers
IF (LBLAZE .AND. NNBSMOKETRACER .GT. 0) THEN
  NSV_FIRE_A(KMI)    = NNBSMOKETRACER
  NSV_FIREBEG_A(KMI) = ISV+1
  NSV_FIREEND_A(KMI) = ISV+NSV_FIRE_A(KMI)
  ISV              = NSV_FIREEND_A(KMI)
ELSE
  NSV_FIRE_A(KMI)   = 0
! force First index to be superior to last index
! in order to create a null section
  NSV_FIREBEG_A(KMI)= 1
  NSV_FIREEND_A(KMI)= 0
END IF
!
! Conditional sampling variables  
IF (LCONDSAMP) THEN
  NSV_CS_A(KMI)   = NCONDSAMP
  NSV_CSBEG_A(KMI)= ISV+1
  NSV_CSEND_A(KMI)= ISV+NSV_CS_A(KMI)
  ISV               = NSV_CSEND_A(KMI)
ELSE
  NSV_CS_A(KMI)   = 0
! force First index to be superior to last index
! in order to create a null section
  NSV_CSBEG_A(KMI)= 1
  NSV_CSEND_A(KMI)= 0
END IF
!
! scalar variables used in chemical core system
!
IF (LUSECHEM) THEN
  CALL CH_INIT_SCHEME_n(KMI,LUSECHAQ,LUSECHIC,LCH_PH,ILUOUT,NVERB)
  IF (LORILAM) CALL CH_AER_INIT_SOA(ILUOUT, NVERB)
END IF

IF (LUSECHEM .AND.(NEQ .GT. 0)) THEN
  NSV_CHEM_A(KMI)   = NEQ
  NSV_CHEMBEG_A(KMI)= ISV+1
  NSV_CHEMEND_A(KMI)= ISV+NSV_CHEM_A(KMI)
  ISV               = NSV_CHEMEND_A(KMI)
  NSV_CHEM_LIST_A(KMI) = NSV_CHEM_LIST_A(KMI) + NSV_CHEM_A(KMI)
ELSE
  NSV_CHEM_A(KMI)   = 0
! force First index to be superior to last index
! in order to create a null section
  NSV_CHEMBEG_A(KMI)= 1
  NSV_CHEMEND_A(KMI)= 0
END IF
!
! aqueous chemistry (part of the "chem" variables)       
!                                                        
IF ((LUSECHAQ .OR. LCHAQDIAG).AND.(NEQ .GT. 0)) THEN     
  NSV_CHGS_A(KMI) = NEQ-NEQAQ                            
  NSV_CHGSBEG_A(KMI)= NSV_CHEMBEG_A(KMI)                 
  NSV_CHGSEND_A(KMI)= NSV_CHEMBEG_A(KMI)+(NEQ-NEQAQ)-1   
  NSV_CHAC_A(KMI) = NEQAQ                                
  NSV_CHACBEG_A(KMI)= NSV_CHGSEND_A(KMI)+1               
  NSV_CHACEND_A(KMI)= NSV_CHEMEND_A(KMI)                 
!  ice phase chemistry
  IF (LUSECHIC) THEN
    NSV_CHIC_A(KMI) = NEQAQ/2. -1.
    NSV_CHICBEG_A(KMI)= ISV+1
    NSV_CHICEND_A(KMI)= ISV+NSV_CHIC_A(KMI)
    ISV               = NSV_CHICEND_A(KMI)
    NSV_CHEM_LIST_A(KMI) = NSV_CHEM_LIST_A(KMI) + NSV_CHIC_A(KMI)
  ELSE
    NSV_CHIC_A(KMI) = 0
    NSV_CHICBEG_A(KMI)= 1
    NSV_CHICEND_A(KMI)= 0
  ENDIF
ELSE                                                     
  IF (NEQ .GT. 0) THEN
    NSV_CHGS_A(KMI) = NEQ-NEQAQ                            
    NSV_CHGSBEG_A(KMI)= NSV_CHEMBEG_A(KMI)                 
    NSV_CHGSEND_A(KMI)= NSV_CHEMBEG_A(KMI)+(NEQ-NEQAQ)-1   
    NSV_CHAC_A(KMI) = 0                                    
    NSV_CHACBEG_A(KMI)= 1                                  
    NSV_CHACEND_A(KMI)= 0                                  
    NSV_CHIC_A(KMI) = 0
    NSV_CHICBEG_A(KMI)= 1
    NSV_CHICEND_A(KMI)= 0
  ELSE
    NSV_CHGS_A(KMI) = 0
    NSV_CHGSBEG_A(KMI)= 1
    NSV_CHGSEND_A(KMI)= 0
    NSV_CHAC_A(KMI) = 0
    NSV_CHACBEG_A(KMI)= 1
    NSV_CHACEND_A(KMI)= 0   
    NSV_CHIC_A(KMI) = 0
    NSV_CHICBEG_A(KMI)= 1
    NSV_CHICEND_A(KMI)= 0    
  ENDIF
END IF
! aerosol variables
IF (LORILAM.AND.(NEQ .GT. 0)) THEN
  NM6_AER = 0
  IF (LVARSIGI) NM6_AER = 1
  IF (LVARSIGJ) NM6_AER = NM6_AER + 1
  NSV_AER_A(KMI)   = (NSP+NCARB+NSOA+1)*JPMODE + NM6_AER
  NSV_AERBEG_A(KMI)= ISV+1
  NSV_AEREND_A(KMI)= ISV+NSV_AER_A(KMI)
  ISV              = NSV_AEREND_A(KMI)
  NSV_CHEM_LIST_A(KMI) = NSV_CHEM_LIST_A(KMI) + NSV_AER_A(KMI)

  ALLOCATE( YAEROLONGNAMES(NSV_AER_A(KMI)) )
ELSE
  NSV_AER_A(KMI)   = 0
! force First index to be superior to last index
! in order to create a null section
  NSV_AERBEG_A(KMI)= 1
  NSV_AEREND_A(KMI)= 0
END IF
IF (LORILAM .AND. LDEPOS_AER(KMI)) THEN
  NSV_AERDEP_A(KMI)   = JPMODE*2
  NSV_AERDEPBEG_A(KMI)= ISV+1
  NSV_AERDEPEND_A(KMI)= ISV+NSV_AERDEP_A(KMI)
  ISV                  = NSV_AERDEPEND_A(KMI)       
  NSV_CHEM_LIST_A(KMI) = NSV_CHEM_LIST_A(KMI) + NSV_AERDEP_A(KMI)
ELSE
  NSV_AERDEP_A(KMI)   = 0
! force First index to be superior to last index
! in order to create a null section
  NSV_AERDEPBEG_A(KMI)= 1
  NSV_AERDEPEND_A(KMI)= 0       
! force First index to be superior to last index
! in order to create a null section
END IF
!
! scalar variables used in dust model
!
IF (LDUST) THEN
  IF (ALLOCATED(XT_LS).AND. .NOT.(LDSTPRES)) LDSTINIT=.TRUE.
  IF (CPROGRAM == 'IDEAL ') LVARSIG = .TRUE.
  IF ((CPROGRAM == 'REAL  ').AND.LDSTINIT) LVARSIG = .TRUE.
  !Determine number of moments
  IF ( LRGFIX_DST ) THEN
    INMOMENTS_DST = 1
    IF ( LVARSIG ) CALL Print_msg( NVERB_WARNING, 'GEN', 'INI_NSV', 'LVARSIG forced to FALSE because LRGFIX_DST is TRUE' )
    LVARSIG = .FALSE.
  ELSE IF ( LVARSIG ) THEN
    INMOMENTS_DST = 3
  ELSE
    INMOMENTS_DST = 2
  END IF
  !Number of entries = number of moments multiplied by number of modes
  NSV_DST_A(KMI) = NMODE_DST * INMOMENTS_DST
  NSV_DSTBEG_A(KMI)= ISV+1
  NSV_DSTEND_A(KMI)= ISV+NSV_DST_A(KMI)
  ISV              = NSV_DSTEND_A(KMI)
  NSV_CHEM_LIST_A(KMI) = NSV_CHEM_LIST_A(KMI) + NSV_DST_A(KMI)
ELSE
  NSV_DST_A(KMI)   = 0
! force First index to be superior to last index
! in order to create a null section
  NSV_DSTBEG_A(KMI)= 1
  NSV_DSTEND_A(KMI)= 0
END IF
IF ( LDUST .AND. LDEPOS_DST(KMI) ) THEN
  NSV_DSTDEP_A(KMI)   = NMODE_DST*2
  NSV_DSTDEPBEG_A(KMI)= ISV+1
  NSV_DSTDEPEND_A(KMI)= ISV+NSV_DSTDEP_A(KMI)
  ISV                  = NSV_DSTDEPEND_A(KMI)       
  NSV_CHEM_LIST_A(KMI) = NSV_CHEM_LIST_A(KMI) + NSV_DSTDEP_A(KMI)
ELSE
  NSV_DSTDEP_A(KMI)   = 0
! force First index to be superior to last index
! in order to create a null section
  NSV_DSTDEPBEG_A(KMI)= 1
  NSV_DSTDEPEND_A(KMI)= 0       
! force First index to be superior to last index
! in order to create a null section

 END IF
! scalar variables used in sea salt model
!
IF (LSALT) THEN
  IF (ALLOCATED(XT_LS).AND. .NOT.(LSLTPRES)) LSLTINIT=.TRUE.
  IF (CPROGRAM == 'IDEAL ') LVARSIG_SLT = .TRUE.
  IF ((CPROGRAM == 'REAL  ').AND. LSLTINIT ) LVARSIG_SLT = .TRUE.
  !Determine number of moments
  IF ( LRGFIX_SLT ) THEN
    INMOMENTS_SLT = 1
    IF ( LVARSIG_SLT ) CALL Print_msg( NVERB_WARNING, 'GEN', 'INI_NSV', 'LVARSIG_SLT forced to FALSE because LRGFIX_SLT is TRUE' )
    LVARSIG_SLT = .FALSE.
  ELSE IF ( LVARSIG_SLT ) THEN
    INMOMENTS_SLT = 3
  ELSE
    INMOMENTS_SLT = 2
  END IF
  !Number of entries = number of moments multiplied by number of modes
  NSV_SLT_A(KMI) = NMODE_SLT * INMOMENTS_SLT
  NSV_SLTBEG_A(KMI)= ISV+1
  NSV_SLTEND_A(KMI)= ISV+NSV_SLT_A(KMI)
  ISV              = NSV_SLTEND_A(KMI)
  NSV_CHEM_LIST_A(KMI) = NSV_CHEM_LIST_A(KMI) + NSV_SLT_A(KMI)
ELSE
  NSV_SLT_A(KMI)   = 0
! force First index to be superior to last index
! in order to create a null section
  NSV_SLTBEG_A(KMI)= 1
  NSV_SLTEND_A(KMI)= 0
END IF
IF ( LSALT .AND. LDEPOS_SLT(KMI) ) THEN
  NSV_SLTDEP_A(KMI)   = NMODE_SLT*2
  NSV_SLTDEPBEG_A(KMI)= ISV+1
  NSV_SLTDEPEND_A(KMI)= ISV+NSV_SLTDEP_A(KMI)
  ISV                  = NSV_SLTDEPEND_A(KMI)       
  NSV_CHEM_LIST_A(KMI) = NSV_CHEM_LIST_A(KMI) + NSV_SLTDEP_A(KMI)
ELSE
  NSV_SLTDEP_A(KMI)   = 0
! force First index to be superior to last index
! in order to create a null section
  NSV_SLTDEPBEG_A(KMI)= 1
  NSV_SLTDEPEND_A(KMI)= 0       
! force First index to be superior to last index
! in order to create a null section
END IF
!
! scalar variables used in blowing snow model
!
IF (LBLOWSNOW) THEN
  NSV_SNW_A(KMI)   = NBLOWSNOW3D
  NSV_SNWBEG_A(KMI)= ISV+1
  NSV_SNWEND_A(KMI)= ISV+NSV_SNW_A(KMI)
  ISV              = NSV_SNWEND_A(KMI)
ELSE
  NSV_SNW_A(KMI)   = 0
! force First index to be superior to last index
! in order to create a null section
  NSV_SNWBEG_A(KMI)= 1
  NSV_SNWEND_A(KMI)= 0
END IF
!
! scalar variables used as LiNOX passive tracer
!
! In case without chemistry
IF (.NOT.(LUSECHEM.OR.LCHEMDIAG) .AND. (LCH_CONV_LINOX.OR.LLNOX_EXPLICIT)) THEN
  NSV_LNOX_A(KMI)   = 1
  NSV_LNOXBEG_A(KMI)= ISV+1
  NSV_LNOXEND_A(KMI)= ISV+NSV_LNOX_A(KMI)
  ISV               = NSV_LNOXEND_A(KMI)
  NSV_CHEM_LIST_A(KMI) = NSV_CHEM_LIST_A(KMI) + NSV_LNOX_A(KMI)
ELSE
  NSV_LNOX_A(KMI)   = 0
! force First index to be superior to last index
! in order to create a null section
  NSV_LNOXBEG_A(KMI)= 1
  NSV_LNOXEND_A(KMI)= 0
END IF
!
! Final number of NSV variables
!
NSV_A(KMI) = ISV
!
!
!*        Update LHORELAX_SV,CGETSVM,CGETSVT for NON USER SV 
!
! C2R2  or KHKO SV case
!*BUG*JPC*MAR2006
! IF (CCLOUD == 'C2R2'  .OR. CCLOUD == 'KHKO' ) &
IF (CCLOUD == 'C2R2' .OR. CCLOUD == 'C3R5' .OR.  CCLOUD == 'KHKO' ) &
!*BUG*JPC*MAR2006
LHORELAX_SV(NSV_C2R2BEG_A(KMI):NSV_C2R2END_A(KMI))=LHORELAX_SVC2R2
! C3R5 SV case
IF (CCLOUD == 'C3R5') &
LHORELAX_SV(NSV_C1R3BEG_A(KMI):NSV_C1R3END_A(KMI))=LHORELAX_SVC1R3
! LIMA SV case
IF (CCLOUD == 'LIMA') &
LHORELAX_SV(NSV_LIMA_BEG_A(KMI):NSV_LIMA_END_A(KMI))=LHORELAX_SVLIMA
! Electrical SV case
IF (CELEC /= 'NONE') &
LHORELAX_SV(NSV_ELECBEG_A(KMI):NSV_ELECEND_A(KMI))=LHORELAX_SVELEC
! Chemical SV case
IF (LUSECHEM .OR. LCHEMDIAG) &
LHORELAX_SV(NSV_CHEMBEG_A(KMI):NSV_CHEMEND_A(KMI))=LHORELAX_SVCHEM
! Ice phase Chemical SV case
IF (LUSECHIC) &
LHORELAX_SV(NSV_CHICBEG_A(KMI):NSV_CHICEND_A(KMI))=LHORELAX_SVCHIC
! LINOX SV case
IF (.NOT.(LUSECHEM .OR. LCHEMDIAG) .AND. LCH_CONV_LINOX) &
LHORELAX_SV(NSV_LNOXBEG_A(KMI):NSV_LNOXEND_A(KMI))=LHORELAX_SVCHEM
! Dust SV case
IF (LDUST) &
LHORELAX_SV(NSV_DSTBEG_A(KMI):NSV_DSTEND_A(KMI))=LHORELAX_SVDST
! Sea Salt SV case
IF (LSALT) &
LHORELAX_SV(NSV_SLTBEG_A(KMI):NSV_SLTEND_A(KMI))=LHORELAX_SVSLT
! Aerosols SV case
IF (LORILAM) &
LHORELAX_SV(NSV_AERBEG_A(KMI):NSV_AEREND_A(KMI))=LHORELAX_SVAER
! Lagrangian variables
IF (LLG) &
LHORELAX_SV(NSV_LGBEG_A(KMI):NSV_LGEND_A(KMI))=LHORELAX_SVLG
! Passive pollutants  
IF (LPASPOL) &
LHORELAX_SV(NSV_PPBEG_A(KMI):NSV_PPEND_A(KMI))=LHORELAX_SVPP
#ifdef MNH_FOREFIRE
! Fire pollutants
IF (LFOREFIRE) &
LHORELAX_SV(NSV_FFBEG_A(KMI):NSV_FFEND_A(KMI))=LHORELAX_SVFF
#endif
! Blaze Fire pollutants
IF (LBLAZE) &
LHORELAX_SV(NSV_FIREBEG_A(KMI):NSV_FIREEND_A(KMI))=LHORELAX_SVFIRE
! Conditional sampling
IF (LCONDSAMP) &
LHORELAX_SV(NSV_CSBEG_A(KMI):NSV_CSEND_A(KMI))=LHORELAX_SVCS
! Blowing snow case
IF (LBLOWSNOW) &
LHORELAX_SV(NSV_SNWBEG_A(KMI):NSV_SNWEND_A(KMI))=LHORELAX_SVSNW
! Update NSV* variables for model KMI
CALL UPDATE_NSV(KMI)
!
!  SET MINIMUN VALUE FOR DIFFERENT SV GROUPS
!
XSVMIN(1:NSV_USER_A(KMI))=0.
IF (CCLOUD == 'C2R2' .OR. CCLOUD == 'C3R5' .OR.  CCLOUD == 'KHKO' ) &
XSVMIN(NSV_C2R2BEG_A(KMI):NSV_C2R2END_A(KMI))=0.
IF (CCLOUD == 'C3R5') &
XSVMIN(NSV_C1R3BEG_A(KMI):NSV_C1R3END_A(KMI))=0.
IF (CCLOUD == 'LIMA') &
XSVMIN(NSV_LIMA_BEG_A(KMI):NSV_LIMA_END_A(KMI))=0.
IF (CELEC /= 'NONE') &
XSVMIN(NSV_ELECBEG_A(KMI):NSV_ELECEND_A(KMI))=0.
IF (LUSECHEM .OR. LCHEMDIAG) &
XSVMIN(NSV_CHEMBEG_A(KMI):NSV_CHEMEND_A(KMI))=0.
IF (LUSECHIC) &
XSVMIN(NSV_CHICBEG_A(KMI):NSV_CHICEND_A(KMI))=0.
IF (.NOT.(LUSECHEM .OR. LCHEMDIAG) .AND. LCH_CONV_LINOX) &
XSVMIN(NSV_LNOXBEG_A(KMI):NSV_LNOXEND_A(KMI))=0.
IF (LORILAM .OR. LCHEMDIAG) &
XSVMIN(NSV_AERBEG_A(KMI):NSV_AEREND_A(KMI))=0.
IF (LDUST) XSVMIN(NSV_DSTBEG_A(KMI):NSV_DSTEND_A(KMI))=XMNH_TINY
IF ((LDUST).AND.(LDEPOS_DST(KMI))) &
XSVMIN(NSV_DSTDEPBEG_A(KMI):NSV_DSTDEPEND_A(KMI))=XMNH_TINY
IF (LSALT) XSVMIN(NSV_SLTBEG_A(KMI):NSV_SLTEND_A(KMI))=XMNH_TINY
IF (LLG) THEN
  XSVMIN(NSV_LGBEG_A(KMI))  =XLG1MIN
  XSVMIN(NSV_LGBEG_A(KMI)+1)=XLG2MIN
  XSVMIN(NSV_LGEND_A(KMI))  =XLG3MIN
ENDIF
IF ((LSALT).AND.(LDEPOS_SLT(KMI))) &
XSVMIN(NSV_SLTDEPBEG_A(KMI):NSV_SLTDEPEND_A(KMI))=XMNH_TINY
IF ((LORILAM).AND.(LDEPOS_AER(KMI))) &
XSVMIN(NSV_AERDEPBEG_A(KMI):NSV_AERDEPEND_A(KMI))=XMNH_TINY
IF (LPASPOL) XSVMIN(NSV_PPBEG_A(KMI):NSV_PPEND_A(KMI))=0.    
#ifdef MNH_FOREFIRE      
IF (LFOREFIRE) XSVMIN(NSV_FFBEG_A(KMI):NSV_FFEND_A(KMI))=0.
#endif
! Blaze smoke
IF (LBLAZE) XSVMIN(NSV_FIREBEG_A(KMI):NSV_FIREEND_A(KMI))=0.
!
IF (LCONDSAMP) XSVMIN(NSV_CSBEG_A(KMI):NSV_CSEND_A(KMI))=0.   
IF (LBLOWSNOW) XSVMIN(NSV_SNWBEG_A(KMI):NSV_SNWEND_A(KMI))=XMNH_TINY
!
!  NAME OF THE SCALAR VARIABLES IN THE DIFFERENT SV GROUPS
!
CSV_A(:, KMI) = '      '
IF (LLG) THEN
  CSV_A(NSV_LGBEG_A(KMI),   KMI) = 'X0     '
  CSV_A(NSV_LGBEG_A(KMI)+1, KMI) = 'Y0     '
  CSV_A(NSV_LGEND_A(KMI),   KMI) = 'Z0     '
ENDIF

! Initialize scalar variable names for dust
IF ( LDUST ) THEN
  IF ( NMODE_DST < 1 .OR. NMODE_DST > 3 ) CALL Print_msg( NVERB_FATAL, 'GEN', 'INI_NSV', 'NMODE_DST must in the 1 to 3 interval' )

  ! Initialization of dust names
  ! Was allocated for previous KMI
  ! We assume that if LDUST=T on a model, NSV_DST_A(KMI) is the same for all
  IF( .NOT. ALLOCATED( CDUSTNAMES ) ) THEN
    ALLOCATE( CDUSTNAMES(NSV_DST_A(KMI)) )
  ELSE IF ( SIZE( CDUSTNAMES ) /= NSV_DST_A(KMI) ) THEN
    CALL Print_msg( NVERB_ERROR, 'GEN', 'INI_NSV', 'NSV_DST not the same for different model (if LDUST=T)' )
    DEALLOCATE( CDUSTNAMES )
    ALLOCATE( CDUSTNAMES(NSV_DST_A(KMI)) )
  END IF
  ALLOCATE( YDUSTLONGNAMES(NSV_DST_A(KMI)) )
  !Loop on all dust modes
  IF ( INMOMENTS_DST == 1 ) THEN
    DO JMODE = 1, NMODE_DST
      IMODEIDX = JPDUSTORDER(JMODE)
      JSV_NAME = ( IMODEIDX - 1 ) * 3 + 2
      CDUSTNAMES(JMODE) = YPDUST_INI(JSV_NAME)
      !Add meaning of the ppv unit (here for moment 3)
      YDUSTLONGNAMES(JMODE) = TRIM( YPDUST_INI(JSV_NAME) ) // ' [molec_{aer}/molec_{air}]'
    END DO
  ELSE
    DO JMODE = 1,NMODE_DST
      !Find which mode we are dealing with
      IMODEIDX = JPDUSTORDER(JMODE)
      DO JMOM = 1, INMOMENTS_DST
        !Find which number this is of the list of scalars
        JSV = ( JMODE - 1 ) * INMOMENTS_DST + JMOM
        !Find what name this corresponds to, always 3 moments assumed in YPDUST_INI
        JSV_NAME = ( IMODEIDX - 1) * 3 + JMOM
        !Get the right CDUSTNAMES which should follow the list of scalars transported in XSVM/XSVT
        CDUSTNAMES(JSV) = YPDUST_INI(JSV_NAME)
        !Add meaning of the ppv unit
        IF ( JMOM == 1 ) THEN !Corresponds to moment 0
          YDUSTLONGNAMES(JSV) = TRIM( YPDUST_INI(JSV_NAME) ) // ' [nb_aerosols/molec_{air}]'
        ELSE IF ( JMOM == 2 ) THEN !Corresponds to moment 3
          YDUSTLONGNAMES(JSV) = TRIM( YPDUST_INI(JSV_NAME) ) // ' [molec_{aer}/molec_{air}]'
        ELSE IF ( JMOM == 3 ) THEN !Corresponds to moment 6
          YDUSTLONGNAMES(JSV) = TRIM( YPDUST_INI(JSV_NAME) ) // ' [um6/molec_{air}*(cm3/m3)]'
        ELSE
          CALL Print_msg( NVERB_WARNING, 'GEN', 'INI_NSV', 'unknown moment for DUST' )
          YDUSTLONGNAMES(JMODE) = TRIM( YPDUST_INI(JSV_NAME) )
        END IF
      ENDDO ! Loop on moments
    ENDDO    ! Loop on dust modes
  END IF

  ! Initialization of deposition scheme names
  IF ( LDEPOS_DST(KMI) ) THEN
    IF( .NOT. ALLOCATED( CDEDSTNAMES ) ) THEN
      ALLOCATE( CDEDSTNAMES(NMODE_DST * 2) )
      DO JMODE = 1, NMODE_DST
        IMODEIDX = JPDUSTORDER(JMODE)
        CDEDSTNAMES(JMODE)             = YPDEDST_INI(IMODEIDX)
        CDEDSTNAMES(NMODE_DST + JMODE) = YPDEDST_INI(NMODE_DST + IMODEIDX)
      ENDDO
    END IF
  END IF
END IF

! Initialize scalar variable names for salt
IF ( LSALT ) THEN
  IF ( NMODE_SLT < 1 .OR. NMODE_SLT > 8 ) CALL Print_msg( NVERB_FATAL, 'GEN', 'INI_NSV', 'NMODE_SLT must in the 1 to 8 interval' )

  ! Was allocated for previous KMI
  ! We assume that if LSALT=T on a model, NSV_SLT_A(KMI) is the same for all
  IF( .NOT. ALLOCATED( CSALTNAMES ) ) THEN
    ALLOCATE( CSALTNAMES(NSV_SLT_A(KMI)) )
  ELSE IF ( SIZE( CSALTNAMES ) /= NSV_SLT_A(KMI) ) THEN
    CALL Print_msg( NVERB_ERROR, 'GEN', 'INI_NSV', 'NSV_SLT not the same for different model (if LSALT=T)' )
    DEALLOCATE( CSALTNAMES )
    ALLOCATE( CSALTNAMES(NSV_SLT_A(KMI)) )
  END IF
  ALLOCATE( YSALTLONGNAMES(NSV_SLT_A(KMI)) )
  !Loop on all dust modes
  IF ( INMOMENTS_SLT == 1 ) THEN
    DO JMODE = 1, NMODE_SLT
      IMODEIDX = JPSALTORDER(JMODE)
      JSV_NAME = ( IMODEIDX - 1 ) * 3 + 2
      CSALTNAMES(JMODE) = YPSALT_INI(JSV_NAME)
      !Add meaning of the ppv unit (here for moment 3)
      YSALTLONGNAMES(JMODE) = TRIM( YPSALT_INI(JSV_NAME) ) // ' [molec_{aer}/molec_{air}]'
    END DO
  ELSE
    DO JMODE = 1, NMODE_SLT
      !Find which mode we are dealing with
      IMODEIDX = JPSALTORDER(JMODE)
      DO JMOM = 1, INMOMENTS_SLT
        !Find which number this is of the list of scalars
        JSV = ( JMODE - 1 ) * INMOMENTS_SLT + JMOM
        !Find what name this corresponds to, always 3 moments assumed in YPSALT_INI
        JSV_NAME = ( IMODEIDX - 1 ) * 3 + JMOM
        !Get the right CSALTNAMES which should follow the list of scalars transported in XSVM/XSVT
        CSALTNAMES(JSV) = YPSALT_INI(JSV_NAME)
        !Add meaning of the ppv unit
        IF ( JMOM == 1 ) THEN !Corresponds to moment 0
          YSALTLONGNAMES(JSV) = TRIM( YPSALT_INI(JSV_NAME) ) // ' [nb_aerosols/molec_{air}]'
        ELSE IF ( JMOM == 2 ) THEN !Corresponds to moment 3
          YSALTLONGNAMES(JSV) = TRIM( YPSALT_INI(JSV_NAME) ) // ' [molec_{aer}/molec_{air}]'
        ELSE IF ( JMOM == 3 ) THEN !Corresponds to moment 6
          YSALTLONGNAMES(JSV) = TRIM( YPSALT_INI(JSV_NAME) ) // ' [um6/molec_{air}*(cm3/m3)]'
        ELSE
          CALL Print_msg( NVERB_WARNING, 'GEN', 'INI_NSV', 'unknown moment for SALT' )
          YSALTLONGNAMES(JMODE) = TRIM( YPSALT_INI(JSV_NAME) )
        END IF
      ENDDO ! Loop on moments
    ENDDO    ! Loop on dust modes
  END IF

  ! Initialization of deposition scheme
  IF ( LDEPOS_SLT(KMI) ) THEN
    IF( .NOT. ALLOCATED( CDESLTNAMES ) ) THEN
      ALLOCATE( CDESLTNAMES(NMODE_SLT * 2) )
      DO JMODE = 1, NMODE_SLT
        IMODEIDX = JPSALTORDER(JMODE)
        CDESLTNAMES(JMODE)             = YPDESLT_INI(IMODEIDX)
        CDESLTNAMES(NMODE_SLT + JMODE) = YPDESLT_INI(NMODE_SLT + IMODEIDX)
      ENDDO
    ENDIF
  ENDIF
END IF

! Initialize scalar variable names for snow
IF ( LBLOWSNOW ) THEN
  IF( .NOT. ALLOCATED( CSNOWNAMES ) ) THEN
    ALLOCATE( CSNOWNAMES(NSV_SNW_A(KMI)) )
    DO JMOM = 1, NSV_SNW_A(KMI)
      CSNOWNAMES(JMOM) = YPSNOW_INI(JMOM)
    END DO
  END IF
END IF

!Fill metadata for model KMI
DO JSV = 1, NSV_USER_A(KMI)
  WRITE( YNUM3, '( I3.3 )' ) JSV

  TSVLIST_A(JSV, KMI) = TFIELDMETADATA(         &
    CMNHNAME   = 'SVUSER' // YNUM3,             &
    CSTDNAME   = '',                            &
    CLONGNAME  = 'SVUSER' // YNUM3,             &
    CUNITS     = 'kg kg-1',                     &
    CDIR       = 'XY',                          &
    CCOMMENT   = 'X_Y_Z_' // 'SVUSER' // YNUM3, &
    NGRID      = 1,                             &
    NTYPE      = TYPEREAL,                      &
    NDIMS      = 3,                             &
    LTIMEDEP   = .TRUE.                         )
END DO

DO JSV = NSV_C2R2BEG_A(KMI), NSV_C2R2END_A(KMI)
  WRITE( YNUM3, '( I3.3 )' ) JSV

  TSVLIST_A(JSV, KMI) = TFIELDMETADATA(                       &
    CMNHNAME   = TRIM( C2R2NAMES(JSV-NSV_C2R2BEG_A(KMI)+1) ), &
    CSTDNAME   = '',                                          &
    CLONGNAME  = TRIM( C2R2NAMES(JSV-NSV_C2R2BEG_A(KMI)+1) ), &
    CUNITS     = 'm-3',                                       &
    CDIR       = 'XY',                                        &
    CCOMMENT   = 'X_Y_Z_' // 'SVT' // YNUM3,                  &
    NGRID      = 1,                                           &
    NTYPE      = TYPEREAL,                                    &
    NDIMS      = 3,                                           &
    LTIMEDEP   = .TRUE.                                       )
END DO

DO JSV = NSV_C1R3BEG_A(KMI), NSV_C1R3END_A(KMI)
  WRITE( YNUM3, '( I3.3 )' ) JSV

  TSVLIST_A(JSV, KMI) = TFIELDMETADATA(                       &
    CMNHNAME   = TRIM( C1R3NAMES(JSV-NSV_C2R2BEG_A(KMI)+1) ), &
    CSTDNAME   = '',                                          &
    CLONGNAME  = TRIM( C1R3NAMES(JSV-NSV_C2R2BEG_A(KMI)+1) ), &
    CUNITS     = 'm-3',                                       &
    CDIR       = 'XY',                                        &
    CCOMMENT   = 'X_Y_Z_' // 'SVT' // YNUM3,                  &
    NGRID      = 1,                                           &
    NTYPE      = TYPEREAL,                                    &
    NDIMS      = 3,                                           &
    LTIMEDEP   = .TRUE.                                       )
END DO

DO JSV = NSV_LIMA_BEG_A(KMI), NSV_LIMA_END_A(KMI)
  WRITE( YNUM3, '( I3.3 )' ) JSV

  TSVLIST_A(JSV, KMI) = TFIELDMETADATA(      &
    CMNHNAME   = 'SV LIMA ' // YNUM3,        &
    CSTDNAME   = '',                         &
    CLONGNAME  = '',                         &
    CUNITS     = 'kg-1',                     &
    CDIR       = 'XY',                       &
    CCOMMENT   = 'X_Y_Z_' // 'SVT' // YNUM3, &
    NGRID      = 1,                          &
    NTYPE      = TYPEREAL,                   &
    NDIMS      = 3,                          &
    LTIMEDEP   = .TRUE.                      )

  IF ( JSV == NSV_LIMA_NC_A(KMI) ) THEN
    TSVLIST_A(JSV, KMI)%CMNHNAME = TRIM( CLIMA_WARM_NAMES(1) )
  ELSE IF ( JSV == NSV_LIMA_NR_A(KMI) ) THEN
    TSVLIST_A(JSV, KMI)%CMNHNAME = TRIM( CLIMA_WARM_NAMES(2) )
  ELSE IF ( JSV >= NSV_LIMA_CCN_FREE_A(KMI) .AND. JSV < NSV_LIMA_CCN_ACTI_A(KMI) ) THEN
    WRITE( YNUM2, '( I2.2 )' ) JSV - NSV_LIMA_CCN_FREE_A(KMI) + 1
    TSVLIST_A(JSV, KMI)%CMNHNAME = TRIM( CLIMA_WARM_NAMES(3) ) // YNUM2
  ELSE IF (JSV >= NSV_LIMA_CCN_ACTI_A(KMI) .AND. JSV < ( NSV_LIMA_CCN_ACTI_A(KMI) + NMOD_CCN ) ) THEN
    WRITE( YNUM2, '( I2.2 )' ) JSV - NSV_LIMA_CCN_ACTI_A(KMI) + 1
    TSVLIST_A(JSV, KMI)%CMNHNAME = TRIM( CLIMA_WARM_NAMES(4) ) // YNUM2
  ELSE IF ( JSV == NSV_LIMA_SCAVMASS_A(KMI) ) THEN
    TSVLIST_A(JSV, KMI)%CMNHNAME = TRIM( CAERO_MASS(1) )
    TSVLIST_A(JSV, KMI)%CUNITS = 'kg kg-1'
  ELSE IF ( JSV == NSV_LIMA_NI_A(KMI) ) THEN
    TSVLIST_A(JSV, KMI)%CMNHNAME = TRIM( CLIMA_COLD_NAMES(1) )
  ELSE IF ( JSV == NSV_LIMA_NS_A(KMI) ) THEN
    TSVLIST_A(JSV, KMI)%CMNHNAME = TRIM( CLIMA_COLD_NAMES(2) )
  ELSE IF ( JSV == NSV_LIMA_NG_A(KMI) ) THEN
    TSVLIST_A(JSV, KMI)%CMNHNAME = TRIM( CLIMA_COLD_NAMES(3) )
  ELSE IF ( JSV == NSV_LIMA_NH_A(KMI) ) THEN
    TSVLIST_A(JSV, KMI)%CMNHNAME = TRIM( CLIMA_COLD_NAMES(4) )
  ELSE IF ( JSV >= NSV_LIMA_IFN_FREE_A(KMI) .AND. JSV < NSV_LIMA_IFN_NUCL_A(KMI) ) THEN
    WRITE( YNUM2, '( I2.2 )' ) JSV - NSV_LIMA_IFN_FREE_A(KMI) + 1
    TSVLIST_A(JSV, KMI)%CMNHNAME = TRIM( CLIMA_COLD_NAMES(5) ) // YNUM2
  ELSE IF ( JSV >= NSV_LIMA_IFN_NUCL_A(KMI) .AND. JSV < ( NSV_LIMA_IFN_NUCL_A(KMI) + NMOD_IFN ) ) THEN
    WRITE( YNUM2, '( I2.2 )' ) JSV - NSV_LIMA_IFN_NUCL_A(KMI) + 1
    TSVLIST_A(JSV, KMI)%CMNHNAME = TRIM( CLIMA_COLD_NAMES(6) ) // YNUM2
  ELSE IF ( JSV >= NSV_LIMA_IMM_NUCL_A(KMI) .AND. JSV < ( NSV_LIMA_IMM_NUCL_A(KMI) + NMOD_IMM ) ) THEN
    WRITE( YNUM2, '( I2.2 )' ) NINDICE_CCN_IMM(JSV-NSV_LIMA_IMM_NUCL_A(KMI)+1)
    TSVLIST_A(JSV, KMI)%CMNHNAME = TRIM( CLIMA_COLD_NAMES(7) ) // YNUM2
  ELSE IF ( JSV == NSV_LIMA_HOM_HAZE_A(KMI) ) THEN
    TSVLIST_A(JSV, KMI)%CMNHNAME = TRIM( CLIMA_COLD_NAMES(8) )
  ELSE IF ( JSV == NSV_LIMA_SPRO_A(KMI) ) THEN
    TSVLIST_A(JSV, KMI)%CUNITS = '1'
    TSVLIST_A(JSV, KMI)%CMNHNAME = TRIM( CLIMA_WARM_NAMES(5) )
  ELSE
    CALL Print_msg( NVERB_FATAL, 'GEN', 'INI_NSV', 'invalid index for LIMA' )
  END IF

   TSVLIST_A(JSV, KMI)%CLONGNAME = TRIM( TSVLIST_A(JSV, KMI)%CMNHNAME )
END DO

DO JSV = NSV_ELECBEG_A(KMI), NSV_ELECEND_A(KMI)
  IF ( JSV > NSV_ELECBEG .AND. JSV < NSV_ELECEND ) THEN
    YUNITS = 'C kg-1'
    WRITE( YCOMMENT, '( A6, A3, I3.3 )' ) 'X_Y_Z_', 'SVT', JSV
  ELSE
    YUNITS = 'kg-1'
    WRITE( YCOMMENT, '( A6, A3, I3.3, A8 )' ) 'X_Y_Z_', 'SVT', JSV, ' (nb ions/kg)'
  END IF

  TSVLIST_A(JSV, KMI) = TFIELDMETADATA(                        &
    CMNHNAME   = TRIM( CELECNAMES(JSV-NSV_ELECBEG_A(KMI)+1) ), &
    CSTDNAME   = '',                                           &
    CLONGNAME  = TRIM( CELECNAMES(JSV-NSV_ELECBEG_A(KMI)+1) ), &
    CUNITS     = TRIM( YUNITS ),                               &
    CDIR       = 'XY',                                         &
    CCOMMENT   = TRIM( YCOMMENT ),                             &
    NGRID      = 1,                                            &
    NTYPE      = TYPEREAL,                                     &
    NDIMS      = 3,                                            &
    LTIMEDEP   = .TRUE.                                        )
END DO

DO JSV = NSV_LGBEG_A(KMI), NSV_LGEND_A(KMI)
  WRITE( YNUM3, '( I3.3 )' ) JSV

  TSVLIST_A(JSV, KMI) = TFIELDMETADATA(                    &
    CMNHNAME   = TRIM( CLGNAMES(JSV-NSV_LGBEG_A(KMI)+1) ), &
    CSTDNAME   = '',                                       &
    CLONGNAME  = TRIM( CLGNAMES(JSV-NSV_LGBEG_A(KMI)+1) ), &
    CUNITS     = 'm',                                      &
    CDIR       = 'XY',                                     &
    CCOMMENT   = 'X_Y_Z_' // 'SVT' // YNUM3,               &
    NGRID      = 1,                                        &
    NTYPE      = TYPEREAL,                                 &
    NDIMS      = 3,                                        &
    LTIMEDEP   = .TRUE.                                    )
END DO

DO JSV = NSV_PPBEG_A(KMI), NSV_PPEND_A(KMI)
  WRITE( YNUM3, '( I3.3 )' ) JSV-NSV_PPBEG_A(KMI)+1

  TSVLIST_A(JSV, KMI) = TFIELDMETADATA(      &
    CMNHNAME   = 'SVPP' // YNUM3,            &
    CSTDNAME   = '',                         &
    CLONGNAME  = 'SVPP' // YNUM3,            &
    CUNITS     = 'kg kg-1',                  &
    CDIR       = 'XY',                       &
    CCOMMENT   = 'X_Y_Z_' // 'SVT' // YNUM3, &
    NGRID      = 1,                          &
    NTYPE      = TYPEREAL,                   &
    NDIMS      = 3,                          &
    LTIMEDEP   = .TRUE.                      )
END DO

#ifdef MNH_FOREFIRE
DO JSV = NSV_FFBEG_A(KMI), NSV_FFEND_A(KMI)
  WRITE( YNUM3, '( I3.3 )' ) JSV-NSV_FFBEG_A(KMI)+1

  TSVLIST_A(JSV, KMI) = TFIELDMETADATA(      &
    CMNHNAME   = 'SVFF' // YNUM3,            &
    CSTDNAME   = '',                         &
    CLONGNAME  = 'SVFF' // YNUM3,            &
    CUNITS     = 'kg kg-1',                  &
    CDIR       = 'XY',                       &
    CCOMMENT   = 'X_Y_Z_' // 'SVT' // YNUM3, &
    NGRID      = 1,                          &
    NTYPE      = TYPEREAL,                   &
    NDIMS      = 3,                          &
    LTIMEDEP   = .TRUE.                      )
END DO
#endif

DO JSV = NSV_FIREBEG_A(KMI), NSV_FIREEND_A(KMI)
  WRITE( YNUM3, '( I3.3 )' ) JSV-NSV_FIREBEG_A(KMI)+1

  TSVLIST_A(JSV, KMI) = TFIELDMETADATA(      &
    CMNHNAME   = 'SVFIRE' // YNUM3,          &
    CSTDNAME   = '',                         &
    CLONGNAME  = 'SVFIRE' // YNUM3,          &
    CUNITS     = 'kg kg-1',                  &
    CDIR       = 'XY',                       &
    CCOMMENT   = 'X_Y_Z_' // 'SVT' // YNUM3, &
    NGRID      = 1,                          &
    NTYPE      = TYPEREAL,                   &
    NDIMS      = 3,                          &
    LTIMEDEP   = .TRUE.                      )
END DO

DO JSV = NSV_CSBEG_A(KMI), NSV_CSEND_A(KMI)
  WRITE( YNUM3, '( I3.3 )' ) JSV-NSV_CSBEG_A(KMI)

  TSVLIST_A(JSV, KMI) = TFIELDMETADATA(      &
    CMNHNAME   = 'SVCS' // YNUM3,            &
    CSTDNAME   = '',                         &
    CLONGNAME  = 'SVCS' // YNUM3,            &
    CUNITS     = 'kg kg-1',                  &
    CDIR       = 'XY',                       &
    CCOMMENT   = 'X_Y_Z_' // 'SVT' // YNUM3, &
    NGRID      = 1,                          &
    NTYPE      = TYPEREAL,                   &
    NDIMS      = 3,                          &
    LTIMEDEP   = .TRUE.                      )
END DO

DO JSV = NSV_CHEMBEG_A(KMI), NSV_CHEMEND_A(KMI)
  ICHIDX = ICHIDX + 1
  CSV_CHEM_LIST_A(ICHIDX, KMI) = TRIM( CNAMES(JSV-NSV_CHEMBEG_A(KMI)+1) )

  WRITE( YNUM3, '( I3.3 )' ) JSV

  TSVLIST_A(JSV, KMI) = TFIELDMETADATA(                    &
    CMNHNAME   = TRIM( CNAMES(JSV-NSV_CHEMBEG_A(KMI)+1) ), &
    CSTDNAME   = '',                                       &
    CLONGNAME  = TRIM( CNAMES(JSV-NSV_CHEMBEG_A(KMI)+1) ), &
    CUNITS     = 'ppv',                                    &
    CDIR       = 'XY',                                     &
    CCOMMENT   = 'X_Y_Z_' // 'SVT' // YNUM3,               &
    NGRID      = 1,                                        &
    NTYPE      = TYPEREAL,                                 &
    NDIMS      = 3,                                        &
    LTIMEDEP   = .TRUE.                                    )
END DO

DO JSV = NSV_CHICBEG_A(KMI), NSV_CHICEND_A(KMI)
  ICHIDX = ICHIDX + 1
  CSV_CHEM_LIST_A(ICHIDX, KMI) = TRIM( CICNAMES(JSV-NSV_CHICBEG_A(KMI)+1) )

  WRITE( YNUM3, '( I3.3 )' ) JSV

  TSVLIST_A(JSV, KMI) = TFIELDMETADATA(                      &
    CMNHNAME   = TRIM( CICNAMES(JSV-NSV_CHICBEG_A(KMI)+1) ), &
    CSTDNAME   = '',                                         &
    CLONGNAME  = TRIM( CICNAMES(JSV-NSV_CHICBEG_A(KMI)+1) ), &
    CUNITS     = 'ppv',                                      &
    CDIR       = 'XY',                                       &
    CCOMMENT   = 'X_Y_Z_' // 'SVT' // YNUM3,                 &
    NGRID      = 1,                                          &
    NTYPE      = TYPEREAL,                                   &
    NDIMS      = 3,                                          &
    LTIMEDEP   = .TRUE.                                      )
END DO

DO JSV = NSV_AERBEG_A(KMI), NSV_AEREND_A(KMI)
  ICHIDX = ICHIDX + 1
  CSV_CHEM_LIST_A(ICHIDX, KMI) = TRIM( CAERONAMES(JSV-NSV_AERBEG_A(KMI)+1) )

  WRITE( YNUM3, '( I3.3 )' ) JSV

  !Determine moment to add meaning of the ppv unit
  JAER = JSV - NSV_AERBEG_A(KMI) + 1
  IF ( ANY( JAER == [JP_CH_M0i, JP_CH_M0j] ) ) THEN
    !Moment 0
    YAEROLONGNAMES = TRIM( CAERONAMES(JAER) ) // ' [nb_aerosols/molec_{air}]'
  ELSE IF ( ANY( JAER == [ JP_CH_SO4i, JP_CH_SO4j, JP_CH_NO3i, JP_CH_NO3j, JP_CH_H2Oi, JP_CH_H2Oj, JP_CH_NH3i, JP_CH_NH3j,   &
                           JP_CH_OCi,  JP_CH_OCj,  JP_CH_BCi,  JP_CH_BCj,  JP_CH_DSTi, JP_CH_DSTj ] )                        &
            .OR. ( NSOA == 10 .AND.                                                                                          &
                   ANY( JAER == [ JP_CH_SOA1i, JP_CH_SOA1j, JP_CH_SOA2i, JP_CH_SOA2j, JP_CH_SOA3i, JP_CH_SOA3j, JP_CH_SOA4i, &
                                  JP_CH_SOA4j, JP_CH_SOA5i, JP_CH_SOA5j, JP_CH_SOA6i, JP_CH_SOA6j, JP_CH_SOA7i, JP_CH_SOA7j, &
                                  JP_CH_SOA8i, JP_CH_SOA8j, JP_CH_SOA9i, JP_CH_SOA9j, JP_CH_SOA10i, JP_CH_SOA10j ] )       ) ) THEN
    !Moment 3
    YAEROLONGNAMES = TRIM( CAERONAMES(JAER) ) // ' [molec_{aer}/molec_{air}]'
  ELSE IF ( ( LVARSIGI .AND. JAER == JP_CH_M6i ) .OR. ( LVARSIGJ .AND. JAER == JP_CH_M6j ) ) THEN
    !Moment 6
    YAEROLONGNAMES = TRIM( CAERONAMES(JAER) ) // ' [um6/molec_{air}*(cm3/m3)]'
  ELSE
    CALL Print_msg( NVERB_WARNING, 'GEN', 'INI_NSV', 'unknown moment for AER' )
    YAEROLONGNAMES = TRIM( CAERONAMES(JAER) )
  END IF

  TSVLIST_A(JSV, KMI) = TFIELDMETADATA(                           &
    CMNHNAME   = TRIM( CAERONAMES(JSV-NSV_AERBEG_A(KMI)+1) ),     &
    CSTDNAME   = '',                                              &
    CLONGNAME  = TRIM( YAEROLONGNAMES(JSV-NSV_AERBEG_A(KMI)+1) ), &
    CUNITS     = 'ppv',                                           &
    CDIR       = 'XY',                                            &
    CCOMMENT   = 'X_Y_Z_' // 'SVT' // YNUM3,                      &
    NGRID      = 1,                                               &
    NTYPE      = TYPEREAL,                                        &
    NDIMS      = 3,                                               &
    LTIMEDEP   = .TRUE.                                           )
END DO

DO JSV = NSV_AERDEPBEG_A(KMI), NSV_AERDEPEND_A(KMI)
  ICHIDX = ICHIDX + 1
  CSV_CHEM_LIST_A(ICHIDX, KMI) = TRIM( CDEAERNAMES(JSV-NSV_AERDEPBEG_A(KMI)+1) )

  WRITE( YNUM3, '( I3.3 )' ) JSV

  TSVLIST_A(JSV, KMI) = TFIELDMETADATA(                           &
    CMNHNAME   = TRIM( CDEAERNAMES(JSV-NSV_AERDEPBEG_A(KMI)+1) ), &
    CSTDNAME   = '',                                              &
    CLONGNAME  = TRIM( CDEAERNAMES(JSV-NSV_AERDEPBEG_A(KMI)+1) ), &
    CUNITS     = 'ppv',                                           &
    CDIR       = 'XY',                                            &
    CCOMMENT   = 'X_Y_Z_' // 'SVT' // YNUM3,                      &
    NGRID      = 1,                                               &
    NTYPE      = TYPEREAL,                                        &
    NDIMS      = 3,                                               &
    LTIMEDEP   = .TRUE.                                           )
END DO

DO JSV = NSV_DSTBEG_A(KMI), NSV_DSTEND_A(KMI)
  ICHIDX = ICHIDX + 1
  CSV_CHEM_LIST_A(ICHIDX, KMI) = TRIM( CDUSTNAMES(JSV-NSV_DSTBEG_A(KMI)+1) )

  WRITE( YNUM3, '( I3.3 )' ) JSV

  TSVLIST_A(JSV, KMI) = TFIELDMETADATA(                           &
    CMNHNAME   = TRIM( CDUSTNAMES(JSV-NSV_DSTBEG_A(KMI)+1) ),     &
    CSTDNAME   = '',                                              &
    CLONGNAME  = TRIM( YDUSTLONGNAMES(JSV-NSV_DSTBEG_A(KMI)+1) ), &
    CUNITS     = 'ppv',                                           &
    CDIR       = 'XY',                                            &
    CCOMMENT   = 'X_Y_Z_' // 'SVT' // YNUM3,                      &
    NGRID      = 1,                                               &
    NTYPE      = TYPEREAL,                                        &
    NDIMS      = 3,                                               &
    LTIMEDEP   = .TRUE.                                           )
END DO

DO JSV = NSV_DSTDEPBEG_A(KMI), NSV_DSTDEPEND_A(KMI)
  ICHIDX = ICHIDX + 1
  CSV_CHEM_LIST_A(ICHIDX, KMI) = TRIM( CDEDSTNAMES(JSV-NSV_DSTDEPBEG_A(KMI)+1) )

  WRITE( YNUM3, '( I3.3 )' ) JSV

  TSVLIST_A(JSV, KMI) = TFIELDMETADATA(                           &
    CMNHNAME   = TRIM( CDEDSTNAMES(JSV-NSV_DSTDEPBEG_A(KMI)+1) ), &
    CSTDNAME   = '',                                              &
    CLONGNAME  = TRIM( CDEDSTNAMES(JSV-NSV_DSTDEPBEG_A(KMI)+1) ), &
    CUNITS     = 'ppv',                                           &
    CDIR       = 'XY',                                            &
    CCOMMENT   = 'X_Y_Z_' // 'SVT' // YNUM3,                      &
    NGRID      = 1,                                               &
    NTYPE      = TYPEREAL,                                        &
    NDIMS      = 3,                                               &
    LTIMEDEP   = .TRUE.                                           )
END DO

DO JSV = NSV_SLTBEG_A(KMI), NSV_SLTEND_A(KMI)
  ICHIDX = ICHIDX + 1
  CSV_CHEM_LIST_A(ICHIDX, KMI) = TRIM( CSALTNAMES(JSV-NSV_SLTBEG_A(KMI)+1) )

  WRITE( YNUM3, '( I3.3 )' ) JSV

  TSVLIST_A(JSV, KMI) = TFIELDMETADATA(                           &
    CMNHNAME   = TRIM( CSALTNAMES(JSV-NSV_SLTBEG_A(KMI)+1) ),     &
    CSTDNAME   = '',                                              &
    CLONGNAME  = TRIM( YSALTLONGNAMES(JSV-NSV_SLTBEG_A(KMI)+1) ), &
    CUNITS     = 'ppv',                                           &
    CDIR       = 'XY',                                            &
    CCOMMENT   = 'X_Y_Z_' // 'SVT' // YNUM3,                      &
    NGRID      = 1,                                               &
    NTYPE      = TYPEREAL,                                        &
    NDIMS      = 3,                                               &
    LTIMEDEP   = .TRUE.                                           )
END DO

DO JSV = NSV_SLTDEPBEG_A(KMI), NSV_SLTDEPEND_A(KMI)
  ICHIDX = ICHIDX + 1
  CSV_CHEM_LIST_A(ICHIDX, KMI) = TRIM( CDESLTNAMES(JSV-NSV_SLTDEPBEG_A(KMI)+1) )

  WRITE( YNUM3, '( I3.3 )' ) JSV

  TSVLIST_A(JSV, KMI) = TFIELDMETADATA(                           &
    CMNHNAME   = TRIM( CDESLTNAMES(JSV-NSV_SLTDEPBEG_A(KMI)+1) ), &
    CSTDNAME   = '',                                              &
    CLONGNAME  = TRIM( CDESLTNAMES(JSV-NSV_SLTDEPBEG_A(KMI)+1) ), &
    CUNITS     = 'ppv',                                           &
    CDIR       = 'XY',                                            &
    CCOMMENT   = 'X_Y_Z_' // 'SVT' // YNUM3,                      &
    NGRID      = 1,                                               &
    NTYPE      = TYPEREAL,                                        &
    NDIMS      = 3,                                               &
    LTIMEDEP   = .TRUE.                                           )
END DO

DO JSV = NSV_SNWBEG_A(KMI), NSV_SNWEND_A(KMI)
  WRITE( YNUM3, '( I3.3 )' ) JSV

  TSVLIST_A(JSV, KMI) = TFIELDMETADATA(                       &
    CMNHNAME   = TRIM( CSNOWNAMES(JSV-NSV_SNWBEG_A(KMI)+1) ), &
    CSTDNAME   = '',                                          &
    CLONGNAME  = TRIM( CSNOWNAMES(JSV-NSV_SNWBEG_A(KMI)+1) ), &
    CUNITS     = 'kg kg-1',                                   &
    CDIR       = 'XY',                                        &
    CCOMMENT   = 'X_Y_Z_' // 'SVT' // YNUM3,                  &
    NGRID      = 1,                                           &
    NTYPE      = TYPEREAL,                                    &
    NDIMS      = 3,                                           &
    LTIMEDEP   = .TRUE.                                       )
END DO

!Check if there is at most 1 LINOX scalar variable
!if not, the name must be modified and different for all of them
IF ( NSV_LNOX_A(KMI) > 1 ) &
  CALL Print_msg( NVERB_ERROR, 'GEN', 'INI_NSV', 'NSV_LNOX_A>1: problem with the names of the corresponding scalar variables' )

DO JSV = NSV_LNOXBEG_A(KMI), NSV_LNOXEND_A(KMI)
  ICHIDX = ICHIDX + 1
  CSV_CHEM_LIST_A(ICHIDX, KMI) = 'LINOX'

  WRITE( YNUM3, '( I3.3 )' ) JSV

  TSVLIST_A(JSV, KMI) = TFIELDMETADATA(      &
    CMNHNAME   = 'LINOX',                    &
    CSTDNAME   = '',                         &
    CLONGNAME  = 'LINOX',                    &
    CUNITS     = 'ppv',                      &
    CDIR       = 'XY',                       &
    CCOMMENT   = 'X_Y_Z_' // 'SVT' // YNUM3, &
    NGRID      = 1,                          &
    NTYPE      = TYPEREAL,                   &
    NDIMS      = 3,                          &
    LTIMEDEP   = .TRUE.                      )
END DO

IF ( ICHIDX /= NSV_CHEM_LIST_A(KMI) ) &
  CALL Print_msg( NVERB_ERROR, 'GEN', 'INI_NSV', 'ICHIDX /= NSV_CHEM_LIST_A(KMI)' )

END SUBROUTINE INI_NSV
