!MNH_LIC Copyright 2001-2023 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-------------------------------------------------------------------------------
!        ###############
         MODULE MODD_NSV
!        ###############
!
!!****  *MODD_NSV* - declaration of scalar variables numbers
!!
!!    PURPOSE
!!    -------
!!       Arrays to store the per-model NSV_* values number (suffix _A denote an array)
!!
!!    AUTHOR
!!    ------
!!      D. Gazen   L.A.
!!
!!    MODIFICATIONS
!!    -------------
!!      Original  01/02/01
!!       J.-P. Pinty   29/11/02 add C3R5, ELEC
!!       V. Masson     01/2004  add scalar names
!!       M. Leriche    12/04/07 add aqueous chemistry
!!       M. Leriche    08/07/10 add ice phase chemistry
!!       C.Lac         07/11    add conditional sampling
!!       Pialat/Tulet  15/02/12 add ForeFire
!!      Modification    01/2016  (JP Pinty) Add LIMA
!!       V. Vionnet     07/17   add blowing snow
!  B. Vie         06/2021: add prognostic supersaturation for LIMA
!  P. Wautelet 26/11/2021: add TSVLIST and TSVLIST_A to store the metadata of all the scalar variables
!  A. Costes      12/2021: add Blaze fire model smoke
!  P. Wautelet 14/01/2022: add CSV_CHEM_LIST(_A) to store the list of all chemical variables
!                          + NSV_CHEM_LIST(_A) the size of the list
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_FIELD,      ONLY: tfieldmetadata
USE MODD_PARAMETERS, ONLY: JPMODELMAX, &     ! Maximum allowed number of nested models
                           JPSVMAX,    &     ! Maximum number of scalar variables
                           JPSVNAMELGTMAX, & ! Maximum length of a scalar variable name
                           NMNHNAMELGTMAX
!
IMPLICIT NONE
SAVE
!
REAL,DIMENSION(JPSVMAX) :: XSVMIN ! minimum value for SV variables
!
LOGICAL :: LINI_NSV(JPMODELMAX) = .FALSE. ! becomes True when routine INI_NSV is called
!
CHARACTER(LEN=NMNHNAMELGTMAX), DIMENSION(:,:), ALLOCATABLE, TARGET :: CSV_CHEM_LIST_A !Names of all the chemical variables
CHARACTER(LEN=6),              DIMENSION(:,:), ALLOCATABLE, TARGET :: CSV_A           !Names of the scalar variables
TYPE(tfieldmetadata), DIMENSION(:,:), ALLOCATABLE, TARGET :: TSVLIST_A !Metadata of all the scalar variables

INTEGER,DIMENSION(JPMODELMAX)::NSV_A = 0 ! total number of scalar variables
                                         ! NSV_A = NSV_USER_A+NSV_C2R2_A+NSV_CHEM_A+..
INTEGER,DIMENSION(JPMODELMAX)::NSV_CHEM_LIST_A = 0 ! total number of chemical variables (including dust, salt...)
INTEGER,DIMENSION(JPMODELMAX)::NSV_USER_A = 0  ! number of user scalar variables with 
                                               ! indices in the range : 1...NSV_USER_A
!
INTEGER,DIMENSION(JPMODELMAX)::NSV_C2R2_A = 0    ! number of liq scalar in C2R2
                                                 !                  and in C3R5
INTEGER,DIMENSION(JPMODELMAX)::NSV_C2R2BEG_A = 0 ! with indices in the range :
INTEGER,DIMENSION(JPMODELMAX)::NSV_C2R2END_A = 0 ! NSV_C2R2BEG_A...NSV_C2R2END_A
!
INTEGER,DIMENSION(JPMODELMAX)::NSV_C1R3_A = 0    ! number of ice scalar in C3R5
INTEGER,DIMENSION(JPMODELMAX)::NSV_C1R3BEG_A = 0 ! with indices in the range :
INTEGER,DIMENSION(JPMODELMAX)::NSV_C1R3END_A = 0 ! NSV_C1R3BEG_A...NSV_C1R3END_A
!
INTEGER,DIMENSION(JPMODELMAX)::NSV_ELEC_A = 0    ! number of scalar in ELEC
INTEGER,DIMENSION(JPMODELMAX)::NSV_ELECBEG_A = 0 ! with indices in the range :
INTEGER,DIMENSION(JPMODELMAX)::NSV_ELECEND_A = 0 ! NSV_ELECBEG_A...NSV_ELECEND_A
!
INTEGER,DIMENSION(JPMODELMAX)::NSV_CHEM_A = 0    ! number of chemical scalar
INTEGER,DIMENSION(JPMODELMAX)::NSV_CHEMBEG_A = 0 ! with indices in the range :
INTEGER,DIMENSION(JPMODELMAX)::NSV_CHEMEND_A = 0 ! NSV_CHEMBEG_A...NSV_CHEMEND_A
!
INTEGER,DIMENSION(JPMODELMAX)::NSV_CHGS_A = 0    ! number of gaseous chemcial species
INTEGER,DIMENSION(JPMODELMAX)::NSV_CHGSBEG_A = 0 ! with indices
INTEGER,DIMENSION(JPMODELMAX)::NSV_CHGSEND_A = 0 ! NSV_CHGSBEG_
!
INTEGER,DIMENSION(JPMODELMAX)::NSV_CHAC_A = 0    ! number of aqueous chemical species
INTEGER,DIMENSION(JPMODELMAX)::NSV_CHACBEG_A = 0 ! with indices
INTEGER,DIMENSION(JPMODELMAX)::NSV_CHACEND_A = 0 ! NSV_CHACBEG
!
INTEGER,DIMENSION(JPMODELMAX)::NSV_CHIC_A = 0    ! number of ice phase chemical species
INTEGER,DIMENSION(JPMODELMAX)::NSV_CHICBEG_A = 0 ! with indices
INTEGER,DIMENSION(JPMODELMAX)::NSV_CHICEND_A = 0 ! NSV_CHICBEG
!
INTEGER,DIMENSION(JPMODELMAX)::NSV_LG_A = 0    ! number of LaGrangian 
INTEGER,DIMENSION(JPMODELMAX)::NSV_LGBEG_A = 0 ! with indices in the range :
INTEGER,DIMENSION(JPMODELMAX)::NSV_LGEND_A = 0 ! NSV_LGBEG_A...NSV_LGEND_A
!
INTEGER,DIMENSION(JPMODELMAX)::NSV_LNOX_A = 0    ! number of lightning NOx
INTEGER,DIMENSION(JPMODELMAX)::NSV_LNOXBEG_A = 0 ! with indices in the range :
INTEGER,DIMENSION(JPMODELMAX)::NSV_LNOXEND_A = 0 ! NSV_LNOXBEG_A...NSV_LNOXEND_A
INTEGER,DIMENSION(JPMODELMAX)::NSV_DST_A = 0    ! number of dust scalar
INTEGER,DIMENSION(JPMODELMAX)::NSV_DSTBEG_A = 0 ! with indices in the range :
INTEGER,DIMENSION(JPMODELMAX)::NSV_DSTEND_A = 0 ! NSV_DSTBEG_A...NSV_DSTEND_A
!
INTEGER,DIMENSION(JPMODELMAX)::NSV_SLT_A = 0    ! number of sea salt scalar
INTEGER,DIMENSION(JPMODELMAX)::NSV_SLTBEG_A = 0 ! with indices in the range :
INTEGER,DIMENSION(JPMODELMAX)::NSV_SLTEND_A = 0 ! NSV_SLTBEG_A...NSV_SLTEND_A
!
INTEGER,DIMENSION(JPMODELMAX)::NSV_AER_A = 0    ! number of aerosol scalar
INTEGER,DIMENSION(JPMODELMAX)::NSV_AERBEG_A = 0 ! with indices in the range :
INTEGER,DIMENSION(JPMODELMAX)::NSV_AEREND_A = 0 ! NSV_AERBEG_A...NSV_AEREND_A
!
INTEGER,DIMENSION(JPMODELMAX)::NSV_DSTDEP_A    = 0    ! number of aerosol scalar
INTEGER,DIMENSION(JPMODELMAX)::NSV_DSTDEPBEG_A = 0 ! with indices in the range :
INTEGER,DIMENSION(JPMODELMAX)::NSV_DSTDEPEND_A = 0 ! NSV_AERBEG_A...NSV_AEREND_A
!
INTEGER,DIMENSION(JPMODELMAX)::NSV_AERDEP_A    = 0    ! number of aerosol scalar
INTEGER,DIMENSION(JPMODELMAX)::NSV_AERDEPBEG_A = 0 ! with indices in the range :
INTEGER,DIMENSION(JPMODELMAX)::NSV_AERDEPEND_A = 0 ! NSV_AERBEG_A...NSV_AEREND_A
!
INTEGER,DIMENSION(JPMODELMAX)::NSV_SLTDEP_A    = 0    ! number of aerosol scalar
INTEGER,DIMENSION(JPMODELMAX)::NSV_SLTDEPBEG_A = 0 ! with indices in the range :
INTEGER,DIMENSION(JPMODELMAX)::NSV_SLTDEPEND_A = 0 ! NSV_SLTBEG_A...NSV_SLTEND_A
!
INTEGER,DIMENSION(JPMODELMAX)::NSV_PP_A = 0    ! number of passive pol.
INTEGER,DIMENSION(JPMODELMAX)::NSV_PPBEG_A = 0 ! with indices in the range :
INTEGER,DIMENSION(JPMODELMAX)::NSV_PPEND_A = 0 ! NSV_PPBEG_A...NSV_PPEND_A
!
INTEGER,DIMENSION(JPMODELMAX)::NSV_CS_A = 0    ! number of condit.samplings
INTEGER,DIMENSION(JPMODELMAX)::NSV_CSBEG_A = 0 ! with indices in the range :
INTEGER,DIMENSION(JPMODELMAX)::NSV_CSEND_A = 0 ! NSV_CSBEG_A...NSV_CSEND_A
!
INTEGER,DIMENSION(JPMODELMAX)::NSV_LIMA_A = 0     ! number of scalar in LIMA
INTEGER,DIMENSION(JPMODELMAX)::NSV_LIMA_BEG_A = 0 ! with indices in the range :
INTEGER,DIMENSION(JPMODELMAX)::NSV_LIMA_END_A = 0 ! NSV_LIMA_BEG_A...NSV_LIMA_END_A
INTEGER,DIMENSION(JPMODELMAX)::NSV_LIMA_NC_A = 0       ! First Nc variable
INTEGER,DIMENSION(JPMODELMAX)::NSV_LIMA_NR_A = 0       ! First Nr variable
INTEGER,DIMENSION(JPMODELMAX)::NSV_LIMA_CCN_FREE_A = 0 ! First Free CCN conc.
INTEGER,DIMENSION(JPMODELMAX)::NSV_LIMA_CCN_ACTI_A = 0 ! First Acti. CNN conc.
INTEGER,DIMENSION(JPMODELMAX)::NSV_LIMA_SCAVMASS_A = 0 ! Scavenged mass variable
INTEGER,DIMENSION(JPMODELMAX)::NSV_LIMA_NI_A = 0       ! First Ni var.
INTEGER,DIMENSION(JPMODELMAX)::NSV_LIMA_NS_A = 0       ! First Ns var.
INTEGER,DIMENSION(JPMODELMAX)::NSV_LIMA_NG_A = 0       ! First Ng var.
INTEGER,DIMENSION(JPMODELMAX)::NSV_LIMA_NH_A = 0       ! First Nh var.
INTEGER,DIMENSION(JPMODELMAX)::NSV_LIMA_IFN_FREE_A = 0 ! First Free IFN conc.
INTEGER,DIMENSION(JPMODELMAX)::NSV_LIMA_IFN_NUCL_A = 0 ! First Nucl. IFN conc.
INTEGER,DIMENSION(JPMODELMAX)::NSV_LIMA_IMM_NUCL_A = 0 ! First Nucl. IMM conc.
INTEGER,DIMENSION(JPMODELMAX)::NSV_LIMA_HOM_HAZE_A = 0 ! Hom. freezing of CCN
INTEGER,DIMENSION(JPMODELMAX)::NSV_LIMA_SPRO_A = 0     ! Supersaturation
!
#ifdef MNH_FOREFIRE
INTEGER,DIMENSION(JPMODELMAX)::NSV_FF_A = 0    ! number of ForeFire scalar variables
INTEGER,DIMENSION(JPMODELMAX)::NSV_FFBEG_A = 0 ! with indices in the range :
INTEGER,DIMENSION(JPMODELMAX)::NSV_FFEND_A = 0 ! NSV_FFBEG_A...NSV_FFEND_A
!
#endif
! Blaze smoke indexes
INTEGER,DIMENSION(JPMODELMAX)::NSV_FIRE_A = 0    ! number of Blaze smoke scalar variables
INTEGER,DIMENSION(JPMODELMAX)::NSV_FIREBEG_A = 0 ! with indices in the range :
INTEGER,DIMENSION(JPMODELMAX)::NSV_FIREEND_A = 0 ! NSV_FIREBEG_A...NSV_FIREEND_A
!
INTEGER,DIMENSION(JPMODELMAX)::NSV_SNW_A = 0    ! number of blowing snow scalar
INTEGER,DIMENSION(JPMODELMAX)::NSV_SNWBEG_A = 0 ! with indices in the range :
INTEGER,DIMENSION(JPMODELMAX)::NSV_SNWEND_A = 0 ! NSV_SNWBEG_A...NSV_SNWEND_A
!
!###############################################################################
!
! variables updated for the current model
!
CHARACTER(LEN=NMNHNAMELGTMAX), DIMENSION(:), POINTER :: CSV_CHEM_LIST !Names of all the chemical variables
CHARACTER(LEN=6),              DIMENSION(:), POINTER :: CSV           !Names of the scalar variables

TYPE(tfieldmetadata), DIMENSION(:), POINTER :: TSVLIST !Metadata of all the scalar variables

INTEGER :: NSV         = 0 ! total number of user scalar variables
!
INTEGER :: NSV_CHEM_LIST = 0 ! total number of chemical variables (including dust, salt...)
!
INTEGER :: NSV_USER    = 0 ! number of user scalar variables with indices
                           ! in the range : 1...NSV_USER
INTEGER :: NSV_C2R2    = 0 ! number of liq scalar used in C2R2 and in C3R5
INTEGER :: NSV_C2R2BEG = 0 ! with indices in the range :
INTEGER :: NSV_C2R2END = 0 ! NSV_C2R2BEG...NSV_C2R2END
!
INTEGER :: NSV_C1R3    = 0 ! number of ice scalar used in C3R5
INTEGER :: NSV_C1R3BEG = 0 ! with indices in the range :
INTEGER :: NSV_C1R3END = 0 ! NSV_C1R3BEG...NSV_C1R3END
!
INTEGER :: NSV_ELEC    = 0 ! number of scalar variables used in ELEC
INTEGER :: NSV_ELECBEG = 0 ! with indices in the range :
INTEGER :: NSV_ELECEND = 0 ! NSV_ELECBEG...NSV_ELECEND
!
INTEGER :: NSV_CHEM    = 0 ! number of chemical scalar variables
INTEGER :: NSV_CHEMBEG = 0 ! with indices in the range :
INTEGER :: NSV_CHEMEND = 0 ! NSV_CHEMBEG...NSV_CHEMEND
!
INTEGER :: NSV_CHGS    = 0 ! number of gas-phase chemicals 
INTEGER :: NSV_CHGSBEG = 0 ! with indices in the range :
INTEGER :: NSV_CHGSEND = 0 ! NSV_CHGSBEG...NSV_CHGSEND         
!
INTEGER :: NSV_CHAC    = 0 ! number of aqueous-phase chemicals
INTEGER :: NSV_CHACBEG = 0 ! with indices in the range :
INTEGER :: NSV_CHACEND = 0 ! NSV_CHACBEG...NSV_CHACEND
!
INTEGER :: NSV_CHIC    = 0 ! number of ice-phase chemicals
INTEGER :: NSV_CHICBEG = 0 ! with indices in the range :
INTEGER :: NSV_CHICEND = 0 ! NSV_CHICBEG...NSV_CHICEND
!
INTEGER :: NSV_LG    = 0 ! number of lagrangian
INTEGER :: NSV_LGBEG = 0 ! with indices in the range :
INTEGER :: NSV_LGEND = 0 ! NSV_LGBEG...NSV_LGEND
!
INTEGER :: NSV_LNOX    = 0 ! number of lightning NOx variables
INTEGER :: NSV_LNOXBEG = 0 ! with indices in the range :
INTEGER :: NSV_LNOXEND = 0 ! NSV_LNOXBEG...NSV_LNOXEND
!
INTEGER :: NSV_DST     = 0 ! number of dust scalar variables
INTEGER :: NSV_DSTBEG  = 0 ! with indices in the range :
INTEGER :: NSV_DSTEND  = 0 ! NSV_DSTBEG...NSV_DSTEND

INTEGER :: NSV_SLT     = 0 ! number of sea salt scalar variables
INTEGER :: NSV_SLTBEG  = 0 ! with indices in the range :
INTEGER :: NSV_SLTEND  = 0 ! NSV_SLTBEG...NSV_SLTEND

INTEGER :: NSV_AER     = 0 ! number of aerosol scalar variables
INTEGER :: NSV_AERBEG  = 0 ! with indices in the range :
INTEGER :: NSV_AEREND  = 0 ! NSV_AERBEG...NSV_AEREND

INTEGER :: NSV_DSTDEP  = 0 ! number of aerosol scalar variables
INTEGER :: NSV_DSTDEPBEG  = 0 ! with indices in the range :
INTEGER :: NSV_DSTDEPEND  = 0 ! NSV_AERBEG...NSV_AEREND
!
INTEGER :: NSV_AERDEP  = 0 ! number of aerosol scalar variables
INTEGER :: NSV_AERDEPBEG  = 0 ! with indices in the range :
INTEGER :: NSV_AERDEPEND  = 0 ! NSV_AERBEG...NSV_AEREND

INTEGER :: NSV_SLTDEP  = 0 ! number of aerosol scalar variables
INTEGER :: NSV_SLTDEPBEG  = 0 ! with indices in the range :
INTEGER :: NSV_SLTDEPEND  = 0 ! NSV_AERBEG...NSV_AEREND
!
INTEGER :: NSV_PP    = 0 ! number of passive pollutants       
INTEGER :: NSV_PPBEG = 0 ! with indices in the range :
INTEGER :: NSV_PPEND = 0 ! NSV_PPBEG...NSV_PPEND
!
INTEGER :: NSV_CS    = 0 ! number of condit.samplings         
INTEGER :: NSV_CSBEG = 0 ! with indices in the range :
INTEGER :: NSV_CSEND = 0 ! NSV_CSBEG...NSV_CSEND
!
INTEGER :: NSV_LIMA     ! number of scalar in LIMA
INTEGER :: NSV_LIMA_BEG ! with indices in the range :
INTEGER :: NSV_LIMA_END ! NSV_LIMA_BEG_A...NSV_LIMA_END_A
INTEGER :: NSV_LIMA_NC       !
INTEGER :: NSV_LIMA_NR       !
INTEGER :: NSV_LIMA_CCN_FREE !
INTEGER :: NSV_LIMA_CCN_ACTI !
INTEGER :: NSV_LIMA_SCAVMASS !
INTEGER :: NSV_LIMA_NI       !
INTEGER :: NSV_LIMA_NS       !
INTEGER :: NSV_LIMA_NG       !
INTEGER :: NSV_LIMA_NH       !
INTEGER :: NSV_LIMA_IFN_FREE !
INTEGER :: NSV_LIMA_IFN_NUCL !
INTEGER :: NSV_LIMA_IMM_NUCL !
INTEGER :: NSV_LIMA_HOM_HAZE !
INTEGER :: NSV_LIMA_SPRO     !
!
#ifdef MNH_FOREFIRE
INTEGER :: NSV_FF    = 0 ! number of ForeFire scalar variables
INTEGER :: NSV_FFBEG = 0 ! with indices in the range :
INTEGER :: NSV_FFEND = 0 ! NSV_FFBEG...NSV_FFEND
!
#endif
! Blaze smoke
INTEGER :: NSV_FIRE    = 0 ! number of Blaze smoke scalar variables
INTEGER :: NSV_FIREBEG = 0 ! with indices in the range :
INTEGER :: NSV_FIREEND = 0 ! NSV_FIREBEG...NSV_FIREEND
!
INTEGER :: NSV_SNW     = 0 ! number of blowing snow scalar variables
INTEGER :: NSV_SNWBEG  = 0 ! with indices in the range :
INTEGER :: NSV_SNWEND  = 0 ! NSV_SNWBEG...NSV_SNWEND
!
INTEGER :: NSV_CO2     = 0  ! index for CO2
!
END MODULE MODD_NSV
