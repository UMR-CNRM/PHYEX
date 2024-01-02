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
TYPE NSV_t
!
REAL,DIMENSION(JPSVMAX) :: XSVMIN ! minimum value for SV variables
!
LOGICAL :: LINI_NSV(JPMODELMAX) = .FALSE. ! becomes True when routine INI_NSV is called
!
CHARACTER(LEN=NMNHNAMELGTMAX), DIMENSION(:,:), ALLOCATABLE :: CSV_CHEM_LIST_A !Names of all the chemical variables
CHARACTER(LEN=6),              DIMENSION(:,:), ALLOCATABLE :: CSV_A           !Names of the scalar variables
TYPE(tfieldmetadata), DIMENSION(:,:), ALLOCATABLE :: TSVLIST_A !Metadata of all the scalar variables

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
END TYPE NSV_t
!
TYPE(NSV_t), TARGET, SAVE :: TNSV
!

REAL, POINTER, DIMENSION(:) :: XSVMIN => NULL()

LOGICAL, POINTER :: LINI_NSV(:) => NULL()
!
CHARACTER(LEN=NMNHNAMELGTMAX), DIMENSION(:,:), POINTER :: CSV_CHEM_LIST_A => NULL()
CHARACTER(LEN=6),              DIMENSION(:,:), POINTER :: CSV_A => NULL()
TYPE(tfieldmetadata), DIMENSION(:,:), POINTER :: TSVLIST_A => NULL()

INTEGER, DIMENSION(:), POINTER ::NSV_A => NULL(), &
                                 NSV_CHEM_LIST_A => NULL(), &
                                 NSV_USER_A => NULL(), &
                                 NSV_C2R2_A => NULL(), &
                                 NSV_C2R2BEG_A => NULL(), &
                                 NSV_C2R2END_A => NULL(), &
                                 NSV_C1R3_A => NULL(), &
                                 NSV_C1R3BEG_A => NULL(), &
                                 NSV_C1R3END_A => NULL(), &
                                 NSV_ELEC_A => NULL(), &
                                 NSV_ELECBEG_A => NULL(), &
                                 NSV_ELECEND_A => NULL(), &
                                 NSV_CHEM_A => NULL(), &
                                 NSV_CHEMBEG_A => NULL(), &
                                 NSV_CHEMEND_A => NULL(), &
                                 NSV_CHGS_A => NULL(), &
                                 NSV_CHGSBEG_A => NULL(), &
                                 NSV_CHGSEND_A => NULL(), &
                                 NSV_CHAC_A => NULL(), &
                                 NSV_CHACBEG_A => NULL(), &
                                 NSV_CHACEND_A => NULL(), &
                                 NSV_CHIC_A => NULL(), &
                                 NSV_CHICBEG_A => NULL(), &
                                 NSV_CHICEND_A => NULL(), &
                                 NSV_LG_A => NULL(), &
                                 NSV_LGBEG_A => NULL(), &
                                 NSV_LGEND_A => NULL(), &
                                 NSV_LNOX_A => NULL(), &
                                 NSV_LNOXBEG_A => NULL(), &
                                 NSV_LNOXEND_A => NULL(), &
                                 NSV_DST_A => NULL(), &
                                 NSV_DSTBEG_A => NULL(), &
                                 NSV_DSTEND_A => NULL(), &
                                 NSV_SLT_A => NULL(), &
                                 NSV_SLTBEG_A => NULL(), &
                                 NSV_SLTEND_A => NULL(), &
                                 NSV_AER_A => NULL(), &
                                 NSV_AERBEG_A => NULL(), &
                                 NSV_AEREND_A => NULL(), &
                                 NSV_DSTDEP_A    => NULL(), &
                                 NSV_DSTDEPBEG_A => NULL(), &
                                 NSV_DSTDEPEND_A => NULL(), &
                                 NSV_AERDEP_A    => NULL(), &
                                 NSV_AERDEPBEG_A => NULL(), &
                                 NSV_AERDEPEND_A => NULL(), &
                                 NSV_SLTDEP_A    => NULL(), &
                                 NSV_SLTDEPBEG_A => NULL(), &
                                 NSV_SLTDEPEND_A => NULL(), &
                                 NSV_PP_A => NULL(), &
                                 NSV_PPBEG_A => NULL(), &
                                 NSV_PPEND_A => NULL(), &
                                 NSV_CS_A => NULL(), &
                                 NSV_CSBEG_A => NULL(), &
                                 NSV_CSEND_A => NULL(), &
                                 NSV_LIMA_A => NULL(), &
                                 NSV_LIMA_BEG_A => NULL(), &
                                 NSV_LIMA_END_A => NULL(), &
                                 NSV_LIMA_NC_A => NULL(), &
                                 NSV_LIMA_NR_A => NULL(), &
                                 NSV_LIMA_CCN_FREE_A => NULL(), &
                                 NSV_LIMA_CCN_ACTI_A => NULL(), &
                                 NSV_LIMA_SCAVMASS_A => NULL(), &
                                 NSV_LIMA_NI_A => NULL(), &
                                 NSV_LIMA_NS_A => NULL(), &
                                 NSV_LIMA_NG_A => NULL(), &
                                 NSV_LIMA_NH_A => NULL(), &
                                 NSV_LIMA_IFN_FREE_A => NULL(), &
                                 NSV_LIMA_IFN_NUCL_A => NULL(), &
                                 NSV_LIMA_IMM_NUCL_A => NULL(), &
                                 NSV_LIMA_HOM_HAZE_A => NULL(), &
                                 NSV_LIMA_SPRO_A => NULL(), &
#ifdef MNH_FOREFIRE
                                 NSV_FF_A => NULL(), &
                                 NSV_FFBEG_A => NULL(), &
                                 NSV_FFEND_A => NULL(), &
#endif
                                 NSV_FIRE_A => NULL(), &
                                 NSV_FIREBEG_A => NULL(), &
                                 NSV_FIREEND_A => NULL(), &
                                 NSV_SNW_A => NULL(), &
                                 NSV_SNWBEG_A => NULL(), &
                                 NSV_SNWEND_A => NULL()

CHARACTER(LEN=NMNHNAMELGTMAX), DIMENSION(:), POINTER :: CSV_CHEM_LIST => NULL()
CHARACTER(LEN=6),              DIMENSION(:), POINTER :: CSV           => NULL()

TYPE(tfieldmetadata), DIMENSION(:), POINTER :: TSVLIST => NULL()

INTEGER, POINTER :: NSV         => NULL(), &
                    NSV_CHEM_LIST => NULL(), &
                    NSV_USER    => NULL(), &
                    NSV_C2R2    => NULL(), &
                    NSV_C2R2BEG => NULL(), &
                    NSV_C2R2END => NULL(), &
                    NSV_C1R3    => NULL(), &
                    NSV_C1R3BEG => NULL(), &
                    NSV_C1R3END => NULL(), &
                    NSV_ELEC    => NULL(), &
                    NSV_ELECBEG => NULL(), &
                    NSV_ELECEND => NULL(), &
                    NSV_CHEM    => NULL(), &
                    NSV_CHEMBEG => NULL(), &
                    NSV_CHEMEND => NULL(), &
                    NSV_CHGS    => NULL(), &
                    NSV_CHGSBEG => NULL(), &
                    NSV_CHGSEND => NULL(), &
                    NSV_CHAC    => NULL(), &
                    NSV_CHACBEG => NULL(), &
                    NSV_CHACEND => NULL(), &
                    NSV_CHIC    => NULL(), &
                    NSV_CHICBEG => NULL(), &
                    NSV_CHICEND => NULL(), &
                    NSV_LG    => NULL(), &
                    NSV_LGBEG => NULL(), &
                    NSV_LGEND => NULL(), &
                    NSV_LNOX    => NULL(), &
                    NSV_LNOXBEG => NULL(), &
                    NSV_LNOXEND => NULL(), &
                    NSV_DST     => NULL(), &
                    NSV_DSTBEG  => NULL(), &
                    NSV_DSTEND  => NULL(), &
                    NSV_SLT     => NULL(), &
                    NSV_SLTBEG  => NULL(), &
                    NSV_SLTEND  => NULL(), &
                    NSV_AER     => NULL(), &
                    NSV_AERBEG  => NULL(), &
                    NSV_AEREND  => NULL(), &
                    NSV_DSTDEP  => NULL(), &
                    NSV_DSTDEPBEG  => NULL(), &
                    NSV_DSTDEPEND  => NULL(), &
                    NSV_AERDEP  => NULL(), &
                    NSV_AERDEPBEG  => NULL(), &
                    NSV_AERDEPEND  => NULL(), &
                    NSV_SLTDEP  => NULL(), &
                    NSV_SLTDEPBEG  => NULL(), &
                    NSV_SLTDEPEND  => NULL(), &
                    NSV_PP    => NULL(), &
                    NSV_PPBEG => NULL(), &
                    NSV_PPEND => NULL(), &
                    NSV_CS    => NULL(), &
                    NSV_CSBEG => NULL(), &
                    NSV_CSEND => NULL(), &
                    NSV_LIMA  => NULL(), &
                    NSV_LIMA_BEG => NULL(), &
                    NSV_LIMA_END => NULL(), &
                    NSV_LIMA_NC       => NULL(), &
                    NSV_LIMA_NR       => NULL(), &
                    NSV_LIMA_CCN_FREE => NULL(), &
                    NSV_LIMA_CCN_ACTI => NULL(), &
                    NSV_LIMA_SCAVMASS => NULL(), &
                    NSV_LIMA_NI       => NULL(), &
                    NSV_LIMA_NS       => NULL(), &
                    NSV_LIMA_NG       => NULL(), &
                    NSV_LIMA_NH       => NULL(), &
                    NSV_LIMA_IFN_FREE => NULL(), &
                    NSV_LIMA_IFN_NUCL => NULL(), &
                    NSV_LIMA_IMM_NUCL => NULL(), &
                    NSV_LIMA_HOM_HAZE => NULL(), &
                    NSV_LIMA_SPRO     => NULL(), &
#ifdef MNH_FOREFIRE
                    NSV_FF    => NULL(), &
                    NSV_FFBEG => NULL(), &
                    NSV_FFEND => NULL(), &
#endif
                    NSV_FIRE    => NULL(), &
                    NSV_FIREBEG => NULL(), &
                    NSV_FIREEND => NULL(), &
                    NSV_SNW     => NULL(), &
                    NSV_SNWBEG  => NULL(), &
                    NSV_SNWEND  => NULL(), &
                    NSV_CO2     => NULL()
!
CONTAINS
!
SUBROUTINE NSV_ASSOCIATE()
IMPLICIT NONE

IF(.NOT. ASSOCIATED(NSV)) THEN
  XSVMIN              =>  TNSV%XSVMIN                
  LINI_NSV            =>  TNSV%LINI_NSV

  NSV_A               =>  TNSV%NSV_A 
  NSV_CHEM_LIST_A     =>  TNSV%NSV_CHEM_LIST_A 
  NSV_USER_A          =>  TNSV%NSV_USER_A 
  NSV_C2R2_A          =>  TNSV%NSV_C2R2_A 
  NSV_C2R2BEG_A       =>  TNSV%NSV_C2R2BEG_A 
  NSV_C2R2END_A       =>  TNSV%NSV_C2R2END_A 
  NSV_C1R3_A          =>  TNSV%NSV_C1R3_A 
  NSV_C1R3BEG_A       =>  TNSV%NSV_C1R3BEG_A 
  NSV_C1R3END_A       =>  TNSV%NSV_C1R3END_A 
  NSV_ELEC_A          =>  TNSV%NSV_ELEC_A 
  NSV_ELECBEG_A       =>  TNSV%NSV_ELECBEG_A 
  NSV_ELECEND_A       =>  TNSV%NSV_ELECEND_A 
  NSV_CHEM_A          =>  TNSV%NSV_CHEM_A 
  NSV_CHEMBEG_A       =>  TNSV%NSV_CHEMBEG_A 
  NSV_CHEMEND_A       =>  TNSV%NSV_CHEMEND_A 
  NSV_CHGS_A          =>  TNSV%NSV_CHGS_A 
  NSV_CHGSBEG_A       =>  TNSV%NSV_CHGSBEG_A 
  NSV_CHGSEND_A       =>  TNSV%NSV_CHGSEND_A 
  NSV_CHAC_A          =>  TNSV%NSV_CHAC_A 
  NSV_CHACBEG_A       =>  TNSV%NSV_CHACBEG_A 
  NSV_CHACEND_A       =>  TNSV%NSV_CHACEND_A 
  NSV_CHIC_A          =>  TNSV%NSV_CHIC_A 
  NSV_CHICBEG_A       =>  TNSV%NSV_CHICBEG_A 
  NSV_CHICEND_A       =>  TNSV%NSV_CHICEND_A 
  NSV_LG_A            =>  TNSV%NSV_LG_A 
  NSV_LGBEG_A         =>  TNSV%NSV_LGBEG_A 
  NSV_LGEND_A         =>  TNSV%NSV_LGEND_A 
  NSV_LNOX_A          =>  TNSV%NSV_LNOX_A 
  NSV_LNOXBEG_A       =>  TNSV%NSV_LNOXBEG_A 
  NSV_LNOXEND_A       =>  TNSV%NSV_LNOXEND_A 
  NSV_DST_A           =>  TNSV%NSV_DST_A 
  NSV_DSTBEG_A        =>  TNSV%NSV_DSTBEG_A 
  NSV_DSTEND_A        =>  TNSV%NSV_DSTEND_A 
  NSV_SLT_A           =>  TNSV%NSV_SLT_A 
  NSV_SLTBEG_A        =>  TNSV%NSV_SLTBEG_A 
  NSV_SLTEND_A        =>  TNSV%NSV_SLTEND_A 
  NSV_AER_A           =>  TNSV%NSV_AER_A 
  NSV_AERBEG_A        =>  TNSV%NSV_AERBEG_A 
  NSV_AEREND_A        =>  TNSV%NSV_AEREND_A 
  NSV_DSTDEP_A        =>  TNSV%NSV_DSTDEP_A    
  NSV_DSTDEPBEG_A     =>  TNSV%NSV_DSTDEPBEG_A 
  NSV_DSTDEPEND_A     =>  TNSV%NSV_DSTDEPEND_A 
  NSV_AERDEP_A        =>  TNSV%NSV_AERDEP_A    
  NSV_AERDEPBEG_A     =>  TNSV%NSV_AERDEPBEG_A 
  NSV_AERDEPEND_A     =>  TNSV%NSV_AERDEPEND_A 
  NSV_SLTDEP_A        =>  TNSV%NSV_SLTDEP_A    
  NSV_SLTDEPBEG_A     =>  TNSV%NSV_SLTDEPBEG_A 
  NSV_SLTDEPEND_A     =>  TNSV%NSV_SLTDEPEND_A 
  NSV_PP_A            =>  TNSV%NSV_PP_A 
  NSV_PPBEG_A         =>  TNSV%NSV_PPBEG_A 
  NSV_PPEND_A         =>  TNSV%NSV_PPEND_A 
  NSV_CS_A            =>  TNSV%NSV_CS_A 
  NSV_CSBEG_A         =>  TNSV%NSV_CSBEG_A 
  NSV_CSEND_A         =>  TNSV%NSV_CSEND_A 
  NSV_LIMA_A          =>  TNSV%NSV_LIMA_A 
  NSV_LIMA_BEG_A      =>  TNSV%NSV_LIMA_BEG_A 
  NSV_LIMA_END_A      =>  TNSV%NSV_LIMA_END_A 
  NSV_LIMA_NC_A       =>  TNSV%NSV_LIMA_NC_A 
  NSV_LIMA_NR_A       =>  TNSV%NSV_LIMA_NR_A 
  NSV_LIMA_CCN_FREE_A =>  TNSV%NSV_LIMA_CCN_FREE_A 
  NSV_LIMA_CCN_ACTI_A =>  TNSV%NSV_LIMA_CCN_ACTI_A 
  NSV_LIMA_SCAVMASS_A =>  TNSV%NSV_LIMA_SCAVMASS_A 
  NSV_LIMA_NI_A       =>  TNSV%NSV_LIMA_NI_A 
  NSV_LIMA_NS_A       =>  TNSV%NSV_LIMA_NS_A 
  NSV_LIMA_NG_A       =>  TNSV%NSV_LIMA_NG_A 
  NSV_LIMA_NH_A       =>  TNSV%NSV_LIMA_NH_A 
  NSV_LIMA_IFN_FREE_A =>  TNSV%NSV_LIMA_IFN_FREE_A 
  NSV_LIMA_IFN_NUCL_A =>  TNSV%NSV_LIMA_IFN_NUCL_A 
  NSV_LIMA_IMM_NUCL_A =>  TNSV%NSV_LIMA_IMM_NUCL_A 
  NSV_LIMA_HOM_HAZE_A =>  TNSV%NSV_LIMA_HOM_HAZE_A 
  NSV_LIMA_SPRO_A     =>  TNSV%NSV_LIMA_SPRO_A 
#ifdef MNH_FOREFIRE
  NSV_FF_A            =>  TNSV%NSV_FF_A 
  NSV_FFBEG_A         =>  TNSV%NSV_FFBEG_A 
  NSV_FFEND_A         =>  TNSV%NSV_FFEND_A 
#endif
  NSV_FIRE_A          =>  TNSV%NSV_FIRE_A 
  NSV_FIREBEG_A       =>  TNSV%NSV_FIREBEG_A 
  NSV_FIREEND_A       =>  TNSV%NSV_FIREEND_A 
  NSV_SNW_A           =>  TNSV%NSV_SNW_A 
  NSV_SNWBEG_A        =>  TNSV%NSV_SNWBEG_A 
  NSV_SNWEND_A        =>  TNSV%NSV_SNWEND_A 

  CSV_CHEM_LIST       =>  TNSV%CSV_CHEM_LIST 
  CSV                 =>  TNSV%CSV           
  TSVLIST             =>  TNSV%TSVLIST 

  NSV                 =>  TNSV%NSV         
  NSV_CHEM_LIST       =>  TNSV%NSV_CHEM_LIST 
  NSV_USER            =>  TNSV%NSV_USER    
  NSV_C2R2            =>  TNSV%NSV_C2R2    
  NSV_C2R2BEG         =>  TNSV%NSV_C2R2BEG 
  NSV_C2R2END         =>  TNSV%NSV_C2R2END 
  NSV_C1R3            =>  TNSV%NSV_C1R3    
  NSV_C1R3BEG         =>  TNSV%NSV_C1R3BEG 
  NSV_C1R3END         =>  TNSV%NSV_C1R3END 
  NSV_ELEC            =>  TNSV%NSV_ELEC    
  NSV_ELECBEG         =>  TNSV%NSV_ELECBEG 
  NSV_ELECEND         =>  TNSV%NSV_ELECEND 
  NSV_CHEM            =>  TNSV%NSV_CHEM    
  NSV_CHEMBEG         =>  TNSV%NSV_CHEMBEG 
  NSV_CHEMEND         =>  TNSV%NSV_CHEMEND 
  NSV_CHGS            =>  TNSV%NSV_CHGS    
  NSV_CHGSBEG         =>  TNSV%NSV_CHGSBEG 
  NSV_CHGSEND         =>  TNSV%NSV_CHGSEND 
  NSV_CHAC            =>  TNSV%NSV_CHAC    
  NSV_CHACBEG         =>  TNSV%NSV_CHACBEG 
  NSV_CHACEND         =>  TNSV%NSV_CHACEND 
  NSV_CHIC            =>  TNSV%NSV_CHIC    
  NSV_CHICBEG         =>  TNSV%NSV_CHICBEG 
  NSV_CHICEND         =>  TNSV%NSV_CHICEND 
  NSV_LG              =>  TNSV%NSV_LG    
  NSV_LGBEG           =>  TNSV%NSV_LGBEG 
  NSV_LGEND           =>  TNSV%NSV_LGEND 
  NSV_LNOX            =>  TNSV%NSV_LNOX    
  NSV_LNOXBEG         =>  TNSV%NSV_LNOXBEG 
  NSV_LNOXEND         =>  TNSV%NSV_LNOXEND 
  NSV_DST             =>  TNSV%NSV_DST     
  NSV_DSTBEG          =>  TNSV%NSV_DSTBEG  
  NSV_DSTEND          =>  TNSV%NSV_DSTEND  
  NSV_SLT             =>  TNSV%NSV_SLT     
  NSV_SLTBEG          =>  TNSV%NSV_SLTBEG  
  NSV_SLTEND          =>  TNSV%NSV_SLTEND  
  NSV_AER             =>  TNSV%NSV_AER     
  NSV_AERBEG          =>  TNSV%NSV_AERBEG  
  NSV_AEREND          =>  TNSV%NSV_AEREND  
  NSV_DSTDEP          =>  TNSV%NSV_DSTDEP  
  NSV_DSTDEPBEG       =>  TNSV%NSV_DSTDEPBEG  
  NSV_DSTDEPEND       =>  TNSV%NSV_DSTDEPEND  
  NSV_AERDEP          =>  TNSV%NSV_AERDEP  
  NSV_AERDEPBEG       =>  TNSV%NSV_AERDEPBEG  
  NSV_AERDEPEND       =>  TNSV%NSV_AERDEPEND  
  NSV_SLTDEP          =>  TNSV%NSV_SLTDEP  
  NSV_SLTDEPBEG       =>  TNSV%NSV_SLTDEPBEG  
  NSV_SLTDEPEND       =>  TNSV%NSV_SLTDEPEND  
  NSV_PP              =>  TNSV%NSV_PP    
  NSV_PPBEG           =>  TNSV%NSV_PPBEG 
  NSV_PPEND           =>  TNSV%NSV_PPEND 
  NSV_CS              =>  TNSV%NSV_CS    
  NSV_CSBEG           =>  TNSV%NSV_CSBEG 
  NSV_CSEND           =>  TNSV%NSV_CSEND 
  NSV_LIMA            =>  TNSV%NSV_LIMA  
  NSV_LIMA_BEG        =>  TNSV%NSV_LIMA_BEG 
  NSV_LIMA_END        =>  TNSV%NSV_LIMA_END 
  NSV_LIMA_NC         =>  TNSV%NSV_LIMA_NC       
  NSV_LIMA_NR         =>  TNSV%NSV_LIMA_NR       
  NSV_LIMA_CCN_FREE   =>  TNSV%NSV_LIMA_CCN_FREE 
  NSV_LIMA_CCN_ACTI   =>  TNSV%NSV_LIMA_CCN_ACTI 
  NSV_LIMA_SCAVMASS   =>  TNSV%NSV_LIMA_SCAVMASS 
  NSV_LIMA_NI         =>  TNSV%NSV_LIMA_NI       
  NSV_LIMA_NS         =>  TNSV%NSV_LIMA_NS       
  NSV_LIMA_NG         =>  TNSV%NSV_LIMA_NG       
  NSV_LIMA_NH         =>  TNSV%NSV_LIMA_NH       
  NSV_LIMA_IFN_FREE   =>  TNSV%NSV_LIMA_IFN_FREE 
  NSV_LIMA_IFN_NUCL   =>  TNSV%NSV_LIMA_IFN_NUCL 
  NSV_LIMA_IMM_NUCL   =>  TNSV%NSV_LIMA_IMM_NUCL 
  NSV_LIMA_HOM_HAZE   =>  TNSV%NSV_LIMA_HOM_HAZE 
  NSV_LIMA_SPRO       =>  TNSV%NSV_LIMA_SPRO     
#ifdef MNH_FOREFIRE
  NSV_FF              =>  TNSV%NSV_FF    
  NSV_FFBEG           =>  TNSV%NSV_FFBEG 
  NSV_FFEND           =>  TNSV%NSV_FFEND 
#endif
  NSV_FIRE            =>  TNSV%NSV_FIRE    
  NSV_FIREBEG         =>  TNSV%NSV_FIREBEG 
  NSV_FIREEND         =>  TNSV%NSV_FIREEND 
  NSV_SNW             =>  TNSV%NSV_SNW     
  NSV_SNWBEG          =>  TNSV%NSV_SNWBEG  
  NSV_SNWEND          =>  TNSV%NSV_SNWEND  
  NSV_CO2             =>  TNSV%NSV_CO2     
ENDIF
!
END SUBROUTINE NSV_ASSOCIATE
!
END MODULE MODD_NSV
