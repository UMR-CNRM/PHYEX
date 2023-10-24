!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ####################
      MODULE MODD_ELEC_n
!     ####################
!
!!****  *MODD_ELEC$n* - declaration of electric fields
!!
!!    PURPOSE
!!    -------
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!          
!!    AUTHOR
!!    ------
!!
!!    MODIFICATIONS
!!    -------------
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_PARAMETERS, ONLY: JPMODELMAX
IMPLICIT NONE

TYPE ELEC_t
!
!  REAL, DIMENSION(:,:,:), POINTER :: XNI_SDRYG=>NULL(), XNI_IDRYG=>NULL(),  &
!     XNI_IAGGS=>NULL(),                 & 
!     XEFIELDU=>NULL(), & ! The 3 components of the electric field
!     XEFIELDV=>NULL(), XEFIELDW=>NULL(), &
  REAL, DIMENSION(:,:,:), POINTER ::      XESOURCEFW=>NULL(),  & ! Fair weather electric charge (C m^-3)
!     XIND_RATE=>NULL(), XEW=>NULL(),  & 
     XEW=>NULL(),  & 
     XIONSOURCEFW =>NULL(), & ! Fair weather ionic source
                              !  (ion pairs m-3 s-1) hold constant in time
     XCION_POS_FW =>NULL(), XCION_NEG_FW =>NULL(), &  !Positive and Negative ion mixing ratio
     XMOBIL_POS =>NULL(), XMOBIL_NEG=>NULL() ! m2/V/s
!
!  Parameters for flat lapalcian operator to solve the electric field
!            (see MODD_DYN_n)
  REAL, DIMENSION(:), POINTER :: XRHOM_E =>NULL(), XAF_E =>NULL(), XCF_E =>NULL()
  REAL, DIMENSION(:,:,:), POINTER :: XBFY_E =>NULL(), &
                                     XBFB_E =>NULL(), XBF_SXP2_YP1_Z_E =>NULL() 
                               ! Z_Splitting
     
!
END TYPE ELEC_t

TYPE(ELEC_t), DIMENSION(JPMODELMAX), TARGET, SAVE :: ELEC_MODEL

REAL, DIMENSION(:,:,:), POINTER :: XNI_SDRYG=>NULL(), XNI_IDRYG=>NULL(),  &
                 XNI_IAGGS=>NULL(), XEFIELDU=>NULL(),                     &
                 XESOURCEFW=>NULL(), XEFIELDV=>NULL(), XEFIELDW=>NULL(),  &
                 XIND_RATE=>NULL(), XIONSOURCEFW =>NULL(), XEW=>NULL(),   &                
                 XCION_POS_FW =>NULL(), XCION_NEG_FW =>NULL(),            &  
                 XMOBIL_POS =>NULL(), XMOBIL_NEG=>NULL(),  XBFY_E =>NULL(), &
                 XBFB_E =>NULL(), XBF_SXP2_YP1_Z_E =>NULL()
REAL, DIMENSION(:), POINTER :: XRHOM_E =>NULL(), XAF_E =>NULL(), XCF_E =>NULL()

CONTAINS

SUBROUTINE ELEC_GOTO_MODEL(KFROM, KTO)
INTEGER, INTENT(IN) :: KFROM, KTO
!
! Save current state for allocated arrays
!ELEC_MODEL(KFROM)%XNI_SDRYG=>XNI_SDRYG !Done in FIELDLIST_GOTO_MODEL
!ELEC_MODEL(KFROM)%XNI_IDRYG=>XNI_IDRYG !Done in FIELDLIST_GOTO_MODEL
!ELEC_MODEL(KFROM)%XNI_IAGGS=>XNI_IAGGS !Done in FIELDLIST_GOTO_MODEL
!ELEC_MODEL(KFROM)%XIND_RATE=>XIND_RATE !Done in FIELDLIST_GOTO_MODEL
ELEC_MODEL(KFROM)%XEW=>XEW
!ELEC_MODEL(KFROM)%XEFIELDU=>XEFIELDU !Done in FIELDLIST_GOTO_MODEL
!ELEC_MODEL(KFROM)%XEFIELDV=>XEFIELDV !Done in FIELDLIST_GOTO_MODEL
!ELEC_MODEL(KFROM)%XEFIELDW=>XEFIELDW !Done in FIELDLIST_GOTO_MODEL
ELEC_MODEL(KFROM)%XESOURCEFW=>XESOURCEFW
ELEC_MODEL(KFROM)%XIONSOURCEFW=>XIONSOURCEFW
ELEC_MODEL(KFROM)%XCION_POS_FW=>XCION_POS_FW
ELEC_MODEL(KFROM)%XCION_NEG_FW=>XCION_NEG_FW
ELEC_MODEL(KFROM)%XMOBIL_POS=>XMOBIL_POS  
ELEC_MODEL(KFROM)%XMOBIL_NEG=>XMOBIL_NEG  
ELEC_MODEL(KFROM)%XBFY_E=>XBFY_E
ELEC_MODEL(KFROM)%XBFB_E=>XBFB_E
ELEC_MODEL(KFROM)%XBF_SXP2_YP1_Z_E=>XBF_SXP2_YP1_Z_E
ELEC_MODEL(KFROM)%XRHOM_E=>XRHOM_E
ELEC_MODEL(KFROM)%XAF_E=>XAF_E
ELEC_MODEL(KFROM)%XCF_E=>XCF_E
!
! Current model is set to model KTO
!XNI_SDRYG=>ELEC_MODEL(KTO)%XNI_SDRYG !Done in FIELDLIST_GOTO_MODEL
!XNI_IDRYG=>ELEC_MODEL(KTO)%XNI_IDRYG !Done in FIELDLIST_GOTO_MODEL
!XNI_IAGGS=>ELEC_MODEL(KTO)%XNI_IAGGS !Done in FIELDLIST_GOTO_MODEL
!XIND_RATE=>ELEC_MODEL(KTO)%XIND_RATE !Done in FIELDLIST_GOTO_MODEL
XEW=>ELEC_MODEL(KTO)%XEW
!XEFIELDU=>ELEC_MODEL(KTO)%XEFIELDU !Done in FIELDLIST_GOTO_MODEL
!XEFIELDV=>ELEC_MODEL(KTO)%XEFIELDV !Done in FIELDLIST_GOTO_MODEL
!XEFIELDW=>ELEC_MODEL(KTO)%XEFIELDW !Done in FIELDLIST_GOTO_MODEL
XESOURCEFW=>ELEC_MODEL(KTO)%XESOURCEFW
XIONSOURCEFW=>ELEC_MODEL(KTO)%XIONSOURCEFW
XCION_POS_FW=>ELEC_MODEL(KTO)%XCION_POS_FW
XCION_NEG_FW=>ELEC_MODEL(KTO)%XCION_NEG_FW
XMOBIL_POS=>ELEC_MODEL(KTO)%XMOBIL_POS
XMOBIL_NEG=>ELEC_MODEL(KTO)%XMOBIL_NEG
XBFY_E=>ELEC_MODEL(KTO)%XBFY_E
XBFB_E=>ELEC_MODEL(KTO)%XBFB_E
XBF_SXP2_YP1_Z_E=>ELEC_MODEL(KTO)%XBF_SXP2_YP1_Z_E
XRHOM_E=>ELEC_MODEL(KTO)%XRHOM_E
XAF_E=>ELEC_MODEL(KTO)%XAF_E
XCF_E=>ELEC_MODEL(KTO)%XCF_E
END SUBROUTINE ELEC_GOTO_MODEL

END MODULE MODD_ELEC_n
