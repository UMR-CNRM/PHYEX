!     ######spl
      MODULE MODN_CH_ORILAM
!!    #####################
!!
!!*** *MODN_CH_ORILAM*
!!
!!    PURPOSE
!!    -------
!       Namelist for ORILAM aerosol scheme parameters 
!!
!!**  AUTHOR
!!    ------
!!    P. Tulet      *CNRM*
!
!!    MODIFICATIONS
!!    -------------
!!    Original 24/02/05
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
USE MODD_CH_AEROSOL, ONLY: LORILAM, XN0IMIN, XN0JMIN, LSEDIMAERO, LAERINIT, &
                          LHETEROSO4, CNUCLEATION, XINISIGI, XINISIGJ,  & 
                          XINIRADIUSI, XINIRADIUSJ, LVARSIGI,&
                          LVARSIGJ, CMINERAL, CORGANIC,&
                          XSIGIMIN, XSIGIMAX,XSIGJMIN, XSIGJMAX,  & 
                          XCOEFRADIMAX, XCOEFRADIMIN, XCOEFRADJMAX, XCOEFRADJMIN,&
                          CRGUNIT, LRGFIX
!!
!-----------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!        -----------------
IMPLICIT NONE
SAVE
NAMELIST /NAM_CH_ORILAM/  LORILAM, XN0IMIN, XN0JMIN, LSEDIMAERO, LAERINIT, &
                          LHETEROSO4, CNUCLEATION, XINISIGI, XINISIGJ,  & 
                          XINIRADIUSI, XINIRADIUSJ, LVARSIGI,&
                          LVARSIGJ, CMINERAL, CORGANIC, &
                          XSIGIMIN, XSIGIMAX,XSIGJMIN, XSIGJMAX,  & 
                          XCOEFRADIMAX, XCOEFRADIMIN, XCOEFRADJMAX, XCOEFRADJMIN,&
                          CRGUNIT, LRGFIX

!
END MODULE MODN_CH_ORILAM
