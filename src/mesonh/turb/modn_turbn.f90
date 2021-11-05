!MNH_LIC Copyright 1995-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ###################
      MODULE MODN_TURB_n
!     ###################
!
!!****  *MODN_TURB$n* - declaration of namelist NAM_TURBn
!!
!!    PURPOSE
!!    -------
!       The purpose of this module is to specify  the namelist NAM_TURBn
!     which concern the parameters of the turbulence scheme for one nested
!     model.   
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_TURB$n : contains declaration of turbulence scheme
!!    variables entering by a namelist
!! 
!!          XIMPL,CTURBLEN,CTURBDIM,LTURB_FLX
!!          LTURB_DIAG,LSUBG_COND,LTGT_FLX
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (module MODD_TURBn)
!!          
!!    AUTHOR
!!    ------
!!	    J. Cuxart and J. Stein     * I.N.M. and Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    January 9, 1995               
!!      J.Cuxart    February 15, 1995 add the switches for diagnostic storages
!!      J. Stein    June 14, 1995   add the subgrid condensation switch
!!      J. Stein    October, 1999   add the tangential fluxes switch
!!      M. Tomasini Jul  05, 2001   add the subgrid autoconversion
!!      P. Bechtold Feb 11, 2002    add switch for Sigma_s computation
!!      P. Jabouille Apr 4, 2002    add switch for Sigma_s convection
!!      V. Masson    Nov 13 2002    add switch for SBL lengths
!!      D. Ricard    May, 2021      add switch for Leonard Terms
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_TURB_n, ONLY: &
         XIMPL_n => XIMPL, &
         XKEMIN_n => XKEMIN, &
         XCEDIS_n => XCEDIS, &
         XCADAP_n => XCADAP, &
         CTURBLEN_n => CTURBLEN, &
         CTURBDIM_n => CTURBDIM, &
         LTURB_FLX_n => LTURB_FLX, &
         LTURB_DIAG_n => LTURB_DIAG, &
         LSUBG_COND_n => LSUBG_COND, &
         LSIGMAS_n => LSIGMAS, &
         LSIG_CONV_n => LSIG_CONV, &
         LRMC01_n => LRMC01, &
         CTOM_n => CTOM, &
         CSUBG_AUCV_n => CSUBG_AUCV, &
         VSIGQSAT_n => VSIGQSAT, &
         CSUBG_AUCV_RI_n => CSUBG_AUCV_RI, &
         CCONDENS_n => CCONDENS, &
         CLAMBDA3_n => CLAMBDA3, &
         CSUBG_MF_PDF_n => CSUBG_MF_PDF, &
         LHGRAD_n => LHGRAD, &
         XCOEFHGRADTHL_n => XCOEFHGRADTHL, &
         XCOEFHGRADRM_n => XCOEFHGRADRM, &
         XALTHGRAD_n => XALTHGRAD, &
         XCLDTHOLD_n => XCLDTHOLD
!
IMPLICIT NONE
!
REAL,SAVE  :: XIMPL
REAL,SAVE :: XKEMIN
REAL,SAVE :: XCEDIS
REAL,SAVE :: XCADAP
CHARACTER (LEN=4),SAVE  :: CTURBLEN
CHARACTER (LEN=4),SAVE  :: CTURBDIM
LOGICAL,SAVE  :: LTURB_FLX
LOGICAL,SAVE  :: LTURB_DIAG
LOGICAL,SAVE  :: LSUBG_COND
LOGICAL,SAVE  :: LSIGMAS
LOGICAL,SAVE  :: LSIG_CONV
LOGICAL,SAVE  :: LRMC01
CHARACTER (LEN=4),SAVE  :: CTOM
CHARACTER (LEN=4),SAVE  :: CSUBG_AUCV
CHARACTER (LEN=80),SAVE  :: CSUBG_AUCV_RI
CHARACTER (LEN=80),SAVE  :: CCONDENS
CHARACTER (LEN=4),SAVE  :: CLAMBDA3
CHARACTER (LEN=80),SAVE  :: CSUBG_MF_PDF
REAL,SAVE :: VSIGQSAT
LOGICAL,SAVE  :: LHGRAD
REAL,SAVE :: XCOEFHGRADTHL
REAL,SAVE :: XCOEFHGRADRM
REAL,SAVE :: XALTHGRAD
REAL,SAVE :: XCLDTHOLD
!
NAMELIST/NAM_TURBn/XIMPL,CTURBLEN,CTURBDIM,LTURB_FLX,LTURB_DIAG,  &
                   LSUBG_COND,LSIGMAS,LSIG_CONV,LRMC01,CTOM,CSUBG_AUCV,&
                   XKEMIN,VSIGQSAT,XCEDIS,XCADAP,CSUBG_AUCV_RI,CCONDENS,&
                   CLAMBDA3,CSUBG_MF_PDF,LHGRAD,XCOEFHGRADTHL, XCOEFHGRADRM, &
                   XALTHGRAD, XCLDTHOLD

!
CONTAINS
!
SUBROUTINE INIT_NAM_TURBn
  XIMPL = XIMPL_n
  XKEMIN = XKEMIN_n
  XCEDIS = XCEDIS_n
  XCADAP = XCADAP_n
  CTURBLEN = CTURBLEN_n
  CTURBDIM = CTURBDIM_n
  LTURB_FLX = LTURB_FLX_n
  LTURB_DIAG = LTURB_DIAG_n
  LSUBG_COND = LSUBG_COND_n
  LSIGMAS = LSIGMAS_n
  LSIG_CONV = LSIG_CONV_n
  LRMC01 = LRMC01_n
  CTOM = CTOM_n
  CSUBG_AUCV = CSUBG_AUCV_n
  VSIGQSAT = VSIGQSAT_n  
  CSUBG_AUCV_RI = CSUBG_AUCV_RI_n
  CCONDENS = CCONDENS_n
  CLAMBDA3 = CLAMBDA3_n
  CSUBG_MF_PDF = CSUBG_MF_PDF_n
  LHGRAD = LHGRAD_n
  XCOEFHGRADTHL = XCOEFHGRADTHL_n
  XCOEFHGRADRM = XCOEFHGRADRM_n
  XALTHGRAD = XALTHGRAD_n
  XCLDTHOLD = XCLDTHOLD_n
END SUBROUTINE INIT_NAM_TURBn

SUBROUTINE UPDATE_NAM_TURBn
  XIMPL_n = XIMPL
  XKEMIN_n = XKEMIN
  XCEDIS_n = XCEDIS
  XCADAP_n = XCADAP
  CTURBLEN_n = CTURBLEN
  CTURBDIM_n = CTURBDIM
  LTURB_FLX_n = LTURB_FLX
  LTURB_DIAG_n = LTURB_DIAG
  LSUBG_COND_n = LSUBG_COND
  LSIGMAS_n = LSIGMAS
  LSIG_CONV_n = LSIG_CONV
  LRMC01_n = LRMC01
  CTOM_n = CTOM
  CSUBG_AUCV_n = CSUBG_AUCV
  VSIGQSAT_n = VSIGQSAT
  CSUBG_AUCV_RI_n = CSUBG_AUCV_RI
  CCONDENS_n = CCONDENS
  CLAMBDA3_n = CLAMBDA3
  CSUBG_MF_PDF_n = CSUBG_MF_PDF
  LHGRAD_n = LHGRAD
  XCOEFHGRADTHL_n = XCOEFHGRADTHL
  XCOEFHGRADRM_n = XCOEFHGRADRM
  XALTHGRAD_n = XALTHGRAD
  XCLDTHOLD_n = XCLDTHOLD
END SUBROUTINE UPDATE_NAM_TURBn

END MODULE MODN_TURB_n
