!MNH_LIC Copyright 2001-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-------------------------------------------------------------------------------
!     ######################
      MODULE MODN_PARAM_LIMA
!     ######################
!
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_PARAM_LIMA
!
IMPLICIT NONE
!
!
NAMELIST/NAM_PARAM_LIMA/LNUCL, LSEDI, LHHONI, LMEYERS,                     &
                        NMOM_I, NMOM_S, NMOM_G, NMOM_H,                    &
                        NMOD_IFN, XIFN_CONC, LIFN_HOM,                     &
                        CIFN_SPECIES, CINT_MIXING, NMOD_IMM, NIND_SPECIE,  &
                        LSNOW_T, CPRISTINE_ICE_LIMA, CHEVRIMED_ICE_LIMA,   &
                        XALPHAI, XNUI, XALPHAS, XNUS, XALPHAG, XNUG,       &
                        XFACTNUC_DEP, XFACTNUC_CON, NPHILLIPS,             &
                        LCIBU, XNDEBRIS_CIBU, LRDSF, LMURAKAMI,            &
                        LACTI, LSEDC, LACTIT, LBOUND, LSPRO,               &
                        LADJ, LKHKO, LKESSLERAC, NMOM_C, NMOM_R,           &
                        NMOD_CCN, XCCN_CONC,                               &
                        LCCN_HOM, CCCN_MODES, HINI_CCN, HTYPE_CCN,         &
                        XALPHAC, XNUC, XALPHAR, XNUR,                      &
                        XFSOLUB_CCN, XACTEMP_CCN, XAERDIFF, XAERHEIGHT,    &
                        LSCAV, LAERO_MASS, LDEPOC, XVDEPOC, LACTTKE,       &
                        LPTSPLIT, LFEEDBACKT, NMAXITER, XMRSTEP, XTSTEP_TS
!
END MODULE MODN_PARAM_LIMA
