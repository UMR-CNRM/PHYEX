!ORILAM_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!ORILAM_LIC This is part of the ORILAM software governed by the CeCILL-C licence
!ORILAM_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!ORILAM_LIC for details.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 modd 2006/05/18 13:07:25
!-----------------------------------------------------------------
!      ######################
        MODULE MODD_CSTS_DUST
!      ######################
!!
!!     PURPOSE
!!     -------
!!
!!     Declaration of dust   constants                               
!!
!!     METHOD
!!     ------
!!
!!
!!     REFERENCE
!!     ---------
!!     none
!!
!!
!!     AUTHOR
!!     ------
!!     P.Tulet (GMEI)               
!!
!!
!!     MODIFICATIONS
!!     -------------
!!
!!--------------------------------------------------------------------
!!     DECLARATIONS
!!     ------------
!
!
IMPLICIT NONE
!
REAL, PARAMETER  :: XDENSITY_DUST = 2.5e3         ![kg/m3] density of dust
REAL, PARAMETER  :: XMOLARWEIGHT_DUST = 100.e-3   ![kg/mol] molar weight dust
REAL, PARAMETER  :: XM3TOUM3          = 1.d18     ![um3/m3] conversion factor
REAL, PARAMETER  :: XUM3TOM3          = 1.d-18    ![m3/um3] conversion factor
REAL, PARAMETER  :: XSIXTH            = 1./6.     ![-] one sixth
!
END MODULE MODD_CSTS_DUST
