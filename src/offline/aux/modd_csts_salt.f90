!ORILAM_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!ORILAM_LIC This is part of the ORILAM software governed by the CeCILL-C licence
!ORILAM_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!ORILAM_LIC for details.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 modd 2006/09/15 12:11:19
!-----------------------------------------------------------------
!      ######################
        MODULE MODD_CSTS_SALT
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
!!      Bielli S. 02/2019  Sea salt : significant sea wave height influences salt emission; 5 salt modes
!!
!!--------------------------------------------------------------------
!!     DECLARATIONS
!!     ------------
!
!
IMPLICIT NONE
!
!densité salt a introduire
! ++ PIERRE / MARINE SSA DUST - MODIF ++
REAL, PARAMETER  :: XDENSITY_DRYSALT  = 2.160e3     ![kg/m3] density of sea salt (dry NaCl 2.160E3) 
REAL, PARAMETER  :: XDENSITY_SALT     = 1.173e3     ![kg/m3] density of wet sea salt (Saltwater at RH80: 1.17e3) 
! -- PIERRE / MARINE SSA DUST - MODIF --
REAL, PARAMETER  :: XMOLARWEIGHT_SALT = 58.e-3   ![kg/mol] molar weight dust
REAL, PARAMETER  :: XM3TOUM3_SALT     = 1.d18     ![um3/m3] conversion factor
REAL, PARAMETER  :: XUM3TOM3_SALT     = 1.d-18    ![m3/um3] conversion factor
REAL, PARAMETER  :: XSIXTH_SALT       = 1./6.     ![-] one sixth
!
END MODULE MODD_CSTS_SALT
