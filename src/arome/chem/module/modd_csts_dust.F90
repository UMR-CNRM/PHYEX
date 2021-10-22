!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$ $Date$
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
USE MODD_CST, ONLY :    &
       XPI              & !Definition of pi
      ,XBOLTZ           & ! Boltzman constant 
      ,XAVOGADRO        & ![molec/mol] avogadros number
      ,XG               & ! Gravity constant
      ,XP00             & ! Reference pressure
      ,XMD              & ![kg/mol] molar weight of air
      ,XRD              & ! Gaz constant for dry air
      ,XCPD               !  Cpd (dry air)
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
