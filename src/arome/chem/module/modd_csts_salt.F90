!     ######spl
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
!densité salt a introduire
REAL, PARAMETER  :: XDENSITY_SALT     = 2.1e3     ![kg/m3] density of dust
REAL, PARAMETER  :: XMOLARWEIGHT_SALT = 58.e-3   ![kg/mol] molar weight dust
REAL, PARAMETER  :: XM3TOUM3_SALT     = 1.d18     ![um3/m3] conversion factor
REAL, PARAMETER  :: XUM3TOM3_SALT     = 1.d-18    ![m3/um3] conversion factor
REAL, PARAMETER  :: XSIXTH_SALT       = 1./6.     ![-] one sixth
!
END MODULE MODD_CSTS_SALT
