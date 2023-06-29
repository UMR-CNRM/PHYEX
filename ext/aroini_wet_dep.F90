!     ######spl
SUBROUTINE AROINI_WET_DEP
!**** *INI_MICRO*   - Initialize common meso_NH MODD_ used in microphysics for AROME

!     Purpose.
!     --------
!           Initialize 
!           MODD_WET_DEP_DESCR, MODD_WET-DEP_PARAM  
!           parameters used in ALADIN_DUST 

!**   Interface.
!     ----------
!        *CALL* *INI_MICRO (KULOUT,KSTEP,KSPLITR)

!        Explicit arguments :
!        --------------------
!        KULOUT : Logical unit for the output
!        PTSTEP  : Time step
!        KSPLITR : Number of small time step interation for rain sedimentation 
!        LDWARM : value assigned to LWARM       

!        Implicit arguments :
!        --------------------
!        

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        Documentation AROME 

!     Author.
!     -------
!        Y. Seity 

!     Modifications.
!     --------------
!        Original : 03-12-12
!        05-08-25 Kovacic  added LDWARM
!        01-02-2011: M. Mokhtari Adaptation of aroini_micro for Aladin
!     ------------------------------------------------------------------

USE MODD_WET_DEP_DESCR
USE MODD_WET_DEP_PARAM

USE MODI_INI_WET_DEP
USE MODE_INI_CST, ONLY: INI_CST
IMPLICIT NONE
! -----------------------------------------------------------------------
!     DUMMY INTEGER SCALARS
! -----------------------------------------------------------------------
!        1.1 Set implicit default values for MODD_PARAMETERS
!       les variables sont initialiso?=es dans le module lui mo?=me 
!        1.2 Set implicit default values for MODD_CST
CALL INI_CST

!        1. Set implicit default values for MODD_PARAM_ICE

!        2. Set implicit default values for MODD_RAIN_ICE_DESCR 
!                     et MODD_RAIN_ICE_PARAM

CALL INI_WET_DEP 
! -----------------------------------------------------------------------

RETURN
END SUBROUTINE AROINI_WET_DEP
