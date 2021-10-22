!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source: /sxturb1/data1/mesonh/CODEMNH/modd_salt.f90,v $ $Revision: 1.1 $
! MASDEV4_7 modd 2007/01/12 14:42:16
!-----------------------------------------------------------------
!!     ######################
       MODULE MODD_SALT
!!     ######################
!!
!!     PURPOSE
!!     -------
!!
!!     declaration of variables and types for the sea salt scheme
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
!!     Pierre Tulet (CNRM)
!!
!!
!!     MODIFICATIONS
!!     -------------
!!
USE MODD_PARAMETERS, ONLY: JPMODELMAX
!!--------------------------------------------------------------------
!!     DECLARATIONS
!!     ------------
IMPLICIT NONE
!
LOGICAL      :: LSALT     = .FALSE.   ! switch to active pronostic sea salts
LOGICAL      :: LSLTINIT  = .FALSE.   ! switch to initialize pronostic sea salts
LOGICAL,DIMENSION(JPMODELMAX)  :: LDEPOS_SLT = .FALSE.    ! switch to SLT wet depositon
INTEGER      :: NMODE_SLT= 3  ! number of sea salt modes (max 3; default = 3)
!
CHARACTER(LEN=9),DIMENSION(:),ALLOCATABLE :: CDESLTNAMES
CHARACTER(LEN=9),DIMENSION(6), PARAMETER  :: YPDESLT_INI = &
     (/'DESLTM31C','DESLTM32C','DESLTM33C' &
      ,'DESLTM31R','DESLTM32R','DESLTM33R' /)
CHARACTER(LEN=6),DIMENSION(:),ALLOCATABLE :: CSALTNAMES
CHARACTER(LEN=6),DIMENSION(9), PARAMETER  :: YPSALT_INI = &
     (/'SLTM01','SLTM31','SLTM61' &
      ,'SLTM02','SLTM32','SLTM62' &
      ,'SLTM03','SLTM33','SLTM63' /)
! Set the order of the loops sorted by importance
!This means that if a user choses 1 mode it will have characteristics of mode 2
!2 modes will be mode 2 & 3, whereas 3 modes will modes 1, 2 and 3
INTEGER, DIMENSION(3),PARAMETER  :: JPSALTORDER = (/3, 2, 1/)
! 
REAL, ALLOCATABLE :: XSLTMSS(:,:,:)   ! [kg/m3] total mass concentration of sea salt
!
! aerosol lognormal parameterization
CHARACTER(LEN=4)  :: CRGUNITS   = 'MASS'  ! type of log-normal geometric mean radius
!                                         !given in namelist (mass on number)
!
LOGICAL      :: LRGFIX_SLT   = .FALSE.    ! switch to fix RG (sedimentation)
LOGICAL      :: LVARSIG_SLT  = .FALSE.    ! switch to active pronostic dispersion for all modes
LOGICAL      :: LSEDIMSALT   = .FALSE.    ! switch to active aerosol sedimentation
REAL         :: XSIGMIN_SLT   = 1.20      ! minimum dispersion value for sea salt mode
REAL         :: XSIGMAX_SLT   = 3.60      ! maximum dispersion value for sea salt mode
REAL         :: XCOEFRADMAX_SLT  = 10.    ! maximum increasement for Rg mode sea salt
REAL         :: XCOEFRADMIN_SLT  = 0.1    ! maximum decreasement for Rg mode sea salt
!

!Initial dry number median radius (um) from Vignati et al., 2001
!REAL, DIMENSION(3)          :: XINIRADIUS_SLT= (/0.2, 2., 12./)
!Initial, standard deviation from Vignati et al., 2001
!REAL, DIMENSION(3)          :: XINISIG_SLT =  (/1.9, 2., 3./)
!Initial dry number median radius (um) from Schultz et al., 2004
REAL, DIMENSION(3)          :: XINIRADIUS_SLT= 0.5*(/0.28, 2.25, 15.28/)
!Initial, standard deviation from Vignati et al., 2001
REAL, DIMENSION(3)          :: XINISIG_SLT =  (/1.9, 2., 2./)
!Minimum allowed number concentration for any mode (#/m3)
REAL, DIMENSION(3)          :: XN0MIN_SLT  = (/1.e4 , 1.e2 , 1.e-1 /)
!
END MODULE MODD_SALT
