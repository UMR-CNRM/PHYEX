!ORILAM_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!ORILAM_LIC This is part of the ORILAM software governed by the CeCILL-C licence
!ORILAM_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!ORILAM_LIC for details.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
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
!!     2014 P.Tulet modif XINIRADIUS_SLT and XN0MIN_SLT
!!
USE MODD_PARAMETERS, ONLY: JPMODELMAX
!USE MODD_SLT_n, ONLY : SLT_t
!!--------------------------------------------------------------------
!!     DECLARATIONS
!!     ------------
IMPLICIT NONE
!
! ++ PIERRE / MARINE SSA DUST - MODIF ++
LOGICAL      :: LSLTCAMS  = .FALSE.   ! switch to active pronostic sea salts  from CAMS
LOGICAL      :: LSALT     = .FALSE.   ! switch to active pronostic sea salts
LOGICAL      :: LONLY     = .FALSE.
LOGICAL      :: LREAD_ONLY_HS_MACC     = .FALSE.
LOGICAL      :: LSLTINIT  = .FALSE.   ! switch to initialize pronostic sea salts
LOGICAL      :: LSLTPRES  = .FALSE.   ! switch to know if pronostic salts exist            
LOGICAL,DIMENSION(JPMODELMAX)  :: LDEPOS_SLT = .FALSE.    ! switch to SLT wet depositon

CHARACTER(LEN=9),DIMENSION(:),ALLOCATABLE :: CDESLTNAMES
CHARACTER(LEN=6),DIMENSION(:), ALLOCATABLE :: CSALTNAMES
CHARACTER(LEN=9),DIMENSION(16), PARAMETER  :: YPDESLT_INI = &
    (/'DESLTM31C','DESLTM32C','DESLTM33C','DESLTM34C', 'DESLTM35C','DESLTM36C', &
       'DESLTM37C','DESLTM38C',& 
      'DESLTM31R','DESLTM32R','DESLTM33R', 'DESLTM34R','DESLTM35R','DESLTM36R','DESLTM37R','DESLTM38R' /)

CHARACTER(LEN=6),DIMENSION(24), PARAMETER  :: YPSALT_INI = &
                  (/'SLTM01','SLTM31','SLTM61',&
                    'SLTM02','SLTM32','SLTM62',&
                    'SLTM03','SLTM33','SLTM63',&
                    'SLTM04','SLTM34','SLTM64',&
                    'SLTM05','SLTM35','SLTM65',&
                    'SLTM06','SLTM36','SLTM66',&
                    'SLTM07','SLTM37','SLTM67',&
                    'SLTM08','SLTM38','SLTM68'/)


INTEGER, DIMENSION(8),PARAMETER  :: JPSALTORDER = (/1,2,3,4,5,6,7,8/)
INTEGER :: NMODE_SLT= 8  ! number of sea salt modes (default = 8)

!Test Thomas (definir rayons et sigma ici si on veut desactiver initialisation MACC)

!REAL, DIMENSION(5) :: XINIRADIUS_SLT,XINISIG_SLT,XN0MIN_SLT

!Initial dry number median radius (um) from Ova et al., 2014
REAL,DIMENSION(8),TARGET    :: XINIRADIUS_SLT=  (/0.009, 0.021, 0.045, 0.115,0.415, 2.5, 7.0, 25.0/)
!Initial, standard deviation from  Ova et al., 2014
REAL,DIMENSION(8),TARGET      :: XINISIG_SLT =  (/ 1.37, 1.5, 1.42, 1.53, 1.85, 1.55, 1.8, 2.1 /)

!Minimum allowed number concentration for any mode (#/m3)
REAL,DIMENSION(8)  :: XN0MIN_SLT  = (/100. , 100., 100., 10., 1., 1E-2 , 1E-3, 1E-4 /)

!Test Thomas
REAL, DIMENSION(:,:,:), ALLOCATABLE :: XSLTMSS   ! [kg/m3] total mass concentration of sea salt
!
! aerosol lognormal parameterization
CHARACTER(LEN=4)  :: CRGUNITS   = 'NUMB'  ! type of log-normal geometric mean radius
!                                         !given in namelist (mass on number)
!
LOGICAL      :: LRGFIX_SLT   = .FALSE.    ! switch to fix RG (sedimentation)
LOGICAL      :: LVARSIG_SLT  = .FALSE.    ! switch to active pronostic dispersion for all modes
LOGICAL      :: LSEDIMSALT   = .FALSE.    ! switch to active aerosol sedimentation
REAL         :: XSIGMIN_SLT   = 1.20      ! minimum dispersion value for sea salt mode
!REAL         :: XSIGMIN_SLT   = 0.      ! minimum dispersion value for sea salt mode
REAL         :: XSIGMAX_SLT   = 3.60      ! maximum dispersion value for sea salt mode
REAL         :: XCOEFRADMAX_SLT  = 10.    ! maximum increasement for Rg mode sea salt
REAL         :: XCOEFRADMIN_SLT  = 0.1    ! minimum decreasement for Rg mode sea salt
!REAL         :: XCOEFRADMIN_SLT  = 0.    ! minimum decreasement for Rg mode sea salt


!
! -- PIERRE / MARINE SSA DUST - MODIF --
!
END MODULE MODD_SALT
