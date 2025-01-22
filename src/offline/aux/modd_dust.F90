!ORILAM_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!ORILAM_LIC This is part of the ORILAM software governed by the CeCILL-C licence
!ORILAM_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!ORILAM_LIC for details.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 modd 2007/02/22 10:05:39
!-----------------------------------------------------------------
!!     ######################
       MODULE MODD_DUST
!!     ######################
!!
!!     PURPOSE
!!     -------
!!
!!     declaration of variables and types for the dust scheme
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
!!     T. Hoarau  03/2019   add a switch for initialisation from MACC
!!--------------------------------------------------------------------
!!     DECLARATIONS
!!     ------------
USE MODD_PARAMETERS, ONLY: JPMODELMAX
!
IMPLICIT NONE
!
LOGICAL      :: LDUST     = .FALSE.   ! switch to active pronostic dusts
LOGICAL      :: LDSTCAMS  = .FALSE.   ! switch to active pronostic dusts from MACC
LOGICAL      :: LDSTINIT  = .FALSE.   ! switch to initialize pronostic dusts
LOGICAL      :: LDSTPRES  = .FALSE.   ! switch to know if pronostic dusts exist
LOGICAL,DIMENSION(JPMODELMAX)  :: LDEPOS_DST = .FALSE.    ! switch to DST wet deposition
INTEGER      :: NMODE_DST= 3  ! number of dust modes (max 3; default = 3)
!
CHARACTER(LEN=6),DIMENSION(:),ALLOCATABLE :: CDUSTNAMES

CHARACTER(LEN=6),DIMENSION(9), PARAMETER  :: YPDUST_INI = &
     (/'DSTM01','DSTM31','DSTM61' &
      ,'DSTM02','DSTM32','DSTM62' &
      ,'DSTM03','DSTM33','DSTM63' /)
! Set the order of the loops sorted by importance
!This means that if a user choses 1 mode it will have characteristics of mode 2
!2 modes will be mode 2 & 3, whereas 3 modes will modes 1, 2 and 3
INTEGER, DIMENSION(3),PARAMETER  :: JPDUSTORDER = (/3, 2, 1/)
REAL         :: XRADMIN = 0.001 ! minimum reasonable value for median radius
! 
REAL, ALLOCATABLE :: XDSTMSS(:,:,:)   ! [kg/m3] total mass concentration of dust
!
! aerosol lognormal parameterization
CHARACTER(LEN=4)  :: CRGUNITD   = 'NUMB'  ! type of log-normal geometric mean radius
!                                          !given in namelist (mass on number)
!
LOGICAL      :: LRGFIX_DST = .FALSE.  ! switch to fix RG (sedimentation)
LOGICAL      :: LVARSIG    = .FALSE.  ! switch to active pronostic dispersion for all modes
LOGICAL      :: LSEDIMDUST = .FALSE.  ! switch to active aerosol sedimentation
REAL         :: XSIGMIN   = 1.20      ! minimum dispersion value for dust mode
REAL         :: XSIGMAX   = 3.60      ! maximum dispersion value for dust mode
REAL         :: XCOEFRADMAX  = 10.    ! maximum increasement for Rg mode dust
REAL         :: XCOEFRADMIN  = 0.1    ! maximum decreasement for Rg mode dust
!
! Alf considers it is better to use initial values as  Schultz et al 1998
! whereas Pierre consider to keep as close as possible initialization
! values close to default emissions; so as you want !!!
!Initial dry mass median radius (um) from Schultz et al 1998
!REAL, DIMENSION(3)          :: XINIRADIUS= (/ 0.0055, 1.26, 21.65 /)
!Initial, standard deviation from Schultz et al 1998
!REAL, DIMENSION(3)          :: XINISIG =  (/2.13, 2.00, 1.89 /)
!Initial dry mass median radius (um) from D'Almeida, 1987 emission fluxes
!REAL, DIMENSION(3)          :: XINIRADIUS= 0.5*(/ 0.832 ,  4.82 , 19.38 /)
!Initial, standard deviation from from D'Almeida, 1987 emission fluxes
!REAL, DIMENSION(3)          :: XINISIG =  (/2.10     ,  1.90    ,  1.60 /)
!Minimum allowed number concentration for any mode (#/m3)
!REAL, DIMENSION(3)          :: XN0MIN  = (/1.e4 , 1.e3 , 1.e-1 /)
!Initial dry mass median radius (um) from Alfaro et al 1998
!REAL, DIMENSION(3)          :: XINIRADIUS= 0.5*(/1.5, 6.7, 14.2/)
!Initial, standard deviation from Alfaro et al 1998
!REAL, DIMENSION(3)          :: XINISIG =  (/1.70, 1.60, 1.50/)
!Minimum allowed number concentration for any mode (#/m3)
!REAL, DIMENSION(3)          :: XN0MIN  = (/1.e2 , 1.e1 , 1.e-2 /)
!
! NEW PARAMETERIZATION FROM AMMA, defalut
!Initial dry number median radius (um) 
REAL, DIMENSION(3)          :: XINIRADIUS= 0.5*(/0.078, 0.641, 5.00 /)
!Initial, standard deviation from Alfaro et al 1998
REAL, DIMENSION(3)          :: XINISIG =  (/1.75, 1.76, 1.70/)
!Minimum allowed number concentration for any mode (#/m3)
REAL, DIMENSION(3)          :: XN0MIN  = (/1.e1 , 1.e-1 , 1.e-4 /)
CHARACTER(LEN=9),DIMENSION(:),ALLOCATABLE :: CDEDSTNAMES
CHARACTER(LEN=9),DIMENSION(6), PARAMETER  :: YPDEDST_INI = &
     (/'DEDSTM31C','DEDSTM32C','DEDSTM33C' &
      ,'DEDSTM31R','DEDSTM32R','DEDSTM33R' /)
!
END MODULE MODD_DUST
