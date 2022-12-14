!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$ $Date$
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!     ########################## 
      SUBROUTINE SWITCH_SBG_LES_n
!     ##########################
!
!!****  *SWITCH_SBG_LESn* - moves LES subgrid quantities from modd_les
!!                          to modd_lesn or the contrary.
!!
!!    PURPOSE
!!    -------
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------ 
!!
!!    REFERENCE
!!    ---------
!!      
!!
!!    AUTHOR
!!    ------
!!	   V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    June 14, 2002                
!!
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_LES
USE MODD_LES_n
USE MODD_CONF_n
USE MODD_NSV
!
USE MODI_SECOND_MNH
!
IMPLICIT NONE
!
REAL :: ZTIME1, ZTIME2
!-------------------------------------------------------------------------------
!
!*      7.4  interactions of resolved and subgrid quantities
!            -----------------------------------------------
!
CALL SECOND_MNH(ZTIME1)
!
IF (.NOT.   ASSOCIATED (X_LES_RES_W_SBG_WThl) ) THEN
!                                                                           ______
  CALL LES_ALLOCATE('X_LES_RES_W_SBG_WThl',(/NLES_K,NLES_TIMES,NLES_MASKS/))!  <w'w'Thl'>
!                                                                        _____
  CALL LES_ALLOCATE('X_LES_RES_W_SBG_Thl2',(/NLES_K,NLES_TIMES,NLES_MASKS/))!  <w'Thl'2>
!                                                                              _____
  CALL LES_ALLOCATE('X_LES_RES_ddxa_U_SBG_UaU',(/NLES_K,NLES_TIMES,NLES_MASKS/))!  <du'/dxa ua'u'>
!                                                                              _____
  CALL LES_ALLOCATE('X_LES_RES_ddxa_V_SBG_UaV',(/NLES_K,NLES_TIMES,NLES_MASKS/))!  <dv'/dxa ua'v'>
!                                                                             _____
  CALL LES_ALLOCATE('X_LES_RES_ddxa_W_SBG_UaW',(/NLES_K,NLES_TIMES,NLES_MASKS/))!  <dw'/dxa ua'w'>
!                                                                              _______
  CALL LES_ALLOCATE('X_LES_RES_ddxa_W_SBG_UaThl',(/NLES_K,NLES_TIMES,NLES_MASKS/))!  <dw'/dxa ua'Thl'>
!                                                                                _____
  CALL LES_ALLOCATE('X_LES_RES_ddxa_Thl_SBG_UaW',(/NLES_K,NLES_TIMES,NLES_MASKS/))!  <dThl'/dxa ua'w'>
!                                                                                  ___
  CALL LES_ALLOCATE('X_LES_RES_ddz_Thl_SBG_W2',(/NLES_K,NLES_TIMES,NLES_MASKS/))!  <dThl'/dz w'2>
!                                                                                _______
  CALL LES_ALLOCATE('X_LES_RES_ddxa_Thl_SBG_UaThl',(/NLES_K,NLES_TIMES,NLES_MASKS/))!  <dThl'/dxa ua'Thl'>
!
  IF (LUSERV) THEN
!                                                                          _____
    CALL LES_ALLOCATE('X_LES_RES_W_SBG_WRt',(/NLES_K,NLES_TIMES,NLES_MASKS/))!  <w'w'Rt'>
!                                                                           ____
    CALL LES_ALLOCATE('X_LES_RES_W_SBG_Rt2',(/NLES_K,NLES_TIMES,NLES_MASKS/))!  <w'Rt'2>
!                                                                          _______
    CALL LES_ALLOCATE('X_LES_RES_W_SBG_ThlRt',(/NLES_K,NLES_TIMES,NLES_MASKS/))!  <w'Thl'Rt'>
!                                                                                ______
    CALL LES_ALLOCATE('X_LES_RES_ddxa_W_SBG_UaRt',(/NLES_K,NLES_TIMES,NLES_MASKS/))!  <dw'/dxa ua'Rt'>
!                                                                                 _____
    CALL LES_ALLOCATE('X_LES_RES_ddxa_Rt_SBG_UaW',(/NLES_K,NLES_TIMES,NLES_MASKS/))!  <dRt'/dxa ua'w'>
!                                                                                   ___
    CALL LES_ALLOCATE('X_LES_RES_ddz_Rt_SBG_W2',(/NLES_K,NLES_TIMES,NLES_MASKS/))!  <dRt'/dz w'2>
!                                                                                  ______
    CALL LES_ALLOCATE('X_LES_RES_ddxa_Thl_SBG_UaRt',(/NLES_K,NLES_TIMES,NLES_MASKS/))!  <dThl'/dxa ua'Rt'>
!                                                                                 _______
    CALL LES_ALLOCATE('X_LES_RES_ddxa_Rt_SBG_UaThl',(/NLES_K,NLES_TIMES,NLES_MASKS/))!  <dRt'/dxa ua'Thl'>
!                                                                                  ______
    CALL LES_ALLOCATE('X_LES_RES_ddxa_Rt_SBG_UaRt',(/NLES_K,NLES_TIMES,NLES_MASKS/)) !  <dRt'/dxa ua'Rt'>
  ELSE
    CALL LES_ALLOCATE('X_LES_RES_W_SBG_WRt',(/0,0,0/))
    CALL LES_ALLOCATE('X_LES_RES_W_SBG_Rt2',(/0,0,0/))
    CALL LES_ALLOCATE('X_LES_RES_W_SBG_ThlRt',(/0,0,0/))
    CALL LES_ALLOCATE('X_LES_RES_ddxa_W_SBG_UaRt',(/0,0,0/))
    CALL LES_ALLOCATE('X_LES_RES_ddxa_Rt_SBG_UaW',(/0,0,0/))
    CALL LES_ALLOCATE('X_LES_RES_ddz_Rt_SBG_W2',(/0,0,0/))
    CALL LES_ALLOCATE('X_LES_RES_ddxa_Thl_SBG_UaRt',(/0,0,0/))
    CALL LES_ALLOCATE('X_LES_RES_ddxa_Rt_SBG_UaThl',(/0,0,0/))
    CALL LES_ALLOCATE('X_LES_RES_ddxa_Rt_SBG_UaRt',(/0,0,0/))
  END IF
!                                                                                   ______
CALL LES_ALLOCATE('X_LES_RES_ddxa_W_SBG_UaSv',(/NLES_K,NLES_TIMES,NLES_MASKS,NSV/))  !  <dw'/dxa ua'Sv'>
!                                                                                    _____
CALL LES_ALLOCATE('X_LES_RES_ddxa_Sv_SBG_UaW',(/NLES_K,NLES_TIMES,NLES_MASKS,NSV/))  !  <dSv'/dxa ua'w'>
!                                                                                    ___
CALL LES_ALLOCATE('X_LES_RES_ddz_Sv_SBG_W2',(/NLES_K,NLES_TIMES,NLES_MASKS,NSV/) ) !  <dSv'/dz w'2>
!                                                                                   ______
CALL LES_ALLOCATE('X_LES_RES_ddxa_Sv_SBG_UaSv',(/NLES_K,NLES_TIMES,NLES_MASKS,NSV/))  !  <dSv'/dxa ua'Sv'>
!                                                                             _____
CALL LES_ALLOCATE('X_LES_RES_W_SBG_WSv',(/NLES_K,NLES_TIMES,NLES_MASKS,NSV/))  !  <w'w'Sv'>
!                                                                             ____
CALL LES_ALLOCATE('X_LES_RES_W_SBG_Sv2',(/NLES_K,NLES_TIMES,NLES_MASKS,NSV/))  !  <w'Sv'2>
!
!
  X_LES_RES_W_SBG_WThl         = XLES_RES_W_SBG_WThl
  X_LES_RES_W_SBG_Thl2         = XLES_RES_W_SBG_Thl2
  X_LES_RES_ddxa_U_SBG_UaU     = XLES_RES_ddxa_U_SBG_UaU
  X_LES_RES_ddxa_V_SBG_UaV     = XLES_RES_ddxa_V_SBG_UaV
  X_LES_RES_ddxa_W_SBG_UaW     = XLES_RES_ddxa_W_SBG_UaW
  X_LES_RES_ddxa_W_SBG_UaThl   = XLES_RES_ddxa_W_SBG_UaThl
  X_LES_RES_ddxa_Thl_SBG_UaW   = XLES_RES_ddxa_Thl_SBG_UaW
  X_LES_RES_ddz_Thl_SBG_W2     = XLES_RES_ddz_Thl_SBG_W2
  X_LES_RES_ddxa_Thl_SBG_UaThl = XLES_RES_ddxa_Thl_SBG_UaThl
  IF (LUSERV) THEN
    X_LES_RES_W_SBG_WRt        = XLES_RES_W_SBG_WRt
    X_LES_RES_W_SBG_Rt2        = XLES_RES_W_SBG_Rt2
    X_LES_RES_W_SBG_ThlRt      = XLES_RES_W_SBG_ThlRt
    X_LES_RES_ddxa_W_SBG_UaRt  = XLES_RES_ddxa_W_SBG_UaRt
    X_LES_RES_ddxa_Rt_SBG_UaW  = XLES_RES_ddxa_Rt_SBG_UaW
    X_LES_RES_ddz_Rt_SBG_W2    = XLES_RES_ddz_Rt_SBG_W2
    X_LES_RES_ddxa_Thl_SBG_UaRt= XLES_RES_ddxa_Thl_SBG_UaRt
    X_LES_RES_ddxa_Rt_SBG_UaThl= XLES_RES_ddxa_Rt_SBG_UaThl
    X_LES_RES_ddxa_Rt_SBG_UaRt = XLES_RES_ddxa_Rt_SBG_UaRt
  END IF
  IF (NSV>0) THEN
    X_LES_RES_ddxa_W_SBG_UaSv  = XLES_RES_ddxa_W_SBG_UaSv
    X_LES_RES_ddxa_Sv_SBG_UaW  = XLES_RES_ddxa_Sv_SBG_UaW
    X_LES_RES_ddz_Sv_SBG_W2    = XLES_RES_ddz_Sv_SBG_W2
    X_LES_RES_ddxa_Sv_SBG_UaSv = XLES_RES_ddxa_Sv_SBG_UaSv
    X_LES_RES_W_SBG_WSv        = XLES_RES_W_SBG_WSv
    X_LES_RES_W_SBG_Sv2        = XLES_RES_W_SBG_Sv2
  END IF
!
!
  CALL LES_ALLOCATE('X_LES_SUBGRID_U2',(/NLES_K,NLES_TIMES,NLES_MASKS/))     ! <u'2>
  CALL LES_ALLOCATE('X_LES_SUBGRID_V2',(/NLES_K,NLES_TIMES,NLES_MASKS/))     ! <v'2>
  CALL LES_ALLOCATE('X_LES_SUBGRID_W2',(/NLES_K,NLES_TIMES,NLES_MASKS/))     ! <w'2>
  CALL LES_ALLOCATE('X_LES_SUBGRID_Thl2',(/NLES_K,NLES_TIMES,NLES_MASKS/))     ! <Thl'2>
  CALL LES_ALLOCATE('X_LES_SUBGRID_UV',(/NLES_K,NLES_TIMES,NLES_MASKS/))     ! <u'v'>
  CALL LES_ALLOCATE('X_LES_SUBGRID_WU',(/NLES_K,NLES_TIMES,NLES_MASKS/))     ! <w'u'>
  CALL LES_ALLOCATE('X_LES_SUBGRID_WV',(/NLES_K,NLES_TIMES,NLES_MASKS/))     ! <w'v'>
  CALL LES_ALLOCATE('X_LES_SUBGRID_UThl',(/NLES_K,NLES_TIMES,NLES_MASKS/))     ! <u'Thl'>
  CALL LES_ALLOCATE('X_LES_SUBGRID_VThl',(/NLES_K,NLES_TIMES,NLES_MASKS/))     ! <v'Thl'>
  CALL LES_ALLOCATE('X_LES_SUBGRID_WThl',(/NLES_K,NLES_TIMES,NLES_MASKS/))     ! <w'Thl'>
  CALL LES_ALLOCATE('X_LES_SUBGRID_WThv',(/NLES_K,NLES_TIMES,NLES_MASKS/))  ! <w'Thv'>
  CALL LES_ALLOCATE('X_LES_SUBGRID_ThlThv',(/NLES_K,NLES_TIMES,NLES_MASKS/))  ! <Thl'Thv'>
  CALL LES_ALLOCATE('X_LES_SUBGRID_W2Thl',(/NLES_K,NLES_TIMES,NLES_MASKS/))  ! <w'2Thl>
  CALL LES_ALLOCATE('X_LES_SUBGRID_WThl2',(/NLES_K,NLES_TIMES,NLES_MASKS/))  ! <w'Thl'2>
  CALL LES_ALLOCATE('X_LES_SUBGRID_DISS_Tke',(/NLES_K,NLES_TIMES,NLES_MASKS/))  ! <epsilon>
  CALL LES_ALLOCATE('X_LES_SUBGRID_DISS_Thl2',(/NLES_K,NLES_TIMES,NLES_MASKS/))  ! <epsilon_Thl2>
  CALL LES_ALLOCATE('X_LES_SUBGRID_WP',(/NLES_K,NLES_TIMES,NLES_MASKS/))  ! <w'p'>
  CALL LES_ALLOCATE('X_LES_SUBGRID_PHI3',(/NLES_K,NLES_TIMES,NLES_MASKS/))  !  phi3
  CALL LES_ALLOCATE('X_LES_SUBGRID_LMix',(/NLES_K,NLES_TIMES,NLES_MASKS/))  !  Lmix
  CALL LES_ALLOCATE('X_LES_SUBGRID_LDiss',(/NLES_K,NLES_TIMES,NLES_MASKS/))  !  Ldiss
  CALL LES_ALLOCATE('X_LES_SUBGRID_Km',(/NLES_K,NLES_TIMES,NLES_MASKS/))  !  Km
  CALL LES_ALLOCATE('X_LES_SUBGRID_Kh',(/NLES_K,NLES_TIMES,NLES_MASKS/))  !  Kh
  CALL LES_ALLOCATE('X_LES_SUBGRID_ThlPz',(/NLES_K,NLES_TIMES,NLES_MASKS/))  ! <Thl'dp'/dz>
  CALL LES_ALLOCATE('X_LES_SUBGRID_UTke',(/NLES_K,NLES_TIMES,NLES_MASKS/))     ! <u'Tke>
  CALL LES_ALLOCATE('X_LES_SUBGRID_VTke',(/NLES_K,NLES_TIMES,NLES_MASKS/))     ! <v'Tke>
  CALL LES_ALLOCATE('X_LES_SUBGRID_WTke',(/NLES_K,NLES_TIMES,NLES_MASKS/))     ! <w'Tke>
  CALL LES_ALLOCATE('X_LES_SUBGRID_ddz_WTke',(/NLES_K,NLES_TIMES,NLES_MASKS/)) ! <dw'Tke/dz>

  CALL LES_ALLOCATE('X_LES_SUBGRID_THLUP_MF',(/NLES_K,NLES_TIMES,NLES_MASKS/))  ! Thl of the Updraft
  CALL LES_ALLOCATE('X_LES_SUBGRID_RTUP_MF',(/NLES_K,NLES_TIMES,NLES_MASKS/))  ! Rt of the Updraft
  CALL LES_ALLOCATE('X_LES_SUBGRID_RVUP_MF',(/NLES_K,NLES_TIMES,NLES_MASKS/))  ! Rv of the Updraft
  CALL LES_ALLOCATE('X_LES_SUBGRID_RCUP_MF',(/NLES_K,NLES_TIMES,NLES_MASKS/))  ! Rc of the Updraft
  CALL LES_ALLOCATE('X_LES_SUBGRID_RIUP_MF',(/NLES_K,NLES_TIMES,NLES_MASKS/))  ! Ri of the Updraft
  CALL LES_ALLOCATE('X_LES_SUBGRID_WUP_MF',(/NLES_K,NLES_TIMES,NLES_MASKS/))  ! Thl of the Updraft
  CALL LES_ALLOCATE('X_LES_SUBGRID_MASSFLUX',(/NLES_K,NLES_TIMES,NLES_MASKS/))  ! Mass Flux
  CALL LES_ALLOCATE('X_LES_SUBGRID_DETR',(/NLES_K,NLES_TIMES,NLES_MASKS/))  ! Detrainment
  CALL LES_ALLOCATE('X_LES_SUBGRID_ENTR',(/NLES_K,NLES_TIMES,NLES_MASKS/))  ! Entrainment
  CALL LES_ALLOCATE('X_LES_SUBGRID_FRACUP',(/NLES_K,NLES_TIMES,NLES_MASKS/))  ! Updraft Fraction 
  CALL LES_ALLOCATE('X_LES_SUBGRID_THVUP_MF',(/NLES_K,NLES_TIMES,NLES_MASKS/))  ! Thv of the Updraft
  CALL LES_ALLOCATE('X_LES_SUBGRID_WTHLMF',(/NLES_K,NLES_TIMES,NLES_MASKS/))! Flux of thl   
  CALL LES_ALLOCATE('X_LES_SUBGRID_WRTMF',(/NLES_K,NLES_TIMES,NLES_MASKS/)) ! Flux of rt
  CALL LES_ALLOCATE('X_LES_SUBGRID_WTHVMF',(/NLES_K,NLES_TIMES,NLES_MASKS/)) ! Flux of thv 
  CALL LES_ALLOCATE('X_LES_SUBGRID_WUMF',(/NLES_K,NLES_TIMES,NLES_MASKS/))! Flux of u
  CALL LES_ALLOCATE('X_LES_SUBGRID_WVMF',(/NLES_K,NLES_TIMES,NLES_MASKS/))! Flux of v
  
  IF (LUSERV ) THEN
    CALL LES_ALLOCATE('X_LES_SUBGRID_Rt2',(/NLES_K,NLES_TIMES,NLES_MASKS/))     ! <Rt'2>
    CALL LES_ALLOCATE('X_LES_SUBGRID_ThlRt',(/NLES_K,NLES_TIMES,NLES_MASKS/))     ! <Thl'Rt'>
    CALL LES_ALLOCATE('X_LES_SUBGRID_URt',(/NLES_K,NLES_TIMES,NLES_MASKS/))     ! <u'Rt'>
    CALL LES_ALLOCATE('X_LES_SUBGRID_VRt',(/NLES_K,NLES_TIMES,NLES_MASKS/))     ! <v'Rt'>
    CALL LES_ALLOCATE('X_LES_SUBGRID_WRt',(/NLES_K,NLES_TIMES,NLES_MASKS/))     ! <w'Rt'>
    CALL LES_ALLOCATE('X_LES_SUBGRID_RtThv',(/NLES_K,NLES_TIMES,NLES_MASKS/))     ! <Rt'Thv'>
    CALL LES_ALLOCATE('X_LES_SUBGRID_W2Rt',(/NLES_K,NLES_TIMES,NLES_MASKS/))     ! <w'2Rt'>
    CALL LES_ALLOCATE('X_LES_SUBGRID_WThlRt',(/NLES_K,NLES_TIMES,NLES_MASKS/))     ! <w'Thl'Rt'>
    CALL LES_ALLOCATE('X_LES_SUBGRID_WRt2',(/NLES_K,NLES_TIMES,NLES_MASKS/))     ! <w'Rt'2>
    CALL LES_ALLOCATE('X_LES_SUBGRID_DISS_Rt2',(/NLES_K,NLES_TIMES,NLES_MASKS/))  ! <epsilon_Rt2>
    CALL LES_ALLOCATE('X_LES_SUBGRID_DISS_ThlRt',(/NLES_K,NLES_TIMES,NLES_MASKS/)) ! <epsilon_ThlRt>
    CALL LES_ALLOCATE('X_LES_SUBGRID_RtPz',(/NLES_K,NLES_TIMES,NLES_MASKS/))     ! <Rt'dp'/dz>
    CALL LES_ALLOCATE('X_LES_SUBGRID_PSI3',(/NLES_K,NLES_TIMES,NLES_MASKS/))     !  psi3  
  ELSE
    CALL LES_ALLOCATE('X_LES_SUBGRID_Rt2',(/0,0,0/))
    CALL LES_ALLOCATE('X_LES_SUBGRID_ThlRt',(/0,0,0/))
    CALL LES_ALLOCATE('X_LES_SUBGRID_URt',(/0,0,0/))
    CALL LES_ALLOCATE('X_LES_SUBGRID_VRt',(/0,0,0/))
    CALL LES_ALLOCATE('X_LES_SUBGRID_WRt',(/0,0,0/))
    CALL LES_ALLOCATE('X_LES_SUBGRID_RtThv',(/0,0,0/))
    CALL LES_ALLOCATE('X_LES_SUBGRID_W2Rt',(/0,0,0/))
    CALL LES_ALLOCATE('X_LES_SUBGRID_WThlRt',(/0,0,0/))
    CALL LES_ALLOCATE('X_LES_SUBGRID_WRt2',(/0,0,0/))
    CALL LES_ALLOCATE('X_LES_SUBGRID_DISS_Rt2',(/0,0,0/))
    CALL LES_ALLOCATE('X_LES_SUBGRID_DISS_ThlRt',(/0,0,0/))
    CALL LES_ALLOCATE('X_LES_SUBGRID_RtPz',(/0,0,0/))
    CALL LES_ALLOCATE('X_LES_SUBGRID_PSI3',(/0,0,0/))    
  END IF
  IF (LUSERC ) THEN
    CALL LES_ALLOCATE('X_LES_SUBGRID_Rc2',(/NLES_K,NLES_TIMES,NLES_MASKS/))     ! <Rc'2>
    CALL LES_ALLOCATE('X_LES_SUBGRID_URc',(/NLES_K,NLES_TIMES,NLES_MASKS/))     ! <u'Rc'>
    CALL LES_ALLOCATE('X_LES_SUBGRID_VRc',(/NLES_K,NLES_TIMES,NLES_MASKS/))     ! <v'Rc'>
    CALL LES_ALLOCATE('X_LES_SUBGRID_WRc',(/NLES_K,NLES_TIMES,NLES_MASKS/))     ! <w'Rc'>
  ELSE
    CALL LES_ALLOCATE('X_LES_SUBGRID_Rc2',(/0,0,0/))
    CALL LES_ALLOCATE('X_LES_SUBGRID_URc',(/0,0,0/))
    CALL LES_ALLOCATE('X_LES_SUBGRID_VRc',(/0,0,0/))
    CALL LES_ALLOCATE('X_LES_SUBGRID_WRc',(/0,0,0/))
  END IF
  IF (LUSERI ) THEN
    CALL LES_ALLOCATE('X_LES_SUBGRID_Ri2',(/NLES_K,NLES_TIMES,NLES_MASKS/))     ! <Ri'2>
  ELSE
    CALL LES_ALLOCATE('X_LES_SUBGRID_Ri2',(/0,0,0/))
  END IF
  IF (NSV>0  ) THEN
    CALL LES_ALLOCATE('X_LES_SUBGRID_USv',(/NLES_K,NLES_TIMES,NLES_MASKS,NSV/)) ! <u'Sv'>
    CALL LES_ALLOCATE('X_LES_SUBGRID_VSv',(/NLES_K,NLES_TIMES,NLES_MASKS,NSV/)) ! <v'Sv'>
    CALL LES_ALLOCATE('X_LES_SUBGRID_WSv',(/NLES_K,NLES_TIMES,NLES_MASKS,NSV/)) ! <w'Sv'>
    CALL LES_ALLOCATE('X_LES_SUBGRID_Sv2',(/NLES_K,NLES_TIMES,NLES_MASKS,NSV/))  ! <Sv'2>
    CALL LES_ALLOCATE('X_LES_SUBGRID_SvThv',(/NLES_K,NLES_TIMES,NLES_MASKS,NSV/))  ! <Sv'Thv'>
    CALL LES_ALLOCATE('X_LES_SUBGRID_W2Sv',(/NLES_K,NLES_TIMES,NLES_MASKS,NSV/))  ! <w'2Sv'>
    CALL LES_ALLOCATE('X_LES_SUBGRID_WSv2',(/NLES_K,NLES_TIMES,NLES_MASKS,NSV/))  ! <w'Sv'2>
    CALL LES_ALLOCATE('X_LES_SUBGRID_DISS_Sv2',(/NLES_K,NLES_TIMES,NLES_MASKS,NSV/))  ! <epsilon_Sv2>
    CALL LES_ALLOCATE('X_LES_SUBGRID_SvPz',(/NLES_K,NLES_TIMES,NLES_MASKS,NSV/))  ! <Sv'dp'/dz>
  ELSE
    CALL LES_ALLOCATE('X_LES_SUBGRID_USv',(/0,0,0,0/))
    CALL LES_ALLOCATE('X_LES_SUBGRID_VSv',(/0,0,0,0/))
    CALL LES_ALLOCATE('X_LES_SUBGRID_WSv',(/0,0,0,0/))
    CALL LES_ALLOCATE('X_LES_SUBGRID_Sv2',(/0,0,0,0/))
    CALL LES_ALLOCATE('X_LES_SUBGRID_SvThv',(/0,0,0,0/))
    CALL LES_ALLOCATE('X_LES_SUBGRID_W2Sv',(/0,0,0,0/))
    CALL LES_ALLOCATE('X_LES_SUBGRID_WSv2',(/0,0,0,0/))
    CALL LES_ALLOCATE('X_LES_SUBGRID_DISS_Sv2',(/0,0,0,0/))
    CALL LES_ALLOCATE('X_LES_SUBGRID_SvPz',(/0,0,0,0/))
  END IF
!
  X_LES_SUBGRID_U2  = XLES_SUBGRID_U2 
  X_LES_SUBGRID_V2  = XLES_SUBGRID_V2 
  X_LES_SUBGRID_W2  = XLES_SUBGRID_W2 
  X_LES_SUBGRID_Thl2= XLES_SUBGRID_Thl2
  X_LES_SUBGRID_UV  = XLES_SUBGRID_UV 
  X_LES_SUBGRID_WU  = XLES_SUBGRID_WU
  X_LES_SUBGRID_WV  = XLES_SUBGRID_WV
  X_LES_SUBGRID_UThl= XLES_SUBGRID_UThl
  X_LES_SUBGRID_VThl= XLES_SUBGRID_VThl
  X_LES_SUBGRID_WThl= XLES_SUBGRID_WThl
  X_LES_SUBGRID_WThv     = XLES_SUBGRID_WThv
  X_LES_SUBGRID_ThlThv   = XLES_SUBGRID_ThlThv
  X_LES_SUBGRID_W2Thl    = XLES_SUBGRID_W2Thl
  X_LES_SUBGRID_WThl2    = XLES_SUBGRID_WThl2
  X_LES_SUBGRID_DISS_Tke = XLES_SUBGRID_DISS_Tke
  X_LES_SUBGRID_DISS_Thl2= XLES_SUBGRID_DISS_Thl2
  X_LES_SUBGRID_WP       = XLES_SUBGRID_WP
  X_LES_SUBGRID_PHI3     = XLES_SUBGRID_PHI3
  X_LES_SUBGRID_LMix     = XLES_SUBGRID_LMix
  X_LES_SUBGRID_LDiss    = XLES_SUBGRID_LDiss
  X_LES_SUBGRID_Km       = XLES_SUBGRID_Km
  X_LES_SUBGRID_Kh       = XLES_SUBGRID_Kh
  X_LES_SUBGRID_ThlPz    = XLES_SUBGRID_ThlPz
  X_LES_SUBGRID_UTke= XLES_SUBGRID_UTke
  X_LES_SUBGRID_VTke= XLES_SUBGRID_VTke
  X_LES_SUBGRID_WTke= XLES_SUBGRID_WTke
  X_LES_SUBGRID_ddz_WTke  =XLES_SUBGRID_ddz_WTke
  
  X_LES_SUBGRID_THLUP_MF = XLES_SUBGRID_THLUP_MF
  X_LES_SUBGRID_RTUP_MF = XLES_SUBGRID_RTUP_MF 
  X_LES_SUBGRID_RVUP_MF = XLES_SUBGRID_RVUP_MF 
  X_LES_SUBGRID_RCUP_MF = XLES_SUBGRID_RCUP_MF 
  X_LES_SUBGRID_RIUP_MF = XLES_SUBGRID_RIUP_MF 
  X_LES_SUBGRID_WUP_MF = XLES_SUBGRID_WUP_MF  
  X_LES_SUBGRID_MASSFLUX = XLES_SUBGRID_MASSFLUX
  X_LES_SUBGRID_DETR = XLES_SUBGRID_DETR    
  X_LES_SUBGRID_ENTR = XLES_SUBGRID_ENTR    
  X_LES_SUBGRID_FRACUP = XLES_SUBGRID_FRACUP   
  X_LES_SUBGRID_THVUP_MF = XLES_SUBGRID_THVUP_MF
  X_LES_SUBGRID_WTHLMF = XLES_SUBGRID_WTHLMF     
  X_LES_SUBGRID_WRTMF = XLES_SUBGRID_WRTMF   
  X_LES_SUBGRID_WTHVMF = XLES_SUBGRID_WTHVMF   
  X_LES_SUBGRID_WUMF = XLES_SUBGRID_WUMF    
  X_LES_SUBGRID_WVMF = XLES_SUBGRID_WVMF    

  IF (LUSERV ) THEN
    X_LES_SUBGRID_Rt2  = XLES_SUBGRID_Rt2
    X_LES_SUBGRID_ThlRt= XLES_SUBGRID_ThlRt
    X_LES_SUBGRID_URt  = XLES_SUBGRID_URt
    X_LES_SUBGRID_VRt  = XLES_SUBGRID_VRt
    X_LES_SUBGRID_WRt  = XLES_SUBGRID_WRt
    X_LES_SUBGRID_RtThv   = XLES_SUBGRID_RtThv
    X_LES_SUBGRID_W2Rt    = XLES_SUBGRID_W2Rt
    X_LES_SUBGRID_WThlRt  = XLES_SUBGRID_WThlRt
    X_LES_SUBGRID_WRt2    = XLES_SUBGRID_WRt2
    X_LES_SUBGRID_DISS_Rt2= XLES_SUBGRID_DISS_Rt2
    X_LES_SUBGRID_DISS_ThlRt= XLES_SUBGRID_DISS_ThlRt
    X_LES_SUBGRID_RtPz    = XLES_SUBGRID_RtPz
    X_LES_SUBGRID_PSI3    = XLES_SUBGRID_PSI3  
  END IF
  IF (LUSERC ) THEN
    X_LES_SUBGRID_Rc2 = XLES_SUBGRID_Rc2
    X_LES_SUBGRID_URc = XLES_SUBGRID_URc
    X_LES_SUBGRID_VRc = XLES_SUBGRID_VRc
    X_LES_SUBGRID_WRc = XLES_SUBGRID_WRc
  END IF
  IF (LUSERI ) THEN
    X_LES_SUBGRID_Ri2 = XLES_SUBGRID_Ri2
  END IF
  IF (NSV>0  ) THEN
    X_LES_SUBGRID_USv = XLES_SUBGRID_USv
    X_LES_SUBGRID_VSv = XLES_SUBGRID_VSv
    X_LES_SUBGRID_WSv = XLES_SUBGRID_WSv
    X_LES_SUBGRID_Sv2      = XLES_SUBGRID_Sv2 
    X_LES_SUBGRID_SvThv    = XLES_SUBGRID_SvThv
    X_LES_SUBGRID_W2Sv     = XLES_SUBGRID_W2Sv 
    X_LES_SUBGRID_WSv2     = XLES_SUBGRID_WSv2    
    X_LES_SUBGRID_DISS_Sv2 = XLES_SUBGRID_DISS_Sv2 
    X_LES_SUBGRID_SvPz     = XLES_SUBGRID_SvPz  
  END IF
!
!
  CALL LES_ALLOCATE('X_LES_UW0',(/NLES_TIMES/))
  CALL LES_ALLOCATE('X_LES_VW0',(/NLES_TIMES/))
  CALL LES_ALLOCATE('X_LES_USTAR',(/NLES_TIMES/))
  CALL LES_ALLOCATE('X_LES_Q0',(/NLES_TIMES/))
  CALL LES_ALLOCATE('X_LES_E0',(/NLES_TIMES/))
  CALL LES_ALLOCATE('X_LES_SV0',(/NLES_TIMES,NSV/))
!
  X_LES_UW0       = XLES_UW0
  X_LES_VW0       = XLES_VW0
  X_LES_USTAR     = XLES_USTAR
  X_LES_Q0        = XLES_Q0
  X_LES_E0        = XLES_E0
  IF (NSV>0) X_LES_SV0       = XLES_SV0

ELSE
!
  XLES_RES_W_SBG_WThl         = X_LES_RES_W_SBG_WThl
  XLES_RES_W_SBG_Thl2         = X_LES_RES_W_SBG_Thl2
  XLES_RES_ddxa_U_SBG_UaU     = X_LES_RES_ddxa_U_SBG_UaU
  XLES_RES_ddxa_V_SBG_UaV     = X_LES_RES_ddxa_V_SBG_UaV
  XLES_RES_ddxa_W_SBG_UaW     = X_LES_RES_ddxa_W_SBG_UaW
  XLES_RES_ddxa_W_SBG_UaThl   = X_LES_RES_ddxa_W_SBG_UaThl
  XLES_RES_ddxa_Thl_SBG_UaW   = X_LES_RES_ddxa_Thl_SBG_UaW
  XLES_RES_ddz_Thl_SBG_W2     = X_LES_RES_ddz_Thl_SBG_W2
  XLES_RES_ddxa_Thl_SBG_UaThl = X_LES_RES_ddxa_Thl_SBG_UaThl
  IF (LUSERV) THEN
    XLES_RES_W_SBG_WRt        = X_LES_RES_W_SBG_WRt
    XLES_RES_W_SBG_Rt2        = X_LES_RES_W_SBG_Rt2
    XLES_RES_W_SBG_ThlRt      = X_LES_RES_W_SBG_ThlRt
    XLES_RES_ddxa_W_SBG_UaRt  = X_LES_RES_ddxa_W_SBG_UaRt
    XLES_RES_ddxa_Rt_SBG_UaW  = X_LES_RES_ddxa_Rt_SBG_UaW
    XLES_RES_ddz_Rt_SBG_W2    = X_LES_RES_ddz_Rt_SBG_W2
    XLES_RES_ddxa_Thl_SBG_UaRt= X_LES_RES_ddxa_Thl_SBG_UaRt
    XLES_RES_ddxa_Rt_SBG_UaThl= X_LES_RES_ddxa_Rt_SBG_UaThl
    XLES_RES_ddxa_Rt_SBG_UaRt = X_LES_RES_ddxa_Rt_SBG_UaRt
  END IF
  IF (NSV>0) THEN
    XLES_RES_ddxa_W_SBG_UaSv  = X_LES_RES_ddxa_W_SBG_UaSv
    XLES_RES_ddxa_Sv_SBG_UaW  = X_LES_RES_ddxa_Sv_SBG_UaW
    XLES_RES_ddz_Sv_SBG_W2    = X_LES_RES_ddz_Sv_SBG_W2
    XLES_RES_ddxa_Sv_SBG_UaSv = X_LES_RES_ddxa_Sv_SBG_UaSv
    XLES_RES_W_SBG_WSv        = X_LES_RES_W_SBG_WSv
    XLES_RES_W_SBG_Sv2        = X_LES_RES_W_SBG_Sv2
  END IF
  XLES_SUBGRID_U2  = X_LES_SUBGRID_U2 
  XLES_SUBGRID_V2  = X_LES_SUBGRID_V2 
  XLES_SUBGRID_W2  = X_LES_SUBGRID_W2 
  XLES_SUBGRID_Thl2= X_LES_SUBGRID_Thl2
  XLES_SUBGRID_UV  = X_LES_SUBGRID_UV 
  XLES_SUBGRID_WU  = X_LES_SUBGRID_WU
  XLES_SUBGRID_WV  = X_LES_SUBGRID_WV
  XLES_SUBGRID_UThl= X_LES_SUBGRID_UThl
  XLES_SUBGRID_VThl= X_LES_SUBGRID_VThl
  XLES_SUBGRID_WThl= X_LES_SUBGRID_WThl
  XLES_SUBGRID_WThv     = X_LES_SUBGRID_WThv
  XLES_SUBGRID_ThlThv   = X_LES_SUBGRID_ThlThv
  XLES_SUBGRID_W2Thl    = X_LES_SUBGRID_W2Thl
  XLES_SUBGRID_WThl2    = X_LES_SUBGRID_WThl2
  XLES_SUBGRID_DISS_Tke = X_LES_SUBGRID_DISS_Tke
  XLES_SUBGRID_DISS_Thl2= X_LES_SUBGRID_DISS_Thl2
  XLES_SUBGRID_WP       = X_LES_SUBGRID_WP
  XLES_SUBGRID_PHI3     = X_LES_SUBGRID_PHI3
  XLES_SUBGRID_LMix     = X_LES_SUBGRID_LMix
  XLES_SUBGRID_LDiss    = X_LES_SUBGRID_LDiss
  XLES_SUBGRID_Km       = X_LES_SUBGRID_Km
  XLES_SUBGRID_Kh       = X_LES_SUBGRID_Kh
  XLES_SUBGRID_ThlPz    = X_LES_SUBGRID_ThlPz
  XLES_SUBGRID_UTke= X_LES_SUBGRID_UTke
  XLES_SUBGRID_VTke= X_LES_SUBGRID_VTke
  XLES_SUBGRID_WTke= X_LES_SUBGRID_WTke
  XLES_SUBGRID_ddz_WTke  =X_LES_SUBGRID_ddz_WTke

  XLES_SUBGRID_THLUP_MF = X_LES_SUBGRID_THLUP_MF
  XLES_SUBGRID_RTUP_MF = X_LES_SUBGRID_RTUP_MF 
  XLES_SUBGRID_RVUP_MF = X_LES_SUBGRID_RVUP_MF 
  XLES_SUBGRID_RCUP_MF = X_LES_SUBGRID_RCUP_MF 
  XLES_SUBGRID_RIUP_MF = X_LES_SUBGRID_RIUP_MF 
  XLES_SUBGRID_WUP_MF = X_LES_SUBGRID_WUP_MF  
  XLES_SUBGRID_MASSFLUX = X_LES_SUBGRID_MASSFLUX
  XLES_SUBGRID_DETR = X_LES_SUBGRID_DETR    
  XLES_SUBGRID_ENTR = X_LES_SUBGRID_ENTR    
  XLES_SUBGRID_FRACUP = X_LES_SUBGRID_FRACUP   
  XLES_SUBGRID_THVUP_MF = X_LES_SUBGRID_THVUP_MF
  XLES_SUBGRID_WTHLMF = X_LES_SUBGRID_WTHLMF     
  XLES_SUBGRID_WRTMF = X_LES_SUBGRID_WRTMF   
  XLES_SUBGRID_WTHVMF = X_LES_SUBGRID_WTHVMF   
  XLES_SUBGRID_WUMF = X_LES_SUBGRID_WUMF    
  XLES_SUBGRID_WVMF = X_LES_SUBGRID_WVMF    
  
  IF (LUSERV ) THEN
    XLES_SUBGRID_Rt2  = X_LES_SUBGRID_Rt2
    XLES_SUBGRID_ThlRt= X_LES_SUBGRID_ThlRt
    XLES_SUBGRID_URt  = X_LES_SUBGRID_URt
    XLES_SUBGRID_VRt  = X_LES_SUBGRID_VRt
    XLES_SUBGRID_WRt  = X_LES_SUBGRID_WRt
    XLES_SUBGRID_RtThv   = X_LES_SUBGRID_RtThv
    XLES_SUBGRID_W2Rt    = X_LES_SUBGRID_W2Rt
    XLES_SUBGRID_WThlRt  = X_LES_SUBGRID_WThlRt
    XLES_SUBGRID_WRt2    = X_LES_SUBGRID_WRt2
    XLES_SUBGRID_DISS_Rt2= X_LES_SUBGRID_DISS_Rt2
    XLES_SUBGRID_DISS_ThlRt= X_LES_SUBGRID_DISS_ThlRt
    XLES_SUBGRID_RtPz    = X_LES_SUBGRID_RtPz
    XLES_SUBGRID_PSI3    = X_LES_SUBGRID_PSI3
  END IF
  IF (LUSERC ) THEN
    XLES_SUBGRID_Rc2 = X_LES_SUBGRID_Rc2
    XLES_SUBGRID_URc = X_LES_SUBGRID_URc
    XLES_SUBGRID_VRc = X_LES_SUBGRID_VRc
    XLES_SUBGRID_WRc = X_LES_SUBGRID_WRc
  END IF
  IF (LUSERI ) THEN
    XLES_SUBGRID_Ri2 = X_LES_SUBGRID_Ri2
  END IF
  IF (NSV>0  ) THEN
    XLES_SUBGRID_USv = X_LES_SUBGRID_USv
    XLES_SUBGRID_VSv = X_LES_SUBGRID_VSv
    XLES_SUBGRID_WSv = X_LES_SUBGRID_WSv
    XLES_SUBGRID_Sv2      = X_LES_SUBGRID_Sv2 
    XLES_SUBGRID_SvThv    = X_LES_SUBGRID_SvThv
    XLES_SUBGRID_W2Sv     = X_LES_SUBGRID_W2Sv 
    XLES_SUBGRID_WSv2     = X_LES_SUBGRID_WSv2    
    XLES_SUBGRID_DISS_Sv2 = X_LES_SUBGRID_DISS_Sv2 
    XLES_SUBGRID_SvPz     = X_LES_SUBGRID_SvPz  
  END IF
  XLES_UW0       = X_LES_UW0
  XLES_VW0       = X_LES_VW0
  XLES_USTAR     = X_LES_USTAR
  XLES_Q0        = X_LES_Q0
  XLES_E0        = X_LES_E0
  IF (NSV>0) XLES_SV0       = X_LES_SV0
!
  CALL LES_DEALLOCATE('X_LES_RES_W_SBG_WThl')
  CALL LES_DEALLOCATE('X_LES_RES_W_SBG_Thl2')
  CALL LES_DEALLOCATE('X_LES_RES_ddxa_U_SBG_UaU')
  CALL LES_DEALLOCATE('X_LES_RES_ddxa_V_SBG_UaV')
  CALL LES_DEALLOCATE('X_LES_RES_ddxa_W_SBG_UaW')
  CALL LES_DEALLOCATE('X_LES_RES_ddxa_W_SBG_UaThl')
  CALL LES_DEALLOCATE('X_LES_RES_ddxa_Thl_SBG_UaW')
  CALL LES_DEALLOCATE('X_LES_RES_ddz_Thl_SBG_W2')
  CALL LES_DEALLOCATE('X_LES_RES_ddxa_Thl_SBG_UaThl')
  CALL LES_DEALLOCATE('X_LES_RES_W_SBG_WRt')
  CALL LES_DEALLOCATE('X_LES_RES_W_SBG_Rt2')
  CALL LES_DEALLOCATE('X_LES_RES_W_SBG_ThlRt')
  CALL LES_DEALLOCATE('X_LES_RES_ddxa_W_SBG_UaRt')
  CALL LES_DEALLOCATE('X_LES_RES_ddxa_Rt_SBG_UaW')
  CALL LES_DEALLOCATE('X_LES_RES_ddz_Rt_SBG_W2')
  CALL LES_DEALLOCATE('X_LES_RES_ddxa_Thl_SBG_UaRt')
  CALL LES_DEALLOCATE('X_LES_RES_ddxa_Rt_SBG_UaThl')
  CALL LES_DEALLOCATE('X_LES_RES_ddxa_Rt_SBG_UaRt')
  CALL LES_DEALLOCATE('X_LES_RES_ddxa_W_SBG_UaSv')
  CALL LES_DEALLOCATE('X_LES_RES_ddxa_Sv_SBG_UaW')
  CALL LES_DEALLOCATE('X_LES_RES_ddz_Sv_SBG_W2')
  CALL LES_DEALLOCATE('X_LES_RES_ddxa_Sv_SBG_UaSv')
  CALL LES_DEALLOCATE('X_LES_RES_W_SBG_WSv')
  CALL LES_DEALLOCATE('X_LES_RES_W_SBG_Sv2')
!
  CALL LES_DEALLOCATE('X_LES_SUBGRID_U2')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_V2')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_W2')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_Thl2')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_UV')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_WU')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_WV')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_UThl')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_VThl')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_WThl')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_WThv')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_ThlThv')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_W2Thl')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_WThl2')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_DISS_Tke')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_DISS_Thl2')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_WP')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_PHI3')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_LMix')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_LDiss')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_Km')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_Kh')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_ThlPz')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_UTke')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_VTke')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_WTke')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_ddz_WTke')

  CALL LES_DEALLOCATE('X_LES_SUBGRID_THLUP_MF')  
  CALL LES_DEALLOCATE('X_LES_SUBGRID_RTUP_MF')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_RVUP_MF')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_RCUP_MF')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_RIUP_MF')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_WUP_MF')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_MASSFLUX')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_DETR')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_ENTR')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_FRACUP')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_THVUP_MF')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_WTHLMF')  
  CALL LES_DEALLOCATE('X_LES_SUBGRID_WRTMF')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_WTHVMF')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_WUMF')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_WVMF')
  
  CALL LES_DEALLOCATE('X_LES_SUBGRID_Rt2')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_ThlRt')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_URt')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_VRt')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_WRt')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_RtThv')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_W2Rt')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_WThlRt')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_WRt2')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_DISS_Rt2')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_DISS_ThlRt')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_RtPz')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_PSI3')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_Rc2')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_URc')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_VRc')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_WRc')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_Ri2')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_USv')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_VSv')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_WSv')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_Sv2')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_SvThv')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_W2Sv')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_WSv2')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_DISS_Sv2')
  CALL LES_DEALLOCATE('X_LES_SUBGRID_SvPz')
  !
  CALL LES_DEALLOCATE('X_LES_UW0')
  CALL LES_DEALLOCATE('X_LES_VW0')
  CALL LES_DEALLOCATE('X_LES_USTAR')
  CALL LES_DEALLOCATE('X_LES_Q0')
  CALL LES_DEALLOCATE('X_LES_E0')
  CALL LES_DEALLOCATE('X_LES_SV0')
!
END IF
!
CALL SECOND_MNH(ZTIME2)
!
XTIME_LES = XTIME_LES + ZTIME2 - ZTIME1
!
END SUBROUTINE SWITCH_SBG_LES_n
