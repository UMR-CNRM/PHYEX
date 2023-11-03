!MNH_LIC Copyright 2012-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ######spl
     MODULE MODI_ICE_ADJUST_BIS
!    ###############################
!
INTERFACE
!
!     #################################################################
      SUBROUTINE ICE_ADJUST_BIS(PP,PTH,PR)
!     #################################################################
!
!!               
!*               1.1  Declaration of Arguments
!! 

REAL, DIMENSION(:,:,:),  INTENT(IN)     :: PP     ! Pressure
REAL, DIMENSION(:,:,:),  INTENT(INOUT)  :: PTH    ! thetal to transform into th
REAL, DIMENSION(:,:,:,:),INTENT(INOUT)  :: PR     ! Total mixing ratios to transform into rv,rc and ri
!
END SUBROUTINE ICE_ADJUST_BIS

END INTERFACE
!
END MODULE MODI_ICE_ADJUST_BIS
!     ######spl
      SUBROUTINE ICE_ADJUST_BIS(PP,PTH,PR)
!     #################################################################
!
!
!!****  *ICE_ADJUST_BIS* - computes an adjusted state of thermodynamical variables
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      Valery Masson & C. Lac * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original         09/2012
!!      M.Moge           08/2015 UPDATE_HALO_ll on PTH, ZRV, ZRC, ZRI
!  P. Wautelet 20/05/2019: add name argument to ADDnFIELD_ll + new ADD4DFIELD_ll subroutine
!!
!! --------------------------------------------------------------------------
!
!*      0. DECLARATIONS
!          ------------
!
USE MODD_CST,   ONLY: XCPD, XRD, XP00, CST
USE MODD_NEB_n, ONLY: NEBN
!
USE MODI_COMPUTE_FUNCTION_THERMO
USE MODI_THLRT_FROM_THRVRCRI
!
USE MODE_ll
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
REAL, DIMENSION(:,:,:),  INTENT(IN)     :: PP    ! Pressure
REAL, DIMENSION(:,:,:),  INTENT(INOUT)  :: PTH   ! thetal to transform into th
REAL, DIMENSION(:,:,:,:),INTENT(INOUT)  :: PR    ! Total mixing ratios to transform into rv,rc and ri
!
!-------------------------------------------------------------------------------
!
!       0.2  declaration of local variables
REAL, DIMENSION(SIZE(PTH,1),SIZE(PTH,2),SIZE(PTH,3)) :: ZTHL, ZRW, ZRV, ZRC, &
                                                        ZRI, ZWORK
REAL, DIMENSION(SIZE(PTH,1),SIZE(PTH,2),SIZE(PTH,3)) :: ZFRAC_ICE, ZRSATW, ZRSATI
REAL, DIMENSION(SIZE(PTH,1),SIZE(PTH,2),SIZE(PTH,3)) :: ZT, ZEXN, ZLVOCPEXN,ZLSOCPEXN
REAL, DIMENSION(SIZE(PTH,1),SIZE(PTH,2),SIZE(PTH,3), 16) :: ZBUF
INTEGER :: IRR
CHARACTER(LEN=1) :: YFRAC_ICE
!
INTEGER :: IINFO_ll
TYPE(LIST_ll), POINTER :: TZFIELDS_ll=>NULL()   ! list of fields to exchange
!----------------------------------------------------------------------------
!
!*      1 Initialisation
!         --------------
!
IRR = SIZE(PR,4)
!
ZRV(:,:,:)=0.
IF (IRR>=1) &
ZRV(:,:,:)=PR(:,:,:,1)
ZRC(:,:,:)=0.
IF (IRR>=2) &
ZRC(:,:,:)=PR(:,:,:,2)
ZRI(:,:,:)=0.
IF (IRR>=4) &
ZRI(:,:,:)=PR(:,:,:,4)
!
YFRAC_ICE='T'
ZFRAC_ICE(:,:,:) = 0.
!
!*      2 Computation
!         -----------
!
ZEXN(:,:,:)=(PP(:,:,:)/XP00)**(XRD/XCPD)
!
CALL COMPUTE_FUNCTION_THERMO( IRR,                                &
                              PTH, PR, ZEXN, PP,                  &
                              ZT,ZLVOCPEXN,ZLSOCPEXN              )

!
CALL THLRT_FROM_THRVRCRI( IRR, PTH, PR, ZLVOCPEXN, ZLSOCPEXN,&
                          ZTHL, ZRW                          )
!
CALL TH_R_FROM_THL_RT(CST, NEBN, SIZE(ZFRAC_ICE), YFRAC_ICE,ZFRAC_ICE(:,:,:),PP(:,:,:), &
                         ZTHL(:,:,:), ZRW(:,:,:), PTH(:,:,:),  &
                         ZRV(:,:,:), ZRC(:,:,:), ZRI(:,:,:),   &
                         ZRSATW(:,:,:), ZRSATI(:,:,:),OOCEAN=.FALSE.,&
                         PBUF=ZBUF)
CALL ADD3DFIELD_ll( TZFIELDS_ll, PTH, 'ICE_ADJUST_BIS::PTH')
IF (IRR>=1) THEN
  CALL ADD3DFIELD_ll( TZFIELDS_ll, ZRV, 'ICE_ADJUST_BIS::ZRV' )
ENDIF
IF (IRR>=2) THEN
  CALL ADD3DFIELD_ll( TZFIELDS_ll, ZRC, 'ICE_ADJUST_BIS::ZRC' )
ENDIF
IF (IRR>=4) THEN
  CALL ADD3DFIELD_ll( TZFIELDS_ll, ZRI, 'ICE_ADJUST_BIS::ZRI' )
ENDIF
CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
CALL CLEANLIST_ll(TZFIELDS_ll)
!

IF (IRR>=1) &
PR(:,:,:,1) = ZRV(:,:,:)
IF (IRR>=2) &
PR(:,:,:,2) = ZRC(:,:,:)
IF (IRR>=4) &
PR(:,:,:,4) = ZRI(:,:,:)
!
CONTAINS
INCLUDE "th_r_from_thl_rt.func.h"
INCLUDE "compute_frac_ice.func.h"
END SUBROUTINE ICE_ADJUST_BIS
