!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!    ################################ 
     MODULE MODI_TURB_HOR_THERMO_FLUX
!    ################################ 
!
INTERFACE  
!
      SUBROUTINE TURB_HOR_THERMO_FLUX(KSPLT, KRR, KRRL, KRRI,        &
                      OTURB_FLX,OSUBG_COND,                          &
                      TPFILE,                                        &
                      PK,PINV_PDXX,PINV_PDYY,PINV_PDZZ,PMZM_PRHODJ,  &
                      PDXX,PDYY,PDZZ,PDZX,PDZY,                      &
                      PDIRCOSXW,PDIRCOSYW,                           &
                      PRHODJ,                                        &
                      PSFTHM,PSFRM,                                  &
                      PWM,PTHLM,PRM,                                 &
                      PATHETA,PAMOIST,PSRCM,PFRAC_ICE,               &
                      PRTHLS,PRRS                                    )
!
USE MODD_IO, ONLY: TFILEDATA
!
INTEGER,                  INTENT(IN)    :: KSPLT         ! split process index
INTEGER,                  INTENT(IN)    :: KRR           ! number of moist var.
INTEGER,                  INTENT(IN)    :: KRRL          ! number of liquid water var.
INTEGER,                  INTENT(IN)    :: KRRI          ! number of ice water var.
LOGICAL,                  INTENT(IN)    ::  OTURB_FLX    ! switch to write the
                                 ! turbulent fluxes in the syncronous FM-file
LOGICAL,                 INTENT(IN)  ::   OSUBG_COND ! Switch for sub-grid 
!                                                    condensation
TYPE(TFILEDATA),          INTENT(IN)    ::  TPFILE       ! Output file
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PK          ! Turbulent diffusion doef.
                                                        ! PK = PLM * SQRT(PTKEM)
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PINV_PDXX   ! 1./PDXX
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PINV_PDYY   ! 1./PDYY
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PINV_PDZZ   ! 1./PDZZ
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PMZM_PRHODJ ! MZM(PRHODJ)
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PDXX, PDYY, PDZZ, PDZX, PDZY 
                                                         ! Metric coefficients
REAL, DIMENSION(:,:),     INTENT(IN)    ::  PDIRCOSXW, PDIRCOSYW
! Director Cosinus along x, y and z directions at surface w-point
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PRHODJ       ! density * grid volume
!
REAL, DIMENSION(:,:),     INTENT(IN)    ::  PSFTHM,PSFRM
!
! Variables at t-1
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PWM 
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PTHLM 
REAL, DIMENSION(:,:,:,:), INTENT(IN)    ::  PRM          ! mixing ratios at t-1,
                              !  where PRM(:,:,:,1) = conservative mixing ratio
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PATHETA      ! coefficients between 
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PAMOIST      ! s and Thetal and Rnp

REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PSRCM
                                  ! normalized 2nd-order flux
                                  ! s'r'c/2Sigma_s2 at t-1 multiplied by Lambda_3
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PFRAC_ICE    ! ri fraction of rc+ri
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) ::  PRTHLS
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) ::  PRRS         ! var. at t+1 -split-
!
!
END SUBROUTINE TURB_HOR_THERMO_FLUX
!
END INTERFACE
!
END MODULE MODI_TURB_HOR_THERMO_FLUX
!     ################################################################
      SUBROUTINE TURB_HOR_THERMO_FLUX(KSPLT, KRR, KRRL, KRRI,        &
                      OTURB_FLX,OSUBG_COND,                          &
                      TPFILE,                                        &
                      PK,PINV_PDXX,PINV_PDYY,PINV_PDZZ,PMZM_PRHODJ,  &
                      PDXX,PDYY,PDZZ,PDZX,PDZY,                      &
                      PDIRCOSXW,PDIRCOSYW,                           &
                      PRHODJ,                                        &
                      PSFTHM,PSFRM,                                  &
                      PWM,PTHLM,PRM,                                 &
                      PATHETA,PAMOIST,PSRCM,PFRAC_ICE,               &
                      PRTHLS,PRRS                                    )
!     ################################################################
!
!
!!****  *TURB_HOR* -routine to compute the source terms in the meso-NH
!!               model equations due to the non-vertical turbulent fluxes.
!!
!!    PURPOSE
!!    -------
!!
!!     see TURB_HOR
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!
!!      Joan Cuxart             * INM and Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!                     Aug    , 1997 (V. Saravane) spliting of TURB_HOR
!!                     Nov  27, 1997 (V. Masson) clearing of the routine
!!                     Feb. 18, 1998 (J. Stein) bug for v'RC'
!!                     Oct  18, 2000 (V. Masson) LES computations + LFLAT switch
!!                     Nov  06, 2002 (V. Masson) LES budgets
!!                     Feb  20, 2003 (JP Pinty)  Add PFRAC_ICE
!!                     October 2009 (G. Tanguy) add ILENCH=LEN(YCOMMENT) after
!!                                              change of YCOMMENT
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!! --------------------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
!
USE MODD_CST
USE MODD_CONF
USE MODD_CTURB
use modd_field,          only: tfielddata, TYPEREAL
USE MODD_IO,             ONLY: TFILEDATA
USE MODD_PARAMETERS
USE MODD_LES
!
USE MODE_IO_FIELD_WRITE, only: IO_Field_write
!
USE MODI_GRADIENT_M
USE MODI_GRADIENT_U
USE MODI_GRADIENT_V
USE MODI_GRADIENT_W
USE MODI_SHUMAN 
USE MODI_LES_MEAN_SUBGRID
!!USE MODI_EMOIST
!!USE MODI_ETHETA
!
USE MODI_SECOND_MNH
!
IMPLICIT NONE
!
!
!*       0.1  declaration of arguments
!
!
!
INTEGER,                  INTENT(IN)    :: KSPLT         ! split process index
INTEGER,                  INTENT(IN)    :: KRR           ! number of moist var.
INTEGER,                  INTENT(IN)    :: KRRL          ! number of liquid water var.
INTEGER,                  INTENT(IN)    :: KRRI          ! number of ice water var.
LOGICAL,                  INTENT(IN)    ::  OTURB_FLX    ! switch to write the
                                 ! turbulent fluxes in the syncronous FM-file
LOGICAL,                 INTENT(IN)  ::   OSUBG_COND ! Switch for sub-grid 
!                                                    condensation
TYPE(TFILEDATA),          INTENT(IN)    ::  TPFILE       ! Output file
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PK          ! Turbulent diffusion doef.
                                                        ! PK = PLM * SQRT(PTKEM)
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PINV_PDXX   ! 1./PDXX
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PINV_PDYY   ! 1./PDYY
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PINV_PDZZ   ! 1./PDZZ
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PMZM_PRHODJ ! MZM(PRHODJ)
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PDXX, PDYY, PDZZ, PDZX, PDZY 
                                                         ! Metric coefficients
REAL, DIMENSION(:,:),     INTENT(IN)    ::  PDIRCOSXW, PDIRCOSYW
! Director Cosinus along x, y and z directions at surface w-point
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PRHODJ       ! density * grid volume
!
REAL, DIMENSION(:,:),     INTENT(IN)    ::  PSFTHM,PSFRM
!
! Variables at t-1
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PWM 
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PTHLM 
REAL, DIMENSION(:,:,:,:), INTENT(IN)    ::  PRM          ! mixing ratios at t-1,
                              !  where PRM(:,:,:,1) = conservative mixing ratio
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PATHETA      ! coefficients between 
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PAMOIST      ! s and Thetal and Rnp

REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PSRCM
                                  ! normalized 2nd-order flux
                                  ! s'r'c/2Sigma_s2 at t-1 multiplied by Lambda_3
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PFRAC_ICE    ! ri fraction of rc+ri
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) ::  PRTHLS
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) ::  PRRS         ! var. at t+1 -split-
!
!
!
!*       0.2  declaration of local variables
!
REAL, DIMENSION(SIZE(PTHLM,1),SIZE(PTHLM,2),SIZE(PTHLM,3))       &
                                     :: ZFLX,ZFLXC
    ! work arrays
!   
!! REAL, DIMENSION(SIZE(PTHLM,1),SIZE(PTHLM,2),SIZE(PTHLM,3))  :: ZVPTV
INTEGER             :: IKB,IKE,IKU
                                    ! Index values for the Beginning and End
                                    ! mass points of the domain  
REAL, DIMENSION(SIZE(PDZZ,1),SIZE(PDZZ,2),1+JPVEXT:3+JPVEXT) :: ZCOEFF 
                                    ! coefficients for the uncentred gradient 
                                    ! computation near the ground
!
REAL :: ZTIME1, ZTIME2
TYPE(TFIELDDATA) :: TZFIELD
! ---------------------------------------------------------------------------
!
!*       1.   PRELIMINARY COMPUTATIONS
!             ------------------------
!
IKB = 1+JPVEXT               
IKE = SIZE(PTHLM,3)-JPVEXT    
IKU = SIZE(PTHLM,3)
!
!
!  compute the coefficients for the uncentred gradient computation near the 
!  ground
ZCOEFF(:,:,IKB+2)= - PDZZ(:,:,IKB+1) /      &
       ( (PDZZ(:,:,IKB+2)+PDZZ(:,:,IKB+1)) * PDZZ(:,:,IKB+2) )
ZCOEFF(:,:,IKB+1)=   (PDZZ(:,:,IKB+2)+PDZZ(:,:,IKB+1)) /      &
       ( PDZZ(:,:,IKB+1) * PDZZ(:,:,IKB+2) )
ZCOEFF(:,:,IKB)= - (PDZZ(:,:,IKB+2)+2.*PDZZ(:,:,IKB+1)) /      &
       ( (PDZZ(:,:,IKB+2)+PDZZ(:,:,IKB+1)) * PDZZ(:,:,IKB+1) )
!
!*       2.   < U' THETA'l >
!             --------------
!
! 
ZFLX(:,:,:)     = -XCSHF * MXM( PK ) * GX_M_U(1,IKU,1,PTHLM,PDXX,PDZZ,PDZX)
ZFLX(:,:,IKE+1) = ZFLX(:,:,IKE) 
!
! Compute the flux at the first inner U-point with an uncentred vertical  
! gradient
ZFLX(:,:,IKB:IKB) = -XCSHF * MXM( PK(:,:,IKB:IKB) ) *          &
  ( DXM(PTHLM(:,:,IKB:IKB)) * PINV_PDXX(:,:,IKB:IKB)           &
   -MXM( ZCOEFF(:,:,IKB+2:IKB+2)*PTHLM(:,:,IKB+2:IKB+2)        &
         +ZCOEFF(:,:,IKB+1:IKB+1)*PTHLM(:,:,IKB+1:IKB+1)       &
         +ZCOEFF(:,:,IKB  :IKB  )*PTHLM(:,:,IKB  :IKB  ))      &
        *0.5* ( PDZX(:,:,IKB+1:IKB+1)+PDZX(:,:,IKB:IKB))       &
        * PINV_PDXX(:,:,IKB:IKB) )
! extrapolates the flux under the ground so that the vertical average with 
! the IKB flux gives the ground value  ( warning the tangential surface
! flux has been set to 0 for the moment !!  to be improved )
ZFLX(:,:,IKB-1:IKB-1) = 2. * MXM(  SPREAD( PSFTHM(:,:)* PDIRCOSXW(:,:), 3,1) )  &
                       - ZFLX(:,:,IKB:IKB)
!
! Add this source to the Theta_l sources
!
IF (.NOT. LFLAT) THEN
  PRTHLS(:,:,:) =  PRTHLS                                                   &
                - DXF( MXM(PRHODJ) * ZFLX * PINV_PDXX )                          &
                + DZF( PMZM_PRHODJ *MXF(PDZX*(MZM(ZFLX * PINV_PDXX))) * PINV_PDZZ )
ELSE
  PRTHLS(:,:,:) =  PRTHLS - DXF( MXM(PRHODJ) * ZFLX * PINV_PDXX )
END IF
!
! Compute the equivalent tendancy for Rc and Ri
!
IF ( KRRL >= 1 ) THEN
  IF (.NOT. LFLAT) THEN
    ZFLXC = 2.*( MXF( MXM( PRHODJ*PATHETA*PSRCM )*ZFLX )                       &
                +MZF( MZM( PRHODJ*PATHETA*PSRCM )*MXF(                         &
                                               PDZX*(MZM( ZFLX*PINV_PDXX )) ) )&
               )
    IF ( KRRI >= 1 ) THEN
      PRRS(:,:,:,2) = PRRS(:,:,:,2) +  2. *                                    &
        (- DXF( MXM( PRHODJ*PATHETA*PSRCM )*ZFLX*PINV_PDXX )                   &
         + DZF( MZM( PRHODJ*PATHETA*PSRCM )*MXF( PDZX*(MZM( ZFLX*PINV_PDXX )) )&
                                           *PINV_PDZZ )                        &
        )*(1.0-PFRAC_ICE(:,:,:))
      PRRS(:,:,:,4) = PRRS(:,:,:,4) +  2. *                                    &
        (- DXF( MXM( PRHODJ*PATHETA*PSRCM )*ZFLX*PINV_PDXX )                   &
         + DZF( MZM( PRHODJ*PATHETA*PSRCM )*MXF( PDZX*(MZM( ZFLX*PINV_PDXX )) )&
                                           *PINV_PDZZ )                        &
        )*PFRAC_ICE(:,:,:)
    ELSE
      PRRS(:,:,:,2) = PRRS(:,:,:,2) +  2. *                                    &
        (- DXF( MXM( PRHODJ*PATHETA*PSRCM )*ZFLX*PINV_PDXX )                   &
         + DZF( MZM( PRHODJ*PATHETA*PSRCM )*MXF( PDZX*(MZM( ZFLX*PINV_PDXX )) )&
                                           *PINV_PDZZ )                        &
        )
    END IF
  ELSE
    ZFLXC = 2.*MXF( MXM( PRHODJ*PATHETA*PSRCM )*ZFLX )
    IF ( KRRI >= 1 ) THEN
      PRRS(:,:,:,2) = PRRS(:,:,:,2) -  2. *                                    &
        DXF( MXM( PRHODJ*PATHETA*PSRCM )*ZFLX*PINV_PDXX )*(1.0-PFRAC_ICE(:,:,:))
      PRRS(:,:,:,4) = PRRS(:,:,:,4) -  2. *                                    &
        DXF( MXM( PRHODJ*PATHETA*PSRCM )*ZFLX*PINV_PDXX )*PFRAC_ICE(:,:,:)
    ELSE
      PRRS(:,:,:,2) = PRRS(:,:,:,2) -  2. *                                    &
        DXF( MXM( PRHODJ*PATHETA*PSRCM )*ZFLX*PINV_PDXX )
    END IF
  END IF
END IF
!
!! stores this flux in ZWORK to compute later <U' VPT'>
!!ZWORK(:,:,:) = ZFLX(:,:,:) 
!
! stores the horizontal  <U THl>
IF ( tpfile%lopened .AND. OTURB_FLX ) THEN
  TZFIELD%CMNHNAME   = 'UTHL_FLX'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'UTHL_FLX'
  TZFIELD%CUNITS     = 'K m s-1'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Z_UTHL_FLX'
  TZFIELD%NGRID      = 2
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZFLX)
END IF
!
IF (KSPLT==1 .AND. LLES_CALL) THEN
  CALL SECOND_MNH(ZTIME1)
  CALL LES_MEAN_SUBGRID( MXF(ZFLX), X_LES_SUBGRID_UThl ) 
  CALL LES_MEAN_SUBGRID( MZF(MXF(GX_W_UW(PWM,PDXX,PDZZ,PDZX)*MZM(ZFLX))),&
                         X_LES_RES_ddxa_W_SBG_UaThl , .TRUE. )
  CALL LES_MEAN_SUBGRID( GX_M_M(PTHLM,PDXX,PDZZ,PDZX)*MXF(ZFLX),&
                         X_LES_RES_ddxa_Thl_SBG_UaThl , .TRUE. )
  IF (KRR>=1) THEN
    CALL LES_MEAN_SUBGRID( GX_M_M(PRM(:,:,:,1),PDXX,PDZZ,PDZX)*MXF(ZFLX), &
                           X_LES_RES_ddxa_Rt_SBG_UaThl , .TRUE. )
  END IF
  CALL SECOND_MNH(ZTIME2)
  XTIME_LES = XTIME_LES + ZTIME2 - ZTIME1
END IF
!
!*       3.   < U' R'np >
!             -----------
IF (KRR/=0) THEN
  !
  ZFLX(:,:,:)     = -XCHF * MXM( PK ) * GX_M_U(1,IKU,1,PRM(:,:,:,1),PDXX,PDZZ,PDZX)
  ZFLX(:,:,IKE+1) = ZFLX(:,:,IKE) 
!
! Compute the flux at the first inner U-point with an uncentred vertical  
! gradient
  ZFLX(:,:,IKB:IKB) = -XCHF * MXM( PK(:,:,IKB:IKB) ) *           &
    ( DXM(PRM(:,:,IKB:IKB,1)) * PINV_PDXX(:,:,IKB:IKB)           &
     -MXM( ZCOEFF(:,:,IKB+2:IKB+2)*PRM(:,:,IKB+2:IKB+2,1)        &
           +ZCOEFF(:,:,IKB+1:IKB+1)*PRM(:,:,IKB+1:IKB+1,1)       &
           +ZCOEFF(:,:,IKB  :IKB  )*PRM(:,:,IKB  :IKB  ,1))      &
          *0.5* ( PDZX(:,:,IKB+1:IKB+1)+PDZX(:,:,IKB:IKB))       &
          * PINV_PDXX(:,:,IKB:IKB) )
! extrapolates the flux under the ground so that the vertical average with 
! the IKB flux gives the ground value  ( warning the tangential surface
! flux has been set to 0 for the moment !!  to be improved )
  ZFLX(:,:,IKB-1:IKB-1) = 2. * MXM(  SPREAD( PSFRM(:,:)* PDIRCOSXW(:,:), 3,1) ) &
                       - ZFLX(:,:,IKB:IKB)
  !
  ! Add this source to the conservative mixing ratio sources
  !
  IF (.NOT. LFLAT) THEN
    PRRS(:,:,:,1) = PRRS(:,:,:,1)                                             &
                  - DXF( MXM(PRHODJ) * ZFLX * PINV_PDXX )                          &
                  + DZF( PMZM_PRHODJ *MXF(PDZX*(MZM(ZFLX * PINV_PDXX))) * PINV_PDZZ )
  ELSE
    PRRS(:,:,:,1) = PRRS(:,:,:,1) - DXF( MXM(PRHODJ) * ZFLX * PINV_PDXX )
  END IF
  !
  ! Compute the equivalent tendancy for Rc and Ri
  !
  IF ( KRRL >= 1 ) THEN
    IF (.NOT. LFLAT) THEN
      ZFLXC = ZFLXC            &
            + 2.*( MXF( MXM( PRHODJ*PAMOIST*PSRCM )*ZFLX )                     &
                  +MZF( MZM( PRHODJ*PAMOIST*PSRCM )*MXF(                       &
                                               PDZX*(MZM( ZFLX*PINV_PDXX )) ) )&
                 )
      IF ( KRRI >= 1 ) THEN
        PRRS(:,:,:,2) = PRRS(:,:,:,2) +  2. *                                  &
        (- DXF( MXM( PRHODJ*PAMOIST*PSRCM )*ZFLX*PINV_PDXX )                   &
         + DZF( MZM( PRHODJ*PAMOIST*PSRCM )*MXF( PDZX*(MZM( ZFLX*PINV_PDXX )) )&
                                           *PINV_PDZZ )                        &
        )*(1.0-PFRAC_ICE(:,:,:))
        PRRS(:,:,:,2) = PRRS(:,:,:,2) +  2. *                                  &
        (- DXF( MXM( PRHODJ*PAMOIST*PSRCM )*ZFLX*PINV_PDXX )                   &
         + DZF( MZM( PRHODJ*PAMOIST*PSRCM )*MXF( PDZX*(MZM( ZFLX*PINV_PDXX )) )&
                                           *PINV_PDZZ )                        &
        )*PFRAC_ICE(:,:,:)
      ELSE
        PRRS(:,:,:,2) = PRRS(:,:,:,2) +  2. *                                  &
        (- DXF( MXM( PRHODJ*PAMOIST*PSRCM )*ZFLX*PINV_PDXX )                   &
         + DZF( MZM( PRHODJ*PAMOIST*PSRCM )*MXF( PDZX*(MZM( ZFLX*PINV_PDXX )) )&
                                           *PINV_PDZZ )                        &
        )
      END IF
    ELSE
      ZFLXC = ZFLXC + 2.*MXF( MXM( PRHODJ*PAMOIST*PSRCM )*ZFLX )
      IF ( KRRI >= 1 ) THEN
        PRRS(:,:,:,2) = PRRS(:,:,:,2) -  2. *                                  &
        DXF( MXM( PRHODJ*PAMOIST*PSRCM )*ZFLX*PINV_PDXX )*(1.0-PFRAC_ICE(:,:,:))
        PRRS(:,:,:,4) = PRRS(:,:,:,4) -  2. *                                  &
        DXF( MXM( PRHODJ*PAMOIST*PSRCM )*ZFLX*PINV_PDXX )*PFRAC_ICE(:,:,:)
      ELSE
        PRRS(:,:,:,2) = PRRS(:,:,:,2) -  2. *                                  &
        DXF( MXM( PRHODJ*PAMOIST*PSRCM )*ZFLX*PINV_PDXX )
      END IF
    END IF
  END IF
  !
  ! stores the horizontal  <U Rnp>
  IF ( tpfile%lopened .AND. OTURB_FLX ) THEN
    TZFIELD%CMNHNAME   = 'UR_FLX'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'UR_FLX'
    TZFIELD%CUNITS     = 'kg kg-1 m s-1'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'X_Y_Z_UR_FLX'
    TZFIELD%NGRID      = 2
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,ZFLX)
  END IF
  !
  IF (KSPLT==1 .AND. LLES_CALL) THEN
    CALL SECOND_MNH(ZTIME1)
    CALL LES_MEAN_SUBGRID( MXF(ZFLX), X_LES_SUBGRID_URt ) 
    CALL LES_MEAN_SUBGRID( MZF(MXF(GX_W_UW(PWM,PDXX,PDZZ,PDZX)*MZM(ZFLX))),&
                           X_LES_RES_ddxa_W_SBG_UaRt , .TRUE. )
    CALL LES_MEAN_SUBGRID( GX_M_M(PTHLM,PDXX,PDZZ,PDZX)*MXF(ZFLX),&
                           X_LES_RES_ddxa_Thl_SBG_UaRt , .TRUE. )
    CALL LES_MEAN_SUBGRID( GX_M_M(PRM(:,:,:,1),PDXX,PDZZ,PDZX)*MXF(ZFLX),&
                           X_LES_RES_ddxa_Rt_SBG_UaRt , .TRUE. )
    CALL SECOND_MNH(ZTIME2)
    XTIME_LES = XTIME_LES + ZTIME2 - ZTIME1
  END IF
!
  !
  IF (KRRL>0 .AND. KSPLT==1 .AND. LLES_CALL) THEN
    CALL SECOND_MNH(ZTIME1)
    CALL LES_MEAN_SUBGRID(MXF(ZFLXC), X_LES_SUBGRID_URc )
    CALL SECOND_MNH(ZTIME2)
    XTIME_LES = XTIME_LES + ZTIME2 - ZTIME1
  END IF
!
END IF 
!
!*       4.   < U' TPV' >
!             -----------
!
!! to be tested later
!!IF (KRR/=0) THEN
!!  ! here ZFLX= <U'Rnp'> and ZWORK= <U'Thetal'>
!!  !
!!  ZVPTU(:,:,:) =                                                        &
!!    ZWORK(:,:,:)*MXM(ETHETA(KRR,KRRI,PTHLT,PEXNREF,PRT,PLOCPT,PSRCM)) +       &
!!     ZFLX(:,:,:)*MXM(EMOIST(KRR,KRRI,PTHLT,PEXNREF,PRT,PLOCPT,PSRCM))
!!  !
!!  ! stores the horizontal  <U VPT>
!!  IF ( tpfile%lopened .AND. OTURB_FLX ) THEN
!!    TZFIELD%CMNHNAME   = 'UVPT_FLX'
!!    TZFIELD%CSTDNAME   = ''
!!    TZFIELD%CLONGNAME  = 'UVPT_FLX'
!!    TZFIELD%CUNITS     = 'K m s-1'
!!    TZFIELD%CDIR       = 'XY'
!!    TZFIELD%CCOMMENT   = 'X_Y_Z_UVPT_FLX'
!!    TZFIELD%NGRID      = 2
!!    TZFIELD%NTYPE      = TYPEREAL
!!    TZFIELD%NDIMS      = 3
!!    TZFIELD%LTIMEDEP   = .TRUE.
!!    CALL IO_Field_write(TPFILE,TZFIELD,ZVPTU)
!!  END IF
!!!
!!ELSE
!!  ZVPTU(:,:,:)=ZWORK(:,:,:)
!!END IF
!
!
!*       5.   < V' THETA'l >
!             --------------
!
!
IF (.NOT. L2D) THEN
  ZFLX(:,:,:)     = -XCSHF * MYM( PK ) * GY_M_V(1,IKU,1,PTHLM,PDYY,PDZZ,PDZY)
  ZFLX(:,:,IKE+1) = ZFLX(:,:,IKE) 
ELSE
  ZFLX(:,:,:)     = 0.
END IF
!
!
! Compute the flux at the first inner U-point with an uncentred vertical  
! gradient
ZFLX(:,:,IKB:IKB) = -XCSHF * MYM( PK(:,:,IKB:IKB) ) *          &
  ( DYM(PTHLM(:,:,IKB:IKB)) * PINV_PDYY(:,:,IKB:IKB)           &
   -MYM( ZCOEFF(:,:,IKB+2:IKB+2)*PTHLM(:,:,IKB+2:IKB+2)        &
         +ZCOEFF(:,:,IKB+1:IKB+1)*PTHLM(:,:,IKB+1:IKB+1)       &
         +ZCOEFF(:,:,IKB  :IKB  )*PTHLM(:,:,IKB  :IKB  ) )     &
        *0.5* ( PDZY(:,:,IKB+1:IKB+1)+PDZY(:,:,IKB:IKB))       &
        * PINV_PDYY(:,:,IKB:IKB) )
! extrapolates the flux under the ground so that the vertical average with 
! the IKB flux gives the ground value  ( warning the tangential surface
! flux has been set to 0 for the moment !!  to be improved )
ZFLX(:,:,IKB-1:IKB-1) = 2. * MYM(  SPREAD( PSFTHM(:,:)* PDIRCOSYW(:,:), 3,1) ) &
                       - ZFLX(:,:,IKB:IKB)
!
! Add this source to the Theta_l sources
!
IF (.NOT. L2D) THEN 
  IF (.NOT. LFLAT) THEN
    PRTHLS(:,:,:) =  PRTHLS                                                         &
                  - DYF( MYM(PRHODJ) * ZFLX * PINV_PDYY )                           &
                  + DZF( PMZM_PRHODJ *MYF(PDZY*(MZM(ZFLX * PINV_PDYY))) * PINV_PDZZ )
  ELSE
    PRTHLS(:,:,:) =  PRTHLS - DYF( MYM(PRHODJ) * ZFLX * PINV_PDYY )
  END IF
END IF
!
! Compute the equivalent tendancy for Rc and Ri
!
!IF ( OSUBG_COND .AND. KRRL > 0 .AND. .NOT. L2D) THEN
IF ( KRRL >= 1 .AND. .NOT. L2D) THEN
  IF (.NOT. LFLAT) THEN
    ZFLXC = 2.*( MYF( MYM( PRHODJ*PATHETA*PSRCM )*ZFLX )                       &
                +MZF( MZM( PRHODJ*PATHETA*PSRCM )*MYF(                         &
                                               PDZY*(MZM( ZFLX*PINV_PDYY )) ) )&
               )
    IF ( KRRI >= 1 ) THEN
      PRRS(:,:,:,2) = PRRS(:,:,:,2) + 2. *                                     &
        (- DYF( MYM( PRHODJ*PATHETA*PSRCM )*ZFLX*PINV_PDYY )                   &
         + DZF( MZM( PRHODJ*PATHETA*PSRCM )*MYF( PDZY*(MZM( ZFLX*PINV_PDYY )) )&
                                           *PINV_PDZZ )                        &
        )*(1.0-PFRAC_ICE(:,:,:))
      PRRS(:,:,:,4) = PRRS(:,:,:,4) + 2. *                                     &
        (- DYF( MYM( PRHODJ*PATHETA*PSRCM )*ZFLX*PINV_PDYY )                   &
         + DZF( MZM( PRHODJ*PATHETA*PSRCM )*MYF( PDZY*(MZM( ZFLX*PINV_PDYY )) )&
                                           *PINV_PDZZ )                        &
        )*PFRAC_ICE(:,:,:)
    ELSE
      PRRS(:,:,:,2) = PRRS(:,:,:,2) + 2. *                                     &
        (- DYF( MYM( PRHODJ*PATHETA*PSRCM )*ZFLX*PINV_PDYY )                   &
         + DZF( MZM( PRHODJ*PATHETA*PSRCM )*MYF( PDZY*(MZM( ZFLX*PINV_PDYY )) )&
                                           *PINV_PDZZ )                        &
        )
    END IF
  ELSE
    ZFLXC = 2.*MYF( MYM( PRHODJ*PATHETA*PSRCM )*ZFLX )
    IF ( KRRI >= 1 ) THEN
      PRRS(:,:,:,2) = PRRS(:,:,:,2) - 2. *                                     &
        DYF( MYM( PRHODJ*PATHETA*PSRCM )*ZFLX*PINV_PDYY )*(1.0-PFRAC_ICE(:,:,:))
      PRRS(:,:,:,4) = PRRS(:,:,:,4) - 2. *                                     &
        DYF( MYM( PRHODJ*PATHETA*PSRCM )*ZFLX*PINV_PDYY )*PFRAC_ICE(:,:,:)
    ELSE
      PRRS(:,:,:,2) = PRRS(:,:,:,2) - 2. *                                     &
        DYF( MYM( PRHODJ*PATHETA*PSRCM )*ZFLX*PINV_PDYY )
    END IF
  END IF
END IF
!
!! stores this flux in ZWORK to compute later <V' VPT'>
!!ZWORK(:,:,:) = ZFLX(:,:,:) 
!
! stores the horizontal  <V THl>
IF ( tpfile%lopened .AND. OTURB_FLX ) THEN
  TZFIELD%CMNHNAME   = 'VTHL_FLX'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'VTHL_FLX'
  TZFIELD%CUNITS     = 'K m s-1'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Z_VTHL_FLX'
  TZFIELD%NGRID      = 3
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZFLX)
END IF
!
IF (KSPLT==1 .AND. LLES_CALL) THEN
  CALL SECOND_MNH(ZTIME1)
  CALL LES_MEAN_SUBGRID( MYF(ZFLX), X_LES_SUBGRID_VThl ) 
  CALL LES_MEAN_SUBGRID( MZF(MYF(GY_W_VW(PWM,PDYY,PDZZ,PDZY)*MZM(ZFLX))),&
                         X_LES_RES_ddxa_W_SBG_UaThl , .TRUE. )
  CALL LES_MEAN_SUBGRID( GY_M_M(PTHLM,PDYY,PDZZ,PDZY)*MYF(ZFLX),&
                         X_LES_RES_ddxa_Thl_SBG_UaThl , .TRUE. )
  IF (KRR>=1) THEN
    CALL LES_MEAN_SUBGRID( GY_M_M(PRM(:,:,:,1),PDYY,PDZZ,PDZY)*MYF(ZFLX),&
                           X_LES_RES_ddxa_Rt_SBG_UaThl , .TRUE. )
  END IF
  CALL SECOND_MNH(ZTIME2)
  XTIME_LES = XTIME_LES + ZTIME2 - ZTIME1
END IF
!
!
!*       6.   < V' R'np >
!             -----------
!
IF (KRR/=0) THEN
  !
  IF (.NOT. L2D) THEN
    ZFLX(:,:,:)     = -XCHF * MYM( PK ) * GY_M_V(1,IKU,1,PRM(:,:,:,1),PDYY,PDZZ,PDZY)
    ZFLX(:,:,IKE+1) = ZFLX(:,:,IKE) 
  ELSE
    ZFLX(:,:,:)     = 0.
  END IF
!
! Compute the flux at the first inner U-point with an uncentred vertical  
! gradient
  ZFLX(:,:,IKB:IKB) = -XCHF * MYM( PK(:,:,IKB:IKB) ) *           &
    ( DYM(PRM(:,:,IKB:IKB,1)) * PINV_PDYY(:,:,IKB:IKB)           &
     -MYM( ZCOEFF(:,:,IKB+2:IKB+2)*PRM(:,:,IKB+2:IKB+2,1)        &
           +ZCOEFF(:,:,IKB+1:IKB+1)*PRM(:,:,IKB+1:IKB+1,1)       &
           +ZCOEFF(:,:,IKB  :IKB  )*PRM(:,:,IKB  :IKB  ,1) )     &
           *0.5* ( PDZY(:,:,IKB+1:IKB+1)+PDZY(:,:,IKB:IKB))      &
          * PINV_PDYY(:,:,IKB:IKB) )
! extrapolates the flux under the ground so that the vertical average with 
! the IKB flux gives the ground value  ( warning the tangential surface
! flux has been set to 0 for the moment !!  to be improved )
  ZFLX(:,:,IKB-1:IKB-1) = 2. * MYM(  SPREAD( PSFRM(:,:)* PDIRCOSYW(:,:), 3,1) ) &
                       - ZFLX(:,:,IKB:IKB)
  !
  ! Add this source to the conservative mixing ratio sources
  !
  IF (.NOT. L2D) THEN 
    IF (.NOT. LFLAT) THEN
      PRRS(:,:,:,1) = PRRS(:,:,:,1)                                              &
                    - DYF( MYM(PRHODJ) * ZFLX * PINV_PDYY )                           &

                    + DZF( PMZM_PRHODJ *MYF(PDZY*(MZM(ZFLX * PINV_PDYY))) * PINV_PDZZ )
    ELSE
      PRRS(:,:,:,1) = PRRS(:,:,:,1) - DYF( MYM(PRHODJ) * ZFLX * PINV_PDYY )
    END IF
  END IF
  !
  ! Compute the equivalent tendancy for Rc and Ri
  !
  IF ( KRRL >= 1 .AND. .NOT. L2D) THEN   ! Sub-grid condensation
    IF (.NOT. LFLAT) THEN
      ZFLXC = ZFLXC            &
            + 2.*( MXF( MYM( PRHODJ*PAMOIST*PSRCM )*ZFLX )                     &
                +  MZF( MZM( PRHODJ*PAMOIST*PSRCM )*MYF(                       &
                                               PDZY*(MZM( ZFLX*PINV_PDYY )) ) )&
                 )
      IF ( KRRI >= 1 ) THEN
        PRRS(:,:,:,2) = PRRS(:,:,:,2) +  2. *                                  &
        (- DYF( MYM( PRHODJ*PAMOIST*PSRCM )*ZFLX/PDYY )                        &
         + DZF( MZM( PRHODJ*PAMOIST*PSRCM )*MYF( PDZY*(MZM( ZFLX*PINV_PDYY )) )&
                                           * PINV_PDZZ )                       &
        )*(1.0-PFRAC_ICE(:,:,:))
        PRRS(:,:,:,4) = PRRS(:,:,:,4) +  2. *                                  &
        (- DYF( MYM( PRHODJ*PAMOIST*PSRCM )*ZFLX/PDYY )                        &
         + DZF( MZM( PRHODJ*PAMOIST*PSRCM )*MYF( PDZY*(MZM( ZFLX*PINV_PDYY )) )&
                                           * PINV_PDZZ )                       &
        )*PFRAC_ICE(:,:,:)
      ELSE
        PRRS(:,:,:,2) = PRRS(:,:,:,2) +  2. *                                  &
        (- DYF( MYM( PRHODJ*PAMOIST*PSRCM )*ZFLX/PDYY )                        &
         + DZF( MZM( PRHODJ*PAMOIST*PSRCM )*MYF( PDZY*(MZM( ZFLX*PINV_PDYY )) )&
                                           * PINV_PDZZ )                       &
        )
      END IF
    ELSE
      ZFLXC = ZFLXC + 2.*MXF( MYM( PRHODJ*PAMOIST*PSRCM )*ZFLX )
      IF ( KRRI >= 1 ) THEN
        PRRS(:,:,:,2) = PRRS(:,:,:,2) - 2. *                                   &
        DYF( MYM( PRHODJ*PAMOIST*PSRCM )*ZFLX/PDYY )*(1.0-PFRAC_ICE(:,:,:))
        PRRS(:,:,:,4) = PRRS(:,:,:,4) - 2. *                                   &
        DYF( MYM( PRHODJ*PAMOIST*PSRCM )*ZFLX/PDYY )*PFRAC_ICE(:,:,:)
      ELSE
        PRRS(:,:,:,2) = PRRS(:,:,:,2) - 2. *                                   &
        DYF( MYM( PRHODJ*PAMOIST*PSRCM )*ZFLX/PDYY )
      END IF
    END IF
  END IF
  !
  ! stores the horizontal  <V Rnp>
  IF ( tpfile%lopened .AND. OTURB_FLX ) THEN
    TZFIELD%CMNHNAME   = 'VR_FLX'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'VR_FLX'
    TZFIELD%CUNITS     = 'kg kg-1 m s-1'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'X_Y_Z_VR_FLX'
    TZFIELD%NGRID      = 3
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,ZFLX)
  END IF
  !
  IF (KSPLT==1 .AND. LLES_CALL) THEN
    CALL SECOND_MNH(ZTIME1)
    CALL LES_MEAN_SUBGRID( MYF(ZFLX), X_LES_SUBGRID_VRt ) 
    CALL LES_MEAN_SUBGRID( MZF(MYF(GY_W_VW(PWM,PDYY,PDZZ,PDZY)*MZM(ZFLX))),&
                           X_LES_RES_ddxa_W_SBG_UaRt , .TRUE. )
    CALL LES_MEAN_SUBGRID( GY_M_M(PTHLM,PDYY,PDZZ,PDZY)*MYF(ZFLX), &
                           X_LES_RES_ddxa_Thl_SBG_UaRt , .TRUE. )
    CALL LES_MEAN_SUBGRID( GY_M_M(PRM(:,:,:,1),PDYY,PDZZ,PDZY)*MYF(ZFLX), &
                           X_LES_RES_ddxa_Rt_SBG_UaRt , .TRUE. )
    CALL SECOND_MNH(ZTIME2)
    XTIME_LES = XTIME_LES + ZTIME2 - ZTIME1
  END IF
!
  !
  IF (KRRL>0 .AND. KSPLT==1 .AND. LLES_CALL) THEN
    CALL SECOND_MNH(ZTIME1)
    CALL LES_MEAN_SUBGRID(MYF(ZFLXC), X_LES_SUBGRID_VRc )
    CALL SECOND_MNH(ZTIME2)
    XTIME_LES = XTIME_LES + ZTIME2 - ZTIME1
  END IF
  !
END IF
!
!*       7.   < V' TPV' >
!             -----------
!
!! to be tested later
!!IF (KRR/=0) THEN
!!  ! here ZFLX= <V'R'np> and ZWORK= <V'Theta'l>
!!  !
!!  IF (.NOT. L2D) THEN        &
!!    ZVPTV(:,:,:) =                                                         &
!!        ZWORK(:,:,:)*MYM(ETHETA(KRR,KRRI,PTHLT,PEXNREF,PRT,PLOCPT,PSRCM)) +       &
!!         ZFLX(:,:,:)*MYM(EMOIST(KRR,KRRI,PTHLT,PEXNREF,PRT,PLOCPT,PSRCM))
!!  ELSE
!!    ZVPTV(:,:,:) = 0.
!!  END IF
!!  !
!!  ! stores the horizontal  <V VPT>
!!  IF ( tpfile%lopened .AND. OTURB_FLX ) THEN
!!    TZFIELD%CMNHNAME   = 'VVPT_FLX'
!!    TZFIELD%CSTDNAME   = ''
!!    TZFIELD%CLONGNAME  = 'VVPT_FLX'
!!    TZFIELD%CUNITS     = 'K m s-1'
!!    TZFIELD%CDIR       = 'XY'
!!    TZFIELD%CCOMMENT   = 'X_Y_Z_VVPT_FLX'
!!    TZFIELD%NGRID      = 3
!!    TZFIELD%NTYPE      = TYPEREAL
!!    TZFIELD%NDIMS      = 3
!!    TZFIELD%LTIMEDEP   = .TRUE.
!!    CALL IO_Field_write(TPFILE,TZFIELD,ZVPTV)
!!  END IF
!!!
!!ELSE
!!  ZVPTV(:,:,:)=ZWORK(:,:,:)
!!END IF
!
!
END SUBROUTINE TURB_HOR_THERMO_FLUX
