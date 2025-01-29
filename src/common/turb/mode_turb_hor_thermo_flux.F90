!MNH_LIC Copyright 1994-2024 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_TURB_HOR_THERMO_FLUX
IMPLICIT NONE
CONTAINS
!     ################################################################
      SUBROUTINE TURB_HOR_THERMO_FLUX(D,TURBN, TLES,KSPLT, KRR, KRRL, KRRI, &
                      TPFILE,OFLAT,O2D,                              &
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
!!                     Oct  18, 2000 (V. Masson) LES computations + OFLAT switch
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
USE MODD_TURB_n, ONLY: TURB_t
USE MODD_DIMPHYEX,   ONLY: DIMPHYEX_t
!
USE MODD_CST
USE MODD_CTURB
USE MODD_FIELD,          ONLY: TFIELDMETADATA, TYPEREAL
USE MODD_IO,             ONLY: TFILEDATA
USE MODD_PARAMETERS
USE MODD_LES, ONLY: TLES_t
!
USE MODE_IO_FIELD_WRITE, ONLY: IO_FIELD_WRITE
!
USE MODI_GRADIENT_M
USE MODI_GRADIENT_U
USE MODI_GRADIENT_V
USE MODI_GRADIENT_W
USE MODE_GRADIENT_M_PHY, ONLY : GX_M_U_PHY, GY_M_V_PHY, GY_M_M_PHY, GX_M_M_PHY
USE MODE_GRADIENT_W_PHY, ONLY: GX_W_UW_PHY, GY_W_VW_PHY
USE MODI_SHUMAN 
USE MODE_SHUMAN_PHY, ONLY: MZM_PHY, MZF_PHY, MXF_PHY, MXM_PHY, MYM_PHY, MYF_PHY, &
                           DYM_PHY, DYF_PHY, DXM_PHY, DXF_PHY, DZF_PHY, &
                           MXM2D_PHY, MYM2D_PHY, DXM2D_PHY, DYM2D_PHY
USE MODI_LES_MEAN_SUBGRID
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
TYPE(DIMPHYEX_t),         INTENT(IN)    :: D
TYPE(TURB_t),             INTENT(IN)    :: TURBN
TYPE(TLES_t),             INTENT(INOUT) :: TLES          ! modd_les structure
INTEGER,                  INTENT(IN)    :: KSPLT         ! split process index
INTEGER,                  INTENT(IN)    :: KRR           ! number of moist var.
INTEGER,                  INTENT(IN)    :: KRRL          ! number of liquid water var.
INTEGER,                  INTENT(IN)    :: KRRI          ! number of ice water var.
LOGICAL,                  INTENT(IN)    ::  OFLAT        ! Logical for zero ororography
LOGICAL,                  INTENT(IN)    ::  O2D          ! Logical for 2D model version (modd_conf)
TYPE(TFILEDATA),          INTENT(INOUT)    ::  TPFILE       ! Output file
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PK          ! Turbulent diffusion doef.
                                                        ! PK = PLM * SQRT(PTKEM)
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PINV_PDXX   ! 1./PDXX
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PINV_PDYY   ! 1./PDYY
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PINV_PDZZ   ! 1./PDZZ
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PMZM_PRHODJ ! MZM(PRHODJ)
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PDXX, PDYY, PDZZ, PDZX, PDZY 
                                                         ! Metric coefficients
REAL, DIMENSION(D%NIT,D%NJT),     INTENT(IN)    ::  PDIRCOSXW, PDIRCOSYW
! Director Cosinus along x, y and z directions at surface w-point
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PRHODJ       ! density * grid volume
!
REAL, DIMENSION(D%NIT,D%NJT),     INTENT(IN)    ::  PSFTHM,PSFRM
!
! Variables at t-1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PWM 
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PTHLM 
REAL, DIMENSION(D%NIT,D%NJT,D%NKT,KRR), INTENT(IN)    ::  PRM          ! mixing ratios at t-1,
                              !  where PRM(:,:,:,1) = conservative mixing ratio
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PATHETA      ! coefficients between 
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PAMOIST      ! s and Thetal and Rnp

REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PSRCM
                                  ! normalized 2nd-order flux
                                  ! s'r'c/2Sigma_s2 at t-1 multiplied by Lambda_3
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PFRAC_ICE    ! ri fraction of rc+ri
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(INOUT) ::  PRTHLS
REAL, DIMENSION(D%NIT,D%NJT,D%NKT,KRR), INTENT(INOUT) ::  PRRS         ! var. at t+1 -split-
!
!
!
!*       0.2  declaration of local variables
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: ZFLX,ZFLXC ! work arrays
!
INTEGER             :: IKB,IKE,IKU, IKT, IIT, IJT, JI,JJ,JK
                                    ! Index values for the Beginning and End
                                    ! mass points of the domain  
REAL, DIMENSION(D%NIT,D%NJT,1+JPVEXT:3+JPVEXT) :: ZCOEFF 
                                    ! coefficients for the uncentred gradient 
                                    ! computation near the ground
!
REAL :: ZTIME1, ZTIME2
TYPE(TFIELDMETADATA) :: TZFIELD
! ---------------------------------------------------------------------------
!
!*       1.   PRELIMINARY COMPUTATIONS
!             ------------------------
!
IKB = 1+JPVEXT               
IKE = SIZE(PTHLM,3)-JPVEXT    
IKU = SIZE(PTHLM,3)
IIT=D%NIT
IJT=D%NJT
IKT=D%NKT
!
!
!  compute the coefficients for the uncentred gradient computation near the 
!  ground
!$acc kernels present_cr(zcoeff)
ZCOEFF(:,:,IKB+2)= - PDZZ(:,:,IKB+1) /      &
       ( (PDZZ(:,:,IKB+2)+PDZZ(:,:,IKB+1)) * PDZZ(:,:,IKB+2) )
ZCOEFF(:,:,IKB+1)=   (PDZZ(:,:,IKB+2)+PDZZ(:,:,IKB+1)) /      &
       ( PDZZ(:,:,IKB+1) * PDZZ(:,:,IKB+2) )
ZCOEFF(:,:,IKB)= - (PDZZ(:,:,IKB+2)+2.*PDZZ(:,:,IKB+1)) /      &
       ( (PDZZ(:,:,IKB+2)+PDZZ(:,:,IKB+1)) * PDZZ(:,:,IKB+1) )
!$acc end kernels
!
!*       2.   < U' THETA'l >
!             --------------
!
!
#define UTHETA UTHETA_ON
!!!!!!! UTHETA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef UTHETA
ZFLX(:,:,:)     = -TURBN%XCSHF * MXM( PK ) * GX_M_U(OFLAT,PTHLM,PDXX,PDZZ,PDZX)
!$acc kernels present_cr(ZFLX)
ZFLX(:,:,IKE+1) = ZFLX(:,:,IKE) 
!$acc end kernels
!
! Compute the flux at the first inner U-point with an uncentred vertical  
! gradient
ZFLX(:,:,IKB) = -TURBN%XCSHF * MXM( PK(:,:,IKB) ) *          &
  ( DXM(PTHLM(:,:,IKB)) * PINV_PDXX(:,:,IKB)           &
   -MXM( ZCOEFF(:,:,IKB+2)*PTHLM(:,:,IKB+2)        &
         +ZCOEFF(:,:,IKB+1)*PTHLM(:,:,IKB+1)       &
         +ZCOEFF(:,:,IKB)*PTHLM(:,:,IKB))      &
        *0.5* ( PDZX(:,:,IKB+1)+PDZX(:,:,IKB))       &
        * PINV_PDXX(:,:,IKB) )
! extrapolates the flux under the ground so that the vertical average with 
! the IKB flux gives the ground value  ( warning the tangential surface
! flux has been set to 0 for the moment !!  to be improved )
ZFLX(:,:,IKB-1) = 2. * MXM(PSFTHM(:,:)* PDIRCOSXW(:,:))  &
                       - ZFLX(:,:,IKB)
!
! Add this source to the Theta_l sources
!
IF (.NOT. OFLAT) THEN
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
  IF (.NOT. OFLAT) THEN
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
IF ( TPFILE%LOPENED .AND. TURBN%LTURB_FLX ) THEN
  TZFIELD = TFIELDMETADATA(        &
    CMNHNAME   = 'UTHL_FLX',       &
    CSTDNAME   = '',               &
    CLONGNAME  = 'UTHL_FLX',       &
    CUNITS     = 'K m s-1',        &
    CDIR       = 'XY',             &
    CCOMMENT   = 'X_Y_Z_UTHL_FLX', &
    NGRID      = 2,                &
    NTYPE      = TYPEREAL,         &
    NDIMS      = 3,                &
    LTIMEDEP   = .TRUE.            )
!$acc update self(ZFLX)
  CALL IO_FIELD_WRITE(TPFILE,TZFIELD,ZFLX)
END IF
!
IF (KSPLT==1 .AND. TLES%LLES_CALL) THEN
  CALL SECOND_MNH(ZTIME1)
  CALL LES_MEAN_SUBGRID(MXF(ZFLX), TLES%X_LES_SUBGRID_UThl ) 
  CALL LES_MEAN_SUBGRID(MZF(MXF(GX_W_UW(OFLAT, PWM,PDXX,PDZZ,PDZX)*MZM(ZFLX))),&
                         TLES%X_LES_RES_ddxa_W_SBG_UaThl , .TRUE. )
  CALL LES_MEAN_SUBGRID(GX_M_M(OFLAT,PTHLM,PDXX,PDZZ,PDZX)*MXF(ZFLX),&
                         TLES%X_LES_RES_ddxa_Thl_SBG_UaThl , .TRUE. )
  IF (KRR>=1) THEN
    CALL LES_MEAN_SUBGRID(GX_M_M(OFLAT,PRM(:,:,:,1),PDXX,PDZZ,PDZX)*MXF(ZFLX), &
                           TLES%X_LES_RES_ddxa_Rt_SBG_UaThl , .TRUE. )
  END IF
  CALL SECOND_MNH(ZTIME2)
  TLES%XTIME_LES = TLES%XTIME_LES + ZTIME2 - ZTIME1
END IF
#endif
!
!*       3.   < U' R'np >
!             -----------
#define URNP URNP_ON
!!!!! URNP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef URNP
IF (KRR/=0) THEN
  !
  ZFLX(:,:,:)     = -TURBN%XCHF * MXM( PK ) * GX_M_U(OFLAT,PRM(:,:,:,1),PDXX,PDZZ,PDZX)
!$acc kernels present_cr(ZFLX)
  ZFLX(:,:,IKE+1) = ZFLX(:,:,IKE)
!$acc end kernels
!
! Compute the flux at the first inner U-point with an uncentred vertical  
! gradient
  ZFLX(:,:,IKB) = -TURBN%XCHF * MXM( PK(:,:,IKB) ) *           &
    ( DXM(PRM(:,:,IKB,1)) * PINV_PDXX(:,:,IKB)           &
     -MXM( ZCOEFF(:,:,IKB+2)*PRM(:,:,IKB+2,1)        &
           +ZCOEFF(:,:,IKB+1)*PRM(:,:,IKB+1,1)       &
           +ZCOEFF(:,:,IKB)*PRM(:,:,IKB,1))      &
          *0.5* ( PDZX(:,:,IKB+1)+PDZX(:,:,IKB))       &
          * PINV_PDXX(:,:,IKB) )
! extrapolates the flux under the ground so that the vertical average with 
! the IKB flux gives the ground value  ( warning the tangential surface
! flux has been set to 0 for the moment !!  to be improved )
  ZFLX(:,:,IKB-1) = 2. * MXM(PSFRM(:,:)* PDIRCOSXW(:,:)) &
                       - ZFLX(:,:,IKB)
  !
  ! Add this source to the conservative mixing ratio sources
  !
  IF (.NOT. OFLAT) THEN
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
    IF (.NOT. OFLAT) THEN
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
  IF ( TPFILE%LOPENED .AND. TURBN%LTURB_FLX ) THEN
    TZFIELD = TFIELDMETADATA(       &
      CMNHNAME   = 'UR_FLX',        &
      CSTDNAME   = '',              &
      CLONGNAME  = 'UR_FLX',        &
      CUNITS     = 'kg kg-1 m s-1', &
      CDIR       = 'XY',            &
      CCOMMENT   = 'X_Y_Z_UR_FLX',  &
      NGRID      = 2,               &
      NTYPE      = TYPEREAL,        &
      NDIMS      = 3,               &
      LTIMEDEP   = .TRUE.           )
!$acc update self(ZFLX)
    CALL IO_FIELD_WRITE(TPFILE,TZFIELD,ZFLX)
  END IF
  !
  IF (KSPLT==1 .AND. TLES%LLES_CALL) THEN
    CALL SECOND_MNH(ZTIME1)
    CALL LES_MEAN_SUBGRID(MXF(ZFLX), TLES%X_LES_SUBGRID_URt ) 
    CALL LES_MEAN_SUBGRID(MZF(MXF(GX_W_UW(OFLAT, PWM,PDXX,PDZZ,PDZX)*MZM(ZFLX))),&
                           TLES%X_LES_RES_ddxa_W_SBG_UaRt , .TRUE. )
    CALL LES_MEAN_SUBGRID(GX_M_M(OFLAT,PTHLM,PDXX,PDZZ,PDZX)*MXF(ZFLX),&
                           TLES%X_LES_RES_ddxa_Thl_SBG_UaRt , .TRUE. )
    CALL LES_MEAN_SUBGRID(GX_M_M(OFLAT,PRM(:,:,:,1),PDXX,PDZZ,PDZX)*MXF(ZFLX),&
                           TLES%X_LES_RES_ddxa_Rt_SBG_UaRt , .TRUE. )
    CALL SECOND_MNH(ZTIME2)
    TLES%XTIME_LES = TLES%XTIME_LES + ZTIME2 - ZTIME1
  END IF
!
  !
  IF (KRRL>0 .AND. KSPLT==1 .AND. TLES%LLES_CALL) THEN
    CALL SECOND_MNH(ZTIME1)
    CALL LES_MEAN_SUBGRID(MXF(ZFLXC), TLES%X_LES_SUBGRID_URc )
    CALL SECOND_MNH(ZTIME2)
    TLES%XTIME_LES = TLES%XTIME_LES + ZTIME2 - ZTIME1
  END IF
!
END IF
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!!  IF ( TPFILE%LOPENED .AND. TURBN%LTURB_FLX ) THEN
!!    TZFIELD = TFIELDMETADATA(        &
!!      CMNHNAME   = 'UVPT_FLX',       &
!!      CSTDNAME   = '',               &
!!      CLONGNAME  = 'UVPT_FLX',       &
!!      CUNITS     = 'K m s-1',        &
!!      CDIR       = 'XY',             &
!!      CCOMMENT   = 'X_Y_Z_UVPT_FLX', &
!!      NGRID      = 2,                &
!!      NTYPE      = TYPEREAL,         &
!!      NDIMS      = 3,                &
!!      LTIMEDEP   = .TRUE.            )
!!    CALL IO_FIELD_WRITE(TPFILE,TZFIELD,ZVPTU)
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
#define VTHETA VTHETA_ON
!!!!!! VTHETA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef VTHETA
IF (.NOT. O2D) THEN
  ZFLX(:,:,:)     = -TURBN%XCSHF * MYM( PK ) * GY_M_V(OFLAT,PTHLM,PDYY,PDZZ,PDZY)
!$acc kernels present_cr(ZFLX)
  ZFLX(:,:,IKE+1) = ZFLX(:,:,IKE)
!$acc end kernels
ELSE
!$acc kernels present_cr(ZFLX)
  ZFLX(:,:,:)     = 0.
!$acc end kernels
END IF
!
!
! Compute the flux at the first inner U-point with an uncentred vertical  
! gradient
ZFLX(:,:,IKB) = -TURBN%XCSHF * MYM( PK(:,:,IKB) ) *          &
  ( DYM(PTHLM(:,:,IKB)) * PINV_PDYY(:,:,IKB)           &
   -MYM( ZCOEFF(:,:,IKB+2)*PTHLM(:,:,IKB+2)        &
         +ZCOEFF(:,:,IKB+1)*PTHLM(:,:,IKB+1)       &
         +ZCOEFF(:,:,IKB)*PTHLM(:,:,IKB) )     &
        *0.5* ( PDZY(:,:,IKB+1)+PDZY(:,:,IKB))       &
        * PINV_PDYY(:,:,IKB) )
! extrapolates the flux under the ground so that the vertical average with 
! the IKB flux gives the ground value  ( warning the tangential surface
! flux has been set to 0 for the moment !!  to be improved )
ZFLX(:,:,IKB-1) = 2. * MYM( PSFTHM(:,:)* PDIRCOSYW(:,:) ) &
                       - ZFLX(:,:,IKB)
!
! Add this source to the Theta_l sources
!
IF (.NOT. O2D) THEN 
  IF (.NOT. OFLAT) THEN
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
!IF ( TURBN%LSUBG_COND .AND. KRRL > 0 .AND. .NOT. O2D) THEN
IF ( KRRL >= 1 .AND. .NOT. O2D) THEN
  IF (.NOT. OFLAT) THEN
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
IF ( TPFILE%LOPENED .AND. TURBN%LTURB_FLX ) THEN
  TZFIELD = TFIELDMETADATA(        &
    CMNHNAME   = 'VTHL_FLX',       &
    CSTDNAME   = '',               &
    CLONGNAME  = 'VTHL_FLX',       &
    CUNITS     = 'K m s-1',        &
    CDIR       = 'XY',             &
    CCOMMENT   = 'X_Y_Z_VTHL_FLX', &
    NGRID      = 3,                &
    NTYPE      = TYPEREAL,         &
    NDIMS      = 3,                &
    LTIMEDEP   = .TRUE.            )
!$acc update self(ZFLX)
  CALL IO_FIELD_WRITE(TPFILE,TZFIELD,ZFLX)
END IF
!
IF (KSPLT==1 .AND. TLES%LLES_CALL) THEN
  CALL SECOND_MNH(ZTIME1)
  CALL LES_MEAN_SUBGRID(MYF(ZFLX), TLES%X_LES_SUBGRID_VThl ) 
  CALL LES_MEAN_SUBGRID(MZF(MYF(GY_W_VW(OFLAT, PWM,PDYY,PDZZ,PDZY)*MZM(ZFLX))),&
                         TLES%X_LES_RES_ddxa_W_SBG_UaThl , .TRUE. )
  CALL LES_MEAN_SUBGRID(GY_M_M(OFLAT,PTHLM,PDYY,PDZZ,PDZY)*MYF(ZFLX),&
                         TLES%X_LES_RES_ddxa_Thl_SBG_UaThl , .TRUE. )
  IF (KRR>=1) THEN
    CALL LES_MEAN_SUBGRID(GY_M_M(OFLAT,PRM(:,:,:,1),PDYY,PDZZ,PDZY)*MYF(ZFLX),&
                           TLES%X_LES_RES_ddxa_Rt_SBG_UaThl , .TRUE. )
  END IF
  CALL SECOND_MNH(ZTIME2)
  TLES%XTIME_LES = TLES%XTIME_LES + ZTIME2 - ZTIME1
END IF
#endif
!
!
!*       6.   < V' R'np >
!             -----------
#define VRNP VRNP_ON
!!!!!!! VRNP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef VRNP
IF (KRR/=0) THEN
  !
  IF (.NOT. O2D) THEN
    ZFLX(:,:,:)     = -TURBN%XCHF * MYM( PK ) * GY_M_V(OFLAT,PRM(:,:,:,1),PDYY,PDZZ,PDZY)
!$acc kernels present_cr(ZFLX)
    ZFLX(:,:,IKE+1) = ZFLX(:,:,IKE)
!$acc end kernels
  ELSE
!$acc kernels present_cr(ZFLX)
    ZFLX(:,:,:)     = 0.
!$acc end kernels
  END IF
!
! Compute the flux at the first inner U-point with an uncentred vertical  
! gradient
  ZFLX(:,:,IKB) = -TURBN%XCHF * MYM( PK(:,:,IKB) ) *           &
    ( DYM(PRM(:,:,IKB,1)) * PINV_PDYY(:,:,IKB)           &
     -MYM( ZCOEFF(:,:,IKB+2)*PRM(:,:,IKB+2,1)        &
           +ZCOEFF(:,:,IKB+1)*PRM(:,:,IKB+1,1)       &
           +ZCOEFF(:,:,IKB)*PRM(:,:,IKB,1) )     &
           *0.5* ( PDZY(:,:,IKB+1)+PDZY(:,:,IKB))      &
          * PINV_PDYY(:,:,IKB) )
! extrapolates the flux under the ground so that the vertical average with 
! the IKB flux gives the ground value  ( warning the tangential surface
! flux has been set to 0 for the moment !!  to be improved )
  ZFLX(:,:,IKB-1) = 2. * MYM(PSFRM(:,:)* PDIRCOSYW(:,:)) &
                       - ZFLX(:,:,IKB)
  !
  ! Add this source to the conservative mixing ratio sources
  !
  IF (.NOT. O2D) THEN 
    IF (.NOT. OFLAT) THEN
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
  IF ( KRRL >= 1 .AND. .NOT. O2D) THEN   ! Sub-grid condensation
    IF (.NOT. OFLAT) THEN
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
  IF ( TPFILE%LOPENED .AND. TURBN%LTURB_FLX ) THEN
    TZFIELD = TFIELDMETADATA(       &
      CMNHNAME   = 'VR_FLX',        &
      CSTDNAME   = '',              &
      CLONGNAME  = 'VR_FLX',        &
      CUNITS     = 'kg kg-1 m s-1', &
      CDIR       = 'XY',            &
      CCOMMENT   = 'X_Y_Z_VR_FLX',  &
      NGRID      = 3,               &
      NTYPE      = TYPEREAL,        &
      NDIMS      = 3,               &
      LTIMEDEP   = .TRUE.           )
!$acc update self(ZFLX)
    CALL IO_FIELD_WRITE(TPFILE,TZFIELD,ZFLX)
  END IF
  !
  IF (KSPLT==1 .AND. TLES%LLES_CALL) THEN
    CALL SECOND_MNH(ZTIME1)
    CALL LES_MEAN_SUBGRID(MYF(ZFLX), TLES%X_LES_SUBGRID_VRt ) 
    CALL LES_MEAN_SUBGRID(MZF(MYF(GY_W_VW(OFLAT, PWM,PDYY,PDZZ,PDZY)*MZM(ZFLX))),&
                           TLES%X_LES_RES_ddxa_W_SBG_UaRt , .TRUE. )
    CALL LES_MEAN_SUBGRID(GY_M_M(OFLAT,PTHLM,PDYY,PDZZ,PDZY)*MYF(ZFLX), &
                           TLES%X_LES_RES_ddxa_Thl_SBG_UaRt , .TRUE. )
    CALL LES_MEAN_SUBGRID(GY_M_M(OFLAT,PRM(:,:,:,1),PDYY,PDZZ,PDZY)*MYF(ZFLX), &
                           TLES%X_LES_RES_ddxa_Rt_SBG_UaRt , .TRUE. )
    CALL SECOND_MNH(ZTIME2)
    TLES%XTIME_LES = TLES%XTIME_LES + ZTIME2 - ZTIME1
  END IF
!
  !
  IF (KRRL>0 .AND. KSPLT==1 .AND. TLES%LLES_CALL) THEN
    CALL SECOND_MNH(ZTIME1)
    CALL LES_MEAN_SUBGRID(MYF(ZFLXC), TLES%X_LES_SUBGRID_VRc )
    CALL SECOND_MNH(ZTIME2)
    TLES%XTIME_LES = TLES%XTIME_LES + ZTIME2 - ZTIME1
  END IF
  !
END IF
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!*       7.   < V' TPV' >
!             -----------
!
!! to be tested later
!!IF (KRR/=0) THEN
!!  ! here ZFLX= <V'R'np> and ZWORK= <V'Theta'l>
!!  !
!!  IF (.NOT. O2D) THEN        &
!!    ZVPTV(:,:,:) =                                                         &
!!        ZWORK(:,:,:)*MYM(ETHETA(KRR,KRRI,PTHLT,PEXNREF,PRT,PLOCPT,PSRCM)) +       &
!!         ZFLX(:,:,:)*MYM(EMOIST(KRR,KRRI,PTHLT,PEXNREF,PRT,PLOCPT,PSRCM))
!!  ELSE
!!    ZVPTV(:,:,:) = 0.
!!  END IF
!!  !
!!  ! stores the horizontal  <V VPT>
!!  IF ( TPFILE%LOPENED .AND. TURBN%LTURB_FLX ) THEN
!!    TZFIELD = TFIELDMETADATA(        &
!!      CMNHNAME   = 'VVPT_FLX',       &
!!      CSTDNAME   = '',               &
!!      CLONGNAME  = 'VVPT_FLX',       &
!!      CUNITS     = 'K m s-1',        &
!!      CDIR       = 'XY',             &
!!      CCOMMENT   = 'X_Y_Z_VVPT_FLX', &
!!      NGRID      = 3,                &
!!      NTYPE      = TYPEREAL,         &
!!      NDIMS      = 3,                &
!!      LTIMEDEP   = .TRUE.            )
!!    CALL IO_FIELD_WRITE(TPFILE,TZFIELD,ZVPTV)
!!  END IF
!!!
!!ELSE
!!  ZVPTV(:,:,:)=ZWORK(:,:,:)
!!END IF
!
!
END SUBROUTINE TURB_HOR_THERMO_FLUX
END MODULE MODE_TURB_HOR_THERMO_FLUX
