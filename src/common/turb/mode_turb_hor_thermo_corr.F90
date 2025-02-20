!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_TURB_HOR_THERMO_CORR
IMPLICIT NONE
CONTAINS
      SUBROUTINE TURB_HOR_THERMO_CORR(D,CST,TURBN,NEBN,TLES,         &
                      KRR, KRRL, KRRI,                               &
                      OOCEAN,OCOMPUTE_SRC,O2D, OFLAT,                &
                      TPFILE,                                        &
                      PINV_PDXX,PINV_PDYY,                           &
                      PDXX,PDYY,PDZZ,PDZX,PDZY,                      &
                      PTHVREF,                                       &
                      PWM,PTHLM,PRM,                                 &
                      PTKEM,PLM,PLEPS,                               &
                      PLOCPEXNM,PATHETA,PAMOIST,PSRCM,               &
                      PSIGS                                          )
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
!!                     Nov  06, 2002 (V. Masson) LES budgets
!!                     Feb  20, 2003 (JP Pinty) Add PFRAC_ICE
!!                     October 2009 (G. Tanguy) add ILENCH=LEN(YCOMMENT) after
!!                                              change of YCOMMENT
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!! --------------------------------------------------------------------------
!
!*      0. DECLARATIONS
!          ------------
!
USE MODD_CST,           ONLY : CST_t
USE MODD_CTURB
USE MODD_DIMPHYEX,       ONLY: DIMPHYEX_t
USE MODD_FIELD,          ONLY: TFIELDMETADATA, TYPEREAL
USE MODD_NEB_n,          ONLY: NEB_t
USE MODD_IO,             ONLY: TFILEDATA
USE MODD_LES, ONLY: TLES_t
USE MODD_PARAMETERS
USE MODD_TURB_n,         ONLY: TURB_t
!
USE MODE_IO_FIELD_WRITE, ONLY: IO_FIELD_WRITE
!
USE MODI_GRADIENT_M
USE MODI_GRADIENT_U
USE MODI_GRADIENT_V
USE MODI_GRADIENT_W
USE MODI_SHUMAN 
USE MODE_SHUMAN_PHY
USE MODE_GRADIENT_M_PHY, ONLY : GX_M_U_PHY, GY_M_V_PHY, GY_M_M_PHY, GX_M_M_PHY
USE MODE_GRADIENT_W_PHY, ONLY : GX_W_UW_PHY, GY_W_VW_PHY
USE MODI_LES_MEAN_SUBGRID
!
USE MODE_EMOIST, ONLY: EMOIST
USE MODE_ETHETA, ONLY: ETHETA
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
TYPE(DIMPHYEX_t),       INTENT(IN)   :: D
TYPE(CST_t),            INTENT(IN)   :: CST
TYPE(TURB_t),           INTENT(IN)   :: TURBN
TYPE(NEB_t),            INTENT(IN)   :: NEBN
TYPE(TLES_t),           INTENT(INOUT):: TLES          ! modd_les structure
INTEGER,                INTENT(IN)   :: KRR           ! number of moist var.
INTEGER,                INTENT(IN)   :: KRRL          ! number of liquid water var.
INTEGER,                INTENT(IN)   :: KRRI          ! number of ice water var.
LOGICAL,                INTENT(IN)   ::  OOCEAN       ! switch for Ocean model version
LOGICAL,                INTENT(IN)   ::  OCOMPUTE_SRC ! flag to define dimensions of SIGS and SRCT variables
LOGICAL,                INTENT(IN)   ::  O2D          ! Logical for 2D model version (modd_conf)
LOGICAL,                INTENT(IN)   ::  OFLAT        ! Logical for zero ororography
TYPE(TFILEDATA),          INTENT(INOUT)    ::  TPFILE       ! Output file
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PINV_PDXX   ! 1./PDXX
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PINV_PDYY   ! 1./PDYY
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PDXX, PDYY, PDZZ, PDZX, PDZY 
                                                         ! Metric coefficients
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PTHVREF      ! ref. state Virtual 
                                                      ! Potential Temperature
!
! Variables at t-1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PWM 
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PTHLM 
REAL, DIMENSION(D%NIT,D%NJT,D%NKT,KRR), INTENT(IN)    ::  PRM          ! mixing ratios at t-1,
                              !  where PRM(:,:,:,1) = conservative mixing ratio
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PTKEM        ! Turb. Kin. Energy
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PLM          ! Turb. mixing length
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PLEPS        ! dissipative length
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PLOCPEXNM    ! Lv(T)/Cp/Exnref at time t-1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PATHETA      ! coefficients between 
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PAMOIST      ! s and Thetal and Rnp
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PSRCM        ! normalized 
!
!
!
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(INOUT) ::  PSIGS
                                  ! IN: Vertical part of Sigma_s at t
                                  ! OUT: Total Sigma_s at t
!
!*       0.2  declaration of local variables
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT)   :: ZFLX,ZWORK,ZA ! work arrays
!   
INTEGER             :: IKB,IKE,IIT,IJT,IKT
INTEGER             :: JI,JJ,JK
                                    ! Index values for the Beginning and End
                                    ! mass points of the domain  
REAL, DIMENSION(D%NIT,D%NJT,1+JPVEXT:3+JPVEXT) :: ZCOEFF 
                                    ! coefficients for the uncentred gradient 
                                    ! computation near the ground
REAL :: ZTIME1, ZTIME2
TYPE(TFIELDMETADATA) :: TZFIELD
!
! ---------------------------------------------------------------------------
!
!*       1.   PRELIMINARY COMPUTATIONS
!             ------------------------
!
IKB = 1+JPVEXT               
IKE = SIZE(PTHLM,3)-JPVEXT
IIT=D%NIT
IJT=D%NJT
IKT=D%NKT
!
!
!
!  compute the coefficients for the uncentred gradient computation near the 
!  ground
!$acc kernels present_cr(ZCOEFF)
ZCOEFF(:,:,IKB+2)= - PDZZ(:,:,IKB+1) /      &
       ( (PDZZ(:,:,IKB+2)+PDZZ(:,:,IKB+1)) * PDZZ(:,:,IKB+2) )
ZCOEFF(:,:,IKB+1)=   (PDZZ(:,:,IKB+2)+PDZZ(:,:,IKB+1)) /      &
       ( PDZZ(:,:,IKB+1) * PDZZ(:,:,IKB+2) )
ZCOEFF(:,:,IKB)= - (PDZZ(:,:,IKB+2)+2.*PDZZ(:,:,IKB+1)) /      &
       ( (PDZZ(:,:,IKB+2)+PDZZ(:,:,IKB+1)) * PDZZ(:,:,IKB+1) )
!$acc end kernels
!
!
!*       8.   TURBULENT CORRELATIONS : <THl THl>, <THl Rnp>, <Rnp Rnp>, Sigma_s
!             -----------------------------------------------------------------
!
!
!
IF ( ( KRRL > 0 .AND. NEBN%LSUBG_COND) .OR. ( TURBN%LTURB_FLX .AND. TPFILE%LOPENED ) &
                                  .OR. ( TLES%LLES_CALL )                  ) THEN
!
!*       8.1  <THl THl>
!
  ! Computes the horizontal variance <THl THl>
  IF (.NOT. O2D) THEN
    ZFLX(:,:,:) = TURBN%XCTV * PLM(:,:,:) * PLEPS(:,:,:) *                           &
       ( GX_M_M(OFLAT,PTHLM,PDXX,PDZZ,PDZX)**2 + GY_M_M(OFLAT,PTHLM,PDYY,PDZZ,PDZY)**2 )
  ELSE
    ZFLX(:,:,:) = TURBN%XCTV * PLM(:,:,:) * PLEPS(:,:,:) *                           &
         GX_M_M(OFLAT,PTHLM,PDXX,PDZZ,PDZX)**2
  END IF
!
! Compute the flux at the first inner U-point with an uncentred vertical  
! gradient
!
  ZFLX(:,:,IKB) = TURBN%XCTV * PLM(:,:,IKB)                  &
  * PLEPS(:,:,IKB) *  (                                    &
  ( MXF(DXM(PTHLM(:,:,IKB)) * PINV_PDXX(:,:,IKB))      &
   - ( ZCOEFF(:,:,IKB+2)*PTHLM(:,:,IKB+2)          &
      +ZCOEFF(:,:,IKB+1)*PTHLM(:,:,IKB+1)          &
      +ZCOEFF(:,:,IKB)*PTHLM(:,:,IKB)          &
     ) * 0.5 * ( PDZX(:,:,IKB+1)+PDZX(:,:,IKB) )     &
     / MXF(PDXX(:,:,IKB))                                  &
  ) ** 2 +                                                     &
  ( MYF(DYM(PTHLM(:,:,IKB)) * PINV_PDYY(:,:,IKB))      &
   - ( ZCOEFF(:,:,IKB+2)*PTHLM(:,:,IKB+2)          &
      +ZCOEFF(:,:,IKB+1)*PTHLM(:,:,IKB+1)          &
      +ZCOEFF(:,:,IKB)*PTHLM(:,:,IKB)          &
     ) * 0.5 * ( PDZY(:,:,IKB+1)+PDZY(:,:,IKB) )     &
     / MYF(PDYY(:,:,IKB))                                  &
  ) ** 2                                             )
  !
!$acc kernels present_cr(ZFLX)
  !$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
  ZFLX(:,:,IKB-1) = ZFLX(:,:,IKB)
  !$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)
  !
  IF ( KRRL > 0 ) THEN
!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
    ZWORK(:,:,:) = ZFLX(:,:,:) * PATHETA(:,:,:) * PATHETA(:,:,:)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
  END IF
!$acc end kernels
  !
  ! stores <THl THl>
  IF ( TURBN%LTURB_FLX .AND. TPFILE%LOPENED ) THEN
    TZFIELD = TFIELDMETADATA(        &
      CMNHNAME   = 'THL_HVAR',       &
      CSTDNAME   = '',               &
      CLONGNAME  = 'THL_HVAR',       &
      CUNITS     = 'K2',             &
      CDIR       = 'XY',             &
      CCOMMENT   = 'X_Y_Z_THL_HVAR', &
      NGRID      = 1,                &
      NTYPE      = TYPEREAL,         &
      NDIMS      = 3,                &
      LTIMEDEP   = .TRUE.            )
!$acc update self(ZFLX)
    CALL IO_FIELD_WRITE(TPFILE,TZFIELD,ZFLX)
  END IF
!
! Storage in the LES configuration (addition to TURB_VER computation)
!
  IF (TLES%LLES_CALL) THEN
    CALL SECOND_MNH(ZTIME1)
    CALL LES_MEAN_SUBGRID( ZFLX, TLES%X_LES_SUBGRID_Thl2, .TRUE. ) 
    CALL LES_MEAN_SUBGRID( MZF(PWM)*ZFLX, TLES%X_LES_RES_W_SBG_Thl2, .TRUE. )
    CALL LES_MEAN_SUBGRID( -2.*XCTD*SQRT(PTKEM)*ZFLX/PLEPS ,TLES%X_LES_SUBGRID_DISS_Thl2, .TRUE. )
    CALL ETHETA(D,CST,KRR,KRRI,PTHLM,PRM,PLOCPEXNM,PATHETA,PSRCM,OOCEAN,OCOMPUTE_SRC,ZA)
    CALL LES_MEAN_SUBGRID( ZA*ZFLX, TLES%X_LES_SUBGRID_ThlThv, .TRUE. ) 
    CALL LES_MEAN_SUBGRID( -CST%XG/PTHVREF/3.*ZA*ZFLX, TLES%X_LES_SUBGRID_ThlPz, .TRUE. )
    CALL SECOND_MNH(ZTIME2)
    TLES%XTIME_LES = TLES%XTIME_LES + ZTIME2 - ZTIME1
  END IF
!
  IF ( KRR /= 0 ) THEN
!
!*       8.3  <THl Rnp>
!
    ! Computes the horizontal correlation <THl Rnp>
    IF (.NOT. O2D) THEN
      ZFLX(:,:,:)=                                                               &
            PLM(:,:,:) * PLEPS(:,:,:) *                                          &
            (GX_M_M(OFLAT,PTHLM,PDXX,PDZZ,PDZX) * GX_M_M(OFLAT,PRM(:,:,:,1),PDXX,PDZZ,PDZX)  &
           + GY_M_M(OFLAT,PTHLM,PDYY,PDZZ,PDZY) * GY_M_M(OFLAT,PRM(:,:,:,1),PDYY,PDZZ,PDZY)  &
            ) * (TURBN%XCHT1+TURBN%XCHT2)
    ELSE
      ZFLX(:,:,:)=                                                               &
            PLM(:,:,:) * PLEPS(:,:,:) *                                          &
            (GX_M_M(OFLAT,PTHLM,PDXX,PDZZ,PDZX) * GX_M_M(OFLAT,PRM(:,:,:,1),PDXX,PDZZ,PDZX)  &
            ) * (TURBN%XCHT1+TURBN%XCHT2)

    END IF
!
! Compute the flux at the first inner U-point with an uncentred vertical  
! gradient
    ZFLX(:,:,IKB) = (TURBN%XCHT1+TURBN%XCHT2) * PLM(:,:,IKB)         &
    * PLEPS(:,:,IKB)  *  (                                   &
    ( MXF(DXM(PTHLM(:,:,IKB)) * PINV_PDXX(:,:,IKB))      &
     - ( ZCOEFF(:,:,IKB+2)*PTHLM(:,:,IKB+2)          &
        +ZCOEFF(:,:,IKB+1)*PTHLM(:,:,IKB+1)          &
        +ZCOEFF(:,:,IKB)*PTHLM(:,:,IKB)          &
       ) * 0.5 * ( PDZX(:,:,IKB+1)+PDZX(:,:,IKB) )     &
       / MXF(PDXX(:,:,IKB))                                  &
    ) *                                                          &
    ( MXF(DXM(PRM(:,:,IKB,1)) * PINV_PDXX(:,:,IKB))      &
     - ( ZCOEFF(:,:,IKB+2)*PRM(:,:,IKB+2,1)          &
        +ZCOEFF(:,:,IKB+1)*PRM(:,:,IKB+1,1)          &
        +ZCOEFF(:,:,IKB)*PRM(:,:,IKB,1)          &
       ) * 0.5 * ( PDZX(:,:,IKB+1)+PDZX(:,:,IKB) )     &
       / MXF(PDXX(:,:,IKB))                                  &
    ) +                                                          &
    ( MYF(DYM(PTHLM(:,:,IKB)) * PINV_PDYY(:,:,IKB))      &
     - ( ZCOEFF(:,:,IKB+2)*PTHLM(:,:,IKB+2)          &
        +ZCOEFF(:,:,IKB+1)*PTHLM(:,:,IKB+1)          &
        +ZCOEFF(:,:,IKB)*PTHLM(:,:,IKB)          &
       ) * 0.5 * ( PDZY(:,:,IKB+1)+PDZY(:,:,IKB) )     &
       / MYF(PDYY(:,:,IKB))                                  &
    ) *                                                          &
    ( MYF(DYM(PRM(:,:,IKB,1)) * PINV_PDYY(:,:,IKB))           &
     - ( ZCOEFF(:,:,IKB+2)*PRM(:,:,IKB+2,1)          &
        +ZCOEFF(:,:,IKB+1)*PRM(:,:,IKB+1,1)          &
        +ZCOEFF(:,:,IKB)*PRM(:,:,IKB,1)          &
       ) * 0.5 * ( PDZY(:,:,IKB+1)+PDZY(:,:,IKB) )     &
       / MYF(PDYY(:,:,IKB))                                  &
    )                                                          )
    !
!$acc kernels present_cr(ZFLX)
    !$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
    ZFLX(:,:,IKB-1) = ZFLX(:,:,IKB)
    !$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)
    !
    IF ( KRRL > 0 )  THEN
!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
      ZWORK(:,:,:) = ZWORK(:,:,:) +       &
                     2. * PATHETA(:,:,:) * PAMOIST(:,:,:) * ZFLX(:,:,:)    
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
    END IF
!$acc end kernels
    ! stores <THl Rnp>
    IF ( TURBN%LTURB_FLX .AND. TPFILE%LOPENED ) THEN
      TZFIELD = TFIELDMETADATA(         &
        CMNHNAME   = 'THLR_HCOR',       &
        CSTDNAME   = '',                &
        CLONGNAME  = 'THLR_HCOR',       &
        CUNITS     = 'K kg kg-1',       &
        CDIR       = 'XY',              &
        CCOMMENT   = 'X_Y_Z_THLR_HCOR', &
        NGRID      = 1,                 &
        NTYPE      = TYPEREAL,          &
        NDIMS      = 3,                 &
        LTIMEDEP   = .TRUE.             )
      CALL IO_FIELD_WRITE(TPFILE,TZFIELD,ZFLX)
    END IF
!
!   Storage in the LES configuration (addition to TURB_VER computation)
!
    IF (TLES%LLES_CALL) THEN
      CALL SECOND_MNH(ZTIME1)
      CALL LES_MEAN_SUBGRID( ZFLX, TLES%X_LES_SUBGRID_ThlRt, .TRUE. ) 
      CALL LES_MEAN_SUBGRID( MZF(PWM)*ZFLX, TLES%X_LES_RES_W_SBG_ThlRt, .TRUE. )
      CALL LES_MEAN_SUBGRID( -XCTD*SQRT(PTKEM)*ZFLX/PLEPS ,TLES%X_LES_SUBGRID_DISS_ThlRt, .TRUE. )
      CALL LES_MEAN_SUBGRID( ZA*ZFLX, TLES%X_LES_SUBGRID_RtThv, .TRUE. ) 
      CALL LES_MEAN_SUBGRID( -CST%XG/PTHVREF/3.*ZA*ZFLX, TLES%X_LES_SUBGRID_RtPz,.TRUE.)
      CALL EMOIST(D,CST,KRR,KRRI,PTHLM,PRM,PLOCPEXNM,PAMOIST,PSRCM,OOCEAN,ZA)
      CALL LES_MEAN_SUBGRID( ZA*ZFLX, TLES%X_LES_SUBGRID_ThlThv, .TRUE. ) 
      CALL LES_MEAN_SUBGRID( -CST%XG/PTHVREF/3.*ZA*ZFLX, TLES%X_LES_SUBGRID_ThlPz,.TRUE.)
      CALL SECOND_MNH(ZTIME2)
      TLES%XTIME_LES = TLES%XTIME_LES + ZTIME2 - ZTIME1
    END IF
!
!*       8.4  <Rnp Rnp>
!
    ! Computes the horizontal variance <Rnp Rnp>
    IF (.NOT. O2D) THEN
      ZFLX(:,:,:) = TURBN%XCHV * PLM(:,:,:) * PLEPS(:,:,:) *                      &
           ( GX_M_M(OFLAT,PRM(:,:,:,1),PDXX,PDZZ,PDZX)**2 +                       &
             GY_M_M(OFLAT,PRM(:,:,:,1),PDYY,PDZZ,PDZY)**2 )
    ELSE
      ZFLX(:,:,:) = TURBN%XCHV * PLM(:,:,:) * PLEPS(:,:,:) *                      &
           ( GX_M_M(OFLAT,PRM(:,:,:,1),PDXX,PDZZ,PDZX)**2  )
    END IF
!
! Compute the flux at the first inner U-point with an uncentred vertical  
! gradient
    ZFLX(:,:,IKB) = TURBN%XCHV * PLM(:,:,IKB)                  &
    * PLEPS(:,:,IKB) *  (                                    &
    ( MXF(DXM(PRM(:,:,IKB,1)) * PINV_PDXX(:,:,IKB))      &
     - ( ZCOEFF(:,:,IKB+2)*PRM(:,:,IKB+2,1)          &
        +ZCOEFF(:,:,IKB+1)*PRM(:,:,IKB+1,1)          &
        +ZCOEFF(:,:,IKB)*PRM(:,:,IKB,1)          &
       ) * 0.5 * ( PDZX(:,:,IKB+1)+PDZX(:,:,IKB) )     &
       / MXF(PDXX(:,:,IKB:))                                  &
    ) ** 2 +                                                     &
    ( MYF(DYM(PRM(:,:,IKB,1)) * PINV_PDYY(:,:,IKB))           &
     - ( ZCOEFF(:,:,IKB+2)*PRM(:,:,IKB+2,1)          &
        +ZCOEFF(:,:,IKB+1)*PRM(:,:,IKB+1,1)          &
        +ZCOEFF(:,:,IKB)*PRM(:,:,IKB,1)          &
       ) * 0.5 * ( PDZY(:,:,IKB+1)+PDZY(:,:,IKB) )     &
       / MYF(PDYY(:,:,IKB))                                  &
    ) ** 2                                             )
!
!$acc kernels present_cr(ZFLX)
    !$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
    ZFLX(:,:,IKB-1) = ZFLX(:,:,IKB)
    !$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)
    !
    IF ( KRRL > 0 ) THEN       
    !$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
      ZWORK(:,:,:) = ZWORK(:,:,:)+ PAMOIST(:,:,:) * PAMOIST(:,:,:) * ZFLX(:,:,:)
    !$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
    END IF
!$acc end kernels
    ! stores <Rnp Rnp>
    IF ( TURBN%LTURB_FLX .AND. TPFILE%LOPENED ) THEN
      TZFIELD = TFIELDMETADATA(      &
        CMNHNAME   = 'R_HVAR',       &
        CSTDNAME   = '',             &
        CLONGNAME  = 'R_HVAR',       &
        CUNITS     = 'kg2 kg-2',     &
        CDIR       = 'XY',           &
        CCOMMENT   = 'X_Y_Z_R_HVAR', &
        NGRID      = 1,              &
        NTYPE      = TYPEREAL,       &
        NDIMS      = 3,              &
        LTIMEDEP   = .TRUE.          )
!$acc update self(ZFLX)
      CALL IO_FIELD_WRITE(TPFILE,TZFIELD,ZFLX)
    END IF
    !
    !   Storage in the LES configuration (addition to TURB_VER computation)
    !
    IF (TLES%LLES_CALL) THEN
      CALL SECOND_MNH(ZTIME1)
      CALL LES_MEAN_SUBGRID( ZFLX, TLES%X_LES_SUBGRID_Rt2, .TRUE. ) 
      CALL LES_MEAN_SUBGRID( MZF(PWM)*ZFLX, TLES%X_LES_RES_W_SBG_Rt2, .TRUE. )
      CALL LES_MEAN_SUBGRID( ZA*ZFLX, TLES%X_LES_SUBGRID_RtThv, .TRUE. ) 
      CALL LES_MEAN_SUBGRID( -CST%XG/PTHVREF/3.*ZA*ZFLX, TLES%X_LES_SUBGRID_RtPz,.TRUE.)
      CALL LES_MEAN_SUBGRID( -2.*XCTD*SQRT(PTKEM)*ZFLX/PLEPS, TLES%X_LES_SUBGRID_DISS_Rt2, .TRUE. )
      CALL SECOND_MNH(ZTIME2)
      TLES%XTIME_LES = TLES%XTIME_LES + ZTIME2 - ZTIME1
    END IF
  !
  END IF
!
!        8.5 Complete the Sigma_s computation:
!
  IF ( KRRL > 0 ) THEN   
    !
    !$acc kernels present_cr(PSIGS)
    PSIGS(:,:,:)=PSIGS(:,:,:)*PSIGS(:,:,:) + ZWORK(:,:,:)
    ! Extrapolate PSIGS at the ground and at the top
    PSIGS(:,:,IKB-1) = PSIGS(:,:,IKB)
    PSIGS(:,:,IKE+1) = PSIGS(:,:,IKE)
    PSIGS(:,:,:) = SQRT(MAX ( PSIGS(:,:,:),1.E-12) ) 
    !$acc end kernels
  END IF       
!
END IF
!
!
!
END SUBROUTINE TURB_HOR_THERMO_CORR
END MODULE MODE_TURB_HOR_THERMO_CORR
