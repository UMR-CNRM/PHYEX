!MNH_LIC Copyright 1994-2025 CNRS, Meteo-France and Universite Paul Sabatier
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
USE MODD_CTURB, ONLY: XCTD
USE MODD_DIMPHYEX,       ONLY: DIMPHYEX_t
USE MODD_FIELD,          ONLY: TFIELDMETADATA, TYPEREAL
USE MODD_NEB_n,          ONLY: NEB_t
USE MODD_IO,             ONLY: TFILEDATA
USE MODD_LES, ONLY: TLES_t
USE MODD_PARAMETERS, ONLY: JPVEXT
USE MODD_TURB_n,         ONLY: TURB_t
!
USE MODE_IO_FIELD_WRITE, ONLY: IO_FIELD_WRITE
!

 
USE MODI_LES_MEAN_SUBGRID, ONLY: LES_MEAN_SUBGRID
!
USE MODE_EMOIST, ONLY: EMOIST
USE MODE_ETHETA, ONLY: ETHETA
!
USE MODI_SECOND_MNH, ONLY: SECOND_MNH
USE MODE_SHUMAN_PHY, ONLY:DXM2D_PHY
USE MODE_SHUMAN_PHY, ONLY:DYM2D_PHY
USE MODE_GRADIENT_M_PHY, ONLY:GX_M_M_PHY
USE MODE_GRADIENT_M_PHY, ONLY:GY_M_M_PHY
USE MODE_SHUMAN_PHY, ONLY:MXF2D_PHY
USE MODE_SHUMAN_PHY, ONLY:MYF2D_PHY
USE MODE_SHUMAN_PHY, ONLY:MZF_PHY
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
REAL, DIMENSION(D%NIT,D%NJT,D%NKT)   :: ZFLX,ZWORK,ZA,ZWKLES ! work arrays
!   
INTEGER             :: IKB,IKE,IIT,IJT,IKT
INTEGER             :: JI,JJ,JK
                                    ! Index values for the Beginning and End
                                    ! mass points of the domain  
REAL, DIMENSION(D%NIT,D%NJT,1+JPVEXT:3+JPVEXT) :: ZCOEFF 
                                    ! coefficients for the uncentred gradient 
                                    ! computation near the ground
REAL :: ZTIME1, ZTIME2
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZGX_M_M3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZGY_M_M3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT) ::ZSHUGRADWK1_2D
REAL, DIMENSION(D%NIT,D%NJT) ::ZMXF2D_WORK1
REAL, DIMENSION(D%NIT,D%NJT) ::ZDXM2D_WORK1
REAL, DIMENSION(D%NIT,D%NJT) ::ZMXF2D_WORK2
REAL, DIMENSION(D%NIT,D%NJT) ::ZMYF2D_WORK1
REAL, DIMENSION(D%NIT,D%NJT) ::ZDYM2D_WORK1
REAL, DIMENSION(D%NIT,D%NJT) ::ZMYF2D_WORK2
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZMZF3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZGX_M_M3D_WORK2
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZGY_M_M3D_WORK2
REAL, DIMENSION(D%NIT,D%NJT) ::ZMXF2D_WORK3
REAL, DIMENSION(D%NIT,D%NJT) ::ZDXM2D_WORK2
REAL, DIMENSION(D%NIT,D%NJT) ::ZMXF2D_WORK4
REAL, DIMENSION(D%NIT,D%NJT) ::ZMYF2D_WORK3
REAL, DIMENSION(D%NIT,D%NJT) ::ZDYM2D_WORK2
REAL, DIMENSION(D%NIT,D%NJT) ::ZMYF2D_WORK4
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

ZCOEFF(:,:,IKB+2)= - PDZZ(:,:,IKB+1) /      &
       ( (PDZZ(:,:,IKB+2)+PDZZ(:,:,IKB+1)) * PDZZ(:,:,IKB+2) )
ZCOEFF(:,:,IKB+1)=   (PDZZ(:,:,IKB+2)+PDZZ(:,:,IKB+1)) /      &
       ( PDZZ(:,:,IKB+1) * PDZZ(:,:,IKB+2) )
ZCOEFF(:,:,IKB)= - (PDZZ(:,:,IKB+2)+2.*PDZZ(:,:,IKB+1)) /      &
       ( (PDZZ(:,:,IKB+2)+PDZZ(:,:,IKB+1)) * PDZZ(:,:,IKB+1) )

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
    CALL GX_M_M_PHY(D, OFLAT,PTHLM,PDXX,PDZZ,PDZX, ZGX_M_M3D_WORK1)
CALL GY_M_M_PHY(D, OFLAT,PTHLM,PDYY,PDZZ,PDZY, ZGY_M_M3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZFLX(:,:,:) = TURBN%XCTV * PLM(:,:,:) * PLEPS(:,:,:) *                           &
       ( ZGX_M_M3D_WORK1(:, :, :)**2 + ZGY_M_M3D_WORK1(:, :, :)**2 )  
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
ELSE
    CALL GX_M_M_PHY(D, OFLAT,PTHLM,PDXX,PDZZ,PDZX, ZGX_M_M3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZFLX(:,:,:) = TURBN%XCTV * PLM(:,:,:) * PLEPS(:,:,:) *                           &
         ZGX_M_M3D_WORK1(:, :, :)**2  
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
END IF
!
! Compute the flux at the first inner U-point with an uncentred vertical  
! gradient
!
  CALL DXM2D_PHY(D, PTHLM(:,:,IKB), ZDXM2D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
ZSHUGRADWK1_2D(:, :) = ZDXM2D_WORK1(:, :) * PINV_PDXX(:,:,IKB)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)

!
CALL MXF2D_PHY(D, ZSHUGRADWK1_2D, ZMXF2D_WORK1)
CALL MXF2D_PHY(D, PDXX(:,:,IKB), ZMXF2D_WORK2)
CALL DYM2D_PHY(D, PTHLM(:,:,IKB), ZDYM2D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
ZSHUGRADWK1_2D(:, :) = ZDYM2D_WORK1(:, :) * PINV_PDYY(:,:,IKB)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)

!
CALL MYF2D_PHY(D, ZSHUGRADWK1_2D, ZMYF2D_WORK1)
CALL MYF2D_PHY(D, PDYY(:,:,IKB), ZMYF2D_WORK2)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
ZFLX(:,:,IKB) = TURBN%XCTV * PLM(:,:,IKB)                  &
  * PLEPS(:,:,IKB) *  (                                    &
  ( ZMXF2D_WORK1(:, :)      &
   - ( ZCOEFF(:,:,IKB+2)*PTHLM(:,:,IKB+2)          &
      +ZCOEFF(:,:,IKB+1)*PTHLM(:,:,IKB+1)          &
      +ZCOEFF(:,:,IKB)*PTHLM(:,:,IKB)          &
     ) * 0.5 * ( PDZX(:,:,IKB+1)+PDZX(:,:,IKB) )     &
     / ZMXF2D_WORK2(:, :)                                  &
  ) ** 2 +                                                     &
  ( ZMYF2D_WORK1(:, :)      &
   - ( ZCOEFF(:,:,IKB+2)*PTHLM(:,:,IKB+2)          &
      +ZCOEFF(:,:,IKB+1)*PTHLM(:,:,IKB+1)          &
      +ZCOEFF(:,:,IKB)*PTHLM(:,:,IKB)          &
     ) * 0.5 * ( PDZY(:,:,IKB+1)+PDZY(:,:,IKB) )     &
     / ZMYF2D_WORK2(:, :)                                  &
  ) ** 2                                             )  
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)

!
!

  !$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
  ZFLX(:,:,IKB-1) = ZFLX(:,:,IKB)
  !$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)
  !
  IF ( KRRL > 0 ) THEN
!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
    ZWORK(:,:,:) = ZFLX(:,:,:) * PATHETA(:,:,:) * PATHETA(:,:,:)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
  END IF

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

    CALL IO_FIELD_WRITE(TPFILE,TZFIELD,ZFLX)
  END IF
!
! Storage in the LES configuration (addition to TURB_VER computation)
!
  IF (TLES%LLES_CALL) THEN
    CALL SECOND_MNH(ZTIME1)
    CALL LES_MEAN_SUBGRID( ZFLX, TLES%X_LES_SUBGRID_Thl2, .TRUE. ) 
    CALL MZF_PHY(D, PWM, ZMZF3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZWKLES(:, :, :) = ZMZF3D_WORK1(:, :, :)*ZFLX(:, :, :)    
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
CALL LES_MEAN_SUBGRID( ZWKLES, TLES%X_LES_RES_W_SBG_Thl2, .TRUE. )
    ZWKLES = -2.*XCTD*SQRT(PTKEM)*ZFLX/PLEPS
    CALL LES_MEAN_SUBGRID( ZWKLES, TLES%X_LES_SUBGRID_DISS_Thl2, .TRUE. )
    CALL ETHETA(D,CST,KRR,KRRI,PTHLM,PRM,PLOCPEXNM,PATHETA,PSRCM,OOCEAN,OCOMPUTE_SRC,ZA)
    ZWKLES = ZA*ZFLX
    CALL LES_MEAN_SUBGRID( ZWKLES, TLES%X_LES_SUBGRID_ThlThv, .TRUE. ) 
    ZWKLES = -CST%XG/PTHVREF/3.*ZA*ZFLX
    CALL LES_MEAN_SUBGRID( ZWKLES, TLES%X_LES_SUBGRID_ThlPz, .TRUE. )
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
      CALL GX_M_M_PHY(D, OFLAT,PTHLM,PDXX,PDZZ,PDZX, ZGX_M_M3D_WORK1)
CALL GX_M_M_PHY(D, OFLAT,PRM(:,:,:,1),PDXX,PDZZ,PDZX, ZGX_M_M3D_WORK2)
CALL GY_M_M_PHY(D, OFLAT,PTHLM,PDYY,PDZZ,PDZY, ZGY_M_M3D_WORK1)
CALL GY_M_M_PHY(D, OFLAT,PRM(:,:,:,1),PDYY,PDZZ,PDZY, ZGY_M_M3D_WORK2)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZFLX(:,:,:)=                                                               &
            PLM(:,:,:) * PLEPS(:,:,:) *                                          &
            (ZGX_M_M3D_WORK1(:, :, :) * ZGX_M_M3D_WORK2(:, :, :)  &
           + ZGY_M_M3D_WORK1(:, :, :) * ZGY_M_M3D_WORK2(:, :, :)  &
            ) * (TURBN%XCHT1+TURBN%XCHT2)    
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
ELSE
      CALL GX_M_M_PHY(D, OFLAT,PTHLM,PDXX,PDZZ,PDZX, ZGX_M_M3D_WORK1)
CALL GX_M_M_PHY(D, OFLAT,PRM(:,:,:,1),PDXX,PDZZ,PDZX, ZGX_M_M3D_WORK2)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZFLX(:,:,:)=                                                               &
            PLM(:,:,:) * PLEPS(:,:,:) *                                          &
            (ZGX_M_M3D_WORK1(:, :, :) * ZGX_M_M3D_WORK2(:, :, :)  &
            ) * (TURBN%XCHT1+TURBN%XCHT2)    
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
END IF
!
! Compute the flux at the first inner U-point with an uncentred vertical  
! gradient
    CALL DXM2D_PHY(D, PTHLM(:,:,IKB), ZDXM2D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
ZSHUGRADWK1_2D(:, :) = ZDXM2D_WORK1(:, :) * PINV_PDXX(:,:,IKB)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)

!
CALL MXF2D_PHY(D, ZSHUGRADWK1_2D, ZMXF2D_WORK1)
CALL MXF2D_PHY(D, PDXX(:,:,IKB), ZMXF2D_WORK2)
CALL DXM2D_PHY(D, PRM(:,:,IKB,1), ZDXM2D_WORK2)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
ZSHUGRADWK1_2D(:, :) = ZDXM2D_WORK2(:, :) * PINV_PDXX(:,:,IKB)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)

!
CALL MXF2D_PHY(D, ZSHUGRADWK1_2D, ZMXF2D_WORK3)
CALL MXF2D_PHY(D, PDXX(:,:,IKB), ZMXF2D_WORK4)
CALL DYM2D_PHY(D, PTHLM(:,:,IKB), ZDYM2D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
ZSHUGRADWK1_2D(:, :) = ZDYM2D_WORK1(:, :) * PINV_PDYY(:,:,IKB)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)

!
CALL MYF2D_PHY(D, ZSHUGRADWK1_2D, ZMYF2D_WORK1)
CALL MYF2D_PHY(D, PDYY(:,:,IKB), ZMYF2D_WORK2)
CALL DYM2D_PHY(D, PRM(:,:,IKB,1), ZDYM2D_WORK2)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
ZSHUGRADWK1_2D(:, :) = ZDYM2D_WORK2(:, :) * PINV_PDYY(:,:,IKB)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)

!
CALL MYF2D_PHY(D, ZSHUGRADWK1_2D, ZMYF2D_WORK3)
CALL MYF2D_PHY(D, PDYY(:,:,IKB), ZMYF2D_WORK4)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
ZFLX(:,:,IKB) = (TURBN%XCHT1+TURBN%XCHT2) * PLM(:,:,IKB)         &
    * PLEPS(:,:,IKB)  *  (                                   &
    ( ZMXF2D_WORK1(:, :)      &
     - ( ZCOEFF(:,:,IKB+2)*PTHLM(:,:,IKB+2)          &
        +ZCOEFF(:,:,IKB+1)*PTHLM(:,:,IKB+1)          &
        +ZCOEFF(:,:,IKB)*PTHLM(:,:,IKB)          &
       ) * 0.5 * ( PDZX(:,:,IKB+1)+PDZX(:,:,IKB) )     &
       / ZMXF2D_WORK2(:, :)                                  &
    ) *                                                          &
    ( ZMXF2D_WORK3(:, :)      &
     - ( ZCOEFF(:,:,IKB+2)*PRM(:,:,IKB+2,1)          &
        +ZCOEFF(:,:,IKB+1)*PRM(:,:,IKB+1,1)          &
        +ZCOEFF(:,:,IKB)*PRM(:,:,IKB,1)          &
       ) * 0.5 * ( PDZX(:,:,IKB+1)+PDZX(:,:,IKB) )     &
       / ZMXF2D_WORK4(:, :)                                  &
    ) +                                                          &
    ( ZMYF2D_WORK1(:, :)      &
     - ( ZCOEFF(:,:,IKB+2)*PTHLM(:,:,IKB+2)          &
        +ZCOEFF(:,:,IKB+1)*PTHLM(:,:,IKB+1)          &
        +ZCOEFF(:,:,IKB)*PTHLM(:,:,IKB)          &
       ) * 0.5 * ( PDZY(:,:,IKB+1)+PDZY(:,:,IKB) )     &
       / ZMYF2D_WORK2(:, :)                                  &
    ) *                                                          &
    ( ZMYF2D_WORK3(:, :)           &
     - ( ZCOEFF(:,:,IKB+2)*PRM(:,:,IKB+2,1)          &
        +ZCOEFF(:,:,IKB+1)*PRM(:,:,IKB+1,1)          &
        +ZCOEFF(:,:,IKB)*PRM(:,:,IKB,1)          &
       ) * 0.5 * ( PDZY(:,:,IKB+1)+PDZY(:,:,IKB) )     &
       / ZMYF2D_WORK4(:, :)                                  &
    )                                                          )    
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)

!
!

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
      CALL MZF_PHY(D, PWM, ZMZF3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZWKLES(:, :, :) = ZMZF3D_WORK1(:, :, :)*ZFLX(:, :, :)      
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
CALL LES_MEAN_SUBGRID( ZWKLES, TLES%X_LES_RES_W_SBG_ThlRt, .TRUE. )
      ZWKLES = -XCTD*SQRT(PTKEM)*ZFLX/PLEPS
      CALL LES_MEAN_SUBGRID( ZWKLES, TLES%X_LES_SUBGRID_DISS_ThlRt, .TRUE. )
      ZWKLES = ZA*ZFLX
      CALL LES_MEAN_SUBGRID( ZWKLES, TLES%X_LES_SUBGRID_RtThv, .TRUE. ) 
      ZWKLES = -CST%XG/PTHVREF/3.*ZA*ZFLX
      CALL LES_MEAN_SUBGRID( ZWKLES, TLES%X_LES_SUBGRID_RtPz,.TRUE.)
      CALL EMOIST(D,CST,KRR,KRRI,PTHLM,PRM,PLOCPEXNM,PAMOIST,PSRCM,OOCEAN,ZA)
      ZWKLES = ZA*ZFLX
      CALL LES_MEAN_SUBGRID( ZWKLES, TLES%X_LES_SUBGRID_ThlThv, .TRUE. ) 
      ZWKLES = -CST%XG/PTHVREF/3.*ZA*ZFLX
      CALL LES_MEAN_SUBGRID( ZWKLES, TLES%X_LES_SUBGRID_ThlPz,.TRUE.)
      CALL SECOND_MNH(ZTIME2)
      TLES%XTIME_LES = TLES%XTIME_LES + ZTIME2 - ZTIME1
    END IF
!
!*       8.4  <Rnp Rnp>
!
    ! Computes the horizontal variance <Rnp Rnp>
    IF (.NOT. O2D) THEN
      CALL GX_M_M_PHY(D, OFLAT,PRM(:,:,:,1),PDXX,PDZZ,PDZX, ZGX_M_M3D_WORK1)
CALL GY_M_M_PHY(D, OFLAT,PRM(:,:,:,1),PDYY,PDZZ,PDZY, ZGY_M_M3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZFLX(:,:,:) = TURBN%XCHV * PLM(:,:,:) * PLEPS(:,:,:) *                      &
           ( ZGX_M_M3D_WORK1(:, :, :)**2 +                       &
             ZGY_M_M3D_WORK1(:, :, :)**2 )    
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
ELSE
      CALL GX_M_M_PHY(D, OFLAT,PRM(:,:,:,1),PDXX,PDZZ,PDZX, ZGX_M_M3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZFLX(:,:,:) = TURBN%XCHV * PLM(:,:,:) * PLEPS(:,:,:) *                      &
           ( ZGX_M_M3D_WORK1(:, :, :)**2  )    
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
END IF
!
! Compute the flux at the first inner U-point with an uncentred vertical  
! gradient
    CALL DXM2D_PHY(D, PRM(:,:,IKB,1), ZDXM2D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
ZSHUGRADWK1_2D(:, :) = ZDXM2D_WORK1(:, :) * PINV_PDXX(:,:,IKB)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)

!
CALL MXF2D_PHY(D, ZSHUGRADWK1_2D, ZMXF2D_WORK1)
CALL MXF2D_PHY(D, PDXX(:,:,IKB:), ZMXF2D_WORK2)
CALL DYM2D_PHY(D, PRM(:,:,IKB,1), ZDYM2D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
ZSHUGRADWK1_2D(:, :) = ZDYM2D_WORK1(:, :) * PINV_PDYY(:,:,IKB)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)

!
CALL MYF2D_PHY(D, ZSHUGRADWK1_2D, ZMYF2D_WORK1)
CALL MYF2D_PHY(D, PDYY(:,:,IKB), ZMYF2D_WORK2)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
ZFLX(:,:,IKB) = TURBN%XCHV * PLM(:,:,IKB)                  &
    * PLEPS(:,:,IKB) *  (                                    &
    ( ZMXF2D_WORK1(:, :)      &
     - ( ZCOEFF(:,:,IKB+2)*PRM(:,:,IKB+2,1)          &
        +ZCOEFF(:,:,IKB+1)*PRM(:,:,IKB+1,1)          &
        +ZCOEFF(:,:,IKB)*PRM(:,:,IKB,1)          &
       ) * 0.5 * ( PDZX(:,:,IKB+1)+PDZX(:,:,IKB) )     &
       / ZMXF2D_WORK2(:, :)                                  &
    ) ** 2 +                                                     &
    ( ZMYF2D_WORK1(:, :)           &
     - ( ZCOEFF(:,:,IKB+2)*PRM(:,:,IKB+2,1)          &
        +ZCOEFF(:,:,IKB+1)*PRM(:,:,IKB+1,1)          &
        +ZCOEFF(:,:,IKB)*PRM(:,:,IKB,1)          &
       ) * 0.5 * ( PDZY(:,:,IKB+1)+PDZY(:,:,IKB) )     &
       / ZMYF2D_WORK2(:, :)                                  &
    ) ** 2                                             )
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)

!
!

    !$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
    ZFLX(:,:,IKB-1) = ZFLX(:,:,IKB)
    !$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)
    !
    IF ( KRRL > 0 ) THEN       
    !$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
      ZWORK(:,:,:) = ZWORK(:,:,:)+ PAMOIST(:,:,:) * PAMOIST(:,:,:) * ZFLX(:,:,:)
    !$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
    END IF

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

      CALL IO_FIELD_WRITE(TPFILE,TZFIELD,ZFLX)
    END IF
    !
    !   Storage in the LES configuration (addition to TURB_VER computation)
    !
    IF (TLES%LLES_CALL) THEN
      CALL SECOND_MNH(ZTIME1)
      CALL LES_MEAN_SUBGRID( ZFLX, TLES%X_LES_SUBGRID_Rt2, .TRUE. ) 
      CALL MZF_PHY(D, PWM, ZMZF3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZWKLES(:, :, :) = ZMZF3D_WORK1(:, :, :)*ZFLX(:, :, :)      
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
CALL LES_MEAN_SUBGRID( ZWKLES, TLES%X_LES_RES_W_SBG_Rt2, .TRUE. )
      ZWKLES = ZA*ZFLX
      CALL LES_MEAN_SUBGRID( ZWKLES, TLES%X_LES_SUBGRID_RtThv, .TRUE. ) 
      ZWKLES = -CST%XG/PTHVREF/3.*ZA*ZFLX
      CALL LES_MEAN_SUBGRID( ZWKLES, TLES%X_LES_SUBGRID_RtPz,.TRUE.)
      ZWKLES = -2.*XCTD*SQRT(PTKEM)*ZFLX/PLEPS
      CALL LES_MEAN_SUBGRID( ZWKLES, TLES%X_LES_SUBGRID_DISS_Rt2, .TRUE. )
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
    
    PSIGS(:,:,:)=PSIGS(:,:,:)*PSIGS(:,:,:) + ZWORK(:,:,:)
    ! Extrapolate PSIGS at the ground and at the top
    PSIGS(:,:,IKB-1) = PSIGS(:,:,IKB)
    PSIGS(:,:,IKE+1) = PSIGS(:,:,IKE)
    PSIGS(:,:,:) = SQRT(MAX ( PSIGS(:,:,:),1.E-12) ) 
    
  END IF       
!
END IF
!
!
!
END SUBROUTINE TURB_HOR_THERMO_CORR
END MODULE MODE_TURB_HOR_THERMO_CORR
