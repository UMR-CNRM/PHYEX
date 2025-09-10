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
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZGX_M_M3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZGY_M_M3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT) ::ZSHUGRADWK1_2D
REAL, DIMENSION(D%NIT,D%NJT) ::ZMXF2D_WORK1
REAL, DIMENSION(D%NIT,D%NJT) ::ZDXM2D_WORK1
REAL, DIMENSION(D%NIT,D%NJT) ::ZMXF2D_WORK2
REAL, DIMENSION(D%NIT,D%NJT) ::ZMYF2D_WORK1
REAL, DIMENSION(D%NIT,D%NJT) ::ZDYM2D_WORK1
REAL, DIMENSION(D%NIT,D%NJT) ::ZMYF2D_WORK2
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZGX_M_M3D_WORK2
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZGY_M_M3D_WORK2
REAL, DIMENSION(D%NIT,D%NJT) ::ZMXF2D_WORK3
REAL, DIMENSION(D%NIT,D%NJT) ::ZDXM2D_WORK2
REAL, DIMENSION(D%NIT,D%NJT) ::ZMXF2D_WORK4
REAL, DIMENSION(D%NIT,D%NJT) ::ZMYF2D_WORK3
REAL, DIMENSION(D%NIT,D%NJT) ::ZDYM2D_WORK2
REAL, DIMENSION(D%NIT,D%NJT) ::ZMYF2D_WORK4
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZMZF3D_WORK1
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
    CALL GX_M_M_PHY(D, OFLAT,PTHLM(:, :, :),PDXX(:, :, :),PDZZ(:, :, :),PDZX(:, :, :), ZGX_M_M3D_WORK1)
CALL GY_M_M_PHY(D, OFLAT,PTHLM(:, :, :),PDYY(:, :, :),PDZZ(:, :, :),PDZY(:, :, :), ZGY_M_M3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZFLX(JI, JJ, JK) = TURBN%XCTV * PLM(JI, JJ, JK) * PLEPS(JI, JJ, JK) *                           &
             ( ZGX_M_M3D_WORK1(JI, JJ, JK)**2 + ZGY_M_M3D_WORK1(JI, JJ, JK)**2 )  
    END DO
  END DO
END DO

!
ELSE
    CALL GX_M_M_PHY(D, OFLAT,PTHLM(:, :, :),PDXX(:, :, :),PDZZ(:, :, :),PDZX(:, :, :), ZGX_M_M3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZFLX(JI, JJ, JK) = TURBN%XCTV * PLM(JI, JJ, JK) * PLEPS(JI, JJ, JK) *                           &
               ZGX_M_M3D_WORK1(JI, JJ, JK)**2  
    END DO
  END DO
END DO

!
END IF
!
! Compute the flux at the first inner U-point with an uncentred vertical  
! gradient
!
  CALL DXM2D_PHY(D, PTHLM(:,:,IKB), ZDXM2D_WORK1)

DO JJ=1, IJT
  DO JI=1, IIT
    ZSHUGRADWK1_2D(JI, JJ) = ZDXM2D_WORK1(JI, JJ) * PINV_PDXX(JI, JJ, IKB)
  END DO
END DO

!
CALL MXF2D_PHY(D, ZSHUGRADWK1_2D, ZMXF2D_WORK1)
CALL MXF2D_PHY(D, PDXX(:,:,IKB), ZMXF2D_WORK2)
CALL DYM2D_PHY(D, PTHLM(:,:,IKB), ZDYM2D_WORK1)

DO JJ=1, IJT
  DO JI=1, IIT
    ZSHUGRADWK1_2D(JI, JJ) = ZDYM2D_WORK1(JI, JJ) * PINV_PDYY(JI, JJ, IKB)
  END DO
END DO

!
CALL MYF2D_PHY(D, ZSHUGRADWK1_2D, ZMYF2D_WORK1)
CALL MYF2D_PHY(D, PDYY(:,:,IKB), ZMYF2D_WORK2)

DO JJ=1, IJT
  DO JI=1, IIT
    ZFLX(JI, JJ, IKB) = TURBN%XCTV * PLM(JI, JJ, IKB)                  &
      * PLEPS(JI, JJ, IKB) *  (                                    &
      ( ZMXF2D_WORK1(JI, JJ)      &
       - ( ZCOEFF(JI, JJ, IKB+2)*PTHLM(JI, JJ, IKB+2)          &
          +ZCOEFF(JI, JJ, IKB+1)*PTHLM(JI, JJ, IKB+1)          &
          +ZCOEFF(JI, JJ, IKB)*PTHLM(JI, JJ, IKB)          &
         ) * 0.5 * ( PDZX(JI, JJ, IKB+1)+PDZX(JI, JJ, IKB) )     &
         / ZMXF2D_WORK2(JI, JJ)                                  &
      ) ** 2 +                                                     &
      ( ZMYF2D_WORK1(JI, JJ)      &
       - ( ZCOEFF(JI, JJ, IKB+2)*PTHLM(JI, JJ, IKB+2)          &
          +ZCOEFF(JI, JJ, IKB+1)*PTHLM(JI, JJ, IKB+1)          &
          +ZCOEFF(JI, JJ, IKB)*PTHLM(JI, JJ, IKB)          &
         ) * 0.5 * ( PDZY(JI, JJ, IKB+1)+PDZY(JI, JJ, IKB) )     &
         / ZMYF2D_WORK2(JI, JJ)                                  &
      ) ** 2                                             )  
  END DO
END DO

!
!

  DO JJ=1, IJT
    DO JI=1, IIT
      ZFLX(JI, JJ, IKB-1) = ZFLX(JI, JJ, IKB)
    END DO
  END DO
  !
  IF ( KRRL > 0 ) THEN
DO JK=1, IKT
      DO JJ=1, IJT
        DO JI=1, IIT
          ZWORK(JI, JJ, JK) = ZFLX(JI, JJ, JK) * PATHETA(JI, JJ, JK) * PATHETA(JI, JJ, JK)
    END DO
      END DO
    END DO
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
    CALL MZF_PHY(D, PWM(:, :, :), ZMZF3D_WORK1)
    CALL LES_MEAN_SUBGRID( ZMZF3D_WORK1(:, :, :)*ZFLX(:, :, :), TLES%X_LES_RES_W_SBG_Thl2, .TRUE. )    
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
      CALL GX_M_M_PHY(D, OFLAT,PTHLM(:, :, :),PDXX(:, :, :),PDZZ(:, :, :),PDZX(:, :, :), ZGX_M_M3D_WORK1)
CALL GX_M_M_PHY(D, OFLAT,PRM(:,:,:,1),PDXX(:, :, :),PDZZ(:, :, :),PDZX(:, :, :), ZGX_M_M3D_WORK2)
CALL GY_M_M_PHY(D, OFLAT,PTHLM(:, :, :),PDYY(:, :, :),PDZZ(:, :, :),PDZY(:, :, :), ZGY_M_M3D_WORK1)
CALL GY_M_M_PHY(D, OFLAT,PRM(:,:,:,1),PDYY(:, :, :),PDZZ(:, :, :),PDZY(:, :, :), ZGY_M_M3D_WORK2)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZFLX(JI, JJ, JK)=                                                               &
                  PLM(JI, JJ, JK) * PLEPS(JI, JJ, JK) *                                          &
                  (ZGX_M_M3D_WORK1(JI, JJ, JK) * ZGX_M_M3D_WORK2(JI, JJ, JK)  &
                 + ZGY_M_M3D_WORK1(JI, JJ, JK) * ZGY_M_M3D_WORK2(JI, JJ, JK)  &
                  ) * (TURBN%XCHT1+TURBN%XCHT2)    
    END DO
  END DO
END DO

!
ELSE
      CALL GX_M_M_PHY(D, OFLAT,PTHLM(:, :, :),PDXX(:, :, :),PDZZ(:, :, :),PDZX(:, :, :), ZGX_M_M3D_WORK1)
CALL GX_M_M_PHY(D, OFLAT,PRM(:,:,:,1),PDXX(:, :, :),PDZZ(:, :, :),PDZX(:, :, :), ZGX_M_M3D_WORK2)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZFLX(JI, JJ, JK)=                                                               &
                  PLM(JI, JJ, JK) * PLEPS(JI, JJ, JK) *                                          &
                  (ZGX_M_M3D_WORK1(JI, JJ, JK) * ZGX_M_M3D_WORK2(JI, JJ, JK)  &
                  ) * (TURBN%XCHT1+TURBN%XCHT2)    
    END DO
  END DO
END DO

!
END IF
!
! Compute the flux at the first inner U-point with an uncentred vertical  
! gradient
    CALL DXM2D_PHY(D, PTHLM(:,:,IKB), ZDXM2D_WORK1)

DO JJ=1, IJT
  DO JI=1, IIT
    ZSHUGRADWK1_2D(JI, JJ) = ZDXM2D_WORK1(JI, JJ) * PINV_PDXX(JI, JJ, IKB)
  END DO
END DO

!
CALL MXF2D_PHY(D, ZSHUGRADWK1_2D, ZMXF2D_WORK1)
CALL MXF2D_PHY(D, PDXX(:,:,IKB), ZMXF2D_WORK2)
CALL DXM2D_PHY(D, PRM(:,:,IKB,1), ZDXM2D_WORK2)

DO JJ=1, IJT
  DO JI=1, IIT
    ZSHUGRADWK1_2D(JI, JJ) = ZDXM2D_WORK2(JI, JJ) * PINV_PDXX(JI, JJ, IKB)
  END DO
END DO

!
CALL MXF2D_PHY(D, ZSHUGRADWK1_2D, ZMXF2D_WORK3)
CALL MXF2D_PHY(D, PDXX(:,:,IKB), ZMXF2D_WORK4)
CALL DYM2D_PHY(D, PTHLM(:,:,IKB), ZDYM2D_WORK1)

DO JJ=1, IJT
  DO JI=1, IIT
    ZSHUGRADWK1_2D(JI, JJ) = ZDYM2D_WORK1(JI, JJ) * PINV_PDYY(JI, JJ, IKB)
  END DO
END DO

!
CALL MYF2D_PHY(D, ZSHUGRADWK1_2D, ZMYF2D_WORK1)
CALL MYF2D_PHY(D, PDYY(:,:,IKB), ZMYF2D_WORK2)
CALL DYM2D_PHY(D, PRM(:,:,IKB,1), ZDYM2D_WORK2)

DO JJ=1, IJT
  DO JI=1, IIT
    ZSHUGRADWK1_2D(JI, JJ) = ZDYM2D_WORK2(JI, JJ) * PINV_PDYY(JI, JJ, IKB)
  END DO
END DO

!
CALL MYF2D_PHY(D, ZSHUGRADWK1_2D, ZMYF2D_WORK3)
CALL MYF2D_PHY(D, PDYY(:,:,IKB), ZMYF2D_WORK4)

DO JJ=1, IJT
  DO JI=1, IIT
    ZFLX(JI, JJ, IKB) = (TURBN%XCHT1+TURBN%XCHT2) * PLM(JI, JJ, IKB)         &
        * PLEPS(JI, JJ, IKB)  *  (                                   &
        ( ZMXF2D_WORK1(JI, JJ)      &
         - ( ZCOEFF(JI, JJ, IKB+2)*PTHLM(JI, JJ, IKB+2)          &
            +ZCOEFF(JI, JJ, IKB+1)*PTHLM(JI, JJ, IKB+1)          &
            +ZCOEFF(JI, JJ, IKB)*PTHLM(JI, JJ, IKB)          &
           ) * 0.5 * ( PDZX(JI, JJ, IKB+1)+PDZX(JI, JJ, IKB) )     &
           / ZMXF2D_WORK2(JI, JJ)                                  &
        ) *                                                          &
        ( ZMXF2D_WORK3(JI, JJ)      &
         - ( ZCOEFF(JI, JJ, IKB+2)*PRM(JI, JJ, IKB+2, 1)          &
            +ZCOEFF(JI, JJ, IKB+1)*PRM(JI, JJ, IKB+1, 1)          &
            +ZCOEFF(JI, JJ, IKB)*PRM(JI, JJ, IKB, 1)          &
           ) * 0.5 * ( PDZX(JI, JJ, IKB+1)+PDZX(JI, JJ, IKB) )     &
           / ZMXF2D_WORK4(JI, JJ)                                  &
        ) +                                                          &
        ( ZMYF2D_WORK1(JI, JJ)      &
         - ( ZCOEFF(JI, JJ, IKB+2)*PTHLM(JI, JJ, IKB+2)          &
            +ZCOEFF(JI, JJ, IKB+1)*PTHLM(JI, JJ, IKB+1)          &
            +ZCOEFF(JI, JJ, IKB)*PTHLM(JI, JJ, IKB)          &
           ) * 0.5 * ( PDZY(JI, JJ, IKB+1)+PDZY(JI, JJ, IKB) )     &
           / ZMYF2D_WORK2(JI, JJ)                                  &
        ) *                                                          &
        ( ZMYF2D_WORK3(JI, JJ)           &
         - ( ZCOEFF(JI, JJ, IKB+2)*PRM(JI, JJ, IKB+2, 1)          &
            +ZCOEFF(JI, JJ, IKB+1)*PRM(JI, JJ, IKB+1, 1)          &
            +ZCOEFF(JI, JJ, IKB)*PRM(JI, JJ, IKB, 1)          &
           ) * 0.5 * ( PDZY(JI, JJ, IKB+1)+PDZY(JI, JJ, IKB) )     &
           / ZMYF2D_WORK4(JI, JJ)                                  &
        )                                                          )    
  END DO
END DO

!
!

    DO JJ=1, IJT
      DO JI=1, IIT
        ZFLX(JI, JJ, IKB-1) = ZFLX(JI, JJ, IKB)
      END DO
    END DO
    !
    IF ( KRRL > 0 )  THEN
DO JK=1, IKT
        DO JJ=1, IJT
          DO JI=1, IIT
            ZWORK(JI, JJ, JK) = ZWORK(JI, JJ, JK) +       &
                           2. * PATHETA(JI, JJ, JK) * PAMOIST(JI, JJ, JK) * ZFLX(JI, JJ, JK)    
    END DO
        END DO
      END DO
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
      CALL MZF_PHY(D, PWM(:, :, :), ZMZF3D_WORK1)
CALL LES_MEAN_SUBGRID( ZMZF3D_WORK1(:, :, :)*ZFLX(:, :, :), TLES%X_LES_RES_W_SBG_ThlRt, .TRUE. )      
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
      CALL GX_M_M_PHY(D, OFLAT,PRM(:,:,:,1),PDXX(:, :, :),PDZZ(:, :, :),PDZX(:, :, :), ZGX_M_M3D_WORK1)
CALL GY_M_M_PHY(D, OFLAT,PRM(:,:,:,1),PDYY(:, :, :),PDZZ(:, :, :),PDZY(:, :, :), ZGY_M_M3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZFLX(JI, JJ, JK) = TURBN%XCHV * PLM(JI, JJ, JK) * PLEPS(JI, JJ, JK) *                      &
                 ( ZGX_M_M3D_WORK1(JI, JJ, JK)**2 +                       &
                   ZGY_M_M3D_WORK1(JI, JJ, JK)**2 )    
    END DO
  END DO
END DO

!
ELSE
      CALL GX_M_M_PHY(D, OFLAT,PRM(:,:,:,1),PDXX(:, :, :),PDZZ(:, :, :),PDZX(:, :, :), ZGX_M_M3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZFLX(JI, JJ, JK) = TURBN%XCHV * PLM(JI, JJ, JK) * PLEPS(JI, JJ, JK) *                      &
                 ( ZGX_M_M3D_WORK1(JI, JJ, JK)**2  )    
    END DO
  END DO
END DO

!
END IF
!
! Compute the flux at the first inner U-point with an uncentred vertical  
! gradient
    CALL DXM2D_PHY(D, PRM(:,:,IKB,1), ZDXM2D_WORK1)

DO JJ=1, IJT
  DO JI=1, IIT
    ZSHUGRADWK1_2D(JI, JJ) = ZDXM2D_WORK1(JI, JJ) * PINV_PDXX(JI, JJ, IKB)
  END DO
END DO

!
CALL MXF2D_PHY(D, ZSHUGRADWK1_2D, ZMXF2D_WORK1)
CALL MXF2D_PHY(D, PDXX(:,:,IKB:), ZMXF2D_WORK2)
CALL DYM2D_PHY(D, PRM(:,:,IKB,1), ZDYM2D_WORK1)

DO JJ=1, IJT
  DO JI=1, IIT
    ZSHUGRADWK1_2D(JI, JJ) = ZDYM2D_WORK1(JI, JJ) * PINV_PDYY(JI, JJ, IKB)
  END DO
END DO

!
CALL MYF2D_PHY(D, ZSHUGRADWK1_2D, ZMYF2D_WORK1)
CALL MYF2D_PHY(D, PDYY(:,:,IKB), ZMYF2D_WORK2)

DO JJ=1, IJT
  DO JI=1, IIT
    ZFLX(JI, JJ, IKB) = TURBN%XCHV * PLM(JI, JJ, IKB)                  &
        * PLEPS(JI, JJ, IKB) *  (                                    &
        ( ZMXF2D_WORK1(JI, JJ)      &
         - ( ZCOEFF(JI, JJ, IKB+2)*PRM(JI, JJ, IKB+2, 1)          &
            +ZCOEFF(JI, JJ, IKB+1)*PRM(JI, JJ, IKB+1, 1)          &
            +ZCOEFF(JI, JJ, IKB)*PRM(JI, JJ, IKB, 1)          &
           ) * 0.5 * ( PDZX(JI, JJ, IKB+1)+PDZX(JI, JJ, IKB) )     &
           / ZMXF2D_WORK2(JI, JJ)                                  &
        ) ** 2 +                                                     &
        ( ZMYF2D_WORK1(JI, JJ)           &
         - ( ZCOEFF(JI, JJ, IKB+2)*PRM(JI, JJ, IKB+2, 1)          &
            +ZCOEFF(JI, JJ, IKB+1)*PRM(JI, JJ, IKB+1, 1)          &
            +ZCOEFF(JI, JJ, IKB)*PRM(JI, JJ, IKB, 1)          &
           ) * 0.5 * ( PDZY(JI, JJ, IKB+1)+PDZY(JI, JJ, IKB) )     &
           / ZMYF2D_WORK2(JI, JJ)                                  &
        ) ** 2                                             )
  END DO
END DO

!
!

    DO JJ=1, IJT
      DO JI=1, IIT
        ZFLX(JI, JJ, IKB-1) = ZFLX(JI, JJ, IKB)
      END DO
    END DO
    !
    IF ( KRRL > 0 ) THEN       
    DO JK=1, IKT
        DO JJ=1, IJT
          DO JI=1, IIT
            ZWORK(JI, JJ, JK) = ZWORK(JI, JJ, JK)+ PAMOIST(JI, JJ, JK) * PAMOIST(JI, JJ, JK) * ZFLX(JI, JJ, JK)
        END DO
        END DO
      END DO
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
      CALL MZF_PHY(D, PWM(:, :, :), ZMZF3D_WORK1)
CALL LES_MEAN_SUBGRID( ZMZF3D_WORK1(:, :, :)*ZFLX(:, :, :), TLES%X_LES_RES_W_SBG_Rt2, .TRUE. )      
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
