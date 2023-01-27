!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_TURB_HOR_THERMO_CORR
IMPLICIT NONE
CONTAINS
      SUBROUTINE TURB_HOR_THERMO_CORR(D,CST,TURBN,TLES,              &
                      KRR, KRRL, KRRI,                               &
                      OOCEAN,OCOMPUTE_SRC,O2D,                       &
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
TYPE(TLES_t),           INTENT(INOUT):: TLES          ! modd_les structure
INTEGER,                INTENT(IN)   :: KRR           ! number of moist var.
INTEGER,                INTENT(IN)   :: KRRL          ! number of liquid water var.
INTEGER,                INTENT(IN)   :: KRRI          ! number of ice water var.
LOGICAL,                INTENT(IN)   ::  OOCEAN       ! switch for Ocean model version
LOGICAL,                INTENT(IN)   ::  OCOMPUTE_SRC ! flag to define dimensions of SIGS and SRCT variables
LOGICAL,                INTENT(IN)   ::  O2D          ! Logical for 2D model version (modd_conf)
TYPE(TFILEDATA),          INTENT(IN)    ::  TPFILE       ! Output file
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PINV_PDXX   ! 1./PDXX
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PINV_PDYY   ! 1./PDYY
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PDXX, PDYY, PDZZ, PDZX, PDZY 
                                                         ! Metric coefficients
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PTHVREF      ! ref. state Virtual 
                                                      ! Potential Temperature
!
! Variables at t-1
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PWM 
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PTHLM 
REAL, DIMENSION(:,:,:,:), INTENT(IN)    ::  PRM          ! mixing ratios at t-1,
                              !  where PRM(:,:,:,1) = conservative mixing ratio
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PTKEM        ! Turb. Kin. Energy
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PLM          ! Turb. mixing length
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PLEPS        ! dissipative length
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PLOCPEXNM    ! Lv(T)/Cp/Exnref at time t-1
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PATHETA      ! coefficients between 
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PAMOIST      ! s and Thetal and Rnp
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PSRCM        ! normalized 
!
!
!
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) ::  PSIGS
                                  ! IN: Vertical part of Sigma_s at t
                                  ! OUT: Total Sigma_s at t
!
!*       0.2  declaration of local variables
!
REAL, DIMENSION(SIZE(PTHLM,1),SIZE(PTHLM,2),SIZE(PTHLM,3))       &
                                     :: ZFLX,ZWORK,ZA
    ! work arrays
!   
INTEGER             :: IKB,IKE
                                    ! Index values for the Beginning and End
                                    ! mass points of the domain  
REAL, DIMENSION(SIZE(PDZZ,1),SIZE(PDZZ,2),1+JPVEXT:3+JPVEXT) :: ZCOEFF 
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
IF ( ( KRRL > 0 .AND. TURBN%LSUBG_COND) .OR. ( TURBN%LTURB_FLX .AND. TPFILE%LOPENED ) &
                                  .OR. ( TLES%LLES_CALL )                  ) THEN
!
!*       8.1  <THl THl>
!
  ! Computes the horizontal variance <THl THl>
  IF (.NOT. O2D) THEN
    ZFLX(:,:,:) = XCTV * PLM(:,:,:) * PLEPS(:,:,:) *                           &
       ( GX_M_M(PTHLM,PDXX,PDZZ,PDZX)**2 + GY_M_M(PTHLM,PDYY,PDZZ,PDZY)**2 )
  ELSE
    ZFLX(:,:,:) = XCTV * PLM(:,:,:) * PLEPS(:,:,:) *                           &
         GX_M_M(PTHLM,PDXX,PDZZ,PDZX)**2
  END IF
!
! Compute the flux at the first inner U-point with an uncentred vertical  
! gradient
!
  ZFLX(:,:,IKB:IKB) = XCTV * PLM(:,:,IKB:IKB)                  &
  * PLEPS(:,:,IKB:IKB) *  (                                    &
  ( MXF(DXM(PTHLM(:,:,IKB:IKB)) * PINV_PDXX(:,:,IKB:IKB))      &
   - ( ZCOEFF(:,:,IKB+2:IKB+2)*PTHLM(:,:,IKB+2:IKB+2)          &
      +ZCOEFF(:,:,IKB+1:IKB+1)*PTHLM(:,:,IKB+1:IKB+1)          &
      +ZCOEFF(:,:,IKB  :IKB  )*PTHLM(:,:,IKB  :IKB  )          &
     ) * 0.5 * ( PDZX(:,:,IKB+1:IKB+1)+PDZX(:,:,IKB:IKB) )     &
     / MXF(PDXX(:,:,IKB:IKB))                                  &
  ) ** 2 +                                                     &
  ( MYF(DYM(PTHLM(:,:,IKB:IKB)) * PINV_PDYY(:,:,IKB:IKB))      &
   - ( ZCOEFF(:,:,IKB+2:IKB+2)*PTHLM(:,:,IKB+2:IKB+2)          &
      +ZCOEFF(:,:,IKB+1:IKB+1)*PTHLM(:,:,IKB+1:IKB+1)          &
      +ZCOEFF(:,:,IKB  :IKB  )*PTHLM(:,:,IKB  :IKB  )          &
     ) * 0.5 * ( PDZY(:,:,IKB+1:IKB+1)+PDZY(:,:,IKB:IKB) )     &
     / MYF(PDYY(:,:,IKB:IKB))                                  &
  ) ** 2                                             )
  !
  ZFLX(:,:,IKB-1) = ZFLX(:,:,IKB)
  !
  IF ( KRRL > 0 ) THEN
    ZWORK(:,:,:) = ZFLX(:,:,:) * PATHETA(:,:,:) * PATHETA(:,:,:)
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
            (GX_M_M(PTHLM,PDXX,PDZZ,PDZX) * GX_M_M(PRM(:,:,:,1),PDXX,PDZZ,PDZX)  &
           + GY_M_M(PTHLM,PDYY,PDZZ,PDZY) * GY_M_M(PRM(:,:,:,1),PDYY,PDZZ,PDZY)  &
            ) * (XCHT1+XCHT2)
    ELSE
      ZFLX(:,:,:)=                                                               &
            PLM(:,:,:) * PLEPS(:,:,:) *                                          &
            (GX_M_M(PTHLM,PDXX,PDZZ,PDZX) * GX_M_M(PRM(:,:,:,1),PDXX,PDZZ,PDZX)  &
            ) * (XCHT1+XCHT2)

    END IF
!
! Compute the flux at the first inner U-point with an uncentred vertical  
! gradient
    ZFLX(:,:,IKB:IKB) = (XCHT1+XCHT2) * PLM(:,:,IKB:IKB)         &
    * PLEPS(:,:,IKB:IKB)  *  (                                   &
    ( MXF(DXM(PTHLM(:,:,IKB:IKB)) * PINV_PDXX(:,:,IKB:IKB))      &
     - ( ZCOEFF(:,:,IKB+2:IKB+2)*PTHLM(:,:,IKB+2:IKB+2)          &
        +ZCOEFF(:,:,IKB+1:IKB+1)*PTHLM(:,:,IKB+1:IKB+1)          &
        +ZCOEFF(:,:,IKB  :IKB  )*PTHLM(:,:,IKB  :IKB  )          &
       ) * 0.5 * ( PDZX(:,:,IKB+1:IKB+1)+PDZX(:,:,IKB:IKB) )     &
       / MXF(PDXX(:,:,IKB:IKB))                                  &
    ) *                                                          &
    ( MXF(DXM(PRM(:,:,IKB:IKB,1)) * PINV_PDXX(:,:,IKB:IKB))      &
     - ( ZCOEFF(:,:,IKB+2:IKB+2)*PRM(:,:,IKB+2:IKB+2,1)          &
        +ZCOEFF(:,:,IKB+1:IKB+1)*PRM(:,:,IKB+1:IKB+1,1)          &
        +ZCOEFF(:,:,IKB  :IKB  )*PRM(:,:,IKB  :IKB  ,1)          &
       ) * 0.5 * ( PDZX(:,:,IKB+1:IKB+1)+PDZX(:,:,IKB:IKB) )     &
       / MXF(PDXX(:,:,IKB:IKB))                                  &
    ) +                                                          &
    ( MYF(DYM(PTHLM(:,:,IKB:IKB)) * PINV_PDYY(:,:,IKB:IKB))      &
     - ( ZCOEFF(:,:,IKB+2:IKB+2)*PTHLM(:,:,IKB+2:IKB+2)          &
        +ZCOEFF(:,:,IKB+1:IKB+1)*PTHLM(:,:,IKB+1:IKB+1)          &
        +ZCOEFF(:,:,IKB  :IKB  )*PTHLM(:,:,IKB  :IKB  )          &
       ) * 0.5 * ( PDZY(:,:,IKB+1:IKB+1)+PDZY(:,:,IKB:IKB) )     &
       / MYF(PDYY(:,:,IKB:IKB))                                  &
    ) *                                                          &
    ( MYF(DYM(PRM(:,:,IKB:IKB,1)) * PINV_PDYY(:,:,IKB:IKB))           &
     - ( ZCOEFF(:,:,IKB+2:IKB+2)*PRM(:,:,IKB+2:IKB+2,1)          &
        +ZCOEFF(:,:,IKB+1:IKB+1)*PRM(:,:,IKB+1:IKB+1,1)          &
        +ZCOEFF(:,:,IKB  :IKB  )*PRM(:,:,IKB  :IKB  ,1)          &
       ) * 0.5 * ( PDZY(:,:,IKB+1:IKB+1)+PDZY(:,:,IKB:IKB) )     &
       / MYF(PDYY(:,:,IKB:IKB))                                  &
    )                                                          )
    !
    ZFLX(:,:,IKB-1) = ZFLX(:,:,IKB)
    !
    IF ( KRRL > 0 )  THEN
      ZWORK(:,:,:) = ZWORK(:,:,:) +       &
                     2. * PATHETA(:,:,:) * PAMOIST(:,:,:) * ZFLX(:,:,:)    
    END IF                    
    !
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
      ZFLX(:,:,:) = XCHV * PLM(:,:,:) * PLEPS(:,:,:) *                      &
           ( GX_M_M(PRM(:,:,:,1),PDXX,PDZZ,PDZX)**2 +                       &
             GY_M_M(PRM(:,:,:,1),PDYY,PDZZ,PDZY)**2 )
    ELSE
      ZFLX(:,:,:) = XCHV * PLM(:,:,:) * PLEPS(:,:,:) *                      &
           ( GX_M_M(PRM(:,:,:,1),PDXX,PDZZ,PDZX)**2  )
    END IF
!
! Compute the flux at the first inner U-point with an uncentred vertical  
! gradient
    ZFLX(:,:,IKB:IKB) = XCHV * PLM(:,:,IKB:IKB)                  &
    * PLEPS(:,:,IKB:IKB) *  (                                    &
    ( MXF(DXM(PRM(:,:,IKB:IKB,1)) * PINV_PDXX(:,:,IKB:IKB))      &
     - ( ZCOEFF(:,:,IKB+2:IKB+2)*PRM(:,:,IKB+2:IKB+2,1)          &
        +ZCOEFF(:,:,IKB+1:IKB+1)*PRM(:,:,IKB+1:IKB+1,1)          &
        +ZCOEFF(:,:,IKB  :IKB  )*PRM(:,:,IKB  :IKB  ,1)          &
       ) * 0.5 * ( PDZX(:,:,IKB+1:IKB+1)+PDZX(:,:,IKB:IKB) )     &
       / MXF(PDXX(:,:,IKB:IKB))                                  &
    ) ** 2 +                                                     &
    ( MYF(DYM(PRM(:,:,IKB:IKB,1)) * PINV_PDYY(:,:,IKB:IKB))           &
     - ( ZCOEFF(:,:,IKB+2:IKB+2)*PRM(:,:,IKB+2:IKB+2,1)          &
        +ZCOEFF(:,:,IKB+1:IKB+1)*PRM(:,:,IKB+1:IKB+1,1)          &
        +ZCOEFF(:,:,IKB  :IKB  )*PRM(:,:,IKB  :IKB  ,1)          &
       ) * 0.5 * ( PDZY(:,:,IKB+1:IKB+1)+PDZY(:,:,IKB:IKB) )     &
       / MYF(PDYY(:,:,IKB:IKB))                                  &
    ) ** 2                                             )
!
    ZFLX(:,:,IKB-1) = ZFLX(:,:,IKB)
    !
    IF ( KRRL > 0 ) THEN       
      ZWORK(:,:,:) = ZWORK(:,:,:)+ PAMOIST(:,:,:) * PAMOIST(:,:,:) * ZFLX(:,:,:)
    END IF
    !
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
