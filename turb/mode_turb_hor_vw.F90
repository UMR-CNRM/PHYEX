!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_TURB_HOR_VW
IMPLICIT NONE
CONTAINS
      SUBROUTINE TURB_HOR_VW(TURBN,TLES,KSPLT,                       &
                      KRR,KSV,OFLAT,O2D,                             &
                      TPFILE,                                        &
                      PK,PINV_PDYY,PINV_PDZZ,PMZM_PRHODJ,            &
                      PDYY,PDZZ,PDZY,                                &
                      PRHODJ,PTHVREF,                                &
                      PVM,PWM,PTHLM,PRM,PSVM,                        &
                      PTKEM,PLM,                                     &
                      PDP,                                           &
                      PRVS,PRWS                                      )
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
!!                     Oct  18, 2000 (V. Masson) LES computations + OFLAT switch
!!                     Feb  14, 2001 (V. Masson and J. Stein) DZF bug on PRWS
!!                                + remove the use of W=0 at the ground
!!                                + extrapolataion under the ground
!!                     Nov  06, 2002 (V. Masson) LES budgets
!!                     October 2009 (G. Tanguy) add ILENCH=LEN(YCOMMENT) after
!!                                              change of YCOMMENT
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!! --------------------------------------------------------------------------
!
!*      0. DECLARATIONS
!          ------------
!
USE MODD_TURB_n, ONLY: TURB_t
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
USE MODI_SHUMAN 
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
TYPE(TURB_t),             INTENT(IN)    :: TURBN
TYPE(TLES_t),             INTENT(INOUT) :: TLES          ! modd_les structure
INTEGER,                  INTENT(IN)    ::  KSPLT        ! split process index
INTEGER,                  INTENT(IN)    ::  KRR          ! number of moist var.
INTEGER,                  INTENT(IN)    ::  KSV          ! number of sv var.
LOGICAL,                  INTENT(IN)    ::  OFLAT        ! Logical for zero ororography
LOGICAL,                  INTENT(IN)    ::  O2D          ! Logical for 2D model version (modd_conf)
TYPE(TFILEDATA),          INTENT(IN)    ::  TPFILE       ! Output file
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PK          ! Turbulent diffusion doef.
                                                        ! PK = PLM * SQRT(PTKEM)
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PINV_PDYY   ! 1./PDYY
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PINV_PDZZ   ! 1./PDZZ
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PMZM_PRHODJ ! MZM(PRHODJ)
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PDYY, PDZZ, PDZY 
                                                         ! Metric coefficients
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PRHODJ       ! density * grid volume
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PTHVREF      ! ref. state VPT       
!
! Variables at t-1
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PVM,PWM,PTHLM
REAL, DIMENSION(:,:,:,:), INTENT(IN)    ::  PRM
REAL, DIMENSION(:,:,:,:), INTENT(IN)    ::  PSVM
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PTKEM        ! TKE at time t- dt
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PLM          ! Turb. mixing length
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) ::  PRVS, PRWS   ! var. at t+1 -split-
REAL, DIMENSION(:,:,:),   INTENT(INOUT) ::  PDP          ! TKE production terms
!
!
!
!*       0.2  declaration of local variables
!
REAL, DIMENSION(SIZE(PWM,1),SIZE(PWM,2),SIZE(PWM,3))       &
                                     :: ZFLX,ZWORK
    ! work arrays
!   
!! REAL, DIMENSION(SIZE(PWM,1),SIZE(PWM,2),SIZE(PWM,3))  :: ZVPTV
INTEGER             :: IKB,IKE,IKU
                                    ! Index values for the Beginning and End
                                    ! mass points of the domain  
INTEGER             :: JSV          ! scalar loop counter
!
REAL, DIMENSION(SIZE(PWM,1),SIZE(PWM,2),SIZE(PWM,3))  :: GY_W_VW_PWM
!
REAL :: ZTIME1, ZTIME2
TYPE(TFIELDMETADATA) :: TZFIELD
! ---------------------------------------------------------------------------
!
!*       1.   PRELIMINARY COMPUTATIONS
!             ------------------------
!
IKB = 1+JPVEXT               
IKE = SIZE(PWM,3)-JPVEXT    
IKU = SIZE(PWM,3)
!
!
IF (.NOT. O2D) GY_W_VW_PWM = GY_W_VW(PWM,PDYY,PDZZ,PDZY)
!
!
!*      14.   < V'W'>
!             -------
!
! residual part of < V'W'> depending on dw/dy
!
IF (.NOT. O2D) THEN
  ZFLX(:,:,:) =                                                      &
    - XCMFS * MYM(MZM(PK)) * GY_W_VW_PWM
  !! &  to be tested
  !!  - (2./3.) * XCMFB * MZM( ZVPTV * MYM( PLM / SQRT(PTKEM) * XG / PTHVREF ) )
ELSE
  ZFLX(:,:,:) = 0.
  !! &  to be tested
  !!  - (2./3.) * XCMFB * MZM( ZVPTV * MYM( PLM / SQRT(PTKEM) * XG / PTHVREF ) )
END IF
!
ZFLX(:,:,IKE+1) = 0.  ! rigid wall condition => no turbulent flux
!
!
! Nullify the flux at the ground level because it has been fully taken into
! account in turb_ver and extrapolate the flux under the ground 
ZFLX(:,:,IKB) = 0.
ZFLX(:,:,IKB-1)= 2.*ZFLX(:,:,IKB) - ZFLX(:,:,IKB+1)
!
! stores  <V W>
IF ( TPFILE%LOPENED .AND. TURBN%LTURB_FLX ) THEN
  TZFIELD = TFIELDMETADATA(       &
    CMNHNAME   = 'VW_HFLX',       &
    CSTDNAME   = '',              &
    CLONGNAME  = 'VW_HFLX',       &
    CUNITS     = 'm2 s-2',        &
    CDIR       = 'XY',            &
    CCOMMENT   = 'X_Y_Z_VW_HFLX', &
    NGRID      = 7,               &
    NTYPE      = TYPEREAL,        &
    NDIMS      = 3,               &
    LTIMEDEP   = .TRUE.           )
  CALL IO_FIELD_WRITE(TPFILE,TZFIELD,ZFLX)
END IF
!
! compute the source for rho*V due to this residual flux ( the other part is
! taken into account in TURB_VER)
IF (.NOT. O2D) &
PRVS(:,:,:) = PRVS(:,:,:) - DZF( ZFLX* MYM( PMZM_PRHODJ ) / MYM ( PDZZ ) )
!
!computation of the source for rho*W due to this flux
IF (.NOT. O2D) THEN 
  IF (.NOT. OFLAT) THEN
    PRWS(:,:,:) = PRWS(:,:,:)                              &
          -DYF( MZM( MYM(PRHODJ) * PINV_PDYY) * ZFLX)           &
          +DZM( PRHODJ * MYF( MZF( ZFLX*PDZY ) * PINV_PDYY ) / MZF(PDZZ) )
  ELSE
    PRWS(:,:,:) = PRWS(:,:,:) - DYF( MZM( MYM(PRHODJ) * PINV_PDYY) * ZFLX)
  END IF
END IF
!
IF (KSPLT==1) THEN
  ! 
  !Contribution to the dynamic production of TKE:
  !
  IF (.NOT. O2D) THEN
    ZWORK(:,:,:) =-MZF( MYF( ZFLX *( GZ_V_VW(PVM,PDZZ) + GY_W_VW_PWM ) ) )
  !
  !
  ! evaluate the dynamic production at w(IKB+1) in PDP(IKB)
    ZWORK(:,:,IKB:IKB) = - MYF (                                               &
       ZFLX(:,:,IKB+1:IKB+1) *                                                 &
     (   (PVM(:,:,IKB+1:IKB+1)-PVM(:,:,IKB:IKB)) / MYM(PDZZ(:,:,IKB+1:IKB+1))  &
       + ( DYM( PWM(:,:,IKB+1:IKB+1) )                                         &
          -MYM(  (PWM(:,:,IKB+2:IKB+2)-PWM(:,:,IKB+1:IKB+1))                   &
                  /(PDZZ(:,:,IKB+2:IKB+2)+PDZZ(:,:,IKB+1:IKB+1))               &
                +(PWM(:,:,IKB+1:IKB+1)-PWM(:,:,IKB  :IKB  ))                   &
                  /(PDZZ(:,:,IKB+1:IKB+1)+PDZZ(:,:,IKB  :IKB  ))               &
              ) * PDZY(:,:,IKB+1:IKB+1)                                        &
         ) / (0.5*(PDYY(:,:,IKB+1:IKB+1)+PDYY(:,:,IKB:IKB)))                 &
     )                        )  
  ENDIF
  !
  ! dynamic production computation
  IF (.NOT. O2D) &
  PDP(:,:,:) = PDP(:,:,:) + ZWORK(:,:,:)  
  !
END IF
!
! Storage in the LES configuration (addition to TURB_VER computation)
!
IF (TLES%LLES_CALL .AND. KSPLT==1) THEN
  CALL SECOND_MNH(ZTIME1)
  CALL LES_MEAN_SUBGRID( MZF(MYF(ZFLX)), TLES%X_LES_SUBGRID_WV , .TRUE. )
  CALL LES_MEAN_SUBGRID( MZF(MYF(GZ_V_VW(PVM,PDZZ)*ZFLX)),&
                         TLES%X_LES_RES_ddxa_V_SBG_UaV , .TRUE.)
  CALL LES_MEAN_SUBGRID( MZF(MYF(GY_W_VW(PWM,PDYY,PDZZ,PDZY)*ZFLX)),&
                         TLES%X_LES_RES_ddxa_W_SBG_UaW , .TRUE.)
  CALL LES_MEAN_SUBGRID( MXF(GY_M_V(1,IKU,1,PTHLM,PDYY,PDZZ,PDZY)*MZF(ZFLX)),&
                         TLES%X_LES_RES_ddxa_Thl_SBG_UaW , .TRUE.)
  IF (KRR>=1) THEN
    CALL LES_MEAN_SUBGRID( MXF(GY_M_V(1,IKU,1,PRM(:,:,:,1),PDYY,PDZZ,PDZY)*MZF(ZFLX)), &
                           TLES%X_LES_RES_ddxa_Rt_SBG_UaW , .TRUE.)
  END IF
  DO JSV=1,KSV
    CALL LES_MEAN_SUBGRID( MXF(GY_M_V(1,IKU,1,PSVM(:,:,:,JSV),PDYY,PDZZ,PDZY)*MZF(ZFLX)), &
                           TLES%X_LES_RES_ddxa_Sv_SBG_UaW(:,:,:,JSV), .TRUE.)
  END DO
  CALL SECOND_MNH(ZTIME2)
  TLES%XTIME_LES = TLES%XTIME_LES + ZTIME2 - ZTIME1
END IF
!
!
!
END SUBROUTINE TURB_HOR_VW
END MODULE MODE_TURB_HOR_VW
