!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_TURB_HOR_SV_FLUX
IMPLICIT NONE
CONTAINS
      SUBROUTINE TURB_HOR_SV_FLUX(D,TURBN,TLES,KSPLT,OBLOWSNOW,OFLAT,&
                      TPFILE,KSV,KSV_LGBEG,KSV_LGEND,O2D,ONOMIXLG,   &
                      PK,PINV_PDXX,PINV_PDYY,PINV_PDZZ,PMZM_PRHODJ,  &
                      PDXX,PDYY,PDZZ,PDZX,PDZY,PRSNOW,               &
                      PDIRCOSXW,PDIRCOSYW,                           &
                      PRHODJ,PWM,                                    &
                      PSFSVM,                                        &
                      PSVM,                                          &
                      PRSVS                                          )
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
!!                     Oct  18, 2000 (V. Masson) LES computations + OFLAT swith
!!                                              + bug on Y scalar flux
!!                     Jun  20, 2001 (J Stein) case of lagragian variables
!!                     Nov  06, 2002 (V. Masson) LES budgets
!!                     October 2009 (G. Tanguy) add ILENCH=LEN(YCOMMENT) after
!!                                              change of YCOMMENT
!!                     04/2016 (M.Moge) Use openACC directives to port the TURB part of Meso-NH on GPU
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!! --------------------------------------------------------------------------
!
!*      0. DECLARATIONS
!          ------------
!
USE MODD_TURB_n, ONLY: TURB_t
USE MODD_DIMPHYEX,       ONLY: DIMPHYEX_t
!
USE MODD_CST
USE MODD_CTURB
use modd_field,          only: tfieldmetadata, TYPEREAL
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
USE MODE_GRADIENT_M_PHY, ONLY: GX_M_U_PHY, GY_M_V_PHY, GX_M_M_PHY, GY_M_M_PHY
USE MODE_GRADIENT_W_PHY, ONLY: GX_W_UW_PHY, GY_W_VW_PHY
USE MODI_SHUMAN
USE MODE_SHUMAN_PHY
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
INTEGER,                  INTENT(IN)    ::  KSPLT        ! split process index
TYPE(TFILEDATA),          INTENT(INOUT)    ::  TPFILE       ! Output file
INTEGER,                  INTENT(IN)    ::  KSV,KSV_LGBEG,KSV_LGEND ! number of sv var.
LOGICAL,                  INTENT(IN)    ::  OFLAT        ! Logical for zero ororography
LOGICAL,                  INTENT(IN)    ::  ONOMIXLG     ! to use turbulence for lagrangian variables (modd_conf)
LOGICAL,                  INTENT(IN)    ::  O2D          ! Logical for 2D model version (modd_conf)
LOGICAL,                  INTENT(IN)    ::  OBLOWSNOW    ! switch to activate pronostic blowing snow
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PK          ! Turbulent diffusion doef.
                                                        ! PK = PLM * SQRT(PTKEM)
REAL,                     INTENT(IN)    ::  PRSNOW       ! Ratio for diffusion coeff. scalar (blowing snow)
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PINV_PDXX   ! 1./PDXX
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PINV_PDYY   ! 1./PDYY
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PINV_PDZZ   ! 1./PDZZ
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PMZM_PRHODJ ! MZM(PRHODJ)
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PDXX, PDYY, PDZZ, PDZX, PDZY 
                                                         ! Metric coefficients
REAL, DIMENSION(:,:),     INTENT(IN)    ::  PDIRCOSXW, PDIRCOSYW
! Director Cosinus along x and y  directions at surface w-point
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PRHODJ       ! density * grid volume
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PWM          ! vertical wind
!
REAL, DIMENSION(D%NIT,D%NJT,KSV),   INTENT(IN)    ::  PSFSVM       ! surface fluxes
!
!
! Variables at t-1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT,KSV), INTENT(IN)    ::  PSVM         ! scalar var. at t-1
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT,KSV), INTENT(INOUT) ::  PRSVS        ! var. at t+1 -split-
!
!
!
!*       0.2  declaration of local variables
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT)       &
                                     :: ZFLXX,ZFLXY
    ! work arrays
REAL, DIMENSION(D%NIT,D%NJT) :: ZWORK2D
!
REAL :: ZCSV          !constant for the scalar flux

INTEGER             :: IKB,IKE
                                    ! Index values for the Beginning and End
                                    ! mass points of the domain  
INTEGER             :: JSV          ! loop counter
INTEGER             :: ISV          ! number of scalar var.
REAL, DIMENSION(SIZE(PDZZ,1),SIZE(PDZZ,2),1+JPVEXT:3+JPVEXT) :: ZCOEFF 
                                    ! coefficients for the uncentred gradient 
                                    ! computation near the ground
!
CHARACTER(LEN=NMNHNAMELGTMAX) :: YMNHNAME
INTEGER :: IKU, IIT, IJT, IKT, JI, JJ, JK
TYPE(TFIELDMETADATA) :: TZFIELD
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZMXM3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZGX_M_U3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT) ::ZMXM2D_WORK1
REAL, DIMENSION(D%NIT,D%NJT) ::ZDXM2D_WORK1
REAL, DIMENSION(D%NIT,D%NJT) ::ZSHUGRADWK1_2D
REAL, DIMENSION(D%NIT,D%NJT) ::ZMXM2D_WORK2
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZMYM3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZGY_M_V3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT) ::ZMYM2D_WORK1
REAL, DIMENSION(D%NIT,D%NJT) ::ZDYM2D_WORK1
REAL, DIMENSION(D%NIT,D%NJT) ::ZMYM2D_WORK2
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZSHUGRADWK1_3D
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZDXF3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZDYF3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZDZF3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZSHUGRADWK2_3D
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZMXF3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZSHUGRADWK3_3D
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZMZM3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZMYF3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZMZM3D_WORK2
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZMZF3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZGX_W_UW3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZGX_M_M3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZGY_W_VW3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZGY_M_M3D_WORK1
REAL :: ZTIME1, ZTIME2
! ---------------------------------------------------------------------------
!
!*       1.   PRELIMINARY COMPUTATIONS
!             ------------------------
!
IKB = 1+JPVEXT               
IKE = SIZE(PSVM,3)-JPVEXT   
IKU = SIZE(PSVM,3)
IIT=D%NIT
IJT=D%NJT
IKT=D%NKT
!
ISV = SIZE(PSVM,4)
!
IF(OBLOWSNOW) THEN
! See Vionnet (PhD, 2012) for a complete discussion around the value of the Schmidt number for blowing snow variables              
   ZCSV= TURBN%XCHF/PRSNOW
ELSE
   ZCSV= TURBN%XCHF
ENDIF
!
!  compute the coefficients for the uncentred gradient computation near the 
!  ground

DO JJ=1, IJT
  DO JI=1, IIT
    ZCOEFF(JI, JJ, IKB+2)= - PDZZ(JI, JJ, IKB+1) /      &
           ( (PDZZ(JI, JJ, IKB+2)+PDZZ(JI, JJ, IKB+1)) * PDZZ(JI, JJ, IKB+2) )
    ZCOEFF(JI, JJ, IKB+1)=   (PDZZ(JI, JJ, IKB+2)+PDZZ(JI, JJ, IKB+1)) /      &
           ( PDZZ(JI, JJ, IKB+1) * PDZZ(JI, JJ, IKB+2) )
    ZCOEFF(JI, JJ, IKB)= - (PDZZ(JI, JJ, IKB+2)+2.*PDZZ(JI, JJ, IKB+1)) /      &
           ( (PDZZ(JI, JJ, IKB+2)+PDZZ(JI, JJ, IKB+1)) * PDZZ(JI, JJ, IKB+1) )
  END DO
END DO

!
!
!*      15.   HORIZONTAL FLUXES OF PASSIVE SCALARS
!             ------------------------------------
!
!
DO JSV=1,ISV
!
  IF (ONOMIXLG .AND. JSV >= KSV_LGBEG .AND. JSV<= KSV_LGEND) CYCLE
!
!       15.1   <U' SVth'>
!              ----------
!
  ! Computes the flux in the X direction
  CALL MXM_PHY(D, PK(:, :, :), ZMXM3D_WORK1)
CALL GX_M_U_PHY(D, OFLAT,PSVM(:,:,:,JSV),PDXX(:, :, :),PDZZ(:, :, :),PDZX(:, :, :), ZGX_M_U3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZFLXX(JI, JJ, JK) = -ZCSV * ZMXM3D_WORK1(JI, JJ, JK) * ZGX_M_U3D_WORK1(JI, JJ, JK)
    END DO
  END DO
END DO

!

DO JJ=1, IJT
    DO JI=1, IIT
      ZFLXX(JI, JJ, IKE+1) = ZFLXX(JI, JJ, IKE)
  END DO
  END DO

!
! Compute the flux at the first inner U-point with an uncentred vertical  
! gradient
  CALL MXM2D_PHY(D, PK(:,:,IKB), ZMXM2D_WORK1)
CALL DXM2D_PHY(D, PSVM(:,:,IKB,JSV), ZDXM2D_WORK1)

DO JJ=1, IJT
  DO JI=1, IIT
    ZSHUGRADWK1_2D(JI, JJ) = ZCOEFF(JI, JJ, 2)*PSVM(JI, JJ, 2, JSV)        +ZCOEFF(JI, JJ, IKB+1)*PSVM(JI, JJ, IKB+1, JSV)        +ZCOEFF(JI, JJ, IKB)*PSVM(JI, JJ, IKB, JSV)
  END DO
END DO

!
CALL MXM2D_PHY(D, ZSHUGRADWK1_2D, ZMXM2D_WORK2)

DO JJ=1, IJT
  DO JI=1, IIT
    ZFLXX(JI, JJ, IKB) = -ZCSV * ZMXM2D_WORK1(JI, JJ) *             &
        ( ZDXM2D_WORK1(JI, JJ) * PINV_PDXX(JI, JJ, IKB)           &
         -ZMXM2D_WORK2(JI, JJ) * 0.5 * ( PDZX(JI, JJ, IKB+1)+PDZX(JI, JJ, IKB) )     &
                * PINV_PDXX(JI, JJ, IKB)                                &
        ) 
  END DO
END DO

!
! extrapolates the flux under the ground so that the vertical average with 
! the IKB flux gives the ground value

DO JJ=1, IJT
    DO JI=1, IIT
      ZWORK2D(JI, JJ)=PSFSVM(JI, JJ, JSV) * PDIRCOSXW(JI, JJ)
  END DO
  END DO

  CALL MXM2D_PHY(D, ZWORK2D(:,:), ZMXM2D_WORK1)

DO JJ=1, IJT
  DO JI=1, IIT
    ZFLXX(JI, JJ, IKB-1) = 2. * ZMXM2D_WORK1(JI, JJ) - ZFLXX(JI, JJ, IKB)  
  END DO
END DO

!
!
  ! stores  <U SVth>
  IF ( TPFILE%LOPENED .AND. TURBN%LTURB_FLX ) THEN
    WRITE(YMNHNAME,'("USV_FLX_",I3.3)') JSV
    TZFIELD = TFIELDMETADATA(                    &
      CMNHNAME   = TRIM( YMNHNAME ),             &
      CSTDNAME   = '',                           &
      CLONGNAME  = TRIM( YMNHNAME ),             &
      CUNITS     = 'SVUNIT m s-1',               &
      CDIR       = 'XY',                         &
      CCOMMENT   = 'X_Y_Z_' // TRIM( YMNHNAME ), &
      NGRID      = 2,                            &
      NTYPE      = TYPEREAL,                     &
      NDIMS      = 3,                            &
      LTIMEDEP   = .TRUE.                        )

    CALL IO_Field_write(TPFILE,TZFIELD,ZFLXX)
  END IF
!
  IF (TLES%LLES_CALL .AND. KSPLT==1) THEN
    CALL SECOND_MNH(ZTIME1)
    CALL MXF_PHY(D, ZFLXX(:, :, :), ZMXF3D_WORK1)
CALL LES_MEAN_SUBGRID( ZMXF3D_WORK1(:, :, :), TLES%X_LES_SUBGRID_USv(:,:,:,JSV) )     
CALL GX_W_UW_PHY(D, OFLAT,PWM(:, :, :),PDXX(:, :, :),PDZZ(:, :, :),PDZX(:, :, :), ZGX_W_UW3D_WORK1)
CALL MZM_PHY(D, ZFLXX(:, :, :), ZMZM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = ZGX_W_UW3D_WORK1(JI, JJ, JK)*ZMZM3D_WORK1(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MXF_PHY(D, ZSHUGRADWK2_3D, ZMXF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMXF3D_WORK1(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MZF_PHY(D, ZSHUGRADWK1_3D, ZMZF3D_WORK1)
CALL LES_MEAN_SUBGRID( ZMZF3D_WORK1(:, :, :), &
                           TLES%X_LES_RES_ddxa_W_SBG_UaSv(:,:,:,JSV) , .TRUE. )    
CALL GX_M_M_PHY(D, OFLAT,PSVM(:,:,:,JSV),PDXX(:, :, :),PDZZ(:, :, :),PDZX(:, :, :), ZGX_M_M3D_WORK1)
CALL MXF_PHY(D, ZFLXX(:, :, :), ZMXF3D_WORK1)
CALL LES_MEAN_SUBGRID( ZGX_M_M3D_WORK1(:, :, :)*ZMXF3D_WORK1(:, :, :), &
                           TLES%X_LES_RES_ddxa_Sv_SBG_UaSv(:,:,:,JSV), .TRUE. )    
CALL SECOND_MNH(ZTIME2)
    TLES%XTIME_LES = TLES%XTIME_LES + ZTIME2 - ZTIME1
  END IF
!
!       15.2   <V' SVth'>
!              ----------
!
  IF (.NOT. O2D) THEN
!
! Computes the flux in the Y direction
    CALL MYM_PHY(D, PK(:, :, :), ZMYM3D_WORK1)
CALL GY_M_V_PHY(D, OFLAT,PSVM(:,:,:,JSV),PDYY(:, :, :),PDZZ(:, :, :),PDZY(:, :, :), ZGY_M_V3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZFLXY(JI, JJ, JK)=-ZCSV * ZMYM3D_WORK1(JI, JJ, JK) * ZGY_M_V3D_WORK1(JI, JJ, JK)
    END DO
  END DO
END DO

!

DO JJ=1, IJT
      DO JI=1, IIT
        ZFLXY(JI, JJ, IKE+1) = ZFLXY(JI, JJ, IKE)
  END DO
    END DO

!
! Compute the flux at the first inner V-point with an uncentred vertical  
! gradient
!
    CALL MYM2D_PHY(D, PK(:,:,IKB), ZMYM2D_WORK1)
CALL DYM2D_PHY(D, PSVM(:,:,IKB,JSV), ZDYM2D_WORK1)

DO JJ=1, IJT
  DO JI=1, IIT
    ZSHUGRADWK1_2D(JI, JJ) = ZCOEFF(JI, JJ, 2)*PSVM(JI, JJ, 2, JSV)        +ZCOEFF(JI, JJ, IKB+1)*PSVM(JI, JJ, IKB+1, JSV)        +ZCOEFF(JI, JJ, IKB)*PSVM(JI, JJ, IKB, JSV)
  END DO
END DO

!
CALL MYM2D_PHY(D, ZSHUGRADWK1_2D, ZMYM2D_WORK2)

DO JJ=1, IJT
  DO JI=1, IIT
    ZFLXY(JI, JJ, IKB) = -ZCSV * ZMYM2D_WORK1(JI, JJ) *             &
          ( ZDYM2D_WORK1(JI, JJ) * PINV_PDYY(JI, JJ, IKB)           &
           -ZMYM2D_WORK2(JI, JJ) * 0.5 * ( PDZY(JI, JJ, IKB+1)+PDZY(JI, JJ, IKB) )     &
                  * PINV_PDYY(JI, JJ, IKB)                                &
          ) 
  END DO
END DO

!
! extrapolates the flux under the ground so that the vertical average with 
! the IKB flux gives the ground value

DO JJ=1, IJT
      DO JI=1, IIT
        ZWORK2D(JI, JJ)=PSFSVM(JI, JJ, JSV) * PDIRCOSYW(JI, JJ)
  END DO
    END DO

    CALL MYM2D_PHY(D, ZWORK2D(:,:), ZMYM2D_WORK1)

DO JJ=1, IJT
  DO JI=1, IIT
    ZFLXY(JI, JJ, IKB-1) = 2. * ZMYM2D_WORK1(JI, JJ) - ZFLXY(JI, JJ, IKB)  
  END DO
END DO

!
!
  ! stores  <V SVth>
    IF ( TPFILE%LOPENED .AND. TURBN%LTURB_FLX ) THEN
      WRITE(YMNHNAME,'("VSV_FLX_",I3.3)') JSV
    TZFIELD = TFIELDMETADATA(                        &
      CMNHNAME   = TRIM( YMNHNAME ),                 &
      CSTDNAME   = '',                               &
      CLONGNAME  = TRIM(TZFIELD%CMNHNAME),           &
      CUNITS     = 'SVUNIT m s-1',                   &
      CDIR       = 'XY',                             &
      CCOMMENT   = 'X_Y_Z_'//TRIM(TZFIELD%CMNHNAME), &
      NGRID      = 3,                                &
      NTYPE      = TYPEREAL,                         &
      NDIMS      = 3,                                &
      LTIMEDEP   = .TRUE.                            )

      CALL IO_Field_write(TPFILE,TZFIELD,ZFLXY)
    END IF
!
  ELSE
    
    DO JK=1, IKT
      DO JJ=1, IJT
        DO JI=1, IIT
          ZFLXY(JI, JJ, JK)=0.
        END DO
      END DO
    END DO
    
  END IF
!
  IF (TLES%LLES_CALL .AND. KSPLT==1) THEN
    CALL SECOND_MNH(ZTIME1)
    CALL MYF_PHY(D, ZFLXY(:, :, :), ZMYF3D_WORK1)
CALL LES_MEAN_SUBGRID( ZMYF3D_WORK1(:, :, :), TLES%X_LES_SUBGRID_VSv(:,:,:,JSV) )     
CALL GY_W_VW_PHY(D, OFLAT,PWM(:, :, :),PDYY(:, :, :),PDZZ(:, :, :),PDZY(:, :, :), ZGY_W_VW3D_WORK1)
CALL MZM_PHY(D, ZFLXY(:, :, :), ZMZM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = ZGY_W_VW3D_WORK1(JI, JJ, JK)*ZMZM3D_WORK1(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MYF_PHY(D, ZSHUGRADWK2_3D, ZMYF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMYF3D_WORK1(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MZF_PHY(D, ZSHUGRADWK1_3D, ZMZF3D_WORK1)
CALL LES_MEAN_SUBGRID( ZMZF3D_WORK1(:, :, :), &
                           TLES%X_LES_RES_ddxa_W_SBG_UaSv(:,:,:,JSV) , .TRUE. )    
CALL GY_M_M_PHY(D, OFLAT,PSVM(:,:,:,JSV),PDYY(:, :, :),PDZZ(:, :, :),PDZY(:, :, :), ZGY_M_M3D_WORK1)
CALL MYF_PHY(D, ZFLXY(:, :, :), ZMYF3D_WORK1)
CALL LES_MEAN_SUBGRID( ZGY_M_M3D_WORK1(:, :, :)*ZMYF3D_WORK1(:, :, :), &
                           TLES%X_LES_RES_ddxa_Sv_SBG_UaSv(:,:,:,JSV) , .TRUE. )    
CALL SECOND_MNH(ZTIME2)
    TLES%XTIME_LES = TLES%XTIME_LES + ZTIME2 - ZTIME1
  END IF
!
!
!       15.3   Horizontal source terms
!              -----------------------
!
  IF (.NOT. O2D) THEN
    IF (.NOT. OFLAT) THEN
      CALL MXM_PHY(D, PRHODJ(:, :, :), ZMXM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMXM3D_WORK1(JI, JJ, JK) * ZFLXX(JI, JJ, JK) * PINV_PDXX(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL DXF_PHY(D, ZSHUGRADWK1_3D, ZDXF3D_WORK1)
CALL MYM_PHY(D, PRHODJ(:, :, :), ZMYM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMYM3D_WORK1(JI, JJ, JK) * ZFLXY(JI, JJ, JK) * PINV_PDYY(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL DYF_PHY(D, ZSHUGRADWK1_3D, ZDYF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK3_3D(JI, JJ, JK) = ZFLXX(JI, JJ, JK) * PINV_PDXX(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MZM_PHY(D, ZSHUGRADWK3_3D, ZMZM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = ZMZM3D_WORK1(JI, JJ, JK) * PDZX(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MXF_PHY(D, ZSHUGRADWK2_3D, ZMXF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = ZFLXY(JI, JJ, JK) * PINV_PDYY(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MZM_PHY(D, ZSHUGRADWK2_3D, ZMZM3D_WORK2)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMZM3D_WORK2(JI, JJ, JK) * PDZY(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MYF_PHY(D, ZSHUGRADWK1_3D, ZMYF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = PMZM_PRHODJ(JI, JJ, JK) * PINV_PDZZ(JI, JJ, JK) *                                              ( ZMXF3D_WORK1(JI, JJ, JK) + ZMYF3D_WORK1(JI, JJ, JK) )
    END DO
  END DO
END DO

!
CALL DZF_PHY(D, ZSHUGRADWK1_3D, ZDZF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      PRSVS(JI, JJ, JK, JSV)=   PRSVS(JI, JJ, JK, JSV)                                          &
              -ZDXF3D_WORK1(JI, JJ, JK)                                    &
              -ZDYF3D_WORK1(JI, JJ, JK)                                    &
              +ZDZF3D_WORK1(JI, JJ, JK)    
    END DO
  END DO
END DO

!
ELSE
      CALL MXM_PHY(D, PRHODJ(:, :, :), ZMXM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMXM3D_WORK1(JI, JJ, JK) * ZFLXX(JI, JJ, JK) * PINV_PDXX(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL DXF_PHY(D, ZSHUGRADWK1_3D, ZDXF3D_WORK1)
CALL MYM_PHY(D, PRHODJ(:, :, :), ZMYM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMYM3D_WORK1(JI, JJ, JK) * ZFLXY(JI, JJ, JK) * PINV_PDYY(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL DYF_PHY(D, ZSHUGRADWK1_3D, ZDYF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      PRSVS(JI, JJ, JK, JSV)=   PRSVS(JI, JJ, JK, JSV)                                          &
              -ZDXF3D_WORK1(JI, JJ, JK)                                    &
              -ZDYF3D_WORK1(JI, JJ, JK)    
    END DO
  END DO
END DO

!
END IF
  ELSE
    IF (.NOT. OFLAT) THEN
      CALL MXM_PHY(D, PRHODJ(:, :, :), ZMXM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMXM3D_WORK1(JI, JJ, JK) * ZFLXX(JI, JJ, JK) * PINV_PDXX(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL DXF_PHY(D, ZSHUGRADWK1_3D, ZDXF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK3_3D(JI, JJ, JK) = ZFLXX(JI, JJ, JK) * PINV_PDXX(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MZM_PHY(D, ZSHUGRADWK3_3D, ZMZM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = ZMZM3D_WORK1(JI, JJ, JK) * PDZX(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MXF_PHY(D, ZSHUGRADWK2_3D, ZMXF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = PMZM_PRHODJ(JI, JJ, JK) * PINV_PDZZ(JI, JJ, JK) *                                              ( ZMXF3D_WORK1(JI, JJ, JK) )
    END DO
  END DO
END DO

!
CALL DZF_PHY(D, ZSHUGRADWK1_3D, ZDZF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      PRSVS(JI, JJ, JK, JSV)=   PRSVS(JI, JJ, JK, JSV)                                          &
              -ZDXF3D_WORK1(JI, JJ, JK)                                    &
              +ZDZF3D_WORK1(JI, JJ, JK)    
    END DO
  END DO
END DO

!
ELSE
      CALL MXM_PHY(D, PRHODJ(:, :, :), ZMXM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMXM3D_WORK1(JI, JJ, JK) *ZFLXX(JI, JJ, JK) * PINV_PDXX(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL DXF_PHY(D, ZSHUGRADWK1_3D, ZDXF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      PRSVS(JI, JJ, JK, JSV)=   PRSVS(JI, JJ, JK, JSV)                                          &
              -ZDXF3D_WORK1(JI, JJ, JK)    
    END DO
  END DO
END DO

!
END IF
  END IF
!
!
END DO    ! end loop JSV
!
!
END SUBROUTINE TURB_HOR_SV_FLUX
END MODULE MODE_TURB_HOR_SV_FLUX
