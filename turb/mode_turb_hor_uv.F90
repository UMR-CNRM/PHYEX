!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_TURB_HOR_UV
IMPLICIT NONE
CONTAINS
!     ################################################################
      SUBROUTINE TURB_HOR_UV(D, TURBN,TLES,KSPLT,OFLAT,O2D,          &
                      TPFILE,                                        &
                      PK,PINV_PDXX,PINV_PDYY,PINV_PDZZ,PMZM_PRHODJ,  &
                      PDXX,PDYY,PDZZ,PDZX,PDZY,                      &
                      PDIRCOSZW,                                     &
                      PCOSSLOPE,PSINSLOPE,                           &
                      PRHODJ,                                        &
                      PCDUEFF,PTAU11M,PTAU12M,PTAU22M,PTAU33M,       &
                      PUM,PVM,PUSLOPEM,PVSLOPEM,                     &
                      PDP,                                           &
                      PRUS,PRVS                                      )
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
!!                     Nov  06, 2002 (V. Masson) LES budgets
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
USE MODE_SHUMAN_PHY
USE MODI_LES_MEAN_SUBGRID
!
USE MODI_SECOND_MNH
!
#ifdef MNH_COMPILER_CCE
!$mnh_undef(LOOP)
!$mnh_undef(OPENACC)
#endif
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
TYPE(TFILEDATA),          INTENT(INOUT) ::  TPFILE       ! Output file
LOGICAL,                  INTENT(IN)    ::  OFLAT        ! Logical for zero ororography
LOGICAL,                  INTENT(IN)    ::  O2D          ! Logical for 2D model version (modd_conf)
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PK          ! Turbulent diffusion doef.
                                                        ! PK = PLM * SQRT(PTKEM)
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PINV_PDXX   ! 1./PDXX
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PINV_PDYY   ! 1./PDYY
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PINV_PDZZ   ! 1./PDZZ
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PMZM_PRHODJ ! MZM(PRHODJ)
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PDXX, PDYY, PDZZ, PDZX, PDZY 
                                                         ! Metric coefficients
REAL, DIMENSION(D%NIT,D%NJT),     INTENT(IN)    ::  PDIRCOSZW
! Director Cosinus along z directions at surface w-point
REAL, DIMENSION(D%NIT,D%NJT),   INTENT(IN)   ::  PCOSSLOPE       ! cosinus of the angle 
                                      ! between i and the slope vector
REAL, DIMENSION(D%NIT,D%NJT),   INTENT(IN)   ::  PSINSLOPE       ! sinus of the angle 
                                      ! between i and the slope vector
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PRHODJ       ! density * grid volume
!
REAL, DIMENSION(D%NIT,D%NJT),   INTENT(IN)   ::  PCDUEFF      ! Cd * || u || at time t
REAL, DIMENSION(D%NIT,D%NJT),   INTENT(IN)   ::  PTAU11M      ! <uu> in the axes linked 
       ! to the maximum slope direction and the surface normal and the binormal 
       ! at time t - dt
REAL, DIMENSION(D%NIT,D%NJT),   INTENT(IN)   ::  PTAU12M      ! <uv> in the same axes
REAL, DIMENSION(D%NIT,D%NJT),   INTENT(IN)   ::  PTAU22M      ! <vv> in the same axes
REAL, DIMENSION(D%NIT,D%NJT),   INTENT(IN)   ::  PTAU33M      ! <ww> in the same axes
!
! Variables at t-1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PUM,PVM
REAL, DIMENSION(D%NIT,D%NJT),      INTENT(IN)   ::  PUSLOPEM     ! wind component along the 
                                     ! maximum slope direction
REAL, DIMENSION(D%NIT,D%NJT),      INTENT(IN)   ::  PVSLOPEM     ! wind component along the 
                                     ! direction normal to the maximum slope one
!
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(INOUT) ::  PRUS, PRVS   ! var. at t+1 -split-
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(INOUT) ::  PDP          ! TKE production terms
!
!
!
!*       0.2  declaration of local variables
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: ZFLX,ZWORK ! work arrays
!   
REAL, DIMENSION(D%NIT,D%NJT) :: ZDIRSINZW 
      ! sinus of the angle between the vertical and the normal to the orography
INTEGER             :: IKB,IKE,IIT,IJT,IKT, JI,JJ,JK
                                    ! Index values for the Beginning and End
                                    ! mass points of the domain  
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT)  :: GY_U_UV_PUM
REAL, DIMENSION(D%NIT,D%NJT,D%NKT)  :: GX_V_UV_PVM
!
REAL :: ZTIME1, ZTIME2
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZGX_V_UV3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZGY_U_UV3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZMYM3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZMXM3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT) ::ZMYM2D_WORK1
REAL, DIMENSION(D%NIT,D%NJT) ::ZMXM2D_WORK1
REAL, DIMENSION(D%NIT,D%NJT) ::ZDYM2D_WORK1
REAL, DIMENSION(D%NIT,D%NJT) ::ZSHUGRADWK1_2D
REAL, DIMENSION(D%NIT,D%NJT) ::ZMYM2D_WORK2
REAL, DIMENSION(D%NIT,D%NJT) ::ZSHUGRADWK2_2D
REAL, DIMENSION(D%NIT,D%NJT) ::ZMXM2D_WORK2
REAL, DIMENSION(D%NIT,D%NJT) ::ZMXM2D_WORK3
REAL, DIMENSION(D%NIT,D%NJT) ::ZMXM2D_WORK4
REAL, DIMENSION(D%NIT,D%NJT) ::ZMXM2D_WORK5
REAL, DIMENSION(D%NIT,D%NJT) ::ZDXM2D_WORK1
REAL, DIMENSION(D%NIT,D%NJT) ::ZMXM2D_WORK6
REAL, DIMENSION(D%NIT,D%NJT) ::ZMYM2D_WORK3
REAL, DIMENSION(D%NIT,D%NJT) ::ZMYM2D_WORK4
REAL, DIMENSION(D%NIT,D%NJT) ::ZMYM2D_WORK5
REAL, DIMENSION(D%NIT,D%NJT) ::ZMYM2D_WORK6
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZSHUGRADWK1_3D
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZDYF3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZSHUGRADWK2_3D
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZDZF3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZMYF3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZMZM3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZMXM3D_WORK2
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZMXM3D_WORK3
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZMZM3D_WORK2
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZDXF3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZMXF3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZMYM3D_WORK2
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZMYM3D_WORK3
REAL, DIMENSION(D%NIT,D%NJT) ::ZMXF2D_WORK1
REAL, DIMENSION(D%NIT,D%NJT) ::ZMYF2D_WORK1
REAL, DIMENSION(D%NIT,D%NJT) ::ZMXF2D_WORK2
REAL, DIMENSION(D%NIT,D%NJT) ::ZMYF2D_WORK2
REAL, DIMENSION(D%NIT,D%NJT) ::ZSHUGRADWK3_2D
REAL, DIMENSION(D%NIT,D%NJT) ::ZMXF2D_WORK3
REAL, DIMENSION(D%NIT,D%NJT) ::ZMYF2D_WORK3
REAL, DIMENSION(D%NIT,D%NJT) ::ZMYF2D_WORK4
REAL, DIMENSION(D%NIT,D%NJT) ::ZMXF2D_WORK4
TYPE(TFIELDMETADATA) :: TZFIELD
! ---------------------------------------------------------------------------
!
!*       1.   PRELIMINARY COMPUTATIONS
!             ------------------------
!
IKB = 1+JPVEXT               
IKE = SIZE(PUM,3)-JPVEXT
IIT=D%NIT
IJT=D%NJT
IKT=D%NKT  
!

DO JJ=1, IJT
  DO JI=1, IIT
    ZDIRSINZW(JI, JJ) = SQRT( 1. - PDIRCOSZW(JI, JJ)**2 )
  END DO
END DO

!
CALL GX_V_UV_DEVICE(PVM(:, :, :),PDXX(:, :, :),PDZZ(:, :, :),PDZX(:, :, :), ZGX_V_UV3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
GX_V_UV_PVM(:, :, :) = ZGX_V_UV3D_WORK1(:, :, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
IF (.NOT. O2D) THEN
  CALL GY_U_UV_DEVICE(PUM(:, :, :),PDYY(:, :, :),PDZZ(:, :, :),PDZY(:, :, :), ZGY_U_UV3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
GY_U_UV_PUM(:, :, :) = ZGY_U_UV3D_WORK1(:, :, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
END IF
!
!
!*      12.   < U'V'>
!             -------
!
!
IF (.NOT. O2D) THEN
  CALL MXM_PHY(D, PK(:, :, :), ZMXM3D_WORK1)
CALL MYM_PHY(D, ZMXM3D_WORK1, ZMYM3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZFLX(:,:,:)= - XCMFS * ZMYM3D_WORK1(:, :, :) *                           &
          (GY_U_UV_PUM(:, :, :) + GX_V_UV_PVM(:, :, :))
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
ELSE
  CALL MXM_PHY(D, PK(:, :, :), ZMXM3D_WORK1)
CALL MYM_PHY(D, ZMXM3D_WORK1, ZMYM3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZFLX(:,:,:)= - XCMFS * ZMYM3D_WORK1(:, :, :) *                           &
          (GX_V_UV_PVM(:, :, :))
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
END IF
!

DO JJ=1, IJT
  DO JI=1, IIT
    ZFLX(JI, JJ, IKE+1)= ZFLX(JI, JJ, IKE)
  END DO
END DO

!
!
! Compute the correlation at the first physical level with the following 
! hypothesis du/dz vary in 1/z and w=0 at the ground
CALL MXM2D_PHY(D, PK(:,:,IKB), ZMXM2D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
ZSHUGRADWK2_2D(:, :) = PDZZ(:,:,IKB+1)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)

!
CALL MXM2D_PHY(D, ZSHUGRADWK2_2D, ZMXM2D_WORK2)
CALL MYM2D_PHY(D, ZMXM2D_WORK1, ZMYM2D_WORK1)
CALL MXM2D_PHY(D, PDZZ(:,:,IKB), ZMXM2D_WORK4)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
ZSHUGRADWK1_2D(:, :) = (PUM(:,:,IKB+1)-PUM(:,:,IKB))                   &
        *(1./ZMXM2D_WORK2(:, :)+1./ZMXM2D_WORK4(:, :))
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)

!
CALL MYM2D_PHY(D, ZSHUGRADWK1_2D, ZMYM2D_WORK2)
CALL DYM2D_PHY(D, PUM(:,:,IKB), ZDYM2D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
ZSHUGRADWK1_2D(:, :) = (PDZY(:,:,IKB+1)+PDZY(:,:,IKB))
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)

!
CALL MXM2D_PHY(D, ZSHUGRADWK1_2D, ZMXM2D_WORK3)
CALL MXM2D_PHY(D, PDYY(:,:,IKB), ZMXM2D_WORK5)
CALL DXM2D_PHY(D, PVM(:,:,IKB), ZDXM2D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
ZSHUGRADWK2_2D(:, :) = PDZZ(:,:,IKB+1)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)

!
CALL MYM2D_PHY(D, ZSHUGRADWK2_2D, ZMYM2D_WORK3)
CALL MYM2D_PHY(D, PDZZ(:,:,IKB), ZMYM2D_WORK5)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
ZSHUGRADWK1_2D(:, :) = (PVM(:,:,IKB+1)-PVM(:,:,IKB))                   &
        *(1./ZMYM2D_WORK3(:, :)+1./ZMYM2D_WORK5(:, :))
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)

!
CALL MXM2D_PHY(D, ZSHUGRADWK1_2D, ZMXM2D_WORK6)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
ZSHUGRADWK1_2D(:, :) = (PDZX(:,:,IKB+1)+PDZX(:,:,IKB))
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)

!
CALL MYM2D_PHY(D, ZSHUGRADWK1_2D, ZMYM2D_WORK4)
CALL MYM2D_PHY(D, PDXX(:,:,IKB), ZMYM2D_WORK6)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
ZFLX(:,:,IKB)   = - XCMFS * ZMYM2D_WORK1(:, :) *  (     &
  ( ZDYM2D_WORK1(:, :)                                        &
   -ZMYM2D_WORK2(:, :)&
    *0.5*ZMXM2D_WORK3(:, :)            &
  ) / ZMXM2D_WORK5(:, :)                                       &
 +( ZDXM2D_WORK1(:, :)                                        &
   -ZMXM2D_WORK6(:, :)&
    *0.5*ZMYM2D_WORK4(:, :)            &
  ) / ZMYM2D_WORK6(:, :)                                 ) 
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)

!
! 
! extrapolates this flux under the ground with the surface flux

DO JJ=1, IJT
  DO JI=1, IIT
    ZFLX(JI, JJ, IKB-1) =                                                           &
       PTAU11M(JI, JJ) * PCOSSLOPE(JI, JJ) * PSINSLOPE(JI, JJ) * PDIRCOSZW(JI, JJ)**2         &
      +PTAU12M(JI, JJ) * (PCOSSLOPE(JI, JJ)**2 - PSINSLOPE(JI, JJ)**2) *                   &
                      PDIRCOSZW(JI, JJ)**2                                           &
      -PTAU22M(JI, JJ) * PCOSSLOPE(JI, JJ) * PSINSLOPE(JI, JJ)                             &
      +PTAU33M(JI, JJ) * PCOSSLOPE(JI, JJ) * PSINSLOPE(JI, JJ) * ZDIRSINZW(JI, JJ)**2         &
      -PCDUEFF(JI, JJ) * (                                                           &
        2. * PUSLOPEM(JI, JJ) * PCOSSLOPE(JI, JJ) * PSINSLOPE(JI, JJ) *                    &
              PDIRCOSZW(JI, JJ) * ZDIRSINZW(JI, JJ)                                     &
        +PVSLOPEM(JI, JJ) * (PCOSSLOPE(JI, JJ)**2 - PSINSLOPE(JI, JJ)**2) * ZDIRSINZW(JI, JJ) &
                       )
    !
  END DO
END DO

!

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
ZSHUGRADWK2_2D(:, :) = ZFLX(:,:,IKB-1)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)

!
CALL MYM2D_PHY(D, ZSHUGRADWK2_2D, ZMYM2D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
ZSHUGRADWK1_2D(:, :) = ZMYM2D_WORK1(:, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)

!
CALL MXM2D_PHY(D, ZSHUGRADWK1_2D, ZMXM2D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
ZFLX(:,:,IKB-1) = 2. * ZMXM2D_WORK1(:, :)  &
                   - ZFLX(:,:,IKB)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)

!
!     
! stores  <U V>
IF ( TPFILE%LOPENED .AND. TURBN%LTURB_FLX ) THEN
  TZFIELD = TFIELDMETADATA(      &
    CMNHNAME   = 'UV_FLX',       &
    CSTDNAME   = '',             &
    CLONGNAME  = 'UV_FLX',       &
    CUNITS     = 'm2 s-2',       &
    CDIR       = 'XY',           &
    CCOMMENT   = 'X_Y_Z_UV_FLX', &
    NGRID      = 5,              &
    NTYPE      = TYPEREAL,       &
    NDIMS      = 3,              &
    LTIMEDEP   = .TRUE.          )

  CALL IO_FIELD_WRITE(TPFILE,TZFIELD,ZFLX)
END IF
!
!
!
!computation of the source for rho*V due to this flux
IF (.NOT. OFLAT) THEN
  CALL MYM_PHY(D, PRHODJ(:, :, :), ZMYM3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZSHUGRADWK2_3D(:, :, :) = ZMYM3D_WORK1(:, :, :) * PINV_PDYY(:, :, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
CALL MXM_PHY(D, ZSHUGRADWK2_3D, ZMXM3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZSHUGRADWK1_3D(:, :, :) = ZFLX(:, :, :) * ZMXM3D_WORK1(:, :, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
CALL DYF_PHY(D, ZSHUGRADWK1_3D, ZDYF3D_WORK1)
CALL MZM_PHY(D, ZFLX(:, :, :), ZMZM3D_WORK1)
CALL MZM_PHY(D, PDYY(:, :, :), ZMZM3D_WORK2)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZSHUGRADWK2_3D(:, :, :) = PDZY(:, :, :)/ZMZM3D_WORK2(:, :, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
CALL MXM_PHY(D, ZSHUGRADWK2_3D, ZMXM3D_WORK3)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZSHUGRADWK2_3D(:, :, :) = ZMZM3D_WORK1(:, :, :)*ZMXM3D_WORK3(:, :, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
CALL MYF_PHY(D, ZSHUGRADWK2_3D, ZMYF3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZSHUGRADWK1_3D(:, :, :) = PMZM_PRHODJ(:, :, :) * PINV_PDZZ(:, :, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
CALL MXM_PHY(D, ZSHUGRADWK1_3D, ZMXM3D_WORK2)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZSHUGRADWK1_3D(:, :, :) = ZMYF3D_WORK1(:, :, :)   &
                    * ZMXM3D_WORK2(:, :, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
CALL DZF_PHY(D, ZSHUGRADWK1_3D, ZDZF3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
PRUS(:,:,:) = PRUS(:,:,:)                                &
              - ZDYF3D_WORK1(:, :, :)         &
              + ZDZF3D_WORK1(:, :, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
ELSE
  CALL MYM_PHY(D, PRHODJ(:, :, :), ZMYM3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZSHUGRADWK2_3D(:, :, :) = ZMYM3D_WORK1(:, :, :) * PINV_PDYY(:, :, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
CALL MXM_PHY(D, ZSHUGRADWK2_3D, ZMXM3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZSHUGRADWK1_3D(:, :, :) = ZFLX(:, :, :) * ZMXM3D_WORK1(:, :, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
CALL DYF_PHY(D, ZSHUGRADWK1_3D, ZDYF3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
PRUS(:,:,:) = PRUS(:,:,:) - ZDYF3D_WORK1(:, :, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
END IF
!
!computation of the source for rho*V due to this flux
IF (.NOT. OFLAT) THEN
  CALL MXM_PHY(D, PRHODJ(:, :, :), ZMXM3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZSHUGRADWK2_3D(:, :, :) = ZMXM3D_WORK1(:, :, :) * PINV_PDXX(:, :, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
CALL MYM_PHY(D, ZSHUGRADWK2_3D, ZMYM3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZSHUGRADWK1_3D(:, :, :) = ZFLX(:, :, :) * ZMYM3D_WORK1(:, :, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
CALL DXF_PHY(D, ZSHUGRADWK1_3D, ZDXF3D_WORK1)
CALL MZM_PHY(D, ZFLX(:, :, :), ZMZM3D_WORK1)
CALL MZM_PHY(D, PDXX(:, :, :), ZMZM3D_WORK2)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZSHUGRADWK2_3D(:, :, :) = PDZX(:, :, :)/ZMZM3D_WORK2(:, :, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
CALL MYM_PHY(D, ZSHUGRADWK2_3D, ZMYM3D_WORK3)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZSHUGRADWK2_3D(:, :, :) = ZMZM3D_WORK1(:, :, :)*ZMYM3D_WORK3(:, :, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
CALL MXF_PHY(D, ZSHUGRADWK2_3D, ZMXF3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZSHUGRADWK1_3D(:, :, :) = PMZM_PRHODJ(:, :, :) * PINV_PDZZ(:, :, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
CALL MYM_PHY(D, ZSHUGRADWK1_3D, ZMYM3D_WORK2)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZSHUGRADWK1_3D(:, :, :) = ZMXF3D_WORK1(:, :, :) &
                      * ZMYM3D_WORK2(:, :, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
CALL DZF_PHY(D, ZSHUGRADWK1_3D, ZDZF3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
PRVS(:,:,:) = PRVS(:,:,:)                             &
                - ZDXF3D_WORK1(:, :, :)    &
                + ZDZF3D_WORK1(:, :, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
ELSE
  CALL MXM_PHY(D, PRHODJ(:, :, :), ZMXM3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZSHUGRADWK2_3D(:, :, :) = ZMXM3D_WORK1(:, :, :) * PINV_PDXX(:, :, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
CALL MYM_PHY(D, ZSHUGRADWK2_3D, ZMYM3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZSHUGRADWK1_3D(:, :, :) = ZFLX(:, :, :) * ZMYM3D_WORK1(:, :, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
CALL DXF_PHY(D, ZSHUGRADWK1_3D, ZDXF3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
PRVS(:,:,:) = PRVS(:,:,:) - ZDXF3D_WORK1(:, :, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
END IF
!
IF (KSPLT==1) THEN
  !
  !Contribution to the dynamic production of TKE:
  !
  IF (.NOT. O2D) THEN
    
!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZSHUGRADWK2_3D(:, :, :) = ZFLX(:, :, :) *                                &
       (GY_U_UV_PUM(:, :, :) + GX_V_UV_PVM(:, :, :))
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
CALL MYF_PHY(D, ZSHUGRADWK2_3D, ZMYF3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZSHUGRADWK1_3D(:, :, :) = ZMYF3D_WORK1(:, :, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
CALL MXF_PHY(D, ZSHUGRADWK1_3D, ZMXF3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZWORK(:,:,:) = - ZMXF3D_WORK1(:, :, :)   
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
ELSE
    
!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZSHUGRADWK2_3D(:, :, :) = ZFLX(:, :, :) *                                &
       (GX_V_UV_PVM(:, :, :))
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
CALL MYF_PHY(D, ZSHUGRADWK2_3D, ZMYF3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZSHUGRADWK1_3D(:, :, :) = ZMYF3D_WORK1(:, :, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
CALL MXF_PHY(D, ZSHUGRADWK1_3D, ZMXF3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZWORK(:,:,:) = - ZMXF3D_WORK1(:, :, :)    
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
ENDIF
  !
  ! evaluate the dynamic production at w(IKB+1) in PDP(IKB)
  !
  
!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
ZSHUGRADWK2_2D(:, :) = 0.5 * (ZFLX(:,:,IKB+1)+ZFLX(:,:,IKB))
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)

!
CALL MYF2D_PHY(D, ZSHUGRADWK2_2D, ZMYF2D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
ZSHUGRADWK1_2D(:, :) = ZMYF2D_WORK1(:, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)

!
CALL MXF2D_PHY(D, ZSHUGRADWK1_2D, ZMXF2D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
ZSHUGRADWK3_2D(:, :) = 0.5 * (PUM(:,:,IKB+1)+PUM(:,:,IKB))
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)

!
CALL DYM2D_PHY(D, ZSHUGRADWK3_2D, ZDYM2D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
ZSHUGRADWK2_2D(:, :) = 0.5*(PDYY(:,:,IKB)+PDYY(:,:,IKB+1))
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)

!
CALL MXM2D_PHY(D, ZSHUGRADWK2_2D, ZMXM2D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
ZSHUGRADWK2_2D(:, :) = 0.5 * (PVM(:,:,IKB+1)+PVM(:,:,IKB))
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)

!
CALL DXM2D_PHY(D, ZSHUGRADWK2_2D, ZDXM2D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
ZSHUGRADWK2_2D(:, :) = 0.5*(PDXX(:,:,IKB)+PDXX(:,:,IKB+1))
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)

!
CALL MYM2D_PHY(D, ZSHUGRADWK2_2D, ZMYM2D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
ZSHUGRADWK2_2D(:, :) = ZDYM2D_WORK1(:, :)           / ZMXM2D_WORK1(:, :)           +ZDXM2D_WORK1(:, :)           / ZMYM2D_WORK1(:, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)

!
CALL MYF2D_PHY(D, ZSHUGRADWK2_2D, ZMYF2D_WORK2)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
ZSHUGRADWK1_2D(:, :) = ZMYF2D_WORK2(:, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)

!
CALL MXF2D_PHY(D, ZSHUGRADWK1_2D, ZMXF2D_WORK2)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
ZSHUGRADWK3_2D(:, :) = PDZZ(:,:,IKB+1)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)

!
CALL MXM2D_PHY(D, ZSHUGRADWK3_2D, ZMXM2D_WORK2)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
ZSHUGRADWK1_2D(:, :) = (PUM(:,:,IKB+1)-PUM(:,:,IKB)) /                    ZMXM2D_WORK2(:, :) * PDZY(:,:,IKB+1)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)

!
CALL MXF2D_PHY(D, ZSHUGRADWK1_2D, ZMXF2D_WORK3)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
ZSHUGRADWK3_2D(:, :) = 0.5*(PDYY(:,:,IKB)+PDYY(:,:,IKB+1))
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)

!
CALL MXM2D_PHY(D, ZSHUGRADWK3_2D, ZMXM2D_WORK3)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
ZSHUGRADWK1_2D(:, :) = ZMXM2D_WORK3(:, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)

!
CALL MYF2D_PHY(D, ZSHUGRADWK1_2D, ZMYF2D_WORK3)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
ZSHUGRADWK3_2D(:, :) = PDZZ(:,:,IKB+1)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)

!
CALL MYM2D_PHY(D, ZSHUGRADWK3_2D, ZMYM2D_WORK2)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
ZSHUGRADWK1_2D(:, :) = (PVM(:,:,IKB+1)-PVM(:,:,IKB)) /                    ZMYM2D_WORK2(:, :) * PDZX(:,:,IKB+1)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)

!
CALL MYF2D_PHY(D, ZSHUGRADWK1_2D, ZMYF2D_WORK4)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
ZSHUGRADWK2_2D(:, :) = 0.5*(PDXX(:,:,IKB)+PDXX(:,:,IKB+1))
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)

!
CALL MYM2D_PHY(D, ZSHUGRADWK2_2D, ZMYM2D_WORK3)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
ZSHUGRADWK1_2D(:, :) = ZMYM2D_WORK3(:, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)

!
CALL MXF2D_PHY(D, ZSHUGRADWK1_2D, ZMXF2D_WORK4)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
ZWORK(:,:,IKB) =  -                                             &
     ZMXF2D_WORK1(:, :)   &
   *(ZMXF2D_WORK2(:, :)                                                       &  
    -ZMXF2D_WORK3(:, :) / ZMYF2D_WORK3(:, :)&
    -ZMYF2D_WORK4(:, :) / ZMXF2D_WORK4(:, :)&
    )   
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)

!
!
  ! dynamic production 
  
  DO JK=1, IKT
    DO JJ=1, IJT
      DO JI=1, IIT
        PDP(JI, JJ, JK) = PDP(JI, JJ, JK) + ZWORK(JI, JJ, JK)
      END DO
    END DO
  END DO
  
  ! 
END IF
!
! Storage in the LES configuration
!
IF (TLES%LLES_CALL .AND. KSPLT==1) THEN
  CALL SECOND_MNH(ZTIME1)
  CALL MYF_PHY(D, ZFLX(:, :, :), ZMYF3D_WORK1)
CALL MXF_PHY(D, ZMYF3D_WORK1, ZMXF3D_WORK1)
CALL LES_MEAN_SUBGRID( ZMXF3D_WORK1(:, :, :), TLES%X_LES_SUBGRID_UV )   
CALL GY_U_UV_DEVICE(PUM(:, :, :),PDYY(:, :, :),PDZZ(:, :, :),PDZY(:, :, :), ZGY_U_UV3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZSHUGRADWK2_3D(:, :, :) = ZGY_U_UV3D_WORK1(:, :, :)*ZFLX(:, :, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
CALL MYF_PHY(D, ZSHUGRADWK2_3D, ZMYF3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZSHUGRADWK1_3D(:, :, :) = ZMYF3D_WORK1(:, :, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
CALL MXF_PHY(D, ZSHUGRADWK1_3D, ZMXF3D_WORK1)
CALL LES_MEAN_SUBGRID( ZMXF3D_WORK1(:, :, :), TLES%X_LES_RES_ddxa_U_SBG_UaU , .TRUE.)  
CALL GX_V_UV_DEVICE(PVM(:, :, :),PDXX(:, :, :),PDZZ(:, :, :),PDZX(:, :, :), ZGX_V_UV3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZSHUGRADWK2_3D(:, :, :) = ZGX_V_UV3D_WORK1(:, :, :)*ZFLX(:, :, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
CALL MYF_PHY(D, ZSHUGRADWK2_3D, ZMYF3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZSHUGRADWK1_3D(:, :, :) = ZMYF3D_WORK1(:, :, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
CALL MXF_PHY(D, ZSHUGRADWK1_3D, ZMXF3D_WORK1)
CALL LES_MEAN_SUBGRID( ZMXF3D_WORK1(:, :, :), TLES%X_LES_RES_ddxa_V_SBG_UaV , .TRUE.)  
CALL SECOND_MNH(ZTIME2)
  TLES%XTIME_LES = TLES%XTIME_LES + ZTIME2 - ZTIME1
END IF
!
!
END SUBROUTINE TURB_HOR_UV
END MODULE MODE_TURB_HOR_UV
