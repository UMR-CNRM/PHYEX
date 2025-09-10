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

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      GX_V_UV_PVM(JI, JJ, JK) = ZGX_V_UV3D_WORK1(JI, JJ, JK)
    END DO
  END DO
END DO

!
IF (.NOT. O2D) THEN
  CALL GY_U_UV_DEVICE(PUM(:, :, :),PDYY(:, :, :),PDZZ(:, :, :),PDZY(:, :, :), ZGY_U_UV3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      GY_U_UV_PUM(JI, JJ, JK) = ZGY_U_UV3D_WORK1(JI, JJ, JK)
    END DO
  END DO
END DO

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

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZFLX(JI, JJ, JK)= - XCMFS * ZMYM3D_WORK1(JI, JJ, JK) *                           &
                (GY_U_UV_PUM(JI, JJ, JK) + GX_V_UV_PVM(JI, JJ, JK))
    END DO
  END DO
END DO

!
ELSE
  CALL MXM_PHY(D, PK(:, :, :), ZMXM3D_WORK1)
CALL MYM_PHY(D, ZMXM3D_WORK1, ZMYM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZFLX(JI, JJ, JK)= - XCMFS * ZMYM3D_WORK1(JI, JJ, JK) *                           &
                (GX_V_UV_PVM(JI, JJ, JK))
    END DO
  END DO
END DO

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

DO JJ=1, IJT
  DO JI=1, IIT
    ZSHUGRADWK2_2D(JI, JJ) = PDZZ(JI, JJ, IKB+1)
  END DO
END DO

!
CALL MXM2D_PHY(D, ZSHUGRADWK2_2D, ZMXM2D_WORK2)
CALL MYM2D_PHY(D, ZMXM2D_WORK1, ZMYM2D_WORK1)
CALL MXM2D_PHY(D, PDZZ(:,:,IKB), ZMXM2D_WORK4)

DO JJ=1, IJT
  DO JI=1, IIT
    ZSHUGRADWK1_2D(JI, JJ) = (PUM(JI, JJ, IKB+1)-PUM(JI, JJ, IKB))                   &
            *(1./ZMXM2D_WORK2(JI, JJ)+1./ZMXM2D_WORK4(JI, JJ))
  END DO
END DO

!
CALL MYM2D_PHY(D, ZSHUGRADWK1_2D, ZMYM2D_WORK2)
CALL DYM2D_PHY(D, PUM(:,:,IKB), ZDYM2D_WORK1)

DO JJ=1, IJT
  DO JI=1, IIT
    ZSHUGRADWK1_2D(JI, JJ) = (PDZY(JI, JJ, IKB+1)+PDZY(JI, JJ, IKB))
  END DO
END DO

!
CALL MXM2D_PHY(D, ZSHUGRADWK1_2D, ZMXM2D_WORK3)
CALL MXM2D_PHY(D, PDYY(:,:,IKB), ZMXM2D_WORK5)
CALL DXM2D_PHY(D, PVM(:,:,IKB), ZDXM2D_WORK1)

DO JJ=1, IJT
  DO JI=1, IIT
    ZSHUGRADWK2_2D(JI, JJ) = PDZZ(JI, JJ, IKB+1)
  END DO
END DO

!
CALL MYM2D_PHY(D, ZSHUGRADWK2_2D, ZMYM2D_WORK3)
CALL MYM2D_PHY(D, PDZZ(:,:,IKB), ZMYM2D_WORK5)

DO JJ=1, IJT
  DO JI=1, IIT
    ZSHUGRADWK1_2D(JI, JJ) = (PVM(JI, JJ, IKB+1)-PVM(JI, JJ, IKB))                   &
            *(1./ZMYM2D_WORK3(JI, JJ)+1./ZMYM2D_WORK5(JI, JJ))
  END DO
END DO

!
CALL MXM2D_PHY(D, ZSHUGRADWK1_2D, ZMXM2D_WORK6)

DO JJ=1, IJT
  DO JI=1, IIT
    ZSHUGRADWK1_2D(JI, JJ) = (PDZX(JI, JJ, IKB+1)+PDZX(JI, JJ, IKB))
  END DO
END DO

!
CALL MYM2D_PHY(D, ZSHUGRADWK1_2D, ZMYM2D_WORK4)
CALL MYM2D_PHY(D, PDXX(:,:,IKB), ZMYM2D_WORK6)

DO JJ=1, IJT
  DO JI=1, IIT
    ZFLX(JI, JJ, IKB)   = - XCMFS * ZMYM2D_WORK1(JI, JJ) *  (     &
      ( ZDYM2D_WORK1(JI, JJ)                                        &
       -ZMYM2D_WORK2(JI, JJ)&
        *0.5*ZMXM2D_WORK3(JI, JJ)            &
      ) / ZMXM2D_WORK5(JI, JJ)                                       &
     +( ZDXM2D_WORK1(JI, JJ)                                        &
       -ZMXM2D_WORK6(JI, JJ)&
        *0.5*ZMYM2D_WORK4(JI, JJ)            &
      ) / ZMYM2D_WORK6(JI, JJ)                                 ) 
  END DO
END DO

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

DO JJ=1, IJT
  DO JI=1, IIT
    ZSHUGRADWK2_2D(JI, JJ) = ZFLX(JI, JJ, IKB-1)
  END DO
END DO

!
CALL MYM2D_PHY(D, ZSHUGRADWK2_2D, ZMYM2D_WORK1)

DO JJ=1, IJT
  DO JI=1, IIT
    ZSHUGRADWK1_2D(JI, JJ) = ZMYM2D_WORK1(JI, JJ)
  END DO
END DO

!
CALL MXM2D_PHY(D, ZSHUGRADWK1_2D, ZMXM2D_WORK1)

DO JJ=1, IJT
  DO JI=1, IIT
    ZFLX(JI, JJ, IKB-1) = 2. * ZMXM2D_WORK1(JI, JJ)  &
                       - ZFLX(JI, JJ, IKB)
  END DO
END DO

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

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = ZMYM3D_WORK1(JI, JJ, JK) * PINV_PDYY(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MXM_PHY(D, ZSHUGRADWK2_3D, ZMXM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZFLX(JI, JJ, JK) * ZMXM3D_WORK1(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL DYF_PHY(D, ZSHUGRADWK1_3D, ZDYF3D_WORK1)
CALL MZM_PHY(D, ZFLX(:, :, :), ZMZM3D_WORK1)
CALL MZM_PHY(D, PDYY(:, :, :), ZMZM3D_WORK2)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = PDZY(JI, JJ, JK)/ZMZM3D_WORK2(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MXM_PHY(D, ZSHUGRADWK2_3D, ZMXM3D_WORK3)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = ZMZM3D_WORK1(JI, JJ, JK)*ZMXM3D_WORK3(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MYF_PHY(D, ZSHUGRADWK2_3D, ZMYF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = PMZM_PRHODJ(JI, JJ, JK) * PINV_PDZZ(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MXM_PHY(D, ZSHUGRADWK1_3D, ZMXM3D_WORK2)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMYF3D_WORK1(JI, JJ, JK)   &
                          * ZMXM3D_WORK2(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL DZF_PHY(D, ZSHUGRADWK1_3D, ZDZF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      PRUS(JI, JJ, JK) = PRUS(JI, JJ, JK)                                &
                    - ZDYF3D_WORK1(JI, JJ, JK)         &
                    + ZDZF3D_WORK1(JI, JJ, JK)
    END DO
  END DO
END DO

!
ELSE
  CALL MYM_PHY(D, PRHODJ(:, :, :), ZMYM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = ZMYM3D_WORK1(JI, JJ, JK) * PINV_PDYY(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MXM_PHY(D, ZSHUGRADWK2_3D, ZMXM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZFLX(JI, JJ, JK) * ZMXM3D_WORK1(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL DYF_PHY(D, ZSHUGRADWK1_3D, ZDYF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      PRUS(JI, JJ, JK) = PRUS(JI, JJ, JK) - ZDYF3D_WORK1(JI, JJ, JK)
    END DO
  END DO
END DO

!
END IF
!
!computation of the source for rho*V due to this flux
IF (.NOT. OFLAT) THEN
  CALL MXM_PHY(D, PRHODJ(:, :, :), ZMXM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = ZMXM3D_WORK1(JI, JJ, JK) * PINV_PDXX(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MYM_PHY(D, ZSHUGRADWK2_3D, ZMYM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZFLX(JI, JJ, JK) * ZMYM3D_WORK1(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL DXF_PHY(D, ZSHUGRADWK1_3D, ZDXF3D_WORK1)
CALL MZM_PHY(D, ZFLX(:, :, :), ZMZM3D_WORK1)
CALL MZM_PHY(D, PDXX(:, :, :), ZMZM3D_WORK2)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = PDZX(JI, JJ, JK)/ZMZM3D_WORK2(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MYM_PHY(D, ZSHUGRADWK2_3D, ZMYM3D_WORK3)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = ZMZM3D_WORK1(JI, JJ, JK)*ZMYM3D_WORK3(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MXF_PHY(D, ZSHUGRADWK2_3D, ZMXF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = PMZM_PRHODJ(JI, JJ, JK) * PINV_PDZZ(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MYM_PHY(D, ZSHUGRADWK1_3D, ZMYM3D_WORK2)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMXF3D_WORK1(JI, JJ, JK) &
                            * ZMYM3D_WORK2(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL DZF_PHY(D, ZSHUGRADWK1_3D, ZDZF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      PRVS(JI, JJ, JK) = PRVS(JI, JJ, JK)                             &
                      - ZDXF3D_WORK1(JI, JJ, JK)    &
                      + ZDZF3D_WORK1(JI, JJ, JK)
    END DO
  END DO
END DO

!
ELSE
  CALL MXM_PHY(D, PRHODJ(:, :, :), ZMXM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = ZMXM3D_WORK1(JI, JJ, JK) * PINV_PDXX(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MYM_PHY(D, ZSHUGRADWK2_3D, ZMYM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZFLX(JI, JJ, JK) * ZMYM3D_WORK1(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL DXF_PHY(D, ZSHUGRADWK1_3D, ZDXF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      PRVS(JI, JJ, JK) = PRVS(JI, JJ, JK) - ZDXF3D_WORK1(JI, JJ, JK)
    END DO
  END DO
END DO

!
END IF
!
IF (KSPLT==1) THEN
  !
  !Contribution to the dynamic production of TKE:
  !
  IF (.NOT. O2D) THEN
    
DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = ZFLX(JI, JJ, JK) *                                &
             (GY_U_UV_PUM(JI, JJ, JK) + GX_V_UV_PVM(JI, JJ, JK))
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
CALL MXF_PHY(D, ZSHUGRADWK1_3D, ZMXF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZWORK(JI, JJ, JK) = - ZMXF3D_WORK1(JI, JJ, JK)   
    END DO
  END DO
END DO

!
ELSE
    
DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = ZFLX(JI, JJ, JK) *                                &
             (GX_V_UV_PVM(JI, JJ, JK))
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
CALL MXF_PHY(D, ZSHUGRADWK1_3D, ZMXF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZWORK(JI, JJ, JK) = - ZMXF3D_WORK1(JI, JJ, JK)    
    END DO
  END DO
END DO

!
ENDIF
  !
  ! evaluate the dynamic production at w(IKB+1) in PDP(IKB)
  !
  
DO JJ=1, IJT
  DO JI=1, IIT
    ZSHUGRADWK2_2D(JI, JJ) = 0.5 * (ZFLX(JI, JJ, IKB+1)+ZFLX(JI, JJ, IKB))
  END DO
END DO

!
CALL MYF2D_PHY(D, ZSHUGRADWK2_2D, ZMYF2D_WORK1)

DO JJ=1, IJT
  DO JI=1, IIT
    ZSHUGRADWK1_2D(JI, JJ) = ZMYF2D_WORK1(JI, JJ)
  END DO
END DO

!
CALL MXF2D_PHY(D, ZSHUGRADWK1_2D, ZMXF2D_WORK1)

DO JJ=1, IJT
  DO JI=1, IIT
    ZSHUGRADWK3_2D(JI, JJ) = 0.5 * (PUM(JI, JJ, IKB+1)+PUM(JI, JJ, IKB))
  END DO
END DO

!
CALL DYM2D_PHY(D, ZSHUGRADWK3_2D, ZDYM2D_WORK1)

DO JJ=1, IJT
  DO JI=1, IIT
    ZSHUGRADWK2_2D(JI, JJ) = 0.5*(PDYY(JI, JJ, IKB)+PDYY(JI, JJ, IKB+1))
  END DO
END DO

!
CALL MXM2D_PHY(D, ZSHUGRADWK2_2D, ZMXM2D_WORK1)

DO JJ=1, IJT
  DO JI=1, IIT
    ZSHUGRADWK2_2D(JI, JJ) = 0.5 * (PVM(JI, JJ, IKB+1)+PVM(JI, JJ, IKB))
  END DO
END DO

!
CALL DXM2D_PHY(D, ZSHUGRADWK2_2D, ZDXM2D_WORK1)

DO JJ=1, IJT
  DO JI=1, IIT
    ZSHUGRADWK2_2D(JI, JJ) = 0.5*(PDXX(JI, JJ, IKB)+PDXX(JI, JJ, IKB+1))
  END DO
END DO

!
CALL MYM2D_PHY(D, ZSHUGRADWK2_2D, ZMYM2D_WORK1)

DO JJ=1, IJT
  DO JI=1, IIT
    ZSHUGRADWK2_2D(JI, JJ) = ZDYM2D_WORK1(JI, JJ)           / ZMXM2D_WORK1(JI, JJ)           +ZDXM2D_WORK1(JI, JJ)           / ZMYM2D_WORK1(JI, JJ)
  END DO
END DO

!
CALL MYF2D_PHY(D, ZSHUGRADWK2_2D, ZMYF2D_WORK2)

DO JJ=1, IJT
  DO JI=1, IIT
    ZSHUGRADWK1_2D(JI, JJ) = ZMYF2D_WORK2(JI, JJ)
  END DO
END DO

!
CALL MXF2D_PHY(D, ZSHUGRADWK1_2D, ZMXF2D_WORK2)

DO JJ=1, IJT
  DO JI=1, IIT
    ZSHUGRADWK3_2D(JI, JJ) = PDZZ(JI, JJ, IKB+1)
  END DO
END DO

!
CALL MXM2D_PHY(D, ZSHUGRADWK3_2D, ZMXM2D_WORK2)

DO JJ=1, IJT
  DO JI=1, IIT
    ZSHUGRADWK1_2D(JI, JJ) = (PUM(JI, JJ, IKB+1)-PUM(JI, JJ, IKB)) /                    ZMXM2D_WORK2(JI, JJ) * PDZY(JI, JJ, IKB+1)
  END DO
END DO

!
CALL MXF2D_PHY(D, ZSHUGRADWK1_2D, ZMXF2D_WORK3)

DO JJ=1, IJT
  DO JI=1, IIT
    ZSHUGRADWK3_2D(JI, JJ) = 0.5*(PDYY(JI, JJ, IKB)+PDYY(JI, JJ, IKB+1))
  END DO
END DO

!
CALL MXM2D_PHY(D, ZSHUGRADWK3_2D, ZMXM2D_WORK3)

DO JJ=1, IJT
  DO JI=1, IIT
    ZSHUGRADWK1_2D(JI, JJ) = ZMXM2D_WORK3(JI, JJ)
  END DO
END DO

!
CALL MYF2D_PHY(D, ZSHUGRADWK1_2D, ZMYF2D_WORK3)

DO JJ=1, IJT
  DO JI=1, IIT
    ZSHUGRADWK3_2D(JI, JJ) = PDZZ(JI, JJ, IKB+1)
  END DO
END DO

!
CALL MYM2D_PHY(D, ZSHUGRADWK3_2D, ZMYM2D_WORK2)

DO JJ=1, IJT
  DO JI=1, IIT
    ZSHUGRADWK1_2D(JI, JJ) = (PVM(JI, JJ, IKB+1)-PVM(JI, JJ, IKB)) /                    ZMYM2D_WORK2(JI, JJ) * PDZX(JI, JJ, IKB+1)
  END DO
END DO

!
CALL MYF2D_PHY(D, ZSHUGRADWK1_2D, ZMYF2D_WORK4)

DO JJ=1, IJT
  DO JI=1, IIT
    ZSHUGRADWK2_2D(JI, JJ) = 0.5*(PDXX(JI, JJ, IKB)+PDXX(JI, JJ, IKB+1))
  END DO
END DO

!
CALL MYM2D_PHY(D, ZSHUGRADWK2_2D, ZMYM2D_WORK3)

DO JJ=1, IJT
  DO JI=1, IIT
    ZSHUGRADWK1_2D(JI, JJ) = ZMYM2D_WORK3(JI, JJ)
  END DO
END DO

!
CALL MXF2D_PHY(D, ZSHUGRADWK1_2D, ZMXF2D_WORK4)

DO JJ=1, IJT
  DO JI=1, IIT
    ZWORK(JI, JJ, IKB) =  -                                             &
         ZMXF2D_WORK1(JI, JJ)   &
       *(ZMXF2D_WORK2(JI, JJ)                                                       &  
        -ZMXF2D_WORK3(JI, JJ) / ZMYF2D_WORK3(JI, JJ)&
        -ZMYF2D_WORK4(JI, JJ) / ZMXF2D_WORK4(JI, JJ)&
        )   
  END DO
END DO

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

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = ZGY_U_UV3D_WORK1(JI, JJ, JK)*ZFLX(JI, JJ, JK)
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
CALL MXF_PHY(D, ZSHUGRADWK1_3D, ZMXF3D_WORK1)
CALL LES_MEAN_SUBGRID( ZMXF3D_WORK1(:, :, :), TLES%X_LES_RES_ddxa_U_SBG_UaU , .TRUE.)  
CALL GX_V_UV_DEVICE(PVM(:, :, :),PDXX(:, :, :),PDZZ(:, :, :),PDZX(:, :, :), ZGX_V_UV3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = ZGX_V_UV3D_WORK1(JI, JJ, JK)*ZFLX(JI, JJ, JK)
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
CALL MXF_PHY(D, ZSHUGRADWK1_3D, ZMXF3D_WORK1)
CALL LES_MEAN_SUBGRID( ZMXF3D_WORK1(:, :, :), TLES%X_LES_RES_ddxa_V_SBG_UaV , .TRUE.)  
CALL SECOND_MNH(ZTIME2)
  TLES%XTIME_LES = TLES%XTIME_LES + ZTIME2 - ZTIME1
END IF
!
!
END SUBROUTINE TURB_HOR_UV
END MODULE MODE_TURB_HOR_UV
