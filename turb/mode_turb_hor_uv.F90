!MNH_LIC Copyright 1994-2025 CNRS, Meteo-France and Universite Paul Sabatier
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
USE MODD_CTURB,          ONLY: CTURB_XCMFS => XCMFS
USE MODD_FIELD,          ONLY: TFIELDMETADATA, TYPEREAL
USE MODD_IO,             ONLY: TFILEDATA
USE MODD_PARAMETERS, ONLY: JPVEXT
USE MODD_LES, ONLY: TLES_t
!
USE MODE_IO_FIELD_WRITE_PHY, ONLY: IO_FIELD_WRITE_PHY
!
USE MODI_LES_MEAN_SUBGRID_PHY, ONLY: LES_MEAN_SUBGRID_PHY
!
USE MODI_SECOND_MNH, ONLY: SECOND_MNH
USE MODE_SHUMAN_PHY, ONLY: DXF_PHY, DXM2D_PHY, DYF_PHY, DYM2D_PHY, DZF_PHY, MXF2D_PHY, MXF_PHY, MXM2D_PHY, MXM_PHY, MYF2D_PHY, MYF_PHY, MYM2D_PHY, MYM_PHY, MZM_PHY
USE MODE_GRADIENT_V_PHY, ONLY: GX_V_UV_PHY
USE MODE_GRADIENT_U_PHY, ONLY: GY_U_UV_PHY
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
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    ::  PK          ! Turbulent diffusion doef.
                                                        ! PK = PLM * SQRT(PTKEM)
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    ::  PINV_PDXX   ! 1./PDXX
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    ::  PINV_PDYY   ! 1./PDYY
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    ::  PINV_PDZZ   ! 1./PDZZ
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    ::  PMZM_PRHODJ ! MZM(PRHODJ)
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    ::  PDXX, PDYY, PDZZ, PDZX, PDZY 
                                                         ! Metric coefficients
REAL, DIMENSION(D%NIJT),     INTENT(IN)    ::  PDIRCOSZW
! Director Cosinus along z directions at surface w-point
REAL, DIMENSION(D%NIJT),   INTENT(IN)   ::  PCOSSLOPE       ! cosinus of the angle 
                                      ! between i and the slope vector
REAL, DIMENSION(D%NIJT),   INTENT(IN)   ::  PSINSLOPE       ! sinus of the angle 
                                      ! between i and the slope vector
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    ::  PRHODJ       ! density * grid volume
!
REAL, DIMENSION(D%NIJT),   INTENT(IN)   ::  PCDUEFF      ! Cd * || u || at time t
REAL, DIMENSION(D%NIJT),   INTENT(IN)   ::  PTAU11M      ! <uu> in the axes linked 
       ! to the maximum slope direction and the surface normal and the binormal 
       ! at time t - dt
REAL, DIMENSION(D%NIJT),   INTENT(IN)   ::  PTAU12M      ! <uv> in the same axes
REAL, DIMENSION(D%NIJT),   INTENT(IN)   ::  PTAU22M      ! <vv> in the same axes
REAL, DIMENSION(D%NIJT),   INTENT(IN)   ::  PTAU33M      ! <ww> in the same axes
!
! Variables at t-1
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    ::  PUM,PVM
REAL, DIMENSION(D%NIJT),      INTENT(IN)   ::  PUSLOPEM     ! wind component along the 
                                     ! maximum slope direction
REAL, DIMENSION(D%NIJT),      INTENT(IN)   ::  PVSLOPEM     ! wind component along the 
                                     ! direction normal to the maximum slope one
!
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) ::  PRUS, PRVS   ! var. at t+1 -split-
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) ::  PDP          ! TKE production terms
!
!
!
!*       0.2  declaration of local variables
!
REAL, DIMENSION(D%NIJT,D%NKT) :: ZFLX,ZWORK ! work arrays
!   
REAL, DIMENSION(D%NIJT) :: ZDIRSINZW 
      ! sinus of the angle between the vertical and the normal to the orography
INTEGER             :: IKB,IKE,IIJT,IKT, JIJ,JK
                                    ! Index values for the Beginning and End
                                    ! mass points of the domain  
!
REAL, DIMENSION(D%NIJT,D%NKT)  :: GY_U_UV_PUM
REAL, DIMENSION(D%NIJT,D%NKT)  :: GX_V_UV_PVM
!
REAL :: ZTIME1, ZTIME2
TYPE(TFIELDMETADATA) :: TZFIELD
!
REAL, DIMENSION(D%NIJT) :: ZWORKSHUGRAD_2D_1
REAL, DIMENSION(D%NIJT) :: ZWORKSHUGRAD_2D_2
REAL, DIMENSION(D%NIJT) :: ZWORKSHUGRAD_2D_3
REAL, DIMENSION(D%NIJT) :: ZWORKSHUGRAD_2D_4
REAL, DIMENSION(D%NIJT) :: ZWORKSHUGRAD_2D_5
REAL, DIMENSION(D%NIJT) :: ZWORKSHUGRAD_2D_6
REAL, DIMENSION(D%NIJT) :: ZWORKSHUGRAD_2D_7
REAL, DIMENSION(D%NIJT) :: ZWORKSHUGRAD_2D_8
REAL, DIMENSION(D%NIJT) :: ZWORKSHUGRAD_2D_9
REAL, DIMENSION(D%NIJT) :: ZWORKSHUGRAD_2D_10
REAL, DIMENSION(D%NIJT) :: ZWORKSHUGRAD_2D_11
REAL, DIMENSION(D%NIJT) :: ZWORKSHUGRAD_2D_12
REAL, DIMENSION(D%NIJT) :: ZWORKSHUGRAD_2D_13
REAL, DIMENSION(D%NIJT) :: ZWORKSHUGRAD_2D_14
REAL, DIMENSION(D%NIJT) :: ZWORKSHUGRAD_2D_15
REAL, DIMENSION(D%NIJT) :: ZWORKSHUGRAD_2D_16
REAL, DIMENSION(D%NIJT) :: ZWORKSHUGRAD_2D_17
REAL, DIMENSION(D%NIJT) :: ZWORKSHUGRAD_2D_18
REAL, DIMENSION(D%NIJT) :: ZWORKSHUGRAD_2D_19
REAL, DIMENSION(D%NIJT) :: ZWORKSHUGRAD_2D_20
REAL, DIMENSION(D%NIJT) :: ZWORKSHUGRAD_2D_21
REAL, DIMENSION(D%NIJT) :: ZWORKSHUGRAD_2D_22
REAL, DIMENSION(D%NIJT) :: ZWORKSHUGRAD_2D_23
REAL, DIMENSION(D%NIJT) :: ZWORKSHUGRAD_2D_24
REAL, DIMENSION(D%NIJT) :: ZWORKSHUGRAD_2D_25
REAL, DIMENSION(D%NIJT) :: ZWORKSHUGRAD_2D_26
REAL, DIMENSION(D%NIJT, D%NKT) :: ZWORKSHUGRAD_3D_1
REAL, DIMENSION(D%NIJT, D%NKT) :: ZWORKSHUGRAD_3D_2
REAL, DIMENSION(D%NIJT, D%NKT) :: ZWORKSHUGRAD_3D_3
REAL, DIMENSION(D%NIJT, D%NKT) :: ZWORKSHUGRAD_3D_4
REAL, DIMENSION(D%NIJT, D%NKT) :: ZWORKSHUGRAD_3D_5
REAL, DIMENSION(D%NIJT, D%NKT) :: ZWORKSHUGRAD_3D_6
REAL, DIMENSION(D%NIJT, D%NKT) :: ZWORKSHUGRAD_3D_7
REAL, DIMENSION(D%NIJT, D%NKT) :: ZWORKSHUGRAD_3D_8
REAL, DIMENSION(D%NIJT, D%NKT) :: ZWORKSHUGRAD_3D_9
REAL, DIMENSION(D%NIJT, D%NKT) :: ZWORKSHUGRAD_3D_10
REAL, DIMENSION(D%NIJT, D%NKT) :: ZWORKSHUGRAD_3D_11
REAL, DIMENSION(D%NIJT, D%NKT) :: ZWORKSHUGRAD_3D_12
REAL, DIMENSION(D%NIJT, D%NKT) :: ZWORKSHUGRAD_3D_13
REAL, DIMENSION(D%NIJT, D%NKT) :: ZWORKSHUGRAD_3D_14
REAL, DIMENSION(D%NIJT, D%NKT) :: ZWORKSHUGRAD_3D_15
REAL :: XCMFS
!
! ---------------------------------------------------------------------------
!
!*       1.   PRELIMINARY COMPUTATIONS
!             ------------------------
!
XCMFS = CTURB_XCMFS
!
IKB = 1+JPVEXT               
IKE = SIZE(PUM,2)-JPVEXT
IIJT=D%NIJT
IKT=D%NKT  
!

DO JIJ=1, IIJT
  ZDIRSINZW(JIJ) = SQRT( 1. - PDIRCOSZW(JIJ)**2 )
END DO

!
CALL GX_V_UV_PHY(D, OFLAT, PVM, PDXX, PDZZ, PDZX, GX_V_UV_PVM)
IF (.NOT. O2D) CALL GY_U_UV_PHY(D, OFLAT, PUM, PDYY, PDZZ, PDZY, GY_U_UV_PUM)
!
!
!*      12.   < U'V'>
!             -------
!
!
IF (.NOT. O2D) THEN
  CALL MXM_PHY(D, PK, ZWORKSHUGRAD_3D_1)
  CALL MYM_PHY(D, ZWORKSHUGRAD_3D_1, ZWORKSHUGRAD_3D_2)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      ZFLX(JIJ, JK)= - XCMFS * ZWORKSHUGRAD_3D_2(JIJ, JK) *                           &
              (GY_U_UV_PUM(JIJ, JK) + GX_V_UV_PVM(JIJ, JK))
    END DO
  END DO
  
ELSE
  CALL MXM_PHY(D, PK, ZWORKSHUGRAD_3D_1)
  CALL MYM_PHY(D, ZWORKSHUGRAD_3D_1, ZWORKSHUGRAD_3D_2)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      ZFLX(JIJ, JK)= - XCMFS * ZWORKSHUGRAD_3D_2(JIJ, JK) *                           &
              (GX_V_UV_PVM(JIJ, JK))
    END DO
  END DO
  
END IF
!

DO JIJ=1, IIJT
  ZFLX(JIJ, IKE+1)= ZFLX(JIJ, IKE)
END DO

!
!
! Compute the correlation at the first physical level with the following 
! hypothesis du/dz vary in 1/z and w=0 at the ground
CALL MXM2D_PHY(D, PDZZ(:,IKB+1), ZWORKSHUGRAD_2D_1)
CALL MXM2D_PHY(D, PDZZ(:,IKB), ZWORKSHUGRAD_2D_2)
CALL MYM2D_PHY(D, PDZZ(:,IKB+1), ZWORKSHUGRAD_2D_3)
CALL MYM2D_PHY(D, PDZZ(:,IKB), ZWORKSHUGRAD_2D_4)

DO JIJ=1, D%NIJT
  ZWORKSHUGRAD_2D_5(JIJ)=PDZY(JIJ, IKB+1)+PDZY(JIJ, IKB)
END DO

CALL MXM2D_PHY(D, ZWORKSHUGRAD_2D_5, ZWORKSHUGRAD_2D_6)

DO JIJ=1, D%NIJT
  ZWORKSHUGRAD_2D_7(JIJ)=PDZX(JIJ, IKB+1)+PDZX(JIJ, IKB)
END DO

CALL MYM2D_PHY(D, ZWORKSHUGRAD_2D_7, ZWORKSHUGRAD_2D_8)

DO JIJ=1, D%NIJT
  ZWORKSHUGRAD_2D_9(JIJ)=(PUM(JIJ, IKB+1)-PUM(JIJ, IKB))                   &
          *(1./ZWORKSHUGRAD_2D_1(JIJ)+1./ZWORKSHUGRAD_2D_2(JIJ))
END DO

CALL MYM2D_PHY(D, ZWORKSHUGRAD_2D_9, ZWORKSHUGRAD_2D_10)

DO JIJ=1, D%NIJT
  ZWORKSHUGRAD_2D_11(JIJ)=(PVM(JIJ, IKB+1)-PVM(JIJ, IKB))                   &
          *(1./ZWORKSHUGRAD_2D_3(JIJ)+1./ZWORKSHUGRAD_2D_4(JIJ))
END DO

CALL MXM2D_PHY(D, ZWORKSHUGRAD_2D_11, ZWORKSHUGRAD_2D_12)
CALL MXM2D_PHY(D, PK(:,IKB), ZWORKSHUGRAD_2D_13)
CALL DYM2D_PHY(D, PUM(:,IKB), ZWORKSHUGRAD_2D_14)
CALL DXM2D_PHY(D, PVM(:,IKB), ZWORKSHUGRAD_2D_15)
CALL MXM2D_PHY(D, PDYY(:,IKB), ZWORKSHUGRAD_2D_16)
CALL MYM2D_PHY(D, PDXX(:,IKB), ZWORKSHUGRAD_2D_17)
CALL MYM2D_PHY(D, ZWORKSHUGRAD_2D_13, ZWORKSHUGRAD_2D_18)

DO JIJ=1, D%NIJT
  ZFLX(JIJ, IKB)   = - XCMFS * ZWORKSHUGRAD_2D_18(JIJ) *  (     &
    ( ZWORKSHUGRAD_2D_14(JIJ)                                        &
     -ZWORKSHUGRAD_2D_10(JIJ)&
      *0.5*ZWORKSHUGRAD_2D_6(JIJ)            &
    ) / ZWORKSHUGRAD_2D_16(JIJ)                                       &
   +( ZWORKSHUGRAD_2D_15(JIJ)                                        &
     -ZWORKSHUGRAD_2D_12(JIJ)&
      *0.5*ZWORKSHUGRAD_2D_8(JIJ)            &
    ) / ZWORKSHUGRAD_2D_17(JIJ)                                 )
END DO

! 
! extrapolates this flux under the ground with the surface flux

DO JIJ=1, IIJT
  ZFLX(JIJ, IKB-1) =                                                           &
     PTAU11M(JIJ) * PCOSSLOPE(JIJ) * PSINSLOPE(JIJ) * PDIRCOSZW(JIJ)**2         &
    +PTAU12M(JIJ) * (PCOSSLOPE(JIJ)**2 - PSINSLOPE(JIJ)**2) *                   &
                    PDIRCOSZW(JIJ)**2                                           &
    -PTAU22M(JIJ) * PCOSSLOPE(JIJ) * PSINSLOPE(JIJ)                             &
    +PTAU33M(JIJ) * PCOSSLOPE(JIJ) * PSINSLOPE(JIJ) * ZDIRSINZW(JIJ)**2         &
    -PCDUEFF(JIJ) * (                                                           &
      2. * PUSLOPEM(JIJ) * PCOSSLOPE(JIJ) * PSINSLOPE(JIJ) *                    &
            PDIRCOSZW(JIJ) * ZDIRSINZW(JIJ)                                     &
      +PVSLOPEM(JIJ) * (PCOSSLOPE(JIJ)**2 - PSINSLOPE(JIJ)**2) * ZDIRSINZW(JIJ) &
                     )
  !
END DO

!
CALL MYM2D_PHY(D, ZFLX(:,IKB-1), ZWORKSHUGRAD_2D_1)
CALL MXM2D_PHY(D, ZWORKSHUGRAD_2D_1, ZWORKSHUGRAD_2D_2)

DO JIJ=1, D%NIJT
  ZFLX(JIJ, IKB-1) = 2. * ZWORKSHUGRAD_2D_2(JIJ)  &
                     - ZFLX(JIJ, IKB)
END DO

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

  CALL IO_FIELD_WRITE_PHY(D,TPFILE,TZFIELD,ZFLX)
END IF
!
!
!
!computation of the source for rho*V due to this flux
IF (.NOT. OFLAT) THEN
  CALL MZM_PHY(D, PDYY, ZWORKSHUGRAD_3D_1)
  CALL MYM_PHY(D, PRHODJ, ZWORKSHUGRAD_3D_2)
  CALL MZM_PHY(D, ZFLX, ZWORKSHUGRAD_3D_3)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      ZWORKSHUGRAD_3D_4(JIJ, JK)=PDZY(JIJ, JK)/ZWORKSHUGRAD_3D_1(JIJ, JK)
    END DO
  END DO
  
  CALL MXM_PHY(D, ZWORKSHUGRAD_3D_4, ZWORKSHUGRAD_3D_5)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      ZWORKSHUGRAD_3D_6(JIJ, JK)=ZWORKSHUGRAD_3D_2(JIJ, JK) * PINV_PDYY(JIJ, JK)
    END DO
  END DO
  
  CALL MXM_PHY(D, ZWORKSHUGRAD_3D_6, ZWORKSHUGRAD_3D_7)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      ZWORKSHUGRAD_3D_8(JIJ, JK)=ZWORKSHUGRAD_3D_3(JIJ, JK)*ZWORKSHUGRAD_3D_5(JIJ, JK)
    END DO
  END DO
  
  CALL MYF_PHY(D, ZWORKSHUGRAD_3D_8, ZWORKSHUGRAD_3D_9)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      ZWORKSHUGRAD_3D_10(JIJ, JK)=PMZM_PRHODJ(JIJ, JK) * PINV_PDZZ(JIJ, JK)
    END DO
  END DO
  
  CALL MXM_PHY(D, ZWORKSHUGRAD_3D_10, ZWORKSHUGRAD_3D_11)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      ZWORKSHUGRAD_3D_12(JIJ, JK)=ZFLX(JIJ, JK) * ZWORKSHUGRAD_3D_7(JIJ, JK)
    END DO
  END DO
  
  CALL DYF_PHY(D, ZWORKSHUGRAD_3D_12, ZWORKSHUGRAD_3D_13)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      ZWORKSHUGRAD_3D_14(JIJ, JK)=ZWORKSHUGRAD_3D_9(JIJ, JK)   &
                        * ZWORKSHUGRAD_3D_11(JIJ, JK)
    END DO
  END DO
  
  CALL DZF_PHY(D, ZWORKSHUGRAD_3D_14, ZWORKSHUGRAD_3D_15)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      PRUS(JIJ, JK) = PRUS(JIJ, JK)                                &
                  - ZWORKSHUGRAD_3D_13(JIJ, JK)         &
                  + ZWORKSHUGRAD_3D_15(JIJ, JK)
    END DO
  END DO
  
ELSE
  CALL MYM_PHY(D, PRHODJ, ZWORKSHUGRAD_3D_1)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      ZWORKSHUGRAD_3D_2(JIJ, JK)=ZWORKSHUGRAD_3D_1(JIJ, JK) * PINV_PDYY(JIJ, JK)
    END DO
  END DO
  
  CALL MXM_PHY(D, ZWORKSHUGRAD_3D_2, ZWORKSHUGRAD_3D_3)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      ZWORKSHUGRAD_3D_4(JIJ, JK)=ZFLX(JIJ, JK) * ZWORKSHUGRAD_3D_3(JIJ, JK)
    END DO
  END DO
  
  CALL DYF_PHY(D, ZWORKSHUGRAD_3D_4, ZWORKSHUGRAD_3D_5)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      PRUS(JIJ, JK) = PRUS(JIJ, JK) - ZWORKSHUGRAD_3D_5(JIJ, JK)
    END DO
  END DO
  
END IF
!
!computation of the source for rho*V due to this flux
IF (.NOT. OFLAT) THEN
  CALL MZM_PHY(D, PDXX, ZWORKSHUGRAD_3D_1)
  CALL MXM_PHY(D, PRHODJ, ZWORKSHUGRAD_3D_2)
  CALL MZM_PHY(D, ZFLX, ZWORKSHUGRAD_3D_3)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      ZWORKSHUGRAD_3D_4(JIJ, JK)=PDZX(JIJ, JK)/ZWORKSHUGRAD_3D_1(JIJ, JK)
    END DO
  END DO
  
  CALL MYM_PHY(D, ZWORKSHUGRAD_3D_4, ZWORKSHUGRAD_3D_5)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      ZWORKSHUGRAD_3D_6(JIJ, JK)=ZWORKSHUGRAD_3D_2(JIJ, JK) * PINV_PDXX(JIJ, JK)
    END DO
  END DO
  
  CALL MYM_PHY(D, ZWORKSHUGRAD_3D_6, ZWORKSHUGRAD_3D_7)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      ZWORKSHUGRAD_3D_8(JIJ, JK)=ZWORKSHUGRAD_3D_3(JIJ, JK)*ZWORKSHUGRAD_3D_5(JIJ, JK)
    END DO
  END DO
  
  CALL MXF_PHY(D, ZWORKSHUGRAD_3D_8, ZWORKSHUGRAD_3D_9)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      ZWORKSHUGRAD_3D_10(JIJ, JK)=PMZM_PRHODJ(JIJ, JK) * PINV_PDZZ(JIJ, JK)
    END DO
  END DO
  
  CALL MYM_PHY(D, ZWORKSHUGRAD_3D_10, ZWORKSHUGRAD_3D_11)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      ZWORKSHUGRAD_3D_12(JIJ, JK)=ZFLX(JIJ, JK) * ZWORKSHUGRAD_3D_7(JIJ, JK)
    END DO
  END DO
  
  CALL DXF_PHY(D, ZWORKSHUGRAD_3D_12, ZWORKSHUGRAD_3D_13)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      ZWORKSHUGRAD_3D_14(JIJ, JK)=ZWORKSHUGRAD_3D_9(JIJ, JK) &
                          * ZWORKSHUGRAD_3D_11(JIJ, JK)
    END DO
  END DO
  
  CALL DZF_PHY(D, ZWORKSHUGRAD_3D_14, ZWORKSHUGRAD_3D_15)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      PRVS(JIJ, JK) = PRVS(JIJ, JK)                             &
                    - ZWORKSHUGRAD_3D_13(JIJ, JK)    &
                    + ZWORKSHUGRAD_3D_15(JIJ, JK)
    END DO
  END DO
  
ELSE
  CALL MXM_PHY(D, PRHODJ, ZWORKSHUGRAD_3D_1)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      ZWORKSHUGRAD_3D_2(JIJ, JK)=ZWORKSHUGRAD_3D_1(JIJ, JK) * PINV_PDXX(JIJ, JK)
    END DO
  END DO
  
  CALL MYM_PHY(D, ZWORKSHUGRAD_3D_2, ZWORKSHUGRAD_3D_3)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      ZWORKSHUGRAD_3D_4(JIJ, JK)=ZFLX(JIJ, JK) * ZWORKSHUGRAD_3D_3(JIJ, JK)
    END DO
  END DO
  
  CALL DXF_PHY(D, ZWORKSHUGRAD_3D_4, ZWORKSHUGRAD_3D_5)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      PRVS(JIJ, JK) = PRVS(JIJ, JK) - ZWORKSHUGRAD_3D_5(JIJ, JK)
    END DO
  END DO
  
END IF
!
IF (KSPLT==1) THEN
  !
  !Contribution to the dynamic production of TKE:
  !
  IF (.NOT. O2D) THEN
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWORKSHUGRAD_3D_1(JIJ, JK)=ZFLX(JIJ, JK) *                                &
           (GY_U_UV_PUM(JIJ, JK) + GX_V_UV_PVM(JIJ, JK))
      END DO
    END DO
    
    CALL MYF_PHY(D, ZWORKSHUGRAD_3D_1, ZWORKSHUGRAD_3D_2)
    CALL MXF_PHY(D, ZWORKSHUGRAD_3D_2, ZWORKSHUGRAD_3D_3)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWORK(JIJ, JK) = - ZWORKSHUGRAD_3D_3(JIJ, JK)
      END DO
    END DO
    
  ELSE
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWORKSHUGRAD_3D_1(JIJ, JK)=ZFLX(JIJ, JK) *                                &
           (GX_V_UV_PVM(JIJ, JK))
      END DO
    END DO
    
    CALL MYF_PHY(D, ZWORKSHUGRAD_3D_1, ZWORKSHUGRAD_3D_2)
    CALL MXF_PHY(D, ZWORKSHUGRAD_3D_2, ZWORKSHUGRAD_3D_3)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWORK(JIJ, JK) = - ZWORKSHUGRAD_3D_3(JIJ, JK)
      END DO
    END DO
    
  ENDIF
  !
  ! evaluate the dynamic production at w(IKB+1) in PDP(IKB)
  !
  
  DO JIJ=1, D%NIJT
    ZWORKSHUGRAD_2D_1(JIJ)=0.5 * (PUM(JIJ, IKB+1)+PUM(JIJ, IKB))
  END DO
  
  CALL DYM2D_PHY(D, ZWORKSHUGRAD_2D_1, ZWORKSHUGRAD_2D_2)
  
  DO JIJ=1, D%NIJT
    ZWORKSHUGRAD_2D_3(JIJ)=0.5*(PDYY(JIJ, IKB)+PDYY(JIJ, IKB+1))
  END DO
  
  CALL MXM2D_PHY(D, ZWORKSHUGRAD_2D_3, ZWORKSHUGRAD_2D_4)
  
  DO JIJ=1, D%NIJT
    ZWORKSHUGRAD_2D_5(JIJ)=0.5 * (PVM(JIJ, IKB+1)+PVM(JIJ, IKB))
  END DO
  
  CALL DXM2D_PHY(D, ZWORKSHUGRAD_2D_5, ZWORKSHUGRAD_2D_6)
  
  DO JIJ=1, D%NIJT
    ZWORKSHUGRAD_2D_7(JIJ)=0.5*(PDXX(JIJ, IKB)+PDXX(JIJ, IKB+1))
  END DO
  
  CALL MYM2D_PHY(D, ZWORKSHUGRAD_2D_7, ZWORKSHUGRAD_2D_8)
  CALL MXM2D_PHY(D, PDZZ(:,IKB+1), ZWORKSHUGRAD_2D_9)
  CALL MYM2D_PHY(D, PDZZ(:,IKB+1), ZWORKSHUGRAD_2D_10)
  
  DO JIJ=1, D%NIJT
    ZWORKSHUGRAD_2D_11(JIJ)=0.5*(PDYY(JIJ, IKB)+PDYY(JIJ, IKB+1))
  END DO
  
  CALL MXM2D_PHY(D, ZWORKSHUGRAD_2D_11, ZWORKSHUGRAD_2D_12)
  
  DO JIJ=1, D%NIJT
    ZWORKSHUGRAD_2D_13(JIJ)=ZWORKSHUGRAD_2D_2(JIJ)          &
        / ZWORKSHUGRAD_2D_4(JIJ)          &
        +ZWORKSHUGRAD_2D_6(JIJ)          &
        / ZWORKSHUGRAD_2D_8(JIJ)
  END DO
  
  CALL MYF2D_PHY(D, ZWORKSHUGRAD_2D_13, ZWORKSHUGRAD_2D_14)
  
  DO JIJ=1, D%NIJT
    ZWORKSHUGRAD_2D_15(JIJ)=0.5*(PDXX(JIJ, IKB)+PDXX(JIJ, IKB+1))
  END DO
  
  CALL MYM2D_PHY(D, ZWORKSHUGRAD_2D_15, ZWORKSHUGRAD_2D_16)
  
  DO JIJ=1, D%NIJT
    ZWORKSHUGRAD_2D_17(JIJ)=0.5 * (ZFLX(JIJ, IKB+1)+ZFLX(JIJ, IKB))
  END DO
  
  CALL MYF2D_PHY(D, ZWORKSHUGRAD_2D_17, ZWORKSHUGRAD_2D_18)
  
  DO JIJ=1, D%NIJT
    ZWORKSHUGRAD_2D_19(JIJ)=(PUM(JIJ, IKB+1)-PUM(JIJ, IKB)) /                   &
                ZWORKSHUGRAD_2D_9(JIJ) * PDZY(JIJ, IKB+1)
  END DO
  
  CALL MXF2D_PHY(D, ZWORKSHUGRAD_2D_19, ZWORKSHUGRAD_2D_20)
  CALL MYF2D_PHY(D, ZWORKSHUGRAD_2D_12, ZWORKSHUGRAD_2D_21)
  CALL MXF2D_PHY(D, ZWORKSHUGRAD_2D_14, ZWORKSHUGRAD_2D_22)
  
  DO JIJ=1, D%NIJT
    ZWORKSHUGRAD_2D_23(JIJ)=(PVM(JIJ, IKB+1)-PVM(JIJ, IKB)) /                   &
                ZWORKSHUGRAD_2D_10(JIJ) * PDZX(JIJ, IKB+1)
  END DO
  
  CALL MYF2D_PHY(D, ZWORKSHUGRAD_2D_23, ZWORKSHUGRAD_2D_24)
  CALL MXF2D_PHY(D, ZWORKSHUGRAD_2D_16, ZWORKSHUGRAD_2D_25)
  CALL MXF2D_PHY(D, ZWORKSHUGRAD_2D_18, ZWORKSHUGRAD_2D_26)
  
  DO JIJ=1, D%NIJT
    ZWORK(JIJ, IKB) =  -                                             &
       ZWORKSHUGRAD_2D_26(JIJ)   &
     *(ZWORKSHUGRAD_2D_22(JIJ)                                                       &  
      -ZWORKSHUGRAD_2D_20(JIJ) / ZWORKSHUGRAD_2D_21(JIJ)&
      -ZWORKSHUGRAD_2D_24(JIJ) / ZWORKSHUGRAD_2D_25(JIJ)&
      )
  END DO
  
  !
  ! dynamic production 
  
  DO JK=1, IKT
    DO JIJ=1, IIJT
      PDP(JIJ, JK) = PDP(JIJ, JK) + ZWORK(JIJ, JK)
    END DO
  END DO
  
  ! 
END IF
!
! Storage in the LES configuration
!
IF (TLES%LLES_CALL .AND. KSPLT==1) THEN
  CALL SECOND_MNH(ZTIME1)
  CALL MYF_PHY(D, ZFLX, ZWORKSHUGRAD_3D_1)
  CALL MXF_PHY(D, ZWORKSHUGRAD_3D_1, ZWORKSHUGRAD_3D_2)
  CALL LES_MEAN_SUBGRID_PHY(D, TLES, ZWORKSHUGRAD_3D_2(1:D%NIJT, 1:D%NKT), TLES%X_LES_SUBGRID_UV ) 
  CALL GY_U_UV_PHY(D, OFLAT, PUM, PDYY, PDZZ, PDZY, ZWORKSHUGRAD_3D_1)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      ZWORKSHUGRAD_3D_2(JIJ, JK)=ZWORKSHUGRAD_3D_1(JIJ, JK)*ZFLX(JIJ, JK)
    END DO
  END DO
  
  CALL MYF_PHY(D, ZWORKSHUGRAD_3D_2, ZWORKSHUGRAD_3D_3)
  CALL MXF_PHY(D, ZWORKSHUGRAD_3D_3, ZWORKSHUGRAD_3D_4)
  CALL LES_MEAN_SUBGRID_PHY(D, TLES, ZWORKSHUGRAD_3D_4(1:D%NIJT, 1:D%NKT), TLES%X_LES_RES_ddxa_U_SBG_UaU , .TRUE.)
  CALL GX_V_UV_PHY(D, OFLAT, PVM, PDXX, PDZZ, PDZX, ZWORKSHUGRAD_3D_1)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      ZWORKSHUGRAD_3D_2(JIJ, JK)=ZWORKSHUGRAD_3D_1(JIJ, JK)*ZFLX(JIJ, JK)
    END DO
  END DO
  
  CALL MYF_PHY(D, ZWORKSHUGRAD_3D_2, ZWORKSHUGRAD_3D_3)
  CALL MXF_PHY(D, ZWORKSHUGRAD_3D_3, ZWORKSHUGRAD_3D_4)
  CALL LES_MEAN_SUBGRID_PHY(D, TLES, ZWORKSHUGRAD_3D_4(1:D%NIJT, 1:D%NKT), TLES%X_LES_RES_ddxa_V_SBG_UaV , .TRUE.)
  CALL SECOND_MNH(ZTIME2)
  TLES%XTIME_LES = TLES%XTIME_LES + ZTIME2 - ZTIME1
END IF
!
!
END SUBROUTINE TURB_HOR_UV
END MODULE MODE_TURB_HOR_UV
