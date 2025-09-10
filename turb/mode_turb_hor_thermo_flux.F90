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
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZMXM3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZGX_M_U3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT) ::ZMXM2D_WORK1
REAL, DIMENSION(D%NIT,D%NJT) ::ZDXM2D_WORK1
REAL, DIMENSION(D%NIT,D%NJT) ::ZSHUGRADWK1_2D
REAL, DIMENSION(D%NIT,D%NJT) ::ZMXM2D_WORK2
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZSHUGRADWK1_3D
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZDXF3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZDZF3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZSHUGRADWK2_3D
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZMXF3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZSHUGRADWK3_3D
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZMZM3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZMZF3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZMXF3D_WORK2
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZMZM3D_WORK2
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZMYM3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZGY_M_V3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT) ::ZMYM2D_WORK1
REAL, DIMENSION(D%NIT,D%NJT) ::ZDYM2D_WORK1
REAL, DIMENSION(D%NIT,D%NJT) ::ZMYM2D_WORK2
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZDYF3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZMYF3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZMYF3D_WORK2
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZGX_W_UW3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZGX_M_M3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZGY_W_VW3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZGY_M_M3D_WORK1
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

ZCOEFF(:,:,IKB+2)= - PDZZ(:,:,IKB+1) /      &
       ( (PDZZ(:,:,IKB+2)+PDZZ(:,:,IKB+1)) * PDZZ(:,:,IKB+2) )
ZCOEFF(:,:,IKB+1)=   (PDZZ(:,:,IKB+2)+PDZZ(:,:,IKB+1)) /      &
       ( PDZZ(:,:,IKB+1) * PDZZ(:,:,IKB+2) )
ZCOEFF(:,:,IKB)= - (PDZZ(:,:,IKB+2)+2.*PDZZ(:,:,IKB+1)) /      &
       ( (PDZZ(:,:,IKB+2)+PDZZ(:,:,IKB+1)) * PDZZ(:,:,IKB+1) )

!
!*       2.   < U' THETA'l >
!             --------------
!
!
#define UTHETA UTHETA_ON
!!!!!!! UTHETA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef UTHETA
CALL MXM_PHY(D, PK(:, :, :), ZMXM3D_WORK1)
CALL GX_M_U_PHY(D, OFLAT,PTHLM(:, :, :),PDXX(:, :, :),PDZZ(:, :, :),PDZX(:, :, :), ZGX_M_U3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZFLX(JI, JJ, JK)     = -TURBN%XCSHF * ZMXM3D_WORK1(JI, JJ, JK) * ZGX_M_U3D_WORK1(JI, JJ, JK)
    END DO
  END DO
END DO

!

ZFLX(:,:,IKE+1) = ZFLX(:,:,IKE) 

!
! Compute the flux at the first inner U-point with an uncentred vertical  
! gradient
CALL MXM2D_PHY(D, PK(:,:,IKB), ZMXM2D_WORK1)

DO JJ=1, IJT
  DO JI=1, IIT
    ZSHUGRADWK1_2D(JI, JJ) = ZCOEFF(JI, JJ, IKB+2)*PTHLM(JI, JJ, IKB+2)        &
             +ZCOEFF(JI, JJ, IKB+1)*PTHLM(JI, JJ, IKB+1)       &
             +ZCOEFF(JI, JJ, IKB)*PTHLM(JI, JJ, IKB)
  END DO
END DO

!
CALL MXM2D_PHY(D, ZSHUGRADWK1_2D, ZMXM2D_WORK2)
CALL DXM2D_PHY(D, PTHLM(:,:,IKB), ZDXM2D_WORK1)

DO JJ=1, IJT
  DO JI=1, IIT
    ZFLX(JI, JJ, IKB) = -TURBN%XCSHF * ZMXM2D_WORK1(JI, JJ) *          &
      ( ZDXM2D_WORK1(JI, JJ) * PINV_PDXX(JI, JJ, IKB)           &
       -ZMXM2D_WORK2(JI, JJ)      &
            *0.5* ( PDZX(JI, JJ, IKB+1)+PDZX(JI, JJ, IKB))       &
            * PINV_PDXX(JI, JJ, IKB) )
  END DO
END DO

!
! extrapolates the flux under the ground so that the vertical average with 
! the IKB flux gives the ground value  ( warning the tangential surface
! flux has been set to 0 for the moment !!  to be improved )

DO JJ=1, IJT
  DO JI=1, IIT
    ZSHUGRADWK1_2D(JI, JJ) = PSFTHM(JI, JJ)* PDIRCOSXW(JI, JJ)
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
! Add this source to the Theta_l sources
!
IF (.NOT. OFLAT) THEN
  CALL MXM_PHY(D, PRHODJ(:, :, :), ZMXM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMXM3D_WORK1(JI, JJ, JK) * ZFLX(JI, JJ, JK) * PINV_PDXX(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL DXF_PHY(D, ZSHUGRADWK1_3D, ZDXF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK3_3D(JI, JJ, JK) = ZFLX(JI, JJ, JK) * PINV_PDXX(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MZM_PHY(D, ZSHUGRADWK3_3D, ZMZM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = PDZX(JI, JJ, JK)*(ZMZM3D_WORK1(JI, JJ, JK))
    END DO
  END DO
END DO

!
CALL MXF_PHY(D, ZSHUGRADWK2_3D, ZMXF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = PMZM_PRHODJ(JI, JJ, JK) *ZMXF3D_WORK1(JI, JJ, JK) * PINV_PDZZ(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL DZF_PHY(D, ZSHUGRADWK1_3D, ZDZF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      PRTHLS(JI, JJ, JK) =  PRTHLS(JI, JJ, JK)                                                   &
                      - ZDXF3D_WORK1(JI, JJ, JK)                          &
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
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMXM3D_WORK1(JI, JJ, JK) * ZFLX(JI, JJ, JK) * PINV_PDXX(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL DXF_PHY(D, ZSHUGRADWK1_3D, ZDXF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      PRTHLS(JI, JJ, JK) =  PRTHLS(JI, JJ, JK) - ZDXF3D_WORK1(JI, JJ, JK)
    END DO
  END DO
END DO

!
END IF
!
! Compute the equivalent tendancy for Rc and Ri
!
IF ( KRRL >= 1 ) THEN
  IF (.NOT. OFLAT) THEN
    
DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = PRHODJ(JI, JJ, JK)*PATHETA(JI, JJ, JK)*PSRCM(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MXM_PHY(D, ZSHUGRADWK2_3D, ZMXM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMXM3D_WORK1(JI, JJ, JK)*ZFLX(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MXF_PHY(D, ZSHUGRADWK1_3D, ZMXF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = PRHODJ(JI, JJ, JK)*PATHETA(JI, JJ, JK)*PSRCM(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MZM_PHY(D, ZSHUGRADWK2_3D, ZMZM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = ZFLX(JI, JJ, JK)*PINV_PDXX(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MZM_PHY(D, ZSHUGRADWK2_3D, ZMZM3D_WORK2)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = PDZX(JI, JJ, JK)*(ZMZM3D_WORK2(JI, JJ, JK))
    END DO
  END DO
END DO

!
CALL MXF_PHY(D, ZSHUGRADWK1_3D, ZMXF3D_WORK2)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMZM3D_WORK1(JI, JJ, JK)*ZMXF3D_WORK2(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MZF_PHY(D, ZSHUGRADWK1_3D, ZMZF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZFLXC(JI, JJ, JK) = 2.*( ZMXF3D_WORK1(JI, JJ, JK)                       &
                      +ZMZF3D_WORK1(JI, JJ, JK)&
                     )    
    END DO
  END DO
END DO

!
IF ( KRRI >= 1 ) THEN
      
DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = PRHODJ(JI, JJ, JK)*PATHETA(JI, JJ, JK)*PSRCM(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MXM_PHY(D, ZSHUGRADWK2_3D, ZMXM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMXM3D_WORK1(JI, JJ, JK)*ZFLX(JI, JJ, JK)*PINV_PDXX(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL DXF_PHY(D, ZSHUGRADWK1_3D, ZDXF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = PRHODJ(JI, JJ, JK)*PATHETA(JI, JJ, JK)*PSRCM(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MZM_PHY(D, ZSHUGRADWK2_3D, ZMZM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = ZFLX(JI, JJ, JK)*PINV_PDXX(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MZM_PHY(D, ZSHUGRADWK2_3D, ZMZM3D_WORK2)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = PDZX(JI, JJ, JK)*(ZMZM3D_WORK2(JI, JJ, JK))
    END DO
  END DO
END DO

!
CALL MXF_PHY(D, ZSHUGRADWK1_3D, ZMXF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMZM3D_WORK1(JI, JJ, JK)*ZMXF3D_WORK1(JI, JJ, JK)&
                                                 *PINV_PDZZ(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL DZF_PHY(D, ZSHUGRADWK1_3D, ZDZF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      PRRS(JI, JJ, JK, 2) = PRRS(JI, JJ, JK, 2) +  2. *                                    &
              (- ZDXF3D_WORK1(JI, JJ, JK)                   &
               + ZDZF3D_WORK1(JI, JJ, JK)                        &
              )*(1.0-PFRAC_ICE(JI, JJ, JK))      
    END DO
  END DO
END DO

!

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = PRHODJ(JI, JJ, JK)*PATHETA(JI, JJ, JK)*PSRCM(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MXM_PHY(D, ZSHUGRADWK2_3D, ZMXM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMXM3D_WORK1(JI, JJ, JK)*ZFLX(JI, JJ, JK)*PINV_PDXX(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL DXF_PHY(D, ZSHUGRADWK1_3D, ZDXF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = PRHODJ(JI, JJ, JK)*PATHETA(JI, JJ, JK)*PSRCM(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MZM_PHY(D, ZSHUGRADWK2_3D, ZMZM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = ZFLX(JI, JJ, JK)*PINV_PDXX(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MZM_PHY(D, ZSHUGRADWK2_3D, ZMZM3D_WORK2)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = PDZX(JI, JJ, JK)*(ZMZM3D_WORK2(JI, JJ, JK))
    END DO
  END DO
END DO

!
CALL MXF_PHY(D, ZSHUGRADWK1_3D, ZMXF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMZM3D_WORK1(JI, JJ, JK)*ZMXF3D_WORK1(JI, JJ, JK)&
                                                 *PINV_PDZZ(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL DZF_PHY(D, ZSHUGRADWK1_3D, ZDZF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      PRRS(JI, JJ, JK, 4) = PRRS(JI, JJ, JK, 4) +  2. *                                    &
              (- ZDXF3D_WORK1(JI, JJ, JK)                   &
               + ZDZF3D_WORK1(JI, JJ, JK)                        &
              )*PFRAC_ICE(JI, JJ, JK)    
    END DO
  END DO
END DO

!
ELSE
      
DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = PRHODJ(JI, JJ, JK)*PATHETA(JI, JJ, JK)*PSRCM(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MXM_PHY(D, ZSHUGRADWK2_3D, ZMXM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMXM3D_WORK1(JI, JJ, JK)*ZFLX(JI, JJ, JK)*PINV_PDXX(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL DXF_PHY(D, ZSHUGRADWK1_3D, ZDXF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = PRHODJ(JI, JJ, JK)*PATHETA(JI, JJ, JK)*PSRCM(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MZM_PHY(D, ZSHUGRADWK2_3D, ZMZM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = ZFLX(JI, JJ, JK)*PINV_PDXX(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MZM_PHY(D, ZSHUGRADWK2_3D, ZMZM3D_WORK2)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = PDZX(JI, JJ, JK)*(ZMZM3D_WORK2(JI, JJ, JK))
    END DO
  END DO
END DO

!
CALL MXF_PHY(D, ZSHUGRADWK1_3D, ZMXF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMZM3D_WORK1(JI, JJ, JK)*ZMXF3D_WORK1(JI, JJ, JK)&
                                                 *PINV_PDZZ(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL DZF_PHY(D, ZSHUGRADWK1_3D, ZDZF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      PRRS(JI, JJ, JK, 2) = PRRS(JI, JJ, JK, 2) +  2. *                                    &
              (- ZDXF3D_WORK1(JI, JJ, JK)                   &
               + ZDZF3D_WORK1(JI, JJ, JK)                        &
              )    
    END DO
  END DO
END DO

!
END IF
  ELSE
    
DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = PRHODJ(JI, JJ, JK)*PATHETA(JI, JJ, JK)*PSRCM(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MXM_PHY(D, ZSHUGRADWK2_3D, ZMXM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMXM3D_WORK1(JI, JJ, JK)*ZFLX(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MXF_PHY(D, ZSHUGRADWK1_3D, ZMXF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZFLXC(JI, JJ, JK) = 2.*ZMXF3D_WORK1(JI, JJ, JK)    
    END DO
  END DO
END DO

!
IF ( KRRI >= 1 ) THEN
      
DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = PRHODJ(JI, JJ, JK)*PATHETA(JI, JJ, JK)*PSRCM(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MXM_PHY(D, ZSHUGRADWK2_3D, ZMXM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMXM3D_WORK1(JI, JJ, JK)*ZFLX(JI, JJ, JK)*PINV_PDXX(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL DXF_PHY(D, ZSHUGRADWK1_3D, ZDXF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      PRRS(JI, JJ, JK, 2) = PRRS(JI, JJ, JK, 2) -  2. *                                    &
              ZDXF3D_WORK1(JI, JJ, JK)*(1.0-PFRAC_ICE(JI, JJ, JK))      
    END DO
  END DO
END DO

!

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = PRHODJ(JI, JJ, JK)*PATHETA(JI, JJ, JK)*PSRCM(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MXM_PHY(D, ZSHUGRADWK2_3D, ZMXM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMXM3D_WORK1(JI, JJ, JK)*ZFLX(JI, JJ, JK)*PINV_PDXX(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL DXF_PHY(D, ZSHUGRADWK1_3D, ZDXF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      PRRS(JI, JJ, JK, 4) = PRRS(JI, JJ, JK, 4) -  2. *                                    &
              ZDXF3D_WORK1(JI, JJ, JK)*PFRAC_ICE(JI, JJ, JK)    
    END DO
  END DO
END DO

!
ELSE
      
DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = PRHODJ(JI, JJ, JK)*PATHETA(JI, JJ, JK)*PSRCM(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MXM_PHY(D, ZSHUGRADWK2_3D, ZMXM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMXM3D_WORK1(JI, JJ, JK)*ZFLX(JI, JJ, JK)*PINV_PDXX(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL DXF_PHY(D, ZSHUGRADWK1_3D, ZDXF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      PRRS(JI, JJ, JK, 2) = PRRS(JI, JJ, JK, 2) -  2. *                                    &
              ZDXF3D_WORK1(JI, JJ, JK)    
    END DO
  END DO
END DO

!
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

  CALL IO_FIELD_WRITE(TPFILE,TZFIELD,ZFLX)
END IF
!
IF (KSPLT==1 .AND. TLES%LLES_CALL) THEN
  CALL SECOND_MNH(ZTIME1)
  CALL MXF_PHY(D, ZFLX(:, :, :), ZMXF3D_WORK1)
CALL LES_MEAN_SUBGRID(ZMXF3D_WORK1(:, :, :), TLES%X_LES_SUBGRID_UThl )   
CALL GX_W_UW_PHY(D, OFLAT, PWM(:, :, :),PDXX(:, :, :),PDZZ(:, :, :),PDZX(:, :, :), ZGX_W_UW3D_WORK1)
CALL MZM_PHY(D, ZFLX(:, :, :), ZMZM3D_WORK1)

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
CALL LES_MEAN_SUBGRID(ZMZF3D_WORK1(:, :, :),&
                         TLES%X_LES_RES_ddxa_W_SBG_UaThl , .TRUE. )  
CALL GX_M_M_PHY(D, OFLAT,PTHLM(:, :, :),PDXX(:, :, :),PDZZ(:, :, :),PDZX(:, :, :), ZGX_M_M3D_WORK1)
CALL MXF_PHY(D, ZFLX(:, :, :), ZMXF3D_WORK1)
CALL LES_MEAN_SUBGRID(ZGX_M_M3D_WORK1(:, :, :)*ZMXF3D_WORK1(:, :, :),&
                         TLES%X_LES_RES_ddxa_Thl_SBG_UaThl , .TRUE. )  
IF (KRR>=1) THEN
    CALL GX_M_M_PHY(D, OFLAT,PRM(:,:,:,1),PDXX(:, :, :),PDZZ(:, :, :),PDZX(:, :, :), ZGX_M_M3D_WORK1)
CALL MXF_PHY(D, ZFLX(:, :, :), ZMXF3D_WORK1)
CALL LES_MEAN_SUBGRID(ZGX_M_M3D_WORK1(:, :, :)*ZMXF3D_WORK1(:, :, :), &
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
  CALL MXM_PHY(D, PK(:, :, :), ZMXM3D_WORK1)
CALL GX_M_U_PHY(D, OFLAT,PRM(:,:,:,1),PDXX(:, :, :),PDZZ(:, :, :),PDZX(:, :, :), ZGX_M_U3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZFLX(JI, JJ, JK)     = -TURBN%XCHF * ZMXM3D_WORK1(JI, JJ, JK) * ZGX_M_U3D_WORK1(JI, JJ, JK)
    END DO
  END DO
END DO

!

  ZFLX(:,:,IKE+1) = ZFLX(:,:,IKE)

!
! Compute the flux at the first inner U-point with an uncentred vertical  
! gradient
  CALL MXM2D_PHY(D, PK(:,:,IKB), ZMXM2D_WORK1)
CALL DXM2D_PHY(D, PRM(:,:,IKB,1), ZDXM2D_WORK1)

DO JJ=1, IJT
  DO JI=1, IIT
    ZSHUGRADWK1_2D(JI, JJ) = ZCOEFF(JI, JJ, IKB+2)*PRM(JI, JJ, IKB+2, 1)        &
               +ZCOEFF(JI, JJ, IKB+1)*PRM(JI, JJ, IKB+1, 1)       &
               +ZCOEFF(JI, JJ, IKB)*PRM(JI, JJ, IKB, 1)
  END DO
END DO

!
CALL MXM2D_PHY(D, ZSHUGRADWK1_2D, ZMXM2D_WORK2)

DO JJ=1, IJT
  DO JI=1, IIT
    ZFLX(JI, JJ, IKB) = -TURBN%XCHF * ZMXM2D_WORK1(JI, JJ) *           &
        ( ZDXM2D_WORK1(JI, JJ) * PINV_PDXX(JI, JJ, IKB)           &
         -ZMXM2D_WORK2(JI, JJ)      &
              *0.5* ( PDZX(JI, JJ, IKB+1)+PDZX(JI, JJ, IKB))       &
              * PINV_PDXX(JI, JJ, IKB) )
  END DO
END DO

!
! extrapolates the flux under the ground so that the vertical average with 
! the IKB flux gives the ground value  ( warning the tangential surface
! flux has been set to 0 for the moment !!  to be improved )
  
DO JJ=1, IJT
  DO JI=1, IIT
    ZSHUGRADWK1_2D(JI, JJ) = PSFRM(JI, JJ)* PDIRCOSXW(JI, JJ)
  END DO
END DO

!
CALL MXM2D_PHY(D, ZSHUGRADWK1_2D, ZMXM2D_WORK1)

DO JJ=1, IJT
  DO JI=1, IIT
    ZFLX(JI, JJ, IKB-1) = 2. * ZMXM2D_WORK1(JI, JJ) &
                           - ZFLX(JI, JJ, IKB)  
  END DO
END DO

!
!
  ! Add this source to the conservative mixing ratio sources
  !
  IF (.NOT. OFLAT) THEN
    CALL MXM_PHY(D, PRHODJ(:, :, :), ZMXM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMXM3D_WORK1(JI, JJ, JK) * ZFLX(JI, JJ, JK) * PINV_PDXX(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL DXF_PHY(D, ZSHUGRADWK1_3D, ZDXF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK3_3D(JI, JJ, JK) = ZFLX(JI, JJ, JK) * PINV_PDXX(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MZM_PHY(D, ZSHUGRADWK3_3D, ZMZM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = PDZX(JI, JJ, JK)*(ZMZM3D_WORK1(JI, JJ, JK))
    END DO
  END DO
END DO

!
CALL MXF_PHY(D, ZSHUGRADWK2_3D, ZMXF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = PMZM_PRHODJ(JI, JJ, JK) *ZMXF3D_WORK1(JI, JJ, JK) * PINV_PDZZ(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL DZF_PHY(D, ZSHUGRADWK1_3D, ZDZF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      PRRS(JI, JJ, JK, 1) = PRRS(JI, JJ, JK, 1)                                             &
                        - ZDXF3D_WORK1(JI, JJ, JK)                          &
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
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMXM3D_WORK1(JI, JJ, JK) * ZFLX(JI, JJ, JK) * PINV_PDXX(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL DXF_PHY(D, ZSHUGRADWK1_3D, ZDXF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      PRRS(JI, JJ, JK, 1) = PRRS(JI, JJ, JK, 1) - ZDXF3D_WORK1(JI, JJ, JK)  
    END DO
  END DO
END DO

!
END IF
  !
  ! Compute the equivalent tendancy for Rc and Ri
  !
  IF ( KRRL >= 1 ) THEN
    IF (.NOT. OFLAT) THEN
      
DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = PRHODJ(JI, JJ, JK)*PAMOIST(JI, JJ, JK)*PSRCM(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MXM_PHY(D, ZSHUGRADWK2_3D, ZMXM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMXM3D_WORK1(JI, JJ, JK)*ZFLX(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MXF_PHY(D, ZSHUGRADWK1_3D, ZMXF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = PRHODJ(JI, JJ, JK)*PAMOIST(JI, JJ, JK)*PSRCM(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MZM_PHY(D, ZSHUGRADWK2_3D, ZMZM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = ZFLX(JI, JJ, JK)*PINV_PDXX(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MZM_PHY(D, ZSHUGRADWK2_3D, ZMZM3D_WORK2)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = PDZX(JI, JJ, JK)*(ZMZM3D_WORK2(JI, JJ, JK))
    END DO
  END DO
END DO

!
CALL MXF_PHY(D, ZSHUGRADWK1_3D, ZMXF3D_WORK2)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMZM3D_WORK1(JI, JJ, JK)*ZMXF3D_WORK2(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MZF_PHY(D, ZSHUGRADWK1_3D, ZMZF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZFLXC(JI, JJ, JK) = ZFLXC(JI, JJ, JK)            &
                  + 2.*( ZMXF3D_WORK1(JI, JJ, JK)                     &
                        +ZMZF3D_WORK1(JI, JJ, JK)&
                       )      
    END DO
  END DO
END DO

!
IF ( KRRI >= 1 ) THEN
        
DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = PRHODJ(JI, JJ, JK)*PAMOIST(JI, JJ, JK)*PSRCM(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MXM_PHY(D, ZSHUGRADWK2_3D, ZMXM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMXM3D_WORK1(JI, JJ, JK)*ZFLX(JI, JJ, JK)*PINV_PDXX(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL DXF_PHY(D, ZSHUGRADWK1_3D, ZDXF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = PRHODJ(JI, JJ, JK)*PAMOIST(JI, JJ, JK)*PSRCM(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MZM_PHY(D, ZSHUGRADWK2_3D, ZMZM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = ZFLX(JI, JJ, JK)*PINV_PDXX(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MZM_PHY(D, ZSHUGRADWK2_3D, ZMZM3D_WORK2)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = PDZX(JI, JJ, JK)*(ZMZM3D_WORK2(JI, JJ, JK))
    END DO
  END DO
END DO

!
CALL MXF_PHY(D, ZSHUGRADWK1_3D, ZMXF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMZM3D_WORK1(JI, JJ, JK)*ZMXF3D_WORK1(JI, JJ, JK)&
                                                 *PINV_PDZZ(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL DZF_PHY(D, ZSHUGRADWK1_3D, ZDZF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      PRRS(JI, JJ, JK, 2) = PRRS(JI, JJ, JK, 2) +  2. *                                  &
              (- ZDXF3D_WORK1(JI, JJ, JK)                   &
               + ZDZF3D_WORK1(JI, JJ, JK)                        &
              )*(1.0-PFRAC_ICE(JI, JJ, JK))        
    END DO
  END DO
END DO

!

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = PRHODJ(JI, JJ, JK)*PAMOIST(JI, JJ, JK)*PSRCM(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MXM_PHY(D, ZSHUGRADWK2_3D, ZMXM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMXM3D_WORK1(JI, JJ, JK)*ZFLX(JI, JJ, JK)*PINV_PDXX(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL DXF_PHY(D, ZSHUGRADWK1_3D, ZDXF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = PRHODJ(JI, JJ, JK)*PAMOIST(JI, JJ, JK)*PSRCM(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MZM_PHY(D, ZSHUGRADWK2_3D, ZMZM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = ZFLX(JI, JJ, JK)*PINV_PDXX(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MZM_PHY(D, ZSHUGRADWK2_3D, ZMZM3D_WORK2)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = PDZX(JI, JJ, JK)*(ZMZM3D_WORK2(JI, JJ, JK))
    END DO
  END DO
END DO

!
CALL MXF_PHY(D, ZSHUGRADWK1_3D, ZMXF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMZM3D_WORK1(JI, JJ, JK)*ZMXF3D_WORK1(JI, JJ, JK)&
                                                 *PINV_PDZZ(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL DZF_PHY(D, ZSHUGRADWK1_3D, ZDZF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      PRRS(JI, JJ, JK, 2) = PRRS(JI, JJ, JK, 2) +  2. *                                  &
              (- ZDXF3D_WORK1(JI, JJ, JK)                   &
               + ZDZF3D_WORK1(JI, JJ, JK)                        &
              )*PFRAC_ICE(JI, JJ, JK)      
    END DO
  END DO
END DO

!
ELSE
        
DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = PRHODJ(JI, JJ, JK)*PAMOIST(JI, JJ, JK)*PSRCM(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MXM_PHY(D, ZSHUGRADWK2_3D, ZMXM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMXM3D_WORK1(JI, JJ, JK)*ZFLX(JI, JJ, JK)*PINV_PDXX(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL DXF_PHY(D, ZSHUGRADWK1_3D, ZDXF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = PRHODJ(JI, JJ, JK)*PAMOIST(JI, JJ, JK)*PSRCM(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MZM_PHY(D, ZSHUGRADWK2_3D, ZMZM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = ZFLX(JI, JJ, JK)*PINV_PDXX(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MZM_PHY(D, ZSHUGRADWK2_3D, ZMZM3D_WORK2)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = PDZX(JI, JJ, JK)*(ZMZM3D_WORK2(JI, JJ, JK))
    END DO
  END DO
END DO

!
CALL MXF_PHY(D, ZSHUGRADWK1_3D, ZMXF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMZM3D_WORK1(JI, JJ, JK)*ZMXF3D_WORK1(JI, JJ, JK)&
                                                 *PINV_PDZZ(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL DZF_PHY(D, ZSHUGRADWK1_3D, ZDZF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      PRRS(JI, JJ, JK, 2) = PRRS(JI, JJ, JK, 2) +  2. *                                  &
              (- ZDXF3D_WORK1(JI, JJ, JK)                   &
               + ZDZF3D_WORK1(JI, JJ, JK)                        &
              )      
    END DO
  END DO
END DO

!
END IF
    ELSE
      
DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = PRHODJ(JI, JJ, JK)*PAMOIST(JI, JJ, JK)*PSRCM(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MXM_PHY(D, ZSHUGRADWK2_3D, ZMXM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMXM3D_WORK1(JI, JJ, JK)*ZFLX(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MXF_PHY(D, ZSHUGRADWK1_3D, ZMXF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZFLXC(JI, JJ, JK) = ZFLXC(JI, JJ, JK) + 2.*ZMXF3D_WORK1(JI, JJ, JK)      
    END DO
  END DO
END DO

!
IF ( KRRI >= 1 ) THEN
        
DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = PRHODJ(JI, JJ, JK)*PAMOIST(JI, JJ, JK)*PSRCM(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MXM_PHY(D, ZSHUGRADWK2_3D, ZMXM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMXM3D_WORK1(JI, JJ, JK)*ZFLX(JI, JJ, JK)*PINV_PDXX(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL DXF_PHY(D, ZSHUGRADWK1_3D, ZDXF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      PRRS(JI, JJ, JK, 2) = PRRS(JI, JJ, JK, 2) -  2. *                                  &
              ZDXF3D_WORK1(JI, JJ, JK)*(1.0-PFRAC_ICE(JI, JJ, JK))        
    END DO
  END DO
END DO

!

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = PRHODJ(JI, JJ, JK)*PAMOIST(JI, JJ, JK)*PSRCM(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MXM_PHY(D, ZSHUGRADWK2_3D, ZMXM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMXM3D_WORK1(JI, JJ, JK)*ZFLX(JI, JJ, JK)*PINV_PDXX(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL DXF_PHY(D, ZSHUGRADWK1_3D, ZDXF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      PRRS(JI, JJ, JK, 4) = PRRS(JI, JJ, JK, 4) -  2. *                                  &
              ZDXF3D_WORK1(JI, JJ, JK)*PFRAC_ICE(JI, JJ, JK)      
    END DO
  END DO
END DO

!
ELSE
        
DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = PRHODJ(JI, JJ, JK)*PAMOIST(JI, JJ, JK)*PSRCM(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MXM_PHY(D, ZSHUGRADWK2_3D, ZMXM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMXM3D_WORK1(JI, JJ, JK)*ZFLX(JI, JJ, JK)*PINV_PDXX(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL DXF_PHY(D, ZSHUGRADWK1_3D, ZDXF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      PRRS(JI, JJ, JK, 2) = PRRS(JI, JJ, JK, 2) -  2. *                                  &
              ZDXF3D_WORK1(JI, JJ, JK)      
    END DO
  END DO
END DO

!
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

    CALL IO_FIELD_WRITE(TPFILE,TZFIELD,ZFLX)
  END IF
  !
  IF (KSPLT==1 .AND. TLES%LLES_CALL) THEN
    CALL SECOND_MNH(ZTIME1)
    CALL MXF_PHY(D, ZFLX(:, :, :), ZMXF3D_WORK1)
CALL LES_MEAN_SUBGRID(ZMXF3D_WORK1(:, :, :), TLES%X_LES_SUBGRID_URt )     
CALL GX_W_UW_PHY(D, OFLAT, PWM(:, :, :),PDXX(:, :, :),PDZZ(:, :, :),PDZX(:, :, :), ZGX_W_UW3D_WORK1)
CALL MZM_PHY(D, ZFLX(:, :, :), ZMZM3D_WORK1)

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
CALL LES_MEAN_SUBGRID(ZMZF3D_WORK1(:, :, :),&
                           TLES%X_LES_RES_ddxa_W_SBG_UaRt , .TRUE. )    
CALL GX_M_M_PHY(D, OFLAT,PTHLM(:, :, :),PDXX(:, :, :),PDZZ(:, :, :),PDZX(:, :, :), ZGX_M_M3D_WORK1)
CALL MXF_PHY(D, ZFLX(:, :, :), ZMXF3D_WORK1)
CALL LES_MEAN_SUBGRID(ZGX_M_M3D_WORK1(:, :, :)*ZMXF3D_WORK1(:, :, :),&
                           TLES%X_LES_RES_ddxa_Thl_SBG_UaRt , .TRUE. )    
CALL GX_M_M_PHY(D, OFLAT,PRM(:,:,:,1),PDXX(:, :, :),PDZZ(:, :, :),PDZX(:, :, :), ZGX_M_M3D_WORK1)
CALL MXF_PHY(D, ZFLX(:, :, :), ZMXF3D_WORK1)
CALL LES_MEAN_SUBGRID(ZGX_M_M3D_WORK1(:, :, :)*ZMXF3D_WORK1(:, :, :),&
                           TLES%X_LES_RES_ddxa_Rt_SBG_UaRt , .TRUE. )    
CALL SECOND_MNH(ZTIME2)
    TLES%XTIME_LES = TLES%XTIME_LES + ZTIME2 - ZTIME1
  END IF
!
  !
  IF (KRRL>0 .AND. KSPLT==1 .AND. TLES%LLES_CALL) THEN
    CALL SECOND_MNH(ZTIME1)
    CALL MXF_PHY(D, ZFLXC(:, :, :), ZMXF3D_WORK1)
CALL LES_MEAN_SUBGRID(ZMXF3D_WORK1(:, :, :), TLES%X_LES_SUBGRID_URc )    
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
  CALL MYM_PHY(D, PK(:, :, :), ZMYM3D_WORK1)
CALL GY_M_V_PHY(D, OFLAT,PTHLM(:, :, :),PDYY(:, :, :),PDZZ(:, :, :),PDZY(:, :, :), ZGY_M_V3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZFLX(JI, JJ, JK)     = -TURBN%XCSHF * ZMYM3D_WORK1(JI, JJ, JK) * ZGY_M_V3D_WORK1(JI, JJ, JK)
    END DO
  END DO
END DO

!

  ZFLX(:,:,IKE+1) = ZFLX(:,:,IKE)

ELSE

  ZFLX(:,:,:)     = 0.

END IF
!
!
! Compute the flux at the first inner U-point with an uncentred vertical  
! gradient
CALL MYM2D_PHY(D, PK(:,:,IKB), ZMYM2D_WORK1)
CALL DYM2D_PHY(D, PTHLM(:,:,IKB), ZDYM2D_WORK1)

DO JJ=1, IJT
  DO JI=1, IIT
    ZSHUGRADWK1_2D(JI, JJ) = ZCOEFF(JI, JJ, IKB+2)*PTHLM(JI, JJ, IKB+2)        &
             +ZCOEFF(JI, JJ, IKB+1)*PTHLM(JI, JJ, IKB+1)       &
             +ZCOEFF(JI, JJ, IKB)*PTHLM(JI, JJ, IKB)
  END DO
END DO

!
CALL MYM2D_PHY(D, ZSHUGRADWK1_2D, ZMYM2D_WORK2)

DO JJ=1, IJT
  DO JI=1, IIT
    ZFLX(JI, JJ, IKB) = -TURBN%XCSHF * ZMYM2D_WORK1(JI, JJ) *          &
      ( ZDYM2D_WORK1(JI, JJ) * PINV_PDYY(JI, JJ, IKB)           &
       -ZMYM2D_WORK2(JI, JJ)     &
            *0.5* ( PDZY(JI, JJ, IKB+1)+PDZY(JI, JJ, IKB))       &
            * PINV_PDYY(JI, JJ, IKB) )
  END DO
END DO

!
! extrapolates the flux under the ground so that the vertical average with 
! the IKB flux gives the ground value  ( warning the tangential surface
! flux has been set to 0 for the moment !!  to be improved )

DO JJ=1, IJT
  DO JI=1, IIT
    ZSHUGRADWK1_2D(JI, JJ) = PSFTHM(JI, JJ)* PDIRCOSYW(JI, JJ)
  END DO
END DO

!
CALL MYM2D_PHY(D, ZSHUGRADWK1_2D, ZMYM2D_WORK1)

DO JJ=1, IJT
  DO JI=1, IIT
    ZFLX(JI, JJ, IKB-1) = 2. * ZMYM2D_WORK1(JI, JJ) &
                           - ZFLX(JI, JJ, IKB)
  END DO
END DO

!
!
! Add this source to the Theta_l sources
!
IF (.NOT. O2D) THEN 
  IF (.NOT. OFLAT) THEN
    CALL MYM_PHY(D, PRHODJ(:, :, :), ZMYM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMYM3D_WORK1(JI, JJ, JK) * ZFLX(JI, JJ, JK) * PINV_PDYY(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL DYF_PHY(D, ZSHUGRADWK1_3D, ZDYF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK3_3D(JI, JJ, JK) = ZFLX(JI, JJ, JK) * PINV_PDYY(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MZM_PHY(D, ZSHUGRADWK3_3D, ZMZM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = PDZY(JI, JJ, JK)*(ZMZM3D_WORK1(JI, JJ, JK))
    END DO
  END DO
END DO

!
CALL MYF_PHY(D, ZSHUGRADWK2_3D, ZMYF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = PMZM_PRHODJ(JI, JJ, JK) *ZMYF3D_WORK1(JI, JJ, JK) * PINV_PDZZ(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL DZF_PHY(D, ZSHUGRADWK1_3D, ZDZF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      PRTHLS(JI, JJ, JK) =  PRTHLS(JI, JJ, JK)                                                         &
                        - ZDYF3D_WORK1(JI, JJ, JK)                           &
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
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMYM3D_WORK1(JI, JJ, JK) * ZFLX(JI, JJ, JK) * PINV_PDYY(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL DYF_PHY(D, ZSHUGRADWK1_3D, ZDYF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      PRTHLS(JI, JJ, JK) =  PRTHLS(JI, JJ, JK) - ZDYF3D_WORK1(JI, JJ, JK)  
    END DO
  END DO
END DO

!
END IF
END IF
!
! Compute the equivalent tendancy for Rc and Ri
!
!IF ( TURBN%LSUBG_COND .AND. KRRL > 0 .AND. .NOT. O2D) THEN
IF ( KRRL >= 1 .AND. .NOT. O2D) THEN
  IF (.NOT. OFLAT) THEN
    
DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = PRHODJ(JI, JJ, JK)*PATHETA(JI, JJ, JK)*PSRCM(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MYM_PHY(D, ZSHUGRADWK2_3D, ZMYM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMYM3D_WORK1(JI, JJ, JK)*ZFLX(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MYF_PHY(D, ZSHUGRADWK1_3D, ZMYF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = PRHODJ(JI, JJ, JK)*PATHETA(JI, JJ, JK)*PSRCM(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MZM_PHY(D, ZSHUGRADWK2_3D, ZMZM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = ZFLX(JI, JJ, JK)*PINV_PDYY(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MZM_PHY(D, ZSHUGRADWK2_3D, ZMZM3D_WORK2)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = PDZY(JI, JJ, JK)*(ZMZM3D_WORK2(JI, JJ, JK))
    END DO
  END DO
END DO

!
CALL MYF_PHY(D, ZSHUGRADWK1_3D, ZMYF3D_WORK2)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMZM3D_WORK1(JI, JJ, JK)*ZMYF3D_WORK2(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MZF_PHY(D, ZSHUGRADWK1_3D, ZMZF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZFLXC(JI, JJ, JK) = 2.*( ZMYF3D_WORK1(JI, JJ, JK)                       &
                      +ZMZF3D_WORK1(JI, JJ, JK)&
                     )    
    END DO
  END DO
END DO

!
IF ( KRRI >= 1 ) THEN
      
DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = PRHODJ(JI, JJ, JK)*PATHETA(JI, JJ, JK)*PSRCM(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MYM_PHY(D, ZSHUGRADWK2_3D, ZMYM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMYM3D_WORK1(JI, JJ, JK)*ZFLX(JI, JJ, JK)*PINV_PDYY(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL DYF_PHY(D, ZSHUGRADWK1_3D, ZDYF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = PRHODJ(JI, JJ, JK)*PATHETA(JI, JJ, JK)*PSRCM(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MZM_PHY(D, ZSHUGRADWK2_3D, ZMZM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = ZFLX(JI, JJ, JK)*PINV_PDYY(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MZM_PHY(D, ZSHUGRADWK2_3D, ZMZM3D_WORK2)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = PDZY(JI, JJ, JK)*(ZMZM3D_WORK2(JI, JJ, JK))
    END DO
  END DO
END DO

!
CALL MYF_PHY(D, ZSHUGRADWK1_3D, ZMYF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMZM3D_WORK1(JI, JJ, JK)*ZMYF3D_WORK1(JI, JJ, JK)&
                                                 *PINV_PDZZ(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL DZF_PHY(D, ZSHUGRADWK1_3D, ZDZF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      PRRS(JI, JJ, JK, 2) = PRRS(JI, JJ, JK, 2) + 2. *                                     &
              (- ZDYF3D_WORK1(JI, JJ, JK)                   &
               + ZDZF3D_WORK1(JI, JJ, JK)                        &
              )*(1.0-PFRAC_ICE(JI, JJ, JK))      
    END DO
  END DO
END DO

!

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = PRHODJ(JI, JJ, JK)*PATHETA(JI, JJ, JK)*PSRCM(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MYM_PHY(D, ZSHUGRADWK2_3D, ZMYM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMYM3D_WORK1(JI, JJ, JK)*ZFLX(JI, JJ, JK)*PINV_PDYY(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL DYF_PHY(D, ZSHUGRADWK1_3D, ZDYF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = PRHODJ(JI, JJ, JK)*PATHETA(JI, JJ, JK)*PSRCM(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MZM_PHY(D, ZSHUGRADWK2_3D, ZMZM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = ZFLX(JI, JJ, JK)*PINV_PDYY(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MZM_PHY(D, ZSHUGRADWK2_3D, ZMZM3D_WORK2)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = PDZY(JI, JJ, JK)*(ZMZM3D_WORK2(JI, JJ, JK))
    END DO
  END DO
END DO

!
CALL MYF_PHY(D, ZSHUGRADWK1_3D, ZMYF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMZM3D_WORK1(JI, JJ, JK)*ZMYF3D_WORK1(JI, JJ, JK)&
                                                 *PINV_PDZZ(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL DZF_PHY(D, ZSHUGRADWK1_3D, ZDZF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      PRRS(JI, JJ, JK, 4) = PRRS(JI, JJ, JK, 4) + 2. *                                     &
              (- ZDYF3D_WORK1(JI, JJ, JK)                   &
               + ZDZF3D_WORK1(JI, JJ, JK)                        &
              )*PFRAC_ICE(JI, JJ, JK)    
    END DO
  END DO
END DO

!
ELSE
      
DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = PRHODJ(JI, JJ, JK)*PATHETA(JI, JJ, JK)*PSRCM(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MYM_PHY(D, ZSHUGRADWK2_3D, ZMYM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMYM3D_WORK1(JI, JJ, JK)*ZFLX(JI, JJ, JK)*PINV_PDYY(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL DYF_PHY(D, ZSHUGRADWK1_3D, ZDYF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = PRHODJ(JI, JJ, JK)*PATHETA(JI, JJ, JK)*PSRCM(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MZM_PHY(D, ZSHUGRADWK2_3D, ZMZM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = ZFLX(JI, JJ, JK)*PINV_PDYY(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MZM_PHY(D, ZSHUGRADWK2_3D, ZMZM3D_WORK2)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = PDZY(JI, JJ, JK)*(ZMZM3D_WORK2(JI, JJ, JK))
    END DO
  END DO
END DO

!
CALL MYF_PHY(D, ZSHUGRADWK1_3D, ZMYF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMZM3D_WORK1(JI, JJ, JK)*ZMYF3D_WORK1(JI, JJ, JK)&
                                                 *PINV_PDZZ(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL DZF_PHY(D, ZSHUGRADWK1_3D, ZDZF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      PRRS(JI, JJ, JK, 2) = PRRS(JI, JJ, JK, 2) + 2. *                                     &
              (- ZDYF3D_WORK1(JI, JJ, JK)                   &
               + ZDZF3D_WORK1(JI, JJ, JK)                        &
              )    
    END DO
  END DO
END DO

!
END IF
  ELSE
    
DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = PRHODJ(JI, JJ, JK)*PATHETA(JI, JJ, JK)*PSRCM(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MYM_PHY(D, ZSHUGRADWK2_3D, ZMYM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMYM3D_WORK1(JI, JJ, JK)*ZFLX(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MYF_PHY(D, ZSHUGRADWK1_3D, ZMYF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZFLXC(JI, JJ, JK) = 2.*ZMYF3D_WORK1(JI, JJ, JK)    
    END DO
  END DO
END DO

!
IF ( KRRI >= 1 ) THEN
      
DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = PRHODJ(JI, JJ, JK)*PATHETA(JI, JJ, JK)*PSRCM(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MYM_PHY(D, ZSHUGRADWK2_3D, ZMYM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMYM3D_WORK1(JI, JJ, JK)*ZFLX(JI, JJ, JK)*PINV_PDYY(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL DYF_PHY(D, ZSHUGRADWK1_3D, ZDYF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      PRRS(JI, JJ, JK, 2) = PRRS(JI, JJ, JK, 2) - 2. *                                     &
              ZDYF3D_WORK1(JI, JJ, JK)*(1.0-PFRAC_ICE(JI, JJ, JK))      
    END DO
  END DO
END DO

!

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = PRHODJ(JI, JJ, JK)*PATHETA(JI, JJ, JK)*PSRCM(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MYM_PHY(D, ZSHUGRADWK2_3D, ZMYM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMYM3D_WORK1(JI, JJ, JK)*ZFLX(JI, JJ, JK)*PINV_PDYY(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL DYF_PHY(D, ZSHUGRADWK1_3D, ZDYF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      PRRS(JI, JJ, JK, 4) = PRRS(JI, JJ, JK, 4) - 2. *                                     &
              ZDYF3D_WORK1(JI, JJ, JK)*PFRAC_ICE(JI, JJ, JK)    
    END DO
  END DO
END DO

!
ELSE
      
DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = PRHODJ(JI, JJ, JK)*PATHETA(JI, JJ, JK)*PSRCM(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MYM_PHY(D, ZSHUGRADWK2_3D, ZMYM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMYM3D_WORK1(JI, JJ, JK)*ZFLX(JI, JJ, JK)*PINV_PDYY(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL DYF_PHY(D, ZSHUGRADWK1_3D, ZDYF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      PRRS(JI, JJ, JK, 2) = PRRS(JI, JJ, JK, 2) - 2. *                                     &
              ZDYF3D_WORK1(JI, JJ, JK)    
    END DO
  END DO
END DO

!
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

  CALL IO_FIELD_WRITE(TPFILE,TZFIELD,ZFLX)
END IF
!
IF (KSPLT==1 .AND. TLES%LLES_CALL) THEN
  CALL SECOND_MNH(ZTIME1)
  CALL MYF_PHY(D, ZFLX(:, :, :), ZMYF3D_WORK1)
CALL LES_MEAN_SUBGRID(ZMYF3D_WORK1(:, :, :), TLES%X_LES_SUBGRID_VThl )   
CALL GY_W_VW_PHY(D, OFLAT, PWM(:, :, :),PDYY(:, :, :),PDZZ(:, :, :),PDZY(:, :, :), ZGY_W_VW3D_WORK1)
CALL MZM_PHY(D, ZFLX(:, :, :), ZMZM3D_WORK1)

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
CALL LES_MEAN_SUBGRID(ZMZF3D_WORK1(:, :, :),&
                         TLES%X_LES_RES_ddxa_W_SBG_UaThl , .TRUE. )  
CALL GY_M_M_PHY(D, OFLAT,PTHLM(:, :, :),PDYY(:, :, :),PDZZ(:, :, :),PDZY(:, :, :), ZGY_M_M3D_WORK1)
CALL MYF_PHY(D, ZFLX(:, :, :), ZMYF3D_WORK1)
CALL LES_MEAN_SUBGRID(ZGY_M_M3D_WORK1(:, :, :)*ZMYF3D_WORK1(:, :, :),&
                         TLES%X_LES_RES_ddxa_Thl_SBG_UaThl , .TRUE. )  
IF (KRR>=1) THEN
    CALL GY_M_M_PHY(D, OFLAT,PRM(:,:,:,1),PDYY(:, :, :),PDZZ(:, :, :),PDZY(:, :, :), ZGY_M_M3D_WORK1)
CALL MYF_PHY(D, ZFLX(:, :, :), ZMYF3D_WORK1)
CALL LES_MEAN_SUBGRID(ZGY_M_M3D_WORK1(:, :, :)*ZMYF3D_WORK1(:, :, :),&
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
    CALL MYM_PHY(D, PK(:, :, :), ZMYM3D_WORK1)
CALL GY_M_V_PHY(D, OFLAT,PRM(:,:,:,1),PDYY(:, :, :),PDZZ(:, :, :),PDZY(:, :, :), ZGY_M_V3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZFLX(JI, JJ, JK)     = -TURBN%XCHF * ZMYM3D_WORK1(JI, JJ, JK) * ZGY_M_V3D_WORK1(JI, JJ, JK)
    END DO
  END DO
END DO

!

    ZFLX(:,:,IKE+1) = ZFLX(:,:,IKE)

  ELSE

    ZFLX(:,:,:)     = 0.

  END IF
!
! Compute the flux at the first inner U-point with an uncentred vertical  
! gradient
  CALL MYM2D_PHY(D, PK(:,:,IKB), ZMYM2D_WORK1)
CALL DYM2D_PHY(D, PRM(:,:,IKB,1), ZDYM2D_WORK1)

DO JJ=1, IJT
  DO JI=1, IIT
    ZSHUGRADWK1_2D(JI, JJ) = ZCOEFF(JI, JJ, IKB+2)*PRM(JI, JJ, IKB+2, 1)        &
               +ZCOEFF(JI, JJ, IKB+1)*PRM(JI, JJ, IKB+1, 1)       &
               +ZCOEFF(JI, JJ, IKB)*PRM(JI, JJ, IKB, 1)
  END DO
END DO

!
CALL MYM2D_PHY(D, ZSHUGRADWK1_2D, ZMYM2D_WORK2)

DO JJ=1, IJT
  DO JI=1, IIT
    ZFLX(JI, JJ, IKB) = -TURBN%XCHF * ZMYM2D_WORK1(JI, JJ) *           &
        ( ZDYM2D_WORK1(JI, JJ) * PINV_PDYY(JI, JJ, IKB)           &
         -ZMYM2D_WORK2(JI, JJ)     &
               *0.5* ( PDZY(JI, JJ, IKB+1)+PDZY(JI, JJ, IKB))      &
              * PINV_PDYY(JI, JJ, IKB) )
  END DO
END DO

!
! extrapolates the flux under the ground so that the vertical average with 
! the IKB flux gives the ground value  ( warning the tangential surface
! flux has been set to 0 for the moment !!  to be improved )
  
DO JJ=1, IJT
  DO JI=1, IIT
    ZSHUGRADWK1_2D(JI, JJ) = PSFRM(JI, JJ)* PDIRCOSYW(JI, JJ)
  END DO
END DO

!
CALL MYM2D_PHY(D, ZSHUGRADWK1_2D, ZMYM2D_WORK1)

DO JJ=1, IJT
  DO JI=1, IIT
    ZFLX(JI, JJ, IKB-1) = 2. * ZMYM2D_WORK1(JI, JJ) &
                           - ZFLX(JI, JJ, IKB)  
  END DO
END DO

!
!
  ! Add this source to the conservative mixing ratio sources
  !
  IF (.NOT. O2D) THEN 
    IF (.NOT. OFLAT) THEN
      CALL MYM_PHY(D, PRHODJ(:, :, :), ZMYM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMYM3D_WORK1(JI, JJ, JK) * ZFLX(JI, JJ, JK) * PINV_PDYY(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL DYF_PHY(D, ZSHUGRADWK1_3D, ZDYF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK3_3D(JI, JJ, JK) = ZFLX(JI, JJ, JK) * PINV_PDYY(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MZM_PHY(D, ZSHUGRADWK3_3D, ZMZM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = PDZY(JI, JJ, JK)*(ZMZM3D_WORK1(JI, JJ, JK))
    END DO
  END DO
END DO

!
CALL MYF_PHY(D, ZSHUGRADWK2_3D, ZMYF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = PMZM_PRHODJ(JI, JJ, JK) *ZMYF3D_WORK1(JI, JJ, JK) * PINV_PDZZ(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL DZF_PHY(D, ZSHUGRADWK1_3D, ZDZF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      PRRS(JI, JJ, JK, 1) = PRRS(JI, JJ, JK, 1)                                              &
                          - ZDYF3D_WORK1(JI, JJ, JK)                           &

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
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMYM3D_WORK1(JI, JJ, JK) * ZFLX(JI, JJ, JK) * PINV_PDYY(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL DYF_PHY(D, ZSHUGRADWK1_3D, ZDYF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      PRRS(JI, JJ, JK, 1) = PRRS(JI, JJ, JK, 1) - ZDYF3D_WORK1(JI, JJ, JK)    
    END DO
  END DO
END DO

!
END IF
  END IF
  !
  ! Compute the equivalent tendancy for Rc and Ri
  !
  IF ( KRRL >= 1 .AND. .NOT. O2D) THEN   ! Sub-grid condensation
    IF (.NOT. OFLAT) THEN
      
DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = PRHODJ(JI, JJ, JK)*PAMOIST(JI, JJ, JK)*PSRCM(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MYM_PHY(D, ZSHUGRADWK2_3D, ZMYM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMYM3D_WORK1(JI, JJ, JK)*ZFLX(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MXF_PHY(D, ZSHUGRADWK1_3D, ZMXF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = PRHODJ(JI, JJ, JK)*PAMOIST(JI, JJ, JK)*PSRCM(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MZM_PHY(D, ZSHUGRADWK2_3D, ZMZM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = ZFLX(JI, JJ, JK)*PINV_PDYY(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MZM_PHY(D, ZSHUGRADWK2_3D, ZMZM3D_WORK2)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = PDZY(JI, JJ, JK)*(ZMZM3D_WORK2(JI, JJ, JK))
    END DO
  END DO
END DO

!
CALL MYF_PHY(D, ZSHUGRADWK1_3D, ZMYF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMZM3D_WORK1(JI, JJ, JK)*ZMYF3D_WORK1(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MZF_PHY(D, ZSHUGRADWK1_3D, ZMZF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZFLXC(JI, JJ, JK) = ZFLXC(JI, JJ, JK)            &
                  + 2.*( ZMXF3D_WORK1(JI, JJ, JK)                     &
                      +  ZMZF3D_WORK1(JI, JJ, JK)&
                       )      
    END DO
  END DO
END DO

!
IF ( KRRI >= 1 ) THEN
        
DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = PRHODJ(JI, JJ, JK)*PAMOIST(JI, JJ, JK)*PSRCM(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MYM_PHY(D, ZSHUGRADWK2_3D, ZMYM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMYM3D_WORK1(JI, JJ, JK)*ZFLX(JI, JJ, JK)/PDYY(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL DYF_PHY(D, ZSHUGRADWK1_3D, ZDYF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = PRHODJ(JI, JJ, JK)*PAMOIST(JI, JJ, JK)*PSRCM(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MZM_PHY(D, ZSHUGRADWK2_3D, ZMZM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = ZFLX(JI, JJ, JK)*PINV_PDYY(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MZM_PHY(D, ZSHUGRADWK2_3D, ZMZM3D_WORK2)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = PDZY(JI, JJ, JK)*(ZMZM3D_WORK2(JI, JJ, JK))
    END DO
  END DO
END DO

!
CALL MYF_PHY(D, ZSHUGRADWK1_3D, ZMYF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMZM3D_WORK1(JI, JJ, JK)*ZMYF3D_WORK1(JI, JJ, JK)&
                                                 * PINV_PDZZ(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL DZF_PHY(D, ZSHUGRADWK1_3D, ZDZF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      PRRS(JI, JJ, JK, 2) = PRRS(JI, JJ, JK, 2) +  2. *                                  &
              (- ZDYF3D_WORK1(JI, JJ, JK)                        &
               + ZDZF3D_WORK1(JI, JJ, JK)                       &
              )*(1.0-PFRAC_ICE(JI, JJ, JK))        
    END DO
  END DO
END DO

!

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = PRHODJ(JI, JJ, JK)*PAMOIST(JI, JJ, JK)*PSRCM(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MYM_PHY(D, ZSHUGRADWK2_3D, ZMYM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMYM3D_WORK1(JI, JJ, JK)*ZFLX(JI, JJ, JK)/PDYY(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL DYF_PHY(D, ZSHUGRADWK1_3D, ZDYF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = PRHODJ(JI, JJ, JK)*PAMOIST(JI, JJ, JK)*PSRCM(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MZM_PHY(D, ZSHUGRADWK2_3D, ZMZM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = ZFLX(JI, JJ, JK)*PINV_PDYY(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MZM_PHY(D, ZSHUGRADWK2_3D, ZMZM3D_WORK2)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = PDZY(JI, JJ, JK)*(ZMZM3D_WORK2(JI, JJ, JK))
    END DO
  END DO
END DO

!
CALL MYF_PHY(D, ZSHUGRADWK1_3D, ZMYF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMZM3D_WORK1(JI, JJ, JK)*ZMYF3D_WORK1(JI, JJ, JK)&
                                                 * PINV_PDZZ(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL DZF_PHY(D, ZSHUGRADWK1_3D, ZDZF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      PRRS(JI, JJ, JK, 4) = PRRS(JI, JJ, JK, 4) +  2. *                                  &
              (- ZDYF3D_WORK1(JI, JJ, JK)                        &
               + ZDZF3D_WORK1(JI, JJ, JK)                       &
              )*PFRAC_ICE(JI, JJ, JK)      
    END DO
  END DO
END DO

!
ELSE
        
DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = PRHODJ(JI, JJ, JK)*PAMOIST(JI, JJ, JK)*PSRCM(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MYM_PHY(D, ZSHUGRADWK2_3D, ZMYM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMYM3D_WORK1(JI, JJ, JK)*ZFLX(JI, JJ, JK)/PDYY(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL DYF_PHY(D, ZSHUGRADWK1_3D, ZDYF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = PRHODJ(JI, JJ, JK)*PAMOIST(JI, JJ, JK)*PSRCM(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MZM_PHY(D, ZSHUGRADWK2_3D, ZMZM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = ZFLX(JI, JJ, JK)*PINV_PDYY(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MZM_PHY(D, ZSHUGRADWK2_3D, ZMZM3D_WORK2)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = PDZY(JI, JJ, JK)*(ZMZM3D_WORK2(JI, JJ, JK))
    END DO
  END DO
END DO

!
CALL MYF_PHY(D, ZSHUGRADWK1_3D, ZMYF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMZM3D_WORK1(JI, JJ, JK)*ZMYF3D_WORK1(JI, JJ, JK)&
                                                 * PINV_PDZZ(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL DZF_PHY(D, ZSHUGRADWK1_3D, ZDZF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      PRRS(JI, JJ, JK, 2) = PRRS(JI, JJ, JK, 2) +  2. *                                  &
              (- ZDYF3D_WORK1(JI, JJ, JK)                        &
               + ZDZF3D_WORK1(JI, JJ, JK)                       &
              )      
    END DO
  END DO
END DO

!
END IF
    ELSE
      
DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = PRHODJ(JI, JJ, JK)*PAMOIST(JI, JJ, JK)*PSRCM(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MYM_PHY(D, ZSHUGRADWK2_3D, ZMYM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMYM3D_WORK1(JI, JJ, JK)*ZFLX(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MXF_PHY(D, ZSHUGRADWK1_3D, ZMXF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZFLXC(JI, JJ, JK) = ZFLXC(JI, JJ, JK) + 2.*ZMXF3D_WORK1(JI, JJ, JK)      
    END DO
  END DO
END DO

!
IF ( KRRI >= 1 ) THEN
        
DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = PRHODJ(JI, JJ, JK)*PAMOIST(JI, JJ, JK)*PSRCM(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MYM_PHY(D, ZSHUGRADWK2_3D, ZMYM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMYM3D_WORK1(JI, JJ, JK)*ZFLX(JI, JJ, JK)/PDYY(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL DYF_PHY(D, ZSHUGRADWK1_3D, ZDYF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      PRRS(JI, JJ, JK, 2) = PRRS(JI, JJ, JK, 2) - 2. *                                   &
              ZDYF3D_WORK1(JI, JJ, JK)*(1.0-PFRAC_ICE(JI, JJ, JK))        
    END DO
  END DO
END DO

!

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = PRHODJ(JI, JJ, JK)*PAMOIST(JI, JJ, JK)*PSRCM(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MYM_PHY(D, ZSHUGRADWK2_3D, ZMYM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMYM3D_WORK1(JI, JJ, JK)*ZFLX(JI, JJ, JK)/PDYY(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL DYF_PHY(D, ZSHUGRADWK1_3D, ZDYF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      PRRS(JI, JJ, JK, 4) = PRRS(JI, JJ, JK, 4) - 2. *                                   &
              ZDYF3D_WORK1(JI, JJ, JK)*PFRAC_ICE(JI, JJ, JK)      
    END DO
  END DO
END DO

!
ELSE
        
DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = PRHODJ(JI, JJ, JK)*PAMOIST(JI, JJ, JK)*PSRCM(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MYM_PHY(D, ZSHUGRADWK2_3D, ZMYM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMYM3D_WORK1(JI, JJ, JK)*ZFLX(JI, JJ, JK)/PDYY(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL DYF_PHY(D, ZSHUGRADWK1_3D, ZDYF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      PRRS(JI, JJ, JK, 2) = PRRS(JI, JJ, JK, 2) - 2. *                                   &
              ZDYF3D_WORK1(JI, JJ, JK)      
    END DO
  END DO
END DO

!
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

    CALL IO_FIELD_WRITE(TPFILE,TZFIELD,ZFLX)
  END IF
  !
  IF (KSPLT==1 .AND. TLES%LLES_CALL) THEN
    CALL SECOND_MNH(ZTIME1)
    CALL MYF_PHY(D, ZFLX(:, :, :), ZMYF3D_WORK1)
CALL LES_MEAN_SUBGRID(ZMYF3D_WORK1(:, :, :), TLES%X_LES_SUBGRID_VRt )     
CALL GY_W_VW_PHY(D, OFLAT, PWM(:, :, :),PDYY(:, :, :),PDZZ(:, :, :),PDZY(:, :, :), ZGY_W_VW3D_WORK1)
CALL MZM_PHY(D, ZFLX(:, :, :), ZMZM3D_WORK1)

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
CALL LES_MEAN_SUBGRID(ZMZF3D_WORK1(:, :, :),&
                           TLES%X_LES_RES_ddxa_W_SBG_UaRt , .TRUE. )    
CALL GY_M_M_PHY(D, OFLAT,PTHLM(:, :, :),PDYY(:, :, :),PDZZ(:, :, :),PDZY(:, :, :), ZGY_M_M3D_WORK1)
CALL MYF_PHY(D, ZFLX(:, :, :), ZMYF3D_WORK1)
CALL LES_MEAN_SUBGRID(ZGY_M_M3D_WORK1(:, :, :)*ZMYF3D_WORK1(:, :, :), &
                           TLES%X_LES_RES_ddxa_Thl_SBG_UaRt , .TRUE. )    
CALL GY_M_M_PHY(D, OFLAT,PRM(:,:,:,1),PDYY(:, :, :),PDZZ(:, :, :),PDZY(:, :, :), ZGY_M_M3D_WORK1)
CALL MYF_PHY(D, ZFLX(:, :, :), ZMYF3D_WORK1)
CALL LES_MEAN_SUBGRID(ZGY_M_M3D_WORK1(:, :, :)*ZMYF3D_WORK1(:, :, :), &
                           TLES%X_LES_RES_ddxa_Rt_SBG_UaRt , .TRUE. )    
CALL SECOND_MNH(ZTIME2)
    TLES%XTIME_LES = TLES%XTIME_LES + ZTIME2 - ZTIME1
  END IF
!
  !
  IF (KRRL>0 .AND. KSPLT==1 .AND. TLES%LLES_CALL) THEN
    CALL SECOND_MNH(ZTIME1)
    CALL MYF_PHY(D, ZFLXC(:, :, :), ZMYF3D_WORK1)
CALL LES_MEAN_SUBGRID(ZMYF3D_WORK1(:, :, :), TLES%X_LES_SUBGRID_VRc )    
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
