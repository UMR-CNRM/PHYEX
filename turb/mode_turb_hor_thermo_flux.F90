!MNH_LIC Copyright 1994-2025 CNRS, Meteo-France and Universite Paul Sabatier
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
USE MODE_SHUMAN_PHY, ONLY: DXF_PHY, DXM2D_PHY, DYF_PHY, DYM2D_PHY, DZF_PHY, MXF_PHY, MXM2D_PHY, MXM_PHY, MYF_PHY, MYM2D_PHY, MYM_PHY, MZF_PHY, MZM_PHY
USE MODE_GRADIENT_M_PHY, ONLY: GX_M_M_PHY, GX_M_U_PHY, GY_M_M_PHY, GY_M_V_PHY
USE MODE_GRADIENT_W_PHY, ONLY: GX_W_UW_PHY, GY_W_VW_PHY
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
TYPE(TFILEDATA),          INTENT(INOUT) ::  TPFILE       ! Output file
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    ::  PK          ! Turbulent diffusion doef.
                                                        ! PK = PLM * SQRT(PTKEM)
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    ::  PINV_PDXX   ! 1./PDXX
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    ::  PINV_PDYY   ! 1./PDYY
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    ::  PINV_PDZZ   ! 1./PDZZ
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    ::  PMZM_PRHODJ ! MZM(PRHODJ)
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    ::  PDXX, PDYY, PDZZ, PDZX, PDZY 
                                                         ! Metric coefficients
REAL, DIMENSION(D%NIJT),     INTENT(IN)    ::  PDIRCOSXW, PDIRCOSYW
! Director Cosinus along x, y and z directions at surface w-point
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    ::  PRHODJ       ! density * grid volume
!
REAL, DIMENSION(D%NIJT),     INTENT(IN)    ::  PSFTHM,PSFRM
!
! Variables at t-1
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    ::  PWM 
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    ::  PTHLM 
REAL, DIMENSION(D%NIJT,D%NKT,KRR), INTENT(IN)    ::  PRM          ! mixing ratios at t-1,
                              !  where PRM(:,:,1) = conservative mixing ratio
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    ::  PATHETA      ! coefficients between 
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    ::  PAMOIST      ! s and Thetal and Rnp

REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    ::  PSRCM
                                  ! normalized 2nd-order flux
                                  ! s'r'c/2Sigma_s2 at t-1 multiplied by Lambda_3
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    ::  PFRAC_ICE    ! ri fraction of rc+ri
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) ::  PRTHLS
REAL, DIMENSION(D%NIJT,D%NKT,KRR), INTENT(INOUT) ::  PRRS         ! var. at t+1 -split-
!
!
!
!*       0.2  declaration of local variables
!
REAL, DIMENSION(D%NIJT,D%NKT) :: ZFLX,ZFLXC,ZWKLES ! work arrays
!
INTEGER             :: IKB,IKE,IKU, IKT, IIJT
                                    ! Index values for the Beginning and End
                                    ! mass points of the domain  
REAL, DIMENSION(D%NIJT,1+JPVEXT:3+JPVEXT) :: ZCOEFF 
                                    ! coefficients for the uncentred gradient 
                                    ! computation near the ground
!
REAL :: ZTIME1, ZTIME2
INTEGER :: JIJ
INTEGER :: JK
REAL, DIMENSION(D%NIJT) :: ZWORKSHUGRAD_2D_1
REAL, DIMENSION(D%NIJT) :: ZWORKSHUGRAD_2D_2
REAL, DIMENSION(D%NIJT) :: ZWORKSHUGRAD_2D_3
REAL, DIMENSION(D%NIJT) :: ZWORKSHUGRAD_2D_4
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
TYPE(TFIELDMETADATA) :: TZFIELD
! ---------------------------------------------------------------------------
!
!*       1.   PRELIMINARY COMPUTATIONS
!             ------------------------
!
IKB = 1+JPVEXT               
IKE = SIZE(PTHLM,2)-JPVEXT    
IKU = SIZE(PTHLM,2)
IIJT=D%NIJT
IKT=D%NKT
!
!
!  compute the coefficients for the uncentred gradient computation near the 
!  ground

ZCOEFF(:,IKB+2)= - PDZZ(:,IKB+1) /      &
       ( (PDZZ(:,IKB+2)+PDZZ(:,IKB+1)) * PDZZ(:,IKB+2) )
ZCOEFF(:,IKB+1)=   (PDZZ(:,IKB+2)+PDZZ(:,IKB+1)) /      &
       ( PDZZ(:,IKB+1) * PDZZ(:,IKB+2) )
ZCOEFF(:,IKB)= - (PDZZ(:,IKB+2)+2.*PDZZ(:,IKB+1)) /      &
       ( (PDZZ(:,IKB+2)+PDZZ(:,IKB+1)) * PDZZ(:,IKB+1) )

!
!*       2.   < U' THETA'l >
!             --------------
!
!
#define UTHETA UTHETA_ON
!!!!!!! UTHETA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef UTHETA
CALL MXM_PHY(D, PK, ZWORKSHUGRAD_3D_1)
CALL GX_M_U_PHY(D, OFLAT, PTHLM, PDXX, PDZZ, PDZX, ZWORKSHUGRAD_3D_2)

DO JK=1, D%NKT
  DO JIJ=1, D%NIJT
    ZFLX(JIJ, JK)     = -TURBN%XCSHF * ZWORKSHUGRAD_3D_1(JIJ, JK) * ZWORKSHUGRAD_3D_2(JIJ, JK)
  END DO
END DO


ZFLX(:,IKE+1) = ZFLX(:,IKE) 

!
! Compute the flux at the first inner U-point with an uncentred vertical  
! gradient
CALL DXM2D_PHY(D, PTHLM(:,IKB), ZWORKSHUGRAD_2D_1)

DO JIJ=1, D%NIJT
  ZWORKSHUGRAD_2D_2(JIJ)=ZCOEFF(JIJ, IKB+2)*PTHLM(JIJ, IKB+2)        &
           +ZCOEFF(JIJ, IKB+1)*PTHLM(JIJ, IKB+1)       &
           +ZCOEFF(JIJ, IKB)*PTHLM(JIJ, IKB)
END DO

CALL MXM2D_PHY(D, ZWORKSHUGRAD_2D_2, ZWORKSHUGRAD_2D_3)
CALL MXM2D_PHY(D, PK(:,IKB), ZWORKSHUGRAD_2D_4)

DO JIJ=1, D%NIJT
  ZFLX(JIJ, IKB) = -TURBN%XCSHF * ZWORKSHUGRAD_2D_4(JIJ) *          &
    ( ZWORKSHUGRAD_2D_1(JIJ) * PINV_PDXX(JIJ, IKB)           &
     -ZWORKSHUGRAD_2D_3(JIJ)      &
          *0.5* ( PDZX(JIJ, IKB+1)+PDZX(JIJ, IKB))       &
          * PINV_PDXX(JIJ, IKB) )
END DO

! extrapolates the flux under the ground so that the vertical average with 
! the IKB flux gives the ground value  ( warning the tangential surface
! flux has been set to 0 for the moment !!  to be improved )

DO JIJ=1, D%NIJT
  ZWORKSHUGRAD_2D_1(JIJ)=PSFTHM(JIJ)* PDIRCOSXW(JIJ)
END DO

CALL MXM2D_PHY(D, ZWORKSHUGRAD_2D_1, ZWORKSHUGRAD_2D_2)

DO JIJ=1, D%NIJT
  ZFLX(JIJ, IKB-1) = 2. * ZWORKSHUGRAD_2D_2(JIJ)  &
                         - ZFLX(JIJ, IKB)
END DO

!
! Add this source to the Theta_l sources
!
IF (.NOT. OFLAT) THEN
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      ZWORKSHUGRAD_3D_1(JIJ, JK)=ZFLX(JIJ, JK) * PINV_PDXX(JIJ, JK)
    END DO
  END DO
  
  CALL MZM_PHY(D, ZWORKSHUGRAD_3D_1, ZWORKSHUGRAD_3D_2)
  CALL MXM_PHY(D, PRHODJ, ZWORKSHUGRAD_3D_3)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      ZWORKSHUGRAD_3D_4(JIJ, JK)=PDZX(JIJ, JK)*(ZWORKSHUGRAD_3D_2(JIJ, JK))
    END DO
  END DO
  
  CALL MXF_PHY(D, ZWORKSHUGRAD_3D_4, ZWORKSHUGRAD_3D_5)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      ZWORKSHUGRAD_3D_6(JIJ, JK)=ZWORKSHUGRAD_3D_3(JIJ, JK) * ZFLX(JIJ, JK) * PINV_PDXX(JIJ, JK)
    END DO
  END DO
  
  CALL DXF_PHY(D, ZWORKSHUGRAD_3D_6, ZWORKSHUGRAD_3D_7)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      ZWORKSHUGRAD_3D_8(JIJ, JK)=PMZM_PRHODJ(JIJ, JK) *ZWORKSHUGRAD_3D_5(JIJ, JK) * PINV_PDZZ(JIJ, JK)
    END DO
  END DO
  
  CALL DZF_PHY(D, ZWORKSHUGRAD_3D_8, ZWORKSHUGRAD_3D_9)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      PRTHLS(JIJ, JK) =  PRTHLS(JIJ, JK)                                                   &
                    - ZWORKSHUGRAD_3D_7(JIJ, JK)                          &
                    + ZWORKSHUGRAD_3D_9(JIJ, JK)
    END DO
  END DO
  
ELSE
  CALL MXM_PHY(D, PRHODJ, ZWORKSHUGRAD_3D_1)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      ZWORKSHUGRAD_3D_2(JIJ, JK)=ZWORKSHUGRAD_3D_1(JIJ, JK) * ZFLX(JIJ, JK) * PINV_PDXX(JIJ, JK)
    END DO
  END DO
  
  CALL DXF_PHY(D, ZWORKSHUGRAD_3D_2, ZWORKSHUGRAD_3D_3)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      PRTHLS(JIJ, JK) =  PRTHLS(JIJ, JK) - ZWORKSHUGRAD_3D_3(JIJ, JK)
    END DO
  END DO
  
END IF
!
! Compute the equivalent tendancy for Rc and Ri
!
IF ( KRRL >= 1 ) THEN
  IF (.NOT. OFLAT) THEN
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWORKSHUGRAD_3D_1(JIJ, JK)=ZFLX(JIJ, JK)*PINV_PDXX(JIJ, JK)
      END DO
    END DO
    
    CALL MZM_PHY(D, ZWORKSHUGRAD_3D_1, ZWORKSHUGRAD_3D_2)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWORKSHUGRAD_3D_3(JIJ, JK)=PRHODJ(JIJ, JK)*PATHETA(JIJ, JK)*PSRCM(JIJ, JK)
      END DO
    END DO
    
    CALL MXM_PHY(D, ZWORKSHUGRAD_3D_3, ZWORKSHUGRAD_3D_4)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWORKSHUGRAD_3D_5(JIJ, JK)=PRHODJ(JIJ, JK)*PATHETA(JIJ, JK)*PSRCM(JIJ, JK)
      END DO
    END DO
    
    CALL MZM_PHY(D, ZWORKSHUGRAD_3D_5, ZWORKSHUGRAD_3D_6)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWORKSHUGRAD_3D_7(JIJ, JK)=PDZX(JIJ, JK)*(ZWORKSHUGRAD_3D_2(JIJ, JK))
      END DO
    END DO
    
    CALL MXF_PHY(D, ZWORKSHUGRAD_3D_7, ZWORKSHUGRAD_3D_8)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWORKSHUGRAD_3D_9(JIJ, JK)=ZWORKSHUGRAD_3D_4(JIJ, JK)*ZFLX(JIJ, JK)
      END DO
    END DO
    
    CALL MXF_PHY(D, ZWORKSHUGRAD_3D_9, ZWORKSHUGRAD_3D_10)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWORKSHUGRAD_3D_11(JIJ, JK)=ZWORKSHUGRAD_3D_6(JIJ, JK)*ZWORKSHUGRAD_3D_8(JIJ, JK)
      END DO
    END DO
    
    CALL MZF_PHY(D, ZWORKSHUGRAD_3D_11, ZWORKSHUGRAD_3D_12)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZFLXC(JIJ, JK) = 2.*( ZWORKSHUGRAD_3D_10(JIJ, JK)                       &
                    +ZWORKSHUGRAD_3D_12(JIJ, JK)&
                   )
      END DO
    END DO
    
    IF ( KRRI >= 1 ) THEN
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_1(JIJ, JK)=ZFLX(JIJ, JK)*PINV_PDXX(JIJ, JK)
        END DO
      END DO
      
      CALL MZM_PHY(D, ZWORKSHUGRAD_3D_1, ZWORKSHUGRAD_3D_2)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_3(JIJ, JK)=PRHODJ(JIJ, JK)*PATHETA(JIJ, JK)*PSRCM(JIJ, JK)
        END DO
      END DO
      
      CALL MXM_PHY(D, ZWORKSHUGRAD_3D_3, ZWORKSHUGRAD_3D_4)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_5(JIJ, JK)=PDZX(JIJ, JK)*(ZWORKSHUGRAD_3D_2(JIJ, JK))
        END DO
      END DO
      
      CALL MXF_PHY(D, ZWORKSHUGRAD_3D_5, ZWORKSHUGRAD_3D_6)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_7(JIJ, JK)=PRHODJ(JIJ, JK)*PATHETA(JIJ, JK)*PSRCM(JIJ, JK)
        END DO
      END DO
      
      CALL MZM_PHY(D, ZWORKSHUGRAD_3D_7, ZWORKSHUGRAD_3D_8)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_9(JIJ, JK)=ZWORKSHUGRAD_3D_4(JIJ, JK)*ZFLX(JIJ, JK)*PINV_PDXX(JIJ, JK)
        END DO
      END DO
      
      CALL DXF_PHY(D, ZWORKSHUGRAD_3D_9, ZWORKSHUGRAD_3D_10)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_11(JIJ, JK)=ZWORKSHUGRAD_3D_8(JIJ, JK)*ZWORKSHUGRAD_3D_6(JIJ, JK)&
                                               *PINV_PDZZ(JIJ, JK)
        END DO
      END DO
      
      CALL DZF_PHY(D, ZWORKSHUGRAD_3D_11, ZWORKSHUGRAD_3D_12)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          PRRS(JIJ, JK, 2) = PRRS(JIJ, JK, 2) +  2. *                                    &
            (- ZWORKSHUGRAD_3D_10(JIJ, JK)                   &
             + ZWORKSHUGRAD_3D_12(JIJ, JK)                        &
            )*(1.0-PFRAC_ICE(JIJ, JK))
        END DO
      END DO
      
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_1(JIJ, JK)=ZFLX(JIJ, JK)*PINV_PDXX(JIJ, JK)
        END DO
      END DO
      
      CALL MZM_PHY(D, ZWORKSHUGRAD_3D_1, ZWORKSHUGRAD_3D_2)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_3(JIJ, JK)=PRHODJ(JIJ, JK)*PATHETA(JIJ, JK)*PSRCM(JIJ, JK)
        END DO
      END DO
      
      CALL MXM_PHY(D, ZWORKSHUGRAD_3D_3, ZWORKSHUGRAD_3D_4)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_5(JIJ, JK)=PDZX(JIJ, JK)*(ZWORKSHUGRAD_3D_2(JIJ, JK))
        END DO
      END DO
      
      CALL MXF_PHY(D, ZWORKSHUGRAD_3D_5, ZWORKSHUGRAD_3D_6)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_7(JIJ, JK)=PRHODJ(JIJ, JK)*PATHETA(JIJ, JK)*PSRCM(JIJ, JK)
        END DO
      END DO
      
      CALL MZM_PHY(D, ZWORKSHUGRAD_3D_7, ZWORKSHUGRAD_3D_8)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_9(JIJ, JK)=ZWORKSHUGRAD_3D_4(JIJ, JK)*ZFLX(JIJ, JK)*PINV_PDXX(JIJ, JK)
        END DO
      END DO
      
      CALL DXF_PHY(D, ZWORKSHUGRAD_3D_9, ZWORKSHUGRAD_3D_10)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_11(JIJ, JK)=ZWORKSHUGRAD_3D_8(JIJ, JK)*ZWORKSHUGRAD_3D_6(JIJ, JK)&
                                               *PINV_PDZZ(JIJ, JK)
        END DO
      END DO
      
      CALL DZF_PHY(D, ZWORKSHUGRAD_3D_11, ZWORKSHUGRAD_3D_12)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          PRRS(JIJ, JK, 4) = PRRS(JIJ, JK, 4) +  2. *                                    &
            (- ZWORKSHUGRAD_3D_10(JIJ, JK)                   &
             + ZWORKSHUGRAD_3D_12(JIJ, JK)                        &
            )*PFRAC_ICE(JIJ, JK)
        END DO
      END DO
      
    ELSE
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_1(JIJ, JK)=ZFLX(JIJ, JK)*PINV_PDXX(JIJ, JK)
        END DO
      END DO
      
      CALL MZM_PHY(D, ZWORKSHUGRAD_3D_1, ZWORKSHUGRAD_3D_2)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_3(JIJ, JK)=PRHODJ(JIJ, JK)*PATHETA(JIJ, JK)*PSRCM(JIJ, JK)
        END DO
      END DO
      
      CALL MXM_PHY(D, ZWORKSHUGRAD_3D_3, ZWORKSHUGRAD_3D_4)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_5(JIJ, JK)=PDZX(JIJ, JK)*(ZWORKSHUGRAD_3D_2(JIJ, JK))
        END DO
      END DO
      
      CALL MXF_PHY(D, ZWORKSHUGRAD_3D_5, ZWORKSHUGRAD_3D_6)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_7(JIJ, JK)=PRHODJ(JIJ, JK)*PATHETA(JIJ, JK)*PSRCM(JIJ, JK)
        END DO
      END DO
      
      CALL MZM_PHY(D, ZWORKSHUGRAD_3D_7, ZWORKSHUGRAD_3D_8)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_9(JIJ, JK)=ZWORKSHUGRAD_3D_4(JIJ, JK)*ZFLX(JIJ, JK)*PINV_PDXX(JIJ, JK)
        END DO
      END DO
      
      CALL DXF_PHY(D, ZWORKSHUGRAD_3D_9, ZWORKSHUGRAD_3D_10)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_11(JIJ, JK)=ZWORKSHUGRAD_3D_8(JIJ, JK)*ZWORKSHUGRAD_3D_6(JIJ, JK)&
                                               *PINV_PDZZ(JIJ, JK)
        END DO
      END DO
      
      CALL DZF_PHY(D, ZWORKSHUGRAD_3D_11, ZWORKSHUGRAD_3D_12)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          PRRS(JIJ, JK, 2) = PRRS(JIJ, JK, 2) +  2. *                                    &
            (- ZWORKSHUGRAD_3D_10(JIJ, JK)                   &
             + ZWORKSHUGRAD_3D_12(JIJ, JK)                        &
            )
        END DO
      END DO
      
    END IF
  ELSE
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWORKSHUGRAD_3D_1(JIJ, JK)=PRHODJ(JIJ, JK)*PATHETA(JIJ, JK)*PSRCM(JIJ, JK)
      END DO
    END DO
    
    CALL MXM_PHY(D, ZWORKSHUGRAD_3D_1, ZWORKSHUGRAD_3D_2)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWORKSHUGRAD_3D_3(JIJ, JK)=ZWORKSHUGRAD_3D_2(JIJ, JK)*ZFLX(JIJ, JK)
      END DO
    END DO
    
    CALL MXF_PHY(D, ZWORKSHUGRAD_3D_3, ZWORKSHUGRAD_3D_4)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZFLXC(JIJ, JK) = 2.*ZWORKSHUGRAD_3D_4(JIJ, JK)
      END DO
    END DO
    
    IF ( KRRI >= 1 ) THEN
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_1(JIJ, JK)=PRHODJ(JIJ, JK)*PATHETA(JIJ, JK)*PSRCM(JIJ, JK)
        END DO
      END DO
      
      CALL MXM_PHY(D, ZWORKSHUGRAD_3D_1, ZWORKSHUGRAD_3D_2)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_3(JIJ, JK)=ZWORKSHUGRAD_3D_2(JIJ, JK)*ZFLX(JIJ, JK)*PINV_PDXX(JIJ, JK)
        END DO
      END DO
      
      CALL DXF_PHY(D, ZWORKSHUGRAD_3D_3, ZWORKSHUGRAD_3D_4)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          PRRS(JIJ, JK, 2) = PRRS(JIJ, JK, 2) -  2. *                                    &
            ZWORKSHUGRAD_3D_4(JIJ, JK)*(1.0-PFRAC_ICE(JIJ, JK))
        END DO
      END DO
      
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_1(JIJ, JK)=PRHODJ(JIJ, JK)*PATHETA(JIJ, JK)*PSRCM(JIJ, JK)
        END DO
      END DO
      
      CALL MXM_PHY(D, ZWORKSHUGRAD_3D_1, ZWORKSHUGRAD_3D_2)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_3(JIJ, JK)=ZWORKSHUGRAD_3D_2(JIJ, JK)*ZFLX(JIJ, JK)*PINV_PDXX(JIJ, JK)
        END DO
      END DO
      
      CALL DXF_PHY(D, ZWORKSHUGRAD_3D_3, ZWORKSHUGRAD_3D_4)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          PRRS(JIJ, JK, 4) = PRRS(JIJ, JK, 4) -  2. *                                    &
            ZWORKSHUGRAD_3D_4(JIJ, JK)*PFRAC_ICE(JIJ, JK)
        END DO
      END DO
      
    ELSE
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_1(JIJ, JK)=PRHODJ(JIJ, JK)*PATHETA(JIJ, JK)*PSRCM(JIJ, JK)
        END DO
      END DO
      
      CALL MXM_PHY(D, ZWORKSHUGRAD_3D_1, ZWORKSHUGRAD_3D_2)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_3(JIJ, JK)=ZWORKSHUGRAD_3D_2(JIJ, JK)*ZFLX(JIJ, JK)*PINV_PDXX(JIJ, JK)
        END DO
      END DO
      
      CALL DXF_PHY(D, ZWORKSHUGRAD_3D_3, ZWORKSHUGRAD_3D_4)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          PRRS(JIJ, JK, 2) = PRRS(JIJ, JK, 2) -  2. *                                    &
            ZWORKSHUGRAD_3D_4(JIJ, JK)
        END DO
      END DO
      
    END IF
  END IF
END IF
!
!! stores this flux in ZWORK to compute later <U' VPT'>
!!ZWORK(:,:) = ZFLX(:,:) 
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

  CALL IO_FIELD_WRITE_PHY(D,TPFILE,TZFIELD,ZFLX)
END IF
!
IF (KSPLT==1 .AND. TLES%LLES_CALL) THEN
  CALL SECOND_MNH(ZTIME1)
  CALL MXF_PHY(D, ZFLX, ZWORKSHUGRAD_3D_1)
  CALL LES_MEAN_SUBGRID_PHY(D, TLES,ZWORKSHUGRAD_3D_1(1:D%NIJT, 1:D%NKT), TLES%X_LES_SUBGRID_UThl ) 
  CALL GX_W_UW_PHY(D, OFLAT, PWM, PDXX, PDZZ, PDZX, ZWORKSHUGRAD_3D_1)
  CALL MZM_PHY(D, ZFLX, ZWORKSHUGRAD_3D_2)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      ZWORKSHUGRAD_3D_3(JIJ, JK)=ZWORKSHUGRAD_3D_1(JIJ, JK)*ZWORKSHUGRAD_3D_2(JIJ, JK)
    END DO
  END DO
  
  CALL MXF_PHY(D, ZWORKSHUGRAD_3D_3, ZWORKSHUGRAD_3D_4)
  CALL MZF_PHY(D, ZWORKSHUGRAD_3D_4, ZWKLES)
  CALL LES_MEAN_SUBGRID_PHY(D, TLES,ZWKLES, TLES%X_LES_RES_ddxa_W_SBG_UaThl , .TRUE. )
  CALL GX_M_M_PHY(D, OFLAT, PTHLM, PDXX, PDZZ, PDZX, ZWORKSHUGRAD_3D_1)
  CALL MXF_PHY(D, ZFLX, ZWORKSHUGRAD_3D_2)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      ZWKLES(JIJ, JK) = ZWORKSHUGRAD_3D_1(JIJ, JK)*ZWORKSHUGRAD_3D_2(JIJ, JK)
    END DO
  END DO
  
  CALL LES_MEAN_SUBGRID_PHY(D, TLES,ZWKLES, TLES%X_LES_RES_ddxa_Thl_SBG_UaThl , .TRUE. )
  IF (KRR>=1) THEN
    CALL GX_M_M_PHY(D, OFLAT, PRM(:,:,1), PDXX, PDZZ, PDZX, ZWORKSHUGRAD_3D_1)
    CALL MXF_PHY(D, ZFLX, ZWORKSHUGRAD_3D_2)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWKLES(JIJ, JK) = ZWORKSHUGRAD_3D_1(JIJ, JK)*ZWORKSHUGRAD_3D_2(JIJ, JK)
      END DO
    END DO
    
    CALL LES_MEAN_SUBGRID_PHY(D, TLES,ZWKLES, TLES%X_LES_RES_ddxa_Rt_SBG_UaThl , .TRUE. )
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
  CALL MXM_PHY(D, PK, ZWORKSHUGRAD_3D_1)
  CALL GX_M_U_PHY(D, OFLAT, PRM(:,:,1), PDXX, PDZZ, PDZX, ZWORKSHUGRAD_3D_2)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      ZFLX(JIJ, JK)     = -TURBN%XCHF * ZWORKSHUGRAD_3D_1(JIJ, JK) * ZWORKSHUGRAD_3D_2(JIJ, JK)
    END DO
  END DO
  

  ZFLX(:,IKE+1) = ZFLX(:,IKE)

!
! Compute the flux at the first inner U-point with an uncentred vertical  
! gradient
  CALL DXM2D_PHY(D, PRM(:,IKB,1), ZWORKSHUGRAD_2D_1)
  
  DO JIJ=1, D%NIJT
    ZWORKSHUGRAD_2D_2(JIJ)=ZCOEFF(JIJ, IKB+2)*PRM(JIJ, IKB+2, 1)        &
             +ZCOEFF(JIJ, IKB+1)*PRM(JIJ, IKB+1, 1)       &
             +ZCOEFF(JIJ, IKB)*PRM(JIJ, IKB, 1)
  END DO
  
  CALL MXM2D_PHY(D, ZWORKSHUGRAD_2D_2, ZWORKSHUGRAD_2D_3)
  CALL MXM2D_PHY(D, PK(:,IKB), ZWORKSHUGRAD_2D_4)
  
  DO JIJ=1, D%NIJT
    ZFLX(JIJ, IKB) = -TURBN%XCHF * ZWORKSHUGRAD_2D_4(JIJ) *           &
      ( ZWORKSHUGRAD_2D_1(JIJ) * PINV_PDXX(JIJ, IKB)           &
       -ZWORKSHUGRAD_2D_3(JIJ)      &
            *0.5* ( PDZX(JIJ, IKB+1)+PDZX(JIJ, IKB))       &
            * PINV_PDXX(JIJ, IKB) )
  END DO
  
! extrapolates the flux under the ground so that the vertical average with 
! the IKB flux gives the ground value  ( warning the tangential surface
! flux has been set to 0 for the moment !!  to be improved )
  
  DO JIJ=1, D%NIJT
    ZWORKSHUGRAD_2D_1(JIJ)=PSFRM(JIJ)* PDIRCOSXW(JIJ)
  END DO
  
  CALL MXM2D_PHY(D, ZWORKSHUGRAD_2D_1, ZWORKSHUGRAD_2D_2)
  
  DO JIJ=1, D%NIJT
    ZFLX(JIJ, IKB-1) = 2. * ZWORKSHUGRAD_2D_2(JIJ) &
                         - ZFLX(JIJ, IKB)
  END DO
  
  !
  ! Add this source to the conservative mixing ratio sources
  !
  IF (.NOT. OFLAT) THEN
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWORKSHUGRAD_3D_1(JIJ, JK)=ZFLX(JIJ, JK) * PINV_PDXX(JIJ, JK)
      END DO
    END DO
    
    CALL MZM_PHY(D, ZWORKSHUGRAD_3D_1, ZWORKSHUGRAD_3D_2)
    CALL MXM_PHY(D, PRHODJ, ZWORKSHUGRAD_3D_3)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWORKSHUGRAD_3D_4(JIJ, JK)=PDZX(JIJ, JK)*(ZWORKSHUGRAD_3D_2(JIJ, JK))
      END DO
    END DO
    
    CALL MXF_PHY(D, ZWORKSHUGRAD_3D_4, ZWORKSHUGRAD_3D_5)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWORKSHUGRAD_3D_6(JIJ, JK)=ZWORKSHUGRAD_3D_3(JIJ, JK) * ZFLX(JIJ, JK) * PINV_PDXX(JIJ, JK)
      END DO
    END DO
    
    CALL DXF_PHY(D, ZWORKSHUGRAD_3D_6, ZWORKSHUGRAD_3D_7)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWORKSHUGRAD_3D_8(JIJ, JK)=PMZM_PRHODJ(JIJ, JK) *ZWORKSHUGRAD_3D_5(JIJ, JK) * PINV_PDZZ(JIJ, JK)
      END DO
    END DO
    
    CALL DZF_PHY(D, ZWORKSHUGRAD_3D_8, ZWORKSHUGRAD_3D_9)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        PRRS(JIJ, JK, 1) = PRRS(JIJ, JK, 1)                                             &
                      - ZWORKSHUGRAD_3D_7(JIJ, JK)                          &
                      + ZWORKSHUGRAD_3D_9(JIJ, JK)
      END DO
    END DO
    
  ELSE
    CALL MXM_PHY(D, PRHODJ, ZWORKSHUGRAD_3D_1)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWORKSHUGRAD_3D_2(JIJ, JK)=ZWORKSHUGRAD_3D_1(JIJ, JK) * ZFLX(JIJ, JK) * PINV_PDXX(JIJ, JK)
      END DO
    END DO
    
    CALL DXF_PHY(D, ZWORKSHUGRAD_3D_2, ZWORKSHUGRAD_3D_3)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        PRRS(JIJ, JK, 1) = PRRS(JIJ, JK, 1) - ZWORKSHUGRAD_3D_3(JIJ, JK)
      END DO
    END DO
    
  END IF
  !
  ! Compute the equivalent tendancy for Rc and Ri
  !
  IF ( KRRL >= 1 ) THEN
    IF (.NOT. OFLAT) THEN
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_1(JIJ, JK)=ZFLX(JIJ, JK)*PINV_PDXX(JIJ, JK)
        END DO
      END DO
      
      CALL MZM_PHY(D, ZWORKSHUGRAD_3D_1, ZWORKSHUGRAD_3D_2)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_3(JIJ, JK)=PRHODJ(JIJ, JK)*PAMOIST(JIJ, JK)*PSRCM(JIJ, JK)
        END DO
      END DO
      
      CALL MXM_PHY(D, ZWORKSHUGRAD_3D_3, ZWORKSHUGRAD_3D_4)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_5(JIJ, JK)=PRHODJ(JIJ, JK)*PAMOIST(JIJ, JK)*PSRCM(JIJ, JK)
        END DO
      END DO
      
      CALL MZM_PHY(D, ZWORKSHUGRAD_3D_5, ZWORKSHUGRAD_3D_6)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_7(JIJ, JK)=PDZX(JIJ, JK)*(ZWORKSHUGRAD_3D_2(JIJ, JK))
        END DO
      END DO
      
      CALL MXF_PHY(D, ZWORKSHUGRAD_3D_7, ZWORKSHUGRAD_3D_8)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_9(JIJ, JK)=ZWORKSHUGRAD_3D_4(JIJ, JK)*ZFLX(JIJ, JK)
        END DO
      END DO
      
      CALL MXF_PHY(D, ZWORKSHUGRAD_3D_9, ZWORKSHUGRAD_3D_10)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_11(JIJ, JK)=ZWORKSHUGRAD_3D_6(JIJ, JK)*ZWORKSHUGRAD_3D_8(JIJ, JK)
        END DO
      END DO
      
      CALL MZF_PHY(D, ZWORKSHUGRAD_3D_11, ZWORKSHUGRAD_3D_12)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZFLXC(JIJ, JK) = ZFLXC(JIJ, JK)            &
                + 2.*( ZWORKSHUGRAD_3D_10(JIJ, JK)                     &
                      +ZWORKSHUGRAD_3D_12(JIJ, JK)&
                     )
        END DO
      END DO
      
      IF ( KRRI >= 1 ) THEN
        
        DO JK=1, D%NKT
          DO JIJ=1, D%NIJT
            ZWORKSHUGRAD_3D_1(JIJ, JK)=ZFLX(JIJ, JK)*PINV_PDXX(JIJ, JK)
          END DO
        END DO
        
        CALL MZM_PHY(D, ZWORKSHUGRAD_3D_1, ZWORKSHUGRAD_3D_2)
        
        DO JK=1, D%NKT
          DO JIJ=1, D%NIJT
            ZWORKSHUGRAD_3D_3(JIJ, JK)=PRHODJ(JIJ, JK)*PAMOIST(JIJ, JK)*PSRCM(JIJ, JK)
          END DO
        END DO
        
        CALL MXM_PHY(D, ZWORKSHUGRAD_3D_3, ZWORKSHUGRAD_3D_4)
        
        DO JK=1, D%NKT
          DO JIJ=1, D%NIJT
            ZWORKSHUGRAD_3D_5(JIJ, JK)=PDZX(JIJ, JK)*(ZWORKSHUGRAD_3D_2(JIJ, JK))
          END DO
        END DO
        
        CALL MXF_PHY(D, ZWORKSHUGRAD_3D_5, ZWORKSHUGRAD_3D_6)
        
        DO JK=1, D%NKT
          DO JIJ=1, D%NIJT
            ZWORKSHUGRAD_3D_7(JIJ, JK)=PRHODJ(JIJ, JK)*PAMOIST(JIJ, JK)*PSRCM(JIJ, JK)
          END DO
        END DO
        
        CALL MZM_PHY(D, ZWORKSHUGRAD_3D_7, ZWORKSHUGRAD_3D_8)
        
        DO JK=1, D%NKT
          DO JIJ=1, D%NIJT
            ZWORKSHUGRAD_3D_9(JIJ, JK)=ZWORKSHUGRAD_3D_4(JIJ, JK)*ZFLX(JIJ, JK)*PINV_PDXX(JIJ, JK)
          END DO
        END DO
        
        CALL DXF_PHY(D, ZWORKSHUGRAD_3D_9, ZWORKSHUGRAD_3D_10)
        
        DO JK=1, D%NKT
          DO JIJ=1, D%NIJT
            ZWORKSHUGRAD_3D_11(JIJ, JK)=ZWORKSHUGRAD_3D_8(JIJ, JK)*ZWORKSHUGRAD_3D_6(JIJ, JK)&
                                               *PINV_PDZZ(JIJ, JK)
          END DO
        END DO
        
        CALL DZF_PHY(D, ZWORKSHUGRAD_3D_11, ZWORKSHUGRAD_3D_12)
        
        DO JK=1, D%NKT
          DO JIJ=1, D%NIJT
            PRRS(JIJ, JK, 2) = PRRS(JIJ, JK, 2) +  2. *                                  &
            (- ZWORKSHUGRAD_3D_10(JIJ, JK)                   &
             + ZWORKSHUGRAD_3D_12(JIJ, JK)                        &
            )*(1.0-PFRAC_ICE(JIJ, JK))
          END DO
        END DO
        
        
        DO JK=1, D%NKT
          DO JIJ=1, D%NIJT
            ZWORKSHUGRAD_3D_1(JIJ, JK)=ZFLX(JIJ, JK)*PINV_PDXX(JIJ, JK)
          END DO
        END DO
        
        CALL MZM_PHY(D, ZWORKSHUGRAD_3D_1, ZWORKSHUGRAD_3D_2)
        
        DO JK=1, D%NKT
          DO JIJ=1, D%NIJT
            ZWORKSHUGRAD_3D_3(JIJ, JK)=PRHODJ(JIJ, JK)*PAMOIST(JIJ, JK)*PSRCM(JIJ, JK)
          END DO
        END DO
        
        CALL MXM_PHY(D, ZWORKSHUGRAD_3D_3, ZWORKSHUGRAD_3D_4)
        
        DO JK=1, D%NKT
          DO JIJ=1, D%NIJT
            ZWORKSHUGRAD_3D_5(JIJ, JK)=PDZX(JIJ, JK)*(ZWORKSHUGRAD_3D_2(JIJ, JK))
          END DO
        END DO
        
        CALL MXF_PHY(D, ZWORKSHUGRAD_3D_5, ZWORKSHUGRAD_3D_6)
        
        DO JK=1, D%NKT
          DO JIJ=1, D%NIJT
            ZWORKSHUGRAD_3D_7(JIJ, JK)=PRHODJ(JIJ, JK)*PAMOIST(JIJ, JK)*PSRCM(JIJ, JK)
          END DO
        END DO
        
        CALL MZM_PHY(D, ZWORKSHUGRAD_3D_7, ZWORKSHUGRAD_3D_8)
        
        DO JK=1, D%NKT
          DO JIJ=1, D%NIJT
            ZWORKSHUGRAD_3D_9(JIJ, JK)=ZWORKSHUGRAD_3D_4(JIJ, JK)*ZFLX(JIJ, JK)*PINV_PDXX(JIJ, JK)
          END DO
        END DO
        
        CALL DXF_PHY(D, ZWORKSHUGRAD_3D_9, ZWORKSHUGRAD_3D_10)
        
        DO JK=1, D%NKT
          DO JIJ=1, D%NIJT
            ZWORKSHUGRAD_3D_11(JIJ, JK)=ZWORKSHUGRAD_3D_8(JIJ, JK)*ZWORKSHUGRAD_3D_6(JIJ, JK)&
                                               *PINV_PDZZ(JIJ, JK)
          END DO
        END DO
        
        CALL DZF_PHY(D, ZWORKSHUGRAD_3D_11, ZWORKSHUGRAD_3D_12)
        
        DO JK=1, D%NKT
          DO JIJ=1, D%NIJT
            PRRS(JIJ, JK, 2) = PRRS(JIJ, JK, 2) +  2. *                                  &
            (- ZWORKSHUGRAD_3D_10(JIJ, JK)                   &
             + ZWORKSHUGRAD_3D_12(JIJ, JK)                        &
            )*PFRAC_ICE(JIJ, JK)
          END DO
        END DO
        
      ELSE
        
        DO JK=1, D%NKT
          DO JIJ=1, D%NIJT
            ZWORKSHUGRAD_3D_1(JIJ, JK)=ZFLX(JIJ, JK)*PINV_PDXX(JIJ, JK)
          END DO
        END DO
        
        CALL MZM_PHY(D, ZWORKSHUGRAD_3D_1, ZWORKSHUGRAD_3D_2)
        
        DO JK=1, D%NKT
          DO JIJ=1, D%NIJT
            ZWORKSHUGRAD_3D_3(JIJ, JK)=PRHODJ(JIJ, JK)*PAMOIST(JIJ, JK)*PSRCM(JIJ, JK)
          END DO
        END DO
        
        CALL MXM_PHY(D, ZWORKSHUGRAD_3D_3, ZWORKSHUGRAD_3D_4)
        
        DO JK=1, D%NKT
          DO JIJ=1, D%NIJT
            ZWORKSHUGRAD_3D_5(JIJ, JK)=PDZX(JIJ, JK)*(ZWORKSHUGRAD_3D_2(JIJ, JK))
          END DO
        END DO
        
        CALL MXF_PHY(D, ZWORKSHUGRAD_3D_5, ZWORKSHUGRAD_3D_6)
        
        DO JK=1, D%NKT
          DO JIJ=1, D%NIJT
            ZWORKSHUGRAD_3D_7(JIJ, JK)=PRHODJ(JIJ, JK)*PAMOIST(JIJ, JK)*PSRCM(JIJ, JK)
          END DO
        END DO
        
        CALL MZM_PHY(D, ZWORKSHUGRAD_3D_7, ZWORKSHUGRAD_3D_8)
        
        DO JK=1, D%NKT
          DO JIJ=1, D%NIJT
            ZWORKSHUGRAD_3D_9(JIJ, JK)=ZWORKSHUGRAD_3D_4(JIJ, JK)*ZFLX(JIJ, JK)*PINV_PDXX(JIJ, JK)
          END DO
        END DO
        
        CALL DXF_PHY(D, ZWORKSHUGRAD_3D_9, ZWORKSHUGRAD_3D_10)
        
        DO JK=1, D%NKT
          DO JIJ=1, D%NIJT
            ZWORKSHUGRAD_3D_11(JIJ, JK)=ZWORKSHUGRAD_3D_8(JIJ, JK)*ZWORKSHUGRAD_3D_6(JIJ, JK)&
                                               *PINV_PDZZ(JIJ, JK)
          END DO
        END DO
        
        CALL DZF_PHY(D, ZWORKSHUGRAD_3D_11, ZWORKSHUGRAD_3D_12)
        
        DO JK=1, D%NKT
          DO JIJ=1, D%NIJT
            PRRS(JIJ, JK, 2) = PRRS(JIJ, JK, 2) +  2. *                                  &
            (- ZWORKSHUGRAD_3D_10(JIJ, JK)                   &
             + ZWORKSHUGRAD_3D_12(JIJ, JK)                        &
            )
          END DO
        END DO
        
      END IF
    ELSE
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_1(JIJ, JK)=PRHODJ(JIJ, JK)*PAMOIST(JIJ, JK)*PSRCM(JIJ, JK)
        END DO
      END DO
      
      CALL MXM_PHY(D, ZWORKSHUGRAD_3D_1, ZWORKSHUGRAD_3D_2)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_3(JIJ, JK)=ZWORKSHUGRAD_3D_2(JIJ, JK)*ZFLX(JIJ, JK)
        END DO
      END DO
      
      CALL MXF_PHY(D, ZWORKSHUGRAD_3D_3, ZWORKSHUGRAD_3D_4)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZFLXC(JIJ, JK) = ZFLXC(JIJ, JK) + 2.*ZWORKSHUGRAD_3D_4(JIJ, JK)
        END DO
      END DO
      
      IF ( KRRI >= 1 ) THEN
        
        DO JK=1, D%NKT
          DO JIJ=1, D%NIJT
            ZWORKSHUGRAD_3D_1(JIJ, JK)=PRHODJ(JIJ, JK)*PAMOIST(JIJ, JK)*PSRCM(JIJ, JK)
          END DO
        END DO
        
        CALL MXM_PHY(D, ZWORKSHUGRAD_3D_1, ZWORKSHUGRAD_3D_2)
        
        DO JK=1, D%NKT
          DO JIJ=1, D%NIJT
            ZWORKSHUGRAD_3D_3(JIJ, JK)=ZWORKSHUGRAD_3D_2(JIJ, JK)*ZFLX(JIJ, JK)*PINV_PDXX(JIJ, JK)
          END DO
        END DO
        
        CALL DXF_PHY(D, ZWORKSHUGRAD_3D_3, ZWORKSHUGRAD_3D_4)
        
        DO JK=1, D%NKT
          DO JIJ=1, D%NIJT
            PRRS(JIJ, JK, 2) = PRRS(JIJ, JK, 2) -  2. *                                  &
            ZWORKSHUGRAD_3D_4(JIJ, JK)*(1.0-PFRAC_ICE(JIJ, JK))
          END DO
        END DO
        
        
        DO JK=1, D%NKT
          DO JIJ=1, D%NIJT
            ZWORKSHUGRAD_3D_1(JIJ, JK)=PRHODJ(JIJ, JK)*PAMOIST(JIJ, JK)*PSRCM(JIJ, JK)
          END DO
        END DO
        
        CALL MXM_PHY(D, ZWORKSHUGRAD_3D_1, ZWORKSHUGRAD_3D_2)
        
        DO JK=1, D%NKT
          DO JIJ=1, D%NIJT
            ZWORKSHUGRAD_3D_3(JIJ, JK)=ZWORKSHUGRAD_3D_2(JIJ, JK)*ZFLX(JIJ, JK)*PINV_PDXX(JIJ, JK)
          END DO
        END DO
        
        CALL DXF_PHY(D, ZWORKSHUGRAD_3D_3, ZWORKSHUGRAD_3D_4)
        
        DO JK=1, D%NKT
          DO JIJ=1, D%NIJT
            PRRS(JIJ, JK, 4) = PRRS(JIJ, JK, 4) -  2. *                                  &
            ZWORKSHUGRAD_3D_4(JIJ, JK)*PFRAC_ICE(JIJ, JK)
          END DO
        END DO
        
      ELSE
        
        DO JK=1, D%NKT
          DO JIJ=1, D%NIJT
            ZWORKSHUGRAD_3D_1(JIJ, JK)=PRHODJ(JIJ, JK)*PAMOIST(JIJ, JK)*PSRCM(JIJ, JK)
          END DO
        END DO
        
        CALL MXM_PHY(D, ZWORKSHUGRAD_3D_1, ZWORKSHUGRAD_3D_2)
        
        DO JK=1, D%NKT
          DO JIJ=1, D%NIJT
            ZWORKSHUGRAD_3D_3(JIJ, JK)=ZWORKSHUGRAD_3D_2(JIJ, JK)*ZFLX(JIJ, JK)*PINV_PDXX(JIJ, JK)
          END DO
        END DO
        
        CALL DXF_PHY(D, ZWORKSHUGRAD_3D_3, ZWORKSHUGRAD_3D_4)
        
        DO JK=1, D%NKT
          DO JIJ=1, D%NIJT
            PRRS(JIJ, JK, 2) = PRRS(JIJ, JK, 2) -  2. *                                  &
            ZWORKSHUGRAD_3D_4(JIJ, JK)
          END DO
        END DO
        
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

    CALL IO_FIELD_WRITE_PHY(D,TPFILE,TZFIELD,ZFLX)
  END IF
  !
  IF (KSPLT==1 .AND. TLES%LLES_CALL) THEN
    CALL SECOND_MNH(ZTIME1)
    CALL MXF_PHY(D, ZFLX, ZWORKSHUGRAD_3D_1)
    CALL LES_MEAN_SUBGRID_PHY(D, TLES,ZWORKSHUGRAD_3D_1(1:D%NIJT, 1:D%NKT), TLES%X_LES_SUBGRID_URt ) 
    CALL GX_W_UW_PHY(D, OFLAT, PWM, PDXX, PDZZ, PDZX, ZWORKSHUGRAD_3D_1)
    CALL MZM_PHY(D, ZFLX, ZWORKSHUGRAD_3D_2)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWORKSHUGRAD_3D_3(JIJ, JK)=ZWORKSHUGRAD_3D_1(JIJ, JK)*ZWORKSHUGRAD_3D_2(JIJ, JK)
      END DO
    END DO
    
    CALL MXF_PHY(D, ZWORKSHUGRAD_3D_3, ZWORKSHUGRAD_3D_4)
    CALL MZF_PHY(D, ZWORKSHUGRAD_3D_4, ZWKLES)
    CALL LES_MEAN_SUBGRID_PHY(D, TLES,ZWKLES, TLES%X_LES_RES_ddxa_W_SBG_UaRt , .TRUE. )
    CALL GX_M_M_PHY(D, OFLAT, PTHLM, PDXX, PDZZ, PDZX, ZWORKSHUGRAD_3D_1)
    CALL MXF_PHY(D, ZFLX, ZWORKSHUGRAD_3D_2)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWKLES(JIJ, JK) = ZWORKSHUGRAD_3D_1(JIJ, JK)*ZWORKSHUGRAD_3D_2(JIJ, JK)
      END DO
    END DO
    
    CALL LES_MEAN_SUBGRID_PHY(D, TLES,ZWKLES, TLES%X_LES_RES_ddxa_Thl_SBG_UaRt , .TRUE. )
    CALL GX_M_M_PHY(D, OFLAT, PRM(:,:,1), PDXX, PDZZ, PDZX, ZWORKSHUGRAD_3D_1)
    CALL MXF_PHY(D, ZFLX, ZWORKSHUGRAD_3D_2)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWKLES(JIJ, JK) = ZWORKSHUGRAD_3D_1(JIJ, JK)*ZWORKSHUGRAD_3D_2(JIJ, JK)
      END DO
    END DO
    
    CALL LES_MEAN_SUBGRID_PHY(D, TLES,ZWKLES, TLES%X_LES_RES_ddxa_Rt_SBG_UaRt , .TRUE. )
    CALL SECOND_MNH(ZTIME2)
    TLES%XTIME_LES = TLES%XTIME_LES + ZTIME2 - ZTIME1
  END IF
!
  !
  IF (KRRL>0 .AND. KSPLT==1 .AND. TLES%LLES_CALL) THEN
    CALL SECOND_MNH(ZTIME1)
    CALL MXF_PHY(D, ZFLXC, ZWORKSHUGRAD_3D_1)
    CALL LES_MEAN_SUBGRID_PHY(D, TLES,ZWORKSHUGRAD_3D_1(1:D%NIJT, 1:D%NKT), TLES%X_LES_SUBGRID_URc )
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
!!  ZVPTU(:,:) =                                                        &
!!    ZWORK(:,:)*MXM(ETHETA(KRR,KRRI,PTHLT,PEXNREF,PRT,PLOCPT,PSRCM)) +       &
!!     ZFLX(:,:)*MXM(EMOIST(KRR,KRRI,PTHLT,PEXNREF,PRT,PLOCPT,PSRCM))
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
!!    CALL IO_FIELD_WRITE_PHY(D,TPFILE,TZFIELD,ZVPTU)
!!  END IF
!!!
!!ELSE
!!  ZVPTU(:,:)=ZWORK(:,:)
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
  CALL MYM_PHY(D, PK, ZWORKSHUGRAD_3D_1)
  CALL GY_M_V_PHY(D, OFLAT, PTHLM, PDYY, PDZZ, PDZY, ZWORKSHUGRAD_3D_2)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      ZFLX(JIJ, JK)     = -TURBN%XCSHF * ZWORKSHUGRAD_3D_1(JIJ, JK) * ZWORKSHUGRAD_3D_2(JIJ, JK)
    END DO
  END DO
  

  ZFLX(:,IKE+1) = ZFLX(:,IKE)

ELSE

  ZFLX(:,:)     = 0.

END IF
!
!
! Compute the flux at the first inner U-point with an uncentred vertical  
! gradient
CALL DYM2D_PHY(D, PTHLM(:,IKB), ZWORKSHUGRAD_2D_1)

DO JIJ=1, D%NIJT
  ZWORKSHUGRAD_2D_2(JIJ)=ZCOEFF(JIJ, IKB+2)*PTHLM(JIJ, IKB+2)        &
           +ZCOEFF(JIJ, IKB+1)*PTHLM(JIJ, IKB+1)       &
           +ZCOEFF(JIJ, IKB)*PTHLM(JIJ, IKB)
END DO

CALL MYM2D_PHY(D, ZWORKSHUGRAD_2D_2, ZWORKSHUGRAD_2D_3)
CALL MYM2D_PHY(D, PK(:,IKB), ZWORKSHUGRAD_2D_4)

DO JIJ=1, D%NIJT
  ZFLX(JIJ, IKB) = -TURBN%XCSHF * ZWORKSHUGRAD_2D_4(JIJ) *          &
    ( ZWORKSHUGRAD_2D_1(JIJ) * PINV_PDYY(JIJ, IKB)           &
     -ZWORKSHUGRAD_2D_3(JIJ)     &
          *0.5* ( PDZY(JIJ, IKB+1)+PDZY(JIJ, IKB))       &
          * PINV_PDYY(JIJ, IKB) )
END DO

! extrapolates the flux under the ground so that the vertical average with 
! the IKB flux gives the ground value  ( warning the tangential surface
! flux has been set to 0 for the moment !!  to be improved )

DO JIJ=1, D%NIJT
  ZWORKSHUGRAD_2D_1(JIJ)=PSFTHM(JIJ)* PDIRCOSYW(JIJ)
END DO

CALL MYM2D_PHY(D, ZWORKSHUGRAD_2D_1, ZWORKSHUGRAD_2D_2)

DO JIJ=1, D%NIJT
  ZFLX(JIJ, IKB-1) = 2. * ZWORKSHUGRAD_2D_2(JIJ) &
                         - ZFLX(JIJ, IKB)
END DO

!
! Add this source to the Theta_l sources
!
IF (.NOT. O2D) THEN 
  IF (.NOT. OFLAT) THEN
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWORKSHUGRAD_3D_1(JIJ, JK)=ZFLX(JIJ, JK) * PINV_PDYY(JIJ, JK)
      END DO
    END DO
    
    CALL MZM_PHY(D, ZWORKSHUGRAD_3D_1, ZWORKSHUGRAD_3D_2)
    CALL MYM_PHY(D, PRHODJ, ZWORKSHUGRAD_3D_3)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWORKSHUGRAD_3D_4(JIJ, JK)=PDZY(JIJ, JK)*(ZWORKSHUGRAD_3D_2(JIJ, JK))
      END DO
    END DO
    
    CALL MYF_PHY(D, ZWORKSHUGRAD_3D_4, ZWORKSHUGRAD_3D_5)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWORKSHUGRAD_3D_6(JIJ, JK)=ZWORKSHUGRAD_3D_3(JIJ, JK) * ZFLX(JIJ, JK) * PINV_PDYY(JIJ, JK)
      END DO
    END DO
    
    CALL DYF_PHY(D, ZWORKSHUGRAD_3D_6, ZWORKSHUGRAD_3D_7)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWORKSHUGRAD_3D_8(JIJ, JK)=PMZM_PRHODJ(JIJ, JK) *ZWORKSHUGRAD_3D_5(JIJ, JK) * PINV_PDZZ(JIJ, JK)
      END DO
    END DO
    
    CALL DZF_PHY(D, ZWORKSHUGRAD_3D_8, ZWORKSHUGRAD_3D_9)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        PRTHLS(JIJ, JK) =  PRTHLS(JIJ, JK)                                                         &
                      - ZWORKSHUGRAD_3D_7(JIJ, JK)                           &
                      + ZWORKSHUGRAD_3D_9(JIJ, JK)
      END DO
    END DO
    
  ELSE
    CALL MYM_PHY(D, PRHODJ, ZWORKSHUGRAD_3D_1)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWORKSHUGRAD_3D_2(JIJ, JK)=ZWORKSHUGRAD_3D_1(JIJ, JK) * ZFLX(JIJ, JK) * PINV_PDYY(JIJ, JK)
      END DO
    END DO
    
    CALL DYF_PHY(D, ZWORKSHUGRAD_3D_2, ZWORKSHUGRAD_3D_3)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        PRTHLS(JIJ, JK) =  PRTHLS(JIJ, JK) - ZWORKSHUGRAD_3D_3(JIJ, JK)
      END DO
    END DO
    
  END IF
END IF
!
! Compute the equivalent tendancy for Rc and Ri
!
!IF ( TURBN%LSUBG_COND .AND. KRRL > 0 .AND. .NOT. O2D) THEN
IF ( KRRL >= 1 .AND. .NOT. O2D) THEN
  IF (.NOT. OFLAT) THEN
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWORKSHUGRAD_3D_1(JIJ, JK)=ZFLX(JIJ, JK)*PINV_PDYY(JIJ, JK)
      END DO
    END DO
    
    CALL MZM_PHY(D, ZWORKSHUGRAD_3D_1, ZWORKSHUGRAD_3D_2)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWORKSHUGRAD_3D_3(JIJ, JK)=PRHODJ(JIJ, JK)*PATHETA(JIJ, JK)*PSRCM(JIJ, JK)
      END DO
    END DO
    
    CALL MYM_PHY(D, ZWORKSHUGRAD_3D_3, ZWORKSHUGRAD_3D_4)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWORKSHUGRAD_3D_5(JIJ, JK)=PRHODJ(JIJ, JK)*PATHETA(JIJ, JK)*PSRCM(JIJ, JK)
      END DO
    END DO
    
    CALL MZM_PHY(D, ZWORKSHUGRAD_3D_5, ZWORKSHUGRAD_3D_6)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWORKSHUGRAD_3D_7(JIJ, JK)=PDZY(JIJ, JK)*(ZWORKSHUGRAD_3D_2(JIJ, JK))
      END DO
    END DO
    
    CALL MYF_PHY(D, ZWORKSHUGRAD_3D_7, ZWORKSHUGRAD_3D_8)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWORKSHUGRAD_3D_9(JIJ, JK)=ZWORKSHUGRAD_3D_4(JIJ, JK)*ZFLX(JIJ, JK)
      END DO
    END DO
    
    CALL MYF_PHY(D, ZWORKSHUGRAD_3D_9, ZWORKSHUGRAD_3D_10)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWORKSHUGRAD_3D_11(JIJ, JK)=ZWORKSHUGRAD_3D_6(JIJ, JK)*ZWORKSHUGRAD_3D_8(JIJ, JK)
      END DO
    END DO
    
    CALL MZF_PHY(D, ZWORKSHUGRAD_3D_11, ZWORKSHUGRAD_3D_12)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZFLXC(JIJ, JK) = 2.*( ZWORKSHUGRAD_3D_10(JIJ, JK)                       &
                    +ZWORKSHUGRAD_3D_12(JIJ, JK)&
                   )
      END DO
    END DO
    
    IF ( KRRI >= 1 ) THEN
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_1(JIJ, JK)=ZFLX(JIJ, JK)*PINV_PDYY(JIJ, JK)
        END DO
      END DO
      
      CALL MZM_PHY(D, ZWORKSHUGRAD_3D_1, ZWORKSHUGRAD_3D_2)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_3(JIJ, JK)=PRHODJ(JIJ, JK)*PATHETA(JIJ, JK)*PSRCM(JIJ, JK)
        END DO
      END DO
      
      CALL MYM_PHY(D, ZWORKSHUGRAD_3D_3, ZWORKSHUGRAD_3D_4)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_5(JIJ, JK)=PDZY(JIJ, JK)*(ZWORKSHUGRAD_3D_2(JIJ, JK))
        END DO
      END DO
      
      CALL MYF_PHY(D, ZWORKSHUGRAD_3D_5, ZWORKSHUGRAD_3D_6)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_7(JIJ, JK)=PRHODJ(JIJ, JK)*PATHETA(JIJ, JK)*PSRCM(JIJ, JK)
        END DO
      END DO
      
      CALL MZM_PHY(D, ZWORKSHUGRAD_3D_7, ZWORKSHUGRAD_3D_8)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_9(JIJ, JK)=ZWORKSHUGRAD_3D_4(JIJ, JK)*ZFLX(JIJ, JK)*PINV_PDYY(JIJ, JK)
        END DO
      END DO
      
      CALL DYF_PHY(D, ZWORKSHUGRAD_3D_9, ZWORKSHUGRAD_3D_10)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_11(JIJ, JK)=ZWORKSHUGRAD_3D_8(JIJ, JK)*ZWORKSHUGRAD_3D_6(JIJ, JK)&
                                               *PINV_PDZZ(JIJ, JK)
        END DO
      END DO
      
      CALL DZF_PHY(D, ZWORKSHUGRAD_3D_11, ZWORKSHUGRAD_3D_12)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          PRRS(JIJ, JK, 2) = PRRS(JIJ, JK, 2) + 2. *                                     &
            (- ZWORKSHUGRAD_3D_10(JIJ, JK)                   &
             + ZWORKSHUGRAD_3D_12(JIJ, JK)                        &
            )*(1.0-PFRAC_ICE(JIJ, JK))
        END DO
      END DO
      
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_1(JIJ, JK)=ZFLX(JIJ, JK)*PINV_PDYY(JIJ, JK)
        END DO
      END DO
      
      CALL MZM_PHY(D, ZWORKSHUGRAD_3D_1, ZWORKSHUGRAD_3D_2)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_3(JIJ, JK)=PRHODJ(JIJ, JK)*PATHETA(JIJ, JK)*PSRCM(JIJ, JK)
        END DO
      END DO
      
      CALL MYM_PHY(D, ZWORKSHUGRAD_3D_3, ZWORKSHUGRAD_3D_4)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_5(JIJ, JK)=PDZY(JIJ, JK)*(ZWORKSHUGRAD_3D_2(JIJ, JK))
        END DO
      END DO
      
      CALL MYF_PHY(D, ZWORKSHUGRAD_3D_5, ZWORKSHUGRAD_3D_6)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_7(JIJ, JK)=PRHODJ(JIJ, JK)*PATHETA(JIJ, JK)*PSRCM(JIJ, JK)
        END DO
      END DO
      
      CALL MZM_PHY(D, ZWORKSHUGRAD_3D_7, ZWORKSHUGRAD_3D_8)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_9(JIJ, JK)=ZWORKSHUGRAD_3D_4(JIJ, JK)*ZFLX(JIJ, JK)*PINV_PDYY(JIJ, JK)
        END DO
      END DO
      
      CALL DYF_PHY(D, ZWORKSHUGRAD_3D_9, ZWORKSHUGRAD_3D_10)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_11(JIJ, JK)=ZWORKSHUGRAD_3D_8(JIJ, JK)*ZWORKSHUGRAD_3D_6(JIJ, JK)&
                                               *PINV_PDZZ(JIJ, JK)
        END DO
      END DO
      
      CALL DZF_PHY(D, ZWORKSHUGRAD_3D_11, ZWORKSHUGRAD_3D_12)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          PRRS(JIJ, JK, 4) = PRRS(JIJ, JK, 4) + 2. *                                     &
            (- ZWORKSHUGRAD_3D_10(JIJ, JK)                   &
             + ZWORKSHUGRAD_3D_12(JIJ, JK)                        &
            )*PFRAC_ICE(JIJ, JK)
        END DO
      END DO
      
    ELSE
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_1(JIJ, JK)=ZFLX(JIJ, JK)*PINV_PDYY(JIJ, JK)
        END DO
      END DO
      
      CALL MZM_PHY(D, ZWORKSHUGRAD_3D_1, ZWORKSHUGRAD_3D_2)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_3(JIJ, JK)=PRHODJ(JIJ, JK)*PATHETA(JIJ, JK)*PSRCM(JIJ, JK)
        END DO
      END DO
      
      CALL MYM_PHY(D, ZWORKSHUGRAD_3D_3, ZWORKSHUGRAD_3D_4)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_5(JIJ, JK)=PDZY(JIJ, JK)*(ZWORKSHUGRAD_3D_2(JIJ, JK))
        END DO
      END DO
      
      CALL MYF_PHY(D, ZWORKSHUGRAD_3D_5, ZWORKSHUGRAD_3D_6)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_7(JIJ, JK)=PRHODJ(JIJ, JK)*PATHETA(JIJ, JK)*PSRCM(JIJ, JK)
        END DO
      END DO
      
      CALL MZM_PHY(D, ZWORKSHUGRAD_3D_7, ZWORKSHUGRAD_3D_8)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_9(JIJ, JK)=ZWORKSHUGRAD_3D_4(JIJ, JK)*ZFLX(JIJ, JK)*PINV_PDYY(JIJ, JK)
        END DO
      END DO
      
      CALL DYF_PHY(D, ZWORKSHUGRAD_3D_9, ZWORKSHUGRAD_3D_10)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_11(JIJ, JK)=ZWORKSHUGRAD_3D_8(JIJ, JK)*ZWORKSHUGRAD_3D_6(JIJ, JK)&
                                               *PINV_PDZZ(JIJ, JK)
        END DO
      END DO
      
      CALL DZF_PHY(D, ZWORKSHUGRAD_3D_11, ZWORKSHUGRAD_3D_12)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          PRRS(JIJ, JK, 2) = PRRS(JIJ, JK, 2) + 2. *                                     &
            (- ZWORKSHUGRAD_3D_10(JIJ, JK)                   &
             + ZWORKSHUGRAD_3D_12(JIJ, JK)                        &
            )
        END DO
      END DO
      
    END IF
  ELSE
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWORKSHUGRAD_3D_1(JIJ, JK)=PRHODJ(JIJ, JK)*PATHETA(JIJ, JK)*PSRCM(JIJ, JK)
      END DO
    END DO
    
    CALL MYM_PHY(D, ZWORKSHUGRAD_3D_1, ZWORKSHUGRAD_3D_2)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWORKSHUGRAD_3D_3(JIJ, JK)=ZWORKSHUGRAD_3D_2(JIJ, JK)*ZFLX(JIJ, JK)
      END DO
    END DO
    
    CALL MYF_PHY(D, ZWORKSHUGRAD_3D_3, ZWORKSHUGRAD_3D_4)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZFLXC(JIJ, JK) = 2.*ZWORKSHUGRAD_3D_4(JIJ, JK)
      END DO
    END DO
    
    IF ( KRRI >= 1 ) THEN
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_1(JIJ, JK)=PRHODJ(JIJ, JK)*PATHETA(JIJ, JK)*PSRCM(JIJ, JK)
        END DO
      END DO
      
      CALL MYM_PHY(D, ZWORKSHUGRAD_3D_1, ZWORKSHUGRAD_3D_2)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_3(JIJ, JK)=ZWORKSHUGRAD_3D_2(JIJ, JK)*ZFLX(JIJ, JK)*PINV_PDYY(JIJ, JK)
        END DO
      END DO
      
      CALL DYF_PHY(D, ZWORKSHUGRAD_3D_3, ZWORKSHUGRAD_3D_4)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          PRRS(JIJ, JK, 2) = PRRS(JIJ, JK, 2) - 2. *                                     &
            ZWORKSHUGRAD_3D_4(JIJ, JK)*(1.0-PFRAC_ICE(JIJ, JK))
        END DO
      END DO
      
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_1(JIJ, JK)=PRHODJ(JIJ, JK)*PATHETA(JIJ, JK)*PSRCM(JIJ, JK)
        END DO
      END DO
      
      CALL MYM_PHY(D, ZWORKSHUGRAD_3D_1, ZWORKSHUGRAD_3D_2)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_3(JIJ, JK)=ZWORKSHUGRAD_3D_2(JIJ, JK)*ZFLX(JIJ, JK)*PINV_PDYY(JIJ, JK)
        END DO
      END DO
      
      CALL DYF_PHY(D, ZWORKSHUGRAD_3D_3, ZWORKSHUGRAD_3D_4)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          PRRS(JIJ, JK, 4) = PRRS(JIJ, JK, 4) - 2. *                                     &
            ZWORKSHUGRAD_3D_4(JIJ, JK)*PFRAC_ICE(JIJ, JK)
        END DO
      END DO
      
    ELSE
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_1(JIJ, JK)=PRHODJ(JIJ, JK)*PATHETA(JIJ, JK)*PSRCM(JIJ, JK)
        END DO
      END DO
      
      CALL MYM_PHY(D, ZWORKSHUGRAD_3D_1, ZWORKSHUGRAD_3D_2)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_3(JIJ, JK)=ZWORKSHUGRAD_3D_2(JIJ, JK)*ZFLX(JIJ, JK)*PINV_PDYY(JIJ, JK)
        END DO
      END DO
      
      CALL DYF_PHY(D, ZWORKSHUGRAD_3D_3, ZWORKSHUGRAD_3D_4)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          PRRS(JIJ, JK, 2) = PRRS(JIJ, JK, 2) - 2. *                                     &
            ZWORKSHUGRAD_3D_4(JIJ, JK)
        END DO
      END DO
      
    END IF
  END IF
END IF
!
!! stores this flux in ZWORK to compute later <V' VPT'>
!!ZWORK(:,:) = ZFLX(:,:) 
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

  CALL IO_FIELD_WRITE_PHY(D,TPFILE,TZFIELD,ZFLX)
END IF
!
IF (KSPLT==1 .AND. TLES%LLES_CALL) THEN
  CALL SECOND_MNH(ZTIME1)
  CALL MYF_PHY(D, ZFLX, ZWORKSHUGRAD_3D_1)
  CALL LES_MEAN_SUBGRID_PHY(D, TLES,ZWORKSHUGRAD_3D_1(1:D%NIJT, 1:D%NKT), TLES%X_LES_SUBGRID_VThl ) 
  CALL GY_W_VW_PHY(D, OFLAT, PWM, PDYY, PDZZ, PDZY, ZWORKSHUGRAD_3D_1)
  CALL MZM_PHY(D, ZFLX, ZWORKSHUGRAD_3D_2)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      ZWORKSHUGRAD_3D_3(JIJ, JK)=ZWORKSHUGRAD_3D_1(JIJ, JK)*ZWORKSHUGRAD_3D_2(JIJ, JK)
    END DO
  END DO
  
  CALL MYF_PHY(D, ZWORKSHUGRAD_3D_3, ZWORKSHUGRAD_3D_4)
  CALL MZF_PHY(D, ZWORKSHUGRAD_3D_4, ZWKLES)
  CALL LES_MEAN_SUBGRID_PHY(D, TLES,ZWKLES, TLES%X_LES_RES_ddxa_W_SBG_UaThl , .TRUE. )
  CALL GY_M_M_PHY(D, OFLAT, PTHLM, PDYY, PDZZ, PDZY, ZWORKSHUGRAD_3D_1)
  CALL MYF_PHY(D, ZFLX, ZWORKSHUGRAD_3D_2)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      ZWKLES(JIJ, JK) = ZWORKSHUGRAD_3D_1(JIJ, JK)*ZWORKSHUGRAD_3D_2(JIJ, JK)
    END DO
  END DO
  
  CALL LES_MEAN_SUBGRID_PHY(D, TLES,ZWKLES, TLES%X_LES_RES_ddxa_Thl_SBG_UaThl , .TRUE. )
  IF (KRR>=1) THEN
    CALL GY_M_M_PHY(D, OFLAT, PRM(:,:,1), PDYY, PDZZ, PDZY, ZWORKSHUGRAD_3D_1)
    CALL MYF_PHY(D, ZFLX, ZWORKSHUGRAD_3D_2)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWKLES(JIJ, JK) = ZWORKSHUGRAD_3D_1(JIJ, JK)*ZWORKSHUGRAD_3D_2(JIJ, JK)
      END DO
    END DO
    
    CALL LES_MEAN_SUBGRID_PHY(D, TLES,ZWKLES, TLES%X_LES_RES_ddxa_Rt_SBG_UaThl , .TRUE. )
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
    CALL MYM_PHY(D, PK, ZWORKSHUGRAD_3D_1)
    CALL GY_M_V_PHY(D, OFLAT, PRM(:,:,1), PDYY, PDZZ, PDZY, ZWORKSHUGRAD_3D_2)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZFLX(JIJ, JK)     = -TURBN%XCHF * ZWORKSHUGRAD_3D_1(JIJ, JK) * ZWORKSHUGRAD_3D_2(JIJ, JK)
      END DO
    END DO
    

    ZFLX(:,IKE+1) = ZFLX(:,IKE)

  ELSE

    ZFLX(:,:)     = 0.

  END IF
!
! Compute the flux at the first inner U-point with an uncentred vertical  
! gradient
  CALL DYM2D_PHY(D, PRM(:,IKB,1), ZWORKSHUGRAD_2D_1)
  
  DO JIJ=1, D%NIJT
    ZWORKSHUGRAD_2D_2(JIJ)=ZCOEFF(JIJ, IKB+2)*PRM(JIJ, IKB+2, 1)        &
             +ZCOEFF(JIJ, IKB+1)*PRM(JIJ, IKB+1, 1)       &
             +ZCOEFF(JIJ, IKB)*PRM(JIJ, IKB, 1)
  END DO
  
  CALL MYM2D_PHY(D, ZWORKSHUGRAD_2D_2, ZWORKSHUGRAD_2D_3)
  CALL MYM2D_PHY(D, PK(:,IKB), ZWORKSHUGRAD_2D_4)
  
  DO JIJ=1, D%NIJT
    ZFLX(JIJ, IKB) = -TURBN%XCHF * ZWORKSHUGRAD_2D_4(JIJ) *           &
      ( ZWORKSHUGRAD_2D_1(JIJ) * PINV_PDYY(JIJ, IKB)           &
       -ZWORKSHUGRAD_2D_3(JIJ)     &
             *0.5* ( PDZY(JIJ, IKB+1)+PDZY(JIJ, IKB))      &
            * PINV_PDYY(JIJ, IKB) )
  END DO
  
! extrapolates the flux under the ground so that the vertical average with 
! the IKB flux gives the ground value  ( warning the tangential surface
! flux has been set to 0 for the moment !!  to be improved )
  
  DO JIJ=1, D%NIJT
    ZWORKSHUGRAD_2D_1(JIJ)=PSFRM(JIJ)* PDIRCOSYW(JIJ)
  END DO
  
  CALL MYM2D_PHY(D, ZWORKSHUGRAD_2D_1, ZWORKSHUGRAD_2D_2)
  
  DO JIJ=1, D%NIJT
    ZFLX(JIJ, IKB-1) = 2. * ZWORKSHUGRAD_2D_2(JIJ) &
                         - ZFLX(JIJ, IKB)
  END DO
  
  !
  ! Add this source to the conservative mixing ratio sources
  !
  IF (.NOT. O2D) THEN 
    IF (.NOT. OFLAT) THEN
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_1(JIJ, JK)=ZFLX(JIJ, JK) * PINV_PDYY(JIJ, JK)
        END DO
      END DO
      
      CALL MZM_PHY(D, ZWORKSHUGRAD_3D_1, ZWORKSHUGRAD_3D_2)
      CALL MYM_PHY(D, PRHODJ, ZWORKSHUGRAD_3D_3)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_4(JIJ, JK)=PDZY(JIJ, JK)*(ZWORKSHUGRAD_3D_2(JIJ, JK))
        END DO
      END DO
      
      CALL MYF_PHY(D, ZWORKSHUGRAD_3D_4, ZWORKSHUGRAD_3D_5)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_6(JIJ, JK)=ZWORKSHUGRAD_3D_3(JIJ, JK) * ZFLX(JIJ, JK) * PINV_PDYY(JIJ, JK)
        END DO
      END DO
      
      CALL DYF_PHY(D, ZWORKSHUGRAD_3D_6, ZWORKSHUGRAD_3D_7)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_8(JIJ, JK)=PMZM_PRHODJ(JIJ, JK) *ZWORKSHUGRAD_3D_5(JIJ, JK) * PINV_PDZZ(JIJ, JK)
        END DO
      END DO
      
      CALL DZF_PHY(D, ZWORKSHUGRAD_3D_8, ZWORKSHUGRAD_3D_9)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          PRRS(JIJ, JK, 1) = PRRS(JIJ, JK, 1)                                              &
                        - ZWORKSHUGRAD_3D_7(JIJ, JK)                           &

                        + ZWORKSHUGRAD_3D_9(JIJ, JK)
        END DO
      END DO
      
    ELSE
      CALL MYM_PHY(D, PRHODJ, ZWORKSHUGRAD_3D_1)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_2(JIJ, JK)=ZWORKSHUGRAD_3D_1(JIJ, JK) * ZFLX(JIJ, JK) * PINV_PDYY(JIJ, JK)
        END DO
      END DO
      
      CALL DYF_PHY(D, ZWORKSHUGRAD_3D_2, ZWORKSHUGRAD_3D_3)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          PRRS(JIJ, JK, 1) = PRRS(JIJ, JK, 1) - ZWORKSHUGRAD_3D_3(JIJ, JK)
        END DO
      END DO
      
    END IF
  END IF
  !
  ! Compute the equivalent tendancy for Rc and Ri
  !
  IF ( KRRL >= 1 .AND. .NOT. O2D) THEN   ! Sub-grid condensation
    IF (.NOT. OFLAT) THEN
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_1(JIJ, JK)=ZFLX(JIJ, JK)*PINV_PDYY(JIJ, JK)
        END DO
      END DO
      
      CALL MZM_PHY(D, ZWORKSHUGRAD_3D_1, ZWORKSHUGRAD_3D_2)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_3(JIJ, JK)=PRHODJ(JIJ, JK)*PAMOIST(JIJ, JK)*PSRCM(JIJ, JK)
        END DO
      END DO
      
      CALL MYM_PHY(D, ZWORKSHUGRAD_3D_3, ZWORKSHUGRAD_3D_4)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_5(JIJ, JK)=PRHODJ(JIJ, JK)*PAMOIST(JIJ, JK)*PSRCM(JIJ, JK)
        END DO
      END DO
      
      CALL MZM_PHY(D, ZWORKSHUGRAD_3D_5, ZWORKSHUGRAD_3D_6)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_7(JIJ, JK)=PDZY(JIJ, JK)*(ZWORKSHUGRAD_3D_2(JIJ, JK))
        END DO
      END DO
      
      CALL MYF_PHY(D, ZWORKSHUGRAD_3D_7, ZWORKSHUGRAD_3D_8)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_9(JIJ, JK)=ZWORKSHUGRAD_3D_4(JIJ, JK)*ZFLX(JIJ, JK)
        END DO
      END DO
      
      CALL MXF_PHY(D, ZWORKSHUGRAD_3D_9, ZWORKSHUGRAD_3D_10)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_11(JIJ, JK)=ZWORKSHUGRAD_3D_6(JIJ, JK)*ZWORKSHUGRAD_3D_8(JIJ, JK)
        END DO
      END DO
      
      CALL MZF_PHY(D, ZWORKSHUGRAD_3D_11, ZWORKSHUGRAD_3D_12)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZFLXC(JIJ, JK) = ZFLXC(JIJ, JK)            &
                + 2.*( ZWORKSHUGRAD_3D_10(JIJ, JK)                     &
                    +  ZWORKSHUGRAD_3D_12(JIJ, JK)&
                     )
        END DO
      END DO
      
      IF ( KRRI >= 1 ) THEN
        
        DO JK=1, D%NKT
          DO JIJ=1, D%NIJT
            ZWORKSHUGRAD_3D_1(JIJ, JK)=ZFLX(JIJ, JK)*PINV_PDYY(JIJ, JK)
          END DO
        END DO
        
        CALL MZM_PHY(D, ZWORKSHUGRAD_3D_1, ZWORKSHUGRAD_3D_2)
        
        DO JK=1, D%NKT
          DO JIJ=1, D%NIJT
            ZWORKSHUGRAD_3D_3(JIJ, JK)=PRHODJ(JIJ, JK)*PAMOIST(JIJ, JK)*PSRCM(JIJ, JK)
          END DO
        END DO
        
        CALL MYM_PHY(D, ZWORKSHUGRAD_3D_3, ZWORKSHUGRAD_3D_4)
        
        DO JK=1, D%NKT
          DO JIJ=1, D%NIJT
            ZWORKSHUGRAD_3D_5(JIJ, JK)=PDZY(JIJ, JK)*(ZWORKSHUGRAD_3D_2(JIJ, JK))
          END DO
        END DO
        
        CALL MYF_PHY(D, ZWORKSHUGRAD_3D_5, ZWORKSHUGRAD_3D_6)
        
        DO JK=1, D%NKT
          DO JIJ=1, D%NIJT
            ZWORKSHUGRAD_3D_7(JIJ, JK)=PRHODJ(JIJ, JK)*PAMOIST(JIJ, JK)*PSRCM(JIJ, JK)
          END DO
        END DO
        
        CALL MZM_PHY(D, ZWORKSHUGRAD_3D_7, ZWORKSHUGRAD_3D_8)
        
        DO JK=1, D%NKT
          DO JIJ=1, D%NIJT
            ZWORKSHUGRAD_3D_9(JIJ, JK)=ZWORKSHUGRAD_3D_4(JIJ, JK)*ZFLX(JIJ, JK)/PDYY(JIJ, JK)
          END DO
        END DO
        
        CALL DYF_PHY(D, ZWORKSHUGRAD_3D_9, ZWORKSHUGRAD_3D_10)
        
        DO JK=1, D%NKT
          DO JIJ=1, D%NIJT
            ZWORKSHUGRAD_3D_11(JIJ, JK)=ZWORKSHUGRAD_3D_8(JIJ, JK)*ZWORKSHUGRAD_3D_6(JIJ, JK)&
                                               * PINV_PDZZ(JIJ, JK)
          END DO
        END DO
        
        CALL DZF_PHY(D, ZWORKSHUGRAD_3D_11, ZWORKSHUGRAD_3D_12)
        
        DO JK=1, D%NKT
          DO JIJ=1, D%NIJT
            PRRS(JIJ, JK, 2) = PRRS(JIJ, JK, 2) +  2. *                                  &
            (- ZWORKSHUGRAD_3D_10(JIJ, JK)                        &
             + ZWORKSHUGRAD_3D_12(JIJ, JK)                       &
            )*(1.0-PFRAC_ICE(JIJ, JK))
          END DO
        END DO
        
        
        DO JK=1, D%NKT
          DO JIJ=1, D%NIJT
            ZWORKSHUGRAD_3D_1(JIJ, JK)=ZFLX(JIJ, JK)*PINV_PDYY(JIJ, JK)
          END DO
        END DO
        
        CALL MZM_PHY(D, ZWORKSHUGRAD_3D_1, ZWORKSHUGRAD_3D_2)
        
        DO JK=1, D%NKT
          DO JIJ=1, D%NIJT
            ZWORKSHUGRAD_3D_3(JIJ, JK)=PRHODJ(JIJ, JK)*PAMOIST(JIJ, JK)*PSRCM(JIJ, JK)
          END DO
        END DO
        
        CALL MYM_PHY(D, ZWORKSHUGRAD_3D_3, ZWORKSHUGRAD_3D_4)
        
        DO JK=1, D%NKT
          DO JIJ=1, D%NIJT
            ZWORKSHUGRAD_3D_5(JIJ, JK)=PDZY(JIJ, JK)*(ZWORKSHUGRAD_3D_2(JIJ, JK))
          END DO
        END DO
        
        CALL MYF_PHY(D, ZWORKSHUGRAD_3D_5, ZWORKSHUGRAD_3D_6)
        
        DO JK=1, D%NKT
          DO JIJ=1, D%NIJT
            ZWORKSHUGRAD_3D_7(JIJ, JK)=PRHODJ(JIJ, JK)*PAMOIST(JIJ, JK)*PSRCM(JIJ, JK)
          END DO
        END DO
        
        CALL MZM_PHY(D, ZWORKSHUGRAD_3D_7, ZWORKSHUGRAD_3D_8)
        
        DO JK=1, D%NKT
          DO JIJ=1, D%NIJT
            ZWORKSHUGRAD_3D_9(JIJ, JK)=ZWORKSHUGRAD_3D_4(JIJ, JK)*ZFLX(JIJ, JK)/PDYY(JIJ, JK)
          END DO
        END DO
        
        CALL DYF_PHY(D, ZWORKSHUGRAD_3D_9, ZWORKSHUGRAD_3D_10)
        
        DO JK=1, D%NKT
          DO JIJ=1, D%NIJT
            ZWORKSHUGRAD_3D_11(JIJ, JK)=ZWORKSHUGRAD_3D_8(JIJ, JK)*ZWORKSHUGRAD_3D_6(JIJ, JK)&
                                               * PINV_PDZZ(JIJ, JK)
          END DO
        END DO
        
        CALL DZF_PHY(D, ZWORKSHUGRAD_3D_11, ZWORKSHUGRAD_3D_12)
        
        DO JK=1, D%NKT
          DO JIJ=1, D%NIJT
            PRRS(JIJ, JK, 4) = PRRS(JIJ, JK, 4) +  2. *                                  &
            (- ZWORKSHUGRAD_3D_10(JIJ, JK)                        &
             + ZWORKSHUGRAD_3D_12(JIJ, JK)                       &
            )*PFRAC_ICE(JIJ, JK)
          END DO
        END DO
        
      ELSE
        
        DO JK=1, D%NKT
          DO JIJ=1, D%NIJT
            ZWORKSHUGRAD_3D_1(JIJ, JK)=ZFLX(JIJ, JK)*PINV_PDYY(JIJ, JK)
          END DO
        END DO
        
        CALL MZM_PHY(D, ZWORKSHUGRAD_3D_1, ZWORKSHUGRAD_3D_2)
        
        DO JK=1, D%NKT
          DO JIJ=1, D%NIJT
            ZWORKSHUGRAD_3D_3(JIJ, JK)=PRHODJ(JIJ, JK)*PAMOIST(JIJ, JK)*PSRCM(JIJ, JK)
          END DO
        END DO
        
        CALL MYM_PHY(D, ZWORKSHUGRAD_3D_3, ZWORKSHUGRAD_3D_4)
        
        DO JK=1, D%NKT
          DO JIJ=1, D%NIJT
            ZWORKSHUGRAD_3D_5(JIJ, JK)=PDZY(JIJ, JK)*(ZWORKSHUGRAD_3D_2(JIJ, JK))
          END DO
        END DO
        
        CALL MYF_PHY(D, ZWORKSHUGRAD_3D_5, ZWORKSHUGRAD_3D_6)
        
        DO JK=1, D%NKT
          DO JIJ=1, D%NIJT
            ZWORKSHUGRAD_3D_7(JIJ, JK)=PRHODJ(JIJ, JK)*PAMOIST(JIJ, JK)*PSRCM(JIJ, JK)
          END DO
        END DO
        
        CALL MZM_PHY(D, ZWORKSHUGRAD_3D_7, ZWORKSHUGRAD_3D_8)
        
        DO JK=1, D%NKT
          DO JIJ=1, D%NIJT
            ZWORKSHUGRAD_3D_9(JIJ, JK)=ZWORKSHUGRAD_3D_4(JIJ, JK)*ZFLX(JIJ, JK)/PDYY(JIJ, JK)
          END DO
        END DO
        
        CALL DYF_PHY(D, ZWORKSHUGRAD_3D_9, ZWORKSHUGRAD_3D_10)
        
        DO JK=1, D%NKT
          DO JIJ=1, D%NIJT
            ZWORKSHUGRAD_3D_11(JIJ, JK)=ZWORKSHUGRAD_3D_8(JIJ, JK)*ZWORKSHUGRAD_3D_6(JIJ, JK)&
                                               * PINV_PDZZ(JIJ, JK)
          END DO
        END DO
        
        CALL DZF_PHY(D, ZWORKSHUGRAD_3D_11, ZWORKSHUGRAD_3D_12)
        
        DO JK=1, D%NKT
          DO JIJ=1, D%NIJT
            PRRS(JIJ, JK, 2) = PRRS(JIJ, JK, 2) +  2. *                                  &
            (- ZWORKSHUGRAD_3D_10(JIJ, JK)                        &
             + ZWORKSHUGRAD_3D_12(JIJ, JK)                       &
            )
          END DO
        END DO
        
      END IF
    ELSE
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_1(JIJ, JK)=PRHODJ(JIJ, JK)*PAMOIST(JIJ, JK)*PSRCM(JIJ, JK)
        END DO
      END DO
      
      CALL MYM_PHY(D, ZWORKSHUGRAD_3D_1, ZWORKSHUGRAD_3D_2)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_3(JIJ, JK)=ZWORKSHUGRAD_3D_2(JIJ, JK)*ZFLX(JIJ, JK)
        END DO
      END DO
      
      CALL MXF_PHY(D, ZWORKSHUGRAD_3D_3, ZWORKSHUGRAD_3D_4)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZFLXC(JIJ, JK) = ZFLXC(JIJ, JK) + 2.*ZWORKSHUGRAD_3D_4(JIJ, JK)
        END DO
      END DO
      
      IF ( KRRI >= 1 ) THEN
        
        DO JK=1, D%NKT
          DO JIJ=1, D%NIJT
            ZWORKSHUGRAD_3D_1(JIJ, JK)=PRHODJ(JIJ, JK)*PAMOIST(JIJ, JK)*PSRCM(JIJ, JK)
          END DO
        END DO
        
        CALL MYM_PHY(D, ZWORKSHUGRAD_3D_1, ZWORKSHUGRAD_3D_2)
        
        DO JK=1, D%NKT
          DO JIJ=1, D%NIJT
            ZWORKSHUGRAD_3D_3(JIJ, JK)=ZWORKSHUGRAD_3D_2(JIJ, JK)*ZFLX(JIJ, JK)/PDYY(JIJ, JK)
          END DO
        END DO
        
        CALL DYF_PHY(D, ZWORKSHUGRAD_3D_3, ZWORKSHUGRAD_3D_4)
        
        DO JK=1, D%NKT
          DO JIJ=1, D%NIJT
            PRRS(JIJ, JK, 2) = PRRS(JIJ, JK, 2) - 2. *                                   &
            ZWORKSHUGRAD_3D_4(JIJ, JK)*(1.0-PFRAC_ICE(JIJ, JK))
          END DO
        END DO
        
        
        DO JK=1, D%NKT
          DO JIJ=1, D%NIJT
            ZWORKSHUGRAD_3D_1(JIJ, JK)=PRHODJ(JIJ, JK)*PAMOIST(JIJ, JK)*PSRCM(JIJ, JK)
          END DO
        END DO
        
        CALL MYM_PHY(D, ZWORKSHUGRAD_3D_1, ZWORKSHUGRAD_3D_2)
        
        DO JK=1, D%NKT
          DO JIJ=1, D%NIJT
            ZWORKSHUGRAD_3D_3(JIJ, JK)=ZWORKSHUGRAD_3D_2(JIJ, JK)*ZFLX(JIJ, JK)/PDYY(JIJ, JK)
          END DO
        END DO
        
        CALL DYF_PHY(D, ZWORKSHUGRAD_3D_3, ZWORKSHUGRAD_3D_4)
        
        DO JK=1, D%NKT
          DO JIJ=1, D%NIJT
            PRRS(JIJ, JK, 4) = PRRS(JIJ, JK, 4) - 2. *                                   &
            ZWORKSHUGRAD_3D_4(JIJ, JK)*PFRAC_ICE(JIJ, JK)
          END DO
        END DO
        
      ELSE
        
        DO JK=1, D%NKT
          DO JIJ=1, D%NIJT
            ZWORKSHUGRAD_3D_1(JIJ, JK)=PRHODJ(JIJ, JK)*PAMOIST(JIJ, JK)*PSRCM(JIJ, JK)
          END DO
        END DO
        
        CALL MYM_PHY(D, ZWORKSHUGRAD_3D_1, ZWORKSHUGRAD_3D_2)
        
        DO JK=1, D%NKT
          DO JIJ=1, D%NIJT
            ZWORKSHUGRAD_3D_3(JIJ, JK)=ZWORKSHUGRAD_3D_2(JIJ, JK)*ZFLX(JIJ, JK)/PDYY(JIJ, JK)
          END DO
        END DO
        
        CALL DYF_PHY(D, ZWORKSHUGRAD_3D_3, ZWORKSHUGRAD_3D_4)
        
        DO JK=1, D%NKT
          DO JIJ=1, D%NIJT
            PRRS(JIJ, JK, 2) = PRRS(JIJ, JK, 2) - 2. *                                   &
            ZWORKSHUGRAD_3D_4(JIJ, JK)
          END DO
        END DO
        
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

    CALL IO_FIELD_WRITE_PHY(D,TPFILE,TZFIELD,ZFLX)
  END IF
  !
  IF (KSPLT==1 .AND. TLES%LLES_CALL) THEN
    CALL SECOND_MNH(ZTIME1)
    CALL MYF_PHY(D, ZFLX, ZWORKSHUGRAD_3D_1)
    CALL LES_MEAN_SUBGRID_PHY(D, TLES,ZWORKSHUGRAD_3D_1(1:D%NIJT, 1:D%NKT), TLES%X_LES_SUBGRID_VRt ) 
    CALL GY_W_VW_PHY(D, OFLAT, PWM, PDYY, PDZZ, PDZY, ZWORKSHUGRAD_3D_1)
    CALL MZM_PHY(D, ZFLX, ZWORKSHUGRAD_3D_2)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWORKSHUGRAD_3D_3(JIJ, JK)=ZWORKSHUGRAD_3D_1(JIJ, JK)*ZWORKSHUGRAD_3D_2(JIJ, JK)
      END DO
    END DO
    
    CALL MYF_PHY(D, ZWORKSHUGRAD_3D_3, ZWORKSHUGRAD_3D_4)
    CALL MZF_PHY(D, ZWORKSHUGRAD_3D_4, ZWKLES)
    CALL LES_MEAN_SUBGRID_PHY(D, TLES,ZWKLES, TLES%X_LES_RES_ddxa_W_SBG_UaRt , .TRUE. )
    CALL GY_M_M_PHY(D, OFLAT, PTHLM, PDYY, PDZZ, PDZY, ZWORKSHUGRAD_3D_1)
    CALL MYF_PHY(D, ZFLX, ZWORKSHUGRAD_3D_2)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWKLES(JIJ, JK) = ZWORKSHUGRAD_3D_1(JIJ, JK)*ZWORKSHUGRAD_3D_2(JIJ, JK)
      END DO
    END DO
    
    CALL LES_MEAN_SUBGRID_PHY(D, TLES,ZWKLES, TLES%X_LES_RES_ddxa_Thl_SBG_UaRt , .TRUE. )
    CALL GY_M_M_PHY(D, OFLAT, PRM(:,:,1), PDYY, PDZZ, PDZY, ZWORKSHUGRAD_3D_1)
    CALL MYF_PHY(D, ZFLX, ZWORKSHUGRAD_3D_2)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWKLES(JIJ, JK) = ZWORKSHUGRAD_3D_1(JIJ, JK)*ZWORKSHUGRAD_3D_2(JIJ, JK)
      END DO
    END DO
    
    CALL LES_MEAN_SUBGRID_PHY(D, TLES,ZWKLES, TLES%X_LES_RES_ddxa_Rt_SBG_UaRt , .TRUE. )
    CALL SECOND_MNH(ZTIME2)
    TLES%XTIME_LES = TLES%XTIME_LES + ZTIME2 - ZTIME1
  END IF
!
  !
  IF (KRRL>0 .AND. KSPLT==1 .AND. TLES%LLES_CALL) THEN
    CALL SECOND_MNH(ZTIME1)
    CALL MYF_PHY(D, ZFLXC, ZWORKSHUGRAD_3D_1)
    CALL LES_MEAN_SUBGRID_PHY(D, TLES,ZWORKSHUGRAD_3D_1(1:D%NIJT, 1:D%NKT), TLES%X_LES_SUBGRID_VRc )
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
!!    ZVPTV(:,:) =                                                         &
!!        ZWORK(:,:)*MYM(ETHETA(KRR,KRRI,PTHLT,PEXNREF,PRT,PLOCPT,PSRCM)) +       &
!!         ZFLX(:,:)*MYM(EMOIST(KRR,KRRI,PTHLT,PEXNREF,PRT,PLOCPT,PSRCM))
!!  ELSE
!!    ZVPTV(:,:) = 0.
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
!!    CALL IO_FIELD_WRITE_PHY(D,TPFILE,TZFIELD,ZVPTV)
!!  END IF
!!!
!!ELSE
!!  ZVPTV(:,:)=ZWORK(:,:)
!!END IF
!
!
END SUBROUTINE TURB_HOR_THERMO_FLUX
END MODULE MODE_TURB_HOR_THERMO_FLUX
