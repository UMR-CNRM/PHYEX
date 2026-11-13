!MNH_LIC Copyright 1994-2025 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_TURB_HOR_TKE
IMPLICIT NONE
CONTAINS
      SUBROUTINE TURB_HOR_TKE(D,KSPLT,TLES,OFLAT,O2D,                &
                      PDXX, PDYY, PDZZ,PDZX,PDZY,                    &
                      PINV_PDXX, PINV_PDYY, PINV_PDZZ, PMZM_PRHODJ,  &
                      PK, PRHODJ, PTKEM,                             &
                      PTRH                                           )
!     ################################################################
!
!
!!****  *TURB_HOR_TKE* computes the horizontal turbulant transports of Tke
!!
!!    PURPOSE
!!    -------

!!**  METHOD
!!    ------
!!
!!
!!
!!    EXTERNAL
!!    --------
!!       
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!           
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      Joan Cuxart             * INM and Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       Aug 29, 1994
!!                     Mar 07  2001 (V. Masson and J. Stein) new routine
!!                     Nov 06, 2002 (V. Masson) LES budgets
!!                     04/2016 (M.Moge) Use openACC directives to port the TURB part of Meso-NH on GPU
!! --------------------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
!
USE MODD_DIMPHYEX,       ONLY: DIMPHYEX_t
USE MODD_CTURB,          ONLY: CTURB_XCET => XCET
USE MODD_PARAMETERS, ONLY: JPVEXT
USE MODD_LES, ONLY: TLES_t
!
!
USE MODI_LES_MEAN_SUBGRID_PHY, ONLY: LES_MEAN_SUBGRID_PHY
!
USE MODI_SECOND_MNH, ONLY: SECOND_MNH
USE MODE_SHUMAN_PHY, ONLY: DXF_PHY, DXM2D_PHY, DYF_PHY, DYM2D_PHY, DZF_PHY, MXF_PHY, MXM2D_PHY, MXM_PHY, MYF_PHY, MYM2D_PHY, &
  & MYM_PHY, MZM_PHY
USE MODE_GRADIENT_M_PHY, ONLY: GX_M_U_PHY, GY_M_V_PHY
!
IMPLICIT NONE
!
!
!*       0.1  declaration of arguments
!
!
TYPE(DIMPHYEX_t),         INTENT(IN)    :: D
TYPE(TLES_t),             INTENT(INOUT) :: TLES          ! modd_les structure
INTEGER,                  INTENT(IN) :: KSPLT        ! current split index
LOGICAL,                  INTENT(IN) ::  OFLAT       ! Logical for zero ororography
LOGICAL,                  INTENT(IN) ::  O2D         ! Logical for 2D model version (modd_conf)
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN) :: PDXX, PDYY, PDZZ, PDZX, PDZY 
                                                     ! Metric coefficients
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN) :: PK           ! Turbulent diffusion doef.
                                                     ! PK = PLM * SQRT(PTKEM)
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN) :: PINV_PDXX    ! 1./PDXX
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN) :: PINV_PDYY    ! 1./PDYY
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN) :: PINV_PDZZ    ! 1./PDZZ
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN) :: PMZM_PRHODJ  ! MZM(PRHODJ)
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN) :: PRHODJ       ! density * grid volume
!
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    ::  PTKEM    ! TKE at time t- dt
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(OUT)   ::  PTRH     ! horizontal transport of Tke
!
!
!
!*       0.2  declaration of local variables
!
INTEGER :: IKB, IKU, IIJT, IKT, JIJ
!
REAL, DIMENSION(D%NIJT,1+JPVEXT:3+JPVEXT) :: ZCOEFF 
                                    ! coefficients for the uncentred gradient 
                                    ! computation near the ground
!
REAL, DIMENSION(D%NIJT,D%NKT):: ZFLX
!
REAL :: ZTIME1, ZTIME2
!
INTEGER :: JK
REAL, DIMENSION(D%NIJT) :: ZWORKSHUGRAD_2D_1
REAL, DIMENSION(D%NIJT) :: ZWORKSHUGRAD_2D_2
REAL, DIMENSION(D%NIJT) :: ZWORKSHUGRAD_2D_3
REAL, DIMENSION(D%NIJT, D%NKT) :: ZWORKSHUGRAD_3D_1
REAL, DIMENSION(D%NIJT, D%NKT) :: ZWORKSHUGRAD_3D_2
REAL, DIMENSION(D%NIJT, D%NKT) :: ZWORKSHUGRAD_3D_3
REAL, DIMENSION(D%NIJT, D%NKT) :: ZWORKSHUGRAD_3D_4
REAL, DIMENSION(D%NIJT, D%NKT) :: ZWORKSHUGRAD_3D_5
REAL, DIMENSION(D%NIJT, D%NKT) :: ZWORKSHUGRAD_3D_6
REAL, DIMENSION(D%NIJT, D%NKT) :: ZWORKSHUGRAD_3D_7
REAL, DIMENSION(D%NIJT, D%NKT) :: ZWORKSHUGRAD_3D_8
REAL, DIMENSION(D%NIJT, D%NKT) :: ZWORKSHUGRAD_3D_9
REAL :: XCET
! ---------------------------------------------------------------------------
!
!*       1.   PRELIMINARY COMPUTATIONS
!             ------------------------
!
XCET = CTURB_XCET
!
IKB = 1.+JPVEXT
IKU = SIZE(PTKEM,2)
IIJT=D%NIJT
IKT=D%NKT 
!
!  compute the coefficients for the uncentred gradient computation near the 
!  ground
!

DO JIJ=1, IIJT
  ZCOEFF(JIJ, IKB+2)= - PDZZ(JIJ, IKB+1) /      &
       ( (PDZZ(JIJ, IKB+2)+PDZZ(JIJ, IKB+1)) * PDZZ(JIJ, IKB+2) )
  ZCOEFF(JIJ, IKB+1)=   (PDZZ(JIJ, IKB+2)+PDZZ(JIJ, IKB+1)) /      &
       ( PDZZ(JIJ, IKB+1) * PDZZ(JIJ, IKB+2) )
  ZCOEFF(JIJ, IKB)= - (PDZZ(JIJ, IKB+2)+2.*PDZZ(JIJ, IKB+1)) /      &
       ( (PDZZ(JIJ, IKB+2)+PDZZ(JIJ, IKB+1)) * PDZZ(JIJ, IKB+1) )
END DO

!
!--------------------------------------------------------------------
!
!*       2.   horizontal transport of Tke u'e
!             -------------------------------
!
!
CALL MXM_PHY(D, PK, ZWORKSHUGRAD_3D_1)
CALL GX_M_U_PHY(D, OFLAT, PTKEM, PDXX, PDZZ, PDZX, ZWORKSHUGRAD_3D_2)

DO JK=1, D%NKT
  DO JIJ=1, D%NIJT
    ZFLX(JIJ, JK) = -XCET * ZWORKSHUGRAD_3D_1(JIJ, JK) * ZWORKSHUGRAD_3D_2(JIJ, JK)
  END DO
END DO

 ! < u'e >
!
! special case near the ground ( uncentred gradient )
!

DO JIJ=1, IIJT
  ZFLX(JIJ, IKB) =  ZCOEFF(JIJ, IKB+2)*PTKEM(JIJ, IKB+2)                     &
                 + ZCOEFF(JIJ, IKB+1)*PTKEM(JIJ, IKB+1)                     &
                 + ZCOEFF(JIJ, IKB)*PTKEM(JIJ, IKB)     
END DO

!
CALL DXM2D_PHY(D, PTKEM(:,IKB), ZWORKSHUGRAD_2D_1)
CALL MXM2D_PHY(D, ZFLX (:,IKB), ZWORKSHUGRAD_2D_2)
CALL MXM2D_PHY(D, PK(:,IKB), ZWORKSHUGRAD_2D_3)

DO JIJ=1, D%NIJT
  ZFLX(JIJ, IKB) =                                                      &
     - XCET * ZWORKSHUGRAD_2D_3(JIJ)                           *  (        &
         ZWORKSHUGRAD_2D_1(JIJ) * PINV_PDXX(JIJ, IKB)               &
        -ZWORKSHUGRAD_2D_2(JIJ) * PINV_PDXX(JIJ, IKB)               &
         * 0.5 * ( PDZX(JIJ, IKB+1) + PDZX(JIJ, IKB) )     )
END DO

!
! extrapolate the fluxes to obtain < u'e > = 0 at the ground
!

DO JIJ=1, IIJT
  ZFLX(JIJ, IKB-1) = - ZFLX(JIJ, IKB)
  !
  ! let the same flux at IKU-1 and IKU level
  !
  ZFLX(JIJ, IKU) =  ZFLX(JIJ, IKU-1)
END DO

!
IF (.NOT. OFLAT) THEN
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      ZWORKSHUGRAD_3D_1(JIJ, JK)=ZFLX(JIJ, JK)*PINV_PDXX(JIJ, JK)
    END DO
  END DO
  
  CALL MZM_PHY(D, ZWORKSHUGRAD_3D_1, ZWORKSHUGRAD_3D_2)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      ZWORKSHUGRAD_3D_3(JIJ, JK)=PDZX(JIJ, JK) * ZWORKSHUGRAD_3D_2(JIJ, JK)
    END DO
  END DO
  
  CALL MXF_PHY(D, ZWORKSHUGRAD_3D_3, ZWORKSHUGRAD_3D_4)
  CALL MXM_PHY(D, PRHODJ, ZWORKSHUGRAD_3D_5)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      ZWORKSHUGRAD_3D_6(JIJ, JK)=ZWORKSHUGRAD_3D_5(JIJ, JK) * ZFLX(JIJ, JK)                             * PINV_PDXX(JIJ, JK)
    END DO
  END DO
  
  CALL DXF_PHY(D, ZWORKSHUGRAD_3D_6, ZWORKSHUGRAD_3D_7)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      ZWORKSHUGRAD_3D_8(JIJ, JK)=PMZM_PRHODJ(JIJ, JK) * ZWORKSHUGRAD_3D_4(JIJ, JK) * PINV_PDZZ(JIJ, JK)
    END DO
  END DO
  
  CALL DZF_PHY(D, ZWORKSHUGRAD_3D_8, ZWORKSHUGRAD_3D_9)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      PTRH(JIJ, JK) =-(  ZWORKSHUGRAD_3D_7(JIJ, JK)&
              - ZWORKSHUGRAD_3D_9(JIJ, JK)&
             ) /PRHODJ(JIJ, JK)
    END DO
  END DO
  
ELSE
  CALL MXM_PHY(D, PRHODJ, ZWORKSHUGRAD_3D_1)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      ZWORKSHUGRAD_3D_2(JIJ, JK)=ZWORKSHUGRAD_3D_1(JIJ, JK) * ZFLX(JIJ, JK)                             * PINV_PDXX(JIJ, JK)
    END DO
  END DO
  
  CALL DXF_PHY(D, ZWORKSHUGRAD_3D_2, ZWORKSHUGRAD_3D_3)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      PTRH(JIJ, JK) =-(  ZWORKSHUGRAD_3D_3(JIJ, JK)&
             ) /PRHODJ(JIJ, JK)
    END DO
  END DO
  
END IF
!
IF (TLES%LLES_CALL .AND. KSPLT==1) THEN
  CALL SECOND_MNH(ZTIME1)
  CALL MXF_PHY(D, ZFLX, ZWORKSHUGRAD_3D_1)
  CALL LES_MEAN_SUBGRID_PHY(D, TLES, ZWORKSHUGRAD_3D_1(1:D%NIJT, 1:D%NKT), TLES%X_LES_SUBGRID_UTke ) 
  CALL SECOND_MNH(ZTIME2)
  TLES%XTIME_LES = TLES%XTIME_LES + ZTIME2 - ZTIME1
END IF
!
!
!--------------------------------------------------------------------
!
!*       3.   horizontal transport of Tke v'e
!             -------------------------------
!
IF (.NOT. O2D) THEN
  CALL MYM_PHY(D, PK, ZWORKSHUGRAD_3D_1)
  CALL GY_M_V_PHY(D, OFLAT, PTKEM, PDYY, PDZZ, PDZY, ZWORKSHUGRAD_3D_2)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      ZFLX(JIJ, JK) =-XCET * ZWORKSHUGRAD_3D_1(JIJ, JK) * ZWORKSHUGRAD_3D_2(JIJ, JK)
    END DO
  END DO
  
 ! < v'e >
!
! special case near the ground ( uncentred gradient )
!

DO JIJ=1, IIJT
    ZFLX(JIJ, IKB) =  ZCOEFF(JIJ, IKB+2)*PTKEM(JIJ, IKB+2)                     &
                   + ZCOEFF(JIJ, IKB+1)*PTKEM(JIJ, IKB+1)                     &
                   + ZCOEFF(JIJ, IKB)*PTKEM(JIJ, IKB)     
  !
END DO

  CALL DYM2D_PHY(D, PTKEM(:,IKB), ZWORKSHUGRAD_2D_1)
  CALL MYM2D_PHY(D, ZFLX (:,IKB), ZWORKSHUGRAD_2D_2)
  CALL MYM2D_PHY(D, PK(:,IKB), ZWORKSHUGRAD_2D_3)
  
  DO JIJ=1, D%NIJT
    ZFLX(JIJ, IKB) =                                                      &
       - XCET * ZWORKSHUGRAD_2D_3(JIJ)                        *  (           &
         ZWORKSHUGRAD_2D_1(JIJ) * PINV_PDYY(JIJ, IKB)                 &
       - ZWORKSHUGRAD_2D_2(JIJ) * PINV_PDYY(JIJ, IKB)                 &
           * 0.5 * ( PDZY(JIJ, IKB+1) + PDZY(JIJ, IKB) )  )
  END DO
  
!
!    extrapolate the fluxes to obtain < v'e > = 0 at the ground
!

DO JIJ=1, IIJT
    ZFLX(JIJ, IKB-1) = - ZFLX(JIJ, IKB)
  !
  !   let the same flux at IKU-1 and IKU level
  !
    ZFLX(JIJ, IKU) =  ZFLX(JIJ, IKU-1)
END DO

!
! complete the explicit turbulent transport
!
  IF (.NOT. OFLAT) THEN
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWORKSHUGRAD_3D_1(JIJ, JK)=ZFLX(JIJ, JK)*PINV_PDYY(JIJ, JK)
      END DO
    END DO
    
    CALL MZM_PHY(D, ZWORKSHUGRAD_3D_1, ZWORKSHUGRAD_3D_2)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWORKSHUGRAD_3D_3(JIJ, JK)=PDZY(JIJ, JK) * ZWORKSHUGRAD_3D_2(JIJ, JK)
      END DO
    END DO
    
    CALL MYF_PHY(D, ZWORKSHUGRAD_3D_3, ZWORKSHUGRAD_3D_4)
    CALL MYM_PHY(D, PRHODJ, ZWORKSHUGRAD_3D_5)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWORKSHUGRAD_3D_6(JIJ, JK)=ZWORKSHUGRAD_3D_5(JIJ, JK) * ZFLX(JIJ, JK)                              * PINV_PDYY(JIJ, JK)
      END DO
    END DO
    
    CALL DYF_PHY(D, ZWORKSHUGRAD_3D_6, ZWORKSHUGRAD_3D_7)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWORKSHUGRAD_3D_8(JIJ, JK)=PMZM_PRHODJ(JIJ, JK) * ZWORKSHUGRAD_3D_4(JIJ, JK) * PINV_PDZZ(JIJ, JK)
      END DO
    END DO
    
    CALL DZF_PHY(D, ZWORKSHUGRAD_3D_8, ZWORKSHUGRAD_3D_9)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        PTRH(JIJ, JK) = PTRH(JIJ, JK) - (  ZWORKSHUGRAD_3D_7(JIJ, JK)  &
                       - ZWORKSHUGRAD_3D_9(JIJ, JK)  &
                      ) /PRHODJ(JIJ, JK)
      END DO
    END DO
    
  ELSE
    CALL MYM_PHY(D, PRHODJ, ZWORKSHUGRAD_3D_1)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWORKSHUGRAD_3D_2(JIJ, JK)=ZWORKSHUGRAD_3D_1(JIJ, JK) * ZFLX(JIJ, JK)                              * PINV_PDYY(JIJ, JK)
      END DO
    END DO
    
    CALL DYF_PHY(D, ZWORKSHUGRAD_3D_2, ZWORKSHUGRAD_3D_3)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        PTRH(JIJ, JK) = PTRH(JIJ, JK) - (  ZWORKSHUGRAD_3D_3(JIJ, JK)  &
                      ) /PRHODJ(JIJ, JK)
      END DO
    END DO
    
  END IF
!
  IF (TLES%LLES_CALL .AND. KSPLT==1) THEN
    CALL SECOND_MNH(ZTIME1)
    CALL MYF_PHY(D, ZFLX, ZWORKSHUGRAD_3D_1)
    CALL LES_MEAN_SUBGRID_PHY(D, TLES, ZWORKSHUGRAD_3D_1(1:D%NIJT, 1:D%NKT), TLES%X_LES_SUBGRID_VTke )
    CALL SECOND_MNH(ZTIME2)
    TLES%XTIME_LES = TLES%XTIME_LES + ZTIME2 - ZTIME1
  END IF
!
END IF
!
!----------------------------------------------------------------------------
!
END SUBROUTINE TURB_HOR_TKE
END MODULE MODE_TURB_HOR_TKE
