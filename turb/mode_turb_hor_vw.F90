!MNH_LIC Copyright 1994-2025 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_TURB_HOR_VW
IMPLICIT NONE
CONTAINS
      SUBROUTINE TURB_HOR_VW(D,TURBN,TLES,KSPLT,                     &
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
!!                                   + remove the use of W=0 at the ground
!!                                   + extrapolation under the ground
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
USE MODE_SHUMAN_PHY, ONLY: DYF_PHY, DYM2D_PHY, DZF_PHY, DZM_PHY, MXF_PHY, MYF2D_PHY, MYF_PHY, MYM2D_PHY, MYM_PHY, MZF_PHY, MZM_PHY
USE MODE_GRADIENT_M_PHY, ONLY: GY_M_V_PHY
USE MODE_GRADIENT_W_PHY, ONLY: GY_W_VW_PHY
USE MODE_GRADIENT_V_PHY, ONLY: GZ_V_VW_PHY
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
INTEGER,                  INTENT(IN)    ::  KRR          ! number of moist var.
INTEGER,                  INTENT(IN)    ::  KSV          ! number of sv var.
LOGICAL,                  INTENT(IN)    ::  OFLAT        ! Logical for zero ororography
LOGICAL,                  INTENT(IN)    ::  O2D          ! Logical for 2D model version (modd_conf)
TYPE(TFILEDATA),          INTENT(INOUT) ::  TPFILE       ! Output file
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    ::  PK          ! Turbulent diffusion doef.
                                                        ! PK = PLM * SQRT(PTKEM)
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    ::  PINV_PDYY   ! 1./PDYY
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    ::  PINV_PDZZ   ! 1./PDZZ
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    ::  PMZM_PRHODJ ! MZM(PRHODJ)
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    ::  PDYY, PDZZ, PDZY 
                                                         ! Metric coefficients
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    ::  PRHODJ       ! density * grid volume
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    ::  PTHVREF      ! ref. state VPT       
!
! Variables at t-1
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    ::  PVM,PWM,PTHLM
REAL, DIMENSION(D%NIJT,D%NKT,KRR), INTENT(IN)    ::  PRM
REAL, DIMENSION(D%NIJT,D%NKT,KSV), INTENT(IN)    ::  PSVM
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    ::  PTKEM        ! TKE at time t- dt
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    ::  PLM          ! Turb. mixing length
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) ::  PRVS, PRWS   ! var. at t+1 -split-
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) ::  PDP          ! TKE production terms
!
!
!
!*       0.2  declaration of local variables
!
REAL, DIMENSION(D%NIJT,D%NKT)  :: ZFLX,ZWORK
    ! work arrays
!   
INTEGER             :: IKB,IKE,IKU, IIJT, IKT
                                    ! Index values for the Beginning and End
                                    ! mass points of the domain  
INTEGER             :: JSV,JIJ,JK          ! scalar loop counter
!
REAL, DIMENSION(D%NIJT,D%NKT)  :: GY_W_VW_PWM
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
IKE = SIZE(PWM,2)-JPVEXT    
IKU = SIZE(PWM,2)
IIJT=D%NIJT
IKT=D%NKT 
!
!
IF (.NOT. O2D) CALL GY_W_VW_PHY(D, OFLAT, PWM, PDYY, PDZZ, PDZY, GY_W_VW_PWM)
!
!
!*      14.   < V'W'>
!             -------
!
! residual part of < V'W'> depending on dw/dy
!
IF (.NOT. O2D) THEN
  CALL MZM_PHY(D, PK, ZWORKSHUGRAD_3D_1)
  CALL MYM_PHY(D, ZWORKSHUGRAD_3D_1, ZWORKSHUGRAD_3D_2)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      ZFLX(JIJ, JK) =                                                      &
        - XCMFS * ZWORKSHUGRAD_3D_2(JIJ, JK) * GY_W_VW_PWM(JIJ, JK)
    END DO
  END DO
  
ELSE

DO JK=1, IKT
    DO JIJ=1, IIJT
      ZFLX(JIJ, JK) = 0.
  END DO
  END DO

END IF
!

DO JIJ=1, IIJT
  ZFLX(JIJ, IKE+1) = 0.  ! rigid wall condition => no turbulent flux
  !
  !
  ! Nullify the flux at the ground level because it has been fully taken into
  ! account in turb_ver and extrapolate the flux under the ground 
  ZFLX(JIJ, IKB) = 0.
  ZFLX(JIJ, IKB-1)= 2.*ZFLX(JIJ, IKB) - ZFLX(JIJ, IKB+1)
END DO

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

  CALL IO_FIELD_WRITE_PHY(D,TPFILE,TZFIELD,ZFLX)
END IF
!
! compute the source for rho*V due to this residual flux ( the other part is
! taken into account in TURB_VER)
IF (.NOT. O2D) THEN
  CALL MYM_PHY(D, PMZM_PRHODJ, ZWORKSHUGRAD_3D_1)
  CALL MYM_PHY(D, PDZZ, ZWORKSHUGRAD_3D_2)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      ZWORKSHUGRAD_3D_3(JIJ, JK)=ZFLX(JIJ, JK)* ZWORKSHUGRAD_3D_1(JIJ, JK) / ZWORKSHUGRAD_3D_2(JIJ, JK)
    END DO
  END DO
  
  CALL DZF_PHY(D, ZWORKSHUGRAD_3D_3, ZWORKSHUGRAD_3D_4)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      PRVS(JIJ, JK) = PRVS(JIJ, JK) - ZWORKSHUGRAD_3D_4(JIJ, JK)
    END DO
  END DO
  
END IF
!
!computation of the source for rho*W due to this flux
IF (.NOT. O2D) THEN 
  IF (.NOT. OFLAT) THEN
    CALL MYM_PHY(D, PRHODJ, ZWORKSHUGRAD_3D_1)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWORKSHUGRAD_3D_2(JIJ, JK)=ZFLX(JIJ, JK)*PDZY(JIJ, JK)
      END DO
    END DO
    
    CALL MZF_PHY(D, ZWORKSHUGRAD_3D_2, ZWORKSHUGRAD_3D_3)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWORKSHUGRAD_3D_4(JIJ, JK)=ZWORKSHUGRAD_3D_1(JIJ, JK) * PINV_PDYY(JIJ, JK)
      END DO
    END DO
    
    CALL MZM_PHY(D, ZWORKSHUGRAD_3D_4, ZWORKSHUGRAD_3D_5)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWORKSHUGRAD_3D_6(JIJ, JK)=ZWORKSHUGRAD_3D_3(JIJ, JK) * PINV_PDYY(JIJ, JK)
      END DO
    END DO
    
    CALL MYF_PHY(D, ZWORKSHUGRAD_3D_6, ZWORKSHUGRAD_3D_7)
    CALL MZF_PHY(D, PDZZ, ZWORKSHUGRAD_3D_8)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWORKSHUGRAD_3D_9(JIJ, JK)=ZWORKSHUGRAD_3D_5(JIJ, JK) * ZFLX(JIJ, JK)
      END DO
    END DO
    
    CALL DYF_PHY(D, ZWORKSHUGRAD_3D_9, ZWORKSHUGRAD_3D_10)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWORKSHUGRAD_3D_11(JIJ, JK)=PRHODJ(JIJ, JK) * ZWORKSHUGRAD_3D_7(JIJ, JK) / ZWORKSHUGRAD_3D_8(JIJ, JK)
      END DO
    END DO
    
    CALL DZM_PHY(D, ZWORKSHUGRAD_3D_11, ZWORKSHUGRAD_3D_12)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        PRWS(JIJ, JK) = PRWS(JIJ, JK)                              &
              -ZWORKSHUGRAD_3D_10(JIJ, JK)           &
              +ZWORKSHUGRAD_3D_12(JIJ, JK)
      END DO
    END DO
    
  ELSE
    CALL MYM_PHY(D, PRHODJ, ZWORKSHUGRAD_3D_1)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWORKSHUGRAD_3D_2(JIJ, JK)=ZWORKSHUGRAD_3D_1(JIJ, JK) * PINV_PDYY(JIJ, JK)
      END DO
    END DO
    
    CALL MZM_PHY(D, ZWORKSHUGRAD_3D_2, ZWORKSHUGRAD_3D_3)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWORKSHUGRAD_3D_4(JIJ, JK)=ZWORKSHUGRAD_3D_3(JIJ, JK) * ZFLX(JIJ, JK)
      END DO
    END DO
    
    CALL DYF_PHY(D, ZWORKSHUGRAD_3D_4, ZWORKSHUGRAD_3D_5)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        PRWS(JIJ, JK) = PRWS(JIJ, JK) - ZWORKSHUGRAD_3D_5(JIJ, JK)
      END DO
    END DO
    
  END IF
END IF
!
IF (KSPLT==1) THEN
  ! 
  !Contribution to the dynamic production of TKE:
  !
  IF (.NOT. O2D) THEN
    CALL GZ_V_VW_PHY(D, PVM, PDZZ, ZWORKSHUGRAD_3D_1)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWORKSHUGRAD_3D_2(JIJ, JK)=ZFLX(JIJ, JK) *( ZWORKSHUGRAD_3D_1(JIJ, JK) + GY_W_VW_PWM(JIJ, JK) )
      END DO
    END DO
    
    CALL MYF_PHY(D, ZWORKSHUGRAD_3D_2, ZWORKSHUGRAD_3D_3)
    CALL MZF_PHY(D, ZWORKSHUGRAD_3D_3, ZWORKSHUGRAD_3D_4)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWORK(JIJ, JK) =-ZWORKSHUGRAD_3D_4(JIJ, JK)
      END DO
    END DO
    
  !
  !
  ! evaluate the dynamic production at w(IKB+1) in PDP(IKB)
    
    DO JIJ=1, D%NIJT
      ZWORKSHUGRAD_2D_1(JIJ)=(PWM(JIJ, IKB+2)-PWM(JIJ, IKB+1))                   &
                    /(PDZZ(JIJ, IKB+2)+PDZZ(JIJ, IKB+1))               &
                  +(PWM(JIJ, IKB+1)-PWM(JIJ, IKB))                   &
                    /(PDZZ(JIJ, IKB+1)+PDZZ(JIJ, IKB))
    END DO
    
    CALL MYM2D_PHY(D, ZWORKSHUGRAD_2D_1, ZWORKSHUGRAD_2D_2)
    CALL DYM2D_PHY(D, PWM(:,IKB+1), ZWORKSHUGRAD_2D_3)
    CALL MYM2D_PHY(D, PDZZ(:,IKB+1), ZWORKSHUGRAD_2D_4)
    
    DO JIJ=1, D%NIJT
      ZWORKSHUGRAD_2D_5(JIJ)=ZFLX(JIJ, IKB+1) *                                                 &
       (   (PVM(JIJ, IKB+1)-PVM(JIJ, IKB)) / ZWORKSHUGRAD_2D_4(JIJ)  &
         + ( ZWORKSHUGRAD_2D_3(JIJ)                                         &
            -ZWORKSHUGRAD_2D_2(JIJ) * PDZY(JIJ, IKB+1)                                        &
           ) / (0.5*(PDYY(JIJ, IKB+1)+PDYY(JIJ, IKB)))                 &
       )
    END DO
    
    CALL MYF2D_PHY(D, ZWORKSHUGRAD_2D_5, ZWORKSHUGRAD_2D_6)
    
    DO JIJ=1, D%NIJT
      ZWORK(JIJ, IKB) = - ZWORKSHUGRAD_2D_6(JIJ)
    END DO
    
  ENDIF
  !
  ! dynamic production computation
  IF (.NOT. O2D) THEN
    
    DO JK=1, IKT
      DO JIJ=1, IIJT
        PDP(JIJ, JK) = PDP(JIJ, JK) + ZWORK(JIJ, JK)  
      END DO
    END DO
    
  ENDIF
  !
END IF
!
! Storage in the LES configuration (addition to TURB_VER computation)
!
IF (TLES%LLES_CALL .AND. KSPLT==1) THEN
  CALL SECOND_MNH(ZTIME1)
  CALL MYF_PHY(D, ZFLX, ZWORKSHUGRAD_3D_1)
  CALL MZF_PHY(D, ZWORKSHUGRAD_3D_1, ZWORKSHUGRAD_3D_2)
  CALL LES_MEAN_SUBGRID_PHY(D, TLES, ZWORKSHUGRAD_3D_2(1:D%NIJT, 1:D%NKT), TLES%X_LES_SUBGRID_WV , .TRUE. )
  CALL GZ_V_VW_PHY(D, PVM, PDZZ, ZWORKSHUGRAD_3D_1)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      ZWORKSHUGRAD_3D_2(JIJ, JK)=ZWORKSHUGRAD_3D_1(JIJ, JK)*ZFLX(JIJ, JK)
    END DO
  END DO
  
  CALL MYF_PHY(D, ZWORKSHUGRAD_3D_2, ZWORKSHUGRAD_3D_3)
  CALL MZF_PHY(D, ZWORKSHUGRAD_3D_3, ZWORKSHUGRAD_3D_4)
  CALL LES_MEAN_SUBGRID_PHY(D, TLES, ZWORKSHUGRAD_3D_4(1:D%NIJT, 1:D%NKT),&
                         TLES%X_LES_RES_ddxa_V_SBG_UaV , .TRUE.)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      ZWORKSHUGRAD_3D_1(JIJ, JK)=GY_W_VW_PWM(JIJ, JK)*ZFLX(JIJ, JK)
    END DO
  END DO
  
  CALL MYF_PHY(D, ZWORKSHUGRAD_3D_1, ZWORKSHUGRAD_3D_2)
  CALL MZF_PHY(D, ZWORKSHUGRAD_3D_2, ZWORKSHUGRAD_3D_3)
  CALL LES_MEAN_SUBGRID_PHY(D, TLES, ZWORKSHUGRAD_3D_3(1:D%NIJT, 1:D%NKT),&
                         TLES%X_LES_RES_ddxa_W_SBG_UaW , .TRUE.)
  CALL GY_M_V_PHY(D, OFLAT, PTHLM, PDYY, PDZZ, PDZY, ZWORKSHUGRAD_3D_1)
  CALL MZF_PHY(D, ZFLX, ZWORKSHUGRAD_3D_2)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      ZWORKSHUGRAD_3D_3(JIJ, JK)=ZWORKSHUGRAD_3D_1(JIJ, JK)*ZWORKSHUGRAD_3D_2(JIJ, JK)
    END DO
  END DO
  
  CALL MXF_PHY(D, ZWORKSHUGRAD_3D_3, ZWORKSHUGRAD_3D_4)
  CALL LES_MEAN_SUBGRID_PHY(D, TLES, ZWORKSHUGRAD_3D_4(1:D%NIJT, 1:D%NKT),&
                         TLES%X_LES_RES_ddxa_Thl_SBG_UaW , .TRUE.)
  IF (KRR>=1) THEN
    CALL GY_M_V_PHY(D, OFLAT, PRM(:,:,1), PDYY, PDZZ, PDZY, ZWORKSHUGRAD_3D_1)
    CALL MZF_PHY(D, ZFLX, ZWORKSHUGRAD_3D_2)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWORKSHUGRAD_3D_3(JIJ, JK)=ZWORKSHUGRAD_3D_1(JIJ, JK)*ZWORKSHUGRAD_3D_2(JIJ, JK)
      END DO
    END DO
    
    CALL MXF_PHY(D, ZWORKSHUGRAD_3D_3, ZWORKSHUGRAD_3D_4)
    CALL LES_MEAN_SUBGRID_PHY(D, TLES, ZWORKSHUGRAD_3D_4(1:D%NIJT, 1:D%NKT), &
                           TLES%X_LES_RES_ddxa_Rt_SBG_UaW , .TRUE.)
  END IF
  DO JSV=1,KSV
    CALL GY_M_V_PHY(D, OFLAT, PSVM(:,:,JSV), PDYY, PDZZ, PDZY, ZWORKSHUGRAD_3D_1)
    CALL MZF_PHY(D, ZFLX, ZWORKSHUGRAD_3D_2)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWORKSHUGRAD_3D_3(JIJ, JK)=ZWORKSHUGRAD_3D_1(JIJ, JK)*ZWORKSHUGRAD_3D_2(JIJ, JK)
      END DO
    END DO
    
    CALL MXF_PHY(D, ZWORKSHUGRAD_3D_3, ZWORKSHUGRAD_3D_4)
    CALL LES_MEAN_SUBGRID_PHY(D, TLES, ZWORKSHUGRAD_3D_4(1:D%NIJT, 1:D%NKT), &
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
