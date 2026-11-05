!MNH_LIC Copyright 1994-2025 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_TURB_HOR_UW
IMPLICIT NONE
CONTAINS
!     ################################################################
      SUBROUTINE TURB_HOR_UW(D,TURBN,TLES,KSPLT,                       &
                      KRR,KSV,OFLAT,                                 &
                      TPFILE,                                        &
                      PK,PINV_PDXX,PINV_PDZZ,PMZM_PRHODJ,            &
                      PDXX,PDZZ,PDZX,                                &
                      PRHODJ,PTHVREF,                                &
                      PUM,PWM,PTHLM,PRM,PSVM,                        &
                      PTKEM,PLM,                                     &
                      PDP,                                           &
                      PRUS,PRWS                                      )
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
USE MODE_SHUMAN_PHY, ONLY: DXF_PHY, DXM2D_PHY, DZF_PHY, DZM_PHY, MXF2D_PHY, MXF_PHY, MXM2D_PHY, MXM_PHY, MZF_PHY, MZM_PHY
USE MODE_GRADIENT_M_PHY, ONLY: GX_M_U_PHY
USE MODE_GRADIENT_W_PHY, ONLY: GX_W_UW_PHY
USE MODE_GRADIENT_U_PHY, ONLY: GZ_U_UW_PHY
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
TYPE(TFILEDATA),          INTENT(INOUT) ::  TPFILE       ! Output file
!

REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    ::  PK          ! Turbulent diffusion doef.
                                                        ! PK = PLM * SQRT(PTKEM)
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    ::  PINV_PDXX   ! 1./PDXX
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    ::  PINV_PDZZ   ! 1./PDZZ
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    ::  PMZM_PRHODJ ! MZM(PRHODJ)
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    ::  PDXX, PDZZ, PDZX
                                                         ! Metric coefficients
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    ::  PRHODJ       ! density * grid volume
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    ::  PTHVREF      ! ref. state VPT       
!
! Variables at t-1
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    ::  PUM,PWM,PTHLM
REAL, DIMENSION(D%NIJT,D%NKT,KRR), INTENT(IN)    ::  PRM
REAL, DIMENSION(D%NIJT,D%NKT,KSV), INTENT(IN)    ::  PSVM
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    ::  PTKEM        ! TKE at time t- dt
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    ::  PLM          ! Turb. mixing length
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) ::  PRUS, PRWS
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) ::  PDP          ! TKE production terms
!
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
REAL, DIMENSION(D%NIJT,D%NKT)  :: GX_W_UW_PWM
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
CALL GX_W_UW_PHY(D, OFLAT, PWM, PDXX, PDZZ, PDZX, GX_W_UW_PWM)
!
!
!*      13.   < U'W'>
!             -------
! 
! residual part of < U'W'> depending on dw/dx
!
CALL MZM_PHY(D, PK, ZWORKSHUGRAD_3D_1)
CALL MXM_PHY(D, ZWORKSHUGRAD_3D_1, ZWORKSHUGRAD_3D_2)

DO JK=1, D%NKT
  DO JIJ=1, D%NIJT
    ZFLX(JIJ, JK) =                                                      &
      - XCMFS * ZWORKSHUGRAD_3D_2(JIJ, JK) * GX_W_UW_PWM(JIJ, JK)
  END DO
END DO

!!         &  to be tested
!!  - (2./3.) * XCMFB * MZM( ZVPTU * MXM( PLM / SQRT(PTKEM) * XG / PTHVREF ) )
!

DO JIJ=1, IIJT
  ZFLX(JIJ, IKE+1) = 0.  ! rigid wall condition => no turbulent flux
  !
  ! Nullify the flux at the ground level because it has been fully taken into
  ! account in turb_ver and extrapolate the flux under the ground 
  ZFLX(JIJ, IKB) = 0.
  ZFLX(JIJ, IKB-1)=2.*ZFLX(JIJ, IKB)- ZFLX(JIJ, IKB+1)
END DO

!
! stores  <U W>
IF ( TPFILE%LOPENED .AND. TURBN%LTURB_FLX ) THEN
  TZFIELD = TFIELDMETADATA(     &
  CMNHNAME   = 'UW_HFLX',       &
  CSTDNAME   = '',              &
  CLONGNAME  = 'UW_HFLX',       &
  CUNITS     = 'm2 s-2',        &
  CDIR       = 'XY',            &
  CCOMMENT   = 'X_Y_Z_UW_HFLX', &
  NGRID      = 6,               &
  NTYPE      = TYPEREAL,        &
  NDIMS      = 3,               &
  LTIMEDEP   = .TRUE.           )

  CALL IO_FIELD_WRITE_PHY(D,TPFILE,TZFIELD,ZFLX)
END IF
!
!
! compute the source for rho*U due to this residual flux ( the other part is
! taken into account in TURB_VER)
CALL MXM_PHY(D, PMZM_PRHODJ, ZWORKSHUGRAD_3D_1)
CALL MXM_PHY(D, PDZZ, ZWORKSHUGRAD_3D_2)

DO JK=1, D%NKT
  DO JIJ=1, D%NIJT
    ZWORKSHUGRAD_3D_3(JIJ, JK)=ZFLX(JIJ, JK)* ZWORKSHUGRAD_3D_1(JIJ, JK) / ZWORKSHUGRAD_3D_2(JIJ, JK)
  END DO
END DO

CALL DZF_PHY(D, ZWORKSHUGRAD_3D_3, ZWORKSHUGRAD_3D_4)

DO JK=1, D%NKT
  DO JIJ=1, D%NIJT
    PRUS(JIJ, JK) = PRUS(JIJ, JK) - ZWORKSHUGRAD_3D_4(JIJ, JK)
  END DO
END DO

!
!computation of the source for rho*W due to this flux
IF (.NOT. OFLAT) THEN
  CALL MXM_PHY(D, PRHODJ, ZWORKSHUGRAD_3D_1)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      ZWORKSHUGRAD_3D_2(JIJ, JK)=ZFLX(JIJ, JK)*PDZX(JIJ, JK)
    END DO
  END DO
  
  CALL MZF_PHY(D, ZWORKSHUGRAD_3D_2, ZWORKSHUGRAD_3D_3)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      ZWORKSHUGRAD_3D_4(JIJ, JK)=ZWORKSHUGRAD_3D_1(JIJ, JK) * PINV_PDXX(JIJ, JK)
    END DO
  END DO
  
  CALL MZM_PHY(D, ZWORKSHUGRAD_3D_4, ZWORKSHUGRAD_3D_5)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      ZWORKSHUGRAD_3D_6(JIJ, JK)=ZWORKSHUGRAD_3D_3(JIJ, JK) * PINV_PDXX(JIJ, JK)
    END DO
  END DO
  
  CALL MXF_PHY(D, ZWORKSHUGRAD_3D_6, ZWORKSHUGRAD_3D_7)
  CALL MZF_PHY(D, PDZZ, ZWORKSHUGRAD_3D_8)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      ZWORKSHUGRAD_3D_9(JIJ, JK)=ZWORKSHUGRAD_3D_5(JIJ, JK) * ZFLX(JIJ, JK)
    END DO
  END DO
  
  CALL DXF_PHY(D, ZWORKSHUGRAD_3D_9, ZWORKSHUGRAD_3D_10)
  
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
  CALL MXM_PHY(D, PRHODJ, ZWORKSHUGRAD_3D_1)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      ZWORKSHUGRAD_3D_2(JIJ, JK)=ZWORKSHUGRAD_3D_1(JIJ, JK) * PINV_PDXX(JIJ, JK)
    END DO
  END DO
  
  CALL MZM_PHY(D, ZWORKSHUGRAD_3D_2, ZWORKSHUGRAD_3D_3)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      ZWORKSHUGRAD_3D_4(JIJ, JK)=ZWORKSHUGRAD_3D_3(JIJ, JK) * ZFLX(JIJ, JK)
    END DO
  END DO
  
  CALL DXF_PHY(D, ZWORKSHUGRAD_3D_4, ZWORKSHUGRAD_3D_5)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      PRWS(JIJ, JK) = PRWS(JIJ, JK) -ZWORKSHUGRAD_3D_5(JIJ, JK)
    END DO
  END DO
  
END IF
! 
IF (KSPLT==1) THEN
  !
  !Contribution to the dynamic production of TKE:
  !
  CALL GZ_U_UW_PHY(D, PUM, PDZZ, ZWORKSHUGRAD_3D_1)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      ZWORKSHUGRAD_3D_2(JIJ, JK)=ZFLX(JIJ, JK) *( ZWORKSHUGRAD_3D_1(JIJ, JK) + GX_W_UW_PWM(JIJ, JK) )
    END DO
  END DO
  
  CALL MXF_PHY(D, ZWORKSHUGRAD_3D_2, ZWORKSHUGRAD_3D_3)
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
    ZWORKSHUGRAD_2D_1(JIJ)=(PWM(JIJ, IKB+2)-PWM(JIJ, IKB+1))                 &
                  /(PDZZ(JIJ, IKB+2)+PDZZ(JIJ, IKB+1))             &
                +(PWM(JIJ, IKB+1)-PWM(JIJ, IKB))                 &
                  /(PDZZ(JIJ, IKB+1)+PDZZ(JIJ, IKB))
  END DO
  
  CALL MXM2D_PHY(D, ZWORKSHUGRAD_2D_1, ZWORKSHUGRAD_2D_2)
  CALL DXM2D_PHY(D, PWM(:,IKB+1), ZWORKSHUGRAD_2D_3)
  CALL MXM2D_PHY(D, PDZZ(:,IKB+1), ZWORKSHUGRAD_2D_4)
  
  DO JIJ=1, D%NIJT
    ZWORKSHUGRAD_2D_5(JIJ)=ZFLX(JIJ, IKB+1) *                                               &
     (   (PUM(JIJ, IKB+1)-PUM(JIJ, IKB)) / ZWORKSHUGRAD_2D_4(JIJ)&
       + ( ZWORKSHUGRAD_2D_3(JIJ)                                       &
          -ZWORKSHUGRAD_2D_2(JIJ)                                                              &
            * PDZX(JIJ, IKB+1)                                          &
         ) / (0.5*(PDXX(JIJ, IKB+1)+PDXX(JIJ, IKB)))                 &
     )
  END DO
  
  CALL MXF2D_PHY(D, ZWORKSHUGRAD_2D_5, ZWORKSHUGRAD_2D_6)
  
  DO JIJ=1, D%NIJT
    ZWORK(JIJ, IKB) = - ZWORKSHUGRAD_2D_6(JIJ)
  END DO
  
  !
  ! dynamic production computation
  
  DO JK=1, IKT
    DO JIJ=1, IIJT
      PDP(JIJ, JK) = PDP(JIJ, JK) +  ZWORK(JIJ, JK)  
    END DO
  END DO
  
  !
END IF
!
! Storage in the LES configuration (addition to TURB_VER computation)
!
IF (TLES%LLES_CALL .AND. KSPLT==1) THEN
  CALL SECOND_MNH(ZTIME1)
  CALL MXF_PHY(D, ZFLX, ZWORKSHUGRAD_3D_1)
  CALL MZF_PHY(D, ZWORKSHUGRAD_3D_1, ZWORKSHUGRAD_3D_2)
  CALL LES_MEAN_SUBGRID_PHY(D, TLES, ZWORKSHUGRAD_3D_2(1:D%NIJT, 1:D%NKT), TLES%X_LES_SUBGRID_WU , .TRUE. )
  CALL GZ_U_UW_PHY(D, PUM, PDZZ, ZWORKSHUGRAD_3D_1)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      ZWORKSHUGRAD_3D_2(JIJ, JK)=ZWORKSHUGRAD_3D_1(JIJ, JK)*ZFLX(JIJ, JK)
    END DO
  END DO
  
  CALL MXF_PHY(D, ZWORKSHUGRAD_3D_2, ZWORKSHUGRAD_3D_3)
  CALL MZF_PHY(D, ZWORKSHUGRAD_3D_3, ZWORKSHUGRAD_3D_4)
  CALL LES_MEAN_SUBGRID_PHY(D, TLES, ZWORKSHUGRAD_3D_4(1:D%NIJT, 1:D%NKT), TLES%X_LES_RES_ddxa_U_SBG_UaU , .TRUE.)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      ZWORKSHUGRAD_3D_1(JIJ, JK)=GX_W_UW_PWM(JIJ, JK)*ZFLX(JIJ, JK)
    END DO
  END DO
  
  CALL MXF_PHY(D, ZWORKSHUGRAD_3D_1, ZWORKSHUGRAD_3D_2)
  CALL MZF_PHY(D, ZWORKSHUGRAD_3D_2, ZWORKSHUGRAD_3D_3)
  CALL LES_MEAN_SUBGRID_PHY(D, TLES, ZWORKSHUGRAD_3D_3(1:D%NIJT, 1:D%NKT), TLES%X_LES_RES_ddxa_W_SBG_UaW , .TRUE.)
  CALL GX_M_U_PHY(D, OFLAT, PTHLM, PDXX, PDZZ, PDZX, ZWORKSHUGRAD_3D_1)
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
    CALL GX_M_U_PHY(D, OFLAT, PRM(:,:,1), PDXX, PDZZ, PDZX, ZWORKSHUGRAD_3D_1)
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
    CALL GX_M_U_PHY(D, OFLAT, PSVM(:,:,JSV), PDXX, PDZZ, PDZX, ZWORKSHUGRAD_3D_1)
    CALL MZF_PHY(D, ZFLX, ZWORKSHUGRAD_3D_2)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWORKSHUGRAD_3D_3(JIJ, JK)=ZWORKSHUGRAD_3D_1(JIJ, JK)*ZWORKSHUGRAD_3D_2(JIJ, JK)
      END DO
    END DO
    
    CALL MXF_PHY(D, ZWORKSHUGRAD_3D_3, ZWORKSHUGRAD_3D_4)
    CALL LES_MEAN_SUBGRID_PHY(D, TLES, ZWORKSHUGRAD_3D_4(1:D%NIJT, 1:D%NKT), &
                           TLES%X_LES_RES_ddxa_Sv_SBG_UaW(:,:,:,JSV) , .TRUE.)
  END DO
  CALL SECOND_MNH(ZTIME2)
  TLES%XTIME_LES = TLES%XTIME_LES + ZTIME2 - ZTIME1
END IF

!
END SUBROUTINE TURB_HOR_UW
END MODULE MODE_TURB_HOR_UW
