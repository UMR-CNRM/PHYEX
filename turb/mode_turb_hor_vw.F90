!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
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
USE MODE_GRADIENT_V_PHY, ONLY: GZ_V_VW_PHY
USE MODE_GRADIENT_M_PHY, ONLY: GY_M_V_PHY
USE MODE_GRADIENT_W_PHY, ONLY: GY_W_VW_PHY
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
INTEGER,                  INTENT(IN)    ::  KRR          ! number of moist var.
INTEGER,                  INTENT(IN)    ::  KSV          ! number of sv var.
LOGICAL,                  INTENT(IN)    ::  OFLAT        ! Logical for zero ororography
LOGICAL,                  INTENT(IN)    ::  O2D          ! Logical for 2D model version (modd_conf)
TYPE(TFILEDATA),          INTENT(INOUT)    ::  TPFILE       ! Output file
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PK          ! Turbulent diffusion doef.
                                                        ! PK = PLM * SQRT(PTKEM)
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PINV_PDYY   ! 1./PDYY
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PINV_PDZZ   ! 1./PDZZ
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PMZM_PRHODJ ! MZM(PRHODJ)
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PDYY, PDZZ, PDZY 
                                                         ! Metric coefficients
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PRHODJ       ! density * grid volume
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PTHVREF      ! ref. state VPT       
!
! Variables at t-1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PVM,PWM,PTHLM
REAL, DIMENSION(D%NIT,D%NJT,D%NKT,KRR), INTENT(IN)    ::  PRM
REAL, DIMENSION(D%NIT,D%NJT,D%NKT,KSV), INTENT(IN)    ::  PSVM
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PTKEM        ! TKE at time t- dt
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PLM          ! Turb. mixing length
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(INOUT) ::  PRVS, PRWS   ! var. at t+1 -split-
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(INOUT) ::  PDP          ! TKE production terms
!
!
!
!*       0.2  declaration of local variables
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT)  :: ZFLX,ZWORK
    ! work arrays
!   
INTEGER             :: IKB,IKE,IKU, IIT, IJT, IKT
                                    ! Index values for the Beginning and End
                                    ! mass points of the domain  
INTEGER             :: JSV,JI,JJ,JK          ! scalar loop counter
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT)  :: GY_W_VW_PWM
!
REAL :: ZTIME1, ZTIME2
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZGY_W_VW3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZMYM3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZMZM3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZSHUGRADWK1_3D
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZDZF3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZMYM3D_WORK2
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZDYF3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZSHUGRADWK2_3D
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZDZM3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZMYF3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZSHUGRADWK3_3D
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZMZF3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZMZF3D_WORK2
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZGZ_V_VW3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT) ::ZSHUGRADWK1_2D
REAL, DIMENSION(D%NIT,D%NJT) ::ZMYF2D_WORK1
REAL, DIMENSION(D%NIT,D%NJT) ::ZSHUGRADWK2_2D
REAL, DIMENSION(D%NIT,D%NJT) ::ZMYM2D_WORK1
REAL, DIMENSION(D%NIT,D%NJT) ::ZDYM2D_WORK1
REAL, DIMENSION(D%NIT,D%NJT) ::ZMYM2D_WORK2
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZMXF3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZGY_M_V3D_WORK1
TYPE(TFIELDMETADATA) :: TZFIELD
! ---------------------------------------------------------------------------
!
!*       1.   PRELIMINARY COMPUTATIONS
!             ------------------------
!
IKB = 1+JPVEXT               
IKE = SIZE(PWM,3)-JPVEXT    
IKU = SIZE(PWM,3)
IIT=D%NIT
IJT=D%NJT
IKT=D%NKT 
!
!
IF (.NOT. O2D) THEN
  CALL GY_W_VW_PHY(D, OFLAT,PWM(:, :, :),PDYY(:, :, :),PDZZ(:, :, :),PDZY(:, :, :), ZGY_W_VW3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      GY_W_VW_PWM(JI, JJ, JK) = ZGY_W_VW3D_WORK1(JI, JJ, JK)
    END DO
  END DO
END DO

!
END IF
!
!
!*      14.   < V'W'>
!             -------
!
! residual part of < V'W'> depending on dw/dy
!
IF (.NOT. O2D) THEN
  CALL MZM_PHY(D, PK(:, :, :), ZMZM3D_WORK1)
CALL MYM_PHY(D, ZMZM3D_WORK1, ZMYM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZFLX(JI, JJ, JK) =                                                      &
          - XCMFS * ZMYM3D_WORK1(JI, JJ, JK) * GY_W_VW_PWM(JI, JJ, JK)
    END DO
  END DO
END DO

!
ELSE

DO JK=1, IKT
    DO JJ=1, IJT
      DO JI=1, IIT
        ZFLX(JI, JJ, JK) = 0.
    END DO
    END DO
  END DO

END IF
!

DO JJ=1, IJT
  DO JI=1, IIT
    ZFLX(JI, JJ, IKE+1) = 0.  ! rigid wall condition => no turbulent flux
    !
    !
    ! Nullify the flux at the ground level because it has been fully taken into
    ! account in turb_ver and extrapolate the flux under the ground 
    ZFLX(JI, JJ, IKB) = 0.
    ZFLX(JI, JJ, IKB-1)= 2.*ZFLX(JI, JJ, IKB) - ZFLX(JI, JJ, IKB+1)
  END DO
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

  CALL IO_FIELD_WRITE(TPFILE,TZFIELD,ZFLX)
END IF
!
! compute the source for rho*V due to this residual flux ( the other part is
! taken into account in TURB_VER)
IF (.NOT. O2D) THEN
  CALL MYM_PHY(D, PMZM_PRHODJ(:, :, :), ZMYM3D_WORK1)
CALL MYM_PHY(D, PDZZ(:, :, :), ZMYM3D_WORK2)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZFLX(JI, JJ, JK)* ZMYM3D_WORK1(JI, JJ, JK) / ZMYM3D_WORK2(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL DZF_PHY(D, ZSHUGRADWK1_3D, ZDZF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      PRVS(JI, JJ, JK) = PRVS(JI, JJ, JK) - ZDZF3D_WORK1(JI, JJ, JK)
    END DO
  END DO
END DO

!
END IF
!
!computation of the source for rho*W due to this flux
IF (.NOT. O2D) THEN 
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
CALL MZM_PHY(D, ZSHUGRADWK2_3D, ZMZM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMZM3D_WORK1(JI, JJ, JK) * ZFLX(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL DYF_PHY(D, ZSHUGRADWK1_3D, ZDYF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK3_3D(JI, JJ, JK) = ZFLX(JI, JJ, JK)*PDZY(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MZF_PHY(D, ZSHUGRADWK3_3D, ZMZF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = ZMZF3D_WORK1(JI, JJ, JK) * PINV_PDYY(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MYF_PHY(D, ZSHUGRADWK2_3D, ZMYF3D_WORK1)
CALL MZF_PHY(D, PDZZ(:, :, :), ZMZF3D_WORK2)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = PRHODJ(JI, JJ, JK) * ZMYF3D_WORK1(JI, JJ, JK) / ZMZF3D_WORK2(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL DZM_PHY(D, ZSHUGRADWK1_3D, ZDZM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      PRWS(JI, JJ, JK) = PRWS(JI, JJ, JK)                              &
                -ZDYF3D_WORK1(JI, JJ, JK)           &
                +ZDZM3D_WORK1(JI, JJ, JK)  
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
CALL MZM_PHY(D, ZSHUGRADWK2_3D, ZMZM3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZMZM3D_WORK1(JI, JJ, JK) * ZFLX(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL DYF_PHY(D, ZSHUGRADWK1_3D, ZDYF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      PRWS(JI, JJ, JK) = PRWS(JI, JJ, JK) - ZDYF3D_WORK1(JI, JJ, JK)  
    END DO
  END DO
END DO

!
END IF
END IF
!
IF (KSPLT==1) THEN
  ! 
  !Contribution to the dynamic production of TKE:
  !
  IF (.NOT. O2D) THEN
    CALL GZ_V_VW_PHY(D, PVM(:, :, :),PDZZ(:, :, :), ZGZ_V_VW3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = ZFLX(JI, JJ, JK) *( ZGZ_V_VW3D_WORK1(JI, JJ, JK) + GY_W_VW_PWM(JI, JJ, JK) )
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

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZWORK(JI, JJ, JK) =-ZMZF3D_WORK1(JI, JJ, JK)  
    END DO
  END DO
END DO

!
!
  !
  ! evaluate the dynamic production at w(IKB+1) in PDP(IKB)
    
DO JJ=1, IJT
  DO JI=1, IIT
    ZSHUGRADWK2_2D(JI, JJ) = PDZZ(JI, JJ, IKB+1)
  END DO
END DO

!
CALL MYM2D_PHY(D, ZSHUGRADWK2_2D, ZMYM2D_WORK1)

DO JJ=1, IJT
  DO JI=1, IIT
    ZSHUGRADWK1_2D(JI, JJ) = PWM(JI, JJ, IKB+1)
  END DO
END DO

!
CALL DYM2D_PHY(D, ZSHUGRADWK1_2D, ZDYM2D_WORK1)

DO JJ=1, IJT
  DO JI=1, IIT
    ZSHUGRADWK1_2D(JI, JJ) = (PWM(JI, JJ, IKB+2)-PWM(JI, JJ, IKB+1))                    /(PDZZ(JI, JJ, IKB+2)+PDZZ(JI, JJ, IKB+1))                +(PWM(JI, JJ, IKB+1)-PWM(JI, JJ, IKB))                    /(PDZZ(JI, JJ, IKB+1)+PDZZ(JI, JJ, IKB))
  END DO
END DO

!
CALL MYM2D_PHY(D, ZSHUGRADWK1_2D, ZMYM2D_WORK2)

DO JJ=1, IJT
  DO JI=1, IIT
    ZSHUGRADWK1_2D(JI, JJ) = ZFLX(JI, JJ, IKB+1) *                                                  (   (PVM(JI, JJ, IKB+1)-PVM(JI, JJ, IKB)) / ZMYM2D_WORK1(JI, JJ)   + ( ZDYM2D_WORK1(JI, JJ)                                          -ZMYM2D_WORK2(JI, JJ)  * PDZY(JI, JJ, IKB+1)                                        )  / (0.5*(PDYY(JI, JJ, IKB+1)+PDYY(JI, JJ, IKB)))                 ) 
  END DO
END DO

!
CALL MYF2D_PHY(D, ZSHUGRADWK1_2D, ZMYF2D_WORK1)

DO JJ=1, IJT
  DO JI=1, IIT
    ZWORK(JI, JJ, IKB) = - ZMYF2D_WORK1(JI, JJ)    
  END DO
END DO

!
ENDIF
  !
  ! dynamic production computation
  IF (.NOT. O2D) THEN
    
    DO JK=1, IKT
      DO JJ=1, IJT
        DO JI=1, IIT
          PDP(JI, JJ, JK) = PDP(JI, JJ, JK) + ZWORK(JI, JJ, JK)  
        END DO
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
  CALL MYF_PHY(D, ZFLX(:, :, :), ZMYF3D_WORK1)
CALL MZF_PHY(D, ZMYF3D_WORK1, ZMZF3D_WORK1)
CALL LES_MEAN_SUBGRID( ZMZF3D_WORK1(:, :, :), TLES%X_LES_SUBGRID_WV , .TRUE. )  
CALL GZ_V_VW_PHY(D, PVM(:, :, :),PDZZ(:, :, :), ZGZ_V_VW3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = ZGZ_V_VW3D_WORK1(JI, JJ, JK)*ZFLX(JI, JJ, JK)
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
CALL LES_MEAN_SUBGRID( ZMZF3D_WORK1(:, :, :),&
                         TLES%X_LES_RES_ddxa_V_SBG_UaV , .TRUE.)  

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK2_3D(JI, JJ, JK) = GY_W_VW_PWM(JI, JJ, JK)*ZFLX(JI, JJ, JK)
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
CALL LES_MEAN_SUBGRID( ZMZF3D_WORK1(:, :, :),&
                         TLES%X_LES_RES_ddxa_W_SBG_UaW , .TRUE.)  
CALL GY_M_V_PHY(D, OFLAT,PTHLM(:, :, :),PDYY(:, :, :),PDZZ(:, :, :),PDZY(:, :, :), ZGY_M_V3D_WORK1)
CALL MZF_PHY(D, ZFLX(:, :, :), ZMZF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZGY_M_V3D_WORK1(JI, JJ, JK)*ZMZF3D_WORK1(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MXF_PHY(D, ZSHUGRADWK1_3D, ZMXF3D_WORK1)
CALL LES_MEAN_SUBGRID( ZMXF3D_WORK1(:, :, :),&
                         TLES%X_LES_RES_ddxa_Thl_SBG_UaW , .TRUE.)  
IF (KRR>=1) THEN
    CALL GY_M_V_PHY(D, OFLAT,PRM(:,:,:,1),PDYY(:, :, :),PDZZ(:, :, :),PDZY(:, :, :), ZGY_M_V3D_WORK1)
CALL MZF_PHY(D, ZFLX(:, :, :), ZMZF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZGY_M_V3D_WORK1(JI, JJ, JK)*ZMZF3D_WORK1(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MXF_PHY(D, ZSHUGRADWK1_3D, ZMXF3D_WORK1)
CALL LES_MEAN_SUBGRID( ZMXF3D_WORK1(:, :, :), &
                           TLES%X_LES_RES_ddxa_Rt_SBG_UaW , .TRUE.)  
END IF
  DO JSV=1,KSV
    CALL GY_M_V_PHY(D, OFLAT,PSVM(:,:,:,JSV),PDYY(:, :, :),PDZZ(:, :, :),PDZY(:, :, :), ZGY_M_V3D_WORK1)
CALL MZF_PHY(D, ZFLX(:, :, :), ZMZF3D_WORK1)

DO JK=1, IKT
  DO JJ=1, IJT
    DO JI=1, IIT
      ZSHUGRADWK1_3D(JI, JJ, JK) = ZGY_M_V3D_WORK1(JI, JJ, JK)*ZMZF3D_WORK1(JI, JJ, JK)
    END DO
  END DO
END DO

!
CALL MXF_PHY(D, ZSHUGRADWK1_3D, ZMXF3D_WORK1)
CALL LES_MEAN_SUBGRID( ZMXF3D_WORK1(:, :, :), &
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
