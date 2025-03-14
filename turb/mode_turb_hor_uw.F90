!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
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
USE MODE_GRADIENT_U_PHY, ONLY: GZ_U_UW_PHY
USE MODE_GRADIENT_M_PHY, ONLY: GX_M_U_PHY
USE MODE_GRADIENT_W_PHY, ONLY: GX_W_UW_PHY
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
TYPE(TFILEDATA),          INTENT(INOUT)    ::  TPFILE       ! Output file
!

REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PK          ! Turbulent diffusion doef.
                                                        ! PK = PLM * SQRT(PTKEM)
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PINV_PDXX   ! 1./PDXX
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PINV_PDZZ   ! 1./PDZZ
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PMZM_PRHODJ ! MZM(PRHODJ)
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PDXX, PDZZ, PDZX
                                                         ! Metric coefficients
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PRHODJ       ! density * grid volume
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PTHVREF      ! ref. state VPT       
!
! Variables at t-1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PUM,PWM,PTHLM
REAL, DIMENSION(D%NIT,D%NJT,D%NKT,KRR), INTENT(IN)    ::  PRM
REAL, DIMENSION(D%NIT,D%NJT,D%NKT,KSV), INTENT(IN)    ::  PSVM
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PTKEM        ! TKE at time t- dt
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PLM          ! Turb. mixing length
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(INOUT) ::  PRUS, PRWS
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(INOUT) ::  PDP          ! TKE production terms
!
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
REAL, DIMENSION(D%NIT,D%NJT,D%NKT)  :: GX_W_UW_PWM
!
REAL :: ZTIME1, ZTIME2
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZGX_W_UW3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZMXM3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZMZM3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZSHUGRADWK1_3D
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZDZF3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZMXM3D_WORK2
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZDXF3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZSHUGRADWK2_3D
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZDZM3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZMXF3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZSHUGRADWK3_3D
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZMZF3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZMZF3D_WORK2
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZGZ_U_UW3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT) ::ZSHUGRADWK1_2D
REAL, DIMENSION(D%NIT,D%NJT) ::ZMXF2D_WORK1
REAL, DIMENSION(D%NIT,D%NJT) ::ZSHUGRADWK2_2D
REAL, DIMENSION(D%NIT,D%NJT) ::ZMXM2D_WORK1
REAL, DIMENSION(D%NIT,D%NJT) ::ZDXM2D_WORK1
REAL, DIMENSION(D%NIT,D%NJT) ::ZMXM2D_WORK2
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZGX_M_U3D_WORK1
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
CALL GX_W_UW_PHY(D, OFLAT,PWM(:, :, :),PDXX(:, :, :),PDZZ(:, :, :),PDZX(:, :, :), ZGX_W_UW3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
GX_W_UW_PWM(:, :, :) = ZGX_W_UW3D_WORK1(:, :, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
!
!
!*      13.   < U'W'>
!             -------
! 
! residual part of < U'W'> depending on dw/dx
!
CALL MZM_PHY(D, PK(:, :, :), ZMZM3D_WORK1)
CALL MXM_PHY(D, ZMZM3D_WORK1, ZMXM3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZFLX(:,:,:) =                                                      &
  - XCMFS * ZMXM3D_WORK1(:, :, :) * GX_W_UW_PWM(:, :, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
!!         &  to be tested
!!  - (2./3.) * XCMFB * MZM( ZVPTU * MXM( PLM / SQRT(PTKEM) * XG / PTHVREF ) )
!

DO JJ=1, IJT
  DO JI=1, IIT
    ZFLX(JI, JJ, IKE+1) = 0.  ! rigid wall condition => no turbulent flux
    !
    ! Nullify the flux at the ground level because it has been fully taken into
    ! account in turb_ver and extrapolate the flux under the ground 
    ZFLX(JI, JJ, IKB) = 0.
    ZFLX(JI, JJ, IKB-1)=2.*ZFLX(JI, JJ, IKB)- ZFLX(JI, JJ, IKB+1)
  END DO
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

  CALL IO_FIELD_WRITE(TPFILE,TZFIELD,ZFLX)
END IF
!
!
! compute the source for rho*U due to this residual flux ( the other part is
CALL MXM_PHY(D, PMZM_PRHODJ(:, :, :), ZMXM3D_WORK1)
CALL MXM_PHY(D, PDZZ(:, :, :), ZMXM3D_WORK2)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZSHUGRADWK1_3D(:, :, :) = ZFLX(:, :, :)* ZMXM3D_WORK1(:, :, :) / ZMXM3D_WORK2(:, :, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
CALL DZF_PHY(D, ZSHUGRADWK1_3D, ZDZF3D_WORK1)
! taken into account in TURB_VER)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
PRUS(:,:,:) = PRUS(:,:,:) - ZDZF3D_WORK1(:, :, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
!
!computation of the source for rho*W due to this flux
IF (.NOT. OFLAT) THEN
  CALL MXM_PHY(D, PRHODJ(:, :, :), ZMXM3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZSHUGRADWK2_3D(:, :, :) = ZMXM3D_WORK1(:, :, :) * PINV_PDXX(:, :, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
CALL MZM_PHY(D, ZSHUGRADWK2_3D, ZMZM3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZSHUGRADWK1_3D(:, :, :) = ZMZM3D_WORK1(:, :, :) * ZFLX(:, :, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
CALL DXF_PHY(D, ZSHUGRADWK1_3D, ZDXF3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZSHUGRADWK3_3D(:, :, :) = ZFLX(:, :, :)*PDZX(:, :, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
CALL MZF_PHY(D, ZSHUGRADWK3_3D, ZMZF3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZSHUGRADWK2_3D(:, :, :) = ZMZF3D_WORK1(:, :, :) * PINV_PDXX(:, :, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
CALL MXF_PHY(D, ZSHUGRADWK2_3D, ZMXF3D_WORK1)
CALL MZF_PHY(D, PDZZ(:, :, :), ZMZF3D_WORK2)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZSHUGRADWK1_3D(:, :, :) = PRHODJ(:, :, :) * ZMXF3D_WORK1(:, :, :) / ZMZF3D_WORK2(:, :, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
CALL DZM_PHY(D, ZSHUGRADWK1_3D, ZDZM3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
PRWS(:,:,:) = PRWS(:,:,:)                              &
        -ZDXF3D_WORK1(:, :, :)           &
        +ZDZM3D_WORK1(:, :, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
ELSE
  CALL MXM_PHY(D, PRHODJ(:, :, :), ZMXM3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZSHUGRADWK2_3D(:, :, :) = ZMXM3D_WORK1(:, :, :) * PINV_PDXX(:, :, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
CALL MZM_PHY(D, ZSHUGRADWK2_3D, ZMZM3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZSHUGRADWK1_3D(:, :, :) = ZMZM3D_WORK1(:, :, :) * ZFLX(:, :, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
CALL DXF_PHY(D, ZSHUGRADWK1_3D, ZDXF3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
PRWS(:,:,:) = PRWS(:,:,:) -ZDXF3D_WORK1(:, :, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
END IF
! 
IF (KSPLT==1) THEN
  !
  !Contribution to the dynamic production of TKE:
  !
  CALL GZ_U_UW_PHY(D, PUM(:, :, :),PDZZ(:, :, :), ZGZ_U_UW3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZSHUGRADWK2_3D(:, :, :) = ZFLX(:, :, :) *( ZGZ_U_UW3D_WORK1(:, :, :) + GX_W_UW_PWM(:, :, :) )
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
CALL MXF_PHY(D, ZSHUGRADWK2_3D, ZMXF3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZSHUGRADWK1_3D(:, :, :) = ZMXF3D_WORK1(:, :, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
CALL MZF_PHY(D, ZSHUGRADWK1_3D, ZMZF3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZWORK(:,:,:) =-ZMZF3D_WORK1(:, :, :)  
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
!
  !
  ! evaluate the dynamic production at w(IKB+1) in PDP(IKB)
  
!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
ZSHUGRADWK2_2D(:, :) = PDZZ(:,:,IKB+1)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)

!
CALL MXM2D_PHY(D, ZSHUGRADWK2_2D, ZMXM2D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
ZSHUGRADWK1_2D(:, :) = PWM(:,:,IKB+1)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)

!
CALL DXM2D_PHY(D, ZSHUGRADWK1_2D, ZDXM2D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
ZSHUGRADWK1_2D(:, :) = (PWM(:,:,IKB+2)-PWM(:,:,IKB+1))                  /(PDZZ(:,:,IKB+2)+PDZZ(:,:,IKB+1))              +(PWM(:,:,IKB+1)-PWM(:,:,IKB))                  /(PDZZ(:,:,IKB+1)+PDZZ(:,:,IKB))
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)

!
CALL MXM2D_PHY(D, ZSHUGRADWK1_2D, ZMXM2D_WORK2)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
ZSHUGRADWK1_2D(:, :) = ZFLX(:,:,IKB+1) *                                                (   (PUM(:,:,IKB+1)-PUM(:,:,IKB)) / ZMXM2D_WORK1(:, :) + ( ZDXM2D_WORK1(:, :)                                        -ZMXM2D_WORK2(:, :)                                                                * PDZX(:,:,IKB+1)                                          )  / (0.5*(PDXX(:,:,IKB+1)+PDXX(:,:,IKB)))                 ) 
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)

!
CALL MXF2D_PHY(D, ZSHUGRADWK1_2D, ZMXF2D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
ZWORK(:,:,IKB) = - ZMXF2D_WORK1(:, :)    
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)

!
!
  ! dynamic production computation
  
  DO JK=1, IKT
    DO JJ=1, IJT
      DO JI=1, IIT
        PDP(JI, JJ, JK) = PDP(JI, JJ, JK) +  ZWORK(JI, JJ, JK)  
      END DO
    END DO
  END DO
  
  !
END IF
!
! Storage in the LES configuration (addition to TURB_VER computation)
!
IF (TLES%LLES_CALL .AND. KSPLT==1) THEN
  CALL SECOND_MNH(ZTIME1)
  CALL MXF_PHY(D, ZFLX(:, :, :), ZMXF3D_WORK1)
CALL MZF_PHY(D, ZMXF3D_WORK1, ZMZF3D_WORK1)
CALL LES_MEAN_SUBGRID( ZMZF3D_WORK1(:, :, :), TLES%X_LES_SUBGRID_WU , .TRUE. )  
CALL GZ_U_UW_PHY(D, PUM(:, :, :),PDZZ(:, :, :), ZGZ_U_UW3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZSHUGRADWK2_3D(:, :, :) = ZGZ_U_UW3D_WORK1(:, :, :)*ZFLX(:, :, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
CALL MXF_PHY(D, ZSHUGRADWK2_3D, ZMXF3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZSHUGRADWK1_3D(:, :, :) = ZMXF3D_WORK1(:, :, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
CALL MZF_PHY(D, ZSHUGRADWK1_3D, ZMZF3D_WORK1)
CALL LES_MEAN_SUBGRID( ZMZF3D_WORK1(:, :, :), TLES%X_LES_RES_ddxa_U_SBG_UaU , .TRUE.)  

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZSHUGRADWK2_3D(:, :, :) = GX_W_UW_PWM(:, :, :)*ZFLX(:, :, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
CALL MXF_PHY(D, ZSHUGRADWK2_3D, ZMXF3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZSHUGRADWK1_3D(:, :, :) = ZMXF3D_WORK1(:, :, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
CALL MZF_PHY(D, ZSHUGRADWK1_3D, ZMZF3D_WORK1)
CALL LES_MEAN_SUBGRID( ZMZF3D_WORK1(:, :, :), TLES%X_LES_RES_ddxa_W_SBG_UaW , .TRUE.)  
CALL GX_M_U_PHY(D, OFLAT,PTHLM(:, :, :),PDXX(:, :, :),PDZZ(:, :, :),PDZX(:, :, :), ZGX_M_U3D_WORK1)
CALL MZF_PHY(D, ZFLX(:, :, :), ZMZF3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZSHUGRADWK1_3D(:, :, :) = ZGX_M_U3D_WORK1(:, :, :)*ZMZF3D_WORK1(:, :, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
CALL MXF_PHY(D, ZSHUGRADWK1_3D, ZMXF3D_WORK1)
CALL LES_MEAN_SUBGRID( ZMXF3D_WORK1(:, :, :),&
                         TLES%X_LES_RES_ddxa_Thl_SBG_UaW , .TRUE.)  
IF (KRR>=1) THEN
    CALL GX_M_U_PHY(D, OFLAT,PRM(:,:,:,1),PDXX(:, :, :),PDZZ(:, :, :),PDZX(:, :, :), ZGX_M_U3D_WORK1)
CALL MZF_PHY(D, ZFLX(:, :, :), ZMZF3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZSHUGRADWK1_3D(:, :, :) = ZGX_M_U3D_WORK1(:, :, :)*ZMZF3D_WORK1(:, :, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
CALL MXF_PHY(D, ZSHUGRADWK1_3D, ZMXF3D_WORK1)
CALL LES_MEAN_SUBGRID( ZMXF3D_WORK1(:, :, :), &
                           TLES%X_LES_RES_ddxa_Rt_SBG_UaW , .TRUE.)  
END IF
  DO JSV=1,KSV
    CALL GX_M_U_PHY(D, OFLAT,PSVM(:,:,:,JSV),PDXX(:, :, :),PDZZ(:, :, :),PDZX(:, :, :), ZGX_M_U3D_WORK1)
CALL MZF_PHY(D, ZFLX(:, :, :), ZMZF3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZSHUGRADWK1_3D(:, :, :) = ZGX_M_U3D_WORK1(:, :, :)*ZMZF3D_WORK1(:, :, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
CALL MXF_PHY(D, ZSHUGRADWK1_3D, ZMXF3D_WORK1)
CALL LES_MEAN_SUBGRID( ZMXF3D_WORK1(:, :, :), &
                           TLES%X_LES_RES_ddxa_Sv_SBG_UaW(:,:,:,JSV) , .TRUE.)  
END DO
  CALL SECOND_MNH(ZTIME2)
  TLES%XTIME_LES = TLES%XTIME_LES + ZTIME2 - ZTIME1
END IF

!
END SUBROUTINE TURB_HOR_UW
END MODULE MODE_TURB_HOR_UW
