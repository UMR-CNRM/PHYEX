!MNH_LIC Copyright 1994-2020 CNRS, Meteo-France and Universite Paul Sabatier
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
USE MODD_CST
USE MODD_CTURB
USE MODD_PARAMETERS
USE MODD_LES, ONLY: TLES_t
!
!
USE MODI_SHUMAN
USE MODE_SHUMAN_PHY
USE MODE_GRADIENT_M_PHY, ONLY: GX_M_U_PHY, GY_M_V_PHY
USE MODI_GRADIENT_M
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
TYPE(DIMPHYEX_t),         INTENT(IN)    :: D
TYPE(TLES_t),             INTENT(INOUT) :: TLES          ! modd_les structure
INTEGER,                  INTENT(IN) :: KSPLT        ! current split index
LOGICAL,                  INTENT(IN) ::  OFLAT       ! Logical for zero ororography
LOGICAL,                  INTENT(IN) ::  O2D         ! Logical for 2D model version (modd_conf)
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN) :: PDXX, PDYY, PDZZ, PDZX, PDZY 
                                                     ! Metric coefficients
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN) :: PK           ! Turbulent diffusion doef.
                                                     ! PK = PLM * SQRT(PTKEM)
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN) :: PINV_PDXX    ! 1./PDXX
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN) :: PINV_PDYY    ! 1./PDYY
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN) :: PINV_PDZZ    ! 1./PDZZ
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN) :: PMZM_PRHODJ  ! MZM(PRHODJ)
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN) :: PRHODJ       ! density * grid volume
!
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PTKEM    ! TKE at time t- dt
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(OUT)   ::  PTRH     ! horizontal transport of Tke
!
!
!
!*       0.2  declaration of local variables
!
INTEGER :: IKB, IKU, IIT, IJT, IKT, JI, JJ, JK
!
REAL, DIMENSION(D%NIT,D%NJT,1+JPVEXT:3+JPVEXT) :: ZCOEFF 
                                    ! coefficients for the uncentred gradient 
                                    ! computation near the ground
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT):: ZFLX
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZMXM3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZGX_M_U3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT) ::ZMXM2D_WORK1
REAL, DIMENSION(D%NIT,D%NJT) ::ZDXM2D_WORK1
REAL, DIMENSION(D%NIT,D%NJT) ::ZMXM2D_WORK2
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZSHUGRADWK1_3D
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZDXF3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZDZF3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZSHUGRADWK2_3D
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZMXF3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZSHUGRADWK3_3D
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZMZM3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZMYM3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZGY_M_V3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT) ::ZMYM2D_WORK1
REAL, DIMENSION(D%NIT,D%NJT) ::ZDYM2D_WORK1
REAL, DIMENSION(D%NIT,D%NJT) ::ZMYM2D_WORK2
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZDYF3D_WORK1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::ZMYF3D_WORK1
REAL :: ZTIME1, ZTIME2
! ---------------------------------------------------------------------------
!
!*       1.   PRELIMINARY COMPUTATIONS
!             ------------------------
!
IKB = 1.+JPVEXT
IKU = SIZE(PTKEM,3)
IIT=D%NIT
IJT=D%NJT
IKT=D%NKT 
!
!  compute the coefficients for the uncentred gradient computation near the 
!  ground
!

DO JJ=1, IJT
  DO JI=1, IIT
    ZCOEFF(JI, JJ, IKB+2)= - PDZZ(JI, JJ, IKB+1) /      &
         ( (PDZZ(JI, JJ, IKB+2)+PDZZ(JI, JJ, IKB+1)) * PDZZ(JI, JJ, IKB+2) )
    ZCOEFF(JI, JJ, IKB+1)=   (PDZZ(JI, JJ, IKB+2)+PDZZ(JI, JJ, IKB+1)) /      &
         ( PDZZ(JI, JJ, IKB+1) * PDZZ(JI, JJ, IKB+2) )
    ZCOEFF(JI, JJ, IKB)= - (PDZZ(JI, JJ, IKB+2)+2.*PDZZ(JI, JJ, IKB+1)) /      &
         ( (PDZZ(JI, JJ, IKB+2)+PDZZ(JI, JJ, IKB+1)) * PDZZ(JI, JJ, IKB+1) )
  END DO
END DO

!
!--------------------------------------------------------------------
!
!*       2.   horizontal transport of Tke u'e
!             -------------------------------
!
!
CALL MXM_PHY(D, PK(:, :, :), ZMXM3D_WORK1)
CALL GX_M_U_PHY(D, OFLAT,PTKEM(:, :, :),PDXX(:, :, :),PDZZ(:, :, :),PDZX(:, :, :), ZGX_M_U3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZFLX(:, :, :) = -XCET * ZMXM3D_WORK1(:, :, :) * ZGX_M_U3D_WORK1(:, :, :) 
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
! < u'e >
!
! special case near the ground ( uncentred gradient )
!

DO JJ=1, IJT
  DO JI=1, IIT
    ZFLX(JI, JJ, IKB) =  ZCOEFF(JI, JJ, IKB+2)*PTKEM(JI, JJ, IKB+2)                     &
                   + ZCOEFF(JI, JJ, IKB+1)*PTKEM(JI, JJ, IKB+1)                     &
                   + ZCOEFF(JI, JJ, IKB)*PTKEM(JI, JJ, IKB)     
  END DO
END DO

!
CALL MXM2D_PHY(D, PK(:,:,IKB), ZMXM2D_WORK1)
CALL DXM2D_PHY(D, PTKEM(:,:,IKB), ZDXM2D_WORK1)
CALL MXM2D_PHY(D, ZFLX (:,:,IKB), ZMXM2D_WORK2)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
ZFLX(:,:,IKB) =                                                      &
   - XCET * ZMXM2D_WORK1(:, :)                           *  (        &
       ZDXM2D_WORK1(:, :) * PINV_PDXX(:,:,IKB)               &
      -ZMXM2D_WORK2(:, :) * PINV_PDXX(:,:,IKB)               &
       * 0.5 * ( PDZX(:,:,IKB+1) + PDZX(:,:,IKB) )     ) 
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)

!
!
! extrapolate the fluxes to obtain < u'e > = 0 at the ground
!

DO JJ=1, IJT
  DO JI=1, IIT
    ZFLX(JI, JJ, IKB-1) = - ZFLX(JI, JJ, IKB)
    !
    ! let the same flux at IKU-1 and IKU level
    !
    ZFLX(JI, JJ, IKU) =  ZFLX(JI, JJ, IKU-1)
  END DO
END DO

!
IF (.NOT. OFLAT) THEN
  CALL MXM_PHY(D, PRHODJ(:, :, :), ZMXM3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZSHUGRADWK1_3D(:, :, :) = ZMXM3D_WORK1(:, :, :) * ZFLX(:, :, :)                             * PINV_PDXX(:, :, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
CALL DXF_PHY(D, ZSHUGRADWK1_3D, ZDXF3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZSHUGRADWK3_3D(:, :, :) = ZFLX(:, :, :)*PINV_PDXX(:, :, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
CALL MZM_PHY(D, ZSHUGRADWK3_3D, ZMZM3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZSHUGRADWK2_3D(:, :, :) = PDZX(:, :, :) * ZMZM3D_WORK1(:, :, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
CALL MXF_PHY(D, ZSHUGRADWK2_3D, ZMXF3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZSHUGRADWK1_3D(:, :, :) = PMZM_PRHODJ(:, :, :) * ZMXF3D_WORK1(:, :, :) * PINV_PDZZ(:, :, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
CALL DZF_PHY(D, ZSHUGRADWK1_3D, ZDZF3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
PTRH(:, :, :) =-(  ZDXF3D_WORK1(:, :, :)&
          - ZDZF3D_WORK1(:, :, :)&
         ) /PRHODJ(:, :, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
ELSE
  CALL MXM_PHY(D, PRHODJ(:, :, :), ZMXM3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZSHUGRADWK1_3D(:, :, :) = ZMXM3D_WORK1(:, :, :) * ZFLX(:, :, :)                             * PINV_PDXX(:, :, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
CALL DXF_PHY(D, ZSHUGRADWK1_3D, ZDXF3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
PTRH(:, :, :) =-(  ZDXF3D_WORK1(:, :, :)&
         ) /PRHODJ(:, :, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
END IF
!
IF (TLES%LLES_CALL .AND. KSPLT==1) THEN
  CALL SECOND_MNH(ZTIME1)
  CALL MXF_PHY(D, ZFLX(:, :, :), ZMXF3D_WORK1)
CALL LES_MEAN_SUBGRID( ZMXF3D_WORK1(:, :, :), TLES%X_LES_SUBGRID_UTke )   
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
  CALL MYM_PHY(D, PK(:, :, :), ZMYM3D_WORK1)
CALL GY_M_V_PHY(D, OFLAT,PTKEM(:, :, :),PDYY(:, :, :),PDZZ(:, :, :),PDZY(:, :, :), ZGY_M_V3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZFLX(:, :, :) =-XCET * ZMYM3D_WORK1(:, :, :) * ZGY_M_V3D_WORK1(:, :, :) 
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
! < v'e >
!
! special case near the ground ( uncentred gradient )
!

DO JJ=1, IJT
    DO JI=1, IIT
      ZFLX(JI, JJ, IKB) =  ZCOEFF(JI, JJ, IKB+2)*PTKEM(JI, JJ, IKB+2)                     &
                     + ZCOEFF(JI, JJ, IKB+1)*PTKEM(JI, JJ, IKB+1)                     &
                     + ZCOEFF(JI, JJ, IKB)*PTKEM(JI, JJ, IKB)     
    !
  END DO
  END DO

  CALL MYM2D_PHY(D, PK(:,:,IKB), ZMYM2D_WORK1)
CALL DYM2D_PHY(D, PTKEM(:,:,IKB), ZDYM2D_WORK1)
CALL MYM2D_PHY(D, ZFLX (:,:,IKB), ZMYM2D_WORK2)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
ZFLX(:,:,IKB) =                                                      &
     - XCET * ZMYM2D_WORK1(:, :)                        *  (           &
       ZDYM2D_WORK1(:, :) * PINV_PDYY(:,:,IKB)                 &
     - ZMYM2D_WORK2(:, :) * PINV_PDYY(:,:,IKB)                 &
         * 0.5 * ( PDZY(:,:,IKB+1) + PDZY(:,:,IKB) )  )
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)

!
!
!    extrapolate the fluxes to obtain < v'e > = 0 at the ground
!

DO JJ=1, IJT
    DO JI=1, IIT
      ZFLX(JI, JJ, IKB-1) = - ZFLX(JI, JJ, IKB)
    !
    !   let the same flux at IKU-1 and IKU level
    !
      ZFLX(JI, JJ, IKU) =  ZFLX(JI, JJ, IKU-1)
  END DO
  END DO

!
! complete the explicit turbulent transport
!
  IF (.NOT. OFLAT) THEN
    CALL MYM_PHY(D, PRHODJ(:, :, :), ZMYM3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZSHUGRADWK1_3D(:, :, :) = ZMYM3D_WORK1(:, :, :) * ZFLX(:, :, :)                              * PINV_PDYY(:, :, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
CALL DYF_PHY(D, ZSHUGRADWK1_3D, ZDYF3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZSHUGRADWK3_3D(:, :, :) = ZFLX(:, :, :)*PINV_PDYY(:, :, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
CALL MZM_PHY(D, ZSHUGRADWK3_3D, ZMZM3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZSHUGRADWK2_3D(:, :, :) = PDZY(:, :, :) * ZMZM3D_WORK1(:, :, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
CALL MYF_PHY(D, ZSHUGRADWK2_3D, ZMYF3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZSHUGRADWK1_3D(:, :, :) = PMZM_PRHODJ(:, :, :) * ZMYF3D_WORK1(:, :, :) * PINV_PDZZ(:, :, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
CALL DZF_PHY(D, ZSHUGRADWK1_3D, ZDZF3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
PTRH(:, :, :) = PTRH(:, :, :) - (  ZDYF3D_WORK1(:, :, :)  &
                   - ZDZF3D_WORK1(:, :, :)  &
                  ) /PRHODJ(:, :, :)  
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
ELSE
    CALL MYM_PHY(D, PRHODJ(:, :, :), ZMYM3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZSHUGRADWK1_3D(:, :, :) = ZMYM3D_WORK1(:, :, :) * ZFLX(:, :, :)                              * PINV_PDYY(:, :, :)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
CALL DYF_PHY(D, ZSHUGRADWK1_3D, ZDYF3D_WORK1)

!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
PTRH(:, :, :) = PTRH(:, :, :) - (  ZDYF3D_WORK1(:, :, :)  &
                  ) /PRHODJ(:, :, :)  
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)

!
END IF
!
  IF (TLES%LLES_CALL .AND. KSPLT==1) THEN
    CALL SECOND_MNH(ZTIME1)
    CALL MYF_PHY(D, ZFLX(:, :, :), ZMYF3D_WORK1)
CALL LES_MEAN_SUBGRID( ZMYF3D_WORK1(:, :, :), TLES%X_LES_SUBGRID_VTke )    
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
