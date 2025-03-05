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
!$acc kernels
!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
ZCOEFF(:,:,IKB+2)= - PDZZ(:,:,IKB+1) /      &
     ( (PDZZ(:,:,IKB+2)+PDZZ(:,:,IKB+1)) * PDZZ(:,:,IKB+2) )
ZCOEFF(:,:,IKB+1)=   (PDZZ(:,:,IKB+2)+PDZZ(:,:,IKB+1)) /      &
     ( PDZZ(:,:,IKB+1) * PDZZ(:,:,IKB+2) )
ZCOEFF(:,:,IKB)= - (PDZZ(:,:,IKB+2)+2.*PDZZ(:,:,IKB+1)) /      &
     ( (PDZZ(:,:,IKB+2)+PDZZ(:,:,IKB+1)) * PDZZ(:,:,IKB+1) )
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)
!$acc end kernels
!
!--------------------------------------------------------------------
!
!*       2.   horizontal transport of Tke u'e
!             -------------------------------
!
!
ZFLX = -XCET * MXM(PK) * GX_M_U(OFLAT,PTKEM,PDXX,PDZZ,PDZX) ! < u'e >
!
! special case near the ground ( uncentred gradient )
!
!$acc kernels
!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
ZFLX(:,:,IKB) =  ZCOEFF(:,:,IKB+2)*PTKEM(:,:,IKB+2)                     &
               + ZCOEFF(:,:,IKB+1)*PTKEM(:,:,IKB+1)                     &
               + ZCOEFF(:,:,IKB  )*PTKEM(:,:,IKB  )     
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)
!$acc end kernels
!
ZFLX(:,:,IKB) =                                                      &
   - XCET * MXM( PK(:,:,IKB) )                           *  (        &
       DXM ( PTKEM(:,:,IKB) ) * PINV_PDXX(:,:,IKB)               &
      -MXM ( ZFLX (:,:,IKB) ) * PINV_PDXX(:,:,IKB)               &
       * 0.5 * ( PDZX(:,:,IKB+1) + PDZX(:,:,IKB) )     ) 
!
! extrapolate the fluxes to obtain < u'e > = 0 at the ground
!
!$acc kernels
!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
ZFLX(:,:,IKB-1) = - ZFLX(:,:,IKB)
!
! let the same flux at IKU-1 and IKU level
!
ZFLX(:,:,IKU) =  ZFLX(:,:,IKU-1)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)
!$acc end kernels 
!
IF (.NOT. OFLAT) THEN
  PTRH =-(  DXF( MXM(PRHODJ) * ZFLX                             * PINV_PDXX)&
          - DZF( PMZM_PRHODJ * MXF( PDZX * MZM(ZFLX*PINV_PDXX)) * PINV_PDZZ)&
         ) /PRHODJ
ELSE
  PTRH =-(  DXF( MXM(PRHODJ) * ZFLX                             * PINV_PDXX)&
         ) /PRHODJ
END IF
!
IF (TLES%LLES_CALL .AND. KSPLT==1) THEN
  CALL SECOND_MNH(ZTIME1)
  CALL LES_MEAN_SUBGRID( MXF(ZFLX), TLES%X_LES_SUBGRID_UTke ) 
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
  ZFLX =-XCET * MYM(PK) * GY_M_V(OFLAT,PTKEM,PDYY,PDZZ,PDZY) ! < v'e >
!
! special case near the ground ( uncentred gradient )
!
!$acc kernels
!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
  ZFLX(:,:,IKB) =  ZCOEFF(:,:,IKB+2)*PTKEM(:,:,IKB+2)                     &
                 + ZCOEFF(:,:,IKB+1)*PTKEM(:,:,IKB+1)                     &
                 + ZCOEFF(:,:,IKB  )*PTKEM(:,:,IKB  )     
!
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)
!$acc end kernels
  ZFLX(:,:,IKB) =                                                      &
     - XCET * MYM( PK(:,:,IKB) )                        *  (           &
       DYM ( PTKEM(:,:,IKB) ) * PINV_PDYY(:,:,IKB)                 &
     - MYM ( ZFLX (:,:,IKB) ) * PINV_PDYY(:,:,IKB)                 &
         * 0.5 * ( PDZY(:,:,IKB+1) + PDZY(:,:,IKB) )  )
!
!    extrapolate the fluxes to obtain < v'e > = 0 at the ground
!
!$acc kernels
!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
  ZFLX(:,:,IKB-1) = - ZFLX(:,:,IKB)
!
!   let the same flux at IKU-1 and IKU level
!
  ZFLX(:,:,IKU) =  ZFLX(:,:,IKU-1)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)
!$acc end kernels
!
! complete the explicit turbulent transport
!
  IF (.NOT. OFLAT) THEN
    PTRH = PTRH - (  DYF( MYM(PRHODJ) * ZFLX                              * PINV_PDYY )  &
                   - DZF( PMZM_PRHODJ * MYF( PDZY * MZM(ZFLX*PINV_PDYY) ) * PINV_PDZZ )  &
                  ) /PRHODJ
  ELSE
    PTRH = PTRH - (  DYF( MYM(PRHODJ) * ZFLX                              * PINV_PDYY )  &
                  ) /PRHODJ
  END IF
!
  IF (TLES%LLES_CALL .AND. KSPLT==1) THEN
    CALL SECOND_MNH(ZTIME1)
    CALL LES_MEAN_SUBGRID( MYF(ZFLX), TLES%X_LES_SUBGRID_VTke )
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
