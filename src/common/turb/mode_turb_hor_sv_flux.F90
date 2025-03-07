!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_TURB_HOR_SV_FLUX
IMPLICIT NONE
CONTAINS
      SUBROUTINE TURB_HOR_SV_FLUX(D,TURBN,TLES,KSPLT,OBLOWSNOW,OFLAT,&
                      TPFILE,KSV,KSV_LGBEG,KSV_LGEND,O2D,ONOMIXLG,   &
                      PK,PINV_PDXX,PINV_PDYY,PINV_PDZZ,PMZM_PRHODJ,  &
                      PDXX,PDYY,PDZZ,PDZX,PDZY,PRSNOW,               &
                      PDIRCOSXW,PDIRCOSYW,                           &
                      PRHODJ,PWM,                                    &
                      PSFSVM,                                        &
                      PSVM,                                          &
                      PRSVS                                          )
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
!!                     Oct  18, 2000 (V. Masson) LES computations + OFLAT swith
!!                                              + bug on Y scalar flux
!!                     Jun  20, 2001 (J Stein) case of lagragian variables
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
use modd_field,          only: tfieldmetadata, TYPEREAL
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
USE MODE_GRADIENT_M_PHY, ONLY: GX_M_U_PHY, GY_M_V_PHY, GX_M_M_PHY, GY_M_M_PHY
USE MODE_GRADIENT_W_PHY, ONLY: GX_W_UW_PHY, GY_W_VW_PHY
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
TYPE(TFILEDATA),          INTENT(INOUT)    ::  TPFILE       ! Output file
INTEGER,                  INTENT(IN)    ::  KSV,KSV_LGBEG,KSV_LGEND ! number of sv var.
LOGICAL,                  INTENT(IN)    ::  OFLAT        ! Logical for zero ororography
LOGICAL,                  INTENT(IN)    ::  ONOMIXLG     ! to use turbulence for lagrangian variables (modd_conf)
LOGICAL,                  INTENT(IN)    ::  O2D          ! Logical for 2D model version (modd_conf)
LOGICAL,                  INTENT(IN)    ::  OBLOWSNOW    ! switch to activate pronostic blowing snow
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PK          ! Turbulent diffusion doef.
                                                        ! PK = PLM * SQRT(PTKEM)
REAL,                     INTENT(IN)    ::  PRSNOW       ! Ratio for diffusion coeff. scalar (blowing snow)
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PINV_PDXX   ! 1./PDXX
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PINV_PDYY   ! 1./PDYY
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PINV_PDZZ   ! 1./PDZZ
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PMZM_PRHODJ ! MZM(PRHODJ)
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PDXX, PDYY, PDZZ, PDZX, PDZY 
                                                         ! Metric coefficients
REAL, DIMENSION(:,:),     INTENT(IN)    ::  PDIRCOSXW, PDIRCOSYW
! Director Cosinus along x and y  directions at surface w-point
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PRHODJ       ! density * grid volume
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PWM          ! vertical wind
!
REAL, DIMENSION(D%NIT,D%NJT,KSV),   INTENT(IN)    ::  PSFSVM       ! surface fluxes
!
!
! Variables at t-1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT,KSV), INTENT(IN)    ::  PSVM         ! scalar var. at t-1
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT,KSV), INTENT(INOUT) ::  PRSVS        ! var. at t+1 -split-
!
!
!
!*       0.2  declaration of local variables
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT)       &
                                     :: ZFLXX,ZFLXY
    ! work arrays
REAL, DIMENSION(D%NIT,D%NJT) :: ZWORK2D
!
REAL :: ZCSV          !constant for the scalar flux

INTEGER             :: IKB,IKE
                                    ! Index values for the Beginning and End
                                    ! mass points of the domain  
INTEGER             :: JSV          ! loop counter
INTEGER             :: ISV          ! number of scalar var.
REAL, DIMENSION(SIZE(PDZZ,1),SIZE(PDZZ,2),1+JPVEXT:3+JPVEXT) :: ZCOEFF 
                                    ! coefficients for the uncentred gradient 
                                    ! computation near the ground
!
CHARACTER(LEN=NMNHNAMELGTMAX) :: YMNHNAME
INTEGER :: IKU, IIT, IJT, IKT, JI, JJ, JK
TYPE(TFIELDMETADATA) :: TZFIELD
REAL :: ZTIME1, ZTIME2
! ---------------------------------------------------------------------------
!
!*       1.   PRELIMINARY COMPUTATIONS
!             ------------------------
!
IKB = 1+JPVEXT               
IKE = SIZE(PSVM,3)-JPVEXT   
IKU = SIZE(PSVM,3)
IIT=D%NIT
IJT=D%NJT
IKT=D%NKT
!
ISV = SIZE(PSVM,4)
!
IF(OBLOWSNOW) THEN
! See Vionnet (PhD, 2012) for a complete discussion around the value of the Schmidt number for blowing snow variables              
   ZCSV= TURBN%XCHF/PRSNOW
ELSE
   ZCSV= TURBN%XCHF
ENDIF
!
!  compute the coefficients for the uncentred gradient computation near the 
!  ground
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
!
!*      15.   HORIZONTAL FLUXES OF PASSIVE SCALARS
!             ------------------------------------
!
!
DO JSV=1,ISV
!
  IF (ONOMIXLG .AND. JSV >= KSV_LGBEG .AND. JSV<= KSV_LGEND) CYCLE
!
!       15.1   <U' SVth'>
!              ----------
!
  ! Computes the flux in the X direction
  ZFLXX(:,:,:) = -ZCSV * MXM(PK) * GX_M_U(OFLAT,PSVM(:,:,:,JSV),PDXX,PDZZ,PDZX)

!$acc kernels
!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
  ZFLXX(:,:,IKE+1) = ZFLXX(:,:,IKE)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)
!$acc end kernels
!
! Compute the flux at the first inner U-point with an uncentred vertical  
! gradient
  ZFLXX(:,:,IKB) = -ZCSV * MXM( PK(:,:,IKB) ) *             &
    ( DXM(PSVM(:,:,IKB,JSV)) * PINV_PDXX(:,:,IKB)           &
     -MXM ( ZCOEFF(:,:,2)*PSVM(:,:,2,JSV)       &
           +ZCOEFF(:,:,IKB+1)*PSVM(:,:,IKB+1,JSV)       &
           +ZCOEFF(:,:,IKB  )*PSVM(:,:,IKB  ,JSV)       &
          ) * 0.5 * ( PDZX(:,:,IKB+1)+PDZX(:,:,IKB) )     &
            * PINV_PDXX(:,:,IKB)                                &
    ) 
! extrapolates the flux under the ground so that the vertical average with 
! the IKB flux gives the ground value
!$acc kernels
!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
  ZWORK2D(:,:)=PSFSVM(:,:,JSV) * PDIRCOSXW(:,:)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)
!$acc end kernels
  ZFLXX(:,:,IKB-1) = 2. * MXM( ZWORK2D(:,:) ) - ZFLXX(:,:,IKB)
  !
  ! stores  <U SVth>
  IF ( TPFILE%LOPENED .AND. TURBN%LTURB_FLX ) THEN
    WRITE(YMNHNAME,'("USV_FLX_",I3.3)') JSV
    TZFIELD = TFIELDMETADATA(                    &
      CMNHNAME   = TRIM( YMNHNAME ),             &
      CSTDNAME   = '',                           &
      CLONGNAME  = TRIM( YMNHNAME ),             &
      CUNITS     = 'SVUNIT m s-1',               &
      CDIR       = 'XY',                         &
      CCOMMENT   = 'X_Y_Z_' // TRIM( YMNHNAME ), &
      NGRID      = 2,                            &
      NTYPE      = TYPEREAL,                     &
      NDIMS      = 3,                            &
      LTIMEDEP   = .TRUE.                        )
!$acc update self(ZFLXX)
    CALL IO_Field_write(TPFILE,TZFIELD,ZFLXX)
  END IF
!
  IF (TLES%LLES_CALL .AND. KSPLT==1) THEN
    CALL SECOND_MNH(ZTIME1)
    CALL LES_MEAN_SUBGRID( MXF(ZFLXX), TLES%X_LES_SUBGRID_USv(:,:,:,JSV) ) 
    CALL LES_MEAN_SUBGRID( MZF(MXF(GX_W_UW(OFLAT,PWM,PDXX,PDZZ,PDZX)*MZM(ZFLXX))), &
                           TLES%X_LES_RES_ddxa_W_SBG_UaSv(:,:,:,JSV) , .TRUE. )
    CALL LES_MEAN_SUBGRID( GX_M_M(OFLAT,PSVM(:,:,:,JSV),PDXX,PDZZ,PDZX)*MXF(ZFLXX), &
                           TLES%X_LES_RES_ddxa_Sv_SBG_UaSv(:,:,:,JSV), .TRUE. )
    CALL SECOND_MNH(ZTIME2)
    TLES%XTIME_LES = TLES%XTIME_LES + ZTIME2 - ZTIME1
  END IF
!
!       15.2   <V' SVth'>
!              ----------
!
  IF (.NOT. O2D) THEN
!
! Computes the flux in the Y direction
    ZFLXY(:,:,:)=-ZCSV * MYM(PK) * GY_M_V(OFLAT,PSVM(:,:,:,JSV),PDYY,PDZZ,PDZY)
!$acc kernels
!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
    ZFLXY(:,:,IKE+1) = ZFLXY(:,:,IKE)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)
!$acc end kernels
!
! Compute the flux at the first inner V-point with an uncentred vertical  
! gradient
!
    ZFLXY(:,:,IKB) = -ZCSV * MYM( PK(:,:,IKB) ) *             &
      ( DYM(PSVM(:,:,IKB,JSV)) * PINV_PDYY(:,:,IKB)           &
       -MYM ( ZCOEFF(:,:,2)*PSVM(:,:,2,JSV)       &
             +ZCOEFF(:,:,IKB+1)*PSVM(:,:,IKB+1,JSV)       &
             +ZCOEFF(:,:,IKB  )*PSVM(:,:,IKB  ,JSV)       &
            ) * 0.5 * ( PDZY(:,:,IKB+1)+PDZY(:,:,IKB) )     &
              * PINV_PDYY(:,:,IKB)                                &
      ) 
! extrapolates the flux under the ground so that the vertical average with 
! the IKB flux gives the ground value
!$acc kernels
!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
    ZWORK2D(:,:)=PSFSVM(:,:,JSV) * PDIRCOSYW(:,:)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)
!$acc end kernels
    ZFLXY(:,:,IKB-1) = 2. * MYM( ZWORK2D(:,:) ) - ZFLXY(:,:,IKB)
  !
  ! stores  <V SVth>
    IF ( TPFILE%LOPENED .AND. TURBN%LTURB_FLX ) THEN
      WRITE(YMNHNAME,'("VSV_FLX_",I3.3)') JSV
    TZFIELD = TFIELDMETADATA(                        &
      CMNHNAME   = TRIM( YMNHNAME ),                 &
      CSTDNAME   = '',                               &
      CLONGNAME  = TRIM(TZFIELD%CMNHNAME),           &
      CUNITS     = 'SVUNIT m s-1',                   &
      CDIR       = 'XY',                             &
      CCOMMENT   = 'X_Y_Z_'//TRIM(TZFIELD%CMNHNAME), &
      NGRID      = 3,                                &
      NTYPE      = TYPEREAL,                         &
      NDIMS      = 3,                                &
      LTIMEDEP   = .TRUE.                            )
!$acc update self(ZFLXY)
      CALL IO_Field_write(TPFILE,TZFIELD,ZFLXY)
    END IF
!
  ELSE
    !$acc kernels
    !$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
    ZFLXY(:,:,:)=0.
    !$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
    !$acc end kernels
  END IF
!
  IF (TLES%LLES_CALL .AND. KSPLT==1) THEN
    CALL SECOND_MNH(ZTIME1)
    CALL LES_MEAN_SUBGRID( MYF(ZFLXY), TLES%X_LES_SUBGRID_VSv(:,:,:,JSV) ) 
    CALL LES_MEAN_SUBGRID( MZF(MYF(GY_W_VW(OFLAT,PWM,PDYY,PDZZ,PDZY)*MZM(ZFLXY))), &
                           TLES%X_LES_RES_ddxa_W_SBG_UaSv(:,:,:,JSV) , .TRUE. )
    CALL LES_MEAN_SUBGRID( GY_M_M(OFLAT,PSVM(:,:,:,JSV),PDYY,PDZZ,PDZY)*MYF(ZFLXY), &
                           TLES%X_LES_RES_ddxa_Sv_SBG_UaSv(:,:,:,JSV) , .TRUE. )
    CALL SECOND_MNH(ZTIME2)
    TLES%XTIME_LES = TLES%XTIME_LES + ZTIME2 - ZTIME1
  END IF
!
!
!       15.3   Horizontal source terms
!              -----------------------
!
  IF (.NOT. O2D) THEN
    IF (.NOT. OFLAT) THEN
      PRSVS(:,:,:,JSV)=   PRSVS(:,:,:,JSV)                                          &
        -DXF( MXM(PRHODJ) * ZFLXX * PINV_PDXX  )                                    &
        -DYF( MYM(PRHODJ) * ZFLXY * PINV_PDYY  )                                    &
        +DZF( PMZM_PRHODJ * PINV_PDZZ *                                             &
              ( MXF( MZM(ZFLXX * PINV_PDXX) * PDZX ) + MYF( MZM(ZFLXY * PINV_PDYY) * PDZY ) ) &
            )
    ELSE
      PRSVS(:,:,:,JSV)=   PRSVS(:,:,:,JSV)                                          &
        -DXF( MXM(PRHODJ) * ZFLXX * PINV_PDXX  )                                    &
        -DYF( MYM(PRHODJ) * ZFLXY * PINV_PDYY  )
    END IF
  ELSE
    IF (.NOT. OFLAT) THEN
      PRSVS(:,:,:,JSV)=   PRSVS(:,:,:,JSV)                                          &
        -DXF( MXM(PRHODJ) * ZFLXX * PINV_PDXX  )                                    &
        +DZF( PMZM_PRHODJ * PINV_PDZZ *                                             &
              ( MXF( MZM(ZFLXX * PINV_PDXX) * PDZX ) )                              &
            )
    ELSE
      PRSVS(:,:,:,JSV)=   PRSVS(:,:,:,JSV)                                          &
        -DXF( MXM(PRHODJ) *ZFLXX * PINV_PDXX  )
    END IF
  END IF
!
!
END DO    ! end loop JSV
!
!
END SUBROUTINE TURB_HOR_SV_FLUX
END MODULE MODE_TURB_HOR_SV_FLUX
