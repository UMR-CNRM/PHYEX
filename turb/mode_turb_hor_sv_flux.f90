!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_TURB_HOR_SV_FLUX
IMPLICIT NONE
CONTAINS
      SUBROUTINE TURB_HOR_SV_FLUX(TURBN,TLES,KSPLT,OBLOWSNOW,OFLAT,  &
                      TPFILE,KSV_LGBEG,KSV_LGEND,O2D,ONOMIXLG,       &
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
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!! --------------------------------------------------------------------------
!
!*      0. DECLARATIONS
!          ------------
!
USE MODD_TURB_n, ONLY: TURB_t
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
USE MODI_SHUMAN 
USE MODE_COEFJ, ONLY: COEFJ
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
TYPE(TURB_t),             INTENT(IN)    :: TURBN
TYPE(TLES_t),             INTENT(INOUT) :: TLES          ! modd_les structure
INTEGER,                  INTENT(IN)    ::  KSPLT        ! split process index
TYPE(TFILEDATA),          INTENT(IN)    ::  TPFILE       ! Output file
INTEGER,                  INTENT(IN)    ::  KSV_LGBEG,KSV_LGEND ! number of sv var.
LOGICAL,                  INTENT(IN)    ::  OFLAT        ! Logical for zero ororography
LOGICAL,                  INTENT(IN)    ::  ONOMIXLG     ! to use turbulence for lagrangian variables (modd_conf)
LOGICAL,                  INTENT(IN)    ::  O2D          ! Logical for 2D model version (modd_conf)
LOGICAL,                  INTENT(IN)    ::  OBLOWSNOW    ! switch to activate pronostic blowing snow
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PK          ! Turbulent diffusion doef.
                                                        ! PK = PLM * SQRT(PTKEM)
REAL,                     INTENT(IN)    ::  PRSNOW       ! Ratio for diffusion coeff. scalar (blowing snow)
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PINV_PDXX   ! 1./PDXX
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PINV_PDYY   ! 1./PDYY
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PINV_PDZZ   ! 1./PDZZ
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PMZM_PRHODJ ! MZM(PRHODJ)
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PDXX, PDYY, PDZZ, PDZX, PDZY 
                                                         ! Metric coefficients
REAL, DIMENSION(:,:),     INTENT(IN)    ::  PDIRCOSXW, PDIRCOSYW
! Director Cosinus along x and y  directions at surface w-point
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PRHODJ       ! density * grid volume
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PWM          ! vertical wind
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PSFSVM       ! surface fluxes
!
!
! Variables at t-1
REAL, DIMENSION(:,:,:,:), INTENT(IN)    ::  PSVM         ! scalar var. at t-1
!
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) ::  PRSVS        ! var. at t+1 -split-
!
!
!
!*       0.2  declaration of local variables
!
REAL, DIMENSION(SIZE(PSVM,1),SIZE(PSVM,2),SIZE(PSVM,3))       &
                                     :: ZFLXX,ZFLXY
    ! work arrays
REAL, DIMENSION(SIZE(PSVM,1),SIZE(PSVM,2),1) :: ZWORK2D
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
INTEGER :: IKU
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
!
ISV = SIZE(PSVM,4)
!
IF(OBLOWSNOW) THEN
! See Vionnet (PhD, 2012) for a complete discussion around the value of the Schmidt number for blowing snow variables              
   ZCSV= XCHF/PRSNOW
ELSE
   ZCSV= XCHF
ENDIF
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
  ZFLXX(:,:,:) = -ZCSV * MXM(PK) * GX_M_U(1,IKU,1,PSVM(:,:,:,JSV),PDXX,PDZZ,PDZX)
  ZFLXX(:,:,IKE+1) = ZFLXX(:,:,IKE) 
!
! Compute the flux at the first inner U-point with an uncentred vertical  
! gradient
  ZFLXX(:,:,IKB:IKB) = -ZCSV * MXM( PK(:,:,IKB:IKB) ) *             &
    ( DXM(PSVM(:,:,IKB:IKB,JSV)) * PINV_PDXX(:,:,IKB:IKB)           &
     -MXM ( ZCOEFF(:,:,IKB+2:IKB+2)*PSVM(:,:,IKB+2:IKB+2,JSV)       &
           +ZCOEFF(:,:,IKB+1:IKB+1)*PSVM(:,:,IKB+1:IKB+1,JSV)       &
           +ZCOEFF(:,:,IKB  :IKB  )*PSVM(:,:,IKB  :IKB  ,JSV)       &
          ) * 0.5 * ( PDZX(:,:,IKB+1:IKB+1)+PDZX(:,:,IKB:IKB) )     &
            * PINV_PDXX(:,:,IKB:IKB)                                &
    ) 
! extrapolates the flux under the ground so that the vertical average with 
! the IKB flux gives the ground value
  ZWORK2D(:,:,1)=PSFSVM(:,:,JSV) * PDIRCOSXW(:,:)
  ZFLXX(:,:,IKB-1:IKB-1) = 2. * MXM( ZWORK2D(:,:,1:1) ) - ZFLXX(:,:,IKB:IKB)
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
    CALL IO_Field_write(TPFILE,TZFIELD,ZFLXX)
  END IF
!
  IF (TLES%LLES_CALL .AND. KSPLT==1) THEN
    CALL SECOND_MNH(ZTIME1)
    CALL LES_MEAN_SUBGRID( MXF(ZFLXX), TLES%X_LES_SUBGRID_USv(:,:,:,JSV) ) 
    CALL LES_MEAN_SUBGRID( MZF(MXF(GX_W_UW(PWM,PDXX,PDZZ,PDZX)*MZM(ZFLXX))), &
                           TLES%X_LES_RES_ddxa_W_SBG_UaSv(:,:,:,JSV) , .TRUE. )
    CALL LES_MEAN_SUBGRID( GX_M_M(PSVM(:,:,:,JSV),PDXX,PDZZ,PDZX)*MXF(ZFLXX), &
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
    ZFLXY(:,:,:)=-ZCSV * MYM(PK) * GY_M_V(1,IKU,1,PSVM(:,:,:,JSV),PDYY,PDZZ,PDZY)
    ZFLXY(:,:,IKE+1) = ZFLXY(:,:,IKE) 
!
! Compute the flux at the first inner V-point with an uncentred vertical  
! gradient
!
    ZFLXY(:,:,IKB:IKB) = -ZCSV * MYM( PK(:,:,IKB:IKB) ) *             &
      ( DYM(PSVM(:,:,IKB:IKB,JSV)) * PINV_PDYY(:,:,IKB:IKB)           &
       -MYM ( ZCOEFF(:,:,IKB+2:IKB+2)*PSVM(:,:,IKB+2:IKB+2,JSV)       &
             +ZCOEFF(:,:,IKB+1:IKB+1)*PSVM(:,:,IKB+1:IKB+1,JSV)       &
             +ZCOEFF(:,:,IKB  :IKB  )*PSVM(:,:,IKB  :IKB  ,JSV)       &
            ) * 0.5 * ( PDZY(:,:,IKB+1:IKB+1)+PDZY(:,:,IKB:IKB) )     &
              * PINV_PDYY(:,:,IKB:IKB)                                &
      ) 
! extrapolates the flux under the ground so that the vertical average with 
! the IKB flux gives the ground value
    ZWORK2D(:,:,1)=PSFSVM(:,:,JSV) * PDIRCOSYW(:,:)
    ZFLXY(:,:,IKB-1:IKB-1) = 2. * MYM( ZWORK2D(:,:,1:1) ) - ZFLXY(:,:,IKB:IKB)
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
      CALL IO_Field_write(TPFILE,TZFIELD,ZFLXY)
    END IF
!
  ELSE
    ZFLXY=0.
  END IF
!
  IF (TLES%LLES_CALL .AND. KSPLT==1) THEN
    CALL SECOND_MNH(ZTIME1)
    CALL LES_MEAN_SUBGRID( MYF(ZFLXY), TLES%X_LES_SUBGRID_VSv(:,:,:,JSV) ) 
    CALL LES_MEAN_SUBGRID( MZF(MYF(GY_W_VW(PWM,PDYY,PDZZ,PDZY)*MZM(ZFLXY))), &
                           TLES%X_LES_RES_ddxa_W_SBG_UaSv(:,:,:,JSV) , .TRUE. )
    CALL LES_MEAN_SUBGRID( GY_M_M(PSVM(:,:,:,JSV),PDYY,PDZZ,PDZY)*MYF(ZFLXY), &
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
