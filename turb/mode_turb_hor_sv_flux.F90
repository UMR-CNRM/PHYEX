!MNH_LIC Copyright 1994-2025 CNRS, Meteo-France and Universite Paul Sabatier
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
use modd_field,          only: tfieldmetadata, TYPEREAL
USE MODD_IO,             ONLY: TFILEDATA
USE MODD_PARAMETERS, ONLY: JPVEXT, NMNHNAMELGTMAX
USE MODD_LES, ONLY: TLES_t
!
USE MODE_IO_FIELD_WRITE_PHY, ONLY: IO_FIELD_WRITE_PHY
!
USE MODI_LES_MEAN_SUBGRID_PHY, ONLY: LES_MEAN_SUBGRID_PHY
!
USE MODI_SECOND_MNH, ONLY: SECOND_MNH
USE MODE_SHUMAN_PHY, ONLY: DXF_PHY, DXM2D_PHY, DYF_PHY, DYM2D_PHY, DZF_PHY, MXF_PHY, MXM2D_PHY, MXM_PHY, MYF_PHY, MYM2D_PHY, &
  & MYM_PHY, MZF_PHY, MZM_PHY
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
INTEGER,                  INTENT(IN)    ::  KSPLT        ! split process index
TYPE(TFILEDATA),          INTENT(INOUT)    ::  TPFILE       ! Output file
INTEGER,                  INTENT(IN)    ::  KSV,KSV_LGBEG,KSV_LGEND ! number of sv var.
LOGICAL,                  INTENT(IN)    ::  OFLAT        ! Logical for zero ororography
LOGICAL,                  INTENT(IN)    ::  ONOMIXLG     ! to use turbulence for lagrangian variables (modd_conf)
LOGICAL,                  INTENT(IN)    ::  O2D          ! Logical for 2D model version (modd_conf)
LOGICAL,                  INTENT(IN)    ::  OBLOWSNOW    ! switch to activate pronostic blowing snow
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    ::  PK          ! Turbulent diffusion doef.
                                                        ! PK = PLM * SQRT(PTKEM)
REAL,                     INTENT(IN)    ::  PRSNOW       ! Ratio for diffusion coeff. scalar (blowing snow)
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    ::  PINV_PDXX   ! 1./PDXX
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    ::  PINV_PDYY   ! 1./PDYY
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    ::  PINV_PDZZ   ! 1./PDZZ
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    ::  PMZM_PRHODJ ! MZM(PRHODJ)
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    ::  PDXX, PDYY, PDZZ, PDZX, PDZY 
                                                         ! Metric coefficients
REAL, DIMENSION(D%NIJT),     INTENT(IN)    ::  PDIRCOSXW, PDIRCOSYW
! Director Cosinus along x and y  directions at surface w-point
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    ::  PRHODJ       ! density * grid volume
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    ::  PWM          ! vertical wind
!
REAL, DIMENSION(D%NIJT,KSV),   INTENT(IN)    ::  PSFSVM       ! surface fluxes
!
!
! Variables at t-1
REAL, DIMENSION(D%NIJT,D%NKT,KSV), INTENT(IN)    ::  PSVM         ! scalar var. at t-1
!
REAL, DIMENSION(D%NIJT,D%NKT,KSV), INTENT(INOUT) ::  PRSVS        ! var. at t+1 -split-
!
!
!
!*       0.2  declaration of local variables
!
REAL, DIMENSION(D%NIJT,D%NKT)       &
                                     :: ZFLXX,ZFLXY,ZWKLES
    ! work arrays
REAL, DIMENSION(D%NIJT) :: ZWORK2D
!
REAL :: ZCSV          !constant for the scalar flux

INTEGER             :: IKB,IKE
                                    ! Index values for the Beginning and End
                                    ! mass points of the domain  
INTEGER             :: JSV          ! loop counter
INTEGER             :: ISV          ! number of scalar var.
REAL, DIMENSION(D%NIJT,1+JPVEXT:3+JPVEXT) :: ZCOEFF 
                                    ! coefficients for the uncentred gradient 
                                    ! computation near the ground
!
CHARACTER(LEN=NMNHNAMELGTMAX) :: YMNHNAME
INTEGER :: IKU, IIJT, IKT, JIJ, JK
TYPE(TFIELDMETADATA) :: TZFIELD
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
REAL, DIMENSION(D%NIJT, D%NKT) :: ZWORKSHUGRAD_3D_13
REAL, DIMENSION(D%NIJT, D%NKT) :: ZWORKSHUGRAD_3D_14
REAL, DIMENSION(D%NIJT, D%NKT) :: ZWORKSHUGRAD_3D_15
REAL, DIMENSION(D%NIJT, D%NKT) :: ZWORKSHUGRAD_3D_16
REAL :: ZTIME1, ZTIME2
! ---------------------------------------------------------------------------
!
!*       1.   PRELIMINARY COMPUTATIONS
!             ------------------------
!
IKB = 1+JPVEXT               
IKE = SIZE(PSVM,2)-JPVEXT   
IKU = SIZE(PSVM,2)
IIJT=D%NIJT
IKT=D%NKT
!
ISV = SIZE(PSVM,3)
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

DO JIJ=1, IIJT
  ZCOEFF(JIJ, IKB+2)= - PDZZ(JIJ, IKB+1) /      &
         ( (PDZZ(JIJ, IKB+2)+PDZZ(JIJ, IKB+1)) * PDZZ(JIJ, IKB+2) )
  ZCOEFF(JIJ, IKB+1)=   (PDZZ(JIJ, IKB+2)+PDZZ(JIJ, IKB+1)) /      &
         ( PDZZ(JIJ, IKB+1) * PDZZ(JIJ, IKB+2) )
  ZCOEFF(JIJ, IKB)= - (PDZZ(JIJ, IKB+2)+2.*PDZZ(JIJ, IKB+1)) /      &
         ( (PDZZ(JIJ, IKB+2)+PDZZ(JIJ, IKB+1)) * PDZZ(JIJ, IKB+1) )
END DO

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
  CALL MXM_PHY(D, PK, ZWORKSHUGRAD_3D_1)
  CALL GX_M_U_PHY(D, OFLAT, PSVM(:,:,JSV), PDXX, PDZZ, PDZX, ZWORKSHUGRAD_3D_2)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      ZFLXX(JIJ, JK) = -ZCSV * ZWORKSHUGRAD_3D_1(JIJ, JK) * ZWORKSHUGRAD_3D_2(JIJ, JK)
    END DO
  END DO
  

DO JIJ=1, IIJT
    ZFLXX(JIJ, IKE+1) = ZFLXX(JIJ, IKE)
END DO

!
! Compute the flux at the first inner U-point with an uncentred vertical  
! gradient
  CALL DXM2D_PHY(D, PSVM(:,IKB,JSV), ZWORKSHUGRAD_2D_1)
  
  DO JIJ=1, D%NIJT
    ZWORKSHUGRAD_2D_2(JIJ)=ZCOEFF(JIJ, 2)*PSVM(JIJ, 2, JSV)       &
             +ZCOEFF(JIJ, IKB+1)*PSVM(JIJ, IKB+1, JSV)       &
             +ZCOEFF(JIJ, IKB)*PSVM(JIJ, IKB, JSV)
  END DO
  
  CALL MXM2D_PHY(D, ZWORKSHUGRAD_2D_2, ZWORKSHUGRAD_2D_3)
  CALL MXM2D_PHY(D, PK(:,IKB), ZWORKSHUGRAD_2D_4)
  
  DO JIJ=1, D%NIJT
    ZFLXX(JIJ, IKB) = -ZCSV * ZWORKSHUGRAD_2D_4(JIJ) *             &
      ( ZWORKSHUGRAD_2D_1(JIJ) * PINV_PDXX(JIJ, IKB)           &
       -ZWORKSHUGRAD_2D_3(JIJ) * 0.5 * ( PDZX(JIJ, IKB+1)+PDZX(JIJ, IKB) )     &
              * PINV_PDXX(JIJ, IKB)                                &
      )
  END DO
  
! extrapolates the flux under the ground so that the vertical average with 
! the IKB flux gives the ground value

DO JIJ=1, IIJT
    ZWORK2D(JIJ)=PSFSVM(JIJ, JSV) * PDIRCOSXW(JIJ)
END DO

  CALL MXM2D_PHY(D, ZWORK2D, ZWORKSHUGRAD_2D_1)
  
  DO JIJ=1, D%NIJT
    ZFLXX(JIJ, IKB-1) = 2. * ZWORKSHUGRAD_2D_1(JIJ) - ZFLXX(JIJ, IKB)
  END DO
  
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

    CALL IO_FIELD_WRITE_PHY(D, TPFILE,TZFIELD,ZFLXX)
  END IF
!
  IF (TLES%LLES_CALL .AND. KSPLT==1) THEN
    CALL SECOND_MNH(ZTIME1)
    CALL MXF_PHY(D, ZFLXX, ZWORKSHUGRAD_3D_1)
    CALL LES_MEAN_SUBGRID_PHY(D, TLES, ZWORKSHUGRAD_3D_1(1:D%NIJT, 1:D%NKT), TLES%X_LES_SUBGRID_USv(:,:,:,JSV) ) 
    CALL GX_W_UW_PHY(D, OFLAT, PWM, PDXX, PDZZ, PDZX, ZWORKSHUGRAD_3D_1)
    CALL MZM_PHY(D, ZFLXX, ZWORKSHUGRAD_3D_2)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWORKSHUGRAD_3D_3(JIJ, JK)=ZWORKSHUGRAD_3D_1(JIJ, JK)*ZWORKSHUGRAD_3D_2(JIJ, JK)
      END DO
    END DO
    
    CALL MXF_PHY(D, ZWORKSHUGRAD_3D_3, ZWORKSHUGRAD_3D_4)
    CALL MZF_PHY(D, ZWORKSHUGRAD_3D_4, ZWKLES)
    CALL LES_MEAN_SUBGRID_PHY(D, TLES, ZWKLES, TLES%X_LES_RES_ddxa_W_SBG_UaSv(:,:,:,JSV) , .TRUE. )
    CALL GX_M_M_PHY(D, OFLAT, PSVM(:,:,JSV), PDXX, PDZZ, PDZX, ZWORKSHUGRAD_3D_1)
    CALL MXF_PHY(D, ZFLXX, ZWORKSHUGRAD_3D_2)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWKLES(JIJ, JK) = ZWORKSHUGRAD_3D_1(JIJ, JK)*ZWORKSHUGRAD_3D_2(JIJ, JK)
      END DO
    END DO
    
    CALL LES_MEAN_SUBGRID_PHY(D, TLES, ZWKLES, TLES%X_LES_RES_ddxa_Sv_SBG_UaSv(:,:,:,JSV), .TRUE. )
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
    CALL MYM_PHY(D, PK, ZWORKSHUGRAD_3D_1)
    CALL GY_M_V_PHY(D, OFLAT, PSVM(:,:,JSV), PDYY, PDZZ, PDZY, ZWORKSHUGRAD_3D_2)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZFLXY(JIJ, JK)=-ZCSV * ZWORKSHUGRAD_3D_1(JIJ, JK) * ZWORKSHUGRAD_3D_2(JIJ, JK)
      END DO
    END DO
    

DO JIJ=1, IIJT
      ZFLXY(JIJ, IKE+1) = ZFLXY(JIJ, IKE)
END DO

!
! Compute the flux at the first inner V-point with an uncentred vertical  
! gradient
!
    CALL DYM2D_PHY(D, PSVM(:,IKB,JSV), ZWORKSHUGRAD_2D_1)
    
    DO JIJ=1, D%NIJT
      ZWORKSHUGRAD_2D_2(JIJ)=ZCOEFF(JIJ, 2)*PSVM(JIJ, 2, JSV)       &
               +ZCOEFF(JIJ, IKB+1)*PSVM(JIJ, IKB+1, JSV)       &
               +ZCOEFF(JIJ, IKB)*PSVM(JIJ, IKB, JSV)
    END DO
    
    CALL MYM2D_PHY(D, ZWORKSHUGRAD_2D_2, ZWORKSHUGRAD_2D_3)
    CALL MYM2D_PHY(D, PK(:,IKB), ZWORKSHUGRAD_2D_4)
    
    DO JIJ=1, D%NIJT
      ZFLXY(JIJ, IKB) = -ZCSV * ZWORKSHUGRAD_2D_4(JIJ) *             &
        ( ZWORKSHUGRAD_2D_1(JIJ) * PINV_PDYY(JIJ, IKB)           &
         -ZWORKSHUGRAD_2D_3(JIJ) * 0.5 * ( PDZY(JIJ, IKB+1)+PDZY(JIJ, IKB) )     &
                * PINV_PDYY(JIJ, IKB)                                &
        )
    END DO
    
! extrapolates the flux under the ground so that the vertical average with 
! the IKB flux gives the ground value

DO JIJ=1, IIJT
      ZWORK2D(JIJ)=PSFSVM(JIJ, JSV) * PDIRCOSYW(JIJ)
END DO

    CALL MYM2D_PHY(D, ZWORK2D, ZWORKSHUGRAD_2D_1)
    
    DO JIJ=1, D%NIJT
      ZFLXY(JIJ, IKB-1) = 2. * ZWORKSHUGRAD_2D_1(JIJ) - ZFLXY(JIJ, IKB)
    END DO
    
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

      CALL IO_FIELD_WRITE_PHY(D,TPFILE,TZFIELD,ZFLXY)
    END IF
!
  ELSE
    
    DO JK=1, IKT
      DO JIJ=1, IIJT
        ZFLXY(JIJ, JK)=0.
      END DO
    END DO
    
  END IF
!
  IF (TLES%LLES_CALL .AND. KSPLT==1) THEN
    CALL SECOND_MNH(ZTIME1)
    CALL MYF_PHY(D, ZFLXY, ZWORKSHUGRAD_3D_1)
    CALL LES_MEAN_SUBGRID_PHY(D, TLES, ZWORKSHUGRAD_3D_1(1:D%NIJT, 1:D%NKT), TLES%X_LES_SUBGRID_VSv(:,:,:,JSV) ) 
    CALL GY_W_VW_PHY(D, OFLAT, PWM, PDYY, PDZZ, PDZY, ZWORKSHUGRAD_3D_1)
    CALL MZM_PHY(D, ZFLXY, ZWORKSHUGRAD_3D_2)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWORKSHUGRAD_3D_3(JIJ, JK)=ZWORKSHUGRAD_3D_1(JIJ, JK)*ZWORKSHUGRAD_3D_2(JIJ, JK)
      END DO
    END DO
    
    CALL MYF_PHY(D, ZWORKSHUGRAD_3D_3, ZWORKSHUGRAD_3D_4)
    CALL MZF_PHY(D, ZWORKSHUGRAD_3D_4, ZWKLES)
    CALL LES_MEAN_SUBGRID_PHY(D, TLES, ZWKLES, TLES%X_LES_RES_ddxa_W_SBG_UaSv(:,:,:,JSV) , .TRUE. )
    CALL GY_M_M_PHY(D, OFLAT, PSVM(:,:,JSV), PDYY, PDZZ, PDZY, ZWORKSHUGRAD_3D_1)
    CALL MYF_PHY(D, ZFLXY, ZWORKSHUGRAD_3D_2)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWKLES(JIJ, JK) = ZWORKSHUGRAD_3D_1(JIJ, JK)*ZWORKSHUGRAD_3D_2(JIJ, JK)
      END DO
    END DO
    
    CALL LES_MEAN_SUBGRID_PHY(D, TLES, ZWKLES, TLES%X_LES_RES_ddxa_Sv_SBG_UaSv(:,:,:,JSV) , .TRUE. )
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
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_1(JIJ, JK)=ZFLXX(JIJ, JK) * PINV_PDXX(JIJ, JK)
        END DO
      END DO
      
      CALL MZM_PHY(D, ZWORKSHUGRAD_3D_1, ZWORKSHUGRAD_3D_2)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_3(JIJ, JK)=ZFLXY(JIJ, JK) * PINV_PDYY(JIJ, JK)
        END DO
      END DO
      
      CALL MZM_PHY(D, ZWORKSHUGRAD_3D_3, ZWORKSHUGRAD_3D_4)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_5(JIJ, JK)=ZWORKSHUGRAD_3D_2(JIJ, JK) * PDZX(JIJ, JK)
        END DO
      END DO
      
      CALL MXF_PHY(D, ZWORKSHUGRAD_3D_5, ZWORKSHUGRAD_3D_6)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_7(JIJ, JK)=ZWORKSHUGRAD_3D_4(JIJ, JK) * PDZY(JIJ, JK)
        END DO
      END DO
      
      CALL MYF_PHY(D, ZWORKSHUGRAD_3D_7, ZWORKSHUGRAD_3D_8)
      CALL MXM_PHY(D, PRHODJ, ZWORKSHUGRAD_3D_9)
      CALL MYM_PHY(D, PRHODJ, ZWORKSHUGRAD_3D_10)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_11(JIJ, JK)=ZWORKSHUGRAD_3D_9(JIJ, JK) * ZFLXX(JIJ, JK) * PINV_PDXX(JIJ, JK)
        END DO
      END DO
      
      CALL DXF_PHY(D, ZWORKSHUGRAD_3D_11, ZWORKSHUGRAD_3D_12)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_13(JIJ, JK)=ZWORKSHUGRAD_3D_10(JIJ, JK) * ZFLXY(JIJ, JK) * PINV_PDYY(JIJ, JK)
        END DO
      END DO
      
      CALL DYF_PHY(D, ZWORKSHUGRAD_3D_13, ZWORKSHUGRAD_3D_14)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_15(JIJ, JK)=PMZM_PRHODJ(JIJ, JK) * PINV_PDZZ(JIJ, JK) *                                             &
                  ( ZWORKSHUGRAD_3D_6(JIJ, JK) + ZWORKSHUGRAD_3D_8(JIJ, JK) )
        END DO
      END DO
      
      CALL DZF_PHY(D, ZWORKSHUGRAD_3D_15, ZWORKSHUGRAD_3D_16)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          PRSVS(JIJ, JK, JSV)=   PRSVS(JIJ, JK, JSV)                                          &
            -ZWORKSHUGRAD_3D_12(JIJ, JK)                                    &
            -ZWORKSHUGRAD_3D_14(JIJ, JK)                                    &
            +ZWORKSHUGRAD_3D_16(JIJ, JK)
        END DO
      END DO
      
    ELSE
      CALL MXM_PHY(D, PRHODJ, ZWORKSHUGRAD_3D_1)
      CALL MYM_PHY(D, PRHODJ, ZWORKSHUGRAD_3D_2)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_3(JIJ, JK)=ZWORKSHUGRAD_3D_1(JIJ, JK) * ZFLXX(JIJ, JK) * PINV_PDXX(JIJ, JK)
        END DO
      END DO
      
      CALL DXF_PHY(D, ZWORKSHUGRAD_3D_3, ZWORKSHUGRAD_3D_4)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_5(JIJ, JK)=ZWORKSHUGRAD_3D_2(JIJ, JK) * ZFLXY(JIJ, JK) * PINV_PDYY(JIJ, JK)
        END DO
      END DO
      
      CALL DYF_PHY(D, ZWORKSHUGRAD_3D_5, ZWORKSHUGRAD_3D_6)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          PRSVS(JIJ, JK, JSV)=   PRSVS(JIJ, JK, JSV)                                          &
            -ZWORKSHUGRAD_3D_4(JIJ, JK)                                    &
            -ZWORKSHUGRAD_3D_6(JIJ, JK)
        END DO
      END DO
      
    END IF
  ELSE
    IF (.NOT. OFLAT) THEN
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_1(JIJ, JK)=ZFLXX(JIJ, JK) * PINV_PDXX(JIJ, JK)
        END DO
      END DO
      
      CALL MZM_PHY(D, ZWORKSHUGRAD_3D_1, ZWORKSHUGRAD_3D_2)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_3(JIJ, JK)=ZWORKSHUGRAD_3D_2(JIJ, JK) * PDZX(JIJ, JK)
        END DO
      END DO
      
      CALL MXF_PHY(D, ZWORKSHUGRAD_3D_3, ZWORKSHUGRAD_3D_4)
      CALL MXM_PHY(D, PRHODJ, ZWORKSHUGRAD_3D_5)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_6(JIJ, JK)=ZWORKSHUGRAD_3D_5(JIJ, JK) * ZFLXX(JIJ, JK) * PINV_PDXX(JIJ, JK)
        END DO
      END DO
      
      CALL DXF_PHY(D, ZWORKSHUGRAD_3D_6, ZWORKSHUGRAD_3D_7)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_8(JIJ, JK)=PMZM_PRHODJ(JIJ, JK) * PINV_PDZZ(JIJ, JK) *                                             &
                  ( ZWORKSHUGRAD_3D_4(JIJ, JK) )
        END DO
      END DO
      
      CALL DZF_PHY(D, ZWORKSHUGRAD_3D_8, ZWORKSHUGRAD_3D_9)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          PRSVS(JIJ, JK, JSV)=   PRSVS(JIJ, JK, JSV)                                          &
            -ZWORKSHUGRAD_3D_7(JIJ, JK)                                    &
            +ZWORKSHUGRAD_3D_9(JIJ, JK)
        END DO
      END DO
      
    ELSE
      CALL MXM_PHY(D, PRHODJ, ZWORKSHUGRAD_3D_1)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          ZWORKSHUGRAD_3D_2(JIJ, JK)=ZWORKSHUGRAD_3D_1(JIJ, JK) *ZFLXX(JIJ, JK) * PINV_PDXX(JIJ, JK)
        END DO
      END DO
      
      CALL DXF_PHY(D, ZWORKSHUGRAD_3D_2, ZWORKSHUGRAD_3D_3)
      
      DO JK=1, D%NKT
        DO JIJ=1, D%NIJT
          PRSVS(JIJ, JK, JSV)=   PRSVS(JIJ, JK, JSV)                                          &
            -ZWORKSHUGRAD_3D_3(JIJ, JK)
        END DO
      END DO
      
    END IF
  END IF
!
!
END DO    ! end loop JSV
!
!
END SUBROUTINE TURB_HOR_SV_FLUX
END MODULE MODE_TURB_HOR_SV_FLUX
