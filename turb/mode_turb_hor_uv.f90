!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_TURB_HOR_UV
IMPLICIT NONE
CONTAINS
!     ################################################################
      SUBROUTINE TURB_HOR_UV(TURBN,TLES,KSPLT,OFLAT,O2D,             &
                      TPFILE,                                        &
                      PK,PINV_PDXX,PINV_PDYY,PINV_PDZZ,PMZM_PRHODJ,  &
                      PDXX,PDYY,PDZZ,PDZX,PDZY,                      &
                      PDIRCOSZW,                                     &
                      PCOSSLOPE,PSINSLOPE,                           &
                      PRHODJ,                                        &
                      PCDUEFF,PTAU11M,PTAU12M,PTAU22M,PTAU33M,       &
                      PUM,PVM,PUSLOPEM,PVSLOPEM,                     &
                      PDP,                                           &
                      PRUS,PRVS                                      )
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
!!                     Nov  06, 2002 (V. Masson) LES budgets
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
USE MODI_SHUMAN 
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
LOGICAL,                  INTENT(IN)    ::  OFLAT        ! Logical for zero ororography
LOGICAL,                  INTENT(IN)    ::  O2D          ! Logical for 2D model version (modd_conf)
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PK          ! Turbulent diffusion doef.
                                                        ! PK = PLM * SQRT(PTKEM)
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PINV_PDXX   ! 1./PDXX
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PINV_PDYY   ! 1./PDYY
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PINV_PDZZ   ! 1./PDZZ
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PMZM_PRHODJ ! MZM(PRHODJ)
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PDXX, PDYY, PDZZ, PDZX, PDZY 
                                                         ! Metric coefficients
REAL, DIMENSION(:,:),     INTENT(IN)    ::  PDIRCOSZW
! Director Cosinus along z directions at surface w-point
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PCOSSLOPE       ! cosinus of the angle 
                                      ! between i and the slope vector
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PSINSLOPE       ! sinus of the angle 
                                      ! between i and the slope vector
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PRHODJ       ! density * grid volume
!
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PCDUEFF      ! Cd * || u || at time t
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PTAU11M      ! <uu> in the axes linked 
       ! to the maximum slope direction and the surface normal and the binormal 
       ! at time t - dt
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PTAU12M      ! <uv> in the same axes
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PTAU22M      ! <vv> in the same axes
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PTAU33M      ! <ww> in the same axes
!
! Variables at t-1
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PUM,PVM
REAL, DIMENSION(:,:),      INTENT(IN)   ::  PUSLOPEM     ! wind component along the 
                                     ! maximum slope direction
REAL, DIMENSION(:,:),      INTENT(IN)   ::  PVSLOPEM     ! wind component along the 
                                     ! direction normal to the maximum slope one
!
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) ::  PRUS, PRVS   ! var. at t+1 -split-
REAL, DIMENSION(:,:,:),   INTENT(INOUT) ::  PDP          ! TKE production terms
!
!
!
!*       0.2  declaration of local variables
!
REAL, DIMENSION(SIZE(PUM,1),SIZE(PUM,2),SIZE(PUM,3))       &
                                     :: ZFLX,ZWORK
    ! work arrays
!   
REAL, DIMENSION(SIZE(PUM,1),SIZE(PUM,2)) ::ZDIRSINZW 
      ! sinus of the angle between the vertical and the normal to the orography
INTEGER             :: IKB,IKE
                                    ! Index values for the Beginning and End
                                    ! mass points of the domain  
!
REAL, DIMENSION(SIZE(PUM,1),SIZE(PUM,2),SIZE(PUM,3))  :: GY_U_UV_PUM
REAL, DIMENSION(SIZE(PVM,1),SIZE(PVM,2),SIZE(PVM,3))  :: GX_V_UV_PVM
!
REAL :: ZTIME1, ZTIME2
TYPE(TFIELDMETADATA) :: TZFIELD
! ---------------------------------------------------------------------------
!
!*       1.   PRELIMINARY COMPUTATIONS
!             ------------------------
!
IKB = 1+JPVEXT               
IKE = SIZE(PUM,3)-JPVEXT    
!
ZDIRSINZW(:,:) = SQRT( 1. - PDIRCOSZW(:,:)**2 )
!
GX_V_UV_PVM = GX_V_UV(PVM,PDXX,PDZZ,PDZX)
IF (.NOT. O2D) GY_U_UV_PUM = GY_U_UV(PUM,PDYY,PDZZ,PDZY)
!
!
!*      12.   < U'V'>
!             -------
!
!
IF (.NOT. O2D) THEN
  ZFLX(:,:,:)= - XCMFS * MYM(MXM(PK)) *                           &
          (GY_U_UV_PUM + GX_V_UV_PVM)
ELSE
  ZFLX(:,:,:)= - XCMFS * MYM(MXM(PK)) *                           &
          (GX_V_UV_PVM)
END IF
!
ZFLX(:,:,IKE+1)= ZFLX(:,:,IKE)
!
!
! Compute the correlation at the first physical level with the following 
! hypothesis du/dz vary in 1/z and w=0 at the ground
ZFLX(:,:,IKB:IKB)   = - XCMFS * MYM(MXM(PK(:,:,IKB:IKB))) *  (     &
  ( DYM( PUM(:,:,IKB:IKB) )                                        &
   -MYM( (PUM(:,:,IKB+1:IKB+1)-PUM(:,:,IKB:IKB))                   &
        *(1./MXM(PDZZ(:,:,IKB+1:IKB+1))+1./MXM(PDZZ(:,:,IKB:IKB))))&
    *0.5*MXM((PDZY(:,:,IKB+1:IKB+1)+PDZY(:,:,IKB:IKB)))            &
  ) / MXM(PDYY(:,:,IKB:IKB))                                       &
 +( DXM( PVM(:,:,IKB:IKB) )                                        &
   -MXM( (PVM(:,:,IKB+1:IKB+1)-PVM(:,:,IKB:IKB))                   &
        *(1./MYM(PDZZ(:,:,IKB+1:IKB+1))+1./MYM(PDZZ(:,:,IKB:IKB))))&
    *0.5*MYM((PDZX(:,:,IKB+1:IKB+1)+PDZX(:,:,IKB:IKB)))            &
  ) / MYM(PDXX(:,:,IKB:IKB))                                 ) 
! 
! extrapolates this flux under the ground with the surface flux
ZFLX(:,:,IKB-1) =                                                           &
   PTAU11M(:,:) * PCOSSLOPE(:,:) * PSINSLOPE(:,:) * PDIRCOSZW(:,:)**2         &
  +PTAU12M(:,:) * (PCOSSLOPE(:,:)**2 - PSINSLOPE(:,:)**2) *                   &
                  PDIRCOSZW(:,:)**2                                           &
  -PTAU22M(:,:) * PCOSSLOPE(:,:) * PSINSLOPE(:,:)                             &
  +PTAU33M(:,:) * PCOSSLOPE(:,:) * PSINSLOPE(:,:) * ZDIRSINZW(:,:)**2         &
  -PCDUEFF(:,:) * (                                                           &
    2. * PUSLOPEM(:,:) * PCOSSLOPE(:,:) * PSINSLOPE(:,:) *                    &
          PDIRCOSZW(:,:) * ZDIRSINZW(:,:)                                     &
    +PVSLOPEM(:,:) * (PCOSSLOPE(:,:)**2 - PSINSLOPE(:,:)**2) * ZDIRSINZW(:,:) &
                   )
!  
ZFLX(:,:,IKB-1:IKB-1) = 2. * MXM( MYM( ZFLX(:,:,IKB-1:IKB-1) ) )  &
                   - ZFLX(:,:,IKB:IKB)
!     
! stores  <U V>
IF ( TPFILE%LOPENED .AND. TURBN%LTURB_FLX ) THEN
  TZFIELD = TFIELDMETADATA(      &
    CMNHNAME   = 'UV_FLX',       &
    CSTDNAME   = '',             &
    CLONGNAME  = 'UV_FLX',       &
    CUNITS     = 'm2 s-2',       &
    CDIR       = 'XY',           &
    CCOMMENT   = 'X_Y_Z_UV_FLX', &
    NGRID      = 5,              &
    NTYPE      = TYPEREAL,       &
    NDIMS      = 3,              &
    LTIMEDEP   = .TRUE.          )
  CALL IO_FIELD_WRITE(TPFILE,TZFIELD,ZFLX)
END IF
!
!
!
!computation of the source for rho*V due to this flux
IF (.NOT. OFLAT) THEN
  PRUS(:,:,:) = PRUS(:,:,:)                                &
              - DYF(ZFLX * MXM(MYM(PRHODJ) * PINV_PDYY) )         &
              + DZF( MYF( MZM(ZFLX)*MXM(PDZY/MZM(PDYY)))   &
                    * MXM(PMZM_PRHODJ * PINV_PDZZ) )
ELSE
  PRUS(:,:,:) = PRUS(:,:,:) - DYF(ZFLX * MXM(MYM(PRHODJ) * PINV_PDYY) )
END IF
!
!computation of the source for rho*V due to this flux
IF (.NOT. OFLAT) THEN
  PRVS(:,:,:) = PRVS(:,:,:)                             &
                - DXF(ZFLX * MYM(MXM(PRHODJ) * PINV_PDXX) )    &
                + DZF( MXF( MZM(ZFLX)*MYM(PDZX/MZM(PDXX))) &
                      * MYM(PMZM_PRHODJ * PINV_PDZZ) )
ELSE
  PRVS(:,:,:) = PRVS(:,:,:) - DXF(ZFLX * MYM(MXM(PRHODJ) * PINV_PDXX) )
END IF
!
IF (KSPLT==1) THEN
  !
  !Contribution to the dynamic production of TKE:
  !
  IF (.NOT. O2D) THEN
    ZWORK(:,:,:) = - MXF( MYF( ZFLX *                                &
       (GY_U_UV_PUM + GX_V_UV_PVM) ) ) 
  ELSE
    ZWORK(:,:,:) = - MXF( MYF( ZFLX *                                &
       (GX_V_UV_PVM) ) )  
  ENDIF
  !
  ! evaluate the dynamic production at w(IKB+1) in PDP(IKB)
  !
  ZWORK(:,:,IKB:IKB) =  -                                             &
     MXF ( MYF( 0.5 * (ZFLX(:,:,IKB+1:IKB+1)+ZFLX(:,:,IKB:IKB)) ) )   &
   *(MXF ( MYF(                                                       &
       DYM( 0.5 * (PUM(:,:,IKB+1:IKB+1)+PUM(:,:,IKB:IKB))  )          &
      / MXM( 0.5*(PDYY(:,:,IKB:IKB)+PDYY(:,:,IKB+1:IKB+1)) )          &
      +DXM( 0.5 * (PVM(:,:,IKB+1:IKB+1)+PVM(:,:,IKB:IKB))  )          &
      / MYM( 0.5*(PDXX(:,:,IKB:IKB)+PDXX(:,:,IKB+1:IKB+1)) )          &
         )    )                                                       &  
    -MXF( (PUM(:,:,IKB+1:IKB+1)-PUM(:,:,IKB:IKB)) /                   &
              MXM(PDZZ(:,:,IKB+1:IKB+1)) * PDZY(:,:,IKB+1:IKB+1)      &
        ) / MYF(MXM( 0.5*(PDYY(:,:,IKB:IKB)+PDYY(:,:,IKB+1:IKB+1)) ) )&
    -MYF( (PVM(:,:,IKB+1:IKB+1)-PVM(:,:,IKB:IKB)) /                   &
              MYM(PDZZ(:,:,IKB+1:IKB+1)) * PDZX(:,:,IKB+1:IKB+1)      &
        ) / MXF(MYM( 0.5*(PDXX(:,:,IKB:IKB)+PDXX(:,:,IKB+1:IKB+1)) ) )&
    ) 
  !
  ! dynamic production 
  PDP(:,:,:) = PDP(:,:,:) + ZWORK(:,:,:)
  ! 
END IF
!
! Storage in the LES configuration
!
IF (TLES%LLES_CALL .AND. KSPLT==1) THEN
  CALL SECOND_MNH(ZTIME1)
  CALL LES_MEAN_SUBGRID( MXF(MYF(ZFLX)), TLES%X_LES_SUBGRID_UV ) 
  CALL LES_MEAN_SUBGRID( MXF(MYF(GY_U_UV(PUM,PDYY,PDZZ,PDZY)*ZFLX)), TLES%X_LES_RES_ddxa_U_SBG_UaU , .TRUE.)
  CALL LES_MEAN_SUBGRID( MXF(MYF(GX_V_UV(PVM,PDXX,PDZZ,PDZX)*ZFLX)), TLES%X_LES_RES_ddxa_V_SBG_UaV , .TRUE.)
  CALL SECOND_MNH(ZTIME2)
  TLES%XTIME_LES = TLES%XTIME_LES + ZTIME2 - ZTIME1
END IF
!
!
END SUBROUTINE TURB_HOR_UV
END MODULE MODE_TURB_HOR_UV
