!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_TURB_HOR_UV
IMPLICIT NONE
CONTAINS
!     ################################################################
      SUBROUTINE TURB_HOR_UV(D, TURBN,TLES,KSPLT,OFLAT,O2D,          &
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
TYPE(TFILEDATA),          INTENT(INOUT) ::  TPFILE       ! Output file
LOGICAL,                  INTENT(IN)    ::  OFLAT        ! Logical for zero ororography
LOGICAL,                  INTENT(IN)    ::  O2D          ! Logical for 2D model version (modd_conf)
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PK          ! Turbulent diffusion doef.
                                                        ! PK = PLM * SQRT(PTKEM)
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PINV_PDXX   ! 1./PDXX
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PINV_PDYY   ! 1./PDYY
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PINV_PDZZ   ! 1./PDZZ
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PMZM_PRHODJ ! MZM(PRHODJ)
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PDXX, PDYY, PDZZ, PDZX, PDZY 
                                                         ! Metric coefficients
REAL, DIMENSION(D%NIT,D%NJT),     INTENT(IN)    ::  PDIRCOSZW
! Director Cosinus along z directions at surface w-point
REAL, DIMENSION(D%NIT,D%NJT),   INTENT(IN)   ::  PCOSSLOPE       ! cosinus of the angle 
                                      ! between i and the slope vector
REAL, DIMENSION(D%NIT,D%NJT),   INTENT(IN)   ::  PSINSLOPE       ! sinus of the angle 
                                      ! between i and the slope vector
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PRHODJ       ! density * grid volume
!
REAL, DIMENSION(D%NIT,D%NJT),   INTENT(IN)   ::  PCDUEFF      ! Cd * || u || at time t
REAL, DIMENSION(D%NIT,D%NJT),   INTENT(IN)   ::  PTAU11M      ! <uu> in the axes linked 
       ! to the maximum slope direction and the surface normal and the binormal 
       ! at time t - dt
REAL, DIMENSION(D%NIT,D%NJT),   INTENT(IN)   ::  PTAU12M      ! <uv> in the same axes
REAL, DIMENSION(D%NIT,D%NJT),   INTENT(IN)   ::  PTAU22M      ! <vv> in the same axes
REAL, DIMENSION(D%NIT,D%NJT),   INTENT(IN)   ::  PTAU33M      ! <ww> in the same axes
!
! Variables at t-1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PUM,PVM
REAL, DIMENSION(D%NIT,D%NJT),      INTENT(IN)   ::  PUSLOPEM     ! wind component along the 
                                     ! maximum slope direction
REAL, DIMENSION(D%NIT,D%NJT),      INTENT(IN)   ::  PVSLOPEM     ! wind component along the 
                                     ! direction normal to the maximum slope one
!
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(INOUT) ::  PRUS, PRVS   ! var. at t+1 -split-
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(INOUT) ::  PDP          ! TKE production terms
!
!
!
!*       0.2  declaration of local variables
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: ZFLX,ZWORK ! work arrays
!   
REAL, DIMENSION(D%NIT,D%NJT) :: ZDIRSINZW 
      ! sinus of the angle between the vertical and the normal to the orography
INTEGER             :: IKB,IKE,IIT,IJT,IKT, JI,JJ,JK
                                    ! Index values for the Beginning and End
                                    ! mass points of the domain  
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT)  :: GY_U_UV_PUM
REAL, DIMENSION(D%NIT,D%NJT,D%NKT)  :: GX_V_UV_PVM
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
IIT=D%NIT
IJT=D%NJT
IKT=D%NKT  
!
!$acc kernels
!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
ZDIRSINZW(:,:) = SQRT( 1. - PDIRCOSZW(:,:)**2 )
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)
!$acc end kernels
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
!$acc kernels
!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
ZFLX(:,:,IKE+1)= ZFLX(:,:,IKE)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)
!$acc end kernels
!
!
! Compute the correlation at the first physical level with the following 
! hypothesis du/dz vary in 1/z and w=0 at the ground
ZFLX(:,:,IKB)   = - XCMFS * MYM(MXM(PK(:,:,IKB))) *  (     &
  ( DYM( PUM(:,:,IKB) )                                        &
   -MYM( (PUM(:,:,IKB+1)-PUM(:,:,IKB))                   &
        *(1./MXM(PDZZ(:,:,IKB+1))+1./MXM(PDZZ(:,:,IKB))))&
    *0.5*MXM((PDZY(:,:,IKB+1)+PDZY(:,:,IKB)))            &
  ) / MXM(PDYY(:,:,IKB))                                       &
 +( DXM( PVM(:,:,IKB) )                                        &
   -MXM( (PVM(:,:,IKB+1)-PVM(:,:,IKB))                   &
        *(1./MYM(PDZZ(:,:,IKB+1))+1./MYM(PDZZ(:,:,IKB))))&
    *0.5*MYM((PDZX(:,:,IKB+1)+PDZX(:,:,IKB)))            &
  ) / MYM(PDXX(:,:,IKB))                                 ) 
! 
! extrapolates this flux under the ground with the surface flux
!$acc kernels present_cr(ZFLX,ZDIRSINZW)
!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
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
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)
!$acc end kernels
!
ZFLX(:,:,IKB-1) = 2. * MXM( MYM( ZFLX(:,:,IKB-1) ) )  &
                   - ZFLX(:,:,IKB)
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
!$acc update self(ZFLX)
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
  ZWORK(:,:,IKB) =  -                                             &
     MXF ( MYF( 0.5 * (ZFLX(:,:,IKB+1)+ZFLX(:,:,IKB)) ) )   &
   *(MXF ( MYF(                                                       &
       DYM( 0.5 * (PUM(:,:,IKB+1)+PUM(:,:,IKB))  )          &
      / MXM( 0.5*(PDYY(:,:,IKB)+PDYY(:,:,IKB+1)) )          &
      +DXM( 0.5 * (PVM(:,:,IKB+1)+PVM(:,:,IKB))  )          &
      / MYM( 0.5*(PDXX(:,:,IKB)+PDXX(:,:,IKB+1)) )          &
         )    )                                                       &  
    -MXF( (PUM(:,:,IKB+1)-PUM(:,:,IKB)) /                   &
              MXM(PDZZ(:,:,IKB+1)) * PDZY(:,:,IKB+1)      &
        ) / MYF(MXM( 0.5*(PDYY(:,:,IKB)+PDYY(:,:,IKB+1)) ) )&
    -MYF( (PVM(:,:,IKB+1)-PVM(:,:,IKB)) /                   &
              MYM(PDZZ(:,:,IKB+1)) * PDZX(:,:,IKB+1)      &
        ) / MXF(MYM( 0.5*(PDXX(:,:,IKB)+PDXX(:,:,IKB+1)) ) )&
    ) 
  !
  ! dynamic production 
  !$acc kernels present_crm(PDP)
  !$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
  PDP(:,:,:) = PDP(:,:,:) + ZWORK(:,:,:)
  !$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
  !$acc end kernels
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
