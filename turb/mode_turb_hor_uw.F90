!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_TURB_HOR_UW
IMPLICIT NONE
CONTAINS
!     ################################################################
      SUBROUTINE TURB_HOR_UW(TURBN,TLES,KSPLT,                       &
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
INTEGER,                  INTENT(IN)    ::  KRR          ! number of moist var.
INTEGER,                  INTENT(IN)    ::  KSV          ! number of sv var.
LOGICAL,                  INTENT(IN)    ::  OFLAT        ! Logical for zero ororography
TYPE(TFILEDATA),          INTENT(IN)    ::  TPFILE       ! Output file
!

REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PK          ! Turbulent diffusion doef.
                                                        ! PK = PLM * SQRT(PTKEM)
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PINV_PDXX   ! 1./PDXX
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PINV_PDZZ   ! 1./PDZZ
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PMZM_PRHODJ ! MZM(PRHODJ)
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PDXX, PDZZ, PDZX
                                                         ! Metric coefficients
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PRHODJ       ! density * grid volume
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PTHVREF      ! ref. state VPT       
!
! Variables at t-1
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PUM,PWM,PTHLM
REAL, DIMENSION(:,:,:,:), INTENT(IN)    ::  PRM
REAL, DIMENSION(:,:,:,:), INTENT(IN)    ::  PSVM
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PTKEM        ! TKE at time t- dt
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PLM          ! Turb. mixing length
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) ::  PRUS, PRWS
REAL, DIMENSION(:,:,:),   INTENT(INOUT) ::  PDP          ! TKE production terms
!
!
!
!
!*       0.2  declaration of local variables
!
REAL, DIMENSION(SIZE(PWM,1),SIZE(PWM,2),SIZE(PWM,3))       &
                                     :: ZFLX,ZWORK
    ! work arrays
!   
INTEGER             :: IKB,IKE,IKU
                                    ! Index values for the Beginning and End
                                    ! mass points of the domain  
INTEGER             :: JSV          ! scalar loop counter
!
REAL, DIMENSION(SIZE(PWM,1),SIZE(PWM,2),SIZE(PWM,3))  :: GX_W_UW_PWM
!
REAL :: ZTIME1, ZTIME2
TYPE(TFIELDMETADATA) :: TZFIELD
! ---------------------------------------------------------------------------
!
!*       1.   PRELIMINARY COMPUTATIONS
!             ------------------------
!
IKB = 1+JPVEXT               
IKE = SIZE(PWM,3)-JPVEXT    
IKU = SIZE(PWM,3)
!
!
GX_W_UW_PWM = GX_W_UW(PWM,PDXX,PDZZ,PDZX)
!
!
!*      13.   < U'W'>
!             -------
! 
! residual part of < U'W'> depending on dw/dx
!
ZFLX(:,:,:) =                                                      &
  - XCMFS * MXM(MZM(PK)) * GX_W_UW_PWM
!!         &  to be tested
!!  - (2./3.) * XCMFB * MZM( ZVPTU * MXM( PLM / SQRT(PTKEM) * XG / PTHVREF ) )
!
ZFLX(:,:,IKE+1) = 0.  ! rigid wall condition => no turbulent flux
!
! Nullify the flux at the ground level because it has been fully taken into
! account in turb_ver and extrapolate the flux under the ground 
ZFLX(:,:,IKB) = 0.
ZFLX(:,:,IKB-1)=2.*ZFLX(:,:,IKB)- ZFLX(:,:,IKB+1)
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
! taken into account in TURB_VER)
PRUS(:,:,:) = PRUS(:,:,:) - DZF( ZFLX* MXM( PMZM_PRHODJ ) / MXM( PDZZ ) )
!
!computation of the source for rho*W due to this flux
IF (.NOT. OFLAT) THEN
  PRWS(:,:,:) = PRWS(:,:,:)                              &
        -DXF( MZM( MXM(PRHODJ) * PINV_PDXX) * ZFLX)           &
        +DZM( PRHODJ * MXF( MZF( ZFLX*PDZX ) * PINV_PDXX ) / MZF(PDZZ) )
ELSE
  PRWS(:,:,:) = PRWS(:,:,:) -DXF( MZM( MXM(PRHODJ) * PINV_PDXX) * ZFLX)
END IF
! 
IF (KSPLT==1) THEN
  !
  !Contribution to the dynamic production of TKE:
  !
  ZWORK(:,:,:) =-MZF( MXF(                               &
     ZFLX *( GZ_U_UW(PUM,PDZZ) + GX_W_UW_PWM ) ) )
  !
  !
  ! evaluate the dynamic production at w(IKB+1) in PDP(IKB)
  ZWORK(:,:,IKB:IKB) = - MXF (                                             &
     ZFLX(:,:,IKB+1:IKB+1) *                                               &
   (   (PUM(:,:,IKB+1:IKB+1)-PUM(:,:,IKB:IKB)) / MXM(PDZZ(:,:,IKB+1:IKB+1))&
     + ( DXM( PWM(:,:,IKB+1:IKB+1) )                                       &
        -MXM(  (PWM(:,:,IKB+2:IKB+2)-PWM(:,:,IKB+1:IKB+1))                 &
                /(PDZZ(:,:,IKB+2:IKB+2)+PDZZ(:,:,IKB+1:IKB+1))             &
              +(PWM(:,:,IKB+1:IKB+1)-PWM(:,:,IKB  :IKB  ))                 &
                /(PDZZ(:,:,IKB+1:IKB+1)+PDZZ(:,:,IKB  :IKB  ))             &
            )                                                              &
          * PDZX(:,:,IKB+1:IKB+1)                                          &
       ) / (0.5*(PDXX(:,:,IKB+1:IKB+1)+PDXX(:,:,IKB:IKB)))                 &
   )                        )  
  !
  ! dynamic production computation
  PDP(:,:,:) = PDP(:,:,:) +  ZWORK(:,:,:)  
  !
END IF
!
! Storage in the LES configuration (addition to TURB_VER computation)
!
IF (TLES%LLES_CALL .AND. KSPLT==1) THEN
  CALL SECOND_MNH(ZTIME1)
  CALL LES_MEAN_SUBGRID( MZF(MXF(ZFLX)), TLES%X_LES_SUBGRID_WU , .TRUE. )
  CALL LES_MEAN_SUBGRID( MZF(MXF(GZ_U_UW(PUM,PDZZ)*ZFLX)), TLES%X_LES_RES_ddxa_U_SBG_UaU , .TRUE.)
  CALL LES_MEAN_SUBGRID( MZF(MXF(GX_W_UW_PWM*ZFLX)), TLES%X_LES_RES_ddxa_W_SBG_UaW , .TRUE.)
  CALL LES_MEAN_SUBGRID( MXF(GX_M_U(1,IKU,1,PTHLM,PDXX,PDZZ,PDZX)*MZF(ZFLX)),&
                         TLES%X_LES_RES_ddxa_Thl_SBG_UaW , .TRUE.)
  IF (KRR>=1) THEN
    CALL LES_MEAN_SUBGRID( MXF(GX_M_U(1,IKU,1,PRM(:,:,:,1),PDXX,PDZZ,PDZX)*MZF(ZFLX)), &
                           TLES%X_LES_RES_ddxa_Rt_SBG_UaW , .TRUE.)
  END IF
  DO JSV=1,KSV
    CALL LES_MEAN_SUBGRID( MXF(GX_M_U(1,IKU,1,PSVM(:,:,:,JSV),PDXX,PDZZ,PDZX)*MZF(ZFLX)), &
                           TLES%X_LES_RES_ddxa_Sv_SBG_UaW(:,:,:,JSV) , .TRUE.)
  END DO
  CALL SECOND_MNH(ZTIME2)
  TLES%XTIME_LES = TLES%XTIME_LES + ZTIME2 - ZTIME1
END IF

!
END SUBROUTINE TURB_HOR_UW
END MODULE MODE_TURB_HOR_UW
