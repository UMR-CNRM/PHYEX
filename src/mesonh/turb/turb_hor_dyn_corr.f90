!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODI_TURB_HOR_DYN_CORR
!
INTERFACE
!
      SUBROUTINE TURB_HOR_DYN_CORR(KSPLT, PTSTEP,                    &
                      OTURB_FLX,KRR,                                 &
                      TPFILE,                                        &
                      PK,PINV_PDZZ,                                  &
                      PDXX,PDYY,PDZZ,PDZX,PDZY,PZZ,                  &
                      PDIRCOSZW,                                     &
                      PCOSSLOPE,PSINSLOPE,                           &
                      PRHODJ,                                        &
                      PCDUEFF,PTAU11M,PTAU12M,PTAU22M,PTAU33M,       &
                      PUM,PVM,PWM,PUSLOPEM,PVSLOPEM,                 &
                      PTHLM,PRM,PSVM,                                &
                      PTKEM,PLM,                                     &
                      PDP,PTP,                                       &
                      PRUS,PRVS,PRWS                                 )
!
USE MODD_IO, ONLY: TFILEDATA
!
INTEGER,                  INTENT(IN)    ::  KSPLT        ! split process index
REAL,                     INTENT(IN)    ::  PTSTEP       ! timestep
LOGICAL,                  INTENT(IN)    ::  OTURB_FLX    ! switch to write the
                                 ! turbulent fluxes in the syncronous FM-file
INTEGER,                  INTENT(IN)    ::  KRR          ! number of moist var.
TYPE(TFILEDATA),          INTENT(IN)    ::  TPFILE       ! Output file
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PK          ! Turbulent diffusion doef.
                                                        ! PK = PLM * SQRT(PTKEM)
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PINV_PDZZ   ! 1./PDZZ
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PDXX, PDYY, PDZZ, PDZX, PDZY 
                                                         ! Metric coefficients
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PZZ          ! vertical grid
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
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PUM,PVM,PWM,PTHLM
REAL, DIMENSION(:,:,:,:), INTENT(IN)    ::  PRM
REAL, DIMENSION(:,:,:,:), INTENT(IN)    ::  PSVM
REAL, DIMENSION(:,:),      INTENT(IN)   ::  PUSLOPEM     ! wind component along the 
                                     ! maximum slope direction
REAL, DIMENSION(:,:),      INTENT(IN)   ::  PVSLOPEM     ! wind component along the 
                                     ! direction normal to the maximum slope one
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PTKEM        ! TKE at time t- dt
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PLM          ! Turb. mixing length
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) ::  PRUS, PRVS, PRWS
REAL, DIMENSION(:,:,:),   INTENT(INOUT) ::  PDP,PTP      ! TKE production terms
!
!
!
END SUBROUTINE TURB_HOR_DYN_CORR
!
END INTERFACE
!
END MODULE MODI_TURB_HOR_DYN_CORR
!     ################################################################
      SUBROUTINE TURB_HOR_DYN_CORR(KSPLT, PTSTEP,                    &
                      OTURB_FLX,KRR,                                 &
                      TPFILE,                                        &
                      PK,PINV_PDZZ,                                  &
                      PDXX,PDYY,PDZZ,PDZX,PDZY,PZZ,                  &
                      PDIRCOSZW,                                     &
                      PCOSSLOPE,PSINSLOPE,                           &
                      PRHODJ,                                        &
                      PCDUEFF,PTAU11M,PTAU12M,PTAU22M,PTAU33M,       &
                      PUM,PVM,PWM,PUSLOPEM,PVSLOPEM,                 &
                      PTHLM,PRM,PSVM,                                &
                      PTKEM,PLM,                                     &
                      PDP,PTP,                                       &
                      PRUS,PRVS,PRWS                                 )
!     ################################################################
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
!!                     Oct  18, 2000 (V. Masson) LES computations + LFLAT switch
!!                     Feb  15, 2001 (J. Stein)  remove the use of w=0 at the
!!                                               ground   
!!                     Mar  12, 2001 (V. Masson and J. Stein) major bugs 
!!                                 + change of discretization at the surface
!!                     Nov  06, 2002 (V. Masson) LES budgets
!!                     October 2009 (G. Tanguy) add ILENCH=LEN(YCOMMENT) after
!!                                              change of YCOMMENT
!!                     July 2012     (V.Masson) Implicitness of W
!!                     March 2014    (V.Masson) tridiag_w : bug between
!!                                               mass and flux position
!!                     J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!  P. Wautelet 20/05/2019: add name argument to ADDnFIELD_ll + new ADD4DFIELD_ll subroutine
!! --------------------------------------------------------------------------
!
!*      0. DECLARATIONS
!          ------------
!
USE MODD_ARGSLIST_ll,    ONLY: LIST_ll
USE MODD_CST
USE MODD_CONF
USE MODD_CTURB
use modd_field,          only: tfielddata, TYPEREAL
USE MODD_IO,             ONLY: TFILEDATA
USE MODD_PARAMETERS
USE MODD_LES
USE MODD_NSV
!
USE MODE_ll
USE MODE_IO_FIELD_WRITE, only: IO_Field_write
!
USE MODI_GRADIENT_M
USE MODI_GRADIENT_U
USE MODI_GRADIENT_V
USE MODI_GRADIENT_W
USE MODI_SHUMAN 
USE MODI_COEFJ
USE MODI_LES_MEAN_SUBGRID
USE MODI_TRIDIAG_W
!
USE MODI_SECOND_MNH
USE MODE_MPPDB
!
IMPLICIT NONE
!
!
!*       0.1  declaration of arguments
!
!
!
INTEGER,                  INTENT(IN)    ::  KSPLT        ! split process index
REAL,                     INTENT(IN)    ::  PTSTEP       ! timestep
LOGICAL,                  INTENT(IN)    ::  OTURB_FLX    ! switch to write the
                                 ! turbulent fluxes in the syncronous FM-file
INTEGER,                  INTENT(IN)    ::  KRR          ! number of moist var.
TYPE(TFILEDATA),          INTENT(IN)    ::  TPFILE       ! Output file
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PK          ! Turbulent diffusion doef.
                                                        ! PK = PLM * SQRT(PTKEM)
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PINV_PDZZ   ! 1./PDZZ
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PDXX, PDYY, PDZZ, PDZX, PDZY 
                                                         ! Metric coefficients
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PZZ          ! vertical grid
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
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PUM,PVM,PWM,PTHLM
REAL, DIMENSION(:,:,:,:), INTENT(IN)    ::  PRM
REAL, DIMENSION(:,:,:,:), INTENT(IN)    ::  PSVM
REAL, DIMENSION(:,:),      INTENT(IN)   ::  PUSLOPEM     ! wind component along the 
                                     ! maximum slope direction
REAL, DIMENSION(:,:),      INTENT(IN)   ::  PVSLOPEM     ! wind component along the 
                                     ! direction normal to the maximum slope one
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PTKEM        ! TKE at time t- dt
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PLM          ! Turb. mixing length
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) ::  PRUS, PRVS, PRWS
REAL, DIMENSION(:,:,:),   INTENT(INOUT) ::  PDP,PTP      ! TKE production terms
!
!
!
!*       0.2  declaration of local variables
!
REAL, DIMENSION(SIZE(PUM,1),SIZE(PUM,2),SIZE(PUM,3))       &
                                     :: ZFLX,ZWORK
    ! work arrays, PK is the turb. mixing coef.
!   
REAL, DIMENSION(SIZE(PUM,1),SIZE(PUM,2)) ::ZDIRSINZW 
      ! sinus of the angle between the vertical and the normal to the orography
INTEGER             :: IKB,IKE
                                    ! Index values for the Beginning and End
                                    ! mass points of the domain  
INTEGER             :: IKU                                   
INTEGER             :: JSV          ! scalar loop counter
!
REAL, DIMENSION(SIZE(PUM,1),SIZE(PUM,2),SIZE(PUM,3))  :: GX_U_M_PUM
REAL, DIMENSION(SIZE(PVM,1),SIZE(PVM,2),SIZE(PVM,3))  :: GY_V_M_PVM
REAL, DIMENSION(SIZE(PWM,1),SIZE(PWM,2),SIZE(PWM,3))  :: GZ_W_M_PWM
REAL, DIMENSION(SIZE(PWM,1),SIZE(PWM,2),SIZE(PWM,3))  :: GZ_W_M_ZWP
REAL, DIMENSION(SIZE(PWM,1),SIZE(PWM,2),SIZE(PWM,3))  :: ZMZF_DZZ   ! MZF(PDZZ)
REAL, DIMENSION(SIZE(PWM,1),SIZE(PWM,2),SIZE(PWM,3))  :: ZDFDDWDZ   ! formal derivative of the 
!                                                                   ! flux (variable: dW/dz)
REAL, DIMENSION(SIZE(PWM,1),SIZE(PWM,2),SIZE(PWM,3))  :: ZWP        ! W at future   time-step
!
REAL, DIMENSION(SIZE(PWM,1),SIZE(PWM,2),1) :: ZDU_DZ_DZS_DX ! du/dz*dzs/dx surf
REAL, DIMENSION(SIZE(PWM,1),SIZE(PWM,2),1) :: ZDV_DZ_DZS_DY ! dv/dz*dzs/dy surf
REAL, DIMENSION(SIZE(PWM,1),SIZE(PWM,2),1) :: ZDU_DX        ! du/dx        surf
REAL, DIMENSION(SIZE(PWM,1),SIZE(PWM,2),1) :: ZDV_DY        ! dv/dy        surf
REAL, DIMENSION(SIZE(PWM,1),SIZE(PWM,2),1) :: ZDW_DZ        ! dw/dz        surf
!
INTEGER                :: IINFO_ll      ! return code of parallel routine
TYPE(LIST_ll), POINTER :: TZFIELDS_ll   ! list of fields to exchange

REAL :: ZTIME1, ZTIME2


REAL, DIMENSION(SIZE(PDZZ,1),SIZE(PDZZ,2),1+JPVEXT:3+JPVEXT) :: ZCOEFF , ZDZZ
                                    ! coefficients for the uncentred gradient 
                                    ! computation near the ground
TYPE(TFIELDDATA) :: TZFIELD
! --------------------------------------------------------------------------
!
!*       1.   PRELIMINARY COMPUTATIONS
!             ------------------------
NULLIFY(TZFIELDS_ll)
!
IKB = 1+JPVEXT               
IKE = SIZE(PUM,3)-JPVEXT    
IKU = SIZE(PUM,3)
!
!
ZDIRSINZW(:,:) = SQRT( 1. - PDIRCOSZW(:,:)**2 )
!
GX_U_M_PUM = GX_U_M(PUM,PDXX,PDZZ,PDZX)
IF (.NOT. L2D) GY_V_M_PVM = GY_V_M(PVM,PDYY,PDZZ,PDZY)
GZ_W_M_PWM = GZ_W_M(PWM,PDZZ)
!
ZMZF_DZZ = MZF(PDZZ)
!
CALL ADD3DFIELD_ll( TZFIELDS_ll, ZFLX, 'TURB_HOR_DYN_CORR::ZFLX' )


!  compute the coefficients for the uncentred gradient computation near the 
!  ground
!
!*       9.   < U'U'>
!             -------
!
! Computes the U variance
IF (.NOT. L2D) THEN
  ZFLX(:,:,:)= (2./3.) * PTKEM                                  &
    - XCMFS * PK *( (4./3.) * GX_U_M_PUM                        &
                   -(2./3.) * ( GY_V_M_PVM                      &
                               +GZ_W_M_PWM                ) ) 
  !!  &   to be tested later
  !!  + XCMFB *  PLM / SQRT(PTKEM) * (-2./3.) * PTP 
ELSE
  ZFLX(:,:,:)= (2./3.) * PTKEM                                  &
    - XCMFS * PK *( (4./3.) * GX_U_M_PUM                        &
                   -(2./3.) * ( GZ_W_M_PWM                ) ) 
  !!  &   to be tested later
  !!  + XCMFB *  PLM / SQRT(PTKEM) * (-2./3.) * PTP 
END IF
!
ZFLX(:,:,IKE+1) = ZFLX(:,:,IKE) 
!
!* prescription of du/dz and dv/dz with uncentered gradient at the surface
!  prescription of dw/dz at Dz/2 above ground using the continuity equation
!  using a Boussinesq hypothesis to remove the z dependance of rhod_ref
!  (div u = 0)
!
ZDZZ(:,:,:) = MXM(PDZZ(:,:,IKB:IKB+2))
ZCOEFF(:,:,IKB+2)= - ZDZZ(:,:,2) /      &
       ( (ZDZZ(:,:,3)+ZDZZ(:,:,2)) * ZDZZ(:,:,3) )
ZCOEFF(:,:,IKB+1)=   (ZDZZ(:,:,3)+ZDZZ(:,:,2)) /      &
       ( ZDZZ(:,:,2) * ZDZZ(:,:,3) )
ZCOEFF(:,:,IKB)= - (ZDZZ(:,:,3)+2.*ZDZZ(:,:,2)) /      &
       ( (ZDZZ(:,:,3)+ZDZZ(:,:,2)) * ZDZZ(:,:,2) )
!
ZDU_DZ_DZS_DX(:,:,:)=MXF ((ZCOEFF(:,:,IKB+2:IKB+2)*PUM(:,:,IKB+2:IKB+2)       &
                          +ZCOEFF(:,:,IKB+1:IKB+1)*PUM(:,:,IKB+1:IKB+1)       &
                          +ZCOEFF(:,:,IKB  :IKB  )*PUM(:,:,IKB  :IKB  )       &
                          )* 0.5 * ( PDZX(:,:,IKB+1:IKB+1)+PDZX(:,:,IKB:IKB)) &
                         )/ MXF(PDXX(:,:,IKB:IKB))
!
ZDZZ(:,:,:) = MYM(PDZZ(:,:,IKB:IKB+2))
ZCOEFF(:,:,IKB+2)= - ZDZZ(:,:,2) /      &
       ( (ZDZZ(:,:,3)+ZDZZ(:,:,2)) * ZDZZ(:,:,3) )
ZCOEFF(:,:,IKB+1)=   (ZDZZ(:,:,3)+ZDZZ(:,:,2)) /      &
       ( ZDZZ(:,:,2) * ZDZZ(:,:,3) )
ZCOEFF(:,:,IKB)= - (ZDZZ(:,:,3)+2.*ZDZZ(:,:,2)) /      &
       ( (ZDZZ(:,:,3)+ZDZZ(:,:,2)) * ZDZZ(:,:,2) )
!

ZDV_DZ_DZS_DY(:,:,:)=MYF ((ZCOEFF(:,:,IKB+2:IKB+2)*PVM(:,:,IKB+2:IKB+2)       &
                          +ZCOEFF(:,:,IKB+1:IKB+1)*PVM(:,:,IKB+1:IKB+1)       &
                          +ZCOEFF(:,:,IKB  :IKB  )*PVM(:,:,IKB  :IKB  )       &
                          )* 0.5 * ( PDZY(:,:,IKB+1:IKB+1)+PDZY(:,:,IKB:IKB)) &
                         )/ MYF(PDYY(:,:,IKB:IKB))
!
!
ZDU_DX(:,:,:)=  DXF(PUM(:,:,IKB:IKB)) / MXF(PDXX(:,:,IKB:IKB))  &
              - ZDU_DZ_DZS_DX(:,:,:)

ZDV_DY(:,:,:)=  DYF(PVM(:,:,IKB:IKB)) / MYF(PDYY(:,:,IKB:IKB)) &
              - ZDV_DZ_DZS_DY(:,:,:)
!
ZDW_DZ(:,:,:)=-ZDU_DX(:,:,:)-ZDV_DY(:,:,:)
!
!* computation 
!
ZFLX(:,:,IKB)   = (2./3.) * PTKEM(:,:,IKB)                           &
  - XCMFS * PK(:,:,IKB) * 2. * ZDU_DX(:,:,1)


!!  &  to be tested later
!! + XCMFB * PLM(:,:,IKB:IKB) /SQRT(PTKEM(:,:,IKB:IKB)) *        &
!!   (-2./3.) * PTP(:,:,IKB:IKB)
!
! extrapolates this flux under the ground with the surface flux
ZFLX(:,:,IKB-1) =                                                            &
        PTAU11M(:,:) * PCOSSLOPE(:,:)**2 * PDIRCOSZW(:,:)**2                 &
  -2. * PTAU12M(:,:) * PCOSSLOPE(:,:)* PSINSLOPE(:,:) * PDIRCOSZW(:,:)       &
  +     PTAU22M(:,:) * PSINSLOPE(:,:)**2                                     &
  +     PTAU33M(:,:) * PCOSSLOPE(:,:)**2 * ZDIRSINZW(:,:)**2                 &
  +2. * PCDUEFF(:,:) *      (                                                &
      PVSLOPEM(:,:) * PCOSSLOPE(:,:)    * PSINSLOPE(:,:) * ZDIRSINZW(:,:)    &
    - PUSLOPEM(:,:) * PCOSSLOPE(:,:)**2 * ZDIRSINZW(:,:) * PDIRCOSZW(:,:)    )
! 
ZFLX(:,:,IKB-1) = 2. * ZFLX(:,:,IKB-1) -  ZFLX(:,:,IKB)
!
CALL UPDATE_HALO_ll(TZFIELDS_ll, IINFO_ll)
IF ( tpfile%lopened .AND. OTURB_FLX ) THEN
  ! stores <U U>  
  TZFIELD%CMNHNAME   = 'U_VAR'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'U_VAR'
  TZFIELD%CUNITS     = 'm2 s-2'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Z_U_VAR'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZFLX)
END IF
!
! Complete the U tendency
IF (.NOT. LFLAT) THEN
CALL MPPDB_CHECK3DM("before turb_corr:PRUS,PRHODJ,ZFLX,PDXX,PDZX,PINV_PDZZ",PRECISION,&
                   & PRUS,PRHODJ,ZFLX,PDXX,PDZX,PINV_PDZZ )

  PRUS(:,:,:)=PRUS                                            &
              -DXM(PRHODJ * ZFLX / MXF(PDXX) )                &
              +DZF( PDZX / MZM(PDXX) * MXM( MZM(PRHODJ*ZFLX) * PINV_PDZZ ) )
CALL MPPDB_CHECK3DM("after  turb_corr:PRUS,PRHODJ,ZFLX,PDXX,PDZX,PINV_PDZZ",PRECISION,&
                   & PRUS,PRHODJ,ZFLX,PDXX,PDZX,PINV_PDZZ )
ELSE
  PRUS(:,:,:)=PRUS -DXM(PRHODJ * ZFLX / MXF(PDXX) )
END IF
!
IF (KSPLT==1) THEN
  ! Contribution to the dynamic production of TKE:
  ZWORK(:,:,:)     = - ZFLX(:,:,:) * GX_U_M_PUM
  !
  ! evaluate the dynamic production at w(IKB+1) in PDP(IKB)
  !
  ZWORK(:,:,IKB) = 0.5* ( -ZFLX(:,:,IKB)*ZDU_DX(:,:,1) + ZWORK(:,:,IKB+1) )
  !
  PDP(:,:,:) = PDP(:,:,:) + ZWORK(:,:,:)
END IF
!
! Storage in the LES configuration
!
IF (LLES_CALL .AND. KSPLT==1) THEN
  CALL SECOND_MNH(ZTIME1)
  CALL LES_MEAN_SUBGRID( ZFLX, X_LES_SUBGRID_U2 ) 
  CALL LES_MEAN_SUBGRID( -ZWORK, X_LES_RES_ddxa_U_SBG_UaU , .TRUE.)
  CALL SECOND_MNH(ZTIME2)
  XTIME_LES = XTIME_LES + ZTIME2 - ZTIME1
END IF

!
!*      10.   < V'V'>
!             -------
!
! Computes the V variance
IF (.NOT. L2D) THEN
  ZFLX(:,:,:)= (2./3.) * PTKEM                                  &
    - XCMFS * PK *( (4./3.) * GY_V_M_PVM                        &
                   -(2./3.) * ( GX_U_M_PUM                      &
                               +GZ_W_M_PWM                ) )  
  !! &  to be tested
  !!  + XCMFB *  PLM / SQRT(PTKEM) * (-2./3.) * PTP 
  !
ELSE
  ZFLX(:,:,:)= (2./3.) * PTKEM                                  &
    - XCMFS * PK *(-(2./3.) * ( GX_U_M_PUM                      &
                               +GZ_W_M_PWM                ) )  
  !! &  to be tested
  !!  + XCMFB *  PLM / SQRT(PTKEM) * (-2./3.) * PTP 
  !
END IF
!
ZFLX(:,:,IKE+1) = ZFLX(:,:,IKE) 
!
ZFLX(:,:,IKB)   = (2./3.) * PTKEM(:,:,IKB)                           &
  - XCMFS * PK(:,:,IKB) * 2. * ZDV_DY(:,:,1)

!!           & to be tested
!! + XCMFB * PLM(:,:,IKB:IKB) /SQRT(PTKEM(:,:,IKB:IKB)) *         &
!!   (-2./3.) * PTP(:,:,IKB:IKB)
!
! extrapolates this flux under the ground with the surface flux
ZFLX(:,:,IKB-1) =                                                            &
        PTAU11M(:,:) * PSINSLOPE(:,:)**2 * PDIRCOSZW(:,:)**2                 &         
  +2. * PTAU12M(:,:) * PCOSSLOPE(:,:)* PSINSLOPE(:,:) * PDIRCOSZW(:,:)       &
  +     PTAU22M(:,:) * PCOSSLOPE(:,:)**2                                     &
  +     PTAU33M(:,:) * PSINSLOPE(:,:)**2 * ZDIRSINZW(:,:)**2                 &
  -2. * PCDUEFF(:,:)*       (                                                &
      PUSLOPEM(:,:) * PSINSLOPE(:,:)**2 * ZDIRSINZW(:,:) * PDIRCOSZW(:,:)    &
    + PVSLOPEM(:,:) * PCOSSLOPE(:,:)    * PSINSLOPE(:,:) * ZDIRSINZW(:,:)    )
! 
ZFLX(:,:,IKB-1) = 2. * ZFLX(:,:,IKB-1) -  ZFLX(:,:,IKB)
!
CALL UPDATE_HALO_ll(TZFIELDS_ll, IINFO_ll)
!
IF ( tpfile%lopened .AND. OTURB_FLX ) THEN
  ! stores <V V>  
  TZFIELD%CMNHNAME   = 'V_VAR'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'V_VAR'
  TZFIELD%CUNITS     = 'm2 s-2'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Z_V_VAR'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZFLX)
END IF
!
! Complete the V tendency
IF (.NOT. L2D) THEN
  IF (.NOT. LFLAT) THEN
    PRVS(:,:,:)=PRVS                                          &
                -DYM(PRHODJ * ZFLX / MYF(PDYY) )              &
                +DZF( PDZY / MZM(PDYY) *                      &
                MYM( MZM(PRHODJ*ZFLX) * PINV_PDZZ ) )
  ELSE
    PRVS(:,:,:)=PRVS -DYM(PRHODJ * ZFLX / MYF(PDYY) )
  END IF
!
! Contribution to the dynamic production of TKE:
  IF (KSPLT==1) ZWORK(:,:,:)     = - ZFLX(:,:,:) * GY_V_M_PVM
ELSE
  ZWORK(:,:,:)     = 0.
END IF
!
IF (KSPLT==1) THEN
  !
  ! evaluate the dynamic production at w(IKB+1) in PDP(IKB)
  !
  ZWORK(:,:,IKB) = 0.5* ( -ZFLX(:,:,IKB)*ZDV_DY(:,:,1) + ZWORK(:,:,IKB+1) )
  !
  PDP(:,:,:) = PDP(:,:,:) + ZWORK(:,:,:)
END IF
!
! Storage in the LES configuration
!
IF (LLES_CALL .AND. KSPLT==1) THEN
  CALL SECOND_MNH(ZTIME1)
  CALL LES_MEAN_SUBGRID( ZFLX, X_LES_SUBGRID_V2 ) 
  CALL LES_MEAN_SUBGRID( -ZWORK, X_LES_RES_ddxa_V_SBG_UaV , .TRUE.)
  CALL SECOND_MNH(ZTIME2)
  XTIME_LES = XTIME_LES + ZTIME2 - ZTIME1
END IF
!
!*      11.   < W'W'>
!             -------
!
! Computes the W variance
IF (.NOT. L2D) THEN
  ZFLX(:,:,:)= (2./3.) * PTKEM                                  &
    - XCMFS * PK *( (4./3.) * GZ_W_M_PWM                        &
                   -(2./3.) * ( GX_U_M_PUM                      &
                               +GY_V_M_PVM                ) ) 
  !!  &  to be tested
  !!    -2.* XCMFB *  PLM / SQRT(PTKEM) * (-2./3.) * PTP 
ELSE
  ZFLX(:,:,:)= (2./3.) * PTKEM                                  &
    - XCMFS * PK *( (4./3.) * GZ_W_M_PWM                        &
                   -(2./3.) * ( GX_U_M_PUM                ) ) 
  !!  &  to be tested
  !!    -2.* XCMFB *  PLM / SQRT(PTKEM) * (-2./3.) * PTP 
END IF
!
ZFLX(:,:,IKE+1)= ZFLX(:,:,IKE)
!
ZFLX(:,:,IKB)   = (2./3.) * PTKEM(:,:,IKB)                           &
  - XCMFS * PK(:,:,IKB) * 2. * ZDW_DZ(:,:,1)

!             &  to be tested
!   - 2.* XCMFB * PLM(:,:,IKB:IKB) /SQRT(PTKEM(:,:,IKB:IKB)) *             &
!  (-2./3.) * PTP(:,:,IKB:IKB)
!
! extrapolates this flux under the ground with the surface flux
ZFLX(:,:,IKB-1) =                                                     &
        PTAU11M(:,:) * ZDIRSINZW(:,:)**2                                &
  +     PTAU33M(:,:) * PDIRCOSZW(:,:)**2                                &
  +2. * PCDUEFF(:,:)* PUSLOPEM(:,:)  * ZDIRSINZW(:,:) * PDIRCOSZW(:,:) 
  ! 
ZFLX(:,:,IKB-1) = 2. * ZFLX(:,:,IKB-1) - ZFLX(:,:,IKB)
!
IF ( tpfile%lopened .AND. OTURB_FLX ) THEN
  ! stores <W W>  
  TZFIELD%CMNHNAME   = 'W_VAR'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'W_VAR'
  TZFIELD%CUNITS     = 'm2 s-2'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Z_W_VAR'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZFLX)
END IF
!
! Complete the W tendency
!
!PRWS(:,:,:)=PRWS(:,:,:) - DZM( PRHODJ*ZFLX/MZF(PDZZ) )
ZDFDDWDZ(:,:,:)    = - XCMFS * PK(:,:,:) * (4./3.)
ZDFDDWDZ(:,:,:IKB) = 0.
!
CALL TRIDIAG_W(PWM,ZFLX,ZDFDDWDZ,PTSTEP,ZMZF_DZZ,PRHODJ,ZWP)
!
PRWS = PRWS(:,:,:) + MZM(PRHODJ(:,:,:))*(ZWP(:,:,:)-PWM(:,:,:))/PTSTEP
!
!* recomputes flux using guess of W
!
GZ_W_M_ZWP = GZ_W_M(ZWP,PDZZ)
ZFLX(:,:,IKB+1:)=ZFLX(:,:,IKB+1:) &
  - XCMFS * PK(:,:,IKB+1:) * (4./3.) * (GZ_W_M_ZWP(:,:,IKB+1:) - GZ_W_M_PWM(:,:,IKB+1:))
!
IF (KSPLT==1) THEN
  !Contribution to the dynamic production of TKE:
! ZWORK(:,:,:) = - ZFLX(:,:,:) * GZ_W_M_PWM
  ZWORK(:,:,:) = - ZFLX(:,:,:) * GZ_W_M_ZWP
  !
  ! evaluate the dynamic production at w(IKB+1) in PDP(IKB)
  !
  ZWORK(:,:,IKB) = 0.5* ( -ZFLX(:,:,IKB)*ZDW_DZ(:,:,1) + ZWORK(:,:,IKB+1) )
  !
  PDP(:,:,:) = PDP(:,:,:) + ZWORK(:,:,:)
END IF
!
! Storage in the LES configuration
!
!
IF (LLES_CALL .AND. KSPLT==1) THEN
  CALL SECOND_MNH(ZTIME1)
  CALL LES_MEAN_SUBGRID( ZFLX, X_LES_SUBGRID_W2 ) 
  CALL LES_MEAN_SUBGRID( -ZWORK, X_LES_RES_ddxa_W_SBG_UaW , .TRUE.)
  CALL LES_MEAN_SUBGRID( GZ_M_M(PTHLM,PDZZ)*ZFLX, X_LES_RES_ddxa_Thl_SBG_UaW , .TRUE.)
  CALL LES_MEAN_SUBGRID(ZFLX*MZF(GZ_M_W(1,IKU,1,PTHLM,PDZZ)),X_LES_RES_ddz_Thl_SBG_W2)
  IF (KRR>=1) THEN
    CALL LES_MEAN_SUBGRID( GZ_M_M(PRM(:,:,:,1),PDZZ)*ZFLX, &
                           X_LES_RES_ddxa_Rt_SBG_UaW , .TRUE.)
    CALL LES_MEAN_SUBGRID(ZFLX*MZF(GZ_M_W(1,IKU,1,PRM(:,:,:,1),PDZZ)), &
                           X_LES_RES_ddz_Rt_SBG_W2)
  END IF
  DO JSV=1,NSV
    CALL LES_MEAN_SUBGRID( GZ_M_M(PSVM(:,:,:,JSV),PDZZ)*ZFLX, &
                           X_LES_RES_ddxa_Sv_SBG_UaW(:,:,:,JSV) , .TRUE.)
    CALL LES_MEAN_SUBGRID(ZFLX*MZF(GZ_M_W(1,IKU,1,PSVM(:,:,:,JSV),PDZZ)), &
                           X_LES_RES_ddz_Sv_SBG_W2(:,:,:,JSV))
  END DO
  CALL SECOND_MNH(ZTIME2)
  XTIME_LES = XTIME_LES + ZTIME2 - ZTIME1
END IF
!
CALL CLEANLIST_ll(TZFIELDS_ll)
!
!
END SUBROUTINE TURB_HOR_DYN_CORR
