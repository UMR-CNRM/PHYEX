!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_TURB_HOR_DYN_CORR
IMPLICIT NONE
CONTAINS
      SUBROUTINE TURB_HOR_DYN_CORR(D,TURBN,TLES,KSPLT, PTSTEP,       &
                      KRR, KSV,OFLAT,O2D,                            &
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
!!                     Oct  18, 2000 (V. Masson) LES computations + OFLAT switch
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
!!                      M.Moge  04/2016 Use openACC directives to port the TURB part of Meso-NH on GPU
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!  P. Wautelet 20/05/2019: add name argument to ADDnFIELD_ll + new ADD4DFIELD_ll subroutine
!  J.Escobar 13/08/2020: PGI/NVHPC BUG , extend DO CONCURRENT to 3D indexes        
!! --------------------------------------------------------------------------
!
!*      0. DECLARATIONS
!          ------------
!
USE MODD_TURB_n, ONLY: TURB_t
USE MODD_DIMPHYEX,   ONLY: DIMPHYEX_t
!
USE MODD_ARGSLIST_ll,    ONLY: LIST_ll
USE MODD_CST
USE MODD_CTURB
use modd_field,          only: tfieldmetadata, TYPEREAL
USE MODD_IO,             ONLY: TFILEDATA
USE MODD_PARAMETERS
USE MODD_LES, ONLY: TLES_t
!
USE MODE_ll
USE MODE_IO_FIELD_WRITE, ONLY: IO_FIELD_WRITE
!
USE MODI_GRADIENT_M
USE MODI_GRADIENT_U
USE MODI_GRADIENT_V
USE MODI_GRADIENT_W
USE MODI_SHUMAN
USE MODE_SHUMAN_PHY
USE MODE_GRADIENT_M_PHY, ONLY: GZ_M_W_PHY
USE MODE_GRADIENT_W_PHY, ONLY: GZ_W_M_PHY, GZ_W_M_PHY
USE MODE_GRADIENT_V_PHY, ONLY: GY_V_M_PHY
USE MODE_GRADIENT_U_PHY, ONLY: GX_U_M_PHY
USE MODE_COEFJ, ONLY: COEFJ
#ifdef MNH_OPENACC
USE MODI_GET_HALO
#endif
USE MODI_LES_MEAN_SUBGRID
USE MODE_TRIDIAG_W, ONLY: TRIDIAG_W
!
USE MODI_SECOND_MNH
USE MODE_MPPDB
!
! These macro are handled by pft_tool.py --craybyPassDOCONCURRENT applied on Cray Rules
#ifdef MNH_COMPILER_CCE
!$mnh_undef(LOOP)
!$mnh_undef(OPENACC)
#endif
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
REAL,                     INTENT(IN)    ::  PTSTEP       ! timestep
INTEGER,                  INTENT(IN)    ::  KRR          ! number of moist var.
INTEGER,                  INTENT(IN)    ::  KSV          ! number of sv var.
LOGICAL,                  INTENT(IN)    ::  OFLAT        ! Logical for zero ororography
LOGICAL,                  INTENT(IN)    ::  O2D          ! Logical for 2D model version (modd_conf)
TYPE(TFILEDATA),          INTENT(INOUT)    ::  TPFILE       ! Output file
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PK          ! Turbulent diffusion doef.
                                                        ! PK = PLM * SQRT(PTKEM)
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PINV_PDZZ   ! 1./PDZZ
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PDXX, PDYY, PDZZ, PDZX, PDZY 
                                                         ! Metric coefficients
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PZZ          ! vertical grid
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
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PUM,PVM,PWM,PTHLM
REAL, DIMENSION(D%NIT,D%NJT,D%NKT,KRR), INTENT(IN)    ::  PRM
REAL, DIMENSION(D%NIT,D%NJT,D%NKT,KSV), INTENT(IN)    ::  PSVM
REAL, DIMENSION(D%NIT,D%NJT),      INTENT(IN)   ::  PUSLOPEM     ! wind component along the 
                                     ! maximum slope direction
REAL, DIMENSION(D%NIT,D%NJT),      INTENT(IN)   ::  PVSLOPEM     ! wind component along the 
                                     ! direction normal to the maximum slope one
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PTKEM        ! TKE at time t- dt
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PLM          ! Turb. mixing length
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(INOUT) ::  PRUS, PRVS, PRWS
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(INOUT) ::  PDP,PTP      ! TKE production terms
!
!
!
!*       0.2  declaration of local variables
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: ZFLX,ZWORK
    ! work arrays, PK is the turb. mixing coef.
!   
REAL, DIMENSION(D%NIT,D%NJT) ::ZDIRSINZW 
      ! sinus of the angle between the vertical and the normal to the orography
INTEGER             :: IKB,IKE
                                    ! Index values for the Beginning and End
                                    ! mass points of the domain  
INTEGER             :: IKU,IKT,IIT,IJT                                   
INTEGER             :: JSV,JI,JJ,JK          ! scalar loop counter
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT)  :: GX_U_M_PUM
REAL, DIMENSION(D%NIT,D%NJT,D%NKT)  :: GY_V_M_PVM
REAL, DIMENSION(D%NIT,D%NJT,D%NKT)  :: GZ_W_M_PWM
REAL, DIMENSION(D%NIT,D%NJT,D%NKT)  :: GZ_W_M_ZWP
REAL, DIMENSION(D%NIT,D%NJT,D%NKT)  :: ZMZF_DZZ   ! MZF(PDZZ)
REAL, DIMENSION(D%NIT,D%NJT,D%NKT)  :: ZDFDDWDZ   ! formal derivative of the 
!                                                                   ! flux (variable: dW/dz)
REAL, DIMENSION(D%NIT,D%NJT,D%NKT)  :: ZWP        ! W at future   time-step
!
REAL, DIMENSION(D%NIT,D%NJT) :: ZDU_DZ_DZS_DX ! du/dz*dzs/dx surf
REAL, DIMENSION(D%NIT,D%NJT) :: ZDV_DZ_DZS_DY ! dv/dz*dzs/dy surf
REAL, DIMENSION(D%NIT,D%NJT) :: ZDU_DX        ! du/dx        surf
REAL, DIMENSION(D%NIT,D%NJT) :: ZDV_DY        ! dv/dy        surf
REAL, DIMENSION(D%NIT,D%NJT) :: ZDW_DZ        ! dw/dz        surf
!
INTEGER                :: IINFO_ll      ! return code of parallel routine
TYPE(LIST_ll), POINTER :: TZFIELDS_ll   ! list of fields to exchange

REAL :: ZTIME1, ZTIME2


REAL, DIMENSION(D%NIT,D%NJT,1+JPVEXT:3+JPVEXT) :: ZCOEFF , ZDZZ
                                    ! coefficients for the uncentred gradient 
                                    ! computation near the ground
TYPE(TFIELDMETADATA) :: TZFIELD
! --------------------------------------------------------------------------
!
!*       1.   PRELIMINARY COMPUTATIONS
!             ------------------------
NULLIFY(TZFIELDS_ll)
!
IKB = 1+JPVEXT               
IKE = SIZE(PUM,3)-JPVEXT    
IKU = SIZE(PUM,3)
IIT=D%NIT
IJT=D%NJT
IKT=D%NKT
!
!
!$acc kernels ! async(1)
!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
ZDIRSINZW(:,:) = SQRT( 1. - PDIRCOSZW(:,:)**2 )
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)
!$acc end kernels
!
GX_U_M_PUM = GX_U_M(OFLAT,PUM,PDXX,PDZZ,PDZX)
IF (.NOT. O2D) GY_V_M_PVM = GY_V_M(OFLAT,PVM,PDYY,PDZZ,PDZY)
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
IF (.NOT. O2D) THEN
   !$acc kernels present_cr(zflx,gz_w_m_pwm) ! async(2)
   !$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
      ZFLX(:,:,:)= (2./3.) * PTKEM(:,:,:)                            &
           - XCMFS * PK(:,:,:) *( (4./3.) * GX_U_M_PUM(:,:,:)        &
           -(2./3.) * ( GY_V_M_PVM(:,:,:)                     &
           +GZ_W_M_PWM(:,:,:)                ) )
   !$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
   !$acc end kernels
   !!  &   to be tested later
  !!  + XCMFB *  PLM / SQRT(PTKEM) * (-2./3.) * PTP 
ELSE
  !$acc kernels ! async(2)
  !$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
  ZFLX(:,:,:)= (2./3.) * PTKEM(:,:,:)                                  &
    - XCMFS * PK(:,:,:) *( (4./3.) * GX_U_M_PUM(:,:,:)                 &
                   -(2./3.) * ( GZ_W_M_PWM(:,:,:)             ) ) 
  !$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
  !$acc end kernels
  !!  &   to be tested later
  !!  + XCMFB *  PLM / SQRT(PTKEM) * (-2./3.) * PTP 
END IF
!
!$acc kernels ! async(2)
!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
   ZFLX(:,:,IKE+1) = ZFLX(:,:,IKE) 
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)
!$acc end kernels
!
!* prescription of du/dz and dv/dz with uncentered gradient at the surface
!  prescription of dw/dz at Dz/2 above ground using the continuity equation
!  using a Boussinesq hypothesis to remove the z dependance of rhod_ref
!  (div u = 0)
!
! ZDZZ(:,:,:) = MXM(PDZZ(:,:,IKB:IKB+2)) ! case not handled yet by pyft --expandAllArraysConcurrent and shumanFUNCtoCALL
ZDZZ(:,:,IKB) = MXM(PDZZ(:,:,IKB))
ZDZZ(:,:,IKB+1) = MXM(PDZZ(:,:,IKB+1))
ZDZZ(:,:,IKB+2) = MXM(PDZZ(:,:,IKB+2))
!$acc kernels present_cr(zdzz,zcoeff) ! async(3)
!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
   ZCOEFF(:,:,IKB+2)= - ZDZZ(:,:,2) /      &
        ( (ZDZZ(:,:,3)+ZDZZ(:,:,2)) * ZDZZ(:,:,3) )
   ZCOEFF(:,:,IKB+1)=   (ZDZZ(:,:,3)+ZDZZ(:,:,2)) /      &
        ( ZDZZ(:,:,2) * ZDZZ(:,:,3) )
   ZCOEFF(:,:,IKB)= - (ZDZZ(:,:,3)+2.*ZDZZ(:,:,2)) /      &
        ( (ZDZZ(:,:,3)+ZDZZ(:,:,2)) * ZDZZ(:,:,2) )
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)
!$acc end kernels
!
ZDU_DZ_DZS_DX(:,:)=MXF ((ZCOEFF(:,:,IKB+2)*PUM(:,:,IKB+2)       &
                          +ZCOEFF(:,:,IKB+1)*PUM(:,:,IKB+1)       &
                          +ZCOEFF(:,:,IKB)*PUM(:,:,IKB)       &
                          )* 0.5 * ( PDZX(:,:,IKB+1)+PDZX(:,:,IKB)) &
                         )/ MXF(PDXX(:,:,IKB))
!
ZDZZ(:,:,IKB) = MYM(PDZZ(:,:,IKB))
ZDZZ(:,:,IKB+1) = MYM(PDZZ(:,:,IKB+1))
ZDZZ(:,:,IKB+2) = MYM(PDZZ(:,:,IKB+2))
!
!$acc kernels present_cr(zdzz,zcoeff) ! async(4)
!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
   ZCOEFF(:,:,IKB+2)= - ZDZZ(:,:,2) /      &
        ( (ZDZZ(:,:,3)+ZDZZ(:,:,2)) * ZDZZ(:,:,3) )
   ZCOEFF(:,:,IKB+1)=   (ZDZZ(:,:,3)+ZDZZ(:,:,2)) /      &
        ( ZDZZ(:,:,2) * ZDZZ(:,:,3) )
   ZCOEFF(:,:,IKB)= - (ZDZZ(:,:,3)+2.*ZDZZ(:,:,2)) /      &
        ( (ZDZZ(:,:,3)+ZDZZ(:,:,2)) * ZDZZ(:,:,2) )
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)
!$acc end kernels

ZDV_DZ_DZS_DY(:,:)=MYF ((ZCOEFF(:,:,IKB+2)*PVM(:,:,IKB+2)       &
                          +ZCOEFF(:,:,IKB+1)*PVM(:,:,IKB+1)       &
                          +ZCOEFF(:,:,IKB)*PVM(:,:,IKB)       &
                          )* 0.5 * ( PDZY(:,:,IKB+1)+PDZY(:,:,IKB)) &
                         )/ MYF(PDYY(:,:,IKB))
!
!
ZDU_DX(:,:)=  DXF(PUM(:,:,IKB)) / MXF(PDXX(:,:,IKB))  &
              - ZDU_DZ_DZS_DX(:,:)

ZDV_DY(:,:)=  DYF(PVM(:,:,IKB)) / MYF(PDYY(:,:,IKB)) &
              - ZDV_DZ_DZS_DY(:,:)
!
!$acc kernels present_cr(zdv_dy,zdw_dz) ! async(4)
!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
ZDW_DZ(:,:)=-ZDU_DX(:,:)-ZDV_DY(:,:)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)
!$acc end kernels
!
!* computation 
!
!!! wait for the computation of ZFLX
!!$acc ! wait(2) async(4)
!!! wait for the computation of ZDW_DZ
!!$acc ! wait(4)
!
! ! !!! we can launch the update of ZFLX on the part that has already been computed
! ! !$acc update self(ZFLX(:,:,IKB+1:)) async(10)
!attention !!!!! je ne comprends pas pourquoi mais ce update plante à l'execution...
! du coup je ne peux pas faire de update self asynchrone...
!
!$acc kernels present_cr(zdu_dx,zflx) ! async(3)
!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
   ZFLX(:,:,IKB)   = (2./3.) * PTKEM(:,:,IKB)                           &
        - XCMFS * PK(:,:,IKB) * 2. * ZDU_DX(:,:)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)
!$acc end kernels

!!  &  to be tested later
!! + XCMFB * PLM(:,:,IKB:IKB) /SQRT(PTKEM(:,:,IKB:IKB)) *        &
!!   (-2./3.) * PTP(:,:,IKB:IKB)
!
! extrapolates this flux under the ground with the surface flux
!
!
!!! wait for the computation of ZDIRSINZW
!!$acc ! wait(1)
!
!$acc kernels present_cr(ZFLX,ZDIRSINZW) ! async(4)
!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
ZFLX(:,:,IKB-1) =                                                            &
        PTAU11M(:,:) * PCOSSLOPE(:,:)**2 * PDIRCOSZW(:,:)**2                 &
  -2. * PTAU12M(:,:) * PCOSSLOPE(:,:)* PSINSLOPE(:,:) * PDIRCOSZW(:,:)       &
  +     PTAU22M(:,:) * PSINSLOPE(:,:)**2                                     &
  +     PTAU33M(:,:) * PCOSSLOPE(:,:)**2 * ZDIRSINZW(:,:)**2                 &
  +2. * PCDUEFF(:,:) *      (                                                &
      PVSLOPEM(:,:) * PCOSSLOPE(:,:)    * PSINSLOPE(:,:) * ZDIRSINZW(:,:)    &
    - PUSLOPEM(:,:) * PCOSSLOPE(:,:)**2 * ZDIRSINZW(:,:) * PDIRCOSZW(:,:)    )
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)
!$acc end kernels
! 
!!! wait for the computation of ZFLX(:,:,IKB) and ZFLX(:,:,IKB-1)
!!$acc ! wait(3) async(4)
!
!$acc kernels ! async(4)
!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
   ZFLX(:,:,IKB-1) = 2. * ZFLX(:,:,IKB-1) -  ZFLX(:,:,IKB)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)
!$acc end kernels
!
!
!!! wait for the computation of ZFLX(:,:,IKB-1)
!!$acc ! wait(4)
!


! ! !!! we can launch the update of ZFLX on the rest
! ! !$acc update self(ZFLX(:,:,1:IKB)) async(11)
! ! !
! ! !!! and wait for the update self(ZFLX(...)) to complete
! ! !$acc wait(10)
! ! !$acc wait(11)
!attention !!!!! je ne comprends pas pourquoi mais le update self(ZFLX(:,:,IKB+1:)) plante à l'execution...
! du coup je ne peux pas faire de update self asynchrone...


!
!!! at this point there are no more async operations running
!!! to be absolutely sure, we do a wait
!!$acc ! wait
!
#ifndef MNH_OPENACC
CALL UPDATE_HALO_ll(TZFIELDS_ll, IINFO_ll)
#else
CALL GET_HALO_D(ZFLX,HNAME='TURB_HOR_DYN_CORR::ZFLX')
#endif
!
IF ( TPFILE%LOPENED .AND. TURBN%LTURB_FLX ) THEN
  !$acc update self( ZFLX )
  ! stores <U U>  
  TZFIELD = TFIELDMETADATA(     &
    CMNHNAME   = 'U_VAR',       &
    CSTDNAME   = '',            &
    CLONGNAME  = 'U_VAR',       &
    CUNITS     = 'm2 s-2',      &
    CDIR       = 'XY',          &
    CCOMMENT   = 'X_Y_Z_U_VAR', &
    NGRID      = 1,             &
    NTYPE      = TYPEREAL,      &
    NDIMS      = 3,             &
    LTIMEDEP   = .TRUE.         )
  CALL IO_Field_write(TPFILE,TZFIELD,ZFLX)
END IF
!
! Complete the U tendency
IF (.NOT. OFLAT) THEN
  PRUS(:,:,:)=PRUS                                            &
              -DXM(PRHODJ * ZFLX / MXF(PDXX) )                &
              +DZF( PDZX / MZM(PDXX) * MXM( MZM(PRHODJ*ZFLX) * PINV_PDZZ ) )
ELSE
  PRUS(:,:,:)=PRUS -DXM(PRHODJ * ZFLX / MXF(PDXX) )
END IF
!
IF (KSPLT==1) THEN
  ! Contribution to the dynamic production of TKE:
   !$acc kernels present_cr(gx_u_m_pum,zwork) ! async(2)
   !$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
      ZWORK(:,:,:)     = - ZFLX(:,:,:) * GX_U_M_PUM(:,:,:)
   !$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
  !$acc end kernels
  !
  ! evaluate the dynamic production at w(IKB+1) in PDP(IKB)
  !
  !$acc kernels present_cr(zdu_dx,zwork) ! async(2)
   !$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
      ZWORK(:,:,IKB) = 0.5* ( -ZFLX(:,:,IKB)*ZDU_DX(:,:) + ZWORK(:,:,IKB+1) )
   !$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)
   !$acc end kernels
  !
  !$acc kernels present_cr(ZWORK) ! async(2)
  PDP(:,:,:) = PDP(:,:,:) + ZWORK(:,:,:)
  !$acc end kernels
END IF
!
! Storage in the LES configuration
!
IF (TLES%LLES_CALL .AND. KSPLT==1) THEN
  CALL SECOND_MNH(ZTIME1)
  CALL LES_MEAN_SUBGRID( ZFLX, TLES%X_LES_SUBGRID_U2 ) 
  CALL LES_MEAN_SUBGRID( -ZWORK, TLES%X_LES_RES_ddxa_U_SBG_UaU , .TRUE.)
  CALL SECOND_MNH(ZTIME2)
  TLES%XTIME_LES = TLES%XTIME_LES + ZTIME2 - ZTIME1
END IF

!
!*      10.   < V'V'>
!             -------
!
!!! wait for the computation of ZWORK and PDP (that uses ZFLX)
!!$acc ! wait(2)
!
! Computes the V variance
IF (.NOT. O2D) THEN
   !$acc kernels present_cr(gz_w_m_pwm,zflx) ! async(3)
   !$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
      ZFLX(:,:,:)= (2./3.) * PTKEM(:,:,:)                                  &
           - XCMFS * PK(:,:,:) *( (4./3.) * GY_V_M_PVM(:,:,:)                        &
           -(2./3.) * ( GX_U_M_PUM(:,:,:)                      &
           +GZ_W_M_PWM(:,:,:)                ) )
   !$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
   !$acc end kernels
  !! &  to be tested
  !!  + XCMFB *  PLM / SQRT(PTKEM) * (-2./3.) * PTP 
  !
ELSE
   !$acc kernels present_cr(gz_w_m_pwm,zflx) ! async(3)
   !$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
      ZFLX(:,:,:)= (2./3.) * PTKEM(:,:,:)                           &
           - XCMFS * PK(:,:,:) *(-(2./3.) * ( GX_U_M_PUM(:,:,:)        &
                                      +GZ_W_M_PWM(:,:,:)     ) )  
   !$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
   !$acc end kernels
  !! &  to be tested
  !!  + XCMFB *  PLM / SQRT(PTKEM) * (-2./3.) * PTP 
  !
END IF
!
!$acc kernels ! async(3)
!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
ZFLX(:,:,IKE+1) = ZFLX(:,:,IKE) 
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)
!$acc end kernels
!
! ! !!! wait for the computation of ZFLX to begin the update
! ! !$acc wait(3)
! ! !$acc update self(ZFLX(:,:,IKB+1:)) async(10)
!
!$acc kernels present_cr(zdv_dy,zflx) ! async(3)
!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
   ZFLX(:,:,IKB)   = (2./3.) * PTKEM(:,:,IKB)                           &
        - XCMFS * PK(:,:,IKB) * 2. * ZDV_DY(:,:)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)
!$acc end kernels

!!           & to be tested
!! + XCMFB * PLM(:,:,IKB:IKB) /SQRT(PTKEM(:,:,IKB:IKB)) *         &
!!   (-2./3.) * PTP(:,:,IKB:IKB)
!
! extrapolates this flux under the ground with the surface flux
!$acc kernels present_cr(ZFLX,ZDIRSINZW)
!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
ZFLX(:,:,IKB-1) =                                                            &
        PTAU11M(:,:) * PSINSLOPE(:,:)**2 * PDIRCOSZW(:,:)**2                 &         
  +2. * PTAU12M(:,:) * PCOSSLOPE(:,:)* PSINSLOPE(:,:) * PDIRCOSZW(:,:)       &
  +     PTAU22M(:,:) * PCOSSLOPE(:,:)**2                                     &
  +     PTAU33M(:,:) * PSINSLOPE(:,:)**2 * ZDIRSINZW(:,:)**2                 &
  -2. * PCDUEFF(:,:)*       (                                                &
      PUSLOPEM(:,:) * PSINSLOPE(:,:)**2 * ZDIRSINZW(:,:) * PDIRCOSZW(:,:)    &
    + PVSLOPEM(:,:) * PCOSSLOPE(:,:)    * PSINSLOPE(:,:) * ZDIRSINZW(:,:)    )
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)
!$acc end kernels
!
!$acc kernels ! async(3)
ZFLX(:,:,IKB-1) = 2. * ZFLX(:,:,IKB-1) -  ZFLX(:,:,IKB)
!$acc end kernels
!
!
! ! !!! wait for the computation of ZFLX(:,:,1:IKB) to begin the update
! ! !$acc update self(ZFLX(:,:,IKB+1:)) async(3)
! ! !
! ! !!! and wait for the update self(ZFLX(...)) to complete
! ! !$acc wait(10)
! ! !$acc wait(3)
!
!!$acc ! wait(3)
#ifndef MNH_OPENACC
CALL UPDATE_HALO_ll(TZFIELDS_ll, IINFO_ll)
#else
CALL GET_HALO_D(ZFLX,HNAME='TURB_HOR_DYN_CORR::ZFLX')
#endif
!
IF ( TPFILE%LOPENED .AND. TURBN%LTURB_FLX ) THEN
  !$acc update self( ZFLX )
  ! stores <V V>  
  TZFIELD = TFIELDMETADATA(     &
    CMNHNAME   = 'V_VAR',       &
    CSTDNAME   = '',            &
    CLONGNAME  = 'V_VAR',       &
    CUNITS     = 'm2 s-2',      &
    CDIR       = 'XY',          &
    CCOMMENT   = 'X_Y_Z_V_VAR', &
    NGRID      = 1,             &
    NTYPE      = TYPEREAL,      &
    NDIMS      = 3,             &
    LTIMEDEP   = .TRUE.         )
  CALL IO_Field_write(TPFILE,TZFIELD,ZFLX)
END IF
!
!!! wait for the computation of PRUS (that uses temporary variables)
!!$acc ! wait(1)
!
!
!
! Complete the V tendency
IF (.NOT. O2D) THEN
  IF (.NOT. OFLAT) THEN
    PRVS(:,:,:)=PRVS                                          &
                -DYM(PRHODJ * ZFLX / MYF(PDYY) )              &
                +DZF( PDZY / MZM(PDYY) *                      &
                MYM( MZM(PRHODJ*ZFLX) * PINV_PDZZ ) )
  ELSE
    PRVS(:,:,:)=PRVS -DYM(PRHODJ * ZFLX / MYF(PDYY) )
  END IF
!
! Contribution to the dynamic production of TKE:
  IF (KSPLT==1) THEN
     !$acc kernels
     !$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
     ZWORK(:,:,:)     = - ZFLX(:,:,:) * GY_V_M_PVM(:,:,:)
     !$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
     !$acc end kernels
  END IF
ELSE
  !$acc kernels ! async(2)
  !$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
  ZWORK(:,:,:)     = 0.
  !$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
  !$acc end kernels
END IF
!
IF (KSPLT==1) THEN
  !
  ! evaluate the dynamic production at w(IKB+1) in PDP(IKB)
  !
   !$acc kernels present_cr(zdv_dy,zwork) ! async(2)
   !$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
      ZWORK(:,:,IKB) = 0.5* ( -ZFLX(:,:,IKB)*ZDV_DY(:,:) + ZWORK(:,:,IKB+1) )
   !$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)
   !$acc end kernels
  !
  !$acc kernels ! async(2)
  !$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
  PDP(:,:,:) = PDP(:,:,:) + ZWORK(:,:,:)
  !$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
  !$acc end kernels
END IF
!
! Storage in the LES configuration
!
IF (TLES%LLES_CALL .AND. KSPLT==1) THEN
  CALL SECOND_MNH(ZTIME1)
  CALL LES_MEAN_SUBGRID( ZFLX, TLES%X_LES_SUBGRID_V2 ) 
  CALL LES_MEAN_SUBGRID( -ZWORK, TLES%X_LES_RES_ddxa_V_SBG_UaV , .TRUE.)
  CALL SECOND_MNH(ZTIME2)
  TLES%XTIME_LES = TLES%XTIME_LES + ZTIME2 - ZTIME1
END IF
!
!*      11.   < W'W'>
!             -------
!
! Computes the W variance
IF (.NOT. O2D) THEN
   !$acc kernels present_cr(gy_v_m_pvm,zflx) ! async(2)
   !$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
      ZFLX(:,:,:) = (2./3.) * PTKEM(:,:,:)                                  &
           - XCMFS * PK(:,:,:) *( (4./3.) * GZ_W_M_PWM(:,:,:)                        &
           -(2./3.) * ( GX_U_M_PUM(:,:,:)                      &
           +GY_V_M_PVM(:,:,:)                ) )
   !$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
  !$acc end kernels
  !!  &  to be tested
  !!    -2.* XCMFB *  PLM / SQRT(PTKEM) * (-2./3.) * PTP 
ELSE
   !$acc kernels present_cr(gx_u_m_pum,zflx) ! async(2)
   !$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
      ZFLX(:,:,:)= (2./3.) * PTKEM(:,:,:)                           &
           - XCMFS * PK(:,:,:) *( (4./3.) * GZ_W_M_PWM(:,:,:)          &
           -(2./3.) * ( GX_U_M_PUM(:,:,:)           ) ) 
   !$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
   !$acc end kernels
  !!  &  to be tested
  !!    -2.* XCMFB *  PLM / SQRT(PTKEM) * (-2./3.) * PTP 
END IF
!
!$acc kernels ! async(2)
!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
ZFLX(:,:,IKE+1)= ZFLX(:,:,IKE)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)
!$acc end kernels
!!! wait for the computation of ZWORK, PDP and ZFLX
!!$acc ! wait(2)
!
!
!$acc kernels present_cr(zdw_dz,zflx) ! async(2)
!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
   ZFLX(:,:,IKB)   = (2./3.) * PTKEM(:,:,IKB)                           &
        - XCMFS * PK(:,:,IKB) * 2. * ZDW_DZ(:,:)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)

!$acc end kernels
!

!             &  to be tested
!   - 2.* XCMFB * PLM(:,:,IKB:IKB) /SQRT(PTKEM(:,:,IKB:IKB)) *             &
!  (-2./3.) * PTP(:,:,IKB:IKB)
! extrapolates this flux under the ground with the surface flux
!$acc kernels present_cr(ZFLX) ! async(3)
!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
ZFLX(:,:,IKB-1) =                                                     &
        PTAU11M(:,:) * ZDIRSINZW(:,:)**2                                &
  +     PTAU33M(:,:) * PDIRCOSZW(:,:)**2                                &
  +2. * PCDUEFF(:,:)* PUSLOPEM(:,:)  * ZDIRSINZW(:,:) * PDIRCOSZW(:,:) 
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)        
!$acc end kernels
! 
!
!!! wait for the computation of ZFLX(:,:,IKB-1) and ZFLX(:,:,IKB)
!!$acc ! wait(2) ! async(3)
!
!$acc kernels ! async(3)
!$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
ZFLX(:,:,IKB-1) = 2. * ZFLX(:,:,IKB-1) - ZFLX(:,:,IKB)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)
!$acc end kernels
!
IF ( TPFILE%LOPENED .AND. TURBN%LTURB_FLX ) THEN
  !!$acc ! wait(3)
  !$acc update self(ZFLX)
  ! stores <W W>  
  TZFIELD = TFIELDMETADATA(     &
    CMNHNAME   = 'W_VAR',       &
    CSTDNAME   = '',            &
    CLONGNAME  = 'W_VAR',       &
    CUNITS     = 'm2 s-2',      &
    CDIR       = 'XY',          &
    CCOMMENT   = 'X_Y_Z_W_VAR', &
    NGRID      = 1,             &
    NTYPE      = TYPEREAL,      &
    NDIMS      = 3,             &
    LTIMEDEP   = .TRUE.         )
  CALL IO_Field_write(TPFILE,TZFIELD,ZFLX)
END IF
!
!
!!! wait for the computation of PRVS (that uses temporary variables)
!!$acc ! wait(1)
!

!
! Complete the W tendency
!
!PRWS(:,:,:)=PRWS(:,:,:) - DZM( PRHODJ*ZFLX/MZF(PDZZ) )
!$acc kernels ! async(2)
!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
ZDFDDWDZ(:,:,:)    = - XCMFS * PK(:,:,:) * (4./3.)
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
!$acc end kernels
!$acc kernels ! async(2)
!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKB)
ZDFDDWDZ(:,:,1:IKB) = 0.
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKB)
!$acc end kernels
!
!!! wait for the computation of ZFLX(:,:,IKB-1) and ZDFDDWDZ
!!$acc ! wait(3) async(2)
!!$acc ! wait(2)
!
CALL TRIDIAG_W(PWM,ZFLX,ZDFDDWDZ,PTSTEP,ZMZF_DZZ,PRHODJ,ZWP)
!
PRWS = PRWS(:,:,:) + MZM(PRHODJ(:,:,:))*(ZWP(:,:,:)-PWM(:,:,:))/PTSTEP
!
!* recomputes flux using guess of W
!
GZ_W_M_ZWP = GZ_W_M(ZWP,PDZZ)
!$acc kernels present_cr(gz_w_m_pwm,zflx) ! async(2)
!$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=IKB+1:IKT)
   ZFLX(:,:,:)=ZFLX(:,:,:) &
        - XCMFS * PK(:,:,:) * (4./3.) * (GZ_W_M_ZWP(:,:,:) - GZ_W_M_PWM(:,:,:))
!$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=IKB+1:IKT)
!$acc end kernels
!
IF (KSPLT==1) THEN
   !Contribution to the dynamic production of TKE:
   !$acc kernels present_cr(gz_w_m_zwp,zwork) ! async(2)
   !$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
      ZWORK(:,:,:) = - ZFLX(:,:,:) * GZ_W_M_ZWP(:,:,:)
   !$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)   
   !$acc end kernels
  !
  ! evaluate the dynamic production at w(IKB+1) in PDP(IKB)
  !
   !$acc kernels present_cr(zdw_dz,zwork) ! async(2)
   !$mnh_expand_array(JI=1:IIT,JJ=1:IJT)
      ZWORK(:,:,IKB) = 0.5* ( -ZFLX(:,:,IKB)*ZDW_DZ(:,:) + ZWORK(:,:,IKB+1) )
   !$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT)
   !$acc end kernels
  !
  !$acc kernels ! async(2)
  !$mnh_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
  PDP(:,:,:) = PDP(:,:,:) + ZWORK(:,:,:)
  !$mnh_end_expand_array(JI=1:IIT,JJ=1:IJT,JK=1:IKT)
  !$acc end kernels
END IF
!
! Storage in the LES configuration
!
!
IF (TLES%LLES_CALL .AND. KSPLT==1) THEN
  CALL SECOND_MNH(ZTIME1)
  CALL LES_MEAN_SUBGRID( ZFLX, TLES%X_LES_SUBGRID_W2 ) 
  CALL LES_MEAN_SUBGRID( -ZWORK, TLES%X_LES_RES_ddxa_W_SBG_UaW , .TRUE.)
  CALL LES_MEAN_SUBGRID( GZ_M_M(PTHLM,PDZZ)*ZFLX, TLES%X_LES_RES_ddxa_Thl_SBG_UaW , .TRUE.)
  CALL LES_MEAN_SUBGRID(ZFLX*MZF(GZ_M_W(PTHLM,PDZZ)),TLES%X_LES_RES_ddz_Thl_SBG_W2)
  IF (KRR>=1) THEN
    CALL LES_MEAN_SUBGRID( GZ_M_M(PRM(:,:,:,1),PDZZ)*ZFLX, &
                           TLES%X_LES_RES_ddxa_Rt_SBG_UaW , .TRUE.)
    CALL LES_MEAN_SUBGRID(ZFLX*MZF(GZ_M_W(PRM(:,:,:,1),PDZZ)), &
                           TLES%X_LES_RES_ddz_Rt_SBG_W2)
  END IF
  DO JSV=1,KSV
    CALL LES_MEAN_SUBGRID( GZ_M_M(PSVM(:,:,:,JSV),PDZZ)*ZFLX, &
                           TLES%X_LES_RES_ddxa_Sv_SBG_UaW(:,:,:,JSV) , .TRUE.)
    CALL LES_MEAN_SUBGRID(ZFLX*MZF(GZ_M_W(PSVM(:,:,:,JSV),PDZZ)), &
                           TLES%X_LES_RES_ddz_Sv_SBG_W2(:,:,:,JSV))
  END DO
  CALL SECOND_MNH(ZTIME2)
  TLES%XTIME_LES = TLES%XTIME_LES + ZTIME2 - ZTIME1
END IF
!
CALL CLEANLIST_ll(TZFIELDS_ll)
!
!
END SUBROUTINE TURB_HOR_DYN_CORR
END MODULE MODE_TURB_HOR_DYN_CORR
