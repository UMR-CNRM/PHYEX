!MNH_LIC Copyright 1994-2025 CNRS, Meteo-France and Universite Paul Sabatier
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
USE MODD_CTURB,          ONLY: CTURB_XCMFS => XCMFS
use modd_field,          only: tfieldmetadata, TYPEREAL
USE MODD_IO,             ONLY: TFILEDATA
USE MODD_PARAMETERS, ONLY: JPVEXT
USE MODD_LES, ONLY: TLES_t
!
USE MODE_ll, ONLY: CLEANLIST_LL, UPDATE_HALO_LL
USE MODE_ARGSLIST_ll_PHY, ONLY: ADD3DFIELD_ll_PHY
USE MODE_IO_FIELD_WRITE_PHY,      ONLY: IO_FIELD_WRITE_PHY
!
#ifdef MNH_OPENACC

#endif
USE MODI_LES_MEAN_SUBGRID_PHY, ONLY: LES_MEAN_SUBGRID_PHY
USE MODE_TRIDIAG_W, ONLY: TRIDIAG_W
!
USE MODI_SECOND_MNH, ONLY: SECOND_MNH
USE MODE_SHUMAN_PHY, ONLY: DXF2D_PHY, DXM_PHY, DYF2D_PHY, DYM_PHY, DZF_PHY, MXF2D_PHY, MXF_PHY, MXM2D_PHY, MXM_PHY, MYF2D_PHY, &
  & MYF_PHY, MYM2D_PHY, MYM_PHY, MZF_PHY, MZM_PHY
USE MODE_GRADIENT_U_PHY, ONLY: GX_U_M_PHY
USE MODE_GRADIENT_V_PHY, ONLY: GY_V_M_PHY
USE MODE_GRADIENT_M_PHY, ONLY: GZ_M_M_PHY, GZ_M_W_PHY
USE MODE_GRADIENT_W_PHY, ONLY: GZ_W_M_PHY
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
TYPE(TFILEDATA),          INTENT(INOUT) ::  TPFILE       ! Output file
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    ::  PK          ! Turbulent diffusion doef.
                                                        ! PK = PLM * SQRT(PTKEM)
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    ::  PINV_PDZZ   ! 1./PDZZ
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    ::  PDXX, PDYY, PDZZ, PDZX, PDZY 
                                                         ! Metric coefficients
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    ::  PZZ          ! vertical grid
REAL, DIMENSION(D%NIJT),     INTENT(IN)    ::  PDIRCOSZW
! Director Cosinus along z directions at surface w-point
REAL, DIMENSION(D%NIJT),   INTENT(IN)   ::  PCOSSLOPE       ! cosinus of the angle 
                                      ! between i and the slope vector
REAL, DIMENSION(D%NIJT),   INTENT(IN)   ::  PSINSLOPE       ! sinus of the angle 
                                      ! between i and the slope vector
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    ::  PRHODJ       ! density * grid volume
!
REAL, DIMENSION(D%NIJT),   INTENT(IN)   ::  PCDUEFF      ! Cd * || u || at time t
REAL, DIMENSION(D%NIJT),   INTENT(IN)   ::  PTAU11M      ! <uu> in the axes linked 
       ! to the maximum slope direction and the surface normal and the binormal 
       ! at time t - dt
REAL, DIMENSION(D%NIJT),   INTENT(IN)   ::  PTAU12M      ! <uv> in the same axes
REAL, DIMENSION(D%NIJT),   INTENT(IN)   ::  PTAU22M      ! <vv> in the same axes
REAL, DIMENSION(D%NIJT),   INTENT(IN)   ::  PTAU33M      ! <ww> in the same axes
!
! Variables at t-1
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    ::  PUM,PVM,PWM,PTHLM
REAL, DIMENSION(D%NIJT,D%NKT,KRR), INTENT(IN)    ::  PRM
REAL, DIMENSION(D%NIJT,D%NKT,KSV), INTENT(IN)    ::  PSVM
REAL, DIMENSION(D%NIJT),      INTENT(IN)   ::  PUSLOPEM     ! wind component along the 
                                     ! maximum slope direction
REAL, DIMENSION(D%NIJT),      INTENT(IN)   ::  PVSLOPEM     ! wind component along the 
                                     ! direction normal to the maximum slope one
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    ::  PTKEM        ! TKE at time t- dt
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    ::  PLM          ! Turb. mixing length
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) ::  PRUS, PRVS, PRWS
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) ::  PDP,PTP      ! TKE production terms
!
!
!
!*       0.2  declaration of local variables
!
REAL, DIMENSION(D%NIJT,D%NKT) :: ZFLX,ZWORK,ZWKLES
    ! work arrays, PK is the turb. mixing coef.
!   
REAL, DIMENSION(D%NIJT) ::ZDIRSINZW 
      ! sinus of the angle between the vertical and the normal to the orography
INTEGER             :: IKB,IKE
                                    ! Index values for the Beginning and End
                                    ! mass points of the domain  
INTEGER             :: IKU,IKT,IIJT
INTEGER             :: JSV,JIJ,JK          ! scalar loop counter
!
REAL, DIMENSION(D%NIJT,D%NKT)  :: GX_U_M_PUM
REAL, DIMENSION(D%NIJT,D%NKT)  :: GY_V_M_PVM
REAL, DIMENSION(D%NIJT,D%NKT)  :: GZ_W_M_PWM
REAL, DIMENSION(D%NIJT,D%NKT)  :: GZ_W_M_ZWP
REAL, DIMENSION(D%NIJT,D%NKT)  :: ZMZF_DZZ   ! MZF(PDZZ)
REAL, DIMENSION(D%NIJT,D%NKT)  :: ZDFDDWDZ   ! formal derivative of the 
!                                                                   ! flux (variable: dW/dz)
REAL, DIMENSION(D%NIJT,D%NKT)  :: ZWP        ! W at future   time-step
!
REAL, DIMENSION(D%NIJT) :: ZDU_DZ_DZS_DX ! du/dz*dzs/dx surf
REAL, DIMENSION(D%NIJT) :: ZDV_DZ_DZS_DY ! dv/dz*dzs/dy surf
REAL, DIMENSION(D%NIJT) :: ZDU_DX        ! du/dx        surf
REAL, DIMENSION(D%NIJT) :: ZDV_DY        ! dv/dy        surf
REAL, DIMENSION(D%NIJT) :: ZDW_DZ        ! dw/dz        surf
!
INTEGER                :: IINFO_ll      ! return code of parallel routine
TYPE(LIST_ll), POINTER :: TZFIELDS_ll   ! list of fields to exchange

REAL :: ZTIME1, ZTIME2


REAL, DIMENSION(D%NIJT,1+JPVEXT:3+JPVEXT) :: ZCOEFF , ZDZZ
                                    ! coefficients for the uncentred gradient 
                                    ! computation near the ground
TYPE(TFIELDMETADATA) :: TZFIELD
!
REAL, DIMENSION(D%NIJT) :: ZWORKSHUGRAD_2D_1
REAL, DIMENSION(D%NIJT) :: ZWORKSHUGRAD_2D_2
REAL, DIMENSION(D%NIJT) :: ZWORKSHUGRAD_2D_3
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
REAL :: XCMFS
!
! --------------------------------------------------------------------------
!
!*       1.   PRELIMINARY COMPUTATIONS
!             ------------------------
!
XCMFS = CTURB_XCMFS
!
NULLIFY(TZFIELDS_ll)
!
IKB = 1+JPVEXT               
IKE = SIZE(PUM,2)-JPVEXT    
IKU = SIZE(PUM,2)
IIJT=D%NIJT
IKT=D%NKT
!
!

DO JIJ=1, IIJT
  ZDIRSINZW(JIJ) = SQRT( 1. - PDIRCOSZW(JIJ)**2 )
END DO

!
CALL GX_U_M_PHY(D, OFLAT, PUM, PDXX, PDZZ, PDZX, GX_U_M_PUM)
IF (.NOT. O2D) CALL GY_V_M_PHY(D, OFLAT, PVM, PDYY, PDZZ, PDZY, GY_V_M_PVM)
CALL GZ_W_M_PHY(D, PWM, PDZZ, GZ_W_M_PWM)
!
CALL MZF_PHY(D, PDZZ, ZMZF_DZZ)
!
CALL ADD3DFIELD_ll_PHY(D, TZFIELDS_ll, ZFLX, 'TURB_HOR_DYN_CORR::ZFLX' )


!  compute the coefficients for the uncentred gradient computation near the 
!  ground
!
!*       9.   < U'U'>
!             -------
!
! Computes the U variance
IF (.NOT. O2D) THEN
   
   DO JK=1, IKT
        DO JIJ=1, IIJT
          ZFLX(JIJ, JK)= (2./3.) * PTKEM(JIJ, JK)                            &
               - XCMFS * PK(JIJ, JK) *( (4./3.) * GX_U_M_PUM(JIJ, JK)        &
               -(2./3.) * ( GY_V_M_PVM(JIJ, JK)                     &
               +GZ_W_M_PWM(JIJ, JK)                ) )
     END DO
      END DO
   
   !!  &   to be tested later
  !!  + XCMFB *  PLM / SQRT(PTKEM) * (-2./3.) * PTP 
ELSE
  
  DO JK=1, IKT
    DO JIJ=1, IIJT
      ZFLX(JIJ, JK)= (2./3.) * PTKEM(JIJ, JK)                                  &
        - XCMFS * PK(JIJ, JK) *( (4./3.) * GX_U_M_PUM(JIJ, JK)                 &
                       -(2./3.) * ( GZ_W_M_PWM(JIJ, JK)             ) ) 
    END DO
  END DO
  
  !!  &   to be tested later
  !!  + XCMFB *  PLM / SQRT(PTKEM) * (-2./3.) * PTP 
END IF
!

DO JIJ=1, IIJT
     ZFLX(JIJ, IKE+1) = ZFLX(JIJ, IKE) 
END DO

!
!* prescription of du/dz and dv/dz with uncentered gradient at the surface
!  prescription of dw/dz at Dz/2 above ground using the continuity equation
!  using a Boussinesq hypothesis to remove the z dependance of rhod_ref
!  (div u = 0)
!
! ZDZZ(:,:) = MXM(PDZZ(:,IKB:IKB+2)) ! case not handled yet by pyft --expandAllArraysConcurrent and shumanFUNCtoCALL
CALL MXM2D_PHY(D, PDZZ(:,IKB), ZDZZ(:,IKB))
CALL MXM2D_PHY(D, PDZZ(:,IKB+1), ZDZZ(:,IKB+1))
CALL MXM2D_PHY(D, PDZZ(:,IKB+2), ZDZZ(:,IKB+2))

DO JIJ=1, IIJT
     ZCOEFF(JIJ, IKB+2)= - ZDZZ(JIJ, 2) /      &
          ( (ZDZZ(JIJ, 3)+ZDZZ(JIJ, 2)) * ZDZZ(JIJ, 3) )
     ZCOEFF(JIJ, IKB+1)=   (ZDZZ(JIJ, 3)+ZDZZ(JIJ, 2)) /      &
          ( ZDZZ(JIJ, 2) * ZDZZ(JIJ, 3) )
     ZCOEFF(JIJ, IKB)= - (ZDZZ(JIJ, 3)+2.*ZDZZ(JIJ, 2)) /      &
          ( (ZDZZ(JIJ, 3)+ZDZZ(JIJ, 2)) * ZDZZ(JIJ, 2) )
END DO

!

DO JIJ=1, D%NIJT
  ZWORKSHUGRAD_2D_1(JIJ)=(ZCOEFF(JIJ, IKB+2)*PUM(JIJ, IKB+2)       &
                            +ZCOEFF(JIJ, IKB+1)*PUM(JIJ, IKB+1)       &
                            +ZCOEFF(JIJ, IKB)*PUM(JIJ, IKB)       &
                            )* 0.5 * ( PDZX(JIJ, IKB+1)+PDZX(JIJ, IKB))
END DO

CALL MXF2D_PHY(D, ZWORKSHUGRAD_2D_1, ZWORKSHUGRAD_2D_2)
CALL MXF2D_PHY(D, PDXX(:,IKB), ZWORKSHUGRAD_2D_3)

DO JIJ=1, D%NIJT
  ZDU_DZ_DZS_DX(JIJ)=ZWORKSHUGRAD_2D_2(JIJ)/ ZWORKSHUGRAD_2D_3(JIJ)
END DO

!
CALL MYM2D_PHY(D, PDZZ(:,IKB), ZDZZ(:,IKB))
CALL MYM2D_PHY(D, PDZZ(:,IKB+1), ZDZZ(:,IKB+1))
CALL MYM2D_PHY(D, PDZZ(:,IKB+2), ZDZZ(:,IKB+2))
!

DO JIJ=1, IIJT
     ZCOEFF(JIJ, IKB+2)= - ZDZZ(JIJ, 2) /      &
          ( (ZDZZ(JIJ, 3)+ZDZZ(JIJ, 2)) * ZDZZ(JIJ, 3) )
     ZCOEFF(JIJ, IKB+1)=   (ZDZZ(JIJ, 3)+ZDZZ(JIJ, 2)) /      &
          ( ZDZZ(JIJ, 2) * ZDZZ(JIJ, 3) )
     ZCOEFF(JIJ, IKB)= - (ZDZZ(JIJ, 3)+2.*ZDZZ(JIJ, 2)) /      &
          ( (ZDZZ(JIJ, 3)+ZDZZ(JIJ, 2)) * ZDZZ(JIJ, 2) )
END DO



DO JIJ=1, D%NIJT
  ZWORKSHUGRAD_2D_1(JIJ)=(ZCOEFF(JIJ, IKB+2)*PVM(JIJ, IKB+2)       &
                            +ZCOEFF(JIJ, IKB+1)*PVM(JIJ, IKB+1)       &
                            +ZCOEFF(JIJ, IKB)*PVM(JIJ, IKB)       &
                            )* 0.5 * ( PDZY(JIJ, IKB+1)+PDZY(JIJ, IKB))
END DO

CALL MYF2D_PHY(D, ZWORKSHUGRAD_2D_1, ZWORKSHUGRAD_2D_2)
CALL MYF2D_PHY(D, PDYY(:,IKB), ZWORKSHUGRAD_2D_3)

DO JIJ=1, D%NIJT
  ZDV_DZ_DZS_DY(JIJ)=ZWORKSHUGRAD_2D_2(JIJ)/ ZWORKSHUGRAD_2D_3(JIJ)
END DO

!
!
CALL DXF2D_PHY(D, PUM(:,IKB), ZWORKSHUGRAD_2D_1)
CALL MXF2D_PHY(D, PDXX(:,IKB), ZWORKSHUGRAD_2D_2)

DO JIJ=1, D%NIJT
  ZDU_DX(JIJ)=  ZWORKSHUGRAD_2D_1(JIJ) / ZWORKSHUGRAD_2D_2(JIJ)  &
                - ZDU_DZ_DZS_DX(JIJ)
END DO

CALL DYF2D_PHY(D, PVM(:,IKB), ZWORKSHUGRAD_2D_1)
CALL MYF2D_PHY(D, PDYY(:,IKB), ZWORKSHUGRAD_2D_2)

DO JIJ=1, D%NIJT
  ZDV_DY(JIJ)=  ZWORKSHUGRAD_2D_1(JIJ) / ZWORKSHUGRAD_2D_2(JIJ) &
                - ZDV_DZ_DZS_DY(JIJ)
END DO

!

DO JIJ=1, IIJT
  ZDW_DZ(JIJ)=-ZDU_DX(JIJ)-ZDV_DY(JIJ)
END DO

!
!* computation 
!
!!! wait for the computation of ZFLX
!!$acc ! wait(2) async(4)
!!! wait for the computation of ZDW_DZ
!!$acc ! wait(4)
!
! ! !!! we can launch the update of ZFLX on the part that has already been computed
! ! !$acc update self(ZFLX(:,IKB+1:)) async(10)
!attention !!!!! je ne comprends pas pourquoi mais ce update plante à l'execution...
! du coup je ne peux pas faire de update self asynchrone...
!

DO JIJ=1, IIJT
     ZFLX(JIJ, IKB)   = (2./3.) * PTKEM(JIJ, IKB)                           &
          - XCMFS * PK(JIJ, IKB) * 2. * ZDU_DX(JIJ)
END DO


!!  &  to be tested later
!! + XCMFB * PLM(:,IKB:IKB) /SQRT(PTKEM(:,IKB:IKB)) *        &
!!   (-2./3.) * PTP(:,IKB:IKB)
!
! extrapolates this flux under the ground with the surface flux
!
!
!!! wait for the computation of ZDIRSINZW
!!$acc ! wait(1)
!

DO JIJ=1, IIJT
  ZFLX(JIJ, IKB-1) =                                                            &
          PTAU11M(JIJ) * PCOSSLOPE(JIJ)**2 * PDIRCOSZW(JIJ)**2                 &
    -2. * PTAU12M(JIJ) * PCOSSLOPE(JIJ)* PSINSLOPE(JIJ) * PDIRCOSZW(JIJ)       &
    +     PTAU22M(JIJ) * PSINSLOPE(JIJ)**2                                     &
    +     PTAU33M(JIJ) * PCOSSLOPE(JIJ)**2 * ZDIRSINZW(JIJ)**2                 &
    +2. * PCDUEFF(JIJ) *      (                                                &
        PVSLOPEM(JIJ) * PCOSSLOPE(JIJ)    * PSINSLOPE(JIJ) * ZDIRSINZW(JIJ)    &
      - PUSLOPEM(JIJ) * PCOSSLOPE(JIJ)**2 * ZDIRSINZW(JIJ) * PDIRCOSZW(JIJ)    )
END DO

! 
!!! wait for the computation of ZFLX(:,IKB) and ZFLX(:,IKB-1)
!!$acc ! wait(3) async(4)
!

DO JIJ=1, IIJT
     ZFLX(JIJ, IKB-1) = 2. * ZFLX(JIJ, IKB-1) -  ZFLX(JIJ, IKB)
END DO

!
!
!!! wait for the computation of ZFLX(:,IKB-1)
!!$acc ! wait(4)
!


! ! !!! we can launch the update of ZFLX on the rest
! ! !$acc update self(ZFLX(:,1:IKB)) async(11)
! ! !
! ! !!! and wait for the update self(ZFLX(...)) to complete
! ! !$acc wait(10)
! ! !$acc wait(11)
!attention !!!!! je ne comprends pas pourquoi mais le update self(ZFLX(:,IKB+1:)) plante à l'execution...
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
  CALL IO_FIELD_WRITE_PHY(D, TPFILE,TZFIELD,ZFLX)
END IF
!
! Complete the U tendency
IF (.NOT. OFLAT) THEN
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      ZWORKSHUGRAD_3D_1(JIJ, JK)=PRHODJ(JIJ, JK)*ZFLX(JIJ, JK)
    END DO
  END DO
  
  CALL MZM_PHY(D, ZWORKSHUGRAD_3D_1, ZWORKSHUGRAD_3D_2)
  CALL MXF_PHY(D, PDXX, ZWORKSHUGRAD_3D_3)
  CALL MZM_PHY(D, PDXX, ZWORKSHUGRAD_3D_4)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      ZWORKSHUGRAD_3D_5(JIJ, JK)=ZWORKSHUGRAD_3D_2(JIJ, JK) * PINV_PDZZ(JIJ, JK)
    END DO
  END DO
  
  CALL MXM_PHY(D, ZWORKSHUGRAD_3D_5, ZWORKSHUGRAD_3D_6)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      ZWORKSHUGRAD_3D_7(JIJ, JK)=PRHODJ(JIJ, JK) * ZFLX(JIJ, JK) / ZWORKSHUGRAD_3D_3(JIJ, JK)
    END DO
  END DO
  
  CALL DXM_PHY(D, ZWORKSHUGRAD_3D_7, ZWORKSHUGRAD_3D_8)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      ZWORKSHUGRAD_3D_9(JIJ, JK)=PDZX(JIJ, JK) / ZWORKSHUGRAD_3D_4(JIJ, JK) * ZWORKSHUGRAD_3D_6(JIJ, JK)
    END DO
  END DO
  
  CALL DZF_PHY(D, ZWORKSHUGRAD_3D_9, ZWORKSHUGRAD_3D_10)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      PRUS(JIJ, JK)=PRUS(JIJ, JK)                                            &
                  -ZWORKSHUGRAD_3D_8(JIJ, JK)                &
                  +ZWORKSHUGRAD_3D_10(JIJ, JK)
    END DO
  END DO
  
ELSE
  CALL MXF_PHY(D, PDXX, ZWORKSHUGRAD_3D_1)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      ZWORKSHUGRAD_3D_2(JIJ, JK)=PRHODJ(JIJ, JK) * ZFLX(JIJ, JK) / ZWORKSHUGRAD_3D_1(JIJ, JK)
    END DO
  END DO
  
  CALL DXM_PHY(D, ZWORKSHUGRAD_3D_2, ZWORKSHUGRAD_3D_3)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      PRUS(JIJ, JK)=PRUS(JIJ, JK) -ZWORKSHUGRAD_3D_3(JIJ, JK)
    END DO
  END DO
  
END IF
!
IF (KSPLT==1) THEN
  ! Contribution to the dynamic production of TKE:
   
   DO JK=1, IKT
        DO JIJ=1, IIJT
          ZWORK(JIJ, JK)     = - ZFLX(JIJ, JK) * GX_U_M_PUM(JIJ, JK)
     END DO
      END DO
  
  !
  ! evaluate the dynamic production at w(IKB+1) in PDP(IKB)
  !
  
   DO JIJ=1, IIJT
        ZWORK(JIJ, IKB) = 0.5* ( -ZFLX(JIJ, IKB)*ZDU_DX(JIJ) + ZWORK(JIJ, IKB+1) )
   END DO
   
  !
  
  PDP(:,:) = PDP(:,:) + ZWORK(:,:)
  
END IF
!
! Storage in the LES configuration
!
IF (TLES%LLES_CALL .AND. KSPLT==1) THEN
  CALL SECOND_MNH(ZTIME1)
  CALL LES_MEAN_SUBGRID_PHY(D, TLES, ZFLX, TLES%X_LES_SUBGRID_U2 ) 
  ZWKLES = -ZWORK
  CALL LES_MEAN_SUBGRID_PHY(D, TLES, ZWKLES, TLES%X_LES_RES_ddxa_U_SBG_UaU , .TRUE.)
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
   
   DO JK=1, IKT
        DO JIJ=1, IIJT
          ZFLX(JIJ, JK)= (2./3.) * PTKEM(JIJ, JK)                                  &
               - XCMFS * PK(JIJ, JK) *( (4./3.) * GY_V_M_PVM(JIJ, JK)                        &
               -(2./3.) * ( GX_U_M_PUM(JIJ, JK)                      &
               +GZ_W_M_PWM(JIJ, JK)                ) )
     END DO
      END DO
   
  !! &  to be tested
  !!  + XCMFB *  PLM / SQRT(PTKEM) * (-2./3.) * PTP 
  !
ELSE
   
   DO JK=1, IKT
        DO JIJ=1, IIJT
          ZFLX(JIJ, JK)= (2./3.) * PTKEM(JIJ, JK)                           &
               - XCMFS * PK(JIJ, JK) *(-(2./3.) * ( GX_U_M_PUM(JIJ, JK)        &
                                          +GZ_W_M_PWM(JIJ, JK)     ) )  
     END DO
      END DO
   
  !! &  to be tested
  !!  + XCMFB *  PLM / SQRT(PTKEM) * (-2./3.) * PTP 
  !
END IF
!

DO JIJ=1, IIJT
  ZFLX(JIJ, IKE+1) = ZFLX(JIJ, IKE) 
END DO

!
! ! !!! wait for the computation of ZFLX to begin the update
! ! !$acc wait(3)
! ! !$acc update self(ZFLX(:,IKB+1:)) async(10)
!

DO JIJ=1, IIJT
     ZFLX(JIJ, IKB)   = (2./3.) * PTKEM(JIJ, IKB)                           &
          - XCMFS * PK(JIJ, IKB) * 2. * ZDV_DY(JIJ)
END DO


!!           & to be tested
!! + XCMFB * PLM(:,IKB:IKB) /SQRT(PTKEM(:,IKB:IKB)) *         &
!!   (-2./3.) * PTP(:,IKB:IKB)
!
! extrapolates this flux under the ground with the surface flux

DO JIJ=1, IIJT
  ZFLX(JIJ, IKB-1) =                                                            &
          PTAU11M(JIJ) * PSINSLOPE(JIJ)**2 * PDIRCOSZW(JIJ)**2                 &         
    +2. * PTAU12M(JIJ) * PCOSSLOPE(JIJ)* PSINSLOPE(JIJ) * PDIRCOSZW(JIJ)       &
    +     PTAU22M(JIJ) * PCOSSLOPE(JIJ)**2                                     &
    +     PTAU33M(JIJ) * PSINSLOPE(JIJ)**2 * ZDIRSINZW(JIJ)**2                 &
    -2. * PCDUEFF(JIJ)*       (                                                &
        PUSLOPEM(JIJ) * PSINSLOPE(JIJ)**2 * ZDIRSINZW(JIJ) * PDIRCOSZW(JIJ)    &
      + PVSLOPEM(JIJ) * PCOSSLOPE(JIJ)    * PSINSLOPE(JIJ) * ZDIRSINZW(JIJ)    )
END DO

!

ZFLX(:,IKB-1) = 2. * ZFLX(:,IKB-1) -  ZFLX(:,IKB)

!
!
! ! !!! wait for the computation of ZFLX(:,1:IKB) to begin the update
! ! !$acc update self(ZFLX(:,IKB+1:)) async(3)
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
  CALL IO_FIELD_WRITE_PHY(D, TPFILE,TZFIELD,ZFLX)
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
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWORKSHUGRAD_3D_1(JIJ, JK)=PRHODJ(JIJ, JK)*ZFLX(JIJ, JK)
      END DO
    END DO
    
    CALL MZM_PHY(D, ZWORKSHUGRAD_3D_1, ZWORKSHUGRAD_3D_2)
    CALL MYF_PHY(D, PDYY, ZWORKSHUGRAD_3D_3)
    CALL MZM_PHY(D, PDYY, ZWORKSHUGRAD_3D_4)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWORKSHUGRAD_3D_5(JIJ, JK)=ZWORKSHUGRAD_3D_2(JIJ, JK) * PINV_PDZZ(JIJ, JK)
      END DO
    END DO
    
    CALL MYM_PHY(D, ZWORKSHUGRAD_3D_5, ZWORKSHUGRAD_3D_6)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWORKSHUGRAD_3D_7(JIJ, JK)=PRHODJ(JIJ, JK) * ZFLX(JIJ, JK) / ZWORKSHUGRAD_3D_3(JIJ, JK)
      END DO
    END DO
    
    CALL DYM_PHY(D, ZWORKSHUGRAD_3D_7, ZWORKSHUGRAD_3D_8)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWORKSHUGRAD_3D_9(JIJ, JK)=PDZY(JIJ, JK) / ZWORKSHUGRAD_3D_4(JIJ, JK) *                      &
                    ZWORKSHUGRAD_3D_6(JIJ, JK)
      END DO
    END DO
    
    CALL DZF_PHY(D, ZWORKSHUGRAD_3D_9, ZWORKSHUGRAD_3D_10)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        PRVS(JIJ, JK)=PRVS(JIJ, JK)                                          &
                    -ZWORKSHUGRAD_3D_8(JIJ, JK)              &
                    +ZWORKSHUGRAD_3D_10(JIJ, JK)
      END DO
    END DO
    
  ELSE
    CALL MYF_PHY(D, PDYY, ZWORKSHUGRAD_3D_1)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWORKSHUGRAD_3D_2(JIJ, JK)=PRHODJ(JIJ, JK) * ZFLX(JIJ, JK) / ZWORKSHUGRAD_3D_1(JIJ, JK)
      END DO
    END DO
    
    CALL DYM_PHY(D, ZWORKSHUGRAD_3D_2, ZWORKSHUGRAD_3D_3)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        PRVS(JIJ, JK)=PRVS(JIJ, JK) -ZWORKSHUGRAD_3D_3(JIJ, JK)
      END DO
    END DO
    
  END IF
!
! Contribution to the dynamic production of TKE:
  IF (KSPLT==1) THEN
     
     DO JK=1, IKT
       DO JIJ=1, IIJT
         ZWORK(JIJ, JK)     = - ZFLX(JIJ, JK) * GY_V_M_PVM(JIJ, JK)
       END DO
     END DO
     
  END IF
ELSE
  
  DO JK=1, IKT
    DO JIJ=1, IIJT
      ZWORK(JIJ, JK)     = 0.
    END DO
  END DO
  
END IF
!
IF (KSPLT==1) THEN
  !
  ! evaluate the dynamic production at w(IKB+1) in PDP(IKB)
  !
   
   DO JIJ=1, IIJT
        ZWORK(JIJ, IKB) = 0.5* ( -ZFLX(JIJ, IKB)*ZDV_DY(JIJ) + ZWORK(JIJ, IKB+1) )
   END DO
   
  !
  
  DO JK=1, IKT
    DO JIJ=1, IIJT
      PDP(JIJ, JK) = PDP(JIJ, JK) + ZWORK(JIJ, JK)
    END DO
  END DO
  
END IF
!
! Storage in the LES configuration
!
IF (TLES%LLES_CALL .AND. KSPLT==1) THEN
  CALL SECOND_MNH(ZTIME1)
  CALL LES_MEAN_SUBGRID_PHY(D, TLES, ZFLX, TLES%X_LES_SUBGRID_V2 ) 
  ZWKLES = -ZWORK
  CALL LES_MEAN_SUBGRID_PHY(D, TLES, ZWKLES, TLES%X_LES_RES_ddxa_V_SBG_UaV , .TRUE.)
  CALL SECOND_MNH(ZTIME2)
  TLES%XTIME_LES = TLES%XTIME_LES + ZTIME2 - ZTIME1
END IF
!
!*      11.   < W'W'>
!             -------
!
! Computes the W variance
IF (.NOT. O2D) THEN
   
   DO JK=1, IKT
        DO JIJ=1, IIJT
          ZFLX(JIJ, JK) = (2./3.) * PTKEM(JIJ, JK)                                  &
               - XCMFS * PK(JIJ, JK) *( (4./3.) * GZ_W_M_PWM(JIJ, JK)                        &
               -(2./3.) * ( GX_U_M_PUM(JIJ, JK)                      &
               +GY_V_M_PVM(JIJ, JK)                ) )
     END DO
      END DO
  
  !!  &  to be tested
  !!    -2.* XCMFB *  PLM / SQRT(PTKEM) * (-2./3.) * PTP 
ELSE
   
   DO JK=1, IKT
        DO JIJ=1, IIJT
          ZFLX(JIJ, JK)= (2./3.) * PTKEM(JIJ, JK)                           &
               - XCMFS * PK(JIJ, JK) *( (4./3.) * GZ_W_M_PWM(JIJ, JK)          &
               -(2./3.) * ( GX_U_M_PUM(JIJ, JK)           ) ) 
     END DO
      END DO
   
  !!  &  to be tested
  !!    -2.* XCMFB *  PLM / SQRT(PTKEM) * (-2./3.) * PTP 
END IF
!

DO JIJ=1, IIJT
  ZFLX(JIJ, IKE+1)= ZFLX(JIJ, IKE)
END DO

!!! wait for the computation of ZWORK, PDP and ZFLX
!!$acc ! wait(2)
!
!

DO JIJ=1, IIJT
     ZFLX(JIJ, IKB)   = (2./3.) * PTKEM(JIJ, IKB)                           &
          - XCMFS * PK(JIJ, IKB) * 2. * ZDW_DZ(JIJ)
END DO


!

!             &  to be tested
!   - 2.* XCMFB * PLM(:,IKB:IKB) /SQRT(PTKEM(:,IKB:IKB)) *             &
!  (-2./3.) * PTP(:,IKB:IKB)
! extrapolates this flux under the ground with the surface flux

DO JIJ=1, IIJT
  ZFLX(JIJ, IKB-1) =                                                     &
          PTAU11M(JIJ) * ZDIRSINZW(JIJ)**2                                &
    +     PTAU33M(JIJ) * PDIRCOSZW(JIJ)**2                                &
    +2. * PCDUEFF(JIJ)* PUSLOPEM(JIJ)  * ZDIRSINZW(JIJ) * PDIRCOSZW(JIJ) 
END DO

! 
!
!!! wait for the computation of ZFLX(:,IKB-1) and ZFLX(:,IKB)
!!$acc ! wait(2) ! async(3)
!

DO JIJ=1, IIJT
  ZFLX(JIJ, IKB-1) = 2. * ZFLX(JIJ, IKB-1) - ZFLX(JIJ, IKB)
END DO

!
IF ( TPFILE%LOPENED .AND. TURBN%LTURB_FLX ) THEN
  !!$acc ! wait(3)
  
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
  CALL IO_FIELD_WRITE_PHY(D, TPFILE,TZFIELD,ZFLX)
END IF
!
!
!!! wait for the computation of PRVS (that uses temporary variables)
!!$acc ! wait(1)
!

!
! Complete the W tendency
!
!PRWS(:,:)=PRWS(:,:) - DZM( PRHODJ*ZFLX/MZF(PDZZ) )

DO JK=1, IKT
  DO JIJ=1, IIJT
    ZDFDDWDZ(JIJ, JK)    = - XCMFS * PK(JIJ, JK) * (4./3.)
  END DO
END DO


DO JK=1, IKB
  DO JIJ=1, IIJT
    ZDFDDWDZ(JIJ, JK) = 0.
  END DO
END DO

!
!!! wait for the computation of ZFLX(:,IKB-1) and ZDFDDWDZ
!!$acc ! wait(3) async(2)
!!$acc ! wait(2)
!
CALL TRIDIAG_W(D,PWM,ZFLX,ZDFDDWDZ,PTSTEP,ZMZF_DZZ,PRHODJ,ZWP)
!
CALL MZM_PHY(D, PRHODJ, ZWORKSHUGRAD_3D_1)

DO JK=1, D%NKT
  DO JIJ=1, D%NIJT
    PRWS(JIJ, JK) = PRWS(JIJ, JK) + ZWORKSHUGRAD_3D_1(JIJ, JK)*(ZWP(JIJ, JK)-PWM(JIJ, JK))/PTSTEP
  END DO
END DO

!
!* recomputes flux using guess of W
!
CALL GZ_W_M_PHY(D, ZWP, PDZZ, GZ_W_M_ZWP)

DO JK=IKB+1, IKT
     DO JIJ=1, IIJT
       ZFLX(JIJ, JK)=ZFLX(JIJ, JK) &
            - XCMFS * PK(JIJ, JK) * (4./3.) * (GZ_W_M_ZWP(JIJ, JK) - GZ_W_M_PWM(JIJ, JK))
  END DO
   END DO

!
IF (KSPLT==1) THEN
   !Contribution to the dynamic production of TKE:
   
   DO JK=1, IKT
        DO JIJ=1, IIJT
          ZWORK(JIJ, JK) = - ZFLX(JIJ, JK) * GZ_W_M_ZWP(JIJ, JK)
     END DO
      END DO
   
  !
  ! evaluate the dynamic production at w(IKB+1) in PDP(IKB)
  !
   
   DO JIJ=1, IIJT
        ZWORK(JIJ, IKB) = 0.5* ( -ZFLX(JIJ, IKB)*ZDW_DZ(JIJ) + ZWORK(JIJ, IKB+1) )
   END DO
   
  !
  
  DO JK=1, IKT
    DO JIJ=1, IIJT
      PDP(JIJ, JK) = PDP(JIJ, JK) + ZWORK(JIJ, JK)
    END DO
  END DO
  
END IF
!
! Storage in the LES configuration
!
!
IF (TLES%LLES_CALL .AND. KSPLT==1) THEN
  CALL SECOND_MNH(ZTIME1)
  CALL LES_MEAN_SUBGRID_PHY(D, TLES, ZFLX, TLES%X_LES_SUBGRID_W2 ) 
  ZWKLES = -ZWORK
  CALL LES_MEAN_SUBGRID_PHY(D, TLES, ZWKLES, TLES%X_LES_RES_ddxa_W_SBG_UaW , .TRUE.)
  CALL GZ_M_M_PHY(D, PTHLM, PDZZ, ZWORKSHUGRAD_3D_1)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      ZWKLES(JIJ, JK) = ZWORKSHUGRAD_3D_1(JIJ, JK)*ZFLX(JIJ, JK)
    END DO
  END DO
  
  CALL LES_MEAN_SUBGRID_PHY(D, TLES, ZWKLES, TLES%X_LES_RES_ddxa_Thl_SBG_UaW , .TRUE.)
  CALL GZ_M_W_PHY(D, PTHLM, PDZZ, ZWORKSHUGRAD_3D_1)
  CALL MZF_PHY(D, ZWORKSHUGRAD_3D_1, ZWORKSHUGRAD_3D_2)
  
  DO JK=1, D%NKT
    DO JIJ=1, D%NIJT
      ZWKLES(JIJ, JK) = ZFLX(JIJ, JK)*ZWORKSHUGRAD_3D_2(JIJ, JK)
    END DO
  END DO
  
  CALL LES_MEAN_SUBGRID_PHY(D, TLES, ZWKLES, TLES%X_LES_RES_ddz_Thl_SBG_W2)
  IF (KRR>=1) THEN
    CALL GZ_M_M_PHY(D, PRM(:,:,1), PDZZ, ZWORKSHUGRAD_3D_1)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWKLES(JIJ, JK) = ZWORKSHUGRAD_3D_1(JIJ, JK)*ZFLX(JIJ, JK)
      END DO
    END DO
    
    CALL LES_MEAN_SUBGRID_PHY(D, TLES, ZWKLES, &
                           TLES%X_LES_RES_ddxa_Rt_SBG_UaW , .TRUE.)
    CALL GZ_M_W_PHY(D, PRM(:,:,1), PDZZ, ZWORKSHUGRAD_3D_1)
    CALL MZF_PHY(D, ZWORKSHUGRAD_3D_1, ZWORKSHUGRAD_3D_2)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWKLES(JIJ, JK) = ZFLX(JIJ, JK)*ZWORKSHUGRAD_3D_2(JIJ, JK)
      END DO
    END DO
    
    CALL LES_MEAN_SUBGRID_PHY(D, TLES, ZWKLES, &
                           TLES%X_LES_RES_ddz_Rt_SBG_W2)
  END IF
  DO JSV=1,KSV
    CALL GZ_M_M_PHY(D, PSVM(:,:,JSV), PDZZ, ZWORKSHUGRAD_3D_1)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWKLES(JIJ, JK) = ZWORKSHUGRAD_3D_1(JIJ, JK)*ZFLX(JIJ, JK)
      END DO
    END DO
    
    CALL LES_MEAN_SUBGRID_PHY(D, TLES, ZWKLES, &
                           TLES%X_LES_RES_ddxa_Sv_SBG_UaW(:,:,:,JSV) , .TRUE.)
    CALL GZ_M_W_PHY(D, PSVM(:,:,JSV), PDZZ, ZWORKSHUGRAD_3D_1)
    CALL MZF_PHY(D, ZWORKSHUGRAD_3D_1, ZWORKSHUGRAD_3D_2)
    
    DO JK=1, D%NKT
      DO JIJ=1, D%NIJT
        ZWKLES(JIJ, JK) = ZFLX(JIJ, JK)*ZWORKSHUGRAD_3D_2(JIJ, JK)
      END DO
    END DO
    
    CALL LES_MEAN_SUBGRID_PHY(D, TLES, ZWKLES, &
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
