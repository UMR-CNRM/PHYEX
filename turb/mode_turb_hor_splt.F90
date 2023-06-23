!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
MODULE MODE_TURB_HOR_SPLT
IMPLICIT NONE
CONTAINS
           SUBROUTINE TURB_HOR_SPLT(D,CST,CSTURB,TURBN,NEBN,TLES,    &
                      KSPLIT, KRR,KRRL,KRRI,KSV, KSV_LGBEG,KSV_LGEND,&
                      PTSTEP,HLBCX,HLBCY, OFLAT, O2D, ONOMIXLG,      &
                      OOCEAN,OCOMPUTE_SRC,OBLOWSNOW,PRSNOW,          &
                      TPFILE, KHALO,                                 &
                      PDXX,PDYY,PDZZ,PDZX,PDZY,PZZ,                  &
                      PDIRCOSXW,PDIRCOSYW,PDIRCOSZW,                 &
                      PCOSSLOPE,PSINSLOPE,                           &
                      PRHODJ,PTHVREF,                                & 
                      PSFTHM,PSFRM,PSFSVM,                           &
                      PCDUEFF,PTAU11M,PTAU12M,PTAU22M,PTAU33M,       &
                      PUM,PVM,PWM,PUSLOPEM,PVSLOPEM,PTHLM,PRM,PSVM,  &
                      PTKEM,PLM,PLEPS,                               &
                      PLOCPEXNM,PATHETA,PAMOIST,PSRCM,PFRAC_ICE,     &
                      PDP,PTP,PSIGS,                                 &
                      PTRH,                                          &
                      PRUS,PRVS,PRWS,PRTHLS,PRRS,PRSVS               )
!     ################################################################
!
!
!!****  *TURB_HOR* -routine to compute the source terms in the meso-NH
!!               model equations due to the non-vertical turbulent fluxes.
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to compute the non-vertical
!     turbulent fluxes of the evolutive variables and give back the 
!     source terms to the main program.
!
!!**  METHOD
!!    ------
!!     Complementary 3D calculations when running at high resolution;
!!    The non-vertical turbulent fluxes are computed explicitly. The 
!!    contributions are cumulated in PRvarS and in DP and TP of TKE
!
! d(rho*T) = -d(rho*u'T'/dxx) -d(-rho*u'T'*dzx/dxx/dzz)
! / dt        / dx             /dz
!!    
!!
!!      Near the bottom of the model, uncentred evaluation of vertical 
!!    gradients are required because no field values are available under 
!!    the level where the gradient must be evaluated. In this case, the 
!!    gradient is computed with a second order accurate uncentred scheme 
!!    according to:
!!
!!        D FF           dzz3                       (dzz3+dzz4)   
!!        ----  = -  ----------------- FF(4)  +  ----------------- FF(3)   
!!        D z         (dzz3+dzz4) dzz4              dzz3 dzz4 
!!  
!!                    dzz4 + 2 dzz3          
!!                -  ----------------- FF(2)
!!                    (dzz3+dzz4) dzz3
!!
!!      where the values are taken from:
!!
!!                  -----    FF(5)
!!                    | 
!!                    |   dzz5
!!                    |    
!!                  -----    FF(4)
!!                    | 
!!                    |   dzz4
!!                    |    
!!                  -----    FF(3)
!!                    | 
!!                    |   dzz3
!!                    |    
!!                  -----    FF(2)    , (D FF / DZ)
!!                    |   dzz2 * 0.5
!!                  -----    ground
!!
!!
!!
!!    EXTERNAL
!!    --------
!!      GX_M_U, GY_M_V
!!      GX_M_M, GY_M_M, GZ_M_M
!!      GY_U_UV,GX_V_UV
!!      GX_U_M, GY_V_M, GZ_W_M
!!      GX_W_UW,GY_W_UW
!!                             :  Cartesian vertical gradient operators 
!!                               
!!
!!      MXM,MXF,MYM,MYF,MZM,MZF
!!                             :  Shuman functions (mean operators)     
!!      DXM,DXF.DYM,DYF,DZM,DZF
!!                             :  Shuman functions (difference operators)     
!!
!!       
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST : contains physical constants
!!
!!           XG         : gravity constant
!!
!!      Module MODD_CTURB: contains the set of constants for
!!                        the turbulence scheme
!!
!!           XCMFS,XCMFB : cts for the momentum flux
!!           XCSHF       : ct for the sensible heat flux
!!           XCHF        : ct for the moisture flux
!!           XCTV,XCHV   : cts for the T and moisture variances
!!
!!      Module MODD_PARAMETERS
!!
!!           JPVEXT     : number of vertical external points
!!
!!      Module MODD_CONF
!!
!!           HPROGRAM
!!           
!!
!!    REFERENCE
!!    ---------
!!      Book 2 of documentation (routine TURB_HOR)
!!      Book 1 of documentation (Chapter: Turbulence)
!!
!!    AUTHOR
!!    ------
!!      Joan Cuxart             * INM and Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       Aug 29, 1994
!!      Modifications: Feb 14, 1995 (J.Cuxart and J.Stein) 
!!                                  Doctorization and Optimization
!!                     March 21, 1995 (J.M. Carriere) 
!!                                  Introduction of cloud water
!!                     June  14, 1995 (J. Stein) 
!!                                  rm the ZVTPV computation + bug in the all 
!!                                  or nothing condens. case
!!                     June 28, 1995 (J.Cuxart)  Add the LES tools 
!!                     Sept 19, 1995 (J. Stein) change the surface flux
!!               computations
!!                     Nov  13, 1995 (J. Stein) include the tangential fluxes
!!               bug in <u'w'> at the surface
!!                     Nov  27, 1997 (V. Saravane) spliting of the routine
!!                     Nov  27, 1997 (V. Masson) clearing of the routine
!!                     Mar  07, 2001 (V. Masson and J. Stein) time splitting
!!                                   + major bugs correction for slopes
!!                     Nov  06, 2002 (V. Masson) LES budgets
!!                     Feb  20, 2003 (JP Pinty)  Add PFRAC_ICE
!!                     Oct.2009  (C.Lac) Introduction of different PTSTEP according to the
!!                              advection schemes
!!                     J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!  P. Wautelet 20/05/2019: add name argument to ADDnFIELD_ll + new ADD4DFIELD_ll subroutine
!! --------------------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
!
USE MODD_CST, ONLY: CST_t
USE MODD_CTURB, ONLY: CSTURB_t
USE MODD_DIMPHYEX, ONLY : DIMPHYEX_t
USE MODD_LES, ONLY: TLES_t
USE MODD_TURB_n, ONLY: TURB_t
USE MODD_NEB_n, ONLY: NEB_t
!
USE MODD_IO, ONLY: TFILEDATA
USE MODD_PARAMETERS
!
!
USE MODI_SHUMAN 
USE MODE_TURB_HOR
USE MODE_TURB_HOR_TKE
!
USE MODE_ll
!
IMPLICIT NONE
!
!
!*       0.1  declaration of arguments
!
!
TYPE(DIMPHYEX_t),       INTENT(IN)   :: D
TYPE(CST_t),            INTENT(IN)   :: CST
TYPE(CSTURB_t),         INTENT(IN)   :: CSTURB
TYPE(TURB_t),           INTENT(IN)   :: TURBN
TYPE(NEB_t),            INTENT(IN)   :: NEBN
TYPE(TLES_t),           INTENT(INOUT):: TLES          ! modd_les structure
INTEGER,                INTENT(IN)   :: KSPLIT        ! number of time splitting
INTEGER,                INTENT(IN)   :: KRR           ! number of moist var.
INTEGER,                INTENT(IN)   :: KRRL          ! number of liquid water var.
INTEGER,                INTENT(IN)   :: KRRI          ! number of ice water var.
INTEGER,                INTENT(IN)   :: KSV,KSV_LGBEG,KSV_LGEND ! number of sv var.
REAL,                   INTENT(IN)   ::  PTSTEP       ! timestep 
CHARACTER (LEN=*), DIMENSION(:), INTENT(IN)       ::  HLBCX,HLBCY
LOGICAL,                INTENT(IN)   ::  OFLAT        ! Logical for zero ororography
LOGICAL,                INTENT(IN)   ::  ONOMIXLG     ! to use turbulence for lagrangian variables (modd_conf)
LOGICAL,                INTENT(IN)   ::  O2D          ! Logical for 2D model version (modd_conf)
LOGICAL,                INTENT(IN)   ::  OOCEAN       ! switch for Ocean model version
LOGICAL,                INTENT(IN)   ::  OCOMPUTE_SRC ! flag to define dimensions of SIGS and SRCT variables
LOGICAL,                INTENT(IN)   ::  OBLOWSNOW    ! switch to activate pronostic blowing snow
INTEGER,                INTENT(IN)   ::  KHALO        ! Size of the halo for parallel distribution
REAL,                   INTENT(IN)   ::  PRSNOW       ! Ratio for diffusion coeff. scalar (blowing snow)
TYPE(TFILEDATA),          INTENT(IN)    ::  TPFILE       ! Output file
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PDXX, PDYY, PDZZ, PDZX, PDZY 
                                                         ! Metric coefficients
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PZZ          ! vertical grid
REAL, DIMENSION(D%NIT,D%NJT),     INTENT(IN)    ::  PDIRCOSXW, PDIRCOSYW, PDIRCOSZW
! Director Cosinus along x, y and z directions at surface w-point
REAL, DIMENSION(D%NIT,D%NJT),   INTENT(IN)   ::  PCOSSLOPE       ! cosinus of the angle 
                                      ! between i and the slope vector
REAL, DIMENSION(D%NIT,D%NJT),   INTENT(IN)   ::  PSINSLOPE       ! sinus of the angle 
                                      ! between i and the slope vector
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PRHODJ       ! density * grid volume
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PTHVREF      ! ref. state VPT       
!
REAL, DIMENSION(D%NIT,D%NJT),     INTENT(IN)    ::  PSFTHM,PSFRM
REAL, DIMENSION(D%NIT,D%NJT,KSV),   INTENT(IN)    ::  PSFSVM       ! surface fluxes
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
REAL, DIMENSION(D%NIT,D%NJT,D%NKT,KRR), INTENT(IN)    ::  PRM          ! mixing ratios at t-1,
                              !  where PRM(:,:,:,1) = conservative mixing ratio
REAL, DIMENSION(D%NIT,D%NJT,D%NKT,KSV), INTENT(IN)    ::  PSVM         ! scalar var. at t-1
REAL, DIMENSION(D%NIT,D%NJT),      INTENT(IN)   ::  PUSLOPEM     ! wind component along the 
                                     ! maximum slope direction
REAL, DIMENSION(D%NIT,D%NJT),      INTENT(IN)   ::  PVSLOPEM     ! wind component along the 
                                     ! direction normal to the maximum slope one
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PTKEM        ! TKE at time t- dt
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PLM          ! Turb. mixing length
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PLEPS        ! dissipative length
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PLOCPEXNM    ! Lv(T)/Cp/Exner at time t-1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PATHETA      ! coefficients between 
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PAMOIST      ! s and Thetal and Rnp

REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PSRCM
                                  ! normalized 2nd-order flux
                                  ! s'r'c/2Sigma_s2 at t-1 multiplied by Lambda_3
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    ::  PFRAC_ICE    ! ri fraction of rc+ri
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(INOUT) ::  PRUS, PRVS, PRWS, PRTHLS
REAL, DIMENSION(D%NIT,D%NJT,D%NKT,KRR), INTENT(INOUT) ::  PRRS    ! var. at t+1 -split-
REAL, DIMENSION(D%NIT,D%NJT,D%NKT,KSV), INTENT(INOUT) ::  PRSVS   ! var. at t+1 -split-
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(INOUT) ::  PDP,PTP      ! TKE production terms
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(OUT)   ::  PTRH
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(INOUT) ::  PSIGS
                                  ! IN: Vertical part of Sigma_s at t
                                  ! OUT: Total Sigma_s at t
!
!
!
!*       0.2  declaration of local variables
!
REAL,ALLOCATABLE,DIMENSION(:,:,:) :: ZK           ! Turbulent diffusion doef.
                                                  ! ZK = PLM * SQRT(PTKEM)
REAL,ALLOCATABLE,DIMENSION(:,:,:) :: ZINV_PDXX    ! 1./PDXX
REAL,ALLOCATABLE,DIMENSION(:,:,:) :: ZINV_PDYY    ! 1./PDYY
REAL,ALLOCATABLE,DIMENSION(:,:,:) :: ZINV_PDZZ    ! 1./PDZZ
REAL,ALLOCATABLE,DIMENSION(:,:,:) :: ZMZM_PRHODJ  ! MZM(PRHODJ)
!
INTEGER :: JSPLT ! current split
!
INTEGER :: IKB, IKE, IIB, IIE, IJB, IJE
INTEGER :: JRR, JSV
!
INTEGER :: ISV
INTEGER :: IINFO_ll
!
REAL,ALLOCATABLE,DIMENSION(:,:,:)   :: ZUM, ZVM, ZWM, ZTHLM, ZTKEM
REAL,ALLOCATABLE,DIMENSION(:,:,:,:) :: ZRM, ZSVM
REAL,ALLOCATABLE,DIMENSION(:,:,:)   :: ZRUS, ZRVS, ZRWS, ZRTHLS
REAL,ALLOCATABLE,DIMENSION(:,:,:,:) :: ZRRS, ZRSVS
!
!
TYPE(LIST_ll), POINTER, SAVE :: TZFIELDS_ll
!
! ---------------------------------------------------------------------------
!
!*       1.   PRELIMINARY COMPUTATIONS
!             ------------------------
!
IKB = 1.+JPVEXT
IKE = SIZE(PUM,3) - JPVEXT
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
ISV=SIZE(PSVM,4)
!
ALLOCATE(ZK(SIZE(PTHLM,1),SIZE(PTHLM,2),SIZE(PTHLM,3)))
ALLOCATE(ZINV_PDXX(SIZE(PDXX,1),SIZE(PDXX,2),SIZE(PDXX,3)))
ALLOCATE(ZINV_PDYY(SIZE(PDYY,1),SIZE(PDYY,2),SIZE(PDYY,3)))
ALLOCATE(ZINV_PDZZ(SIZE(PDZZ,1),SIZE(PDZZ,2),SIZE(PDZZ,3)))
ALLOCATE(ZMZM_PRHODJ(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3)))
!
ZINV_PDXX = 1./PDXX
ZINV_PDYY = 1./PDYY
ZINV_PDZZ = 1./PDZZ
ZMZM_PRHODJ = MZM(PRHODJ)
!
ZK(:,:,:)         = PLM(:,:,:) * SQRT(PTKEM(:,:,:))
!
NULLIFY(TZFIELDS_ll)
!
!--------------------------------------------------------------------
!
!*       2.   SPLIT PROCESS LOOP
!             ------------------
!
IF (KSPLIT>1) THEN
!
!*       2.1  allocations
!             -----------
!
  ALLOCATE(ZUM(SIZE(PUM,1),SIZE(PUM,2),SIZE(PUM,3)))
  ALLOCATE(ZVM(SIZE(PVM,1),SIZE(PVM,2),SIZE(PVM,3)))
  ALLOCATE(ZWM(SIZE(PWM,1),SIZE(PWM,2),SIZE(PWM,3)))
  ALLOCATE(ZSVM(SIZE(PSVM,1),SIZE(PSVM,2),SIZE(PSVM,3),SIZE(PSVM,4)))
  ALLOCATE(ZTHLM(SIZE(PTHLM,1),SIZE(PTHLM,2),SIZE(PTHLM,3)))
  ALLOCATE(ZTKEM(SIZE(PTKEM,1),SIZE(PTKEM,2),SIZE(PTKEM,3)))
  ALLOCATE(ZRM(SIZE(PRM,1),SIZE(PRM,2),SIZE(PRM,3),SIZE(PRM,4)))
  ALLOCATE(ZRUS(SIZE(PRUS,1),SIZE(PRUS,2),SIZE(PRUS,3)))
  ALLOCATE(ZRVS(SIZE(PRVS,1),SIZE(PRVS,2),SIZE(PRVS,3)))
  ALLOCATE(ZRWS(SIZE(PRWS,1),SIZE(PRWS,2),SIZE(PRWS,3)))
  ALLOCATE(ZRSVS(SIZE(PRSVS,1),SIZE(PRSVS,2),SIZE(PRSVS,3),SIZE(PRSVS,4)))
  ALLOCATE(ZRTHLS(SIZE(PRTHLS,1),SIZE(PRTHLS,2),SIZE(PRTHLS,3)))
  ALLOCATE(ZRRS(SIZE(PRRS,1),SIZE(PRRS,2),SIZE(PRRS,3),SIZE(PRRS,4)))
!
!
!*       2.2  list for parallel exchanges
!             ---------------------------
!
  CALL ADD3DFIELD_ll( TZFIELDS_ll, ZUM,               'TURB_HOR_SPLT::ZUM'               )
  CALL ADD3DFIELD_ll( TZFIELDS_ll, ZVM,               'TURB_HOR_SPLT::ZVM'               )
  CALL ADD3DFIELD_ll( TZFIELDS_ll, ZWM,               'TURB_HOR_SPLT::ZWM'               )
  CALL ADD3DFIELD_ll( TZFIELDS_ll, ZTHLM,             'TURB_HOR_SPLT::ZTHLM'             )
  CALL ADD3DFIELD_ll( TZFIELDS_ll, ZTKEM,             'TURB_HOR_SPLT::ZTKEM'             )
  CALL ADD4DFIELD_ll( TZFIELDS_ll, ZSVM(:,:,:,1:ISV), 'TURB_HOR_SPLT::ZSVM(:,:,:,1:ISV)' )
  CALL ADD4DFIELD_ll( TZFIELDS_ll, ZRM(:,:,:,1:KRR),  'TURB_HOR_SPLT::ZRM(:,:,:,1:KRR)'  )
!
!
!*       2.3  initializations
!             ---------------
!
!
  ZUM=PUM
  ZVM=PVM
  ZWM=PWM
  IF (ISV>0) ZSVM=PSVM
  ZTHLM=PTHLM
  ZTKEM=PTKEM
  IF (KRR>0) ZRM=PRM
  !
  ZRUS=PRUS*KSPLIT
  ZRVS=PRVS*KSPLIT
  ZRWS=PRWS*KSPLIT
  IF (ISV>0) ZRSVS=PRSVS*KSPLIT
  ZRTHLS=PRTHLS*KSPLIT
  IF (KRR>0) ZRRS=PRRS*KSPLIT

!
!*       2.4  split process
!             -------------
!
  DO JSPLT=1,KSPLIT
!
! compute the turbulent tendencies for the small time step
    CALL TURB_HOR(D,CST,CSTURB,TURBN,NEBN,TLES,                   &
                   JSPLT, KRR, KRRL, KRRI, PTSTEP,                &
                   KSV, KSV_LGBEG, KSV_LGEND, OFLAT,O2D, ONOMIXLG,&
                   OOCEAN,OCOMPUTE_SRC,OBLOWSNOW,                 &
                   TPFILE,                                        &
                   PDXX,PDYY,PDZZ,PDZX,PDZY,PZZ,PRSNOW,           &
                   PDIRCOSXW,PDIRCOSYW,PDIRCOSZW,                 &
                   PCOSSLOPE,PSINSLOPE,                           &
                   ZINV_PDXX, ZINV_PDYY, ZINV_PDZZ, ZMZM_PRHODJ,  &
                   ZK,                                            &
                   PRHODJ,PTHVREF,                                & 
                   PSFTHM,PSFRM,PSFSVM,                           &
                   PCDUEFF,PTAU11M,PTAU12M,PTAU22M,PTAU33M,       &
                   ZUM,ZVM,ZWM,PUSLOPEM,PVSLOPEM,ZTHLM,ZRM,ZSVM,  &
                   PTKEM,PLM,PLEPS,                               &
                   PLOCPEXNM,PATHETA,PAMOIST,PSRCM,PFRAC_ICE,     &
                   PDP,PTP,PSIGS,                                 &
                   ZRUS,ZRVS,ZRWS,ZRTHLS,ZRRS,ZRSVS               )
!
! horizontal transport of Tke
!
  CALL   TURB_HOR_TKE(JSPLT,TLES,OFLAT,O2D,                          &
                      PDXX,PDYY,PDZZ,PDZX,PDZY,                      &
                      ZINV_PDXX, ZINV_PDYY, ZINV_PDZZ, ZMZM_PRHODJ,  &
                      ZK, PRHODJ, ZTKEM,                             &
                      PTRH                                           )
!
!
! split temporal advance

    ZUM=PUM+(ZRUS/KSPLIT-PRUS)/MXM(PRHODJ)*PTSTEP
    ZVM=PVM+(ZRVS/KSPLIT-PRVS)/MYM(PRHODJ)*PTSTEP
    ZWM=PWM+(ZRWS/KSPLIT-PRWS)/ZMZM_PRHODJ*PTSTEP
    DO JSV=1,ISV
      ZSVM(:,:,:,JSV)=PSVM(:,:,:,JSV)+   &
        (ZRSVS(:,:,:,JSV)/KSPLIT-PRSVS(:,:,:,JSV))/PRHODJ*PTSTEP
    END DO
    ZTHLM=PTHLM+(ZRTHLS/KSPLIT-PRTHLS)/PRHODJ*PTSTEP
    ZTKEM=ZTKEM+PTRH*PTSTEP/KSPLIT
    DO JRR=1,KRR
      ZRM(:,:,:,JRR)=PRM(:,:,:,JRR)+   &
       (ZRRS(:,:,:,JRR)/KSPLIT-PRRS(:,:,:,JRR))/PRHODJ*PTSTEP
    END DO
!
! reinforce boundary conditions
!
    IF (JSPLT<KSPLIT-KHALO+1) CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
    !
    IF ( HLBCX(1) /= "CYCL" .AND. LWEST_ll()) THEN
       ZUM(IIB  ,:,:)=PUM(IIB  ,:,:)
       ZVM(IIB-1,:,:)=PVM(IIB-1,:,:)
       ZWM(IIB-1,:,:)=PWM(IIB-1,:,:)
       ZTHLM(IIB-1,:,:)=PTHLM(IIB-1,:,:)
       ZTKEM(IIB-1,:,:)=PTKEM(IIB-1,:,:)
       IF (ISV>0) ZSVM(IIB-1,:,:,:)=PSVM(IIB-1,:,:,:)
       IF (KRR>0) ZRM (IIB-1,:,:,:)=PRM (IIB-1,:,:,:)
     ENDIF
     !
     IF ( HLBCX(2) /= "CYCL" .AND. LEAST_ll()) THEN
       ZUM(IIE+1,:,:)=PUM(IIE+1,:,:)
       ZVM(IIE+1,:,:)=PVM(IIE+1,:,:)
       ZWM(IIE+1,:,:)=PWM(IIE+1,:,:)
       ZTHLM(IIE+1,:,:)=PTHLM(IIE+1,:,:)
       ZTKEM(IIE+1,:,:)=PTKEM(IIE+1,:,:)
       IF (ISV>0) ZSVM(IIE+1,:,:,:)=PSVM(IIE+1,:,:,:)
       IF (KRR>0) ZRM (IIE+1,:,:,:)=PRM(IIE+1,:,:,:)
     ENDIF
     !
     IF ( HLBCY(1) /= "CYCL" .AND. LSOUTH_ll()) THEN
       ZUM(:,IJB-1,:)=PUM(:,IJB-1,:)
       ZVM(:,IJB  ,:)=PVM(:,IJB  ,:)
       ZWM(:,IJB-1,:)=PWM(:,IJB-1,:)
       ZTHLM(:,IJB-1,:)=PTHLM(:,IJB-1,:)
       ZTKEM(:,IJB-1,:)=PTKEM(:,IJB-1,:)
       IF (ISV>0) ZSVM(:,IJB-1,:,:)=PSVM(:,IJB-1,:,:)
       IF (KRR>0) ZRM (:,IJB-1,:,:)=PRM (:,IJB-1,:,:)
     ENDIF
     !
     IF ( HLBCY(2) /= "CYCL" .AND. LNORTH_ll()) THEN
       ZUM(:,IJE+1,:)=PUM(:,IJE+1,:)
       ZVM(:,IJE+1,:)=PVM(:,IJE+1,:)
       ZWM(:,IJE+1,:)=PWM(:,IJE+1,:)
       ZTHLM(:,IJE+1,:)=PTHLM(:,IJE+1,:)
       ZTKEM(:,IJE+1,:)=PTKEM(:,IJE+1,:)
       IF (ISV>0) ZSVM(:,IJE+1,:,:)=PSVM(:,IJE+1,:,:)
       IF (KRR>0) ZRM (:,IJE+1,:,:)=PRM(:,IJE+1,:,:)
     ENDIF
     !
     ZUM(:,:,IKB-1)=ZUM(:,:,IKB)
     ZVM(:,:,IKB-1)=ZVM(:,:,IKB)
     ZWM(:,:,IKB-1)=ZWM(:,:,IKB)
     ZTHLM(:,:,IKB-1)=ZTHLM(:,:,IKB)
     ZTKEM(:,:,IKB-1)=ZTKEM(:,:,IKB)
     IF (ISV>0) ZSVM(:,:,IKB-1,:)=ZSVM(:,:,IKB,:)
     IF (KRR>0) ZRM (:,:,IKB-1,:)=ZRM (:,:,IKB,:)
     !
     ZUM(:,:,IKE+1)=ZUM(:,:,IKE)
     ZVM(:,:,IKE+1)=ZVM(:,:,IKE)
     ZWM(:,:,IKE+1)=ZWM(:,:,IKE)
     ZTHLM(:,:,IKE+1)=ZTHLM(:,:,IKE)
     ZTKEM(:,:,IKE+1)=ZTKEM(:,:,IKE)
     IF (ISV>0) ZSVM(:,:,IKE+1,:)=ZSVM(:,:,IKE,:)
     IF (KRR>0) ZRM (:,:,IKE+1,:)=ZRM (:,:,IKE,:)
     !
  END DO
!
!*       2.5  update the complete tendencies
!             ------------------------------
!
  PRUS=ZRUS/KSPLIT
  PRVS=ZRVS/KSPLIT
  PRWS=ZRWS/KSPLIT
  IF (ISV>0) PRSVS=ZRSVS/KSPLIT
  PRTHLS=ZRTHLS/KSPLIT
  IF (KRR>0) PRRS=ZRRS/KSPLIT
  PTRH=(ZTKEM-PTKEM)/PTSTEP
!
!*       2.6  deallocations
!             -------------
!
  DEALLOCATE(ZUM)
  DEALLOCATE(ZVM)
  DEALLOCATE(ZWM)
  DEALLOCATE(ZSVM)
  DEALLOCATE(ZTHLM)
  DEALLOCATE(ZTKEM)
  DEALLOCATE(ZRM)
  DEALLOCATE(ZRUS)
  DEALLOCATE(ZRVS)
  DEALLOCATE(ZRWS)
  DEALLOCATE(ZRSVS)
  DEALLOCATE(ZRTHLS)
  DEALLOCATE(ZRRS)
  !
  CALL CLEANLIST_ll(TZFIELDS_ll)
!
!-------------------------------------------------------------------
!
!*       4.   NO SPLIT PROCESS CASE
!             ---------------------
!
ELSE
!
  CALL TURB_HOR(D,CST,CSTURB,TURBN,NEBN,TLES,                  &
                1, KRR, KRRL, KRRI,  PTSTEP,                   &
                KSV, KSV_LGBEG, KSV_LGEND, OFLAT,O2D, ONOMIXLG,&                
                OOCEAN,OCOMPUTE_SRC,OBLOWSNOW,                 &
                TPFILE,                                        &
                PDXX,PDYY,PDZZ,PDZX,PDZY,PZZ,PRSNOW,           &
                PDIRCOSXW,PDIRCOSYW,PDIRCOSZW,                 &
                PCOSSLOPE,PSINSLOPE,                           &
                ZINV_PDXX, ZINV_PDYY, ZINV_PDZZ, ZMZM_PRHODJ,  &
                ZK,                                            &
                PRHODJ,PTHVREF,                                & 
                PSFTHM,PSFRM,PSFSVM,                           &
                PCDUEFF,PTAU11M,PTAU12M,PTAU22M,PTAU33M,       &
                PUM,PVM,PWM,PUSLOPEM,PVSLOPEM,PTHLM,PRM,PSVM,  &
                PTKEM,PLM,PLEPS,                               &
                PLOCPEXNM,PATHETA,PAMOIST,PSRCM,PFRAC_ICE,     &
                PDP,PTP,PSIGS,                                 &
                PRUS,PRVS,PRWS,PRTHLS,PRRS,PRSVS               )

! horizontal transport of Tke
!

  CALL   TURB_HOR_TKE(1,TLES,OFLAT,O2D,                              &
                      PDXX,PDYY,PDZZ,PDZX,PDZY,                      &
                      ZINV_PDXX, ZINV_PDYY, ZINV_PDZZ, ZMZM_PRHODJ,  &
                      ZK, PRHODJ, PTKEM,                             &
                      PTRH                                           )
!
END IF
!--------------------------------------------------------------------
!
DEALLOCATE(ZK)
DEALLOCATE(ZINV_PDXX)
DEALLOCATE(ZINV_PDYY)
DEALLOCATE(ZINV_PDZZ)
DEALLOCATE(ZMZM_PRHODJ)
!
END SUBROUTINE TURB_HOR_SPLT
END MODULE MODE_TURB_HOR_SPLT
