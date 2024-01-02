!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 conv 2006/05/18 13:07:25
!-----------------------------------------------------------------
!     #################
      MODULE MODI_CONVECT_CLOSURE
!     #################
!
INTERFACE
!
       SUBROUTINE CONVECT_CLOSURE( KLON, KLEV,                                 &
                                 PPRES, PDPRES, PZ, PDXDY, PLMASS,           &
                                 PTHL, PTH, PRW, PRC, PRI, OTRIG1,           &
                                 PTHC, PRWC, PRCC, PRIC, PWSUB,              &
                                 KLCL, KDPL, KPBL, KLFS, KCTL, KML,          &
                                 PUMF, PUER, PUDR, PUTHL, PURW,              &
                                 PURC, PURI, PUPR,                           &
                                 PDMF, PDER, PDDR, PDTHL, PDRW,              &
                                 PTPR, PSPR, PDTEVR,                         &
                                 PCAPE, PTIMEC,                              &
                                 KFTSTEPS,                                   &
                                 PDTEVRF, PPRLFLX, PPRSFLX                   )
!
INTEGER,                   INTENT(IN) :: KLON   ! horizontal dimension
INTEGER,                   INTENT(IN) :: KLEV   ! vertical dimension
INTEGER, DIMENSION(KLON),  INTENT(IN) :: KLFS   ! index for level of free sink
INTEGER, DIMENSION(KLON),  INTENT(IN) :: KLCL   ! index lifting condens. level
INTEGER, DIMENSION(KLON),  INTENT(IN) :: KCTL   ! index for cloud top level
INTEGER, DIMENSION(KLON),  INTENT(IN) :: KDPL   ! index for departure level 
INTEGER, DIMENSION(KLON),  INTENT(IN) :: KPBL   ! index for top of source layer
INTEGER, DIMENSION(KLON),  INTENT(IN) :: KML    ! index for melting level
REAL, DIMENSION(KLON),  INTENT(INOUT) :: PTIMEC ! convection time step 
REAL, DIMENSION(KLON),     INTENT(IN) :: PDXDY  ! grid area (m^2)
REAL, DIMENSION(KLON,KLEV),INTENT(IN) :: PTHL   ! grid scale enthalpy (J/kg)
REAL, DIMENSION(KLON,KLEV),INTENT(IN) :: PTH    ! grid scale theta        
REAL, DIMENSION(KLON,KLEV),INTENT(IN) :: PRW    ! grid scale total water  
                                                ! mixing ratio 
REAL, DIMENSION(KLON,KLEV),INTENT(IN) :: PRC    ! grid scale r_c 
REAL, DIMENSION(KLON,KLEV),INTENT(IN) :: PRI    ! grid scale r_i 
LOGICAL, DIMENSION(KLON),  INTENT(IN) :: OTRIG1 ! logical to keep trace of 
                                                ! convective arrays modified in UPDRAFT
!   
!
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PPRES  ! pressure (P)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PDPRES ! pressure difference between 
                                                 ! bottom and top of layer (Pa)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PLMASS ! mass of model layer (kg)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PZ     ! height of model layer (m) 
REAL, DIMENSION(KLON),     INTENT(IN)  :: PCAPE  ! available potent. energy
INTEGER,                INTENT(OUT)   :: KFTSTEPS! maximum of fract time steps
                                                 ! only used for chemical tracers
!
!
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PUMF  ! updraft mass flux (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PUER  ! updraft entrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PUDR  ! updraft detrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PUPR  ! updraft precipitation in
                                                  ! flux units (kg water / s)
REAL, DIMENSION(KLON,KLEV), INTENT(IN)  :: PUTHL  ! updraft enthalpy (J/kg)
REAL, DIMENSION(KLON,KLEV), INTENT(IN)  :: PURW   ! updraft total water (kg/kg)
REAL, DIMENSION(KLON,KLEV), INTENT(IN)  :: PURC   ! updraft cloud water (kg/kg)
REAL, DIMENSION(KLON,KLEV), INTENT(IN)  :: PURI   ! updraft cloud ice   (kg/kg)
!
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PDMF  ! downdraft mass flux (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PDER  ! downdraft entrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PDDR  ! downdraft detrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(IN)   :: PDTHL ! downdraft enthalpy (J/kg)
REAL, DIMENSION(KLON,KLEV), INTENT(IN)   :: PDRW  ! downdraft total water (kg/kg)
REAL, DIMENSION(KLON),      INTENT(INOUT):: PTPR  ! total surf precipitation (kg/s)
REAL, DIMENSION(KLON),      INTENT(OUT)  :: PSPR  ! solid surf precipitation (kg/s)
REAL, DIMENSION(KLON),      INTENT(INOUT):: PDTEVR! donwndraft evapor. (kg/s)
!
REAL, DIMENSION(KLON,KLEV), INTENT(OUT)  :: PTHC  ! conv. adj. grid scale theta
REAL, DIMENSION(KLON,KLEV), INTENT(OUT)  :: PRWC  ! conv. adj. grid scale r_w 
REAL, DIMENSION(KLON,KLEV), INTENT(OUT)  :: PRCC  ! conv. adj. grid scale r_c 
REAL, DIMENSION(KLON,KLEV), INTENT(OUT)  :: PRIC  ! conv. adj. grid scale r_i 
REAL, DIMENSION(KLON,KLEV), INTENT(OUT)  :: PWSUB ! envir. compensating subsidence(Pa/s)
!
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PDTEVRF! downdraft evaporation rate
REAL, DIMENSION(KLON,KLEV), INTENT(OUT)  :: PPRLFLX! liquid precip flux
REAL, DIMENSION(KLON,KLEV), INTENT(OUT)  :: PPRSFLX! solid  precip flux
!
END SUBROUTINE CONVECT_CLOSURE
!
END INTERFACE
!
END MODULE MODI_CONVECT_CLOSURE
!    #########################################################################
     SUBROUTINE CONVECT_CLOSURE( KLON, KLEV,                                 &
                                 PPRES, PDPRES, PZ, PDXDY, PLMASS,           &
                                 PTHL, PTH, PRW, PRC, PRI, OTRIG1,           &
                                 PTHC, PRWC, PRCC, PRIC, PWSUB,              &
                                 KLCL, KDPL, KPBL, KLFS, KCTL, KML,          &
                                 PUMF, PUER, PUDR, PUTHL, PURW,              &
                                 PURC, PURI, PUPR,                           &
                                 PDMF, PDER, PDDR, PDTHL, PDRW,              &
                                 PTPR, PSPR, PDTEVR,                         &
                                 PCAPE, PTIMEC,                              &
                                 KFTSTEPS,                                   &
                                 PDTEVRF, PPRLFLX, PPRSFLX                   )
!    #########################################################################
!
!!**** Uses modified Fritsch-Chappell closure
!!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to determine the final adjusted 
!!     (over a time step PTIMEC) environmental values of THETA_l, R_w, R_c, R_i
!!      The final convective tendencies can then be evaluated in the main
!!      routine DEEP_CONVECT by (PTHC-PTH)/PTIMEC
!!
!!
!!**  METHOD
!!    ------
!!      Computations are done at every model level starting from bottom.
!!      The use of masks allows to optimise the inner loops (horizontal loops).
!!      
!!     
!!
!!    EXTERNAL
!!    --------
!!     
!!    CONVECT_CLOSURE_THRVLCL
!!    CONVECT_CLOSURE_ADJUST
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST
!!          XG                 ! gravity constant
!!          XP00               ! reference pressure
!!          XRD, XRV           ! gaz  constants for dry air and water vapor
!!          XCPD, XCPV         ! specific heat for dry air and water vapor
!!          XCL, XCI           ! specific heat for liquid water and ice
!!          XTT                ! triple point temperature
!!          XLVTT, XLSTT       ! vaporization, sublimation heat constant
!!
!!      Module MODD_CONVPAR
!!          XA25               ! reference grid area
!!          XSTABT             ! stability factor in time integration 
!!          XSTABC             ! stability factor in CAPE adjustment
!!          XMELDPTH           ! allow melting over specific pressure depth
!!
!!     Module MODD_CONVPAREXT
!!          JCVEXB, JCVEXT     ! extra levels on the vertical boundaries
!!
!!
!!    REFERENCE
!!    ---------
!!
!!      Book1,2 of documentation ( routine CONVECT_CLOSURE)
!!      Fritsch and Chappell, 1980, J. Atmos. Sci.
!!      Kain and Fritsch, 1993, Meteor. Monographs, Vol.
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    26/03/96 
!!   Peter Bechtold 04/10/97 change for enthalpie, r_c + r_i tendencies
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
USE MODD_CONVPAR
USE MODD_CONVPAREXT
!
USE MODI_CONVECT_SATMIXRATIO
USE MODI_CONVECT_CLOSURE_THRVLCL
USE MODI_CONVECT_CLOSURE_ADJUST
!
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER,                   INTENT(IN) :: KLON   ! horizontal dimension
INTEGER,                   INTENT(IN) :: KLEV   ! vertical dimension
INTEGER, DIMENSION(KLON),  INTENT(IN) :: KLFS   ! index for level of free sink
INTEGER, DIMENSION(KLON),  INTENT(IN) :: KLCL   ! index lifting condens. level
INTEGER, DIMENSION(KLON),  INTENT(IN) :: KCTL   ! index for cloud top level
INTEGER, DIMENSION(KLON),  INTENT(IN) :: KDPL   ! index for departure level 
INTEGER, DIMENSION(KLON),  INTENT(IN) :: KPBL   ! index for top of source layer
INTEGER, DIMENSION(KLON),  INTENT(IN) :: KML    ! index for melting level
REAL, DIMENSION(KLON),  INTENT(INOUT) :: PTIMEC ! convection time step 
REAL, DIMENSION(KLON),     INTENT(IN) :: PDXDY  ! grid area (m^2)
REAL, DIMENSION(KLON,KLEV),INTENT(IN) :: PTHL   ! grid scale enthalpy (J/kg)
REAL, DIMENSION(KLON,KLEV),INTENT(IN) :: PTH    ! grid scale theta        
REAL, DIMENSION(KLON,KLEV),INTENT(IN) :: PRW    ! grid scale total water  
			                        ! mixing ratio 
REAL, DIMENSION(KLON,KLEV),INTENT(IN) :: PRC    ! grid scale r_c 
REAL, DIMENSION(KLON,KLEV),INTENT(IN) :: PRI    ! grid scale r_i 
LOGICAL, DIMENSION(KLON),  INTENT(IN) :: OTRIG1 ! logical to keep trace of 
                                                ! convective arrays modified in UPDRAFT
!   
!
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PPRES  ! pressure (P)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PDPRES ! pressure difference between 
                                                 ! bottom and top of layer (Pa)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PLMASS ! mass of model layer (kg)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PZ     ! height of model layer (m) 
REAL, DIMENSION(KLON),     INTENT(IN)  :: PCAPE  ! available potent. energy
INTEGER,                INTENT(OUT)   :: KFTSTEPS! maximum of fract time steps
                                                 ! only used for chemical tracers
!
!
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PUMF  ! updraft mass flux (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PUER  ! updraft entrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PUDR  ! updraft detrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PUPR  ! updraft precipitation in
                                                  ! flux units (kg water / s)
REAL, DIMENSION(KLON,KLEV), INTENT(IN)  :: PUTHL  ! updraft enthalpy (J/kg)
REAL, DIMENSION(KLON,KLEV), INTENT(IN)  :: PURW   ! updraft total water (kg/kg)
REAL, DIMENSION(KLON,KLEV), INTENT(IN)  :: PURC   ! updraft cloud water (kg/kg)
REAL, DIMENSION(KLON,KLEV), INTENT(IN)  :: PURI   ! updraft cloud ice   (kg/kg)
!
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PDMF  ! downdraft mass flux (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PDER  ! downdraft entrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PDDR  ! downdraft detrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(IN)   :: PDTHL ! downdraft enthalpy (J/kg)
REAL, DIMENSION(KLON,KLEV), INTENT(IN)   :: PDRW  ! downdraft total water (kg/kg)
REAL, DIMENSION(KLON),      INTENT(INOUT):: PTPR  ! total surf precipitation (kg/s)
REAL, DIMENSION(KLON),      INTENT(OUT)  :: PSPR  ! solid surf precipitation (kg/s)
REAL, DIMENSION(KLON),      INTENT(INOUT):: PDTEVR! donwndraft evapor. (kg/s)
!
REAL, DIMENSION(KLON,KLEV), INTENT(OUT)  :: PTHC  ! conv. adj. grid scale theta
REAL, DIMENSION(KLON,KLEV), INTENT(OUT)  :: PRWC  ! conv. adj. grid scale r_w 
REAL, DIMENSION(KLON,KLEV), INTENT(OUT)  :: PRCC  ! conv. adj. grid scale r_c 
REAL, DIMENSION(KLON,KLEV), INTENT(OUT)  :: PRIC  ! conv. adj. grid scale r_i 
REAL, DIMENSION(KLON,KLEV), INTENT(OUT)  :: PWSUB ! envir. compensating subsidence(Pa/s)
!
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PDTEVRF! downdraft evaporation rate
REAL, DIMENSION(KLON,KLEV), INTENT(OUT)  :: PPRLFLX! liquid precip flux
REAL, DIMENSION(KLON,KLEV), INTENT(OUT)  :: PPRSFLX! solid  precip flux
!
!*       0.2   Declarations of local variables :
!
INTEGER :: IIE, IKB, IKE  ! horizontal + vertical loop bounds
INTEGER :: IKS            ! vertical dimension
INTEGER :: JK, JKP, JKMAX ! vertical loop index
INTEGER :: JI             ! horizontal loop index
INTEGER :: JITER          ! iteration loop index
INTEGER :: JSTEP          ! fractional time loop index
REAL    :: ZCPORD, ZRDOCP ! C_pd / R_d, R_d / C_pd
!
REAL, DIMENSION(KLON,KLEV) :: ZTHLC       ! convectively adjusted 
                                          ! grid scale enthalpy
REAL, DIMENSION(KLON,KLEV) :: ZOMG        ! conv. environm. subsidence (Pa/s)
REAL, DIMENSION(KLON,KLEV) :: ZUMF        ! non-adjusted updraft mass flux
REAL, DIMENSION(KLON,KLEV) :: ZUER        !   "     updraft  entrainm. rate
REAL, DIMENSION(KLON,KLEV) :: ZUDR        !   "     updraft  detrainm. rate
REAL, DIMENSION(KLON,KLEV) :: ZDMF        !   "   downdraft mass flux
REAL, DIMENSION(KLON,KLEV) :: ZDER        !   "   downdraft  entrainm. rate
REAL, DIMENSION(KLON,KLEV) :: ZDDR        !   "   downdraft  detrainm. rate
REAL, DIMENSION(KLON)     :: ZTPR         !   "   total precipitation
REAL, DIMENSION(KLON)     :: ZDTEVR       !   "   total downdraft evapor. 
REAL, DIMENSION(KLON,KLEV):: ZPRLFLX      !   "   liquid precip flux
REAL, DIMENSION(KLON,KLEV):: ZPRSFLX      !   "   solid  precip flux
REAL, DIMENSION(KLON)     :: ZPRMELT      ! melting of precipitation
REAL, DIMENSION(KLON)     :: ZPRMELTO     ! non-adjusted  "
REAL, DIMENSION(KLON)     :: ZADJ         ! mass adjustment factor
REAL, DIMENSION(KLON)     :: ZADJMAX      ! limit value for ZADJ
REAL, DIMENSION(KLON)     :: ZCAPE        ! new CAPE after adjustment
REAL, DIMENSION(KLON)     :: ZTIMEC       ! fractional convective time step
REAL, DIMENSION(KLON,KLEV):: ZTIMC        ! 2D work array for ZTIMEC
!
REAL, DIMENSION(KLON)     :: ZTHLCL       ! new  theta at LCL
REAL, DIMENSION(KLON)     :: ZRVLCL       ! new  r_v at LCL
REAL, DIMENSION(KLON)     :: ZZLCL        ! height of LCL
REAL, DIMENSION(KLON)     :: ZTLCL        ! temperature at LCL
REAL, DIMENSION(KLON)     :: ZTELCL       ! envir. temper. at LCL
REAL, DIMENSION(KLON)     :: ZTHEUL       ! theta_e for undilute ascent
REAL, DIMENSION(KLON)     :: ZTHES1, ZTHES2! saturation environm. theta_e
REAL, DIMENSION(KLON,KLEV) :: ZTHMFIN, ZTHMFOUT, ZRWMFIN, ZRWMFOUT
REAL, DIMENSION(KLON,KLEV) :: ZRCMFIN, ZRCMFOUT, ZRIMFIN, ZRIMFOUT
                                    ! work arrays for environm. compensat. mass flux
REAL, DIMENSION(KLON)     :: ZPI          ! (P/P00)**R_d/C_pd 
REAL, DIMENSION(KLON)     :: ZLV          ! latent heat of vaporisation
REAL, DIMENSION(KLON)     :: ZLS          ! latent heat of sublimation 
REAL, DIMENSION(KLON)     :: ZLM          ! latent heat of melting
REAL, DIMENSION(KLON)     :: ZCPH         ! specific heat C_ph
REAL, DIMENSION(KLON)     :: ZMELDPTH     ! actual depth of melting layer 
INTEGER, DIMENSION(KLON)  :: ITSTEP       ! fractional convective time step
INTEGER, DIMENSION(KLON)  :: ICOUNT       ! timestep counter 
INTEGER, DIMENSION(KLON)  :: ILCL         ! index lifting condens. level
INTEGER, DIMENSION(KLON)  :: IWORK1       ! work array
REAL, DIMENSION(KLON)     :: ZWORK1, ZWORK2, ZWORK3, ZWORK4, ZWORK5
REAL, DIMENSION(KLON,KLEV):: ZWORK6
LOGICAL, DIMENSION(KLON)  :: GWORK1, GWORK3! work arrays
LOGICAL, DIMENSION(KLON,KLEV) :: GWORK4    ! work array
!
!
!-------------------------------------------------------------------------------
!
!*       0.2    Initialize  local variables
!               ----------------------------
!
!
PSPR(:)   = 0.
ZTIMC(:,:) = 0.
ZTHES2(:) = 0.
ZWORK1(:) = 0. 
ZWORK2(:) = 0. 
ZWORK3(:) = 0. 
ZWORK4(:) = 0. 
ZWORK5(:) = 0. 
GWORK1(:) = .FALSE.
GWORK3(:) = .FALSE.  
GWORK4(:,:) = .FALSE.  
ILCL(:)   = KLCL(:)
!
ZCPORD    = XCPD / XRD
ZRDOCP    = XRD / XCPD
!
ZADJ(:)   = 1.
ZWORK5(:) = 1.
WHERE( .NOT. OTRIG1(:) ) ZWORK5(:) = 0. 
!
!
!*       0.3   Compute loop bounds
!              ------------------- 
!
IIE    = KLON
IKB    = 1 + JCVEXB 
IKS    = KLEV
IKE    = KLEV - JCVEXT 
JKMAX  = MAXVAL( KCTL(:) )
!
!
!*       2.     Save initial mass flux values to be used in adjustment procedure
!               ---------------------------------------------------------------
!
ZUMF(:,:)  = PUMF(:,:)
ZUER(:,:)  = PUER(:,:)
ZUDR(:,:)  = PUDR(:,:)
ZDMF(:,:)  = PDMF(:,:)
ZDER(:,:)  = PDER(:,:)
ZDDR(:,:)  = PDDR(:,:)
ZTPR(:)    = PTPR(:)
ZDTEVR(:)  = PDTEVR(:)
ZOMG(:,:)  = 0.
PWSUB(:,:) = 0. 
ZPRMELT(:) = 0.
PPRLFLX(:,:) = 0.
ZPRLFLX(:,:) = 0.
PPRSFLX(:,:) = 0.
ZPRSFLX(:,:) = 0.
!
!
!*       2.1    Some preliminar computations for melting of precipitation
!               used later in section 9 and computation of precip fluxes
!               Precipitation fluxes are updated for melting and evaporation
!               ---------------------------------------------------------
!
!
ZWORK1(:) = 0.
ZMELDPTH(:) = 0.
ZWORK6(:,:) = 0.
DO JK = JKMAX + 1, IKB + 1, -1
   ! Nota: PUPR is total precipitation flux, but the solid, liquid
   !       precipitation is stored in units kg/kg; therefore we compute here
   !       the solid fraction of the total precipitation flux.
  DO JI = 1, IIE
     ZWORK2(JI)    = PUPR(JI,JK) / ( PURC(JI,JK) + PURI(JI,JK) + 1.E-8 )
     ZPRMELT(JI)   = ZPRMELT(JI) + PURI(JI,JK) * ZWORK2(JI)
     ZWORK1(JI)    = ZWORK1(JI) + PURC(JI,JK) * ZWORK2(JI) - PDTEVRF(JI,JK)
     ZPRLFLX(JI,JK)= MAX( 0., ZWORK1(JI) )
     ZPRMELT(JI)   = ZPRMELT(JI) + MIN( 0., ZWORK1(JI) )
     ZPRSFLX(JI,JK)= ZPRMELT(JI) 
     IF ( KML(JI) >= JK .AND. ZMELDPTH(JI) <= XMELDPTH ) THEN                 
          ZPI(JI)    = ( PPRES(JI,JK) / XP00 ) ** ZRDOCP 
          ZWORK3(JI) = PTH(JI,JK) * ZPI(JI)            ! temperature estimate
          ZLM(JI)    = XLSTT + ( XCPV - XCI ) * ( ZWORK3(JI) - XTT ) -       &
               ( XLVTT + ( XCPV - XCL ) * ( ZWORK3(JI) - XTT ) ) ! L_s - L_v
          ZCPH(JI)   = XCPD + XCPV * PRW(JI,JK)
          ZMELDPTH(JI) = ZMELDPTH(JI) + PDPRES(JI,JK)
          ZWORK6(JI,JK)= ZLM(JI) * PTIMEC(JI) / PLMASS(JI,JK) * PDPRES(JI,JK)
          ZOMG(JI,JK)= 1. ! at this place only used as work variable
     END IF
  END DO
!
END DO
!
ZWORK2(:) = 0.
DO JK = JKMAX, IKB + 1, -1
    ZWORK1(:) = ZPRMELT(:) * PDPRES(:,JK) / MAX( XMELDPTH, ZMELDPTH(:) )
    ZWORK2(:) = ZWORK2(:) + ZWORK1(:) * ZOMG(:,JK)
    ZPRLFLX(:,JK) = ZPRLFLX(:,JK) + ZWORK2(:) 
    ZPRSFLX(:,JK) = ZPRSFLX(:,JK) - ZWORK2(:)
END DO 
WHERE( ZPRSFLX(:,:) < 1. ) ZPRSFLX(:,:)=0.
ZPRMELTO(:) = ZPRMELT(:)
!
!
!*       3.     Compute limits on the closure adjustment factor so that the
!               inflow in convective drafts from a given layer can't be larger 
!               than the mass contained in this layer initially.
!               ---------------------------------------------------------------
!
ZADJMAX(:) = 1000.
IWORK1(:) = MAX( ILCL(:), KLFS(:) )
JKP = MINVAL( KDPL(:) )
DO JK = JKP, IKE
  DO JI = 1, IIE
    IF( JK > KDPL(JI) .AND. JK <= IWORK1(JI) ) THEN
        ZWORK1(JI)  = PLMASS(JI,JK) /                                      &
                  ( ( PUER(JI,JK) + PDER(JI,JK) + 1.E-5 ) * PTIMEC(JI) )
        ZADJMAX(JI) = MIN( ZADJMAX(JI), ZWORK1(JI) )
    END IF
  END DO
END DO
!
!
GWORK1(:) = OTRIG1(:)  ! logical array to limit adjustment to not definitively
                       ! adjusted columns
!
DO JK = IKB, IKE
  ZTHLC(:,JK) = PTHL(:,JK) ! initialize adjusted envir. values 
  PRWC(:,JK)  = PRW(:,JK)
  PRCC(:,JK)  = PRC(:,JK)
  PRIC(:,JK)  = PRI(:,JK)
  PTHC(:,JK)  = PTH(:,JK)
END DO
!
!
!
DO JITER = 1, 6  ! Enter adjustment loop to assure that all CAPE is
                 ! removed within the advective time interval TIMEC
!
     ZTIMEC(:) = PTIMEC(:)
     GWORK4(:,:)   = SPREAD( GWORK1(:), DIM=2, NCOPIES=IKS )
     WHERE( GWORK4(:,:) ) PWSUB(:,:) = 0.
     ZOMG(:,:)=0.
!
     DO JK = IKB + 1, JKMAX 
           JKP = MAX( IKB + 1, JK - 1 )
           WHERE ( GWORK1(:) .AND. JK <= KCTL(:) )
!
!
!*       4.     Determine vertical velocity at top and bottom of each layer
!               to satisfy mass continuity.
!               ---------------------------------------------------------------
              ! we compute here Domega/Dp = - g rho Dw/Dz = 1/Dt
!
             ZWORK1(:)   = - ( PUER(:,JKP) + PDER(:,JKP) -                   &
                           PUDR(:,JKP) - PDDR(:,JKP) ) / PLMASS(:,JKP)
!    
             PWSUB(:,JK) = PWSUB(:,JKP) - PDPRES(:,JK-1) * ZWORK1(:)
              ! we use PDPRES(JK-1) and not JKP in order to have zero subsidence
              ! at the first layer
!
!   
!*       5.     Compute fractional time step. For stability or 
!               mass conservation reasons one must split full time step PTIMEC)
!               ---------------------------------------------------------------
!
             ZWORK1(:) = XSTABT * PDPRES(:,JKP) / ( ABS( PWSUB(:,JK) ) + 1.E-10 )
              ! the factor XSTABT is used for stability reasons
             ZTIMEC(:) = MIN( ZTIMEC(:), ZWORK1(:) ) 
!
              ! transform vertical velocity in mass flux units
             ZOMG(:,JK) = PWSUB(:,JK) * PDXDY(:) / XG 
         END WHERE
     END DO
!
!
     WHERE( GWORK4(:,:) )
           ZTHLC(:,:) = PTHL(:,:) ! reinitialize adjusted envir. values 
           PRWC(:,:)  = PRW(:,:)  ! when iteration criterium not attained
           PRCC(:,:)  = PRC(:,:)
           PRIC(:,:)  = PRI(:,:)
           PTHC(:,:)  = PTH(:,:)
     END WHERE
!
! 
!        6. Check for mass conservation, i.e. ZWORK1 > 1.E-2
!           If mass is not conserved, the convective tendencies
!           automatically become zero.
!           ----------------------------------------------------
!
    DO JI = 1, IIE
       JK=KCTL(JI)
       ZWORK1(JI) = PUDR(JI,JK) * PDPRES(JI,JK) / ( PLMASS(JI,JK) + .1 )    &
                                                            - PWSUB(JI,JK)
    END DO
    WHERE( GWORK1(:) .AND. ABS( ZWORK1(:) ) - .01 > 0. )
        GWORK1(:) = .FALSE.
        PTIMEC(:) = 1.E-1
        ZTPR(:)   = 0.
        ZWORK5(:) = 0.
    END WHERE
    DO JK = IKB, IKE
        PWSUB(:,JK) = PWSUB(:,JK) * ZWORK5(:)
        ZPRLFLX(:,JK) = ZPRLFLX(:,JK) * ZWORK5(:)
        ZPRSFLX(:,JK) = ZPRSFLX(:,JK) * ZWORK5(:)
    END DO
    GWORK4(:,1:IKB) = .FALSE.
    GWORK4(:,IKE:IKS)   = .FALSE.
!
    ITSTEP(:) = INT( PTIMEC(:) / ZTIMEC(:) ) + 1 
    ZTIMEC(:) = PTIMEC(:) / REAL( ITSTEP(:) ) ! adjust  fractional time step
                                           ! to be an integer multiple of PTIMEC
    ZTIMC(:,:)= SPREAD( ZTIMEC(:), DIM=2, NCOPIES=IKS )
    ICOUNT(:) = 0
!
!
!
    KFTSTEPS = MAXVAL( ITSTEP(:) )
    DO JSTEP = 1, KFTSTEPS ! Enter the fractional time step loop here
!
     	 ICOUNT(:) = ICOUNT(:) + 1
!
	     GWORK3(:) =  ITSTEP(:) >= ICOUNT(:) .AND. GWORK1(:) 
!
!
!*       7.     Assign enthalpy and r_w values at the top and bottom of each
!               layer based on the sign of w
!               ------------------------------------------------------------
!
             ZTHMFIN(:,:)   = 0.
             ZRWMFIN(:,:)   = 0.
             ZRCMFIN(:,:)   = 0.
             ZRIMFIN(:,:)   = 0.
             ZTHMFOUT(:,:)  = 0.
             ZRWMFOUT(:,:)  = 0.
             ZRCMFOUT(:,:)  = 0.
             ZRIMFOUT(:,:)  = 0.
!
         DO JK = IKB + 1, JKMAX
           DO JI = 1, IIE
              GWORK4(JI,JK) = GWORK3(JI) .AND. JK <= KCTL(JI)
           END DO
           JKP = MAX( IKB + 1, JK - 1 )
           DO JI = 1, IIE
           IF ( GWORK3(JI) ) THEN
!
               ZWORK1(JI)       = SIGN( 1., ZOMG(JI,JK) )
               ZWORK2(JI)       = 0.5 * ( 1. + ZWORK1(JI) )
               ZWORK1(JI)       = 0.5 * ( 1. - ZWORK1(JI) )
               ZTHMFIN(JI,JK)   = - ZOMG(JI,JK) * ZTHLC(JI,JKP) * ZWORK1(JI)
               ZTHMFOUT(JI,JK)  =   ZOMG(JI,JK) * ZTHLC(JI,JK)  * ZWORK2(JI)
               ZRWMFIN(JI,JK)   = - ZOMG(JI,JK) * PRWC(JI,JKP) * ZWORK1(JI)
               ZRWMFOUT(JI,JK)  =   ZOMG(JI,JK) * PRWC(JI,JK)  * ZWORK2(JI)
               ZRCMFIN(JI,JK)   = - ZOMG(JI,JK) * PRCC(JI,JKP) * ZWORK1(JI)
               ZRCMFOUT(JI,JK)  =   ZOMG(JI,JK) * PRCC(JI,JK)  * ZWORK2(JI)
               ZRIMFIN(JI,JK)   = - ZOMG(JI,JK) * PRIC(JI,JKP) * ZWORK1(JI)
               ZRIMFOUT(JI,JK)  =   ZOMG(JI,JK) * PRIC(JI,JK)  * ZWORK2(JI)
           END IF
           END DO
           DO JI = 1, IIE
           IF ( GWORK3(JI) ) THEN
               ZTHMFIN(JI,JKP)  = ZTHMFIN(JI,JKP)  + ZTHMFOUT(JI,JK) * ZWORK2(JI)
               ZTHMFOUT(JI,JKP) = ZTHMFOUT(JI,JKP) + ZTHMFIN(JI,JK)  * ZWORK1(JI)
               ZRWMFIN(JI,JKP)  = ZRWMFIN(JI,JKP)  + ZRWMFOUT(JI,JK) * ZWORK2(JI)
               ZRWMFOUT(JI,JKP) = ZRWMFOUT(JI,JKP) + ZRWMFIN(JI,JK)  * ZWORK1(JI)
               ZRCMFIN(JI,JKP)  = ZRCMFIN(JI,JKP)  + ZRCMFOUT(JI,JK) * ZWORK2(JI)
               ZRCMFOUT(JI,JKP) = ZRCMFOUT(JI,JKP) + ZRCMFIN(JI,JK)  * ZWORK1(JI)
               ZRIMFIN(JI,JKP)  = ZRIMFIN(JI,JKP)  + ZRIMFOUT(JI,JK) * ZWORK2(JI)
               ZRIMFOUT(JI,JKP) = ZRIMFOUT(JI,JKP) + ZRIMFIN(JI,JK)  * ZWORK1(JI)
!
           END IF
           END DO
         END DO
!
         WHERE ( GWORK4(:,:) )
!
!******************************************************************************
!
!*       8.     Update the environmental values of enthalpy and r_w at each level
!               NOTA: These are the MAIN EQUATIONS of the scheme
!               -----------------------------------------------------------------
!
!
           ZTHLC(:,:) = ZTHLC(:,:) + ZTIMC(:,:) / PLMASS(:,:) * (      &
                          ZTHMFIN(:,:) + PUDR(:,:) * PUTHL(:,:)  +     &
                          PDDR(:,:) * PDTHL(:,:) - ZTHMFOUT(:,:) -     &
                        ( PUER(:,:) + PDER(:,:) ) * PTHL(:,:)   )
           PRWC(:,:)  = PRWC(:,:) + ZTIMC(:,:) / PLMASS(:,:) *  (      &
                         ZRWMFIN(:,:) + PUDR(:,:) * PURW(:,:)  +       &
                         PDDR(:,:) * PDRW(:,:) - ZRWMFOUT(:,:) -       &
                        ( PUER(:,:) + PDER(:,:) ) * PRW(:,:)    )    
           PRCC(:,:)  = PRCC(:,:) + ZTIMC(:,:) / PLMASS(:,:) *  (      &
               ZRCMFIN(:,:) + PUDR(:,:) * PURC(:,:) - ZRCMFOUT(:,:) -  &
                        ( PUER(:,:) + PDER(:,:) ) * PRC(:,:)    )    
           PRIC(:,:)  = PRIC(:,:) + ZTIMC(:,:) / PLMASS(:,:) *  (      &
               ZRIMFIN(:,:) + PUDR(:,:) * PURI(:,:) - ZRIMFOUT(:,:) -  & 
                        ( PUER(:,:) + PDER(:,:) ) * PRI(:,:)    )    
!
!
!******************************************************************************
!
         END WHERE
!
    END DO ! Exit the fractional time step loop
!
! 
!*           9.    Allow frozen precipitation to melt over a 200 mb deep layer
!                  -----------------------------------------------------------
!
      DO JK = JKMAX, IKB + 1, -1
            ZTHLC(:,JK) = ZTHLC(:,JK) -                                &
               ZPRMELT(:) * ZWORK6(:,JK) / MAX( XMELDPTH, ZMELDPTH(:) )
      END DO
!
!
!*          10.    Compute final linearized value of theta envir.
!                  ----------------------------------------------
!
      DO JK = IKB + 1, JKMAX
         DO JI = 1, IIE
         IF( GWORK1(JI) .AND. JK <= KCTL(JI) ) THEN
           ZPI(JI)    = ( XP00 / PPRES(JI,JK) ) ** ZRDOCP
           ZCPH(JI)   = XCPD + PRWC(JI,JK) * XCPV
           ZWORK2(JI) = PTH(JI,JK) / ZPI(JI)  ! first temperature estimate
           ZLV(JI)    = XLVTT + ( XCPV - XCL ) * ( ZWORK2(JI) - XTT )
           ZLS(JI)    = XLVTT + ( XCPV - XCI ) * ( ZWORK2(JI) - XTT )
             ! final linearized temperature
           ZWORK2(JI) = ( ZTHLC(JI,JK) + ZLV(JI) * PRCC(JI,JK) + ZLS(JI) * PRIC(JI,JK) &
                       - (1. + PRWC(JI,JK) ) * XG * PZ(JI,JK) ) / ZCPH(JI)
           ZWORK2(JI) = MAX( 180., MIN( 340., ZWORK2(JI) ) )
           PTHC(JI,JK)= ZWORK2(JI) * ZPI(JI) ! final adjusted envir. theta
         END IF
         END DO
      END DO
!
!
!*         11.     Compute new cloud ( properties at new LCL )
!                     NOTA: The computations are very close to
!                           that in routine TRIGGER_FUNCT
!                  ---------------------------------------------
!
      CALL CONVECT_CLOSURE_THRVLCL(  KLON, KLEV,                           &
                                     PPRES, PTHC, PRWC, PZ, GWORK1,        &
                                     ZTHLCL, ZRVLCL, ZZLCL, ZTLCL, ZTELCL, &
                                     ILCL, KDPL, KPBL )
!
!
       ZTLCL(:)  = MAX( 230., MIN( 335., ZTLCL(:) ) )  ! set some overflow bounds
       ZTELCL(:) = MAX( 230., MIN( 335., ZTELCL(:) ) )
       ZTHLCL(:) = MAX( 230., MIN( 345., ZTHLCL(:) ) )
       ZRVLCL(:) = MAX(   0., MIN(   1., ZRVLCL(:) ) )
!
!
!*         12.    Compute adjusted CAPE
!                 ---------------------
!
       ZCAPE(:)  = 0.
       ZPI(:)    = ZTHLCL(:) / ZTLCL(:)
       ZPI(:)    = MAX( 0.95, MIN( 1.5, ZPI(:) ) )
       ZWORK1(:) = XP00 / ZPI(:) ** ZCPORD ! pressure at LCL
!
       CALL CONVECT_SATMIXRATIO( KLON, ZWORK1, ZTELCL, ZWORK3, ZLV, ZLS, ZCPH )
       ZWORK3(:) = MIN(   .1, MAX(   0., ZWORK3(:) ) )
!
	        ! compute theta_e updraft undilute
       ZTHEUL(:) = ZTLCL(:) * ZPI(:) ** ( 1. - 0.28 * ZRVLCL(:) )            &
                                  * EXP( ( 3374.6525 / ZTLCL(:) - 2.5403 )   &
                                  * ZRVLCL(:) * ( 1. + 0.81 * ZRVLCL(:) ) )
!
	        ! compute theta_e saturated environment at LCL
       ZTHES1(:) = ZTELCL(:) * ZPI(:) ** ( 1. - 0.28 * ZWORK3(:) )           &
                                  * EXP( ( 3374.6525 / ZTELCL(:) - 2.5403 )  &
                                  * ZWORK3(:) * ( 1. + 0.81 * ZWORK3(:) ) )
!
      DO JK = MINVAL( ILCL(:) ), JKMAX
        JKP = JK - 1
        DO JI = 1, IIE
          ZWORK4(JI) = 1.
          IF ( JK == ILCL(JI) ) ZWORK4(JI) = 0.
!
           ! compute theta_e saturated environment and adjusted values
           ! of theta
!
          GWORK3(JI)  = JK >= ILCL(JI) .AND. JK <= KCTL(JI) .AND. GWORK1(JI) 
!
          ZPI(JI)     = ( XP00 / PPRES(JI,JK) ) ** ZRDOCP
          ZWORK2(JI)  = PTHC(JI,JK) / ZPI(JI)
        END DO
!
        CALL CONVECT_SATMIXRATIO( KLON, PPRES(:,JK), ZWORK2, ZWORK3, ZLV, ZLS, ZCPH )
!
!
        DO JI = 1, IIE
          IF ( GWORK3(JI) ) THEN
              ZTHES2(JI)  = ZWORK2(JI) * ZPI(JI) ** ( 1. - 0.28 * ZWORK3(JI) )   &
                                   * EXP( ( 3374.6525 / ZWORK2(JI) - 2.5403 ) &
                                   * ZWORK3(JI) * ( 1. + 0.81 * ZWORK3(JI) ) )
!
              ZWORK3(JI)  = PZ(JI,JK) - PZ(JI,JKP) * ZWORK4(JI) -                &
                           ( 1. - ZWORK4(JI) ) * ZZLCL(JI)    ! level thickness
              ZWORK1(JI)  = ( 2. * ZTHEUL(JI) ) / ( ZTHES1(JI) + ZTHES2(JI) ) - 1.
              ZCAPE(JI)   = ZCAPE(JI) + XG * ZWORK3(JI) * MAX( 0., ZWORK1(JI) )
              ZTHES1(JI)  = ZTHES2(JI)
          END IF
        END DO
      END DO
!
!                                                          
!*         13.     Determine mass adjustment factor knowing how much
!                  CAPE has been removed.
!                  -------------------------------------------------
!
       WHERE ( GWORK1(:) )
           ZWORK1(:) = MAX( PCAPE(:) - ZCAPE(:), 0.1 * PCAPE(:) )
           ZWORK2(:) = ZCAPE(:) / ( PCAPE(:) + 1.E-8 )
!       
           GWORK1(:) = ZWORK2(:) > 0.1 .OR. ZCAPE(:) == 0. ! mask for adjustment
       END WHERE
!
       WHERE ( ZCAPE(:) == 0. .AND. GWORK1(:) )  ZADJ(:) = ZADJ(:) * 0.5
       WHERE ( ZCAPE(:) /= 0. .AND. GWORK1(:) )                              &
               ZADJ(:) = ZADJ(:) * XSTABC * PCAPE(:) / ( ZWORK1(:) + 1.E-8 )
       ZADJ(:) = MIN( ZADJ(:), ZADJMAX(:) )  
!
!
!*         13.     Adjust mass flux by the factor ZADJ to converge to
!                  specified degree of stabilization
!                 ----------------------------------------------------
!
       CALL CONVECT_CLOSURE_ADJUST( KLON, KLEV, ZADJ,                     &
                                    PUMF, ZUMF, PUER, ZUER, PUDR, ZUDR,   &
                                    PDMF, ZDMF, PDER, ZDER, PDDR, ZDDR,   &
                                    ZPRMELT, ZPRMELTO, PDTEVR, ZDTEVR,    &
                                    PTPR, ZTPR,                           &
                                    PPRLFLX, ZPRLFLX, PPRSFLX, ZPRSFLX    )
!
!
      IF ( COUNT( GWORK1(:) ) == 0 ) EXIT ! exit big adjustment iteration loop
                                          ! when all columns have reached 
                                          ! desired degree of stabilization.
!
END DO  ! end of big adjustment iteration loop
!
!
        ! skip adj. total water array  to water vapor
DO JK = IKB, IKE
  PRWC(:,JK) = MAX( 0., PRWC(:,JK) - PRCC(:,JK) - PRIC(:,JK) )
END DO
!
        ! compute surface solid (ice) precipitation 
PSPR(:) = ZPRMELT(:) * ( 1. - ZMELDPTH(:) / XMELDPTH )
PSPR(:) = MAX( 0., PSPR(:) )
!
!
END SUBROUTINE CONVECT_CLOSURE
