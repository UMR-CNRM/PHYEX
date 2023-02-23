!MNH_LIC Copyright 1995-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     #################
      MODULE MODI_CONVECT_UPDRAFT
!     #################
!
INTERFACE
!
       SUBROUTINE CONVECT_UPDRAFT( KLON, KLEV,                                      &
                              KICE, PPRES, PDPRES, PZ, PTHL, PTHV, PTHES, PRW, &
                              PTHLCL, PTLCL, PRVLCL, PWLCL, PZLCL, PTHVELCL,   &
                              PMFLCL, OTRIG, KLCL, KDPL, KPBL,                 &
                              PUMF, PUER, PUDR, PUTHL, PUTHV, PURW,            &
                              PURC, PURI, PURR, PURS, PUPR,                    &
                              PUTPR, PCAPE, KCTL, KETL, PUTT )
!
INTEGER, INTENT(IN)                    :: KLON  ! horizontal dimension
INTEGER, INTENT(IN)                    :: KLEV  ! vertical dimension
INTEGER, INTENT(IN)                    :: KICE  ! flag for ice ( 1 = yes,
                                                !                0 = no ice )
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PTHL  ! grid scale enthalpy (J/kg)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PTHV  ! grid scale theta_v     
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PTHES ! grid scale saturated theta_e 
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PRW   ! grid scale total water  
                                                ! mixing ratio 
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PPRES ! pressure (P)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PDPRES! pressure difference between 
                                                ! bottom and top of layer (Pa) 
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PZ    ! height of model layer (m) 
REAL, DIMENSION(KLON),     INTENT(IN) :: PTHLCL ! theta at LCL
REAL, DIMENSION(KLON),     INTENT(IN) :: PTLCL  ! temp. at LCL
REAL, DIMENSION(KLON),     INTENT(IN) :: PRVLCL ! vapor mixing ratio at  LCL
REAL, DIMENSION(KLON),     INTENT(IN) :: PWLCL  ! parcel velocity at LCL (m/s)
REAL, DIMENSION(KLON),     INTENT(IN) :: PMFLCL ! cloud  base unit mass flux
                                                ! (kg/s)
REAL, DIMENSION(KLON),     INTENT(IN) :: PZLCL  ! height at LCL (m)
REAL, DIMENSION(KLON),     INTENT(IN) :: PTHVELCL  ! environm. theta_v at LCL (K)
LOGICAL, DIMENSION(KLON),  INTENT(INOUT):: OTRIG! logical mask for convection 
INTEGER, DIMENSION(KLON),  INTENT(IN) :: KLCL   ! contains vert. index of LCL
INTEGER, DIMENSION(KLON),  INTENT(IN) :: KDPL   ! contains vert. index of DPL 
INTEGER, DIMENSION(KLON),  INTENT(IN) :: KPBL   !  " vert. index of source layertop
!
!
INTEGER, DIMENSION(KLON),  INTENT(OUT):: KCTL   ! contains vert. index of CTL 
INTEGER, DIMENSION(KLON),  INTENT(OUT):: KETL   ! contains vert. index of        &
                                                !equilibrium (zero buoyancy) level 
REAL, DIMENSION(KLON,KLEV), INTENT(OUT):: PUMF  ! updraft mass flux (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(OUT):: PUER  ! updraft entrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(OUT):: PUDR  ! updraft detrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(OUT):: PUTHL ! updraft enthalpy (J/kg)
REAL, DIMENSION(KLON,KLEV), INTENT(OUT):: PUTHV ! updraft theta_v (K)
REAL, DIMENSION(KLON,KLEV), INTENT(OUT):: PUTT  ! updraft temperature(K)
REAL, DIMENSION(KLON,KLEV), INTENT(OUT):: PURW  ! updraft total water (kg/kg)
REAL, DIMENSION(KLON,KLEV), INTENT(OUT):: PURC  ! updraft cloud water (kg/kg)
REAL, DIMENSION(KLON,KLEV), INTENT(OUT):: PURI  ! updraft cloud ice   (kg/kg)
REAL, DIMENSION(KLON,KLEV), INTENT(OUT):: PURR  ! liquid precipit. (kg/kg)
                                                ! produced in  model layer
REAL, DIMENSION(KLON,KLEV),   INTENT(OUT)::PURS ! solid precipit. (kg/kg)
                                                ! produced in  model layer
REAL, DIMENSION(KLON,KLEV),   INTENT(OUT)::PUPR ! updraft precipitation in
                                                ! flux units (kg water / s)
REAL, DIMENSION(KLON),     INTENT(OUT):: PUTPR  ! total updraft precipitation
                                                ! in flux units (kg water / s)
REAL, DIMENSION(KLON),     INTENT(OUT):: PCAPE  ! available potent. energy
!
END SUBROUTINE CONVECT_UPDRAFT
!
END INTERFACE
!
END MODULE MODI_CONVECT_UPDRAFT
!     ##########################################################################
  SUBROUTINE CONVECT_UPDRAFT( KLON, KLEV,                                      &
                              KICE, PPRES, PDPRES, PZ, PTHL, PTHV, PTHES, PRW, &
                              PTHLCL, PTLCL, PRVLCL, PWLCL, PZLCL, PTHVELCL,   &
                              PMFLCL, OTRIG, KLCL, KDPL, KPBL,                 &
                              PUMF, PUER, PUDR, PUTHL, PUTHV, PURW,            &
                              PURC, PURI, PURR, PURS, PUPR,                    &
                              PUTPR, PCAPE, KCTL, KETL, PUTT )
!     ##########################################################################
!
!!**** Compute updraft properties from DPL to CTL. 
!!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to determine updraft properties
!!      ( mass flux, thermodynamics, precipitation ) 
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
!!     Routine CONVECT_MIXING_FUNCT
!!     Routine CONVECT_CONDENS
!!     
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST
!!          XG                 ! gravity constant
!!          XP00               ! reference pressure
!!          XRD, XRV           ! gaz  constants for dry air and water vapor
!!          XCPD, XCPV, XCL    ! Cp of dry air, water vapor and liquid water
!!          XTT                ! triple point temperature
!!          XLVTT              ! vaporisation heat at XTT
!!        
!!
!!      Module MODD_CONVPAR
!!          XA25               ! reference grid area
!!          XCRAD              ! cloud radius
!!          XCDEPTH            ! minimum necessary cloud depth
!!          XENTR              ! entrainment constant
!!          XRCONV             ! constant in precipitation conversion 
!!          XNHGAM             ! coefficient for buoyancy term in w eq.
!!                             ! accounting for nh-pressure
!!          XTFRZ1             ! begin of freezing interval
!!          XTFRZ2             ! begin of freezing interval
!!
!!     Module MODD_CONVPAREXT
!!          JCVEXB, JCVEXT     ! extra levels on the vertical boundaries
!!
!!    REFERENCE
!!    ---------
!!
!!      Book1,2 of documentation ( routine CONVECT_UPDRAFT)
!!      Kain and Fritsch, 1990, J. Atmos. Sci., Vol.
!!      Kain and Fritsch, 1993, Meteor. Monographs, Vol.
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    07/11/95 
!!   Last modified  10/12/97
!!   V.Masson, C.Lac, Sept. 2010 : Correction of a loop for reproducibility
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
USE MODD_CONVPAR
USE MODD_CONVPAREXT
!
USE MODI_CONVECT_CONDENS
USE MODI_CONVECT_MIXING_FUNCT
!
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER, INTENT(IN)                    :: KLON  ! horizontal dimension
INTEGER, INTENT(IN)                    :: KLEV  ! vertical dimension
INTEGER, INTENT(IN)                    :: KICE  ! flag for ice ( 1 = yes,
                                                !                0 = no ice )
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PTHL  ! grid scale enthalpy (J/kg)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PTHV  ! grid scale theta_v     
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PTHES ! grid scale saturated theta_e 
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PRW   ! grid scale total water  
                                                ! mixing ratio 
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PPRES ! pressure (P)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PDPRES! pressure difference between 
                                                ! bottom and top of layer (Pa) 
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PZ    ! height of model layer (m) 
REAL, DIMENSION(KLON),     INTENT(IN) :: PTHLCL ! theta at LCL
REAL, DIMENSION(KLON),     INTENT(IN) :: PTLCL  ! temp. at LCL
REAL, DIMENSION(KLON),     INTENT(IN) :: PRVLCL ! vapor mixing ratio at  LCL
REAL, DIMENSION(KLON),     INTENT(IN) :: PWLCL  ! parcel velocity at LCL (m/s)
REAL, DIMENSION(KLON),     INTENT(IN) :: PMFLCL ! cloud  base unit mass flux
                                                ! (kg/s)
REAL, DIMENSION(KLON),     INTENT(IN) :: PZLCL  ! height at LCL (m)
REAL, DIMENSION(KLON),     INTENT(IN) :: PTHVELCL  ! environm. theta_v at LCL (K)
LOGICAL, DIMENSION(KLON),  INTENT(INOUT):: OTRIG! logical mask for convection 
INTEGER, DIMENSION(KLON),  INTENT(IN) :: KLCL   ! contains vert. index of LCL
INTEGER, DIMENSION(KLON),  INTENT(IN) :: KDPL   ! contains vert. index of DPL 
INTEGER, DIMENSION(KLON),  INTENT(IN) :: KPBL   !  " vert. index of source layertop
!
!
INTEGER, DIMENSION(KLON),  INTENT(OUT):: KCTL   ! contains vert. index of CTL 
INTEGER, DIMENSION(KLON),  INTENT(OUT):: KETL   ! contains vert. index of        &
                                                !equilibrium (zero buoyancy) level 
REAL, DIMENSION(KLON,KLEV), INTENT(OUT):: PUMF  ! updraft mass flux (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(OUT):: PUER  ! updraft entrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(OUT):: PUDR  ! updraft detrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(OUT):: PUTHL ! updraft enthalpy (J/kg)
REAL, DIMENSION(KLON,KLEV), INTENT(OUT):: PUTHV ! updraft theta_v (K)
REAL, DIMENSION(KLON,KLEV), INTENT(OUT):: PUTT  ! updraft temperature(K)
REAL, DIMENSION(KLON,KLEV), INTENT(OUT):: PURW  ! updraft total water (kg/kg)
REAL, DIMENSION(KLON,KLEV), INTENT(OUT):: PURC  ! updraft cloud water (kg/kg)
REAL, DIMENSION(KLON,KLEV), INTENT(OUT):: PURI  ! updraft cloud ice   (kg/kg)
REAL, DIMENSION(KLON,KLEV), INTENT(OUT):: PURR  ! liquid precipit. (kg/kg)
                                                ! produced in  model layer
REAL, DIMENSION(KLON,KLEV),   INTENT(OUT)::PURS ! solid precipit. (kg/kg)
                                                ! produced in  model layer
REAL, DIMENSION(KLON,KLEV),   INTENT(OUT)::PUPR ! updraft precipitation in
                                                ! flux units (kg water / s)
REAL, DIMENSION(KLON),     INTENT(OUT):: PUTPR  ! total updraft precipitation
                                                ! in flux units (kg water / s)
REAL, DIMENSION(KLON),     INTENT(OUT):: PCAPE  ! available potent. energy
!
!*       0.2   Declarations of local variables :
!
INTEGER :: IIE, IKB, IKE  ! horizontal and vertical loop bounds
INTEGER :: JI             ! horizontal loop index
INTEGER :: JK, JKP, JKM, JK1, JK2, JKMIN   ! vertical loop index
REAL    :: ZEPSA          ! R_v / R_d, C_pv / C_pd 
REAL    :: ZRDOCP         ! C_pd / R_d, R_d / C_pd
!
REAL, DIMENSION(KLON)    :: ZUT             ! updraft temperature (K)
REAL, DIMENSION(KLON)    :: ZUW1, ZUW2      ! square of updraft vert.
                                            ! velocity at levels k and k+1
REAL, DIMENSION(KLON)    :: ZE1,ZE2,ZD1,ZD2 ! fractional entrainm./detrain
                                            ! rates at levels k and k+1
REAL, DIMENSION(KLON)    :: ZMIXF           ! critical mixed fraction  
REAL, DIMENSION(KLON)    :: ZCPH            ! specific heat C_ph 
REAL, DIMENSION(KLON)    :: ZLV, ZLS        ! latent heat of vaporis., sublim.       
REAL, DIMENSION(KLON)    :: ZURV            ! updraft water vapor at level k+1
REAL, DIMENSION(KLON)    :: ZPI             ! Pi=(P0/P)**(Rd/Cpd)  
REAL, DIMENSION(KLON)    :: ZTHEUL          ! theta_e for undilute ascent
REAL, DIMENSION(KLON)    :: ZWORK1, ZWORK2, ZWORK3, ZWORK4, ZWORK5,   &
                            ZWORK6          ! work arrays
INTEGER, DIMENSION(KLON) :: IWORK           ! wok array
LOGICAL, DIMENSION(KLON) :: GWORK1, GWORK2, GWORK4
                                            ! work arrays
LOGICAL, DIMENSION(KLON,KLEV) :: GWORK6     ! work array
!
!
!-------------------------------------------------------------------------------
!
!        0.3   Set loop bounds
!              ---------------
!
IKB = 1 + JCVEXB 
IKE = KLEV - JCVEXT 
IIE = KLON
!
!
!*       1.     Initialize updraft properties and local variables
!               -------------------------------------------------
!
ZEPSA      = XRV / XRD 
ZRDOCP     = XRD / XCPD
!
PUMF(:,:)  = 0.
PUER(:,:)  = 0.
PUDR(:,:)  = 0.
PUTHL(:,:) = 0.
PUTHV(:,:) = 0.
PUTT(:,:)  = 0.
PURW(:,:)  = 0.
PURC(:,:)  = 0.
PURI(:,:)  = 0.
PUPR(:,:)  = 0.
PURR(:,:)  = 0.
PURS(:,:)  = 0.
PUTPR(:)   = 0.
ZUW1(:)    = PWLCL(:) * PWLCL(:)
ZUW2(:)    = 0.
ZE1(:)     = 1.
ZD1(:)     = 0.
PCAPE(:)   = 0.
KCTL(:)    = IKB
KETL(:)    = KLCL(:)
GWORK2(:)  = .TRUE.
ZPI(:)     = 1.
ZWORK3(:)  = 0.
ZWORK4(:)  = 0.
ZWORK5(:)  = 0.
ZWORK6(:)  = 0.
GWORK1(:)  = .FALSE.
GWORK4(:)  = .FALSE.
!
!
!*       1.1    Compute undilute updraft theta_e for CAPE computations
!               Bolton (1980) formula.
!               Define accurate enthalpy for updraft
!               -----------------------------------------------------
!
ZTHEUL(:) = PTLCL(:) * ( PTHLCL(:) / PTLCL(:) ) ** ( 1. - 0.28 * PRVLCL(:) )  &
            * EXP( ( 3374.6525 / PTLCL(:) - 2.5403 ) *                        &
                                   PRVLCL(:) * ( 1. + 0.81 * PRVLCL(:) ) )
!
!
ZWORK1(:) = ( XCPD + PRVLCL(:) * XCPV ) * PTLCL(:)                            &
            + ( 1. + PRVLCL(:) ) * XG * PZLCL(:)
!
!
!*       2.     Set updraft properties between DPL and LCL
!               ------------------------------------------
!
JKP = MAXVAL( KLCL(:) )
JKM = MINVAL( KDPL(:) )
DO JK = JKM, JKP
   DO JI = 1, IIE
   IF ( JK >= KDPL(JI) .AND. JK < KLCL(JI) ) THEN
        PUMF(JI,JK)  = PMFLCL(JI)
        PUTHL(JI,JK) = ZWORK1(JI) 
        PUTHV(JI,JK) = PTHLCL(JI) * ( 1. + ZEPSA * PRVLCL(JI) ) /             &
                                  ( 1. + PRVLCL(JI) )
        PURW(JI,JK)  = PRVLCL(JI) 
   END IF
   END DO
END DO                        
!
!
!*       3.     Enter loop for updraft computations
!               ------------------------------------
!
! Correction for reproduciblity
!JKMIN = MINVAL( KLCL(:) ) - 1
JKMIN = MINVAL( KLCL(:) ) - 2
DO JK = MAX( IKB + 1, JKMIN ), IKE - 1
  ZWORK6(:) = 1.
  JKP = JK + 1  
!
  GWORK4(:) = JK >= KLCL(:) - 1 
  GWORK1(:) = GWORK4(:) .AND. GWORK2(:) ! this mask is used to confine
                           ! updraft computations between the LCL and the CTL
!                                                         
  WHERE( JK == KLCL(:) - 1 ) ZWORK6(:) = 0. ! factor that is used in buoyancy
                                        ! computation at first level above LCL
!
!
!*       4.     Estimate condensate, L_v L_i, Cph and theta_v at level k+1   
!               ----------------------------------------------------------
!
    ZWORK1(:) = PURC(:,JK) + PURR(:,JK)
    ZWORK2(:) = PURI(:,JK) + PURS(:,JK)
    CALL CONVECT_CONDENS( KLON, KICE, PPRES(:,JKP), PUTHL(:,JK), PURW(:,JK),&
                          ZWORK1, ZWORK2, PZ(:,JKP), GWORK1, ZUT, ZURV,     &
                          PURC(:,JKP), PURI(:,JKP), ZLV, ZLS, ZCPH )
!
!
  ZPI(:) = ( XP00 / PPRES(:,JKP) ) ** ZRDOCP   
  WHERE ( GWORK1(:) )
!
    PUTHV(:,JKP) = ZPI(:) * ZUT(:) * ( 1. + ZEPSA * ZURV(:) )           &  
                         / ( 1. + PURW(:,JK) )     
    PUTT(:,JKP) = ZUT(:)
!
!
!*       5.     Compute square of vertical velocity using entrainment   
!               at level k
!               -----------------------------------------------------
!    
    ZWORK3(:) = PZ(:,JKP) - PZ(:,JK) * ZWORK6(:) -         &
                     ( 1. - ZWORK6(:) ) * PZLCL(:)          ! level thickness  
    ZWORK4(:) = PTHV(:,JK) * ZWORK6(:) +                   &
                 ( 1. - ZWORK6(:) ) * PTHVELCL(:)
    ZWORK5(:) = 2. * ZUW1(:) * PUER(:,JK) / MAX( .1, PUMF(:,JK) )
    ZUW2(:)   = ZUW1(:) + ZWORK3(:) * XNHGAM * XG *        & 
                  ( ( PUTHV(:,JK) + PUTHV(:,JKP) ) /       &
                  ( ZWORK4(:) + PTHV(:,JKP) ) - 1. )       & ! buoyancy term
                - ZWORK5(:)                                  ! entrainment term
!
!
!*       6.     Update total precipitation: dr_r=(r_c+r_i)*exp(-rate*dz)  
!               --------------------------------------------------------
!
!                    compute level mean vertical velocity  
    ZWORK2(:)   = 0.5 *                                                    &
                       ( SQRT( MAX( 1.E-2, ZUW2(:) ) ) +                   &
                         SQRT( MAX( 1.E-2, ZUW1(:) ) ) )          
    PURR(:,JKP) = 0.5 * ( PURC(:,JK) + PURC(:,JKP) + PURI(:,JK) + PURI(:,JKP) )&
                      * ( 1. - EXP( - XRCONV  * ZWORK3(:) / ZWORK2(:) ) )
    PUPR(:,JKP) = PURR(:,JKP) * PUMF(:,JK) ! precipitation rate ( kg water / s)
    PUTPR(:)    = PUTPR(:) + PUPR(:,JKP)   ! total precipitation rate
    ZWORK2(:)   = PURR(:,JKP) / MAX( 1.E-8, PURC(:,JKP) + PURI(:,JKP) )
    PURR(:,JKP) = ZWORK2(:) * PURC(:,JKP)          ! liquid precipitation
    PURS(:,JKP) = ZWORK2(:) * PURI(:,JKP)          ! solid precipitation
!
!
!*       7.     Update r_c, r_i, enthalpy, r_w  for precipitation 
!               -------------------------------------------------------
!
    PURW(:,JKP)  = PURW(:,JK) - PURR(:,JKP) - PURS(:,JKP) 
    PURC(:,JKP)  = PURC(:,JKP) - PURR(:,JKP)
    PURI(:,JKP)  = PURI(:,JKP) - PURS(:,JKP)       
    PUTHL(:,JKP) = ( XCPD + PURW(:,JKP) * XCPV ) * ZUT(:)                     &
                   + ( 1. + PURW(:,JKP) ) * XG * PZ(:,JKP)                    &
                   - ZLV(:) * PURC(:,JKP) - ZLS(:) * PURI(:,JKP)             
!    
    ZUW1(:)      = ZUW2(:)       
!
  END WHERE
!
!
!*       8.     Compute entrainment and detrainment using conservative
!               variables adjusted for precipitation ( not for entrainment)
!               -----------------------------------------------------------
!
!*       8.1    Compute critical mixed fraction by estimating unknown  
!               T^mix r_c^mix and r_i^mix from enthalpy^mix and r_w^mix
!               We determine the zero crossing of the linear curve
!               evaluating the derivative using ZMIXF=0.1.
!               -----------------------------------------------------
!    
    ZMIXF(:)  = 0.1   ! starting value for critical mixed fraction
    ZWORK1(:) = ZMIXF(:) * PTHL(:,JKP)                                     &
                     + ( 1. - ZMIXF(:) ) * PUTHL(:,JKP) ! mixed enthalpy
    ZWORK2(:) = ZMIXF(:) * PRW(:,JKP)                                      &
                     + ( 1. - ZMIXF(:) ) * PURW(:,JKP)  ! mixed r_w
!
    CALL CONVECT_CONDENS( KLON, KICE, PPRES(:,JKP), ZWORK1, ZWORK2,        &
                          PURC(:,JKP), PURI(:,JKP), PZ(:,JKP), GWORK1, ZUT,&
                          ZWORK3, ZWORK4, ZWORK5, ZLV, ZLS, ZCPH )
!        put in enthalpy and r_w and get T r_c, r_i (ZUT, ZWORK4-5)
!        
     ! compute theta_v of mixture
    ZWORK3(:) = ZUT(:) * ZPI(:) * ( 1. + ZEPSA * (                         &
                ZWORK2(:) - ZWORK4(:) - ZWORK5(:) ) ) / ( 1. + ZWORK2(:) )
     ! compute final value of critical mixed fraction using theta_v
     ! of mixture, grid-scale and updraft
    ZMIXF(:) = MAX( 0., PUTHV(:,JKP) - PTHV(:,JKP) ) * ZMIXF(:) /          &
                              ( PUTHV(:,JKP) - ZWORK3(:) + 1.E-10 )
    ZMIXF(:) = MAX( 0., MIN( 1., ZMIXF(:) ) )
!    
!
!*       8.2     Compute final midlevel values for entr. and detrainment    
!                after call of distribution function
!                -------------------------------------------------------
!    
!
    CALL CONVECT_MIXING_FUNCT ( KLON, ZMIXF, 1, ZE2, ZD2 )
!       Note: routine MIXING_FUNCT returns fractional entrainm/detrainm. rates
!
! ZWORK1(:) = XENTR * PMFLCL(:) * PDPRES(:,JKP) / XCRAD ! rate of env. inflow
!*MOD
  zwork1(:) = xentr * xg / xcrad * pumf(:,jk) * ( pz(:,jkp) - pz(:,jk) )
! ZWORK1(:) = XENTR * pumf(:,jk) * PDPRES(:,JKP) / XCRAD ! rate of env. inflow
!*MOD
  ZWORK2(:) = 0.
  WHERE ( GWORK1(:) ) ZWORK2(:) = 1.
  ZE2(:) = .5; ZD2(:) = .6 ! set entrainment=detrainment for better
                           ! mass flux profiles in deep continental convection
  WHERE ( PUTHV(:,JKP) > PTHV(:,JKP) )
    PUER(:,JKP) = 0.5 * ZWORK1(:) * ( ZE1(:) + ZE2(:) ) * ZWORK2(:)
    PUDR(:,JKP) = 0.5 * ZWORK1(:) * ( ZD1(:) + ZD2(:) ) * ZWORK2(:)
  ELSEWHERE
    PUER(:,JKP) = 0.
    PUDR(:,JKP) = ZWORK1(:) * ZWORK2(:)
  END WHERE
!
!*       8.3     Determine equilibrium temperature level
!                --------------------------------------
!
   WHERE ( PUTHV(:,JKP) > PTHV(:,JKP) .AND. JK > KLCL(:) + 1 &   
           .AND. GWORK1(:) )
         KETL(:) = JKP            ! equilibrium temperature level 
   END WHERE
!
!*       8.4     If the calculated detrained mass flux is greater than    
!                the total updraft mass flux, or vertical velocity is
!                negative, all cloud mass detrains at previous model level,
!                exit updraft calculations - CTL is attained
!                -------------------------------------------------------
!
  WHERE( GWORK1(:) )                                                   &
        GWORK2(:) = PUMF(:,JK) - PUDR(:,JKP) > 10. .AND. ZUW2(:) > 0.        
  WHERE ( GWORK2(:) ) KCTL(:) = JKP   ! cloud top level
!!!! Correction Bug C.Lac 30/10/08 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  KCTL(:) = MIN( KCTL(:), IKE-1 )
  GWORK1(:) = GWORK2(:) .AND. GWORK4(:)
!
  IF ( COUNT( GWORK2(:) ) == 0 ) EXIT           
!
!
!*       9.   Compute CAPE for undilute ascent using theta_e and 
!             theta_es instead of theta_v. This estimation produces 
!             a significantly larger value for CAPE than the actual one.
!             ----------------------------------------------------------
!
  WHERE ( GWORK1(:) )
!
    ZWORK3(:)   = PZ(:,JKP) - PZ(:,JK) * ZWORK6(:) -                      &
                  ( 1. - ZWORK6(:) ) *  PZLCL(:)              ! level thickness
    ZWORK2(:)   = PTHES(:,JK) + ( 1. - ZWORK6(:) ) *                      &
     ( PTHES(:,JKP) - PTHES(:,JK) ) / ( PZ(:,JKP) - PZ(:,JK) ) *          &
     ( PZLCL(:) - PZ(:,JK) ) ! linear interpolation for theta_es at LCL
                            ! ( this is only done for model level just above LCL
!
    ZWORK1(:) = ( 2. * ZTHEUL(:) ) / ( ZWORK2(:) + PTHES(:,JKP) ) - 1.   
    PCAPE(:)  = PCAPE(:) + XG * ZWORK3(:) * MAX( 0., ZWORK1(:) )
!
!
!*       10.   Compute final values of updraft mass flux, enthalpy, r_w 
!              at level k+1    
!              --------------------------------------------------------
!    
    PUMF(:,JKP)  = PUMF(:,JK) - PUDR(:,JKP) + PUER(:,JKP) 
    PUMF(:,JKP)  = MAX( PUMF(:,JKP), 0.1 )
    PUTHL(:,JKP) = ( PUMF(:,JK) * PUTHL(:,JK) +                              &
                     PUER(:,JKP) * PTHL(:,JK) - PUDR(:,JKP) * PUTHL(:,JK) )  &
                    / PUMF(:,JKP) + PUTHL(:,JKP) - PUTHL(:,JK)
    PURW(:,JKP)  = ( PUMF(:,JK) * PURW(:,JK) +                               &
                     PUER(:,JKP) * PRW(:,JK) - PUDR(:,JKP) * PURW(:,JK) )    &
                    / PUMF(:,JKP) - PURR(:,JKP) - PURS(:,JKP)
!    
!
    ZE1(:) = ZE2(:) ! update fractional entrainment/detrainment
    ZD1(:) = ZD2(:)
!
  END WHERE
!
END DO
!
!*       12.1    Set OTRIG to False if cloud thickness < XCDEPTH
!                or CAPE < 1
!                ------------------------------------------------
!
    DO JI = 1, IIE
          JK  = KCTL(JI)
          OTRIG(JI) = PZ(JI,JK) - PZLCL(JI) >= XCDEPTH               &
                     .AND. PCAPE(JI) > 1. 
    END DO
    WHERE( .NOT. OTRIG(:) )
          KCTL(:) = IKB 
    END WHERE
KETL(:) = MAX( KETL(:), KLCL(:) + 2 )
KETL(:) = MIN( KETL(:), KCTL(:) )
!
!
!*       12.2    If the ETL and CTL are the same detrain updraft mass   
!                flux at this level
!                ------------------------------------------------------- 
!
ZWORK1(:) = 0.
WHERE ( KETL(:) == KCTL(:) ) ZWORK1(:) = 1.
!
DO JI = 1, IIE
    JK = KETL(JI) 
    PUDR(JI,JK)   = PUDR(JI,JK) +                                    &
                          ( PUMF(JI,JK) - PUER(JI,JK) )  * ZWORK1(JI)  
    PUER(JI,JK)   = PUER(JI,JK) * ( 1. - ZWORK1(JI) )
    PUMF(JI,JK)   = PUMF(JI,JK) * ( 1. - ZWORK1(JI) )
    JKP = KCTL(JI) + 1
    PUER(JI,JKP)  = 0. ! entrainm/detr rates have been already computed
    PUDR(JI,JKP)  = 0. ! at level KCTL+1, set them to zero
    PURW(JI,JKP)  = 0.
    PURC(JI,JKP)  = 0.
    PURI(JI,JKP)  = 0.
    PUTHL(JI,JKP) = 0.
    PURI(JI,JKP+1)= 0.
    PURC(JI,JKP+1)= 0.
END DO
!    
!*       12.3    Adjust mass flux profiles, detrainment rates, and   
!                precipitation fallout rates to reflect linear decrease
!                in mass flux between the ETL and CTL
!                -------------------------------------------------------        
! 
ZWORK1(:) = 0.
JK1 = MINVAL( KETL(:) )
JK2 = MAXVAL( KCTL(:) )
DO JK = JK1, JK2
    DO JI = 1, IIE
    IF( JK > KETL(JI) .AND. JK <= KCTL(JI) ) THEN
        ZWORK1(JI) = ZWORK1(JI) + PDPRES(JI,JK)
    END IF
    END DO
END DO
!
DO JI = 1, IIE
    JK = KETL(JI) 
    ZWORK1(JI) = PUMF(JI,JK) / MAX( 1., ZWORK1(JI) )
END DO
!
DO JK = JK1 + 1, JK2
    JKP = JK - 1
    DO JI = 1, IIE
    IF ( JK > KETL(JI) .AND. JK <= KCTL(JI) ) THEN
      ! PUTPR(JI)    = PUTPR(JI) - ( PURR(JI,JK) + PURS(JI,JK) ) * PUMF(JI,JKP)      
        PUTPR(JI)    = PUTPR(JI) - PUPR(JI,JK)
        PUDR(JI,JK)  = PDPRES(JI,JK) * ZWORK1(JI)
        PUMF(JI,JK)  = PUMF(JI,JKP) - PUDR(JI,JK)
        PUPR(JI,JK)  = PUMF(JI,JKP) * ( PURR(JI,JK) + PURS(JI,JK) )
        PUTPR(JI)    = PUTPR(JI) + PUPR(JI,JK)
    END IF
    END DO
END DO
!
!         12.4   Set mass flux and entrainment in the source layer.
!                Linear increase throughout the source layer.
!                -------------------------------------------------------
!
!IWORK(:) = MIN( KPBL(:), KLCL(:) - 1 )
IWORK(:) = KPBL(:)
DO JI = 1, IIE
     JK  = KDPL(JI)
     JKP = IWORK(JI)
!          mixed layer depth
     ZWORK2(JI) = PPRES(JI,JK) - PPRES(JI,JKP) + PDPRES(JI,JK)
END DO
!
JKP = MAXVAL( IWORK(:) )
DO JK = JKM, JKP
   DO JI = 1, IIE
   IF ( JK >= KDPL(JI)  .AND. JK <= IWORK(JI) ) THEN
       PUER(JI,JK) = PUER(JI,JK) + PMFLCL(JI) * PDPRES(JI,JK) / ( ZWORK2(JI) + 0.1 )
       PUMF(JI,JK) = PUMF(JI,JK-1) + PUER(JI,JK)
   END IF
   END DO
END DO
!
!
!*       13.   If cloud thickness is smaller than  3 km, no
!              convection is allowed
!              Nota: For technical reasons, we stop the convection
!                    computations in this case and do not go back to
!                    TRIGGER_FUNCT to look for the next unstable LCL
!                    which could produce a thicker cloud.
!              ---------------------------------------------------
!
GWORK6(:,:) = SPREAD( OTRIG(:), DIM=2, NCOPIES=KLEV )
WHERE ( .NOT. OTRIG(:) ) PUTPR(:) = 0.
WHERE ( .NOT. GWORK6(:,:) )
    PUMF(:,:)  = 0.
    PUDR(:,:)  = 0.
    PUER(:,:)  = 0.
    PUTHL(:,:) = PTHL(:,:)
    PURW(:,:)  = PRW(:,:)
    PUPR(:,:)  = 0.
    PURC(:,:)  = 0.
    PURI(:,:)  = 0.
    PURR(:,:)  = 0.
    PURS(:,:)  = 0.
END WHERE
!
END SUBROUTINE CONVECT_UPDRAFT
