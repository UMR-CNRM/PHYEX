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
      MODULE MODI_CONVECT_DOWNDRAFT
!     #################
!
INTERFACE
!
       SUBROUTINE CONVECT_DOWNDRAFT( KLON, KLEV,                                &
                                   KICE, PPRES, PDPRES, PZ, PTH, PTHES,       &
                                   PRW, PRC, PRI,                             &
                                   PPREF, KLCL, KCTL, KETL,                   &
                                   PUTHL, PURW, PURC, PURI,                   &
                                   PDMF, PDER, PDDR, PDTHL, PDRW,             &
                                   PMIXF, PDTEVR, KLFS, KDBL, KML,            &
                                   PDTEVRF )
!
INTEGER,                    INTENT(IN) :: KLON  ! horizontal dimension
INTEGER,                    INTENT(IN) :: KLEV  ! vertical dimension
INTEGER,                    INTENT(IN) :: KICE  ! flag for ice ( 1 = yes,
                                                !                0 = no ice )
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PTH   ! grid scale theta        
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PTHES ! grid scale saturated theta_e 
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PRW   ! grid scale total water  
                                                ! mixing ratio 
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PRC   ! grid scale r_c (cloud water) 
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PRI   ! grid scale r_i (cloud ice) 
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PPRES ! pressure (Pa)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PDPRES! pressure difference between 
						! bottom and top of layer (Pa) 
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PZ    ! level height (m)
INTEGER, DIMENSION(KLON),   INTENT(IN) :: KLCL  ! contains vert. index of LCL
INTEGER, DIMENSION(KLON),   INTENT(IN) :: KCTL  ! contains vert. index of CTL 
INTEGER, DIMENSION(KLON),   INTENT(IN) :: KETL  ! contains vert. index of 
						! equilibrium (zero buoyancy) level 
INTEGER, DIMENSION(KLON),   INTENT(IN) :: KML   ! " vert. index of melting level
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PUTHL ! updraft enthalpy (J/kg)      
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PURW  ! updraft total water (kg/kg)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PURC  ! updraft r_c (kg/kg)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PURI  ! updraft r_i (kg/kg)
REAL, DIMENSION(KLON),      INTENT(IN) :: PPREF ! precipitation efficiency
!
!
REAL, DIMENSION(KLON,KLEV), INTENT(OUT):: PDMF   ! downdraft mass flux (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(OUT):: PDER   ! downdraft entrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(OUT):: PDDR   ! downdraft detrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(OUT):: PDTHL  ! downdraft enthalpy (J/kg)
REAL, DIMENSION(KLON,KLEV), INTENT(OUT):: PDRW   ! downdraft total water (kg/kg)
REAL, DIMENSION(KLON),      INTENT(OUT):: PMIXF  ! mixed fraction at LFS
REAL, DIMENSION(KLON),      INTENT(OUT):: PDTEVR ! total downdraft evaporation
                                                 ! rate at LFS (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(OUT):: PDTEVRF! downdraft evaporation rate
INTEGER, DIMENSION(KLON),  INTENT(OUT):: KLFS    ! contains vert. index of LFS 
INTEGER, DIMENSION(KLON),  INTENT(OUT):: KDBL    ! contains vert. index of DBL   
!
END SUBROUTINE CONVECT_DOWNDRAFT
!
END INTERFACE
!
END MODULE MODI_CONVECT_DOWNDRAFT
!    ##########################################################################
     SUBROUTINE CONVECT_DOWNDRAFT( KLON, KLEV,                                &
                                   KICE, PPRES, PDPRES, PZ, PTH, PTHES,       &
                                   PRW, PRC, PRI,                             &
                                   PPREF, KLCL, KCTL, KETL,                   &
                                   PUTHL, PURW, PURC, PURI,                   &
                                   PDMF, PDER, PDDR, PDTHL, PDRW,             &
                                   PMIXF, PDTEVR, KLFS, KDBL, KML,            &
                                   PDTEVRF )
!    ##########################################################################
!
!!**** Compute downdraft properties from LFS to DBL. 
!!
!!
!!    PDRPOSE                                                       
!!    -------
!!      The purpose of this routine is to determine downdraft properties
!!      ( mass flux, thermodynamics ) 
!!
!!
!!**  METHOD
!!    ------
!!      Computations are done at every model level starting from top.
!!      The use of masks allows to optimise the inner loops (horizontal loops).
!!      
!!     
!!
!!    EXTERNAL
!!    --------
!!     Routine CONVECT_SATMIXRATIO        
!!                
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      Module MODD_CST
!!          XG                 ! gravity constant
!!          XPI                ! Pi
!!          XP00               ! reference pressure
!!          XRD, XRV           ! gaz  constants for dry air and water vapor
!!          XCPD               ! Cpd (dry air)
!!          XCPV, XCL, XCI     ! Cp of water vapor, liquid water and ice
!!          XTT                ! triple point temperature
!!          XLVTT, XLSTT       ! vaporisation/sublimation heat at XTT
!!
!!      Module MODD_CONVPAR
!!          XCRAD              ! cloud radius
!!          XZPBL              ! thickness of downdraft detrainment layer
!!          XENTR              ! entrainment constant in pressure coordinates
!!          XRHDBC             ! relative humidity in downdraft below cloud
!!
!!     Module MODD_CONVPAREXT
!!          JCVEXB, JCVEXT     ! extra levels on the vertical boundaries
!!
!!    REFERENCE
!!    ---------
!!
!!      Book1,2 of documentation ( routine CONVECT_DOWNDRAFT)
!!      Kain and Fritsch, 1993, Meteor. Monographs, Vol.
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    07/11/95 
!!   Last modified  04/10/97
!!   C.Lac          27/09/10 modification loop index for reproducibility
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
!
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
INTEGER,                    INTENT(IN) :: KLON  ! horizontal dimension
INTEGER,                    INTENT(IN) :: KLEV  ! vertical dimension
INTEGER,                    INTENT(IN) :: KICE  ! flag for ice ( 1 = yes,
                                                !                0 = no ice )
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PTH   ! grid scale theta        
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PTHES ! grid scale saturated theta_e 
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PRW   ! grid scale total water  
                                                ! mixing ratio 
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PRC   ! grid scale r_c (cloud water) 
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PRI   ! grid scale r_i (cloud ice) 
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PPRES ! pressure (Pa)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PDPRES! pressure difference between 
						! bottom and top of layer (Pa) 
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PZ    ! level height (m)
INTEGER, DIMENSION(KLON),   INTENT(IN) :: KLCL  ! contains vert. index of LCL
INTEGER, DIMENSION(KLON),   INTENT(IN) :: KCTL  ! contains vert. index of CTL 
INTEGER, DIMENSION(KLON),   INTENT(IN) :: KETL  ! contains vert. index of 
						! equilibrium (zero buoyancy) level 
INTEGER, DIMENSION(KLON),   INTENT(IN) :: KML   ! " vert. index of melting level
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PUTHL ! updraft enthalpy (J/kg)      
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PURW  ! updraft total water (kg/kg)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PURC  ! updraft r_c (kg/kg)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PURI  ! updraft r_i (kg/kg)
REAL, DIMENSION(KLON),      INTENT(IN) :: PPREF ! precipitation efficiency
!
!
REAL, DIMENSION(KLON,KLEV), INTENT(OUT):: PDMF   ! downdraft mass flux (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(OUT):: PDER   ! downdraft entrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(OUT):: PDDR   ! downdraft detrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(OUT):: PDTHL  ! downdraft enthalpy (J/kg)
REAL, DIMENSION(KLON,KLEV), INTENT(OUT):: PDRW   ! downdraft total water (kg/kg)
REAL, DIMENSION(KLON),      INTENT(OUT):: PMIXF  ! mixed fraction at LFS
REAL, DIMENSION(KLON),      INTENT(OUT):: PDTEVR ! total downdraft evaporation
                                                 ! rate at LFS (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(OUT):: PDTEVRF! downdraft evaporation rate
INTEGER, DIMENSION(KLON),  INTENT(OUT):: KLFS    ! contains vert. index of LFS 
INTEGER, DIMENSION(KLON),  INTENT(OUT):: KDBL    ! contains vert. index of DBL   
!
!*       0.2   Declarations of local variables :
!
INTEGER :: IIE, IKB, IKE     ! horizontal + vertical loop bounds
INTEGER :: JK, JKP, JKM, JKT ! vertical loop index
INTEGER :: JI, JL            ! horizontal loop index
INTEGER :: JITER          ! iteration loop index
REAL    :: ZRDOCP         ! R_d / C_pd
REAL    :: ZEPS           ! R_d / R_v
!
INTEGER, DIMENSION(KLON) :: IDDT      ! top level of detrainm. layer
REAL, DIMENSION(KLON)    :: ZTHE      ! environm. theta_e (K)
REAL, DIMENSION(KLON)    :: ZDT, ZDTP ! downdraft temperature (K)
REAL, DIMENSION(KLON)    :: ZCPH      ! specific heat C_ph 
REAL, DIMENSION(KLON)    :: ZLV, ZLS  ! latent heat of vaporis., sublim.       
REAL, DIMENSION(KLON)    :: ZDDT      ! thickness (hPa) of detrainm. layer
REAL, DIMENSION(KLON)    :: ZPI       ! Pi=(P0/P)**(Rd/Cpd)  
REAL, DIMENSION(KLON)    :: ZWORK1, ZWORK2, ZWORK3, ZWORK4 ! work arrays
LOGICAL, DIMENSION(KLON) :: GWORK1                         ! work array
!
!
!-------------------------------------------------------------------------------
!
!        0.3    Set loop bounds
!               ---------------
!
IIE = KLON
IKB = 1 + JCVEXB 
IKE = KLEV - JCVEXT 
!
!
!*       1.     Initialize downdraft properties
!               -------------------------------
!
ZRDOCP     = XRD / XCPD
ZEPS       = XRD / XRV
PDMF(:,:)  = 0.
PDER(:,:)  = 0.
PDDR(:,:)  = 0.
PDRW(:,:)  = 0.
PDTHL(:,:) = 0.
PDTEVR(:)  = 0.
PMIXF(:)   = 0.
ZTHE(:)    = 0.
ZDDT(:)    = PDPRES(:,IKB+2)
KDBL(:)    = IKB + 1
KLFS(:)    = IKB + 1
IDDT(:)    = KDBL(:) + 1
!  
!
!*       2.     Determine the LFS by looking for minimum of environmental 
!               saturated theta_e 
!               ----------------------------------------------------------
!
ZWORK1(:) = 900.   ! starting value for search of minimum envir. theta_e
DO JK = MINVAL( KLCL(:) ) + 2, MAXVAL( KETL(:) )
   DO JI = 1, IIE
      GWORK1(JI) = JK >= KLCL(JI) + 2 .AND. JK < KETL(JI)  
      IF ( GWORK1(JI) .AND. ZWORK1(JI) > PTHES(JI,JK) ) THEN
         KLFS(JI)   = JK
         ZWORK1(JI) = MIN( ZWORK1(JI), PTHES(JI,JK) )
      END IF
   END DO
END DO
!
!
!*       3.     Determine the mixed fraction using environmental and updraft
!               values of theta_e at LFS
!               ---------------------------------------------------------   
!
DO JI = 1, IIE
    JK = KLFS(JI)
    ZPI(JI)    = ( XP00 / PPRES(JI,JK) ) ** ZRDOCP
      ! compute updraft theta_e
    ZWORK3(JI) = PURW(JI,JK) - PURC(JI,JK) - PURI(JI,JK)
    ZDT(JI)    = PTH(JI,JK) / ZPI(JI) 
    ZLV(JI)    = XLVTT + ( XCPV - XCL ) * ( ZDT(JI) - XTT )                   
    ZLS(JI)    = XLSTT + ( XCPV - XCI ) * ( ZDT(JI) - XTT )                   
    ZCPH(JI)   = XCPD + XCPV * PURW(JI,JK)
    ZDT(JI)    = ( PUTHL(JI,JK) - ( 1. + PURW(JI,JK) ) * XG * PZ(JI,JK)       &
                 + ZLV(JI) * PURC(JI,JK) + ZLS(JI) * PURI(JI,JK) ) / ZCPH(JI)           
    ZWORK1(JI) = ZDT(JI) * ZPI(JI) ** ( 1. - 0.28 * ZWORK3(JI) )              &
                                  * EXP( ( 3374.6525 / ZDT(JI) - 2.5403 )     &
                                  * ZWORK3(JI) * ( 1. + 0.81 * ZWORK3(JI) ) )
      ! compute environmental theta_e
    ZDT(JI)    = PTH(JI,JK) / ZPI(JI)
    ZLV(JI)    = XLVTT + ( XCPV - XCL ) * ( ZDT(JI) - XTT )                   
    ZLS(JI)    = XLSTT + ( XCPV - XCI ) * ( ZDT(JI) - XTT )                   
    ZWORK3(JI) = PRW(JI,JK) - PRC(JI,JK) - PRI(JI,JK)
    ZCPH(JI)   = XCPD + XCPV * PRW(JI,JK)
    ZWORK2(JI) = ZDT(JI) * ZPI(JI) ** ( 1. - 0.28 * ZWORK3(JI) )              &
                                  * EXP( ( 3374.6525 / ZDT(JI) - 2.5403 )     &
                                  * ZWORK3(JI) * ( 1. + 0.81 * ZWORK3(JI) ) )
      ! compute mixed fraction
    PMIXF(JI)  = MAX( 0., ( ZWORK1(JI) - PTHES(JI,JK) ) )                   &
                  / ( ZWORK1(JI) - ZWORK2(JI) + 1.E-10 )
    PMIXF(JI)  = MAX(0., MIN( 1., PMIXF(JI) ) )
    ZWORK4(JI) = PPRES(JI,JK)
END DO
!
!
!*       4.     Estimate the effect of melting on the downdraft  
!               ---------------------------------------------
!
ZWORK1(:) = 0.
      ! use total solid precipitation
!DO JK = IKB + 1, IKE
!    ZWORK1(:) = ZWORK1(:) + PURS(:,JK) ! total snow/hail content
!END DO
!
DO JI = 1, IIE
     JK  = KLCL(JI)
     JKP = KCTL(JI)
     ZWORK1(JI) = 0.5 * ( PURW(JI,JK) - PURW(JI,JKP) )
END DO
!
      ! temperature perturbation due to melting at LFS
ZWORK3(:) = 0.
WHERE( KML(:) > IKB + 2 )
	  ZWORK3(:) = ZWORK1(:) * ( ZLS(:) - ZLV(:) ) / ZCPH(:)
	  ZDT(:)    = ZDT(:) - ZWORK3(:) * REAL(KICE)
END WHERE
!
!
!*       5.     Initialize humidity at LFS as a saturated mixture of
!               updraft and environmental air
!               -----------------------------------------------------    
!
DO JI = 1, IIE
     JK = KLFS(JI)
     PDRW(JI,JK)  = PMIXF(JI) * PRW(JI,JK) + ( 1. - PMIXF(JI) ) * PURW(JI,JK)
     ZWORK2(JI)   = PDRW(JI,JK) - ( 1. - PMIXF(JI) )                          &
                                     * ( PURC(JI,JK) + PURI(JI,JK) )
END DO
!
!
!*       6.1    Determine the DBL by looking for level where the envir.
!               theta_es at the LFS corrected by melting effects  becomes
!               larger than envir. value
!               ---------------------------------------------------------
!
      ! compute satur. mixing ratio for melting corrected temperature
CALL CONVECT_SATMIXRATIO( KLON, ZWORK4, ZDT, ZWORK3, ZLV, ZLS, ZCPH )  
!
      ! compute envir. saturated theta_e for melting corrected temperature
    ZWORK1(:) = MIN( ZWORK2(:), ZWORK3(:) )
    ZWORK3(:) = ZWORK3(:) * ZWORK4(:) / ( ZWORK3(:) + ZEPS ) ! sat. pressure
    ZWORK3(:) = ALOG( ZWORK3(:) / 613.3 )
              ! dewp point temperature
    ZWORK3(:) = ( 4780.8 - 32.19 * ZWORK3(:) ) / ( 17.502 - ZWORK3(:) )
              ! adiabatic saturation temperature
    ZWORK3(:) = ZWORK3(:) - ( .212 + 1.571E-3 * ( ZWORK3(:) - XTT )          &
                  - 4.36E-4 * ( ZDT(:) - XTT ) ) * ( ZDT(:) - ZWORK3(:) )
    ZWORK4(:) = SIGN(0.5, ZWORK2(:) - ZWORK3(:) )
    ZDT(:)    = ZDT(:) * ( .5 + ZWORK4(:) ) + ( .5 - ZWORK4(:) ) * ZWORK3(:) 
    ZWORK2(:) = ZDT(:) * ZPI(:) ** ( 1. - 0.28 * ZWORK2(:) )                 &
                                  * EXP( ( 3374.6525 / ZDT(:) - 2.5403 )     &
                                  * ZWORK1(:) * ( 1. + 0.81 * ZWORK1(:) ) )
!
GWORK1(:) = .TRUE.
JKM = MAXVAL( KLFS(:) )
DO JK = JKM - 1, IKB + 1, -1
  DO JI = 1, IIE
     IF ( JK < KLFS(JI) .AND. ZWORK2(JI) > PTHES(JI,JK) .AND. GWORK1(JI) ) THEN
	  KDBL(JI) = JK
          GWORK1(JI) = .FALSE.
     END IF
  END DO
END DO
!
!
!*       7.     Define mass flux and entr/detr. rates at LFS
!               -------------------------------------------
!
DO JI = 1, IIE
     JK = KLFS(JI)
     ZWORK1(JI)  = PPRES(JI,JK) /                                            &
                   ( XRD * ZDT(JI) * ( 1. + ZEPS * ZWORK1(JI) ) ) ! density
     PDMF(JI,JK) = - ( 1. - PPREF(JI) ) * ZWORK1(JI) * XPI * XCRAD * XCRAD
     PDTHL(JI,JK)= ZWORK2(JI)   ! theta_l is here actually theta_e
     ZWORK2(JI)  = PDMF(JI,JK)
     PDDR(JI,JK) = 0.
     PDER(JI,JK) = - PMIXF(JI) * PDMF(JI,JK)
END DO
!
!
!         7.1   Downdraft detrainment is assumed to occur in a layer
!               of 60 hPa, determine top level IDDT of this layer
!               ---------------------------------------------------------
!
ZWORK1(:) = 0.
DO JK = IKB + 2, JKM
      ZWORK1(:) = ZWORK1(:) + PDPRES(:,JK)
      !WHERE ( JK > KDBL(:) .AND. ZWORK1(:) <= XZPBL )
      WHERE ( JK > KDBL(:) .AND. JK <= KLCL(:) )
           ZDDT(:) = ZWORK1(:) 
           IDDT(:) = JK
      END WHERE
END DO
!
!
!*       8.     Enter loop for downdraft computations. Make a first guess
!               of initial downdraft mass flux. 
!               In the downdraft computations we use theta_es instead of 
!               enthalpy as it allows to better take into account evaporation
!               effects. As the downdraft detrainment rate is zero apart 
!               from the detrainment layer, we just compute enthalpy 
!               downdraft from theta_es in this layer.
!               ----------------------------------------------------------
!
!
!
DO JK =  JKM - 1, IKB + 1, -1
  JKP = JK + 1
  DO JI = 1, IIE
    IF ( JK < KLFS(JI) .AND. JK >= IDDT(JI) )  THEN
      PDER(JI,JK)  = - ZWORK2(JI) * XENTR * PDPRES(JI,JKP) / XCRAD 
                                               ! DER and DPRES are positive
      PDMF(JI,JK)  = PDMF(JI,JKP) - PDER(JI,JK) 
      ZPI(JI)      = ( XP00 / PPRES(JI,JK) ) ** ZRDOCP
      ZDT(JI)      = PTH(JI,JK) / ZPI(JI)
      ZWORK1(JI)   = PRW(JI,JK) - PRC(JI,JK) - PRI(JI,JK)
      ZTHE(JI)     = ZDT(JI) * ZPI(JI) ** ( 1. - 0.28 * ZWORK1(JI) )           &
                               * EXP( ( 3374.6525 / ZDT(JI) - 2.5403 )         &
                               * ZWORK1(JI) * ( 1. + 0.81 * ZWORK1(JI) ) )
         ! PDTHL is here theta_es, later on in this routine this table is
         ! reskipped to enthalpy 
      PDTHL(JI,JK) = ( PDTHL(JI,JKP) * PDMF(JI,JKP) - ZTHE(JI) * PDER(JI,JK)    &
                    ) / ( PDMF(JI,JK) - 1.E-7 )      
      PDRW(JI,JK)  = ( PDRW(JI,JKP) * PDMF(JI,JKP) - PRW(JI,JK) * PDER(JI,JK)   &
                    ) / ( PDMF(JI,JK) - 1.E-7 )       
    END IF
    IF ( JK < IDDT(JI) .AND. JK >= KDBL(JI) )   THEN
      JL = IDDT(JI)
      PDDR(JI,JK)  = - PDMF(JI,JL) * PDPRES(JI,JKP) / ZDDT(JI) 
      PDMF(JI,JK)  = PDMF(JI,JKP) + PDDR(JI,JK) 
      PDTHL(JI,JK) = PDTHL(JI,JKP)
      PDRW(JI,JK)  = PDRW(JI,JKP)
    END IF
  END DO
END DO
!
!
!*       9.     Calculate total downdraft evaporation 
!               rate for given mass flux (between DBL and IDDT)
!               -----------------------------------------------
!
PDTEVRF(:,:) = 0.
! Reproducibility
!JKT = MAXVAL( IDDT(:) )
!DO JK = IKB + 1, JKT
DO JK = IKB + 1, IKE
!
       ZPI(:) = ( XP00 / PPRES(:,JK) ) ** ZRDOCP
       ZDT(:) = PTH(:,JK) / ZPI(:)
!
!*       9.1    Determine wet bulb temperature at DBL from theta_e.
!               The iteration algoritm is similar to that used in
!               routine CONVECT_CONDENS
!               --------------------------------------------------
!
   DO JITER = 1, 4
       CALL CONVECT_SATMIXRATIO( KLON, PPRES(:,JK), ZDT, ZWORK1, ZLV, ZLS, ZCPH )  
       ZDTP(:) = PDTHL(:,JK) / ( ZPI(:) ** ( 1. - 0.28 * ZWORK1(:) )         &
                      * EXP( ( 3374.6525 / ZDT(:) - 2.5403 )                 &
                             * ZWORK1(:) * ( 1. + 0.81 * ZWORK1(:) ) ) )
       ZDT(:)  = 0.4 * ZDTP(:) + 0.6 * ZDT(:) ! force convergence
   END DO
!
!
!*       9.2    Sum total downdraft evaporation rate. No evaporation
!               if actual humidity is larger than specified one.
!               -----------------------------------------------------
!
   ZWORK2(:) = ZWORK1(:) / ZDT(:) * ( XBETAW / ZDT(:) - XGAMW ) ! dr_sat/dT
   ZWORK2(:) = ZLV(:) / ZCPH(:) * ZWORK1(:) * ( 1. - XRHDBC ) /              &
                    ( 1. + ZLV(:) / ZCPH(:) * ZWORK2(:) ) ! temperature perturb                                                           ! due to evaporation
   ZDT(:)    = ZDT(:) + ZWORK2(:)
!
   CALL CONVECT_SATMIXRATIO( KLON, PPRES(:,JK), ZDT, ZWORK3, ZLV, ZLS, ZCPH )
!
   ZWORK3(:)    = ZWORK3(:) * XRHDBC
   ZWORK1(:)    = MAX( 0., ZWORK3(:) - PDRW(:,JK) ) 
   PDTEVR(:)    = PDTEVR(:) + ZWORK1(:) * PDDR(:,JK) 
   PDTEVRF(:,JK)= PDTEVRF(:,JK) + ZWORK1(:) * PDDR(:,JK) 
      ! compute enthalpie and humidity in the detrainment layer
   PDRW(:,JK)   = MAX( PDRW(:,JK), ZWORK3(:) ) 
   PDTHL(:,JK)  = ( ( XCPD + PDRW(:,JK) * XCPV ) * ZDT(:)                    &
                    + ( 1. + PDRW(:,JK) ) * XG * PZ(:,JK) ) 
!
END DO
!
!
!*      12.     If downdraft does not evaporate any water for specified 
!               relative humidity, no downdraft is allowed
!               ---------------------------------------------------------
!
ZWORK2(:) = 1.
WHERE ( PDTEVR(:) < 1. .OR. KLFS(:) == IKB + 1 ) ZWORK2(:) = 0.
DO JK = IKB, JKM
      KDBL(:)     = KDBL(:) * INT( ZWORK2(:) ) + ( 1 - INT( ZWORK2(:) ) ) * IKB
      KLFS(:)     = KLFS(:) * INT( ZWORK2(:) ) + ( 1 - INT( ZWORK2(:) ) ) * IKB
      PDMF(:,JK)  = PDMF(:,JK)  * ZWORK2(:)
      PDER(:,JK)  = PDER(:,JK)  * ZWORK2(:) 
      PDDR(:,JK)  = PDDR(:,JK)  * ZWORK2(:) 
      ZWORK1(:)   = REAL( KLFS(:) - JK )         ! use this to reset thl_d
      ZWORK1(:)   = MAX( 0.,MIN(1.,ZWORK1(:) ) ) ! and rv_d to zero above LFS
      PDTHL(:,JK) = PDTHL(:,JK) * ZWORK2(:) * ZWORK1(:)
      PDRW(:,JK)  = PDRW(:,JK)  * ZWORK2(:) * ZWORK1(:)
      PDTEVR(:)   = PDTEVR(:)   * ZWORK2(:)
      PDTEVRF(:,JK)= PDTEVRF(:,JK) * ZWORK2(:)
END DO
!
END SUBROUTINE CONVECT_DOWNDRAFT
