!     ######spl
     SUBROUTINE CONVECT_CLOSURE_SHAL( CVP_SHAL, CVPEXT, CST, D,         &
                                      PPRES, PDPRES, PZ, PLMASS,           &
                                      PTHL, PTH, PRW, PRC, PRI, OTRIG1,           &
                                      PTHC, PRWC, PRCC, PRIC, PWSUB,              &
                                      KLCL, KDPL, KPBL, KCTL,                     &
                                      PUMF, PUER, PUDR, PUTHL, PURW,              &
                                      PURC, PURI, PCAPE, PTIMEC, KFTSTEPS         )
     USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
!    ##############################################################################
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
!!    CONVECT_CLOSURE_ADJUST_SHAL
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
!!      Module MODD_CONVPAR_SHAL
!!          XA25               ! reference grid area
!!          XSTABT             ! stability factor in time integration
!!          XSTABC             ! stability factor in CAPE adjustment
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
!!   Peter Bechtold 15/11/96 change for enthalpie, r_c + r_i tendencies
!!      Tony Dore   14/10/96 Initialise local variables
!!      F Bouyssel  08/11/13 Modifications for reproductibility
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST, ONLY : CST_T
USE MODD_CONVPAR_SHAL, ONLY : CONVPAR_SHAL
USE MODD_CONVPAREXT, ONLY : CONVPAREXT
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_T
!
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(CONVPAR_SHAL)                         ,INTENT(IN)     :: CVP_SHAL
TYPE(CONVPAREXT)                           ,INTENT(IN)     :: CVPEXT
TYPE(CST_T)                                ,INTENT(IN)     :: CST
TYPE(DIMPHYEX_T)                           ,INTENT(IN)     :: D
REAL               ,DIMENSION(D%NIT,D%NKT) ,INTENT(IN)     :: PPRES  ! pressure (P)
REAL               ,DIMENSION(D%NIT,D%NKT) ,INTENT(IN)     :: PDPRES ! pressure difference between
                                                 ! bottom and top of layer (Pa)
REAL               ,DIMENSION(D%NIT,D%NKT) ,INTENT(IN)     :: PZ     ! height of model layer (m)
REAL               ,DIMENSION(D%NIT,D%NKT) ,INTENT(IN)     :: PLMASS ! mass of model layer (kg)
REAL               ,DIMENSION(D%NIT,D%NKT) ,INTENT(IN)     :: PTHL   ! grid scale enthalpy (J/kg)
REAL               ,DIMENSION(D%NIT,D%NKT) ,INTENT(IN)     :: PTH    ! grid scale theta
REAL               ,DIMENSION(D%NIT,D%NKT) ,INTENT(IN)     :: PRW    ! grid scale total water
                                                ! mixing ratio
REAL               ,DIMENSION(D%NIT,D%NKT) ,INTENT(IN)     :: PRC    ! grid scale r_c
REAL               ,DIMENSION(D%NIT,D%NKT) ,INTENT(IN)     :: PRI    ! grid scale r_i
LOGICAL            ,DIMENSION(D%NIT)       ,INTENT(IN)     :: OTRIG1 ! logical to keep trace of
REAL               ,DIMENSION(D%NIT,D%NKT) ,INTENT(OUT)    :: PTHC  ! conv. adj. grid scale theta
REAL               ,DIMENSION(D%NIT,D%NKT) ,INTENT(OUT)    :: PRWC  ! conv. adj. grid scale r_w
REAL               ,DIMENSION(D%NIT,D%NKT) ,INTENT(OUT)    :: PRCC  ! conv. adj. grid scale r_c
REAL               ,DIMENSION(D%NIT,D%NKT) ,INTENT(OUT)    :: PRIC  ! conv. adj. grid scale r_i
REAL               ,DIMENSION(D%NIT,D%NKT) ,INTENT(OUT)    :: PWSUB ! envir. compensating subsidence(Pa/s)
INTEGER            ,DIMENSION(D%NIT)       ,INTENT(IN)     :: KLCL   ! index lifting condens. level
INTEGER            ,DIMENSION(D%NIT)       ,INTENT(IN)     :: KDPL   ! index for departure level
INTEGER            ,DIMENSION(D%NIT)       ,INTENT(IN)     :: KPBL   ! index for top of source layer
INTEGER            ,DIMENSION(D%NIT)       ,INTENT(IN)     :: KCTL   ! index for cloud top level
REAL               ,DIMENSION(D%NIT,D%NKT) ,INTENT(INOUT)  :: PUMF  ! updraft mass flux (kg/s)
REAL               ,DIMENSION(D%NIT,D%NKT) ,INTENT(INOUT)  :: PUER  ! updraft entrainment (kg/s)
REAL               ,DIMENSION(D%NIT,D%NKT) ,INTENT(INOUT)  :: PUDR  ! updraft detrainment (kg/s)
REAL               ,DIMENSION(D%NIT,D%NKT) ,INTENT(IN)     :: PUTHL  ! updraft enthalpy (J/kg)
REAL               ,DIMENSION(D%NIT,D%NKT) ,INTENT(IN)     :: PURW   ! updraft total water (kg/kg)
REAL               ,DIMENSION(D%NIT,D%NKT) ,INTENT(IN)     :: PURC   ! updraft cloud water (kg/kg)
REAL               ,DIMENSION(D%NIT,D%NKT) ,INTENT(IN)     :: PURI   ! updraft cloud ice   (kg/kg)
REAL               ,DIMENSION(D%NIT)       ,INTENT(IN)     :: PCAPE  ! available potent. energy
REAL               ,DIMENSION(D%NIT)       ,INTENT(INOUT)  :: PTIMEC ! convection time step
                                                ! convective arrays modified in UPDRAFT
!
!
INTEGER                                    ,INTENT(OUT)    :: KFTSTEPS! maximum of fract time steps
                                                 ! only used for chemical tracers
!
!
!
!
!*       0.2   Declarations of local variables :
!
INTEGER :: IKB, IKE  ! horizontal + vertical loop bounds
INTEGER :: IKS            ! vertical dimension
INTEGER :: JK, JKP, JKMAX ! vertical loop index
INTEGER :: JI             ! horizontal loop index
INTEGER :: JITER          ! iteration loop index
INTEGER :: JSTEP          ! fractional time loop index
REAL    :: ZCPORD, ZRDOCP ! C_pd / R_d, R_d / C_pd
!
REAL, DIMENSION(D%NIT,D%NKT) :: ZTHLC       ! convectively adjusted
                                          ! grid scale enthalpy
REAL, DIMENSION(D%NIT,D%NKT) :: ZOMG        ! conv. environm. subsidence (Pa/s)
REAL, DIMENSION(D%NIT,D%NKT) :: ZUMF        ! non-adjusted updraft mass flux
REAL, DIMENSION(D%NIT,D%NKT) :: ZUER        !   "     updraft  entrainm. rate
REAL, DIMENSION(D%NIT,D%NKT) :: ZUDR        !   "     updraft  detrainm. rate
REAL, DIMENSION(D%NIT)     :: ZADJ         ! mass adjustment factor
REAL, DIMENSION(D%NIT)     :: ZADJMAX      ! limit value for ZADJ
REAL, DIMENSION(D%NIT)     :: ZCAPE        ! new CAPE after adjustment
REAL, DIMENSION(D%NIT)     :: ZTIMEC       ! fractional convective time step
REAL, DIMENSION(D%NIT,D%NKT):: ZTIMC        ! 2D work array for ZTIMEC
!
REAL, DIMENSION(D%NIT)     :: ZTHLCL       ! new  theta at LCL
REAL, DIMENSION(D%NIT)     :: ZRVLCL       ! new  r_v at LCL
REAL, DIMENSION(D%NIT)     :: ZZLCL        ! height of LCL
REAL, DIMENSION(D%NIT)     :: ZTLCL        ! temperature at LCL
REAL, DIMENSION(D%NIT)     :: ZTELCL       ! envir. temper. at LCL
REAL, DIMENSION(D%NIT)     :: ZTHEUL       ! theta_e for undilute ascent
REAL, DIMENSION(D%NIT)     :: ZTHES1, ZTHES2! saturation environm. theta_e
REAL, DIMENSION(D%NIT,D%NKT) :: ZTHMFIN, ZTHMFOUT, ZRWMFIN, ZRWMFOUT
REAL, DIMENSION(D%NIT,D%NKT) :: ZRCMFIN, ZRCMFOUT, ZRIMFIN, ZRIMFOUT
                                    ! work arrays for environm. compensat. mass flux
REAL, DIMENSION(D%NIT)     :: ZPI          ! (P/P00)**R_d/C_pd
REAL, DIMENSION(D%NIT)     :: ZLV          ! latent heat of vaporisation
REAL, DIMENSION(D%NIT)     :: ZLS          ! latent heat of sublimation
REAL, DIMENSION(D%NIT)     :: ZCPH         ! specific heat C_ph
INTEGER, DIMENSION(D%NIT)  :: ITSTEP       ! fractional convective time step
INTEGER, DIMENSION(D%NIT)  :: ICOUNT       ! timestep counter
INTEGER, DIMENSION(D%NIT)  :: ILCL         ! index lifting condens. level
INTEGER, DIMENSION(D%NIT)  :: IWORK1       ! work array
REAL, DIMENSION(D%NIT)     :: ZWORK1, ZWORK2, ZWORK3, ZWORK4, ZWORK5
LOGICAL, DIMENSION(D%NIT)  :: GWORK1, GWORK3! work arrays
LOGICAL, DIMENSION(D%NIT,D%NKT) :: GWORK4    ! work array
!
!
!-------------------------------------------------------------------------------
!
!*       0.2    Initialize  local variables
!               ----------------------------
!
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "convect_closure_adjust_shal.h"

IF (LHOOK) CALL DR_HOOK('CONVECT_CLOSURE_SHAL',0,ZHOOK_HANDLE)
ZTIMC(:,:) = 0.
ZTHES2(:) = 0.
ZWORK1(:) = 0.
ZWORK2(:) = 0.
ZWORK3(:) = 0.
ZWORK4(:) = 0.
GWORK1(:) = .FALSE.
GWORK3(:) = .FALSE.
GWORK4(:,:) = .FALSE.
ILCL(:)   = KLCL(:)
!
ZCPORD    = CST%XCPD / CST%XRD
ZRDOCP    = CST%XRD / CST%XCPD
!
ZADJ(:)   = 1.
ZWORK5(:) = 1.
DO JI=D%NIB,D%NIE
  IF (.NOT. OTRIG1(JI)) ZWORK5 = 0
ENDDO
!
!
!*       0.3   Compute loop bounds
!              -------------------
!
IKB    = 1 + CVPEXT%JCVEXB
IKS    = D%NKT
IKE    = D%NKT - CVPEXT%JCVEXT
JKMAX=IKE
!
!
!*       2.     Save initial mass flux values to be used in adjustment procedure
!               ---------------------------------------------------------------
!
ZUMF(:,:)  = PUMF(:,:)
ZUER(:,:)  = PUER(:,:)
ZUDR(:,:)  = PUDR(:,:)
ZOMG(:,:)  = 0.
PWSUB(:,:) = 0.
!
!
!*       3.     Compute limits on the closure adjustment factor so that the
!               inflow in convective drafts from a given layer can't be larger
!               than the mass contained in this layer initially.
!               ---------------------------------------------------------------
!
ZADJMAX(:) = 1000.
IWORK1(:) = ILCL(:)
JKP=IKB
DO JK = JKP, IKE
  DO JI = D%NIB, D%NIE
    IF( JK > KDPL(JI))THEN
            IF(JK <= IWORK1(JI) ) THEN
        ZWORK1(JI)  = PLMASS(JI,JK) / ( ( PUER(JI,JK) + 1.E-5 ) * PTIMEC(JI) )
        ZADJMAX(JI) = MIN( ZADJMAX(JI), ZWORK1(JI) )
    END IF
    END IF
  END DO
END DO
!
!
GWORK1(D%NIB:D%NIE) = OTRIG1(D%NIB:D%NIE)  ! logical array to limit adjustment to not definitively
                                           ! adjusted columns
!
DO JK = IKB, IKE
  DO JI=D%NIB, D%NIE
    ZTHLC(JI,JK) = PTHL(JI,JK) ! initialize adjusted envir. values
    PRWC(JI,JK)  = PRW(JI,JK)
    PRCC(JI,JK)  = MAX(0., PRC(JI,JK))
    PRIC(JI,JK)  = MAX(0., PRI(JI,JK))
    PTHC(JI,JK)  = PTH(JI,JK)
  ENDDO
END DO
!
!
!
DO JITER = 1, 4  ! Enter adjustment loop to assure that all CAPE is
                 ! removed within the advective time interval TIMEC
!
     ZTIMEC(D%NIB:D%NIE) = PTIMEC(D%NIB:D%NIE)
     DO JI=1, D%NIE
     DO JK=1, IKS
       GWORK4(JI,JK) = GWORK1(JI)
     ENDDO
     ENDDO
     DO JK = IKB, IKE
       DO JI=D%NIB, D%NIE
       IF(GWORK4(JI,JK))PWSUB(JI,JK) = 0.
       ENDDO
     END DO
     ZOMG(D%NIB:D%NIE,1:D%NKT)=0.
!
     DO JK = IKB + 1, JKMAX
           JKP = MAX( IKB + 1, JK - 1 )
           DO JI=D%NIB,D%NIE
             IF(GWORK1(JI) .AND. JK <= KCTL(JI)) THEN
!
!
!*       4.     Determine vertical velocity at top and bottom of each layer
!               to satisfy mass continuity.
!               ---------------------------------------------------------------
              ! we compute here Domega/Dp = - g rho Dw/Dz = 1/Dt
!
               ZWORK1(JI)   = - ( PUER(JI,JKP) - PUDR(JI,JKP) ) / PLMASS(JI,JKP)
!
               PWSUB(JI,JK) = PWSUB(JI,JKP) - PDPRES(JI,JK-1) * ZWORK1(JI)
              ! we use PDPRES(JK-1) and not JKP in order to have zero subsidence
              ! at the first layer
!
!
!*       5.     Compute fractional time step. For stability or
!               mass conservation reasons one must split full time step PTIMEC)
!               ---------------------------------------------------------------
!
               ZWORK1(JI) = CVP_SHAL%XSTABT * PDPRES(JI,JKP) / ( ABS( PWSUB(JI,JK) ) + 1.E-10 )
              ! the factor XSTABT is used for stability reasons
               ZTIMEC(JI) = MIN( ZTIMEC(JI), ZWORK1(JI) )
!
              ! transform vertical velocity in mass flux units
               ZOMG(JI,JK) = PWSUB(JI,JK) * CVP_SHAL%XA25 / CST%XG
             ENDIF
           ENDDO
     END DO
!
!
     DO JK = IKB, IKE
       DO JI=D%NIB, D%NIE
         IF(GWORK4(JI,JK)) THEN
           ZTHLC(JI,JK) = PTHL(JI,JK) ! reinitialize adjusted envir. values
           PRWC(JI,JK)  = PRW(JI,JK)  ! when iteration criterium not attained
           PRCC(JI,JK)  = MAX(0., PRC(JI,JK))
           PRIC(JI,JK)  = MAX(0., PRI(JI,JK))
           PTHC(JI,JK)  = PTH(JI,JK)
         ENDIF
       ENDDO
     ENDDO
!
!
!        6. Check for mass conservation, i.e. ZWORK1 > 1.E-2
!           If mass is not conserved, the convective tendencies
!           automatically become zero.
!           ----------------------------------------------------
!
    DO JI = D%NIB, D%NIE
       JK=KCTL(JI)
       ZWORK1(JI) = PUDR(JI,JK) * PDPRES(JI,JK) / ( PLMASS(JI,JK) + .1 )    &
                                                            - PWSUB(JI,JK)
    END DO
    DO JI = D%NIB, D%NIE
      IF(GWORK1(JI) .AND. ABS( ZWORK1(JI) ) - .01 > 0. ) THEN
        GWORK1(JI) = .FALSE.
        PTIMEC(JI) = 1.E-1
        ZWORK5(JI) = 0.
      ENDIF
    ENDDO
    DO JK = IKB, IKE
        PWSUB(:,JK) = PWSUB(:,JK) * ZWORK5(:)
    END DO
    GWORK4(D%NIB:D%NIE,1:IKB) = .FALSE.
    GWORK4(D%NIB:D%NIE,IKE:IKS)   = .FALSE.
!
    DO JI=D%NIB,D%NIE
      ITSTEP(JI) = INT( PTIMEC(JI) / ZTIMEC(JI) ) + 1
    ENDDO
    DO JI=D%NIB,D%NIE
      ZTIMEC(JI) = PTIMEC(JI) / REAL( ITSTEP(JI) ) ! adjust  fractional time step
    ENDDO
                                           ! to be an integer multiple of PTIMEC
    DO JI=1, D%NIE
    DO JK=1, IKS
      ZTIMC(JI,JK) = ZTIMEC(JI)
    ENDDO
    ENDDO
    ICOUNT(D%NIB:D%NIE) = 0
!
!
!
    KFTSTEPS = 0
    DO JI=D%NIB,D%NIE
      KFTSTEPS = MAX(KFTSTEPS, ITSTEP(JI))
    ENDDO
    DO JSTEP = 1, KFTSTEPS ! Enter the fractional time step loop here
!
            ICOUNT(D%NIB:D%NIE) = ICOUNT(D%NIB:D%NIE) + 1
!
            GWORK3(D%NIB:D%NIE) =  ITSTEP(D%NIB:D%NIE) >= ICOUNT(D%NIB:D%NIE) .AND. GWORK1(D%NIB:D%NIE)
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
           DO JI = D%NIB, D%NIE
              GWORK4(JI,JK) = GWORK3(JI) .AND. JK <= KCTL(JI)
           END DO
           JKP = MAX( IKB + 1, JK - 1 )
           DO JI = D%NIB, D%NIE
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
           DO JI = D%NIB, D%NIE
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
         DO JK = IKB, IKE
           DO JI=D%NIB, D%NIE
             IF(GWORK4(JI,JK)) THEN
!
!******************************************************************************
!
!*       8.     Update the environmental values of enthalpy and r_w at each level
!               NOTA: These are the MAIN EQUATIONS of the scheme
!               -----------------------------------------------------------------
!
!
           ZTHLC(JI,JK) = ZTHLC(JI,JK) + ZTIMC(JI,JK) / PLMASS(JI,JK) * (      &
                          ZTHMFIN(JI,JK) + PUDR(JI,JK) * PUTHL(JI,JK)        &
                      - ZTHMFOUT(JI,JK) - PUER(JI,JK) * PTHL(JI,JK)   )
           PRWC(JI,JK)  = PRWC(JI,JK) + ZTIMC(JI,JK) / PLMASS(JI,JK) *  (      &
                          ZRWMFIN(JI,JK) + PUDR(JI,JK) * PURW(JI,JK)          &
                      - ZRWMFOUT(JI,JK) - PUER(JI,JK) * PRW(JI,JK)    )
           PRCC(JI,JK)  = PRCC(JI,JK) + ZTIMC(JI,JK) / PLMASS(JI,JK) *  (      &
               ZRCMFIN(JI,JK) + PUDR(JI,JK) * PURC(JI,JK) - ZRCMFOUT(JI,JK) -  &
                         PUER(JI,JK) * MAX(0., PRC(JI,JK))    )
           PRIC(JI,JK)  = PRIC(JI,JK) + ZTIMC(JI,JK) / PLMASS(JI,JK) *  (      &
               ZRIMFIN(JI,JK) + PUDR(JI,JK) * PURI(JI,JK) - ZRIMFOUT(JI,JK) -  &
                         PUER(JI,JK) * MAX(0., PRI(JI,JK))    )
!
!
!******************************************************************************
!
             ENDIF
           ENDDO
         ENDDO
!
    END DO ! Exit the fractional time step loop
!
!
!*          10.    Compute final linearized value of theta envir.
!                  ----------------------------------------------
!
      DO JK = IKB + 1, JKMAX
         DO JI = D%NIB, D%NIE
         IF( GWORK1(JI) .AND. JK <= KCTL(JI) ) THEN
           ZPI(JI)    = ( CST%XP00 / PPRES(JI,JK) ) ** ZRDOCP
           ZCPH(JI)   = CST%XCPD + PRWC(JI,JK) * CST%XCPV
           ZWORK2(JI) = PTH(JI,JK) / ZPI(JI)  ! first temperature estimate
           ZLV(JI)    = CST%XLVTT + ( CST%XCPV - CST%XCL ) * ( ZWORK2(JI) - CST%XTT )
           ZLS(JI)    = CST%XLVTT + ( CST%XCPV - CST%XCI ) * ( ZWORK2(JI) - CST%XTT )
             ! final linearized temperature
           ZWORK2(JI) = ( ZTHLC(JI,JK) + ZLV(JI) * PRCC(JI,JK) + ZLS(JI) * PRIC(JI,JK) &
                       - (1. + PRWC(JI,JK) ) * CST%XG * PZ(JI,JK) ) / ZCPH(JI)
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
      CALL CONVECT_CLOSURE_THRVLCL(  CVPEXT, CST, D,     &
                                     PPRES, PTHC, PRWC, PZ, GWORK1,        &
                                     ZTHLCL, ZRVLCL, ZZLCL, ZTLCL, ZTELCL, &
                                     ILCL, KDPL, KPBL )
!
!
       ZTLCL(D%NIB:D%NIE)  = MAX( 230., MIN( 335., ZTLCL(D%NIB:D%NIE) ) )  ! set some overflow bounds
       ZTELCL(D%NIB:D%NIE) = MAX( 230., MIN( 335., ZTELCL(D%NIB:D%NIE) ) )
       ZTHLCL(D%NIB:D%NIE) = MAX( 230., MIN( 345., ZTHLCL(D%NIB:D%NIE) ) )
       ZRVLCL(D%NIB:D%NIE) = MAX(   0., MIN(   1., ZRVLCL(D%NIB:D%NIE) ) )
!
!
!*         12.    Compute adjusted CAPE
!                 ---------------------
!
       ZCAPE(D%NIB:D%NIE)  = 0.
       ZPI(D%NIB:D%NIE)    = ZTHLCL(D%NIB:D%NIE) / ZTLCL(D%NIB:D%NIE)
       ZPI(D%NIB:D%NIE)    = MAX( 0.95, MIN( 1.5, ZPI(D%NIB:D%NIE) ) )
       ZWORK1(D%NIB:D%NIE) = CST%XP00 / ZPI(D%NIB:D%NIE) ** ZCPORD ! pressure at LCL
!
       CALL CONVECT_SATMIXRATIO( CST, D, ZWORK1, ZTELCL, ZWORK3, ZLV, ZLS, ZCPH )
       ZWORK3(D%NIB:D%NIE) = MIN(   .1, MAX(   0., ZWORK3(D%NIB:D%NIE) ) )
!
                ! compute theta_e updraft undilute
       ZTHEUL(D%NIB:D%NIE) = ZTLCL(D%NIB:D%NIE) * ZPI(D%NIB:D%NIE) ** ( 1. - 0.28 * ZRVLCL(D%NIB:D%NIE) )            &
                                  * EXP( ( 3374.6525 / ZTLCL(D%NIB:D%NIE) - 2.5403 )   &
                                  * ZRVLCL(D%NIB:D%NIE) * ( 1. + 0.81 * ZRVLCL(D%NIB:D%NIE) ) )
!
                ! compute theta_e saturated environment at LCL
       ZTHES1(D%NIB:D%NIE) = ZTELCL(D%NIB:D%NIE) * ZPI(D%NIB:D%NIE) ** ( 1. - 0.28 * ZWORK3(D%NIB:D%NIE) )           &
                                  * EXP( ( 3374.6525 / ZTELCL(D%NIB:D%NIE) - 2.5403 )  &
                                  * ZWORK3(D%NIB:D%NIE) * ( 1. + 0.81 * ZWORK3(D%NIB:D%NIE) ) )
!
      DO JK = IKB, JKMAX
        JKP = JK - 1
        DO JI = D%NIB, D%NIE
          ZWORK4(JI) = 1.
          IF ( JK == ILCL(JI) ) ZWORK4(JI) = 0.
!
           ! compute theta_e saturated environment and adjusted values
           ! of theta
!
          GWORK3(JI)  = JK >= ILCL(JI) .AND. JK <= KCTL(JI) .AND. GWORK1(JI)
!
          ZPI(JI)     = ( CST%XP00 / PPRES(JI,JK) ) ** ZRDOCP
          ZWORK2(JI)  = PTHC(JI,JK) / ZPI(JI)
        END DO
!
        CALL CONVECT_SATMIXRATIO( CST, D, PPRES(:,JK), ZWORK2, ZWORK3, ZLV, ZLS, ZCPH )
!
!
        DO JI = D%NIB, D%NIE
          IF ( GWORK3(JI) ) THEN
              ZTHES2(JI)  = ZWORK2(JI) * ZPI(JI) ** ( 1. - 0.28 * ZWORK3(JI) )   &
                                   * EXP( ( 3374.6525 / ZWORK2(JI) - 2.5403 ) &
                                   * ZWORK3(JI) * ( 1. + 0.81 * ZWORK3(JI) ) )
!
              ZWORK3(JI)  = PZ(JI,JK) - PZ(JI,JKP) * ZWORK4(JI) -                &
                           ( 1. - ZWORK4(JI) ) * ZZLCL(JI)    ! level thickness
              ZWORK1(JI)  = ( 2. * ZTHEUL(JI) ) / ( ZTHES1(JI) + ZTHES2(JI) ) - 1.
              ZCAPE(JI)   = ZCAPE(JI) + CST%XG * ZWORK3(JI) * MAX( 0., ZWORK1(JI) )
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
       DO JI=D%NIB,D%NIE
         IF ( GWORK1(JI) ) THEN
           ZWORK1(JI) = MAX( PCAPE(JI) - ZCAPE(JI), 0.2 * PCAPE(JI) )
           ZWORK2(JI) = ZCAPE(JI) / ( PCAPE(JI) + 1.E-8 )
!
           GWORK1(JI) = ZWORK2(JI) > 0.2 .OR. ZCAPE(JI) == 0. ! mask for adjustment
         END IF
       ENDDO
!
       DO JI=D%NIB,D%NIE
         IF( ZCAPE(JI) == 0. .AND. GWORK1(JI) )  ZADJ(JI) = ZADJ(JI) * 0.5
       ENDDO
       DO JI=D%NIB,D%NIE
         IF ( ZCAPE(JI) /= 0. .AND. GWORK1(JI) ) THEN
               ZADJ(JI) = ZADJ(JI) * CVP_SHAL%XSTABC * PCAPE(JI) / ( ZWORK1(JI) + 1.E-8 )
         ENDIF
       ENDDO
       DO JI=D%NIB,D%NIE
         ZADJ(JI) = MIN( ZADJ(JI), ZADJMAX(JI) )
       ENDDO
!
!
!*         13.     Adjust mass flux by the factor ZADJ to converge to
!                  specified degree of stabilization
!                 ----------------------------------------------------
!
       CALL CONVECT_CLOSURE_ADJUST_SHAL( CVPEXT, D, ZADJ,&
                                         PUMF, ZUMF, PUER, ZUER, PUDR, ZUDR    )
!
!
      !IF ( COUNT( GWORK1(:) ) == 0 ) EXIT ! exit big adjustment iteration loop
                                          ! when all columns have reached
                                          ! desired degree of stabilization.
!
END DO  ! end of big adjustment iteration loop
!
!
        ! skip adj. total water array  to water vapor
DO JK = IKB, IKE
  DO JI=D%NIB,D%NIE
    PRWC(JI,JK) = MAX( 0., PRWC(JI,JK) - PRCC(JI,JK) - PRIC(JI,JK) )
  END DO
END DO
!
!
IF (LHOOK) CALL DR_HOOK('CONVECT_CLOSURE_SHAL',1,ZHOOK_HANDLE)
END SUBROUTINE CONVECT_CLOSURE_SHAL

