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
TYPE(CONVPAR_SHAL),        INTENT(IN) :: CVP_SHAL
TYPE(CONVPAREXT),          INTENT(IN) :: CVPEXT
TYPE(CST_T),               INTENT(IN) :: CST
TYPE(DIMPHYEX_T),          INTENT(IN) :: D
INTEGER, DIMENSION(D%NIT),  INTENT(IN) :: KLCL   ! index lifting condens. level
INTEGER, DIMENSION(D%NIT),  INTENT(IN) :: KCTL   ! index for cloud top level
INTEGER, DIMENSION(D%NIT),  INTENT(IN) :: KDPL   ! index for departure level
INTEGER, DIMENSION(D%NIT),  INTENT(IN) :: KPBL   ! index for top of source layer
REAL, DIMENSION(D%NIT),  INTENT(INOUT) :: PTIMEC ! convection time step
REAL, DIMENSION(D%NIT,D%NKT),INTENT(IN) :: PTHL   ! grid scale enthalpy (J/kg)
REAL, DIMENSION(D%NIT,D%NKT),INTENT(IN) :: PTH    ! grid scale theta
REAL, DIMENSION(D%NIT,D%NKT),INTENT(IN) :: PRW    ! grid scale total water
                                                ! mixing ratio
REAL, DIMENSION(D%NIT,D%NKT),INTENT(IN) :: PRC    ! grid scale r_c
REAL, DIMENSION(D%NIT,D%NKT),INTENT(IN) :: PRI    ! grid scale r_i
LOGICAL, DIMENSION(D%NIT),  INTENT(IN) :: OTRIG1 ! logical to keep trace of
                                                ! convective arrays modified in UPDRAFT
!
!
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN) :: PPRES  ! pressure (P)
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN) :: PDPRES ! pressure difference between
                                                 ! bottom and top of layer (Pa)
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN) :: PLMASS ! mass of model layer (kg)
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN) :: PZ     ! height of model layer (m)
REAL, DIMENSION(D%NIT),     INTENT(IN)  :: PCAPE  ! available potent. energy
INTEGER,                INTENT(OUT)   :: KFTSTEPS! maximum of fract time steps
                                                 ! only used for chemical tracers
!

REAL, DIMENSION(D%NIT,D%NKT), INTENT(INOUT):: PUMF  ! updraft mass flux (kg/s)
REAL, DIMENSION(D%NIT,D%NKT), INTENT(INOUT):: PUER  ! updraft entrainment (kg/s)
REAL, DIMENSION(D%NIT,D%NKT), INTENT(INOUT):: PUDR  ! updraft detrainment (kg/s)
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)  :: PUTHL  ! updraft enthalpy (J/kg)
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)  :: PURW   ! updraft total water (kg/kg)
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)  :: PURC   ! updraft cloud water (kg/kg)
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)  :: PURI   ! updraft cloud ice   (kg/kg)
REAL, DIMENSION(D%NIT,D%NKT), INTENT(OUT)  :: PTHC  ! conv. adj. grid scale theta
REAL, DIMENSION(D%NIT,D%NKT), INTENT(OUT)  :: PRWC  ! conv. adj. grid scale r_w
REAL, DIMENSION(D%NIT,D%NKT), INTENT(OUT)  :: PRCC  ! conv. adj. grid scale r_c
REAL, DIMENSION(D%NIT,D%NKT), INTENT(OUT)  :: PRIC  ! conv. adj. grid scale r_i
REAL, DIMENSION(D%NIT,D%NKT), INTENT(OUT)  :: PWSUB ! envir. compensating subsidence(Pa/s)
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
REAL    :: ZEPS  
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
REAL, DIMENSION(D%NIT)     :: ZWORK1, ZWORK2, ZWORK3, ZWORK5
LOGICAL, DIMENSION(D%NIT)  :: GWORK1, GWORK3! work arrays
LOGICAL, DIMENSION(D%NIT,D%NKT) :: GWORK4    ! work array
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
INTEGER           :: JLCLMIN,JDPLMIN !MIN value of LCL and DPL on all grid points,
                                     !used to remove unnecessary computations
                                     !for smaller vertical levels
INTEGER           :: JCTLMAX,JLCLMAX !MAX value of CTL and LCL on all grid points,
                                     !used to remove unnecessary computations 
                                     !for greater vertical levels
INTEGER           :: JKPP
REAL, DIMENSION(D%NIT,D%NKT) :: ZPIAUX 
!
#include "convect_closure_adjust_shal.h"
#include "convect_closure_thrvlcl.h"
!-------------------------------------------------------------------------------
!
!*       0.2    Initialize  local variables
!               ----------------------------
!
!


IF (LHOOK) CALL DR_HOOK('CONVECT_CLOSURE_SHAL',0,ZHOOK_HANDLE)
ZTIMC(:,:) = 0.
ZTHES2(:) = 0.
ZWORK1(:) = 0.
ZWORK2(:) = 0.
ZWORK3(:) = 0.
GWORK1(:) = .FALSE.
GWORK3(:) = .FALSE.
GWORK4(:,:) = .FALSE.
ILCL(:)   = KLCL(:)
!
ZEPS      = CST%XRD / CST%XRV
ZCPORD    = CST%XCPD / CST%XRD
ZRDOCP    = CST%XRD / CST%XCPD
!
ZADJ(:)   = 1.
ZWORK5(:) = 1.
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
PWSUB(:,:) = 0.
!
!
!*       3.     Compute limits on the closure adjustment factor so that the
!               inflow in convective drafts from a given layer can't be larger
!               than the mass contained in this layer initially.
!               ---------------------------------------------------------------
!
ZADJMAX(:) = 1000.
JCTLMAX=KCTL(D%NIB)   
JDPLMIN=KDPL(D%NIB)
JLCLMAX=ILCL(D%NIB)
DO JI=D%NIB,D%NIE
  IF (.NOT. OTRIG1(JI)) ZWORK5(JI) = 0
  JCTLMAX=MAX(JCTLMAX,KCTL(JI)) 
  JDPLMIN=MIN(JDPLMIN,KDPL(JI))
  JLCLMAX=MAX(JLCLMAX,ILCL(JI))
ENDDO
JCTLMAX=MIN(JCTLMAX,IKE)
JKP=IKB
DO JK = MAX(JKP,JDPLMIN+1) , MIN(IKE,JLCLMAX)
  DO JI = D%NIB, D%NIE
    IF( JK > KDPL(JI))THEN
           IF(JK <= ILCL(JI) ) THEN  
        ZWORK1(JI)  = PLMASS(JI,JK) / ( ( PUER(JI,JK) + 1.E-5 ) * PTIMEC(JI) )
        ZADJMAX(JI) = MIN( ZADJMAX(JI), ZWORK1(JI) )
    END IF
    END IF
  END DO
END DO
!
!
GWORK1(:) = OTRIG1(:)  ! logical array to limit adjustment to not definitively
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
DO JK = IKB , JCTLMAX 
  DO JI = D%NIB, D%NIE
    ZPIAUX(JI,JK)    = ( CST%XP00 / PPRES(JI,JK) ) ** ZRDOCP
  ENDDO
ENDDO

DO JITER = 1, 4  ! Enter adjustment loop to assure that all CAPE is
                 ! removed within the advective time interval TIMEC
!
     ZTIMEC(:) = PTIMEC(:)
     DO JK=1, IKS
       GWORK4(1:D%NIE,JK) = GWORK1(1:D%NIE)
     ENDDO
     DO JK = IKB, IKE
       DO JI=D%NIB, D%NIE
         IF(GWORK4(JI,JK)) PWSUB(JI,JK) = 0.
       ENDDO
     END DO
     ZOMG(:,1:D%NKT)=0.
!
     DO JK = IKB + 1, JCTLMAX 
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
!
    DO JI=D%NIB,D%NIE
      ITSTEP(JI) = INT( PTIMEC(JI) / ZTIMEC(JI) ) + 1
      ZTIMEC(JI) = PTIMEC(JI) / REAL( ITSTEP(JI) ) ! adjust  fractional time step
    ENDDO
                                           ! to be an integer multiple of PTIMEC
    DO JK=1, JCTLMAX !ZTIMC is only used for JK <= JCTLMAX in the rest of the code 
      ZTIMC(1:D%NIE,JK) = ZTIMEC(1:D%NIE)
    ENDDO
    ICOUNT(:) = 0
!
!
!
    KFTSTEPS = 0
    DO JI=D%NIB,D%NIE
      KFTSTEPS = MAX(KFTSTEPS, ITSTEP(JI))
    ENDDO
    DO JSTEP = 1, KFTSTEPS ! Enter the fractional time step loop here
!
      DO JI=D%NIB,D%NIE
        ICOUNT(JI) = ICOUNT(JI) + 1
        GWORK3(JI) = ITSTEP(JI) >= ICOUNT(JI) .AND. GWORK1(JI)
      ENDDO
!
!
!*       7.     Assign enthalpy and r_w values at the top and bottom of each
!               layer based on the sign of w
!               ------------------------------------------------------------
!
      DO JK = 1, JCTLMAX+1
        ZTHMFIN(:,JK)   = 0.
        ZRWMFIN(:,JK)   = 0.
        ZRCMFIN(:,JK)   = 0.
        ZRIMFIN(:,JK)   = 0.
        ZTHMFOUT(:,JK)  = 0.
        ZRWMFOUT(:,JK)  = 0.
        ZRCMFOUT(:,JK)  = 0.
        ZRIMFOUT(:,JK)  = 0.
      ENDDO

      IF (JKMAX-IKB <= 4) THEN
        !!!original code for the two loops
        DO JK = IKB + 1, JKMAX
          DO JI = D%NIB, D%NIE
             GWORK4(JI,JK) = GWORK3(JI) .AND. JK <= KCTL(JI)
          ENDDO
          JKP = MAX( IKB + 1, JK - 1 )
          DO JI = D%NIB, D%NIE !original loop 1
            IF ( GWORK3(JI) ) THEN  
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
            ENDIF
          ENDDO
          DO JI = D%NIB, D%NIE  !original loop 2
            IF ( GWORK3(JI) ) THEN
              ZTHMFIN(JI,JKP)  = ZTHMFIN(JI,JKP)  + ZTHMFOUT(JI,JK) * ZWORK2(JI)
              ZTHMFOUT(JI,JKP) = ZTHMFOUT(JI,JKP) + ZTHMFIN(JI,JK)  * ZWORK1(JI)
              ZRWMFIN(JI,JKP)  = ZRWMFIN(JI,JKP)  + ZRWMFOUT(JI,JK) * ZWORK2(JI)
              ZRWMFOUT(JI,JKP) = ZRWMFOUT(JI,JKP) + ZRWMFIN(JI,JK)  * ZWORK1(JI)
              ZRCMFIN(JI,JKP)  = ZRCMFIN(JI,JKP)  + ZRCMFOUT(JI,JK) * ZWORK2(JI)
              ZRCMFOUT(JI,JKP) = ZRCMFOUT(JI,JKP) + ZRCMFIN(JI,JK)  * ZWORK1(JI)
              ZRIMFIN(JI,JKP)  = ZRIMFIN(JI,JKP)  + ZRIMFOUT(JI,JK) * ZWORK2(JI)
              ZRIMFOUT(JI,JKP) = ZRIMFOUT(JI,JKP) + ZRIMFIN(JI,JK)  * ZWORK1(JI)
            ENDIF
          ENDDO
        ENDDO

      ELSE
        !!!the original two loops are rewritten to allow vectorization
        !!!if there are enough levels (JKMAX-IBK>4). The values are reported from
        !!!the original two loops to have direct expressions for 
        !!!arrays ZTHMFIN etc from arrays ZTHLC etc.,
        !!! with special cases for the first vertical
        !!!levels and for the last vertical level. The arrays ZWORK1 and ZWORK2
        !!!times ZOMG are replaced with equivalent expressions based on MAX and MIN.
        DO JK = IKB + 1, JKMAX
          DO JI = D%NIB, D%NIE
             GWORK4(JI,JK) = GWORK3(JI) .AND. JK <= KCTL(JI)
          ENDDO
        ENDDO
        IF (JKMAX .LE. JCTLMAX) THEN
          JK= JKMAX
          JKP = MAX( IKB + 1, JK - 1 )
          DO JI = D%NIB, D%NIE !if JK=JKMAX, the column JK of the arrays is only modified by 
                               !index JKMAX of original loop 1, and is not used
                               !in the remaining code if JKMAX is greater than JCTLMAX
            IF ( GWORK3(JI) ) THEN
              ZTHMFIN(JI,JK)   = - MIN(ZOMG(JI,JK),0.) * ZTHLC(JI,JKP)
              ZTHMFOUT(JI,JK)  =   MAX(ZOMG(JI,JK),0.) * ZTHLC(JI,JK)  
              ZRWMFIN(JI,JK)   = - MIN(ZOMG(JI,JK),0.) * PRWC(JI,JKP) 
              ZRWMFOUT(JI,JK)  =   MAX(ZOMG(JI,JK),0.) * PRWC(JI,JK)  
              ZRCMFIN(JI,JK)   = - MIN(ZOMG(JI,JK),0.) * PRCC(JI,JKP) 
              ZRCMFOUT(JI,JK)  =   MAX(ZOMG(JI,JK),0.) * PRCC(JI,JK)  
              ZRIMFIN(JI,JK)   = - MIN(ZOMG(JI,JK),0.) * PRIC(JI,JKP) 
              ZRIMFOUT(JI,JK)  =   MAX(ZOMG(JI,JK),0.) * PRIC(JI,JK)  
            ENDIF
          ENDDO
        ENDIF
        JK = IKB+1
        DO JI = D%NIB, D%NIE !if JK=IKB+1, the column JK of the arrays is modified by 
                             !index IKB+1 of both original loops
          IF ( GWORK3(JI) ) THEN
            ZTHMFIN(JI,JK)  =  - MIN(ZOMG(JI,JK),0.) * ZTHLC(JI,JK)   + &
               &  MAX(ZOMG(JI,JK),0.) * ZTHLC(JI,JK)               
            ZTHMFOUT(JI,JK) =  MAX(ZOMG(JI,JK),0.) * ZTHLC(JI,JK) -  &
               &  MIN(ZOMG(JI,JK),0.) * ZTHLC(JI,JK) 
            ZRWMFIN(JI,JK)  =  - MIN(ZOMG(JI,JK),0.) * PRWC(JI,JK)  + &
               &  MAX(ZOMG(JI,JK),0.) * PRWC(JI,JK)  
            ZRWMFOUT(JI,JK) =  MAX(ZOMG(JI,JK),0.) * PRWC(JI,JK)  - &
               &    MIN(ZOMG(JI,JK),0.) * PRWC(JI,JK) 
            ZRCMFIN(JI,JKP)  =  - MIN(ZOMG(JI,JK),0.) * PRCC(JI,JK)  + &
               &  MAX(ZOMG(JI,JK),0.) * PRCC(JI,JK)    
            ZRCMFOUT(JI,JK) =  MAX(ZOMG(JI,JK),0.) * PRCC(JI,JK)   - &
               &  MIN(ZOMG(JI,JK),0.) * PRCC(JI,JK)  
            ZRIMFIN(JI,JKP)  =  - MIN(ZOMG(JI,JK),0.) * PRIC(JI,JK) + &
               &  MAX(ZOMG(JI,JK),0.) * PRIC(JI,JK)  
            ZRIMFOUT(JI,JK) =  MAX(ZOMG(JI,JK),0.) * PRIC(JI,JK)  - &
               &  MIN(ZOMG(JI,JK),0.) * PRIC(JI,JK) 
          ENDIF
        ENDDO
        JK = IKB + 2
        JKP= JK - 1  !JKP = MAX (IKB + 1 , JK - 1)
        DO JI = D%NIB, D%NIE !if JK=IBK+2, the column JKB+1 of the arrays is modified 
                             !again by index IKB+2 of original loop 2
          IF ( GWORK3(JI) ) THEN
            ZTHMFIN(JI,JKP)  = ZTHMFIN(JI,JKP)  + &
              & MAX(ZOMG(JI,JK),0.) * ZTHLC(JI,JK)   
            ZTHMFOUT(JI,JKP) = ZTHMFOUT(JI,JKP) - &
              & MIN(ZOMG(JI,JK),0.) * ZTHLC(JI,JKP)  
            ZRWMFIN(JI,JKP)  = ZRWMFIN(JI,JKP)  + &
              & MAX(ZOMG(JI,JK),0.) * PRWC(JI,JK)  
            ZRWMFOUT(JI,JKP) = ZRWMFOUT(JI,JKP) - &
              & MIN(ZOMG(JI,JK),0.) * PRWC(JI,JKP) 
            ZRCMFIN(JI,JKP)  = ZRCMFIN(JI,JKP)  + &
              & MAX(ZOMG(JI,JK),0.) * PRCC(JI,JK)  
            ZRCMFOUT(JI,JKP) = ZRCMFOUT(JI,JKP) - &
              & MIN(ZOMG(JI,JK),0.) * PRCC(JI,JKP) 
            ZRIMFIN(JI,JKP)  = ZRIMFIN(JI,JKP)  + &
              & MAX(ZOMG(JI,JK),0.) * PRIC(JI,JK)  
            ZRIMFOUT(JI,JKP) = ZRIMFOUT(JI,JKP) - &
              & MIN(ZOMG(JI,JK),0.) * PRIC(JI,JKP) 
          ENDIF
        ENDDO
        DO JK = IKB + 3, MIN(JKMAX-1,JCTLMAX+1) 
          JKP= JK - 1 !JKP = MAX( IKB + 1, JK - 1 )
          JKPP=JK - 2 !JKPP= MAX( IKB + 1, JK - 2 )
          DO JI = D%NIB, D%NIE !if JK is greater than IKB+3, the column JK of the arrays is modified 
                               !by index JK of original loop 1 and index JK+1 of original loop 2
            IF ( GWORK3(JI) ) THEN
              ZTHMFIN(JI,JKP)  =  - MIN(ZOMG(JI,JKP),0.) * ZTHLC(JI,JKPP) + &
                 &  MAX(ZOMG(JI,JK),0.) * ZTHLC(JI,JK)   
              ZTHMFOUT(JI,JKP) =  MAX(ZOMG(JI,JKP),0.) * ZTHLC(JI,JKP)   - &
                 &   MIN(ZOMG(JI,JK),0.) * ZTHLC(JI,JKP)  
              ZRWMFIN(JI,JKP)  =  - MIN(ZOMG(JI,JKP),0.) * PRWC(JI,JKPP) + &
                 &  MAX(ZOMG(JI,JK),0.) * PRWC(JI,JK)   
              ZRWMFOUT(JI,JKP) =  MAX(ZOMG(JI,JKP),0.) * PRWC(JI,JKP)  - &
                 &  MIN(ZOMG(JI,JK),0.) * PRWC(JI,JKP) 
              ZRCMFIN(JI,JKP)  =  - MIN(ZOMG(JI,JKP),0.) * PRCC(JI,JKPP) + &
                 &  MAX(ZOMG(JI,JK),0.) * PRCC(JI,JK)   
              ZRCMFOUT(JI,JKP) =  MAX(ZOMG(JI,JKP),0.) * PRCC(JI,JKP)    - &
                 &  MIN(ZOMG(JI,JK),0.) * PRCC(JI,JKP) 
              ZRIMFIN(JI,JKP)  =  - MIN(ZOMG(JI,JKP),0.) * PRIC(JI,JKPP) + &
                 &  MAX(ZOMG(JI,JK),0.) * PRIC(JI,JK)   
              ZRIMFOUT(JI,JKP) =  MAX(ZOMG(JI,JKP),0.) * PRIC(JI,JKP) - &
                 &  MIN(ZOMG(JI,JK),0.) * PRIC(JI,JKP)  
            ENDIF
          ENDDO
        ENDDO

      ENDIF

!
!******************************************************************************
!
!*       8.     Update the environmental values of enthalpy and r_w at each level
!               NOTA: These are the MAIN EQUATIONS of the scheme
!               -----------------------------------------------------------------
!
!
      DO JK = IKB, JCTLMAX !GWORK4 is .FALSE. if JK is greater than JCTLMAX 
        DO JI=D%NIB, D%NIE
          IF(GWORK4(JI,JK)) THEN
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
          ENDIF
        ENDDO
      ENDDO
!
!
!******************************************************************************
!
!
    ENDDO ! Exit the fractional time step loop
!
!
!*          10.    Compute final linearized value of theta envir.
!                  ----------------------------------------------
!
    DO JK = IKB + 1,JCTLMAX 
      DO JI = D%NIB, D%NIE
        IF( GWORK1(JI) .AND. JK <= KCTL(JI) ) THEN

          ZCPH(JI)   = CST%XCPD + PRWC(JI,JK) * CST%XCPV
          ZWORK2(JI) = PTH(JI,JK) / ZPIAUX(JI,JK)  ! first temperature estimate
          ZLV(JI)    = CST%XLVTT + ( CST%XCPV - CST%XCL ) * ( ZWORK2(JI) - CST%XTT )
          ZLS(JI)    = CST%XLVTT + ( CST%XCPV - CST%XCI ) * ( ZWORK2(JI) - CST%XTT )
          ! final linearized temperature
          ZWORK2(JI) = ( ZTHLC(JI,JK) + ZLV(JI) * PRCC(JI,JK) + ZLS(JI) * PRIC(JI,JK) &
                      - (1. + PRWC(JI,JK) ) * CST%XG * PZ(JI,JK) ) / ZCPH(JI)
          ZWORK2(JI) = MAX( 180., MIN( 340., ZWORK2(JI) ) )
          PTHC(JI,JK)= ZWORK2(JI) * ZPIAUX(JI,JK) ! final adjusted envir. theta
        ENDIF
      ENDDO
    ENDDO
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
JLCLMIN=ILCL(D%NIB)  

!DIR$ IVDEP
    DO JI = D%NIB, D%NIE
      JLCLMIN=MIN(JLCLMIN,ILCL(JI)) 
      ZTLCL (JI)  = MAX( 230., MIN( 335., ZTLCL(JI) ) )  ! set some overflow bounds
      ZTELCL(JI) = MAX( 230., MIN( 335., ZTELCL(JI) ) )
      ZTHLCL(JI) = MAX( 230., MIN( 345., ZTHLCL(JI) ) )
      ZRVLCL(JI) = MAX(   0., MIN(   1., ZRVLCL(JI) ) )
      ZCAPE (JI) = 0.
      ZPI   (JI) = MAX( 0.95, MIN( 1.5, ZTHLCL(JI)/ZTLCL(JI) ) )
      ZWORK1(JI) = CST%XP00 / ZPI(JI) ** ZCPORD ! pressure at LCL
      CALL CONVECT_SATMIXRATIO( ZWORK1(JI), ZTELCL(JI), ZEPS, ZWORK3(JI), ZLV(JI), ZLS(JI), ZCPH(JI) )
!
!*         12.    Compute adjusted CAPE
!                 ---------------------
      ZWORK3(JI) = MIN(   .1, MAX(   0., ZWORK3(JI) ) )


      ! compute theta_e updraft undilute
      ZTHEUL(JI) = ZTLCL(JI) * ZPI(JI) ** ( 1. - 0.28 * ZRVLCL(JI) )            &
       & * EXP( ( 3374.6525 / ZTLCL(JI) - 2.5403 ) * ZRVLCL(JI) * ( 1. + 0.81 * ZRVLCL(JI) ) )
      ! compute theta_e saturated environment at LCL
      ZTHES1(JI) = ZTELCL(JI) * ZPI(JI) ** ( 1. - 0.28 * ZWORK3(JI) )           &
       & * EXP( ( 3374.6525 / ZTELCL(JI) - 2.5403 ) * ZWORK3(JI) * ( 1. + 0.81 * ZWORK3(JI) ) )
    ENDDO
!
    JLCLMIN=MAX(JLCLMIN,IKB)        
    DO JK = JLCLMIN, JCTLMAX  
      JKP = JK - 1
!DIR$ IVDEP
      DO JI = D%NIB, D%NIE
        ! compute theta_e saturated environment and adjusted values of theta
        GWORK3(JI)  = JK >= ILCL(JI) .AND. JK <= KCTL(JI) .AND. GWORK1(JI)
        ZWORK2(JI)  = PTHC(JI,JK) / ZPIAUX(JI,JK)
        CALL CONVECT_SATMIXRATIO( PPRES(JI,JK), ZWORK2(JI), ZEPS, ZWORK3(JI), ZLV(JI), ZLS(JI), ZCPH(JI) )
      END DO
!
      DO JI = D%NIB, D%NIE 
        IF ( GWORK3(JI) ) THEN
           ZTHES2(JI)  = ZWORK2(JI) * ZPIAUX(JI,JK) ** ( 1. - 0.28 * ZWORK3(JI) )   &
           & * EXP( ( 3374.6525 / ZWORK2(JI) - 2.5403 ) * ZWORK3(JI) * ( 1. + 0.81 * ZWORK3(JI) ) )
          IF ( JK == ILCL(JI) ) THEN
            ZWORK3(JI)  = PZ(JI,JK) - ZZLCL(JI)  ! level thickness
          ELSE
            ZWORK3(JI)  = PZ(JI,JK) - PZ(JI,JKP) ! level thickness
          ENDIF
          ZWORK1(JI)  = ( 2. * ZTHEUL(JI) ) / ( ZTHES1(JI) + ZTHES2(JI) ) - 1.
          ZCAPE(JI)   = ZCAPE(JI) + CST%XG * ZWORK3(JI) * MAX( 0., ZWORK1(JI) )
          ZTHES1(JI)  = ZTHES2(JI)
        ENDIF
      ENDDO

    ENDDO
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
      ENDIF
    ENDDO
!
    DO JI=D%NIB,D%NIE
      IF( ZCAPE(JI) == 0. .AND. GWORK1(JI) ) THEN
        ZADJ(JI) = ZADJ(JI) * 0.5
      ENDIF
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
    CALL CONVECT_CLOSURE_ADJUST_SHAL( CVPEXT, D, ZADJ, PUMF, ZUMF, PUER, ZUER, PUDR, ZUDR    )
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
CONTAINS
INCLUDE "convect_satmixratio.h"
!
END SUBROUTINE CONVECT_CLOSURE_SHAL
