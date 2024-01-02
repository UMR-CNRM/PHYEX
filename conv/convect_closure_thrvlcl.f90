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
      MODULE MODI_CONVECT_CLOSURE_THRVLCL
!     #################
!
INTERFACE
!
       SUBROUTINE CONVECT_CLOSURE_THRVLCL( KLON, KLEV,                         &
                                          PPRES, PTH, PRV, PZ, OWORK1,        &
                                         PTHLCL, PRVLCL, PZLCL, PTLCL, PTELCL,&
                                          KLCL, KDPL, KPBL )
!
INTEGER,                    INTENT(IN) :: KLON  ! horizontal dimension
INTEGER,                    INTENT(IN) :: KLEV  ! vertical dimension
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PTH   ! theta
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PRV   ! vapor mixing ratio 
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PPRES ! pressure
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PZ    ! height of grid point (m)
INTEGER, DIMENSION(KLON),   INTENT(IN) :: KDPL  ! contains vert. index of DPL
INTEGER, DIMENSION(KLON),   INTENT(IN) :: KPBL  ! " vert. index of source layer top
LOGICAL, DIMENSION(KLON),   INTENT(IN) :: OWORK1! logical mask 
!
REAL, DIMENSION(KLON),     INTENT(OUT):: PTHLCL ! theta at LCL
REAL, DIMENSION(KLON),     INTENT(OUT):: PRVLCL ! vapor mixing ratio at  LCL
REAL, DIMENSION(KLON),     INTENT(OUT):: PZLCL  ! height at LCL (m)
REAL, DIMENSION(KLON),     INTENT(OUT):: PTLCL  ! temperature at LCL (m)
REAL, DIMENSION(KLON),     INTENT(OUT):: PTELCL ! environm. temp. at LCL (K)
INTEGER, DIMENSION(KLON),  INTENT(OUT):: KLCL   ! contains vert. index of LCL
!
END SUBROUTINE CONVECT_CLOSURE_THRVLCL
!
END INTERFACE
!
END MODULE MODI_CONVECT_CLOSURE_THRVLCL
!     #########################################################################
      SUBROUTINE CONVECT_CLOSURE_THRVLCL( KLON, KLEV,                         &
                                          PPRES, PTH, PRV, PZ, OWORK1,        &
                                         PTHLCL, PRVLCL, PZLCL, PTLCL, PTELCL,&
                                          KLCL, KDPL, KPBL )
!     #########################################################################
!
!!**** Determine thermodynamic properties at new LCL
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to determine the thermodynamic
!!      properties at the new lifting condensation level LCL
!!   
!!
!!
!!**  METHOD
!!    ------
!!    see CONVECT_TRIGGER_FUNCT
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
!!      Module MODD_CST
!!          XG                 ! gravity constant
!!          XP00               ! Reference pressure
!!          XRD, XRV           ! Gaz  constants for dry air and water vapor
!!          XCPD               ! Cpd (dry air)
!!          XTT                ! triple point temperature
!!          XBETAW, XGAMW      ! constants for vapor saturation pressure
!!
!!      Module MODD_CONVPAR
!!          XA25               ! reference grid area
!!          XZLCL              ! lowest allowed pressure difference between
!!                             ! surface and LCL
!!          XZPBL              ! minimum mixed layer depth to sustain convection
!!          XWTRIG             ! constant in vertical velocity trigger
!!
!!      Module MODD_CONVPAREXT
!!          JCVEXB, JCVEXT     ! extra levels on the vertical boundaries
!!
!!    REFERENCE
!!    ---------
!!
!!      Book2 of documentation ( routine TRIGGER_FUNCT)
!!      Fritsch and Chappell (1980), J. Atm. Sci., Vol. 37, 1722-1761.
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    07/11/95 
!!   Last modified  04/10/97
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
USE MODD_CONVPAR
USE MODD_CONVPAREXT
USE MODI_CONVECT_SATMIXRATIO
!
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER,                    INTENT(IN) :: KLON  ! horizontal dimension
INTEGER,                    INTENT(IN) :: KLEV  ! vertical dimension
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PTH   ! theta
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PRV   ! vapor mixing ratio 
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PPRES ! pressure
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PZ    ! height of grid point (m)
INTEGER, DIMENSION(KLON),   INTENT(IN) :: KDPL  ! contains vert. index of DPL
INTEGER, DIMENSION(KLON),   INTENT(IN) :: KPBL  ! " vert. index of source layer top
LOGICAL, DIMENSION(KLON),   INTENT(IN) :: OWORK1! logical mask 
!
REAL, DIMENSION(KLON),     INTENT(OUT):: PTHLCL ! theta at LCL
REAL, DIMENSION(KLON),     INTENT(OUT):: PRVLCL ! vapor mixing ratio at  LCL
REAL, DIMENSION(KLON),     INTENT(OUT):: PZLCL  ! height at LCL (m)
REAL, DIMENSION(KLON),     INTENT(OUT):: PTLCL  ! temperature at LCL (m)
REAL, DIMENSION(KLON),     INTENT(OUT):: PTELCL ! environm. temp. at LCL (K)
INTEGER, DIMENSION(KLON),  INTENT(OUT):: KLCL   ! contains vert. index of LCL
!
!*       0.2   Declarations of local variables :
!
INTEGER :: JK, JKM, JKMIN, JKMAX      ! vertical loop index
INTEGER :: JI                         ! horizontal loop index 
INTEGER :: IIE, IKB, IKE              ! horizontal + vertical loop bounds
REAL    :: ZEPS           ! R_d / R_v
REAL    :: ZCPORD, ZRDOCP ! C_pd / R_d, R_d / C_pd
!
REAL, DIMENSION(KLON) :: ZPLCL    ! pressure at LCL
REAL, DIMENSION(KLON) :: ZTMIX    ! mixed layer temperature
REAL, DIMENSION(KLON) :: ZEVMIX   ! mixed layer water vapor pressure 
REAL, DIMENSION(KLON) :: ZDPTHMIX, ZPRESMIX ! mixed layer depth and pressure
REAL, DIMENSION(KLON) :: ZLV, ZCPH! specific heats of vaporisation, dry air
REAL, DIMENSION(KLON) :: ZDP      ! pressure between LCL and model layer
REAL, DIMENSION(KLON) :: ZWORK1, ZWORK2     ! work arrays
!
!
!-------------------------------------------------------------------------------
!
!*       0.3    Compute array bounds
!               --------------------
!
IIE = KLON
IKB = 1 + JCVEXB 
IKE = KLEV - JCVEXT 
!
!
!*       1.     Initialize local variables
!               --------------------------
!
ZEPS      = XRD / XRV
ZCPORD    = XCPD / XRD
ZRDOCP    = XRD / XCPD
!
ZDPTHMIX(:) = 0.
ZPRESMIX(:) = 0.
PTHLCL(:)   = 300.
PTLCL(:)    = 300.
PTELCL(:)   = 300.
PRVLCL(:)   = 0.
PZLCL(:)    = PZ(:,IKB)
ZTMIX(:)    = 230.
ZPLCL(:)    = 1.E4 
KLCL(:)     = IKB + 1
!
!
!*       2.     Construct a mixed layer as in TRIGGER_FUNCT
!               -------------------------------------------
!
     JKMAX = MAXVAL( KPBL(:) )
     JKMIN = MINVAL( KDPL(:) )
     DO JK = IKB + 1, JKMAX
        JKM = JK + 1
        DO JI = 1, IIE
        IF ( JK >= KDPL(JI) .AND. JK <= KPBL(JI) ) THEN
!           
            ZWORK1(JI)   = PPRES(JI,JK) - PPRES(JI,JKM)
            ZDPTHMIX(JI) = ZDPTHMIX(JI) + ZWORK1(JI)
            ZPRESMIX(JI) = ZPRESMIX(JI) + PPRES(JI,JK) * ZWORK1(JI)
            PTHLCL(JI)   = PTHLCL(JI)   + PTH(JI,JK)   * ZWORK1(JI)
            PRVLCL(JI)   = PRVLCL(JI)   + PRV(JI,JK)   * ZWORK1(JI)
!
        END IF
        END DO
     END DO
!
!
WHERE ( OWORK1(:) )
!
        ZPRESMIX(:) = ZPRESMIX(:) / ZDPTHMIX(:)
        PTHLCL(:)   = PTHLCL(:)   / ZDPTHMIX(:)
        PRVLCL(:)   = PRVLCL(:)   / ZDPTHMIX(:)
!
!*       3.1    Use an empirical direct solution ( Bolton formula )
!               to determine temperature and pressure at LCL.
!               Nota: the adiabatic saturation temperature is not
!                     equal to the dewpoint temperature
!               --------------------------------------------------
!
!
        ZTMIX(:)  = PTHLCL(:) * ( ZPRESMIX(:) / XP00 ) ** ZRDOCP
        ZEVMIX(:) = PRVLCL(:) * ZPRESMIX(:) / ( PRVLCL(:) + ZEPS )
        ZEVMIX(:) = MAX( 1.E-8, ZEVMIX(:) )
        ZWORK1(:) = ALOG( ZEVMIX(:) / 613.3 )
              ! dewpoint temperature
        ZWORK1(:) = ( 4780.8 - 32.19 * ZWORK1(:) ) / ( 17.502 - ZWORK1(:) ) 
              ! adiabatic saturation temperature
        PTLCL(:)  = ZWORK1(:) - ( .212 + 1.571E-3 * ( ZWORK1(:) - XTT )      &
                  - 4.36E-4 * ( ZTMIX(:) - XTT ) ) * ( ZTMIX(:) - ZWORK1(:) )
        PTLCL(:)  = MIN( PTLCL(:), ZTMIX(:) )
        ZPLCL(:)  = XP00 * ( PTLCL(:) / PTHLCL(:) ) ** ZCPORD
!
END WHERE
!
     ZPLCL(:) = MIN( 2.E5, MAX( 10., ZPLCL(:) ) ) ! bound to avoid overflow
!
!
!*       3.2    Correct PTLCL in order to be completely consistent
!               with MNH saturation formula
!               --------------------------------------------------
!
     CALL CONVECT_SATMIXRATIO( KLON, ZPLCL, PTLCL, ZWORK1, ZLV, ZWORK2, ZCPH )
     WHERE( OWORK1(:) )
        ZWORK2(:) = ZWORK1(:) / PTLCL(:) * ( XBETAW / PTLCL(:) - XGAMW ) ! dr_sat/dT
        ZWORK2(:) = ( ZWORK1(:) - PRVLCL(:) ) /                              &
                        ( 1. + ZLV(:) / ZCPH(:) * ZWORK2(:) ) 
        PTLCL(:)  = PTLCL(:) - ZLV(:) / ZCPH(:) * ZWORK2(:)
!
     END WHERE
!
!
!*       3.3    If PRVLCL is oversaturated set humidity and temperature
!               to saturation values.
!               -------------------------------------------------------
!
    CALL CONVECT_SATMIXRATIO( KLON, ZPRESMIX, ZTMIX, ZWORK1, ZLV, ZWORK2, ZCPH )
     WHERE( OWORK1(:) .AND. PRVLCL(:) > ZWORK1(:) )
        ZWORK2(:) = ZWORK1(:) / ZTMIX(:) * ( XBETAW / ZTMIX(:) - XGAMW ) ! dr_sat/dT
        ZWORK2(:) = ( ZWORK1(:) - PRVLCL(:) ) /                              &
                        ( 1. + ZLV(:) / ZCPH(:) * ZWORK2(:) )
        PTLCL(:)  = ZTMIX(:) + ZLV(:) / ZCPH(:) * ZWORK2(:)
        PRVLCL(:) = PRVLCL(:) - ZWORK2(:)
        ZPLCL(:)  = ZPRESMIX(:)
        PTHLCL(:) = PTLCL(:) * ( XP00 / ZPLCL(:) ) ** ZRDOCP
     END WHERE
!
!
!*        4.1   Determine  vertical loop index at the LCL 
!               -----------------------------------------
!
     DO JK = JKMIN, IKE - 1
        DO JI = 1, IIE
        IF ( ZPLCL(JI) <= PPRES(JI,JK) .AND. OWORK1(JI) ) THEN
            KLCL(JI)  = JK + 1
            PZLCL(JI) = PZ(JI,JK+1)
        END IF
        END DO
     END DO
!
!
!*        4.2   Estimate height and environmental temperature at LCL
!               ----------------------------------------------------
!
    DO JI = 1, IIE
        JK   = KLCL(JI)
        JKM  = JK - 1
        ZDP(JI)     = ALOG( ZPLCL(JI) / PPRES(JI,JKM) ) /                     &
                      ALOG( PPRES(JI,JK) / PPRES(JI,JKM) )
        ZWORK1(JI)  = PTH(JI,JK)  * ( PPRES(JI,JK) / XP00 ) ** ZRDOCP
        ZWORK2(JI)  = PTH(JI,JKM) * ( PPRES(JI,JKM) / XP00 ) ** ZRDOCP
        ZWORK1(JI)  = ZWORK2(JI) + ( ZWORK1(JI) - ZWORK2(JI) ) * ZDP(JI) 
           ! we compute the precise value of the LCL
           ! The precise height is between the levels KLCL and KLCL-1.
        ZWORK2(JI) = PZ(JI,JKM) + ( PZ(JI,JK) - PZ(JI,JKM) ) * ZDP(JI)
    END DO
    WHERE( OWORK1(:) )
       PTELCL(:) = ZWORK1(:)
       PZLCL(:)  = ZWORK2(:)
    END WHERE
!        
!
!
END SUBROUTINE CONVECT_CLOSURE_THRVLCL
