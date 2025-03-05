!MNH_LIC Copyright 2018-2024 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_LIMA_GRAUPEL
  IMPLICIT NONE
CONTAINS
!     #################################################################################
  SUBROUTINE LIMA_GRAUPEL (CST, LIMAP, LIMAC, LIMAM, KSIZE, PTSTEP, ODCOMPUTE,                              &
                           PRHODREF, PPRES, PT, PKA, PDV, PCJ,                    &
                           PRVT, PRCT, PRRT, PRIT, PRST, PRGT,                    &
                           PCCT, PCRT, PCIT, PCST, PCGT,                          &
                           PLBDC, PLBDR, PLBDS, PLBDG,                            &
                           PLVFACT, PLSFACT,                                      &
                           P_TH_WETG, P_RC_WETG, P_CC_WETG, P_RR_WETG, P_CR_WETG, &
                           P_RI_WETG, P_CI_WETG, P_RS_WETG, P_CS_WETG, P_RG_WETG, P_CG_WETG, P_RH_WETG, &
                           P_TH_DRYG, P_RC_DRYG, P_CC_DRYG, P_RR_DRYG, P_CR_DRYG, &
                           P_RI_DRYG, P_CI_DRYG, P_RS_DRYG, P_CS_DRYG, P_RG_DRYG, &
                           P_RI_HMG, P_CI_HMG, P_RG_HMG,                          &
                           P_TH_GMLT, P_RR_GMLT, P_CR_GMLT, P_CG_GMLT,            &
                           PA_TH, PA_RC, PA_CC, PA_RR, PA_CR,                     &
                           PA_RI, PA_CI, PA_RS, PA_CS, PA_RG, PA_CG, PA_RH, PA_CH )
!     #################################################################################
!
!!    PURPOSE
!!    -------
!!      Compute the wet/dry growth of graupel, associated Hallett-Mossop ice production,
!!    and graupel melting rates
!!
!!    AUTHOR
!!    ------
!!      J.-M. Cohard     * Laboratoire d'Aerologie*
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!      S.    Berthet    * Laboratoire d'Aerologie*
!!      B.    Vié        * CNRM *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original             15/03/2018
!  P. Wautelet 26/04/2019: replace non-standard FLOAT function by REAL function
!  J. Wurtz       03/2022: new snow characteristics
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAM_LIMA_MIXED, ONLY:PARAM_LIMA_MIXED_T
USE MODD_PARAM_LIMA_COLD, ONLY:PARAM_LIMA_COLD_T
USE MODD_PARAM_LIMA, ONLY:PARAM_LIMA_T
USE MODD_CST, ONLY:CST_T
USE YOMHOOK, ONLY:LHOOK, DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(PARAM_LIMA_MIXED_T),INTENT(IN)::LIMAM
TYPE(PARAM_LIMA_COLD_T),INTENT(IN)::LIMAC
TYPE(PARAM_LIMA_T),INTENT(IN)::LIMAP
TYPE(CST_T),INTENT(IN)::CST
INTEGER, INTENT(IN) :: KSIZE
REAL,                 INTENT(IN)    :: PTSTEP 
LOGICAL, DIMENSION(KSIZE),INTENT(IN)    :: ODCOMPUTE
!
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRHODREF    ! 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PPRES    ! 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PT   ! 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PKA   ! 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PDV   ! 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PCJ   ! 
!
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRVT    ! 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRCT    ! 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRRT    ! 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRIT    ! 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRST    ! 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRGT    !
!
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PCCT    !
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PCRT    !
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PCIT    !
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PCST    !
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PCGT    !
!
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLBDC   ! 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLBDR   ! 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLBDS   ! 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLBDG   ! 
!
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLVFACT ! 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLSFACT ! 
!
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_TH_WETG
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_RC_WETG
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_CC_WETG
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_RR_WETG
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_CR_WETG
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_RI_WETG
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_CI_WETG
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_RS_WETG
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_CS_WETG
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_RG_WETG
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_CG_WETG
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_RH_WETG
!
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_TH_DRYG
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_RC_DRYG
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_CC_DRYG
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_RR_DRYG
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_CR_DRYG
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_RI_DRYG
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_CI_DRYG
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_RS_DRYG
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_CS_DRYG
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_RG_DRYG
!
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_RI_HMG
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_CI_HMG
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_RG_HMG
!
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_TH_GMLT
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_RR_GMLT
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_CR_GMLT
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_CG_GMLT
!
REAL, DIMENSION(KSIZE),   INTENT(INOUT) :: PA_TH
REAL, DIMENSION(KSIZE),   INTENT(INOUT) :: PA_RC
REAL, DIMENSION(KSIZE),   INTENT(INOUT) :: PA_CC
REAL, DIMENSION(KSIZE),   INTENT(INOUT) :: PA_RR
REAL, DIMENSION(KSIZE),   INTENT(INOUT) :: PA_CR
REAL, DIMENSION(KSIZE),   INTENT(INOUT) :: PA_RI
REAL, DIMENSION(KSIZE),   INTENT(INOUT) :: PA_CI
REAL, DIMENSION(KSIZE),   INTENT(INOUT) :: PA_RS
REAL, DIMENSION(KSIZE),   INTENT(INOUT) :: PA_CS
REAL, DIMENSION(KSIZE),   INTENT(INOUT) :: PA_RG
REAL, DIMENSION(KSIZE),   INTENT(INOUT) :: PA_CG
REAL, DIMENSION(KSIZE),   INTENT(INOUT) :: PA_RH
REAL, DIMENSION(KSIZE),   INTENT(INOUT) :: PA_CH
!
!*       0.2   Declarations of local variables :
!
LOGICAL, DIMENSION(SIZE(PRCT))  :: GDRY
!
REAL,    DIMENSION(SIZE(PRCT))  :: Z1, Z2, Z3, Z4
REAL,    DIMENSION(SIZE(PRCT))  :: ZZX, ZZW, ZZW1, ZZW2, ZZW3, ZZW4, ZZW5, ZZW6, ZZW7
REAL,    DIMENSION(SIZE(PRCT))  :: ZZW3N, ZZW4N, ZZW6N 
REAL,    DIMENSION(SIZE(PRCT))  :: ZRDRYG, ZRWETG
REAL,    DIMENSION(SIZE(PRCT))  :: ZSIGMOIDE
!
INTEGER, DIMENSION(SIZE(PRCT))  :: IVEC1,IVEC2        ! Vectors of indices
REAL,    DIMENSION(SIZE(PRCT))  :: ZVEC1,ZVEC2, ZVEC3 ! Work vectors
!
INTEGER                         :: IHAIL
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
!
IF (LHOOK) CALL DR_HOOK('LIMA_GRAUPEL', 0, ZHOOK_HANDLE)
P_TH_WETG(:) = 0.
P_RC_WETG(:) = 0.
P_CC_WETG(:) = 0.
P_RR_WETG(:) = 0.
P_CR_WETG(:) = 0.
P_RI_WETG(:) = 0.
P_CI_WETG(:) = 0.
P_RS_WETG(:) = 0.
P_CS_WETG(:) = 0.
P_RG_WETG(:) = 0.
P_CG_WETG(:) = 0.
P_RH_WETG(:) = 0.
!
P_TH_DRYG(:) = 0.
P_RC_DRYG(:) = 0.
P_CC_DRYG(:) = 0.
P_RR_DRYG(:) = 0.
P_CR_DRYG(:) = 0.
P_RI_DRYG(:) = 0.
P_CI_DRYG(:) = 0.
P_RS_DRYG(:) = 0.
P_CS_DRYG(:) = 0.
P_RG_DRYG(:) = 0.
!
P_RI_HMG(:) = 0.
P_CI_HMG(:) = 0.
P_RG_HMG(:) = 0.
!
P_TH_GMLT(:) = 0.
P_RR_GMLT(:) = 0.
P_CR_GMLT(:) = 0.
P_CG_GMLT(:) = 0.
!
ZZW1(:) = 0. ! RCDRYG
ZZW2(:) = 0. ! RIDRYG
ZZW3(:) = 0. ! RSDRYG
ZZW3N(:) = 0.! NSDRYG
ZZW4(:) = 0. ! RRDRYG
ZZW4N(:) = 0.! NRDRYG
ZZW5(:) = 0. ! RIWETG
ZZW6(:) = 0. ! RSWETG
ZZW7(:) = 0. ! 
!
ZRDRYG(:) = 0.
ZRWETG(:) = 0.
!
!
!*       1.  Graupel growth by collection (dry or wet case)
!        --------------------------------------------------
!
!            1.a Collection of rc and ri in the dry mode
!            --------------------------------------------
!
WHERE( PRGT(:)>LIMAP%XRTMIN(6) .AND. PCGT(:)>LIMAP%XCTMIN(6) .AND. ODCOMPUTE(:) )
   ZZW(:) = PCGT(:) * PLBDG(:)**(-LIMAM%XDG-2.0) * PRHODREF(:)**(1-LIMAP%XCEXVT)
   ZZW1(:) = LIMAM%XFCDRYG * PRCT(:) * ZZW(:)                               ! RCDRYG - rc collected by graupel in dry mode 
   ZZW2(:) = LIMAM%XFIDRYG * EXP( LIMAM%XCOLEXIG*(PT(:)-CST%XTT) ) * PRIT(:) * ZZW(:) ! RIDRYG - ri collected by graupel in dry mode
END WHERE
!
!*           1.b Collection of rs in the dry mode
!            ------------------------------------
!
GDRY(:) = PRST(:)>LIMAP%XRTMIN(5) .AND. PCST(:)>LIMAP%XCTMIN(5) .AND. PRGT(:)>LIMAP%XRTMIN(6) &
     .AND. PCGT(:)>LIMAP%XCTMIN(6) .AND. ODCOMPUTE(:)
!
WHERE( GDRY )
!
!*       Select the (ZLBDAG,ZLBDAS) couplet
!
   ZVEC1(:) = PLBDG(:)
   ZVEC2(:) = PLBDS(:)
!
!*       find the next lower indice for the ZLBDAG and for the ZLBDAS
!        in the geometrical set of (Lbda_g,Lbda_s) couplet use to
!        tabulate the SDRYG-kernel
!
   ZVEC1(:) = MAX( 1.0001, MIN( REAL(LIMAM%NDRYLBDAG)-0.0001,           &
                         LIMAM%XDRYINTP1G * LOG( ZVEC1(:) ) + LIMAM%XDRYINTP2G ) )
   IVEC1(:) = INT( ZVEC1(:) )
   ZVEC1(:) = ZVEC1(:) - REAL( IVEC1(:) )
!
   ZVEC2(:) = MAX( 1.0001, MIN( REAL(LIMAM%NDRYLBDAS)-0.0001,           &
                         LIMAM%XDRYINTP1S * LOG( ZVEC2(:) ) + LIMAM%XDRYINTP2S ) )
   IVEC2(:) = INT( ZVEC2(:) )
   ZVEC2(:) = ZVEC2(:) - REAL( IVEC2(:) )
!
!*       perform the bilinear interpolation of the normalized
!        SDRYG-kernel
   !
   Z1(:) = GET_XKER_SDRYG(KSIZE,IVEC1(:)+1,IVEC2(:)+1)
   Z2(:) = GET_XKER_SDRYG(KSIZE,IVEC1(:)+1,IVEC2(:)  )
   Z3(:) = GET_XKER_SDRYG(KSIZE,IVEC1(:)  ,IVEC2(:)+1)
   Z4(:) = GET_XKER_SDRYG(KSIZE,IVEC1(:)  ,IVEC2(:)  )
   ZVEC3(:) =  (      Z1(:)* ZVEC2(:)          &
                    - Z2(:)*(ZVEC2(:) - 1.0) ) &
                                               *  ZVEC1(:)    &
                 - (  Z3(:)* ZVEC2(:)          &
                    - Z4(:)*(ZVEC2(:) - 1.0) ) &
                                                    * (ZVEC1(:) - 1.0)
   ZZW(:) = ZVEC3(:)
   ZZW3(:) = LIMAM%XFSDRYG * ZZW(:) * EXP( LIMAM%XCOLEXSG*(PT(:)-CST%XTT) )       & ! RSDRYG - rs collected by graupel in dry mode
                    *  PRST(:) * PCGT(:)                          &
                    *  PRHODREF(:)**(1-LIMAP%XCEXVT)                    &
                    *( LIMAM%XLBSDRYG1/( PLBDG(:)**2                ) + &
                       LIMAM%XLBSDRYG2/( PLBDG(:)    * PLBDS(:)     ) + &
                       LIMAM%XLBSDRYG3/(               PLBDS(:)**2) )
!
   Z1(:) = GET_XKER_N_SDRYG(KSIZE,IVEC1(:)+1,IVEC2(:)+1)
   Z2(:) = GET_XKER_N_SDRYG(KSIZE,IVEC1(:)+1,IVEC2(:)  )
   Z3(:) = GET_XKER_N_SDRYG(KSIZE,IVEC1(:)  ,IVEC2(:)+1)
   Z4(:) = GET_XKER_N_SDRYG(KSIZE,IVEC1(:)  ,IVEC2(:)  )
   ZVEC3(:) =  (      Z1(:)* ZVEC2(:)          &
                    - Z2(:)*(ZVEC2(:) - 1.0) ) &
                                               *  ZVEC1(:)    &
                 - (  Z3(:)* ZVEC2(:)          &
                    - Z4(:)*(ZVEC2(:) - 1.0) ) &
                                                    * (ZVEC1(:) - 1.0)
   ZZW(:) = ZVEC3(:)
   ZZW3N(:) = LIMAM%XFNSDRYG * ZZW(:) * EXP( LIMAM%XCOLEXSG*(PT(:)-CST%XTT) )        & ! NSDRYG - Ns collected by graupel in dry mode
                      *  PCST(:) * PCGT(:)                           &
                      *  PRHODREF(:)**(1-LIMAP%XCEXVT)                     &
                      *( LIMAM%XLBNSDRYG1/( PLBDG(:)**2                ) + &
                         LIMAM%XLBNSDRYG2/( PLBDG(:)    * PLBDS(:)     ) + &
                         LIMAM%XLBNSDRYG3/(               PLBDS(:)**2) )
END WHERE
!
!*           1.c  Collection of rr in the dry mode
!            -------------------------------------
!
GDRY(:) = PRRT(:)>LIMAP%XRTMIN(3) .AND. PCRT(:)>LIMAP%XCTMIN(3) .AND. PRGT(:)>LIMAP%XRTMIN(6) &
     .AND. PCGT(:)>LIMAP%XCTMIN(6) .AND. ODCOMPUTE(:)
!
WHERE( GDRY )
!
!*       Select the (ZLBDAG,ZLBDAR) couplet
!
   ZVEC1(:) = PLBDG(:)
   ZVEC2(:) = PLBDR(:)
!
!*       Find the next lower indice for the ZLBDAG and for the ZLBDAR
!        in the geometrical set of (Lbda_g,Lbda_r) couplet use to
!        tabulate the RDRYG-kernel
!
   ZVEC1(:) = MAX( 1.0001, MIN( REAL(LIMAM%NDRYLBDAG)-0.0001,           &
                         LIMAM%XDRYINTP1G * LOG( ZVEC1(:) ) + LIMAM%XDRYINTP2G ) )
   IVEC1(:) = INT( ZVEC1(:) )
   ZVEC1(:) = ZVEC1(:) - REAL( IVEC1(:) )
!
   ZVEC2(:) = MAX( 1.0001, MIN( REAL(LIMAM%NDRYLBDAR)-0.0001,           &
                         LIMAM%XDRYINTP1R * LOG( ZVEC2(:) ) + LIMAM%XDRYINTP2R ) )
   IVEC2(:) = INT( ZVEC2(:) )
   ZVEC2(:) = ZVEC2(:) - REAL( IVEC2(:) )
!
!*       Perform the bilinear interpolation of the normalized
!        RDRYG-kernel
!
   Z1(:) = GET_XKER_RDRYG(KSIZE,IVEC1(:)+1,IVEC2(:)+1)
   Z2(:) = GET_XKER_RDRYG(KSIZE,IVEC1(:)+1,IVEC2(:)  )
   Z3(:) = GET_XKER_RDRYG(KSIZE,IVEC1(:)  ,IVEC2(:)+1)
   Z4(:) = GET_XKER_RDRYG(KSIZE,IVEC1(:)  ,IVEC2(:)  )
      ZVEC3(:) =  (   Z1(:)* ZVEC2(:)          &
                    - Z2(:)*(ZVEC2(:) - 1.0) ) &
                                               *  ZVEC1(:)   &
                 - (  Z3(:)* ZVEC2(:)          &
                    - Z4(:)*(ZVEC2(:) - 1.0) ) &
                                               * (ZVEC1(:) - 1.0)
   ZZW(:) = ZVEC3(:)
   ZZW4(:) = LIMAM%XFRDRYG * ZZW(:)                                     & ! RRDRYG
                     * PRRT(:) * PCGT(:)                          &
                     * PRHODREF(:)**(1-LIMAP%XCEXVT)                    &
                     *( LIMAM%XLBRDRYG1/( PLBDG(:)**2               ) + &
                        LIMAM%XLBRDRYG2/( PLBDG(:)   * PLBDR(:)     ) + &
                        LIMAM%XLBRDRYG3/(              PLBDR(:)**2) )
!
   Z1(:) = GET_XKER_N_RDRYG(KSIZE,IVEC1(:)+1,IVEC2(:)+1)
   Z2(:) = GET_XKER_N_RDRYG(KSIZE,IVEC1(:)+1,IVEC2(:)  )
   Z3(:) = GET_XKER_N_RDRYG(KSIZE,IVEC1(:)  ,IVEC2(:)+1)
   Z4(:) = GET_XKER_N_RDRYG(KSIZE,IVEC1(:)  ,IVEC2(:)  )
      ZVEC3(:) =  (   Z1(:)* ZVEC2(:)          &
                    - Z2(:)*(ZVEC2(:) - 1.0) ) &
                                               *  ZVEC1(:)   &
                 - (  Z3(:)* ZVEC2(:)          &
                    - Z4(:)*(ZVEC2(:) - 1.0) ) &
                                               * (ZVEC1(:) - 1.0)
   ZZW(:) = ZVEC3(:)
   ZZW4N(:) = LIMAM%XFNRDRYG * ZZW(:)                                      & ! NRDRYG
                       * PCRT(:) * PCGT(:)                           &
                       * PRHODREF(:)**(1-LIMAP%XCEXVT)                     &
                       *( LIMAM%XLBNRDRYG1/( PLBDG(:)**2               ) + &
                          LIMAM%XLBNRDRYG2/( PLBDG(:)   * PLBDR(:)     ) + &
                          LIMAM%XLBNRDRYG3/(              PLBDR(:)**2) )
   
END WHERE
!
!            1.d Total collection in the dry mode
!            ------------------------------------
!
ZRDRYG(:) = ZZW1(:) + ZZW2(:) + ZZW3(:) + ZZW4(:)
!
!            1.e Collection in the wet mode
!            ------------------------------
!
ZZW(:) = 0.0
WHERE( PRGT(:)>LIMAP%XRTMIN(6) .AND. PCGT(:)>LIMAP%XCTMIN(6) .AND. ODCOMPUTE(:) )
   ZZW5(:) = ZZW2(:) / (LIMAM%XCOLIG*EXP(LIMAM%XCOLEXIG*(PT(:)-CST%XTT)) ) ! RIWETG
   ZZW6(:) = ZZW3(:) / (LIMAM%XCOLSG*EXP(LIMAM%XCOLEXSG*(PT(:)-CST%XTT)) ) ! RSWETG
   ZZW6N(:)= ZZW3N(:)/ (LIMAM%XCOLSG*EXP(LIMAM%XCOLEXSG*(PT(:)-CST%XTT)) ) ! NSWETG
!
   ZZW(:) = PRVT(:)*PPRES(:)/((CST%XMV/CST%XMD)+PRVT(:)) ! Vapor pressure
   ZZW(:) = PKA(:)*(CST%XTT-PT(:)) +                                  &
              ( PDV(:)*(CST%XLVTT + ( CST%XCPV - CST%XCL ) * ( PT(:) - CST%XTT ))   &
                          *(CST%XESTT-ZZW(:))/(CST%XRV*PT(:))             )
!
! Total mass gained by graupel in wet mode
   ZRWETG(:)  = MAX( 0.0,                                                        &
                   ( ZZW(:) * PCGT(:) * ( LIMAM%X0DEPG*       PLBDG(:)**LIMAM%XEX0DEPG +     &
                                          LIMAM%X1DEPG*PCJ(:)*PLBDG(:)**LIMAM%XEX1DEPG ) +   &
                   ( ZZW5(:)+ZZW6(:) ) * ( CST%XLMTT + (CST%XCI-CST%XCL)*(CST%XTT-PT(:)) )   )   &
                   / (CST%XLMTT-CST%XCL*(CST%XTT-PT(:)))                                     )
  !We must agregate, at least, the cold species
   ZRWETG(:)=MAX(ZRWETG(:), ZZW5(:)+ZZW6(:))
END WHERE
!
!            1.f Wet mode and partial conversion to hail
!            -------------------------------------------
!
ZZW(:) = 0.0
IHAIL = 0.
IF (LIMAP%NMOM_H.GE.1) IHAIL = 1. 
WHERE( ODCOMPUTE(:) .AND. PRGT(:)>LIMAP%XRTMIN(6) .AND. PCGT(:)>LIMAP%XCTMIN(6) .AND. PT(:)<CST%XTT .AND. &
       (ZRDRYG(:)-ZZW2(:)-ZZW3(:))>=(ZRWETG(:)-ZZW5(:)-ZZW6(:)) .AND. ZRWETG(:)-ZZW5(:)-ZZW6(:)>0.0 ) 
!
! Mass of rain and cloud droplets frozen by graupel in wet mode : RCWETG + RRWETG = RWETG - RIWETG - RSWETG
   ZZW7(:) = ZRWETG(:) - ZZW5(:) - ZZW6(:)
!
! assume a linear percent of conversion of graupel into hail
! ZZW = percentage of graupel transformed
!
   ZZW(:) = ZRDRYG(:)*IHAIL/(ZRWETG(:)+ZRDRYG(:)) 
!
   P_RC_WETG(:) = - ZZW1(:)
   P_CC_WETG(:) = P_RC_WETG(:) * PCCT(:)/MAX(PRCT(:),LIMAP%XRTMIN(2))
   P_RR_WETG(:) = - ZZW7(:) + ZZW1(:)
   P_CR_WETG(:) = P_RR_WETG(:) * PCRT(:)/MAX(PRRT(:),LIMAP%XRTMIN(3))
   P_RI_WETG(:) = - ZZW5(:)
   P_CI_WETG(:) = P_RI_WETG(:) * PCIT(:)/MAX(PRIT(:),LIMAP%XRTMIN(4))
   P_RS_WETG(:) = - ZZW6(:)
   P_CS_WETG(:) = - ZZW6N(:)
   P_RG_WETG(:) = - PRGT(:)/PTSTEP * ZZW(:) + ZRWETG(:) * (1.-ZZW(:))
   P_CG_WETG(:) = - PCGT(:)/PTSTEP * ZZW(:)
   P_RH_WETG(:) =   PRGT(:)/PTSTEP * ZZW(:) + ZRWETG(:) *     ZZW(:)
   !
   P_TH_WETG(:) = ZZW7(:) * (PLSFACT(:)-PLVFACT(:))
END WHERE
!
!            1.g Dry mode
!            ------------
!
WHERE( ODCOMPUTE(:) .AND. PRGT(:)>LIMAP%XRTMIN(6) .AND. PCGT(:)>LIMAP%XCTMIN(6) .AND. PT(:)<CST%XTT .AND.                  &
       (ZRDRYG(:)-ZZW2(:)-ZZW3(:))<(ZRWETG(:)-ZZW5(:)-ZZW6(:)) .AND. ZRDRYG(:)>0.0 )
   !
   P_RC_DRYG(:) = - ZZW1(:)
   P_CC_DRYG(:) = P_RC_DRYG(:) * PCCT(:)/MAX(PRCT(:),LIMAP%XRTMIN(2))
   P_RR_DRYG(:) = - ZZW4(:)
   P_CR_DRYG(:) = - ZZW4N(:)
   P_RI_DRYG(:) = - ZZW2(:)
   P_CI_DRYG(:) = P_RI_DRYG(:) * PCIT(:)/MAX(PRIT(:),LIMAP%XRTMIN(4))
   P_RS_DRYG(:) = - ZZW3(:)
   P_CS_DRYG(:) = - ZZW3N(:)
   P_RG_DRYG(:) =   ZRDRYG(:)
   !
   P_TH_DRYG(:) = (ZZW1(:) + ZZW4(:)) * (PLSFACT(:)-PLVFACT(:))
END WHERE
!
!            1.h Prevent graupel growth for too small graupel diameter
!            --------------------------------------------------------
!
IF (LIMAP%LSIGMOIDE_G) THEN
!
   ZSIGMOIDE(:) = 1/(1 + exp(-LIMAP%XSIGMOIDE_G*(PRGT(:)-LIMAM%XMINDG/PRHODREF(:))))
!
   P_RC_DRYG(:) = P_RC_DRYG(:) * ZSIGMOIDE(:)
   P_CC_DRYG(:) = P_CC_DRYG(:) * ZSIGMOIDE(:)
   P_RI_DRYG(:) = P_RI_DRYG(:) * ZSIGMOIDE(:)
   P_CI_DRYG(:) = P_CI_DRYG(:) * ZSIGMOIDE(:)
   P_RR_DRYG(:) = P_RR_DRYG(:) * ZSIGMOIDE(:)
   P_CR_DRYG(:) = P_CR_DRYG(:) * ZSIGMOIDE(:)
   P_RS_DRYG(:) = P_RS_DRYG(:) * ZSIGMOIDE(:)
   P_CS_DRYG(:) = P_CS_DRYG(:) * ZSIGMOIDE(:)
   P_RG_DRYG(:) = P_RG_DRYG(:) * ZSIGMOIDE(:)
   P_TH_DRYG(:) = P_TH_DRYG(:) * ZSIGMOIDE(:)
!
   P_RC_WETG(:) = P_RC_WETG(:) * ZSIGMOIDE(:)
   P_CC_WETG(:) = P_CC_WETG(:) * ZSIGMOIDE(:)
   P_RI_WETG(:) = P_RI_WETG(:) * ZSIGMOIDE(:)
   P_CI_WETG(:) = P_CI_WETG(:) * ZSIGMOIDE(:)
   P_RR_WETG(:) = P_RR_WETG(:) * ZSIGMOIDE(:)
   P_CR_WETG(:) = P_CR_WETG(:) * ZSIGMOIDE(:)
   P_RS_WETG(:) = P_RS_WETG(:) * ZSIGMOIDE(:)
   P_CS_WETG(:) = P_CS_WETG(:) * ZSIGMOIDE(:)
   P_RG_WETG(:) = P_RG_WETG(:) * ZSIGMOIDE(:)
   P_TH_WETG(:) = P_TH_WETG(:) * ZSIGMOIDE(:)
!
END IF 
!
!
!*       2.  Hallett-Mossop process (HMG)
!        --------------------------------
!
! BVIE test ZRDRYG<ZZW ?????????????????????????
!GDRY(:) = (PT(:)<LIMAM%XHMTMAX) .AND. (PT(:)>LIMAM%XHMTMIN)    .AND. (ZRDRYG(:)<ZZW(:))&
GDRY(:) = PT(:)<LIMAM%XHMTMAX .AND. PT(:)>LIMAM%XHMTMIN .AND. PRGT(:)>LIMAP%XRTMIN(6) &
     .AND. PRCT(:)>LIMAP%XRTMIN(2) .AND. ODCOMPUTE(:) .AND. &
     PCGT(:)>LIMAP%XCTMIN(6) .AND. PCCT(:)>LIMAP%XCTMIN(2) &
     .AND. (ZRDRYG(:)-ZZW2(:)-ZZW3(:))<(ZRWETG(:)-ZZW5(:)-ZZW6(:))

ZZX(:)=9999.
ZVEC1(:)=0.
ZVEC2(:)=0.
IVEC1(:)=0
IVEC2(:)=0
WHERE( GDRY(:) )
   ZVEC1(:) = PLBDC(:)
   ZVEC2(:) = MAX( 1.0001, MIN( REAL(LIMAM%NGAMINC)-0.0001,           &
                         LIMAM%XHMLINTP1 * LOG( ZVEC1(:) ) + LIMAM%XHMLINTP2 ) )
   IVEC2(:) = INT( ZVEC2(:) )
   ZVEC2(:) = ZVEC2(:) - REAL( IVEC2(:) )
   ZVEC1(:) =   LIMAM%XGAMINC_HMC( IVEC2(:)+1 )* ZVEC2(:)      &
                    - LIMAM%XGAMINC_HMC( IVEC2(:)   )*(ZVEC2(:) - 1.0)
   ZZX(:) = ZVEC1(:) ! Large droplets
END WHERE
WHERE ( GDRY(:) .AND. ZZX(:)<0.99 ) ! Dry case
   P_CI_HMG(:) = ZZW1(:)*(PCCT(:)/PRCT(:))*(1.0-ZZX(:))*LIMAM%XHM_FACTG*  &
        MAX( 0.0, MIN( (PT(:)-LIMAM%XHMTMIN)/3.0,(LIMAM%XHMTMAX-PT(:))/2.0 ) )
   P_RI_HMG(:) = P_CI_HMG(:) * LIMAC%XMNU0
   P_RG_HMG(:) = - P_RI_HMG(:)
END WHERE
!
!            2.bis  Prevent Hallett-Mossop mechanism if graupel too small (no collection)
!            ----------------------------------------------------------------------------
!
IF (LIMAP%LSIGMOIDE_G) THEN
   P_CI_HMG(:) = P_CI_HMG(:) * ZSIGMOIDE(:)
   P_RI_HMG(:) = P_RI_HMG(:) * ZSIGMOIDE(:)
   P_RG_HMG(:) = P_RG_HMG(:) * ZSIGMOIDE(:)
END IF 
!
!
!*       3.  Graupel Melting
!        -------------------
!
ZZX(:) = 0.0
WHERE( PRGT(:)>LIMAP%XRTMIN(6) .AND. PCGT(:)>LIMAP%XCTMIN(6) .AND. PT(:)>CST%XTT .AND. ODCOMPUTE(:) )
   ZZX(:) = PRVT(:)*PPRES(:)/((CST%XMV/CST%XMD)+PRVT(:)) ! Vapor pressure
   ZZX(:) = PKA(:)*(CST%XTT-PT(:)) +                                 &
              ( PDV(:)*(CST%XLVTT + ( CST%XCPV - CST%XCL ) * ( PT(:) - CST%XTT )) &
                          *(CST%XESTT-ZZX(:))/(CST%XRV*PT(:))             )
!
! compute RGMLTR
!
   ZZX(:)  = MAX( 0.0,( -ZZX(:) * PCGT(:) *                        &
                          ( LIMAM%X0DEPG*       PLBDG(:)**LIMAM%XEX0DEPG +     &
                            LIMAM%X1DEPG*PCJ(:)*PLBDG(:)**LIMAM%XEX1DEPG ) -   &
                        ( ZZW1(:)+ZZW4(:) ) * ( CST%XCL*(CST%XTT-PT(:))) ) &
                        / CST%XLMTT                                    )
   P_RR_GMLT(:) = ZZX(:)
   P_CR_GMLT(:) = ZZX(:) * 5.0E6  ! obtained after averaging, Dshed=1mm and 500 microns 
   P_CG_GMLT(:) = - ZZX(:) * PCGT(:)/PRGT(:)
   !
   P_TH_GMLT(:) = - P_RR_GMLT(:) * (PLSFACT(:)-PLVFACT(:))
END WHERE
!
!
!
!
PA_RC(:) = PA_RC(:) + P_RC_WETG(:) + P_RC_DRYG(:) 
PA_CC(:) = PA_CC(:) + P_CC_WETG(:) + P_CC_DRYG(:) 
PA_RR(:) = PA_RR(:) + P_RR_WETG(:) + P_RR_DRYG(:)               + P_RR_GMLT(:)
PA_CR(:) = PA_CR(:) + P_CR_WETG(:) + P_CR_DRYG(:)               + P_CR_GMLT(:)
PA_RI(:) = PA_RI(:) + P_RI_WETG(:) + P_RI_DRYG(:) + P_RI_HMG(:) 
PA_CI(:) = PA_CI(:) + P_CI_WETG(:) + P_CI_DRYG(:) + P_CI_HMG(:) 
PA_RS(:) = PA_RS(:) + P_RS_WETG(:) + P_RS_DRYG(:) 
PA_CS(:) = PA_CS(:) + P_CS_WETG(:) + P_CS_DRYG(:) 
PA_RG(:) = PA_RG(:) + P_RG_WETG(:) + P_RG_DRYG(:) + P_RG_HMG(:) - P_RR_GMLT(:)
PA_CG(:) = PA_CG(:) + P_CG_WETG(:)                              + P_CG_GMLT(:)
PA_RH(:) = PA_RH(:) + P_RH_WETG(:) 
PA_CH(:) = PA_CH(:) - P_CG_WETG(:) 
PA_TH(:) = PA_TH(:) + P_TH_WETG(:) + P_TH_DRYG(:)               + P_TH_GMLT(:) 
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('LIMA_GRAUPEL', 1, ZHOOK_HANDLE)
CONTAINS
  FUNCTION GET_XKER_SDRYG(KSIZE,KG,KS) RESULT(RET)
    USE YOMHOOK, ONLY:LHOOK, DR_HOOK, JPHOOK
    INTEGER, INTENT(IN) :: KSIZE
    INTEGER, DIMENSION(KSIZE), INTENT(IN) :: KG
    INTEGER, DIMENSION(KSIZE), INTENT(IN) :: KS
    REAL, DIMENSION(KSIZE) :: RET
    !
    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
    INTEGER I
    !
    IF (LHOOK) CALL DR_HOOK('GET_XKER_SDRYG', 0, ZHOOK_HANDLE)
    DO I=1,SIZE(KG)
       RET(I) = LIMAM%XKER_SDRYG(MAX(MIN(KG(I),SIZE(LIMAM%XKER_SDRYG,1)),1),&
            MAX(MIN(KS(I),SIZE(LIMAM%XKER_SDRYG,2)),1))
    END DO
    IF (LHOOK) CALL DR_HOOK('GET_XKER_SDRYG', 1, ZHOOK_HANDLE)
  END FUNCTION GET_XKER_SDRYG
!
!-------------------------------------------------------------------------------
!
  FUNCTION GET_XKER_N_SDRYG(KSIZE,KG,KS) RESULT(RET)
  USE YOMHOOK, ONLY:LHOOK, DR_HOOK, JPHOOK
    INTEGER, INTENT(IN) :: KSIZE
    INTEGER, DIMENSION(KSIZE), INTENT(IN) :: KG
    INTEGER, DIMENSION(KSIZE), INTENT(IN) :: KS
    REAL, DIMENSION(KSIZE) :: RET
    !
    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
    INTEGER I
    !
    IF (LHOOK) CALL DR_HOOK('GET_XKER_N_SDRYG', 0, ZHOOK_HANDLE)
    DO I=1,SIZE(KG)
       RET(I) = LIMAM%XKER_N_SDRYG(MAX(MIN(KG(I),SIZE(LIMAM%XKER_N_SDRYG,1)),1),&
            MAX(MIN(KS(I),SIZE(LIMAM%XKER_N_SDRYG,2)),1))
    END DO
    IF (LHOOK) CALL DR_HOOK('GET_XKER_N_SDRYG', 1, ZHOOK_HANDLE)
  END FUNCTION GET_XKER_N_SDRYG
!
!-------------------------------------------------------------------------------
!
  FUNCTION GET_XKER_RDRYG(KSIZE,KG,KR) RESULT(RET)
    USE YOMHOOK, ONLY:LHOOK, DR_HOOK, JPHOOK
    INTEGER, INTENT(IN) :: KSIZE
    INTEGER, DIMENSION(KSIZE), INTENT(IN) :: KG
    INTEGER, DIMENSION(KSIZE), INTENT(IN) :: KR
    REAL, DIMENSION(KSIZE) :: RET
    !
    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
    INTEGER I
    !
    IF (LHOOK) CALL DR_HOOK('GET_XKER_RDRYG', 0, ZHOOK_HANDLE)
    DO I=1,SIZE(KG)
       RET(I) = LIMAM%XKER_RDRYG(MAX(MIN(KG(I),SIZE(LIMAM%XKER_RDRYG,1)),1),&
            MAX(MIN(KR(I),SIZE(LIMAM%XKER_RDRYG,2)),1))
    END DO
  IF (LHOOK) CALL DR_HOOK('GET_XKER_RDRYG', 1, ZHOOK_HANDLE)
  END FUNCTION GET_XKER_RDRYG
!
!-------------------------------------------------------------------------------
!
  FUNCTION GET_XKER_N_RDRYG(KSIZE,KG,KR) RESULT(RET)
    USE YOMHOOK, ONLY:LHOOK, DR_HOOK, JPHOOK
    INTEGER, INTENT(IN) :: KSIZE
    INTEGER, DIMENSION(KSIZE), INTENT(IN) :: KG
    INTEGER, DIMENSION(KSIZE), INTENT(IN) :: KR
    REAL, DIMENSION(KSIZE) :: RET
    !
    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
    INTEGER I
    !
    IF (LHOOK) CALL DR_HOOK('GET_XKER_N_RDRYG', 0, ZHOOK_HANDLE)
    DO I=1,SIZE(KG)
       RET(I) = LIMAM%XKER_N_RDRYG(MAX(MIN(KG(I),SIZE(LIMAM%XKER_N_RDRYG,1)),1),&
            MAX(MIN(KR(I),SIZE(LIMAM%XKER_N_RDRYG,2)),1))
    END DO
  IF (LHOOK) CALL DR_HOOK('GET_XKER_N_RDRYG', 1, ZHOOK_HANDLE)
END FUNCTION GET_XKER_N_RDRYG
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_GRAUPEL
END MODULE MODE_LIMA_GRAUPEL
