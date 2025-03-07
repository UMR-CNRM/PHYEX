!MNH_LIC Copyright 2018-2025 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_LIMA_HAIL
  IMPLICIT NONE
CONTAINS
!     #################################################################################
  SUBROUTINE LIMA_HAIL (CST, LIMAP, LIMAM, KSIZE, PTSTEP, ODCOMPUTE,                              &
                        PRHODREF, PPRES, PT, PKA, PDV, PCJ,                    &
                        PRVT, PRCT, PRRT, PRIT, PRST, PRGT, PRHT,              &
                        PCCT, PCRT, PCIT, PCST, PCGT, PCHT,                    &
                        PLBDC, PLBDR, PLBDS, PLBDG, PLBDH,                     &
                        PLVFACT, PLSFACT,                                      &
                        P_TH_WETH, P_RC_WETH, P_CC_WETH, P_RR_WETH, P_CR_WETH, &
                        P_RI_WETH, P_CI_WETH, P_RS_WETH, P_CS_WETH, P_RG_WETH, P_CG_WETH, P_RH_WETH, &
                        P_RG_COHG, P_CG_COHG,                                  &
                        P_TH_HMLT, P_RR_HMLT, P_CR_HMLT, P_CH_HMLT,            &
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
USE MODD_PARAM_LIMA, ONLY:PARAM_LIMA_T
USE MODD_CST, ONLY:CST_T
USE YOMHOOK, ONLY:LHOOK, DR_HOOK, JPHOOK

!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(PARAM_LIMA_MIXED_T),INTENT(IN)::LIMAM
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
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRHT    !
!
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PCCT    !
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PCRT    !
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PCIT    !
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PCST    !
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PCGT    !
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PCHT    !
!
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLBDC   ! 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLBDR   ! 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLBDS   ! 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLBDG   ! 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLBDH   ! 
!
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLVFACT ! 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLSFACT ! 
!
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_TH_WETH
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_RC_WETH
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_CC_WETH
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_RR_WETH
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_CR_WETH
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_RI_WETH
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_CI_WETH
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_RS_WETH
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_CS_WETH
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_RG_WETH
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_CG_WETH
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_RH_WETH
!
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_RG_COHG
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_CG_COHG
!
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_TH_HMLT
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_RR_HMLT
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_CR_HMLT
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_CH_HMLT
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
LOGICAL, DIMENSION(SIZE(PRCT))  :: GWET
!
REAL,    DIMENSION(SIZE(PRCT))  :: Z1, Z2, Z3, Z4
REAL,    DIMENSION(SIZE(PRCT))  :: ZZW, ZZW1, ZZW2, ZZW3, ZZW4, ZZW5, ZZW6
REAL,    DIMENSION(SIZE(PRCT))  :: ZZW3N, ZZW4N
REAL,    DIMENSION(SIZE(PRCT))  :: ZRWETH
!
INTEGER, DIMENSION(SIZE(PRCT))  :: IVEC1,IVEC2        ! Vectors of indices
REAL,    DIMENSION(SIZE(PRCT))  :: ZVEC1,ZVEC2, ZVEC3 ! Work vectors
!
REAL                            :: ZTHRH, ZTHRC
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
!
IF (LHOOK) CALL DR_HOOK('LIMA_HAIL', 0, ZHOOK_HANDLE)
P_TH_WETH(:) = 0.
P_RC_WETH(:) = 0.
P_CC_WETH(:) = 0.
P_RR_WETH(:) = 0.
P_CR_WETH(:) = 0.
P_RI_WETH(:) = 0.
P_CI_WETH(:) = 0.
P_RS_WETH(:) = 0.
P_CS_WETH(:) = 0.
P_RG_WETH(:) = 0.
P_CG_WETH(:) = 0.
P_RH_WETH(:) = 0.
!
P_RG_COHG(:) = 0.
P_CG_COHG(:) = 0.
!
P_TH_HMLT(:) = 0.
P_RR_HMLT(:) = 0.
P_CR_HMLT(:) = 0.
P_CH_HMLT(:) = 0.
!
ZZW1(:) = 0. ! RCWETH
ZZW2(:) = 0. ! RIWETH
ZZW3(:) = 0. ! RSWETH
ZZW3N(:) = 0.! NSWETH
ZZW4(:) = 0. ! RGWETH
ZZW4N(:) = 0.! NGWETH
ZZW5(:) = 0. ! RCWETH+RRWETH
ZZW6(:) = 0. ! RSWETH
!
ZRWETH(:) = 0.
!
!
!*       1.  Hail growth by collection (wet case only ?)
!        -----------------------------------------------
!
!            1.a Collection of rc and ri
!            ---------------------------
!
WHERE( PRHT(:)>LIMAP%XRTMIN(7) .AND. PCHT(:)>LIMAP%XCTMIN(7) .AND. ODCOMPUTE(:) )
   ZZW(:) = PCHT(:) * PLBDH(:)**(-LIMAM%XDH-2.0) * PRHODREF(:)**(1-LIMAP%XCEXVT)
   ZZW1(:) = LIMAM%XFWETH * PRCT(:) * ZZW(:) ! RCWETH
   ZZW2(:) = LIMAM%XFWETH * PRIT(:) * ZZW(:) ! RIWETH
END WHERE
!
!*           1.b Collection of rs
!            --------------------
!
GWET(:) = PRST(:)>LIMAP%XRTMIN(5) .AND. PRHT(:)>LIMAP%XRTMIN(7) .AND. ODCOMPUTE(:) .AND. &
          PCST(:)>LIMAP%XCTMIN(5) .AND. PCHT(:)>LIMAP%XCTMIN(7)
!
WHERE( GWET )
!
!*       Select the (ZLBDAG,ZLBDAS) couplet
!
   ZVEC1(:) = PLBDH(:)
   ZVEC2(:) = PLBDS(:)
!
!*       find the next lower indice for the ZLBDAG and for the ZLBDAS
!        in the geometrical set of (Lbda_g,Lbda_s) couplet use to
!        tabulate the SDRYG-kernel
!
   ZVEC1(:) = MAX( 1.0001, MIN( REAL(LIMAM%NWETLBDAH)-0.0001,           &
                         LIMAM%XWETINTP1H * LOG( ZVEC1(:) ) + LIMAM%XWETINTP2H ) )
   IVEC1(:) = INT( ZVEC1(:) )
   ZVEC1(:) = ZVEC1(:) - REAL( IVEC1(:) )
!
   ZVEC2(:) = MAX( 1.0001, MIN( REAL(LIMAM%NWETLBDAS)-0.0001,           &
                         LIMAM%XWETINTP1S * LOG( ZVEC2(:) ) + LIMAM%XWETINTP2S ) )
   IVEC2(:) = INT( ZVEC2(:) )
   ZVEC2(:) = ZVEC2(:) - REAL( IVEC2(:) )
!
!*       perform the bilinear interpolation of the normalized
!        SDRYG-kernel
   !
   Z1(:) = GET_XKER_SWETH(KSIZE,IVEC1(:)+1,IVEC2(:)+1)
   Z2(:) = GET_XKER_SWETH(KSIZE,IVEC1(:)+1,IVEC2(:)  )
   Z3(:) = GET_XKER_SWETH(KSIZE,IVEC1(:)  ,IVEC2(:)+1)
   Z4(:) = GET_XKER_SWETH(KSIZE,IVEC1(:)  ,IVEC2(:)  )
   ZVEC3(:) =  (      Z1(:)* ZVEC2(:)          &
                    - Z2(:)*(ZVEC2(:) - 1.0) ) &
                                               *  ZVEC1(:)    &
                 - (  Z3(:)* ZVEC2(:)          &
                    - Z4(:)*(ZVEC2(:) - 1.0) ) &
                                                    * (ZVEC1(:) - 1.0)
   ZZW(:) = ZVEC3(:)
   ZZW3(:) = LIMAM%XFSWETH * ZZW(:)                                     & ! RSWETH
                    *  PRST(:) * PCHT(:)                          &
                    *  PRHODREF(:)**(1-LIMAP%XCEXVT)                    &
                    *( LIMAM%XLBSWETH1/( PLBDH(:)**2                ) + &
                       LIMAM%XLBSWETH2/( PLBDH(:)    * PLBDS(:)     ) + &
                       LIMAM%XLBSWETH3/(               PLBDS(:)**2) )
!
   Z1(:) = GET_XKER_N_SWETH(KSIZE,IVEC1(:)+1,IVEC2(:)+1)
   Z2(:) = GET_XKER_N_SWETH(KSIZE,IVEC1(:)+1,IVEC2(:)  )
   Z3(:) = GET_XKER_N_SWETH(KSIZE,IVEC1(:)  ,IVEC2(:)+1)
   Z4(:) = GET_XKER_N_SWETH(KSIZE,IVEC1(:)  ,IVEC2(:)  )
   ZVEC3(:) =  (      Z1(:)* ZVEC2(:)          &
                    - Z2(:)*(ZVEC2(:) - 1.0) ) &
                                               *  ZVEC1(:)    &
                 - (  Z3(:)* ZVEC2(:)          &
                    - Z4(:)*(ZVEC2(:) - 1.0) ) &
                                                    * (ZVEC1(:) - 1.0)
   ZZW(:) = ZVEC3(:)
   ZZW3N(:) = LIMAM%XFNSWETH * ZZW(:)                                      & ! NSWETH
                      *  PCST(:) * PCHT(:)                           &
                      *  PRHODREF(:)**(1-LIMAP%XCEXVT)                     &
                      *( LIMAM%XLBNSWETH1/( PLBDH(:)**2                ) + &
                         LIMAM%XLBNSWETH2/( PLBDH(:)    * PLBDS(:)     ) + &
                         LIMAM%XLBNSWETH3/(               PLBDS(:)**2) )
END WHERE
!
!*           1.c  Collection of rg
!            ---------------------
!
GWET(:) = PRGT(:)>LIMAP%XRTMIN(6) .AND. PRHT(:)>LIMAP%XRTMIN(7) .AND. ODCOMPUTE(:) .AND. &
          PCGT(:)>LIMAP%XCTMIN(6) .AND. PCHT(:)>LIMAP%XCTMIN(7)
!
WHERE( GWET )
!
!*       Select the (ZLBDAG,ZLBDAR) couplet
!
   ZVEC1(:) = PLBDH(:)
   ZVEC2(:) = PLBDG(:)
!
!*       Find the next lower indice for the ZLBDAG and for the ZLBDAR
!        in the geometrical set of (Lbda_g,Lbda_r) couplet use to
!        tabulate the RDRYG-kernel
!
   ZVEC1(:) = MAX( 1.0001, MIN( REAL(LIMAM%NWETLBDAH)-0.0001,           &
                         LIMAM%XWETINTP1H * LOG( ZVEC1(:) ) + LIMAM%XWETINTP2H ) )
   IVEC1(:) = INT( ZVEC1(:) )
   ZVEC1(:) = ZVEC1(:) - REAL( IVEC1(:) )
!
   ZVEC2(:) = MAX( 1.0001, MIN( REAL(LIMAM%NWETLBDAG)-0.0001,           &
                         LIMAM%XWETINTP1G * LOG( ZVEC2(:) ) + LIMAM%XWETINTP2G ) )
   IVEC2(:) = INT( ZVEC2(:) )
   ZVEC2(:) = ZVEC2(:) - REAL( IVEC2(:) )
!
!*       Perform the bilinear interpolation of the normalized
!        RDRYG-kernel
!
   Z1(:) = GET_XKER_GWETH(KSIZE,IVEC1(:)+1,IVEC2(:)+1)
   Z2(:) = GET_XKER_GWETH(KSIZE,IVEC1(:)+1,IVEC2(:)  )
   Z3(:) = GET_XKER_GWETH(KSIZE,IVEC1(:)  ,IVEC2(:)+1)
   Z4(:) = GET_XKER_GWETH(KSIZE,IVEC1(:)  ,IVEC2(:)  )
      ZVEC3(:) =  (   Z1(:)* ZVEC2(:)          &
                    - Z2(:)*(ZVEC2(:) - 1.0) ) &
                                               *  ZVEC1(:)   &
                 - (  Z3(:)* ZVEC2(:)          &
                    - Z4(:)*(ZVEC2(:) - 1.0) ) &
                                               * (ZVEC1(:) - 1.0)
   ZZW(:) = ZVEC3(:)
   ZZW4(:) = LIMAM%XFGWETH * ZZW(:)                                     & ! RGWETH
                     * PRGT(:) * PCHT(:)                          &
                     * PRHODREF(:)**(1-LIMAP%XCEXVT)                    &
                     *( LIMAM%XLBGWETH1/( PLBDH(:)**2               ) + &
                        LIMAM%XLBGWETH2/( PLBDH(:)   * PLBDG(:)     ) + &
                        LIMAM%XLBGWETH3/(              PLBDG(:)**2) )
!
   Z1(:) = GET_XKER_N_GWETH(KSIZE,IVEC1(:)+1,IVEC2(:)+1)
   Z2(:) = GET_XKER_N_GWETH(KSIZE,IVEC1(:)+1,IVEC2(:)  )
   Z3(:) = GET_XKER_N_GWETH(KSIZE,IVEC1(:)  ,IVEC2(:)+1)
   Z4(:) = GET_XKER_N_GWETH(KSIZE,IVEC1(:)  ,IVEC2(:)  )
      ZVEC3(:) =  (   Z1(:)* ZVEC2(:)          &
                    - Z2(:)*(ZVEC2(:) - 1.0) ) &
                                               *  ZVEC1(:)   &
                 - (  Z3(:)* ZVEC2(:)          &
                    - Z4(:)*(ZVEC2(:) - 1.0) ) &
                                               * (ZVEC1(:) - 1.0)
   ZZW(:) = ZVEC3(:)
   ZZW4N(:) = LIMAM%XFNGWETH * ZZW(:)                                      & ! NGWETH
                       * PCGT(:) * PCHT(:)                           &
                       * PRHODREF(:)**(1-LIMAP%XCEXVT)                     &
                       *( LIMAM%XLBNGWETH1/( PLBDH(:)**2               ) + &
                          LIMAM%XLBNGWETH2/( PLBDH(:)   * PLBDG(:)     ) + &
                          LIMAM%XLBNGWETH3/(              PLBDG(:)**2) )
   
END WHERE
!
!            1.e Wet growth of hail
!            ----------------------
!
ZZW(:) = 0.0
WHERE( PRHT(:)>LIMAP%XRTMIN(6) .AND. PCHT(:) >LIMAP%XCTMIN(6) .AND. ODCOMPUTE(:) )
   ZZW(:) = PRVT(:)*PPRES(:)/((CST%XMV/CST%XMD)+PRVT(:)) ! Vapor pressure
   ZZW(:) = PKA(:)*(CST%XTT-PT(:)) +                                  &
              ( PDV(:)*(CST%XLVTT + ( CST%XCPV - CST%XCL ) * ( PT(:) - CST%XTT ))   &
                          *(CST%XESTT-ZZW(:))/(CST%XRV*PT(:))             )
!
! Total mass gained by hail in wet mode
   ZRWETH(:)  = MAX( 0.0,                                                                &
                   ( ZZW(:) * PCHT(:) * ( LIMAM%X0DEPH*       PLBDH(:)**LIMAM%XEX0DEPH +             &
                                          LIMAM%X1DEPH*PCJ(:)*PLBDH(:)**LIMAM%XEX1DEPH ) +           &
                   ( ZZW2(:)+ZZW3(:)+ZZW4(:) ) * ( CST%XLMTT + (CST%XCI-CST%XCL)*(CST%XTT-PT(:)) )   )   &
                   / (CST%XLMTT-CST%XCL*(CST%XTT-PT(:)))                                     )
   ! We must agregate, at least, the cold species
   ZRWETH(:)=MAX(ZRWETH(:), ZZW2(:)+ZZW3(:)+ZZW4(:))
END WHERE
WHERE( PRHT(:)>LIMAP%XRTMIN(6) .AND. PCHT(:) >LIMAP%XCTMIN(6) .AND. PRRT(:)>LIMAP%XRTMIN(3) &
     .AND. PCRT(:) >LIMAP%XCTMIN(3) .AND. ODCOMPUTE(:) )
   ! Mass of rain frozen by hail RRWETH
   ZZW5(:) = ZRWETH(:) - ZZW2(:) - ZZW3(:) - ZZW4(:) - ZZW1(:)
END WHERE
!
ZZW(:) = 0.0
WHERE( ODCOMPUTE(:) .AND. PT(:)<CST%XTT .AND. ZZW5(:)>0.0 ) 
   P_RC_WETH(:) = - ZZW1(:)
   P_CC_WETH(:) = P_RC_WETH(:) * PCCT(:)/MAX(PRCT(:),LIMAP%XRTMIN(2))
   P_RR_WETH(:) = - ZZW5(:)
   P_CR_WETH(:) = P_RR_WETH(:) * PCRT(:)/MAX(PRRT(:),LIMAP%XRTMIN(3))
   P_RI_WETH(:) = - ZZW2(:)
   P_CI_WETH(:) = P_RI_WETH(:) * PCIT(:)/MAX(PRIT(:),LIMAP%XRTMIN(4))
   P_RS_WETH(:) = - ZZW3(:)
   P_CS_WETH(:) = - ZZW3N(:)
   P_RG_WETH(:) = - ZZW4(:)
   P_CG_WETH(:) = - ZZW4N(:)
   P_RH_WETH(:) =   ZRWETH(:)
   !
   P_TH_WETH(:) = (ZZW5(:)+ZZW1(:)) * (PLSFACT(:)-PLVFACT(:))
END WHERE
!
!
!*       2.  Hail -> graupel conversion (COHG)
!        -------------------------------------
!
ZTHRH=0.01E-3
ZTHRC=0.001E-3
ZZW(:) = 0.0
GWET(:) = PRHT(:)<ZTHRH .AND. PRCT(:)<ZTHRC .AND. PT(:)<CST%XTT 
WHERE( GWET(:) )
   P_RG_COHG(:) = PRHT(:) * MIN( 1.0,MAX( 0.0,1.0-(PRCT(:)/ZTHRC) ) )
   P_CG_COHG(:) = P_RG_COHG(:) * PCHT(:)/MAX(PRHT(:),LIMAP%XRTMIN(7))
END WHERE
!
!
!*       3.  Hail Melting (HMLT)
!        -----------------------
!
ZZW(:) = 0.0
WHERE( PRHT(:)>LIMAP%XRTMIN(6) .AND. PCHT(:)>LIMAP%XCTMIN(7) .AND. PT(:)>CST%XTT .AND. ODCOMPUTE(:) )
   ZZW(:) = PRVT(:)*PPRES(:)/((CST%XMV/CST%XMD)+PRVT(:)) ! Vapor pressure
   ZZW(:) = PKA(:)*(CST%XTT-PT(:)) +                                 &
              ( PDV(:)*(CST%XLVTT + ( CST%XCPV - CST%XCL ) * ( PT(:) - CST%XTT )) &
                          *(CST%XESTT-ZZW(:))/(CST%XRV*PT(:))             )
!
! compute RHMLTR
!
   ZZW(:)  = MAX( 0.0,( -ZZW(:) * PCHT(:) *                        &
                          ( LIMAM%X0DEPH*       PLBDH(:)**LIMAM%XEX0DEPH +     &
                            LIMAM%X1DEPH*PCJ(:)*PLBDH(:)**LIMAM%XEX1DEPH ) -   &
                        ( ZZW5(:) ) * ( CST%XCL*(CST%XTT-PT(:))) ) &
                        / CST%XLMTT                                    )
   P_RR_HMLT(:) = ZZW(:)
   P_CR_HMLT(:) = ZZW(:) * 5.0E6  ! obtained after averaging, Dshed=1mm and 500 microns 
   P_CH_HMLT(:) = - ZZW(:) * PCHT(:)/PRHT(:)
   !
   P_TH_HMLT(:) = - P_RR_HMLT(:) * (PLSFACT(:)-PLVFACT(:))
END WHERE
!
!
!
!
PA_RC(:) = PA_RC(:) + P_RC_WETH(:)
PA_CC(:) = PA_CC(:) + P_CC_WETH(:)
PA_RR(:) = PA_RR(:) + P_RR_WETH(:)                + P_RR_HMLT(:)
PA_CR(:) = PA_CR(:) + P_CR_WETH(:)                + P_CR_HMLT(:)
PA_RI(:) = PA_RI(:) + P_RI_WETH(:) 
PA_CI(:) = PA_CI(:) + P_CI_WETH(:)
PA_RS(:) = PA_RS(:) + P_RS_WETH(:)
PA_CS(:) = PA_CS(:) + P_CS_WETH(:)
PA_RG(:) = PA_RG(:) + P_RG_WETH(:) + P_RG_COHG(:)
PA_CG(:) = PA_CG(:) + P_CG_WETH(:) + P_CG_COHG(:)
PA_RH(:) = PA_RH(:) + P_RH_WETH(:) - P_RG_COHG(:) - P_RR_HMLT(:)
PA_CH(:) = PA_CH(:)                - P_CG_COHG(:) + P_CH_HMLT(:)
PA_TH(:) = PA_TH(:) + P_TH_WETH(:)                + P_TH_HMLT(:) 
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('LIMA_HAIL', 1, ZHOOK_HANDLE)
CONTAINS
  FUNCTION GET_XKER_SWETH(KSIZE,KH,KS) RESULT(RET)
USE YOMHOOK, ONLY:LHOOK, DR_HOOK, JPHOOK
    INTEGER, INTENT(IN) :: KSIZE
    INTEGER, DIMENSION(KSIZE), INTENT(IN) :: KH
    INTEGER, DIMENSION(KSIZE), INTENT(IN) :: KS
    REAL, DIMENSION(KSIZE) :: RET
    !
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
    INTEGER I
    !
    IF (LHOOK) CALL DR_HOOK('GET_XKER_SWETH', 0, ZHOOK_HANDLE)
    DO I=1,SIZE(KH)
       RET(I) = LIMAM%XKER_SWETH(MAX(MIN(KH(I),SIZE(LIMAM%XKER_SWETH,1)),1),&
            MAX(MIN(KS(I),SIZE(LIMAM%XKER_SWETH,2)),1))
    END DO
  IF (LHOOK) CALL DR_HOOK('GET_XKER_SWETH', 1, ZHOOK_HANDLE)
  END FUNCTION GET_XKER_SWETH
!
!-------------------------------------------------------------------------------
!
  FUNCTION GET_XKER_N_SWETH(KSIZE,KH,KS) RESULT(RET)
USE YOMHOOK, ONLY:LHOOK, DR_HOOK, JPHOOK
    INTEGER, INTENT(IN) :: KSIZE
    INTEGER, DIMENSION(KSIZE), INTENT(IN) :: KH
    INTEGER, DIMENSION(KSIZE), INTENT(IN) :: KS
    REAL, DIMENSION(KSIZE) :: RET
    !
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
    INTEGER I
    !
    IF (LHOOK) CALL DR_HOOK('GET_XKER_N_SWETH', 0, ZHOOK_HANDLE)
    DO I=1,SIZE(KH)
       RET(I) = LIMAM%XKER_N_SWETH(MAX(MIN(KH(I),SIZE(LIMAM%XKER_N_SWETH,1)),1),&
            MAX(MIN(KS(I),SIZE(LIMAM%XKER_N_SWETH,2)),1))
    END DO
  IF (LHOOK) CALL DR_HOOK('GET_XKER_N_SWETH', 1, ZHOOK_HANDLE)
  END FUNCTION GET_XKER_N_SWETH
!
!-------------------------------------------------------------------------------
!
  FUNCTION GET_XKER_GWETH(KSIZE,KH,KG) RESULT(RET)
USE YOMHOOK, ONLY:LHOOK, DR_HOOK, JPHOOK
    INTEGER, INTENT(IN) :: KSIZE
    INTEGER, DIMENSION(KSIZE), INTENT(IN) :: KH
    INTEGER, DIMENSION(KSIZE), INTENT(IN) :: KG
    REAL, DIMENSION(KSIZE) :: RET
    !
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
    INTEGER I
    !
    IF (LHOOK) CALL DR_HOOK('GET_XKER_GWETH', 0, ZHOOK_HANDLE)
    DO I=1,SIZE(KH)
       RET(I) = LIMAM%XKER_GWETH(MAX(MIN(KH(I),SIZE(LIMAM%XKER_GWETH,1)),1),&
            MAX(MIN(KG(I),SIZE(LIMAM%XKER_GWETH,2)),1))
    END DO
  IF (LHOOK) CALL DR_HOOK('GET_XKER_GWETH', 1, ZHOOK_HANDLE)
END FUNCTION GET_XKER_GWETH
!
!-------------------------------------------------------------------------------
!
  FUNCTION GET_XKER_N_GWETH(KSIZE,KH,KG) RESULT(RET)
USE YOMHOOK, ONLY:LHOOK, DR_HOOK, JPHOOK
    INTEGER, INTENT(IN) :: KSIZE
    INTEGER, DIMENSION(KSIZE), INTENT(IN) :: KH
    INTEGER, DIMENSION(KSIZE), INTENT(IN) :: KG
    REAL, DIMENSION(KSIZE) :: RET
    !
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
    INTEGER I
    !
    IF (LHOOK) CALL DR_HOOK('GET_XKER_N_GWETH', 0, ZHOOK_HANDLE)
    DO I=1,SIZE(KH)
       RET(I) = LIMAM%XKER_N_GWETH(MAX(MIN(KH(I),SIZE(LIMAM%XKER_N_GWETH,1)),1),&
            MAX(MIN(KG(I),SIZE(LIMAM%XKER_N_GWETH,2)),1))
    END DO
  IF (LHOOK) CALL DR_HOOK('GET_XKER_N_GWETH', 1, ZHOOK_HANDLE)
END FUNCTION GET_XKER_N_GWETH
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_HAIL
END MODULE MODE_LIMA_HAIL
