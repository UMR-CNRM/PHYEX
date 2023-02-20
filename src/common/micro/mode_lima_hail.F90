!MNH_LIC Copyright 2018-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_LIMA_HAIL
  IMPLICIT NONE
CONTAINS
!     #################################################################################
  SUBROUTINE LIMA_HAIL (PTSTEP, LDCOMPUTE,                                     &
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
USE MODD_CST,              ONLY : XTT, XMD, XMV, XRD, XRV, XLVTT, XLMTT, XESTT, XCL, XCI, XCPV
USE MODD_PARAM_LIMA,       ONLY : XRTMIN, XCTMIN, XCEXVT
USE MODD_PARAM_LIMA_MIXED, ONLY : NWETLBDAG, XWETINTP1G, XWETINTP2G, &
                                  NWETLBDAH, X0DEPH, X1DEPH, XDH, XEX0DEPH, XEX1DEPH, &
                                  XFWETH, XWETINTP1H, XWETINTP2H, &
                                  NWETLBDAS, XWETINTP1S, XWETINTP2S, &
                                  XKER_N_GWETH, XKER_GWETH, XKER_N_SWETH, XKER_SWETH, &
                                  XFSWETH, XLBSWETH1, XLBSWETH2, XLBSWETH3, &
                                  XFNSWETH, XLBNSWETH1, XLBNSWETH2, XLBNSWETH3, &
                                  XFGWETH, XLBGWETH1, XLBGWETH2, XLBGWETH3, &
                                  XFNGWETH, XLBNGWETH1, XLBNGWETH2, XLBNGWETH3

USE MODD_PARAM_LIMA_COLD,  ONLY : XMNU0, XCXS, XBS
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
REAL,                 INTENT(IN)    :: PTSTEP 
LOGICAL, DIMENSION(:),INTENT(IN)    :: LDCOMPUTE
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRHODREF    ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PPRES    ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PT   ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PKA   ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PDV   ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PCJ   ! 
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRVT    ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PRCT    ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PRRT    ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PRIT    ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PRST    ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PRGT    !
REAL, DIMENSION(:),   INTENT(IN)    :: PRHT    !
!
REAL, DIMENSION(:),   INTENT(IN)    :: PCCT    !
REAL, DIMENSION(:),   INTENT(IN)    :: PCRT    !
REAL, DIMENSION(:),   INTENT(IN)    :: PCIT    !
REAL, DIMENSION(:),   INTENT(IN)    :: PCST    !
REAL, DIMENSION(:),   INTENT(IN)    :: PCGT    !
REAL, DIMENSION(:),   INTENT(IN)    :: PCHT    !
!
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDC   ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDR   ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDS   ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDG   ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDH   ! 
!
REAL, DIMENSION(:),   INTENT(IN)    :: PLVFACT ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PLSFACT ! 
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_TH_WETH
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RC_WETH
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CC_WETH
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RR_WETH
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CR_WETH
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RI_WETH
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CI_WETH
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RS_WETH
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CS_WETH
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RG_WETH
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CG_WETH
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RH_WETH
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RG_COHG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CG_COHG
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_TH_HMLT
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RR_HMLT
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CR_HMLT
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CH_HMLT
!
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_TH
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_RC
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_CC
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_RR
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_CR
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_RI
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_CI
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_RS
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_CS
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_RG
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_CG
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_RH
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_CH
!
!*       0.2   Declarations of local variables :
!
LOGICAL, DIMENSION(SIZE(PRCT))  :: GWET
INTEGER                         :: JJ
!
REAL,    DIMENSION(SIZE(PRCT))  :: Z1, Z2, Z3, Z4
REAL,    DIMENSION(SIZE(PRCT))  :: ZZX, ZZW, ZZW1, ZZW2, ZZW3, ZZW4, ZZW5, ZZW6
REAL,    DIMENSION(SIZE(PRCT))  :: ZZW3N, ZZW4N, ZZW6N 
REAL,    DIMENSION(SIZE(PRCT))  :: ZRWETH
!
INTEGER, DIMENSION(SIZE(PRCT))  :: IVEC1,IVEC2        ! Vectors of indices
REAL,    DIMENSION(SIZE(PRCT))  :: ZVEC1,ZVEC2, ZVEC3 ! Work vectors
!
INTEGER                         :: NHAIL
REAL                            :: ZTHRH, ZTHRC
!
!-------------------------------------------------------------------------------
!
!
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
WHERE( PRHT(:)>XRTMIN(7) .AND. PCHT(:)>XCTMIN(7) .AND. LDCOMPUTE(:) )
   ZZW(:) = PCHT(:) * PLBDH(:)**(-XDH-2.0) * PRHODREF(:)**(1-XCEXVT)
   ZZW1(:) = XFWETH * PRCT(:) * ZZW(:) ! RCWETH
   ZZW2(:) = XFWETH * PRIT(:) * ZZW(:) ! RIWETH
END WHERE
!
!*           1.b Collection of rs
!            --------------------
!
GWET(:) = PRST(:)>XRTMIN(5) .AND. PRHT(:)>XRTMIN(7) .AND. LDCOMPUTE(:) .AND. &
          PCST(:)>XCTMIN(5) .AND. PCHT(:)>XCTMIN(7)
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
   ZVEC1(:) = MAX( 1.0001, MIN( REAL(NWETLBDAH)-0.0001,           &
                         XWETINTP1H * LOG( ZVEC1(:) ) + XWETINTP2H ) )
   IVEC1(:) = INT( ZVEC1(:) )
   ZVEC1(:) = ZVEC1(:) - REAL( IVEC1(:) )
!
   ZVEC2(:) = MAX( 1.0001, MIN( REAL(NWETLBDAS)-0.0001,           &
                         XWETINTP1S * LOG( ZVEC2(:) ) + XWETINTP2S ) )
   IVEC2(:) = INT( ZVEC2(:) )
   ZVEC2(:) = ZVEC2(:) - REAL( IVEC2(:) )
!
!*       perform the bilinear interpolation of the normalized
!        SDRYG-kernel
   !
   Z1(:) = GET_XKER_SWETH(IVEC1(:)+1,IVEC2(:)+1)
   Z2(:) = GET_XKER_SWETH(IVEC1(:)+1,IVEC2(:)  )
   Z3(:) = GET_XKER_SWETH(IVEC1(:)  ,IVEC2(:)+1)
   Z4(:) = GET_XKER_SWETH(IVEC1(:)  ,IVEC2(:)  )
   ZVEC3(:) =  (      Z1(:)* ZVEC2(:)          &
                    - Z2(:)*(ZVEC2(:) - 1.0) ) &
      			 	                            *  ZVEC1(:)    &
                 - (  Z3(:)* ZVEC2(:)          &
                    - Z4(:)*(ZVEC2(:) - 1.0) ) &
       			                                    * (ZVEC1(:) - 1.0)
   ZZW(:) = ZVEC3(:)
   ZZW3(:) = XFSWETH * ZZW(:)                                     & ! RSWETH
                    *  PRST(:) * PCHT(:)                          &
                    *  PRHODREF(:)**(1-XCEXVT)                    &
                    *( XLBSWETH1/( PLBDH(:)**2                ) + &
                       XLBSWETH2/( PLBDH(:)    * PLBDS(:)     ) + &
                       XLBSWETH3/(               PLBDS(:)**2) )
!
   Z1(:) = GET_XKER_N_SWETH(IVEC1(:)+1,IVEC2(:)+1)
   Z2(:) = GET_XKER_N_SWETH(IVEC1(:)+1,IVEC2(:)  )
   Z3(:) = GET_XKER_N_SWETH(IVEC1(:)  ,IVEC2(:)+1)
   Z4(:) = GET_XKER_N_SWETH(IVEC1(:)  ,IVEC2(:)  )
   ZVEC3(:) =  (      Z1(:)* ZVEC2(:)          &
                    - Z2(:)*(ZVEC2(:) - 1.0) ) &
      			 	                            *  ZVEC1(:)    &
                 - (  Z3(:)* ZVEC2(:)          &
                    - Z4(:)*(ZVEC2(:) - 1.0) ) &
       			                                    * (ZVEC1(:) - 1.0)
   ZZW(:) = ZVEC3(:)
   ZZW3N(:) = XFNSWETH * ZZW(:)                                      & ! NSWETH
                      *  PCST(:) * PCHT(:)                           &
                      *  PRHODREF(:)**(1-XCEXVT)                     &
                      *( XLBNSWETH1/( PLBDH(:)**2                ) + &
                         XLBNSWETH2/( PLBDH(:)    * PLBDS(:)     ) + &
                         XLBNSWETH3/(               PLBDS(:)**2) )
END WHERE
!
!*           1.c  Collection of rg
!            ---------------------
!
GWET(:) = PRGT(:)>XRTMIN(6) .AND. PRHT(:)>XRTMIN(7) .AND. LDCOMPUTE(:) .AND. &
          PCGT(:)>XCTMIN(6) .AND. PCHT(:)>XCTMIN(7)
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
   ZVEC1(:) = MAX( 1.0001, MIN( REAL(NWETLBDAH)-0.0001,           &
                         XWETINTP1H * LOG( ZVEC1(:) ) + XWETINTP2H ) )
   IVEC1(:) = INT( ZVEC1(:) )
   ZVEC1(:) = ZVEC1(:) - REAL( IVEC1(:) )
!
   ZVEC2(:) = MAX( 1.0001, MIN( REAL(NWETLBDAG)-0.0001,           &
                         XWETINTP1G * LOG( ZVEC2(:) ) + XWETINTP2G ) )
   IVEC2(:) = INT( ZVEC2(:) )
   ZVEC2(:) = ZVEC2(:) - REAL( IVEC2(:) )
!
!*       Perform the bilinear interpolation of the normalized
!        RDRYG-kernel
!
   Z1(:) = GET_XKER_GWETH(IVEC1(:)+1,IVEC2(:)+1)
   Z2(:) = GET_XKER_GWETH(IVEC1(:)+1,IVEC2(:)  )
   Z3(:) = GET_XKER_GWETH(IVEC1(:)  ,IVEC2(:)+1)
   Z4(:) = GET_XKER_GWETH(IVEC1(:)  ,IVEC2(:)  )
      ZVEC3(:) =  (   Z1(:)* ZVEC2(:)          &
                    - Z2(:)*(ZVEC2(:) - 1.0) ) &
                     			 	             *  ZVEC1(:)   &
                 - (  Z3(:)* ZVEC2(:)          &
                    - Z4(:)*(ZVEC2(:) - 1.0) ) &
                                 			     * (ZVEC1(:) - 1.0)
   ZZW(:) = ZVEC3(:)
   ZZW4(:) = XFGWETH * ZZW(:)                                     & ! RGWETH
                     * PRGT(:) * PCHT(:)                          &
                     * PRHODREF(:)**(1-XCEXVT)                    &
                     *( XLBGWETH1/( PLBDH(:)**2               ) + &
                        XLBGWETH2/( PLBDH(:)   * PLBDG(:)     ) + &
                        XLBGWETH3/(              PLBDG(:)**2) )
!
   Z1(:) = GET_XKER_N_GWETH(IVEC1(:)+1,IVEC2(:)+1)
   Z2(:) = GET_XKER_N_GWETH(IVEC1(:)+1,IVEC2(:)  )
   Z3(:) = GET_XKER_N_GWETH(IVEC1(:)  ,IVEC2(:)+1)
   Z4(:) = GET_XKER_N_GWETH(IVEC1(:)  ,IVEC2(:)  )
      ZVEC3(:) =  (   Z1(:)* ZVEC2(:)          &
                    - Z2(:)*(ZVEC2(:) - 1.0) ) &
                     			 	             *  ZVEC1(:)   &
                 - (  Z3(:)* ZVEC2(:)          &
                    - Z4(:)*(ZVEC2(:) - 1.0) ) &
                                 			     * (ZVEC1(:) - 1.0)
   ZZW(:) = ZVEC3(:)
   ZZW4N(:) = XFNGWETH * ZZW(:)                                      & ! NGWETH
                       * PCGT(:) * PCHT(:)                           &
                       * PRHODREF(:)**(1-XCEXVT)                     &
                       *( XLBNGWETH1/( PLBDH(:)**2               ) + &
                          XLBNGWETH2/( PLBDH(:)   * PLBDG(:)     ) + &
                          XLBNGWETH3/(              PLBDG(:)**2) )
   
END WHERE
!
!            1.e Wet growth of hail
!            ----------------------
!
ZZW(:) = 0.0
WHERE( PRHT(:)>XRTMIN(6) .AND. PCHT(:) >XCTMIN(6) .AND. LDCOMPUTE(:) )
   ZZW(:) = PRVT(:)*PPRES(:)/((XMV/XMD)+PRVT(:)) ! Vapor pressure
   ZZW(:) = PKA(:)*(XTT-PT(:)) +                                  &
              ( PDV(:)*(XLVTT + ( XCPV - XCL ) * ( PT(:) - XTT ))   &
                          *(XESTT-ZZW(:))/(XRV*PT(:))             )
!
! Total mass gained by hail in wet mode
   ZRWETH(:)  = MAX( 0.0,                                                                &
                   ( ZZW(:) * PCHT(:) * ( X0DEPH*       PLBDH(:)**XEX0DEPH +             &
                                          X1DEPH*PCJ(:)*PLBDH(:)**XEX1DEPH ) +           &
                   ( ZZW2(:)+ZZW3(:)+ZZW4(:) ) * ( XLMTT + (XCI-XCL)*(XTT-PT(:)) )   )   &
                   / (XLMTT-XCL*(XTT-PT(:)))                                     )
   ! We must agregate, at least, the cold species
   ZRWETH(:)=MAX(ZRWETH(:), ZZW2(:)+ZZW3(:)+ZZW4(:))
END WHERE
WHERE( PRHT(:)>XRTMIN(6) .AND. PCHT(:) >XCTMIN(6) .AND. PRRT(:)>XRTMIN(3) .AND. PCRT(:) >XCTMIN(3) .AND. LDCOMPUTE(:) )
   ! Mass of rain frozen by hail RRWETH
   ZZW5(:) = ZRWETH(:) - ZZW2(:) - ZZW3(:) - ZZW4(:) - ZZW1(:)
END WHERE
!
ZZW(:) = 0.0
WHERE( LDCOMPUTE(:) .AND. PT(:)<XTT .AND. ZZW5(:)>0.0 ) 
   P_RC_WETH(:) = - ZZW1(:)
   P_CC_WETH(:) = P_RC_WETH(:) * PCCT(:)/MAX(PRCT(:),XRTMIN(2))
   P_RR_WETH(:) = - ZZW5(:)
   P_CR_WETH(:) = P_RR_WETH(:) * PCRT(:)/MAX(PRRT(:),XRTMIN(3))
   P_RI_WETH(:) = - ZZW2(:)
   P_CI_WETH(:) = P_RI_WETH(:) * PCIT(:)/MAX(PRIT(:),XRTMIN(4))
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
GWET(:) = PRHT(:)<ZTHRH .AND. PRCT(:)<ZTHRC .AND. PT(:)<XTT 
WHERE( GWET(:) )
   P_RG_COHG(:) = PRHT * MIN( 1.0,MAX( 0.0,1.0-(PRCT(:)/ZTHRC) ) )
   P_CG_COHG(:) = P_RG_COHG(:) * PCHT(:)/MAX(PRHT(:),XRTMIN(7))
END WHERE
!
!
!*       3.  Hail Melting (HMLT)
!        -----------------------
!
ZZW(:) = 0.0
WHERE( PRHT(:)>XRTMIN(6) .AND. PCHT(:)>XCTMIN(7) .AND. PT(:)>XTT .AND. LDCOMPUTE(:) )
   ZZW(:) = PRVT(:)*PPRES(:)/((XMV/XMD)+PRVT(:)) ! Vapor pressure
   ZZW(:) = PKA(:)*(XTT-PT(:)) +                                 &
              ( PDV(:)*(XLVTT + ( XCPV - XCL ) * ( PT(:) - XTT )) &
                          *(XESTT-ZZW(:))/(XRV*PT(:))             )
!
! compute RHMLTR
!
   ZZW(:)  = MAX( 0.0,( -ZZW(:) * PCHT(:) *                        &
                          ( X0DEPH*       PLBDH(:)**XEX0DEPH +     &
                            X1DEPH*PCJ(:)*PLBDH(:)**XEX1DEPH ) -   &
                        ( ZZW5(:) ) * ( XCL*(XTT-PT(:))) ) &
                        / XLMTT                                    )
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
CONTAINS
  FUNCTION GET_XKER_SWETH(GRAUPEL,SNOW) RESULT(RET)
    INTEGER, DIMENSION(:) :: GRAUPEL
    INTEGER, DIMENSION(:) :: SNOW
    REAL, DIMENSION(SIZE(SNOW)) :: RET
    !
    INTEGER I
    !
    DO I=1,SIZE(GRAUPEL)
       RET(I) = XKER_SWETH(MAX(MIN(GRAUPEL(I),SIZE(XKER_SWETH,1)),1),MAX(MIN(SNOW(I),SIZE(XKER_SWETH,2)),1))
    END DO
  END FUNCTION GET_XKER_SWETH
!
!-------------------------------------------------------------------------------
!
  FUNCTION GET_XKER_N_SWETH(GRAUPEL,SNOW) RESULT(RET)
    INTEGER, DIMENSION(:) :: GRAUPEL
    INTEGER, DIMENSION(:) :: SNOW
    REAL, DIMENSION(SIZE(SNOW)) :: RET
    !
    INTEGER I
    !
    DO I=1,SIZE(GRAUPEL)
       RET(I) = XKER_N_SWETH(MAX(MIN(GRAUPEL(I),SIZE(XKER_N_SWETH,1)),1),MAX(MIN(SNOW(I),SIZE(XKER_N_SWETH,2)),1))
    END DO
  END FUNCTION GET_XKER_N_SWETH
!
!-------------------------------------------------------------------------------
!
  FUNCTION GET_XKER_GWETH(GRAUPEL,SNOW) RESULT(RET)
    INTEGER, DIMENSION(:) :: GRAUPEL
    INTEGER, DIMENSION(:) :: SNOW
    REAL, DIMENSION(SIZE(SNOW)) :: RET
    !
    INTEGER I
    !
    DO I=1,SIZE(GRAUPEL)
       RET(I) = XKER_GWETH(MAX(MIN(GRAUPEL(I),SIZE(XKER_GWETH,1)),1),MAX(MIN(SNOW(I),SIZE(XKER_GWETH,2)),1))
    END DO
  END FUNCTION GET_XKER_GWETH
!
!-------------------------------------------------------------------------------
!
  FUNCTION GET_XKER_N_GWETH(GRAUPEL,SNOW) RESULT(RET)
    INTEGER, DIMENSION(:) :: GRAUPEL
    INTEGER, DIMENSION(:) :: SNOW
    REAL, DIMENSION(SIZE(SNOW)) :: RET
    !
    INTEGER I
    !
    DO I=1,SIZE(GRAUPEL)
       RET(I) = XKER_N_GWETH(MAX(MIN(GRAUPEL(I),SIZE(XKER_N_GWETH,1)),1),MAX(MIN(SNOW(I),SIZE(XKER_N_GWETH,2)),1))
    END DO
  END FUNCTION GET_XKER_N_GWETH
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_HAIL
END MODULE MODE_LIMA_HAIL
