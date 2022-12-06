!MNH_LIC Copyright 2018-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!      #################################
       MODULE MODI_LIMA_GRAUPEL
!      #################################
!
INTERFACE
   SUBROUTINE LIMA_GRAUPEL (PTSTEP, LDCOMPUTE,                                     &
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
!
REAL, DIMENSION(:),   INTENT(IN)    :: PCCT    !
REAL, DIMENSION(:),   INTENT(IN)    :: PCRT    !
REAL, DIMENSION(:),   INTENT(IN)    :: PCIT    !
REAL, DIMENSION(:),   INTENT(IN)    :: PCST    !
REAL, DIMENSION(:),   INTENT(IN)    :: PCGT    !
!
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDC   ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDR   ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDS   ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDG   ! 
!
REAL, DIMENSION(:),   INTENT(IN)    :: PLVFACT ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PLSFACT ! 
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_TH_WETG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RC_WETG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CC_WETG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RR_WETG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CR_WETG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RI_WETG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CI_WETG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RS_WETG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CS_WETG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RG_WETG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CG_WETG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RH_WETG
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_TH_DRYG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RC_DRYG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CC_DRYG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RR_DRYG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CR_DRYG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RI_DRYG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CI_DRYG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RS_DRYG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CS_DRYG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RG_DRYG
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RI_HMG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CI_HMG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RG_HMG
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_TH_GMLT
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RR_GMLT
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CR_GMLT
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CG_GMLT
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
END SUBROUTINE LIMA_GRAUPEL
END INTERFACE
END MODULE MODI_LIMA_GRAUPEL
!
!     #################################################################################
      SUBROUTINE LIMA_GRAUPEL (PTSTEP, LDCOMPUTE,                                     &
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
USE MODD_CST,              ONLY : XTT, XMD, XMV, XRD, XRV, XLVTT, XLMTT, XESTT, XCL, XCI, XCPV
USE MODD_PARAM_LIMA,       ONLY : XRTMIN, XCEXVT, LHAIL
USE MODD_PARAM_LIMA_MIXED, ONLY : XCXG, XDG, X0DEPG, X1DEPG, NGAMINC,                             &
                                  XFCDRYG, XFIDRYG, XCOLIG, XCOLSG, XCOLEXIG, XCOLEXSG,           &
                                  XFSDRYG, XLBSDRYG1, XLBSDRYG2, XLBSDRYG3, XKER_SDRYG,           &
                                  XFNSDRYG, XLBNSDRYG1, XLBNSDRYG2, XLBNSDRYG3, XKER_N_SDRYG,     &
                                  XFRDRYG, XLBRDRYG1, XLBRDRYG2, XLBRDRYG3, XKER_RDRYG,           &
                                  XFNRDRYG, XLBNRDRYG1, XLBNRDRYG2, XLBNRDRYG3, XKER_N_RDRYG,     &
                                  XHMTMIN, XHMTMAX, XHMLINTP1, XHMLINTP2, XHM_FACTG, XGAMINC_HMC, &
                                  XEX0DEPG, XEX1DEPG,                                             &
                                  XDRYINTP1R, XDRYINTP1S, XDRYINTP1G,                             &
                                  XDRYINTP2R, XDRYINTP2S, XDRYINTP2G,                             &
                                  NDRYLBDAR, NDRYLBDAS, NDRYLBDAG
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
!
REAL, DIMENSION(:),   INTENT(IN)    :: PCCT    !
REAL, DIMENSION(:),   INTENT(IN)    :: PCRT    !
REAL, DIMENSION(:),   INTENT(IN)    :: PCIT    !
REAL, DIMENSION(:),   INTENT(IN)    :: PCST    !
REAL, DIMENSION(:),   INTENT(IN)    :: PCGT    !
!
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDC   ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDR   ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDS   ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDG   ! 
!
REAL, DIMENSION(:),   INTENT(IN)    :: PLVFACT ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PLSFACT ! 
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_TH_WETG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RC_WETG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CC_WETG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RR_WETG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CR_WETG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RI_WETG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CI_WETG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RS_WETG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CS_WETG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RG_WETG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CG_WETG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RH_WETG
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_TH_DRYG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RC_DRYG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CC_DRYG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RR_DRYG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CR_DRYG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RI_DRYG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CI_DRYG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RS_DRYG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CS_DRYG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RG_DRYG
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RI_HMG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CI_HMG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RG_HMG
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_TH_GMLT
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RR_GMLT
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CR_GMLT
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CG_GMLT
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
LOGICAL, DIMENSION(SIZE(PRCT))  :: GDRY
INTEGER                         :: IGDRY
INTEGER                         :: JJ
!
REAL,    DIMENSION(SIZE(PRCT))  :: Z1, Z2, Z3, Z4
REAL,    DIMENSION(SIZE(PRCT))  :: ZZX, ZZW, ZZW1, ZZW2, ZZW3, ZZW4, ZZW5, ZZW6, ZZW7
REAL,    DIMENSION(SIZE(PRCT))  :: ZZW3N, ZZW4N, ZZW6N 
REAL,    DIMENSION(SIZE(PRCT))  :: ZRDRYG, ZRWETG
!
INTEGER, DIMENSION(SIZE(PRCT))  :: IVEC1,IVEC2        ! Vectors of indices
REAL,    DIMENSION(SIZE(PRCT))  :: ZVEC1,ZVEC2, ZVEC3 ! Work vectors
!
INTEGER                         :: NHAIL
!
!-------------------------------------------------------------------------------
!
!
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
WHERE( PRGT(:)>XRTMIN(6) .AND. LDCOMPUTE(:) )
   ZZW(:) = PCGT(:) * PLBDG(:)**(-XDG-2.0) * PRHODREF(:)**(1-XCEXVT)
   ZZW1(:) = XFCDRYG * PRCT(:) * ZZW(:)                               ! RCDRYG - rc collected by graupel in dry mode 
   ZZW2(:) = XFIDRYG * EXP( XCOLEXIG*(PT(:)-XTT) ) * PRIT(:) * ZZW(:) ! RIDRYG - ri collected by graupel in dry mode
END WHERE
!
!*           1.b Collection of rs in the dry mode
!            ------------------------------------
!
GDRY(:) = PRST(:)>XRTMIN(5) .AND. PRGT(:)>XRTMIN(6) .AND. LDCOMPUTE(:)
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
   ZVEC1(:) = MAX( 1.0001, MIN( REAL(NDRYLBDAG)-0.0001,           &
                         XDRYINTP1G * LOG( ZVEC1(:) ) + XDRYINTP2G ) )
   IVEC1(:) = INT( ZVEC1(:) )
   ZVEC1(:) = ZVEC1(:) - REAL( IVEC1(:) )
!
   ZVEC2(:) = MAX( 1.0001, MIN( REAL(NDRYLBDAS)-0.0001,           &
                         XDRYINTP1S * LOG( ZVEC2(:) ) + XDRYINTP2S ) )
   IVEC2(:) = INT( ZVEC2(:) )
   ZVEC2(:) = ZVEC2(:) - REAL( IVEC2(:) )
!
!*       perform the bilinear interpolation of the normalized
!        SDRYG-kernel
   !
   Z1(:) = GET_XKER_SDRYG(IVEC1(:)+1,IVEC2(:)+1)
   Z2(:) = GET_XKER_SDRYG(IVEC1(:)+1,IVEC2(:)  )
   Z3(:) = GET_XKER_SDRYG(IVEC1(:)  ,IVEC2(:)+1)
   Z4(:) = GET_XKER_SDRYG(IVEC1(:)  ,IVEC2(:)  )
   ZVEC3(:) =  (      Z1(:)* ZVEC2(:)          &
                    - Z2(:)*(ZVEC2(:) - 1.0) ) &
      			 	                            *  ZVEC1(:)    &
                 - (  Z3(:)* ZVEC2(:)          &
                    - Z4(:)*(ZVEC2(:) - 1.0) ) &
       			                                    * (ZVEC1(:) - 1.0)
   ZZW(:) = ZVEC3(:)
   ZZW3(:) = XFSDRYG * ZZW(:) * EXP( XCOLEXSG*(PT(:)-XTT) )       & ! RSDRYG - rs collected by graupel in dry mode
                    *  PRST(:) * PCGT(:)                          &
                    *  PRHODREF(:)**(1-XCEXVT)                    &
                    *( XLBSDRYG1/( PLBDG(:)**2                ) + &
                       XLBSDRYG2/( PLBDG(:)    * PLBDS(:)     ) + &
                       XLBSDRYG3/(               PLBDS(:)**2) )
!
   Z1(:) = GET_XKER_N_SDRYG(IVEC1(:)+1,IVEC2(:)+1)
   Z2(:) = GET_XKER_N_SDRYG(IVEC1(:)+1,IVEC2(:)  )
   Z3(:) = GET_XKER_N_SDRYG(IVEC1(:)  ,IVEC2(:)+1)
   Z4(:) = GET_XKER_N_SDRYG(IVEC1(:)  ,IVEC2(:)  )
   ZVEC3(:) =  (      Z1(:)* ZVEC2(:)          &
                    - Z2(:)*(ZVEC2(:) - 1.0) ) &
      			 	                            *  ZVEC1(:)    &
                 - (  Z3(:)* ZVEC2(:)          &
                    - Z4(:)*(ZVEC2(:) - 1.0) ) &
       			                                    * (ZVEC1(:) - 1.0)
   ZZW(:) = ZVEC3(:)
   ZZW3N(:) = XFNSDRYG * ZZW(:) * EXP( XCOLEXSG*(PT(:)-XTT) )        & ! NSDRYG - Ns collected by graupel in dry mode
                      *  PCST(:) * PCGT(:)                           &
                      *  PRHODREF(:)**(1-XCEXVT)                     &
                      *( XLBNSDRYG1/( PLBDG(:)**2                ) + &
                         XLBNSDRYG2/( PLBDG(:)    * PLBDS(:)     ) + &
                         XLBNSDRYG3/(               PLBDS(:)**2) )
END WHERE
!
!*           1.c  Collection of rr in the dry mode
!            -------------------------------------
!
GDRY(:) = PRRT(:)>XRTMIN(3) .AND. PRGT(:)>XRTMIN(6) .AND. LDCOMPUTE(:)
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
   ZVEC1(:) = MAX( 1.0001, MIN( REAL(NDRYLBDAG)-0.0001,           &
                         XDRYINTP1G * LOG( ZVEC1(:) ) + XDRYINTP2G ) )
   IVEC1(:) = INT( ZVEC1(:) )
   ZVEC1(:) = ZVEC1(:) - REAL( IVEC1(:) )
!
   ZVEC2(:) = MAX( 1.0001, MIN( REAL(NDRYLBDAR)-0.0001,           &
                         XDRYINTP1R * LOG( ZVEC2(:) ) + XDRYINTP2R ) )
   IVEC2(:) = INT( ZVEC2(:) )
   ZVEC2(:) = ZVEC2(:) - REAL( IVEC2(:) )
!
!*       Perform the bilinear interpolation of the normalized
!        RDRYG-kernel
!
   Z1(:) = GET_XKER_RDRYG(IVEC1(:)+1,IVEC2(:)+1)
   Z2(:) = GET_XKER_RDRYG(IVEC1(:)+1,IVEC2(:)  )
   Z3(:) = GET_XKER_RDRYG(IVEC1(:)  ,IVEC2(:)+1)
   Z4(:) = GET_XKER_RDRYG(IVEC1(:)  ,IVEC2(:)  )
      ZVEC3(:) =  (   Z1(:)* ZVEC2(:)          &
                    - Z2(:)*(ZVEC2(:) - 1.0) ) &
                     			 	             *  ZVEC1(:)   &
                 - (  Z3(:)* ZVEC2(:)          &
                    - Z4(:)*(ZVEC2(:) - 1.0) ) &
                                 			     * (ZVEC1(:) - 1.0)
   ZZW(:) = ZVEC3(:)
   ZZW4(:) = XFRDRYG * ZZW(:)                                     & ! RRDRYG
                     * PRRT(:) * PCGT(:)                          &
                     * PRHODREF(:)**(1-XCEXVT)                    &
                     *( XLBRDRYG1/( PLBDG(:)**2               ) + &
                        XLBRDRYG2/( PLBDG(:)   * PLBDR(:)     ) + &
                        XLBRDRYG3/(              PLBDR(:)**2) )
!
   Z1(:) = GET_XKER_N_RDRYG(IVEC1(:)+1,IVEC2(:)+1)
   Z2(:) = GET_XKER_N_RDRYG(IVEC1(:)+1,IVEC2(:)  )
   Z3(:) = GET_XKER_N_RDRYG(IVEC1(:)  ,IVEC2(:)+1)
   Z4(:) = GET_XKER_N_RDRYG(IVEC1(:)  ,IVEC2(:)  )
      ZVEC3(:) =  (   Z1(:)* ZVEC2(:)          &
                    - Z2(:)*(ZVEC2(:) - 1.0) ) &
                     			 	             *  ZVEC1(:)   &
                 - (  Z3(:)* ZVEC2(:)          &
                    - Z4(:)*(ZVEC2(:) - 1.0) ) &
                                 			     * (ZVEC1(:) - 1.0)
   ZZW(:) = ZVEC3(:)
   ZZW4N(:) = XFNRDRYG * ZZW(:)                                      & ! NRDRYG
                       * PCRT(:) * PCGT(:)                           &
                       * PRHODREF(:)**(1-XCEXVT)                     &
                       *( XLBNRDRYG1/( PLBDG(:)**2               ) + &
                          XLBNRDRYG2/( PLBDG(:)   * PLBDR(:)     ) + &
                          XLBNRDRYG3/(              PLBDR(:)**2) )
   
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
WHERE( PRGT(:)>XRTMIN(6) .AND. LDCOMPUTE(:) )
   ZZW5(:) = ZZW2(:) / (XCOLIG*EXP(XCOLEXIG*(PT(:)-XTT)) ) ! RIWETG
   ZZW6(:) = ZZW3(:) / (XCOLSG*EXP(XCOLEXSG*(PT(:)-XTT)) ) ! RSWETG
   ZZW6N(:)= ZZW3N(:)/ (XCOLSG*EXP(XCOLEXSG*(PT(:)-XTT)) ) ! NSWETG
!
   ZZW(:) = PRVT(:)*PPRES(:)/((XMV/XMD)+PRVT(:)) ! Vapor pressure
   ZZW(:) = PKA(:)*(XTT-PT(:)) +                                  &
              ( PDV(:)*(XLVTT + ( XCPV - XCL ) * ( PT(:) - XTT ))   &
                          *(XESTT-ZZW(:))/(XRV*PT(:))             )
!
! Total mass gained by graupel in wet mode
   ZRWETG(:)  = MAX( 0.0,                                                        &
                   ( ZZW(:) * PCGT(:) * ( X0DEPG*       PLBDG(:)**XEX0DEPG +     &
                                          X1DEPG*PCJ(:)*PLBDG(:)**XEX1DEPG ) +   &
                   ( ZZW5(:)+ZZW6(:) ) * ( XLMTT + (XCI-XCL)*(XTT-PT(:)) )   )   &
                   / (XLMTT-XCL*(XTT-PT(:)))                                     )
  !We must agregate, at least, the cold species
   ZRWETG(:)=MAX(ZRWETG(:), ZZW5(:)+ZZW6(:))
END WHERE
!
!            1.f Wet mode and partial conversion to hail
!            -------------------------------------------
!
ZZW(:) = 0.0
NHAIL = 0.
IF (LHAIL) NHAIL = 1. 
WHERE( LDCOMPUTE(:) .AND. PRGT(:)>XRTMIN(6) .AND. PT(:)<XTT .AND. &
       (ZRDRYG(:)-ZZW2(:)-ZZW3(:))>=(ZRWETG(:)-ZZW5(:)-ZZW6(:)) .AND. ZRWETG(:)-ZZW5(:)-ZZW6(:)>0.0 ) 
!
! Mass of rain and cloud droplets frozen by graupel in wet mode : RCWETG + RRWETG = RWETG - RIWETG - RSWETG
   ZZW7(:) = ZRWETG(:) - ZZW5(:) - ZZW6(:)
!
! assume a linear percent of conversion of graupel into hail
! ZZW = percentage of graupel transformed
!
   ZZW(:) = ZRDRYG(:)*NHAIL/(ZRWETG(:)+ZRDRYG(:)) 
!
   P_RC_WETG(:) = - ZZW1(:)
   P_CC_WETG(:) = P_RC_WETG(:) * PCCT(:)/MAX(PRCT(:),XRTMIN(2))
   P_RR_WETG(:) = - ZZW7(:) + ZZW1(:)
   P_CR_WETG(:) = P_RR_WETG(:) * PCRT(:)/MAX(PRRT(:),XRTMIN(3))
   P_RI_WETG(:) = - ZZW5(:)
   P_CI_WETG(:) = P_RI_WETG(:) * PCIT(:)/MAX(PRIT(:),XRTMIN(4))
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
WHERE( LDCOMPUTE(:) .AND. PRGT(:)>XRTMIN(6) .AND. PT(:)<XTT .AND.                  &
       (ZRDRYG(:)-ZZW2(:)-ZZW3(:))<(ZRWETG(:)-ZZW5(:)-ZZW6(:)) .AND. ZRDRYG(:)>0.0 )
   !
   P_RC_DRYG(:) = - ZZW1(:)
   P_CC_DRYG(:) = P_RC_DRYG(:) * PCCT(:)/MAX(PRCT(:),XRTMIN(2))
   P_RR_DRYG(:) = - ZZW4(:)
   P_CR_DRYG(:) = - ZZW4N(:)
   P_RI_DRYG(:) = - ZZW2(:)
   P_CI_DRYG(:) = P_RI_DRYG(:) * PCIT(:)/MAX(PRIT(:),XRTMIN(4))
   P_RS_DRYG(:) = - ZZW3(:)
   P_CS_DRYG(:) = - ZZW3N(:)
   P_RG_DRYG(:) =   ZRDRYG(:)
   !
   P_TH_DRYG(:) = (ZZW1(:) + ZZW4(:)) * (PLSFACT(:)-PLVFACT(:))
END WHERE
!
!
!*       2.  Hallett-Mossop process (HMG)
!        --------------------------------
!
! BVIE test ZRDRYG<ZZW ?????????????????????????
!GDRY(:) = (PT(:)<XHMTMAX) .AND. (PT(:)>XHMTMIN)    .AND. (ZRDRYG(:)<ZZW(:))&
GDRY(:) = PT(:)<XHMTMAX .AND. PT(:)>XHMTMIN .AND. PRGT(:)>XRTMIN(6) .AND. PRCT(:)>XRTMIN(2) .AND. LDCOMPUTE(:) .AND. &
          (ZRDRYG(:)-ZZW2(:)-ZZW3(:))<(ZRWETG(:)-ZZW5(:)-ZZW6(:))

ZZX(:)=9999.
ZVEC1(:)=0.
ZVEC2(:)=0.
IVEC1(:)=0
IVEC2(:)=0
WHERE( GDRY(:) )
!
   ZVEC1(:) = PLBDC(:)
   ZVEC2(:) = MAX( 1.0001, MIN( REAL(NGAMINC)-0.0001,           &
                         XHMLINTP1 * LOG( ZVEC1(:) ) + XHMLINTP2 ) )
   IVEC2(:) = INT( ZVEC2(:) )
   ZVEC2(:) = ZVEC2(:) - REAL( IVEC2(:) )
   ZVEC1(:) =   XGAMINC_HMC( IVEC2(:)+1 )* ZVEC2(:)      &
                    - XGAMINC_HMC( IVEC2(:)   )*(ZVEC2(:) - 1.0)
   ZZX(:) = ZVEC1(:) ! Large droplets
!
   WHERE ( ZZX(:)<0.99 ) ! Dry case
      P_CI_HMG(:) = ZZW1(:)*(PCCT(:)/PRCT(:))*(1.0-ZZX(:))*XHM_FACTG*  &
                        MAX( 0.0, MIN( (PT(:)-XHMTMIN)/3.0,(XHMTMAX-PT(:))/2.0 ) )
      P_RI_HMG(:) = P_CI_HMG(:) * XMNU0
      P_RG_HMG(:) = - P_RI_HMG(:)
   END WHERE
END WHERE
!
!
!*       3.  Graupel Melting
!        -------------------
!
ZZX(:) = 0.0
WHERE( (PRGT(:)>XRTMIN(6)) .AND. (PT(:)>XTT) .AND. LDCOMPUTE(:) )
   ZZX(:) = PRVT(:)*PPRES(:)/((XMV/XMD)+PRVT(:)) ! Vapor pressure
   ZZX(:) = PKA(:)*(XTT-PT(:)) +                                 &
              ( PDV(:)*(XLVTT + ( XCPV - XCL ) * ( PT(:) - XTT )) &
                          *(XESTT-ZZX(:))/(XRV*PT(:))             )
!
! compute RGMLTR
!
   ZZX(:)  = MAX( 0.0,( -ZZX(:) * PCGT(:) *                        &
                          ( X0DEPG*       PLBDG(:)**XEX0DEPG +     &
                            X1DEPG*PCJ(:)*PLBDG(:)**XEX1DEPG ) -   &
                        ( ZZW1(:)+ZZW4(:) ) * ( XCL*(XTT-PT(:))) ) &
                        / XLMTT                                    )
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
CONTAINS
  FUNCTION GET_XKER_SDRYG(GRAUPEL,SNOW) RESULT(RET)
    INTEGER, DIMENSION(:) :: GRAUPEL
    INTEGER, DIMENSION(:) :: SNOW
    REAL, DIMENSION(SIZE(SNOW)) :: RET
    !
    INTEGER I
    !
    DO I=1,SIZE(GRAUPEL)
       RET(I) = XKER_SDRYG(MAX(MIN(GRAUPEL(I),SIZE(XKER_SDRYG,1)),1),MAX(MIN(SNOW(I),SIZE(XKER_SDRYG,2)),1))
    END DO
  END FUNCTION GET_XKER_SDRYG
!
!-------------------------------------------------------------------------------
!
  FUNCTION GET_XKER_N_SDRYG(GRAUPEL,SNOW) RESULT(RET)
    INTEGER, DIMENSION(:) :: GRAUPEL
    INTEGER, DIMENSION(:) :: SNOW
    REAL, DIMENSION(SIZE(SNOW)) :: RET
    !
    INTEGER I
    !
    DO I=1,SIZE(GRAUPEL)
       RET(I) = XKER_N_SDRYG(MAX(MIN(GRAUPEL(I),SIZE(XKER_N_SDRYG,1)),1),MAX(MIN(SNOW(I),SIZE(XKER_N_SDRYG,2)),1))
    END DO
  END FUNCTION GET_XKER_N_SDRYG
!
!-------------------------------------------------------------------------------
!
  FUNCTION GET_XKER_RDRYG(GRAUPEL,RAIN) RESULT(RET)
    INTEGER, DIMENSION(:) :: GRAUPEL
    INTEGER, DIMENSION(:) :: RAIN
    REAL, DIMENSION(SIZE(RAIN)) :: RET
    !
    INTEGER I
    !
    DO I=1,SIZE(GRAUPEL)
       RET(I) = XKER_RDRYG(MAX(MIN(GRAUPEL(I),SIZE(XKER_RDRYG,1)),1),MAX(MIN(RAIN(I),SIZE(XKER_RDRYG,2)),1))
    END DO
  END FUNCTION GET_XKER_RDRYG
!
!-------------------------------------------------------------------------------
!
  FUNCTION GET_XKER_N_RDRYG(GRAUPEL,RAIN) RESULT(RET)
    INTEGER, DIMENSION(:) :: GRAUPEL
    INTEGER, DIMENSION(:) :: RAIN
    REAL, DIMENSION(SIZE(RAIN)) :: RET
    !
    INTEGER I
    !
    DO I=1,SIZE(GRAUPEL)
       RET(I) = XKER_N_RDRYG(MAX(MIN(GRAUPEL(I),SIZE(XKER_N_RDRYG,1)),1),MAX(MIN(RAIN(I),SIZE(XKER_N_RDRYG,2)),1))
    END DO
  END FUNCTION GET_XKER_N_RDRYG
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_GRAUPEL
