!MNH_LIC Copyright 2018-2024 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_LIMA_RAIN_ACCR_SNOW
  IMPLICIT NONE
CONTAINS
!     ######################################################################################
  SUBROUTINE LIMA_RAIN_ACCR_SNOW (CST, LIMAP, LIMAW, LIMAC, LIMAM, KSIZE, PTSTEP, ODCOMPUTE,&
                                  PRHODREF, PT,                                             &
                                  PRCT, PCCT, PRRT, PCRT, PRST, PCST, PLBDC, PLBDR, PLBDS,  &
                                  PCF1D, PPF1D, PLVFACT, PLSFACT,                           &
                                  P_TH_RIM, P_CC_RIM, P_CS_RIM,                             &
                                  P_RC_RIMSS, P_RC_RIMSG, P_RS_RIMCG,                       &
                                  P_RI_HMS, P_CI_HMS, P_RS_HMS,                             &
                                  P_TH_ACC, P_CR_ACC, P_CS_ACC,                             &
                                  P_RR_ACCSS, P_RR_ACCSG, P_RS_ACCRG                        )
!--cb--
!     ######################################################################################
!
!!    PURPOSE
!!    -------
!!      Compute the rain drops accretion on aggregates
!!
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
!  C. Barthe   04/07/2022: modify the microphysics terms to save to simplify the merging
!                          with the electrification scheme
! QUESTION : ne fonctionne pas si NMOM_R=1 ???
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAM_LIMA_MIXED, ONLY:PARAM_LIMA_MIXED_T
USE MODD_PARAM_LIMA_COLD, ONLY:PARAM_LIMA_COLD_T
USE MODD_PARAM_LIMA_WARM, ONLY:PARAM_LIMA_WARM_T
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
TYPE(PARAM_LIMA_WARM_T),INTENT(IN)::LIMAW
TYPE(PARAM_LIMA_T),INTENT(IN)::LIMAP
TYPE(CST_T),INTENT(IN)::CST
INTEGER, INTENT(IN) :: KSIZE
REAL,                 INTENT(IN)    :: PTSTEP 
LOGICAL, DIMENSION(KSIZE),INTENT(IN)    :: ODCOMPUTE
!
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRHODREF    ! 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PT   ! 
!
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRCT    ! Cloud water mr at t
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PCCT    ! Cloud water C. at t
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRRT    ! Rain mr at t
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PCRT    ! Rain C. at t
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRST    ! Snow mr at t
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PCST    ! Snow C. at t
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLBDC   ! 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLBDR   ! 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLBDS   ! 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PCF1D   ! 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PPF1D   ! 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLVFACT ! 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLSFACT ! 
!
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_CC_RIM
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_CS_RIM
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_RC_RIMSS
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_RC_RIMSG
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_RS_RIMCG
!
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_TH_RIM
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_RI_HMS
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_CI_HMS
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_RS_HMS
!
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_TH_ACC
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_CR_ACC
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_CS_ACC
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_RR_ACCSS
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_RR_ACCSG
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_RS_ACCRG
!
!*       0.2   Declarations of local variables :
!
LOGICAL, DIMENSION(SIZE(PRRT))  :: GACC
!
REAL,    DIMENSION(SIZE(PRRT))  :: Z1, Z2, Z3, Z4
REAL,    DIMENSION(SIZE(PRRT))  :: ZZW1, ZZW2, ZZW3, ZZW4, ZZW5
REAL,    DIMENSION(SIZE(PRRT))  :: ZZWC1, ZZWC2, ZZWC3, ZZWC4, ZZWC5
REAL,    DIMENSION(SIZE(PRCT))  :: Z_RC_RIM
!
INTEGER                         :: II
INTEGER, DIMENSION(SIZE(PRRT))  :: IVEC1,IVEC2              ! Vectors of indices
REAL,    DIMENSION(SIZE(PRRT))  :: ZVEC1,ZVEC1W,ZVEC2,ZVEC3 ! Work vectors
REAL,    DIMENSION(SIZE(PRRT))  :: Z_RR_ACC  ! ++cb-- for elec
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('LIMA_RAIN_ACCR_SNOW', 0, ZHOOK_HANDLE)
!
P_TH_RIM(:) = 0.
P_RC_RIMSS(:) = 0.
P_RC_RIMSG(:) = 0.
P_RS_RIMCG(:) = 0.
Z_RC_RIM(:) = 0.
P_CC_RIM(:) = 0.
P_CS_RIM(:) = 0.
P_RI_HMS(:) = 0.
P_CI_HMS(:) = 0.
P_RS_HMS(:) = 0.
!
P_TH_ACC(:) = 0.
P_CR_ACC(:) = 0.
P_CS_ACC(:) = 0.
P_RR_ACCSS(:) = 0.
P_RR_ACCSG(:) = 0.
P_RS_ACCRG(:) = 0.
Z_RR_ACC(:) = 0.
!
ZZW1(:) = 0.
ZZW2(:) = 0.
ZZW3(:) = 0.
ZZW4(:) = 0.
ZZW5(:) = 0.
!
ZZWC1(:) = 0.
ZZWC2(:) = 0.
ZZWC3(:) = 0.
ZZWC4(:) = 0.
ZZWC5(:) = 0.
!
IVEC1(:) = 0
IVEC2(:) = 0
ZVEC1(:) = 0.
ZVEC1W(:)= 0.
ZVEC2(:) = 0.
ZVEC3(:) = 0.
!
!
!*       Cloud droplet riming of the aggregates  
!        --------------------------------------
!
!
DO II = 1, SIZE(PRCT)
   IF ( PRCT(II)/PCF1D(II)>LIMAP%XRTMIN(2) .AND. PRST(II)/PPF1D(II)>LIMAP%XRTMIN(5) .AND. PT(II)<CST%XTT .AND. &
        PCCT(II)/PCF1D(II)>LIMAP%XCTMIN(2) .AND. PCST(II/PPF1D(II))>LIMAP%XCTMIN(5) .AND. ODCOMPUTE(II) ) THEN
!
      ZVEC1(II) = PLBDS(II)
      ZVEC1W(II)= ( LIMAC%XFVELOS**LIMAP%XALPHAS + PLBDS(II)**LIMAP%XALPHAS ) ** (1./LIMAP%XALPHAS) ! modified equivalent lambda
!
!        2.     perform the linear interpolation of the normalized
!               "2+LIMAC%XDS"-moment of the incomplete gamma function using the modified equivalent lambda
!
      ZVEC2(II) = MAX( 1.0001, MIN( REAL(LIMAM%NGAMINC)-0.0001,           &
                       LIMAM%XRIMINTP1 * LOG( ZVEC1W(II) ) + LIMAM%XRIMINTP2 ) )
      IVEC2(II) = INT( ZVEC2(II) )
      ZVEC2(II) = ZVEC2(II) - REAL( IVEC2(II) )
!
      ZZW1(II)  =  LIMAM%XGAMINC_RIM1( IVEC2(II)+1 )* ZVEC2(II)      &
                - LIMAM%XGAMINC_RIM1( IVEC2(II)   )*(ZVEC2(II) - 1.0)
!
!        3.     perform the linear interpolation of the normalized
!               "LIMAC%XBS"-moment of the incomplete gamma function
!
      ZVEC2(II) = MAX( 1.0001, MIN( REAL(LIMAM%NGAMINC)-0.0001,           &
                       LIMAM%XRIMINTP1 * LOG( ZVEC1W(II) ) + LIMAM%XRIMINTP2 ) )
      IVEC2(II) = INT( ZVEC2(II) )
      ZVEC2(II) = ZVEC2(II) - REAL( IVEC2(II) )
!
      ZZW2(II)  = LIMAM%XGAMINC_RIM2( IVEC2(II)+1 )* ZVEC2(II)      &
                - LIMAM%XGAMINC_RIM2( IVEC2(II)   )*(ZVEC2(II) - 1.0)
!
!        4.     riming
!
   ! Cloud droplets collected
   ! total mass loss of cloud droplets, < 0
      Z_RC_RIM(II) = - LIMAM%XCRIMSS  * PRCT(II)/PCF1D(II) * PCST(II)/PPF1D(II)*(1+(LIMAC%XFVELOS/PLBDS(II))&
           **LIMAP%XALPHAS)**(-LIMAP%XNUS+LIMAM%XEXCRIMSS/LIMAP%XALPHAS) &
           * PRHODREF(II)**(-LIMAP%XCEXVT+1) * PLBDS(II)**LIMAM%XEXCRIMSS
      P_CC_RIM(II) = Z_RC_RIM(II) * (PCCT(II) / PRCT(II)) ! Lambda_c**3
    !
    ! Cloud droplets collected on small aggregates add to snow
!      P_RS_RIM(II) = - P_RC_RIM(II) * ZZW1(II)
      P_RC_RIMSS(II) = Z_RC_RIM(II) * ZZW1(II)  ! < 0, loss of mass for rc
    !
    ! Cloud droplets collected on large aggregates add to graupel
      P_RC_RIMSG(II) = Z_RC_RIM(II) - P_RC_RIMSS(II) ! < 0, loss of mass for rc
    !
      IF (LIMAP%LMURAKAMI) THEN
    ! Graupel formation based on Murakami
         ZVEC1(II) = LIMAM%XGAMINC_RIM4( IVEC2(II)+1 )* ZVEC2(II)      &
                   - LIMAM%XGAMINC_RIM4( IVEC2(II)   )*(ZVEC2(II) - 1.0)
         ZZW5(II) = ZVEC1(II)
         ZZW3(II) = LIMAM%XSRIMCG * PRHODREF(II) * PCST(II)/PPF1D(II) * PLBDS(II)**LIMAM%XEXSRIMCG * (1.0 - ZZW2(II))!/(PTSTEP*PRHODREF(II))
         ZZW3(II) = - P_RC_RIMSG(II)*ZZW3(II)/ &
                    MAX(1.E-20, & !-20
                        LIMAM%XSRIMCG3*LIMAM%XSRIMCG2*PCST(II)/PPF1D(II)*PRHODREF(II)*PLBDS(II)**(LIMAM%XEXSRIMCG2)*(1.-ZZW5(II))- &
                        LIMAM%XSRIMCG3*ZZW3(II))
      ELSE
    ! Large aggregates collecting droplets add to graupel (instant process ???)
         ZZW3(II) = PRST(II)/PPF1D(II)*(1.0 - ZZW2(II))/PTSTEP
      END IF
    !
      P_RS_RIMCG(II) = ZZW3(II) !!! RSRIMCG is a loss for snow, but must be positive value for electricity !!!
      P_CS_RIM(II) = -ZZW3(II) * PCST(II)/PRST(II)
    !
      P_TH_RIM(II) = - Z_RC_RIM(II)*(PLSFACT(II)-PLVFACT(II))
   END IF
!
!
!*       Hallett-Mossop ice production (HMS)  
!        -----------------------------------
!
   IF ( PRST(II)/PPF1D(II)>LIMAP%XRTMIN(5) .AND. PRCT(II)/PCF1D(II)>LIMAP%XRTMIN(2) .AND. &
        PT(II)<LIMAM%XHMTMAX .AND. PT(II)>LIMAM%XHMTMIN .AND. &
        PCST(II)/PPF1D(II)>LIMAP%XCTMIN(5) .AND. PCCT(II)/PCF1D(II)>LIMAP%XCTMIN(2) .AND. ODCOMPUTE(II) ) THEN
!
      ZVEC1(II) = PLBDC(II)
      ZVEC2(II) = MAX( 1.0001, MIN( REAL(LIMAM%NGAMINC)-0.0001,           &
                       LIMAM%XHMLINTP1 * LOG( ZVEC1(II) ) + LIMAM%XHMLINTP2 ) )
      IVEC2(II) = INT( ZVEC2(II) )
      ZVEC2(II) = ZVEC2(II) - REAL( IVEC2(II) )
      ZVEC1(II) =  LIMAM%XGAMINC_HMC( IVEC2(II)+1 )* ZVEC2(II)      &
                - LIMAM%XGAMINC_HMC( IVEC2(II)   )*(ZVEC2(II) - 1.0)
      ZZW4(II) = ZVEC1(II) ! Large droplets
!
      IF ( ZZW4(II)<0.99 ) THEN
         P_CI_HMS(II) = - Z_RC_RIM(II) * (PCCT(II)/PRCT(II)) * (1.0-ZZW4(II)) * LIMAM%XHM_FACTS * &
                       MAX( 0.0, MIN( (PT(II)-LIMAM%XHMTMIN)/3.0,(LIMAM%XHMTMAX-PT(II))/2.0 ) ) ! CCHMSI
!
         P_RI_HMS(II) = P_CI_HMS(II) * LIMAC%XMNU0                                     ! RCHMSI
         P_RS_HMS(II) = - P_RI_HMS(II)
      ELSE
         P_RI_HMS(II) = 0.
         P_CI_HMS(II) = 0.
         P_RS_HMS(II) = 0.
      END IF
   END IF
END DO
!
!
!*       Rain drops accretion on the aggregates  
!        --------------------------------------
!
!
ZZW1(:) = 0.
ZZW2(:) = 0.
ZZW3(:) = 0.
ZZW4(:) = 0.
ZZW5(:) = 0.
IVEC1(:) = 0
IVEC2(:) = 0
ZVEC1(:) = 0.
ZVEC1W(:)= 0.
ZVEC2(:) = 0.
ZVEC3(:) = 0.
!
GACC(:) = .FALSE.
GACC(:) = (PRRT(:)/PPF1D(:)>LIMAP%XRTMIN(3)) .AND. (PRST(:)/PPF1D(:)>LIMAP%XRTMIN(5)) .AND. &
           (PT(:)<CST%XTT) .AND. ODCOMPUTE(:) .AND. &
          (PCRT(:)/PPF1D(:)>LIMAP%XCTMIN(3)) .AND. (PCST(:)/PPF1D(:)>LIMAP%XCTMIN(5))
!
WHERE( GACC )
!
!        1.3.1  select the (ZLBDAS,ZLBDAR) couplet
   !
   ZVEC1(:) = MAX(MIN(PLBDS(:),5.E5*LIMAC%XTRANS_MP_GAMMAS),5.E1*LIMAC%XTRANS_MP_GAMMAS)
   ZVEC2(:) = PLBDR(:)
!
!        1.3.2  find the next lower indice for the ZLBDAS and for the ZLBDAR
!               in the geometrical set of (Lbda_s,Lbda_r) couplet use to
!               tabulate the RACCSS-kernel
!
   ZVEC1(:) = MAX( 1.0001, MIN( REAL(LIMAM%NACCLBDAS)-0.0001,           &
                         LIMAM%XACCINTP1S * LOG( ZVEC1(:) ) + LIMAM%XACCINTP2S ) )
   IVEC1(:) = INT( ZVEC1(:) )
   ZVEC1(:) = ZVEC1(:) - REAL( IVEC1(:) )
!
   ZVEC2(:) = MAX( 1.0001, MIN( REAL(LIMAM%NACCLBDAR)-0.0001,           &
                         LIMAM%XACCINTP1R * LOG( ZVEC2(:) ) + LIMAM%XACCINTP2R ) )
   IVEC2(:) = INT( ZVEC2(:) )
   ZVEC2(:) = ZVEC2(:) - REAL( IVEC2(:) )
!
!        1.3.3  perform the bilinear interpolation of the normalized
!               RACCSS-kernel : for small rain drops transformed into snow
!
   Z1(:) = GET_XKER_RACCSS(KSIZE,IVEC1(:)+1,IVEC2(:)+1)
   Z2(:) = GET_XKER_RACCSS(KSIZE,IVEC1(:)+1,IVEC2(:)  )
   Z3(:) = GET_XKER_RACCSS(KSIZE,IVEC1(:)  ,IVEC2(:)+1)
   Z4(:) = GET_XKER_RACCSS(KSIZE,IVEC1(:)  ,IVEC2(:)  )
      ZVEC3(:) =  (   Z1(:)* ZVEC2(:)          &
                    - Z2(:)*(ZVEC2(:) - 1.0) ) &
                                               *  ZVEC1(:)    &
                 - (  Z3(:)* ZVEC2(:)          &
                    - Z4(:)*(ZVEC2(:) - 1.0) ) &
                                                       * (ZVEC1(:) - 1.0)
   ZZW1(:) = ZVEC3(:)
!
!        1.3.3b perform the bilinear interpolation of the normalized
!               RACCSS-kernel for concentration : for small rain drops transformed into snow
!
!!$   Z1(:) = GET_XKER_N_RACCSS(IVEC1(:)+1,IVEC2(:)+1)
!!$   Z2(:) = GET_XKER_N_RACCSS(IVEC1(:)+1,IVEC2(:)  )
!!$   Z3(:) = GET_XKER_N_RACCSS(IVEC1(:)  ,IVEC2(:)+1)
!!$   Z4(:) = GET_XKER_N_RACCSS(IVEC1(:)  ,IVEC2(:)  )
!!$      ZVEC3(:) =  (   Z1(:)* ZVEC2(:)          &
!!$                    - Z2(:)*(ZVEC2(:) - 1.0) ) &
!!$                                               *  ZVEC1(:)    &
!!$                 - (  Z3(:)* ZVEC2(:)          &
!!$                    - Z4(:)*(ZVEC2(:) - 1.0) ) &
!!$                                                       * (ZVEC1(:) - 1.0)
!!$   ZZWC1(:) = ZVEC3(:)
!
!        1.3.4  perform the bilinear interpolation of the normalized
!               RACCS-kernel : total frozen rain drops
!
   Z1(:) = GET_XKER_RACCS(KSIZE,IVEC1(:)+1,IVEC2(:)+1)
   Z2(:) = GET_XKER_RACCS(KSIZE,IVEC1(:)+1,IVEC2(:)  )
   Z3(:) = GET_XKER_RACCS(KSIZE,IVEC1(:)  ,IVEC2(:)+1)
   Z4(:) = GET_XKER_RACCS(KSIZE,IVEC1(:)  ,IVEC2(:)  )
      ZVEC3(:) =  (    Z1(:)* ZVEC2(:)          &
                    -  Z2(:)*(ZVEC2(:) - 1.0) ) &
                                                           *  ZVEC1(:)      &
                 - (   Z3(:)* ZVEC2(:)          &
                    -  Z4(:)*(ZVEC2(:) - 1.0) ) &
                                                           * (ZVEC1(:) - 1.0)
   ZZW2(:) = ZVEC3(:)
!
!        1.3.4b perform the bilinear interpolation of the normalized
!               RACCS-kernel for concentration : total frozen rain drops
!
   Z1(:) = GET_XKER_N_RACCS(KSIZE,IVEC1(:)+1,IVEC2(:)+1)
   Z2(:) = GET_XKER_N_RACCS(KSIZE,IVEC1(:)+1,IVEC2(:)  )
   Z3(:) = GET_XKER_N_RACCS(KSIZE,IVEC1(:)  ,IVEC2(:)+1)
   Z4(:) = GET_XKER_N_RACCS(KSIZE,IVEC1(:)  ,IVEC2(:)  )
      ZVEC3(:) =  (    Z1(:)* ZVEC2(:)          &
                    -  Z2(:)*(ZVEC2(:) - 1.0) ) &
                                                           *  ZVEC1(:)      &
                 - (   Z3(:)* ZVEC2(:)          &
                    -  Z4(:)*(ZVEC2(:) - 1.0) ) &
                                                           * (ZVEC1(:) - 1.0)
   ZZWC2(:) = ZVEC3(:)
!
! Correction of ZZW1 to ensure that ZZW1 <= ZZW2
! ie                coll. of small drops <= coll. of all drops
!
   ZZW1(:) = MIN(ZZW1(:),ZZW2(:))
!!$   ZZWC1(:)= MIN(ZZWC1(:),ZZWC2(:))
!
!        1.3.5  perform the bilinear interpolation of the normalized
!               SACCRG-kernel : snow transformed into graupel
!
   Z1(:) = GET_XKER_SACCRG(KSIZE,IVEC2(:)+1,IVEC1(:)+1)
   Z2(:) = GET_XKER_SACCRG(KSIZE,IVEC2(:)+1,IVEC1(:)  )
   Z3(:) = GET_XKER_SACCRG(KSIZE,IVEC2(:)  ,IVEC1(:)+1)
   Z4(:) = GET_XKER_SACCRG(KSIZE,IVEC2(:)  ,IVEC1(:)  )
      ZVEC3(:) =  (   Z1(:)* ZVEC1(:)          &
                    - Z2(:)*(ZVEC1(:) - 1.0) ) &
                                                        *  ZVEC2(:)    &
                 - (  Z3(:)* ZVEC1(:)          &
                    - Z4(:)*(ZVEC1(:) - 1.0) ) &
                                                    * (ZVEC2(:) - 1.0)
   ZZW3(:) = ZVEC3(:)
!
!        1.3.5b perform the bilinear interpolation of the normalized
!               SACCRG-kernel for concentration : snow transformed into graupel
!
   Z1(:) = GET_XKER_N_SACCRG(KSIZE,IVEC2(:)+1,IVEC1(:)+1)
   Z2(:) = GET_XKER_N_SACCRG(KSIZE,IVEC2(:)+1,IVEC1(:)  )
   Z3(:) = GET_XKER_N_SACCRG(KSIZE,IVEC2(:)  ,IVEC1(:)+1)
   Z4(:) = GET_XKER_N_SACCRG(KSIZE,IVEC2(:)  ,IVEC1(:)  )
      ZVEC3(:) =  (   Z1(:)* ZVEC1(:)          &
                    - Z2(:)*(ZVEC1(:) - 1.0) ) &
                                                        *  ZVEC2(:)    &
                 - (  Z3(:)* ZVEC1(:)          &
                    - Z4(:)*(ZVEC1(:) - 1.0) ) &
                                                    * (ZVEC2(:) - 1.0)
   ZZWC3(:) = ZVEC3(:)
!
!        1.3.4  raindrop accretion on the small sized aggregates
!      
   ZZW4(:) = PCRT(:)/PPF1D(:) *                                                      & !! coef of RRACCS and RRACCS
            LIMAM%XFRACCSS * PCST(:)/PPF1D(:) * PRHODREF(:)**(1-LIMAP%XCEXVT)                    &
         *( LIMAM%XLBRACCS1/( PLBDS(:)**2               ) +                       &
            LIMAM%XLBRACCS2/( PLBDS(:)    * PLBDR(:)    ) +                       &
            LIMAM%XLBRACCS3/(               PLBDR(:)**2 ) ) / PLBDR(:)**LIMAW%XBR
!
   ZZWC4(:)= PCRT(:)/PPF1D(:) *                                                      & !! coef of RRACCS and RRACCS
            LIMAM%XFNRACCSS * PCST(:)/PPF1D(:) * PRHODREF(:)**(1-LIMAP%XCEXVT)                    &
         *( LIMAM%XLBNRACCS1/( PLBDS(:)**2               ) +                       &
            LIMAM%XLBNRACCS2/( PLBDS(:)    * PLBDR(:)    ) +                       &
            LIMAM%XLBNRACCS3/(               PLBDR(:)**2 ) )

!
!        1.3.6  raindrop accretion-conversion of the large sized aggregates
!               into graupeln
!
   ZZW5(:) = LIMAM%XFSACCRG * ZZW3(:) * PCRT(:)/PPF1D(:) *             & ! RSACCRG
             PCST(:)/PPF1D(:) * PLBDS(:)**(-LIMAC%XBS) * PRHODREF(:)**(1-LIMAP%XCEXVT) * &
            ( LIMAM%XLBSACCR1/( PLBDR(:)**2               ) + &
              LIMAM%XLBSACCR2/( PLBDR(:)    * PLBDS(:)    ) + &
              LIMAM%XLBSACCR3/(               PLBDS(:)**2 ) )
!
   ZZWC5(:)= LIMAM%XFNSACCRG * ZZWC3(:) * PCRT(:)/PPF1D(:) *            & ! RSACCRG
             PCST(:)/PPF1D(:) * PRHODREF(:)**(1-LIMAP%XCEXVT) * &
            ( LIMAM%XLBNSACCR1/( PLBDR(:)**2               ) + &
              LIMAM%XLBNSACCR2/( PLBDR(:)    * PLBDS(:)    ) + &
              LIMAM%XLBNSACCR3/(               PLBDS(:)**2 ) )
!
   Z_RR_ACC(:) = - ZZW4(:) *  ZZW2(:) ! < 0
   P_CR_ACC(:) = - ZZWC4(:) * ZZWC2(:)
   P_CS_ACC(:) = - ZZWC5(:)
   P_RR_ACCSS(:) = ZZW4(:) *  ZZW1(:)  ! perte pour rr, > 0
   P_RR_ACCSG(:) = ZZW4(:) * ( ZZW2(:) - ZZW1(:) )  ! rraccsg = rraccs - rraccss
   P_RS_ACCRG(:) = ZZW5(:)             ! perte pour rs, > 0
   P_TH_ACC(:) = - Z_RR_ACC(:) * (PLSFACT(:)-PLVFACT(:))
!
END WHERE
!
!
!
P_CC_RIM(:) = P_CC_RIM(:) * PCF1D(:)
P_CS_RIM(:) = P_CS_RIM(:) * PCF1D(:)
P_RC_RIMSS(:) = P_RC_RIMSS(:) * PCF1D(:)
P_RC_RIMSG(:) = P_RC_RIMSG(:) * PCF1D(:)
P_RS_RIMCG(:) = P_RS_RIMCG(:) * PCF1D(:)
P_TH_RIM(:) = - (P_RC_RIMSS(:) + P_RC_RIMSG(:)) * (PLSFACT(:)-PLVFACT(:))
P_RI_HMS(:) = P_RI_HMS(:) * PCF1D(:)
P_CI_HMS(:) = P_CI_HMS(:) * PCF1D(:)
P_RS_HMS(:) = P_RS_HMS(:) * PCF1D(:)
!
P_CR_ACC(:) = P_CR_ACC(:) * PPF1D(:)
P_CS_ACC(:) = P_CS_ACC(:) * PPF1D(:)
P_RR_ACCSS(:) = P_RR_ACCSS(:) * PPF1D(:)
P_RR_ACCSG(:) = P_RR_ACCSG(:) * PPF1D(:)
P_RS_ACCRG(:) = P_RS_ACCRG(:) * PPF1D(:)
P_TH_ACC(:)   = (P_RR_ACCSS(:) + P_RR_ACCSG(:)) * (PLSFACT(:)-PLVFACT(:))
!
!
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('LIMA_RAIN_ACCR_SNOW', 1, ZHOOK_HANDLE)
CONTAINS
  FUNCTION GET_XKER_RACCSS(KSIZE,K1,K2) RESULT(RET)
USE YOMHOOK, ONLY:LHOOK, DR_HOOK, JPHOOK
    INTEGER, INTENT(IN) :: KSIZE
    INTEGER, DIMENSION(KSIZE), INTENT(IN) :: K1
    INTEGER, DIMENSION(KSIZE), INTENT(IN) :: K2
    REAL, DIMENSION(KSIZE) :: RET
    !
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
    INTEGER I
    !
    IF (LHOOK) CALL DR_HOOK('GET_XKER_RACCSS', 0, ZHOOK_HANDLE)
    DO I=1,KSIZE
       RET(I) = LIMAM%XKER_RACCSS(MAX(MIN(K1(I),SIZE(LIMAM%XKER_RACCSS,1)),1),MAX(MIN(K2(I),SIZE(LIMAM%XKER_RACCSS,2)),1))
    END DO
  IF (LHOOK) CALL DR_HOOK('GET_XKER_RACCSS', 1, ZHOOK_HANDLE)
END FUNCTION GET_XKER_RACCSS
!
!-------------------------------------------------------------------------------
!
  FUNCTION GET_XKER_N_RACCSS(KSIZE,K1,K2) RESULT(RET)
USE YOMHOOK, ONLY:LHOOK, DR_HOOK, JPHOOK
    INTEGER, INTENT(IN) :: KSIZE
    INTEGER, DIMENSION(KSIZE), INTENT(IN) :: K1
    INTEGER, DIMENSION(KSIZE), INTENT(IN) :: K2
    REAL, DIMENSION(KSIZE) :: RET
    !
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
    INTEGER I
    !
    IF (LHOOK) CALL DR_HOOK('GET_XKER_N_RACCSS', 0, ZHOOK_HANDLE)
    DO I=1,KSIZE
       RET(I) = LIMAM%XKER_N_RACCSS(MAX(MIN(K1(I),SIZE(LIMAM%XKER_N_RACCSS,1)),1),MAX(MIN(K2(I),SIZE(LIMAM%XKER_N_RACCSS,2)),1))
    END DO
  IF (LHOOK) CALL DR_HOOK('GET_XKER_N_RACCSS', 1, ZHOOK_HANDLE)
END FUNCTION GET_XKER_N_RACCSS
!
!-------------------------------------------------------------------------------
!
  FUNCTION GET_XKER_RACCS(KSIZE,K1,K2) RESULT(RET)
USE YOMHOOK, ONLY:LHOOK, DR_HOOK, JPHOOK
    INTEGER, INTENT(IN) :: KSIZE
    INTEGER, DIMENSION(KSIZE), INTENT(IN) :: K1
    INTEGER, DIMENSION(KSIZE), INTENT(IN) :: K2
    REAL, DIMENSION(KSIZE) :: RET
    !
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
    INTEGER I
    !
    IF (LHOOK) CALL DR_HOOK('GET_XKER_RACCS', 0, ZHOOK_HANDLE)
    DO I=1,KSIZE
       RET(I) = LIMAM%XKER_RACCS(MAX(MIN(K1(I),SIZE(LIMAM%XKER_RACCS,1)),1),MAX(MIN(K2(I),SIZE(LIMAM%XKER_RACCS,2)),1))
    END DO
  IF (LHOOK) CALL DR_HOOK('GET_XKER_RACCS', 1, ZHOOK_HANDLE)
END FUNCTION GET_XKER_RACCS
!
!-------------------------------------------------------------------------------
!
  FUNCTION GET_XKER_N_RACCS(KSIZE,K1,K2) RESULT(RET)
USE YOMHOOK, ONLY:LHOOK, DR_HOOK, JPHOOK
    INTEGER, INTENT(IN) :: KSIZE
    INTEGER, DIMENSION(KSIZE), INTENT(IN) :: K1
    INTEGER, DIMENSION(KSIZE), INTENT(IN) :: K2
    REAL, DIMENSION(KSIZE) :: RET
    !
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
    INTEGER I
    !
    IF (LHOOK) CALL DR_HOOK('GET_XKER_N_RACCS', 0, ZHOOK_HANDLE)
    DO I=1,KSIZE
       RET(I) = LIMAM%XKER_N_RACCS(MAX(MIN(K1(I),SIZE(LIMAM%XKER_N_RACCS,1)),1),MAX(MIN(K2(I),SIZE(LIMAM%XKER_N_RACCS,2)),1))
    END DO
  IF (LHOOK) CALL DR_HOOK('GET_XKER_N_RACCS', 1, ZHOOK_HANDLE)
END FUNCTION GET_XKER_N_RACCS
!
!-------------------------------------------------------------------------------
!
  FUNCTION GET_XKER_SACCRG(KSIZE,K1,K2) RESULT(RET)
USE YOMHOOK, ONLY:LHOOK, DR_HOOK, JPHOOK
    INTEGER, INTENT(IN) :: KSIZE
    INTEGER, DIMENSION(KSIZE), INTENT(IN) :: K1
    INTEGER, DIMENSION(KSIZE), INTENT(IN) :: K2
    REAL, DIMENSION(KSIZE) :: RET
    !
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
    INTEGER I
    !
    IF (LHOOK) CALL DR_HOOK('GET_XKER_SACCRG', 0, ZHOOK_HANDLE)
    DO I=1,KSIZE
       RET(I) = LIMAM%XKER_SACCRG(MAX(MIN(K1(I),SIZE(LIMAM%XKER_SACCRG,1)),1),MAX(MIN(K2(I),SIZE(LIMAM%XKER_SACCRG,2)),1))
    END DO
  IF (LHOOK) CALL DR_HOOK('GET_XKER_SACCRG', 1, ZHOOK_HANDLE)
END FUNCTION GET_XKER_SACCRG
!
!-------------------------------------------------------------------------------
!
  FUNCTION GET_XKER_N_SACCRG(KSIZE,K1,K2) RESULT(RET)
USE YOMHOOK, ONLY:LHOOK, DR_HOOK, JPHOOK
    INTEGER, INTENT(IN) :: KSIZE
    INTEGER, DIMENSION(KSIZE), INTENT(IN) :: K1
    INTEGER, DIMENSION(KSIZE), INTENT(IN) :: K2
    REAL, DIMENSION(KSIZE) :: RET
    !
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
    INTEGER I
    !
    IF (LHOOK) CALL DR_HOOK('GET_XKER_N_SACCRG', 0, ZHOOK_HANDLE)
    DO I=1,KSIZE
       RET(I) = LIMAM%XKER_N_SACCRG(MAX(MIN(K1(I),SIZE(LIMAM%XKER_N_SACCRG,1)),1),MAX(MIN(K2(I),SIZE(LIMAM%XKER_N_SACCRG,2)),1))
    END DO
  IF (LHOOK) CALL DR_HOOK('GET_XKER_N_SACCRG', 1, ZHOOK_HANDLE)
END FUNCTION GET_XKER_N_SACCRG
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_RAIN_ACCR_SNOW
END MODULE MODE_LIMA_RAIN_ACCR_SNOW
