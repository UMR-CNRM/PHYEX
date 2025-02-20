!MNH_LIC Copyright 2018-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-------------------------------------------------------------------------------
MODULE MODE_LIMA_DROPLETS_RIMING_SNOW
  IMPLICIT NONE
CONTAINS
!     ###########################################################################################
  SUBROUTINE LIMA_DROPLETS_RIMING_SNOW (CST, LIMAP, LIMAC, LIMAM, KSIZE, PTSTEP, ODCOMPUTE,                               &
                                        PRHODREF, PT,                                           &
                                        PRCT, PCCT, PRST, PCST, PLBDC, PLBDS, PLVFACT, PLSFACT, &
!++cb++
!                                        P_TH_RIM, P_RC_RIM, P_CC_RIM, P_RS_RIM, P_CS_RIM, P_RG_RIM,       &
                                        P_TH_RIM, P_CC_RIM, P_CS_RIM,                           &
                                        P_RC_RIMSS, P_RC_RIMSG, P_RS_RIMCG,                     &
                                        P_RI_HMS, P_CI_HMS, P_RS_HMS                            )
!     ###########################################################################################
!
!!    PURPOSE
!!    -------
!!      Compute the cloud droplets riming of the aggregates rate, and the associated
!!    Hallett-Mossop ice production rate          
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
!  C. Barthe      06/2022: modify the microphysics terms to save to simplify the merging with the electrification scheme
!                          (same terms as in ICE3)
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
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PT   ! 
!
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRCT    ! Cloud water mr at t
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PCCT    ! Cloud water C. at t
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRST    ! Snow mr at t
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PCST    ! Snow C. at t
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLBDC   ! 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLBDS   ! 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLVFACT ! 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLSFACT ! 
!
!++cb++
!REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_RC_RIM
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_CC_RIM
!REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_RS_RIM
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_CS_RIM
!REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_RG_RIM
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_RC_RIMSS
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_RC_RIMSG
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_RS_RIMCG
!--cb--
!
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_TH_RIM
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_RI_HMS
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_CI_HMS
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_RS_HMS
!
!*       0.2   Declarations of local variables :
!
REAL,    DIMENSION(SIZE(PRCT))  :: ZZW1, ZZW2, ZZW3, ZZW4, ZZW5
REAL,    DIMENSION(SIZE(PRCT))  :: Z_RC_RIM, Z_RS_RIM, Z_RG_RIM    !++cb--
!
INTEGER, DIMENSION(SIZE(PRCT))  :: IVEC2              ! Vector of indices
REAL,    DIMENSION(SIZE(PRCT))  :: ZVEC1,ZVEC2,ZVEC1W ! Work vectors
INTEGER                         :: II
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!
IF (LHOOK) CALL DR_HOOK('LIMA_DROPLETS_RIMING_SNOW', 0, ZHOOK_HANDLE)
DO II = 1, SIZE(PRCT)
!
!*       Cloud droplet riming of the aggregates  
!        --------------------------------------
!
  IF ( PRCT(II)>LIMAP%XRTMIN(2) .AND. PRST(II)>LIMAP%XRTMIN(5) .AND. PT(II)<CST%XTT .AND. &
       PCCT(II)>LIMAP%XCTMIN(2) .AND. PCST(II)>LIMAP%XCTMIN(5) .AND. ODCOMPUTE(II) ) THEN
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
                     LIMAM%XRIMINTP1 * LOG( ZVEC1(II) ) + LIMAM%XRIMINTP2 ) )
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
    Z_RC_RIM(II) = - LIMAM%XCRIMSS  * PRCT(II) * PCST(II)*(1+(LIMAC%XFVELOS/PLBDS(II))&
         **LIMAP%XALPHAS)**(-LIMAP%XNUS+LIMAM%XEXCRIMSS/LIMAP%XALPHAS) &
                                * PRHODREF(II)**(-LIMAP%XCEXVT+1) * PLBDS(II)**LIMAM%XEXCRIMSS
    P_CC_RIM(II) = Z_RC_RIM(II) * (PCCT(II) / PRCT(II)) ! Lambda_c**3
    !
    ! Cloud droplets collected on small aggregates add to snow
!      P_RS_RIM(II) = - P_RC_RIM(II) * ZZW1(II)
    Z_RS_RIM(II)   = -Z_RC_RIM(II) * ZZW1(II)
    P_RC_RIMSS(II) = Z_RC_RIM(II) * ZZW1(II)  ! < 0, loss of mass for rc
    !
    ! Cloud droplets collected on large aggregates add to graupel
    Z_RG_RIM(II)   = -Z_RC_RIM(II) - Z_RS_RIM(II)
    P_RC_RIMSG(II) = Z_RC_RIM(II) - P_RC_RIMSS(II) ! < 0, loss of mass for rc
    !
    IF (LIMAP%LMURAKAMI) THEN
    ! Graupel formation based on Murakami
      ZVEC1(II) = LIMAM%XGAMINC_RIM4( IVEC2(II)+1 )* ZVEC2(II)      &
                - LIMAM%XGAMINC_RIM4( IVEC2(II)   )*(ZVEC2(II) - 1.0)
      ZZW5(II) = ZVEC1(II)
      ZZW3(II) = LIMAM%XSRIMCG * PRHODREF(II) * PCST(II) * PLBDS(II)**LIMAM%XEXSRIMCG * (1.0 - ZZW2(II))!/(PTSTEP*PRHODREF(II))
      ZZW3(II) = Z_RG_RIM(II)*ZZW3(II)/ &
                    MAX(1.E-10, & !-20
                        LIMAM%XSRIMCG3*LIMAM%XSRIMCG2*PCST(II)*PRHODREF(II)*PLBDS(II)**(LIMAM%XEXSRIMCG2)*(1.-ZZW5(II))- &
                        LIMAM%XSRIMCG3*ZZW3(II))
    ELSE
    ! Large aggregates collecting droplets add to graupel (instant process ???)
      ZZW3(II) = PRST(II)*(1.0 - ZZW2(II))/PTSTEP
    END IF
    !
    P_RS_RIMCG(II) = ZZW3(II)
    P_CS_RIM(II) = -ZZW3(II) * PCST(II)/PRST(II)
    !
    P_TH_RIM(II) = - Z_RC_RIM(II)*(PLSFACT(II)-PLVFACT(II))
  ELSE
    P_TH_RIM(II) = 0.
    P_RC_RIMSS(II) = 0.
    P_RC_RIMSG(II) = 0.
    P_RS_RIMCG(II) = 0.
    Z_RC_RIM(II) = 0.
    P_CC_RIM(II) = 0.
    Z_RS_RIM(II) = 0.
    P_CS_RIM(II) = 0.
    Z_RG_RIM(II) = 0.
  END IF
!
!*       Hallett-Mossop ice production (HMS)  
!        -----------------------------------
!
  IF ( PRST(II)>LIMAP%XRTMIN(5) .AND. PRCT(II)>LIMAP%XRTMIN(2) .AND. PT(II)<LIMAM%XHMTMAX .AND. PT(II)>LIMAM%XHMTMIN .AND. &
       PCST(II)>LIMAP%XCTMIN(5) .AND. PCCT(II)>LIMAP%XCTMIN(2) .AND. ODCOMPUTE(II) ) THEN
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
  ELSE
    P_RI_HMS(II) = 0.
    P_CI_HMS(II) = 0.
    P_RS_HMS(II) = 0.
  END IF
END DO
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('LIMA_DROPLETS_RIMING_SNOW', 1, ZHOOK_HANDLE)
END SUBROUTINE LIMA_DROPLETS_RIMING_SNOW
END MODULE MODE_LIMA_DROPLETS_RIMING_SNOW
