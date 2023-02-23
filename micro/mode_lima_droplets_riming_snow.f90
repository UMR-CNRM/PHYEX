!MNH_LIC Copyright 2018-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-------------------------------------------------------------------------------
MODULE MODE_LIMA_DROPLETS_RIMING_SNOW
  IMPLICIT NONE
CONTAINS
!     #########################################################################################
  SUBROUTINE LIMA_DROPLETS_RIMING_SNOW (PTSTEP, LDCOMPUTE,                                                &
                                        PRHODREF, PT,                                                     &
                                        PRCT, PCCT, PRST, PCST, PLBDC, PLBDS, PLVFACT, PLSFACT,           &
                                        P_TH_RIM, P_RC_RIM, P_CC_RIM, P_RS_RIM, P_CS_RIM, P_RG_RIM,       &
                                        P_RI_HMS, P_CI_HMS, P_RS_HMS                                      )
!     #########################################################################################
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
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST,              ONLY : XTT
USE MODD_PARAM_LIMA,       ONLY : XRTMIN, XCTMIN, XCEXVT, XNUS, XALPHAS, LMURAKAMI
USE MODD_PARAM_LIMA_MIXED, ONLY : NGAMINC, XRIMINTP1, XRIMINTP2, XGAMINC_RIM1, XGAMINC_RIM2, XGAMINC_RIM4, &
                                  XCRIMSS, XEXCRIMSS, XSRIMCG, XEXSRIMCG, XSRIMCG2, XSRIMCG3, XEXSRIMCG2, &
                                  XHMLINTP1, XHMLINTP2, XGAMINC_HMC, XHM_FACTS, XHMTMIN, XHMTMAX
USE MODD_PARAM_LIMA_COLD,  ONLY : XMNU0, XBS, XFVELOS
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
REAL,                 INTENT(IN)    :: PTSTEP 
LOGICAL, DIMENSION(:),INTENT(IN)    :: LDCOMPUTE
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRHODREF    ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PT   ! 
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRCT    ! Cloud water mr at t
REAL, DIMENSION(:),   INTENT(IN)    :: PCCT    ! Cloud water C. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PRST    ! Snow mr at t
REAL, DIMENSION(:),   INTENT(IN)    :: PCST    ! Snow C. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDC   ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDS   ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PLVFACT ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PLSFACT ! 
!
REAL, DIMENSION(:),   INTENT(OUT)   :: P_RC_RIM
REAL, DIMENSION(:),   INTENT(OUT)   :: P_CC_RIM
REAL, DIMENSION(:),   INTENT(OUT)   :: P_RS_RIM
REAL, DIMENSION(:),   INTENT(OUT)   :: P_CS_RIM
REAL, DIMENSION(:),   INTENT(OUT)   :: P_RG_RIM
!
REAL, DIMENSION(:),   INTENT(OUT)   :: P_TH_RIM
REAL, DIMENSION(:),   INTENT(OUT)   :: P_RI_HMS
REAL, DIMENSION(:),   INTENT(OUT)   :: P_CI_HMS
REAL, DIMENSION(:),   INTENT(OUT)   :: P_RS_HMS
!
!*       0.2   Declarations of local variables :
!
REAL,    DIMENSION(SIZE(PRCT))  :: ZZW1, ZZW2, ZZW3, ZZW4, ZZW5
!
INTEGER, DIMENSION(SIZE(PRCT))  :: IVEC1,IVEC2        ! Vectors of indices
REAL,    DIMENSION(SIZE(PRCT))  :: ZVEC1,ZVEC2,ZVEC1W ! Work vectors
INTEGER                         :: JI
!
!-------------------------------------------------------------------------------
!
!
!
DO JI = 1, SIZE(PRCT)
!
!*       Cloud droplet riming of the aggregates  
!        --------------------------------------
!
   IF ( PRCT(JI)>XRTMIN(2) .AND. PRST(JI)>XRTMIN(5) .AND. PT(JI)<XTT .AND. &
        PCCT(JI)>XCTMIN(2) .AND. PCST(JI)>XCTMIN(5) .AND. LDCOMPUTE(JI) ) THEN
!
      ZVEC1(JI) = PLBDS(JI)
      ZVEC1W(JI)= ( XFVELOS**XALPHAS + PLBDS(JI)**XALPHAS ) ** (1./XALPHAS) ! modified equivalent lambda
!
!        2.     perform the linear interpolation of the normalized
!               "2+XDS"-moment of the incomplete gamma function using the modified equivalent lambda
!
      ZVEC2(JI) = MAX( 1.0001, MIN( REAL(NGAMINC)-0.0001,           &
                       XRIMINTP1 * LOG( ZVEC1W(JI) ) + XRIMINTP2 ) )
      IVEC2(JI) = INT( ZVEC2(JI) )
      ZVEC2(JI) = ZVEC2(JI) - REAL( IVEC2(JI) )
!
      ZZW1(JI)  =  XGAMINC_RIM1( IVEC2(JI)+1 )* ZVEC2(JI)      &
                 - XGAMINC_RIM1( IVEC2(JI)   )*(ZVEC2(JI) - 1.0)
!
!        3.     perform the linear interpolation of the normalized
!               "XBS"-moment of the incomplete gamma function
!
      ZVEC2(JI) = MAX( 1.0001, MIN( REAL(NGAMINC)-0.0001,           &
                       XRIMINTP1 * LOG( ZVEC1(JI) ) + XRIMINTP2 ) )
      IVEC2(JI) = INT( ZVEC2(JI) )
      ZVEC2(JI) = ZVEC2(JI) - REAL( IVEC2(JI) )
!
      ZZW2(JI)  = XGAMINC_RIM2( IVEC2(JI)+1 )* ZVEC2(JI)      &
                - XGAMINC_RIM2( IVEC2(JI)   )*(ZVEC2(JI) - 1.0)
!
!        4.     riming
!
   ! Cloud droplets collected
      P_RC_RIM(JI) = - XCRIMSS  * PRCT(JI) * PCST(JI)*(1+(XFVELOS/PLBDS(JI))**XALPHAS)**(-XNUS+XEXCRIMSS/XALPHAS) &
                                * PRHODREF(JI)**(-XCEXVT+1) * PLBDS(JI)**XEXCRIMSS
      P_CC_RIM(JI) = P_RC_RIM(JI) * PCCT(JI)/PRCT(JI) ! Lambda_c**3
   !
   ! Cloud droplets collected on small aggregates add to snow
      P_RS_RIM(JI) = - P_RC_RIM(JI) * ZZW1(JI)
   !
   ! Cloud droplets collected on large aggregates add to graupel
      P_RG_RIM(JI) = - P_RC_RIM(JI) - P_RS_RIM(JI) 
   !
      IF (LMURAKAMI) THEN
   ! Graupel formation based on Murakami
         ZVEC1(JI) = XGAMINC_RIM4( IVEC2(JI)+1 )* ZVEC2(JI)      &
                   - XGAMINC_RIM4( IVEC2(JI)   )*(ZVEC2(JI) - 1.0)
         ZZW5(JI) = ZVEC1(JI)
         ZZW3(JI) = XSRIMCG * PRHODREF(JI) * PCST(JI) * PLBDS(JI)**XEXSRIMCG * (1.0 - ZZW2(JI))!/(PTSTEP*PRHODREF(JI))
         ZZW3(JI) = P_RG_RIM(JI)*ZZW3(JI)/ &
                       MAX(1.E-10, & !-20
                           XSRIMCG3*XSRIMCG2*PCST(JI)*PRHODREF(JI)*PLBDS(JI)**(XEXSRIMCG2)*(1.-ZZW5(JI))- &
                           XSRIMCG3*ZZW3(JI))
      ELSE
      ! Large aggregates collecting droplets add to graupel (instant process ???)
         ZZW3(JI) = PRST(JI)*(1.0 - ZZW2(JI))/PTSTEP
      END IF
      P_RS_RIM(JI) = P_RS_RIM(JI) - ZZW3(JI) 
      P_CS_RIM(JI) = -ZZW3(JI) * PCST(JI)/PRST(JI)
      P_RG_RIM(JI) = P_RG_RIM(JI) + ZZW3(JI) 
   !
      P_TH_RIM(JI) = - P_RC_RIM(JI)*(PLSFACT(JI)-PLVFACT(JI))
   ELSE
      P_TH_RIM(JI) = 0.
      P_RC_RIM(JI) = 0.
      P_CC_RIM(JI) = 0.
      P_RS_RIM(JI) = 0.
      P_CS_RIM(JI) = 0.
      P_RG_RIM(JI) = 0.
   END IF
!
!*       Hallett-Mossop ice production (HMS)  
!        -----------------------------------
!
   IF ( PRST(JI)>XRTMIN(5) .AND. PRCT(JI)>XRTMIN(2) .AND. PT(JI)<XHMTMAX .AND. PT(JI)>XHMTMIN .AND. &
        PCST(JI)>XCTMIN(5) .AND. PCCT(JI)>XCTMIN(2) .AND. LDCOMPUTE(JI) ) THEN
!
      ZVEC1(JI) = PLBDC(JI)
      ZVEC2(JI) = MAX( 1.0001, MIN( REAL(NGAMINC)-0.0001,           &
                         XHMLINTP1 * LOG( ZVEC1(JI) ) + XHMLINTP2 ) )
      IVEC2(JI) = INT( ZVEC2(JI) )
      ZVEC2(JI) = ZVEC2(JI) - REAL( IVEC2(JI) )
      ZVEC1(JI) =  XGAMINC_HMC( IVEC2(JI)+1 )* ZVEC2(JI)      &
                 - XGAMINC_HMC( IVEC2(JI)   )*(ZVEC2(JI) - 1.0)
      ZZW4(JI) = ZVEC1(JI) ! Large droplets
!
      IF ( ZZW4(JI)<0.99 ) THEN
         P_CI_HMS(JI) = - P_RC_RIM(JI) * (PCCT(JI)/PRCT(JI)) * (1.0-ZZW4(JI)) * XHM_FACTS * &
                          MAX( 0.0, MIN( (PT(JI)-XHMTMIN)/3.0,(XHMTMAX-PT(JI))/2.0 ) ) ! CCHMSI
!
         P_RI_HMS(JI) = P_CI_HMS(JI) * XMNU0                                     ! RCHMSI
         P_RS_HMS(JI) = - P_RI_HMS(JI)
      ELSE
         P_RI_HMS(JI) = 0.
         P_CI_HMS(JI) = 0.
         P_RS_HMS(JI) = 0.
      END IF
   ELSE
      P_RI_HMS(JI) = 0.
      P_CI_HMS(JI) = 0.
      P_RS_HMS(JI) = 0.
   END IF
END DO
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_DROPLETS_RIMING_SNOW
END MODULE MODE_LIMA_DROPLETS_RIMING_SNOW
