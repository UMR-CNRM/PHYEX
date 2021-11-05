!MNH_LIC Copyright 2018-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-------------------------------------------------------------------------------
!      #################################
       MODULE MODI_LIMA_DROPLETS_RIMING_SNOW
!      #################################
!
INTERFACE
   SUBROUTINE LIMA_DROPLETS_RIMING_SNOW (PTSTEP,  LDCOMPUTE,                               &
                                         PRHODREF, PT,                                     &
                                         PRCT, PCCT, PRST, PLBDC, PLBDS, PLVFACT, PLSFACT, &
                                         P_TH_RIM, P_RC_RIM, P_CC_RIM, P_RS_RIM, P_RG_RIM, &
                                         P_RI_HMS, P_CI_HMS, P_RS_HMS                      )
!
REAL,                 INTENT(IN)    :: PTSTEP 
LOGICAL, DIMENSION(:),INTENT(IN)    :: LDCOMPUTE
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRHODREF    ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PT   ! 
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRCT    ! Cloud water C. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PCCT    ! Cloud water C. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PRST    ! Cloud water C. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDC   ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDS   ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PLVFACT ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PLSFACT ! 
!
REAL, DIMENSION(:),   INTENT(OUT)   :: P_TH_RIM
REAL, DIMENSION(:),   INTENT(OUT)   :: P_RC_RIM
REAL, DIMENSION(:),   INTENT(OUT)   :: P_CC_RIM
REAL, DIMENSION(:),   INTENT(OUT)   :: P_RS_RIM
REAL, DIMENSION(:),   INTENT(OUT)   :: P_RG_RIM
!
REAL, DIMENSION(:),   INTENT(OUT)   :: P_RI_HMS
REAL, DIMENSION(:),   INTENT(OUT)   :: P_CI_HMS
REAL, DIMENSION(:),   INTENT(OUT)   :: P_RS_HMS
!
END SUBROUTINE LIMA_DROPLETS_RIMING_SNOW
END INTERFACE
END MODULE MODI_LIMA_DROPLETS_RIMING_SNOW
!
!     #########################################################################################
      SUBROUTINE LIMA_DROPLETS_RIMING_SNOW (PTSTEP, LDCOMPUTE,                                &
                                            PRHODREF, PT,                                     &
                                            PRCT, PCCT, PRST, PLBDC, PLBDS, PLVFACT, PLSFACT, &
                                            P_TH_RIM, P_RC_RIM, P_CC_RIM, P_RS_RIM, P_RG_RIM, &
                                            P_RI_HMS, P_CI_HMS, P_RS_HMS                      )
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
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST,              ONLY : XTT
USE MODD_PARAM_LIMA,       ONLY : XRTMIN, XCEXVT
USE MODD_PARAM_LIMA_MIXED, ONLY : NGAMINC, XRIMINTP1, XRIMINTP2, XGAMINC_RIM1, XGAMINC_RIM2, &
                                  XCRIMSS, XEXCRIMSS, XSRIMCG, XEXSRIMCG, &
                                  XHMLINTP1, XHMLINTP2, XGAMINC_HMC, XHM_FACTS, XHMTMIN, XHMTMAX
USE MODD_PARAM_LIMA_COLD,  ONLY : XMNU0
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
REAL, DIMENSION(:),   INTENT(IN)    :: PRCT    ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PCCT    ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PRST    ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDC   ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDS   ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PLVFACT ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PLSFACT ! 
!
REAL, DIMENSION(:),   INTENT(OUT)   :: P_RC_RIM
REAL, DIMENSION(:),   INTENT(OUT)   :: P_CC_RIM
REAL, DIMENSION(:),   INTENT(OUT)   :: P_RS_RIM
REAL, DIMENSION(:),   INTENT(OUT)   :: P_RG_RIM
!
REAL, DIMENSION(:),   INTENT(OUT)   :: P_TH_RIM
REAL, DIMENSION(:),   INTENT(OUT)   :: P_RI_HMS
REAL, DIMENSION(:),   INTENT(OUT)   :: P_CI_HMS
REAL, DIMENSION(:),   INTENT(OUT)   :: P_RS_HMS
!
!*       0.2   Declarations of local variables :
!
LOGICAL, DIMENSION(SIZE(PRCT))  :: GRIM
!
REAL,    DIMENSION(SIZE(PRCT))  :: ZZW1, ZZW2, ZZW3, ZZW4
!
INTEGER, DIMENSION(SIZE(PRCT))  :: IVEC1,IVEC2      ! Vectors of indices
REAL,    DIMENSION(SIZE(PRCT))  :: ZVEC1,ZVEC2      ! Work vectors
!
!-------------------------------------------------------------------------------
!
!
P_TH_RIM(:) = 0.
P_RC_RIM(:) = 0.
P_CC_RIM(:) = 0.
P_RS_RIM(:) = 0.
P_RG_RIM(:) = 0.
!
P_RI_HMS(:) = 0.
P_CI_HMS(:) = 0.
P_RS_HMS(:) = 0.
!
ZZW1(:) = 0.
ZZW2(:) = 0.
ZZW3(:) = 0.
ZZW4(:) = 0.
!
!*       Cloud droplet riming of the aggregates  
!        --------------------------------------
!
!
GRIM(:) = .False.
GRIM(:) = (PRCT(:)>XRTMIN(2)) .AND. (PRST(:)>XRTMIN(5)) .AND. (PT(:)<XTT) .AND. LDCOMPUTE(:)
!
WHERE( GRIM )
!
   ZVEC1(:) = PLBDS(:)
!
!        1.     find the next lower indice for the ZLBDAS in the geometrical
!               set of Lbda_s used to tabulate some moments of the incomplete 
!               gamma function
!
   ZVEC2(:) = MAX( 1.0001, MIN( REAL(NGAMINC)-0.0001,           &
                       XRIMINTP1 * LOG( ZVEC1(:) ) + XRIMINTP2 ) )
   IVEC2(:) = INT( ZVEC2(:) )
   ZVEC2(:) = ZVEC2(:) - REAL( IVEC2(:) )
!
!        2.     perform the linear interpolation of the normalized
!               "2+XDS"-moment of the incomplete gamma function
!
   ZVEC1(:) =   XGAMINC_RIM1( IVEC2(:)+1 )* ZVEC2(:)      &
                    - XGAMINC_RIM1( IVEC2(:)   )*(ZVEC2(:) - 1.0)
   ZZW1(:) = ZVEC1(:)
!
!        3.     perform the linear interpolation of the normalized
!               "XBS"-moment of the incomplete gamma function
!
   ZVEC1(:) =  XGAMINC_RIM2( IVEC2(:)+1 )* ZVEC2(:)      &
                   - XGAMINC_RIM2( IVEC2(:)   )*(ZVEC2(:) - 1.0)
   ZZW2(:) = ZVEC1(:)
!
!        4.     riming
!
   ! Cloud droplets collected
   P_RC_RIM(:) = - XCRIMSS * PRCT(:) * PLBDS(:)**XEXCRIMSS * PRHODREF(:)**(-XCEXVT)
   P_CC_RIM(:) = P_RC_RIM(:) *(PCCT(:)/PRCT(:)) ! Lambda_c**3
   !
   ! Cloud droplets collected on small aggregates add to snow
   P_RS_RIM(:) = - P_RC_RIM(:) * ZZW1(:)
   !
   ! Cloud droplets collected on large aggregates add to graupel
   P_RG_RIM(:) = - P_RC_RIM(:) - P_RS_RIM(:) 
   !
   ! Large aggregates collecting droplets add to graupel (instant process ???)
   ZZW3(:) = XSRIMCG * PLBDS(:)**XEXSRIMCG * (1.0 - ZZW2(:))/(PTSTEP*PRHODREF(:))
   P_RS_RIM(:) = P_RS_RIM(:) - ZZW3(:) 
   P_RG_RIM(:) = P_RG_RIM(:) + ZZW3(:) 
   !
   P_TH_RIM(:) = - P_RC_RIM(:)*(PLSFACT(:)-PLVFACT(:))
END WHERE
!
!
!*       Hallett-Mossop ice production (HMS)  
!        -----------------------------------
!
!
GRIM(:) = .False.
GRIM(:) = (PT(:)<XHMTMAX)     .AND. (PT(:)>XHMTMIN)     .AND. &
              (PRST(:)>XRTMIN(5)) .AND. (PRCT(:)>XRTMIN(2)) .AND. &
               LDCOMPUTE(:)
!
WHERE ( GRIM )
!
   ZVEC1(:) = PLBDC(:)
   ZVEC2(:) = MAX( 1.0001, MIN( REAL(NGAMINC)-0.0001,           &
                         XHMLINTP1 * LOG( ZVEC1(:) ) + XHMLINTP2 ) )
   IVEC2(:) = INT( ZVEC2(:) )
   ZVEC2(:) = ZVEC2(:) - REAL( IVEC2(:) )
   ZVEC1(:) =   XGAMINC_HMC( IVEC2(:)+1 )* ZVEC2(:)      &
                    - XGAMINC_HMC( IVEC2(:)   )*(ZVEC2(:) - 1.0)
   ZZW4(:) = ZVEC1(:) ! Large droplets
!
   WHERE ( ZZW4(:)<0.99 )
      P_CI_HMS(:) = - P_RC_RIM(:) * (PCCT(:)/PRCT(:)) * (1.0-ZZW4(:)) * XHM_FACTS * &
                          MAX( 0.0, MIN( (PT(:)-XHMTMIN)/3.0,(XHMTMAX-PT(:))/2.0 ) ) ! CCHMSI
!
      P_RI_HMS(:) = P_CI_HMS(:) * XMNU0                                     ! RCHMSI
      P_RS_HMS(:) = - P_RI_HMS(:) 
   END WHERE

END WHERE
!
!
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_DROPLETS_RIMING_SNOW
