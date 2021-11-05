!MNH_LIC Copyright 2018-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!      #################################
       MODULE MODI_LIMA_RAIN_ACCR_SNOW
!      #################################
!
INTERFACE
   SUBROUTINE LIMA_RAIN_ACCR_SNOW (PTSTEP, LDCOMPUTE,                                &
                                   PRHODREF, PT,                                     &
                                   PRRT, PCRT, PRST, PLBDR, PLBDS, PLVFACT, PLSFACT, &
                                   P_TH_ACC, P_RR_ACC, P_CR_ACC, P_RS_ACC, P_RG_ACC  )
!
REAL,                 INTENT(IN)    :: PTSTEP 
LOGICAL, DIMENSION(:),INTENT(IN)    :: LDCOMPUTE
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRHODREF    ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PT   ! 
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRRT    ! Cloud water C. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PCRT    ! Cloud water C. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PRST    ! Cloud water C. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDR   ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDS   ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PLVFACT ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PLSFACT ! 
!
REAL, DIMENSION(:),   INTENT(OUT)   :: P_TH_ACC
REAL, DIMENSION(:),   INTENT(OUT)   :: P_RR_ACC
REAL, DIMENSION(:),   INTENT(OUT)   :: P_CR_ACC
REAL, DIMENSION(:),   INTENT(OUT)   :: P_RS_ACC
REAL, DIMENSION(:),   INTENT(OUT)   :: P_RG_ACC
!
END SUBROUTINE LIMA_RAIN_ACCR_SNOW
END INTERFACE
END MODULE MODI_LIMA_RAIN_ACCR_SNOW
!
!     ###################################################################################
      SUBROUTINE LIMA_RAIN_ACCR_SNOW (PTSTEP, LDCOMPUTE,                                &
                                      PRHODREF, PT,                                     &
                                      PRRT, PCRT, PRST, PLBDR, PLBDS, PLVFACT, PLSFACT, &
                                      P_TH_ACC, P_RR_ACC, P_CR_ACC, P_RS_ACC, P_RG_ACC  )
!     ###################################################################################
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
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST,              ONLY : XTT
USE MODD_PARAM_LIMA,       ONLY : XRTMIN, XCEXVT
USE MODD_PARAM_LIMA_COLD,  ONLY : XBS, XCXS
USE MODD_PARAM_LIMA_MIXED, ONLY : NACCLBDAS, XACCINTP1S, XACCINTP2S,         &
                                  NACCLBDAR, XACCINTP1R, XACCINTP2R,         &
                                  XKER_RACCSS, XKER_RACCS, XKER_SACCRG,      &
                                  XFRACCSS, XLBRACCS1, XLBRACCS2, XLBRACCS3, &
                                  XFSACCRG, XLBSACCR1, XLBSACCR2, XLBSACCR3
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
REAL, DIMENSION(:),   INTENT(IN)    :: PRRT    ! Cloud water C. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PCRT    ! Cloud water C. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PRST    ! Cloud water C. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDR   ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDS   ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PLVFACT ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PLSFACT ! 
!
REAL, DIMENSION(:),   INTENT(OUT)   :: P_TH_ACC
REAL, DIMENSION(:),   INTENT(OUT)   :: P_RR_ACC
REAL, DIMENSION(:),   INTENT(OUT)   :: P_CR_ACC
REAL, DIMENSION(:),   INTENT(OUT)   :: P_RS_ACC
REAL, DIMENSION(:),   INTENT(OUT)   :: P_RG_ACC
!
!*       0.2   Declarations of local variables :
!
LOGICAL, DIMENSION(SIZE(PRRT))  :: GACC
!
REAL,    DIMENSION(SIZE(PRRT))  :: Z1, Z2, Z3, Z4
REAL,    DIMENSION(SIZE(PRRT))  :: ZZW1, ZZW2, ZZW3, ZZW4, ZZW5
!
INTEGER, DIMENSION(SIZE(PRRT))  :: IVEC1,IVEC2       ! Vectors of indices
REAL,    DIMENSION(SIZE(PRRT))  :: ZVEC1,ZVEC2,ZVEC3 ! Work vectors
!
!-------------------------------------------------------------------------------
!
!
P_TH_ACC(:) = 0.
P_RR_ACC(:) = 0.
P_CR_ACC(:) = 0.
P_RS_ACC(:) = 0.
P_RG_ACC(:) = 0.
!
ZZW1(:) = 0.
ZZW2(:) = 0.
ZZW3(:) = 0.
ZZW4(:) = 0.
ZZW5(:) = 0.
!
IVEC1(:) = 0
IVEC2(:) = 0
ZVEC1(:) = 0.
ZVEC2(:) = 0.
ZVEC3(:) = 0.
!
!*       Cloud droplet riming of the aggregates  
!        -------------------------------------------
!
!
GACC(:) = .False.
GACC(:) = (PRRT(:)>XRTMIN(3)) .AND. (PRST(:)>XRTMIN(5)) .AND. (PT(:)<XTT) .AND. LDCOMPUTE(:)
!
WHERE( GACC )
!
!        1.3.1  select the (ZLBDAS,ZLBDAR) couplet
   !
   ZVEC1(:) = MAX(MIN(PLBDS(:),5.E5),5.E1)
   ZVEC2(:) = PLBDR(:)
!
!        1.3.2  find the next lower indice for the ZLBDAS and for the ZLBDAR
!               in the geometrical set of (Lbda_s,Lbda_r) couplet use to
!               tabulate the RACCSS-kernel
!
   ZVEC1(:) = MAX( 1.0001, MIN( REAL(NACCLBDAS)-0.0001,           &
                         XACCINTP1S * LOG( ZVEC1(:) ) + XACCINTP2S ) )
   IVEC1(:) = INT( ZVEC1(:) )
   ZVEC1(:) = ZVEC1(:) - REAL( IVEC1(:) )
!
   ZVEC2(:) = MAX( 1.0001, MIN( REAL(NACCLBDAR)-0.0001,           &
                         XACCINTP1R * LOG( ZVEC2(:) ) + XACCINTP2R ) )
   IVEC2(:) = INT( ZVEC2(:) )
   ZVEC2(:) = ZVEC2(:) - REAL( IVEC2(:) )
!
!        1.3.3  perform the bilinear interpolation of the normalized
!               RACCSS-kernel : for small rain drops transformed into snow
   !
   Z1(:) = GET_XKER_RACCSS(IVEC1(:)+1,IVEC2(:)+1)
   Z2(:) = GET_XKER_RACCSS(IVEC1(:)+1,IVEC2(:)  )
   Z3(:) = GET_XKER_RACCSS(IVEC1(:)  ,IVEC2(:)+1)
   Z4(:) = GET_XKER_RACCSS(IVEC1(:)  ,IVEC2(:)  )
      ZVEC3(:) =  (   Z1(:)* ZVEC2(:)          &
                    - Z2(:)*(ZVEC2(:) - 1.0) ) &
                				 	     *  ZVEC1(:)    &
                 - (  Z3(:)* ZVEC2(:)          &
                    - Z4(:)*(ZVEC2(:) - 1.0) ) &
  	                    			             * (ZVEC1(:) - 1.0)
   ZZW1(:) = ZVEC3(:)
!
!        1.3.4b perform the bilinear interpolation of the normalized
!               RACCS-kernel : total frozen rain drops
!
   Z1(:) = GET_XKER_RACCS(IVEC1(:)+1,IVEC2(:)+1)
   Z2(:) = GET_XKER_RACCS(IVEC1(:)+1,IVEC2(:)  )
   Z3(:) = GET_XKER_RACCS(IVEC1(:)  ,IVEC2(:)+1)
   Z4(:) = GET_XKER_RACCS(IVEC1(:)  ,IVEC2(:)  )
      ZVEC3(:) =  (    Z1(:)* ZVEC2(:)          &
                    -  Z2(:)*(ZVEC2(:) - 1.0) ) &
                                                           *  ZVEC1(:)      &
                 - (   Z3(:)* ZVEC2(:)          &
                    -  Z4(:)*(ZVEC2(:) - 1.0) ) &
                                                           * (ZVEC1(:) - 1.0)
   ZZW2(:) = ZVEC3(:)
!
! Correction of ZZW1 to ensure that ZZW1 <= ZZW2
! ie                coll. of small drops <= coll. of all drops
!
   ZZW1(:) = MIN(ZZW1(:),ZZW2(:))
!
!        1.3.5  perform the bilinear interpolation of the normalized
!               SACCRG-kernel : snow transformed into graupel
!
   Z1(:) = GET_XKER_SACCRG(IVEC2(:)+1,IVEC1(:)+1)
   Z2(:) = GET_XKER_SACCRG(IVEC2(:)+1,IVEC1(:)  )
   Z3(:) = GET_XKER_SACCRG(IVEC2(:)  ,IVEC1(:)+1)
   Z4(:) = GET_XKER_SACCRG(IVEC2(:)  ,IVEC1(:)  )
      ZVEC3(:) =  (   Z1(:)* ZVEC1(:)          &
                    - Z2(:)*(ZVEC1(:) - 1.0) ) &
      			 	                             *  ZVEC2(:)    &
                 - (  Z3(:)* ZVEC1(:)          &
                    - Z4(:)*(ZVEC1(:) - 1.0) ) &
			                                     * (ZVEC2(:) - 1.0)
   ZZW3(:) = ZVEC3(:)
!
!        1.3.4  raindrop accretion on the small sized aggregates
!      
! BVIE manque PCRT ???????????????????????????????????
!      ZZW4(:) =                                            & !! coef of RRACCS and RRACCS
   ZZW4(:) = PCRT(:)                                     & !! coef of RRACCS and RRACCS
         *  XFRACCSS *( PLBDS(:)**XCXS )*( PRHODREF(:)**(-XCEXVT-1.) ) &
         *( XLBRACCS1/( PLBDS(:)**2 )                +                     &
            XLBRACCS2/( PLBDS(:) * PLBDR(:)    ) +                     &
            XLBRACCS3/(                PLBDR(:)**2 ) ) / PLBDR(:)**3

!      ZRRS(:) = ZRRS(:) - ZZW1(:,4)
!      ZRSS(:) = ZRSS(:) + ZZW1(:,4)
!      ZTHS(:) = ZTHS(:) + ZZW1(:,4)*(ZLSFACT(:)-ZLVFACT(:)) ! f(L_f*(RRACCSS))
!
!      ZCRS(:) = MAX( ZCRS(:)-ZZW1(:,4)*(ZCRT(:)/ZRRT(:)),0.0 ) ! Lambda_r**3 
!
!        1.3.6  raindrop accretion-conversion of the large sized aggregates
!               into graupeln
!
   ZZW5(:) = XFSACCRG*ZZW3(:) *                             & ! RSACCRG
            ( PLBDS(:)**(XCXS-XBS) )*( PRHODREF(:)**(-XCEXVT-1.) ) &
           *( XLBSACCR1/((PLBDR(:)**2)               ) +           &
              XLBSACCR2/( PLBDR(:)    * PLBDS(:) ) +           &
              XLBSACCR3/(                  (PLBDS(:)**2)) )
      !
!      P_RR_ACC(:) = - ZZW4(:) * ZZW1(:)           ! RRACCSS
!      P_CR_ACC(:) = P_RR_ACC(:) * PCRT(:)/PRRT(:) ! Lambda_r**3 
!      P_RS_ACC(:) = - P_RR_ACC(:) 
      !
!      P_RR_ACC(:) = P_RR_ACC(:) - ( ZZW2(:)-P_RS_ACC(:) )
!      P_CR_ACC(:) = P_CR_ACC(:) - ( ZZW2(:)-P_RS_ACC(:) ) * PCRT(:)/PRRT(:) ! Lambda_r**3
!      P_RS_ACC(:) = P_RS_ACC(:) - ZZW5(:)
!      P_RG_ACC(:) = ( ZZW2(:)-P_RS_ACC(:) ) + ZZW5(:)
      !
   P_RR_ACC(:) = - ZZW4(:) *  ZZW2(:)
   P_CR_ACC(:) = P_RR_ACC(:) * PCRT(:)/PRRT(:)
   P_RS_ACC(:) = ZZW4(:) *  ZZW1(:) - ZZW5(:)
   P_RG_ACC(:) = ZZW4(:) * ( ZZW2(:) - ZZW1(:) ) + ZZW5(:)
   P_TH_ACC(:) = - P_RR_ACC(:) * (PLSFACT(:)-PLVFACT(:))
   !
END WHERE
!
!
!-------------------------------------------------------------------------------
!
CONTAINS
  FUNCTION GET_XKER_RACCSS(I1,I2) RESULT(RET)
    INTEGER, DIMENSION(:) :: I1
    INTEGER, DIMENSION(:) :: I2
    REAL, DIMENSION(SIZE(I1)) :: RET
    !
    INTEGER I
    !
    DO I=1,SIZE(I1)
       RET(I) = XKER_RACCSS(MAX(MIN(I1(I),SIZE(XKER_RACCSS,1)),1),MAX(MIN(I2(I),SIZE(XKER_RACCSS,2)),1))
    END DO
  END FUNCTION GET_XKER_RACCSS
!
!-------------------------------------------------------------------------------
!
  FUNCTION GET_XKER_RACCS(I1,I2) RESULT(RET)
    INTEGER, DIMENSION(:) :: I1
    INTEGER, DIMENSION(:) :: I2
    REAL, DIMENSION(SIZE(I1)) :: RET
    !
    INTEGER I
    !
    DO I=1,SIZE(I1)
       RET(I) = XKER_RACCS(MAX(MIN(I1(I),SIZE(XKER_RACCS,1)),1),MAX(MIN(I2(I),SIZE(XKER_RACCS,2)),1))
    END DO
  END FUNCTION GET_XKER_RACCS
!
!-------------------------------------------------------------------------------
!
  FUNCTION GET_XKER_SACCRG(I1,I2) RESULT(RET)
    INTEGER, DIMENSION(:) :: I1
    INTEGER, DIMENSION(:) :: I2
    REAL, DIMENSION(SIZE(I1)) :: RET
    !
    INTEGER I
    !
    DO I=1,SIZE(I1)
       RET(I) = XKER_SACCRG(MAX(MIN(I1(I),SIZE(XKER_SACCRG,1)),1),MAX(MIN(I2(I),SIZE(XKER_SACCRG,2)),1))
    END DO
  END FUNCTION GET_XKER_SACCRG
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_RAIN_ACCR_SNOW
