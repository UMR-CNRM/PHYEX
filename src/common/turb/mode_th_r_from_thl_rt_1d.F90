!MNH_LIC Copyright 2006-2022 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
MODULE MODE_TH_R_FROM_THL_RT_1D
IMPLICIT NONE
CONTAINS
      SUBROUTINE TH_R_FROM_THL_RT_1D(HFRAC_ICE,PFRAC_ICE,PP,             &
                                  PTHL, PRT, PTH, PRV, PRL, PRI,         &
                                  PRSATW, PRSATI, PRR, PRS, PRG, PRH,OOCEAN)
!     #################################################################
!
!
!!****  *TH_R_FROM_THL_RT_1D* - computes the non-conservative variables
!!                          from conservative variables
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!    
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      Julien PERGAUD      * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original         13/03/06
!!      S. Riette April 2011 : ice added, allow ZRLTEMP to be negative
!!                             we use dQsat/dT to help convergence
!!                             use of optional PRR, PRS, PRG, PRH
!!      S. Riette Nov 2016: support for HFRAC_ICE='S'
!!
!! --------------------------------------------------------------------------
!
!*      0. DECLARATIONS
!          ------------
!
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
USE MODD_CST !, ONLY: XP00, XRD, XCPD, XCPV, XCL, XCI, XLVTT, XTT, XLSTT
USE MODD_NEB, ONLY: NEB
USE MODE_THERMO
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
CHARACTER(LEN=1),   INTENT(IN) :: HFRAC_ICE
LOGICAL,            INTENT(IN)   ::  OOCEAN       ! switch for Ocean model version
REAL, DIMENSION(:), INTENT(INOUT) :: PFRAC_ICE
REAL, DIMENSION(:), INTENT(IN) :: PP          ! Pressure
REAL, DIMENSION(:), INTENT(IN) :: PTHL    ! thetal to transform into th
REAL, DIMENSION(:),INTENT(IN)  :: PRT    ! Total mixing ratios to transform into rv,rc and ri
REAL, DIMENSION(:),OPTIONAL,INTENT(IN) :: PRR, PRS, PRG, PRH
REAL, DIMENSION(:), INTENT(OUT):: PTH    ! th
REAL, DIMENSION(:), INTENT(OUT):: PRV    ! vapor mixing ratio
REAL, DIMENSION(:), INTENT(INOUT):: PRL    ! vapor mixing ratio
REAL, DIMENSION(:), INTENT(INOUT):: PRI    ! vapor mixing ratio
REAL, DIMENSION(:), INTENT(OUT)  :: PRSATW ! estimated mixing ration at saturation over water
REAL, DIMENSION(:), INTENT(OUT)  :: PRSATI ! estimated mixing ration at saturation over ice
!
!-------------------------------------------------------------------------------
!
!       0.2  declaration of local variables
INTEGER                       :: II ! Loop control
INTEGER                       :: JITER ! number of iterations
INTEGER                       :: J
REAL, DIMENSION(SIZE(PTHL,1))   :: ZEXN
REAL, DIMENSION(SIZE(PTHL,1)) :: ZRVSAT,ZCPH,ZRLTEMP,ZCPH2
REAL, DIMENSION(SIZE(PTHL,1)) :: ZT,ZLVOCPEXN,ZLSOCPEXN
REAL, DIMENSION(SIZE(PTHL,1)) :: ZDRSATODT,ZDRSATODTW,ZDRSATODTI
REAL, DIMENSION(SIZE(PTHL,1)) :: ZFOESW, ZFOESI
REAL, DIMENSION(SIZE(PTHL,1)) :: ZLOGT, Z99PP, Z1PRT
REAL(KIND=JPRB) :: ZVAR1, ZVAR2, ZTPOW2, ZDELT
INTEGER, DIMENSION(SIZE(PTHL,1)) :: IERR

REAL(KIND=JPRB) :: ZHOOK_HANDLE
!----------------------------------------------------------------------------
!
!*      1 Initialisation
!         --------------
!
!
IF (LHOOK) CALL DR_HOOK('TH_R_FROM_THL_RT_1D',0,ZHOOK_HANDLE)
!
!Number of iterations
JITER=2
!
!Computation of ZCPH2 depending on dummy arguments received
ZCPH2(:)=0
IF(PRESENT(PRR)) ZCPH2(:)=ZCPH2(:) + XCL*PRR(:)
IF(PRESENT(PRS)) ZCPH2(:)=ZCPH2(:) + XCI*PRS(:)
IF(PRESENT(PRG)) ZCPH2(:)=ZCPH2(:) + XCI*PRG(:)
IF(PRESENT(PRH)) ZCPH2(:)=ZCPH2(:) + XCI*PRH(:)
!
!Computation of an approximate state thanks to PRL and PRI guess
ZEXN(:)=(PP(:)/XP00) ** RDSCPD

DO J=1,SIZE(PTHL,1)
Z99PP(J)=0.99*PP(J)
PRV(J)=PRT(J)-PRL(J)-PRI(J)
ZCPH(J)=XCPD+ XCPV * PRV(J)+ XCL * PRL(J) + XCI * PRI(J) + ZCPH2(J)
ZVAR2=ZCPH(J)*ZEXN(J)
ZDELT=(PTHL(J)*ZEXN(J))-XTT
ZLVOCPEXN(J) = (XLVTT + (XCPV-XCL) * ZDELT) /ZVAR2
ZLSOCPEXN(J) = (XLSTT + (XCPV-XCI) * ZDELT) /ZVAR2 
PTH(J)=PTHL(J)+ZLVOCPEXN(J)*PRL(J)+ZLSOCPEXN(J)*PRI(J)
Z1PRT(J)=1+PRT(J)
ENDDO
!
!
!       2 Iteration
!         ---------

DO II=1,JITER
  IF (OOCEAN) THEN
    ZT=PTH                  
  ELSE
    ZT(:)=PTH(:)*ZEXN(:)
  END IF
  !Computation of liquid/ice fractions
  PFRAC_ICE(:) = 0.
  DO J=1, SIZE(PFRAC_ICE, 1)
    IF(PRL(J)+PRI(J) > 1.E-20) THEN
      PFRAC_ICE(J) = PRI(J) / (PRL(J)+PRI(J))
    ENDIF
  ENDDO
  CALL COMPUTE_FRAC_ICE(HFRAC_ICE,NEB,PFRAC_ICE(:),ZT(:), IERR(:))

  !Computation of Rvsat and dRsat/dT
  !In this version QSAT, QSATI, DQSAT and DQASATI functions are not used
  !due to performance issue

  ! Log does not vectorize on all compilers:
  ZLOGT(:)=LOG(ZT(:))

  DO J=1,SIZE(PTHL,1)

  ZFOESW(J) = MIN(EXP( XALPW - XBETAW/ZT(J) - XGAMW*ZLOGT(J)  ), Z99PP(J))
  ZFOESI(J) = MIN(EXP( XALPI - XBETAI/ZT(J) - XGAMI*ZLOGT(J)  ), Z99PP(J))
  PRSATW(J) = XRD/XRV*ZFOESW(J)/PP(J) / (1.+(XRD/XRV-1.)*ZFOESW(J)/PP(J))
  PRSATI(J) = XRD/XRV*ZFOESI(J)/PP(J) / (1.+(XRD/XRV-1.)*ZFOESI(J)/PP(J))
  ZTPOW2=ZT(J)**2
  ZDRSATODTW(J) = PRSATW(J) / (1.+(XRD/XRV-1.)*ZFOESW(J)/PP(J) ) &
                   * (XBETAW/ZTPOW2 - XGAMW/ZT(J))*Z1PRT(J)
  ZDRSATODTI(J) = PRSATI(J) / (1.+(XRD/XRV-1.)*ZFOESI(J)/PP(J) ) &
                   * (XBETAI/ZTPOW2 - XGAMI/ZT(J))*Z1PRT(J)
  !PRSATW(J) =  QSAT(ZT(J),PP(J)) !qsatw
  !PRSATI(J) = QSATI(ZT(J),PP(J)) !qsati
  !ZDRSATODTW(J) =  DQSAT(ZT(J),PP(J),PRSATW(J))*Z1PRT(J)
  !ZDRSATODTI(J) = DQSATI(ZT(J),PP(J),PRSATI(J))*Z1PRT(J)
  PRSATW(J) = PRSATW(J)*Z1PRT(J)
  PRSATI(J) = PRSATI(J)*Z1PRT(J)
  ZRVSAT(J) = PRSATW(J)*(1-PFRAC_ICE(J)) + PRSATI(J)*PFRAC_ICE(J)
  ZDRSATODT(J) = (ZDRSATODTW(J)*(1-PFRAC_ICE(J))+ &
            & ZDRSATODTI(J)*PFRAC_ICE(J))

  !Computation of new PRL, PRI and PRV
  !Correction term applied to (PRV(J)-ZRVSAT(J)) is computed assuming that
  !ZLVOCPEXN, ZLSOCPEXN and ZCPH don't vary to much with T. It takes into account
  !the variation (estimated linear) of Qsat with T
  ZRLTEMP(J)=(PRV(J)-ZRVSAT(J))/ &
                &(1 + ZDRSATODT(J)*ZEXN(J)* &
                &     (ZLVOCPEXN(J)*(1-PFRAC_ICE(J))+ZLSOCPEXN(J)*PFRAC_ICE(J)))
  ZRLTEMP(J)=MIN(MAX(-PRL(J)-PRI(J), ZRLTEMP(J)),PRV(J))
  PRV(J)=PRV(J)-ZRLTEMP(J)
  PRL(J)=PRL(J)+PRI(J)+ZRLTEMP(J)
  PRI(J)=PFRAC_ICE(J)     * (PRL(J))
  PRL(J)=(1-PFRAC_ICE(J)) * (PRT(J) - PRV(J))

  !Computation of Cph (as defined in Meso-NH doc, equation 2.2, to be used with mixing ratios)
  ZCPH(J)=XCPD+ XCPV * PRV(J)+ XCL * PRL(J) + XCI * PRI(J) + ZCPH2(J)

  !Computation of L/Cph/EXN, then new PTH
  ZVAR2=ZCPH(J)*ZEXN(J)
  ZLVOCPEXN(J) = (XLVTT + (XCPV-XCL) * (ZT(J)-XTT)) /ZVAR2
  ZLSOCPEXN(J) = (XLSTT + (XCPV-XCI) * (ZT(J)-XTT)) /ZVAR2
  PTH(J)=PTHL(J)+ZLVOCPEXN(J)*PRL(J)+ZLSOCPEXN(J)*PRI(J)

  !Computation of estimated mixing ration at saturation
  !To compute the adjustement a first order development was used
  ZVAR1=PTH(J)*ZEXN(J)-ZT(J)
  PRSATW(J)=PRSATW(J) + ZDRSATODTW(J)*ZVAR1
  PRSATI(J)=PRSATI(J) + ZDRSATODTI(J)*ZVAR1

  ENDDO
ENDDO

IF (LHOOK) CALL DR_HOOK('TH_R_FROM_THL_RT_1D',1,ZHOOK_HANDLE)

!
CONTAINS
INCLUDE "compute_frac_ice.func.h"
!
END SUBROUTINE TH_R_FROM_THL_RT_1D
END MODULE MODE_TH_R_FROM_THL_RT_1D
