!MNH_LIC Copyright 2006-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ######spl
      MODULE MODI_TH_R_FROM_THL_RT_1D
!     ###############################
!
      INTERFACE
      SUBROUTINE TH_R_FROM_THL_RT_1D(HFRAC_ICE,PFRAC_ICE,PP,             &
                                  PTHL, PRT, PTH, PRV, PRL, PRI,         &
                                  PRSATW, PRSATI, PRR, PRS, PRG, PRH     )
CHARACTER(len=1),   INTENT(IN)    :: HFRAC_ICE
REAL, DIMENSION(:), INTENT(INOUT) :: PFRAC_ICE
REAL, DIMENSION(:), INTENT(IN) :: PP     ! Pressure
REAL, DIMENSION(:), INTENT(IN) :: PTHL   ! Liquid pot. temp.
REAL, DIMENSION(:), INTENT(IN) :: PRT    ! Total mixing ratios 
REAL, DIMENSION(:),OPTIONAL,INTENT(IN) :: PRR, PRS, PRG, PRH
REAL, DIMENSION(:), INTENT(OUT):: PTH    ! Potential temp.
REAL, DIMENSION(:), INTENT(OUT):: PRV    ! vapor mixing ratio
REAL, DIMENSION(:), INTENT(INOUT):: PRL  ! cloud mixing ratio
REAL, DIMENSION(:), INTENT(INOUT):: PRI  ! ice   mixing ratio
REAL, DIMENSION(:), INTENT(OUT)  :: PRSATW ! estimated mixing ration at saturation over water
REAL, DIMENSION(:), INTENT(OUT)  :: PRSATI ! estimated mixing ration at saturation over ice

      END SUBROUTINE TH_R_FROM_THL_RT_1D
      END INTERFACE
      END MODULE MODI_TH_R_FROM_THL_RT_1D
!     ######spl
      SUBROUTINE TH_R_FROM_THL_RT_1D(HFRAC_ICE,PFRAC_ICE,PP,             &
                                  PTHL, PRT, PTH, PRV, PRL, PRI,         &
                                  PRSATW, PRSATI, PRR, PRS, PRG, PRH     )
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
USE MODI_COMPUTE_FRAC_ICE
USE MODD_CST
USE MODD_DYN_n, ONLY : LOCEAN
USE MODE_THERMO
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
CHARACTER(len=1),   INTENT(IN) :: HFRAC_ICE
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
REAL, DIMENSION(SIZE(PTHL,1))   :: ZEXN
REAL, DIMENSION(SIZE(PTHL,1)) :: ZRVSAT,ZCPH,ZRLTEMP,ZCPH2
REAL, DIMENSION(SIZE(PTHL,1)) :: ZT,ZLVOCPEXN,ZLSOCPEXN
REAL, DIMENSION(SIZE(PTHL,1)) :: ZDRSATODT,ZDRSATODTW,ZDRSATODTI
REAL, DIMENSION(SIZE(PTHL,1)) :: ZFOESW, ZFOESI
!----------------------------------------------------------------------------
!
!*      1 Initialisation
!         --------------
!
!
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
ZEXN(:)=(PP(:)/XP00) ** (XRD/XCPD)
ZT(:)=PTHL(:)*ZEXN(:)
PRV(:)=PRT(:)-PRL(:)-PRI(:)
ZCPH(:)=XCPD+ XCPV * PRV(:)+ XCL * PRL(:) + XCI * PRI(:) + ZCPH2(:)
ZLVOCPEXN(:) = (XLVTT + (XCPV-XCL) * (ZT(:)-XTT)) &
                        /(ZCPH(:)*ZEXN(:))
ZLSOCPEXN(:) = (XLSTT + (XCPV-XCI) * (ZT(:)-XTT)) &
                        /(ZCPH(:)*ZEXN(:))
PTH(:)=PTHL(:)+ZLVOCPEXN(:)*PRL(:)+ZLSOCPEXN(:)*PRI(:)
!
!
!       2 Iteration
!         ---------

DO II=1,JITER
  IF (LOCEAN) THEN
    ZT=PTH                  
  ELSE
    ZT(:)=PTH(:)*ZEXN(:)
  END IF
  !Computation of liquid/ice fractions
  PFRAC_ICE(:) = 0.
  WHERE(PRL(:)+PRI(:) > 1.E-20)
    PFRAC_ICE(:) = PRI(:) / (PRL(:)+PRI(:))
  ENDWHERE
  CALL COMPUTE_FRAC_ICE(HFRAC_ICE,PFRAC_ICE(:),ZT(:))

  !Computation of Rvsat and dRsat/dT
  !In this version QSAT, QSATI, DQSAT and DQASATI functions are not used
  !due to performance issue
  ZFOESW(:) = MIN(EXP( XALPW - XBETAW/ZT(:) - XGAMW*LOG(ZT(:))  ), 0.99*PP(:))
  ZFOESI(:) = MIN(EXP( XALPI - XBETAI/ZT(:) - XGAMI*LOG(ZT(:))  ), 0.99*PP(:))
  PRSATW(:) = XRD/XRV*ZFOESW(:)/PP(:) / (1.+(XRD/XRV-1.)*ZFOESW(:)/PP(:))
  PRSATI(:) = XRD/XRV*ZFOESI(:)/PP(:) / (1.+(XRD/XRV-1.)*ZFOESI(:)/PP(:))
  ZDRSATODTW(:) = PRSATW(:) / (1.+(XRD/XRV-1.)*ZFOESW(:)/PP(:) ) &
                   * (XBETAW/ZT(:)**2 - XGAMW/ZT(:))*(1+PRT(:))
  ZDRSATODTI(:) = PRSATI(:) / (1.+(XRD/XRV-1.)*ZFOESI(:)/PP(:) ) &
                   * (XBETAI/ZT(:)**2 - XGAMI/ZT(:))*(1+PRT(:))
  !PRSATW(:) =  QSAT(ZT(:),PP(:)) !qsatw
  !PRSATI(:) = QSATI(ZT(:),PP(:)) !qsati
  !ZDRSATODTW(:) =  DQSAT(ZT(:),PP(:),PRSATW(:))*(1+PRT(:))
  !ZDRSATODTI(:) = DQSATI(ZT(:),PP(:),PRSATI(:))*(1+PRT(:))
  PRSATW(:) = PRSATW(:)*(1+PRT(:))
  PRSATI(:) = PRSATI(:)*(1+PRT(:))
  ZRVSAT(:) = PRSATW(:)*(1-PFRAC_ICE(:)) + PRSATI(:)*PFRAC_ICE(:)
  ZDRSATODT(:) = (ZDRSATODTW(:)*(1-PFRAC_ICE(:))+ &
            & ZDRSATODTI(:)*PFRAC_ICE(:))

  !Computation of new PRL, PRI and PRV
  !Correction term applied to (PRV(:)-ZRVSAT(:)) is computed assuming that
  !ZLVOCPEXN, ZLSOCPEXN and ZCPH don't vary to much with T. It takes into account
  !the variation (estimated linear) of Qsat with T
  ZRLTEMP(:)=(PRV(:)-ZRVSAT(:))/ &
                &(1 + ZDRSATODT(:)*ZEXN(:)* &
                &     (ZLVOCPEXN(:)*(1-PFRAC_ICE(:))+ZLSOCPEXN(:)*PFRAC_ICE(:)))
  ZRLTEMP(:)=MIN(MAX(-PRL(:)-PRI(:), ZRLTEMP(:)),PRV(:))
  PRV(:)=PRV(:)-ZRLTEMP(:)
  PRL(:)=PRL(:)+PRI(:)+ZRLTEMP(:)
  PRI(:)=PFRAC_ICE(:)     * (PRL(:))
  PRL(:)=(1-PFRAC_ICE(:)) * (PRT(:) - PRV(:))

  !Computation of Cph (as defined in Meso-NH doc, equation 2.2, to be used with mixing ratios)
  ZCPH(:)=XCPD+ XCPV * PRV(:)+ XCL * PRL(:) + XCI * PRI(:) + ZCPH2(:)

  !Computation of L/Cph/EXN, then new PTH
  ZLVOCPEXN(:) = (XLVTT + (XCPV-XCL) * (ZT(:)-XTT)) &
                    /(ZCPH(:)*ZEXN(:))
  ZLSOCPEXN(:) = (XLSTT + (XCPV-XCI) * (ZT(:)-XTT)) &
                    /(ZCPH(:)*ZEXN(:))
  PTH(:)=PTHL(:)+ZLVOCPEXN(:)*PRL(:)+ZLSOCPEXN(:)*PRI(:)

  !Computation of estimated mixing ration at saturation
  !To compute the adjustement a first order development was used
  PRSATW(:)=PRSATW(:) + ZDRSATODTW(:)*(PTH(:)*ZEXN(:)-ZT(:))
  PRSATI(:)=PRSATI(:) + ZDRSATODTI(:)*(PTH(:)*ZEXN(:)-ZT(:))
ENDDO


END SUBROUTINE TH_R_FROM_THL_RT_1D
