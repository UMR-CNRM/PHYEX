!MNH_LIC Copyright 2006-2023 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
      SUBROUTINE TH_R_FROM_THL_RT(CST, NEBN, KT, HFRAC_ICE,PFRAC_ICE,PP, &
                                  PTHL, PRT, PTH, PRV, PRL, PRI,           &
                                  PRSATW, PRSATI, PRR, PRS, PRG, PRH, OOCEAN,&
                                  PBUF, KB, KE)
! ******* TO BE INCLUDED IN THE *CONTAINS* OF A SUBROUTINE, IN ORDER TO EASE AUTOMATIC INLINING ******
! => Don't use drHook !!!
! "compute_frac_ice.func.h" must be included at the same time
!     #################################################################
!
!
!!****  *TH_R_FROM_THL_RT* - computes the non-conservative variables
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
!!      P. Wautelet Mar 2023: bugfix: protect operations on optional dummy arguments
!!
!! --------------------------------------------------------------------------
!
!*      0. DECLARATIONS
!          ------------
!
USE MODD_CST, ONLY : CST_t
USE MODD_NEB_n, ONLY : NEB_t
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
TYPE(CST_t),        INTENT(IN) :: CST
TYPE(NEB_t),        INTENT(IN) :: NEBN
INTEGER,            INTENT(IN) :: KT
CHARACTER(LEN=1),   INTENT(IN) :: HFRAC_ICE
LOGICAL,            INTENT(IN)   ::  OOCEAN       ! switch for Ocean model version
REAL, DIMENSION(KT), INTENT(INOUT) :: PFRAC_ICE
REAL, DIMENSION(KT), INTENT(IN) :: PP          ! Pressure
REAL, DIMENSION(KT), INTENT(IN) :: PTHL    ! thetal to transform into th
REAL, DIMENSION(KT),INTENT(IN)  :: PRT    ! Total mixing ratios to transform into rv,rc and ri
REAL, DIMENSION(KT),OPTIONAL,INTENT(IN) :: PRR, PRS, PRG, PRH
REAL, DIMENSION(KT), INTENT(OUT):: PTH    ! th
REAL, DIMENSION(KT), INTENT(OUT):: PRV    ! vapor mixing ratio
REAL, DIMENSION(KT), INTENT(INOUT):: PRL    ! vapor mixing ratio
REAL, DIMENSION(KT), INTENT(INOUT):: PRI    ! vapor mixing ratio
REAL, DIMENSION(KT), INTENT(OUT)  :: PRSATW ! estimated mixing ration at saturation over water
REAL, DIMENSION(KT), INTENT(OUT)  :: PRSATI ! estimated mixing ration at saturation over ice
REAL, DIMENSION(KT, 16), INTENT(OUT) :: PBUF ! buffer to replace automatic arrays
INTEGER, OPTIONAL, INTENT(IN) :: KB !first index to deal with (default is 1)
INTEGER, OPTIONAL, INTENT(IN) :: KE !last index to deal with (default if KT)
!
!-------------------------------------------------------------------------------
!
!       0.2  declaration of local variables
INTEGER                       :: II ! Loop control
INTEGER                       :: JITER ! number of iterations
INTEGER                       :: JIJ, IIJB, IIJE
INTEGER, PARAMETER :: IEXN=1, IRVSAT=2, ICPH=3, IRLTEMP=4, ICPH2=5, IT=6, ILVOCPEXN=7, ILSOCPEXN=8, &
                    & IDRSATODT=9, IDRSATODTW=10, IDRSATODTI=11, IFOESW=12, IFOESI=13, &
                    & ILOGT=14, I99PP=15, I1PRT=16
REAL :: ZVAR1, ZVAR2, ZTPOW2, ZDELT

!----------------------------------------------------------------------------
!
!*      1 Initialisation
!         --------------
!
!
!
IF ( PRESENT(KB) ) THEN
  IIJB = KB
ELSE
  IIJB = 1
END IF

IF ( PRESENT(KE) ) THEN
  IIJE = KE
ELSE
  IIJE = KT
END IF

!Number of iterations
JITER=2
!
!Computation of PBUF(:, ICPH2) depending on dummy arguments received
PBUF(:, ICPH2)=0
IF(PRESENT(PRR)) PBUF(:, ICPH2)=PBUF(:, ICPH2) + CST%XCL*PRR(:)
IF(PRESENT(PRS)) PBUF(:, ICPH2)=PBUF(:, ICPH2) + CST%XCI*PRS(:)
IF(PRESENT(PRG)) PBUF(:, ICPH2)=PBUF(:, ICPH2) + CST%XCI*PRG(:)
IF(PRESENT(PRH)) PBUF(:, ICPH2)=PBUF(:, ICPH2) + CST%XCI*PRH(:)
!
!Computation of an approximate state thanks to PRL and PRI guess
PBUF(:, IEXN)=(PP(:)/CST%XP00) ** CST%RDSCPD

DO JIJ=IIJB,IIJE
  PBUF(JIJ, I99PP)=0.99*PP(JIJ)
  PRV(JIJ)=PRT(JIJ)-PRL(JIJ)-PRI(JIJ)
  PBUF(JIJ, ICPH)=CST%XCPD+ CST%XCPV * PRV(JIJ)+ CST%XCL * PRL(JIJ) + CST%XCI * PRI(JIJ) + PBUF(JIJ, ICPH2)
  ZVAR2=PBUF(JIJ, ICPH)*PBUF(JIJ, IEXN)
  ZDELT=(PTHL(JIJ)*PBUF(JIJ, IEXN))-CST%XTT
  PBUF(JIJ, ILVOCPEXN) = (CST%XLVTT + (CST%XCPV-CST%XCL) * ZDELT) /ZVAR2
  PBUF(JIJ, ILSOCPEXN) = (CST%XLSTT + (CST%XCPV-CST%XCI) * ZDELT) /ZVAR2 
  PTH(JIJ)=PTHL(JIJ)+PBUF(JIJ, ILVOCPEXN)*PRL(JIJ)+PBUF(JIJ, ILSOCPEXN)*PRI(JIJ)
  PBUF(JIJ, I1PRT)=1+PRT(JIJ)
ENDDO
!
!
!       2 Iteration
!         ---------

DO II=1,JITER
  IF (OOCEAN) THEN
    PBUF(:, IT)=PTH(:)
  ELSE
    PBUF(:, IT)=PTH(:)*PBUF(:, IEXN)
  END IF
  !Computation of liquid/ice fractions
  PFRAC_ICE(:) = 0.
  DO JIJ=IIJB, IIJE
    IF(PRL(JIJ)+PRI(JIJ) > 1.E-20) THEN
      PFRAC_ICE(JIJ) = PRI(JIJ) / (PRL(JIJ)+PRI(JIJ))
    ENDIF
  ENDDO
  CALL COMPUTE_FRAC_ICE(CST, HFRAC_ICE,NEBN,PFRAC_ICE(:),PBUF(:, IT))

  !Computation of Rvsat and dRsat/dT
  !In this version QSAT, QSATI, DQSAT and DQASATI functions are not used
  !due to performance issue

  ! Log does not vectorize on all compilers:
  PBUF(:, ILOGT)=LOG(PBUF(:, IT))

  DO JIJ=IIJB, IIJE
    PBUF(JIJ, IFOESW) = MIN(EXP( CST%XALPW - CST%XBETAW/PBUF(JIJ, IT) - CST%XGAMW*PBUF(JIJ, ILOGT)  ), PBUF(JIJ, I99PP))
    PBUF(JIJ, IFOESI) = MIN(EXP( CST%XALPI - CST%XBETAI/PBUF(JIJ, IT) - CST%XGAMI*PBUF(JIJ, ILOGT)  ), PBUF(JIJ, I99PP))
    PRSATW(JIJ) = CST%XRD/CST%XRV*PBUF(JIJ, IFOESW)/PP(JIJ) / (1.+(CST%XRD/CST%XRV-1.)*PBUF(JIJ, IFOESW)/PP(JIJ))
    PRSATI(JIJ) = CST%XRD/CST%XRV*PBUF(JIJ, IFOESI)/PP(JIJ) / (1.+(CST%XRD/CST%XRV-1.)*PBUF(JIJ, IFOESI)/PP(JIJ))
    ZTPOW2=PBUF(JIJ, IT)**2
    PBUF(JIJ, IDRSATODTW) = PRSATW(JIJ) / (1.+(CST%XRD/CST%XRV-1.)*PBUF(JIJ, IFOESW)/PP(JIJ) ) &
                     * (CST%XBETAW/ZTPOW2 - CST%XGAMW/PBUF(JIJ, IT))*PBUF(JIJ, I1PRT)
    PBUF(JIJ, IDRSATODTI) = PRSATI(JIJ) / (1.+(CST%XRD/CST%XRV-1.)*PBUF(JIJ, IFOESI)/PP(JIJ) ) &
                     * (CST%XBETAI/ZTPOW2 - CST%XGAMI/PBUF(JIJ, IT))*PBUF(JIJ, I1PRT)
    !PRSATW(JIJ) =  QSAT(PBUF(JIJ, IT),PP(JIJ)) !qsatw
    !PRSATI(JIJ) = QSATI(PBUF(JIJ, IT),PP(JIJ)) !qsati
    !PBUF(JIJ, IDRSATODTW) =  DQSAT(PBUF(JIJ, IT),PP(JIJ),PRSATW(JIJ))*PBUF(JIJ, I1PRT)
    !PBUF(JIJ, IDRSATODTI) = DQSATI(PBUF(JIJ, IT),PP(JIJ),PRSATI(JIJ))*PBUF(JIJ, I1PRT)
    PRSATW(JIJ) = PRSATW(JIJ)*PBUF(JIJ, I1PRT)
    PRSATI(JIJ) = PRSATI(JIJ)*PBUF(JIJ, I1PRT)
    PBUF(JIJ, IRVSAT) = PRSATW(JIJ)*(1-PFRAC_ICE(JIJ)) + PRSATI(JIJ)*PFRAC_ICE(JIJ)
    PBUF(JIJ, IDRSATODT) = (PBUF(JIJ, IDRSATODTW)*(1-PFRAC_ICE(JIJ))+ &
              & PBUF(JIJ, IDRSATODTI)*PFRAC_ICE(JIJ))

    !Computation of new PRL, PRI and PRV
    !Correction term applied to (PRV(JIJ)-PBUF(JIJ, IRVSAT)) is computed assuming that
    !PBUF(JIJ, ILVOCPEXN), PBUF(JIJ, ILSOCPEXN) and PBUF(JIJ, ICPH) don't vary too much with T. It takes into account
    !the variation (estimated linear) of Qsat with T
    PBUF(JIJ, IRLTEMP)=(PRV(JIJ)-PBUF(JIJ, IRVSAT))/ &
                  &(1 + PBUF(JIJ, IDRSATODT)*PBUF(JIJ, IEXN)* &
                  &     (PBUF(JIJ, ILVOCPEXN)*(1-PFRAC_ICE(JIJ))+PBUF(JIJ, ILSOCPEXN)*PFRAC_ICE(JIJ)))
    PBUF(JIJ, IRLTEMP)=MIN(MAX(-PRL(JIJ)-PRI(JIJ), PBUF(JIJ, IRLTEMP)),PRV(JIJ))
    PRV(JIJ)=PRV(JIJ)-PBUF(JIJ, IRLTEMP)
    PRL(JIJ)=PRL(JIJ)+PRI(JIJ)+PBUF(JIJ, IRLTEMP)
    PRI(JIJ)=PFRAC_ICE(JIJ)     * (PRL(JIJ))
    PRL(JIJ)=(1-PFRAC_ICE(JIJ)) * (PRT(JIJ) - PRV(JIJ))

    !Computation of Cph (as defined in Meso-NH doc, equation 2.2, to be used with mixing ratios)
    PBUF(JIJ, ICPH)=CST%XCPD+ CST%XCPV * PRV(JIJ)+ CST%XCL * PRL(JIJ) + CST%XCI * PRI(JIJ) + PBUF(JIJ, ICPH2)

    !Computation of L/Cph/EXN, then new PTH
    ZVAR2=PBUF(JIJ, ICPH)*PBUF(JIJ, IEXN)
    PBUF(JIJ, ILVOCPEXN) = (CST%XLVTT + (CST%XCPV-CST%XCL) * (PBUF(JIJ, IT)-CST%XTT)) /ZVAR2
    PBUF(JIJ, ILSOCPEXN) = (CST%XLSTT + (CST%XCPV-CST%XCI) * (PBUF(JIJ, IT)-CST%XTT)) /ZVAR2
    PTH(JIJ)=PTHL(JIJ)+PBUF(JIJ, ILVOCPEXN)*PRL(JIJ)+PBUF(JIJ, ILSOCPEXN)*PRI(JIJ)

    !Computation of estimated mixing ration at saturation
    !To compute the adjustement a first order development was used
    ZVAR1=PTH(JIJ)*PBUF(JIJ, IEXN)-PBUF(JIJ, IT)
    PRSATW(JIJ)=PRSATW(JIJ) + PBUF(JIJ, IDRSATODTW)*ZVAR1
    PRSATI(JIJ)=PRSATI(JIJ) + PBUF(JIJ, IDRSATODTI)*ZVAR1
  ENDDO
ENDDO

END SUBROUTINE TH_R_FROM_THL_RT
