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
INTEGER                       :: J, IB, IE
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
  IB = KB
ELSE
  IB = 1
END IF

IF ( PRESENT(KE) ) THEN
  IE = KE
ELSE
  IE = KT
END IF

!Number of iterations
JITER=2
!
!Computation of PBUF(IB:IE, ICPH2) depending on dummy arguments received
PBUF(IB:IE, ICPH2)=0
IF(PRESENT(PRR)) PBUF(IB:IE, ICPH2)=PBUF(IB:IE, ICPH2) + CST%XCL*PRR(IB:IE)
IF(PRESENT(PRS)) PBUF(IB:IE, ICPH2)=PBUF(IB:IE, ICPH2) + CST%XCI*PRS(IB:IE)
IF(PRESENT(PRG)) PBUF(IB:IE, ICPH2)=PBUF(IB:IE, ICPH2) + CST%XCI*PRG(IB:IE)
IF(PRESENT(PRH)) PBUF(IB:IE, ICPH2)=PBUF(IB:IE, ICPH2) + CST%XCI*PRH(IB:IE)
!
!Computation of an approximate state thanks to PRL and PRI guess
PBUF(IB:IE, IEXN)=(PP(IB:IE)/CST%XP00) ** CST%RDSCPD

DO J=IB,IE
  PBUF(J, I99PP)=0.99*PP(J)
  PRV(J)=PRT(J)-PRL(J)-PRI(J)
  PBUF(J, ICPH)=CST%XCPD+ CST%XCPV * PRV(J)+ CST%XCL * PRL(J) + CST%XCI * PRI(J) + PBUF(J, ICPH2)
  ZVAR2=PBUF(J, ICPH)*PBUF(J, IEXN)
  ZDELT=(PTHL(J)*PBUF(J, IEXN))-CST%XTT
  PBUF(J, ILVOCPEXN) = (CST%XLVTT + (CST%XCPV-CST%XCL) * ZDELT) /ZVAR2
  PBUF(J, ILSOCPEXN) = (CST%XLSTT + (CST%XCPV-CST%XCI) * ZDELT) /ZVAR2 
  PTH(J)=PTHL(J)+PBUF(J, ILVOCPEXN)*PRL(J)+PBUF(J, ILSOCPEXN)*PRI(J)
  PBUF(J, I1PRT)=1+PRT(J)
ENDDO
!
!
!       2 Iteration
!         ---------

DO II=1,JITER
  IF (OOCEAN) THEN
    PBUF(IB:IE, IT)=PTH(IB:IE)
  ELSE
    PBUF(IB:IE, IT)=PTH(IB:IE)*PBUF(IB:IE, IEXN)
  END IF
  !Computation of liquid/ice fractions
  PFRAC_ICE(IB:IE) = 0.
  DO J=IB, IE
    IF(PRL(J)+PRI(J) > 1.E-20) THEN
      PFRAC_ICE(J) = PRI(J) / (PRL(J)+PRI(J))
    ENDIF
  ENDDO
  CALL COMPUTE_FRAC_ICE(HFRAC_ICE,NEBN,PFRAC_ICE(IB:IE),PBUF(IB:IE, IT))

  !Computation of Rvsat and dRsat/dT
  !In this version QSAT, QSATI, DQSAT and DQASATI functions are not used
  !due to performance issue

  ! Log does not vectorize on all compilers:
  PBUF(IB:IE, ILOGT)=LOG(PBUF(IB:IE, IT))

  DO J=IB, IE
    PBUF(J, IFOESW) = MIN(EXP( CST%XALPW - CST%XBETAW/PBUF(J, IT) - CST%XGAMW*PBUF(J, ILOGT)  ), PBUF(J, I99PP))
    PBUF(J, IFOESI) = MIN(EXP( CST%XALPI - CST%XBETAI/PBUF(J, IT) - CST%XGAMI*PBUF(J, ILOGT)  ), PBUF(J, I99PP))
    PRSATW(J) = CST%XRD/CST%XRV*PBUF(J, IFOESW)/PP(J) / (1.+(CST%XRD/CST%XRV-1.)*PBUF(J, IFOESW)/PP(J))
    PRSATI(J) = CST%XRD/CST%XRV*PBUF(J, IFOESI)/PP(J) / (1.+(CST%XRD/CST%XRV-1.)*PBUF(J, IFOESI)/PP(J))
    ZTPOW2=PBUF(J, IT)**2
    PBUF(J, IDRSATODTW) = PRSATW(J) / (1.+(CST%XRD/CST%XRV-1.)*PBUF(J, IFOESW)/PP(J) ) &
                     * (CST%XBETAW/ZTPOW2 - CST%XGAMW/PBUF(J, IT))*PBUF(J, I1PRT)
    PBUF(J, IDRSATODTI) = PRSATI(J) / (1.+(CST%XRD/CST%XRV-1.)*PBUF(J, IFOESI)/PP(J) ) &
                     * (CST%XBETAI/ZTPOW2 - CST%XGAMI/PBUF(J, IT))*PBUF(J, I1PRT)
    !PRSATW(J) =  QSAT(PBUF(J, IT),PP(J)) !qsatw
    !PRSATI(J) = QSATI(PBUF(J, IT),PP(J)) !qsati
    !PBUF(J, IDRSATODTW) =  DQSAT(PBUF(J, IT),PP(J),PRSATW(J))*PBUF(J, I1PRT)
    !PBUF(J, IDRSATODTI) = DQSATI(PBUF(J, IT),PP(J),PRSATI(J))*PBUF(J, I1PRT)
    PRSATW(J) = PRSATW(J)*PBUF(J, I1PRT)
    PRSATI(J) = PRSATI(J)*PBUF(J, I1PRT)
    PBUF(J, IRVSAT) = PRSATW(J)*(1-PFRAC_ICE(J)) + PRSATI(J)*PFRAC_ICE(J)
    PBUF(J, IDRSATODT) = (PBUF(J, IDRSATODTW)*(1-PFRAC_ICE(J))+ &
              & PBUF(J, IDRSATODTI)*PFRAC_ICE(J))

    !Computation of new PRL, PRI and PRV
    !Correction term applied to (PRV(J)-PBUF(J, IRVSAT)) is computed assuming that
    !PBUF(J, ILVOCPEXN), PBUF(J, ILSOCPEXN) and PBUF(J, ICPH) don't vary too much with T. It takes into account
    !the variation (estimated linear) of Qsat with T
    PBUF(J, IRLTEMP)=(PRV(J)-PBUF(J, IRVSAT))/ &
                  &(1 + PBUF(J, IDRSATODT)*PBUF(J, IEXN)* &
                  &     (PBUF(J, ILVOCPEXN)*(1-PFRAC_ICE(J))+PBUF(J, ILSOCPEXN)*PFRAC_ICE(J)))
    PBUF(J, IRLTEMP)=MIN(MAX(-PRL(J)-PRI(J), PBUF(J, IRLTEMP)),PRV(J))
    PRV(J)=PRV(J)-PBUF(J, IRLTEMP)
    PRL(J)=PRL(J)+PRI(J)+PBUF(J, IRLTEMP)
    PRI(J)=PFRAC_ICE(J)     * (PRL(J))
    PRL(J)=(1-PFRAC_ICE(J)) * (PRT(J) - PRV(J))

    !Computation of Cph (as defined in Meso-NH doc, equation 2.2, to be used with mixing ratios)
    PBUF(J, ICPH)=CST%XCPD+ CST%XCPV * PRV(J)+ CST%XCL * PRL(J) + CST%XCI * PRI(J) + PBUF(J, ICPH2)

    !Computation of L/Cph/EXN, then new PTH
    ZVAR2=PBUF(J, ICPH)*PBUF(J, IEXN)
    PBUF(J, ILVOCPEXN) = (CST%XLVTT + (CST%XCPV-CST%XCL) * (PBUF(J, IT)-CST%XTT)) /ZVAR2
    PBUF(J, ILSOCPEXN) = (CST%XLSTT + (CST%XCPV-CST%XCI) * (PBUF(J, IT)-CST%XTT)) /ZVAR2
    PTH(J)=PTHL(J)+PBUF(J, ILVOCPEXN)*PRL(J)+PBUF(J, ILSOCPEXN)*PRI(J)

    !Computation of estimated mixing ration at saturation
    !To compute the adjustement a first order development was used
    ZVAR1=PTH(J)*PBUF(J, IEXN)-PBUF(J, IT)
    PRSATW(J)=PRSATW(J) + PBUF(J, IDRSATODTW)*ZVAR1
    PRSATI(J)=PRSATI(J) + PBUF(J, IDRSATODTI)*ZVAR1
  ENDDO
ENDDO

END SUBROUTINE TH_R_FROM_THL_RT
