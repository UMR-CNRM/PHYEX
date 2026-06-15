!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_DEAR

!$ACDC singlecolumn --dummy

IMPLICIT NONE
CONTAINS
  SUBROUTINE DEAR (D, CST, TURBN, KRR, KRRI, O2D, OCOMPUTE_SRC, OOCEAN, &
                & PLM,PRT, PDZZ, PZZ, PTKET, PTHVREF, &
                & PTHLT, PLOCPEXNM, PSRCT, PAMOIST, PDIRCOSZW, &
                & PDXX, PDYY, PATHETA)
    !     ####################
    !!
    !!****  *DEAR* routine to compute mixing length for DEARdorff case
    !
    !!    AUTHOR
    !!    ------
    !!
    !!     M Tomasini      *Meteo-France
    !!
    !!    MODIFICATIONS
    !!    -------------
    !!      Original   01/05
    !!      I.Sandu (Sept.2006) : Modification of the stability criterion
    !!                            (theta_v -> theta_l)
    !!
    !-------------------------------------------------------------------------------
    !
    !*       0.    DECLARATIONS
    !              ------------
    USE YOMHOOK,         ONLY: LHOOK, DR_HOOK, JPHOOK
    USE MODD_CST,        ONLY: CST_t
    USE MODD_DIMPHYEX,   ONLY: DIMPHYEX_t
    USE MODD_TURB_n,     ONLY: TURB_t
    USE MODE_SHUMAN_PHY, ONLY: MXF_PHY,MYF_PHY
    USE MODE_ETHETA,     ONLY: ETHETA
    USE MODE_EMOIST,     ONLY: EMOIST

    IMPLICIT NONE
    !
    !*       0.1   Declarations of dummy arguments
    !
    TYPE(DIMPHYEX_t), INTENT(IN) :: D
    TYPE(CST_t),      INTENT(IN) :: CST
    TYPE(TURB_t),     INTENT(IN) :: TURBN
    LOGICAL, INTENT(IN) :: OCOMPUTE_SRC
    LOGICAL, INTENT(IN) :: O2D
    LOGICAL, INTENT(IN) :: OOCEAN
    INTEGER, INTENT(IN) :: KRRI, KRR
    REAL, INTENT(IN) :: PDZZ(D%NIJT, D%NKT)
    REAL, INTENT(IN) :: PZZ(D%NIJT, D%NKT)
    REAL, INTENT(IN) :: PTKET(D%NIJT, D%NKT)
    REAL, INTENT(IN) :: PDIRCOSZW(D%NIJT)
    REAL, INTENT(IN) :: PDXX(D%NIJT, D%NKT)
    REAL, INTENT(IN) :: PDYY(D%NIJT, D%NKT)
    REAL, INTENT(IN) :: PTHVREF(D%NIJT, D%NKT)
    REAL, INTENT(IN) :: PSRCT(MERGE(D%NIJT, 0, OCOMPUTE_SRC), MERGE(D%NKT, 0, OCOMPUTE_SRC))
    REAL, INTENT(IN) :: PRT(D%NIJT, D%NKT, KRR)
    REAL, INTENT(IN) :: PTHLT(D%NIJT, D%NKT)
    REAL, INTENT(IN) :: PLOCPEXNM(D%NIJT, D%NKT)
    REAL, INTENT(INOUT) :: PATHETA(D%NIJT, D%NKT)
    REAL, INTENT(INOUT) :: PAMOIST(D%NIJT, D%NKT)
    REAL, INTENT(OUT), DIMENSION(D%NIJT, D%NKT) :: PLM
    REAL :: ZALPHA,ZD,ZVAR
    LOGICAL :: GZD
    REAL :: PWORK1(D%NIJT, D%NKT)
    REAL :: PWORK2(D%NIJT, D%NKT)
    REAL :: PWORK2D(D%NIJT)
    REAL :: ZETHETA(D%NIJT, D%NKT)
    REAL :: ZEMOIST(D%NIJT, D%NKT)
    REAL :: ZDRTDZ(D%NIJT, D%NKT) !drt_dz used for computing the stablity criterion
    REAL :: ZDTHLDZ(D%NIJT, D%NKT) !dtheta_l/dz used for computing the stablity criterion

    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE2
    INTEGER :: IKB, IKE, IKU, IKA, IKT, IKL
    INTEGER :: IIJB, IKTB, IIJE, IKTE
    INTEGER :: JK, JIJ
    !-------------------------------------------------------------------------------
    !
    !   initialize the mixing length with the mesh grid
    IF (LHOOK) CALL DR_HOOK('TURB:DEAR', 0, ZHOOK_HANDLE2)
    IKT = D%NKT
    IKTB = D%NKTB
    IKTE = D%NKTE
    IKB = D%NKB
    IKE = D%NKE
    IKA = D%NKA
    IKU = D%NKU
    IKL = D%NKL
    IIJE = D%NIJE
    IIJB = D%NIJB
    IF (TURBN%CTURBDIM /= '1DIM') THEN
      CALL MXF_PHY(D, PDXX, PWORK1)
      IF (.not.O2D) THEN
        CALL MYF_PHY(D, PDYY, PWORK2)
      END IF
    END IF
    ! 1D turbulence scheme

DO JK=IKTB, IKTE
      DO JIJ=IIJB, IIJE
        PLM(JIJ, JK) = PZZ(JIJ, IKL + JK) - PZZ(JIJ, JK)
  END DO
    END DO
DO JIJ=IIJB, IIJE
      PLM(JIJ, IKU) = PLM(JIJ, IKE)
      PLM(JIJ, IKA) = PZZ(JIJ, IKB) - PZZ(JIJ, IKA)
END DO

    !
    IF (TURBN%CTURBDIM /= '1DIM') THEN
      ! 3D turbulence scheme
      IF (O2D) THEN

DO JK=1, IKT
          DO JIJ=IIJB, IIJE
            PLM(JIJ, JK) = SQRT(PLM(JIJ, JK)*PWORK1(JIJ, JK))
  END DO
        END DO

      ELSE

DO JK=1, IKT
          DO JIJ=IIJB, IIJE
            PLM(JIJ, JK) = (PLM(JIJ, JK)*PWORK1(JIJ, JK)*PWORK2(JIJ, JK))**(1. / 3.)
  END DO
        END DO

      END IF
    END IF
    !   compute a mixing length limited by the stability
    !
    CALL ETHETA(D, CST, KRR, KRRI, PTHLT, PRT, PLOCPEXNM, PATHETA, PSRCT, OOCEAN, OCOMPUTE_SRC, ZETHETA)
    CALL EMOIST(D, CST, KRR, KRRI, PTHLT, PRT, PLOCPEXNM, PAMOIST, PSRCT, OOCEAN, ZEMOIST)
    !
    IF (KRR > 0) THEN


      DO JK=IKTB + 1,IKTE - 1
        DO JIJ=IIJB,IIJE
          ZDTHLDZ(JIJ, JK) = 0.5*((PTHLT(JIJ, JK + IKL) - PTHLT(JIJ, JK)) / PDZZ(JIJ, JK + IKL) + (PTHLT(JIJ, JK) - &
                             PTHLT(JIJ,JK-IKL)) / PDZZ(JIJ, JK))
          ZDRTDZ(JIJ, JK) = 0.5*((PRT(JIJ, JK + IKL, 1) - PRT(JIJ, JK, 1)) / PDZZ(JIJ, JK + IKL) + (PRT(JIJ, JK, 1) - &
                             PRT(JIJ,JK-IKL, 1)) / PDZZ(JIJ, JK))
        END DO
      END DO




      DO JK=IKTB + 1,IKTE - 1
        DO JIJ=IIJB,IIJE
          IF (OOCEAN) THEN
            ZVAR = CST%XG*(CST%XALPHAOC*ZDTHLDZ(JIJ, JK) - CST%XBETAOC*ZDRTDZ(JIJ, JK))
          ELSE
            ZVAR = CST%XG / PTHVREF(JIJ, JK)*(ZETHETA(JIJ, JK)*ZDTHLDZ(JIJ, JK) + ZEMOIST(JIJ, JK)*ZDRTDZ(JIJ, JK))
          END IF
          !
          IF (ZVAR > 0.) THEN
            PLM(JIJ, JK) = MAX(CST%XMNH_EPSILON, MIN(PLM(JIJ, JK), 0.76*SQRT(PTKET(JIJ, JK) / ZVAR)))
          END IF
        END DO
      END DO


    ELSE
      ! For dry atmos or unsalted ocean runs


      DO JK=IKTB + 1,IKTE - 1
        DO JIJ=IIJB,IIJE
          ZDTHLDZ(JIJ, JK) = 0.5*((PTHLT(JIJ, JK + IKL) - PTHLT(JIJ, JK)) / PDZZ(JIJ, JK + IKL) + (PTHLT(JIJ, JK) - &
                             PTHLT(JIJ,JK-IKL)) / PDZZ(JIJ, JK))
          IF (OOCEAN) THEN
            ZVAR = CST%XG*CST%XALPHAOC*ZDTHLDZ(JIJ, JK)
          ELSE
            ZVAR = CST%XG / PTHVREF(JIJ, JK)*ZETHETA(JIJ, JK)*ZDTHLDZ(JIJ, JK)
          END IF
          !
          IF (ZVAR > 0.) THEN
            PLM(JIJ, JK) = MAX(CST%XMNH_EPSILON, MIN(PLM(JIJ, JK), 0.76*SQRT(PTKET(JIJ, JK) / ZVAR)))
          END IF
        END DO
      END DO

    END IF

    !  special case near the surface
DO JIJ=IIJB, IIJE
      ZDTHLDZ(JIJ, IKB) = (PTHLT(JIJ, IKB + IKL) - PTHLT(JIJ, IKB)) / PDZZ(JIJ, IKB + IKL)
END DO
    ! For dry simulations
    IF (KRR > 0) THEN
DO JIJ=IIJB, IIJE
        ZDRTDZ(JIJ, IKB) = (PRT(JIJ, IKB + IKL, 1) - PRT(JIJ, IKB, 1)) / PDZZ(JIJ, IKB + IKL)
END DO
    ELSE
      ZDRTDZ(:, IKB) = 0
    END IF
    !
    IF (OOCEAN) THEN
DO JIJ=IIJB, IIJE
        PWORK2D(JIJ) = CST%XG*(CST%XALPHAOC*ZDTHLDZ(JIJ, IKB) - CST%XBETAOC*ZDRTDZ(JIJ, IKB))
END DO
    ELSE
DO JIJ=IIJB, IIJE
        PWORK2D(JIJ) = CST%XG / PTHVREF(JIJ, IKB)*(ZETHETA(JIJ, IKB)*ZDTHLDZ(JIJ, IKB) + &
                             ZEMOIST(JIJ, IKB)*ZDRTDZ(JIJ, IKB))
END DO
    END IF
DO JIJ=IIJB, IIJE
      IF (PWORK2D(JIJ) > 0.) THEN
        PLM(JIJ, IKB) =  &
        & MAX(CST%XMNH_EPSILON, MIN(PLM(JIJ, IKB), 0.76*SQRT(PTKET(JIJ, IKB) / PWORK2D(JIJ))))
      END IF
END DO
    !
    !  mixing length limited by the distance normal to the surface (with the same factor as for BL89)
    !
    IF (.NOT. TURBN%LRMC01) THEN
      ZALPHA = 0.5**(-1.5)
      !

      DO JIJ=IIJB,IIJE
        GZD = .TRUE.
        IF (OOCEAN) THEN

          DO JK=IKTE,IKTB,-1
            ZD = ZALPHA*(PZZ(JIJ, IKTE + 1) - PZZ(JIJ, JK))
            IF (PLM(JIJ, JK) > ZD .AND. GZD) THEN
              PLM(JIJ, JK) = ZD
            ELSE
              GZD = .FALSE.
            END IF
          END DO
        ELSE
          DO JK=IKTB,IKTE
            ZD = ZALPHA*(0.5*(PZZ(JIJ, JK) + PZZ(JIJ, JK + IKL)) - PZZ(JIJ, IKB))*PDIRCOSZW(JIJ)
            IF (PLM(JIJ, JK) > ZD .AND. GZD) THEN
              PLM(JIJ, JK) = ZD
            ELSE
              GZD = .FALSE.
            END IF
          END DO
        END IF
      END DO
    END IF
    !
DO JIJ=IIJB, IIJE
      PLM(JIJ, IKA) = PLM(JIJ, IKB)
      PLM(JIJ, IKE) = PLM(JIJ, IKE - IKL)
      PLM(JIJ, IKU) = PLM(JIJ, IKU - IKL)
END DO
    !

    IF (LHOOK) CALL DR_HOOK('TURB:DEAR', 1, ZHOOK_HANDLE2)
END SUBROUTINE DEAR
END MODULE MODE_DEAR
