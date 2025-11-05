!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_DELT
IMPLICIT NONE
CONTAINS
  SUBROUTINE DELT (D,TURBN, O2D, ODZ, OOCEAN, PZZ, PDYY, PDXX, PDIRCOSZW, PLM)
    !     ####################
    !!
    !!****  *DELT* routine to compute mixing length for DELT case
    !
    !!    AUTHOR
    !!    ------
    !!
    !!     M Tomasini      *Meteo-France
    !!
    !!    MODIFICATIONS
    !!    -------------
    !!      Original   01/05
    !!
    !-------------------------------------------------------------------------------
    !
    !*       0.    DECLARATIONS
    !              ------------
    !
    USE YOMHOOK, ONLY: LHOOK, DR_HOOK, JPHOOK
    USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
    USE MODD_TURB_n, ONLY: TURB_t
    USE MODE_SHUMAN_PHY, ONLY: MXF_PHY,MYF_PHY

    IMPLICIT NONE
    !*       0.1   Declarations of dummy arguments
    !
    TYPE(DIMPHYEX_t), INTENT(IN) :: D
    TYPE(TURB_t), INTENT(IN) :: TURBN
    REAL, INTENT(IN) :: PDXX(D%NIJT, D%NKT)
    LOGICAL, INTENT(IN) :: ODZ
    LOGICAL, INTENT(IN) :: O2D
    LOGICAL, INTENT(IN) :: OOCEAN
    REAL, INTENT(IN) :: PZZ(D%NIJT, D%NKT)
    REAL, INTENT(IN) :: PDYY(D%NIJT, D%NKT)
    REAL, INTENT(IN) :: PDIRCOSZW(D%NIJT)
    REAL, INTENT(OUT), DIMENSION(D%NIJT, D%NKT) :: PLM
    REAL :: ZWORK2(D%NIJT, D%NKT)
    REAL :: ZWORK1(D%NIJT, D%NKT)
    REAL :: ZALPHA
    REAL :: ZD

    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE2
    INTEGER :: IKB, IKE, IKU, IKA, IKT, IKL
    INTEGER :: IIJB, IKTB, IIJE, IKTE
    INTEGER :: JK, JIJ
    !-------------------------------------------------------------------------------
    !
    IF (LHOOK) CALL DR_HOOK('TURB:DELT', 0, ZHOOK_HANDLE2)
    !
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
    CALL MXF_PHY(D, PDXX, ZWORK1)
    IF (.NOT. O2D) THEN
      CALL MYF_PHY(D, PDYY, ZWORK2)
    END IF
    !
    IF (ODZ) THEN
!$acc kernels present_cr( PLM )
      ! Dz is take into account in the computation
      DO JK=IKTB,IKTE
        ! 1D turbulence scheme
!$mnh_expand_array ( JIJ=IIJB:IIJE )
        PLM(IIJB:IIJE, JK) = PZZ(IIJB:IIJE, JK + IKL) - PZZ(IIJB:IIJE, JK)
!$mnh_end_expand_array ( JIJ=IIJB:IIJE )
      END DO
!$mnh_expand_array ( JIJ=IIJB:IIJE )
      PLM(IIJB:IIJE, IKU) = PLM(IIJB:IIJE, IKE)
      PLM(IIJB:IIJE, IKA) = PZZ(IIJB:IIJE, IKB) - PZZ(IIJB:IIJE, IKA)
!$mnh_end_expand_array ( JIJ=IIJB:IIJE )
!$acc end kernels
      IF (TURBN%CTURBDIM /= '1DIM') THEN
        ! 3D turbulence scheme
        IF (O2D) THEN
!$acc kernels present_cr( PLM )
!$mnh_expand_array ( JIJ=IIJB:IIJE,JK=1:IKT )
          PLM(IIJB:IIJE, 1:IKT) = SQRT(PLM(IIJB:IIJE, 1:IKT)*ZWORK1(IIJB:IIJE, 1:IKT))
!$mnh_end_expand_array ( JIJ=IIJB:IIJE,JK=1:IKT )
!$acc end kernels
        ELSE
!$acc kernels present_cr( PLM )
!$mnh_expand_array ( JIJ=IIJB:IIJE,JK=1:IKT )
          PLM(IIJB:IIJE, 1:IKT) = (PLM(IIJB:IIJE, 1:IKT)*ZWORK1(IIJB:IIJE, 1:IKT)*ZWORK2(IIJB:IIJE, 1:IKT))**(1. / 3.)
!$mnh_end_expand_array ( JIJ=IIJB:IIJE,JK=1:IKT )
!$acc end kernels
        END IF
      END IF
    ELSE
      ! Dz not taken into account in computation to assure invariability with vertical grid mesh
!$acc kernels present_cr( PLM )
!$mnh_expand_array ( JIJ=IIJB:IIJE,JK=1:IKT )
      PLM(IIJB:IIJE, 1:IKT) = 1.E10
!$mnh_end_expand_array ( JIJ=IIJB:IIJE,JK=1:IKT )
!$acc end kernels
      IF (TURBN%CTURBDIM /= '1DIM') THEN
        ! 3D turbulence scheme
        IF (O2D) THEN
!$acc kernels present_cr( PLM )
!$mnh_expand_array ( JIJ=IIJB:IIJE,JK=1:IKT )
          PLM(:, :) = ZWORK1(:, :)
!$mnh_end_expand_array ( JIJ=IIJB:IIJE,JK=1:IKT )
!$acc end kernels
        ELSE
!$acc kernels present_cr( PLM )
!$mnh_expand_array ( JIJ=IIJB:IIJE,JK=1:IKT )
          PLM(IIJB:IIJE, 1:IKT) = (ZWORK1(IIJB:IIJE, 1:IKT)*ZWORK2(IIJB:IIJE, 1:IKT))**(1. / 2.)
!$mnh_end_expand_array ( JIJ=IIJB:IIJE,JK=1:IKT )
!$acc end kernels
        END IF
      END IF
    END IF
    !
    !  mixing length limited by the distance normal to the surface
    !  (with the same factor as for BL89)
    !
    IF (.NOT. TURBN%LRMC01) THEN
      ZALPHA = 0.5**(-1.5)
      !
!$acc kernels
      DO JIJ=IIJB,IIJE
        IF (OOCEAN) THEN
          DO JK=IKTE,IKTB,-1
            ZD = ZALPHA*(PZZ(JIJ, IKTE + 1) - PZZ(JIJ, JK))
            IF (PLM(JIJ, JK) > ZD) THEN
              PLM(JIJ, JK) = ZD
            ELSE
              EXIT
            END IF
          END DO
        ELSE
          DO JK=IKTB,IKTE
            ZD = ZALPHA*(0.5*(PZZ(JIJ, JK) + PZZ(JIJ, JK + IKL)) - PZZ(JIJ, IKB))*PDIRCOSZW(JIJ)
            IF (PLM(JIJ, JK) > ZD) THEN
              PLM(JIJ, JK) = ZD
            ELSE
              EXIT
            END IF
          END DO
        END IF
      END DO
!$acc end kernels
    END IF
    !
!$acc kernels
!$mnh_expand_array ( JIJ=IIJB:IIJE )
    PLM(IIJB:IIJE, IKA) = PLM(IIJB:IIJE, IKB)
    PLM(IIJB:IIJE, IKU) = PLM(IIJB:IIJE, IKE)
!$mnh_end_expand_array ( JIJ=IIJB:IIJE )
!$acc end kernels
    !
    IF (LHOOK) CALL DR_HOOK('TURB:DELT', 1, ZHOOK_HANDLE2)
  END SUBROUTINE DELT
END MODULE MODE_DELT
