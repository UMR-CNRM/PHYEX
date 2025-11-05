!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_CLOUD_MODIF_LM
IMPLICIT NONE
CONTAINS
  SUBROUTINE CLOUD_MODIF_LM (D, CST, CSTURB, TURBN, TPFILE, TZFIELD, KRR, KRRI, &
                             & OCLOUDMODIFLM, OOCEAN, OCOMPUTE_SRC, O2D, HTURBLEN_CL, &
                             & PDZZ, PDXX, PDYY, PZZ, &
                             & PRT, PTKET, PTHLT, PTHLM, PRM, PTHVREF, &
                             & PLOCPEXNM, PSRCT, PCOEF_AMPL_SAT, PAMOIST, PATHETA, PDIRCOSZW, &
                             & PCEI, PCEI_MIN, PCEI_MAX, PLM)
    !     #########################
    !!
    !!*****CLOUD_MODIF_LM routine to:
    !!       1/ change the mixing length in the clouds
    !!       2/ emphasize the mixing length in the cloud
    !!           by the coefficient ZCOEF_AMPL calculated here
    !!             when the CEI index is above ZCEI_MIN.
    !!
    !!
    !!      ZCOEF_AMPL ^
    !!                 |
    !!                 |
    !!  PCOEF_AMPL_SAT -                       ---------- Saturation
    !!    (XDUMMY1)    |                      -
    !!                 |                     -
    !!                 |                    -
    !!                 |                   -
    !!                 |                  - Amplification
    !!                 |                 - straight
    !!                 |                - line
    !!                 |               -
    !!                 |              -
    !!                 |             -
    !!                 |            -
    !!                 |           -
    !!               1 ------------
    !!                 |
    !!                 |
    !!               0 -----------|------------|----------> PCEI
    !!                 0      ZCEI_MIN     ZCEI_MAX
    !!                        (XDUMMY2)    (XDUMMY3)
    !!
    !!
    !!
    !!    AUTHOR
    !!    ------
    !!     M. Tomasini   *CNRM METEO-FRANCE
    !!
    !!    MODIFICATIONS
    !!    -------------
    !!     Original   09/07/04
    !!
    !-------------------------------------------------------------------------------
    !
    !*       0.    DECLARATIONS
    !              ------------
    USE YOMHOOK, ONLY: LHOOK, DR_HOOK, JPHOOK
    USE MODD_CST, ONLY: CST_t
    USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
    USE MODD_IO, ONLY: TFILEDATA
    USE MODD_FIELD, ONLY: TFIELDMETADATA, TYPEREAL
    USE MODD_TURB_n, ONLY: TURB_t
    USE MODD_CTURB, ONLY: CSTURB_t

    USE MODE_DELT, ONLY: DELT
    USE MODE_DEAR, ONLY: DEAR
    USE MODE_BL89, ONLY: BL89
    USE MODE_IO_FIELD_WRITE_PHY,      ONLY: IO_FIELD_WRITE_PHY
    !
    IMPLICIT NONE
    !
    !*       0.1   Declarations of dummy arguments
    !
    TYPE(DIMPHYEX_t), INTENT(IN) :: D
    TYPE(CST_t), INTENT(IN) :: CST
    TYPE(CSTURB_t), INTENT(IN) :: CSTURB
    TYPE(TURB_t), INTENT(IN) :: TURBN
    TYPE(TFIELDMETADATA), INTENT(INOUT) :: TZFIELD
    TYPE(TFILEDATA), INTENT(INOUT) :: TPFILE
    INTEGER, INTENT(IN) :: KRRI, KRR
    LOGICAL, INTENT(IN) :: OCLOUDMODIFLM  
    LOGICAL, INTENT(IN) :: OOCEAN
    LOGICAL, INTENT(IN) :: O2D
    LOGICAL, INTENT(IN) :: OCOMPUTE_SRC
    CHARACTER(LEN=4), INTENT(IN) :: HTURBLEN_CL
    REAL, INTENT(IN) :: PTKET(D%NIJT, D%NKT)
    REAL, INTENT(IN) :: PDZZ(D%NIJT, D%NKT)
    REAL, INTENT(IN) :: PZZ(D%NIJT, D%NKT)
    REAL, INTENT(IN) :: PTHVREF(D%NIJT, D%NKT)
    REAL, INTENT(IN) :: PSRCT(MERGE(D%NIJT, 0, OCOMPUTE_SRC), MERGE(D%NKT, 0, OCOMPUTE_SRC))
    REAL, INTENT(IN) :: PCOEF_AMPL_SAT
    REAL, INTENT(IN) :: PDIRCOSZW(D%NIJT)
    REAL, INTENT(IN) :: PCEI(MERGE(D%NIJT, 0, OCLOUDMODIFLM), MERGE(D%NKT, 0, OCLOUDMODIFLM))
    REAL, INTENT(IN) :: PCEI_MAX
    REAL, INTENT(IN) :: PCEI_MIN
    REAL, INTENT(IN) :: PDXX(D%NIJT, D%NKT)
    REAL, INTENT(IN) :: PDYY(D%NIJT, D%NKT)
    REAL, INTENT(IN) :: PTHLT(D%NIJT, D%NKT)
    REAL, INTENT(IN) :: PTHLM(D%NIJT, D%NKT)
    REAL, INTENT(IN) :: PRM(D%NIJT, D%NKT, KRR)
    REAL, INTENT(IN) :: PLOCPEXNM(D%NIJT, D%NKT)
    REAL, INTENT(IN) :: PRT(D%NIJT, D%NKT, KRR)
    REAL, INTENT(INOUT) :: PLM(D%NIJT, D%NKT)
    REAL, INTENT(INOUT) :: PAMOIST(D%NIJT, D%NKT)
    REAL, INTENT(INOUT) :: PATHETA(D%NIJT, D%NKT)
    REAL :: ZCOEF_AMPL(D%NIJT, D%NKT)
    REAL :: ZSHEAR(D%NIJT, D%NKT)
    REAL :: PCOEF_AMPL_CEI_NUL
    REAL :: ZPENTE
    REAL :: PLM_CLOUD(D%NIJT, D%NKT) ! Turbulent mixing length in the clouds
    INTEGER :: IKE, IKU, IKA, IKT, IKL
    INTEGER :: IIJB, IKTB, IIJE, IKTE
    INTEGER :: JK, JIJ
    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE2
    !
    !-------------------------------------------------------------------------------
    !
    !*       1.    INITIALISATION
    !              --------------
    !
    IF (LHOOK) CALL DR_HOOK('TURB:CLOUD_MODIF_LM', 0, ZHOOK_HANDLE2)
    IKT = D%NKT
    IKTB = D%NKTB
    IKTE = D%NKTE
    IKE = D%NKE
    IKA = D%NKA
    IKU = D%NKU
    IKL = D%NKL
    IIJE = D%NIJE
    IIJB = D%NIJB
    ZPENTE = (PCOEF_AMPL_SAT - 1.) / (PCEI_MAX - PCEI_MIN)
    PCOEF_AMPL_CEI_NUL = 1. - ZPENTE*PCEI_MIN
    !
!$acc kernels
!$mnh_expand_array( JIJ=IIJB:IIJE,JK=1:IKT )
    ZCOEF_AMPL(IIJB:IIJE, 1:IKT) = 1.
!$mnh_end_expand_array ( JIJ=IIJB:IIJE,JK=1:IKT )
!$acc end kernels
    !
    !*       2.    CALCULATION OF THE AMPLIFICATION COEFFICIENT
    !              --------------------------------------------
    !
    ! Saturation
    !
!$acc kernels
!$mnh_expand_where( JIJ=IIJB:IIJE,JK=1:IKT )
    WHERE (PCEI(IIJB:IIJE, 1:IKT) >= PCEI_MAX)
      ZCOEF_AMPL(IIJB:IIJE, 1:IKT) = PCOEF_AMPL_SAT
    END WHERE
!$mnh_end_expand_where ( JIJ=IIJB:IIJE,JK=1:IKT )
!$acc end kernels
    !
    ! Between the min and max limits of CEI index, linear variation of the
    ! amplification coefficient ZCOEF_AMPL as a function of CEI
    !
!$acc kernels
!$mnh_expand_where( JIJ=IIJB:IIJE,JK=1:IKT )
    WHERE (PCEI(IIJB:IIJE, 1:IKT) < PCEI_MAX .and. PCEI(IIJB:IIJE, 1:IKT) > PCEI_MIN)
      ZCOEF_AMPL(IIJB:IIJE, 1:IKT) = ZPENTE*PCEI(IIJB:IIJE, 1:IKT) + PCOEF_AMPL_CEI_NUL
    END WHERE
!$mnh_end_expand_where ( JIJ=IIJB:IIJE,JK=1:IKT )
!$acc end kernels
    !
    !
    !*       3.    CALCULATION OF THE MIXING LENGTH IN CLOUDS
    !              ------------------------------------------
    !
    IF (HTURBLEN_CL == TURBN%CTURBLEN) THEN
!$acc kernels
!$mnh_expand_array( JIJ=IIJB:IIJE,JK=1:IKT )
      PLM_CLOUD(:, :) = PLM(:, :)
!$mnh_end_expand_array ( JIJ=IIJB:IIJE,JK=1:IKT )
!$acc end kernels
    ELSE
      !
      !*         3.1 BL89 mixing length
      !           ------------------
      SELECT CASE (HTURBLEN_CL)
      CASE ('BL89', 'RM17', 'HM21')
!$acc kernels
!$mnh_expand_array( JIJ=IIJB:IIJE,JK=1:IKT )
        ZSHEAR(:, :) = 0.
!$mnh_end_expand_array ( JIJ=IIJB:IIJE,JK=1:IKT )
!$acc end kernels
        CALL BL89(D, CST, CSTURB, TURBN, PZZ, PDZZ, PTHVREF, PTHLM, KRR, PRM, PTKET, ZSHEAR, PLM_CLOUD, OOCEAN)
        !
        !*         3.2 Delta mixing length
        !           -------------------
      CASE ('DELT')
        CALL DELT (D,TURBN, O2D, .TRUE., OOCEAN, PZZ, PDYY, PDXX, PDIRCOSZW, PLM_CLOUD)
        !
        !*         3.3 Deardorff mixing length
        !           -----------------------
      CASE ('DEAR')
        CALL DEAR(D, CST, TURBN, KRR, KRRI, O2D, OCOMPUTE_SRC, OOCEAN, &
                  & PLM_CLOUD, PRT, PDZZ, PZZ, PTKET, PTHVREF, &
                  & PTHLT, PLOCPEXNM, PSRCT, PAMOIST, PDIRCOSZW, &
                  & PDXX, PDYY, PATHETA)
        !
      END SELECT
    END IF
    !
    !*       4.    MODIFICATION OF THE MIXING LENGTH IN THE CLOUDS
    !              -----------------------------------------------
    !
    ! Impression before modification of the mixing length
IF ( TURBN%LTURB_DIAG .AND. TPFILE%LOPENED ) THEN
  TZFIELD = TFIELDMETADATA(            &
    CMNHNAME   = 'LM_CLEAR_SKY',       &
    CSTDNAME   = '',                   &
    CLONGNAME  = 'LM_CLEAR_SKY',       &
    CUNITS     = 'm',                  &
    CDIR       = 'XY',                 &
    CCOMMENT   = 'X_Y_Z_LM CLEAR SKY', &
    NGRID      = 1,                    &
    NTYPE      = TYPEREAL,             &
    NDIMS      = 3,                    &
    LTIMEDEP   = .TRUE.                )
!$acc update self(PLM)
  CALL IO_FIELD_WRITE_PHY(D,TPFILE,TZFIELD,PLM)
ENDIF
    !
    ! Amplification of the mixing length when the criteria are verified
    !
!$acc kernels
!$mnh_expand_where( JIJ=IIJB:IIJE,JK=1:IKT )
    WHERE (ZCOEF_AMPL(IIJB:IIJE, 1:IKT) /= 1.)
      PLM(IIJB:IIJE, 1:IKT) = ZCOEF_AMPL(IIJB:IIJE, 1:IKT)*PLM_CLOUD(IIJB:IIJE, 1:IKT)
    END WHERE
!$mnh_end_expand_where ( JIJ=IIJB:IIJE,JK=1:IKT )
!$acc end kernels
    !
    ! Cloud mixing length in the clouds at the points which do not verified the CEI
    !
!$acc kernels
!$mnh_expand_where( JIJ=IIJB:IIJE,JK=1:IKT )
    WHERE (PCEI(IIJB:IIJE, 1:IKT) == -1.)
      PLM(IIJB:IIJE, 1:IKT) = PLM_CLOUD(IIJB:IIJE, 1:IKT)
    END WHERE
!$mnh_end_expand_where ( JIJ=IIJB:IIJE,JK=1:IKT )
!$acc end kernels
    !
    !
    !*       5.    IMPRESSION
    !              ----------
    !
IF ( TURBN%LTURB_DIAG .AND. TPFILE%LOPENED ) THEN
  TZFIELD = TFIELDMETADATA(         &
    CMNHNAME   = 'COEF_AMPL',       &
    CSTDNAME   = '',                &
    CLONGNAME  = 'COEF_AMPL',       &
    CUNITS     = '1',               &
    CDIR       = 'XY',              &
    CCOMMENT   = 'X_Y_Z_COEF AMPL', &
    NGRID      = 1,                 &
    NTYPE      = TYPEREAL,          &
    NDIMS      = 3,                 &
    LTIMEDEP   = .TRUE.             )
!$acc update self(ZCOEF_AMPL)
  CALL IO_FIELD_WRITE_PHY(D,TPFILE,TZFIELD,ZCOEF_AMPL)
  !
  TZFIELD = TFIELDMETADATA(        &
    CMNHNAME   = 'LM_CLOUD',       &
    CSTDNAME   = '',               &
    CLONGNAME  = 'LM_CLOUD',       &
    CUNITS     = 'm',              &
    CDIR       = 'XY',             &
    CCOMMENT   = 'X_Y_Z_LM CLOUD', &
    NGRID      = 1,                &
    NTYPE      = TYPEREAL,         &
    NDIMS      = 3,                &
    LTIMEDEP   = .TRUE.            )
!$acc update self( PLM_CLOUD )
   CALL IO_FIELD_WRITE_PHY(D, TPFILE, TZFIELD, PLM_CLOUD)
   !
END IF
    !
    IF (LHOOK) CALL DR_HOOK('TURB:CLOUD_MODIF_LM', 1, ZHOOK_HANDLE2)
  END SUBROUTINE CLOUD_MODIF_LM
END MODULE MODE_CLOUD_MODIF_LM
