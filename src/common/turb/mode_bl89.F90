!MNH_LIC Copyright 1994-2022 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
MODULE MODE_BL89
IMPLICIT NONE
CONTAINS
!     ######spl
      SUBROUTINE BL89(D,CST,CSTURB,PZZ,PDZZ,PTHVREF,PTHLM,KRR,PRM,PTKEM,PSHEAR,PLM,OOCEAN,HPROGRAM)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!     #########################################################
!
!!****  *BL89* -
!!
!!    PURPOSE
!!    -------
!!    This routine computes the mixing length from Bougeault-Lacarrere 89
!!    formula.
!!
!!**  METHOD
!!    ------
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
!!      Book 2
!!
!!    AUTHOR
!!    ------
!!
!!      J. Cuxart  INM and Meteo-France
!!
!!    MODIFICATIONS
!!    -------------
!!     Original     27/04/97 (V. Masson) separation from turb.f90
!!                                       and optimization
!!                  06/01/98 (V. Masson and P. Jabouille) optimization
!!                  15/03/99 (V. Masson) new lup ldown averaging
!!                  21/02/01 (P. Jabouille) improve vectorization
!!                  2012-02 (Y. Seity) add possibility to run with
!!                            reversed vertical levels
!!  Philippe 13/02/2018: use ifdef MNH_REAL to prevent problems with intrinsics on Blue Gene/Q
!!                  01/2019 (Q. Rodier) support for RM17 mixing length
!!                  03/2021 (JL Redelsperger) Ocean model case 
!!                  06/2021 (P. Marquet) correction of exponent on final length according to Lemarié et al. 2021
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST, ONLY: CST_t
USE MODD_CTURB, ONLY: CSTURB_t
USE MODD_DIMPHYEX,   ONLY: DIMPHYEX_t
USE MODD_PARAMETERS, ONLY: JPVEXT_TURB
USE MODD_PRECISION, ONLY: MNHREAL
!
!
IMPLICIT NONE
!
!*       0.1   Declaration of arguments
!              ------------------------
!
TYPE(DIMPHYEX_t),         INTENT(IN)  :: D
TYPE(CST_t),              INTENT(IN)  :: CST
TYPE(CSTURB_t),           INTENT(IN)  :: CSTURB
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN),TARGET  :: PZZ
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN),TARGET  :: PDZZ
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN),TARGET  :: PTHVREF
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN),TARGET  :: PTHLM       ! conservative pot. temp.
INTEGER,                  INTENT(IN)  :: KRR
REAL, DIMENSION(D%NIT,D%NJT,D%NKT,KRR), INTENT(IN),TARGET :: PRM       ! water var.
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN),TARGET  :: PTKEM     ! TKE
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN),TARGET  :: PSHEAR
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(OUT),TARGET :: PLM       ! Mixing length
LOGICAL,                  INTENT(IN)  ::  OOCEAN       ! switch for Ocean model version
CHARACTER(LEN=6),         INTENT(IN)  ::  HPROGRAM     ! CPROGRAM is the program currently running (modd_conf)
!   thermodynamical variables PTHLM=Theta at the begining
!
!*       0.2   Declaration of local variables
!              ------------------------------
!
INTEGER :: IKB,IKE
INTEGER :: IKTB,IKTE    ! start, end of k loops in physical domain

REAL, DIMENSION(D%NIT*D%NJT,D%NKT) :: ZVPT  ! Virtual Potential Temp at half levels
REAL, DIMENSION(D%NIT*D%NJT,D%NKT) :: ZDELTVPT
            ! Increment of Virtual Potential Temp between two following levels
REAL, DIMENSION(D%NIT*D%NJT,D%NKT) :: ZHLVPT
            ! Virtual Potential Temp at half levels
REAL, DIMENSION(D%NIT*D%NJT) ::  ZLWORK,ZINTE
!           ! downwards then upwards vertical displacement,
!           ! residual internal energy,
!           ! residual potential energy
REAL, POINTER, DIMENSION(:,:) :: ZZZ,ZDZZ,       &
                                 ZTHM,ZTKEM,ZLM, &
                                 ZSHEAR,ZTHVREF
!           ! input and output arrays packed according one horizontal coord.
REAL, POINTER, DIMENSION(:,:,:)  :: ZRM
!           ! input array packed according one horizontal coord.
REAL, DIMENSION(D%NIT*D%NJT,D%NKT) :: ZSUM ! to replace SUM function
REAL, DIMENSION(D%NIT*D%NJT,D%NKT) :: ZG_O_THVREF
REAL, DIMENSION(D%NIT*D%NJT,D%NKT) :: ZSQRT_TKE
REAL, DIMENSION(D%NIT*D%NJT,D%NKT) :: ZLMDN
!
INTEGER :: IIU,IJU,IPROMA
INTEGER :: J1D        ! horizontal loop counter
INTEGER :: JK,JKK     ! loop counters
INTEGER :: JRR        ! moist loop counter
REAL    :: ZRVORD     ! Rv/Rd
REAL    :: ZPOTE,ZLWORK1,ZLWORK2
REAL    :: ZTEST,ZTEST0,ZTESTM ! test for vectorization
REAL    :: Z2SQRT2,ZUSRBL89,ZBL89EXP
!-------------------------------------------------------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('BL89',0,ZHOOK_HANDLE)
Z2SQRT2=2.*SQRT(2.)
!
ZRVORD = CST%XRV / CST%XRD
IPROMA = D%NIT*D%NJT
!
!-------------------------------------------------------------------------------
!
!*       1.    pack the horizontal dimensions into one
!              ---------------------------------------
!
! Pointer remapping instead of RESHAPE (contiguous memory)
! 2D array => 3D array
ZZZ(1:IPROMA,1:D%NKT) => PZZ
ZDZZ(1:IPROMA,1:D%NKT) => PDZZ
ZTHM(1:IPROMA,1:D%NKT) => PTHLM
ZSHEAR(1:IPROMA,1:D%NKT) => PSHEAR
ZTKEM(1:IPROMA,1:D%NKT) => PTKEM
ZTHVREF(1:IPROMA,1:D%NKT) => PTHVREF
ZRM(1:IPROMA,1:D%NKT,1:KRR) =>PRM
ZLM(1:IPROMA,1:D%NKT) =>PLM
!
IF (OOCEAN) THEN
  DO JK=1,D%NKT
    DO J1D=1,IPROMA
      ZG_O_THVREF(J1D,JK) = CST%XG * CST%XALPHAOC
    END DO
  END DO
ELSE !Atmosphere case
  DO JK=1,D%NKT
    DO J1D=1,IPROMA
      ZG_O_THVREF(J1D,JK) = CST%XG / ZTHVREF(J1D,JK)
    END DO
  END DO
END IF
!
!$mnh_expand_array(J1D=1:IPROMA,JK=1:D%NKT)
ZSQRT_TKE(:,:) = SQRT(ZTKEM(:,:))
!$mnh_end_expand_array(J1D=1:IPROMA,JK=1:D%NKT)
!
!ZBL89EXP is defined here because (and not in ini_cturb) because CSTURB%XCED is defined in read_exseg (depending on BL89/RM17)
ZBL89EXP = LOG(16.)/(4.*LOG(CST%XKARMAN)+LOG(CSTURB%XCED)-3.*LOG(CSTURB%XCMFS))
ZUSRBL89 = 1./ZBL89EXP
!-------------------------------------------------------------------------------
!
!*       2.    Virtual potential temperature on the model grid
!              -----------------------------------------------
!
IF(KRR /= 0) THEN
  ZSUM(:,:) = 0.
  DO JRR=1,KRR
    !$mnh_expand_array(J1D=1:IPROMA,JK=1:D%NKT)
    ZSUM(:,:) = ZSUM(:,:)+ZRM(:,:,JRR)
    !$mnh_end_expand_array(J1D=1:IPROMA,JK=1:D%NKT)
  ENDDO
  !$mnh_expand_array(J1D=1:IPROMA,JK=1:D%NKT)
  ZVPT(:,:)=ZTHM(:,:) * ( 1. + ZRVORD*ZRM(:,:,1) )  &
                           / ( 1. + ZSUM(:,:) )
  !$mnh_end_expand_array(J1D=1:IPROMA,JK=1:D%NKT)
ELSE
  ZVPT(:,:)=ZTHM(:,:)
END IF
!
!!!!!!!!!!!!
!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!! WARNING !!
!!!!!!!!!!!!
!!!!!!!!!!!!
!Any modification done to the following lines and to the sections 4 and
!6 must be copied in compute_bl89_ml routine.
!We do not call directly this routine for numerical performance reasons
!but algorithm must remain the same.
!!!!!!!!!!!!
!
DO JK=D%NKTB,D%NKTE
  DO J1D=1,IPROMA
    ZDELTVPT(J1D,JK) = ZVPT(J1D,JK) - ZVPT(J1D,JK-D%NKL)
    ZHLVPT(J1D,JK) = 0.5 * ( ZVPT(J1D,JK) + ZVPT(J1D,JK-D%NKL) )
  END DO
END DO
!
DO J1D=1,IPROMA
  ZDELTVPT(J1D,D%NKU) = ZVPT(J1D,D%NKU) - ZVPT(J1D,D%NKU-D%NKL)
  ZDELTVPT(J1D,D%NKA) = 0.
  ZHLVPT(J1D,D%NKU) = 0.5 * ( ZVPT(J1D,D%NKU) + ZVPT(J1D,D%NKU-D%NKL) )
  ZHLVPT(J1D,D%NKA) = ZVPT(J1D,D%NKA)
END DO
!
DO JK=1,D%NKT
  DO J1D=1,IPROMA
    IF(ABS(ZDELTVPT(J1D,JK))<CSTURB%XLINF) THEN
      ZDELTVPT(J1D,JK)=CSTURB%XLINF
    END IF
  END DO
END DO
!
!-------------------------------------------------------------------------------
!
!*       3.  loop on model levels
!            --------------------
!
DO JK=D%NKTB,D%NKTE
!
!-------------------------------------------------------------------------------
!
!*       4.  mixing length for a downwards displacement
!            ------------------------------------------
  ZINTE(:)=ZTKEM(:,JK)
  ZLWORK=0.
  ZTESTM=1.
  DO JKK=JK,D%NKB,-D%NKL
    IF(ZTESTM > 0.) THEN
      ZTESTM=0.
      DO J1D=1,IPROMA
        ZTEST0=0.5+SIGN(0.5,ZINTE(J1D))
        !--------- SHEAR + STABILITY -----------
        ZPOTE = ZTEST0* &
                (-ZG_O_THVREF(J1D,JK)*(ZHLVPT(J1D,JKK)-ZVPT(J1D,JK)) &
                + CSTURB%XRM17*ZSHEAR(J1D,JKK)*ZSQRT_TKE(J1D,JK) &
                )*ZDZZ(J1D,JKK)

        ZTEST =0.5+SIGN(0.5,ZINTE(J1D)-ZPOTE)
        ZTESTM=ZTESTM+ZTEST0
        ZLWORK1=ZDZZ(J1D,JKK)

        !--------- SHEAR + STABILITY ----------- 
        ZLWORK2 = (ZG_O_THVREF(J1D,JK) *(ZVPT(J1D,JKK) - ZVPT(J1D,JK))  & 
          -CSTURB%XRM17*ZSHEAR(J1D,JKK)*ZSQRT_TKE(J1D,JK) &
          + sqrt(abs( (CSTURB%XRM17*ZSHEAR(J1D,JKK)*ZSQRT_TKE(J1D,JK) &
          + ( -ZG_O_THVREF(J1D,JK) * (ZVPT(J1D,JKK) - ZVPT(J1D,JK)) ))**2.0 + &
          2. * ZINTE(J1D) * &
#ifdef REPRO48
          ZG_O_THVREF(J1D,JK) * ZDELTVPT(J1D,JKK)/ ZDZZ(J1D,JKK)))) / &
#else
          (ZG_O_THVREF(J1D,JK) * ZDELTVPT(J1D,JKK)/ ZDZZ(J1D,JKK))))) / &
#endif
          (ZG_O_THVREF(J1D,JK) * ZDELTVPT(J1D,JKK) / ZDZZ(J1D,JKK))
        ZLWORK(J1D)=ZLWORK(J1D)+ZTEST0*(ZTEST*ZLWORK1+(1-ZTEST)*ZLWORK2)
        ZINTE(J1D) = ZINTE(J1D) - ZPOTE
      END DO
    ENDIF
  END DO
!-------------------------------------------------------------------------------
!
!*       5.  intermediate storage of the final mixing length
!            -----------------------------------------------
!
  DO J1D=1,IPROMA
    ZLMDN(J1D,JK)=MIN(ZLWORK(J1D),0.5*(ZZZ(J1D,JK)+ZZZ(J1D,JK+D%NKL))-ZZZ(J1D,D%NKB))
  END DO
!
!-------------------------------------------------------------------------------
!
!*       6.  mixing length for an upwards displacement
!            -----------------------------------------
!
  ZINTE(:)=ZTKEM(:,JK)
  ZLWORK(:)=0.
  ZTESTM=1.
!
  DO JKK=JK+D%NKL,D%NKE,D%NKL
    IF(ZTESTM > 0.) THEN
      ZTESTM=0.
      DO J1D=1,D%NIT*D%NJT
        ZTEST0=0.5+SIGN(0.5,ZINTE(J1D))
        !--------- SHEAR + STABILITY -----------
        ZPOTE = ZTEST0* &
                (ZG_O_THVREF(J1D,JK)*(ZHLVPT(J1D,JKK)-ZVPT(J1D,JK)) &
                +CSTURB%XRM17*ZSHEAR(J1D,JKK)*ZSQRT_TKE(J1D,JK) &
                )*ZDZZ(J1D,JKK)
        ZTEST =0.5+SIGN(0.5,ZINTE(J1D)-ZPOTE)
        ZTESTM=ZTESTM+ZTEST0
        ZLWORK1=ZDZZ(J1D,JKK)
        !--------- SHEAR + STABILITY ----------- 
        ZLWORK2= ( - ZG_O_THVREF(J1D,JK) *(ZVPT(J1D,JKK-D%NKL) - ZVPT(J1D,JK) )  &
                   - CSTURB%XRM17*ZSHEAR(J1D,JKK)*ZSQRT_TKE(J1D,JK)  &
          + SQRT (ABS(                                                       &
          (CSTURB%XRM17*ZSHEAR(J1D,JKK)*ZSQRT_TKE(J1D,JK)   &
            + ( ZG_O_THVREF(J1D,JK) * (ZVPT(J1D,JKK-D%NKL) - ZVPT(J1D,JK))) )**2    &
            + 2. * ZINTE(J1D) * &
#ifdef REPRO48
             ZG_O_THVREF(J1D,JK)* ZDELTVPT(J1D,JKK)/ZDZZ(J1D,JKK)))) / &
#else
            (ZG_O_THVREF(J1D,JK)* ZDELTVPT(J1D,JKK)/ZDZZ(J1D,JKK))))) / &
#endif
            (ZG_O_THVREF(J1D,JK) * ZDELTVPT(J1D,JKK) / ZDZZ(J1D,JKK))
        ZLWORK(J1D)=ZLWORK(J1D)+ZTEST0*(ZTEST*ZLWORK1+(1-ZTEST)*ZLWORK2)
        ZINTE(J1D) = ZINTE(J1D) - ZPOTE
      END DO
    ENDIF
  END DO
!
!-------------------------------------------------------------------------------
!
!*       7.  final mixing length
!
  DO J1D=1,IPROMA
    ZLWORK1=MAX(ZLMDN(J1D,JK),1.E-10_MNHREAL)
    ZLWORK2=MAX(ZLWORK(J1D),1.E-10_MNHREAL)
    ZPOTE = ZLWORK1 / ZLWORK2
#ifdef REPRO48
    ZLWORK2=1.d0 + ZPOTE**(2./3.)
    ZLM(J1D,JK) = Z2SQRT2*ZLWORK1/(ZLWORK2*SQRT(ZLWORK2))
#else
    ZLWORK2=1.d0 + ZPOTE**ZBL89EXP
    ZLM(J1D,JK) = ZLWORK1*(2./ZLWORK2)**ZUSRBL89
#endif
    ZLM(J1D,JK)=MAX(ZLM(J1D,JK),CSTURB%XLINI)
  END DO


!-------------------------------------------------------------------------------
!*       8.  end of the loop on the vertical levels
!            --------------------------------------
!
END DO
!
!-------------------------------------------------------------------------------
!
!*       9.  boundaries
!            ----------
!
ZLM(:,D%NKA)=ZLM(:,D%NKB)
ZLM(:,D%NKE)=ZLM(:,D%NKE-D%NKL)
ZLM(:,D%NKU)=ZLM(:,D%NKE-D%NKL)
!
!-------------------------------------------------------------------------------
!
!*      10.  retrieve output array in model coordinates
!            ------------------------------------------
! Not needed anymore because of the use of Pointer remapping (see 1.)
! PLM (3D array) is the target of ZLM (2D array) in a contiguous way
!
IF (LHOOK) CALL DR_HOOK('BL89',1,ZHOOK_HANDLE)
END SUBROUTINE BL89
END MODULE MODE_BL89
