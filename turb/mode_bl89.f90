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
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN),TARGET  :: PZZ
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN),TARGET  :: PDZZ
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN),TARGET  :: PTHVREF
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN),TARGET  :: PTHLM       ! conservative pot. temp.
INTEGER,                  INTENT(IN)  :: KRR
REAL, DIMENSION(D%NIJT,D%NKT,KRR), INTENT(IN),TARGET :: PRM       ! water var.
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN),TARGET  :: PTKEM     ! TKE
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN),TARGET  :: PSHEAR
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(OUT),TARGET :: PLM       ! Mixing length
LOGICAL,                  INTENT(IN)  ::  OOCEAN       ! switch for Ocean model version
CHARACTER(LEN=6),         INTENT(IN)  ::  HPROGRAM     ! CPROGRAM is the program currently running (modd_conf)
!   thermodynamical variables PTHLM=Theta at the begining
!
!*       0.2   Declaration of local variables
!              ------------------------------
!
INTEGER :: IKT,IKB,IKA,IKU

REAL, DIMENSION(D%NIJT,D%NKT) :: ZVPT  ! Virtual Potential Temp at half levels
REAL, DIMENSION(D%NIJT,D%NKT) :: ZDELTVPT
            ! Increment of Virtual Potential Temp between two following levels
REAL, DIMENSION(D%NIJT,D%NKT) :: ZHLVPT
            ! Virtual Potential Temp at half levels
REAL, DIMENSION(D%NIJT) ::  ZLWORK,ZINTE
!           ! downwards then upwards vertical displacement,
!           ! residual internal energy,
!           ! residual potential energy
!           ! input and output arrays packed according one horizontal coord.
!           ! input array packed according one horizontal coord.
REAL, DIMENSION(D%NIJT,D%NKT) :: ZSUM ! to replace SUM function
REAL, DIMENSION(D%NIJT,D%NKT) :: ZG_O_THVREF
REAL, DIMENSION(D%NIJT,D%NKT) :: ZSQRT_TKE
REAL, DIMENSION(D%NIJT,D%NKT) :: PLMDN
!
INTEGER :: IIJB, IIJE
INTEGER :: IKTB, IKTE, IKE,IKL
INTEGER :: JIJ        ! horizontal loop counter
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
!
IIJB=D%NIJB
IIJE=D%NIJE
IKTB=D%NKTB
IKTE=D%NKTE
IKT=D%NKT
IKB=D%NKB
IKE=D%NKE
IKA=D%NKA
IKU=D%NKU
IKL=D%NKL
!-------------------------------------------------------------------------------
!
!*       1.    pack the horizontal dimensions into one
!              ---------------------------------------
!
! Pointer remapping instead of RESHAPE (contiguous memory)
! 2D array => 3D array
!
IF (OOCEAN) THEN
  DO JK=1,IKT
    DO JIJ=IIJB,IIJE
      ZG_O_THVREF(JIJ,JK) = CST%XG * CST%XALPHAOC
    END DO
  END DO
ELSE !Atmosphere case
  DO JK=1,IKT
    DO JIJ=IIJB,IIJE
      ZG_O_THVREF(JIJ,JK) = CST%XG / PTHVREF(JIJ,JK)
    END DO
  END DO
END IF
!
!$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
ZSQRT_TKE(:,:) = SQRT(PTKEM(:,:))
!$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
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
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ZSUM(:,:) = ZSUM(:,:)+PRM(:,:,JRR)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  ENDDO
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  ZVPT(:,:)=PTHLM(:,:) * ( 1. + ZRVORD*PRM(:,:,1) )  &
                           / ( 1. + ZSUM(:,:) )
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
ELSE
  ZVPT(:,:)=PTHLM(:,:)
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
DO JK=IKTB,IKTE
  DO JIJ=IIJB,IIJE
    ZDELTVPT(JIJ,JK) = ZVPT(JIJ,JK) - ZVPT(JIJ,JK-IKL)
    ZHLVPT(JIJ,JK) = 0.5 * ( ZVPT(JIJ,JK) + ZVPT(JIJ,JK-IKL) )
  END DO
END DO
!
DO JIJ=IIJB,IIJE
  ZDELTVPT(JIJ,IKU) = ZVPT(JIJ,IKU) - ZVPT(JIJ,IKU-IKL)
  ZDELTVPT(JIJ,IKA) = 0.
  ZHLVPT(JIJ,IKU) = 0.5 * ( ZVPT(JIJ,IKU) + ZVPT(JIJ,IKU-IKL) )
  ZHLVPT(JIJ,IKA) = ZVPT(JIJ,IKA)
END DO
!
DO JK=1,IKT
  DO JIJ=IIJB,IIJE
    IF(ABS(ZDELTVPT(JIJ,JK))<CSTURB%XLINF) THEN
      ZDELTVPT(JIJ,JK)=CSTURB%XLINF
    END IF
  END DO
END DO
!
!-------------------------------------------------------------------------------
!
!*       3.  loop on model levels
!            --------------------
!
DO JK=IKTB,IKTE
!
!-------------------------------------------------------------------------------
!
!*       4.  mixing length for a downwards displacement
!            ------------------------------------------
  ZINTE(:)=PTKEM(:,JK)
  ZLWORK=0.
  ZTESTM=1.
  DO JKK=JK,IKB,-IKL
    IF(ZTESTM > 0.) THEN
      ZTESTM=0.
      DO JIJ=IIJB,IIJE
        ZTEST0=0.5+SIGN(0.5,ZINTE(JIJ))
        !--------- SHEAR + STABILITY -----------
        ZPOTE = ZTEST0* &
                (-ZG_O_THVREF(JIJ,JK)*(ZHLVPT(JIJ,JKK)-ZVPT(JIJ,JK)) &
                + CSTURB%XRM17*PSHEAR(JIJ,JKK)*ZSQRT_TKE(JIJ,JK) &
                )*PDZZ(JIJ,JKK)

        ZTEST =0.5+SIGN(0.5,ZINTE(JIJ)-ZPOTE)
        ZTESTM=ZTESTM+ZTEST0
        ZLWORK1=PDZZ(JIJ,JKK)

        !--------- SHEAR + STABILITY ----------- 
        ZLWORK2 = (ZG_O_THVREF(JIJ,JK) *(ZVPT(JIJ,JKK) - ZVPT(JIJ,JK))  & 
          -CSTURB%XRM17*PSHEAR(JIJ,JKK)*ZSQRT_TKE(JIJ,JK) &
          + sqrt(abs( (CSTURB%XRM17*PSHEAR(JIJ,JKK)*ZSQRT_TKE(JIJ,JK) &
          + ( -ZG_O_THVREF(JIJ,JK) * (ZVPT(JIJ,JKK) - ZVPT(JIJ,JK)) ))**2.0 + &
          2. * ZINTE(JIJ) * &
#ifdef REPRO48
          ZG_O_THVREF(JIJ,JK) * ZDELTVPT(JIJ,JKK)/ PDZZ(JIJ,JKK)))) / &
#else
          (ZG_O_THVREF(JIJ,JK) * ZDELTVPT(JIJ,JKK)/ PDZZ(JIJ,JKK))))) / &
#endif
          (ZG_O_THVREF(JIJ,JK) * ZDELTVPT(JIJ,JKK) / PDZZ(JIJ,JKK))
        ZLWORK(JIJ)=ZLWORK(JIJ)+ZTEST0*(ZTEST*ZLWORK1+(1-ZTEST)*ZLWORK2)
        ZINTE(JIJ) = ZINTE(JIJ) - ZPOTE
      END DO
    ENDIF
  END DO
!-------------------------------------------------------------------------------
!
!*       5.  intermediate storage of the final mixing length
!            -----------------------------------------------
!
  DO JIJ=IIJB,IIJE
    PLMDN(JIJ,JK)=MIN(ZLWORK(JIJ),0.5*(PZZ(JIJ,JK)+PZZ(JIJ,JK+IKL))-PZZ(JIJ,IKB))
  END DO
!
!-------------------------------------------------------------------------------
!
!*       6.  mixing length for an upwards displacement
!            -----------------------------------------
!
  ZINTE(:)=PTKEM(:,JK)
  ZLWORK(:)=0.
  ZTESTM=1.
!
  DO JKK=JK+IKL,IKE,IKL
    IF(ZTESTM > 0.) THEN
      ZTESTM=0.
      DO JIJ=IIJB,IIJE
        ZTEST0=0.5+SIGN(0.5,ZINTE(JIJ))
        !--------- SHEAR + STABILITY -----------
        ZPOTE = ZTEST0* &
                (ZG_O_THVREF(JIJ,JK)*(ZHLVPT(JIJ,JKK)-ZVPT(JIJ,JK)) &
                +CSTURB%XRM17*PSHEAR(JIJ,JKK)*ZSQRT_TKE(JIJ,JK) &
                )*PDZZ(JIJ,JKK)
        ZTEST =0.5+SIGN(0.5,ZINTE(JIJ)-ZPOTE)
        ZTESTM=ZTESTM+ZTEST0
        ZLWORK1=PDZZ(JIJ,JKK)
        !--------- SHEAR + STABILITY ----------- 
        ZLWORK2= ( - ZG_O_THVREF(JIJ,JK) *(ZVPT(JIJ,JKK-IKL) - ZVPT(JIJ,JK) )  &
                   - CSTURB%XRM17*PSHEAR(JIJ,JKK)*ZSQRT_TKE(JIJ,JK)  &
          + SQRT (ABS(                                                       &
          (CSTURB%XRM17*PSHEAR(JIJ,JKK)*ZSQRT_TKE(JIJ,JK)   &
            + ( ZG_O_THVREF(JIJ,JK) * (ZVPT(JIJ,JKK-IKL) - ZVPT(JIJ,JK))) )**2    &
            + 2. * ZINTE(JIJ) * &
#ifdef REPRO48
             ZG_O_THVREF(JIJ,JK)* ZDELTVPT(JIJ,JKK)/PDZZ(JIJ,JKK)))) / &
#else
            (ZG_O_THVREF(JIJ,JK)* ZDELTVPT(JIJ,JKK)/PDZZ(JIJ,JKK))))) / &
#endif
            (ZG_O_THVREF(JIJ,JK) * ZDELTVPT(JIJ,JKK) / PDZZ(JIJ,JKK))
        ZLWORK(JIJ)=ZLWORK(JIJ)+ZTEST0*(ZTEST*ZLWORK1+(1-ZTEST)*ZLWORK2)
        ZINTE(JIJ) = ZINTE(JIJ) - ZPOTE
      END DO
    ENDIF
  END DO
!
!-------------------------------------------------------------------------------
!
!*       7.  final mixing length
!
  DO JIJ=IIJB,IIJE
    ZLWORK1=MAX(PLMDN(JIJ,JK),1.E-10_MNHREAL)
    ZLWORK2=MAX(ZLWORK(JIJ),1.E-10_MNHREAL)
    ZPOTE = ZLWORK1 / ZLWORK2
#ifdef REPRO48
    ZLWORK2=1.d0 + ZPOTE**(2./3.)
    PLM(JIJ,JK) = Z2SQRT2*ZLWORK1/(ZLWORK2*SQRT(ZLWORK2))
#else
    ZLWORK2=1.d0 + ZPOTE**ZBL89EXP
    PLM(JIJ,JK) = ZLWORK1*(2./ZLWORK2)**ZUSRBL89
#endif
    PLM(JIJ,JK)=MAX(PLM(JIJ,JK),CSTURB%XLINI)
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
PLM(:,IKA)=PLM(:,IKB)
PLM(:,IKE)=PLM(:,IKE-IKL)
PLM(:,IKU)=PLM(:,IKE-IKL)
!
!-------------------------------------------------------------------------------
!
!*      10.  retrieve output array in model coordinates
!            ------------------------------------------
! Not needed anymore because of the use of Pointer remapping (see 1.)
! PLM (3D array) is the target of PLM (2D array) in a contiguous way
!
IF (LHOOK) CALL DR_HOOK('BL89',1,ZHOOK_HANDLE)
END SUBROUTINE BL89
END MODULE MODE_BL89
