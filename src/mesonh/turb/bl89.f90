!MNH_LIC Copyright 1997-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ################
      MODULE MODI_BL89
!     ################
INTERFACE
      SUBROUTINE BL89(KKA,KKU,KKL,PZZ,PDZZ,PTHVREF,PTHLM,KRR,PRM,PTKEM,PSHEAR,PLM)
!
INTEGER,                  INTENT(IN)  :: KKA 
INTEGER,                  INTENT(IN)  :: KKU 
INTEGER,                  INTENT(IN)  :: KKL 
REAL, DIMENSION(:,:,:),   INTENT(IN)  :: PZZ
REAL, DIMENSION(:,:,:),   INTENT(IN)  :: PDZZ
REAL, DIMENSION(:,:,:),   INTENT(IN)  :: PTHVREF
REAL, DIMENSION(:,:,:),   INTENT(IN)  :: PTHLM
INTEGER,                  INTENT(IN)  :: KRR
REAL, DIMENSION(:,:,:,:), INTENT(IN)  :: PRM
REAL, DIMENSION(:,:,:),   INTENT(IN)  :: PTKEM
REAL, DIMENSION(:,:,:),   INTENT(IN)  :: PSHEAR
REAL, DIMENSION(:,:,:),   INTENT(OUT) :: PLM

END SUBROUTINE BL89
END INTERFACE
END MODULE MODI_BL89
!
!     #########################################################
      SUBROUTINE BL89(KKA,KKU,KKL,PZZ,PDZZ,PTHVREF,PTHLM,KRR,PRM,PTKEM,PSHEAR,PLM)
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
USE MODD_CONF,      ONLY: CPROGRAM
USE MODD_CST
USE MODD_CTURB
USE MODD_DYN_n,     ONLY: LOCEAN
USE MODD_PARAMETERS
use modd_precision, only: MNHREAL
!
!
IMPLICIT NONE
!
!*       0.1   Declaration of arguments
!              ------------------------
!
INTEGER,                  INTENT(IN)  :: KKA      !near ground array index  
INTEGER,                  INTENT(IN)  :: KKU      !uppest atmosphere array index
INTEGER,                  INTENT(IN)  :: KKL      !vert. levels type 1=MNH -1=ARO
REAL, DIMENSION(:,:,:),   INTENT(IN)  :: PZZ
REAL, DIMENSION(:,:,:),   INTENT(IN)  :: PDZZ
REAL, DIMENSION(:,:,:),   INTENT(IN)  :: PTHVREF
REAL, DIMENSION(:,:,:),   INTENT(IN)  :: PTHLM       ! conservative pot. temp.
INTEGER,                  INTENT(IN)  :: KRR
REAL, DIMENSION(:,:,:,:), INTENT(IN)  :: PRM       ! water var.
REAL, DIMENSION(:,:,:),   INTENT(IN)  :: PTKEM     ! TKE
REAL, DIMENSION(:,:,:),   INTENT(IN)  :: PSHEAR
REAL, DIMENSION(:,:,:),   INTENT(OUT) :: PLM       ! Mixing length
!   thermodynamical variables PTHLM=Theta at the begining
!
!*       0.2   Declaration of local variables
!              ------------------------------
!
INTEGER :: IKB,IKE
INTEGER :: IKT          ! array size in k direction
INTEGER :: IKTB,IKTE    ! start, end of k loops in physical domain 

REAL, DIMENSION(SIZE(PTKEM,1)*SIZE(PTKEM,2),SIZE(PTKEM,3)) :: ZVPT  ! Virtual Potential Temp at half levels
REAL, DIMENSION(SIZE(PTKEM,1)*SIZE(PTKEM,2),SIZE(PTKEM,3)) :: ZDELTVPT
            ! Increment of Virtual Potential Temp between two following levels
REAL, DIMENSION(SIZE(PTKEM,1)*SIZE(PTKEM,2),SIZE(PTKEM,3)) :: ZHLVPT
            ! Virtual Potential Temp at half levels
REAL, DIMENSION(SIZE(PTKEM,1)*SIZE(PTKEM,2)) ::  ZLWORK,ZINTE
!           ! downwards then upwards vertical displacement,
!           ! residual internal energy,
!           ! residual potential energy
REAL, DIMENSION(SIZE(PTKEM,1)*SIZE(PTKEM,2),SIZE(PTKEM,3)) :: ZZZ,ZDZZ,       &
                                                              ZG_O_THVREF,    &
                                                              ZTHM,ZTKEM,ZLM, &
                                                              ZLMDN,ZSHEAR,   &
                                                              ZSQRT_TKE
!           ! input and output arrays packed according one horizontal coord.
REAL, DIMENSION(SIZE(PRM,1)*SIZE(PRM,2),SIZE(PRM,3),SIZE(PRM,4)) :: ZRM
!           ! input array packed according one horizontal coord.
REAL, DIMENSION(SIZE(PRM,1)*SIZE(PRM,2),SIZE(PRM,3)) :: ZSUM ! to replace SUM function
!
INTEGER :: IIU,IJU
INTEGER :: J1D        ! horizontal loop counter
INTEGER :: JK,JKK,J3RD     ! loop counters
INTEGER :: JRR        ! moist loop counter
REAL    :: ZRVORD     ! Rv/Rd
REAL    :: ZPOTE,ZLWORK1,ZLWORK2
REAL    :: ZTEST,ZTEST0,ZTESTM ! test for vectorization
REAL    :: Z2SQRT2,ZUSRBL89,ZBL89EXP
!-------------------------------------------------------------------------------
!
Z2SQRT2=2.*SQRT(2.)
IIU=SIZE(PTKEM,1)
IJU=SIZE(PTKEM,2)
!
IKB=KKA+JPVEXT_TURB*KKL
IKE=KKU-JPVEXT_TURB*KKL

IKTB = JPVEXT_TURB + 1
IKT = SIZE(PTKEM,3)
IKTE = IKT-JPVEXT_TURB
ZRVORD = XRV / XRD
!
!-------------------------------------------------------------------------------
!
!*       1.    pack the horizontal dimensions into one
!              ---------------------------------------
!
IF (CPROGRAM=='AROME ') THEN
  DO JK=1,IKT
    ZZZ    (:,JK)   = PZZ    (:,1,JK)
    ZDZZ   (:,JK)   = PDZZ   (:,1,JK)
    ZTHM   (:,JK)   = PTHLM  (:,1,JK)
    ZTKEM  (:,JK)   = PTKEM  (:,1,JK)
    ZG_O_THVREF(:,JK)   = XG/PTHVREF(:,1,JK)
  END DO
  DO JK=1,IKT
    DO JRR=1,KRR
      ZRM  (:,JK,JRR) = PRM    (:,1,JK,JRR)
    END DO
  END DO
ELSE
  DO JK=1,IKT
    ZZZ    (:,JK)   = RESHAPE(PZZ    (:,:,JK),(/ IIU*IJU /) )
    ZDZZ   (:,JK)   = RESHAPE(PDZZ   (:,:,JK),(/ IIU*IJU /) )
    ZTHM   (:,JK)   = RESHAPE(PTHLM  (:,:,JK),(/ IIU*IJU /) )
    ZSHEAR   (:,JK)   = RESHAPE(PSHEAR  (:,:,JK),(/ IIU*IJU /) )    
    ZTKEM  (:,JK)   = RESHAPE(PTKEM  (:,:,JK),(/ IIU*IJU /) )
    ZG_O_THVREF(:,JK)   = RESHAPE(XG/PTHVREF(:,:,JK),(/ IIU*IJU /) )
    IF (LOCEAN) ZG_O_THVREF(:,JK)   = XG * XALPHAOC
    DO JRR=1,KRR
      ZRM  (:,JK,JRR) = RESHAPE(PRM    (:,:,JK,JRR),(/ IIU*IJU /) )
    END DO
  END DO
END IF
!
ZSQRT_TKE = SQRT(ZTKEM)
!
!ZBL89EXP is defined here because (and not in ini_cturb) because XCED is defined in read_exseg (depending on BL89/RM17)
ZBL89EXP = LOG(16.)/(4.*LOG(XKARMAN)+LOG(XCED)-3.*LOG(XCMFS))
ZUSRBL89 = 1./ZBL89EXP
!-------------------------------------------------------------------------------
!
!*       2.    Virtual potential temperature on the model grid
!              -----------------------------------------------
!
IF(KRR /= 0) THEN
  ZSUM(:,:) = 0.
  DO JRR=1,KRR
    ZSUM(:,:) = ZSUM(:,:)+ZRM(:,:,JRR)
  ENDDO
  ZVPT(:,1:)=ZTHM(:,:) * ( 1. + ZRVORD*ZRM(:,:,1) )  &
                           / ( 1. + ZSUM(:,:) )
ELSE
  ZVPT(:,1:)=ZTHM(:,:)
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

ZDELTVPT(:,IKTB:IKTE)=ZVPT(:,IKTB:IKTE)-ZVPT(:,IKTB-KKL:IKTE-KKL)
ZDELTVPT(:,KKU)=ZVPT(:,KKU)-ZVPT(:,KKU-KKL)
ZDELTVPT(:,KKA)=0.
WHERE (ABS(ZDELTVPT(:,:))<XLINF)
  ZDELTVPT(:,:)=XLINF
END WHERE
!
ZHLVPT(:,IKTB:IKTE)= 0.5 * ( ZVPT(:,IKTB:IKTE)+ZVPT(:,IKTB-KKL:IKTE-KKL) )
ZHLVPT(:,KKU)= 0.5 * ( ZVPT(:,KKU)+ZVPT(:,KKU-KKL) )
ZHLVPT(:,KKA)    =         ZVPT(:,KKA)
!-------------------------------------------------------------------------------
!
!*       3.  loop on model levels
!            --------------------
DO JK=IKTB,IKTE
!
!-------------------------------------------------------------------------------
!
!

!*       4.  mixing length for a downwards displacement
!            ------------------------------------------
  ZINTE(:)=ZTKEM(:,JK)
  ZLWORK=0.
  ZTESTM=1.
  DO JKK=JK,IKB,-KKL  
    IF(ZTESTM > 0.) THEN
      ZTESTM=0.
      DO J1D=1,IIU*IJU
        
       ZTEST0=0.5+SIGN(0.5,ZINTE(J1D))
        
        !--------- SHEAR + STABILITY -----------
        ZPOTE = ZTEST0* &
                (-ZG_O_THVREF(J1D,JK)*(ZHLVPT(J1D,JKK)-ZVPT(J1D,JK)) &
                + XRM17*ZSHEAR(J1D,JKK)*ZSQRT_TKE(J1D,JK) &
                )*ZDZZ(J1D,JKK)

        ZTEST =0.5+SIGN(0.5,ZINTE(J1D)-ZPOTE)
        ZTESTM=ZTESTM+ZTEST0
        ZLWORK1=ZDZZ(J1D,JKK)
        
        !-------- ORIGINAL -------------        
!        ZLWORK2= ( + ZG_O_THVREF(J1D,JK) *                                   &
!            (  ZVPT(J1D,JKK) - ZVPT(J1D,JK) )                                &
!          + SQRT (ABS(                                                       &
!            ( ZG_O_THVREF(J1D,JK) * (ZVPT(J1D,JKK) - ZVPT(J1D,JK)) )**2      &
!            + 2. * ZINTE(J1D) * ZG_O_THVREF(J1D,JK)                          &
!                 * ZDELTVPT(J1D,JKK) / ZDZZ(J1D,JKK)   ))) /                 &
!        ( ZG_O_THVREF(J1D,JK) * ZDELTVPT(J1D,JKK) / ZDZZ(J1D,JKK)) 
        

        !--------- SHEAR + STABILITY ----------- 
    ZLWORK2 = (ZG_O_THVREF(J1D,JK) *(ZVPT(J1D,JKK) - ZVPT(J1D,JK))  & 
          -XRM17*ZSHEAR(J1D,JKK)*ZSQRT_TKE(J1D,JK) &
          + sqrt(abs( (XRM17*ZSHEAR(J1D,JKK)*ZSQRT_TKE(J1D,JK) &
          + ( -ZG_O_THVREF(J1D,JK) * (ZVPT(J1D,JKK) - ZVPT(J1D,JK)) ))**2.0 + &
          2. * ZINTE(J1D) * &
          (ZG_O_THVREF(J1D,JK) * ZDELTVPT(J1D,JKK)/ ZDZZ(J1D,JKK))))) / &
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
  ZLMDN(:,JK)=MIN(ZLWORK(:),0.5*(ZZZ(:,JK)+ZZZ(:,JK+KKL))-ZZZ(:,IKB))
!
!-------------------------------------------------------------------------------
!
!*       6.  mixing length for an upwards displacement
!            -----------------------------------------
!
  ZINTE(:)=ZTKEM(:,JK)
  ZLWORK=0.
  ZTESTM=1.
!
  DO JKK=JK+KKL,IKE,KKL 
    IF(ZTESTM > 0.) THEN
      ZTESTM=0.
      DO J1D=1,IIU*IJU
        ZTEST0=0.5+SIGN(0.5,ZINTE(J1D))

        !-------- ORIGINAL -------------        
       !ZPOTE = ZTEST0*ZG_O_THVREF(J1D,JK)      *      &
       !        (ZHLVPT(J1D,JKK) - ZVPT(J1D,JK) ) *ZDZZ(J1D,JKK)

        !--------- SHEAR + STABILITY -----------
        ZPOTE = ZTEST0* &
                (ZG_O_THVREF(J1D,JK)*(ZHLVPT(J1D,JKK)-ZVPT(J1D,JK)) &
                +XRM17*ZSHEAR(J1D,JKK)*ZSQRT_TKE(J1D,JK) &
                )*ZDZZ(J1D,JKK)

        ZTEST =0.5+SIGN(0.5,ZINTE(J1D)-ZPOTE)
        ZTESTM=ZTESTM+ZTEST0
        ZLWORK1=ZDZZ(J1D,JKK)

        !-------- ORIGINAL -------------        
       ! ZLWORK2= ( - ZG_O_THVREF(J1D,JK) *                                   &
       !     (  ZVPT(J1D,JKK-KKL) - ZVPT(J1D,JK) )                              &
       !   + SQRT (ABS(                                                       &
       !     ( ZG_O_THVREF(J1D,JK) * (ZVPT(J1D,JKK-KKL) - ZVPT(J1D,JK)) )**2    &
       !     + 2. * ZINTE(J1D) * ZG_O_THVREF(J1D,JK)                          &
       !          * ZDELTVPT(J1D,JKK) / ZDZZ(J1D,JKK) ))    ) /               &
       ! ( ZG_O_THVREF(J1D,JK) * ZDELTVPT(J1D,JKK) / ZDZZ(J1D,JKK) ) 
        
        !--------- SHEAR + STABILITY ----------- 
       ZLWORK2= ( - ZG_O_THVREF(J1D,JK) *(ZVPT(J1D,JKK-KKL) - ZVPT(J1D,JK) )  &
                   - XRM17*ZSHEAR(J1D,JKK)*ZSQRT_TKE(J1D,JK)  &
          + SQRT (ABS(                                                       &
          (XRM17*ZSHEAR(J1D,JKK)*ZSQRT_TKE(J1D,JK)   &
            + ( ZG_O_THVREF(J1D,JK) * (ZVPT(J1D,JKK-KKL) - ZVPT(J1D,JK))) )**2    &
            + 2. * ZINTE(J1D) * &
            ( ZG_O_THVREF(J1D,JK)* ZDELTVPT(J1D,JKK)/ZDZZ(J1D,JKK))))) / &
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
  DO J1D=1,IIU*IJU
    ZLWORK1=MAX(ZLMDN(J1D,JK),1.E-10_MNHREAL)
    ZLWORK2=MAX(ZLWORK(J1D),1.E-10_MNHREAL)
    ZPOTE = ZLWORK1 / ZLWORK2
    ZLWORK2=1.d0 + ZPOTE**ZBL89EXP
    ZLM(J1D,JK) = ZLWORK1*(2./ZLWORK2)**ZUSRBL89
  END DO 

ZLM(:,JK)=MAX(ZLM(:,JK),XLINI)

!
!
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
ZLM(:,KKA)=ZLM(:,IKB)
ZLM(:,IKE)=ZLM(:,IKE-KKL)
ZLM(:,KKU)=ZLM(:,IKE-KKL)
!
!-------------------------------------------------------------------------------
!
!*      10.  retrieve output array in model coordinates
!            ------------------------------------------
!
IF (CPROGRAM=='AROME ') THEN
  DO JK=1,IKT
    PLM  (:,1,JK)   = ZLM  (:,JK)
  END DO
ELSE
  DO JK=1,IKT
    PLM  (:,:,JK)   = RESHAPE(ZLM  (:,JK), (/ IIU,IJU /) )
  END DO
END IF

!
END SUBROUTINE BL89
