!MNH_LIC Copyright 2006-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!    ###########################
     MODULE MODI_COMPUTE_BL89_ML
!    ###########################

INTERFACE

!     ###################################################################
      SUBROUTINE COMPUTE_BL89_ML(KKA,KKB,KKE,KKU,KKL,PDZZ2D, &
             PTKEM_DEP,PG_O_THVREF,PVPT,KK,OUPORDN,OFLUX,PSHEAR,PLWORK)
!     ###################################################################

!*               1.1  Declaration of Arguments

INTEGER,                INTENT(IN)   :: KKA          ! near ground array index
INTEGER,                INTENT(IN)   :: KKB          ! near ground physical index
INTEGER,                INTENT(IN)   :: KKE          ! uppest atmosphere physical index
INTEGER,                INTENT(IN)   :: KKU          ! uppest atmosphere array index
INTEGER,                INTENT(IN)   :: KKL          ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(:,:),   INTENT(IN)  :: PDZZ2D        ! height difference between two mass levels
REAL, DIMENSION(:),     INTENT(IN)  :: PTKEM_DEP     ! TKE to consume
REAL, DIMENSION(:),     INTENT(IN)  :: PG_O_THVREF   ! g/ThetaVRef at the departure point
REAL, DIMENSION(:,:),   INTENT(IN)  :: PVPT          ! ThetaV on mass levels
INTEGER,                INTENT(IN)  :: KK            ! index of departure level
LOGICAL,                INTENT(IN)  :: OUPORDN       ! switch to compute upward (true) or
                                                     !   downward (false) mixing length
LOGICAL,                INTENT(IN)  :: OFLUX         ! Computation must be done from flux level
REAL, DIMENSION(:),     INTENT(OUT) :: PLWORK        ! Resulting mixing length
REAL, DIMENSION(:,:),   INTENT(IN)  :: PSHEAR        ! vertical wind shear for RM17 mixing length

END SUBROUTINE COMPUTE_BL89_ML

END INTERFACE
!
END MODULE MODI_COMPUTE_BL89_ML
!     ######spl
      SUBROUTINE COMPUTE_BL89_ML(KKA,KKB,KKE,KKU,KKL,PDZZ2D, &
             PTKEM_DEP,PG_O_THVREF,PVPT,KK,OUPORDN,OFLUX,PSHEAR,PLWORK)

!     ###################################################################
!!
!!     COMPUTE_BL89_ML routine to:
!!       1/ compute upward or downward mixing length with BL89 formulation
!!
!!    AUTHOR
!!    ------
!!     J. PERGAUD
!!
!!    MODIFICATIONS
!!    -------------
!!     Original   19/01/06
!!     S. Riette Jan 2012: support for both order of vertical levels and cleaning
!!     R.Honnert Oct 2016 : Update with AROME
!!     Q.Rodier  01/2019 : support RM17 mixing length as in bl89.f90 
!  P. Wautelet 10/04/2019: replace ABORT and STOP calls by Print_msg
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
!              ------------
!
!!!!!!!!!!!!
!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!! WARNING !!
!!!!!!!!!!!!
!!!!!!!!!!!!
!Any modification done to this routine must be copied in bl89.f90.
!This routine was inlined in bl89 for numerical performance reasons
!but algorithm must remain the same.
!!!!!!!!!!!!
!
USE MODD_CTURB
USE MODD_PARAMETERS, ONLY: JPVEXT
!
use mode_msg
!
USE MODI_SHUMAN_MF
!
IMPLICIT NONE
!
!          0.1 arguments
!
INTEGER,                INTENT(IN)   :: KKA          ! near ground array index
INTEGER,                INTENT(IN)   :: KKB          ! near ground physical index
INTEGER,                INTENT(IN)   :: KKE          ! uppest atmosphere physical index
INTEGER,                INTENT(IN)   :: KKU          ! uppest atmosphere array index
INTEGER,                INTENT(IN)   :: KKL          ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(:,:),   INTENT(IN)  :: PDZZ2D        ! height difference between two mass levels
REAL, DIMENSION(:),     INTENT(IN)  :: PTKEM_DEP     ! TKE to consume
REAL, DIMENSION(:),     INTENT(IN)  :: PG_O_THVREF   ! g/ThetaVRef at the departure point
REAL, DIMENSION(:,:),   INTENT(IN)  :: PVPT          ! ThetaV on mass levels
INTEGER,                INTENT(IN)  :: KK            ! index of departure level
LOGICAL,                INTENT(IN)  :: OUPORDN       ! switch to compute upward (true) or
                                                     !   downward (false) mixing length
LOGICAL,                INTENT(IN)  :: OFLUX         ! Computation must be done from flux level
REAL, DIMENSION(:),     INTENT(OUT) :: PLWORK        ! Resulting mixing length
REAL, DIMENSION(:,:),   INTENT(IN)  :: PSHEAR        ! vertical wind shear for RM17 mixing length

!          0.2 Local variable
!
REAL, DIMENSION(SIZE(PVPT,1)) :: ZLWORK1,ZLWORK2 ! Temporary mixing length
REAL, DIMENSION(SIZE(PVPT,1)) :: ZINTE,ZPOTE     ! TKE and potential energy
                                                 !   between 2 levels
REAL, DIMENSION(SIZE(PVPT,1)) :: ZVPT_DEP        ! Thetav on departure point
!
REAL, DIMENSION(SIZE(PVPT,1),SIZE(PVPT,2)) :: ZDELTVPT,ZHLVPT                                
                      !Virtual Potential Temp at Half level and DeltaThv between
                      !2 mass levels

INTEGER :: IIJU                 !Internal Domain
INTEGER :: J1D                  !horizontal loop counter
INTEGER :: JKK                  !loop counters
REAL    :: ZTEST,ZTEST0,ZTESTM  !test for vectorization
!-------------------------------------------------------------------------------------
!
!*       1.    INITIALISATION
!              --------------
IIJU=SIZE(PVPT,1)
!
ZDELTVPT(:,:)=DZM_MF(KKA,KKU,KKL,PVPT(:,:))
ZDELTVPT(:,KKA)=0.
WHERE (ABS(ZDELTVPT(:,:))<XLINF)
  ZDELTVPT(:,:)=XLINF
END WHERE
!
ZHLVPT(:,:)=MZM_MF(KKA,KKU,KKL,PVPT(:,:))
!
!We consider that gradient between mass levels KKB and KKB+KKL is the same as
!the gradient between flux level KKB and mass level KKB
ZDELTVPT(:,KKB)=PDZZ2D(:,KKB)*ZDELTVPT(:,KKB+KKL)/PDZZ2D(:,KKB+KKL)
ZHLVPT(:,KKB)=PVPT(:,KKB)-ZDELTVPT(:,KKB)*0.5
!
!
!
!*       2.    CALCULATION OF THE UPWARD MIXING LENGTH
!              ---------------------------------------
!

IF (OUPORDN.EQV..TRUE.) THEN 
 ZINTE(:)=PTKEM_DEP(:)
 PLWORK=0.
 ZTESTM=1.
 IF(OFLUX)THEN
   ZVPT_DEP(:)=ZHLVPT(:,KK) ! departure point is on flux level
   !We must compute what happens between flux level KK and mass level KK
   DO J1D=1,IIJU
     ZTEST0=0.5+SIGN(0.5,ZINTE(J1D)) ! test if there's energy to consume
     ! Energy consumed if parcel cross the entire layer
     ZPOTE(J1D) = ZTEST0*(PG_O_THVREF(J1D)      *      &
         (0.5*(ZHLVPT(J1D,KK)+ PVPT(J1D,KK)) - ZVPT_DEP(J1D)) + &
         XRM17*PSHEAR(J1D,KK)*SQRT(ABS(PTKEM_DEP(J1D))))  * &
         PDZZ2D(J1D,KK)*0.5
     ! Test if it rests some energy to consume
     ZTEST =0.5+SIGN(0.5,ZINTE(J1D)-ZPOTE(J1D))
     ! Length travelled by parcel if it rests energy to consume
     ZLWORK1(J1D)=PDZZ2D(J1D,KK)*0.5
     ! Lenght travelled by parcel to nullify energy
      ZLWORK2(J1D)=        ( - PG_O_THVREF(J1D) *                     &
            (  ZHLVPT(J1D,KK) - ZVPT_DEP(J1D) )                              &
            - XRM17*PSHEAR(J1D,JKK)*sqrt(abs(PTKEM_DEP(J1D))) &
          + SQRT (ABS(                                                       &
            (XRM17*PSHEAR(J1D,JKK)*sqrt(abs(PTKEM_DEP(J1D))) +  &
             PG_O_THVREF(J1D) * (ZHLVPT(J1D,KK) - ZVPT_DEP(J1D)) )**2  &
            + 2. * ZINTE(J1D) * PG_O_THVREF(J1D)                        &
                 * ZDELTVPT(J1D,KK) / PDZZ2D(J1D,KK) ))    ) /             &
        ( PG_O_THVREF(J1D) * ZDELTVPT(J1D,KK) / PDZZ2D(J1D,KK) )  
     ! Effective length travelled by parcel
      PLWORK(J1D)=PLWORK(J1D)+ZTEST0*(ZTEST*ZLWORK1(J1D)+  &
                                    (1-ZTEST)*ZLWORK2(J1D))
      ! Rest of energy to consume
      ZINTE(J1D) = ZINTE(J1D) - ZPOTE(J1D)
   ENDDO
 ELSE
   ZVPT_DEP(:)=PVPT(:,KK) ! departure point is on mass level
 ENDIF

 DO JKK=KK+KKL,KKE,KKL
    IF(ZTESTM > 0.) THEN
      ZTESTM=0
      DO J1D=1,IIJU
        ZTEST0=0.5+SIGN(0.5,ZINTE(J1D))
        ZPOTE(J1D) = ZTEST0*(PG_O_THVREF(J1D)      *      &
            (ZHLVPT(J1D,JKK) - ZVPT_DEP(J1D))   &
           + XRM17*PSHEAR(J1D,JKK)*SQRT(ABS(PTKEM_DEP(J1D))))* PDZZ2D(J1D,JKK) 
        ZTEST =0.5+SIGN(0.5,ZINTE(J1D)-ZPOTE(J1D))
        ZTESTM=ZTESTM+ZTEST0
        ZLWORK1(J1D)=PDZZ2D(J1D,JKK)
        !ZLWORK2 jump of the last reached level
        ZLWORK2(J1D)=        ( - PG_O_THVREF(J1D) *                     &
            (  PVPT(J1D,JKK-KKL) - ZVPT_DEP(J1D) )                              &
            - XRM17*PSHEAR(J1D,JKK)*sqrt(abs(PTKEM_DEP(J1D))) &
          + SQRT (ABS(                                                   &
            (XRM17*PSHEAR(J1D,JKK)*sqrt(abs(PTKEM_DEP(J1D))) +  &
             PG_O_THVREF(J1D) * (PVPT(J1D,JKK-KKL) - ZVPT_DEP(J1D)) )**2  &
            + 2. * ZINTE(J1D) * PG_O_THVREF(J1D)                        &
                 * ZDELTVPT(J1D,JKK) / PDZZ2D(J1D,JKK) ))    ) /             &
        ( PG_O_THVREF(J1D) * ZDELTVPT(J1D,JKK) / PDZZ2D(J1D,JKK) ) 
!
        PLWORK(J1D)=PLWORK(J1D)+ZTEST0*(ZTEST*ZLWORK1(J1D)+  &
                                    (1-ZTEST)*ZLWORK2(J1D))
        ZINTE(J1D) = ZINTE(J1D) - ZPOTE(J1D)
      END DO 
    ENDIF
  END DO 
ENDIF
!!
!*       2.    CALCULATION OF THE DOWNWARD MIXING LENGTH
!              ---------------------------------------
!

IF (OUPORDN.EQV..FALSE.) THEN 
 IF(OFLUX) call Print_msg(NVERB_FATAL,'GEN','COMPUTE_BL89_ML','OFLUX option not coded for downward mixing length')
 ZINTE(:)=PTKEM_DEP(:)
 PLWORK=0.
 ZTESTM=1.
 DO JKK=KK,KKB,-KKL
    IF(ZTESTM > 0.) THEN
      ZTESTM=0
      DO J1D=1,IIJU
        ZTEST0=0.5+SIGN(0.5,ZINTE(J1D))
         ZPOTE(J1D) = ZTEST0*(-PG_O_THVREF(J1D)      *      &
            (ZHLVPT(J1D,JKK) - PVPT(J1D,KK)) &
         + XRM17*PSHEAR(J1D,JKK)*SQRT(ABS(PTKEM_DEP(J1D))))* PDZZ2D(J1D,JKK) 
        ZTEST =0.5+SIGN(0.5,ZINTE(J1D)-ZPOTE(J1D))
        ZTESTM=ZTESTM+ZTEST0
        ZLWORK1(J1D)=PDZZ2D(J1D,JKK)
      ZLWORK2(J1D)=        ( + PG_O_THVREF(J1D) *                     &
            (  PVPT(J1D,JKK) - PVPT(J1D,KK) )                              &
             -XRM17*PSHEAR(J1D,JKK)*sqrt(abs(PTKEM_DEP(J1D))) &
          + SQRT (ABS(                                                       &
            (XRM17*PSHEAR(J1D,JKK)*sqrt(abs(PTKEM_DEP(J1D))) - &
             PG_O_THVREF(J1D) * (PVPT(J1D,JKK) - PVPT(J1D,KK)) )**2  &
            + 2. * ZINTE(J1D) * PG_O_THVREF(J1D)                        &
                 * ZDELTVPT(J1D,JKK) / PDZZ2D(J1D,JKK) ))    ) /             &
        ( PG_O_THVREF(J1D) * ZDELTVPT(J1D,JKK) / PDZZ2D(J1D,JKK) ) 
!
        PLWORK(J1D)=PLWORK(J1D)+ZTEST0*(ZTEST*ZLWORK1(J1D)+  &
                                    (1-ZTEST)*ZLWORK2(J1D)) 
        ZINTE(J1D) = ZINTE(J1D) - ZPOTE(J1D)
      END DO 
    ENDIF
  END DO 
ENDIF
 
END SUBROUTINE COMPUTE_BL89_ML
