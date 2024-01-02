MODULE MODE_COMPUTE_BL89_ML
IMPLICIT NONE
CONTAINS
!     ######spl
      SUBROUTINE COMPUTE_BL89_ML(D, CST, CSTURB,PDZZ2D, &
             PTKEM_DEP,PG_O_THVREF,PVPT,KK,OUPORDN,OFLUX,PSHEAR,PLWORK)

      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
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
!  P. Wautelet 17/12/2021: bugfix: KK instead of JKK
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
USE MODD_DIMPHYEX,        ONLY: DIMPHYEX_t
USE MODD_CST, ONLY: CST_t
USE MODD_CTURB, ONLY: CSTURB_t
!
USE MODE_MSG
!
USE MODI_SHUMAN_MF, ONLY: DZM_MF, MZM_MF
!
IMPLICIT NONE
!
!          0.1 arguments
!
TYPE(DIMPHYEX_t),       INTENT(IN)   :: D
TYPE(CST_t),            INTENT(IN)   :: CST
TYPE(CSTURB_t),         INTENT(IN)   :: CSTURB
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)  :: PDZZ2D        ! height difference between two mass levels
REAL, DIMENSION(D%NIJT),     INTENT(IN)  :: PTKEM_DEP     ! TKE to consume
REAL, DIMENSION(D%NIJT),     INTENT(IN)  :: PG_O_THVREF   ! g/ThetaVRef at the departure point
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)  :: PVPT          ! ThetaV on mass levels
INTEGER,                INTENT(IN)  :: KK            ! index of departure level
LOGICAL,                INTENT(IN)  :: OUPORDN       ! switch to compute upward (true) or
                                                     !   downward (false) mixing length
LOGICAL,                INTENT(IN)  :: OFLUX         ! Computation must be done from flux level
REAL, DIMENSION(D%NIJT),     INTENT(OUT) :: PLWORK        ! Resulting mixing length
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)  :: PSHEAR        ! vertical wind shear for RM17 mixing length

!          0.2 Local variable
!
REAL, DIMENSION(D%NIJT) :: ZLWORK1,ZLWORK2 ! Temporary mixing length
REAL, DIMENSION(D%NIJT) :: ZINTE,ZPOTE     ! TKE and potential energy
                                                 !   between 2 levels
REAL, DIMENSION(D%NIJT) :: ZVPT_DEP        ! Thetav on departure point
!
REAL, DIMENSION(D%NIJT,D%NKT) :: ZDELTVPT,ZHLVPT                                
                      !Virtual Potential Temp at Half level and DeltaThv between
                      !2 mass levels

INTEGER :: JKK                  !loop counters
INTEGER :: JIJ, JK
INTEGER :: IIJB,IIJE ! physical horizontal domain indices
INTEGER :: IKT,IKB,IKA,IKE,IKL
REAL    :: ZTEST,ZTEST0,ZTESTM  !test for vectorization
!-------------------------------------------------------------------------------------
!
!*       1.    INITIALISATION
!              --------------
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('COMPUTE_BL89_ML',0,ZHOOK_HANDLE)
!
IIJE=D%NIJE
IIJB=D%NIJB
IKT=D%NKT
IKB=D%NKB
IKA=D%NKA
IKE=D%NKE
IKL=D%NKL
!
CALL DZM_MF(D, PVPT(:,:), ZDELTVPT(:,:))
ZDELTVPT(:,IKA)=0.
!$mnh_expand_where(JIJ=IIJB:IIJE,JK=1:IKT)
WHERE (ABS(ZDELTVPT(:,:))<CSTURB%XLINF)
  ZDELTVPT(:,:)=CSTURB%XLINF
END WHERE
!$mnh_end_expand_where(JIJ=IIJB:IIJE,JK=1:IKT)
!
CALL MZM_MF(D, PVPT(:,:), ZHLVPT(:,:))
!
!We consider that gradient between mass levels KKB and KKB+KKL is the same as
!the gradient between flux level KKB and mass level KKB
!$mnh_expand_array(JIJ=IIJB:IIJE)
ZDELTVPT(:,IKB)=PDZZ2D(:,IKB)*ZDELTVPT(:,IKB+IKL)/PDZZ2D(:,IKB+IKL)
ZHLVPT(:,IKB)=PVPT(:,IKB)-ZDELTVPT(:,IKB)*0.5
!$mnh_end_expand_array(JIJ=IIJB:IIJE)
!
!
!
!*       2.    CALCULATION OF THE UPWARD MIXING LENGTH
!              ---------------------------------------
!

IF (OUPORDN.EQV..TRUE.) THEN 
 !$mnh_expand_array(JIJ=IIJB:IIJE)
 ZINTE(:)=PTKEM_DEP(:)
 !$mnh_end_expand_array(JIJ=IIJB:IIJE)
 PLWORK=0.
 ZTESTM=1.
 IF(OFLUX)THEN
   !$mnh_expand_array(JIJ=IIJB:IIJE)
   ZVPT_DEP(:)=ZHLVPT(:,KK) ! departure point is on flux level
   !$mnh_end_expand_array(JIJ=IIJB:IIJE)
   !We must compute what happens between flux level KK and mass level KK
   DO JIJ=IIJB,IIJE
     ZTEST0=0.5+SIGN(0.5,ZINTE(JIJ)) ! test if there's energy to consume
     ! Energy consumed if parcel cross the entire layer
     ZPOTE(JIJ) = ZTEST0*(PG_O_THVREF(JIJ)      *      &
         (0.5*(ZHLVPT(JIJ,KK)+ PVPT(JIJ,KK)) - ZVPT_DEP(JIJ)) + &
         CSTURB%XRM17*PSHEAR(JIJ,KK)*SQRT(ABS(PTKEM_DEP(JIJ))))  * &
         PDZZ2D(JIJ,KK)*0.5
     ! Test if it rests some energy to consume
     ZTEST =0.5+SIGN(0.5,ZINTE(JIJ)-ZPOTE(JIJ))
     ! Length travelled by parcel if it rests energy to consume
     ZLWORK1(JIJ)=PDZZ2D(JIJ,KK)*0.5
     ! Lenght travelled by parcel to nullify energy
     ZLWORK2(JIJ)=        ( - PG_O_THVREF(JIJ) *                     &
            (  ZHLVPT(JIJ,KK) - ZVPT_DEP(JIJ) )                              &
            - CSTURB%XRM17*PSHEAR(JIJ,KK)*SQRT(ABS(PTKEM_DEP(JIJ))) &
          + SQRT (ABS(                                                       &
            (CSTURB%XRM17*PSHEAR(JIJ,KK)*SQRT(ABS(PTKEM_DEP(JIJ))) +  &
             PG_O_THVREF(JIJ) * (ZHLVPT(JIJ,KK) - ZVPT_DEP(JIJ)) )**2  &
            + 2. * ZINTE(JIJ) * PG_O_THVREF(JIJ)                        &
                 * ZDELTVPT(JIJ,KK) / PDZZ2D(JIJ,KK) ))    ) /             &
        ( PG_O_THVREF(JIJ) * ZDELTVPT(JIJ,KK) / PDZZ2D(JIJ,KK) ) 
      ! Effective length travelled by parcel
      PLWORK(JIJ)=PLWORK(JIJ)+ZTEST0*(ZTEST*ZLWORK1(JIJ)+  &
                                    (1-ZTEST)*ZLWORK2(JIJ))
      ! Rest of energy to consume
      ZINTE(JIJ) = ZINTE(JIJ) - ZPOTE(JIJ)
   ENDDO
 ELSE
   !$mnh_expand_array(JIJ=IIJB:IIJE)
   ZVPT_DEP(:)=PVPT(:,KK) ! departure point is on mass level
   !$mnh_end_expand_array(JIJ=IIJB:IIJE)
 ENDIF

 DO JKK=KK+IKL,IKE,IKL
    IF(ZTESTM > 0.) THEN
      ZTESTM=0
      DO JIJ=IIJB,IIJE
        ZTEST0=0.5+SIGN(0.5,ZINTE(JIJ))
        ZPOTE(JIJ) = ZTEST0*(PG_O_THVREF(JIJ)      *      &
            (ZHLVPT(JIJ,JKK) - ZVPT_DEP(JIJ))   &
           + CSTURB%XRM17*PSHEAR(JIJ,JKK)*SQRT(ABS(PTKEM_DEP(JIJ))))* PDZZ2D(JIJ,JKK) 
        ZTEST =0.5+SIGN(0.5,ZINTE(JIJ)-ZPOTE(JIJ))
        ZTESTM=ZTESTM+ZTEST0
        ZLWORK1(JIJ)=PDZZ2D(JIJ,JKK)
        !ZLWORK2 jump of the last reached level
        ZLWORK2(JIJ)=        ( - PG_O_THVREF(JIJ) *                     &
            (  PVPT(JIJ,JKK-IKL) - ZVPT_DEP(JIJ) )                              &
            - CSTURB%XRM17*PSHEAR(JIJ,JKK)*sqrt(abs(PTKEM_DEP(JIJ))) &
          + SQRT (ABS(                                                   &
            (CSTURB%XRM17*PSHEAR(JIJ,JKK)*sqrt(abs(PTKEM_DEP(JIJ))) +  &
             PG_O_THVREF(JIJ) * (PVPT(JIJ,JKK-IKL) - ZVPT_DEP(JIJ)) )**2  &
            + 2. * ZINTE(JIJ) * PG_O_THVREF(JIJ)                        &
                 * ZDELTVPT(JIJ,JKK) / PDZZ2D(JIJ,JKK) ))    ) /             &
        ( PG_O_THVREF(JIJ) * ZDELTVPT(JIJ,JKK) / PDZZ2D(JIJ,JKK) ) 
      !
        PLWORK(JIJ)=PLWORK(JIJ)+ZTEST0*(ZTEST*ZLWORK1(JIJ)+  &
                                    (1-ZTEST)*ZLWORK2(JIJ))
        ZINTE(JIJ) = ZINTE(JIJ) - ZPOTE(JIJ)
      END DO 
    ENDIF
  END DO 
ENDIF
!!
!*       2.    CALCULATION OF THE DOWNWARD MIXING LENGTH
!              ---------------------------------------
!

IF (OUPORDN.EQV..FALSE.) THEN 
 IF(OFLUX) CALL PRINT_MSG(NVERB_FATAL,'GEN','COMPUTE_BL89_ML','OFLUX option not coded for downward mixing length')
 !$mnh_expand_array(JIJ=IIJB:IIJE)
 ZINTE(:)=PTKEM_DEP(:)
 !$mnh_end_expand_array(JIJ=IIJB:IIJE)
 PLWORK=0.
 ZTESTM=1.
 DO JKK=KK,IKB,-IKL
    IF(ZTESTM > 0.) THEN
      ZTESTM=0
      DO JIJ=IIJB,IIJE
        ZTEST0=0.5+SIGN(0.5,ZINTE(JIJ))
         ZPOTE(JIJ) = ZTEST0*(-PG_O_THVREF(JIJ)      *      &
            (ZHLVPT(JIJ,JKK) - PVPT(JIJ,KK)) &
         + CSTURB%XRM17*PSHEAR(JIJ,JKK)*SQRT(ABS(PTKEM_DEP(JIJ))))* PDZZ2D(JIJ,JKK) 
        ZTEST =0.5+SIGN(0.5,ZINTE(JIJ)-ZPOTE(JIJ))
        ZTESTM=ZTESTM+ZTEST0
        ZLWORK1(JIJ)=PDZZ2D(JIJ,JKK)
        ZLWORK2(JIJ)=        ( + PG_O_THVREF(JIJ) *                     &
            (  PVPT(JIJ,JKK) - PVPT(JIJ,KK) )                              &
             -CSTURB%XRM17*PSHEAR(JIJ,JKK)*sqrt(abs(PTKEM_DEP(JIJ))) &
          + SQRT (ABS(                                                       &
            (CSTURB%XRM17*PSHEAR(JIJ,JKK)*sqrt(abs(PTKEM_DEP(JIJ))) - &
             PG_O_THVREF(JIJ) * (PVPT(JIJ,JKK) - PVPT(JIJ,KK)) )**2  &
            + 2. * ZINTE(JIJ) * PG_O_THVREF(JIJ)                        &
                 * ZDELTVPT(JIJ,JKK) / PDZZ2D(JIJ,JKK) ))    ) /             &
        ( PG_O_THVREF(JIJ) * ZDELTVPT(JIJ,JKK) / PDZZ2D(JIJ,JKK) ) 
      !
        PLWORK(JIJ)=PLWORK(JIJ)+ZTEST0*(ZTEST*ZLWORK1(JIJ)+  &
                                    (1-ZTEST)*ZLWORK2(JIJ)) 
        ZINTE(JIJ) = ZINTE(JIJ) - ZPOTE(JIJ)
      END DO 
    ENDIF
  END DO 
ENDIF
  
IF (LHOOK) CALL DR_HOOK('COMPUTE_BL89_ML',1,ZHOOK_HANDLE)
END SUBROUTINE COMPUTE_BL89_ML
END MODULE MODE_COMPUTE_BL89_ML
