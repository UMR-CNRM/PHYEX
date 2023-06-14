!     ######spl
      SUBROUTINE CONVECT_CHEM_TRANSPORT( CVPEXT, D, NSV, KCH, PCH1, PCH1C, &
                                         KDPL, KPBL, KLCL, KCTL, KLFS, KDBL, &
                                         PUMF, PUER, PUDR, PDMF, PDER, PDDR, &
                                         PTIMEC, PDXDY, PMIXF, PLMASS, PWSUB,&
                                         KFTSTEPS )
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
!     ########################################################################
!
!!**** Compute  modified chemical tracer values due to convective event
!!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to determine the final adjusted
!!      environmental values of the chemical tracers
!!      The final convective tendencies can then be evaluated in the main
!!      routine DEEP_CONVECT by (PCH1C-PCH1)/PTIMEC
!!
!!
!!**  METHOD
!!    ------
!!      Identical to the computation of the conservative variables in the
!!      main deep convection code
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST
!!          XG                 ! gravity constant
!!
!!     Module MODD_CONVPAREXT
!!          JCVEXB, JCVEXT     ! extra levels on the vertical boundaries
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      Original    11/12/97
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST, ONLY : XG
USE MODD_CONVPAREXT, ONLY : CONVPAREXT
USE MODD_NSV,  ONLY : NSV_T
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_T
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(CONVPAREXT),       INTENT(IN) :: CVPEXT
TYPE(DIMPHYEX_T),       INTENT(IN) :: D
TYPE(NSV_T),            INTENT(IN) :: NSV
INTEGER,                INTENT(IN) :: KCH      ! number of passive tracers
!
REAL,DIMENSION(D%NIT,D%NKT,KCH),INTENT(IN) :: PCH1 ! grid scale tracer concentr.
REAL,DIMENSION(D%NIT,D%NKT,KCH),INTENT(OUT):: PCH1C! conv adjusted tracer concntr.
!
INTEGER, DIMENSION(D%NIT), INTENT(IN) :: KDPL   ! index for departure level
INTEGER, DIMENSION(D%NIT), INTENT(IN) :: KPBL   ! index for top of source layer
INTEGER, DIMENSION(D%NIT), INTENT(IN) :: KLCL   ! index lifting condens. level
INTEGER, DIMENSION(D%NIT), INTENT(IN) :: KCTL   ! index for cloud top level
INTEGER, DIMENSION(D%NIT), INTENT(IN) :: KLFS   ! index for level of free sink
INTEGER, DIMENSION(D%NIT), INTENT(IN) :: KDBL   ! index for downdraft base level
!
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN) :: PUMF ! updraft mass flux (kg/s)
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN) :: PUER ! updraft entrainment (kg/s)
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN) :: PUDR ! updraft detrainment (kg/s)
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN) :: PDMF ! downdraft mass flux (kg/s)
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN) :: PDER ! downdraft entrainment (kg/s)
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN) :: PDDR ! downdraft detrainment (kg/s)
!
REAL, DIMENSION(D%NIT),     INTENT(IN) :: PTIMEC! convection time step
REAL,                      INTENT(IN) :: PDXDY ! grid area (m^2)
REAL, DIMENSION(D%NIT),     INTENT(IN) :: PMIXF ! mixed fraction at LFS
REAL, DIMENSION(D%NIT,D%NKT),INTENT(IN) :: PLMASS! mass of model layer (kg)
REAL, DIMENSION(D%NIT,D%NKT),INTENT(IN) :: PWSUB ! envir. compensating subsidence(Pa/s)
INTEGER,                INTENT(IN) :: KFTSTEPS  ! maximum fractional time steps
!
!
!*       0.2   Declarations of local variables :
!
INTEGER :: INCH1          ! number of chemical tracers
INTEGER :: IKB, IKE  ! horizontal + vertical loop bounds
INTEGER :: IKS            ! vertical dimension
INTEGER :: JI             ! horizontal loop index
INTEGER :: JK, JKP        ! vertical loop index
INTEGER :: JN             ! chemical tracer loop index
INTEGER :: JSTEP          ! fractional time loop index
INTEGER :: JKLD, JKLP, JKMIN, JKMAX, JKMAX2 ! loop index for levels
!
REAL, DIMENSION(D%NIT,D%NKT)     :: ZOMG ! compensat. subsidence (Pa/s)
REAL, DIMENSION(D%NIT,D%NKT,KCH) :: ZUCH1, ZDCH1 ! updraft/downdraft values
REAL, DIMENSION(D%NIT)          :: ZTIMEC  ! fractional convective time step
REAL, DIMENSION(D%NIT,D%NKT)     :: ZTIMC! 2D work array for ZTIMEC
REAL, DIMENSION(D%NIT,D%NKT,KCH) :: ZCH1MFIN, ZCH1MFOUT
                                   ! work arrays for environm. compensat. mass
REAL, DIMENSION(D%NIT,KCH)      :: ZWORK1, ZWORK2, ZWORK3
!
!-------------------------------------------------------------------------------
!
!*       0.3   Compute loop bounds
!              -------------------
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('CONVECT_CHEM_TRANSPORT',0,ZHOOK_HANDLE)
INCH1  = KCH
IKB    = 1 + CVPEXT%JCVEXB
IKS    = D%NKT
IKE    = D%NKT - CVPEXT%JCVEXT
JKMAX  = 0
JKMIN  = 999999999
DO JI=D%NIB, D%NIE
  JKMIN = MIN(JKMIN, KDPL(JI))
  JKMAX = MAX(JKMAX, KCTL(JI))
ENDDO

!
!
!*      2.      Updraft computations
!               --------------------
!
ZUCH1(:,:,:) = 0.
!
!*      2.1     Initialization  at LCL
!               ----------------------------------
!
DO JI = D%NIB, D%NIE
    DO JN = 1, INCH1
      JKLD = KDPL(JI)
      JKLP = KPBL(JI)
      ZWORK1(JI,JN) = .5 * ( PCH1(JI,JKLD,JN) + PCH1(JI,JKLP,JN) )
    ENDDO
END DO
!
!*      2.2     Final updraft loop
!               ------------------
!
DO JK = JKMIN, JKMAX
JKP = JK + 1
!
    DO JN = 1, INCH1
     DO JI = D%NIB, D%NIE
       IF ( KDPL(JI) <= JK .AND. MIN(KLCL(JI), KCTL(JI)) > JK )                             &
            ZUCH1(JI,JK,JN) = ZWORK1(JI,JN)
!
       IF ( MIN(KLCL(JI), KCTL(JI)) - 1 <= JK .AND. KCTL(JI) > JK ) THEN
           ZUCH1(JI,JKP,JN) = ZUCH1(JI,JK,JN)
                            !if you have reactive i.e. non-passive tracers
                            ! update their values here and add the corresponding
                            ! sink term in the following equation
           ZUCH1(JI,JKP,JN) = ( PUMF(JI,JK) * ZUCH1(JI,JK,JN) +              &
                                PUER(JI,JKP) * PCH1(JI,JK,JN) )  /           &
                              ( PUMF(JI,JKP) + PUDR(JI,JKP) + 1.E-7 )
       END IF
     END DO
   END DO
!
END DO
!
!*      3.      Downdraft computations
!               ----------------------
!
ZDCH1(:,:,:) = 0.
!
!*      3.1     Initialization at the LFS
!               -------------------------
!
DO JI=1,D%NIE
DO JK=1,INCH1
ZWORK1(JI,JK) = PMIXF(JI)
ENDDO
ENDDO
DO JI = D%NIB, D%NIE
    DO JN = 1, INCH1
     JK = KLFS(JI)
     ZDCH1(JI,JK,JN) = ZWORK1(JI,JN) * PCH1(JI,JK,JN) +                          &
                                       ( 1. - ZWORK1(JI,JN) ) * ZUCH1(JI,JK,JN)
    ENDDO
END DO
!
!*      3.2     Final downdraft loop
!               --------------------
!
JKMAX2 = 0
DO JI=D%NIB, D%NIE
  JKMAX2 = MAX(JKMAX2, KLFS(JI))
ENDDO
DO JK = JKMAX2, IKB + 1, -1
JKP = JK - 1
    DO JN = 1, INCH1
    DO JI = D%NIB, D%NIE
      IF ( JK <= KLFS(JI) .AND. JKP >= KDBL(JI) ) THEN
       ZDCH1(JI,JKP,JN) = ( ZDCH1(JI,JK,JN) * PDMF(JI,JK) -              &
                            PCH1(JI,JK,JN) *  PDER(JI,JKP) ) /           &
                          ( PDMF(JI,JKP) - PDDR(JI,JKP) - 1.E-7 )
      END IF
    END DO
    END DO
END DO
!
!
!*      4.      Final closure (environmental) computations
!               ------------------------------------------
!
PCH1C(D%NIB:D%NIE,IKB:IKE,1:KCH) = PCH1(D%NIB:D%NIE,IKB:IKE,1:KCH) ! initialize adjusted envir. values
!
DO JK = IKB, IKE
  DO JI=D%NIB,D%NIE
    ZOMG(JI,JK) = PWSUB(JI,JK) * PDXDY / XG ! environmental subsidence
  ENDDO
END DO
!
DO JI=D%NIB,D%NIE
  ZTIMEC(JI) = PTIMEC(JI) / REAL( KFTSTEPS ) ! adjust  fractional time step
ENDDO                                    ! to be an integer multiple of PTIMEC

DO JI=D%NIB,D%NIE
  IF(PTIMEC(JI) < 1.) ZTIMEC(JI) = 0
ENDDO
DO JI=1,D%NIE
DO JK=1,IKS
  ZTIMC(JI,JK) = ZTIMEC(JI)
ENDDO
ENDDO
!
ZCH1MFIN(D%NIB:D%NIE,1:D%NKT,1:KCH)   = 0.
ZCH1MFOUT(D%NIB:D%NIE,1:D%NKT,1:KCH)  = 0.
!
DO JSTEP = 1, KFTSTEPS ! Enter the fractional time step loop
!
      DO JK = IKB + 1, JKMAX
          JKP = MAX( IKB + 1, JK - 1 )
          DO JI=1,D%NIE
          DO JN=1,INCH1
            ZWORK3(JI,JN) = ZOMG(JI,JK)
          ENDDO
          ENDDO
          ZWORK1(:,:) = SIGN( 1., ZWORK3(:,:) )
          ZWORK2(:,:) = 0.5 * ( 1. + ZWORK1(:,:) )
          ZWORK1(:,:) = 0.5 * ( 1. - ZWORK1(:,:) )
          ZCH1MFIN(:,JK,:)  = - ZWORK3(:,:) * PCH1C(:,JKP,:) * ZWORK1(:,:)
          ZCH1MFOUT(:,JK,:) =   ZWORK3(:,:) * PCH1C(:,JK,:)  * ZWORK2(:,:)
         ZCH1MFIN(:,JKP,:) = ZCH1MFIN(:,JKP,:) + ZCH1MFOUT(:,JK,:) * ZWORK2(:,:)
         ZCH1MFOUT(:,JKP,:)= ZCH1MFOUT(:,JKP,:) + ZCH1MFIN(:,JK,:) * ZWORK1(:,:)
      END DO
!
      DO JN = 1, INCH1
       DO JK = IKB + 1, JKMAX
         DO JI=D%NIB,D%NIE
           PCH1C(JI,JK,JN) = PCH1C(JI,JK,JN) + ZTIMC(JI,JK) / PLMASS(JI,JK) *  (    &
                        ZCH1MFIN(JI,JK,JN) + PUDR(JI,JK) * ZUCH1(JI,JK,JN) +       &
                        PDDR(JI,JK) * ZDCH1(JI,JK,JN) - ZCH1MFOUT(JI,JK,JN) -      &
                        ( PUER(JI,JK) + PDER(JI,JK) ) * PCH1(JI,JK,JN)    )
           IF(JN < NSV%NSV_LGBEG .OR. JN>NSV%NSV_LGEND-1) THEN
             PCH1C(JI,JK,JN) = MAX( 0., PCH1C(JI,JK,JN) )
           ELSE
             ! no tendency for horizontal Lagrangian variables
             PCH1C(JI,JK,JN) = PCH1(JI,JK,JN)
           END IF
         ENDDO
       END DO
      END DO
!
END DO ! Exit the fractional time step loop
!
!
IF (LHOOK) CALL DR_HOOK('CONVECT_CHEM_TRANSPORT',1,ZHOOK_HANDLE)
END SUBROUTINE CONVECT_CHEM_TRANSPORT
