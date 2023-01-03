!      ###################################
       MODULE MODI_LIMA_COLD_SEDIMENTATION
!      ###################################
!
INTERFACE
      SUBROUTINE LIMA_COLD_SEDIMENTATION (OSEDI, KSPLITG, PTSTEP, KMI,     &
                                          HFMFILE, HLUOUT, OCLOSE_OUT,     &
                                          PZZ, PRHODJ, PRHODREF,           &
                                          PRIT, PCIT,                      &
                                          PRIS, PRSS, PRGS, PRHS, PCIS,    &
                                          PINPRS, PINPRG, PINPRH )
!
LOGICAL,                  INTENT(IN)    :: OSEDI      ! switch to activate the 
                                                      ! cloud ice sedimentation
INTEGER,                  INTENT(IN)    :: KSPLITG    ! Number of small time step 
                                                      ! for ice sedimendation
REAL,                     INTENT(IN)    :: PTSTEP     ! Time step          
INTEGER,                  INTENT(IN)    :: KMI        ! Model index 
CHARACTER(LEN=*),         INTENT(IN)    :: HFMFILE    ! Name of the output FM-file
CHARACTER(LEN=*),         INTENT(IN)    :: HLUOUT     ! Output-listing name for
                                                      ! model n
LOGICAL,                  INTENT(IN)    :: OCLOSE_OUT ! Conditional closure of 
                                                      ! the FM file output
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PZZ        ! Height (z)
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ     ! Dry density * Jacobian (Budgets)
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF   ! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRIT       ! Cloud ice m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCIT       ! Ice crystal C. at t
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRIS       ! Pristine ice m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRSS       ! Snow/aggregate m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRGS       ! Graupel m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRHS       ! Hail m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCIS       ! Ice crystal C. source
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PINPRS  ! Snow instant precip
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PINPRG  ! Graupel instant precip
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PINPRH  ! Hail instant precip

!
      END SUBROUTINE LIMA_COLD_SEDIMENTATION
END INTERFACE
END MODULE MODI_LIMA_COLD_SEDIMENTATION
!
!
!     ######################################################################
      SUBROUTINE LIMA_COLD_SEDIMENTATION (OSEDI, KSPLITG, PTSTEP, KMI,     &
                                          HFMFILE, HLUOUT, OCLOSE_OUT,     &
                                          PZZ, PRHODJ, PRHODREF,           &
                                          PRIT, PCIT,                      &
                                          PRIS, PRSS, PRGS, PRHS, PCIS,    &
                                          PINPRS,PINPRG,PINPRH                 )
!     ######################################################################
!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the sediimentation
!!    of primary ice, snow and graupel.
!!
!!    METHOD
!!    ------
!!      The sedimentation rates are computed with a time spliting technique: 
!!    an upstream scheme, written as a difference of non-advective fluxes. 
!!    This source term is added to the next coming time step (split-implicit 
!!    process).
!!
!!    AUTHOR
!!    ------
!!      J.-M. Cohard     * Laboratoire d'Aerologie*
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!      S.    Berthet    * Laboratoire d'Aerologie*
!!      B.    Vié        * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original             ??/??/13 
!!      C. Barthe  * LACy *  jan. 2014   add budgets
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAM_LIMA_COLD,  ONLY : XLBEXI, XLBI, XDI,                 &
                                  XFSEDRI, XFSEDCI, XFSEDS, XEXSEDS
USE MODD_PARAM_LIMA_MIXED, ONLY : XFSEDG, XEXSEDG, XFSEDH, XEXSEDH
USE MODD_PARAM_LIMA,       ONLY : XCEXVT, XRTMIN, XCTMIN
USE MODD_CST,              ONLY : XRHOLW
USE MODD_PARAMETERS,       ONLY : JPHEXT, JPVEXT
USE MODI_LIMA_FUNCTIONS,   ONLY : COUNTJV
!
USE MODD_NSV

IMPLICIT NONE

!
!*       0.1   Declarations of dummy arguments :
!
LOGICAL,                  INTENT(IN)    :: OSEDI      ! switch to activate the 
                                                      ! cloud ice sedimentation
INTEGER,                  INTENT(IN)    :: KSPLITG    ! Number of small time step 
                                                      ! for ice sedimendation
REAL,                     INTENT(IN)    :: PTSTEP     ! Time step          
INTEGER,                  INTENT(IN)    :: KMI        ! Model index 
CHARACTER(LEN=*),         INTENT(IN)    :: HFMFILE    ! Name of the output FM-file
CHARACTER(LEN=*),         INTENT(IN)    :: HLUOUT     ! Output-listing name for
                                                      ! model n
LOGICAL,                  INTENT(IN)    :: OCLOSE_OUT ! Conditional closure of 
                                                      ! the FM file output
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PZZ        ! Height (z)
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ     ! Dry density * Jacobian (Budgets)
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF   ! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRIT       ! Cloud ice m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCIT       ! Ice crystal C. at t
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRIS       ! Pristine ice m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRSS       ! Snow/aggregate m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRGS       ! Graupel m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRHS       ! Hail m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCIS       ! Ice crystal C. source
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PINPRS  ! Snow instant precip
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PINPRG  ! Graupel instant precip
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PINPRH  ! Hail instant precip

!
!*       0.2   Declarations of local variables :
!
INTEGER :: JK, JL, JN                     ! Loop index
INTEGER :: IIB, IIE, IJB, IJE, IKB, IKE   ! Physical domain
INTEGER :: ISEDIM                         ! Case number of sedimentation
!
LOGICAL, DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3)) &
                           :: GSEDIM      ! Test where to compute the SED processes
REAL,    DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3)) &
                           :: ZW          ! Work array
REAL,    DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3)+1) &
                           :: ZWSEDR,   & ! Sedimentation of MMR
                              ZWSEDC      ! Sedimentation of number conc.
!
REAL, DIMENSION(:), ALLOCATABLE         &
                           :: ZRIS,     & ! Pristine ice m.r. source
                              ZCIS,     & ! Pristine ice conc. source
                              ZRSS,     & ! Snow/aggregate m.r. source
                              ZRGS,     & ! Graupel/hail m.r. source
                              ZRHS,     & ! Graupel/hail m.r. source
                              ZRIT,     & ! Pristine ice m.r. at t
                              ZCIT,     & ! Pristine ice conc. at t
                              ZRHODREF, & ! RHO Dry REFerence
                              ZRHODJ,   & ! RHO times Jacobian
                              ZZW,      & ! Work array
                              ZZX,      & ! Work array
                              ZZY,      & ! Work array
                              ZLBDAI,   & ! Slope parameter of the ice crystal distr.
                              ZRTMIN 
!
INTEGER , DIMENSION(SIZE(PRHODREF)) :: I1,I2,I3 ! Indexes for PACK replacement
!
REAL    :: ZTSPLITG                       ! Small time step for rain sedimentation
!
INTEGER :: IKMAX 
!
!!!!!!!!!!! Entiers pour niveaux inversés dans AROME !!!!!!!!!!!!!!!!!!!!!!!!!!!
INTEGER :: IBOTTOM, INVLVL
!
!-------------------------------------------------------------------------------
!
! Physical domain
!
IIB=1+JPHEXT
IIE=SIZE(PZZ,1) - JPHEXT
IJB=1+JPHEXT
IJE=SIZE(PZZ,2) - JPHEXT
IKB=1+JPVEXT
IKE=SIZE(PZZ,3) - JPVEXT
!
!!!!!!!!!!! Entiers pour niveaux inversés dans AROME !!!!!!!!!!!!!!!!!!!!!!!!!!!
IBOTTOM=IKE
INVLVL=-1
!
ZWSEDR(:,:,:)=0.
ZWSEDC(:,:,:)=0.
IKMAX=SIZE(PRHODREF,3)
!
! Time splitting and ZRTMIN
!
ALLOCATE(ZRTMIN(SIZE(XRTMIN)))
ZRTMIN(:) = XRTMIN(:) / PTSTEP
!
ZTSPLITG= PTSTEP / FLOAT(KSPLITG)
!
PINPRS(:,:) = 0.
PINPRG(:,:) = 0.
PINPRH(:,:) = 0.
!
! ################################
! Compute the sedimentation fluxes
! ################################
!
DO JN = 1 , KSPLITG 
  ! Computation only where enough ice, snow, graupel or hail
   GSEDIM(:,:,:) = .FALSE.
   GSEDIM(IIB:IIE,IJB:IJE,IKB:IKE) = PRSS(IIB:IIE,IJB:IJE,IKB:IKE)>ZRTMIN(5) &
                                .OR. PRGS(IIB:IIE,IJB:IJE,IKB:IKE)>ZRTMIN(6) &
                                .OR. PRHS(IIB:IIE,IJB:IJE,IKB:IKE)>ZRTMIN(7)
   IF( OSEDI ) THEN
      GSEDIM(IIB:IIE,IJB:IJE,IKB:IKE) = GSEDIM(IIB:IIE,IJB:IJE,IKB:IKE) &
                                .OR. PRIS(IIB:IIE,IJB:IJE,IKB:IKE)>ZRTMIN(4)
   END IF
!
   ISEDIM = COUNTJV( GSEDIM(:,:,:),I1(:),I2(:),I3(:))
   IF( ISEDIM >= 1 ) THEN
!
      IF( JN==1 ) THEN
         IF( OSEDI ) THEN
            PCIS(:,:,:) = PCIS(:,:,:) * PTSTEP
            PRIS(:,:,:) = PRIS(:,:,:) * PTSTEP
         END IF
         PRSS(:,:,:) = PRSS(:,:,:) * PTSTEP
         PRGS(:,:,:) = PRGS(:,:,:) * PTSTEP
         PRHS(:,:,:) = PRHS(:,:,:) * PTSTEP
         DO JK = IKB , IKE
!Dans AROME, PZZ = épaisseur de la couche
!            ZW(:,:,JK)=ZTSPLITG/(PZZ(:,:,JK+1)-PZZ(:,:,JK))
            ZW(:,:,JK)=ZTSPLITG/(PZZ(:,:,JK))
         END DO
      END IF
!
      ALLOCATE(ZRHODREF(ISEDIM))
      DO JL = 1,ISEDIM
         ZRHODREF(JL) = PRHODREF(I1(JL),I2(JL),I3(JL))
      END DO
!
      ALLOCATE(ZZW(ISEDIM)) ; ZZW(:) = 0.0
      ALLOCATE(ZZX(ISEDIM)) ; ZZX(:) = 0.0
      ALLOCATE(ZZY(ISEDIM)) ; ZZY(:) = 0.0
!
!*       2.21   for pristine ice
!
      IF( OSEDI.AND.MAXVAL(PRIS(:,:,:))>ZRTMIN(4) ) THEN
         ALLOCATE(ZRIS(ISEDIM))
         ALLOCATE(ZCIS(ISEDIM))
         ALLOCATE(ZRIT(ISEDIM))
         ALLOCATE(ZCIT(ISEDIM))
         ALLOCATE(ZLBDAI(ISEDIM))
         DO JL = 1,ISEDIM
            ZRIS(JL) = PRIS(I1(JL),I2(JL),I3(JL))
            ZCIS(JL) = PCIS(I1(JL),I2(JL),I3(JL))
            ZRIT(JL) = PRIT(I1(JL),I2(JL),I3(JL))
            ZCIT(JL) = PCIT(I1(JL),I2(JL),I3(JL))
         END DO
         ZLBDAI(:)  = 1.E10
         WHERE (ZRIT(:)>XRTMIN(4) .AND. ZCIT(:)>XCTMIN(4))
            ZLBDAI(:) = ( XLBI*ZCIT(:) / ZRIT(:) )**XLBEXI
         END WHERE
         WHERE( ZRIS(:)>ZRTMIN(4) )
            ZZY(:) = ZRHODREF(:)**(-XCEXVT) * ZLBDAI(:)**(-XDI)
            ZZW(:) = XFSEDRI * ZRIS(:) * ZZY(:) * ZRHODREF(:)
            ZZX(:) = XFSEDCI * ZCIS(:) * ZZY(:) * ZRHODREF(:)
         END WHERE
         ZWSEDR(:,:,1:IKMAX) = UNPACK( ZZW(:),MASK=GSEDIM(:,:,:),FIELD=0.0 )
         ZWSEDC(:,:,1:IKMAX) = UNPACK( ZZX(:),MASK=GSEDIM(:,:,:),FIELD=0.0 )
         DO JK = IKB+1 , IKE
            PRIS(:,:,JK) = PRIS(:,:,JK) + ZW(:,:,JK)*    &
                 (ZWSEDR(:,:,JK+1*INVLVL)-ZWSEDR(:,:,JK))/PRHODREF(:,:,JK)
            PCIS(:,:,JK) = PCIS(:,:,JK) + ZW(:,:,JK)*    &
                 (ZWSEDC(:,:,JK+1*INVLVL)-ZWSEDC(:,:,JK))/PRHODREF(:,:,JK)
         END DO
            PRIS(:,:,1) = PRIS(:,:,1) + ZW(:,:,1)*    &
                 (0.-ZWSEDR(:,:,1))/PRHODREF(:,:,1)
            PCIS(:,:,1) = PCIS(:,:,1) + ZW(:,:,1)*    &
                 (0.-ZWSEDC(:,:,1))/PRHODREF(:,:,1)
         DEALLOCATE(ZRIS)
         DEALLOCATE(ZCIS)
         DEALLOCATE(ZRIT)
         DEALLOCATE(ZCIT)
         DEALLOCATE(ZLBDAI)
      END IF
!
!*       2.22   for aggregates
!
      ZZW(:) = 0.
      IF( MAXVAL(PRSS(:,:,:))>ZRTMIN(5) ) THEN
         ALLOCATE(ZRSS(ISEDIM)) 
         DO JL = 1,ISEDIM
            ZRSS(JL) = PRSS(I1(JL),I2(JL),I3(JL))
         END DO
         WHERE( ZRSS(:)>ZRTMIN(5) )
! Correction BVIE ZRHODREF
!            ZZW(:) = XFSEDS * ZRSS(:)**XEXSEDS * ZRHODREF(:)**(XEXSEDS-XCEXVT)
            ZZW(:) = XFSEDS * ZRSS(:)**XEXSEDS * ZRHODREF(:)**(-XCEXVT) * ZRHODREF(:)
         END WHERE
         ZWSEDR(:,:,1:IKMAX) = UNPACK( ZZW(:),MASK=GSEDIM(:,:,:),FIELD=0.0 )
         DO JK = IKB+1 , IKE
            PRSS(:,:,JK) = PRSS(:,:,JK) + ZW(:,:,JK)* &
                 (ZWSEDR(:,:,JK+1*INVLVL)-ZWSEDR(:,:,JK))/PRHODREF(:,:,JK)
         END DO
            PRSS(:,:,1) = PRSS(:,:,1) + ZW(:,:,1)* &
                 (0.-ZWSEDR(:,:,1))/PRHODREF(:,:,1)
         DEALLOCATE(ZRSS)
      ELSE
         ZWSEDR(:,:,IBOTTOM) = 0.0
      END IF
!    
      PINPRS(:,:) = PINPRS(:,:) + ZWSEDR(:,:,IBOTTOM)/XRHOLW/KSPLITG                          ! in m/s
!
!*       2.23   for graupeln
!
      ZZW(:) = 0.
      IF( MAXVAL(PRGS(:,:,:))>ZRTMIN(6) ) THEN
         ALLOCATE(ZRGS(ISEDIM)) 
         DO JL = 1,ISEDIM
            ZRGS(JL) = PRGS(I1(JL),I2(JL),I3(JL))
         END DO
         WHERE( ZRGS(:)>ZRTMIN(6) )
! Correction BVIE ZRHODREF
!            ZZW(:) = XFSEDG * ZRGS(:)**XEXSEDG * ZRHODREF(:)**(XEXSEDG-XCEXVT)
            ZZW(:) = XFSEDG * ZRGS(:)**XEXSEDG * ZRHODREF(:)**(-XCEXVT) * ZRHODREF(:)
         END WHERE
         ZWSEDR(:,:,1:IKMAX) = UNPACK( ZZW(:),MASK=GSEDIM(:,:,:),FIELD=0.0 )
         DO JK = IKB+1 , IKE
            PRGS(:,:,JK) = PRGS(:,:,JK) + ZW(:,:,JK)* &
                 (ZWSEDR(:,:,JK+1*INVLVL)-ZWSEDR(:,:,JK))/PRHODREF(:,:,JK)
         END DO
            PRGS(:,:,1) = PRGS(:,:,1) + ZW(:,:,1)* &
                 (0.-ZWSEDR(:,:,1))/PRHODREF(:,:,1)
         DEALLOCATE(ZRGS)
      ELSE
         ZWSEDR(:,:,IBOTTOM) = 0.0
      END IF
!    
      PINPRG(:,:) = PINPRG(:,:) + ZWSEDR(:,:,IBOTTOM)/XRHOLW/KSPLITG                        ! in m/s
!
!*       2.23   for hail
!
      ZZW(:) = 0.
      IF( MAXVAL(PRHS(:,:,:))>ZRTMIN(7) ) THEN
         ALLOCATE(ZRHS(ISEDIM)) 
         DO JL = 1,ISEDIM
            ZRHS(JL) = PRHS(I1(JL),I2(JL),I3(JL))
         END DO
         WHERE( ZRHS(:)>ZRTMIN(7) )
! Correction BVIE ZRHODREF
!            ZZW(:) = XFSEDH * ZRHS(:)**XEXSEDH * ZRHODREF(:)**(XEXSEDH-XCEXVT)
            ZZW(:) = XFSEDH * ZRHS(:)**XEXSEDH * ZRHODREF(:)**(-XCEXVT) * ZRHODREF(:)
         END WHERE
         ZWSEDR(:,:,1:IKMAX) = UNPACK( ZZW(:),MASK=GSEDIM(:,:,:),FIELD=0.0 )
         DO JK = IKB+1 , IKE
            PRHS(:,:,JK) = PRHS(:,:,JK) + ZW(:,:,JK)* &
                 (ZWSEDR(:,:,JK+1*INVLVL)-ZWSEDR(:,:,JK))/PRHODREF(:,:,JK)
         END DO
            PRHS(:,:,1) = PRHS(:,:,1) + ZW(:,:,1)* &
                 (0.-ZWSEDR(:,:,1))/PRHODREF(:,:,1)
         DEALLOCATE(ZRHS)
      ELSE
         ZWSEDR(:,:,IBOTTOM) = 0.0
      END IF
!    
      PINPRH(:,:) = PINPRH(:,:) + ZWSEDR(:,:,IBOTTOM)/XRHOLW/KSPLITG                        ! in m/s
!
!*       2.24 End of sedimentation  
!
      DEALLOCATE(ZRHODREF)
      DEALLOCATE(ZZW)
      DEALLOCATE(ZZX)
      DEALLOCATE(ZZY)
      IF( JN==KSPLITG ) THEN
         IF( OSEDI ) THEN
            PRIS(:,:,:) = PRIS(:,:,:) / PTSTEP
            PCIS(:,:,:) = PCIS(:,:,:) / PTSTEP
         END IF
         PRSS(:,:,:) = PRSS(:,:,:) / PTSTEP
         PRGS(:,:,:) = PRGS(:,:,:) / PTSTEP
         PRHS(:,:,:) = PRHS(:,:,:) / PTSTEP
      END IF
   END IF
END DO
!
DEALLOCATE(ZRTMIN)
!
END SUBROUTINE LIMA_COLD_SEDIMENTATION
!
!-------------------------------------------------------------------------------
