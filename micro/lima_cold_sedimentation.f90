!MNH_LIC Copyright 2013-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!      ###################################
       MODULE MODI_LIMA_COLD_SEDIMENTATION
!      ###################################
!
INTERFACE
      SUBROUTINE LIMA_COLD_SEDIMENTATION (OSEDI, KSPLITG, PTSTEP, KMI,     &
                                          PZZ, PRHODJ, PRHODREF,           &
                                          PRIT, PCIT,                      &
                                          PRIS, PRSS, PRGS, PRHS, PCIS,    &
                                          PINPRS, PINPRG, PINPRH,          &
                                          PCSS, PCGS, PCHS)
!
LOGICAL,                  INTENT(IN)    :: OSEDI      ! switch to activate the 
                                                      ! cloud ice sedimentation
INTEGER,                  INTENT(IN)    :: KSPLITG    ! Number of small time step 
                                                      ! for ice sedimendation
REAL,                     INTENT(IN)    :: PTSTEP     ! Time step          
INTEGER,                  INTENT(IN)    :: KMI        ! Model index 
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
REAL, DIMENSION(:,:,:),OPTIONAL,   INTENT(INOUT) :: PCSS       ! Snow/aggregate C. source   
REAL, DIMENSION(:,:,:),OPTIONAL,   INTENT(INOUT) :: PCGS       ! Graupel C. source       
REAL, DIMENSION(:,:,:),OPTIONAL,   INTENT(INOUT) :: PCHS       ! Hail C. source
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
                                          PZZ, PRHODJ, PRHODREF,           &
                                          PRIT, PCIT,                      &
                                          PRIS, PRSS, PRGS, PRHS, PCIS,    &
                                          PINPRS,PINPRG,PINPRH,            &
                                          PCSS, PCGS, PCHS)
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
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!  P. Wautelet 26/04/2019: replace non-standard FLOAT function by REAL function
!  P. Wautelet 28/05/2019: move COUNTJV function to tools.f90
!  M. Taufour     07/2022: add concentration for snow, graupel, hail
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST,              ONLY : XRHOLW
USE MODD_NSV
USE MODD_PARAMETERS,       ONLY : JPHEXT, JPVEXT
USE MODD_PARAM_LIMA,       ONLY : XCEXVT, XRTMIN, XCTMIN, &
                                  NMOM_S, NMOM_G, NMOM_H
USE MODD_PARAM_LIMA_COLD,  ONLY : XLBEXI, XLBI, XDI,                 &
                                  XFSEDRI, XFSEDCI, XFSEDS, XEXSEDS, &
                                  XLBEXS, XLBS, XDS,                 &
                                  XFSEDRS, XFSEDCS 
                               
USE MODD_PARAM_LIMA_MIXED, ONLY : XFSEDG, XEXSEDG, XFSEDH, XEXSEDH, &
                                  XLBEXG, XLBG, XDG,                &
                                  XLBEXH, XLBH, XDH,                &
                                  XFSEDRG, XFSEDCG, XFSEDRH, XFSEDCH
!
use mode_tools,           only: Countjv
!
IMPLICIT NONE
!

!
!*       0.1   Declarations of dummy arguments :
!
LOGICAL,                  INTENT(IN)    :: OSEDI      ! switch to activate the 
                                                      ! cloud ice sedimentation
INTEGER,                  INTENT(IN)    :: KSPLITG    ! Number of small time step 
                                                      ! for ice sedimendation
REAL,                     INTENT(IN)    :: PTSTEP     ! Time step          
INTEGER,                  INTENT(IN)    :: KMI        ! Model index 
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
REAL, DIMENSION(:,:,:),OPTIONAL,   INTENT(INOUT) :: PCSS       ! Snow/aggregate C. source  
REAL, DIMENSION(:,:,:),OPTIONAL,   INTENT(INOUT) :: PCGS       ! Graupel C. source       
REAL, DIMENSION(:,:,:),OPTIONAL,   INTENT(INOUT) :: PCHS       ! Hail C. source
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
                           :: ZW,       & ! Work array
                              ZWSEDR,   & ! Sedimentation of MMR
                              ZWSEDC      ! Sedimentation of number conc.
!
REAL, DIMENSION(:), ALLOCATABLE         &
                           :: ZRIS,     & ! Pristine ice m.r. source
                              ZCIS,     & ! Pristine ice conc. source
                              ZRSS,     & ! Snow/aggregate m.r. source
                              ZCSS,     & ! Snow/aggregate conc. source
                              ZRGS,     & ! Graupel/hail m.r. source
                              ZCGS,     & ! Graupel/hail conc. source
                              ZRHS,     & ! Graupel/hail m.r. source
                              ZCHS,     & ! Graupel/hail conc. source 
                              ZRIT,     & ! Pristine ice m.r. at t
                              ZCIT,     & ! Pristine ice conc. at t
                              ZRHODREF, & ! RHO Dry REFerence
                              ZRHODJ,   & ! RHO times Jacobian
                              ZZW,      & ! Work array
                              ZZX,      & ! Work array
                              ZZY,      & ! Work array
                              ZLBDAI, ZLBDAS, ZLBDAG, ZLBDAH,   & ! Slope parameter of the ice crystal distr.
                              ZRTMIN, ZCTMIN 
!
INTEGER , DIMENSION(SIZE(PRHODREF)) :: I1,I2,I3 ! Indexes for PACK replacement
!
REAL    :: ZTSPLITG                       ! Small time step for rain sedimentation
!
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
! Time splitting and ZRTMIN
!
ALLOCATE(ZRTMIN(SIZE(XRTMIN)))
ALLOCATE(ZCTMIN(SIZE(XCTMIN)))
ZRTMIN(:) = XRTMIN(:) / PTSTEP
ZCTMIN(:) = XCTMIN(:) / PTSTEP
!
ZTSPLITG= PTSTEP / REAL(KSPLITG)
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
IF (NMOM_S.EQ.2) THEN
   GSEDIM(IIB:IIE,IJB:IJE,IKB:IKE) = PRSS(IIB:IIE,IJB:IJE,IKB:IKE)>ZRTMIN(5) .AND. &
                                     PCSS(IIB:IIE,IJB:IJE,IKB:IKE)>ZCTMIN(5)
ELSE
   GSEDIM(IIB:IIE,IJB:IJE,IKB:IKE) = PRSS(IIB:IIE,IJB:IJE,IKB:IKE)>ZRTMIN(5)
END IF
IF (NMOM_G.EQ.2) THEN
   GSEDIM(IIB:IIE,IJB:IJE,IKB:IKE) = GSEDIM(IIB:IIE,IJB:IJE,IKB:IKE) .OR. &
                                   ( PRGS(IIB:IIE,IJB:IJE,IKB:IKE)>ZRTMIN(6) .AND. &
                                     PCGS(IIB:IIE,IJB:IJE,IKB:IKE)>ZCTMIN(6) )
ELSE
   GSEDIM(IIB:IIE,IJB:IJE,IKB:IKE) = GSEDIM(IIB:IIE,IJB:IJE,IKB:IKE) .OR. &
                                     PRGS(IIB:IIE,IJB:IJE,IKB:IKE)>ZRTMIN(6)
END IF
IF (NMOM_H.EQ.2) THEN
   GSEDIM(IIB:IIE,IJB:IJE,IKB:IKE) = GSEDIM(IIB:IIE,IJB:IJE,IKB:IKE) .OR. &
                                   ( PRSS(IIB:IIE,IJB:IJE,IKB:IKE)>ZRTMIN(7) .AND. &
                                     PCSS(IIB:IIE,IJB:IJE,IKB:IKE)>ZCTMIN(7) )
ELSE
   GSEDIM(IIB:IIE,IJB:IJE,IKB:IKE) = GSEDIM(IIB:IIE,IJB:IJE,IKB:IKE) .OR. &
                                     PRSS(IIB:IIE,IJB:IJE,IKB:IKE)>ZRTMIN(7)
END IF
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
         IF(NMOM_S.EQ.2) PCSS(:,:,:) = PCSS(:,:,:) * PTSTEP
         IF(NMOM_G.EQ.2) PCGS(:,:,:) = PCGS(:,:,:) * PTSTEP
         IF(NMOM_H.EQ.2) PCHS(:,:,:) = PCHS(:,:,:) * PTSTEP
         DO JK = IKB , IKE
            ZW(:,:,JK)=ZTSPLITG/(PZZ(:,:,JK+1)-PZZ(:,:,JK))
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
         WHERE (ZRIS(:)>XRTMIN(4) .AND. ZCIS(:)>XCTMIN(4))
            ZLBDAI(:) = ( XLBI*ZCIS(:) / ZRIS(:) )**XLBEXI
            ZZY(:) = ZRHODREF(:)**(-XCEXVT) * ZLBDAI(:)**(-XDI)
            ZZW(:) = XFSEDRI * ZRIS(:) * ZZY(:) * ZRHODREF(:)
            ZZX(:) = XFSEDCI * ZCIS(:) * ZZY(:) * ZRHODREF(:)
         END WHERE
         ZWSEDR(:,:,:) = UNPACK( ZZW(:),MASK=GSEDIM(:,:,:),FIELD=0.0 )
         ZWSEDR(:,:,IKB:IKE) = MIN( ZWSEDR(:,:,IKB:IKE), PRIS(:,:,IKB:IKE) * PRHODREF(:,:,IKB:IKE) / ZW(:,:,IKB:IKE) )
         ZWSEDC(:,:,:) = UNPACK( ZZX(:),MASK=GSEDIM(:,:,:),FIELD=0.0 )
         ZWSEDC(:,:,IKB:IKE) = MIN( ZWSEDC(:,:,IKB:IKE), PCIS(:,:,IKB:IKE) * PRHODREF(:,:,IKB:IKE) / ZW(:,:,IKB:IKE) )
         DO JK = IKB , IKE
            PRIS(:,:,JK) = PRIS(:,:,JK) + ZW(:,:,JK)*    &
                 (ZWSEDR(:,:,JK+1)-ZWSEDR(:,:,JK))/PRHODREF(:,:,JK)
            PCIS(:,:,JK) = PCIS(:,:,JK) + ZW(:,:,JK)*    &
                 (ZWSEDC(:,:,JK+1)-ZWSEDC(:,:,JK))/PRHODREF(:,:,JK)
         END DO
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
      ZZX(:) = 0.0                                                
      ZZY(:) = 0.0
      IF( MAXVAL(PRSS(:,:,:))>XRTMIN(5) ) THEN
         ALLOCATE(ZRSS(ISEDIM)) 
         IF(NMOM_S.GE.2) THEN
            ALLOCATE(ZCSS(ISEDIM))  
            ALLOCATE(ZLBDAS(ISEDIM))
            DO JL = 1,ISEDIM
               ZRSS(JL) = PRSS(I1(JL),I2(JL),I3(JL))
               ZCSS(JL) = PCSS(I1(JL),I2(JL),I3(JL))       
            END DO
            ZLBDAS(:)  = 1.E10         
            WHERE( ZRSS(:)>XRTMIN(5) .AND. ZCSS(:)>XCTMIN(5) )            
               ZLBDAS(:) = ( XLBS*ZCSS(:) / ZRSS(:) )**XLBEXS    
               ZZY(:) = ZRHODREF(:)**(-XCEXVT) * (ZLBDAS(:)**(-XDS))
               ZZW(:) = XFSEDRS * ZRSS(:) * ZZY(:) * ZRHODREF(:)
               ZZX(:) = XFSEDCS * ZCSS(:) * ZZY(:) * ZRHODREF(:)
            END WHERE
            ZWSEDR(:,:,:) = UNPACK( ZZW(:),MASK=GSEDIM(:,:,:),FIELD=0.0 )
            ZWSEDR(:,:,IKB:IKE) = MIN( ZWSEDR(:,:,IKB:IKE), PRSS(:,:,IKB:IKE) * PRHODREF(:,:,IKB:IKE) / ZW(:,:,IKB:IKE) )
            ZWSEDC(:,:,:) = UNPACK( ZZX(:),MASK=GSEDIM(:,:,:),FIELD=0.0 )
            ZWSEDC(:,:,IKB:IKE) = MIN( ZWSEDC(:,:,IKB:IKE), PCSS(:,:,IKB:IKE) * PRHODREF(:,:,IKB:IKE) / ZW(:,:,IKB:IKE) )
            DO JK = IKB , IKE
               PRSS(:,:,JK) = PRSS(:,:,JK) + ZW(:,:,JK)*                      &
                    (ZWSEDR(:,:,JK+1)-ZWSEDR(:,:,JK))/PRHODREF(:,:,JK)
               PCSS(:,:,JK) = PCSS(:,:,JK) + ZW(:,:,JK)*                      &
                    (ZWSEDC(:,:,JK+1)-ZWSEDC(:,:,JK))/PRHODREF(:,:,JK)
            END DO
            DEALLOCATE(ZRSS)
            DEALLOCATE(ZCSS)
            DEALLOCATE(ZLBDAS) 
         ELSE
            DO JL = 1,ISEDIM
               ZRSS(JL) = PRSS(I1(JL),I2(JL),I3(JL))
            END DO
            WHERE( ZRSS(:)>XRTMIN(5) )
               ZZW(:) = XFSEDS * (ZRSS(:)*ZRHODREF(:))**XEXSEDS * ZRHODREF(:)**(-XCEXVT)
            END WHERE
            ZWSEDR(:,:,:) = UNPACK( ZZW(:),MASK=GSEDIM(:,:,:),FIELD=0.0 )
            ZWSEDR(:,:,IKB:IKE) = MIN( ZWSEDR(:,:,IKB:IKE), PRSS(:,:,IKB:IKE) * PRHODREF(:,:,IKB:IKE) / ZW(:,:,IKB:IKE) )
            DO JK = IKB , IKE
               PRSS(:,:,JK) = PRSS(:,:,JK) + ZW(:,:,JK)* &
                    (ZWSEDR(:,:,JK+1)-ZWSEDR(:,:,JK))/PRHODREF(:,:,JK)
            END DO
            DEALLOCATE(ZRSS)
         END IF
      ELSE
         ZWSEDR(:,:,IKB) = 0.0
      END IF
!    
      PINPRS(:,:) = PINPRS(:,:) + ZWSEDR(:,:,IKB)/XRHOLW/KSPLITG                          ! in m/s
!
!*       2.23   for graupeln
!
      ZZW(:) = 0.
      ZZX(:) = 0.0                                          
      ZZY(:) = 0.0
      IF( MAXVAL(PRGS(:,:,:))>XRTMIN(6) ) THEN
         ALLOCATE(ZRGS(ISEDIM)) 
         IF(NMOM_G.GE.2) THEN
                 ALLOCATE(ZCGS(ISEDIM))                                   
                 ALLOCATE(ZLBDAG(ISEDIM))                  
                 DO JL = 1,ISEDIM
                    ZRGS(JL) = PRGS(I1(JL),I2(JL),I3(JL))
                    ZCGS(JL) = PCGS(I1(JL),I2(JL),I3(JL))            
                 END DO 
                 ZLBDAG(:)  = 1.E10         
                 WHERE( ZRGS(:)>XRTMIN(6) .AND. ZCGS(:)>XCTMIN(6) )                          
                    ZLBDAG(:) = ( XLBS*ZCGS(:) / ZRGS(:) )**XLBEXG    
                    ZZY(:) = ZRHODREF(:)**(-XCEXVT) * (ZLBDAG(:)**(-XDG))
                    ZZW(:) = XFSEDRG * ZRGS(:) * ZZY(:) * ZRHODREF(:)
                    ZZX(:) = XFSEDCG * ZCGS(:) * ZZY(:) * ZRHODREF(:)
                 END WHERE
                 ZWSEDR(:,:,:) = UNPACK( ZZW(:),MASK=GSEDIM(:,:,:),FIELD=0.0 )
                 ZWSEDR(:,:,IKB:IKE) = MIN( ZWSEDR(:,:,IKB:IKE), PRGS(:,:,IKB:IKE) * PRHODREF(:,:,IKB:IKE) / ZW(:,:,IKB:IKE) )
                 ZWSEDC(:,:,:) = UNPACK( ZZX(:),MASK=GSEDIM(:,:,:),FIELD=0.0 )
                 ZWSEDC(:,:,IKB:IKE) = MIN( ZWSEDC(:,:,IKB:IKE), PCGS(:,:,IKB:IKE) * PRHODREF(:,:,IKB:IKE) / ZW(:,:,IKB:IKE) )
                 DO JK = IKB , IKE
                    PRGS(:,:,JK) = PRGS(:,:,JK) + ZW(:,:,JK)*                      &
                         (ZWSEDR(:,:,JK+1)-ZWSEDR(:,:,JK))/PRHODREF(:,:,JK)
                    PCGS(:,:,JK) = PCGS(:,:,JK) + ZW(:,:,JK)*                      &
                         (ZWSEDC(:,:,JK+1)-ZWSEDC(:,:,JK))/PRHODREF(:,:,JK)
                 END DO          
                 DEALLOCATE(ZRGS)
                 DEALLOCATE(ZCGS)           
                 DEALLOCATE(ZLBDAG)
         ELSE
                 DO JL = 1,ISEDIM
                    ZRGS(JL) = PRGS(I1(JL),I2(JL),I3(JL))
                 END DO
                 WHERE( ZRGS(:)>XRTMIN(6) )
                    ZZW(:) = XFSEDG * (ZRGS(:)*ZRHODREF(:))**XEXSEDG * ZRHODREF(:)**(-XCEXVT)
                 END WHERE
                 ZWSEDR(:,:,:) = UNPACK( ZZW(:),MASK=GSEDIM(:,:,:),FIELD=0.0 )
                 ZWSEDR(:,:,IKB:IKE) = MIN( ZWSEDR(:,:,IKB:IKE), PRGS(:,:,IKB:IKE) * PRHODREF(:,:,IKB:IKE) / ZW(:,:,IKB:IKE) )
                 DO JK = IKB , IKE
                    PRGS(:,:,JK) = PRGS(:,:,JK) + ZW(:,:,JK)* &
                         (ZWSEDR(:,:,JK+1)-ZWSEDR(:,:,JK))/PRHODREF(:,:,JK)
                 END DO
                 DEALLOCATE(ZRGS)
         END IF
      ELSE
         ZWSEDR(:,:,IKB) = 0.0
      END IF
!    
      PINPRG(:,:) = PINPRG(:,:) + ZWSEDR(:,:,IKB)/XRHOLW/KSPLITG                        ! in m/s
!
!*       2.23   for hail
!
      ZZW(:) = 0.
      ZZX(:) = 0.
      ZZY(:) = 0.      
      IF( MAXVAL(PRHS(:,:,:))>XRTMIN(7) ) THEN
         ALLOCATE(ZRHS(ISEDIM)) 
         IF(NMOM_H.GE.2) THEN
                ALLOCATE(ZCHS(ISEDIM))          
                ALLOCATE(ZLBDAH(ISEDIM))                  
                 DO JL = 1,ISEDIM
                    ZRHS(JL) = PRHS(I1(JL),I2(JL),I3(JL))
                    ZCHS(JL) = PCHS(I1(JL),I2(JL),I3(JL))            
                 END DO
                 ZLBDAH(:)  = 1.E10         
                 WHERE( ZRHS(:)>XRTMIN(7) .AND. ZCHS(:)>XCTMIN(7) )
                    ZLBDAH(:) = ( XLBH*ZCHS(:) / ZRHS(:) )**XLBEXH    
                    ZZY(:) = ZRHODREF(:)**(-XCEXVT) * (ZLBDAH(:)**(-XDH))
                    ZZW(:) = XFSEDRH * ZRHS(:) * ZZY(:) * ZRHODREF(:)
                    ZZX(:) = XFSEDCH * ZCHS(:) * ZZY(:) * ZRHODREF(:)
                 END WHERE
                 ZWSEDR(:,:,:) = UNPACK( ZZW(:),MASK=GSEDIM(:,:,:),FIELD=0.0 )
                 ZWSEDR(:,:,IKB:IKE) = MIN( ZWSEDR(:,:,IKB:IKE), PRHS(:,:,IKB:IKE) * PRHODREF(:,:,IKB:IKE) / ZW(:,:,IKB:IKE) )
                 DO JK = IKB , IKE
                    PRHS(:,:,JK) = PRHS(:,:,JK) + ZW(:,:,JK)*                      &
                         (ZWSEDR(:,:,JK+1)-ZWSEDR(:,:,JK))/PRHODREF(:,:,JK)
                 END DO
                 DEALLOCATE(ZRHS)
                 DEALLOCATE(ZLBDAH)
                 DEALLOCATE(ZCHS)
         ELSE
                 DO JL = 1,ISEDIM
                    ZRHS(JL) = PRHS(I1(JL),I2(JL),I3(JL))
                 END DO
                 WHERE( ZRHS(:)>XRTMIN(7) )
                    ZZW(:) = XFSEDH * (ZRHS(:)*ZRHODREF(:))**XEXSEDH * ZRHODREF(:)**(-XCEXVT)
                 END WHERE
                 ZWSEDR(:,:,:) = UNPACK( ZZW(:),MASK=GSEDIM(:,:,:),FIELD=0.0 )
                 ZWSEDR(:,:,IKB:IKE) = MIN( ZWSEDR(:,:,IKB:IKE), PRHS(:,:,IKB:IKE) * PRHODREF(:,:,IKB:IKE) / ZW(:,:,IKB:IKE) )
                 DO JK = IKB , IKE
                    PRHS(:,:,JK) = PRHS(:,:,JK) + ZW(:,:,JK)* &
                         (ZWSEDR(:,:,JK+1)-ZWSEDR(:,:,JK))/PRHODREF(:,:,JK)
                 END DO
                 DEALLOCATE(ZRHS)
         END IF
      ELSE
         ZWSEDR(:,:,IKB) = 0.0
      END IF
!    
      PINPRH(:,:) = PINPRH(:,:) + ZWSEDR(:,:,IKB)/XRHOLW/KSPLITG                        ! in m/s
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
         IF(NMOM_S.GE.2) PCSS(:,:,:) = PCSS(:,:,:) / PTSTEP
         IF(NMOM_G.GE.2) PCGS(:,:,:) = PCGS(:,:,:) / PTSTEP
         IF(NMOM_H.GE.2) PCHS(:,:,:) = PCHS(:,:,:) / PTSTEP 
      END IF
   END IF
END DO
!++cb++
DEALLOCATE(ZRTMIN)
DEALLOCATE(ZCTMIN)
!--cb--
!
END SUBROUTINE LIMA_COLD_SEDIMENTATION
!
!-------------------------------------------------------------------------------
