!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 conv 2006/05/18 13:07:25
!-----------------------------------------------------------------
!     #################
      MODULE MODI_CONVECT_CHEM_TRANSPORT
!     #################
!
INTERFACE
!
       SUBROUTINE CONVECT_CHEM_TRANSPORT( KLON, KLEV, KCH, PCH1, PCH1C,       &
                                         KDPL, KPBL, KLCL, KCTL, KLFS, KDBL, &
                                         PUMF, PUER, PUDR, PDMF, PDER, PDDR, &
                                         PTIMEC, PDXDY, PMIXF, PLMASS, PWSUB,&
                                         KFTSTEPS )

!
INTEGER,                INTENT(IN) :: KLON     ! horizontal dimension
INTEGER,                INTENT(IN) :: KLEV     ! vertical dimension
INTEGER,                INTENT(IN) :: KCH      ! number of passive tracers
!
REAL,DIMENSION(KLON,KLEV,KCH),INTENT(IN) :: PCH1 ! grid scale tracer concentr.
REAL,DIMENSION(KLON,KLEV,KCH),INTENT(OUT):: PCH1C! conv adjusted tracer concntr.
!
INTEGER, DIMENSION(KLON), INTENT(IN) :: KDPL   ! index for departure level
INTEGER, DIMENSION(KLON), INTENT(IN) :: KPBL   ! index for top of source layer
INTEGER, DIMENSION(KLON), INTENT(IN) :: KLCL   ! index lifting condens. level
INTEGER, DIMENSION(KLON), INTENT(IN) :: KCTL   ! index for cloud top level
INTEGER, DIMENSION(KLON), INTENT(IN) :: KLFS   ! index for level of free sink
INTEGER, DIMENSION(KLON), INTENT(IN) :: KDBL   ! index for downdraft base level
!
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PUMF ! updraft mass flux (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PUER ! updraft entrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PUDR ! updraft detrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PDMF ! downdraft mass flux (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PDER ! downdraft entrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PDDR ! downdraft detrainment (kg/s)
!
REAL, DIMENSION(KLON),     INTENT(IN) :: PTIMEC! convection time step
REAL, DIMENSION(KLON),     INTENT(IN) :: PDXDY ! grid area (m^2)
REAL, DIMENSION(KLON),     INTENT(IN) :: PMIXF ! mixed fraction at LFS
REAL, DIMENSION(KLON,KLEV),INTENT(IN) :: PLMASS! mass of model layer (kg)
REAL, DIMENSION(KLON,KLEV),INTENT(IN) :: PWSUB ! envir. compensating subsidence(Pa/s)
INTEGER,                INTENT(IN) :: KFTSTEPS  ! maximum fractional time steps
!
END SUBROUTINE CONVECT_CHEM_TRANSPORT
!
END INTERFACE
!
END MODULE MODI_CONVECT_CHEM_TRANSPORT
!     ########################################################################
      SUBROUTINE CONVECT_CHEM_TRANSPORT( KLON, KLEV, KCH, PCH1, PCH1C,       &
                                         KDPL, KPBL, KLCL, KCTL, KLFS, KDBL, &
                                         PUMF, PUER, PUDR, PDMF, PDER, PDDR, &
                                         PTIMEC, PDXDY, PMIXF, PLMASS, PWSUB,&
                                         KFTSTEPS )
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
USE MODD_CST
USE MODD_CONVPAREXT
USE MODD_NSV,  ONLY : NSV_LGBEG,NSV_LGEND
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER,                INTENT(IN) :: KLON     ! horizontal dimension
INTEGER,                INTENT(IN) :: KLEV     ! vertical dimension
INTEGER,                INTENT(IN) :: KCH      ! number of passive tracers
!
REAL,DIMENSION(KLON,KLEV,KCH),INTENT(IN) :: PCH1 ! grid scale tracer concentr.
REAL,DIMENSION(KLON,KLEV,KCH),INTENT(OUT):: PCH1C! conv adjusted tracer concntr.
!
INTEGER, DIMENSION(KLON), INTENT(IN) :: KDPL   ! index for departure level
INTEGER, DIMENSION(KLON), INTENT(IN) :: KPBL   ! index for top of source layer
INTEGER, DIMENSION(KLON), INTENT(IN) :: KLCL   ! index lifting condens. level
INTEGER, DIMENSION(KLON), INTENT(IN) :: KCTL   ! index for cloud top level
INTEGER, DIMENSION(KLON), INTENT(IN) :: KLFS   ! index for level of free sink
INTEGER, DIMENSION(KLON), INTENT(IN) :: KDBL   ! index for downdraft base level
!
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PUMF ! updraft mass flux (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PUER ! updraft entrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PUDR ! updraft detrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PDMF ! downdraft mass flux (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PDER ! downdraft entrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PDDR ! downdraft detrainment (kg/s)
!
REAL, DIMENSION(KLON),     INTENT(IN) :: PTIMEC! convection time step
REAL, DIMENSION(KLON),     INTENT(IN) :: PDXDY ! grid area (m^2)
REAL, DIMENSION(KLON),     INTENT(IN) :: PMIXF ! mixed fraction at LFS
REAL, DIMENSION(KLON,KLEV),INTENT(IN) :: PLMASS! mass of model layer (kg)
REAL, DIMENSION(KLON,KLEV),INTENT(IN) :: PWSUB ! envir. compensating subsidence(Pa/s)
INTEGER,                INTENT(IN) :: KFTSTEPS  ! maximum fractional time steps
!
!
!*       0.2   Declarations of local variables :
!
INTEGER :: INCH1          ! number of chemical tracers
INTEGER :: IIE, IKB, IKE  ! horizontal + vertical loop bounds
INTEGER :: IKS            ! vertical dimension
INTEGER :: JI             ! horizontal loop index
INTEGER :: JK, JKP        ! vertical loop index
INTEGER :: JN             ! chemical tracer loop index
INTEGER :: JSTEP          ! fractional time loop index
INTEGER :: JKLD, JKLP, JKMAX ! loop index for levels
!
REAL, DIMENSION(KLON,KLEV)     :: ZOMG ! compensat. subsidence (Pa/s)
REAL, DIMENSION(KLON,KLEV,KCH) :: ZUCH1, ZDCH1 ! updraft/downdraft values
REAL, DIMENSION(KLON)          :: ZTIMEC  ! fractional convective time step
REAL, DIMENSION(KLON,KLEV)     :: ZTIMC! 2D work array for ZTIMEC
REAL, DIMENSION(KLON,KLEV,KCH) :: ZCH1MFIN, ZCH1MFOUT
                                   ! work arrays for environm. compensat. mass
REAL, DIMENSION(KLON,KCH)      :: ZWORK1, ZWORK2, ZWORK3
!
!-------------------------------------------------------------------------------
!
!*       0.3   Compute loop bounds
!              -------------------
!
INCH1  = KCH
IIE    = KLON
IKB    = 1 + JCVEXB
IKS    = KLEV
IKE    = KLEV - JCVEXT
JKMAX  = MAXVAL( KCTL(:) )
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
DO JI = 1, IIE
    JKLD = KDPL(JI)
    JKLP = KPBL(JI)
    ZWORK1(JI,:) = .5 * ( PCH1(JI,JKLD,:) + PCH1(JI,JKLP,:) )
END DO
!
!*      2.2     Final updraft loop
!               ------------------
!
DO JK = MINVAL( KDPL(:) ), JKMAX
JKP = JK + 1
!
    DO JN = 1, INCH1
     DO JI = 1, IIE
       IF ( KDPL(JI) <= JK .AND. KLCL(JI) > JK )                             &
            ZUCH1(JI,JK,JN) = ZWORK1(JI,JN)
!
       IF ( KLCL(JI) - 1 <= JK .AND. KCTL(JI) > JK ) THEN
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
ZWORK1(:,:) = SPREAD( PMIXF(:), DIM=2, NCOPIES=INCH1 )
DO JI = 1, IIE
     JK = KLFS(JI)
     ZDCH1(JI,JK,:) = ZWORK1(JI,:) * PCH1(JI,JK,:) +                          &
                                       ( 1. - ZWORK1(JI,:) ) * ZUCH1(JI,JK,:)
END DO
!
!*      3.2     Final downdraft loop
!               --------------------
!
DO JK = MAXVAL( KLFS(:) ), IKB + 1, -1
JKP = JK - 1
    DO JN = 1, INCH1
    DO JI = 1, IIE
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
PCH1C(:,IKB:IKE,:) = PCH1(:,IKB:IKE,:) ! initialize adjusted envir. values
!
DO JK = IKB, IKE
   ZOMG(:,JK) = PWSUB(:,JK) * PDXDY(:) / XG ! environmental subsidence
END DO
!
ZTIMEC(:) = PTIMEC(:) / REAL( KFTSTEPS ) ! adjust  fractional time step
                                         ! to be an integer multiple of PTIMEC
WHERE ( PTIMEC(:) < 1. ) ZTIMEC(:) = 0.
ZTIMC(:,:)= SPREAD( ZTIMEC(:), DIM=2, NCOPIES=IKS )
!
ZCH1MFIN(:,:,:)   = 0.
ZCH1MFOUT(:,:,:)  = 0.
!
DO JSTEP = 1, KFTSTEPS ! Enter the fractional time step loop
!
      DO JK = IKB + 1, JKMAX
          JKP = MAX( IKB + 1, JK - 1 )
          ZWORK3(:,:) = SPREAD( ZOMG(:,JK), DIM=2, NCOPIES=INCH1 )
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
         PCH1C(:,JK,JN) = PCH1C(:,JK,JN) + ZTIMC(:,JK) / PLMASS(:,JK) *  (    &
                      ZCH1MFIN(:,JK,JN) + PUDR(:,JK) * ZUCH1(:,JK,JN) +       &
                      PDDR(:,JK) * ZDCH1(:,JK,JN) - ZCH1MFOUT(:,JK,JN) -      &
                      ( PUER(:,JK) + PDER(:,JK) ) * PCH1(:,JK,JN)    )
         IF(JN < NSV_LGBEG .OR. JN>NSV_LGEND-1) THEN
           PCH1C(:,JK,JN) = MAX( 0., PCH1C(:,JK,JN) )
         ELSE
           ! no tendency for horizontal Lagrangian variables
           PCH1C(:,JK,JN) = PCH1(:,JK,JN)
         END IF
       END DO
      END DO
!
END DO ! Exit the fractional time step loop
!
!
END SUBROUTINE CONVECT_CHEM_TRANSPORT
