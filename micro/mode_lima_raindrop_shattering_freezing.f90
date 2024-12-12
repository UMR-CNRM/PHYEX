!MNH_LIC Copyright 2018-2025 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-------------------------------------------------------------------------------
MODULE MODE_LIMA_RAINDROP_SHATTERING_FREEZING
  IMPLICIT NONE
CONTAINS
!     #######################################################################
  SUBROUTINE LIMA_RAINDROP_SHATTERING_FREEZING (LDCOMPUTE,                    &
                                                PRHODREF,                     &
                                                PRRT, PCRT, PRIT, PCIT, PRGT, &
                                                PLBDR,                        &
                                                P_RI_RDSF, P_CI_RDSF, PT      )
!     #######################################################################
!
!!    PURPOSE
!!    -------
!!      Compute the raindrop shattering when freezing (secondary ice production process)
!!
!!    AUTHOR
!!    ------
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/2022  duplicate from original for LIMA_SPLIT
!!    I. Vongapseut 01/2024  Add a dependance on temperature
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST,              ONLY : XPI
USE MODD_PARAM_LIMA,       ONLY : XRTMIN, XCTMIN, LRDSF, XCEXVT, XPSH_MAX_RDSF
USE MODD_PARAM_LIMA_COLD,  ONLY : XMNU0
USE MODD_PARAM_LIMA_WARM,  ONLY : XDR
USE MODD_PARAM_LIMA_MIXED, ONLY : NGAMINC, XGAMINC_RDSF_R, &
                                  XFACTOR_RDSF_NI, XMOMGR_RDSF, XRDSFINTP1_R, XRDSFINTP_R, &
                                  XTM_PSH, XSIG_PSH
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
LOGICAL, DIMENSION(:),INTENT(IN)    :: LDCOMPUTE
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRHODREF
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRRT
REAL, DIMENSION(:),   INTENT(IN)    :: PCRT
REAL, DIMENSION(:),   INTENT(IN)    :: PRIT
REAL, DIMENSION(:),   INTENT(IN)    :: PCIT
REAL, DIMENSION(:),   INTENT(IN)    :: PRGT
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDR 
REAL, DIMENSION(:),  INTENT(IN)     :: PT              ! Temperature
!
REAL, DIMENSION(:),   INTENT(OUT)   :: P_RI_RDSF
REAL, DIMENSION(:),   INTENT(OUT)   :: P_CI_RDSF
!
!
!*       0.2   Declarations of local variables :
!
LOGICAL, DIMENSION(SIZE(PRRT))     :: GRDSF              ! Test where to compute collision process
INTEGER :: IRDSF
REAL,    DIMENSION(:), ALLOCATABLE :: ZVEC1_R            ! Work vectors for rain
REAL,    DIMENSION(:), ALLOCATABLE :: ZVEC1_R1           ! Work vectors for rain
REAL,    DIMENSION(:), ALLOCATABLE :: ZVEC2_R            ! Work vectors for rain
INTEGER, DIMENSION(:), ALLOCATABLE :: IVEC2_R            ! Rain indice vector
REAL,    DIMENSION(SIZE(PRRT))     :: ZINTG_RAIN         ! incomplete gamma function for rain
REAL,    DIMENSION(SIZE(PRRT))     :: ZNI_RDSF,ZRI_RDSF  ! RDSF rates
!
REAL,    DIMENSION(:), ALLOCATABLE :: ZT_RDSF            ! Temperature where RDSF
REAL,    DIMENSION(:), ALLOCATABLE :: ZNORMT             ! Normal distribution temperature shattering probability
REAL,    DIMENSION(:), ALLOCATABLE :: ZPSH_R             ! Shattering probability work vectors
REAL,    DIMENSION(SIZE(PRRT))     :: ZPSH               !Shattering probability
LOGICAL, DIMENSION(:), ALLOCATABLE :: LTEMP              ! Define the mask where ZNORMT is negligible
integer :: ii
!-------------------------------------------------------------------------------
P_RI_RDSF(:)=0.
P_CI_RDSF(:)=0.
!
GRDSF(:) = LRDSF .AND. LDCOMPUTE(:) .AND. (PRIT(:)>XRTMIN(4)) .AND. (PRRT(:)>XRTMIN(3)) &
                                    .AND. (PCIT(:)>XCTMIN(4)) .AND. (PCRT(:)>XCTMIN(3)) &
                                    .AND. (PRGT(:)>XRTMIN(6))

IRDSF    = COUNT( GRDSF(:) )
!
IF (IRDSF > 0) THEN
!
  ALLOCATE(ZVEC1_R(IRDSF))
  ALLOCATE(ZVEC1_R1(IRDSF))
  ALLOCATE(ZVEC2_R(IRDSF))
  ALLOCATE(IVEC2_R(IRDSF))
  ALLOCATE(ZT_RDSF(IRDSF))
  ALLOCATE(ZNORMT(IRDSF))
  ALLOCATE(ZPSH_R(IRDSF))
  ALLOCATE(LTEMP(IRDSF))
!
!
!*       2.2.1  select the ZLBDAR
!
  ZVEC1_R(:) = PACK( PLBDR(:),MASK=GRDSF(:) )
  ZT_RDSF(:) = PACK( PT(:),MASK=GRDSF(:) )
!
!
!*       2.2.2  find the next lower indice for the ZLBDAR in the
!               geometrical set of Lbda_r used to tabulate some moments of the
!               incomplete gamma function, for the lower boundary (0.1 mm)
!
  ZVEC2_R(1:IRDSF) = MAX( 1.00001, MIN( FLOAT(NGAMINC)-0.00001,XRDSFINTP_R  &
                        * LOG( ZVEC1_R(1:IRDSF) ) + XRDSFINTP1_R  ) )
  IVEC2_R(1:IRDSF) = INT( ZVEC2_R(1:IRDSF) )
  ZVEC2_R(1:IRDSF) = ZVEC2_R(1:IRDSF) - FLOAT( IVEC2_R(1:IRDSF) )
!
!*       2.2.3  perform the linear interpolation of the
!               normalized "2+XDR"-moment of the incomplete gamma function
!
  ZVEC1_R1(1:IRDSF) = XGAMINC_RDSF_R(IVEC2_R(1:IRDSF)+1) *  ZVEC2_R(1:IRDSF)    &
                    - XGAMINC_RDSF_R(IVEC2_R(1:IRDSF))   * (ZVEC2_R(1:IRDSF) - 1.0)
!
!  From 0.1 mm to infinity we need
  ZVEC1_R1(1:IRDSF) = XMOMGR_RDSF * (1.0 - ZVEC1_R1(1:IRDSF))
!
  ZINTG_RAIN(:) = UNPACK ( VECTOR=ZVEC1_R1(:),MASK=GRDSF,FIELD=0.0 )
!
!*       2.2.4  To compute final "RDSF" contributions
!
  LTEMP(1:IRDSF) = .FALSE.
!  Add shattering probability 
! 
  ZNORMT(1:IRDSF) =  1. / SQRT(2.*XSIG_PSH*XPI) * EXP(-(ZT_RDSF(1:IRDSF)-XTM_PSH)**2. / (2.*XSIG_PSH**2.))
  LTEMP(1:IRDSF) = ZNORMT(1:IRDSF) >= 0. .AND. ZNORMT(1:IRDSF) < 0.1
!
  ZPSH_R(1:IRDSF) = XPSH_MAX_RDSF / MAXVAL(ZNORMT(1:IRDSF)) * ZNORMT(1:IRDSF)
  !
  ZPSH(:) = UNPACK ( VECTOR=ZPSH_R(:),MASK=GRDSF,FIELD=0.0 )
!
  ZNI_RDSF(:) = (XFACTOR_RDSF_NI * ZPSH(:) / (PRHODREF(:)**(XCEXVT-1.0))) * (  &
                 PCIT(:) * PCRT(:) * ZINTG_RAIN(:) * PLBDR(:)**(-(XDR+6.0)) )
!
  P_CI_RDSF(:) = MAX(ZNI_RDSF(:),0.)
!
! The value of rg removed by RDSF is determined by the mean mass of pristine ice
!++cb++
!  ZRI_RDSF(:) = MAX( ZNI_RDSF(:)*XMNU0,XRTMIN(5) )
  ZRI_RDSF(:) = MAX( ZNI_RDSF(:)*XMNU0,0. )
!--cb--
!
  P_RI_RDSF(:) = ZRI_RDSF(:)
!
  DEALLOCATE(ZVEC1_R)
  DEALLOCATE(ZVEC1_R1)
  DEALLOCATE(ZVEC2_R)
  DEALLOCATE(IVEC2_R)
  DEALLOCATE(ZT_RDSF)
  DEALLOCATE(ZNORMT)
  DEALLOCATE(ZPSH_R)
  DEALLOCATE(LTEMP)
!
ENDIF
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_RAINDROP_SHATTERING_FREEZING
END MODULE MODE_LIMA_RAINDROP_SHATTERING_FREEZING
