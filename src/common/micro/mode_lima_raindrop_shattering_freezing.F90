!MNH_LIC Copyright 2018-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-------------------------------------------------------------------------------
MODULE MODE_LIMA_RAINDROP_SHATTERING_FREEZING
  IMPLICIT NONE
CONTAINS
!     #######################################################################
  SUBROUTINE LIMA_RAINDROP_SHATTERING_FREEZING (LIMAP, LIMAW, LIMAC, LIMAM, KSIZE, ODCOMPUTE,             &
                                                PRHODREF,                     &
                                                PRRT, PCRT, PRIT, PCIT, PRGT, &
                                                PLBDR,                        &
                                                P_RI_RDSF, P_CI_RDSF          )
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
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAM_LIMA_MIXED, ONLY:PARAM_LIMA_MIXED_T
USE MODD_PARAM_LIMA_COLD, ONLY:PARAM_LIMA_COLD_T
USE MODD_PARAM_LIMA_WARM, ONLY:PARAM_LIMA_WARM_T
USE MODD_PARAM_LIMA, ONLY:PARAM_LIMA_T
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
INTEGER, INTENT(IN) :: KSIZE
LOGICAL, DIMENSION(KSIZE),INTENT(IN)    :: ODCOMPUTE
!
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRHODREF
!
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRRT
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PCRT
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRIT
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PCIT
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRGT
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLBDR 
!
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_RI_RDSF
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_CI_RDSF
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
TYPE(PARAM_LIMA_MIXED_T),INTENT(IN)::LIMAM
TYPE(PARAM_LIMA_COLD_T),INTENT(IN)::LIMAC
TYPE(PARAM_LIMA_WARM_T),INTENT(IN)::LIMAW
TYPE(PARAM_LIMA_T),INTENT(IN)::LIMAP
REAL,    DIMENSION(SIZE(PRRT))     :: ZNI_RDSF,ZRI_RDSF  ! RDSF rates
!
!-------------------------------------------------------------------------------
P_RI_RDSF(:)=0.
P_CI_RDSF(:)=0.
!
GRDSF(:) = LIMAP%LRDSF .AND. ODCOMPUTE .AND. (PRIT(:)>LIMAP%XRTMIN(4)) .AND. (PRRT(:)>LIMAP%XRTMIN(3)) &
                                 .AND. (PCIT(:)>LIMAP%XCTMIN(4)) .AND. (PCRT(:)>LIMAP%XCTMIN(3)) &
                                 .AND. (PRGT(:)>LIMAP%XRTMIN(6))

IRDSF    = COUNT( GRDSF(:) )
!
IF (IRDSF > 0) THEN
!
  ALLOCATE(ZVEC1_R(IRDSF))
  ALLOCATE(ZVEC1_R1(IRDSF))
  ALLOCATE(ZVEC2_R(IRDSF))
  ALLOCATE(IVEC2_R(IRDSF))
!
!*       2.2.1  select the ZLBDAR
!
  ZVEC1_R(:) = PACK( PLBDR(:),MASK=GRDSF(:) )
!
!*       2.2.2  find the next lower indice for the ZLBDAR in the
!               geometrical set of Lbda_r used to tabulate some moments of the
!               incomplete gamma function, for the lower boundary (0.1 mm)
!
  ZVEC2_R(1:IRDSF) = MAX( 1.00001, MIN( FLOAT(LIMAM%NGAMINC)-0.00001,LIMAM%XRDSFINTP_R  &
                        * LOG( ZVEC1_R(1:IRDSF) ) + LIMAM%XRDSFINTP1_R  ) )
  IVEC2_R(1:IRDSF) = INT( ZVEC2_R(1:IRDSF) )
  ZVEC2_R(1:IRDSF) = ZVEC2_R(1:IRDSF) - FLOAT( IVEC2_R(1:IRDSF) )
!
!*       2.2.3  perform the linear interpolation of the
!               normalized "2+XDR"-moment of the incomplete gamma function
!
  ZVEC1_R1(1:IRDSF) = LIMAM%XGAMINC_RDSF_R(IVEC2_R(1:IRDSF)+1) *  ZVEC2_R(1:IRDSF)    &
                    - LIMAM%XGAMINC_RDSF_R(IVEC2_R(1:IRDSF))   * (ZVEC2_R(1:IRDSF) - 1.0)
!
!  From 0.1 mm to infinity we need
  ZVEC1_R1(1:IRDSF) = LIMAM%XMOMGR_RDSF * (1.0 - ZVEC1_R1(1:IRDSF))
!
  ZINTG_RAIN(:) = UNPACK ( VECTOR=ZVEC1_R1(:),MASK=GRDSF,FIELD=0.0 )
!
!*       2.2.4  To compute final "RDSF" contributions
!
  ZNI_RDSF(:) = (LIMAM%XFACTOR_RDSF_NI / (PRHODREF(:)**(LIMAP%XCEXVT-1.0))) * (  &
                 PCIT(:) * PCRT(:) * ZINTG_RAIN(:) * PLBDR(:)**(-(LIMAW%XDR+6.0)) )
!
  P_CI_RDSF(:) = MAX(ZNI_RDSF(:),0.)
!
! The value of rg removed by RDSF is determined by the mean mass of pristine ice
  ZRI_RDSF(:) = MAX( ZNI_RDSF(:)*LIMAC%XMNU0,LIMAP%XRTMIN(5) )
!
  P_RI_RDSF(:) = ZRI_RDSF(:)
!
  DEALLOCATE(ZVEC1_R)
  DEALLOCATE(ZVEC1_R1)
  DEALLOCATE(ZVEC2_R)
  DEALLOCATE(IVEC2_R)
  !
ENDIF
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_RAINDROP_SHATTERING_FREEZING
END MODULE MODE_LIMA_RAINDROP_SHATTERING_FREEZING
