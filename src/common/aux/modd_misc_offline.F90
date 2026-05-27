MODULE MODD_MISC_OFFLINE
!
USE MODD_BUDGET,     ONLY: NBUDGET_RH, TBUDGETDATA_PTR, TBUDGETCONF_t
USE MODD_IO,         ONLY: TFILEDATA
IMPLICIT NONE
!
!> @file  
!!    MODD_MISC_OFFLINE - definition of a structure containing all the control parameters for the testprogs
!!
!!    This is a structure specifically built for the available testprogs.
!!    This would certainly be a bad idea to take this structure for an example on how to plug PHYEX
!!    in a real model.
!!
!!    The idea is to put here all the constants needed to call the parametrisations in order to reduce
!!    the number of objects to deal with in the calling loop.
TYPE :: MISC_OFFLINE_t
  REAL :: PTSTEP
  CHARACTER(LEN=4)         :: HBUNAME
  LOGICAL                  :: LMFCONV
  INTEGER                  :: KRR, KRRL, KRRI, KSV
  LOGICAL                  :: OCOMPUTE_SRC
  TYPE(TBUDGETDATA_PTR), DIMENSION(NBUDGET_RH) :: YLBUDGET
  INTEGER                  :: NBUDGET
  TYPE(TBUDGETCONF_t)      :: TBUCONF
  LOGICAL                  :: ONOMIXLG
  INTEGER                  :: KSV_LGBEG, KSV_LGEND
  REAL                     :: PDX, PDY
  INTEGER                  :: KGRADIENTSLEO, KGRADIENTSGOG, KHALO
  CHARACTER(LEN=4),DIMENSION(2)  :: HLBCX, HLBCY
  CHARACTER(LEN=6)         :: CPROGRAM
  INTEGER                  :: KSV_LIMA_NR, KSV_LIMA_NS, KSV_LIMA_NG, KSV_LIMA_NH
  LOGICAL                  :: O2D, OFLAT, OCOUPLES, OBLOWSNOW, OOCEAN, ODEEPOC
  LOGICAL                  :: OIBM, OFLYER
  TYPE(TFILEDATA)          :: ZTFILE
  REAL                     :: PRSNOW
  LOGICAL                  :: ODIAG_IN_RUN
  CHARACTER(LEN=4)         :: CMICRO
  LOGICAL                  :: OELEC=.FALSE.        !< Lightning prognostic scheme
  CHARACTER(LEN=4)         :: CELEC='NONE'         !< Name of the electricity scheme
  LOGICAL                  :: OSEDIM_BEARD=.FALSE. !< Switch for effect of electrical forces on sedim.
  TYPE(TFILEDATA)          :: TPFILE
  REAL,DIMENSION(8)        :: XINIRADIUS_SLT=(/0.009, 0.021, 0.045, 0.115,0.415, 2.5, 7.0, 25.0/)
  REAL,DIMENSION(8)        :: XINISIG_SLT =(/ 1.37, 1.5, 1.42, 1.53, 1.85, 1.55, 1.8, 2.1 /)
  LOGICAL                  :: LSALT=.FALSE.
  INTEGER                  :: NMODE_SLT=8
  CHARACTER(LEN=4)         :: CRGUNITS='NUMB'
  REAL                     :: XDENSITY_SALT=1.173e3
END TYPE MISC_OFFLINE_t
END MODULE MODD_MISC_OFFLINE
