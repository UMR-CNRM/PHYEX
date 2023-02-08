MODULE YOMPARAR

USE PARKIND1  ,ONLY : JPIM     ,JPRB
!
USE MODD_CST, ONLY: CST_t
USE MODD_PARAM_ICE, ONLY: PARAM_ICE_t
USE MODD_RAIN_ICE_DESCR, ONLY: RAIN_ICE_DESCR_t
USE MODD_RAIN_ICE_PARAM, ONLY: RAIN_ICE_PARAM_t
USE MODD_CLOUDPAR_N, ONLY: CLOUDPAR_t
USE MODD_PARAM_MFSHALL_N, ONLY: PARAM_MFSHALL_t

IMPLICIT NONE

SAVE

TYPE TPARAR
!*
!     ------------------------------------------------------------------

!     VARIABLES pour utiliser la PHYSIQUE de meso_NH :
!     VARIABLES to use the MESO-NH physics:

CHARACTER(LEN=4)   :: CMICRO      ! Microphysics scheme ('ICE3', 'ICE4' or 'LIMA')
TYPE(CST_t)            :: CST
TYPE(PARAM_ICE_t)      :: PARAM_ICE
TYPE(RAIN_ICE_DESCR_t) :: RAIN_ICE_DESCR
TYPE(RAIN_ICE_PARAM_t) :: RAIN_ICE_PARAM
TYPE(CLOUDPAR_t)       :: CLOUDPARN
TYPE(PARAM_MFSHALL_t)  :: PARAM_MFSHALLN
!
INTEGER(KIND=JPIM) :: NSPLITR     ! Time splitting for Eulerian sedimentation
INTEGER(KIND=JPIM) :: NSPLITG     ! Time splitting for Eulerian sedimentation
INTEGER(KIND=JPIM) :: NRR, NRRL, NRRI   !number of microphysical species
INTEGER(KIND=JPIM) :: NSV         !number of passiv variables in MesoNH,
                                  ! always 0 in AROME
INTEGER(KIND=JPIM) :: NSWB_MNH    !number of SW bands for surface
                                  ! (must be equal to NSW !!)
INTEGER(KIND=JPIM) :: NGPAR       !number of fields in the buffer containing
                                  ! the 2D pseudo-historical variables.
INTEGER(KIND=JPIM) :: MINPRR      !pointer on INPRR
INTEGER(KIND=JPIM) :: MINPRS      !pointer on INPRS
INTEGER(KIND=JPIM) :: MINPRG      !pointer on INPRG
INTEGER(KIND=JPIM) :: MACPRR      !pointer on ACPRR 
INTEGER(KIND=JPIM) :: MACPRS      !pointer on ACPRS
INTEGER(KIND=JPIM) :: MACPRG      !pointer on ACPRG
INTEGER(KIND=JPIM) :: MALBDIR     !pointer on ALBDIR
INTEGER(KIND=JPIM) :: MALBSCA     !pointer on ALBSCA
INTEGER(KIND=JPIM) :: MRAIN       !pointer on surface rain
INTEGER(KIND=JPIM) :: MSNOW       !pointer on surface snow
INTEGER(KIND=JPIM) :: MGZ0        !pointer on GZ0
INTEGER(KIND=JPIM) :: MGZ0H       !pointer on GZ0H
INTEGER(KIND=JPIM) :: MVQS        !pointer on surface moisture
INTEGER(KIND=JPIM) :: MVTS        !pointer on surface temperature
INTEGER(KIND=JPIM) :: MVEMIS      !pointer on surface emissivity
INTEGER(KIND=JPIM) :: MSWDIR      !pointer on SW direct surface flux
INTEGER(KIND=JPIM) :: MSWDIF      !pointer on SW surface diffuse flux
INTEGER(KIND=JPIM) :: MLSM        !pointer on land-sea mask
INTEGER(KIND=JPIM) :: MCD         !pointer on drag coefficient

REAL(KIND=JPRB), DIMENSION(:), ALLOCATABLE  :: XSW_BANDS  !SW spectral bands
! for ext. surface scheme
LOGICAL :: LOSUBG_COND ! switch to activate subgrid condensation
LOGICAL :: LOSIGMAS    ! activate calculation of variance of departure to 
                       ! saturation in turb scheme (to be used in subgrid condensation)
LOGICAL :: LOLSMC      ! Land/sea mask for cloud droplet number conc.
LOGICAL :: LOTOWNC     ! Town mask for cloud droplet number conc.

LOGICAL :: LKOGAN      ! Use Kogan autocoversion of liquid
LOGICAL :: LMODICEDEP  ! Logical switch for alternative dep/evap of ice
LOGICAL :: LICERAD     ! Assume higher fraction of condensate for 
                       ! ice/snow/graupel than the actual cloud cover in
                       ! radiation
REAL(KIND=JPRB) :: RADGR  ! Tuning of ice for radiation, TO BE REMOVED
REAL(KIND=JPRB) :: RADSN  ! Tuning of ice for radiation, TO BE REMOVED

REAL(KIND=JPRB) :: VSIGQSAT ! coeff applied to qsat variance contribution
                            ! for subgrid condensation

LOGICAL            :: LLCRIT    ! True if temperature dependent 
                                ! critical condensation in EDMFm 
LOGICAL            :: LTOTPREC  ! True if precipitation tendencies
                                ! from the sub-grid scheme are
                                ! added to the total precip tendencies.
LOGICAL            :: LTOTPRECL ! As LTOTPREC but updraft fraction untouched
LOGICAL            :: LHGT_QS   ! Switch for height dependent VQSIGSAT

! * for the squall line case:
LOGICAL :: LSQUALL ! use for the squall line case
INTEGER(KIND=JPIM) :: NREFROI1 !starting point for cooling (lsquall case) 
INTEGER(KIND=JPIM) :: NREFROI2 !end point for cooling (lsquall case)
REAL(KIND=JPRB) :: VSQUALL ! mean velocity displacement of the squall line (lsquall case)

! * for the MESO-NH physics printings:
INTEGER(KIND=JPIM) :: NPTP ! index in NPROMA paquet where the print will be done
INTEGER(KIND=JPIM) :: NPRINTFR !frequency of physical prints in apl_arome

!* for other diagnostics
!  wmax per vertical level
LOGICAL :: LDIAGWMAX !activate print of WMAX in apl_arome
INTEGER(KIND=JPIM) :: NDIAGWMAX ! frequency of preceding prints (in time step)

!* for chemical scheme
! time step factor
INTEGER(KIND=JPIM) :: NDTCHEM ! time step factor for chemical scheme
!* for MNH budget anlysis
LOGICAL :: LAROBU_ENABLE ! for MNH budget anlysis
!* for turbulence scheme
REAL(KIND=JPRB) :: XLINI ! minimum bl89 mixing length
LOGICAL :: LSTATNW !  updated full statistical cloud scheme
                   !  (yet only to be used in combination with EDMFm convection (DUAL))
LOGICAL :: LHARATU !  if true RACMO turbulence is used
                   !  (yet only to be used in combination with EDMFm convection (DUAL))
!* Subgrid precipitation scheme
CHARACTER (LEN=4) :: CSUBG_AUCV_RC ! type of rc->rr autoconversion scheme
                                   ! 'CLFR', 'PDF' or 'NONE'
CHARACTER(LEN=80) :: CSUBG_AUCV_RI ! type of ri->rs autoconversion scheme
                                   ! 'NONE', 'CLFR' or 'ADJU'
CHARACTER(LEN=80) :: CSUBG_MF_PDF  ! PDF to use on MF cloud to retrieve low and high cloud parts
                                   ! 'NONE' or 'TRIANGLE'

!
!* For total cumulative 3D prec flux for MOCAGE
LOGICAL :: LFPREC3D ! Switch on total cumulative 3D prec flux output (for MOCAGE use)
!* For radiation :
REAL(KIND=JPRB) :: XCQVR ! reduction factor of water vapour used for radiation computation.
REAL(KIND=JPRB) :: GQVPLIM ! pressure value over which qv is damped towards 0 for radiation.
REAL(KIND=JPRB) :: GQVTOP ! qv value at the top of the atmopshere.
LOGICAL :: LQVTOP ! to activate modification of qv in input to radiation.

CHARACTER(LEN=80) :: CCONDENS !condensation formulation. 'GAUS' or 'CB02'
CHARACTER(LEN=4) :: CLAMBDA3 !formulation for the lambda3 coeff used with s'r'. 'CB' or 'NONE'

END TYPE TPARAR

!!TYPE(TPARAR), POINTER :: YRPARAR => NULL()


!     ------------------------------------------------------------------
END MODULE YOMPARAR
