MODULE MODD_FIELD
  USE MODD_PARAMETERS, ONLY: NGRIDUNKNOWN, NMNHNAMELGTMAX, NSTDNAMELGTMAX
  INTEGER, PARAMETER :: NMNHDIM_UNKNOWN             = -2
  INTEGER, PARAMETER :: NMNHMAXDIMS = 6 ! Cannot be less than 6
  INTEGER,PARAMETER :: TYPEUNDEF = -1, TYPEINT = 1, TYPELOG = 2, TYPEREAL = 3, TYPECHAR = 4, TYPEDATE = 5
!
TYPE TFIELDMETADATA
  CHARACTER(LEN=NMNHNAMELGTMAX) :: CMNHNAME  = '' !Name of the field (for MesoNH, non CF convention)
  CHARACTER(LEN=NSTDNAMELGTMAX) :: CSTDNAME  = '' !Standard name (CF convention)
  CHARACTER(LEN=32)  :: CLONGNAME = '' !Long name (CF convention)
  CHARACTER(LEN=40)  :: CUNITS    = '' !Canonical units (CF convention)
  CHARACTER(LEN=100) :: CCOMMENT  = '' !Comment (for MesoNH, non CF convention)
  INTEGER            :: NGRID     = NGRIDUNKNOWN !Localization on the model grid
  INTEGER            :: NTYPE     = TYPEUNDEF !Datatype
  INTEGER            :: NDIMS     = 0  !Number of dimensions
  INTEGER, DIMENSION(NMNHMAXDIMS) :: NDIMLIST = NMNHDIM_UNKNOWN ! List of dimensions of the data field
  !
  INTEGER            :: NFILLVALUE =  -2147483647            !Fill value for integer fields
  REAL               :: XFILLVALUE =  9.9692099683868690e+36 !Fill value for real fields
  INTEGER            :: NVALIDMIN  = -2147483646 !Minimum valid value for integer fields
  INTEGER            :: NVALIDMAX  =  2147483647 !Maximum valid value for integer fields
  REAL               :: XVALIDMIN  = -1.E36 !Minimum valid value for real fields
  REAL               :: XVALIDMAX  =  1.E36 !Maximum valid value for real fields
  CHARACTER(LEN=2)   :: CDIR      = '' !Type of the data field (XX,XY,--...)
  CHARACTER(LEN=4)   :: CLBTYPE   = 'NONE' !Type of the lateral boundary (LBX,LBY,LBXU,LBYV)
  LOGICAL            :: LTIMEDEP  = .FALSE. !Is the field time-dependent?
END TYPE TFIELDMETADATA
END MODULE MODD_FIELD
