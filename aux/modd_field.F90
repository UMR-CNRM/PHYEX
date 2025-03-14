MODULE MODD_FIELD
  USE MODD_PARAMETERS, ONLY: NGRIDUNKNOWN, NMNHNAMELGTMAX, NSTDNAMELGTMAX
  IMPLICIT NONE
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
  !
  CHARACTER(LEN=2)   :: CDIR      = '' !Type of the data field (XX,XY,--...)
  LOGICAL            :: LTIMEDEP  = .FALSE. !Is the field time-dependent?
END TYPE TFIELDMETADATA
END MODULE MODD_FIELD
