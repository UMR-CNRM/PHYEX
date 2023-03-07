!     ##################
      MODULE DDH_MIX
!     ##################


!!    PURPOSE
!!    -------
!       module for type structure in  DDH

!!**  IMPLICIT ARGUMENTS
!!    ------------------

!!    AUTHOR
!!    ------
!!    O.Riviere   *Meteo France*
!!
!!    Modifications :
!!    --------------
!!    F.Voitus     18/09/17 : New DDH flexible structure for OpenMP
  
!-------------------------------------------------------------------------------


USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK,DR_HOOK
USE YOMLUN_IFSAUX, ONLY : NULOUT

IMPLICIT NONE

SAVE


!*       1.   TYPE DEFINITION


! structure for vertical profile
TYPE TYP_FIELD3D
  REAL(KIND=JPRB), POINTER  :: RVAL(:,:)=>NULL() ! values (domains x vertical)
  CHARACTER(LEN=16)         :: CNAME             ! name
  CHARACTER(LEN=1)          :: CTYPE             ! indicate wether it's a flux, variable or tendency
END TYPE TYP_FIELD3D

TYPE TYP_FIELD2D
  REAL(KIND=JPRB), POINTER :: RVAL(:)=>NULL() ! values (domains)
  CHARACTER(LEN=16)        :: CNAME           ! name
  CHARACTER(LEN=1)         :: CTYPE           ! indicate wether it's a flux, variable or tendency
END TYPE TYP_FIELD2D

TYPE TYP_BUVAR
  REAL(KIND=JPRB), POINTER  :: RVAL(:,:)=>NULL() ! values (domains)
  CHARACTER(LEN=2)          :: CNAME             ! name
END TYPE TYP_BUVAR


! structure containing fields, weights for horizontal averaging, etc.
TYPE TYP_DDH
  INTEGER(KIND=JPIM)          :: NLEV                  ! vertical dimension
  INTEGER(KIND=JPIM)          :: NPROMA                ! horizontal dimension                                                    
  INTEGER(KIND=JPIM)          :: KST                   ! first point
  INTEGER(KIND=JPIM)          :: KEND                  ! last point

  ! 3D field info
  INTEGER(KIND=JPIM)          :: NFIELDS3D             ! number of fields
  INTEGER(KIND=JPIM)          :: NFIELDS3D_OFFSET      ! offset value for the fields when 
                                                       ! putting in the global arrays
  INTEGER(KIND=JPIM)          :: NFIELDS3D_AUTO        ! number of fields that should be 
                                                       ! automatically allocated 
                                                       ! (maximum of cpg and cpglag)
  TYPE(TYP_FIELD3D), POINTER  :: YFIELD3D(:)=>NULL()   ! array of fields
  REAL(KIND=JPRB), POINTER    :: RVAL3D(:,:,:)=>NULL() ! auxiliary array for fast 
                                                       ! (automatic) allocation
                                                       ! (vertical x fields x domains)
  ! 2D field info
  INTEGER(KIND=JPIM)          :: NFIELDS2D             ! number of fields
  INTEGER(KIND=JPIM)          :: NFIELDS2D_OFFSET      ! offset value for the fields when 
                                                       ! putting in the global arrays
  INTEGER(KIND=JPIM)          :: NFIELDS2D_AUTO        ! number of fields that should be 
                                                       ! automatically allocated 
                                                       ! (maximum of cpg and cpglag)
  TYPE(TYP_FIELD2D), POINTER  :: YFIELD2D(:)=>NULL()   ! array of fields
  REAL(KIND=JPRB), POINTER    :: RVAL2D(:,:)=>NULL()   ! auxiliary array for fast 
                                                       ! (automatic) allocation
                                                       ! (vertical x fields x domains)
  ! horizontal info
  REAL(KIND=JPRB), POINTER    :: WEIGHT(:)=>NULL()     ! weights inside one NPROMA block
  INTEGER(KIND=JPIM), POINTER :: NDDHI(:)=>NULL()      ! cfr. KDDHI in cpg_dia


  TYPE(TYP_BUVAR), POINTER    :: YVARMULT(:)=>NULL()    ! array of fields
  REAL(KIND=JPRB), POINTER    :: RVARSM(:,:,:,:)=>NULL()! auxiliary array for fast 
                                                        ! (automatic) allocation

END TYPE TYP_DDH



END MODULE DDH_MIX

