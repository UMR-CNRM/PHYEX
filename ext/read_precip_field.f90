!MNH_LIC Copyright 1996-2020 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     #############################
      MODULE MODI_READ_PRECIP_FIELD
!     #############################
!
!
!
INTERFACE
!
      SUBROUTINE READ_PRECIP_FIELD(TPINIFILE,HPROGRAM,HCONF,                       &
                              HGETRCT,HGETRRT,HGETRST,HGETRGT,HGETRHT,             &
                              PINPRC,PACPRC,PINDEP,PACDEP,PINPRR,PINPRR3D,PEVAP3D, &
                              PACPRR,PINPRS,PACPRS,PINPRG,PACPRG,PINPRH,PACPRH     )
!
USE MODD_IO, ONLY : TFILEDATA
!
!*       0.1   declarations of arguments
!
TYPE(TFILEDATA),        INTENT(IN)    :: TPINIFILE ! Initial file
CHARACTER (LEN=*),      INTENT(IN)    :: HPROGRAM  !
CHARACTER (LEN=*),      INTENT(IN)    :: HCONF     !
!                    
CHARACTER (LEN=*),      INTENT(IN)    :: HGETRCT, HGETRRT, HGETRST, HGETRGT, HGETRHT
                                                  ! Get indicator RCT,RRT,RST,RGT,RHT
!
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PINPRC   ! Droplet instant precip
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PACPRC   ! Droplet accumulated precip
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PINDEP   ! Droplet instant deposition
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PACDEP   ! Droplet accumulated dep
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PINPRR   ! Rain instant precip
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PINPRR3D ! Rain precipitation flux 3D
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PEVAP3D  ! Rain evaporation flux 3D
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PACPRR   ! Rain accumulated precip
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PINPRS   ! Snow instant precip
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PACPRS   ! Snow accumulated precip
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PINPRG   ! Graupel instant precip
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PACPRG   ! Graupel accumulated precip
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PINPRH   ! Hail instant precip
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PACPRH   ! Hail accumulated precip
!
END SUBROUTINE READ_PRECIP_FIELD
!
END INTERFACE
!
END MODULE MODI_READ_PRECIP_FIELD 
!
!     ##############################################################################
      SUBROUTINE READ_PRECIP_FIELD(TPINIFILE,HPROGRAM,HCONF,                       &
                              HGETRCT,HGETRRT,HGETRST,HGETRGT,HGETRHT,             &
                              PINPRC,PACPRC,PINDEP,PACDEP,PINPRR,PINPRR3D,PEVAP3D, &
                              PACPRR,PINPRS,PACPRS,PINPRG,PACPRG,PINPRH,PACPRH     )
!     ##############################################################################
!
!!****  *READ_PRECIP_FIELD* - routine to read precipitation surface fields
!!
!!    PURPOSE
!!    -------
!       Initialize precipitation fields by reading their value in an initial
!     MNH file.
!
!!**  METHOD
!!    ------
!!    
!!    
!!
!!    EXTERNAL
!!    --------
!!      FMREAD   : to read data in LFIFM file
!!       
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!      None
!!
!!    REFERENCE
!!    ---------
!!      Book2 of the documentation (routine READ_PRECIP_FIELD)
!!
!!    AUTHOR
!!    ------
!!  	J.-P. Pinty     *Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       13/06/96 
!!      (J. Viviand)   04/02/97  convert precipitation rates in m/s
!!      (V. Ducrocq)   14/08/98  // remove KIINF,KJINF,KISUP,KJSUP
!!      (JP Pinty)     29/11/02  add C3R5, ICE2, ICE4
!!      (C.Lac)        04/03/13  add YGETxxx for FIT scheme
!!                    10/2016 (C.Lac) Add droplet deposition
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!  P. Wautelet 14/02/2019: remove CLUOUT/CLUOUT0 and associated variables
!!
!-----------------------------------------------------------------------------
!
!*       0.    DECLARATIONS

use modd_field,         only: tfieldmetadata, tfieldlist
USE MODD_IO,            ONLY: TFILEDATA
USE MODD_PARAM_ICE_n,     ONLY: LDEPOSC
USE MODD_PARAM_C2R2,    ONLY: LDEPOC
USE MODD_PARAM_LIMA,    ONLY: MDEPOC=>LDEPOC
!
use mode_field,         only: Find_field_id_from_mnhname
USE MODE_IO_FIELD_READ, only: IO_Field_read
!
IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
TYPE(TFILEDATA),        INTENT(IN)    :: TPINIFILE ! Initial file
CHARACTER (LEN=*),      INTENT(IN)    :: HPROGRAM  !
CHARACTER (LEN=*),      INTENT(IN)    :: HCONF     !
!                    
CHARACTER (LEN=*),      INTENT(IN)    :: HGETRCT, HGETRRT, HGETRST, HGETRGT, HGETRHT
                                                  ! Get indicator RCT,RRT,RST,RGT,RHT
!
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PINPRC   ! Droplet instant precip
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PACPRC   ! Droplet accumulated precip
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PINDEP   ! Droplet instant deposition
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PACDEP   ! Droplet accumulated dep
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PINPRR   ! Rain instant precip
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PINPRR3D ! Rain precipitation flux 3D
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PEVAP3D  ! Rain evaporation flux 3D
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PACPRR   ! Rain accumulated precip
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PINPRS   ! Snow instant precip
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PACPRS   ! Snow accumulated precip
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PINPRG   ! Graupel instant precip
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PACPRG   ! Graupel accumulated precip
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PINPRH   ! Hail instant precip
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PACPRH   ! Hail accumulated precip
!
!*       0.2   declarations of local variables
!
REAL, DIMENSION(SIZE(PINPRR,1),SIZE(PINPRR,2)) :: Z2D ! 2D array to read  data
REAL, DIMENSION(SIZE(PINPRR3D,1),SIZE(PINPRR3D,2),SIZE(PINPRR3D,3)) :: Z3D ! 3D array to read  data
                                                  ! in initial file 
INTEGER              :: IID
INTEGER              :: IRESP
CHARACTER(LEN=4)     :: YGETRCT,YGETRRT,YGETRST,YGETRGT,YGETRHT
TYPE(TFIELDMETADATA) :: TZFIELD
!
!-------------------------------------------------------------------------------
!
!*       1..    INITIALIZATION
!              ----------------
!
IF ((HPROGRAM == 'MESONH') .AND. (HCONF == 'START')) THEN
  YGETRCT = 'INIT'
  YGETRRT = 'INIT'
  YGETRST = 'INIT'
  YGETRGT = 'INIT'
  YGETRHT = 'INIT'
ELSE
  YGETRCT = HGETRCT
  YGETRRT = HGETRRT
  YGETRST = HGETRST
  YGETRGT = HGETRGT
  YGETRHT = HGETRHT
END IF
!-------------------------------------------------------------------------------
!
!*       2..    READ PROGNOSTIC VARIABLES
!              -------------------------
!
IF (SIZE(PINPRC) /= 0 ) THEN
  SELECT CASE(YGETRCT)
  CASE ('READ')
    CALL FIND_FIELD_ID_FROM_MNHNAME('INPRC',IID,IRESP)
    TZFIELD = TFIELDMETADATA( TFIELDLIST(IID) )
    TZFIELD%CUNITS = 'mm hour-1'
    CALL IO_Field_read(TPINIFILE,TZFIELD,Z2D,IRESP)
    IF (IRESP == 0) PINPRC(:,:)=Z2D(:,:)/(1000.*3600.)
    !
    CALL FIND_FIELD_ID_FROM_MNHNAME('ACPRC',IID,IRESP)
    TZFIELD = TFIELDMETADATA( TFIELDLIST(IID) )
    TZFIELD%CUNITS = 'mm'
    CALL IO_Field_read(TPINIFILE,TZFIELD,Z2D,IRESP)
    IF (IRESP == 0) PACPRC(:,:)=Z2D(:,:)/(1000.)
  CASE ('INIT')
    PINPRC(:,:) = 0.0
    PACPRC(:,:) = 0.0
  END SELECT
END IF
!
IF (SIZE(PINDEP) /= 0 ) THEN
  SELECT CASE(YGETRCT)
  CASE ('READ')
    CALL FIND_FIELD_ID_FROM_MNHNAME('INDEP',IID,IRESP)
    TZFIELD = TFIELDMETADATA( TFIELDLIST(IID) )
    TZFIELD%CUNITS = 'mm hour-1'
    CALL IO_Field_read(TPINIFILE,TZFIELD,Z2D,IRESP)
    IF (IRESP == 0) PINDEP(:,:)=Z2D(:,:)/(1000.*3600.)
    !
    CALL FIND_FIELD_ID_FROM_MNHNAME('ACDEP',IID,IRESP)
    TZFIELD = TFIELDMETADATA( TFIELDLIST(IID) )
    TZFIELD%CUNITS = 'mm'
    CALL IO_Field_read(TPINIFILE,TZFIELD,Z2D,IRESP)
    IF (IRESP == 0) PACDEP(:,:)=Z2D(:,:)/(1000.)
  CASE ('INIT')
    PINDEP(:,:) = 0.0
    PACDEP(:,:) = 0.0
  END SELECT
END IF
!
IF (SIZE(PINPRR) /= 0 ) THEN
  SELECT CASE(YGETRRT)
  CASE ('READ')
    CALL FIND_FIELD_ID_FROM_MNHNAME('INPRR',IID,IRESP)
    TZFIELD = TFIELDMETADATA( TFIELDLIST(IID) )
    TZFIELD%CUNITS = 'mm hour-1'
    CALL IO_Field_read(TPINIFILE,TZFIELD,Z2D,IRESP)
    IF (IRESP == 0) PINPRR(:,:)=Z2D(:,:)/(1000.*3600.)
    !
    CALL IO_Field_read(TPINIFILE,'INPRR3D',Z3D,IRESP)
    IF (IRESP == 0) PINPRR3D(:,:,:)=Z3D(:,:,:)
    !
    CALL IO_Field_read(TPINIFILE,'EVAP3D',Z3D,IRESP)
    IF (IRESP == 0) PEVAP3D(:,:,:)=Z3D(:,:,:)
    !
    CALL FIND_FIELD_ID_FROM_MNHNAME('ACPRR',IID,IRESP)
    TZFIELD = TFIELDMETADATA( TFIELDLIST(IID) )
    TZFIELD%CUNITS = 'mm'
    CALL IO_Field_read(TPINIFILE,TZFIELD,Z2D,IRESP)
    IF (IRESP == 0) PACPRR(:,:)=Z2D(:,:)/(1000.)
  CASE ('INIT')
    PINPRR(:,:) = 0.0
    PINPRR3D(:,:,:) = 0.0
    PEVAP3D(:,:,:) = 0.0
    PACPRR(:,:) = 0.0
  END SELECT
END IF
!
IF (SIZE(PINPRS) /= 0 ) THEN
  SELECT CASE(YGETRST)
  CASE ('READ')
    CALL FIND_FIELD_ID_FROM_MNHNAME('INPRS',IID,IRESP)
    TZFIELD = TFIELDMETADATA( TFIELDLIST(IID) )
    TZFIELD%CUNITS = 'mm hour-1'
    CALL IO_Field_read(TPINIFILE,TZFIELD,Z2D,IRESP)
    IF (IRESP == 0) PINPRS(:,:)=Z2D(:,:)/(1000.*3600.)
    !
    CALL FIND_FIELD_ID_FROM_MNHNAME('ACPRS',IID,IRESP)
    TZFIELD = TFIELDMETADATA( TFIELDLIST(IID) )
    TZFIELD%CUNITS = 'mm'
    CALL IO_Field_read(TPINIFILE,TZFIELD,Z2D,IRESP)
    IF (IRESP == 0) PACPRS(:,:)=Z2D(:,:)/(1000.)
  CASE ('INIT')
    PINPRS(:,:) = 0.0
    PACPRS(:,:) = 0.0
  END SELECT
END IF
!
IF (SIZE(PINPRG) /= 0 ) THEN
  SELECT CASE(YGETRGT)
  CASE ('READ')
    CALL FIND_FIELD_ID_FROM_MNHNAME('INPRG',IID,IRESP)
    TZFIELD = TFIELDMETADATA( TFIELDLIST(IID) )
    TZFIELD%CUNITS = 'mm hour-1'
    CALL IO_Field_read(TPINIFILE,TZFIELD,Z2D,IRESP)
    IF (IRESP == 0) PINPRG(:,:)=Z2D(:,:)/(1000.*3600.)
    !
    CALL FIND_FIELD_ID_FROM_MNHNAME('ACPRG',IID,IRESP)
    TZFIELD = TFIELDMETADATA( TFIELDLIST(IID) )
    TZFIELD%CUNITS = 'mm'
    CALL IO_Field_read(TPINIFILE,TZFIELD,Z2D,IRESP)
    IF (IRESP == 0) PACPRG(:,:)=Z2D(:,:)/(1000.)
  CASE ('INIT')
    PINPRG(:,:) = 0.0
    PACPRG(:,:) = 0.0
  END SELECT
END IF
!
IF (SIZE(PINPRH) /= 0 ) THEN
  SELECT CASE(YGETRHT)
  CASE ('READ')
    CALL FIND_FIELD_ID_FROM_MNHNAME('INPRH',IID,IRESP)
    TZFIELD = TFIELDMETADATA( TFIELDLIST(IID) )
    TZFIELD%CUNITS = 'mm hour-1'
    CALL IO_Field_read(TPINIFILE,TZFIELD,Z2D,IRESP)
    IF (IRESP == 0) PINPRH(:,:)=Z2D(:,:)/(1000.*3600.)
    !
    CALL FIND_FIELD_ID_FROM_MNHNAME('ACPRH',IID,IRESP)
    TZFIELD = TFIELDMETADATA( TFIELDLIST(IID) )
    TZFIELD%CUNITS = 'mm'
    CALL IO_Field_read(TPINIFILE,TZFIELD,Z2D,IRESP)
    IF (IRESP == 0) PACPRH(:,:)=Z2D(:,:)/(1000.)
  CASE ('INIT')
    PINPRH(:,:) = 0.0
    PACPRH(:,:) = 0.0
  END SELECT
END IF
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE READ_PRECIP_FIELD
