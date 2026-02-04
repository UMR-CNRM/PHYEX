!MNH_LIC Copyright 2022-2026 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
! Modifications:
!  Q.Rodier 02/2023 Creation call to mode_io_field_write inside PHYEX
!-----------------------------------------------------------------
!
MODULE MODE_IO_FIELD_WRITE_PHY
 USE MODD_IO,         ONLY: TFILEDATA
 USE MODD_FIELD,      ONLY: TFIELDLIST, TFIELDMETADATA

 USE MODE_FIELD,      ONLY: FIND_FIELD_ID_FROM_MNHNAME

 IMPLICIT NONE

 INTERFACE IO_Field_write_phy
    MODULE PROCEDURE IO_Field_write_phy_byname_X2,  IO_Field_write_phy_byname_X1, &
                     IO_Field_write_phy_byfield_X2, IO_Field_write_phy_byfield_X1
 END INTERFACE
CONTAINS
  SUBROUTINE IO_Field_write_phy_byname_X2(D, TPFILE, HNAME, PFIELD, KRESP, koffset )
    USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
    !
    IMPLICIT NONE
    !
    !*      0.1   Declarations of arguments
    !
    TYPE(DIMPHYEX_t),                     INTENT(IN)  :: D
    TYPE(TFILEDATA),                      INTENT(INOUT)  :: TPFILE
    CHARACTER(LEN=*),                     INTENT(IN)  :: HNAME    ! name of the field to write
    REAL,DIMENSION(D%NIJT,D%NKT),TARGET,  INTENT(IN)  :: PFIELD   ! array containing the data field
    INTEGER,                    OPTIONAL, INTENT(OUT) :: KRESP    ! return-code
    integer, dimension(3),      optional, intent(in)  :: koffset
    !
    !*      0.2   Declarations of local variables
    !
    INTEGER :: ID    ! Index of the field
    INTEGER :: IRESP ! return-code
    !
    CALL FIND_FIELD_ID_FROM_MNHNAME( HNAME, ID, IRESP )
    !
    IF( IRESP == 0 ) CALL IO_FIELD_WRITE_PHY_UNPACK2D( D, TPFILE, TFIELDLIST(ID), PFIELD, IRESP, KOFFSET )
    !
    IF (PRESENT(KRESP)) KRESP = IRESP
    !
  END SUBROUTINE IO_Field_write_phy_byname_X2
!
  SUBROUTINE IO_Field_write_phy_byfield_X2(D, TPFILE, TPFIELD, PFIELD, KRESP, koffset )
    USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
    !
    IMPLICIT NONE
    !
    !*      0.1   Declarations of arguments
    !
    TYPE(DIMPHYEX_t),                     INTENT(IN)  :: D
    TYPE(TFILEDATA),                      INTENT(INOUT)  :: TPFILE
    CLASS(TFIELDMETADATA),                INTENT(INOUT)  :: TPFIELD
    REAL,DIMENSION(D%NIJT,D%NKT),TARGET,  INTENT(IN)  :: PFIELD   ! array containing the data field
    INTEGER,                    OPTIONAL, INTENT(OUT) :: KRESP    ! return-code
    integer, dimension(3),      optional, intent(in)  :: koffset
    !
    CALL IO_Field_write_phy_unpack2D(D,TPFILE, TPFIELD, PFIELD, KRESP, koffset )
    !
  END SUBROUTINE IO_Field_write_phy_byfield_X2
!
  SUBROUTINE IO_Field_write_phy_unpack2D(D, TPFILE, TPFIELD, PFIELD, KRESP, koffset )
    USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
    USE MODE_IO_FIELD_WRITE, ONLY: IO_Field_write
    !
    IMPLICIT NONE
    !
    !*      0.1   Declarations of arguments
    !
    TYPE(DIMPHYEX_t),                     INTENT(IN)  :: D
    TYPE(TFILEDATA),                      INTENT(INOUT)  :: TPFILE
    CLASS(TFIELDMETADATA),                INTENT(INOUT)  :: TPFIELD
    REAL,DIMENSION(D%NIT,D%NJT,D%NKT),TARGET,  INTENT(IN)  :: PFIELD   ! array containing the data field
    INTEGER,                    OPTIONAL, INTENT(OUT) :: KRESP    ! return-code
    integer, dimension(3),      optional, intent(in)  :: koffset
    !
    CALL IO_Field_write(TPFILE, TPFIELD, PFIELD, KRESP, koffset )
    !
  END SUBROUTINE IO_Field_write_phy_unpack2D
!
  SUBROUTINE IO_Field_write_phy_byname_X1(D, TPFILE, HNAME, PFIELD, KRESP, koffset )
    USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
    !
    IMPLICIT NONE
    !
    !*      0.1   Declarations of arguments
    !
    TYPE(DIMPHYEX_t),                     INTENT(IN)  :: D
    TYPE(TFILEDATA),                      INTENT(INOUT)  :: TPFILE
    CHARACTER(LEN=*),                     INTENT(IN)  :: HNAME    ! name of the field to write
    REAL,DIMENSION(D%NIJT),TARGET,  INTENT(IN)  :: PFIELD   ! array containing the data field
    INTEGER,                    OPTIONAL, INTENT(OUT) :: KRESP    ! return-code
    integer, dimension(3),      optional, intent(in)  :: koffset
    !
    !*      0.2   Declarations of local variables
    !
    INTEGER :: ID    ! Index of the field
    INTEGER :: IRESP ! return-code
    !
    CALL FIND_FIELD_ID_FROM_MNHNAME( HNAME, ID, IRESP )
    !
    IF( IRESP == 0 ) CALL IO_FIELD_WRITE_PHY_UNPACK1D( D, TPFILE, TFIELDLIST(ID), PFIELD, IRESP, KOFFSET )
    !
    IF (PRESENT(KRESP)) KRESP = IRESP
    !
  END SUBROUTINE IO_Field_write_phy_byname_X1
!
  SUBROUTINE IO_Field_write_phy_byfield_X1(D, TPFILE, TPFIELD, PFIELD, KRESP, koffset )
    USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
    !
    IMPLICIT NONE
    !
    !*      0.1   Declarations of arguments
    !
    TYPE(DIMPHYEX_t),                     INTENT(IN)  :: D
    TYPE(TFILEDATA),                      INTENT(INOUT)  :: TPFILE
    CLASS(TFIELDMETADATA),                INTENT(INOUT)  :: TPFIELD
    REAL,DIMENSION(D%NIJT),TARGET,  INTENT(IN)  :: PFIELD   ! array containing the data field
    INTEGER,                    OPTIONAL, INTENT(OUT) :: KRESP    ! return-code
    integer, dimension(3),      optional, intent(in)  :: koffset
    !
    CALL IO_Field_write_phy_unpack1D(D,TPFILE, TPFIELD, PFIELD, KRESP, koffset )
    !
  END SUBROUTINE IO_Field_write_phy_byfield_X1
!
  SUBROUTINE IO_Field_write_phy_unpack1D(D, TPFILE, TPFIELD, PFIELD, KRESP, koffset )
    USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
    USE MODE_IO_FIELD_WRITE, ONLY: IO_Field_write
    !
    IMPLICIT NONE
    !
    !*      0.1   Declarations of arguments
    !
    TYPE(DIMPHYEX_t),                    INTENT(IN)  :: D
    TYPE(TFILEDATA),                     INTENT(INOUT)  :: TPFILE
    CLASS(TFIELDMETADATA),               INTENT(INOUT)  :: TPFIELD
    REAL,DIMENSION(D%NIT,D%NJT),TARGET,  INTENT(IN)  :: PFIELD   ! array containing the data field
    INTEGER,                    OPTIONAL, INTENT(OUT) :: KRESP    ! return-code
    integer, dimension(3),      optional, intent(in)  :: koffset
    !
    CALL IO_Field_write(TPFILE, TPFIELD, PFIELD, KRESP, koffset )
    !
  END SUBROUTINE IO_Field_write_phy_unpack1D
!
!
END MODULE MODE_IO_FIELD_WRITE_PHY
