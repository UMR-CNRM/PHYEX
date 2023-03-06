!MNH_LIC Copyright 2022-2023 CNRS, Meteo-France and Universite Paul Sabatier
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
 USE MODD_FIELD, ONLY: TFIELDMETADATA
 IMPLICIT NONE
 INTERFACE IO_Field_write_phy
    MODULE PROCEDURE IO_Field_write_phy_byfield_X2, IO_Field_write_phy_byfield_X1
 END INTERFACE
CONTAINS
  SUBROUTINE IO_Field_write_phy_byfield_X2(D, TPFILE, TPFIELD, PFIELD, KRESP, koffset )
    USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
    !
    IMPLICIT NONE
    !
    !*      0.1   Declarations of arguments
    !
    TYPE(DIMPHYEX_t),                     INTENT(IN)  :: D
    TYPE(TFILEDATA),                      INTENT(IN)  :: TPFILE
    TYPE(TFIELDMETADATA),                 INTENT(IN)  :: TPFIELD
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
    TYPE(TFILEDATA),                      INTENT(IN)  :: TPFILE
    TYPE(TFIELDMETADATA),                 INTENT(IN)  :: TPFIELD
    REAL,DIMENSION(D%NIT,D%NJT,D%NKT),TARGET,  INTENT(IN)  :: PFIELD   ! array containing the data field
    INTEGER,                    OPTIONAL, INTENT(OUT) :: KRESP    ! return-code
    integer, dimension(3),      optional, intent(in)  :: koffset
    !
    CALL IO_Field_write(TPFILE, TPFIELD, PFIELD, KRESP, koffset )
    !
  END SUBROUTINE IO_Field_write_phy_unpack2D
!
  SUBROUTINE IO_Field_write_phy_byfield_X1(D, TPFILE, TPFIELD, PFIELD, KRESP, koffset )
    USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
    !
    IMPLICIT NONE
    !
    !*      0.1   Declarations of arguments
    !
    TYPE(DIMPHYEX_t),                     INTENT(IN)  :: D
    TYPE(TFILEDATA),                      INTENT(IN)  :: TPFILE
    TYPE(TFIELDMETADATA),                 INTENT(IN)  :: TPFIELD
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
    TYPE(TFILEDATA),                     INTENT(IN)  :: TPFILE
    TYPE(TFIELDMETADATA),                INTENT(IN)  :: TPFIELD
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
