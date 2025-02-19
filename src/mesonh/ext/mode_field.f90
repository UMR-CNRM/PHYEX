!MNH_LIC Copyright 2016-2023 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
! Original version:
!  P. Wautelet: 05/2016-04/2018: new data structures and calls for I/O
! Modifications:
!  P. Wautelet 29/01/2019: small bug correction (null pointers) in FIELDLIST_GOTO_MODEL if NESPGD or PGD
!  P. Wautelet 01/02/2019: bug correction in case XRT is not associated
!  C. Lac         02/2019: add rain fraction as an output field
!  S. Bielli      02/2019: sea salt: significant sea wave height influences salt emission; 5 salt modes
!  P. Wautelet 06/03/2019: correct ZWS entry
!  P. Wautelet 06/06/2019: bug correction in FIELDLIST_GOTO_MODEL (XLSTHM was overwritten if LUSERV=.FALSE. due to wrong IF block)
!  P. Wautelet 19/06/2019: add Fieldlist_nmodel_resize subroutine + provide KMODEL to INI_FIELD_LIST when known
!  P. Wautelet 23/01/2020: split in modd_field.f90 and mode_field.f90
!  JL Redelsperger 03/2021: add variables for Ocean LES and auto-coupled version
!  P. Wautelet 08/10/2021: add Goto_model_1field + Add_field2list procedures + remove Fieldlist_nmodel_resize
!  P. Wautelet 14/10/2021: dynamically allocate tfieldlist (+ reallocate if necessary)
!  A. Costes      12/2021: add Blaze fire model variables
!  A. Marcel Jan 2025: EDMF contribution to dynamic TKE production
!-----------------------------------------------------------------
module mode_field

use modd_conf,       only: cprogram
use modd_field
use modd_io,         only: NVERB_DEBUG, NVERB_INFO, NVERB_WARNING, NVERB_ERROR, NVERB_FATAL
use modd_parameters, only: JPMODELMAX
use modd_fire_n,     only: NREFINX, NREFINY

use mode_msg

implicit none

private

public :: Ini_field_list
public :: Find_field_id_from_mnhname
public :: Alloc_field_scalars
public :: Fieldlist_goto_model
public :: Ini_field_scalars

interface Goto_model_1field
  module procedure :: Goto_model_1field_c0d
  module procedure :: Goto_model_1field_c1d
  module procedure :: Goto_model_1field_l0d
  module procedure :: Goto_model_1field_l1d
  module procedure :: Goto_model_1field_n0d
  module procedure :: Goto_model_1field_n1d
  module procedure :: Goto_model_1field_n2d
  module procedure :: Goto_model_1field_n3d
  module procedure :: Goto_model_1field_t0d
  module procedure :: Goto_model_1field_t1d
  module procedure :: Goto_model_1field_x0d
  module procedure :: Goto_model_1field_x1d
  module procedure :: Goto_model_1field_x2d
  module procedure :: Goto_model_1field_x3d
  module procedure :: Goto_model_1field_x4d
  module procedure :: Goto_model_1field_x5d
  module procedure :: Goto_model_1field_x6d
end interface


contains

SUBROUTINE INI_FIELD_LIST()
! Modif
!  J.Escobar 25/04/2018: missing def of FRC
!------------------------------------------------

CHARACTER(LEN=64) :: YMSG

CALL PRINT_MSG(NVERB_DEBUG,'GEN','INI_FIELD_LIST','called')
IF (LFIELDLIST_ISINIT) THEN
  CALL PRINT_MSG(NVERB_ERROR,'GEN','INI_FIELD_LIST','already called')
  RETURN
END IF

LFIELDLIST_ISINIT = .TRUE.

Allocate( tfieldlist(NMAXFIELDINIT) )
NMAXFIELDS = NMAXFIELDINIT

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'MNHVERSION',     &
  CSTDNAME   = '',               &
  CLONGNAME  = 'MesoNH version', &
  CUNITS     = '',               &
  CDIR       = '--',             &
  CCOMMENT   = '',               &
  NGRID      = 0,                &
  NTYPE      = TYPEINT,          &
  NDIMS      = 1,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA(                  &
  CMNHNAME   = 'MASDEV',                          &
  CSTDNAME   = '',                                &
  CLONGNAME  = 'MesoNH version (without bugfix)', &
  CUNITS     = '',                                &
  CDIR       = '--',                              &
  CCOMMENT   = '',                                &
  NGRID      = 0,                                 &
  NTYPE      = TYPEINT,                           &
  NDIMS      = 0,                                 &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA(       &
  CMNHNAME   = 'BUGFIX',               &
  CSTDNAME   = '',                     &
  CLONGNAME  = 'MesoNH bugfix number', &
  CUNITS     = '',                     &
  CDIR       = '--',                   &
  CCOMMENT   = '',                     &
  NGRID      = 0,                      &
  NTYPE      = TYPEINT,                &
  NDIMS      = 0,                      &
  LTIMEDEP   = .FALSE.                 ) )

call Add_field2list( TFIELDDATA(              &
  CMNHNAME   = 'BIBUSER',                     &
  CSTDNAME   = '',                            &
  CLONGNAME  = 'MesoNH: user binary library', &
  CUNITS     = '',                            &
  CDIR       = '--',                          &
  CCOMMENT   = '',                            &
  NGRID      = 0,                             &
  NTYPE      = TYPECHAR,                      &
  NDIMS      = 0,                             &
  LTIMEDEP   = .FALSE.                        ) )

call Add_field2list( TFIELDDATA(               &
  CMNHNAME   = 'VERSION',                      &
  CSTDNAME   = '',                             &
  CLONGNAME  = 'SURFEX version (without BUG)', &
  CUNITS     = '',                             &
  CDIR       = '--',                           &
  CCOMMENT   = '',                             &
  NGRID      = 0,                              &
  NTYPE      = TYPEINT,                        &
  NDIMS      = 0,                              &
  LTIMEDEP   = .FALSE.                         ) )

call Add_field2list( TFIELDDATA(       &
  CMNHNAME   = 'BUG',                  &
  CSTDNAME   = '',                     &
  CLONGNAME  = 'SURFEX bugfix number', &
  CUNITS     = '',                     &
  CDIR       = '--',                   &
  CCOMMENT   = '',                     &
  NGRID      = 0,                      &
  NTYPE      = TYPEINT,                &
  NDIMS      = 0,                      &
  LTIMEDEP   = .FALSE.                 ) )

call Add_field2list( TFIELDDATA(              &
  CMNHNAME   = 'PROGRAM',                     &
  CSTDNAME   = '',                            &
  CLONGNAME  = 'MesoNH family: used program', &
  CUNITS     = '',                            &
  CDIR       = '--',                          &
  CCOMMENT   = '',                            &
  NGRID      = 0,                             &
  NTYPE      = TYPECHAR,                      &
  NDIMS      = 0,                             &
  LTIMEDEP   = .FALSE.                        ) )

call Add_field2list( TFIELDDATA(    &
  CMNHNAME   = 'FILETYPE',          &
  CSTDNAME   = '',                  &
  CLONGNAME  = 'type of this file', &
  CUNITS     = '',                  &
  CDIR       = '--',                &
  CCOMMENT   = '',                  &
  NGRID      = 0,                   &
  NTYPE      = TYPECHAR,            &
  NDIMS      = 0,                   &
  LTIMEDEP   = .FALSE.              ) )

call Add_field2list( TFIELDDATA(          &
  CMNHNAME   = 'MY_NAME',                 &
  CSTDNAME   = '',                        &
  CLONGNAME  = 'filename (no extension)', &
  CUNITS     = '',                        &
  CDIR       = '--',                      &
  CCOMMENT   = '',                        &
  NGRID      = 0,                         &
  NTYPE      = TYPECHAR,                  &
  NDIMS      = 0,                         &
  LTIMEDEP   = .FALSE.                    ) )

call Add_field2list( TFIELDDATA(           &
  CMNHNAME   = 'DAD_NAME',                 &
  CSTDNAME   = '',                         &
  CLONGNAME  = 'filename of the dad file', &
  CUNITS     = '',                         &
  CDIR       = '--',                       &
  CCOMMENT   = '',                         &
  NGRID      = 0,                          &
  NTYPE      = TYPECHAR,                   &
  NDIMS      = 0,                          &
  LTIMEDEP   = .FALSE.                     ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'DXRATIO',        &
  CSTDNAME   = '',               &
  CLONGNAME  = 'DXRATIO',        &
  CUNITS     = '',               &
  CDIR       = '--',             &
  CCOMMENT   = 'Resolution ratio between this mesh and its father in x-direction', &
  NGRID      = 0,                &
  NTYPE      = TYPEINT,          &
  NDIMS      = 0,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'DYRATIO',        &
  CSTDNAME   = '',               &
  CLONGNAME  = 'DYRATIO',        &
  CUNITS     = '',               &
  CDIR       = '--',             &
  CCOMMENT   = 'Resolution ratio between this mesh and its father in y-direction', &
  NGRID      = 0,                &
  NTYPE      = TYPEINT,          &
  NDIMS      = 0,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'XSIZE',          &
  CSTDNAME   = '',               &
  CLONGNAME  = 'XSIZE',          &
  CUNITS     = '',               &
  CDIR       = '--',             &
  CCOMMENT   = 'Number of model 1 grid points in x-direction in the model 2 physical domain', &
  NGRID      = 0,                &
  NTYPE      = TYPEINT,          &
  NDIMS      = 0,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'YSIZE',          &
  CSTDNAME   = '',               &
  CLONGNAME  = 'YSIZE',          &
  CUNITS     = '',               &
  CDIR       = '--',             &
  CCOMMENT   = 'Number of model 1 grid points in y-direction in the model 2 physical domain', &
  NGRID      = 0,                &
  NTYPE      = TYPEINT,          &
  NDIMS      = 0,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'XOR',            &
  CSTDNAME   = '',               &
  CLONGNAME  = 'XOR',            &
  CUNITS     = '',               &
  CDIR       = '--',             &
  CCOMMENT   = 'Horizontal position of this mesh relative to its father', &
  NGRID      = 0,                &
  NTYPE      = TYPEINT,          &
  NDIMS      = 0,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'YOR',            &
  CSTDNAME   = '',               &
  CLONGNAME  = 'YOR',            &
  CUNITS     = '',               &
  CDIR       = '--',             &
  CCOMMENT   = 'Vertical position of this mesh relative to its father', &
  NGRID      = 0,                &
  NTYPE      = TYPEINT,          &
  NDIMS      = 0,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'STORAGE_TYPE',   &
  CSTDNAME   = '',               &
  CLONGNAME  = 'STORAGE_TYPE',   &
  CUNITS     = '',               &
  CDIR       = '--',             &
  CCOMMENT   = 'Storage type for the information written in the FM files', &
  NGRID      = 0,                &
  NTYPE      = TYPECHAR,         &
  NDIMS      = 0,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'IMAX',           &
  CSTDNAME   = '',               &
  CLONGNAME  = 'IMAX',           &
  CUNITS     = '',               &
  CDIR       = '--',             &
  CCOMMENT   = 'x-dimension of the physical domain', &
  NGRID      = 0,                &
  NTYPE      = TYPEINT,          &
  NDIMS      = 0,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'JMAX',           &
  CSTDNAME   = '',               &
  CLONGNAME  = 'JMAX',           &
  CUNITS     = '',               &
  CDIR       = '--',             &
  CCOMMENT   = 'y-dimension of the physical domain', &
  NGRID      = 0,                &
  NTYPE      = TYPEINT,          &
  NDIMS      = 0,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'KMAX',           &
  CSTDNAME   = '',               &
  CLONGNAME  = 'KMAX',           &
  CUNITS     = '',               &
  CDIR       = '--',             &
  CCOMMENT   = 'z-dimension of the physical domain', &
  NGRID      = 0,                &
  NTYPE      = TYPEINT,          &
  NDIMS      = 0,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'JPHEXT',         &
  CSTDNAME   = '',               &
  CLONGNAME  = 'JPHEXT',         &
  CUNITS     = '',               &
  CDIR       = '--',             &
  CCOMMENT   = 'Number of horizontal external points on each side', &
  NGRID      = 0,                &
  NTYPE      = TYPEINT,          &
  NDIMS      = 0,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'RPK',            &
  CSTDNAME   = '',               &
  CLONGNAME  = 'RPK',            &
  CUNITS     = '',               &
  CDIR       = '--',             &
  CCOMMENT   = 'Projection parameter for conformal projection', &
  NGRID      = 0,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 0,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'LONORI',         &
  CSTDNAME   = '',               &
  CLONGNAME  = 'LONORI',         &
  CUNITS     = 'degree',         &
  CDIR       = '--',             &
  CCOMMENT   = 'Longitude of the point of coordinates x=0, y=0 for conformal projection', &
  NGRID      = 0,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 0,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'LATORI',         &
  CSTDNAME   = '',               &
  CLONGNAME  = 'LATORI',         &
  CUNITS     = 'degree',         &
  CDIR       = '--',             &
  CCOMMENT   = 'Latitude of the point of coordinates x=0, y=0 for conformal projection', &
  NGRID      = 0,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 0,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'LONOR',          &
  CSTDNAME   = '',               &
  CLONGNAME  = 'LONOR',          &
  CUNITS     = 'degree',         &
  CDIR       = '--',             &
  CCOMMENT   = 'Longitude of 1st mass point', &
  NGRID      = 0,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 0,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'LATOR',          &
  CSTDNAME   = '',               &
  CLONGNAME  = 'LATOR',          &
  CUNITS     = 'degree',         &
  CDIR       = '--',             &
  CCOMMENT   = 'Latitude of 1st mass point', &
  NGRID      = 0,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 0,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'THINSHELL',      &
  CSTDNAME   = '',               &
  CLONGNAME  = 'THINSHELL',      &
  CUNITS     = '',               &
  CDIR       = '--',             &
  CCOMMENT   = 'Logical for thinshell approximation', &
  NGRID      = 0,                &
  NTYPE      = TYPELOG,          &
  NDIMS      = 0,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'LAT0',           &
  CSTDNAME   = '',               &
  CLONGNAME  = 'LAT0',           &
  CUNITS     = 'degree',         &
  CDIR       = '--',             &
  CCOMMENT   = 'Reference latitude for conformal projection', &
  NGRID      = 0,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 0,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'LON0',           &
  CSTDNAME   = '',               &
  CLONGNAME  = 'LON0',           &
  CUNITS     = 'degree',         &
  CDIR       = '--',             &
  CCOMMENT   = 'Reference longitude for conformal projection', &
  NGRID      = 0,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 0,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'BETA',           &
  CSTDNAME   = '',               &
  CLONGNAME  = 'BETA',           &
  CUNITS     = 'degree',         &
  CDIR       = '--',             &
  CCOMMENT   = 'Rotation angle for conformal projection', &
  NGRID      = 0,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 0,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'XHAT',           &
  CSTDNAME   = 'projection_x_coordinate', &
  CLONGNAME  = 'XHAT',           &
  CUNITS     = 'm',              &
  CDIR       = 'XX',             &
  CCOMMENT   = 'Position x in the conformal or cartesian plane', &
  NGRID      = 2,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 1,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'YHAT',           &
  CSTDNAME   = 'projection_y_coordinate', &
  CLONGNAME  = 'YHAT',           &
  CUNITS     = 'm',              &
  CDIR       = 'YY',             &
  CCOMMENT   = 'Position y in the conformal or cartesian plane', &
  NGRID      = 3,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 1,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'XHATM',          &
  CSTDNAME   = 'projection_x_coordinate', &
  CLONGNAME  = 'XHATM',          &
  CUNITS     = 'm',              &
  CDIR       = 'XX',             &
  CCOMMENT   = 'Position x in the conformal or cartesian plane at mass points', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 1,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'YHATM',          &
  CSTDNAME   = 'projection_y_coordinate', &
  CLONGNAME  = 'YHATM',          &
  CUNITS     = 'm',              &
  CDIR       = 'YY',             &
  CCOMMENT   = 'Position y in the conformal or cartesian plane at mass points', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 1,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'ZHAT',           &
!TODO: check stdname
  CSTDNAME   = '',               &
  CLONGNAME  = 'ZHAT',           &
  CUNITS     = 'm',              &
  CDIR       = 'ZZ',             &
  CCOMMENT   = 'Height level without orography', &
  NGRID      = 4,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 1,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'ZHATM',          &
!TODO: check stdname
  CSTDNAME   = '',               &
  CLONGNAME  = 'ZHATM',          &
  CUNITS     = 'm',              &
  CDIR       = 'ZZ',             &
  CCOMMENT   = 'Height level without orography at mass point', &
  NGRID      = 4,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 1,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'HAT_BOUND', &
  CSTDNAME   = '',               &
  CLONGNAME  = 'HAT_BOUND', &
  CUNITS     = 'm',              &
  CDIR       = '--',             &
  CCOMMENT   = 'Boundaries of domain in the conformal or cartesian plane at u and v points', &
  NGRID      = 0,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 1,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA(  &
  CMNHNAME   = 'HATM_BOUND', &
  CSTDNAME   = '',                &
  CLONGNAME  = 'HATM_BOUND', &
  CUNITS     = 'm',               &
  CDIR       = '--',              &
  CCOMMENT   = 'Boundaries of domain in the conformal or cartesian plane at mass points', &
  NGRID      = 0,                 &
  NTYPE      = TYPEREAL,          &
  NDIMS      = 1,                 &
  LTIMEDEP   = .FALSE.            ) )

call Add_field2list( TFIELDDATA(      &
  CMNHNAME   = 'ZTOP',                &
  CSTDNAME   = 'altitude_at_top_of_atmosphere_model', &
  CLONGNAME  = 'ZTOP',                &
  CUNITS     = 'm',                   &
  CDIR       = '--',                  &
  CCOMMENT   = 'Height of top level', &
  NGRID      = 4,                     &
  NTYPE      = TYPEREAL,              &
  NDIMS      = 0,                     &
  LTIMEDEP   = .FALSE.                ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'DXHAT',          &
!TODO: check stdname
  CSTDNAME   = '',               &
  CLONGNAME  = 'DXHAT',          &
  CUNITS     = 'm',              &
  CDIR       = 'XX',             &
  CCOMMENT   = 'Horizontal stretching in x', &
  NGRID      = 2,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 1,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'DYHAT',          &
!TODO: check stdname
  CSTDNAME   = '',               &
  CLONGNAME  = 'DYHAT',          &
  CUNITS     = 'm',              &
  CDIR       = 'YY',             &
  CCOMMENT   = 'Horizontal stretching in y', &
  NGRID      = 3,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 1,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'ALT',            &
  CSTDNAME   = 'altitude',       &
  CLONGNAME  = 'ALT',            &
  CUNITS     = 'm',              &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Z_ALTitude', &
  NGRID      = 4,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'DIRCOSXW',       &
  CSTDNAME   = '',               &
  CLONGNAME  = 'DIRCOSXW',       &
  CUNITS     = '',               &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X director cosinus of the normal to the ground surface', &
  NGRID      = 4,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 2,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'DIRCOSYW',       &
  CSTDNAME   = '',               &
  CLONGNAME  = 'DIRCOSYW',       &
  CUNITS     = '',               &
  CDIR       = 'XY',             &
  CCOMMENT   = 'Y director cosinus of the normal to the ground surface', &
  NGRID      = 4,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 2,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'DIRCOSZW',       &
  CSTDNAME   = '',               &
  CLONGNAME  = 'DIRCOSZW',       &
  CUNITS     = '',               &
  CDIR       = 'XY',             &
  CCOMMENT   = 'Z director cosinus of the normal to the ground surface', &
  NGRID      = 4,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 2,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'COSSLOPE',       &
  CSTDNAME   = '',               &
  CLONGNAME  = 'COSSLOPE',       &
  CUNITS     = '',               &
  CDIR       = 'XY',             &
  CCOMMENT   = 'cosinus of the angle between i and the slope vector', &
  NGRID      = 4,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 2,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'SINSLOPE',       &
  CSTDNAME   = '',               &
  CLONGNAME  = 'SINSLOPE',       &
  CUNITS     = '',               &
  CDIR       = 'XY',             &
  CCOMMENT   = 'sinus of the angle between i and the slope vector', &
  NGRID      = 4,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 2,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'MAP',            &
!TODO: check stdname
  CSTDNAME   = '',               &
  CLONGNAME  = 'MAP',            &
  CUNITS     = '',               &
  CDIR       = 'XY',             &
  CCOMMENT   = 'Map factor',     &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 2,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'latitude',       &
  CSTDNAME   = 'latitude',       &
  CLONGNAME  = 'latitude',       &
  CUNITS     = 'degrees_north',  &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_latitude at mass point', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 2,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'longitude',      &
  CSTDNAME   = 'longitude',      &
  CLONGNAME  = 'longitude',      &
  CUNITS     = 'degrees_east',   &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_longitude at mass point', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 2,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA(          &
  CMNHNAME   = 'latitude_u',              &
  CSTDNAME   = 'latitude_at_u_location',  &
  CLONGNAME  = 'latitude at u location',  &
  CUNITS     = 'degrees_north',           &
  CDIR       = 'XY',                      &
  CCOMMENT   = 'X_Y_latitude at u point', &
  NGRID      = 2,                         &
  NTYPE      = TYPEREAL,                  &
  NDIMS      = 2,                         &
  LTIMEDEP   = .FALSE.                    ) )

call Add_field2list( TFIELDDATA(           &
  CMNHNAME   = 'longitude_u',              &
  CSTDNAME   = 'longitude_at_u_location',  &
  CLONGNAME  = 'longitude at u location',  &
  CUNITS     = 'degrees_east',             &
  CDIR       = 'XY',                       &
  CCOMMENT   = 'X_Y_longitude at u point', &
  NGRID      = 2,                          &
  NTYPE      = TYPEREAL,                   &
  NDIMS      = 2,                          &
  LTIMEDEP   = .FALSE.                     ) )

call Add_field2list( TFIELDDATA(          &
  CMNHNAME   = 'latitude_v',              &
  CSTDNAME   = 'latitude_at_v_location',  &
  CLONGNAME  = 'latitude at v location',  &
  CUNITS     = 'degrees_north',           &
  CDIR       = 'XY',                      &
  CCOMMENT   = 'X_Y_latitude at v point', &
  NGRID      = 3,                         &
  NTYPE      = TYPEREAL,                  &
  NDIMS      = 2,                         &
  LTIMEDEP   = .FALSE.                    ) )

call Add_field2list( TFIELDDATA(           &
  CMNHNAME   = 'longitude_v',              &
  CSTDNAME   = 'longitude_at_v_location',  &
  CLONGNAME  = 'longitude at v location',  &
  CUNITS     = 'degrees_east',             &
  CDIR       = 'XY',                       &
  CCOMMENT   = 'X_Y_longitude at v point', &
  NGRID      = 3,                          &
  NTYPE      = TYPEREAL,                   &
  NDIMS      = 2,                          &
  LTIMEDEP   = .FALSE.                     ) )

call Add_field2list( TFIELDDATA(          &
  CMNHNAME   = 'latitude_f',              &
  CSTDNAME   = 'latitude_at_f_location',  &
  CLONGNAME  = 'latitude at f location',  &
  CUNITS     = 'degrees_north',           &
  CDIR       = 'XY',                      &
  CCOMMENT   = 'X_Y_latitude at f point', &
  NGRID      = 5,                         &
  NTYPE      = TYPEREAL,                  &
  NDIMS      = 2,                         &
  LTIMEDEP   = .FALSE.                    ) )

call Add_field2list( TFIELDDATA(           &
  CMNHNAME   = 'longitude_f',              &
  CSTDNAME   = 'longitude_at_f_location',  &
  CLONGNAME  = 'longitude at f location',  &
  CUNITS     = 'degrees_east',             &
  CDIR       = 'XY',                       &
  CCOMMENT   = 'X_Y_longitude at f point', &
  NGRID      = 5,                          &
  NTYPE      = TYPEREAL,                   &
  NDIMS      = 2,                          &
  LTIMEDEP   = .FALSE.                     ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'LAT',            &
!   CSTDNAME   = 'latitude',        &
  CSTDNAME   = '',               &
  CLONGNAME  = 'LAT',            &
  CUNITS     = 'degrees_north',  &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_latitude',   &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 2,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'LON',            &
!   CSTDNAME   = 'longitude',       &
  CSTDNAME   = '',               &
  CLONGNAME  = 'LON',            &
  CUNITS     = 'degrees_east',   &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_longitude',  &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 2,                &
  LTIMEDEP   = .FALSE.           ) )

!Note: do not use XHAT_ll in I/O (use XHAT instead)
call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'XHAT_ll',        &
  CSTDNAME   = 'projection_x_coordinate', &
  CLONGNAME  = 'XHAT_ll',        &
  CUNITS     = 'm',              &
!PW:TODO?: create a new field to say if the variable is distributed? and how (X,Y,XY...)?
  CDIR       = 'XX',             &
  CCOMMENT   = 'Position x in the conformal or cartesian plane (all domain)', &
  NGRID      = 2,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 1,                &
  LTIMEDEP   = .FALSE.           ) )

!Note: do not use YHAT_ll in I/O (use YHAT instead)
call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'YHAT_ll',        &
  CSTDNAME   = 'projection_y_coordinate', &
  CLONGNAME  = 'YHAT_ll',        &
  CUNITS     = 'm',              &
  CDIR       = 'YY',             &
  CCOMMENT   = 'Position y in the conformal or cartesian plane (all domain)', &
  NGRID      = 3,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 1,                &
  LTIMEDEP   = .FALSE.           ) )

!Note: do not use XHATM_ll in I/O (use XHATM instead)
call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'XHATM_ll',        &
  CSTDNAME   = 'projection_x_coordinate', &
  CLONGNAME  = 'XHATL_ll',        &
  CUNITS     = 'm',              &
  CDIR       = 'XX',             &
  CCOMMENT   = 'Position x in the conformal or cartesian plane at mass points (all domain)', &
  NGRID      = 2,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 1,                &
  LTIMEDEP   = .FALSE.           ) )

!Note: do not use YHATM_ll in I/O (use YHATM instead)
call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'YHATM_ll',        &
  CSTDNAME   = 'projection_y_coordinate', &
  CLONGNAME  = 'YHATM_ll',        &
  CUNITS     = 'm',              &
  CDIR       = 'YY',             &
  CCOMMENT   = 'Position y in the conformal or cartesian plane at mass points (all domain)', &
  NGRID      = 3,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 1,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA(   &
  CMNHNAME   = 'ZS',               &
  CSTDNAME   = 'surface_altitude', &
  CLONGNAME  = 'ZS',               &
  CUNITS     = 'm',                &
  CDIR       = 'XY',               &
  CCOMMENT   = 'orography',        &
  NGRID      = 4,                  &
  NTYPE      = TYPEREAL,           &
  NDIMS      = 2,                  &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA(  &
  CMNHNAME   = 'ZWS',             &
  CSTDNAME   = 'sea_surface_wave_significant_height', &
  CLONGNAME  = 'ZWS',             &
  CUNITS     = 'm',               &
  CDIR       = 'XY',              &
  CCOMMENT   = 'sea wave height', &
  NGRID      = 4,                 &
  NTYPE      = TYPEREAL,          &
  NDIMS      = 2,                 &
  LTIMEDEP   = .TRUE.             ) )

call Add_field2list( TFIELDDATA(   &
  CMNHNAME   = 'ZSMT',             &
  CSTDNAME   = '',                 &
  CLONGNAME  = 'ZSMT',             &
  CUNITS     = 'm',                &
  CDIR       = 'XY',               &
  CCOMMENT   = 'smooth orography', &
  NGRID      = 4,                  &
  NTYPE      = TYPEREAL,           &
  NDIMS      = 2,                  &
  LTIMEDEP   = .FALSE.             ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'SLEVE',          &
  CSTDNAME   = '',               &
  CLONGNAME  = 'SLEVE',          &
  CUNITS     = '',               &
  CDIR       = '--',             &
  CCOMMENT   = 'Logical for SLEVE coordinate', &
  NGRID      = 4,                &
  NTYPE      = TYPELOG,          &
  NDIMS      = 0,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'LEN1',           &
  CSTDNAME   = '',               &
  CLONGNAME  = 'LEN1',           &
  CUNITS     = '',               &
  CDIR       = '--',             &
  CCOMMENT   = 'Decay scale for smooth topography', &
  NGRID      = 4,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 0,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'LEN2',           &
  CSTDNAME   = '',               &
  CLONGNAME  = 'LEN2',           &
  CUNITS     = '',               &
  CDIR       = '--',             &
  CCOMMENT   = 'Decay scale for small-scale topography deviation', &
  NGRID      = 4,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 0,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA(                      &
  CMNHNAME   = 'DTMOD',                               &
  CSTDNAME   = '',                                    &
  CLONGNAME  = 'DTMOD',                               &
  CUNITS     = 'seconds since YYYY-MM-DD HH:MM:SS.S', &
  CDIR       = '--',                                  &
  CCOMMENT   = 'Time and date of model beginning',    &
  NGRID      = 0,                                     &
  NTYPE      = TYPEDATE,                              &
  NDIMS      = 0,                                     &
  LTIMEDEP   = .FALSE.                                ) )

call Add_field2list( TFIELDDATA(                      &
  CMNHNAME   = 'DTCUR',                               &
  CSTDNAME   = 'time',                                &
  CLONGNAME  = 'DTCUR',                               &
  CUNITS     = 'seconds since YYYY-MM-DD HH:MM:SS.S', &
  CDIR       = '--',                                  &
  CCOMMENT   = 'Current time and date',               &
  NGRID      = 0,                                     &
  NTYPE      = TYPEDATE,                              &
  NDIMS      = 0,                                     &
  LTIMEDEP   = .FALSE.                                ) )

call Add_field2list( TFIELDDATA(                            &
  CMNHNAME   = 'DTRAD_FULL',                                &
  CSTDNAME   = '',                                          &
  CLONGNAME  = 'DTRAD_FULL',                                &
  CUNITS     = 'seconds since YYYY-MM-DD HH:MM:SS.S',       &
  CDIR       = '--',                                        &
  CCOMMENT   = 'Time and date of last full radiation call', &
  NGRID      = 0,                                           &
  NTYPE      = TYPEDATE,                                    &
  NDIMS      = 0,                                           &
  LTIMEDEP   = .FALSE.                                      ) )

call Add_field2list( TFIELDDATA(                      &
  CMNHNAME   = 'DTRAD_CLLY',                          &
  CSTDNAME   = '',                                    &
  CLONGNAME  = 'DTRAD_CLLY',                          &
  CUNITS     = 'seconds since YYYY-MM-DD HH:MM:SS.S', &
  CDIR       = '--',                                  &
  CCOMMENT   = 'Time and date of last radiation call for only cloudy verticals', &
  NGRID      = 0,                                     &
  NTYPE      = TYPEDATE,                              &
  NDIMS      = 0,                                     &
  LTIMEDEP   = .FALSE.                                ) )

call Add_field2list( TFIELDDATA(                                 &
  CMNHNAME   = 'DTDCONV',                                        &
  CSTDNAME   = '',                                               &
  CLONGNAME  = 'DTDCONV',                                        &
  CUNITS     = 'seconds since YYYY-MM-DD HH:MM:SS.S',            &
  CDIR       = '--',                                             &
  CCOMMENT   = 'Time and date of the last deep convection call', &
  NGRID      = 0,                                                &
  NTYPE      = TYPEDATE,                                         &
  NDIMS      = 0,                                                &
  LTIMEDEP   = .FALSE.                                           ) )

call Add_field2list( TFIELDDATA(                        &
  CMNHNAME   = 'DTEXP',                                 &
  CSTDNAME   = '',                                      &
  CLONGNAME  = 'DTEXP',                                 &
  CUNITS     = 'seconds since YYYY-MM-DD HH:MM:SS.S',   &
  CDIR       = '--',                                    &
  CCOMMENT   = 'Time and date of experiment beginning', &
  NGRID      = 0,                                       &
  NTYPE      = TYPEDATE,                                &
  NDIMS      = 0,                                       &
  LTIMEDEP   = .FALSE.                                  ) )

call Add_field2list( TFIELDDATA(                      &
  CMNHNAME   = 'DTSEG',                               &
  CSTDNAME   = '',                                    &
  CLONGNAME  = 'DTSEG',                               &
  CUNITS     = 'seconds since YYYY-MM-DD HH:MM:SS.S', &
  CDIR       = '--',                                  &
  CCOMMENT   = 'Time and date of segment beginning',  &
  NGRID      = 0,                                     &
  NTYPE      = TYPEDATE,                              &
  NDIMS      = 0,                                     &
  LTIMEDEP   = .FALSE.                                ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'L1D',            &
  CSTDNAME   = '',               &
  CLONGNAME  = 'L1D',            &
  CUNITS     = '',               &
  CDIR       = '--',             &
  CCOMMENT   = 'Logical for 1D model version', &
  NGRID      = 0,                &
  NTYPE      = TYPELOG,          &
  NDIMS      = 0,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'L2D',            &
  CSTDNAME   = '',               &
  CLONGNAME  = 'L2D',            &
  CUNITS     = '',               &
  CDIR       = '--',             &
  CCOMMENT   = 'Logical for 2D model version', &
  NGRID      = 0,                &
  NTYPE      = TYPELOG,          &
  NDIMS      = 0,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'PACK',           &
  CSTDNAME   = '',               &
  CLONGNAME  = 'PACK',           &
  CUNITS     = '',               &
  CDIR       = '--',             &
  CCOMMENT   = 'Logical to compress 1D or 2D FM files', &
  NGRID      = 0,                &
  NTYPE      = TYPELOG,          &
  NDIMS      = 0,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'CARTESIAN',      &
  CSTDNAME   = '',               &
  CLONGNAME  = 'CARTESIAN',      &
  CUNITS     = '',               &
  CDIR       = '--',             &
  CCOMMENT   = 'Logical for cartesian geometry', &
  NGRID      = 0,                &
  NTYPE      = TYPELOG,          &
  NDIMS      = 0,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'LBOUSS',         &
  CSTDNAME   = '',               &
  CLONGNAME  = 'LBOUSS',         &
  CUNITS     = '',               &
  CDIR       = '--',             &
  CCOMMENT   = 'Logical for Boussinesq approximation', &
  NGRID      = 0,                &
  NTYPE      = TYPELOG,          &
  NDIMS      = 0,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'LOCEAN',         &
  CSTDNAME   = '',               &
  CLONGNAME  = 'LOCEAN',         &
  CUNITS     = '',               &
  CDIR       = '--',             &
  CCOMMENT   = 'Logical for Ocean MesoNH', &
  NGRID      = 0,                &
  NTYPE      = TYPELOG,          &
  NDIMS      = 0,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'LCOUPLES',       &
  CSTDNAME   = '',               &
  CLONGNAME  = 'LCOUPLES',       &
  CUNITS     = '',               &
  CDIR       = '--',             &
  CCOMMENT   = 'Logical for coupling O-A LES', &
  NGRID      = 0,                &
  NTYPE      = TYPELOG,          &
  NDIMS      = 0,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'SURF',           &
  CSTDNAME   = '',               &
  CLONGNAME  = 'SURF',           &
  CUNITS     = '',               &
  CDIR       = '--',             &
  CCOMMENT   = 'Kind of surface processes parameterization', &
  NGRID      = 0,                &
  NTYPE      = TYPECHAR,         &
  NDIMS      = 0,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'CPL_AROME',      &
  CSTDNAME   = '',               &
  CLONGNAME  = 'CPL_AROME',      &
  CUNITS     = '',               &
  CDIR       = '--',             &
  CCOMMENT   = 'Logical for AROME coupling file', &
  NGRID      = 0,                &
  NTYPE      = TYPELOG,          &
  NDIMS      = 0,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'COUPLING',       &
  CSTDNAME   = '',               &
  CLONGNAME  = 'COUPLING',       &
  CUNITS     = '',               &
  CDIR       = '--',             &
  CCOMMENT   = 'Logical for coupling file', &
  NGRID      = 0,                &
  NTYPE      = TYPELOG,          &
  NDIMS      = 0,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA(            &
  CMNHNAME   = 'UT',                        &
  CSTDNAME   = 'x_wind',                    &
  CLONGNAME  = 'UT',                        &
  CUNITS     = 'm s-1',                     &
  CDIR       = 'XY',                        &
  CCOMMENT   = 'X_Y_Z_U component of wind', &
  NGRID      = 2,                           &
  NTYPE      = TYPEREAL,                    &
  NDIMS      = 3,                           &
  LTIMEDEP   = .TRUE.                       ) )

call Add_field2list( TFIELDDATA(            &
  CMNHNAME   = 'VT',                        &
  CSTDNAME   = 'y_wind',                    &
  CLONGNAME  = 'VT',                        &
  CUNITS     = 'm s-1',                     &
  CDIR       = 'XY',                        &
  CCOMMENT   = 'X_Y_Z_V component of wind', &
  NGRID      = 3,                           &
  NTYPE      = TYPEREAL,                    &
  NDIMS      = 3,                           &
  LTIMEDEP   = .TRUE.                       ) )

call Add_field2list( TFIELDDATA(      &
  CMNHNAME   = 'WT',                  &
  CSTDNAME   = 'upward_air_velocity', &
  CLONGNAME  = 'WT',                  &
  CUNITS     = 'm s-1',               &
  CDIR       = 'XY',                  &
  CCOMMENT   = 'X_Y_Z_vertical wind', &
  NGRID      = 4,                     &
  NTYPE      = TYPEREAL,              &
  NDIMS      = 3,                     &
  LTIMEDEP   = .TRUE.                 ) )

call Add_field2list( TFIELDDATA(              &
  CMNHNAME   = 'THT',                         &
  CSTDNAME   = 'air_potential_temperature',   &
  CLONGNAME  = 'THT',                         &
  CUNITS     = 'K',                           &
  CDIR       = 'XY',                          &
  CCOMMENT   = 'X_Y_Z_potential temperature', &
  NGRID      = 1,                             &
  NTYPE      = TYPEREAL,                      &
  NDIMS      = 3,                             &
  LTIMEDEP   = .TRUE.                         ) )

call Add_field2list( TFIELDDATA(            &
  CMNHNAME   = 'UM',                        &
  CSTDNAME   = 'x_wind',                    &
  CLONGNAME  = 'UM',                        &
  CUNITS     = 'm s-1',                     &
  CDIR       = 'XY',                        &
  CCOMMENT   = 'X_Y_Z_U component of wind', &
  NGRID      = 2,                           &
  NTYPE      = TYPEREAL,                    &
  NDIMS      = 3,                           &
  LTIMEDEP   = .TRUE.                       ) )

call Add_field2list( TFIELDDATA(            &
  CMNHNAME   = 'VM',                        &
  CSTDNAME   = 'y_wind',                    &
  CLONGNAME  = 'VM',                        &
  CUNITS     = 'm s-1',                     &
  CDIR       = 'XY',                        &
  CCOMMENT   = 'X_Y_Z_V component of wind', &
  NGRID      = 3,                           &
  NTYPE      = TYPEREAL,                    &
  NDIMS      = 3,                           &
  LTIMEDEP   = .TRUE.                       ) )

call Add_field2list( TFIELDDATA(      &
  CMNHNAME   = 'WM',                  &
  CSTDNAME   = 'upward_air_velocity', &
  CLONGNAME  = 'WM',                  &
  CUNITS     = 'm s-1',               &
  CDIR       = 'XY',                  &
  CCOMMENT   = 'X_Y_Z_vertical wind', &
  NGRID      = 4,                     &
  NTYPE      = TYPEREAL,              &
  NDIMS      = 3,                     &
  LTIMEDEP   = .TRUE.                 ) )

call Add_field2list( TFIELDDATA(            &
  CMNHNAME   = 'DUM',                       &
  CSTDNAME   = 'x_wind',                    &
  CLONGNAME  = 'DUM',                       &
  CUNITS     = 'm s-1',                     &
  CDIR       = 'XY',                        &
  CCOMMENT   = 'X_Y_Z_U component of wind', &
  NGRID      = 2,                           &
  NTYPE      = TYPEREAL,                    &
  NDIMS      = 3,                           &
  LTIMEDEP   = .TRUE.                       ) )

call Add_field2list( TFIELDDATA(            &
  CMNHNAME   = 'DVM',                       &
  CSTDNAME   = 'y_wind',                    &
  CLONGNAME  = 'DVM',                       &
  CUNITS     = 'm s-1',                     &
  CDIR       = 'XY',                        &
  CCOMMENT   = 'X_Y_Z_V component of wind', &
  NGRID      = 3,                           &
  NTYPE      = TYPEREAL,                    &
  NDIMS      = 3,                           &
  LTIMEDEP   = .TRUE.                       ) )

call Add_field2list( TFIELDDATA(      &
  CMNHNAME   = 'DWM',                 &
  CSTDNAME   = 'upward_air_velocity', &
  CLONGNAME  = 'DWM',                 &
  CUNITS     = 'm s-1',               &
  CDIR       = 'XY',                  &
  CCOMMENT   = 'X_Y_Z_vertical wind', &
  NGRID      = 4,                     &
  NTYPE      = TYPEREAL,              &
  NDIMS      = 3,                     &
  LTIMEDEP   = .TRUE.                 ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'TKET',           &
  CSTDNAME   = '',               &
  CLONGNAME  = 'TKET',           &
  CUNITS     = 'm2 s-2',         &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Z_Turbulent Kinetic Energy', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'TKEMS',          &
  CSTDNAME   = '',               &
  CLONGNAME  = 'TKEMS',          &
  CUNITS     = 'm2 s-3',         &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Z_Turbulent Kinetic Energy adv source', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA(          &
  CMNHNAME   = 'PABST',                   &
  CSTDNAME   = 'air_pressure',            &
  CLONGNAME  = 'PABST',                   &
  CUNITS     = 'Pa',                      &
  CDIR       = 'XY',                      &
  CCOMMENT   = 'X_Y_Z_ABSolute Pressure', &
  NGRID      = 1,                         &
  NTYPE      = TYPEREAL,                  &
  NDIMS      = 3,                         &
  LTIMEDEP   = .TRUE.                     ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'PHIT',           &
  CSTDNAME   = '',               &
  CLONGNAME  = 'PHIT',           &
  CUNITS     = 'Pa',             &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Z_Reduced Pressure Oce/Shallow conv', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'RT',             &
  CSTDNAME   = '',               &
  CLONGNAME  = 'RT',             &
  CUNITS     = 'kg kg-1',        &
  CDIR       = 'XY',             &
  CCOMMENT   = 'Moist variables (rho Rn)', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 4,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA(           &
  CMNHNAME   = 'RVT',                      &
!TODO: check stdname
  CSTDNAME   = 'humidity_mixing_ratio',    &
  CLONGNAME  = 'RVT',                      &
  CUNITS     = 'kg kg-1',                  &
  CDIR       = 'XY',                       &
  CCOMMENT   = 'X_Y_Z_Vapor mixing Ratio', &
  NGRID      = 1,                          &
  NTYPE      = TYPEREAL,                   &
  NDIMS      = 3,                          &
  LTIMEDEP   = .TRUE.                      ) )

call Add_field2list( TFIELDDATA(           &
  CMNHNAME   = 'RCT',                      &
  CSTDNAME   = '',                         &
  CLONGNAME  = 'RCT',                      &
  CUNITS     = 'kg kg-1',                  &
  CDIR       = 'XY',                       &
  CCOMMENT   = 'X_Y_Z_Cloud mixing Ratio', &
  NGRID      = 1,                          &
  NTYPE      = TYPEREAL,                   &
  NDIMS      = 3,                          &
  LTIMEDEP   = .TRUE.                      ) )

call Add_field2list( TFIELDDATA(          &
  CMNHNAME   = 'RRT',                     &
  CSTDNAME   = '',                        &
  CLONGNAME  = 'RRT',                     &
  CUNITS     = 'kg kg-1',                 &
  CDIR       = 'XY',                      &
  CCOMMENT   = 'X_Y_Z_Rain mixing Ratio', &
  NGRID      = 1,                         &
  NTYPE      = TYPEREAL,                  &
  NDIMS      = 3,                         &
  LTIMEDEP   = .TRUE.                     ) )

call Add_field2list( TFIELDDATA(         &
  CMNHNAME   = 'RIT',                    &
!TODO: check stdname
  CSTDNAME   = 'cloud_ice_mixing_ratio', &
  CLONGNAME  = 'RIT',                    &
  CUNITS     = 'kg kg-1',                &
  CDIR       = 'XY',                     &
  CCOMMENT   = 'X_Y_Z_Ice mixing Ratio', &
  NGRID      = 1,                        &
  NTYPE      = TYPEREAL,                 &
  NDIMS      = 3,                        &
  LTIMEDEP   = .TRUE.                    ) )

call Add_field2list( TFIELDDATA(          &
  CMNHNAME   = 'RST',                     &
  CSTDNAME   = '',                        &
  CLONGNAME  = 'RST',                     &
  CUNITS     = 'kg kg-1',                 &
  CDIR       = 'XY',                      &
  CCOMMENT   = 'X_Y_Z_Snow mixing Ratio', &
  NGRID      = 1,                         &
  NTYPE      = TYPEREAL,                  &
  NDIMS      = 3,                         &
  LTIMEDEP   = .TRUE.                     ) )

call Add_field2list( TFIELDDATA(             &
  CMNHNAME   = 'RGT',                        &
  CSTDNAME   = '',                           &
  CLONGNAME  = 'RGT',                        &
  CUNITS     = 'kg kg-1',                    &
  CDIR       = 'XY',                         &
  CCOMMENT   = 'X_Y_Z_Graupel mixing Ratio', &
  NGRID      = 1,                            &
  NTYPE      = TYPEREAL,                     &
  NDIMS      = 3,                            &
  LTIMEDEP   = .TRUE.                        ) )

call Add_field2list( TFIELDDATA(          &
  CMNHNAME   = 'RHT',                     &
  CSTDNAME   = '',                        &
  CLONGNAME  = 'RHT',                     &
  CUNITS     = 'kg kg-1',                 &
  CDIR       = 'XY',                      &
  CCOMMENT   = 'X_Y_Z_Hail mixing Ratio', &
  NGRID      = 1,                         &
  NTYPE      = TYPEREAL,                  &
  NDIMS      = 3,                         &
  LTIMEDEP   = .TRUE.                     ) )

call Add_field2list( TFIELDDATA(        &
  CMNHNAME   = 'SUPSATMAX',             &
  CSTDNAME   = '',                      &
  CLONGNAME  = 'SUPSATMAX',             &
  CUNITS     = '',                      &
  CDIR       = 'XY',                    &
  CCOMMENT   = 'X_Y_Z_Supersaturation', &
  NGRID      = 1,                       &
  NTYPE      = TYPEREAL,                &
  NDIMS      = 3,                       &
  LTIMEDEP   = .TRUE.                   ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'NACT',           &
  CSTDNAME   = '',               &
  CLONGNAME  = 'NACT',           &
  CUNITS     = '',               &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Z_Nact',     &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'SSPRO',          &
  CSTDNAME   = '',               &
  CLONGNAME  = 'SSPRO',          &
  CUNITS     = '',               &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Z_Supersaturation', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'NPRO',           &
  CSTDNAME   = '',               &
  CLONGNAME  = 'NPRO',           &
  CUNITS     = '',               &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Z_NPRO',     &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'INPAP',          &
  CSTDNAME   = '',               &
  CLONGNAME  = 'INPAP',          &
  CUNITS     = 'kg m-2 s-1',     &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_INstantaneous Precipitating Aerosol Rate', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 2,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'ACPAP',          &
  CSTDNAME   = '',               &
  CLONGNAME  = 'ACPAP',          &
  CUNITS     = 'kg m-2',         &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_ACcumulated Precipitating Aerosol Rate', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 2,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'EFIELDU',        &
  CSTDNAME   = '',               &
  CLONGNAME  = 'EFIELDU',        &
  CUNITS     = 'V m-1',          &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Z_EFIELDU',  &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'EFIELDV',        &
  CSTDNAME   = '',               &
  CLONGNAME  = 'EFIELDV',        &
  CUNITS     = 'V m-1',          &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Z_EFIELDV',  &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'EFIELDW',        &
  CSTDNAME   = '',               &
  CLONGNAME  = 'EFIELDW',        &
  CUNITS     = 'V m-1',          &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Z_EFIELDW',  &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'NI_IAGGS',       &
  CSTDNAME   = '',               &
  CLONGNAME  = 'NI_IAGGS',       &
  CUNITS     = 'C m-3 s-1',      &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Z_NI_IAGGS', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'NI_IDRYG',       &
  CSTDNAME   = '',               &
  CLONGNAME  = 'NI_IDRYG',       &
  CUNITS     = 'C m-3 s-1',      &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Z_NI_IDRYG', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'NI_SDRYG',       &
  CSTDNAME   = '',               &
  CLONGNAME  = 'NI_SDRYG',       &
  CUNITS     = 'C m-3 s-1',      &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Z_NI_SDRYG', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'INDUC_CG',       &
  CSTDNAME   = '',               &
  CLONGNAME  = 'INDUC_CG',       &
  CUNITS     = 'C m-3 s-1',      &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Z_INDUC_CG', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'TRIG_IC',        &
  CSTDNAME   = '',               &
  CLONGNAME  = 'TRIG_IC',        &
  CUNITS     = '1',              &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_FLASH_MAP_TRIG_IC', &
  NGRID      = 1,                &
  NTYPE      = TYPEINT,          &
  NDIMS      = 2,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'IMPACT_CG',      &
  CSTDNAME   = '',               &
  CLONGNAME  = 'IMPACT_CG',      &
  CUNITS     = '1',              &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_FLASH_MAP_IMPACT_CG', &
  NGRID      = 1,                &
  NTYPE      = TYPEINT,          &
  NDIMS      = 2,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'AREA_CG',        &
  CSTDNAME   = '',               &
  CLONGNAME  = 'AREA_CG',        &
  CUNITS     = '1',              &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_FLASH_MAP_2DAREA_CG', &
  NGRID      = 1,                &
  NTYPE      = TYPEINT,          &
  NDIMS      = 2,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'AREA_IC',        &
  CSTDNAME   = '',               &
  CLONGNAME  = 'AREA_IC',        &
  CUNITS     = '1',              &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_FLASH_MAP_2DAREA_IC', &
  NGRID      = 1,                &
  NTYPE      = TYPEINT,          &
  NDIMS      = 2,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'FLASH_3DCG',     &
  CSTDNAME   = '',               &
  CLONGNAME  = 'FLASH_3DCG',     &
  CUNITS     = '1',              &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Z_FLASH_MAP_3DCG', &
  NGRID      = 1,                &
  NTYPE      = TYPEINT,          &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'FLASH_3DIC',     &
  CSTDNAME   = '',               &
  CLONGNAME  = 'FLASH_3DIC',     &
  CUNITS     = '1',              &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Z_FLASH_MAP_3DIC', &
  NGRID      = 1,                &
  NTYPE      = TYPEINT,          &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'PHC',            &
  CSTDNAME   = '',               &
  CLONGNAME  = 'PHC',            &
  CUNITS     = '1',              &
  CDIR       = 'XY',             &
  CCOMMENT   = 'pH in cloud',    &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'PHR',            &
  CSTDNAME   = '',               &
  CLONGNAME  = 'PHR',            &
  CUNITS     = '1',              &
  CDIR       = 'XY',             &
  CCOMMENT   = 'pH in rain',     &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'LSUM',           &
  CSTDNAME   = '',               &
  CLONGNAME  = 'LSUM',           &
  CUNITS     = 'm s-1',          &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Z_Large Scale U component', &
  NGRID      = 2,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'LSVM',           &
  CSTDNAME   = '',               &
  CLONGNAME  = 'LSVM',           &
  CUNITS     = 'm s-1',          &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Z_Large Scale V component', &
  NGRID      = 3,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'LSWM',           &
  CSTDNAME   = '',               &
  CLONGNAME  = 'LSWM',           &
  CUNITS     = 'm s-1',          &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Z_Large Scale vertical wind', &
  NGRID      = 4,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'LSTHM',          &
  CSTDNAME   = '',               &
  CLONGNAME  = 'LSTHM',          &
  CUNITS     = 'K',              &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Z_Large Scale potential Temperature', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'LSRVM',          &
  CSTDNAME   = '',               &
  CLONGNAME  = 'LSRVM',          &
  CUNITS     = 'kg kg-1',        &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Z_Large Scale Vapor Mixing Ratio', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'RIMX',           &
  CSTDNAME   = '',               &
  CLONGNAME  = 'RIMX',           &
  CUNITS     = '',               &
  CDIR       = '--',             &
  CCOMMENT   = 'Number of points in the lateral absorbing layer in the x direction', &
  NGRID      = 1,                &
  NTYPE      = TYPEINT,          &
  NDIMS      = 0,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'RIMY',           &
  CSTDNAME   = '',               &
  CLONGNAME  = 'RIMY',           &
  CUNITS     = '',               &
  CDIR       = '--',             &
  CCOMMENT   = 'Number of points in the lateral absorbing layer in the y direction', &
  NGRID      = 1,                &
  NTYPE      = TYPEINT,          &
  NDIMS      = 0,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'HORELAX_UVWTH',  &
  CSTDNAME   = '',               &
  CLONGNAME  = 'HORELAX_UVWTH',  &
  CUNITS     = '',               &
  CDIR       = '--',             &
  CCOMMENT   = 'Switch to activate the HOrizontal RELAXation', &
  NGRID      = 1,                &
  NTYPE      = TYPELOG,          &
  NDIMS      = 0,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'LBXUM',          &
  CSTDNAME   = '',               &
  CLONGNAME  = 'LBXUM',          &
  CUNITS     = 'm s-1',          &
!   CDIR       = ''
  CLBTYPE    = 'LBXU',           &
  CCOMMENT   = '2_Y_Z_LBXUM',    &
  NGRID      = 2,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'LBXVM',          &
  CSTDNAME   = '',               &
  CLONGNAME  = 'LBXVM',          &
  CUNITS     = 'm s-1',          &
!   CDIR       = ''
  CLBTYPE    = 'LBX',            &
  CCOMMENT   = '2_Y_Z_LBXVM',    &
  NGRID      = 3,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'LBXWM',          &
  CSTDNAME   = '',               &
  CLONGNAME  = 'LBXWM',          &
  CUNITS     = 'm s-1',          &
!   CDIR       = ''
  CLBTYPE    = 'LBX',            &
  CCOMMENT   = '2_Y_Z_LBXWM',    &
  NGRID      = 4,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'LBYUM',          &
  CSTDNAME   = '',               &
  CLONGNAME  = 'LBYUM',          &
  CUNITS     = 'm s-1',          &
!   CDIR       = ''
  CLBTYPE    = 'LBY',            &
  CCOMMENT   = '2_Y_Z_LBYUM',    &
  NGRID      = 2,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'LBYVM',          &
  CSTDNAME   = '',               &
  CLONGNAME  = 'LBYVM',          &
  CUNITS     = 'm s-1',          &
!   CDIR       = ''
  CLBTYPE    = 'LBYV',           &
  CCOMMENT   = '2_Y_Z_LBYVM',    &
  NGRID      = 3,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'LBYWM',          &
  CSTDNAME   = '',               &
  CLONGNAME  = 'LBYWM',          &
  CUNITS     = 'm s-1',          &
!   CDIR       = ''
  CLBTYPE    = 'LBY',            &
  CCOMMENT   = '2_Y_Z_LBYWM',    &
  NGRID      = 4,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'LBXTHM',         &
  CSTDNAME   = '',               &
  CLONGNAME  = 'LBXTHM',         &
  CUNITS     = 'K',              &
!   CDIR       = ''
  CLBTYPE    = 'LBX',            &
  CCOMMENT   = '2_Y_Z_LBXTHM',   &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'LBYTHM',         &
  CSTDNAME   = '',               &
  CLONGNAME  = 'LBYTHM',         &
  CUNITS     = 'K',              &
!   CDIR       = ''
  CLBTYPE    = 'LBY',            &
  CCOMMENT   = '2_Y_Z_LBYTHM',   &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'HORELAX_TKE',    &
  CSTDNAME   = '',               &
  CLONGNAME  = 'HORELAX_TKE',    &
  CUNITS     = '',               &
  CDIR       = '--',             &
  CCOMMENT   = 'Switch to activate the HOrizontal RELAXation', &
  NGRID      = 1,                &
  NTYPE      = TYPELOG,          &
  NDIMS      = 0,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'LBXTKEM',        &
  CSTDNAME   = '',               &
  CLONGNAME  = 'LBXTKEM',        &
  CUNITS     = 'm2 s-2',         &
!   CDIR       = ''
  CLBTYPE    = 'LBX',            &
  CCOMMENT   = '2_Y_Z_LBXTKEM',  &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'LBYTKEM',        &
  CSTDNAME   = '',               &
  CLONGNAME  = 'LBYTKEM',        &
  CUNITS     = 'm2 s-2',         &
!   CDIR       = ''
  CLBTYPE    = 'LBY',            &
  CCOMMENT   = '2_Y_Z_LBYTKEM',  &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'DRYMASST',       &
  CSTDNAME   = '',               &
  CLONGNAME  = 'DRYMASST',       &
  CUNITS     = 'kg',             &
  CDIR       = '--',             &
  CCOMMENT   = 'Total Dry Mass', &
  NGRID      = 0,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 0,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA(        &
  CMNHNAME   = 'DRYMASSS',              &
  CSTDNAME   = '',                      &
  CLONGNAME  = 'DRYMASSS',              &
  CUNITS     = 'kg',                    &
  CDIR       = '--',                    &
  CCOMMENT   = 'Total Dry Mass Source', &
  NGRID      = 0,                       &
  NTYPE      = TYPEREAL,                &
  NDIMS      = 0,                       &
  LTIMEDEP   = .TRUE.                   ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'BL_DEPTH',       &
  CSTDNAME   = '',               &
  CLONGNAME  = 'BL_DEPTH',       &
  CUNITS     = 'm',              &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_BL_DEPTH',   &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 2,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'SBL_DEPTH',      &
  CSTDNAME   = '',               &
  CLONGNAME  = 'SBL_DEPTH',      &
  CUNITS     = 'm',              &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_BL_SDEPTH',  &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 2,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'WTHVMF',         &
  CSTDNAME   = '',               &
  CLONGNAME  = 'WTHVMF',         &
  CUNITS     = 'm K s-1',        &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_WTHVMF',     &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'WUMF',           &
  CSTDNAME   = '',               &
  CLONGNAME  = 'WUMF',           &
  CUNITS     = 'm K s-1',        &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_WUMF',       &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'WVMF',           &
  CSTDNAME   = '',               &
  CLONGNAME  = 'WVMF',           &
  CUNITS     = 'm K s-1',        &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_WVMF',       &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'SRCT',           &
  CSTDNAME   = '',               &
  CLONGNAME  = 'SRCT',           &
  CUNITS     = 'kg kg-2',        &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Z_normalized 2nd_order moment s_r_c/2Sigma_s2', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'SIGS',           &
  CSTDNAME   = '',               &
  CLONGNAME  = 'SIGS',           &
  CUNITS     = 'kg kg-2',        &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Z_Sigma_s from turbulence scheme', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'RHOREFZ',        &
  CSTDNAME   = '',               &
  CLONGNAME  = 'RHOREFZ',        &
  CUNITS     = 'kg m-3',         &
  CDIR       = 'ZZ',             &
  CCOMMENT   = 'rhodz for reference state without orography', &
  NGRID      = 4,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 1,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'THVREFZ',        &
  CSTDNAME   = '',               &
  CLONGNAME  = 'THVREFZ',        &
  CUNITS     = 'K',              &
  CDIR       = 'ZZ',             &
  CCOMMENT   = 'thetavz for reference state without orography', &
  NGRID      = 4,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 1,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'EXNTOP',         &
  CSTDNAME   = 'dimensionless_exner_function', &
  CLONGNAME  = 'EXNTOP',         &
  CUNITS     = '',               &
  CDIR       = '--',             &
  CCOMMENT   = 'Exner function at model top', &
  NGRID      = 4,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 0,                &
  LTIMEDEP   = .FALSE.           ) )

  IF (TRIM(CPROGRAM) == 'MESONH' .OR. TRIM(CPROGRAM) == 'DIAG' .OR. TRIM(CPROGRAM) == 'LFICDF') THEN

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'US_PRES',        &
  CSTDNAME   = '',               &
  CLONGNAME  = 'US_PRES',        &
!TODO: units?
  CUNITS     = '',               &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Z_US_PRES',  &
  NGRID      = 2,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'VS_PRES',        &
  CSTDNAME   = '',               &
  CLONGNAME  = 'VS_PRES',        &
!TODO: units?
  CUNITS     = '',               &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Z_VS_PRES',  &
  NGRID      = 3,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'WS_PRES',        &
  CSTDNAME   = '',               &
  CLONGNAME  = 'WS_PRES',        &
!TODO: units?
  CUNITS     = '',               &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Z_WS_PRES',  &
  NGRID      = 4,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'THS_CLD',        &
  CSTDNAME   = '',               &
  CLONGNAME  = 'THS_CLD',        &
!TODO: units?
  CUNITS     = '',               &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Z_THS_CLD',  &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'RS_CLD',         &
  CSTDNAME   = '',               &
  CLONGNAME  = 'RS_CLD',         &
  CUNITS     = '',               &
  CDIR       = 'XY',             &
  CCOMMENT   = 'Source of Moist variables', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 4,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'RVS_CLD',        &
  CSTDNAME   = '',               &
  CLONGNAME  = 'RVS_CLD',        &
  CUNITS     = '',               &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Z_RVS_CLD',  &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'RCS_CLD',        &
  CSTDNAME   = '',               &
  CLONGNAME  = 'RCS_CLD',        &
  CUNITS     = '',               &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Z_RCS_CLD',  &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'RRS_CLD',        &
  CSTDNAME   = '',               &
  CLONGNAME  = 'RRS_CLD',        &
  CUNITS     = '',               &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Z_RRS_CLD',  &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'RIS_CLD',        &
  CSTDNAME   = '',               &
  CLONGNAME  = 'RIS_CLD',        &
  CUNITS     = '',               &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Z_RIS_CLD',  &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'RSS_CLD',        &
  CSTDNAME   = '',               &
  CLONGNAME  = 'RSS_CLD',        &
  CUNITS     = '',               &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Z_RSS_CLD',  &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'RGS_CLD',        &
  CSTDNAME   = '',               &
  CLONGNAME  = 'RGS_CLD',        &
  CUNITS     = '',               &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Z_RGS_CLD',  &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'RHS_CLD',        &
  CSTDNAME   = '',               &
  CLONGNAME  = 'RHS_CLD',        &
  CUNITS     = '',               &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Z_RHS_CLD',  &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'CLDFR',          &
  CSTDNAME   = '',               &
  CLONGNAME  = 'CLDFR',          &
  CUNITS     = '1',              &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Z_CLouD FRaction', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'ICEFR',          &
  CSTDNAME   = '',               &
  CLONGNAME  = 'ICEFR',          &
  CUNITS     = '1',              &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Z_ICE cloud FRaction', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'CIT',            &
  CSTDNAME   = '',               &
  CLONGNAME  = 'CIT',            &
  CUNITS     = 'm-3',            &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Z_Cloud Ice concentration', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA(      &
  CMNHNAME   = 'RAINFR',              &
  CSTDNAME   = '',                    &
  CLONGNAME  = 'RAINFR',              &
  CUNITS     = '1',                   &
  CDIR       = 'XY',                  &
  CCOMMENT   = 'X_Y_Z_Rain FRaction', &
  NGRID      = 1,                     &
  NTYPE      = TYPEREAL,              &
  NDIMS      = 3,                     &
  LTIMEDEP   = .TRUE.                 ) )
!
END IF ! CPROGRAM=MESONH .OR. DIAG .OR. LFICDF
!
!
call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'RHODREF',        &
  CSTDNAME   = '',               &
  CLONGNAME  = 'RHODREF',        &
  CUNITS     = 'kg m-3',         &
  CDIR       = 'XY',             &
  CCOMMENT   = 'Dry density for reference state with orography', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'THVREF',         &
  CSTDNAME   = '',               &
  CLONGNAME  = 'THVREF',         &
  CUNITS     = 'K',              &
  CDIR       = 'XY',             &
  CCOMMENT   = 'Thetav for reference state with orography', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .FALSE.           ) )
!
!
IF (     TRIM(CPROGRAM) =='MESONH' .OR. TRIM(CPROGRAM) == 'DIAG'  &
    .OR. TRIM(CPROGRAM) == 'LFICDF'.OR. TRIM(CPROGRAM) == 'SPAWN' ) THEN
!
call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'DTHRAD',         &
  CSTDNAME   = '',               &
  CLONGNAME  = 'DTHRAD',         &
  CUNITS     = 'K s-1',          &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Z_RADiative heating/cooling rate', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'FLALWD',         &
  CSTDNAME   = '',               &
  CLONGNAME  = 'FLALWD',         &
  CUNITS     = 'W m-2',          &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Downward Long Waves on FLAT surface', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 2,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'DIRFLASWD',      &
  CSTDNAME   = '',               &
  CLONGNAME  = 'DIRFLASWD',      &
  CUNITS     = 'W m-2',          &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_DIRect Downward Short Waves on FLAT surface', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  NDIMLIST   = [ NMNHDIM_NI, NMNHDIM_NJ, NMNHDIM_NSWB ], &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'SCAFLASWD',      &
  CSTDNAME   = '',               &
  CLONGNAME  = 'SCAFLASWD',      &
  CUNITS     = 'W m-2',          &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_SCAttered Downward Short Waves on FLAT surface', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  NDIMLIST   = [ NMNHDIM_NI, NMNHDIM_NJ, NMNHDIM_NSWB ], &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'DIRSRFSWD',      &
  CSTDNAME   = '',               &
  CLONGNAME  = 'DIRSRFSWD',      &
  CUNITS     = 'W m-2',          &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_DIRect Downward Short Waves', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  NDIMLIST   = [ NMNHDIM_NI, NMNHDIM_NJ, NMNHDIM_NSWB ], &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'CLEARCOL_TM1',   &
  CSTDNAME   = '',               &
  CLONGNAME  = 'CLEARCOL_TM1',   &
  CUNITS     = '',               &
  CDIR       = 'XY',             &
  CCOMMENT   = 'TRACE OF CLOUD', &
  NGRID      = 1,                &
  NTYPE      = TYPEINT,          &
  NDIMS      = 2,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'ZENITH',         &
  CSTDNAME   = 'zenith_angle',   &
  CLONGNAME  = 'ZENITH',         &
  CUNITS     = 'rad',            &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_ZENITH',     &
  NGRID      = 4,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 2,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'AZIM',           &
  CSTDNAME   = '',               &
  CLONGNAME  = 'AZIM',           &
  CUNITS     = 'rad',            &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_AZIMuth',    &
  NGRID      = 4,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 2,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA(    &
  CMNHNAME   = 'DIR_ALB',           &
  CSTDNAME   = '',                  &
  CLONGNAME  = 'DIR_ALB',           &
  CUNITS     = '1',                 &
  CDIR       = 'XY',                &
  CCOMMENT   = 'X_Y_DIRect ALBedo', &
  NGRID      = 1,                   &
  NTYPE      = TYPEREAL,            &
  NDIMS      = 3,                   &
  NDIMLIST   = [ NMNHDIM_NI, NMNHDIM_NJ, NMNHDIM_NSWB ], &
  LTIMEDEP   = .TRUE.               ) )

call Add_field2list( TFIELDDATA(       &
  CMNHNAME   = 'SCA_ALB',              &
  CSTDNAME   = '',                     &
  CLONGNAME  = 'SCA_ALB',              &
  CUNITS     = '1',                    &
  CDIR       = 'XY',                   &
  CCOMMENT   = 'X_Y_SCAttered ALBedo', &
  NGRID      = 1,                      &
  NTYPE      = TYPEREAL,               &
  NDIMS      = 3,                      &
  NDIMLIST   = [ NMNHDIM_NI, NMNHDIM_NJ, NMNHDIM_NSWB ], &
  LTIMEDEP   = .TRUE.                  ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'EMIS',           &
  CSTDNAME   = '',               &
  CLONGNAME  = 'EMIS',           &
  CUNITS     = '1',              &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_EMISsivity', &
  NGRID      = 4,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  NDIMLIST   = [ NMNHDIM_NI, NMNHDIM_NJ, NMNHDIM_NLWB ], &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'TSRAD',          &
  CSTDNAME   = '',               &
  CLONGNAME  = 'TSRAD',          &
  CUNITS     = 'K',              &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_RADiative Surface Temperature', &
  NGRID      = 4,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 2,                &
  LTIMEDEP   = .TRUE.            ) )

END IF !CPROGRAM=MESONH .OR. DIAG .OR. LFICDF .OR. SPAWN

IF (     TRIM(CPROGRAM) == 'MESONH' .OR. TRIM(CPROGRAM) == 'DIAG'  .OR. TRIM(CPROGRAM) == 'REAL' &
    .OR. TRIM(CPROGRAM) == 'LFICDF' .OR. TRIM(CPROGRAM) == 'SPAWN'                               ) THEN
!
! Blaze fire model fields
!
call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'FMREFINRATIOX',  &
  CLONGNAME  = 'FMREFINRATIOX',  &
  CSTDNAME   = '',               &
  CUNITS     = '1',              &
  CDIR       = '--',             &
  CCOMMENT   = 'Blaze fire model: grid refinement ratio (x direction)', &
  NGRID      = 0,                &
  NTYPE      = TYPEINT,          &
  NDIMS      = 0,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'FMREFINRATIOY',  &
  CLONGNAME  = 'FMREFINRATIOY',  &
  CSTDNAME   = '',               &
  CUNITS     = '1',              &
  CDIR       = '--',             &
  CCOMMENT   = 'Blaze fire model: grid refinement ratio (y direction)', &
  NGRID      = 0,                &
  NTYPE      = TYPEINT,          &
  NDIMS      = 0,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA(     &
  CMNHNAME   = 'FMPHI',              &
  CSTDNAME   = '',                   &
  CLONGNAME  = 'level set function', &
  CUNITS     = '',                   &
  CDIR       = 'XY',                 &
  CCOMMENT   = 'X_Y_F Blaze fire model level set function', &
  NGRID      = 1,                    &
  NTYPE      = TYPEREAL,             &
  NDIMS      = 3,                    &
  NDIMLIST   = [ NMNHDIM_NI, NMNHDIM_NJ, NMNHDIM_NOTLISTED ], &
  LTIMEDEP   = .TRUE.                ) )

call Add_field2list( TFIELDDATA(   &
  CMNHNAME   = 'FMBMAP',           &
  CSTDNAME   = '',                 &
  CLONGNAME  = 'fire burning map', &
  CUNITS     = 's',                &
  CDIR       = 'XY',               &
  CCOMMENT   = 'X_Y_F Blaze fire model burning map, i.e. arrival time matrix', &
  NGRID      = 1,                  &
  NTYPE      = TYPEREAL,           &
  NDIMS      = 3,                  &
  NDIMLIST   = [ NMNHDIM_NI, NMNHDIM_NJ, NMNHDIM_NOTLISTED ], &
  LTIMEDEP   = .TRUE.              ) )

call Add_field2list( TFIELDDATA(            &
  CMNHNAME   = 'FMASE',                     &
  CSTDNAME   = '',                          &
  CLONGNAME  = 'available sensible energy', &
  CUNITS     = 'kJ m-2',                    &
  CDIR       = 'XY',                        &
  CCOMMENT   = 'X_Y_F Blaze fire model available sensible energy of vegetation', &
  NGRID      = 1,                           &
  NTYPE      = TYPEREAL,                    &
  NDIMS      = 3,                           &
  NDIMLIST   = [ NMNHDIM_NI, NMNHDIM_NJ, NMNHDIM_NOTLISTED ], &
  LTIMEDEP   = .TRUE.                       ) )

call Add_field2list( TFIELDDATA(          &
  CMNHNAME   = 'FMAWC',                   &
  CSTDNAME   = '',                        &
  CLONGNAME  = 'available water content', &
  CUNITS     = 'kg m-2',                  &
  CDIR       = 'XY',                      &
  CCOMMENT   = 'X_Y_F Blaze fire model available liquid water of vegetation', &
  NGRID      = 1,                         &
  NTYPE      = TYPEREAL,                  &
  NDIMS      = 3,                         &
  NDIMLIST   = [ NMNHDIM_NI, NMNHDIM_NJ, NMNHDIM_NOTLISTED ], &
  LTIMEDEP   = .TRUE.                     ) )

call Add_field2list( TFIELDDATA(             &
  CMNHNAME   = 'FMWINDU',                    &
  CSTDNAME   = '',                           &
  CLONGNAME  = 'fire model filtered wind u', &
  CUNITS     = 'm s-1',                      &
  CDIR       = 'XY',                         &
  CCOMMENT   = 'X_Y_F Blaze fire model EWAM filtered u wind', &
  NGRID      = 1,                            &
  NTYPE      = TYPEREAL,                     &
  NDIMS      = 3,                            &
  NDIMLIST   = [ NMNHDIM_NI, NMNHDIM_NJ, NMNHDIM_NOTLISTED ], &
  LTIMEDEP   = .TRUE.                        ) )

call Add_field2list( TFIELDDATA(             &
  CMNHNAME   = 'FMWINDV',                    &
  CSTDNAME   = '',                           &
  CLONGNAME  = 'fire model filtered wind v', &
  CUNITS     = 'm s-1',                      &
  CDIR       = 'XY',                         &
  CCOMMENT   = 'X_Y_F Blaze fire model EWAM filtered v wind', &
  NGRID      = 1,                            &
  NTYPE      = TYPEREAL,                     &
  NDIMS      = 3,                            &
  NDIMLIST   = [ NMNHDIM_NI, NMNHDIM_NJ, NMNHDIM_NOTLISTED ], &
  LTIMEDEP   = .TRUE.                        ) )

call Add_field2list( TFIELDDATA(             &
  CMNHNAME   = 'FMWINDW',                    &
  CSTDNAME   = '',                           &
  CLONGNAME  = 'fire model filtered wind w', &
  CUNITS     = 'm s-1',                      &
  CDIR       = 'XY',                         &
  CCOMMENT   = 'X_Y_F Blaze fire model EWAM filtered w wind', &
  NGRID      = 1,                            &
  NTYPE      = TYPEREAL,                     &
  NDIMS      = 3,                            &
  NDIMLIST   = [ NMNHDIM_NI, NMNHDIM_NJ, NMNHDIM_NOTLISTED ], &
  LTIMEDEP   = .TRUE.                        ) )

call Add_field2list( TFIELDDATA(                 &
  CMNHNAME   = 'FMHWS',                          &
  CSTDNAME   = '',                               &
  CLONGNAME  = 'filtered horizontal wind speed', &
  CUNITS     = 'm s-1',                          &
  CDIR       = 'XY',                             &
  CCOMMENT   = 'X_Y_F Blaze filtered horizontal wind speed', &
  NGRID      = 1,                                &
  NTYPE      = TYPEREAL,                         &
  NDIMS      = 3,                                &
  NDIMLIST   = [ NMNHDIM_NI, NMNHDIM_NJ, NMNHDIM_NOTLISTED ], &
  LTIMEDEP   = .TRUE.                            ) )

call Add_field2list( TFIELDDATA(      &
  CMNHNAME   = 'FMROS',               &
  CSTDNAME   = '',                    &
  CLONGNAME  = 'fire rate of spread', &
  CUNITS     = 'm s-1',               &
  CDIR       = 'XY',                  &
  CCOMMENT   = 'X_Y_F Blaze fire model rate of spread', &
  NGRID      = 1,                     &
  NTYPE      = TYPEREAL,              &
  NDIMS      = 3,                     &
  NDIMLIST   = [ NMNHDIM_NI, NMNHDIM_NJ, NMNHDIM_NOTLISTED ], &
  LTIMEDEP   = .TRUE.                 ) )

call Add_field2list( TFIELDDATA(              &
  CMNHNAME   = 'FMROS0',                      &
  CSTDNAME   = '',                            &
  CLONGNAME  = 'fire rate of spread no wind', &
  CUNITS     = 'm s-1',                       &
  CDIR       = 'XY',                          &
  CCOMMENT   = 'X_Y_F Blaze fire model rate of spread without wind and slope', &
  NGRID      = 1,                             &
  NTYPE      = TYPEREAL,                      &
  NDIMS      = 3,                             &
  NDIMLIST   = [ NMNHDIM_NI, NMNHDIM_NJ, NMNHDIM_NOTLISTED ], &
  LTIMEDEP   = .TRUE.                         ) )

call Add_field2list( TFIELDDATA(          &
  CMNHNAME   = 'FMFLUXHDH',               &
  CSTDNAME   = '',                        &
  CLONGNAME  = 'fire sensible heat flux', &
  CUNITS     = 'W m-2',                   &
  CDIR       = 'XY',                      &
  CCOMMENT   = 'X_Y_F Blaze fire model sensible heat flux', &
  NGRID      = 1,                         &
  NTYPE      = TYPEREAL,                  &
  NDIMS      = 3,                         &
  NDIMLIST   = [ NMNHDIM_NI, NMNHDIM_NJ, NMNHDIM_NOTLISTED ], &
  LTIMEDEP   = .TRUE.                     ) )

call Add_field2list( TFIELDDATA(        &
  CMNHNAME   = 'FMFLUXHDW',             &
  CSTDNAME   = '',                      &
  CLONGNAME  = 'fire latent heat flux', &
  CUNITS     = 'kg m-2 s-1',            &
  CDIR       = 'XY',                    &
  CCOMMENT   = 'X_Y_F Blaze fire model latent heat flux', &
  NGRID      = 1,                       &
  NTYPE      = TYPEREAL,                &
  NDIMS      = 3,                       &
  NDIMLIST   = [ NMNHDIM_NI, NMNHDIM_NJ, NMNHDIM_NOTLISTED ], &
  LTIMEDEP   = .TRUE.                   ) )

call Add_field2list( TFIELDDATA(        &
  CMNHNAME   = 'FMGRADOROX',            &
  CSTDNAME   = '',                      &
  CLONGNAME  = 'orographic x-gradient', &
  CUNITS     = '',                      &
  CDIR       = 'XY',                    &
  CCOMMENT   = 'X_Y_F Blaze fire model orographic gradient on x direction on fire mesh', &
  NGRID      = 1,                       &
  NTYPE      = TYPEREAL,                &
  NDIMS      = 3,                       &
  NDIMLIST   = [ NMNHDIM_NI, NMNHDIM_NJ, NMNHDIM_NOTLISTED ], &
  LTIMEDEP   = .TRUE.                   ) )

call Add_field2list( TFIELDDATA(        &
  CMNHNAME   = 'FMGRADOROY',            &
  CSTDNAME   = '',                      &
  CLONGNAME  = 'orographic y-gradient', &
  CUNITS     = '',                      &
  CDIR       = 'XY',                    &
  CCOMMENT   = 'X_Y_F Blaze fire model orographic gradient on y direction on fire mesh', &
  NGRID      = 1,                       &
  NTYPE      = TYPEREAL,                &
  NDIMS      = 3,                       &
  NDIMLIST   = [ NMNHDIM_NI, NMNHDIM_NJ, NMNHDIM_NOTLISTED ], &
  LTIMEDEP   = .TRUE.                   ) )
!
! end of Blaze fields
!
END IF !CPROGRAM=MESONH .OR. DIAG .OR. LFICDF .OR. SPAWN .OR. REAL
!
!
IF ( TRIM(CPROGRAM) /= 'PGD' .AND. TRIM(CPROGRAM) /= 'NESPGD' ) THEN
!
call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'COUNTCONV',      &
  CSTDNAME   = '',               &
  CLONGNAME  = 'COUNTCONV',      &
  CUNITS     = '',               &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_COUNTCONV',  &
  NGRID      = 1,                &
  NTYPE      = TYPEINT,          &
  NDIMS      = 2,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'DTHCONV',        &
  CSTDNAME   = '',               &
  CLONGNAME  = 'DTHCONV',        &
  CUNITS     = 'K s-1',          &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Z_CONVective heating/cooling rate', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'DRVCONV',        &
  CSTDNAME   = '',               &
  CLONGNAME  = 'DRVCONV',        &
  CUNITS     = 's-1',            &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Z_CONVective R_v tendency', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'DRCCONV',        &
  CSTDNAME   = '',               &
  CLONGNAME  = 'DRCCONV',        &
  CUNITS     = 's-1',            &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Z_CONVective R_c tendency', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'DRICONV',        &
  CSTDNAME   = '',               &
  CLONGNAME  = 'DRICONV',        &
  CUNITS     = 's-1',            &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Z_CONVective R_i tendency', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'PRCONV',         &
  CSTDNAME   = '',               &
  CLONGNAME  = 'PRCONV',         &
  CUNITS     = 'm s-1',          &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_CONVective instantaneous Precipitation Rate', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 2,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'PACCONV',        &
  CSTDNAME   = '',               &
  CLONGNAME  = 'PACCONV',        &
  CUNITS     = 'm',              &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_CONVective ACcumulated Precipitation rate', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 2,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'PRSCONV',        &
  CSTDNAME   = '',               &
  CLONGNAME  = 'PRSCONV',        &
  CUNITS     = 'm s-1',          &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_CONVective instantaneous Precipitation Rate for Snow', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 2,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA(    &
  CMNHNAME   = 'DSVCONV',           &
  CSTDNAME   = '',                  &
  CLONGNAME  = 'DSVCONV',           &
  CUNITS     = 's-1',               &
  CDIR       = 'XY',                &
  CCOMMENT   = 'Tracer tendencies', &
  NGRID      = 1,                   &
  NTYPE      = TYPEREAL,            &
  NDIMS      = 4,                   &
  LTIMEDEP   = .TRUE.               ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'PRLFLXCONV',     &
  CSTDNAME   = '',               &
  CLONGNAME  = 'PRLFLXCONV',     &
  CUNITS     = 'm s-1',          &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Liquid Precipitation Convective Flux', &
  NGRID      = 4,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'PRSFLXCONV',     &
  CSTDNAME   = '',               &
  CLONGNAME  = 'PRSFLXCONV',     &
  CUNITS     = 'm s-1',          &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Solid Precipitation Convective Flux', &
  NGRID      = 4,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'UMFCONV',        &
  CSTDNAME   = '',               &
  CLONGNAME  = 'UMFCONV',        &
  CUNITS     = 'kg s-1 m-2',     &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Updraft Convective Mass Flux', &
  NGRID      = 4,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'DMFCONV',        &
  CSTDNAME   = '',               &
  CLONGNAME  = 'DMFCONV',        &
  CUNITS     = 'kg s-1 m-2',     &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Downdraft Convective Mass Flux', &
  NGRID      = 4,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'MFCONV',         &
  CSTDNAME   = '',               &
  CLONGNAME  = 'MFCONV',         &
  CUNITS     = 'kg s-1 m-2',     &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Convective Mass Flux', &
  NGRID      = 4,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'CAPE',           &
  CSTDNAME   = '',               &
  CLONGNAME  = 'CAPE',           &
  CUNITS     = 'J kg-1',         &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Convective Available Potentiel Energy', &
  NGRID      = 4,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 2,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'CLTOPCONV_LVL',  &
  CSTDNAME   = '',               &
  CLONGNAME  = 'CLTOPCONV_LVL',  &
  CUNITS     = '1',              &
  CDIR       = 'XY',             &
  CCOMMENT   = 'Convective cloud top level', &
  NGRID      = 1,                &
  NTYPE      = TYPEINT,          &
  NDIMS      = 2,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'CLBASCONV_LVL',  &
  CSTDNAME   = '',               &
  CLONGNAME  = 'CLBASCONV_LVL',  &
  CUNITS     = '1',              &
  CDIR       = 'XY',             &
  CCOMMENT   = 'Convective cloud base level', &
  NGRID      = 1,                &
  NTYPE      = TYPEINT,          &
  NDIMS      = 2,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'IC_RATE',        &
  CSTDNAME   = '',               &
  CLONGNAME  = 'IC_RATE',        &
  CUNITS     = 's-1',            &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_IntraCloud lightning Rate', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 2,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'CG_RATE',        &
  CSTDNAME   = '',               &
  CLONGNAME  = 'CG_RATE',        &
  CUNITS     = 's-1',            &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_CloudGround lightning Rate', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 2,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'IC_TOTAL_NB',    &
  CSTDNAME   = '',               &
  CLONGNAME  = 'IC_TOTAL_NB',    &
  CUNITS     = '1',              &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_IntraCloud lightning Number', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 2,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'CG_TOTAL_NB',    &
  CSTDNAME   = '',               &
  CLONGNAME  = 'CG_TOTAL_NB',    &
  CUNITS     = '1',              &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_CloudGround lightning Number', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 2,                &
  LTIMEDEP   = .TRUE.            ) )
!
END IF !CPROGRAM/=PGD , NESPGD
!
!
call Add_field2list( TFIELDDATA(     &
  CMNHNAME   = 'SSO_ANIS',           &
  CSTDNAME   = '',                   &
  CLONGNAME  = 'SSO_ANIS',           &
  CUNITS     = 'm',                  &
  CDIR       = 'XY',                 &
  CCOMMENT   = 'X_Y_SSO_ANISOTROPY', &
  NGRID      = 4,                    &
  NTYPE      = TYPEREAL,             &
  NDIMS      = 2,                    &
  LTIMEDEP   = .FALSE.               ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'SSO_SLOPE',      &
  CSTDNAME   = '',               &
  CLONGNAME  = 'SSO_SLOPE',      &
  CUNITS     = '',               &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_SSO_SLOPE',  &
  NGRID      = 4,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 2,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'SSO_DIR',        &
  CSTDNAME   = '',               &
  CLONGNAME  = 'SSO_DIR',        &
  CUNITS     = 'degree',         &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_SSO_DIR',    &
  NGRID      = 4,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 2,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'AVG_ZS',         &
  CSTDNAME   = '',               &
  CLONGNAME  = 'AVG_ZS',         &
  CUNITS     = 'm',              &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_AVG_ZS',     &
  NGRID      = 4,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 2,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'SIL_ZS',         &
  CSTDNAME   = '',               &
  CLONGNAME  = 'SIL_ZS',         &
  CUNITS     = 'm',              &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_SIL_ZS',     &
  NGRID      = 4,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 2,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'MAX_ZS',         &
  CSTDNAME   = '',               &
  CLONGNAME  = 'MAX_ZS',         &
  CUNITS     = 'm',              &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_MAX_ZS',     &
  NGRID      = 4,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 2,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'MIN_ZS',         &
  CSTDNAME   = '',               &
  CLONGNAME  = 'MIN_ZS',         &
  CUNITS     = 'm',              &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_MIN_ZS',     &
  NGRID      = 4,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 2,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'SSO_STDEV',      &
  CSTDNAME   = '',               &
  CLONGNAME  = 'SSO_STDEV',      &
  CUNITS     = 'm',              &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_SSO_STDEV',  &
  NGRID      = 4,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 2,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'INPRC',          &
  CSTDNAME   = '',               &
  CLONGNAME  = 'INPRC',          &
  CUNITS     = 'm s-1',          &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_INstantaneous Cloud Precipitation Rate', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 2,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'ACPRC',          &
  CSTDNAME   = '',               &
  CLONGNAME  = 'ACPRC',          &
  CUNITS     = 'm',              &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_ACcumulated Cloud Precipitation Rate', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 2,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'INDEP',          &
  CSTDNAME   = '',               &
  CLONGNAME  = 'INDEP',          &
  CUNITS     = 'm s-1',          &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_INstantaneous Cloud Deposition Rate', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 2,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'ACDEP',          &
  CSTDNAME   = '',               &
  CLONGNAME  = 'ACDEP',          &
  CUNITS     = 'm',              &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_ACcumulated Cloud Deposition Rate', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 2,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'INPRR',          &
  CSTDNAME   = 'rainfall_rate',  &
  CLONGNAME  = 'INPRR',          &
  CUNITS     = 'm s-1',          &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_INstantaneous Precipitation Rain Rate', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 2,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'INPRR3D',        &
  CSTDNAME   = '',               &
  CLONGNAME  = 'INPRR3D',        &
  CUNITS     = 'm s-1',          &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_INstantaneous 3D Rain Precipitation flux', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'EVAP3D',         &
  CSTDNAME   = '',               &
  CLONGNAME  = 'EVAP3D',         &
  CUNITS     = 'kg kg-1 s-1',    &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_INstantaneous 3D Rain Evaporation flux', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA(                          &
  CMNHNAME   = 'ACPRR',                                   &
  CSTDNAME   = 'thickness_of_rainfall_amount',            &
  CLONGNAME  = 'ACPRR',                                   &
  CUNITS     = 'm',                                       &
  CDIR       = 'XY',                                      &
  CCOMMENT   = 'X_Y_ACcumulated Precipitation Rain Rate', &
  NGRID      = 1,                                         &
  NTYPE      = TYPEREAL,                                  &
  NDIMS      = 2,                                         &
  LTIMEDEP   = .TRUE.                                     ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'INPRS',          &
  CSTDNAME   = '',               &
  CLONGNAME  = 'INPRS',          &
  CUNITS     = 'm s-1',          &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_INstantaneous PRecipitation Snow Rate', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 2,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'ACPRS',          &
  CSTDNAME   = '',               &
  CLONGNAME  = 'ACPRS',          &
  CUNITS     = 'm',              &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_ACcumulated PRecipitation Snow Rate', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 2,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'INPRG',          &
  CSTDNAME   = '',               &
  CLONGNAME  = 'INPRG',          &
  CUNITS     = 'm s-1',          &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_INstantaneous PRecipitation Graupel Rate', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 2,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'ACPRG',          &
  CSTDNAME   = '',               &
  CLONGNAME  = 'ACPRG',          &
  CUNITS     = 'm',              &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_ACcumulated PRecipitation Graupel Rate', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 2,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'INPRH',          &
  CSTDNAME   = '',               &
  CLONGNAME  = 'INPRH',          &
  CUNITS     = 'm s-1',          &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_INstantaneous PRecipitation Hail Rate', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 2,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'ACPRH',          &
  CSTDNAME   = '',               &
  CLONGNAME  = 'ACPRH',          &
  CUNITS     = 'm',              &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_ACcumulated PRecipitation Hail Rate', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 2,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'INPRT',          &
  CSTDNAME   = '',               &
  CLONGNAME  = 'INPRT',          &
  CUNITS     = 'm s-1',          &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Total INstantaneaous PRecipitation rate', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 2,                &
  LTIMEDEP   = .TRUE.            ) )
!No permanent variable associated to this field

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'ACPRT',          &
  CSTDNAME   = '',               &
  CLONGNAME  = 'ACPRT',          &
  CUNITS     = 'm',              &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Total ACcumulated PRecipitation rate', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 2,                &
  LTIMEDEP   = .TRUE.            ) )
!No permanent variable associated to this field

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'VT_FLX',         &
  CSTDNAME   = '',               &
  CLONGNAME  = 'VT_FLX',         &
  CUNITS     = 'K m s-1',        &
  CDIR       = 'XY',             &
  CCOMMENT   = '',               &
  NGRID      = 2,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'WT_FLX',         &
  CSTDNAME   = '',               &
  CLONGNAME  = 'WT_FLX',         &
  CUNITS     = 'K m s-1',        &
  CDIR       = 'XY',             &
  CCOMMENT   = '',               &
  NGRID      = 4,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'RTHS_EDDY_FLUX', &
  CSTDNAME   = '',               &
  CLONGNAME  = 'RTHS_EDDY_FLUX', &
  CUNITS     = '',               &
  CDIR       = 'XY',             &
  CCOMMENT   = '',               &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'VU_FLX',         &
  CSTDNAME   = '',               &
  CLONGNAME  = 'VU_FLX',         &
  CUNITS     = 'm s-2',          &
  CDIR       = 'XY',             &
  CCOMMENT   = '',               &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'RVS_EDDY_FLUX',  &
  CSTDNAME   = '',               &
  CLONGNAME  = 'RVS_EDDY_FLUX',  &
  CUNITS     = '',               &
  CDIR       = 'XY',             &
  CCOMMENT   = '',               &
  NGRID      = 3,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .TRUE.            ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'FRC',            &
  CSTDNAME   = '',               &
  CLONGNAME  = 'FRC',            &
  CUNITS     = '',               &
  CDIR       = '--',             &
  CCOMMENT   = 'Number of forcing profiles', &
  NGRID      = 0,                &
  NTYPE      = TYPEINT,          &
  NDIMS      = 0,                &
  LTIMEDEP   = .FALSE.           ) )
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
IF (TRIM(CPROGRAM)=='REAL' .OR. TRIM(CPROGRAM) == 'LFICDF') THEN
!
call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'UT15',           &
  CSTDNAME   = '',               &
  CLONGNAME  = 'UT15',           &
  CUNITS     = 'm s-1',          &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Z_U component of Total wind', &
  NGRID      = 2,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'VT15',           &
  CSTDNAME   = '',               &
  CLONGNAME  = 'VT15',           &
  CUNITS     = 'm s-1',          &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Z_V component of Total wind', &
  NGRID      = 3,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'TEMPTOT',        &
  CSTDNAME   = '',               &
  CLONGNAME  = 'TEMPTOT',        &
  CUNITS     = 'K',              &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Z_TOTal TEMPerature', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'PRESTOT',        &
  CSTDNAME   = '',               &
  CLONGNAME  = 'PRESTOT',        &
  CUNITS     = 'Pa',             &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Z_TOTal PRESsure', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 2,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'HUMTOT',         &
  CSTDNAME   = '',               &
  CLONGNAME  = 'HUMTOT',         &
  CUNITS     = 'kg kg-1',        &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Z_TOTal specific HUMidity', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'UT16',           &
  CSTDNAME   = '',               &
  CLONGNAME  = 'UT16',           &
  CUNITS     = 'm s-1',          &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Z_U component of Environmental wind', &
  NGRID      = 2,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'VT16',           &
  CSTDNAME   = '',               &
  CLONGNAME  = 'VT16',           &
  CUNITS     = 'm s-1',          &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Z_V component of Environmental wind', &
  NGRID      = 3,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'TEMPENV',        &
  CSTDNAME   = '',               &
  CLONGNAME  = 'TEMPENV',        &
  CUNITS     = 'K',              &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Z_ENVironmental TEMPerature', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'PRESENV',        &
  CSTDNAME   = '',               &
  CLONGNAME  = 'PRESENV',        &
  CUNITS     = 'Pa',             &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Z_ENVironmental PRESsure', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 2,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'HUMENV',         &
  CSTDNAME   = '',               &
  CLONGNAME  = 'HUMENV',         &
  CUNITS     = 'kg kg-1',        &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Z_ENVironmental specific HUMidity', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'UT17',           &
  CSTDNAME   = '',               &
  CLONGNAME  = 'UT17',           &
  CUNITS     = 'm s-1',          &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Z_U component of Basic wind', &
  NGRID      = 2,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'VT17',           &
  CSTDNAME   = '',               &
  CLONGNAME  = 'VT17',           &
  CUNITS     = 'm s-1',          &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Z_V component of Basic wind', &
  NGRID      = 3,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'TEMPBAS',        &
  CSTDNAME   = '',               &
  CLONGNAME  = 'TEMPBAS',        &
  CUNITS     = 'K',              &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Z_BASic TEMPerature', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'PRESBAS',        &
  CSTDNAME   = '',               &
  CLONGNAME  = 'PRESBAS',        &
  CUNITS     = 'Pa',             &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Z_BASic PRESsure', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'HUMBAS',         &
  CSTDNAME   = '',               &
  CLONGNAME  = 'HUMBAS',         &
  CUNITS     = 'kg kg-1',        &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Z_BASic specific HUMidity', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'VTDIS',          &
  CSTDNAME   = '',               &
  CLONGNAME  = 'VTDIS',          &
  CUNITS     = 'm s-1',          &
  CDIR       = 'XY',             &
  CCOMMENT   = 'X_Y_Z_Total disturbance tangential wind', &
  NGRID      = 1,                &
  NTYPE      = TYPEREAL,         &
  NDIMS      = 3,                &
  LTIMEDEP   = .FALSE.           ) )
!
!END IF !LFILTERING
END IF !CPROGRAM==REAL .OR. LFICDF
!
call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'NFRCLT',         &
  CSTDNAME   = '',               &
  CLONGNAME  = 'NFRCLT',         &
  CUNITS     = '1',              &
  CDIR       = '--',             &
  CCOMMENT   = 'number of sea surface forcings + 1', &
  NGRID      = 0,                &
  NTYPE      = TYPEINT,          &
  NDIMS      = 0,                &
  LTIMEDEP   = .FALSE.           ) )

call Add_field2list( TFIELDDATA( &
  CMNHNAME   = 'NINFRT',         &
  CSTDNAME   = '',               &
  CLONGNAME  = 'NINFRT',         &
  CUNITS     = 's',              &
  CDIR       = '--',             &
  CCOMMENT   = 'Interval in seconds between forcings', &
  NGRID      = 0,                &
  NTYPE      = TYPEINT,          &
  NDIMS      = 0,                &
  LTIMEDEP   = .FALSE.           ) )
!
!
WRITE(YMSG,'("number of used fields=",I4," out of ",I4)') nfields_used-1,NMAXFIELDS
CALL PRINT_MSG(NVERB_INFO,'GEN','INI_FIELD_LIST',TRIM(YMSG))
!
#if 0
!
call Add_field2list( TFIELDDATA( &
  CMNHNAME   = '',               &
  CSTDNAME   = '',               &
  CLONGNAME  = '',               &
  CUNITS     = '',               &
  CDIR       = '',               &
  CLBTYPE    = '',               &
  CCOMMENT   = '',               &
  NGRID      = ,                 &
  NTYPE      = ,                 &
  NDIMS      = ,                 &
  LTIMEDEP   = ,                 ) )
#endif
!
END SUBROUTINE INI_FIELD_LIST
!
SUBROUTINE FIND_FIELD_ID_FROM_MNHNAME(HMNHNAME,KID,KRESP,ONOWARNING)
!
CHARACTER(LEN=*),            INTENT(IN) :: HMNHNAME !Name of the field to find
INTEGER,                     INTENT(OUT):: KID      !Index of the field
INTEGER,                     INTENT(OUT):: KRESP    !Return-code 
LOGICAL, OPTIONAL,           INTENT(IN) :: ONOWARNING !If true, do not print warning
!
INTEGER :: IDX,JI
INTEGER :: ICOUNT
INTEGER,SAVE :: IFIRSTGUESS=1 !Store first field to test
CHARACTER(LEN=64) :: YMSG
LOGICAL :: GNOWARNING
!
!
IF (.NOT.LFIELDLIST_ISINIT) THEN
  CALL PRINT_MSG(NVERB_FATAL,'GEN','FIND_FIELD_ID_FROM_MNHNAME','TFIELDLIST not yet initialized')
END IF
!
KID = 0
KRESP = 0
ICOUNT = 0
IDX = IFIRSTGUESS
!
IF (PRESENT(ONOWARNING)) THEN
  GNOWARNING = ONOWARNING
ELSE
  GNOWARNING = .FALSE.
END IF
!
DO
  ICOUNT = ICOUNT + 1
  IF (TRIM(TFIELDLIST(IDX)%CMNHNAME)==TRIM(HMNHNAME)) THEN
    KID = IDX
    EXIT
  ELSE
    IDX = IDX + 1
    IF ( IDX > nfields_used ) IDX = 1
  END IF
  IF (IDX == IFIRSTGUESS) EXIT !All entries have been tested
END DO
!
IF (KID==0) THEN
  !Field not found
  KRESP = -1
  IF (.NOT.GNOWARNING) THEN
    CALL PRINT_MSG(NVERB_WARNING,'GEN','FIND_FIELD_ID_FROM_MNHNAME','field '//TRIM(HMNHNAME)//' not known')
  ELSE
    CALL PRINT_MSG(NVERB_DEBUG,'GEN','FIND_FIELD_ID_FROM_MNHNAME','field '//TRIM(HMNHNAME)//' not known (not unexpected)')
  END IF
ELSE
  IFIRSTGUESS = IDX + 1
  IF ( IFIRSTGUESS > nfields_used ) IFIRSTGUESS = 1
  WRITE(YMSG,'( "field ",A16," found after ",I4," attempt(s)" )') TRIM(HMNHNAME),ICOUNT
  CALL PRINT_MSG(NVERB_DEBUG,'GEN','FIND_FIELD_ID_FROM_MNHNAME',TRIM(YMSG))
END IF
!
END SUBROUTINE FIND_FIELD_ID_FROM_MNHNAME
!
!
SUBROUTINE ALLOC_FIELD_SCALARS
!
USE MODD_DYN_n
USE MODD_PARAM_n
!
CALL PRINT_MSG(NVERB_DEBUG,'GEN','ALLOC_FIELD_SCALARS','called')
!
IF (LFIELDLIST_ISINIT) THEN
  CALL PRINT_MSG(NVERB_FATAL,'GEN','ALLOC_FIELD_SCALARS','TFIELDLIST already initialized')
END IF
!
!
IF (.NOT.ASSOCIATED(NRIMX)) THEN
  CALL PRINT_MSG(NVERB_DEBUG,'GEN','ALLOC_FIELD_SCALARS',' NRIMX was not associated')
  ALLOCATE(NRIMX)
END IF
IF (.NOT.ASSOCIATED(NRIMY)) THEN
  CALL PRINT_MSG(NVERB_DEBUG,'GEN','ALLOC_FIELD_SCALARS',' NRIMY was not associated')
  ALLOCATE(NRIMY)
END IF
IF (.NOT.ASSOCIATED(LHORELAX_UVWTH)) THEN
  CALL PRINT_MSG(NVERB_DEBUG,'GEN','ALLOC_FIELD_SCALARS',' LHORELAX_UVWTH was not associated')
  ALLOCATE(LHORELAX_UVWTH)
END IF
IF (.NOT.ASSOCIATED(LHORELAX_TKE)) THEN
  CALL PRINT_MSG(NVERB_DEBUG,'GEN','ALLOC_FIELD_SCALARS',' LHORELAX_TKE was not associated')
  ALLOCATE(LHORELAX_TKE)
END IF
IF (.NOT.ASSOCIATED(CSURF)) THEN
  CALL PRINT_MSG(NVERB_DEBUG,'GEN','ALLOC_FIELD_SCALARS',' CSURF was not associated')
  ALLOCATE(CHARACTER(LEN=4) :: CSURF)
  CSURF = ''
END IF
!
END SUBROUTINE ALLOC_FIELD_SCALARS
!
!
SUBROUTINE INI_FIELD_SCALARS
!
USE MODD_DYN_n
USE MODD_FIELD_n
USE MODD_GRID_n
USE MODD_TIME_n
USE MODD_PARAM_n
!
INTEGER :: IID,IRESP
!
CALL PRINT_MSG(NVERB_DEBUG,'GEN','INI_FIELD_SCALARS','called')
!
IF (.NOT.LFIELDLIST_ISINIT) THEN
  CALL PRINT_MSG(NVERB_FATAL,'GEN','FIND_FIELD_ID_FROM_MNHNAME','TFIELDLIST not yet initialized')
END IF
!
!
IF (.NOT.ASSOCIATED(XZTOP)) THEN
  CALL PRINT_MSG(NVERB_DEBUG,'GEN','INI_FIELD_SCALARS',' XZTOP was not associated')
  ALLOCATE(XZTOP)
  call Goto_model_1field( 'ZTOP', 1, 1, XZTOP )
END IF
!
IF (.NOT.ASSOCIATED(LSLEVE)) THEN
  CALL PRINT_MSG(NVERB_DEBUG,'GEN','INI_FIELD_SCALARS',' LSLEVE was not associated')
  ALLOCATE(LSLEVE)
  call Goto_model_1field( 'SLEVE', 1, 1, LSLEVE )
END IF
!
IF (.NOT.ASSOCIATED(XLEN1)) THEN
  CALL PRINT_MSG(NVERB_DEBUG,'GEN','INI_FIELD_SCALARS',' XLEN1 was not associated')
  ALLOCATE(XLEN1)
  call Goto_model_1field( 'LEN1', 1, 1, XLEN1 )
END IF
!
IF (.NOT.ASSOCIATED(XLEN2)) THEN
  CALL PRINT_MSG(NVERB_DEBUG,'GEN','INI_FIELD_SCALARS',' XLEN2 was not associated')
  ALLOCATE(XLEN2)
  call Goto_model_1field( 'LEN2', 1, 1, XLEN2 )
END IF
!
IF (.NOT.ASSOCIATED(TDTMOD)) THEN
  CALL PRINT_MSG(NVERB_DEBUG,'GEN','INI_FIELD_SCALARS',' TDTMOD was not associated')
  ALLOCATE(TDTMOD)
  call Goto_model_1field( 'DTMOD', 1, 1, TDTMOD )
END IF
!
IF (.NOT.ASSOCIATED(TDTCUR)) THEN
  CALL PRINT_MSG(NVERB_DEBUG,'GEN','INI_FIELD_SCALARS',' TDTCUR was not associated')
  ALLOCATE(TDTCUR)
  call Goto_model_1field( 'DTCUR', 1, 1, TDTCUR )
END IF
!
IF (.NOT.ASSOCIATED(TDTRAD_FULL)) THEN
  CALL PRINT_MSG(NVERB_DEBUG,'GEN','INI_FIELD_SCALARS',' TDTRAD_FULL was not associated')
  ALLOCATE(TDTRAD_FULL)
  call Goto_model_1field( 'DTRAD_FULL', 1, 1, TDTRAD_FULL )
END IF
!
IF (.NOT.ASSOCIATED(TDTRAD_CLONLY)) THEN
  CALL PRINT_MSG(NVERB_DEBUG,'GEN','INI_FIELD_SCALARS',' TDTRAD_CLONLY was not associated')
  ALLOCATE(TDTRAD_CLONLY)
  call Goto_model_1field( 'DTRAD_CLLY', 1, 1, TDTRAD_CLONLY )
END IF
!
IF (.NOT.ASSOCIATED(TDTDCONV)) THEN
  CALL PRINT_MSG(NVERB_DEBUG,'GEN','INI_FIELD_SCALARS',' TDTDCONV was not associated')
  ALLOCATE(TDTDCONV)
  call Goto_model_1field( 'DTDCONV', 1, 1, TDTDCONV )
END IF
!
IF (.NOT.ASSOCIATED(CSURF)) THEN
  CALL PRINT_MSG(NVERB_DEBUG,'GEN','INI_FIELD_SCALARS',' CSURF was not associated')
  ALLOCATE(CHARACTER(LEN=4) :: CSURF)
  CSURF = ''
  call Goto_model_1field( 'SURF', 1, 1, CSURF )
END IF
!
IF (.NOT.ASSOCIATED(XDRYMASST)) THEN
  CALL PRINT_MSG(NVERB_DEBUG,'GEN','INI_FIELD_SCALARS',' XDRYMASST was not associated')
  ALLOCATE(XDRYMASST)
  call Goto_model_1field( 'DRYMASST', 1, 1, XDRYMASST )
END IF
!
IF (.NOT.ASSOCIATED(XDRYMASSS)) THEN
  CALL PRINT_MSG(NVERB_DEBUG,'GEN','INI_FIELD_SCALARS',' XDRYMASSS was not associated')
  ALLOCATE(XDRYMASSS)
  call Goto_model_1field( 'DRYMASSS', 1, 1, XDRYMASSS )
END IF
!
IF (.NOT.ASSOCIATED(NRIMX)) THEN
  CALL PRINT_MSG(NVERB_DEBUG,'GEN','INI_FIELD_SCALARS',' NRIMX was not associated')
  ALLOCATE(NRIMX)
  call Goto_model_1field( 'RIMX', 1, 1, NRIMX )
END IF
!
IF (.NOT.ASSOCIATED(NRIMY)) THEN
  CALL PRINT_MSG(NVERB_DEBUG,'GEN','INI_FIELD_SCALARS',' NRIMY was not associated')
  ALLOCATE(NRIMY)
  call Goto_model_1field( 'RIMY', 1, 1, NRIMY )
END IF
!
IF (.NOT.ASSOCIATED(LHORELAX_UVWTH)) THEN
  CALL PRINT_MSG(NVERB_DEBUG,'GEN','INI_FIELD_SCALARS',' LHORELAX_UVWTH was not associated')
  ALLOCATE(LHORELAX_UVWTH)
  call Goto_model_1field( 'HORELAX_UVWTH', 1, 1, LHORELAX_UVWTH )
END IF
!
IF (.NOT.ASSOCIATED(LHORELAX_TKE)) THEN
  CALL PRINT_MSG(NVERB_DEBUG,'GEN','INI_FIELD_SCALARS',' LHORELAX_TKE was not associated')
  ALLOCATE(LHORELAX_TKE)
  call Goto_model_1field( 'HORELAX_TKE', 1, 1, LHORELAX_TKE )
END IF
!
END SUBROUTINE INI_FIELD_SCALARS
!
!
SUBROUTINE FIELDLIST_GOTO_MODEL(KFROM, KTO)
!
USE MODD_REF
!
USE MODD_ADV_n
USE MODD_CONF_n
USE MODD_DEEP_CONVECTION_n
USE MODD_DEF_EDDY_FLUX_n
USE MODD_DEF_EDDYUV_FLUX_n
USE MODD_DYN_n
USE MODD_ELEC_n
USE MODD_FIELD_n
USE MODD_FIRE_n
USE MODD_GR_FIELD_n
USE MODD_GRID_n
USE MODD_HURR_FIELD_n
USE MODD_LIMA_PRECIP_SCAVENGING_n
USE MODD_LSFIELD_n
USE MODD_OCEANH
USE MODD_PARAM_n
USE MODD_PAST_FIELD_n
USE MODD_CH_PH_n
USE MODD_PRECIP_n
USE MODD_RADIATIONS_n
USE MODD_REF_n
USE MODD_TIME_n
USE MODD_TURB_n
!
INTEGER, INTENT(IN) :: KFROM, KTO
!
!LOGICAL,SAVE :: GFIRST_CALL=.TRUE.
INTEGER :: IID,IID2,IRESP
CHARACTER(LEN=64) :: YMSG
!
WRITE(YMSG,'( I4,"->",I4 )') KFROM,KTO
CALL PRINT_MSG(NVERB_DEBUG,'GEN','FIELDLIST_GOTO_MODEL',TRIM(YMSG))
!
! IF (GFIRST_CALL) THEN
!   !This is necessary because the first time this subroutine is called
!   !the TFIELDLIST is not yet initialized.
!   !The use of this subroutine is not useful the first timebecause the
!   !data for the fields has not yet been allocated.
!   GFIRST_CALL = .FALSE.
!   RETURN
! END IF
!
IF (.NOT.LFIELDLIST_ISINIT) THEN
  CALL PRINT_MSG(NVERB_WARNING,'GEN','FIELDLIST_GOTO_MODEL','TFIELDLIST not yet initialized')
  RETURN
END IF
!
! MODD_FIELD_n variables
!
call Goto_model_1field( 'ZWS',   kfrom, kto, xzws   )
call Goto_model_1field( 'UT',    kfrom, kto, xut    )
call Goto_model_1field( 'VT',    kfrom, kto, xvt    )
call Goto_model_1field( 'WT',    kfrom, kto, xwt    )
call Goto_model_1field( 'THT',   kfrom, kto, xtht   )
call Goto_model_1field( 'TKET',  kfrom, kto, xtket  )
call Goto_model_1field( 'PABST', kfrom, kto, xpabst )
call Goto_model_1field( 'PHIT',  kfrom, kto, xphit  )
call Goto_model_1field( 'RT',    kfrom, kto, xrt    )
!
CALL FIND_FIELD_ID_FROM_MNHNAME( 'RT', IID2, IRESP )

IF (CONF_MODEL(KFROM)%IDX_RVT>0) THEN
  CALL FIND_FIELD_ID_FROM_MNHNAME( 'RVT', IID, IRESP )
  call Extend_1field_x3d( tfieldlist(iid), Max( kfrom, kto ) )
  if ( Associated( TFIELDLIST(IID2)%TFIELD_X4D(KFROM)%DATA ) ) &
  TFIELDLIST(IID)%TFIELD_X3D(KFROM)%DATA => TFIELDLIST(IID2)%TFIELD_X4D(KFROM)%DATA(:,:,:,CONF_MODEL(KFROM)%IDX_RVT)
  if ( kfrom /= kto .and. Associated( TFIELDLIST(IID2)%TFIELD_X4D(KTO)%DATA ) ) &
  TFIELDLIST(IID)%TFIELD_X3D(KTO)%DATA   => TFIELDLIST(IID2)%TFIELD_X4D(KTO)%DATA(:,:,:,CONF_MODEL(KTO)%IDX_RVT)
END IF
IF (CONF_MODEL(KFROM)%IDX_RCT>0) THEN
  CALL FIND_FIELD_ID_FROM_MNHNAME( 'RCT', IID, IRESP )
  call Extend_1field_x3d( tfieldlist(iid), Max( kfrom, kto ) )
  if ( Associated( TFIELDLIST(IID2)%TFIELD_X4D(KFROM)%DATA ) ) &
  TFIELDLIST(IID)%TFIELD_X3D(KFROM)%DATA => TFIELDLIST(IID2)%TFIELD_X4D(KFROM)%DATA(:,:,:,CONF_MODEL(KFROM)%IDX_RCT)
  if ( kfrom /= kto .and. Associated( TFIELDLIST(IID2)%TFIELD_X4D(KTO)%DATA ) ) &
  TFIELDLIST(IID)%TFIELD_X3D(KTO)%DATA   => TFIELDLIST(IID2)%TFIELD_X4D(KTO)%DATA(:,:,:,CONF_MODEL(KTO)%IDX_RCT)
END IF
IF (CONF_MODEL(KFROM)%IDX_RRT>0) THEN
  CALL FIND_FIELD_ID_FROM_MNHNAME( 'RRT', IID, IRESP )
  call Extend_1field_x3d( tfieldlist(iid), Max( kfrom, kto ) )
  if ( Associated( TFIELDLIST(IID2)%TFIELD_X4D(KFROM)%DATA ) ) &
  TFIELDLIST(IID)%TFIELD_X3D(KFROM)%DATA => TFIELDLIST(IID2)%TFIELD_X4D(KFROM)%DATA(:,:,:,CONF_MODEL(KFROM)%IDX_RRT)
  if ( kfrom /= kto .and. Associated( TFIELDLIST(IID2)%TFIELD_X4D(KTO)%DATA ) ) &
  TFIELDLIST(IID)%TFIELD_X3D(KTO)%DATA   => TFIELDLIST(IID2)%TFIELD_X4D(KTO)%DATA(:,:,:,CONF_MODEL(KTO)%IDX_RRT)
END IF
IF (CONF_MODEL(KFROM)%IDX_RIT>0) THEN
  CALL FIND_FIELD_ID_FROM_MNHNAME( 'RIT', IID, IRESP )
  call Extend_1field_x3d( tfieldlist(iid), Max( kfrom, kto ) )
  if ( Associated( TFIELDLIST(IID2)%TFIELD_X4D(KFROM)%DATA ) ) &
  TFIELDLIST(IID)%TFIELD_X3D(KFROM)%DATA => TFIELDLIST(IID2)%TFIELD_X4D(KFROM)%DATA(:,:,:,CONF_MODEL(KFROM)%IDX_RIT)
  if ( kfrom /= kto .and. Associated( TFIELDLIST(IID2)%TFIELD_X4D(KTO)%DATA ) ) &
  TFIELDLIST(IID)%TFIELD_X3D(KTO)%DATA   => TFIELDLIST(IID2)%TFIELD_X4D(KTO)%DATA(:,:,:,CONF_MODEL(KTO)%IDX_RIT)
END IF
IF (CONF_MODEL(KFROM)%IDX_RST>0) THEN
  CALL FIND_FIELD_ID_FROM_MNHNAME( 'RST', IID, IRESP )
  call Extend_1field_x3d( tfieldlist(iid), Max( kfrom, kto ) )
  if ( Associated( TFIELDLIST(IID2)%TFIELD_X4D(KFROM)%DATA ) ) &
  TFIELDLIST(IID)%TFIELD_X3D(KFROM)%DATA => TFIELDLIST(IID2)%TFIELD_X4D(KFROM)%DATA(:,:,:,CONF_MODEL(KFROM)%IDX_RST)
  if ( kfrom /= kto .and. Associated( TFIELDLIST(IID2)%TFIELD_X4D(KTO)%DATA ) ) &
  TFIELDLIST(IID)%TFIELD_X3D(KTO)%DATA   => TFIELDLIST(IID2)%TFIELD_X4D(KTO)%DATA(:,:,:,CONF_MODEL(KTO)%IDX_RST)
END IF
IF (CONF_MODEL(KFROM)%IDX_RGT>0) THEN
  CALL FIND_FIELD_ID_FROM_MNHNAME( 'RGT', IID, IRESP )
  call Extend_1field_x3d( tfieldlist(iid), Max( kfrom, kto ) )
  if ( Associated( TFIELDLIST(IID2)%TFIELD_X4D(KFROM)%DATA ) ) &
  TFIELDLIST(IID)%TFIELD_X3D(KFROM)%DATA => TFIELDLIST(IID2)%TFIELD_X4D(KFROM)%DATA(:,:,:,CONF_MODEL(KFROM)%IDX_RGT)
  if ( kfrom /= kto .and. Associated( TFIELDLIST(IID2)%TFIELD_X4D(KTO)%DATA ) ) &
  TFIELDLIST(IID)%TFIELD_X3D(KTO)%DATA   => TFIELDLIST(IID2)%TFIELD_X4D(KTO)%DATA(:,:,:,CONF_MODEL(KTO)%IDX_RGT)
END IF
IF (CONF_MODEL(KFROM)%IDX_RHT>0) THEN
  CALL FIND_FIELD_ID_FROM_MNHNAME( 'RHT', IID, IRESP )
  call Extend_1field_x3d( tfieldlist(iid), Max( kfrom, kto ) )
  if ( Associated( TFIELDLIST(IID2)%TFIELD_X4D(KFROM)%DATA ) ) &
  TFIELDLIST(IID)%TFIELD_X3D(KFROM)%DATA => TFIELDLIST(IID2)%TFIELD_X4D(KFROM)%DATA(:,:,:,CONF_MODEL(KFROM)%IDX_RHT)
  if ( kfrom /= kto .and. Associated( TFIELDLIST(IID2)%TFIELD_X4D(KTO)%DATA ) ) &
  TFIELDLIST(IID)%TFIELD_X3D(KTO)%DATA   => TFIELDLIST(IID2)%TFIELD_X4D(KTO)%DATA(:,:,:,CONF_MODEL(KTO)%IDX_RHT)
END IF
!
call Goto_model_1field( 'SUPSATMAX', kfrom, kto, XSUPSAT )
call Goto_model_1field( 'NACT',      kfrom, kto, XNACT   )
call Goto_model_1field( 'SSPRO',     kfrom, kto, XSSPRO  )
call Goto_model_1field( 'NPRO',      kfrom, kto, XNPRO   )
call Goto_model_1field( 'SRCT',      kfrom, kto, XSRCT   )
call Goto_model_1field( 'SIGS',      kfrom, kto, XSIGS   )
!
IF (CPROGRAM == 'MESONH') THEN
  call Goto_model_1field( 'US_PRES', kfrom, kto, XRUS_PRES )
  call Goto_model_1field( 'VS_PRES', kfrom, kto, XRVS_PRES )
  call Goto_model_1field( 'WS_PRES', kfrom, kto, XRWS_PRES )
  call Goto_model_1field( 'THS_CLD', kfrom, kto, XRTHS_CLD )

  call Goto_model_1field( 'RS_CLD', kfrom, kto, XRRS_CLD )
  !
  CALL FIND_FIELD_ID_FROM_MNHNAME( 'RS_CLD', IID2, IRESP )

  IF (CONF_MODEL(KFROM)%IDX_RVT>0) THEN
    CALL FIND_FIELD_ID_FROM_MNHNAME('RVS_CLD',IID,IRESP)
    call Extend_1field_x3d( tfieldlist(iid), Max( kfrom, kto ) )
    if ( Associated( TFIELDLIST(IID2)%TFIELD_X4D(KFROM)%DATA ) ) &
    TFIELDLIST(IID)%TFIELD_X3D(KFROM)%DATA => TFIELDLIST(IID2)%TFIELD_X4D(KFROM)%DATA(:,:,:,CONF_MODEL(KFROM)%IDX_RVT)
    if ( kfrom /= kto .and. Associated( TFIELDLIST(IID2)%TFIELD_X4D(KTO)%DATA ) ) &
    TFIELDLIST(IID)%TFIELD_X3D(KTO)%DATA => TFIELDLIST(IID2)%TFIELD_X4D(KTO)%DATA(:,:,:,CONF_MODEL(KTO)%IDX_RVT)
  END IF
  IF (CONF_MODEL(KFROM)%IDX_RCT>0) THEN
    CALL FIND_FIELD_ID_FROM_MNHNAME('RCS_CLD',IID,IRESP)
    call Extend_1field_x3d( tfieldlist(iid), Max( kfrom, kto ) )
    if ( Associated( TFIELDLIST(IID2)%TFIELD_X4D(KFROM)%DATA ) ) &
    TFIELDLIST(IID)%TFIELD_X3D(KFROM)%DATA => TFIELDLIST(IID2)%TFIELD_X4D(KFROM)%DATA(:,:,:,CONF_MODEL(KFROM)%IDX_RCT)
    if ( kfrom /= kto .and. Associated( TFIELDLIST(IID2)%TFIELD_X4D(KTO)%DATA ) ) &
    TFIELDLIST(IID)%TFIELD_X3D(KTO)%DATA => TFIELDLIST(IID2)%TFIELD_X4D(KTO)%DATA(:,:,:,CONF_MODEL(KTO)%IDX_RCT)
  END IF
  IF (CONF_MODEL(KFROM)%IDX_RRT>0) THEN
    CALL FIND_FIELD_ID_FROM_MNHNAME('RRS_CLD',IID,IRESP)
    call Extend_1field_x3d( tfieldlist(iid), Max( kfrom, kto ) )
    if ( Associated( TFIELDLIST(IID2)%TFIELD_X4D(KFROM)%DATA ) ) &
    TFIELDLIST(IID)%TFIELD_X3D(KFROM)%DATA => TFIELDLIST(IID2)%TFIELD_X4D(KFROM)%DATA(:,:,:,CONF_MODEL(KFROM)%IDX_RRT)
    if ( kfrom /= kto .and. Associated( TFIELDLIST(IID2)%TFIELD_X4D(KTO)%DATA ) ) &
    TFIELDLIST(IID)%TFIELD_X3D(KTO)%DATA => TFIELDLIST(IID2)%TFIELD_X4D(KTO)%DATA(:,:,:,CONF_MODEL(KTO)%IDX_RRT)
  END IF
  IF (CONF_MODEL(KFROM)%IDX_RIT>0) THEN
    CALL FIND_FIELD_ID_FROM_MNHNAME('RIS_CLD',IID,IRESP)
    call Extend_1field_x3d( tfieldlist(iid), Max( kfrom, kto ) )
    if ( Associated( TFIELDLIST(IID2)%TFIELD_X4D(KFROM)%DATA ) ) &
    TFIELDLIST(IID)%TFIELD_X3D(KFROM)%DATA => TFIELDLIST(IID2)%TFIELD_X4D(KFROM)%DATA(:,:,:,CONF_MODEL(KFROM)%IDX_RIT)
    if ( kfrom /= kto .and. Associated( TFIELDLIST(IID2)%TFIELD_X4D(KTO)%DATA ) ) &
    TFIELDLIST(IID)%TFIELD_X3D(KTO)%DATA => TFIELDLIST(IID2)%TFIELD_X4D(KTO)%DATA(:,:,:,CONF_MODEL(KTO)%IDX_RIT)
  END IF
  IF (CONF_MODEL(KFROM)%IDX_RST>0) THEN
    CALL FIND_FIELD_ID_FROM_MNHNAME('RSS_CLD',IID,IRESP)
    call Extend_1field_x3d( tfieldlist(iid), Max( kfrom, kto ) )
    if ( Associated( TFIELDLIST(IID2)%TFIELD_X4D(KFROM)%DATA ) ) &
    TFIELDLIST(IID)%TFIELD_X3D(KFROM)%DATA => TFIELDLIST(IID2)%TFIELD_X4D(KFROM)%DATA(:,:,:,CONF_MODEL(KFROM)%IDX_RST)
    if ( kfrom /= kto .and. Associated( TFIELDLIST(IID2)%TFIELD_X4D(KTO)%DATA ) ) &
    TFIELDLIST(IID)%TFIELD_X3D(KTO)%DATA => TFIELDLIST(IID2)%TFIELD_X4D(KTO)%DATA(:,:,:,CONF_MODEL(KTO)%IDX_RST)
  END IF
  IF (CONF_MODEL(KFROM)%IDX_RGT>0) THEN
    CALL FIND_FIELD_ID_FROM_MNHNAME('RGS_CLD',IID,IRESP)
    call Extend_1field_x3d( tfieldlist(iid), Max( kfrom, kto ) )
    if ( Associated( TFIELDLIST(IID2)%TFIELD_X4D(KFROM)%DATA ) ) &
    TFIELDLIST(IID)%TFIELD_X3D(KFROM)%DATA => TFIELDLIST(IID2)%TFIELD_X4D(KFROM)%DATA(:,:,:,CONF_MODEL(KFROM)%IDX_RGT)
    if ( kfrom /= kto .and. Associated( TFIELDLIST(IID2)%TFIELD_X4D(KTO)%DATA ) ) &
    TFIELDLIST(IID)%TFIELD_X3D(KTO)%DATA => TFIELDLIST(IID2)%TFIELD_X4D(KTO)%DATA(:,:,:,CONF_MODEL(KTO)%IDX_RGT)
  END IF
  IF (CONF_MODEL(KFROM)%IDX_RHT>0) THEN
    CALL FIND_FIELD_ID_FROM_MNHNAME('RHS_CLD',IID,IRESP)
    call Extend_1field_x3d( tfieldlist(iid), Max( kfrom, kto ) )
    if ( Associated( TFIELDLIST(IID2)%TFIELD_X4D(KFROM)%DATA ) ) &
    TFIELDLIST(IID)%TFIELD_X3D(KFROM)%DATA => TFIELDLIST(IID2)%TFIELD_X4D(KFROM)%DATA(:,:,:,CONF_MODEL(KFROM)%IDX_RHT)
    if ( kfrom /= kto .and. Associated( TFIELDLIST(IID2)%TFIELD_X4D(KTO)%DATA ) ) &
    TFIELDLIST(IID)%TFIELD_X3D(KTO)%DATA => TFIELDLIST(IID2)%TFIELD_X4D(KTO)%DATA(:,:,:,CONF_MODEL(KTO)%IDX_RHT)
  END IF

  call Goto_model_1field( 'CLDFR',  kfrom, kto, XCLDFR  )
  call Goto_model_1field( 'ICEFR',  kfrom, kto, XICEFR  )
  call Goto_model_1field( 'CIT',    kfrom, kto, XCIT    )
  call Goto_model_1field( 'RAINFR', kfrom, kto, XRAINFR )
END IF
!
! MODD_PAST_FIELD_n variables
!
call Goto_model_1field( 'UM',  kfrom, kto, XUM  )
call Goto_model_1field( 'VM',  kfrom, kto, XVM  )
call Goto_model_1field( 'WM',  kfrom, kto, XWM  )
call Goto_model_1field( 'DUM', kfrom, kto, XDUM )
call Goto_model_1field( 'DVM', kfrom, kto, XDVM )
call Goto_model_1field( 'DWM', kfrom, kto, XDWM )
!
! MODD_LIMA_PRECIP_SCAVENGING_n variables
!
call Goto_model_1field( 'INPAP', kfrom, kto, XINPAP )
call Goto_model_1field( 'ACPAP', kfrom, kto, XACPAP )
!
! MODD_ELEC_n variables
!
call Goto_model_1field( 'EFIELDU',  kfrom, kto, XEFIELDU  )
call Goto_model_1field( 'EFIELDV',  kfrom, kto, XEFIELDV  )
call Goto_model_1field( 'EFIELDW',  kfrom, kto, XEFIELDW  )
call Goto_model_1field( 'NI_IAGGS', kfrom, kto, XNI_IAGGS )
call Goto_model_1field( 'NI_IDRYG', kfrom, kto, XNI_IDRYG )
call Goto_model_1field( 'NI_SDRYG', kfrom, kto, XNI_SDRYG )
call Goto_model_1field( 'INDUC_CG', kfrom, kto, XIND_RATE )
!
! MODD_CH_PH_n variables
!
call Goto_model_1field( 'PHC', kfrom, kto, XPHC )
call Goto_model_1field( 'PHR', kfrom, kto, XPHR )
!
! MODD_LSFIELD_n variables
!
call Goto_model_1field( 'LSUM',  kfrom, kto, XLSUM  )
call Goto_model_1field( 'LSVM',  kfrom, kto, XLSVM  )
call Goto_model_1field( 'LSWM',  kfrom, kto, XLSWM  )
call Goto_model_1field( 'LSTHM', kfrom, kto, XLSTHM )
IF (LUSERV) THEN
  call Goto_model_1field( 'LSRVM', kfrom, kto, XLSRVM )
END IF
call Goto_model_1field( 'LBXUM',  kfrom, kto, XLBXUM  )
call Goto_model_1field( 'LBXVM',  kfrom, kto, XLBXVM  )
call Goto_model_1field( 'LBXWM',  kfrom, kto, XLBXWM  )
call Goto_model_1field( 'LBYUM',  kfrom, kto, XLBYUM  )
call Goto_model_1field( 'LBYVM',  kfrom, kto, XLBYVM  )
call Goto_model_1field( 'LBYWM',  kfrom, kto, XLBYWM  )
call Goto_model_1field( 'LBXTHM', kfrom, kto, XLBXTHM )
call Goto_model_1field( 'LBYTHM', kfrom, kto, XLBYTHM )

call Goto_model_1field( 'DRYMASST', kfrom, kto, XDRYMASST )
call Goto_model_1field( 'DRYMASSS', kfrom, kto, XDRYMASSS )
!
! MODD_DYN_n variables
!
call Goto_model_1field( 'RIMX', kfrom, kto, NRIMX )
call Goto_model_1field( 'RIMY', kfrom, kto, NRIMY )
call Goto_model_1field( 'HORELAX_UVWTH', kfrom, kto, LHORELAX_UVWTH )
call Goto_model_1field( 'HORELAX_TKE',   kfrom, kto, LHORELAX_TKE   )
!
! MODD_ADV_n variables
!
call Goto_model_1field( 'TKEMS', kfrom, kto, XRTKEMS  )
!
! MODD_GRID_n variables
!
call Goto_model_1field( 'ZS'  ,  kfrom, kto, XZS    )
call Goto_model_1field( 'ZSMT',  kfrom, kto, XZSMT  )
call Goto_model_1field( 'XHAT',  kfrom, kto, XXHAT  )
call Goto_model_1field( 'YHAT',  kfrom, kto, XYHAT  )
call Goto_model_1field( 'XHATM', kfrom, kto, XXHATM )
call Goto_model_1field( 'YHATM', kfrom, kto, XYHATM )
call Goto_model_1field( 'ZHAT',  kfrom, kto, XZHAT  )
call Goto_model_1field( 'ZHATM', kfrom, kto, XZHATM )
call Goto_model_1field( 'HAT_BOUND',  kfrom, kto, XHAT_BOUND  )
call Goto_model_1field( 'HATM_BOUND', kfrom, kto, XHATM_BOUND )
call Goto_model_1field( 'ZTOP',  kfrom, kto, XZTOP  )
call Goto_model_1field( 'DXHAT', kfrom, kto, XDXHAT )
call Goto_model_1field( 'DYHAT', kfrom, kto, XDYHAT )
call Goto_model_1field( 'SLEVE', kfrom, kto, LSLEVE )
call Goto_model_1field( 'LEN1',  kfrom, kto, XLEN1  )
call Goto_model_1field( 'LEN2',  kfrom, kto, XLEN2  )
call Goto_model_1field( 'ALT',   kfrom, kto, XZZ    )
call Goto_model_1field( 'DIRCOSXW', kfrom, kto, XDIRCOSXW )
call Goto_model_1field( 'DIRCOSYW', kfrom, kto, XDIRCOSYW )
call Goto_model_1field( 'DIRCOSZW', kfrom, kto, XDIRCOSZW )
call Goto_model_1field( 'COSSLOPE', kfrom, kto, XCOSSLOPE )
call Goto_model_1field( 'SINSLOPE', kfrom, kto, XSINSLOPE )
call Goto_model_1field( 'MAP',      kfrom, kto, XMAP )
call Goto_model_1field( 'LAT',      kfrom, kto, XLAT )
call Goto_model_1field( 'LON',      kfrom, kto, XLON )
call Goto_model_1field( 'XHAT_ll',  kfrom, kto, XXHAT_ll  )
call Goto_model_1field( 'YHAT_ll',  kfrom, kto, XYHAT_ll  )
call Goto_model_1field( 'XHATM_ll', kfrom, kto, XXHATM_ll )
call Goto_model_1field( 'YHATM_ll', kfrom, kto, XYHATM_ll )
!
! MODD_TIME_n variables
!
call Goto_model_1field( 'DTMOD',      kfrom, kto, TDTMOD        )
call Goto_model_1field( 'DTCUR',      kfrom, kto, TDTCUR        )
call Goto_model_1field( 'DTRAD_FULL', kfrom, kto, TDTRAD_FULL   )
call Goto_model_1field( 'DTRAD_CLLY', kfrom, kto, TDTRAD_CLONLY )
call Goto_model_1field( 'DTDCONV',    kfrom, kto, TDTDCONV      )
!
! MODD_PARAM_n variables
!
call Goto_model_1field( 'SURF', kfrom, kto, CSURF )
!
! MODD_TURB_n variables
!
call Goto_model_1field( 'BL_DEPTH',  kfrom, kto, XBL_DEPTH  )
call Goto_model_1field( 'SBL_DEPTH', kfrom, kto, XSBL_DEPTH )
call Goto_model_1field( 'WTHVMF',    kfrom, kto, XWTHVMF    )
call Goto_model_1field( 'WUMF',      kfrom, kto, XWUMF      )
call Goto_model_1field( 'WVMF',      kfrom, kto, XWVMF      )
!
! MODD_REF_n variables
!
call Goto_model_1field( 'RHODREF', kfrom, kto, XRHODREF )
call Goto_model_1field( 'THVREF',  kfrom, kto, XTHVREF  )
!
! MODD_RADIATIONS_n variables
!
IF (CPROGRAM=='MESONH') THEN
  call Goto_model_1field( 'DTHRAD',       kfrom, kto, XDTHRAD       )
  call Goto_model_1field( 'FLALWD',       kfrom, kto, XFLALWD       )
  call Goto_model_1field( 'DIRFLASWD',    kfrom, kto, XDIRFLASWD    )
  call Goto_model_1field( 'SCAFLASWD',    kfrom, kto, XSCAFLASWD    )
  call Goto_model_1field( 'DIRSRFSWD',    kfrom, kto, XDIRSRFSWD    )
  call Goto_model_1field( 'CLEARCOL_TM1', kfrom, kto, NCLEARCOL_TM1 )
  call Goto_model_1field( 'ZENITH',       kfrom, kto, XZENITH       )
  call Goto_model_1field( 'AZIM',         kfrom, kto, XAZIM         )
  call Goto_model_1field( 'DIR_ALB',      kfrom, kto, XDIR_ALB      )
  call Goto_model_1field( 'SCA_ALB',      kfrom, kto, XSCA_ALB      )
  call Goto_model_1field( 'EMIS',         kfrom, kto, XEMIS         )
  call Goto_model_1field( 'TSRAD',        kfrom, kto, XTSRAD        )
END IF
!
! MODD_FIRE_n variables
!
IF (     TRIM(CPROGRAM) == 'MESONH' .OR. TRIM(CPROGRAM) == 'DIAG'  .OR. TRIM(CPROGRAM) == 'REAL' &
    .OR. TRIM(CPROGRAM) == 'LFICDF' .OR. TRIM(CPROGRAM) == 'SPAWN'                               ) THEN
  call Goto_model_1field( 'FMREFINRATIOX', kfrom, kto, NREFINX )
  call Goto_model_1field( 'FMREFINRATIOY', kfrom, kto, NREFINY )
  call Goto_model_1field( 'FMPHI',  kfrom, kto, XLSPHI )
  call Goto_model_1field( 'FMBMAP', kfrom, kto, XBMAP )
  call Goto_model_1field( 'FMASE',  kfrom, kto, XFMASE )
  call Goto_model_1field( 'FMAWC',  kfrom, kto, XFMAWC )
  call Goto_model_1field( 'FMWINDU', kfrom, kto, XFMWINDU )
  call Goto_model_1field( 'FMWINDV', kfrom, kto, XFMWINDV )
  call Goto_model_1field( 'FMWINDW', kfrom, kto, XFMWINDW )
  call Goto_model_1field( 'FMHWS',  kfrom, kto, XFMHWS )
  call Goto_model_1field( 'FMROS',  kfrom, kto, XFIRERW )
  call Goto_model_1field( 'FMROS0', kfrom, kto, XFMR0 )
  call Goto_model_1field( 'FMFLUXHDH', kfrom, kto, XFMFLUXHDH )
  call Goto_model_1field( 'FMFLUXHDW', kfrom, kto, XFMFLUXHDW )
  call Goto_model_1field( 'FMGRADOROX', kfrom, kto, XFMGRADOROX )
  call Goto_model_1field( 'FMGRADOROY', kfrom, kto, XFMGRADOROY )
END IF
!
! MODD_DEEP_CONVECTION_n variables
!
IF (TRIM(CPROGRAM) /= 'PGD' .AND. TRIM(CPROGRAM) /= 'NESPGD' .AND. TRIM(CPROGRAM) /= 'SPAWN') THEN
  call Goto_model_1field( 'COUNTCONV',     kfrom, kto, NCOUNTCONV       )
  call Goto_model_1field( 'DTHCONV',       kfrom, kto, XDTHCONV         )
  call Goto_model_1field( 'DRVCONV',       kfrom, kto, XDRVCONV         )
  call Goto_model_1field( 'DRCCONV',       kfrom, kto, XDRCCONV         )
  call Goto_model_1field( 'DRICONV',       kfrom, kto, XDRICONV         )
  call Goto_model_1field( 'PRCONV',        kfrom, kto, XPRCONV          )
  call Goto_model_1field( 'PACCONV',       kfrom, kto, XPACCONV         )
  call Goto_model_1field( 'PRSCONV',       kfrom, kto, XPRSCONV         )
  call Goto_model_1field( 'DSVCONV',       kfrom, kto, XDSVCONV         )
  call Goto_model_1field( 'PRLFLXCONV',    kfrom, kto, XPRLFLXCONV      )
  call Goto_model_1field( 'PRSFLXCONV',    kfrom, kto, XPRSFLXCONV      )
  call Goto_model_1field( 'UMFCONV',       kfrom, kto, XUMFCONV         )
  call Goto_model_1field( 'DMFCONV',       kfrom, kto, XDMFCONV         )
  call Goto_model_1field( 'MFCONV',        kfrom, kto, XMFCONV          )
  call Goto_model_1field( 'CAPE',          kfrom, kto, XCAPE            )
  call Goto_model_1field( 'CLTOPCONV_LVL', kfrom, kto, NCLTOPCONV       )
  call Goto_model_1field( 'CLBASCONV_LVL', kfrom, kto, NCLBASCONV       )
  call Goto_model_1field( 'IC_RATE',       kfrom, kto, XIC_RATE         )
  call Goto_model_1field( 'CG_RATE',       kfrom, kto, XCG_RATE         )
  call Goto_model_1field( 'IC_TOTAL_NB',   kfrom, kto, XIC_TOTAL_NUMBER )
  call Goto_model_1field( 'CG_TOTAL_NB',   kfrom, kto, XCG_TOTAL_NUMBER )
END IF
!
! MODD_GR_FIELD_n variables
!
call Goto_model_1field( 'SSO_ANIS',  kfrom, kto, XSSO_ANISOTROPY )
call Goto_model_1field( 'SSO_SLOPE', kfrom, kto, XSSO_SLOPE      )
call Goto_model_1field( 'SSO_DIR',   kfrom, kto, XSSO_DIRECTION  )
call Goto_model_1field( 'AVG_ZS',    kfrom, kto, XAVG_ZS         )
call Goto_model_1field( 'SIL_ZS',    kfrom, kto, XSIL_ZS         )
call Goto_model_1field( 'MAX_ZS',    kfrom, kto, XMAX_ZS         )
call Goto_model_1field( 'MIN_ZS',    kfrom, kto, XMIN_ZS         )
call Goto_model_1field( 'SSO_STDEV', kfrom, kto, XSSO_STDEV      )
!
! MODD_PRECIP_n variables
!
call Goto_model_1field( 'INPRC',   kfrom, kto, XINPRC   )
call Goto_model_1field( 'ACPRC',   kfrom, kto, XACPRC   )
call Goto_model_1field( 'INDEP',   kfrom, kto, XINDEP   )
call Goto_model_1field( 'ACDEP',   kfrom, kto, XACDEP   )
call Goto_model_1field( 'INPRR',   kfrom, kto, XINPRR   )
call Goto_model_1field( 'INPRR3D', kfrom, kto, XINPRR3D )
call Goto_model_1field( 'EVAP3D',  kfrom, kto, XEVAP3D  )
call Goto_model_1field( 'ACPRR',   kfrom, kto, XACPRR   )
call Goto_model_1field( 'INPRS',   kfrom, kto, XINPRS   )
call Goto_model_1field( 'ACPRS',   kfrom, kto, XACPRS   )
call Goto_model_1field( 'INPRG',   kfrom, kto, XINPRG   )
call Goto_model_1field( 'ACPRG',   kfrom, kto, XACPRG   )
call Goto_model_1field( 'INPRH',   kfrom, kto, XINPRH   )
call Goto_model_1field( 'ACPRH',   kfrom, kto, XACPRH   )
!
! MODD_DEF_EDDY_FLUX_n variables
!
call Goto_model_1field( 'VT_FLX',         kfrom, kto, XVTH_FLUX_M      )
call Goto_model_1field( 'WT_FLX',         kfrom, kto, XWTH_FLUX_M      )
call Goto_model_1field( 'RTHS_EDDY_FLUX', kfrom, kto, XRTHS_EDDY_FLUX  )
!
! MODD_DEF_EDDYUV_FLUX_n variables
!
call Goto_model_1field( 'VU_FLX',        kfrom, kto, XVU_FLUX_M     )
call Goto_model_1field( 'RVS_EDDY_FLUX', kfrom, kto, XRVS_EDDY_FLUX )
!
! MODD_HURR_FIELD_n variables
!
IF (CPROGRAM=='REAL') THEN
  call Goto_model_1field( 'UT15',    kfrom, kto, XUTOT   )
  call Goto_model_1field( 'VT15',    kfrom, kto, XVTOT   )
  call Goto_model_1field( 'TEMPTOT', kfrom, kto, XTTOT   )
  call Goto_model_1field( 'PRESTOT', kfrom, kto, XPTOT   )
  call Goto_model_1field( 'HUMTOT',  kfrom, kto, XQTOT   )
  call Goto_model_1field( 'UT16',    kfrom, kto, XUENV   )
  call Goto_model_1field( 'VT16',    kfrom, kto, XVENV   )
  call Goto_model_1field( 'TEMPENV', kfrom, kto, XTENV   )
  call Goto_model_1field( 'PRESENV', kfrom, kto, XPENV   )
  call Goto_model_1field( 'HUMENV',  kfrom, kto, XQENV   )
  call Goto_model_1field( 'UT17',    kfrom, kto, XUBASIC )
  call Goto_model_1field( 'VT17',    kfrom, kto, XVBASIC )
  call Goto_model_1field( 'TEMPBAS', kfrom, kto, XTBASIC )
  call Goto_model_1field( 'PRESBAS', kfrom, kto, XPBASIC )
  call Goto_model_1field( 'HUMBAS',  kfrom, kto, XQBASIC )
  call Goto_model_1field( 'VTDIS',   kfrom, kto, XVTDIS  )
END IF
!
END SUBROUTINE FIELDLIST_GOTO_MODEL


subroutine Add_field2list( tpfield )

implicit none

type(tfielddata) :: tpfield

character(len=64) :: ymsg
type(tfielddata), allocatable, dimension(:) :: tzfieldlistnew

!Check if tfieldlist big enough and enlarge it if necessary
if ( nfields_used >= NMAXFIELDS ) then
  Allocate( tzfieldlistnew(nmaxfields + nmaxfieldstep) )
  tzfieldlistnew(1 : nmaxfields) = tfieldlist(1 : nmaxfields)
  call Move_alloc( from = tzfieldlistnew, to = tfieldlist )
  nmaxfields = nmaxfields + nmaxfieldstep
  Write( ymsg, '( "nmaxfields increased from ", i5, " to ", i5 )') nmaxfields - nmaxfieldstep, nmaxfields
  call Print_msg( NVERB_DEBUG, 'GEN', 'Add_field2list', Trim( ymsg ) )
end if

nfields_used = nfields_used + 1

tfieldlist(nfields_used) = tpfield

end subroutine Add_field2list


subroutine Goto_model_1field_c0d( hname, kfrom, kto, pdata )

implicit none

character(len=*),          intent(in)    :: hname
integer,                   intent(in)    :: kfrom
integer,                   intent(in)    :: kto
character(len=*), pointer, intent(inout) :: pdata

integer :: iid
integer :: iresp

call Find_field_id_from_mnhname( hname, iid, iresp )

call Extend_1field_c0d( tfieldlist(iid), Max( kfrom, kto ) )

tfieldlist(iid)%tfield_c0d(kfrom)%data => pdata
if ( kfrom /= kto ) then
  if ( .not. Associated( tfieldlist(iid)%tfield_c0d(kto)%data ) ) then
    Allocate( character(len=Len(pdata)) :: tfieldlist(iid)%tfield_c0d(kto)%data )
    tfieldlist(iid)%tfield_c0d(kto)%data(:) = ''
  end if
  pdata => tfieldlist(iid)%tfield_c0d(kto)%data
end if

end subroutine Goto_model_1field_c0d


subroutine Goto_model_1field_c1d( hname, kfrom, kto, pdata )

implicit none

character(len=*),                        intent(in)    :: hname
integer,                                 intent(in)    :: kfrom
integer,                                 intent(in)    :: kto
character(len=*), dimension(:), pointer, intent(inout) :: pdata

integer :: iid
integer :: iresp
integer :: ji

call Find_field_id_from_mnhname( hname, iid, iresp )

call Extend_1field_c1d( tfieldlist(iid), Max( kfrom, kto ) )

tfieldlist(iid)%tfield_c1d(kfrom)%data => pdata
if ( kfrom /= kto ) then
  if ( .not. Associated( tfieldlist(iid)%tfield_c1d(kto)%data ) ) then
    Allocate( character(len=Len(pdata)) :: tfieldlist(iid)%tfield_c1d(kto)%data(Size(pdata)) )
    do ji = 1, Size(pdata)
      tfieldlist(iid)%tfield_c1d(kto)%data(ji) = ''
    end do
  end if
  pdata => tfieldlist(iid)%tfield_c1d(kto)%data
end if

end subroutine Goto_model_1field_c1d


subroutine Goto_model_1field_l0d( hname, kfrom, kto, pdata )

implicit none

character(len=*),          intent(in)    :: hname
integer,                   intent(in)    :: kfrom
integer,                   intent(in)    :: kto
logical,          pointer, intent(inout) :: pdata

integer :: iid
integer :: iresp

call Find_field_id_from_mnhname( hname, iid, iresp )

call Extend_1field_l0d( tfieldlist(iid), Max( kfrom, kto ) )

tfieldlist(iid)%tfield_l0d(kfrom)%data => pdata
if ( kfrom /= kto ) pdata => tfieldlist(iid)%tfield_l0d(kto)%data

end subroutine Goto_model_1field_l0d


subroutine Goto_model_1field_l1d( hname, kfrom, kto, pdata )

implicit none

character(len=*),               intent(in)    :: hname
integer,                        intent(in)    :: kfrom
integer,                        intent(in)    :: kto
logical, dimension(:), pointer, intent(inout) :: pdata

integer :: iid
integer :: iresp

call Find_field_id_from_mnhname( hname, iid, iresp )

call Extend_1field_l1d( tfieldlist(iid), Max( kfrom, kto ) )

tfieldlist(iid)%tfield_l1d(kfrom)%data => pdata
if ( kfrom /= kto ) pdata => tfieldlist(iid)%tfield_l1d(kto)%data

end subroutine Goto_model_1field_l1d


subroutine Goto_model_1field_n0d( hname, kfrom, kto, pdata )

implicit none

character(len=*),          intent(in)    :: hname
integer,                   intent(in)    :: kfrom
integer,                   intent(in)    :: kto
integer,          pointer, intent(inout) :: pdata

integer :: iid
integer :: iresp

call Find_field_id_from_mnhname( hname, iid, iresp )

call Extend_1field_n0d( tfieldlist(iid), Max( kfrom, kto ) )

tfieldlist(iid)%tfield_n0d(kfrom)%data => pdata
if ( kfrom /= kto ) pdata => tfieldlist(iid)%tfield_n0d(kto)%data

end subroutine Goto_model_1field_n0d


subroutine Goto_model_1field_n1d( hname, kfrom, kto, pdata )

implicit none

character(len=*),               intent(in)    :: hname
integer,                        intent(in)    :: kfrom
integer,                        intent(in)    :: kto
integer, dimension(:), pointer, intent(inout) :: pdata

integer :: iid
integer :: iresp

call Find_field_id_from_mnhname( hname, iid, iresp )

call Extend_1field_n1d( tfieldlist(iid), Max( kfrom, kto ) )

tfieldlist(iid)%tfield_n1d(kfrom)%data => pdata
if ( kfrom /= kto ) pdata => tfieldlist(iid)%tfield_n1d(kto)%data

end subroutine Goto_model_1field_n1d


subroutine Goto_model_1field_n2d( hname, kfrom, kto, pdata )

implicit none

character(len=*),                 intent(in)    :: hname
integer,                          intent(in)    :: kfrom
integer,                          intent(in)    :: kto
integer, dimension(:,:), pointer, intent(inout) :: pdata

integer :: iid
integer :: iresp

call Find_field_id_from_mnhname( hname, iid, iresp )

call Extend_1field_n2d( tfieldlist(iid), Max( kfrom, kto ) )

tfieldlist(iid)%tfield_n2d(kfrom)%data => pdata
if ( kfrom /= kto ) pdata => tfieldlist(iid)%tfield_n2d(kto)%data

end subroutine Goto_model_1field_n2d


subroutine Goto_model_1field_n3d( hname, kfrom, kto, pdata )

implicit none

character(len=*),                   intent(in)    :: hname
integer,                            intent(in)    :: kfrom
integer,                            intent(in)    :: kto
integer, dimension(:,:,:), pointer, intent(inout) :: pdata

integer :: iid
integer :: iresp

call Find_field_id_from_mnhname( hname, iid, iresp )

call Extend_1field_n3d( tfieldlist(iid), Max( kfrom, kto ) )

tfieldlist(iid)%tfield_n3d(kfrom)%data => pdata
if ( kfrom /= kto ) pdata => tfieldlist(iid)%tfield_n3d(kto)%data

end subroutine Goto_model_1field_n3d


subroutine Goto_model_1field_t0d( hname, kfrom, kto, pdata )

use modd_type_date, only: date_time

implicit none

character(len=*),          intent(in)    :: hname
integer,                   intent(in)    :: kfrom
integer,                   intent(in)    :: kto
type(date_time),  pointer, intent(inout) :: pdata

integer :: iid
integer :: iresp

call Find_field_id_from_mnhname( hname, iid, iresp )

call Extend_1field_t0d( tfieldlist(iid), Max( kfrom, kto ) )

tfieldlist(iid)%tfield_t0d(kfrom)%data => pdata
if ( kfrom /= kto ) pdata => tfieldlist(iid)%tfield_t0d(kto)%data

end subroutine Goto_model_1field_t0d


subroutine Goto_model_1field_t1d( hname, kfrom, kto, pdata )

use modd_type_date, only: date_time

implicit none

character(len=*),                       intent(in)    :: hname
integer,                                intent(in)    :: kfrom
integer,                                intent(in)    :: kto
type(date_time), dimension(:), pointer, intent(inout) :: pdata

integer :: iid
integer :: iresp

call Find_field_id_from_mnhname( hname, iid, iresp )

call Extend_1field_t1d( tfieldlist(iid), Max( kfrom, kto ) )

tfieldlist(iid)%tfield_t1d(kfrom)%data => pdata
if ( kfrom /= kto ) pdata => tfieldlist(iid)%tfield_t1d(kto)%data

end subroutine Goto_model_1field_t1d


subroutine Goto_model_1field_x0d( hname, kfrom, kto, pdata )

implicit none

character(len=*),          intent(in)    :: hname
integer,                   intent(in)    :: kfrom
integer,                   intent(in)    :: kto
real,             pointer, intent(inout) :: pdata

integer :: iid
integer :: iresp

call Find_field_id_from_mnhname( hname, iid, iresp )

call Extend_1field_x0d( tfieldlist(iid), Max( kfrom, kto ) )

tfieldlist(iid)%tfield_x0d(kfrom)%data => pdata
if ( kfrom /= kto ) pdata => tfieldlist(iid)%tfield_x0d(kto)%data

end subroutine Goto_model_1field_x0d


subroutine Goto_model_1field_x1d( hname, kfrom, kto, pdata )

implicit none

character(len=*),               intent(in)    :: hname
integer,                        intent(in)    :: kfrom
integer,                        intent(in)    :: kto
real,    dimension(:), pointer, intent(inout) :: pdata

integer :: iid
integer :: iresp

call Find_field_id_from_mnhname( hname, iid, iresp )

call Extend_1field_x1d( tfieldlist(iid), Max( kfrom, kto ) )

tfieldlist(iid)%tfield_x1d(kfrom)%data => pdata
if ( kfrom /= kto ) pdata => tfieldlist(iid)%tfield_x1d(kto)%data

end subroutine Goto_model_1field_x1d


subroutine Goto_model_1field_x2d( hname, kfrom, kto, pdata )

implicit none

character(len=*),                 intent(in)    :: hname
integer,                          intent(in)    :: kfrom
integer,                          intent(in)    :: kto
real,    dimension(:,:), pointer, intent(inout) :: pdata

integer :: iid
integer :: iresp

call Find_field_id_from_mnhname( hname, iid, iresp )

call Extend_1field_x2d( tfieldlist(iid), Max( kfrom, kto ) )

tfieldlist(iid)%tfield_x2d(kfrom)%data => pdata
if ( kfrom /= kto ) pdata => tfieldlist(iid)%tfield_x2d(kto)%data

end subroutine Goto_model_1field_x2d


subroutine Goto_model_1field_x3d( hname, kfrom, kto, pdata )

implicit none

character(len=*),                   intent(in)    :: hname
integer,                            intent(in)    :: kfrom
integer,                            intent(in)    :: kto
real,    dimension(:,:,:), pointer, intent(inout) :: pdata

integer :: iid
integer :: iresp

call Find_field_id_from_mnhname( hname, iid, iresp )

call Extend_1field_x3d( tfieldlist(iid), Max( kfrom, kto ) )

tfieldlist(iid)%tfield_x3d(kfrom)%data => pdata
if ( kfrom /= kto ) pdata => tfieldlist(iid)%tfield_x3d(kto)%data

end subroutine Goto_model_1field_x3d


subroutine Goto_model_1field_x4d( hname, kfrom, kto, pdata )

implicit none

character(len=*),                     intent(in)    :: hname
integer,                              intent(in)    :: kfrom
integer,                              intent(in)    :: kto
real,    dimension(:,:,:,:), pointer, intent(inout) :: pdata

integer :: iid
integer :: iresp

call Find_field_id_from_mnhname( hname, iid, iresp )

call Extend_1field_x4d( tfieldlist(iid), Max( kfrom, kto ) )

tfieldlist(iid)%tfield_x4d(kfrom)%data => pdata
if ( kfrom /= kto ) pdata => tfieldlist(iid)%tfield_x4d(kto)%data

end subroutine Goto_model_1field_x4d


subroutine Goto_model_1field_x5d( hname, kfrom, kto, pdata )

implicit none

character(len=*),                       intent(in)    :: hname
integer,                                intent(in)    :: kfrom
integer,                                intent(in)    :: kto
real,    dimension(:,:,:,:,:), pointer, intent(inout) :: pdata

integer :: iid
integer :: iresp

call Find_field_id_from_mnhname( hname, iid, iresp )

call Extend_1field_x5d( tfieldlist(iid), Max( kfrom, kto ) )

tfieldlist(iid)%tfield_x5d(kfrom)%data => pdata
if ( kfrom /= kto ) pdata => tfieldlist(iid)%tfield_x5d(kto)%data

end subroutine Goto_model_1field_x5d


subroutine Goto_model_1field_x6d( hname, kfrom, kto, pdata )

implicit none

character(len=*),                         intent(in)    :: hname
integer,                                  intent(in)    :: kfrom
integer,                                  intent(in)    :: kto
real,    dimension(:,:,:,:,:,:), pointer, intent(inout) :: pdata

integer :: iid
integer :: iresp

call Find_field_id_from_mnhname( hname, iid, iresp )

call Extend_1field_x6d( tfieldlist(iid), Max( kfrom, kto ) )

tfieldlist(iid)%tfield_x6d(kfrom)%data => pdata
if ( kfrom /= kto ) pdata => tfieldlist(iid)%tfield_x6d(kto)%data

end subroutine Goto_model_1field_x6d


subroutine Extend_1field_c0d( tpfield, ksize )

implicit none

type(tfielddata), intent(inout) :: tpfield
integer,          intent(in)    :: ksize

integer :: ji
type(tfieldptr_c0d), dimension(:), allocatable :: tzfield_c0d

if ( tpfield%nmodelmax < 0 ) then
  !nmodelmax is < 0 if the allocation of the field has been done by hand
  !(not using a constructor, default value of nmodelmax)
  !The correct value of nmodelmax is hence computed here
  tpfield%nmodelmax = Size( tpfield%tfield_c0d )
end if

if ( ksize > tpfield%nmodelmax ) then
  if ( tpfield%nmodelmax == 0 ) then
    Allocate( tpfield%tfield_c0d(ksize) )
    do ji = 1, ksize
      tpfield%tfield_c0d(ji)%data => null()
    end do
  else
    Allocate( tzfield_c0d(ksize) )
    do ji = 1, Size( tpfield%tfield_c0d)
      tzfield_c0d(ji)%data => tpfield%tfield_c0d(ji)%data
    end do
    do ji = Size( tpfield%tfield_c0d) + 1, ksize
      tzfield_c0d(ji)%data => null()
    end do
    call Move_alloc( from = tzfield_c0d, to = tpfield%tfield_c0d )
  end if
  tpfield%nmodelmax = ksize
end if

end subroutine Extend_1field_c0d


subroutine Extend_1field_c1d( tpfield, ksize )

implicit none

type(tfielddata), intent(inout) :: tpfield
integer,          intent(in)    :: ksize

integer :: ji
type(tfieldptr_c1d), dimension(:), allocatable :: tzfield_c1d

if ( tpfield%nmodelmax < 0 ) then
  !nmodelmax is < 0 if the allocation of the field has been done by hand
  !(not using a constructor, default value of nmodelmax)
  !The correct value of nmodelmax is hence computed here
  tpfield%nmodelmax = Size( tpfield%tfield_c1d )
end if

if ( ksize > tpfield%nmodelmax ) then
  if ( tpfield%nmodelmax == 0 ) then
    Allocate( tpfield%tfield_c1d(ksize) )
    do ji = 1, ksize
      tpfield%tfield_c1d(ji)%data => null()
    end do
  else
    Allocate( tzfield_c1d(ksize) )
    do ji = 1, Size( tpfield%tfield_c1d)
      tzfield_c1d(ji)%data => tpfield%tfield_c1d(ji)%data
    end do
    do ji = Size( tpfield%tfield_c1d) + 1, ksize
      tzfield_c1d(ji)%data => null()
    end do
    call Move_alloc( from = tzfield_c1d, to = tpfield%tfield_c1d )
  end if
  tpfield%nmodelmax = ksize
end if

end subroutine Extend_1field_c1d


subroutine Extend_1field_l0d( tpfield, ksize )

implicit none

type(tfielddata), intent(inout) :: tpfield
integer,          intent(in)    :: ksize

integer :: ji
type(tfieldptr_l0d), dimension(:), allocatable :: tzfield_l0d

if ( tpfield%nmodelmax < 0 ) then
  !nmodelmax is < 0 if the allocation of the field has been done by hand
  !(not using a constructor, default value of nmodelmax)
  !The correct value of nmodelmax is hence computed here
  tpfield%nmodelmax = Size( tpfield%tfield_l0d )
end if

if ( ksize > tpfield%nmodelmax ) then
  if ( tpfield%nmodelmax == 0 ) then
    Allocate( tpfield%tfield_l0d(ksize) )
    do ji = 1, ksize
      ! tpfield%tfield_l0d(ji)%data => null()
      Allocate( tpfield%tfield_l0d(ji)%data )
      tpfield%tfield_l0d(ji)%data = .false.
    end do
  else
    Allocate( tzfield_l0d(ksize) )
    do ji = 1, Size( tpfield%tfield_l0d)
      tzfield_l0d(ji)%data => tpfield%tfield_l0d(ji)%data
    end do
    do ji = Size( tpfield%tfield_l0d) + 1, ksize
      ! tzfield_l0d(ji)%data => null()
      Allocate( tzfield_l0d(ji)%data )
      tzfield_l0d(ji)%data = .false.
    end do
    call Move_alloc( from = tzfield_l0d, to = tpfield%tfield_l0d )
  end if
  tpfield%nmodelmax = ksize
end if

end subroutine Extend_1field_l0d


subroutine Extend_1field_l1d( tpfield, ksize )

implicit none

type(tfielddata), intent(inout) :: tpfield
integer,          intent(in)    :: ksize

integer :: ji
type(tfieldptr_l1d), dimension(:), allocatable :: tzfield_l1d

if ( tpfield%nmodelmax < 0 ) then
  !nmodelmax is < 0 if the allocation of the field has been done by hand
  !(not using a constructor, default value of nmodelmax)
  !The correct value of nmodelmax is hence computed here
  tpfield%nmodelmax = Size( tpfield%tfield_l1d )
end if

if ( ksize > tpfield%nmodelmax ) then
  if ( tpfield%nmodelmax == 0 ) then
    Allocate( tpfield%tfield_l1d(ksize) )
    do ji = 1, ksize
      tpfield%tfield_l1d(ji)%data => null()
    end do
  else
    Allocate( tzfield_l1d(ksize) )
    do ji = 1, Size( tpfield%tfield_l1d)
      tzfield_l1d(ji)%data => tpfield%tfield_l1d(ji)%data
    end do
    do ji = Size( tpfield%tfield_l1d) + 1, ksize
      tzfield_l1d(ji)%data => null()
    end do
    call Move_alloc( from = tzfield_l1d, to = tpfield%tfield_l1d )
  end if
  tpfield%nmodelmax = ksize
end if

end subroutine Extend_1field_l1d


subroutine Extend_1field_n0d( tpfield, ksize )

use modd_parameters, only: NUNDEF

implicit none

type(tfielddata), intent(inout) :: tpfield
integer,          intent(in)    :: ksize

integer :: ji
type(tfieldptr_n0d), dimension(:), allocatable :: tzfield_n0d

if ( tpfield%nmodelmax < 0 ) then
  !nmodelmax is < 0 if the allocation of the field has been done by hand
  !(not using a constructor, default value of nmodelmax)
  !The correct value of nmodelmax is hence computed here
  tpfield%nmodelmax = Size( tpfield%tfield_n0d )
end if

if ( ksize > tpfield%nmodelmax ) then
  if ( tpfield%nmodelmax == 0 ) then
    Allocate( tpfield%tfield_n0d(ksize) )
    do ji = 1, ksize
      ! tpfield%tfield_n0d(ji)%data => null()
      Allocate( tpfield%tfield_n0d(ji)%data )
      tpfield%tfield_n0d(ji)%data = NUNDEF
    end do
  else
    Allocate( tzfield_n0d(ksize) )
    do ji = 1, Size( tpfield%tfield_n0d)
      tzfield_n0d(ji)%data => tpfield%tfield_n0d(ji)%data
    end do
    do ji = Size( tpfield%tfield_n0d) + 1, ksize
      ! tzfield_n0d(ji)%data => null()
      Allocate( tzfield_n0d(ji)%data )
      tzfield_n0d(ji)%data = NUNDEF
    end do
    call Move_alloc( from = tzfield_n0d, to = tpfield%tfield_n0d )
  end if
  tpfield%nmodelmax = ksize
end if

end subroutine Extend_1field_n0d


subroutine Extend_1field_n1d( tpfield, ksize )

implicit none

type(tfielddata), intent(inout) :: tpfield
integer,          intent(in)    :: ksize

integer :: ji
type(tfieldptr_n1d), dimension(:), allocatable :: tzfield_n1d

if ( tpfield%nmodelmax < 0 ) then
  !nmodelmax is < 0 if the allocation of the field has been done by hand
  !(not using a constructor, default value of nmodelmax)
  !The correct value of nmodelmax is hence computed here
  tpfield%nmodelmax = Size( tpfield%tfield_n1d )
end if

if ( ksize > tpfield%nmodelmax ) then
  if ( tpfield%nmodelmax == 0 ) then
    Allocate( tpfield%tfield_n1d(ksize) )
    do ji = 1, ksize
      tpfield%tfield_n1d(ji)%data => null()
    end do
  else
    Allocate( tzfield_n1d(ksize) )
    do ji = 1, Size( tpfield%tfield_n1d)
      tzfield_n1d(ji)%data => tpfield%tfield_n1d(ji)%data
    end do
    do ji = Size( tpfield%tfield_n1d) + 1, ksize
      tzfield_n1d(ji)%data => null()
    end do
    call Move_alloc( from = tzfield_n1d, to = tpfield%tfield_n1d )
  end if
  tpfield%nmodelmax = ksize
end if

end subroutine Extend_1field_n1d


subroutine Extend_1field_n2d( tpfield, ksize )

implicit none

type(tfielddata), intent(inout) :: tpfield
integer,          intent(in)    :: ksize

integer :: ji
type(tfieldptr_n2d), dimension(:), allocatable :: tzfield_n2d

if ( tpfield%nmodelmax < 0 ) then
  !nmodelmax is < 0 if the allocation of the field has been done by hand
  !(not using a constructor, default value of nmodelmax)
  !The correct value of nmodelmax is hence computed here
  tpfield%nmodelmax = Size( tpfield%tfield_n2d )
end if

if ( ksize > tpfield%nmodelmax ) then
  if ( tpfield%nmodelmax == 0 ) then
    Allocate( tpfield%tfield_n2d(ksize) )
    do ji = 1, ksize
      tpfield%tfield_n2d(ji)%data => null()
    end do
  else
    Allocate( tzfield_n2d(ksize) )
    do ji = 1, Size( tpfield%tfield_n2d)
      tzfield_n2d(ji)%data => tpfield%tfield_n2d(ji)%data
    end do
    do ji = Size( tpfield%tfield_n2d) + 1, ksize
      tzfield_n2d(ji)%data => null()
    end do
    call Move_alloc( from = tzfield_n2d, to = tpfield%tfield_n2d )
  end if
  tpfield%nmodelmax = ksize
end if

end subroutine Extend_1field_n2d


subroutine Extend_1field_n3d( tpfield, ksize )

implicit none

type(tfielddata), intent(inout) :: tpfield
integer,          intent(in)    :: ksize

integer :: ji
type(tfieldptr_n3d), dimension(:), allocatable :: tzfield_n3d

if ( tpfield%nmodelmax < 0 ) then
  !nmodelmax is < 0 if the allocation of the field has been done by hand
  !(not using a constructor, default value of nmodelmax)
  !The correct value of nmodelmax is hence computed here
  tpfield%nmodelmax = Size( tpfield%tfield_n3d )
end if

if ( ksize > tpfield%nmodelmax ) then
  if ( tpfield%nmodelmax == 0 ) then
    Allocate( tpfield%tfield_n3d(ksize) )
    do ji = 1, ksize
      tpfield%tfield_n3d(ji)%data => null()
    end do
  else
    Allocate( tzfield_n3d(ksize) )
    do ji = 1, Size( tpfield%tfield_n3d)
      tzfield_n3d(ji)%data => tpfield%tfield_n3d(ji)%data
    end do
    do ji = Size( tpfield%tfield_n3d) + 1, ksize
      tzfield_n3d(ji)%data => null()
    end do
    call Move_alloc( from = tzfield_n3d, to = tpfield%tfield_n3d )
  end if
  tpfield%nmodelmax = ksize
end if

end subroutine Extend_1field_n3d


subroutine Extend_1field_t0d( tpfield, ksize )

implicit none

type(tfielddata), intent(inout) :: tpfield
integer,          intent(in)    :: ksize

integer :: ji
type(tfieldptr_t0d), dimension(:), allocatable :: tzfield_t0d

if ( tpfield%nmodelmax < 0 ) then
  !nmodelmax is < 0 if the allocation of the field has been done by hand
  !(not using a constructor, default value of nmodelmax)
  !The correct value of nmodelmax is hence computed here
  tpfield%nmodelmax = Size( tpfield%tfield_t0d )
end if

if ( ksize > tpfield%nmodelmax ) then
  if ( tpfield%nmodelmax == 0 ) then
    Allocate( tpfield%tfield_t0d(ksize) )
    do ji = 1, ksize
      ! tpfield%tfield_t0d(ji)%data => null()
      Allocate( tpfield%tfield_t0d(ji)%data )
    end do
  else
    Allocate( tzfield_t0d(ksize) )
    do ji = 1, Size( tpfield%tfield_t0d)
      tzfield_t0d(ji)%data => tpfield%tfield_t0d(ji)%data
    end do
    do ji = Size( tpfield%tfield_t0d) + 1, ksize
      ! tzfield_t0d(ji)%data => null()
      Allocate( tzfield_t0d(ji)%data )
    end do
    call Move_alloc( from = tzfield_t0d, to = tpfield%tfield_t0d )
  end if
  tpfield%nmodelmax = ksize
end if

end subroutine Extend_1field_t0d


subroutine Extend_1field_t1d( tpfield, ksize )

implicit none

type(tfielddata), intent(inout) :: tpfield
integer,          intent(in)    :: ksize

integer :: ji
type(tfieldptr_t1d), dimension(:), allocatable :: tzfield_t1d

if ( tpfield%nmodelmax < 0 ) then
  !nmodelmax is < 0 if the allocation of the field has been done by hand
  !(not using a constructor, default value of nmodelmax)
  !The correct value of nmodelmax is hence computed here
  tpfield%nmodelmax = Size( tpfield%tfield_t1d )
end if

if ( ksize > tpfield%nmodelmax ) then
  if ( tpfield%nmodelmax == 0 ) then
    Allocate( tpfield%tfield_t1d(ksize) )
    do ji = 1, ksize
      tpfield%tfield_t1d(ji)%data => null()
    end do
  else
    Allocate( tzfield_t1d(ksize) )
    do ji = 1, Size( tpfield%tfield_t1d)
      tzfield_t1d(ji)%data => tpfield%tfield_t1d(ji)%data
    end do
    do ji = Size( tpfield%tfield_t1d) + 1, ksize
      tzfield_t1d(ji)%data => null()
    end do
    call Move_alloc( from = tzfield_t1d, to = tpfield%tfield_t1d )
  end if
  tpfield%nmodelmax = ksize
end if

end subroutine Extend_1field_t1d


subroutine Extend_1field_x0d( tpfield, ksize )

use modd_parameters, only: XUNDEF

implicit none

type(tfielddata), intent(inout) :: tpfield
integer,          intent(in)    :: ksize

integer :: ji
type(tfieldptr_x0d), dimension(:), allocatable :: tzfield_x0d

if ( tpfield%nmodelmax < 0 ) then
  !nmodelmax is < 0 if the allocation of the field has been done by hand
  !(not using a constructor, default value of nmodelmax)
  !The correct value of nmodelmax is hence computed here
  tpfield%nmodelmax = Size( tpfield%tfield_x0d )
end if

if ( ksize > tpfield%nmodelmax ) then
  if ( tpfield%nmodelmax == 0 ) then
    Allocate( tpfield%tfield_x0d(ksize) )
    do ji = 1, ksize
      ! tpfield%tfield_x0d(ji)%data => null()
      Allocate( tpfield%tfield_x0d(ji)%data )
      tpfield%tfield_x0d(ji)%data = XUNDEF
    end do
  else
    Allocate( tzfield_x0d(ksize) )
    do ji = 1, Size( tpfield%tfield_x0d)
      tzfield_x0d(ji)%data => tpfield%tfield_x0d(ji)%data
    end do
    do ji = Size( tpfield%tfield_x0d) + 1, ksize
      ! tzfield_x0d(ji)%data => null()
      Allocate( tzfield_x0d(ji)%data )
      tzfield_x0d(ji)%data = XUNDEF
    end do
    call Move_alloc( from = tzfield_x0d, to = tpfield%tfield_x0d )
  end if
  tpfield%nmodelmax = ksize
end if

end subroutine Extend_1field_x0d


subroutine Extend_1field_x1d( tpfield, ksize )

implicit none

type(tfielddata), intent(inout) :: tpfield
integer,          intent(in)    :: ksize

integer :: ji
type(tfieldptr_x1d), dimension(:), allocatable :: tzfield_x1d

if ( tpfield%nmodelmax < 0 ) then
  !nmodelmax is < 0 if the allocation of the field has been done by hand
  !(not using a constructor, default value of nmodelmax)
  !The correct value of nmodelmax is hence computed here
  tpfield%nmodelmax = Size( tpfield%tfield_x1d )
end if

if ( ksize > tpfield%nmodelmax ) then
  if ( tpfield%nmodelmax == 0 ) then
    Allocate( tpfield%tfield_x1d(ksize) )
    do ji = 1, ksize
      tpfield%tfield_x1d(ji)%data => null()
    end do
  else
    Allocate( tzfield_x1d(ksize) )
    do ji = 1, Size( tpfield%tfield_x1d)
      tzfield_x1d(ji)%data => tpfield%tfield_x1d(ji)%data
    end do
    do ji = Size( tpfield%tfield_x1d) + 1, ksize
      tzfield_x1d(ji)%data => null()
    end do
    call Move_alloc( from = tzfield_x1d, to = tpfield%tfield_x1d )
  end if
  tpfield%nmodelmax = ksize
end if

end subroutine Extend_1field_x1d


subroutine Extend_1field_x2d( tpfield, ksize )

implicit none

type(tfielddata), intent(inout) :: tpfield
integer,          intent(in)    :: ksize

integer :: ji
type(tfieldptr_x2d), dimension(:), allocatable :: tzfield_x2d

if ( tpfield%nmodelmax < 0 ) then
  !nmodelmax is < 0 if the allocation of the field has been done by hand
  !(not using a constructor, default value of nmodelmax)
  !The correct value of nmodelmax is hence computed here
  tpfield%nmodelmax = Size( tpfield%tfield_x2d )
end if

if ( ksize > tpfield%nmodelmax ) then
  if ( tpfield%nmodelmax == 0 ) then
    Allocate( tpfield%tfield_x2d(ksize) )
    do ji = 1, ksize
      tpfield%tfield_x2d(ji)%data => null()
    end do
  else
    Allocate( tzfield_x2d(ksize) )
    do ji = 1, Size( tpfield%tfield_x2d)
      tzfield_x2d(ji)%data => tpfield%tfield_x2d(ji)%data
    end do
    do ji = Size( tpfield%tfield_x2d) + 1, ksize
      tzfield_x2d(ji)%data => null()
    end do
    call Move_alloc( from = tzfield_x2d, to = tpfield%tfield_x2d )
  end if
  tpfield%nmodelmax = ksize
end if

end subroutine Extend_1field_x2d


subroutine Extend_1field_x3d( tpfield, ksize )

implicit none

type(tfielddata), intent(inout) :: tpfield
integer,          intent(in)    :: ksize

integer :: ji
type(tfieldptr_x3d), dimension(:), allocatable :: tzfield_x3d

if ( tpfield%nmodelmax < 0 ) then
  !nmodelmax is < 0 if the allocation of the field has been done by hand
  !(not using a constructor, default value of nmodelmax)
  !The correct value of nmodelmax is hence computed here
  tpfield%nmodelmax = Size( tpfield%tfield_x3d )
end if

if ( ksize > tpfield%nmodelmax ) then
  if ( tpfield%nmodelmax == 0 ) then
    Allocate( tpfield%tfield_x3d(ksize) )
    do ji = 1, ksize
      tpfield%tfield_x3d(ji)%data => null()
    end do
  else
    Allocate( tzfield_x3d(ksize) )
    do ji = 1, Size( tpfield%tfield_x3d)
      tzfield_x3d(ji)%data => tpfield%tfield_x3d(ji)%data
    end do
    do ji = Size( tpfield%tfield_x3d) + 1, ksize
      tzfield_x3d(ji)%data => null()
    end do
    call Move_alloc( from = tzfield_x3d, to = tpfield%tfield_x3d )
  end if
  tpfield%nmodelmax = ksize
end if

end subroutine Extend_1field_x3d


subroutine Extend_1field_x4d( tpfield, ksize )

implicit none

type(tfielddata), intent(inout) :: tpfield
integer,          intent(in)    :: ksize

integer :: ji
type(tfieldptr_x4d), dimension(:), allocatable :: tzfield_x4d

if ( tpfield%nmodelmax < 0 ) then
  !nmodelmax is < 0 if the allocation of the field has been done by hand
  !(not using a constructor, default value of nmodelmax)
  !The correct value of nmodelmax is hence computed here
  tpfield%nmodelmax = Size( tpfield%tfield_x4d )
end if

if ( ksize > tpfield%nmodelmax ) then
  if ( tpfield%nmodelmax == 0 ) then
    Allocate( tpfield%tfield_x4d(ksize) )
    do ji = 1, ksize
      tpfield%tfield_x4d(ji)%data => null()
    end do
  else
    Allocate( tzfield_x4d(ksize) )
    do ji = 1, Size( tpfield%tfield_x4d)
      tzfield_x4d(ji)%data => tpfield%tfield_x4d(ji)%data
    end do
    do ji = Size( tpfield%tfield_x4d) + 1, ksize
      tzfield_x4d(ji)%data => null()
    end do
    call Move_alloc( from = tzfield_x4d, to = tpfield%tfield_x4d )
  end if
  tpfield%nmodelmax = ksize
end if

end subroutine Extend_1field_x4d


subroutine Extend_1field_x5d( tpfield, ksize )

implicit none

type(tfielddata), intent(inout) :: tpfield
integer,          intent(in)    :: ksize

integer :: ji
type(tfieldptr_x5d), dimension(:), allocatable :: tzfield_x5d

if ( tpfield%nmodelmax < 0 ) then
  !nmodelmax is < 0 if the allocation of the field has been done by hand
  !(not using a constructor, default value of nmodelmax)
  !The correct value of nmodelmax is hence computed here
  tpfield%nmodelmax = Size( tpfield%tfield_x5d )
end if

if ( ksize > tpfield%nmodelmax ) then
  if ( tpfield%nmodelmax == 0 ) then
    Allocate( tpfield%tfield_x5d(ksize) )
    do ji = 1, ksize
      tpfield%tfield_x5d(ji)%data => null()
    end do
  else
    Allocate( tzfield_x5d(ksize) )
    do ji = 1, Size( tpfield%tfield_x5d)
      tzfield_x5d(ji)%data => tpfield%tfield_x5d(ji)%data
    end do
    do ji = Size( tpfield%tfield_x5d) + 1, ksize
      tzfield_x5d(ji)%data => null()
    end do
    call Move_alloc( from = tzfield_x5d, to = tpfield%tfield_x5d )
  end if
  tpfield%nmodelmax = ksize
end if

end subroutine Extend_1field_x5d


subroutine Extend_1field_x6d( tpfield, ksize )

implicit none

type(tfielddata), intent(inout) :: tpfield
integer,          intent(in)    :: ksize

integer :: ji
type(tfieldptr_x6d), dimension(:), allocatable :: tzfield_x6d

if ( tpfield%nmodelmax < 0 ) then
  !nmodelmax is < 0 if the allocation of the field has been done by hand
  !(not using a constructor, default value of nmodelmax)
  !The correct value of nmodelmax is hence computed here
  tpfield%nmodelmax = Size( tpfield%tfield_x6d )
end if

if ( ksize > tpfield%nmodelmax ) then
  if ( tpfield%nmodelmax == 0 ) then
    Allocate( tpfield%tfield_x6d(ksize) )
    do ji = 1, ksize
      tpfield%tfield_x6d(ji)%data => null()
    end do
  else
    Allocate( tzfield_x6d(ksize) )
    do ji = 1, Size( tpfield%tfield_x6d)
      tzfield_x6d(ji)%data => tpfield%tfield_x6d(ji)%data
    end do
    do ji = Size( tpfield%tfield_x6d) + 1, ksize
      tzfield_x6d(ji)%data => null()
    end do
    call Move_alloc( from = tzfield_x6d, to = tpfield%tfield_x6d )
  end if
  tpfield%nmodelmax = ksize
end if

end subroutine Extend_1field_x6d

end module mode_field
