!MNH_LIC Copyright 1998-2020 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
! Modifications:
!  J. Escobar  15/09/2015: WENO5 & JPHEXT <> 1
!  P. Wautelet 10/04/2019: replace ABORT and STOP calls by Print_msg
!  P. Wautelet 26/04/2019: use modd_precision parameters for datatypes of MPI communications
!-----------------------------------------------------------------

!     ####################
      MODULE MODE_TOOLS_ll
!     ####################
!
!!    Purpose
!!    -------
!
!     The Purpose of this module is to provide subroutines and functions for
!     for the initialization of parallel data variables
!
!!    Routines Of The User Interface
!!    ------------------------------
!
!     FUNCTIONS   : LNORTH_ll, LSOUTH_ll, LEAST_ll, LWEST_ll
!     SUBROUTINES : GET_DIM_EXT_ll, GET_OR_ll, GET_INDICE_ll, GET_PHYSICAL_ll
!                   GET_INTERSECTION_ll,
!                   GET_GLOBALSLICE_ll
!                     (GET_1DGLOBALSLICE_ll, GET_2DGLOBALSLICE_ll),
!                   GET_SLICE_ll
!                     (GET_1DSLICE_ll, GET_2DSLICE_ll)
!                   GET_L2_NORM_ll
!
!!    Reference
!!    ---------
!
!     User Interface for Meso-NH parallel package
!     Ph. Kloos, L. Giraud, R. Guivarch, D. Lugato
!
!!    Authors
!!    -------
!
!     R. Guivarch, D. Lugato    * CERFACS *
!     Ph. Kloos                 * CERFACS - CNRM *
!
!     Modif
!     Juan/Didier 12/03/2009: array bound bug correction with 1proc/MPIVIDE
!     J. Escobar  27/06/2011  correction for gridnesting with different SHAPE 
! 
USE MODD_MPIF
use modd_precision, only: MNHINT_MPI, MNHREAL_MPI
USE MODD_STRUCTURE_ll
USE MODD_VAR_ll, ONLY : NMNH_COMM_WORLD

use mode_msg

implicit none

interface GET_GLOBALSLICE_ll
  module procedure GET_1DGLOBALSLICE_ll, GET_2DGLOBALSLICE_ll
end interface

interface GET_SLICE_ll
  module procedure GET_1DSLICE_ll, GET_2DSLICE_ll
end interface

CONTAINS

  SUBROUTINE SLIDE_COORD(KDIM_DATA,KDIM_PROC,THIS_PROC,KOR,KEND)

    !!    Purpose
    !
    !   Compute for the processor=THIS_PROC the origine/end of slide in decomposing
    !   an array of data of dimension=KDIM_DATA on KDIM_PROC
    !
    !!    Author
    !!    ------
    !     J. ESCOBAR    * LA *

    IMPLICIT NONE
    !
    !*       0.1   declarations of arguments
    !
    INTEGER, INTENT(IN)  :: KDIM_DATA ! dimension of data to split
    INTEGER, INTENT(IN)  :: KDIM_PROC ! numbers of processor to use in splitting
    INTEGER, INTENT(IN)  :: THIS_PROC ! processor id from 1..NB_PROC
    INTEGER, INTENT(OUT) :: KOR,KEND  ! Origine/End coordonate
    !
    !*       0.2   declarations of local variables
    !
    INTEGER             :: IDIM_SLIDE ! slide dimension ( without rest/delta )
    INTEGER             :: IREST      ! number of point in surabondance to distribut 
    INTEGER             :: IDELTAOR,IDELTAEND     ! offset in origine to apply 

    IDIM_SLIDE   = KDIM_DATA/KDIM_PROC
    IREST        = MOD(KDIM_DATA,KDIM_PROC)
    IDELTAOR     = MIN(IREST,THIS_PROC-1)
    IDELTAEND    = MIN(IREST,THIS_PROC)

    KOR   = ( THIS_PROC - 1 ) * IDIM_SLIDE + 1 + IDELTAOR
    KEND  =   THIS_PROC       * IDIM_SLIDE     + IDELTAEND

  END SUBROUTINE SLIDE_COORD

!
!     ########################################
      LOGICAL FUNCTION LNORTH_ll(K,HSPLITTING)
!     ########################################
!
!!****  *LNORTH_ll* - function which returns the position on to the boundaries
!                  of the subdomain K according to the splitting
! 
!!    Purpose
!!    -------
!     the Purpose of this routine is to offer a transparent way to obtain
!     the position of a subdomain
!
!!**  Method
!!    ------
!     if the argument HSPLITTING is omitted the 2Way splitting is considered
!     if the argument K is omitted, the local subdomain is considered
! 
!!    External
!!    --------
! 
!!    Implicit Arguments
!!    ------------------
!     Module MODD_VAR_ll
!        TCRRT_PROCONF - Current configuration for current model
!        IP - Number of local processor=subdomain
!        NPROC - Number of processors
! 
!!    Reference
!!    ---------
! 
!!    Author
!!    ------
!     R. Guivarch
! 
!!    Modifications
!!    -------------
!     Original 01/05/98
! 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_VAR_ll, ONLY : TCRRT_PROCONF, IP, NPROC
!
  IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
  INTEGER, INTENT(IN), OPTIONAL :: K ! number of the subdomain
  CHARACTER(len=1), INTENT(IN), OPTIONAL :: HSPLITTING ! kind of splitting
!
!*       0.2   declarations of local variables
!
  INTEGER :: IT ! number of the tested subdomain
!
!-------------------------------------------------------------------------------
!
!*       1.    COMPUTATION OF THE RESULT :
!              -------------------------
!
  IF( PRESENT(K) ) THEN
     IT = K
  ELSE
     IT = IP
  ENDIF
!
  LNORTH_ll = .FALSE.
!
  IF(.NOT.PRESENT(HSPLITTING)) THEN
    LNORTH_ll = TCRRT_PROCONF%TBOUND(IT)%NORTH
  ELSEIF(HSPLITTING .EQ. 'B') THEN
    LNORTH_ll = TCRRT_PROCONF%TBOUND(IT)%NORTH
  ELSEIF(HSPLITTING .EQ. 'X') THEN
    LNORTH_ll = (IT.EQ.NPROC)
  ELSEIF(HSPLITTING .EQ. 'Y') THEN
    LNORTH_ll = .TRUE.
  ENDIF
!
!-------------------------------------------------------------------------------
!
      END FUNCTION LNORTH_ll
!
!     #######################################
      LOGICAL FUNCTION LWEST_ll(K,HSPLITTING)
!     #######################################
!
!!****  *LWEST_ll* - function which returns the position on to the boundaries
!                  of the subdomain K according to the splitting
! 
!!    Purpose
!!    -------
!     the Purpose of this routine is to offer a transparent way to obtain
!     the position of a subdomain
!
!!**  Method
!!    ------
!     if the argument HSPLITTING is omitted the 2Way splitting is considered
!     if the argument K is omitted, the local subdomain is considered
! 
!!    External
!!    --------
! 
!!    Implicit Arguments
!!    ------------------
!     Module MODD_VAR_ll
!        TCRRT_PROCONF - Current configuration for current model
!        IP - Number of local processor=subdomain
! 
!!    Reference
!!    ---------
! 
!!    Author
!!    ------
!     R. Guivarch
! 
!!    Modifications
!!    -------------
!     Original 01/05/98
! 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_VAR_ll, ONLY : TCRRT_PROCONF, IP
!
  IMPLICIT NONE
!
!
!*       0.1   declarations of arguments
!
  INTEGER, INTENT(IN), OPTIONAL :: K ! number of the subdomain
  CHARACTER(len=1), INTENT(IN), OPTIONAL :: HSPLITTING ! kind of splitting

!!
!*       0.2   declarations of local variables
!
  INTEGER :: IT ! number of the tested subdomain
!
!-------------------------------------------------------------------------------
!
!*       1.    COMPUTATION OF THE RESULT :
!              -------------------------
!
  IF( PRESENT(K) ) THEN
     IT = K
  ELSE
     IT = IP
  ENDIF
!
  LWEST_ll = .FALSE.
!
  IF(.NOT.PRESENT(HSPLITTING)) THEN
    LWEST_ll = TCRRT_PROCONF%TBOUND(IT)%WEST
  ELSEIF(HSPLITTING .EQ. 'B') THEN
    LWEST_ll = TCRRT_PROCONF%TBOUND(IT)%WEST
  ELSEIF(HSPLITTING .EQ. 'X') THEN
    LWEST_ll = .TRUE.
  ELSEIF(HSPLITTING .EQ. 'Y') THEN
    LWEST_ll = (IT.EQ.1)
  ENDIF
!
!-------------------------------------------------------------------------------
!
      END FUNCTION LWEST_ll
!
!     ########################################
      LOGICAL FUNCTION LSOUTH_ll(K,HSPLITTING)
!     ########################################
!
!!****  *LSOUTH_ll* - function which returns the position on to the boundaries
!!                 of the subdomain K according to the splitting
!!
!!    Purpose
!!    -------
!     the Purpose of this routine is to offer a transparent way to obtain
!     the position of a subdomain
!
!!**  Method
!!    ------
!     if the argument HSPLITTING is omitted the 2Way splitting is considered
!     if the argument K is omitted, the local subdomain is considered
! 
!!    External
!!    --------
! 
!!    Implicit Arguments
!!    ------------------
!     Module MODD_VAR_ll
!        TCRRT_PROCONF - Current configuration for current model
!        IP - Number of local processor=subdomain
! 
!!    Reference
!!    ---------
! 
!!    Author
!!    ------
!     R. Guivarch
! 
!!    Modifications
!!    -------------
!     Original 01/05/98
! 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_VAR_ll, ONLY : TCRRT_PROCONF, IP
!
  IMPLICIT NONE
!
!
!*       0.1   declarations of arguments
!
  INTEGER, INTENT(IN), OPTIONAL :: K ! number of the subdomain
  CHARACTER(len=1), INTENT(IN), OPTIONAL :: HSPLITTING ! kind of splitting
 
!!
!*       0.2   declarations of local variables
!
  INTEGER :: IT ! number of the tested subdomain
!
!-------------------------------------------------------------------------------
!
!*       1.    COMPUTATION OF THE RESULT :
!              -------------------------
!
  IF( PRESENT(K) ) THEN
     IT = K
  ELSE
     IT = IP
  ENDIF
!
  LSOUTH_ll = .FALSE.
!
  IF(.NOT.PRESENT(HSPLITTING)) THEN
    LSOUTH_ll = TCRRT_PROCONF%TBOUND(IT)%SOUTH
  ELSEIF(HSPLITTING .EQ. 'B') THEN
    LSOUTH_ll = TCRRT_PROCONF%TBOUND(IT)%SOUTH
  ELSEIF(HSPLITTING .EQ. 'X') THEN
    LSOUTH_ll = (IT.EQ.1)
  ELSEIF(HSPLITTING .EQ. 'Y') THEN
     LSOUTH_ll = .TRUE.
  ENDIF
!
!-------------------------------------------------------------------------------
!
      END FUNCTION LSOUTH_ll
!
!     #######################################
      LOGICAL FUNCTION LEAST_ll(K,HSPLITTING)
!     #######################################
!
!!****  *LEAST_ll* - function which returns the position on to the boundaries
!!                 of the subdomain K according to the splitting
!!
!!    Purpose
!!    -------
!     the Purpose of this routine is to offer a transparent way to obtain
!     the position of a subdomain
!
!!**  Method
!!    ------
!     if the argument HSPLITTING is omitted the 2Way splitting is considered
!     if the argument K is omitted, the local subdomain is considered
! 
!!    External
!!    --------
! 
!!    Implicit Arguments
!!    ------------------
!     Module MODD_VAR_ll
!        TCRRT_PROCONF - Current configuration for current model
!        IP - Number of local processor=subdomain
!        NPROC - Number of processors
! 
!!    Reference
!!    ---------
! 
!!    Author
!!    ------
!     R. Guivarch
! 
!!    Modifications
!!    -------------
!     Original 01/05/98
! 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_VAR_ll, ONLY : TCRRT_PROCONF, IP, NPROC
!
  IMPLICIT NONE
!
!
!*       0.1   declarations of arguments
!
  INTEGER, INTENT(IN), OPTIONAL :: K ! number of the subdomain
  CHARACTER(len=1), INTENT(IN), OPTIONAL :: HSPLITTING ! kind of splitting
 
!!
!*       0.2   declarations of local variables
!
  INTEGER :: IT ! number of the tested subdomain
!
!-------------------------------------------------------------------------------
!
!*       1.    COMPUTATION OF THE RESULT :
!              -------------------------
!
  IF( PRESENT(K) ) THEN
     IT = K
  ELSE
     IT = IP
  ENDIF
!
  LEAST_ll = .FALSE.
!
  IF(.NOT.PRESENT(HSPLITTING)) THEN
    LEAST_ll = TCRRT_PROCONF%TBOUND(IT)%EAST
  ELSEIF(HSPLITTING .EQ. 'B') THEN
    LEAST_ll = TCRRT_PROCONF%TBOUND(IT)%EAST
  ELSEIF(HSPLITTING .EQ. 'X') THEN
    LEAST_ll = .TRUE.
  ELSEIF(HSPLITTING .EQ. 'Y') THEN
    LEAST_ll = (IT.EQ.NPROC)
  ENDIF
!
!-------------------------------------------------------------------------------
!
      END FUNCTION LEAST_ll
!
!     #################################################
      SUBROUTINE GET_DIM_EXT_ll( HSPLIT, KXDIM, KYDIM )
!     #################################################
!
!!****  *GET_DIM_EXT_ll* - returns the dimensions of the extended 2way subdomain
!                   or of the x-slices subdomain or of the y-slices
!                   subdomain of the local processor
! 
!!    Purpose
!!    -------
!     the Purpose of this routine is to give subdomain dimension
!
!!**  Method
!!    ------
!     if HSPLIT='B', the dimensions of the extended 2way subdomain are returned
!     if HSPLIT='X', the dimensions of x-slices subdomain are returned
!     if HSPLIT='Y', the dimensions of y-slices subdomain are returned
! 
!!    External
!!    --------
! 
!!    Implicit Arguments
!!    ------------------
!     Module MODD_VAR_ll
!        TCRRT_COMDATA - Current communication data structure for current model
!                        and local processor
! 
!!    Reference
!!    ---------
! 
!!    Author
!!    ------
!     R. Guivarch
! 
!!    Modifications
!!    -------------
!     Original 01/05/98
! 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_VAR_ll, ONLY : TCRRT_COMDATA
!
  IMPLICIT NONE
!
!
!*       0.1   declarations of arguments
!
  CHARACTER(len=1), INTENT(IN) :: HSPLIT
!
  INTEGER, INTENT(OUT) :: KXDIM, KYDIM
!
!
!*       0.2   declarations of local variables
!
!-------------------------------------------------------------------------------
!
!*       1.    Return the dimensions
!
  IF( HSPLIT .EQ. 'B' ) THEN
    KXDIM = TCRRT_COMDATA%TSPLIT_B%NDIMXE
    KYDIM = TCRRT_COMDATA%TSPLIT_B%NDIMYE
  ELSEIF ( HSPLIT .EQ. 'X' ) THEN
    KXDIM = TCRRT_COMDATA%TSPLIT_X%NDIMXP
    KYDIM = TCRRT_COMDATA%TSPLIT_X%NDIMYP
  ELSE
    KXDIM = TCRRT_COMDATA%TSPLIT_Y%NDIMXP
    KYDIM = TCRRT_COMDATA%TSPLIT_Y%NDIMYP
  ENDIF
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE GET_DIM_EXT_ll
!
!     ##################################################
      SUBROUTINE GET_DIM_PHYS_ll( HSPLIT, KXDIM, KYDIM )
!     ##################################################
!
!!****  *GET_DIM_PHYS_ll* - returns the dimensions of the physical
!                           2way subdomain or of the x-slices subdomain or
!                           of the y-slices subdomain of the local processor
! 
!!    Purpose
!!    -------
!     the Purpose of this routine is to give subdomain dimension
!
!!**  Method
!!    ------
!     if HSPLIT='B', the dimensions of the physical 2way subdomain are returned
!     if HSPLIT='X', the dimensions of x-slices subdomain are returned
!     if HSPLIT='Y', the dimensions of y-slices subdomain are returned
! 
!!    External
!!    --------
! 
!!    Implicit Arguments
!!    ------------------
!     Module MODD_VAR_ll
!        TCRRT_COMDATA - Current communication data structure for current model
!                        and local processor
! 
!!    Reference
!!    ---------
! 
!!    Author
!!    ------
!     R. Guivarch
! 
!!    Modifications
!!    -------------
!     Original 01/05/98
! 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_VAR_ll, ONLY : TCRRT_COMDATA
!
  IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
  CHARACTER(len=1), INTENT(IN) :: HSPLIT
!
  INTEGER, INTENT(OUT) :: KXDIM, KYDIM
!
!*       0.2   declarations of local variables
!
!-------------------------------------------------------------------------------
!
!*       1.    Return the dimensions
!
  IF( HSPLIT .EQ. 'B' ) THEN
    KXDIM = TCRRT_COMDATA%TSPLIT_B%NDIMXP
    KYDIM = TCRRT_COMDATA%TSPLIT_B%NDIMYP
  ELSEIF ( HSPLIT .EQ. 'X' ) THEN
    KXDIM = TCRRT_COMDATA%TSPLIT_X%NDIMXP
    KYDIM = TCRRT_COMDATA%TSPLIT_X%NDIMYP
  ELSE
    KXDIM = TCRRT_COMDATA%TSPLIT_Y%NDIMXP
    KYDIM = TCRRT_COMDATA%TSPLIT_Y%NDIMYP
  ENDIF
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE GET_DIM_PHYS_ll
!
!     ##########################################
      SUBROUTINE GET_OR_ll( HSPLIT, KXOR, KYOR )
!     ##########################################
!
!!****  *GET_OR_ll* - returns the origin'coordinates of the extended
!                     2way subdomain or of the x-slices subdomain
!                     or of the y-slices
!                     subdomain of the local processor (global indices)
! 
!!    Purpose
!!    -------
!     the Purpose of this routine is to give subdomain origin
!
!!**  Method
!!    ------
!     if HSPLIT = 'B', the origin of the extended subdomain are returned
!     if HSPLIT = 'X', the origin of x-slices subdomain are returned
!     if HSPLIT = 'Y', the origin of y-slices subdomain are returned
! 
!!    External
!!    --------
!     
!!    Implicit Arguments
!!    ------------------
!     Module MODD_VAR_ll
!        TCRRT_COMDATA - Current communication data structure for current model
!                        and local processor
!     
!!    Reference
!!    ---------
!     
!!    Author
!!    ------
!     R. Guivarch
!     
!!    Modifications
!!    -------------
!     Original 01/05/98
! 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_VAR_ll, ONLY : TCRRT_COMDATA
!
  IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
  CHARACTER(len=1), INTENT(IN) :: HSPLIT
!
  INTEGER, INTENT(OUT) :: KXOR, KYOR
!
!*       0.2   declarations of local variables
!
!-------------------------------------------------------------------------------
!
!*       1.    
!
  IF( HSPLIT .EQ. 'B' ) THEN
    KXOR = TCRRT_COMDATA%TSPLIT_B%NXORE
    KYOR = TCRRT_COMDATA%TSPLIT_B%NYORE
  ELSEIF ( HSPLIT .EQ. 'X' ) THEN
    KXOR = TCRRT_COMDATA%TSPLIT_X%NXORP
    KYOR = TCRRT_COMDATA%TSPLIT_X%NYORP
  ELSE
    KXOR = TCRRT_COMDATA%TSPLIT_Y%NXORP
    KYOR = TCRRT_COMDATA%TSPLIT_Y%NYORP
  ENDIF
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE GET_OR_ll
!
!     ####################################################
      SUBROUTINE GET_INDICE_ll( KXOR, KYOR, KXEND, KYEND, KSIZE1, KSIZE2 )
!     ####################################################
!
!!****  *GET_INDICE_ll* - returns the origin's coordinates and the end's
!                           coordinates of the local physical
!                           subdomain (in local indices)
! 
!!    Purpose
!!    -------
!
!!**  Method
!!    ------
! 
!!    External
!!    --------
!     Module MODE_TOOLS_ll
!       LWEST_ll, LSOUTH_ll, LEAST_ll, LNORTH_ll
! 
!!    Implicit Arguments
!!    ------------------
!     Module MODD_VAR_ll
!        TCRRT_COMDATA - Current communication data structure for current model
!                        and local processor
!        JPHALO- halo size
! 
!     Module MODD_PARAMETERS_ll
!        JPHEXT - halo size
! 
!!    Reference
!!    ---------
! 
!!    Author
!!    ------
!     R. Guivarch
! 
!!    Modifications
!!    -------------
!     Original 08/07/98
! 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_PARAMETERS_ll, ONLY : JPHEXT
  USE MODD_VAR_ll, ONLY : TCRRT_COMDATA, JPHALO
!
  IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
  INTEGER, INTENT(OUT) :: KXOR, KYOR, KXEND, KYEND
  INTEGER, INTENT(IN),OPTIONAL  :: KSIZE1, KSIZE2
!
!*       0.2   declarations of local variables
!
!-------------------------------------------------------------------------------
!
  IF(LWEST_ll()) THEN
    KXOR = 1 + JPHEXT
  ELSE
    KXOR = 1 + JPHALO
  ENDIF
!
  IF(LSOUTH_ll()) THEN
    KYOR = 1 + JPHEXT
  ELSE
    KYOR = 1 + JPHALO
  ENDIF
!
  IF(LEAST_ll()) THEN
    KXEND = TCRRT_COMDATA%TSPLIT_B%NDIMXE - JPHEXT
  ELSE
    KXEND = TCRRT_COMDATA%TSPLIT_B%NDIMXE - JPHALO
  ENDIF
!
  IF(LNORTH_ll()) THEN
    KYEND = TCRRT_COMDATA%TSPLIT_B%NDIMYE - JPHEXT
  ELSE
    KYEND = TCRRT_COMDATA%TSPLIT_B%NDIMYE - JPHALO
  ENDIF
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE GET_INDICE_ll
!
!     ##########################################
      SUBROUTINE GET_GLOBALDIMS_ll(KIMAX, KJMAX, KMODEL)
!     ##########################################
!
!!****  *GET_GLOBALDIMS_ll* - returns the global horizontal dimensions
!                             of the current model (External halo excluded)
! 
!!    Purpose
!!    -------
!
!!**  Method
!!    ------
! 
!!    External
!!    --------
! 
!!    Implicit Arguments
!!    ------------------
!     Module MODD_VAR_ll
!        TCRRT_PROCONF - Current configuration for current model
!
!     Module MODD_DIM_ll
!       NDXRATIO_ALL, NDYRATIO_ALL, NXOR_ALL, NYOR_ALL,
!       NXEND_ALL, NYEND_ALL
! 
!!    Reference
!!    ---------
! 
!!    Author
!!    ------
!     P. Kloos (CERFACS)
! 
!     Modifications
!!    -------------
!     Original 15 september 1998
! 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_DIM_ll, ONLY : NDXRATIO_ALL, NDYRATIO_ALL, &
                          NXOR_ALL, NYOR_ALL,         &
                          NXEND_ALL, NYEND_ALL
  USE MODD_VAR_ll, ONLY : TCRRT_PROCONF

  USE MODD_PARAMETERS_ll, ONLY : JPHEXT
!
  IMPLICIT NONE
!
!
!*       0.1   declarations of arguments
!
  INTEGER, INTENT(OUT) :: KIMAX, KJMAX ! current model dimensions
  INTEGER, OPTIONAL, INTENT(IN) :: KMODEL  ! number of the current model
!
!*       0.2   declarations of local variables
!
  INTEGER :: IMODEL ! number of the current model
!
!-------------------------------------------------------------------------------
!
!*       1.    Extract the number of the current model.
!
IF ( PRESENT(KMODEL) ) THEN
  IMODEL = KMODEL
ELSE
  IMODEL = TCRRT_PROCONF%NUMBER
ENDIF
!
!*       2.    Compute the dimensions of the model
!
  KIMAX = NDXRATIO_ALL(IMODEL) * (NXEND_ALL(IMODEL)-NXOR_ALL(IMODEL) -2*JPHEXT + 1)
  KJMAX = NDYRATIO_ALL(IMODEL) * (NYEND_ALL(IMODEL)-NYOR_ALL(IMODEL) -2*JPHEXT + 1)
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE GET_GLOBALDIMS_ll
! 
!     ######################################################
      SUBROUTINE GET_PHYSICAL_ll( KXOR, KYOR, KXEND, KYEND )
!     ######################################################
!
!!****  *GET_PHYSICAL_ll* - returns the origin's coordinates and the end's
!                           coordinates of the intersection
!                           of the physical global domain with the local
!                           extended subdomain (in local indices)
!
!!    Purpose
!!    -------
!
!!**  Method
!!    ------
!
!!    External
!!    --------
!     Module MODE_TOOLS_ll
!       LWEST_ll, LSOUTH_ll, LEAST_ll, LNORTH_ll
! 
!!    Implicit Arguments
!!    ------------------
!     Module MODD_VAR_ll
!        TCRRT_COMDATA - Current communication data structure for current model
!                        and local processor
!
!     Module MODD_PARAMETERS_ll
!       JPHEXT
! 
!!    Reference
!!    ---------
! 
!!    Author
!!    ------
!     R. Guivarch
! 
!!    Modifications
!!    -------------
!     Original 01/05/98
! 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_PARAMETERS_ll, ONLY : JPHEXT
  USE MODD_VAR_ll, ONLY : TCRRT_COMDATA
!
  IMPLICIT NONE
!
!
!*       0.1   declarations of arguments
!
  INTEGER, INTENT(OUT) :: KXOR, KYOR, KXEND, KYEND
!
!
!*       0.2   declarations of local variables
!
!-------------------------------------------------------------------------------
!
  KXOR = 1
  IF(LWEST_ll()) KXOR = KXOR + JPHEXT
!
  KYOR = 1
  IF(LSOUTH_ll()) KYOR = KYOR + JPHEXT
!
  KXEND = TCRRT_COMDATA%TSPLIT_B%NDIMXE
  IF(LEAST_ll()) KXEND = KXEND - JPHEXT
!
  KYEND = TCRRT_COMDATA%TSPLIT_B%NDIMYE
  IF(LNORTH_ll()) KYEND = KYEND - JPHEXT
!
!
      END SUBROUTINE GET_PHYSICAL_ll
!
!     ###################################################################
      SUBROUTINE GET_INTERSECTION_ll( KXOR,KYOR,KXEND,KYEND,KXORI,KYORI,&
                                       KXENDI,KYENDI,HDOM,KINFO,KIP)
!     ###################################################################
!
!!****  *GET_INTERSECTION_ll* - routine to get the indices of the intersection
!                               between a geographic region and the EXTENDED or
!                               PHYSICAL subdomain of the KIP or current
!                               processor.
!                               The input indices are global.
!                               The output indices are local.
!                               If the returned indices are null
!                               the intersection is void.
! 
!!    Purpose
!!    -------
!
!!**  Method
!!    ------
!     The processor computes the intersection of a sub-domain
!     with the geographic region.
! 
!!    External
!!    --------
! 
!!    Implicit Arguments
!!    ------------------
!     Module MODD_VAR_ll 
!       TCRRT_PROCONF - Current configuration for current model
!       IP - Number of the local processor
!
!     Module MODD_PARAMETERS_ll
!       JPHEXT, JPVEXT - halo sizes
!
!     Module MODD_STRUCTURE_ll
!       type MODELSPLITTING_ll
!
!!    Reference
!!    ---------
! 
!!    Author
!!    ------
!     R. Guivarch, D. Gazen
! 
!!    Modifications
!!    -------------
!     Original 20/01/99
! 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_STRUCTURE_ll, ONLY : MODELSPLITTING_ll
  USE MODD_PARAMETERS_ll, ONLY : JPHEXT, JPVEXT
  USE MODD_VAR_ll, ONLY : TCRRT_PROCONF, IP
!
  IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
  CHARACTER(LEN=4), INTENT(IN) :: HDOM     ! 'EXTE' for extended subdomain
                                            ! 'PHYS' for physical subdomain
  INTEGER, INTENT(IN)  :: KXOR, KYOR, KXEND, KYEND   ! Coordinates of the region
  INTEGER, INTENT(OUT) :: KXORI, KYORI, KXENDI, KYENDI ! Global Coordinates of
                                                       ! the intersection
  INTEGER, INTENT(OUT) :: KINFO             ! Returned Info
  INTEGER, INTENT(IN), OPTIONAL:: KIP       ! Processor number
                                            ! (or subdomain number)
!
!*       0.2   declarations of local variables
!
  INTEGER :: IIMAX_ll, IJMAX_ll           ! size of the total physical model
  INTEGER :: IXORI, IYORI, IXENDI, IYENDI ! Local coordinates of intersection
  TYPE(MODELSPLITTING_ll), POINTER :: TZSPLIT ! 2way-splitting
  INTEGER :: IXOR, IYOR, IXEND, IYEND     ! global coordinate of KIP subdomain
  INTEGER :: IIP ! subdomain or processor number
!
!-------------------------------------------------------------------------------
!
  KINFO = 0
  IF (PRESENT(KIP)) THEN
    IIP = KIP
  ELSE
    IIP = IP
  END IF
!
!-------------------------------------------------------------------------------
!
!*       1.    INITIALISE BOUNDARY VALUES
!              --------------------------
!
!        1.1   'EXTE'nded sub-domain

  TZSPLIT => TCRRT_PROCONF%TSPLITS_B(IIP)
!
  IXOR  = TZSPLIT%NXORE
  IYOR  = TZSPLIT%NYORE
  IXEND = TZSPLIT%NXENDE
  IYEND = TZSPLIT%NYENDE
!
!        1.2  'PHYS'ical sub-domain
!
  IF (HDOM=='PHYS') THEN
    IF (.NOT. LWEST_ll(IIP))  IXOR  = TZSPLIT%NXORP
    IF (.NOT. LEAST_ll(IIP))  IXEND = TZSPLIT%NXENDP
    IF (.NOT. LSOUTH_ll(IIP)) IYOR  = TZSPLIT%NYORP
    IF (.NOT. LNORTH_ll(IIP)) IYEND = TZSPLIT%NYENDP
  END IF
!
!-------------------------------------------------------------------------------
!
!*       2.    THE COORDINATES ARE NOT IN THE PHYSICAL DOMAIN -> ERROR :
!              -------------------------------------------------------
!
  CALL GET_GLOBALDIMS_ll(IIMAX_ll, IJMAX_ll)
!
  IF((KXOR < 1 ) .OR. (KYOR < 1 ) .OR. &
     (KXEND > IIMAX_ll + 2*JPHEXT) .OR. (KYEND > IJMAX_ll + 2*JPHEXT)) THEN
!
!   Error
!
    KINFO = -1
!
    RETURN
!
!-------------------------------------------------------------------------------
!
!*       3.    THE COORDINATES ARE IN THE PHYSICAL DOMAIN :
!              ------------------------------------------
!
  ELSE
!
    IXORI = MAX( IXOR, KXOR )
    IYORI = MAX( IYOR, KYOR )
!
    IXENDI = MIN( IXEND, KXEND )
    IYENDI = MIN( IYEND, KYEND )
!
!*       3.1   The intersection is empty
!
    IF((IXORI > IXENDI) .OR. (IYORI > IYENDI)) THEN
!
      KXORI = 0
      KYORI = 0
!
      KXENDI = 0
      KYENDI = 0
!
      KINFO = 1
!
    ELSE
!
!*       3.2   Switch to local coordinates
!
      KXORI = IXORI - TZSPLIT%NXORE + 1
      KYORI = IYORI - TZSPLIT%NYORE + 1
      KXENDI = IXENDI - TZSPLIT%NXORE + 1
      KYENDI = IYENDI - TZSPLIT%NYORE + 1
    ENDIF
!
  ENDIF
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE GET_INTERSECTION_ll
!
!      ##################################################################
!       SUBROUTINE GET_GLOBALSLICE_ll(PARRAY, HDIR, KLOC, PGLOBALSLICE, &
!                                     KB, KE, KERR)
!      ##################################################################
!
!!     Purpose
!!     -------  
!         The Purpose of this routine is to get a slice of the
!      domain in the x or y direction.
!
!      ###################################################################
       SUBROUTINE GET_1DGLOBALSLICE_ll(PARRAY, HDIR, KLOC, PGLOBALSLICE, &
                                       KB, KE, KERR)
!      ###################################################################
!
!!     Purpose
!!     -------  
!         The Purpose of this routine is to extract a slice of an horizontal
!      field PARRAY along the x or y direction
!
!!**   Method
!!     ------
! 
!         An MPI communicator with the processes corresponding to the
!      subdomains intersecting with the slice is built. This
!      communicator is then used to gather the whole slice (i.e. 
!      the global slice) on these procs. The global slice is then
!      broadcasted on all procs that are not on the slice.
! 
!!     External
!!     --------
!      Module MODE_TOOLS_ll
!        LWEST_ll, LSOUTH_ll, LNORTH_ll, LEAST_ll
!        GET_GLOBALDIMS_ll
!
!!    Implicit Arguments
!!    ------------------
!      Module MODD_STRUCTURE_ll
!        type MODELSPLITTING_ll
!
!      Module MODD_VAR_ll
!        NPROC - Number of processors
!        TCRRT_PROCONF -  Current configuration for current model
!        IP - Number of the local processor
!
!      Module MODD_PARAMETERS_ll
!        JPHEXT - halo size
!
!!    Reference
!!    ---------
!
!     User Interface for the MesoNH Parallel Package
! 
!!    Author
!!    ------
!     P. Kloos (CERFACS)
! 
!!    Modifications
!!    -------------
!     Original 14 August 1998
! 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_PARAMETERS_ll, ONLY: JPHEXT
  USE MODD_STRUCTURE_ll,  ONLY: MODELSPLITTING_ll
  USE MODD_VAR_ll,        ONLY: NPROC, TCRRT_PROCONF, IP
!
  IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
  REAL, DIMENSION(:,:), TARGET, INTENT(IN) :: PARRAY ! horizontal field
  CHARACTER(LEN=1), INTENT(IN)     :: HDIR   ! direction ("X" or "Y")
  INTEGER, INTENT(IN)              :: KLOC   ! coordinate of the slice
                                             ! to extract 
                                             ! (in global coordinates)
  REAL, DIMENSION(:), INTENT(OUT)  :: PGLOBALSLICE ! output slice
  INTEGER, OPTIONAL    :: KB, KE             ! begin and end positions of the
                                             ! extracted slices
                                             ! (global coordinates)
  INTEGER, OPTIONAL    :: KERR               ! error code
!
!*       0.2   declarations of local variables
!
  TYPE(MODELSPLITTING_ll), POINTER :: TZSPLIT
  INTEGER :: IDIM, IDIM2
  INTEGER :: INUMPROC, IERR
  INTEGER :: ILOC   ! local location of the slice
  INTEGER :: ICOUNT ! number of relevant procs (= those which are on the slice)
  INTEGER :: ISIZE  ! length of the local slice
  INTEGER, DIMENSION(:), ALLOCATABLE :: IPROCS ! array of procs that are on 
                                               ! the slice 
  INTEGER, DIMENSION(:), ALLOCATABLE :: ISIZES ! length of the local slice
                                               ! on all procs 
  INTEGER, DIMENSION(:), ALLOCATABLE :: IDISPL ! array of locations of the 
                                               ! procs in the slice
                                               ! (for MPI_ALLGATHERV)
  INTEGER :: IWRLD_GROUP                ! world group
  INTEGER :: IGROUP_GLOBALSLICE         ! group for the proc on the slice
  INTEGER :: ICOMM_GLOBALSLICE          ! communicator for the proc on the slice
  INTEGER :: IGROUP, ICOMM
  INTEGER, DIMENSION(NPROC) :: IGLOBALSLICEPROC
  INTEGER, DIMENSION(2) :: IOR, IEND, IORP, IENDP, IORE, IENDE  ! splitting
  INTEGER :: IWEST, IEAST, INORTH, ISOUTH
  INTEGER :: IB, IE ! beginning and end of the local slice
  INTEGER :: IDISPL1  ! beginning of the slice (global)
  INTEGER :: J          ! loop index
  INTEGER :: IGLOBALSLICELENGTH ! length of the global slice
  INTEGER, DIMENSION(2) :: IMAX ! maximum dimensions
  REAL, DIMENSION(:), POINTER :: ZPTR
!
!-------------------------------------------------------------------------------
!
!*       1.    INITIALISATIONS
!              ---------------
  IWEST=0; IEAST=0; INORTH=0; ISOUTH=0
  IB=0; IE=0
  IDISPL1=0
!
!*       1.1   Get current splitting
!
  IF (LWEST_ll())  IWEST=-1
  IF (LEAST_ll())  IEAST=1
  IF (LNORTH_ll()) INORTH=1
  IF (LSOUTH_ll()) ISOUTH=-1
  TZSPLIT => TCRRT_PROCONF%TSPLITS_B(IP)
  IOR(1)  = TZSPLIT%NXORP+IWEST
  IOR(2)  = TZSPLIT%NYORP+ISOUTH 
  IEND(1) = TZSPLIT%NXENDP+IEAST
  IEND(2) = TZSPLIT%NYENDP+INORTH
!
  IORP(1)  = TZSPLIT%NXORP
  IORP(2)  = TZSPLIT%NYORP
  IENDP(1) = TZSPLIT%NXENDP
  IENDP(2) = TZSPLIT%NYENDP
!
  IORE(1)  = TZSPLIT%NXORE
  IORE(2)  = TZSPLIT%NYORE
  IENDE(1) = TZSPLIT%NXENDE
  IENDE(2) = TZSPLIT%NYENDE
!
  IGLOBALSLICEPROC = 0
  CALL GET_GLOBALDIMS_ll(IMAX(1),IMAX(2))
  IMAX(:) = IMAX(:) + 2*JPHEXT
!
!*       1.2   Set dimension (1 for X and 2 for Y)
!
  IDIM  = IACHAR(HDIR)-IACHAR('X')+1
  IDIM2 = 3-IDIM ! IDIM's inverse 
!
!*       1.4   Set beginning and end of local slice
!
  IF (PRESENT(KB) .AND. PRESENT(KE)) THEN
!
!*      Test the ranges
!
    IF((KB < 1 ) .OR. (KE > IMAX(IDIM))) THEN
!
!   Error
!
      KERR = -1
      RETURN
!
    ENDIF
!
    IB = MAX(KB, IOR(IDIM))  - IORE(IDIM) + 1
    IE = MIN(KE, IEND(IDIM)) - IORE(IDIM) + 1
!    IDISPL1 = KB-1  ! old version
    IDISPL1 = 0
  ELSE ! default : physical domain
    IB = 1+JPHEXT
    IE = IENDP(IDIM)-IORE(IDIM)+1
!    IDISPL1 = JPHEXT ! old version
    IDISPL1 = 0
  ENDIF
!
!-------------------------------------------------------------------------------
!
!*       2.    CREATE MPI COMMUNICATOR WITH THE PROCS ON THE SLICE
!              ---------------------------------------------------
!
!*       2.1   Test if i am on the slice
!              if so, INUMPROC = my MPI rank
!              if not, INUMPROC = MPI_PROC_NULL
!
  IF (KLOC >= IOR(IDIM2) .AND. KLOC <= IEND(IDIM2) .AND. IB<=IE) THEN
!
! Set local location
!
    ILOC = KLOC-IORE(IDIM2)+1
!
! Set relevant procs
!
    INUMPROC = IP-1
!
! Set lenght of the local slice
!
    ISIZE = IE - IB + 1
!
!*        2.2   Have ZPTR point to the slice
!
    SELECT CASE(HDIR)
    CASE("X")
      ZPTR => PARRAY(IB:IE,ILOC)
!
    CASE("Y")
      ZPTR => PARRAY(ILOC,IB:IE)
!
    CASE DEFAULT
      call Print_msg( NVERB_FATAL, 'GEN', 'GET_1DGLOBALSLICE_ll', 'invalid HDIR dummy argument ('//hdir//')' )
!
    END SELECT 
!
  ELSE
!
    INUMPROC = MPI_PROC_NULL
!
  ENDIF
!
!*        2.3    Gather values of INUMPROC 
!
  CALL MPI_ALLGATHER( (/ INUMPROC /) , 1, MNHINT_MPI, IGLOBALSLICEPROC, 1, &
                     MNHINT_MPI, NMNH_COMM_WORLD, IERR)
!
!*        2.4     Get MPI world group
!
  CALL MPI_COMM_GROUP(NMNH_COMM_WORLD, IWRLD_GROUP, IERR)
!
!*        2.5     Count number of proc that contain the slice
!
  ICOUNT = COUNT(IGLOBALSLICEPROC.NE.MPI_PROC_NULL)
!
!*        2.6     Create MPI group with the procs that contain the slice
!
  ALLOCATE(IPROCS(ICOUNT))
  IPROCS = PACK(IGLOBALSLICEPROC, MASK=IGLOBALSLICEPROC.NE.MPI_PROC_NULL)
  CALL MPI_GROUP_INCL(IWRLD_GROUP, ICOUNT, IPROCS, IGROUP_GLOBALSLICE, IERR)
!
!*        2.7     Create MPI communicator associated to new group
!
  CALL MPI_COMM_CREATE(NMNH_COMM_WORLD, IGROUP_GLOBALSLICE, &
                       ICOMM_GLOBALSLICE, IERR)
!
!-------------------------------------------------------------------------------
!
!*        3.      GATHER THE LOCAL SLICES ON ALL PROCS THAT CONTAIN THE SLICE
!                 -----------------------------------------------------------
!
!*        3.1     Have the length of the local slice on each proc known
!                 by all procs on the global slice
!
  IF (ICOMM_GLOBALSLICE .NE. MPI_COMM_NULL) THEN
!
    ALLOCATE(ISIZES(ICOUNT))
    ISIZES = 0
    CALL MPI_ALLGATHER( (/ ISIZE /) , 1, MNHINT_MPI, ISIZES, 1, MNHINT_MPI, &
                     ICOMM_GLOBALSLICE, IERR)
!
!*        3.2     Compute array of displacements in the slice relative to the 
!                 origin of the global domain
!
    ALLOCATE(IDISPL(ICOUNT+1))
    IDISPL(1) = IDISPL1
    DO J=2, ICOUNT+1
      IDISPL(J) = IDISPL(J-1)+ISIZES(J-1)
    ENDDO
    IGLOBALSLICELENGTH = IDISPL(ICOUNT+1) - IDISPL(1)
!
!*        3.3     Have the values of the local slice on each proc known
!                 by all procs on the global slice
!
    CALL MPI_ALLGATHERV(ZPTR, ISIZE, MNHREAL_MPI, PGLOBALSLICE, &
                        ISIZES, IDISPL, MNHREAL_MPI, ICOMM_GLOBALSLICE, IERR)
!
!*        3.4     Delete slice communicator
!
    CALL MPI_COMM_FREE(ICOMM_GLOBALSLICE, IERR)
!
    DEALLOCATE(ISIZES, IDISPL)
!
  ENDIF
!
!*        3.5     Delete slice group 
!
    CALL MPI_GROUP_FREE(IGROUP_GLOBALSLICE, IERR)
!
!-------------------------------------------------------------------------------
!
!*        4.      BROADCAST THE SLICE ON ALL PROCS THAT ARE NOT ON THE SLICE
!                 ----------------------------------------------------------
!
!*        4.1     Create communicator with the first proc
!*                on the slice and the procs that are not on  the
!*                slice
!
  CALL MPI_GROUP_EXCL(IWRLD_GROUP, ICOUNT-1, IPROCS(2:2), IGROUP, IERR)
  CALL MPI_COMM_CREATE(NMNH_COMM_WORLD, IGROUP, ICOMM, IERR)
!
!*        4.2     Broadcast the slice
!
  IF (ICOMM .NE. MPI_COMM_NULL) THEN
!
    CALL MPI_BCAST(IGLOBALSLICELENGTH, 1, MNHINT_MPI, IPROCS(1), ICOMM, IERR)
    CALL MPI_BCAST(PGLOBALSLICE(IDISPL1+1), IGLOBALSLICELENGTH, MNHREAL_MPI, &
                    IPROCS(1), ICOMM, IERR)
!
    CALL MPI_COMM_FREE(ICOMM, IERR)
    CALL MPI_GROUP_FREE(IGROUP, IERR)
!
  ENDIF
!
  CALL MPI_GROUP_FREE(IWRLD_GROUP, IERR)
!
  IF (PRESENT(KERR)) KERR=IERR
!
  DEALLOCATE(IPROCS)
!
!-------------------------------------------------------------------------------
!
       END SUBROUTINE GET_1DGLOBALSLICE_ll
!
!      ###################################################################
       SUBROUTINE GET_2DGLOBALSLICE_ll(PARRAY, HDIR, KLOC, PGLOBALSLICE, &
                                       KB, KE, KKB, KKE, KERR)
!      ###################################################################
!
!!     Purpose
!!     -------
!         The Purpose of this routine is to extract a slice of 
!      3D field PARRAY along the x or y direction
!
!!**   Method
!!     ------
! 
!         An MPI communicator with the processes corresponding to the
!      subdomains intersecting with the slice is built. This
!      communicator is then used to gather the whole slice (i.e. 
!      the global slice) on these procs. The global slice is then
!      broadcasted on all procs that are not on the slice.
! 
!!     External
!!     --------
!      Module MODE_TOOLS_ll
!        LWEST_ll, LSOUTH_ll, LNORTH_ll, LEAST_ll
!        GET_GLOBALDIMS_ll
!
!!    Implicit Arguments
!!    ------------------
!      Module MODD_STRUCTURE_ll
!        type MODELSPLITTING_ll
!
!      Module MODD_VAR_ll
!        NPROC - Number of processors
!        TCRRT_PROCONF -  Current configuration for current model
!        IP - Number of the local processor
!
!      Module MODD_PARAMETERS_ll
!        JPHEXT, JPVEXT - halo size
!
!!    Reference
!!    ---------
!
!     User Interface for the MesoNH Parallel Package
! 
!!    Author
!!    ------
!     P. Kloos (CERFACS)
! 
!!    Modifications
!!    -------------
!     Original 14 August 1998
! 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_PARAMETERS_ll, ONLY: JPHEXT, JPVEXT
  USE MODD_STRUCTURE_ll,  ONLY: MODELSPLITTING_ll
  USE MODD_VAR_ll,        ONLY: NPROC, TCRRT_PROCONF, IP
!
  IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
  REAL, DIMENSION(:,:,:), TARGET, INTENT(IN) :: PARRAY ! horizontal field
  CHARACTER(LEN=1), INTENT(IN)     :: HDIR   ! direction ("X" or "Y")
  INTEGER, INTENT(IN)              :: KLOC   ! coordinate of the slice
                                             ! to extract
                                             ! (in global coordinates)
  REAL, DIMENSION(:,:), INTENT(OUT)  :: PGLOBALSLICE ! output slice
  INTEGER, OPTIONAL    :: KB, KE             ! begin and end positions of the
                                             ! extracted slices in the HDIR
                                             ! direction
  INTEGER, OPTIONAL    :: KKB, KKE           ! begin and end positions of the
                                             ! extracted slices in the vertical
                                             ! direction
  INTEGER, OPTIONAL    :: KERR               ! error code
!
!*       0.2   declarations of local variables
!
  TYPE(MODELSPLITTING_ll), POINTER :: TZSPLIT
  INTEGER :: IDIM, IDIM2
  INTEGER :: INUMPROC, IERR
  INTEGER :: ILOC   ! local location of the slice
  INTEGER :: ICOUNT ! number of relevant procs (= those which are on the slice)
  INTEGER :: ISIZE  ! length of the local slice
  INTEGER, DIMENSION(:), ALLOCATABLE :: IPROCS ! array of procs that are on
                                               ! the slice
  INTEGER, DIMENSION(:), ALLOCATABLE :: ISIZES ! length of the local slice
                                               ! on all procs
  INTEGER, DIMENSION(:), ALLOCATABLE :: IDISPL ! array of locations of the
                                               ! procs in the slice
                                               ! (for MPI_ALLGATHERV)
  INTEGER :: IWRLD_GROUP                ! world group
  INTEGER :: IGROUP_GLOBALSLICE         ! group for the proc on the slice
  INTEGER :: ICOMM_GLOBALSLICE          ! communicator for the proc on the slice
  INTEGER :: IGROUP, ICOMM
  INTEGER, DIMENSION(NPROC) :: IGLOBALSLICEPROC
  INTEGER, DIMENSION(2) :: IOR, IEND, IORP, IENDP, IORE, IENDE  ! splitting
  INTEGER :: IWEST, IEAST, INORTH, ISOUTH
  INTEGER :: IB, IE ! beginning and end of the local slice
  INTEGER :: IDISPL1  ! beginning of the slice (global)
  INTEGER :: J          ! loop index
  INTEGER :: IGLOBALSLICELENGTH ! length of the global slice
  INTEGER :: IGLOBALSLICEHEIGHT
  INTEGER :: JK
  INTEGER, DIMENSION(2) :: IMAX ! maximum dimensions
  REAL, DIMENSION(:,:), ALLOCATABLE :: ZPTR
!
!-------------------------------------------------------------------------------
!
!*       1.    INITIALISATIONS
!              ---------------
 IWEST=0; IEAST=0; INORTH=0; ISOUTH=0
 IB=0; IE=0
 IDISPL1=0  
!
!*       1.1   Get current splitting
!
  IF (LWEST_ll())  IWEST=-1
  IF (LEAST_ll())  IEAST=1
  IF (LNORTH_ll()) INORTH=1
  IF (LSOUTH_ll()) ISOUTH=-1
  TZSPLIT => TCRRT_PROCONF%TSPLITS_B(IP)
  IOR(1)  = TZSPLIT%NXORP+IWEST
  IOR(2)  = TZSPLIT%NYORP+ISOUTH
  IEND(1) = TZSPLIT%NXENDP+IEAST
  IEND(2) = TZSPLIT%NYENDP+INORTH
!
  IORP(1)  = TZSPLIT%NXORP
  IORP(2)  = TZSPLIT%NYORP
  IENDP(1) = TZSPLIT%NXENDP
  IENDP(2) = TZSPLIT%NYENDP
!
  IORE(1)  = TZSPLIT%NXORE
  IORE(2)  = TZSPLIT%NYORE
  IENDE(1) = TZSPLIT%NXENDE
  IENDE(2) = TZSPLIT%NYENDE
!
  IGLOBALSLICEPROC = 0
  CALL GET_GLOBALDIMS_ll(IMAX(1),IMAX(2))
  IMAX(:) = IMAX(:) + 2*JPHEXT
!
!*       1.2   Set dimension (1 for X and 2 for Y)
!
  IDIM  = IACHAR(HDIR)-IACHAR('X')+1
  IDIM2 = 3-IDIM ! IDIM's inverse
!
  IGLOBALSLICEHEIGHT = KKE - KKB + 1
!
!*       1.4   Set beginning and end of local slice
!
  IF (.NOT.PRESENT(KKB) .AND. .NOT.PRESENT(KKE)) THEN
    KKB = 1 + JPVEXT
    KKE = SIZE(PARRAY,3) - JPVEXT
  ENDIF
!
  IF (PRESENT(KB) .AND. PRESENT(KE)) THEN
!
!*      Test the ranges
!
    IF((KB < 1 ) .OR. (KE > IMAX(IDIM))) THEN
!
!   Error
!
      KERR = -1
      RETURN
!
    ENDIF
!
    IB = MAX(KB, IOR(IDIM))  - IORE(IDIM) + 1
    IE = MIN(KE, IEND(IDIM)) - IORE(IDIM) + 1
!    IDISPL1 = KB-1  ! old version
    IDISPL1 = 0
  ELSE ! default : physical domain
    IB = 1+JPHEXT
    IE = IENDP(IDIM)-IORE(IDIM)+1
!    IDISPL1 = JPHEXT ! old version
    IDISPL1 = 0
  ENDIF
!
!-------------------------------------------------------------------------------
!
!*       2.    CREATE MPI COMMUNICATOR WITH THE PROCS ON THE SLICE
!              ---------------------------------------------------
!
!*       2.1   Test if i am on the slice
!              if so, INUMPROC = my MPI rank
!              if not, INUMPROC = MPI_PROC_NULL
!
  IF (KLOC >= IOR(IDIM2) .AND. KLOC <= IEND(IDIM2) .AND. IB<=IE) THEN
!
! Set local location
!
    ILOC = KLOC-IORE(IDIM2)+1
!
! Set relevant procs
!
    INUMPROC = IP-1
!
! Set lenght of the local slice
!
    ISIZE = IE - IB + 1
!
!*        2.2   Have ZPTR point to the slice
!
    SELECT CASE(HDIR)
    CASE("X")
      ALLOCATE(ZPTR(IE-IB+1, IGLOBALSLICEHEIGHT))
      ZPTR = PARRAY(IB:IE,ILOC,KKB:KKE)
!
    CASE("Y")
      ALLOCATE(ZPTR(IE-IB+1, IGLOBALSLICEHEIGHT))
      ZPTR = PARRAY(ILOC,IB:IE,KKB:KKE)
!
    CASE DEFAULT
      call Print_msg( NVERB_FATAL, 'GEN', 'GET_2DGLOBALSLICE_ll', 'invalid HDIR dummy argument ('//hdir//')' )
!
    END SELECT
!
  ELSE
!
    INUMPROC = MPI_PROC_NULL
!
  ENDIF
!
!*        2.3    Gather values of INUMPROC
!
  CALL MPI_ALLGATHER( (/ INUMPROC /) , 1, MNHINT_MPI, IGLOBALSLICEPROC, 1, &
                     MNHINT_MPI, NMNH_COMM_WORLD, IERR)
!
!*        2.4     Get MPI world group
!
  CALL MPI_COMM_GROUP(NMNH_COMM_WORLD, IWRLD_GROUP, IERR)
!
!*        2.5     Count number of proc that contain the slice
!
  ICOUNT = COUNT(IGLOBALSLICEPROC.NE.MPI_PROC_NULL)
!
!*        2.6     Create MPI group with the procs that contain the slice
!
  ALLOCATE(IPROCS(ICOUNT+1))
  IPROCS(1:ICOUNT) = PACK(IGLOBALSLICEPROC, MASK=IGLOBALSLICEPROC.NE.MPI_PROC_NULL)
  IPROCS(ICOUNT+1) = MPI_PROC_NULL
  PRINT *,'SIZE(IPROCS) = ',SIZE(IPROCS)
  CALL MPI_GROUP_INCL(IWRLD_GROUP, ICOUNT, IPROCS, IGROUP_GLOBALSLICE, IERR)
!
!*        2.7     Create MPI communicator associated to new group
!
  CALL MPI_COMM_CREATE(NMNH_COMM_WORLD, IGROUP_GLOBALSLICE, &
                       ICOMM_GLOBALSLICE, IERR)
! 
!-------------------------------------------------------------------------------
!
!*        3.      GATHER THE LOCAL SLICES ON ALL PROCS THAT CONTAIN THE SLICE
!                 -----------------------------------------------------------
!
!*        3.1     Have the length of the local slice on each proc known
!                 by all procs on the global slice
!
  IF (ICOMM_GLOBALSLICE .NE. MPI_COMM_NULL) THEN
!
    ALLOCATE(ISIZES(ICOUNT))
    ISIZES = 0
    CALL MPI_ALLGATHER( (/ ISIZE /) , 1, MNHINT_MPI, ISIZES, 1, MNHINT_MPI, &
                     ICOMM_GLOBALSLICE, IERR)
!
!*        3.2     Compute array of displacements in the slice relative to the
!                 origin of the global domain
!
    ALLOCATE(IDISPL(ICOUNT+1))
    IDISPL(1) = IDISPL1
    DO J=2, ICOUNT+1
      IDISPL(J) = IDISPL(J-1)+ISIZES(J-1)
    ENDDO
    IGLOBALSLICELENGTH = IDISPL(ICOUNT+1) - IDISPL(1)
!
!
!*        3.3     Have the values of the local slice on each proc known
!                 by all procs on the global slice
!
    DO JK = 1, IGLOBALSLICEHEIGHT
      CALL MPI_ALLGATHERV(ZPTR(1,JK), ISIZE, MNHREAL_MPI, &
                          PGLOBALSLICE(1,JK), ISIZES, IDISPL, &
                          MNHREAL_MPI, ICOMM_GLOBALSLICE, IERR)
    ENDDO
!
!*        3.4     Delete slice communicator
!
    CALL MPI_COMM_FREE(ICOMM_GLOBALSLICE, IERR)
!
    DEALLOCATE(ISIZES, IDISPL)
    DEALLOCATE(ZPTR)
!
  ENDIF
!
!*        3.5     Delete slice group
!
    CALL MPI_GROUP_FREE(IGROUP_GLOBALSLICE, IERR)
!
!-------------------------------------------------------------------------------
!
!*        4.      BROADCAST THE SLICE ON ALL PROCS THAT ARE NOT ON THE SLICE
!                 ----------------------------------------------------------
!
!*        4.1     Create communicator with the first proc
!*                on the slice and the procs that are not on  the
!*                slice
!
  CALL MPI_GROUP_EXCL(IWRLD_GROUP, ICOUNT-1, IPROCS(2:2), IGROUP, IERR)
  CALL MPI_COMM_CREATE(NMNH_COMM_WORLD, IGROUP, ICOMM, IERR)
!
!*        4.2     Broadcast the slice
!
  IF (ICOMM .NE. MPI_COMM_NULL) THEN
!
    CALL MPI_BCAST(IGLOBALSLICELENGTH, 1, MNHINT_MPI, IPROCS(1), ICOMM, IERR)
    DO JK = 1, IGLOBALSLICEHEIGHT
      CALL MPI_BCAST(PGLOBALSLICE(1,JK), IGLOBALSLICELENGTH, MNHREAL_MPI, &
                    IPROCS(1), ICOMM, IERR)
    ENDDO
!
    CALL MPI_COMM_FREE(ICOMM, IERR)
!
  ENDIF
!
  CALL MPI_GROUP_FREE(IGROUP, IERR)
!
!  CALL MPI_GROUP_FREE(IWRLD_GROUP, IERR)
!
  IF (PRESENT(KERR)) KERR=IERR
!
  DEALLOCATE(IPROCS)
!
       END SUBROUTINE GET_2DGLOBALSLICE_ll
!
!      ###################################################################
!       SUBROUTINE GET_SLICE_ll(PARRAY, HDIR, KLOC, PSLICE, KB, KE, KERR)
!      ###################################################################
!
!!     PURPOSE
!!     -------
!         The purpose of this routine is to get a slice of the
!      domain in the x or y direction.
!
!      ###################################################################
       SUBROUTINE GET_1DSLICE_ll(PARRAY, HDIR, KLOC, PSLICE, KB, KE, KERR)
!      ###################################################################
!
!!     Purpose
!!     -------  
!         The Purpose of this routine is to extract a slice of an horizontal
!      field PARRAY along the x or y direction
!
!!**   Method
!!     ------
!!
!         An MPI communicator with the processes corresponding to the
!      subdomains intersecting with the slice is built. This
!      communicator is then used to gather the whole slice (i.e. 
!      the global slice) on these procs. The global slice is then
!      broadcasted on all procs that are not on the slice.
! 
!!     External
!!     --------
!      Module MODE_TOOLS_ll
!        LWEST_ll, LSOUTH_ll, LNORTH_ll, LEAST_ll
!        GET_GLOBALDIMS_ll
!
!!    Implicit Arguments
!!    ------------------
!      Module MODD_STRUCTURE_ll
!        type MODELSPLITTING_ll
!
!      Module MODD_VAR_ll
!        NPROC - Number of processors
!        TCRRT_PROCONF -  Current configuration for current model
!        IP - Number of the local processor
!
!      Module MODD_PARAMETERS_ll
!        JPHEXT, JPVEXT - halo size
!
!!    Reference
!!    ---------
!     User Interface for the MesoNH Parallel Package
! 
!!    Author
!!    ------
!     P. Kloos (CERFACS)
! 
!!    Modifications
!!    -------------
!     Original 14 August 1998
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_PARAMETERS_ll, ONLY: JPHEXT
  USE MODD_STRUCTURE_ll,  ONLY: MODELSPLITTING_ll
  USE MODD_VAR_ll,        ONLY: NPROC, TCRRT_PROCONF, IP
!
  IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
  REAL, DIMENSION(:,:), TARGET, INTENT(IN) :: PARRAY ! horizontal field
  CHARACTER(LEN=1), INTENT(IN)     :: HDIR   ! direction ("X" or "Y")
  INTEGER, INTENT(IN)              :: KLOC   ! coordinate of the slice to
                                             ! extract (in global coordinates)
  REAL, DIMENSION(:), INTENT(OUT)  :: PSLICE ! output slice
  INTEGER, OPTIONAL    :: KB, KE             ! begin and end positions of the
                                             ! extracted slices
                                             ! (local coordinates)
  INTEGER, OPTIONAL    :: KERR               ! error code

!
!*       0.2   declarations of local variables
!
  TYPE(MODELSPLITTING_ll), POINTER :: TZSPLIT
  INTEGER :: IDIM, IDIM2
  INTEGER :: INUMPROC, IERR
  INTEGER :: ILOC   ! local location of the slice
  INTEGER :: ICOUNT ! number of relevant procs (= those which are on the slice)
  INTEGER :: ISIZE  ! length of the local slice
  INTEGER, DIMENSION(:), ALLOCATABLE :: IPROCS ! array of procs that are on 
                                               ! the slice 
  INTEGER, DIMENSION(:), ALLOCATABLE :: ISIZES ! length of the local slice
                                               ! on all procs 
  INTEGER, DIMENSION(:), ALLOCATABLE :: IDISPL ! array of locations of the 
                                               ! procs in the slice
                                               ! (for MPI_ALLGATHERV)
  
  REAL, DIMENSION(:), ALLOCATABLE :: ITOTALSLICE

  INTEGER :: IWRLD_GROUP                ! world group
  INTEGER :: IGROUP_SLICE               ! group for the proc on the slice
  INTEGER :: ICOMM_SLICE                ! communicator for the proc on the slice
  INTEGER :: IGROUP, ICOMM
  INTEGER, DIMENSION(NPROC) :: ISLICEPROC
  INTEGER, DIMENSION(2) :: IOR, IEND, IORP, IENDP, IORE, IENDE  ! splitting
  INTEGER :: IWEST, IEAST, INORTH, ISOUTH
  INTEGER :: IB, IE ! beginning and end of the local slice
  INTEGER :: IDISPL1  ! beginning of the slice (global)
  INTEGER :: J          ! loop index
  INTEGER :: ISLICELENGTH ! length of the global slice
  INTEGER, DIMENSION(2) :: IMAX ! maximum dimensions
  REAL, DIMENSION(:), POINTER :: ZPTR
  INTEGER :: IRES
!
  INTEGER :: IIB, IIE, IJB, IJE
!
!-------------------------------------------------------------------------------
!
!*       1.    INITIALISATIONS
!              ---------------
 IWEST=0; IEAST=0; INORTH=0; ISOUTH=0
 IB=0; IE=0
 IDISPL1=0  
!
!*       1.1   Get current splitting
!
  IF (LWEST_ll())  IWEST=-JPHEXT      ! -1
  IF (LEAST_ll())  IEAST=JPHEXT       ! 1
  IF (LNORTH_ll()) INORTH=JPHEXT      ! 1
  IF (LSOUTH_ll()) ISOUTH=-JPHEXT     ! -1
  TZSPLIT => TCRRT_PROCONF%TSPLITS_B(IP)
  IOR(1)  = TZSPLIT%NXORP+IWEST
  IOR(2)  = TZSPLIT%NYORP+ISOUTH 
  IEND(1) = TZSPLIT%NXENDP+IEAST
  IEND(2) = TZSPLIT%NYENDP+INORTH
!
  CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
!
  IIB = IIB+IWEST
  IJB = IJB+ISOUTH
  IIE = IIE+IEAST
  IJE = IJE+INORTH
!
  IORP(1)  = TZSPLIT%NXORP
  IORP(2)  = TZSPLIT%NYORP
  IENDP(1) = TZSPLIT%NXENDP
  IENDP(2) = TZSPLIT%NYENDP
!
  IORE(1)  = TZSPLIT%NXORE
  IORE(2)  = TZSPLIT%NYORE
  IENDE(1) = TZSPLIT%NXENDE
  IENDE(2) = TZSPLIT%NYENDE
!
  ISLICEPROC = 0
  CALL GET_GLOBALDIMS_ll(IMAX(1),IMAX(2))
  IMAX(:) = IMAX(:) + 2*JPHEXT
!
!*       1.2   Set dimension (1 for X and 2 for Y)
!
  IDIM  = IACHAR(HDIR)-IACHAR('X')+1
  IDIM2 = 3-IDIM ! IDIM's inverse 

  ALLOCATE(ITOTALSLICE(IMAX(IDIM)))

!
!*       1.4   Set beginning and end of local slice
!
  IF (PRESENT(KB) .AND. PRESENT(KE)) THEN
!
!*      Test the ranges
!
    IF((KB < 1 ) .OR. (KE > IENDE(IDIM))) THEN
!
!   Error
!
      KERR = -1
      RETURN
!
    ENDIF
!
    IB = KB
    IE = KE
!
!*      Correction in case the user has specified the 
!*      extended subdomain (so that the halo overlapping
!*      between two processors won't be gathered twice).
!
    IF (IB == 1) IB = IOR(IDIM)-IORE(IDIM)+1
    IF (IE == IENDE(IDIM)-IORE(IDIM)+1) IE = IEND(IDIM)-IORE(IDIM)+1
!
!    IDISPL1 = KB-1  ! old version
    IDISPL1 = 0
  ELSE ! default : physical domain
    IB = 1+JPHEXT
    IE = IENDP(IDIM)-IORE(IDIM)+1
!    IDISPL1 = JPHEXT ! old version
    IDISPL1 = 0
  ENDIF
!
!-------------------------------------------------------------------------------
!
!*       2.    CREATE MPI COMMUNICATOR WITH THE PROCS ON THE SLICE
!              ---------------------------------------------------
!
!*       2.1   Test if i am on the slice
!              if so, INUMPROC = my MPI rank
!              if not, INUMPROC = MPI_PROC_NULL
!
!
  IF (KLOC >= IOR(IDIM2) .AND. KLOC <= IEND(IDIM2) .AND. IB<=IE) THEN
!
! Set local location
!
    ILOC = KLOC-IORE(IDIM2)+1
!
! Set relevant procs
!
    INUMPROC = IP-1
!
! Set lenght of the local slice
!
!    ISIZE = IENDE(IDIM) - IORE(IDIM) + 1
!
!*        2.2   Have ZPTR point to the slice
!
    SELECT CASE(HDIR)
    CASE("X")
      ISIZE = IIE - IIB + 1
      ZPTR => PARRAY(IIB:IIE,ILOC)
!
    CASE("Y")
      ISIZE = IJE -IJB + 1
      ZPTR => PARRAY(ILOC,IJB:IJE)
!
    CASE DEFAULT
      call Print_msg( NVERB_FATAL, 'GEN', 'GET_1DSLICE_ll', 'invalid HDIR dummy argument ('//hdir//')' )
!
    END SELECT 
!
  ELSE
!
    INUMPROC = MPI_PROC_NULL
!
  ENDIF
!
!*        2.3    Gather values of INUMPROC 
!
  CALL MPI_ALLGATHER( (/ INUMPROC /) , 1, MNHINT_MPI, ISLICEPROC, 1, MNHINT_MPI, &
                     NMNH_COMM_WORLD, IERR)
!
!*        2.4     Get MPI world group
!
  CALL MPI_COMM_GROUP(NMNH_COMM_WORLD, IWRLD_GROUP, IERR)
!
!*        2.5     Count number of proc that contain the slice
!
  ICOUNT = COUNT(ISLICEPROC.NE.MPI_PROC_NULL)
!
!*        2.6     Create MPI group with the procs that contain the slice
!
  ALLOCATE(IPROCS(ICOUNT+1))
  IPROCS(1:ICOUNT) = PACK(ISLICEPROC, MASK=ISLICEPROC.NE.MPI_PROC_NULL)
  IPROCS(ICOUNT+1) = MPI_PROC_NULL
  CALL MPI_GROUP_INCL(IWRLD_GROUP, ICOUNT, IPROCS, IGROUP_SLICE, IERR)
!
!*        2.7     Create MPI communicator associated to new group
!
  CALL MPI_COMM_CREATE(NMNH_COMM_WORLD, IGROUP_SLICE, ICOMM_SLICE, IERR)
!
!-------------------------------------------------------------------------------
!
!*        3.      GATHER THE LOCAL SLICES ON ALL PROCS THAT CONTAIN THE SLICE
!                 -----------------------------------------------------------
!
!*        3.1     Have the length of the local slice on each proc known
!                 by all procs on the global slice
!
  IF (ICOMM_SLICE .NE. MPI_COMM_NULL) THEN
!
    ALLOCATE(ISIZES(ICOUNT))
    ISIZES = 0
    CALL MPI_ALLGATHER( (/ ISIZE /) , 1, MNHINT_MPI, ISIZES, 1, MNHINT_MPI, &
                     ICOMM_SLICE, IERR)
!
!*        3.2     Compute array of displacements in the slice relative to the 
!                 origin of the global domain
!
    ALLOCATE(IDISPL(ICOUNT+1))
    IDISPL(1) = IDISPL1
    DO J=2, ICOUNT+1
      IDISPL(J) = IDISPL(J-1)+ISIZES(J-1)
    ENDDO
    ISLICELENGTH = IDISPL(ICOUNT+1) - IDISPL(1)
!
!*        3.3     Have the values of the local slice on each proc known
!                 by all procs on the global slice
!
    CALL MPI_ALLGATHERV(ZPTR, ISIZE, MNHREAL_MPI, ITOTALSLICE, ISIZES, &
                        IDISPL, MNHREAL_MPI, ICOMM_SLICE, IERR)
!
    DEALLOCATE(ISIZES, IDISPL)
!
!*        3.4     Delete slice communicator
!
    CALL MPI_COMM_FREE(ICOMM_SLICE, IERR)
!
  ENDIF
!
!*        3.5     Delete slice group 
!
    CALL MPI_GROUP_FREE(IGROUP_SLICE, IERR)
!
!-------------------------------------------------------------------------------
!
!*        4.      BROADCAST THE SLICE ON ALL PROCS THAT ARE NOT ON THE SLICE
!                 ----------------------------------------------------------
!
!*        4.1     Create communicator with the first proc
!*                on the slice and the procs that are not on  the
!*                slice
!
  CALL MPI_GROUP_EXCL(IWRLD_GROUP, ICOUNT-1, IPROCS(2:2), IGROUP, IERR)
  
  CALL MPI_COMM_CREATE(NMNH_COMM_WORLD, IGROUP, ICOMM, IERR)
! CALL MPI_COMM_COMPARE(ICOMM, MPI_COMM_NULL, IRES, IERR)
!
!*        4.2     Broadcast the slice
!
  IF (ICOMM .NE. MPI_COMM_NULL) THEN
!
    CALL MPI_BCAST(ISLICELENGTH, 1, MNHINT_MPI, IPROCS(1), ICOMM, IERR)
    CALL MPI_BCAST(ITOTALSLICE, ISLICELENGTH, MNHREAL_MPI, &
                    IPROCS(1), ICOMM, IERR)
    CALL MPI_COMM_FREE(ICOMM, IERR)
  ENDIF
!
  PSLICE(1:(KE-KB+1)) = &
                       ITOTALSLICE((IORE(IDIM) + KB - 1): (IORE(IDIM) + KE - 1))

  CALL MPI_GROUP_FREE(IGROUP, IERR)
  CALL MPI_GROUP_FREE(IWRLD_GROUP, IERR)
!
  IF (PRESENT(KERR)) KERR=IERR
!
  DEALLOCATE(IPROCS)
  DEALLOCATE(ITOTALSLICE)
!
!-------------------------------------------------------------------------------
!
       END SUBROUTINE GET_1DSLICE_ll
!
!      #######################################################
       SUBROUTINE GET_2DSLICE_ll(PARRAY, HDIR, KLOC, PSLICE, &
                                 KB, KE, KKB, KKE, KERR)
!      #######################################################
!
!!     Purpose
!!     -------
!         The Purpose of this routine is to extract a slice of 
!      3D field PARRAY along the x or y direction
!
!!**   Method
!!     ------
!!
!         An MPI communicator with the processes corresponding to the
!      subdomains intersecting with the slice is built. This
!      communicator is then used to gather the whole slice (i.e. 
!      the global slice) on these procs. The global slice is then
!      broadcasted on all procs that are not on the slice.
! 
!!     External
!!     --------
!      Module MODE_TOOLS_ll
!        LWEST_ll, LSOUTH_ll, LNORTH_ll, LEAST_ll
!        GET_GLOBALDIMS_ll
!
!!    Implicit Arguments
!!    ------------------
!      Module MODD_STRUCTURE_ll
!        type MODELSPLITTING_ll
!
!      Module MODD_VAR_ll
!        NPROC - Number of processors
!        TCRRT_PROCONF -  Current configuration for current model
!        IP - Number of the local processor
!
!      Module MODD_PARAMETERS_ll
!        JPHEXT, JPVEXT - halo size
!
!!    Reference
!!    ---------
!     User Interface for the MesoNH Parallel Package
! 
!!    Author
!!    ------
!     P. Kloos (CERFACS)
! 
!!    Modifications
!!    -------------
!     Original 14 August 1998
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_PARAMETERS_ll, ONLY: JPHEXT, JPVEXT
  USE MODD_STRUCTURE_ll,  ONLY: MODELSPLITTING_ll
  USE MODD_VAR_ll,        ONLY: NPROC, TCRRT_PROCONF, IP
!
  IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
  REAL, DIMENSION(:,:,:), TARGET, INTENT(IN) :: PARRAY ! horizontal field
  CHARACTER(LEN=1), INTENT(IN)     :: HDIR   ! direction ("X" or "Y")
  INTEGER, INTENT(IN)              :: KLOC   ! coordinate of the slice to
                                             ! extract (in global coordinates)
  REAL, DIMENSION(:,:), INTENT(OUT)  :: PSLICE ! output slice
  INTEGER, OPTIONAL    :: KB, KE             ! begin and end positions of the
                                             ! extracted slices in the HDIR 
                                             ! direction
  INTEGER, OPTIONAL    :: KKB, KKE           ! begin and end positions of the
                                             ! extracted slices in the vertical                                              ! direction
  INTEGER, OPTIONAL    :: KERR               ! error code

!
!*       0.2   declarations of local variables
!
  TYPE(MODELSPLITTING_ll), POINTER :: TZSPLIT
  INTEGER :: IDIM, IDIM2
  INTEGER :: INUMPROC, IERR
  INTEGER :: ILOC   ! local location of the slice
  INTEGER :: ICOUNT ! number of relevant procs (= those which are on the slice)
  INTEGER :: ISIZE  ! length of the local slice
  INTEGER, DIMENSION(:), ALLOCATABLE :: IPROCS ! array of procs that are on
                                               ! the slice
  INTEGER, DIMENSION(:), ALLOCATABLE :: ISIZES ! length of the local slice
                                               ! on all procs
  INTEGER, DIMENSION(:), ALLOCATABLE :: IDISPL ! array of locations of the
                                               ! procs in the slice
                                               ! (for MPI_ALLGATHERV)
!
  REAL, DIMENSION(:,:), ALLOCATABLE :: ITOTALSLICE
!
  INTEGER :: IWRLD_GROUP                ! world group
  INTEGER :: IGROUP_SLICE               ! group for the proc on the slice
  INTEGER :: ICOMM_SLICE                ! communicator for the proc on the slice
  INTEGER :: IGROUP, ICOMM
  INTEGER, DIMENSION(NPROC) :: ISLICEPROC
  INTEGER, DIMENSION(2) :: IOR, IEND, IORP, IENDP, IORE, IENDE  ! splitting
  INTEGER :: IWEST, IEAST, INORTH, ISOUTH
  INTEGER :: IB, IE ! beginning and end of the local slice
  INTEGER :: IDISPL1  ! beginning of the slice (global)
  INTEGER :: J          ! loop index
  INTEGER :: ISLICELENGTH ! length of the global slice
  INTEGER :: ISLICEHEIGHT
  INTEGER :: JK
  INTEGER, DIMENSION(2) :: IMAX ! maximum dimensions
  REAL, DIMENSION(:,:), ALLOCATABLE :: ZPTR
!
  INTEGER :: IIB, IIE, IJB, IJE
!
!-------------------------------------------------------------------------------
!
!*       1.    INITIALISATIONS
!              ---------------
 IWEST=0; IEAST=0; INORTH=0; ISOUTH=0
 IB=0; IE=0
 IDISPL1=0  
!
!*       1.1   Get current splitting
!
  IF (LWEST_ll())  IWEST=-1
  IF (LEAST_ll())  IEAST=1
  IF (LNORTH_ll()) INORTH=1
  IF (LSOUTH_ll()) ISOUTH=-1

  TZSPLIT => TCRRT_PROCONF%TSPLITS_B(IP)
  IOR(1)  = TZSPLIT%NXORP+IWEST
  IOR(2)  = TZSPLIT%NYORP+ISOUTH
  IEND(1) = TZSPLIT%NXENDP+IEAST
  IEND(2) = TZSPLIT%NYENDP+INORTH
!
  IORP(1)  = TZSPLIT%NXORP
  IORP(2)  = TZSPLIT%NYORP
  IENDP(1) = TZSPLIT%NXENDP
  IENDP(2) = TZSPLIT%NYENDP
!
  IORE(1)  = TZSPLIT%NXORE
  IORE(2)  = TZSPLIT%NYORE
  IENDE(1) = TZSPLIT%NXENDE
  IENDE(2) = TZSPLIT%NYENDE
!
  CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
!
  IIB = IIB+IWEST
  IJB = IJB+ISOUTH
  IIE = IIE+IEAST
  IJE = IJE+INORTH
!
  ISLICEPROC = 0
  CALL GET_GLOBALDIMS_ll(IMAX(1),IMAX(2))
  IMAX(:) = IMAX(:) + 2*JPHEXT
!
!*       1.2   Set dimension (1 for X and 2 for Y)
!
  IDIM  = IACHAR(HDIR)-IACHAR('X')+1
  IDIM2 = 3-IDIM ! IDIM's inverse
!
  ISLICEHEIGHT = KKE - KKB + 1
  ALLOCATE(ITOTALSLICE(IMAX(IDIM),ISLICEHEIGHT))
!
!*       1.4   Set beginning and end of local slice
!
  IF (.NOT.PRESENT(KKB) .AND. .NOT.PRESENT(KKE)) THEN
    KKB = 1 + JPVEXT
    KKE = SIZE(PARRAY,3) - JPVEXT
  ENDIF
!
  IF (PRESENT(KB) .AND. PRESENT(KE)) THEN
!
!*      Test the ranges
!
    IF((KB < 1 ) .OR. (KE > IMAX(IDIM))) THEN
!
!   Error
!
      KERR = -1
      RETURN
!
    ENDIF
!
    IB = KB
    IE = KE
!
!*      Correction in case the user has specified the
!*      extended subdomain (so that the halo overlapping
!*      between two processors won't be gathered twice).
!
    IF (IB == 1) IB = IOR(IDIM)-IORE(IDIM)+1
    IF (IE == IENDE(IDIM)-IORE(IDIM)+1) IE = IEND(IDIM)-IORE(IDIM)+1
!
!    IDISPL1 = KB-1  ! old version
    IDISPL1 = 0
  ELSE ! default : physical domain
    IB = 1+JPHEXT
    IE = IENDP(IDIM)-IORE(IDIM)+1
!    IDISPL1 = JPHEXT ! old version
    IDISPL1 = 0
  ENDIF
!
!-------------------------------------------------------------------------------
!
!*       2.    CREATE MPI COMMUNICATOR WITH THE PROCS ON THE SLICE
!              ---------------------------------------------------
!
!*       2.1   Test if i am on the slice
!              if so, INUMPROC = my MPI rank
!              if not, INUMPROC = MPI_PROC_NULL
!
  IF (KLOC >= IOR(IDIM2) .AND. KLOC <= IEND(IDIM2) .AND. IB<=IE) THEN
!
! Set local location
!
    ILOC = KLOC-IORE(IDIM2)+1
!
! Set relevant procs
!
    INUMPROC = IP-1
!
! Set lenght of the local slice
!
!    ISIZE = IE - IB + 1
!
!*        2.2   Have ZPTR point to the slice
!
    SELECT CASE(HDIR)
    CASE("X")
      ISIZE = IIE - IIB + 1
      ALLOCATE(ZPTR(ISIZE, ISLICEHEIGHT))
      ZPTR = PARRAY(IIB:IIE,ILOC,KKB:KKE)
!
    CASE("Y")
      ISIZE = IJE - IJB + 1
      ALLOCATE(ZPTR(ISIZE, ISLICEHEIGHT))
      ZPTR = PARRAY(ILOC,IJB:IJE,KKB:KKE)
!
    CASE DEFAULT
      call Print_msg( NVERB_FATAL, 'GEN', 'GET_2DSLICE_ll', 'invalid HDIR dummy argument ('//hdir//')' )
!
    END SELECT
!
  ELSE
!
    INUMPROC = MPI_PROC_NULL
!
  ENDIF
!
!*        2.3    Gather values of INUMPROC
!
  CALL MPI_ALLGATHER( (/ INUMPROC /) , 1, MNHINT_MPI, ISLICEPROC, 1, MNHINT_MPI, &
                     NMNH_COMM_WORLD, IERR)
!
!*        2.4     Get MPI world group
!
  CALL MPI_COMM_GROUP(NMNH_COMM_WORLD, IWRLD_GROUP, IERR)
!
!*        2.5     Count number of proc that contain the slice
!
  ICOUNT = COUNT(ISLICEPROC.NE.MPI_PROC_NULL)
!
!*        2.6     Create MPI group with the procs that contain the slice
!
  ALLOCATE(IPROCS(ICOUNT+1))
  IPROCS(1:ICOUNT) = PACK(ISLICEPROC, MASK=ISLICEPROC.NE.MPI_PROC_NULL)
  IPROCS(ICOUNT+1) = MPI_PROC_NULL
  CALL MPI_GROUP_INCL(IWRLD_GROUP, ICOUNT, IPROCS, IGROUP_SLICE, IERR)
!
!*        2.7     Create MPI communicator associated to new group
!
  CALL MPI_COMM_CREATE(NMNH_COMM_WORLD, IGROUP_SLICE, ICOMM_SLICE, IERR)
!
!-------------------------------------------------------------------------------
!
!*        3.      GATHER THE LOCAL SLICES ON ALL PROCS THAT CONTAIN THE SLICE
!                 -----------------------------------------------------------
!
!*        3.1     Have the length of the local slice on each proc known
!                 by all procs on the global slice
!
  IF (ICOMM_SLICE .NE. MPI_COMM_NULL) THEN
!
    ALLOCATE(ISIZES(ICOUNT))
    ISIZES = 0
    CALL MPI_ALLGATHER( (/ ISIZE /) , 1, MNHINT_MPI, ISIZES, 1, MNHINT_MPI, &
                     ICOMM_SLICE, IERR)
!
!*        3.2     Compute array of displacements in the slice relative to the
!                 origin of the global domain
!
    ALLOCATE(IDISPL(ICOUNT+1))
    IDISPL(1) = IDISPL1
    DO J=2, ICOUNT+1
      IDISPL(J) = IDISPL(J-1)+ISIZES(J-1)
    ENDDO
    ISLICELENGTH = IDISPL(ICOUNT+1) - IDISPL(1)
!
!
!*        3.3     Have the values of the local slice on each proc known
!                 by all procs on the global slice
!
    DO JK = 1, ISLICEHEIGHT
      CALL MPI_ALLGATHERV(ZPTR(1,JK), ISIZE, MNHREAL_MPI, &
                          ITOTALSLICE(1,JK), &
                          ISIZES, IDISPL, MNHREAL_MPI, ICOMM_SLICE, IERR)
    ENDDO
!
!*        3.4     Delete slice communicator
!
    CALL MPI_COMM_FREE(ICOMM_SLICE, IERR)
!
    DEALLOCATE(ISIZES, IDISPL)
    DEALLOCATE(ZPTR)
!
  ENDIF
!
!*        3.5     Delete slice group
!
    CALL MPI_GROUP_FREE(IGROUP_SLICE, IERR)
!
!-------------------------------------------------------------------------------
!
!*        4.      BROADCAST THE SLICE ON ALL PROCS THAT ARE NOT ON THE SLICE
!                 ----------------------------------------------------------
!
!*        4.1     Create communicator with the first proc
!*                on the slice and the procs that are not on  the
!*                slice
!
  CALL MPI_GROUP_EXCL(IWRLD_GROUP, ICOUNT-1, IPROCS(2:2), IGROUP, IERR)
  CALL MPI_COMM_CREATE(NMNH_COMM_WORLD, IGROUP, ICOMM, IERR)
!
!*        4.2     Broadcast the slice
!
  IF (ICOMM .NE. MPI_COMM_NULL) THEN
!
    CALL MPI_BCAST(ISLICELENGTH, 1, MNHINT_MPI, IPROCS(1), ICOMM, IERR)
    DO JK = 1, ISLICEHEIGHT
      CALL MPI_BCAST(ITOTALSLICE(1,JK), ISLICELENGTH, MNHREAL_MPI, &
                    IPROCS(1), ICOMM, IERR)
    ENDDO
!
    CALL MPI_COMM_FREE(ICOMM, IERR)
!
  ENDIF
    PSLICE(1:(KE-KB+1),KKB:KKE) = &
          ITOTALSLICE((IORE(IDIM) + KB - 1):(IORE(IDIM) + KE - 1),&
                      1:ISLICEHEIGHT)
!
  CALL MPI_GROUP_FREE(IGROUP, IERR)
!
  CALL MPI_GROUP_FREE(IWRLD_GROUP, IERR)
!
  IF (PRESENT(KERR)) KERR=IERR
!
  DEALLOCATE(IPROCS)
  DEALLOCATE(ITOTALSLICE)
!
      END SUBROUTINE GET_2DSLICE_ll
!
!     ####################################################
      SUBROUTINE INTERSECTION( TPSPLIT, K, TPZONE, TPRES )
!     ####################################################
!
!!****  *INTERSECTION* - routine to compute the intersections of a zone TPZONE,
!                        in a domain D with the subsets of D resulting 
!                        from a splitting TPSPLIT of D
!                        (TPSPLIT is a any splitting of D in K parts with
!                         or without overlapping)
!
!                        the result TPRES is a splitting
!!
!!    Purpose
!!    -------
!     To compute the correspondants of a processor for each kind of
!     exchange, the exchange's zones have to be computed ; these zones
!     can be seen as intersections between zones of the domain.
!     
!!**  Method
!!    ------
!     For each element of the splitting TPSPLIT, its intersection with
!     the ZONE_ll TPZONE is computed by a max and a min computations on the
!     coordinates of the TPSPLIT'element and TPZONE
! 
!!    External
!!    --------   
!     MAX, MIN - functions which compute the MIN and the MAX between 2 elements
! 
!!    Implicit Arguments
!!    ------------------ 
!         
!!    Reference
!!    ---------
! 
!!    Author
!!    ------
!     R. Guivarch
!!
!!    Modifications
!!    -------------
!     Original 01/05/98
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_STRUCTURE_ll, ONLY : ZONE_ll
  USE MODD_VAR_ll, ONLY : DIMZ
!
  IMPLICIT NONE
!
!
!*       0.1   declarations of arguments
!
  TYPE(ZONE_ll), DIMENSION(:), INTENT(IN) :: TPSPLIT ! Splitting of the domain
!
  INTEGER, INTENT(IN) :: K ! Number of elements of TPSPLIT
!
  TYPE(ZONE_ll), INTENT(IN) :: TPZONE ! Zone to be split
!
  TYPE(ZONE_ll), DIMENSION(:), INTENT(OUT) :: TPRES ! Splitting of the zone
!
!*       0.2   declarations of local variables
!
  INTEGER :: J ! loop control variable
!
!-------------------------------------------------------------------------------
!
!*       1.    LIST AND COMPUTE INTERSECTION BETWEEN TPSPLIT(J) AND TPZONE :
!              -----------------------------------------------------------
!
  DO J = 1, K
! 
!   Which subdomain is the owner of TPSPLIT(J)
    TPRES(J)%NUMBER = TPSPLIT(J)%NUMBER
!
!   Computation of the origin coordinate
    TPRES(J)%NXOR = MAX( TPZONE%NXOR, TPSPLIT(J)%NXOR )
    TPRES(J)%NYOR = MAX( TPZONE%NYOR, TPSPLIT(J)%NYOR )
!
!   Computation of the last coordinate
    TPRES(J)%NXEND = MIN( TPZONE%NXEND, TPSPLIT(J)%NXEND )
    TPRES(J)%NYEND = MIN( TPZONE%NYEND, TPSPLIT(J)%NYEND )
!
!   for z-direction all the domain is considered
    TPRES(J)%NZOR = 1
    TPRES(J)%NZEND = DIMZ
! 
!   if the intersection is void, the result is nullified
    IF((TPRES(J)%NXOR > TPRES(J)%NXEND) .OR. &
       (TPRES(J)%NYOR > TPRES(J)%NYEND) ) &
      TPRES(J) = ZONE_ll ( 0, 0, 0, 0, 0, 0, 0, 0 )
! 
  ENDDO
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE INTERSECTION
!
!     ####################################
      SUBROUTINE ADD_ZONE( TPHEAD, TPELT )
!     ####################################
!
!!****  *ADD_ZONE* - routine to add a zone at the end of a correspondant
! 
!!    Purpose
!!    -------
!     the Purpose of this routine is to add a element of type ZONE_ll to
!     a variable of type CRSPD_ll which is a list of ZONE_ll
!
!!**  Method
!!    ------
!     if the list is void, we create the list and put the element
!     as the first element else we add the element at the end of the list
! 
!!    External
!!    --------   
! 
!!    Implicit Arguments
!!    ------------------ 
!     Module MODD_STRUCTURE_ll
!       types CRSPD_ll, ZONE_ll
!
!     Module MODD_VAR_ll
!       IP - Number of the local processor
!
!!    Reference
!!    ---------
! 
!!    Author
!     ------
!     R. Guivarch
!
!!    Modifications
!!    -------------
!     Original 01/05/98
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!         
  USE MODD_STRUCTURE_ll, ONLY : CRSPD_ll, ZONE_ll
  USE MODD_VAR_ll, ONLY : IP
!
  IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
  TYPE(CRSPD_ll), POINTER :: TPHEAD ! head of the list
!
  TYPE(ZONE_ll), INTENT(IN) :: TPELT ! element to be added
!
!*       0.2   declarations of local variables
!
  TYPE(CRSPD_ll), POINTER :: TZCURRENT, TZNEW ! intermediate variables
!
!-------------------------------------------------------------------------------
!
!*       1.    ADD THE ELEMENT TPELT :
!              ---------------------
!
  IF(.NOT.ASSOCIATED(TPHEAD)) THEN
! 
! first element of the list
!
    ALLOCATE(TPHEAD)
    TPHEAD%TELT = TPELT
    NULLIFY( TPHEAD%TNEXT )
    TPHEAD%NCARD = 1
    IF (TPELT%NUMBER /= IP) THEN
      TPHEAD%NCARDDIF = 1
    ELSE
      TPHEAD%NCARDDIF = 0
    ENDIF
! 
  ELSE
!
! others elements
!
    TZCURRENT => TPHEAD
!
!   Go to the end of the list
    DO WHILE(ASSOCIATED(TZCURRENT%TNEXT))
      TZCURRENT => TZCURRENT%TNEXT
    ENDDO
!
!   Add the element
    ALLOCATE(TZNEW)
    TZNEW%TELT = TPELT
    NULLIFY(TZNEW%TNEXT)
! 
    TZCURRENT%TNEXT => TZNEW
! 
    TPHEAD%NCARD = TPHEAD%NCARD + 1
    IF (TPELT%NUMBER /= IP) TPHEAD%NCARDDIF = TPHEAD%NCARDDIF + 1
! 
  ENDIF
!
!-------------------------------------------------------------------------------
! 
      END SUBROUTINE ADD_ZONE
!
!     ##################################
      INTEGER FUNCTION SIZE_ZONE(TPZONE)
!     ##################################
!
!*       0.    DECLARATIONS
!
  USE MODD_STRUCTURE_ll, ONLY : ZONE_ll
!
  IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
  TYPE(ZONE_ll) :: TPZONE
!
!-------------------------------------------------------------------------------
!
  SIZE_ZONE = &
             (TPZONE%NXEND - TPZONE%NXOR + 1) * (TPZONE%NYEND - TPZONE%NYOR + 1)
!
!-------------------------------------------------------------------------------
!
      END FUNCTION SIZE_ZONE
!
!     #####################################
      INTEGER FUNCTION GET_MAX_SIZE(TPLIST)
!     #####################################
!
!*       0.    DECLARATIONS
!
  USE MODD_STRUCTURE_ll, ONLY : ZONE_ll
!
  IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
  TYPE(CRSPD_ll), POINTER :: TPLIST
!
!*       0.2   declarations of local variables
!
  TYPE(CRSPD_ll), POINTER :: TZLIST
  INTEGER :: KCURSIZE, KMAXSIZE
!
!-------------------------------------------------------------------------------
!
  KMAXSIZE = 0
  TZLIST => TPLIST
  DO WHILE(ASSOCIATED(TZLIST))
    KCURSIZE = SIZE_ZONE(TZLIST%TELT)
    IF (KMAXSIZE < KCURSIZE) KMAXSIZE = KCURSIZE
    TZLIST => TZLIST%TNEXT
  ENDDO
!
  GET_MAX_SIZE = KMAXSIZE
!
!-------------------------------------------------------------------------------
!
      END FUNCTION GET_MAX_SIZE 
!
!     #################################################
      SUBROUTINE EXTRACT_ZONE( TPSPLITS, TPPZS, TPEZS )
!     #################################################
!
!!****  *EXTRACT_ZONE* - routine to construct two splittings variables
!!                       from a MODELSPLITTING_ll variable
! 
!!    Purpose
!!    -------
!     the Purpose of this routine is to extract two splittings TPPZS,
!     physical zone splitting and TPEZS, extended zone splitting
!     from a MODELSPLITTING_ll TPSPLITS
!
!!**  Method
!!    ------
! 
!!    External
!!    --------
! 
!!    Implicit Arguments
!!    ------------------
!     Module MODD_STRUCTURE_ll
!       types MODELSPLITTING_ll, ZONE_ll
!
!     Module MODD_VAR_ll
!        NPROC - Number of processors
! 
!!    Reference
!!    ---------
! 
!!    Author
!!    ------
!     R. Guivarch
! 
!!    Modifications
!!    -------------
!     Original 01/05/98
! 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_STRUCTURE_ll, ONLY : MODELSPLITTING_ll, ZONE_ll
  USE MODD_VAR_ll, ONLY : NPROC
!
  IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
  TYPE(MODELSPLITTING_ll), DIMENSION(:), POINTER :: TPSPLITS
!
  TYPE(ZONE_ll), DIMENSION(:), INTENT(OUT) :: TPPZS, TPEZS
!
!*       0.2   declarations of local variables
!
    INTEGER :: J ! loop control variable
!
!-------------------------------------------------------------------------------
!
!*       1.    FILL TPPZS AND TPEZS FOR EACH J :
!              -------------------------------
!
  DO J = 1, NPROC
!
    TPPZS(J) = ZONE_ll( 0, 0, 0, 0, 0, 0, 0, 0 )
    TPEZS(J) = ZONE_ll( 0, 0, 0, 0, 0, 0, 0, 0 )
!
    TPPZS(J)%NUMBER = TPSPLITS(J)%NUMBER
    TPPZS(J)%NXOR   = TPSPLITS(J)%NXORP
    TPPZS(J)%NYOR   = TPSPLITS(J)%NYORP
    TPPZS(J)%NXEND  = TPSPLITS(J)%NXENDP
    TPPZS(J)%NYEND  = TPSPLITS(J)%NYENDP
!
    TPEZS(J)%NUMBER = TPSPLITS(J)%NUMBER
    TPEZS(J)%NXOR   = TPSPLITS(J)%NXORE
    TPEZS(J)%NYOR   = TPSPLITS(J)%NYORE
    TPEZS(J)%NXEND  = TPSPLITS(J)%NXENDE
    TPEZS(J)%NYEND  = TPSPLITS(J)%NYENDE
!
  ENDDO
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE EXTRACT_ZONE
!
!     #################################################
      SUBROUTINE EXTRACT_ZONE_EXTENDED( TPSPLITS, TPPZS, TPEZS_EXTENDED, HALOSIZE )
!     #################################################
!
!!****  *EXTRACT_ZONE* - routine to construct two splittings variables
!!                       from a MODELSPLITTING_ll variable
!
!!    Purpose
!!    -------
!     the Purpose of this routine is to extract two splittings TPPZS,
!     physical zone splitting and TPEZS_EXTENDED, extended zone splitting with halo of size HALOSIZE
!     from a MODELSPLITTING_ll TPSPLITS
!
!!**  Method
!!    ------
!
!!    External
!!    --------
!
!!    Implicit Arguments
!!    ------------------
!     Module MODD_STRUCTURE_ll
!       types MODELSPLITTING_ll, ZONE_ll
!
!     Module MODD_VAR_ll
!        NPROC - Number of processors
!
!!    Reference
!!    ---------
!
!!    Author
!!    ------
!     R. Guivarch
!
!!    Modifications
!!    -------------
!     Original 01/05/98
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_STRUCTURE_ll, ONLY : MODELSPLITTING_ll, ZONE_ll
  USE MODD_VAR_ll, ONLY : NPROC
!
  IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
  TYPE(MODELSPLITTING_ll), DIMENSION(:), POINTER :: TPSPLITS
!
  TYPE(ZONE_ll), DIMENSION(:), INTENT(OUT) :: TPPZS, TPEZS_EXTENDED
!
  INTEGER, INTENT(IN) :: HALOSIZE
!
!*       0.2   declarations of local variables
!
    INTEGER :: J ! loop control variable
!
!-------------------------------------------------------------------------------
!
!*       1.    FILL TPPZS AND TPEZS FOR EACH J :
!              -------------------------------
!
  DO J = 1, NPROC
!
    TPPZS(J) = ZONE_ll( 0, 0, 0, 0, 0, 0, 0, 0 )
    TPEZS_EXTENDED(J) = ZONE_ll( 0, 0, 0, 0, 0, 0, 0, 0 )
!
    TPPZS(J)%NUMBER = TPSPLITS(J)%NUMBER
    TPPZS(J)%NXOR   = TPSPLITS(J)%NXORP+1
    TPPZS(J)%NYOR   = TPSPLITS(J)%NYORP+1
    TPPZS(J)%NXEND  = TPSPLITS(J)%NXENDP+1
    TPPZS(J)%NYEND  = TPSPLITS(J)%NYENDP+1
!
    IF (  TPSPLITS(J)%NDIMXP < HALOSIZE .OR. TPSPLITS(J)%NDIMYP < HALOSIZE ) THEN
      WRITE(*,*) "WARNING : HALOSIZE is greater than model dimension"
      WRITE(*,*) "HALOSIZE = ", HALOSIZE
      WRITE(*,*) "model dimensions : ", TPSPLITS(J)%NDIMXP, "x", TPSPLITS(J)%NDIMYP
    ENDIF
!
    TPEZS_EXTENDED(J)%NUMBER = TPSPLITS(J)%NUMBER
    TPEZS_EXTENDED(J)%NXOR   = TPSPLITS(J)%NXORP+1-HALOSIZE
    TPEZS_EXTENDED(J)%NYOR   = TPSPLITS(J)%NYORP+1-HALOSIZE
    TPEZS_EXTENDED(J)%NXEND  = TPSPLITS(J)%NXENDP+1+HALOSIZE
    TPEZS_EXTENDED(J)%NYEND  = TPSPLITS(J)%NYENDP+1+HALOSIZE
!
  ENDDO
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE EXTRACT_ZONE_EXTENDED
!
!     ###########################################
      SUBROUTINE GLOBAL2LOCAL(TPPROCONF, TPCRSPD)
!     ###########################################
!
!!****  *GLOBAL2LOCAL* - routine to switch from global coordinates to local ones
! 
!!    Purpose
!!    -------
!     the Purpose of this routine is to compute the coordinates of the 
!     subdomains included in a variable of type CRSPD_ll in the local
!     referential for 2way splitting of the local subdomain
!
!!**  Method
!!    ------
!     the coordinates of the elements of TPCRSPD are global ;
!     the coordinates of the local subdomain in 2way splitting 
!     are in the variable TPPROCONF%TSPLITS_B(IP) ;
! 
!     we substract the coordinates of each element by the origin
!     of the extended subdomain to obtain local coordinates
! 
!!    External
!!    --------   
!!
!!    Implicit Arguments
!!    ------------------ 
!     Module MODD_STRUCTURE_ll
!       types PROCONF_ll, CRSPD_ll, ZONE_ll
!
!     Module MODD_VAR_ll
!       IP - Number of the local processor
!         
!!    Reference
!!    ---------
!!
!!    Author
!!    ------
!!      Ph. Kloos
!!
!!    Modifications
!!    -------------
!!      Original 01/05/98
!!               03/02/99  change declaration of TPPROCONF
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_STRUCTURE_ll, ONLY : PROCONF_ll, CRSPD_ll, ZONE_ll
  USE MODD_VAR_ll, ONLY : IP
!
  IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
  TYPE(PROCONF_ll), POINTER :: TPPROCONF
  TYPE(CRSPD_ll), POINTER :: TPCRSPD ! CRSPD_ll to be switched
!
!
!*       0.2   declarations of local variables
!
! intermediate variables to describe the list and the elements
  TYPE(ZONE_ll), POINTER :: TZZONE
  TYPE(CRSPD_ll), POINTER :: TZCRSPD
!
! 2way-splitting
  TYPE(MODELSPLITTING_ll), POINTER :: TPSPLIT
!
!-------------------------------------------------------------------------------
!
!*       1.    SWITCH :
!              ------
! 
! we point to the structure which contains informations 
! of the 2way local subdomain
  TPSPLIT => TPPROCONF%TSPLITS_B(IP)
!
! we list the variable TPCRSPD of type CRSPD_ll
  TZCRSPD => TPCRSPD
  DO WHILE (ASSOCIATED(TZCRSPD))
    TZZONE => TZCRSPD%TELT
!   
!   we substract the origin of the local subdomain (extended subdomain)
    TZZONE%NXOR = TZZONE%NXOR - TPSPLIT%NXORE + 1
    TZZONE%NXEND = TZZONE%NXEND - TPSPLIT%NXORE + 1
    TZZONE%NYOR = TZZONE%NYOR - TPSPLIT%NYORE + 1
    TZZONE%NYEND = TZZONE%NYEND - TPSPLIT%NYORE + 1
!
    TZCRSPD => TZCRSPD%TNEXT
  ENDDO
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE GLOBAL2LOCAL
!
!     ################################
      SUBROUTINE G2LX(TPSPLIT,TPCRSPD)
!     ################################
!
!!****  *G2LX* - routine to switch from global coordinates to local ones
! 
!!    Purpose
!!    -------
!     the Purpose of this routine is to compute the coordinates of the
!     subdomains included in a variable of type CRSPD_ll in the local
!     referentiel for x-slices or y-slices splitting of the local subdomain
!
!!**  Method
!!    ------
!     the coordinates of the elements of TPCRSPD are global ;
!     the coordinates of the local subdomain in 2way splitting
!     are in the variable TPSPLIT ;
! 
!     we substract the coordinates of each element by the origin
!     of the extended subdomain to obtain local coordinates
! 
!!    External
!!    --------
! 
!!    Implicit Arguments
!!    ------------------
!     Module MODD_STRUCTURE_ll
!       types MODELSPLITTING_ll, CRSPD_ll
!
!!    Reference
!!    ---------
! 
!!    Author
!!    ------
!     R. Guivarch
!
!!    Modifications
!!    -------------
!     Original 01/05/98
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_STRUCTURE_ll, ONLY : MODELSPLITTING_ll, CRSPD_ll
!
  IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
  TYPE(MODELSPLITTING_ll), INTENT(IN) :: TPSPLIT ! x-slices or y-slices
                                                 ! splitting
!
  TYPE(CRSPD_ll), POINTER :: TPCRSPD ! CRSPD_ll to be switch
!
!*       0.2   declarations of local variables
!
! intermediate variables to describe the list and the elements
  TYPE(ZONE_ll), POINTER :: TZZONE
  TYPE(CRSPD_ll), POINTER :: TZCRSPD
!
!-------------------------------------------------------------------------------
!
!*       1.    SWITCH :
!              ------
! we list the variable TPCRSPD of type CRSPD_ll
  TZCRSPD => TPCRSPD
  DO WHILE (ASSOCIATED(TZCRSPD))
    TZZONE => TZCRSPD%TELT
!
!   we substract the origin of the local subdomain (physical subdomain)
    TZZONE%NXOR = TZZONE%NXOR - TPSPLIT%NXORP + 1
    TZZONE%NXEND = TZZONE%NXEND - TPSPLIT%NXORP + 1
    TZZONE%NYOR = TZZONE%NYOR - TPSPLIT%NYORP + 1
    TZZONE%NYEND = TZZONE%NYEND - TPSPLIT%NYORP + 1
!
    TZCRSPD => TZCRSPD%TNEXT
  ENDDO
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE G2LX
!
!     #################################################
      SUBROUTINE GET_OR_SURFEX_ll( HSPLIT, KOR )
!     #################################################
!
!!****  *GET_LOCAL_PORTION_OF_SURFEX_FIELD2D* - returns the origin index of the extended
!                     2way subdomain or of the x-slices subdomain
!                     or of the y-slices
!                     subdomain of the local processor in a surfex field (global indices)
!
!!    Purpose
!!    -------
!!     returns the origin index of the extended
!!                     2way subdomain or of the x-slices subdomain
!!                     or of the y-slices
!!                     subdomain of the local processor in a surfex field (global indices)
!
!!**  Method
!!    ------
!
!!    External
!!    --------
!
!!    Implicit Arguments
!!    ------------------
!
!!    Reference
!!    ---------
!
!!    Author
!!    ------
!     M.Moge
!
!!    Modifications
!!    -------------
!     Original 16/12/14
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_PARAMETERS, ONLY : JPHEXT
!
  IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
  CHARACTER(len=1), INTENT(IN) :: HSPLIT
  INTEGER, INTENT(OUT) :: KOR
!
!*       0.2   declarations of local variables
!
  INTEGER :: IXOR_ll, IYOR_ll ! beginning of local subdomain in global coordinates
!
!-------------------------------------------------------------------------------
!
  CALL GET_OR_ll( HSPLIT, IXOR_ll, IYOR_ll )
  KOR = (IXOR_ll-JPHEXT)*(IYOR_ll-JPHEXT)
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE GET_OR_SURFEX_ll
!
!
!     #################################################
      SUBROUTINE GET_LOCAL_PORTION_OF_SURFEX_FIELD2D( PSURFEXFIELDGLB, POUTPUTFIELDLCL )
!     #################################################
!
!!****  *GET_LOCAL_PORTION_OF_SURFEX_FIELD2D* - extracts local portion of a global
!!                       surfex field (2D field stored in 1D array)
!
!!    Purpose
!!    -------
!     extract local portion of a global
!!    surfex field (2D field stored in 1D array)
!
!!**  Method
!!    ------
!
!!    External
!!    --------
!
!!    Implicit Arguments
!!    ------------------
!
!!    Reference
!!    ---------
!
!!    Author
!!    ------
!     M.Moge
!
!!    Modifications
!!    -------------
!     Original 08/12/14
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_DIM_n, ONLY : NIMAX_ll, NJMAX_ll
  USE MODD_PARAMETERS, ONLY : JPHEXT
!
  IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
  REAL, DIMENSION(:), INTENT(IN) :: PSURFEXFIELDGLB
!
  REAL, DIMENSION(:), INTENT(OUT) :: POUTPUTFIELDLCL
!
!*       0.2   declarations of local variables
!
  INTEGER :: JI,JJ ! loop control variables
  INTEGER :: IXOR, IYOR, IXEND, IYEND ! beginning and end of local subdomain in local coordinates
  INTEGER :: IXOR_ll, IYOR_ll ! beginning of local subdomain in global coordinates
  INTEGER :: ICOUNT
!
!-------------------------------------------------------------------------------
!
  CALL GET_INDICE_ll( IXOR, IYOR, IXEND, IYEND )
  CALL GET_OR_ll( 'B', IXOR_ll, IYOR_ll )
!
  ICOUNT = 1
  DO JJ=IYOR_ll+IYOR-1-JPHEXT,IYOR_ll+IYEND-1-JPHEXT
    DO JI=IXOR_ll+IXOR-1-JPHEXT,IXOR_ll+IXEND-1-JPHEXT
      POUTPUTFIELDLCL(ICOUNT) = PSURFEXFIELDGLB(JI+(NIMAX_ll)*(JJ-1))
      ICOUNT = ICOUNT+1
    ENDDO
  ENDDO
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE GET_LOCAL_PORTION_OF_SURFEX_FIELD2D
!
!
!     #################################################
      SUBROUTINE SET_LOCAL_PORTION_OF_SURFEX_FIELD2D( PFIELDLCL, PSURFEXFIELDGLB )
!     #################################################
!
!!****  *GET_LOCAL_PORTION_OF_SURFEX_FIELD2D* - sets values of local portion of a global
!!                       surfex field (2D field stored in 1D array)
!
!!    Purpose
!!    -------
!     sets values of local portion of a global
!!    surfex field (2D field stored in 1D array)
!
!!**  Method
!!    ------
!
!!    External
!!    --------
!
!!    Implicit Arguments
!!    ------------------
!
!!    Reference
!!    ---------
!
!!    Author
!!    ------
!     M.Moge
!
!!    Modifications
!!    -------------
!     Original 09/12/14
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_DIM_n, ONLY : NIMAX_ll, NJMAX_ll
  USE MODD_PARAMETERS, ONLY : JPHEXT
!
  IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
  REAL, DIMENSION(:), INTENT(IN) :: PFIELDLCL
!
  REAL, DIMENSION(:), INTENT(OUT) :: PSURFEXFIELDGLB
!
!*       0.2   declarations of local variables
!
  INTEGER :: JI,JJ ! loop control variables
  INTEGER :: IXOR, IYOR, IXEND, IYEND ! beginning and end of local subdomain in local coordinates
  INTEGER :: IXOR_ll, IYOR_ll ! beginning of local subdomain in global coordinates
  INTEGER :: ICOUNT
!
!-------------------------------------------------------------------------------
!
  CALL GET_INDICE_ll( IXOR, IYOR, IXEND, IYEND )
  CALL GET_OR_ll( 'B', IXOR_ll, IYOR_ll )
!
  ICOUNT = 1
  DO JJ=IYOR_ll+IYOR-1-JPHEXT,IYOR_ll+IYEND-1-JPHEXT
    DO JI=IXOR_ll+IXOR-1-JPHEXT,IXOR_ll+IXEND-1-JPHEXT
      PSURFEXFIELDGLB(JI+(NIMAX_ll)*(JJ-1)) = PFIELDLCL(ICOUNT)
      ICOUNT = ICOUNT+1
    ENDDO
  ENDDO
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE SET_LOCAL_PORTION_OF_SURFEX_FIELD2D
!
!
!     #################################################
      SUBROUTINE GET_MEAN_OF_COORD_SQRT_ll(PARRAY,KSIZELOC,KSIZEGLB,PMEANSQRT)
!     #################################################
!
!!****  *GET_L2_NORM_ll* - computes the L2 norm of 1D array PARRAY accross all processes
!
!!    Purpose
!!    -------
!     computes the L2 norm of 1D array PARRAY accross all processes
!
!!**  Method
!!    ------
!
!!    External
!!    --------
!
!!    Implicit Arguments
!!    ------------------
!
!!    Reference
!!    ---------
!
!!    Author
!!    ------
!     M.Moge
!
!!    Modifications
!!    -------------
!     Original 10/12/14
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
  REAL, DIMENSION(:), INTENT(IN)  :: PARRAY
  INTEGER,            INTENT(IN)  :: KSIZELOC
  INTEGER,            INTENT(IN)  :: KSIZEGLB
!
  REAL,               INTENT(OUT) :: PMEANSQRT
!
!*       0.2   declarations of local variables
!
  REAL    :: IMEANSQRTLOC
  INTEGER :: IINFO
!
!-------------------------------------------------------------------------------
!
IMEANSQRTLOC = SUM(SQRT(PARRAY))
CALL MPI_ALLREDUCE(IMEANSQRTLOC, PMEANSQRT, 1, MNHREAL_MPI, MPI_SUM, NMNH_COMM_WORLD,IINFO)
PMEANSQRT = PMEANSQRT / KSIZEGLB
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE GET_MEAN_OF_COORD_SQRT_ll
!
!     ##########################################################################
      FUNCTION SPREAD_X_ll(HSPLIT, PSOURCE, KDIM, KX, KCOPIES) RESULT(PSPREAD_X)
!     ##########################################################################
!
!!****  *SPREAD_X_ll* - perform the spread in y-direction of the local
!!                   part of a x-vector for the local processor
! 
!!    Purpose
!!    -------
!
!!**  Method
!!    ------
!     In function of the position of the local subdomain, we extract
!     the local part of the x-vector and spread it in y-direction.
!     the spread is done with the dimension of the local subdomain in
!     y-direction

!!    External
!!    --------
!     Module MODE_TOOLS_ll
!       GET_DIM_EXT_ll
!       GET_OR_ll
!
!!    Implicit Arguments
!!    ------------------
!
!!    Reference
!!    ---------
!
!!    Author
!!    ------
!     R. Guivarch
! 
!!    Modifications
!!    -------------
!     Original 01/05/98
! 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  IMPLICIT NONE
!
!
!*       0.0   declarations of arguments
!
  CHARACTER(len=1), INTENT(IN) :: HSPLIT ! Splitting flag (B, X or Y)
!
  REAL, DIMENSION(:), INTENT(IN) :: PSOURCE ! x-vector
!
  INTEGER, INTENT(IN) :: KDIM ! direction of the spread
!
  INTEGER, INTENT(IN) :: KX, KCOPIES ! dimension of the local subdomain
!
!*       0.1   declaration of returned variable
!
  REAL, DIMENSION( KX , KCOPIES ) :: PSPREAD_X

!*       0.2   declarations of local variables
!
  REAL, ALLOCATABLE :: ZSUBSET(:) ! Intermediate buffer
!
  INTEGER :: INDEX 
  INTEGER :: JI ! loop control variable
!
  INTEGER :: IXOR, IYOR, IDIMX, IDIMY ! origin of the local subdoamin
!
!-------------------------------------------------------------------------------
!
!*       1.    GET THE ORIGIN AND THE DIMENSION OF THE LOCAL SUBDOMAIN :
!              -------------------------------------------------------
!
  CALL GET_OR_ll( HSPLIT, IXOR, IYOR )
  CALL GET_DIM_EXT_ll( HSPLIT, IDIMX, IDIMY )
!
!-------------------------------------------------------------------------------
!
!*       2.    ALLOCATION OF THE INTERMEDIATE BUFFER :
!              -------------------------------------
!
  ALLOCATE(ZSUBSET(IDIMX))
!
!-------------------------------------------------------------------------------
!
!*       3.    FILL THE INTERMEDIATE BUFFER :
!              ----------------------------
!
  INDEX = 0
  DO JI = IXOR, IXOR + IDIMX - 1
!
    INDEX = INDEX + 1
    ZSUBSET(INDEX) = PSOURCE(JI)
!
  ENDDO
!
!-------------------------------------------------------------------------------
!
!*       4.    SPREAD THE BUFFER IN Y-DIRECTION WITH KCOPIES
!
  PSPREAD_X = SPREAD(ZSUBSET, KDIM, KCOPIES )
!
!-------------------------------------------------------------------------------
!
!*       5.    DEALLOCATION OF THE INTERMEDIATE BUFFER :
!              -------------------------------------
!
  DEALLOCATE(ZSUBSET)
!
!-------------------------------------------------------------------------------
!
      END FUNCTION SPREAD_X_ll
!
!     ##########################################################################
      FUNCTION SPREAD_Y_ll(HSPLIT, PSOURCE, KDIM, KY, KCOPIES) RESULT(PSPREAD_Y)
!     ##########################################################################
!
!!****  *SPREAD_Y_ll* - perform the spread in x-direction of the local
!!                   part of a y-vector for the local processor
! 
!!    Purpose
!!    -------
!
!!**  Method
!!    ------
!     In function of the position of the local subdomain, we extract
!     the local part of the y-vector and spread it in x-direction.
!     the spread is done with the dimension of the local subdomain in
!     x-direction
! 
!!    External
!!    --------
!     Module MODE_TOOLS_ll
!       GET_DIM_EXT_ll, GET_OR_ll
!
!!    Implicit Arguments
!!    ------------------
! 
!!    Reference
!!    ---------
! 
!!    Author
!!    ------
!     R. Guivarch
! 
!!    Modifications
!!    -------------
!     Original 01/05/98
! 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  IMPLICIT NONE
!
!
!*       0.0   declarations of arguments
!
  CHARACTER(len=1), INTENT(IN) :: HSPLIT ! Splitting flag (B, X or Y)
!
  REAL, DIMENSION(:), INTENT(IN) :: PSOURCE ! x-vector
!
  INTEGER, INTENT(IN) :: KDIM ! direction of the spread
!
  INTEGER, INTENT(IN) :: KY, KCOPIES ! dimension of the local subdomain
!
!*       0.1   declaration of returned variable
!
  REAL, DIMENSION( KCOPIES , KY ) :: PSPREAD_Y
!

!*       0.2   declarations of local variables
!
  REAL, ALLOCATABLE :: ZSUBSET(:) ! Intermediate buffer
!
  INTEGER :: INDEX
  INTEGER :: JJ ! loop control variable
!
  INTEGER :: IXOR, IYOR, IDIMX, IDIMY ! origin of the local subdoamin
!
!-------------------------------------------------------------------------------
!
!*       1.    GET THE ORIGIN AND THE DIMENSION OF THE LOCAL SUBDOMAIN :
!              -------------------------------------------------------
!
  CALL GET_OR_ll( HSPLIT, IXOR, IYOR )
  CALL GET_DIM_EXT_ll( HSPLIT, IDIMX, IDIMY )
!
!-------------------------------------------------------------------------------
!
!*       2.    ALLOCATION OF THE INTERMEDIATE BUFFER :
!              -------------------------------------
!
  ALLOCATE(ZSUBSET(IDIMY))
!
!-------------------------------------------------------------------------------
!
!*       3.    FILL THE INTERMEDIATE BUFFER :
!              ----------------------------
!
  INDEX = 0
  DO JJ = IYOR, IYOR + IDIMY - 1
!
    INDEX = INDEX + 1
    ZSUBSET(INDEX) = PSOURCE(JJ)
!
  ENDDO
!
!-------------------------------------------------------------------------------
!
!*       4.    SPREAD THE BUFFER IN X-DIRECTION WITH KCOPIES
!
  PSPREAD_Y = SPREAD(ZSUBSET, KDIM, KCOPIES )
!
!-------------------------------------------------------------------------------
!
!*       5.    DEALLOCATION OF THE INTERMEDIATE BUFFER :
!              -------------------------------------
!
  DEALLOCATE(ZSUBSET)
!
!-------------------------------------------------------------------------------
!
      END FUNCTION SPREAD_Y_ll
!
!     #################################################################
      FUNCTION SPREAD_XY_ll( HSPLIT, PSOURCE, KDIM, KX, KY, KCOPIES ) &
      RESULT( PSPREAD_XY )
!     #################################################################
!
!!****  *SPREAD_XY_ll* - perform the spread in z-direction of the local
!!                    part of a 2D-array for the local processor
! 
!!    Purpose
!!    -------
!
!!**  Method
!!    ------
!     In function of the position of the local subdomain, we extract
!     the local part of the 2D-array and spread it in z-direction.
!     the spread is done with KCOPIES
 
!!    External
!!    --------
!     Module MODE_TOOLS_ll
!       GET_DIM_EXT_ll, GET_OR_ll
!
!!    Implicit Arguments
!!    ------------------
!
!!    Reference
!!    ---------
!
!!    Author
!!    ------
!     R. Guivarch
! 
!!    Modifications
!!    -------------
!     Original 01/05/98
! 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  IMPLICIT NONE
!
!
!*       0.0   declarations of arguments
!
  CHARACTER(len=1), INTENT(IN) :: HSPLIT ! Splitting flag (B, X or Y)
!
  REAL, DIMENSION(:,:), INTENT(IN) :: PSOURCE ! x-vector
!
  INTEGER, INTENT(IN) :: KDIM ! direction of the spread
!
  INTEGER, INTENT(IN) :: KX, KY ! dimension of the local subdomain
!
  INTEGER, INTENT(IN) :: KCOPIES ! number of spread
!
!*       0.1   declaration of returned variable
!
  REAL, DIMENSION( KX , KY, KCOPIES ) :: PSPREAD_XY

!*       0.2   declarations of local variables
!
  REAL, ALLOCATABLE :: ZSUBSET(:,:) ! Intermediate buffer
!
  INTEGER :: INDEXX, INDEXY
  INTEGER :: JI, JJ ! loop control variable
!
  INTEGER :: IXOR, IYOR, IDIMX, IDIMY ! origin of the local subdoamin
!
!-------------------------------------------------------------------------------
!
!*       1.    GET THE ORIGIN AND THE DIMENSION OF THE LOCAL SUBDOMAIN :
!              -------------------------------------------------------
!
  CALL GET_OR_ll( HSPLIT, IXOR, IYOR )
  CALL GET_DIM_EXT_ll( HSPLIT, IDIMX, IDIMY )
!
!-------------------------------------------------------------------------------
!
!*       2.    ALLOCATION OF THE INTERMEDIATE BUFFER :
!              -------------------------------------
!
  ALLOCATE(ZSUBSET(IDIMX,IDIMY))
!
!-------------------------------------------------------------------------------
!
!*       3.    FILL THE INTERMEDIATE BUFFER :
!              ----------------------------
  INDEXY = 0
!
  DO JJ = IYOR, IYOR + IDIMY - 1
!
    INDEXY = INDEXY + 1
    INDEXX = 0
!
    DO JI = IXOR, IXOR + IDIMX - 1
!
      INDEXX = INDEXX + 1
      ZSUBSET(INDEXX, INDEXY) = PSOURCE(JI,JJ)
!
    ENDDO
!
  ENDDO
!
!-------------------------------------------------------------------------------
!
!*       4.    SPREAD THE BUFFER IN Z-DIRECTION WITH KCOPIES
!
  PSPREAD_XY = SPREAD(ZSUBSET, KDIM, KCOPIES )
!
!-------------------------------------------------------------------------------
!
!*       5.    DEALLOCATION OF THE INTERMEDIATE BUFFER :
!              -----------------------------------------
!
  DEALLOCATE(ZSUBSET)
!
!-------------------------------------------------------------------------------
!
      END FUNCTION SPREAD_XY_ll
!
END MODULE MODE_TOOLS_ll
