!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$ $Date$
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!     ################
      MODULE MODD_DYN
!     ################
!
!!****  *MODD_DYN* - declaration of dynamic control variables 
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to declare the dynamic
!     control variables for all models.    
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_PARAMETERS : contains the maximum number of coupling files
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (module MODD_DYN)
!!      Technical Specifications Report of the Meso-NH (chapters 2 and 3)
!!          
!!    AUTHOR
!!    ------
!!	V. Ducrocq   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original      06/05/94
!!      Modifications 17/10/94   (Stein)  For LCORIO              
!!      Modifications 12/12/94   (Stein)  remove LABSLAYER + add XALZBOT
!!                                        and XALKBOT     
!!      Modifications 10/03/95  (I.Mallet) add coupling variables
!!                    04/05/07  (C.Lac) Separation of num.diffusion
!!                                      between variables
!---------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_PARAMETERS
!
IMPLICIT NONE
!
REAL   ,SAVE   :: XSEGLEN    ! Duration of segment (in seconds)
REAL   ,SAVE   :: XASSELIN   ! Asselin coefficient
REAL   ,SAVE   :: XASSELIN_SV! Asselin coefficient for tracer variables
LOGICAL,SAVE   :: LCORIO     ! Coriolis flag
LOGICAL,SAVE   :: LNUMDIFU   ! logical switch for the NUMerical DIFFusion:
                             ! .TRUE. active .FALSE. unactive
                             ! for momentum
LOGICAL,SAVE   :: LNUMDIFTH  ! For theta and mixing ratio
LOGICAL,SAVE   :: LNUMDIFSV  ! For scalar variables           
LOGICAL,SAVE   :: LSTEADYLS  ! logical switch to remove all Larger Scale fields
                             ! evolution during the segment (the LS fields are 
                             ! STEADY).TRUE. LS steady .FALSE. LS unsteady
REAL   ,SAVE   :: XALKTOP    ! Damping coef. at the top of the absorbing
                             ! layer
REAL   ,SAVE   :: XALZBOT    ! Height of the absorbing layer base
!
REAL   ,SAVE   :: XALKGRD    ! Damping coef. at the top of the absorbing
                             ! layer (bottom layer)
REAL   ,SAVE   :: XALZBAS    ! Height of the absorbing layer base  (bottom layer)
!
INTEGER,SAVE  :: NCPL_NBR    ! NumBeR of CouPLing files
INTEGER,SAVE  :: NCPL_CUR    ! Number of CURrent CouPLing file
!
LOGICAL, SAVE, DIMENSION(JPMODELMAX) :: LUSERV_G, LUSERC_G, LUSERR_G, LUSERI_G
LOGICAL, SAVE, DIMENSION(JPMODELMAX) :: LUSERS_G, LUSERH_G, LUSERG_G
LOGICAL, SAVE, DIMENSION(JPMODELMAX) :: LUSETKE
LOGICAL, SAVE, DIMENSION(JPSVMAX,JPMODELMAX) :: LUSESV
REAL,    SAVE, DIMENSION(JPCPLFILEMAX,JPMODELMAX)          :: NCPL_TIMES   ! array of
                ! the number for the coupling instants of every nested model
REAL,    SAVE                       :: XTSTEP_MODEL1  ! time step of the
                                                      ! outermost model
END MODULE MODD_DYN
