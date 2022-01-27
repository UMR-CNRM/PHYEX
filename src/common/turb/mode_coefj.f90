!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 turb 2006/05/18 13:07:25
!-----------------------------------------------------------------
!################
MODULE MODE_COEFJ
IMPLICIT NONE
CONTAINS
!    #######################################################
     FUNCTION COEFJ(PTHL,PEXNREF,PFRAC_ICE)   RESULT(PCOEFJ)
!    #######################################################
!
!      PURPOSE
!!     -------
!      COEFJ computes the coefficient J of the documentation.
!
!!**   METHOD
!!     ------                                 rvs(Tl) Lv(Tl)
!!       The value of this coefficient is J = --------------, for rc only
!!                                             Rv Tl THETAl
!!
!!           rvsw(Tl) Lv(Tl)                rvsi(Tl) Ls(Tl)
!!       or  --------------- (1-Pfrac_ri) + --------------- Pfrac_ri, for rc+ri.
!!            Rv Tl THETAl                   Rv Tl THETAl
!!
!!     EXTERNAL
!!     --------
!!       None.
!!
!!     IMPLICIT ARGUMENTS
!!     ------------------
!!       Module MODD_CST : contains physical constants.
!!   
!!     REFERENCE
!!     ---------
!!       Book 1 of documentation of Meso-NH
!!       Book 2 of documentation of Meso-NH
!!
!!
!!     AUTHOR
!!     ------
!!       Jean-Marie Carriere      * Meteo-France *
!!
!!     MODIFICATIONS
!!     -------------
!!       Original       20/03/95
!!       J.-P. Pinty    20/02/03 add non-precipitating ice
!!
!! ----------------------------------------------------------------------
!
!*       0. DECLARATIONS
!           ------------
USE MODD_CST
!
IMPLICIT NONE
!
!*       0.1 declarations of arguments and result
!
REAL, DIMENSION(:,:,:),  INTENT(IN)  ::   PTHL     ! Temperature variable
REAL, DIMENSION(:,:,:),  INTENT(IN)  ::   PEXNREF  ! Exner function of the 
!                                                          reference state
REAL, DIMENSION(:,:,:),  INTENT(IN), OPTIONAL ::   PFRAC_ICE 
                                                   ! Fraction of ri in the 
                                                   !   non-precipating
                                                   !  "rc+ri" condensate
!
REAL,DIMENSION(SIZE(PTHL,1),SIZE(PTHL,2),SIZE(PTHL,3)):: PCOEFJ ! result
!
!*       0.2 declarations of local variables
!
REAL,DIMENSION(SIZE(PTHL,1),SIZE(PTHL,2),SIZE(PTHL,3)) ::       &
                                          ZTL, ZL, ZES, ZRVS, ZP
!                ZTL = Tl, ZL = Lv(Tl) or Ls(Tl), ZES = esw(Tl) or esi(Tl)
!                ZRVS = rvsw(Tl) or rvsi(Tl), ZP = p
!
REAL                                 ::   ZEPS     ! = Mv/Md
!---------------------------------------------------------------------------
!
!*       1. COMPUTATION OF Tl
!           -----------------
!
ZTL(:,:,:) = PTHL(:,:,:) * PEXNREF(:,:,:)
!
!*       2. COMPUTATION OF Lv(Tl)
!           ---------------------
!
ZL(:,:,:) = XLVTT + ( XCPV - XCL ) * ( ZTL(:,:,:) -XTT )
!
!*       3. COMPUTATION OF rvs(Tl)
!           ----------------------
!
ZEPS      = XMV/XMD
ZP(:,:,:) = (PEXNREF(:,:,:)**(XCPD/XRD))*XP00
ZES(:,:,:)  = EXP( XALPW - XBETAW/ZTL(:,:,:) - XGAMW*ALOG(ZTL(:,:,:) ) )
ZRVS(:,:,:) =  ZES(:,:,:) * ZEPS / ( ZP(:,:,:) - ZES(:,:,:) )             
!
!        4. RESULT FOR rc only
!           ------------------
!
PCOEFJ(:,:,:) =    ZRVS(:,:,:)*ZL(:,:,:)/   &
             (  XRV*ZTL(:,:,:)*PTHL(:,:,:)  )
!
! Add case when rc+ri
!
IF(PRESENT(PFRAC_ICE)) THEN
!
!*       5. COMPUTATION OF Ls(Tl)
!           ---------------------
!
  ZL(:,:,:) = XLSTT + ( XCPV - XCI ) * ( ZTL(:,:,:) -XTT )
!
!*       6. COMPUTATION OF rvs(Tl)
!           ----------------------
!
  ZES(:,:,:)  = EXP( XALPI - XBETAI/ZTL(:,:,:) - XGAMI*ALOG(ZTL(:,:,:) ) )
  ZRVS(:,:,:) =  ZES(:,:,:) * ZEPS / ( ZP(:,:,:) - ZES(:,:,:) )
!
!        7. RESULT FOR rc and ri
!           --------------------
!
  PCOEFJ(:,:,:) = (1.0 - PFRAC_ICE(:,:,:))*PCOEFJ(:,:,:)            &
                       + PFRAC_ICE(:,:,:) *ZRVS(:,:,:)*ZL(:,:,:)/   &
                                     (  XRV*ZTL(:,:,:)*PTHL(:,:,:)  )
END IF
!
!---------------------------------------------------------------------------
!
END FUNCTION COEFJ
END MODULE MODE_COEFJ
