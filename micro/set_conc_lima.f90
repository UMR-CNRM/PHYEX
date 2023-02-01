!MNH_LIC Copyright 2000-2020 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!#######################################
module mode_set_conc_lima
!#######################################

implicit none

contains

!     ###########################################################################
      SUBROUTINE SET_CONC_LIMA( kmi, HGETCLOUD, PRHODREF, PRT, PSVT )
!     ###########################################################################
!
!!****  *SET_CONC_LIMA * - initialize droplet, raindrop and ice
!!                   concentration for a RESTArt simulation of the LIMA scheme
!!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to initialize cloud droplet and rain drop
!!    concentrations when the cloud droplet and rain drop mixing ratios are
!!    only available (generally from a previous run using the Kessler scheme).
!!      This routine is used to initialize the droplet/drop concentrations
!!    using the r_c and r_r of a previous REVE or KESS run but also to compute
!!    the LB tendencies in ONE_WAY$n in case of grid-nesting when the optional
!!    argument PTIME is set (a LIMA run embedded in a KESS or REVE run).
!!
!!**  METHOD
!!    ------
!!      The method assumes a Csk law for the activation of aerososl with "s"
!!    the supersaturation (here 0.05 % is chosen). A Marshall-Palmer law with
!!    N_o=10**(-7) m**(-4) is assumed for the rain drop concentration.
!!      The initialization of the PSVT is straightforward for the cloud droplets
!!    while N_r=N_0/Lambda_r with Rho*r_r=Pi*Rho_w*N_0/(Lambda_r**4) is used for
!!    the rain drops. The HGETCLOUD test is used to discriminate between the
!!    'REVE' and 'KESS' options for CCLOUD in the previous run (from which
!!     PRT was calculated).
!!
!!    EXTERNAL
!!    --------
!!      None
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_RAIN_C2R2_DESCR, ONLY : XRTMIN, XCTMIN
!!      Module MODD_RAIN_C2R2_KHKO_PARAM, ONLY : XCONCC_INI, XCONCR_PARAM_INI
!!      Module MODD_CONF,            ONLY : NVERB
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation ( routine SET_CONC_RAIN_C2R2 )
!!
!!    AUTHOR
!!    ------
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!      P. Jabouille     * CNRM/GMME *
!!      B. Vié           * CNRM/GMME *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    15/11/00
!!                        2014 G.Delautier : remplace MODD_RAIN_C2R2_PARAM par MODD_RAIN_C2R2_KHKO_PARAM        *
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!  B. Vié      03/03/2020: secure physical tests
!  P. Wautelet 04/06/2020: correct array start for microphys. concentrations + add kmi dummy argument
!                          (this subroutine is also called for other models)
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAM_LIMA,      ONLY : XRTMIN, XCTMIN, LCOLD, LWARM, LRAIN, NMOD_CCN, NMOD_IFN, &
                                 NMOM_C, NMOM_R, NMOM_I
USE MODD_PARAM_LIMA_COLD, ONLY : XAI, XBI, XAS, XBS
USE MODD_PARAM_LIMA_MIXED,ONLY : XAG, XBG, XAH, XBH
USE MODD_NSV,             ONLY : NSV_LIMA_BEG_A, NSV_LIMA_NC_A, NSV_LIMA_NR_A, NSV_LIMA_CCN_ACTI_A, &
                                 NSV_LIMA_NI_A, NSV_LIMA_NS_A, NSV_LIMA_NG_A, NSV_LIMA_NH_A, NSV_LIMA_IFN_NUCL_A
USE MODD_CST,             ONLY : XPI, XRHOLW, XRHOLI
USE MODD_CONF,            ONLY : NVERB
USE MODD_CONF_n,          ONLY : NRR
USE MODD_LUNIT_n,         ONLY : TLUOUT
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
integer,                   intent(in) :: kmi        ! Model number
CHARACTER (LEN=4),         INTENT(IN) :: HGETCLOUD  ! Get indicator
REAL, DIMENSION(:,:,:),    INTENT(IN) :: PRHODREF   ! Reference density
!
REAL, DIMENSION(:,:,:,:),  INTENT(INOUT) :: PRT     ! microphysical mixing ratios
!
REAL,  DIMENSION(:,:,:,NSV_LIMA_BEG_A(kmi):), INTENT(INOUT):: PSVT     ! microphys. concentrations
!
!
!*       0.2   Declarations of local variables :
!
INTEGER    :: IRESP   ! Return code of FM routines
INTEGER    :: ILUOUT  ! Logical unit number of output-listing
REAL       :: ZCONC
!
!-------------------------------------------------------------------------------
!*       1.    RETRIEVE LOGICAL UNIT NUMBER
!              ----------------------------
!
ILUOUT = TLUOUT%NLU
!
!*       2.    INITIALIZATION
!              --------------
!
IF (LWARM .AND. NRR.GE.2 .AND. NMOM_C.GE.2) THEN
!
!  droplets
!
   ZCONC = 300.E6 ! droplet concentration set at 300 cm-3
   WHERE ( PRT(:,:,:,2) > 1.E-11 )
      PSVT(:,:,:,NSV_LIMA_NC_A(kmi)) = ZCONC
   END WHERE
   WHERE ( PRT(:,:,:,2) <= 1.E-11 )
      PRT(:,:,:,2)  = 0.0
      PSVT(:,:,:,NSV_LIMA_NC_A(kmi)) = 0.0
   END WHERE
   
   IF (NMOD_CCN .GE. 1) THEN
      WHERE ( PRT(:,:,:,2) > 1.E-11 )
         PSVT(:,:,:,NSV_LIMA_CCN_ACTI_A(kmi)) = ZCONC
      END WHERE
      WHERE ( PRT(:,:,:,2) <= 1.E-11 )
         PSVT(:,:,:,NSV_LIMA_CCN_ACTI_A(kmi)) = 0.0
      END WHERE
   END IF
   
!   IF( NVERB >= 5 ) THEN
!      WRITE (UNIT=ILUOUT,FMT=*) "!INI_MODEL$n: The droplet concentration has "
!      WRITE (UNIT=ILUOUT,FMT=*) "been roughly initialised"
!   END IF
END IF
!
IF (LWARM .AND. LRAIN .AND. NRR.GE.3 .AND. NMOM_R.GE.2) THEN
!
!  drops
!
   ZCONC = (1.E7)**3/(XPI*XRHOLW) ! cf XCONCR_PARAM_INI in ini_rain_c2r2.f90
   IF (HGETCLOUD == 'INI1') THEN ! init from REVE scheme
      PSVT(:,:,:,NSV_LIMA_NR_A(kmi)) = 0.0
   ELSE ! init from KESS, ICE3...
      WHERE ( PRT(:,:,:,3) > 1.E-11 )
         PSVT(:,:,:,NSV_LIMA_NR_A(kmi)) = MAX( SQRT(SQRT(PRHODREF(:,:,:)*PRT(:,:,:,3) &
              *ZCONC)),1. )
      END WHERE
      WHERE ( PRT(:,:,:,3) <= 1.E-11 )
         PRT(:,:,:,3)  = 0.0
         PSVT(:,:,:,NSV_LIMA_NR_A(kmi)) = 0.0
      END WHERE
!      IF( NVERB >= 5 ) THEN
!         WRITE (UNIT=ILUOUT,FMT=*) "!INI_MODEL$n: The raindrop concentration has "
!         WRITE (UNIT=ILUOUT,FMT=*) "been roughly initialised"
!      END IF
   END IF
END IF
!
IF (LCOLD .AND. NRR.GE.4 .AND. NMOM_I.GE.2) THEN
!
! ice crystals
!
   ZCONC = 100.E3 ! maximum ice concentration set at 100/L
   WHERE ( PRT(:,:,:,4) > 1.E-11 )
!
!      PSVT(:,:,:,NSV_LIMA_NI_A(kmi)) = MIN( PRHODREF(:,:,:) /                                     &
!           ( XRHOLI * XAI*(10.E-06)**XBI * PRT(:,:,:,4) ), &
!           ZCONC )
! Correction
      PSVT(:,:,:,NSV_LIMA_NI_A(kmi)) = MIN(PRT(:,:,:,4)/(0.82*(10.E-06)**2.5),ZCONC )
   END WHERE
   WHERE ( PRT(:,:,:,4) <= 1.E-11 )
      PRT(:,:,:,4)  = 0.0
      PSVT(:,:,:,NSV_LIMA_NI_A(kmi)) = 0.0
   END WHERE

   IF (NMOD_IFN .GE. 1) THEN
      WHERE ( PRT(:,:,:,4) > 1.E-11 )
         PSVT(:,:,:,NSV_LIMA_IFN_NUCL_A(kmi)) = PSVT(:,:,:,NSV_LIMA_NI_A(kmi))
      END WHERE
      WHERE ( PRT(:,:,:,4) <= 1.E-11 )
         PSVT(:,:,:,NSV_LIMA_IFN_NUCL_A(kmi)) = 0.0
      END WHERE
   END IF

!   IF( NVERB >= 5 ) THEN
!      WRITE (UNIT=ILUOUT,FMT=*) "!INI_MODEL$n: The cloud ice concentration has "
!      WRITE (UNIT=ILUOUT,FMT=*) "been roughly initialised"
!   END IF
!
END IF
!
IF (NSV_LIMA_NS_A(KMI).GE.1) THEN
!
!  snow
!
   ZCONC = 1./ (XAS*0.001**XBS) ! 1mm particle size
   WHERE ( PRT(:,:,:,5) > 1.E-11 )
      PSVT(:,:,:,NSV_LIMA_NS_A(KMI)) = PRT(:,:,:,5) * ZCONC
   ELSEWHERE
      PRT(:,:,:,5)  = 0.0
      PSVT(:,:,:,NSV_LIMA_NS_A(KMI)) = 0.0
   END WHERE
END IF
!
IF (NSV_LIMA_NG_A(KMI).GE.1) THEN
!
!  graupel
!
   ZCONC = 1./ (XAG*0.001**XBG) ! 1mm particle size
   WHERE ( PRT(:,:,:,6) > 1.E-11 )
      PSVT(:,:,:,NSV_LIMA_NG_A(KMI)) = PRT(:,:,:,6) * ZCONC
   ELSEWHERE
      PRT(:,:,:,6)  = 0.0
      PSVT(:,:,:,NSV_LIMA_NG_A(KMI)) = 0.0
   END WHERE
END IF
!
IF (NSV_LIMA_NH_A(KMI).GE.1) THEN
!
!  hail
!
   ZCONC = 1./ (XAH*0.001**XBH) ! 1mm particle size
   WHERE ( PRT(:,:,:,7) > 1.E-11 )
      PSVT(:,:,:,NSV_LIMA_NH_A(KMI)) = PRT(:,:,:,7) * ZCONC
   ELSEWHERE
      PRT(:,:,:,7)  = 0.0
      PSVT(:,:,:,NSV_LIMA_NH_A(KMI)) = 0.0
   END WHERE
END IF
!
END SUBROUTINE SET_CONC_LIMA

end module mode_set_conc_lima
