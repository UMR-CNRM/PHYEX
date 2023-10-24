!MNH_LIC Copyright 2002-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!       #######################
        MODULE  MODD_ELEC_DESCR
!       #######################
!
!!****  *MODD_ELEC_DESCR* - declaration of the electrical descriptive constants
!!
!!	PURPOSE
!!	-------
!
!!**	IMPLICIT ARGUMENTS
!!	------------------
!!	  None
!!
!!	REFERENCE
!!	---------
!!
!!	AUTHOR
!!	------
!!       Gilles Molinie    * Laboratoire d'Aerologie *
!!
!!	MODIFICATIONS
!!	-------------
!!	  Original	14/11/02
!!        M. Chong      26/01/10  Small ions parameters
!!                               +Option for Fair weather field from
!!                               Helsdon-Farley (JGR, 1987, 5661-5675)
!!                               Add "Beard" effect via sedimentation process
!!        J.-P. Pinty   25/10/13 Add "Latham" effect via aggregation process
!!        C. Barthe     05/07/23 New data structures for PHYEX - for sedimentation in ICE3
!!                               + Remove unused variables
!!
!-------------------------------------------------------------------------------
!
!*	0.	DECLARATIONS
!		------------
!
USE MODD_PARAMETERS, ONLY: JPSVNAMELGTMAX

IMPLICIT NONE
!
! Namelist
LOGICAL :: LOCG=.FALSE.          ! .T.: only charge generation
LOGICAL :: LELEC_FIELD=.TRUE.    ! .T.: the electric field is computed
LOGICAL :: LFW_HELFA=.FALSE.  ! .T.: Helsdon-Farley Fair Weather field
LOGICAL :: LCOSMIC_APPROX=.FALSE.  ! .T.: Neglecting height variations of fair
                                   ! weather ion current in calculating ion
                                   ! source (XIONSOURCEFW) from cosmic rays
LOGICAL :: LION_ATTACH = .TRUE.    ! .T.: Ion attachment to hydrometeors
CHARACTER (LEN=3) :: CDRIFT = 'PPM' ! PPM (advection) or DIV (Divergence form) 
LOGICAL :: LRELAX2FW_ION = .FALSE. ! .T.= Relaxation to fair weather ion
                                   ! concentration in rim zone and top absorbing
                                   ! layer
LOGICAL :: LFLASH_GEOM=.TRUE.    ! .T.: the 'geometric' flash scheme is used
LOGICAL :: LSAVE_COORD=.FALSE.   ! .T.: the flash coord are written in an ascii file
INTEGER :: NFLASH_WRITE = 1000   ! Number of flashes to be saved before writing
                                 ! the diag and/or coordinates in ascii files
LOGICAL :: LINDUCTIVE=.FALSE.      ! .T.: inductive process is taken into account
LOGICAL :: LLNOX_EXPLICIT=.FALSE.  ! .T.: lnox production is computed
LOGICAL :: LSERIES_ELEC=.FALSE.    ! .T.: looking for flash rate proxies
INTEGER :: NTSAVE_SERIES = 60      ! time interval at which data from 
                                   ! series_cloud_elec are written in an ascii file
!
CHARACTER (LEN=5) :: CNI_CHARGING='TAKAH'  ! Choice of the charging process
REAL :: XQTC=263.   ! temperature charge reversal for 'HELFA'
REAL :: XLIM_NI_IS=10.E-15   ! max magnitude of dq for I-S non-inductive charging (C)
REAL :: XLIM_NI_IG=30.E-15   ! max magnitude of dq for I-G non-inductive charging (C)
REAL :: XLIM_NI_SG=100.E-15  ! max magnitude of dq for S-G non-inductive charging (C)
!
CHARACTER (LEN=5) :: CLSOL='RICHA'    ! Choice of the Laplace equation solver
INTEGER :: NLAPITR_ELEC=4       ! Nb of iteration for the elec field solveur
REAL    :: XRELAX_ELEC=1        ! Relaxation factor for the elec field solveur
!
REAL :: XETRIG=200.E3       ! E threshold for lightning triggering
REAL :: XEBALANCE=0.1       ! Proportion of XETRIG that the lightning must reduce
REAL :: XEPROP=15.E3        ! E threshold for lightning propagation
!
REAL :: XQEXCES=2.E-10      ! Charge in excess of qexces => pt available for cell detection
REAL :: XQNEUT=1.E-10       ! Charge in excess of qneut is neutralized
REAL :: XDFRAC_ECLAIR=2.3   ! Fractal dimension of lightning flashes
REAL :: XDFRAC_L=1500.      ! Linear coefficient for the branch number
!
REAL :: XWANG_A = 0.34E21 ! Wang eta al. parameters of
REAL :: XWANG_B = 1.3E16  ! LNOX production
!
!
REAL, DIMENSION(:), SAVE, ALLOCATABLE :: XQTMIN ! Min values allowed for the 
                                                ! volumetric charge
REAL, DIMENSION(:) ,      ALLOCATABLE :: XRTMIN_ELEC    ! Limit value of R where charge is available
!
REAL       :: XEPSILON        ! Dielectric permittivity of air (F/m)
REAL       :: XECHARGE        ! Elementary charge (C)
!
!
! parameters relative to electrification
!
REAL, SAVE :: XLBDAR_MAXE, &  ! Max values allowed for the shape
              XLBDAS_MAXE, &  ! when computation of charge separation
              XLBDAG_MAXE, &  ! and of lightning neutralisation
              XLBDAH_MAXE     ! 
REAL :: XALPHACQ, XNUCQ, XLBDACQ
!
!
! parameters relative to the electric field
!
REAL :: XE_0, XKEF            ! Constant for the fair weather electric field

REAL, SAVE :: XE0_HF, XA1_HF, XB1_HF, XA2_HF, XB2_HF, XA3_HF, XB3_HF ! Coeffs.
                                   ! Helsdon-Farley Fair Weather Electric Field
REAL, SAVE :: XIONCOMB             ! Ionic recombination coefficient (m3/s)
REAL, SAVE :: XF_POS, XF_NEG      ! Constant for positive/negative ion mobility
                                   ! law (m2/V/s)
REAL, SAVE :: XEXPMOB              ! Exponent of ion mobility law (m-1)

REAL, SAVE :: XFCORONA             ! Factor for corona current (A m /V3)
REAL, SAVE :: XECORONA             ! Electric field threshold for corona (V/m)

! Fair Weather electric property (Chiu, JGR 1978, 5025-5049)
!
REAL, SAVE :: XJCURR_FW            ! Air-earth conduction current (A/m2)
!
! Lightning flashes
!
INTEGER :: NMAX_CELL  ! max number of electrified cells in the domain
INTEGER :: NBRANCH_MAX  ! max number of branches per flash
INTEGER :: NLEADER_MAX ! max number of segments in the bi-leader
REAL    :: XE_THRESH    ! electric field threshold for cell detection
!
INTEGER, PARAMETER  :: NLGHTMAX = 5000, &            ! Nb max of lightnings 
                       NSEGMAX  = 500                ! Nb max of segments
!
! Parameters relative to the lightning
! 
INTEGER :: NNBLIGHT=0              ! Nb of lightning flashes
!
REAL, DIMENSION(:), ALLOCATABLE :: XNEUT_POS, XNEUT_NEG
INTEGER :: NNB_CG          ! Nb of CG flashes
INTEGER :: NNB_CG_POS      ! Nb of positive CG flashes
REAL    :: XALT_CG         ! Altitude (m) at which CG are detected
!
CHARACTER(LEN=JPSVNAMELGTMAX), DIMENSION(8) &
         :: CELECNAMES=(/'QNIONP','QCELEC','QRELEC','QIELEC','QSELEC',   &
                         'QGELEC','QHELEC','QNIONN'/)
! QNIONP (QNIONN): Positive (Negative) ion concentration
! basenames of the SV articles stored in the binary files
!
REAL :: XLNOX_ECLAIR
!
REAL, DIMENSION(:,:),   ALLOCATABLE :: XEPOTFW_TOP
!
! Parameters relative to the "Beard" effect ELEC=>MICROPHYS
!
LOGICAL :: LSEDIM_BEARD=.FALSE.    ! .T.: to enable ELEC=>MICROPHYS via
!                                  ! particule sedimentation rate
LOGICAL :: LIAGGS_LATHAM=.FALSE.   ! .T.: to enable ELEC=>MICROPHYS via
!                                  ! ice aggregation rate
!
! The following variables must be declared with a derived type to match with PHYEX requirements
TYPE ELEC_DESCR_t
  REAL :: XFC, XFR, XFI, XFS, XFG, XFH ! f_x in q_x = e_x D^f_x
  REAL :: XCXR            ! Exponent in the concentration-slope
END TYPE ELEC_DESCR_t
!
TYPE(ELEC_DESCR_t), SAVE, TARGET :: ELEC_DESCR
!
REAL, POINTER :: XFC => NULL(), &
                 XFR => NULL(), &
                 XFI => NULL(), &
                 XFS => NULL(), &
                 XFG => NULL(), &
                 XFH => NULL(), &
                 XCXR => NULL()
!
CONTAINS
!
SUBROUTINE ELEC_DESCR_ASSOCIATE()
  IMPLICIT NONE
  !
  XFC => ELEC_DESCR%XFC
  XFR => ELEC_DESCR%XFR
  XFI => ELEC_DESCR%XFI
  XFS => ELEC_DESCR%XFS
  XFG => ELEC_DESCR%XFG
  XFH => ELEC_DESCR%XFH
  XCXR => ELEC_DESCR%XCXR
END SUBROUTINE ELEC_DESCR_ASSOCIATE
!
END MODULE MODD_ELEC_DESCR
