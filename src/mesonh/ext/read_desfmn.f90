!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ######################
      MODULE MODI_READ_DESFM_n
!     ######################
!
INTERFACE
!
      SUBROUTINE READ_DESFM_n(KMI,TPDATAFILE,HCONF,OFLAT,OUSERV,                 &
                   OUSERC,OUSERR,OUSERI,OUSECI,OUSERS,OUSERG,OUSERH,             &
                   OUSECHEM,OUSECHAQ,OUSECHIC,OCH_PH,OCH_CONV_LINOX,OSALT,       &
                   ODEPOS_SLT,ODUST,ODEPOS_DST, OCHTRANS,                        &
                   OORILAM,ODEPOS_AER,OLG,OPASPOL,                               &
#ifdef MNH_FOREFIRE
                   OFOREFIRE,                                                    &
#endif
                   OLNOX_EXPLICIT,                                               &
                   OCONDSAMP,OBLOWSNOW,                                          &
                   KRIMX,KRIMY,KSV_USER,                                         &
                   HTURB,HTOM,ORMC01,HRAD,HDCONV,HSCONV,HCLOUD,HELEC,HEQNSYS     )
!
USE MODD_IO, ONLY: TFILEDATA
USE MODD_PARAMETERS
!
INTEGER,            INTENT(IN)  :: KMI    ! Model index
TYPE(TFILEDATA),    INTENT(IN)  :: TPDATAFILE ! Datafile
CHARACTER (LEN=5),  INTENT(OUT) :: HCONF  ! configuration var. linked to FMfile
LOGICAL,            INTENT(OUT) :: OFLAT  ! Logical for zero orography 
LOGICAL,            INTENT(OUT) :: OUSERV ! use Rv mixing ratio               
LOGICAL,            INTENT(OUT) :: OUSERC ! use Rc mixing ratio               
LOGICAL,            INTENT(OUT) :: OUSERR ! use Rr mixing ratio               
LOGICAL,            INTENT(OUT) :: OUSERI ! use Ri mixing ratio               
LOGICAL,            INTENT(OUT) :: OUSECI ! use Ci concentration of Ice cristals
LOGICAL,            INTENT(OUT) :: OUSERS ! use Rs mixing ratio               
LOGICAL,            INTENT(OUT) :: OUSERG ! use Rg mixing ratio               
LOGICAL,            INTENT(OUT) :: OUSERH ! use Rh mixing ratio
LOGICAL,            INTENT(OUT) :: OUSECHEM ! Chemical flag
LOGICAL,            INTENT(OUT) :: OUSECHAQ ! Aqueous Chemical flag
LOGICAL,            INTENT(OUT) :: OUSECHIC ! Ice phase Chemical flag
LOGICAL,            INTENT(OUT) :: OCH_PH   ! pH flag
LOGICAL,            INTENT(OUT) :: OCH_CONV_LINOX ! LiNOX flag
LOGICAL,            INTENT(OUT) :: OLG      ! lagrangian flag
LOGICAL,            INTENT(OUT) :: OSALT    ! Sea Salt flag
LOGICAL,            INTENT(OUT) :: ODUST    ! Dust flag
LOGICAL,            INTENT(OUT) :: OPASPOL  ! Passive pollutant flag
#ifdef MNH_FOREFIRE
LOGICAL,            INTENT(OUT) :: OFOREFIRE! ForeFire flag
#endif
LOGICAL,            INTENT(OUT) :: OLNOX_EXPLICIT ! explicit LNOx flag
LOGICAL,            INTENT(OUT) :: OCONDSAMP! Conditional sampling flag
LOGICAL,            INTENT(OUT) :: OBLOWSNOW ! Blowing snow flag
LOGICAL,            INTENT(OUT) :: OORILAM  ! Orilam flag
LOGICAL,            INTENT(OUT) :: OCHTRANS ! Deep convection on scalar
LOGICAL,DIMENSION(JPMODELMAX),INTENT(OUT) :: ODEPOS_DST    ! Dust Wet Deposition flag
LOGICAL,DIMENSION(JPMODELMAX),INTENT(OUT) :: ODEPOS_SLT    ! Sea Salt Wet Deposition flag
LOGICAL,DIMENSION(JPMODELMAX),INTENT(OUT) :: ODEPOS_AER    ! Aerosols Wet Deposition flag
INTEGER,            INTENT(OUT) :: KRIMX, KRIMY        ! number of points for the
                       ! horizontal relaxation for the outermost verticals
INTEGER,            INTENT(OUT) :: KSV_USER    ! number of additional scalar
                                          ! variables in FMfile 
CHARACTER (LEN=4),  INTENT(OUT) :: HTURB  ! Kind of turbulence parameterization
                                          ! used to produce the FMfile
CHARACTER (LEN=4),  INTENT(OUT) :: HTOM   ! Kind of third order moment
LOGICAL,            INTENT(OUT) :: ORMC01 ! flag for RMC01 SBL computations
CHARACTER (LEN=4),  INTENT(OUT) :: HRAD   ! Kind of radiation scheme
CHARACTER (LEN=4),  INTENT(OUT) :: HDCONV ! Kind of deep convection scheme
CHARACTER (LEN=4),  INTENT(OUT) :: HSCONV ! Kind of shallow convection scheme
CHARACTER (LEN=4),  INTENT(OUT) :: HCLOUD ! Kind of microphysical scheme
CHARACTER (LEN=4),  INTENT(OUT) :: HELEC  ! Kind of electrical scheme
CHARACTER (LEN=*),  INTENT(OUT) :: HEQNSYS! type of equations' system
END SUBROUTINE READ_DESFM_n
!
END INTERFACE
!
END MODULE MODI_READ_DESFM_n
!     #########################################################################
      SUBROUTINE READ_DESFM_n(KMI,TPDATAFILE,HCONF,OFLAT,OUSERV,                 &
                   OUSERC,OUSERR,OUSERI,OUSECI,OUSERS,OUSERG,OUSERH,             &
                   OUSECHEM,OUSECHAQ,OUSECHIC,OCH_PH,OCH_CONV_LINOX,OSALT,       &
                   ODEPOS_SLT,ODUST,ODEPOS_DST, OCHTRANS,                        &
                   OORILAM,ODEPOS_AER,OLG,OPASPOL,                               &
#ifdef MNH_FOREFIRE
                   OFOREFIRE,                                                    &
#endif
                   OLNOX_EXPLICIT,                                               &
                   OCONDSAMP,OBLOWSNOW,                                          &
                   KRIMX,KRIMY,KSV_USER,                                         &
                   HTURB,HTOM,ORMC01,HRAD,HDCONV,HSCONV,HCLOUD,HELEC,HEQNSYS     )
!     #########################################################################
!
!!****  *READ_DESFM_n * - routine to read  the descriptor file DESFM
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to read the descriptor file called 
!     DESFM. 
!       
!!
!!**  METHOD
!!    ------
!!      The descriptor file is read. Namelists (NAMXXXn) which contain
!!    informations linked to one nested model are at the beginning of the file.
!!    Namelists (NAMXXX) which contain variables common to all models
!!    are at the end of the file. When the  model index is different from 1, 
!!    the end of the file (namelists NAMXXX) is not read. 
!!      Some attributes of the FMfile are saved in order to check coherence
!!    between initial file and the segment to perform (description given by
!!    EXSEG file), i.e. :
!!      - the configuration which has been used to produce the initial file
!!        (CCONF)
!!      - logical switch for flat configuration (zero orography) in initial file
!!        (LFLAT)
!!      - kind of moist variables in initial file (LUSERV,LUSERC,LUSERR,
!!        LUSERI,LUSERS,LUSERG,LUSERH)
!!      - number of additional scalar variables in initial file (NSV_USER)
!!      - kind of turbulence parameterization  used to produce the initial
!!        file (CTURB)
!!      - kind of mixing length  used to produce the initial
!!        file (CTURBLEN)
!!      - time step of each model stored in PTSTEP_OLD, to correct the initial 
!!        field at t-dt in routine READ_FIELD in case of time step change
!!      - type of equation system in order to verify that the anelastic is the
!!        same for the initila file generation and the run  
!!       
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!     Module MODN_CONF : CCONF,LFLAT,CEQNSYS
!!
!!     Module MODN_CONF1 : LUSERV,LUSERC,LUSERR,LUSERI,LUSECI,
!!                 LUSERS,LUSERG,LUSERH
!!
!!     Module MODN_PARAM1 : CTURB,CRAD,CDCONV,CSCONV
!!                   
!!     Module MODN_TURB$n : CTURBLEN
!!
!!     Module MODN_DYN$n : NRIMX,NRIMY
!!
!!    REFERENCE
!!    ---------
!!      Book2 of the documentation (routine READ_DESFM_n)
!!      
!!
!!    AUTHOR
!!    ------
!!  	V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       07/06/94
!!      Modifications  17/10/94  (Stein)  For LCORIO
!!      Modifications  26/10/94  (Stein)  remove NAM_GET from the Namelists
!!      present in DESFM + change the namelist names
!!      Modifications  09/01/95  (Stein)  add the turbulence scheme
!!      Modifications  09/01/95  (Stein)  add the 1D switch
!!      Modifications  13/02/95  (Stein)  save HTURBLEN
!!      Modifications  30/06/95  (Stein)  add new namelists
!!      Modifications  18/08/95  (Lafore) time step change
!!      Modifications  15/09/95  (Pinty)  add the radiations
!!      Modifications  06/02/96  (J.Vila) add the new scalar advection scheme
!!      Modifications  20/02/96  (Stein)  add the LES namelist + cleaning
!!      Modifications  25/04/96  (Suhre)  add NAM_BLANK
!!      Modifications  25/04/96  (Suhre)  add NAM_FRC
!!      Modifications  25/04/96  (Suhre)  add NAM_CH_MNHCn and NAM_CH_SOLVER
!!      Modifications  11/04/96  (Pinty)  add the ice concentration
!!      Modifications  11/01/97  (Pinty)  add the deep convection
!!      Modifications  22/07/96  (Lafore) gridnesting implementation
!!      Modifications  22/06/97  (Stein ) save the equations' system+ cleaning
!!      Modifications  09/07/97  (Masson) add NAM_PARAM_GROUND
!!      Modifications  25/08/97  (Masson) add HGROUND
!!      Modifications  25/10/97  (Stein ) new namelists
!!      Modification   04/06/00  (Pinty)  add C2R2 scheme 
!!      Modification   22/01/01  (Gazen)  Add  OUSECHEM and OLG
!!      Modification   15/10/01  (Mallet) allow namelists in different orders
!!      Modification   29/11/02  (Pinty)  add C3R5, ICE2, ICE4, ELEC
!!      Modification   01/2004   (Masson) removes surface (externalization)
!!      Modification   01/2005   (Masson) removes 1D and 2D switches
!!      Modification   03/2005   (Tulet)  add dust, aerosols
!!      Modification   03/2006   (O.Geoffroy) Add KHKO scheme
!!      Modification   04/2010   (M. Leriche) Add aqueous + ice chemistry
!!      Modification   07/2013   (Bosseur & Filippi) Adds Forefire
!!      Modification   01/2015   (C. Barthe) Add explicit LNOx
!!      Modification   2016      (B.VIE) LIMA
!!      Modification   11/2016   (Ph. Wautelet) Allocate/initialise some output/backup structures
!!      Modification   02/2018   (Q.Libois) ECRAD
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!      Modification   07/2017   (V. Vionnet) Add blowing snow scheme
!!      Modification   02/2021   (F.Auguste)  add IBM
!!                               (T.Nagel)    add turbulence recycling
!!                               (E.Jezequel) add stations read from CSV file
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
USE MODD_IO,      ONLY: TFILEDATA
USE MODD_LUNIT_n, ONLY: TLUOUT
USE MODD_PARAMETERS
!
USE MODN_BACKUP
USE MODN_BUDGET
USE MODN_CONF
USE MODN_DYN
USE MODN_NESTING
USE MODN_OUTPUT
USE MODN_LES
USE MODN_CONF_n
USE MODN_DYN_n
USE MODN_ADV_n
USE MODN_PARAM_n
USE MODN_PARAM_RAD_n
USE MODN_PARAM_ECRAD_n
USE MODN_PARAM_KAFR_n
USE MODN_PARAM_MFSHALL_n
USE MODN_PARAM_ICE, ONLY : NAM_PARAM_ICE, ZWARM=>LWARM, ZSEDIC=>LSEDIC, &
                           ZPRISTINE_ICE=>CPRISTINE_ICE, ZSEDIM=>CSEDIM
USE MODN_LUNIT_n
USE MODN_LBC_n
USE MODN_NUDGING_n
USE MODN_TURB_n
USE MODN_FRC
USE MODN_BLANK_n
USE MODN_CH_SOLVER_n
USE MODN_CH_MNHC_n
USE MODN_PARAM_C2R2, ONLY : HPARAM_CCN_C2R2=>HPARAM_CCN,HINI_CCN_C2R2=>HINI_CCN, &
                            HTYPE_CCN_C2R2=>HTYPE_CCN,LRAIN_C2R2=>LRAIN, &
                            LSEDC_C2R2=>LSEDC,LACTIT_C2R2=>LACTIT,XCHEN_C2R2=>XCHEN, &
                            XKHEN_C2R2=>XKHEN,XMUHEN_C2R2=>XMUHEN, &
                            XBETAHEN_C2R2=>XBETAHEN,XCONC_CCN_C2R2=>XCONC_CCN, &
                            XR_MEAN_CCN_C2R2=>XR_MEAN_CCN,XLOGSIG_CCN_C2R2=>XLOGSIG_CCN, &
                            XFSOLUB_CCN_C2R2=>XFSOLUB_CCN,XACTEMP_CCN_C2R2=>XACTEMP_CCN, &
                            XALPHAC_C2R2=>XALPHAC,XNUC_C2R2=>XNUC,XALPHAR_C2R2=>XALPHAR, &
                            XNUR_C2R2=>XNUR,XAERDIFF_C2R2=>XAERDIFF, &
                            XAERHEIGHT_C2R2=>XAERHEIGHT,NAM_PARAM_C2R2
USE MODN_PARAM_C1R3, ONLY : XALPHAI_C1R3=>XALPHAI,XNUI_C1R3=>XNUI,XALPHAS_C1R3=>XALPHAS, &
                            XNUS_C1R3=>XNUS,XALPHAG_C1R3=>XALPHAG,XNUG_C1R3=>XNUG, &
                            XFACTNUC_DEP_C1R3=>XFACTNUC_DEP, &
                            XFACTNUC_CON_C1R3=>XFACTNUC_CON,LSEDI_C1R3=>LSEDI, &
                            LHHONI_C1R3=>LHHONI,CPRISTINE_ICE_C1R3,CHEVRIMED_ICE_C1R3, &
                            NAM_PARAM_C1R3
USE MODN_ELEC
USE MODN_SERIES
USE MODN_SERIES_n
USE MODN_TURB_CLOUD
USE MODN_TURB
USE MODN_CH_ORILAM
USE MODN_DUST
USE MODN_SALT
USE MODN_PASPOL
USE MODN_VISCOSITY
USE MODN_DRAG_n
#ifdef MNH_FOREFIRE
USE MODN_FOREFIRE
#endif
USE MODN_CONDSAMP
USE MODN_LATZ_EDFLX
USE MODN_2D_FRC
USE MODN_BLOWSNOW_n
USE MODN_BLOWSNOW
USE MODN_STATION_n
!
USE MODN_PARAM_LIMA
!
USE MODE_MSG
USE MODE_POS
USE MODN_RECYCL_PARAM_n
USE MODN_IBM_PARAM_n
USE MODD_IBM_LSF, ONLY: LIBM_LSF
!
IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
!
!
INTEGER,            INTENT(IN)  :: KMI    ! Model index
TYPE(TFILEDATA),    INTENT(IN)  :: TPDATAFILE ! Datafile
CHARACTER (LEN=5),  INTENT(OUT) :: HCONF  ! configuration var. linked to FMfile
LOGICAL,            INTENT(OUT) :: OFLAT  ! Logical for zero orography 
LOGICAL,            INTENT(OUT) :: OUSERV ! use Rv mixing ratio               
LOGICAL,            INTENT(OUT) :: OUSERC ! use Rc mixing ratio               
LOGICAL,            INTENT(OUT) :: OUSERR ! use Rr mixing ratio               
LOGICAL,            INTENT(OUT) :: OUSERI ! use Ri mixing ratio               
LOGICAL,            INTENT(OUT) :: OUSECI ! use Ci concentration of Ice cristals
LOGICAL,            INTENT(OUT) :: OUSERS ! use Rs mixing ratio               
LOGICAL,            INTENT(OUT) :: OUSERG ! use Rg mixing ratio               
LOGICAL,            INTENT(OUT) :: OUSERH ! use Rh mixing ratio               
LOGICAL,            INTENT(OUT) :: OUSECHEM ! Chemical flag
LOGICAL,            INTENT(OUT) :: OUSECHAQ ! Aqueous Chemical flag
LOGICAL,            INTENT(OUT) :: OUSECHIC ! Ice phase Chemical flag
LOGICAL,            INTENT(OUT) :: OCH_PH   ! pH flag
LOGICAL,            INTENT(OUT) :: OCH_CONV_LINOX ! LiNOX flag
LOGICAL,            INTENT(OUT) :: OLG      ! lagrangian flag
INTEGER,            INTENT(OUT) :: KRIMX, KRIMY        ! number of points for the
                       ! horizontal relaxation for the outermost verticals
INTEGER,            INTENT(OUT) :: KSV_USER    ! number of additional scalar
                                          ! variables in FMfile 
CHARACTER (LEN=4),  INTENT(OUT) :: HTURB  ! Kind of turbulence parameterization
                                          ! used to produce the FMfile
CHARACTER (LEN=4),  INTENT(OUT) :: HTOM   ! Kind of third order moment
LOGICAL,            INTENT(OUT) :: ORMC01 ! flag for RMC01 SBL computations
CHARACTER (LEN=4),  INTENT(OUT) :: HRAD   ! Kind of radiation scheme
CHARACTER (LEN=4),  INTENT(OUT) :: HDCONV ! Kind of deep convection scheme
CHARACTER (LEN=4),  INTENT(OUT) :: HSCONV ! Kind of shallow convection scheme
CHARACTER (LEN=4),  INTENT(OUT) :: HCLOUD ! Kind of microphysical scheme
CHARACTER (LEN=4),  INTENT(OUT) :: HELEC  ! Kind of electrical scheme
CHARACTER (LEN=*),  INTENT(OUT) :: HEQNSYS! type of equations' system
LOGICAL,            INTENT(OUT) :: OSALT    ! Sea Salt flag
LOGICAL,            INTENT(OUT) :: OPASPOL  ! Passive pollutant flag
#ifdef MNH_FOREFIRE
LOGICAL,            INTENT(OUT) :: OFOREFIRE ! ForeFire flag
#endif
LOGICAL,            INTENT(OUT) :: OLNOX_EXPLICIT ! explicit LNOx flag
LOGICAL,            INTENT(OUT) :: OCONDSAMP! Conditional sampling flag
LOGICAL,            INTENT(OUT) :: OBLOWSNOW! Blowing snow flag
LOGICAL,            INTENT(OUT) :: ODUST    ! Dust flag
LOGICAL,            INTENT(OUT) :: OORILAM  ! Dust flag
LOGICAL,            INTENT(OUT) :: OCHTRANS ! Deep convection on scalar
                                            ! variables flag 
LOGICAL,DIMENSION(JPMODELMAX),INTENT(OUT) :: ODEPOS_DST    ! Dust Wet Deposition flag
LOGICAL,DIMENSION(JPMODELMAX),INTENT(OUT) :: ODEPOS_SLT    ! Sea Salt Wet Deposition flag
LOGICAL,DIMENSION(JPMODELMAX),INTENT(OUT) :: ODEPOS_AER    ! Aerosols Wet Deposition flag
!
!*       0.2   declarations of local variables 
!
INTEGER :: ILUDES, & ! logical unit numbers of
           ILUOUT    ! DESFM file and output listing
LOGICAL :: GFOUND          ! Return code when searching namelist
LOGICAL,DIMENSION(JPMODELMAX),SAVE :: LTEMPDEPOS_DST ! Dust Moist flag
LOGICAL,DIMENSION(JPMODELMAX),SAVE :: LTEMPDEPOS_SLT ! Sea Salt Moist flag
LOGICAL,DIMENSION(JPMODELMAX),SAVE :: LTEMPDEPOS_AER ! Orilam Moist flag
!
!-------------------------------------------------------------------------------
!
!*       1.    READ DESFM FILE
!              ---------------
!
CALL PRINT_MSG(NVERB_DEBUG,'IO','READ_DESFM_n','called for '//TRIM(TPDATAFILE%CNAME))
!
IF (.NOT.ASSOCIATED(TPDATAFILE%TDESFILE)) &
  CALL PRINT_MSG(NVERB_FATAL,'IO','READ_DESFM_n','TDESFILE not associated for '//TRIM(TPDATAFILE%CNAME))
!
ILUDES = TPDATAFILE%TDESFILE%NLU
ILUOUT = TLUOUT%NLU
!
CALL POSNAM(ILUDES,'NAM_LUNITN',GFOUND)
CALL INIT_NAM_LUNITN
IF (GFOUND) THEN 
  READ(UNIT=ILUDES,NML=NAM_LUNITn)
  CALL UPDATE_NAM_LUNITN
END IF
CALL POSNAM(ILUDES,'NAM_CONFN',GFOUND)
CALL INIT_NAM_CONFN
IF (GFOUND) THEN 
  READ(UNIT=ILUDES,NML=NAM_CONFn)
  CALL UPDATE_NAM_CONFN
END IF
CALL POSNAM(ILUDES,'NAM_DYNN',GFOUND)
CALL INIT_NAM_DYNN
IF (GFOUND) THEN 
  READ(UNIT=ILUDES,NML=NAM_DYNn)
  CALL UPDATE_NAM_DYNN
END IF
CALL POSNAM(ILUDES,'NAM_ADVN',GFOUND)
CALL INIT_NAM_ADVN
IF (GFOUND) THEN 
  READ(UNIT=ILUDES,NML=NAM_ADVn)
  CALL UPDATE_NAM_ADVN
END IF
CALL POSNAM(ILUDES,'NAM_PARAMN',GFOUND)
CALL INIT_NAM_PARAMn
IF (GFOUND) THEN 
  READ(UNIT=ILUDES,NML=NAM_PARAMn)
  CALL UPDATE_NAM_PARAMn
END IF
CALL POSNAM(ILUDES,'NAM_PARAM_RADN',GFOUND)
CALL INIT_NAM_PARAM_RADn
IF (GFOUND) THEN 
  READ(UNIT=ILUDES,NML=NAM_PARAM_RADn)
  CALL UPDATE_NAM_PARAM_RADn
END IF
#ifdef MNH_ECRAD
CALL POSNAM(ILUDES,'NAM_PARAM_ECRADN',GFOUND)
CALL INIT_NAM_PARAM_ECRADn
IF (GFOUND) THEN 
  READ(UNIT=ILUDES,NML=NAM_PARAM_ECRADn)
  CALL UPDATE_NAM_PARAM_ECRADn
END IF
#endif
CALL POSNAM(ILUDES,'NAM_PARAM_KAFRN',GFOUND)
CALL INIT_NAM_PARAM_KAFRn
IF (GFOUND) THEN 
  READ(UNIT=ILUDES,NML=NAM_PARAM_KAFRn)
  CALL UPDATE_NAM_PARAM_KAFRn
END IF
CALL POSNAM(ILUDES,'NAM_PARAM_MFSHALLN',GFOUND)
CALL INIT_NAM_PARAM_MFSHALLn
IF (GFOUND) THEN 
  READ(UNIT=ILUDES,NML=NAM_PARAM_MFSHALLn)
  CALL UPDATE_NAM_PARAM_MFSHALLn
END IF
CALL POSNAM(ILUDES,'NAM_LBCN',GFOUND)
CALL INIT_NAM_LBCn
IF (GFOUND) THEN 
  READ(UNIT=ILUDES,NML=NAM_LBCn)
  CALL UPDATE_NAM_LBCn
END IF
CALL POSNAM(ILUDES,'NAM_NUDGINGN',GFOUND)
CALL INIT_NAM_NUDGINGn
IF (GFOUND) THEN 
  READ(UNIT=ILUDES,NML=NAM_NUDGINGn)
  CALL UPDATE_NAM_NUDGINGn
END IF
CALL POSNAM(ILUDES,'NAM_TURBN',GFOUND)
CALL INIT_NAM_TURBn
IF (GFOUND) THEN 
  READ(UNIT=ILUDES,NML=NAM_TURBn)
  CALL UPDATE_NAM_TURBn
END IF
CALL POSNAM(ILUDES,'NAM_CH_MNHCN',GFOUND)
CALL INIT_NAM_CH_MNHCn
IF (GFOUND) THEN 
  READ(UNIT=ILUDES,NML=NAM_CH_MNHCn)
  CALL UPDATE_NAM_CH_MNHCn
END IF
CALL POSNAM(ILUDES,'NAM_CH_SOLVERn',GFOUND)
CALL INIT_NAM_CH_SOLVERn
IF (GFOUND) THEN 
  READ(UNIT=ILUDES,NML=NAM_CH_SOLVERn)
  CALL UPDATE_NAM_CH_SOLVERn
END IF
CALL POSNAM(ILUDES,'NAM_DRAGN',GFOUND)
CALL INIT_NAM_DRAGn
IF (GFOUND) THEN 
  READ(UNIT=ILUDES,NML=NAM_DRAGn)
  CALL UPDATE_NAM_DRAGn
END IF
CALL POSNAM(ILUDES,'NAM_IBM_PARAMN',GFOUND,ILUOUT)
CALL INIT_NAM_IBM_PARAMn
IF (GFOUND) THEN
  READ(UNIT=ILUDES,NML=NAM_IBM_PARAMn)
  CALL UPDATE_NAM_IBM_PARAMn
END IF
CALL POSNAM(ILUDES,'NAM_RECYCL_PARAMN',GFOUND,ILUOUT)
CALL INIT_NAM_RECYCL_PARAMn
IF (GFOUND) THEN
  READ(UNIT=ILUDES,NML=NAM_RECYCL_PARAMn)
  CALL UPDATE_NAM_RECYCL_PARAMn
END IF
CALL POSNAM(ILUDES,'NAM_SERIESN',GFOUND,ILUOUT)
CALL INIT_NAM_SERIESn
IF (GFOUND) THEN 
  READ(UNIT=ILUDES,NML=NAM_SERIESn)
  CALL UPDATE_NAM_SERIESn
END IF
CALL POSNAM(ILUDES,'NAM_BLOWSNOWn',GFOUND,ILUOUT)
CALL INIT_NAM_BLOWSNOWn
IF (GFOUND) THEN
  READ(UNIT=ILUDES,NML=NAM_BLOWSNOWn)
  CALL UPDATE_NAM_BLOWSNOWn
END IF
CALL POSNAM(ILUDES,'NAM_BLANKN',GFOUND,ILUOUT)
CALL INIT_NAM_BLANKn
IF (GFOUND) THEN
  READ(UNIT=ILUDES,NML=NAM_BLANKn)
  CALL UPDATE_NAM_BLANKn
END IF
CALL POSNAM(ILUDES,'NAM_STATIONN',GFOUND,ILUOUT)
CALL INIT_NAM_STATIONn
IF (GFOUND) THEN
  READ(UNIT=ILUDES,NML=NAM_STATIONn)
  CALL UPDATE_NAM_STATIONn
END IF
!
!
IF (KMI == 1) THEN
  CALL POSNAM(ILUDES,'NAM_CONF',GFOUND)
  IF (GFOUND) READ(UNIT=ILUDES,NML=NAM_CONF)
  CALL POSNAM(ILUDES,'NAM_DYN',GFOUND)
  IF (GFOUND) READ(UNIT=ILUDES,NML=NAM_DYN)
  CALL POSNAM(ILUDES,'NAM_NESTING',GFOUND)
  IF (GFOUND) READ(UNIT=ILUDES,NML=NAM_NESTING)
  CALL POSNAM(ILUDES,'NAM_BACKUP',GFOUND)
  IF (GFOUND) THEN
    IF (.NOT.ALLOCATED(XBAK_TIME)) THEN
      ALLOCATE(XBAK_TIME(NMODEL,JPOUTMAX))
      XBAK_TIME(:,:) = XNEGUNDEF
    END IF
    IF (.NOT.ALLOCATED(XOUT_TIME)) THEN
      ALLOCATE(XOUT_TIME(NMODEL,JPOUTMAX)) !Allocate *OUT* variables to prevent
      XOUT_TIME(:,:) = XNEGUNDEF
    END IF
    IF (.NOT.ALLOCATED(NBAK_STEP)) THEN
      ALLOCATE(NBAK_STEP(NMODEL,JPOUTMAX))
      NBAK_STEP(:,:) = NNEGUNDEF
    END IF
    IF (.NOT.ALLOCATED(NOUT_STEP)) THEN
      ALLOCATE(NOUT_STEP(NMODEL,JPOUTMAX)) !problems if NAM_OUTPUT does not exist
      NOUT_STEP(:,:) = NNEGUNDEF
    END IF
    IF (.NOT.ALLOCATED(COUT_VAR)) THEN
      ALLOCATE(COUT_VAR (NMODEL,JPOUTVARMAX))
      COUT_VAR(:,:)  = ''
    END IF
    READ(UNIT=ILUDES,NML=NAM_BACKUP)
  ELSE
    CALL POSNAM(ILUDES,'NAM_FMOUT',GFOUND)
    IF (GFOUND) CALL PRINT_MSG(NVERB_FATAL,'IO','READ_DESFM_n','use namelist NAM_BACKUP instead of namelist NAM_FMOUT')
  END IF
  CALL POSNAM(ILUDES,'NAM_OUTPUT',GFOUND)
  IF (GFOUND) THEN
    IF (.NOT.ALLOCATED(XBAK_TIME)) THEN
      ALLOCATE(XBAK_TIME(NMODEL,JPOUTMAX)) !Allocate *BAK* variables to prevent
      XBAK_TIME(:,:) = XNEGUNDEF
    END IF
    IF (.NOT.ALLOCATED(XOUT_TIME)) THEN
      ALLOCATE(XOUT_TIME(NMODEL,JPOUTMAX))
      XOUT_TIME(:,:) = XNEGUNDEF
    END IF
    IF (.NOT.ALLOCATED(NBAK_STEP)) THEN
      ALLOCATE(NBAK_STEP(NMODEL,JPOUTMAX)) !problems if NAM_BACKUP does not exist
      NBAK_STEP(:,:) = NNEGUNDEF
    END IF
    IF (.NOT.ALLOCATED(NOUT_STEP)) THEN
      ALLOCATE(NOUT_STEP(NMODEL,JPOUTMAX))
      NOUT_STEP(:,:) = NNEGUNDEF
    END IF
    IF (.NOT.ALLOCATED(COUT_VAR)) THEN
      ALLOCATE(COUT_VAR (NMODEL,JPOUTVARMAX))
      COUT_VAR(:,:)  = ''
    END IF
    READ(UNIT=ILUDES,NML=NAM_OUTPUT)
  END IF
! Note: it is not useful to read the budget namelists in the .des files
! The value here (if present in file) don't need to be compared with the ones in the EXSEGn files
!   CALL POSNAM(ILUDES,'NAM_BUDGET',GFOUND)
!   IF (GFOUND) READ(UNIT=ILUDES,NML=NAM_BUDGET)
!   CALL POSNAM(ILUDES,'NAM_BU_RU',GFOUND)
!   IF (GFOUND) READ(UNIT=ILUDES,NML=NAM_BU_RU)
!   CALL POSNAM(ILUDES,'NAM_BU_RV',GFOUND)
!   IF (GFOUND) READ(UNIT=ILUDES,NML=NAM_BU_RV)
!   CALL POSNAM(ILUDES,'NAM_BU_RW',GFOUND)
!   IF (GFOUND) READ(UNIT=ILUDES,NML=NAM_BU_RW)
!   CALL POSNAM(ILUDES,'NAM_BU_RTH',GFOUND)
!   IF (GFOUND) READ(UNIT=ILUDES,NML=NAM_BU_RTH)
!   CALL POSNAM(ILUDES,'NAM_BU_RTKE',GFOUND)
!   IF (GFOUND) READ(UNIT=ILUDES,NML=NAM_BU_RTKE)
!   CALL POSNAM(ILUDES,'NAM_BU_RRV',GFOUND)
!   IF (GFOUND) READ(UNIT=ILUDES,NML=NAM_BU_RRV)
!   CALL POSNAM(ILUDES,'NAM_BU_RRC',GFOUND)
!   IF (GFOUND) READ(UNIT=ILUDES,NML=NAM_BU_RRC)
!   CALL POSNAM(ILUDES,'NAM_BU_RRR',GFOUND)
!   IF (GFOUND) READ(UNIT=ILUDES,NML=NAM_BU_RRR)
!   CALL POSNAM(ILUDES,'NAM_BU_RRI',GFOUND)
!   IF (GFOUND) READ(UNIT=ILUDES,NML=NAM_BU_RRI)
!   CALL POSNAM(ILUDES,'NAM_BU_RRS',GFOUND)
!   IF (GFOUND) READ(UNIT=ILUDES,NML=NAM_BU_RRS)
!   CALL POSNAM(ILUDES,'NAM_BU_RRG',GFOUND)
!   IF (GFOUND) READ(UNIT=ILUDES,NML=NAM_BU_RRG)
!   CALL POSNAM(ILUDES,'NAM_BU_RRH',GFOUND)
!   IF (GFOUND) READ(UNIT=ILUDES,NML=NAM_BU_RRH)
!   CALL POSNAM(ILUDES,'NAM_BU_RSV',GFOUND)
!   IF (GFOUND) READ(UNIT=ILUDES,NML=NAM_BU_RSV)
  CALL POSNAM(ILUDES,'NAM_LES',GFOUND)
  IF (GFOUND) READ(UNIT=ILUDES,NML=NAM_LES)
  CALL POSNAM(ILUDES,'NAM_PDF',GFOUND)
  IF (GFOUND) READ(UNIT=ILUDES,NML=NAM_PDF)
  CALL POSNAM(ILUDES,'NAM_FRC',GFOUND)
  IF (GFOUND) READ(UNIT=ILUDES,NML=NAM_FRC)
  CALL POSNAM(ILUDES,'NAM_PARAM_ICE',GFOUND)
  IF (GFOUND) READ(UNIT=ILUDES,NML=NAM_PARAM_ICE)
  CALL POSNAM(ILUDES,'NAM_PARAM_C2R2',GFOUND)
  IF (GFOUND) READ(UNIT=ILUDES,NML=NAM_PARAM_C2R2)
  CALL POSNAM(ILUDES,'NAM_PARAM_C1R3',GFOUND)
  IF (GFOUND) READ(UNIT=ILUDES,NML=NAM_PARAM_C1R3)
  CALL POSNAM(ILUDES,'NAM_PARAM_LIMA',GFOUND)
  IF (GFOUND) READ(UNIT=ILUDES,NML=NAM_PARAM_LIMA)
  CALL POSNAM(ILUDES,'NAM_ELEC',GFOUND)
  IF (GFOUND) READ(UNIT=ILUDES,NML=NAM_ELEC)
  CALL POSNAM(ILUDES,'NAM_SERIES',GFOUND,ILUOUT)
  IF (GFOUND) READ(UNIT=ILUDES,NML=NAM_SERIES)
  CALL POSNAM(ILUDES,'NAM_TURB_CLOUD',GFOUND,ILUOUT)
  IF (GFOUND) READ(UNIT=ILUDES,NML=NAM_TURB_CLOUD)
  CALL POSNAM(ILUDES,'NAM_TURB',GFOUND,ILUOUT)
  IF (GFOUND) READ(UNIT=ILUDES,NML=NAM_TURB)
  CALL POSNAM(ILUDES,'NAM_CH_ORILAM',GFOUND,ILUOUT)
  IF (GFOUND) READ(UNIT=ILUDES,NML=NAM_CH_ORILAM)
  CALL POSNAM(ILUDES,'NAM_DUST',GFOUND,ILUOUT)
  IF (GFOUND) READ(UNIT=ILUDES,NML=NAM_DUST)
  CALL POSNAM(ILUDES,'NAM_SALT',GFOUND,ILUOUT)
  IF (GFOUND) READ(UNIT=ILUDES,NML=NAM_SALT)
  CALL POSNAM(ILUDES,'NAM_PASPOL',GFOUND,ILUOUT)
  IF (GFOUND) READ(UNIT=ILUDES,NML=NAM_PASPOL)
#ifdef MNH_FOREFIRE
  CALL POSNAM(ILUDES,'NAM_FOREFIRE',GFOUND,ILUOUT)
  IF (GFOUND) READ(UNIT=ILUDES,NML=NAM_FOREFIRE)
#endif
  CALL POSNAM(ILUDES,'NAM_CONDSAMP',GFOUND,ILUOUT)
  IF (GFOUND) READ(UNIT=ILUDES,NML=NAM_CONDSAMP)
  CALL POSNAM(ILUDES,'NAM_BLOWSNOW',GFOUND,ILUOUT)
  IF (GFOUND) READ(UNIT=ILUDES,NML=NAM_BLOWSNOW)
  CALL POSNAM(ILUDES,'NAM_2D_FRC',GFOUND,ILUOUT)
  IF (GFOUND) READ(UNIT=ILUDES,NML=NAM_2D_FRC)
  LTEMPDEPOS_DST(:) = LDEPOS_DST(:)
  LTEMPDEPOS_SLT(:) = LDEPOS_SLT(:)
  LTEMPDEPOS_AER(:) = LDEPOS_AER(:)
  CALL POSNAM(ILUDES,'NAM_LATZ_EDFLX',GFOUND)
  IF (GFOUND) READ(UNIT=ILUDES,NML=NAM_LATZ_EDFLX)
  CALL POSNAM(ILUDES,'NAM_VISC',GFOUND,ILUOUT)
  IF (GFOUND) READ(UNIT=ILUDES,NML=NAM_VISC)
END IF                                                       
!
!-------------------------------------------------------------------------------
!
!*       2.    SAVE SOME FMFILE ATTRIBUTES
!              ---------------------------
HCONF  = CCONF
OFLAT  = LFLAT
OUSERV = LUSERV
OUSERC = LUSERC
OUSERR = LUSERR
OUSERI = LUSERI
OUSECI = LUSECI
OUSERS = LUSERS
OUSERG = LUSERG
OUSERH = LUSERH
OUSECHEM = LUSECHEM
OUSECHAQ = LUSECHAQ
OUSECHIC = LUSECHIC
OCH_PH = LCH_PH
OCH_CONV_LINOX = LCH_CONV_LINOX
ODUST    = LDUST
ODEPOS_DST(KMI)    = LTEMPDEPOS_DST(KMI)
ODEPOS_SLT(KMI)    = LTEMPDEPOS_SLT(KMI)
ODEPOS_AER(KMI)    = LTEMPDEPOS_AER(KMI)
OCHTRANS = LCHTRANS
OSALT    = LSALT
OORILAM  = LORILAM
OLG      = LLG
OPASPOL  = LPASPOL
#ifdef MNH_FOREFIRE
OFOREFIRE  = LFOREFIRE
#endif
OLNOX_EXPLICIT = LLNOX_EXPLICIT
OCONDSAMP= LCONDSAMP
OBLOWSNOW= LBLOWSNOW
! Initially atmosphere free of blowing snow particles
IF(KMI>1) OBLOWSNOW=.FALSE.
KRIMX  = NRIMX
KRIMY  = NRIMY
KSV_USER = NSV_USER
HTURB  = CTURB
HTOM = CTOM
ORMC01 = LRMC01
HRAD   = CRAD
HDCONV = CDCONV
HSCONV = CSCONV
HCLOUD = CCLOUD
HELEC  = CELEC
HEQNSYS = CEQNSYS
!
!-------------------------------------------------------------------------------
!
!*       3.    WRITE DESFM ON OUTPUT LISTING
!              ------------------------------
!
IF (NVERB >= 10) THEN
  WRITE(UNIT=ILUOUT,FMT="(/,'DESCRIPTOR OF INITIAL FILE FOR MODEL ',I2)") KMI
  WRITE(UNIT=ILUOUT,FMT="(  '------------------------------------ '   )")
!  
  WRITE(UNIT=ILUOUT,FMT="('********** LOGICAL UNITSn **********')")
  WRITE(UNIT=ILUOUT,NML=NAM_LUNITn)
!  
  WRITE(UNIT=ILUOUT,FMT="('********** CONFIGURATIONn **********')")
  WRITE(UNIT=ILUOUT,NML=NAM_CONFn)
!  
  WRITE(UNIT=ILUOUT,FMT="('********** DYNAMICn ****************')")
  WRITE(UNIT=ILUOUT,NML=NAM_DYNn)
!  
  WRITE(UNIT=ILUOUT,FMT="('********** ADVECTIONn **************')")
  WRITE(UNIT=ILUOUT,NML=NAM_ADVn)
!  
  WRITE(UNIT=ILUOUT,FMT="('********** PARAMETERIZATIONSn ******')")
  WRITE(UNIT=ILUOUT,NML=NAM_PARAMn)
!  
  WRITE(UNIT=ILUOUT,FMT="('********** RADIATIONSn *************')")
  WRITE(UNIT=ILUOUT,NML=NAM_PARAM_RADn)
!  
#ifdef MNH_ECRAD
  WRITE(UNIT=ILUOUT,FMT="('********** ECRAD RADIATIONSn *************')")
  WRITE(UNIT=ILUOUT,NML=NAM_PARAM_ECRADn)
#endif
!  
  WRITE(UNIT=ILUOUT,FMT="('********** DEEP CONVECTIONn ********')")
  WRITE(UNIT=ILUOUT,NML=NAM_PARAM_KAFRn)
!  
  WRITE(UNIT=ILUOUT,FMT="('*** MASS FLUX SHALLOW CONVECTION ***')")
  WRITE(UNIT=ILUOUT,NML=NAM_PARAM_MFSHALLn)
!
  WRITE(UNIT=ILUOUT,FMT="('********** LBCn ********************')")
  WRITE(UNIT=ILUOUT,NML=NAM_LBCn)
!  
  WRITE(UNIT=ILUOUT,FMT="('********** TURBn *******************')")  
  WRITE(UNIT=ILUOUT,NML=NAM_TURBn)
!
  WRITE(UNIT=ILUOUT,FMT="('********** DRAGn *******************')")  
  WRITE(UNIT=ILUOUT,NML=NAM_DRAGn)
!
  WRITE(UNIT=ILUOUT,FMT="('********** IBM FORCING *************')")
  WRITE(UNIT=ILUOUT,NML=NAM_IBM_PARAMn)  
!
  WRITE(UNIT=ILUOUT,FMT="('********** RECYLING *************')")
  WRITE(UNIT=ILUOUT,NML=NAM_RECYCL_PARAMn)
!
  WRITE(UNIT=ILUOUT,FMT="('********** NUDGINGn ****************')")  
  WRITE(UNIT=ILUOUT,NML=NAM_NUDGINGn)
!  
  WRITE(UNIT=ILUOUT,FMT="('********** CHEMICAL MONITORn *******')")  
  WRITE(UNIT=ILUOUT,NML=NAM_CH_MNHCn)
!
  WRITE(UNIT=ILUOUT,FMT="('********** CHEMICAL SOLVER *********')")
  WRITE(UNIT=ILUOUT,NML=NAM_CH_SOLVERn)
!
  WRITE(UNIT=ILUOUT,FMT="('********** BLOWSNOWn ***************')")
  WRITE(UNIT=ILUOUT,NML=NAM_BLOWSNOWn)
!
  WRITE(UNIT=ILUOUT,FMT="('********** BLANKn ******************')")
  WRITE(UNIT=ILUOUT,NML=NAM_BLANKn)
!
  WRITE(UNIT=ILUOUT,FMT="('********** STATIONn ******************')")
  WRITE(UNIT=ILUOUT,NML=NAM_STATIONn)
!
  IF (KMI==1) THEN
    WRITE(UNIT=ILUOUT,FMT="(/,'PART OF INITIAL FILE COMMON TO ALL THE MODELS')")
    WRITE(UNIT=ILUOUT,FMT="(  '---------------------------------------------')")
!
    WRITE(UNIT=ILUOUT,FMT="('************ CONFIGURATION ********************')")
    WRITE(UNIT=ILUOUT,NML=NAM_CONF)
!    
    WRITE(UNIT=ILUOUT,FMT="('************ DYNAMIC **************************')")
    WRITE(UNIT=ILUOUT,NML=NAM_DYN)
!    
! Budget namelists not read anymore in READ_DESFM_n
!     WRITE(UNIT=ILUOUT,FMT="('************ BUDGET ***************************')")
!     WRITE(UNIT=ILUOUT,NML=NAM_BUDGET)
! !
!     WRITE(UNIT=ILUOUT,FMT="('************ U BUDGET *************************')")
!     WRITE(UNIT=ILUOUT,NML=NAM_BU_RU)
! !
!     WRITE(UNIT=ILUOUT,FMT="('************ V BUDGET *************************')")
!     WRITE(UNIT=ILUOUT,NML=NAM_BU_RV)
! !
!     WRITE(UNIT=ILUOUT,FMT="('************ W BUDGET *************************')")
!     WRITE(UNIT=ILUOUT,NML=NAM_BU_RW)
! !
!     WRITE(UNIT=ILUOUT,FMT="('************ TH BUDGET ************************')")
!     WRITE(UNIT=ILUOUT,NML=NAM_BU_RTH)
! !
!     WRITE(UNIT=ILUOUT,FMT="('************ TKE BUDGET ***********************')")
!     WRITE(UNIT=ILUOUT,NML=NAM_BU_RTKE)
! !
!     WRITE(UNIT=ILUOUT,FMT="('************ RV BUDGET ************************')")
!     WRITE(UNIT=ILUOUT,NML=NAM_BU_RRV)
! !
!     WRITE(UNIT=ILUOUT,FMT="('************ RC BUDGET ************************')")
!     WRITE(UNIT=ILUOUT,NML=NAM_BU_RRC)
! !
!     WRITE(UNIT=ILUOUT,FMT="('************ RR BUDGET ************************')")
!     WRITE(UNIT=ILUOUT,NML=NAM_BU_RRR)
! !
!     WRITE(UNIT=ILUOUT,FMT="('************ RI BUDGET ************************')")
!     WRITE(UNIT=ILUOUT,NML=NAM_BU_RRI)
! !
!     WRITE(UNIT=ILUOUT,FMT="('************ RS BUDGET ************************')")
!     WRITE(UNIT=ILUOUT,NML=NAM_BU_RRS)
! !
!     WRITE(UNIT=ILUOUT,FMT="('************ RG BUDGET ************************')")
!     WRITE(UNIT=ILUOUT,NML=NAM_BU_RRG)
! !
!     WRITE(UNIT=ILUOUT,FMT="('************ RH BUDGET ************************')")
!     WRITE(UNIT=ILUOUT,NML=NAM_BU_RRH)
! !
!     WRITE(UNIT=ILUOUT,FMT="('************ SVx BUDGET ***********************')")
!     WRITE(UNIT=ILUOUT,NML=NAM_BU_RSV)
!    
    WRITE(UNIT=ILUOUT,FMT="('************ LES ******************************')")
    WRITE(UNIT=ILUOUT,NML=NAM_LES)
!
    WRITE(UNIT=ILUOUT,FMT="('************ PDF ******************************')")
    WRITE(UNIT=ILUOUT,NML=NAM_PDF)
!    
    WRITE(UNIT=ILUOUT,FMT="('************ FORCING **************************')")
    WRITE(UNIT=ILUOUT,NML=NAM_FRC)
!    
    WRITE(UNIT=ILUOUT,FMT="('************ ICE SCHEME ***********************')")
    WRITE(UNIT=ILUOUT,NML=NAM_PARAM_ICE)
!    
    WRITE(UNIT=ILUOUT,FMT="('************ ORILAM SCHEME ********************')")
    WRITE(UNIT=ILUOUT,NML=NAM_CH_ORILAM)
!
    WRITE(UNIT=ILUOUT,FMT="('************ SALT SCHEME **********************')")
    WRITE(UNIT=ILUOUT,NML=NAM_SALT)
!
    WRITE(UNIT=ILUOUT,FMT="('************ DUST SCHEME **********************')")
    WRITE(UNIT=ILUOUT,NML=NAM_DUST)
!    
    WRITE(UNIT=ILUOUT,FMT="('************ PASSIVE POLLUTANT  ***************')")
    WRITE(UNIT=ILUOUT,NML=NAM_PASPOL)
!
    WRITE(UNIT=ILUOUT,FMT="('************ VISCOSITY  ***************')")
    WRITE(UNIT=ILUOUT,NML=NAM_VISC)
!
#ifdef MNH_FOREFIRE
	WRITE(UNIT=ILUOUT,FMT="('************ FOREFIRE  ***************')")
	WRITE(UNIT=ILUOUT,NML=NAM_FOREFIRE)
!
#endif		
    WRITE(UNIT=ILUOUT,FMT="('************ CONDITIONAL SAMPLING *************')")
    WRITE(UNIT=ILUOUT,NML=NAM_CONDSAMP)
 !
    WRITE(UNIT=ILUOUT,FMT="('********** BLOWING SNOW SCHEME******************')")
    WRITE(UNIT=ILUOUT,NML=NAM_BLOWSNOW)
!
    IF( CCLOUD == 'C2R2' ) THEN
      WRITE(UNIT=ILUOUT,FMT="('************ C2R2 SCHEME **********************')")
      WRITE(UNIT=ILUOUT,NML=NAM_PARAM_C2R2)
    END IF
!
    IF( CCLOUD == 'KHKO' ) THEN                                                    !modif
      WRITE(UNIT=ILUOUT,FMT="('************ KHKO SCHEME **********************')")
      WRITE(UNIT=ILUOUT,NML=NAM_PARAM_C2R2)
    END IF
!
    IF( CCLOUD == 'C3R5' ) THEN
      WRITE(UNIT=ILUOUT,FMT="('************ C3R5 SCHEME **********************')")
      WRITE(UNIT=ILUOUT,NML=NAM_PARAM_C2R2)
      WRITE(UNIT=ILUOUT,NML=NAM_PARAM_C1R3)
    END IF
!
    IF( CCLOUD == 'LIMA' ) THEN
      WRITE(UNIT=ILUOUT,FMT="('************ LIMA SCHEME **********************')")
      WRITE(UNIT=ILUOUT,NML=NAM_PARAM_LIMA)
    END IF
!
   IF (CELEC /= 'NONE') THEN
     WRITE(UNIT=ILUOUT,FMT="('************ ELEC SCHEME **********************')")
     WRITE(UNIT=ILUOUT,NML=NAM_ELEC)                                             
   END IF
!    
  END IF
!
END IF  
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE READ_DESFM_n
