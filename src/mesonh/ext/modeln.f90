!MNH_LIC Copyright 1994-2025 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ###################
      MODULE MODI_MODEL_n
!     ###################
!
INTERFACE
!
       SUBROUTINE MODEL_n( KTCOUNT, TPBAKFILE, TPDTMODELN, OEXIT )
!
USE MODD_IO,        ONLY: TFILEDATA
USE MODD_TYPE_DATE, ONLY: DATE_TIME
!
INTEGER,                  INTENT(IN)    :: KTCOUNT    ! Temporal loop index of model KMODEL
TYPE(TFILEDATA), POINTER, INTENT(OUT)   :: TPBAKFILE  ! Pointer for backup file
TYPE(DATE_TIME),          INTENT(OUT)   :: TPDTMODELN ! Time of current model computation
LOGICAL,                  INTENT(INOUT) :: OEXIT      ! Switch for the end of the temporal loop
!
END SUBROUTINE MODEL_n
!
END INTERFACE
!
END MODULE MODI_MODEL_n

!     ################################### 
      SUBROUTINE MODEL_n( KTCOUNT, TPBAKFILE, TPDTMODELN, OEXIT )
!     ###################################
!
!!****  *MODEL_n * -monitor of the model version _n 
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to build up a typical model version
!     by sequentially calling the specialized routines.
!
!!**  METHOD
!!    ------
!!      Some preliminary initializations are performed in the first section.
!!    Then, specialized routines are called to update the guess of the future
!!    instant XRxxS of the variable xx by adding the effects of all the
!!    different sources of evolution.
!!
!!              (guess of xx at t+dt) * Rhod_ref * Jacobian
!!      XRxxS = -------------------------------------------
!!                           2 dt
!!
!!      At this level, the informations are transferred with a USE association
!!    from the INIT step, where the modules have been previously filled. The
!!    transfer to the subroutines computing each source term is performed by
!!    argument in order to avoid repeated compilations of these subroutines.
!!      This monitor model_n, must therefore be duplicated for each model,
!!    model1 corresponds in this case to the outermost model, model2 is used
!!    for the first level of gridnesting,....  
!!      The effect of all parameterizations is computed in PHYS_PARAM_n, which
!!    is itself a monitor. This is due to a possible large number of
!!    parameterizations, which can be activated and therefore, will require a
!!    very large list of arguments. To circumvent this problem, we transfer by
!!    a USE association, the necessary informations in this monitor, which will
!!    dispatch the pertinent information to every parametrization.
!!      Some elaborated diagnostics, LES tools, budget storages are also called
!!    at this level because they require informations about the fields at every
!!    timestep.
!!
!!
!!    EXTERNAL
!!    --------
!!      Subroutine IO_File_open: to open a file
!!      Subroutine WRITE_DESFM: to write the descriptive part of a FMfile
!!      Subroutine WRITE_LFIFM: to write the binary part of a FMfile
!!      Subroutine SET_MASK   : to compute all the masks selected for budget
!!                         computations
!!      Subroutine BOUNDARIES   : set the fields at the marginal points in every
!!                         directions according the selected boundary conditions
!!      Subroutine INITIAL_GUESS: initializes the guess of the future instant
!!      Subroutine LES_FLX_SPECTRA: computes the resolved fluxes and the
!!                     spectra of some quantities when running in LES mode.
!!      Subroutine ADVECTION: computes the advection terms.
!!      Subroutine DYN_SOURCES: computes the curvature, Coriolis, gravity terms.
!!      Subroutine NUM_DIFF: applies the fourth order numerical diffusion.
!!      Subroutine RELAXATION: performs the relaxation to Larger Scale fields
!!                             in the upper levels and outermost vertical planes
!!      Subroutine PHYS_PARAM_n : computes the parameterized physical terms
!!      Subroutine RAD_BOUND: prepares the velocity normal components for the bc.
!!      Subroutine RESOLVED_CLOUD : computes the sources terms for water in any
!!                                  form
!!      Subroutine PRESSURE : computes the pressure gradient term and the
!!                            absolute pressure
!!      Subroutine EXCHANGE : updates the halo of each subdomains
!!      Subroutine ENDSTEP : advances in time the  fields.
!!      Subroutines UVW_LS_COUPLING and SCALAR_LS_COUPLING:
!!                                 compute the large scale fields, used to
!!                                 couple Model_n with outer informations.
!!      Subroutine ENDSTEP_BUDGET: writes the budget informations.
!!      Subroutine IO_File_close: closes a file
!!      Subroutine DATETIME_CORRECTDATE: transform the current time in GMT
!!      Subroutine FORCING : computes forcing terms
!!      Subroutine ADD3DFIELD_ll : add a field to 3D-list
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!          MODD_DYN
!!          MODD_CONF
!!          MODD_NESTING
!!          MODD_BUDGET
!!          MODD_PARAMETERS
!!          MODD_CONF_n
!!          MODD_CURVCOR_n
!!          MODD_DYN_n
!!          MODD_DIM_n
!!          MODD_ADV_n
!!          MODD_FIELD_n
!!          MODD_LSFIELD_n
!!          MODD_GRID_n
!!          MODD_METRICS_n
!!          MODD_LBC_n
!!          MODD_PARAM_n
!!          MODD_REF_n
!!          MODD_LUNIT_n
!!          MODD_OUT_n
!!          MODD_TIME_n
!!          MODD_TURB_n
!!          MODD_CLOUDPAR_n
!!          MODD_TIME
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      J.-P. Pinty                  * LA *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    15/09/94
!!      Modification 20/10/94  (J.Stein) for the outputs and abs_layers routines
!!      Modification 10/11/94  (J.Stein) change ABS_LAYER_FIELDS call
!!      Modification 16/11/94  (J.Stein) add call to the renormalization
!!      Modification 17/11/94  (J.-P. Lafore and J.-P. Pinty) call NUM_DIFF
!!      Modification 08/12/94  (J.Stein) cleaning + remove (RENORM + ABS_LAYER..
!!                             ..) + add RELAXATION + LS fiels in the arguments
!!      Modification 19/12/94  (J.Stein) switch for the num diff
!!      Modification 22/12/94  (J.Stein) update tdtcur + change dyn_source call
!!      Modification 05/01/95  (J.Stein) add the parameterization monitor
!!      Modification 09/01/95  (J.Stein) add the 1D switch
!!      Modification 10/01/95  (J.Stein) displace the TDTCUR computation
!!      Modification 03/01/95  (J.-P. Lafore) Absolute pressure diagnosis
!!      Modification Jan 19, 1995 (J. Cuxart) Shunt the DYN_SOURCES in 1D cases.
!!      Modification Jan 24, 1995 (J. Stein)  Interchange Boundaries and
!!                           Initial_guess to correct a bug in 2D configuration
!!      Modification Feb 02, 1995 (I.Mallet) update BOUNDARIES and RAD_BOUND
!!                                           calls
!!      Modification Mar 10, 1995 (I.Mallet) add call to SET_COUPLING
!!                   March,21, 1995 (J. Stein) remove R from the historical var.
!!                   March,26, 1995 (J. Stein) add the EPS variable
!!                   April 18, 1995 (J. Cuxart) add the LES call
!!                   Sept 20,1995 (Lafore) coupling for the dry mass Md
!!                   Nov   2,1995 (Stein) displace the temporal counter increase
!!                   Jan   2,1996 (Stein) rm the test on the temporal counter
!!      Modification Feb   5,1996 (J. Vila) implementation new advection
!!                                          schemes for scalars
!!      Modification Feb  20,1996 (J.Stein) doctor norm
!!                   Dec95 - Jul96 (Georgelin, Pinty, Mari, Suhre) FORCING
!!                   June 17,1996 (Vincent, Lafore, Jabouille)
!!                                        statistics of computing time
!!                   Aug 8, 1996 (K. Suhre) add chemistry
!!                   October 12, 1996 (J. Stein) save the PSRC value
!!                   Sept 05,1996 (V.Masson) print of loop index for debugging
!!                                           purposes
!!                   July 22,1996 (Lafore) improve write of computing time statistics
!!                   July 29,1996 (Lafore) nesting introduction
!!                   Aug.  1,1996 (Lafore) synchronization between models
!!                   Sept. 4,1996 (Lafore) modification of call to routine SET_COUPLING
!!                                         now split in 2 routines
!!                                         (UVW_LS_COUPLING and SCALAR_LS_COUPLING)
!!                   Sept  5,1996 (V.Masson) print of loop index for debugging
!!                                           purposes
!!                   Sept 25,1996 (V.Masson) test for coupling performed here
!!                   Oct. 29,1996 (Lafore)   one-way nesting implementation
!!                   Oct. 12,1996 (J. Stein) save the PSRC value
!!                   Dec. 12,1996 (Lafore)   change call to RAD_BOUND
!!                   Dec. 21,1996 (Lafore)   two-way nesting implementation
!!                   Mar. 12,1997 (Lafore)   introduction of "surfacic" LS fields
!!                   Nov 18, 1996 (J.-P. Pinty) FORCING revisited (translation)
!!                   Dec 04, 1996 (J.-P. Pinty) include mixed-phase clouds
!!                   Dec 20, 1996 (J.-P. Pinty) update the budgets
!!                   Dec 23, 1996 (J.-P. Pinty) add the diachronic file control
!!                   Jan 11, 1997 (J.-P. Pinty) add the deep convection control
!!                   Dec  20,1996 (V.Masson) call boundaries before the writing
!!                   Fev 25, 1997 (P.Jabouille) modify the LES tools
!!                   April 3,1997 (Lafore)      merging of the nesting
!!                                              developments on MASTER3
!!                   Jul.  8,1997 (Lafore)  print control for nesting (NVERB>=7)
!!                   Jul. 28,1997 (Masson)  supress LSTEADY_DMASS
!!                   Aug. 19,1997 (Lafore)  full Clark's formulation introduction
!!                   Sept 26,1997 (Lafore)  LS source calculation at restart
!!                                          (temporarily test to have LS at instant t)
!!                   Jan. 28,1998 (Bechtold) add SST forcing
!!                   fev. 10,1998 (Lafore)  RHODJ computation and storage for budget
!!                   Jul. 10,1998 (Stein )  sequentiel loop for nesting
!!                   Apr. 07,1999 (Stein )  cleaning of the nesting subroutines
!!                   oct. 20,1998 (Jabouille) //
!!                   oct. 20,2000 (J.-P. Pinty) add the C2R2 scheme
!!                   fev. 01,2001 (D.Gazen) add module MODD_NSV for NSV variables
!!                   mar,  4,2002 (V.Ducrocq) call to temporal series
!!                   mar, 8, 2001 (V. Masson) advection of perturbation of theta in neutral cases.
!!                   Nov, 6, 2002 (V. Masson) time counters for budgets & LES
!!                   mars 20,2001 (Pinty)   add ICE4 and C3R5 options
!!                   jan. 2004    (Masson)  surface externalization
!!                   sept 2004 (M. Tomasini) Cloud mixing length modification
!!                   june 2005 (P. Tulet)  add aerosols / dusts
!!                   Jul. 2005 (N. Asencio)  two_way and phys_param calls: 
!!                             Add the surface parameters : precipitating 
!!                             hydrometeors, Short and Long Wave , MASKkids array 
!!                   Fev. 2006 (M. Leriche) add aqueous phase chemistry
!!                   april 2006 (T.Maric) Add halo related to 4th order advection scheme
!!                   May 2006 Remove KEPS
!!                   Oct 2008 (C.Lac) FIT for variables advected with PPM
!!                   July 2009 : Displacement of surface diagnostics call to be
!!                               coherent with  surface diagnostics obtained with DIAG
!!                   10/11/2009 (P. Aumond) Add mean moments
!!                   Nov, 12, 2009 (C. Barthe) add cloud electrification and lightning flashes
!!                   July 2010 (M. Leriche) add ice phase chemical species
!!                   April 2011 (C.Lac) : Remove instant M 
!!                   April 2011 (C.Lac, V.Masson) : Time splitting for advection
!!      J.Escobar 21/03/2013: for HALOK comment all NHALO=1 test
!!       P. Tulet      Nov 2014 accumulated moles of aqueous species that fall at the surface   
!!                   Dec 2014 (C.Lac) : For reproducibility START/RESTA
!!      J.Escobar 20/04/2015: missing UPDATE_HALO before UPDATE_HALO2
!!              July, 2015 (O.Nuissier/F.Duffourg) Add microphysics diagnostic for
!!                                      aircraft, ballon and profiler
!!       C.Lac    11/09/2015: correction of the budget due to FIT temporal scheme
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1
!!                   Sep 2015 (S. Bielli) : Remove YDADFILE from argument call 
!                              of write_phys_param
!!      J.Escobar : 19/04/2016 : Pb IOZ/NETCDF , missing OPARALLELIO=.FALSE. for PGD files
!!      M.Mazoyer : 04/2016      DTHRAD used for radiative cooling when LACTIT
!!!      Modification    01/2016  (JP Pinty) Add LIMA
!!  06/2016     (G.Delautier) phasage surfex 8
!!      M.Leriche : 03/2016 Move computation of accumulated chem. in rain to ch_monitor
!!                  09/2016 Add filter on negative values on AERDEP SV before relaxation
!!                  10/2016  (C.Lac) _ Correction on the flag for Strang splitting
!!                                  to insure reproducibility between START and RESTA
!!                                  _  Add OSPLIT_WENO
!!                                  _ Add droplet deposition 
!!                   10/2016 (M.Mazoyer) New KHKO output fields
!!      P.Wautelet : 11/07/2016 : removed MNH_NCWRIT define
!!                   09/2017 Q.Rodier add LTEND_UV_FRC
!!                   10/2017 (C.Lac) Necessity to have chemistry processes as
!!                            the las process modifying XRSVS
!!  01/2018      (G.Delautier) SURFEX 8.1
!!  03/2018     (P.Wautelet)   replace ADD_FORECAST_TO_DATE by DATETIME_CORRECTDATE
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!                  07/2017  (V. Vionnet) : Add blowing snow scheme 
!!      S. Riette : 11/2016 Add ZPABST to keep pressure constant during timestep
!!                   01/2018 (C.Lac) Add VISCOSITY
!!  Philippe Wautelet: 21/01/2019: add LIO_ALLOW_NO_BACKUP and LIO_NO_WRITE to modd_io_ll
!                                  to allow to disable writes (for bench purposes)
!  P. Wautelet 07/02/2019: remove OPARALLELIO argument from open and close files subroutines
!                          (nsubfiles_ioz is now determined in IO_File_add2list)
!!                   02/2019 C.Lac add rain fraction as an output field
!!      Bielli S. 02/2019  Sea salt : significant sea wave height influences salt emission; 5 salt modes
!  P. Wautelet 28/03/2019: use MNHTIME for time measurement variables
!  P. Wautelet 28/03/2019: use TFILE instead of unit number for set_iluout_timing
!  P. Wautelet 19/04/2019: removed unused dummy arguments and variables
!  P. Wautelet 26/04/2019: replace non-standard FLOAT function by REAL function
!  P. Wautelet 20/05/2019: add name argument to ADDnFIELD_ll + new ADD4DFIELD_ll subroutine
!  J. Escobar  09/07/2019: norme Doctor -> Rename Module Type variable TZ -> T
!  J. Escobar  09/07/2019: for bug in management of XLSZWSM variable, add/use specific 2D TLSFIELD2D_ll pointer
!  P. Wautelet 23/07/2019: OpenACC: move data creations from resolved_cloud to modeln and optimize updates
!  P. Wautelet 13/09/2019: budget: simplify and modernize date/time management
!  J. Escobar  27/09/2019: add missing report timing of RESOLVED_ELEC
!  P. Wautelet 02-03/2020: use the new data structures and subroutines for budgets
!  P. Wautelet 12/10/2020: Write_les_n: remove HLES_AVG dummy argument and group all 4 calls
!  F. Auguste  01/02/2021: add IBM
!  T. Nagel    01/02/2021: add turbulence recycling
!  P. Wautelet 19/02/2021: add NEGA2 term for SV budgets
!  J.L. Redelsperger 03/2021: add Call NHOA_COUPLN (coupling O & A LES version)
!  R. Schoetter   12/2021:  multi-level coupling between MesoNH and SURFEX
!  A. Costes      12/2021: add Blaze fire model
!  C. Barthe   07/04/2022: deallocation of ZSEA
!  P. Wautelet 08/12/2022: bugfix if no TDADFILE
!  P. Wautelet 13/01/2023: manage close of backup files outside of MODEL_n
!                          (useful to close them in reverse model order (child before parent, needed by WRITE_BALLOON_n)
!  J. Wurtz       01/2023: correction for mean in SURFEX outputs
!  C. Barthe   03/02/2022: cloud electrification is now called from resolved_cloud to avoid duplicated routines
!  P. Wautelet 02/02/2024: restructure backups/outputs lists
!  A. Marcel      01/2025: relaxation of the small fraction assumption
!!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_2D_FRC
USE MODD_ADV_n
USE MODD_AIRCRAFT_BALLOON
USE MODD_ARGSLIST_ll,     ONLY : LIST_ll
USE MODD_BIKHARDT_n
USE MODD_BLANK_n
USE MODD_BLOWSNOW
USE MODD_BLOWSNOW_n
use modd_budget,          only: cbutype, lbu_ru, lbu_rv, lbu_rw, lbudget_u, lbudget_v, lbudget_w, lbudget_sv, lbu_enable, &
                                NBUDGET_U, NBUDGET_V, NBUDGET_W, NBUDGET_SV1, nbumod, nbutime,                            &
                                tbudgets, tbudgets_ptr, tbuconf, tburhodj,                                                                       &
                                xtime_bu, xtime_bu_process
USE MODD_CH_AERO_n,      ONLY: XSOLORG, XMI
USE MODD_CH_MNHC_n,      ONLY: LUSECHEM,LCH_CONV_LINOX,LUSECHAQ,LUSECHIC, &
                               LCH_INIT_FIELD
USE MODD_CLOUD_MF_n
USE MODD_CLOUDPAR_n
USE MODD_CONF
USE MODD_CONF_n
USE MODD_CST,            ONLY: CST
USE MODD_CURVCOR_n
USE MODD_DEEP_CONVECTION_n
USE MODD_DIM_n
USE MODD_DIMPHYEX,       ONLY: DIMPHYEX_t
USE MODD_DRAG_n
USE MODD_DUST,           ONLY: LDUST
USE MODD_DYN
USE MODD_DYN_n
USE MODD_DYNZD
USE MODD_DYNZD_n
USE MODD_ELEC_DESCR
USE MODD_EOL_MAIN
USE MODD_FIELD_n
USE MODD_FIRE_n
USE MODD_FRC
USE MODD_FRC_n
USE MODD_GET_n
USE MODD_GRID,           ONLY: XLONORI,XLATORI
USE MODD_GRID_n
USE MODD_IBM_PARAM_n,    ONLY: CIBM_ADV, LIBM, LIBM_TROUBLE, XIBM_LS
USE MODD_ICE_C1R3_DESCR, ONLY: XRTMIN_C1R3=>XRTMIN
USE MODD_IO,             ONLY: LIO_NO_WRITE, TFILEDATA, TFILE_SURFEX, TFILE_DUMMY
USE MODD_LBC_n
USE MODD_LES
USE MODD_LES_BUDGET
USE MODD_LIMA_PRECIP_SCAVENGING_n
USE MODD_LSFIELD_n
USE MODD_LUNIT,          ONLY: TOUTDATAFILE
USE MODD_LUNIT_n,        ONLY: TDIAFILE,TINIFILE,TINIFILEPGD,TLUOUT
USE MODD_MEAN_FIELD
USE MODD_MEAN_FIELD_n
USE MODD_METRICS_n
USE MODD_MNH_SURFEX_n
#ifndef MNH_COMPILER_CCE
USE MODD_NEB_n,          ONLY: LSIGMAS, LSUBG_COND, VSIGQSAT
#else
USE MODD_NEB_n,          ONLY: LSIGMAS, LSUBG_COND, VSIGQSAT_MODD => VSIGQSAT
#endif
USE MODD_NESTING
USE MODD_NSV
USE MODD_NUDGING_n
USE MODD_OUT_n
USE MODD_PARAM_C1R3,     ONLY: NSEDI => LSEDI, NHHONI => LHHONI
USE MODD_PARAM_C2R2,     ONLY: NSEDC => LSEDC, NRAIN => LRAIN, NACTIT => LACTIT,LACTTKE,LDEPOC
USE MODD_PARAMETERS
USE MODD_PARAM_ICE_n,    ONLY: LWARM,LSEDIC,LCONVHG,LDEPOSC, CSUBG_AUCV_RC
USE MODD_PARAM_LIMA,     ONLY: MSEDC => LSEDC, NMOM_C, NMOM_R,                           &
                               MACTIT => LACTIT, LSCAV, NMOM_I,                          &
                               MSEDI => LSEDI, MHHONI => LHHONI, NMOM_S, NMOM_G, NMOM_H, &
                               XRTMIN_LIMA=>XRTMIN, MACTTKE=>LACTTKE, LPTSPLIT
USE MODD_PARAM_MFSHALL_n
USE MODD_PARAM_n
USE MODD_PAST_FIELD_n
USE MODD_PRECIP_n
use modd_precision,      only: MNHTIME
USE MODD_PROFILER_n
USE MODD_RADIATIONS_n,   ONLY: XTSRAD,XSCAFLASWD,XDIRFLASWD,XDIRSRFSWD, XAER, XDTHRAD
USE MODD_RAIN_ICE_DESCR_n, ONLY: XRTMIN
USE MODD_RECYCL_PARAM_n, ONLY: LRECYCL
USE MODD_REF,            ONLY: LCOUPLES
USE MODD_REF_n
USE MODD_SALT,           ONLY: LSALT
USE MODD_SERIES,         ONLY: LSERIES
USE MODD_SERIES_n,       ONLY: NFREQSERIES
USE MODD_STATION_n
USE MODD_SUB_MODEL_n
USE MODD_TIME
USE MODD_TIME_n
USE MODD_TIMEZ
USE MODD_TURB_n
USE MODD_TYPE_DATE,      ONLY: DATE_TIME
USE MODD_VISCOSITY
!
USE MODE_AIRCRAFT_BALLOON
use mode_budget,           only: Budget_store_init, Budget_store_end
USE MODE_DATETIME
USE MODE_ELEC_ll
USE MODE_FILL_DIMPHYEX,    ONLY: FILL_DIMPHYEX
USE MODE_GRIDCART
USE MODE_GRIDPROJ
USE MODE_IO_FIELD_WRITE,   only: IO_Field_user_write, IO_Fieldlist_write, IO_Header_write
USE MODE_IO_FILE,          only: IO_File_close, IO_File_open
USE MODE_IO_MANAGE_STRUCT, only: IO_File_add2list, IO_File_remove_from_list, &
                                 IO_Is_backup_time, IO_Is_output_time, IO_BakOut_file_create
USE MODE_ll
USE MODE_MNH_TIMING
#ifdef MNH_OPENACC
USE MODE_MNH_ZWORK,      ONLY: MNH_MEM_GET, MNH_MEM_POSITION_PIN, MNH_MEM_RELEASE
#endif
USE MODE_MODELN_HANDLER
USE MODE_MPPDB
USE MODE_MSG
USE MODE_ONE_WAY_n
USE MODE_WRITE_AIRCRAFT_BALLOON
use mode_write_les_n,               only: Write_les_n
use mode_write_lfifmn_fordiachro_n, only: WRITE_LFIFMN_FORDIACHRO_n
USE MODE_WRITE_STATPROF_n,          ONLY: WRITE_STATPROF_n
!
USE MODI_ADDFLUCTUATIONS
USE MODI_ADVECTION_METSV
USE MODI_ADVECTION_UVW
USE MODI_ADVECTION_UVW_CEN
USE MODI_ADV_FORCING_n
USE MODI_AER_MONITOR_n
USE MODI_BLOWSNOW
USE MODI_BOUNDARIES
USE MODI_BUDGET_FLAGS
USE MODI_CART_COMPRESS
USE MODI_CH_MONITOR_n
USE MODI_DIAG_SURF_ATM_N
USE MODI_DYN_SOURCES
USE MODI_END_DIAG_IN_RUN
USE MODI_ENDSTEP
USE MODI_ENDSTEP_BUDGET
USE MODI_EXCHANGE
USE MODI_FORCING
USE MODI_FORC_SQUALL_LINE
USE MODI_FORC_WIND
USE MODI_GET_HALO
USE MODI_GRAVITY_IMPL
USE MODI_IBM_INIT
USE MODI_IBM_FORCING
USE MODI_IBM_FORCING_TR
USE MODI_IBM_FORCING_ADV
USE MODI_INI_DIAG_IN_RUN
USE MODI_INI_LG
USE MODI_INI_MEAN_FIELD
USE MODI_INITIAL_GUESS
USE MODI_LES_INI_TIMESTEP_n
USE MODI_LES_N
USE MODI_LIMA_PRECIP_SCAVENGING
USE MODI_LS_COUPLING
USE MODI_MASK_COMPRESS
USE MODI_MEAN_FIELD
USE MODI_MNHGET_SURF_PARAM_n
USE MODI_MNHWRITE_ZS_DUMMY_n
USE MODI_NUDGING
USE MODI_NUM_DIFF
USE MODI_PHYS_PARAM_n
USE MODI_PRESSUREZ
USE MODI_PROFILER_n
USE MODI_RAD_BOUND
USE MODI_RECYCLING
USE MODI_RELAX2FW_ION 
USE MODI_RELAXATION
USE MODI_REL_FORCING_n
USE MODI_RESOLVED_CLOUD
USE MODI_RESOLVED_ELEC_n
USE MODI_SERIES_N
USE MODI_SETLB_LG
USE MODI_SET_MASK
USE MODI_SHUMAN
USE MODI_SPAWN_LS_n
USE MODI_STATION_n
USE MODI_TURB_CLOUD_INDEX
USE MODI_TWO_WAY
USE MODI_UPDATE_NSV
USE MODI_VISCOSITY
USE MODI_WRITE_DESFM_n
USE MODI_WRITE_DIAG_SURF_ATM_N
USE MODI_WRITE_LFIFM_n
USE MODI_WRITE_SERIES_n
USE MODI_WRITE_SURF_ATM_N
#ifdef MNH_BITREP_OMP
USE MODI_BITREPZ
#endif
!
IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
!
!
INTEGER,                  INTENT(IN)    :: KTCOUNT    ! Temporal loop index of model KMODEL
TYPE(TFILEDATA), POINTER, INTENT(OUT)   :: TPBAKFILE  ! Pointer for backup file
TYPE(DATE_TIME),          INTENT(OUT)   :: TPDTMODELN ! Time of current model computation
LOGICAL,                  INTENT(INOUT) :: OEXIT      ! Switch for the end of the temporal loop
!
!*       0.2   declarations of local variables
!
LOGICAL, PARAMETER :: GFULLSTAT_PRESS_SLV = .FALSE. ! To write detailed stats of pressure solver in files
!
INTEGER :: ILUOUT      ! Logical unit number for the output listing
INTEGER :: IIU,IJU,IKU ! array size in first, second and third dimensions
INTEGER :: IIB,IIE,IJB,IJE ! index values for the physical subdomain
INTEGER :: JSV,JRR     ! Loop index for scalar and moist variables
INTEGER  :: INBVAR              ! number of HALO2_lls to allocate
INTEGER  :: IINFO_ll            ! return code of parallel routine
INTEGER :: IVERB                ! LFI verbosity level
LOGICAL :: GSTEADY_DMASS        ! conditional call to mass computation
!
                                ! for computing time analysis
REAL(kind=MNHTIME), DIMENSION(2) :: ZTIME, ZTIME1, ZTIME2, ZEND, ZTOT, ZALL, ZTOT_PT, ZBLAZETOT
REAL(kind=MNHTIME), DIMENSION(2) :: ZTIME_STEP,ZTIME_STEP_PTS
CHARACTER                 :: YMI
INTEGER                   :: IPOINTS
CHARACTER(len=16)         :: YTCOUNT,YPOINTS
CHARACTER(LEN=:), ALLOCATABLE :: YDADNAME
!
INTEGER :: IBAKSTEP ! Number of timesteps since previous backup
INTEGER :: ISYNCHRO          ! model synchronic index relative to its father
                             ! = 1  for the first time step in phase with DAD
                             ! = 0  for the last  time step (out of phase)
INTEGER      :: IMI ! Current model index
REAL, DIMENSION(:,:), POINTER, CONTIGUOUS :: ZSEA
REAL, DIMENSION(:,:), POINTER, CONTIGUOUS :: ZTOWN
INTEGER      :: INUMBAK ! Current backup number
INTEGER      :: INUMOUT ! Current output number
! Dummy pointers needed to correct an ifort Bug
REAL, DIMENSION(:), POINTER :: DPTR_XZHAT
REAL, DIMENSION(:), POINTER :: DPTR_XBMX1,DPTR_XBMX2,DPTR_XBMX3,DPTR_XBMX4
REAL, DIMENSION(:), POINTER :: DPTR_XBMY1,DPTR_XBMY2,DPTR_XBMY3,DPTR_XBMY4
REAL, DIMENSION(:), POINTER :: DPTR_XBFX1,DPTR_XBFX2,DPTR_XBFX3,DPTR_XBFX4
REAL, DIMENSION(:), POINTER :: DPTR_XBFY1,DPTR_XBFY2,DPTR_XBFY3,DPTR_XBFY4
CHARACTER(LEN=4), DIMENSION(:), POINTER :: DPTR_CLBCX,DPTR_CLBCY
INTEGER, DIMENSION(:,:,:), POINTER :: DPTR_NKLIN_LBXU,DPTR_NKLIN_LBYU,DPTR_NKLIN_LBXV,DPTR_NKLIN_LBYV
INTEGER, DIMENSION(:,:,:), POINTER :: DPTR_NKLIN_LBXW,DPTR_NKLIN_LBYW,DPTR_NKLIN_LBXM,DPTR_NKLIN_LBYM
REAL, DIMENSION(:,:,:), POINTER :: DPTR_XCOEFLIN_LBXU,DPTR_XCOEFLIN_LBYU
REAL, DIMENSION(:,:,:), POINTER :: DPTR_XCOEFLIN_LBXV,DPTR_XCOEFLIN_LBYV
REAL, DIMENSION(:,:,:), POINTER :: DPTR_XCOEFLIN_LBXW,DPTR_XCOEFLIN_LBYW
REAL, DIMENSION(:,:,:), POINTER :: DPTR_XCOEFLIN_LBXM,DPTR_XCOEFLIN_LBYM
REAL, DIMENSION(:,:,:),   POINTER :: DPTR_XLBXUM,DPTR_XLBYUM,DPTR_XLBXVM,DPTR_XLBYVM
REAL, DIMENSION(:,:,:),   POINTER :: DPTR_XLBXWM,DPTR_XLBYWM,DPTR_XLBXTHM,DPTR_XLBYTHM
REAL, DIMENSION(:,:,:),   POINTER :: DPTR_XLBXTKEM,DPTR_XLBYTKEM
REAL, DIMENSION(:,:,:,:),   POINTER :: DPTR_XLBXSVM,DPTR_XLBYSVM                
REAL, DIMENSION(:,:,:,:), POINTER :: DPTR_XLBXRM,DPTR_XLBYRM
REAL, DIMENSION(:,:,:),   POINTER ::  DPTR_XZZ
REAL, DIMENSION(:,:,:), POINTER :: DPTR_XLSUM,DPTR_XLSVM,DPTR_XLSWM,DPTR_XLSTHM,DPTR_XLSRVM
REAL, DIMENSION(:,:,:), POINTER :: DPTR_XLSUS,DPTR_XLSVS,DPTR_XLSWS,DPTR_XLSTHS,DPTR_XLSRVS
REAL, DIMENSION(:,:),   POINTER :: DPTR_XLSZWSM,DPTR_XLSZWSS
REAL, DIMENSION(:,:,:), POINTER :: DPTR_XLBXUS,DPTR_XLBYUS,DPTR_XLBXVS,DPTR_XLBYVS
REAL, DIMENSION(:,:,:), POINTER :: DPTR_XLBXWS,DPTR_XLBYWS,DPTR_XLBXTHS,DPTR_XLBYTHS
REAL, DIMENSION(:,:,:), POINTER :: DPTR_XLBXTKES,DPTR_XLBYTKES
REAL, DIMENSION(:,:,:,:), POINTER :: DPTR_XLBXRS,DPTR_XLBYRS,DPTR_XLBXSVS,DPTR_XLBYSVS
!
REAL, DIMENSION(:,:,:), POINTER :: DPTR_XRHODJ,DPTR_XUM,DPTR_XVM,DPTR_XWM,DPTR_XTHM
REAL, DIMENSION(:,:,:), POINTER :: DPTR_XTKEM,DPTR_XRUS,DPTR_XRVS,DPTR_XRWS,DPTR_XRTHS
REAL, DIMENSION(:,:,:), POINTER :: DPTR_XRTKES,DPTR_XDIRFLASWD,DPTR_XSCAFLASWD,DPTR_XDIRSRFSWD
REAL, DIMENSION(:,:,:,:), POINTER :: DPTR_XRM,DPTR_XSVM,DPTR_XRRS,DPTR_XRSVS
REAL, DIMENSION(:,:), POINTER :: DPTR_XINPRC,DPTR_XINPRR,DPTR_XINPRS,DPTR_XINPRG
REAL, DIMENSION(:,:), POINTER :: DPTR_XINPRH,DPTR_XPRCONV,DPTR_XPRSCONV
LOGICAL, DIMENSION(:,:),POINTER :: DPTR_GMASKkids
!
REAL, DIMENSION(SIZE(XTHT,1),SIZE(XTHT,2),SIZE(XTHT,3))           :: ZSPEEDC
REAL, DIMENSION(SIZE(XTHT,1),SIZE(XTHT,2),SIZE(XTHT,3))           :: ZSPEEDR
REAL, DIMENSION(SIZE(XTHT,1),SIZE(XTHT,2),SIZE(XTHT,3))           :: ZSPEEDS
REAL, DIMENSION(SIZE(XTHT,1),SIZE(XTHT,2),SIZE(XTHT,3))           :: ZSPEEDG
REAL, DIMENSION(SIZE(XTHT,1),SIZE(XTHT,2),SIZE(XTHT,3))           :: ZSPEEDH
REAL, DIMENSION(SIZE(XTHT,1),SIZE(XTHT,2),SIZE(XTHT,3))           :: ZINPRC3D
REAL, DIMENSION(SIZE(XTHT,1),SIZE(XTHT,2),SIZE(XTHT,3))           :: ZINPRS3D
REAL, DIMENSION(SIZE(XTHT,1),SIZE(XTHT,2),SIZE(XTHT,3))           :: ZINPRG3D
REAL, DIMENSION(SIZE(XTHT,1),SIZE(XTHT,2),SIZE(XTHT,3))           :: ZINPRH3D
!
LOGICAL :: KWARM
LOGICAL :: KRAIN
LOGICAL :: KSEDC
LOGICAL :: KACTIT
LOGICAL :: KSEDI
LOGICAL :: KHHONI
!
#ifndef MNH_OPENACC
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZRUS,ZRVS
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZRWS
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZPABST !To give pressure at t
                                                     ! (and not t+1) to resolved_cloud
#else
REAL, DIMENSION(:,:,:), POINTER, CONTIGUOUS :: ZRUS,ZRVS
REAL, DIMENSION(:,:,:), POINTER, CONTIGUOUS :: ZRWS
REAL, DIMENSION(:,:,:), POINTER, CONTIGUOUS :: ZPABST !To give pressure at t
                                                     ! (and not t+1) to resolved_cloud
#endif
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZJ
!
TYPE(LIST_ll), POINTER :: TZFIELDC_ll   ! list of fields to exchange
TYPE(HALO2LIST_ll), POINTER :: TZHALO2C_ll   ! list of fields to exchange
#ifdef MNH_OPENACC
TYPE(HALO2LIST_ll), SAVE , POINTER :: TZHALO2_UT,TZHALO2_VT,TZHALO2_WT
LOGICAL , SAVE :: GFIRST_CALL_MODELN = .TRUE.
#endif
LOGICAL :: GCLD                     ! conditionnal call for dust wet deposition
LOGICAL :: GCLOUD_ONLY              ! conditionnal radiation computations for
                                !      the only cloudy columns
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZWETDEPAER

LOGICAL :: GNVS_NO_0
!
TYPE(TFILEDATA),POINTER :: TZOUTFILE
! TYPE(TFILEDATA),SAVE    :: TZDIACFILE
TYPE(DIMPHYEX_t) :: YLDIMPHYEX
!-------------------------------------------------------------------------------
!
!
#ifdef MNH_COMPILER_CCE
!Bypass cray bug with scalar pointer
REAL :: VSIGQSAT
VSIGQSAT = VSIGQSAT_MODD
#endif

!-------------------------------------------------------------------------------

TPBAKFILE=> NULL()
TZOUTFILE=> NULL()

#ifndef MNH_OPENACC
allocate( ZRUS  (SIZE(XTHT,1),SIZE(XTHT,2),SIZE(XTHT,3)) )
allocate( ZRVS  (SIZE(XTHT,1),SIZE(XTHT,2),SIZE(XTHT,3)) )
allocate( ZRWS  (SIZE(XTHT,1),SIZE(XTHT,2),SIZE(XTHT,3)) )
allocate( ZPABST(SIZE(XTHT,1),SIZE(XTHT,2),SIZE(XTHT,3)) )
#else
!Pin positions in the pools of MNH memory
CALL MNH_MEM_POSITION_PIN( 'MODEL_n 1')
CALL MNH_MEM_GET( ZRUS,   SIZE( XTHT, 1 ), SIZE( XTHT, 2 ), SIZE( XTHT, 3 ) )
CALL MNH_MEM_GET( ZRVS,   SIZE( XTHT, 1 ), SIZE( XTHT, 2 ), SIZE( XTHT, 3 ) )
CALL MNH_MEM_GET( ZRWS,   SIZE( XTHT, 1 ), SIZE( XTHT, 2 ), SIZE( XTHT, 3 ) )
CALL MNH_MEM_GET( ZPABST, SIZE( XTHT, 1 ), SIZE( XTHT, 2 ), SIZE( XTHT, 3 ) )
#endif
allocate( ZJ    (SIZE(XTHT,1),SIZE(XTHT,2),SIZE(XTHT,3)) )
allocate( ZWETDEPAER(SIZE(XRSVS,1), SIZE(XRSVS,2), SIZE(XRSVS,3), NSV_AER) )
!
TPDTMODELN = TDTCUR
!
!*       0.    MICROPHYSICAL SCHEME
!              ------------------- 
SELECT CASE(CCLOUD)
CASE('C2R2','KHKO','C3R5')
  KWARM  = .TRUE.        
  KRAIN  = NRAIN
  KSEDC  = NSEDC
  KACTIT = NACTIT
!
  KSEDI  = NSEDI
  KHHONI = NHHONI
CASE('LIMA')
  KRAIN  = NMOM_R.GE.1
  KWARM  = NMOM_C.GE.1
  KSEDC  = MSEDC
  KACTIT = MACTIT
!
  KSEDI  = MSEDI
  KHHONI = MHHONI
CASE('ICE3','ICE4') !default values
  KWARM  = LWARM        
  KRAIN  = .TRUE.
  KSEDC  = .TRUE.
  KACTIT = .FALSE.
!
  KSEDI  = .TRUE.
  KHHONI = .FALSE.
CASE DEFAULT
  ! Set default values to prevent false positives with debugging tools
  ! (ie KSEDC is used in an IF statement and could be correctly detected as not-initialised) (valgrind)
  KWARM  = .FALSE.
  KRAIN  = .FALSE.
  KSEDC  = .FALSE.
  KACTIT = .FALSE.
END SELECT
!
!
!*        1    PRELIMINARY
!              ------------
IMI = GET_CURRENT_MODEL_INDEX()
!
!*       1.0   update NSV_* variables for current model
!              ----------------------------------------
!
CALL UPDATE_NSV(IMI)
!
GNVS_NO_0 = ( NSV /= 0 )
!
!*       1.1   RECOVER THE LOGICAL UNIT NUMBER FOR THE OUTPUT PRINTS
!
ILUOUT = TLUOUT%NLU
!
!*       1.2   SET ARRAY SIZE
!
CALL GET_DIM_EXT_ll('B',IIU,IJU)
IKU=NKMAX+2*JPVEXT
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
!
CALL FILL_DIMPHYEX(YLDIMPHYEX, SIZE(XWT,1), SIZE(XWT,2), SIZE(XWT,3))
!
IF (IMI==1) THEN
  GSTEADY_DMASS=LSTEADYLS
ELSE
  GSTEADY_DMASS=.FALSE.
END IF
!
!*       1.3   OPEN THE DIACHRONIC FILE
!
IF (KTCOUNT == 1) THEN
!
  NULLIFY(TFIELDS_ll,TLSFIELD_ll,TFIELDT_ll)
  NULLIFY(TLSFIELD2D_ll)
  NULLIFY(THALO2T_ll)
  NULLIFY(TLSHALO2_ll)
  NULLIFY(TFIELDSC_ll)
!
  ALLOCATE(XWT_ACT_NUC(SIZE(XWT,1),SIZE(XWT,2),SIZE(XWT,3)))
  ALLOCATE(GMASKkids(SIZE(XWT,1),SIZE(XWT,2)))
!
  IF ( .NOT. LIO_NO_WRITE ) THEN
    CALL IO_File_open(TDIAFILE)
!
    CALL IO_Header_write(TDIAFILE)
    CALL WRITE_DESFM_n(IMI,TDIAFILE)
    CALL WRITE_LFIFMN_FORDIACHRO_n(TDIAFILE)
  END IF
!
!*       1.4   Initialization of the list of fields for the halo updates
!
!                 a) Sources terms
!
  CALL ADD3DFIELD_ll( TFIELDS_ll, XRUS,      'MODEL_n::XRUS' )
  CALL ADD3DFIELD_ll( TFIELDS_ll, XRVS,      'MODEL_n::XRVS' )
  CALL ADD3DFIELD_ll( TFIELDS_ll, XRWS,      'MODEL_n::XRWS' )
  CALL ADD3DFIELD_ll( TFIELDS_ll, XRTHS,     'MODEL_n::XRTHS' )
  CALL ADD3DFIELD_ll( TFIELDS_ll, XRUS_PRES, 'MODEL_n::XRUS_PRES' )
  CALL ADD3DFIELD_ll( TFIELDS_ll, XRVS_PRES, 'MODEL_n::XRVS_PRES' )
  CALL ADD3DFIELD_ll( TFIELDS_ll, XRWS_PRES, 'MODEL_n::XRWS_PRES' )
  CALL ADD3DFIELD_ll( TFIELDS_ll, XRTHS_CLD, 'MODEL_n::XRTHS_CLD' )
  IF (SIZE(XRTKES,1) /= 0) CALL ADD3DFIELD_ll( TFIELDS_ll, XRTKES, 'MODEL_n::XRTKES' )
  CALL ADD4DFIELD_ll( TFIELDS_ll, XRRS     (:,:,:,1:NRR), 'MODEL_n::XRRS' )
  CALL ADD4DFIELD_ll( TFIELDS_ll, XRRS_CLD (:,:,:,1:NRR), 'MODEL_n::XRRS_CLD' )
  IF ( GNVS_NO_0 ) CALL ADD4DFIELD_ll( TFIELDS_ll, XRSVS    (:,:,:,1:NSV), 'MODEL_n::XRSVS')
  IF ( GNVS_NO_0 ) CALL ADD4DFIELD_ll( TFIELDS_ll, XRSVS_CLD(:,:,:,1:NSV), 'MODEL_n::XRSVS_CLD')
  IF (SIZE(XSRCT,1) /= 0) CALL ADD3DFIELD_ll( TFIELDS_ll, XSRCT, 'MODEL_n::XSRCT' )
  ! Fire model parallel setup
  IF (LBLAZE) THEN
    CALL ADD3DFIELD_ll( TFIELDS_ll, XLSPHI,     'MODEL_n::XLSPHI')
    CALL ADD3DFIELD_ll( TFIELDS_ll, XBMAP,      'MODEL_n::XBMAP')
    CALL ADD3DFIELD_ll( TFIELDS_ll, XFMRFA,     'MODEL_n::XFMRFA')
    CALL ADD3DFIELD_ll( TFIELDS_ll, XFMWF0,     'MODEL_n::XFMWF0')
    CALL ADD3DFIELD_ll( TFIELDS_ll, XFMR0,      'MODEL_n::XFMR0')
    CALL ADD3DFIELD_ll( TFIELDS_ll, XFMR00,     'MODEL_n::XFMR00')
    CALL ADD3DFIELD_ll( TFIELDS_ll, XFMIGNITION, 'MODEL_n::XFMIGNITION')
    CALL ADD3DFIELD_ll( TFIELDS_ll, XFMFUELTYPE, 'MODEL_n::XFMFUELTYPE')
    CALL ADD3DFIELD_ll( TFIELDS_ll, XFIRETAU,   'MODEL_n::XFIRETAU')
    CALL ADD4DFIELD_ll( TFIELDS_ll, XFLUXPARAMH(:,:,:,1:SIZE(XFLUXPARAMH,4)), 'MODEL_n::XFLUXPARAMH')
    CALL ADD4DFIELD_ll( TFIELDS_ll, XFLUXPARAMW(:,:,:,1:SIZE(XFLUXPARAMW,4)), 'MODEL_n::XFLUXPARAMW')
    CALL ADD3DFIELD_ll( TFIELDS_ll, XFIRERW,    'MODEL_n::XFIRERW')
    CALL ADD3DFIELD_ll( TFIELDS_ll, XFMASE,     'MODEL_n::XFMASE')
    CALL ADD3DFIELD_ll( TFIELDS_ll, XFMAWC,     'MODEL_n::XFMAWC')
    CALL ADD3DFIELD_ll( TFIELDS_ll, XFMWALKIG,  'MODEL_n::XFMWALKIG')
    CALL ADD3DFIELD_ll( TFIELDS_ll, XFMFLUXHDH, 'MODEL_n::XFMFLUXHDH')
    CALL ADD3DFIELD_ll( TFIELDS_ll, XFMFLUXHDW, 'MODEL_n::XFMFLUXHDW')
    CALL ADD3DFIELD_ll( TFIELDS_ll, XFMHWS,     'MODEL_n::XFMHWS')
    CALL ADD3DFIELD_ll( TFIELDS_ll, XFMWINDU,   'MODEL_n::XFMWINDU')
    CALL ADD3DFIELD_ll( TFIELDS_ll, XFMWINDV,   'MODEL_n::XFMWINDV')
    CALL ADD3DFIELD_ll( TFIELDS_ll, XFMWINDW,   'MODEL_n::XFMWINDW')
    CALL ADD3DFIELD_ll( TFIELDS_ll, XFMGRADOROX, 'MODEL_n::XFMGRADOROX')
    CALL ADD3DFIELD_ll( TFIELDS_ll, XFMGRADOROY, 'MODEL_n::XFMGRADOROY')
  END IF
  !
  IF ((LNUMDIFU .OR. LNUMDIFTH .OR. LNUMDIFSV) ) THEN
  !
  !                 b) LS fields
  !
    CALL ADD3DFIELD_ll( TLSFIELD_ll,   XLSUM,   'MODEL_n::XLSUM'   )
    CALL ADD3DFIELD_ll( TLSFIELD_ll,   XLSVM,   'MODEL_n::XLSVM'   )
    CALL ADD3DFIELD_ll( TLSFIELD_ll,   XLSWM,   'MODEL_n::XLSWM'   )
    CALL ADD3DFIELD_ll( TLSFIELD_ll,   XLSTHM,  'MODEL_n::XLSTHM'  )
    CALL ADD2DFIELD_ll( TLSFIELD2D_ll, XLSZWSM, 'MODEL_n::XLSZWSM' )
    IF (NRR >= 1) THEN
      CALL ADD3DFIELD_ll( TLSFIELD_ll, XLSRVM, 'MODEL_n::XLSRVM' )
    ENDIF
  !
  !                 c) Fields at t
  !
    CALL ADD3DFIELD_ll( TFIELDT_ll, XUT,  'MODEL_n::XUT'  )
    CALL ADD3DFIELD_ll( TFIELDT_ll, XVT,  'MODEL_n::XVT'  )
    CALL ADD3DFIELD_ll( TFIELDT_ll, XWT,  'MODEL_n::XWT'  )
    CALL ADD3DFIELD_ll( TFIELDT_ll, XTHT, 'MODEL_n::XTHT' )
    IF (SIZE(XRTKES,1) /= 0) CALL ADD3DFIELD_ll( TFIELDT_ll, XTKET, 'MODEL_n::XTKET' )
    CALL ADD4DFIELD_ll(TFIELDT_ll, XRT (:,:,:,1:NRR), 'MODEL_n::XSV'  )
    CALL ADD4DFIELD_ll(TFIELDT_ll, XSVT(:,:,:,1:NSV), 'MODEL_n::XSVT' )
  !
  !*       1.5   Initialize the list of fields for the halo updates (2nd layer)
  !
    INBVAR = 4+NRR+NSV
    IF (SIZE(XRTKES,1) /= 0) INBVAR=INBVAR+1
    CALL INIT_HALO2_ll (THALO2T_ll,  INBVAR,       IIU, IJU, IKU, 'MODEL_n::HALO2T'  )
    CALL INIT_HALO2_ll (TLSHALO2_ll, 4+MIN(1,NRR), IIU, IJU, IKU, 'MODEL_n::LSHALO2' )
  !
  !*       1.6   Initialise the 2nd layer of the halo of the LS fields
  !
    IF ( LSTEADYLS ) THEN
       CALL UPDATE_HALO_ll(TLSFIELD_ll, IINFO_ll)
       CALL UPDATE_HALO_ll(TLSFIELD2D_ll,IINFO_ll) 
       CALL UPDATE_HALO2_ll(TLSFIELD_ll, TLSHALO2_ll, IINFO_ll)
    END IF
  END IF
  !
!
  !
  XT_START     = 0.0_MNHTIME
  !
  XT_STORE     = 0.0_MNHTIME
  XT_BOUND     = 0.0_MNHTIME
  XT_GUESS     = 0.0_MNHTIME
  XT_FORCING   = 0.0_MNHTIME
  XT_NUDGING   = 0.0_MNHTIME
  XT_ADV       = 0.0_MNHTIME
  XT_ADVUVW    = 0.0_MNHTIME
  XT_GRAV      = 0.0_MNHTIME
  XT_SOURCES   = 0.0_MNHTIME
  !
  XT_DIFF      = 0.0_MNHTIME
  XT_RELAX     = 0.0_MNHTIME
  XT_PARAM     = 0.0_MNHTIME
  XT_SPECTRA   = 0.0_MNHTIME
  XT_HALO      = 0.0_MNHTIME
  XT_VISC      = 0.0_MNHTIME
  XT_RAD_BOUND = 0.0_MNHTIME
  XT_PRESS     = 0.0_MNHTIME
  !
  XT_CLOUD     = 0.0_MNHTIME
  XT_STEP_SWA  = 0.0_MNHTIME
  XT_STEP_MISC = 0.0_MNHTIME
  XT_COUPL     = 0.0_MNHTIME
  XT_1WAY      = 0.0_MNHTIME
  XT_STEP_BUD  = 0.0_MNHTIME
  !
  XT_RAD       = 0.0_MNHTIME
  XT_DCONV     = 0.0_MNHTIME
  XT_GROUND    = 0.0_MNHTIME
  XT_TURB      = 0.0_MNHTIME
  XT_MAFL      = 0.0_MNHTIME
  XT_DRAG      = 0.0_MNHTIME
  XT_EOL       = 0.0_MNHTIME
  XT_TRACER    = 0.0_MNHTIME
  XT_SHADOWS   = 0.0_MNHTIME
  XT_ELEC      = 0.0_MNHTIME
  XT_CHEM      = 0.0_MNHTIME
  XT_2WAY      = 0.0_MNHTIME
  !
  XT_IBM_FORC  = 0.0_MNHTIME
  ! Blaze fire model
  XFIREPERF    = 0.0_MNHTIME
  !
END IF
!
!*       1.7   Allocation of arrays for observation diagnostics
!
CALL INI_DIAG_IN_RUN(IIU,IJU,IKU,LFLYER,LSTATION,LPROFILER)
!
!
CALL SECOND_MNH2(ZEND)
!
!-------------------------------------------------------------------------------
!
!*       2.    ONE-WAY NESTING AND LARGE SCALE FIELD REFRESH
!              ---------------------------------------------
!
!
CALL SECOND_MNH2(ZTIME1)
!
ISYNCHRO = MODULO (KTCOUNT, NDTRATIO(IMI) )      ! test of synchronisation
!
!
IF (LCOUPLES.AND.LOCEAN) THEN
   CALL NHOA_COUPL_n(NDAD(IMI),XTSTEP,IMI,KTCOUNT,IKU)
END IF
! No Gridnest in coupled OA LES for now
IF (.NOT. LCOUPLES .AND. IMI/=1 .AND. NDAD(IMI)/=IMI .AND. (ISYNCHRO==1 .OR. NDTRATIO(IMI) == 1) ) THEN     
!                                                                         
! Use dummy pointers to correct an ifort BUG
  DPTR_XBMX1=>XBMX1
  DPTR_XBMX2=>XBMX2
  DPTR_XBMX3=>XBMX3
  DPTR_XBMX4=>XBMX4
  DPTR_XBMY1=>XBMY1
  DPTR_XBMY2=>XBMY2
  DPTR_XBMY3=>XBMY3
  DPTR_XBMY4=>XBMY4
  DPTR_XBFX1=>XBFX1
  DPTR_XBFX2=>XBFX2
  DPTR_XBFX3=>XBFX3
  DPTR_XBFX4=>XBFX4
  DPTR_XBFY1=>XBFY1
  DPTR_XBFY2=>XBFY2
  DPTR_XBFY3=>XBFY3
  DPTR_XBFY4=>XBFY4
  DPTR_CLBCX=>CLBCX
  DPTR_CLBCY=>CLBCY
  !
  DPTR_XZZ=>XZZ
  DPTR_XZHAT=>XZHAT
  DPTR_XCOEFLIN_LBXM=>XCOEFLIN_LBXM
  DPTR_XLSTHM=>XLSTHM
  DPTR_XLSRVM=>XLSRVM
  DPTR_XLSUM=>XLSUM
  DPTR_XLSVM=>XLSVM
  DPTR_XLSWM=>XLSWM
  DPTR_XLSZWSM=>XLSZWSM
  DPTR_XLSTHS=>XLSTHS
  DPTR_XLSRVS=>XLSRVS
  DPTR_XLSUS=>XLSUS
  DPTR_XLSVS=>XLSVS
  DPTR_XLSWS=>XLSWS
  DPTR_XLSZWSS=>XLSZWSS
  !
  IF ( LSTEADYLS                     ) THEN
    NCPL_CUR=0
  ELSE
    IF (NCPL_CUR/=1) THEN
      IF ( KTCOUNT+1 == NCPL_TIMES(NCPL_CUR-1,IMI)  ) THEN
        !
        !  LS sources are interpolated from the LS field 
        ! values of model DAD(IMI)
        CALL SPAWN_LS_n(NDAD(IMI),XTSTEP,IMI,                        &
             DPTR_XBMX1,DPTR_XBMX2,DPTR_XBMX3,DPTR_XBMX4,DPTR_XBMY1,DPTR_XBMY2,DPTR_XBMY3,DPTR_XBMY4,        &
             DPTR_XBFX1,DPTR_XBFX2,DPTR_XBFX3,DPTR_XBFX4,DPTR_XBFY1,DPTR_XBFY2,DPTR_XBFY3,DPTR_XBFY4,        &
             NDXRATIO_ALL(IMI),NDYRATIO_ALL(IMI),                    &
             DPTR_CLBCX,DPTR_CLBCY,DPTR_XZZ,DPTR_XZHAT,LSLEVE,XLEN1,XLEN2,DPTR_XCOEFLIN_LBXM, &
             DPTR_XLSTHM,DPTR_XLSRVM,DPTR_XLSUM,DPTR_XLSVM,DPTR_XLSWM,DPTR_XLSZWSM,                        &
             DPTR_XLSTHS,DPTR_XLSRVS,DPTR_XLSUS,DPTR_XLSVS,DPTR_XLSWS, DPTR_XLSZWSS                         )
      END IF
    END IF
    !
  END IF
  !
  DPTR_NKLIN_LBXU=>NKLIN_LBXU
  DPTR_XCOEFLIN_LBXU=>XCOEFLIN_LBXU
  DPTR_NKLIN_LBYU=>NKLIN_LBYU
  DPTR_XCOEFLIN_LBYU=>XCOEFLIN_LBYU
  DPTR_NKLIN_LBXV=>NKLIN_LBXV
  DPTR_XCOEFLIN_LBXV=>XCOEFLIN_LBXV
  DPTR_NKLIN_LBYV=>NKLIN_LBYV
  DPTR_XCOEFLIN_LBYV=>XCOEFLIN_LBYV
  DPTR_NKLIN_LBXW=>NKLIN_LBXW
  DPTR_XCOEFLIN_LBXW=>XCOEFLIN_LBXW
  DPTR_NKLIN_LBYW=>NKLIN_LBYW
  DPTR_XCOEFLIN_LBYW=>XCOEFLIN_LBYW
  !
  DPTR_NKLIN_LBXM=>NKLIN_LBXM
  DPTR_XCOEFLIN_LBXM=>XCOEFLIN_LBXM
  DPTR_NKLIN_LBYM=>NKLIN_LBYM
  DPTR_XCOEFLIN_LBYM=>XCOEFLIN_LBYM
  !
  DPTR_XLBXUM=>XLBXUM
  DPTR_XLBYUM=>XLBYUM
  DPTR_XLBXVM=>XLBXVM
  DPTR_XLBYVM=>XLBYVM
  DPTR_XLBXWM=>XLBXWM
  DPTR_XLBYWM=>XLBYWM
  DPTR_XLBXTHM=>XLBXTHM
  DPTR_XLBYTHM=>XLBYTHM
  DPTR_XLBXTKEM=>XLBXTKEM
  DPTR_XLBYTKEM=>XLBYTKEM
  DPTR_XLBXRM=>XLBXRM
  DPTR_XLBYRM=>XLBYRM
  DPTR_XLBXSVM=>XLBXSVM
  DPTR_XLBYSVM=>XLBYSVM
  !  
  DPTR_XLBXUS=>XLBXUS
  DPTR_XLBYUS=>XLBYUS
  DPTR_XLBXVS=>XLBXVS
  DPTR_XLBYVS=>XLBYVS
  DPTR_XLBXWS=>XLBXWS
  DPTR_XLBYWS=>XLBYWS
  DPTR_XLBXTHS=>XLBXTHS
  DPTR_XLBYTHS=>XLBYTHS
  DPTR_XLBXTKES=>XLBXTKES
  DPTR_XLBYTKES=>XLBYTKES
  DPTR_XLBXRS=>XLBXRS
  DPTR_XLBYRS=>XLBYRS
  DPTR_XLBXSVS=>XLBXSVS
  DPTR_XLBYSVS=>XLBYSVS
  !
  CALL ONE_WAY_n(NDAD(IMI),XTSTEP,IMI,KTCOUNT,                   &
       DPTR_XBMX1,DPTR_XBMX2,DPTR_XBMX3,DPTR_XBMX4,DPTR_XBMY1,DPTR_XBMY2,DPTR_XBMY3,DPTR_XBMY4,        &
       DPTR_XBFX1,DPTR_XBFX2,DPTR_XBFX3,DPTR_XBFX4,DPTR_XBFY1,DPTR_XBFY2,DPTR_XBFY3,DPTR_XBFY4,        &
       NDXRATIO_ALL(IMI),NDYRATIO_ALL(IMI),NDTRATIO(IMI),         &
       DPTR_CLBCX,DPTR_CLBCY,NRIMX,NRIMY,                                &
       DPTR_NKLIN_LBXU,DPTR_XCOEFLIN_LBXU,DPTR_NKLIN_LBYU,DPTR_XCOEFLIN_LBYU,      &
       DPTR_NKLIN_LBXV,DPTR_XCOEFLIN_LBXV,DPTR_NKLIN_LBYV,DPTR_XCOEFLIN_LBYV,      &
       DPTR_NKLIN_LBXW,DPTR_XCOEFLIN_LBXW,DPTR_NKLIN_LBYW,DPTR_XCOEFLIN_LBYW,      &
       DPTR_NKLIN_LBXM,DPTR_XCOEFLIN_LBXM,DPTR_NKLIN_LBYM,DPTR_XCOEFLIN_LBYM,      &
       GSTEADY_DMASS,CCLOUD,LUSECHAQ,LUSECHIC,                           &
       DPTR_XLBXUM,DPTR_XLBYUM,DPTR_XLBXVM,DPTR_XLBYVM,DPTR_XLBXWM,DPTR_XLBYWM,              &
       DPTR_XLBXTHM,DPTR_XLBYTHM,                                        &
       DPTR_XLBXTKEM,DPTR_XLBYTKEM,                                      &
       DPTR_XLBXRM,DPTR_XLBYRM,DPTR_XLBXSVM,DPTR_XLBYSVM,                          &
       XDRYMASST,XDRYMASSS,                                    &
       DPTR_XLBXUS,DPTR_XLBYUS,DPTR_XLBXVS,DPTR_XLBYVS,DPTR_XLBXWS,DPTR_XLBYWS,              &
       DPTR_XLBXTHS,DPTR_XLBYTHS,                                        &
       DPTR_XLBXTKES,DPTR_XLBYTKES,                                      &
       DPTR_XLBXRS,DPTR_XLBYRS,DPTR_XLBXSVS,DPTR_XLBYSVS                           )
  !
END IF
!
CALL SECOND_MNH2(ZTIME2)                                                  
XT_1WAY = XT_1WAY + ZTIME2 - ZTIME1 
!
!*       2.1    RECYCLING TURBULENCE
!              ---- 
IF (CTURB /= 'NONE' .AND. LRECYCL) THEN
  CALL RECYCLING(XFLUCTUNW,XFLUCTVNN,XFLUCTUTN,XFLUCTVTW,XFLUCTWTW,XFLUCTWTN, &
                 XFLUCTUNE,XFLUCTVNS,XFLUCTUTS,XFLUCTVTE,XFLUCTWTE,XFLUCTWTS, &
                 KTCOUNT)
ENDIF
!
!*       2.2    IBM
!              ----
!
IF (LIBM .AND. KTCOUNT==1) THEN
  !
  IF (.NOT.LCARTESIAN) THEN
    CALL PRINT_MSG(NVERB_FATAL, 'GEN', 'MODELN', 'IBM can only be used in combination with cartesian coordinates')
  ENDIF
  !
  CALL IBM_INIT(XIBM_LS)
  !
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       3.    LATERAL BOUNDARY CONDITIONS EXCEPT FOR NORMAL VELOCITY
!              ------------------------------------------------------
!
ZTIME1=ZTIME2
!
!*       3.1   Set the lagragian variables values at the LB
!
IF( LLG .AND. IMI==1 ) CALL SETLB_LG
!
IF (CCONF == "START" .OR. (CCONF == "RESTA" .AND. KTCOUNT /= 1 )) THEN
CALL BOUNDARIES (                                                   &
            XTSTEP,CLBCX,CLBCY,NRR,NSV,KTCOUNT,                     &
            XLBXUM,XLBXVM,XLBXWM,XLBXTHM,XLBXTKEM,XLBXRM,XLBXSVM,   &
            XLBYUM,XLBYVM,XLBYWM,XLBYTHM,XLBYTKEM,XLBYRM,XLBYSVM,   &
            XLBXUS,XLBXVS,XLBXWS,XLBXTHS,XLBXTKES,XLBXRS,XLBXSVS,   &
            XLBYUS,XLBYVS,XLBYWS,XLBYTHS,XLBYTKES,XLBYRS,XLBYSVS,   &
            XRHODJ,XRHODREF,                                        &
            XUT, XVT, XWT, XTHT, XTKET, XRT, XSVT, XSRCT            )
END IF
!
CALL SECOND_MNH2(ZTIME2)
!
XT_BOUND = XT_BOUND + ZTIME2 - ZTIME1
!
!
! For START/RESTART MPPDB_CHECK use 
!IF ( (IMI==1) .AND. (CCONF == "START") .AND. (KTCOUNT == 2) ) THEN
!   CALL MPPDB_START_DEBUG()
!ENDIF
!IF ( (IMI==1) .AND. (CCONF == "RESTA") .AND. (KTCOUNT == 1) ) THEN
!   CALL MPPDB_START_DEBUG()
!ENDIF
!-------------------------------------------------------------------------------
!* initializes surface number
IF (CSURF=='EXTE') CALL GOTO_SURFEX(IMI)
!-------------------------------------------------------------------------------
!
!*       4.    STORAGE IN A SYNCHRONOUS FILE
!              -----------------------------
!
ZTIME1 = ZTIME2
!
IF ( IO_IS_BACKUP_TIME( IMI, KTCOUNT, INUMBAK ) ) THEN
  CALL IO_BAKOUT_FILE_CREATE( TPBAKFILE, 'MNHBACKUP', IMI, KTCOUNT, INUMBAK )
  !
  IVERB    = TPBAKFILE%NLFIVERB
  !
  CALL IO_File_open(TPBAKFILE)
  !
  CALL WRITE_DESFM_n(IMI,TPBAKFILE)
  CALL IO_Header_write( TPBAKFILE )
  IF ( ASSOCIATED( TPBAKFILE%TDADFILE ) ) THEN
    YDADNAME = TPBAKFILE%TDADFILE%CNAME
  ELSE
    ! Set a dummy name for the dad file. Its non-zero size will allow the writing of some data in the backup file
    YDADNAME = 'DUMMY'
  END IF
  CALL WRITE_LFIFM_n( TPBAKFILE, TRIM( YDADNAME ) )
  TOUTDATAFILE => TPBAKFILE
  CALL MNHWRITE_ZS_DUMMY_n(TPBAKFILE)
  IF (CSURF=='EXTE') THEN
    TFILE_SURFEX => TPBAKFILE
    CALL GOTO_SURFEX(IMI)
    CALL WRITE_SURF_ATM_n(YSURF_CUR,'MESONH','ALL')
    IF ( KTCOUNT > 1) THEN
      CALL DIAG_SURF_ATM_n(YSURF_CUR,'MESONH')
      ! Determine number of timestep since previous backup
      IBAKSTEP = KTCOUNT - NBAK_PREV_STEP
      NBAK_PREV_STEP = KTCOUNT

      CALL WRITE_DIAG_SURF_ATM_n( YSURF_CUR, 'MESONH', 'ALL', IBAKSTEP )
    END IF
    NULLIFY(TFILE_SURFEX)
  END IF
  !
  ! Reinitialize Lagragian variables at every model backup
  IF (LLG .AND. LINIT_LG .AND. CINIT_LG=='FMOUT') THEN
    CALL INI_LG( XXHATM, XYHATM, XZZ, XSVT, XLBXSVM, XLBYSVM )
    IF (IVERB>=5) THEN
      WRITE(UNIT=ILUOUT,FMT=*) '************************************'
      WRITE(UNIT=ILUOUT,FMT=*) '*** Lagrangian variables refreshed after ',TRIM(TPBAKFILE%CNAME),' backup'
      WRITE(UNIT=ILUOUT,FMT=*) '************************************'
    END IF
  END IF
  ! Reinitialise mean variables
  IF (LMEAN_FIELD) THEN
     CALL INI_MEAN_FIELD
  END IF
ELSE
  !Necessary to have a 'valid' CNAME when calling some subroutines
  TPBAKFILE => TFILE_DUMMY
END IF
!
IF ( IO_IS_OUTPUT_TIME( IMI, KTCOUNT, INUMOUT ) ) THEN
  CALL IO_BAKOUT_FILE_CREATE( TZOUTFILE, 'MNHOUTPUT', IMI, KTCOUNT, INUMOUT )
  !
  CALL IO_File_open(TZOUTFILE)
  !
  CALL IO_Header_write(TZOUTFILE)
  CALL IO_Fieldlist_write( TZOUTFILE )
  CALL IO_Field_user_write( TZOUTFILE )
  !
  CALL IO_File_close(TZOUTFILE)
  CALL IO_File_remove_from_list( TZOUTFILE )
END IF
!
CALL SECOND_MNH2(ZTIME2)
!
XT_STORE = XT_STORE + ZTIME2 - ZTIME1
!
!-------------------------------------------------------------------------------
!
!*       4.BIS    IBM and Fluctuations application
!              -----------------------------
!
!*       4.B1   Add fluctuations at the domain boundaries
!
IF (LRECYCL) THEN
  CALL ADDFLUCTUATIONS (                                    &
           CLBCX,CLBCY,                                     &
           XUT, XVT, XWT, XTHT, XTKET, XRT, XSVT, XSRCT,    &
           XFLUCTUTN,XFLUCTVTW,XFLUCTUTS,XFLUCTVTE,         &
           XFLUCTWTW,XFLUCTWTN,XFLUCTWTS,XFLUCTWTE          )
ENDIF
!
!*       4.B2   Immersed boundaries
!
IF (LIBM) THEN
  !
  ZTIME1=ZTIME2
  !
  IF (.NOT.LCARTESIAN) THEN
    CALL PRINT_MSG(NVERB_FATAL, 'GEN', 'MODELN', 'IBM can only be used in combination with cartesian coordinates')
  ENDIF
  !
  CALL IBM_FORCING(XUT,XVT,XWT,XTHT,XRT,XSVT,XTKET)
  !
  IF (LIBM_TROUBLE) THEN
     CALL IBM_FORCING_TR(XUT,XVT,XWT,XTHT,XRT,XSVT,XTKET)
  ENDIF
  !
  CALL SECOND_MNH2(ZTIME2)
  !
  XT_IBM_FORC = XT_IBM_FORC + ZTIME2 - ZTIME1
  !
ENDIF
!-------------------------------------------------------------------------------
!
!*       5.    INITIALIZATION OF THE BUDGET VARIABLES
!              --------------------------------------
!
IF (NBUMOD==IMI) THEN
  LBU_ENABLE = CBUTYPE /='NONE'.AND. CBUTYPE /='SKIP' 
ELSE
  LBU_ENABLE = .FALSE.
END IF
!
IF (NBUMOD==IMI .AND. CBUTYPE=='MASK' ) THEN
  CALL SET_MASK()
  if ( lbu_ru ) then
    tbudgets(NBUDGET_U)%trhodj%xdata(:, nbutime, :) = tbudgets(NBUDGET_U)%trhodj%xdata(:, nbutime, :) &
                                                      + Mask_compress( Mxm( xrhodj(:, :, :) ) )
  end if
  if ( lbu_rv ) then
    tbudgets(NBUDGET_V)%trhodj%xdata(:, nbutime, :) = tbudgets(NBUDGET_V)%trhodj%xdata(:, nbutime, :) &
                                                      + Mask_compress( Mym( xrhodj(:, :, :) ) )
  end if
  if ( lbu_rw ) then
    tbudgets(NBUDGET_W)%trhodj%xdata(:, nbutime, :) = tbudgets(NBUDGET_W)%trhodj%xdata(:, nbutime, :) &
                                                      + Mask_compress( Mzm( xrhodj(:, :, :) ) )
  end if
  if ( associated( tburhodj ) ) tburhodj%xdata(:, nbutime, :) = tburhodj%xdata(:, nbutime, :) + Mask_compress( xrhodj(:, :, :) )
END IF
!
IF (NBUMOD==IMI .AND. CBUTYPE=='CART' ) THEN
  if ( lbu_ru ) then
    tbudgets(NBUDGET_U)%trhodj%xdata(:, :, :) = tbudgets(NBUDGET_U)%trhodj%xdata(:, :, :) + Cart_compress( Mxm( xrhodj(:, :, :) ) )
  end if
  if ( lbu_rv ) then
    tbudgets(NBUDGET_V)%trhodj%xdata(:, :, :) = tbudgets(NBUDGET_V)%trhodj%xdata(:, :, :) + Cart_compress( Mym( xrhodj(:, :, :) ) )
  end if
  if ( lbu_rw ) then
    tbudgets(NBUDGET_W)%trhodj%xdata(:, :, :) = tbudgets(NBUDGET_W)%trhodj%xdata(:, :, :) &
                                                + Cart_compress( Mzm( xrhodj(:, :, :) ) )
  end if
  if ( associated( tburhodj ) ) tburhodj%xdata(:, :, :) = tburhodj%xdata(:, :, :) + Cart_compress( xrhodj(:, :, :) )
END IF
!
CALL BUDGET_FLAGS(LUSERV, LUSERC, LUSERR,         &
                  LUSERI, LUSERS, LUSERG, LUSERH  )
!
XTIME_BU   = 0.0
!
!-------------------------------------------------------------------------------
!
!*       6.    INITIALIZATION OF THE FIELD TENDENCIES
!              --------------------------------------
!
ZTIME1 = ZTIME2
XTIME_BU_PROCESS = 0.
!
CALL INITIAL_GUESS ( NRR, NSV, KTCOUNT, XRHODJ,IMI, XTSTEP,                 &
                     XRUS, XRVS, XRWS, XRTHS, XRRS, XRTKES, XRSVS,          &
                     XUT, XVT, XWT, XTHT, XRT, XTKET, XSVT )
!
CALL SECOND_MNH2(ZTIME2)
!
XT_GUESS = XT_GUESS + ZTIME2 - ZTIME1 - XTIME_BU_PROCESS
!
!-------------------------------------------------------------------------------
!
!*       7.    INITIALIZATION OF THE LES FOR CURRENT TIME-STEP
!              -----------------------------------------------
!
XTIME_LES_BU   = 0.0
XTIME_LES      = 0.0
IF (LLES) CALL LES_INI_TIMESTEP_n(KTCOUNT)
!
!-------------------------------------------------------------------------------
!
!*       8.    TWO-WAY INTERACTIVE GRID-NESTING
!              --------------------------------
!
!
CALL SECOND_MNH2(ZTIME1)
XTIME_BU_PROCESS = 0.
XTIME_LES_BU_PROCESS = 0.
!
GMASKkids(:,:)=.FALSE.
!
IF (NMODEL>1) THEN
  ! correct an ifort bug
  DPTR_XRHODJ=>XRHODJ
  DPTR_XUM=>XUT
  DPTR_XVM=>XVT
  DPTR_XWM=>XWT
  DPTR_XTHM=>XTHT
  DPTR_XRM=>XRT
  DPTR_XTKEM=>XTKET
  DPTR_XSVM=>XSVT
  DPTR_XRUS=>XRUS
  DPTR_XRVS=>XRVS
  DPTR_XRWS=>XRWS
  DPTR_XRTHS=>XRTHS
  DPTR_XRRS=>XRRS
  DPTR_XRTKES=>XRTKES
  DPTR_XRSVS=>XRSVS
  DPTR_XINPRC=>XINPRC
  DPTR_XINPRR=>XINPRR
  DPTR_XINPRS=>XINPRS
  DPTR_XINPRG=>XINPRG
  DPTR_XINPRH=>XINPRH
  DPTR_XPRCONV=>XPRCONV
  DPTR_XPRSCONV=>XPRSCONV
  DPTR_XDIRFLASWD=>XDIRFLASWD
  DPTR_XSCAFLASWD=>XSCAFLASWD
  DPTR_XDIRSRFSWD=>XDIRSRFSWD
  DPTR_GMASKkids=>GMASKkids
  !
  CALL TWO_WAY(     NRR,NSV,KTCOUNT,DPTR_XRHODJ,IMI,XTSTEP,                                    &
       DPTR_XUM ,DPTR_XVM ,DPTR_XWM , DPTR_XTHM, DPTR_XRM,DPTR_XSVM,                           &
       DPTR_XRUS,DPTR_XRVS,DPTR_XRWS,DPTR_XRTHS,DPTR_XRRS,DPTR_XRSVS,                          &
       DPTR_XINPRC,DPTR_XINPRR,DPTR_XINPRS,DPTR_XINPRG,DPTR_XINPRH,DPTR_XPRCONV,DPTR_XPRSCONV, &
       DPTR_XDIRFLASWD,DPTR_XSCAFLASWD,DPTR_XDIRSRFSWD,DPTR_GMASKkids           )
END IF
!
CALL SECOND_MNH2(ZTIME2)
XT_2WAY = XT_2WAY + ZTIME2 - ZTIME1 - XTIME_LES_BU_PROCESS - XTIME_BU_PROCESS
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
!*       10.    FORCING
!               -------
!
!
ZTIME1 = ZTIME2
XTIME_BU_PROCESS = 0.
XTIME_LES_BU_PROCESS = 0.
IF (LCARTESIAN) THEN
  CALL SM_GRIDCART(XXHAT,XYHAT,XZHAT,XZS,LSLEVE,XLEN1,XLEN2,XZSMT,XDXHAT,XDYHAT,XZZ,ZJ)
  XMAP=1.
ELSE
  CALL SM_GRIDPROJ( XXHAT, XYHAT, XZHAT, XXHATM, XYHATM, XZS,      &
                    LSLEVE, XLEN1, XLEN2, XZSMT, XLATORI, XLONORI, &
                    XMAP, XLAT, XLON, XDXHAT, XDYHAT, XZZ, ZJ      )
END IF
!
IF ( LFORCING ) THEN
  CALL FORCING(XTSTEP,LUSERV,XRHODJ,XCORIOZ,XZHAT,XZZ,TDTCUR,&
               XUFRC_PAST, XVFRC_PAST,XWTFRC,         &
               XUT,XVT,XWT,XTHT,XTKET,XRT,XSVT,       &
               XRUS,XRVS,XRWS,XRTHS,XRTKES,XRRS,XRSVS,IMI,ZJ)
END IF
!
IF ( L2D_ADV_FRC ) THEN 
  CALL ADV_FORCING_n(XRHODJ,TDTCUR,XTHT,XRT,XZZ,XRTHS,XRRS)
END IF
IF ( L2D_REL_FRC ) THEN 
  CALL REL_FORCING_n(XRHODJ,TDTCUR,XTHT,XRT,XZZ,XRTHS,XRRS)
END IF
!
CALL SECOND_MNH2(ZTIME2)
!
XT_FORCING = XT_FORCING + ZTIME2 - ZTIME1 &
             - XTIME_LES_BU_PROCESS - XTIME_BU_PROCESS
!
!-------------------------------------------------------------------------------
!
!*       11.    NUDGING
!               -------
!
!
ZTIME1 = ZTIME2
XTIME_BU_PROCESS = 0.
XTIME_LES_BU_PROCESS = 0.
!
IF ( LNUDGING ) THEN
  CALL NUDGING(LUSERV,XRHODJ,XTNUDGING,         &
               XUT,XVT,XWT,XTHT,XRT,            &
               XLSUM,XLSVM,XLSWM,XLSTHM,XLSRVM, &
               XRUS,XRVS,XRWS,XRTHS,XRRS)

END IF
!
CALL SECOND_MNH2(ZTIME2)
!
XT_NUDGING = XT_NUDGING + ZTIME2 - ZTIME1 &
             - XTIME_LES_BU_PROCESS - XTIME_BU_PROCESS
!
!-------------------------------------------------------------------------------
!
!*       12.    DYNAMICAL SOURCES
!               -----------------
!
ZTIME1 = ZTIME2
XTIME_BU_PROCESS = 0.
XTIME_LES_BU_PROCESS = 0.
!
IF( LTRANS ) THEN
  XUT(:,:,:) = XUT(:,:,:) + XUTRANS
  XVT(:,:,:) = XVT(:,:,:) + XVTRANS
END IF
!
CALL DYN_SOURCES( NRR,NRRL, NRRI,                              &
                  XUT, XVT, XWT, XTHT, XRT,                    &
                  XCORIOX, XCORIOY, XCORIOZ, XCURVX, XCURVY,   &
                  XRHODJ, XZZ, XTHVREF, XEXNREF,               &
                  XRUS, XRVS, XRWS, XRTHS                      )
!
IF( LTRANS ) THEN
  XUT(:,:,:) = XUT(:,:,:) - XUTRANS
  XVT(:,:,:) = XVT(:,:,:) - XVTRANS
END IF
!
CALL SECOND_MNH2(ZTIME2)
!
XT_SOURCES = XT_SOURCES + ZTIME2 - ZTIME1 &
             - XTIME_LES_BU_PROCESS - XTIME_BU_PROCESS
!
!-------------------------------------------------------------------------------
!
!*       13.    NUMERICAL DIFFUSION
!               -------------------
!
ZTIME1 = ZTIME2
XTIME_BU_PROCESS = 0.
XTIME_LES_BU_PROCESS = 0.
!
IF ( LNUMDIFU .OR. LNUMDIFTH .OR. LNUMDIFSV ) THEN
!
  CALL UPDATE_HALO_ll(TFIELDT_ll, IINFO_ll)
  CALL UPDATE_HALO2_ll(TFIELDT_ll, THALO2T_ll, IINFO_ll)
  IF ( .NOT. LSTEADYLS ) THEN
     CALL UPDATE_HALO_ll(TLSFIELD_ll, IINFO_ll)
     CALL UPDATE_HALO_ll(TLSFIELD2D_ll,IINFO_ll) 
     CALL UPDATE_HALO2_ll(TLSFIELD_ll, TLSHALO2_ll, IINFO_ll)
  END IF
  CALL NUM_DIFF ( CLBCX, CLBCY, NRR, NSV,                               &
                  XDK2U, XDK4U, XDK2TH, XDK4TH, XDK2SV, XDK4SV, IMI,    &
                  XUT, XVT, XWT, XTHT, XTKET, XRT, XSVT,                &
                  XLSUM,XLSVM,XLSWM,XLSTHM,XLSRVM,XRHODJ,               &
                  XRUS, XRVS, XRWS, XRTHS, XRTKES, XRRS, XRSVS,         &
                  LZDIFFU,LNUMDIFU, LNUMDIFTH, LNUMDIFSV,               &
                  THALO2T_ll, TLSHALO2_ll,XZDIFFU_HALO2      )
END IF

if ( lbudget_sv ) then
  do jsv = 1, nsv
    call Budget_store_init( tbudgets(jsv + NBUDGET_SV1 - 1), 'NEGA2', xrsvs(:, :, :, jsv) )
  end do
end if

DO JSV = NSV_CHEMBEG,NSV_CHEMEND
  XRSVS(:,:,:,JSV) = MAX(XRSVS(:,:,:,JSV),0.)
END DO
DO JSV = NSV_CHICBEG,NSV_CHICEND
  XRSVS(:,:,:,JSV) = MAX(XRSVS(:,:,:,JSV),0.)
END DO
DO JSV = NSV_AERBEG,NSV_AEREND
  XRSVS(:,:,:,JSV) = MAX(XRSVS(:,:,:,JSV),0.)
END DO
DO JSV = NSV_LNOXBEG,NSV_LNOXEND
  XRSVS(:,:,:,JSV) = MAX(XRSVS(:,:,:,JSV),0.)
END DO
DO JSV = NSV_DSTBEG,NSV_DSTEND
  XRSVS(:,:,:,JSV) = MAX(XRSVS(:,:,:,JSV),0.)
END DO
DO JSV = NSV_SLTBEG,NSV_SLTEND
  XRSVS(:,:,:,JSV) = MAX(XRSVS(:,:,:,JSV),0.)
END DO
DO JSV = NSV_PPBEG,NSV_PPEND
  XRSVS(:,:,:,JSV) = MAX(XRSVS(:,:,:,JSV),0.)
END DO
#ifdef MNH_FOREFIRE
DO JSV = NSV_FFBEG,NSV_FFEND
  XRSVS(:,:,:,JSV) = MAX(XRSVS(:,:,:,JSV),0.)
END DO
#endif
! Blaze smoke
DO JSV = NSV_FIREBEG,NSV_FIREEND
  XRSVS(:,:,:,JSV) = MAX(XRSVS(:,:,:,JSV),0.)
END DO
DO JSV = NSV_CSBEG,NSV_CSEND
  XRSVS(:,:,:,JSV) = MAX(XRSVS(:,:,:,JSV),0.)
END DO
DO JSV = NSV_DSTDEPBEG,NSV_DSTDEPEND
  XRSVS(:,:,:,JSV) = MAX(XRSVS(:,:,:,JSV),0.)
END DO
DO JSV = NSV_SLTDEPBEG,NSV_SLTDEPEND
  XRSVS(:,:,:,JSV) = MAX(XRSVS(:,:,:,JSV),0.)
END DO
DO JSV = NSV_AERDEPBEG,NSV_AERDEPEND
  XRSVS(:,:,:,JSV) = MAX(XRSVS(:,:,:,JSV),0.)
END DO
DO JSV = NSV_SNWBEG,NSV_SNWEND
  XRSVS(:,:,:,JSV) = MAX(XRSVS(:,:,:,JSV),0.)
END DO
IF (CELEC .NE. 'NONE') THEN
  XRSVS(:,:,:,NSV_ELECBEG) = MAX(XRSVS(:,:,:,NSV_ELECBEG),0.)
  XRSVS(:,:,:,NSV_ELECEND) = MAX(XRSVS(:,:,:,NSV_ELECEND),0.)
END IF

if ( lbudget_sv ) then
  do jsv = 1, nsv
    call Budget_store_end( tbudgets(jsv + NBUDGET_SV1 - 1), 'NEGA2', xrsvs(:, :, :, jsv) )
  end do
end if
!
CALL SECOND_MNH2(ZTIME2)
!
XT_DIFF = XT_DIFF + ZTIME2 - ZTIME1 &
          - XTIME_LES_BU_PROCESS - XTIME_BU_PROCESS
!
!-------------------------------------------------------------------------------
!
!*       14.    UPPER AND LATERAL RELAXATION
!               ----------------------------
!
ZTIME1 = ZTIME2
XTIME_BU_PROCESS = 0.
XTIME_LES_BU_PROCESS = 0.
!
IF(LVE_RELAX .OR. LVE_RELAX_GRD .OR. LHORELAX_UVWTH .OR. LHORELAX_RV .OR.&
   LHORELAX_RC .OR. LHORELAX_RR .OR. LHORELAX_RI .OR. LHORELAX_RS .OR.   &
   LHORELAX_RG .OR. LHORELAX_RH .OR. LHORELAX_TKE .OR.                   &
   ANY(LHORELAX_SV)) THEN
  CALL RELAXATION (LVE_RELAX,LVE_RELAX_GRD,LHORELAX_UVWTH,LHORELAX_RV,LHORELAX_RC,   &
                   LHORELAX_RR,LHORELAX_RI,LHORELAX_RS,LHORELAX_RG,    &
                   LHORELAX_RH,LHORELAX_TKE,LHORELAX_SV,               &
                   LHORELAX_SVC2R2,LHORELAX_SVC1R3,                    &
                   LHORELAX_SVELEC,LHORELAX_SVLG,                      &
                   LHORELAX_SVCHEM,LHORELAX_SVCHIC,LHORELAX_SVAER,     &
                   LHORELAX_SVDST,LHORELAX_SVSLT,LHORELAX_SVPP,        &
                   LHORELAX_SVCS,LHORELAX_SVSNW,LHORELAX_SVFIRE,       &
#ifdef MNH_FOREFIRE
                   LHORELAX_SVFF,                                      &
#endif
                   KTCOUNT,NRR,NSV,XTSTEP,XRHODJ,                      &
                   XUT, XVT, XWT, XTHT, XRT, XSVT, XTKET,              &
                   XLSUM, XLSVM, XLSWM, XLSTHM,                        &
                   XLBXUM, XLBXVM, XLBXWM, XLBXTHM,                    &
                   XLBXRM, XLBXSVM, XLBXTKEM,                          &
                   XLBYUM, XLBYVM, XLBYWM, XLBYTHM,                    &
                   XLBYRM, XLBYSVM, XLBYTKEM,                          &
                   NALBOT, XALK, XALKW,                                &
                   NALBAS, XALKBAS, XALKWBAS,                          &
                   LMASK_RELAX,XKURELAX, XKVRELAX, XKWRELAX,           &
                   NRIMX,NRIMY,                                        &
                   XRUS, XRVS, XRWS, XRTHS, XRRS, XRSVS, XRTKES        )
END IF

IF (CELEC.NE.'NONE' .AND. LRELAX2FW_ION) THEN
   CALL RELAX2FW_ION (KTCOUNT, IMI, XTSTEP, XRHODJ, XSVT, NALBOT,      &
                      XALK, LMASK_RELAX, XKWRELAX, XRSVS )   
END IF
!
CALL SECOND_MNH2(ZTIME2)
!
XT_RELAX = XT_RELAX + ZTIME2 - ZTIME1 &
           - XTIME_LES_BU_PROCESS - XTIME_BU_PROCESS
!
!-------------------------------------------------------------------------------
!
!*       15.    PARAMETRIZATIONS' MONITOR
!               -------------------------
!
ZTIME1 = ZTIME2
!
CALL PHYS_PARAM_n( KTCOUNT, TPBAKFILE,                            &
                   XT_RAD,  XT_SHADOWS, XT_DCONV, XT_GROUND,      &
                   XT_MAFL, XT_DRAG, XT_EOL, XT_TURB,  XT_TRACER, &
                   ZTIME, ZWETDEPAER, GMASKkids, GCLOUD_ONLY      )
!
IF (CDCONV/='NONE') THEN
  XPACCONV = XPACCONV + XPRCONV * XTSTEP
  IF (LCH_CONV_LINOX) THEN
    XIC_TOTAL_NUMBER = XIC_TOTAL_NUMBER + XIC_RATE * XTSTEP
    XCG_TOTAL_NUMBER = XCG_TOTAL_NUMBER + XCG_RATE * XTSTEP
  END IF
END IF
!
!
CALL SECOND_MNH2(ZTIME2)
!
XT_PARAM = XT_PARAM + ZTIME2 - ZTIME1 - XTIME_LES - ZTIME
!
!-------------------------------------------------------------------------------
!
!*       16.    TEMPORAL SERIES
!               ---------------
!
ZTIME1 = ZTIME2
!
IF (LSERIES) THEN
  IF ( MOD (KTCOUNT-1,NFREQSERIES) == 0 ) CALL SERIES_n
END IF
!
CALL SECOND_MNH2(ZTIME2)
!
XT_STEP_MISC = XT_STEP_MISC + ZTIME2 - ZTIME1
!
!-------------------------------------------------------------------------------
!
!*       17.    LARGE SCALE FIELD REFRESH
!               -------------------------
!
ZTIME1 = ZTIME2
!
IF (.NOT. LSTEADYLS) THEN
  IF (  IMI==1                             .AND.      &
    NCPL_CUR < NCPL_NBR                              ) THEN
    IF (KTCOUNT+1 == NCPL_TIMES(NCPL_CUR,1)          ) THEN
                                  ! The next current time reachs a
      NCPL_CUR=NCPL_CUR+1         ! coupling one, LS sources are refreshed
      !
      CALL LS_COUPLING(XTSTEP,GSTEADY_DMASS,CCONF,                          &
             CGETTKET,                                                      &
             CGETRVT,CGETRCT,CGETRRT,CGETRIT,                               &
             CGETRST,CGETRGT,CGETRHT,CGETSVT,LCH_INIT_FIELD, NSV,           &
             NIMAX_ll,NJMAX_ll,                                             &
             NSIZELBX_ll,NSIZELBXU_ll,NSIZELBY_ll,NSIZELBYV_ll,             &
             NSIZELBXTKE_ll,NSIZELBYTKE_ll,                                 &
             NSIZELBXR_ll,NSIZELBYR_ll,NSIZELBXSV_ll,NSIZELBYSV_ll,         &
             XLSUM,XLSVM,XLSWM,XLSTHM,XLSRVM,XLSZWSM,XDRYMASST,             &
             XLBXUM,XLBXVM,XLBXWM,XLBXTHM,XLBXTKEM,XLBXRM,XLBXSVM,          &
             XLBYUM,XLBYVM,XLBYWM,XLBYTHM,XLBYTKEM,XLBYRM,XLBYSVM,          &
             XLSUS,XLSVS,XLSWS,XLSTHS,XLSRVS,XLSZWSS,XDRYMASSS,             &
             XLBXUS,XLBXVS,XLBXWS,XLBXTHS,XLBXTKES,XLBXRS,XLBXSVS,          &
             XLBYUS,XLBYVS,XLBYWS,XLBYTHS,XLBYTKES,XLBYRS,XLBYSVS           )
      !
      DO JSV=NSV_CHEMBEG,NSV_CHEMEND
        XLBXSVS(:,:,:,JSV)=MAX(XLBXSVS(:,:,:,JSV),0.)
        XLBYSVS(:,:,:,JSV)=MAX(XLBYSVS(:,:,:,JSV),0.)
      ENDDO
      !
      DO JSV=NSV_LNOXBEG,NSV_LNOXEND
        XLBXSVS(:,:,:,JSV)=MAX(XLBXSVS(:,:,:,JSV),0.)
        XLBYSVS(:,:,:,JSV)=MAX(XLBYSVS(:,:,:,JSV),0.)
      ENDDO
      !
      DO JSV=NSV_AERBEG,NSV_AEREND
        XLBXSVS(:,:,:,JSV)=MAX(XLBXSVS(:,:,:,JSV),0.)
        XLBYSVS(:,:,:,JSV)=MAX(XLBYSVS(:,:,:,JSV),0.)
      ENDDO
      !
      DO JSV=NSV_DSTBEG,NSV_DSTEND
        XLBXSVS(:,:,:,JSV)=MAX(XLBXSVS(:,:,:,JSV),0.)
        XLBYSVS(:,:,:,JSV)=MAX(XLBYSVS(:,:,:,JSV),0.)
      ENDDO
      !
      DO JSV=NSV_DSTDEPBEG,NSV_DSTDEPEND
        XLBXSVS(:,:,:,JSV)=MAX(XLBXSVS(:,:,:,JSV),0.)
        XLBYSVS(:,:,:,JSV)=MAX(XLBYSVS(:,:,:,JSV),0.)
      ENDDO
      !
      DO JSV=NSV_SLTBEG,NSV_SLTEND
        XLBXSVS(:,:,:,JSV)=MAX(XLBXSVS(:,:,:,JSV),0.)
        XLBYSVS(:,:,:,JSV)=MAX(XLBYSVS(:,:,:,JSV),0.)
      ENDDO
      !
      DO JSV=NSV_SLTDEPBEG,NSV_SLTDEPEND
        XLBXSVS(:,:,:,JSV)=MAX(XLBXSVS(:,:,:,JSV),0.)
        XLBYSVS(:,:,:,JSV)=MAX(XLBYSVS(:,:,:,JSV),0.)
      ENDDO
      !
      DO JSV=NSV_PPBEG,NSV_PPEND
        XLBXSVS(:,:,:,JSV)=MAX(XLBXSVS(:,:,:,JSV),0.)
        XLBYSVS(:,:,:,JSV)=MAX(XLBYSVS(:,:,:,JSV),0.)
      ENDDO
      !
#ifdef MNH_FOREFIRE
      DO JSV=NSV_FFBEG,NSV_FFEND
        XLBXSVS(:,:,:,JSV)=MAX(XLBXSVS(:,:,:,JSV),0.)
        XLBYSVS(:,:,:,JSV)=MAX(XLBYSVS(:,:,:,JSV),0.)
      ENDDO
      !
#endif
      DO JSV=NSV_FIREBEG,NSV_FIREEND
        XLBXSVS(:,:,:,JSV)=MAX(XLBXSVS(:,:,:,JSV),0.)
        XLBYSVS(:,:,:,JSV)=MAX(XLBYSVS(:,:,:,JSV),0.)
      ENDDO
      !
      DO JSV=NSV_CSBEG,NSV_CSEND
        XLBXSVS(:,:,:,JSV)=MAX(XLBXSVS(:,:,:,JSV),0.)
        XLBYSVS(:,:,:,JSV)=MAX(XLBYSVS(:,:,:,JSV),0.)
      ENDDO
      !
      DO JSV=NSV_SNWBEG,NSV_SNWEND
        XLBXSVS(:,:,:,JSV)=MAX(XLBXSVS(:,:,:,JSV),0.)
        XLBYSVS(:,:,:,JSV)=MAX(XLBYSVS(:,:,:,JSV),0.)
      ENDDO
      !
     END IF
  END IF
END IF
!
CALL SECOND_MNH2(ZTIME2)
!
XT_COUPL = XT_COUPL + ZTIME2 - ZTIME1
!
!-------------------------------------------------------------------------------
!
!
!
!*      8 Bis . Blowing snow scheme
!              ---------
!
IF ( LBLOWSNOW ) THEN
 CALL BLOWSNOW( XTSTEP, NRR, XPABST, XTHT, XRT, XZZ, XRHODREF,        &
                XRHODJ, XEXNREF, XRRS, XRTHS, XSVT, XRSVS, XSNWSUBL3D )
ENDIF
!
!-----------------------------------------------------------------------
!
!*       8 Ter  VISCOSITY (no-slip condition inside)
!              ---------
!
!
IF ( LVISC ) THEN
!
ZTIME1 = ZTIME2
!
   CALL VISCOSITY(CLBCX, CLBCY, NRR, NSV, XMU_V,XPRANDTL,         &
                  LVISC_UVW,LVISC_TH,LVISC_SV,LVISC_R,            &
                  LDRAG,    &
                  XUT, XVT, XWT, XTHT, XRT, XSVT,                 &
                  XRHODJ, XDXX, XDYY, XDZZ, XDZX, XDZY,           &
                  XRUS, XRVS, XRWS, XRTHS, XRRS, XRSVS,XDRAG )
!
ENDIF
!
CALL SECOND_MNH2(ZTIME2)
!
XT_VISC = XT_VISC + ZTIME2 - ZTIME1
!!
!-------------------------------------------------------------------------------
!
!*       9.    ADVECTION
!              ---------
!
ZTIME1 = ZTIME2
XTIME_BU_PROCESS = 0.
XTIME_LES_BU_PROCESS = 0.
!
!
!$acc data create ( XUT, XVT, XWT )                                             &
!$acc &    present( XTHT, XRT, XPABST, XTHVREF, XRHODJ )                        &
!$acc &    present( XDXX, XDYY, XDZZ, XDZX, XDZY, XRUS, XRVS, XRWS, XRWS_PRES ) &
!$acc &    present( XRRS, XRSVS, XRTHS_CLD, XRRS_CLD, XRSVS_CLD )               &
!$acc &    copyin ( XTKET, XSVT )
!
CALL MPPDB_CHECK3DM("before ADVEC_METSV:XU/V/W/TH/TKE/T,XRHODJ",PRECISION,&
                   &  XUT, XVT, XWT, XTHT, XTKET,XRHODJ)
!$acc update device( XRUS, XRVS, XRWS, XRRS, XRSVS )
!$acc update device( XRUS_PRES, XRVS_PRES, XRWS_PRES )
!$acc update device( XUT, XVT, XWT, XTHT, XRT, XRHODJ, XRTHS )
!$acc update device( XRTHS_CLD, XRRS_CLD, XRSVS_CLD )
!
!$acc data copy   (XRTKES)  &
!$acc &    copyout(XRTKEMS)
#ifdef MNH_BITREP_OMP
CALL SBR_FZ(XRRS_CLD)
CALL SBR_FZ(XRT)
#endif
CALL ADVECTION_METSV ( TPBAKFILE, CUVW_ADV_SCHEME,                     &
                 CMET_ADV_SCHEME, CSV_ADV_SCHEME, CCLOUD, CELEC,       &
                 NSPLIT, LSPLIT_CFL, XSPLIT_CFL, LCFL_WRIT,            &
                 CLBCX, CLBCY, NRR, NSV, TDTCUR, XTSTEP,               &
                 XUT, XVT, XWT, XTHT, XRT, XTKET, XSVT, XPABST,        &
                 XTHVREF, XRHODJ, XDXX, XDYY, XDZZ, XDZX, XDZY,        &
                 XRTHS, XRRS, XRTKES, XRSVS,                           &
                 XRTHS_CLD, XRRS_CLD, XRSVS_CLD, XRTKEMS               )
!$acc end data
!
!
CALL SECOND_MNH2(ZTIME2)
!
XT_ADV = XT_ADV + ZTIME2 - ZTIME1 - XTIME_LES_BU_PROCESS - XTIME_BU_PROCESS
!
!-------------------------------------------------------------------------------
ZTIME1 = ZTIME2
XTIME_BU_PROCESS = 0.
XTIME_LES_BU_PROCESS = 0.
!
!$acc kernels present( ZRWS )
ZRWS(:,:,:) = XRWS(:,:,:)
!$acc end kernels
!
CALL GRAVITY_IMPL ( CLBCX, CLBCY, NRR, NRRL, NRRI,XTSTEP,            &
                 XTHT, XRT, XTHVREF, XRHODJ, XRWS, XRTHS, XRRS,      &
                 XRTHS_CLD, XRRS_CLD                                 )
!
! At the initial instant the difference with the ref state creates a
! vertical velocity production that must not be advected as it is
! compensated by the pressure gradient
!
IF (KTCOUNT == 1 .AND. CCONF=='START') THEN
!$acc kernels present( ZRWS,XRWS_PRES )
  XRWS_PRES(:,:,:) = ZRWS(:,:,:) - XRWS(:,:,:)
!$acc end kernels
END IF
!
CALL SECOND_MNH2(ZTIME2)
!
XT_GRAV = XT_GRAV + ZTIME2 - ZTIME1 - XTIME_LES_BU_PROCESS - XTIME_BU_PROCESS
! 
IF ( LIBM .AND. CIBM_ADV=='FORCIN' ) THEN
  !
  ZTIME1=ZTIME2
  !
  CALL IBM_FORCING_ADV   (XRUS,XRVS,XRWS)
  !
  CALL SECOND_MNH2(ZTIME2)
  !
  XT_IBM_FORC = XT_IBM_FORC + ZTIME2 - ZTIME1
  !
ENDIF
!
ZTIME1 = ZTIME2
XTIME_BU_PROCESS = 0.
XTIME_LES_BU_PROCESS = 0.
!
IF ((CUVW_ADV_SCHEME(1:3)=='CEN') .AND. (CTEMP_SCHEME == 'LEFR')) THEN
  IF (CUVW_ADV_SCHEME=='CEN4TH') THEN
    NULLIFY(TZFIELDC_ll)
    NULLIFY(TZHALO2C_ll)
      CALL ADD3DFIELD_ll( TZFIELDC_ll, XUT, 'MODEL_n::XUT' )
      CALL ADD3DFIELD_ll( TZFIELDC_ll, XVT, 'MODEL_n::XVT' )
      CALL ADD3DFIELD_ll( TZFIELDC_ll, XWT, 'MODEL_n::XWT' )
#ifndef MNH_OPENACC
      CALL INIT_HALO2_ll( TZHALO2C_ll, 3, IIU, IJU, IKU, 'MODEL_n::HALO2C' )
      CALL UPDATE_HALO_ll(TZFIELDC_ll,IINFO_ll)
      CALL UPDATE_HALO2_ll(TZFIELDC_ll, TZHALO2C_ll, IINFO_ll)
#else
  IF (GFIRST_CALL_MODELN) THEN
    GFIRST_CALL_MODELN = .FALSE.
    NULLIFY(TZHALO2_UT,TZHALO2_VT,TZHALO2_WT)
    CALL INIT_HALO2_ll( TZHALO2_UT, 1, IIU, IJU, IKU, 'MODEL_n::XUT' )
    CALL INIT_HALO2_ll( TZHALO2_VT, 1, IIU, IJU, IKU, 'MODEL_n::XVT' )
    CALL INIT_HALO2_ll( TZHALO2_WT, 1, IIU, IJU, IKU, 'MODEL_n::XWT' )
  END IF

  CALL GET_HALO2_DF(XUT,TZHALO2_UT,HNAME='XUT')
  CALL GET_HALO2_DF(XVT,TZHALO2_VT,HNAME='XVT')
  CALL GET_HALO2_DF(XWT,TZHALO2_WT,HNAME='XWT')
#endif
!$acc update device(XUT, XVT, XWT)
  END IF
!$acc data copyin(XUM, XVM, XWM)    &
!$acc &    copy  (XDUM, XDVM, XDWM)
 CALL ADVECTION_UVW_CEN(CUVW_ADV_SCHEME,                           &
                           CLBCX, CLBCY,                           &
                           XTSTEP, KTCOUNT,                        &
                           XUM, XVM, XWM, XDUM, XDVM, XDWM,        &
                           XUT, XVT, XWT,                          &
                           XRHODJ, XDXX, XDYY, XDZZ, XDZX, XDZY,   &
                           XRUS,XRVS, XRWS,                        &
#ifndef MNH_OPENACC
                           TZHALO2C_ll                             )
#else
                           TZHALO2_UT,TZHALO2_VT,TZHALO2_WT        )
#endif
!$acc end data
  IF (CUVW_ADV_SCHEME=='CEN4TH') THEN
    CALL CLEANLIST_ll(TZFIELDC_ll)
    NULLIFY(TZFIELDC_ll)
#ifndef MNH_OPENACC
    CALL  DEL_HALO2_ll(TZHALO2C_ll)
    NULLIFY(TZHALO2C_ll)
#endif
  END IF
ELSE

  CALL ADVECTION_UVW(CUVW_ADV_SCHEME, CTEMP_SCHEME,                  &
                 NWENO_ORDER, LSPLIT_WENO,                           &
                 CLBCX, CLBCY, XTSTEP,                               &
                 XUT, XVT, XWT,                                      &
                 XRHODJ, XDXX, XDYY, XDZZ, XDZX, XDZY,               &
                 XRUS, XRVS, XRWS,                                   &
                 XRUS_PRES, XRVS_PRES, XRWS_PRES                     )
END IF
!
!
CALL SECOND_MNH2(ZTIME2)
!
XT_ADVUVW = XT_ADVUVW + ZTIME2 - ZTIME1 - XTIME_LES_BU_PROCESS - XTIME_BU_PROCESS
!
!-------------------------------------------------------------------------------
!
IF (LCLOUDMODIFLM) THEN
  CALL TURB_CLOUD_INDEX( XTSTEP, TPBAKFILE,                               &
                         LTURB_DIAG, NRRI,                                &
                         XRRS, XRT, XRHODJ, XDXX, XDYY, XDZZ, XDZX, XDZY, &
                         XCEI                                             )
END IF
!
!-------------------------------------------------------------------------------
!
!*       18.    LATERAL BOUNDARY CONDITION FOR THE NORMAL VELOCITY
!               --------------------------------------------------
!
ZTIME1 = ZTIME2
!
!$acc kernels present(ZRUS,ZRVS,ZRWS)
ZRUS(:,:,:)=XRUS(:,:,:)
ZRVS(:,:,:)=XRVS(:,:,:)
ZRWS(:,:,:)=XRWS(:,:,:)
!$acc end kernels
!
if ( .not. l1d ) then
  if ( lbudget_u ) call Budget_store_init( tbudgets(NBUDGET_U), 'PRES', xrus(:, :, :) )
  if ( lbudget_v ) call Budget_store_init( tbudgets(NBUDGET_V), 'PRES', xrvs(:, :, :) )
  if ( lbudget_w ) call Budget_store_init( tbudgets(NBUDGET_W), 'PRES', xrws(:, :, :) )
end if
!
CALL MPPDB_CHECK3DM("before RAD_BOUND : other var",PRECISION,XUT,XVT,XRHODJ,XTKET)
CALL MPPDB_CHECKLB(XLBXUM,"modeln XLBXUM",PRECISION,'LBXU',NRIMX)
CALL MPPDB_CHECKLB(XLBYVM,"modeln XLBYVM",PRECISION,'LBYV',NRIMY)
CALL MPPDB_CHECKLB(XLBXUS,"modeln XLBXUS",PRECISION,'LBXU',NRIMX)
CALL MPPDB_CHECKLB(XLBYVS,"modeln XLBYVS",PRECISION,'LBYV',NRIMY)
!
!$acc data copyin( XDXHAT, XDYHAT, XZHAT, XLBXUM, XLBYVM, XLBXUS, XLBYVS, &
!$acc &            XFLUCTUNW, XFLUCTVNN, XFLUCTUNE, XFLUCTVNS )
CALL RAD_BOUND (CLBCX,CLBCY,CTURB,XCARPKMAX,             &
                XTSTEP,                                  &
                XDXHAT, XDYHAT, XZHAT,                   &
                XUT, XVT,                                &
                XLBXUM, XLBYVM, XLBXUS, XLBYVS,          &
                XFLUCTUNW,XFLUCTVNN,XFLUCTUNE,XFLUCTVNS, &
                XCPHASE, XCPHASE_PBL, XRHODJ,            &
                XTKET, XRUS, XRVS, XRWS                  )
!$acc end data
!$acc kernels present(ZRUS,ZRVS,ZRWS)
ZRUS(:,:,:)=XRUS(:,:,:)-ZRUS(:,:,:)
ZRVS(:,:,:)=XRVS(:,:,:)-ZRVS(:,:,:)
ZRWS(:,:,:)=XRWS(:,:,:)-ZRWS(:,:,:)
!$acc end kernels
!
CALL SECOND_MNH2(ZTIME2)
!
XT_RAD_BOUND = XT_RAD_BOUND + ZTIME2 - ZTIME1
!
!$acc end data
!-------------------------------------------------------------------------------
!
!*       19.    PRESSURE COMPUTATION
!               --------------------
!
ZTIME1 = ZTIME2
XTIME_BU_PROCESS = 0.
XTIME_LES_BU_PROCESS = 0.
!
!$acc kernels present( XPABST, ZPABST )
ZPABST(:,:,:) = XPABST(:,:,:)
!$acc end kernels
!
IF(.NOT. L1D) THEN
!
!$acc kernels
  XRUS_PRES(:,:,:) = XRUS(:,:,:)
  XRVS_PRES(:,:,:) = XRVS(:,:,:)
  XRWS_PRES(:,:,:) = XRWS(:,:,:)
!$acc end kernels
!
!$acc data present( XRHODJ, XDXX, XDYY, XDZZ, XDZX, XDZY )                    &
!$acc &    present( XRHOM, XAF, XBFY, XCF, XTRIGSX, XTRIGSY, NIFAXX, NIFAXY ) &
!$acc &    present( XTHT, XRT, XRHODREF, XTHVREF, XRVREF, XEXNREF )           &
!$acc &    present( XRUS, XRVS, XRWS, XPABST, XBFB, XBF_SXP2_YP1_Z )
  CALL PRESSUREZ( CLBCX,CLBCY,CPRESOPT,NITR,LITRADJ,KTCOUNT, XRELAX,IMI, &
                  XRHODJ,XDXX,XDYY,XDZZ,XDZX,XDZY,XDXHATM,XDYHATM,XRHOM, &
                  XAF,XBFY,XCF,XTRIGSX,XTRIGSY,NIFAXX,NIFAXY,            &
                  NRR,NRRL,NRRI,XDRYMASST,XREFMASS,XMASS_O_PHI0,         &
                  XTHT,XRT,XRHODREF,XTHVREF,XRVREF,XEXNREF, XLINMASS,    &
                  XRUS, XRVS, XRWS, XPABST,                              &
                  XBFB,                                                  &
                  XBF_SXP2_YP1_Z) !JUAN Z_SPLITING
!$acc end data
!$acc update self( XRUS, XRVS, XRWS, XPABST )
!
!$acc kernels
  XRUS_PRES(:,:,:) = XRUS(:,:,:) - XRUS_PRES(:,:,:) + ZRUS(:,:,:)
  XRVS_PRES(:,:,:) = XRVS(:,:,:) - XRVS_PRES(:,:,:) + ZRVS(:,:,:)
  XRWS_PRES(:,:,:) = XRWS(:,:,:) - XRWS_PRES(:,:,:) + ZRWS(:,:,:)
!$acc end kernels
!$acc update self( XRUS_PRES, XRVS_PRES, XRWS_PRES )
  CALL MPPDB_CHECK3DM("after pressurez:XRU/V/WS",PRECISION,XRUS,XRVS,XRWS)
!
END IF
!
CALL SECOND_MNH2(ZTIME2)
!
XT_PRESS = XT_PRESS + ZTIME2 - ZTIME1 &
           - XTIME_LES_BU_PROCESS - XTIME_BU_PROCESS
!
!-------------------------------------------------------------------------------
!
!*       20.    CHEMISTRY/AEROSOLS
!               ------------------
!
ZTIME1 = ZTIME2
XTIME_BU_PROCESS = 0.
XTIME_LES_BU_PROCESS = 0.
!
IF (LUSECHEM) THEN
  CALL CH_MONITOR_n(ZWETDEPAER,KTCOUNT,XTSTEP, ILUOUT, NVERB)
END IF
!
! For inert aerosol (dust and sea salt) => aer_monitor_n
IF ((LDUST).OR.(LSALT)) THEN
!
! tests to see if any cloud exists
!   
    GCLD=.TRUE.
    IF (GCLD .AND. NRR.LE.3 ) THEN 
      IF( MAX(MAXVAL(XCLDFR(:,:,:)),MAXVAL(XICEFR(:,:,:))).LE. 1.E-10 .AND. GCLOUD_ONLY ) THEN
          GCLD = .FALSE.                ! only the cloudy verticals would be 
                                        ! refreshed but there is no clouds 
      END IF
    END IF
!
    IF (GCLD .AND. NRR.GE.4 ) THEN 
      IF( CCLOUD(1:3)=='ICE' )THEN
        IF( MAXVAL(XRT(:,:,:,2)).LE.XRTMIN(2) .AND.             &
            MAXVAL(XRT(:,:,:,4)).LE.XRTMIN(4) .AND. GCLOUD_ONLY ) THEN
            GCLD = .FALSE.            ! only the cloudy verticals would be 
                                      ! refreshed but there is no cloudwater and ice
        END IF
      END IF
      IF( CCLOUD=='C3R5' )THEN
        IF( MAXVAL(XRT(:,:,:,2)).LE.XRTMIN_C1R3(2) .AND.             &
            MAXVAL(XRT(:,:,:,4)).LE.XRTMIN_C1R3(4) .AND. GCLOUD_ONLY ) THEN
            GCLD = .FALSE.            ! only the cloudy verticals would be 
                                      ! refreshed but there is no cloudwater and ice
        END IF
      END IF
      IF( CCLOUD=='LIMA' )THEN
        IF( MAXVAL(XRT(:,:,:,2)).LE.XRTMIN_LIMA(2) .AND.             &
            MAXVAL(XRT(:,:,:,4)).LE.XRTMIN_LIMA(4) .AND. GCLOUD_ONLY ) THEN
            GCLD = .FALSE.            ! only the cloudy verticals would be 
                                      ! refreshed but there is no cloudwater and ice
        END IF
      END IF
    END IF

!
        CALL AER_MONITOR_n(KTCOUNT,XTSTEP, ILUOUT, NVERB, GCLD)
END IF
!
!
CALL SECOND_MNH2(ZTIME2)
!
XT_CHEM = XT_CHEM + ZTIME2 - ZTIME1 &
      - XTIME_LES_BU_PROCESS - XTIME_BU_PROCESS
!
ZTIME = ZTIME + XTIME_LES_BU_PROCESS + XTIME_BU_PROCESS

!-------------------------------------------------------------------------------
!
!*       20.    WATER MICROPHYSICS
!               ------------------
!
ZTIME1 = ZTIME2
XTIME_BU_PROCESS = 0.
XTIME_LES_BU_PROCESS = 0.
!
IF (CCLOUD /= 'NONE') THEN
!
  IF (CCLOUD == 'C2R2' .OR. CCLOUD == 'KHKO' .OR. CCLOUD == 'C3R5' &
                                             .OR. CCLOUD == "LIMA" ) THEN
    IF ( LFORCING ) THEN
      XWT_ACT_NUC(:,:,:) = XWT(:,:,:) + XWTFRC(:,:,:)
    ELSE
      XWT_ACT_NUC(:,:,:) = XWT(:,:,:)
    END IF
    IF (CTURB /= 'NONE' ) THEN
     IF ( ((CCLOUD=='C2R2'.OR.CCLOUD=='KHKO').AND.LACTTKE) .OR. (CCLOUD=='LIMA'.AND.MACTTKE) ) THEN 
       XWT_ACT_NUC(:,:,:) = XWT_ACT_NUC(:,:,:) +  (2./3. * XTKET(:,:,:))**0.5
     ELSE
       XWT_ACT_NUC(:,:,:) = XWT_ACT_NUC(:,:,:) 
     ENDIF
    ENDIF
  ELSE
    XWT_ACT_NUC(:,:,:) = 0.
  END IF
!
!$acc kernels present_cr(XRTHS_CLD,XRRS_CLD)
  XRTHS_CLD(:, :, : )    = XRTHS(:, :, : )
  XRRS_CLD (:, :, :, : ) = XRRS(:, :, :, : )
!$acc end kernels
  IF ( GNVS_NO_0 ) THEN
!$acc kernels present_cr(XRSVS_CLD)
     XRSVS_CLD(:, :, :, : ) = XRSVS(:, :, :, : )
!$acc end kernels
  END IF

!$acc data present( XRHODJ, XRHODREF, XEXNREF, XRTHS, XRRS, XRSVS, ZPABST, XTHT )                &
!$acc &    copyin (XZZ, XSIGS, VSIGQSAT, XMFCONV, XTHM, XPABSM,                                  &
!$acc &            XRCM, XWT_ACT_NUC, XDTHRAD, XCF_MF, XRC_MF, XRI_MF,                           &
!$acc &            XSOLORG, XMI)                                                                 &
!$acc &    copy   (XSUPSAT, XNACT, XNPRO, XSSPRO,                                                &
!$acc &            XRT, XSVT, XCLDFR, XCIT, XINPRR3D, XEVAP3D,                                   &
!$acc &            XINPRC, XINPRR, XINPRS, XINPRG, XINPRH, XINDEP,                               &
!$acc &            XHLC_HRC, XHLC_HCF, XHLI_HRI, XHLI_HCF)                                       &
!$acc &    copyout(XSRCT, XRAINFR)

  IF (CSURF=='EXTE') THEN
#ifndef MNH_OPENACC
    ALLOCATE (ZSEA(SIZE(XRHODJ,1),SIZE(XRHODJ,2)))
    ALLOCATE (ZTOWN(SIZE(XRHODJ,1),SIZE(XRHODJ,2)))
#else
    CALL MNH_MEM_POSITION_PIN( 'MODEL_n 2')

    CALL MNH_MEM_GET( ZSEA,  SIZE( XRHODJ, 1 ), SIZE( XRHODJ, 2 ) )
    CALL MNH_MEM_GET( ZTOWN, SIZE( XRHODJ, 1 ), SIZE( XRHODJ, 2 ) )
#endif
    CALL MNHGET_SURF_PARAM_n (PSEA=ZSEA(:,:),PTOWN=ZTOWN(:,:))
!$acc update device( ZSEA, ZTOWN )
    CALL RESOLVED_CLOUD ( CCLOUD, CELEC,CACTCCN, CSCONV, CMF_CLOUD, NRR, NSPLITR,   &
                          NSPLITG, IMI, KTCOUNT,                                    &
                          CLBCX,CLBCY,TPBAKFILE, CRAD, CTURBDIM,                    &
                          LSUBG_COND,LSIGMAS,CSUBG_AUCV_RC,XTSTEP,                  &
                          XZZ, XRHODJ, XRHODREF, XEXNREF,                           &
                          ZPABST, XTHT,XRT,XSIGS,VSIGQSAT,XMFCONV,XTHM,XRCM,        &
                          XPABST, XWT_ACT_NUC,XDTHRAD, XRTHS, XRRS,                 &
                          XSVT, XRSVS,                                              &
                          XSRCT, XCLDFR,XICEFR, XCIT,                               &
                          LSEDIC,KACTIT, KSEDC, KSEDI, KRAIN, KWARM, KHHONI,        &
                          LCONVHG, XCF_MF,XRC_MF, XRI_MF, XWEIGHT_MF_CLOUD,         &
                          XHLC_HRC_MF, XHLC_HCF_MF, XHLI_HRI_MF, XHLI_HCF_MF,       &
                          XINPRC,ZINPRC3D,XINPRR, XINPRR3D, XEVAP3D,                &
                          XINPRS,ZINPRS3D, XINPRG,ZINPRG3D, XINPRH,ZINPRH3D,        &
                          XSOLORG, XMI,ZSPEEDC, ZSPEEDR, ZSPEEDS, ZSPEEDG, ZSPEEDH, &
                          XINDEP, XSUPSAT, XNACT, XNPRO,XSSPRO, XRAINFR,            &
                          XHLC_HRC, XHLC_HCF, XHLI_HRI, XHLI_HCF,                   &
                          ZSEA, ZTOWN                                          )
#ifndef MNH_OPENACC
    IF (CELEC == 'NONE') DEALLOCATE(ZTOWN)
    IF (CELEC == 'NONE') DEALLOCATE(ZSEA)
#else
    IF (CELEC == 'NONE') CALL MNH_MEM_RELEASE( 'MODEL_n 2' ) !Release ZSEA and ZTOWN
#endif
  ELSE
    CALL RESOLVED_CLOUD ( CCLOUD, CELEC, CACTCCN, CSCONV, CMF_CLOUD, NRR, NSPLITR,  &
                          NSPLITG, IMI, KTCOUNT,                                    &
                          CLBCX,CLBCY,TPBAKFILE, CRAD, CTURBDIM,                    &
                          LSUBG_COND,LSIGMAS,CSUBG_AUCV_RC,                         &
                          XTSTEP,XZZ, XRHODJ, XRHODREF, XEXNREF,                    &
                          ZPABST, XTHT,XRT,XSIGS,VSIGQSAT,XMFCONV,XTHM,XRCM,        &
                          XPABST, XWT_ACT_NUC,XDTHRAD, XRTHS, XRRS,                 &
                          XSVT, XRSVS,                                              &
                          XSRCT, XCLDFR, XICEFR, XCIT,                              &
                          LSEDIC,KACTIT, KSEDC, KSEDI, KRAIN, KWARM, KHHONI,        &
                          LCONVHG, XCF_MF,XRC_MF, XRI_MF, XWEIGHT_MF_CLOUD,         &
                          XHLC_HRC_MF, XHLC_HCF_MF, XHLI_HRI_MF, XHLI_HCF_MF,       &
                          XINPRC,ZINPRC3D,XINPRR, XINPRR3D, XEVAP3D,                &
                          XINPRS,ZINPRS3D, XINPRG,ZINPRG3D, XINPRH,ZINPRH3D,        &
                          XSOLORG, XMI,ZSPEEDC, ZSPEEDR, ZSPEEDS, ZSPEEDG, ZSPEEDH, &
                          XINDEP, XSUPSAT, XNACT, XNPRO,XSSPRO, XRAINFR,            &
                          XHLC_HRC, XHLC_HCF, XHLI_HRI, XHLI_HCF                    )
  END IF
!$acc end data

!$acc kernels present_cr(XRTHS_CLD,XRRS_CLD)
  XRTHS_CLD(:, :, : )    = XRTHS(:, :, : )    - XRTHS_CLD(:, :, : )
  XRRS_CLD (:, :, :, : ) = XRRS (:, :, :, : ) - XRRS_CLD (:, :, :, : )
!$acc end kernels
  IF ( GNVS_NO_0 ) THEN
!$acc kernels present_cr(XRSVS_CLD)  
     XRSVS_CLD(:, :, :, : ) = XRSVS(:, :, :, : ) - XRSVS_CLD(:, :, :, : )
!$acc end kernels
  END IF
!
!$acc update host( XRTHS_CLD, XRRS_CLD, XRSVS_CLD )
!
  IF (CCLOUD /= 'REVE' ) THEN
    XACPRR = XACPRR + XINPRR * XTSTEP
          IF ( (CCLOUD(1:3) == 'ICE' .AND. LSEDIC ) .OR.                     &
        ((CCLOUD == 'C2R2' .OR. CCLOUD == 'C3R5' .OR. CCLOUD == 'KHKO' &
                           .OR. CCLOUD == 'LIMA' ) .AND. KSEDC ) )    THEN
      XACPRC = XACPRC + XINPRC * XTSTEP
      IF (LDEPOSC .OR. LDEPOC) XACDEP = XACDEP + XINDEP * XTSTEP
    END IF
    IF (CCLOUD(1:3) == 'ICE' .OR. CCLOUD == 'C3R5' .OR. &
                                 (CCLOUD == 'LIMA' .AND. NMOM_I.GE.1 ) ) THEN
      XACPRS = XACPRS + XINPRS * XTSTEP
      XACPRG = XACPRG + XINPRG * XTSTEP
      IF (CCLOUD == 'ICE4' .OR. (CCLOUD == 'LIMA' .AND. NMOM_H.GE.1)) XACPRH = XACPRH + XINPRH * XTSTEP          
    END IF
!
! Lessivage des CCN et IFN nucléables par Slinn
!
    IF (LSCAV .AND. (CCLOUD == 'LIMA')) THEN
       CALL FILL_DIMPHYEX( YLDIMPHYEX, SIZE(XTHT,1), SIZE(XTHT,2), SIZE(XTHT,3) )
       CALL LIMA_PRECIP_SCAVENGING( TNSV, YLDIMPHYEX,CST,TBUCONF,TBUDGETS_PTR,SIZE(TBUDGETS_PTR), &
                                    CCLOUD, CCONF, ILUOUT, KTCOUNT,XTSTEP,XRT(:,:,:,3),    &
                                    XRHODREF, XRHODJ, XZZ, XPABST, XTHT,            &
                                    XSVT(:,:,:,NSV_LIMA_BEG:NSV_LIMA_END),          &
                                    XRSVS(:,:,:,NSV_LIMA_BEG:NSV_LIMA_END), XINPAP  )
!
      XACPAP(:,:) = XACPAP(:,:) + XINPAP(:,:) * XTSTEP
    END IF
  END IF
!
! It is necessary that SV_C2R2 and SV_C1R3 are contiguous in the preceeding CALL
!
END IF
!$acc update host( XRTHS, XRRS, XRSVS )
!
CALL SECOND_MNH2(ZTIME2)
!
XT_CLOUD = XT_CLOUD + ZTIME2 - ZTIME1 &
           - XTIME_LES_BU_PROCESS - XTIME_BU_PROCESS
!
!-------------------------------------------------------------------------------
!
!*       21.    CLOUD ELECTRIFICATION AND LIGHTNING FLASHES
!               -------------------------------------------
!
! Cloud electrification is now called directly from resolved_cloud
! It avoids duplicating microphysics routines.
! Resolved_elec solves the ion recombination and attachement, and
! lightning flash triggering and propagation
!
ZTIME1 = ZTIME2
XTIME_BU_PROCESS = 0.
XTIME_LES_BU_PROCESS = 0.
!
IF (CELEC /= 'NONE') THEN  !++cb-- ATTENTION : le cas rain_ice_elec n'est pas traite !!!
  IF (CCLOUD(1:3) == 'ICE') THEN
    IF (CSURF == 'EXTE') THEN
      IF (LLNOX_EXPLICIT) THEN
        CALL RESOLVED_ELEC_n (CCLOUD, NRR, IMI, KTCOUNT, OEXIT,       &
                              XTSTEP, XZZ, XRHODJ, XRHODREF, XEXNREF, &
                              ZPABST, XTHT, XWT, XRT, XRRS,           &
                              XSVT(:,:,:,NSV_ELECBEG:NSV_ELECEND),    &
                              XRSVS(:,:,:,NSV_ELECBEG:NSV_ELECEND),   &
                              XCIT, XINPRR,                           &
                              PSEA=ZSEA, PTOWN=ZTOWN,                 &
                              PSVS_LNOX=XRSVS(:,:,:,NSV_LNOXBEG)      )
      ELSE
        CALL RESOLVED_ELEC_n (CCLOUD, NRR, IMI, KTCOUNT, OEXIT,       &
                              XTSTEP, XZZ, XRHODJ, XRHODREF, XEXNREF, &
                              ZPABST, XTHT, XWT, XRT, XRRS,           &
                              XSVT(:,:,:,NSV_ELECBEG:NSV_ELECEND),    &
                              XRSVS(:,:,:,NSV_ELECBEG:NSV_ELECEND),   &
                              XCIT, XINPRR,                           &
                              PSEA=ZSEA, PTOWN=ZTOWN                  )
      END IF
    ELSE
      IF (LLNOX_EXPLICIT) THEN
        CALL RESOLVED_ELEC_n (CCLOUD, NRR, IMI, KTCOUNT, OEXIT,       &
                              XTSTEP, XZZ, XRHODJ, XRHODREF, XEXNREF, &
                              ZPABST, XTHT, XWT, XRT, XRRS,           &
                              XSVT(:,:,:,NSV_ELECBEG:NSV_ELECEND),    &
                              XRSVS(:,:,:,NSV_ELECBEG:NSV_ELECEND),   &
                              XCIT, XINPRR,                           &
                              PSVS_LNOX=XRSVS(:,:,:,NSV_LNOXBEG)      )
      ELSE
        CALL RESOLVED_ELEC_n (CCLOUD, NRR, IMI, KTCOUNT, OEXIT,       &
                              XTSTEP, XZZ, XRHODJ, XRHODREF, XEXNREF, &
                              ZPABST, XTHT, XWT, XRT, XRRS,           &
                              XSVT(:,:,:,NSV_ELECBEG:NSV_ELECEND),    &
                              XRSVS(:,:,:,NSV_ELECBEG:NSV_ELECEND),   &
                              XCIT, XINPRR                            )
      END IF
    END IF
  ELSE IF (CCLOUD == 'LIMA' .AND. LPTSPLIT) THEN
    IF (LLNOX_EXPLICIT) THEN
      IF ((NRR == 6 .AND. NMOM_S == 1 .AND. NMOM_G == 1) .OR. &
          (NRR == 7 .AND. NMOM_S == 1 .AND. NMOM_G == 1 .AND. NMOM_H == 1)) THEN
        CALL RESOLVED_ELEC_n (CCLOUD, NRR, IMI, KTCOUNT, OEXIT,       &
                              XTSTEP, XZZ, XRHODJ, XRHODREF, XEXNREF, &
                              ZPABST, XTHT, XWT, XRT, XRRS,           &
                              XSVT(:,:,:,NSV_ELECBEG:NSV_ELECEND),    &
                              XRSVS(:,:,:,NSV_ELECBEG:NSV_ELECEND),   &
                              XRSVS(:,:,:,NSV_LIMA_NI), XINPRR,       &
                              PCCS=XRSVS(:,:,:,NSV_LIMA_NC),          &
                              PCRS=XRSVS(:,:,:,NSV_LIMA_NR),          &
                              PSVS_LNOX=XRSVS(:,:,:,NSV_LNOXBEG)      )
      ELSE IF (NRR == 6 .AND. NMOM_S == 2 .AND. NMOM_G == 2) THEN
        CALL RESOLVED_ELEC_n (CCLOUD, NRR, IMI, KTCOUNT, OEXIT,       &
                              XTSTEP, XZZ, XRHODJ, XRHODREF, XEXNREF, &
                              ZPABST, XTHT, XWT, XRT, XRRS,           &
                              XSVT(:,:,:,NSV_ELECBEG:NSV_ELECEND),    &
                              XRSVS(:,:,:,NSV_ELECBEG:NSV_ELECEND),   &
                              XRSVS(:,:,:,NSV_LIMA_NI), XINPRR,       &
                              PCCS=XRSVS(:,:,:,NSV_LIMA_NC),          &
                              PCRS=XRSVS(:,:,:,NSV_LIMA_NR),          &
                              PCSS=XRSVS(:,:,:,NSV_LIMA_NS),          &
                              PCGS=XRSVS(:,:,:,NSV_LIMA_NG),          &
                              PSVS_LNOX=XRSVS(:,:,:,NSV_LNOXBEG)      )
      ELSE IF (NRR == 7 .AND. NMOM_S == 2 .AND. NMOM_G == 2 .AND. NMOM_H == 2) THEN
        CALL RESOLVED_ELEC_n (CCLOUD, NRR, IMI, KTCOUNT, OEXIT,       &
                              XTSTEP, XZZ, XRHODJ, XRHODREF, XEXNREF, &
                              ZPABST, XTHT, XWT, XRT, XRRS,           &
                              XSVT(:,:,:,NSV_ELECBEG:NSV_ELECEND),    &
                              XRSVS(:,:,:,NSV_ELECBEG:NSV_ELECEND),   &
                              XRSVS(:,:,:,NSV_LIMA_NI), XINPRR,       &
                              PCCS=XRSVS(:,:,:,NSV_LIMA_NC),          &
                              PCRS=XRSVS(:,:,:,NSV_LIMA_NR),          &
                              PCSS=XRSVS(:,:,:,NSV_LIMA_NS),          &
                              PCGS=XRSVS(:,:,:,NSV_LIMA_NG),          &
                              PCHS=XRSVS(:,:,:,NSV_LIMA_NH),          &
                              PSVS_LNOX=XRSVS(:,:,:,NSV_LNOXBEG)      )
      END IF
    ELSE
      IF ((NRR == 6 .AND. NMOM_S == 1 .AND. NMOM_G == 1) .OR. &
          (NRR == 7 .AND. NMOM_S == 1 .AND. NMOM_G == 1 .AND. NMOM_H == 1)) THEN
        CALL RESOLVED_ELEC_n (CCLOUD, NRR, IMI, KTCOUNT, OEXIT,       &
                              XTSTEP, XZZ, XRHODJ, XRHODREF, XEXNREF, &
                              ZPABST, XTHT, XWT, XRT, XRRS,           &
                              XSVT(:,:,:,NSV_ELECBEG:NSV_ELECEND),    &
                              XRSVS(:,:,:,NSV_ELECBEG:NSV_ELECEND),   &
                              XRSVS(:,:,:,NSV_LIMA_NI), XINPRR,       &
                              PCCS=XRSVS(:,:,:,NSV_LIMA_NC),          &
                              PCRS=XRSVS(:,:,:,NSV_LIMA_NR))
      ELSE IF (NRR == 6 .AND. NMOM_S == 2 .AND. NMOM_G == 2) THEN
        CALL RESOLVED_ELEC_n (CCLOUD, NRR, IMI, KTCOUNT, OEXIT,       &
                              XTSTEP, XZZ, XRHODJ, XRHODREF, XEXNREF, &
                              ZPABST, XTHT, XWT, XRT, XRRS,           &
                              XSVT(:,:,:,NSV_ELECBEG:NSV_ELECEND),    &
                              XRSVS(:,:,:,NSV_ELECBEG:NSV_ELECEND),   &
                              XRSVS(:,:,:,NSV_LIMA_NI), XINPRR,       &
                              PCCS=XRSVS(:,:,:,NSV_LIMA_NC),          &
                              PCRS=XRSVS(:,:,:,NSV_LIMA_NR),          &
                              PCSS=XRSVS(:,:,:,NSV_LIMA_NS),          &
                              PCGS=XRSVS(:,:,:,NSV_LIMA_NG))
      ELSE IF (NRR == 7 .AND. NMOM_S == 2 .AND. NMOM_G == 2 .AND. NMOM_H == 2) THEN
        CALL RESOLVED_ELEC_n (CCLOUD, NRR, IMI, KTCOUNT, OEXIT,       &
                              XTSTEP, XZZ, XRHODJ, XRHODREF, XEXNREF, &
                              ZPABST, XTHT, XWT, XRT, XRRS,           &
                              XSVT(:,:,:,NSV_ELECBEG:NSV_ELECEND),    &
                              XRSVS(:,:,:,NSV_ELECBEG:NSV_ELECEND),   &
                              XRSVS(:,:,:,NSV_LIMA_NI), XINPRR,       &
                              PCCS=XRSVS(:,:,:,NSV_LIMA_NC),          &
                              PCRS=XRSVS(:,:,:,NSV_LIMA_NR),          &
                              PCSS=XRSVS(:,:,:,NSV_LIMA_NS),          &
                              PCGS=XRSVS(:,:,:,NSV_LIMA_NG),          &
                              PCHS=XRSVS(:,:,:,NSV_LIMA_NH))
      END IF
    END IF
  END IF
  !
  IF ( CSURF == 'EXTE' ) THEN
#ifndef MNH_OPENACC
    DEALLOCATE(ZTOWN)
    DEALLOCATE(ZSEA)
#else
    CALL MNH_MEM_RELEASE( 'MODEL_n 2' ) !Release ZSEA and ZTOWN
#endif
  END IF
END IF
!
CALL SECOND_MNH2(ZTIME2)
!
XT_ELEC = XT_ELEC + ZTIME2 - ZTIME1 &
           - XTIME_LES_BU_PROCESS - XTIME_BU_PROCESS
!
!-------------------------------------------------------------------------------
!
!*       21.    L.E.S. COMPUTATIONS
!               -------------------
!
ZTIME1 = ZTIME2
!
CALL LES_n
!
CALL SECOND_MNH2(ZTIME2)
!
XT_SPECTRA = XT_SPECTRA + ZTIME2 - ZTIME1 + XTIME_LES_BU + XTIME_LES
!
!-------------------------------------------------------------------------------
!
!*       21. bis    MEAN_UM
!               --------------------
!
IF (LMEAN_FIELD) THEN
   CALL MEAN_FIELD(XUT, XVT, XWT, XTHT, XTKET, XPABST, XRT(:,:,:,1), XSVT(:,:,:,1))
END IF
!
!-------------------------------------------------------------------------------
!
!*       22.    UPDATE HALO OF EACH SUBDOMAINS FOR TIME T+DT
!               --------------------------------------------
!
ZTIME1 = ZTIME2
!
CALL EXCHANGE (XTSTEP,NRR,NSV,XRHODJ,TFIELDS_ll,     &
               XRUS, XRVS,XRWS,XRTHS,XRRS,XRTKES,XRSVS)
!
CALL SECOND_MNH2(ZTIME2)
!
XT_HALO = XT_HALO + ZTIME2 - ZTIME1
!
!-------------------------------------------------------------------------------
!
!*       23.    TEMPORAL SWAPPING
!               -----------------
!
ZTIME1 = ZTIME2
XTIME_BU_PROCESS = 0.
!
CALL ENDSTEP  ( XTSTEP,NRR,NSV,KTCOUNT,IMI,               &
                CUVW_ADV_SCHEME,CTEMP_SCHEME,XRHODJ,      &
                XRUS,XRVS,XRWS,XDRYMASSS,                 &
                XRTHS,XRRS,XRTKES,XRSVS,                  &
                XLSUS,XLSVS,XLSWS,                        &
                XLSTHS,XLSRVS,XLSZWSS,                    &
                XLBXUS,XLBXVS,XLBXWS,                     &
                XLBXTHS,XLBXRS,XLBXTKES,XLBXSVS,          &
                XLBYUS,XLBYVS,XLBYWS,                     &
                XLBYTHS,XLBYRS,XLBYTKES,XLBYSVS,          &
                XUM,XVM,XWM,XZWS,                         &
                XUT,XVT,XWT,XPABST,XDRYMASST,             &
                XTHT, XRT, XTHM, XRCM, XPABSM,XTKET, XSVT,&
                XLSUM,XLSVM,XLSWM,                        &
                XLSTHM,XLSRVM,XLSZWSM,                    &
                XLBXUM,XLBXVM,XLBXWM,                     &
                XLBXTHM,XLBXRM,XLBXTKEM,XLBXSVM,          &
                XLBYUM,XLBYVM,XLBYWM,                     &
                XLBYTHM,XLBYRM,XLBYTKEM,XLBYSVM           )
!
CALL SECOND_MNH2(ZTIME2)
!
XT_STEP_SWA = XT_STEP_SWA + ZTIME2 - ZTIME1 - XTIME_BU_PROCESS
!
!-------------------------------------------------------------------------------
!
!*       24.1    BALLOON and AIRCRAFT
!               --------------------
!
ZTIME1 = ZTIME2
!
IF (LFLYER) THEN
  IF (CSURF=='EXTE') THEN
#ifndef MNH_OPENACC
    ALLOCATE(ZSEA(IIU,IJU))
#else
    CALL MNH_MEM_POSITION_PIN( 'MODEL_n 4' )

    CALL MNH_MEM_GET( ZSEA, IIU, IJU )
#endif
    ZSEA(:,:) = 0.
    CALL MNHGET_SURF_PARAM_n (PSEA=ZSEA(:,:))
    CALL AIRCRAFT_BALLOON( XTSTEP, XZZ, XMAP, XLONORI, XLATORI,                   &
                           XUT, XVT, XWT, XPABST, XTHT, XRT, XSVT, XTKET, XTSRAD, &
                           XRHODREF, XCIT, PSEA = ZSEA(:,:)                       )
#ifndef MNH_OPENACC
    DEALLOCATE(ZSEA)
#else
    CALL MNH_MEM_RELEASE( 'MODEL_n 4' ) !Release ZSEA
#endif
  ELSE
    CALL AIRCRAFT_BALLOON( XTSTEP, XZZ, XMAP, XLONORI, XLATORI,                   &
                           XUT, XVT, XWT, XPABST, XTHT, XRT, XSVT, XTKET, XTSRAD, &
                           XRHODREF, XCIT                                         )
  END IF
END IF

!-------------------------------------------------------------------------------
!
!*       24.2    STATION (observation diagnostic)
!               --------------------------------
!
IF ( LSTATION ) &
  CALL STATION_n( XZZ, XRHODREF, XUT, XVT, XWT, XTHT, XRT, XSVT, XTKET, XTSRAD, XPABST )
!
!---------------------------------------------------------
!
!*       24.3    PROFILER (observation diagnostic)
!               ---------------------------------
!
IF (LPROFILER)  THEN
  IF (CSURF=='EXTE') THEN
#ifndef MNH_OPENACC
    ALLOCATE(ZSEA(IIU,IJU))
#else
    CALL MNH_MEM_POSITION_PIN( 'MODEL_n 5' )

    CALL MNH_MEM_GET( ZSEA, IIU, IJU )
#endif
    ZSEA(:,:) = 0.
    CALL MNHGET_SURF_PARAM_n (PSEA=ZSEA(:,:))
    CALL PROFILER_n( XZZ, XRHODREF,                             &
                     XUT, XVT, XWT, XTHT, XRT, XSVT, XTKET,     &
                     XTSRAD, XPABST, XAER, XCIT, PSEA=ZSEA(:,:) )
#ifndef MNH_OPENACC
    DEALLOCATE(ZSEA)
#else
    CALL MNH_MEM_RELEASE( 'MODEL_n 5' ) !Release ZSEA
#endif
  ELSE
    CALL PROFILER_n( XZZ, XRHODREF,                         &
                     XUT, XVT, XWT, XTHT, XRT, XSVT, XTKET, &
                     XTSRAD, XPABST, XAER, XCIT             )
  END IF
END IF
!
!
CALL SECOND_MNH2(ZTIME2)
!
XT_STEP_MISC = XT_STEP_MISC + ZTIME2 - ZTIME1
!
!-------------------------------------------------------------------------------
!
!*       24.4   deallocation of observation diagnostics
!               ---------------------------------------
!
CALL END_DIAG_IN_RUN
!
!-------------------------------------------------------------------------------
!
!
!*       25.    STORAGE OF BUDGET FIELDS
!               ------------------------
!
ZTIME1 = ZTIME2
!
IF ( .NOT. LIO_NO_WRITE ) THEN
  IF (NBUMOD==IMI .AND. CBUTYPE/='NONE') THEN
    CALL ENDSTEP_BUDGET(TDIAFILE,KTCOUNT,TDTCUR,XTSTEP,NSV)
  END IF
END IF
!
CALL SECOND_MNH2(ZTIME2)
!
XT_STEP_BUD = XT_STEP_BUD + ZTIME2 - ZTIME1 + XTIME_BU
!
!-------------------------------------------------------------------------------
!
!*       27.    CURRENT TIME REFRESH
!               --------------------
!
TDTCUR = TDTCUR + XTSTEP
!
!-------------------------------------------------------------------------------
!
!*       28.    CPU ANALYSIS
!               ------------
!
CALL SECOND_MNH2(ZTIME2)
XT_START=XT_START+ZTIME2-ZEND
!
!
IF ( KTCOUNT == NSTOP .AND. IMI==1) THEN
  OEXIT=.TRUE.
END IF
!
IF (OEXIT) THEN
!
  IF ( .NOT. LIO_NO_WRITE ) THEN
    IF ( LSERIES )   CALL WRITE_SERIES_n(TDIAFILE)
    IF ( LFLYER )    CALL WRITE_AIRCRAFT_BALLOON(TDIAFILE)
    IF ( LSTATION )  CALL WRITE_STATPROF_n( TDIAFILE, TSTATIONS )
    IF ( LPROFILER ) CALL WRITE_STATPROF_n( TDIAFILE, TPROFILERS )
    call Write_les_n( tdiafile )
    CALL IO_File_close(TDIAFILE)
    ! Free memory of flyer that is not present on the master process of the file (was allocated in WRITE_AIRCRAFT_BALLOON)
    IF ( LFLYER ) CALL AIRCRAFT_BALLOON_FREE_NONLOCAL( TDIAFILE )
  END IF
  !
  CALL IO_File_close(TINIFILE)
  IF ( CSURF == "EXTE" .AND. ASSOCIATED( TINIFILEPGD ) ) CALL IO_File_close( TINIFILEPGD )
!
!*       28.1   print statistics!
!
  ! Set File Timing OUTPUT
  !
  CALL SET_ILUOUT_TIMING(TLUOUT)
  !
  ! Compute global time
  !
  CALL TIME_STAT_ll(XT_START,ZTOT)
  !
  CALL TIME_HEADER_ll(IMI)
  !
  CALL TIME_STAT_ll(XT_1WAY,ZTOT,       ' ONE WAY','=')
  CALL TIME_STAT_ll(XT_BOUND,ZTOT,      ' BOUNDARIES','=')
  CALL TIME_STAT_ll(XT_STORE,ZTOT,      ' STORE-FIELDS','=')
    CALL TIME_STAT_ll(TIMEZ%T_WRIT3D_SEND,ZTOT,    '   W3D_SEND ','-')
    CALL TIME_STAT_ll(TIMEZ%T_WRIT3D_RECV,ZTOT,    '   W3D_RECV ','-')
    CALL TIME_STAT_ll(TIMEZ%T_WRIT3D_WRIT,ZTOT,    '   W3D_WRIT ','-')
    CALL TIME_STAT_ll(TIMEZ%T_WRIT3D_WAIT,ZTOT,    '   W3D_WAIT ','-')
    CALL TIME_STAT_ll(TIMEZ%T_WRIT3D_ALL ,ZTOT,    '   W3D_ALL ','-')
    CALL TIME_STAT_ll(TIMEZ%T_WRIT2D_GATH,ZTOT,    '   W2D_GATH ','-')
    CALL TIME_STAT_ll(TIMEZ%T_WRIT2D_WRIT,ZTOT,    '   W2D_WRIT ','-')
    CALL TIME_STAT_ll(TIMEZ%T_WRIT2D_ALL ,ZTOT,    '   W2D_ALL ','-')
  CALL TIME_STAT_ll(XT_GUESS,ZTOT,      ' INITIAL_GUESS','=')
  CALL TIME_STAT_ll(XT_2WAY,ZTOT,       ' TWO WAY','=')
  CALL TIME_STAT_ll(XT_ADV,ZTOT,        ' ADVECTION MET','=')
  CALL TIME_STAT_ll(XT_ADVUVW,ZTOT,     ' ADVECTION UVW','=')
  CALL TIME_STAT_ll(XT_GRAV,ZTOT,       ' GRAVITY','=')
  CALL TIME_STAT_ll(XT_FORCING,ZTOT,    ' FORCING','=')
  CALL TIME_STAT_ll(XT_IBM_FORC,ZTOT,   ' IBM','=')  
  CALL TIME_STAT_ll(XT_NUDGING,ZTOT,    ' NUDGING','=')
  CALL TIME_STAT_ll(XT_SOURCES,ZTOT,    ' DYN_SOURCES','=')
  CALL TIME_STAT_ll(XT_DIFF,ZTOT,       ' NUM_DIFF','=')
  CALL TIME_STAT_ll(XT_RELAX,ZTOT,      ' RELAXATION','=')
  !
  CALL  TIMING_LEGEND() 
  !
  CALL TIME_STAT_ll(XT_PARAM,ZTOT,      ' PHYS_PARAM','=')
    CALL TIME_STAT_ll(XT_RAD,ZTOT,      '   RAD       = '//CRAD  ,'-')
    CALL TIME_STAT_ll(XT_SHADOWS,ZTOT,  '   SHADOWS'             ,'-')
    CALL TIME_STAT_ll(XT_DCONV,ZTOT,    '   DEEP CONV = '//CDCONV,'-')
    CALL TIME_STAT_ll(XT_GROUND,ZTOT,   '   GROUND'              ,'-')
    ! Blaze perf
    IF (LBLAZE) THEN
      CALL TIME_STAT_ll(XFIREPERF,ZBLAZETOT)
      CALL TIME_STAT_ll(XFIREPERF,ZTOT,           '     BLAZE'          ,'~')
      CALL TIME_STAT_ll(XGRADPERF,ZBLAZETOT,      '       GRAD(PHI)'    ,' ')
      CALL TIME_STAT_ll(XROSWINDPERF,ZBLAZETOT,   '       ROS & WIND'   ,' ')
      CALL TIME_STAT_ll(XPROPAGPERF,ZBLAZETOT,    '       PROPAGATION'  ,' ')
      CALL TIME_STAT_ll(XFLUXPERF,ZBLAZETOT,      '       HEAT FLUXES'  ,' ')
    END IF
    CALL TIME_STAT_ll(XT_TURB,ZTOT,     '   TURB      = '//CTURB ,'-')
    CALL TIME_STAT_ll(XT_MAFL,ZTOT,     '   MAFL      = '//CSCONV,'-')
    CALL TIME_STAT_ll(XT_CHEM,ZTOT,     '   CHIMIE'              ,'-')
    CALL TIME_STAT_ll(XT_EOL,ZTOT,      '   WIND TURBINE'        ,'-')
  CALL  TIMING_LEGEND()
  CALL TIME_STAT_ll(XT_COUPL,ZTOT,      ' SET_COUPLING','=')
  CALL TIME_STAT_ll(XT_RAD_BOUND,ZTOT,  ' RAD_BOUND','=')
  !
  CALL  TIMING_LEGEND()
  ! 
  CALL TIME_STAT_ll( XT_PRESS, ZTOT,      ' PRESSURE ', '=', GFULLSTAT_PRESS_SLV )
  !JUAN Z_SPLITTING
    CALL TIME_STAT_ll( TIMEZ%T_MAP_B_SX_YP2_ZP1,          ZTOT, '   REMAP       B=>FFTXZ'  , '-', GFULLSTAT_PRESS_SLV )
    CALL TIME_STAT_ll( TIMEZ%T_MAP_SX_YP2_ZP1_SXP2_Y_ZP1, ZTOT, '   REMAP   FFTXZ=>FFTYZ'  , '-', GFULLSTAT_PRESS_SLV )
    CALL TIME_STAT_ll( TIMEZ%T_MAP_SXP2_Y_ZP1_B,          ZTOT, '   REMAP   FTTYZ=>B'      , '-', GFULLSTAT_PRESS_SLV )
    CALL TIME_STAT_ll( TIMEZ%T_MAP_SXP2_Y_ZP1_SXP2_YP1_Z, ZTOT, '   REMAP   FFTYZ=>SUBZ'   , '-', GFULLSTAT_PRESS_SLV )
    CALL TIME_STAT_ll( TIMEZ%T_MAP_B_SXP2_Y_ZP1,          ZTOT, '   REMAP       B=>FFTYZ-1', '-', GFULLSTAT_PRESS_SLV )
    CALL TIME_STAT_ll( TIMEZ%T_MAP_SXP2_YP1_Z_SXP2_Y_ZP1, ZTOT, '   REMAP    SUBZ=>FFTYZ-1', '-', GFULLSTAT_PRESS_SLV )
    CALL TIME_STAT_ll( TIMEZ%T_MAP_SXP2_Y_ZP1_SX_YP2_ZP1, ZTOT, '   REMAP FFTYZ-1=>FFTXZ-1', '-', GFULLSTAT_PRESS_SLV )
    CALL TIME_STAT_ll( TIMEZ%T_MAP_SX_YP2_ZP1_B,          ZTOT, '   REMAP FFTXZ-1=>B     ' , '-', GFULLSTAT_PRESS_SLV )
  ! JUAN P1/P2
  CALL TIME_STAT_ll(XT_CLOUD,ZTOT,      ' RESOLVED_CLOUD','=')
  CALL TIME_STAT_ll(XT_ELEC,ZTOT,      ' RESOLVED_ELEC','=')
  CALL TIME_STAT_ll(XT_HALO,ZTOT,       ' EXCHANGE_HALO','=')
  CALL TIME_STAT_ll(XT_STEP_SWA,ZTOT,   ' ENDSTEP','=')
  CALL TIME_STAT_ll(XT_STEP_BUD,ZTOT,   ' BUDGETS','=')
  CALL TIME_STAT_ll(XT_SPECTRA,ZTOT,    ' LES','=')
  CALL TIME_STAT_ll(XT_STEP_MISC,ZTOT,  ' MISCELLANEOUS','=')
  IF (LIBM) CALL TIME_STAT_ll(XT_IBM_FORC,ZTOT,' IBM FORCING','=') 
  !
  ! sum of call subroutine
  !
  ZALL   = XT_1WAY + XT_BOUND   + XT_STORE   + XT_GUESS    +  XT_2WAY   + &
           XT_ADV  + XT_FORCING + XT_NUDGING + XT_SOURCES  +  XT_DIFF   + &
           XT_ADVUVW  + XT_GRAV + XT_IBM_FORC                           + &
           XT_RELAX+ XT_PARAM   + XT_COUPL   + XT_RAD_BOUND+XT_PRESS    + &
           XT_CLOUD+ XT_ELEC    + XT_HALO    + XT_SPECTRA + XT_STEP_SWA + &
           XT_STEP_MISC+ XT_STEP_BUD
  CALL TIME_STAT_ll(ZALL,ZTOT,          ' SUM(CALL)','=')
  CALL  TIMING_SEPARATOR('=')
  !
  ! Gobale Stat
  !
  WRITE(ILUOUT,FMT=*)
  WRITE(ILUOUT,FMT=*)
  CALL  TIMING_LEGEND() 
  !
  ! MODELN all included
  !
  CALL  TIMING_SEPARATOR('+')
  CALL  TIMING_SEPARATOR('+')  
  WRITE(YMI,FMT="(I0)") IMI
  CALL TIME_STAT_ll(XT_START,ZTOT,      ' MODEL'//YMI,'+')
  CALL  TIMING_SEPARATOR('+')
  CALL  TIMING_SEPARATOR('+')
  CALL  TIMING_SEPARATOR('+')
  !
  ! Timing/ Steps
  !
  ZTIME_STEP     =  XT_START / REAL(KTCOUNT)
  WRITE(YTCOUNT,FMT="(I0)") KTCOUNT
  CALL TIME_STAT_ll(ZTIME_STEP,ZTOT,     ' SECOND/STEP='//YTCOUNT,'=')
  !
  ! Timing/Step/Points
  !
  IPOINTS = NIMAX_ll*NJMAX_ll*NKMAX
  WRITE(YPOINTS,FMT="(I0)") IPOINTS
  ZTIME_STEP_PTS =  ZTIME_STEP / REAL(IPOINTS) * 1e6
  CALL TIME_STAT_ll(ZTIME_STEP_PTS,ZTOT_PT)
  CALL TIME_STAT_ll(ZTIME_STEP_PTS,ZTOT_PT,  ' MICROSEC/STP/PT='//YPOINTS,'-')
  !
  CALL  TIMING_SEPARATOR('=')
  !
END IF

#ifndef MNH_OPENACC
DEALLOCATE( ZRUS,ZRVS,ZRWS )
DEALLOCATE( ZPABST )
#else
CALL MNH_MEM_RELEASE( 'MODEL_n 1' )
#endif

END SUBROUTINE MODEL_n
