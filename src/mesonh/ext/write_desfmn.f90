!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     #########################
      MODULE MODI_WRITE_DESFM_n
!     #########################
!
INTERFACE
!
SUBROUTINE WRITE_DESFM_n(KMI,TPDATAFILE)
!
USE MODD_IO, ONLY: TFILEDATA
!
INTEGER,            INTENT(IN)  :: KMI        ! Model index
TYPE(TFILEDATA),    INTENT(IN)  :: TPDATAFILE ! Datafile
!
END SUBROUTINE WRITE_DESFM_n
!
END INTERFACE
!
END MODULE MODI_WRITE_DESFM_n
!
!
!     ###################################################
      SUBROUTINE WRITE_DESFM_n(KMI,TPDATAFILE)
!     ###################################################
!
!!****  *WRITE_DESFM_n * - routine to write a descriptor file ( DESFM )
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to write the descriptive part of a Mesonh
!     file (FM-file). The resulting file is called  DESFM. 
!       
!!
!!**  METHOD
!!    ------
!!     
!!     This routine writes in the file HDESFM, previously opened, the group of
!!     all the namelists used to specify a Mesonh simulation.
!!     If verbose option is high enough : NVERB>=5, the variables in descriptor
!!     file are printed on the right output-listing corresponding tomodel _n.
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      Module MODN_LUNIT_n : contains declarations of namelist NAM_LUNITn
!!                           and module MODD_LUNIT_n 
!!
!!
!!      Module MODN_CONF_n : contains declaration of namelist NAM_CONFn and
!!                          uses module MODD_CONF1 (configuration variables
!!                             for  model _n )
!!
!!      Module MODN_DYN_n : contains declaration of namelist NAM_DYNn and
!!                         uses module MODD_DYN_n (dynamic control variables
!!                             for   model _n )
!!
!!      Module MODN_ADV_n : contains declaration of namelist NAM_ADVn and
!!                         uses module MODD_ADV_n (control variables for the
!!                         advection scheme for model _n )
!!
!!      Module MODN_PARAM_n : contains declaration of namelist NAM_PARAMn and
!!                         uses module MODD_PARAM_n (names of the physical
!!                         parameterizations for model _n )
!!
!!      Module MODN_PARAM_RAD_n : contains declaration of the control parameters
!!                                for calling the radiation scheme
!!
!!      Module MODN_PARAM_KAFR_n : contains declaration of control parameters
!!                                    for calling the deep convection scheme
!! 
!!      Module MODN_LBC_n : contains declaration of namelis NAM_LBCn and
!!                         uses module MODD_LBC_n (lateral boundary conditions)
!!
!!
!!      Module MODN_TURB_n : contains declaration of turbulence scheme options
!!                        present in the namelist
!!
!!      Module MODN_CONF : contains declaration of namelist NAM_CONF and 
!!                        uses module MODD_CONF (configuration variables) 
!!
!!      Module MODN_DYN : contains the  declaration of namelist NAM_DYN and
!!                        uses module MODD_DYN (dynamic control variables)
!!                         
!!      Module MODN_BUDGET : contains declaration of all the namelists
!!                        related to the budget computations
!!
!!      Module MODN_LES  : contains declaration of the control parameters
!!                                for Large Eddy Simulations' storages
!!      Module MODN_BLANK_n : contains declaration of MesoNH developper variables
!!                                for test and debugging purposes.
!!           
!!
!!    REFERENCE
!!    ---------
!!      None
!!      
!!
!!    AUTHOR
!!    ------
!!	V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original                      07/06/94 
!!      Updated      V.Ducrocq        06/09/94
!!      Updated      J.Stein          20/10/94 to include NAM_OUTn
!!      Updated      J.Stein          24/10/94 change routine name
!!      Updated      J.Stein          26/10/94 add the OWRIGET argument
!!      Updated      J.Stein          06/12/94 add the LS fields       
!!      Updated      J.Stein          09/01/95 add the turbulence scheme     
!!      Updated      J.Stein          09/01/95 add the 1D switch
!!      Updated      J.Stein          20/03/95 remove R from the historical var.
!!      Updated      Ph.Hereil        20/06/95 add the budgets
!!      Updated      J.-P. Pinty      15/09/95 add the radiations
!!      Updated      J.Vila           06/02/96 implementation of scalar
!!                                             advection schemes
!!      Updated      J.Stein          20/02/96 cleaning + add the LES namelist
!!      Modifications  25/04/96  (Suhre)  add NAM_BLANK
!!      Modifications  25/04/96  (Suhre)  add NAM_FRC
!!      Modifications  25/04/96  (Suhre)  add NAM_CH_MNHCn and NAM_CH_SOLVER
!!      Modifications  11/04/96  (Pinty)  add the ice concentration
!!      Modifications  11/01/97  (Pinty)  add the deep convection
!!      Temporary Modification (Masson 06/09/96) manual write of the first and
!!                              third namelists because of compiler version.
!!      Modifications  J.-P. Lafore   22/07/96 gridnesting implementation
!!      Modifications  J.-P. Lafore   29/07/96 add NAM_FMOUT (renamed in NAM_OUTPUT/NAM_BACKUP)
!!      Modifications  V. Masson      10/07/97 add NAM_PARAM_GROUNDn
!!      Modifications  V. Masson      28/07/97 supress LSTEADY_DMASS
!!      Modifications  P. Jabouille   03/10/01 LHORELAX_ modifications
!!      Modifications  P. Jabouille   12/03/02 conditional writing of namelists
!!      Modifications  J.-P. Pinty    29/11/02 add C3R5, ICE2, ICE4, CELEC
!!      Modification   V. Masson      01/2004  removes surface (externalization)
!!      Modification   P. Tulet       01/2005  add dust, orilam
!!      Modification                  05/2006  Remove EPS and OWRIGET
!!      Modification    01/2016  (JP Pinty) Add LIMA
!!                   02/2018 Q.Libois ECRAD
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!      Modification   V. Vionnet     07/2017  add blowing snow variables
!!      Modification   F.Auguste      02/2021  add IBM
!!                     E.Jezequel     02/2021  add stations read from CSV file
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
USE MODD_CONF
USE MODD_DYN_n,   ONLY: LHORELAX_SVLIMA
USE MODD_IO,      ONLY: TFILEDATA
USE MODD_LUNIT_n, ONLY: TLUOUT
USE MODD_PARAMETERS
!
USE MODE_MSG
!
USE MODN_BACKUP
USE MODN_CONF
USE MODN_DYN
USE MODN_NESTING
USE MODN_OUTPUT
USE MODN_BUDGET
USE MODN_LES
USE MODN_DYN_n
USE MODN_ADV_n
USE MODN_PARAM_n
USE MODN_PARAM_RAD_n
USE MODN_PARAM_ECRAD_n
USE MODN_PARAM_KAFR_n
USE MODN_PARAM_MFSHALL_n
USE MODN_PARAM_ICE
USE MODN_CONF_n
USE MODN_LUNIT_n
USE MODN_LBC_n
USE MODN_NUDGING_n
USE MODN_TURB_n
USE MODN_BLANK_n
USE MODN_FRC
USE MODN_CH_MNHC_n
USE MODN_CH_SOLVER_n
USE MODN_PARAM_C2R2
USE MODN_PARAM_C1R3
USE MODN_PARAM_LIMA
USE MODN_ELEC
USE MODN_SERIES
USE MODN_SERIES_n 
USE MODN_TURB_CLOUD
USE MODN_TURB
USE MODN_CH_ORILAM
USE MODN_DUST
USE MODN_SALT
USE MODN_PASPOL
USE MODN_CONDSAMP
USE MODN_2D_FRC
USE MODN_LATZ_EDFLX
#ifdef MNH_FOREFIRE
USE MODN_FOREFIRE
USE MODD_FOREFIRE_n, ONLY : FFCOUPLING
#endif
USE MODN_BLOWSNOW_n
USE MODN_BLOWSNOW
USE MODN_IBM_PARAM_n
USE MODN_RECYCL_PARAM_n
USE MODD_IBM_LSF, ONLY: LIBM_LSF
USE MODN_STATION_n
!
IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
INTEGER,            INTENT(IN)  :: KMI     ! Model index
TYPE(TFILEDATA),    INTENT(IN)  :: TPDATAFILE ! Datafile
!
!*       0.2   declarations of local variables
!
INTEGER :: ILUSEG  ! logical unit number of EXSEG file
INTEGER :: ILUOUT  ! Logical unit number for output-listing TLUOUT file
!
LOGICAL                     ::  GHORELAX_UVWTH,                               &
                                GHORELAX_RV,  GHORELAX_RC, GHORELAX_RR,       &
                                GHORELAX_RI,  GHORELAX_RS, GHORELAX_RG,       &
                                GHORELAX_TKE, GHORELAX_SVC2R2, GHORELAX_SVPP, &
                                GHORELAX_SVCS, GHORELAX_SVCHIC,               &
#ifdef MNH_FOREFIRE
                                GHORELAX_SVFF,                                &
#endif
                                GHORELAX_SVCHEM, GHORELAX_SVC1R3,             &
                                GHORELAX_SVELEC, GHORELAX_SVLIMA,GHORELAX_SVSNW
LOGICAL                     ::  GHORELAX_SVDST, GHORELAX_SVSLT,  GHORELAX_SVAER
LOGICAL, DIMENSION(JPSVMAX) ::  GHORELAX_SV
!
!-------------------------------------------------------------------------------
!
!*       1.    UPDATE DESFM FILE
!              -----------------
!
CALL PRINT_MSG(NVERB_DEBUG,'IO','WRITE_DESFM_n','called for '//TRIM(TPDATAFILE%CNAME))
!
IF (.NOT.ASSOCIATED(TPDATAFILE%TDESFILE)) &
  CALL PRINT_MSG(NVERB_FATAL,'IO','WRITE_DESFM_n','TDESFILE not associated for '//TRIM(TPDATAFILE%CNAME))
!
ILUSEG = TPDATAFILE%TDESFILE%NLU
!
CALL INIT_NAM_LUNITn
WRITE(UNIT=ILUSEG,NML=NAM_LUNITn)
IF (CPROGRAM/='MESONH') THEN
  LUSECI=.FALSE.
  NSV_USER = 0      
ENDIF
CALL INIT_NAM_CONFn
WRITE(UNIT=ILUSEG,NML=NAM_CONFn)
!
!
CALL INIT_NAM_DYNn
IF (CPROGRAM/='MESONH') THEN   ! impose default value for next simulation
  GHORELAX_UVWTH = LHORELAX_UVWTH
  GHORELAX_RV    = LHORELAX_RV
  GHORELAX_RC    = LHORELAX_RC
  GHORELAX_RR    = LHORELAX_RR
  GHORELAX_RI    = LHORELAX_RI
  GHORELAX_RS    = LHORELAX_RS
  GHORELAX_RG    = LHORELAX_RG
  GHORELAX_TKE   = LHORELAX_TKE
  GHORELAX_SV(:) = LHORELAX_SV(:)
  GHORELAX_SVC2R2= LHORELAX_SVC2R2
  GHORELAX_SVC1R3= LHORELAX_SVC1R3
  GHORELAX_SVLIMA= LHORELAX_SVLIMA
  GHORELAX_SVELEC= LHORELAX_SVELEC
  GHORELAX_SVCHEM= LHORELAX_SVCHEM
  GHORELAX_SVCHIC= LHORELAX_SVCHIC
  GHORELAX_SVDST = LHORELAX_SVDST
  GHORELAX_SVSLT = LHORELAX_SVSLT
  GHORELAX_SVPP  = LHORELAX_SVPP 
#ifdef MNH_FOREFIRE
  GHORELAX_SVFF  = LHORELAX_SVFF
#endif
  GHORELAX_SVCS  = LHORELAX_SVCS 
  GHORELAX_SVAER = LHORELAX_SVAER
  GHORELAX_SVSNW = LHORELAX_SVSNW    
!
  LHORELAX_UVWTH = .FALSE.
  LHORELAX_RV    = .FALSE.
  LHORELAX_RC    = .FALSE.
  LHORELAX_RR    = .FALSE.
  LHORELAX_RI    = .FALSE.
  LHORELAX_RS    = .FALSE.
  LHORELAX_RG    = .FALSE.
  LHORELAX_TKE   = .FALSE.
  LHORELAX_SV(:) = .FALSE.
  LHORELAX_SVC2R2= .FALSE.
  LHORELAX_SVC1R3= .FALSE.
  LHORELAX_SVLIMA= .FALSE.
  LHORELAX_SVELEC= .FALSE.
  LHORELAX_SVCHEM= .FALSE.
  LHORELAX_SVCHIC= .FALSE.
  LHORELAX_SVLG  = .FALSE.
  LHORELAX_SVPP  = .FALSE.
#ifdef MNH_FOREFIRE
  LHORELAX_SVFF  = .FALSE.
#endif
  LHORELAX_SVCS  = .FALSE.
  LHORELAX_SVDST= .FALSE.
  LHORELAX_SVSLT= .FALSE.
  LHORELAX_SVAER= .FALSE.
  LHORELAX_SVSNW= .FALSE.    
ELSE  !return to namelist meaning of LHORELAX_SV
  GHORELAX_SV(:) = LHORELAX_SV(:)
  LHORELAX_SV(NSV_USER+1:)=.FALSE. 
END IF
WRITE(UNIT=ILUSEG,NML=NAM_DYNn)
!
IF (LIBM_LSF) THEN
  !
  CALL INIT_NAM_IBM_PARAMn
  !
  WRITE(UNIT=ILUSEG,NML=NAM_IBM_PARAMn)
  !
  IF (CPROGRAM/='MESONH') THEN
    LIBM         = .FALSE.
    LIBM_TROUBLE = .FALSE.  
    CIBM_ADV     = 'NOTHIN' 
  END IF
  !
END IF
!
CALL INIT_NAM_ADVn
WRITE(UNIT=ILUSEG,NML=NAM_ADVn)
IF (CPROGRAM/='MESONH') THEN
  CTURB   = 'NONE'
  CRAD    = 'NONE'
  CCLOUD  = 'NONE'
  CDCONV  = 'NONE'
  CSCONV  = 'NONE'
  CELEC   = 'NONE'
  CACTCCN = 'NONE'
END IF
CALL INIT_NAM_PARAMn
WRITE(UNIT=ILUSEG,NML=NAM_PARAMn)
!
CALL INIT_NAM_PARAM_RADn
IF(CRAD /= 'NONE') WRITE(UNIT=ILUSEG,NML=NAM_PARAM_RADn)
#ifdef MNH_ECRAD
CALL INIT_NAM_PARAM_ECRADn
IF(CRAD /= 'NONE') WRITE(UNIT=ILUSEG,NML=NAM_PARAM_ECRADn)
#endif
!
CALL INIT_NAM_PARAM_KAFRn
IF(CDCONV /= 'NONE' .OR. CSCONV == 'KAFR') &
       WRITE(UNIT=ILUSEG,NML=NAM_PARAM_KAFRn)
!
CALL INIT_NAM_PARAM_MFSHALLn
IF (CSCONV == 'EDKF' ) WRITE(UNIT=ILUSEG,NML=NAM_PARAM_MFSHALLn)
!
CALL INIT_NAM_LBCn
WRITE(UNIT=ILUSEG,NML=NAM_LBCn)
!
CALL INIT_NAM_NUDGINGn
WRITE(UNIT=ILUSEG,NML=NAM_NUDGINGn)
!
CALL INIT_NAM_TURBn
IF(CTURB /= 'NONE') WRITE(UNIT=ILUSEG,NML=NAM_TURBn)
!
CALL INIT_NAM_BLANKn
WRITE(UNIT=ILUSEG,NML=NAM_BLANKn)
!
!IF (CPROGRAM/='MESONH') THEN
!  LUSECHEM   = .FALSE.
!  LORILAM    = .FALSE.
!  LDEPOS_AER = .FALSE.
!  LDUST      = .FALSE.
!  LDEPOS_DST = .FALSE.
!  LSALT      = .FALSE.
!  LDEPOS_SLT = .FALSE.
!  LPASPOL    = .FALSE.
!  LCONDSAMP  = .FALSE.
!END IF
CALL INIT_NAM_CH_MNHCn
IF(LUSECHEM .OR. LCH_CONV_LINOX .OR. LCH_CONV_SCAV) &
 WRITE(UNIT=ILUSEG,NML=NAM_CH_MNHCn)
!
CALL INIT_NAM_CH_SOLVERn
IF(LUSECHEM) WRITE(UNIT=ILUSEG,NML=NAM_CH_SOLVERn)
!
CALL INIT_NAM_BLOWSNOWn
IF(LBLOWSNOW) WRITE(UNIT=ILUSEG,NML=NAM_BLOWSNOWn)
IF(LBLOWSNOW) WRITE(UNIT=ILUSEG,NML=NAM_BLOWSNOW)
!
CALL INIT_NAM_STATIONn
IF(LSTATION) WRITE(UNIT=ILUSEG,NML=NAM_STATIONn)
!
IF(LDUST) WRITE(UNIT=ILUSEG,NML=NAM_DUST)
IF(LSALT) WRITE(UNIT=ILUSEG,NML=NAM_SALT)
IF(LPASPOL) WRITE(UNIT=ILUSEG,NML=NAM_PASPOL)
#ifdef MNH_FOREFIRE
IF(FFCOUPLING) WRITE(UNIT=ILUSEG,NML=NAM_FOREFIRE)
#endif
IF(LCONDSAMP) WRITE(UNIT=ILUSEG,NML=NAM_CONDSAMP)
IF(LORILAM.AND.LUSECHEM) WRITE(UNIT=ILUSEG,NML=NAM_CH_ORILAM)
!
CALL INIT_NAM_SERIESn
IF(LSERIES) WRITE(UNIT=ILUSEG,NML=NAM_SERIESn)
IF(L2D_ADV_FRC .OR. L2D_REL_FRC) WRITE(UNIT=ILUSEG,NML=NAM_2D_FRC)
!
IF (LUV_FLX .OR. LTH_FLX) WRITE(UNIT=ILUSEG,NML=NAM_LATZ_EDFLX)
!
IF (CPROGRAM/='MESONH') THEN
  LLG  = .FALSE.
END IF
WRITE(UNIT=ILUSEG,NML=NAM_CONF)
WRITE(UNIT=ILUSEG,NML=NAM_DYN)
WRITE(UNIT=ILUSEG,NML=NAM_NESTING)
!WRITE(UNIT=ILUSEG,NML=NAM_BACKUP)
!WRITE(UNIT=ILUSEG,NML=NAM_OUTPUT)
IF(CBUTYPE /= 'NONE') THEN
  IF(CBUTYPE=='SKIP') CBUTYPE='CART'
  WRITE(UNIT=ILUSEG,NML=NAM_BUDGET)
END IF
IF(LBU_RU) WRITE(UNIT=ILUSEG,NML=NAM_BU_RU)
IF(LBU_RV) WRITE(UNIT=ILUSEG,NML=NAM_BU_RV)
IF(LBU_RW) WRITE(UNIT=ILUSEG,NML=NAM_BU_RW)
IF(LBU_RTH) WRITE(UNIT=ILUSEG,NML=NAM_BU_RTH)
IF(LBU_RTKE) WRITE(UNIT=ILUSEG,NML=NAM_BU_RTKE)
IF(LBU_RRV) WRITE(UNIT=ILUSEG,NML=NAM_BU_RRV)
IF(LBU_RRC) WRITE(UNIT=ILUSEG,NML=NAM_BU_RRC)
IF(LBU_RRR) WRITE(UNIT=ILUSEG,NML=NAM_BU_RRR)
IF(LBU_RRI) WRITE(UNIT=ILUSEG,NML=NAM_BU_RRI)
IF(LBU_RRS) WRITE(UNIT=ILUSEG,NML=NAM_BU_RRS)
IF(LBU_RRG) WRITE(UNIT=ILUSEG,NML=NAM_BU_RRG)
IF(LBU_RRH) WRITE(UNIT=ILUSEG,NML=NAM_BU_RRH)
IF(LBU_RSV) WRITE(UNIT=ILUSEG,NML=NAM_BU_RSV)
IF(LLES_MEAN .OR. LLES_RESOLVED .OR. LLES_SUBGRID .OR. LLES_UPDRAFT  &
.OR. LLES_DOWNDRAFT .OR. LLES_SPECTRA) WRITE(UNIT=ILUSEG,NML=NAM_LES)
WRITE(UNIT=ILUSEG,NML=NAM_BLANKn)
IF(LFORCING .OR. LTRANS) WRITE(UNIT=ILUSEG,NML=NAM_FRC)
IF(CCLOUD(1:3) == 'ICE')  WRITE(UNIT=ILUSEG,NML=NAM_PARAM_ICE)
IF(CCLOUD == 'C2R2' .OR. CCLOUD == 'C3R5' .OR. CCLOUD == 'KHKO') &
                     WRITE(UNIT=ILUSEG,NML=NAM_PARAM_C2R2)
IF(CCLOUD == 'C3R5' ) WRITE(UNIT=ILUSEG,NML=NAM_PARAM_C1R3) 
IF(CCLOUD == 'LIMA' ) WRITE(UNIT=ILUSEG,NML=NAM_PARAM_LIMA) 
IF(CELEC /= 'NONE') WRITE(UNIT=ILUSEG,NML=NAM_ELEC) 
IF(LSERIES) WRITE(UNIT=ILUSEG,NML=NAM_SERIES)
IF(NMODEL_CLOUD/=NUNDEF) WRITE(UNIT=ILUSEG,NML=NAM_TURB_CLOUD)
IF(CTURB /= 'NONE') WRITE(UNIT=ILUSEG,NML=NAM_TURB)
!
!
!
!-------------------------------------------------------------------------------
!
!*       2.    WRITE UPDATED DESFM ON OUTPUT LISTING
!              -------------------------------------
!
IF (NVERB >= 5) THEN
!
  ILUOUT = TLUOUT%NLU
!
  WRITE(UNIT=ILUOUT,FMT="(/,'DESCRIPTOR OF SEGMENT FOR MODEL ',I2)") KMI
  WRITE(UNIT=ILUOUT,FMT="(  '------------------------------- '   )")
!
  WRITE(UNIT=ILUOUT,FMT="('********** LOGICAL UNITSn **********')")
  WRITE(UNIT=ILUOUT,NML=NAM_LUNITn)
!  
  WRITE(UNIT=ILUOUT,FMT="('********** CONFIGURATIONn **********')")
  WRITE(UNIT=ILUOUT,NML=NAM_CONFn)
!  
!  
  WRITE(UNIT=ILUOUT,FMT="('********** DYNAMICn ****************')")
  WRITE(UNIT=ILUOUT,NML=NAM_DYNn)
!  
  WRITE(UNIT=ILUOUT,FMT="('********** ADVECTIONn **************')")
  WRITE(UNIT=ILUOUT,NML=NAM_ADVn)
  !  
  IF (LIBM_LSF) THEN
    WRITE(UNIT=ILUOUT,FMT="('********** IBM_PARAMn **************')")                           
    WRITE(UNIT=ILUOUT,NML=NAM_IBM_PARAMn)
  ENDIF
  !
  IF (LRECYCL) THEN
    WRITE(UNIT=ILUOUT,FMT="('********** RECYCL_PARAMn **************')")
    WRITE(UNIT=ILUOUT,NML=NAM_RECYCL_PARAMn)
  ENDIF
  !  
  WRITE(UNIT=ILUOUT,FMT="('********** PARAMETERIZATIONSn ******')")
  WRITE(UNIT=ILUOUT,NML=NAM_PARAMn)
!  
  WRITE(UNIT=ILUOUT,FMT="('********** RADIATIONn **************')")
  WRITE(UNIT=ILUOUT,NML=NAM_PARAM_RADn)
#ifdef MNH_ECRAD
  WRITE(UNIT=ILUOUT,FMT="('********** ECRADn **************')")
  WRITE(UNIT=ILUOUT,NML=NAM_PARAM_ECRADn)
#endif
!  
  WRITE(UNIT=ILUOUT,FMT="('********** CONVECTIONn *************')")
  WRITE(UNIT=ILUOUT,NML=NAM_PARAM_KAFRn)
!  
  WRITE(UNIT=ILUOUT,FMT="('************ PARAM_MFSHALLn  *******')")
  WRITE(UNIT=ILUOUT,NML=NAM_PARAM_MFSHALLn)
!
  WRITE(UNIT=ILUOUT,FMT="('********** LBCn ********************')")
  WRITE(UNIT=ILUOUT,NML=NAM_LBCn)
!
  WRITE(UNIT=ILUOUT,FMT="('********** NUDGINGn*****************')")  
  WRITE(UNIT=ILUOUT,NML=NAM_NUDGINGn)
!  
  WRITE(UNIT=ILUOUT,FMT="('********** TURBn *******************')")  
  WRITE(UNIT=ILUOUT,NML=NAM_TURBn)
!  
  WRITE(UNIT=ILUOUT,FMT="('********** CHEMICAL MONITORn *******')")  
  WRITE(UNIT=ILUOUT,NML=NAM_CH_MNHCn)
!
  WRITE(UNIT=ILUOUT,FMT="('************ CHEMICAL SOLVERn ******************')")
  WRITE(UNIT=ILUOUT,NML=NAM_CH_SOLVERn)
!
  WRITE(UNIT=ILUOUT,FMT="('************ TEMPORAL SERIESn ******************')")
  WRITE(UNIT=ILUOUT,NML=NAM_SERIESn)
!  
  WRITE(UNIT=ILUOUT,FMT="('********** BLOWING SNOW SCHEME ****************')")
  WRITE(UNIT=ILUOUT,NML=NAM_BLOWSNOWn)
!
  WRITE(UNIT=ILUOUT,FMT="('********** BLANKn *****************************')")
  WRITE(UNIT=ILUOUT,NML=NAM_BLANKn)
!
  IF (KMI==1) THEN
    WRITE(UNIT=ILUOUT,FMT="(/,'PART OF SEGMENT FILE COMMON TO ALL THE MODELS')")
    WRITE(UNIT=ILUOUT,FMT="(  '---------------------------------------------')")
!
    WRITE(UNIT=ILUOUT,FMT="('************ CONFIGURATION ********************')")
    WRITE(UNIT=ILUOUT,NML=NAM_CONF)
!    
    WRITE(UNIT=ILUOUT,FMT="('************ DYNAMIC **************************')")
    WRITE(UNIT=ILUOUT,NML=NAM_DYN)
!
    WRITE(UNIT=ILUOUT,FMT="(/,'********** NESTING **************************')")
    WRITE(UNIT=ILUOUT,NML=NAM_NESTING)
!
!    WRITE(UNIT=ILUOUT,FMT="(/,'********** BACKUP ***************************')")
!    WRITE(UNIT=ILUOUT,NML=NAM_BACKUP)
!
!    WRITE(UNIT=ILUOUT,FMT="(/,'********** OUTPUT ***************************')")
!    WRITE(UNIT=ILUOUT,NML=NAM_OUTPUT)
!
    WRITE(UNIT=ILUOUT,FMT="('************ BUDGET ***************************')")
    WRITE(UNIT=ILUOUT,NML=NAM_BUDGET)
!    
    IF ( .NOT. ALLOCATED( CBULIST_RU ) ) ALLOCATE( CHARACTER(LEN=NBULISTMAXLEN) :: CBULIST_RU(0) )
    WRITE(UNIT=ILUOUT,FMT="('************ U BUDGET *************************')")
    WRITE(UNIT=ILUOUT,NML=NAM_BU_RU)
!
    IF ( .NOT. ALLOCATED( CBULIST_RV ) ) ALLOCATE( CHARACTER(LEN=NBULISTMAXLEN) :: CBULIST_RV(0) )
    WRITE(UNIT=ILUOUT,FMT="('************ V BUDGET *************************')")
    WRITE(UNIT=ILUOUT,NML=NAM_BU_RV)
!    
    IF ( .NOT. ALLOCATED( CBULIST_RW ) ) ALLOCATE( CHARACTER(LEN=NBULISTMAXLEN) :: CBULIST_RW(0) )
    WRITE(UNIT=ILUOUT,FMT="('************ W BUDGET *************************')")
    WRITE(UNIT=ILUOUT,NML=NAM_BU_RW)
!    
    IF ( .NOT. ALLOCATED( CBULIST_RTH ) ) ALLOCATE( CHARACTER(LEN=NBULISTMAXLEN) :: CBULIST_RTH(0) )
    WRITE(UNIT=ILUOUT,FMT="('************ TH BUDGET ************************')")
    WRITE(UNIT=ILUOUT,NML=NAM_BU_RTH)
!    
    IF ( .NOT. ALLOCATED( CBULIST_RTKE ) ) ALLOCATE( CHARACTER(LEN=NBULISTMAXLEN) :: CBULIST_RTKE(0) )
    WRITE(UNIT=ILUOUT,FMT="('************ TKE BUDGET ***********************')")
    WRITE(UNIT=ILUOUT,NML=NAM_BU_RTKE)
!    
    IF ( .NOT. ALLOCATED( CBULIST_RRV ) ) ALLOCATE( CHARACTER(LEN=NBULISTMAXLEN) :: CBULIST_RRV(0) )
    WRITE(UNIT=ILUOUT,FMT="('************ RV BUDGET ************************')")
    WRITE(UNIT=ILUOUT,NML=NAM_BU_RRV)
!    
    IF ( .NOT. ALLOCATED( CBULIST_RRC ) ) ALLOCATE( CHARACTER(LEN=NBULISTMAXLEN) :: CBULIST_RRC(0) )
    WRITE(UNIT=ILUOUT,FMT="('************ RC BUDGET ************************')")
    WRITE(UNIT=ILUOUT,NML=NAM_BU_RRC)
!    
    IF ( .NOT. ALLOCATED( CBULIST_RRR ) ) ALLOCATE( CHARACTER(LEN=NBULISTMAXLEN) :: CBULIST_RRR(0) )
    WRITE(UNIT=ILUOUT,FMT="('************ RR BUDGET ************************')")
    WRITE(UNIT=ILUOUT,NML=NAM_BU_RRR)
!    
    IF ( .NOT. ALLOCATED( CBULIST_RRI ) ) ALLOCATE( CHARACTER(LEN=NBULISTMAXLEN) :: CBULIST_RRI(0) )
    WRITE(UNIT=ILUOUT,FMT="('************ RI BUDGET ************************')")
    WRITE(UNIT=ILUOUT,NML=NAM_BU_RRI)
!    
    IF ( .NOT. ALLOCATED( CBULIST_RRS ) ) ALLOCATE( CHARACTER(LEN=NBULISTMAXLEN) :: CBULIST_RRS(0) )
    WRITE(UNIT=ILUOUT,FMT="('************ RS BUDGET ************************')")
    WRITE(UNIT=ILUOUT,NML=NAM_BU_RRS)
!    
    IF ( .NOT. ALLOCATED( CBULIST_RRG ) ) ALLOCATE( CHARACTER(LEN=NBULISTMAXLEN) :: CBULIST_RRG(0) )
    WRITE(UNIT=ILUOUT,FMT="('************ RG BUDGET ************************')")
    WRITE(UNIT=ILUOUT,NML=NAM_BU_RRG)
!    
    IF ( .NOT. ALLOCATED( CBULIST_RRH ) ) ALLOCATE( CHARACTER(LEN=NBULISTMAXLEN) :: CBULIST_RRH(0) )
    WRITE(UNIT=ILUOUT,FMT="('************ RH BUDGET ************************')")
    WRITE(UNIT=ILUOUT,NML=NAM_BU_RRH)
!
    IF ( .NOT. ALLOCATED( CBULIST_RSV ) ) ALLOCATE( CHARACTER(LEN=NBULISTMAXLEN) :: CBULIST_RSV(0) )
    WRITE(UNIT=ILUOUT,FMT="('************ SVx BUDGET ***********************')")
    WRITE(UNIT=ILUOUT,NML=NAM_BU_RSV)
!    
    WRITE(UNIT=ILUOUT,FMT="('************ LES ******************************')")
    WRITE(UNIT=ILUOUT,NML=NAM_LES)
!    
    WRITE(UNIT=ILUOUT,FMT="('************ FORCING **************************')")
    WRITE(UNIT=ILUOUT,NML=NAM_FRC)
!    
    WRITE(UNIT=ILUOUT,FMT="('************ ICE SCHEME ***********************')")
    WRITE(UNIT=ILUOUT,NML=NAM_PARAM_ICE)
!    
    WRITE(UNIT=ILUOUT,FMT="('********** DUST SCHEME ************************')")
    WRITE(UNIT=ILUOUT,NML=NAM_DUST)
!
    WRITE(UNIT=ILUOUT,FMT="('********** PASPOL *****************************')")
    WRITE(UNIT=ILUOUT,NML=NAM_PASPOL)
!
#ifdef MNH_FOREFIRE
    WRITE(UNIT=ILUOUT,FMT="('********** FOREFIRE *****************************')")
    WRITE(UNIT=ILUOUT,NML=NAM_FOREFIRE)
!
#endif
    WRITE(UNIT=ILUOUT,FMT="('********** CONDSAMP****************************')")
    WRITE(UNIT=ILUOUT,NML=NAM_CONDSAMP)
!
    WRITE(UNIT=ILUOUT,FMT="('********** SALT SCHEME ************************')")
    WRITE(UNIT=ILUOUT,NML=NAM_SALT)
!
    WRITE(UNIT=ILUOUT,FMT="('********** BLOWING SNOW SCHEME ****************')")
    WRITE(UNIT=ILUOUT,NML=NAM_BLOWSNOW)    
!
    WRITE(UNIT=ILUOUT,FMT="('************ ORILAM SCHEME ********************')")
    WRITE(UNIT=ILUOUT,NML=NAM_CH_ORILAM)
!    
    IF( CCLOUD == 'C2R2' .OR. CCLOUD == 'C3R5') THEN
      WRITE(UNIT=ILUOUT,FMT="('*********** C2R2 SCHEME *********************')")
      WRITE(UNIT=ILUOUT,NML=NAM_PARAM_C2R2)
      IF( CCLOUD == 'C3R5' ) THEN
        WRITE(UNIT=ILUOUT,FMT="('*********** C1R3 SCHEME *********************')")
        WRITE(UNIT=ILUOUT,NML=NAM_PARAM_C1R3)
      END IF
    END IF
!
    IF( CCLOUD == 'LIMA' ) THEN
      WRITE(UNIT=ILUOUT,FMT="('*********** LIMA SCHEME *********************')")
      WRITE(UNIT=ILUOUT,NML=NAM_PARAM_LIMA)
    END IF
!
    IF( CCLOUD == 'KHKO' ) THEN
      WRITE(UNIT=ILUOUT,FMT="('*********** KHKO SCHEME *********************')")
      WRITE(UNIT=ILUOUT,NML=NAM_PARAM_C2R2)
    END IF
!
    IF( CELEC /= 'NONE' ) THEN
      WRITE(UNIT=ILUOUT,FMT="('*********** ELEC SCHEME *********************')")
      WRITE(UNIT=ILUOUT,NML=NAM_ELEC)
    END IF
!
    WRITE(UNIT=ILUOUT,FMT="('************ TEMPORAL SERIES ****************')")
    WRITE(UNIT=ILUOUT,NML=NAM_SERIES)
!
    WRITE(UNIT=ILUOUT,FMT="('************ MIXING LENGTH FOR CLOUD ***********')")
    WRITE(UNIT=ILUOUT,NML=NAM_TURB_CLOUD)
!
  END IF
!
END IF  
!
IF (CPROGRAM /='MESONH') THEN !return to previous LHORELAX_
  LHORELAX_UVWTH = GHORELAX_UVWTH
  LHORELAX_RV    = GHORELAX_RV
  LHORELAX_RC    = GHORELAX_RC
  LHORELAX_RR    = GHORELAX_RR
  LHORELAX_RI    = GHORELAX_RI
  LHORELAX_RS    = GHORELAX_RS
  LHORELAX_RG    = GHORELAX_RG
  LHORELAX_TKE   = GHORELAX_TKE
  LHORELAX_SV(:) = GHORELAX_SV(:)
  LHORELAX_SVC2R2= GHORELAX_SVC2R2
  LHORELAX_SVC1R3= GHORELAX_SVC1R3
  LHORELAX_SVLIMA= GHORELAX_SVLIMA
  LHORELAX_SVELEC= GHORELAX_SVELEC
  LHORELAX_SVCHEM= GHORELAX_SVCHEM
  LHORELAX_SVCHIC= GHORELAX_SVCHIC
  LHORELAX_SVLG  = .FALSE.
  LHORELAX_SVDST = GHORELAX_SVDST
  LHORELAX_SVSLT = GHORELAX_SVSLT
  LHORELAX_SVPP  = GHORELAX_SVPP 
#ifdef MNH_FOREFIRE
  LHORELAX_SVFF  = GHORELAX_SVFF
#endif
  LHORELAX_SVCS  = GHORELAX_SVCS 
  LHORELAX_SVAER = GHORELAX_SVAER
  LHORELAX_SVSNW = GHORELAX_SVSNW  
ELSE
  LHORELAX_SV(:) = GHORELAX_SV(:)
ENDIF
CALL UPDATE_NAM_DYNn
!------------------------------------------------------------------------------
!
END SUBROUTINE WRITE_DESFM_n
