!MNH_LIC Copyright 1994-2023 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ###########################
      MODULE MODI_ADVECTION_METSV
!     ###########################
!
INTERFACE
      SUBROUTINE ADVECTION_METSV (TPFILE, HUVW_ADV_SCHEME,                     &
                            HMET_ADV_SCHEME,HSV_ADV_SCHEME, HCLOUD, KSPLIT,    &
                            OSPLIT_CFL, PSPLIT_CFL, OCFL_WRIT,                 &
                            HLBCX, HLBCY, KRR, KSV, TPDTCUR, PTSTEP,           &
                            PUT, PVT, PWT, PTHT, PRT, PTKET, PSVT, PPABST,     &
                            PTHVREF, PRHODJ, PDXX, PDYY, PDZZ, PDZX, PDZY,     &
                            PRTHS, PRRS, PRTKES, PRSVS,                        &
                            PRTHS_CLD, PRRS_CLD, PRSVS_CLD, PRTKES_ADV         )
!
USE MODD_IO,        ONLY: TFILEDATA
USE MODD_TYPE_DATE, ONLY: DATE_TIME
!
TYPE(TFILEDATA),        INTENT(IN)   ::  TPFILE       ! Output file
CHARACTER(LEN=6),       INTENT(IN)   :: HMET_ADV_SCHEME, & ! Control of the 
                                        HSV_ADV_SCHEME, &  ! scheme applied 
                                        HUVW_ADV_SCHEME
CHARACTER (LEN=4),      INTENT(IN)   :: HCLOUD      ! Kind of cloud parameterization                                
!
INTEGER,                INTENT(INOUT):: KSPLIT       ! Number of time splitting
                                                     ! for PPM advection
LOGICAL,                INTENT(IN)   :: OSPLIT_CFL   ! flag to automatically chose number of iterations
REAL,                   INTENT(IN)   :: PSPLIT_CFL   ! maximum CFL to automatically chose number of iterations
LOGICAL,                INTENT(IN)   :: OCFL_WRIT    ! flag to write CFL fields in output files            
!
CHARACTER(LEN=4),DIMENSION(2),INTENT(IN):: HLBCX, HLBCY  ! X- and Y-direc LBC
!
INTEGER,                  INTENT(IN)    :: KRR  ! Number of moist variables
INTEGER,                  INTENT(IN)    :: KSV  ! Number of Scalar Variables
!
TYPE (DATE_TIME),         INTENT(IN)    :: TPDTCUR ! current date and time
REAL,                     INTENT(IN)    :: PTSTEP
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PUT , PVT  , PWT
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT, PTKET, PRHODJ
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABST                 
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRT , PSVT
                                                  ! Variables at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHVREF   ! Virtual Temperature
                                          ! of the reference state
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PDXX,PDYY,PDZZ
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PDZX,PDZY
                                                  !  metric coefficients
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRTHS, PRTKES
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRRS , PRSVS
                                                  ! Sources terms 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRTHS_CLD
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRRS_CLD,PRSVS_CLD
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PRTKES_ADV  ! Advection TKE source term 
!
END SUBROUTINE ADVECTION_METSV
!
END INTERFACE
!
END MODULE MODI_ADVECTION_METSV
!     ##########################################################################
      SUBROUTINE ADVECTION_METSV (TPFILE, HUVW_ADV_SCHEME,                     &
                            HMET_ADV_SCHEME,HSV_ADV_SCHEME, HCLOUD, KSPLIT,    &
                            OSPLIT_CFL, PSPLIT_CFL, OCFL_WRIT,                 &
                            HLBCX, HLBCY, KRR, KSV, TPDTCUR, PTSTEP,           &
                            PUT, PVT, PWT, PTHT, PRT, PTKET, PSVT, PPABST,     &
                            PTHVREF, PRHODJ, PDXX, PDYY, PDZZ, PDZX, PDZY,     &
                            PRTHS, PRRS, PRTKES, PRSVS,                        &
                            PRTHS_CLD, PRRS_CLD, PRSVS_CLD, PRTKES_ADV         )
!     ##########################################################################
!
!!****  *ADVECTION_METSV * - routine to call the specialized advection routines
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to control the advection routines.
!!    For that, it is first necessary to compute the metric coefficients
!!    and the contravariant components of the momentum.
!!
!!**  METHOD
!!    ------
!!      Once the scheme is selected, it is applied to the following group of
!!    variables: METeorologicals (temperature, water substances, TKE,
!!    dissipation TKE) and Scalar Variables. It is possible to select different
!!    advection schemes for each group of variables.
!!
!!    EXTERNAL
!!    --------
!!      CONTRAV              : computes the contravariant components.
!!      ADVECUVW             : computes the advection terms for momentum.
!!      ADVECSCALAR          : computes the advection terms for scalar fields.
!!      ADD3DFIELD_ll        : add a field to 3D-list
!!      ADVEC_4TH_ORDER      : 4th order advection scheme
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      NONE
!!
!!    REFERENCE
!!    ---------
!!      Book1 and book2 ( routine ADVECTION )
!!
!!    AUTHOR
!!    ------
!!	J.-P. Pinty      * Laboratoire d'Aerologie*
!!	J.-P. Lafore     * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    06/07/94 
!!                  01/04/95 (Ph. Hereil J. Nicolau) add the model number
!!                  23/10/95 (J. Vila and JP Lafore) advection schemes scalar
!!                  16/01/97 (JP Pinty)              change presentation 
!!                  30/04/98 (J. Stein P Jabouille)  extrapolation for the cyclic
!!                                                   case and parallelisation
!!                  24/06/99 (P Jabouille)           case of NHALO>1
!!                  25/10/05 (JP Pinty)              4th order scheme
!!                  24/04/06 (C.Lac)                 Split scalar and passive
!!                                                   tracer routines
!!                  08/06    (T.Maric)               PPM scheme
!!                  04/2011  (V.Masson & C. Lac)     splits the routine and add time splitting
!!                  04/2014  (C.Lac)                 adaptation of time
!!                                                   splitting for L1D and L2D
!!                  09/2014  (G.Delautier)              close OUTPUT_LISTING before STOP
!!                  04/2015  (J.Escobar) remove/commente some NHALO=1 test
!!                  J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!!                  J.Escobar : 01/10/2015 : add computation of CFL for L1D case
!!                  04/2016  (C.Lac)       : correction of negativity for KHKO
!!                  10/2016  (C.Lac) Correction on the flag for Strang splitting
!!                                  to insure reproducibility between START and RESTA
!  V. Vionnet     07/2017: add advection of 2D variables at the surface for the blowing snow scheme
!  P. Wautelet 05/2016-04/2018: new data structures and calls for I/O
!  P. Wautelet    02/2020: use the new data structures and subroutines for budgets
!  B. Vie         03/2020: LIMA negativity checks after turbulence, advection and microphysics budgets
!  P. Wautelet 20/05/2019: add name argument to ADDnFIELD_ll + new ADD4DFIELD_ll subroutine
!  P. Wautelet 11/06/2020: bugfix: correct PRSVS array indices
!  P. Wautelet + Benoît Vié 06/2020: improve removal of negative scalar variables + adapt the corresponding budgets
!  P. Wautelet 30/06/2020: move removal of negative scalar variables to Sources_neg_correct
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
use modd_budget,     only: lbudget_th, lbudget_tke, lbudget_rv, lbudget_rc,                          &
                           lbudget_rr, lbudget_ri,  lbudget_rs, lbudget_rg, lbudget_rh, lbudget_sv,  &
                           NBUDGET_TH, NBUDGET_TKE, NBUDGET_RV, NBUDGET_RC,                          &
                           NBUDGET_RR, NBUDGET_RI,  NBUDGET_RS, NBUDGET_RG, NBUDGET_RH, NBUDGET_SV1, &
                           tbudgets
USE MODD_CST
USE MODD_TURB_n,         ONLY: XTKEMIN
USE MODD_CONF,           ONLY: LNEUTRAL,NHALO,L1D, L2D
use modd_field,          only: tfieldmetadata, TYPEREAL
USE MODD_IBM_PARAM_n,    ONLY: LIBM,XIBM_LS,XIBM_EPSI
USE MODD_IO,             ONLY: TFILEDATA
USE MODD_LUNIT_n,        ONLY: TLUOUT
USE MODD_NSV
USE MODD_PARAM_LIMA
USE MODD_PARAM_n
USE MODD_TYPE_DATE,      ONLY: DATE_TIME
USE MODD_BLOWSNOW
USE MODD_BLOWSNOW_n
USE MODD_PARAMETERS
USE MODD_REF_n,          ONLY: XRHODJ,XRHODREF
!
use mode_budget,         only: Budget_store_init, Budget_store_end
USE MODE_IO_FIELD_WRITE, only: IO_Field_write
USE MODE_ll
USE MODE_MSG
use mode_sources_neg_correct, only: Sources_neg_correct
!
USE MODI_ADV_BOUNDARIES
USE MODI_CONTRAV
USE MODI_GET_HALO
USE MODI_PPM_RHODJ
USE MODI_PPM_MET
USE MODI_PPM_SCALAR
!
!
!-------------------------------------------------------------------------------
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(TFILEDATA),        INTENT(IN)   ::  TPFILE       ! Output file
CHARACTER(LEN=6),       INTENT(IN)   :: HMET_ADV_SCHEME, & ! Control of the 
                                        HSV_ADV_SCHEME, &  ! scheme applied 
                                        HUVW_ADV_SCHEME
CHARACTER (LEN=4),      INTENT(IN)   :: HCLOUD      ! Kind of cloud parameterization                                
!
INTEGER,                INTENT(INOUT):: KSPLIT       ! Number of time splitting
                                                     ! for PPM advection
LOGICAL,                INTENT(IN)   :: OSPLIT_CFL   ! flag to automatically chose number of iterations
REAL,                   INTENT(IN)   :: PSPLIT_CFL   ! maximum CFL to automatically chose number of iterations
LOGICAL,                INTENT(IN)   :: OCFL_WRIT    ! flag to write CFL fields in output files            
!
CHARACTER(LEN=4),DIMENSION(2),INTENT(IN):: HLBCX, HLBCY  ! X- and Y-direc LBC
!
INTEGER,                  INTENT(IN)    :: KRR  ! Number of moist variables
INTEGER,                  INTENT(IN)    :: KSV  ! Number of Scalar Variables
!
TYPE (DATE_TIME),         INTENT(IN)    :: TPDTCUR ! current date and time
REAL,                     INTENT(IN)    :: PTSTEP
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PUT , PVT  , PWT
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT, PTKET, PRHODJ
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABST                 
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRT , PSVT
                                                  ! Variables at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHVREF   ! Virtual Temperature
                                          ! of the reference state
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PDXX,PDYY,PDZZ
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PDZX,PDZY
                                                  !  metric coefficients
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRTHS, PRTKES
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRRS , PRSVS
                                                  ! Sources terms 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRTHS_CLD
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRRS_CLD, PRSVS_CLD
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PRTKES_ADV  ! Advection TKE source term 
!
!
!*       0.2   declarations of local variables
!
!
REAL, DIMENSION(SIZE(PUT,1),SIZE(PUT,2),SIZE(PUT,3)) :: ZRUCPPM
REAL, DIMENSION(SIZE(PUT,1),SIZE(PUT,2),SIZE(PUT,3)) :: ZRVCPPM
REAL, DIMENSION(SIZE(PUT,1),SIZE(PUT,2),SIZE(PUT,3)) :: ZRWCPPM
                                                  ! contravariant
                                                  ! components
                                                  ! of momentum
REAL, DIMENSION(SIZE(PUT,1),SIZE(PUT,2),SIZE(PUT,3)) :: ZCFLU
REAL, DIMENSION(SIZE(PUT,1),SIZE(PUT,2),SIZE(PUT,3)) :: ZCFLV
REAL, DIMENSION(SIZE(PUT,1),SIZE(PUT,2),SIZE(PUT,3)) :: ZCFLW
!                                                 ! CFL numbers on each direction
REAL, DIMENSION(SIZE(PUT,1),SIZE(PUT,2),SIZE(PUT,3)) :: ZCFL
!                                                 ! CFL number
!
REAL :: ZCFLU_MAX, ZCFLV_MAX, ZCFLW_MAX, ZCFL_MAX ! maximum CFL numbers
!
REAL, DIMENSION(SIZE(PTHT,1), SIZE(PTHT,2), SIZE(PTHT,3) ) :: ZTH
REAL, DIMENSION(SIZE(PTKET,1),SIZE(PTKET,2),SIZE(PTKET,3)) :: ZTKE
REAL, DIMENSION(SIZE(PTHT,1), SIZE(PTHT,2), SIZE(PTHT,3) ) :: ZRTHS_OTHER
REAL, DIMENSION(SIZE(PTKET,1),SIZE(PTKET,2),SIZE(PTKET,3)) :: ZRTKES_OTHER
REAL, DIMENSION(SIZE(PTHT,1), SIZE(PTHT,2), SIZE(PTHT,3) ) :: ZRTHS_PPM
REAL, DIMENSION(SIZE(PTKET,1),SIZE(PTKET,2),SIZE(PTKET,3)) :: ZRTKES_PPM
REAL, DIMENSION(SIZE(PRT,1), SIZE(PRT,2), SIZE(PRT,3), SIZE(PRT,4) ) :: ZR
REAL, DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3),SIZE(PSVT,4)) :: ZSV
REAL, DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3), NBLOWSNOW_2D) :: ZSNWC
REAL, DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3), NBLOWSNOW_2D) :: ZSNWC_INIT
REAL, DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3), NBLOWSNOW_2D) :: ZRSNWCS
! Guess at the sub time step
REAL, DIMENSION(SIZE(PRT,1), SIZE(PRT,2), SIZE(PRT,3), SIZE(PRT,4) ) :: ZRRS_OTHER
REAL, DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3),SIZE(PSVT,4)) :: ZRSVS_OTHER
REAL, DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3),NBLOWSNOW_2D) ::  ZRSNWCS_OTHER
! Tendencies since the beginning of the time step
REAL, DIMENSION(SIZE(PRT,1), SIZE(PRT,2), SIZE(PRT,3), SIZE(PRT,4) ) :: ZRRS_PPM
REAL, DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3),SIZE(PSVT,4)) :: ZRSVS_PPM
REAL, DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3),NBLOWSNOW_2D) :: ZRSNWCS_PPM
! Guess at the end of the sub time step
REAL, DIMENSION(SIZE(PUT,1),SIZE(PUT,2),SIZE(PUT,3)) :: ZRHOX1,ZRHOX2
REAL, DIMENSION(SIZE(PUT,1),SIZE(PUT,2),SIZE(PUT,3)) :: ZRHOY1,ZRHOY2
REAL, DIMENSION(SIZE(PUT,1),SIZE(PUT,2),SIZE(PUT,3)) :: ZRHOZ1,ZRHOZ2
REAL, DIMENSION(SIZE(PUT,1),SIZE(PUT,2),SIZE(PUT,3)):: ZT,ZEXN,ZLV,ZLS,ZCPH,ZCOR
! Temporary advected rhodj for PPM routines
!
INTEGER :: JS,JR,JSV,JSPL, JI, JJ  ! Loop index
REAL    :: ZTSTEP_PPM ! Sub Time step 
LOGICAL :: GTKE
!
INTEGER                     :: IINFO_ll    ! return code of parallel routine
TYPE(LIST_ll), POINTER      :: TZFIELDS0_ll ! list of fields to exchange
TYPE(LIST_ll), POINTER      :: TZFIELDS1_ll ! list of fields to exchange
!
!
INTEGER              :: IRESP        ! Return code of FM routines
INTEGER              :: ILUOUT       ! logical unit
INTEGER              :: ISPLIT_PPM   ! temporal time splitting
INTEGER              :: IIB, IIE, IJB, IJE,IKB,IKE
TYPE(TFIELDMETADATA) :: TZFIELD
!-------------------------------------------------------------------------------
!
!*       0.     INITIALIZATION                        
!	        --------------

GTKE=(SIZE(PTKET)/=0)

if ( lbudget_th  ) call Budget_store_init( tbudgets(NBUDGET_TH ), 'ADV', prths (:, :, :)    )
if ( lbudget_tke ) call Budget_store_init( tbudgets(NBUDGET_TKE), 'ADV', prtkes(:, :, :)    )
if ( lbudget_rv  ) call Budget_store_init( tbudgets(NBUDGET_RV ), 'ADV', prrs  (:, :, :, 1) )
if ( lbudget_rc  ) call Budget_store_init( tbudgets(NBUDGET_RC ), 'ADV', prrs  (:, :, :, 2) )
if ( lbudget_rr  ) call Budget_store_init( tbudgets(NBUDGET_RR ), 'ADV', prrs  (:, :, :, 3) )
if ( lbudget_ri  ) call Budget_store_init( tbudgets(NBUDGET_RI ), 'ADV', prrs  (:, :, :, 4) )
if ( lbudget_rs  ) call Budget_store_init( tbudgets(NBUDGET_RS ), 'ADV', prrs  (:, :, :, 5) )
if ( lbudget_rg  ) call Budget_store_init( tbudgets(NBUDGET_RG ), 'ADV', prrs  (:, :, :, 6) )
if ( lbudget_rh  ) call Budget_store_init( tbudgets(NBUDGET_RH ), 'ADV', prrs  (:, :, :, 7) )
if ( lbudget_sv) then
  do jsv = 1, ksv
    call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + jsv ), 'ADV', prsvs(:, :, :, jsv) )
  end do
end if

ILUOUT = TLUOUT%NLU
!
CALL GET_INDICE_ll(IIB,IJB,IIE,IJE)
IKB=1+JPVEXT
IKE=SIZE(PSVT,3) - JPVEXT
!
IF(LBLOWSNOW) THEN    ! Put 2D Canopy blowing snow variables into a 3D array for advection
  ZSNWC_INIT = 0.
  ZRSNWCS = 0.

  DO JSV=1,(NBLOWSNOW_2D)
     ZSNWC_INIT(:,:,IKB,JSV) = XSNWCANO(:,:,JSV)
     ZRSNWCS(:,:,IKB,JSV)    = XRSNWCANOS(:,:,JSV)
  END DO
ENDIF
!
!
!-------------------------------------------------------------------------------
!
!*       2.     COMPUTES THE CONTRAVARIANT COMPONENTS (FOR PPM ONLY)
!	        --------------------------------------
!
!*       2.1 computes contravariant components
!
IF (HUVW_ADV_SCHEME=='CEN2ND' ) THEN
 CALL CONTRAV (HLBCX,HLBCY,PUT,PVT,PWT,PDXX,PDYY,PDZZ,PDZX,PDZY,ZRUCPPM,ZRVCPPM,ZRWCPPM,2)
ELSE
 CALL CONTRAV (HLBCX,HLBCY,PUT,PVT,PWT,PDXX,PDYY,PDZZ,PDZX,PDZY,ZRUCPPM,ZRVCPPM,ZRWCPPM,4)
END IF
!
!
!*       2.2 computes CFL numbers
!

IF (.NOT. L1D) THEN
  ZCFLU = 0.0 ; ZCFLV = 0.0 ;  ZCFLW = 0.0
  ZCFLU(IIB:IIE,IJB:IJE,:) = ABS(ZRUCPPM(IIB:IIE,IJB:IJE,:) * PTSTEP)
  ZCFLV(IIB:IIE,IJB:IJE,:) = ABS(ZRVCPPM(IIB:IIE,IJB:IJE,:) * PTSTEP)
  ZCFLW(IIB:IIE,IJB:IJE,:) = ABS(ZRWCPPM(IIB:IIE,IJB:IJE,:) * PTSTEP)
   IF (LIBM) THEN
    ZCFLU(IIB:IIE,IJB:IJE,:) = ZCFLU(IIB:IIE,IJB:IJE,:)*(1.-exp(-(XIBM_LS(IIB:IIE,IJB:IJE,:,2)/&
                                                        (XRHODJ(IIB:IIE,IJB:IJE,:)/XRHODREF(IIB:IIE,IJB:IJE,:))**(1./3.))**2.))
    ZCFLV(IIB:IIE,IJB:IJE,:) = ZCFLV(IIB:IIE,IJB:IJE,:)*(1.-exp(-(XIBM_LS(IIB:IIE,IJB:IJE,:,3)/&
                                                        (XRHODJ(IIB:IIE,IJB:IJE,:)/XRHODREF(IIB:IIE,IJB:IJE,:))**(1./3.))**2.))
    ZCFLW(IIB:IIE,IJB:IJE,:) = ZCFLW(IIB:IIE,IJB:IJE,:)*(1.-exp(-(XIBM_LS(IIB:IIE,IJB:IJE,:,4)/&
                                                        (XRHODJ(IIB:IIE,IJB:IJE,:)/XRHODREF(IIB:IIE,IJB:IJE,:))**(1./3.))**2.))
    WHERE (XIBM_LS(IIB:IIE,IJB:IJE,:,2).GT.(-XIBM_EPSI)) ZCFLU(IIB:IIE,IJB:IJE,:)=0.
    WHERE (XIBM_LS(IIB:IIE,IJB:IJE,:,3).GT.(-XIBM_EPSI)) ZCFLV(IIB:IIE,IJB:IJE,:)=0.
    WHERE (XIBM_LS(IIB:IIE,IJB:IJE,:,4).GT.(-XIBM_EPSI)) ZCFLW(IIB:IIE,IJB:IJE,:)=0.
  ENDIF 
  IF (.NOT. L2D) THEN
    ZCFL  = SQRT(ZCFLU**2+ZCFLV**2+ZCFLW**2)
  ELSE
    ZCFL  = SQRT(ZCFLU**2+ZCFLW**2)
  END IF
ELSE
   ZCFLU = 0.0 ; ZCFLV = 0.0 ;  ZCFLW = 0.0 
   ZCFLW(IIB:IIE,IJB:IJE,:) = ABS(ZRWCPPM(IIB:IIE,IJB:IJE,:) * PTSTEP)
   ZCFL = SQRT(ZCFLW**2)
END IF
!
!* prints in the file the 3D Courant numbers (one should flag this)
!
IF ( tpfile%lopened .AND. OCFL_WRIT .AND. (.NOT. L1D) ) THEN
    TZFIELD = TFIELDMETADATA(    &
      CMNHNAME   = 'CFLU',       &
      CSTDNAME   = '',           &
      CLONGNAME  = 'CFLU',       &
      CUNITS     = '1',          &
      CDIR       = 'XY',         &
      CCOMMENT   = 'X_Y_Z_CFLU', &
      NGRID      = 1,            &
      NTYPE      = TYPEREAL,     &
      NDIMS      = 3,            &
      LTIMEDEP   = .TRUE.        )
    CALL IO_Field_write(TPFILE,TZFIELD,ZCFLU)
!
  IF (.NOT. L2D) THEN
    TZFIELD = TFIELDMETADATA(    &
      CMNHNAME   = 'CFLV',       &
      CSTDNAME   = '',           &
      CLONGNAME  = 'CFLV',       &
      CUNITS     = '1',          &
      CDIR       = 'XY',         &
      CCOMMENT   = 'X_Y_Z_CFLV', &
      NGRID      = 1,            &
      NTYPE      = TYPEREAL,     &
      NDIMS      = 3,            &
      LTIMEDEP   = .TRUE.        )
    CALL IO_Field_write(TPFILE,TZFIELD,ZCFLV)
  END IF
!
    TZFIELD = TFIELDMETADATA(    &
      CMNHNAME   = 'CFLW',       &
      CSTDNAME   = '',           &
      CLONGNAME  = 'CFLW',       &
      CUNITS     = '1',          &
      CDIR       = 'XY',         &
      CCOMMENT   = 'X_Y_Z_CFLW', &
      NGRID      = 1,            &
      NTYPE      = TYPEREAL,     &
      NDIMS      = 3,            &
      LTIMEDEP   = .TRUE.        )
    CALL IO_Field_write(TPFILE,TZFIELD,ZCFLW)
!
    TZFIELD = TFIELDMETADATA(   &
      CMNHNAME   = 'CFL',       &
      CSTDNAME   = '',          &
      CLONGNAME  = 'CFL',       &
      CUNITS     = '1',         &
      CDIR       = 'XY',        &
      CCOMMENT   = 'X_Y_Z_CFL', &
      NGRID      = 1,           &
      NTYPE      = TYPEREAL,    &
      NDIMS      = 3,           &
      LTIMEDEP   = .TRUE.       )
    CALL IO_Field_write(TPFILE,TZFIELD,ZCFL)
END IF
!
!* prints in the output file the maximum CFL
!
ZCFLU_MAX = MAX_ll(ZCFLU,IINFO_ll)
ZCFLV_MAX = MAX_ll(ZCFLV,IINFO_ll)
ZCFLW_MAX = MAX_ll(ZCFLW,IINFO_ll)
ZCFL_MAX  = MAX_ll(ZCFL,IINFO_ll)
!
WRITE(ILUOUT,FMT='(A24,F10.2,A5,F10.2,A5,F10.2,A9,F10.2)') &
                'Max. CFL number for U : ',ZCFLU_MAX,  &
                '  V : ',ZCFLV_MAX,'  W : ', ZCFLW_MAX,&
                'global : ',ZCFL_MAX
!
!
!*       2.3 updates time step splitting loop
!
IF (OSPLIT_CFL .AND. (.NOT.L1D)  ) THEN
!
 ISPLIT_PPM = INT(ZCFL_MAX/PSPLIT_CFL)+1
 IF ( KSPLIT /= ISPLIT_PPM )                                    &
 WRITE(ILUOUT,FMT='(A37,I2,A4,I2,A11)')                         &
                  'PPM  time spliting loop changed from ',      &
                  KSPLIT,' to ',ISPLIT_PPM, ' iterations'
!
 KSPLIT =     ISPLIT_PPM                      
!
END IF
! ---------------------------------------------------------------
IF (( (ZCFLU_MAX>=3.) .AND. (.NOT.L1D) ) .OR. &
    ( (ZCFLV_MAX>=3.) .AND. (.NOT.L1D) .AND. (.NOT.L2D) ) .OR. &
    ( (ZCFLW_MAX>=8.) .AND. (.NOT.L1D) ) ) THEN
  WRITE(ILUOUT,*) ' '
  WRITE(ILUOUT,*) ' +---------------------------------------------------+'
  WRITE(ILUOUT,*) ' |                   MODEL ERROR                     |'
  WRITE(ILUOUT,*) ' +---------------------------------------------------+'
  WRITE(ILUOUT,*) ' |                                                   |'
  WRITE(ILUOUT,*) ' |       The model wind speed becomes too high       |'
  WRITE(ILUOUT,*) ' |                                                   |'
  IF ( ZCFLU_MAX>=3. .OR. ZCFLV_MAX>=3. ) &
  WRITE(ILUOUT,*) ' |    The  horizontal CFL value reaches 3. or more   |'
  IF ( ZCFLW_MAX>=8.                    ) &
  WRITE(ILUOUT,*) ' |    The  vertical   CFL value reaches 8. or more   |'
  WRITE(ILUOUT,*) ' |                                                   |'
  WRITE(ILUOUT,*) ' |    This can be due either to :                    |'
  WRITE(ILUOUT,*) ' |     - a numerical explosion of the model          |'
  WRITE(ILUOUT,*) ' |     - or a too high wind speed for an             |'
  WRITE(ILUOUT,*) ' |       acceptable accuracy of the advection        |'
  WRITE(ILUOUT,*) ' |                                                   |'
  WRITE(ILUOUT,*) ' |        Please decrease your time-step             |'
  WRITE(ILUOUT,*) ' |                                                   |'
  WRITE(ILUOUT,*) ' +---------------------------------------------------+'
  WRITE(ILUOUT,*) ' '
  WRITE(ILUOUT,*) ' +---------------------------------------------------+'
  WRITE(ILUOUT,*) ' |                   MODEL STOPS                     |'
  WRITE(ILUOUT,*) ' +---------------------------------------------------+'
  CALL PRINT_MSG(NVERB_FATAL,'GEN','ADVECTION_METSV','')
END IF
!
!
ZTSTEP_PPM = PTSTEP / REAL(KSPLIT)
!
!
!*      2.4 normalized contravariant components for split PPM time-step
!
ZRUCPPM = ZRUCPPM*ZTSTEP_PPM
ZRVCPPM = ZRVCPPM*ZTSTEP_PPM
ZRWCPPM = ZRWCPPM*ZTSTEP_PPM
!
!
!-------------------------------------------------------------------------------
!
!
!*       3.     COMPUTES THE TENDENCIES SINCE THE BEGINNING OF THE TIME STEP
!	        ------------------------------------------------------------
!
!* This represent the effects of all OTHER processes
!  Clouds    related processes from previous time-step are     taken into account in PRTHS_CLD
!  Advection related processes from previous time-step will be taken into account in ZRTHS_PPM
!
ZRTHS_OTHER = PRTHS - PTHT * PRHODJ / PTSTEP
IF (GTKE) ZRTKES_OTHER = PRTKES - PTKET * PRHODJ / PTSTEP                      
DO JR = 1, KRR
 ZRRS_OTHER(:,:,:,JR) = PRRS(:,:,:,JR) - PRT(:,:,:,JR) * PRHODJ(:,:,:) / PTSTEP
END DO
DO JSV = 1, KSV
 ZRSVS_OTHER(:,:,:,JSV) = PRSVS(:,:,:,JSV) - PSVT(:,:,:,JSV) * PRHODJ / PTSTEP
END DO
IF(LBLOWSNOW) THEN
   DO JSV = 1, (NBLOWSNOW_2D)        
     ZRSNWCS_OTHER(:,:,:,JSV) = ZRSNWCS(:,:,:,JSV) - ZSNWC_INIT(:,:,:,JSV) * PRHODJ / PTSTEP
   END DO   
ENDIF
!
! Top and bottom Boundaries 
!
CALL ADV_BOUNDARIES (HLBCX, HLBCY, ZRTHS_OTHER)
IF (GTKE) CALL ADV_BOUNDARIES (HLBCX, HLBCY, ZRTKES_OTHER)
DO JR = 1, KRR
  CALL ADV_BOUNDARIES (HLBCX, HLBCY, ZRRS_OTHER(:,:,:,JR))
END DO
DO JSV = 1, KSV
  CALL ADV_BOUNDARIES (HLBCX, HLBCY, ZRSVS_OTHER(:,:,:,JSV))
END DO
IF(LBLOWSNOW) THEN
  DO JSV = 1, (NBLOWSNOW_2D)        
    CALL ADV_BOUNDARIES (HLBCX, HLBCY, ZRSNWCS_OTHER(:,:,:,JSV))
  END DO 
END IF
!
! Exchanges on processors
!
NULLIFY(TZFIELDS0_ll)
!!$IF(NHALO == 1) THEN
  CALL ADD3DFIELD_ll( TZFIELDS0_ll, ZRTHS_OTHER, 'ADVECTION_METSV::ZRTHS_OTHER' )
  IF (GTKE) CALL ADD3DFIELD_ll( TZFIELDS0_ll, ZRTKES_OTHER, 'ADVECTION_METSV::ZRTKES_OTHER' )
  IF ( KRR>0 )  CALL ADD4DFIELD_ll( TZFIELDS0_ll, ZRRS_OTHER(:,:,:,1:KRR),             'ADVECTION_METSV::ZRRS_OTHER'    )
  IF ( KSV>0 )  CALL ADD4DFIELD_ll( TZFIELDS0_ll, ZRSVS_OTHER(:,:,:,1:KSV),            'ADVECTION_METSV::ZRSVS_OTHER'   )
  IF(LBLOWSNOW) CALL ADD4DFIELD_ll( TZFIELDS0_ll, ZRSNWCS_OTHER(:,:,:,1:NBLOWSNOW_2D), 'ADVECTION_METSV::ZRSNWCS_OTHER' )
  CALL UPDATE_HALO_ll(TZFIELDS0_ll,IINFO_ll)
  CALL CLEANLIST_ll(TZFIELDS0_ll)
!!$END IF
!
!

!-------------------------------------------------------------------------------
!
!*       4.     CALLS THE PPM ADVECTION INSIDE A TIME SPLITTING         
!	        --------------------------------------
!
CALL PPM_RHODJ(HLBCX,HLBCY, ZRUCPPM, ZRVCPPM, ZRWCPPM,              &
               ZTSTEP_PPM, PRHODJ, ZRHOX1, ZRHOX2, ZRHOY1, ZRHOY2,  &
               ZRHOZ1, ZRHOZ2                                       )
!
!* values of the fields at the beginning of the time splitting loop
ZTH   = PTHT
ZTKE   = PTKET
IF (KRR /=0 ) ZR    = PRT
IF (KSV /=0 ) ZSV   = PSVT
IF(LBLOWSNOW) THEN
    DO JSV = 1, (NBLOWSNOW_2D)
        ZSNWC(:,:,:,JSV) = ZRSNWCS(:,:,:,JSV)* PTSTEP/ PRHODJ
        CALL ADV_BOUNDARIES (HLBCX, HLBCY, ZSNWC(:,:,:,JSV))
    END DO
    ZSNWC_INIT=ZSNWC
ENDIF
!
IF (GTKE)    PRTKES_ADV(:,:,:)  = 0.              
!
!* time splitting loop
DO JSPL=1,KSPLIT
!
   !ZRTHS_PPM(:,:,:)   = 0.
   !ZRTKES_PPM(:,:,:)   = 0.
   !IF (KRR /=0) ZRRS_PPM(:,:,:,:)   = 0.
   !IF (KSV /=0) ZRSVS_PPM(:,:,:,:)   = 0.
!
   IF (LNEUTRAL) ZTH=ZTH-PTHVREF  !* To be removed with the new PPM scheme ?
   CALL PPM_MET (HLBCX,HLBCY, KRR, TPDTCUR,ZRUCPPM, ZRVCPPM, ZRWCPPM, PTSTEP,ZTSTEP_PPM, &
              PRHODJ,  ZRHOX1, ZRHOX2, ZRHOY1, ZRHOY2,  ZRHOZ1, ZRHOZ2,                  &
              ZTH, ZTKE, ZR, ZRTHS_PPM, ZRTKES_PPM, ZRRS_PPM, HMET_ADV_SCHEME)
   IF (LNEUTRAL) ZTH=ZTH+PTHVREF  !* To be removed with the new PPM scheme ?
!
   CALL PPM_SCALAR (HLBCX,HLBCY, KSV, TPDTCUR, ZRUCPPM, ZRVCPPM, ZRWCPPM, PTSTEP,     &
                 ZTSTEP_PPM, PRHODJ, ZRHOX1, ZRHOX2, ZRHOY1, ZRHOY2,  ZRHOZ1, ZRHOZ2, &
                 ZSV, ZRSVS_PPM, HSV_ADV_SCHEME                                       )
!
! Tendencies of PPM
!
   PRTHS(:,:,:)                      = PRTHS     (:,:,:)   + ZRTHS_PPM (:,:,:)   / KSPLIT
   IF (GTKE)     PRTKES_ADV(:,:,:)   = PRTKES_ADV(:,:,:)   + ZRTKES_PPM(:,:,:)   / KSPLIT
   IF (KRR /=0)  PRRS      (:,:,:,:) = PRRS      (:,:,:,:) + ZRRS_PPM  (:,:,:,:) / KSPLIT
   IF (KSV /=0 ) PRSVS     (:,:,:,:) = PRSVS     (:,:,:,:) + ZRSVS_PPM (:,:,:,:) / KSPLIT
!
   IF (JSPL<KSPLIT) THEN
!
!  Guesses of the field inside the time splitting loop
!
   ZTH = ZTH + ( ZRTHS_PPM(:,:,:) + ZRTHS_OTHER(:,:,:) + PRTHS_CLD(:,:,:)) * &
           ZTSTEP_PPM / PRHODJ(:,:,:)
   IF (GTKE) ZTKE = ZTKE + ( ZRTKES_PPM(:,:,:) + ZRTKES_OTHER(:,:,:) ) * ZTSTEP_PPM / PRHODJ(:,:,:)
   DO JR = 1, KRR
    ZR(:,:,:,JR) = ZR(:,:,:,JR) + ( ZRRS_PPM(:,:,:,JR) + ZRRS_OTHER(:,:,:,JR) + PRRS_CLD(:,:,:,JR) ) &
                    * ZTSTEP_PPM / PRHODJ(:,:,:)
   END DO
   DO JSV = 1, KSV
    ZSV(:,:,:,JSV) = ZSV(:,:,:,JSV) + ( ZRSVS_PPM(:,:,:,JSV) + ZRSVS_OTHER(:,:,:,JSV) +  &
                     PRSVS_CLD(:,:,:,JSV) ) * ZTSTEP_PPM / PRHODJ(:,:,:)
   END DO
!
! Top and bottom Boundaries and LBC for the guesses
!
   CALL ADV_BOUNDARIES (HLBCX, HLBCY, ZTH, PTHT )    
    IF (GTKE) CALL ADV_BOUNDARIES (HLBCX, HLBCY, ZTKE, PTKET)
   DO JR = 1, KRR
     CALL ADV_BOUNDARIES (HLBCX, HLBCY, ZR(:,:,:,JR), PRT(:,:,:,JR))
   END DO
   DO JSV = 1, KSV
     CALL ADV_BOUNDARIES (HLBCX, HLBCY, ZSV(:,:,:,JSV), PSVT(:,:,:,JSV))
   END DO

   IF(LBLOWSNOW) THEN ! Advection of Canopy mass at the 1st atmospheric level
      ZRSNWCS_PPM(:,:,:,:) = 0.
   !

      CALL PPM_SCALAR (HLBCX,HLBCY, NBLOWSNOW_2D, TPDTCUR, ZRUCPPM, ZRVCPPM, ZRWCPPM,PTSTEP,    &
                 ZTSTEP_PPM, PRHODJ, ZRHOX1, ZRHOX2, ZRHOY1, ZRHOY2,  ZRHOZ1, ZRHOZ2,          &
                 ZSNWC, ZRSNWCS_PPM, HSV_ADV_SCHEME)


! Tendencies of PPM
      ZRSNWCS(:,:,:,:)        =    ZRSNWCS(:,:,:,:)  + ZRSNWCS_PPM (:,:,:,:)   / KSPLIT
!  Guesses of the field inside the time splitting loop 
      DO JSV = 1, ( NBLOWSNOW_2D)
          ZSNWC(:,:,:,JSV) = ZSNWC(:,:,:,JSV) + ZRSNWCS_PPM(:,:,:,JSV)*ZTSTEP_PPM/ PRHODJ(:,:,:)  
      END DO

! Top and bottom Boundaries and LBC for the guesses      
      DO JSV = 1, (NBLOWSNOW_2D)
          CALL ADV_BOUNDARIES (HLBCX, HLBCY, ZSNWC(:,:,:,JSV), ZSNWC_INIT(:,:,:,JSV))
      END DO
   END IF   
!
!  Exchanges fields between processors
!
   NULLIFY(TZFIELDS1_ll)
!!$   IF(NHALO == 1) THEN
    CALL ADD3DFIELD_ll( TZFIELDS1_ll, ZTH, 'ZTH' )
    IF (GTKE)        CALL ADD3DFIELD_ll( TZFIELDS1_ll, ZTKE,                        'ADVECTION_METSV::ZTKE' )
    IF ( KRR>0 )     CALL ADD4DFIELD_ll( TZFIELDS1_ll, ZR (:,:,:,1:KRR),            'ADVECTION_METSV::ZR'    )
    IF ( KSV>0 )     CALL ADD4DFIELD_ll( TZFIELDS1_ll, ZSV(:,:,:,1:KSV),            'ADVECTION_METSV::ZSV'   )
    IF ( LBLOWSNOW ) CALL ADD4DFIELD_ll( TZFIELDS1_ll, ZSNWC(:,:,:,1:NBLOWSNOW_2D), 'ADVECTION_METSV::ZSNWC' )
    CALL UPDATE_HALO_ll(TZFIELDS1_ll,IINFO_ll)
    CALL CLEANLIST_ll(TZFIELDS1_ll)
!!$   END IF
   END IF
!
END DO
!
!-------------------------------------------------------------------------------
!
!  TKE special case: advection is the last process for TKE
!
! TKE must be greater than its minimum value
! (previously done in tke_eps_sources)
!
IF (GTKE) THEN
   PRTKES(:,:,:)  = PRTKES(:,:,:) + PRTKES_ADV(:,:,:)
   PRTKES(:,:,:) = MAX (PRTKES(:,:,:) , XTKEMIN * PRHODJ(:,:,:) / PTSTEP )
END IF
!
!
!-------------------------------------------------------------------------------
! Update tendency for cano variables : from 3D to 2D
! 
IF(LBLOWSNOW) THEN

    DO JSV=1,(NBLOWSNOW_2D) 
       DO JI=1,SIZE(PSVT,1)
         DO JJ=1,SIZE(PSVT,2)
             XRSNWCANOS(JI,JJ,JSV) = SUM(ZRSNWCS(JI,JJ,IKB:IKE,JSV))
         END DO
       END DO
    END DO
IF(LWEST_ll())  XRSNWCANOS(IIB,:,:)  = ZRSNWCS(IIB,:,IKB,:)
IF(LEAST_ll())  XRSNWCANOS(IIE,:,:)  = ZRSNWCS(IIE,:,IKB,:)
IF(LSOUTH_ll()) XRSNWCANOS(:,IJB,:)  = ZRSNWCS(:,IJB,IKB,:)
IF(LNORTH_ll()) XRSNWCANOS(:,IJE,:)  = ZRSNWCS(:,IJE,IKB,:)

END IF
!-------------------------------------------------------------------------------
!
!*       5.     BUDGETS                                                 
!	        -------
!
if ( lbudget_th  ) call Budget_store_end( tbudgets(NBUDGET_TH ), 'ADV', prths (:, :, :)    )
if ( lbudget_tke ) call Budget_store_end( tbudgets(NBUDGET_TKE), 'ADV', prtkes(:, :, :)    )
if ( lbudget_rv  ) call Budget_store_end( tbudgets(NBUDGET_RV ), 'ADV', prrs  (:, :, :, 1) )
if ( lbudget_rc  ) call Budget_store_end( tbudgets(NBUDGET_RC ), 'ADV', prrs  (:, :, :, 2) )
if ( lbudget_rr  ) call Budget_store_end( tbudgets(NBUDGET_RR ), 'ADV', prrs  (:, :, :, 3) )
if ( lbudget_ri  ) call Budget_store_end( tbudgets(NBUDGET_RI ), 'ADV', prrs  (:, :, :, 4) )
if ( lbudget_rs  ) call Budget_store_end( tbudgets(NBUDGET_RS ), 'ADV', prrs  (:, :, :, 5) )
if ( lbudget_rg  ) call Budget_store_end( tbudgets(NBUDGET_RG ), 'ADV', prrs  (:, :, :, 6) )
if ( lbudget_rh  ) call Budget_store_end( tbudgets(NBUDGET_RH ), 'ADV', prrs  (:, :, :, 7) )
if ( lbudget_sv) then
  do jsv = 1, ksv
    call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + jsv ), 'ADV', prsvs(:, :, :, jsv) )
  end do
end if

! Remove non-physical negative values (unnecessary in a perfect world) + corresponding budgets
call Sources_neg_correct( hcloud, 'NEADV', krr, ptstep, ppabst, ptht, prt, prths, prrs, prsvs )

!-------------------------------------------------------------------------------
!
END SUBROUTINE ADVECTION_METSV
