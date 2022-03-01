!     ######spl
      SUBROUTINE INI_BUDGET(KLUOUT, HLUOUT,PTSTEP,KSV,KRR,            &
      ONUMDIFU,ONUMDIFTH,ONUMDIFSV,                                   &
      OHORELAX_UVWTH,OHORELAX_RV,OHORELAX_RC,OHORELAX_RR,             &
      OHORELAX_RI,OHORELAX_RS, OHORELAX_RG, OHORELAX_RH,OHORELAX_TKE, &
      OHORELAX_SV,OVE_RELAX,OCHTRANS,ONUDGING,ODRAGTREE,              &
      HRAD,HDCONV,HSCONV,HTURB,HTURBDIM,HCLOUD,                       &
      HMET_ADV_SCHEME,HSV_ADV_SCHEME)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!     #################################################################
!
!!****  *INI_BUDGET* - routine to initialize the parameters for the budgets
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to set or compute the parameters used
!     by the MESONH budgets. Names of files for budget recording are processed
!     and storage arrays are initialized.
!
!!**  METHOD
!!    ------
!!      The essential of information is passed by modules. The choice of budgets
!!    and processes set by the user as integers is converted in "actions"
!!    readable  by the subroutine BUDGET under the form of string characters.
!!    For each complete process composed of several elementary processes, names
!!    of elementary processes are concatenated in order to have an explicit name
!!    in the comment of the recording file for budget.
!!
!!
!!    EXTERNAL
!!    --------
!!      None
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_PARAMETERS: JPBUMAX,JPBUPROMAX
!!
!!      Module MODD_CONF: CCONF
!!
!!      Module MODD_DYN:  XSEGLEN
!!
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation (routine INI_BUDGET)
!!
!!
!!    AUTHOR
!!    ------
!!      P. Hereil      * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original        01/03/95
!!      J. Stein        25/06/95  put the sources in phase with the code
!!      J. Stein        20/07/95  reset to FALSE of all the switches when
!!                                CBUTYPE /= MASK or CART
!!      J. Stein        26/06/96  add the new sources + add the increment between
!!                                2 active processes
!!      J.-P. Pinty     13/12/96  Allowance of multiple SVs
!!      J.-P. Pinty     11/01/97  Includes deep convection ice and forcing processes
!!      J.-P. Lafore    10/02/98  Allocation of the RHODJs for budget
!!      V. Ducrocq      04/06/99  //
!!      N. Asencio      18/06/99  // MASK case : delete KIMAX and KJMAX arguments,
!!                                GET_DIM_EXT_ll initializes the dimensions of the
!!                                extended local domain.
!!                                LBU_MASK and XBUSURF are allocated on the extended
!!                                local domain.
!!                                add 3 local variables IBUDIM1,IBUDIM2,IBUDIM3
!!                                to define the dimensions of the budget arrays
!!                                in the different cases CART and MASK
!!      J.-P. Pinty     23/09/00  add budget for C2R2
!!      V. Masson       18/11/02  add budget for 2way nesting
!!      O.Geoffroy      03/2006   Add KHKO scheme
!!      J.-P. Pinty     22/04/97  add the explicit hail processes
!!      C.Lac           10/08/07  Add ADV for PPM without contribution
!!                                of each direction
!!      C. Barthe       19/11/09  Add atmospheric electricity
!!      C.Lac           01/07/11  Add vegetation drag
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAMETERS
USE MODD_BUDGET
USE MODD_DYN
USE MODD_CONF
USE MODD_PARAM_ICE
USE MODD_PARAM_C2R2
USE MODD_ELEC_DESCR, ONLY : LINDUCTIVE
!
!USE MODE_ll
!USE MODE_IO_ll
!
IMPLICIT NONE
!
!*       0.1   declarations of argument
!
!
INTEGER,         INTENT(IN) :: KLUOUT   ! Logical unit number for prints
CHARACTER (LEN=*), INTENT(IN) :: HLUOUT ! name of output listing
REAL, INTENT(IN) :: PTSTEP              ! time step
INTEGER, INTENT(IN) :: KSV              ! number of scalar variables
INTEGER, INTENT(IN) :: KRR              ! number of moist variables
LOGICAL, INTENT(IN) :: ONUMDIFU         ! switch to activate the numerical
                                        ! diffusion for momentum
LOGICAL, INTENT(IN) :: ONUMDIFTH        ! for meteorological scalar variables
LOGICAL, INTENT(IN) :: ONUMDIFSV        ! for tracer scalar variables
LOGICAL, INTENT(IN) :: OHORELAX_UVWTH  ! switch for the
                       ! horizontal relaxation for U,V,W,TH
LOGICAL, INTENT(IN) :: OHORELAX_RV     ! switch for the
                       ! horizontal relaxation for Rv
LOGICAL, INTENT(IN) :: OHORELAX_RC     ! switch for the
                       ! horizontal relaxation for Rc
LOGICAL, INTENT(IN) :: OHORELAX_RR     ! switch for the
                       ! horizontal relaxation for Rr
LOGICAL, INTENT(IN) :: OHORELAX_RI     ! switch for the
                       ! horizontal relaxation for Ri
LOGICAL, INTENT(IN) :: OHORELAX_RS     ! switch for the
                       ! horizontal relaxation for Rs
LOGICAL, INTENT(IN) :: OHORELAX_RG     ! switch for the
                       ! horizontal relaxation for Rg
LOGICAL, INTENT(IN) :: OHORELAX_RH     ! switch for the
                       ! horizontal relaxation for Rh
LOGICAL, INTENT(IN) :: OHORELAX_TKE    ! switch for the
                       ! horizontal relaxation for tke
LOGICAL,DIMENSION(:),INTENT(IN):: OHORELAX_SV     ! switch for the
                       ! horizontal relaxation for scalar variables
LOGICAL, INTENT(IN) :: OVE_RELAX        ! switch to activate the vertical
                                        ! relaxation
LOGICAL, INTENT(IN) :: OCHTRANS         ! switch to activate convective
                                        !transport for SV
LOGICAL, INTENT(IN) :: ONUDGING         ! switch to activate nudging
LOGICAL, INTENT(IN) :: ODRAGTREE        ! switch to activate vegetation drag
CHARACTER (LEN=*), INTENT(IN) :: HRAD   ! type of the radiation scheme
CHARACTER (LEN=*), INTENT(IN) :: HDCONV ! type of the deep convection scheme
CHARACTER (LEN=*), INTENT(IN) :: HSCONV ! type of the shallow convection scheme
CHARACTER (LEN=*), INTENT(IN) :: HTURB  ! type of the turbulence scheme
CHARACTER (LEN=*), INTENT(IN) :: HTURBDIM! dimensionnality of the turbulence
                                        ! scheme
CHARACTER (LEN=*), INTENT(IN) :: HCLOUD ! type of microphysical scheme
CHARACTER (LEN=*), INTENT(IN) :: HMET_ADV_SCHEME ! type of advection scheme
                                        ! for meteorological scalar variables
CHARACTER (LEN=*), INTENT(IN) :: HSV_ADV_SCHEME ! type of advection scheme
                                        ! for tracer scalar variables
!
!*       0.2   declarations of local variables
!
INTEGER, DIMENSION(JPBUMAX,JPBUPROMAX+1) :: IPROACTV      ! switches set by the
                                                          ! user for process
                                                          ! activation
INTEGER :: JI, JJ, JK , JJJ                               ! loop indices
INTEGER :: IIMAX_ll, IJMAX_ll ! size of the physical global domain
INTEGER :: ITEN                                           ! tens for CBURECORD
INTEGER :: IPROC                                          ! counter for processes
INTEGER :: IIU, IJU                                       ! size along x and y directions
                                                          ! of the extended subdomain
INTEGER :: IBUDIM1                                        ! first dimension of the budget arrays
                                                          ! = NBUIMAX in CART case
                                                          ! = NBUKMAX in MASK case
INTEGER :: IBUDIM2                                        ! second dimension of the budget arrays
                                                          ! = NBUJMAX in CART case
                                                          ! = NBUWRNB in MASK case
INTEGER :: IBUDIM3                                        ! third dimension of the budget arrays
                                                          ! = NBUKMAX in CART case
                                                          ! = NBUMASK in MASK case
LOGICAL :: GERROR                                         ! switch for error in
                                                          ! budget specifcation
CHARACTER(LEN=7), DIMENSION(JPBUMAX) :: YEND_COMMENT      ! last part of comment
                                                          ! for budgets records
CHARACTER(LEN=6), DIMENSION(JPBUMAX,JPBUPROMAX) :: YWORK2 ! used for
                                                          ! concatenattion of
                                                          ! comments for budgets
CHARACTER(LEN=40)                  :: YSTRING
INTEGER                            :: ILEN
INTEGER :: JSV               ! loop indice for the SVs
INTEGER :: IBUPROCNBR_SV_MAX ! Max number of processes for the SVs
INTEGER :: ILAST_PROC_NBR    ! Index of the last process number
INTEGER :: IINFO_ll ! return status of the interface routine
INTEGER :: IRESP   ! Return code of FM-routines
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! the lines below must be update as soon as MODD_BUDGET is updated
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-------------------------------------------------------------------------------
!
!*       1.    COMPUTE BUDGET VARIABLES
!              ------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('INI_BUDGET',0,ZHOOK_HANDLE)
NBUSTEP = NINT (XBULEN / PTSTEP)
NBUTSHIFT=0
!
!  common dimension for all CBUTYPE values
!
IF (LBU_KCP) THEN
  NBUKMAX = 1
ELSE
  NBUKMAX = NBUKH - NBUKL +1
END IF
!
IF (CBUTYPE=='CART') THEN              ! cartesian case only
!
  NBUWRNB = 1        ! after every budget period, we write the result on the
                     ! FM_FILE
! Following lines are not used in AROME
!  IF (LBU_ICP) THEN
!    NBUIMAX_ll = 1
!  ELSE
!    NBUIMAX_ll = NBUIH - NBUIL +1
!  END IF
!  IF (LBU_JCP) THEN
!    NBUJMAX_ll = 1
!  ELSE
!    NBUJMAX_ll = NBUJH - NBUJL +1
!  END IF
!
!  CALL GET_INTERSECTION_ll(NBUIL,NBUJL,NBUIH,NBUJH, &
!      NBUSIL,NBUSJL,NBUSIH,NBUSJH,"PHYS",IINFO_ll)
! Call to GET_INTERSECTION_ll is replaced by:
  NBUSIL=NBUIL
  NBUSJL=NBUJL
  NBUSIH=NBUIH
  NBUSJH=NBUJH
  IINFO_ll=0
  IF ( IINFO_ll /= 1 ) THEN !
    IF (LBU_ICP) THEN
      NBUIMAX = 1
    ELSE
      NBUIMAX = NBUSIH - NBUSIL +1
    END IF
    IF (LBU_JCP) THEN
      NBUJMAX = 1
    ELSE
      NBUJMAX =  NBUSJH - NBUSJL +1
    END IF
  ELSE ! the intersection is void
    CBUTYPE='SKIP'  ! no budget on this processor
    NBUIMAX = 0     ! in order to allocate void arrays
    NBUJMAX = 0
  ENDIF
! three first dimensions of budget arrays in cart and skip cases
   IBUDIM1=NBUIMAX
   IBUDIM2=NBUJMAX
   IBUDIM3=NBUKMAX
! these variables are not be used
   NBUWRNB=-1
   NBUMASK=-1
!
ELSEIF (CBUTYPE=='MASK') THEN          ! mask case only
!
  LBU_ENABLE=.TRUE.
  NBUWRNB = NINT (XBUWRI / XBULEN)  ! only after NBUWRNB budget periods, we write the
                                    ! result on the FM_FILE
  NBUTIME = 1

!  CALL GET_DIM_EXT_ll ('B', IIU,IJU)
  ALLOCATE( LBU_MASK( IIU ,IJU, NBUMASK) )
  LBU_MASK(:,:,:)=.FALSE.
  ALLOCATE( XBUSURF( IIU, IJU, NBUMASK, NBUWRNB) )
  XBUSURF(:,:,:,:) = 0.
!
! three first dimensions of budget arrays in mask case
!  the order of the dimensions are the order expected in WRITE_DIACHRO routine:
!  x,y,z,time,mask,processus  and in this case x and y are missing
!  first dimension of the arrays : dimension along K
!  second dimension of the arrays : number of the budget time period
!  third dimension of the arrays : number of the budget masks zones
  IBUDIM1=NBUKMAX
  IBUDIM2=NBUWRNB
  IBUDIM3=NBUMASK
! these variables are not used in this case
  NBUIMAX=-1
  NBUJMAX=-1
! the beginning and the end along x and y direction : global extended domain
 ! get dimensions of the physical global domain
 !  CALL GET_GLOBALDIMS_ll (IIMAX_ll,IJMAX_ll)
   NBUIL=1
   NBUIH=IIMAX_ll + 2 * JPHEXT
   NBUJL=1
   NBUJH=IJMAX_ll + 2 * JPHEXT

!
ELSE                      ! default case
!
  LBU_ENABLE=.FALSE.
  NBUIMAX = -1
  NBUJMAX = -1
  LBU_RU = .FALSE.
  LBU_RV = .FALSE.
  LBU_RW = .FALSE.
  LBU_RTH= .FALSE.
  LBU_RTKE= .FALSE.
  LBU_RRV= .FALSE.
  LBU_RRC= .FALSE.
  LBU_RRR= .FALSE.
  LBU_RRI= .FALSE.
  LBU_RRS= .FALSE.
  LBU_RRG= .FALSE.
  LBU_RRH= .FALSE.
  LBU_RSV= .FALSE.
!
! three first dimensions of budget arrays in default case
  IBUDIM1=0
  IBUDIM2=0
  IBUDIM3=0
!
END IF
!
!
!-------------------------------------------------------------------------------
!
!*       2.    ALLOCATE MEMORY FOR BUDGET ARRAYS AND INITIALIZE
!              ------------------------------------------------
!
ALLOCATE( NBUPROCNBR(JPBUMAX) )
ALLOCATE( NBUPROCCTR(JPBUMAX) )
ALLOCATE( CBUACTION(JPBUMAX, JPBUPROMAX) )
ALLOCATE( CBUCOMMENT(JPBUMAX, JPBUPROMAX) )
ALLOCATE( CBURECORD(JPBUMAX, JPBUPROMAX) )
NBUPROCCTR(:) = 0
NBUCTR_ACTV(:) = 0
NBUPROCNBR(:) = 0
CBUACTION(:,:) = 'OF'
CBURECORD(:,:) = ' '
CBUCOMMENT(:,:) = ' '
LBU_BEG =.TRUE.
!
!-------------------------------------------------------------------------------
!
!*       3.    INITALIZE VARIABLES
!              -------------------
!
IPROACTV(:,:) = 3
IPROACTV(:,4) = 1
IPROACTV(:,JPBUPROMAX+1) = 0
GERROR=.FALSE.
YWORK2(:,:) = ' '
YEND_COMMENT(:) = ' '
!
!                        Budget of RU
IF (LBU_RU) THEN
  IPROC=4
  IPROACTV(1,IPROC) = NASSEU
  IPROC=IPROC+1
  IF( NMODEL>1 ) IPROACTV(1,IPROC) = NNESTU
  IPROC=IPROC+1
  IPROACTV(1,IPROC) = NADVXU
  IPROC=IPROC+1
  IPROACTV(1,IPROC) = NADVYU
  IPROC=IPROC+1
  IPROACTV(1,IPROC) = NADVZU
  IPROC=IPROC+1
  IF( LFORCING ) IPROACTV(1,IPROC)  = NFRCU
  IPROC=IPROC+1
  IF( ONUDGING ) IPROACTV(1,IPROC)  = NNUDU
  IPROC=IPROC+1
  IF ( .NOT. LCARTESIAN ) THEN
    IPROACTV(1,IPROC) = NCURVU
  ELSE
    IPROACTV(1,IPROC) = 4
  END IF
  IPROC=IPROC+1
  IF ( LCORIO ) THEN
    IPROACTV(1,IPROC) = NCORU
  ELSE
    IPROACTV(1,IPROC) = 4
  END IF
  IPROC=IPROC+1
  IF ( ONUMDIFU ) IPROACTV(1,IPROC) = NDIFU
  IPROC=IPROC+1
  IF ( OHORELAX_UVWTH .OR. OVE_RELAX ) THEN
    IPROACTV(1,IPROC) = NRELU
  ELSE
    IF(OVE_RELAX .OR. OHORELAX_UVWTH .OR. OHORELAX_RV .OR.                 &
     OHORELAX_RC .OR. OHORELAX_RR .OR. OHORELAX_RI .OR. OHORELAX_RS .OR.   &
     OHORELAX_RG .OR. OHORELAX_RH .OR. OHORELAX_TKE .OR. ANY(OHORELAX_SV)) THEN
      IPROACTV(1,IPROC) = 4
    ELSE
      IPROACTV(1,IPROC) = 3
    END IF
  END IF
  IPROC=IPROC+1
  IF( ODRAGTREE ) IPROACTV(1,IPROC)  = NDRAGU
  IPROC=IPROC+1
  IF ( HTURB /= 'NONE' ) IPROACTV(1,IPROC) = NVTURBU
  IPROC=IPROC+1
  IF ( HTURB /= 'NONE' .AND. HTURBDIM == '3DIM' ) THEN
    IPROACTV(1,IPROC) = NHTURBU
  ELSE
    IF ( HTURB /= 'NONE' ) THEN
      IPROACTV(1,IPROC) = 4
    ELSE
      IPROACTV(1,IPROC) = 3
    END IF
  END IF
  IPROC=IPROC+1
  IF ( HSCONV == 'EDKF' ) IPROACTV(1,IPROC) = NMAFLU
  IPROC=IPROC+1
  IPROACTV(1,IPROC) = NPRESU
!
  YWORK2(1,1) = 'INIF_'
  YWORK2(1,2) = 'ENDF_'
  YWORK2(1,3) = 'AVEF_'
  IPROC=4
  YWORK2(1,IPROC) = 'ASSE_'
  IPROC=IPROC+1
  YWORK2(1,IPROC) = 'NEST_'
  IPROC=IPROC+1
  YWORK2(1,IPROC) = 'ADVX_'
  IPROC=IPROC+1
  YWORK2(1,IPROC) = 'ADVY_'
  IPROC=IPROC+1
  YWORK2(1,IPROC) = 'ADVZ_'
  IPROC=IPROC+1
  YWORK2(1,IPROC) = 'FRC_'
  IPROC=IPROC+1
  YWORK2(1,IPROC) = 'NUD_'
  IPROC=IPROC+1
  YWORK2(1,IPROC) = 'CURV_'
  IPROC=IPROC+1
  YWORK2(1,IPROC) = 'COR_'
  IPROC=IPROC+1
  YWORK2(1,IPROC) = 'DIF_'
  IPROC=IPROC+1
  YWORK2(1,IPROC) = 'REL_'
  IPROC=IPROC+1
  YWORK2(1,IPROC) = 'DRAG_'
  IPROC=IPROC+1
  YWORK2(1,IPROC) = 'VTURB_'
  IPROC=IPROC+1
  YWORK2(1,IPROC) = 'HTURB_'
  IPROC=IPROC+1
  YWORK2(1,IPROC) = 'MAFL_'
  IPROC=IPROC+1
  YWORK2(1,IPROC) = 'PRES_'
!
  YEND_COMMENT(1) = 'BU_RU'
  NBUPROCNBR(1) = 3
!
  CBUACTION(1,1) = 'IG'
  CBUACTION(1,2) = 'CC'
  CBUACTION(1,3) = 'ES'
!
  DO JJ=1,3
    CBUCOMMENT(1,JJ) = ADJUSTL( ADJUSTR( YWORK2(1,JJ) ) // &
                                ADJUSTL( YEND_COMMENT(1) ) )
  END DO
!
END IF
!
!                        Budget of RV
IF (LBU_RV) THEN
  IPROC=4
  IPROACTV(2,IPROC) = NASSEV
  IPROC=IPROC+1
  IF( NMODEL>1 ) IPROACTV(2,IPROC) = NNESTV
  IPROC=IPROC+1
  IPROACTV(2,IPROC) = NADVXV
  IPROC=IPROC+1
  IPROACTV(2,IPROC) = NADVYV
  IPROC=IPROC+1
  IPROACTV(2,IPROC) = NADVZV
  IPROC=IPROC+1
  IF( LFORCING ) IPROACTV(2,IPROC)  = NFRCV
  IPROC=IPROC+1
  IF( ONUDGING ) IPROACTV(2,IPROC)  = NNUDV
  IPROC=IPROC+1
  IF ( .NOT. LCARTESIAN ) THEN
    IPROACTV(2,IPROC) = NCURVV
  ELSE
    IPROACTV(2,IPROC) = 4
  END IF
  IPROC=IPROC+1
  IF ( LCORIO ) THEN
    IPROACTV(2,IPROC) = NCORV
  ELSE
    IPROACTV(2,IPROC) = 4
  END IF
  IPROC=IPROC+1
  IF ( ONUMDIFU ) IPROACTV(2,IPROC) = NDIFV
  IPROC=IPROC+1
  IF ( OHORELAX_UVWTH .OR. OVE_RELAX ) THEN
    IPROACTV(2,IPROC) = NRELV
  ELSE
    IF(OVE_RELAX .OR. OHORELAX_UVWTH .OR. OHORELAX_RV .OR.                 &
     OHORELAX_RC .OR. OHORELAX_RR .OR. OHORELAX_RI .OR. OHORELAX_RS .OR.   &
     OHORELAX_RG .OR. OHORELAX_RH .OR. OHORELAX_TKE .OR. ANY(OHORELAX_SV)) THEN
      IPROACTV(2,IPROC) = 4
    ELSE
      IPROACTV(2,IPROC) = 3
    END IF
  END IF
  IPROC=IPROC+1
  IF( ODRAGTREE ) IPROACTV(2,IPROC)  = NDRAGV
  IPROC=IPROC+1
  IF ( HTURB /= 'NONE' ) IPROACTV(2,IPROC) = NVTURBV
  IPROC=IPROC+1
  IF ( HTURB /= 'NONE' .AND. HTURBDIM == '3DIM' ) THEN
    IPROACTV(2,IPROC) = NHTURBV
  ELSE
    IF ( HTURB /= 'NONE' ) THEN
      IPROACTV(2,IPROC) = 4
    ELSE
      IPROACTV(2,IPROC) = 3
    END IF
  END IF
  IPROC=IPROC+1
  IF ( HSCONV == 'EDKF' ) IPROACTV(2,IPROC) = NMAFLV
  IPROC=IPROC+1
  IPROACTV(2,IPROC) = NPRESV
!
  YWORK2(2,1) = 'INIF_'
  YWORK2(2,2) = 'ENDF_'
  YWORK2(2,3) = 'AVEF_'
  IPROC=4
  YWORK2(2,IPROC) = 'ASSE_'
  IPROC=IPROC+1
  YWORK2(2,IPROC) = 'NEST_'
  IPROC=IPROC+1
  YWORK2(2,IPROC) = 'ADVX_'
  IPROC=IPROC+1
  YWORK2(2,IPROC) = 'ADVY_'
  IPROC=IPROC+1
  YWORK2(2,IPROC) = 'ADVZ_'
  IPROC=IPROC+1
  YWORK2(2,IPROC) = 'FRC_'
  IPROC=IPROC+1
  YWORK2(2,IPROC) = 'NUD_'
  IPROC=IPROC+1
  YWORK2(2,IPROC) = 'CURV_'
  IPROC=IPROC+1
  YWORK2(2,IPROC) = 'COR_'
  IPROC=IPROC+1
  YWORK2(2,IPROC) = 'DIF_'
  IPROC=IPROC+1
  YWORK2(2,IPROC) = 'REL_'
  IPROC=IPROC+1
  YWORK2(2,IPROC) = 'DRAG_'
  IPROC=IPROC+1
  YWORK2(2,IPROC) = 'VTURB_'
  IPROC=IPROC+1
  YWORK2(2,IPROC) = 'HTURB_'
  IPROC=IPROC+1
  YWORK2(2,IPROC) = 'MAFL_'
  IPROC=IPROC+1
  YWORK2(2,IPROC) = 'PRES_'
!
  YEND_COMMENT(2) = 'BU_RV'
  NBUPROCNBR(2) = 3
!
  CBUACTION(2,1) = 'IG'
  CBUACTION(2,2) = 'CC'
  CBUACTION(2,3) = 'ES'
!
  DO JJ=1,3
    CBUCOMMENT(2,JJ) = ADJUSTL( ADJUSTR( YWORK2(2,JJ) ) // &
                                ADJUSTL( YEND_COMMENT(2) ) )
  END DO
!
END IF
!
!                        Budget of RW
IF (LBU_RW) THEN
  IPROC=4
  IPROACTV(3,IPROC) = NASSEW
  IPROC=IPROC+1
  IF( NMODEL>1 ) IPROACTV(3,IPROC) = NNESTW
  IPROC=IPROC+1
  IPROACTV(3,IPROC) = NADVXW
  IPROC=IPROC+1
  IPROACTV(3,IPROC) = NADVYW
  IPROC=IPROC+1
  IPROACTV(3,IPROC) = NADVZW
  IPROC=IPROC+1
  IF( LFORCING ) IPROACTV(3,IPROC)  = NFRCW
  IPROC=IPROC+1
  IF( ONUDGING ) IPROACTV(3,IPROC)  = NNUDW
  IPROC=IPROC+1
  IF ( .NOT. LCARTESIAN ) THEN
    IPROACTV(3,IPROC) = NCURVW
  ELSE
    IPROACTV(3,IPROC) = 4
  END IF
  IPROC=IPROC+1
  IF ( LCORIO ) THEN
    IPROACTV(3,IPROC) = NCORW
  ELSE
    IPROACTV(3,IPROC) = 4
  END IF
  IPROC=IPROC+1
  IPROACTV(3,IPROC) = NGRAVW
  IPROC=IPROC+1
  IF ( ONUMDIFU ) IPROACTV(3,IPROC) = NDIFW
  IPROC=IPROC+1
  IF ( OHORELAX_UVWTH .OR. OVE_RELAX ) THEN
    IPROACTV(3,IPROC) = NRELW
  ELSE
    IF(OVE_RELAX .OR. OHORELAX_UVWTH .OR. OHORELAX_RV .OR.                 &
     OHORELAX_RC .OR. OHORELAX_RR .OR. OHORELAX_RI .OR. OHORELAX_RS .OR.   &
     OHORELAX_RG .OR. OHORELAX_RH .OR. OHORELAX_TKE .OR. ANY(OHORELAX_SV)) THEN
      IPROACTV(3,IPROC) = 4
    ELSE
      IPROACTV(3,IPROC) = 3
    END IF
  END IF
  IPROC=IPROC+1
  IF ( HTURB /= 'NONE' ) IPROACTV(3,IPROC) = NVTURBW
  IPROC=IPROC+1
  IF ( HTURB /= 'NONE' .AND. HTURBDIM == '3DIM' ) THEN
    IPROACTV(3,IPROC) = NHTURBW
  ELSE
    IF ( HTURB /= 'NONE' ) THEN
      IPROACTV(3,IPROC) = 4
    ELSE
      IPROACTV(3,IPROC) = 3
    END IF
  END IF
  IPROC=IPROC+1
  IPROACTV(3,IPROC) = NPRESW
!
  YWORK2(3,1) = 'INIF_'
  YWORK2(3,2) = 'ENDF_'
  YWORK2(3,3) = 'AVEF_'
  IPROC=4
  YWORK2(3,IPROC) = 'ASSE_'
  IPROC=IPROC+1
  YWORK2(3,IPROC) = 'NEST_'
  IPROC=IPROC+1
  YWORK2(3,IPROC) = 'ADVX_'
  IPROC=IPROC+1
  YWORK2(3,IPROC) = 'ADVY_'
  IPROC=IPROC+1
  YWORK2(3,IPROC) = 'ADVZ_'
  IPROC=IPROC+1
  YWORK2(3,IPROC) = 'FRC_'
  IPROC=IPROC+1
  YWORK2(3,IPROC) = 'NUD_'
  IPROC=IPROC+1
  YWORK2(3,IPROC) = 'CURV_'
  IPROC=IPROC+1
  YWORK2(3,IPROC) = 'COR_'
  IPROC=IPROC+1
  YWORK2(3,IPROC) = 'GRAV_'
  IPROC=IPROC+1
  YWORK2(3,IPROC) = 'DIF_'
  IPROC=IPROC+1
  YWORK2(3,IPROC) = 'REL_'
  IPROC=IPROC+1
  YWORK2(3,IPROC) = 'VTURB_'
  IPROC=IPROC+1
  YWORK2(3,IPROC) = 'HTURB_'
  IPROC=IPROC+1
  YWORK2(3,IPROC) = 'PRES_'
!
  YEND_COMMENT(3) = 'BU_RW'
  NBUPROCNBR(3) = 3
!
  CBUACTION(3,1) = 'IG'
  CBUACTION(3,2) = 'CC'
  CBUACTION(3,3) = 'ES'
!
  DO JJ=1,3
    CBUCOMMENT(3,JJ) = ADJUSTL( ADJUSTR( YWORK2(3,JJ) ) // &
                                ADJUSTL( YEND_COMMENT(3) ) )
  END DO
!
END IF
!
!                        Budget of RTH
IF (LBU_RTH) THEN
  IPROC=4
  IPROACTV(4,IPROC) = NASSETH
  IPROC=IPROC+1
  IF( NMODEL>1 ) IPROACTV(4,IPROC) = NNESTTH
  IPROC=IPROC+1
  IF ( HMET_ADV_SCHEME(1:3) == 'PPM' ) &
     IPROACTV(4,IPROC) = NADVTH
  IPROC=IPROC+1
  IF ( HMET_ADV_SCHEME(1:3) /= 'PPM' ) &
     IPROACTV(4,IPROC) = NADVXTH
  IPROC=IPROC+1
  IF ( HMET_ADV_SCHEME(1:3) /= 'PPM' ) &
       IPROACTV(4,IPROC) = NADVYTH
  IPROC=IPROC+1
  IF ( HMET_ADV_SCHEME(1:3) /= 'PPM' ) &
       IPROACTV(4,IPROC) = NADVZTH
  IPROC=IPROC+1
  IF( LFORCING ) IPROACTV(4,IPROC)  = NFRCTH
  IPROC=IPROC+1
  IF( ONUDGING ) IPROACTV(4,IPROC)  = NNUDTH
  IPROC=IPROC+1
  IF ( KRR > 0 ) THEN
    IF(.NOT.L1D) IPROACTV(4,IPROC) = NPREFTH
  ELSE
    IPROACTV(4,IPROC) = 4
  END IF
  IPROC=IPROC+1
  IF ( ONUMDIFTH ) IPROACTV(4,IPROC) = NDIFTH
  IPROC=IPROC+1
  IF ( OHORELAX_UVWTH .OR. OVE_RELAX ) THEN
    IPROACTV(4,IPROC) = NRELTH
  ELSE
    IF(OVE_RELAX .OR. OHORELAX_UVWTH .OR. OHORELAX_RV .OR.                 &
     OHORELAX_RC .OR. OHORELAX_RR .OR. OHORELAX_RI .OR. OHORELAX_RS .OR.   &
     OHORELAX_RG .OR. OHORELAX_RH .OR. OHORELAX_TKE .OR. ANY(OHORELAX_SV)) THEN
      IPROACTV(4,IPROC) = 4
    ELSE
      IPROACTV(4,IPROC) = 3
    END IF
  END IF
  IPROC=IPROC+1
  IF ( HRAD /= 'NONE' ) IPROACTV(4,IPROC) = NRADTH
  IPROC=IPROC+1
  IF ( HDCONV /= 'NONE' .OR. HSCONV == 'KAFR') IPROACTV(4,IPROC) = NDCONVTH
  IPROC=IPROC+1
  IF ( HTURB /= 'NONE' ) IPROACTV(4,IPROC) = NVTURBTH
  IPROC=IPROC+1
  IF ( HTURB /= 'NONE' .AND. HTURBDIM == '3DIM' ) THEN
    IPROACTV(4,IPROC) = NHTURBTH
  ELSE
    IF ( HTURB /= 'NONE' ) THEN
      IPROACTV(4,IPROC) = 4
    ELSE
      IPROACTV(4,IPROC) = 3
    END IF
  END IF
  IPROC=IPROC+1
  IF (HTURB /= 'NONE')     IPROACTV(4,IPROC) = NDISSHTH
  IPROC=IPROC+1
  IF ( HSCONV == 'EDKF' ) IPROACTV(4,IPROC) = NMAFLTH
  IPROC=IPROC+1
  IF (HCLOUD /= 'NONE')     IPROACTV(4,IPROC) = NNEGATH
  IPROC=IPROC+1
  IF (HCLOUD(1:3) == 'ICE' .OR. HCLOUD == 'C2R2' .OR. HCLOUD == 'KHKO' ) &
      IPROACTV(4,IPROC) = NHENUTH
  IPROC=IPROC+1
  IF (HCLOUD(1:3) == 'ICE') IPROACTV(4,IPROC) = NHONTH
  IPROC=IPROC+1
  IF (HCLOUD(1:3) == 'ICE') IPROACTV(4,IPROC) = NSFRTH
  IPROC=IPROC+1
  IF (HCLOUD(1:3) == 'ICE') IPROACTV(4,IPROC) = NDEPSTH
  IPROC=IPROC+1
  IF (HCLOUD(1:3) == 'ICE') IPROACTV(4,IPROC) = NDEPGTH
  IPROC=IPROC+1
  IF (((HCLOUD(1:3) == 'ICE') .AND. LWARM) .OR. ((HCLOUD == 'C2R2' &
   .OR. HCLOUD == 'KHKO') .AND. LRAIN) .OR. HCLOUD(1:3) == 'KES')             &
        IPROACTV(4,IPROC) = NREVATH
  IPROC=IPROC+1
  IF (HCLOUD(1:3) == 'ICE') IPROACTV(4,IPROC) = NRIMTH
  IPROC=IPROC+1
  IF (HCLOUD(1:3) == 'ICE') IPROACTV(4,IPROC) = NACCTH
  IPROC=IPROC+1
  IF (HCLOUD(1:3) == 'ICE') IPROACTV(4,IPROC) = NCFRZTH
  IPROC=IPROC+1
  IF (HCLOUD(1:3) == 'ICE') IPROACTV(4,IPROC) = NWETGTH
  IPROC=IPROC+1
  IF (HCLOUD(1:3) == 'ICE') IPROACTV(4,IPROC) = NDRYGTH
  IPROC=IPROC+1
  IF (HCLOUD(1:3) == 'ICE') IPROACTV(4,IPROC) = NGMLTTH
  IPROC=IPROC+1
  IF (HCLOUD == 'ICE4') IPROACTV(4,IPROC) = NWETHTH
  IPROC=IPROC+1
  IF (HCLOUD == 'ICE4') IPROACTV(4,IPROC) = NHMLTTH
  IPROC=IPROC+1
  IF (HCLOUD(1:3) == 'ICE') IPROACTV(4,IPROC) = NIMLTTH
  IPROC=IPROC+1
  IF (HCLOUD(1:3) == 'ICE') IPROACTV(4,IPROC) = NBERFITH
  IPROC=IPROC+1
  IF (HCLOUD(1:3) == 'ICE') IPROACTV(4,IPROC) = NCDEPITH
  IPROC=IPROC+1
  IF (HCLOUD == 'C2R2' .OR. HCLOUD == 'KHKO' .OR. HCLOUD(1:3) == 'KES' .OR.   &
      HCLOUD == 'REVE')   IPROACTV(4,IPROC) = NCONDTH
!
  YWORK2(4,1) = 'INIF_'
  YWORK2(4,2) = 'ENDF_'
  YWORK2(4,3) = 'AVEF_'
  IPROC=4
  YWORK2(4,IPROC) = 'ASSE_'
  IPROC=IPROC+1
  YWORK2(4,IPROC) = 'NEST_'
  IPROC=IPROC+1
  YWORK2(4,IPROC) = 'ADV_'
  IPROC=IPROC+1
  YWORK2(4,IPROC) = 'ADVX_'
  IPROC=IPROC+1
  YWORK2(4,IPROC) = 'ADVY_'
  IPROC=IPROC+1
  YWORK2(4,IPROC) = 'ADVZ_'
  IPROC=IPROC+1
  YWORK2(4,IPROC) = 'FRC_'
  IPROC=IPROC+1
  YWORK2(4,IPROC) = 'NUD_'
  IPROC=IPROC+1
  YWORK2(4,IPROC) = 'PREF_'
  IPROC=IPROC+1
  YWORK2(4,IPROC) = 'DIF_'
  IPROC=IPROC+1
  YWORK2(4,IPROC) = 'REL_'
  IPROC=IPROC+1
  YWORK2(4,IPROC) = 'RAD_'
  IPROC=IPROC+1
  YWORK2(4,IPROC) = 'DCONV_'
  IPROC=IPROC+1
  YWORK2(4,IPROC) = 'VTURB_'
  IPROC=IPROC+1
  YWORK2(4,IPROC) = 'HTURB_'
  IPROC=IPROC+1
  YWORK2(4,IPROC) = 'DISSH_'
  IPROC=IPROC+1
  YWORK2(4,IPROC) = 'MAFL_'
  IPROC=IPROC+1
  YWORK2(4,IPROC) = 'NEGA_'
  IPROC=IPROC+1
  YWORK2(4,IPROC) = 'HENU_'
  IPROC=IPROC+1
  YWORK2(4,IPROC) = 'HON_'
  IPROC=IPROC+1
  YWORK2(4,IPROC) = 'SFR_'
  IPROC=IPROC+1
  YWORK2(4,IPROC) = 'DEPS_'
  IPROC=IPROC+1
  YWORK2(4,IPROC) = 'DEPG_'
  IPROC=IPROC+1
  YWORK2(4,IPROC) = 'REVA_'
  IPROC=IPROC+1
  YWORK2(4,IPROC) = 'RIM_'
  IPROC=IPROC+1
  YWORK2(4,IPROC) = 'ACC_'
  IPROC=IPROC+1
  YWORK2(4,IPROC) = 'CFRZ_'
  IPROC=IPROC+1
  YWORK2(4,IPROC) = 'WETG_'
  IPROC=IPROC+1
  YWORK2(4,IPROC) = 'DRYG_'
  IPROC=IPROC+1
  YWORK2(4,IPROC) = 'GMLT_'
  IPROC=IPROC+1
  YWORK2(4,IPROC) = 'WETH_'
  IPROC=IPROC+1
  YWORK2(4,IPROC) = 'HMLT_'
  IPROC=IPROC+1
  YWORK2(4,IPROC) = 'IMLT_'
  IPROC=IPROC+1
  YWORK2(4,IPROC) = 'BERFI_'
  IPROC=IPROC+1
  YWORK2(4,IPROC) = 'CDEPI_'
  IPROC=IPROC+1
  YWORK2(4,IPROC) = 'COND_'
!
  YEND_COMMENT(4) = 'BU_RTH'
  NBUPROCNBR(4) = 3
!
  CBUACTION(4,1) = 'IG'
  CBUACTION(4,2) = 'CC'
  CBUACTION(4,3) = 'ES'
!
  DO JJ=1,3
    CBUCOMMENT(4,JJ) = ADJUSTL( ADJUSTR( YWORK2(4,JJ) ) // &
                                ADJUSTL( YEND_COMMENT(4) ) )
  END DO
!
END IF
!
!                        Budget of RTKE
IF (LBU_RTKE) THEN
  IPROC=4
  IPROACTV(5,IPROC) = NASSETKE
  IPROC=IPROC+1
  IF ( HMET_ADV_SCHEME(1:3) == 'PPM' ) &
     IPROACTV(5,IPROC) = NADVTKE
  IPROC=IPROC+1
  IF ( HMET_ADV_SCHEME(1:3) /= 'PPM' ) &
     IPROACTV(5,IPROC) = NADVXTKE
  IPROC=IPROC+1
  IF ( HMET_ADV_SCHEME(1:3) /= 'PPM' ) &
     IPROACTV(5,IPROC) = NADVYTKE
  IPROC=IPROC+1
  IF ( HMET_ADV_SCHEME(1:3) /= 'PPM' ) &
     IPROACTV(5,IPROC) = NADVZTKE
  IPROC=IPROC+1
  IF( LFORCING ) IPROACTV(5,IPROC)  = NFRCTKE
  IPROC=IPROC+1
  IF ( ONUMDIFTH ) IPROACTV(5,IPROC) = NDIFTKE
  IPROC=IPROC+1
  IF ( OHORELAX_TKE ) THEN
    IPROACTV(5,IPROC) = NRELTKE
  ELSE
    IF(OVE_RELAX .OR. OHORELAX_UVWTH .OR. OHORELAX_RV .OR.                 &
     OHORELAX_RC .OR. OHORELAX_RR .OR. OHORELAX_RI .OR. OHORELAX_RS .OR.   &
     OHORELAX_RG .OR. OHORELAX_RH .OR. OHORELAX_TKE .OR. ANY(OHORELAX_SV)) THEN
      IPROACTV(5,IPROC) = 4
    ELSE
      IPROACTV(5,IPROC) = 3
    END IF
  END IF
  IPROC=IPROC+1
  IF( ODRAGTREE ) IPROACTV(5,IPROC)  = NDRAGTKE
  IPROC=IPROC+1
  IPROACTV(5,IPROC) = NDPTKE
  IPROC=IPROC+1
  IPROACTV(5,IPROC) = NTPTKE
  IPROC=IPROC+1
  IPROACTV(5,IPROC) = NDISSTKE
  IPROC=IPROC+1
  IPROACTV(5,IPROC) = NTRTKE
!
  YWORK2(5,1) = 'INIF_'
  YWORK2(5,2) = 'ENDF_'
  YWORK2(5,3) = 'AVEF_'
  IPROC=4
  YWORK2(5,IPROC) = 'ASSE_'
  IPROC=IPROC+1
  YWORK2(5,IPROC) = 'ADV_'
  IPROC=IPROC+1
  YWORK2(5,IPROC) = 'ADVX_'
  IPROC=IPROC+1
  YWORK2(5,IPROC) = 'ADVY_'
  IPROC=IPROC+1
  YWORK2(5,IPROC) = 'ADVZ_'
  IPROC=IPROC+1
  YWORK2(5,IPROC) = 'FRC_'
  IPROC=IPROC+1
  YWORK2(5,IPROC) = 'DIF_'
  IPROC=IPROC+1
  YWORK2(5,IPROC) = 'REL_'
  IPROC=IPROC+1
  YWORK2(5,IPROC) = 'DRAG_'
  IPROC=IPROC+1
  YWORK2(5,IPROC) = 'DP_'
  IPROC=IPROC+1
  YWORK2(5,IPROC) = 'TP_'
  IPROC=IPROC+1
  YWORK2(5,IPROC) = 'DISS_'
  IPROC=IPROC+1
  YWORK2(5,IPROC) = 'TR_'
!
  YEND_COMMENT(5) = 'BU_RTKE'
  NBUPROCNBR(5) = 3
!
  CBUACTION(5,1) = 'IG'
  CBUACTION(5,2) = 'CC'
  CBUACTION(5,3) = 'ES'
!
  DO JJ=1,3
    CBUCOMMENT(5,JJ) = ADJUSTL( ADJUSTR( YWORK2(5,JJ) ) // &
                                ADJUSTL( YEND_COMMENT(5) ) )
  END DO
!
END IF
!
!                        Budget of RRV
IF (LBU_RRV) THEN
  IPROC=4
  IPROACTV(6,IPROC) = NASSERV
  IPROC=IPROC+1
  IF( NMODEL>1 ) IPROACTV(6,IPROC) = NNESTRV
  IPROC=IPROC+1
  IF ( HMET_ADV_SCHEME(1:3) == 'PPM' ) &
     IPROACTV(6,IPROC) = NADVRV
  IPROC=IPROC+1
  IF ( HMET_ADV_SCHEME(1:3) /= 'PPM' ) &
     IPROACTV(6,IPROC) = NADVXRV
  IPROC=IPROC+1
  IF ( HMET_ADV_SCHEME(1:3) /= 'PPM' ) &
     IPROACTV(6,IPROC) = NADVYRV
  IPROC=IPROC+1
  IF ( HMET_ADV_SCHEME(1:3) /= 'PPM' ) &
     IPROACTV(6,IPROC) = NADVZRV
  IPROC=IPROC+1
  IF( LFORCING ) IPROACTV(6,IPROC)  = NFRCRV
  IPROC=IPROC+1
  IF( ONUDGING ) IPROACTV(6,IPROC)  = NNUDRV
  IPROC=IPROC+1
  IF ( ONUMDIFTH ) IPROACTV(6,IPROC) = NDIFRV
  IPROC=IPROC+1
  IF ( OHORELAX_RV .OR. OVE_RELAX ) THEN
    IPROACTV(6,IPROC) = NRELRV
  ELSE
    IF(OVE_RELAX .OR. OHORELAX_UVWTH .OR. OHORELAX_RV .OR.                 &
     OHORELAX_RC .OR. OHORELAX_RR .OR. OHORELAX_RI .OR. OHORELAX_RS .OR.   &
     OHORELAX_RG .OR. OHORELAX_RH .OR. OHORELAX_TKE .OR. ANY(OHORELAX_SV)) THEN
      IPROACTV(6,IPROC) = 4
    ELSE
      IPROACTV(6,IPROC) = 3
    END IF
  END IF
  IPROC=IPROC+1
  IF ( HDCONV /= 'NONE' .OR. HSCONV == 'KAFR') IPROACTV(6,IPROC) = NDCONVRV
  IPROC=IPROC+1
  IF ( HTURB /= 'NONE' ) IPROACTV(6,IPROC) = NVTURBRV
  IPROC=IPROC+1
  IF ( HTURB /= 'NONE' .AND. HTURBDIM == '3DIM' ) THEN
    IPROACTV(6,IPROC) = NHTURBRV
  ELSE
    IF ( HTURB /= 'NONE' ) THEN
      IPROACTV(6,IPROC) = 4
    ELSE
      IPROACTV(6,IPROC) = 3
    END IF
  END IF
  IPROC=IPROC+1
  IF ( HSCONV == 'EDKF' ) IPROACTV(6,IPROC) = NMAFLRV
  IPROC=IPROC+1
  IF ( HCLOUD /= 'NONE' ) IPROACTV(6,IPROC) = NNEGARV
  IPROC=IPROC+1
  IF (HCLOUD == 'C2R2' .OR. HCLOUD == 'KHKO' .OR. HCLOUD(1:3) == 'ICE')  &
        IPROACTV(6,IPROC) = NHENURV
  IPROC=IPROC+1
  IF (HCLOUD(1:3) == 'ICE') IPROACTV(6,IPROC) = NDEPSRV
  IPROC=IPROC+1
  IF (HCLOUD(1:3) == 'ICE') IPROACTV(6,IPROC) = NDEPGRV
  IPROC=IPROC+1
  IF (HCLOUD(1:3) == 'KES' .OR. ((HCLOUD == 'C2R2' .OR. HCLOUD == 'KHKO') .AND. LRAIN) .OR. &
      ((HCLOUD(1:3) == 'ICE') .AND. LWARM)) IPROACTV(6,IPROC) = NREVARV
  IPROC=IPROC+1
  IF (HCLOUD(1:3) == 'KES'.OR. HCLOUD == 'REVE' .OR. &
      HCLOUD == 'C2R2' .OR. HCLOUD == 'KHKO' ) IPROACTV(6,IPROC) = NCONDRV
  IPROC=IPROC+1
  IF (HCLOUD(1:3) == 'ICE') IPROACTV(6,IPROC) = NCDEPIRV
  IPROC=IPROC+1
!
  YWORK2(6,1) = 'INIF_'
  YWORK2(6,2) = 'ENDF_'
  YWORK2(6,3) = 'AVEF_'
  IPROC=4
  YWORK2(6,IPROC) = 'ASSE_'
  IPROC=IPROC+1
  YWORK2(6,IPROC) = 'NEST_'
  IPROC=IPROC+1
  YWORK2(6,IPROC) = 'ADV_'
  IPROC=IPROC+1
  YWORK2(6,IPROC) = 'ADVX_'
  IPROC=IPROC+1
  YWORK2(6,IPROC) = 'ADVY_'
  IPROC=IPROC+1
  YWORK2(6,IPROC) = 'ADVZ_'
  IPROC=IPROC+1
  YWORK2(6,IPROC) = 'FRC_'
  IPROC=IPROC+1
  YWORK2(6,IPROC) = 'NUD_'
  IPROC=IPROC+1
  YWORK2(6,IPROC) = 'DIF_'
  IPROC=IPROC+1
  YWORK2(6,IPROC) = 'REL_'
  IPROC=IPROC+1
  YWORK2(6,IPROC) = 'DCONV_'
  IPROC=IPROC+1
  YWORK2(6,IPROC) = 'VTURB_'
  IPROC=IPROC+1
  YWORK2(6,IPROC) = 'HTURB_'
  IPROC=IPROC+1
  YWORK2(6,IPROC) = 'MAFL_'
  IPROC=IPROC+1
  YWORK2(6,IPROC) = 'NEGA_'
  IPROC=IPROC+1
  YWORK2(6,IPROC) = 'HENU_'
  IPROC=IPROC+1
  YWORK2(6,IPROC) = 'DEPS_'
  IPROC=IPROC+1
  YWORK2(6,IPROC) = 'DEPG_'
  IPROC=IPROC+1
  YWORK2(6,IPROC) = 'REVA_'
  IPROC=IPROC+1
  YWORK2(6,IPROC) = 'COND_'
  IPROC=IPROC+1
  YWORK2(6,IPROC) = 'CDEPI_'
!
  YEND_COMMENT(6) = 'BU_RRV'
  NBUPROCNBR(6) = 3
!
  CBUACTION(6,1) = 'IG'
  CBUACTION(6,2) = 'CC'
  CBUACTION(6,3) = 'ES'
!
  DO JJ=1,3
    CBUCOMMENT(6,JJ) = ADJUSTL( ADJUSTR( YWORK2(6,JJ) ) // &
                                ADJUSTL( YEND_COMMENT(6) ) )
  END DO
!
END IF
!
!                        Budget of RRC
IF (LBU_RRC) THEN
  IPROC=4
  IPROACTV(7,IPROC) = NASSERC
  IPROC=IPROC+1
  IF( NMODEL>1 ) IPROACTV(7,IPROC) = NNESTRC
  IPROC=IPROC+1
  IF ( HMET_ADV_SCHEME(1:3) == 'PPM' ) &
     IPROACTV(7,IPROC) = NADVRC
  IPROC=IPROC+1
  IF ( HMET_ADV_SCHEME(1:3) /= 'PPM' ) &
     IPROACTV(7,IPROC) = NADVXRC
  IPROC=IPROC+1
  IF ( HMET_ADV_SCHEME(1:3) /= 'PPM' ) &
     IPROACTV(7,IPROC) = NADVYRC
  IPROC=IPROC+1
  IF ( HMET_ADV_SCHEME(1:3) /= 'PPM' ) &
     IPROACTV(7,IPROC) = NADVZRC
  IPROC=IPROC+1
  IF( LFORCING ) IPROACTV(7,IPROC)  = NFRCRC
  IPROC=IPROC+1
  IF ( ONUMDIFTH ) IPROACTV(7,IPROC) = NDIFRC
  IPROC=IPROC+1
  IF ( OHORELAX_RC ) THEN
    IPROACTV(7,IPROC) = NRELRC
  ELSE
    IF(OVE_RELAX .OR. OHORELAX_UVWTH .OR. OHORELAX_RV .OR.                 &
     OHORELAX_RC .OR. OHORELAX_RR .OR. OHORELAX_RI .OR. OHORELAX_RS .OR.   &
     OHORELAX_RG .OR. OHORELAX_RH .OR. OHORELAX_TKE .OR. ANY(OHORELAX_SV)) THEN
      IPROACTV(7,IPROC) = 4
    ELSE
      IPROACTV(7,IPROC) = 3
    END IF
  END IF
  IPROC=IPROC+1
  IF( HDCONV /= 'NONE' .OR. HSCONV == 'KAFR') IPROACTV(7,IPROC)  = NDCONVRC
  IPROC=IPROC+1
  IF ( HTURB /= 'NONE' ) IPROACTV(7,IPROC) = NVTURBRC
  IPROC=IPROC+1
  IF ( HTURB /= 'NONE' .AND. HTURBDIM == '3DIM' ) THEN
    IPROACTV(7,IPROC) = NHTURBRC
  ELSE
    IF ( HTURB /= 'NONE' ) THEN
      IPROACTV(7,IPROC) = 4
    ELSE
      IPROACTV(7,IPROC) = 3
    END IF
  END IF
  IPROC=IPROC+1
  IF( HCLOUD /= 'NONE' ) IPROACTV(7,IPROC)  = NNEGARC
  IPROC=IPROC+1
  IF (HCLOUD(1:3) == 'KES' )                      IPROACTV(7,IPROC  ) = NACCRRC
  IPROC=IPROC+1
  IF (HCLOUD(1:3) == 'KES' )                      IPROACTV(7,IPROC) = NAUTORC
  IPROC=IPROC+1
  IF (HCLOUD == 'C2R2' .OR. HCLOUD == 'KHKO') IPROACTV(7,IPROC) = NHENURC
  IPROC=IPROC+1
!
  IF (HCLOUD(1:3) == 'ICE') IPROACTV(7,IPROC) = NHONRC
  IPROC=IPROC+1
  IF (((HCLOUD(1:3) == 'ICE' ) .AND. LWARM) .OR.  ((HCLOUD == 'C2R2' .OR. &
        HCLOUD == 'KHKO') .AND. LRAIN)) IPROACTV(7,IPROC) = NAUTORC
  IPROC=IPROC+1
  !modif
  IF (((HCLOUD(1:3) == 'ICE' ) .AND. LWARM) .OR.  ((HCLOUD == 'C2R2' .OR. &
        HCLOUD == 'KHKO') .AND. LRAIN)) IPROACTV(7,IPROC) = NACCRRC
  IPROC=IPROC+1
  IF (HCLOUD(1:3) == 'ICE') IPROACTV(7,IPROC) = NRIMRC
  IPROC=IPROC+1
  IF (HCLOUD(1:3) == 'ICE') IPROACTV(7,IPROC) = NWETGRC
  IPROC=IPROC+1
  IF (HCLOUD(1:3) == 'ICE') IPROACTV(7,IPROC) = NDRYGRC
  IPROC=IPROC+1
  IF (HCLOUD == 'ICE4')  IPROACTV(7,IPROC) = NWETHRC
  IPROC=IPROC+1
  IF (HCLOUD(1:3) == 'ICE' ) IPROACTV(7,IPROC) = NIMLTRC
  IPROC=IPROC+1
  IF (HCLOUD(1:3) == 'ICE') IPROACTV(7,IPROC) = NBERFIRC
  IPROC=IPROC+1
  IF (((HCLOUD == 'C2R2' .OR. HCLOUD == 'KHKO').AND. LSEDC) .OR.   &
     (HCLOUD(1:3) == 'ICE' .AND. LSEDIC)) IPROACTV(7,IPROC) = NSEDIRC
  IPROC=IPROC+1
  IF (HCLOUD(1:3) == 'ICE') IPROACTV(7,IPROC) = NCDEPIRC
  IPROC=IPROC+1
  IF (HCLOUD == 'C2R2'.OR. HCLOUD == 'KHKO' .OR. &
      HCLOUD(1:3) == 'KES'.OR. HCLOUD == 'REVE') IPROACTV(7,IPROC) = NCONDRC
  IPROC=IPROC+1
!
  YWORK2(7,1) = 'INIF_'
  YWORK2(7,2) = 'ENDF_'
  YWORK2(7,3) = 'AVEF_'
  IPROC=4
  YWORK2(7,IPROC) = 'ASSE_'
  IPROC=IPROC+1
  YWORK2(7,IPROC) = 'NEST_'
  IPROC=IPROC+1
  YWORK2(7,IPROC) = 'ADV_'
  IPROC=IPROC+1
  YWORK2(7,IPROC) = 'ADVX_'
  IPROC=IPROC+1
  YWORK2(7,IPROC) = 'ADVY_'
  IPROC=IPROC+1
  YWORK2(7,IPROC) = 'ADVZ_'
  IPROC=IPROC+1
  YWORK2(7,IPROC) = 'FRC_'
  IPROC=IPROC+1
  YWORK2(7,IPROC) = 'DIF_'
  IPROC=IPROC+1
  YWORK2(7,IPROC) = 'REL_'
  IPROC=IPROC+1
  YWORK2(7,IPROC) = 'DCONV_'
  IPROC=IPROC+1
  YWORK2(7,IPROC) = 'VTURB_'
  IPROC=IPROC+1
  YWORK2(7,IPROC) = 'HTURB_'
  IPROC=IPROC+1
  YWORK2(7,IPROC) = 'NEGA_'
  IPROC=IPROC+1
  YWORK2(7,IPROC) = 'ACCR_'
  IPROC=IPROC+1
  YWORK2(7,IPROC) = 'AUTO_'
  IPROC=IPROC+1
  YWORK2(7,IPROC) = 'HENU_'
  IPROC=IPROC+1
  YWORK2(7,IPROC) = 'HON_'
  IPROC=IPROC+1
  YWORK2(7,IPROC) = 'AUTO_'
  IPROC=IPROC+1
  YWORK2(7,IPROC) = 'ACCR_'
  IPROC=IPROC+1
  YWORK2(7,IPROC) = 'RIM_'
  IPROC=IPROC+1
  YWORK2(7,IPROC) = 'WETG_'
  IPROC=IPROC+1
  YWORK2(7,IPROC) = 'DRYG_'
  IPROC=IPROC+1
  YWORK2(7,IPROC) = 'WETH_'
  IPROC=IPROC+1
  YWORK2(7,IPROC) = 'IMLT_'
  IPROC=IPROC+1
  YWORK2(7,IPROC) = 'BERFI_'
  IPROC=IPROC+1
  YWORK2(7,IPROC) = 'SEDI_'
  IPROC=IPROC+1
  YWORK2(7,IPROC) = 'CDEPI_'
  IPROC=IPROC+1
  YWORK2(7,IPROC) = 'COND_'
!
  YEND_COMMENT(7) = 'BU_RRC'
  NBUPROCNBR(7) = 3
!
  CBUACTION(7,1) = 'IG'
  CBUACTION(7,2) = 'CC'
  CBUACTION(7,3) = 'ES'
!
  DO JJ=1,3
    CBUCOMMENT(7,JJ) = ADJUSTL( ADJUSTR( YWORK2(7,JJ) ) // &
                                ADJUSTL( YEND_COMMENT(7) ) )
  END DO
!
END IF
!
!                        Budget of RRR
IF (LBU_RRR) THEN
  IPROC=4
  IPROACTV(8,IPROC) = NASSERR
  IPROC=IPROC+1
  IF( NMODEL>1 ) IPROACTV(8,IPROC) = NNESTRR
  IPROC=IPROC+1
  IF ( HMET_ADV_SCHEME(1:3) == 'PPM' ) &
     IPROACTV(8,IPROC) = NADVRR
  IPROC=IPROC+1
  IF ( HMET_ADV_SCHEME(1:3) /= 'PPM' ) &
     IPROACTV(8,IPROC) = NADVXRR
  IPROC=IPROC+1
  IF ( HMET_ADV_SCHEME(1:3) /= 'PPM' ) &
     IPROACTV(8,IPROC) = NADVYRR
  IPROC=IPROC+1
  IF ( HMET_ADV_SCHEME(1:3) /= 'PPM' ) &
     IPROACTV(8,IPROC) = NADVZRR
  IPROC=IPROC+1
  IF( LFORCING ) IPROACTV(8,IPROC)  = NFRCRR
  IPROC=IPROC+1
  IF ( ONUMDIFTH ) IPROACTV(8,IPROC) = NDIFRR
  IPROC=IPROC+1
  IF ( OHORELAX_RR ) THEN
    IPROACTV(8,IPROC) = NRELRR
  ELSE
    IF(OVE_RELAX .OR. OHORELAX_UVWTH .OR. OHORELAX_RV .OR.                 &
     OHORELAX_RC .OR. OHORELAX_RR .OR. OHORELAX_RI .OR. OHORELAX_RS .OR.   &
     OHORELAX_RG .OR. OHORELAX_RH .OR. OHORELAX_TKE .OR. ANY(OHORELAX_SV)) THEN
      IPROACTV(8,IPROC) = 4
    ELSE
      IPROACTV(8,IPROC) = 3
    END IF
  END IF
  IPROC=IPROC+1
  IF ( HCLOUD /= 'NONE' ) IPROACTV(8,IPROC) = NNEGARR
  IPROC=IPROC+1
  IF (HCLOUD(1:3) == 'KES' )   IPROACTV(8,IPROC) = NSEDIRR
  IPROC=IPROC+1
  IF (HCLOUD(1:3) == 'KES' )   IPROACTV(8,IPROC) = NACCRRR
  IPROC=IPROC+1
  IF (HCLOUD(1:3) == 'KES' )   IPROACTV(8,IPROC) = NAUTORR
  IPROC=IPROC+1
  IF (HCLOUD(1:3) == 'KES' )   IPROACTV(8,IPROC) = NREVARR
  IPROC=IPROC+1
!
  IF (HCLOUD(1:3) == 'ICE') IPROACTV(8,IPROC) = NSFRRR
  IPROC=IPROC+1
  IF ((HCLOUD(1:3) == 'ICE' ) .AND. LWARM) &
                                               IPROACTV(8,IPROC) = NAUTORR
  IPROC=IPROC+1
  IF ((HCLOUD(1:3) == 'ICE' ) .AND. LWARM) &
                                               IPROACTV(8,IPROC) = NACCRRR
  IPROC=IPROC+1
  IF ((HCLOUD(1:3) == 'ICE' ) .AND. LWARM) &
                                               IPROACTV(8,IPROC) = NREVARR
  IPROC=IPROC+1
  IF (HCLOUD(1:3) == 'ICE') IPROACTV(8,IPROC) = NACCRR
  IPROC=IPROC+1
  IF (HCLOUD(1:3) == 'ICE') IPROACTV(8,IPROC) = NCFRZRR
  IPROC=IPROC+1
  IF (HCLOUD(1:3) == 'ICE') IPROACTV(8,IPROC) = NWETGRR
  IPROC=IPROC+1
  IF (HCLOUD(1:3) == 'ICE') IPROACTV(8,IPROC) = NDRYGRR
  IPROC=IPROC+1
  IF (HCLOUD(1:3) == 'ICE') IPROACTV(8,IPROC) = NGMLTRR
  IPROC=IPROC+1
  IF (HCLOUD == 'ICE4') IPROACTV(8,IPROC) = NWETHRR
  IPROC=IPROC+1
  IF (HCLOUD == 'ICE4') IPROACTV(8,IPROC) = NHMLTRR
  IPROC=IPROC+1
  IF ((HCLOUD == 'C2R2' .OR. HCLOUD == 'KHKO') .AND. LRAIN)  IPROACTV(8,IPROC) = NAUTORR
  IPROC=IPROC+1
  IF ((HCLOUD == 'C2R2' .OR. HCLOUD == 'KHKO') .AND. LRAIN)  IPROACTV(8,IPROC) = NACCRRR
  IPROC=IPROC+1
  IF ((HCLOUD == 'C2R2'.OR. HCLOUD == 'KHKO') .AND. LRAIN)   IPROACTV(8,IPROC) = NREVARR
  IPROC=IPROC+1
  IF (HCLOUD == 'C2R2' .OR. HCLOUD == 'KHKO' .OR. HCLOUD(1:3) == 'ICE')  &
   IPROACTV(8,IPROC) = NSEDIRR
  IPROC=IPROC+1
!
  YWORK2(8,1) = 'INIF_'
  YWORK2(8,2) = 'ENDF_'
  YWORK2(8,3) = 'AVEF_'
  IPROC=4
  YWORK2(8,IPROC) = 'ASSE_'
  IPROC=IPROC+1
  YWORK2(8,IPROC) = 'NEST_'
  IPROC=IPROC+1
  YWORK2(8,IPROC) = 'ADV_'
  IPROC=IPROC+1
  YWORK2(8,IPROC) = 'ADVX_'
  IPROC=IPROC+1
  YWORK2(8,IPROC) = 'ADVY_'
  IPROC=IPROC+1
  YWORK2(8,IPROC) = 'ADVZ_'
  IPROC=IPROC+1
  YWORK2(8,IPROC) = 'FRC_'
  IPROC=IPROC+1
  YWORK2(8,IPROC) = 'DIF_'
  IPROC=IPROC+1
  YWORK2(8,IPROC) = 'REL_'
  IPROC=IPROC+1
  YWORK2(8,IPROC) = 'NEGA_'
  IPROC=IPROC+1
  YWORK2(8,IPROC) = 'SEDI_'
  IPROC=IPROC+1
  YWORK2(8,IPROC) = 'ACCR_'
  IPROC=IPROC+1
  YWORK2(8,IPROC) = 'AUTO_'
  IPROC=IPROC+1
  YWORK2(8,IPROC) = 'REVA_'
  IPROC=IPROC+1
  YWORK2(8,IPROC) = 'SFR_'
  IPROC=IPROC+1
  YWORK2(8,IPROC) = 'AUTO_'
  IPROC=IPROC+1
  YWORK2(8,IPROC) = 'ACCR_'
  IPROC=IPROC+1
  YWORK2(8,IPROC) = 'REVA_'
  IPROC=IPROC+1
  YWORK2(8,IPROC) = 'ACC_'
  IPROC=IPROC+1
  YWORK2(8,IPROC) = 'CFRZ_'
  IPROC=IPROC+1
  YWORK2(8,IPROC) = 'WETG_'
  IPROC=IPROC+1
  YWORK2(8,IPROC) = 'DRYG_'
  IPROC=IPROC+1
  YWORK2(8,IPROC) = 'GMLT_'
  IPROC=IPROC+1
  YWORK2(8,IPROC) = 'WETH_'
  IPROC=IPROC+1
  YWORK2(8,IPROC) = 'HMLT_'
  IPROC=IPROC+1
  YWORK2(8,IPROC) = 'AUTO_'
  IPROC=IPROC+1
  YWORK2(8,IPROC) = 'ACCR_'
  IPROC=IPROC+1
  YWORK2(8,IPROC) = 'REVA_'
  IPROC=IPROC+1
  YWORK2(8,IPROC) = 'SEDI_'
!
  YEND_COMMENT(8) = 'BU_RRR'
  NBUPROCNBR(8) = 3
!
  CBUACTION(8,1) = 'IG'
  CBUACTION(8,2) = 'CC'
  CBUACTION(8,3) = 'ES'
!
  DO JJ=1,3
    CBUCOMMENT(8,JJ) = ADJUSTL( ADJUSTR( YWORK2(8,JJ) ) // &
                                ADJUSTL( YEND_COMMENT(8) ) )
  END DO
!
END IF
!
!                        Budget of RRI
IF (LBU_RRI) THEN
  IPROC=4
  IPROACTV(9,IPROC) = NASSERI
  IPROC=IPROC+1
  IF( NMODEL>1 ) IPROACTV(9,IPROC) = NNESTRI
  IPROC=IPROC+1
  IF ( HMET_ADV_SCHEME(1:3) == 'PPM' ) &
     IPROACTV(9,IPROC) = NADVRI
  IPROC=IPROC+1
  IF ( HMET_ADV_SCHEME(1:3) /= 'PPM' ) &
     IPROACTV(9,IPROC) = NADVXRI
  IPROC=IPROC+1
  IF ( HMET_ADV_SCHEME(1:3) /= 'PPM' ) &
     IPROACTV(9,IPROC) = NADVYRI
  IPROC=IPROC+1
  IF ( HMET_ADV_SCHEME(1:3) /= 'PPM' ) &
     IPROACTV(9,IPROC) = NADVZRI
  IPROC=IPROC+1
  IF( LFORCING ) IPROACTV(9,IPROC)  = NFRCRI
  IPROC=IPROC+1
  IF( ONUMDIFTH ) IPROACTV(9,IPROC) = NDIFRI
  IPROC=IPROC+1
  IF ( OHORELAX_RI ) THEN
    IPROACTV(9,IPROC) = NRELRI
  ELSE
    IF(OVE_RELAX .OR. OHORELAX_UVWTH .OR. OHORELAX_RV .OR.                 &
     OHORELAX_RC .OR. OHORELAX_RR .OR. OHORELAX_RI .OR. OHORELAX_RS .OR.   &
     OHORELAX_RG .OR. OHORELAX_RH .OR. OHORELAX_TKE .OR. ANY(OHORELAX_SV)) THEN
      IPROACTV(9,IPROC) = 4
    ELSE
      IPROACTV(9,IPROC) = 3
    END IF
  END IF
  IPROC=IPROC+1
  IF( HDCONV /= 'NONE' .OR. HSCONV == 'KAFR') IPROACTV(9,IPROC) = NDCONVRI
  IPROC=IPROC+1
  IF ( HTURB /= 'NONE' ) IPROACTV(9,IPROC) = NVTURBRI
  IPROC=IPROC+1
  IF ( HTURB /= 'NONE' .AND. HTURBDIM == '3DIM' ) THEN
    IPROACTV(9,IPROC) = NHTURBRI
  ELSE
    IF ( HTURB /= 'NONE' ) THEN
      IPROACTV(9,IPROC) = 4
    ELSE
      IPROACTV(9,IPROC) = 3
    END IF
  END IF
  IPROC=IPROC+1
  IF( HCLOUD /= 'NONE' ) IPROACTV(9,IPROC) = NNEGARI
  IPROC=IPROC+1
  IF( HCLOUD(1:3) == 'ICE') IPROACTV(9,IPROC) = NHENURI
  IPROC=IPROC+1
  IF( HCLOUD(1:3) == 'ICE') IPROACTV(9,IPROC) = NHONRI
  IPROC=IPROC+1
  IF( HCLOUD(1:3) == 'ICE') IPROACTV(9,IPROC) = NAGGSRI
  IPROC=IPROC+1
  IF( HCLOUD(1:3) == 'ICE') IPROACTV(9,IPROC) = NAUTSRI
  IPROC=IPROC+1
  IF( HCLOUD(1:3) == 'ICE') IPROACTV(9,IPROC) = NCFRZRI
  IPROC=IPROC+1
  IF( HCLOUD(1:3) == 'ICE') IPROACTV(9,IPROC) = NWETGRI
  IPROC=IPROC+1
  IF( HCLOUD(1:3) == 'ICE') IPROACTV(9,IPROC) = NDRYGRI
  IPROC=IPROC+1
  IF( HCLOUD == 'ICE4') IPROACTV(9,IPROC) = NWETHRI
  IPROC=IPROC+1
  IF( HCLOUD(1:3) == 'ICE') IPROACTV(9,IPROC) = NIMLTRI
  IPROC=IPROC+1
  IF( HCLOUD(1:3) == 'ICE') IPROACTV(9,IPROC) = NBERFIRI
  IPROC=IPROC+1
  IF( HCLOUD(1:3) == 'ICE') IPROACTV(9,IPROC) = NSEDIRI
  IPROC=IPROC+1
  IF( HCLOUD(1:3) == 'ICE') IPROACTV(9,IPROC) = NCDEPIRI
  IPROC=IPROC+1
!
  YWORK2(9,1) = 'INIF_'
  YWORK2(9,2) = 'ENDF_'
  YWORK2(9,3) = 'AVEF_'
  IPROC=4
  YWORK2(9,IPROC) = 'ASSE_'
  IPROC=  IPROC+1
  YWORK2(9,IPROC) = 'NEST_'
  IPROC=  IPROC+1
  YWORK2(9,IPROC) = 'ADV_'
  IPROC=  IPROC+1
  YWORK2(9,IPROC) = 'ADVX_'
  IPROC=  IPROC+1
  YWORK2(9,IPROC) = 'ADVY_'
  IPROC=  IPROC+1
  YWORK2(9,IPROC) = 'ADVZ_'
  IPROC=  IPROC+1
  YWORK2(9,IPROC) = 'FRC_'
  IPROC=  IPROC+1
  YWORK2(9,IPROC) = 'DIF_'
  IPROC=  IPROC+1
  YWORK2(9,IPROC) = 'REL_'
  IPROC=IPROC+1
  YWORK2(9,IPROC) = 'DCONV_'
  IPROC=IPROC+1
  YWORK2(9,IPROC) = 'VTURB_'
  IPROC=IPROC+1
  YWORK2(9,IPROC) = 'HTURB_'
  IPROC=IPROC+1
  YWORK2(9,IPROC) = 'NEGA_'
  IPROC=  IPROC+1
  YWORK2(9,IPROC) = 'HENU_'
  IPROC=  IPROC+1
  YWORK2(9,IPROC) = 'HON_'
  IPROC=  IPROC+1
  YWORK2(9,IPROC) = 'AGGS_'
  IPROC=  IPROC+1
  YWORK2(9,IPROC) = 'AUTS_'
  IPROC=  IPROC+1
  YWORK2(9,IPROC) = 'CFRZ_'
  IPROC=  IPROC+1
  YWORK2(9,IPROC) = 'WETG_'
  IPROC=  IPROC+1
  YWORK2(9,IPROC) = 'DRYG_'
  IPROC=  IPROC+1
  YWORK2(9,IPROC) = 'WETH_'
  IPROC=  IPROC+1
  YWORK2(9,IPROC) = 'IMLT_'
  IPROC=  IPROC+1
  YWORK2(9,IPROC) = 'BERFI_'
  IPROC=  IPROC+1
  YWORK2(9,IPROC) = 'SEDI_'
  IPROC=  IPROC+1
  YWORK2(9,IPROC) = 'CDEPI_'
!
  YEND_COMMENT(9) = 'BU_RRI'
  NBUPROCNBR(9) = 3
!
  CBUACTION(9,1) = 'IG'
  CBUACTION(9,2) = 'CC'
  CBUACTION(9,3) = 'ES'
!
  DO JJ=1,3
    CBUCOMMENT(9,JJ) = ADJUSTL( ADJUSTR( YWORK2(9,JJ) ) // &
                                ADJUSTL( YEND_COMMENT(9) ) )
  END DO
!
END IF
!
!                        Budget of RRS
IF (LBU_RRS) THEN
  IPROC=4
  IPROACTV(10,IPROC) = NASSERS
  IPROC=  IPROC+1
  IF( NMODEL>1 ) IPROACTV(10,IPROC) = NNESTRS
  IPROC=  IPROC+1
  IF ( HMET_ADV_SCHEME(1:3) == 'PPM' ) &
     IPROACTV(10,IPROC) = NADVRS
  IPROC=IPROC+1
  IF ( HMET_ADV_SCHEME(1:3) /= 'PPM' ) &
     IPROACTV(10,IPROC) = NADVXRS
  IPROC=IPROC+1
  IF ( HMET_ADV_SCHEME(1:3) /= 'PPM' ) &
     IPROACTV(10,IPROC) = NADVYRS
  IPROC=IPROC+1
  IF ( HMET_ADV_SCHEME(1:3) /= 'PPM' ) &
     IPROACTV(10,IPROC) = NADVZRS
  IPROC=  IPROC+1
  IF( LFORCING )  IPROACTV(10,IPROC)  = NFRCRS
  IPROC=  IPROC+1
  IF( ONUMDIFTH ) IPROACTV(10,IPROC) = NDIFRS
  IPROC=IPROC+1
  IF ( OHORELAX_RS ) THEN
    IPROACTV(10,IPROC) = NRELRS
  ELSE
    IF(OVE_RELAX .OR. OHORELAX_UVWTH .OR. OHORELAX_RV .OR.                 &
     OHORELAX_RC .OR. OHORELAX_RR .OR. OHORELAX_RI .OR. OHORELAX_RS .OR.   &
     OHORELAX_RG .OR. OHORELAX_RH .OR. OHORELAX_TKE .OR. ANY(OHORELAX_SV)) THEN
      IPROACTV(10,IPROC) = 4
    ELSE
      IPROACTV(10,IPROC) = 3
    END IF
  END IF
  IPROC=  IPROC+1
  IF( HCLOUD /= 'NONE' ) IPROACTV(10,IPROC) = NNEGARS
  IPROC=IPROC+1
  IF( HCLOUD(1:3) == 'ICE') IPROACTV(10,IPROC) = NDEPSRS
  IPROC=IPROC+1
  IF( HCLOUD(1:3) == 'ICE') IPROACTV(10,IPROC) = NAGGSRS
  IPROC=IPROC+1
  IF( HCLOUD(1:3) == 'ICE') IPROACTV(10,IPROC) = NAUTSRS
  IPROC=IPROC+1
  IF( HCLOUD(1:3) == 'ICE') IPROACTV(10,IPROC) = NRIMRS
  IPROC=IPROC+1
  IF( HCLOUD(1:3) == 'ICE') IPROACTV(10,IPROC) = NACCRS
  IPROC=IPROC+1
  IF( HCLOUD(1:3) == 'ICE') IPROACTV(10,IPROC) = NCMELRS
  IPROC=IPROC+1
  IF( HCLOUD(1:3) == 'ICE') IPROACTV(10,IPROC) = NWETGRS
  IPROC=IPROC+1
  IF( HCLOUD(1:3) == 'ICE') IPROACTV(10,IPROC) = NDRYGRS
  IPROC=IPROC+1
  IF( HCLOUD == 'ICE4') IPROACTV(10,IPROC) = NWETHRS
  IPROC=IPROC+1
  IF( HCLOUD(1:3) == 'ICE') IPROACTV(10,IPROC) = NSEDIRS
  IPROC=IPROC+1
!
  YWORK2(10,1) = 'INIF_'
  YWORK2(10,2) = 'ENDF_'
  YWORK2(10,3) = 'AVEF_'
  IPROC= 4
  YWORK2(10,IPROC) = 'ASSE_'
  IPROC=  IPROC+1
  YWORK2(10,IPROC) = 'NEST_'
  IPROC=  IPROC+1
  YWORK2(10,IPROC) = 'ADV_'
  IPROC=  IPROC+1
  YWORK2(10,IPROC) = 'ADVX_'
  IPROC=  IPROC+1
  YWORK2(10,IPROC) = 'ADVY_'
  IPROC=  IPROC+1
  YWORK2(10,IPROC) = 'ADVZ_'
  IPROC=  IPROC+1
  YWORK2(10,IPROC) = 'FRC_'
  IPROC=  IPROC+1
  YWORK2(10,IPROC) = 'DIF_'
  IPROC=  IPROC+1
  YWORK2(10,IPROC) = 'REL_'
  IPROC=IPROC+1
  YWORK2(10,IPROC) = 'NEGA_'
  IPROC=  IPROC+1
  YWORK2(10,IPROC) = 'DEPS_'
  IPROC=  IPROC+1
  YWORK2(10,IPROC) = 'AGGS_'
  IPROC=  IPROC+1
  YWORK2(10,IPROC) = 'AUTS_'
  IPROC=  IPROC+1
  YWORK2(10,IPROC) = 'RIM_'
  IPROC=  IPROC+1
  YWORK2(10,IPROC) = 'ACC_'
  IPROC=  IPROC+1
  YWORK2(10,IPROC) = 'CMEL_'
  IPROC=  IPROC+1
  YWORK2(10,IPROC) = 'WETG_'
  IPROC=  IPROC+1
  YWORK2(10,IPROC) = 'DRYG_'
  IPROC=  IPROC+1
  YWORK2(10,IPROC) = 'WETH_'
  IPROC=IPROC+1
  YWORK2(10,IPROC) = 'SEDI_'
!
  YEND_COMMENT(10) = 'BU_RRS'
  NBUPROCNBR(10) = 3
!
  CBUACTION(10,1) = 'IG'
  CBUACTION(10,2) = 'CC'
  CBUACTION(10,3) = 'ES'
!
  DO JJ=1,3
    CBUCOMMENT(10,JJ) = ADJUSTL( ADJUSTR( YWORK2(10,JJ) ) // &
                                 ADJUSTL( YEND_COMMENT(10) ) )
  END DO
!
END IF
!
!                        Budget of RRG
IF (LBU_RRG) THEN
  IPROC=4
  IPROACTV(11,IPROC) = NASSERG
  IPROC=IPROC+1
  IF( NMODEL>1 ) IPROACTV(11,IPROC) = NNESTRG
  IPROC=IPROC+1
  IF ( HMET_ADV_SCHEME(1:3) == 'PPM' ) &
     IPROACTV(11,IPROC) = NADVRG
  IPROC=IPROC+1
  IF ( HMET_ADV_SCHEME(1:3) /= 'PPM' ) &
     IPROACTV(11,IPROC) = NADVXRG
  IPROC=IPROC+1
  IF ( HMET_ADV_SCHEME(1:3) /= 'PPM' ) &
     IPROACTV(11,IPROC) = NADVYRG
  IPROC=IPROC+1
  IF ( HMET_ADV_SCHEME(1:3) /= 'PPM' ) &
     IPROACTV(11,IPROC) = NADVZRG
  IPROC=IPROC+1
  IF( LFORCING ) IPROACTV(11,IPROC)  = NFRCRG
  IPROC=IPROC+1
  IF( ONUMDIFTH ) IPROACTV(11,IPROC) = NDIFRG
  IPROC=IPROC+1
  IF ( OHORELAX_RG ) THEN
    IPROACTV(11,IPROC) = NRELRG
  ELSE
    IF(OVE_RELAX .OR. OHORELAX_UVWTH .OR. OHORELAX_RV .OR.                 &
     OHORELAX_RC .OR. OHORELAX_RR .OR. OHORELAX_RI .OR. OHORELAX_RS .OR.   &
     OHORELAX_RG .OR. OHORELAX_RH .OR. OHORELAX_TKE .OR. ANY(OHORELAX_SV)) THEN
      IPROACTV(11,IPROC) = 4
    ELSE
      IPROACTV(11,IPROC) = 3
    END IF
  END IF
  IPROC=IPROC+1
  IF( HCLOUD /= 'NONE'  ) IPROACTV(11,IPROC) = NNEGARG
  IPROC=IPROC+1
  IF( HCLOUD(1:3) == 'ICE') IPROACTV(11,IPROC) = NSFRRG
  IPROC=IPROC+1
  IF( HCLOUD(1:3) == 'ICE') IPROACTV(11,IPROC) = NDEPGRG
  IPROC=IPROC+1
  IF( HCLOUD(1:3) == 'ICE') IPROACTV(11,IPROC) = NRIMRG
  IPROC=IPROC+1
  IF( HCLOUD(1:3) == 'ICE') IPROACTV(11,IPROC) = NACCRG
  IPROC=IPROC+1
  IF( HCLOUD(1:3) == 'ICE') IPROACTV(11,IPROC) = NCMELRG
  IPROC=IPROC+1
  IF( HCLOUD(1:3) == 'ICE') IPROACTV(11,IPROC) = NCFRZRG
  IPROC=IPROC+1
  IF( HCLOUD(1:3) == 'ICE') IPROACTV(11,IPROC) = NWETGRG
  IPROC=IPROC+1
  IF( HCLOUD(1:3) == 'ICE') IPROACTV(11,IPROC) = NDRYGRG
  IPROC=IPROC+1
  IF( HCLOUD(1:3) == 'ICE') IPROACTV(11,IPROC) = NGMLTRG
  IPROC=IPROC+1
  IF( HCLOUD == 'ICE4') IPROACTV(11,IPROC) = NWETHRG
  IPROC=IPROC+1
  IF( HCLOUD(1:3) == 'ICE')IPROACTV(11,IPROC) = NSEDIRG
  IPROC=IPROC+1
!
  YWORK2(11,1) = 'INIF_'
  YWORK2(11,2) = 'ENDF_'
  YWORK2(11,3) = 'AVEF_'
  IPROC=4
  YWORK2(11,IPROC) = 'ASSE_'
  IPROC=IPROC+1
  YWORK2(11,IPROC) = 'NEST_'
  IPROC=IPROC+1
  YWORK2(11,IPROC) = 'ADV_'
  IPROC=IPROC+1
  YWORK2(11,IPROC) = 'ADVX_'
  IPROC=IPROC+1
  YWORK2(11,IPROC) = 'ADVY_'
  IPROC=IPROC+1
  YWORK2(11,IPROC) = 'ADVZ_'
  IPROC=IPROC+1
  YWORK2(11,IPROC) = 'FRC_'
  IPROC=IPROC+1
  YWORK2(11,IPROC) = 'DIF_'
  IPROC=IPROC+1
  YWORK2(11,IPROC) = 'REL_'
  IPROC=IPROC+1
  YWORK2(11,IPROC) = 'NEGA_'
  IPROC=IPROC+1
  YWORK2(11,IPROC)= 'SFR_'
  IPROC=IPROC+1
  YWORK2(11,IPROC) = 'DEPG_'
  IPROC=IPROC+1
  YWORK2(11,IPROC) = 'RIM_'
  IPROC=IPROC+1
  YWORK2(11,IPROC) = 'ACC_'
  IPROC=IPROC+1
  YWORK2(11,IPROC) = 'CMEL_'
  IPROC=IPROC+1
  YWORK2(11,IPROC) = 'CFRZ_'
  IPROC=IPROC+1
  YWORK2(11,IPROC) = 'WETG_'
  IPROC=IPROC+1
  YWORK2(11,IPROC) = 'DRYG_'
  IPROC=IPROC+1
  YWORK2(11,IPROC) = 'GMLT_'
  IPROC=IPROC+1
  YWORK2(11,IPROC) = 'WETH_'
  IPROC=IPROC+1
  YWORK2(11,IPROC)= 'SEDI_'
!
  YEND_COMMENT(11) = 'BU_RRG'
  NBUPROCNBR(11) = 3
!
  CBUACTION(11,1) = 'IG'
  CBUACTION(11,2) = 'CC'
  CBUACTION(11,3) = 'ES'
!
  DO JJ=1,3
    CBUCOMMENT(11,JJ) = ADJUSTL( ADJUSTR( YWORK2(11,JJ) ) // &
                                 ADJUSTL( YEND_COMMENT(11) ) )
  END DO
!
END IF
!
!                        Budget of RRH
IF (LBU_RRH) THEN
  IPROC=4
  IPROACTV(12,IPROC) = NASSERH
  IPROC=IPROC+1
  IF( NMODEL>1 ) THEN
    IPROACTV(12,IPROC) = NNESTRH
  ELSE
    IPROACTV(12,IPROC) = 3
  END IF
  IPROC=IPROC+1
  IF ( HMET_ADV_SCHEME(1:3) == 'PPM' ) &
     IPROACTV(12,IPROC) = NADVRH
  IPROC=IPROC+1
  IF ( HMET_ADV_SCHEME(1:3) /= 'PPM' ) &
     IPROACTV(12,IPROC) = NADVXRH
  IPROC=IPROC+1
  IF ( HMET_ADV_SCHEME(1:3) /= 'PPM' ) &
     IPROACTV(12,IPROC) = NADVYRH
  IPROC=IPROC+1
  IF ( HMET_ADV_SCHEME(1:3) /= 'PPM' ) &
     IPROACTV(12,IPROC) = NADVZRH
  IPROC=IPROC+1
   IF( LFORCING ) THEN
    IPROACTV(12,IPROC)  = NFRCRH
  ELSE
    IPROACTV(12,IPROC)  = 3
  END IF
  IPROC=IPROC+1
  IF( ONUMDIFTH ) THEN
    IPROACTV(12,IPROC) = NDIFRH
  ELSE
    IPROACTV(12,IPROC) = 3
  END IF
  IPROC=IPROC+1
  IF ( OHORELAX_RH ) THEN
    IPROACTV(12,IPROC) = NRELRH
  ELSE
    IF(OVE_RELAX .OR. OHORELAX_UVWTH .OR. OHORELAX_RV .OR.                 &
     OHORELAX_RC .OR. OHORELAX_RR .OR. OHORELAX_RI .OR. OHORELAX_RS .OR.   &
     OHORELAX_RG .OR. OHORELAX_RH .OR. OHORELAX_TKE .OR. ANY(OHORELAX_SV)) THEN
      IPROACTV(12,IPROC) = 4
    ELSE
      IPROACTV(12,IPROC) = 3
    END IF
  END IF
  IPROC=IPROC+1
  IF( HCLOUD /= 'NONE' ) THEN
    IPROACTV(12,IPROC) = NNEGARH
  ELSE
    IPROACTV(12,IPROC) = 3
  END IF
  IPROC=IPROC+1
  IF( HCLOUD == 'ICE4') IPROACTV(12,IPROC) = NWETGRH
  IPROC=IPROC+1
  IF( HCLOUD == 'ICE4') IPROACTV(12,IPROC) = NWETHRH
  IPROC=IPROC+1
  IF( HCLOUD == 'ICE4') IPROACTV(12,IPROC) = NHMLTRH
  IPROC=IPROC+1
  IF( HCLOUD == 'ICE4') IPROACTV(12,IPROC) = NSEDIRH
!
  YWORK2(12,1) = 'INIF_'
  YWORK2(12,2) = 'ENDF_'
  YWORK2(12,3) = 'AVEF_'
  IPROC=4
  YWORK2(12,IPROC) = 'ASSE_'
  IPROC=IPROC+1
  YWORK2(12,IPROC) = 'NEST_'
  IPROC=IPROC+1
  YWORK2(12,IPROC) = 'ADV_'
  IPROC=IPROC+1
  YWORK2(12,IPROC) = 'ADVX_'
  IPROC=IPROC+1
  YWORK2(12,IPROC) = 'ADVY_'
  IPROC=IPROC+1
  YWORK2(12,IPROC) = 'ADVZ_'
  IPROC=IPROC+1
  YWORK2(12,IPROC) = 'FRC_'
  IPROC=IPROC+1
  YWORK2(12,IPROC) = 'DIF_'
  IPROC=IPROC+1
  YWORK2(12,IPROC) = 'REL_'
  IPROC=IPROC+1
  YWORK2(12,IPROC) = 'NEGA_'
  IPROC=IPROC+1
  YWORK2(12,IPROC) = 'WETG_'
  IPROC=IPROC+1
  YWORK2(12,IPROC) = 'WETH_'
  IPROC=IPROC+1
  YWORK2(12,IPROC) = 'HMLT_'
  IPROC=IPROC+1
  YWORK2(12,IPROC) = 'SEDI_'
!
  YEND_COMMENT(12) = 'BU_RRH'
  NBUPROCNBR(12) = 3
!
  CBUACTION(12,1) = 'IG'
  CBUACTION(12,2) = 'CC'
  CBUACTION(12,3) = 'ES'
!
  DO JJ=1,3
    CBUCOMMENT(12,JJ) = ADJUSTL( ADJUSTR( YWORK2(12,JJ) ) // &
                                 ADJUSTL( YEND_COMMENT(12) ) )
  END DO
!

END IF
!
!                        Budget of RSV
IF (LBU_RSV) THEN
  IBUPROCNBR_SV_MAX = 0 ! initialize the Max nunmber of processes for the SVs
  DO JSV = 1,KSV
    IPROC=4
    IPROACTV(12+JSV,IPROC) = NASSESV
    IPROC=IPROC+1
    IF( NMODEL>1 ) IPROACTV(12+JSV,IPROC) = NNESTSV
    IPROC=IPROC+1
    IF ( HSV_ADV_SCHEME(1:3) == 'PPM' ) &
     IPROACTV(12+JSV,IPROC) = NADVSV
    IPROC=IPROC+1
    IF ( HSV_ADV_SCHEME(1:3) /= 'PPM' ) &
     IPROACTV(12+JSV,IPROC) = NADVXSV
    IPROC=IPROC+1
    IF ( HSV_ADV_SCHEME(1:3) /= 'PPM' ) &
     IPROACTV(12+JSV,IPROC) = NADVYSV
    IPROC=IPROC+1
    IF ( HSV_ADV_SCHEME(1:3) /= 'PPM' ) &
     IPROACTV(12+JSV,IPROC) = NADVZSV
    IPROC=IPROC+1
    IF( LFORCING ) IPROACTV(12+JSV,IPROC)  = NFRCSV
    IPROC=IPROC+1
    IF ( ONUMDIFSV ) IPROACTV(12+JSV,IPROC) = NDIFSV
    IPROC=IPROC+1
    IF ( OHORELAX_SV(JSV) ) THEN
      IPROACTV(12+JSV,IPROC) = NRELSV
    ELSE
    IF(OVE_RELAX .OR. OHORELAX_UVWTH .OR. OHORELAX_RV .OR.                 &
     OHORELAX_RC .OR. OHORELAX_RR .OR. OHORELAX_RI .OR. OHORELAX_RS .OR.   &
     OHORELAX_RG .OR. OHORELAX_RH .OR. OHORELAX_TKE .OR. ANY(OHORELAX_SV)) THEN
        IPROACTV(12+JSV,IPROC) = 4
      ELSE
        IPROACTV(12+JSV,IPROC) = 3
      END IF
    END IF
    IPROC=IPROC+1
    IF ( (HDCONV /= 'NONE' .OR. HSCONV == 'KAFR') .AND. OCHTRANS ) &
        IPROACTV(12+JSV,IPROC) = NDCONVSV
    IPROC=IPROC+1
    IF ( HTURB /= 'NONE' ) IPROACTV(12+JSV,IPROC) = NVTURBSV
    IPROC=IPROC+1
    IF ( HTURB /= 'NONE' .AND. HTURBDIM == '3DIM' ) THEN
      IPROACTV(12+JSV,IPROC) = NHTURBSV
    ELSE
      IF ( HTURB /= 'NONE' ) THEN
        IPROACTV(12+JSV,IPROC) = 4
      ELSE
        IPROACTV(12+JSV,IPROC) = 3
      END IF
    END IF
    IPROC=IPROC+1
    IF ( HSCONV == 'EDKF' ) IPROACTV(12+JSV,IPROC)= NMAFLSV
    IPROC=IPROC+1
!
    YWORK2(12+JSV,1) = 'INIF_'
    YWORK2(12+JSV,2) = 'ENDF_'
    YWORK2(12+JSV,3) = 'AVEF_'
    IPROC=4
    YWORK2(12+JSV,IPROC) = 'ASSE_'
    IPROC=IPROC+1
    YWORK2(12+JSV,IPROC) = 'NEST_'
    IPROC=IPROC+1
    YWORK2(12+JSV,IPROC) = 'ADV_'
    IPROC=IPROC+1
    YWORK2(12+JSV,IPROC) = 'ADVX_'
    IPROC=IPROC+1
    YWORK2(12+JSV,IPROC) = 'ADVY_'
    IPROC=IPROC+1
    YWORK2(12+JSV,IPROC) = 'ADVZ_'
    IPROC=IPROC+1
    YWORK2(12+JSV,IPROC) = 'FRC_'
    IPROC=IPROC+1
    YWORK2(12+JSV,IPROC) = 'DIF_'
    IPROC=IPROC+1
    YWORK2(12+JSV,IPROC) = 'REL_'
    IPROC=IPROC+1
    YWORK2(12+JSV,IPROC) = 'DCONV_'
    IPROC=IPROC+1
    YWORK2(12+JSV,IPROC) = 'VTURB_'
    IPROC=IPROC+1
    YWORK2(12+JSV,IPROC) = 'HTURB_'
    IPROC=IPROC+1
    YWORK2(12+JSV,IPROC) = 'MAFL_'
!
! complete with the budget of other processes
!
    ILAST_PROC_NBR = IPROC
    CALL BUDGET_OTHERPROC_SV
!
    YEND_COMMENT(12+JSV) = 'BU_RSV'
    IBUPROCNBR_SV_MAX   = MAX( IBUPROCNBR_SV_MAX, ILAST_PROC_NBR )
    NBUPROCNBR(12+JSV) = 3
!
    CBUACTION(12+JSV,1) = 'IG'
    CBUACTION(12+JSV,2) = 'CC'
    CBUACTION(12+JSV,3) = 'ES'
!
    DO JJ=1,3
      CBUCOMMENT(12+JSV,JJ) = ADJUSTL( ADJUSTR( YWORK2(12+JSV,JJ) ) // &
                                       ADJUSTL( YEND_COMMENT(12+JSV) ) )
    END DO
  END DO
!
END IF

!-------------------------------------------------------------------------------
!*       4.    COMPUTE THE INCREMENT BETWEEN TWO ACTIVE SOURCES
!              ------------------------------------------------
!
NBUINC(:,:) = 1
DO JI = 1, JPBUMAX
  DO JJ = 4, JPBUPROMAX-1
    DO JK = JJ+1,JPBUPROMAX
      IF ( IPROACTV(JI,JK) /= 3 ) EXIT
      NBUINC(JI,JJ) = NBUINC(JI,JJ) +1
    END DO
  END DO
END DO
!
!-------------------------------------------------------------------------------
!*       5.    COMPUTE PROCESSES ACTIONS AND NAMES OF BUDGET OUTPUT ARRAYS
!              -----------------------------------------------------------
!
DO JI=1,JPBUMAX                                ! loop on the allowed budgets
                                               ! names of recording files for:
  CBURECORD(JI,1) = ADJUSTL( CBUCOMMENT(JI,1) )   ! initial guess
  CBURECORD(JI,2) = ADJUSTL( CBUCOMMENT(JI,2) )   ! source cumul
  CBURECORD(JI,3) = ADJUSTL( CBUCOMMENT(JI,3) )   ! end step
!
  IF (IPROACTV(JI,4) >= 2) THEN
    WRITE(UNIT=KLUOUT,FMT= '("Error in budget specification of ",A7,/," &
    & The first source either is the first element of a group of sources or &
    & is not considered")')  YEND_COMMENT(JI)
    WRITE(UNIT=KLUOUT,FMT= '("change this namelist element ")')
    GERROR = .TRUE.
  END IF
!
  DO JJ=4,JPBUPROMAX                           ! loop on the allowed processes
    IF (IPROACTV(JI,JJ) == 0) THEN
      IF(IPROACTV(JI,JJ+NBUINC(JI,JJ)) == 0) THEN
        CBUACTION(JI,JJ)='OF'
      ELSE IF (IPROACTV(JI,JJ+NBUINC(JI,JJ)) == 1) THEN
        CBUACTION(JI,JJ)='CC'
      ELSE IF (IPROACTV(JI,JJ+NBUINC(JI,JJ)) == 2) THEN
        WRITE(UNIT=KLUOUT,FMT= '("Error in budget specification of ",A15)') &
        ADJUSTL( ADJUSTR(YWORK2(JI,JJ+NBUINC(JI,JJ)))//ADJUSTL(YEND_COMMENT(JI)))
        WRITE(UNIT=KLUOUT,FMT= '("change this namelist ")')
        GERROR = .TRUE.
      END IF
    ELSE IF (IPROACTV(JI,JJ) <= 2) THEN
      DO JJJ = JJ+NBUINC(JI,JJ), JPBUPROMAX
         IF(IPROACTV(JI,JJJ) /= 3 .AND. IPROACTV(JI,JJJ) /= 4) EXIT
      END DO
!
      IF (IPROACTV(JI,JJJ) == 1) THEN
        NBUPROCNBR(JI) = NBUPROCNBR(JI)+1
        CBUACTION(JI,JJ) = 'DC'
        CBUCOMMENT(JI,NBUPROCNBR(JI)) =                            ADJUSTL(    &
                                   ADJUSTR( CBUCOMMENT(JI,NBUPROCNBR(JI)) ) // &
                                   ADJUSTL( ADJUSTR( YWORK2(JI,JJ) ) //        &
                                            ADJUSTL( YEND_COMMENT(JI) ) ) )
        ITEN=INT(NBUPROCNBR(JI)/10)
        CBURECORD(JI,NBUPROCNBR(JI)) = 'S' // CHAR( ITEN + 48 )               &
                  // CHAR(  48+ MODULO( NBUPROCNBR(JI),10*MAX(1,ITEN) )  )    &
                  // '_' // ADJUSTL( YEND_COMMENT(JI) )
      ELSE IF (IPROACTV(JI,JJJ) == 0) THEN
        NBUPROCNBR(JI) = NBUPROCNBR(JI)+1
        CBUACTION(JI,JJ) = 'DD'
        CBUCOMMENT(JI,NBUPROCNBR(JI)) =                            ADJUSTL(    &
                                   ADJUSTR( CBUCOMMENT(JI,NBUPROCNBR(JI)) ) // &
                                   ADJUSTL( ADJUSTR( YWORK2(JI,JJ) ) //        &
                                            ADJUSTL( YEND_COMMENT(JI) ) ) )
        ITEN=INT(NBUPROCNBR(JI)/10)
        CBURECORD(JI,NBUPROCNBR(JI)) = 'S' // CHAR( ITEN + 48 )               &
                  // CHAR(  48+ MODULO( NBUPROCNBR(JI),10*MAX(1,ITEN) )  )    &
                  // '_' // ADJUSTL( YEND_COMMENT(JI) )
      ELSE IF (IPROACTV(JI,JJJ) == 2) THEN
        CBUACTION(JI,JJ) = 'NO'
        CBUCOMMENT(JI,NBUPROCNBR(JI)+1) =           ADJUSTL(                   &
                                  ADJUSTR( CBUCOMMENT(JI,NBUPROCNBR(JI)+1)) // &
                                  ADJUSTL( YWORK2(JI,JJ) ) )
      END IF
    ELSEIF (IPROACTV(JI,JJ) == 3) THEN
      CBUACTION(JI,JJ) = 'RM'
    ELSEIF (IPROACTV(JI,JJ) == 4) THEN
      CBUACTION(JI,JJ) = 'OF'
    ELSE
      WRITE(UNIT=KLUOUT,FMT= '("Error in budget specification of ",A7)') &
            YEND_COMMENT(JI)
      WRITE(UNIT=KLUOUT,FMT= '("change this namelist ")')
      GERROR = .TRUE.
    END IF
  END DO
END DO
!   writes on output the explicit chain of sources for all the budgets
DO JI=1,JPBUMAX                            ! loop over the allowed budgets
  YSTRING = ADJUSTL( YEND_COMMENT(JI) )
  ILEN    = LEN_TRIM(YSTRING)
  IF( ILEN /= 0 ) THEN
    IF( JI <= 12 ) THEN
      WRITE (UNIT=KLUOUT,FMT='(/,"budget ",A7," with ",I2," vectors")')        &
                                                  YSTRING(1:ILEN),NBUPROCNBR(JI)
      DO JJ=1,3
        YSTRING = CBUCOMMENT(JI,JJ)
        ILEN    = LEN_TRIM(YSTRING)
        WRITE (UNIT=KLUOUT,FMT='(12X,A40)')         YSTRING(1:ILEN)
      END DO
      DO JJ=4,NBUPROCNBR(JI)                  ! loop over the allowed processes
        YSTRING = CBUCOMMENT(JI,JJ)
        ILEN    = LEN_TRIM(YSTRING)
        WRITE (UNIT=KLUOUT,FMT='(20X,A40)')         YSTRING(1:ILEN)
      END DO
    ELSE
      WRITE (UNIT=KLUOUT,                                                      &
             FMT='(/,"budget ",A7," (number ",I3,") with ",I2," vectors")')    &
                                            YSTRING(1:ILEN),JI-12,NBUPROCNBR(JI)
      DO JJ=1,3
        YSTRING = CBUCOMMENT(JI,JJ)
        ILEN    = LEN_TRIM(YSTRING)
        WRITE (UNIT=KLUOUT,FMT='(12X,A40)')         YSTRING(1:ILEN)
      END DO
      DO JJ=4,NBUPROCNBR(JI)                  ! loop over the allowed processes
        YSTRING = CBUCOMMENT(JI,JJ)
        ILEN    = LEN_TRIM(YSTRING)
        WRITE (UNIT=KLUOUT,FMT='(20X,A40)')         YSTRING(1:ILEN)
      END DO
    END IF
  END IF
END DO
!
IF (CBUTYPE=='CART') THEN
  WRITE(UNIT=KLUOUT, FMT= '(2/,"DESCRIPTION OF THE BUDGET BOX")' )
  WRITE(UNIT=KLUOUT, FMT= '("BUIL = ",I4.4)' ) NBUIL
  WRITE(UNIT=KLUOUT, FMT= '("BUIH = ",I4.4)' ) NBUIH
  WRITE(UNIT=KLUOUT, FMT= '("BUJL = ",I4.4)' ) NBUJL
  WRITE(UNIT=KLUOUT, FMT= '("BUJH = ",I4.4)' ) NBUJH
  WRITE(UNIT=KLUOUT, FMT= '("BUKL = ",I4.4)' ) NBUKL
  WRITE(UNIT=KLUOUT, FMT= '("BUKH = ",I4.4)' ) NBUKH
  WRITE(UNIT=KLUOUT, FMT= '("BUIMAX = ",I4.4)' ) NBUIMAX
  WRITE(UNIT=KLUOUT, FMT= '("BUJMAX = ",I4.4)' ) NBUJMAX
  WRITE(UNIT=KLUOUT, FMT= '("BUKMAX = ",I4.4)' ) NBUKMAX
END IF
IF (CBUTYPE=='MASK') THEN
  WRITE(UNIT=KLUOUT, FMT= '(2/,"DESCRIPTION OF THE BUDGET MASK")' )
  WRITE(UNIT=KLUOUT, FMT= '("BUIL = ",I4.4)' ) NBUIL
  WRITE(UNIT=KLUOUT, FMT= '("BUIH = ",I4.4)' ) NBUIH
  WRITE(UNIT=KLUOUT, FMT= '("BUJL = ",I4.4)' ) NBUJL
  WRITE(UNIT=KLUOUT, FMT= '("BUJH = ",I4.4)' ) NBUJH
  WRITE(UNIT=KLUOUT, FMT= '("BUKL = ",I4.4)' ) NBUKL
  WRITE(UNIT=KLUOUT, FMT= '("BUKH = ",I4.4)' ) NBUKH
  WRITE(UNIT=KLUOUT, FMT= '("BUKMAX = ",I4.4)' ) NBUKMAX
  WRITE(UNIT=KLUOUT, FMT= '("BUWRNB = ",I4.4)' ) NBUWRNB
  WRITE(UNIT=KLUOUT, FMT= '("BUMASK = ",I4.4)' ) NBUMASK
END IF
IF (GERROR) THEN
   !callabortstop
  !CALL CLOSE_ll(HLUOUT,IOSTAT=IRESP)
  CALL ABORT
  STOP
ENDIF
!-------------------------------------------------------------------------------
!*       5.    ALLOCATE MEMORY FOR BUDGET STORAGE ARRAYS
!              -----------------------------------------
IF (LBU_RU) THEN
  ALLOCATE ( XBURU(IBUDIM1, IBUDIM2, IBUDIM3, NBUPROCNBR(1)) )
  XBURU(:,:,:,:)=0.
  ALLOCATE ( XBURHODJU(IBUDIM1, IBUDIM2, IBUDIM3) )
  XBURHODJU(:,:,:)=0.
END IF
!
IF (LBU_RV) THEN
  ALLOCATE ( XBURV(IBUDIM1, IBUDIM2, IBUDIM3, NBUPROCNBR(2)) )
  XBURV(:,:,:,:)=0.
  ALLOCATE ( XBURHODJV(IBUDIM1, IBUDIM2, IBUDIM3) )
  XBURHODJV(:,:,:)=0.
END IF
!
IF (LBU_RW) THEN
  ALLOCATE ( XBURW(IBUDIM1, IBUDIM2, IBUDIM3, NBUPROCNBR(3)) )
  XBURW(:,:,:,:)=0.
  ALLOCATE ( XBURHODJW(IBUDIM1, IBUDIM2, IBUDIM3) )
  XBURHODJW(:,:,:)=0.
END IF
!
IF (LBU_RTH .OR. LBU_RTKE .OR. LBU_RRV .OR. LBU_RRC .OR. LBU_RRR .OR. &
    LBU_RRI .OR. LBU_RRS  .OR. LBU_RRG .OR. LBU_RRH .OR. LBU_RSV      ) THEN
  ALLOCATE ( XBURHODJ(IBUDIM1, IBUDIM2, IBUDIM3) )
  XBURHODJ(:,:,:)=0.
END IF
!
IF (LBU_RTH) THEN
  ALLOCATE ( XBURTH(IBUDIM1, IBUDIM2, IBUDIM3, NBUPROCNBR(4)) )
  XBURTH(:,:,:,:)=0.
END IF
!
IF (LBU_RTKE) THEN
  ALLOCATE ( XBURTKE(IBUDIM1, IBUDIM2, IBUDIM3, NBUPROCNBR(5)) )
  XBURTKE(:,:,:,:)=0.
END IF
!
IF (LBU_RRV) THEN
  ALLOCATE ( XBURRV(IBUDIM1, IBUDIM2, IBUDIM3, NBUPROCNBR(6)) )
  XBURRV(:,:,:,:)=0.
END IF
!
IF (LBU_RRC) THEN
  ALLOCATE ( XBURRC(IBUDIM1, IBUDIM2, IBUDIM3, NBUPROCNBR(7)) )
  XBURRC(:,:,:,:)=0.
END IF
!
IF (LBU_RRR) THEN
  ALLOCATE ( XBURRR(IBUDIM1, IBUDIM2, IBUDIM3, NBUPROCNBR(8)) )
  XBURRR(:,:,:,:)=0.
END IF
!
IF (LBU_RRI) THEN
  ALLOCATE ( XBURRI(IBUDIM1, IBUDIM2, IBUDIM3, NBUPROCNBR(9)) )
  XBURRI(:,:,:,:)=0.
END IF
!
IF (LBU_RRS) THEN
  ALLOCATE ( XBURRS(IBUDIM1, IBUDIM2, IBUDIM3, NBUPROCNBR(10)) )
  XBURRS(:,:,:,:)=0.
END IF
!
IF (LBU_RRG) THEN
  ALLOCATE ( XBURRG(IBUDIM1, IBUDIM2, IBUDIM3, NBUPROCNBR(11)) )
  XBURRG(:,:,:,:)=0.
END IF
!
IF (LBU_RRH) THEN
  ALLOCATE ( XBURRH(IBUDIM1, IBUDIM2, IBUDIM3, NBUPROCNBR(12)) )
  XBURRH(:,:,:,:)=0.
END IF
!
IF (LBU_RSV) THEN
  ALLOCATE ( XBURSV(IBUDIM1, IBUDIM2, IBUDIM3, IBUPROCNBR_SV_MAX, KSV) )
  XBURSV(:,:,:,:,:)=0.
END IF
!
IF (LHOOK) CALL DR_HOOK('INI_BUDGET',1,ZHOOK_HANDLE)
CONTAINS
! ##############################
  SUBROUTINE BUDGET_OTHERPROC_SV
! ##############################
!
!
USE MODD_NSV, ONLY : NSV_USER, NSV_C2R2BEG, NSV_C2R2END, &
                               NSV_CHEMBEG, NSV_CHEMEND, &
                               NSV_ELECBEG, NSV_ELECEND
!
  REAL(KIND=JPRB) :: ZHOOK_HANDLE
  IF (LHOOK) CALL DR_HOOK('INI_BUDGET:BUDGET_OTHERPROC_SV',0,ZHOOK_HANDLE)
  IF (JSV <= NSV_USER) THEN
    ! NSV_USER Case
    SELECT CASE(JSV)
    CASE (1)
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR)= 'PROC1_'
      IPROACTV(12+JSV,ILAST_PROC_NBR) = 1
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR)= 'PROC2_'
      IPROACTV(12+JSV,ILAST_PROC_NBR) = 1
    CASE (2)
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR)= 'PROC3_'
      IPROACTV(12+JSV,ILAST_PROC_NBR) = 1
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR)= 'PROC4_'
      IPROACTV(12+JSV,ILAST_PROC_NBR) = 1
    END SELECT
    !
  ELSEIF (JSV >= NSV_C2R2BEG .AND. JSV <= NSV_C2R2END) THEN
    ! C2R2 or KHKO Case
    SELECT CASE(JSV-NSV_C2R2BEG+1)
    CASE (1)                               ! Concentration of activated nuclei
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR)= 'HENU_'
      IPROACTV(12+JSV,ILAST_PROC_NBR) = 1
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR)= 'CEVA_'
      IPROACTV(12+JSV,ILAST_PROC_NBR) = 1
    CASE (2)                               ! Concentration of cloud droplets
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR)= 'HENU_'
      IPROACTV(12+JSV,ILAST_PROC_NBR) = 1
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR)= 'SELF_'
      IPROACTV(12+JSV,ILAST_PROC_NBR) = 1
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR)= 'ACCR_'
      IPROACTV(12+JSV,ILAST_PROC_NBR) = 1
      IF (LSEDC) THEN
        ILAST_PROC_NBR = ILAST_PROC_NBR + 1
        YWORK2(12+JSV,ILAST_PROC_NBR)= 'SEDI_'
        IPROACTV(12+JSV,ILAST_PROC_NBR) = 1
      END IF
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR)= 'CEVA_'
      IPROACTV(12+JSV,ILAST_PROC_NBR) = 1
    CASE (3)                               ! Concentration of raindrops
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR)= 'AUTO_'
      IPROACTV(12+JSV,ILAST_PROC_NBR) = 1
      IF (HCLOUD /= 'KHKO') THEN
       ILAST_PROC_NBR = ILAST_PROC_NBR + 1
       YWORK2(12+JSV,ILAST_PROC_NBR)= 'SCBU_'
       IPROACTV(12+JSV,ILAST_PROC_NBR) = 1
      END IF
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR)= 'CEVA_'
      IPROACTV(12+JSV,ILAST_PROC_NBR) = 1
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR)= 'BRKU_'
      IPROACTV(12+JSV,ILAST_PROC_NBR) = 1
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR)= 'SEDI_'
      IPROACTV(12+JSV,ILAST_PROC_NBR) = 1
    END SELECT
    !
  ELSEIF (JSV >= NSV_ELECBEG .AND. JSV <= NSV_ELECEND) THEN
    SELECT CASE(JSV-NSV_ELECBEG+1)
    CASE(1)  ! volumetric charge of water vapor
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR) = 'DEPS_'
      IPROACTV(12+JSV,ILAST_PROC_NBR) = NDEPSQV
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR) = 'DEPG_'
      IPROACTV(12+JSV,ILAST_PROC_NBR) = NDEPGQV
      IF (LWARM) THEN
        ILAST_PROC_NBR = ILAST_PROC_NBR + 1
        YWORK2(12+JSV,ILAST_PROC_NBR) = 'REVA_'
        IPROACTV(12+JSV,ILAST_PROC_NBR) = NREVAQV
      END IF
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR) = 'DEPI_'
      IPROACTV(12+JSV,ILAST_PROC_NBR) = NDEPIQV
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR) = 'NEUT_'
      IPROACTV(12+JSV,ILAST_PROC_NBR) = NNEUTQV
    CASE(2)  ! volumetric charge of cloud droplets
      IF (LWARM) THEN
        ILAST_PROC_NBR = ILAST_PROC_NBR + 1
        YWORK2(12+JSV,ILAST_PROC_NBR) = 'AUTO_'
        IPROACTV(12+JSV,ILAST_PROC_NBR) = NAUTOQC
        ILAST_PROC_NBR = ILAST_PROC_NBR + 1
        YWORK2(12+JSV,ILAST_PROC_NBR) = 'ACCR_'
        IPROACTV(12+JSV,ILAST_PROC_NBR) = NACCRQC
      END IF
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR) = 'RIM_'
      IPROACTV(12+JSV,ILAST_PROC_NBR) = NRIMQC
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR) = 'WETG_'
      IPROACTV(12+JSV,ILAST_PROC_NBR) = NWETGQC
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR) = 'DRYG_'
      IPROACTV(12+JSV,ILAST_PROC_NBR) = NDRYGQC
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR) = 'IMLT_'
      IPROACTV(12+JSV,ILAST_PROC_NBR) = NIMLTQC
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR) = 'BERFI_'
      IPROACTV(12+JSV,ILAST_PROC_NBR) = NBERFIQC
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR) = 'DEPI_'
      IPROACTV(12+JSV,ILAST_PROC_NBR) = NDEPIQC
      IF (LINDUCTIVE) THEN
        ILAST_PROC_NBR = ILAST_PROC_NBR + 1
        YWORK2(12+JSV,ILAST_PROC_NBR) = 'IND_'
      END IF
        IPROACTV(12+JSV,ILAST_PROC_NBR) = NINDQC
      IF (LSEDIC) THEN
        ILAST_PROC_NBR = ILAST_PROC_NBR + 1
        YWORK2(12+JSV,ILAST_PROC_NBR) = 'SEDI_'
        IPROACTV(12+JSV,ILAST_PROC_NBR) = NSEDIQC
      END IF
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR) = 'NEUT_'
      IPROACTV(12+JSV,ILAST_PROC_NBR) = NNEUTQC
    CASE(3)  ! volumetric charge of rain drops
      IF (LWARM) THEN
        ILAST_PROC_NBR = ILAST_PROC_NBR + 1
        YWORK2(12+JSV,ILAST_PROC_NBR) = 'AUTO_'
        IPROACTV(12+JSV,ILAST_PROC_NBR) = NAUTOQR
        ILAST_PROC_NBR = ILAST_PROC_NBR + 1
        YWORK2(12+JSV,ILAST_PROC_NBR) = 'ACCR_'
        IPROACTV(12+JSV,ILAST_PROC_NBR) = NACCRQR
        ILAST_PROC_NBR = ILAST_PROC_NBR + 1
        YWORK2(12+JSV,ILAST_PROC_NBR) = 'REVA_'
        IPROACTV(12+JSV,ILAST_PROC_NBR) = NREVAQR
      END IF
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR) = 'ACC_'
      IPROACTV(12+JSV,ILAST_PROC_NBR) = NACCQR
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR) = 'CFRZ_'
      IPROACTV(12+JSV,ILAST_PROC_NBR) = NCFRZQR
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR) = 'WETG_'
      IPROACTV(12+JSV,ILAST_PROC_NBR) = NWETGQR
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR) = 'DRYG_'
      IPROACTV(12+JSV,ILAST_PROC_NBR) = NDRYGQR
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR) = 'GMLT_'
      IPROACTV(12+JSV,ILAST_PROC_NBR) = NGMLTQR
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR) = 'SEDI_'
      IPROACTV(12+JSV,ILAST_PROC_NBR) = NSEDIQR
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR) = 'NEUT_'
      IPROACTV(12+JSV,ILAST_PROC_NBR) = NNEUTQR
    CASE(4)  ! volumetric charge of ice crystals
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR) = 'AGGS_'
      IPROACTV(12+JSV,ILAST_PROC_NBR) = NAGGSQI
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR) = 'AUTS_'
      IPROACTV(12+JSV,ILAST_PROC_NBR) = NAUTSQI
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR) = 'CFRZ_'
      IPROACTV(12+JSV,ILAST_PROC_NBR) = NCFRZQI
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR) = 'WETG_'
      IPROACTV(12+JSV,ILAST_PROC_NBR) = NWETGQI
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR) = 'DRYG_'
      IPROACTV(12+JSV,ILAST_PROC_NBR) = NDRYGQI
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR) = 'IMLT_'
      IPROACTV(12+JSV,ILAST_PROC_NBR) = NIMLTQI
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR) = 'BERFI_'
      IPROACTV(12+JSV,ILAST_PROC_NBR) = NBERFIQI
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR) = 'DEPI_'
      IPROACTV(12+JSV,ILAST_PROC_NBR) = NDEPIQI
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR) = 'NIIS_'
      IPROACTV(12+JSV,ILAST_PROC_NBR) = NNIISQI
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR) = 'SEDI_'
      IPROACTV(12+JSV,ILAST_PROC_NBR) = NSEDIQI
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR) = 'NEUT_'
      IPROACTV(12+JSV,ILAST_PROC_NBR) = NNEUTQI
    CASE(5)  ! volumetric charge of snow
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR) = 'DEPS_'
      IPROACTV(12+JSV,ILAST_PROC_NBR) = NDEPSQS
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR) = 'AGGS_'
      IPROACTV(12+JSV,ILAST_PROC_NBR) = NAGGSQS
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR) = 'AUTS_'
      IPROACTV(12+JSV,ILAST_PROC_NBR) = NAUTSQS
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR) = 'RIM_'
      IPROACTV(12+JSV,ILAST_PROC_NBR) = NRIMQS
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR) = 'ACC_'
      IPROACTV(12+JSV,ILAST_PROC_NBR) = NACCQS
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR) = 'CMEL_'
      IPROACTV(12+JSV,ILAST_PROC_NBR) = NCMELQS
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR) = 'WETG_'
      IPROACTV(12+JSV,ILAST_PROC_NBR) = NWETGQS
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR) = 'DRYG_'
      IPROACTV(12+JSV,ILAST_PROC_NBR) = NDRYGQS
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR) = 'NIIS_'
      IPROACTV(12+JSV,ILAST_PROC_NBR) = NNIISQS
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR) = 'SEDI_'
      IPROACTV(12+JSV,ILAST_PROC_NBR) = NSEDIQS
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR) = 'NEUT_'
    CASE(6)  ! volumetric charge of graupel
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR) = 'DEPG_'
      IPROACTV(12+JSV,ILAST_PROC_NBR) = NDEPGQG
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR) = 'RIM_'
      IPROACTV(12+JSV,ILAST_PROC_NBR) = NRIMQG
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR) = 'ACC_'
      IPROACTV(12+JSV,ILAST_PROC_NBR) = NACCQG
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR) = 'CMEL_'
      IPROACTV(12+JSV,ILAST_PROC_NBR) = NCMELQG
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR) = 'CFRZ_'
      IPROACTV(12+JSV,ILAST_PROC_NBR) = NCFRZQG
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR) = 'WETG_'
      IPROACTV(12+JSV,ILAST_PROC_NBR) = NWETGQG
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR) = 'DRYG_'
      IPROACTV(12+JSV,ILAST_PROC_NBR) = NDRYGQG
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR) = 'GMLT_'
      IPROACTV(12+JSV,ILAST_PROC_NBR) = NGMLTQG
      IF (LINDUCTIVE) THEN
        ILAST_PROC_NBR = ILAST_PROC_NBR + 1
        YWORK2(12+JSV,ILAST_PROC_NBR) = 'IND_'
        IPROACTV(12+JSV,ILAST_PROC_NBR) = NINDQG
      END IF
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR) = 'SEDI_'
      IPROACTV(12+JSV,ILAST_PROC_NBR) = NSEDIQG
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR) = 'NEUT_'
      IPROACTV(12+JSV,ILAST_PROC_NBR) = NNEUTQG
    CASE(7)  ! volumetric charge of hail
! add budget for hail volumetric charge
    END SELECT
!
  ELSEIF (JSV >= NSV_CHEMBEG .AND. JSV <= NSV_CHEMEND) THEN
    ! Chemical Case
      ILAST_PROC_NBR = ILAST_PROC_NBR + 1
      YWORK2(12+JSV,ILAST_PROC_NBR)= 'CHEM_'
      IPROACTV(12+JSV,ILAST_PROC_NBR) = NCHEMSV

  END IF
  !
  IF (LHOOK) CALL DR_HOOK('INI_BUDGET:BUDGET_OTHERPROC_SV',1,ZHOOK_HANDLE)
  END SUBROUTINE BUDGET_OTHERPROC_SV
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE INI_BUDGET
