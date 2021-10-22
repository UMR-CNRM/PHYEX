!     ######spl
      SUBROUTINE CH_UPDATE_JVALUES_n(KLUOUT, PZENITH, PRT, &
           PALB_UV, PZS, PZZ, PLAT0, PLON0,      &
           KLON, KLAT, KLEV, KRR,                &
           KDAY, KMONTH, KYEAR, PTIME,           &
           OCH_TUV_ONLINE,  HCH_TUV_CLOUDS,      &
           PALBNEW, PDOBNEW, PRHODREF, PJVALUES, &
           NIB,NIE,NJB,NJE,NIU,NJU, KVERB )
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!!    #######################################################
!!
!!*** *CH_UPDATE_JVALUES_n* updates photolysis rates
!!
!!    PURPOSE
!!    -------
!!      Update photolysis rates in MODD_CH_JVALUES_n for later extraction
!!    by CH_JVALUES
!!
!!**  METHOD
!!    ------
!!    It converts some input information into the format expected by the
!!    FORTRAN90 subroutine tuvmain, which is then called in order to
!!    calculate the photolysis rates.
!!
!!
!!    REFERENCE
!!    ---------
!!    MesoNH documentation
!!
!!    AUTHOR
!!    ------
!!    Karsten Suhre (LA)
!!
!!    MODIFICATIONS
!!    -------------
!!    Original 05/03/97
!!    05/03/05  P. Tulet (CNRM/GMEI) Update for Arome
!!
!!------------------------------------------------------------------------------
!!
!!    EXTERNAL
!!    --------
USE MODI_CH_INTERP_JVALUES_n
USE MODI_CH_JVALUES_CLOUDS_n
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
USE MODD_CST
USE MODD_PARAMETERS
USE MODD_CH_INIT_JVALUES, ONLY : JPJVMAX
USE MODD_CONF
!!
!!    EXPLICIT ARGUMENTS
!!    ------------------
IMPLICIT NONE

INTEGER,                  INTENT(IN)   :: KLON     !NPROMA under CPG
INTEGER,                  INTENT(IN)   :: KLAT     !NPROMA under CPG
INTEGER,                  INTENT(IN)   :: KLEV     !Number of vertical levels
INTEGER,                  INTENT(IN)   :: KRR      ! Number of moist variables
REAL,                     INTENT(IN)   :: PLAT0, PLON0
REAL, DIMENSION(KLON,KLAT),      INTENT(IN) :: PZENITH, PALB_UV, PZS
REAL, DIMENSION(KLON,KLAT,KLEV), INTENT(IN) :: PZZ
REAL, DIMENSION(KLON,KLAT,KLEV), INTENT(IN) :: PRHODREF
REAL, DIMENSION(KLON,KLAT,KLEV,KRR),  INTENT(IN) :: PRT
INTEGER,                   INTENT(IN) :: KDAY, KMONTH, KYEAR ! current date
REAL,                      INTENT(IN) :: PTIME    ! current time (s)
INTEGER,                   INTENT(IN) :: KLUOUT
LOGICAL,                   INTENT(IN) :: OCH_TUV_ONLINE ! online/lookup table
CHARACTER*4,               INTENT(IN) :: HCH_TUV_CLOUDS ! clouds and radiation
REAL,                      INTENT(IN) :: PALBNEW  ! surface albedo
REAL,                      INTENT(IN) :: PDOBNEW  ! ozone column dobson
REAL,DIMENSION(KLON,KLAT,KLEV,JPJVMAX), INTENT(INOUT) :: PJVALUES    ! Tuv coefficients
INTEGER,                   INTENT(IN)    :: KVERB      ! verbosity level
INTEGER,                   INTENT(IN)    :: NIB,NIE,NJB,NJE,NIU,NJU   !  domain dim
!!
!!    LOCAL VARIABLES
!!    ---------------
INTEGER       :: JVLEVEL       ! loop counter
REAL, DIMENSION(:,:), ALLOCATABLE :: ZSZA
LOGICAL, SAVE :: GSFIRSTCALL = .TRUE. ! flag for initialization on first call
INTEGER       :: IDATE                ! date in format yymmdd
INTEGER       :: JJ                   ! loop counter
!
REAL,DIMENSION(:,:),     ALLOCATABLE  :: ZJOUT1D ! dummy parameter
!REAL,DIMENSION(:,:,:,:), ALLOCATABLE  :: ZJOUT3D ! dummy parameter
REAL,DIMENSION(:,:,:,:), ALLOCATABLE  :: ZFCLD   ! cloud correction
CHARACTER*40, DIMENSION(JPJVMAX)      :: YJLABELOUT ! names of J-reacts.
!!
REAL, DIMENSION(:), ALLOCATABLE       :: ZAZ
REAL, DIMENSION(:,:,:), ALLOCATABLE   :: ZAZ3D, ZLWC3D
INTEGER                               :: IIU,IJU,IKU,IKB,IKE, JPV
REAL                                  :: ZALBNEW, ZDOBNEW
REAL, DIMENSION(:), ALLOCATABLE       :: ZLWC


!
!------------------------------------------------------------------------------
!!
!!    EXECUTABLE STATEMENTS
!!    ---------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('CH_UPDATE_JVALUES_N',0,ZHOOK_HANDLE)
IF (CPROGRAM == "MOD0D") THEN

IIU=1
IJU=1
IKB=1
IKE=1
IKU=1
JPV=0

ELSE

IIU=NIU
IJU=NJU
IKU=SIZE(PZZ,3)
IKB=1+JPVEXT
IKE=IKU-JPVEXT
JPV=JPVEXT

END IF
!
IF (.NOT.ALLOCATED(ZSZA))     ALLOCATE(ZSZA(IIU,IJU))
!
IF (OCH_TUV_ONLINE) THEN
!
   IF ((.NOT.L1D).OR.(CPROGRAM .EQ. "AROME")) THEN
     WRITE(KLUOUT,*)"ERROR in CH_UPDATE_JVALUES_n:                           "
     WRITE(KLUOUT,*)"you want to use ON-LINE calculation of photolysis rates "
     WRITE(KLUOUT,*)"but you are not runnning in 1D                          "
     WRITE(KLUOUT,*)"Program is STOPPED now                                  "
     STOP
  ENDIF

!*        1. TUV 3D ON LINE
!         -------------------
  IDATE = 10000*MOD(KYEAR,100) + 100*KMONTH + KDAY
!
  ZALBNEW          = PALBNEW
  ZDOBNEW          = PDOBNEW
!
   IF (.NOT.ALLOCATED(ZAZ)) THEN
       ALLOCATE(ZAZ(IKE-IKB+1))
  ENDIF
  IF (.NOT.ALLOCATED(ZLWC)) THEN
       ALLOCATE(ZLWC(IKE-IKB+1))
  ENDIF
  ZAZ(1:IKE-IKB+1) = PZZ(1+JPV,1+JPV,1:IKE-IKB+1)
!
  IF (HCH_TUV_CLOUDS .EQ. "NONE") THEN
    ZLWC(1:IKE-IKB+1) = 0.0
  ELSE
    ZLWC(1:IKE-IKB+1) = PRT(1+JPV,1+JPV,1:IKE-IKB+1,2) * PRHODREF(1+JPV,1+JPV,1:IKE-IKB+1)
  END IF
!
  JVLEVEL = IKE-IKB+1
!
  ALLOCATE(ZJOUT1D(JVLEVEL,JPJVMAX))
!
  ZSZA(:,:) = PZENITH(:,:) * 180. / XPI

  CALL TUVMAIN(ZSZA, IDATE, ZALBNEW, ZDOBNEW,  &
               JVLEVEL, ZAZ, ZLWC,             &
               JPJVMAX, ZJOUT1D, YJLABELOUT    )
!
  PJVALUES(:,:,:,:) = SPREAD(SPREAD(ZJOUT1D(1:JVLEVEL,1:JPJVMAX),1,IIU),2,IJU)
  IF (KVERB >= 6) THEN
  WRITE(KLUOUT,*)'PJVALUES',SIZE(PJVALUES,1),SIZE(PJVALUES,2),SIZE(PJVALUES,3),SIZE(PJVALUES,4)
  WRITE(KLUOUT,*)'MAX',MAXVAL(PJVALUES), MAXLOC(PJVALUES)
  WRITE(KLUOUT,*)'MIN',MINVAL(PJVALUES), MINLOC(PJVALUES)
  ENDIF
!
  IF ((GSFIRSTCALL).AND.(KVERB >= 6)) THEN
!
WRITE(KLUOUT,*) "CH_UPDATE_JVALUES_n information on first call"
WRITE(KLUOUT,*) "-------------------------------------------"
WRITE(KLUOUT,*) "You are using ON LINE computation of photolysis rates "
WRITE(KLUOUT,*) "with CH_TUV_CLOUDS ",HCH_TUV_CLOUDS
WRITE(KLUOUT,*) "-------------------------------------------"
WRITE(KLUOUT,*) "J-Values at the bottom, at the center of the domain: Z = " , ZAZ(1)
    DO JJ = 1, JPJVMAX
      WRITE(KLUOUT,*) JJ, YJLABELOUT(JJ), PJVALUES(IIU/2+1,IJU/2+1,1,JJ)
    ENDDO
    WRITE(KLUOUT,*) "J-Values at the top, at the center of the domain: Z = ", ZAZ(JVLEVEL)
    DO JJ = 1, JPJVMAX
      WRITE(KLUOUT,*) JJ, YJLABELOUT(JJ), PJVALUES(IIU/2+1,IJU/2+1,JVLEVEL,JJ)
    ENDDO
    WRITE(KLUOUT,*) "-------------------------------------------"
!
  ENDIF
  IF (ALLOCATED(ZAZ)) THEN
    DEALLOCATE(ZAZ)
  ENDIF
  IF (ALLOCATED(ZLWC)) THEN
    DEALLOCATE(ZLWC)
  ENDIF
!

ELSE
!*        2. TUV 3D OFF LINE
!         -------------------
!
!
  CALL CH_INTERP_JVALUES_n(PJVALUES,PZZ, PZS, PALB_UV, PZENITH,&
                           KLON, KLAT, KLEV, &
                           NIB,NIE,NJB,NJE,NIU,NJU,KVERB, KLUOUT)
  IF (KVERB >= 6) THEN
   WRITE(KLUOUT,*) "-------------------------------------------"
   WRITE(KLUOUT,*) "J-Values at the bottom, at the center of the domain"
    DO JJ = 1, JPJVMAX
      WRITE(KLUOUT,*) JJ, PJVALUES(IIU/2+1,IJU/2+1,IKB,JJ)
    ENDDO
    WRITE(KLUOUT,*) "J-Values at the top, at the center of the domain"
    DO JJ = 1, JPJVMAX
      WRITE(KLUOUT,*) JJ, PJVALUES(IIU/2+1,IJU/2+1,IKE-IKB+1,JJ)
    ENDDO
    WRITE(KLUOUT,*) "-------------------------------------------"
  ENDIF
!
  IF (HCH_TUV_CLOUDS .EQ. "CHAN") THEN
!
    IF (.NOT.ALLOCATED(ZAZ3D))  ALLOCATE(ZAZ3D(1:IIU,1:IJU,1:IKU))
    IF (.NOT.ALLOCATED(ZLWC3D)) ALLOCATE(ZLWC3D(1:IIU,1:IJU,1:IKU))
    IF (.NOT.ALLOCATED(ZFCLD))  ALLOCATE(ZFCLD(1:IIU,1:IJU,1:IKU,JPJVMAX))
!
     ZAZ3D(:,:,:) = PZZ(:,:,:)
     ZLWC3D(:,:,:) = PRT(:,:,:,2) * PRHODREF(:,:,:)
!
    ZFCLD(:,:,:,:) = 1.
!
    CALL CH_JVALUES_CLOUDS_n(PZENITH, JPJVMAX, IKU, ZAZ3D, ZLWC3D, ZFCLD)
!
    PJVALUES(:,:,:,:) = PJVALUES(:,:,:,:) * ZFCLD(:,:,:,:)
    IF (KVERB >= 6) WRITE(KLUOUT,*)" J VALUES CLOUD CORR MAX/MIN", MAXVAL(ZFCLD), MINVAL(ZFCLD)
!
  ENDIF
!
  IF (ALLOCATED(ZAZ3D)) THEN
    DEALLOCATE(ZAZ3D)
  ENDIF
  IF (ALLOCATED(ZLWC3D)) THEN
    DEALLOCATE(ZLWC3D)
  ENDIF
!
END IF
GSFIRSTCALL = .FALSE.
!
IF (LHOOK) CALL DR_HOOK('CH_UPDATE_JVALUES_N',1,ZHOOK_HANDLE)
END SUBROUTINE CH_UPDATE_JVALUES_n
