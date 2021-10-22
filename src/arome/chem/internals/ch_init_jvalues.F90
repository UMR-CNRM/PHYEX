!     ######spl
       SUBROUTINE CH_INIT_JVALUES(KDAY, KMONTH, KYEAR, KLUOUT)
       USE PARKIND1, ONLY : JPRB
       USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!      ###################################
!!
!!***  *CH_INIT_JVALUES*
!!
!!    PURPOSE
!!    -------
!      calculate J-Values used by MesoNH-chemistry
!!
!!**  METHOD
!!    ------
!!    The radiative transfer model TUV is called for a dicrete number
!!    of solar zenith angle (latitude, longitude, time) and altitude.
!!    TUV itself requires
!!    a number of input files that will reside in directories
!!    DATAX, DATA0 and DATA4 (to be created at runtime by tar xvf TUVDATA.tar).
!!    Temporal interpolation is performed in CH_MONITOR.
!!
!!==============================================================================
!!
!!    REFERENCES
!!    ----------
!!    MesoNH-chemistry book 3
!!
!!    AUTHOR
!!    ------
!!    C. Mari
!!
!!    MODIFICATIONS
!!    -------------
!!    Original 11/03/01
!!    Modification
!!          12/03/01 (P. Tulet) Include variation of J with surface albedo
!!          01/02/04 (P. Tulet) externalisation , modification of the albedo UV
!!                              interpolation
!!          01/12/04 (P. Tulet) update for arome
!!
!!    EXTERNAL
!!    --------
!!    TUV39.f (Fortran 77 code from S. Madronich)
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!    None
!
USE MODD_CH_INIT_JVALUES
!
!*       0.   DECLARATIONS
!        -----------------
!
IMPLICIT NONE
!
INTEGER, INTENT(IN)    :: KLUOUT
INTEGER, INTENT(IN)    :: KYEAR       ! Current Year
INTEGER, INTENT(IN)    :: KMONTH      ! Current Month
INTEGER, INTENT(IN)    :: KDAY        ! Current Day

!
! Calculation of J values
!
REAL, DIMENSION(NZZ_JVAL)         :: ZLWC
REAL                              :: ZDZ
REAL                              :: ZDOBNEW
REAL                              :: ZSZALOOP, ZALBLOOP
REAL                              :: ZMAX = 30E3
                ! MAXIMUM HEIGHT FOR WHICH J-VALUES WILL BE COMPUTED
INTEGER                           :: IDATE
!
! J value storage
!
CHARACTER*40, DIMENSION(JPJVMAX)  :: YLABELOUT
REAL, DIMENSION(NZZ_JVAL,JPJVMAX) :: ZJOUT
!
! Parameters for interpolation
!
INTEGER                           :: JALB, JKLEV, JSZA, JJVAL
!
!----------------------------------------------------------------------------
!
!
!*       1.   INITIALISATION
!        -------------------
!
! Nb of interpolation values for albedo; values = 10
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('CH_INIT_JVALUES',0,ZHOOK_HANDLE)
NBALB=10
IF (.NOT.ALLOCATED(XJDATA)) ALLOCATE(XJDATA(NSZA_INCR,NZZ_JVAL,JPJVMAX,NBALB))
  !
  !* Albedo
  ! now denifed by radiation scheme
  !
  !* Solar Zenith Angle
  !
  IF (.NOT. ALLOCATED(XSZA_JVAL)) ALLOCATE(XSZA_JVAL(NSZA_INCR))
  DO JSZA = 1, NSZA_INCR 
   XSZA_JVAL(JSZA) = FLOAT(JSZA-1)
  ENDDO
  !
  !* Ozone Column
  !
  ZDOBNEW = 320.
  !
  !* Date
  !
  IDATE= MOD(KYEAR,100) *10000 + KMONTH * 100 + KDAY
  !
  !* Vertical Levels
  !
  ZDZ = ZMAX / FLOAT(NZZ_JVAL - 1)
  IF(.NOT.ALLOCATED(XZZ_JVAL))  ALLOCATE(XZZ_JVAL(NZZ_JVAL))
  DO JKLEV = 1, NZZ_JVAL
    XZZ_JVAL(JKLEV) = FLOAT(JKLEV-1) * ZDZ
    ZLWC(JKLEV)= 0.0
  ENDDO
  !
  !*       CALL TUV 3.9
  !        ------------
  !
  DO JALB=1,NBALB
    ZALBLOOP=0.02+0.20*FLOAT(JALB-1)/FLOAT(NBALB-1) 
    DO JSZA = 1, NSZA_INCR
      ZSZALOOP = XSZA_JVAL(JSZA)
      CALL TUVMAIN(  ZSZALOOP, IDATE,           &
                     ZALBLOOP, ZDOBNEW,         &
                     NZZ_JVAL,XZZ_JVAL, ZLWC,   &
                     JPJVMAX, ZJOUT, YLABELOUT  )
       DO JKLEV = 1, NZZ_JVAL
         DO JJVAL = 1, JPJVMAX
             XJDATA(JSZA,JKLEV,JJVAL,JALB) = ZJOUT(JKLEV,JJVAL)
         ENDDO
       ENDDO
    ENDDO
  ENDDO
  XJDATA(:,:,:,:)=MAX(0.,XJDATA(:,:,:,:))
  !
WRITE(KLUOUT,*) "There are ", NZZ_JVAL, " vertical levels with a spacing of ",&
                 ZDZ, "m"
WRITE(KLUOUT,*) "The solar zenith angle varies from", XSZA_JVAL(1)," to ", &
                 XSZA_JVAL(NSZA_INCR)
WRITE(KLUOUT,*) "There are the following ", JPJVMAX, " photolysis reactions:"
DO JJVAL = 1, JPJVMAX
    WRITE(KLUOUT,'(I4, A,  A70)') JJVAL, ". ", YLABELOUT(JJVAL)
ENDDO
DO JJVAL = 1, JPJVMAX
 WRITE(KLUOUT,'(A12,1X,I2,E12.4,1X,E12.4)')'J Values',JJVAL, &
                MAXVAL(XJDATA(:,:,JJVAL,:)), MINVAL(XJDATA(:,:,JJVAL,:))
ENDDO
!
IF (LHOOK) CALL DR_HOOK('CH_INIT_JVALUES',1,ZHOOK_HANDLE)
END SUBROUTINE CH_INIT_JVALUES
