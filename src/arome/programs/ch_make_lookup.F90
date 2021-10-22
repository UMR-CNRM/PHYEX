!     ######spl
       PROGRAM CH_MAKE_LOOKUP
       USE PARKIND1, ONLY : JPRB
       USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!      ######################
!!
!!***  *CH_MAKE_LOOKUP*
!!
!!    PURPOSE
!!    -------
!      create a lookup table for J-Values used by MesoNH-chemistry
!!
!!**  METHOD
!!    ------
!!      The radiative transfer model TUV is called for a dicrete number
!!    of times. The lookup table is fixed for one latitude, longitude, date, 
!!    surface albedo and ozone column dobson. TUV itself requires
!!    a number of input files that will reside in directories
!!    DATAX, DATA0 and DATA4 (to be created at runtime by tar xvf TUVDATA.tar).
!!    The program may be piloted via a namelistfile LOOKUP1.nam. 
!!    The format of the lookup table is as follows:
!!
!!==============================================================================
!! COMMENT LINE
!! NPHOTO  (NUMBER OF REACTIONS)
!! NLEVEL, DZ  (NUMBER OF VERTICAL POINTS, HEIGHT INCREMENT IN M)
!! NTIME,  DT  (NUMBER OF TEMPORAL POINTS, TIME INCREMENT IN H)
!! REACTION 1
!!  ...
!! REACTION NPHOTO
!! FORTRAN FORMAT FOR THE FOLLOWING DATA
!! HEIGHT1  TIME1  J1 ... JN 
!! HEIGHT1  TIME2  J1 ... JN 
!! ...
!! HEIGHT1  TIMEN  J1 ... JN 
!! HEIGHT2  TIME1  J1 ... JN 
!! ...
!! HEIGHTM  TIMEN  J1 ... JN 
!!==============================================================================
!!
!!    REFERENCES
!!    ----------
!!    MesoNH-chemistry book 3
!!
!!    AUTHOR
!!    ------
!!    K. Suhre
!!
!!    MODIFICATIONS
!!    -------------
!!    Original 01/03/99
!!
!!    EXTERNAL
!!    --------
!!    TUV39.f (Fortran 77 code from S. Madronich)
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!    None
!!
!------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!        -----------------
!
IMPLICIT NONE
!
REAL    :: ALAT, ALONG    ! LATITUDE AND LONGITUDE
REAL    :: ALBNEW, DOBNEW ! SURFACE ALBEDO AND O3 COLUMN DOBSON
INTEGER :: IDATE          ! DATE IN FORMAT YYMMDD
!
INTEGER, PARAMETER :: NTIME  = 96 + 1 ! TEMPORAL DISCRETIZATION
INTEGER, PARAMETER :: NLEVEL = 30 + 1 ! VERTICAL DISCRETIZATION
REAL, DIMENSION(NLEVEL) :: AZ, LWC
REAL, DIMENSION(NTIME)  :: ATIME
REAL :: DZ, DT
REAL :: ZMAX = 30E3       ! MAXIMUM HEIGHT FOR WHICH J-VALUES WILL BE COMPUTED
!
!      J VALUE STORAGE
INTEGER, PARAMETER                    :: NJOUT = 21
REAL, DIMENSION(NLEVEL,NJOUT)         :: JOUT
REAL, DIMENSION(NLEVEL, NJOUT, NTIME) :: JDATA
CHARACTER*40, DIMENSION(NJOUT)        :: JLABELOUT
!
CHARACTER*120 :: HEADDER
REAL          :: UT
INTEGER       :: I, J, K, NJIO
CHARACTER*40  :: YFMT = '(2F11.2,5E11.4/99(7E11.4/))'
!
! NAMELIST for options
NAMELIST /NAM_TUV/ ALAT, ALONG, IDATE, ALBNEW, DOBNEW
!
!
!*       1.   INITIALISATION
!        -------------------
!
!      initialize az and atime
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('CH_MAKE_LOOKUP',0,ZHOOK_HANDLE)
DZ = ZMAX / FLOAT(NLEVEL - 1)
DO J = 1, NLEVEL
  AZ(J) = FLOAT(J-1) * DZ
  LWC(J)= 0.0
ENDDO
DT = 24.00 / FLOAT(NTIME - 1)
DO I = 1, NTIME
  ATIME(I) = FLOAT(I-1) * DT
ENDDO
!
!      initialize default values
ALAT = 45.
ALONG = 0.
IDATE = 970621
ALBNEW = -1.
DOBNEW = -1.
!
!      read the namelist
WRITE(*,*) "TRYING TO OPEN FILE LOOKUP1.nam ..."
OPEN(UNIT=42,FILE="LOOKUP1.nam",STATUS="OLD",FORM="FORMATTED", ERR=100)
READ(42,NAM_TUV)
WRITE(*,*) "NAMELIST NAM_TUV INITIALIZED TO:"
WRITE(*,NAM_TUV)
CLOSE(42)
GOTO 200
!
100    CONTINUE
!      read from input if no namelist file exists
WRITE(*,*) "NO FILE LOOKUP1.NAM FOUND, TRYING INTERACTIVE ..."
PRINT *, "ENTER LATITUDE (+N, -S)"
READ(*,*) ALAT
PRINT *, "ENTER LONGITUDE (+E, -W)"
READ(*,*) ALONG
PRINT *, "ENTER DATE (YYMMDD)"
READ(*,*) IDATE
PRINT *, "ENTER SURFACE ALBEDO (IF <0 -->DATAX/ALBEDO.DAT)"
READ(*,*) ALBNEW
PRINT *, "ENTER TOTAL COLUMN DOBSON (IF <0 NO SCALING)"
READ(*,*) DOBNEW
!
200    CONTINUE
!
WRITE(HEADDER,'(A,F6.1,A,F6.1,A,I6,A,F6.3,A,F8.2)') &
      "TUV39,LAT=", ALAT, ",LON=", ALONG,           &
      ",IDATE(YYMMDD)=", IDATE,                     &
      ",ALB=", ALBNEW,                              &
      ",DOB=", DOBNEW
PRINT *, HEADDER
!
!*       2.   CALL TUV 3.5
!        -----------------
!
DO I = 1, NTIME
  UT = ATIME(I)
  PRINT *, I, "  FROM  ", NTIME, "  DONE ... UT = ", UT
  NJIO = NJOUT
  CALL TUVMAIN(ALAT, ALONG, IDATE, UT, ALBNEW, DOBNEW, NLEVEL, AZ, LWC, &
              NJIO, JOUT, JLABELOUT)
  DO J = 1, NLEVEL
    DO K = 1, NJIO
      JDATA(J, K, I) = JOUT(J,K) 
    ENDDO
  ENDDO
ENDDO
!
!*       3.   WRITE LOOKUP TABLE FILE
!        -----------------
!
OPEN(66,FILE="PHOTO.TUV39",STATUS="UNKNOWN",FORM="FORMATTED")
WRITE(66,'(A)') HEADDER
WRITE(66,'(I6,A)') NJIO,   "  NUMBER OF PHOTOLYSIS REACTIONS"
WRITE(66,'(I6, F10.0, A)') NLEVEL, DZ, &
      "  NUMBER OF LEVELS, VERTICAL INCREMENT (M)"
WRITE(66,'(I6, F10.4, A)') NTIME,  DT, &
      "  NUMBER OF TEMPORAL RECS, TEMPORAL INCREMENT (H)"
WRITE(66,'(A)') (JLABELOUT(K), K=1, NJIO)
WRITE(66,'(A)') YFMT
DO J = 1, NLEVEL
  DO I = 1, NTIME
    WRITE(66,YFMT) AZ(J), ATIME(I), (JDATA(J,K,I), K=1, NJIO)
  ENDDO
ENDDO
CLOSE(66)
!
PRINT *, 'Lookup table file PHOTO.TUV39 has been generated.'
PRINT *, 'CH_MAKE_LOOKUP ended correctly.'
!
IF (LHOOK) CALL DR_HOOK('CH_MAKE_LOOKUP',1,ZHOOK_HANDLE)
END PROGRAM CH_MAKE_LOOKUP
