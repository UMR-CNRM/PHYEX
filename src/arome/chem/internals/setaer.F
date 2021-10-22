!     ######spl
      SUBROUTINE setaer(nz,z,nw,wl,dtaer,omaer,gaer)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK


      IMPLICIT NONE
!_________________________________________________
!ks:   INCLUDE 'params' has been inserted here
!_________________________________________________
! BROADLY USED PARAMETERS:
!_________________________________________________
! i/o file unit numbers
      INTEGER kout, kin
! output
      PARAMETER(kout=6)
! input
      PARAMETER(kin=78)
!_________________________________________________
! altitude, wavelength grids
      INTEGER kz, kw
! altitude
      PARAMETER(kz=151)
! wavelength
      PARAMETER(kw=150)
!_________________________________________________
! number of weighting functions
      INTEGER ks, kj
!  wavelength dependent
      PARAMETER(ks=40)
!  wavelength and altitude dependent
      PARAMETER(kj=50)

! delta for adding points at beginning or end of data grids
      REAL deltax
      PARAMETER (deltax = 1.E-4)

! some constants...

! pi:
      REAL pi
      PARAMETER(pi=3.1415926535898)

! radius of the earth:
      REAL radius
      PARAMETER(radius=6.371E+3)

! largest number of the machine:
      REAL largest
      PARAMETER(largest=1.E+36)
!_________________________________________________

      INTEGER kdata
      PARAMETER(kdata=51)

! input: (grids)
      REAL wl(kw)
      REAL z(kz)
      INTEGER nz
      INTEGER nw

! output: (on converted grid)
      REAL dtaer(kz,kw), omaer(kz,kw), gaer(kz,kw)

! local:
      REAL zd(kdata), aer(kdata)
      REAL cd(kdata), omd(kdata), gd(kdata)
      REAL womd(kdata), wgd(kdata)

      REAL cz(kz)
      REAL omz(kz)
      REAL gz(kz)

      REAL wc, wscale
      INTEGER i, iw, nd

      LOGICAL aerosl

      REAL fsum
      EXTERNAL fsum

!_______________________________________________________________________

! Aerosol data from Elterman (1968)
! These are verical optical depths per km, in 1 km
! intervals from 0 km to 50 km, at 340 nm.
! This is one option.  User can specify different data set.

      DATA aer/                                                         &
     &     2.40E-01,1.06E-01,4.56E-02,1.91E-02,1.01E-02,7.63E-03,       &
     &     5.38E-03,5.00E-03,5.15E-03,4.94E-03,4.82E-03,4.51E-03,       &
     &     4.74E-03,4.37E-03,4.28E-03,4.03E-03,3.83E-03,3.78E-03,       &
     &     3.88E-03,3.08E-03,2.26E-03,1.64E-03,1.23E-03,9.45E-04,       &
     &     7.49E-04,6.30E-04,5.50E-04,4.21E-04,3.22E-04,2.48E-04,       &
     &     1.90E-04,1.45E-04,1.11E-04,8.51E-05,6.52E-05,5.00E-05,       &
     &     3.83E-05,2.93E-05,2.25E-05,1.72E-05,1.32E-05,1.01E-05,       &
     &     7.72E-06,5.91E-06,4.53E-06,3.46E-06,2.66E-06,2.04E-06,       &
     &     1.56E-06,1.19E-06,9.14E-07/
!_______________________________________________________________________

! initialize

      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('SETAER',0,ZHOOK_HANDLE)
      DO 15, iw = 1, nw - 1
         DO 10, i = 1, nz - 1
            dtaer(i,iw) = 0.
            omaer(i,iw) = 1.
            gaer(i,iw) = 0.
   10    CONTINUE
   15 CONTINUE

! if dont want any aerosols, set AEROSL = .FALSE.

!     aerosl = .FALSE.
      aerosl = .TRUE.
      IF (.NOT. aerosl) THEN
         WRITE(kout,*) 'no aerosols'
         IF (LHOOK) CALL DR_HOOK('SETAER',1,ZHOOK_HANDLE)
         RETURN
      ENDIF

! Altitudes corresponding to Elterman profile, from bottom to top:

      WRITE(kout,*)'aerosols:  Elterman (1968)'
      nd = 51
      DO 22, i = 1, nd
         zd(i) = FLOAT(i-1)
   22 CONTINUE

! assume these are point values (at each level), so find column
! increments

      DO 27, i = 1, nd - 1
         cd(i) = (aer(i+1) + aer(i)) / 2.
         omd(i) = .99
         gd(i) = .61
   27 CONTINUE

!********** end data input.

! Compute integrals and averages over grid layers:
! for g and omega, use averages weigthed by optical depth

      DO 29, i = 1, nd-1
         womd(i) = omd(i) * cd(i)
         wgd(i) = gd(i) * cd(i)
   29 CONTINUE
      CALL inter3(nz,z,cz, nd,zd,cd, 1)
      CALL inter3(nz,z,omz, nd, zd,womd, 1)
      CALL inter3(nz,z,gz , nd, zd,wgd, 1)
      DO 30, i = 1, nz-1
         IF (cz(i) .GT. 0.) THEN
            omz(i) = omz(i)/cz(i)
            gz(i)  = gz(i) /cz(i)
         ELSE
            omz(i) = 1.
            gz(i) = 0.
         ENDIF
   30 CONTINUE

      WRITE(kout,*) '  Total aerosol od at 340 nm = ', fsum(nz-1,cz)

! assign at all wavelengths
! (can move wavelength loop outside if want to vary with wavelength)

      DO 50, iw = 1, nw - 1
         wc = (wl(iw)+wl(iw+1))/2.

! Elterman's data are for 340 nm, so assume optical depth scales
! inversely with first power of wavelength.

         wscale = 340./wc

! optical depths:

         DO 40, i = 1, nz - 1
            dtaer(i,iw) = cz(i)  * wscale
            omaer(i,iw) = omz(i)
            gaer(i,iw) = gz(i)
   40    CONTINUE

   50 CONTINUE
!_______________________________________________________________________

      IF (LHOOK) CALL DR_HOOK('SETAER',1,ZHOOK_HANDLE)
      RETURN
      ENDSUBROUTINE SETAER
