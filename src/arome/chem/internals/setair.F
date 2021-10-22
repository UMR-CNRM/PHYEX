!     ######spl
      SUBROUTINE setair(pmbnew,                                         &
     &     nz,z,nw,wl,                                                  &
     &     airlev,dtrl,cz)
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
      PARAMETER(kdata=150)

! input: (grids)
      REAL wl(kw)
      REAL z(kz)
      INTEGER nw
      INTEGER nz
      REAL  pmbnew

! output:
! air density (molec cm-3) at each grid level
! Rayleigh optical depths

      REAL airlev(kz)
      REAL dtrl(kz,kw)

! local:
      REAL scale
      real airnew(kdata)
      REAL colold, colnew, pmbold
      REAL pconv
      PARAMETER(pconv = 980.665 * 1.E-3 * 28.9644 / 6.022169E23)
! specified data:
      REAL zd(kdata), air(kdata)
      REAL hscale
      REAL cd(kdata)

! other:
      REAL cz(kz)
      REAL srayl(kw)
      REAL deltaz
      REAL colz, pressz
      REAL wc, wmicrn, xx
      INTEGER i, iw, nd


! External functions:
      REAL fsum
      EXTERNAL fsum

!_______________________________________________________________________

! read in air density profile

      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('SETAIR',0,ZHOOK_HANDLE)
      WRITE(kout,*) 'air density: USSA, 1976'

      OPEN(kin,FILE='DATAE1/ATM/ussa.dens',STATUS='old')
      DO i = 1, 3
         READ(kin,*)
      ENDDO
      nd = 1
 4    CONTINUE
        READ(kin,*,END=5) zd(nd), air(nd)
        nd = nd+1
        GOTO 4
 5    CONTINUE
      CLOSE(kin)
      nd = nd-1

! compute column increments (logarithmic integrals)

      DO 6, i = 1, nd - 1
         deltaz = 1.E5 * (zd(i+1)-zd(i))
         cd(i) =  (air(i+1)-air(i)) /LOG(air(i+1)/air(i)) * deltaz
!         cd(i) = (air(i+1)+air(i)) * deltaz / 2.
    6 CONTINUE

! Include exponential tail integral from infinity to 50 km,
! fold tail integral into top layer
! specify scale height near top of data.

      hscale = 8.05e5
      cd(nd-1) = cd(nd-1) + hscale * air(nd)

! alternative input air density data could include, e.g., a read file here:

! If want, can rescale to any total pressure:

      colold = fsum(nd-1,cd)
      pmbold = colold * pconv
      WRITE(kout,100) colold, pmbold
  100 FORMAT(5x,'old sea level air column = ', 1pe11.4,1x,'# cm-2  = ', &
     &     0pf8.2,' mbar')

! assign new sea level pressure

      if (pmbnew .lt. 0.) then
         scale = 1.
      else
         scale = pmbnew/pmbold
      endif

      DO i = 1, nd-1
         cd(i) = cd(i) * scale
         airnew(i) = air(i) * scale
      ENDDO
      airnew(nd) = air(nd) * scale

      colnew = fsum(nd-1,cd)
      WRITE(kout,105) colnew, colnew * pconv
  105 FORMAT(5x,'new sea level air column = ', 1pe11.4,1x,'# cm-2  = ', &
     &     0pf8.2,' mbar')

!********************* end data input.

! Compute air density at each level

      CALL ch_inter1(nz,z,airlev,nd,zd,airnew)

! Compute column increments on standard z-grid.

      CALL inter3(nz,z,cz, nd,zd,cd, 1)

      colz = fsum(nz-1,cz)
      pressz =  colz * pconv
      write(kout,110) colz, pressz
 110  FORMAT(5x,'surface air column = ', 1pe11.4,1x,'# cm-2  = ',       &
     &     0pf8.2,' mbar')


! compute Rayleigh cross sections and depths:

      DO 30, iw = 1, nw - 1
         wc = (wl(iw) + wl(iw+1))/2.

! Rayleigh scattering cross section
! from WMO 1985 (originally from Nicolet)

         wmicrn =  wc/1.E3
         IF( wmicrn .LE. 0.55) THEN
            xx = 3.677 + 0.389*wmicrn + 0.09426/wmicrn
         ELSE
            xx = 4. + 0.04
         ENDIF
         srayl(iw) = 4.02e-28/(wmicrn)**xx

! alternate (older) expression from
! Frohlich and Shaw, Appl.Opt. v.11, p.1773 (1980).
!     xx = 3.916 + 0.074*wmicrn + 0.050/wmicrn
!     srayl(iw) = 3.90e-28/(wmicrn)**xx

         DO 40, i = 1, nz - 1
            dtrl(i,iw) = cz(i)*srayl(iw)
   40    CONTINUE

   30 CONTINUE
!_______________________________________________________________________

      IF (LHOOK) CALL DR_HOOK('SETAIR',1,ZHOOK_HANDLE)
      RETURN
      ENDSUBROUTINE SETAIR
