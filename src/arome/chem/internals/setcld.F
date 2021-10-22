!     ######spl
      SUBROUTINE setcld(nz,z,nw,wl,                                     &
     &                  lwc, nlevel,                                    &
     &                  dtcld,omcld,gcld)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!XXX      SUBROUTINE setcld(nz,z,nw,wl,dtcld,omcld,gcld)
!KS:11===================================================================


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
!KS:12===================================================================
!XXX      PARAMETER(kdata=51)
      PARAMETER(kdata=151)
!KS:12===================================================================


! input: (grids)
      REAL wl(kw)
      REAL z(kz)
      INTEGER nz
      INTEGER nw
!KS:13===================================================================
!     LWC is the liquid water content (!! kg/m3 !!) on the calling model
!     grid (which has NLEVEL points: Z(1:NLEVEL) = AZ(*)
      REAL lwc(*)
      INTEGER nlevel
!KS:13===================================================================


! Output:
      REAL dtcld(kz,kw), omcld(kz,kw), gcld(kz,kw)

! local:

! specified data:
      REAL zd(kdata), cd(kdata), omd(kdata), gd(kdata)
      REAL womd(kdata), wgd(kdata)

! other:
      REAL cz(kz)
      REAL omz(kz)
      REAL gz(kz)
      INTEGER i, iw, n

! External functions:
      REAL fsum
      EXTERNAL fsum
!_______________________________________________________________________


! cloud properties are set for each layer (not each level)

! Set as many clouds as want here:
! First choose a cloud grid, zd(n), in km above sea level
! Can allow altitude variation of omega, g:

!KS:14===================================================================
!XXX      n = 4
!XXX
!XXX      zd(1) = 5.
!XXX      cd(1) = 0.
!XXX      omd(1) = .9999
!XXX      gd(1) = .85
!XXX
!XXX      zd(2) = 7.
!XXX      cd(2) = 0.
!XXX      omd(2) = .5
!XXX      gd(2) = .5
!XXX
!XXX      zd(3) = 9.
!XXX      cd(3) = 0.
!XXX      omd(3) = .9999
!XXX      gd(3) = .85
!XXX
!XXX      zd(4) = 11.
!XXX
          REAL(KIND=JPRB) :: ZHOOK_HANDLE
          IF (LHOOK) CALL DR_HOOK('SETCLD',0,ZHOOK_HANDLE)
          n = nlevel + 1
          if (n .gt. kdata) stop "SETCLD: not enough memory: KDATA"
          zd(1) = 0
          do 110, i = 2, n
            zd(i) = 0.5*( z(i-1) + z(i) )
110       continue
!
!         calculate cloud optical properties
          do 120, i = 1, nlevel
!
!           reference: Fouquart et al., Rev. Geophys., 1990
!           TAU = 3/2 LWC*DZ / (RHOWATER * Reff)
!           RHOWATER = 1E3 kg/m3
!           Reff = (11 w + 4) 1E-6
!           w = LWC * 1E+3 (in g/cm3, since LWC is given in kg/m3)
!
            cd(i)  = 1.5 * ( lwc(i) * 1E3*(zd(i+1) - zd(i)) )           &
     &             / ( 1E3 * (11.*lwc(i)*1E+3+4.) * 1E-6)
            omd(i) = .9999
            gd(i)  = .85
!           print '(A,I5,99E12.5)', "I,TAU,LWC,REFF(um)"
!    +            , i, cd(i), lwc(i)
!    +            , ((11.*lwc(i)*1E+3+4.))
120       continue
!KS:14===================================================================

!*****************

! compute integrals and averages over grid layers:
! for g and omega, use averages weigthed by optical depth

!     DO 11, i = 1, n    !***** CHANGED!!See header!!*****
      DO 11, i = 1, n-1
         womd(i) = omd(i) * cd(i)
         wgd(i) = gd(i) * cd(i)
   11 CONTINUE
      CALL inter3(nz,z,cz,  n, zd,cd, 0)
      CALL inter3(nz,z,omz, n, zd,womd, 0)
      CALL inter3(nz,z,gz , n, zd,wgd, 0)

      DO 15, i = 1, nz-1
         IF (cz(i) .GT. 0.) THEN
            omz(i) = omz(i)/cz(i)
            gz(i)  = gz(i) /cz(i)
         ELSE
            omz(i) = 1.
            gz(i) = 0.
         ENDIF
   15 CONTINUE

!      WRITE(kout,*) 'Cloud: ', n, 'levels, tot opt. dep. = ',
!     $     fsum(nz-1,cz)

! assign at all wavelengths
! (can move wavelength loop outside if want to vary with wavelength)

      DO 17, iw = 1, nw-1
         DO 16, i = 1, nz-1
            dtcld(i,iw) = cz(i)
            omcld(i,iw) = omz(i)
            gcld (i,iw) = gz(i)
   16    CONTINUE
   17 CONTINUE
!_______________________________________________________________________

      IF (LHOOK) CALL DR_HOOK('SETCLD',1,ZHOOK_HANDLE)
      RETURN
      ENDSUBROUTINE SETCLD
