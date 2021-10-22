!     ######spl
      subroutine tuvmain (asza, idate,                                  &
     &           albnew, dobnew,                                        &
     &           nlevel, az, lwc,                                       &
     &           njout, jout, jlabelout)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!_______________________________________________________________________
!     PROGRAM tuv
!     Tropospheric Ultraviolet-Visible radiation model
!     version 3.9
!     Feb. 97 by Madronich et al.
!_______________________________________________________________________
! DATE:    Mon Mar 15 09:50:41 MET 1999
! FILE:    tuv39.f
! MODIFIED BY: Karsten Suhre as follows
! - all modifications are indicated by *KS: or *ks:
! - concatenate all sources in one file
! - insert all include 'params' directly into the code
!   and change the following parameters: iout=6, kin=78, kz=151, kw=150
! - modify tuv.f to become a subroutine called tuvmain.f
! - save ALL variables (on HP) and only concerned variables (on LINUX pgf90)
!
!_______________________________________________________________________
!-----------------------------------------------------------------------------*
!= This program is free software;  you can redistribute it and/or modify     =*
!= it under the terms of the GNU General Public License as published by the  =*
!= Free Software Foundation;  either version 2 of the license, or (at your   =*
!= option) any later version.                                                =*
!= The TUV package is distributed in the hope that it will be useful, but    =*
!= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
!= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
!= License for more details.                                                 =*
!= To obtain a copy of the GNU General Public License, write to:             =*
!= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
!-----------------------------------------------------------------------------*
!= To contact the authors, please mail to:                                   =*
!= Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
!= send email to:  sasha@ucar.edu                                            =*
!-----------------------------------------------------------------------------*
!= Copyright (C) 1994,95,96  University Corporation for Atmospheric Research =*
!-----------------------------------------------------------------------------*

      IMPLICIT NONE
      SAVE

! include parameter file

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

! ___ SECTION 1: VARIABLES AND PARAMETERS ______________________________
!
! interface variables:
!
! INPUT:
!       ALAT  - latitude in decimal degrees
!       ALONG - longitude in decimal degrees
!       IDATE - Date at Greenwich - specify year (yy), month (mm), day (dd)
!               format is six-digit integer:  yymmdd
!               (only a test if yy is a factor of 4 will be perfomed,
!                so there is no year 2000 problem!)
!       UT    - Greenwich mean time - decimal military eg.
!               22.75 = 45 min after ten pm UT
!       ALBNEW- surface albedo, if negative the default value will be used
!       DOBNEW- scale the O3 profile to the given Dobson unit
!                (no scaling is done if this variable is negative)
!       NLEVEL- number of points in AZ
!       AZ    - vector of altitudes (in m) where J-values are calculated
!               (supposed to be ordered in increasing order!), mass level,
!               dimension = NLEVEL
!       LWC   - liquid water content (kg/m3),
!               dimension = NLEVEL
!
! OUTPUT:
!       NJOUT - number of J-values to be returned (actually 21)
!       JOUT  - array of J-values (dimension must be at least 21)
!               1st dimension = NLEVEL, 2nd dimension >= NJOUT
!       JLABELOUT - array of labels for J-values (dimension must be at least 21)
!               dimension >= NJOUT
!
! method:
!
!     on first call, all necessary data will be read from the files
!     in directories DATA0, DATA4, and DATAX;
!     all variables are saved for future calls to TUV
!
      REAL,         INTENT(IN) :: asza
!      REAL,         INTENT(IN) :: ut
      INTEGER,      INTENT(IN) :: idate
      INTEGER,      INTENT(IN) :: nlevel
      REAL,         INTENT(IN) :: dobnew, albnew
      REAL,         INTENT(IN) :: az(nlevel)
      REAL,         INTENT(IN) :: lwc(nlevel)
      INTEGER,      INTENT(IN) :: njout
      REAL,         INTENT(OUT) :: jout(nlevel,njout)
      CHARACTER*40, INTENT(OUT) :: jlabelout(njout)
!
      LOGICAL LFIRSTCALL
      DATA LFIRSTCALL /.TRUE./
!
!
! geometry

      REAL zen

! altitude grid:

      INTEGER nz, iz
      REAL z(kz)

! slant path lengths in spherical geometry

      INTEGER nid(0:kz)
      REAL dsdh(0:kz,kz)

! wavelength grid:

      INTEGER nw, iw
      REAL wl(kw), wc(kw), wu(kw)

! extra terrestrial solar flux and earth-Sun distance ^-2

      REAL f(kw)
      REAL esrm2

! ozone absorption cross section and ozone optical depth:

      REAL xso3(kw), s226(kw), s263(kw), s298(kw)

! O2 absorption cross section

      REAL xso2(kz,kw)

! SO2 absorption cross section

      REAL xsso2(kw)

! NO2 absorption cross section

      REAL xsno2(kw)

! atmospheric optical parameters:

      REAL tlev(kz), tlay(kz)
      REAL airlev(kz), colinc(kz)
      REAL dtrl(kz,kw)
      REAL dto3(kz,kw), dto2(kz,kw), dtso2(kz,kw), dtno2(kz,kw)
      REAL dtcld(kz,kw), omcld(kz,kw), gcld(kz,kw)
      REAL dtaer(kz,kw), omaer(kz,kw), gaer(kz,kw)
      REAL albedo(kw)

! spectral irradiance and actinic flux (scalar irradiance):

      REAL edir(kz), edn(kz), eup(kz)
      REAL fdir(kz), fdn(kz), fup(kz)

!! j-values:

      INTEGER nj, ij
      REAL sj(kj,kz,kw), valj(kj,kz)
      REAL saflux, deltaj
      CHARACTER*40 jlabel(kj)

! User-defined
! new sea level pressure, surface dobson, etc.

      REAL azim

      REAL pmbnew, so2new, no2new

      integer i, j

!________________________________________________________________________
! BEGIN OF INITIALIZATION (only executed on firstcall of tuvmain)
!________________________________________________________________________

      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('TUVMAIN',0,ZHOOK_HANDLE)
      IF (LFIRSTCALL) THEN

        WRITE(kout,*) 'running TUVMAIN, version 3.9'

! ___ SECTION 2: SET GRIDS _________________________________________________

! wavelengths

      CALL gridw(nw,wl,wc,wu)

! altitudes

      CALL gridz(nz,z,az,nlevel)

! ___ SECTION 3: SPECTRAL DATA ____________________________

! read (and grid) extra terrestrial flux data:

      CALL rdetfl(nw,wl,f)

! read cross section data for
!    ozone (temperature-dependent)
!    SO2
!    NO2

      CALL rdo3xs(nw,wl,xso3,s226,s263,s298)
      CALL rdso2xs(nw,wl,xsso2)
      CALL rdno2xs(nw,wl,xsno2)

! ___ SECTION 4: SET MODEL ATMOSPHERE __________________________________

! temperature profile

      CALL settmp(nz,z,                                                 &
     &     tlev,tlay)

!  air profile and rayleigh optical depths

      pmbnew = -999.
      CALL setair(pmbnew,                                               &
     &     nz,z,nw,wl,                                                  &
     &     airlev,dtrl,colinc)


!!  now can select (altitude-dependent) photolysis reaction spectra
!  (cross section)x(quantum yield)
!  (must know airlev, tlev from setair, settmp)
! reserve last slot for O2 photolysis (extra because of Schuman-Runge band
! parameterization -> zenith angle dependence of effective cross section)

      ij = 0
      CALL jspec1(nw,wl,nz,tlev,airlev,                                 &
     &     ij,sj,jlabel)
      CALL o2spec(ij,jlabel)
      nj = ij

! ozone optical depths (must give temperature)

      CALL setozo(dobnew,                                               &
     &     nz,z,nw,wl,                                                  &
     &     xso3,s226,s263,s298,tlay,                                    &
     &     dto3)

! SO2 optical depth (also has temperature correction)
! so2new = new column SO2, in Dobson Units

      so2new =  0.
      CALL setso2(so2new,                                               &
     &     nz,z,nw,wl,                                                  &
     &     xsso2, tlay,                                                 &
     &     dtso2)

! NO2 optical depth (also has temperature correction)
! no2new = new column NO2, in Dobson Units

      no2new =  0.
      CALL setno2(no2new,                                               &
     &     nz,z,nw,wl,                                                  &
     &     xsno2, tlay,                                                 &
     &     dtno2)

!  aerosol optical depths:

      CALL setaer(nz,z,nw,wl,                                           &
     &     dtaer,omaer,gaer)

        LFIRSTCALL = .FALSE.
      ENDIF

!________________________________________________________________________
! END OF INITIALIZATION (only executed on firstcall of tuvmain)
!________________________________________________________________________

! UPG:CM 08/06/01 surface albedo calculation taken out of firstcall loop
! surface albedo:

      if (albnew .ge. 0.) then
!        WRITE(kout,*)'wavelength-independent albedo = ', albnew
        DO iw = 1, nw - 1
          albedo(iw) = albnew
        enddo
      else
        CALL setalb(nw,wl,albedo)
      endif

!  cloud optical depth (may be modified at each call, so it has ben moved
!  outside the if(firstcall) section:

      CALL setcld(nz,z,nw,wl,                                           &
     &     lwc,nlevel,                                                  &
     &     dtcld,omcld,gcld)

! ___ SECTION 5: TIME AND LOCATION _____________________________________

! date, earth-sun distance

      CALL sundis(idate,esrm2)
!      WRITE(kout,*) 'idate = ', idate,' esrm2 = ', esrm2

! solar zenith angle calculation:

!      WRITE(kout,*)'lat = ',alat,' long = ',along
!
!      CALL zenith(alat,along,idate,ut,azim,zen)
      zen = asza
!      WRITE(kout,*) 'ut = ', ut, 'azimuth = ', azim, ' zen = ', zen
!      WRITE(kout,*) ' zen = ', zen

! ____ SECTION 6: CALCULATE ZENITH-ANGLE DEPENDENT QUANTITIES __________

! slant path lengths for spherical geometry

      CALL sphers(nz, z, zen, dsdh, nid)

!! effective O2 optical depth (SR bands, must know zenith angle!)
! reassign O2 cross section to sj(nj,*,*)

       CALL zero2(dto2,kz,kw)
       CALL zero2(xso2,kz,kw)
       CALL seto2(nz,z,nw,wl,colinc,zen,dto2,xso2)
       CALL sto2xs(nz,nw,xso2,nj,sj)

! ____ SECTION 7: WAVELENGTH LOOP ______________________________________

! initialize for wavelength integration

       call zero2(valj,kj,kz)

      DO 10, iw = 1, nw-1

! monochromatic radiative transfer:

         CALL rtlink(nz,z,                                              &
     &        iw, albedo(iw), zen,                                      &
     &        dsdh,nid,                                                 &
     &        dtrl,                                                     &
     &        dto3,                                                     &
     &        dto2,                                                     &
     &        dtso2,                                                    &
     &        dtno2,                                                    &
     &        dtcld, omcld, gcld,                                       &
     &        dtaer,omaer,gaer,                                         &
     &        edir, edn, eup, fdir, fdn, fup)

!**** spherical irradiance (actinic flux)
! as a function of altitude
! convert to quanta s-1 nm-1 cm-2
! ( 1.e-4 * (wc*1e-9) / (hc = 6.62E-34 * 2.998E8) )

         DO 17 iz = 1, nz
            saflux = f(iw) * esrm2 * 5.039e11 * wc(iw) *                &
     &           (fdir(iz) + fdn(iz) + fup(iz))
            DO 16, ij = 1, nj
               deltaj = saflux * sj(ij,iz,iw)
               valj(ij,iz) = valj(ij,iz) + deltaj * (wu(iw) - wl(iw))
 16         CONTINUE
 17      CONTINUE

 10   CONTINUE

!^^^^^^^^^^^^^^^^ end wavelength loop

! ____ SECTION 8: OUTPUT _______________________________________________
!
!     iz = 1
!     DO 36, ij = 1, nj
!        WRITE(kout,99) ij, jlabel(ij), valj(ij,iz)
!36   CONTINUE
!
! copy labels into output array
! (to be compatible with TUV39, photolysis raction 1-20 will be transfered
!  as 2-21, and the first will be set to zero (former user-added CH3I))

         if (njout .ne. 21) then
           WRITE(kout,*) 'There should be 21 J-Values to be updated!'
           WRITE(kout,*) 'We better stop here ... in tuvmain.f'
           STOP 1
         endif

         jlabelout(1) = 'not used'
         DO ij = 2, njout
           jlabelout(ij) = jlabel(ij-1)
         ENDDO
!
! return J values on AZ grid
!
!        for each individual J value do..
         DO iz = 1, nlevel
           jout(iz,1) = 0.0
         ENDDO
         DO ij = 2, njout
           DO iz = 1, nlevel
             jout(iz,ij) = valj(ij-1,iz)
           ENDDO
         ENDDO
!
!_______________________________________________________________________
   99 FORMAT(I4,1X,A40,1X,1PE10.3)
!_______________________________________________________________________

      IF (LHOOK) CALL DR_HOOK('TUVMAIN',1,ZHOOK_HANDLE)
      ENDSUBROUTINE TUVMAIN
