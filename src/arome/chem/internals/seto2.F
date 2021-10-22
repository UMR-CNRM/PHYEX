!     ######spl
      SUBROUTINE seto2(nz,z,nw,wl,cz,zen,dto2,xso2)
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

      REAL wl(kw)
      REAL z(kz), cz(kz)
      INTEGER nz, nw
      REAL zen
      REAL dto2(kz,kw), xso2(kz,kw)

! grid on which Kockarts' parameterization is defined
      INTEGER ngast
      PARAMETER (ngast = 17)
      REAL wlgast(ngast)
      SAVE wlgast


! O2 optical depth and equivalent cross section on Kockarts' grid
      REAL dto2k(kz,ngast-1), xso2k(kz,ngast-1)


! internal grid and O2 cross section on internal grid
      INTEGER kdata
      PARAMETER (kdata = 200)
      REAL wlint(kdata), xso2int(kdata)
      SAVE wlint, xso2int
      INTEGER nwint
      SAVE nwint

! temporary one-dimensional storage for optical depth and cross section values
! XXtmp  - on internal grid
! XXuser - on user defined grid
      REAL dttmp(2*kw), xstmp(2*kw)
      REAL dtuser(kw), xsuser(kw)

      REAL o2col(kz)

      REAL secchi
      REAL fchap
      EXTERNAL fchap

! cross section data for use outside the SR-Bands (combined from
! Brasseur and Solomon and the JPL 1994 recommendation)
      INTEGER nosr
      PARAMETER (nosr = 105)
      REAL x1(nosr), y1(nosr)

! auxiliaries
      REAL x, y
      REAL dr
      PARAMETER (dr = pi/180.)
      INTEGER i, iw, igast, ierr, icount
      INTEGER iz
      INTEGER ifirst, n

      LOGICAL call1
      SAVE call1, icount
      DATA call1/.TRUE./



! check, whether user grid is in the O2 absorption band at all...
! if not, set cross section and optical depth values to zero and return

      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('SETO2',0,ZHOOK_HANDLE)
      IF (wl(1) .GT. 243.) THEN
         DO iw = 1, nw-1
           DO i = 1, nz
             dto2(i,iw) = 0.
             xso2(i,iw) = 0.
           ENDDO
         ENDDO
         IF (LHOOK) CALL DR_HOOK('SETO2',1,ZHOOK_HANDLE)
         RETURN
      ENDIF

! sec Xhi or Chapman calculation
      IF (zen .LE. 75.) THEN
         secchi = 1./COS(zen*dr)
      ELSEIF (zen .LE. 95. ) THEN
         secchi = fchap(zen)
      ELSE
         IF (LHOOK) CALL DR_HOOK('SETO2',1,ZHOOK_HANDLE)
         RETURN
      ENDIF

! O2 overhead columns calculation
      o2col(nz-1) = 0.2095 * cz(nz-1) * secchi
      DO i = nz-2, 1, -1
        o2col(i) = o2col(i+1) + 0.2095*cz(i)*secchi
      END DO

! read O2 cross section data outside SR-bands only in the very first call
      IF (call1) THEN
!************ O2 absorption cross sections:
! from 116 nm to 245 nm, including Schumann-Runge continumm
! from Brasseur and Solomon 1986.

        OPEN(UNIT=kin,FILE='DATAE1/O2/O2_src.abs')
        READ(kin,*) ifirst, n
        DO i = 1, ifirst-2
          READ(kin,*)
        ENDDO
        icount = 0
        DO i = 1, n
          READ(kin,*) x, y
          IF (x .LT. 204.) THEN
            icount = icount+1
            x1(icount) = x
            y1(icount) = y
          ENDIF
        ENDDO
        CLOSE(kin)

! overwrite from 204 to 241 nm (Herzberg continuum)

        OPEN(UNIT=kin,FILE='DATAE1/O2/O2_jpl94.abs',STATUS='old')
        read(kin,*)
        n = 40
        DO i = 1, n
           icount = icount+1
           READ(kin,*) y
           y1(icount) = y*1E-24
           x1(icount) = 204. + FLOAT(i-1)
        END DO
        CLOSE (kin)

! set values to zero outside the wavelength range defined by the data files

        CALL addpnt(x1,y1,nosr,icount,     x1(1)-deltax,0.)
        CALL addpnt(x1,y1,nosr,icount,               0.,0.)
        CALL addpnt(x1,y1,nosr,icount,x1(icount)+deltax,0.)
        CALL addpnt(x1,y1,nosr,icount,            1.e38,0.)


! set up the internal grid, use full resolution of the cross section data
! outside the SR bands, use Kockarts' grid inside the SR bands

! define Kockarts' grid points
        DO iw = 1, ngast
           wlgast(iw) = 1E7/(57000.-(iw-1)*500.)
        ENDDO

        nwint = 0
        DO iw = 1, 47
           nwint = nwint + 1
           wlint(nwint) = x1(iw)
        ENDDO
        DO iw = 1, ngast
           nwint = nwint+1
           wlint(nwint) = wlgast(iw)
        ENDDO
        DO iw = 65, 105
           nwint = nwint+1
           wlint(nwint) = x1(iw)
        ENDDO

! interpolate Brasseur/Solomon and JPL data onto internal grid
        CALL ch_inter2(nwint,wlint,xso2int, icount,x1,y1, ierr)

        IF (call1) call1 = .FALSE.

      ENDIF

! if necessary:
! do Kockarts' parameterization of the SR bands, output values of O2
! optical depth and O2 equivalent cross section are on his grid
      IF ((wl(1) .LT. wlgast(ngast)) .AND.                              &
     &    (wl(nw) .GT. wlgast(1))) THEN
        DO iw = 1, ngast-1
           CALL schu(nz,o2col,iw,secchi,dto2k,xso2k)
        ENDDO
      ENDIF


! loop through the altitude levels
      DO iz = 1, nz

         igast = 0

! loop through the internal wavelength grid
         DO iw = 1, nwint-1

! if outside Kockarts' grid, use the JPL/Brasseur+Solomon data, if inside
! Kockarts' grid, use the parameterized values from the call to SCHU
           IF ((wlint(iw) .LT. wlgast(1)) .OR.                          &
     &         (wlint(iw) .GT. wlgast(ngast-1))) THEN
              IF (iz .EQ. nz) THEN
                dttmp(iw) = 0.
              ELSE
                dttmp(iw) = xso2int(iw) * 0.2095*cz(iz)
              ENDIF
              xstmp(iw) = xso2int(iw)
           ELSE
              igast = igast+1
              dttmp(iw) = dto2k(iz,igast)
              xstmp(iw) = xso2k(iz,igast)
           ENDIF

! compute the area in each bin (for correct interpolation purposes only!)
           dttmp(iw) = dttmp(iw) * (wlint(iw+1)-wlint(iw))
           xstmp(iw) = xstmp(iw) * (wlint(iw+1)-wlint(iw))

         ENDDO

! interpolate O2 optical depth from the internal grid onto the user grid
         CALL inter3(nw,wl,dtuser, nwint,wlint,dttmp, 0)
         DO iw = 1, nw-1
            dto2(iz,iw) = dtuser(iw)/(wl(iw+1)-wl(iw))
         ENDDO

! interpolate O2 cross section from the internal grid onto the user grid
         CALL inter3(nw,wl,xsuser, nwint,wlint,xstmp, 0)

         DO iw = 1, nw-1
            xso2(iz,iw) = xsuser(iw)/(wl(iw+1)-wl(iw))
         ENDDO

      ENDDO

      IF (LHOOK) CALL DR_HOOK('SETO2',1,ZHOOK_HANDLE)
      RETURN
      ENDSUBROUTINE SETO2
