!MNH_LIC Copyright 2009-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------------
!     ######spl
        SUBROUTINE MNH2LPDM_INI(TPFILE1,TPFILE2,TPLOGFILE,TPGRIDFILE,TPDATEFILE)
!--------------------------------------------------------------------------
!* MNH2S2_INI    : INITIALISATION DU COUPLAGE MESO-NH / LPDM.
!
!  Auteur        : Francois BONNARDOT, DP/SERV/ENV
!  Creation      : 04.01.2009 (mnh2s2_ini.f90)
!
!
!	Arguments explicites.
!	---------------------
! TPFILE1,TPFILE2 First and last files to treat
! TPLOGFILE       Log file
! TPGRIDFILE      Grid file
! TPDATEFILE      Date file
!
! Modifications:
!  P. Wautelet 05/2016-04/2018: new data structures and calls for I/O
!  P. Wautelet 26/04/2019: replace non-standard FLOAT function by REAL function
!  P. Wautelet 05/11/2020: correct I/O of MNH2LPDM
!--------------------------------------------------------------------------
!
!
!
!*	0.  INITIALISATION.
!	    ---------------
!
!*	0.1 Modules.
!
USE MODD_CST
USE MODD_DIM_n
use modd_field,         only: tfieldmetadata, TYPEREAL
USE MODD_GRID
USE MODD_GRID_n
USE MODD_IO,            ONLY: TFILEDATA
USE MODD_LUNIT
USE MODD_MNH2LPDM
USE MODD_PARAMETERS
USE MODD_TIME
USE MODD_TIME_n
!
USE MODE_DATETIME
USE MODE_GRIDPROJ
USE MODE_INI_CST,       ONLY: INI_CST
USE MODE_IO_FILE,       only: IO_File_close, IO_File_open
USE MODE_IO_FIELD_READ, only: IO_Field_read
USE MODE_MODELN_HANDLER
!
USE MODI_READ_HGRID
USE MODI_XYTOLATLON
!
!*	0.2 Arguments.
!
IMPLICIT NONE
!
TYPE(TFILEDATA),POINTER,INTENT(INOUT) :: TPFILE1,TPFILE2
TYPE(TFILEDATA),POINTER,INTENT(IN)    :: TPLOGFILE
TYPE(TFILEDATA),POINTER,INTENT(IN)    :: TPGRIDFILE
TYPE(TFILEDATA),POINTER,INTENT(IN)    :: TPDATEFILE
!
!
!*	0.3 Variables locales.
!
CHARACTER(LEN=28)     :: YNAME,YDAD     ! Noms du FM et de son papa.
CHARACTER(LEN=2)      :: YSTORAGE       ! Type de variable.
!
REAL                  :: ZECHEANCE1,ZECHEANCE2     ! dist temp date modele - date courante
INTEGER               :: IHHMDL,IMNMDL,ISSMDL ! h - mn - s du model
INTEGER               :: IHHCUR1,IMNCUR1,ISSCUR1
INTEGER               :: IHHCUR2,IMNCUR2,ISSCUR2
CHARACTER(LEN=14)     :: YDATMDL,YDATCUR1,YDATCUR2
!
REAL                  :: XLATOR,XLONOR,XPTLAT,XPTLON
REAL                  :: XXPTSOMNH,XYPTSOMNH
INTEGER               :: JI,JJ,JK,a
INTEGER               :: b,c,I
INTEGER, DIMENSION(:),   ALLOCATABLE   :: TAB1D
INTEGER, DIMENSION(:,:), ALLOCATABLE   :: TAB2D
TYPE(DATE_TIME)         :: TZDTCUR1,TZDTCUR2,TZDTEXP1
INTEGER                :: IFDAT,IFGRI,IFLOG
type(tfieldmetadata)   :: tzfield
!
!
!
!*	1.  INITIALISATION.
!	    ---------------
!
IFDAT = TPDATEFILE%NLU
IFGRI = TPGRIDFILE%NLU
IFLOG = TPLOGFILE%NLU
!
CALL INI_CST
!
CALL GOTO_MODEL(1)
!
!
!*	2.  DONNEES MESO-NH.
!	    ----------------
!
!*	2.1 Ouverture du fichier Meso-NH.
!
CALL IO_File_open(TPFILE1)
CALL IO_File_open(TPFILE2)
!
!
!*      2.2 Date et heure du modele.
!
CALL IO_Field_read(TPFILE1,'DTEXP',TZDTEXP1)
CALL IO_Field_read(TPFILE1,'DTCUR',TZDTCUR1)
CALL IO_Field_read(TPFILE2,'DTCUR',TZDTCUR2)
!
CALL DATETIME_DISTANCE(TZDTEXP1,TZDTCUR1,ZECHEANCE1)
CALL DATETIME_DISTANCE(TZDTEXP1,TZDTCUR2,ZECHEANCE2)
!
IHHMDL=INT(TZDTEXP1%xtime/3600)
IMNMDL=INT((TZDTEXP1%xtime-IHHMDL*3600)/60)
ISSMDL=INT(TZDTEXP1%xtime-IHHMDL*3600-IMNMDL*60)
IHHCUR1=INT(TZDTCUR1%xtime/3600)
IMNCUR1=INT((TZDTCUR1%xtime-IHHCUR1*3600)/60)
ISSCUR1=INT(TZDTCUR1%xtime-IHHCUR1*3600-IMNCUR1*60)
IHHCUR2=INT(TZDTCUR2%xtime/3600)
IMNCUR2=INT((TZDTCUR2%xtime-IHHCUR2*3600)/60)
ISSCUR2=INT(TZDTCUR2%xtime-IHHCUR2*3600-IMNCUR2*60)
!
WRITE(YDATMDL, '(I4.4,5I2.2)') TZDTEXP1%nyear, TZDTEXP1%nmonth, TZDTEXP1%nday, &
                               IHHMDL, IMNMDL, ISSMDL
WRITE(YDATCUR1,'(I4.4,5I2.2)') TZDTCUR1%nyear, TZDTCUR1%nmonth, TZDTCUR1%nday, &
                               IHHCUR1, IMNCUR1, ISSCUR1
WRITE(YDATCUR2,'(I4.4,5I2.2)') TZDTCUR2%nyear, TZDTCUR2%nmonth, TZDTCUR2%nday, &
                               IHHCUR2, IMNCUR2, ISSCUR2
! 
NMDLAA=MOD( TZDTEXP1%nyear, 100 )  ! Annee arrondi a 2 chiffres.
NMDLMM=TZDTEXP1%nmonth
NMDLJJ=TZDTEXP1%nday
NMDLSS=NINT(TZDTEXP1%xtime)
!
!*	Heure du modele arrondie a 5 minutes pres.
!
NMDLMN = NINT( (REAL(NMDLSS)/60.0)/5.0 )*5
NMDLSS = 0
NMDLHH =NMDLMN/60
NMDLMN =NMDLMN-NMDLHH*60
!
!*	2.3 Grille horizontale.
!
CALL READ_HGRID(1,TPFILE1,YNAME,YDAD,YSTORAGE)
IF (YNAME == YDAD) THEN
IGRILLE=1
ELSE
IGRILLE=2
ENDIF
print*,IGRILLE
!
! Lecture grille horizontale
!
NIU=NIMAX+2*JPHEXT
NJU=NJMAX+2*JPHEXT
NIB=1+JPHEXT
NJB=1+JPHEXT
NIE=NIU-JPHEXT
NJE=NJU-JPHEXT
!
!
!*	2.4 Nombre de niveaux-verticaux.
!
CALL IO_Field_read(TPFILE1,'KMAX',NKMAX)
!WRITE(IFLOG,*) '%%% MNH2S2_INI Lecture du nombre de niveau OK.'
!
NKU = NKMAX+2*JPVEXT
NKB = 1+JPVEXT
NKE = NKU-JPVEXT
!
!
!*	2.5 Allocations Meso-NH.
!
ALLOCATE( XZHAT(NKU) )
ALLOCATE( XZS(NIU,NJU) )
ALLOCATE( XZ0(NIU,NJU) )
ALLOCATE( XUT(NIU,NJU,NKU))
ALLOCATE( XVT(NIU,NJU,NKU))
ALLOCATE( XWT(NIU,NJU,NKU))
ALLOCATE( XTHT(NIU,NJU,NKU))
ALLOCATE( XTKET(NIU,NJU,NKU))
ALLOCATE( XLM(NIU,NJU,NKU))
ALLOCATE( XDISSIP(NIU,NJU,NKU))
ALLOCATE( XWPTHP(NIU,NJU,NKU))
ALLOCATE( XRMVT(NIU,NJU,NKU))
ALLOCATE( XRMCT(NIU,NJU,NKU))
ALLOCATE( XRMRT(NIU,NJU,NKU))
ALLOCATE( XINRT(NIU,NJU))
ALLOCATE( XSFU(NIU,NJU))
ALLOCATE( XSFV(NIU,NJU))
!
!*	2.6 Decoupage vertical.
!
CALL IO_Field_read(TPFILE1,'ZHAT',XZHAT)
CALL IO_Field_read(TPFILE1,'ZTOP',XZTOP)
!
!*	2.7 Orographie. 
!
CALL IO_Field_read(TPFILE1,'ZS',XZS)
!
!*	2.8 Rugosite Z0. 
!
tzfield = tfieldmetadata( &
  cmnhname  = 'Z0',       &
  clongname = '',         &
  cunits    = 'm',        &
  cdir      = 'XY',       &
  ccomment  = 'X_Y_Z0',   &
  ngrid     = 4,          &
  ntype     = TYPEREAL,   &
  ndims     = 2           )
CALL IO_Field_read(TPFILE1,tzfield,XZ0)
!
XXPTSOMNH=XXHAT(1)+(XXHAT(2)-XXHAT(1))/2
XYPTSOMNH=XYHAT(1)+(XYHAT(2)-XYHAT(1))/2
CALL SM_LATLON(XLATORI,XLONORI,XXPTSOMNH,XYPTSOMNH,XLATOR,XLONOR)
!
!*	2.9  DOMAINE D'EXTRACTION.
!	    ---------------------
!
NSIB   = NIB
NSIE   = NIE
NSJB   = NJB
NSJE   = NJE
!
NSIMAX = NSIE-NSIB+1
NSJMAX = NSJE-NSJB+1
!
!
!*	3. Impression de controle Meso-NH.
!	    -------------------------------
!
!           Domaine horizontal Meso-NH.
!modif 12.2014 : passage a 1 seul domaine MesoNH
!           ---------------------------
WRITE(IFLOG,'(I1,a12)') IGRILLE,'      ngrid '
!WRITE(IFLOG,'(a13)') '2      ngrids'
WRITE(IFLOG,'(a13)') '1      ngrids'
WRITE(IFLOG,'(i4,3x,a6)') NSIMAX,'nx    '
WRITE(IFLOG,'(i4,3x,a6)') NSJMAX,'ny    '
WRITE(IFLOG,'(i4,3x,a6)') NKU-2,'nz    '
WRITE(IFLOG,'(i4,3x,a6)') NKU-3,'nzg   ' 
WRITE(IFLOG,'(a13)') '12     npatch'
WRITE(IFLOG,'(a13)') '0      icloud'
WRITE(IFLOG,'(a11)') '0.0  wlon  '
WRITE(IFLOG,'(a11)') '45.0 rnlat '
WRITE(IFLOG,'(f10.1,3x,a6)') XZHAT(NKE),'s     '
WRITE(IFLOG,'(f8.0,a8)') ZECHEANCE1,'  time1 '
WRITE(IFLOG,'(f8.0,a8)') ZECHEANCE2,'  time2 '
WRITE(IFLOG,'(a13)') '3600    dtmet '
WRITE(IFLOG,'(a13)') 'm       tunits'
WRITE(IFLOG,'(a13)') '12     nvout '
WRITE(IFLOG,'(6x,a8,i4)') 'u       ',1
WRITE(IFLOG,'(6x,a8,i4)') 'v       ',1
WRITE(IFLOG,'(6x,a8,i4)') 'w       ',1
WRITE(IFLOG,'(6x,a8,i4)') 'tp      ',1
WRITE(IFLOG,'(6x,a8,i4)') 'tke     ',1
WRITE(IFLOG,'(6x,a8,i4)') 'uu      ',1
WRITE(IFLOG,'(6x,a8,i4)') 'vv      ',1
WRITE(IFLOG,'(6x,a8,i4)') 'ww      ',1
WRITE(IFLOG,'(6x,a8,i4)') 'tlx     ',1
WRITE(IFLOG,'(6x,a8,i4)') 'tly     ',1
WRITE(IFLOG,'(6x,a8,i4)') 'tlz     ',1
WRITE(IFLOG,'(6x,a8,i4)') 'intopr  ',1
WRITE(IFLOG,*) '  grid structure'
!
!*	4.  FICHIER METEO.
!	    --------------
!
!*	4.1 Allocations.
!
ALLOCATE( XSHAUT(NKMAX))
ALLOCATE( XSREL(NSIMAX,NSJMAX) )
ALLOCATE( XSZ0(NSIMAX,NSJMAX) )
ALLOCATE( XSCORIOZ (NSIMAX,NSJMAX) )
ALLOCATE( XSPHI(NSIMAX,NSJMAX,NKMAX) )
ALLOCATE( XSU(NSIMAX,NSJMAX,NKMAX))
ALLOCATE( XSV(NSIMAX,NSJMAX,NKMAX))
ALLOCATE( XSW(NSIMAX,NSJMAX,NKMAX))
ALLOCATE( XSTH(NSIMAX,NSJMAX,NKMAX))
ALLOCATE( XSTKE(NSIMAX,NSJMAX,NKMAX))
ALLOCATE( XSLM(NSIMAX,NSJMAX,NKMAX))
ALLOCATE( XSDISSIP(NSIMAX,NSJMAX,NKMAX))
ALLOCATE( XSWPTHP(NSIMAX,NSJMAX,NKMAX))
ALLOCATE( XSRMV(NSIMAX,NSJMAX,NKMAX))
ALLOCATE( XSRMC(NSIMAX,NSJMAX,NKMAX))
ALLOCATE( XSRMR(NSIMAX,NSJMAX,NKMAX))
ALLOCATE( XSINRT(NSIMAX,NSJMAX))
ALLOCATE( XSSFU(NSIMAX,NSJMAX))
ALLOCATE( XSSFV(NSIMAX,NSJMAX))
ALLOCATE( XSTIMEW(NSIMAX,NSJMAX,NKMAX))
ALLOCATE( XSTIMEU(NSIMAX,NSJMAX,NKMAX))
ALLOCATE( XSSIGW(NSIMAX,NSJMAX,NKMAX))
ALLOCATE( XSSIGU(NSIMAX,NSJMAX,NKMAX))
ALLOCATE( XSUSTAR(NSIMAX,NSJMAX))
ALLOCATE( XSWSTAR(NSIMAX,NSJMAX))
ALLOCATE( XSHMIX(NSIMAX,NSJMAX))
ALLOCATE( XSLMO(NSIMAX,NSJMAX))
ALLOCATE( XSTHETAV(NKMAX))

!
!   4.2.    Nombre de niveaux en Z
!
XSHAUT(1:NKMAX) = (XZHAT(NKB:NKE)+XZHAT(NKB+1:NKE+1))/2.
print*,"niveaux hauteur"
DO JK=1,NKMAX
print*,XSHAUT(JK)
ENDDO
!
!   4.3.    Calcul du tableau contenant les coef. de coriolis de la grille
!
DO JI=NSIB,NSIE ; DO JJ=NSJB,NSJE
   CALL SM_LATLON(XLATORI,XLONORI,XXHAT(JI),XYHAT(JJ),XPTLAT,XPTLON)
   XSCORIOZ(JI-1,JJ-1)=2.*XOMEGA*SIN(XPTLAT*XPI/180.)
ENDDO ; ENDDO
!
!
!* 4.4 Geometrie de la grille et positionnement.
!
!
! On a besoin du point sud-ouest, c'est-a-dire de l'angle inferieur gauche
! du domaine physique de la maille "en bas a gauche". Ca tombe bien, on
! va travailler avec les XXHAT et les XYHAT directement.
!
XPASXM = XXHAT(2)-XXHAT(1)     ! Pas selon X en metres.
XPASYM = XYHAT(2)-XYHAT(1)     ! Pas selon Y en metres.
ZMAILLE = MAX(XPASXM,XPASYM)
!
!* 4.5 Constantes et champs constants.
!
!* Relief.
!
XSREL(:,:)  =  XZS(NSIB:NSIE,NSJB:NSJE)
!
!* Geopotentiel PHI
!
print*,"Geopotentiel"
DO JK=1,NKMAX
XSPHI(:,:,JK) = (XSREL(:,:)+XSHAUT(JK))*XG
print*,MINVAL(XSPHI(:,:,JK)),MAXVAL(XSPHI(:,:,JK))
ENDDO
!
!* Rugosite.
!
XSZ0(:,:)  =  XZ0(NSIB:NSIE,NSJB:NSJE)
print*,"Rugosite"
print*,MINVAL(XSZ0),MAXVAL(XSZ0)
!
!* 5   FICHIER DATES.
!      -------------
!
WRITE(IFDAT,'(A14)') YDATMDL 
WRITE(IFDAT,'(A14)') YDATCUR1 
WRITE(IFDAT,'(A14)') YDATCUR2
!
!* 5.  FICHIER GRILLE.
!      --------------
!
!
!* 5.1 Infos franchement utiles.
!
WRITE(IFGRI,'(F15.8,1X,A)') &
                  XLON0,   'XLON0   Longitude reference (deg.deci.)'
WRITE(IFGRI,'(F15.8,1X,A)') &
                  XLAT0,   'XLAT0   Latitude  reference (deg.deci.)'
WRITE(IFGRI,'(F15.8,1X,A)') &
                  XBETA,   'XBETA   Rotation  grille    (deg.deci.)'
WRITE(IFGRI,'(F15.8,1X,A)') XRPK,    'XRPK    Facteur de conicite'
WRITE(IFGRI,'(F15.8,1X,A)') &
                  XLONOR,  'XLONOR  Longitude origine   (deg.deci.)'
WRITE(IFGRI,'(F15.8,1X,A)') &
                  XLATOR,  'XLATOR  Latitude  origine   (deg.deci.)'
WRITE(IFGRI,'(F15.1,1X,A)') XXHAT(1),'XHAT(1) Coord. Cartesienne  (m)'
WRITE(IFGRI,'(F15.1,1X,A)') XXHAT(2),'XHAT(2) Coord. Cartesienne  (m)'
WRITE(IFGRI,'(F15.1,1X,A)') XYHAT(1),'YHAT(1) Coord. Cartesienne  (m)'
WRITE(IFGRI,'(F15.1,1X,A)') XYHAT(2),'YHAT(2) Coord. Cartesienne  (m)'
!
print*,"GRILLE"
print*,"LON0 : ",XLON0
print*,"LAT0 : ",XLAT0
print*,"BETA : ",XBETA
print*,"RPK  : ",XRPK
print*,"LONOR: ",XLONOR
print*,"LATOR: ",XLATOR
!
!* 5.2  Points de grille x y z zg
!
WRITE(IFLOG,*)NSIMAX,' gridpoints in x direction'
WRITE(IFLOG,'(8f10.0)')XXHAT(NSIB:NSIE)
WRITE(IFLOG,*)NSJMAX,' gridpoints y direction'
WRITE(IFLOG,'(8f10.0)')XYHAT(NSJB:NSJE)
WRITE(IFLOG,*)NKMAX,'  main gridpoints in z direction'
WRITE(IFLOG,'(8f10.2)')XSHAUT(1:NKMAX)
WRITE(IFLOG,'(i4,3x,a38)')NKU-2,'intermediate gridpoints in z direction'
WRITE(IFLOG,'(8f10.2)')XZHAT(2:NKU-1)
WRITE(IFLOG,*)'   =================================================='
!
!     Topographie
!
WRITE(IFLOG,*) 'TERRAIN TOPOGRAPHY'
c=1
a=0
!modif 12/2014 : passage a une grille haute resolution MesoNH, on depasse 99
!300 format(i2,'|',18i4)
300 format(i3,'|',18i5)
!400 format(i2,'|',18(f4.2))
!400 format(i3,'|',18(f5.2))
!301 format(3x,18('__',i2))
301 format(3x,18('__',i3))
ALLOCATE(TAB2D(NSIMAX,NSJMAX))
ALLOCATE(TAB1D(NSIMAX))
DO I=1,NSIMAX
   TAB1D(I)=I
ENDDO
TAB2D(:,:) = NINT(XSREL(:,:))
DO WHILE (c.lt.(NSIMAX+1))
  DO b=NSJB,NSJE
     IF ((c+17).LT.(NSIMAX+1)) then
        a=NSJMAX-b+NSJB
        WRITE(IFLOG,300) a,TAB2D(c:c+17,a)
     ELSE  
        a=NSJMAX-b+NSJB
        WRITE(IFLOG,300) a,TAB2D(c:NSIMAX,a)
     ENDIF
  ENDDO
IF ((c+17).LT.(NSIMAX+1)) then
   WRITE(IFLOG,301) TAB1D(c:c+17)
ELSE
   WRITE(IFLOG,301) TAB1D(c:NSIMAX)
ENDIF

c=c+18
ENDDO
!
DEALLOCATE(TAB2D)
DEALLOCATE(TAB1D)
DEALLOCATE(XZS)
DEALLOCATE(XZ0)
DEALLOCATE(XZHAT)
!
!	 Fermeture du fichier Meso-NH.
!
CALL IO_File_close(TPFILE1)
CALL IO_File_close(TPFILE2)
!
!
!-------------------------------------------'
print*,'          FIN MNH2LPDM_INI'
!-------------------------------------------'
!
!
END SUBROUTINE MNH2LPDM_INI
