!MNH_LIC Copyright 2009-2020 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------------
!     ######spl
        SUBROUTINE MNH2LPDM_ECH(TPFILE,TPMETEOFILE)
!	##################################################
!-----------------------------------------------------------------------
!****	MNH2S2_ECH TRAITEMENT D'UNE ECHEANCE.
!
! Auteur   : Francois Bonnardot, DP/SERV/ENV
! Creation : 07.01.2009
! Modifications:
!  P. Wautelet 05/2016-04/2018: new data structures and calls for I/O
!  P. Wautelet 28/05/2018: corrected truncated integer division (1/3 -> 1./3.)
!  P. Wautelet 26/04/2019: replace non-standard FLOAT function by REAL function
!  P. Wautelet 05/11/2020: correct I/O of MNH2LPDM
!-----------------------------------------------------------------------
!
!*	0.  DECLARATIONS.
!	    -------------
!
!*	0.1 Modules.
!
!
!
USE MODD_DIM_n
USE MODD_IO,               ONLY: TFILEDATA
USE MODD_TIME_n
USE MODD_GRID_n
!
USE MODD_CST
USE MODD_PARAMETERS
USE MODD_TIME
!
USE MODD_MNH2LPDM
!
use modd_field,            only: tfielddata, TYPEREAL
USE MODE_IO_FILE,          only: IO_File_close, IO_File_open
USE MODE_IO_FIELD_READ,    only: IO_Field_read
USE MODE_IO_MANAGE_STRUCT, only: IO_File_add2list
!
USE MODI_INI_CST
!
IMPLICIT NONE
!
!
!*	0.2 Arguments.
!
TYPE(TFILEDATA),POINTER,INTENT(INOUT) :: TPFILE
TYPE(TFILEDATA),POINTER,INTENT(IN)    :: TPMETEOFILE
!
!
!*	0.3 Variables locales.
!
CHARACTER(LEN=100)   :: YFTURB                       ! Stockage champs de turbulence.
INTEGER              :: IFTURB
INTEGER              :: IFMTO,IREP
INTEGER              :: ICURAA,ICURMM,ICURJJ         ! Date  courante.
INTEGER              :: ICURHH,ICURMN,ICURSS         ! Heure courante.
INTEGER              :: JI,JJ,JK
TYPE(DATE_TIME)      :: TZDTCUR
type(tfielddata)        :: tzfield
TYPE(TFILEDATA),POINTER :: TZFILE
!
!
!
!
!*	1.  INITIALISATION.
!	    ---------------
!
!*	1.1 Blabla.
!
TZFILE => NULL()
IFMTO = TPMETEOFILE%NLU
!
!*	2.  LECTURE DES DONNEES MESO-NH DE BASE.
!	    ------------------------------------
!
!*	2.1 Ouverture du fichier Meso-NH.
!
CALL IO_File_open(TPFILE)
!
!*	2.2 Date et heure courante.
!
CALL IO_Field_read(TPFILE,'DTCUR',TZDTCUR)
! 
ICURAA=MOD(TZDTCUR%nyear,100)  ! Annee sur 2 caracteres.
ICURMM=TZDTCUR%nmonth
ICURJJ=TZDTCUR%nday
ICURSS=NINT(TZDTCUR%xtime)
!
ICURMN = NINT( (REAL(ICURSS)/60.0)/5.0 )*5   ! Heure arrondie a 5 minutes pres.
ICURSS = 0
ICURHH =ICURMN/60
ICURMN =ICURMN-ICURHH*60
!
print*, '%%% MNH2LPDM2_ECH Date et heure des donnees :'
print 20300, ICURJJ,ICURMM,ICURAA,ICURHH,ICURMN,ICURSS
20300 FORMAT(I2.2,'/',I2.2,'/',I4.4,' ',I2.1,'h',I2.1,'mn',I2.1,'sec')
!
!
!
!*	2.3 Lecture des champs Meso-NH de base.
!
CALL IO_Field_read(TPFILE,'UT',     XUT)
CALL IO_Field_read(TPFILE,'VT',     XVT)
CALL IO_Field_read(TPFILE,'WT',     XWT)
CALL IO_Field_read(TPFILE,'THT',    XTHT)
CALL IO_Field_read(TPFILE,'TKET',   XTKET)

tzfield%cmnhname  = 'LM'
tzfield%clongname = ''
tzfield%cunits    = 'm'
tzfield%cdir      = 'XY'
tzfield%ccomment  = 'Mixing length'
tzfield%ngrid     = 1
tzfield%ntype     = TYPEREAL
tzfield%ndims     = 3
CALL IO_Field_read(TPFILE, tzfield, XLM)

tzfield%cmnhname  = 'THW_FLX'
tzfield%clongname = ''
tzfield%cunits    = 'K s-1' !correct?
tzfield%cdir      = 'XY'
tzfield%ccomment  = 'Conservative potential temperature vertical flux'
tzfield%ngrid     = 4
tzfield%ntype     = TYPEREAL
tzfield%ndims     = 3
CALL IO_Field_read(TPFILE, tzfield, XWPTHP)

tzfield%cmnhname  = 'DISS'
tzfield%clongname = ''
tzfield%cunits    = '' !TODO: set units
tzfield%cdir      = 'XY'
tzfield%ccomment  = 'X_Y_Z_DISS'
tzfield%ngrid     = 1
tzfield%ntype     = TYPEREAL
tzfield%ndims     = 3
CALL IO_Field_read(TPFILE, tzfield, XDISSIP)

tzfield%cmnhname  = 'FMU'
tzfield%clongname = ''
tzfield%cunits    = 'kg m-1 s-2'
tzfield%cdir      = 'XY'
tzfield%ccomment  = 'X_Y_FMU'
tzfield%ngrid     = 4
tzfield%ntype     = TYPEREAL
tzfield%ndims     = 2
CALL IO_Field_read(TPFILE, tzfield, XSFU)

tzfield%cmnhname  = 'FMV'
tzfield%clongname = ''
tzfield%cunits    = 'kg m-1 s-2'
tzfield%cdir      = 'XY'
tzfield%ccomment  = 'X_Y_FMV'
tzfield%ngrid     = 4
tzfield%ntype     = TYPEREAL
tzfield%ndims     = 2
CALL IO_Field_read(TPFILE, tzfield, XSFV)

CALL IO_Field_read(TPFILE,'INPRT',  XINRT)
CALL IO_Field_read(TPFILE,'RVT',    XRMVT)
CALL IO_Field_read(TPFILE,'RCT',    XRMCT)
CALL IO_Field_read(TPFILE,'RRT',    XRMRT)
!
!              Lecture des donnees Meso-NH terminee.'
!
!*	2.4 Fermeture du fichier Meso-NH.
!
CALL IO_File_close(TPFILE)
!
!
!*	3.  PREPARATION DES DONNEES.
!	    ------------------------
!
!
!*	3.2 Niveaux altitude "hors-sol" (1:NKMAX).
!
XSU(:,:,1:NKMAX) = XUT(NSIB:NSIE,NSJB:NSJE,NKB:NKE)
XSV(:,:,1:NKMAX) = XVT(NSIB:NSIE,NSJB:NSJE,NKB:NKE)
XSW(:,:,1:NKMAX) = XWT(NSIB:NSIE,NSJB:NSJE,NKB:NKE)
XSTH(:,:,1:NKMAX) = XTHT(NSIB:NSIE,NSJB:NSJE,NKB:NKE)
XSTKE(:,:,1:NKMAX) = XTKET(NSIB:NSIE,NSJB:NSJE,NKB:NKE)
XSLM(:,:,1:NKMAX) = XLM(NSIB:NSIE,NSJB:NSJE,NKB:NKE)
XSDISSIP(:,:,1:NKMAX) = XDISSIP(NSIB:NSIE,NSJB:NSJE,NKB:NKE)
XSINRT(:,:) = XINRT(NSIB:NSIE,NSJB:NSJE)
XSWPTHP(:,:,1:NKMAX) = XWPTHP(NSIB:NSIE,NSJB:NSJE,NKB:NKE)
XSRMV(:,:,1:NKMAX) = XRMVT(NSIB:NSIE,NSJB:NSJE,NKB:NKE)
XSRMC(:,:,1:NKMAX) = XRMCT(NSIB:NSIE,NSJB:NSJE,NKB:NKE)
XSRMR(:,:,1:NKMAX) = XRMRT(NSIB:NSIE,NSJB:NSJE,NKB:NKE)
XSSFU(:,:) = XSFU(NSIB:NSIE,NSJB:NSJE)
XSSFV(:,:) = XSFV(NSIB:NSIE,NSJB:NSJE)
!
!
!*	4.  CALCULS DES TEMPS LAGRANGIENS ET VARIANCES DU VENT POUR LPDM.
!	    ------------------------------------------------------------
!
      XRVSRD  = XRV/XRD
!
      XSUSTAR (:,:)   = XUNDEF
      XSLMO   (:,:)   = XUNDEF
      XSHMIX  (:,:)   = XUNDEF
      XSWSTAR (:,:)   = XUNDEF
      XSSIGU  (:,:,:) = XUNDEF
      XSSIGW  (:,:,:) = XUNDEF
      XSTIMEU (:,:,:) = XUNDEF
      XSTIMEW (:,:,:) = XUNDEF
!
      DO JI=1,NSIMAX ; DO JJ=1,NSJMAX
        !
        !* Temperature potentielle virtuelle.
        !
        XSTHETAV(:)=1.0+XSRMV(JI,JJ,:)+XSRMC(JI,JJ,:)+XSRMR(JI,JJ,:)
        XSTHETAV(:) = XSTH(JI,JJ,:)*(1.0+XSRMV(JI,JJ,:)*XRVSRD)/XSTHETAV(:)
        !
        !* ZHMIX Hauteur de melange.
        !
        XTHSOL       = XSTHETAV(1)+0.5
        XSHMIX(JI,JJ) = 0.0
        DO JK=2,NKMAX
           IF ( XSTHETAV(JK).GT.XTHSOL ) THEN
              XSHMIX(JI,JJ)     = XSHAUT  (JK-1)  &
                +( XSHAUT  (JK) - XSHAUT  (JK-1) )  &
                /( XSTHETAV(JK) - XSTHETAV(JK-1) )  &
                *( XTHSOL      - XSTHETAV(JK-1) )
              EXIT
           ENDIF
        END DO
        XSHMIX(JI,JJ)=MAX(XSHMIX(JI,JJ),50.0)

        !
        !* XSUSTAR Vitesse de frottement.
        !
        XSUSTAR(JI,JJ) = XSSFU(JI,JJ)*XSSFU(JI,JJ) &
                        +XSSFV(JI,JJ)*XSSFV(JI,JJ)
        XSUSTAR(JI,JJ) = SQRT(SQRT(XSUSTAR(JI,JJ)))
        !
        !
        !
        !* XSLMO Longueur de Monin-Obukhov.
        !
        IF (XSWPTHP(JI,JJ,1).NE.0.) THEN
         XSLMO(JI,JJ)= -XSTHETAV(1)*(XSUSTAR(JI,JJ)**3) &
                     / (XKARMAN*XG*XSWPTHP(JI,JJ,1))
        ENDIF
        !
        !
        !* XSWSTAR Vitesse Verticale Convective.
        !
        XSWSTAR(JI,JJ)=XG/XSTHETAV(1)*XSWPTHP(JI,JJ,1)*XSHMIX(JI,JJ)
        XSWSTAR(JI,JJ)=SIGN(1.,XSWSTAR(JI,JJ)) &
                               * ( ABS(XSWSTAR(JI,JJ))**(1./3.))
        !
        !
        IF (CTURBPARAM=="HANNA".OR.CTURBPARAM=="HANNABIS") THEN
        !
        IF ((XSLMO(JI,JJ).GT.0).AND.(XSLMO(JI,JJ).LE.300)) THEN
           !
           !* Conditions stables.
           !
           !*	XSSIGU,XSSIGW <u'2>**0.5, <w'2>**0.5
           DO JK=1,NKMAX
           IF (XSHAUT(JK).LT.XSHMIX(JI,JJ)) THEN
           XSSIGU(JI,JJ,JK) = SQRT( 0.5 *                              &
                ((2.0*(1-XSHAUT(JK)/XSHMIX(JI,JJ))*XSUSTAR(JI,JJ))**2) &
              + ((1.3*(1-XSHAUT(JK)/XSHMIX(JI,JJ))*XSUSTAR(JI,JJ))**2) )
           XSSIGW(JI,JJ,JK) = 1.3*(1-XSHAUT(JK)/XSHMIX(JI,JJ))         &
                             *XSUSTAR(JI,JJ)
           ELSE
           XSSIGU(JI,JJ,JK) = 0.001
           XSSIGW(JI,JJ,JK) = 0.001
           ENDIF
           ENDDO
           ! 
           XSSIGU(JI,JJ,:)=MAX(0.001,XSSIGU(JI,JJ,:))
           XSSIGW(JI,JJ,:)=MAX(0.001,XSSIGW(JI,JJ,:))
           !
           !* Lagrangian time scale
           XSTIMEU(JI,JJ,:) = 0.11*XSHMIX(JI,JJ)/XSSIGU(JI,JJ,:)  &
                             *SQRT( XSHAUT(:)/XSHMIX(JI,JJ) )
           XSTIMEW(JI,JJ,:) = 0.10*XSHMIX(JI,JJ)/XSSIGW(JI,JJ,:)  &
                             *( XSHAUT(:)/XSHMIX(JI,JJ) )**0.8
           !
           !
        ENDIF
        !
        !
        IF (ABS(XSLMO(JI,JJ)).GT.300) THEN
           !
           !* Conditions neutres.
           !
           !* XSSIGU,XSSIGW <u'2>**0.5, <w'2>**0.5
           XSSIGU(JI,JJ,:)=SQRT( 0.5 *                                &
        ((2.0*XSUSTAR(JI,JJ)*EXP(-3*XSCORIOZ(JI,JJ)*XSHAUT(:)/XSUSTAR(JI,JJ)))**2)  &
     +  ((1.3*XSUSTAR(JI,JJ)*EXP(-2*XSCORIOZ(JI,JJ)*XSHAUT(:)/XSUSTAR(JI,JJ)))**2) )
           XSSIGW(JI,JJ,:)=1.3*XSUSTAR(JI,JJ)*EXP(-2*XSCORIOZ(JI,JJ)*XSHAUT(:)/XSUSTAR(JI,JJ))
           XSSIGU(JI,JJ,:)=MAX(0.001,XSSIGU(JI,JJ,:))
           XSSIGW(JI,JJ,:)=MAX(0.001,XSSIGW(JI,JJ,:))
           !
           !* lagrangian time scale
           XSTIMEU(JI,JJ,:) = 0.5*XSHAUT(:)/   &
                 (XSSIGW(JI,JJ,:)*(1.+15.0*XSCORIOZ(JI,JJ)*XSHAUT(:)/XSUSTAR(JI,JJ)))
           XSTIMEW(JI,JJ,:) = XSTIMEU(JI,JJ,:)
           !
        ENDIF 
        !
        ! 
        IF ((XSLMO(JI,JJ).LT.0).AND.(XSLMO(JI,JJ).GE.-300)) THEN
           !
           !* Conditions instables.
           !
           !* XSSIGU,XSSIGW <u'2>**0.5, <w'2>**0.5
           ! 
           IF (CTURBPARAM=="HANNA") THEN
           !
             DO JK=1,NKMAX
             IF (XSHAUT(JK).LE.XSHMIX(JI,JJ)) THEN
             XSSIGU(JI,JJ,JK)=XSUSTAR(JI,JJ)         &
                        * (12+0.5*XSHMIX(JI,JJ)/ABS(XSLMO(JI,JJ)))**(1./3.)
             ELSE
             XSSIGU(JI,JJ,JK)=0.001
             ENDIF
             ENDDO
             ! 
             DO JK=1,NKMAX
             !IF (XSHAUT(JK).LE.XSHMIX(JI,JJ)) THEN
             !  XSSIGW(JI,JJ,JK)=SQRT(  1.2*XSWSTAR(JI,JJ)**2 &
             !               *(1-0.9*XSHAUT(JK)/XSHMIX(JI,JJ)) &
             !               *(XSHAUT(JK)/XSHMIX(JI,JJ))**(2/3) &
             !            +   (1.8-1.4*XSHAUT(JK)/XSHMIX(JI,JJ)) &
             !               *XSUSTAR(JI,JJ)**2 )  
             !ELSE
             IF (XSHAUT(JK).LE.0.4*XSHMIX(JI,JJ)) THEN
               XSSIGW(JI,JJ,JK)=0.763*(XSHAUT(JK)/XSHMIX(JI,JJ))**0.175
             ELSE IF (XSHAUT(JK).LE.0.96*XSHMIX(JI,JJ)) THEN
               XSSIGW(JI,JJ,JK)=0.722*XSWSTAR(JI,JJ)* &
                                (1-XSHAUT(JK)/XSHMIX(JI,JJ))**0.207
             ELSE IF (XSHAUT(JK).LE.XSHMIX(JI,JJ)) THEN
               XSSIGW(JI,JJ,JK)=0.37*XSWSTAR(JI,JJ)
             ELSE
               XSSIGW(JI,JJ,JK)=0.001
             ENDIF
             ENDDO
             !
             XSSIGU(JI,JJ,:)=MAX(0.001,XSSIGU(JI,JJ,:))
             XSSIGW(JI,JJ,:)=MAX(0.001,XSSIGW(JI,JJ,:))
             !
             !* Lagrangian time scale
             XSTIMEU(JI,JJ,:) = 0.15*XSHMIX(JI,JJ)/XSSIGU(JI,JJ,:)
             DO JK=1,NKMAX
                 IF (XSHAUT(JK).LE.(0.1*XSHMIX(JI,JJ))) THEN
                   IF ( XSHAUT(JK).LT.(XSZ0(JI,JJ)-XSLMO(JI,JJ)) ) THEN
                     XSTIMEW(JI,JJ,JK) = 0.1*XSHAUT(JK)/XSSIGW(JI,JJ,JK)   &
                      / ( 0.55 - 0.38*(XSHAUT(JK)-XSZ0(JI,JJ))/ABS(XSLMO(JI,JJ)))
                   ELSE
                     XSTIMEW(JI,JJ,JK) = 0.59*XSHAUT(JK)/XSSIGW(JI,JJ,JK)
                   ENDIF
                 ELSE
                     XSTIMEW(JI,JJ,JK) = 0.15*XSHMIX(JI,JJ)/XSSIGW(JI,JJ,JK) &
                   *( 1.-EXP(-5*XSHAUT(JK)/XSHMIX(JI,JJ)) )
                 ENDIF
             END DO
           !
           ELSE IF (CTURBPARAM=="HANNABIS") THEN
             !* sigmas
             XSSIGW(JI,JJ,:) = SQRT(2./3.*XSTKE(JI,JJ,:))
             XSSIGU(JI,JJ,:) = XSSIGW(JI,JJ,:)
             !* Temps Lagrangien
             DO JK=1,NKMAX
               IF (XSHAUT(JK).LE.XSHMIX(JI,JJ)) THEN
               XSTIMEU(JI,JJ,JK)=0.17*XSHMIX(JI,JJ)/XSSIGU(JI,JJ,JK)
               XSTIMEW(JI,JJ,JK)=0.2*XSHMIX(JI,JJ)/XSSIGW(JI,JJ,JK)* &
               (1-EXP(-4*XSHAUT(JK)/XSHMIX(JI,JJ)) &
                 -0.0003*EXP(8*XSHAUT(JK)/XSHMIX(JI,JJ)))
               ELSE IF (XSHAUT(JK).LE.XSHMIX(JI,JJ)*1.2) THEN
               XSTIMEU(JI,JJ,JK)= &
                (1-(XSHAUT(JK)-XSHAUT(JK-1))/(XSHAUT(JK+1)-XSHAUT(JK-1)))* &
                XSTIMEU(JI,JJ,JK-1) &
                +(XSHAUT(JK)-XSHAUT(JK-1))/(XSHAUT(JK+1)-XSHAUT(JK-1))*10000.0
               XSTIMEW(JI,JJ,JK)= &
                (1-(XSHAUT(JK)-XSHAUT(JK-1))/(XSHAUT(JK+1)-XSHAUT(JK-1)))* &
                XSTIMEW(JI,JJ,JK-1)  &
                +(XSHAUT(JK)-XSHAUT(JK-1))/(XSHAUT(JK+1)-XSHAUT(JK-1))*10000.0
               ELSE
               XSTIMEU(JI,JJ,JK)=10000.0
               XSTIMEW(JI,JJ,JK)=10000.0
               ENDIF
             ENDDO
             !
           ENDIF ! CTURBPARAM=HANNA ou HANNABIS
           !  
        ENDIF ! instable
        !
        ELSE ! CTURBPARAM=="ISOTROPE"
          !
          !*	XSSIGU,XSSIGW <u'2>**0.5, <w'2>**0.5
          !
          XSSIGW(JI,JJ,:) = SQRT(2./3.*XSTKE(JI,JJ,:))
          XSSIGU(JI,JJ,:) = XSSIGW(JI,JJ,:)
          !
          !* Lagrangian time scale
          DO JK=1,NKMAX
            IF (XSHAUT(JK).LE.XSHMIX(JI,JJ)) THEN
            XSTIMEU(JI,JJ,JK)=ABS(2*(XSSIGU(JI,JJ,JK)**2)/(3*XSDISSIP(JI,JJ,JK)))
            XSTIMEW(JI,JJ,JK)=ABS(2*(XSSIGW(JI,JJ,JK)**2)/(3*XSDISSIP(JI,JJ,JK)))
            ELSE IF (XSHAUT(JK).LE.XSHMIX(JI,JJ)*1.2) THEN
            XSTIMEU(JI,JJ,JK)= &
            (1-(XSHAUT(JK)-XSHAUT(JK-1))/(XSHAUT(JK+1)-XSHAUT(JK-1)))*XSTIMEU(JI,JJ,JK-1) &
           +(XSHAUT(JK)-XSHAUT(JK-1))/(XSHAUT(JK+1)-XSHAUT(JK-1))*1000.0
            XSTIMEW(JI,JJ,JK)=XSTIMEU(JI,JJ,JK)
            ELSE
            XSTIMEU(JI,JJ,JK)=1000.0
            XSTIMEW(JI,JJ,JK)=1000.0
            ENDIF
          ENDDO
          !
        ENDIF
        !
        !
     END DO
  END DO
  !
     IF (IGRILLE.EQ.2) THEN
     WRITE(YFTURB,'("TURB_LPDM",5I2.2)') ICURAA,ICURMM,ICURJJ,ICURHH,ICURMN
     CALL IO_File_add2list(TZFILE,YFTURB,'TXT','WRITE')
     CALL IO_File_open(TZFILE)
     IFTURB = TZFILE%NLU
     WRITE(UNIT=IFTURB,FMT='(5A12)') "WSTAR       ","USTAR       ", &
                                     "HMIX        ","LMO         ", &
                                     "WPTHP"
     WRITE(UNIT=IFTURB,FMT='(5F12.5)') XSWSTAR(15,15),XSUSTAR(15,15), &
                                       XSHMIX(15,15),XSLMO(15,15),    &
                                       XSWPTHP(15,15,1)


     WRITE(UNIT=IFTURB,FMT='(8A12)') "HAUT          ","TKE           ", &
                                     "DISS          ","THETA         ", &
                                     "SIGU          ","SIGW          ", &
                                     "TIMEU         ","TIMEW         "
     DO JK=1,NKMAX
     WRITE(UNIT=IFTURB,FMT='(6F12.5,2F12.1)') XSHAUT(JK),XSTKE(15,15,JK), &
                                    XSDISSIP(15,15,JK),XSTH(15,15,JK), &
                                    XSSIGU(15,15,JK),XSSIGW(15,15,JK), &
                                    XSTIMEU(15,15,JK),XSTIMEW(15,15,JK)
  
     ENDDO
     CALL IO_File_close(TZFILE)
     ENDIF
!
               

!
!*	5.  ECRITURES FIC MTO.
!	    ------------------
!
!
DO JK = 1,NKMAX
WRITE(IFMTO) XSU(:,:,JK)        ! Composante zonale du vent.
ENDDO
DO JK = 1,NKMAX
WRITE(IFMTO) XSV(:,:,JK)        ! Composante meridienne du vent.
ENDDO
DO JK = 1,NKMAX
WRITE(IFMTO) XSW(:,:,JK)        ! Vitesse verticale.
ENDDO
DO JK = 1,NKMAX
WRITE(IFMTO) XSTH(:,:,JK)     ! Temperature potentielle.
ENDDO
DO JK = 1,NKMAX
WRITE(IFMTO) XSTKE(:,:,JK)       ! Energie cinetique Turbulence
ENDDO
DO JK = 1,NKMAX
WRITE(IFMTO) (XSSIGU(:,:,JK))**2       ! SigmaU
ENDDO
DO JK = 1,NKMAX
WRITE(IFMTO) (XSSIGU(:,:,JK))**2       ! SigmaV
ENDDO
DO JK = 1,NKMAX
WRITE(IFMTO) (XSSIGW(:,:,JK))**2       ! SigmaW
ENDDO
DO JK = 1,NKMAX
WRITE(IFMTO) XSTIMEU(:,:,JK)       ! Temps lagrangien U
ENDDO
DO JK = 1,NKMAX
WRITE(IFMTO) XSTIMEU(:,:,JK)       ! Temps lagrangien V
ENDDO
DO JK = 1,NKMAX
WRITE(IFMTO) XSTIMEW(:,:,JK)       ! Dissipation de TKE
ENDDO
WRITE(IFMTO) XSINRT
!
END SUBROUTINE MNH2LPDM_ECH
