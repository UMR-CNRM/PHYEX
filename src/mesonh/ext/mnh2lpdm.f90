!MNH_LIC Copyright 2002-2020 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------------
!     ######spl
        PROGRAM MNH2LPDM
!	##############
!-----------------------------------------------------------------------------
!****	MNH2DIF COUPLAGE MESO-NH / SPRAY.
!
!	Auteur   :   Michel Bouzom, DP/SERV/ENV
!	Creation :   16.07.2002
!       Modification  : 07.01.2006 (T.LAUVAUX, adaptation LPDM)
!       Modification  : 04.01.2009 (F. BONNARDOT, DP/SER/ENV )
!  P. Wautelet 05/2016-04/2018: new data structures and calls for I/O
!  P. Wautelet 07/02/2019: force TYPE to a known value for IO_File_add2list
!  P. Wautelet 10/04/2019: replace ABORT and STOP calls by Print_msg
!  P. Wautelet 05/11/2020: correct I/O of MNH2LPDM
!
!-----------------------------------------------------------------------------
!
!
!
!*	0.  DECLARATIONS.
!	    -------------
!
!*	0.1 Modules.
!
USE MODD_CONF,             ONLY : CPROGRAM
USE MODD_IO,               ONLY : TFILEDATA, TFILE_OUTPUTLISTING, TPTR2FILE
use modd_lunit,            only: TLUOUT0
use modd_lunit_n,          only: TLUOUT
USE MODD_MNH2LPDM
!
USE MODE_FIELD,            ONLY: INI_FIELD_LIST, INI_FIELD_SCALARS
USE MODE_IO,               ONLY: IO_Init, IO_Config_set
USE MODE_IO_FILE,          ONLY: IO_File_open, IO_File_close
USE MODE_IO_MANAGE_STRUCT, ONLY: IO_File_add2list
USE MODE_MODELN_HANDLER
use mode_msg
USE MODE_POS
!
USE MODI_INI_CST
USE MODI_MNH2LPDM_ECH
USE MODI_MNH2LPDM_INI
USE MODI_VERSION
!
USE MODN_CONFIO
!
!
!*	0.2 Variables locales.
!
IMPLICIT NONE
!
CHARACTER(LEN=*),PARAMETER :: YFLOG = 'METEO.log'     ! Log filename
CHARACTER(LEN=*),PARAMETER :: YFNML = 'MNH2LPDM1.nam' ! Namelist filename
INTEGER,         PARAMETER :: IVERB = 5
!
INTEGER :: IFNML  ! Unit of namelist
INTEGER :: JFIC
LOGICAL :: GFOUND ! Return code when searching namelist
TYPE(TPTR2FILE),DIMENSION(JPMNHMAX) :: TZFMNH  ! MesoNH files
TYPE(TFILEDATA),POINTER :: TZDATEFILE  => NULL() ! Date file
TYPE(TFILEDATA),POINTER :: TZGRIDFILE  => NULL() ! Grid file
TYPE(TFILEDATA),POINTER :: TZMETEOFILE => NULL() ! Meteo file
TYPE(TFILEDATA),POINTER :: TZLOGFILE   => NULL() ! Log file
TYPE(TFILEDATA),POINTER :: TZNMLFILE   => NULL() ! Namelist file
!
!
!
!
!*	1.  INITIALISATION.
!	    ---------------
!
CPROGRAM='M2LPDM'
CALL GOTO_MODEL(1)
CALL VERSION()
CALL IO_Init()
CALL INI_CST()
CALL INI_FIELD_LIST(1)
CALL INI_FIELD_SCALARS()
!
CALL IO_File_add2list(TLUOUT0,'OUTPUT_LISTING1','OUTPUTLISTING','WRITE')
CALL IO_File_open(TLUOUT0)
!Set output files for PRINT_MSG
TLUOUT              => TLUOUT0
TFILE_OUTPUTLISTING => TLUOUT0
!
!*	1.1 Variables generales.
!
 CFMNH(:) = ''
!
!
!*	1.2 Initialisation routines LL.
!
CALL IO_Init()
!
!
!*	1.3 Ouverture du fichier log.
!
CALL IO_File_add2list(TZLOGFILE,YFLOG,'TXT','WRITE')
CALL IO_File_open(TZLOGFILE)
!
!
!*	1.4 Lecture des namelists.
!
CALL IO_File_add2list(TZNMLFILE,YFNML,'NML','READ')
CALL IO_File_open(TZNMLFILE)
IFNML = TZNMLFILE%NLU

READ(UNIT=IFNML,NML=NAM_TURB)
READ(UNIT=IFNML,NML=NAM_FIC)
print *,'Lecture de NAM_FIC OK.'

CALL POSNAM(IFNML,'NAM_CONFIO',GFOUND)
IF (GFOUND) THEN
  READ(UNIT=IFNML,NML=NAM_CONFIO)
END IF
LCDF4 = .FALSE.
LLFIOUT  = .FALSE.
LLFIREAD = .FALSE.
CALL IO_Config_set()
CALL IO_File_close(TZNMLFILE)
!
!
!*	1.5 Comptage des FM a traiter.
!
IF (LEN_TRIM(CFMNH(1))>0) THEN
   NBMNH=1
   CALL IO_File_add2list(TZFMNH(1)%TZFILE,TRIM(CFMNH(1)),'MNH','READ',KLFITYPE=2,KLFIVERB=IVERB)
   DO WHILE (CFMNH(NBMNH+1).NE.'VIDE')
      NBMNH=NBMNH+1
      CALL IO_File_add2list(TZFMNH(NBMNH)%TZFILE,TRIM(CFMNH(NBMNH)),'MNH','READ',KLFITYPE=2,KLFIVERB=IVERB)
   END DO
   print *,NBMNH,' fichiers a traiter.'
ELSE
  call Print_msg( NVERB_FATAL, 'GEN', 'MNH2LPDM', 'no CFMNH file given' )
END IF
!
!
!
!
!*	2.  TRAITEMENTS.
!	    ------------
!
!*	2.1 Ouverture des fichiers METEO et GRILLE et DATE.
!
CALL IO_File_add2list(TZGRIDFILE,CFGRI,'TXT','WRITE')
CALL IO_File_open(TZGRIDFILE)
CALL IO_File_add2list(TZDATEFILE,CFDAT,'TXT','WRITE')
CALL IO_File_open(TZDATEFILE)
!
!
!*	2.2 Preparation du couplage.
!
CALL MNH2LPDM_INI(TZFMNH(1)%TZFILE,TZFMNH(NBMNH)%TZFILE,TZLOGFILE,TZGRIDFILE,TZDATEFILE)
!
!
!*	2.3 Traitement des echeances.
!
DO JFIC=1,NBMNH
   print*,"CFMTO(JFIC)=",CFMTO(JFIC)
   CALL IO_File_add2list(TZMETEOFILE,CFMTO(JFIC),'METEO','WRITE')
   CALL IO_File_open(TZMETEOFILE)
   CALL MNH2LPDM_ECH(TZFMNH(JFIC)%TZFILE,TZMETEOFILE)
   print*,"CLOSE_LL(CFMTO(JFIC)"
   CALL IO_File_close(TZMETEOFILE)
   TZMETEOFILE => NULL()
END DO
!
!
!*	2.4 Fermeture des fichiers, METEO, GRILLE et LOG.
!
CALL IO_File_close(TZGRIDFILE)
CALL IO_File_close(TZDATEFILE)
CALL IO_File_close(TZLOGFILE)
!
!
!
END PROGRAM MNH2LPDM
