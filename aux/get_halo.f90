!MNH_LIC Copyright 1994-2025 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
! Modifications:
!  P. Wautelet 20/05/2019: add name argument to ADDnFIELD_ll + new ADD4DFIELD_ll subroutine
!  P. Wautelet 18/07/2019: add optional dummy argument with name of the field
!-----------------------------------------------------------------
!     ####################
      MODULE MODI_GET_HALO
!     ####################
!
IMPLICIT NONE
INTERFACE
!
SUBROUTINE GET_HALO_PHY(D,PSRC,OONGPU)
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
IMPLICIT NONE
!
TYPE(DIMPHYEX_t),        INTENT(IN)   :: D
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)  :: PSRC    ! variable at t
LOGICAL , OPTIONAL :: OONGPU
!
END SUBROUTINE GET_HALO_PHY
!
   SUBROUTINE GET_HALO2(PSRC, TP_PSRC_HALO2_ll, HNAME)
     !
     USE MODD_ARGSLIST_ll, ONLY : HALO2LIST_ll
     !
     IMPLICIT NONE
     !
     REAL, DIMENSION(:,:,:), INTENT(IN)  :: PSRC    ! variable at t
     TYPE(HALO2LIST_ll), POINTER         :: TP_PSRC_HALO2_ll          ! halo2 for SRC
     character(len=*), optional, intent(in) :: HNAME ! Name of the field to be added
     !
   END SUBROUTINE GET_HALO2
END INTERFACE
!
INTERFACE 
   SUBROUTINE GET_HALO(PSRC, HDIR, HNAME)
     IMPLICIT NONE
     !
     REAL, DIMENSION(:,:,:), INTENT(IN)  :: PSRC    ! variable at t
     CHARACTER(len=4), OPTIONAL :: HDIR ! to send only halo on X or Y direction
     character(len=*), optional, intent(in) :: HNAME ! Name of the field to be added
     !
   END SUBROUTINE GET_HALO
END INTERFACE
! 
#ifdef MNH_OPENACC
INTERFACE
   SUBROUTINE GET_HALO2_D(PSRC, TP_PSRC_HALO2_ll, HNAME)
     !
     USE MODD_ARGSLIST_ll, ONLY : HALO2LIST_ll
     !
     IMPLICIT NONE
     !
     REAL, DIMENSION(:,:,:), INTENT(INOUT)  :: PSRC    ! variable at t
     TYPE(HALO2LIST_ll), POINTER         :: TP_PSRC_HALO2_ll          ! halo2 for SRC
     character(len=*), optional, intent(in) :: HNAME ! Name of the field to be added
     !
   END SUBROUTINE GET_HALO2_D
END INTERFACE

INTERFACE
   SUBROUTINE GET_HALO2_DD(PSRC, TP_PSRC_HALO2_ll, HNAME)
     !
     USE MODD_ARGSLIST_ll, ONLY : HALO2LIST_ll
     !
     IMPLICIT NONE
     !
     REAL, DIMENSION(:,:,:), INTENT(INOUT)  :: PSRC    ! variable at t
     TYPE(HALO2LIST_ll), POINTER         :: TP_PSRC_HALO2_ll          ! halo2 for SRC
     character(len=*), optional, intent(in) :: HNAME ! Name of the field to be added
     !
   END SUBROUTINE GET_HALO2_DD
END INTERFACE

INTERFACE
   SUBROUTINE GET_HALO2_DF(PSRC, TP_PSRC_HALO2_ll, HNAME)
     !
     USE MODD_ARGSLIST_ll, ONLY : HALO2LIST_ll
     !
     IMPLICIT NONE
     !
     REAL, DIMENSION(:,:,:), INTENT(INOUT)  :: PSRC    ! variable at t
     TYPE(HALO2LIST_ll), POINTER         :: TP_PSRC_HALO2_ll          ! halo2 for SRC
     character(len=*), optional, intent(in) :: HNAME ! Name of the field to be added
     !
   END SUBROUTINE GET_HALO2_DF
END INTERFACE

INTERFACE
   SUBROUTINE GET_HALO_D(PSRC, HDIR, HNAME)
     IMPLICIT NONE
     !
     REAL, DIMENSION(:,:,:), INTENT(INOUT)   :: PSRC    ! variable at t
     CHARACTER(len=4), OPTIONAL :: HDIR ! to send only halo on X or Y direction
     character(len=*), optional, intent(in) :: HNAME ! Name of the field to be added
     !
   END SUBROUTINE GET_HALO_D
END INTERFACE
!
INTERFACE
   SUBROUTINE GET_HALO_START_D(PSRC,KNB_REQ,KREQ,&
     PZNORTH_IN , PZSOUTH_IN , PZWEST_IN , PZEAST_IN , &
     PZNORTH_OUT, PZSOUTH_OUT, PZWEST_OUT, PZEAST_OUT, &
     KIIB,KIIE,KIJB,KIJE,KIIU,KIJU,KIKU,KIHALO_1,&
     HDIR)
     
     IMPLICIT NONE
     !
     REAL, DIMENSION(KIIU,KIJU,KIKU), INTENT(INOUT) :: PSRC    ! variable at t
     ! acc declare present (PSRC)
     INTEGER                               :: KNB_REQ ,  KREQ(8)
     REAL            :: PZSOUTH_IN ( KIIB:KIIE   , KIJB:KIJB+KIHALO_1   , KIKU ) ,&
                        PZNORTH_IN ( KIIB:KIIE   , KIJE-KIHALO_1:KIJE   , KIKU ) ,&
                        PZWEST_IN  ( KIIB:KIIB+KIHALO_1   , KIJB:KIJE   , KIKU ) ,&
                        PZEAST_IN  ( KIIE-KIHALO_1:KIIE   , KIJB:KIJE   , KIKU ) ,&
       !
                        PZSOUTH_OUT (   KIIB:KIIE   , 1:KIJB-1 , KIKU ) ,&
                        PZNORTH_OUT (   KIIB:KIIE   , KIJE+1:KIJU , KIKU ) ,&
                        PZWEST_OUT  ( 1:KIIB-1 ,   KIJB:KIJE   , KIKU ) ,&
                        PZEAST_OUT  ( KIIE+1:KIIU ,   KIJB:KIJE   , KIKU ) 
     INTEGER         :: KIIB,KIIE,KIJB,KIJE,KIIU,KIJU,KIKU,KIHALO_1
     CHARACTER(len=4), OPTIONAL :: HDIR ! to send only halo on X or Y direction
     !
   END SUBROUTINE GET_HALO_START_D
END INTERFACE

INTERFACE
   SUBROUTINE GET_HALO_STOP_D(PSRC,KNB_REQ,KREQ,&
     PZNORTH_OUT, PZSOUTH_OUT, PZWEST_OUT, PZEAST_OUT, &
     KIIB,KIIE,KIJB,KIJE,KIIU,KIJU,KIKU,&
     HDIR)
     
     IMPLICIT NONE
     !
     REAL, DIMENSION(KIIU,KIJU,KIKU), INTENT(INOUT) :: PSRC    ! variable at t
     ! acc declare present (PSRC)
     INTEGER                               :: KNB_REQ ,  KREQ(8)
     REAL            :: PZSOUTH_OUT (   KIIB:KIIE   , 1:KIJB-1 , KIKU ) ,&
                        PZNORTH_OUT (   KIIB:KIIE   , KIJE+1:KIJU , KIKU ) ,&
                        PZWEST_OUT  ( 1:KIIB-1 ,   KIJB:KIJE   , KIKU ) ,&
                        PZEAST_OUT  ( KIIE+1:KIIU ,   KIJB:KIJE   , KIKU ) 
     INTEGER         :: KIIB,KIIE,KIJB,KIJE,KIIU,KIJU,KIKU
     CHARACTER(len=4), OPTIONAL :: HDIR ! to send only halo on X or Y direction
     !
   END SUBROUTINE GET_HALO_STOP_D
END INTERFACE

INTERFACE
   SUBROUTINE GET_HALO_DD(PSRC, HDIR, HNAME)
     IMPLICIT NONE
     !
     REAL, DIMENSION(:,:,:), INTENT(INOUT)   :: PSRC    ! variable at t
     CHARACTER(len=4), OPTIONAL :: HDIR ! to send only halo on X or Y direction
     character(len=*), optional, intent(in) :: HNAME ! Name of the field to be added
     !
   END SUBROUTINE GET_HALO_DD
END INTERFACE

INTERFACE
   SUBROUTINE GET_HALO_DDC(PSRC, HDIR, HNAME)
     IMPLICIT NONE
     !
     REAL, DIMENSION(:,:,:), INTENT(INOUT)   :: PSRC    ! variable at t
     CHARACTER(len=4), OPTIONAL :: HDIR ! to send only halo on X or Y direction
     character(len=*), optional, intent(in) :: HNAME ! Name of the field to be added
     !
   END SUBROUTINE GET_HALO_DDC
END INTERFACE

INTERFACE
   SUBROUTINE GET_2D_HALO_DD(PSRC, HDIR, HNAME)
     IMPLICIT NONE
     !
     REAL, DIMENSION(:,:), INTENT(INOUT)   :: PSRC    ! variable at t
     CHARACTER(len=4), OPTIONAL :: HDIR ! to send only halo on X or Y direction
     character(len=*), optional, intent(in) :: HNAME ! Name of the field to be added
     !
   END SUBROUTINE GET_2D_HALO_DD
END INTERFACE

INTERFACE
   SUBROUTINE GET_2D_HALO_DDC(PSRC, HDIR, HNAME)
     IMPLICIT NONE
     !
     REAL, DIMENSION(:,:), INTENT(INOUT)   :: PSRC    ! variable at t
     CHARACTER(len=4), OPTIONAL :: HDIR ! to send only halo on X or Y direction
     character(len=*), optional, intent(in) :: HNAME ! Name of the field to be added
     !
   END SUBROUTINE GET_2D_HALO_DDC
END INTERFACE
#endif
!
INTERFACE
   SUBROUTINE DEL_HALO2_ll(TPHALO2LIST)
     !
     USE MODD_ARGSLIST_ll, ONLY : HALO2LIST_ll
     !
     IMPLICIT NONE
     !
     TYPE(HALO2LIST_ll), POINTER :: TPHALO2LIST ! list of HALO2_lls
     !
   END SUBROUTINE DEL_HALO2_ll
   !
END INTERFACE
!
END MODULE MODI_GET_HALO
!
!     ###################################################
      SUBROUTINE GET_HALO2(PSRC, TP_PSRC_HALO2_ll, HNAME)
!     ###################################################
!
USE MODE_ll
USE MODD_ARGSLIST_ll, ONLY : LIST_ll, HALO2LIST_ll
!
IMPLICIT NONE
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PSRC    ! variable at t
TYPE(HALO2LIST_ll), POINTER         :: TP_PSRC_HALO2_ll          ! halo2 for SRC
character(len=*), optional, intent(in) :: HNAME ! Name of the field to be added
!
character(len=:), allocatable    :: yname
INTEGER                          :: IIU,IJU,IKU            ! domain sizes
TYPE(LIST_ll)     , POINTER      :: TZ_PSRC_ll               ! halo
INTEGER                          :: IERROR                 ! error return code 
!
IIU = SIZE(PSRC,1)
IJU = SIZE(PSRC,2)
IKU = SIZE(PSRC,3)

if ( present ( hname ) ) then
  yname = hname
else
  yname = 'PSRC'
end if

NULLIFY( TZ_PSRC_ll,TP_PSRC_HALO2_ll)
CALL INIT_HALO2_ll( TP_PSRC_HALO2_ll, 1, IIU, IJU, IKU, 'GET_HALO2::' // TRIM( yname ) )
!
CALL ADD3DFIELD_ll( TZ_PSRC_ll, PSRC, 'GET_HALO2::'//trim( yname ) )
CALL UPDATE_HALO_ll(TZ_PSRC_ll,IERROR)
CALL UPDATE_HALO2_ll(TZ_PSRC_ll,TP_PSRC_HALO2_ll,IERROR)
!
!   clean local halo list
!
CALL CLEANLIST_ll(TZ_PSRC_ll)
!
END SUBROUTINE GET_HALO2
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!     ######################################
      SUBROUTINE GET_HALO(PSRC, HDIR, HNAME)
!     ######################################
!
USE MODE_ll
USE MODD_ARGSLIST_ll, ONLY : LIST_ll
!
IMPLICIT NONE
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PSRC    ! variable at t
CHARACTER(len=4), OPTIONAL :: HDIR ! to send only halo on X or Y direction
character(len=*), optional, intent(in) :: HNAME ! Name of the field to be added
!
character(len=:), allocatable    :: yname
TYPE(LIST_ll)     , POINTER      :: TZ_PSRC_ll               ! halo
INTEGER                          :: IERROR                 ! error return code 
!
NULLIFY( TZ_PSRC_ll)

if ( present ( hname ) ) then
  yname = hname
else
  yname = 'PSRC'
end if

CALL ADD3DFIELD_ll( TZ_PSRC_ll, PSRC, 'GET_HALO::'//trim( yname ) )
CALL UPDATE_HALO_ll(TZ_PSRC_ll,IERROR, HDIR=HDIR )
CALL CLEANLIST_ll(TZ_PSRC_ll)
!
END SUBROUTINE GET_HALO
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!     #########################
      SUBROUTINE GET_HALO_PHY(D,PSRC,OONGPU)
!     #########################
!
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
USE MODE_ll
USE MODD_ARGSLIST_ll, ONLY : LIST_ll
!
IMPLICIT NONE
!
TYPE(DIMPHYEX_t),        INTENT(IN)   :: D
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)  :: PSRC    ! variable at t
LOGICAL , OPTIONAL :: OONGPU
!
TYPE(LIST_ll)     , POINTER      :: TZ_PSRC_ll               ! halo
INTEGER                          :: IERROR                 ! error return code 
!
NULLIFY( TZ_PSRC_ll)
!
CALL ADD3DFIELD_ll( TZ_PSRC_ll, PSRC, 'GET_HALO::PSRC' )
CALL UPDATE_HALO_ll(TZ_PSRC_ll,IERROR, OONGPU=OONGPU )
CALL CLEANLIST_ll(TZ_PSRC_ll)
!
END SUBROUTINE GET_HALO_PHY
!-----------------------------------------------------------------------
#ifdef MNH_OPENACC
MODULE MODD_HALO_D

  USE MODD_PRECISION,  ONLY: MNHINT64
  
  IMPLICIT NONE
  
  REAL, SAVE , ALLOCATABLE, DIMENSION(:,:,:)  :: ZNORTH_IN, ZSOUTH_IN, ZWEST_IN, ZEAST_IN
  REAL, SAVE , ALLOCATABLE, DIMENSION(:,:,:)  :: ZNORTH_OUT, ZSOUTH_OUT, ZWEST_OUT, ZEAST_OUT

  REAL, SAVE , ALLOCATABLE, DIMENSION(:,:,:)  :: ZNORTHC_IN, ZSOUTHC_IN, ZWESTC_IN, ZEASTC_IN
  REAL, SAVE , ALLOCATABLE, DIMENSION(:,:,:)  :: ZNORTHC_OUT, ZSOUTHC_OUT, ZWESTC_OUT, ZEASTC_OUT 

  REAL, SAVE , ALLOCATABLE, DIMENSION(:,:)  :: ZNORTH2_IN, ZSOUTH2_IN, ZWEST2_IN, ZEAST2_IN
  REAL, SAVE , ALLOCATABLE, DIMENSION(:,:)  :: ZNORTH2_OUT, ZSOUTH2_OUT, ZWEST2_OUT, ZEAST2_OUT

  REAL, SAVE , ALLOCATABLE, DIMENSION(:,:)  :: ZNORTH_2D_IN, ZSOUTH_2D_IN, ZWEST_2D_IN, ZEAST_2D_IN
  REAL, SAVE , ALLOCATABLE, DIMENSION(:,:)  :: ZNORTH_2D_OUT, ZSOUTH_2D_OUT, ZWEST_2D_OUT, ZEAST_2D_OUT

  REAL, SAVE , ALLOCATABLE, DIMENSION(:,:)  :: ZNORTHC_2D_IN, ZSOUTHC_2D_IN, ZWESTC_2D_IN, ZEASTC_2D_IN
  REAL, SAVE , ALLOCATABLE, DIMENSION(:,:)  :: ZNORTHC_2D_OUT, ZSOUTHC_2D_OUT, ZWESTC_2D_OUT, ZEASTC_2D_OUT

  REAL, SAVE , ALLOCATABLE, DIMENSION(:,:,:)  :: ZNORTH2F_IN, ZSOUTH2F_IN, ZWEST2F_IN, ZEAST2F_IN
  REAL, SAVE , ALLOCATABLE, DIMENSION(:,:,:)  :: ZNORTH2F_OUT, ZSOUTH2F_OUT, ZWEST2F_OUT, ZEAST2F_OUT  
     
  LOGICAL, SAVE                               :: GFIRST_GET_HALO_D = .TRUE.
  
  LOGICAL, SAVE     :: GFIRST_INIT_HALO_D = .TRUE.
  INTEGER, SAVE     :: IHALO_1  
  INTEGER, SAVE     :: NP_NORTH,NP_SOUTH,NP_WEST,NP_EAST
  INTEGER, SAVE     :: IHALO2,IHALO2_1

CONTAINS
  
  SUBROUTINE INIT_HALO_D()

    USE MODE_MNH_ZWORK, ONLY : IIU,IJU,IKU
    USE MODE_MNH_ZWORK, ONLY : IIB,IJB ,IIE,IJE 
    USE MODE_MNH_ZWORK, ONLY : GWEST , GEAST, GSOUTH , GNORTH
    USE MODD_CONF, ONLY      : NHALO
    USE MODE_MNH_ZWORK, ONLY : NPMAX_T1DFLAT_R

    USE MODD_VAR_ll, ONLY    : IP,NPROC,NP1,NP2

    USE MODE_HALO_MANAGED
    
    IMPLICIT NONE

    IF (GFIRST_INIT_HALO_D) THEN 
       !
       IHALO_1  = NHALO-1
       IHALO2   = MAX(2,NHALO)
       IHALO2_1 = IHALO2-1
       !
       !  Init HALO
       !
       ALLOCATE  ( ZSOUTH_IN ( IIB:IIE   , IJB:IJB+IHALO_1   , IKU ) )
       ALLOCATE  ( ZNORTH_IN ( IIB:IIE   , IJE-IHALO_1:IJE   , IKU ) )
       ALLOCATE  ( ZWEST_IN  ( IIB:IIB+IHALO_1   , IJB:IJE   , IKU ) )
       ALLOCATE  ( ZEAST_IN  ( IIE-IHALO_1:IIE   , IJB:IJE   , IKU ) )
       !$acc enter data create (ZNORTH_IN, ZSOUTH_IN, ZWEST_IN, ZEAST_IN)
       !
       ALLOCATE  ( ZSOUTH_OUT (   IIB:IIE   , 1:IJB-1 , IKU ) )
       ALLOCATE  ( ZNORTH_OUT (   IIB:IIE   , IJE+1:IJU , IKU ) )
       ALLOCATE  ( ZWEST_OUT  ( 1:IIB-1 ,   IJB:IJE   , IKU ) )
       ALLOCATE  ( ZEAST_OUT  ( IIE+1:IIU ,   IJB:IJE   , IKU ) )
       !$acc enter data create (ZNORTH_OUT, ZSOUTH_OUT, ZWEST_OUT, ZEAST_OUT)
       !
       !  Init HALO with Corner
       !
       ALLOCATE  ( ZSOUTHC_IN (   1:IIU   , IJB:IJB+IHALO_1   , IKU ) )
       ALLOCATE  ( ZNORTHC_IN (   1:IIU   , IJE-IHALO_1:IJE   , IKU ) )
       ALLOCATE  ( ZWESTC_IN  ( IIB:IIB+IHALO_1   , IJB:IJE   , IKU ) )
       ALLOCATE  ( ZEASTC_IN  ( IIE-IHALO_1:IIE   , IJB:IJE   , IKU ) )
       !$acc enter data create (ZNORTHC_IN, ZSOUTHC_IN, ZWESTC_IN, ZEASTC_IN)
       !
       ALLOCATE  ( ZSOUTHC_OUT (   1:IIU   , 1:IJB-1 , IKU ) )
       ALLOCATE  ( ZNORTHC_OUT (   1:IIU   , IJE+1:IJU , IKU ) )
       ALLOCATE  ( ZWESTC_OUT  ( 1:IIB-1 ,   IJB:IJE   , IKU ) )
       ALLOCATE  ( ZEASTC_OUT  ( IIE+1:IIU ,   IJB:IJE   , IKU ) )
       !$acc enter data create (ZNORTHC_OUT, ZSOUTHC_OUT, ZWESTC_OUT, ZEASTC_OUT)       
       !
       !  Init HALO2
       !
       ALLOCATE  ( ZSOUTH2_IN ( IIU , IKU ) )
       ALLOCATE  ( ZNORTH2_IN ( IIU , IKU ) )
       ALLOCATE  ( ZWEST2_IN  ( IJU , IKU ) )
       ALLOCATE  ( ZEAST2_IN  ( IJU , IKU ) )
       !$acc enter data create (ZNORTH2_IN, ZSOUTH2_IN, ZWEST2_IN, ZEAST2_IN)
       !
       ALLOCATE  ( ZSOUTH2_OUT ( IIU , IKU ) )
       ALLOCATE  ( ZNORTH2_OUT ( IIU , IKU ) )
       ALLOCATE  ( ZWEST2_OUT  ( IJU , IKU ) )
       ALLOCATE  ( ZEAST2_OUT  ( IJU , IKU ) )
       !$acc enter data create (ZNORTH2_OUT, ZSOUTH2_OUT, ZWEST2_OUT, ZEAST2_OUT)
       !
       !  Init HALO_2D
       !
       ALLOCATE  ( ZSOUTH_2D_IN ( IIB:IIE   , IJB:IJB+IHALO_1  ) )
       ALLOCATE  ( ZNORTH_2D_IN ( IIB:IIE   , IJE-IHALO_1:IJE  ) )
       ALLOCATE  ( ZWEST_2D_IN  ( IIB:IIB+IHALO_1   , IJB:IJE  ) )
       ALLOCATE  ( ZEAST_2D_IN  ( IIE-IHALO_1:IIE   , IJB:IJE  ) )
       !$acc enter data create (ZNORTH_2D_IN, ZSOUTH_2D_IN, ZWEST_2D_IN, ZEAST_2D_IN)
       !
       ALLOCATE  ( ZSOUTH_2D_OUT (   IIB:IIE   , 1:IJB-1) )
       ALLOCATE  ( ZNORTH_2D_OUT (   IIB:IIE   , IJE+1:IJU) )
       ALLOCATE  ( ZWEST_2D_OUT  ( 1:IIB-1 ,   IJB:IJE  ) )
       ALLOCATE  ( ZEAST_2D_OUT  ( IIE+1:IIU ,   IJB:IJE  ) )
       !$acc enter data create (ZNORTH_2D_OUT, ZSOUTH_2D_OUT, ZWEST_2D_OUT, ZEAST_2D_OUT)
       !
       !  Init HALO 2D with Corner
       !
       ALLOCATE  ( ZSOUTHC_2D_IN (   1:IIU   , IJB:IJB+IHALO_1   ) )
       ALLOCATE  ( ZNORTHC_2D_IN (   1:IIU   , IJE-IHALO_1:IJE   ) )
       ALLOCATE  ( ZWESTC_2D_IN  ( IIB:IIB+IHALO_1   , IJB:IJE   ) )
       ALLOCATE  ( ZEASTC_2D_IN  ( IIE-IHALO_1:IIE   , IJB:IJE   ) )
       !$acc enter data create (ZNORTHC_2D_IN, ZSOUTHC_2D_IN, ZWESTC_2D_IN, ZEASTC_2D_IN)
       !
       ALLOCATE  ( ZSOUTHC_2D_OUT (   1:IIU   , 1:IJB-1 ) )
       ALLOCATE  ( ZNORTHC_2D_OUT (   1:IIU   , IJE+1:IJU ) )
       ALLOCATE  ( ZWESTC_2D_OUT  ( 1:IIB-1 ,   IJB:IJE   ) )
       ALLOCATE  ( ZEASTC_2D_OUT  ( IIE+1:IIU ,   IJB:IJE   ) )
       !$acc enter data create (ZNORTHC_2D_OUT, ZSOUTHC_2D_OUT, ZWESTC_2D_OUT, ZEASTC_2D_OUT)
       !
       !  Init HALO2 for Full update in 1 time <-> GET_HALO + GET_HALO2 
       !
       ALLOCATE  ( ZSOUTH2F_IN ( IIB:IIE   , IJB:IJB+IHALO2_1   , IKU ) )
       ALLOCATE  ( ZNORTH2F_IN ( IIB:IIE   , IJE-IHALO2_1:IJE   , IKU ) )
       ALLOCATE  ( ZWEST2F_IN  ( IIB:IIB+IHALO2_1   , IJB:IJE   , IKU ) )
       ALLOCATE  ( ZEAST2F_IN  ( IIE-IHALO2_1:IIE   , IJB:IJE   , IKU ) )
       !$acc enter data create (ZNORTH2F_IN, ZSOUTH2F_IN, ZWEST2F_IN, ZEAST2F_IN)
       !
       ALLOCATE  ( ZSOUTH2F_OUT (   IIB:IIE   , IJB-IHALO2:IJB-1 , IKU ) )
       ALLOCATE  ( ZNORTH2F_OUT (   IIB:IIE   , IJE+1:IJE+IHALO2 , IKU ) )
       ALLOCATE  ( ZWEST2F_OUT  ( IIB-IHALO2:IIB-1 ,   IJB:IJE   , IKU ) )
       ALLOCATE  ( ZEAST2F_OUT  ( IIE+1:IIE+IHALO2 ,   IJB:IJE   , IKU ) )
       !$acc enter data create (ZNORTH2F_OUT, ZSOUTH2F_OUT, ZWEST2F_OUT, ZEAST2F_OUT)

       CALL MNH_HALO_MANAGED_INIT()

       IF (.NOT. GWEST ) THEN
          NP_WEST = ( IP-1 -1 ) + 1
       ELSE
          NP_WEST = 0
       ENDIF
       IF (.NOT. GEAST ) THEN
          NP_EAST = ( IP-1 +1 ) + 1
       ELSE
          NP_EAST = 0
       ENDIF
       IF (.NOT. GSOUTH ) THEN
          NP_SOUTH =  ( IP-1 -NP1 ) + 1
       ELSE
          NP_SOUTH = 0
       ENDIF
       IF (.NOT. GNORTH ) THEN
          NP_NORTH =  ( IP-1 +NP1 ) + 1
       ELSE
          NP_NORTH = 0
       ENDIF

       !print*,"PROC=",IP, GWEST,NP_WEST, GEAST,NP_EAST, GSOUTH,NP_SOUTH ,  GNORTH,NP_NORTH
       
       GFIRST_INIT_HALO_D = .FALSE.
       
    END IF

  END SUBROUTINE INIT_HALO_D

END MODULE MODD_HALO_D

!     #########################
      SUBROUTINE GET_HALO_D(PSRC,HDIR,HNAME)
!     #########################
!define MNH_GPUDIRECT
!
USE MODD_HALO_D

!USE MODE_MNH_ZWORK, ONLY : GWEST , GEAST, GSOUTH , GNORTH
USE MODE_MNH_ZWORK, ONLY : IIU,IJU,IKU
USE MODE_MNH_ZWORK, ONLY : IIB,IJB ,IIE,IJE

!!
USE MODE_MPPDB
USE MODI_GET_HALO, ONLY : GET_HALO_START_D,GET_HALO_STOP_D
!
IMPLICIT NONE
!
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PSRC    ! variable at t
! acc declare present (PSRC)
CHARACTER(len=4), OPTIONAL :: HDIR ! to send only halo on X or Y direction
character(len=*), optional, intent(in) :: HNAME ! Name of the field to be added
!
INTEGER            :: INB_REQ , IREQ(8)
!
CALL MPPDB_CHECK(PSRC,"GET_HALO_D big:PSRC")
!
CALL GET_HALO_START_D(PSRC,INB_REQ,IREQ,&
     ZNORTH_IN , ZSOUTH_IN , ZWEST_IN , ZEAST_IN , &
     ZNORTH_OUT, ZSOUTH_OUT, ZWEST_OUT, ZEAST_OUT, &
     IIB,IIE,IJB,IJE,IIU,IJU,IKU,IHALO_1,HDIR)

CALL GET_HALO_STOP_D(PSRC,INB_REQ,IREQ,&
     ZNORTH_OUT, ZSOUTH_OUT, ZWEST_OUT, ZEAST_OUT, &
     IIB,IIE,IJB,IJE,IIU,IJU,IKU,HDIR)
! !   CALL MPPDB_CHECK(PSRC,"UPDATE_HALO_ll::"//Trim(HNAME))
!   CALL MPPDB_CHECK(PSRC,Trim(HNAME))
!
CALL MPPDB_CHECK(PSRC,"GET_HALO_D end:PSRC")
!
END SUBROUTINE GET_HALO_D
!     #########################
SUBROUTINE GET_HALO_START_D(PSRC,KNB_REQ,KREQ,&
     PZNORTH_IN , PZSOUTH_IN , PZWEST_IN , PZEAST_IN , &
     PZNORTH_OUT, PZSOUTH_OUT, PZWEST_OUT, PZEAST_OUT, &
     KIIB,KIIE,KIJB,KIJE,KIIU,KIJU,KIKU,KIHALO_1,&
     HDIR)
!     #########################
!define MNH_GPUDIRECT
!
USE MODD_HALO_D

USE MODE_MNH_ZWORK, ONLY : GWEST , GEAST, GSOUTH , GNORTH
!
USE MODD_VAR_ll, ONLY : NMNH_COMM_WORLD
USE MODD_VAR_ll,    ONLY : MNH_STATUSES_IGNORE
USE MODD_PRECISION,  ONLY : MNHREAL_MPI
!
USE MODE_MPPDB
!
IMPLICIT NONE
!
REAL, DIMENSION(KIIU,KIJU,KIKU), INTENT(INOUT) :: PSRC    ! variable at t
! acc declare present (PSRC)
INTEGER                               :: KNB_REQ ,  KREQ(8)
REAL            :: PZSOUTH_IN ( KIIB:KIIE   , KIJB:KIJB+KIHALO_1   , KIKU ) ,&
                   PZNORTH_IN ( KIIB:KIIE   , KIJE-KIHALO_1:KIJE   , KIKU ) ,&
                   PZWEST_IN  ( KIIB:KIIB+KIHALO_1   , KIJB:KIJE   , KIKU ) ,&
                   PZEAST_IN  ( KIIE-KIHALO_1:KIIE   , KIJB:KIJE   , KIKU ) ,&
       !
                   PZSOUTH_OUT (   KIIB:KIIE   , 1:KIJB-1 , KIKU ) ,&
                   PZNORTH_OUT (   KIIB:KIIE   , KIJE+1:KIJU , KIKU ) ,&
                   PZWEST_OUT  ( 1:KIIB-1 ,   KIJB:KIJE   , KIKU ) ,&
                   PZEAST_OUT  ( KIIE+1:KIIU ,   KIJB:KIJE   , KIKU ) 
INTEGER         :: KIIB,KIIE,KIJB,KIJE,KIIU,KIJU,KIKU,KIHALO_1

CHARACTER(len=4), OPTIONAL :: HDIR ! to send only halo on X or Y direction
!
INTEGER                          :: IERROR                 ! error return code 

INTEGER,PARAMETER :: IS_WEST=1 , IS_EAST=2, IS_SOUTH=3, IS_NORTH=4
LOGICAL      :: LX , LY
INTEGER      :: NB_REQ, IERR
!
INTEGER :: JI,JJ,JK 

CALL INIT_HALO_D()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#if 0
!$acc data present (PSRC) &
!$acc present (PZNORTH_IN, PZSOUTH_IN, PZWEST_IN, PZEAST_IN) &
!$acc present (PZNORTH_OUT, PZSOUTH_OUT, PZWEST_OUT, PZEAST_OUT) 
#else
!PW: note: present for PZ* arrays removed to prevent crash with Cray compiler (CCE 15.0.1)
!$acc data present( PSRC )
#endif

LX = .FALSE.
LY = .FALSE. 

IF (.NOT. PRESENT(HDIR) ) THEN
LX = .TRUE.
LY = .TRUE.
ELSE
!!$LX = ( HDIR == "01_X"  .OR. HDIR == "S0_X" )
!!$LY = ( HDIR == "01_Y"  .OR. HDIR == "S0_Y" )
LX = ( HDIR == "01_X" .OR. HDIR == "S0_X" .OR. HDIR == "S0_Y" )
LY = ( HDIR == "01_Y" .OR. HDIR == "S0_Y" .OR. HDIR == "S0_X" )
!!$print *,"KIIB=",KIIB," HDIR=",HDIR," LX=",LX," LY=",LY ; call flush(6)
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

NB_REQ = 0

!
! Post the recieve of Zxxxx_IN buffer first via MPI(Gpu_direct)
!

IF (LX) THEN
   IF (.NOT. GWEST) THEN
#ifdef MNH_GPUDIRECT
      !$acc host_data use_device(PZWEST_OUT)
#endif
      NB_REQ = NB_REQ + 1
      CALL MPI_IRECV(PZWEST_OUT,SIZE(PZWEST_OUT),MNHREAL_MPI,NP_WEST-1,1000+IS_EAST,&
                     NMNH_COMM_WORLD,KREQ(NB_REQ),IERR)
#ifdef MNH_GPUDIRECT
      !$acc end host_data
#endif
   END IF
   IF (.NOT.GEAST) THEN 
#ifdef MNH_GPUDIRECT
      !$acc host_data use_device(PZEAST_OUT)
#endif
      NB_REQ = NB_REQ + 1
      CALL MPI_IRECV(PZEAST_OUT,SIZE(PZEAST_OUT),MNHREAL_MPI,NP_EAST-1,1000+IS_WEST,&
                     NMNH_COMM_WORLD,KREQ(NB_REQ),IERR)
#ifdef MNH_GPUDIRECT
      !$acc end host_data
#endif
   ENDIF
END IF

IF (LY) THEN
   IF (.NOT.GSOUTH) THEN
#ifdef MNH_GPUDIRECT
      !$acc host_data use_device(PZSOUTH_OUT)
#endif
      NB_REQ = NB_REQ + 1
      CALL MPI_IRECV(PZSOUTH_OUT,SIZE(PZSOUTH_OUT),MNHREAL_MPI,NP_SOUTH-1,1000+IS_NORTH,&
                     NMNH_COMM_WORLD,KREQ(NB_REQ),IERR)
#ifdef MNH_GPUDIRECT
      !$acc end host_data
#endif
   ENDIF
   IF (.NOT.GNORTH) THEN
#ifdef MNH_GPUDIRECT
      !$acc host_data use_device(PZNORTH_OUT)
#endif
      NB_REQ = NB_REQ + 1
      CALL MPI_IRECV(PZNORTH_OUT,SIZE(PZNORTH_OUT),MNHREAL_MPI,NP_NORTH-1,1000+IS_SOUTH,&
                     NMNH_COMM_WORLD,KREQ(NB_REQ),IERR)
#ifdef MNH_GPUDIRECT
      !$acc end host_data
#endif
   ENDIF
END IF

!
! Copy the halo(async) on the device PSRC to Zxxxx_IN buffer 
!

IF (LX) THEN
   IF (.NOT. GWEST) THEN
      !$acc kernels async(IS_WEST)
      !$mnh_expand_array(JI=KIIB:KIIB+KIHALO_1 , JJ=KIJB:KIJE , JK=1:KIKU )
           PZWEST_IN ( KIIB:KIIB+KIHALO_1  ,    KIJB:KIJE  , : )  = PSRC( KIIB:KIIB+KIHALO_1  ,  KIJB:KIJE  , : )
      !$mnh_end_expand_array()
      !$acc end kernels
   END IF
   IF (.NOT.GEAST) THEN
      !$acc kernels async(IS_EAST)
      !$mnh_expand_array(JI=KIIE-KIHALO_1:KIIE , JJ=KIJB:KIJE , JK=1:KIKU)
           PZEAST_IN ( KIIE-KIHALO_1:KIIE  ,    KIJB:KIJE  , : )  = PSRC( KIIE-KIHALO_1:KIIE  ,  KIJB:KIJE  , : )
      !$mnh_end_expand_array()
      !$acc end kernels
   ENDIF
END IF

IF (LY) THEN
   IF (.NOT.GSOUTH) THEN
      !$acc kernels async(IS_SOUTH)
      !$mnh_expand_array(JI=KIIB:KIIE , JJ=KIJB:KIJB+KIHALO_1 , JK=1:KIKU )
           PZSOUTH_IN ( KIIB:KIIE  ,    KIJB:KIJB+KIHALO_1  , : ) = PSRC( KIIB:KIIE  ,    KIJB:KIJB+KIHALO_1  , : )
      !$mnh_end_expand_array()
      !$acc end kernels
   ENDIF
   IF (.NOT.GNORTH) THEN
      !$acc kernels async(IS_NORTH)
      !$mnh_expand_array(JI=KIIB:KIIE , JJ=KIJE-KIHALO_1:KIJE , JK=1:KIKU )
           PZNORTH_IN ( KIIB:KIIE  ,    KIJE-KIHALO_1:KIJE  , : ) = PSRC( KIIB:KIIE  ,    KIJE-KIHALO_1:KIJE  , : )
      !$mnh_end_expand_array()
      !$acc end kernels
   ENDIF   
ENDIF

!$acc wait

!
! Send  Zxxxx_IN buffer via MPI(Gpu_direct)
!
IF (LX) THEN
   IF (.NOT. GWEST) THEN
#ifdef MNH_GPUDIRECT
      !$acc host_data use_device(PZWEST_IN)
#else
      !$acc update host(PZWEST_IN)
#endif
      NB_REQ = NB_REQ + 1
      CALL MPI_ISEND(PZWEST_IN,SIZE(PZWEST_IN)  ,MNHREAL_MPI,NP_WEST-1,1000+IS_WEST,&
                     NMNH_COMM_WORLD,KREQ(NB_REQ),IERR)
#ifdef MNH_GPUDIRECT
      !$acc end host_data
#endif
   END IF
   IF (.NOT.GEAST) THEN
#ifdef MNH_GPUDIRECT
      !$acc host_data use_device(PZEAST_IN)
#else
      !$acc update host(PZEAST_IN)
#endif
      NB_REQ = NB_REQ + 1
      CALL MPI_ISEND(PZEAST_IN,SIZE(PZEAST_IN)  ,MNHREAL_MPI,NP_EAST-1,1000+IS_EAST,&
                     NMNH_COMM_WORLD,KREQ(NB_REQ),IERR)
#ifdef MNH_GPUDIRECT
      !$acc end host_data
#endif
   ENDIF
END IF

IF (LY) THEN
   IF (.NOT.GSOUTH) THEN
#ifdef MNH_GPUDIRECT
      !$acc host_data use_device(PZSOUTH_IN)
#else
      !$acc update host(PZSOUTH_IN)
#endif
      NB_REQ = NB_REQ + 1
      CALL MPI_ISEND(PZSOUTH_IN,SIZE(PZSOUTH_IN)  ,MNHREAL_MPI,NP_SOUTH-1,1000+IS_SOUTH,&
                     NMNH_COMM_WORLD,KREQ(NB_REQ),IERR)
#ifdef MNH_GPUDIRECT
      !$acc end host_data
#endif
   ENDIF
   IF (.NOT.GNORTH) THEN
#ifdef MNH_GPUDIRECT
      !$acc host_data use_device(PZNORTH_IN)
#else
      !$acc update host(PZNORTH_IN)
#endif
      NB_REQ = NB_REQ + 1
      CALL MPI_ISEND(PZNORTH_IN,SIZE(PZNORTH_IN)  ,MNHREAL_MPI,NP_NORTH-1,1000+IS_NORTH,&
                     NMNH_COMM_WORLD,KREQ(NB_REQ),IERR)
#ifdef MNH_GPUDIRECT
      !$acc end host_data
#endif
   ENDIF   
ENDIF

!$acc end data

KNB_REQ = NB_REQ
!
END SUBROUTINE GET_HALO_START_D
!
!     #########################
SUBROUTINE GET_HALO_STOP_D(PSRC,KNB_REQ,KREQ,&
     PZNORTH_OUT, PZSOUTH_OUT, PZWEST_OUT, PZEAST_OUT, &
     KIIB,KIIE,KIJB,KIJE,KIIU,KIJU,KIKU,&
     HDIR)
!     #########################
!define MNH_GPUDIRECT
!
USE MODD_HALO_D

USE MODE_MNH_ZWORK, ONLY : GWEST , GEAST, GSOUTH , GNORTH
!!$USE MODE_MNH_ZWORK, ONLY : IIU,IJU,IKU
!!$USE MODE_MNH_ZWORK, ONLY : IIB,IJB ,IIE,IJE 
!
USE MODD_VAR_ll, ONLY : NMNH_COMM_WORLD
USE MODD_VAR_ll,    ONLY : MNH_STATUSES_IGNORE
USE MODD_PRECISION,  ONLY : MNHREAL_MPI
!
USE MODE_MPPDB
!
IMPLICIT NONE
!
REAL, DIMENSION(KIIU,KIJU,KIKU), INTENT(INOUT) :: PSRC    ! variable at t
! acc declare present (PSRC)
INTEGER                               :: KNB_REQ , KREQ(8)
REAL            :: PZSOUTH_OUT (   KIIB:KIIE   , 1:KIJB-1 , KIKU ) ,&
                   PZNORTH_OUT (   KIIB:KIIE   , KIJE+1:KIJU , KIKU ) ,&
                   PZWEST_OUT  ( 1:KIIB-1 ,   KIJB:KIJE   , KIKU ) ,&
                   PZEAST_OUT  ( KIIE+1:KIIU ,   KIJB:KIJE   , KIKU ) 
INTEGER         :: KIIB,KIIE,KIJB,KIJE,KIIU,KIJU,KIKU

CHARACTER(len=4), OPTIONAL :: HDIR ! to send only halo on X or Y direction
!
INTEGER                          :: IERROR                 ! error return code 

INTEGER,PARAMETER :: IS_WEST=1 , IS_EAST=2, IS_SOUTH=3, IS_NORTH=4
LOGICAL      :: LX , LY
INTEGER      :: NB_REQ, IERR
!
INTEGER :: JI,JJ,JK

CALL INIT_HALO_D()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!$acc data  present (PSRC) &
!$acc & present (PZNORTH_OUT,PZSOUTH_OUT,PZWEST_OUT,PZEAST_OUT)

LX = .FALSE.
LY = .FALSE. 

IF (.NOT. PRESENT(HDIR) ) THEN
LX = .TRUE.
LY = .TRUE.
ELSE
!!$LX = ( HDIR == "01_X"  .OR. HDIR == "S0_X" )
!!$LY = ( HDIR == "01_Y"  .OR. HDIR == "S0_Y" )
LX = ( HDIR == "01_X" .OR. HDIR == "S0_X" .OR. HDIR == "S0_Y" )
LY = ( HDIR == "01_Y" .OR. HDIR == "S0_Y" .OR. HDIR == "S0_X" )
!!$print *,"KIIB=",KIIB," HDIR=",HDIR," LX=",LX," LY=",LY ; call flush(6)
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

NB_REQ = KNB_REQ

CALL MPI_WAITALL(NB_REQ,KREQ,MNH_STATUSES_IGNORE,IERR)

!
!   Copy back the Zxxx_OUT buffer recv via MPI(gpu_direct) to PSRC halo
! 

IF (LX) THEN
   IF (.NOT.GWEST) THEN
#ifndef MNH_GPUDIRECT
   !$acc update device(PZWEST_OUT) async(IS_WEST)
#endif
   !$acc kernels async(IS_WEST)
   !$mnh_expand_array(JI=1:KIIB-1 , JJ=KIJB:KIJE , JK=1:KIKU )
        PSRC( 1:KIIB-1  ,      KIJB:KIJE      , : ) = PZWEST_OUT( 1:KIIB-1  ,   KIJB:KIJE    , : )
   !$mnh_end_expand_array()
   !$acc end kernels
   ENDIF
   IF (.NOT.GEAST) THEN
#ifndef MNH_GPUDIRECT
   !$acc update device(PZEAST_OUT) async(IS_EAST)
#endif
   !$acc kernels async(IS_EAST)
   !$mnh_expand_array(JI=KIIE+1:KIIU , JJ=KIJB:KIJE , JK=1:KIKU )
        PSRC( KIIE+1:KIIU  ,      KIJB:KIJE      , : ) = PZEAST_OUT( KIIE+1:KIIU  ,   KIJB:KIJE    , : )  
   !$mnh_end_expand_array()
   !$acc end kernels
   ENDIF
END IF
IF (LY) THEN
   IF (.NOT.GSOUTH) THEN
#ifndef MNH_GPUDIRECT
   !$acc update device(PZSOUTH_OUT) async(IS_SOUTH)
#endif
   !$acc kernels async(IS_SOUTH)
   !$mnh_expand_array(JI=KIIB:KIIE , JJ=1:KIJB-1 , JK=1:KIKU )
        PSRC(      KIIB:KIIE       ,  1:KIJB-1 , : ) = PZSOUTH_OUT(  KIIB:KIIE     , 1:KIJB-1  , : )
   !$mnh_end_expand_array()
   !$acc end kernels
   ENDIF
   IF (.NOT.GNORTH) THEN
#ifndef MNH_GPUDIRECT
   !$acc update device(PZNORTH_OUT) async(IS_NORTH)
#endif
   !$acc kernels async(IS_NORTH)
   !$mnh_expand_array(JI=KIIB:KIIE , JJ=KIJE+1:KIJU , JK=1:KIKU )
        PSRC(      KIIB:KIIE       , KIJE+1:KIJU , : ) = PZNORTH_OUT (  KIIB:KIIE     , KIJE+1:KIJU  , : )
   !$mnh_end_expand_array()
   !$acc end kernels
   ENDIF
END IF
!$acc wait

!$acc end data
!
END SUBROUTINE GET_HALO_STOP_D
!-------------------------------------------------------------------------------
!     ########################################
      SUBROUTINE GET_HALO_DD(PSRC, HDIR, HNAME)
!     ########################################
!define MNH_GPUDIRECT
!
USE MODD_HALO_D
USE MODE_ll
USE MODD_ARGSLIST_ll, ONLY : LIST_ll
USE MODD_PARAMETERS, ONLY : JPHEXT
!
USE MODD_IO,        ONLY : GSMONOPROC
USE MODE_MNH_ZWORK, ONLY : GWEST , GEAST, GSOUTH , GNORTH
!
USE MODD_CONF, ONLY : NHALO
USE MODE_MPPDB

USE MODD_VAR_ll, ONLY    : IP,NPROC,NP1,NP2
USE MODD_VAR_ll, ONLY : NMNH_COMM_WORLD
USE MODD_VAR_ll,    ONLY : MNH_STATUSES_IGNORE
USE MODD_PRECISION,  ONLY : MNHREAL_MPI
!
USE MODE_MNH_ZWORK, ONLY : IIU,IJU,IKU
USE MODE_MNH_ZWORK, ONLY : IIB,IJB ,IIE,IJE 
!
IMPLICIT NONE
!
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PSRC    ! variable at t
CHARACTER(len=4), OPTIONAL :: HDIR ! to send only halo on X or Y direction
character(len=*), optional, intent(in) :: HNAME ! Name of the field to be added
!
character(len=:), allocatable    :: yname
TYPE(LIST_ll)     , POINTER      :: TZ_PSRC_ll               ! halo
INTEGER                          :: IERROR                 ! error return code 

INTEGER,PARAMETER :: IS_WEST=1 , IS_EAST=2, IS_SOUTH=3, IS_NORTH=4

LOGICAL      :: LX , LY

INTEGER      :: INB_REQ , IREQ(8)
INTEGER      :: IERR

CALL MPPDB_CHECK(PSRC,"GET_HALO_DD big:PSRC")

if ( NPROC == 1 ) THEN
   CALL MPPDB_CHECK(PSRC,"GET_HALO_DD end:PSRC")
   RETURN
end if

CALL INIT_HALO_D()

!$acc data present ( PSRC )

NULLIFY( TZ_PSRC_ll)
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

LX = .FALSE.
LY = .FALSE. 

IF (.NOT. PRESENT(HDIR) ) THEN
LX = .TRUE.
LY = .TRUE.
ELSE
   !
   !  Problem of reproductibility in ppm_s0_x/y if only S0_X or S0_Y
   !  so add S0_X + S0_Y for ppm_s0*
   !
!!$LX = ( HDIR == "01_X" .OR. HDIR == "S0_X" )
!!$LY = ( HDIR == "01_Y" .OR. HDIR == "S0_Y" )
LX = ( HDIR == "01_X" .OR. HDIR == "S0_X" .OR. HDIR == "S0_Y" )
LY = ( HDIR == "01_Y" .OR. HDIR == "S0_Y" .OR. HDIR == "S0_X" )
END IF

!!$LX = .TRUE.
!!$LY = .TRUE.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

INB_REQ = 0

!
! Post the recieve of Zxxxx_IN buffer first via MPI(Gpu_direct)
!

IF (LX) THEN
   IF (.NOT. GWEST) THEN
#ifdef MNH_GPUDIRECT
      !$acc host_data use_device(ZWEST_OUT)
#endif
      INB_REQ = INB_REQ + 1
      CALL MPI_IRECV(ZWEST_OUT,SIZE(ZWEST_OUT),MNHREAL_MPI,NP_WEST-1,1000+IS_EAST,&
                     NMNH_COMM_WORLD,IREQ(INB_REQ),IERR)
#ifdef MNH_GPUDIRECT
      !$acc end host_data
#endif
   END IF
   IF (.NOT.GEAST) THEN 
#ifdef MNH_GPUDIRECT
      !$acc host_data use_device(ZEAST_OUT)
#endif
      INB_REQ = INB_REQ + 1
      CALL MPI_IRECV(ZEAST_OUT,SIZE(ZEAST_OUT),MNHREAL_MPI,NP_EAST-1,1000+IS_WEST,&
                     NMNH_COMM_WORLD,IREQ(INB_REQ),IERR)
#ifdef MNH_GPUDIRECT
      !$acc end host_data
#endif
   ENDIF
END IF

IF (LY) THEN
   IF (.NOT.GSOUTH) THEN
#ifdef MNH_GPUDIRECT
      !$acc host_data use_device(ZSOUTH_OUT)
#endif
      INB_REQ = INB_REQ + 1
      CALL MPI_IRECV(ZSOUTH_OUT,SIZE(ZSOUTH_OUT),MNHREAL_MPI,NP_SOUTH-1,1000+IS_NORTH,&
                     NMNH_COMM_WORLD,IREQ(INB_REQ),IERR)
#ifdef MNH_GPUDIRECT
      !$acc end host_data
#endif
   ENDIF
   IF (.NOT.GNORTH) THEN
#ifdef MNH_GPUDIRECT
      !$acc host_data use_device(ZNORTH_OUT)
#endif
      INB_REQ = INB_REQ + 1
      CALL MPI_IRECV(ZNORTH_OUT,SIZE(ZNORTH_OUT),MNHREAL_MPI,NP_NORTH-1,1000+IS_SOUTH,&
                     NMNH_COMM_WORLD,IREQ(INB_REQ),IERR)
#ifdef MNH_GPUDIRECT
      !$acc end host_data
#endif
   ENDIF
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Copy the halo on the device PSRC to Zxxxx_IN

IF (LX) THEN
   IF (.NOT. GWEST) THEN
   !$acc kernels async(IS_WEST)
   ZWEST_IN ( IIB:IIB+IHALO_1  ,    IJB:IJE  , : )  = PSRC( IIB:IIB+IHALO_1  ,  IJB:IJE  , : )
   !$acc end kernels
      END IF
   IF (.NOT.GEAST) THEN
   !$acc kernels async(IS_EAST)
   ZEAST_IN ( IIE-IHALO_1:IIE  ,    IJB:IJE  , : )  = PSRC( IIE-IHALO_1:IIE  ,  IJB:IJE  , : )
   !$acc end kernels
      ENDIF
END IF
IF (LY) THEN
   IF (.NOT.GSOUTH) THEN
   !$acc kernels async(IS_SOUTH)
   ZSOUTH_IN ( IIB:IIE  ,    IJB:IJB+IHALO_1  , : ) = PSRC( IIB:IIE  ,    IJB:IJB+IHALO_1  , : )
   !$acc end kernels
      ENDIF
   IF (.NOT.GNORTH) THEN
   !$acc kernels async(IS_NORTH)
   ZNORTH_IN ( IIB:IIE  ,    IJE-IHALO_1:IJE  , : ) = PSRC( IIB:IIE  ,    IJE-IHALO_1:IJE  , : )
   !$acc end kernels
   ENDIF
ENDIF
!$acc wait

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Send  Zxxxx_IN buffer via MPI(Gpu_direct) or copy to host
!
IF (LX) THEN
   IF (.NOT. GWEST) THEN
#ifdef MNH_GPUDIRECT
      !$acc host_data use_device(ZWEST_IN)
#else
      !$acc update host(ZWEST_IN)
#endif
      INB_REQ = INB_REQ + 1
      CALL MPI_ISEND(ZWEST_IN,SIZE(ZWEST_IN)  ,MNHREAL_MPI,NP_WEST-1,1000+IS_WEST,&
                     NMNH_COMM_WORLD,IREQ(INB_REQ),IERR)
#ifdef MNH_GPUDIRECT
      !$acc end host_data
#endif
   END IF
   IF (.NOT.GEAST) THEN
#ifdef MNH_GPUDIRECT
      !$acc host_data use_device(ZEAST_IN)
#else
      !$acc update host(ZEAST_IN)
#endif
      INB_REQ = INB_REQ + 1
      CALL MPI_ISEND(ZEAST_IN,SIZE(ZEAST_IN)  ,MNHREAL_MPI,NP_EAST-1,1000+IS_EAST,&
                     NMNH_COMM_WORLD,IREQ(INB_REQ),IERR)
#ifdef MNH_GPUDIRECT
      !$acc end host_data
#endif
   ENDIF
END IF

IF (LY) THEN
   IF (.NOT.GSOUTH) THEN
#ifdef MNH_GPUDIRECT
      !$acc host_data use_device(ZSOUTH_IN)
#else
      !$acc update host(ZSOUTH_IN)
#endif
      INB_REQ = INB_REQ + 1
      CALL MPI_ISEND(ZSOUTH_IN,SIZE(ZSOUTH_IN)  ,MNHREAL_MPI,NP_SOUTH-1,1000+IS_SOUTH,&
                     NMNH_COMM_WORLD,IREQ(INB_REQ),IERR)
#ifdef MNH_GPUDIRECT
      !$acc end host_data
#endif
   ENDIF
   IF (.NOT.GNORTH) THEN
#ifdef MNH_GPUDIRECT
      !$acc host_data use_device(ZNORTH_IN)
#else
      !$acc update host(ZNORTH_IN)
#endif
      INB_REQ = INB_REQ + 1
      CALL MPI_ISEND(ZNORTH_IN,SIZE(ZNORTH_IN)  ,MNHREAL_MPI,NP_NORTH-1,1000+IS_NORTH,&
                     NMNH_COMM_WORLD,IREQ(INB_REQ),IERR)
#ifdef MNH_GPUDIRECT
      !$acc end host_data
#endif
   ENDIF   
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF ( INB_REQ > 0 ) THEN
   CALL MPI_WAITALL(INB_REQ,IREQ,MNH_STATUSES_IGNORE,IERR)
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Is update halo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


IF (LX) THEN
   IF (.NOT.GWEST) THEN
#ifndef MNH_GPUDIRECT
   !$acc update device(ZWEST_OUT) async(IS_WEST)
#endif
   !$acc kernels async(IS_WEST)
   PSRC( 1:IIB-1  ,      IJB:IJE      , : ) = ZWEST_OUT( 1:IIB-1  ,   IJB:IJE    , : )
   !$acc end kernels
   ENDIF
   IF (.NOT.GEAST) THEN
#ifndef MNH_GPUDIRECT
   !$acc update device(ZEAST_OUT) async(IS_EAST)
#endif
   !$acc kernels async(IS_EAST)
   PSRC( IIE+1:IIU  ,      IJB:IJE      , : ) = ZEAST_OUT( IIE+1:IIU  ,   IJB:IJE    , : )  
   !$acc end kernels
   ENDIF
END IF
IF (LY) THEN
   IF (.NOT.GSOUTH) THEN
#ifndef MNH_GPUDIRECT
   !$acc update device(ZSOUTH_OUT) async(IS_SOUTH)
#endif
   !$acc kernels async(IS_SOUTH)
   PSRC(      IIB:IIE       ,  1:IJB-1 , : ) = ZSOUTH_OUT(  IIB:IIE     , 1:IJB-1  , : )
   !$acc end kernels
   ENDIF
   IF (.NOT.GNORTH) THEN
#ifndef MNH_GPUDIRECT
   !$acc update device(ZNORTH_OUT) async(IS_NORTH)
#endif
   !$acc kernels async(IS_NORTH)
   PSRC(      IIB:IIE       , IJE+1:IJU , : ) = ZNORTH_OUT (  IIB:IIE     , IJE+1:IJU  , : )
   !$acc end kernels
   ENDIF
END IF
!$acc wait

!$acc end data

CALL MPPDB_CHECK(PSRC,"GET_HALO_DD end:PSRC")

END SUBROUTINE GET_HALO_DD
!     ########################################
      SUBROUTINE GET_2D_HALO_DD(PSRC, HDIR, HNAME)
!     ########################################
!define MNH_GPUDIRECT
!
USE MODD_HALO_D
USE MODE_ll
USE MODD_ARGSLIST_ll, ONLY : LIST_ll
USE MODD_PARAMETERS, ONLY : JPHEXT
!
USE MODD_IO,        ONLY : GSMONOPROC
USE MODE_MNH_ZWORK, ONLY : GWEST , GEAST, GSOUTH , GNORTH
!
USE MODD_CONF, ONLY : NHALO
USE MODE_MPPDB

USE MODD_VAR_ll, ONLY    : IP,NPROC,NP1,NP2
USE MODD_VAR_ll, ONLY : NMNH_COMM_WORLD
USE MODD_VAR_ll,    ONLY : MNH_STATUSES_IGNORE
USE MODD_PRECISION,  ONLY : MNHREAL_MPI
!
USE MODE_MNH_ZWORK, ONLY : IIU,IJU,IKU
USE MODE_MNH_ZWORK, ONLY : IIB,IJB ,IIE,IJE 
!
IMPLICIT NONE
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSRC    ! variable at t
CHARACTER(len=4), OPTIONAL :: HDIR ! to send only halo on X or Y direction
character(len=*), optional, intent(in) :: HNAME ! Name of the field to be added
!
character(len=:), allocatable    :: yname
TYPE(LIST_ll)     , POINTER      :: TZ_PSRC_ll               ! halo
INTEGER                          :: IERROR                 ! error return code 

INTEGER,PARAMETER :: IS_WEST=1 , IS_EAST=2, IS_SOUTH=3, IS_NORTH=4

LOGICAL      :: LX , LY

INTEGER      :: INB_REQ , IREQ(8)
INTEGER      :: IERR

CALL MPPDB_CHECK(PSRC,"GET_2D_HALO_DD big:PSRC")

if ( NPROC == 1 ) THEN
   CALL MPPDB_CHECK(PSRC,"GET_2D_HALO_DD end:PSRC")
   RETURN
end if

CALL INIT_HALO_D()

!$acc data present ( PSRC )

NULLIFY( TZ_PSRC_ll)
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

LX = .FALSE.
LY = .FALSE. 

IF (.NOT. PRESENT(HDIR) ) THEN
LX = .TRUE.
LY = .TRUE.
ELSE
   !
   !  Problem of reproductibility in ppm_s0_x/y if only S0_X or S0_Y
   !  so add S0_X + S0_Y for ppm_s0*
   !
!!$LX = ( HDIR == "01_X" .OR. HDIR == "S0_X" )
!!$LY = ( HDIR == "01_Y" .OR. HDIR == "S0_Y" )
LX = ( HDIR == "01_X" .OR. HDIR == "S0_X" .OR. HDIR == "S0_Y" )
LY = ( HDIR == "01_Y" .OR. HDIR == "S0_Y" .OR. HDIR == "S0_X" )
END IF

!!$LX = .TRUE.
!!$LY = .TRUE.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

INB_REQ = 0

!
! Post the recieve of Zxxxx_2D_OUT buffer first via MPI(Gpu_direct)
!

IF (LX) THEN
   IF (.NOT. GWEST) THEN
#ifdef MNH_GPUDIRECT
      !$acc host_data use_device(ZWEST_2D_OUT)
#endif
      INB_REQ = INB_REQ + 1
      CALL MPI_IRECV(ZWEST_2D_OUT,SIZE(ZWEST_2D_OUT),MNHREAL_MPI,NP_WEST-1,1000+IS_EAST,&
                     NMNH_COMM_WORLD,IREQ(INB_REQ),IERR)
#ifdef MNH_GPUDIRECT
      !$acc end host_data
#endif
   END IF
   IF (.NOT.GEAST) THEN 
#ifdef MNH_GPUDIRECT
      !$acc host_data use_device(ZEAST_2D_OUT)
#endif
      INB_REQ = INB_REQ + 1
      CALL MPI_IRECV(ZEAST_2D_OUT,SIZE(ZEAST_2D_OUT),MNHREAL_MPI,NP_EAST-1,1000+IS_WEST,&
                     NMNH_COMM_WORLD,IREQ(INB_REQ),IERR)
#ifdef MNH_GPUDIRECT
      !$acc end host_data
#endif
   ENDIF
END IF

IF (LY) THEN
   IF (.NOT.GSOUTH) THEN
#ifdef MNH_GPUDIRECT
      !$acc host_data use_device(ZSOUTH_2D_OUT)
#endif
      INB_REQ = INB_REQ + 1
      CALL MPI_IRECV(ZSOUTH_2D_OUT,SIZE(ZSOUTH_2D_OUT),MNHREAL_MPI,NP_SOUTH-1,1000+IS_NORTH,&
                     NMNH_COMM_WORLD,IREQ(INB_REQ),IERR)
#ifdef MNH_GPUDIRECT
      !$acc end host_data
#endif
   ENDIF
   IF (.NOT.GNORTH) THEN
#ifdef MNH_GPUDIRECT
      !$acc host_data use_device(ZNORTH_2D_OUT)
#endif
      INB_REQ = INB_REQ + 1
      CALL MPI_IRECV(ZNORTH_2D_OUT,SIZE(ZNORTH_2D_OUT),MNHREAL_MPI,NP_NORTH-1,1000+IS_SOUTH,&
                     NMNH_COMM_WORLD,IREQ(INB_REQ),IERR)
#ifdef MNH_GPUDIRECT
      !$acc end host_data
#endif
   ENDIF
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Copy the halo on the device PSRC to Zxxxx_2D_IN

IF (LX) THEN
   IF (.NOT. GWEST) THEN
   !$acc kernels async(IS_WEST)
   ZWEST_2D_IN ( IIB:IIB+IHALO_1  ,    IJB:IJE )  = PSRC( IIB:IIB+IHALO_1  ,  IJB:IJE )
   !$acc end kernels
      END IF
   IF (.NOT.GEAST) THEN
   !$acc kernels async(IS_EAST)
   ZEAST_2D_IN ( IIE-IHALO_1:IIE  ,    IJB:IJE )  = PSRC( IIE-IHALO_1:IIE  ,  IJB:IJE )
   !$acc end kernels
      ENDIF
END IF
IF (LY) THEN
   IF (.NOT.GSOUTH) THEN
   !$acc kernels async(IS_SOUTH)
   ZSOUTH_2D_IN ( IIB:IIE  ,    IJB:IJB+IHALO_1 ) = PSRC( IIB:IIE  ,    IJB:IJB+IHALO_1 )
   !$acc end kernels
      ENDIF
   IF (.NOT.GNORTH) THEN
   !$acc kernels async(IS_NORTH)
   ZNORTH_2D_IN ( IIB:IIE  ,    IJE-IHALO_1:IJE ) = PSRC( IIB:IIE  ,    IJE-IHALO_1:IJE )
   !$acc end kernels
   ENDIF
ENDIF
!$acc wait

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Send  Zxxxx_2D_IN buffer via MPI(Gpu_direct) or copy to host
!
IF (LX) THEN
   IF (.NOT. GWEST) THEN
#ifdef MNH_GPUDIRECT
      !$acc host_data use_device(ZWEST_2D_IN)
#else
      !$acc update host(ZWEST_2D_IN)
#endif
      INB_REQ = INB_REQ + 1
      CALL MPI_ISEND(ZWEST_2D_IN,SIZE(ZWEST_2D_IN)  ,MNHREAL_MPI,NP_WEST-1,1000+IS_WEST,&
                     NMNH_COMM_WORLD,IREQ(INB_REQ),IERR)
#ifdef MNH_GPUDIRECT
      !$acc end host_data
#endif
   END IF
   IF (.NOT.GEAST) THEN
#ifdef MNH_GPUDIRECT
      !$acc host_data use_device(ZEAST_2D_IN)
#else
      !$acc update host(ZEAST_2D_IN)
#endif
      INB_REQ = INB_REQ + 1
      CALL MPI_ISEND(ZEAST_2D_IN,SIZE(ZEAST_2D_IN)  ,MNHREAL_MPI,NP_EAST-1,1000+IS_EAST,&
                     NMNH_COMM_WORLD,IREQ(INB_REQ),IERR)
#ifdef MNH_GPUDIRECT
      !$acc end host_data
#endif
   ENDIF
END IF

IF (LY) THEN
   IF (.NOT.GSOUTH) THEN
#ifdef MNH_GPUDIRECT
      !$acc host_data use_device(ZSOUTH_2D_IN)
#else
      !$acc update host(ZSOUTH_2D_IN)
#endif
      INB_REQ = INB_REQ + 1
      CALL MPI_ISEND(ZSOUTH_2D_IN,SIZE(ZSOUTH_2D_IN)  ,MNHREAL_MPI,NP_SOUTH-1,1000+IS_SOUTH,&
                     NMNH_COMM_WORLD,IREQ(INB_REQ),IERR)
#ifdef MNH_GPUDIRECT
      !$acc end host_data
#endif
   ENDIF
   IF (.NOT.GNORTH) THEN
#ifdef MNH_GPUDIRECT
      !$acc host_data use_device(ZNORTH_2D_IN)
#else
      !$acc update host(ZNORTH_2D_IN)
#endif
      INB_REQ = INB_REQ + 1
      CALL MPI_ISEND(ZNORTH_2D_IN,SIZE(ZNORTH_2D_IN)  ,MNHREAL_MPI,NP_NORTH-1,1000+IS_NORTH,&
                     NMNH_COMM_WORLD,IREQ(INB_REQ),IERR)
#ifdef MNH_GPUDIRECT
      !$acc end host_data
#endif
   ENDIF   
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF ( INB_REQ > 0 ) THEN
   CALL MPI_WAITALL(INB_REQ,IREQ,MNH_STATUSES_IGNORE,IERR)
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Is update halo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


IF (LX) THEN
   IF (.NOT.GWEST) THEN
#ifndef MNH_GPUDIRECT
   !$acc update device(ZWEST_2D_OUT) async(IS_WEST)
#endif
   !$acc kernels async(IS_WEST)
   PSRC( 1:IIB-1  ,      IJB:IJE     ) = ZWEST_2D_OUT( 1:IIB-1  ,   IJB:IJE   )
   !$acc end kernels
   ENDIF
   IF (.NOT.GEAST) THEN
#ifndef MNH_GPUDIRECT
   !$acc update device(ZEAST_2D_OUT) async(IS_EAST)
#endif
   !$acc kernels async(IS_EAST)
   PSRC( IIE+1:IIU  ,      IJB:IJE     ) = ZEAST_2D_OUT( IIE+1:IIU  ,   IJB:IJE   )  
   !$acc end kernels
   ENDIF
END IF
IF (LY) THEN
   IF (.NOT.GSOUTH) THEN
#ifndef MNH_GPUDIRECT
   !$acc update device(ZSOUTH_2D_OUT) async(IS_SOUTH)
#endif
   !$acc kernels async(IS_SOUTH)
   PSRC(      IIB:IIE       ,  1:IJB-1) = ZSOUTH_2D_OUT(  IIB:IIE     , 1:IJB-1 )
   !$acc end kernels
   ENDIF
   IF (.NOT.GNORTH) THEN
#ifndef MNH_GPUDIRECT
   !$acc update device(ZNORTH_2D_OUT) async(IS_NORTH)
#endif
   !$acc kernels async(IS_NORTH)
   PSRC(      IIB:IIE       , IJE+1:IJU) = ZNORTH_2D_OUT (  IIB:IIE     , IJE+1:IJU  )
   !$acc end kernels
   ENDIF
END IF
!$acc wait

!$acc end data

CALL MPPDB_CHECK(PSRC,"GET_2D_HALO_DD end:PSRC")

END SUBROUTINE GET_2D_HALO_DD
!-------------------------------------------------------------------------------
!     ########################################
      SUBROUTINE GET_HALO_DDC(PSRC, HDIR, HNAME)
!     ########################################
!define MNH_GPUDIRECT
!
USE MODD_HALO_D
USE MODE_ll
USE MODD_ARGSLIST_ll, ONLY : LIST_ll
USE MODD_PARAMETERS, ONLY : JPHEXT
!
USE MODD_IO,        ONLY : GSMONOPROC
USE MODE_MNH_ZWORK, ONLY : GWEST , GEAST, GSOUTH , GNORTH
!
USE MODD_CONF, ONLY : NHALO
USE MODE_MPPDB

USE MODD_VAR_ll, ONLY    : IP,NPROC,NP1,NP2
USE MODD_VAR_ll, ONLY : NMNH_COMM_WORLD
USE MODD_VAR_ll,    ONLY : MNH_STATUSES_IGNORE
USE MODD_PRECISION,  ONLY : MNHREAL_MPI
!
USE MODE_MNH_ZWORK, ONLY : IIU,IJU,IKU
USE MODE_MNH_ZWORK, ONLY : IIB,IJB ,IIE,IJE 
!
IMPLICIT NONE
!
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PSRC    ! variable at t
CHARACTER(len=4), OPTIONAL :: HDIR ! to send only halo on X or Y direction
character(len=*), optional, intent(in) :: HNAME ! Name of the field to be added
!
character(len=:), allocatable    :: yname
TYPE(LIST_ll)     , POINTER      :: TZ_PSRC_ll               ! halo
INTEGER                          :: IERROR                 ! error return code 

INTEGER,PARAMETER :: IS_WEST=1 , IS_EAST=2, IS_SOUTH=3, IS_NORTH=4

LOGICAL      :: LX , LY

INTEGER      :: INB_REQEW , IREQEW(4)
INTEGER      :: INB_REQNS , IREQNS(4)
INTEGER      :: IERR

CALL MPPDB_CHECK(PSRC,"GET_HALO_DDC big:PSRC")

if ( NPROC == 1 ) then
!CALL MPPDB_CHECK(PSRC,'UPDATE_HALO_ll::'//TRIM(HNAME))
   CALL MPPDB_CHECK(PSRC,"GET_HALO_DDC end:PSRC")
   RETURN
   endif

CALL INIT_HALO_D()

!$acc data present ( PSRC )

NULLIFY( TZ_PSRC_ll)
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

LX = .FALSE.
LY = .FALSE. 

IF (.NOT. PRESENT(HDIR) ) THEN
LX = .TRUE.
LY = .TRUE.
ELSE
   !
   !  Problem of reproductibility in ppm_s0_x/y if only S0_X or S0_Y
   !  so add S0_X + S0_Y for ppm_s0*
   !
!!$LX = ( HDIR == "01_X" .OR. HDIR == "S0_X" )
!!$LY = ( HDIR == "01_Y" .OR. HDIR == "S0_Y" )
LX = ( HDIR == "01_X" .OR. HDIR == "S0_X" .OR. HDIR == "S0_Y" )
LY = ( HDIR == "01_Y" .OR. HDIR == "S0_Y" .OR. HDIR == "S0_X" )
END IF

!!$LX = .TRUE.
!!$LY = .TRUE.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Post first the recieve of ZxxxxC_OUT buffer via MPI(Gpu_direct)
!
!-------------------------------------------------------------------------------!
!  IRecv  E/W                                                                   !
!-------------------------------------------------------------------------------!
INB_REQEW = 0
IF (LX) THEN
   IF (.NOT. GWEST) THEN
#ifdef MNH_GPUDIRECT
      !$acc host_data use_device(ZWESTC_OUT)
#endif
      INB_REQEW = INB_REQEW + 1
      CALL MPI_IRECV(ZWESTC_OUT,SIZE(ZWESTC_OUT),MNHREAL_MPI,NP_WEST-1,1000+IS_EAST,&
                     NMNH_COMM_WORLD,IREQEW(INB_REQEW),IERR)
#ifdef MNH_GPUDIRECT
      !$acc end host_data
#endif
   END IF
   IF (.NOT.GEAST) THEN 
#ifdef MNH_GPUDIRECT
      !$acc host_data use_device(ZEASTC_OUT)
#endif
      INB_REQEW = INB_REQEW + 1
      CALL MPI_IRECV(ZEASTC_OUT,SIZE(ZEASTC_OUT),MNHREAL_MPI,NP_EAST-1,1000+IS_WEST,&
                     NMNH_COMM_WORLD,IREQEW(INB_REQEW),IERR)
#ifdef MNH_GPUDIRECT
      !$acc end host_data
#endif
   ENDIF
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Copy the halo E/W on the device PSRC to ZxxxxC_IN

IF (LX) THEN
   IF (.NOT. GWEST) THEN
   !$acc kernels async(IS_WEST)
   ZWESTC_IN ( IIB:IIB+IHALO_1  ,    IJB:IJE  , : )  = PSRC( IIB:IIB+IHALO_1  ,  IJB:IJE  , : )
   !$acc end kernels
   END IF
   IF (.NOT.GEAST) THEN
   !$acc kernels async(IS_EAST)
   ZEASTC_IN ( IIE-IHALO_1:IIE  ,    IJB:IJE  , : )  = PSRC( IIE-IHALO_1:IIE  ,  IJB:IJE  , : )
   !$acc end kernels
   ENDIF
   !$acc wait
END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Send  E/W ZxxxxC_IN buffer via MPI(Gpu_direct) or copy to host
!
IF (LX) THEN
   IF (.NOT. GWEST) THEN
#ifdef MNH_GPUDIRECT
      !$acc host_data use_device(ZWESTC_IN)
#else
      !$acc update host(ZWESTC_IN)
#endif
      INB_REQEW = INB_REQEW + 1
      CALL MPI_ISEND(ZWESTC_IN,SIZE(ZWESTC_IN)  ,MNHREAL_MPI,NP_WEST-1,1000+IS_WEST,&
                     NMNH_COMM_WORLD,IREQEW(INB_REQEW),IERR)
#ifdef MNH_GPUDIRECT
      !$acc end host_data
#endif
   END IF
   IF (.NOT.GEAST) THEN
#ifdef MNH_GPUDIRECT
      !$acc host_data use_device(ZEASTC_IN)
#else
      !$acc update host(ZEASTC_IN)
#endif
      INB_REQEW = INB_REQEW + 1
      CALL MPI_ISEND(ZEASTC_IN,SIZE(ZEASTC_IN)  ,MNHREAL_MPI,NP_EAST-1,1000+IS_EAST,&
                     NMNH_COMM_WORLD,IREQEW(INB_REQEW),IERR)
#ifdef MNH_GPUDIRECT
      !$acc end host_data
#endif
   ENDIF
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF ( INB_REQEW > 0 ) THEN
   CALL MPI_WAITALL(INB_REQEW,IREQEW,MNH_STATUSES_IGNORE,IERR)
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Update halo E/W from buffer to PSRC

IF (LX) THEN
   IF (.NOT.GWEST) THEN
#ifndef MNH_GPUDIRECT
   !$acc update device(ZWESTC_OUT) async(IS_WEST)
#endif
   !$acc kernels async(IS_WEST)
   PSRC( 1:IIB-1  ,      IJB:IJE      , : ) = ZWESTC_OUT( 1:IIB-1  ,   IJB:IJE    , : )
   !$acc end kernels
   ENDIF
   IF (.NOT.GEAST) THEN
#ifndef MNH_GPUDIRECT
   !$acc update device(ZEASTC_OUT) async(IS_EAST)
#endif
   !$acc kernels async(IS_EAST)
   PSRC( IIE+1:IIU  ,      IJB:IJE      , : ) = ZEASTC_OUT( IIE+1:IIU  ,   IJB:IJE    , : )  
   !$acc end kernels
   ENDIF
   !$acc wait
END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Post first the recieve of N/S ZxxxxC_OUT buffer via MPI(Gpu_direct)
!
!-------------------------------------------------------------------------------!
!  IRecv  N/S                                                                   !
!-------------------------------------------------------------------------------!
INB_REQNS = 0
IF (LY) THEN
   IF (.NOT.GSOUTH) THEN
#ifdef MNH_GPUDIRECT
      !$acc host_data use_device(ZSOUTHC_OUT)
#endif
      INB_REQNS = INB_REQNS + 1
      CALL MPI_IRECV(ZSOUTHC_OUT,SIZE(ZSOUTHC_OUT),MNHREAL_MPI,NP_SOUTH-1,1000+IS_NORTH,&
                     NMNH_COMM_WORLD,IREQNS(INB_REQNS),IERR)
#ifdef MNH_GPUDIRECT
      !$acc end host_data
#endif
   ENDIF
   IF (.NOT.GNORTH) THEN
#ifdef MNH_GPUDIRECT
      !$acc host_data use_device(ZNORTHC_OUT)
#endif
      INB_REQNS = INB_REQNS + 1
      CALL MPI_IRECV(ZNORTHC_OUT,SIZE(ZNORTHC_OUT),MNHREAL_MPI,NP_NORTH-1,1000+IS_SOUTH,&
                     NMNH_COMM_WORLD,IREQNS(INB_REQNS),IERR)
#ifdef MNH_GPUDIRECT
      !$acc end host_data
#endif
   ENDIF
END IF

!
!Copy the halo N/S on the device PSRC to ZxxxxC_IN
!
IF (LY) THEN
   IF (.NOT.GSOUTH) THEN
   !$acc kernels async(IS_SOUTH)
   ZSOUTHC_IN ( 1:IIU  ,    IJB:IJB+IHALO_1  , : ) = PSRC( 1:IIU  ,    IJB:IJB+IHALO_1  , : )
   !$acc end kernels
      ENDIF
   IF (.NOT.GNORTH) THEN
   !$acc kernels async(IS_NORTH)
   ZNORTHC_IN ( 1:IIU  ,    IJE-IHALO_1:IJE  , : ) = PSRC( 1:IIU  ,    IJE-IHALO_1:IJE  , : )
   !$acc end kernels
   ENDIF
   !$acc wait
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Send N/S ZxxxxC_IN buffer via MPI(Gpu_direct) or copy to host
!
IF (LY) THEN
   IF (.NOT.GSOUTH) THEN
#ifdef MNH_GPUDIRECT
      !$acc host_data use_device(ZSOUTHC_IN)
#else
      !$acc update host(ZSOUTHC_IN)
#endif
      INB_REQNS = INB_REQNS + 1
      CALL MPI_ISEND(ZSOUTHC_IN,SIZE(ZSOUTHC_IN)  ,MNHREAL_MPI,NP_SOUTH-1,1000+IS_SOUTH,&
                     NMNH_COMM_WORLD,IREQNS(INB_REQNS),IERR)
#ifdef MNH_GPUDIRECT
      !$acc end host_data
#endif
   ENDIF
   IF (.NOT.GNORTH) THEN
#ifdef MNH_GPUDIRECT
      !$acc host_data use_device(ZNORTHC_IN)
#else
      !$acc update host(ZNORTHC_IN)
#endif
      INB_REQNS = INB_REQNS + 1
      CALL MPI_ISEND(ZNORTHC_IN,SIZE(ZNORTHC_IN)  ,MNHREAL_MPI,NP_NORTH-1,1000+IS_NORTH,&
                     NMNH_COMM_WORLD,IREQNS(INB_REQNS),IERR)
#ifdef MNH_GPUDIRECT
      !$acc end host_data
#endif
   ENDIF   
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF ( INB_REQNS > 0 ) THEN
   CALL MPI_WAITALL(INB_REQNS,IREQNS,MNH_STATUSES_IGNORE,IERR)
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Update halo N/S

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Update halo N/S/W from buffer to PSRC

IF (LY) THEN
   IF (.NOT.GSOUTH) THEN
#ifndef MNH_GPUDIRECT
   !$acc update device(ZSOUTHC_OUT) async(IS_SOUTH)
#endif
   !$acc kernels async(IS_SOUTH)
   PSRC(      1:IIU       ,  1:IJB-1 , : ) = ZSOUTHC_OUT(  1:IIU     , 1:IJB-1  , : )
   !$acc end kernels
   ENDIF
   IF (.NOT.GNORTH) THEN
#ifndef MNH_GPUDIRECT
   !$acc update device(ZNORTHC_OUT) async(IS_NORTH)
#endif
   !$acc kernels async(IS_NORTH)
   PSRC(      1:IIU       , IJE+1:IJU , : ) = ZNORTHC_OUT (  1:IIU     , IJE+1:IJU  , : )
   !$acc end kernels
   ENDIF
   !$acc wait
END IF

!$acc end data

CALL MPPDB_CHECK(PSRC,"GET_HALO_DDC end:PSRC")
!CALL MPPDB_CHECK(PSRC,'UPDATE_HALO_ll::'//TRIM(HNAME))

END SUBROUTINE GET_HALO_DDC
!-------------------------------------------------------------------------------
!     ########################################
      SUBROUTINE GET_2D_HALO_DDC(PSRC, HDIR, HNAME)
!     ########################################
!define MNH_GPUDIRECT
!
USE MODD_HALO_D
USE MODE_ll
USE MODD_ARGSLIST_ll, ONLY : LIST_ll
USE MODD_PARAMETERS, ONLY : JPHEXT
!
USE MODD_IO,        ONLY : GSMONOPROC
USE MODE_MNH_ZWORK, ONLY : GWEST , GEAST, GSOUTH , GNORTH
!
USE MODD_CONF, ONLY : NHALO
USE MODE_MPPDB

USE MODD_VAR_ll, ONLY    : IP,NPROC,NP1,NP2
USE MODD_VAR_ll, ONLY : NMNH_COMM_WORLD
USE MODD_VAR_ll,    ONLY : MNH_STATUSES_IGNORE
USE MODD_PRECISION,  ONLY : MNHREAL_MPI
!
USE MODE_MNH_ZWORK, ONLY : IIU,IJU,IKU
USE MODE_MNH_ZWORK, ONLY : IIB,IJB ,IIE,IJE 
!
IMPLICIT NONE
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSRC    ! variable at t
CHARACTER(len=4), OPTIONAL :: HDIR ! to send only halo on X or Y direction
character(len=*), optional, intent(in) :: HNAME ! Name of the field to be added
!
character(len=:), allocatable    :: yname
TYPE(LIST_ll)     , POINTER      :: TZ_PSRC_ll               ! halo
INTEGER                          :: IERROR                 ! error return code 

INTEGER,PARAMETER :: IS_WEST=1 , IS_EAST=2, IS_SOUTH=3, IS_NORTH=4

LOGICAL      :: LX , LY

INTEGER      :: INB_REQEW , IREQEW(4)
INTEGER      :: INB_REQNS , IREQNS(4)
INTEGER      :: IERR

CALL MPPDB_CHECK(PSRC,"GET_2D_HALO_DDC big:PSRC")

if ( NPROC == 1 ) then
   CALL MPPDB_CHECK(PSRC,"GET_2D_HALO_DDC end:PSRC")
   !CALL MPPDB_CHECK(PSRC,'UPDATE_HALO_ll::'//TRIM(HNAME))
   RETURN
end if

CALL INIT_HALO_D()

!$acc data present ( PSRC )

NULLIFY( TZ_PSRC_ll)
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

LX = .FALSE.
LY = .FALSE. 

IF (.NOT. PRESENT(HDIR) ) THEN
LX = .TRUE.
LY = .TRUE.
ELSE
   !
   !  Problem of reproductibility in ppm_s0_x/y if only S0_X or S0_Y
   !  so add S0_X + S0_Y for ppm_s0*
   !
!!$LX = ( HDIR == "01_X" .OR. HDIR == "S0_X" )
!!$LY = ( HDIR == "01_Y" .OR. HDIR == "S0_Y" )
LX = ( HDIR == "01_X" .OR. HDIR == "S0_X" .OR. HDIR == "S0_Y" )
LY = ( HDIR == "01_Y" .OR. HDIR == "S0_Y" .OR. HDIR == "S0_X" )
END IF

!!$LX = .TRUE.
!!$LY = .TRUE.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Post first the recieve of ZxxxxC_2D_OUT buffer via MPI(Gpu_direct)
!
!-------------------------------------------------------------------------------!
!  IRecv  E/W                                                                   !
!-------------------------------------------------------------------------------!
INB_REQEW = 0
IF (LX) THEN
   IF (.NOT. GWEST) THEN
#ifdef MNH_GPUDIRECT
      !$acc host_data use_device(ZWESTC_2D_OUT)
#endif
      INB_REQEW = INB_REQEW + 1
      CALL MPI_IRECV(ZWESTC_2D_OUT,SIZE(ZWESTC_2D_OUT),MNHREAL_MPI,NP_WEST-1,1000+IS_EAST,&
                     NMNH_COMM_WORLD,IREQEW(INB_REQEW),IERR)
#ifdef MNH_GPUDIRECT
      !$acc end host_data
#endif
   END IF
   IF (.NOT.GEAST) THEN 
#ifdef MNH_GPUDIRECT
      !$acc host_data use_device(ZEASTC_2D_OUT)
#endif
      INB_REQEW = INB_REQEW + 1
      CALL MPI_IRECV(ZEASTC_2D_OUT,SIZE(ZEASTC_2D_OUT),MNHREAL_MPI,NP_EAST-1,1000+IS_WEST,&
                     NMNH_COMM_WORLD,IREQEW(INB_REQEW),IERR)
#ifdef MNH_GPUDIRECT
      !$acc end host_data
#endif
   ENDIF
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Copy the halo E/W on the device PSRC to ZxxxxC_2D_IN

IF (LX) THEN
   IF (.NOT. GWEST) THEN
   !$acc kernels async(IS_WEST)
   ZWESTC_2D_IN ( IIB:IIB+IHALO_1  ,    IJB:IJE )  = PSRC( IIB:IIB+IHALO_1  ,  IJB:IJE )
   !$acc end kernels
   END IF
   IF (.NOT.GEAST) THEN
   !$acc kernels async(IS_EAST)
   ZEASTC_2D_IN ( IIE-IHALO_1:IIE  ,    IJB:IJE )  = PSRC( IIE-IHALO_1:IIE  ,  IJB:IJE )
   !$acc end kernels
   ENDIF
   !$acc wait
END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Send  E/W ZxxxxC_2D_IN buffer via MPI(Gpu_direct) or copy to host
!
IF (LX) THEN
   IF (.NOT. GWEST) THEN
#ifdef MNH_GPUDIRECT
      !$acc host_data use_device(ZWESTC_2D_IN)
#else
      !$acc update host(ZWESTC_2D_IN)
#endif
      INB_REQEW = INB_REQEW + 1
      CALL MPI_ISEND(ZWESTC_2D_IN,SIZE(ZWESTC_2D_IN)  ,MNHREAL_MPI,NP_WEST-1,1000+IS_WEST,&
                     NMNH_COMM_WORLD,IREQEW(INB_REQEW),IERR)
#ifdef MNH_GPUDIRECT
      !$acc end host_data
#endif
   END IF
   IF (.NOT.GEAST) THEN
#ifdef MNH_GPUDIRECT
      !$acc host_data use_device(ZEASTC_2D_IN)
#else
      !$acc update host(ZEASTC_2D_IN)
#endif
      INB_REQEW = INB_REQEW + 1
      CALL MPI_ISEND(ZEASTC_2D_IN,SIZE(ZEASTC_2D_IN)  ,MNHREAL_MPI,NP_EAST-1,1000+IS_EAST,&
                     NMNH_COMM_WORLD,IREQEW(INB_REQEW),IERR)
#ifdef MNH_GPUDIRECT
      !$acc end host_data
#endif
   ENDIF
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF ( INB_REQEW > 0 ) THEN
   CALL MPI_WAITALL(INB_REQEW,IREQEW,MNH_STATUSES_IGNORE,IERR)
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Update halo E/W from buffer to PSRC

IF (LX) THEN
   IF (.NOT.GWEST) THEN
#ifndef MNH_GPUDIRECT
   !$acc update device(ZWESTC_2D_OUT) async(IS_WEST)
#endif
   !$acc kernels async(IS_WEST)
   PSRC( 1:IIB-1  ,      IJB:IJE     ) = ZWESTC_2D_OUT( 1:IIB-1  ,   IJB:IJE   )
   !$acc end kernels
   ENDIF
   IF (.NOT.GEAST) THEN
#ifndef MNH_GPUDIRECT
   !$acc update device(ZEASTC_2D_OUT) async(IS_EAST)
#endif
   !$acc kernels async(IS_EAST)
   PSRC( IIE+1:IIU  ,      IJB:IJE     ) = ZEASTC_2D_OUT( IIE+1:IIU  ,   IJB:IJE   )  
   !$acc end kernels
   ENDIF
   !$acc wait
END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Post first the recieve of N/S ZxxxxC_2D_OUT buffer via MPI(Gpu_direct)
!
!-------------------------------------------------------------------------------!
!  IRecv  N/S                                                                   !
!-------------------------------------------------------------------------------!
INB_REQNS = 0
IF (LY) THEN
   IF (.NOT.GSOUTH) THEN
#ifdef MNH_GPUDIRECT
      !$acc host_data use_device(ZSOUTHC_2D_OUT)
#endif
      INB_REQNS = INB_REQNS + 1
      CALL MPI_IRECV(ZSOUTHC_2D_OUT,SIZE(ZSOUTHC_2D_OUT),MNHREAL_MPI,NP_SOUTH-1,1000+IS_NORTH,&
                     NMNH_COMM_WORLD,IREQNS(INB_REQNS),IERR)
#ifdef MNH_GPUDIRECT
      !$acc end host_data
#endif
   ENDIF
   IF (.NOT.GNORTH) THEN
#ifdef MNH_GPUDIRECT
      !$acc host_data use_device(ZNORTHC_2D_OUT)
#endif
      INB_REQNS = INB_REQNS + 1
      CALL MPI_IRECV(ZNORTHC_2D_OUT,SIZE(ZNORTHC_2D_OUT),MNHREAL_MPI,NP_NORTH-1,1000+IS_SOUTH,&
                     NMNH_COMM_WORLD,IREQNS(INB_REQNS),IERR)
#ifdef MNH_GPUDIRECT
      !$acc end host_data
#endif
   ENDIF
END IF

!
!Copy the halo N/S on the device PSRC to ZxxxxC_2D_IN
!
IF (LY) THEN
   IF (.NOT.GSOUTH) THEN
   !$acc kernels async(IS_SOUTH)
   ZSOUTHC_2D_IN ( 1:IIU  ,    IJB:IJB+IHALO_1 ) = PSRC( 1:IIU  ,    IJB:IJB+IHALO_1 )
   !$acc end kernels
      ENDIF
   IF (.NOT.GNORTH) THEN
   !$acc kernels async(IS_NORTH)
   ZNORTHC_2D_IN ( 1:IIU  ,    IJE-IHALO_1:IJE ) = PSRC( 1:IIU  ,    IJE-IHALO_1:IJE )
   !$acc end kernels
   ENDIF
   !$acc wait
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Send N/S ZxxxxC_2D_IN buffer via MPI(Gpu_direct) or copy to host
!
IF (LY) THEN
   IF (.NOT.GSOUTH) THEN
#ifdef MNH_GPUDIRECT
      !$acc host_data use_device(ZSOUTHC_2D_IN)
#else
      !$acc update host(ZSOUTHC_2D_IN)
#endif
      INB_REQNS = INB_REQNS + 1
      CALL MPI_ISEND(ZSOUTHC_2D_IN,SIZE(ZSOUTHC_2D_IN)  ,MNHREAL_MPI,NP_SOUTH-1,1000+IS_SOUTH,&
                     NMNH_COMM_WORLD,IREQNS(INB_REQNS),IERR)
#ifdef MNH_GPUDIRECT
      !$acc end host_data
#endif
   ENDIF
   IF (.NOT.GNORTH) THEN
#ifdef MNH_GPUDIRECT
      !$acc host_data use_device(ZNORTHC_2D_IN)
#else
      !$acc update host(ZNORTHC_2D_IN)
#endif
      INB_REQNS = INB_REQNS + 1
      CALL MPI_ISEND(ZNORTHC_2D_IN,SIZE(ZNORTHC_2D_IN)  ,MNHREAL_MPI,NP_NORTH-1,1000+IS_NORTH,&
                     NMNH_COMM_WORLD,IREQNS(INB_REQNS),IERR)
#ifdef MNH_GPUDIRECT
      !$acc end host_data
#endif
   ENDIF   
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF ( INB_REQNS > 0 ) THEN
   CALL MPI_WAITALL(INB_REQNS,IREQNS,MNH_STATUSES_IGNORE,IERR)
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Update halo N/S

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Update halo N/S/W from buffer to PSRC

IF (LY) THEN
   IF (.NOT.GSOUTH) THEN
#ifndef MNH_GPUDIRECT
   !$acc update device(ZSOUTHC_2D_OUT) async(IS_SOUTH)
#endif
   !$acc kernels async(IS_SOUTH)
   PSRC(      1:IIU       ,  1:IJB-1) = ZSOUTHC_2D_OUT(  1:IIU     , 1:IJB-1 )
   !$acc end kernels
   ENDIF
   IF (.NOT.GNORTH) THEN
#ifndef MNH_GPUDIRECT
   !$acc update device(ZNORTHC_2D_OUT) async(IS_NORTH)
#endif
   !$acc kernels async(IS_NORTH)
   PSRC(      1:IIU       , IJE+1:IJU) = ZNORTHC_2D_OUT (  1:IIU     , IJE+1:IJU )
   !$acc end kernels
   ENDIF
   !$acc wait
END IF

!$acc end data

CALL MPPDB_CHECK(PSRC,"GET_2D_HALO_DDC end:PSRC")
!CALL MPPDB_CHECK(PSRC,'UPDATE_HALO_ll::'//TRIM(HNAME))

END SUBROUTINE GET_2D_HALO_DDC
!-------------------------------------------------------------------------------
!     ########################################
      SUBROUTINE GET_HALO2_DD(PSRC, TP_PSRC_HALO2_ll, HNAME)
!     ########################################
!define MNH_GPUDIRECT
!
USE MODD_HALO_D
USE MODE_ll
USE MODD_ARGSLIST_ll, ONLY : LIST_ll
USE MODD_PARAMETERS, ONLY : JPHEXT
!
USE MODD_IO,        ONLY : GSMONOPROC
USE MODE_MNH_ZWORK, ONLY : GWEST , GEAST, GSOUTH , GNORTH
USE MODE_MNH_ZWORK, ONLY : IIU,IJU,IKU
USE MODE_MNH_ZWORK, ONLY : IIB,IJB ,IIE,IJE
!
USE MODD_CONF, ONLY : NHALO
USE MODE_MPPDB
!
USE MODD_VAR_ll,     ONLY : IP,NPROC,NP1,NP2
USE MODD_VAR_ll,     ONLY : NMNH_COMM_WORLD
USE MODD_VAR_ll,     ONLY : MNH_STATUSES_IGNORE
USE MODD_PRECISION,  ONLY : MNHREAL_MPI
!
IMPLICIT NONE
!
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PSRC    ! variable at t
TYPE(HALO2LIST_ll), POINTER  :: TP_PSRC_HALO2_ll  ! halo2 for SRC
character(len=*), optional, intent(in) :: HNAME ! Name of the field to be added
!
character(len=:), allocatable    :: yname
INTEGER                          :: IERROR                 ! error return code 

INTEGER,PARAMETER :: IS_WEST=1 , IS_EAST=2, IS_SOUTH=3, IS_NORTH=4

LOGICAL      :: LX , LY

INTEGER      :: INB_REQ , IREQ(8)
INTEGER      :: IERR

REAL , DIMENSION(:,:) , POINTER , CONTIGUOUS :: ZH2_EAST,ZH2_WEST,ZH2_NORTH,ZH2_SOUTH

CALL MPPDB_CHECK(PSRC,"GET_HALO2_DD big:PSRC")

if ( NPROC == 1 ) then
   CALL MPPDB_CHECK(PSRC,"GET_HALO2_DD end:PSRC")
   RETURN
end if

CALL INIT_HALO_D()

!$acc data present ( PSRC ) &
!$acc present (ZNORTH2_IN, ZSOUTH2_IN, ZWEST2_IN, ZEAST2_IN) &
!$acc present (ZNORTH2_OUT, ZSOUTH2_OUT, ZWEST2_OUT, ZEAST2_OUT)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

LX = .TRUE.
LY = .TRUE.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

INB_REQ = 0

!
! Post the recieve of Zxxxx_IN buffer first via MPI(Gpu_direct)
!

IF (LX) THEN
   IF (.NOT. GWEST) THEN
#ifdef MNH_GPUDIRECT
      !$acc host_data use_device(ZWEST2_OUT)
#endif
      INB_REQ = INB_REQ + 1
      CALL MPI_IRECV(ZWEST2_OUT,SIZE(ZWEST2_OUT),MNHREAL_MPI,NP_WEST-1,1000+IS_EAST,&
                     NMNH_COMM_WORLD,IREQ(INB_REQ),IERR)
#ifdef MNH_GPUDIRECT
      !$acc end host_data
#endif
   END IF
   IF (.NOT.GEAST) THEN 
#ifdef MNH_GPUDIRECT
      !$acc host_data use_device(ZEAST2_OUT)
#endif
      INB_REQ = INB_REQ + 1
      CALL MPI_IRECV(ZEAST2_OUT,SIZE(ZEAST2_OUT),MNHREAL_MPI,NP_EAST-1,1000+IS_WEST,&
                     NMNH_COMM_WORLD,IREQ(INB_REQ),IERR)
#ifdef MNH_GPUDIRECT
      !$acc end host_data
#endif
   ENDIF
END IF

IF (LY) THEN
   IF (.NOT.GSOUTH) THEN
#ifdef MNH_GPUDIRECT
      !$acc host_data use_device(ZSOUTH2_OUT)
#endif
      INB_REQ = INB_REQ + 1
      CALL MPI_IRECV(ZSOUTH2_OUT,SIZE(ZSOUTH2_OUT),MNHREAL_MPI,NP_SOUTH-1,1000+IS_NORTH,&
                     NMNH_COMM_WORLD,IREQ(INB_REQ),IERR)
#ifdef MNH_GPUDIRECT
      !$acc end host_data
#endif
   ENDIF
   IF (.NOT.GNORTH) THEN
#ifdef MNH_GPUDIRECT
      !$acc host_data use_device(ZNORTH2_OUT)
#endif
      INB_REQ = INB_REQ + 1
      CALL MPI_IRECV(ZNORTH2_OUT,SIZE(ZNORTH2_OUT),MNHREAL_MPI,NP_NORTH-1,1000+IS_SOUTH,&
                     NMNH_COMM_WORLD,IREQ(INB_REQ),IERR)
#ifdef MNH_GPUDIRECT
      !$acc end host_data
#endif
   ENDIF
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Copy the halo on the device PSRC to Zxxxx_IN

IF (LX) THEN
   IF (.NOT. GWEST) THEN
   !$acc kernels async(IS_WEST)
!!$      ZWEST2_IN ( IIB:IIB+IHALO_1  ,    IJB:IJE  , : )  = PSRC( IIB:IIB+IHALO_1  ,  IJB:IJE  , : )
      ZWEST2_IN ( : , : )  = PSRC( IIB+1  , : , : )
   !$acc end kernels
      END IF
   IF (.NOT.GEAST) THEN
   !$acc kernels async(IS_EAST)
!!$      ZEAST2_IN ( IIE-IHALO_1:IIE  ,    IJB:IJE  , : )  = PSRC( IIE-IHALO_1:IIE  ,  IJB:IJE  , : )
      ZEAST2_IN ( : , : )  = PSRC( IIE-1 ,  :  , : )
   !$acc end kernels
      ENDIF
END IF
IF (LY) THEN
   IF (.NOT.GSOUTH) THEN
   !$acc kernels async(IS_SOUTH)
!!$   ZSOUTH2_IN ( IIB:IIE  ,    IJB:IJB+IHALO_1  , : ) = PSRC( IIB:IIE  ,    IJB:IJB+IHALO_1  , : )
      ZSOUTH2_IN ( : , : ) = PSRC( : , IJB+1 , : )
   !$acc end kernels
      ENDIF
   IF (.NOT.GNORTH) THEN
   !$acc kernels async(IS_NORTH)
!!$      ZNORTH2_IN ( IIB:IIE  ,    IJE-IHALO_1:IJE  , : ) = PSRC( IIB:IIE  ,    IJE-IHALO_1:IJE  , : )
      ZNORTH2_IN ( : , : ) = PSRC( : , IJE-1  , : )      
   !$acc end kernels
   ENDIF
ENDIF
!$acc wait

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Send  Zxxxx2_IN buffer via MPI(Gpu_direct) or copy to host
!
IF (LX) THEN
   IF (.NOT. GWEST) THEN
#ifdef MNH_GPUDIRECT
      !$acc host_data use_device(ZWEST2_IN)
#else
      !$acc update host(ZWEST2_IN)
#endif
      INB_REQ = INB_REQ + 1
      CALL MPI_ISEND(ZWEST2_IN,SIZE(ZWEST2_IN)  ,MNHREAL_MPI,NP_WEST-1,1000+IS_WEST,&
                     NMNH_COMM_WORLD,IREQ(INB_REQ),IERR)
#ifdef MNH_GPUDIRECT
      !$acc end host_data
#endif
   END IF
   IF (.NOT.GEAST) THEN
#ifdef MNH_GPUDIRECT
      !$acc host_data use_device(ZEAST2_IN)
#else
      !$acc update host(ZEAST2_IN)
#endif
      INB_REQ = INB_REQ + 1
      CALL MPI_ISEND(ZEAST2_IN,SIZE(ZEAST2_IN)  ,MNHREAL_MPI,NP_EAST-1,1000+IS_EAST,&
                     NMNH_COMM_WORLD,IREQ(INB_REQ),IERR)
#ifdef MNH_GPUDIRECT
      !$acc end host_data
#endif
   ENDIF
END IF

IF (LY) THEN
   IF (.NOT.GSOUTH) THEN
#ifdef MNH_GPUDIRECT
      !$acc host_data use_device(ZSOUTH2_IN)
#else
      !$acc update host(ZSOUTH2_IN)
#endif
      INB_REQ = INB_REQ + 1
      CALL MPI_ISEND(ZSOUTH2_IN,SIZE(ZSOUTH2_IN)  ,MNHREAL_MPI,NP_SOUTH-1,1000+IS_SOUTH,&
                     NMNH_COMM_WORLD,IREQ(INB_REQ),IERR)
#ifdef MNH_GPUDIRECT
      !$acc end host_data
#endif
   ENDIF
   IF (.NOT.GNORTH) THEN
#ifdef MNH_GPUDIRECT
      !$acc host_data use_device(ZNORTH2_IN)
#else
      !$acc update host(ZNORTH2_IN)
#endif
      INB_REQ = INB_REQ + 1
      CALL MPI_ISEND(ZNORTH2_IN,SIZE(ZNORTH2_IN)  ,MNHREAL_MPI,NP_NORTH-1,1000+IS_NORTH,&
                     NMNH_COMM_WORLD,IREQ(INB_REQ),IERR)
#ifdef MNH_GPUDIRECT
      !$acc end host_data
#endif
   ENDIF   
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF ( INB_REQ > 0 ) THEN
   CALL MPI_WAITALL(INB_REQ,IREQ,MNH_STATUSES_IGNORE,IERR)
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Is update halo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


IF (LX) THEN
   IF (.NOT.GWEST) THEN
#ifndef MNH_GPUDIRECT
   !$acc update device(ZWEST2_OUT) async(IS_WEST)
#endif
   ZH2_WEST => TP_PSRC_HALO2_ll%HALO2%WEST  
   !$acc kernels async(IS_WEST)
!!$      PSRC( 1:IIB-1  ,      IJB:IJE      , : ) = ZWEST2_OUT( 1:IIB-1  ,   IJB:IJE    , : )
      ZH2_WEST( : , : ) = ZWEST2_OUT( : , : )      
   !$acc end kernels
   ENDIF
   IF (.NOT.GEAST) THEN
#ifndef MNH_GPUDIRECT
   !$acc update device(ZEAST2_OUT) async(IS_EAST)
#endif
      ZH2_EAST => TP_PSRC_HALO2_ll%HALO2%EAST
   !$acc kernels async(IS_EAST)
!!$      PSRC( IIE+1:IIU  ,      IJB:IJE      , : ) = ZEAST2_OUT( IIE+1:IIU  ,   IJB:IJE    , : )
      ZH2_EAST( : , : ) = ZEAST2_OUT( : , : )
   !$acc end kernels
   ENDIF
END IF
IF (LY) THEN
   IF (.NOT.GSOUTH) THEN
#ifndef MNH_GPUDIRECT
   !$acc update device(ZSOUTH2_OUT) async(IS_SOUTH)
#endif
   ZH2_SOUTH => TP_PSRC_HALO2_ll%HALO2%SOUTH    
   !$acc kernels async(IS_SOUTH)
!!$      PSRC(      IIB:IIE       ,  1:IJB-1 , : ) = ZSOUTH2_OUT(  IIB:IIE     , 1:IJB-1  , : )
     ZH2_SOUTH( : , : ) = ZSOUTH2_OUT( : , : )  
   !$acc end kernels
   ENDIF
   IF (.NOT.GNORTH) THEN
#ifndef MNH_GPUDIRECT
   !$acc update device(ZNORTH2_OUT) async(IS_NORTH)
#endif
   ZH2_NORTH => TP_PSRC_HALO2_ll%HALO2%NORTH
   !$acc kernels async(IS_NORTH)
!!$      PSRC(      IIB:IIE       , IJE+1:IJU , : ) = ZNORTH2_OUT (  IIB:IIE     , IJE+1:IJU  , : )
     ZH2_NORTH( : , : ) = ZNORTH2_OUT ( : , : )      
   !$acc end kernels
   ENDIF
END IF
!$acc wait

!$acc end data

CALL MPPDB_CHECK(PSRC,"GET_HALO2_DD end:PSRC")

END SUBROUTINE GET_HALO2_DD
!-------------------------------------------------------------------------------
!     ########################################
SUBROUTINE GET_HALO2_DF(PSRC, TP_PSRC_HALO2F_ll, HNAME)
!     ########################################
!define MNH_GPUDIRECT
  !
USE MODD_HALO_D  
USE MODD_ARGSLIST_ll, ONLY : HALO2LIST_ll
USE MODE_MNH_ZWORK, ONLY : IIU,IJU,IKU
USE MODE_MNH_ZWORK, ONLY : IIB,IJB ,IIE,IJE  
IMPLICIT NONE
!
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PSRC    ! variable at t
TYPE(HALO2LIST_ll), POINTER  :: TP_PSRC_HALO2F_ll  ! halo2 for SRC
character(len=*), optional, intent(in) :: HNAME ! Name of the field to be added
!
!  local var
!
REAL , DIMENSION(:,:) , POINTER , CONTIGUOUS :: ZH2F_EAST,ZH2F_WEST,ZH2F_NORTH,ZH2F_SOUTH
!
CALL INIT_HALO_D()
!
ZH2F_EAST => TP_PSRC_HALO2F_ll%HALO2%EAST
ZH2F_WEST => TP_PSRC_HALO2F_ll%HALO2%WEST
ZH2F_NORTH => TP_PSRC_HALO2F_ll%HALO2%NORTH
ZH2F_SOUTH => TP_PSRC_HALO2F_ll%HALO2%SOUTH
!
CALL  GET_HALO2_DF_DIM(PSRC,&
      ZSOUTH2F_IN, ZNORTH2F_IN, ZWEST2F_IN, ZEAST2F_IN,&
      ZSOUTH2F_OUT, ZNORTH2F_OUT, ZWEST2F_OUT, ZEAST2F_OUT,&
      ZH2F_EAST,ZH2F_WEST,ZH2F_NORTH,ZH2F_SOUTH,&
      IIB,IIE,IJB,IJE,IIU,IJU,IKU,IHALO2_1,IHALO2,&
      HNAME)
!
CONTAINS
!
SUBROUTINE GET_HALO2_DF_DIM(PSRC,&
      PZSOUTH2F_IN, PZNORTH2F_IN, PZWEST2F_IN, PZEAST2F_IN,&
      PZSOUTH2F_OUT, PZNORTH2F_OUT, PZWEST2F_OUT, PZEAST2F_OUT,&
      PZH2F_EAST,PZH2F_WEST,PZH2F_NORTH,PZH2F_SOUTH,&
      KIIB,KIIE,KIJB,KIJE,KIIU,KIJU,KIKU,KIHALO2_1,KIHALO2,&
      HNAME)  
!        
USE MODE_ll
!
USE MODE_MNH_ZWORK, ONLY : GWEST , GEAST, GSOUTH , GNORTH

USE MODE_MPPDB
!
USE MODD_VAR_ll,     ONLY : NPROC
USE MODD_VAR_ll,     ONLY : NMNH_COMM_WORLD
USE MODD_VAR_ll,     ONLY : MNH_STATUSES_IGNORE
USE MODD_PRECISION,  ONLY : MNHREAL_MPI
!
IMPLICIT NONE
!
REAL, DIMENSION(KIIU,KIJU,KIKU), INTENT(INOUT) :: PSRC    ! variable at t
REAL :: PZSOUTH2F_IN ( KIIB:KIIE   , KIJB:KIJB+KIHALO2_1   , KIKU ) ,&
        PZNORTH2F_IN ( KIIB:KIIE   , KIJE-KIHALO2_1:KIJE   , KIKU ) ,&
        PZWEST2F_IN  ( KIIB:KIIB+KIHALO2_1   , KIJB:KIJE   , KIKU ) ,&
        PZEAST2F_IN  ( KIIE-KIHALO2_1:KIIE   , KIJB:KIJE   , KIKU ) ,&
        !
        PZSOUTH2F_OUT (   KIIB:KIIE   , KIJB-KIHALO2:KIJB-1 , KIKU ) ,&
        PZNORTH2F_OUT (   KIIB:KIIE   , KIJE+1:KIJE+KIHALO2 , KIKU ) ,&
        PZWEST2F_OUT  ( KIIB-KIHALO2:KIIB-1 ,   KIJB:KIJE   , KIKU ) ,&
        PZEAST2F_OUT  ( KIIE+1:KIIE+KIHALO2 ,   KIJB:KIJE   , KIKU ) ,&
        !
        PZH2F_EAST ( KIJU , KIKU ) ,&
        PZH2F_WEST ( KIJU , KIKU ) ,&
        PZH2F_NORTH( KIIU , KIKU ) ,&
        PZH2F_SOUTH( KIIU , KIKU )
INTEGER ::  KIIB,KIIE,KIJB,KIJE,KIIU,KIJU,KIKU,KIHALO2_1,KIHALO2
character(len=*), optional, intent(in) :: HNAME ! Name of the field to be added
!
character(len=:), allocatable    :: yname
INTEGER                          :: IERROR                 ! error return code 

INTEGER,PARAMETER :: IS_WEST=1 , IS_EAST=2, IS_SOUTH=3, IS_NORTH=4

LOGICAL      :: LX , LY

INTEGER      :: INB_REQ , IREQ(8)
INTEGER      :: IERR

!!$REAL , DIMENSION(:,:) , POINTER , CONTIGUOUS :: ZH2F_EAST,ZH2F_WEST,ZH2F_NORTH,ZH2F_SOUTH

CALL MPPDB_CHECK(PSRC,"GET_HALO2_DF beg:PSRC")

if ( NPROC == 1 ) then
   CALL MPPDB_CHECK(PSRC,"GET_HALO2_DF end:PSRC")
   !CALL MPPDB_CHECK(PSRC,'UPDATE_HALO_ll::'//TRIM(HNAME))
   !CALL MPPDB_CHECK(PSRC,'UPDATE_HALO2_ll::'//TRIM(HNAME))
   RETURN
end if

CALL INIT_HALO_D()

!$acc data present ( PSRC ) &
!$acc & present (PZNORTH2F_IN, PZSOUTH2F_IN, PZWEST2F_IN, PZEAST2F_IN) &
!$acc & present (PZNORTH2F_OUT, PZSOUTH2F_OUT, PZWEST2F_OUT, PZEAST2F_OUT) &
!$acc & present (PZH2F_EAST,PZH2F_WEST,PZH2F_NORTH,PZH2F_SOUTH)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

LX = .TRUE.
LY = .TRUE.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

INB_REQ = 0

!
! Post the recieve of Zxxxx_IN buffer first via MPI(Gpu_direct)
!

IF (LX) THEN
   IF (.NOT. GWEST) THEN
#ifdef MNH_GPUDIRECT
      !$acc host_data use_device(PZWEST2F_OUT)
#endif
      INB_REQ = INB_REQ + 1
      CALL MPI_IRECV(PZWEST2F_OUT,SIZE(PZWEST2F_OUT),MNHREAL_MPI,NP_WEST-1,1000+IS_EAST,&
                     NMNH_COMM_WORLD,IREQ(INB_REQ),IERR)
#ifdef MNH_GPUDIRECT
      !$acc end host_data
#endif
   END IF
   IF (.NOT.GEAST) THEN 
#ifdef MNH_GPUDIRECT
      !$acc host_data use_device(PZEAST2F_OUT)
#endif
      INB_REQ = INB_REQ + 1
      CALL MPI_IRECV(PZEAST2F_OUT,SIZE(PZEAST2F_OUT),MNHREAL_MPI,NP_EAST-1,1000+IS_WEST,&
                     NMNH_COMM_WORLD,IREQ(INB_REQ),IERR)
#ifdef MNH_GPUDIRECT
      !$acc end host_data
#endif
   ENDIF
END IF

IF (LY) THEN
   IF (.NOT.GSOUTH) THEN
#ifdef MNH_GPUDIRECT
      !$acc host_data use_device(PZSOUTH2F_OUT)
#endif
      INB_REQ = INB_REQ + 1
      CALL MPI_IRECV(PZSOUTH2F_OUT,SIZE(PZSOUTH2F_OUT),MNHREAL_MPI,NP_SOUTH-1,1000+IS_NORTH,&
                     NMNH_COMM_WORLD,IREQ(INB_REQ),IERR)
#ifdef MNH_GPUDIRECT
      !$acc end host_data
#endif
   ENDIF
   IF (.NOT.GNORTH) THEN
#ifdef MNH_GPUDIRECT
      !$acc host_data use_device(PZNORTH2F_OUT)
#endif
      INB_REQ = INB_REQ + 1
      CALL MPI_IRECV(PZNORTH2F_OUT,SIZE(PZNORTH2F_OUT),MNHREAL_MPI,NP_NORTH-1,1000+IS_SOUTH,&
                     NMNH_COMM_WORLD,IREQ(INB_REQ),IERR)
#ifdef MNH_GPUDIRECT
      !$acc end host_data
#endif
   ENDIF
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Copy the halo on the device PSRC to Zxxxx_IN

IF (LX) THEN
   IF (.NOT. GWEST) THEN
   !$acc kernels async(IS_WEST)
      PZWEST2F_IN ( KIIB:KIIB+KIHALO2_1  ,    KIJB:KIJE  , : )  = PSRC( KIIB:KIIB+KIHALO2_1  ,  KIJB:KIJE  , : )
!!$      PZWEST2F_IN ( : , : )  = PSRC( KIIB+1  , : , : )
   !$acc end kernels
      END IF
   IF (.NOT.GEAST) THEN
   !$acc kernels async(IS_EAST)
      PZEAST2F_IN ( KIIE-KIHALO2_1:KIIE  ,    KIJB:KIJE  , : )  = PSRC( KIIE-KIHALO2_1:KIIE  ,  KIJB:KIJE  , : )
!!$      PZEAST2F_IN ( : , : )  = PSRC( KIIE-1 ,  :  , : )
   !$acc end kernels
      ENDIF
END IF
IF (LY) THEN
   IF (.NOT.GSOUTH) THEN
   !$acc kernels async(IS_SOUTH)
   PZSOUTH2F_IN ( KIIB:KIIE  ,    KIJB:KIJB+KIHALO2_1  , : ) = PSRC( KIIB:KIIE  ,    KIJB:KIJB+KIHALO2_1  , : )
!!$      PZSOUTH2F_IN ( : , : ) = PSRC( : , KIJB+1 , : )
   !$acc end kernels
      ENDIF
   IF (.NOT.GNORTH) THEN
   !$acc kernels async(IS_NORTH)
      PZNORTH2F_IN ( KIIB:KIIE  ,    KIJE-KIHALO2_1:KIJE  , : ) = PSRC( KIIB:KIIE  ,    KIJE-KIHALO2_1:KIJE  , : )
!!$      PZNORTH2F_IN ( : , : ) = PSRC( : , KIJE-1  , : )      
   !$acc end kernels
   ENDIF
ENDIF
!$acc wait

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Send  Zxxxx2F_IN buffer via MPI(Gpu_direct) or copy to host
!
IF (LX) THEN
   IF (.NOT. GWEST) THEN
#ifdef MNH_GPUDIRECT
      !$acc host_data use_device(PZWEST2F_IN)
#else
      !$acc update host(PZWEST2F_IN)
#endif
      INB_REQ = INB_REQ + 1
      CALL MPI_ISEND(PZWEST2F_IN,SIZE(PZWEST2F_IN)  ,MNHREAL_MPI,NP_WEST-1,1000+IS_WEST,&
                     NMNH_COMM_WORLD,IREQ(INB_REQ),IERR)
#ifdef MNH_GPUDIRECT
      !$acc end host_data
#endif
   END IF
   IF (.NOT.GEAST) THEN
#ifdef MNH_GPUDIRECT
      !$acc host_data use_device(PZEAST2F_IN)
#else
      !$acc update host(PZEAST2F_IN)
#endif
      INB_REQ = INB_REQ + 1
      CALL MPI_ISEND(PZEAST2F_IN,SIZE(PZEAST2F_IN)  ,MNHREAL_MPI,NP_EAST-1,1000+IS_EAST,&
                     NMNH_COMM_WORLD,IREQ(INB_REQ),IERR)
#ifdef MNH_GPUDIRECT
      !$acc end host_data
#endif
   ENDIF
END IF

IF (LY) THEN
   IF (.NOT.GSOUTH) THEN
#ifdef MNH_GPUDIRECT
      !$acc host_data use_device(PZSOUTH2F_IN)
#else
      !$acc update host(PZSOUTH2F_IN)
#endif
      INB_REQ = INB_REQ + 1
      CALL MPI_ISEND(PZSOUTH2F_IN,SIZE(PZSOUTH2F_IN)  ,MNHREAL_MPI,NP_SOUTH-1,1000+IS_SOUTH,&
                     NMNH_COMM_WORLD,IREQ(INB_REQ),IERR)
#ifdef MNH_GPUDIRECT
      !$acc end host_data
#endif
   ENDIF
   IF (.NOT.GNORTH) THEN
#ifdef MNH_GPUDIRECT
      !$acc host_data use_device(PZNORTH2F_IN)
#else
      !$acc update host(PZNORTH2F_IN)
#endif
      INB_REQ = INB_REQ + 1
      CALL MPI_ISEND(PZNORTH2F_IN,SIZE(PZNORTH2F_IN)  ,MNHREAL_MPI,NP_NORTH-1,1000+IS_NORTH,&
                     NMNH_COMM_WORLD,IREQ(INB_REQ),IERR)
#ifdef MNH_GPUDIRECT
      !$acc end host_data
#endif
   ENDIF   
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF ( INB_REQ > 0 ) THEN
   CALL MPI_WAITALL(INB_REQ,IREQ,MNH_STATUSES_IGNORE,IERR)
END IF

!PW:disabled car KO si on compare avec 1 run 1 prc
! CALL MPPDB_CHECK(PSRC,'UPDATE_HALO_ll::'//TRIM(HNAME))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Update halo in PSRC + %HALO2 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


IF (LX) THEN
   IF (.NOT.GWEST) THEN
#ifndef MNH_GPUDIRECT
   !$acc update device(PZWEST2F_OUT) async(IS_WEST)
#endif
   !$acc kernels async(IS_WEST)
   PSRC( 1:KIIB-1  ,      KIJB:KIJE      , : ) = PZWEST2F_OUT( 1:KIIB-1  ,   KIJB:KIJE    , : )
   PZH2F_WEST( KIJB:KIJE , : ) = PZWEST2F_OUT( KIIB-2, KIJB:KIJE , : )      
   !$acc end kernels
   ENDIF
   IF (.NOT.GEAST) THEN
#ifndef MNH_GPUDIRECT
   !$acc update device(PZEAST2F_OUT) async(IS_EAST)
#endif
   !$acc kernels async(IS_EAST)
   PSRC( KIIE+1:KIIU  ,      KIJB:KIJE      , : ) = PZEAST2F_OUT( KIIE+1:KIIU  ,   KIJB:KIJE    , : )
   PZH2F_EAST( KIJB:KIJE , : ) = PZEAST2F_OUT( KIIE+2 , KIJB:KIJE , : )
   !$acc end kernels
   ENDIF
END IF
IF (LY) THEN
   IF (.NOT.GSOUTH) THEN
#ifndef MNH_GPUDIRECT
   !$acc update device(PZSOUTH2F_OUT) async(IS_SOUTH)
#endif
   !$acc kernels async(IS_SOUTH)
   PSRC(      KIIB:KIIE       ,  1:KIJB-1 , : ) = PZSOUTH2F_OUT(  KIIB:KIIE     , 1:KIJB-1  , : )
   PZH2F_SOUTH( KIIB:KIIE , : ) = PZSOUTH2F_OUT( KIIB:KIIE , KIJB-2 , : )  
   !$acc end kernels
   ENDIF
   IF (.NOT.GNORTH) THEN
#ifndef MNH_GPUDIRECT
   !$acc update device(PZNORTH2F_OUT) async(IS_NORTH)
#endif
   !$acc kernels async(IS_NORTH)
   PSRC(      KIIB:KIIE       , KIJE+1:KIJU , : ) = PZNORTH2F_OUT (  KIIB:KIIE     , KIJE+1:KIJU  , : )
   PZH2F_NORTH( KIIB:KIIE , : ) = PZNORTH2F_OUT ( KIIB:KIIE , KIJE+2 , : )      
   !$acc end kernels
   ENDIF
END IF
!$acc wait

!$acc end data

CALL MPPDB_CHECK(PSRC,"GET_HALO2_DF end:PSRC")
!CALL MPPDB_CHECK(PSRC,'UPDATE_HALO2_ll::'//TRIM(HNAME))

END SUBROUTINE GET_HALO2_DF_DIM

END SUBROUTINE GET_HALO2_DF
!
!     ###################################################
      SUBROUTINE GET_HALO2_D(PSRC, TP_PSRC_HALO2_ll, HNAME)
!     ###################################################
!
USE MODE_ll
USE MODD_ARGSLIST_ll, ONLY : LIST_ll, HALO2LIST_ll
USE MODI_GET_HALO, ONLY : GET_HALO_D,GET_HALO_DD,GET_HALO2_DD
USE MODE_MPPDB
!
IMPLICIT NONE
!
REAL, DIMENSION(:,:,:), INTENT(INOUT)  :: PSRC    ! variable at t
TYPE(HALO2LIST_ll), POINTER         :: TP_PSRC_HALO2_ll          ! halo2 for SRC
character(len=*), optional, intent(in) :: HNAME ! Name of the field to be added
!
character(len=:), allocatable    :: yname
INTEGER                          :: IIU,IJU,IKU            ! domain sizes
TYPE(LIST_ll)     , POINTER      :: TZ_PSRC_ll               ! halo
INTEGER                          :: IERROR                 ! error return code 
!
IIU = SIZE(PSRC,1)
IJU = SIZE(PSRC,2)
IKU = SIZE(PSRC,3)

if ( present ( hname ) ) then
  yname = hname
else
  yname = 'PSRC'
end if

CALL MPPDB_CHECK(PSRC,"GET_HALO2_D big:PSRC")

CALL GET_HALO_DD(PSRC,HNAME=yname)

!!$NULLIFY( TZ_PSRC_ll,TP_PSRC_HALO2_ll)
!!$CALL INIT_HALO2_ll(TP_PSRC_HALO2_ll,1,IIU,IJU,IKU)
!
CALL GET_HALO2_DD(PSRC,TP_PSRC_HALO2_ll,'GET_HALO2_DD::'//trim( yname ) )
!
!   clean local halo list , must be done outside
!
!!$CALL CLEANLIST_ll(TZ_PSRC_ll)
!
CALL MPPDB_CHECK(PSRC,"GET_HALO2_D end:PSRC")
!
END SUBROUTINE GET_HALO2_D
!
#endif
!-----------------------------------------------------------------------
!
!
!     ####################################
      SUBROUTINE DEL_HALO2_ll(TPHALO2LIST)
!     ####################################
!
!!****  *DEL_HALO2_ll* delete the second layer of the halo
!!
!!
!!    Purpose
!!    -------
!       The purpose of this routine is to deallocate the 
!     TPHALO2LIST variable which contains the second layer of the
!     halo for each variable.
! 
!!    Implicit Arguments
!!    ------------------
!     Module MODD_ARGSLIST_ll
!       type HALO2LIST_ll
!!
!!    Reference
!!    ---------
! 
!!    Author
!!    ------
!     J. Escobar                 * LA - CNRS *
!
!     Modification :
!     -------------
!     Juan  11/03/2010 : Memory Leak add DEALLOCATE(TZHALO2LIST%HALO2)     
! 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
  USE MODD_ARGSLIST_ll, ONLY : HALO2LIST_ll
!
  IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
  TYPE(HALO2LIST_ll), POINTER, INTENT(INOUT) :: TPHALO2LIST ! list of HALO2_lls
!
!
!*       0.2   Declarations of local variables :
!
  TYPE(HALO2LIST_ll), POINTER :: TZHALO2LIST
!
!-------------------------------------------------------------------------------
!
!*       1.    Deallocate the list of HALO2_lls
!
  TZHALO2LIST => TPHALO2LIST
!
  DO WHILE(ASSOCIATED(TZHALO2LIST))
!
   !Remark: OpenACC data create done in INIT_HALO2_ll
   !$acc exit data delete( TZHALO2LIST%HALO2%WEST, TZHALO2LIST%HALO2%EAST, TZHALO2LIST%HALO2%SOUTH, TZHALO2LIST%HALO2%NORTH )

   TPHALO2LIST => TZHALO2LIST%NEXT
    DEALLOCATE(TZHALO2LIST%HALO2%WEST)
    DEALLOCATE(TZHALO2LIST%HALO2%EAST)
    DEALLOCATE(TZHALO2LIST%HALO2%SOUTH)
    DEALLOCATE(TZHALO2LIST%HALO2%NORTH)
    DEALLOCATE(TZHALO2LIST%HALO2)
    DEALLOCATE(TZHALO2LIST)
    TZHALO2LIST => TPHALO2LIST
!
  ENDDO
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE DEL_HALO2_ll

MODULE MODE_HALO_MANAGED
  
  USE MODD_PRECISION,  ONLY: MNHINT64
  IMPLICIT NONE
  REAL, SAVE , ALLOCATABLE, DIMENSION(:) , TARGET :: ZHALO1DFLAT
  INTEGER(KIND=MNHINT64) , SAVE                   :: NHALO1DFLAT_MAXSIZE_R
  
  INTERFACE
     MODULE SUBROUTINE MNH_HALO_MANAGED_INIT()
     END SUBROUTINE MNH_HALO_MANAGED_INIT
  END INTERFACE
  INTERFACE
     MODULE SUBROUTINE MNH_HALO_MANAGED_GET( PBUFFER, KBUFFSIZE, KSENDNB_RECVNB )
       IMPLICIT NONE
       REAL, DIMENSION (:,:), POINTER, CONTIGUOUS :: PBUFFER
       INTEGER :: KBUFFSIZE, KSENDNB_RECVNB
     END SUBROUTINE MNH_HALO_MANAGED_GET
  END INTERFACE
  
END MODULE MODE_HALO_MANAGED

SUBMODULE (MODE_HALO_MANAGED) SMODE_HALO_MANAGED
  IMPLICIT NONE
CONTAINS
  
  MODULE SUBROUTINE MNH_HALO_MANAGED_INIT()

    USE MODE_MNH_ZWORK, ONLY : IIU,IJU,IKU,NPMAX_T1DFLAT_R
    USE MODD_CONF, ONLY      : NHALO
    IMPLICIT NONE
    
    NHALO1DFLAT_MAXSIZE_R = INT( 2*IIU*NHALO + 2*IJU*NHALO, KIND=MNHINT64 ) * IKU * NPMAX_T1DFLAT_R
    ALLOCATE (ZHALO1DFLAT(NHALO1DFLAT_MAXSIZE_R))
    !$acc enter data create (ZHALO1DFLAT)
    
  END SUBROUTINE MNH_HALO_MANAGED_INIT
  
  MODULE SUBROUTINE MNH_HALO_MANAGED_GET( PBUFFER, KBUFFSIZE, KSENDNB_RECVNB )
    IMPLICIT NONE
    REAL, DIMENSION (:,:), POINTER, CONTIGUOUS :: PBUFFER
    INTEGER :: KBUFFSIZE, KSENDNB_RECVNB
    
    PBUFFER(1:KBUFFSIZE, 1:KSENDNB_RECVNB ) =>  ZHALO1DFLAT(1:KBUFFSIZE*KSENDNB_RECVNB)
    
  END SUBROUTINE MNH_HALO_MANAGED_GET  
  
 END SUBMODULE SMODE_HALO_MANAGED      
