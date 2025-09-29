SUBROUTINE ARO_CONV_MNH(&
&KIDIA,KFDIA,&
&KLON,KLEV,&
&LDEEP,LSHALLOW,LDIAGCONV,&
&LSETTADJ,PTADJD,PTADJS,PDTCONV,&
&KSETENS,KICE,LREFRESH_ALL,LDOWN,PDXDY,&
&PAPRSF,PZZF,&
&PT,PRV,PRC,PRI,PRHOREF,&
&PU,PV,PW,&
&PTTEN,PRVTEN,PRCTEN,PRITEN,&
&PMF,&
&PPRTEN, PPRSTEN,&
&PPRLFLX,PPRSFLX,&
&PUMF,PDMF,&
&PCAPE,&
&KCLTOP,KCLBAS)
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
USE MODD_CONVPAR, ONLY : INI_CONVPAR
USE MODI_SHALLOW_CONVECTION, ONLY : SHALLOW_CONVECTION
USE MODI_DEEP_CONVECTION, ONLY: DEEP_CONVECTION

! Purpose:
! -------
!*** APL_CONV_MNH -MNH convection call (deep and shallow)
!*** ------------
 
! This routine is an interface routine to call the parametrisation
! of convection of the MNH physical package. It is a compilation of
! the CONVECTION routine of MNH and ACBECHT (E.Bazile and P. Bechtold)
! for ARPEGE/ALADIN. 
 

! Interface:
! ---------
! Input arguments
 
!-  Dimensions
 
!KIDIA,KFDIA : Start/end of horizontal loop
!KLON : Horizontal dimension (NPROMA in CPG)
!KLEV : Vertical dimension
 
!- Controls of callings
 
!LDEEP : Control the deep convection call (if true)
!LSHALLOW : Control the shallow convection call (if true)
!LDIAGCONV : Control the computation of some convective diagnostics (if true
!            cape, updraft mass flux, downdraft mass flux, liquid precip flux
!            and solid precip flux)
!LSETTADJ : logical to set convective adjustment time by user
!KSETENS : use of several calls to deep_conv (0 only one call)
!KICE : Ice cooking (0/1)
!LREFRESH_ALL : refresh or not tendencies at every call
!LDOWN : take or not convective downdrafts into account
 
!- 1D

!PTADJD ! user defined deep adjustement time (hard coded)
!PTADJS ! user defined shallow adjustement tile (hard coded)
!PDTCONV ! interval of time between 2 calls of the deep conv

!- 2D 
!
!PAPRSF : Pressure on full levels 
!PZZF : Altitude of full levels
!PT : Temperature
!PRV : Water vapor (rv)
!PRC : Cloud contain (rc)
!PRI : Ice (ri)
!PRHOREF : Rho
!PU : X-component of wind
!PV : Y-component of wind
!PW : Vertical velocity
!PDXDY : Grid size (m2)

! Output arguments
! ----------------
 
!-2D (1:KLEV)
 
!PTTEN : Convective T tendency (K/s)
!PRVTEN : Convective r_v tendency (1/s)
!PRCTEN : Convective r_c tendency (1/s)
!PRITEN : Convective r_i tendency (1/s)
!PMF : Convective mass flux for subgrid condensation
 
!-1D
!PPRTEN : Total surf precipitation tendency (m/s)
!PPRSTEN : Solid surf precipitation tendency (m/s)
 
!-Convective diagnostics (if required)
!  2D (1:KLEV)
!PPRLFLX : Liquid precip flux
!PPRSFLX : Solid precip flux
!PUMF : Updraft mass flux
!PDMF : Downdraft mass flux
!  1D
!PCAPE : Cape (JOULE)
!KCLTOP : Cloud top level (0 if no convection)
!KCLBAS : Cloud base level (0 if no convection)
 
 
! Method
! ------
! Adapted from the convection meso-nh routine (KFB Scheme)
! with the help of Eric Bazile previous work for Aladin.
 
! Externals
! ---------
! (mnh routines)
! INI_CONVPAR, INI_CONVPAR_E1, INI_CONVPAR_SHAL
! DEEP_CONVECTION, SHALLOW_CONVECTION 
 
! Author
! -----
! G. Hello : 03-10-08 
 
! Modifications
! -------------
!  E. Bazile : 2009/07/20 Local array TKECLS : input of shallow Conv. (KFB)
! End modifications
!--------------------------------------------------------------------
IMPLICIT NONE

!Input variables

INTEGER, INTENT(IN) :: KIDIA, KFDIA, KLON, KLEV

LOGICAL, INTENT(IN) :: LDEEP, LSHALLOW, LDIAGCONV 
LOGICAL, INTENT(IN) :: LSETTADJ
REAL, INTENT(IN) :: PTADJD,PTADJS,PDTCONV 
INTEGER, INTENT(IN) :: KSETENS, KICE
LOGICAL, INTENT(IN) :: LREFRESH_ALL, LDOWN
REAL, DIMENSION(KLON), INTENT(IN) :: PDXDY
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PAPRSF, PZZF 

REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PT, PRV, PRC, PRI, PRHOREF 
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PU, PV, PW


!Output variables

REAL, DIMENSION(KLON,KLEV), INTENT(OUT) :: PTTEN, PRVTEN, PRCTEN, PRITEN
REAL, DIMENSION(KLON,KLEV), INTENT(OUT) :: PMF
REAL, DIMENSION(KLON), INTENT(OUT) :: PPRTEN, PPRSTEN
REAL, DIMENSION(KLON,KLEV), INTENT(OUT) :: PPRLFLX, PPRSFLX, PUMF, PDMF
REAL, DIMENSION(KLON), INTENT(OUT) :: PCAPE
INTEGER, DIMENSION(KLON),INTENT(INOUT) :: KCLTOP,KCLBAS

!Local variables

LOGICAL :: OCHTRANS ! flag to compute convective transport
                        ! for chemical tracers
LOGICAL :: OCH_CONV_SCAV ! ?

INTEGER :: JLON,JLEV,JN ! loop index
INTEGER :: IIDIA, IFDIA ! horizontal loop bounds
INTEGER :: ILEV_MNH, IKU, IKB, IKE ! vertical loop bounds
INTEGER :: IBDIA, ITDIA, IPVEXT ! vertical loop bounds
INTEGER :: IENS ! number of calls to deep convection
INTEGER :: ICH1 ! number of chemical species (dimension)
INTEGER, DIMENSION(KLON) :: ICOUNT ! convective counter
INTEGER, DIMENSION(KLON) :: ICLTOP ! cloud top level
INTEGER, DIMENSION(KLON) :: ICLBAS ! cloud base level
INTEGER, DIMENSION(KLON) :: ICLTOPS ! cloud top level (shallow)
INTEGER, DIMENSION(KLON) :: ICLBASS ! cloud base level (shallow)
REAL, DIMENSION(KLON,KLEV):: ZAPRSF, ZZZF
REAL, DIMENSION(KLON,KLEV):: ZT, ZRV, ZRC, ZRI, ZRHOREF
REAL, DIMENSION(KLON,KLEV):: ZU, ZV, ZW, ZTEMP



REAL, DIMENSION(KLON) :: ZTIMEC ! time adjustement
REAL, DIMENSION(KLON) :: ZCAPE ! cape
REAL, DIMENSION(KLON) :: ZTKECLS ! TKE in CLS


!allocatable arrays

INTEGER, DIMENSION(:), ALLOCATABLE :: IEDUMMY ! field not to be recomputed by ensemble

REAL, DIMENSION(:,:), ALLOCATABLE :: ZTTEN,ZRVTEN,ZRCTEN,ZRITEN
REAL, DIMENSION(:,:), ALLOCATABLE :: ZPRLFLX, ZPRSFLX
REAL, DIMENSION(:,:), ALLOCATABLE :: ZUMF, ZDMF,ZUMFS
REAL, DIMENSION(:), ALLOCATABLE :: ZPRLTEN, ZPRSTEN
!chemical tracers
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZCH1, ZCH1TEN, ZCH1TENS 

! special for shallow convection
REAL, DIMENSION(:,:), ALLOCATABLE :: ZTTENS,ZRVTENS,ZRCTENS,ZRITENS
REAL, DIMENSION(:), ALLOCATABLE :: ZPRLTENS, ZPRSTENS


! Declarations of addiational ensemble fields 

REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZTTENE, ZRVTENE, ZRCTENE
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZRITENE, ZUMFE, ZDMFE
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZCH1TENE
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZPRLFLXE, ZPRSFLXE
REAL, DIMENSION(:,:), ALLOCATABLE :: ZPRLTENE, ZPRSTENE
REAL, DIMENSION(:), ALLOCATABLE :: ZEDUMMY, ZWEIGHT
REAL :: ZSUM


!1-INIT
!1-0 Dimension
!-------------
!a)levels
!Watch in MNH the vertical dimensionning is KLEV !!!
!there is one level underground and one level more upper level.
!Watch in MNH KLEV is upper atmosphere and 2 is ground.
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('ARO_CONV_MNH',0,ZHOOK_HANDLE)
ILEV_MNH=KLEV
IPVEXT=0 ! The MNH additionnal levels corresponds to JPVEXT in MODD_PARAMETERS
IKU=ILEV_MNH ! We mimic there the MNH convection routine
IKB=1+IPVEXT
IKE=IKU-IPVEXT
IBDIA=IKB
ITDIA=1+IPVEXT
! go to MNH level ordering
DO JLEV=1,KLEV
  ZAPRSF(:,JLEV)=PAPRSF(:,KLEV+1-JLEV)
  ZZZF(:,JLEV)=PZZF(:,KLEV+1-JLEV)
  ZT(:,JLEV)=PT(:,KLEV+1-JLEV)
  ZRV(:,JLEV)=PRV(:,KLEV+1-JLEV)
  ZRC(:,JLEV)=PRC(:,KLEV+1-JLEV)
  ZRI(:,JLEV)=PRI(:,KLEV+1-JLEV)
  ZRHOREF(:,JLEV)=PRHOREF(:,KLEV+1-JLEV)
  ZU(:,JLEV)=PU(:,KLEV+1-JLEV)
  ZV(:,JLEV)=PV(:,KLEV+1-JLEV)
  ZW(:,JLEV)=PW(:,KLEV+1-JLEV)
ENDDO

!b)horizontal
IIDIA= KIDIA! start index of horizontal convection computations
IFDIA=KFDIA ! end index of horizontal
! coming from arp world (IEST,IEND, NPROMA)
!1-1 Time 
!--------
IF(LSETTADJ) THEN 
    ZTIMEC(:)=PTADJD
ELSE
    ZTIMEC(:)=0.0
ENDIF

!1-2 Ensemble cooking
!--------------------
IENS=MIN(KSETENS,3)

!1-4 chemical tracers
!--------------------
!For the time being, things are not going to a common
!this is hard-coded.

OCHTRANS=.FALSE.
ICH1=1
OCH_CONV_SCAV=.FALSE.

!1-5 others
!----------
ICOUNT(:)=0 ! what to do with that ?
ZTKECLS(:)=0.0

!1-6 allocations
!---------------
ALLOCATE( ZTTEN(KLON,ILEV_MNH) )
ALLOCATE( ZRVTEN(KLON,ILEV_MNH) )
ALLOCATE( ZRCTEN(KLON,ILEV_MNH) )
ALLOCATE( ZRITEN(KLON,ILEV_MNH) )
ALLOCATE( ZUMF(KLON,ILEV_MNH) )
ALLOCATE( ZDMF(KLON,ILEV_MNH) )
ALLOCATE( ZPRLTEN(KLON) )
ALLOCATE( ZPRSTEN(KLON) )
ALLOCATE( ZPRLFLX(KLON,ILEV_MNH) )
ALLOCATE( ZPRSFLX(KLON,ILEV_MNH) )
ALLOCATE( ZCH1(KLON,ILEV_MNH,ICH1) )
ALLOCATE( ZCH1TEN(KLON,ILEV_MNH,ICH1) )

ALLOCATE( ZTTENS(KLON,ILEV_MNH) )
ALLOCATE( ZRVTENS(KLON,ILEV_MNH) )
ALLOCATE( ZRCTENS(KLON,ILEV_MNH) )
ALLOCATE( ZRITENS(KLON,ILEV_MNH) )
ALLOCATE( ZUMFS(KLON,ILEV_MNH) )
ALLOCATE( ZPRLTENS(KLON) )
ALLOCATE( ZPRSTENS(KLON) )
ALLOCATE( ZCH1TENS(KLON,ILEV_MNH,ICH1) )

!* Allocation of additional ensemble members
IF ( IENS > 0 ) THEN
     ALLOCATE( ZTTENE(KLON,ILEV_MNH,IENS) )
     ALLOCATE( ZRVTENE(KLON,ILEV_MNH,IENS) )
     ALLOCATE( ZRCTENE(KLON,ILEV_MNH,IENS) )
     ALLOCATE( ZRITENE(KLON,ILEV_MNH,IENS) )
     ALLOCATE( ZUMFE(KLON,ILEV_MNH,IENS) )
     ALLOCATE( ZDMFE(KLON,ILEV_MNH,IENS) )
     ALLOCATE( ZCH1TENE(KLON,ILEV_MNH,ICH1,IENS) )
     ALLOCATE( ZPRLFLXE(KLON,ILEV_MNH,IENS) )
     ALLOCATE( ZPRSFLXE(KLON,ILEV_MNH,IENS) )
     ALLOCATE( ZPRLTENE(KLON,IENS) )
     ALLOCATE( ZPRSTENE(KLON,IENS) )
     ALLOCATE( ZEDUMMY(KLON) )
     ALLOCATE( IEDUMMY(KLON) )
     ALLOCATE( ZWEIGHT(IENS) )
     IEDUMMY(:) = 0
     ZCH1TENE(:,:,:,:)=0.0
ENDIF

!2-Interface to MNH world
!------------------------


! Initializing tendencies !?!?!?!?

  PTTEN(:,:) =0.0
  PRVTEN(:,:)= 0.0
  PRCTEN(:,:)=0.0            
  PRITEN(:,:)=0.0


!chemical species (all put to zero)
ZCH1TEN(:,:,:)=0.0
ZCH1TENS(:,:,:)=0.0
ZCH1(:,:,:)=0.0


!3-Deep convection call
!----------------------
IF(LDEEP) THEN

! 3.1. Base Version

  CALL INI_CONVPAR

  CALL ABOR1('FIXME: THE INTERFACE IS WRONG')
  !CALL DEEP_CONVECTION( KLON, ILEV_MNH, IIDIA, IFDIA, IBDIA, ITDIA,&
  !                      &PDTCONV, KICE, LREFRESH_ALL, LDOWN, LSETTADJ,& 
  !                      &PAPRSF, PZZF,&
  !                      &PDXDY,ZTIMEC,&
  !                      &PT, PRV,& 
  !                      &PRC, PRI,&
  !                      &PU,PV,&
  !                      &PW,&
  !                      &ICOUNT, ZTTEN,&
  !                      &ZRVTEN, ZRCTEN,&
  !                      & ZRITEN,&
  !                      &ZPRLTEN, ZPRSTEN,&
  !                      &ICLTOP, ICLBAS,&
  !                      & ZPRLFLX, ZPRSFLX,& 
  !                      &ZUMF, ZDMF,&
  !                      &ZCAPE,&
  !                      &OCHTRANS, ICH1, ZCH1,&
  !                      & ZCH1TEN,&
  !                      &OCH_CONV_SCAV, PRHOREF)
!  3.2. Additional Ensemble members
 
    IF ( IENS > 0 ) THEN
 
    CALL INI_CONVPAR_E1
 
!* first member - changes in MODD_CONVPAR (cloud radius of 500 m)
 
     CALL ABOR1('FIXME: THE INTERFACE IS WRONG')
     !CALL DEEP_CONVECTION( KLON, ILEV_MNH, IIDIA, IFDIA, IBDIA, ITDIA,&
     !                      &PDTCONV, KICE, LREFRESH_ALL, LDOWN, LSETTADJ,&
     !                      &PAPRSF,  PZZF, PDXDY, ZTIMEC,&
 
     !                     &PT, PRV, PRC, PRI, PU, PV, PW,&
     !                     &IEDUMMY, ZTTENE(:,:,1), ZRVTENE(:,:,1), ZRCTENE(:,:,1), ZRITENE(:,:,1),&
     !                     &ZPRLTENE(:,1), ZPRSTENE(:,1),&
     !                     &IEDUMMY, IEDUMMY, ZPRLFLXE(:,:,1), ZPRSFLXE(:,:,1),&
     !                     &ZUMFE(:,:,1), ZDMFE(:,:,1), ZEDUMMY,&
     !                     &OCHTRANS, ICH1, ZCH1, ZCH1TENE(:,:,:,1),&
     !                     &OCH_CONV_SCAV, PRHOREF)                  
    ENDIF
 
    IF (  IENS > 1 ) THEN
 
    CALL INI_CONVPAR
 
!* second member (positive vertical velocity perturb for Trigger)
 
    CALL ABOR1('FIXME: THE INTERFACE IS WRONG')
    !CALL DEEP_CONVECTION( KLON, ILEV_MNH, IIDIA, IFDIA, IBDIA, ITDIA,&
    !                      &PDTCONV, KICE, LREFRESH_ALL, LDOWN, LSETTADJ,&
    !                      &PAPRSF, PZZF, PDXDY, ZTIMEC,&
    !                      &PT, PRV, PRC, PRI, PU, PV, PW*1.5+1.E-4,&
    !                      &IEDUMMY, ZTTENE(:,:,2), ZRVTENE(:,:,2), ZRCTENE(:,:,2), ZRITENE(:,:,2),&
    !                      &ZPRLTENE(:,2), ZPRSTENE(:,2),&
    !                      &IEDUMMY, IEDUMMY, ZPRLFLXE(:,:,2), ZPRSFLXE(:,:,2),&
    !                      &ZUMFE(:,:,2), ZDMFE(:,:,2), ZEDUMMY,&
    !                      &OCHTRANS, ICH1, ZCH1, ZCH1TENE(:,:,:,2),&
    !                      &OCH_CONV_SCAV, PRHOREF)                  
    ENDIF
!
    IF ( IENS > 2 ) THEN
 
!* third member (positive vertical velocity perturb for Trigger)
 
     CALL ABOR1('FIXME: THE INTERFACE IS WRONG')
     !CALL DEEP_CONVECTION( KLON, ILEV_MNH, IIDIA, IFDIA, IBDIA, ITDIA,&
     !                     &PDTCONV, KICE, LREFRESH_ALL, LDOWN, LSETTADJ,&
     !                     &PAPRSF, PZZF, PDXDY, ZTIMEC,&
     !                     &PT, PRV, PRC, PRI, PU, PV, PW*.5-1.E-4,&
     !                     &IEDUMMY, ZTTENE(:,:,3), ZRVTENE(:,:,3), ZRCTENE(:,:,3), ZRITENE(:,:,3),&
     !                     &ZPRLTENE(:,3), ZPRSTENE(:,3),&
     !                     &IEDUMMY, IEDUMMY, ZPRLFLXE(:,:,3), ZPRSFLXE(:,:,3),&
     !                     &ZUMFE(:,:,3), ZDMFE(:,:,3), ZEDUMMY,&
     !                     &OCHTRANS, ICH1, ZCH1, ZCH1TENE(:,:,:,3),&
     !                     &OCH_CONV_SCAV, PRHOREF)                  
    ENDIF
 
ENDIF
IF ( .NOT. LDEEP ) THEN
  ICOUNT(:)=0
  ZTTEN(:,:)=0.0
  ZRVTEN(:,:)=0.0
  ZRCTEN(:,:)=0.0
  ZRITEN(:,:)=0.0
  ZUMF(:,:)=0.
  ZDMF(:,:)=0.
  ICLTOP(:)=0
  ICLBAS(:)=0
  ZCH1TEN(:,:,:)=0.0
  ZPRLTEN(:)=0.0
  ZPRSTEN(:)=0.0
  ZPRLFLX(:,:)=0.0
  ZPRSFLX(:,:)=0.0
  ZCAPE(:)=0.0
ENDIF

!4-Call shallow convective routine
!---------------------------------

IF ( LSHALLOW ) THEN
!
  IF ( .NOT. LDEEP ) CALL INI_CONVPAR 
  CALL INI_CONVPAR_SHAL
!
CALL ABOR1('FIXME:THE INTERFACE IS WRONG')
  !CALL SHALLOW_CONVECTION( KLON, ILEV_MNH, IIDIA, IFDIA, IBDIA, ITDIA,&
                           !PDTCONV, KICE, LSETTADJ, PTADJS,           &
                           !PAPRSF, PZZF, ZTKECLS,                     &
                           !PT, PRV, PRC, PRI, PW,                     &
                           !ZTTENS, ZRVTENS, ZRCTENS, ZRITENS,         &
                           !ICLTOPS, ICLBASS, ZUMFS,                   &
                           !OCHTRANS, ICH1, ZCH1, ZCH1TENS             &
                          !)
ENDIF
IF ( .NOT. LSHALLOW ) THEN

  ZTTENS(:,:)=0.0
  ZRVTENS(:,:)=0.0
  ZRCTENS(:,:)=0.0
  ZRITENS(:,:)=0.0
  ZUMFS(:,:)=0.
  ICLTOPS(:)=0
  ICLBASS(:)=0
  ZCH1TENS(:,:,:)=0.0
END IF

!-------------------------------------------------------------------------
!5-Add deep and shallow tendencies, and - if activated - ensemble members
!-------------------------------------------------------------------------

!5.1-Ensemble members for deep
!-----------------------------
ZSUM = 1.
IF ( IENS > 0 ) THEN
    IF ( IENS == 1 ) ZWEIGHT(:) = .5
    IF ( IENS >  1 ) ZWEIGHT(:) = 1.
    DO JN = 1, IENS
       ZTTEN(:,:)  = ZTTEN(:,:)  + ZWEIGHT(JN) * ZTTENE(:,:,JN)
       ZRVTEN(:,:) = ZRVTEN(:,:) + ZWEIGHT(JN) * ZRVTENE(:,:,JN)
       ZRCTEN(:,:) = ZRCTEN(:,:) + ZWEIGHT(JN) * ZRCTENE(:,:,JN)
       ZRITEN(:,:) = ZRITEN(:,:) + ZWEIGHT(JN) * ZRITENE(:,:,JN)
       ZPRLFLX(:,:)= ZPRLFLX(:,:)+ ZWEIGHT(JN) * ZPRLFLXE(:,:,JN)
       ZPRSFLX(:,:)= ZPRSFLX(:,:)+ ZWEIGHT(JN) * ZPRSFLXE(:,:,JN)
       ZUMF(:,:)   = ZUMF(:,:)   + ZWEIGHT(JN) * ZUMFE(:,:,JN)
       ZDMF(:,:)   = ZDMF(:,:)   + ZWEIGHT(JN) * ZDMFE(:,:,JN)
       ZPRLTEN(:)  = ZPRLTEN(:)  + ZWEIGHT(JN) * ZPRLTENE(:,JN)
       ZPRSTEN(:)  = ZPRSTEN(:)  + ZWEIGHT(JN) * ZPRSTENE(:,JN)
       IF ( OCHTRANS )  &
         & ZCH1TEN(:,:,:) = ZCH1TEN(:,:,:) + ZWEIGHT(JN) * ZCH1TENE(:,:,:,JN)
    ENDDO

    ZSUM = 1. / ( 1. + SUM( ZWEIGHT(:) ) )
ENDIF

!5.2- Add deep and shallow
!-------------------------
 DO JLON = IIDIA,IFDIA
    DO JLEV = 1,KLEV
      PTTEN(JLON,JLEV) = ZTTEN(JLON,JLEV) * ZSUM + ZTTENS(JLON,JLEV)
      PRVTEN(JLON,JLEV) = ZRVTEN(JLON,JLEV)  * ZSUM + ZRVTENS(JLON,JLEV)
      PRCTEN(JLON,JLEV) = ZRCTEN(JLON,JLEV)  * ZSUM + ZRCTENS(JLON,JLEV)
      PRITEN(JLON,JLEV) = ZRITEN(JLON,JLEV)  * ZSUM + ZRITENS(JLON,JLEV)
    ENDDO
ENDDO

DO JLON =  IIDIA,IFDIA
     PPRTEN(JLON)  = (ZPRLTEN(JLON) + ZPRSTEN(JLON) ) * ZSUM 
     PPRSTEN(JLON) = ZPRSTEN(JLON) * ZSUM
ENDDO
!IF ( OCHTRANS ) THEN
!  DO JLON = IIDIA,IFDIA
!    DO JLEV = 1,KLEV
!        PCH1TEN(JLON,JLEV,:) = ZCH1TEN(JLON,JLEV,:) * ZSUM + ZCH1TENS(JLON,JLEV,:)
!    ENDDO
!  ENDDO
!ENDIF
!
!5.3-convective mass flux for subgrid condensation
!-------------------------------------------------
IF(SIZE(PMF) .NE. 0) THEN
  DO JLON = IIDIA,IFDIA
      DO JLEV = 1,KLEV
          PMF(JLON,JLEV)   = (ZUMF(JLON,JLEV) + ZDMF(JLON,JLEV))*ZSUM + ZUMFS(JLON,JLEV) 
      ENDDO
   ENDDO
ELSE
  PMF(:,:)=0.0
ENDIF
!
!5.4-Diagnostics
!--------------------
IF( LDIAGCONV ) THEN
  DO JLON = IIDIA,IFDIA
      DO JLEV = 1,KLEV 
        PPRLFLX(JLON,JLEV)= ZPRLFLX(JLON,JLEV) * ZSUM
        PPRSFLX(JLON,JLEV)= ZPRSFLX(JLON,JLEV) * ZSUM
        PUMF(JLON,JLEV)   = ZUMF(JLON,JLEV) * ZSUM + ZUMFS(JLON,JLEV)
        PDMF(JLON,JLEV)   = ZDMF(JLON,JLEV) * ZSUM
      ENDDO
   ENDDO

  DO JLON = IIDIA,IFDIA
      PCAPE(JLON)   = ZCAPE(JLON)
      KCLTOP(JLON)  = MAX(ICLTOP(JLON), ICLTOPS(JLON))
      KCLBAS(JLON)  = MAX(ICLBAS(JLON), ICLBASS(JLON))
  ENDDO
ELSE
   PPRLFLX(:,:)=0.0
   PPRSFLX(:,:)=0.0
   PUMF(:,:)=0.0
   PDMF(:,:)=0.0
   PCAPE(:)=0.0
   KCLTOP(:)=0
   KCLBAS(:)=0
ENDIF

!--------------------
!   6.  Deallocation 
!--------------------

!* additional ensemble members
!
IF ( IENS > 0 ) THEN
     DEALLOCATE( ZTTENE )
     DEALLOCATE( ZRVTENE )
     DEALLOCATE( ZRCTENE )
     DEALLOCATE( ZRITENE )
     DEALLOCATE( ZUMFE )
     DEALLOCATE( ZDMFE )
     DEALLOCATE( ZCH1TENE )
     DEALLOCATE( ZPRLFLXE )
     DEALLOCATE( ZPRSFLXE )
     DEALLOCATE( ZPRLTENE )
     DEALLOCATE( ZPRSTENE )
     DEALLOCATE( ZEDUMMY )
     DEALLOCATE( IEDUMMY )
     DEALLOCATE( ZWEIGHT )
ENDIF

!* local arrays
 
DEALLOCATE(ZTTEN)
DEALLOCATE(ZRVTEN)
DEALLOCATE(ZRCTEN)
DEALLOCATE(ZRITEN)
DEALLOCATE(ZPRLTEN)
DEALLOCATE(ZPRSTEN)
DEALLOCATE(ZPRLFLX)
DEALLOCATE(ZPRSFLX)
DEALLOCATE(ZUMF)
DEALLOCATE(ZDMF)
DEALLOCATE(ZCH1)
DEALLOCATE(ZCH1TEN)

DEALLOCATE(ZTTENS)
DEALLOCATE(ZRVTENS)
DEALLOCATE(ZRCTENS)
DEALLOCATE(ZRITENS)
DEALLOCATE(ZUMFS)
DEALLOCATE(ZPRLTENS)
DEALLOCATE(ZPRSTENS)
DEALLOCATE(ZCH1TENS)

! go back to ARPEGE level ordering
DO JLEV=1,KLEV
  ZTEMP(:,JLEV)=PTTEN(:,KLEV+1-JLEV)
ENDDO
PTTEN(:,:)=ZTEMP(:,:)
DO JLEV=1,KLEV
  ZTEMP(:,JLEV)=PRVTEN(:,KLEV+1-JLEV)
ENDDO
PRVTEN(:,:)=ZTEMP(:,:)
DO JLEV=1,KLEV
  ZTEMP(:,JLEV)=PRCTEN(:,KLEV+1-JLEV)
ENDDO
PRCTEN(:,:)=ZTEMP(:,:)
DO JLEV=1,KLEV
  ZTEMP(:,JLEV)=PRITEN(:,KLEV+1-JLEV)
ENDDO
PRITEN(:,:)=ZTEMP(:,:)
DO JLEV=1,KLEV
  ZTEMP(:,JLEV)=PMF(:,KLEV+1-JLEV)
ENDDO
PMF(:,:)=ZTEMP(:,:)
DO JLEV=1,KLEV
  ZTEMP(:,JLEV)=PPRLFLX(:,KLEV+1-JLEV)
ENDDO
PPRLFLX(:,:)=ZTEMP(:,:)
DO JLEV=1,KLEV
  ZTEMP(:,JLEV)=PPRSFLX(:,KLEV+1-JLEV)
ENDDO
PPRSFLX(:,:)=ZTEMP(:,:)
DO JLEV=1,KLEV
  ZTEMP(:,JLEV)=PUMF(:,KLEV+1-JLEV)
ENDDO
PUMF(:,:)=ZTEMP(:,:)
DO JLEV=1,KLEV
  ZTEMP(:,JLEV)=PDMF(:,KLEV+1-JLEV)
ENDDO
PDMF(:,:)=ZTEMP(:,:)



IF (LHOOK) CALL DR_HOOK('ARO_CONV_MNH',1,ZHOOK_HANDLE)
END SUBROUTINE ARO_CONV_MNH
