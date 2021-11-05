!MNH_LIC Copyright 2012-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     #######################
      MODULE MODI_PROGNOS_LIMA
!     #######################
!
INTERFACE
!
SUBROUTINE PROGNOS_LIMA(PTSTEP,PDZ,PLV,PCPH,PPRES,PRHOD,PRR,PTT,PRV,PRC,PS0,PNAS,PCCS,PNFS)
!
REAL,                     INTENT(IN)    :: PTSTEP
REAL, DIMENSION(:),       INTENT(IN)    :: PPRES
REAL, DIMENSION(:),       INTENT(IN)    :: PDZ
REAL, DIMENSION(:),       INTENT(IN)    :: PLV
REAL, DIMENSION(:),       INTENT(IN)    :: PCPH
REAL, DIMENSION(:),       INTENT(IN)    :: PRHOD
REAL, DIMENSION(:),       INTENT(IN)    :: PRR
REAL, DIMENSION(:),       INTENT(INOUT) :: PTT ! PTHS
REAL, DIMENSION(:),       INTENT(INOUT) :: PRV ! PRVS
REAL, DIMENSION(:),       INTENT(INOUT) :: PRC ! PRCS
REAL, DIMENSION(:),       INTENT(INOUT) :: PS0 ! PSVS sursat source
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PNAS ! PSVS activated aerosols source
REAL, DIMENSION(:),       INTENT(INOUT) :: PCCS ! PSVS droplet concentration source
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PNFS ! PSVS free aerosol source           
!
END SUBROUTINE PROGNOS_LIMA
!
END INTERFACE
!
END MODULE MODI_PROGNOS_LIMA
!
!     ###################################################################################
      SUBROUTINE PROGNOS_LIMA(PTSTEP,PDZ,PLV,PCPH,PPRES,PRHOD,PRR,PTT,PRV,PRC,PS0,PNAS,PCCS,PNFS)
!     ###################################################################################
!
!!****  * -  compute pseudo-prognostic of supersaturation according to Thouron
!                                                                     et al. 2012
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!
!!    REFERENCE
!!    ---------
!!
!!      Thouron, O., J.-L. Brenguier, and F. Burnet, Supersaturation calculation
!!      in large eddy simulation models for prediction of the droplet number
!!      concentration, Geosci. Model Dev., 5, 761-772, 2012.
!!
!!    AUTHOR
!!    ------
!!     06/2021 B. Vie forked from prognos.f90 
!!
!!    MODIFICATIONS
!!    -------------
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
USE MODD_CST
USE MODD_PARAM_LIMA
USE MODD_PARAM_LIMA_WARM
!
USE MODE_IO
USE MODE_MSG
!
USE MODI_GAMMA
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
!
REAL,                     INTENT(IN)    :: PTSTEP
REAL, DIMENSION(:),       INTENT(IN)    :: PPRES
REAL, DIMENSION(:),       INTENT(IN)    :: PDZ
REAL, DIMENSION(:),       INTENT(IN)    :: PLV
REAL, DIMENSION(:),       INTENT(IN)    :: PCPH
REAL, DIMENSION(:),       INTENT(IN)    :: PRHOD
REAL, DIMENSION(:),       INTENT(IN)    :: PRR
REAL, DIMENSION(:),       INTENT(INOUT) :: PTT ! PTHS
REAL, DIMENSION(:),       INTENT(INOUT) :: PRV ! PRVS
REAL, DIMENSION(:),       INTENT(INOUT) :: PRC ! PRCS
REAL, DIMENSION(:),       INTENT(INOUT) :: PS0 ! PSVS sursat source
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PNAS ! PSVS activated aerosols source
REAL, DIMENSION(:),       INTENT(INOUT) :: PCCS ! PSVS droplet concentration source
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PNFS ! PSVS free aerosol source           
!
!
!*       0.2   Declarations of local variables :
!
!
REAL, DIMENSION(SIZE(PRHOD,1))   ::  ZW1,ZW2,ZDZRC2,ZDZRC,ZCPH
REAL, DIMENSION(SIZE(PRHOD,1))   :: ZA1,ZA2,ZB,ZC,ZG
REAL, DIMENSION(SIZE(PRHOD,1))   :: ZLV,ZTT1,ZRT,ZTL,ZTT1_TEMP,ZTT2_TEMP
REAL, DIMENSION(SIZE(PRHOD,1))   :: ZRMOY,ZRVSAT1,ZRVSAT2
REAL, DIMENSION(SIZE(PRHOD,1))   :: ZVEC2  ! Work vectors forinterpolations
INTEGER, DIMENSION(SIZE(PRHOD,1)):: IVEC2   ! Vectors of indices for interpolations
INTEGER :: J1,J2,JMOD,INUCT,JL
REAL,DIMENSION(SIZE(PS0,1))      ::MEM_PS0,ADJU2
REAL::AER_RAD
REAL, DIMENSION(SIZE(PRHOD,1))   :: ZFLAG_ACT   !Flag for activation
!
INTEGER                          :: IRESP      ! Return code of FM routines
INTEGER                          :: ILUOUT     ! Logical unit of output listing
CHARACTER(LEN=100)               :: YMSG
!
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZCHEN_MULTI,ZTMP
REAL, DIMENSION(:), ALLOCATABLE    :: ZZW1, ZZW2, ZZW6, ZVEC1
INTEGER, DIMENSION(:), ALLOCATABLE :: IVEC1             ! Vectors of indices for
                                                        ! interpolations

!
INUCT = SIZE(PTT,1)
!
!
 ALLOCATE(ZZW1(INUCT))
 ALLOCATE(ZZW2(INUCT))
 ALLOCATE(ZZW6(INUCT))
 ALLOCATE(ZCHEN_MULTI(INUCT,NMOD_CCN))
 ALLOCATE(ZTMP(INUCT,NMOD_CCN))
 ALLOCATE(ZVEC1(INUCT))
 ALLOCATE(IVEC1(INUCT))
!
!
 DO JL=1,INUCT
  DO JMOD = 1,NMOD_CCN
   ZCHEN_MULTI(JL,JMOD) = (PNFS(JL,JMOD)+PNAS(JL,JMOD))*PRHOD(JL)  &
                                                / XLIMIT_FACTOR(JMOD)
  ENDDO
 END DO
!print*,'ZCHEN_MULTI=',MINVAL(ZCHEN_MULTI(:,1)), MAXVAL(ZCHEN_MULTI(:,1)), &
!       'ZCHEN_MULTI(1,1)=',ZCHEN_MULTI(1,1)
!
!*       . Compute the nucleus source
!   	 -----------------------------
!
!
! Modified values for Beta and C (see in init_aerosol_properties) account for that
!
   WHERE ( PS0(:) > 0.)
      ZVEC1(:) = MAX( 1.0001, MIN( REAL(NHYP)-0.0001,  &
                                    XHYPINTP1*LOG(PS0(:))+XHYPINTP2 ) )
      IVEC1(:) = INT( ZVEC1(:) )
      ZVEC1(:) = ZVEC1(:) - REAL( IVEC1(:) )
   END WHERE
!print*,'ZVEC1=',MINVAL(ZVEC1), MAXVAL(ZVEC1)
   ZZW6(:)  = 0. ! initialize the change of cloud droplet concentration
!
   ZTMP(:,:)=0.0
!
! Compute the concentration of activable aerosols for each mode
! based on the supersaturation ( -> ZTMP )
!
   DO JMOD = 1, NMOD_CCN                     ! iteration on mode number
      ZZW1(:) = 0.
   !
      WHERE( PS0(:)>0.0 )
         ZZW1(:) =  XHYPF12( IVEC1(:)+1,JMOD )* ZVEC1(:)      & ! hypergeo function
                  - XHYPF12( IVEC1(:)  ,JMOD )*(ZVEC1(:) - 1.0) ! XHYPF12 is tabulated
   !
         ZTMP(:,JMOD) = (ZCHEN_MULTI(:,JMOD)/PRHOD(:))*PS0(:)**XKHEN_MULTI(JMOD) &
                                                         *ZZW1(:)
   !     ZTMP(:,JMOD) = (ZCHEN_MULTI(:,JMOD)/PRHOD(:))*100*PS0(:)**XKHEN_MULTI(JMOD) &
      ENDWHERE
!print*,'ZZW1=',MINVAL(ZZW1), MAXVAL(ZZW1)
!print*,'ZTMP=',MINVAL(ZTMP), MAXVAL(ZTMP)
   ENDDO
!
! Compute the concentration of aerosols activated at this time step
! as the difference between ZTMP and the aerosols already activated at t-dt (ZZW1)
!
   DO JMOD = 1, NMOD_CCN                     ! iteration on mode number
      ZZW2(:) = 0.
   !
!     WHERE( SUM(ZTMP(:,:),DIM=2)*PTSTEP .GT. 15.E6/PRHOD(:) ) 
         ZZW2(:) = MIN( PNFS(:,JMOD),MAX( ZTMP(:,JMOD)- PNAS(:,JMOD) , 0.0 ) )
!     ENDWHERE
!print*,'ZTMP=',ZTMP(:,1)
!print*,'PNAS=',PNAS(:,1)
!print*,'PNFS=',PNFS(:,1)
!print*,'ZZW2=',ZZW2(:)
   !
   !* update the concentration of activated CCN = Na
   !
      PNAS(:,JMOD) = (PNAS(:,JMOD) +  ZZW2(:))
   !
   !* update the concentration of free CCN = Nf
   !
      PNFS(:,JMOD) = (PNFS(:,JMOD) -  ZZW2(:)) 
   !
   !* prepare to update the cloud water concentration 
   !
      ZZW6(:) = ZZW6(:) + ZZW2(:)
!print*,'ZZW6=',MINVAL(ZZW6), MAXVAL(ZZW6)
   ENDDO
!
!FLAG ACTIVE A TRUE (1.0) si on active pas
ZFLAG_ACT(:)=0.0
DO J2=1,SIZE(PRC,1)
 IF (ZZW2(J2).EQ.0.0) THEN
 ZFLAG_ACT(J2)=1.0
 ENDIF
!print*,'ZFLAG_ACT=',ZFLAG_ACT(J2)
ENDDO
!
! Mean radius
!minimum radius of cloud droplet
AER_RAD=1.0E-6
ZRMOY(:)=0.0
DO J2=1,SIZE(PRC,1)
 IF ((PRC(J2).NE.0.) .AND. (PCCS(J2).NE.0.)) THEN
  ZRMOY(J2)=(MOMG(XALPHAC,XNUC,3.0)*4.0*XPI*PCCS(J2)*XRHOLW/&
        (3.0*PRC(J2)*PRHOD(J2)))**(1.0/3.0)
  ZRMOY(J2)=(PCCS(J2)*MOMG(XALPHAC,XNUC,1.0)/ZRMOY(J2))
 ENDIF
!ZRMOY(J2)=ZRMOY(J2)+(ZZW2(J2)*AER_RAD)
 ZRMOY(J2)=ZRMOY(J2)+(ZZW6(J2)*AER_RAD)
ENDDO
   !print*,'prognos RMOY=',MINVAL(ZRMOY),MAXVAL(ZRMOY)
!
!  PCCS(:) = ZZW6(:) * PTSTEP
   PCCS(:) = PCCS(:) + ZZW6(:) 
   !print*,'prognos PCCS=',MINVAL(PCCS),MAXVAL(PCCS)
!
!CALCUL DE A1 => Estimation de (drs/dt)f
!T(=à determiner) avant forcage; T'(=PTT) apres forcage
!Calcul de ZTT1: calculé en inversant S0(T)jusqu'à T:
! l'erreur faite sur cette inversion est supérieur à la précision
! recherchée, on applique à rs(T') pour cxalculer le DT=T'-T qui
! correspond à la variation rs(T')-rs(T). Permet de recuperer une valeur
! correcte de DT et donc de determiner T comme T=T'-DT         
!ZRVSAT1=rs(T)
!
!print*,'prognos : PS0=',MINVAL(PS0),MAXVAL(PS0)
ZRVSAT1(:)=PRV(:)/(PS0(:)+1.0)
!ZTT1<--es(T) de rs(T)
ZTT1_TEMP(:)=PPRES(:)*((((XMV / XMD)/ZRVSAT1(:))+1.0)**(-1D0))
!ZTT1<--T de es(T)
ZTT1_TEMP(:)=LOG(ZTT1_TEMP(:)/610.8)
ZTT1_TEMP(:)=(31.25*ZTT1_TEMP(:) -17.5688*273.15)/(ZTT1_TEMP(:) - 17.5688)
!es(T')
ZW1(:)=EXP(XALPW-XBETAW/PTT(:)-XGAMW*LOG(PTT(:)))
!ZRVSAT2=rs(T')
ZRVSAT2(:)=(XMV / XMD)*ZW1(:)/(PPRES(:)-ZW1(:))
!ZTT2<--es(T') de rs(T')
ZTT2_TEMP(:)=PPRES(:)*((((XMV / XMD)/ZRVSAT2(:))+1.0)**(-1D0))
!ZTT2<--T' de es(T')
IF (MINVAL(ZTT2_TEMP).LT.0.0) THEN
  WRITE(YMSG,*) 'ZTT2_TEMP',MINVAL(ZTT2_TEMP),MINLOC(ZTT2_TEMP)
  CALL PRINT_MSG(NVERB_FATAL,'GEN','PROGNOS_LIMA',YMSG)
ENDIF
!
ZTT2_TEMP(:)=LOG(ZW1(:)/610.8)
ZTT2_TEMP(:)=(31.25*ZTT2_TEMP(:) -17.5688*273.15)/(ZTT2_TEMP(:) - 17.5688)
!ZTT1=T'-DT
ZTT1(:)=PTT(:)-(ZTT2_TEMP(:)-ZTT1_TEMP(:))
!Lv(T)
ZLV(:) = XLVTT+(XCPV-XCL)*(ZTT1(:)-XTT)                
!
ZA1(:)=-(((PS0(:)+1.0)**2.0)/PRV(:))*(ZRVSAT2(:)-(PRV(:)/(PS0(:)+1.0)))/PTSTEP
!G
ZG(:)= 1.0/(XRHOLW*((XRV*ZTT1(:)/(XDIVA*EXP(XALPW-(XBETAW/ZTT1(:))-(XGAMW*LOG(ZTT1(:))))))     &
+((ZLV(:)/(XTHCO*ZTT1(:)))*((ZLV(:)/(ZTT1(:)*XRV))-1.0))))
!
ZC(:)=4.0*XPI*(XRHOLW/PRHOD(:))*ZG(:)
ZDZRC(:)=0.0
ZDZRC(:)=ZC(:)*PS0(:)*ZRMOY(:)
MEM_PS0(:)=PS0(:)
!CALCUL DE B => Estimation de (drs/dT)ce
!T(=PTT) avant condensation; T'(=à determiner) apres condensation
!Lv(T),Cph(T)
ZLV(:) = XLVTT+(XCPV-XCL)*(PTT(:)-XTT)                
ZCPH(:)= XCPD+XCPV*PRV(:)+XCL*(PRC(:)+PRR(:))
!T'=T+(DT)ce
ZTT1(:)=PTT(:)+(ZDZRC(:)*PTSTEP*ZLV(:)/ZCPH(:))
!es(T')
ZW1(:)=EXP(XALPW-XBETAW/PTT(:)-XGAMW*LOG(PTT(:)))
!rs(T')
ZW1(:)=(XMV / XMD)*ZW1(:)/(PPRES(:)-ZW1(:))
!es(Tcond)
ZW2(:)=EXP(XALPW-XBETAW/ZTT1(:)-XGAMW*LOG(ZTT1(:)))
!rs(Tcond)
ZW2(:)=(XMV / XMD)*ZW2(:)/(PPRES(:)-ZW2(:))
!
WHERE (ZTT1(:).NE.PTT(:))
 ZB(:)=(ZLV(:)/ZCPH(:))*((ZW2(:)-ZW1(:))/(ZTT1(:)-PTT(:)))
ELSEWHERE
 ZB(:)=0.0
 ZDZRC(:)=0.0
ENDWHERE
!Calcul de S+dS
PS0(:)=PS0(:)+((ZA1(:)-(((ZB(:)*(PS0(:)+1.0)+1.0)*ZDZRC(:))/ZRVSAT1(:)))*PTSTEP)
!
PS0=MAX(PS0,-0.98)
!Ajustement tel que rv=(s+1)*rvs
ZTL(:)=PTT(:)-(PLV(:)/PCPH(:))*PRC(:)
ZRT(:)=PRC(:)+PRV(:)
ZDZRC2(:)=PRC(:)
DO J2=1,SIZE(ZDZRC,1)
  IF ((ZDZRC(J2).NE.0.0).OR.(ZDZRC2(J2).NE.0.0)) THEN 
    DO J1=1,5
     ZLV(J2) = XLVTT+(XCPV-XCL)*(PTT(J2)-XTT)                
     ZCPH(J2)=XCPD+XCPV*PRV(J2)+XCL*(PRC(J2)+PRR(J2))
     ZW1(J2)=EXP(XALPW-XBETAW/PTT(J2)-XGAMW*LOG(PTT(J2)))
     ZRVSAT1(J2)=(XMV / XMD)*ZW1(J2)/(PPRES(J2)-ZW1(J2))
     PRV(J2)=MIN(ZRT(J2),(PS0(J2)+1.0)*ZRVSAT1(J2))
     PRC(J2)=MAX(ZRT(J2)-PRV(J2),0.0)
     PTT(J2)=0.5*PTT(J2)+0.5*(ZTL(J2)+(ZLV(J2)*PRC(J2)/ZCPH(J2)))
    ENDDO
    ZLV(J2) = XLVTT+(XCPV-XCL)*(PTT(J2)-XTT)                
    ZCPH(J2)=XCPD+XCPV*PRV(J2)+XCL*(PRC(J2)+PRR(J2))
    PTT(J2)=ZTL(J2)+(ZLV(J2)*PRC(J2)/ZCPH(J2))
 ENDIF
ENDDO
ADJU2(:)=0.0
!
!Correction dans les mailles où ds a été surestimée
ZDZRC2(:)=PRC(:)-ZDZRC2(:)
WHERE ((MEM_PS0(:).LE.0.0).AND.(PS0(:).GT.0.0).AND.(ZDZRC2(:).LT.0.0))
  PS0(:)=0.0
  ADJU2(:)=1.0
ENDWHERE
!
WHERE ((MEM_PS0(:).GE.0.0).AND.(PS0(:).LT.0.0).AND.(ZDZRC2(:).GT.0.0))
  PS0(:)=0.0
  ADJU2(:)=1.0
ENDWHERE
!
DO J2=1,SIZE(ADJU2,1)
  IF (ADJU2(J2)==1) THEN 
   DO J1=1,5
    ZLV(J2) = XLVTT+(XCPV-XCL)*(PTT(J2)-XTT)                
    ZCPH(J2)=XCPD+XCPV*PRV(J2)+XCL*(PRC(J2)+PRR(J2))
    ZW1(J2)=EXP(XALPW-XBETAW/PTT(J2)-XGAMW*LOG(PTT(J2)))
    ZRVSAT1(J2)=(XMV / XMD)*ZW1(J2)/(PPRES(J2)-ZW1(J2))
    PRV(J2)=MIN(ZRT(J2),(PS0(J2)+1.0)*ZRVSAT1(J2))
    PRC(J2)=MAX(ZRT(J2)-PRV(J2),0.0)
    PTT(J2)=0.5*PTT(J2)+0.5*(ZTL(J2)+(ZLV(J2)*PRC(J2)/ZCPH(J2)))
   ENDDO
   ZLV(J2) = XLVTT+(XCPV-XCL)*(PTT(J2)-XTT)                
   ZCPH(J2)=XCPD+XCPV*PRV(J2)+XCL*(PRC(J2)+PRR(J2))
   PTT(J2)=ZTL(J2)+(ZLV(J2)*PRC(J2)/ZCPH(J2))
  ENDIF
ENDDO
!
!Elimination de l'eau liquide dans les mailles où le rayon des gouttelettes est
!inférieur à AER_RAD
ZRMOY(:)=0.0
DO J2=1,SIZE(PRC,1)
 IF ((PRC(J2).NE.0.) .AND. (PCCS(J2).NE.0.)) THEN
  ZRMOY(J2)=(MOMG(XALPHAC,XNUC,3.0)*4.0*XPI*PCCS(J2)*XRHOLW/&
        (3.0*PRC(J2)*PRHOD(J2)))**(1.0/3.0)
  ZRMOY(J2)=MOMG(XALPHAC,XNUC,1.0)/ZRMOY(J2)
  IF ((ZFLAG_ACT(J2).EQ.1.0).AND.(MEM_PS0(J2).LT.0.0).AND.(ZRMOY(J2).LT.AER_RAD)) THEN
   PTT(J2)=ZTL(J2)
   PRV(J2)=ZRT(J2)
   PRC(J2)=0.0
  ENDIF
 ENDIF
ENDDO
!
!Calcul de S au regard de T et rv en fin de pas de temps
ZW1=EXP(XALPW-XBETAW/PTT(:)-XGAMW*LOG(PTT(:)))
 !rvsat
ZRVSAT1(:)=(XMV / XMD)*ZW1(:)/(PPRES-ZW1(:))
!
WHERE (PRC(:)==0.0D0)
 PS0(:)=(PRV(:)/ZRVSAT1(:))-1D0
ENDWHERE
!
 DEALLOCATE(ZZW1,ZZW2,ZZW6,ZCHEN_MULTI,ZTMP,ZVEC1,IVEC1)
!
!
CONTAINS
!
FUNCTION MOMG (PALPHA,PNU,PP) RESULT (PMOMG)
USE MODI_GAMMA
IMPLICIT NONE
REAL     :: PALPHA ! first shape parameter of the DIMENSIONnal distribution
REAL     :: PNU    ! second shape parameter of the DIMENSIONnal distribution
REAL     :: PP     ! order of the moment
REAL     :: PMOMG  ! result: moment of order ZP
PMOMG = GAMMA(PNU+PP/PALPHA)/GAMMA(PNU)
!
END FUNCTION MOMG
!
END SUBROUTINE PROGNOS_LIMA
