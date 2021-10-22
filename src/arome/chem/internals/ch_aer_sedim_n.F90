!     ######spl
     SUBROUTINE CH_AER_SEDIM_n(PDTMONITOR,PSVT,PTHT,&
                       PRHODREF,PPABST, PVDEPAERO, PZZ, PSEDA)
     USE PARKIND1, ONLY : JPRB
     USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!!   #######################################
!!
!!   PURPOSE
!!   -------
!!
!!   REFERENCE
!!   ---------
!!   none
!!
!!   AUTHOR
!!    ------
!!    Pierre TULET (GMEI) / Vincent Crassier (LA)
!!
!!   MODIFICATIONS
!!    -------------
!!   Original
!!
! Entry variables:
!
! PM(IN)       -Array of moments
!
!*************************************************************
! Exit variables:
!
! PDMDEPOS(IN)  -Array of moment variation due to dry deposition
!
!*************************************************************
! Variables used during the deposition velocity calculation
! 
! ZVGK       -Polydisperse settling velocity of the kth moment (m/s)
!************************************************************
!!
!!   IMPLICIT ARGUMENTS
!
USE MODD_CH_AEROSOL
USE MODD_CH_AERO_n
USE MODI_CH_AER_VELGRAV_n
USE MODE_AERO_PSD
!
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!
REAL,    INTENT(IN) :: PDTMONITOR
REAL,  DIMENSION(:,:,:,:),  INTENT(IN)    :: PSVT
REAL,  DIMENSION(:,:,:),    INTENT(IN)    :: PTHT,PRHODREF, PZZ
REAL,  DIMENSION(:,:,:),    INTENT(IN)    :: PPABST
REAL,  DIMENSION(:,:,:),    INTENT(INOUT) :: PVDEPAERO
REAL,  DIMENSION(:,:,:,:),  INTENT(INOUT) :: PSEDA
!
!*       0.2   Declarations of local variables
!
INTEGER :: JT, JN, JJ, JK
INTEGER :: ISPLITA
REAL    :: ZTSPLITR
REAL,    DIMENSION(JPIN) :: ZPMIN
REAL,  DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3),JPMODE)   :: ZRG, ZSIG, ZRHOP
REAL,  DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3),JPMODE*3) :: ZPM, ZFLUXSED, ZFLUXMAX, ZPMOLD
REAL,  DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3)) :: ZW,ZH,ZMU
REAL,  DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3),JPMODE*3) :: ZVGK, ZDPK
REAL,  DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3),JPMODE)   :: ZVG
REAL,  DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3)) :: ZSUM
REAL,  DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3)) :: ZVSNUMMAX
REAL,  DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3),NSP+NCARB+NSOA,JPMODE) :: ZCTOTA, ZCCTOT
REAL :: ZVSMAX, ZHMIN
REAL    :: ZINIRADIUSI, ZINIRADIUSJ, ZRGMIN
INTEGER             :: ILU  ! indice K End       in z direction
!-------------------------------------------------------------------------------
!
!*       1.1   compute dimensions of arrays
!
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('CH_AER_SEDIM_N',0,ZHOOK_HANDLE)
ILU = SIZE(PSVT,3) -1

ZMU(:,:,:)   = 0.
ZH(:,:,:)    = 0.
ZW(:,:,:)    = 0.
ZVGK(:,:,:,:) = 0.
ZVG(:,:,:,:)  = 0.
ZDPK(:,:,:,:) = 0.
ZFLUXSED(:,:,:,:) = 0.
!
IF (CRGUNIT=="MASS") THEN
    ZINIRADIUSI = XINIRADIUSI * EXP(-3.*(LOG(XINISIGI))**2)
    ZINIRADIUSJ = XINIRADIUSJ * EXP(-3.*(LOG(XINISIGJ))**2)
ELSE
    ZINIRADIUSI = XINIRADIUSI 
    ZINIRADIUSJ = XINIRADIUSJ
END IF
!
!Get minimum values possible
ZPMIN(1) = XN0IMIN
ZRGMIN = XCOEFRADIMIN * ZINIRADIUSI
ZPMIN(2) = ZPMIN(1) * (ZRGMIN**3)*EXP(4.5 * LOG(XSIGIMIN)**2) 
ZPMIN(3) = ZPMIN(1) * (ZRGMIN**6)*EXP(18. * LOG(XSIGIMIN)**2)
ZPMIN(4) = XN0JMIN
ZRGMIN = XCOEFRADJMIN * ZINIRADIUSJ
ZPMIN(5) = ZPMIN(4) * (ZRGMIN**3)*EXP(4.5 * LOG(XSIGJMIN)**2) 
ZPMIN(6) = ZPMIN(4) * (ZRGMIN**6)*EXP(18. * LOG(XSIGJMIN)**2)
!
!
CALL  PPP2AERO(PSVT, PRHODREF, &
               PSIG3D=ZSIG, PRG3D=ZRG, PCTOTA=ZCTOTA, PM3D=ZPM)
!
ZPMOLD(:,:,:,:) = ZPM(:,:,:,:)
!
ZRHOP(:,:,:,:)   = 0.
ZCCTOT(:,:,:,:,:)= 0.
DO JN=1,JPMODE
  ZSUM(:,:,:)=0.
  DO JJ=1,NSP+NCARB+NSOA
    ZSUM(:,:,:)=ZSUM(:,:,:)+ZCTOTA(:,:,:,JJ,JN)/XRHOI(JJ)
  END DO
  DO JJ=1,NSP+NCARB+NSOA
    ZCCTOT(:,:,:,JJ,JN) = ZCTOTA(:,:,:,JJ,JN)/XRHOI(JJ)/ZSUM(:,:,:)
    ZRHOP(:,:,:,JN)=ZRHOP(:,:,:,JN)+ZCCTOT(:,:,:,JJ,JN)*XRHOI(JJ)
  ENDDO
ENDDO
!
CALL CH_AER_VELGRAV_n(ZSIG, ZRG, PTHT,  PPABST, PRHODREF, &
                      ZRHOP, ZMU, ZVGK, ZDPK,ZVG)
!
! Compute time-splitting condition
ZH=9999.
ZVSNUMMAX(:,:,:)   = 0.
DO JK=1,ILU
  ZH(:,:,JK)=PZZ(:,:,JK+1)-PZZ(:,:,JK)
  ! Maximum velocity is 10 m/s
  ZVSNUMMAX(:,:,JK) = MIN(10. *ZH(:,:,JK) / PDTMONITOR,20.)
END DO
ZHMIN=MINVAL(ZH(:,:,1:ILU))
!
ZFLUXSED(:,:,:,:) = 0.
ZFLUXMAX(:,:,:,:) = 0.
!
DO JN=1,JPIN
  !
  ZVGK(:,:,1:ILU,JN)  = MIN(ZVGK(:,:,1:ILU,JN), ZVSNUMMAX(:,:,1:ILU))
  !
  !PVDEPAERO(:,:,JN) = MIN(PVDEPAERO(:,:,JN), ZVSNUMMAX(:,:,1))
  PVDEPAERO(:,:,JN) = MIN(ZVGK(:,:,1,JN), ZVSNUMMAX(:,:,1))
  !
  ZVSMAX=MAXVAL(ZVGK(:,:,1:ILU,JN))
  !
  ISPLITA = INT(ZVSMAX*PDTMONITOR/ZHMIN)+1
  ISPLITA = MIN(50, ISPLITA)
  !
  ZTSPLITR  = PDTMONITOR / FLOAT(ISPLITA)   
  !
  DO JT=1,ISPLITA
    ZFLUXSED(:,:,1:ILU+1,JN)= ZVGK(:,:,1:ILU+1,JN)* ZPM(:,:,1:ILU+1,JN)
! first level: aerosol dry deposition velocity (turbulence + gravitation)
    ZFLUXSED(:,:,1,JN) =   PVDEPAERO(:,:,JN)* ZPM(:,:,1,JN)
    DO JK=1,ILU
      ZW(:,:,JK)  = ZTSPLITR /(PZZ(:,:,JK+1)-PZZ(:,:,JK))
    END DO
    !
    ZFLUXMAX(:,:,1:ILU,JN)  =  ZPM(:,:,1:ILU,JN) / ZW(:,:,1:ILU)
    ZFLUXSED(:,:,1:ILU,JN)  = MAX(MIN(ZFLUXSED(:,:,1:ILU,JN), ZFLUXMAX(:,:,1:ILU,JN)),0.)
    !
    DO JK=1,ILU
      ZPM(:,:,JK,JN)= ZPM(:,:,JK,JN) + &
                      ZW(:,:,JK)*(ZFLUXSED(:,:,JK+1,JN)- ZFLUXSED(:,:,JK,JN))

    END DO
  END DO
  !
  ! No variation more than 10 % each time step.
  !
  ZPM(:,:,:,JN) = MIN(MAX(ZPM(:,:,:,JN),0.90*ZPMOLD(:,:,:,JN)),1.10*ZPMOLD(:,:,:,JN))
!
END DO
!
DO JN=1,JPMODE
! Calcul pour maintenir le rayon fixe pour la sedimentation :
 IF (LRGFIX) THEN
    ZPM(:,:,:,NM3(JN)) = ZPM(:,:,:,NM0(JN)) *&
                        (ZRG(:,:,:,JN)**3)*EXP(4.5 * LOG(ZSIG(:,:,:,JN))**2)
 ENDIF

! Calcul pour maintenir la dispersion fixe pour la sedimentation :
! sinon Rg augmente lors de la reconstruction d'une loi log-normale
! a partir des nouveau Mk (sigma diminue plus vite que Rg)

! calcul de M6 en conservant sigma
ZPM(:,:,1:ILU, NM6(JN)) = ZPM(:,:,1:ILU,NM0(JN)) &
       * ( (ZPM(:,:,1:ILU,NM3(JN))/ZPM(:,:,1:ILU,NM0(JN)))**(1./3.) &
       * exp(-(3./2.)*LOG(ZSIG(:,:,1:ILU,JN))**2))**6 &
       * exp(18.*LOG(ZSIG(:,:,1:ILU,JN))**2)

  IF ((GSEDFIX).AND.&
      (((JN .EQ. 1).AND. (LVARSIGI)).OR.&
       ((JN .EQ. 2).AND. (LVARSIGJ)))) THEN
! calcul de M6 en conservant Rg
    ZPM(:,:,2:ILU,NM6(JN)) = ZPM(:,:,2:ILU,NM3(JN)) ** 4 / &
                     (ZRG(:,:,2:ILU,JN)**6 * ZPM(:,:,2:ILU,NM0(JN))**3)

  END IF
END DO
!
IF (GSEDFIX) THEN
  GSEDFIX = .FALSE.
ELSE
  GSEDFIX = .TRUE.
END IF
!
PSEDA(:,:,:,:) = 0.
DO JN=1,JPMODE
  WHERE ((ZPM(:,:,:,NM0(JN)) .LE. ZPMIN(NM0(JN))).OR.&
         (ZPM(:,:,:,NM3(JN)) .LE. ZPMIN(NM3(JN))).OR.& 
         (ZPM(:,:,:,NM6(JN)) .LE. ZPMIN(NM6(JN)))) 
    ZPM(:,:,:,NM0(JN)) = ZPMOLD(:,:,:,NM0(JN)) 
    ZPM(:,:,:,NM3(JN)) = ZPMOLD(:,:,:,NM3(JN)) 
    ZPM(:,:,:,NM6(JN)) = ZPMOLD(:,:,:,NM6(JN)) 
  END WHERE
ENDDO 
!
DO JN=1,JPIN
 PSEDA(:,:,1:ILU,JN) = (ZPM(:,:,1:ILU,JN) -  &
                                ZPMOLD(:,:,1:ILU,JN)) / PDTMONITOR
END DO
!
IF (LHOOK) CALL DR_HOOK('CH_AER_SEDIM_N',1,ZHOOK_HANDLE)
END SUBROUTINE CH_AER_SEDIM_n
