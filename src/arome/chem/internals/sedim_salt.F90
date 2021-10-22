!     ######spl
     SUBROUTINE SEDIM_SALT(PTHT,PDTMONITOR,&
                       PRHODREF,PPABST,PZZ,PSVT)
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
!!    Pierre TULET (GMEI) 
!!
!!   MODIFICATIONS
!!    -------------
!!   Original
!!
! Entry variables:
!
! PSVTS(INOUT)       -Array of moments included in PSVTS
!
!*************************************************************
! Exit variables:
!
!*************************************************************
! Variables used during the deposition velocity calculation
! 
! ZVGK       -Polydisperse settling velocity of the kth moment (m/s)
!************************************************************
!!
!!   IMPLICIT ARGUMENTS
!
USE MODD_SALT
USE MODD_CSTS_SALT
USE MODI_SALT_VELGRAV
USE MODE_SALT_PSD
!
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL,  INTENT(IN) :: PDTMONITOR
REAL,  DIMENSION(:,:,:,:),  INTENT(INOUT)    :: PSVT  !scalar variable 
REAL,  DIMENSION(:,:,:),    INTENT(IN)    :: PTHT,PRHODREF, PZZ
REAL,  DIMENSION(:,:,:),    INTENT(IN)    :: PPABST
!
!*      0.2    declarations of local variables
!
INTEGER :: JK, JN, JT
INTEGER :: IMODEIDX
INTEGER :: ISPLITA
REAL    :: ZTSPLITR
REAL,    DIMENSION(NMODE_SLT*3) :: ZPMIN
REAL,    DIMENSION(NMODE_SLT)   :: ZINIRADIUS
REAL :: ZRGMIN, ZSIGMIN
REAL,  DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3),NMODE_SLT)   :: ZRG, ZSIG
REAL,  DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3),3*NMODE_SLT) :: ZPM, ZPMOLD
REAL,  DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3)+1,3*NMODE_SLT) :: ZFLUXSED, ZFLUXMAX
REAL,  DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3)) :: ZH,ZMU, ZW, ZVSNUMMAX
REAL,  DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3),3*NMODE_SLT) :: ZVGK, ZDPK
REAL,  DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3),NMODE_SLT)   :: ZVG
REAL :: ZVSMAX, ZHMIN, ZRHOI, ZFAC
INTEGER             :: ILU  ! indice K End       in z direction
INTEGER,DIMENSION(NMODE_SLT) :: NM0, NM3, NM6
REAL                :: ZDEN2MOL, ZAVOGADRO, ZMD, ZMI 
LOGICAL, SAVE       :: LSEDFIX = .TRUE.
!
!*       0.3   initialize constant
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('SEDIM_SALT',0,ZHOOK_HANDLE)
ZRHOI = XDENSITY_SALT

ZMU(:,:,:)     = 0.
ZH(:,:,:)      = 0.
ZVGK(:,:,:,:)  = 0.
ZVG(:,:,:,:)   = 0.
ZDPK(:,:,:,:)  = 0.
ZW(:,:,:)      = 0.
ZFLUXSED(:,:,:,:) = 0.
ILU = SIZE(PSVT,3)
ZFAC  = (4./3.)*XPI*ZRHOI*1.e-9
ZAVOGADRO = 6.0221367E+23
ZMD = 28.9644E-3
ZDEN2MOL = 1E-6 * ZAVOGADRO  / ZMD
ZMI   = XMOLARWEIGHT_SALT ! molecular mass in kg/mol
!
!*       1.   compute dimensions of arrays
!
DO JN=1,NMODE_SLT
  IMODEIDX = JPSALTORDER(JN)
  !Calculations here are for one mode only
  IF (CRGUNITS=="MASS") THEN
    ZINIRADIUS(JN) = XINIRADIUS_SLT(IMODEIDX) * EXP(-3.*(LOG(XINISIG_SLT(IMODEIDX)))**2)
  ELSE
    ZINIRADIUS(JN) = XINIRADIUS_SLT(IMODEIDX)
  END IF

  !Set counter for number, M3 and M6
  NM0(JN) = 1+(JN-1)*3
  NM3(JN) = 2+(JN-1)*3
  NM6(JN) = 3+(JN-1)*3
  !Get minimum values possible
  ZPMIN(NM0(JN)) = XN0MIN_SLT(IMODEIDX)
  ZRGMIN     = ZINIRADIUS(JN)
  IF (LVARSIG_SLT) THEN
    ZSIGMIN = XSIGMIN_SLT
  ELSE
    ZSIGMIN = XINISIG_SLT(IMODEIDX)
  ENDIF
  ZPMIN(NM3(JN)) = XN0MIN_SLT(IMODEIDX) * (ZRGMIN**3)*EXP(4.5 * LOG(ZSIGMIN)**2) 
  ZPMIN(NM6(JN)) = XN0MIN_SLT(IMODEIDX) * (ZRGMIN**6)*EXP(18. * LOG(ZSIGMIN)**2)
END DO
!
!*       2.   compute SIG, RG and moments
!
CALL PPP2SALT(PSVT, PRHODREF, PSIG3D=ZSIG, PRG3D=ZRG, PM3D=ZPM)

ZPMOLD(:,:,:,:)=ZPM(:,:,:,:)
!
!*       3.   compute gravitational velocities
!
CALL SALT_VELGRAV(ZSIG(:,:,1:ILU,:), ZRG(:,:,1:ILU,:),    &
                      PTHT(:,:,1:ILU),  PPABST(:,:,1:ILU),&
                      PRHODREF(:,:,1:ILU), ZRHOI,         &
                      ZMU(:,:,1:ILU), ZVGK(:,:,1:ILU,:),  &
                      ZDPK(:,:,1:ILU,:),ZVG(:,:,1:ILU,:))
!
!*       4.   Compute time-splitting condition
!
ZH=9999.
ZVSNUMMAX(:,:,:)   = 0.

DO JK=1,ILU
   ZH(:,:,JK)=PZZ(:,:,JK+1)-PZZ(:,:,JK)
   ! Maximum velocity 
   ZVSNUMMAX(:,:,JK) = MIN(10. *ZH(:,:,JK) / PDTMONITOR,20.)
ENDDO
!
ZHMIN=MINVAL(ZH(:,:,1:ILU))

ZFLUXSED(:,:,:,:) = 0.
ZFLUXMAX(:,:,:,:) = 0.

!Get loss rate [1/s] instead of [m/s]
DO JN=1,NMODE_SLT*3
  !ZVGK(:,:,:,JN)=ZVGK(:,:,:,JN)/ZH(:,:,:)
  ZVGK(:,:,1:ILU,JN)  = MIN(ZVGK(:,:,1:ILU,JN), ZVSNUMMAX(:,:,1:ILU))
  ZVSMAX=MAXVAL(ZVGK(:,:,1:ILU,JN))
  ZVSMAX=MAX(ZVSMAX,5.)
  !
  ISPLITA = INT(ZVSMAX*PDTMONITOR/ZHMIN)+1
  ISPLITA = MIN(20, ISPLITA)
  !
  ZTSPLITR  = PDTMONITOR / FLOAT(ISPLITA)   
  !
  ZFLUXSED(:,:,ILU+1,JN) = 0.
  
  DO JT=1,ISPLITA
    ZFLUXSED(:,:,1:ILU,JN)= ZVGK(:,:,1:ILU,JN)* ZPM(:,:,1:ILU,JN)
    DO JK=1,ILU
      ZW(:,:,JK)  = ZTSPLITR /(PZZ(:,:,JK+1)-PZZ(:,:,JK))
    END DO
    !
    !
!
    DO JK=1,ILU
      ZPM(:,:,JK,JN)= ZPM(:,:,JK,JN) + &
                      ZW(:,:,JK)*(ZFLUXSED(:,:,JK+1,JN)- ZFLUXSED(:,:,JK,JN))

    END DO
  ENDDO

ENDDO
!
DO JN=1,NMODE_SLT
!! Calcul pour maintenir le rayon fixe pour la sedimentation :
  IF (LRGFIX_SLT) THEN
    ZPM(:,:,:,NM3(JN)) = ZPM(:,:,:,NM0(JN)) *&
                        (ZRG(:,:,:,JN)**3)*EXP(4.5 * LOG(ZSIG(:,:,:,JN))**2)
  ENDIF
! Calcul pour maintenir la dispersion fixe pour la sedimentation :
! sinon Rg augmente lors de la reconstruction d'une loi log-normale
! a partir des nouveau Mk (sigma diminue plus vite que Rg)
! calcul de M6 en conservant sigma
  ZPM(:,:,:, NM6(JN)) = ZPM(:,:,:,NM0(JN)) &
       * ( (ZPM(:,:,:,NM3(JN))/ZPM(:,:,:,NM0(JN)))**(1./3.) &
       * exp(-(3./2.)*LOG(ZSIG(:,:,:,JN))**2))**6 &
       * exp(18.*LOG(ZSIG(:,:,:,JN))**2)

  IF ((LSEDFIX).AND.(LVARSIG_SLT)) THEN
  ! calcul de M6 en conservant Rg

  ZPM(:,:,:,NM6(JN)) = ZPM(:,:,:,NM3(JN)) ** 4 / &
                     (ZRG(:,:,:,JN)**6 * ZPM(:,:,:,NM0(JN))**3)

  END IF
ENDDO
!
IF (LSEDFIX) THEN
  LSEDFIX = .FALSE.
ELSE
  LSEDFIX = .TRUE.
END IF
!
!*       5.   Return to concentration in ppp (#/molec_{air})
!
DO JN=1,NMODE_SLT
IF (LVARSIG_SLT) THEN
WHERE ((ZPM(:,:,:,NM0(JN)) .GT. 100.*ZPMIN(NM0(JN))).AND.&
       (ZPM(:,:,:,NM3(JN)) .GT. 100.*ZPMIN(NM3(JN))).AND.&
       (ZPM(:,:,:,NM6(JN)) .GT. 100.*ZPMIN(NM6(JN))))
  PSVT(:,:,:,1+(JN-1)*3) = ZPM(:,:,:,NM0(JN)) * XMD / &
                           (XAVOGADRO*PRHODREF(:,:,:))
  PSVT(:,:,:,2+(JN-1)*3) = ZPM(:,:,:,NM3(JN)) * XMD*XPI * 4./3.*ZRHOI  / &
                           (ZMI*PRHODREF(:,:,:)*XM3TOUM3_SALT)
  PSVT(:,:,:,3+(JN-1)*3) = ZPM(:,:,:,NM6(JN)) * XMD / &
                           (XAVOGADRO*PRHODREF(:,:,:)*1.d-6)
ELSEWHERE
  PSVT(:,:,:,1+(JN-1)*3) = ZPMOLD(:,:,:,NM0(JN)) * XMD / &
                           (XAVOGADRO*PRHODREF(:,:,:))
  PSVT(:,:,:,2+(JN-1)*3) = ZPMOLD(:,:,:,NM3(JN)) * XMD*XPI * 4./3.*ZRHOI  / &
                           (ZMI*PRHODREF(:,:,:)*XM3TOUM3_SALT)
  PSVT(:,:,:,3+(JN-1)*3) = ZPMOLD(:,:,:,NM6(JN)) * XMD / &
                           (XAVOGADRO*PRHODREF(:,:,:)*1.d-6)
ENDWHERE
ELSE IF (LRGFIX_SLT) THEN
WHERE ((ZPM(:,:,:,NM3(JN)) .GT. 100*ZPMIN(NM3(JN))))
  PSVT(:,:,:,JN)     = ZPM(:,:,:,NM3(JN)) * XMD*XPI * 4./3.*ZRHOI  / &
                           (ZMI*PRHODREF(:,:,:)*XM3TOUM3_SALT)
ELSEWHERE
  PSVT(:,:,:,JN)     = ZPMOLD(:,:,:,NM3(JN)) * XMD*XPI * 4./3.*ZRHOI  / &
                           (ZMI*PRHODREF(:,:,:)*XM3TOUM3_SALT)
ENDWHERE

ELSE


WHERE ((ZPM(:,:,:,NM0(JN)) .GT. 100*ZPMIN(NM0(JN))).AND.&
       (ZPM(:,:,:,NM3(JN)) .GT. 100*ZPMIN(NM3(JN))))
  PSVT(:,:,:,1+(JN-1)*2) = ZPM(:,:,:,NM0(JN)) * XMD / &
                           (XAVOGADRO*PRHODREF(:,:,:))
  PSVT(:,:,:,2+(JN-1)*2) = ZPM(:,:,:,NM3(JN)) * XMD*XPI * 4./3.*ZRHOI  / &
                           (ZMI*PRHODREF(:,:,:)*XM3TOUM3_SALT)
ELSEWHERE
  PSVT(:,:,:,1+(JN-1)*2) = ZPMOLD(:,:,:,NM0(JN)) * XMD / &
                           (XAVOGADRO*PRHODREF(:,:,:))
  PSVT(:,:,:,2+(JN-1)*2) = ZPMOLD(:,:,:,NM3(JN)) * XMD*XPI * 4./3.*ZRHOI  / &
                           (ZMI*PRHODREF(:,:,:)*XM3TOUM3_SALT)
ENDWHERE

END IF
ENDDO

CALL PPP2SALT(PSVT, PRHODREF, PSIG3D=ZSIG, PRG3D=ZRG, PM3D=ZPM)

!
IF (LHOOK) CALL DR_HOOK('SEDIM_SALT',1,ZHOOK_HANDLE)
END SUBROUTINE SEDIM_SALT
