!     ######spl
     SUBROUTINE CH_AER_INIT(PCHEM,PAERO, PRHODREF)
     USE PARKIND1, ONLY : JPRB
     USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!!   ############################################################
!!
!!    PURPOSE
!!    -------
!!    Realise l'equilibre entre les moments via la masse contenue 
!!    dans les aerosols, les diametres moyens et la dispersion.
!!
!!    REFERENCE
!!    ---------
!!    none
!!
!!    AUTHOR
!!    ------
!!    Pierre TULET (LA)
!!
!!    MODIFICATIONS
!!    -------------
!!    none
!!
!!    EXTERNAL
!!    --------
!!    None
!!
USE MODD_CH_AEROSOL
USE MODD_CH_AERO_n
USE MODD_CH_M9,     ONLY : CNAMES
USE MODD_CH_MNHC_n, ONLY : LCH_INIT_FIELD
USE MODD_NSV
!!
IMPLICIT NONE
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
!*      0.1    declarations of arguments
!
REAL,   DIMENSION(:,:,:,:),    INTENT(INOUT)   :: PCHEM, PAERO
REAL,   DIMENSION(:,:,:),      INTENT(IN)      :: PRHODREF
!
!
!*      0.2    declarations local variables
!
REAL,DIMENSION(SIZE(PCHEM,1),SIZE(PCHEM,2),SIZE(PCHEM,3),JPMODE) :: ZN, ZRG, ZSIG
REAL,DIMENSION(SIZE(PCHEM,1),SIZE(PCHEM,2),SIZE(PCHEM,3),NSP+NCARB+NSOA,JPMODE) :: ZCTOTA
REAL,DIMENSION(SIZE(PCHEM,1),SIZE(PCHEM,2),SIZE(PCHEM,3)) :: ZSIGMA
REAL,DIMENSION(SIZE(PCHEM,1),SIZE(PCHEM,2),SIZE(PCHEM,3)) :: ZBCMINI, ZBCMINJ, ZOCMINI, ZOCMINJ
REAL,DIMENSION(SIZE(PCHEM,1),SIZE(PCHEM,2),SIZE(PCHEM,3),JPMODE*3) :: ZM, ZPM
REAL,DIMENSION(NSP+NCARB+NSOA) :: ZMI
INTEGER :: JN, JJ,  JK  ! loop counter
REAL    :: ZDEN2MOL, ZRHODREFMIN, ZCOEFAEROBC, ZCOEFAEROOC
REAL    :: ZVALBC, ZVALOC, ZMINRGI, ZMINRGJ
REAL    :: ZINIRADIUSI, ZINIRADIUSJ
!
!-------------------------------------------------------------------------------
!

!*       1.     TRANSFER FROM GAS TO AEROSOL MODULE
!               ------------------------------------
!        1.1    initialisation 
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('CH_AER_INIT',0,ZHOOK_HANDLE)
DO JJ=1,SIZE(CNAMES)
IF (CNAMES(JJ) == "CO")  JP_CH_CO   = JJ
END DO

IF (CRGUNIT=="MASS") THEN
ZINIRADIUSI = XINIRADIUSI * EXP(-3.*(LOG(XINISIGI))**2)
ZINIRADIUSJ = XINIRADIUSJ * EXP(-3.*(LOG(XINISIGJ))**2)
ELSE
ZINIRADIUSI = XINIRADIUSI 
ZINIRADIUSJ = XINIRADIUSJ
END IF
ZMINRGI = ZINIRADIUSI * XCOEFRADIMIN
ZMINRGJ = ZINIRADIUSJ * XCOEFRADJMIN

DO JK=1, size(PRHODREF,3)
WHERE (PRHODREF(:,:,JK) .GT. 0.)
ZPM(:,:,JK,1) = 100.* XN0IMIN*PRHODREF(:,:,JK)/PRHODREF(:,:,size(PRHODREF,3))
ELSEWHERE
ZPM(:,:,JK,1) = XN0IMIN
ENDWHERE
END DO
ZPM(:,:,:,2) = ZPM(:,:,:,1) * (ZINIRADIUSI**3)*EXP(4.5 * LOG(XINISIGI)**2) 

DO JK=1, size(PRHODREF,3)
WHERE (PRHODREF(:,:,JK) .GT. 0.)
ZPM(:,:,JK,4) = 100.* XN0JMIN*PRHODREF(:,:,JK)/PRHODREF(:,:,size(PRHODREF,3))
ELSEWHERE
ZPM(:,:,JK,4) = XN0JMIN
ENDWHERE
END DO
ZPM(:,:,:,5) = ZPM(:,:,:,4) * (ZINIRADIUSJ**3)*EXP(4.5 * LOG(XINISIGJ)**2) 
ZDEN2MOL = 1E-6 * XAVOGADRO / XMD

! Default values of molar mass

ZMI(:) = 250.
ZMI(JP_AER_SO4)  = 98.
ZMI(JP_AER_NO3)  = 63.
ZMI(JP_AER_NH3)  = 17.
ZMI(JP_AER_H2O)  = 18.

IF (NSOA .EQ. 10) THEN
ZMI(JP_AER_SOA1) = 88. 
ZMI(JP_AER_SOA2) = 180.
ZMI(JP_AER_SOA3) = 1.5374857E+02
ZMI(JP_AER_SOA4) = 1.9586780E+02
ZMI(JP_AER_SOA5) = 195.
ZMI(JP_AER_SOA6) = 195.
ZMI(JP_AER_SOA7) = 165.
ZMI(JP_AER_SOA8) = 195.
ZMI(JP_AER_SOA9) = 270.
ZMI(JP_AER_SOA10) = 210.
END IF

! Aerosol Density
! Cf Ackermann (all to black carbon except water)
 XRHOI(:) = 1.8e3
 XRHOI(JP_AER_H2O) = 1.0e3   ! water

PCHEM(:,:,:,:) = MAX(PCHEM(:,:,:,:), 0.)
PAERO(:,:,:,:) = MAX(PAERO(:,:,:,:), 0.)
!
DO JJ=1,NSP+NCARB+NSOA
  XFAC(JJ)=(4./3.)*3.14292654*XRHOI(JJ)*1.e-9
ENDDO
!
!
!*       1.n    transfer aerosol mass from gas to aerosol variables
!               (and conversion of part/part --> microgram/m3)
!
DO JJ=1,NSV_AER
  PAERO(:,:,:,JJ) =  PAERO(:,:,:,JJ) * ZDEN2MOL * PRHODREF(:,:,:)
ENDDO

ZBCMINI(:,:,:) = 0.6 * 6.0221367E+11 * XFAC(JP_AER_BC) * ZPM(:,:,:,2) / ZMI(JP_AER_BC)
ZOCMINI(:,:,:) = 0.4 * 6.0221367E+11 * XFAC(JP_AER_OC) * ZPM(:,:,:,2) / ZMI(JP_AER_OC)
ZBCMINJ(:,:,:) = 0.6 * 6.0221367E+11 * XFAC(JP_AER_BC) * ZPM(:,:,:,5) / ZMI(JP_AER_BC)
ZOCMINJ(:,:,:) = 0.4 * 6.0221367E+11 * XFAC(JP_AER_OC) * ZPM(:,:,:,5) / ZMI(JP_AER_OC)

PAERO(:,:,:,JP_CH_BCi)=MAX(PAERO(:,:,:,JP_CH_BCi), ZBCMINI(:,:,:))
PAERO(:,:,:,JP_CH_BCj)=MAX(PAERO(:,:,:,JP_CH_BCj), ZBCMINJ(:,:,:))

PAERO(:,:,:,JP_CH_OCi)=MAX(PAERO(:,:,:,JP_CH_OCi), ZOCMINI(:,:,:))
PAERO(:,:,:,JP_CH_OCj)=MAX(PAERO(:,:,:,JP_CH_OCj), ZOCMINJ(:,:,:))

!
! mineral phase
  ZCTOTA(:,:,:,JP_AER_SO4,1) = PAERO(:,:,:,JP_CH_SO4i)*ZMI(JP_AER_SO4)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SO4,2) = PAERO(:,:,:,JP_CH_SO4j)*ZMI(JP_AER_SO4)/6.0221367E+11

  ZCTOTA(:,:,:,JP_AER_NO3,1) = PAERO(:,:,:,JP_CH_NO3i)*ZMI(JP_AER_NO3)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_NO3,2) = PAERO(:,:,:,JP_CH_NO3j)*ZMI(JP_AER_NO3)/6.0221367E+11

  ZCTOTA(:,:,:,JP_AER_NH3,1) = PAERO(:,:,:,JP_CH_NH3i)*ZMI(JP_AER_NH3)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_NH3,2) = PAERO(:,:,:,JP_CH_NH3j)*ZMI(JP_AER_NH3)/6.0221367E+11

! water
  ZCTOTA(:,:,:,JP_AER_H2O,1) = PAERO(:,:,:,JP_CH_H2Oi)*ZMI(JP_AER_H2O)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_H2O,2) = PAERO(:,:,:,JP_CH_H2Oj)*ZMI(JP_AER_H2O)/6.0221367E+11

!
! primary organic carbon
  ZCTOTA(:,:,:,JP_AER_OC,1) = PAERO(:,:,:,JP_CH_OCi)*ZMI(JP_AER_OC)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_OC,2) = PAERO(:,:,:,JP_CH_OCj)*ZMI(JP_AER_OC)/6.0221367E+11

! primary black carbon
  ZCTOTA(:,:,:,JP_AER_BC,1) = PAERO(:,:,:,JP_CH_BCi)*ZMI(JP_AER_BC)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_BC,2) = PAERO(:,:,:,JP_CH_BCj)*ZMI(JP_AER_BC)/6.0221367E+11

!
IF (NSOA .EQ. 10) THEN
  ZCTOTA(:,:,:,JP_AER_SOA1,1) = PAERO(:,:,:,JP_CH_SOA1i)*ZMI(JP_AER_SOA1)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA1,2) = PAERO(:,:,:,JP_CH_SOA1j)*ZMI(JP_AER_SOA1)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA2,1) = PAERO(:,:,:,JP_CH_SOA2i)*ZMI(JP_AER_SOA2)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA2,2) = PAERO(:,:,:,JP_CH_SOA2j)*ZMI(JP_AER_SOA2)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA3,1) = PAERO(:,:,:,JP_CH_SOA3i)*ZMI(JP_AER_SOA3)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA3,2) = PAERO(:,:,:,JP_CH_SOA3j)*ZMI(JP_AER_SOA3)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA4,1) = PAERO(:,:,:,JP_CH_SOA4i)*ZMI(JP_AER_SOA4)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA4,2) = PAERO(:,:,:,JP_CH_SOA4j)*ZMI(JP_AER_SOA4)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA5,1) = PAERO(:,:,:,JP_CH_SOA5i)*ZMI(JP_AER_SOA5)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA5,2) = PAERO(:,:,:,JP_CH_SOA5j)*ZMI(JP_AER_SOA5)/6.0221367E+11

  ZCTOTA(:,:,:,JP_AER_SOA6,1) = PAERO(:,:,:,JP_CH_SOA6i)*ZMI(JP_AER_SOA6)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA6,2) = PAERO(:,:,:,JP_CH_SOA6j)*ZMI(JP_AER_SOA6)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA7,1) = PAERO(:,:,:,JP_CH_SOA7i)*ZMI(JP_AER_SOA7)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA7,2) = PAERO(:,:,:,JP_CH_SOA7j)*ZMI(JP_AER_SOA7)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA8,1) = PAERO(:,:,:,JP_CH_SOA8i)*ZMI(JP_AER_SOA8)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA8,2) = PAERO(:,:,:,JP_CH_SOA8j)*ZMI(JP_AER_SOA8)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA9,1) = PAERO(:,:,:,JP_CH_SOA9i)*ZMI(JP_AER_SOA9)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA9,2) = PAERO(:,:,:,JP_CH_SOA9j)*ZMI(JP_AER_SOA9)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA10,1) = PAERO(:,:,:,JP_CH_SOA10i)*ZMI(JP_AER_SOA10)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA10,2) = PAERO(:,:,:,JP_CH_SOA10j)*ZMI(JP_AER_SOA10)/6.0221367E+11
END IF


!*       1.1    calculate moment 3 from mass
    
    ZM(:,:,:,2) = 0.
    ZM(:,:,:,5) = 0.
    ZCTOTA(:,:,:,:,:) = MAX(ZCTOTA(:,:,:,:,:), 0.)
    DO JJ = 1,NSP+NCARB+NSOA
    ZM(:,:,:,2) = ZM(:,:,:,2)+ZCTOTA(:,:,:,JJ,1)/XFAC(JJ)
    ZM(:,:,:,5) = ZM(:,:,:,5)+ZCTOTA(:,:,:,JJ,2)/XFAC(JJ)
    ENDDO
!
!
!*       1.2    calculate moment 0 from dispersion and mean radius
    ZM(:,:,:,1)= ZM(:,:,:,2) / &
               ((ZINIRADIUSI**3)*EXP(4.5 * (LOG(XINISIGI))**2))

    ZM(:,:,:,4)= ZM(:,:,:,5) / &
               ((ZINIRADIUSJ**3)*EXP(4.5 * (LOG(XINISIGJ))**2))

!*       1.3    calculate moment 6 from dispersion and mean radius
   ZM(:,:,:,3) = ZM(:,:,:,1) * (ZINIRADIUSI**6) *EXP(18 *(LOG(XINISIGI))**2)
   ZM(:,:,:,6) = ZM(:,:,:,4) * (ZINIRADIUSJ**6) *EXP(18 *(LOG(XINISIGJ))**2)

IF (LVARSIGI) THEN ! set M6 variable standard deviation
 ZM(:,:,:,3) = MAX(PAERO(:,:,:,JP_CH_M6i), 1E-80)
ELSE ! fixed standard deviation
 ZM(:,:,:,3) = ZM(:,:,:,1) &
          * ( (ZM(:,:,:,2)/ZM(:,:,:,1))**(1./3.)  &
          * exp(-(3./2.)*log(XINISIGI)**2))**6 &
          * exp(18.*log(XINISIGI)**2)
END IF

IF (LVARSIGJ) THEN ! set M6 variable standard deviation
 ZM(:,:,:,6) = MAX(PAERO(:,:,:,JP_CH_M6j), 1E-80)
ELSE ! fixed standard deviation
 ZM(:,:,:,6) = ZM(:,:,:,4) &
          * ( (ZM(:,:,:,5)/ZM(:,:,:,4))**(1./3.)  &
          * exp(-(3./2.)*log(XINISIGJ)**2))**6 &
          * exp(18.*log(XINISIGJ)**2)
END IF


DO JN=1,JPMODE
  IF (JN .EQ. 1) THEN

    IF (LVARSIGI) THEN ! variable dispersion for mode 1

      ZSIGMA(:,:,:)=ZM(:,:,:,NM3(JN))**2/(ZM(:,:,:,NM0(JN))*ZM(:,:,:,NM6(JN)))
      ZSIGMA(:,:,:)=MIN(1-1E-10,ZSIGMA(:,:,:))
      ZSIGMA(:,:,:)=MAX(1E-10,ZSIGMA(:,:,:))
      ZSIGMA(:,:,:)= LOG(ZSIGMA(:,:,:))
      ZSIGMA(:,:,:)= EXP(1./3.*SQRT(-ZSIGMA(:,:,:)))
      WHERE (ZSIGMA(:,:,:) > XSIGIMAX)
      ZSIGMA(:,:,:) =  XSIGIMAX
      END WHERE
      WHERE (ZSIGMA(:,:,:) < XSIGIMIN)
      ZSIGMA(:,:,:) =  XSIGIMIN
      END WHERE

    ELSE ! fixed dispersion for mode 1
      ZSIGMA(:,:,:) = XINISIGI
    END IF
  END IF

!
  IF (JN .EQ. 2) THEN

    IF (LVARSIGJ) THEN ! variable dispersion for mode 2

      ZSIGMA(:,:,:)=ZM(:,:,:,NM3(JN))**2/(ZM(:,:,:,NM0(JN))*ZM(:,:,:,NM6(JN)))
      ZSIGMA(:,:,:)=MIN(1-1E-10,ZSIGMA(:,:,:))
      ZSIGMA(:,:,:)=MAX(1E-10,ZSIGMA(:,:,:))
      ZSIGMA(:,:,:)= LOG(ZSIGMA(:,:,:))
      ZSIGMA(:,:,:)= EXP(1./3.*SQRT(-ZSIGMA(:,:,:)))
      WHERE (ZSIGMA(:,:,:) > XSIGJMAX)
      ZSIGMA(:,:,:) =  XSIGJMAX
      END WHERE
      WHERE (ZSIGMA(:,:,:) < XSIGJMIN)
      ZSIGMA(:,:,:) =  XSIGJMIN
      END WHERE

    ELSE ! fixed dispersion for mode 2
      ZSIGMA(:,:,:) = XINISIGJ
    END IF
  END IF
!


!*       1.4    calculate modal parameters from moments
ZSIG(:,:,:,JN) = ZSIGMA(:,:,:)
ZN(:,:,:,JN) = ZM(:,:,:,NM0(JN))

ZSIGMA(:,:,:)=LOG(ZSIG(:,:,:,JN))**2

ZRG(:,:,:,JN)=(ZM(:,:,:,NM3(JN))/ZN(:,:,:,JN))**(1./3.)*EXP(-1.5*ZSIGMA(:,:,:))

ZM(:,:,:,NM6(JN))=ZN(:,:,:,JN)*ZRG(:,:,:,JN)**6*EXP(18.*ZSIGMA(:,:,:))

!
ENDDO
!
!
PAERO(:,:,:,JP_CH_M0i) = ZM(:,:,:,1) * 1E-6 
PAERO(:,:,:,JP_CH_M0j) = ZM(:,:,:,4) * 1E-6

IF (LVARSIGI) PAERO(:,:,:,JP_CH_M6i) = ZM(:,:,:,3) 
IF (LVARSIGJ) PAERO(:,:,:,JP_CH_M6j) = ZM(:,:,:,6)

!
  DO JJ=1,NSV_AER
  PAERO(:,:,:,JJ) =  PAERO(:,:,:,JJ) / (ZDEN2MOL * PRHODREF(:,:,:))
  ENDDO

IF (LHOOK) CALL DR_HOOK('CH_AER_INIT',1,ZHOOK_HANDLE)
END SUBROUTINE CH_AER_INIT
