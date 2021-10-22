!     ######spl
     SUBROUTINE CH_AER_REALLFI_n(PSV, PCO, PRHODREF)
     USE PARKIND1, ONLY : JPRB
     USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!!   ############################################################
!!
!!    PURPOSE
!!    -------
!!    Initialise des  valeurs constantes pour les mineraux 
!!    Initialise des valeurs  proportionnelles au monoxyde de carbone pour OC et BC 
!!    Realise l'équilibre des moments à partir du sigma et du diametre moyen
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
!USE MODE_ll

USE MODD_CH_AEROSOL
USE MODD_CH_AERO_n
USE MODD_CH_MNHC_n , ONLY : CCH_SCHEME
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
REAL,   DIMENSION(:,:,:,:),    INTENT(INOUT) :: PSV
REAL,   DIMENSION(:,:,:),      INTENT(INOUT) :: PCO
REAL,   DIMENSION(:,:,:),      INTENT(IN) :: PRHODREF
!
!
!*      0.2    declarations local variables
!
REAL,DIMENSION(SIZE(PSV,1),SIZE(PSV,2),SIZE(PSV,3),NSP+NCARB+NSOA,JPMODE) :: ZCTOTA
REAL,DIMENSION(NSP+NCARB+NSOA) :: ZFAC, ZMI, ZRHOI
REAL,DIMENSION(SIZE(PSV,1),SIZE(PSV,2),SIZE(PSV,3),JPIN) :: ZM
REAL,DIMENSION(SIZE(PSV,1),SIZE(PSV,2),SIZE(PSV,3),JPMODE) :: ZRG, ZSIG, ZN
REAL,DIMENSION(SIZE(PSV,1),SIZE(PSV,2),SIZE(PSV,3)) :: ZSIGMA
REAL    :: ZPI, ZCOEFAEROBC, ZCOEFAEROOC, ZSUMAEROCO
INTEGER :: NJAERO, NIAERO
INTEGER :: JJ, JN  ! loop counter
REAL    :: ZDEN2MOL, ZVALBC, ZVALOC
REAL    :: ZINIRADIUSI, ZINIRADIUSJ
!
!-------------------------------------------------------------------------------
!
!*       1.     TRANSFER FROM GAS TO AEROSOL MODULE
!               ------------------------------------
!        1.1    initialisation 

REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('CH_AER_REALLFI_N',0,ZHOOK_HANDLE)
IF (CRGUNIT=="MASS") THEN
  ZINIRADIUSI = XINIRADIUSI * EXP(-3.*(LOG(XINISIGI))**2)
  ZINIRADIUSJ = XINIRADIUSJ * EXP(-3.*(LOG(XINISIGJ))**2)
ELSE
  ZINIRADIUSI = XINIRADIUSI 
  ZINIRADIUSJ = XINIRADIUSJ
END IF
!
! Moments index
NM0(1) = 1
NM3(1) = 2
NM6(1) = 3
NM0(2) = 4
NM3(2) = 5
NM6(2) = 6
!
! Default values of molar mass
ZMI(:) = 12.
ZMI(JP_AER_SO4)  = 98.
ZMI(JP_AER_NO3)  = 63.
ZMI(JP_AER_NH3)  = 17.
ZMI(JP_AER_H2O)  = 18.
!
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
ZRHOI(:) = 1.8e3
ZRHOI(JP_AER_H2O) = 1.0e3   ! water
!
PSV(:,:,:,:)  = MAX(PSV(:,:,:,:), 0.)
PCO(:,:,:)    = MAX(PCO(:,:,:), 10.E-9*PRHODREF(:,:,:))
!
! Special treatment for BC and OC (link to CO gaseous concentration)

IF (LINITPM) THEN
NIAERO=SIZE(PSV,1)-1
NJAERO=SIZE(PSV,2)-1
ZSUMAEROCO=SUM(PCO(2:NIAERO,2:NJAERO,2))/&
           ((NIAERO-1)*(NJAERO-1))


!ZVALBC=1.28058E-9  ! value in kg/m3 (escompte values)
!ZVALOC=2.304978E-9 ! value in kg/m3 (escompte values)
ZVALBC=1.E-9  ! value in kg/m3 (default values)
ZVALOC=2.E-9 ! value in kg/m3  (default values)
ZVALBC= ZVALBC *24.47 / 12. ! conversion into ppp 
ZVALOC= ZVALOC *24.47 / 12. ! conversion into ppp
ZCOEFAEROBC=ZVALBC/ZSUMAEROCO
ZCOEFAEROOC=ZVALOC/ZSUMAEROCO

PSV(:,:,:,JP_CH_BCi)=PCO(:,:,:) * ZCOEFAEROBC / 2.
PSV(:,:,:,JP_CH_BCj)=PCO(:,:,:) * ZCOEFAEROBC / 2.
PSV(:,:,:,JP_CH_OCi)=PCO(:,:,:) * ZCOEFAEROOC / 2.
PSV(:,:,:,JP_CH_OCj)=PCO(:,:,:) * ZCOEFAEROOC / 2.
END IF

!
!
ZDEN2MOL = 1E-6 * XAVOGADRO / XMD
!
! conversion into mol.cm-3
DO JJ=1, NSV_AER
  PSV(:,:,:,JJ) =  PSV(:,:,:,JJ) * ZDEN2MOL * PRHODREF(:,:,:)
ENDDO
!
ZPI = 2.*ASIN(1.)
DO JJ=1,NSP+NCARB+NSOA
  ZFAC(JJ)=(4./3.)*ZPI*ZRHOI(JJ)*1.e-9
ENDDO
!
!
!*       1.n    transfer aerosol mass from gas to aerosol variables
!               (and conversion of mol.cm-3 --> microgram/m3)
!
! aerosol phase
ZCTOTA(:,:,:,JP_AER_SO4,1) = PSV(:,:,:,JP_CH_SO4i)*ZMI(JP_AER_SO4)/6.0221367E+11
ZCTOTA(:,:,:,JP_AER_SO4,2) = PSV(:,:,:,JP_CH_SO4j)*ZMI(JP_AER_SO4)/6.0221367E+11
ZCTOTA(:,:,:,JP_AER_NO3,1) = PSV(:,:,:,JP_CH_NO3i)*ZMI(JP_AER_NO3)/6.0221367E+11
ZCTOTA(:,:,:,JP_AER_NO3,2) = PSV(:,:,:,JP_CH_NO3j)*ZMI(JP_AER_NO3)/6.0221367E+11
ZCTOTA(:,:,:,JP_AER_NH3,1) = PSV(:,:,:,JP_CH_NH3i)*ZMI(JP_AER_NH3)/6.0221367E+11
ZCTOTA(:,:,:,JP_AER_NH3,2) = PSV(:,:,:,JP_CH_NH3j)*ZMI(JP_AER_NH3)/6.0221367E+11
!
! water
ZCTOTA(:,:,:,JP_AER_H2O,1) = PSV(:,:,:,JP_CH_H2Oi)*ZMI(JP_AER_H2O)/6.0221367E+11
ZCTOTA(:,:,:,JP_AER_H2O,2) = PSV(:,:,:,JP_CH_H2Oj)*ZMI(JP_AER_H2O)/6.0221367E+11
!
! primary organic carbon
ZCTOTA(:,:,:,JP_AER_OC,1) = PSV(:,:,:,JP_CH_OCi)*ZMI(JP_AER_OC)/6.0221367E+11
ZCTOTA(:,:,:,JP_AER_OC,2) = PSV(:,:,:,JP_CH_OCj)*ZMI(JP_AER_OC)/6.0221367E+11
!
! primary black carbon
ZCTOTA(:,:,:,JP_AER_BC,1) = PSV(:,:,:,JP_CH_BCi)*ZMI(JP_AER_BC)/6.0221367E+11
ZCTOTA(:,:,:,JP_AER_BC,2) = PSV(:,:,:,JP_CH_BCj)*ZMI(JP_AER_BC)/6.0221367E+11
!
IF (NSOA .EQ. 10) THEN
  ZCTOTA(:,:,:,JP_AER_SOA1,1) = PSV(:,:,:,JP_CH_SOA1i)*ZMI(JP_AER_SOA1)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA1,2) = PSV(:,:,:,JP_CH_SOA1j)*ZMI(JP_AER_SOA1)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA2,1) = PSV(:,:,:,JP_CH_SOA2i)*ZMI(JP_AER_SOA2)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA2,2) = PSV(:,:,:,JP_CH_SOA2j)*ZMI(JP_AER_SOA2)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA3,1) = PSV(:,:,:,JP_CH_SOA3i)*ZMI(JP_AER_SOA3)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA3,2) = PSV(:,:,:,JP_CH_SOA3j)*ZMI(JP_AER_SOA3)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA4,1) = PSV(:,:,:,JP_CH_SOA4i)*ZMI(JP_AER_SOA4)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA4,2) = PSV(:,:,:,JP_CH_SOA4j)*ZMI(JP_AER_SOA4)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA5,1) = PSV(:,:,:,JP_CH_SOA5i)*ZMI(JP_AER_SOA5)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA5,2) = PSV(:,:,:,JP_CH_SOA5j)*ZMI(JP_AER_SOA5)/6.0221367E+11

  ZCTOTA(:,:,:,JP_AER_SOA6,1) = PSV(:,:,:,JP_CH_SOA6i)*ZMI(JP_AER_SOA6)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA6,2) = PSV(:,:,:,JP_CH_SOA6j)*ZMI(JP_AER_SOA6)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA7,1) = PSV(:,:,:,JP_CH_SOA7i)*ZMI(JP_AER_SOA7)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA7,2) = PSV(:,:,:,JP_CH_SOA7j)*ZMI(JP_AER_SOA7)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA8,1) = PSV(:,:,:,JP_CH_SOA8i)*ZMI(JP_AER_SOA8)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA8,2) = PSV(:,:,:,JP_CH_SOA8j)*ZMI(JP_AER_SOA8)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA9,1) = PSV(:,:,:,JP_CH_SOA9i)*ZMI(JP_AER_SOA9)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA9,2) = PSV(:,:,:,JP_CH_SOA9j)*ZMI(JP_AER_SOA9)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA10,1) = PSV(:,:,:,JP_CH_SOA10i)*ZMI(JP_AER_SOA10)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA10,2) = PSV(:,:,:,JP_CH_SOA10j)*ZMI(JP_AER_SOA10)/6.0221367E+11
END IF



!
!*       1.1    calculate moment 3 from total aerosol mass
!
ZM(:,:,:,2) = 0.
ZM(:,:,:,5) = 0.
DO JJ = 1,NSP+NCARB+NSOA
  ZM(:,:,:,2) = ZM(:,:,:,2)+ZCTOTA(:,:,:,JJ,1)/ZFAC(JJ)
  ZM(:,:,:,5) = ZM(:,:,:,5)+ZCTOTA(:,:,:,JJ,2)/ZFAC(JJ)
ENDDO
!
!
!*       1.2    calculate moment 0 from dispersion and mean radius
!
ZM(:,:,:,1)= ZM(:,:,:,2) / &
            ((ZINIRADIUSI**3)*EXP(4.5 * (LOG(XINISIGI))**2))
ZM(:,:,:,4)= ZM(:,:,:,5) / &
            ((ZINIRADIUSJ**3)*EXP(4.5 * (LOG(XINISIGJ))**2))
!
!*       1.3    calculate moment 6 from dispersion and mean radius
!
ZM(:,:,:,3) = ZM(:,:,:,1) * (ZINIRADIUSI**6) *EXP(18 *(LOG(XINISIGI))**2)
ZM(:,:,:,6) = ZM(:,:,:,4) * (ZINIRADIUSJ**6) *EXP(18 *(LOG(XINISIGJ))**2)
!
PSV(:,:,:,JP_CH_M0i) = ZM(:,:,:,1) * 1E-6 
PSV(:,:,:,JP_CH_M0j) = ZM(:,:,:,4) * 1E-6
IF (LVARSIGI) PSV(:,:,:,JP_CH_M6i) = ZM(:,:,:,3) 
IF (LVARSIGJ) PSV(:,:,:,JP_CH_M6j) = ZM(:,:,:,6)
!
DO JN=1,JPMODE
  IF (JN .EQ. 1) THEN
  !
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
  !
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
!
  ZSIG(:,:,:,JN) = ZSIGMA(:,:,:)
  ZN(:,:,:,JN) = ZM(:,:,:,NM0(JN))
!
  ZSIGMA(:,:,:)=LOG(ZSIG(:,:,:,JN))**2
!
  ZRG(:,:,:,JN)=(ZM(:,:,:,NM3(JN))/ZN(:,:,:,JN))**(1./3.)*EXP(-1.5*ZSIGMA(:,:,:))
!
ENDDO
!
!conversion into ppp
DO JJ=1,NSV_AER
  PSV(:,:,:,JJ) =  PSV(:,:,:,JJ) /  (ZDEN2MOL*PRHODREF(:,:,:)) 
ENDDO
!
!
IF (LHOOK) CALL DR_HOOK('CH_AER_REALLFI_N',1,ZHOOK_HANDLE)
END SUBROUTINE CH_AER_REALLFI_n
