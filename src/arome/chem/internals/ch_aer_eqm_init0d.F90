!     ######spl
     SUBROUTINE CH_AER_EQM_INIT0d(PMI, PAERO, PM3D, PRHOP3D, PSIG3D, PRG3D, &
                             PN3D, PCTOTA)
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
USE MODD_CH_M9, ONLY : CNAMES
USE MODD_CH_AERO_n
USE MODD_CH_MNHC_n

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
REAL,   DIMENSION(:,:),    INTENT(INOUT)   :: PAERO, PMI
REAL,   DIMENSION(:,:),    INTENT(INOUT)   :: PM3D, PRHOP3D, PSIG3D, PRG3D, PN3D
REAL,   DIMENSION(:,:,:),  INTENT(INOUT)   :: PCTOTA


!
!
!*      0.2    declarations local variables
!
REAL,DIMENSION(1,NSP+NCARB+NSOA,JPMODE) :: ZCCTOT
REAL,DIMENSION(1) :: ZSUM
REAL,DIMENSION(1) :: ZSIGMA
INTEGER :: JN, JJ ! loop counter
!
!-------------------------------------------------------------------------------
!

!*       1.     TRANSFER FROM GAS TO AEROSOL MODULE
!               ------------------------------------
!        1.1    initialisation 
! Index gas scheme <=> Index Orilam

REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('CH_AER_EQM_INIT0D',0,ZHOOK_HANDLE)
DO JJ=1,SIZE(CNAMES)
IF (CNAMES(JJ) == "CO")  JP_CH_CO   = JJ
END DO


IF (CORGANIC == TRIM("MPMPO") .OR. CORGANIC == TRIM("PUN") .OR. CORGANIC == TRIM("EQSAM2")) THEN
  IF ((CCH_SCHEME .NE. TRIM("CACM")) .AND. (CCH_SCHEME .NE. TRIM("RELACS2"))) THEN
  print*, '**********************************************'
  print*, 'WARNING : NO SOA !!!!'
  print*, 'YOU WANT TO USE SOA GAS PARTICLE BALANCE'
  print*, 'BUT THE SCHEME NEED TO BE CACM or RELACS 2'
  print*, 'CORGANIC HAS BEEN SET TO NONE'
  print*, 'OTHERWISE COMPILE THE CORRECT SCHEME BEFORE'
  print*, '**********************************************'
  CORGANIC = "NONE"
  END IF
END IF

IF (.NOT.(ALLOCATED(XRHOI))) ALLOCATE(XRHOI(NSP+NSOA+NCARB))
IF (.NOT.(ALLOCATED(XFAC))) ALLOCATE(XFAC(NSP+NSOA+NCARB))


! Moments index
 NM0(1) = 1 
 NM3(1) = 2
 NM6(1) = 3 
 NM0(2) = 4 
 NM3(2) = 5 
 NM6(2) = 6 

! Aerosol Density
! Cf Ackermann (all to black carbon except water)

 XRHOI(:) = 1.8e3
 XRHOI(JP_AER_H2O) = 1.0e3   ! water

!
DO JJ=1,NSP+NCARB+NSOA
  XFAC(JJ)=(4./3.)*3.14292654*XRHOI(JJ)*1.e-9
ENDDO

!
!
!*       1.n    transfer aerosol mass from gas to aerosol variables
!               (and conversion of part/part --> microgram/m3)
!

!
! mineral phase
  PCTOTA(:,JP_AER_SO4,1) = PAERO(:,JP_CH_SO4i)*PMI(:,JP_AER_SO4)/6.0221367E+11
  PCTOTA(:,JP_AER_SO4,2) = PAERO(:,JP_CH_SO4j)*PMI(:,JP_AER_SO4)/6.0221367E+11

  PCTOTA(:,JP_AER_NO3,1) = PAERO(:,JP_CH_NO3i)*PMI(:,JP_AER_NO3)/6.0221367E+11
  PCTOTA(:,JP_AER_NO3,2) = PAERO(:,JP_CH_NO3j)*PMI(:,JP_AER_NO3)/6.0221367E+11

  PCTOTA(:,JP_AER_NH3,1) = PAERO(:,JP_CH_NH3i)*PMI(:,JP_AER_NH3)/6.0221367E+11
  PCTOTA(:,JP_AER_NH3,2) = PAERO(:,JP_CH_NH3j)*PMI(:,JP_AER_NH3)/6.0221367E+11

! water
  PCTOTA(:,JP_AER_H2O,1) = PAERO(:,JP_CH_H2Oi)*PMI(:,JP_AER_H2O)/6.0221367E+11
  PCTOTA(:,JP_AER_H2O,2) = PAERO(:,JP_CH_H2Oj)*PMI(:,JP_AER_H2O)/6.0221367E+11

!
! primary organic carbon
  PCTOTA(:,JP_AER_OC,1) = PAERO(:,JP_CH_OCi)*PMI(:,JP_AER_OC)/6.0221367E+11
  PCTOTA(:,JP_AER_OC,2) = PAERO(:,JP_CH_OCj)*PMI(:,JP_AER_OC)/6.0221367E+11

! primary black carbon
  PCTOTA(:,JP_AER_BC,1) = PAERO(:,JP_CH_BCi)*PMI(:,JP_AER_BC)/6.0221367E+11
  PCTOTA(:,JP_AER_BC,2) = PAERO(:,JP_CH_BCj)*PMI(:,JP_AER_BC)/6.0221367E+11

!
  PCTOTA(:,JP_AER_SOA1,1) = PAERO(:,JP_CH_SOA1i)*PMI(:,JP_AER_SOA1)/6.0221367E+11
  PCTOTA(:,JP_AER_SOA1,2) = PAERO(:,JP_CH_SOA1j)*PMI(:,JP_AER_SOA1)/6.0221367E+11
  PCTOTA(:,JP_AER_SOA2,1) = PAERO(:,JP_CH_SOA2i)*PMI(:,JP_AER_SOA2)/6.0221367E+11
  PCTOTA(:,JP_AER_SOA2,2) = PAERO(:,JP_CH_SOA2j)*PMI(:,JP_AER_SOA2)/6.0221367E+11
  PCTOTA(:,JP_AER_SOA3,1) = PAERO(:,JP_CH_SOA3i)*PMI(:,JP_AER_SOA3)/6.0221367E+11
  PCTOTA(:,JP_AER_SOA3,2) = PAERO(:,JP_CH_SOA3j)*PMI(:,JP_AER_SOA3)/6.0221367E+11
  PCTOTA(:,JP_AER_SOA4,1) = PAERO(:,JP_CH_SOA4i)*PMI(:,JP_AER_SOA4)/6.0221367E+11
  PCTOTA(:,JP_AER_SOA4,2) = PAERO(:,JP_CH_SOA4j)*PMI(:,JP_AER_SOA4)/6.0221367E+11
  PCTOTA(:,JP_AER_SOA5,1) = PAERO(:,JP_CH_SOA5i)*PMI(:,JP_AER_SOA5)/6.0221367E+11
  PCTOTA(:,JP_AER_SOA5,2) = PAERO(:,JP_CH_SOA5j)*PMI(:,JP_AER_SOA5)/6.0221367E+11
  PCTOTA(:,JP_AER_SOA6,1) = PAERO(:,JP_CH_SOA6i)*PMI(:,JP_AER_SOA6)/6.0221367E+11
  PCTOTA(:,JP_AER_SOA6,2) = PAERO(:,JP_CH_SOA6j)*PMI(:,JP_AER_SOA6)/6.0221367E+11

  PCTOTA(:,JP_AER_SOA6,1) = PAERO(:,JP_CH_SOA6i)*PMI(:,JP_AER_SOA6)/6.0221367E+11
  PCTOTA(:,JP_AER_SOA6,2) = PAERO(:,JP_CH_SOA6j)*PMI(:,JP_AER_SOA6)/6.0221367E+11
  PCTOTA(:,JP_AER_SOA7,1) = PAERO(:,JP_CH_SOA7i)*PMI(:,JP_AER_SOA7)/6.0221367E+11
  PCTOTA(:,JP_AER_SOA7,2) = PAERO(:,JP_CH_SOA7j)*PMI(:,JP_AER_SOA7)/6.0221367E+11
  PCTOTA(:,JP_AER_SOA8,1) = PAERO(:,JP_CH_SOA8i)*PMI(:,JP_AER_SOA8)/6.0221367E+11
  PCTOTA(:,JP_AER_SOA8,2) = PAERO(:,JP_CH_SOA8j)*PMI(:,JP_AER_SOA8)/6.0221367E+11
  PCTOTA(:,JP_AER_SOA9,1) = PAERO(:,JP_CH_SOA9i)*PMI(:,JP_AER_SOA9)/6.0221367E+11
  PCTOTA(:,JP_AER_SOA9,2) = PAERO(:,JP_CH_SOA9j)*PMI(:,JP_AER_SOA9)/6.0221367E+11
  PCTOTA(:,JP_AER_SOA10,1) = PAERO(:,JP_CH_SOA10i)*PMI(:,JP_AER_SOA10)/6.0221367E+11
  PCTOTA(:,JP_AER_SOA10,2) = PAERO(:,JP_CH_SOA10j)*PMI(:,JP_AER_SOA10)/6.0221367E+11


!
!*       1.1    calculate moment 3 from mass
    
    PM3D(:,2) = 0.
    PM3D(:,5) = 0.
    PCTOTA(:,:,:) = MAX(PCTOTA(:,:,:), 0.)
    DO JJ = 1,NSP+NCARB+NSOA
    PM3D(:,2) = PM3D(:,2)+PCTOTA(:,JJ,1)/XFAC(JJ)
    PM3D(:,5) = PM3D(:,5)+PCTOTA(:,JJ,2)/XFAC(JJ)
    ENDDO
!
!
!*       1.2    calculate moment 0 from dispersion and mean radius
    PM3D(:,1)= PM3D(:,2) / &
               ((XINIRADIUSI**3)*EXP(4.5 * (LOG(XINISIGI))**2))

IF (ANY(PM3D(:,1) < XN0IMIN)) THEN

print*, 'FATAL ERROR '
print*, 'COMPATIBILITY ERROR: Initialization of particle number mode I < XN0IMIN '
print*, ' MINIMAL NUMBER PARTICLE BY m3 is ', MINVAL(PM3D(:,1)),&
'located at ',MINLOC(PM3D(:,1))
print*, 'PLEASE CHANGE MASS OR XN0IMIN INITIALIZATION '
STOP
END IF
    PM3D(:,4)= PM3D(:,5) / &
               ((XINIRADIUSJ**3)*EXP(4.5 * (LOG(XINISIGJ))**2))

IF (ANY(PM3D(:,4) < XN0JMIN)) THEN
print*, 'FATAL ERROR '
print*, 'COMPATIBILITY ERROR: Initialization of particle number mode J < XN0JMIN '
print*, ' MINIMAL NUMBER PARTICLE BY m3 is ',MINVAL(PM3D(:,4)),&
'located at ',MINLOC(PM3D(:,4))
print*, 'PLEASE CHANGE MASS OR XN0JMIN INITIALIZATION '
STOP
END IF

!*       1.3    calculate moment 6 from dispersion and mean radius
   PM3D(:,3) = PM3D(:,1) * (XINIRADIUSI**6) *EXP(18 *(LOG(XINISIGI))**2)
   PM3D(:,6) = PM3D(:,4) * (XINIRADIUSJ**6) *EXP(18 *(LOG(XINISIGJ))**2)



!
!**********************************************
! Calcul de XRHOP3D
!**********************************************

PRHOP3D(:,:)=0.
DO JN=1,JPMODE
  ZSUM(:)=0.
  DO JJ=1,NSP+NCARB+NSOA
   ZSUM(:)=ZSUM(:)+PCTOTA(:,JJ,JN)/XRHOI(JJ)
  ENDDO
  DO JJ=1,NSP+NCARB+NSOA
  ZCCTOT(:,JJ,JN)=PCTOTA(:,JJ,JN)/XRHOI(JJ)/ZSUM(:)
  PRHOP3D(:,JN)=PRHOP3D(:,JN)+ZCCTOT(:,JJ,JN)*XRHOI(JJ)
  ENDDO
ENDDO

DO JN=1,JPMODE

  IF (JN .EQ. 1) THEN

    IF (LVARSIGI) THEN ! variable dispersion for mode 1

      ZSIGMA(:)=PM3D(:,NM3(JN))**2/(PM3D(:,NM0(JN))*PM3D(:,NM6(JN)))
      ZSIGMA(:)=MIN(1-1E-10,ZSIGMA(:))
      ZSIGMA(:)=MAX(1E-10,ZSIGMA(:))
      ZSIGMA(:)= LOG(ZSIGMA(:))
      ZSIGMA(:)= EXP(1./3.*SQRT(-ZSIGMA(:)))
      WHERE (ZSIGMA(:) > XSIGIMAX)
      ZSIGMA(:) =  XSIGIMAX
      END WHERE
      WHERE (ZSIGMA(:) < XSIGIMIN)
      ZSIGMA(:) =  XSIGIMIN
      END WHERE

    ELSE ! fixed dispersion for mode 1
      ZSIGMA(:) = XINISIGI
    END IF
  END IF

!
  IF (JN .EQ. 2) THEN

    IF (LVARSIGJ) THEN ! variable dispersion for mode 2

      ZSIGMA(:)=PM3D(:,NM3(JN))**2/(PM3D(:,NM0(JN))*PM3D(:,NM6(JN)))
      ZSIGMA(:)=MIN(1-1E-10,ZSIGMA(:))
      ZSIGMA(:)=MAX(1E-10,ZSIGMA(:))
      ZSIGMA(:)= LOG(ZSIGMA(:))
      ZSIGMA(:)= EXP(1./3.*SQRT(-ZSIGMA(:)))
      WHERE (ZSIGMA(:) > XSIGJMAX)
      ZSIGMA(:) =  XSIGJMAX
      END WHERE
      WHERE (ZSIGMA(:) < XSIGJMIN)
      ZSIGMA(:) =  XSIGJMIN
      END WHERE

    ELSE ! fixed dispersion for mode 2
      ZSIGMA(:) = XINISIGJ
    END IF
  END IF
!
!*       1.4    calculate modal parameters from moments
PSIG3D(:,JN) = ZSIGMA(:)
PN3D(:,JN) = PM3D(:,NM0(JN))

ZSIGMA(:)=LOG(PSIG3D(:,JN))**2

PRG3D(:,JN)=(PM3D(:,NM3(JN))/PN3D(:,JN))**(1./3.)*EXP(-1.5*ZSIGMA(:))

PM3D(:,NM6(JN))=PN3D(:,JN)*PRG3D(:,JN)**6*EXP(18.*ZSIGMA(:))
!
PSIG3D(:,JN)=LOG(PSIG3D(:,JN))
ENDDO
!
!
PAERO(:,JP_CH_M0i) = PM3D(:,1) * 1E-6 
PAERO(:,JP_CH_M0j) = PM3D(:,4) * 1E-6
IF (LVARSIGI) PAERO(:,JP_CH_M6i) = PM3D(:,3) 
IF (LVARSIGJ) PAERO(:,JP_CH_M6j) = PM3D(:,6)

!
!
IF (LHOOK) CALL DR_HOOK('CH_AER_EQM_INIT0D',1,ZHOOK_HANDLE)
END SUBROUTINE CH_AER_EQM_INIT0d
