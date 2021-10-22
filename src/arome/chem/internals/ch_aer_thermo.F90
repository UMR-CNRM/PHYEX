!     ######spl
SUBROUTINE CH_AER_THERMO(PAER,PRH, PDENAIR, PPRESSURE, PTEMP, PRC)
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!
!
!!    PURPOSE
!!    -------
! Make the thermodynamic equilibrium calculation
! Give gas equilibrium concentration for semi-volatile
! species NH3 and HNO3 ... at that time.
! For SO4, the equilibrium concentration is 0.0

!!    AUTHOR
!!    ------
!!      F. Cousin          * LA *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    05/04/02
!!      Modifications  15/03/03 (P.Tulet)  update to log-normal model
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! INPUT
! PAER(:,1) :: H2SO4    in micrograms / m**3
! PAER(:,2) :: NH3(g)   in micrograms / m**3
! PAER(:,3) :: HNO3(g)  in micrograms / m**3
! PAER(:,4) :: H2O(a)   in micrograms / m**3
! PAER(:,5) :: NO3(a)   in micrograms / m**3
! PAER(:,6) :: NH4(a)   in micrograms / m**3
! PTEMP   : temperature in K
! PPRESSURE   : pressure = RHO * T * R  (in Pa) assuming 20.95 vol% O2
! PPKM   : n_molec (moelc./cm3):      M = 1E-3*RHO(kg/m3) * Navo / m_mol
! PPKH2  : molecular weight of H2O:   m_mol^H2O = 18.0 g/mol
!                        ==>   m_mol^air / m_mol^H2O = 1.6
!
! OUTPUT
!   zout(:,1) :: waterc : ratio [water content]/[mineral species content]
!   zout(:,2)*cnh3(:)  : NH3 gas equilibrimum concentrations (kg/m3)
!   zout(:,3)*cno3(:)  : NO3 gas equilibrimum concentrations (kg/m3)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
USE MODI_CH_AER_INTERMIN
USE MODD_CH_AEROSOL
IMPLICIT NONE
!
!*       0.1   declarations of arguments
REAL, DIMENSION(:,:), INTENT(INOUT) :: PAER
REAL, DIMENSION(:),   INTENT(IN)    :: PRH, PDENAIR, PPRESSURE, PTEMP, PRC

!
!*       0.2   declarations of local variables
DOUBLE PRECISION, DIMENSION(SIZE(PAER,1)) :: waterc,chno3,ch2so4,cnh3
DOUBLE PRECISION, DIMENSION(SIZE(PAER,1)) :: TEMPE, ZTOTNH3, ZTOTNO3
DOUBLE PRECISION, DIMENSION(SIZE(PAER,1),3) :: zout
!
!...........PARAMETERS and their descriptions:

REAL, PARAMETER ::  ZMWH2O = 18.0           ! molecular weight for water
REAL, PARAMETER ::  ZMWNO3 = 62.0049        ! molecular weight for NO3
REAL, PARAMETER ::  ZMH2SO4 = 98.07354      ! molecular weight for H2SO4
REAL, PARAMETER ::  ZMWNH3 = 17.03061       ! molecular weight for NH3
! 
! Conversion micrograms / m3 in  kg/m3
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('CH_AER_THERMO',0,ZHOOK_HANDLE)
PAER(:,:) = PAER(:,:) * 1E-9
!
IF((NSP-1).NE.0) THEN
!
!
!
! Compute total mass concentration for mineral species
!
!!!!ch2so4(:) = PAER(:,1)  + conc(:,numesoc(ih2so4-npri))*fac1*conmw(ih2so4-npri)
ZTOTNO3(:) = PAER(:,5) + PAER(:,3) 
ZTOTNH3(:) = PAER(:,6) + PAER(:,2)
ch2so4(:) = PAER(:,1) ! all so4 is in the aerosol phase: to be update in the future
chno3(:) = ZTOTNO3(:) 
cnh3(:)  = ZTOTNH3(:)

! Temperature
TEMPE=PTEMP(:)

!print*,'avant CH_AER_INTERMIN'
!print*,'ch2so4=', ch2so4(1:5)
!print*,'cnh3=', cnh3(1:5)
!print*,'cno3=', chno3(1:5)
!
!
CALL CH_AER_INTERMIN(PRH,TEMPE,ch2so4,cnh3,chno3,zout,SIZE(PAER,1))
!
!print*,'apres CH_AER_INTERMIN'
!print*,'zout(:1:5,1)=',zout(1:5,1)
!print*,'zout(:1:5,2)=',zout(1:5,2)
!print*,'zout(:1:5,3)=',zout(1:5,3)
waterc(:) = zout(:,1)
PAER(:,2) = zout(:,2)*cnh3(:) 
PAER(:,3) = zout(:,3)*chno3(:)
PAER(:,5) = ZTOTNO3(:) - PAER(:,3)
PAER(:,6) = ZTOTNH3(:) - PAER(:,2)
PAER(:,4) = waterc(:) * ZMWH2O * &
            & (PAER(:,1) / ZMH2SO4 + PAER(:,5) / ZMWNO3 + PAER(:,6) / ZMWNH3)

! 
! Conversion kg/m3 in  micrograms / m3 

PAER(:,:) = PAER(:,:) * 1E9
!
!
ENDIF
!
IF (LHOOK) CALL DR_HOOK('CH_AER_THERMO',1,ZHOOK_HANDLE)
RETURN
END SUBROUTINE CH_AER_THERMO
