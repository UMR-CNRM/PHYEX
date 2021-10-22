!     ######spl
     SUBROUTINE CH_AER_MPMPO(PTOT, PTOTG, PTEMP, PRH, PLWC, PPROTON, PTOTNEW,PTOTGNEW,PSOLORG)
     USE PARKIND1, ONLY : JPRB
     USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!###########################################################################################
!!
!!   PURPOSE
!!   -------
!!   solve the organic thermodynamic balance using MPMP0 (Griffin et al. 2003, J. Atm. Chem, 44, 171-190)) ,
!!   if use CACM our ReLACS2 chemical scheme
!!   Code updated from Griffin 2005.
!!
!!    REFERENCE
!!    ---------
!!    none
!!
!!    AUTHOR
!!    ------
!!    Alf Grini / P. Tulet (CNRM)
!!
!!    MODIFICATIONS
!!    -------------
!!
!!    EXTERNAL
!!    --------
!!    None
!-------------------------------------------------------------------------------!
USE MODD_CH_AEROSOL
USE MODD_GLO
USE MODE_OAMAIN
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!
REAL, DIMENSION(:,:),   INTENT(IN)    :: PTOT  ![ug/m3] aerosol conc
REAL, DIMENSION(:,:),   INTENT(IN)    :: PTOTG ![ug/m3] gas concentration
REAL, DIMENSION(:),     INTENT(IN)    :: PTEMP, PRH ![K/-] temp and RH
REAL, DIMENSION(:),     INTENT(IN)    :: PLWC  ![ug/m3] liquid water
REAL, DIMENSION(:),     INTENT(IN)    :: PPROTON ![mole/g_{water]] H+
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PTOTNEW ![ug/m3] new aerosol conc
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PTOTGNEW ![ug/m3] new gas conc
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PSOLORG ! [%] Solubility of SOA (fraction)
!
!*       0.2   Declarations of local variables
!
INTEGER :: JI
REAL, DIMENSION(SIZE(PTOT,1),NSOA)   :: ZCPART     ![ug/m3] total ORG+AQ+GAS SOA
REAL, DIMENSION(SIZE(PTOT,1),NSOA)   :: ZAERO_ORG   ![ug/m3] organic phase SOA
REAL, DIMENSION(SIZE(PTOT,1),NAAERO) :: ZAERO_AQ    ![ug/m3] aquous phase SOA
REAL, DIMENSION(SIZE(PTOT,1),NSOA)   :: ZPARTORG    ![ug/m3] particle phase organic 
REAL, DIMENSION(SIZE(PTOT,1),NSOA)   :: ZGASORG    ![ug/m3] gas phase organic
REAL, DIMENSION(SIZE(PTOT,1),NBSPOA) :: ZCPT        ![ug/m3] Primary organic aerosol
REAL, DIMENSION(SIZE(PTOT,1)) :: ZORGANION         ![mole/m3] negative charge associated with organics
REAL, DIMENSION(SIZE(PTOT,1)) :: ZDELLWC           ![ug/m3] liquid water assiciated with organics
REAL, DIMENSION(SIZE(PTOT,1)) :: ZHCONC            ![mole/kg_{water}] proton concentration
!-------------------------------------------------------------------------------
!
!cpart = total concentration of secondary species (gas + aero), microgram/m3
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('CH_AER_MPMPO',0,ZHOOK_HANDLE)
DO JI=NSP+NCARB+1,NSP+NCARB+NSOA
   ZCPART(:,JI-NSP-NCARB) = PTOT(:,JI) + PTOTG(:,JI)
END DO
!
!Distribute POA equally between the 8 species..?
DO JI=1,NBSPOA
   ZCPT(:,JI) = 1./dble(NBSPOA)*PTOT(:,JP_AER_OC)
ENDDO

!Save LWC in a local variable (why??)
!ZLWC(:) = PLWC(:)

! Compute relative humidity (why??)
!ZRH(:)= PRH(:)

!Compute H+ concentration in right unit
ZHCONC(:) = 1E3 * PPROTON(:) ! conversion in mol H+ / kg solvent
!
!Calling the mpmpo model

CALL MPMPO(        &
     PTEMP         &  !I [K] temperature
     ,PRH          &  !I [0-1] Relative humidity
     ,ZCPART       &  !I [ug/m3] total aerosol (aq+org+gas)
     ,ZCPT         &  !I [ug/m3] primary organic aerosol
     ,ZGASORG      &  !O [ug/m3] gas phase concentration
     ,ZAERO_ORG    &  !O [ug/m3] aerosol phase concentration of organic phase
     ,ZAERO_AQ     &  !O [ug/m3] aerosol phase concentrations of aquous phase
     ,ZPARTORG     &  !O [ug/m3] New total aerosol phase concentrations
     ,PLWC         &  !I [ug/m3] liquid water content already in the aerosols
     ,ZHCONC       &  !I [mole/kg_{water}] proton concentration
     ,ZDELLWC      &  !I [ug/m3] additional liquid water content because of organics
     ,ZORGANION    &  !O [mole/m3] negative charge from organic anions
     ,PSOLORG      & !IO [%] Solubility of SOA (fraction)
     )
!
! New SOA
DO JI=1,NSOA
   PTOTNEW(:,JI+NSP+NCARB) = ZPARTORG(:,JI) 
   PTOTGNEW(:,JI+NSP+NCARB) = ZGASORG(:,JI) 
END DO
!
! New LWC
PTOTNEW(:,JP_AER_H2O) = PLWC(:) + ZDELLWC(:)
!
IF (LHOOK) CALL DR_HOOK('CH_AER_MPMPO',1,ZHOOK_HANDLE)
END SUBROUTINE CH_AER_MPMPO
