!     ######spl
     SUBROUTINE SALT_VELGRAV(PSIG, PRG, PTHT, PABST, PRHODREF, PRHOP, PMU, PVGK,PDPK, PVGG)
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
!!   P. Tulet (CNRM/GMEI)
!!
!!   MODIFICATIONS
!!    -------------
!!
! Entry variables:
!
! PM(IN)       -Array of moments
!
!*************************************************************
! Exit variables:
!
! PFSED(IN)  -Array of moment variation due to dry deposition
!
!*************************************************************
! Variables used during the deposition velocity calculation
! 
! PDPK       -Polydisperse diffusivity (m2/s)
! PVGK       -Polydisperse settling velocity of the kth moment (m/s)
!************************************************************
!!
!-----------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
!
USE MODD_SALT
USE MODD_CSTS_SALT
!
IMPLICIT NONE
!
!*      0.1  declarations of arguments
!
REAL,                     INTENT(IN) :: PRHOP
REAL, DIMENSION(:,:,:,:), INTENT(IN) :: PSIG, PRG
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PTHT, PABST, PRHODREF
REAL, DIMENSION(:,:,:,:), INTENT(OUT) :: PVGK,PDPK
REAL, DIMENSION(:,:,:),   INTENT(OUT) :: PMU
REAL, DIMENSION(:,:,:,:), INTENT(OUT) :: PVGG
!
!*      0.2  declaration of local variables
!
REAL, DIMENSION(SIZE(PSIG,1),SIZE(PSIG,2),SIZE(PSIG,3)) :: ZTEMP,ZLAMBDA
!
REAL, DIMENSION(SIZE(PSIG,1),SIZE(PSIG,2),SIZE(PSIG,3)) :: ZRG,ZLN2S
!
REAL, DIMENSION(SIZE(PSIG,1),SIZE(PSIG,2),SIZE(PSIG,3)) :: ZKNG
!
REAL, DIMENSION(SIZE(PSIG,1),SIZE(PSIG,2),SIZE(PSIG,3),NMODE_SLT) :: ZDPG
!
!
REAL, PARAMETER :: gasmw=28.9644d0
REAL :: ZK
!
INTEGER :: JI,JJ
!
!-----------------------------------------------------------------
!temperature
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('SALT_VELGRAV',0,ZHOOK_HANDLE)
ZTEMP(:,:,:)=PTHT(:,:,:)*(PABST(:,:,:)/XP00)**(XRD/XCPD)
!
! Sutherland's equation for viscosity
PMU(:,:,:)=1.8325d-5*416.16/(ZTEMP(:,:,:)+120)*(ZTEMP(:,:,:)/296.16)*SQRT(ZTEMP(:,:,:)/296.16)
!
! Mean free path (Seinfeld and Pandis p455)
ZLAMBDA(:,:,:)=PMU(:,:,:)/PRHODREF(:,:,:)*SQRT(1.89d-4*gasmw/ZTEMP(:,:,:))*1.e6
!
!
DO JI=1,NMODE_SLT
  ZRG(:,:,:)=PRG(:,:,:,JI) * 1E-6 
  ZLN2S(:,:,:)=LOG(PSIG(:,:,:,JI))**2 
  !
  ZKNG(:,:,:)=ZLAMBDA(:,:,:) / PRG(:,:,:,JI) 
  !
  PVGG(:,:,:,JI)= 2.*XG*PRHOP*ZRG(:,:,:)**2 /(9.*PMU(:,:,:))
  ZDPG(:,:,:,JI)=XBOLTZ*ZTEMP(:,:,:)/ (6.*XPI* ZRG(:,:,:)*PMU(:,:,:))
  ! 
  ! 
  DO JJ=0,2
    ZK=REAL(3*JJ)
    PDPK(:,:,:,3*JI+JJ-2)=ZDPG(:,:,:,JI)*(exp((-2.*ZK+1.)/2.*ZLN2S(:,:,:))+1.246*ZKNG(:,:,:)*&
              exp((-4.*ZK+4)/2.*ZLN2S(:,:,:)))
    !
    PVGK(:,:,:,3*JI+JJ-2)=PVGG(:,:,:,JI)*&
    (exp((4.*ZK+4.)/2.*ZLN2S(:,:,:)) + 1.246*ZKNG(:,:,:)* exp((2.*ZK+1.)/2.*ZLN2S(:,:,:)))
  ENDDO
  !
ENDDO
!
!
IF (LHOOK) CALL DR_HOOK('SALT_VELGRAV',1,ZHOOK_HANDLE)
END SUBROUTINE SALT_VELGRAV
