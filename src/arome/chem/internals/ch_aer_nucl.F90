!     ######spl
SUBROUTINE CH_AER_NUCL(ZRH,ZT,ZCONC,ZJ,ZAL,KVECNPT)
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!###########################################################
!
!!                   
!!                       
!!
!!    PURPOSE
!!    -------
!!      
!!    compute nucleation rate for binary sulfate/H2O 
!!    (Kulmala 1998)
!!
!!    AUTHOR
!!    ------
!!      F.Cousin      * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!     P. Tulet       05/01/04 
!----------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
! 
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
REAL, DIMENSION(:), INTENT(INOUT) :: ZJ,ZAL
REAL, DIMENSION(:), INTENT(IN)    :: ZRH,ZT
REAL, DIMENSION(:), INTENT(INOUT) :: ZCONC
INTEGER, INTENT(IN)               :: KVECNPT
INTEGER :: II
REAL, DIMENSION(KVECNPT) :: RA,XH2O,PVH2O,PVH2SO4,KHI,SIG,XNSULFC,XNSULF
REAL :: Kb,TC,T0


REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('CH_AER_NUCL',0,ZHOOK_HANDLE)
Kb=1.381E-23
TC = 905.15
T0 = 360.15
!   1. Saturation vapor pressure for water (N/m2, T in K)
!         (Preining et al, 1981)
PVH2O(:) = EXP(77.344913-7235.4247/ZT(:)-8.2*LOG(ZT(:))+0.0057113*ZT(:))
!   2.   Water concentration (molec/cm3) 
XH2O(:) = ZRH(:)*PVH2O(:)/(Kb*ZT(:))/1.E6

ZJ(:)=0.
ZAL(:)=0.
WHERE(((ZT(:)>=223.).OR.(ZT(:)<=298)).AND.(ZRH(:)>=0.1))
!   1. Saturation vapor pressure for water (N/m2, T in K)
!         (Preining et al, 1981)
  PVH2O(:) = EXP(77.344913-7235.4247/ZT(:)-8.2*LOG(ZT(:))+0.0057113*ZT(:))
!   2.   Water concentration (molec/cm3) 
  XH2O(:) = ZRH(:)*PVH2O(:)/(Kb*ZT(:))/1.E6
!   3.   Saturation vapor pressure for H2SO4 
!         (Kulmala et al 1990, Seinfeld 577p)
!PVH2SO4 = 1./T-1./T0+0.38*T0*(1./T-1./T0)/(TC-T0)
!PVH2SO4 = PVH2SO4 + 0.38/(TC-T0)*LOG(T/T0)
!PVH2SO4 = -10156.0*PVH2SO4 - 0.414
!PVH2SO4= EXP(PVH2SO4)
  PVH2SO4(:)=EXP(-10156./T0+16.259+10156.*(-1./ZT(:)+1./T0+0.38/(TC-T0)*&
       (1.+LOG(T0/ZT(:))-T0/ZT(:))))*101325.

!   4.   Relative Acidity
  RA(:)=ZCONC(:)*1.E6*(Kb*ZT(:))/PVH2SO4(:)

END WHERE
!   5.    H2SO4 mole fraction in the critical nucleous 
  DO II=1,SIZE(ZCONC,1)
  IF((ZCONC(II)>0.).AND.(XH2O(II)>0.)) THEN
    ZAL(II)=1.2233-0.0154*RA(II)/(RA(II)+ZRH(II))+0.0102*&
         LOG(ZCONC(II))-0.0415*LOG(XH2O(II))+0.0016*ZT(II)
  END IF
  END DO

WHERE(((ZT(:)>=223.).OR.(ZT(:)<=298)).AND.(ZRH(:)>=0.1))
  WHERE(ZAL(:)/=0.)
    !   6.    Sulfuric nucleation rate (molec/cm3/s) 
    XNSULFC(:)=EXP(-14.5125+0.1335*ZT(:)-10.5462*ZRH(:)+1958.4*ZRH(:)/ZT(:))
    SIG(:) = 1.+(ZT(:)-273.15)/273.15
    XNSULF(:)=LOG(ZCONC(:)/XNSULFC(:))
    KHI(:)=25.1289*XNSULF(:)-4890.8*XNSULF(:)/ZT(:)-1743.3/ZT(:)-2.2479*SIG(:)*XNSULF(:)*ZRH(:)+&
         7643.4*ZAL(:)/ZT(:)-1.9712*ZAL(:)*SIG(:)/ZRH(:)
    ZJ(:)=EXP(KHI(:))
  END WHERE
END WHERE
!WHERE(RA(:)>1)
!  ZJ(:)=0.
!END WHERE

IF (LHOOK) CALL DR_HOOK('CH_AER_NUCL',1,ZHOOK_HANDLE)
RETURN
END SUBROUTINE CH_AER_NUCL
