!     ######spl
     SUBROUTINE CH_AER_PUN(PTOT, PTOTG,PTEMP, PRH, PLWC, PPROTON, PTOTNEW, PTOTGNEW)
     USE PARKIND1, ONLY : JPRB
     USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!###########################################################################################
!!
!!   PURPOSE
!!   -------
!!   solve the organic thermodynamic balance , if use CACM our ReLACS2 chemical
!!   scheme
!!
!!    REFERENCE
!!    ---------
!!    none
!!
!!    AUTHOR
!!    ------
!!    A. Grini / P. Tulet  (GMEI)
!!
!!    MODIFICATIONS
!!    -------------
!!
!!    EXTERNAL
!!    --------
!!    None
!!
USE MODD_CH_AEROSOL
!!
!! ORGANIC AEROSOL MODULE
USE MODD_GLO
USE MODE_OAMAIN

IMPLICIT NONE
!!
!!  Declaration arguments
REAL, DIMENSION(:,:),   INTENT(IN)    :: PTOT, PTOTG
REAL, DIMENSION(:),     INTENT(IN)    :: PTEMP, PRH
REAL, DIMENSION(:),     INTENT(INOUT) :: PLWC, PPROTON
REAL, DIMENSION(:,:),   INTENT(OUT)   :: PTOTNEW, PTOTGNEW


!!
!!  Declaration variables internes
!! 
INTEGER :: IFLAG=1  !Flag to do iterations or not for type B aerosols
REAL, DIMENSION(SIZE(PTOT,1), NBSPA+NBSPB)        :: ZAER
REAL, DIMENSION(SIZE(PTOT,1),NBSPA+NBSPB)         :: ZGAS   
REAL, DIMENSION(SIZE(PTOT,1),NBSPA+NBSPB)         :: ZWORG    ![ug/m3] total (gas+aerosol) SOA
REAL, DIMENSION(SIZE(PTOT,1))                     :: ZPROTON  ![mole/kg_{water}] proton concentration
!The only use of these variables which would make sense, would be
!to have them as intent OUT and send them back to the inorganic module
REAL, DIMENSION(SIZE(PTOT,1))                     :: ZORGANION
REAL, DIMENSION(SIZE(PTOT,1))                     :: ZDELTALWC
REAL, DIMENSION(SIZE(PTOT,1))                     :: ZPAOM

!Fraction of primary aerosol assumed to be butadionic acid
REAL, PARAMETER                                   :: ZPOA_TO_BUTACID=0.01d0 !1%


!  For all variables:
!  First dimension is the data of simulation domain
!  Second dimension represent the aerosol species
!  1 => NSP : mineral aerosol  species
!  NSP+1 => NSP + NCARB : primary aerosol species
!  NSP + NCARB +1 => NSP + NCARB + NSOA : secondary organic aerosol
!                                      Third dimension (only for PCTOTA and PCCTOT): lognormal mode

!*****************************************************************
!*****************************************************************
! SOLVEUR DE L'EQUILIBRE CHIMIQUE PUN
!*****************************************************************
!*****************************************************************
! Option  flag to turn on type B oa module if == 1
!          (type a is iterated, but no need to rerun type b)

!PPROTON(:) = MAX(PPROTON(:), 0.)
!PLWC(:)    = MAX(PLWC(:), 0.)
!ZLWC(:)   = PLWC(:)
!ZRH(:)    = PRH(:)
!ZTEMP(:)  = PTEMP(:)

!Convert proton concentration to right units
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('CH_AER_PUN',0,ZHOOK_HANDLE)
ZPROTON(:)= PPROTON(:)*1.d3  ![mol/g_{water}] ==> [mole/kg_{water}]
!Initialize (don't think it is needed)
!ZORGANION(:) = 0.
!ZDELTALWC(:) = 0.

!Concentration of water soluble species (1-6)
ZAER(:,1)=PTOT(:,JP_AER_SOA1)
ZAER(:,2)=PTOT(:,JP_AER_SOA2)
ZAER(:,3)=PTOT(:,JP_AER_SOA3)
ZAER(:,4)=PTOT(:,JP_AER_SOA4)
ZAER(:,5)=PTOT(:,JP_AER_SOA5)
!assume butadionic acid is a fraction of primary
ZAER(:,6)=ZPOA_TO_BUTACID*PTOT(:,JP_AER_OC)
!Concentration of organic species (7-11)
ZAER(:,7)=PTOT(:,JP_AER_SOA6)
ZAER(:,8)=PTOT(:,JP_AER_SOA7)
ZAER(:,9)=PTOT(:,JP_AER_SOA8)
ZAER(:,10)=PTOT(:,JP_AER_SOA9)
ZAER(:,11)=PTOT(:,JP_AER_SOA10)

ZGAS(:,1)=PTOTG(:,JP_AER_SOA1)
ZGAS(:,2)=PTOTG(:,JP_AER_SOA2)
ZGAS(:,3)=PTOTG(:,JP_AER_SOA3)
ZGAS(:,4)=PTOTG(:,JP_AER_SOA4)
ZGAS(:,5)=PTOTG(:,JP_AER_SOA5)
ZGAS(:,6)=0.

ZGAS(:,7)=PTOTG(:,JP_AER_SOA6)
ZGAS(:,8)=PTOTG(:,JP_AER_SOA7)
ZGAS(:,9)=PTOTG(:,JP_AER_SOA8)
ZGAS(:,10)=PTOTG(:,JP_AER_SOA9)
ZGAS(:,11)=PTOTG(:,JP_AER_SOA10)

!Get total organic aerosol
ZWORG(:,:) = ZGAS(:,:) + ZAER(:,:)


!ZWORG(:,11)   = 0. ! non volatil organic coupounds for absorbing medium (primary aerosol as butanedioic acid)
!Get primary organic aerosol
ZPAOM(:)=PTOT(:,JP_AER_OC)*(1.d0-ZPOA_TO_BUTACID)

CALL PUN(              &
     PTEMP             & !I [K] Temperature
     ,PRH              & !I [-] Relative humidity
     ,ZWORG            & !I [ug/m3] total (aerosol+gas) conc
     ,ZGAS             & !O [ug/m3] gas phase concentrations
     ,ZPAOM            & !I [ug/m3] Primary aerosol organic matter
     ,ZAER             & !O [ug/m3] aerosol phase concentrations
     ,PLWC             & !I [ug/m3] liquid water content
     ,ZPROTON          & !I [mol/kg_{water}] H+ concentrations
     ,ZDELTALWC        & !O [ug/m3] change in LWC
     ,ZORGANION        & !O [mole//m3] anion charge
     ,IFLAG            & !I [1/0] Flag for iteration or not
     )

!aerosol concentration of hydrophilic
PTOTNEW(:,JP_AER_SOA1)=ZAER(:,1)
PTOTNEW(:,JP_AER_SOA2)=ZAER(:,2)
PTOTNEW(:,JP_AER_SOA3)=ZAER(:,3)
PTOTNEW(:,JP_AER_SOA4)=ZAER(:,4)
PTOTNEW(:,JP_AER_SOA5)=ZAER(:,5)

!aerosol concentration of hydrophobic
PTOTNEW(:,JP_AER_SOA6)=ZAER(:,7)
PTOTNEW(:,JP_AER_SOA7)=ZAER(:,8)
PTOTNEW(:,JP_AER_SOA8)=ZAER(:,9)
PTOTNEW(:,JP_AER_SOA9)=ZAER(:,10)
PTOTNEW(:,JP_AER_SOA10)=ZAER(:,11)

!Gasesous concentration of hydrophilic
PTOTGNEW(:,JP_AER_SOA1)=ZGAS(:,1)
PTOTGNEW(:,JP_AER_SOA2)=ZGAS(:,2)
PTOTGNEW(:,JP_AER_SOA3)=ZGAS(:,3)
PTOTGNEW(:,JP_AER_SOA4)=ZGAS(:,4)
PTOTGNEW(:,JP_AER_SOA5)=ZGAS(:,5)

!Gasesous concentration of hydrophobic
PTOTGNEW(:,JP_AER_SOA6)=ZGAS(:,7)
PTOTGNEW(:,JP_AER_SOA7)=ZGAS(:,8)
PTOTGNEW(:,JP_AER_SOA8)=ZGAS(:,9)
PTOTGNEW(:,JP_AER_SOA9)=ZGAS(:,10)
PTOTGNEW(:,JP_AER_SOA10)=ZGAS(:,11)

!New water content
PTOTNEW(:,JP_AER_H2O)=PTOT(:,JP_AER_H2O)+ZDELTALWC(:)

!Convert back ZPROTON to right units
!mole/kg ==> mole/g
PPROTON(:)    = ZPROTON(:)*1.d-3

IF (LHOOK) CALL DR_HOOK('CH_AER_PUN',1,ZHOOK_HANDLE)
END SUBROUTINE CH_AER_PUN
