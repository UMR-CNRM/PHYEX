!     ######spl
      SUBROUTINE INI_AEROSOLS_CAMS
      USE YOMHOOK, ONLY : LHOOK, DR_HOOK, JPHOOK
!     ###########################################################
!
!!****  *INI_AEROSOLS_CAMS* - initialize the constants necessary 
!!                            for cams aerosols
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to initialize the constants used to
!!    resolve the mixed phase microphysical scheme. 
!!      AEROSOLS FROM CAMS: A total number of 14 species.
!!
!!**  METHOD
!!    ------
!!      The constants are initialized to their numerical values and the number
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_AEROSOL_PROP
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    D. Martin-Perez 11/02/2019
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_AEROSOL_PROP
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
!
!*       0.2   Declarations of local variables :
!
!
LOGICAL  :: GFLAG   ! Logical flag for printing the constatnts on the output
                    ! listing
!-------------------------------------------------------------------------------
!
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('INI_AEROSOLS_CAMS',0,ZHOOK_HANDLE)
!
!
!-------------------------------------------------------------------------------
!
!*       1.     CLOUD CONDENSATION NUMBER
!               -------------------------
!
!     XACCNI: Index of aerosol species considered condensation nuclei


XACCNI(:)=0
XIFNI(:)=0
XCOARSE(:)=0
XSEASALT(:)=0
XDUST(:)=0
XDRYDEP(:)=0
XINTSS(:)=0

! Cloud condensation nuclei species
NCCN=9 ! Number of species becoming CCN
XACCNI(1)=1   ! Sea salt 1
XACCNI(2)=2   ! Sea salt 2
XACCNI(3)=3   ! Sea salt 3
XACCNI(4)=8   ! Hydrophilic Organic Matter
XACCNI(5)=10  ! Hydrophilic Black Carbon
XACCNI(6)=11  ! Sulfate
XACCNI(7)=12  ! Nitrate fine mode
XACCNI(8)=13  ! Nitrate coarse mode
XACCNI(9)=14  ! Ammonium

! Ice nuclei species
NIFN=5 ! Number of species becoming IFN
XIFNI(1)=4  ! Desert dust 1
XIFNI(2)=5  ! Desert dust 2
XIFNI(3)=6  ! Desert dust 3
XIFNI(4)=7  ! Hydrophobic OM
XIFNI(5)=9  ! Hydrophobic BC

! For supercoarse sea salt and dust (Remy 2019)
NCOAR=2
XCOARSE(1)=3
XCOARSE(2)=6

! Sea salt
NSEASALT=3
XSEASALT(1)=1
XSEASALT(2)=2
XSEASALT(3)=3
! Coarse Sea salt
NCOARSEASALT=2
XCOARSEASALT(1)=2
XCOARSEASALT(2)=3

! Desert Dust
NDUST=3
XDUST(1)=4
XDUST(2)=5
XDUST(3)=6

! Dry Deposition activated
NDRYDEP=6
XDRYDEP(1)=1
XDRYDEP(2)=2
XDRYDEP(3)=3
XDRYDEP(4)=4
XDRYDEP(5)=5
XDRYDEP(6)=6


!
!*       2.     AEROSOL PHYSICAL PROPERTIES
!               ---------------------------
!

!
! Optical properties
!
!     XMEXT: Mass Extinction (m2 kg-1)
!
! Aerosol distribution: log normal size distribution 
!
!     XMRAE: Modal Radius(micrometers),Rg
!     XSDAE: geometric Standard Ddeviation,sigma 
!
! Limits of the error function
!
!     XBINDOWN: Lower bin limit (micrometers)
!     XBINUP: Upper bin limit (micrometers)
!
!     XERFDOWN(JC)=0.5*erf(log({DOWNlimit}/Rg)/(sqrt(2.0))*log(sigma)
!     XERFUP(JC)=0.5*erf(log({UPlimit}/Rg)/(sqrt(2.0))*log(sigma)
!
! Mean cubic radius of the distribution:
!    <r3> = Rg**3*((3.0/sqrt(2.0)))*log(sigma))**2
!    X    = (3.0/sqrt(2.0)))*log(sigma)
!    {ERF3DOWN,ERF3UP} =
!    0.5*erf(log({DOWNlimit,Uplimit}-X)/Rg)/(sqrt(2.0))*log(sigma)
!    Correction of the mean cubic radius when the bin size limits are considered
!
!    XR3CORR(JC)= <r3>*(ERF3UP-ERF3DOWN)/(XERFUP(JC)-XERFDOWN(JC)
!
! Physical properties
!
!     XRHOAE: mass density (Kg m-3)
!     XKHYGROS: Hygroscopic parameter
!
!     XMEXT: Mass extinction 
!
! Variables for removal processes
!    
!     XVDRY: Dry Deposition Velocity (cm s-1) (Remy et al., 2019)
!     XFRICS: Fraction of aerosols in aqueous phase for in-cloud scavenging
!             (Remy et al., 2019)

XMEXT(:,:)=0.0
XMRAE(:)=0.0
XSDAE(:)=0.0
XBINDOWN(:)=0.0
XBINUP(:)=0.0
XERFDOWN(:)=0.0
XERFUP(:)=0.0
XR3CORR(:)=0.0
XRHOAE(:)=0.0
XKHYGROS(:)=0.0
XVDRY(:,:)=0.0 ! XVDRY= Dry deposition over (sea,land)
XFRICS(:)=0.0

!-------------------------------------------------
! 2.1  Sea salt 
!-------------------------------------------------
!
!      Sea salt aerosol density (Kg m-3)
!          2160 (dry particles)
!          1182 (sea salt production asuming 80% humidity, Mocrette et al. 2009)

! Factor for Sea Salt emision
! Monahan, 1986
! INTEGRAL (1+0.057*r**1.05)10**(1.19*exp(-B**2))
! B=(0.380-log(r))/0.650
! r: aerosol radius in um

XINTSS(1)=2.09913  ! r: 2x0.03 - 2x0.5 um
XINTSS(2)=29.6649  ! r: 2x0.5  - 2x5.0 um
XINTSS(3)=80.3799  ! r: 2x5.0  - 2x20.0 um

!
! Sea salt, size bin limits: 0.03 - 0.5 microm.
! Size distribution parameters:
XMRAE(1)=0.1992
XSDAE(1)= 1.92
XBINDOWN(1)=0.03
XBINUP(1)=0.5
XERFDOWN(1)=ERFPDF(XBINDOWN(1),XMRAE(1),XSDAE(1)) ! -0.4981
XERFUP(1)=ERFPDF(XBINUP(1),XMRAE(1),XSDAE(1))     !  0.4208
! Mean cubic radius of the distribution:  <r3>=5.3640e-02
! Corrected mean cubic radius
XR3CORR(1)=1.7071e-02
! Density (kg/m3):
XRHOAE(1)=  2160.0

XKHYGROS(1)=1.28

XVDRY(1,:)=(/0.05,0.2/)
XFRICS(1)=0.7
 
! Mass extinction
XMEXT(1,:)=(/827.98,755.72,684.17,614.08,2197.81,2612.42,3223.86,4048.34,5381.25,6506.83,8497.96,13330.20/)

!
! Sea salt, size bin limits: 0.5 - 5.0 microm.
!

! Size distribution parameters:
XMRAE(2)=1.9920
XSDAE(2)= 2.00
XBINDOWN(2)=0.5
XBINUP(2)=5.0
XERFDOWN(2)=ERFPDF(XBINDOWN(2),XMRAE(2),XSDAE(2)) ! -0.4769
XERFUP(2)=ERFPDF(XBINUP(2),XMRAE(2),XSDAE(2))     !  0.4079
! Mean cubic radius of the distribution:    <r3>=6.8680e+01
! Corrected mean cubic radius
XR3CORR(2)=1.7549e+01
! Density (kg/m3):
XRHOAE(2)=2160.0

XKHYGROS(2)=1.28

XVDRY(2,:)=(/0.05,0.5/)
XFRICS(2)=0.9

! Mass extinction
XMEXT(2,:)=(/143.79,143.36,143.30,143.72,285.42,330.74,376.77,432.98,520.54,590.13,713.12,1030.59/)

!
! Sea salt, size bin limits: 5.0 - 20.0 microm.
!
! Size distribution parameters:
XMRAE(3)=1.9920
XSDAE(3)=2.00
XBINDOWN(3)=5.0
XBINUP(3)=20.0
XERFDOWN(3)=ERFPDF(XBINDOWN(3),XMRAE(3),XSDAE(3)) ! 0.4079
XERFUP(3)=ERFPDF(XBINUP(3),XMRAE(3),XSDAE(3))     ! 0.4996
! Mean cubic radius of the distribution:   <r3>=6.8680e+01
! Corrected mean cubic radius
XR3CORR(3)=5.0026e+02
! Density (kg/m3):
XRHOAE(3)=2160.0

XKHYGROS(3)=1.28

XVDRY(3,:)=(/1.20,1.50/)
XFRICS(3)=0.9

! Mass extinction
XMEXT(3,:)=(/38.65,38.85,38.85,38.63,79.27,92.17,105.23,122.77,149.16,171.24,209.44,309.58/)

!-------------------------------------------------
! 2.2 Desert Dust
!-------------------------------------------------

!
! Desert Dust, size bin limits: 0.03 - 0.55 microm. 
!
! Size distribution parameters:
XMRAE(4)=0.2900
XSDAE(4)=2.00
XBINDOWN(4)=0.03
XBINUP(4)=0.55
XERFDOWN(4)=ERFPDF(XBINDOWN(4),XMRAE(4),XSDAE(4)) ! -0.4995
XERFUP(4)=ERFPDF(XBINUP(4),XMRAE(4),XSDAE(4))     ! 0.3221
! Mean cubic radius of the distribution:   <r3>=2.1191e-01
! Corrected mean cubic radius
XR3CORR(4)=3.1940e-02
! Density (kg/m3):
XRHOAE(4)=2610.0

XKHYGROS(4)=0.0

XVDRY(4,:)=(/0.02,0.02/)
XFRICS(4)=0.3

! Mass extinction
XMEXT(4,:)=2496.68

!
! Desert Dust, size bin limits: 0.55 - 0.9 microm. 
!

! Size distribution parameters:
XMRAE(5)=0.2900
XSDAE(5)=2.00
XBINDOWN(5)=0.55
XBINUP(5)=0.9
XERFDOWN(5)=ERFPDF(XBINDOWN(5),XMRAE(5),XSDAE(5)) ! 0.3221
XERFUP(5)=ERFPDF(XBINUP(5),XMRAE(5),XSDAE(5))     ! 0.4489
! Mean cubic radius of the distribution:     <r3>=2.1191e-01
! Corrected mean cubic radius
XR3CORR(5)=3.4124e-01
! Density (kg/m3):
XRHOAE(5)=2610.0

XKHYGROS(5)=0.0

XVDRY(5,:)=(/0.2,0.2/)
XFRICS(5)=0.4

! Mass extinction
XMEXT(5,:)=955.078

!
! Desert Dust, size bin limits: 0.9 - 20.0 microm. 
!

! Size distribution parameters:
XMRAE(6)=0.2900
XSDAE(6)=2.00
XBINDOWN(6)=0.9
XBINUP(6)=20.0
XERFDOWN(6)=ERFPDF(XBINDOWN(6),XMRAE(6),XSDAE(6)) ! 0.4489
XERFUP(6)=ERFPDF(XBINUP(6),XMRAE(6),XSDAE(6))     ! 0.5000
! Mean cubic radius of the distribution:    <r3>=2.1191e-01
! Corrected mean cubic radius
XR3CORR(6)=2.7845e+00
! Density (kg/m3):
XRHOAE(6)=2610.0

XKHYGROS(6)=0.0

XVDRY(6,:)=(/1.2,1.2/)
XFRICS(6)=0.4

! Mass extinction
XMEXT(6,:)=406.533

!-------------------------------------------------
! 2.3  Organic matter
!-------------------------------------------------
!

!
! 7:Hydrophobic Organic Matter
!

! Size distribution parameters:
XMRAE(7)=0.02
XSDAE(7)=2.2
XBINDOWN(7)=0.005
XBINUP(7)=20.0
XERFDOWN(7)=ERFPDF(XBINDOWN(7),XMRAE(7),XSDAE(7)) ! -0.4996
XERFUP(7)=ERFPDF(XBINUP(7),XMRAE(7),XSDAE(7))     ! 0.4996
! Mean cubic radius of the distribution:   <r3>=1.0861e-03
! Corrected mean cubic radius
XR3CORR(7)=9.7073e-04
! Density (kg/m3):
XRHOAE(7)=1825.0

XKHYGROS(7)=0.0

XVDRY(7,:)=(/0.1,0.1/)
XFRICS(7)=0.0

! Mass extinction (RH=80)
XMEXT(7,:)=2321.03

!
! 8:Hydrophilic Organic Matter
!
! Size distribution parameters:
XMRAE(8)=0.02
XSDAE(8)=2.2
XBINDOWN(8)=0.005
XBINUP(8)=20.0
XERFDOWN(8)=ERFPDF(XBINDOWN(8),XMRAE(8),XSDAE(8)) ! -0.4996
XERFUP(8)=ERFPDF(XBINUP(8),XMRAE(8),XSDAE(8))     ! 0.4996
! Mean cubic radius of the distribution:   <r3>=1.0861e-03
! Corrected mean cubic radius
XR3CORR(8)=9.7073e-04
! Density (kg/m3):
XRHOAE(8)=1825.0

XKHYGROS(8)=0.3

XVDRY(8,:)=(/0.1,0.1/)
XFRICS(8)=0.7

! Mass extinction
XMEXT(8,:)=(/1942.13,2131.58,2321.03,2510.48,2699.93,2889.38,3185.61,3481.84,4119.88,4909.06,5698.24,8206.75/)

!
!-------------------------------------------------
! 2.4  Black Carbon
!-------------------------------------------------
!

! Hydrophobic Black Carbon size bin limits  0.005 -  0.500 micrometers.

! Size distribution (Oshima et al., 2009)
!    XMRAE(9)=0.053; XSDAE(9)=1.5; XRHOAE(9)=1800.
! Size distribution parameters:
XMRAE(9)=0.0118
XSDAE(9)=2.00
XBINDOWN(9)=0.005
XBINUP(9)=0.5
XERFDOWN(9)=ERFPDF(XBINDOWN(9),XMRAE(9),XSDAE(9)) ! -0.3923
XERFUP(9)=ERFPDF(XBINUP(9),XMRAE(9),XSDAE(9))     !  0.5000
! Mean cubic radius of the distribution:    <r3>=1.4276e-05
! Corrected mean cubic radius
XR3CORR(9)=1.5985e-05
! Density (kg/m3):
XRHOAE(9)=1000.0

XKHYGROS(9)=0.0

XVDRY(9,:)=(/0.1,0.1/)
XFRICS(9)=0.0

XMEXT(9,:)=13487.8

!
! Hydrophilic Black Carbon size bin limits  0.005 -  0.500 micrometers.
!


! Size distribution parameters:
XMRAE(10)=0.0118
XSDAE(10)=2.00
XBINDOWN(10)=0.005
XBINUP(10)=0.5
XERFDOWN(10)=ERFPDF(XBINDOWN(10),XMRAE(10),XSDAE(10)) ! -0.3923
XERFUP(10)=ERFPDF(XBINUP(10),XMRAE(10),XSDAE(10))     !  0.5000
! Mean cubic radius of the distribution:    <r3>=1.4276e-05
! Corrected mean cubic radius
XR3CORR(10)=1.5985e-05
! Density (kg/m3):
XRHOAE(10)=1000.0

XKHYGROS(10)=0.1

XVDRY(10,:)=(/0.1,0.1/)
XFRICS(10)=0.7

! Mass extinction
XMEXT(10,:)=13487.8

!
!-------------------------------------------------
! 2.5. Sulphate properties (0.005-20 micrometers)
!-------------------------------------------------
!
! Sulfate (Vd=0.05 for oceans; Vd=0.25 all other surfaces )

! Size distribution parameters:
XMRAE(11)=0.0355
XSDAE(11)=2.00
XBINDOWN(11)=0.005
XBINUP(11)=20.0
XERFDOWN(11)=ERFPDF(XBINDOWN(11),XMRAE(11),XSDAE(11)) ! -0.4977
XERFUP(11)=ERFPDF(XBINUP(11),XMRAE(11),XSDAE(11))     !  0.5000
! Mean cubic radius of the distribution:   <r3>=3.8873e-04
! Corrected mean cubic radius
XR3CORR(11)=3.8964e-04
! Density (kg/m3):
XRHOAE(11)=1760.0

XKHYGROS(11)=0.6

XVDRY(11,:)=(/0.15,0.25/)
XFRICS(11)=0.7

! Mass extinction
XMEXT(11,:)=(/3119.62,3510.56,3901.51,4292.45,4683.40,5074.35,5685.64,6296.94,7613.58,9242.11,10870.64,16047.13/)

!
!-------------------------------------------------
! 2.6. Nitrate properties (0.005-0.9 micrometers)
!-------------------------------------------------

!
! Nitrate fine mode (NH4NO3), size bin limits: 0.005 - 0.9 microm.
!

! Size distribution parameters:
XMRAE(12)=0.0355
XSDAE(12)=2.0
XBINDOWN(12)=0.005
XBINUP(12)=0.9
XERFDOWN(12)=ERFPDF(XBINDOWN(12),XMRAE(12),XSDAE(12)) ! -0.4977
XERFUP(12)=ERFPDF(XBINUP(12),XMRAE(12),XSDAE(12))     !  0.5000
! Corrected mean cubic radius
XR3CORR(12)=3.8774e-04
! Density (kg/m3):
XRHOAE(12)=1730.0

XKHYGROS(12)=0.64

XVDRY(12,:)=(/0.15,0.15/)
XFRICS(12)=0.4

XMEXT(12,:)=(/3501.201,3501.201,3501.201,3501.201,4757.47,5354.868,6161.467,7361.848,13555.21,15649.4,14113.68,23958.8/)

!
! Nitrate coarse mode ( NaNO3+Ca(NO3)2 ), size bin limits: 0.9 - 20.0 microm.
!

! Size distribution parameters:
XMRAE(13)=1.992
XSDAE(13)=2.0
XBINDOWN(13)=0.9
XBINUP(13)=20.0
XERFDOWN(13)=ERFPDF(XBINDOWN(13),XMRAE(13),XSDAE(13)) !  -0.3741
XERFUP(13)=ERFPDF(XBINUP(13),XMRAE(13),XSDAE(13))     !  0.4996
! Corrected mean cubic radius
XR3CORR(13)=7.0228e+01
! Density (kg/m3):
XRHOAE(13)=1400.0

XKHYGROS(13)=0.9

XVDRY(13,:)=(/0.15,0.15/)
XFRICS(13)=0.4

! Mass extinction
XMEXT(13,:)=(/4295.583,4295.583,4295.583,4295.583,5150.058,6185.985,6780.04,7425.708,8129.261,10575.19,14709.95,26338.47/)

!
!-------------------------------------------------
! 2.7. Ammonium properties (0.005-20 micrometers)
!-------------------------------------------------
!      (Ammonium sulfate for Amonia)
! Ammonium (Vd=0.05 for oceans; Vd=0.25 all other surfaces )

! Size distribution parameters:
XMRAE(14)=0.0355
XSDAE(14)=2.00
XBINDOWN(14)=0.005
XBINUP(14)=20.0
XERFDOWN(14)=ERFPDF(XBINDOWN(14),XMRAE(14),XSDAE(14)) ! -0.4977
XERFUP(14)=ERFPDF(XBINUP(14),XMRAE(14),XSDAE(14))     ! 0.5000
! Corrected mean cubic radius
XR3CORR(14)=3.8964e-04
! Density (kg/m3):
XRHOAE(14)=1760.0

XKHYGROS(14)=0.6

XVDRY(14,:)=(/0.15,0.15/)
XFRICS(14)=0.4

! Mass extinction
XMEXT(14,:)=(/190.3313,232.0307,274.9777,321.6228,346.3498,371.914,425.8234,483.4848,544.9835,610.332,751.631,906.7416/)


!
IF (LHOOK) CALL DR_HOOK('INI_AEROSOL_CAMS',1,ZHOOK_HANDLE)

CONTAINS
!
!----------------------------------------
!

FUNCTION ERFPDF(PRAD,PMODRAD,PSIGMA) RESULT (PERFPDF)
!
! 0.5*erf(log(PRAD/Rg)/(sqrt(2.0)*log(sigma)))
!
!
  IMPLICIT NONE
!
  REAL   :: PRAD
  REAL   :: PMODRAD
  REAL   :: PSIGMA
  REAL   :: PERFPDF
!
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  IF (LHOOK) CALL DR_HOOK('INI_AEROSOLS_CAMS:ERFPDF',0,ZHOOK_HANDLE)
!
  PERFPDF = 0.5*erf(0.7071*LOG(PRAD/PMODRAD)/LOG(PSIGMA))
!
  IF (LHOOK) CALL DR_HOOK('INI_AEROSOLS_CAMS:ERFPDF',1,ZHOOK_HANDLE)
!
END FUNCTION ERFPDF


FUNCTION RAD3CORR(PMODRAD,PSIGMA,PBINDOWN,PBINUP) RESULT (PR3CORR)
! Mean cubic radius of the distribution:
!    <r3> = Rg**3*((3.0/sqrt(2.0)))*log(sigma))**2
!    X    = (3.0/sqrt(2.0)))*log(sigma)
!    {ERF3DOWN,ERF3UP} =
!    0.5*erf(log({DOWNlimit,Uplimit}-X)/Rg)/(sqrt(2.0))*log(sigma)
!    Correction of the mean cubic radius when the bin size limits are considered
!
!    XR3CORR(JC)= <r3>*(ERF3UP-ERF3DOWN)/(XERFUP(JC)-XERFDOWN(JC)
!
  IMPLICIT NONE
!
  REAL   :: PMODRAD,PSIGMA,PBINDOWN,PBINUP
  REAL   :: PR3CORR
  REAL   :: ZR3, ZX,ZWUP,ZWDOWN,ZW3UP,ZW3DOWN
!
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  IF (LHOOK) CALL DR_HOOK('INI_AEROSOLS_CAMS:RAD3CORR',0,ZHOOK_HANDLE)

  ZX=2.12132*LOG(PSIGMA)

  ZR3=(PMODRAD**3)*ZX**2

  ZWUP=ERFPDF(PBINUP,PMODRAD,PSIGMA)
  ZWDOWN=ERFPDF(PBINDOWN,PMODRAD,PSIGMA)
  ZW3UP=ERFPDF(PBINUP-ZX,PMODRAD,PSIGMA)
  ZW3DOWN=ERFPDF(PBINDOWN-ZX,PMODRAD,PSIGMA)

  PR3CORR=ZR3*(ZW3UP-ZW3DOWN)/(ZWUP-ZWDOWN)

  IF (LHOOK) CALL DR_HOOK('INI_AEROSOLS_CAMS:RAD3CORR',1,ZHOOK_HANDLE)

END FUNCTION RAD3CORR

END SUBROUTINE INI_AEROSOLS_CAMS
