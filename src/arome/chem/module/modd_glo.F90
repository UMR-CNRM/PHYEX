!     ######spl
MODULE MODD_GLO

IMPLICIT NONE

!**********************************************
!PARAMETERS USED IN BOTH PUN AND GRIFFIN'S CODES
!***********************************************

! flag to use binary solution and zsr if = 1 */
! zsrflag = 0 calls unifac and solves implicit
! equation for a.c. water = RH using newt1 */ 
INTEGER, PARAMETER :: zsrflag = 1 
     
! saturationflag = 1 means to 
! use saturation to determine particulate-phase
! concentration when inorganic particle is dry. 
! If = 0, use absorption */
INTEGER, PARAMETER  :: saturationflag = 0

! Criterion for setting intital particle 
! conc, compared to VP in torr or VP 
! rated by PAOM */
REAL, PARAMETER :: VPCrit = 1.d-7

!Critical value for LWC above which much aerosol can dissolve
REAL, PARAMETER :: LWCCRIT = 2.0d0 ![ug/m3] 

!Critical value for Henry's law constant above which much aerosol will dissolve
REAL, PARAMETER :: HLCRIT = 1.d-2

REAL, PARAMETER :: TINY = 1.d-8

!Set the universal gas constant in non SI units
REAL, PARAMETER       :: R_UNIV = 8.206e-5  !gas constant in units of (m3*atm)/(mol*K) 

!Molecular weight of water
REAL, PARAMETER                    :: MW_WATER = 18.0  ![g/mol] molar weight of water

!Are we running box model?
LOGICAL                               :: LBOX=.FALSE.

!Do we want a lot of prints           
LOGICAL                               :: LPRINT=.FALSE.
!************************************************************
!******END PARAMETERS USED IN BOTH CODES 
!************************************************************

!*************************************************************
!*******PARAMETERS USED IN GRIFFIN'S CODE
!************************************************************
INTEGER, PARAMETER     :: NBSPOA = 8 !Number of POA species
INTEGER, PARAMETER     :: NBSP = 10  !Number of main SOA species
INTEGER, PARAMETER     :: NAAERO= 17 !Number of SOA species and their ions

!Molar weight of all parameters including ions
REAL, PARAMETER, DIMENSION(NAAERO) ::  MW_SOA = (/  &
     211.0, 210.0       & !Comp #1 (NK=2)
     ,178.0, 177.0      & !Comp #2 (NK=2)
     ,217.0             & !Comp #3 (NK=1)
     ,303.0             & !Comp #4 (NK=1)
     ,217.0             & !Comp #5 (NK=1)
     ,90.0, 89.0, 88.0  & !Comp #6 (NK=3)
     ,184.0, 183.0, 182.0 &!Comp #7 (NK=3)
     ,154.0               &!Comp #8 (NK=1)
     ,186.0, 185.0        &!Comp #9 (NK=2)
     ,186.0               & !Comp #10 (NK=1)
     /)

!Molecular weight of primary organic paerosols
REAL, PARAMETER, DIMENSION(NBSPOA) :: MW_POA = (/ 408., 118., 216., 276., 412., 166., 284., 390. /)

!Parameter needed to get the saturation vapor pressures of organics
!The below is cut and pasted from aerodriv.f recieved from Griffin
REAL, PARAMETER, DIMENSION(NBSP)  :: HBN = (/6.7e-03,5.61e-03,0.0,3.3e-03,4.61e-03, &
     1.57e-02,7.68e-03,6.49e-03,7.6e-03,5.37e-03/)

REAL, PARAMETER, DIMENSION(NBSP)  :: TAUVP=(/3.5,2.5,6.0,13.5,1.0,0.0,2.5,2.0,3.0,5.0/)

REAL, PARAMETER, DIMENSION(NBSP)  :: TBOIL = (/685.3,634.0,645.5,672.5,566.3,560.0,698.0 &
     ,575.0,679.0,615.0 /) 

!partition parameters H and K
!H is in units of m3/ug estimated based on Suzuki et al., 1992
!K is in units of [mol/kg water] (same as {H+}) with
!concentrations of molecules in ions in the same molar units
REAL, PARAMETER, DIMENSION(NAAERO) :: K_298 = (/ &
     1.82e-2, 1.7e-3  & !Comp 1
     , 1.202e-4, 7.33e-5 & !Comp 2
     , 1.38e-6           & !Comp 3
     , 2.455e-7          & !Comp 4
     , 2.884e-5          & !Comp 5 
     , 2.512e-2, 5.4e-2, 5.2e-5 & !Comp 6,
     , 22.01, 3.7e-5, 3.9e-6    & !Comp 7
     , 1.47e-4                  & !Comp 8
     , 0.0489, 6.52e-4          & !Comp 9 
     , 9.55e-4                  & !Comp 10
     /)

!Total number of (main+sub-components (ions)) per main components
!NK=3 means main component can dissociate twice (acid=H2A)
!NK=2 means main component can dissociate once (acid=HA)
!NK=1 means main component can not dissociate to ions
INTEGER, PARAMETER, DIMENSION(NBSP)  :: NKAQ=(/ 2, 2, 1, 1, 1, 3, 3, 1, 2, 1 /)
!**********************************************
!END PARAMETERS USED ONLY IN GRIFFINS CODE
!**********************************************

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

!*****************************************
!PARAMETERS USED ONLY IN PUN'S CODE
!*****************************************
INTEGER, PARAMETER :: NBSPA=6        !number of molecular species in group A
INTEGER, PARAMETER :: NAAEROA=13     !possible no. of particle phase solutes */ 
INTEGER, PARAMETER  ::NBSPB=5        !number of species in group A 
INTEGER, PARAMETER  ::NBSPAO_PUN=5   !number of species in Primary organics

!*****************TYPE B PARAMETERS *************************************

!Saturation vapor pressure of type B species (torr) 
REAL, DIMENSION(NBSPB) :: VPB = (/ 3.e-10, 3.e-6, 7.e-6, 5.e-8, 3.e-8 /)

!Set molecular weight of Primary aerosol organic matter
REAL, PARAMETER       :: MWAOM = 280.0 ![g/mol]
                      
!Set molecular weight of the secondary organic components
REAL, PARAMETER, DIMENSION(NBSPB)  :: MWB=(/197.0,  164.0, 181.0, 301.0, 170.0 /)

!Assumed molar fraction of the primary organics
REAL, PARAMETER, DIMENSION(NBSPB)  :: XAOM=(/0.4, 0.05, 0.15, 0.12, 0.28 /)
!****************END TYPE B PARAMETERS**************************************

!***************TYPE A PARAMETERS ******************************************

!Molar weight of TYPE A species including ions
REAL, PARAMETER, DIMENSION(NAAEROA) ::  MWA = (/    &
     104.0, 103.0, 102.0                           &
     ,184.0, 183.0, 182.0                          &
     ,154.0                                        &
     ,186.0, 185.0                                 &
     ,186.0                                        &
     ,118.0, 117.0, 116.0                          &
     /)

!partition parameters H and K
!H is in units of (m3/ug) estimated based on Suzuki et al., 1992
!K is in units of [mol/kg water] (same as {H+}) with
!concentrations of molecules in ions in the same mass-based units
!K's of malic acid and glyoxalic acid used respectively for
!compounds that dissociate twice and once 
REAL, PARAMETER, DIMENSION(NAAEROA) :: K_A_298 = (/  & 
     3.87e-5, 3.95e-4, 7.70e-6                      &
     , 22.01, 3.95e-4, 7.70e-6                      &
     , 1.47e-4                                      &
     , 0.0489, 6.52e-4                              &
     , 0.0196                                       &
     , 7.34e-3, 3.95e-4, 7.70e-6                    &
     /)

!Saturation vapor pressure in units of ug/m3
REAL, PARAMETER, DIMENSION(NBSPA)  :: VPA_298=(/ 462.22, 0.0127, 1.86e5, 1.11, 71.62, 50.16 /)

!Same as VP, but this time, unit is torr
REAL, PARAMETER, DIMENSION(NBSPA)  :: VPAtorr_298=(/ 8.26e-5, 1.28e-9, 2.24e-2, 1.1e-7, 7.16e-6, 7.9e-6 /)

!Total number of (main+sub-components (ions)) per main components
!NK=2 means main component can dissociate once
!NK=1 means main component can not dissociate
INTEGER, PARAMETER, DIMENSION(NBSPA)  :: NKA=(/ 3, 3, 1, 2, 1, 3 /)

!Value to determine if component soluble or not
REAL, PARAMETER                       :: CRITSOL=100.d0 ![ug_{aer}/ug_{gas} * kg_{water}/umol_{H+}]
!********************END OF TYPE A PARAMETERS ****************************************************

END MODULE MODD_GLO
