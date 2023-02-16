!ORILAM_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!ORILAM_LIC This is part of the ORILAM software governed by the CeCILL-C licence
!ORILAM_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!ORILAM_LIC for details.
!
!     ######spll
      MODULE MODI_AER_EFFIC3D
!!    ########################
!!
!
INTERFACE
!!
SUBROUTINE AER_EFFIC3D(PRG,PVGG,    & !aerosol radius/fall speed (m/s)
                PRHODREF,           & !Air density     
                PMUW, PMU,          & !mu water/air
                PDPG,               & !diffusivity
                PURR,               & ! Rain water m.r. at time t 
                NMODE_DST,          & ! Number of aerosol modes
                PTEMP, PCOR,        & ! air temp, cunningham corr factor
                PDENSITY_AER,       & ! aerosol density
                PEFFIC )              ! scavenging efficiency
!
IMPLICIT NONE
REAL, DIMENSION(:,:,:,:), INTENT(IN) ::  PRG,  PVGG
REAL, DIMENSION(:,:,:),   INTENT(IN) ::  PRHODREF
REAL, DIMENSION(:,:,:,:), INTENT(IN) ::  PDPG
REAL, DIMENSION(:,:,:),   INTENT(IN) ::  PMU, PMUW
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PURR
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PTEMP
REAL, DIMENSION(:,:,:,:), INTENT(IN) :: PCOR
INTEGER, INTENT(IN)                  :: NMODE_DST
REAL, DIMENSION(:,:,:,:), INTENT(IN) :: PDENSITY_AER
REAL, DIMENSION(:,:,:,:), INTENT(INOUT):: PEFFIC



END SUBROUTINE  AER_EFFIC3D 
!!
END INTERFACE
END MODULE MODI_AER_EFFIC3D
!     ######spll
SUBROUTINE AER_EFFIC3D(PRG,PVGG,    & !aerosol radius/fall speed (m/s)
                PRHODREF,           & !Air density     
                PMUW, PMU,          & !mu water/air
                PDPG,               & !diffusivity
                PURR,               & ! Rain water m.r. at time t 
                NMODE_DST,          & ! Number of aerosol modes
                PTEMP, PCOR,        & ! air temp, cunningham corr factor
                PDENSITY_AER,       & ! aerosol density
                PEFFIC )              ! scavenging efficiency
!!   #######################################
!!**********AER_EFFIC3D********** 
!!   PURPOSE
!!   -------
!!  Calculate the collection efficiency of
!   a falling drop interacting with a dust aerosol
!   for use with aer_wet_dep_kmt_warm.f90
!!
!!**  METHOD
!!    ------
!!    Using basic theory, and the one dimensional variables sent 
!!    from aer_wet_dep_kmt_warm.f90, calculation of the average 
!!    fall speed calculations, chapter 17.3.4, MESONH Handbook
!!    droplet number based on the Marshall_Palmer distribution
!!    and Stokes number, Reynolds number, etc. based on theory
!!    (S&P, p.1019)
!!
!!   REFERENCE
!!   ---------
!!   Seinfeld and Pandis p.1019
!!   MESONH Handbook chapter 17.3.4
!!
!!   AUTHOR
!!    ------
!!   K. Crahan Kaku / P. Tulet (CNRM/GMEI)
!!
!!   MODIFICATIONS
!!    -------------
!!  Philippe Wautelet 28/05/2018: corrected truncated integer division (1/12 -> 1./12.)
!!
!-----------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
!
USE MODD_RAIN_ICE_PARAM
USE MODD_RAIN_ICE_DESCR
USE MODD_CST,         ONLY : XPI, XRHOLW, XP00, XRD
USE MODD_PARAMETERS , ONLY : JPVEXT
USE MODD_REF,         ONLY : XTHVREFZ
!
IMPLICIT NONE
!
!*      0.1  declarations of arguments
!
REAL, DIMENSION(:,:,:,:), INTENT(IN) ::  PRG,  PVGG
REAL, DIMENSION(:,:,:),   INTENT(IN) ::  PRHODREF
REAL, DIMENSION(:,:,:,:), INTENT(IN) ::  PDPG
REAL, DIMENSION(:,:,:),   INTENT(IN) ::  PMU, PMUW
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PURR
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PTEMP
REAL, DIMENSION(:,:,:,:), INTENT(IN) :: PCOR
INTEGER, INTENT(IN)                  :: NMODE_DST
REAL, DIMENSION(:,:,:,:), INTENT(IN) :: PDENSITY_AER
REAL, DIMENSION(:,:,:,:), INTENT(INOUT):: PEFFIC
!
!*      0.2  declaration of local variables
!
INTEGER :: IKB                ! Coordinates of the first physical 
                              ! points along z
REAL :: ZRHO00                ! Surface reference air density
!viscosity ratio, Reynolds number
REAL, DIMENSION(SIZE(PRG,1),SIZE(PRG,2),SIZE(PRG,3)) :: ZOMG, ZREY
!rain radius, m, and rain fall speed, m/s; aerosol radius (m),
REAL, DIMENSION(SIZE(PRG,1),SIZE(PRG,2),SIZE(PRG,3)) :: ZRR, ZVR 
!lambda, number concentration according to marshall palmer, 
REAL, DIMENSION(SIZE(PRG,1),SIZE(PRG,2),SIZE(PRG,3)) :: ZNT, ZLBDA
! Rain water m.r. source
REAL, DIMENSION(SIZE(PRG,1),SIZE(PRG,2),SIZE(PRG,3)) :: ZRRS
!RHO_dref*r_r, Rain LWC
REAL, DIMENSION(SIZE(PRG,1),SIZE(PRG,2),SIZE(PRG,3)) :: ZRLWC
! schmidts number
REAL, DIMENSION(SIZE(PRG,1),SIZE(PRG,2),SIZE(PRG,3),NMODE_DST) ::  ZSCH
!
!Stokes number, ratio of diameters,aerosol radius
REAL, DIMENSION(SIZE(PRG,1),SIZE(PRG,2),SIZE(PRG,3),NMODE_DST) :: ZSTO, ZPHI, ZRG
! S Star Term
REAL, DIMENSION(SIZE(PRG,1),SIZE(PRG,2),SIZE(PRG,3)) :: ZSTA, ZDIFF, ZTAU
!
!Term 1, Term 2, Term 3, Term 4 such that
! E = Term1 * Term 2 + Term 3 + Term 4
REAL, DIMENSION(SIZE(PRG,1),SIZE(PRG,2),SIZE(PRG,3),NMODE_DST) :: ZT1, ZT2
REAL, DIMENSION(SIZE(PRG,1),SIZE(PRG,2),SIZE(PRG,3),NMODE_DST) :: ZT3, ZT4      
!
INTEGER :: JI,JK
!
!-----------------------------------------------------------------
ZLBDA = 1E20
ZNT = 1E-20
ZRR = 10E-6
ZRRS(:,:,:)=PURR(:,:,:)
IKB = 1 + JPVEXT
ZRHO00 = XP00/(XRD*XTHVREFZ(IKB))                                        
ZRG(:,:,:,:)=PRG(:,:,:,:)*1.E-6 !change units to meters
!
!Fall Speed calculations
!similar to rain_ice.f90, chapter 17.3.4, MESONH Handbook
!
ZVR (:,:,:)=  XFSEDR  * ZRRS(:,:,:)**(XEXSEDR-1) *   &
          PRHODREF(:,:,:)**(XEXSEDR-XCEXVT-1) 

! Drop Radius calculation in m 
!lbda = pi*No*rho(lwc)/(rho(dref)*rain rate) p.212 MESONH Handbook
! compute the slope parameter Lbda_r

WHERE((ZRRS(:,:,:).GT. 0.).AND.(PRHODREF(:,:,:) .GT. 0.))

ZLBDA(:,:,:)  = XLBR*(PRHODREF(:,:,:)*ZRRS(:,:,:))**XLBEXR
!Number concentration NT=No/lbda   p. 415 Jacobson
ZNT(:,:,:) = XCCR/ZLBDA(:,:,:)  
!rain  lwc (kg/m3) =  rain m.r.(kg/kg) * rho_air(kg/m3)
ZRLWC(:,:,:)=ZRRS(:,:,:)*PRHODREF(:,:,:)
!4/3 *pi *r³*NT*rho_eau(kg/m3) =rho(lwc)=rho(air)* qc(kg/kg)
ZRR(:,:,:) =  (ZRLWC(:,:,:)/(XRHOLW*ZNT(:,:,:)*4./3.*XPI))**(1./3.)
END WHERE

ZRR(:,:,:) = MIN(ZRR(:,:,:), 100.E-6)
!Fall speed cannot be faster than 7 m/s
ZVR (:,:,:)=MIN(ZVR (:,:,:),7.)   

!Ref SEINFELD AND PANDIS p.1019
! Viscosity Ratio      
ZOMG(:,:,:)=PMUW(:,:,:)/PMU(:,:,:)
!!Reynolds number
ZREY(:,:,:)=ZRR(:,:,:)*ZVR(:,:,:)*PRHODREF(:,:,:)/PMU(:,:,:)
ZREY(:,:,:)= MAX(ZREY(:,:,:), 1E-2)


!S Star
ZSTA(:,:,:)=(1.2+(1./12.)*LOG(1.+ZREY(:,:,:)))/(1.+LOG(1.+ZREY(:,:,:)))
PEFFIC(:,:,:,:)=0.0
DO JI=1,NMODE_DST
!
!Scmidts number
  ZSCH(:,:,:,JI)=PMU(:,:,:)/PRHODREF(:,:,:)/PDPG(:,:,:,JI)      
! Rain-Aerosol relative velocity 
  ZDIFF(:,:,:) = MAX(ZVR(:,:,:)-PVGG(:,:,:,JI),0.)
! Relaxation time 
  ZTAU(:,:,:) = (ZRG(:,:,:,JI)*2.)**2. * PDENSITY_AER(:,:,:,JI) * PCOR(:,:,:,JI) / (18.*PMU(:,:,:))
! Stockes number
  ZSTO(:,:,:,JI)= ZTAU(:,:,:) * ZDIFF(:,:,:) / ZRR(:,:,:)
!Ratio of diameters  
  ZPHI(:,:,:,JI)=ZRG(:,:,:,JI)/ZRR(:,:,:)
  ZPHI(:,:,:,JI)=MIN(ZPHI(:,:,:,JI), 1.)
!Term 1
  ZT1(:,:,:,JI)=4.0/ZREY(:,:,:)/ZSCH(:,:,:,JI)
!Term 2
  ZT2(:,:,:,JI)=1.0+(0.4*ZREY(:,:,:)**(0.5)*ZSCH(:,:,:,JI)**(1./3.))+ &
            (0.16*ZREY(:,:,:)**(0.5)*ZSCH(:,:,:,JI)**(0.5))   
     
!Brownian diffusion          
  ZT1(:,:,:,JI)= ZT1(:,:,:,JI)*ZT2(:,:,:,JI)
!Term 3 - interception
  ZT3(:,:,:,JI)=4.*ZPHI(:,:,:,JI)*(1./ZOMG(:,:,:)+        &
             (1.0+(2.0*ZREY(:,:,:)**0.5))*ZPHI(:,:,:,JI))

  ZT4(:,:,:,JI)=0.0           
 WHERE(ZSTO(:,:,:,JI).GT.ZSTA(:,:,:)) 
!Term 4 - impaction
   ZT4(:,:,:,JI)=((ZSTO(:,:,:,JI)-ZSTA(:,:,:))/            &
              (ZSTO(:,:,:,JI)-ZSTA(:,:,:)+2./3.))**(3./2.) &
             *((XRHOLW/PDENSITY_AER(:,:,:,JI))**(1./2.))

  END WHERE
!Collision Efficiancy 
   PEFFIC(:,:,:,JI)=ZT1(:,:,:,JI)+  ZT3(:,:,:,JI)+ZT4(:,:,:,JI)     
! Physical radius of a rain collector droplet up than 20 um
WHERE (ZRR(:,:,:) .LE. 9.9E-6)
  PEFFIC(:,:,:,JI)= 0.
END WHERE
ENDDO
PEFFIC(:,:,:,:)=MIN(PEFFIC(:,:,:,:),1.0)
PEFFIC(:,:,:,:)=MAX(PEFFIC(:,:,:,:),0.0)

END SUBROUTINE AER_EFFIC3D
