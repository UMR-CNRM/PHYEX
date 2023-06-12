!ORILAM_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!ORILAM_LIC This is part of the ORILAM software governed by the CeCILL-C licence
!ORILAM_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!ORILAM_LIC for details.
!     ######spl
      MODULE MODI_AER_EFFIC
!!    ########################
!!
!
INTERFACE
!!
SUBROUTINE AER_EFFIC(PRG,PVGG,      & !aerosol radius/fall speed (m/s)
                PRHODREF,           & !Air density     
                PMUW, PMU,          & !mu water/air
                PDPG, PEFC,         & !diffusivity, efficiency
                PRRS,               & ! Rain water m.r. at time 
                KMODE,              & ! Number of aerosol modes
                PTEMP, PCOR,        & ! air temp, cunningham corr factor
                PDENSITY_AER,       & ! aerosol density
                PRR, PNT            ) ! radius and number of rain drops
!
IMPLICIT NONE
REAL, DIMENSION(:,:), INTENT(IN) ::  PRG,  PVGG
REAL, DIMENSION(:),   INTENT(IN) ::  PRHODREF
REAL, DIMENSION(:,:), INTENT(IN) ::  PDPG
REAL, DIMENSION(:),   INTENT(IN) ::  PMU, PMUW
REAL, DIMENSION(:,:), INTENT(INOUT) :: PEFC
REAL, DIMENSION(:),   INTENT(IN)    :: PRRS
REAL, DIMENSION(:),   INTENT(IN)    :: PTEMP
REAL, DIMENSION(:,:), INTENT(IN)    :: PCOR
REAL, DIMENSION(:),   INTENT(IN)    :: PRR, PNT
INTEGER, INTENT(IN)                 :: KMODE
REAL, DIMENSION(:,:), INTENT(IN)    :: PDENSITY_AER


END SUBROUTINE  AER_EFFIC 
!!
END INTERFACE
END MODULE MODI_AER_EFFIC
!     ######spl
SUBROUTINE AER_EFFIC(PRG,PVGG,      & !aerosol radius/fall speed (m/s)
                PRHODREF,           & !Air density     
                PMUW, PMU,          & !mu water/air
                PDPG, PEFC,         & !diffusivity, efficiency
                PRRT,               & ! Rain water m.r. at time t 
                KMODE,              & ! Number of aerosol modes
                PTEMP, PCOR,        & ! air temp, cunningham corr factor
                PDENSITY_AER,       & ! aerosol density
                PRR, PNT            ) ! radius and number of rain drops
!!   #######################################
!!**********AER_EFFIC********** 
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
!!   -------------
!!   T. Hoarau (LACy) 15/05/17   add LIMA
!!   Philippe Wautelet 28/05/2018: corrected truncated integer division (1/12 -> 1./12.)
!!   P. Tulet and C. Barthe (LAERO) 15/01/22   correction for lima 
!!
!-----------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
!
USE MODD_RAIN_ICE_PARAM_n,  ONLY : YFSEDR => XFSEDR, YEXSEDR => XEXSEDR
!++cb++
!++th++
USE MODD_RAIN_ICE_DESCR_n,  ONLY : YCCR => XCCR, YLBR => XLBR, YLBEXR => XLBEXR, &
                                 YCEXVT => XCEXVT
USE MODD_PARAM_LIMA_WARM, ONLY : WCCR => XCCR, WLBR => XLBR, WLBEXR => XLBEXR, &
                                 XFSEDRR, XFSEDRC 
USE MODD_PARAM_LIMA,      ONLY : WCEXVT => XCEXVT, WFSEDR => XFSEDR, WFSEDC=>XFSEDC, &
                                 XRTMIN
!--cb--
USE MODD_PARAM_n,        ONLY: CCLOUD
!--th--
USE MODD_CST,            ONLY : XPI, XRHOLW, XP00, XRD
USE MODD_PARAMETERS ,    ONLY : JPVEXT
USE MODD_REF,            ONLY : XTHVREFZ
!
IMPLICIT NONE
!
!*      0.1  declarations of arguments
REAL, DIMENSION(:,:), INTENT(IN) ::  PRG,  PVGG
REAL, DIMENSION(:),   INTENT(IN) ::  PRHODREF
REAL, DIMENSION(:,:), INTENT(IN) ::  PDPG
REAL, DIMENSION(:),   INTENT(IN) ::  PMU, PMUW
REAL, DIMENSION(:,:), INTENT(INOUT) :: PEFC
REAL, DIMENSION(:),   INTENT(IN)    :: PRRT
REAL, DIMENSION(:),   INTENT(IN)    :: PTEMP
REAL, DIMENSION(:),   INTENT(IN)    :: PRR, PNT
REAL, DIMENSION(:,:), INTENT(IN)    :: PCOR
INTEGER, INTENT(IN)                 :: KMODE
REAL, DIMENSION(:,:), INTENT(IN)    :: PDENSITY_AER
!
!
!*      0.2  declaration of local variables
!
INTEGER :: IKB                ! Coordinates of the first physical 
                              ! points along z
REAL :: ZRHO00                ! Surface reference air density
!viscosity ratio, Reynolds number
REAL, DIMENSION(SIZE(PRG,1)) :: ZOMG, ZREY
!rain radius, m, and rain fall speed, m/s; aerosol radius (m),
REAL, DIMENSION(SIZE(PRG,1)) :: ZRR, ZVR 
!lambda, number concentration according to marshall palmer, 
REAL, DIMENSION(SIZE(PRG,1)) :: ZNT, ZLBDA1
!RHO_dref*r_r, Rain LWC
REAL, DIMENSION(SIZE(PRG,1)) :: RLWC
! schmidts number
REAL, DIMENSION(SIZE(PRG,1),KMODE) ::  ZSCH
!
!Stokes number, ratio of diameters,aerosol radius
REAL, DIMENSION(SIZE(PRG,1),KMODE) :: ZSTO, ZPHI, ZRG
! S Star Term
REAL, DIMENSION(SIZE(PRG,1)) :: ZSTA, ZDIFF, ZTAU
!
!Term 1, Term 2, Term 3, Term 4 such that
! E = Term1 * Term 2 + Term 3 + Term 4
REAL, DIMENSION(SIZE(PRG,1),KMODE) :: ZT1, ZT2
REAL, DIMENSION(SIZE(PRG,1),KMODE) :: ZT3, ZT4      
!
INTEGER :: JI,JK
!++th++
REAL :: KLBEXR, KLBR, KCEXVT, KCCR, ZFSEDR, ZBR, ZDR, ZEXSEDR
!--th--
!
!-----------------------------------------------------------------
IKB = 1 + JPVEXT
ZRHO00 = XP00 / (XRD * XTHVREFZ(IKB))                                        
ZRG(:,:) = PRG(:,:) * 1.E-6 !change units to meters
ZVR(:) = 0.

SELECT CASE(CCLOUD)
CASE('ICE3')
  KLBEXR  = YLBEXR
  KLBR    = YLBR
  KCEXVT  = YCEXVT 
  KCCR    = YCCR
  ZFSEDR  = YFSEDR
  ZEXSEDR = YEXSEDR

!Fall Speed calculations
!similar to rain_ice.f90, chapter 17.3.4, MESONH Handbook
  ZVR(:) = ZFSEDR * PRRT(:)**(ZEXSEDR-1) *   &
         PRHODREF(:)**(ZEXSEDR-KCEXVT) 

CASE('LIMA')
  KLBEXR  = WLBEXR
  KLBR    = WLBR
  KCEXVT  = WCEXVT     
  KCCR    = WCCR
  ZFSEDR  = XFSEDRR
  ZBR = 3.0
  ZDR = 0.8
  ZEXSEDR = (ZBR + ZDR + 1.0) / (ZBR + 1.0)
  WHERE (PRRT(:) > XRTMIN(3) .AND. PNT(:) > 0.)
    ZLBDA1(:) = (KLBR * PNT(:) / PRRT(:))**KLBEXR
    ZVR(:) = XFSEDRR * PRHODREF(:)**(1.-KCEXVT) * ZLBDA1(:)**(-ZDR)
  END WHERE
END SELECT


!Fall speed cannot be faster than 7 m/s
ZVR(:) = MIN(ZVR(:), 7.)   

KCCR = 8.E6


!Ref SEINFELD AND PANDIS p.1019
! Viscosity Ratio      
ZOMG(:) = PMUW(:) / PMU(:)
!!Reynolds number
ZREY(:) = PRR(:) * ZVR(:) * PRHODREF(:) / PMU(:)
ZREY(:) =  MAX(ZREY(:), 1.E-2)
!
!S Star
ZSTA(:) = (1.2 + 1./12. * LOG(1.+ZREY(:))) / (1. + LOG(1.+ZREY(:)))

PEFC(:,:) = 0.0
!
DO JI = 1, KMODE
!Scmidts number
  ZSCH(:,JI) = PMU(:) / PRHODREF(:) / PDPG(:,JI)      
!
! Rain-Aerosol relative velocity 
  ZDIFF(:) = MAX(ZVR(:)-PVGG(:,JI), 0.)
!
! Relaxation time 
  ZTAU(:) = (ZRG(:,JI)*2.)**2. * PDENSITY_AER(:,JI) * PCOR(:,JI) / (18. * PMU(:))
!
! Stockes number
  ZSTO(:,JI) = ZTAU(:) * ZDIFF(:) / PRR(:)
!
!Ratio of diameters  
  ZPHI(:,JI) = ZRG(:,JI) / PRR(:)
  ZPHI(:,JI) = MIN(ZPHI(:,JI), 1.)
!
!Term 1
  ZT1(:,JI) = 4.0 / ZREY(:) / ZSCH(:,JI)
!
!Term 2
  ZT2(:,JI) = 1.0 + 0.4  * ZREY(:)**(0.5) * ZSCH(:,JI)**(1./3.) + &
                    0.16 * ZREY(:)**(0.5) * ZSCH(:,JI)**(0.5)   
!     
!Brownian diffusion          
  ZT1(:,JI) = ZT1(:,JI) * ZT2(:,JI)
!
!Term 3 - interception
  ZT3(:,JI) = 4. * ZPHI(:,JI) * (1. / ZOMG(:) + &
             (1.0 + 2.0 * ZREY(:)**0.5) * ZPHI(:,JI))
!
  ZT4(:,JI) = 0.0     
!      
  WHERE(ZSTO(:,JI) .GT. ZSTA(:)) 
!Term 4 - impaction
    ZT4(:,JI) = ((ZSTO(:,JI) - ZSTA(:)) /                     &
                 (ZSTO(:,JI) - ZSTA(:) + 2. / 3.))**(3./2.) * &
                 (XRHOLW / PDENSITY_AER(:,JI))**(1./2.)
              
  END WHERE
!
!Collision Efficiancy 
  PEFC(:,JI) = ZT1(:,JI) + ZT3(:,JI) + ZT4(:,JI)     
!
! Physical radius of a rain collector droplet up than 20 um
  WHERE (PRR(:) .LE. 20.E-6)
    PEFC(:,JI) = 0.
  END WHERE
ENDDO
!
PEFC(:,:) = MIN(PEFC(:,:), 1.0)
PEFC(:,:) = MAX(PEFC(:,:), 0.0)

END SUBROUTINE AER_EFFIC
