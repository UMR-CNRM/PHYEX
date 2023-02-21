!MNH_LIC Copyright 2013-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!######################################
 MODULE MODI_INIT_AEROSOL_CONCENTRATION
!######################################
!
INTERFACE INIT_AEROSOL_CONCENTRATION
   SUBROUTINE INIT_AEROSOL_CONCENTRATION(PRHODREF, PSVT, PZZ)
!
     REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF !Air Density [kg/m**3]
     REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PSVT     !Particles Concentration [/m**3]
     REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PZZ      ! Height (z)
!
   END SUBROUTINE INIT_AEROSOL_CONCENTRATION
END INTERFACE INIT_AEROSOL_CONCENTRATION
!
END MODULE MODI_INIT_AEROSOL_CONCENTRATION
!
!     ##########################################################
      SUBROUTINE INIT_AEROSOL_CONCENTRATION(PRHODREF, PSVT, PZZ)
!     ##########################################################
!!
!!    PURPOSE
!!    -------
!!    Define the aerosol distributions
!! 
!!
!!      MODD_BLANKn :
!!      CDUMMY2 : CCN ou IFN pour le panache
!!      NDUMMY1 : hauteur base du panache
!!      NDUMMY2 : hauteur sommet du panache
!!      XDUMMY8 : Concentration du panache (N/cm3 pour des CCN, N/L pour des IFN)
!!
!!
!!    AUTHOR
!!    ------
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!      S.    Berthet    * Laboratoire d'Aerologie*
!!      B.    Vié        * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original             ??/??/13 
!!      Modification    01/2016  (JP Pinty) Add LIMA
!!
!!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_NSV
USE MODD_PARAM_n,    ONLY : CCLOUD
USE MODD_PARAM_LIMA, ONLY : NMOM_C, LACTI, NMOD_CCN, LSCAV, LAERO_MASS,     &
                            XCCN_CONC, LCCN_HOM,                            &
                            NMOM_I, LNUCL, NMOD_IFN, LMEYERS,               &
                            XIFN_CONC, LIFN_HOM
USE MODD_PARAMETERS, ONLY : JPVEXT
USE MODD_BLANK_n,      ONLY : CDUMMY2, NDUMMY1, NDUMMY2, XDUMMY8
!
IMPLICIT NONE
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF !Air Density [kg/m**3]
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PSVT     !Particles Concentration    
                                                    ![particles/kg of dry air]
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PZZ      ! Height (z)
!
! Local variables
INTEGER                                 :: JMOD_IFN
INTEGER                                 :: JSV, JINIT
INTEGER                                 :: IKB, IKE
!
!-------------------------------------------------------------------------------
!
!
!*initialization of N_FREE_CCN/N_ACTIVATED_CCN et N_FREE_IN/N_ACTIVATED_IN
!
!
IF ( NMOM_C.GE.2 .AND. LACTI ) THEN
  DO JSV = NSV_LIMA_CCN_FREE, NSV_LIMA_CCN_ACTI+NMOD_CCN-1
    PSVT(:,:,:,JSV) = 0.0
  ENDDO
  IKB = 1+JPVEXT
  IKE = SIZE(PSVT,3)-JPVEXT
!
! Initialisation des concentrations en CCN
!
!
  IF (LCCN_HOM) THEN
! concentration homogène (en nombre par m3) sur la verticale
    DO JSV = 1, NMOD_CCN
      PSVT(:,:,IKB:IKE,NSV_LIMA_CCN_FREE+JSV-1) = &
        XCCN_CONC(JSV)*1.0E6 / PRHODREF(:,:,IKB:IKE)
    END DO
  ELSE
! concentration décroissante selon z
    DO JSV = 1, NMOD_CCN
      WHERE (PZZ(:,:,:) .LE. 1000.) 
         PSVT(:,:,:,NSV_LIMA_CCN_FREE+JSV-1) = XCCN_CONC(JSV)*1.0E6 / PRHODREF(:,:,:)
      ELSEWHERE (PZZ(:,:,:) .LE. 10000.)
         PSVT(:,:,:,NSV_LIMA_CCN_FREE+JSV-1) = XCCN_CONC(JSV)*1.0E6 &
                 / PRHODREF(:,:,:) * EXP(-LOG(XCCN_CONC(JSV)/0.01)*PZZ(:,:,:)/10000.)
      ELSEWHERE
         PSVT(:,:,:,NSV_LIMA_CCN_FREE+JSV-1) = 0.01*1.0E6 / PRHODREF(:,:,:)
      ENDWHERE
    END DO
  ENDIF
END IF ! LWARM AND LACTI
!
! Initialisation des concentrations en IFN
!
IF ( NMOM_I.GE.2 .AND. LNUCL .AND. (.NOT. LMEYERS) ) THEN
  DO JSV = NSV_LIMA_IFN_FREE, NSV_LIMA_IFN_NUCL+NMOD_IFN-1
    PSVT(:,:,:,JSV) = 0.0
  ENDDO
  IKB = 1+JPVEXT
  IKE = SIZE(PSVT,3)-JPVEXT
!
  IF (LIFN_HOM) THEN
! concentration homogène (en nombre par m3) sur la verticale
    DO JSV = 1, NMOD_IFN
      PSVT(:,:,IKB:IKE,NSV_LIMA_IFN_FREE+JSV-1) = &
             XIFN_CONC(JSV)*1.0E3 / PRHODREF(:,:,IKB:IKE)
     END DO
  ELSE
! concentration décroissante selon z
    DO JSV = 1, NMOD_IFN
      WHERE (PZZ(:,:,:) .LE. 1000.)
         PSVT(:,:,:,NSV_LIMA_IFN_FREE+JSV-1) = XIFN_CONC(JSV)*1.0E3 / PRHODREF(:,:,:)
      ELSEWHERE (PZZ(:,:,:) .LE. 10000.)
         PSVT(:,:,:,NSV_LIMA_IFN_FREE+JSV-1) = XIFN_CONC(JSV)*1.0E3 &
                 / PRHODREF(:,:,:) * EXP(-LOG(XIFN_CONC(JSV)/1.)*PZZ(:,:,:)/10000.)
      ELSEWHERE
         PSVT(:,:,:,NSV_LIMA_IFN_FREE+JSV-1) = 1*1.0E3 / PRHODREF(:,:,:)
      ENDWHERE
    END DO
  ENDIF
END IF ! LCOLD AND LNUCL AND NOT LMEYERS
!
!
! Cas d'un panache de "pollution", concentration homogène dans le panache :
!
SELECT CASE (CDUMMY2)
  CASE ('CCN')
    PSVT(:,:,:,NSV_LIMA_CCN_FREE+NMOD_CCN-1)=0.
    WHERE ( (PZZ(:,:,:) .GE. NDUMMY1) .AND. (PZZ(:,:,:) .LE. NDUMMY2) ) &
      PSVT(:,:,:,NSV_LIMA_CCN_FREE+NMOD_CCN-1)=XDUMMY8*1.0E6 / PRHODREF(:,:,:)
  CASE ('IFN')
    PSVT(:,:,:,NSV_LIMA_IFN_FREE+NMOD_IFN-1)=0.
    WHERE ( (PZZ(:,:,:) .GE. NDUMMY1) .AND. (PZZ(:,:,:) .LE. NDUMMY2) ) &
      PSVT(:,:,:,NSV_LIMA_IFN_FREE+NMOD_IFN-1)=XDUMMY8*1.0E3 / PRHODREF(:,:,:)
END SELECT
!
! 
END SUBROUTINE INIT_AEROSOL_CONCENTRATION
