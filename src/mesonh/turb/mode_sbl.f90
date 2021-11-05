!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 mode 2006/05/18 13:07:25
!-----------------------------------------------------------------
!     ############### 
      MODULE MODE_SBL 
!     ###############
!
!!****  *MODE_SBL * - contains Surface Boundary Layer characteristics functions
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
!!    
!!    
!!
!!    EXTERNAL
!!    --------
!!       
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!
!!    REFERENCE
!!    ---------
!!
!!    Businger et al 1971,   Wyngaard and Cote 1974
!!      
!!
!!    AUTHOR
!!    ------
!!	V. Masson       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original        13/10/99
!!      V. Masson       06/11/02 optimization and add Businger fonction for TKE
!!      V. Masson       01/01/03 use PAULSON_PSIM function
!-----------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
!
INTERFACE BUSINGER_PHIM
  MODULE PROCEDURE BUSINGER_PHIM_0D
  MODULE PROCEDURE BUSINGER_PHIM_1D
  MODULE PROCEDURE BUSINGER_PHIM_2D
  MODULE PROCEDURE BUSINGER_PHIM_3D
END INTERFACE
INTERFACE BUSINGER_PHIH
  MODULE PROCEDURE BUSINGER_PHIH_0D
  MODULE PROCEDURE BUSINGER_PHIH_1D
  MODULE PROCEDURE BUSINGER_PHIH_2D
  MODULE PROCEDURE BUSINGER_PHIH_3D
END INTERFACE
INTERFACE BUSINGER_PHIE
  MODULE PROCEDURE BUSINGER_PHIE_3D
END INTERFACE
INTERFACE PAULSON_PSIM
  MODULE PROCEDURE PAULSON_PSIM_0D
  MODULE PROCEDURE PAULSON_PSIM_1D
  MODULE PROCEDURE PAULSON_PSIM_2D
END INTERFACE
INTERFACE LMO
  MODULE PROCEDURE LMO_0D
  MODULE PROCEDURE LMO_1D
  MODULE PROCEDURE LMO_2D
END INTERFACE
INTERFACE USTAR
  MODULE PROCEDURE USTAR_0D
  MODULE PROCEDURE USTAR_1D
  MODULE PROCEDURE USTAR_2D
END INTERFACE
!
!-------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------
!
FUNCTION BUSINGER_PHIM_3D(PZ_O_LMO)
  REAL, DIMENSION(:,:,:), INTENT(IN)                 :: PZ_O_LMO
  REAL, DIMENSION(SIZE(PZ_O_LMO,1), &
                  SIZE(PZ_O_LMO,2),SIZE(PZ_O_LMO,3)) :: BUSINGER_PHIM_3D
!
  WHERE ( PZ_O_LMO(:,:,:) < 0. )
    BUSINGER_PHIM_3D(:,:,:) = (1.-15.*PZ_O_LMO)**(-0.25)
  ELSEWHERE
    BUSINGER_PHIM_3D(:,:,:) = 1. + 4.7 * PZ_O_LMO
  END WHERE
END FUNCTION BUSINGER_PHIM_3D
!
!-------------------------------------------------------------------------------
!
FUNCTION BUSINGER_PHIM_2D(PZ_O_LMO)
  REAL, DIMENSION(:,:), INTENT(IN)                   :: PZ_O_LMO
  REAL, DIMENSION(SIZE(PZ_O_LMO,1),SIZE(PZ_O_LMO,2)) :: BUSINGER_PHIM_2D
!
  WHERE ( PZ_O_LMO(:,:) < 0. )
    BUSINGER_PHIM_2D(:,:) = (1.-15.*PZ_O_LMO)**(-0.25)
  ELSEWHERE
    BUSINGER_PHIM_2D(:,:) = 1. + 4.7 * PZ_O_LMO
  END WHERE
END FUNCTION BUSINGER_PHIM_2D
!
!-------------------------------------------------------------------------------
!
FUNCTION BUSINGER_PHIM_1D(PZ_O_LMO)
  REAL, DIMENSION(:), INTENT(IN)  :: PZ_O_LMO
  REAL, DIMENSION(SIZE(PZ_O_LMO)) :: BUSINGER_PHIM_1D
!
  WHERE ( PZ_O_LMO(:) < 0. )
    BUSINGER_PHIM_1D(:) = (1.-15.*PZ_O_LMO)**(-0.25)
  ELSEWHERE
    BUSINGER_PHIM_1D(:) = 1. + 4.7 * PZ_O_LMO
  END WHERE
END FUNCTION BUSINGER_PHIM_1D
!
!-------------------------------------------------------------------------------
!
FUNCTION BUSINGER_PHIM_0D(PZ_O_LMO)
  REAL, INTENT(IN)                   :: PZ_O_LMO
  REAL                               :: BUSINGER_PHIM_0D
!
  IF ( PZ_O_LMO < 0. ) THEN
    BUSINGER_PHIM_0D = (1.-15.*PZ_O_LMO)**(-0.25)
  ELSE
    BUSINGER_PHIM_0D = 1. + 4.7 * PZ_O_LMO
  END IF
END FUNCTION BUSINGER_PHIM_0D
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
FUNCTION BUSINGER_PHIH_3D(PZ_O_LMO)
  REAL, DIMENSION(:,:,:), INTENT(IN)                 :: PZ_O_LMO
  REAL, DIMENSION(SIZE(PZ_O_LMO,1), &
                  SIZE(PZ_O_LMO,2),SIZE(PZ_O_LMO,3)) :: BUSINGER_PHIH_3D
!
  WHERE ( PZ_O_LMO(:,:,:) < 0. )
    BUSINGER_PHIH_3D(:,:,:) = 0.74 * (1.-9.*PZ_O_LMO)**(-0.5)
  ELSEWHERE
    BUSINGER_PHIH_3D(:,:,:) = 0.74 + 4.7 * PZ_O_LMO
  END WHERE
END FUNCTION BUSINGER_PHIH_3D
!
!-------------------------------------------------------------------------------
!
FUNCTION BUSINGER_PHIH_2D(PZ_O_LMO)
  REAL, DIMENSION(:,:), INTENT(IN)                   :: PZ_O_LMO
  REAL, DIMENSION(SIZE(PZ_O_LMO,1),SIZE(PZ_O_LMO,2)) :: BUSINGER_PHIH_2D
!
  WHERE ( PZ_O_LMO(:,:) < 0. )
    BUSINGER_PHIH_2D(:,:) = 0.74 * (1.-9.*PZ_O_LMO)**(-0.5)
  ELSEWHERE
    BUSINGER_PHIH_2D(:,:) = 0.74 + 4.7 * PZ_O_LMO
  END WHERE
END FUNCTION BUSINGER_PHIH_2D
!
!-------------------------------------------------------------------------------
!
FUNCTION BUSINGER_PHIH_1D(PZ_O_LMO)
  REAL, DIMENSION(:), INTENT(IN)  :: PZ_O_LMO
  REAL, DIMENSION(SIZE(PZ_O_LMO)) :: BUSINGER_PHIH_1D
!
  WHERE ( PZ_O_LMO(:) < 0. )
    BUSINGER_PHIH_1D(:) = 0.74 * (1.-9.*PZ_O_LMO)**(-0.5)
  ELSEWHERE
    BUSINGER_PHIH_1D(:) = 0.74 + 4.7 * PZ_O_LMO
  END WHERE
END FUNCTION BUSINGER_PHIH_1D
!
!-------------------------------------------------------------------------------
!
FUNCTION BUSINGER_PHIH_0D(PZ_O_LMO)
  REAL, INTENT(IN)                   :: PZ_O_LMO
  REAL                               :: BUSINGER_PHIH_0D
!
  IF ( PZ_O_LMO < 0. ) THEN
    BUSINGER_PHIH_0D = 0.74 * (1.-9.*PZ_O_LMO)**(-0.5)
  ELSE
    BUSINGER_PHIH_0D = 0.74 + 4.7 * PZ_O_LMO
  END IF
END FUNCTION BUSINGER_PHIH_0D
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
FUNCTION BUSINGER_PHIE_3D(PZ_O_LMO)
  USE MODD_CTURB
  REAL, DIMENSION(:,:,:), INTENT(IN)                 :: PZ_O_LMO
  REAL, DIMENSION(SIZE(PZ_O_LMO,1), &
                  SIZE(PZ_O_LMO,2),SIZE(PZ_O_LMO,3)) :: BUSINGER_PHIE_3D
!
  WHERE ( PZ_O_LMO(:,:,:) < 0. )
    BUSINGER_PHIE_3D(:,:,:) =   (1.+(-PZ_O_LMO)**(2./3.)/XALPSBL) &
                              * (1.-15.*PZ_O_LMO)**(0.5)
  ELSEWHERE
    BUSINGER_PHIE_3D(:,:,:) = 1./(1. + 4.7 * PZ_O_LMO)**2
  END WHERE
END FUNCTION BUSINGER_PHIE_3D
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
FUNCTION PAULSON_PSIM_2D(PZ_O_LMO)
  USE MODD_CST
  REAL, DIMENSION(:,:), INTENT(IN)                   :: PZ_O_LMO
  REAL, DIMENSION(SIZE(PZ_O_LMO,1),SIZE(PZ_O_LMO,2)) :: PAULSON_PSIM_2D
!
  REAL, DIMENSION(SIZE(PZ_O_LMO,1),SIZE(PZ_O_LMO,2)) :: ZX

  ZX=1.
  WHERE ( PZ_O_LMO(:,:) < 0. )
    ZX=(1.-15.*PZ_O_LMO)**(0.25)
    PAULSON_PSIM_2D(:,:) = LOG( (1.+ZX**2)*(1+ZX)**2/8. ) - 2.*ATAN(ZX) + XPI/2.
  ELSEWHERE
    PAULSON_PSIM_2D(:,:) = - 4.7 * PZ_O_LMO
  END WHERE
END FUNCTION PAULSON_PSIM_2D
!
!-------------------------------------------------------------------------------
!
FUNCTION PAULSON_PSIM_1D(PZ_O_LMO)
  USE MODD_CST
  REAL, DIMENSION(:), INTENT(IN)    :: PZ_O_LMO
  REAL, DIMENSION(SIZE(PZ_O_LMO,1)) :: PAULSON_PSIM_1D
!
  REAL, DIMENSION(SIZE(PZ_O_LMO,1)) :: ZX

  ZX=1.
  WHERE ( PZ_O_LMO(:) < 0. )
    ZX=(1.-15.*PZ_O_LMO)**(0.25)
    PAULSON_PSIM_1D(:) = LOG( (1.+ZX**2)*(1+ZX)**2/8. ) - 2.*ATAN(ZX) + XPI/2.
  ELSEWHERE
    PAULSON_PSIM_1D(:) = - 4.7 * PZ_O_LMO
  END WHERE
END FUNCTION PAULSON_PSIM_1D
!
!-------------------------------------------------------------------------------
!
FUNCTION PAULSON_PSIM_0D(PZ_O_LMO)
  USE MODD_CST
  REAL, INTENT(IN)    :: PZ_O_LMO
  REAL                :: PAULSON_PSIM_0D
!
  REAL                :: ZX

  ZX=1.
  IF ( PZ_O_LMO < 0. ) THEN
    ZX=(1.-15.*PZ_O_LMO)**(0.25)
    PAULSON_PSIM_0D = LOG( (1.+ZX**2)*(1+ZX)**2/8. ) - 2.*ATAN(ZX) + XPI/2.
  ELSE
    PAULSON_PSIM_0D = - 4.7 * PZ_O_LMO
  END IF
END FUNCTION PAULSON_PSIM_0D
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
FUNCTION LMO_2D(PUSTAR,PTHETA,PRV,PSFTH,PSFRV)
  USE MODD_CST
  USE MODD_PARAMETERS
  REAL, DIMENSION(:,:), INTENT(IN)               :: PUSTAR
  REAL, DIMENSION(:,:), INTENT(IN)               :: PTHETA
  REAL, DIMENSION(:,:), INTENT(IN)               :: PRV
  REAL, DIMENSION(:,:), INTENT(IN)               :: PSFTH
  REAL, DIMENSION(:,:), INTENT(IN)               :: PSFRV
  REAL, DIMENSION(SIZE(PUSTAR,1),SIZE(PUSTAR,2)) :: LMO_2D
!
  REAL, DIMENSION(SIZE(PUSTAR,1),SIZE(PUSTAR,2)) :: ZTHETAV
  REAL, DIMENSION(SIZE(PUSTAR,1),SIZE(PUSTAR,2)) :: ZQ0
  REAL                                           :: ZEPS
!
!
  ZEPS=(XRV-XRD)/XRD
  ZTHETAV(:,:) = PTHETA(:,:) * ( 1. +ZEPS * PRV(:,:))
  ZQ0    (:,:) = PSFTH(:,:) + ZTHETAV(:,:) * ZEPS * PSFRV(:,:) 
!
  LMO_2D(:,:) = XUNDEF
  WHERE ( ZQ0(:,:) /=0.  )                                   &
    LMO_2D(:,:) = - MAX(PUSTAR(:,:),1.E-6)**3                &
                  / ( XKARMAN * XG / ZTHETAV(:,:) *ZQ0(:,:) )

END FUNCTION LMO_2D
!
!-------------------------------------------------------------------------------
!
FUNCTION LMO_1D(PUSTAR,PTHETA,PRV,PSFTH,PSFRV)
  USE MODD_CST
  USE MODD_PARAMETERS
  REAL, DIMENSION(:), INTENT(IN)  :: PUSTAR
  REAL, DIMENSION(:), INTENT(IN)  :: PTHETA
  REAL, DIMENSION(:), INTENT(IN)  :: PRV
  REAL, DIMENSION(:), INTENT(IN)  :: PSFTH
  REAL, DIMENSION(:), INTENT(IN)  :: PSFRV
  REAL, DIMENSION(SIZE(PUSTAR))   :: LMO_1D
!
  REAL, DIMENSION(SIZE(PUSTAR))   :: ZTHETAV
  REAL                                           :: ZEPS
!
!
  ZEPS=(XRV-XRD)/XRD
!
  ZTHETAV(:) = PTHETA(:) * ( 1. +ZEPS * PRV(:))
!
  LMO_1D(:) = XUNDEF
  WHERE ( PSFTH(:)/ZTHETAV(:)+ZEPS*PSFRV(:)/=0. )                  &
    LMO_1D(:) = - MAX(PUSTAR(:),1.E-6)**3                          &
              / ( XKARMAN * (  XG / ZTHETAV(:)    * PSFTH(:)       &
                             + XG * ZEPS * PSFRV(:) )  )
END FUNCTION LMO_1D
!
!-------------------------------------------------------------------------------
!
FUNCTION LMO_0D(PUSTAR,PTHETA,PRV,PSFTH,PSFRV)
  USE MODD_CST
  USE MODD_PARAMETERS
  REAL, INTENT(IN)  :: PUSTAR
  REAL, INTENT(IN)  :: PTHETA
  REAL, INTENT(IN)  :: PRV
  REAL, INTENT(IN)  :: PSFTH
  REAL, INTENT(IN)  :: PSFRV
  REAL              :: LMO_0D
!
  REAL              :: ZTHETAV
  REAL              :: ZEPS
!
!
  ZEPS=(XRV-XRD)/XRD
!
!
  ZTHETAV = PTHETA * ( 1. +ZEPS * PRV)
!
  LMO_0D = XUNDEF
  IF ( PSFTH/ZTHETAV+ZEPS*PSFRV/=0. )                      &
  LMO_0D = - MAX(PUSTAR,1.E-6)**3                          &
           / ( XKARMAN * (  XG / ZTHETAV       * PSFTH     &
                          + XG * ZEPS * PSFRV )  )
END FUNCTION LMO_0D
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
FUNCTION USTAR_2D(PU,PV,PZ,PZ0,PLMO)
  USE MODD_CST
  USE MODD_PARAMETERS
  REAL, DIMENSION(:,:), INTENT(IN)               :: PU
  REAL, DIMENSION(:,:), INTENT(IN)               :: PV
  REAL, DIMENSION(:,:), INTENT(IN)               :: PZ
  REAL, DIMENSION(:,:), INTENT(IN)               :: PZ0
  REAL, DIMENSION(:,:), INTENT(IN)               :: PLMO
  REAL, DIMENSION(SIZE(PU,1),SIZE(PU,2))         :: USTAR_2D

  REAL, DIMENSION(SIZE(PU,1),SIZE(PU,2))         :: ZZ_O_LMO
  REAL, DIMENSION(SIZE(PU,1),SIZE(PU,2))         :: ZZ0_O_LMO
!
!* purely unstable case
  USTAR_2D(:,:) = 0.
  ZZ_O_LMO(:,:) = XUNDEF
  ZZ0_O_LMO(:,:) = XUNDEF
!
!* general case
  WHERE(ABS(PLMO) > 1.E-20 .AND. PLMO/=XUNDEF)
    ZZ_O_LMO  = PZ(:,:)  / PLMO(:,:)
    ZZ0_O_LMO = PZ0(:,:) / PLMO(:,:)
    USTAR_2D(:,:) = SQRT( PU(:,:)**2+PV(:,:)**2 )               &
                  * XKARMAN / ( LOG(PZ(:,:)/PZ0(:,:))           &
                                 - PAULSON_PSIM(ZZ_O_LMO(:,:))  &
                                 + PAULSON_PSIM(ZZ0_O_LMO(:,:)) )
  END WHERE
!
!* purely neutral case
  WHERE(PLMO==XUNDEF)
    ZZ_O_LMO = 0.
    USTAR_2D(:,:) = SQRT( PU(:,:)**2+PV(:,:)**2 )      &
                  * XKARMAN / LOG(PZ(:,:)/PZ0(:,:))
  END WHERE
!
END FUNCTION USTAR_2D
!
!-------------------------------------------------------------------------------
!
FUNCTION USTAR_1D(PU,PV,PZ,PZ0,PLMO)
  USE MODD_CST
  USE MODD_PARAMETERS
  REAL, DIMENSION(:), INTENT(IN)               :: PU
  REAL, DIMENSION(:), INTENT(IN)               :: PV
  REAL, DIMENSION(:), INTENT(IN)               :: PZ
  REAL, DIMENSION(:), INTENT(IN)               :: PZ0
  REAL, DIMENSION(:), INTENT(IN)               :: PLMO
  REAL, DIMENSION(SIZE(PU))                    :: USTAR_1D

  REAL, DIMENSION(SIZE(PU))                    :: ZZ_O_LMO
  REAL, DIMENSION(SIZE(PU))                    :: ZZ0_O_LMO
!
!* purely unstable case
  USTAR_1D(:) = 0.
  ZZ_O_LMO(:) = XUNDEF
  ZZ0_O_LMO(:) = XUNDEF
!
!* general case
  WHERE(ABS(PLMO) > 1.E-20 .AND. PLMO/=XUNDEF)
    ZZ_O_LMO  = PZ(:)  / PLMO(:)
    ZZ0_O_LMO = PZ0(:) / PLMO(:)
    USTAR_1D(:) = SQRT( PU(:)**2+PV(:)**2 )               &
                * XKARMAN / ( LOG(PZ(:)/PZ0(:))           &
                             - PAULSON_PSIM(ZZ_O_LMO(:))  &
                             + PAULSON_PSIM(ZZ0_O_LMO(:)) )
  END WHERE
!
!* purely neutral case
  WHERE(PLMO==XUNDEF)
    ZZ_O_LMO = 0.
    USTAR_1D(:) = SQRT( PU(:)**2+PV(:)**2 )      &
                  * XKARMAN / LOG(PZ(:)/PZ0(:))
  END WHERE
!
END FUNCTION USTAR_1D
!
!-------------------------------------------------------------------------------
!
FUNCTION USTAR_0D(PU,PV,PZ,PZ0,PLMO)
  USE MODD_CST
  USE MODD_PARAMETERS
  REAL, INTENT(IN)               :: PU
  REAL, INTENT(IN)               :: PV
  REAL, INTENT(IN)               :: PZ
  REAL, INTENT(IN)               :: PZ0
  REAL, INTENT(IN)               :: PLMO
  REAL                           :: USTAR_0D
!
!* purely unstable case
  USTAR_0D = 0.
!
!* general case
  IF ( ABS(PLMO) >= 1.E-20 .AND. PLMO/=XUNDEF)    &
  USTAR_0D = SQRT( PU**2+PV**2 )                  &
             * XKARMAN / ( LOG(PZ/PZ0)            &
                          - PAULSON_PSIM(PZ/PLMO) &
                          + PAULSON_PSIM(PZ0/PLMO))
!
!* purely neutral case
  IF (PLMO==XUNDEF)                  &
  USTAR_0D = SQRT( PU**2+PV**2 )     &
             * XKARMAN / LOG(PZ/PZ0)

END FUNCTION USTAR_0D
!
!-------------------------------------------------------------------------------
!
END MODULE MODE_SBL
