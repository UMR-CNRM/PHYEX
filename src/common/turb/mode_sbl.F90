!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
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
!!      V. Masson       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original        13/10/99
!!      V. Masson       06/11/02 optimization and add Businger fonction for TKE
!!      V. Masson       01/01/03 use PAULSON_PSIM function
!-----------------------------------------------------------------------------
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
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
SUBROUTINE BUSINGER_PHIM_3D(PZ_O_LMO,BUSINGER_PHIM3D)
  REAL, DIMENSION(:,:,:), INTENT(IN)                 :: PZ_O_LMO
  REAL, DIMENSION(SIZE(PZ_O_LMO,1), &
                  SIZE(PZ_O_LMO,2),SIZE(PZ_O_LMO,3)),INTENT(OUT) :: BUSINGER_PHIM3D
!
  REAL(KIND=JPRB) :: ZHOOK_HANDLE
  IF (LHOOK) CALL DR_HOOK('MODE_SBL:BUSINGER_PHIM_3D',0,ZHOOK_HANDLE)
  WHERE ( PZ_O_LMO(:,:,:) < 0. )
    BUSINGER_PHIM3D(:,:,:) = (1.-15.*PZ_O_LMO)**(-0.25)
  ELSEWHERE
    BUSINGER_PHIM3D(:,:,:) = 1. + 4.7 * PZ_O_LMO
  END WHERE
  IF (LHOOK) CALL DR_HOOK('MODE_SBL:BUSINGER_PHIM_3D',1,ZHOOK_HANDLE)
END SUBROUTINE BUSINGER_PHIM_3D
!
!-------------------------------------------------------------------------------
!
SUBROUTINE BUSINGER_PHIM_2D(PZ_O_LMO,BUSINGER_PHIM2D)
  REAL, DIMENSION(:,:), INTENT(IN)                   :: PZ_O_LMO
  REAL, DIMENSION(SIZE(PZ_O_LMO,1),SIZE(PZ_O_LMO,2)),INTENT(OUT) :: BUSINGER_PHIM2D
!
  REAL(KIND=JPRB) :: ZHOOK_HANDLE
  IF (LHOOK) CALL DR_HOOK('MODE_SBL:BUSINGER_PHIM_2D',0,ZHOOK_HANDLE)
  WHERE ( PZ_O_LMO(:,:) < 0. )
    BUSINGER_PHIM2D(:,:) = (1.-15.*PZ_O_LMO)**(-0.25)
  ELSEWHERE
    BUSINGER_PHIM2D(:,:) = 1. + 4.7 * PZ_O_LMO
  END WHERE
  IF (LHOOK) CALL DR_HOOK('MODE_SBL:BUSINGER_PHIM_2D',1,ZHOOK_HANDLE)
END SUBROUTINE BUSINGER_PHIM_2D
!
!-------------------------------------------------------------------------------
!
SUBROUTINE BUSINGER_PHIM_1D(PZ_O_LMO,BUSINGER_PHIM1D)
  REAL, DIMENSION(:), INTENT(IN)  :: PZ_O_LMO
  REAL, DIMENSION(SIZE(PZ_O_LMO)),INTENT(OUT) :: BUSINGER_PHIM1D
!
  REAL(KIND=JPRB) :: ZHOOK_HANDLE
  IF (LHOOK) CALL DR_HOOK('MODE_SBL:BUSINGER_PHIM_1D',0,ZHOOK_HANDLE)
  WHERE ( PZ_O_LMO(:) < 0. )
    BUSINGER_PHIM1D(:) = (1.-15.*PZ_O_LMO)**(-0.25)
  ELSEWHERE
    BUSINGER_PHIM1D(:) = 1. + 4.7 * PZ_O_LMO
  END WHERE
  IF (LHOOK) CALL DR_HOOK('MODE_SBL:BUSINGER_PHIM_1D',1,ZHOOK_HANDLE)
END SUBROUTINE BUSINGER_PHIM_1D
!
!-------------------------------------------------------------------------------
!
SUBROUTINE BUSINGER_PHIM_0D(PZ_O_LMO,BUSINGER_PHIM0D)
  REAL, INTENT(IN)                   :: PZ_O_LMO
  REAL,INTENT(OUT)                   :: BUSINGER_PHIM0D
!
  REAL(KIND=JPRB) :: ZHOOK_HANDLE
  IF (LHOOK) CALL DR_HOOK('MODE_SBL:BUSINGER_PHIM_0D',0,ZHOOK_HANDLE)
  IF ( PZ_O_LMO < 0. ) THEN
    BUSINGER_PHIM0D = (1.-15.*PZ_O_LMO)**(-0.25)
  ELSE
    BUSINGER_PHIM0D = 1. + 4.7 * PZ_O_LMO
  END IF
  IF (LHOOK) CALL DR_HOOK('MODE_SBL:BUSINGER_PHIM_0D',1,ZHOOK_HANDLE)
END SUBROUTINE BUSINGER_PHIM_0D
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
SUBROUTINE BUSINGER_PHIH_3D(PZ_O_LMO,BUSINGER_PHIH3D)
  REAL, DIMENSION(:,:,:), INTENT(IN)                 :: PZ_O_LMO
  REAL, DIMENSION(SIZE(PZ_O_LMO,1), &
                  SIZE(PZ_O_LMO,2),SIZE(PZ_O_LMO,3)),INTENT(OUT) :: BUSINGER_PHIH3D
!
  REAL(KIND=JPRB) :: ZHOOK_HANDLE
  IF (LHOOK) CALL DR_HOOK('MODE_SBL:BUSINGER_PHIH_3D',0,ZHOOK_HANDLE)
  WHERE ( PZ_O_LMO(:,:,:) < 0. )
    BUSINGER_PHIH3D(:,:,:) = 0.74 * (1.-9.*PZ_O_LMO)**(-0.5)
  ELSEWHERE
    BUSINGER_PHIH3D(:,:,:) = 0.74 + 4.7 * PZ_O_LMO
  END WHERE
  IF (LHOOK) CALL DR_HOOK('MODE_SBL:BUSINGER_PHIH_3D',1,ZHOOK_HANDLE)
END SUBROUTINE BUSINGER_PHIH_3D
!
!-------------------------------------------------------------------------------
!
SUBROUTINE BUSINGER_PHIH_2D(PZ_O_LMO,BUSINGER_PHIH2D)
  REAL, DIMENSION(:,:), INTENT(IN)                   :: PZ_O_LMO
  REAL, DIMENSION(SIZE(PZ_O_LMO,1),SIZE(PZ_O_LMO,2)),INTENT(OUT) :: BUSINGER_PHIH2D
!
  REAL(KIND=JPRB) :: ZHOOK_HANDLE
  IF (LHOOK) CALL DR_HOOK('MODE_SBL:BUSINGER_PHIH_2D',0,ZHOOK_HANDLE)
  WHERE ( PZ_O_LMO(:,:) < 0. )
    BUSINGER_PHIH2D(:,:) = 0.74 * (1.-9.*PZ_O_LMO)**(-0.5)
  ELSEWHERE
    BUSINGER_PHIH2D(:,:) = 0.74 + 4.7 * PZ_O_LMO
  END WHERE
  IF (LHOOK) CALL DR_HOOK('MODE_SBL:BUSINGER_PHIH_2D',1,ZHOOK_HANDLE)
END SUBROUTINE BUSINGER_PHIH_2D
!
!-------------------------------------------------------------------------------
!
SUBROUTINE BUSINGER_PHIH_1D(PZ_O_LMO,BUSINGER_PHIH1D)
  REAL, DIMENSION(:), INTENT(IN)  :: PZ_O_LMO
  REAL, DIMENSION(SIZE(PZ_O_LMO)),INTENT(OUT) :: BUSINGER_PHIH1D
!
  REAL(KIND=JPRB) :: ZHOOK_HANDLE
  IF (LHOOK) CALL DR_HOOK('MODE_SBL:BUSINGER_PHIH_1D',0,ZHOOK_HANDLE)
  WHERE ( PZ_O_LMO(:) < 0. )
    BUSINGER_PHIH1D(:) = 0.74 * (1.-9.*PZ_O_LMO)**(-0.5)
  ELSEWHERE
    BUSINGER_PHIH1D(:) = 0.74 + 4.7 * PZ_O_LMO
  END WHERE
  IF (LHOOK) CALL DR_HOOK('MODE_SBL:BUSINGER_PHIH_1D',1,ZHOOK_HANDLE)
END SUBROUTINE BUSINGER_PHIH_1D
!
!-------------------------------------------------------------------------------
!
SUBROUTINE BUSINGER_PHIH_0D(PZ_O_LMO,BUSINGER_PHIH0D)
  REAL, INTENT(IN)                   :: PZ_O_LMO
  REAL,INTENT(OUT)                   :: BUSINGER_PHIH0D
!
  REAL(KIND=JPRB) :: ZHOOK_HANDLE
  IF (LHOOK) CALL DR_HOOK('MODE_SBL:BUSINGER_PHIH_0D',0,ZHOOK_HANDLE)
  IF ( PZ_O_LMO < 0. ) THEN
    BUSINGER_PHIH0D = 0.74 * (1.-9.*PZ_O_LMO)**(-0.5)
  ELSE
    BUSINGER_PHIH0D = 0.74 + 4.7 * PZ_O_LMO
  END IF
  IF (LHOOK) CALL DR_HOOK('MODE_SBL:BUSINGER_PHIH_0D',1,ZHOOK_HANDLE)
END SUBROUTINE BUSINGER_PHIH_0D
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
SUBROUTINE BUSINGER_PHIE_3D(PZ_O_LMO,BUSINGER_PHIE3D)
  USE MODD_CTURB
  REAL, DIMENSION(:,:,:), INTENT(IN)                 :: PZ_O_LMO
  REAL, DIMENSION(SIZE(PZ_O_LMO,1), &
                  SIZE(PZ_O_LMO,2),SIZE(PZ_O_LMO,3)),INTENT(OUT) :: BUSINGER_PHIE3D
!
  REAL(KIND=JPRB) :: ZHOOK_HANDLE
  IF (LHOOK) CALL DR_HOOK('MODE_SBL:BUSINGER_PHIE_3D',0,ZHOOK_HANDLE)
  WHERE ( PZ_O_LMO(:,:,:) < 0. )
    BUSINGER_PHIE3D(:,:,:) =   (1.+(-PZ_O_LMO)**(2./3.)/XALPSBL) &
                              * (1.-15.*PZ_O_LMO)**(0.5)
  ELSEWHERE
    BUSINGER_PHIE3D(:,:,:) = 1./(1. + 4.7 * PZ_O_LMO)**2
  END WHERE
  IF (LHOOK) CALL DR_HOOK('MODE_SBL:BUSINGER_PHIE_3D',1,ZHOOK_HANDLE)
END SUBROUTINE BUSINGER_PHIE_3D
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
SUBROUTINE PAULSON_PSIM_2D(PZ_O_LMO,PAULSON_PSIM2D)
  USE MODD_CST
  REAL, DIMENSION(:,:), INTENT(IN)                   :: PZ_O_LMO
  REAL, DIMENSION(SIZE(PZ_O_LMO,1),SIZE(PZ_O_LMO,2)),INTENT(OUT) :: PAULSON_PSIM2D
!
  REAL, DIMENSION(SIZE(PZ_O_LMO,1),SIZE(PZ_O_LMO,2)) :: ZX

  REAL(KIND=JPRB) :: ZHOOK_HANDLE
  IF (LHOOK) CALL DR_HOOK('MODE_SBL:PAULSON_PSIM_2D',0,ZHOOK_HANDLE)
  ZX=1.
  WHERE ( PZ_O_LMO(:,:) < 0. )
    ZX=(1.-15.*PZ_O_LMO)**(0.25)
    PAULSON_PSIM2D(:,:) = LOG( (1.+ZX**2)*(1+ZX)**2/8. ) - 2.*ATAN(ZX) + XPI/2.
  ELSEWHERE
    PAULSON_PSIM2D(:,:) = - 4.7 * PZ_O_LMO
  END WHERE
  IF (LHOOK) CALL DR_HOOK('MODE_SBL:PAULSON_PSIM_2D',1,ZHOOK_HANDLE)
END SUBROUTINE PAULSON_PSIM_2D
!
!-------------------------------------------------------------------------------
!
SUBROUTINE PAULSON_PSIM_1D(PZ_O_LMO,PAULSON_PSIM1D)
  USE MODD_CST
  REAL, DIMENSION(:), INTENT(IN)    :: PZ_O_LMO
  REAL, DIMENSION(SIZE(PZ_O_LMO,1)),INTENT(OUT) :: PAULSON_PSIM1D
!
  REAL, DIMENSION(SIZE(PZ_O_LMO,1)) :: ZX

  REAL(KIND=JPRB) :: ZHOOK_HANDLE
  IF (LHOOK) CALL DR_HOOK('MODE_SBL:PAULSON_PSIM_1D',0,ZHOOK_HANDLE)
  ZX=1.
  WHERE ( PZ_O_LMO(:) < 0. )
    ZX=(1.-15.*PZ_O_LMO)**(0.25)
    PAULSON_PSIM1D(:) = LOG( (1.+ZX**2)*(1+ZX)**2/8. ) - 2.*ATAN(ZX) + XPI/2.
  ELSEWHERE
    PAULSON_PSIM1D(:) = - 4.7 * PZ_O_LMO
  END WHERE
  IF (LHOOK) CALL DR_HOOK('MODE_SBL:PAULSON_PSIM_1D',1,ZHOOK_HANDLE)
END SUBROUTINE PAULSON_PSIM_1D
!
!-------------------------------------------------------------------------------
!
SUBROUTINE PAULSON_PSIM_0D(PZ_O_LMO,PAULSON_PSIM0D)
  USE MODD_CST
  REAL, INTENT(IN)    :: PZ_O_LMO
  REAL,INTENT(OUT)    :: PAULSON_PSIM0D
!
  REAL                :: ZX

  REAL(KIND=JPRB) :: ZHOOK_HANDLE
  IF (LHOOK) CALL DR_HOOK('MODE_SBL:PAULSON_PSIM_0D',0,ZHOOK_HANDLE)
  ZX=1.
  IF ( PZ_O_LMO < 0. ) THEN
    ZX=(1.-15.*PZ_O_LMO)**(0.25)
    PAULSON_PSIM0D = LOG( (1.+ZX**2)*(1+ZX)**2/8. ) - 2.*ATAN(ZX) + XPI/2.
  ELSE
    PAULSON_PSIM0D = - 4.7 * PZ_O_LMO
  END IF
  IF (LHOOK) CALL DR_HOOK('MODE_SBL:PAULSON_PSIM_0D',1,ZHOOK_HANDLE)
END SUBROUTINE PAULSON_PSIM_0D
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
SUBROUTINE LMO_2D(PUSTAR,PTHETA,PRV,PSFTH,PSFRV,LMO2D)
  USE MODD_CST
  USE MODD_PARAMETERS, ONLY: JPVEXT_TURB,XUNDEF
  REAL, DIMENSION(:,:), INTENT(IN)               :: PUSTAR
  REAL, DIMENSION(:,:), INTENT(IN)               :: PTHETA
  REAL, DIMENSION(:,:), INTENT(IN)               :: PRV
  REAL, DIMENSION(:,:), INTENT(IN)               :: PSFTH
  REAL, DIMENSION(:,:), INTENT(IN)               :: PSFRV
  REAL, DIMENSION(SIZE(PUSTAR,1),SIZE(PUSTAR,2)),INTENT(OUT) :: LMO2D
!
  REAL, DIMENSION(SIZE(PUSTAR,1),SIZE(PUSTAR,2)) :: ZTHETAV
  REAL, DIMENSION(SIZE(PUSTAR,1),SIZE(PUSTAR,2)) :: ZQ0
  REAL                                           :: ZEPS
!
!
  REAL(KIND=JPRB) :: ZHOOK_HANDLE
  IF (LHOOK) CALL DR_HOOK('MODE_SBL:LMO_2D',0,ZHOOK_HANDLE)
  ZEPS=(XRV-XRD)/XRD
  ZTHETAV(:,:) = PTHETA(:,:) * ( 1. +ZEPS * PRV(:,:))
  ZQ0    (:,:) = PSFTH(:,:) + ZTHETAV(:,:) * ZEPS * PSFRV(:,:)
!
  LMO2D(:,:) = XUNDEF
  WHERE ( ZQ0(:,:) /=0.  )                                   &
    LMO2D(:,:) = - MAX(PUSTAR(:,:),1.E-6)**3                &
                  / ( XKARMAN * XG / ZTHETAV(:,:) *ZQ0(:,:) )

  IF (LHOOK) CALL DR_HOOK('MODE_SBL:LMO_2D',1,ZHOOK_HANDLE)
END SUBROUTINE LMO_2D
!
!-------------------------------------------------------------------------------
!
SUBROUTINE LMO_1D(PUSTAR,PTHETA,PRV,PSFTH,PSFRV,LMO1D)
  USE MODD_CST
  USE MODD_PARAMETERS, ONLY: JPVEXT_TURB,XUNDEF
  REAL, DIMENSION(:), INTENT(IN)  :: PUSTAR
  REAL, DIMENSION(:), INTENT(IN)  :: PTHETA
  REAL, DIMENSION(:), INTENT(IN)  :: PRV
  REAL, DIMENSION(:), INTENT(IN)  :: PSFTH
  REAL, DIMENSION(:), INTENT(IN)  :: PSFRV
  REAL, DIMENSION(SIZE(PUSTAR)),INTENT(OUT)   :: LMO1D
!
  REAL, DIMENSION(SIZE(PUSTAR))   :: ZTHETAV
  REAL                                           :: ZEPS
!
!
  REAL(KIND=JPRB) :: ZHOOK_HANDLE
  IF (LHOOK) CALL DR_HOOK('MODE_SBL:LMO_1D',0,ZHOOK_HANDLE)
  ZEPS=(XRV-XRD)/XRD
!
  ZTHETAV(:) = PTHETA(:) * ( 1. +ZEPS * PRV(:))
!
  LMO1D(:) = XUNDEF
  WHERE ( PSFTH(:)/ZTHETAV(:)+ZEPS*PSFRV(:)/=0. )                  &
    LMO1D(:) = - MAX(PUSTAR(:),1.E-6)**3                          &
              / ( XKARMAN * (  XG / ZTHETAV(:)    * PSFTH(:)       &
                             + XG * ZEPS * PSFRV(:) )  )
  IF (LHOOK) CALL DR_HOOK('MODE_SBL:LMO_1D',1,ZHOOK_HANDLE)
END SUBROUTINE LMO_1D
!
!-------------------------------------------------------------------------------
!
SUBROUTINE LMO_0D(PUSTAR,PTHETA,PRV,PSFTH,PSFRV,LMO0D)
  USE MODD_CST
  USE MODD_PARAMETERS, ONLY: JPVEXT_TURB,XUNDEF
  REAL, INTENT(IN)  :: PUSTAR
  REAL, INTENT(IN)  :: PTHETA
  REAL, INTENT(IN)  :: PRV
  REAL, INTENT(IN)  :: PSFTH
  REAL, INTENT(IN)  :: PSFRV
  REAL, INTENT(OUT) :: LMO0D
!
  REAL              :: ZTHETAV
  REAL              :: ZEPS
!
!
  REAL(KIND=JPRB) :: ZHOOK_HANDLE
  IF (LHOOK) CALL DR_HOOK('MODE_SBL:LMO_0D',0,ZHOOK_HANDLE)
  ZEPS=(XRV-XRD)/XRD
!
!
  ZTHETAV = PTHETA * ( 1. +ZEPS * PRV)
!
  LMO0D = XUNDEF
  IF ( PSFTH/ZTHETAV+ZEPS*PSFRV/=0. )                      &
  LMO0D = - MAX(PUSTAR,1.E-6)**3                          &
           / ( XKARMAN * (  XG / ZTHETAV       * PSFTH     &
                          + XG * ZEPS * PSFRV )  )
  IF (LHOOK) CALL DR_HOOK('MODE_SBL:LMO_0D',1,ZHOOK_HANDLE)
END SUBROUTINE LMO_0D
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
SUBROUTINE USTAR_2D(PU,PV,PZ,PZ0,PLMO,USTAR2D)
  USE MODD_CST
  USE MODD_PARAMETERS, ONLY: JPVEXT_TURB,XUNDEF
  REAL, DIMENSION(:,:), INTENT(IN)               :: PU
  REAL, DIMENSION(:,:), INTENT(IN)               :: PV
  REAL, DIMENSION(:,:), INTENT(IN)               :: PZ
  REAL, DIMENSION(:,:), INTENT(IN)               :: PZ0
  REAL, DIMENSION(:,:), INTENT(IN)               :: PLMO
  REAL, DIMENSION(SIZE(PU,1),SIZE(PU,2)),INTENT(OUT) :: USTAR2D

  REAL, DIMENSION(SIZE(PU,1),SIZE(PU,2))         :: ZZ_O_LMO
  REAL, DIMENSION(SIZE(PU,1),SIZE(PU,2))         :: ZZ0_O_LMO
  REAL, DIMENSION(SIZE(PU,1),SIZE(PU,2))         :: ZWORK1,ZWORK2
!
!* purely unstable case
  REAL(KIND=JPRB) :: ZHOOK_HANDLE
  IF (LHOOK) CALL DR_HOOK('MODE_SBL:USTAR_2D',0,ZHOOK_HANDLE)
  USTAR2D(:,:) = 0.
  ZZ_O_LMO(:,:) = XUNDEF
  ZZ0_O_LMO(:,:) = XUNDEF
!
!* general case
  CALL PAULSON_PSIM(ZZ_O_LMO,ZWORK1)
  CALL PAULSON_PSIM(ZZ0_O_LMO,ZWORK2)
  WHERE(ABS(PLMO) > 1.E-20 .AND. PLMO/=XUNDEF)
    ZZ_O_LMO  = PZ(:,:)  / PLMO(:,:)
    ZZ0_O_LMO = PZ0(:,:) / PLMO(:,:)      
    USTAR2D(:,:) = SQRT( PU(:,:)**2+PV(:,:)**2 )               &
                  * XKARMAN / ( LOG(PZ(:,:)/PZ0(:,:))           &
                                 - ZWORK1(:,:) + ZWORK2(:,:))
  END WHERE
!
!* purely neutral case
  WHERE(PLMO==XUNDEF)
    ZZ_O_LMO = 0.
    USTAR2D(:,:) = SQRT( PU(:,:)**2+PV(:,:)**2 )      &
                  * XKARMAN / LOG(PZ(:,:)/PZ0(:,:))
  END WHERE
!
  IF (LHOOK) CALL DR_HOOK('MODE_SBL:USTAR_2D',1,ZHOOK_HANDLE)
END SUBROUTINE USTAR_2D
!
!-------------------------------------------------------------------------------
!
SUBROUTINE USTAR_1D(PU,PV,PZ,PZ0,PLMO,USTAR1D)
  USE MODD_CST
  USE MODD_PARAMETERS, ONLY: JPVEXT_TURB,XUNDEF
  REAL, DIMENSION(:), INTENT(IN)               :: PU
  REAL, DIMENSION(:), INTENT(IN)               :: PV
  REAL, DIMENSION(:), INTENT(IN)               :: PZ
  REAL, DIMENSION(:), INTENT(IN)               :: PZ0
  REAL, DIMENSION(:), INTENT(IN)               :: PLMO
  REAL, DIMENSION(SIZE(PU)),INTENT(OUT)        :: USTAR1D

  REAL, DIMENSION(SIZE(PU))                    :: ZZ_O_LMO
  REAL, DIMENSION(SIZE(PU))                    :: ZZ0_O_LMO
  REAL, DIMENSION(SIZE(PU))                    :: ZWORK1,ZWORK2
!
!* purely unstable case
  REAL(KIND=JPRB) :: ZHOOK_HANDLE
  IF (LHOOK) CALL DR_HOOK('MODE_SBL:USTAR_1D',0,ZHOOK_HANDLE)
  USTAR1D(:) = 0.
  ZZ_O_LMO(:) = XUNDEF
  ZZ0_O_LMO(:) = XUNDEF
!
!* general case
  CALL PAULSON_PSIM(ZZ_O_LMO,ZWORK1)
  CALL PAULSON_PSIM(ZZ0_O_LMO,ZWORK2)
  WHERE(ABS(PLMO) > 1.E-20 .AND. PLMO/=XUNDEF)
    ZZ_O_LMO  = PZ(:)  / PLMO(:)
    ZZ0_O_LMO = PZ0(:) / PLMO(:)
    USTAR1D(:) = SQRT( PU(:)**2+PV(:)**2 )               &
                * XKARMAN / ( LOG(PZ(:)/PZ0(:))           &
                             - ZWORK1(:) + ZWORK2(:))
  END WHERE
!
!* purely neutral case
  WHERE(PLMO==XUNDEF)
    ZZ_O_LMO = 0.
    USTAR1D(:) = SQRT( PU(:)**2+PV(:)**2 )      &
                  * XKARMAN / LOG(PZ(:)/PZ0(:))
  END WHERE
!
  IF (LHOOK) CALL DR_HOOK('MODE_SBL:USTAR_1D',1,ZHOOK_HANDLE)
END SUBROUTINE USTAR_1D
!
!-------------------------------------------------------------------------------
!
SUBROUTINE USTAR_0D(PU,PV,PZ,PZ0,PLMO,USTAR0D)
  USE MODD_CST
  USE MODD_PARAMETERS, ONLY: JPVEXT_TURB,XUNDEF
  REAL, INTENT(IN)               :: PU
  REAL, INTENT(IN)               :: PV
  REAL, INTENT(IN)               :: PZ
  REAL, INTENT(IN)               :: PZ0
  REAL, INTENT(IN)               :: PLMO
  REAL, INTENT(OUT)              :: USTAR0D
  REAL :: ZWORK, ZWORK2
!
!* purely unstable case
  REAL(KIND=JPRB) :: ZHOOK_HANDLE
  IF (LHOOK) CALL DR_HOOK('MODE_SBL:USTAR_0D',0,ZHOOK_HANDLE)
  USTAR0D = 0.
!
!* general case
  IF ( ABS(PLMO) >= 1.E-20 .AND. PLMO/=XUNDEF) THEN
    CALL PAULSON_PSIM(PZ/PLMO,ZWORK1)
    CALL PAULSON_PSIM(PZ0/PLMO,ZWORK2)
    USTAR0D = SQRT( PU**2+PV**2 )                  &
             * XKARMAN / ( LOG(PZ/PZ0)            &
                          - ZWORK1 + ZWORK2)
  END IF
!
!* purely neutral case
  IF (PLMO==XUNDEF)                  &
  USTAR0D = SQRT( PU**2+PV**2 )     &
             * XKARMAN / LOG(PZ/PZ0)

  IF (LHOOK) CALL DR_HOOK('MODE_SBL:USTAR_0D',1,ZHOOK_HANDLE)
END SUBROUTINE USTAR_0D
!
!-------------------------------------------------------------------------------
!
END MODULE MODE_SBL
