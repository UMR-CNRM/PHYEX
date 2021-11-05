!MNH_LIC Copyright 2002-2020 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ################
      MODULE MODI_RMC01
!     ################
INTERFACE
      SUBROUTINE RMC01(HTURBLEN,KKA,KKU,KKL,PZZ,PDXX,PDYY,PDZZ,PDIRCOSZW, &
                       PSBL_DEPTH, PLMO, PLK, PLEPS                )
!
CHARACTER(LEN=4),         INTENT(IN)    :: HTURBLEN ! type of mixing length
INTEGER,                  INTENT(IN)    :: KKA           !near ground array index  
INTEGER,                  INTENT(IN)    :: KKU           !uppest atmosphere array index
INTEGER,                  INTENT(IN)    :: KKL           !vert. levels type 1=MNH -1=ARO
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PZZ   ! altitude of flux points
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PDXX  ! width of grid mesh (X dir)
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PDYY  ! width of grid mesh (Y dir)
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PDZZ  ! width of vert. layers
REAL, DIMENSION(:,:),     INTENT(IN)    :: PDIRCOSZW ! Director Cosinus 
REAL, DIMENSION(:,:),     INTENT(IN)    :: PSBL_DEPTH! SBL depth
REAL, DIMENSION(:,:),     INTENT(IN)    :: PLMO  ! Monin Obuhkov length
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PLK   ! Mixing length
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PLEPS ! Dissipative length

END SUBROUTINE RMC01
END INTERFACE
END MODULE MODI_RMC01
!
!     ##############################################################
      SUBROUTINE RMC01(HTURBLEN,KKA, KKU, KKL, PZZ, PDXX, PDYY, PDZZ, PDIRCOSZW, &
                       PSBL_DEPTH, PLMO, PLK, PLEPS                )
!     ##############################################################
!
!!****  *RMC01* -
!! 
!!    PURPOSE
!!    -------
!!    This routine modifies the mixing and dissipative length near the SBL.
!!    (Redelsperger, Mahe and Carlotti, 2001)
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!      Book 2
!!
!!    AUTHOR
!!    ------
!!	
!!      V. Masson  - Meteo-France -
!!
!!    MODIFICATIONS
!!    -------------
!!     Original     14/02/02
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAMETERS
USE MODD_CST
USE MODD_CTURB
!
USE MODE_SBL
!
USE MODI_SHUMAN
!
IMPLICIT NONE
!
!*       0.1   Declaration of arguments
!              ------------------------
!
CHARACTER(LEN=4),         INTENT(IN)    :: HTURBLEN ! type of mixing length
INTEGER,                   INTENT(IN)   :: KKA      !near ground array index  
INTEGER,                   INTENT(IN)   :: KKU      !uppest atmosphere array index
INTEGER,                   INTENT(IN)   :: KKL      !vert. levels type 1=MNH -1=ARO
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PZZ   ! altitude of flux points
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PDXX  ! width of grid mesh (X dir)
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PDYY  ! width of grid mesh (Y dir)
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PDZZ  ! width of vert. layers
REAL, DIMENSION(:,:),     INTENT(IN)    :: PDIRCOSZW ! Director Cosinus 
REAL, DIMENSION(:,:),     INTENT(IN)    :: PSBL_DEPTH! SBL depth
REAL, DIMENSION(:,:),     INTENT(IN)    :: PLMO  ! Monin Obuhkov length
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PLK   ! Mixing length
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PLEPS ! Dissipative length
!
!*       0.2   Declaration of local variables
!              ------------------------------
!
INTEGER :: IKB,IKE    ! first,last physical level
INTEGER :: IKT        ! array size in k direction
INTEGER :: IKTB,IKTE  ! start, end of k loops in physical domain 
INTEGER :: IIU        ! horizontal x boundary
INTEGER :: IJU        ! horizontal y boundary
INTEGER :: JK         ! loop counter
!
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)) :: ZZZ  ! height of mass
                                                             ! points above ground
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)) :: ZZ_O_LMO ! height / LMO
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)) :: ZGAM ! factor controling
                                                             ! transition betw.
                                                             ! SBL and free BL

REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)) :: ZPHIM! MO function
                                                             ! for stress
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)) :: ZPHIE! MO function
                                                             ! for TKE
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)) :: ZDH  ! hor. grid mesh
                                                             ! size
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)) :: ZL   ! SBL length
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)) :: ZZC  ! alt. where
                                                             ! turb. is isotr.
!-------------------------------------------------------------------------------
!
!*     1. Initializations
!         ---------------
!
! horizontal boundaries
IIU=SIZE(PZZ,1)
IJU=SIZE(PZZ,2)
!
! vertical boundaries
IKB=KKA+JPVEXT_TURB*KKL
IKE=KKU-JPVEXT_TURB*KKL

IKTB=1+JPVEXT_TURB
IKT=SIZE(PZZ,3)
IKTE=IKT-JPVEXT_TURB
!
! altitude of mass points
ZZZ=MZF(PZZ)
! replace by height of mass points
DO JK=1,IKT
  ZZZ(:,:,JK) = ZZZ(:,:,JK) - PZZ(:,:,IKB)
END DO
! fill upper level with physical value
ZZZ(:,:,KKU) = 2.*ZZZ(:,:,KKU-KKL) - ZZZ(:,:,KKU-2*KKL)
!
!-------------------------------------------------------------------------------
!
!*     2. MO quantities
!         -------------
!
! z/LMO
DO JK=1,IKT
  WHERE (PLMO(:,:)==XUNDEF)
    ZZ_O_LMO(:,:,JK)=0.
  ELSEWHERE
    ZZ_O_LMO(:,:,JK)=ZZZ(:,:,JK)*PDIRCOSZW(:,:)/PLMO(:,:) 
  END WHERE
END DO
ZZ_O_LMO(:,:,:) = MAX(ZZ_O_LMO(:,:,:),-10.)
ZZ_O_LMO(:,:,:) = MIN(ZZ_O_LMO(:,:,:), 10.)
!
!
! MO function for stress
ZPHIM(:,:,:) = BUSINGER_PHIM(ZZ_O_LMO)
!
! MO function for TKE
ZPHIE(:,:,:) = BUSINGER_PHIE(ZZ_O_LMO)
!
!-------------------------------------------------------------------------------
SELECT CASE (HTURBLEN)
!-------------------------------------------------------------------------------
!
!*     3. altitude where turbulence is isotropic inside a layer of given width (3D case)
!         --------------------------------------------------------------------
!
!
!* LES subgrid mixing (unresolved eddies all below mesh size)
!  For stable cases, the vertical size of eddies is supposed to be given by the
!  same law as in the neutral case (i.e. with Phim = 1).
!
  CASE ('DELT','DEAR')
    ZDH = SQRT(MXF(PDXX)*MYF(PDYY))
    ZDH(IIU,:,:) = ZDH(IIU-1,:,:)
    ZDH(:,IJU,:) = ZDH(:,IJU-1,:)
    DO JK=1,IKT
      ZZC(:,:,JK) = 2.*MIN(ZPHIM(:,:,JK),1.)/XKARMAN    &
                     * MAX( PDZZ(:,:,JK)*PDIRCOSZW(:,:) , ZDH(:,:,JK)/PDIRCOSZW(:,:)/3. )
    END DO
!
!*     4. factor controling the transition between SBL and free isotropic turb. (3D case)
!         --------------------------------------------------------------------
!
    ZGAM(:,:,KKA) = 0.
    DO JK=IKTB,IKTE
      ZGAM(:,:,JK) = 1.  - EXP( -3.*(ZZZ(:,:,JK)-ZZZ(:,:,IKB))/(ZZC(:,:,JK)) )
      WHERE (ZGAM(:,:,JK-KKL)>ZGAM(:,:,JK) .OR. ZGAM(:,:,JK-KKL)>0.99 ) ZGAM(:,:,JK) = 1.
    END DO
    ZGAM(:,:,KKU) = 1.  - EXP( -3.*(ZZZ(:,:,KKU)-ZZZ(:,:,IKB))/(ZZC(:,:,KKU)) )
    WHERE (ZGAM(:,:,KKU-KKL)>ZGAM(:,:,KKU) .OR. ZGAM(:,:,KKU-KKL)>0.99 ) ZGAM(:,:,KKU) = 1.
!
!
!-------------------------------------------------------------------------------
!
!*     5. factor controling the transition between SBL and free isotropic turb.(1D case)
!         --------------------------------------------------------------------
!
  CASE DEFAULT
!* SBL depth is used
    ZGAM(:,:,:) = 1.
    ZGAM(:,:,KKA) = 0.
    DO JK=IKTB,IKTE
      WHERE(PSBL_DEPTH>0.) &
      ZGAM(:,:,JK) = TANH( (ZZZ(:,:,JK)-ZZZ(:,:,IKB))/PSBL_DEPTH(:,:) )
      WHERE (ZGAM(:,:,JK-KKL)>0.99 ) ZGAM(:,:,JK) = 1.
    END DO
    WHERE(PSBL_DEPTH>0.) &
      ZGAM(:,:,KKU) = TANH( (ZZZ(:,:,KKU)-ZZZ(:,:,IKB))/PSBL_DEPTH(:,:) )
    WHERE (ZGAM(:,:,KKU-KKL)>0.99 ) ZGAM(:,:,JK) = 1.
!
!-------------------------------------------------------------------------------
END SELECT
!-------------------------------------------------------------------------------
!
!*     6. Modification of the mixing length
!         ---------------------------------
!
DO JK=1,IKT
  ZL(:,:,JK) =  XKARMAN/SQRT(XALPSBL)/XCMFS                                      &
              * ZZZ(:,:,JK)*PDIRCOSZW(:,:)/(ZPHIM(:,:,JK)**2*SQRT(ZPHIE(:,:,JK)))
END DO
!
PLK(:,:,:)=(1.-ZGAM)*ZL+ZGAM*PLK
!
PLK(:,:,KKA) = PLK(:,:,IKB)
PLK(:,:,KKU) = PLK(:,:,IKE)
!-------------------------------------------------------------------------------
!
!*     7. Modification of the dissipative length
!         --------------------------------------
!
ZL = ZL * (XALPSBL**(3./2.)*XKARMAN*XCED) &
        / (XKARMAN/SQRT(XALPSBL)/XCMFS)
!
WHERE (ZZ_O_LMO<0.)
  ZL = ZL/(1.-1.9*ZZ_O_LMO)
ELSEWHERE
  ZL = ZL/(1.-0.3*SQRT(ZZ_O_LMO))
ENDWHERE
!
PLEPS(:,:,:)=(1.-ZGAM)*ZL+ZGAM*PLEPS
!
PLEPS(:,:,KKA) = PLEPS(:,:,IKB)
PLEPS(:,:,KKU  ) = PLEPS(:,:,IKE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE RMC01
