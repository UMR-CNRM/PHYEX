!MNH_LIC Copyright 1994-2022 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
MODULE MODE_RMC01
IMPLICIT NONE
CONTAINS
SUBROUTINE RMC01(D,CST,CSTURB,TURBN,PZZ,PDXX,PDYY,PDZZ,PDIRCOSZW,PSBL_DEPTH,PLMO,PLK,PLEPS)
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
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
USE MODD_PARAMETERS, ONLY: XUNDEF
USE MODD_CST, ONLY : CST_t
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
USE MODD_CTURB, ONLY: CSTURB_t
USE MODD_TURB_n, ONLY: TURB_t
!
USE MODE_UPDATE_IIJU_PHY, ONLY: UPDATE_IIJU_PHY
USE MODE_SBL_PHY, ONLY: BUSINGER_PHIM, BUSINGER_PHIE
!
USE MODE_SHUMAN_PHY, ONLY: MZF_PHY, MYF_PHY, MXF_PHY
!
IMPLICIT NONE
!
!*       0.1   Declaration of arguments
!              ------------------------
!
TYPE(DIMPHYEX_t),         INTENT(IN)   :: D
TYPE(CST_t),              INTENT(IN)   :: CST
TYPE(CSTURB_t),           INTENT(IN)   :: CSTURB
TYPE(TURB_t),             INTENT(IN)   :: TURBN
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PZZ   ! altitude of flux points
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PDXX  ! width of grid mesh (X dir)
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PDYY  ! width of grid mesh (Y dir)
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PDZZ  ! width of vert. layers
REAL, DIMENSION(D%NIJT),     INTENT(IN)    :: PDIRCOSZW ! Director Cosinus
REAL, DIMENSION(D%NIJT),     INTENT(IN)    :: PSBL_DEPTH! SBL depth
REAL, DIMENSION(D%NIJT),     INTENT(IN)    :: PLMO  ! Monin Obuhkov length
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PLK   ! Mixing length
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PLEPS ! Dissipative length
!
!*       0.2   Declaration of local variables
!              ------------------------------
!
INTEGER :: IKB,IKE    ! first,last physical level
INTEGER :: IKTB,IKTE,IKT,IKA,IKU,IKL  ! start, end of k loops in physical domain
INTEGER :: JK,JIJ   ! loop counter
INTEGER :: IIJB,IIJE
!
REAL, DIMENSION(D%NIJT,D%NKT) :: ZZZ  ! height of mass
                                                             ! points above ground
REAL, DIMENSION(D%NIJT,D%NKT) :: ZZ_O_LMO ! height / LMO
REAL, DIMENSION(D%NIJT,D%NKT) :: ZGAM ! factor controling
                                                             ! transition betw.
                                                             ! SBL and free BL
REAL, DIMENSION(D%NIJT,D%NKT) :: ZPHIM! MO function
                                                             ! for stress
REAL, DIMENSION(D%NIJT,D%NKT) :: ZPHIE! MO function
                                                             ! for TKE
REAL, DIMENSION(D%NIJT,D%NKT) :: ZZC  ! alt. where turb. is isotr.
                                                             ! size
REAL, DIMENSION(D%NIJT,D%NKT) :: ZL   ! SBL length
REAL, DIMENSION(D%NIJT,D%NKT) :: ZWORK1, ZWORK2
REAL, DIMENSION(D%NIJT,D%NKT) :: ZDH  ! hor. grid mesh
!-------------------------------------------------------------------------------
!
!*     1. Initializations
!         ---------------
!
! horizontal boundaries
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('RMC01',0,ZHOOK_HANDLE)
IKTB=D%NKTB          
IKTE=D%NKTE
IKB=D%NKB
IKE=D%NKE
IIJB=D%NIJB
IIJE=D%NIJE
IKT=D%NKT
IKA=D%NKA
IKU=D%NKU
IKL=D%NKL
!
! altitude of mass points
CALL MZF_PHY(D,PZZ,ZZZ)
! replace by height of mass points
DO JK=1,IKT
  !$mnh_expand_array(JIJ=IIJB:IIJE)
  ZZZ(:,JK) = ZZZ(:,JK) - PZZ(:,IKB)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE)
END DO
! fill upper level with physical value
!$mnh_expand_array(JIJ=IIJB:IIJE)
ZZZ(:,IKU) = 2.*ZZZ(:,IKU-IKL) - ZZZ(:,IKU-2*IKL)
!$mnh_end_expand_array(JIJ=IIJB:IIJE)
!
!-------------------------------------------------------------------------------
!
!*     2. MO quantities
!         -------------
!
! z/LMO
DO JK=1,IKT
  !$mnh_expand_where(JIJ=IIJB:IIJE)
  WHERE (PLMO(:)==XUNDEF)
    ZZ_O_LMO(:,JK)=0.
  ELSEWHERE
    ZZ_O_LMO(:,JK)=ZZZ(:,JK)*PDIRCOSZW(:)/PLMO(:)
  END WHERE
  !$mnh_end_expand_where(JIJ=IIJB:IIJE)
END DO
!$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
ZZ_O_LMO(:,:) = MAX(ZZ_O_LMO(:,:),-10.)
ZZ_O_LMO(:,:) = MIN(ZZ_O_LMO(:,:), 10.)
!$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
!
!
! MO function for stress
CALL BUSINGER_PHIM(D,ZZ_O_LMO,ZPHIM)
!
! MO function for TKE
CALL BUSINGER_PHIE(D,CSTURB,ZZ_O_LMO,ZPHIE)
!
!-------------------------------------------------------------------------------
SELECT CASE (TURBN%CTURBLEN)
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
    CALL MXF_PHY(D,PDXX,ZWORK1)
    CALL MYF_PHY(D,PDYY,ZWORK2)
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ZDH(:,:) = SQRT(ZWORK1(:,:)*ZWORK2(:,:))
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    !
    CALL UPDATE_IIJU_PHY(D,ZZC)
    !
    DO JK=1,IKT
      !$mnh_expand_array(JIJ=IIJB:IIJE)
      ZZC(:,JK) = 2.*MIN(ZPHIM(:,JK),1.)/CST%XKARMAN    &
                            * MAX( PDZZ(:,JK)*PDIRCOSZW(:) , & 
                            ZDH(:,JK)/PDIRCOSZW(:)/3. )
      !$mnh_end_expand_array(JIJ=IIJB:IIJE)
    END DO
    !
    !*     4. factor controling the transition between SBL and free isotropic turb. (3D case)
    !         --------------------------------------------------------------------
    !
    ZGAM(:,IKA) = 0.
    DO JK=IKTB,IKTE
      !$mnh_expand_array(JIJ=IIJB:IIJE)
      ZGAM(:,JK) = 1.  - EXP( -3.*(ZZZ(:,JK)-ZZZ(:,IKB))/(ZZC(:,JK)) )
      !$mnh_end_expand_array(JIJ=IIJB:IIJE)
      !$mnh_expand_where(JIJ=IIJB:IIJE)
      WHERE (ZGAM(:,JK-IKL)>ZGAM(:,JK) .OR. ZGAM(:,JK-IKL)>0.99 ) 
        ZGAM(:,JK) = 1.
      END WHERE
     !$mnh_end_expand_where(JIJ=IIJB:IIJE)
    END DO
    !$mnh_expand_array(JIJ=IIJB:IIJE)
    ZGAM(:,IKU) = 1.  - EXP( -3.*(ZZZ(:,IKU)-ZZZ(:,IKB))& 
                                   /(ZZC(:,IKU)) )
    !$mnh_end_expand_array(JIJ=IIJB:IIJE)
    !$mnh_expand_where(JIJ=IIJB:IIJE)
    WHERE (ZGAM(:,IKU-IKL)>ZGAM(:,IKU) .OR. ZGAM(:,IKU-IKL)>0.99 ) 
      ZGAM(:,IKU) = 1.
    END WHERE
    !$mnh_end_expand_where(JIJ=IIJB:IIJE)
  !   
  !
  !-------------------------------------------------------------------------------
  !
  !*     5. factor controling the transition between SBL and free isotropic turb.(1D case)
  !         --------------------------------------------------------------------
  !
  CASE DEFAULT
    !* SBL depth is used
    ZGAM(:,:) = 1.
    ZGAM(:,IKA) = 0.
    DO JK=IKTB,IKTE
      !$mnh_expand_where(JIJ=IIJB:IIJE)
      WHERE(PSBL_DEPTH(:)>0.)
        ZGAM(:,JK) = TANH( (ZZZ(:,JK)-ZZZ(:,IKB))/PSBL_DEPTH(:) )
      END WHERE
      !$mnh_end_expand_where(JIJ=IIJB:IIJE)
      !$mnh_expand_where(JIJ=IIJB:IIJE)
      WHERE (ZGAM(:,JK-IKL)>0.99 ) 
        ZGAM(:,JK) = 1.
      END WHERE
      !$mnh_end_expand_where(JIJ=IIJB:IIJE)
    END DO
    !$mnh_expand_where(JIJ=IIJB:IIJE)
    WHERE(PSBL_DEPTH(:)>0.)
      ZGAM(:,IKU) = TANH( (ZZZ(:,IKU)-ZZZ(:,IKB))/PSBL_DEPTH(:) )
    END WHERE
    !$mnh_end_expand_where(JIJ=IIJB:IIJE)
    !$mnh_expand_where(JIJ=IIJB:IIJE)
    WHERE (ZGAM(:,IKU-IKL)>0.99 ) 
      ZGAM(:,JK) = 1.
    END WHERE
    !$mnh_end_expand_where(JIJ=IIJB:IIJE)
  !
  !-------------------------------------------------------------------------------
END SELECT
!-------------------------------------------------------------------------------
!
!*     6. Modification of the mixing length
!         ---------------------------------
!
DO JK=1,IKT
  !$mnh_expand_array(JIJ=IIJB:IIJE)
  ZL(:,JK) =  CST%XKARMAN/SQRT(CSTURB%XALPSBL)/CSTURB%XCMFS                                      &
              * ZZZ(:,JK)*PDIRCOSZW(:)/(ZPHIM(:,JK)**2*SQRT(ZPHIE(:,JK)))
  !$mnh_end_expand_array(JIJ=IIJB:IIJE)
END DO
!
!$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
PLK(:,:)=(1.-ZGAM(:,:))*ZL(:,:) &
                             +ZGAM(:,:)*PLK(:,:)
!$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
!
PLK(:,IKA) = PLK(:,IKB)
PLK(:,IKU) = PLK(:,IKE)
!-------------------------------------------------------------------------------
!
!*     7. Modification of the dissipative length
!         --------------------------------------
!
!$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
ZL(:,:) = ZL(:,:) * (CSTURB%XALPSBL**(3./2.)*CST%XKARMAN*TURBN%XCED) &
        / (CST%XKARMAN/SQRT(CSTURB%XALPSBL)/CSTURB%XCMFS)
!$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
!
!$mnh_expand_where(JIJ=IIJB:IIJE,JK=1:IKT)
WHERE (ZZ_O_LMO(:,:)<0.)
  ZL(:,:) = ZL(:,:)/(1.-1.9*ZZ_O_LMO(:,:))
ELSEWHERE
  ZL(:,:) = ZL(:,:)/(1.-0.3*SQRT(ZZ_O_LMO(:,:)))
END WHERE
!$mnh_end_expand_where(JIJ=IIJB:IIJE,JK=1:IKT)
!
!$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
PLEPS(:,:)=(1.-ZGAM(:,:))*ZL(:,:) &
                               +ZGAM(:,:)*PLEPS(:,:)
!$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
!
PLEPS(:,IKA) = PLEPS(:,IKB)
PLEPS(:,IKU) = PLEPS(:,IKE)
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('RMC01',1,ZHOOK_HANDLE)
END SUBROUTINE RMC01
END MODULE MODE_RMC01
