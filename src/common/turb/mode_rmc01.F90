!MNH_LIC Copyright 1994-2022 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
MODULE MODE_RMC01
IMPLICIT NONE
CONTAINS
SUBROUTINE RMC01(D,CST,CSTURB,HTURBLEN,PZZ,PDXX,PDYY,PDZZ,PDIRCOSZW,PSBL_DEPTH,PLMO,PLK,PLEPS)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
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
!
USE MODE_SBL_PHY, ONLY: BUSINGER_PHIM, BUSINGER_PHIE
!
USE SHUMAN_PHY, ONLY: MZF_PHY, MYF_PHY, MXF_PHY
!
IMPLICIT NONE
!
!*       0.1   Declaration of arguments
!              ------------------------
!
TYPE(DIMPHYEX_t),         INTENT(IN)   :: D
TYPE(CST_t),              INTENT(IN)   :: CST
TYPE(CSTURB_t),           INTENT(IN)   :: CSTURB
CHARACTER(LEN=4),         INTENT(IN)   :: HTURBLEN ! type of mixing length
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
INTEGER :: IKTB,IKTE  ! start, end of k loops in physical domain
INTEGER :: JK,JIJ   ! loop counter
INTEGER :: IIE,IIB,IJE,IJB
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
REAL, DIMENSION(D%NIJT,D%NKT) :: ZDH  ! hor. grid mesh
                                                             ! size
REAL, DIMENSION(D%NIJT,D%NKT) :: ZL   ! SBL length
REAL, DIMENSION(D%NIJT,D%NKT) :: ZZC  ! alt. where turb. is isotr.
REAL, DIMENSION(D%NIJT,D%NKT) :: ZWORK1, ZWORK2
!-------------------------------------------------------------------------------
!
!*     1. Initializations
!         ---------------
!
! horizontal boundaries
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('RMC01',0,ZHOOK_HANDLE)
IKTB=D%NKTB          
IKTE=D%NKTE
IKB=D%NKB
IKE=D%NKE
IIE=D%NIEC
IIB=D%NIBC
IJE=D%NJEC
IJB=D%NJBC
!
! altitude of mass points
CALL MZF_PHY(D,PZZ,ZZZ)
! replace by height of mass points
DO JK=1,D%NKT
  !$mnh_expand_array(JIJ=D%NIJB:D%NIJE)
  ZZZ(D%NIJB:D%NIJE,JK) = ZZZ(D%NIJB:D%NIJE,JK) - PZZ(D%NIJB:D%NIJE,IKB)
  !$mnh_end_expand_array(JIJ=D%NIJB:D%NIJE)
END DO
! fill upper level with physical value
!$mnh_expand_array(JIJ=D%NIJB:D%NIJE)
ZZZ(D%NIJB:D%NIJE,D%NKU) = 2.*ZZZ(D%NIJB:D%NIJE,D%NKU-D%NKL) - ZZZ(D%NIJB:D%NIJE,D%NKU-2*D%NKL)
!$mnh_end_expand_array(JIJ=D%NIJB:D%NIJE)
!
!-------------------------------------------------------------------------------
!
!*     2. MO quantities
!         -------------
!
! z/LMO
DO JK=1,D%NKT
  !$mnh_expand_where(JIJ=D%NIJB:D%NIJE)
  WHERE (PLMO(D%NIJB:D%NIJE)==XUNDEF)
    ZZ_O_LMO(D%NIJB:D%NIJE,JK)=0.
  ELSEWHERE
    ZZ_O_LMO(D%NIJB:D%NIJE,JK)=ZZZ(D%NIJB:D%NIJE,JK)*PDIRCOSZW(D%NIJB:D%NIJE)/PLMO(D%NIJB:D%NIJE)
  END WHERE
  !$mnh_end_expand_where(JIJ=D%NIJB:D%NIJE)
END DO
!$mnh_expand_array(JIJ=D%NIJB:D%NIJE,JK=1:D%NKT)
ZZ_O_LMO(D%NIJB:D%NIJE,1:D%NKT) = MAX(ZZ_O_LMO(D%NIJB:D%NIJE,1:D%NKT),-10.)
ZZ_O_LMO(D%NIJB:D%NIJE,1:D%NKT) = MIN(ZZ_O_LMO(D%NIJB:D%NIJE,1:D%NKT), 10.)
!$mnh_end_expand_array(JIJ=D%NIJB:D%NIJE,JK=1:D%NKT)
!
!
! MO function for stress
CALL BUSINGER_PHIM(D,ZZ_O_LMO,ZPHIM)
!
! MO function for TKE
CALL BUSINGER_PHIE(D,CSTURB,ZZ_O_LMO,ZPHIE)
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
    CALL MXF_PHY(D,PDXX,ZWORK1)
    CALL MYF_PHY(D,PDYY,ZWORK2)
    !$mnh_expand_array(JIJ=D%NIJB:D%NIJE,JK=1:D%NKT)
    ZDH(D%NIJB:D%NIJE,1:D%NKT) = SQRT(ZWORK1(D%NIJB:D%NIJE,1:D%NKT)*ZWORK2(D%NIJB:D%NIJE,1:D%NKT))
    !$mnh_end_expand_array(JIJ=D%NIJB:D%NIJE,JK=1:D%NKT)
    ZDH(D%NIT*IJB:D%NIT*IJE:D%NIT,1:D%NKT) = ZDH(D%NIT*IJB-1:D%NIT*IJE-1:D%NIT,1:D%NKT)
    ZDH(D%NIJT-IIE+IIB:D%NIJT,1:D%NKT) = ZDH(D%NIJT-2*IIE+IIB:D%NIJT-IIE,1:D%NKT)
    DO JK=1,D%NKT
      !$mnh_expand_array(JIJ=D%NIJB:D%NIJE)
      ZZC(D%NIJB:D%NIJE,JK) = 2.*MIN(ZPHIM(D%NIJB:D%NIJE,JK),1.)/CST%XKARMAN    &
                     * MAX( PDZZ(D%NIJB:D%NIJE,JK)*PDIRCOSZW(D%NIJB:D%NIJE) , & 
                     ZDH(D%NIJB:D%NIJE,JK)/PDIRCOSZW(D%NIJB:D%NIJE)/3. )
      !$mnh_end_expand_array(JIJ=D%NIJB:D%NIJE)
    END DO
!
!*     4. factor controling the transition between SBL and free isotropic turb. (3D case)
!         --------------------------------------------------------------------
!
    ZGAM(D%NIJB:D%NIJE,D%NKA) = 0.
    DO JK=IKTB,IKTE
      !$mnh_expand_array(JIJ=D%NIJB:D%NIJE)
      ZGAM(D%NIJB:D%NIJE,JK) = 1.  - EXP( -3.*(ZZZ(D%NIJB:D%NIJE,JK)-ZZZ(D%NIJB:D%NIJE,IKB))/(ZZC(D%NIJB:D%NIJE,JK)) )
      !$mnh_end_expand_array(JIJ=D%NIJB:D%NIJE)
      !$mnh_expand_where(JIJ=D%NIJB:D%NIJE)
      WHERE (ZGAM(D%NIJB:D%NIJE,JK-D%NKL)>ZGAM(D%NIJB:D%NIJE,JK) .OR. ZGAM(D%NIJB:D%NIJE,JK-D%NKL)>0.99 ) 
        ZGAM(D%NIJB:D%NIJE,JK) = 1.
      END WHERE
     !$mnh_end_expand_where(JIJ=D%NIJB:D%NIJE)
    END DO
    !$mnh_expand_array(JIJ=D%NIJB:D%NIJE)
    ZGAM(D%NIJB:D%NIJE,D%NKU) = 1.  - EXP( -3.*(ZZZ(D%NIJB:D%NIJE,D%NKU)-ZZZ(D%NIJB:D%NIJE,IKB))& 
                                   /(ZZC(D%NIJB:D%NIJE,D%NKU)) )
    !$mnh_end_expand_array(JIJ=D%NIJB:D%NIJE)
    !$mnh_expand_where(JIJ=D%NIJB:D%NIJE)
    WHERE (ZGAM(D%NIJB:D%NIJE,D%NKU-D%NKL)>ZGAM(D%NIJB:D%NIJE,D%NKU) .OR. ZGAM(D%NIJB:D%NIJE,D%NKU-D%NKL)>0.99 ) 
      ZGAM(D%NIJB:D%NIJE,D%NKU) = 1.
    END WHERE
    !$mnh_end_expand_where(JIJ=D%NIJB:D%NIJE)
!   
!
!-------------------------------------------------------------------------------
!
!*     5. factor controling the transition between SBL and free isotropic turb.(1D case)
!         --------------------------------------------------------------------
!
  CASE DEFAULT
!* SBL depth is used
    ZGAM(D%NIJB:D%NIJE,1:D%NKT) = 1.
    ZGAM(D%NIJB:D%NIJE,D%NKA) = 0.
    DO JK=IKTB,IKTE
      !$mnh_expand_where(JIJ=D%NIJB:D%NIJE)
      WHERE(PSBL_DEPTH(D%NIJB:D%NIJE)>0.)
        ZGAM(D%NIJB:D%NIJE,JK) = TANH( (ZZZ(D%NIJB:D%NIJE,JK)-ZZZ(D%NIJB:D%NIJE,IKB))/PSBL_DEPTH(D%NIJB:D%NIJE) )
      END WHERE
      !$mnh_end_expand_where(JIJ=D%NIJB:D%NIJE)
      !$mnh_expand_where(JIJ=D%NIJB:D%NIJE)
      WHERE (ZGAM(D%NIJB:D%NIJE,JK-D%NKL)>0.99 ) 
        ZGAM(D%NIJB:D%NIJE,JK) = 1.
      END WHERE
      !$mnh_end_expand_where(JIJ=D%NIJB:D%NIJE)
    END DO
    !$mnh_expand_where(JIJ=D%NIJB:D%NIJE)
    WHERE(PSBL_DEPTH(D%NIJB:D%NIJE)>0.)
      ZGAM(D%NIJB:D%NIJE,D%NKU) = TANH( (ZZZ(D%NIJB:D%NIJE,D%NKU)-ZZZ(D%NIJB:D%NIJE,IKB))/PSBL_DEPTH(D%NIJB:D%NIJE) )
    END WHERE
   !$mnh_end_expand_where(JIJ=D%NIJB:D%NIJE)
   !$mnh_expand_where(JIJ=D%NIJB:D%NIJE)
    WHERE (ZGAM(D%NIJB:D%NIJE,D%NKU-D%NKL)>0.99 ) 
      ZGAM(D%NIJB:D%NIJE,JK) = 1.
    END WHERE
    !$mnh_end_expand_where(JIJ=D%NIJB:D%NIJE)
!
!-------------------------------------------------------------------------------
END SELECT
!-------------------------------------------------------------------------------
!
!*     6. Modification of the mixing length
!         ---------------------------------
!
DO JK=1,D%NKT
!$mnh_expand_array(JIJ=D%NIJB:D%NIJE)
  ZL(D%NIJB:D%NIJE,JK) =  CST%XKARMAN/SQRT(CSTURB%XALPSBL)/CSTURB%XCMFS                                      &
              * ZZZ(D%NIJB:D%NIJE,JK)*PDIRCOSZW(D%NIJB:D%NIJE)/(ZPHIM(D%NIJB:D%NIJE,JK)**2*SQRT(ZPHIE(D%NIJB:D%NIJE,JK)))
!$mnh_end_expand_array(JIJ=D%NIJB:D%NIJE)
END DO
!
!$mnh_expand_array(JIJ=D%NIJB:D%NIJE,JK=1:D%NKT)
PLK(D%NIJB:D%NIJE,1:D%NKT)=(1.-ZGAM(D%NIJB:D%NIJE,1:D%NKT))*ZL(D%NIJB:D%NIJE,1:D%NKT) &
                             +ZGAM(D%NIJB:D%NIJE,1:D%NKT)*PLK(D%NIJB:D%NIJE,1:D%NKT)
!$mnh_end_expand_array(JIJ=D%NIJB:D%NIJE,JK=1:D%NKT)
!
PLK(D%NIJB:D%NIJE,D%NKA) = PLK(D%NIJB:D%NIJE,IKB)
PLK(D%NIJB:D%NIJE,D%NKU) = PLK(D%NIJB:D%NIJE,IKE)
!-------------------------------------------------------------------------------
!
!*     7. Modification of the dissipative length
!         --------------------------------------
!
!$mnh_expand_array(JIJ=D%NIJB:D%NIJE,JK=1:D%NKT)
ZL(D%NIJB:D%NIJE,1:D%NKT) = ZL(D%NIJB:D%NIJE,1:D%NKT) * (CSTURB%XALPSBL**(3./2.)*CST%XKARMAN*CSTURB%XCED) &
        / (CST%XKARMAN/SQRT(CSTURB%XALPSBL)/CSTURB%XCMFS)
!$mnh_end_expand_array(JIJ=D%NIJB:D%NIJE,JK=1:D%NKT)
!
!$mnh_expand_where(JIJ=D%NIJB:D%NIJE,JK=1:D%NKT)
WHERE (ZZ_O_LMO(D%NIJB:D%NIJE,1:D%NKT)<0.)
  ZL(D%NIJB:D%NIJE,1:D%NKT) = ZL(D%NIJB:D%NIJE,1:D%NKT)/(1.-1.9*ZZ_O_LMO(D%NIJB:D%NIJE,1:D%NKT))
ELSEWHERE
  ZL(D%NIJB:D%NIJE,1:D%NKT) = ZL(D%NIJB:D%NIJE,1:D%NKT)/(1.-0.3*SQRT(ZZ_O_LMO(D%NIJB:D%NIJE,1:D%NKT)))
END WHERE
!$mnh_end_expand_where(JIJ=D%NIJB:D%NIJE,JK=1:D%NKT)
!
!$mnh_expand_array(JIJ=D%NIJB:D%NIJE,JK=1:D%NKT)
PLEPS(D%NIJB:D%NIJE,1:D%NKT)=(1.-ZGAM(D%NIJB:D%NIJE,1:D%NKT))*ZL(D%NIJB:D%NIJE,1:D%NKT) &
                               +ZGAM(D%NIJB:D%NIJE,1:D%NKT)*PLEPS(D%NIJB:D%NIJE,1:D%NKT)
!$mnh_end_expand_array(JIJ=D%NIJB:D%NIJE,JK=1:D%NKT)
!
PLEPS(D%NIJB:D%NIJE,D%NKA) = PLEPS(D%NIJB:D%NIJE,IKB)
PLEPS(D%NIJB:D%NIJE,D%NKU) = PLEPS(D%NIJB:D%NIJE,IKE)
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('RMC01',1,ZHOOK_HANDLE)
END SUBROUTINE RMC01
END MODULE MODE_RMC01
