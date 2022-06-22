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
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    :: PZZ   ! altitude of flux points
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    :: PDXX  ! width of grid mesh (X dir)
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    :: PDYY  ! width of grid mesh (Y dir)
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    :: PDZZ  ! width of vert. layers
REAL, DIMENSION(D%NIT,D%NJT),     INTENT(IN)    :: PDIRCOSZW ! Director Cosinus
REAL, DIMENSION(D%NIT,D%NJT),     INTENT(IN)    :: PSBL_DEPTH! SBL depth
REAL, DIMENSION(D%NIT,D%NJT),     INTENT(IN)    :: PLMO  ! Monin Obuhkov length
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(INOUT) :: PLK   ! Mixing length
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(INOUT) :: PLEPS ! Dissipative length
!
!*       0.2   Declaration of local variables
!              ------------------------------
!
INTEGER :: IKB,IKE    ! first,last physical level
INTEGER :: IKTB,IKTE  ! start, end of k loops in physical domain
INTEGER :: IIU        ! horizontal x boundary
INTEGER :: IJU        ! horizontal y boundary
INTEGER :: JK,JI,JJ   ! loop counter
INTEGER :: IIE,IIB,IJE,IJB
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: ZZZ  ! height of mass
                                                             ! points above ground
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: ZZ_O_LMO ! height / LMO
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: ZGAM ! factor controling
                                                             ! transition betw.
                                                             ! SBL and free BL

REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: ZPHIM! MO function
                                                             ! for stress
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: ZPHIE! MO function
                                                             ! for TKE
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: ZDH  ! hor. grid mesh
                                                             ! size
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: ZL   ! SBL length
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: ZZC  ! alt. where turb. is isotr.
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: ZWORK1, ZWORK2
!-------------------------------------------------------------------------------
!
!*     1. Initializations
!         ---------------
!
! horizontal boundaries
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('RMC01',0,ZHOOK_HANDLE)
IIU=D%NIT
IJU=D%NJT
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
  !$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE)
  ZZZ(IIB:IIE,IJB:IJE,JK) = ZZZ(IIB:IIE,IJB:IJE,JK) - PZZ(IIB:IIE,IJB:IJE,IKB)
  !$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE)
END DO
! fill upper level with physical value
!$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE)
ZZZ(IIB:IIE,IJB:IJE,D%NKU) = 2.*ZZZ(IIB:IIE,IJB:IJE,D%NKU-D%NKL) - ZZZ(IIB:IIE,IJB:IJE,D%NKU-2*D%NKL)
!$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE)
!
!-------------------------------------------------------------------------------
!
!*     2. MO quantities
!         -------------
!
! z/LMO
DO JK=1,D%NKT
  !$mnh_expand_where(JI=IIB:IIE,JJ=IJB:IJE)
  WHERE (PLMO(IIB:IIE,IJB:IJE)==XUNDEF)
    ZZ_O_LMO(IIB:IIE,IJB:IJE,JK)=0.
  ELSEWHERE
    ZZ_O_LMO(IIB:IIE,IJB:IJE,JK)=ZZZ(IIB:IIE,IJB:IJE,JK)*PDIRCOSZW(IIB:IIE,IJB:IJE)/PLMO(IIB:IIE,IJB:IJE)
  END WHERE
  !$mnh_end_expand_where(JI=IIB:IIE,JJ=IJB:IJE)
END DO
!$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)
ZZ_O_LMO(IIB:IIE,IJB:IJE,1:D%NKT) = MAX(ZZ_O_LMO(IIB:IIE,IJB:IJE,1:D%NKT),-10.)
ZZ_O_LMO(IIB:IIE,IJB:IJE,1:D%NKT) = MIN(ZZ_O_LMO(IIB:IIE,IJB:IJE,1:D%NKT), 10.)
!$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)
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
    !$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)
    ZDH(IIB:IIE,IJB:IJE,1:D%NKT) = SQRT(ZWORK1(IIB:IIE,IJB:IJE,1:D%NKT)*ZWORK2(IIB:IIE,IJB:IJE,1:D%NKT))
    !$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)
    ZDH(IIU,IJB:IJE,1:D%NKT) = ZDH(IIU-1,IJB:IJE,1:D%NKT)
    ZDH(IIB:IIE,IJU,1:D%NKT) = ZDH(IIB:IIE,IJU-1,1:D%NKT)
    DO JK=1,D%NKT
      !$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE)
      ZZC(IIB:IIE,IJB:IJE,JK) = 2.*MIN(ZPHIM(IIB:IIE,IJB:IJE,JK),1.)/CST%XKARMAN    &
                     * MAX( PDZZ(IIB:IIE,IJB:IJE,JK)*PDIRCOSZW(IIB:IIE,IJB:IJE) , & 
                     ZDH(IIB:IIE,IJB:IJE,JK)/PDIRCOSZW(IIB:IIE,IJB:IJE)/3. )
      !$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE)
    END DO
!
!*     4. factor controling the transition between SBL and free isotropic turb. (3D case)
!         --------------------------------------------------------------------
!
    ZGAM(IIB:IIE,IJB:IJE,D%NKA) = 0.
    DO JK=IKTB,IKTE
      !$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE)
      ZGAM(IIB:IIE,IJB:IJE,JK) = 1.  - EXP( -3.*(ZZZ(IIB:IIE,IJB:IJE,JK)-ZZZ(IIB:IIE,IJB:IJE,IKB))/(ZZC(IIB:IIE,IJB:IJE,JK)) )
      !$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE)
      !$mnh_expand_where(JI=IIB:IIE,JJ=IJB:IJE)
      WHERE (ZGAM(IIB:IIE,IJB:IJE,JK-D%NKL)>ZGAM(IIB:IIE,IJB:IJE,JK) .OR. ZGAM(IIB:IIE,IJB:IJE,JK-D%NKL)>0.99 ) 
        ZGAM(IIB:IIE,IJB:IJE,JK) = 1.
      END WHERE
     !$mnh_end_expand_where(JI=IIB:IIE,JJ=IJB:IJE)
    END DO
    !$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE)
    ZGAM(IIB:IIE,IJB:IJE,D%NKU) = 1.  - EXP( -3.*(ZZZ(IIB:IIE,IJB:IJE,D%NKU)-ZZZ(IIB:IIE,IJB:IJE,IKB))& 
                                   /(ZZC(IIB:IIE,IJB:IJE,D%NKU)) )
    !$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE)
    !$mnh_expand_where(JI=IIB:IIE,JJ=IJB:IJE)
    WHERE (ZGAM(IIB:IIE,IJB:IJE,D%NKU-D%NKL)>ZGAM(IIB:IIE,IJB:IJE,D%NKU) .OR. ZGAM(IIB:IIE,IJB:IJE,D%NKU-D%NKL)>0.99 ) 
      ZGAM(IIB:IIE,IJB:IJE,D%NKU) = 1.
    END WHERE
    !$mnh_end_expand_where(JI=IIB:IIE,JJ=IJB:IJE)
!   
!
!-------------------------------------------------------------------------------
!
!*     5. factor controling the transition between SBL and free isotropic turb.(1D case)
!         --------------------------------------------------------------------
!
  CASE DEFAULT
!* SBL depth is used
    ZGAM(IIB:IIE,IJB:IJE,1:D%NKT) = 1.
    ZGAM(IIB:IIE,IJB:IJE,D%NKA) = 0.
    DO JK=IKTB,IKTE
      !$mnh_expand_where(JI=IIB:IIE,JJ=IJB:IJE)
      WHERE(PSBL_DEPTH(IIB:IIE,IJB:IJE)>0.)
        ZGAM(IIB:IIE,IJB:IJE,JK) = TANH( (ZZZ(IIB:IIE,IJB:IJE,JK)-ZZZ(IIB:IIE,IJB:IJE,IKB))/PSBL_DEPTH(IIB:IIE,IJB:IJE) )
      END WHERE
      !$mnh_end_expand_where(JI=IIB:IIE,JJ=IJB:IJE)
      !$mnh_expand_where(JI=IIB:IIE,JJ=IJB:IJE)
      WHERE (ZGAM(IIB:IIE,IJB:IJE,JK-D%NKL)>0.99 ) 
        ZGAM(IIB:IIE,IJB:IJE,JK) = 1.
      END WHERE
      !$mnh_end_expand_where(JI=IIB:IIE,JJ=IJB:IJE)
    END DO
    !$mnh_expand_where(JI=IIB:IIE,JJ=IJB:IJE)
    WHERE(PSBL_DEPTH(IIB:IIE,IJB:IJE)>0.)
      ZGAM(IIB:IIE,IJB:IJE,D%NKU) = TANH( (ZZZ(IIB:IIE,IJB:IJE,D%NKU)-ZZZ(IIB:IIE,IJB:IJE,IKB))/PSBL_DEPTH(IIB:IIE,IJB:IJE) )
    END WHERE
   !$mnh_end_expand_where(JI=IIB:IIE,JJ=IJB:IJE)
   !$mnh_expand_where(JI=IIB:IIE,JJ=IJB:IJE)
    WHERE (ZGAM(IIB:IIE,IJB:IJE,D%NKU-D%NKL)>0.99 ) 
      ZGAM(IIB:IIE,IJB:IJE,JK) = 1.
    END WHERE
    !$mnh_end_expand_where(JI=IIB:IIE,JJ=IJB:IJE)
!
!-------------------------------------------------------------------------------
END SELECT
!-------------------------------------------------------------------------------
!
!*     6. Modification of the mixing length
!         ---------------------------------
!
DO JK=1,D%NKT
!$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE)
  ZL(IIB:IIE,IJB:IJE,JK) =  CST%XKARMAN/SQRT(CSTURB%XALPSBL)/CSTURB%XCMFS                                      &
              * ZZZ(IIB:IIE,IJB:IJE,JK)*PDIRCOSZW(IIB:IIE,IJB:IJE)/(ZPHIM(IIB:IIE,IJB:IJE,JK)**2*SQRT(ZPHIE(IIB:IIE,IJB:IJE,JK)))
!$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE)
END DO
!
!$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)
PLK(IIB:IIE,IJB:IJE,1:D%NKT)=(1.-ZGAM(IIB:IIE,IJB:IJE,1:D%NKT))*ZL(IIB:IIE,IJB:IJE,1:D%NKT) &
                             +ZGAM(IIB:IIE,IJB:IJE,1:D%NKT)*PLK(IIB:IIE,IJB:IJE,1:D%NKT)
!$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)
!
PLK(IIB:IIE,IJB:IJE,D%NKA) = PLK(IIB:IIE,IJB:IJE,IKB)
PLK(IIB:IIE,IJB:IJE,D%NKU) = PLK(IIB:IIE,IJB:IJE,IKE)
!-------------------------------------------------------------------------------
!
!*     7. Modification of the dissipative length
!         --------------------------------------
!
!$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)
ZL(IIB:IIE,IJB:IJE,1:D%NKT) = ZL(IIB:IIE,IJB:IJE,1:D%NKT) * (CSTURB%XALPSBL**(3./2.)*CST%XKARMAN*CSTURB%XCED) &
        / (CST%XKARMAN/SQRT(CSTURB%XALPSBL)/CSTURB%XCMFS)
!$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)
!
!$mnh_expand_where(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)
WHERE (ZZ_O_LMO(IIB:IIE,IJB:IJE,1:D%NKT)<0.)
  ZL(IIB:IIE,IJB:IJE,1:D%NKT) = ZL(IIB:IIE,IJB:IJE,1:D%NKT)/(1.-1.9*ZZ_O_LMO(IIB:IIE,IJB:IJE,1:D%NKT))
ELSEWHERE
  ZL(IIB:IIE,IJB:IJE,1:D%NKT) = ZL(IIB:IIE,IJB:IJE,1:D%NKT)/(1.-0.3*SQRT(ZZ_O_LMO(IIB:IIE,IJB:IJE,1:D%NKT)))
END WHERE
!$mnh_end_expand_where(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)
!
!$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)
PLEPS(IIB:IIE,IJB:IJE,1:D%NKT)=(1.-ZGAM(IIB:IIE,IJB:IJE,1:D%NKT))*ZL(IIB:IIE,IJB:IJE,1:D%NKT) &
                               +ZGAM(IIB:IIE,IJB:IJE,1:D%NKT)*PLEPS(IIB:IIE,IJB:IJE,1:D%NKT)
!$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)
!
PLEPS(IIB:IIE,IJB:IJE,D%NKA) = PLEPS(IIB:IIE,IJB:IJE,IKB)
PLEPS(IIB:IIE,IJB:IJE,D%NKU) = PLEPS(IIB:IIE,IJB:IJE,IKE)
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('RMC01',1,ZHOOK_HANDLE)
END SUBROUTINE RMC01
END MODULE MODE_RMC01
