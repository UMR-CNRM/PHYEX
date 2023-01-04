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
USE MODE_UPDATE_IIJU_PHY, ONLY: UPDATE_IIJU_PHY
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
REAL(KIND=JPRB) :: ZHOOK_HANDLE
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
  ZZZ(IIJB:IIJE,JK) = ZZZ(IIJB:IIJE,JK) - PZZ(IIJB:IIJE,IKB)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE)
END DO
! fill upper level with physical value
!$mnh_expand_array(JIJ=IIJB:IIJE)
ZZZ(IIJB:IIJE,IKU) = 2.*ZZZ(IIJB:IIJE,IKU-IKL) - ZZZ(IIJB:IIJE,IKU-2*IKL)
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
  WHERE (PLMO(IIJB:IIJE)==XUNDEF)
    ZZ_O_LMO(IIJB:IIJE,JK)=0.
  ELSEWHERE
    ZZ_O_LMO(IIJB:IIJE,JK)=ZZZ(IIJB:IIJE,JK)*PDIRCOSZW(IIJB:IIJE)/PLMO(IIJB:IIJE)
  END WHERE
  !$mnh_end_expand_where(JIJ=IIJB:IIJE)
END DO
!$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
ZZ_O_LMO(IIJB:IIJE,1:IKT) = MAX(ZZ_O_LMO(IIJB:IIJE,1:IKT),-10.)
ZZ_O_LMO(IIJB:IIJE,1:IKT) = MIN(ZZ_O_LMO(IIJB:IIJE,1:IKT), 10.)
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
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ZDH(IIJB:IIJE,1:IKT) = SQRT(ZWORK1(IIJB:IIJE,1:IKT)*ZWORK2(IIJB:IIJE,1:IKT))
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    !
    CALL UPDATE_IIJU_PHY(D,ZZC)
    !
    DO JK=1,IKT
      !$mnh_expand_array(JIJ=IIJB:IIJE)
      ZZC(IIJB:IIJE,JK) = 2.*MIN(ZPHIM(IIJB:IIJE,JK),1.)/CST%XKARMAN    &
                            * MAX( PDZZ(IIJB:IIJE,JK)*PDIRCOSZW(IIJB:IIJE) , & 
                            ZDH(IIJB:IIJE,JK)/PDIRCOSZW(IIJB:IIJE)/3. )
      !$mnh_end_expand_array(JIJ=IIJB:IIJE)
    END DO
    !
!*     4. factor controling the transition between SBL and free isotropic turb. (3D case)
!         --------------------------------------------------------------------
!
    ZGAM(IIJB:IIJE,IKA) = 0.
    DO JK=IKTB,IKTE
      !$mnh_expand_array(JIJ=IIJB:IIJE)
      ZGAM(IIJB:IIJE,JK) = 1.  - EXP( -3.*(ZZZ(IIJB:IIJE,JK)-ZZZ(IIJB:IIJE,IKB))/(ZZC(IIJB:IIJE,JK)) )
      !$mnh_end_expand_array(JIJ=IIJB:IIJE)
      !$mnh_expand_where(JIJ=IIJB:IIJE)
      WHERE (ZGAM(IIJB:IIJE,JK-IKL)>ZGAM(IIJB:IIJE,JK) .OR. ZGAM(IIJB:IIJE,JK-IKL)>0.99 ) 
        ZGAM(IIJB:IIJE,JK) = 1.
      END WHERE
     !$mnh_end_expand_where(JIJ=IIJB:IIJE)
    END DO
    !$mnh_expand_array(JIJ=IIJB:IIJE)
    ZGAM(IIJB:IIJE,IKU) = 1.  - EXP( -3.*(ZZZ(IIJB:IIJE,IKU)-ZZZ(IIJB:IIJE,IKB))& 
                                   /(ZZC(IIJB:IIJE,IKU)) )
    !$mnh_end_expand_array(JIJ=IIJB:IIJE)
    !$mnh_expand_where(JIJ=IIJB:IIJE)
    WHERE (ZGAM(IIJB:IIJE,IKU-IKL)>ZGAM(IIJB:IIJE,IKU) .OR. ZGAM(IIJB:IIJE,IKU-IKL)>0.99 ) 
      ZGAM(IIJB:IIJE,IKU) = 1.
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
    ZGAM(IIJB:IIJE,1:IKT) = 1.
    ZGAM(IIJB:IIJE,IKA) = 0.
    DO JK=IKTB,IKTE
      !$mnh_expand_where(JIJ=IIJB:IIJE)
      WHERE(PSBL_DEPTH(IIJB:IIJE)>0.)
        ZGAM(IIJB:IIJE,JK) = TANH( (ZZZ(IIJB:IIJE,JK)-ZZZ(IIJB:IIJE,IKB))/PSBL_DEPTH(IIJB:IIJE) )
      END WHERE
      !$mnh_end_expand_where(JIJ=IIJB:IIJE)
      !$mnh_expand_where(JIJ=IIJB:IIJE)
      WHERE (ZGAM(IIJB:IIJE,JK-IKL)>0.99 ) 
        ZGAM(IIJB:IIJE,JK) = 1.
      END WHERE
      !$mnh_end_expand_where(JIJ=IIJB:IIJE)
    END DO
    !$mnh_expand_where(JIJ=IIJB:IIJE)
    WHERE(PSBL_DEPTH(IIJB:IIJE)>0.)
      ZGAM(IIJB:IIJE,IKU) = TANH( (ZZZ(IIJB:IIJE,IKU)-ZZZ(IIJB:IIJE,IKB))/PSBL_DEPTH(IIJB:IIJE) )
    END WHERE
   !$mnh_end_expand_where(JIJ=IIJB:IIJE)
   !$mnh_expand_where(JIJ=IIJB:IIJE)
    WHERE (ZGAM(IIJB:IIJE,IKU-IKL)>0.99 ) 
      ZGAM(IIJB:IIJE,JK) = 1.
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
  ZL(IIJB:IIJE,JK) =  CST%XKARMAN/SQRT(CSTURB%XALPSBL)/CSTURB%XCMFS                                      &
              * ZZZ(IIJB:IIJE,JK)*PDIRCOSZW(IIJB:IIJE)/(ZPHIM(IIJB:IIJE,JK)**2*SQRT(ZPHIE(IIJB:IIJE,JK)))
!$mnh_end_expand_array(JIJ=IIJB:IIJE)
END DO
!
!$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
PLK(IIJB:IIJE,1:IKT)=(1.-ZGAM(IIJB:IIJE,1:IKT))*ZL(IIJB:IIJE,1:IKT) &
                             +ZGAM(IIJB:IIJE,1:IKT)*PLK(IIJB:IIJE,1:IKT)
!$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
!
PLK(IIJB:IIJE,IKA) = PLK(IIJB:IIJE,IKB)
PLK(IIJB:IIJE,IKU) = PLK(IIJB:IIJE,IKE)
!-------------------------------------------------------------------------------
!
!*     7. Modification of the dissipative length
!         --------------------------------------
!
!$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
ZL(IIJB:IIJE,1:IKT) = ZL(IIJB:IIJE,1:IKT) * (CSTURB%XALPSBL**(3./2.)*CST%XKARMAN*CSTURB%XCED) &
        / (CST%XKARMAN/SQRT(CSTURB%XALPSBL)/CSTURB%XCMFS)
!$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
!
!$mnh_expand_where(JIJ=IIJB:IIJE,JK=1:IKT)
WHERE (ZZ_O_LMO(IIJB:IIJE,1:IKT)<0.)
  ZL(IIJB:IIJE,1:IKT) = ZL(IIJB:IIJE,1:IKT)/(1.-1.9*ZZ_O_LMO(IIJB:IIJE,1:IKT))
ELSEWHERE
  ZL(IIJB:IIJE,1:IKT) = ZL(IIJB:IIJE,1:IKT)/(1.-0.3*SQRT(ZZ_O_LMO(IIJB:IIJE,1:IKT)))
END WHERE
!$mnh_end_expand_where(JIJ=IIJB:IIJE,JK=1:IKT)
!
!$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
PLEPS(IIJB:IIJE,1:IKT)=(1.-ZGAM(IIJB:IIJE,1:IKT))*ZL(IIJB:IIJE,1:IKT) &
                               +ZGAM(IIJB:IIJE,1:IKT)*PLEPS(IIJB:IIJE,1:IKT)
!$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
!
PLEPS(IIJB:IIJE,IKA) = PLEPS(IIJB:IIJE,IKB)
PLEPS(IIJB:IIJE,IKU) = PLEPS(IIJB:IIJE,IKE)
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('RMC01',1,ZHOOK_HANDLE)
END SUBROUTINE RMC01
END MODULE MODE_RMC01
