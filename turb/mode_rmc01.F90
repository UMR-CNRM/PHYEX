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
  DO JIJ=IIJB, IIJE
    ZZZ(JIJ, JK) = ZZZ(JIJ, JK) - PZZ(JIJ, IKB)
  END DO
END DO
! fill upper level with physical value
DO JIJ=IIJB, IIJE
  ZZZ(JIJ, IKU) = 2.*ZZZ(JIJ, IKU-IKL) - ZZZ(JIJ, IKU-2*IKL)
END DO
!
!-------------------------------------------------------------------------------
!
!*     2. MO quantities
!         -------------
!
! z/LMO
DO JK=1,IKT
  DO JIJ=IIJB, IIJE
    IF (PLMO(JIJ)==XUNDEF) THEN
      ZZ_O_LMO(JIJ, JK)=0.
    ELSE
      ZZ_O_LMO(JIJ, JK)=ZZZ(JIJ, JK)*PDIRCOSZW(JIJ)/PLMO(JIJ)
    END IF
  END DO
END DO
DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZZ_O_LMO(JIJ, JK) = MAX(ZZ_O_LMO(JIJ, JK),-10.)
    ZZ_O_LMO(JIJ, JK) = MIN(ZZ_O_LMO(JIJ, JK), 10.)
  END DO
END DO
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
    DO JK=1, IKT
      DO JIJ=IIJB, IIJE
        ZDH(JIJ, JK) = SQRT(ZWORK1(JIJ, JK)*ZWORK2(JIJ, JK))
      END DO
    END DO
    !
    CALL UPDATE_IIJU_PHY(D,ZZC)
    !
    DO JK=1,IKT
      DO JIJ=IIJB, IIJE
        ZZC(JIJ, JK) = 2.*MIN(ZPHIM(JIJ, JK),1.)/CST%XKARMAN    &
                              * MAX( PDZZ(JIJ, JK)*PDIRCOSZW(JIJ) , & 
                              ZDH(JIJ, JK)/PDIRCOSZW(JIJ)/3. )
      END DO
    END DO
    !
    !*     4. factor controling the transition between SBL and free isotropic turb. (3D case)
    !         --------------------------------------------------------------------
    !
    ZGAM(IIJB:IIJE,IKA) = 0.
    DO JK=IKTB,IKTE
      DO JIJ=IIJB, IIJE
        ZGAM(JIJ, JK) = 1.  - EXP( -3.*(ZZZ(JIJ, JK)-ZZZ(JIJ, IKB))/(ZZC(JIJ, JK)) )
      END DO
      DO JIJ=IIJB, IIJE
        IF (ZGAM(JIJ, JK-IKL)>ZGAM(JIJ, JK) .OR. ZGAM(JIJ, JK-IKL)>0.99 ) THEN 
          ZGAM(JIJ, JK) = 1.
        END IF
     END DO
    END DO
    DO JIJ=IIJB, IIJE
      ZGAM(JIJ, IKU) = 1.  - EXP( -3.*(ZZZ(JIJ, IKU)-ZZZ(JIJ, IKB))& 
                                     /(ZZC(JIJ, IKU)) )
    END DO
    DO JIJ=IIJB, IIJE
      IF (ZGAM(JIJ, IKU-IKL)>ZGAM(JIJ, IKU) .OR. ZGAM(JIJ, IKU-IKL)>0.99 ) THEN 
        ZGAM(JIJ, IKU) = 1.
      END IF
    END DO
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
      DO JIJ=IIJB, IIJE
        IF (PSBL_DEPTH(JIJ)>0.) THEN
          ZGAM(JIJ, JK) = TANH( (ZZZ(JIJ, JK)-ZZZ(JIJ, IKB))/PSBL_DEPTH(JIJ) )
        END IF
      END DO
      DO JIJ=IIJB, IIJE
        IF (ZGAM(JIJ, JK-IKL)>0.99 ) THEN 
          ZGAM(JIJ, JK) = 1.
        END IF
      END DO
    END DO
    DO JIJ=IIJB, IIJE
      IF (PSBL_DEPTH(JIJ)>0.) THEN
        ZGAM(JIJ, IKU) = TANH( (ZZZ(JIJ, IKU)-ZZZ(JIJ, IKB))/PSBL_DEPTH(JIJ) )
      END IF
    END DO
    DO JIJ=IIJB, IIJE
      IF (ZGAM(JIJ, IKU-IKL)>0.99 ) THEN 
        ZGAM(JIJ, JK) = 1.
      END IF
    END DO
  !
  !-------------------------------------------------------------------------------
END SELECT
!-------------------------------------------------------------------------------
!
!*     6. Modification of the mixing length
!         ---------------------------------
!
DO JK=1,IKT
  DO JIJ=IIJB, IIJE
    ZL(JIJ, JK) =  CST%XKARMAN/SQRT(CSTURB%XALPSBL)/CSTURB%XCMFS                                      &
                * ZZZ(JIJ, JK)*PDIRCOSZW(JIJ)/(ZPHIM(JIJ, JK)**2*SQRT(ZPHIE(JIJ, JK)))
  END DO
END DO
!
DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    PLK(JIJ, JK)=(1.-ZGAM(JIJ, JK))*ZL(JIJ, JK) &
                                 +ZGAM(JIJ, JK)*PLK(JIJ, JK)
  END DO
END DO
!
PLK(IIJB:IIJE,IKA) = PLK(IIJB:IIJE,IKB)
PLK(IIJB:IIJE,IKU) = PLK(IIJB:IIJE,IKE)
!-------------------------------------------------------------------------------
!
!*     7. Modification of the dissipative length
!         --------------------------------------
!
DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZL(JIJ, JK) = ZL(JIJ, JK) * (CSTURB%XALPSBL**(3./2.)*CST%XKARMAN*TURBN%XCED) &
            / (CST%XKARMAN/SQRT(CSTURB%XALPSBL)/CSTURB%XCMFS)
  END DO
END DO
!
DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    IF (ZZ_O_LMO(JIJ, JK)<0.) THEN
      ZL(JIJ, JK) = ZL(JIJ, JK)/(1.-1.9*ZZ_O_LMO(JIJ, JK))
    ELSE
      ZL(JIJ, JK) = ZL(JIJ, JK)/(1.-0.3*SQRT(ZZ_O_LMO(JIJ, JK)))
    END IF
  END DO
END DO
!
DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    PLEPS(JIJ, JK)=(1.-ZGAM(JIJ, JK))*ZL(JIJ, JK) &
                                   +ZGAM(JIJ, JK)*PLEPS(JIJ, JK)
  END DO
END DO
!
PLEPS(IIJB:IIJE,IKA) = PLEPS(IIJB:IIJE,IKB)
PLEPS(IIJB:IIJE,IKU) = PLEPS(IIJB:IIJE,IKE)
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('RMC01',1,ZHOOK_HANDLE)
END SUBROUTINE RMC01
END MODULE MODE_RMC01
