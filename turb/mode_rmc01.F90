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
USE MODE_RMC01_3D, ONLY: RMC01_3D
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
!
! altitude of mass points
CALL MZF_PHY(D,PZZ,ZZZ)
! replace by height of mass points
DO JK=1,D%NKT
  DO JIJ=IIJB,IIJE 
    ZZZ(JIJ,JK) = ZZZ(JIJ,JK) - PZZ(JIJ,IKB)
  ENDDO
END DO
! fill upper level with physical value
DO JIJ=IIJB,IIJE 
  ZZZ(JIJ,D%NKU) = 2.*ZZZ(JIJ,D%NKU-D%NKL) - ZZZ(JIJ,D%NKU-2*D%NKL)
ENDDO
!
!-------------------------------------------------------------------------------
!
!*     2. MO quantities
!         -------------
!
! z/LMO
DO JK=1,D%NKT
  DO JIJ=IIJB,IIJE 
    IF (PLMO(JIJ)==XUNDEF)THEN
      ZZ_O_LMO(JIJ,JK)=0.
    ELSE
      ZZ_O_LMO(JIJ,JK)=ZZZ(JIJ,JK)*PDIRCOSZW(JIJ)/PLMO(JIJ)
    ENDIF
  ENDDO
END DO
DO JK=1,D%NKT 
  DO JIJ=IIJB,IIJE 
    ZZ_O_LMO(JIJ,JK) = MAX(ZZ_O_LMO(JIJ,JK),-10.)
    ZZ_O_LMO(JIJ,JK) = MIN(ZZ_O_LMO(JIJ,JK), 10.)
  ENDDO
ENDDO
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
    CALL RMC01_3D(D,CST,PDXX,PDYY,PDZZ,PDIRCOSZW,ZPHIM,ZZC)
    !
!*     4. factor controling the transition between SBL and free isotropic turb. (3D case)
!         --------------------------------------------------------------------
!
    ZGAM(IIJB:IIJE,D%NKA) = 0.
    DO JK=IKTB,IKTE
  DO JIJ=IIJB,IIJE 
    ZGAM(JIJ,JK) = 1.  - EXP( -3.*(ZZZ(JIJ,JK)-ZZZ(JIJ,IKB))/(ZZC(JIJ,JK)) )
  ENDDO
  DO JIJ=IIJB,IIJE 
    IF (ZGAM(JIJ,JK-D%NKL)>ZGAM(JIJ,JK) .OR. ZGAM(JIJ,JK-D%NKL)>0.99 ) THEN
      ZGAM(JIJ,JK) = 1.
    ENDIF
  ENDDO
    END DO
DO JIJ=IIJB,IIJE 
  ZGAM(JIJ,D%NKU) = 1.  - EXP( -3.*(ZZZ(JIJ,D%NKU)-ZZZ(JIJ,IKB))& 
  /(ZZC(JIJ,D%NKU)) )
ENDDO
DO JIJ=IIJB,IIJE 
  IF (ZGAM(JIJ,D%NKU-D%NKL)>ZGAM(JIJ,D%NKU) .OR. ZGAM(JIJ,D%NKU-D%NKL)>0.99 ) THEN
    ZGAM(JIJ,D%NKU) = 1.
  ENDIF
ENDDO
!   
!
!-------------------------------------------------------------------------------
!
!*     5. factor controling the transition between SBL and free isotropic turb.(1D case)
!         --------------------------------------------------------------------
!
  CASE DEFAULT
!* SBL depth is used
    ZGAM(IIJB:IIJE,1:D%NKT) = 1.
    ZGAM(IIJB:IIJE,D%NKA) = 0.
    DO JK=IKTB,IKTE
  DO JIJ=IIJB,IIJE 
    IF(PSBL_DEPTH(JIJ)>0.)THEN
      ZGAM(JIJ,JK) = TANH( (ZZZ(JIJ,JK)-ZZZ(JIJ,IKB))/PSBL_DEPTH(JIJ) )
    ENDIF
  ENDDO
  DO JIJ=IIJB,IIJE 
    IF (ZGAM(JIJ,JK-D%NKL)>0.99 ) THEN
      ZGAM(JIJ,JK) = 1.
    ENDIF
  ENDDO
    END DO
DO JIJ=IIJB,IIJE 
  IF(PSBL_DEPTH(JIJ)>0.)THEN
    ZGAM(JIJ,D%NKU) = TANH( (ZZZ(JIJ,D%NKU)-ZZZ(JIJ,IKB))/PSBL_DEPTH(JIJ) )
  ENDIF
ENDDO
DO JIJ=IIJB,IIJE 
  IF (ZGAM(JIJ,D%NKU-D%NKL)>0.99 ) THEN
    ZGAM(JIJ,JK) = 1.
  ENDIF
ENDDO
!
!-------------------------------------------------------------------------------
END SELECT
!-------------------------------------------------------------------------------
!
!*     6. Modification of the mixing length
!         ---------------------------------
!
DO JK=1,D%NKT
  DO JIJ=IIJB,IIJE 
    ZL(JIJ,JK) =  CST%XKARMAN/SQRT(CSTURB%XALPSBL)/CSTURB%XCMFS                                      &
    * ZZZ(JIJ,JK)*PDIRCOSZW(JIJ)/(ZPHIM(JIJ,JK)**2*SQRT(ZPHIE(JIJ,JK)))
  ENDDO
END DO
!
DO JK=1,D%NKT 
  DO JIJ=IIJB,IIJE 
    PLK(JIJ,JK)=(1.-ZGAM(JIJ,JK))*ZL(JIJ,JK) &
    +ZGAM(JIJ,JK)*PLK(JIJ,JK)
  ENDDO
ENDDO
!
PLK(IIJB:IIJE,D%NKA) = PLK(IIJB:IIJE,IKB)
PLK(IIJB:IIJE,D%NKU) = PLK(IIJB:IIJE,IKE)
!-------------------------------------------------------------------------------
!
!*     7. Modification of the dissipative length
!         --------------------------------------
!
DO JK=1,D%NKT 
  DO JIJ=IIJB,IIJE 
    ZL(JIJ,JK) = ZL(JIJ,JK) * (CSTURB%XALPSBL**(3./2.)*CST%XKARMAN*CSTURB%XCED) &
    / (CST%XKARMAN/SQRT(CSTURB%XALPSBL)/CSTURB%XCMFS)
  ENDDO
ENDDO
!
DO JK=1,D%NKT 
  DO JIJ=IIJB,IIJE 
    IF (ZZ_O_LMO(JIJ,JK)<0.)THEN
      ZL(JIJ,JK) = ZL(JIJ,JK)/(1.-1.9*ZZ_O_LMO(JIJ,JK))
    ELSE
      ZL(JIJ,JK) = ZL(JIJ,JK)/(1.-0.3*SQRT(ZZ_O_LMO(JIJ,JK)))
    ENDIF
  ENDDO
ENDDO
!
DO JK=1,D%NKT 
  DO JIJ=IIJB,IIJE 
    PLEPS(JIJ,JK)=(1.-ZGAM(JIJ,JK))*ZL(JIJ,JK) &
    +ZGAM(JIJ,JK)*PLEPS(JIJ,JK)
  ENDDO
ENDDO
!
PLEPS(IIJB:IIJE,D%NKA) = PLEPS(IIJB:IIJE,IKB)
PLEPS(IIJB:IIJE,D%NKU) = PLEPS(IIJB:IIJE,IKE)
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('RMC01',1,ZHOOK_HANDLE)
END SUBROUTINE RMC01
END MODULE MODE_RMC01
