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
  DO JJ=IJB,IJE 
 DO JI=IIB,IIE 
  ZZZ(JI,JJ,JK) = ZZZ(JI,JJ,JK) - PZZ(JI,JJ,IKB)
   ENDDO
ENDDO
END DO
! fill upper level with physical value
DO JJ=IJB,IJE 
 DO JI=IIB,IIE 
ZZZ(JI,JJ,D%NKU) = 2.*ZZZ(JI,JJ,D%NKU-D%NKL) - ZZZ(JI,JJ,D%NKU-2*D%NKL)
 ENDDO
ENDDO
!
!-------------------------------------------------------------------------------
!
!*     2. MO quantities
!         -------------
!
! z/LMO
DO JK=1,D%NKT
  DO JJ=IJB,IJE 
 DO JI=IIB,IIE 
  IF (PLMO(JI,JJ)==XUNDEF)THEN
    ZZ_O_LMO(JI,JJ,JK)=0.
  ELSE
    ZZ_O_LMO(JI,JJ,JK)=ZZZ(JI,JJ,JK)*PDIRCOSZW(JI,JJ)/PLMO(JI,JJ)
  ENDIF
   ENDDO
ENDDO
END DO
DO JK=1,D%NKT 
 DO JJ=IJB,IJE 
  DO JI=IIB,IIE 
ZZ_O_LMO(JI,JJ,JK) = MAX(ZZ_O_LMO(JI,JJ,JK),-10.)
ZZ_O_LMO(JI,JJ,JK) = MIN(ZZ_O_LMO(JI,JJ,JK), 10.)
  ENDDO
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
    CALL MXF_PHY(D,PDXX,ZWORK1)
    CALL MYF_PHY(D,PDYY,ZWORK2)
    DO JK=1,D%NKT 
 DO JJ=IJB,IJE 
  DO JI=IIB,IIE 
    ZDH(JI,JJ,JK) = SQRT(ZWORK1(JI,JJ,JK)*ZWORK2(JI,JJ,JK))
      ENDDO
 ENDDO
ENDDO
    ZDH(IIU,IJB:IJE,1:D%NKT) = ZDH(IIU-1,IJB:IJE,1:D%NKT)
    ZDH(IIB:IIE,IJU,1:D%NKT) = ZDH(IIB:IIE,IJU-1,1:D%NKT)
    DO JK=1,D%NKT
      DO JJ=IJB,IJE 
 DO JI=IIB,IIE 
      ZZC(JI,JJ,JK) = 2.*MIN(ZPHIM(JI,JJ,JK),1.)/CST%XKARMAN    &
                     * MAX( PDZZ(JI,JJ,JK)*PDIRCOSZW(JI,JJ) , & 
                     ZDH(JI,JJ,JK)/PDIRCOSZW(JI,JJ)/3. )
       ENDDO
ENDDO
    END DO
!
!*     4. factor controling the transition between SBL and free isotropic turb. (3D case)
!         --------------------------------------------------------------------
!
    ZGAM(IIB:IIE,IJB:IJE,D%NKA) = 0.
    DO JK=IKTB,IKTE
      DO JJ=IJB,IJE 
 DO JI=IIB,IIE 
      ZGAM(JI,JJ,JK) = 1.  - EXP( -3.*(ZZZ(JI,JJ,JK)-ZZZ(JI,JJ,IKB))/(ZZC(JI,JJ,JK)) )
       ENDDO
ENDDO
      DO JJ=IJB,IJE 
 DO JI=IIB,IIE 
      IF (ZGAM(JI,JJ,JK-D%NKL)>ZGAM(JI,JJ,JK) .OR. ZGAM(JI,JJ,JK-D%NKL)>0.99 ) THEN
        ZGAM(JI,JJ,JK) = 1.
      ENDIF
      ENDDO
ENDDO
    END DO
    DO JJ=IJB,IJE 
 DO JI=IIB,IIE 
    ZGAM(JI,JJ,D%NKU) = 1.  - EXP( -3.*(ZZZ(JI,JJ,D%NKU)-ZZZ(JI,JJ,IKB))& 
                                   /(ZZC(JI,JJ,D%NKU)) )
     ENDDO
ENDDO
    DO JJ=IJB,IJE 
 DO JI=IIB,IIE 
    IF (ZGAM(JI,JJ,D%NKU-D%NKL)>ZGAM(JI,JJ,D%NKU) .OR. ZGAM(JI,JJ,D%NKU-D%NKL)>0.99 ) THEN
      ZGAM(JI,JJ,D%NKU) = 1.
    ENDIF
     ENDDO
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
    ZGAM(IIB:IIE,IJB:IJE,1:D%NKT) = 1.
    ZGAM(IIB:IIE,IJB:IJE,D%NKA) = 0.
    DO JK=IKTB,IKTE
      DO JJ=IJB,IJE 
 DO JI=IIB,IIE 
      IF(PSBL_DEPTH(JI,JJ)>0.)THEN
        ZGAM(JI,JJ,JK) = TANH( (ZZZ(JI,JJ,JK)-ZZZ(JI,JJ,IKB))/PSBL_DEPTH(JI,JJ) )
      ENDIF
       ENDDO
ENDDO
      DO JJ=IJB,IJE 
 DO JI=IIB,IIE 
      IF (ZGAM(JI,JJ,JK-D%NKL)>0.99 ) THEN
        ZGAM(JI,JJ,JK) = 1.
      ENDIF
       ENDDO
ENDDO
    END DO
    DO JJ=IJB,IJE 
 DO JI=IIB,IIE 
    IF(PSBL_DEPTH(JI,JJ)>0.)THEN
      ZGAM(JI,JJ,D%NKU) = TANH( (ZZZ(JI,JJ,D%NKU)-ZZZ(JI,JJ,IKB))/PSBL_DEPTH(JI,JJ) )
    ENDIF
    ENDDO
ENDDO
   DO JJ=IJB,IJE 
 DO JI=IIB,IIE 
    IF (ZGAM(JI,JJ,D%NKU-D%NKL)>0.99 ) THEN
      ZGAM(JI,JJ,JK) = 1.
    ENDIF
     ENDDO
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
DO JJ=IJB,IJE 
 DO JI=IIB,IIE 
  ZL(JI,JJ,JK) =  CST%XKARMAN/SQRT(CSTURB%XALPSBL)/CSTURB%XCMFS                                      &
              * ZZZ(JI,JJ,JK)*PDIRCOSZW(JI,JJ)/(ZPHIM(JI,JJ,JK)**2*SQRT(ZPHIE(JI,JJ,JK)))
 ENDDO
ENDDO
END DO
!
DO JK=1,D%NKT 
 DO JJ=IJB,IJE 
  DO JI=IIB,IIE 
PLK(JI,JJ,JK)=(1.-ZGAM(JI,JJ,JK))*ZL(JI,JJ,JK) &
                             +ZGAM(JI,JJ,JK)*PLK(JI,JJ,JK)
  ENDDO
 ENDDO
ENDDO
!
PLK(IIB:IIE,IJB:IJE,D%NKA) = PLK(IIB:IIE,IJB:IJE,IKB)
PLK(IIB:IIE,IJB:IJE,D%NKU) = PLK(IIB:IIE,IJB:IJE,IKE)
!-------------------------------------------------------------------------------
!
!*     7. Modification of the dissipative length
!         --------------------------------------
!
DO JK=1,D%NKT 
 DO JJ=IJB,IJE 
  DO JI=IIB,IIE 
ZL(JI,JJ,JK) = ZL(JI,JJ,JK) * (CSTURB%XALPSBL**(3./2.)*CST%XKARMAN*CSTURB%XCED) &
        / (CST%XKARMAN/SQRT(CSTURB%XALPSBL)/CSTURB%XCMFS)
  ENDDO
 ENDDO
ENDDO
!
DO JK=1,D%NKT 
 DO JJ=IJB,IJE 
  DO JI=IIB,IIE 
IF (ZZ_O_LMO(JI,JJ,JK)<0.)THEN
  ZL(JI,JJ,JK) = ZL(JI,JJ,JK)/(1.-1.9*ZZ_O_LMO(JI,JJ,JK))
ELSE
  ZL(JI,JJ,JK) = ZL(JI,JJ,JK)/(1.-0.3*SQRT(ZZ_O_LMO(JI,JJ,JK)))
ENDIF
  ENDDO
 ENDDO
ENDDO
!
DO JK=1,D%NKT 
 DO JJ=IJB,IJE 
  DO JI=IIB,IIE 
PLEPS(JI,JJ,JK)=(1.-ZGAM(JI,JJ,JK))*ZL(JI,JJ,JK) &
                               +ZGAM(JI,JJ,JK)*PLEPS(JI,JJ,JK)
  ENDDO
 ENDDO
ENDDO
!
PLEPS(IIB:IIE,IJB:IJE,D%NKA) = PLEPS(IIB:IIE,IJB:IJE,IKB)
PLEPS(IIB:IIE,IJB:IJE,D%NKU) = PLEPS(IIB:IIE,IJB:IJE,IKE)
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('RMC01',1,ZHOOK_HANDLE)
END SUBROUTINE RMC01
END MODULE MODE_RMC01
