!MNH_LIC Copyright 1994-2022 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
MODULE MODE_RMC01_3D
IMPLICIT NONE
CONTAINS
SUBROUTINE RMC01_3D(D,CST,PDXX,PDYY,PDZZ,PDIRCOSZW,PPHIM,PZC)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!     ##############################################################
!
!!****  *RMC01* -
!!
!!    PURPOSE
!!    -------
!!    This routine computes 3D parts of the rmc01.f90 routine 
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
!!    AUTHOR
!!    ------
!!
!!      Q. Rodier  - Meteo-France -
!!
!!    MODIFICATIONS
!!    -------------
!!     Original     18/08/022
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST, ONLY : CST_t
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
!
USE SHUMAN_PHY, ONLY: MYF_PHY, MXF_PHY
!
IMPLICIT NONE
!
!*       0.1   Declaration of arguments
TYPE(DIMPHYEX_t),         INTENT(IN)   :: D
TYPE(CST_t),              INTENT(IN)   :: CST
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    :: PDXX  ! width of grid mesh (X dir)
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    :: PDYY  ! width of grid mesh (Y dir)
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    :: PDZZ  ! width of vert. layers
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    :: PPHIM ! MO function
REAL, DIMENSION(D%NIT,D%NJT),         INTENT(IN)    :: PDIRCOSZW ! Director Cosinus
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(OUT)   :: PZC  ! alt. where turb. is isotr.
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: ZWORK1, ZWORK2
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: ZDH  ! hor. grid mesh
!
INTEGER :: IKB,IKE    ! first,last physical level
INTEGER :: IKTB,IKTE  ! start, end of k loops in physical domain
INTEGER :: JK,JI,JJ   ! loop counter
INTEGER :: IIE,IIB,IJE,IJB,IIU,IJU
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('RMC01_3D',0,ZHOOK_HANDLE)
IKTB=D%NKTB          
IKTE=D%NKTE
IKB=D%NKB
IKE=D%NKE
IIE=D%NIEC
IIB=D%NIBC
IJE=D%NJEC
IJB=D%NJBC
IIU=D%NIT
IJU=D%NJT
!
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
      PZC(JI,JJ,JK) = 2.*MIN(PPHIM(JI,JJ,JK),1.)/CST%XKARMAN    &
      * MAX( PDZZ(JI,JJ,JK)*PDIRCOSZW(JI,JJ) , & 
      ZDH(JI,JJ,JK)/PDIRCOSZW(JI,JJ)/3. )
    ENDDO
  ENDDO
END DO
!
IF (LHOOK) CALL DR_HOOK('RMC01_3D',1,ZHOOK_HANDLE)
END SUBROUTINE RMC01_3D
END MODULE MODE_RMC01_3D
