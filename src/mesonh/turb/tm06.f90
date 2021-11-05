!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 newsrc 2006/05/18 13:07:25
!-----------------------------------------------------------------
!    ################ 
     MODULE MODI_TM06  
!    ################ 
!
INTERFACE
!
      SUBROUTINE TM06(KKA,KKU,KKL,PTHVREF,PBL_DEPTH,PZZ,PSFTH,PMWTH,PMTH2)
!
INTEGER,                INTENT(IN)   :: KKA           !near ground array index  
INTEGER,                INTENT(IN)   :: KKU           !uppest atmosphere array index
INTEGER,                INTENT(IN)   :: KKL           !vert. levels type 1=MNH -1=ARO
REAL, DIMENSION(:,:,:), INTENT(IN) :: PTHVREF    ! reference potential temperature
REAL, DIMENSION(:,:),   INTENT(IN) :: PBL_DEPTH ! boundary layer height
REAL, DIMENSION(:,:,:), INTENT(IN) :: PZZ        ! altitude of flux levels
REAL, DIMENSION(:,:),   INTENT(IN) :: PSFTH      ! surface heat flux
REAL, DIMENSION(:,:,:), INTENT(OUT):: PMWTH      ! w'2th'
REAL, DIMENSION(:,:,:), INTENT(OUT):: PMTH2      ! w'th'2
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE TM06
!
END INTERFACE
!
END MODULE MODI_TM06
!
!     #################################################################
      SUBROUTINE TM06(KKA,KKU,KKL,PTHVREF,PBL_DEPTH,PZZ,PSFTH,PMWTH,PMTH2)
!     #################################################################
!
!
!!****  *TM06* - computes the Third Order Moments
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
!!
!!     TOMs are deduced from convective normalized TOMs according to Tomas and
!!     Masson 2006
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
!!      V. MAsson and S. Tomas  * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original         sept. 2005
!!
!! --------------------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
!
USE MODD_PARAMETERS, ONLY : XUNDEF
USE MODD_CST,        ONLY : XG
USE MODD_PARAMETERS, ONLY : JPVEXT_TURB

!
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
INTEGER,                INTENT(IN) :: KKA           !near ground array index  
INTEGER,                INTENT(IN) :: KKU           !uppest atmosphere array index
INTEGER,                INTENT(IN) :: KKL           !vert. levels type 1=MNH -1=ARO
REAL, DIMENSION(:,:,:), INTENT(IN) :: PTHVREF    ! reference potential temperature
REAL, DIMENSION(:,:),   INTENT(IN) :: PBL_DEPTH ! boundary layer height
REAL, DIMENSION(:,:,:), INTENT(IN) :: PZZ        ! altitude of flux levels
REAL, DIMENSION(:,:),   INTENT(IN) :: PSFTH      ! surface heat flux
REAL, DIMENSION(:,:,:), INTENT(OUT):: PMWTH      ! w'2th'
REAL, DIMENSION(:,:,:), INTENT(OUT):: PMTH2      ! w'th'2
!
!-------------------------------------------------------------------------------
!
!       0.2  declaration of local variables
!
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)):: ZZ_O_H ! normalized height z/h (where h=BL height)
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2))            :: ZWSTAR ! normalized convective velocity w*
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2))            :: ZTSTAR ! normalized temperature velocity w*
!
INTEGER                                             :: JK     ! loop counter
INTEGER                                             :: IKT    ! vertical size
INTEGER                                             :: IKTB,IKTE,IKB,IKE ! vertical levels
!----------------------------------------------------------------------------
!
IKT=SIZE(PZZ,3)          
IKTB=1+JPVEXT_TURB              
IKTE=IKT-JPVEXT_TURB
IKB=KKA+JPVEXT_TURB*KKL
IKE=KKU-JPVEXT_TURB*KKL

!
!
!* w* and T*
!
WHERE(PSFTH>0.)
  ZWSTAR = ((XG/PTHVREF(:,:,IKB))*PSFTH*PBL_DEPTH)**(1./3.)
  ZTSTAR = PSFTH / ZWSTAR
ELSEWHERE
  ZWSTAR = 0.
  ZTSTAR = 0.
END WHERE
!
!
!* normalized height
!
ZZ_O_H = XUNDEF
DO JK=1,IKT
  WHERE (PBL_DEPTH/=XUNDEF)
    ZZ_O_H(:,:,JK) = (PZZ(:,:,JK)-PZZ(:,:,IKB)) / PBL_DEPTH(:,:)
  END WHERE
END DO
!
!* w'th'2
!
PMTH2 = 0.
WHERE(ZZ_O_H < 0.95 .AND. ZZ_O_H/=XUNDEF)
  PMTH2(:,:,:) = 4.*(MAX(ZZ_O_H,0.))**0.4*(ZZ_O_H-0.95)**2
END WHERE
DO JK=IKTB+1,IKTE-1
  PMTH2(:,:,JK) = PMTH2(:,:,JK) * ZTSTAR(:,:)**2*ZWSTAR(:,:)
END DO
PMTH2(:,:,IKE)=PMTH2(:,:,IKE) * ZTSTAR(:,:)**2*ZWSTAR(:,:)
PMTH2(:,:,KKU)=PMTH2(:,:,KKU) * ZTSTAR(:,:)**2*ZWSTAR(:,:)

!
!
!* w'2th'
!
PMWTH = 0.
WHERE(ZZ_O_H <0.9 .AND. ZZ_O_H/=XUNDEF)
  PMWTH(:,:,:) = MAX(-7.9*(ABS(ZZ_O_H-0.35))**2.9 * (ABS(ZZ_O_H-1.))**0.58 + 0.37, 0.)
END WHERE
DO JK=IKTB+1,IKTE-1
  PMWTH(:,:,JK) = PMWTH(:,:,JK) * ZWSTAR(:,:)**2*ZTSTAR(:,:)
END DO
PMWTH(:,:,IKE) = PMWTH(:,:,IKE) * ZWSTAR(:,:)**2*ZTSTAR(:,:)
PMWTH(:,:,KKU) = PMWTH(:,:,KKU) * ZWSTAR(:,:)**2*ZTSTAR(:,:)
!
!----------------------------------------------------------------------------
END SUBROUTINE TM06
