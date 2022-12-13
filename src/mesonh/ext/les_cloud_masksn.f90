!MNH_LIC Copyright 2006-2020 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     #######################
      SUBROUTINE  LES_CLOUD_MASKS_n
!     #######################
!
!
!!****  *LES_MASKS_n* initializes the masks for clouds
!!                         
!!
!!    PURPOSE
!!    -------
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      V. Masson
!!
!!    MODIFICATIONS
!!    -------------
!!      Original         07/2006
!!      P. Aumond        10/2009 Add possibility of user maskS
!!       F.Couvreux      06/2011 : Conditional sampling
!!       C.Lac           10/2014 : Correction on user masks   
!!       Q.Rodier        05/2019 : Missing parallelization 
!!
!! --------------------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
!
USE MODD_LES
USE MODD_LES_n
USE MODD_FIELD_n
USE MODD_CONF_n
USE MODD_CST    , ONLY : XRD, XRV
USE MODD_NSV    , ONLY : NSV_CSBEG, NSV_CSEND, NSV_CS
USE MODD_GRID_n , ONLY : XZHAT
USE MODD_CONDSAMP
!
USE MODE_ll
!
USE MODI_LES_VER_INT
USE MODI_LES_MEAN_ll
USE MODI_SHUMAN
!
IMPLICIT NONE
!
!
!       0.2  declaration of local variables
!
INTEGER :: JK                                  ! vertical loop counter
INTEGER :: JI                                  ! loop index on masks  
INTEGER :: IIU, IJU,IIB,IJB,IIE,IJE            ! hor. indices
INTEGER :: IKU, KBASE, KTOP                    ! ver. index
INTEGER :: IRR, IRRC, IRRR, IRRI, IRRS, IRRG   ! moist variables indices
INTEGER :: JSV                                  ! ind of scalars
!
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZRT  ! total water
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZTHV ! Virtual potential temperature
!
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZW_LES    ! W     on LES vertical grid
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZRC_LES   ! Rc    on LES vertical grid
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZRI_LES   ! Ri    on LES vertical grid
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZRT_LES   ! Rt    on LES vertical grid
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZTHV_LES  ! thv   on LES vertical grid
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZSV_LES  ! thv   on LES vertical grid
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZTHV_ANOM ! thv-thv_mean   on LES vertical grid
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZSV_ANOM ! sv-sv_mean
REAL, DIMENSION(:,:),     ALLOCATABLE :: ZSTD_SV
REAL, DIMENSION(:,:),     ALLOCATABLE :: ZSTD_SVTRES ! threshold of sv
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZWORK3D,ZWORK3DB
REAL, DIMENSION(:),       ALLOCATABLE :: ZWORK1D
REAL, DIMENSION(:),       ALLOCATABLE :: ZMEANRC
!
!
!-------------------------------------------------------------------------------
!
CALL GET_DIM_EXT_ll('B',IIU,IJU)
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
!
IKU = SIZE(XTHT,3)
!
!-------------------------------------------------------------------------------
!
!*      1.0  Thermodynamical computations
!            ----------------------------
!
ALLOCATE(ZRT      (IIU,IJU,IKU))
ALLOCATE(ZMEANRC (IKU))
ZRT = 0.
!
IRR=0
IF (LUSERV) THEN
  IRR=IRR+1
  ZRT = ZRT + XRT(:,:,:,1)
END IF
IF (LUSERC) THEN
  IRR=IRR+1
  IRRC=IRR
  ZRT = ZRT + XRT(:,:,:,IRRC)
END IF
IF (LUSERR) THEN
  IRR=IRR+1
  IRRR=IRR
  ZRT = ZRT + XRT(:,:,:,IRRR)
END IF
IF (LUSERI) THEN
  IRR=IRR+1
  IRRI=IRR
  ZRT = ZRT + XRT(:,:,:,IRRI)
END IF
IF (LUSERS) THEN
  IRR=IRR+1
  IRRS=IRR
  ZRT = ZRT + XRT(:,:,:,IRRS)
END IF
IF (LUSERG) THEN
  IRR=IRR+1
  IRRG=IRR
  ZRT = ZRT + XRT(:,:,:,IRRG)
END IF
!
!
!* computes fields on the LES grid in order to compute masks
!
ALLOCATE(ZTHV     (IIU,IJU,IKU))
ZTHV = XTHT
IF (LUSERV) ZTHV=ZTHV*(1.+XRV/XRD*XRT(:,:,:,1))/(1.+ZRT(:,:,:))
!
!-------------------------------------------------------------------------------
!
!*      2.0  Fields on LES grid
!            ------------------
!
!* allocates fields on the LES grid
!
!
ALLOCATE(ZW_LES   (IIU,IJU,NLES_K))
ALLOCATE(ZRC_LES  (IIU,IJU,NLES_K))
ALLOCATE(ZRI_LES  (IIU,IJU,NLES_K))
ALLOCATE(ZRT_LES  (IIU,IJU,NLES_K))
ALLOCATE(ZTHV_LES (IIU,IJU,NLES_K))
ALLOCATE(ZTHV_ANOM(IIU,IJU,NLES_K))
ALLOCATE(ZSV_LES (IIU,IJU,NLES_K,NSV_CS))
ALLOCATE(ZSV_ANOM(IIU,IJU,NLES_K,NSV_CS))
ALLOCATE(ZSTD_SV(NLES_K,NSV_CS))
ALLOCATE(ZSTD_SVTRES(NLES_K,NSV_CS))
ALLOCATE(ZWORK1D(NLES_K))
ALLOCATE(ZWORK3D(IIU,IJU,IKU))
ALLOCATE(ZWORK3DB(IIU,IJU,NLES_K))
!
ZWORK1D=0.
ZWORK3D=0.  
ZWORK3DB=0.
!
CALL LES_VER_INT(MZF(XWT), ZW_LES)
IF (NSV_CS>0) THEN
  DO JSV=NSV_CSBEG, NSV_CSEND
    CALL LES_VER_INT(  XSVT(:,:,:,JSV),  &
                       ZSV_LES(:,:,:,JSV-NSV_CSBEG+1) )
  END DO
END IF
IF (LUSERC) THEN
  CALL LES_VER_INT(XRT(:,:,:,IRRC), ZRC_LES)
ELSE
  ZRC_LES = 0.
END IF
IF (LUSERI) THEN
  CALL LES_VER_INT(XRT(:,:,:,IRRI), ZRI_LES)
ELSE
  ZRI_LES = 0.
END IF
CALL LES_VER_INT(ZRT, ZRT_LES)
CALL LES_VER_INT(ZTHV, ZTHV_LES)
CALL LES_ANOMALY_FIELD(ZTHV,ZTHV_ANOM)
!
IF (NSV_CS>0) THEN
 DO JSV=NSV_CSBEG, NSV_CSEND
  ZWORK3D(:,:,:)=XSVT(:,:,:,JSV)
  CALL LES_ANOMALY_FIELD(ZWORK3D,ZWORK3DB)
  ZSV_ANOM(:,:,:,JSV-NSV_CSBEG+1)=ZWORK3DB(:,:,:)
  CALL LES_STDEV(ZWORK3DB,ZWORK1D)
  ZSTD_SV(:,JSV-NSV_CSBEG+1)=ZWORK1D(:)
  DO JK=1,NLES_K
       ZSTD_SVTRES(JK,JSV-NSV_CSBEG+1)=SUM(ZSTD_SV(1:JK,JSV-NSV_CSBEG+1))/(1.*JK)
  END DO
 END DO
END IF
!
DEALLOCATE(ZTHV     )
DEALLOCATE(ZWORK3D)
DEALLOCATE(ZWORK3DB)
DEALLOCATE(ZWORK1D)
!
!-------------------------------------------------------------------------------
!
!*      3.0  Cloud mask
!            ----------
!
IF (LLES_NEB_MASK) THEN
  CALL LES_ALLOCATE('LLES_CURRENT_NEB_MASK',(/IIU,IJU,NLES_K/))
  LLES_CURRENT_NEB_MASK (:,:,:) = .FALSE.
  WHERE ((ZRC_LES(IIB:IIE,IJB:IJE,:)>1.E-6 .OR. ZRI_LES(IIB:IIE,IJB:IJE,:)>1.E-6) .AND. ZW_LES(IIB:IIE,IJB:IJE,:)>0.)
    LLES_CURRENT_NEB_MASK (IIB:IIE,IJB:IJE,:) = .TRUE.
  END WHERE
END IF
!
!-------------------------------------------------------------------------------
!
!*      4.0  Cloud core mask
!            ---------------
!
IF (LLES_CORE_MASK) THEN
  CALL LES_ALLOCATE('LLES_CURRENT_CORE_MASK',(/IIU,IJU,NLES_K/))
  LLES_CURRENT_CORE_MASK (:,:,:) = .FALSE.
  WHERE ((ZRC_LES(IIB:IIE,IJB:IJE,:)>1.E-6 .OR. ZRI_LES(IIB:IIE,IJB:IJE,:)>1.E-6) &
         .AND. ZW_LES(IIB:IIE,IJB:IJE,:)>0. .AND. ZTHV_ANOM(IIB:IIE,IJB:IJE,:)>0.)  
    LLES_CURRENT_CORE_MASK (IIB:IIE,IJB:IJE,:) = .TRUE.
  END WHERE
END IF
!
!-------------------------------------------------------------------------------
!
!*      4.0  Conditional sampling mask
!            -------------------------
!
IF (LLES_CS_MASK) THEN
!
  CALL LES_MEAN_ll(ZRC_LES, LLES_CURRENT_CART_MASK, ZMEANRC  )
  CALL LES_ALLOCATE('LLES_CURRENT_CS1_MASK',(/IIU,IJU,NLES_K/))
  LLES_CURRENT_CS1_MASK(:,:,:) = .FALSE.
  IF (NSV_CS >= 2) THEN
    CALL LES_ALLOCATE('LLES_CURRENT_CS2_MASK',(/IIU,IJU,NLES_K/))
    LLES_CURRENT_CS2_MASK(:,:,:) = .FALSE.
    IF (NSV_CS == 3) THEN
      CALL LES_ALLOCATE('LLES_CURRENT_CS3_MASK',(/IIU,IJU,NLES_K/))
      LLES_CURRENT_CS3_MASK (:,:,:) = .FALSE.
    END IF
  END IF

!
! Cloud top and base computation                   
!
  KBASE=2
  KTOP=NLES_K
  DO JK=2,NLES_K
   IF ((ZMEANRC(JK) > 1.E-7) .AND. (KBASE == 2)) KBASE=JK
   IF ((ZMEANRC(JK) < 1.E-7) .AND. (KBASE > 2) .AND. (KTOP == NLES_K)) &
              KTOP=JK-1
  END DO
!
 DO JSV=NSV_CSBEG, NSV_CSEND
  DO JK=2,NLES_K
    IF (ZSTD_SV(JK,JSV-NSV_CSBEG+1) < 0.05*ZSTD_SVTRES(JK,JSV-NSV_CSBEG+1)) &
       ZSTD_SV(JK,JSV-NSV_CSBEG+1)=0.05*ZSTD_SVTRES(JK,JSV-NSV_CSBEG+1)
! case no cloud top and base                    
    IF (JSV == NSV_CSBEG) THEN
     IF ((KBASE ==2) .AND. (KTOP == NLES_K)) THEN
      WHERE (ZW_LES(IIB:IIE,IJB:IJE,JK)>0. .AND. ZSV_ANOM(IIB:IIE,IJB:IJE,JK,JSV-NSV_CSBEG+1) >  &
             XSCAL(JSV-NSV_CSBEG+1) * ZSTD_SV(JK,JSV-NSV_CSBEG+1))
          LLES_CURRENT_CS1_MASK (IIB:IIE,IJB:IJE,JK) = .TRUE.
      END WHERE       
     END IF
!
! case cloud top and base defined                    
!
     IF (XZHAT(JK) < XZHAT(KBASE)+(XZHAT(KTOP)-XZHAT(KBASE))/4.) THEN      
      WHERE (ZW_LES(IIB:IIE,IJB:IJE,JK)>0. .AND. ZSV_ANOM(IIB:IIE,IJB:IJE,JK,JSV-NSV_CSBEG+1) > &
             XSCAL(JSV-NSV_CSBEG+1) *ZSTD_SV(JK,JSV-NSV_CSBEG+1))
          LLES_CURRENT_CS1_MASK (IIB:IIE,IJB:IJE,JK) = .TRUE.
      END WHERE
     END IF     
!
     IF (XZHAT(JK) >= XZHAT(KBASE)+(XZHAT(KTOP)-XZHAT(KBASE))/4.) THEN
      WHERE (ZW_LES(IIB:IIE,IJB:IJE,JK)>0. .AND. ZSV_ANOM(IIB:IIE,IJB:IJE,JK,JSV-NSV_CSBEG+1) > &
             XSCAL(JSV-NSV_CSBEG+1) * ZSTD_SV(JK,JSV-NSV_CSBEG+1) .AND. &
          ZRC_LES(IIB:IIE,IJB:IJE,JK)>1.E-6)
          LLES_CURRENT_CS1_MASK (IIB:IIE,IJB:IJE,JK) = .TRUE.
      END WHERE  
     END IF
    ELSE IF ( JSV == NSV_CSBEG + 1 ) THEN
     IF ((KBASE ==2) .AND. (KTOP == NLES_K)) THEN
      WHERE ( ZSV_ANOM(IIB:IIE,IJB:IJE,JK,JSV-NSV_CSBEG+1) > &
              XSCAL(JSV-NSV_CSBEG+1) * ZSTD_SV(JK,JSV-NSV_CSBEG+1))
          LLES_CURRENT_CS2_MASK (IIB:IIE,IJB:IJE,JK) = .TRUE.
      END WHERE       
     END IF
!
! case cloud top and base defined                    
!
     IF (XZHAT(JK) < XZHAT(KBASE)+(XZHAT(KTOP)-XZHAT(KBASE))/4.) THEN      
      WHERE (ZSV_ANOM(IIB:IIE,IJB:IJE,JK,JSV-NSV_CSBEG+1) >  &
             XSCAL(JSV-NSV_CSBEG+1) *ZSTD_SV(JK,JSV-NSV_CSBEG+1))
          LLES_CURRENT_CS2_MASK (IIB:IIE,IJB:IJE,JK) = .TRUE.
      END WHERE
     END IF     
!
     IF (XZHAT(JK) >= XZHAT(KBASE)+(XZHAT(KTOP)-XZHAT(KBASE))/4.) THEN
      WHERE (ZSV_ANOM(IIB:IIE,IJB:IJE,JK,JSV-NSV_CSBEG+1) > &
             XSCAL(JSV-NSV_CSBEG+1) * ZSTD_SV(JK,JSV-NSV_CSBEG+1)) 
          LLES_CURRENT_CS2_MASK (IIB:IIE,IJB:IJE,JK) = .TRUE.
      END WHERE  
     END IF
!
    ELSE 
     IF ((KBASE ==2) .AND. (KTOP == NLES_K)) THEN
      WHERE ( ZSV_ANOM(IIB:IIE,IJB:IJE,JK,JSV-NSV_CSBEG+1) > &
              XSCAL(JSV-NSV_CSBEG+1) * ZSTD_SV(JK,JSV-NSV_CSBEG+1))
          LLES_CURRENT_CS3_MASK (IIB:IIE,IJB:IJE,JK) = .TRUE.
      END WHERE       
     END IF
!
! case cloud top and base defined                    
!
     IF (XZHAT(JK) < XZHAT(KBASE)+(XZHAT(KTOP)-XZHAT(KBASE))/4.) THEN      
      WHERE (ZSV_ANOM(IIB:IIE,IJB:IJE,JK,JSV-NSV_CSBEG+1) > &
             XSCAL(JSV-NSV_CSBEG+1) *ZSTD_SV(JK,JSV-NSV_CSBEG+1))
          LLES_CURRENT_CS3_MASK (IIB:IIE,IJB:IJE,JK) = .TRUE.
      END WHERE
     END IF     
!
     IF (XZHAT(JK) >= XZHAT(KBASE)+(XZHAT(KTOP)-XZHAT(KBASE))/4.) THEN
      WHERE (ZSV_ANOM(IIB:IIE,IJB:IJE,JK,JSV-NSV_CSBEG+1) > &
             XSCAL(JSV-NSV_CSBEG+1) * ZSTD_SV(JK,JSV-NSV_CSBEG+1)) 
          LLES_CURRENT_CS3_MASK (IIB:IIE,IJB:IJE,JK) = .TRUE.
      END WHERE  
     END IF
    END IF
  END DO
 END DO
END IF
!
!-------------------------------------------------------------------------------
!
!*      5.0  User mask
!            ---------
!
IF (LLES_MY_MASK) THEN
  CALL LES_ALLOCATE('LLES_CURRENT_MY_MASKS',(/IIU,IJU,NLES_K,NLES_MASKS_USER/))
  DO JI=1,NLES_MASKS_USER
    LLES_CURRENT_MY_MASKS (IIB:IIE,IJB:IJE,:,JI) = .FALSE.
  END DO
! WHERE ((ZRC_LES + ZRI_LES) > 1.E-06) 
!    LLES_CURRENT_MY_MASKS (:,:,:,1) = .TRUE.
! END WHERE
!
END IF
!-------------------------------------------------------------------------------
!
DEALLOCATE(ZW_LES   )
DEALLOCATE(ZRC_LES  )
DEALLOCATE(ZRI_LES  )
DEALLOCATE(ZRT_LES  )
DEALLOCATE(ZTHV_LES )
DEALLOCATE(ZSV_LES  )
DEALLOCATE(ZTHV_ANOM)
DEALLOCATE(ZSV_ANOM)
DEALLOCATE(ZSTD_SV)
DEALLOCATE(ZSTD_SVTRES)
!-------------------------------------------------------------------------------
DEALLOCATE(ZRT )
DEALLOCATE(ZMEANRC)
!--------------------------------------------------------------------------------
!
CONTAINS
!
!--------------------------------------------------------------------------------
!
SUBROUTINE LES_ANOMALY_FIELD(PF,PF_ANOM)
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PF
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PF_ANOM

REAL, DIMENSION(SIZE(PF_ANOM,3)) :: ZMEAN
INTEGER :: JI, JJ

CALL LES_VER_INT(PF, PF_ANOM)
CALL LES_MEAN_ll(PF_ANOM, LLES_CURRENT_CART_MASK, ZMEAN  )
DO JJ=1,SIZE(PF_ANOM,2)
  DO JI=1,SIZE(PF_ANOM,1)
    PF_ANOM(JI,JJ,:) = PF_ANOM(JI,JJ,:) - ZMEAN(:)
  END DO
END DO

END SUBROUTINE LES_ANOMALY_FIELD
!--------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------
!
SUBROUTINE LES_STDEV(PF_ANOM,PF_STD)
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PF_ANOM
REAL, DIMENSION(:), INTENT(OUT) :: PF_STD

REAL, DIMENSION(SIZE(PF_ANOM,1),SIZE(PF_ANOM,2),SIZE(PF_ANOM,3)) :: Z2
INTEGER :: JK

Z2(:,:,:)=PF_ANOM(:,:,:)*PF_ANOM(:,:,:)
CALL LES_MEAN_ll(Z2, LLES_CURRENT_CART_MASK, PF_STD  )
DO JK=1,SIZE(PF_ANOM,3)
     PF_STD(JK)=SQRT(PF_STD(JK))
END DO

END SUBROUTINE LES_STDEV
!-------------------------------------------------------------------------------  
!
END SUBROUTINE LES_CLOUD_MASKS_n   
