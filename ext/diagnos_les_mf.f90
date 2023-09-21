!MNH_LIC Copyright 2009-2020 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!    ###########################
     MODULE MODI_DIAGNOS_LES_MF
!    ###########################
!
INTERFACE
!
!     #################################################################
      SUBROUTINE DIAGNOS_LES_MF(KIU,KJU,KKU,PTIME_LES,               &
                                PTHL_UP,PRT_UP,PRV_UP,PRC_UP,PRI_UP, &
                                PU_UP, PV_UP, PTHV_UP, PW_UP,        &
                                PFRAC_UP,PEMF,PDETR,PENTR,           &
                                PWTHMF,PWTHVMF,PWRTMF,               &
                                PWUMF,PWVMF,                         &
                                KKLCL,KKETL,KKCTL)
!     #################################################################
!
!*                    1.1  Declaration of Arguments
!
use modd_precision, only: MNHTIME
!
INTEGER,                INTENT(IN)  :: KIU, KJU, KKU ! 3D grid size
REAL(kind=MNHTIME), DIMENSION(2), INTENT(OUT) :: PTIME_LES
REAL, DIMENSION(:,:,:),   INTENT(IN)  :: PTHL_UP,PRT_UP,PRV_UP,&
                                       PRC_UP,PRI_UP   ! updraft properties
REAL, DIMENSION(:,:,:),   INTENT(IN)  :: PU_UP, PV_UP
REAL, DIMENSION(:,:,:),   INTENT(IN)  :: PTHV_UP,PW_UP,&
                                       PFRAC_UP,PEMF,PDETR,PENTR
REAL, DIMENSION(:,:,:),   INTENT(IN)  :: PWTHMF,PWTHVMF,PWRTMF, &
                                       PWUMF,PWVMF
INTEGER, DIMENSION(:),  INTENT(IN)  :: KKLCL,KKETL,KKCTL
                                           

END SUBROUTINE DIAGNOS_LES_MF

END INTERFACE
!
END MODULE MODI_DIAGNOS_LES_MF
!
!     #################################################################
      SUBROUTINE DIAGNOS_LES_MF(KIU,KJU,KKU,PTIME_LES,               &
                                PTHL_UP,PRT_UP,PRV_UP,PRC_UP,PRI_UP, &
                                PU_UP, PV_UP, PTHV_UP, PW_UP,        &
                                PFRAC_UP,PEMF,PDETR,PENTR,           &
                                PWTHMF,PWTHVMF,PWRTMF,               &
                                PWUMF,PWVMF,                         &
                                KKLCL,KKETL,KKCTL)
!     #################################################################
!!
!!****  *DIAGNOS_LES_MF* - Edit in File the updraft properties as
!!                         LES diagnostics 
!!
!!    PURPOSE
!!    -------
!!****  The purpose of this routine is to write updraft variable as
!!       LES diagnostics         
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!      
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!     REFERENCE
!!     ---------
!!      
!!
!!     AUTHOR
!!     ------
!!     J.pergaud
!
! Modifications:
!  V. Masson      09/2010: Optimization
!  P. Wautelet 28/03/2019: use MNHTIME for time measurement variables
!! --------------------------------------------------------------------------
!
!*      0. DECLARATIONS
!          ------------
!
USE MODD_LES
use modd_precision, only: MNHTIME
!
USE MODE_MNH_TIMING
!
USE MODI_LES_VER_INT
USE MODI_LES_MEAN_ll
USE MODI_SHUMAN
!
IMPLICIT NONE

!*                    0.1  Declaration of Arguments
!
!
INTEGER,                INTENT(IN)  :: KIU, KJU, KKU ! 3D grid size
REAL(kind=MNHTIME), DIMENSION(2), INTENT(OUT) :: PTIME_LES
REAL, DIMENSION(:,:,:),   INTENT(IN)  :: PTHL_UP,PRT_UP,PRV_UP,&
                                       PRC_UP,PRI_UP   ! updraft properties
REAL, DIMENSION(:,:,:),   INTENT(IN)  :: PU_UP, PV_UP
REAL, DIMENSION(:,:,:),   INTENT(IN)  :: PTHV_UP,PW_UP,&
                                       PFRAC_UP,PEMF,PDETR,PENTR
REAL, DIMENSION(:,:,:),   INTENT(IN)  :: PWTHMF,PWTHVMF,PWRTMF, &
                                       PWUMF,PWVMF
INTEGER, DIMENSION(:),  INTENT(IN)  :: KKLCL,KKETL,KKCTL
                                           
!
!
!                     0.2  Declaration of local variables
!
REAL, DIMENSION(:,:,:), ALLOCATABLE   :: ZTHLMFFLX_LES,ZRTMFFLX_LES, &
                                         ZTHVMFFLX_LES,ZUMFFLX_LES, &
                                         ZVMFFLX_LES
REAL, DIMENSION(:,:,:), ALLOCATABLE   :: ZTHLUP_MF_LES,ZRTUP_MF_LES, &
                                         ZRCUP_MF_LES,ZEMF_MF_LES,   &
                                         ZDETR_MF_LES, ZENTR_MF_LES, & 
                                         ZWUP_MF_LES,ZFRACUP_MF_LES, &
                                         ZTHVUP_MF_LES,ZRVUP_MF_LES, &
                                         ZRIUP_MF_LES
REAL(kind=MNHTIME), DIMENSION(2) :: ZTIME1, ZTIME2
!------------------------------------------------------------------------
!

CALL SECOND_MNH2(ZTIME1)

 IF (LLES_CALL) THEN

    ALLOCATE( ZTHLUP_MF_LES(KIU,KJU,NLES_K) )
    ALLOCATE( ZRTUP_MF_LES(KIU,KJU,NLES_K) )
    ALLOCATE( ZRVUP_MF_LES(KIU,KJU,NLES_K) )
    ALLOCATE( ZRCUP_MF_LES(KIU,KJU,NLES_K) )
    ALLOCATE( ZRIUP_MF_LES(KIU,KJU,NLES_K) )
    ALLOCATE( ZEMF_MF_LES(KIU,KJU,NLES_K) )
    ALLOCATE( ZDETR_MF_LES(KIU,KJU,NLES_K) )
    ALLOCATE( ZENTR_MF_LES(KIU,KJU,NLES_K) )
    ALLOCATE( ZWUP_MF_LES(KIU,KJU,NLES_K) )
    ALLOCATE( ZFRACUP_MF_LES(KIU,KJU,NLES_K) )
    ALLOCATE( ZTHVUP_MF_LES(KIU,KJU,NLES_K) )

    ALLOCATE( ZTHLMFFLX_LES(KIU,KJU,NLES_K) )
    ALLOCATE( ZRTMFFLX_LES (KIU,KJU,NLES_K) )
    ALLOCATE( ZTHVMFFLX_LES(KIU,KJU,NLES_K) )
    ALLOCATE( ZUMFFLX_LES  (KIU,KJU,NLES_K) )
    ALLOCATE( ZVMFFLX_LES  (KIU,KJU,NLES_K) )

    
    CALL LES_VER_INT(MZF(PWTHMF)  ,ZTHLMFFLX_LES  )
    CALL LES_MEAN_ll(ZTHLMFFLX_LES,LLES_CURRENT_CART_MASK, &
                    X_LES_SUBGRID_WTHLMF(:,NLES_CURRENT_TCOUNT,1))
    
    CALL LES_VER_INT( MZF(PWRTMF)  ,ZRTMFFLX_LES  )
    CALL LES_MEAN_ll (ZRTMFFLX_LES , LLES_CURRENT_CART_MASK,          &
                    X_LES_SUBGRID_WRTMF(:,NLES_CURRENT_TCOUNT,1)     )
                    
    CALL LES_VER_INT( MZF(PWUMF)  ,ZUMFFLX_LES  )
    CALL LES_MEAN_ll (ZUMFFLX_LES , LLES_CURRENT_CART_MASK,              &
                    X_LES_SUBGRID_WUMF(:,NLES_CURRENT_TCOUNT,1)     )
                    
    CALL LES_VER_INT( MZF(PWVMF)  ,ZVMFFLX_LES  )
    CALL LES_MEAN_ll (ZVMFFLX_LES , LLES_CURRENT_CART_MASK,                   &
                    X_LES_SUBGRID_WVMF(:,NLES_CURRENT_TCOUNT,1)     )
                    
    CALL LES_VER_INT( MZF(PWTHVMF)  ,ZTHVMFFLX_LES  )
    CALL LES_MEAN_ll (ZTHVMFFLX_LES , LLES_CURRENT_CART_MASK,    &
                    X_LES_SUBGRID_WTHVMF(:,NLES_CURRENT_TCOUNT,1)     )
                    

    CALL LES_VER_INT( MZF(PTHL_UP)  ,ZTHLUP_MF_LES  )
    CALL LES_MEAN_ll (ZTHLUP_MF_LES , LLES_CURRENT_CART_MASK,    &
                   X_LES_SUBGRID_THLUP_MF(:,NLES_CURRENT_TCOUNT,1)     )
                   
    CALL LES_VER_INT( MZF(PRT_UP)  ,ZRTUP_MF_LES  )
    CALL LES_MEAN_ll (ZRTUP_MF_LES , LLES_CURRENT_CART_MASK,       &
                    X_LES_SUBGRID_RTUP_MF(:,NLES_CURRENT_TCOUNT,1)     )
                    
    CALL LES_VER_INT( MZF(PRV_UP)  ,ZRVUP_MF_LES  )
    CALL LES_MEAN_ll (ZRVUP_MF_LES , LLES_CURRENT_CART_MASK,       &
                    X_LES_SUBGRID_RVUP_MF(:,NLES_CURRENT_TCOUNT,1)     )
                    
    CALL LES_VER_INT( MZF(PRC_UP)  ,ZRCUP_MF_LES  )
    CALL LES_MEAN_ll (ZRCUP_MF_LES , LLES_CURRENT_CART_MASK,        &
                    X_LES_SUBGRID_RCUP_MF(:,NLES_CURRENT_TCOUNT,1)     )
                    
    CALL LES_VER_INT( MZF(PRI_UP)  ,ZRIUP_MF_LES  )
    CALL LES_MEAN_ll (ZRIUP_MF_LES , LLES_CURRENT_CART_MASK,        &
                    X_LES_SUBGRID_RIUP_MF(:,NLES_CURRENT_TCOUNT,1)     )            
                    
    CALL LES_VER_INT( MZF(PEMF)  ,ZEMF_MF_LES  )
    CALL LES_MEAN_ll (ZEMF_MF_LES , LLES_CURRENT_CART_MASK,       &
                   X_LES_SUBGRID_MASSFLUX(:,NLES_CURRENT_TCOUNT,1)     )
                   
    CALL LES_VER_INT( MZF(PDETR)  ,ZDETR_MF_LES  )
    CALL LES_MEAN_ll (ZDETR_MF_LES , LLES_CURRENT_CART_MASK,       &
                   X_LES_SUBGRID_DETR(:,NLES_CURRENT_TCOUNT,1)         )
                   
    CALL LES_VER_INT( MZF(PENTR)  ,ZENTR_MF_LES  )
    CALL LES_MEAN_ll (ZENTR_MF_LES , LLES_CURRENT_CART_MASK,       &
                   X_LES_SUBGRID_ENTR(:,NLES_CURRENT_TCOUNT,1)     )
                   
    CALL LES_VER_INT( MZF(PW_UP)  ,ZWUP_MF_LES  )
    CALL LES_MEAN_ll (ZWUP_MF_LES , LLES_CURRENT_CART_MASK,       &
                   X_LES_SUBGRID_WUP_MF(:,NLES_CURRENT_TCOUNT,1)     )
                   
    CALL LES_VER_INT( MZF(PFRAC_UP)  ,ZFRACUP_MF_LES  )
    CALL LES_MEAN_ll (ZFRACUP_MF_LES , LLES_CURRENT_CART_MASK,       &
                   X_LES_SUBGRID_FRACUP(:,NLES_CURRENT_TCOUNT,1)     )
                   
    CALL LES_VER_INT( MZF(PTHV_UP)  ,ZTHVUP_MF_LES  )
    CALL LES_MEAN_ll (ZTHVUP_MF_LES , LLES_CURRENT_CART_MASK,       &
                   X_LES_SUBGRID_THVUP_MF(:,NLES_CURRENT_TCOUNT,1)     )
                   
     
        
    DEALLOCATE( ZTHLMFFLX_LES )
    DEALLOCATE( ZRTMFFLX_LES )
    DEALLOCATE( ZTHVMFFLX_LES )
    DEALLOCATE( ZUMFFLX_LES  )
    DEALLOCATE( ZVMFFLX_LES  )

    
    DEALLOCATE( ZTHLUP_MF_LES )
    DEALLOCATE( ZRTUP_MF_LES )
    DEALLOCATE( ZRVUP_MF_LES )    
    DEALLOCATE( ZRCUP_MF_LES )
    DEALLOCATE( ZRIUP_MF_LES )
    DEALLOCATE( ZENTR_MF_LES )
    DEALLOCATE( ZDETR_MF_LES )
    DEALLOCATE( ZEMF_MF_LES )
    DEALLOCATE( ZWUP_MF_LES )
    DEALLOCATE( ZFRACUP_MF_LES )
    DEALLOCATE( ZTHVUP_MF_LES )
    
ENDIF

CALL SECOND_MNH2(ZTIME2)
PTIME_LES =             ZTIME2 - ZTIME1
XTIME_LES = XTIME_LES + ZTIME2 - ZTIME1

END SUBROUTINE DIAGNOS_LES_MF           
