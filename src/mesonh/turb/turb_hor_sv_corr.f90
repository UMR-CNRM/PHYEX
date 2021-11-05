!MNH_LIC Copyright 2002-2020 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!    ############################
     MODULE MODI_TURB_HOR_SV_CORR
!    ############################
!
INTERFACE  
!
      SUBROUTINE TURB_HOR_SV_CORR(KRR,KRRL,KRRI,                     &
                      PDXX,PDYY,PDZZ,PDZX,PDZY,                      &
                      PLM,PLEPS,PTKEM,PTHVREF,                       &
                      PTHLM,PRM,                                     &
                      PLOCPEXNM,PATHETA,PAMOIST,PSRCM,               &
                      PWM,PSVM                                       )
!
INTEGER,                  INTENT(IN)    ::  KRR          ! number of moist var.
INTEGER,                  INTENT(IN)    ::  KRRL         ! number of liquid var.
INTEGER,                  INTENT(IN)    ::  KRRI         ! number of ice var.
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PDXX, PDYY, PDZZ, PDZX, PDZY 
                                                         ! Metric coefficients
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PLM          ! mixing length
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PLEPS        ! dissipative length
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PTKEM        ! tke
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PTHVREF      ! reference Thv
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PTHLM        ! potential temperature at t-Delta t
REAL, DIMENSION(:,:,:,:), INTENT(IN)    ::  PRM          ! Mixing ratios at t-Delta t
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PLOCPEXNM    ! Lv(T)/Cp/Exnref at time t-1
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PATHETA      ! coefficients between 
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PAMOIST      ! s and Thetal and Rnp
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PSRCM        ! normalized 
                  ! 2nd-order flux s'r'c/2Sigma_s2 at t-1 multiplied by Lambda_3
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PWM          ! w at t-1
REAL, DIMENSION(:,:,:,:), INTENT(IN)    ::  PSVM         ! scalar var. at t-1
!
!
END SUBROUTINE TURB_HOR_SV_CORR
!
END INTERFACE
!
END MODULE MODI_TURB_HOR_SV_CORR
!     ################################################################
      SUBROUTINE TURB_HOR_SV_CORR(KRR,KRRL,KRRI,                     &
                      PDXX,PDYY,PDZZ,PDZX,PDZY,                      &
                      PLM,PLEPS,PTKEM,PTHVREF,                       &
                      PTHLM,PRM,                                     &
                      PLOCPEXNM,PATHETA,PAMOIST,PSRCM,               &
                      PWM,PSVM                                       )
!     ################################################################
!
!
!!****  *TURB_HOT_SV_CORR*  computes subgrid Sv2 and SvThv terms
!!
!!    PURPOSE
!!    -------
!!
!!     see TURB_HOR
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
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!
!!      V. Masson               * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!   Original  06/11/02
!!      JP Pinty       Feb 20, 2003 Add PFRAC_ICE
!! --------------------------------------------------------------------------
!
!*      0. DECLARATIONS
!          ------------
!
USE MODD_CST
USE MODD_CONF
USE MODD_CTURB
USE MODD_PARAMETERS
USE MODD_NSV, ONLY : NSV,NSV_LGBEG,NSV_LGEND
USE MODD_LES
USE MODD_BLOWSNOW
!
USE MODI_GRADIENT_M
USE MODI_GRADIENT_U
USE MODI_GRADIENT_V
USE MODI_GRADIENT_W
USE MODI_SHUMAN 
USE MODI_LES_MEAN_SUBGRID
USE MODI_EMOIST
USE MODI_ETHETA
!
USE MODI_SECOND_MNH
!
IMPLICIT NONE
!
!
!*       0.1  declaration of arguments
!
!
!
INTEGER,                  INTENT(IN)    ::  KRR          ! number of moist var.
INTEGER,                  INTENT(IN)    ::  KRRL         ! number of liquid var.
INTEGER,                  INTENT(IN)    ::  KRRI         ! number of ice var.
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PDXX, PDYY, PDZZ, PDZX, PDZY 
                                                         ! Metric coefficients
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PLM          ! mixing length
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PLEPS        ! dissipative length
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PTKEM        ! tke
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PTHVREF      ! reference Thv
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PTHLM        ! potential temperature at t-Delta t
REAL, DIMENSION(:,:,:,:), INTENT(IN)    ::  PRM          ! Mixing ratios at t-Delta t
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PLOCPEXNM    ! Lv(T)/Cp/Exnref at time t-1
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PATHETA      ! coefficients between 
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PAMOIST      ! s and Thetal and Rnp
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PSRCM        ! normalized 
                  ! 2nd-order flux s'r'c/2Sigma_s2 at t-1 multiplied by Lambda_3
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PWM          ! w at t-1
REAL, DIMENSION(:,:,:,:), INTENT(IN)    ::  PSVM         ! scalar var. at t-1
!
!
!
!*       0.2  declaration of local variables
!
REAL, DIMENSION(SIZE(PSVM,1),SIZE(PSVM,2),SIZE(PSVM,3))       &
                                     :: ZFLX, ZA
!
INTEGER             :: JSV          ! loop counter
!
REAL :: ZTIME1, ZTIME2
!
REAL :: ZCSVD  = 1.2  ! constant for scalar variance dissipation
REAL :: ZCTSVD = 2.4  ! constant for temperature - scalar covariance dissipation
REAL :: ZCQSVD = 2.4  ! constant for humidity - scalar covariance dissipation
!
REAL :: ZCSV          !constant for the scalar flux 
! ---------------------------------------------------------------------------
!
CALL SECOND_MNH(ZTIME1)
!
IF(LBLOWSNOW) THEN
! See Vionnet (PhD, 2012) for a complete discussion around the value of the Schmidt number for blowing snow variables        
   ZCSV= XCHF/XRSNOW 
ELSE
   ZCSV= XCHF
ENDIF
!
DO JSV=1,NSV
!
  IF (LNOMIXLG .AND. JSV >= NSV_LGBEG .AND. JSV<= NSV_LGEND) CYCLE
  !
  ! variance Sv2
  !
  IF (LLES_CALL) THEN
    IF (.NOT. L2D) THEN
      ZFLX(:,:,:) =  ZCSV / ZCSVD * PLM(:,:,:) * PLEPS(:,:,:) *   &
         (  GX_M_M(PSVM(:,:,:,JSV),PDXX,PDZZ,PDZX)**2             &
          + GY_M_M(PSVM(:,:,:,JSV),PDYY,PDZZ,PDZY)**2 )
    ELSE
      ZFLX(:,:,:) =  ZCSV / ZCSVD * PLM(:,:,:) * PLEPS(:,:,:) *   &
            GX_M_M(PSVM(:,:,:,JSV),PDXX,PDZZ,PDZX)**2
    END IF
    CALL LES_MEAN_SUBGRID( -2.*ZCSVD*SQRT(PTKEM)*ZFLX/PLEPS, &
                           X_LES_SUBGRID_DISS_Sv2(:,:,:,JSV), .TRUE. )
    CALL LES_MEAN_SUBGRID( MZF(PWM)*ZFLX, X_LES_RES_W_SBG_Sv2(:,:,:,JSV), .TRUE. )
  END IF
  !
  ! covariance SvThv
  !
  IF (LLES_CALL) THEN
    ZA(:,:,:)   =  ETHETA(KRR,KRRI,PTHLM,PRM,PLOCPEXNM,PATHETA,PSRCM)
    IF (.NOT. L2D) THEN
      ZFLX(:,:,:)=  PLM(:,:,:) * PLEPS(:,:,:)                                          &
          *  (  GX_M_M(PTHLM,PDXX,PDZZ,PDZX) * GX_M_M(PSVM(:,:,:,JSV),PDXX,PDZZ,PDZX)  &
              + GY_M_M(PTHLM,PDYY,PDZZ,PDZY) * GY_M_M(PSVM(:,:,:,JSV),PDYY,PDZZ,PDZY)  &
             ) * (XCSHF+ZCSV) / (2.*ZCTSVD)
    ELSE
      ZFLX(:,:,:)=  PLM(:,:,:) * PLEPS(:,:,:)                                          &
              * GX_M_M(PTHLM,PDXX,PDZZ,PDZX) * GX_M_M(PSVM(:,:,:,JSV),PDXX,PDZZ,PDZX)  &
              * (XCSHF+ZCSV) / (2.*ZCTSVD)
    END IF
    CALL LES_MEAN_SUBGRID( ZA*ZFLX, X_LES_SUBGRID_SvThv(:,:,:,JSV) , .TRUE.)
    CALL LES_MEAN_SUBGRID( -XG/PTHVREF/3.*ZA*ZFLX, X_LES_SUBGRID_SvPz(:,:,:,JSV), .TRUE. )
    !
    IF (KRR>=1) THEN
      ZA(:,:,:)   =  EMOIST(KRR,KRRI,PTHLM,PRM,PLOCPEXNM,PAMOIST,PSRCM)
      IF (.NOT. L2D) THEN
        ZFLX(:,:,:)=  PLM(:,:,:) * PLEPS(:,:,:)                                                 &
            *  (  GX_M_M(PRM(:,:,:,1),PDXX,PDZZ,PDZX) * GX_M_M(PSVM(:,:,:,JSV),PDXX,PDZZ,PDZX)  &
                + GY_M_M(PRM(:,:,:,1),PDYY,PDZZ,PDZY) * GY_M_M(PSVM(:,:,:,JSV),PDYY,PDZZ,PDZY)  &
               ) * (XCHF+ZCSV) / (2.*ZCQSVD)
      ELSE
        ZFLX(:,:,:)=  PLM(:,:,:) * PLEPS(:,:,:)                                                 &
                * GX_M_M(PRM(:,:,:,1),PDXX,PDZZ,PDZX) * GX_M_M(PSVM(:,:,:,JSV),PDXX,PDZZ,PDZX)  &
                * (XCHF+ZCSV) / (2.*ZCQSVD)
      END IF
      CALL LES_MEAN_SUBGRID( ZA*ZFLX, X_LES_SUBGRID_SvThv(:,:,:,JSV) , .TRUE.)
      CALL LES_MEAN_SUBGRID( -XG/PTHVREF/3.*ZA*ZFLX, X_LES_SUBGRID_SvPz(:,:,:,JSV), .TRUE. )
    END IF
  END IF
!
END DO    ! end loop JSV
!
CALL SECOND_MNH(ZTIME2)
XTIME_LES = XTIME_LES + ZTIME2 - ZTIME1
!
END SUBROUTINE TURB_HOR_SV_CORR
