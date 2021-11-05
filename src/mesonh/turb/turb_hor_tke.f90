!MNH_LIC Copyright 1994-2020 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!    ####################  
     MODULE MODI_TURB_HOR_TKE
!    ####################  
!
INTERFACE  
!
      SUBROUTINE TURB_HOR_TKE(KSPLT,                                 &
                      PDXX,PDYY,PDZZ,PDZX,PDZY,                      &
                      PINV_PDXX, PINV_PDYY, PINV_PDZZ, PMZM_PRHODJ,  &
                      PK, PRHODJ, PTKEM,                             &
                      PTRH                                           )

!
INTEGER,                  INTENT(IN) :: KSPLT        ! current split index
!
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PDXX, PDYY, PDZZ, PDZX, PDZY 
                                                     ! Metric coefficients
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PK           ! Turbulent diffusion doef.
                                                     ! PK = PLM * SQRT(PTKEM)
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PINV_PDXX    ! 1./PDXX
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PINV_PDYY    ! 1./PDYY
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PINV_PDZZ    ! 1./PDZZ
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PMZM_PRHODJ  ! MZM(PRHODJ)
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PRHODJ       ! density * grid volume
!
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PTKEM    ! TKE at time t- dt
REAL, DIMENSION(:,:,:),   INTENT(OUT)   ::  PTRH     ! horizontal transport of Tke 
!
!
!
END SUBROUTINE TURB_HOR_TKE
!
END INTERFACE
!
END MODULE MODI_TURB_HOR_TKE
!     ################################################################
      SUBROUTINE TURB_HOR_TKE(KSPLT,                                 &
                      PDXX, PDYY, PDZZ,PDZX,PDZY,                    &
                      PINV_PDXX, PINV_PDYY, PINV_PDZZ, PMZM_PRHODJ,  &
                      PK, PRHODJ, PTKEM,                             &
                      PTRH                                           )
!     ################################################################
!
!
!!****  *TURB_HOR_TKE* computes the horizontal turbulant transports of Tke
!!
!!    PURPOSE
!!    -------

!!**  METHOD
!!    ------
!!
!!
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
!!      Joan Cuxart             * INM and Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       Aug 29, 1994
!!                     Mar 07  2001 (V. Masson and J. Stein) new routine
!!                     Nov 06, 2002 (V. Masson) LES budgets
!! --------------------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
!
USE MODD_CONF
USE MODD_CST
USE MODD_CTURB
USE MODD_PARAMETERS
USE MODD_LES
!
!
USE MODI_SHUMAN 
USE MODI_GRADIENT_M
USE MODI_LES_MEAN_SUBGRID
!
USE MODI_SECOND_MNH
!
IMPLICIT NONE
!
!
!*       0.1  declaration of arguments
!
!
INTEGER,                  INTENT(IN) :: KSPLT        ! current split index
!
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PDXX, PDYY, PDZZ, PDZX, PDZY 
                                                     ! Metric coefficients
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PK           ! Turbulent diffusion doef.
                                                     ! PK = PLM * SQRT(PTKEM)
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PINV_PDXX    ! 1./PDXX
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PINV_PDYY    ! 1./PDYY
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PINV_PDZZ    ! 1./PDZZ
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PMZM_PRHODJ  ! MZM(PRHODJ)
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PRHODJ       ! density * grid volume
!
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PTKEM    ! TKE at time t- dt
REAL, DIMENSION(:,:,:),   INTENT(OUT)   ::  PTRH     ! horizontal transport of Tke
!
!
!
!*       0.2  declaration of local variables
!
INTEGER :: IKB, IKU
!
REAL, DIMENSION(SIZE(PDZZ,1),SIZE(PDZZ,2),1+JPVEXT:3+JPVEXT) :: ZCOEFF 
                                    ! coefficients for the uncentred gradient 
                                    ! computation near the ground
!
REAL, DIMENSION(SIZE(PTKEM,1),SIZE(PTKEM,2),SIZE(PTKEM,3)):: ZFLX
!
REAL :: ZTIME1, ZTIME2
! ---------------------------------------------------------------------------
!
!*       1.   PRELIMINARY COMPUTATIONS
!             ------------------------
!
IKB = 1.+JPVEXT
IKU = SIZE(PTKEM,3)
!
!  compute the coefficients for the uncentred gradient computation near the 
!  ground
!
ZCOEFF(:,:,IKB+2)= - PDZZ(:,:,IKB+1) /      &
     ( (PDZZ(:,:,IKB+2)+PDZZ(:,:,IKB+1)) * PDZZ(:,:,IKB+2) )
ZCOEFF(:,:,IKB+1)=   (PDZZ(:,:,IKB+2)+PDZZ(:,:,IKB+1)) /      &
     ( PDZZ(:,:,IKB+1) * PDZZ(:,:,IKB+2) )
ZCOEFF(:,:,IKB)= - (PDZZ(:,:,IKB+2)+2.*PDZZ(:,:,IKB+1)) /      &
     ( (PDZZ(:,:,IKB+2)+PDZZ(:,:,IKB+1)) * PDZZ(:,:,IKB+1) )
!
!--------------------------------------------------------------------
!
!*       2.   horizontal transport of Tke u'e
!             -------------------------------
!
!
ZFLX = -XCET * MXM(PK) * GX_M_U(1,IKU,1,PTKEM,PDXX,PDZZ,PDZX) ! < u'e >
!
! special case near the ground ( uncentred gradient )
!
ZFLX(:,:,IKB) =  ZCOEFF(:,:,IKB+2)*PTKEM(:,:,IKB+2)                     &
               + ZCOEFF(:,:,IKB+1)*PTKEM(:,:,IKB+1)                     &
               + ZCOEFF(:,:,IKB  )*PTKEM(:,:,IKB  )     
!
ZFLX(:,:,IKB:IKB) =                                                      &
   - XCET * MXM( PK(:,:,IKB:IKB) )                           *  (        &
       DXM ( PTKEM(:,:,IKB:IKB) ) * PINV_PDXX(:,:,IKB:IKB)               &
      -MXM ( ZFLX (:,:,IKB:IKB) ) * PINV_PDXX(:,:,IKB:IKB)               &
       * 0.5 * ( PDZX(:,:,IKB+1:IKB+1) + PDZX(:,:,IKB:IKB) )     ) 
!
! extrapolate the fluxes to obtain < u'e > = 0 at the ground
!
ZFLX(:,:,IKB-1) = - ZFLX(:,:,IKB)
!
! let the same flux at IKU-1 and IKU level
!
ZFLX(:,:,IKU) =  ZFLX(:,:,IKU-1)
!
IF (.NOT. LFLAT) THEN
  PTRH =-(  DXF( MXM(PRHODJ) * ZFLX                             * PINV_PDXX)&
          - DZF( PMZM_PRHODJ * MXF( PDZX * MZM(ZFLX*PINV_PDXX)) * PINV_PDZZ)&
         ) /PRHODJ
ELSE
  PTRH =-(  DXF( MXM(PRHODJ) * ZFLX                             * PINV_PDXX)&
         ) /PRHODJ
END IF
!
IF (LLES_CALL .AND. KSPLT==1) THEN
  CALL SECOND_MNH(ZTIME1)
  CALL LES_MEAN_SUBGRID( MXF(ZFLX), X_LES_SUBGRID_UTke ) 
  CALL SECOND_MNH(ZTIME2)
  XTIME_LES = XTIME_LES + ZTIME2 - ZTIME1
END IF
!
!
!--------------------------------------------------------------------
!
!*       3.   horizontal transport of Tke v'e
!             -------------------------------
!
IF (.NOT. L2D) THEN
  ZFLX =-XCET * MYM(PK) * GY_M_V(1,IKU,1,PTKEM,PDYY,PDZZ,PDZY) ! < v'e >
!
! special case near the ground ( uncentred gradient )
!
  ZFLX(:,:,IKB) =  ZCOEFF(:,:,IKB+2)*PTKEM(:,:,IKB+2)                     &
                 + ZCOEFF(:,:,IKB+1)*PTKEM(:,:,IKB+1)                     &
                 + ZCOEFF(:,:,IKB  )*PTKEM(:,:,IKB  )     
!
  ZFLX(:,:,IKB:IKB) =                                                      &
     - XCET * MYM( PK(:,:,IKB:IKB) )                        *  (           &
       DYM ( PTKEM(:,:,IKB:IKB) ) * PINV_PDYY(:,:,IKB:IKB)                 &
     - MYM ( ZFLX (:,:,IKB:IKB) ) * PINV_PDYY(:,:,IKB:IKB)                 &
         * 0.5 * ( PDZY(:,:,IKB+1:IKB+1) + PDZY(:,:,IKB:IKB) )  )
!
!    extrapolate the fluxes to obtain < v'e > = 0 at the ground
!
  ZFLX(:,:,IKB-1) = - ZFLX(:,:,IKB)
!
!   let the same flux at IKU-1 and IKU level
!
  ZFLX(:,:,IKU) =  ZFLX(:,:,IKU-1)
!
! complete the explicit turbulent transport
!
  IF (.NOT. LFLAT) THEN
    PTRH = PTRH - (  DYF( MYM(PRHODJ) * ZFLX                              * PINV_PDYY )  &
                   - DZF( PMZM_PRHODJ * MYF( PDZY * MZM(ZFLX*PINV_PDYY) ) * PINV_PDZZ )  &
                  ) /PRHODJ
  ELSE
    PTRH = PTRH - (  DYF( MYM(PRHODJ) * ZFLX                              * PINV_PDYY )  &
                  ) /PRHODJ
  END IF
!
  IF (LLES_CALL .AND. KSPLT==1) THEN
    CALL SECOND_MNH(ZTIME1)
    CALL LES_MEAN_SUBGRID( MYF(ZFLX), X_LES_SUBGRID_VTke )
    CALL SECOND_MNH(ZTIME2)
    XTIME_LES = XTIME_LES + ZTIME2 - ZTIME1
  END IF
!
END IF
!
!----------------------------------------------------------------------------
!
END SUBROUTINE TURB_HOR_TKE
