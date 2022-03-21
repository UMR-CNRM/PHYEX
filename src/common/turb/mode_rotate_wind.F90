!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!    #######################  
     MODULE MODE_ROTATE_WIND
!    #######################
IMPLICIT NONE
CONTAINS
!     ###########################################################
      SUBROUTINE ROTATE_WIND(PU,PV,PW,                          &
                             PDIRCOSXW, PDIRCOSYW, PDIRCOSZW,   &
                             PCOSSLOPE,PSINSLOPE,               &
                             PDXX,PDYY,PDZZ,                    &
                             PUSLOPE,PVSLOPE                    )
!     ###########################################################
!
!
!!****  *ROTATE_WIND* - computes the wind components along the maximum slope 
!!               direction and its normal direction in the first mass level.
!!
!!    PURPOSE
!!    -------
!!**** 
!        The purpose of this routine is to compute the wind component parallel 
!     to the orography at the first mass level. The exact location where these
!     components are computed is the point of intersection between the normal 
!     to the orography and the first mass-level hyper-plane at PDZZ(:,:,IKB)/2 
!        
!!**  METHOD
!!    ------
!!       The values of the 3 cartesian components of the wind are determined
!!    by a bilinear interpolation between the 4 nearest points in the first 
!!    mass-level hyper-plane. These points are found according to the signs of 
!!    the slopes' sinus and cosinus. For each direction of interpolation, the 
!!    two different localizations (mass or flux grids) are used to avoid 
!!    lateral boundary problems.  
!!       Then, the rotation is performed for the wind components. The rotation 
!!    angle is the angle between the x axe and the maximum slope direction 
!!    defined by the slope vector (dZs/dx , dZs/dy).
!!        Finally, the horizontal components are set at the marginal points 
!!    according to cyclic boundary conditions because this is the only case
!!    where these points can be considered.
!!
!!    EXTERNAL
!!    --------
!!       NONE
!!
!!    IMPLICIT ARGUMENTS 
!!    ------------------
!!
!!       MODD_CONF      : L2D   switch for 2D model version
!!
!!
!!    REFERENCE
!!    ---------
!!      Book 1 of documentation (Chapter: Turbulence)
!!
!!    AUTHOR
!!    ------
!!      Joel Stein              * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original         14/11/95
!!      Modifications:   15/05/96, (N. wood)
!!                                 take into account no slip conditions 
!!                                 at the surface
!!                       14/02/01  (V. Masson)
!!                                 Slip condition at the surface restored
!!
!! --------------------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
USE MODD_PARAMETERS, ONLY: JPVEXT
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PU,PV,PW        ! cartesian components
                                 ! of the wind
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PDIRCOSXW, PDIRCOSYW, PDIRCOSZW
! Director Cosinus along x, y and z directions at surface w-point
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PCOSSLOPE       ! cosinus of the angle 
                                 ! between i and the slope vector
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PSINSLOPE       ! sinus of the angle 
                                 ! between i and the slope vector
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PDXX, PDYY, PDZZ
                                 ! Metric coefficients
REAL, DIMENSION(:,:),   INTENT(OUT)  ::  PUSLOPE         ! wind component along 
                                 ! the maximum slope direction
REAL, DIMENSION(:,:),   INTENT(OUT)  ::  PVSLOPE         ! wind component along
                                 !  the direction normal to the maximum slope one
!
!-------------------------------------------------------------------------------
!
!       0.2  declaration of local variables
!
INTEGER, DIMENSION(SIZE(PDIRCOSXW,1),SIZE(PDIRCOSXW,2)) :: ILOC,JLOC
              ! shift index to find the 4 nearest points in x and y directions
REAL,    DIMENSION(SIZE(PDIRCOSXW,1),SIZE(PDIRCOSXW,2)) :: ZCOEFF,ZCOEFM,     &
              ! interpolation weigths for flux and mass locations
                                                           ZUINT,ZVINT,ZWINT, &
              ! intermediate values of the cartesian components after x interp.
                                                           ZUFIN,ZVFIN,ZWFIN, &
              ! final values of the cartesian components after the 2 interp.
                                                           ZWGROUND
              ! vertical velocity at the surface                                                            
INTEGER     :: IIB,IIE,IJB,IJE,IKB
              ! index values for the Beginning or the End of the physical 
              ! domain in x,y and z directions
INTEGER     :: IIU,IJU
              ! arrays' sizes for i and j indices
INTEGER     :: JI,JJ
!      
!----------------------------------------------------------------------------
!
!*      1.    PRELIMINARIES
!             -------------
!
PUSLOPE=0.
PVSLOPE=0.
!
IIB = 2
IJB = 2
IIU = SIZE(PU,1)
IJU = SIZE(PU,2)
IIE = IIU - 1
IJE = IJU - 1
IKB = 1+JPVEXT
!
ZWGROUND(:,:) = PW(:,:,IKB)
!
!*      2.    INTERPOLATE THE CARTESIAN COMPONENTS
!             ------------------------------------
!
ILOC(:,:)=NINT(SIGN(1.,-PCOSSLOPE(:,:)))
JLOC(:,:)=NINT(SIGN(1.,-PSINSLOPE(:,:)))
!
! interpolation in x direction
!
DO JJ = 1,IJU
  DO JI = IIB,IIE 
    ZCOEFF(JI,JJ) =                                                  &
      (0.5*PDXX(JI,JJ,IKB) + 0.5*PDZZ(JI,JJ,IKB)*PDIRCOSXW(JI,JJ) )  & 
      * 2. / (PDXX(JI,JJ,IKB)+PDXX(JI+1,JJ,IKB))
    ZUINT(JI,JJ) = ZCOEFF(JI,JJ)      * PU(JI+1,JJ,IKB)  +           &
                   (1.-ZCOEFF(JI,JJ)) * PU(JI,JJ,IKB)
    !
    ZCOEFM(JI,JJ) = 1. - 0.5 * PDZZ(JI,JJ,IKB) * ABS(PDIRCOSXW(JI,JJ))      & 
                             / PDXX(JI+(ILOC(JI,JJ)+1)/2,JJ,IKB)
    ZVINT(JI,JJ) = ZCOEFM(JI,JJ)      * PV(JI,JJ,IKB)              +        &
                   (1.-ZCOEFM(JI,JJ)) * PV(JI+ILOC(JI,JJ),JJ,IKB)
    !
    ZWINT(JI,JJ) = ZCOEFM(JI,JJ)  * (PW(JI,JJ,IKB+1)+ZWGROUND(JI,JJ)) * 0.5    &
              + (1.-ZCOEFM(JI,JJ))                                             &
               *(PW(JI+ILOC(JI,JJ),JJ,IKB+1)+ZWGROUND(JI+ILOC(JI,JJ),JJ)) * 0.5
  END DO
END DO
!
! interpolation in y direction
!
DO JJ = IJB,IJE
  DO JI = IIB,IIE
    ZCOEFF(JI,JJ) =                                                     &
      (0.5*PDYY(JI,JJ,IKB) + 0.5*PDZZ(JI,JJ,IKB)*PDIRCOSYW(JI,JJ) )     & 
      * 2. / (PDYY(JI,JJ,IKB)+PDYY(JI+1,JJ,IKB))
    ZVFIN(JI,JJ) = ZCOEFF(JI,JJ)      * ZVINT(JI,JJ+1)  +               &
                   (1.-ZCOEFF(JI,JJ)) * ZVINT(JI,JJ)
    !
    ZCOEFM(JI,JJ) = 1. - 0.5 * PDZZ(JI,JJ,IKB) * ABS(PDIRCOSYW(JI,JJ))   & 
                             / PDYY(JI,JJ+(JLOC(JI,JJ)+1)/2,IKB)
    ZUFIN(JI,JJ) = ZCOEFM(JI,JJ)      * ZUINT(JI,JJ)                +    &
                   (1.-ZCOEFM(JI,JJ)) * ZUINT(JI,JJ+JLOC(JI,JJ))
    ZWFIN(JI,JJ) = ZCOEFM(JI,JJ)      * ZWINT(JI,JJ)                +    &
                   (1.-ZCOEFM(JI,JJ)) * ZWINT(JI,JJ+JLOC(JI,JJ))
  END DO
END DO
!
!*      3.    ROTATE THE WIND
!             ---------------
!
!
DO JJ = IJB,IJE 
  DO JI = IIB,IIE
    PUSLOPE(JI,JJ) = PCOSSLOPE(JI,JJ) * PDIRCOSZW(JI,JJ) * ZUFIN(JI,JJ) +   &
                     PSINSLOPE(JI,JJ) * PDIRCOSZW(JI,JJ) * ZVFIN(JI,JJ) +   &
                            SQRT(1.-PDIRCOSZW(JI,JJ)**2) * ZWFIN(JI,JJ)
    !              
    PVSLOPE(JI,JJ) =-PSINSLOPE(JI,JJ)                    * ZUFIN(JI,JJ) +   &
                     PCOSSLOPE(JI,JJ)                    * ZVFIN(JI,JJ)
    !
  END DO
END DO
!
!
!
!----------------------------------------------------------------------------
!
END SUBROUTINE ROTATE_WIND
END MODULE MODE_ROTATE_WIND
