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
      SUBROUTINE ROTATE_WIND(D,PU,PV,PW,                          &
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
USE MODD_DIMPHYEX,   ONLY: DIMPHYEX_t
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
TYPE(DIMPHYEX_t),       INTENT(IN) :: D
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PU,PV,PW        ! cartesian components
                                 ! of the wind
REAL, DIMENSION(D%NIT,D%NJT),   INTENT(IN)   ::  PDIRCOSXW, PDIRCOSYW, PDIRCOSZW
! Director Cosinus along x, y and z directions at surface w-point
REAL, DIMENSION(D%NIT,D%NJT),   INTENT(IN)   ::  PCOSSLOPE       ! cosinus of the angle 
                                 ! between i and the slope vector
REAL, DIMENSION(D%NIT,D%NJT),   INTENT(IN)   ::  PSINSLOPE       ! sinus of the angle 
                                 ! between i and the slope vector
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PDXX, PDYY, PDZZ
                                 ! Metric coefficients
REAL, DIMENSION(D%NIT,D%NJT),   INTENT(OUT)  ::  PUSLOPE         ! wind component along 
                                 ! the maximum slope direction
REAL, DIMENSION(D%NIT,D%NJT),   INTENT(OUT)  ::  PVSLOPE         ! wind component along
                                 !  the direction normal to the maximum slope one
!
!-------------------------------------------------------------------------------
!
!       0.2  declaration of local variables
!
INTEGER, DIMENSION(D%NIT,D%NJT) :: ILOC,JLOC
              ! shift index to find the 4 nearest points in x and y directions
REAL,    DIMENSION(D%NIT,D%NJT) :: ZCOEFF,ZCOEFM,     &
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
PUSLOPE(:,:)=0.
PVSLOPE(:,:)=0.
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
!
!     ##############################################
      SUBROUTINE UPDATE_ROTATE_WIND(D,PUSLOPE,PVSLOPE,HLBCX,HLBCY)
!     ##############################################
!!
!!****  *UPDATE_ROTATE_WIND* routine to set rotate wind values at the border
!
!!    AUTHOR
!!    ------
!!
!!     P Jabouille   *CNRM METEO-FRANCE
!!
!!    MODIFICATIONS
!!    -------------
!!      Original   24/06/99
!!      J.Escobar 21/03/2013: for HALOK comment all NHALO=1 test
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODE_ll, ONLY: ADD2DFIELD_ll, UPDATE_HALO_ll, CLEANLIST_ll, &
                   LWEST_ll, LEAST_ll, LSOUTH_ll, LNORTH_ll
USE MODD_ARGSLIST_ll, ONLY : LIST_ll
USE MODD_DIMPHYEX,   ONLY: DIMPHYEX_t
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(DIMPHYEX_t),       INTENT(IN) :: D
CHARACTER(LEN=4),DIMENSION(2),INTENT(IN):: HLBCX, HLBCY  ! X- and Y-direc LBC
REAL, DIMENSION(D%NIT,D%NJT), INTENT(INOUT) :: PUSLOPE,PVSLOPE
!
TYPE(LIST_ll), POINTER :: TZFIELDS_ll  ! list of fields to exchange
INTEGER                :: IINFO_ll       ! return code of parallel routine
!
! tangential surface fluxes in the axes following the orography
!
!*        1  PROLOGUE
!
NULLIFY(TZFIELDS_ll)
!
!         2 Update halo if necessary
!
!!$IF (NHALO == 1) THEN
  CALL ADD2DFIELD_ll( TZFIELDS_ll, PUSLOPE, 'UPDATE_ROTATE_WIND::PUSLOPE' )
  CALL ADD2DFIELD_ll( TZFIELDS_ll, PVSLOPE, 'UPDATE_ROTATE_WIND::PVSLOPE' )
  CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
  CALL CLEANLIST_ll(TZFIELDS_ll)
!!$ENDIF
!
!        3 Boundary conditions for non cyclic case
!
IF ( HLBCX(1) /= "CYCL" .AND. LWEST_ll()) THEN
  PUSLOPE(D%NIB-1,:)=PUSLOPE(D%NIB,:)
  PVSLOPE(D%NIB-1,:)=PVSLOPE(D%NIB,:)
END IF
IF ( HLBCX(2) /= "CYCL" .AND. LEAST_ll()) THEN
  PUSLOPE(D%NIE+1,:)=PUSLOPE(D%NIE,:)
  PVSLOPE(D%NIE+1,:)=PVSLOPE(D%NIE,:)
END IF
IF ( HLBCY(1) /= "CYCL" .AND. LSOUTH_ll()) THEN
  PUSLOPE(:,D%NJB-1)=PUSLOPE(:,D%NJB)
  PVSLOPE(:,D%NJB-1)=PVSLOPE(:,D%NJB)
END IF
IF(  HLBCY(2) /= "CYCL" .AND. LNORTH_ll()) THEN
  PUSLOPE(:,D%NJE+1)=PUSLOPE(:,D%NJE)
  PVSLOPE(:,D%NJE+1)=PVSLOPE(:,D%NJE)
END IF
!
END SUBROUTINE UPDATE_ROTATE_WIND
END MODULE MODE_ROTATE_WIND
