!MNH_LIC Copyright 2002-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!      ###########################
       MODULE MODI_TO_ELEC_FIELD_n
!      ###########################
!
INTERFACE
      SUBROUTINE TO_ELEC_FIELD_n(PRT, PSVT, PRHODJ, KTCOUNT, KRR, &
                                 PEFIELDU, PEFIELDV, PEFIELDW, PPHIT)
!
INTEGER,                  INTENT(IN)    :: KTCOUNT ! counter value of the 
                                                   ! model temporal loop
INTEGER,                  INTENT(IN)    :: KRR     ! Number of moist variables
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ  ! Jacobian
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PSVT    ! Scalar variables with
                                                   ! electric charge density
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRT     ! Mixing ratio 
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PEFIELDU !  3 components 
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PEFIELDV !     of the 
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PEFIELDW ! electric field
REAL, DIMENSION(:,:,:),   OPTIONAL, INTENT(OUT) :: PPHIT   ! Electrostatic potential

END SUBROUTINE TO_ELEC_FIELD_n
END INTERFACE
END MODULE MODI_TO_ELEC_FIELD_n
!
!     ###############################################################
      SUBROUTINE TO_ELEC_FIELD_n(PRT, PSVT, PRHODJ, KTCOUNT, KRR, &
                                 PEFIELDU, PEFIELDV, PEFIELDW, PPHIT)
!     ###############################################################
!
!
!!****  * -  compute the electric field
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute...
!!
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!      None
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      C. Barthe, G. Molinie, J.-P. Pinty      *Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    2002
!!      C. Barthe   06/11/09   update to version 4.8.1
!!      M. Chong    26/01/10   Add Small ions 
!  P. Wautelet 20/05/2019: add name argument to ADDnFIELD_ll + new ADD4DFIELD_ll subroutine
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_REF_n, ONLY : XRHODREF
USE MODD_PARAMETERS, ONLY : JPVEXT
USE MODD_RAIN_ICE_DESCR, ONLY : XRTMIN
USE MODD_ELEC_DESCR, ONLY : XRELAX_ELEC, XECHARGE 
USE MODD_ELEC_n, ONLY : XESOURCEFW
!
USE MODI_ELEC_FIELD_n
!
USE MODE_ll
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER,                  INTENT(IN)    :: KTCOUNT ! counter value of the 
                                                   ! model temporal loop
INTEGER,                  INTENT(IN)    :: KRR     ! Number of moist variables
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ  ! Jacobian
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PSVT    ! Scalar variables with
                                                   ! electric charge density
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRT     ! Mixing ratio 
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PEFIELDU !  3 components 
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PEFIELDV !     of the 
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PEFIELDW ! electric field
REAL, DIMENSION(:,:,:),   OPTIONAL, INTENT(OUT) :: PPHIT   ! Electrostatic potential
!
!
!*       0.2   Declarations of local variables :
!
REAL,    DIMENSION(:,:,:),   ALLOCATABLE :: ZW ! work array
!
INTEGER :: IIB     ! Define 
INTEGER :: IIE     !   the 
INTEGER :: IJB     ! physical  
INTEGER :: IJE     !  domain
INTEGER :: IKB     ! 
INTEGER :: IKE     ! 
INTEGER :: IIU, IJU, IKU 
INTEGER :: II
INTEGER :: IINFO_ll  
!
TYPE(LIST_ll), POINTER :: TZFIELDS_ll ! list of fields to exchange
!
!
!-------------------------------------------------------------------------------
!
!*       1.     COMPUTE THE LOOP BOUNDS
!   	        -----------------------
!
NULLIFY(TZFIELDS_ll)
!
! Compute loop bounds
!
CALL GET_PHYSICAL_ll(IIB,IJB,IIE,IJE)
CALL GET_DIM_EXT_ll('B',IIU,IJU)
!
IKB = 1 + JPVEXT
IKU = SIZE(XESOURCEFW,3)
IKE = IKU - JPVEXT
!
! allocations
!
ALLOCATE(ZW(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3)))
ZW(:,:,:) = 0.
!
!
!-------------------------------------------------------------------------------
!
!*       2.     TRANSFORM PSVT from C/kg INTO C/m3 and SUM
!   	        ----------------------------------
!
DO II = 1, KRR+1
  ZW(:,:,:) = ZW(:,:,:) + PSVT(:,:,:,II) * XRHODREF(:,:,:)
END DO
!
!-------------------------------------------------------------------------------
!
!*       3.     BOUNDARY CONDITIONS
!   	        -------------------
!
ZW(:,:,1:IKB-1)        = 0.0 ! Setup to neutralize the computation on the
                             ! first ligne of the tridiagonal system starting
                             ! at IKB-1 
ZW(:,:,IKE:IKE+JPVEXT) = XESOURCEFW(:,:,IKE:IKE+JPVEXT)
!
CALL ADD3DFIELD_ll( TZFIELDS_ll, ZW, 'TO_ELEC_FIELD_n::ZW' )
CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
CALL CLEANLIST_ll(TZFIELDS_ll)
!
!
!-------------------------------------------------------------------------------
!
!*       4.     COMPUTE THE ELECTRIC FIELD
!   	        --------------------------
!
IF (PRESENT(PPHIT)) THEN 
  CALL ELEC_FIELD_n (ZW, KTCOUNT, XRELAX_ELEC, PRHODJ, &
                     PEFIELDU, PEFIELDV, PEFIELDW, PPHIT) 
ELSE 
  CALL ELEC_FIELD_n (ZW, KTCOUNT, XRELAX_ELEC, PRHODJ, &
                     PEFIELDU, PEFIELDV, PEFIELDW) 
ENDIF
!
DEALLOCATE(ZW)
!
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE TO_ELEC_FIELD_n

