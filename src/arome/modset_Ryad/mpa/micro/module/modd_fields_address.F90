!     ######spl
      MODULE MODD_FIELDS_ADDRESS
!     ######################
!
!!****  *MODD_FIELDS_ADDRESS* - declaration of fields adress in arrays, as parameter variables
!!
!!    PURPOSE
!!    -------
!     To share fields adress in arrays
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      Ryad El Khatib   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    24-Aug-2021
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE
!
! Pointer of fields in microphysic species arrays. Microphysics species are
! usually counted as KRR=6 or KRR=7. The extra "zero" adress for potential
! temperature is a trick to improve vectorization when all these fields needs
! the same treatement.
INTEGER, PARAMETER :: & ! pointer of fields in microphysic species arrays :
      & ITH=0,     & ! Potential temperature
      & IRV=1,     & ! Water vapor
      & IRC=2,     & ! Cloud water
      & IRR=3,     & ! Rain water
      & IRI=4,     & ! Pristine ice
      & IRS=5,     & ! Snow/aggregate
      & IRG=6,     & ! Graupel
      & IRH=7        ! Hail
!
END MODULE MODD_FIELDS_ADDRESS
