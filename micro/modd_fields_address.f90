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
! Pointers for tendency arrays
! Tendencies are computed either directly as a tendency or as a mixing ratio change that is transformed, afterwards, in a tendency
! The second type is suffixed by _MR
! Some final tendencies can have two contributions (one from a tendency, one from a mixing ratio change).
! In the following list, order matters:
! - first are the normal tendencies directly computed as tendencies
! - second are the tendencies computed only from a mixing ratio change
! - third are the indexes used to designate the mising ratio change part of double-contribution tendencies
INTEGER, PARAMETER :: IBUNUM=47,    & ! Total number
                      IBUNUM_MR=3,  & ! Number of tendencies computed only from a mixing ratio change
                      IBUNUM_EXTRA=2  ! Extra terms
INTEGER, PARAMETER :: &
                    !normal tendencies directly computed as tendencies
                    & IRCHONI=1,     & ! Homogeneous nucleation
                    & IRVDEPS=2,     & ! Deposition on r_s,
                    & IRIAGGS=3,     & ! Aggregation on r_s
                    & IRIAUTS=4,     & ! Autoconversion of r_i for r_s production
                    & IRVDEPG=5,     & ! Deposition on r_g
                    & IRCAUTR=6,     & ! Autoconversion of r_c for r_r production
                    & IRCACCR=7,     & ! Accretion of r_c for r_r production
                    & IRREVAV=8,     & ! Evaporation of r_r
                    & IRCBERI=9,     & ! Bergeron-Findeisen effect
                    & IRHMLTR=10,    & ! Melting of the hailstones
                    & IRSMLTG=11,    & ! Conversion-Melting of the aggregates
                    & IRCMLTSR=12,   & ! Cloud droplet collection onto aggregates by positive temperature
                    & IRRACCSS=13, IRRACCSG=14, IRSACCRG=15, & ! Rain accretion onto the aggregates
                    & IRCRIMSS=16, IRCRIMSG=17, IRSRIMCG=18, & ! Cloud droplet riming of the aggregates
                    & IRICFRRG=19, IRRCFRIG=20, IRICFRR=21,  & ! Rain contact freezing
                    & IRCWETG=22,  IRIWETG=23,  IRRWETG=24,  IRSWETG=25, &  ! Graupel wet growth
                    & IRCDRYG=26,  IRIDRYG=27,  IRRDRYG=28,  IRSDRYG=29, &  ! Graupel dry growth
                    & IRWETGH=30,    & ! Conversion of graupel into hail
                    & IRGMLTR=31,    & ! Melting of the graupel
                    & IRCWETH=32,  IRIWETH=33,  IRSWETH=34,  IRGWETH=35,  IRRWETH=36, & ! Dry growth of hailstone
                    & IRCDRYH=37,  IRIDRYH=38,  IRSDRYH=39,  IRRDRYH=40,  IRGDRYH=41, & ! Wet growth of hailstone
                    & IRDRYHG=42,    &

                    !tendencies computed only with a mixing ratio change
                    & IRVHENI_MR=43, & ! heterogeneous nucleation mixing ratio change
                    & IRRHONG_MR=44, & ! Spontaneous freezing mixing ratio change
                    & IRIMLTC_MR=45, & ! Cloud ice melting mixing ratio change

                    !Extra term computed as a mixing ratio change, to be added to other term
                    & IRSRIMCG_MR=46,& ! Cloud droplet riming of the aggregates
                    & IRWETGH_MR=47    ! Conversion of graupel into hail
INTEGER, PARAMETER, DIMENSION(IBUNUM-IBUNUM_EXTRA+1:IBUNUM) :: IBUEXTRAIND=(/18, 30/)
!
END MODULE MODD_FIELDS_ADDRESS
