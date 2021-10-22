!     ######spl
MODULE mode_zsrpun
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK

  USE modd_glo

  IMPLICIT NONE

  !**************************************************************************
  !Purpose: ZSR is used to calculate the amount of water associated with Type
  !         A organic molecules only when zsrflag = 1; when zsrflag = 0 use
  !         Newt and UNIFAC to solve for Aw = RH
  !
  !Revision history: 1. Developed by Betty Pun, AER, Feb 99 under EPRI funding
  !                     for the prototype SOA module.
  !               2. Modified November 99 by Betty Pun, AER, under CARB
  !                  funding to adhere to models-3 coding standards
  !               3. Modified to use ZSR or unifac to calculate water
  !                  associated with organics, as specified by zsrflag
  !                  by Betty Pun, Nov, 99.  A file with xi at given Aw
  !                  is included binsolu.h
  !               4. Rewritten to FORTRAN90 by Alf Grini (alf.grini@cnrm.meteo.fr)
  !
  !***************************************************************************/

CONTAINS

  subroutine ZSRPUN(         &
       AERO                  & !I [ug/m3] guess for aerosol concentrations
       ,RH                   & !I [0-1] relative humidity
       ,DELTALWC             & !O [ug/m3] liquid water content
       ,MOLALBIN             & !I [umol/ug_{water}] molality of binary solutions
       ,NSPIN                & !I [nbr] number of species to take into account
       ,NK                   & !I [nbr] number of sub and main component for main component
       ,MW                   & !I [g/mol] molecular weight of species
       )

    USE modd_binsolu

    IMPLICIT NONE

    !INPUT
    REAL, DIMENSION(:,:), INTENT(IN)      :: AERO     !I [ug/m3] guess for aerosol concentrations
    REAL, DIMENSION(:), INTENT(IN)        :: RH       !I [0-1] relative humidity
    REAL, DIMENSION(:,:), INTENT(IN)      :: MOLALBIN !I ![umol/ug_{water}] molality of binary solutions
    INTEGER, INTENT(IN)                   :: NSPIN    !I [nbr] number of (main) species to take into account
    INTEGER, DIMENSION(:), INTENT(IN)     :: NK       !I [nbr] number of sub and main component for main component
    REAL, DIMENSION(:), INTENT(IN)        :: MW       !I [g/mol] molecular weight of species

    !OUTPUT
    REAL, DIMENSION(:), INTENT(OUT)       :: DELTALWC ![ug/m3] Liquid water content associated with organics

    !LOCAL, AUTOMATIC ARRAYS
    REAL, DIMENSION(SIZE(AERO,1))             :: MOLAL1        ![umol_{aer}/ug_{water}] molality of aerosol in binary water mix @ RH1
    REAL, DIMENSION(SIZE(AERO,1))             :: MOLAL2        ![umol_{aer}/ug_{water}] molality of aerosol in binary water mix @ RH2
    REAL, DIMENSION(SIZE(AERO,1))             :: MOLALW        ![umol_{aer}/ug_{water}] molality of aerosol binary water mix @ RH
    INTEGER, DIMENSION(SIZE(AERO,1))          :: IRH1          ![0-10] lower index for relative humidity
    INTEGER, DIMENSION(SIZE(AERO,1))          :: IRH2          ![1-11] upper index for relative humidity

    !LOCAL, ALLOCATABLE ARRAYS
    REAL, DIMENSION(:,:), ALLOCATABLE         :: AEROMOLEC     ![umol/m3] concentration of main components

    !Small, integers, counters
    INTEGER                                   :: COMP_IDX      ![idx] index for components (ions and main comp)
    INTEGER                                   :: I,J           ![idx] counter for main and sub-components
    INTEGER                                   :: IX            ![idx] conter for physical points

    !Allocate memory
    REAL(KIND=JPRB) :: ZHOOK_HANDLE
    IF (LHOOK) CALL DR_HOOK('MODE_ZSRPUN:ZSRPUN',0,ZHOOK_HANDLE)
    ALLOCATE(AEROMOLEC(SIZE(AERO,1),NSPIN))

    !Initialize water associated with organics
    deltaLWC(:)=0.d0
    AEROMOLEC(:,:)=0.d0

    !Start code
    IF (ZSRFLAG.eq.0)THEN
       stop "ZSRFLAG=0 not implemented yet"

    ELSE  !        zsrflag = 1

       !Get the total moles (umole/m3) of the main components,
       !i.e. sum ions into main components
       COMP_IDX=1
       DO I =1,NSPIN
          AEROMOLEC(:,I)=AERO(:,COMP_IDX)/MW(COMP_IDX)
          COMP_IDX=COMP_IDX+1

          !add up the number of moles of ions
          DO J=2,NK(I)
             AEROMOLEC(:,I)=AEROMOLEC(:,I) &
                  + AERO(:,COMP_IDX)/MW(COMP_IDX)

             COMP_IDX=COMP_IDX+1
          ENDDO
       ENDDO

       !Get the molality of component I in mix with water
       !Note that is molality is low, it means RH is high, since we then
       !have much water. So for RH=100%, molality should be very low
       DO I=1,NSPIN

          DO IX=1,SIZE(AERO,1) !Loop on points

             IRH1(IX)=int(RH(IX)/RHGRAD) +1 !==> transforms lower RH into a 1-11 scale, we know binsolu for 0, 1, 2, 3, ..10

             IRH2(IX)=int(RH(IX)/RHGRAD) +2 !==> transforms upper RH into a 1-11 scale, we know binsolu for 0, 1, 2, 3, ..10

             !Get lower part of RH spectre
             MOLAL1(IX)=MOLALBIN(IRH1(IX),I)

             !Get upper part of RH spectre
             MOLAL2(IX)=MOLALBIN(IRH2(IX),I)
          ENDDO !!Loop on points

          !Interpolate ==> molality at current RH for current component
          MOLALW(:) = MOLAL1(:)       &
               +(MOLAL2(:)-MOLAL1(:))/RHgrad *(RH(:) - RHgrad*DBLE(IRH1(:)-1))


          !Use ZSR equation to add total water to deltaLWC
          deltaLWC(:) = deltaLWC(:) +  AEROMOLEC(:,I)/MOLALW(:)

       ENDDO

    ENDIF

    DEALLOCATE(AEROMOLEC)

  IF (LHOOK) CALL DR_HOOK('MODE_ZSRPUN:ZSRPUN',1,ZHOOK_HANDLE)
  END SUBROUTINE ZSRPUN

END MODULE mode_zsrpun
