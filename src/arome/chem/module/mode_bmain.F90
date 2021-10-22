!     ######spl
MODULE MODE_BMAIN
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK

  USE MODD_GLO     !Global definitions

  IMPLICIT NONE
  
CONTAINS

  SUBROUTINE BMAIN(      &
       CB                & !I [ug/m3] total concentration 
       ,GAS              & !I/O [ug/m3] gas phase concentrations
       ,AERO             & !I/O [ug/m3] 
       ,PAOM             & !I [ug/m3] primary aerosol organic matter
       ,TEMPK            & !I [K] Temperature
       )

    USE mode_typeb

    IMPLICIT NONE
    
    !INPUT
    REAL, INTENT(IN), DIMENSION(:,:)    :: CB     ![ug/m3] total concentration
    REAL, INTENT(IN), DIMENSION(:)      :: PAOM   ![ug/m3] Primary aerosol organic matter
    REAL, INTENT(IN), DIMENSION(:)      :: TEMPK  ![K] temperature
    
    !OUTPUT
    REAL, INTENT(OUT), DIMENSION(:,:) :: GAS    ![ug/m3] gas concentratin
    REAL, INTENT(OUT), DIMENSION(:,:) :: AERO   ![ug/m3] aerosol conc

    !LOCAL VARIABLES
    REAL, DIMENSION(SIZE(GAS,1))             :: GUESSSOA  ![ug/m3] guess for soa
    REAL                                     :: GUESSFRC  ![frc] guess for fraction of total in aerosol phase

    INTEGER                                  :: II        ![idx] counters

    !Set the total concentration of type B species
    !CB(:,:) = GAS(:,:) + AERO (:,:)

    REAL(KIND=JPRB) :: ZHOOK_HANDLE
    IF (LHOOK) CALL DR_HOOK('MODE_BMAIN:BMAIN',0,ZHOOK_HANDLE)
    GUESSSOA(:)=0.0 !Used later if SOA >> PAOM
    
    DO II=1,NBSPB
       !Get initial guesses for aerosol part
       !If component has high vapor pressure, compared to POA
       !==>guess 30% of total
       !If component has low vapor pressure, compared to POA
       !==> guess 90% of total
       !This guessing is not really OK for cases where SOA is the 
       !most important solvent species. In those cases, we get longer iterations
       WHERE(VPB(II)/(PAOM(:)+1.d-20) > VPCRIT ) 
          AERO(:,II) = 0.3 * CB(:,II)
          GUESSSOA(:)=GUESSSOA(:)+AERO(:,II)
       ELSEWHERE
          AERO(:,II) = 0.9 * CB(:,II)
          GUESSSOA(:)=GUESSSOA(:)+AERO(:,II)
       ENDWHERE
    ENDDO
    
    !Special calculation for the cases where the SOA part is a lot greater 
    !than the POA part: Need special guess:
    DO II=1,NBSPB       
       IF(VPB(II) < VPCRIT) THEN
          GUESSFRC=0.9d0
       ELSE
          GUESSFRC=0.3d0
       ENDIF
       
       !Since the above test did not work well if POA >> SOA
       !We have to check if the GUESSSOA variable is large compared to POA
       WHERE(10.0*PAOM(:) < GUESSSOA(:))
          
          AERO(:,II) = GUESSFRC * CB(:,II)
          
       ENDWHERE!Check on total SOA larger than POA
       
    ENDDO  !Loop on SOA components

    !Limit the new aerosol values to reasonable values
    !WHERE (AERO(:,:) < 0.0)
    !   AERO(:,:)=0.d0
    !ELSEWHERE(AERO(:,:) > CB(:,:))
    !   AERO(:,:) = CB(:,:)
    !ENDWHERE
    
    !Get the gas phase concentration
    GAS(:,:)=CB(:,:) - AERO (:,:)

    CALL TYPEB(         &
         PAOM           & !I [ug/m3] primary organic aerosol
         ,AERO          & !I/O [ug/m3] aerosol concentrations
         ,GAS           & !I/O [ug/m3] gas phase concentrations
         ,CB            & !I [ug/m3] total concentrations
         ,TEMPK         & !I [K] temperature
         )                
    
  IF (LHOOK) CALL DR_HOOK('MODE_BMAIN:BMAIN',1,ZHOOK_HANDLE)
  END SUBROUTINE BMAIN
  
END MODULE MODE_BMAIN
