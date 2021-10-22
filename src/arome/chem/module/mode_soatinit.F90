!     ######spl
module mode_soatinit
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK

  use modd_glo

  implicit none

  !Purpose: Contains subroutines to initialize all parameters which are
  !only Temperature-dependent

contains

  !**********************************
  subroutine VP_GET(     &
       TEMPK             & !I [K] Temperature
       ,VP               & !I [torr] saturation vapor pressures
       )

    !Purpose: Find the saturation vapor pressures corrected for temperature
    implicit none

    !Input
    REAL, DIMENSION(:), INTENT(IN)     :: TEMPK !I [K] Temperature

    !Output
    REAL, DIMENSION(:,:), INTENT(OUT)  :: VP    !O [torr] saturation vapor pressures

    !Local
    REAL, DIMENSION(SIZE(VP,2))            :: DELCP
    REAL, DIMENSION(SIZE(VP,2))            :: DELSB
    REAL, DIMENSION(SIZE(VP,1),SIZE(VP,2)) :: TERM1
    REAL, DIMENSION(SIZE(VP,1),SIZE(VP,2)) :: TERM2

    INTEGER                            :: I     ![idx] counter for component

    !The following is cut and paste from Griffin's code
    !I have to check the units here!

    REAL(KIND=JPRB) :: ZHOOK_HANDLE
    IF (LHOOK) CALL DR_HOOK('MODE_SOATINIT:VP_GET',0,ZHOOK_HANDLE)
    DO I = 1,NBSP

       !Calculate DELSB
       DELSB(I)= 86.0 + 0.4*TAUVP(I)+1421.*HBN(I)

       !Calculate DELCP
       DELCP(I) = -90.-2.1*TAUVP(I)

       !Calculate Intermediate terms 1
       TERM1(:,I) = DELSB(I)*(TEMPK(:)-TBOIL(I)) &
            /19.1/TEMPK(:)

       !Calcualte intermediate term 2
       TERM2(:,I) = ((TBOIL(I)-TEMPK(:))/TEMPK(:) &
            -LOG(TBOIL(I)/TEMPK(:)))  &  
            *DELCP(I)/19.1

       !Get the temperature-dependent 
       VP(:,I) = 760.*(10.**(TERM1(:,I)+TERM2(:,I)))

    ENDDO

    !Tuning: Correct vp so that they match chamber data
    VP(:,1) = VP(:,1)/1.5
    VP(:,2) = VP(:,2)/1.4
    VP(:,3) = VP(:,3)/2.4
    VP(:,5) = VP(:,5)/66.
    VP(:,8) = VP(:,8)/33.
    VP(:,9) = VP(:,9)/3.3
    VP(:,10) = VP(:,10)/15.


  IF (LHOOK) CALL DR_HOOK('MODE_SOATINIT:VP_GET',1,ZHOOK_HANDLE)
  end subroutine VP_GET

  !*********************************

  subroutine AQCONST_GET(      &
       TEMPK                   & !I [K] temperature
       ,K                      & !O [units] Henry's law and dissociation constants
       ,NK                     & !I [nbr] Number of total and sub-components for one acid
       )

    !Purpose: Calculate 
    implicit none
    !Input
    REAL, INTENT(IN), DIMENSION(:)            :: TEMPK     !I [K] temperature
    INTEGER, INTENT(IN), DIMENSION(:)         :: NK        !I [nbr] Number of total and sub-components for one acid

    !Output
    REAL, INTENT(OUT), DIMENSION(:,:)         :: K         !O [units] temperature dependent value of K_298
    
    !Local 
    INTEGER                            :: I        ![idx] counter for components
    INTEGER                            :: COMP_IDX ![idx] pointer to right component
    INTEGER                            :: J        ![idx] counter for sub components (ions)

    REAL(KIND=JPRB) :: ZHOOK_HANDLE
    IF (LHOOK) CALL DR_HOOK('MODE_SOATINIT:AQCONST_GET',0,ZHOOK_HANDLE)
    COMP_IDX=1
    DO I=1,NBSP
       
       !Adjust Henry's law constant
       K(:,COMP_IDX) = K_298(COMP_IDX)     &
            *EXP(15.*((1./TEMPK(:))-(1./298.15))/0.001987)
       !Prepare for next component
       COMP_IDX=COMP_IDX+1

       !Calculate dissociation constants of sub-components (ions)
       DO J=2,NK(I)
          K(:,COMP_IDX) = K_298(COMP_IDX)  &
               *EXP(0.5*((1./TEMPK(:))-(1./298.15))/0.001987)
          !Prepare for next component
          COMP_IDX=COMP_IDX + 1
       ENDDO

    ENDDO !loop on I
    
  IF (LHOOK) CALL DR_HOOK('MODE_SOATINIT:AQCONST_GET',1,ZHOOK_HANDLE)
  end subroutine AQCONST_GET
     
  !*************************************

  subroutine SI_GET(           &
       TEMPK                    & !I [K] temperature
       ,A                       & !I [units?] term in UNIFAC parameterization
       ,SI                      & !O [units?] term in UNIFAC parameterization
       ,NFUNC                   & !I [nbr] number of functional group in mixture
       )

    !Purpose: Calculate the term in the unifac method which
    !is dependent only on temperature. This term includes
    !exponential of a large matrix, so it presumably takes 
    !up a lot of computer time.
    implicit none
    
    !Input
    REAL, DIMENSION(:), INTENT(IN)       :: TEMPK     !I [K] temperature
    REAL, DIMENSION(:,:), INTENT(IN)     :: A         !I [?] Term in UNIFAC parameterization
    INTEGER, INTENT(IN)                  :: NFUNC     !I [nbr] number of functional groups in mixture

    !Output
    REAL, DIMENSION(:,:,:), INTENT(OUT)  :: SI        !O [?] term in UNIFAC parameterization

    !Local
    REAL, DIMENSION(SIZE(TEMPK))         :: TEMPKINV  ![1/K] inverse of temperature
    INTEGER                              :: J1, J2    ![idx] counter for functional groups
    

    !Get values of SI (temperature dependent), don't need to be in the iteration!
    !eqn 9 in MaP05
    REAL(KIND=JPRB) :: ZHOOK_HANDLE
    IF (LHOOK) CALL DR_HOOK('MODE_SOATINIT:SI_GET',0,ZHOOK_HANDLE)
    TEMPKINV(:) = 1.d0/TEMPK(:)
    DO J1=1,NFUNC
       DO J2=1,NFUNC
          SI(:,J1,J2) = EXP(-A(J1,J2)*TEMPKINV(:))
       ENDDO
    ENDDO

  IF (LHOOK) CALL DR_HOOK('MODE_SOATINIT:SI_GET',1,ZHOOK_HANDLE)
  end subroutine SI_GET

end module mode_soatinit
