!     ######spl
         MODULE MODE_DUSTOPT
         USE PARKIND1, ONLY : JPRB
         USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!        ###################
!
!!
!!    PURPOSE
!!    -------
!!
!!    AUTHOR
!!    ------
!!      Alf Grini (CNRM/GMEI)
!!
!!    MODIFICATIONS
!!    -------------
!!     29 Juin 2007 P. TULET (CNRM/GMEI) / M. Mallet (LA)
!!                  Introduction of new aerosols radiative index updated to the
!!                  west african dusts observed during the AMMA campaign
!!

  IMPLICIT NONE
  PUBLIC
  PRIVATE :: DUSTOPT_LKT

CONTAINS
  
  !****************************************************
  SUBROUTINE DUSTOPT_GET(     &
       PSVT                   & !I [moments/molec_{air}] Transported moments of dust
       ,PZZ                   & !I [m] height of layers
       ,PRHODREF              & !I [kg/m3] density of air
       ,PPIZA_WVL             & !O [-] single scattering albedo of dust layer for all SW wavelengths
       ,PCGA_WVL              & !O [-] assymetry factor for dust layer for all SW wavelengths
       ,PTAUREL_WVL           & !O [-] opt.depth/opt.depth(550) for dust layer for all SW wvl 
       ,PTAU550               & !O [-] opt.depth at 550nm for all dust layer
       ,KSWB                  & !I [nbr] number of shortwave bands
       )
    

    USE MODE_DUST_PSD   !Conversion procedures from moments to radius, ,number, mass and sigma
    USE MODD_DUST, ONLY : NMODE_DST
    IMPLICIT NONE
    
    !INPUT
    REAL, DIMENSION(:,:,:,:),INTENT(IN)      :: PSVT       !I [moments/molec_{air}] transported moments of dust
    REAL, DIMENSION(:,:,:),INTENT(IN)        :: PZZ        !I [m] height of layers
    REAL, DIMENSION(:,:,:),INTENT(IN)        :: PRHODREF   !I [kg/m3] density of air
    INTEGER, INTENT(IN)                      :: KSWB       !I [nbr] number of shortwave wavelengths
    REAL, PARAMETER                          :: EPSILON=1.e-8 !a very low number for optical depth in a layer

    !OUTPUT
    REAL, DIMENSION(:,:,:,:),INTENT(INOUT)     :: PPIZA_WVL   !O [-] single scattering albedo of dust layer for all SW wavelengths
    REAL, DIMENSION(:,:,:,:),INTENT(INOUT)     :: PCGA_WVL    !O [-] assymetry factor for dust layer for all SW wavelengths
    REAL, DIMENSION(:,:,:,:),INTENT(INOUT)     :: PTAUREL_WVL !O [-] opt.depth/opt.depth(550) for dust layer for all SW wvl 
    REAL, DIMENSION(:,:,:), INTENT(INOUT)      :: PTAU550     !O [-] opt.depth at 550nm for all dust layer

    !LOCAL VARIABLES
    REAL, DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3), SIZE(PSVT,4)) :: ZSVT 
    REAL, DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3), NMODE_DST) :: ZMASS         ![kg/m3] mass of one dust mode
    REAL, DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3), NMODE_DST) :: ZRADIUS       ![um] number median radius of one dust mode
    REAL, DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3), NMODE_DST) :: ZSIGMA        ![-] dispersion coefficient one dust mode
    REAL, ALLOCATABLE, DIMENSION(:,:,:,:)                   :: ZTAU550_MDE   ![-] opt.depth 550nm one mode
    REAL, ALLOCATABLE, DIMENSION(:,:,:,:,:)                 :: ZTAU_WVL_MDE  ![-] opt.depth @ wvl, one mode
    REAL, ALLOCATABLE, DIMENSION(:,:,:,:,:)                 :: ZPIZA_WVL_MDE ![-] single scattering albedo @ wvl, one mode
    REAL, ALLOCATABLE, DIMENSION(:,:,:,:,:)                 :: ZCGA_WVL_MDE  ![-] assymetry factor @ wvl, one mode
    INTEGER                                                 :: JMDE        ![idx] counter for modes
    INTEGER                                                 :: JWVL        ![idx] counter for wavelengths
    !Allocate arrays which size depend on number of modes
    REAL(KIND=JPRB) :: ZHOOK_HANDLE
    IF (LHOOK) CALL DR_HOOK('MODE_DUSTOPT:DUSTOPT_GET',0,ZHOOK_HANDLE)
    ALLOCATE(ZTAU550_MDE(SIZE(PTAU550,1),SIZE(PTAU550,2),SIZE(PTAU550,3),NMODE_DST)) 
    ALLOCATE(ZTAU_WVL_MDE(SIZE(PTAU550,1),SIZE(PTAU550,2),SIZE(PTAU550,3),KSWB,NMODE_DST))
    ALLOCATE(ZPIZA_WVL_MDE(SIZE(PTAU550,1),SIZE(PTAU550,2),SIZE(PTAU550,3),KSWB,NMODE_DST))
    ALLOCATE(ZCGA_WVL_MDE(SIZE(PTAU550,1),SIZE(PTAU550,2),SIZE(PTAU550,3),KSWB,NMODE_DST))
    
       
    ZSVT(:,:,:,:) = PSVT(:,:,:,:)

    CALL PPP2DUST(     &
         ZSVT                                   & !I [moments/molec_{air}] moments of dust for all modes
         ,PRHODREF                              & !I [kg/m3] air density
         ,PSIG3D=ZSIGMA                         & !O [-] dispersion coefficient
         ,PRG3D=ZRADIUS                         & !O [um] number median radius
         ,PMASS3D=ZMASS                         & !O [kg/m3] mass of dust
         )
       
    DO JMDE=1,NMODE_DST
       !Get dust optical properties from look up tables
       CALL DUSTOPT_LKT(                     &
            ZRADIUS(:,:,:,JMDE)                   &  !I [um] number median radius for current mode
            ,ZSIGMA(:,:,:,JMDE)                   &  !I [none] dispersion coefficient for current mode
            ,ZMASS(:,:,:, JMDE)                    &  !I [kg/m3] Mass of dust for current mode
            ,ZTAU550_MDE(:,:,:,JMDE)         &  !O [-] optical depth at 550 nm wavelength
            ,ZTAU_WVL_MDE(:,:,:,:,JMDE)      &  !O [-] opt.depth(lambda)/opt.depth(550nm)
            ,ZPIZA_WVL_MDE(:,:,:,:,JMDE)     &  !O [-] single scattering coefficient at any wavelength
            ,ZCGA_WVL_MDE(:,:,:,:,JMDE)      &  !O [-] assymetry factor at any wavelength
            ,PZZ(:,:,:)                      &  !I [m] height of layers
            ,KSWB                            &  !I [nbr] number of shortwave bands
            )
    ENDDO  !Loop on modes

    !Erase earlier value of optical depth at 550 nm
    PTAU550(:,:,:)=0.d0
    
    !Get total at 550 nm from all modes 
    DO JMDE=1,NMODE_DST
       PTAU550(:,:,:) =        &        !Dust optical depth at 550 nm for all dust
            PTAU550(:,:,:)     &        !Dust optical depth at 550 nm for all dust
            + ZTAU550_MDE(:,:,:,JMDE)   !Optical depth for one mode at 550 nm      
    ENDDO
    
    !Initialize output variables
    PTAUREL_WVL(:,:,:,:)=0.d0         !Initialize opt.depth at wvl=lambda
    PCGA_WVL(:,:,:,:)=0.d0           !Initialize assym.factor at wvl=lambda
    PPIZA_WVL(:,:,:,:)=0.d0          !Initialize single scattering albedo at wvl=lambda


    !Find the numerator in the expression for the average of the optical properties   
    DO JMDE=1,NMODE_DST             !Number of modes
       DO JWVL=1,KSWB                    !Number of SW wavelengths

          !Get sum of optical depth from all modes at wvl
          PTAUREL_WVL(:,:,:,JWVL)  =                &  !new opt.depth(lambda) / opt.depth(550)
               PTAUREL_WVL(:,:,:,JWVL)              &  !old sum for all modes at wvl=lambda
               +ZTAU_WVL_MDE(:,:,:,JWVL,JMDE)       !optical depth for one mode at wvl=lambda
          
          !Get sum of all assymmetry factors  from all modes at wvl=lambda 
          PCGA_WVL(:,:,:,JWVL) =                     &  !New sum of assymetry factors
               PCGA_WVL(:,:,:,JWVL)                  &  !old sum of assymetry factors
               +ZCGA_WVL_MDE(:,:,:,JWVL,JMDE)     &  !Assymetry factor for one mode and one wavelength
               *ZTAU_WVL_MDE(:,:,:,JWVL,JMDE)     &  !Optical depth of this wavelength and mode
               *ZPIZA_WVL_MDE(:,:,:,JWVL,JMDE)       !Fraction of radiation scattered

          !Get sum of single scattering albdedo at wvl=lambda
          PPIZA_WVL(:,:,:,JWVL)  =                  & !New sum of single scattering albedo
               PPIZA_WVL(:,:,:,JWVL)             & !Old sum of single scattering albedo
               +ZPIZA_WVL_MDE(:,:,:,JWVL,JMDE)   & !SSA for onen mode and one wavelength
               *ZTAU_WVL_MDE(:,:,:,JWVL,JMDE)        !Optical depth for this wavelength and mode
          
       ENDDO
    ENDDO
  
    !Compute the output values for dust optical properties
    DO JWVL=1,KSWB
       
       !Divide total single scattering albdeo by total optical depth at this wavelength
       !This is needed since we weight all single scattering alebdos by wavelengths just above
       PPIZA_WVL(:,:,:,JWVL) =               &     !The value we want is ....
            PPIZA_WVL(:,:,:,JWVL)            &     !..value weighted by optical depths of all wvl and modes
            /max(epsilon,PTAUREL_WVL(:,:,:,JWVL))            !..divided by the optical depth for all wvl 
       
       
       !Divide total assymetry factor by total optical depth at this wavelength
       !This is needed since we weight all assymetry factors by wavelengths just above
       PCGA_WVL(:,:,:,JWVL) =              &  !The value we want is ....
            PCGA_WVL(:,:,:,JWVL)           &  !..value weighted by optical depths of all wvl and modes
            /                              &
            (max(epsilon,                  &
            (PTAUREL_WVL(:,:,:,JWVL)       &  !..divided scattered fraction of by the optical depth
            *PPIZA_WVL(:,:,:,JWVL))))

       !Finally convert PTAUREL_WVL which was until now an optical depth to a fraction of optical depth
       PTAUREL_WVL(:,:,:,JWVL) =                    &
            PTAUREL_WVL(:,:,:,JWVL)                 &  !Opt.depth at lambda with contr. from all modes
            /max(epsilon,PTAU550(:,:,:))               !Optical depth at 550 contr. from all modes
       
    ENDDO !Loop on wavelenghts
       

    !DEALLOCATE local  arrays which size depend on number of modes
    DEALLOCATE(ZTAU550_MDE)
    DEALLOCATE(ZTAU_WVL_MDE)
    DEALLOCATE(ZPIZA_WVL_MDE)
    DEALLOCATE(ZCGA_WVL_MDE)
    
  IF (LHOOK) CALL DR_HOOK('MODE_DUSTOPT:DUSTOPT_GET',1,ZHOOK_HANDLE)
  END SUBROUTINE DUSTOPT_GET

  !*****************************************************************

  SUBROUTINE DUSTOPT_LKT(              &
        PRG                            & !I [um] number median radius of aerosol mode
       ,PSIGMA                         & !I [-] lognormal dispersion coefficient
       ,PMASS                          & !I [kg/m3] Mass concentration of dust
       ,PTAU550                        & !O [optical depth at 550 nm
       ,PTAU_WVL                       & !O [-] opt.depth(lambda)/opt.depth(550nm)
       ,PPIZA_WVL                      & !O [-] single scattering coefficient at any wavelength
       ,PCGA_WVL                       & !O [-] assymetry factor at any wavelength
       ,PZZ                            & !I [m] height of layers
       ,KSWB                           & !I [nbr] number of short wave bands
       )

    !Purpose: Get optical properties of one dust mode from the mass concentration, 
    !dispersion coefficient and number median radius.

    !Use the module with the dust optical properties look up tables
    USE MODD_DUST_OPT_LKT  

    IMPLICIT NONE
    !INPUT
    REAL, DIMENSION(:,:,:), INTENT(IN)        :: PRG          !I [um] number median radius for one mode
    REAL, DIMENSION(:,:,:), INTENT(IN)        :: PSIGMA       !I [-] dispersion coefficient for one mode
    REAL, DIMENSION(:,:,:), INTENT(IN)        :: PMASS        !I [kg/m3] mass of dust
    REAL, DIMENSION(:,:,:), INTENT(IN)        :: PZZ          !I [m] height of layers
    INTEGER, INTENT(IN)                       :: KSWB         !I [nbr] number of shortwave bands

    !OUTPUT
    REAL, DIMENSION(:,:,:), INTENT(OUT)       :: PTAU550      !O [-] optical depth at 550 nm
    REAL, DIMENSION(:,:,:,:), INTENT(OUT)     :: PTAU_WVL     !O [-] optical depth at wvl
    REAL, DIMENSION(:,:,:,:), INTENT(OUT)     :: PPIZA_WVL    !O [-] single scattering albedo @ wvl
    REAL, DIMENSION(:,:,:,:), INTENT(OUT)     :: PCGA_WVL     !O [-] assymetry factor @ wvl

    !LOCALS
    REAL, DIMENSION(SIZE(PTAU550,1),SIZE(PTAU550,2),SIZE(PTAU550,3),KSWB) :: ZEXT_COEFF_WVL    ![m2/kg] Extinction coefficient at wvl
    REAL, DIMENSION(SIZE(PTAU550,1),SIZE(PTAU550,2),SIZE(PTAU550,3)     ) :: ZEXT_COEFF_550    ![m2/kg] Extinction coefficient at 550nm
    REAL                      :: FACT_SIGMA  ![-] factor needed to get right index in look up table for sigma
    REAL                      :: FACT_RADIUS ![-] factor needed to get right index in look up table for radius
    INTEGER                   :: WVL_IDX     ![idx] counter for wavelengths
    INTEGER                   :: JI, JJ, JK  ![idx] counters for lon, lat and lev
    INTEGER                   :: RG_IDX      ![idx] index for radius to get in look up table
    INTEGER                   :: SG_IDX      ![idx] index for sigma to get in look up table
    REAL,DIMENSION(SIZE(PRG,1),SIZE(PRG,2),SIZE(PRG,3)) :: ZRG     ![um] bounded value for number median radius
    REAL,DIMENSION(SIZE(PRG,1),SIZE(PRG,2),SIZE(PRG,3)) :: ZSIGMA  ![um] bounded value for sigma
    REAL, PARAMETER           :: EPSILON=1.d-8                     ![um] a small number used to avoid zero
    REAL                      :: ZRADIUS_LKT_MAX, ZRADIUS_LKT_MIN  ![um] values limited at midpoint values of bin
    REAL                      :: ZSIGMA_LKT_MAX, ZSIGMA_LKT_MIN    ![-] values limited at midpoint of bin
    INTEGER                   :: JKRAD                             !Index valid for radiation code

    !Limit max and min values to be midpont of bin to avoid 0 or NMAX+1 values
    REAL(KIND=JPRB) :: ZHOOK_HANDLE
    IF (LHOOK) CALL DR_HOOK('MODE_DUSTOPT:DUSTOPT_LKT',0,ZHOOK_HANDLE)
    ZRADIUS_LKT_MAX=exp(log(XRADIUS_LKT_MAX)   &
         - 0.5d0/DBLE(NMAX_RADIUS_LKT)*log(XRADIUS_LKT_MAX/XRADIUS_LKT_MIN))
    ZRADIUS_LKT_MIN=exp(log(XRADIUS_LKT_MIN)   &
         + 0.5d0/DBLE(NMAX_RADIUS_LKT)*log(XRADIUS_LKT_MAX/XRADIUS_LKT_MIN))
    ZSIGMA_LKT_MAX=XSIGMA_LKT_MAX - 0.5d0/DBLE(NMAX_SIGMA_LKT)*(XSIGMA_LKT_MAX-XSIGMA_LKT_MIN)
    ZSIGMA_LKT_MIN=XSIGMA_LKT_MIN + 0.5d0/DBLE(NMAX_SIGMA_LKT)*(XSIGMA_LKT_MAX-XSIGMA_LKT_MIN)

    !Begin code
    FACT_SIGMA = DBLE(NMAX_SIGMA_LKT)/(XSIGMA_LKT_MAX-XSIGMA_LKT_MIN)
    FACT_RADIUS = DBLE(NMAX_RADIUS_LKT)/(LOG(XRADIUS_LKT_MAX/XRADIUS_LKT_MIN))

    !Remove unphysical values for rg
    ZRG(:,:,:) = min( max(ZRADIUS_LKT_MIN,PRG(:,:,:)), ZRADIUS_LKT_MAX)
    ZSIGMA(:,:,:) = min( max(ZSIGMA_LKT_MIN,PSIGMA(:,:,:)), ZSIGMA_LKT_MAX)

    !Initilalize arrays to make sure, they are intent(OUT), 
    !and may be initialized strangely by the computer
    PTAU550(:,:,:)=EPSILON 
    PTAU_WVL(:,:,:,:)=EPSILON     
    PPIZA_WVL(:,:,:,:)=EPSILON
    PCGA_WVL(:,:,:,:)=EPSILON

    DO WVL_IDX = 1,KSWB
       DO JK=2,SIZE(PMASS,3)
          JKRAD = JK - 1  !Index in radiation code
          DO JJ=1,SIZE(PMASS,2)
             DO JI=1,SIZE(PMASS,1)

                !Get the correct indexes for the look up tables
                RG_IDX =  nint(                                  &
                     log(ZRG(JI,JJ,JK)/XRADIUS_LKT_MIN)          &
                     *FACT_RADIUS                                &
                     +0.5)
                   
                SG_IDX = nint((ZSIGMA(JI,JJ,JK)-XSIGMA_LKT_MIN)*FACT_SIGMA + 0.5) 

                !Open the look up tables and get the right values out of them
                !The extinction coefficient for 
                ZEXT_COEFF_WVL(JI,JJ,JK,WVL_IDX)         = XEXT_COEFF_WVL_LKT(RG_IDX,SG_IDX,WVL_IDX)
                ZEXT_COEFF_550(JI,JJ,JK)                 = XEXT_COEFF_550_LKT(RG_IDX,SG_IDX)
                !Switch to radiation code indexes for the output values
                PPIZA_WVL(JI,JJ,JKRAD,WVL_IDX)              = XPIZA_LKT(RG_IDX,SG_IDX,WVL_IDX)
                PCGA_WVL(JI,JJ,JKRAD,WVL_IDX)               = XCGA_LKT(RG_IDX,SG_IDX,WVL_IDX)
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    !Get the optical depth of this mode using the looked up extinction coeffient
    DO JK=2,SIZE(PZZ,3)-1
       JKRAD = JK - 1           !Index in radiation code
       PTAU550(:,:,JKRAD) = ZEXT_COEFF_550(:,:,JK) &
            * PMASS(:,:,JK)                    &
            * (PZZ(:,:,JK+1) - PZZ(:,:,JK))
    ENDDO

    !Get the optical depth of whatever wavelength using the looked up tables
    DO WVL_IDX=1,KSWB
       DO JK=2,SIZE(PZZ,3)-1
          JKRAD = JK -1
          PTAU_WVL(:,:,JKRAD,WVL_IDX) =              &
               PMASS(:,:,JK)                      & ![kg/m3] Mass in this mode
               *ZEXT_COEFF_WVL(:,:,JK,WVL_IDX)    & ![m2/kg] mass exinction coefficient
               *(PZZ(:,:,JK+1) - PZZ(:,:,JK))       ![m] Height of layer
       ENDDO !Loop on levels
    ENDDO    !Loop on wavelengths

    !Avoid unphysical values (which might occur on grid edges) on grid edges
    PTAU550(:,:,:)=max(PTAU550(:,:,:),EPSILON)
    PTAU_WVL(:,:,:,:)=max(PTAU_WVL(:,:,:,:),EPSILON)
    PPIZA_WVL(:,:,:,:)=max(PPIZA_WVL(:,:,:,:),EPSILON)
    PCGA_WVL(:,:,:,:)=max(PCGA_WVL(:,:,:,:),EPSILON)

  IF (LHOOK) CALL DR_HOOK('MODE_DUSTOPT:DUSTOPT_LKT',1,ZHOOK_HANDLE)
  END SUBROUTINE DUSTOPT_LKT

  !***************************************************************************
  SUBROUTINE DUST_OPT_LKT_SET1()

!Purpose: Read the look up tables for dust optical properties
!All variables are defined in the module MODD_DUST_OPT_LKT
! Based upon refractive indexes at wavelength intervals as:
!
! New tables from Mallet (LA) and Tulet (CNRM) :
! 0.185-0.25  RI[1]="(1.448,-0.00292)"
! 0.25-0.44   RI[2]="(1.448,-0.00292)"
! 0.44-0.69   RI[3]="(1.448,-0.00292)"
! 0.69-1.19   RI[4]="(1.44023,-0.00116)"
! 1.19-2.38   RI[5]="(1.41163,-0.00106)"
! 2.38-4.0    RI[6]="(1.41163,-0.00106)"
!             RI550="(1.44412,-0.00204)"

    USE MODD_DUST_OPT_LKT

    IMPLICIT NONE
    
    !Here are the output values from the mie program:
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_DUSTOPT:DUST_OPT_LKT_SET1',0,ZHOOK_HANDLE)
XEXT_COEFF_WVL_LKT(1,1,1:6)=(/ 92.520000,37.760000,21.553000,5.277700,2.711100,1.337700 /)
XPIZA_LKT(1,1,1:6)=(/ 0.431792,0.157096,0.049108,0.030213,0.005366,0.000659 /)
XCGA_LKT(1,1,1:6)=(/ 0.016660,0.006443,0.002727,0.001070,0.000340,0.000083 /)
XEXT_COEFF_550_LKT(1,1)=15.743000 !rg=0.0104084 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(1,2,1:6)=(/ 103.240000,39.382000,21.850000,5.321300,2.715200,1.338000 /)
XPIZA_LKT(1,2,1:6)=(/ 0.487818,0.189906,0.061069,0.037775,0.006754,0.000830 /)
XCGA_LKT(1,2,1:6)=(/ 0.021147,0.008187,0.003463,0.001360,0.000430,0.000107 /)
XEXT_COEFF_550_LKT(1,2)=16.050000 !rg=0.0104084 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(1,3,1:6)=(/ 129.230000,43.318000,22.568000,5.427000,2.725100,1.338600 /)
XPIZA_LKT(1,3,1:6)=(/ 0.585845,0.259826,0.088974,0.055729,0.010124,0.001248 /)
XCGA_LKT(1,3,1:6)=(/ 0.032163,0.012480,0.005287,0.002080,0.000657,0.000163 /)
XEXT_COEFF_550_LKT(1,3)=16.796000 !rg=0.0104084 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(1,4,1:6)=(/ 183.290000,51.589000,24.073000,5.649500,2.745900,1.340000 /)
XPIZA_LKT(1,4,1:6)=(/ 0.702104,0.373129,0.142704,0.091571,0.017184,0.002133 /)
XCGA_LKT(1,4,1:6)=(/ 0.054600,0.021287,0.009037,0.003557,0.001127,0.000277 /)
XEXT_COEFF_550_LKT(1,4)=18.363000 !rg=0.0104084 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(1,5,1:6)=(/ 286.070000,67.814000,27.025000,6.087600,2.786700,1.342500 /)
XPIZA_LKT(1,5,1:6)=(/ 0.803408,0.516645,0.231770,0.154935,0.030865,0.003880 /)
XCGA_LKT(1,5,1:6)=(/ 0.091107,0.035690,0.015190,0.005987,0.001897,0.000467 /)
XEXT_COEFF_550_LKT(1,5)=21.447000 !rg=0.0104084 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(1,6,1:6)=(/ 460.450000,97.279000,32.421000,6.891000,2.861200,1.347200 /)
XPIZA_LKT(1,6,1:6)=(/ 0.872997,0.656423,0.353869,0.250768,0.055064,0.007082 /)
XCGA_LKT(1,6,1:6)=(/ 0.142283,0.055803,0.023800,0.009397,0.002977,0.000733 /)
XEXT_COEFF_550_LKT(1,6)=27.095000 !rg=0.0104084 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(1,7,1:6)=(/ 733.310000,149.150000,42.096000,8.338200,2.995000,1.355600 /)
XPIZA_LKT(1,7,1:6)=(/ 0.916224,0.769783,0.495850,0.377503,0.095802,0.012818 /)
XCGA_LKT(1,7,1:6)=(/ 0.213453,0.083587,0.035690,0.014113,0.004477,0.001103 /)
XEXT_COEFF_550_LKT(1,7)=37.236000 !rg=0.0104084 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(1,8,1:6)=(/ 1130.600000,236.290000,59.060000,10.896000,3.231000,1.370300 /)
XPIZA_LKT(1,8,1:6)=(/ 0.942169,0.849454,0.634032,0.519983,0.159848,0.022828 /)
XCGA_LKT(1,8,1:6)=(/ 0.308897,0.122627,0.052300,0.020717,0.006583,0.001623 /)
XEXT_COEFF_550_LKT(1,8)=55.031000 !rg=0.0104084 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(1,9,1:6)=(/ 1705.100000,376.030000,88.742000,15.448000,3.650500,1.396300 /)
XPIZA_LKT(1,9,1:6)=(/ 0.958589,0.901098,0.750204,0.657642,0.253797,0.040195 /)
XCGA_LKT(1,9,1:6)=(/ 0.411487,0.179057,0.076007,0.030137,0.009593,0.002367 /)
XEXT_COEFF_550_LKT(1,9)=86.152000 !rg=0.0104084 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(1,10,1:6)=(/ 2437.700000,582.320000,139.140000,23.463000,4.390600,1.442100 /)
XPIZA_LKT(1,10,1:6)=(/ 0.968736,0.932546,0.835278,0.771030,0.376436,0.069456 /)
XCGA_LKT(1,10,1:6)=(/ 0.493533,0.260817,0.110187,0.043647,0.013920,0.003440 /)
XEXT_COEFF_550_LKT(1,10)=138.880000 !rg=0.0104084 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(1,11,1:6)=(/ 3314.000000,881.530000,221.820000,37.643000,5.710600,1.523500 /)
XPIZA_LKT(1,11,1:6)=(/ 0.975085,0.952114,0.892184,0.854112,0.517023,0.117531 /)
XCGA_LKT(1,11,1:6)=(/ 0.567660,0.368473,0.160363,0.063207,0.020190,0.004997 /)
XEXT_COEFF_550_LKT(1,11)=224.920000 !rg=0.0104084 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(1,12,1:6)=(/ 4242.300000,1306.600000,347.100000,62.286000,8.055400,1.668200 /)
XPIZA_LKT(1,12,1:6)=(/ 0.979085,0.964921,0.927413,0.909165,0.653952,0.191816 /)
XCGA_LKT(1,12,1:6)=(/ 0.624753,0.465660,0.234643,0.091640,0.029290,0.007260 /)
XEXT_COEFF_550_LKT(1,12)=354.350000 !rg=0.0104084 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(1,13,1:6)=(/ 5110.500000,1812.200000,525.740000,103.700000,12.207000,1.925100 /)
XPIZA_LKT(1,13,1:6)=(/ 0.981483,0.972719,0.948753,0.943296,0.768173,0.296860 /)
XCGA_LKT(1,13,1:6)=(/ 0.668683,0.538647,0.339447,0.133300,0.042487,0.010553 /)
XEXT_COEFF_550_LKT(1,13)=538.780000 !rg=0.0104084 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(1,14,1:6)=(/ 5748.900000,2387.700000,784.190000,168.610000,19.470000,2.379400 /)
XPIZA_LKT(1,14,1:6)=(/ 0.982603,0.977645,0.962617,0.963443,0.851580,0.427755 /)
XCGA_LKT(1,14,1:6)=(/ 0.698187,0.605163,0.448927,0.194953,0.061620,0.015337 /)
XEXT_COEFF_550_LKT(1,14)=805.170000 !rg=0.0104084 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(1,15,1:6)=(/ 6015.700000,2950.500000,1106.000000,262.040000,31.999000,3.182800 /)
XPIZA_LKT(1,15,1:6)=(/ 0.982578,0.980638,0.971414,0.975083,0.907112,0.568583 /)
XCGA_LKT(1,15,1:6)=(/ 0.714507,0.653480,0.524447,0.285813,0.089433,0.022283 /)
XEXT_COEFF_550_LKT(1,15)=1131.400000 !rg=0.0104084 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(1,16,1:6)=(/ 5851.500000,3422.500000,1472.900000,395.050000,52.980000,4.605800 /)
XPIZA_LKT(1,16,1:6)=(/ 0.981412,0.982285,0.976787,0.982146,0.941826,0.698244 /)
XCGA_LKT(1,16,1:6)=(/ 0.718937,0.689083,0.593623,0.402740,0.130147,0.032367 /)
XEXT_COEFF_550_LKT(1,16)=1505.100000 !rg=0.0104084 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(1,17,1:6)=(/ 5366.800000,3697.700000,1840.800000,580.530000,85.671000,7.093600 /)
XPIZA_LKT(1,17,1:6)=(/ 0.979145,0.982774,0.980136,0.986797,0.962418,0.800732 /)
XCGA_LKT(1,17,1:6)=(/ 0.715933,0.710930,0.644890,0.497120,0.190150,0.046947 /)
XEXT_COEFF_550_LKT(1,17)=1875.900000 !rg=0.0104084 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(1,18,1:6)=(/ 4801.300000,3708.500000,2163.800000,789.620000,132.270000,11.407000 /)
XPIZA_LKT(1,18,1:6)=(/ 0.976068,0.982123,0.982064,0.989568,0.974360,0.873176 /)
XCGA_LKT(1,18,1:6)=(/ 0.713493,0.719603,0.683553,0.565050,0.278340,0.068027 /)
XEXT_COEFF_550_LKT(1,18)=2198.200000 !rg=0.0104084 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(1,19,1:6)=(/ 4249.400000,3478.400000,2378.500000,1023.100000,196.690000,18.760000 /)
XPIZA_LKT(1,19,1:6)=(/ 0.972714,0.980276,0.982817,0.991341,0.981545,0.920482 /)
XCGA_LKT(1,19,1:6)=(/ 0.718660,0.717890,0.708233,0.627220,0.394907,0.098600 /)
XEXT_COEFF_550_LKT(1,19)=2408.700000 !rg=0.0104084 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(1,20,1:6)=(/ 3667.400000,3122.100000,2430.400000,1242.500000,287.790000,30.753000 /)
XPIZA_LKT(1,20,1:6)=(/ 0.968698,0.977438,0.982433,0.992399,0.986356,0.949592 /)
XCGA_LKT(1,20,1:6)=(/ 0.723267,0.712220,0.719883,0.671263,0.497200,0.143163 /)
XEXT_COEFF_550_LKT(1,20)=2452.500000 !rg=0.0104084 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(2,1,1:6)=(/ 103.740000,39.460000,21.865000,5.323400,2.715400,1.338000 /)
XPIZA_LKT(2,1,1:6)=(/ 0.489711,0.191161,0.061553,0.038086,0.006812,0.000837 /)
XCGA_LKT(2,1,1:6)=(/ 0.019537,0.007560,0.003200,0.001257,0.000397,0.000097 /)
XEXT_COEFF_550_LKT(2,1)=16.065000 !rg=0.011276 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(2,2,1:6)=(/ 117.330000,41.516000,22.240000,5.378600,2.720600,1.338300 /)
XPIZA_LKT(2,2,1:6)=(/ 0.545773,0.229082,0.076290,0.047519,0.008571,0.001055 /)
XCGA_LKT(2,2,1:6)=(/ 0.024793,0.009603,0.004067,0.001597,0.000507,0.000123 /)
XEXT_COEFF_550_LKT(2,2)=16.454000 !rg=0.011276 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(2,3,1:6)=(/ 150.230000,46.504000,23.149000,5.512700,2.733200,1.339100 /)
XPIZA_LKT(2,3,1:6)=(/ 0.640299,0.307762,0.110299,0.069760,0.012835,0.001586 /)
XCGA_LKT(2,3,1:6)=(/ 0.037687,0.014633,0.006203,0.002440,0.000770,0.000190 /)
XEXT_COEFF_550_LKT(2,3)=17.399000 !rg=0.011276 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(2,4,1:6)=(/ 218.290000,56.986000,25.055000,5.794800,2.759500,1.340800 /)
XPIZA_LKT(2,4,1:6)=(/ 0.746746,0.429423,0.174353,0.113519,0.021742,0.002710 /)
XCGA_LKT(2,4,1:6)=(/ 0.063893,0.024950,0.010600,0.004173,0.001320,0.000327 /)
XEXT_COEFF_550_LKT(2,4)=19.387000 !rg=0.011276 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(2,5,1:6)=(/ 345.760000,77.517000,28.796000,6.350700,2.811200,1.344100 /)
XPIZA_LKT(2,5,1:6)=(/ 0.834754,0.573976,0.276580,0.188852,0.038904,0.004927 /)
XCGA_LKT(2,5,1:6)=(/ 0.106493,0.041807,0.017810,0.007023,0.002223,0.000547 /)
XEXT_COEFF_550_LKT(2,5)=23.299000 !rg=0.011276 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(2,6,1:6)=(/ 556.620000,114.640000,35.633000,7.370400,2.905600,1.350000 /)
XPIZA_LKT(2,6,1:6)=(/ 0.892861,0.705457,0.409267,0.298127,0.068943,0.008986 /)
XCGA_LKT(2,6,1:6)=(/ 0.166257,0.065340,0.027897,0.011020,0.003493,0.000860 /)
XEXT_COEFF_550_LKT(2,6)=30.459000 !rg=0.011276 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(2,7,1:6)=(/ 877.950000,179.370000,47.880000,9.207500,3.075300,1.360600 /)
XPIZA_LKT(2,7,1:6)=(/ 0.928292,0.805941,0.553705,0.434671,0.118628,0.016238 /)
XCGA_LKT(2,7,1:6)=(/ 0.248543,0.097890,0.041817,0.016550,0.005253,0.001297 /)
XEXT_COEFF_550_LKT(2,7)=43.302000 !rg=0.011276 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(2,8,1:6)=(/ 1343.100000,286.000000,69.289000,12.455000,3.374700,1.379200 /)
XPIZA_LKT(2,8,1:6)=(/ 0.949787,0.873449,0.685079,0.578317,0.194580,0.028837 /)
XCGA_LKT(2,8,1:6)=(/ 0.350880,0.143747,0.061267,0.024287,0.007723,0.001903 /)
XEXT_COEFF_550_LKT(2,8)=65.758000 !rg=0.011276 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(2,9,1:6)=(/ 1991.200000,451.820000,106.470000,18.228000,3.907000,1.412200 /)
XPIZA_LKT(2,9,1:6)=(/ 0.963390,0.915917,0.789123,0.708153,0.301475,0.050533 /)
XCGA_LKT(2,9,1:6)=(/ 0.445603,0.210107,0.089063,0.035320,0.011250,0.002777 /)
XEXT_COEFF_550_LKT(2,9)=104.720000 !rg=0.011276 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(2,10,1:6)=(/ 2778.000000,691.380000,168.610000,28.374000,4.846200,1.470200 /)
XPIZA_LKT(2,10,1:6)=(/ 0.971647,0.941603,0.861814,0.809101,0.433519,0.086622 /)
XCGA_LKT(2,10,1:6)=(/ 0.523827,0.303943,0.129253,0.051140,0.016323,0.004033 /)
XEXT_COEFF_550_LKT(2,10)=169.590000 !rg=0.011276 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(2,11,1:6)=(/ 3694.000000,1041.700000,267.610000,46.237000,6.520900,1.573500 /)
XPIZA_LKT(2,11,1:6)=(/ 0.976966,0.958062,0.908812,0.879885,0.575374,0.144689 /)
XCGA_LKT(2,11,1:6)=(/ 0.592367,0.412317,0.188487,0.074063,0.023673,0.005860 /)
XEXT_COEFF_550_LKT(2,11)=272.300000 !rg=0.011276 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(2,12,1:6)=(/ 4618.100000,1504.900000,412.780000,76.935000,9.493100,1.757000 /)
XPIZA_LKT(2,12,1:6)=(/ 0.980232,0.968603,0.937403,0.925366,0.704711,0.231541 /)
XCGA_LKT(2,12,1:6)=(/ 0.643963,0.495827,0.275467,0.107457,0.034327,0.008517 /)
XEXT_COEFF_550_LKT(2,12)=422.070000 !rg=0.011276 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(2,13,1:6)=(/ 5413.500000,2043.900000,621.270000,127.360000,14.744000,2.083100 /)
XPIZA_LKT(2,13,1:6)=(/ 0.982066,0.974988,0.955133,0.952971,0.806551,0.348772 /)
XCGA_LKT(2,13,1:6)=(/ 0.682020,0.567970,0.388163,0.156593,0.049783,0.012377 /)
XEXT_COEFF_550_LKT(2,13)=637.700000 !rg=0.011276 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(2,14,1:6)=(/ 5913.800000,2619.700000,913.080000,203.540000,23.886000,2.659500 /)
XPIZA_LKT(2,14,1:6)=(/ 0.982723,0.979049,0.966827,0.969037,0.877713,0.486431 /)
XCGA_LKT(2,14,1:6)=(/ 0.706313,0.625800,0.483003,0.229497,0.072207,0.017987 /)
XEXT_COEFF_550_LKT(2,14)=935.960000 !rg=0.011276 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(2,15,1:6)=(/ 6000.100000,3158.300000,1249.300000,310.900000,39.482000,3.678800 /)
XPIZA_LKT(2,15,1:6)=(/ 0.982241,0.981474,0.973882,0.978376,0.923657,0.625059 /)
XCGA_LKT(2,15,1:6)=(/ 0.717430,0.669540,0.553157,0.333777,0.104867,0.026123 /)
XEXT_COEFF_550_LKT(2,15)=1278.000000 !rg=0.011276 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(2,16,1:6)=(/ 5684.500000,3562.200000,1628.300000,466.940000,65.029000,5.481900 /)
XPIZA_LKT(2,16,1:6)=(/ 0.980599,0.982639,0.978404,0.984335,0.951778,0.744846 /)
XCGA_LKT(2,16,1:6)=(/ 0.718130,0.699263,0.616980,0.447937,0.152870,0.037933 /)
XEXT_COEFF_550_LKT(2,16)=1660.700000 !rg=0.011276 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(2,17,1:6)=(/ 5141.000000,3734.000000,1986.800000,665.130000,103.390000,8.625700 /)
XPIZA_LKT(2,17,1:6)=(/ 0.977913,0.982642,0.981098,0.988136,0.968216,0.834668 /)
XCGA_LKT(2,17,1:6)=(/ 0.713957,0.715560,0.662640,0.525477,0.223803,0.055013 /)
XEXT_COEFF_550_LKT(2,17)=2021.600000 !rg=0.011276 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(2,18,1:6)=(/ 4574.400000,3638.100000,2271.100000,886.600000,156.520000,14.043000 /)
XPIZA_LKT(2,18,1:6)=(/ 0.974904,0.981497,0.982498,0.990398,0.977764,0.895751 /)
XCGA_LKT(2,18,1:6)=(/ 0.715933,0.719760,0.695253,0.593790,0.325740,0.079730 /)
XEXT_COEFF_550_LKT(2,18)=2303.900000 !rg=0.011276 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(2,19,1:6)=(/ 4020.000000,3334.500000,2421.300000,1115.600000,231.770000,23.157000 /)
XPIZA_LKT(2,19,1:6)=(/ 0.971278,0.979229,0.982793,0.991842,0.983800,0.934593 /)
XCGA_LKT(2,19,1:6)=(/ 0.721507,0.715440,0.714623,0.646460,0.443660,0.115677 /)
XEXT_COEFF_550_LKT(2,19)=2447.500000 !rg=0.011276 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(2,20,1:6)=(/ 3434.300000,2974.300000,2397.500000,1320.700000,331.450000,37.640000 /)
XPIZA_LKT(2,20,1:6)=(/ 0.966326,0.976205,0.981925,0.992685,0.987818,0.958051 /)
XCGA_LKT(2,20,1:6)=(/ 0.724030,0.712603,0.720917,0.685953,0.527733,0.168320 /)
XEXT_COEFF_550_LKT(2,20)=2415.200000 !rg=0.011276 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(3,1,1:6)=(/ 117.960000,41.614000,22.259000,5.381200,2.720900,1.338400 /)
XPIZA_LKT(3,1,1:6)=(/ 0.547621,0.230509,0.076882,0.047907,0.008644,0.001064 /)
XCGA_LKT(3,1,1:6)=(/ 0.022907,0.008867,0.003753,0.001477,0.000467,0.000113 /)
XEXT_COEFF_550_LKT(3,1)=16.473000 !rg=0.0122159 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(3,2,1:6)=(/ 135.180000,44.219000,22.734000,5.451200,2.727400,1.338800 /)
XPIZA_LKT(3,2,1:6)=(/ 0.602197,0.273509,0.094903,0.059618,0.010870,0.001341 /)
XCGA_LKT(3,2,1:6)=(/ 0.029067,0.011263,0.004770,0.001877,0.000593,0.000147 /)
XEXT_COEFF_550_LKT(3,2)=16.966000 !rg=0.0122159 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(3,3,1:6)=(/ 176.800000,50.542000,23.885000,5.621200,2.743400,1.339800 /)
XPIZA_LKT(3,3,1:6)=(/ 0.690964,0.360080,0.135939,0.086990,0.016259,0.002016 /)
XCGA_LKT(3,3,1:6)=(/ 0.044153,0.017160,0.007277,0.002863,0.000907,0.000223 /)
XEXT_COEFF_550_LKT(3,3)=18.164000 !rg=0.0122159 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(3,4,1:6)=(/ 262.170000,63.821000,26.299000,5.979100,2.776700,1.341900 /)
XPIZA_LKT(3,4,1:6)=(/ 0.786163,0.487337,0.211219,0.139896,0.027473,0.003443 /)
XCGA_LKT(3,4,1:6)=(/ 0.074740,0.029237,0.012430,0.004897,0.001550,0.000380 /)
XEXT_COEFF_550_LKT(3,4)=20.684000 !rg=0.0122159 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(3,5,1:6)=(/ 419.010000,89.770000,31.040000,6.684600,2.842200,1.346100 /)
XPIZA_LKT(3,5,1:6)=(/ 0.861225,0.628955,0.326229,0.228145,0.048925,0.006256 /)
XCGA_LKT(3,5,1:6)=(/ 0.124397,0.048957,0.020880,0.008240,0.002610,0.000643 /)
XEXT_COEFF_550_LKT(3,5)=25.646000 !rg=0.0122159 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(3,6,1:6)=(/ 671.210000,136.400000,39.701000,7.978800,2.961900,1.353600 /)
XPIZA_LKT(3,6,1:6)=(/ 0.909216,0.749553,0.466818,0.350170,0.085992,0.011396 /)
XCGA_LKT(3,6,1:6)=(/ 0.193970,0.076497,0.032690,0.012927,0.004100,0.001010 /)
XEXT_COEFF_550_LKT(3,6)=34.722000 !rg=0.0122159 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(3,7,1:6)=(/ 1048.200000,216.630000,55.188000,10.311000,3.177100,1.367000 /)
XPIZA_LKT(3,7,1:6)=(/ 0.938278,0.836859,0.609746,0.493498,0.145993,0.020551 /)
XCGA_LKT(3,7,1:6)=(/ 0.287080,0.114640,0.048987,0.019403,0.006163,0.001520 /)
XEXT_COEFF_550_LKT(3,7)=50.966000 !rg=0.0122159 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(3,8,1:6)=(/ 1588.400000,345.550000,82.140000,14.432000,3.556900,1.390500 /)
XPIZA_LKT(3,8,1:6)=(/ 0.956153,0.893261,0.731466,0.634353,0.234705,0.036368 /)
XCGA_LKT(3,8,1:6)=(/ 0.390013,0.168527,0.071773,0.028467,0.009060,0.002237 /)
XEXT_COEFF_550_LKT(3,8)=79.229000 !rg=0.0122159 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(3,9,1:6)=(/ 2297.900000,539.990000,128.460000,21.750000,4.232500,1.432300 /)
XPIZA_LKT(3,9,1:6)=(/ 0.967263,0.927974,0.822709,0.753756,0.353778,0.063349 /)
XCGA_LKT(3,9,1:6)=(/ 0.476813,0.246147,0.104387,0.041387,0.013197,0.003260 /)
XEXT_COEFF_550_LKT(3,9)=127.710000 !rg=0.0122159 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(3,10,1:6)=(/ 3146.100000,819.440000,204.240000,34.572000,5.424300,1.505900 /)
XPIZA_LKT(3,10,1:6)=(/ 0.974108,0.949179,0.883857,0.841852,0.492284,0.107528 /)
XCGA_LKT(3,10,1:6)=(/ 0.553383,0.349627,0.151713,0.059917,0.019140,0.004733 /)
XEXT_COEFF_550_LKT(3,10)=206.630000 !rg=0.0122159 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(3,11,1:6)=(/ 4072.200000,1221.800000,321.000000,56.990000,7.548300,1.636900 /)
XPIZA_LKT(3,11,1:6)=(/ 0.978476,0.963032,0.922295,0.901309,0.631501,0.176837 /)
XCGA_LKT(3,11,1:6)=(/ 0.613887,0.449493,0.221580,0.086800,0.027750,0.006877 /)
XEXT_COEFF_550_LKT(3,11)=327.410000 !rg=0.0122159 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(3,12,1:6)=(/ 4964.200000,1713.000000,488.770000,94.924000,11.313000,1.869700 /)
XPIZA_LKT(3,12,1:6)=(/ 0.981136,0.971550,0.945595,0.938518,0.750611,0.276619 /)
XCGA_LKT(3,12,1:6)=(/ 0.660627,0.524873,0.321297,0.126077,0.040227,0.009990 /)
XEXT_COEFF_550_LKT(3,12)=500.580000 !rg=0.0122159 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(3,13,1:6)=(/ 5660.800000,2283.500000,732.950000,155.440000,17.944000,2.283600 /)
XPIZA_LKT(3,13,1:6)=(/ 0.982483,0.976910,0.960594,0.960696,0.839616,0.404441 /)
XCGA_LKT(3,13,1:6)=(/ 0.692997,0.594180,0.432363,0.184137,0.058333,0.014517 /)
XEXT_COEFF_550_LKT(3,13)=752.670000 !rg=0.0122159 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(3,14,1:6)=(/ 6003.900000,2852.300000,1046.700000,243.600000,29.407000,3.015000 /)
XPIZA_LKT(3,14,1:6)=(/ 0.982662,0.980191,0.970186,0.973487,0.899469,0.545325 /)
XCGA_LKT(3,14,1:6)=(/ 0.711893,0.644943,0.511987,0.269877,0.084623,0.021087 /)
XEXT_COEFF_550_LKT(3,14)=1071.200000 !rg=0.0122159 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(3,15,1:6)=(/ 5916.800000,3345.300000,1403.900000,368.220000,48.673000,4.307500 /)
XPIZA_LKT(3,15,1:6)=(/ 0.981714,0.982054,0.975973,0.981123,0.937112,0.678113 /)
XCGA_LKT(3,15,1:6)=(/ 0.718627,0.683010,0.582010,0.384037,0.123033,0.030623 /)
XEXT_COEFF_550_LKT(3,15)=1435.400000 !rg=0.0122159 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(3,16,1:6)=(/ 5478.800000,3666.200000,1778.000000,546.840000,79.342000,6.590200 /)
XPIZA_LKT(3,16,1:6)=(/ 0.979616,0.982749,0.979670,0.986165,0.959737,0.786195 /)
XCGA_LKT(3,16,1:6)=(/ 0.716330,0.707423,0.636630,0.484270,0.179723,0.044453 /)
XEXT_COEFF_550_LKT(3,16)=1812.300000 !rg=0.0122159 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(3,17,1:6)=(/ 4907.300000,3725.800000,2114.800000,752.700000,123.630000,10.555000 /)
XPIZA_LKT(3,17,1:6)=(/ 0.976767,0.982311,0.981826,0.989196,0.972823,0.863522 /)
XCGA_LKT(3,17,1:6)=(/ 0.714543,0.718613,0.677650,0.553597,0.263273,0.064463 /)
XEXT_COEFF_550_LKT(3,17)=2149.000000 !rg=0.0122159 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(3,18,1:6)=(/ 4361.000000,3529.900000,2350.200000,984.390000,184.580000,17.328000 /)
XPIZA_LKT(3,18,1:6)=(/ 0.973539,0.980681,0.982758,0.991098,0.980575,0.914388 /)
XCGA_LKT(3,18,1:6)=(/ 0.718843,0.718403,0.704433,0.618487,0.376813,0.093473 /)
XEXT_COEFF_550_LKT(3,18)=2381.100000 !rg=0.0122159 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(3,19,1:6)=(/ 3779.800000,3186.300000,2431.500000,1207.700000,271.840000,28.521000 /)
XPIZA_LKT(3,19,1:6)=(/ 0.969284,0.978079,0.982560,0.992253,0.985726,0.946005 /)
XCGA_LKT(3,19,1:6)=(/ 0.722557,0.714117,0.718290,0.664513,0.484347,0.135813 /)
XEXT_COEFF_550_LKT(3,19)=2455.100000 !rg=0.0122159 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(3,20,1:6)=(/ 3203.500000,2842.300000,2337.800000,1390.300000,376.110000,45.733000 /)
XPIZA_LKT(3,20,1:6)=(/ 0.963890,0.974757,0.981213,0.992874,0.988973,0.964790 /)
XCGA_LKT(3,20,1:6)=(/ 0.725090,0.713283,0.720253,0.698123,0.555967,0.198110 /)
XEXT_COEFF_550_LKT(3,20)=2352.500000 !rg=0.0122159 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(4,1,1:6)=(/ 135.990000,44.342000,22.758000,5.454500,2.727800,1.338800 /)
XPIZA_LKT(4,1,1:6)=(/ 0.603958,0.275097,0.095621,0.060096,0.010963,0.001353 /)
XCGA_LKT(4,1,1:6)=(/ 0.026860,0.010403,0.004403,0.001730,0.000547,0.000133 /)
XEXT_COEFF_550_LKT(4,1)=16.990000 !rg=0.0132342 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(4,2,1:6)=(/ 157.810000,47.644000,23.359000,5.543300,2.736100,1.339300 /)
XPIZA_LKT(4,2,1:6)=(/ 0.655675,0.322804,0.117451,0.074547,0.013777,0.001704 /)
XCGA_LKT(4,2,1:6)=(/ 0.034070,0.013213,0.005597,0.002200,0.000697,0.000170 /)
XEXT_COEFF_550_LKT(4,2)=17.615000 !rg=0.0132342 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(4,3,1:6)=(/ 210.380000,55.660000,24.816000,5.758800,2.756200,1.340600 /)
XPIZA_LKT(4,3,1:6)=(/ 0.736988,0.415754,0.166379,0.107970,0.020577,0.002562 /)
XCGA_LKT(4,3,1:6)=(/ 0.051720,0.020120,0.008537,0.003360,0.001063,0.000260 /)
XEXT_COEFF_550_LKT(4,3)=19.134000 !rg=0.0132342 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(4,4,1:6)=(/ 316.960000,72.477000,27.875000,6.212900,2.798500,1.343300 /)
XPIZA_LKT(4,4,1:6)=(/ 0.820320,0.545288,0.253403,0.171196,0.034661,0.004373 /)
XCGA_LKT(4,4,1:6)=(/ 0.087393,0.034257,0.014577,0.005743,0.001820,0.000447 /)
XEXT_COEFF_550_LKT(4,4)=22.330000 !rg=0.0132342 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(4,5,1:6)=(/ 508.090000,105.220000,33.883000,7.108200,2.881400,1.348500 /)
XPIZA_LKT(4,5,1:6)=(/ 0.883319,0.680299,0.379933,0.272812,0.061359,0.007940 /)
XCGA_LKT(4,5,1:6)=(/ 0.145183,0.057317,0.024473,0.009663,0.003063,0.000753 /)
XEXT_COEFF_550_LKT(4,5)=28.622000 !rg=0.0132342 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(4,6,1:6)=(/ 806.970000,163.530000,44.849000,8.751300,3.033200,1.358000 /)
XPIZA_LKT(4,6,1:6)=(/ 0.922669,0.788356,0.524964,0.405950,0.106764,0.014442 /)
XCGA_LKT(4,6,1:6)=(/ 0.225563,0.089537,0.038300,0.015157,0.004810,0.001187 /)
XEXT_COEFF_550_LKT(4,6)=40.119000 !rg=0.0132342 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(4,7,1:6)=(/ 1248.400000,262.090000,64.405000,11.712000,3.306200,1.375000 /)
XPIZA_LKT(4,7,1:6)=(/ 0.946609,0.862880,0.662571,0.552350,0.178367,0.025978 /)
XCGA_LKT(4,7,1:6)=(/ 0.326890,0.134260,0.057380,0.022750,0.007230,0.001783 /)
XEXT_COEFF_550_LKT(4,7)=60.633000 !rg=0.0132342 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(4,8,1:6)=(/ 1861.600000,415.950000,98.214000,16.939000,3.788200,1.404900 /)
XPIZA_LKT(4,8,1:6)=(/ 0.961374,0.909478,0.772664,0.686759,0.280171,0.045770 /)
XCGA_LKT(4,8,1:6)=(/ 0.425003,0.197520,0.084080,0.033360,0.010627,0.002623 /)
XEXT_COEFF_550_LKT(4,8)=96.065000 !rg=0.0132342 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(4,9,1:6)=(/ 2629.300000,642.570000,155.490000,26.207000,4.645500,1.457800 /)
XPIZA_LKT(4,9,1:6)=(/ 0.970432,0.937854,0.851193,0.794040,0.409728,0.079137 /)
XCGA_LKT(4,9,1:6)=(/ 0.507727,0.286953,0.122390,0.048490,0.015473,0.003823 /)
XEXT_COEFF_550_LKT(4,9)=155.910000 !rg=0.0132342 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(4,10,1:6)=(/ 3525.000000,969.860000,246.670000,42.373000,6.157600,1.551100 /)
XPIZA_LKT(4,10,1:6)=(/ 0.976159,0.955603,0.901947,0.869585,0.551102,0.132729 /)
XCGA_LKT(4,10,1:6)=(/ 0.579427,0.393593,0.178183,0.070200,0.022437,0.005553 /)
XEXT_COEFF_550_LKT(4,10)=250.610000 !rg=0.0132342 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(4,11,1:6)=(/ 4454.200000,1414.100000,382.650000,70.358000,8.850500,1.717300 /)
XPIZA_LKT(4,11,1:6)=(/ 0.979729,0.967053,0.933215,0.918925,0.684060,0.214301 /)
XCGA_LKT(4,11,1:6)=(/ 0.634003,0.480940,0.260127,0.101753,0.032523,0.008067 /)
XEXT_COEFF_550_LKT(4,11)=390.980000 !rg=0.0132342 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(4,12,1:6)=(/ 5284.600000,1938.000000,578.030000,116.740000,13.614000,2.012800 /)
XPIZA_LKT(4,12,1:6)=(/ 0.981802,0.974000,0.952474,0.949113,0.791225,0.326673 /)
XCGA_LKT(4,12,1:6)=(/ 0.674957,0.554537,0.369380,0.148033,0.047137,0.011717 /)
XEXT_COEFF_550_LKT(4,12)=593.020000 !rg=0.0132342 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(4,13,1:6)=(/ 5856.300000,2517.300000,857.220000,188.210000,21.968000,2.538000 /)
XPIZA_LKT(4,13,1:6)=(/ 0.982674,0.978456,0.965157,0.966835,0.867655,0.462561 /)
XCGA_LKT(4,13,1:6)=(/ 0.702090,0.615980,0.468663,0.216657,0.068350,0.017023 /)
XEXT_COEFF_550_LKT(4,13)=879.180000 !rg=0.0132342 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(4,14,1:6)=(/ 6023.400000,3070.400000,1186.100000,289.610000,36.269000,3.465900 /)
XPIZA_LKT(4,14,1:6)=(/ 0.982419,0.981130,0.972867,0.977071,0.917388,0.602799 /)
XCGA_LKT(4,14,1:6)=(/ 0.715883,0.662093,0.540530,0.315857,0.099207,0.024723 /)
XEXT_COEFF_550_LKT(4,14)=1213.400000 !rg=0.0132342 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(4,15,1:6)=(/ 5774.600000,3503.000000,1560.600000,435.830000,59.825000,5.104000 /)
XPIZA_LKT(4,15,1:6)=(/ 0.980989,0.982491,0.977735,0.983468,0.947973,0.726708 /)
XCGA_LKT(4,15,1:6)=(/ 0.718240,0.694237,0.607020,0.431060,0.144450,0.035890 /)
XEXT_COEFF_550_LKT(4,15)=1592.600000 !rg=0.0132342 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(4,16,1:6)=(/ 5251.700000,3722.800000,1927.200000,630.140000,96.041000,7.990200 /)
XPIZA_LKT(4,16,1:6)=(/ 0.978604,0.982714,0.980708,0.987630,0.966071,0.822169 /)
XCGA_LKT(4,16,1:6)=(/ 0.715820,0.713153,0.655047,0.513870,0.211450,0.052093 /)
XEXT_COEFF_550_LKT(4,16)=1962.400000 !rg=0.0132342 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(4,17,1:6)=(/ 4693.900000,3673.500000,2229.200000,847.240000,146.650000,12.976000 /)
XPIZA_LKT(4,17,1:6)=(/ 0.975466,0.981763,0.982317,0.990080,0.976513,0.887717 /)
XCGA_LKT(4,17,1:6)=(/ 0.715887,0.719313,0.690183,0.582740,0.308653,0.075547 /)
XEXT_COEFF_550_LKT(4,17)=2263.200000 !rg=0.0132342 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(4,18,1:6)=(/ 4132.100000,3400.300000,2405.700000,1077.600000,217.610000,21.393000 /)
XPIZA_LKT(4,18,1:6)=(/ 0.971902,0.979688,0.982788,0.991649,0.982970,0.929627 /)
XCGA_LKT(4,18,1:6)=(/ 0.720730,0.716660,0.711720,0.638627,0.426830,0.109637 /)
XEXT_COEFF_550_LKT(4,18)=2433.500000 !rg=0.0132342 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(4,19,1:6)=(/ 3529.600000,3043.700000,2411.300000,1290.300000,314.720000,34.973000 /)
XPIZA_LKT(4,19,1:6)=(/ 0.967404,0.976734,0.982134,0.992580,0.987307,0.955171 /)
XCGA_LKT(4,19,1:6)=(/ 0.723760,0.712587,0.720253,0.680220,0.516730,0.159610 /)
XEXT_COEFF_550_LKT(4,19)=2431.000000 !rg=0.0132342 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(4,20,1:6)=(/ 3001.100000,2703.200000,2258.900000,1444.200000,423.470000,55.068000 /)
XPIZA_LKT(4,20,1:6)=(/ 0.961274,0.973500,0.980258,0.993006,0.989915,0.970138 /)
XCGA_LKT(4,20,1:6)=(/ 0.729503,0.716473,0.717963,0.707840,0.585287,0.233330 /)
XEXT_COEFF_550_LKT(4,20)=2269.400000 !rg=0.0132342 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(5,1,1:6)=(/ 158.850000,47.800000,23.389000,5.547400,2.736500,1.339400 /)
XPIZA_LKT(5,1,1:6)=(/ 0.657319,0.324535,0.118314,0.075136,0.013895,0.001719 /)
XCGA_LKT(5,1,1:6)=(/ 0.031487,0.012203,0.005167,0.002030,0.000643,0.000157 /)
XEXT_COEFF_550_LKT(5,1)=17.645000 !rg=0.0143374 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(5,2,1:6)=(/ 186.470000,51.985000,24.150000,5.659900,2.747000,1.340000 /)
XPIZA_LKT(5,2,1:6)=(/ 0.705103,0.376215,0.144465,0.092836,0.017448,0.002166 /)
XCGA_LKT(5,2,1:6)=(/ 0.039933,0.015497,0.006567,0.002583,0.000817,0.000200 /)
XEXT_COEFF_550_LKT(5,2)=18.437000 !rg=0.0143374 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(5,3,1:6)=(/ 252.710000,62.147000,25.995000,5.933300,2.772500,1.341700 /)
XPIZA_LKT(5,3,1:6)=(/ 0.777890,0.473425,0.201976,0.133252,0.026009,0.003255 /)
XCGA_LKT(5,3,1:6)=(/ 0.060580,0.023587,0.010013,0.003940,0.001247,0.000307 /)
XEXT_COEFF_550_LKT(5,3)=20.363000 !rg=0.0143374 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(5,4,1:6)=(/ 384.940000,83.428000,29.872000,6.509600,2.826000,1.345100 /)
XPIZA_LKT(5,4,1:6)=(/ 0.849434,0.601669,0.300677,0.207771,0.043641,0.005553 /)
XCGA_LKT(5,4,1:6)=(/ 0.102130,0.040123,0.017093,0.006740,0.002133,0.000527 /)
XEXT_COEFF_550_LKT(5,4)=24.418000 !rg=0.0143374 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(5,5,1:6)=(/ 615.490000,124.630000,37.484000,7.646000,2.931200,1.351700 /)
XPIZA_LKT(5,5,1:6)=(/ 0.901607,0.727066,0.436525,0.322503,0.076693,0.010072 /)
XCGA_LKT(5,5,1:6)=(/ 0.169203,0.067080,0.028680,0.011337,0.003593,0.000887 /)
XEXT_COEFF_550_LKT(5,5)=32.394000 !rg=0.0143374 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(5,6,1:6)=(/ 967.410000,197.120000,51.357000,9.731800,3.123700,1.363700 /)
XPIZA_LKT(5,6,1:6)=(/ 0.933776,0.821837,0.582066,0.464140,0.131811,0.018287 /)
XCGA_LKT(5,6,1:6)=(/ 0.260543,0.104777,0.044867,0.017773,0.005643,0.001390 /)
XEXT_COEFF_550_LKT(5,6)=46.944000 !rg=0.0143374 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(5,7,1:6)=(/ 1480.200000,316.900000,75.999000,13.489000,3.470100,1.385200 /)
XPIZA_LKT(5,7,1:6)=(/ 0.953537,0.884495,0.711097,0.609591,0.216067,0.032789 /)
XCGA_LKT(5,7,1:6)=(/ 0.365220,0.157213,0.067207,0.026663,0.008483,0.002093 /)
XEXT_COEFF_550_LKT(5,7)=72.787000 !rg=0.0143374 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(5,8,1:6)=(/ 2158.200000,498.380000,118.200000,20.117000,4.081800,1.423000 /)
XPIZA_LKT(5,8,1:6)=(/ 0.965605,0.922710,0.808504,0.734555,0.330563,0.057454 /)
XCGA_LKT(5,8,1:6)=(/ 0.457370,0.231153,0.098510,0.039093,0.012463,0.003077 /)
XEXT_COEFF_550_LKT(5,8)=116.970000 !rg=0.0143374 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(5,9,1:6)=(/ 2990.000000,762.940000,188.320000,31.836000,5.169600,1.490200 /)
XPIZA_LKT(5,9,1:6)=(/ 0.973098,0.946081,0.874987,0.828945,0.467992,0.098437 /)
XCGA_LKT(5,9,1:6)=(/ 0.538093,0.330833,0.143567,0.056810,0.018143,0.004487 /)
XEXT_COEFF_550_LKT(5,9)=190.070000 !rg=0.0143374 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(5,10,1:6)=(/ 3904.800000,1141.000000,296.460000,52.150000,7.087600,1.608500 /)
XPIZA_LKT(5,10,1:6)=(/ 0.977818,0.960985,0.916668,0.892752,0.608335,0.162736 /)
XCGA_LKT(5,10,1:6)=(/ 0.602080,0.431997,0.209310,0.082257,0.026303,0.006517 /)
XEXT_COEFF_550_LKT(5,10)=302.060000 !rg=0.0143374 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(5,11,1:6)=(/ 4814.900000,1616.600000,453.940000,86.832000,10.500000,1.819400 /)
XPIZA_LKT(5,11,1:6)=(/ 0.980747,0.970274,0.942141,0.933273,0.732061,0.257173 /)
XCGA_LKT(5,11,1:6)=(/ 0.651820,0.510420,0.303843,0.119340,0.038113,0.009463 /)
XEXT_COEFF_550_LKT(5,11)=464.620000 !rg=0.0143374 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(5,12,1:6)=(/ 5559.400000,2175.600000,682.950000,142.790000,16.519000,2.194400 /)
XPIZA_LKT(5,12,1:6)=(/ 0.982312,0.976082,0.958347,0.957592,0.826471,0.380925 /)
XCGA_LKT(5,12,1:6)=(/ 0.687053,0.582050,0.414690,0.173963,0.055230,0.013743 /)
XEXT_COEFF_550_LKT(5,12)=701.270000 !rg=0.0143374 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(5,13,1:6)=(/ 5979.900000,2751.000000,987.860000,225.950000,27.009000,2.860800 /)
XPIZA_LKT(5,13,1:6)=(/ 0.982714,0.979698,0.968830,0.971716,0.891112,0.521568 /)
XCGA_LKT(5,13,1:6)=(/ 0.708773,0.635753,0.498880,0.254793,0.080097,0.019960 /)
XEXT_COEFF_550_LKT(5,13)=1011.400000 !rg=0.0143374 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(5,14,1:6)=(/ 5974.800000,3266.100000,1336.700000,343.290000,44.725000,4.037700 /)
XPIZA_LKT(5,14,1:6)=(/ 0.981970,0.981802,0.975108,0.980031,0.932010,0.657367 /)
XCGA_LKT(5,14,1:6)=(/ 0.717780,0.676500,0.569820,0.365303,0.116357,0.028980 /)
XEXT_COEFF_550_LKT(5,14)=1367.300000 !rg=0.0143374 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(5,15,1:6)=(/ 5588.600000,3622.700000,1711.300000,512.520000,73.151000,6.112100 /)
XPIZA_LKT(5,15,1:6)=(/ 0.980139,0.982692,0.979133,0.985447,0.956679,0.770200 /)
XCGA_LKT(5,15,1:6)=(/ 0.717677,0.703213,0.627647,0.470053,0.169743,0.042060 /)
XEXT_COEFF_550_LKT(5,15)=1744.700000 !rg=0.0143374 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(5,16,1:6)=(/ 5033.900000,3737.000000,2062.300000,715.790000,115.210000,9.755000 /)
XPIZA_LKT(5,16,1:6)=(/ 0.977360,0.982451,0.981540,0.988784,0.971103,0.852940 /)
XCGA_LKT(5,16,1:6)=(/ 0.715070,0.716997,0.671143,0.541850,0.248763,0.061040 /)
XEXT_COEFF_550_LKT(5,16)=2096.400000 !rg=0.0143374 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(5,17,1:6)=(/ 4469.700000,3582.300000,2319.600000,945.420000,173.130000,16.000000 /)
XPIZA_LKT(5,17,1:6)=(/ 0.974161,0.981037,0.982669,0.990833,0.979533,0.907770 /)
XCGA_LKT(5,17,1:6)=(/ 0.718090,0.718913,0.700340,0.609053,0.358580,0.088560 /)
XEXT_COEFF_550_LKT(5,17)=2350.900000 !rg=0.0143374 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(5,18,1:6)=(/ 3890.200000,3257.700000,2427.900000,1170.400000,255.870000,26.371000 /)
XPIZA_LKT(5,18,1:6)=(/ 0.970229,0.978515,0.982655,0.992090,0.985026,0.941988 /)
XCGA_LKT(5,18,1:6)=(/ 0.722673,0.714127,0.716280,0.657153,0.470170,0.128680 /)
XEXT_COEFF_550_LKT(5,18)=2452.900000 !rg=0.0143374 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(5,19,1:6)=(/ 3309.700000,2900.300000,2364.000000,1363.100000,358.710000,42.601000 /)
XPIZA_LKT(5,19,1:6)=(/ 0.964441,0.975520,0.981486,0.992797,0.988562,0.962487 /)
XCGA_LKT(5,19,1:6)=(/ 0.725733,0.714007,0.720157,0.693183,0.545250,0.187777 /)
XEXT_COEFF_550_LKT(5,19)=2379.100000 !rg=0.0143374 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(5,20,1:6)=(/ 2786.500000,2568.900000,2162.700000,1482.500000,474.170000,65.686000 /)
XPIZA_LKT(5,20,1:6)=(/ 0.959361,0.972083,0.979225,0.993033,0.990722,0.974385 /)
XCGA_LKT(5,20,1:6)=(/ 0.731550,0.719600,0.715697,0.715443,0.613660,0.274610 /)
XEXT_COEFF_550_LKT(5,20)=2174.400000 !rg=0.0143374 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(6,1,1:6)=(/ 187.820000,52.181000,24.187000,5.665200,2.747600,1.340100 /)
XPIZA_LKT(6,1,1:6)=(/ 0.706604,0.378045,0.145485,0.093552,0.017596,0.002185 /)
XCGA_LKT(6,1,1:6)=(/ 0.036910,0.014313,0.006063,0.002383,0.000753,0.000187 /)
XEXT_COEFF_550_LKT(6,1)=18.475000 !rg=0.0155325 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(6,2,1:6)=(/ 222.740000,57.488000,25.151000,5.807800,2.760900,1.340900 /)
XPIZA_LKT(6,2,1:6)=(/ 0.749736,0.432606,0.176395,0.115039,0.022074,0.002752 /)
XCGA_LKT(6,2,1:6)=(/ 0.046800,0.018170,0.007703,0.003030,0.000960,0.000237 /)
XEXT_COEFF_550_LKT(6,2)=19.480000 !rg=0.0155325 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(6,3,1:6)=(/ 305.950000,70.366000,27.489000,6.154700,2.793200,1.343000 /)
XPIZA_LKT(6,3,1:6)=(/ 0.813546,0.531531,0.242883,0.163344,0.032826,0.004134 /)
XCGA_LKT(6,3,1:6)=(/ 0.070943,0.027647,0.011743,0.004623,0.001463,0.000360 /)
XEXT_COEFF_550_LKT(6,3)=21.922000 !rg=0.0155325 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(6,4,1:6)=(/ 468.710000,97.265000,32.402000,6.886000,2.861000,1.347300 /)
XPIZA_LKT(6,4,1:6)=(/ 0.873915,0.655052,0.352447,0.249755,0.054811,0.007049 /)
XCGA_LKT(6,4,1:6)=(/ 0.119287,0.046987,0.020040,0.007907,0.002503,0.000617 /)
XEXT_COEFF_550_LKT(6,4)=27.065000 !rg=0.0155325 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(6,5,1:6)=(/ 743.990000,148.930000,42.043000,8.328600,2.994300,1.355600 /)
XPIZA_LKT(6,5,1:6)=(/ 0.916675,0.768693,0.494537,0.376462,0.095456,0.012769 /)
XCGA_LKT(6,5,1:6)=(/ 0.196697,0.078480,0.033603,0.013293,0.004217,0.001040 /)
XEXT_COEFF_550_LKT(6,5)=37.172000 !rg=0.0155325 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(6,6,1:6)=(/ 1156.300000,238.360000,59.572000,10.976000,3.238600,1.370800 /)
XPIZA_LKT(6,6,1:6)=(/ 0.942981,0.850246,0.636612,0.523175,0.161646,0.023131 /)
XCGA_LKT(6,6,1:6)=(/ 0.297457,0.122573,0.052550,0.020833,0.006620,0.001633 /)
XEXT_COEFF_550_LKT(6,6)=55.558000 !rg=0.0155325 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(6,7,1:6)=(/ 1740.900000,382.190000,90.524000,15.743000,3.678000,1.398100 /)
XPIZA_LKT(6,7,1:6)=(/ 0.959219,0.902279,0.754618,0.663754,0.259166,0.041306 /)
XCGA_LKT(6,7,1:6)=(/ 0.400810,0.183997,0.078710,0.031250,0.009950,0.002457 /)
XEXT_COEFF_550_LKT(6,7)=88.003000 !rg=0.0155325 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(6,8,1:6)=(/ 2481.300000,594.660000,142.840000,24.139000,4.454200,1.446000 /)
XPIZA_LKT(6,8,1:6)=(/ 0.969072,0.933556,0.839115,0.777161,0.385073,0.071889 /)
XCGA_LKT(6,8,1:6)=(/ 0.489240,0.269360,0.115440,0.045800,0.014613,0.003610 /)
XEXT_COEFF_550_LKT(6,8)=142.690000 !rg=0.0155325 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(6,9,1:6)=(/ 3366.800000,904.440000,227.630000,38.927000,5.834500,1.531200 /)
XPIZA_LKT(6,9,1:6)=(/ 0.975326,0.953023,0.894608,0.858680,0.526985,0.121806 /)
XCGA_LKT(6,9,1:6)=(/ 0.565413,0.374173,0.168473,0.066553,0.021273,0.005267 /)
XEXT_COEFF_550_LKT(6,9)=230.860000 !rg=0.0155325 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(6,10,1:6)=(/ 4291.300000,1326.400000,354.230000,64.331000,8.266500,1.681300 /)
XPIZA_LKT(6,10,1:6)=(/ 0.979189,0.965366,0.928610,0.911881,0.662526,0.197947 /)
XCGA_LKT(6,10,1:6)=(/ 0.623170,0.464750,0.245617,0.096403,0.030830,0.007647 /)
XEXT_COEFF_550_LKT(6,10)=361.670000 !rg=0.0155325 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(6,11,1:6)=(/ 5151.000000,1834.900000,537.420000,106.900000,12.585000,1.948900 /)
XPIZA_LKT(6,11,1:6)=(/ 0.981509,0.972929,0.949587,0.944865,0.774900,0.305220 /)
XCGA_LKT(6,11,1:6)=(/ 0.667197,0.540353,0.350683,0.140053,0.044660,0.011100 /)
XEXT_COEFF_550_LKT(6,11)=551.050000 !rg=0.0155325 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(6,12,1:6)=(/ 5782.700000,2410.800000,801.660000,173.390000,20.175000,2.424700 /)
XPIZA_LKT(6,12,1:6)=(/ 0.982594,0.977782,0.963288,0.964342,0.856542,0.438201 /)
XCGA_LKT(6,12,1:6)=(/ 0.697103,0.605177,0.452997,0.204563,0.064710,0.016117 /)
XEXT_COEFF_550_LKT(6,12)=822.530000 !rg=0.0155325 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(6,13,1:6)=(/ 6037.500000,2977.700000,1123.800000,269.270000,33.287000,3.270300 /)
XPIZA_LKT(6,13,1:6)=(/ 0.982553,0.980736,0.971753,0.975630,0.910509,0.579808 /)
XCGA_LKT(6,13,1:6)=(/ 0.713803,0.653927,0.527503,0.298613,0.093883,0.023403 /)
XEXT_COEFF_550_LKT(6,13)=1149.800000 !rg=0.0155325 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(6,14,1:6)=(/ 5860.200000,3439.700000,1493.400000,406.700000,55.032000,4.762100 /)
XPIZA_LKT(6,14,1:6)=(/ 0.981366,0.982309,0.977015,0.982544,0.943844,0.707842 /)
XCGA_LKT(6,14,1:6)=(/ 0.718567,0.688750,0.596357,0.413430,0.136560,0.033967 /)
XEXT_COEFF_550_LKT(6,14)=1525.100000 !rg=0.0155325 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(6,15,1:6)=(/ 5383.700000,3701.600000,1861.800000,594.030000,88.804000,7.386500 /)
XPIZA_LKT(6,15,1:6)=(/ 0.979058,0.982737,0.980258,0.987050,0.963620,0.808329 /)
XCGA_LKT(6,15,1:6)=(/ 0.716103,0.710130,0.646627,0.501373,0.199620,0.049287 /)
XEXT_COEFF_550_LKT(6,15)=1897.000000 !rg=0.0155325 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(6,16,1:6)=(/ 4805.200000,3702.400000,2182.500000,807.520000,137.040000,11.972000 /)
XPIZA_LKT(6,16,1:6)=(/ 0.976195,0.982013,0.982108,0.989730,0.975121,0.878873 /)
XCGA_LKT(6,16,1:6)=(/ 0.716263,0.718713,0.684537,0.571013,0.291983,0.071533 /)
XEXT_COEFF_550_LKT(6,16)=2217.100000 !rg=0.0155325 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(6,17,1:6)=(/ 4246.200000,3462.000000,2385.600000,1039.800000,204.180000,19.750000 /)
XPIZA_LKT(6,17,1:6)=(/ 0.972711,0.980127,0.982770,0.991439,0.982083,0.924221 /)
XCGA_LKT(6,17,1:6)=(/ 0.720983,0.717137,0.708377,0.630437,0.409203,0.103850 /)
XEXT_COEFF_550_LKT(6,17)=2415.300000 !rg=0.0155325 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(6,18,1:6)=(/ 3649.500000,3109.400000,2421.400000,1257.200000,297.710000,32.388000 /)
XPIZA_LKT(6,18,1:6)=(/ 0.968018,0.977376,0.982311,0.992455,0.986735,0.951938 /)
XCGA_LKT(6,18,1:6)=(/ 0.723677,0.713797,0.719237,0.673867,0.504823,0.151163 /)
XEXT_COEFF_550_LKT(6,18)=2442.000000 !rg=0.0155325 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(6,19,1:6)=(/ 3082.800000,2769.700000,2291.000000,1423.800000,404.770000,51.452000 /)
XPIZA_LKT(6,19,1:6)=(/ 0.962630,0.974170,0.980667,0.992958,0.989571,0.968299 /)
XCGA_LKT(6,19,1:6)=(/ 0.727313,0.716390,0.718767,0.703837,0.574277,0.221093 /)
XEXT_COEFF_550_LKT(6,19)=2304.700000 !rg=0.0155325 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(6,20,1:6)=(/ 2578.000000,2423.000000,2067.000000,1500.900000,524.570000,77.739000 /)
XPIZA_LKT(6,20,1:6)=(/ 0.957392,0.970310,0.977974,0.993000,0.991395,0.977791 /)
XCGA_LKT(6,20,1:6)=(/ 0.734070,0.721477,0.713567,0.720567,0.637487,0.321877 /)
XEXT_COEFF_550_LKT(6,20)=2071.400000 !rg=0.0155325 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(7,1,1:6)=(/ 224.510000,57.735000,25.198000,5.814400,2.761500,1.341000 /)
XPIZA_LKT(7,1,1:6)=(/ 0.751093,0.434492,0.177583,0.115901,0.022260,0.002776 /)
XCGA_LKT(7,1,1:6)=(/ 0.043263,0.016787,0.007113,0.002797,0.000883,0.000217 /)
XEXT_COEFF_550_LKT(7,1)=19.526000 !rg=0.0168272 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(7,2,1:6)=(/ 268.580000,64.465000,26.418000,5.995400,2.778400,1.342100 /)
XPIZA_LKT(7,2,1:6)=(/ 0.789211,0.490554,0.213548,0.141701,0.027890,0.003496 /)
XCGA_LKT(7,2,1:6)=(/ 0.054850,0.021307,0.009037,0.003557,0.001127,0.000277 /)
XEXT_COEFF_550_LKT(7,2)=20.801000 !rg=0.0168272 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(7,3,1:6)=(/ 372.630000,80.781000,29.382000,6.435500,2.819300,1.344600 /)
XPIZA_LKT(7,3,1:6)=(/ 0.844110,0.588473,0.288963,0.198642,0.041351,0.005250 /)
XCGA_LKT(7,3,1:6)=(/ 0.083077,0.032403,0.013773,0.005427,0.001717,0.000423 /)
XEXT_COEFF_550_LKT(7,3)=23.899000 !rg=0.0168272 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(7,4,1:6)=(/ 571.050000,114.720000,35.608000,7.363700,2.905300,1.350000 /)
XPIZA_LKT(7,4,1:6)=(/ 0.894276,0.704317,0.407731,0.296977,0.068628,0.008944 /)
XCGA_LKT(7,4,1:6)=(/ 0.139217,0.055003,0.023487,0.009273,0.002940,0.000723 /)
XEXT_COEFF_550_LKT(7,4)=30.421000 !rg=0.0168272 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(7,5,1:6)=(/ 896.760000,179.190000,47.811000,9.195100,3.074400,1.360600 /)
XPIZA_LKT(7,5,1:6)=(/ 0.929079,0.804999,0.552360,0.433560,0.118207,0.016176 /)
XCGA_LKT(7,5,1:6)=(/ 0.227547,0.091780,0.039367,0.015587,0.004947,0.001220 /)
XEXT_COEFF_550_LKT(7,5)=43.219000 !rg=0.0168272 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(7,6,1:6)=(/ 1375.400000,288.450000,69.917000,12.556000,3.384200,1.379900 /)
XPIZA_LKT(7,6,1:6)=(/ 0.950585,0.874019,0.687351,0.581396,0.196670,0.029217 /)
XCGA_LKT(7,6,1:6)=(/ 0.334430,0.143333,0.061537,0.024423,0.007767,0.001917 /)
XEXT_COEFF_550_LKT(7,6)=66.406000 !rg=0.0168272 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(7,7,1:6)=(/ 2028.100000,459.200000,108.620000,18.600000,3.941900,1.414400 /)
XPIZA_LKT(7,7,1:6)=(/ 0.963844,0.916846,0.792814,0.713694,0.307422,0.051913 /)
XCGA_LKT(7,7,1:6)=(/ 0.434533,0.215010,0.092180,0.036617,0.011670,0.002883 /)
XEXT_COEFF_550_LKT(7,7)=106.940000 !rg=0.0168272 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(7,8,1:6)=(/ 2834.400000,707.700000,172.890000,29.223000,4.926800,1.475200 /)
XPIZA_LKT(7,8,1:6)=(/ 0.971978,0.942558,0.864843,0.814371,0.442508,0.089596 /)
XCGA_LKT(7,8,1:6)=(/ 0.520550,0.310887,0.135310,0.053657,0.017137,0.004237 /)
XEXT_COEFF_550_LKT(7,8)=173.990000 !rg=0.0168272 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(7,9,1:6)=(/ 3748.700000,1066.700000,274.040000,47.824000,6.677700,1.583200 /)
XPIZA_LKT(7,9,1:6)=(/ 0.977142,0.958838,0.910637,0.883649,0.585052,0.149777 /)
XCGA_LKT(7,9,1:6)=(/ 0.589377,0.413197,0.197720,0.077970,0.024937,0.006177 /)
XEXT_COEFF_550_LKT(7,9)=278.870000 !rg=0.0168272 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(7,10,1:6)=(/ 4665.100000,1523.100000,421.120000,79.385000,9.759700,1.773600 /)
XPIZA_LKT(7,10,1:6)=(/ 0.980321,0.968885,0.938356,0.927517,0.712522,0.238559 /)
XCGA_LKT(7,10,1:6)=(/ 0.642177,0.494943,0.287060,0.113023,0.036130,0.008967 /)
XEXT_COEFF_550_LKT(7,10)=430.760000 !rg=0.0168272 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(7,11,1:6)=(/ 5451.700000,2068.700000,635.750000,130.990000,15.220000,2.113200 /)
XPIZA_LKT(7,11,1:6)=(/ 0.982107,0.975182,0.955914,0.954164,0.812365,0.357836 /)
XCGA_LKT(7,11,1:6)=(/ 0.680453,0.568920,0.396313,0.164477,0.052327,0.013017 /)
XEXT_COEFF_550_LKT(7,11)=652.660000 !rg=0.0168272 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(7,12,1:6)=(/ 5942.100000,2645.100000,928.550000,208.790000,24.761000,2.717100 /)
XPIZA_LKT(7,12,1:6)=(/ 0.982726,0.979143,0.967300,0.969709,0.881832,0.497035 /)
XCGA_LKT(7,12,1:6)=(/ 0.705027,0.625690,0.484700,0.240500,0.075823,0.018897 /)
XEXT_COEFF_550_LKT(7,12)=951.180000 !rg=0.0168272 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(7,13,1:6)=(/ 6020.800000,3182.300000,1270.000000,319.590000,41.048000,3.789500 /)
XPIZA_LKT(7,13,1:6)=(/ 0.982212,0.981514,0.974165,0.978836,0.926391,0.635708 /)
XCGA_LKT(7,13,1:6)=(/ 0.716693,0.669357,0.556960,0.346707,0.110080,0.027433 /)
XEXT_COEFF_550_LKT(7,13)=1299.400000 !rg=0.0168272 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(7,14,1:6)=(/ 5701.000000,3575.100000,1645.600000,479.800000,67.416000,5.679400 /)
XPIZA_LKT(7,14,1:6)=(/ 0.980551,0.982612,0.978554,0.984675,0.953353,0.753421 /)
XCGA_LKT(7,14,1:6)=(/ 0.717917,0.698613,0.618237,0.454850,0.160390,0.039807 /)
XEXT_COEFF_550_LKT(7,14)=1678.200000 !rg=0.0168272 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(7,15,1:6)=(/ 5155.700000,3735.700000,2003.700000,677.880000,106.880000,8.993800 /)
XPIZA_LKT(7,15,1:6)=(/ 0.978020,0.982565,0.981189,0.988315,0.969137,0.841152 /)
XCGA_LKT(7,15,1:6)=(/ 0.716123,0.714690,0.663787,0.529580,0.234803,0.057753 /)
XEXT_COEFF_550_LKT(7,15)=2037.900000 !rg=0.0168272 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(7,16,1:6)=(/ 4589.500000,3629.100000,2284.300000,905.040000,162.030000,14.746000 /)
XPIZA_LKT(7,16,1:6)=(/ 0.974899,0.981374,0.982533,0.990538,0.978383,0.900454 /)
XCGA_LKT(7,16,1:6)=(/ 0.718660,0.718887,0.695700,0.598690,0.340337,0.083840 /)
XEXT_COEFF_550_LKT(7,16)=2316.300000 !rg=0.0168272 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(7,17,1:6)=(/ 4011.500000,3324.200000,2421.100000,1132.700000,240.490000,24.360000 /)
XPIZA_LKT(7,17,1:6)=(/ 0.970830,0.979110,0.982726,0.991916,0.984271,0.937600 /)
XCGA_LKT(7,17,1:6)=(/ 0.722103,0.715917,0.714060,0.649413,0.454863,0.121847 /)
XEXT_COEFF_550_LKT(7,17)=2447.200000 !rg=0.0168272 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(7,18,1:6)=(/ 3403.500000,2974.300000,2384.900000,1333.200000,341.020000,39.544000 /)
XPIZA_LKT(7,18,1:6)=(/ 0.965591,0.976040,0.981755,0.992708,0.988101,0.959896 /)
XCGA_LKT(7,18,1:6)=(/ 0.724093,0.713837,0.719757,0.687703,0.534097,0.177753 /)
XEXT_COEFF_550_LKT(7,18)=2403.000000 !rg=0.0168272 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(7,19,1:6)=(/ 2867.300000,2632.300000,2205.700000,1467.600000,454.450000,61.562000 /)
XPIZA_LKT(7,19,1:6)=(/ 0.961163,0.972642,0.979616,0.993020,0.990426,0.972911 /)
XCGA_LKT(7,19,1:6)=(/ 0.732103,0.718633,0.716407,0.712130,0.603430,0.260270 /)
XEXT_COEFF_550_LKT(7,19)=2213.200000 !rg=0.0168272 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(7,20,1:6)=(/ 2394.700000,2270.700000,1969.700000,1501.500000,573.360000,91.629000 /)
XPIZA_LKT(7,20,1:6)=(/ 0.954022,0.968499,0.976666,0.992866,0.991927,0.980592 /)
XCGA_LKT(7,20,1:6)=(/ 0.735900,0.722903,0.711993,0.723560,0.657597,0.373333 /)
XEXT_COEFF_550_LKT(7,20)=1977.800000 !rg=0.0168272 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(8,1,1:6)=(/ 270.930000,64.778000,26.477000,6.003700,2.779200,1.342100 /)
XPIZA_LKT(8,1,1:6)=(/ 0.790436,0.492444,0.214911,0.142727,0.028123,0.003527 /)
XCGA_LKT(8,1,1:6)=(/ 0.050713,0.019683,0.008347,0.003283,0.001040,0.000257 /)
XEXT_COEFF_550_LKT(8,1)=20.860000 !rg=0.0182298 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(8,2,1:6)=(/ 326.400000,73.309000,28.024000,6.233300,2.800600,1.343500 /)
XPIZA_LKT(8,2,1:6)=(/ 0.823492,0.548471,0.256000,0.173305,0.035180,0.004441 /)
XCGA_LKT(8,2,1:6)=(/ 0.064287,0.024980,0.010600,0.004173,0.001320,0.000327 /)
XEXT_COEFF_550_LKT(8,2)=22.476000 !rg=0.0182298 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(8,3,1:6)=(/ 455.680000,93.966000,31.780000,6.791800,2.852400,1.346700 /)
XPIZA_LKT(8,3,1:6)=(/ 0.869932,0.642772,0.339713,0.239336,0.051967,0.006665 /)
XCGA_LKT(8,3,1:6)=(/ 0.097287,0.037967,0.016153,0.006367,0.002017,0.000497 /)
XEXT_COEFF_550_LKT(8,3)=26.407000 !rg=0.0182298 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(8,4,1:6)=(/ 694.960000,136.680000,39.669000,7.970000,2.961400,1.353600 /)
XPIZA_LKT(8,4,1:6)=(/ 0.911075,0.748717,0.465211,0.348894,0.085605,0.011342 /)
XCGA_LKT(8,4,1:6)=(/ 0.162287,0.064367,0.027527,0.010877,0.003450,0.000850 /)
XEXT_COEFF_550_LKT(8,4)=34.676000 !rg=0.0182298 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(8,5,1:6)=(/ 1076.900000,216.610000,55.100000,10.295000,3.175900,1.366900 /)
XPIZA_LKT(8,5,1:6)=(/ 0.939281,0.836103,0.608399,0.492339,0.145487,0.020473 /)
XCGA_LKT(8,5,1:6)=(/ 0.261123,0.107277,0.046103,0.018277,0.005807,0.001430 /)
XEXT_COEFF_550_LKT(8,5)=50.863000 !rg=0.0182298 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(8,6,1:6)=(/ 1624.400000,348.670000,82.904000,14.559000,3.569000,1.391300 /)
XPIZA_LKT(8,6,1:6)=(/ 0.956816,0.893704,0.733388,0.637234,0.237091,0.036843 /)
XCGA_LKT(8,6,1:6)=(/ 0.370480,0.167473,0.072043,0.028623,0.009110,0.002250 /)
XEXT_COEFF_550_LKT(8,6)=80.014000 !rg=0.0182298 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(8,7,1:6)=(/ 2344.200000,549.690000,131.010000,22.219000,4.276700,1.435100 /)
XPIZA_LKT(8,7,1:6)=(/ 0.967646,0.928804,0.825695,0.758656,0.360206,0.065053 /)
XCGA_LKT(8,7,1:6)=(/ 0.467803,0.250290,0.107953,0.042900,0.013687,0.003383 /)
XEXT_COEFF_550_LKT(8,7)=130.330000 !rg=0.0182298 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(8,8,1:6)=(/ 3208.000000,840.630000,209.080000,35.631000,5.526400,1.512200 /)
XPIZA_LKT(8,8,1:6)=(/ 0.974406,0.950113,0.886179,0.846292,0.501377,0.111129 /)
XCGA_LKT(8,8,1:6)=(/ 0.549200,0.352807,0.158630,0.062853,0.020093,0.004973 /)
XEXT_COEFF_550_LKT(8,8)=211.590000 !rg=0.0182298 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(8,9,1:6)=(/ 4139.800000,1245.000000,328.170000,58.932000,7.746700,1.649200 /)
XPIZA_LKT(8,9,1:6)=(/ 0.978641,0.963591,0.923674,0.904357,0.640645,0.182807 /)
XCGA_LKT(8,9,1:6)=(/ 0.611613,0.447110,0.231837,0.091357,0.029230,0.007247 /)
XEXT_COEFF_550_LKT(8,9)=334.790000 !rg=0.0182298 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(8,10,1:6)=(/ 5016.700000,1735.200000,499.330000,97.797000,11.649000,1.890800 /)
XPIZA_LKT(8,10,1:6)=(/ 0.981189,0.971772,0.946445,0.940188,0.757562,0.284487 /)
XCGA_LKT(8,10,1:6)=(/ 0.658710,0.525257,0.332183,0.132573,0.042337,0.010520 /)
XEXT_COEFF_550_LKT(8,10)=511.710000 !rg=0.0182298 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(8,11,1:6)=(/ 5702.200000,2304.500000,748.440000,159.460000,18.539000,2.321800 /)
XPIZA_LKT(8,11,1:6)=(/ 0.982488,0.977044,0.961257,0.961583,0.844547,0.414011 /)
XCGA_LKT(8,11,1:6)=(/ 0.691550,0.593427,0.436177,0.193280,0.061307,0.015267 /)
XEXT_COEFF_550_LKT(8,11)=768.110000 !rg=0.0182298 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(8,12,1:6)=(/ 6036.600000,2877.900000,1060.900000,249.500000,30.485000,3.088000 /)
XPIZA_LKT(8,12,1:6)=(/ 0.982650,0.980280,0.970500,0.974003,0.902839,0.555788 /)
XCGA_LKT(8,12,1:6)=(/ 0.710960,0.644810,0.513737,0.282057,0.088860,0.022157 /)
XEXT_COEFF_550_LKT(8,12)=1085.800000 !rg=0.0182298 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(8,13,1:6)=(/ 5941.600000,3369.300000,1425.200000,378.940000,50.545000,4.447700 /)
XPIZA_LKT(8,13,1:6)=(/ 0.981685,0.982090,0.976222,0.981536,0.939284,0.687954 /)
XCGA_LKT(8,13,1:6)=(/ 0.718213,0.682593,0.584770,0.395160,0.129147,0.032157 /)
XEXT_COEFF_550_LKT(8,13)=1456.500000 !rg=0.0182298 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(8,14,1:6)=(/ 5500.000000,3675.300000,1796.600000,559.110000,82.055000,6.839400 /)
XPIZA_LKT(8,14,1:6)=(/ 0.979679,0.982725,0.979779,0.986422,0.960951,0.793700 /)
XCGA_LKT(8,14,1:6)=(/ 0.717653,0.706633,0.637790,0.488207,0.188530,0.046647 /)
XEXT_COEFF_550_LKT(8,14)=1831.400000 !rg=0.0182298 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(8,15,1:6)=(/ 4940.000000,3723.800000,2129.500000,766.710000,127.500000,11.016000 /)
XPIZA_LKT(8,15,1:6)=(/ 0.976838,0.982219,0.981855,0.989336,0.973533,0.868968 /)
XCGA_LKT(8,15,1:6)=(/ 0.716987,0.717513,0.678127,0.558557,0.275783,0.067677 /)
XEXT_COEFF_550_LKT(8,15)=2164.200000 !rg=0.0182298 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(8,16,1:6)=(/ 4367.000000,3523.200000,2360.000000,1000.700000,191.190000,18.196000 /)
XPIZA_LKT(8,16,1:6)=(/ 0.973368,0.980565,0.982727,0.991202,0.981110,0.918224 /)
XCGA_LKT(8,16,1:6)=(/ 0.720160,0.718123,0.704523,0.621583,0.390873,0.098297 /)
XEXT_COEFF_550_LKT(8,16)=2390.900000 !rg=0.0182298 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(8,17,1:6)=(/ 3754.200000,3185.700000,2428.300000,1223.000000,281.040000,29.958000 /)
XPIZA_LKT(8,17,1:6)=(/ 0.969080,0.977862,0.982448,0.992314,0.986112,0.948399 /)
XCGA_LKT(8,17,1:6)=(/ 0.723673,0.714067,0.717660,0.667023,0.492053,0.143080 /)
XEXT_COEFF_550_LKT(8,17)=2450.800000 !rg=0.0182298 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(8,18,1:6)=(/ 3177.700000,2833.500000,2324.200000,1400.000000,385.920000,47.902000 /)
XPIZA_LKT(8,18,1:6)=(/ 0.963867,0.974858,0.981001,0.992892,0.989190,0.966227 /)
XCGA_LKT(8,18,1:6)=(/ 0.728737,0.715977,0.719050,0.699277,0.562847,0.209210 /)
XEXT_COEFF_550_LKT(8,18)=2337.200000 !rg=0.0182298 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(8,19,1:6)=(/ 2664.400000,2489.700000,2108.200000,1494.400000,505.160000,73.033000 /)
XPIZA_LKT(8,19,1:6)=(/ 0.958284,0.971106,0.978484,0.993014,0.991151,0.976595 /)
XCGA_LKT(8,19,1:6)=(/ 0.733197,0.721307,0.714010,0.718347,0.628820,0.305523 /)
XEXT_COEFF_550_LKT(8,19)=2119.500000 !rg=0.0182298 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(8,20,1:6)=(/ 2214.500000,2127.900000,1877.100000,1481.400000,622.150000,107.960000 /)
XPIZA_LKT(8,20,1:6)=(/ 0.950658,0.965893,0.975428,0.992657,0.992368,0.982973 /)
XCGA_LKT(8,20,1:6)=(/ 0.737173,0.724117,0.713430,0.724380,0.676287,0.424487 /)
XEXT_COEFF_550_LKT(8,20)=1883.700000 !rg=0.0182298 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(9,1,1:6)=(/ 329.580000,73.707000,28.098000,6.243800,2.801600,1.343500 /)
XPIZA_LKT(9,1,1:6)=(/ 0.824607,0.550320,0.257532,0.174510,0.035472,0.004479 /)
XCGA_LKT(9,1,1:6)=(/ 0.059450,0.023080,0.009790,0.003853,0.001220,0.000300 /)
XEXT_COEFF_550_LKT(9,1)=22.551000 !rg=0.0197494 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(9,2,1:6)=(/ 399.130000,84.521000,30.059000,6.535100,2.828700,1.345300 /)
XPIZA_LKT(9,2,1:6)=(/ 0.852798,0.604787,0.303513,0.210195,0.044288,0.005639 /)
XCGA_LKT(9,2,1:6)=(/ 0.075360,0.029283,0.012433,0.004897,0.001550,0.000380 /)
XEXT_COEFF_550_LKT(9,2)=24.602000 !rg=0.0197494 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(9,3,1:6)=(/ 558.400000,110.650000,34.820000,7.243900,2.894400,1.349400 /)
XPIZA_LKT(9,3,1:6)=(/ 0.891482,0.693236,0.394253,0.285340,0.065117,0.008458 /)
XCGA_LKT(9,3,1:6)=(/ 0.113940,0.044483,0.018940,0.007467,0.002367,0.000583 /)
XEXT_COEFF_550_LKT(9,3)=29.588000 !rg=0.0197494 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(9,4,1:6)=(/ 843.500000,164.210000,44.812000,8.739700,3.032500,1.358000 /)
XPIZA_LKT(9,4,1:6)=(/ 0.924855,0.787888,0.523339,0.404563,0.106289,0.014374 /)
XCGA_LKT(9,4,1:6)=(/ 0.188790,0.075293,0.032253,0.012757,0.004047,0.000997 /)
XEXT_COEFF_550_LKT(9,4)=40.067000 !rg=0.0197494 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(9,5,1:6)=(/ 1286.900000,262.520000,64.296000,11.691000,3.304700,1.374900 /)
XPIZA_LKT(9,5,1:6)=(/ 0.947652,0.862358,0.661261,0.551175,0.177767,0.025879 /)
XCGA_LKT(9,5,1:6)=(/ 0.296533,0.125307,0.053980,0.021423,0.006810,0.001680 /)
XEXT_COEFF_550_LKT(9,5)=60.506000 !rg=0.0197494 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(9,6,1:6)=(/ 1903.100000,420.400000,99.131000,17.099000,3.803500,1.405800 /)
XPIZA_LKT(9,6,1:6)=(/ 0.961910,0.909909,0.774224,0.689379,0.282837,0.046361 /)
XCGA_LKT(9,6,1:6)=(/ 0.405893,0.195370,0.084333,0.033540,0.010687,0.002640 /)
XEXT_COEFF_550_LKT(9,6)=97.004000 !rg=0.0197494 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(9,7,1:6)=(/ 2691.600000,656.200000,158.440000,26.794000,4.701500,1.461300 /)
XPIZA_LKT(9,7,1:6)=(/ 0.970822,0.938705,0.853528,0.798272,0.416494,0.081228 /)
XCGA_LKT(9,7,1:6)=(/ 0.500423,0.288947,0.126427,0.050257,0.016050,0.003967 /)
XEXT_COEFF_550_LKT(9,7)=158.930000 !rg=0.0197494 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(9,8,1:6)=(/ 3592.200000,993.840000,252.080000,43.684000,6.287000,1.559100 /)
XPIZA_LKT(9,8,1:6)=(/ 0.976397,0.956423,0.903693,0.873257,0.560042,0.137042 /)
XCGA_LKT(9,8,1:6)=(/ 0.574737,0.391747,0.185960,0.073623,0.023553,0.005833 /)
XEXT_COEFF_550_LKT(9,8)=256.150000 !rg=0.0197494 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(9,9,1:6)=(/ 4525.400000,1436.000000,391.080000,72.697000,9.101300,1.732900 /)
XPIZA_LKT(9,9,1:6)=(/ 0.979888,0.967426,0.934316,0.921348,0.692475,0.221188 /)
XCGA_LKT(9,9,1:6)=(/ 0.631887,0.478260,0.270930,0.107067,0.034257,0.008500 /)
XEXT_COEFF_550_LKT(9,9)=399.780000 !rg=0.0197494 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(9,10,1:6)=(/ 5341.100000,1964.400000,591.470000,120.020000,14.037000,2.039600 /)
XPIZA_LKT(9,10,1:6)=(/ 0.981877,0.974213,0.953284,0.950381,0.797275,0.335290 /)
XCGA_LKT(9,10,1:6)=(/ 0.673177,0.554703,0.377370,0.155593,0.049603,0.012340 /)
XEXT_COEFF_550_LKT(9,10)=607.000000 !rg=0.0197494 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(9,11,1:6)=(/ 5897.000000,2539.500000,870.910000,192.590000,22.707000,2.586500 /)
XPIZA_LKT(9,11,1:6)=(/ 0.982705,0.978539,0.965629,0.967489,0.871765,0.472388 /)
XCGA_LKT(9,11,1:6)=(/ 0.700720,0.614863,0.469427,0.227120,0.071830,0.017903 /)
XEXT_COEFF_550_LKT(9,11)=892.550000 !rg=0.0197494 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(9,12,1:6)=(/ 6057.400000,3091.800000,1202.700000,296.650000,37.580000,3.558400 /)
XPIZA_LKT(9,12,1:6)=(/ 0.982415,0.981172,0.973113,0.977497,0.920105,0.612838 /)
XCGA_LKT(9,12,1:6)=(/ 0.715080,0.661347,0.543283,0.328377,0.104167,0.025973 /)
XEXT_COEFF_550_LKT(9,12)=1230.800000 !rg=0.0197494 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(9,13,1:6)=(/ 5803.500000,3521.100000,1578.800000,448.140000,62.018000,5.281400 /)
XPIZA_LKT(9,13,1:6)=(/ 0.981021,0.982495,0.977911,0.983831,0.949674,0.735587 /)
XCGA_LKT(9,13,1:6)=(/ 0.718760,0.693463,0.608083,0.438530,0.151617,0.037687 /)
XEXT_COEFF_550_LKT(9,13)=1610.900000 !rg=0.0197494 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(9,14,1:6)=(/ 5293.700000,3729.300000,1943.900000,641.300000,99.059000,8.303800 /)
XPIZA_LKT(9,14,1:6)=(/ 0.978586,0.982655,0.980802,0.987808,0.966996,0.828616 /)
XCGA_LKT(9,14,1:6)=(/ 0.717167,0.712120,0.655910,0.517073,0.221690,0.054660 /)
XEXT_COEFF_550_LKT(9,14)=1978.500000 !rg=0.0197494 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(9,15,1:6)=(/ 4716.900000,3669.900000,2241.400000,862.480000,151.070000,13.550000 /)
XPIZA_LKT(9,15,1:6)=(/ 0.975604,0.981687,0.982344,0.990203,0.977081,0.892230 /)
XCGA_LKT(9,15,1:6)=(/ 0.718277,0.718653,0.690273,0.587230,0.322233,0.079310 /)
XEXT_COEFF_550_LKT(9,15)=2274.400000 !rg=0.0197494 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(9,16,1:6)=(/ 4126.700000,3396.500000,2409.600000,1093.600000,225.460000,22.451000 /)
XPIZA_LKT(9,16,1:6)=(/ 0.971817,0.979551,0.982755,0.991723,0.983443,0.932724 /)
XCGA_LKT(9,16,1:6)=(/ 0.722753,0.716200,0.711347,0.641163,0.438337,0.115303 /)
XEXT_COEFF_550_LKT(9,16)=2436.700000 !rg=0.0197494 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(9,17,1:6)=(/ 3517.500000,3040.000000,2403.600000,1302.500000,323.640000,36.656000 /)
XPIZA_LKT(9,17,1:6)=(/ 0.966296,0.976714,0.981996,0.992610,0.987599,0.957057 /)
XCGA_LKT(9,17,1:6)=(/ 0.725263,0.714523,0.719247,0.681840,0.522587,0.168173 /)
XEXT_COEFF_550_LKT(9,17)=2422.500000 !rg=0.0197494 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(9,18,1:6)=(/ 2959.500000,2700.300000,2242.700000,1449.600000,434.280000,57.497000 /)
XPIZA_LKT(9,18,1:6)=(/ 0.961918,0.973460,0.980097,0.992995,0.990099,0.971253 /)
XCGA_LKT(9,18,1:6)=(/ 0.730633,0.718990,0.717253,0.708347,0.592387,0.246287 /)
XEXT_COEFF_550_LKT(9,18)=2256.000000 !rg=0.0197494 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(9,19,1:6)=(/ 2462.700000,2342.200000,2012.300000,1501.800000,554.150000,86.180000 /)
XPIZA_LKT(9,19,1:6)=(/ 0.955088,0.969070,0.977277,0.992917,0.991731,0.979599 /)
XCGA_LKT(9,19,1:6)=(/ 0.733940,0.722417,0.713330,0.721950,0.649830,0.355740 /)
XEXT_COEFF_550_LKT(9,19)=2019.100000 !rg=0.0197494 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(9,20,1:6)=(/ 2050.600000,1984.300000,1793.600000,1444.900000,667.020000,126.990000 /)
XPIZA_LKT(9,20,1:6)=(/ 0.946824,0.963550,0.973952,0.992360,0.992725,0.985024 /)
XCGA_LKT(9,20,1:6)=(/ 0.741877,0.724877,0.714890,0.723477,0.692750,0.469380 /)
XEXT_COEFF_550_LKT(9,20)=1797.400000 !rg=0.0197494 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(10,1,1:6)=(/ 403.520000,85.028000,30.151000,6.548300,2.830000,1.345300 /)
XPIZA_LKT(10,1,1:6)=(/ 0.853826,0.606542,0.305190,0.211581,0.044651,0.005688 /)
XCGA_LKT(10,1,1:6)=(/ 0.069703,0.027060,0.011483,0.004520,0.001430,0.000353 /)
XEXT_COEFF_550_LKT(10,1)=24.695000 !rg=0.0213956 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(10,2,1:6)=(/ 490.200000,98.727000,32.637000,6.918000,2.864300,1.347500 /)
XPIZA_LKT(10,2,1:6)=(/ 0.877511,0.658097,0.355471,0.252484,0.055611,0.007158 /)
XCGA_LKT(10,2,1:6)=(/ 0.088370,0.034327,0.014583,0.005743,0.001820,0.000447 /)
XEXT_COEFF_550_LKT(10,2)=27.298000 !rg=0.0213956 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(10,3,1:6)=(/ 684.270000,131.720000,38.672000,7.817800,2.947500,1.352700 /)
XPIZA_LKT(10,3,1:6)=(/ 0.909284,0.739026,0.451340,0.336200,0.081302,0.010728 /)
XCGA_LKT(10,3,1:6)=(/ 0.133477,0.052107,0.022203,0.008760,0.002777,0.000683 /)
XEXT_COEFF_550_LKT(10,3)=33.623000 !rg=0.0213956 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(10,4,1:6)=(/ 1019.600000,198.580000,51.321000,9.716700,3.122800,1.363700 /)
XPIZA_LKT(10,4,1:6)=(/ 0.936112,0.821798,0.580494,0.462676,0.131236,0.018201 /)
XCGA_LKT(10,4,1:6)=(/ 0.218873,0.088037,0.037783,0.014960,0.004747,0.001170 /)
XEXT_COEFF_550_LKT(10,4)=46.892000 !rg=0.0213956 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(10,5,1:6)=(/ 1528.000000,318.330000,75.867000,13.462000,3.468000,1.385100 /)
XPIZA_LKT(10,5,1:6)=(/ 0.954499,0.884259,0.709854,0.608423,0.215364,0.032664 /)
XCGA_LKT(10,5,1:6)=(/ 0.333157,0.146237,0.063187,0.025110,0.007990,0.001970 /)
XEXT_COEFF_550_LKT(10,5)=72.637000 !rg=0.0213956 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(10,6,1:6)=(/ 2213.700000,505.340000,119.280000,20.318000,4.101000,1.424200 /)
XPIZA_LKT(10,6,1:6)=(/ 0.966107,0.923239,0.809711,0.736868,0.333468,0.058185 /)
XCGA_LKT(10,6,1:6)=(/ 0.441190,0.227137,0.098697,0.039293,0.012533,0.003097 /)
XEXT_COEFF_550_LKT(10,6)=118.070000 !rg=0.0213956 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(10,7,1:6)=(/ 3063.000000,781.440000,191.660000,32.567000,5.240500,1.494600 /)
XPIZA_LKT(10,7,1:6)=(/ 0.973469,0.946965,0.876754,0.832518,0.474912,0.100979 /)
XCGA_LKT(10,7,1:6)=(/ 0.530703,0.328743,0.148057,0.058863,0.018817,0.004657 /)
XEXT_COEFF_550_LKT(10,7)=193.500000 !rg=0.0213956 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(10,8,1:6)=(/ 3988.600000,1164.100000,302.580000,53.753000,7.251300,1.618600 /)
XPIZA_LKT(10,8,1:6)=(/ 0.978045,0.961593,0.917987,0.895730,0.616876,0.167829 /)
XCGA_LKT(10,8,1:6)=(/ 0.598450,0.426563,0.217793,0.086240,0.027607,0.006843 /)
XEXT_COEFF_550_LKT(10,8)=308.370000 !rg=0.0213956 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(10,9,1:6)=(/ 4893.100000,1642.600000,464.610000,89.592000,10.816000,1.839200 /)
XPIZA_LKT(10,9,1:6)=(/ 0.980867,0.970567,0.943117,0.935163,0.739612,0.264960 /)
XCGA_LKT(10,9,1:6)=(/ 0.649710,0.509183,0.313967,0.125523,0.040140,0.009973 /)
XEXT_COEFF_550_LKT(10,9)=475.880000 !rg=0.0213956 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(10,10,1:6)=(/ 5619.100000,2199.700000,697.960000,146.430000,17.047000,2.228300 /)
XPIZA_LKT(10,10,1:6)=(/ 0.982359,0.976245,0.959058,0.958532,0.831630,0.390118 /)
XCGA_LKT(10,10,1:6)=(/ 0.685443,0.580550,0.418143,0.182703,0.058110,0.014470 /)
XEXT_COEFF_550_LKT(10,10)=716.370000 !rg=0.0213956 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(10,11,1:6)=(/ 6026.600000,2776.500000,999.730000,230.810000,27.919000,2.922300 /)
XPIZA_LKT(10,11,1:6)=(/ 0.982730,0.979783,0.969128,0.972206,0.894484,0.531377 /)
XCGA_LKT(10,11,1:6)=(/ 0.707707,0.634903,0.499200,0.266403,0.084167,0.020990 /)
XEXT_COEFF_550_LKT(10,11)=1023.500000 !rg=0.0213956 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(10,12,1:6)=(/ 6011.900000,3289.900000,1355.100000,352.110000,46.294000,4.154700 /)
XPIZA_LKT(10,12,1:6)=(/ 0.982001,0.981830,0.975333,0.980415,0.934167,0.666738 /)
XCGA_LKT(10,12,1:6)=(/ 0.717633,0.675580,0.572067,0.376397,0.122163,0.030443 /)
XEXT_COEFF_550_LKT(10,12)=1385.700000 !rg=0.0213956 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(10,13,1:6)=(/ 5632.600000,3640.800000,1730.100000,524.710000,75.657000,6.335900 /)
XPIZA_LKT(10,13,1:6)=(/ 0.980134,0.982680,0.979252,0.985729,0.957990,0.778028 /)
XCGA_LKT(10,13,1:6)=(/ 0.718387,0.702470,0.628343,0.474027,0.178123,0.044160 /)
XEXT_COEFF_550_LKT(10,13)=1764.300000 !rg=0.0213956 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(10,14,1:6)=(/ 5070.300000,3739.200000,2075.700000,727.570000,118.550000,10.148000 /)
XPIZA_LKT(10,14,1:6)=(/ 0.977476,0.982402,0.981575,0.988916,0.971809,0.858384 /)
XCGA_LKT(10,14,1:6)=(/ 0.717400,0.716000,0.671283,0.545880,0.260457,0.064047 /)
XEXT_COEFF_550_LKT(10,14)=2110.300000 !rg=0.0213956 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(10,15,1:6)=(/ 4492.800000,3581.000000,2327.400000,958.900000,178.410000,16.708000 /)
XPIZA_LKT(10,15,1:6)=(/ 0.974229,0.980946,0.982643,0.990928,0.980017,0.911460 /)
XCGA_LKT(10,15,1:6)=(/ 0.721093,0.718180,0.699983,0.611750,0.372030,0.092970 /)
XEXT_COEFF_550_LKT(10,15)=2359.100000 !rg=0.0213956 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(10,16,1:6)=(/ 3882.700000,3253.600000,2428.800000,1186.200000,264.430000,27.640000 /)
XPIZA_LKT(10,16,1:6)=(/ 0.969720,0.978547,0.982570,0.992154,0.985422,0.944458 /)
XCGA_LKT(10,16,1:6)=(/ 0.723823,0.715653,0.715667,0.659520,0.478127,0.135343 /)
XEXT_COEFF_550_LKT(10,16)=2453.600000 !rg=0.0213956 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(10,17,1:6)=(/ 3271.700000,2906.400000,2353.000000,1374.100000,367.590000,44.529000 /)
XPIZA_LKT(10,17,1:6)=(/ 0.964811,0.975451,0.981346,0.992816,0.988780,0.963958 /)
XCGA_LKT(10,17,1:6)=(/ 0.726757,0.716223,0.719217,0.694283,0.551263,0.197850 /)
XEXT_COEFF_550_LKT(10,17)=2369.400000 !rg=0.0213956 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(10,18,1:6)=(/ 2745.900000,2561.500000,2153.200000,1484.600000,484.900000,68.401000 /)
XPIZA_LKT(10,18,1:6)=(/ 0.959502,0.971690,0.979000,0.993013,0.990877,0.975257 /)
XCGA_LKT(10,18,1:6)=(/ 0.732290,0.720370,0.715430,0.715557,0.619280,0.289407 /)
XEXT_COEFF_550_LKT(10,18)=2160.000000 !rg=0.0213956 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(10,19,1:6)=(/ 2288.100000,2186.000000,1922.500000,1489.800000,603.120000,101.560000 /)
XPIZA_LKT(10,19,1:6)=(/ 0.951589,0.967004,0.975882,0.992741,0.992201,0.982126 /)
XCGA_LKT(10,19,1:6)=(/ 0.737820,0.723347,0.712667,0.723743,0.669027,0.407307 /)
XEXT_COEFF_550_LKT(10,19)=1927.800000 !rg=0.0213956 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(10,20,1:6)=(/ 1881.900000,1851.800000,1702.700000,1396.100000,707.430000,147.990000 /)
XPIZA_LKT(10,20,1:6)=(/ 0.943790,0.961963,0.972587,0.991972,0.992979,0.986744 /)
XCGA_LKT(10,20,1:6)=(/ 0.744250,0.730417,0.718377,0.721283,0.706687,0.505157 /)
XEXT_COEFF_550_LKT(10,20)=1711.500000 !rg=0.0213956 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(11,1,1:6)=(/ 496.440000,99.380000,32.754000,6.934700,2.865900,1.347600 /)
XPIZA_LKT(11,1,1:6)=(/ 0.878489,0.659733,0.357269,0.254054,0.056062,0.007220 /)
XCGA_LKT(11,1,1:6)=(/ 0.081757,0.031723,0.013470,0.005303,0.001680,0.000413 /)
XEXT_COEFF_550_LKT(11,1)=27.416000 !rg=0.0231791 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(11,2,1:6)=(/ 603.570000,116.720000,35.906000,7.403900,2.909400,1.350400 /)
XPIZA_LKT(11,2,1:6)=(/ 0.898114,0.707311,0.410883,0.299984,0.069613,0.009083 /)
XCGA_LKT(11,2,1:6)=(/ 0.103690,0.040233,0.017103,0.006740,0.002133,0.000527 /)
XEXT_COEFF_550_LKT(11,2)=30.718000 !rg=0.0231791 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(11,3,1:6)=(/ 836.770000,158.300000,43.555000,8.546300,3.014900,1.357000 /)
XPIZA_LKT(11,3,1:6)=(/ 0.923863,0.779683,0.509473,0.391067,0.101062,0.013598 /)
XCGA_LKT(11,3,1:6)=(/ 0.156403,0.061030,0.026027,0.010277,0.003257,0.000803 /)
XEXT_COEFF_550_LKT(11,3)=38.740000 !rg=0.0231791 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(11,4,1:6)=(/ 1226.100000,241.210000,59.550000,10.957000,3.237300,1.370800 /)
XPIZA_LKT(11,4,1:6)=(/ 0.945286,0.850684,0.635170,0.521676,0.160960,0.023022 /)
XCGA_LKT(11,4,1:6)=(/ 0.252533,0.102880,0.044250,0.017540,0.005570,0.001373 /)
XEXT_COEFF_550_LKT(11,4)=55.523000 !rg=0.0231791 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(11,5,1:6)=(/ 1802.000000,385.590000,90.377000,15.709000,3.675300,1.397900 /)
XPIZA_LKT(11,5,1:6)=(/ 0.960107,0.902381,0.753477,0.662624,0.258361,0.041149 /)
XCGA_LKT(11,5,1:6)=(/ 0.370760,0.170417,0.073933,0.029423,0.009370,0.002313 /)
XEXT_COEFF_550_LKT(11,5)=87.840000 !rg=0.0231791 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(11,6,1:6)=(/ 2557.100000,605.750000,144.100000,24.390000,4.478600,1.447600 /)
XPIZA_LKT(11,6,1:6)=(/ 0.969598,0.934249,0.840009,0.779151,0.388162,0.072788 /)
XCGA_LKT(11,6,1:6)=(/ 0.475830,0.262247,0.115477,0.046027,0.014697,0.003633 /)
XEXT_COEFF_550_LKT(11,6)=143.980000 !rg=0.0231791 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(11,7,1:6)=(/ 3451.000000,926.240000,231.420000,39.831000,5.924300,1.536800 /)
XPIZA_LKT(11,7,1:6)=(/ 0.975649,0.953830,0.895924,0.861636,0.533865,0.124864 /)
XCGA_LKT(11,7,1:6)=(/ 0.558250,0.366953,0.173323,0.068937,0.022060,0.005463 /)
XEXT_COEFF_550_LKT(11,7)=234.750000 !rg=0.0231791 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(11,8,1:6)=(/ 4385.700000,1349.000000,361.550000,66.264000,8.473400,1.694200 /)
XPIZA_LKT(11,8,1:6)=(/ 0.979419,0.965783,0.929674,0.914249,0.670468,0.203866 /)
XCGA_LKT(11,8,1:6)=(/ 0.620213,0.458870,0.254327,0.101033,0.032357,0.008030 /)
XEXT_COEFF_550_LKT(11,8)=369.340000 !rg=0.0231791 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(11,9,1:6)=(/ 5239.800000,1867.300000,551.200000,110.080000,12.983000,1.974000 /)
XPIZA_LKT(11,9,1:6)=(/ 0.981648,0.973212,0.950519,0.946308,0.781529,0.313836 /)
XCGA_LKT(11,9,1:6)=(/ 0.665487,0.539460,0.358020,0.147220,0.047030,0.011697 /)
XEXT_COEFF_550_LKT(11,9)=565.470000 !rg=0.0231791 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(11,10,1:6)=(/ 5848.200000,2435.600000,815.450000,177.360000,20.833000,2.467800 /)
XPIZA_LKT(11,10,1:6)=(/ 0.982663,0.977888,0.963811,0.965031,0.860868,0.447757 /)
XCGA_LKT(11,10,1:6)=(/ 0.695877,0.603100,0.452790,0.214537,0.068080,0.016970 /)
XEXT_COEFF_550_LKT(11,10)=836.030000 !rg=0.0231791 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(11,11,1:6)=(/ 6087.200000,2999.500000,1137.200000,275.040000,34.395000,3.348200 /)
XPIZA_LKT(11,11,1:6)=(/ 0.982590,0.980788,0.971971,0.976028,0.913232,0.589325 /)
XCGA_LKT(11,11,1:6)=(/ 0.713037,0.652590,0.528933,0.310673,0.098640,0.024607 /)
XEXT_COEFF_550_LKT(11,11)=1164.100000 !rg=0.0231791 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(11,12,1:6)=(/ 5909.600000,3458.800000,1509.300000,417.180000,56.870000,4.910200 /)
XPIZA_LKT(11,12,1:6)=(/ 0.981400,0.982331,0.977181,0.982890,0.945530,0.716377 /)
XCGA_LKT(11,12,1:6)=(/ 0.718887,0.687560,0.596843,0.421040,0.143353,0.035680 /)
XEXT_COEFF_550_LKT(11,12)=1540.800000 !rg=0.0231791 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(11,13,1:6)=(/ 5425.900000,3715.500000,1881.000000,605.060000,91.603000,7.668100 /)
XPIZA_LKT(11,13,1:6)=(/ 0.979195,0.982710,0.980367,0.987247,0.964619,0.815097 /)
XCGA_LKT(11,13,1:6)=(/ 0.718107,0.709037,0.647310,0.503967,0.209360,0.051747 /)
XEXT_COEFF_550_LKT(11,13)=1915.900000 !rg=0.0231791 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(11,14,1:6)=(/ 4850.800000,3707.300000,2195.700000,820.970000,140.800000,12.462000 /)
XPIZA_LKT(11,14,1:6)=(/ 0.976289,0.981945,0.982131,0.989845,0.975676,0.883401 /)
XCGA_LKT(11,14,1:6)=(/ 0.719303,0.717857,0.684367,0.575153,0.304833,0.075053 /)
XEXT_COEFF_550_LKT(11,14)=2229.800000 !rg=0.0231791 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(11,15,1:6)=(/ 4263.200000,3463.200000,2390.800000,1052.200000,210.560000,20.615000 /)
XPIZA_LKT(11,15,1:6)=(/ 0.972535,0.980103,0.982736,0.991501,0.982514,0.927204 /)
XCGA_LKT(11,15,1:6)=(/ 0.722340,0.717623,0.707853,0.632193,0.420600,0.109030 /)
XEXT_COEFF_550_LKT(11,15)=2419.300000 !rg=0.0231791 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(11,16,1:6)=(/ 3621.900000,3118.300000,2418.200000,1269.800000,306.130000,33.883000 /)
XPIZA_LKT(11,16,1:6)=(/ 0.967654,0.977322,0.982205,0.992494,0.987040,0.953888 /)
XCGA_LKT(11,16,1:6)=(/ 0.724440,0.715267,0.718273,0.675397,0.510377,0.159007 /)
XEXT_COEFF_550_LKT(11,16)=2439.000000 !rg=0.0231791 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(11,17,1:6)=(/ 3047.600000,2769.600000,2282.600000,1430.300000,414.540000,53.619000 /)
XPIZA_LKT(11,17,1:6)=(/ 0.963064,0.974098,0.980489,0.992959,0.989751,0.969438 /)
XCGA_LKT(11,17,1:6)=(/ 0.730273,0.718017,0.717997,0.704273,0.580853,0.232873 /)
XEXT_COEFF_550_LKT(11,17)=2293.500000 !rg=0.0231791 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(11,18,1:6)=(/ 2550.200000,2408.500000,2059.000000,1499.300000,534.270000,80.844000 /)
XPIZA_LKT(11,18,1:6)=(/ 0.956089,0.970076,0.977758,0.992955,0.991510,0.978496 /)
XCGA_LKT(11,18,1:6)=(/ 0.734233,0.722637,0.713247,0.719980,0.641493,0.337997 /)
XEXT_COEFF_550_LKT(11,18)=2067.200000 !rg=0.0231791 sigma=2.75 wvl=0.55
 
IF (LHOOK) CALL DR_HOOK('MODE_DUSTOPT:DUST_OPT_LKT_SET1',1,ZHOOK_HANDLE)
END SUBROUTINE DUST_OPT_LKT_SET1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE DUST_OPT_LKT_SET2()

  USE MODD_DUST_OPT_LKT
  
  IMPLICIT NONE


REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_DUSTOPT:DUST_OPT_LKT_SET2',0,ZHOOK_HANDLE)
XEXT_COEFF_WVL_LKT(11,19,1:6)=(/ 2107.700000,2046.700000,1831.600000,1460.700000,649.760000,119.620000 /)
XPIZA_LKT(11,19,1:6)=(/ 0.948817,0.964413,0.974678,0.992480,0.992594,0.984301 /)
XCGA_LKT(11,19,1:6)=(/ 0.739493,0.726350,0.715220,0.723640,0.686390,0.454440 /)
XEXT_COEFF_550_LKT(11,19)=1839.700000 !rg=0.0231791 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(11,20,1:6)=(/ 1729.500000,1722.900000,1611.700000,1337.300000,742.290000,169.760000 /)
XPIZA_LKT(11,20,1:6)=(/ 0.941121,0.959822,0.971116,0.991518,0.993175,0.988125 /)
XCGA_LKT(11,20,1:6)=(/ 0.748793,0.732190,0.721223,0.718110,0.718213,0.534630 /)
XEXT_COEFF_550_LKT(11,20)=1615.400000 !rg=0.0231791 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(12,1,1:6)=(/ 612.620000,117.570000,36.053000,7.425000,2.911400,1.350500 /)
XPIZA_LKT(12,1,1:6)=(/ 0.899073,0.708806,0.412754,0.301717,0.070167,0.009160 /)
XCGA_LKT(12,1,1:6)=(/ 0.095957,0.037187,0.015797,0.006223,0.001970,0.000487 /)
XEXT_COEFF_550_LKT(12,1)=30.868000 !rg=0.0251112 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(12,2,1:6)=(/ 743.530000,139.490000,40.049000,8.020700,2.966500,1.353900 /)
XPIZA_LKT(12,2,1:6)=(/ 0.915121,0.751708,0.468420,0.352124,0.086806,0.011518 /)
XCGA_LKT(12,2,1:6)=(/ 0.121770,0.047150,0.020053,0.007907,0.002503,0.000617 /)
XEXT_COEFF_550_LKT(12,2)=35.058000 !rg=0.0251112 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(12,3,1:6)=(/ 1019.100000,191.700000,49.741000,9.471100,3.100400,1.362300 /)
XPIZA_LKT(12,3,1:6)=(/ 0.935722,0.815094,0.567041,0.448713,0.124954,0.017222 /)
XCGA_LKT(12,3,1:6)=(/ 0.183307,0.071473,0.030503,0.012053,0.003823,0.000940 /)
XEXT_COEFF_550_LKT(12,3)=45.227000 !rg=0.0251112 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(12,4,1:6)=(/ 1465.500000,293.720000,69.939000,12.531000,3.382500,1.379800 /)
XPIZA_LKT(12,4,1:6)=(/ 0.952755,0.874958,0.686112,0.579903,0.195862,0.029081 /)
XCGA_LKT(12,4,1:6)=(/ 0.289610,0.120157,0.051807,0.020560,0.006537,0.001613 /)
XEXT_COEFF_550_LKT(12,4)=66.419000 !rg=0.0251112 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(12,5,1:6)=(/ 2110.800000,466.030000,108.480000,18.556000,3.938400,1.414200 /)
XPIZA_LKT(12,5,1:6)=(/ 0.964722,0.917313,0.791808,0.712626,0.306518,0.051718 /)
XCGA_LKT(12,5,1:6)=(/ 0.409000,0.198073,0.086477,0.034473,0.010990,0.002713 /)
XEXT_COEFF_550_LKT(12,5)=106.800000 !rg=0.0251112 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(12,6,1:6)=(/ 2928.500000,723.860000,174.370000,29.533000,4.957700,1.477200 /)
XPIZA_LKT(12,6,1:6)=(/ 0.972497,0.943372,0.865485,0.816037,0.445700,0.090694 /)
XCGA_LKT(12,6,1:6)=(/ 0.508517,0.299203,0.135057,0.053903,0.017233,0.004263 /)
XEXT_COEFF_550_LKT(12,6)=175.530000 !rg=0.0251112 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(12,7,1:6)=(/ 3855.400000,1088.900000,278.470000,48.929000,6.791300,1.590300 /)
XPIZA_LKT(12,7,1:6)=(/ 0.977461,0.959458,0.911643,0.886041,0.591700,0.153413 /)
XCGA_LKT(12,7,1:6)=(/ 0.583943,0.402423,0.202693,0.080727,0.025860,0.006410 /)
XEXT_COEFF_550_LKT(12,7)=283.480000 !rg=0.0251112 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(12,8,1:6)=(/ 4770.900000,1550.500000,430.620000,81.672000,10.021000,1.790000 /)
XPIZA_LKT(12,8,1:6)=(/ 0.980519,0.969220,0.939321,0.929361,0.719716,0.245311 /)
XCGA_LKT(12,8,1:6)=(/ 0.639583,0.490733,0.294843,0.118380,0.037913,0.009417 /)
XEXT_COEFF_550_LKT(12,8)=440.820000 !rg=0.0251112 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(12,9,1:6)=(/ 5545.500000,2101.400000,651.800000,134.570000,15.718000,2.145100 /)
XPIZA_LKT(12,9,1:6)=(/ 0.982230,0.975421,0.956754,0.955243,0.818062,0.367134 /)
XCGA_LKT(12,9,1:6)=(/ 0.679067,0.566600,0.398997,0.172727,0.055093,0.013717 /)
XEXT_COEFF_550_LKT(12,9)=668.960000 !rg=0.0251112 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(12,10,1:6)=(/ 6014.000000,2675.600000,940.410000,213.190000,25.573000,2.771700 /)
XPIZA_LKT(12,10,1:6)=(/ 0.982789,0.979248,0.967634,0.970222,0.885395,0.506688 /)
XCGA_LKT(12,10,1:6)=(/ 0.704090,0.624120,0.483530,0.251563,0.079760,0.019897 /)
XEXT_COEFF_550_LKT(12,10)=963.120000 !rg=0.0251112 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(12,11,1:6)=(/ 6079.700000,3207.600000,1286.200000,326.900000,42.374000,3.888200 /)
XPIZA_LKT(12,11,1:6)=(/ 0.982277,0.981540,0.974370,0.979195,0.928554,0.644692 /)
XCGA_LKT(12,11,1:6)=(/ 0.716717,0.667873,0.558447,0.357627,0.115643,0.028843 /)
XEXT_COEFF_550_LKT(12,11)=1315.900000 !rg=0.0251112 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(12,12,1:6)=(/ 5754.100000,3596.100000,1661.000000,490.450000,69.517000,5.866400 /)
XPIZA_LKT(12,12,1:6)=(/ 0.980674,0.982603,0.978657,0.984953,0.954649,0.761015 /)
XCGA_LKT(12,12,1:6)=(/ 0.719310,0.697533,0.618027,0.458533,0.168323,0.041813 /)
XEXT_COEFF_550_LKT(12,12)=1694.400000 !rg=0.0251112 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(12,13,1:6)=(/ 5212.300000,3748.900000,2019.200000,689.020000,109.970000,9.347300 /)
XPIZA_LKT(12,13,1:6)=(/ 0.978163,0.982543,0.981250,0.988455,0.969896,0.846903 /)
XCGA_LKT(12,13,1:6)=(/ 0.719030,0.713993,0.663783,0.532790,0.245963,0.060630 /)
XEXT_COEFF_550_LKT(12,13)=2053.700000 !rg=0.0251112 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(12,14,1:6)=(/ 4631.100000,3634.100000,2292.600000,917.550000,166.500000,15.352000 /)
XPIZA_LKT(12,14,1:6)=(/ 0.974928,0.981335,0.982529,0.990632,0.978848,0.904173 /)
XCGA_LKT(12,14,1:6)=(/ 0.720603,0.718560,0.695083,0.601240,0.353387,0.087967 /)
XEXT_COEFF_550_LKT(12,14)=2324.900000 !rg=0.0251112 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(12,15,1:6)=(/ 4003.200000,3336.400000,2423.000000,1145.800000,247.650000,25.400000 /)
XPIZA_LKT(12,15,1:6)=(/ 0.970896,0.979017,0.982654,0.991968,0.984642,0.939983 /)
XCGA_LKT(12,15,1:6)=(/ 0.724270,0.716210,0.713107,0.651170,0.462793,0.127933 /)
XEXT_COEFF_550_LKT(12,15)=2449.000000 !rg=0.0251112 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(12,16,1:6)=(/ 3374.700000,2977.600000,2380.300000,1345.100000,349.220000,41.268000 /)
XPIZA_LKT(12,16,1:6)=(/ 0.966242,0.976120,0.981642,0.992729,0.988326,0.961419 /)
XCGA_LKT(12,16,1:6)=(/ 0.727943,0.716157,0.719087,0.688697,0.539333,0.186983 /)
XEXT_COEFF_550_LKT(12,16)=2396.800000 !rg=0.0251112 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(12,17,1:6)=(/ 2841.400000,2629.700000,2195.600000,1472.200000,464.630000,63.978000 /)
XPIZA_LKT(12,17,1:6)=(/ 0.960131,0.972599,0.979495,0.993003,0.990582,0.973797 /)
XCGA_LKT(12,17,1:6)=(/ 0.731767,0.721037,0.716007,0.712283,0.609023,0.273813 /)
XEXT_COEFF_550_LKT(12,17)=2208.400000 !rg=0.0251112 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(12,18,1:6)=(/ 2354.900000,2260.000000,1964.400000,1495.500000,583.140000,95.315000 /)
XPIZA_LKT(12,18,1:6)=(/ 0.952918,0.967630,0.976582,0.992816,0.992016,0.981193 /)
XCGA_LKT(12,18,1:6)=(/ 0.734833,0.723830,0.713790,0.722797,0.661197,0.389337 /)
XEXT_COEFF_550_LKT(12,18)=1972.600000 !rg=0.0251112 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(12,19,1:6)=(/ 1935.100000,1906.200000,1746.500000,1416.300000,691.510000,139.960000 /)
XPIZA_LKT(12,19,1:6)=(/ 0.946223,0.962787,0.973304,0.992131,0.992882,0.986144 /)
XCGA_LKT(12,19,1:6)=(/ 0.744860,0.728377,0.718063,0.721793,0.701137,0.492890 /)
XEXT_COEFF_550_LKT(12,19)=1750.700000 !rg=0.0251112 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(12,20,1:6)=(/ 1598.000000,1595.400000,1513.700000,1275.600000,769.040000,192.130000 /)
XPIZA_LKT(12,20,1:6)=(/ 0.936499,0.957164,0.969204,0.991039,0.993281,0.989218 /)
XCGA_LKT(12,20,1:6)=(/ 0.751800,0.733320,0.722627,0.716383,0.727550,0.562997 /)
XEXT_COEFF_550_LKT(12,20)=1520.300000 !rg=0.0251112 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(13,1,1:6)=(/ 756.910000,140.600000,40.235000,8.047300,2.969100,1.354100 /)
XPIZA_LKT(13,1,1:6)=(/ 0.916094,0.753061,0.470318,0.353994,0.087483,0.011616 /)
XCGA_LKT(13,1,1:6)=(/ 0.112737,0.043590,0.018527,0.007300,0.002313,0.000570 /)
XEXT_COEFF_550_LKT(13,1)=35.248000 !rg=0.0272044 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(13,2,1:6)=(/ 914.330000,168.270000,45.302000,8.803700,3.039000,1.358500 /)
XPIZA_LKT(13,2,1:6)=(/ 0.929038,0.790940,0.526542,0.407941,0.107742,0.014596 /)
XCGA_LKT(13,2,1:6)=(/ 0.143187,0.055260,0.023513,0.009277,0.002940,0.000723 /)
XEXT_COEFF_550_LKT(13,2)=40.563000 !rg=0.0272044 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(13,3,1:6)=(/ 1233.900000,233.540000,57.577000,10.645000,3.208800,1.369100 /)
XPIZA_LKT(13,3,1:6)=(/ 0.945319,0.845426,0.622506,0.507631,0.153508,0.021789 /)
XCGA_LKT(13,3,1:6)=(/ 0.214797,0.083697,0.035747,0.014137,0.004483,0.001107 /)
XEXT_COEFF_550_LKT(13,3)=53.447000 !rg=0.0272044 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(13,4,1:6)=(/ 1740.200000,357.850000,83.023000,14.529000,3.566700,1.391300 /)
XPIZA_LKT(13,4,1:6)=(/ 0.958843,0.895141,0.732428,0.635792,0.236155,0.036671 /)
XCGA_LKT(13,4,1:6)=(/ 0.329683,0.140227,0.060633,0.024100,0.007667,0.001890 /)
XEXT_COEFF_550_LKT(13,4)=80.141000 !rg=0.0272044 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(13,5,1:6)=(/ 2455.000000,561.640000,130.940000,22.163000,4.272100,1.434900 /)
XPIZA_LKT(13,5,1:6)=(/ 0.968539,0.929604,0.824863,0.757667,0.359217,0.064812 /)
XCGA_LKT(13,5,1:6)=(/ 0.446957,0.229073,0.101100,0.040383,0.012890,0.003187 /)
XEXT_COEFF_550_LKT(13,5)=130.270000 !rg=0.0272044 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(13,6,1:6)=(/ 3323.200000,860.800000,210.910000,36.014000,5.565400,1.514700 /)
XPIZA_LKT(13,6,1:6)=(/ 0.974894,0.950906,0.886657,0.847649,0.504583,0.112457 /)
XCGA_LKT(13,6,1:6)=(/ 0.538943,0.336147,0.157863,0.063113,0.020203,0.005003 /)
XEXT_COEFF_550_LKT(13,6)=213.510000 !rg=0.0272044 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(13,7,1:6)=(/ 4265.500000,1268.000000,333.790000,60.262000,7.890400,1.658200 /)
XPIZA_LKT(13,7,1:6)=(/ 0.978972,0.964040,0.924525,0.906250,0.646889,0.187064 /)
XCGA_LKT(13,7,1:6)=(/ 0.607593,0.436100,0.236393,0.094530,0.030303,0.007520 /)
XEXT_COEFF_550_LKT(13,7)=340.750000 !rg=0.0272044 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(13,8,1:6)=(/ 5140.200000,1770.500000,511.930000,100.440000,11.978000,1.911600 /)
XPIZA_LKT(13,8,1:6)=(/ 0.981406,0.972102,0.947392,0.941593,0.763928,0.292032 /)
XCGA_LKT(13,8,1:6)=(/ 0.656897,0.521970,0.337047,0.138737,0.044417,0.011047 /)
XEXT_COEFF_550_LKT(13,8)=524.970000 !rg=0.0272044 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(13,9,1:6)=(/ 5807.400000,2338.800000,764.130000,163.420000,19.160000,2.362300 /)
XPIZA_LKT(13,9,1:6)=(/ 0.982624,0.977219,0.961901,0.962379,0.849351,0.423783 /)
XCGA_LKT(13,9,1:6)=(/ 0.690823,0.590450,0.434703,0.202643,0.064537,0.016087 /)
XEXT_COEFF_550_LKT(13,9)=783.650000 !rg=0.0272044 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(13,10,1:6)=(/ 6113.900000,2907.000000,1074.000000,254.700000,31.475000,3.157200 /)
XPIZA_LKT(13,10,1:6)=(/ 0.982747,0.980368,0.970733,0.974412,0.905726,0.565266 /)
XCGA_LKT(13,10,1:6)=(/ 0.710673,0.643003,0.513677,0.293607,0.093453,0.023327 /)
XEXT_COEFF_550_LKT(13,10)=1099.600000 !rg=0.0272044 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(13,11,1:6)=(/ 6009.000000,3392.800000,1439.900000,387.940000,52.099000,4.572500 /)
XPIZA_LKT(13,11,1:6)=(/ 0.981784,0.982131,0.976383,0.981869,0.940972,0.696216 /)
XCGA_LKT(13,11,1:6)=(/ 0.718857,0.681017,0.584600,0.402777,0.135643,0.033807 /)
XEXT_COEFF_550_LKT(13,11)=1470.900000 !rg=0.0272044 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(13,12,1:6)=(/ 5568.900000,3692.900000,1814.000000,568.530000,84.395000,7.074900 /)
XPIZA_LKT(13,12,1:6)=(/ 0.979825,0.982727,0.979872,0.986614,0.961930,0.800313 /)
XCGA_LKT(13,12,1:6)=(/ 0.719997,0.705333,0.637783,0.489833,0.197737,0.048993 /)
XEXT_COEFF_550_LKT(13,12)=1848.900000 !rg=0.0272044 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(13,13,1:6)=(/ 4997.500000,3736.800000,2145.800000,779.710000,130.970000,11.457000 /)
XPIZA_LKT(13,13,1:6)=(/ 0.976915,0.982191,0.981886,0.989455,0.974123,0.873779 /)
XCGA_LKT(13,13,1:6)=(/ 0.719497,0.716783,0.677803,0.562360,0.288160,0.071043 /)
XEXT_COEFF_550_LKT(13,13)=2180.700000 !rg=0.0272044 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(13,14,1:6)=(/ 4385.800000,3534.800000,2367.800000,1011.500000,196.640000,18.937000 /)
XPIZA_LKT(13,14,1:6)=(/ 0.973587,0.980506,0.982693,0.991263,0.981524,0.921238 /)
XCGA_LKT(13,14,1:6)=(/ 0.723640,0.717923,0.703847,0.622823,0.402353,0.103137 /)
XEXT_COEFF_550_LKT(13,14)=2398.200000 !rg=0.0272044 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(13,15,1:6)=(/ 3748.500000,3194.900000,2426.600000,1233.400000,288.140000,31.186000 /)
XPIZA_LKT(13,15,1:6)=(/ 0.968862,0.977958,0.982382,0.992351,0.986401,0.950279 /)
XCGA_LKT(13,15,1:6)=(/ 0.725930,0.715997,0.716790,0.668123,0.497097,0.150237 /)
XEXT_COEFF_550_LKT(13,15)=2449.000000 !rg=0.0272044 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(13,16,1:6)=(/ 3149.600000,2843.100000,2317.900000,1408.100000,394.770000,49.848000 /)
XPIZA_LKT(13,16,1:6)=(/ 0.963858,0.974788,0.980888,0.992905,0.989369,0.967407 /)
XCGA_LKT(13,16,1:6)=(/ 0.729117,0.718597,0.718357,0.699697,0.568777,0.220017 /)
XEXT_COEFF_550_LKT(13,16)=2333.300000 !rg=0.0272044 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(13,17,1:6)=(/ 2624.700000,2484.800000,2102.800000,1494.500000,514.430000,75.774000 /)
XPIZA_LKT(13,17,1:6)=(/ 0.957172,0.970685,0.978429,0.992985,0.991270,0.977303 /)
XCGA_LKT(13,17,1:6)=(/ 0.732213,0.722243,0.715193,0.717757,0.632683,0.320510 /)
XEXT_COEFF_550_LKT(13,17)=2111.800000 !rg=0.0272044 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(13,18,1:6)=(/ 2176.000000,2105.100000,1878.800000,1474.000000,631.160000,112.350000 /)
XPIZA_LKT(13,18,1:6)=(/ 0.949932,0.965679,0.975216,0.992584,0.992442,0.983502 /)
XCGA_LKT(13,18,1:6)=(/ 0.740540,0.724647,0.714787,0.723197,0.679400,0.438190 /)
XEXT_COEFF_550_LKT(13,18)=1883.400000 !rg=0.0272044 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(13,19,1:6)=(/ 1785.000000,1774.300000,1653.600000,1362.200000,729.090000,161.400000 /)
XPIZA_LKT(13,19,1:6)=(/ 0.942242,0.960852,0.971715,0.991716,0.993099,0.987641 /)
XCGA_LKT(13,19,1:6)=(/ 0.747343,0.731450,0.720497,0.719703,0.713570,0.523830 /)
XEXT_COEFF_550_LKT(13,19)=1661.600000 !rg=0.0272044 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(13,20,1:6)=(/ 1467.100000,1483.000000,1414.500000,1217.700000,787.460000,216.090000 /)
XPIZA_LKT(13,20,1:6)=(/ 0.931744,0.953610,0.967342,0.990460,0.993319,0.990123 /)
XCGA_LKT(13,20,1:6)=(/ 0.753397,0.736027,0.723803,0.714267,0.734850,0.592470 /)
XEXT_COEFF_550_LKT(13,20)=1416.500000 !rg=0.0272044 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(14,1,1:6)=(/ 934.340000,169.740000,45.538000,8.837200,3.042200,1.358700 /)
XPIZA_LKT(14,1,1:6)=(/ 0.930051,0.792159,0.528412,0.409902,0.108560,0.014720 /)
XCGA_LKT(14,1,1:6)=(/ 0.132650,0.051093,0.021723,0.008567,0.002713,0.000667 /)
XEXT_COEFF_550_LKT(14,1)=40.805000 !rg=0.0294721 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(14,2,1:6)=(/ 1119.700000,204.560000,51.960000,9.797700,3.130800,1.364200 /)
XPIZA_LKT(14,2,1:6)=(/ 0.940333,0.824980,0.583643,0.466110,0.132970,0.018480 /)
XCGA_LKT(14,2,1:6)=(/ 0.168667,0.064767,0.027567,0.010880,0.003450,0.000850 /)
XEXT_COEFF_550_LKT(14,2)=47.546000 !rg=0.0294721 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(14,3,1:6)=(/ 1483.300000,285.630000,67.493000,12.136000,3.346300,1.377600 /)
XPIZA_LKT(14,3,1:6)=(/ 0.953057,0.871037,0.674551,0.566175,0.187158,0.027533 /)
XCGA_LKT(14,3,1:6)=(/ 0.251383,0.098013,0.041883,0.016577,0.005263,0.001297 /)
XEXT_COEFF_550_LKT(14,3)=63.852000 !rg=0.0294721 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(14,4,1:6)=(/ 2052.400000,435.470000,99.454000,17.064000,3.800400,1.405800 /)
XPIZA_LKT(14,4,1:6)=(/ 0.963822,0.911787,0.773610,0.688029,0.281775,0.046147 /)
XCGA_LKT(14,4,1:6)=(/ 0.371960,0.163453,0.070937,0.028243,0.008993,0.002220 /)
XEXT_COEFF_550_LKT(14,4)=97.364000 !rg=0.0294721 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(14,5,1:6)=(/ 2832.100000,674.360000,158.580000,26.725000,4.695600,1.461000 /)
XPIZA_LKT(14,5,1:6)=(/ 0.971696,0.939714,0.852912,0.797372,0.415433,0.080931 /)
XCGA_LKT(14,5,1:6)=(/ 0.483567,0.262757,0.118127,0.047290,0.015113,0.003737 /)
XEXT_COEFF_550_LKT(14,5)=159.120000 !rg=0.0294721 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(14,6,1:6)=(/ 3739.000000,1016.300000,254.580000,44.149000,6.336200,1.562200 /)
XPIZA_LKT(14,6,1:6)=(/ 0.976892,0.957079,0.904115,0.874329,0.563164,0.138627 /)
XCGA_LKT(14,6,1:6)=(/ 0.567457,0.372157,0.184297,0.073883,0.023680,0.005867 /)
XEXT_COEFF_550_LKT(14,6)=258.830000 !rg=0.0294721 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(14,7,1:6)=(/ 4670.800000,1465.200000,398.840000,74.268000,9.282300,1.744300 /)
XPIZA_LKT(14,7,1:6)=(/ 0.980199,0.967808,0.935150,0.922809,0.698180,0.226082 /)
XCGA_LKT(14,7,1:6)=(/ 0.628890,0.469367,0.273953,0.110690,0.035510,0.008820 /)
XEXT_COEFF_550_LKT(14,7)=408.100000 !rg=0.0294721 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(14,8,1:6)=(/ 5475.300000,2002.900000,606.650000,123.000000,14.449000,2.065900 /)
XPIZA_LKT(14,8,1:6)=(/ 0.982091,0.974509,0.954162,0.951427,0.802782,0.343517 /)
XCGA_LKT(14,8,1:6)=(/ 0.672023,0.550467,0.377430,0.162613,0.052033,0.012957 /)
XEXT_COEFF_550_LKT(14,8)=622.540000 !rg=0.0294721 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(14,9,1:6)=(/ 6009.800000,2581.800000,885.130000,197.020000,23.476000,2.637800 /)
XPIZA_LKT(14,9,1:6)=(/ 0.982848,0.978706,0.966062,0.968084,0.875745,0.482375 /)
XCGA_LKT(14,9,1:6)=(/ 0.700390,0.612620,0.466537,0.237453,0.075600,0.018863 /)
XEXT_COEFF_550_LKT(14,9)=906.870000 !rg=0.0294721 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(14,10,1:6)=(/ 6146.600000,3124.800000,1219.200000,303.280000,38.768000,3.646000 /)
XPIZA_LKT(14,10,1:6)=(/ 0.982529,0.981225,0.973332,0.977864,0.922403,0.621882 /)
XCGA_LKT(14,10,1:6)=(/ 0.715483,0.659437,0.543810,0.338967,0.109530,0.027343 /)
XEXT_COEFF_550_LKT(14,10)=1247.900000 !rg=0.0294721 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(14,11,1:6)=(/ 5883.100000,3547.000000,1592.300000,457.600000,63.791000,5.438900 /)
XPIZA_LKT(14,11,1:6)=(/ 0.981167,0.982500,0.978007,0.984106,0.950964,0.742997 /)
XCGA_LKT(14,11,1:6)=(/ 0.720467,0.692030,0.606900,0.441873,0.159183,0.039617 /)
XEXT_COEFF_550_LKT(14,11)=1624.900000 !rg=0.0294721 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(14,12,1:6)=(/ 5367.700000,3749.700000,1958.300000,650.220000,101.630000,8.599200 /)
XPIZA_LKT(14,12,1:6)=(/ 0.978775,0.982643,0.980866,0.987937,0.967730,0.834268 /)
XCGA_LKT(14,12,1:6)=(/ 0.719963,0.711213,0.655380,0.518980,0.232240,0.057407 /)
XEXT_COEFF_550_LKT(14,12)=1992.800000 !rg=0.0294721 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(14,13,1:6)=(/ 4765.100000,3685.800000,2253.900000,875.570000,155.140000,14.097000 /)
XPIZA_LKT(14,13,1:6)=(/ 0.975749,0.981656,0.982372,0.990306,0.977566,0.896197 /)
XCGA_LKT(14,13,1:6)=(/ 0.722030,0.718273,0.689603,0.589807,0.335070,0.083257 /)
XEXT_COEFF_550_LKT(14,13)=2286.700000 !rg=0.0294721 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(14,14,1:6)=(/ 4139.100000,3409.000000,2413.500000,1105.400000,231.770000,23.344000 /)
XPIZA_LKT(14,14,1:6)=(/ 0.971821,0.979622,0.982711,0.991769,0.983808,0.935134 /)
XCGA_LKT(14,14,1:6)=(/ 0.725013,0.717487,0.710220,0.642383,0.446533,0.120983 /)
XEXT_COEFF_550_LKT(14,14)=2441.100000 !rg=0.0294721 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(14,15,1:6)=(/ 3490.500000,3059.400000,2402.200000,1312.000000,330.330000,38.073000 /)
XPIZA_LKT(14,15,1:6)=(/ 0.967075,0.976809,0.981901,0.992623,0.987808,0.958521 /)
XCGA_LKT(14,15,1:6)=(/ 0.726613,0.717183,0.718253,0.682340,0.526760,0.176580 /)
XEXT_COEFF_550_LKT(14,15)=2421.400000 !rg=0.0294721 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(14,16,1:6)=(/ 2925.900000,2706.600000,2239.800000,1456.000000,443.900000,59.668000 /)
XPIZA_LKT(14,16,1:6)=(/ 0.961426,0.973261,0.979980,0.992983,0.990256,0.972166 /)
XCGA_LKT(14,16,1:6)=(/ 0.730730,0.720127,0.717323,0.708473,0.597887,0.258763 /)
XEXT_COEFF_550_LKT(14,16)=2249.600000 !rg=0.0294721 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(14,17,1:6)=(/ 2433.100000,2324.500000,2013.100000,1499.800000,563.230000,89.408000 /)
XPIZA_LKT(14,17,1:6)=(/ 0.953627,0.968810,0.977145,0.992873,0.991819,0.980192 /)
XCGA_LKT(14,17,1:6)=(/ 0.735777,0.723520,0.714057,0.721393,0.653023,0.371067 /)
XEXT_COEFF_550_LKT(14,17)=2019.400000 !rg=0.0294721 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(14,18,1:6)=(/ 1999.000000,1962.300000,1788.100000,1436.000000,674.460000,131.870000 /)
XPIZA_LKT(14,18,1:6)=(/ 0.946939,0.964339,0.973947,0.992279,0.992769,0.985475 /)
XCGA_LKT(14,18,1:6)=(/ 0.742833,0.728987,0.717470,0.722443,0.695037,0.479370 /)
XEXT_COEFF_550_LKT(14,18)=1798.000000 !rg=0.0294721 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(14,19,1:6)=(/ 1642.200000,1651.100000,1557.000000,1304.200000,758.530000,183.370000 /)
XPIZA_LKT(14,19,1:6)=(/ 0.937739,0.957817,0.970088,0.991203,0.993241,0.988827 /)
XCGA_LKT(14,19,1:6)=(/ 0.749310,0.733090,0.722550,0.716437,0.723580,0.552193 /)
XEXT_COEFF_550_LKT(14,19)=1560.100000 !rg=0.0294721 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(14,20,1:6)=(/ 1350.000000,1370.900000,1326.100000,1159.100000,795.930000,241.410000 /)
XPIZA_LKT(14,20,1:6)=(/ 0.926849,0.950503,0.964539,0.989933,0.993274,0.990899 /)
XCGA_LKT(14,20,1:6)=(/ 0.759170,0.737157,0.725363,0.714997,0.739617,0.619943 /)
XEXT_COEFF_550_LKT(14,20)=1321.300000 !rg=0.0294721 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(15,1,1:6)=(/ 1149.500000,206.560000,52.261000,9.840200,3.134900,1.364500 /)
XPIZA_LKT(15,1,1:6)=(/ 0.941408,0.826090,0.585441,0.468112,0.133950,0.018637 /)
XCGA_LKT(15,1,1:6)=(/ 0.156420,0.059897,0.025470,0.010050,0.003183,0.000783 /)
XEXT_COEFF_550_LKT(15,1)=47.856000 !rg=0.0319288 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(15,2,1:6)=(/ 1361.900000,250.200000,60.400000,11.060000,3.247300,1.371500 /)
XPIZA_LKT(15,2,1:6)=(/ 0.949427,0.854060,0.638245,0.525068,0.163002,0.023373 /)
XCGA_LKT(15,2,1:6)=(/ 0.199117,0.075927,0.032313,0.012763,0.004047,0.000997 /)
XEXT_COEFF_550_LKT(15,2)=56.401000 !rg=0.0319288 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(15,3,1:6)=(/ 1769.300000,350.040000,80.029000,14.028000,3.520700,1.388500 /)
XPIZA_LKT(15,3,1:6)=(/ 0.959290,0.892400,0.722194,0.622745,0.226171,0.034734 /)
XCGA_LKT(15,3,1:6)=(/ 0.293287,0.114793,0.049063,0.019437,0.006173,0.001523 /)
XEXT_COEFF_550_LKT(15,3)=77.008000 !rg=0.0319288 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(15,4,1:6)=(/ 2403.400000,528.460000,120.000000,20.277000,4.097000,1.424100 /)
XPIZA_LKT(15,4,1:6)=(/ 0.967908,0.925437,0.809506,0.735646,0.332291,0.057918 /)
XCGA_LKT(15,4,1:6)=(/ 0.415220,0.190123,0.082957,0.033090,0.010547,0.002603 /)
XEXT_COEFF_550_LKT(15,4)=118.880000 !rg=0.0319288 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(15,5,1:6)=(/ 3239.500000,805.590000,192.290000,32.483000,5.233000,1.494200 /)
XPIZA_LKT(15,5,1:6)=(/ 0.974312,0.948006,0.876404,0.831713,0.473804,0.100615 /)
XCGA_LKT(15,5,1:6)=(/ 0.518313,0.298230,0.137907,0.055363,0.017720,0.004387 /)
XEXT_COEFF_550_LKT(15,5)=194.260000 !rg=0.0319288 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(15,6,1:6)=(/ 4166.200000,1190.200000,306.370000,54.311000,7.313400,1.622600 /)
XPIZA_LKT(15,6,1:6)=(/ 0.978555,0.962126,0.918469,0.896550,0.619831,0.169698 /)
XCGA_LKT(15,6,1:6)=(/ 0.593730,0.407553,0.214610,0.086473,0.027753,0.006883 /)
XEXT_COEFF_550_LKT(15,6)=312.520000 !rg=0.0319288 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(15,7,1:6)=(/ 5064.900000,1681.700000,475.450000,91.408000,11.044000,1.853700 /)
XPIZA_LKT(15,7,1:6)=(/ 0.981200,0.970958,0.943997,0.936261,0.744695,0.270477 /)
XCGA_LKT(15,7,1:6)=(/ 0.648070,0.501913,0.313657,0.129610,0.041600,0.010347 /)
XEXT_COEFF_550_LKT(15,7)=487.390000 !rg=0.0319288 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(15,8,1:6)=(/ 5771.400000,2241.700000,713.430000,149.760000,17.561000,2.261700 /)
XPIZA_LKT(15,8,1:6)=(/ 0.982585,0.976481,0.959752,0.959302,0.836298,0.398863 /)
XCGA_LKT(15,8,1:6)=(/ 0.685290,0.575887,0.413743,0.190560,0.060943,0.015193 /)
XEXT_COEFF_550_LKT(15,8)=731.800000 !rg=0.0319288 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(15,9,1:6)=(/ 6149.000000,2820.700000,1015.200000,236.070000,28.859000,2.987300 /)
XPIZA_LKT(15,9,1:6)=(/ 0.982903,0.979940,0.969439,0.972682,0.897722,0.541294 /)
XCGA_LKT(15,9,1:6)=(/ 0.708313,0.632777,0.497383,0.277153,0.088560,0.022113 /)
XEXT_COEFF_550_LKT(15,9)=1039.800000 !rg=0.0319288 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(15,10,1:6)=(/ 6110.400000,3324.900000,1371.600000,360.490000,47.691000,4.265600 /)
XPIZA_LKT(15,10,1:6)=(/ 0.982158,0.981905,0.975521,0.980762,0.935962,0.675141 /)
XCGA_LKT(15,10,1:6)=(/ 0.718917,0.673797,0.571207,0.383847,0.128410,0.032047 /)
XEXT_COEFF_550_LKT(15,10)=1402.000000 !rg=0.0319288 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(15,11,1:6)=(/ 5722.600000,3666.000000,1746.500000,533.020000,77.629000,6.534500 /)
XPIZA_LKT(15,11,1:6)=(/ 0.980332,0.982711,0.979335,0.985922,0.958958,0.784535 /)
XCGA_LKT(15,11,1:6)=(/ 0.720930,0.701087,0.627480,0.474643,0.186887,0.046420 /)
XEXT_COEFF_550_LKT(15,11)=1781.200000 !rg=0.0319288 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(15,12,1:6)=(/ 5144.400000,3760.400000,2090.800000,738.080000,121.390000,10.517000 /)
XPIZA_LKT(15,12,1:6)=(/ 0.977753,0.982399,0.981598,0.989020,0.972368,0.863133 /)
XCGA_LKT(15,12,1:6)=(/ 0.721503,0.715230,0.670400,0.548707,0.272203,0.067260 /)
XEXT_COEFF_550_LKT(15,12)=2126.100000 !rg=0.0319288 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(15,13,1:6)=(/ 4531.700000,3598.500000,2339.700000,970.210000,183.400000,17.379000 /)
XPIZA_LKT(15,13,1:6)=(/ 0.974390,0.980982,0.982626,0.990999,0.980446,0.914686 /)
XCGA_LKT(15,13,1:6)=(/ 0.723933,0.718707,0.699257,0.612740,0.383740,0.097593 /)
XEXT_COEFF_550_LKT(15,13)=2371.500000 !rg=0.0319288 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(15,14,1:6)=(/ 3870.400000,3279.500000,2432.300000,1196.400000,270.850000,28.699000 /)
XPIZA_LKT(15,14,1:6)=(/ 0.970001,0.978556,0.982514,0.992192,0.985713,0.946360 /)
XCGA_LKT(15,14,1:6)=(/ 0.725747,0.717317,0.714823,0.660340,0.483033,0.142013 /)
XEXT_COEFF_550_LKT(15,14)=2455.900000 !rg=0.0319288 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(15,15,1:6)=(/ 3252.100000,2924.000000,2351.900000,1381.500000,374.550000,46.127000 /)
XPIZA_LKT(15,15,1:6)=(/ 0.964953,0.975493,0.981254,0.992828,0.988938,0.965084 /)
XCGA_LKT(15,15,1:6)=(/ 0.728560,0.718153,0.718627,0.694350,0.556067,0.207687 /)
XEXT_COEFF_550_LKT(15,15)=2366.800000 !rg=0.0319288 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(15,16,1:6)=(/ 2718.900000,2554.100000,2152.900000,1487.000000,493.930000,70.845000 /)
XPIZA_LKT(15,16,1:6)=(/ 0.957989,0.971768,0.978868,0.992998,0.991002,0.975977 /)
XCGA_LKT(15,16,1:6)=(/ 0.732930,0.723013,0.715343,0.715117,0.623103,0.303373 /)
XEXT_COEFF_550_LKT(15,16)=2161.900000 !rg=0.0319288 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(15,17,1:6)=(/ 2239.400000,2171.800000,1920.500000,1484.900000,612.090000,105.440000 /)
XPIZA_LKT(15,17,1:6)=(/ 0.951196,0.966862,0.975971,0.992686,0.992276,0.982648 /)
XCGA_LKT(15,17,1:6)=(/ 0.738167,0.726167,0.715417,0.722743,0.671953,0.421003 /)
XEXT_COEFF_550_LKT(15,17)=1931.000000 !rg=0.0319288 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(15,18,1:6)=(/ 1837.100000,1833.300000,1699.200000,1386.800000,714.140000,152.880000 /)
XPIZA_LKT(15,18,1:6)=(/ 0.943661,0.961747,0.972534,0.991885,0.993010,0.987098 /)
XCGA_LKT(15,18,1:6)=(/ 0.745410,0.730577,0.720457,0.720057,0.708350,0.512257 /)
XEXT_COEFF_550_LKT(15,18)=1702.600000 !rg=0.0319288 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(15,19,1:6)=(/ 1513.400000,1525.300000,1460.300000,1242.700000,780.800000,206.640000 /)
XPIZA_LKT(15,19,1:6)=(/ 0.933103,0.954611,0.967895,0.990715,0.993301,0.989791 /)
XCGA_LKT(15,19,1:6)=(/ 0.754093,0.733407,0.723317,0.715443,0.731823,0.581517 /)
XEXT_COEFF_550_LKT(15,19)=1459.000000 !rg=0.0319288 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(15,20,1:6)=(/ 1234.100000,1264.200000,1237.100000,1106.400000,794.150000,266.150000 /)
XPIZA_LKT(15,20,1:6)=(/ 0.922375,0.947587,0.962192,0.989360,0.993154,0.991536 /)
XCGA_LKT(15,20,1:6)=(/ 0.762373,0.742967,0.726490,0.717157,0.742547,0.642613 /)
XEXT_COEFF_550_LKT(15,20)=1239.300000 !rg=0.0319288 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(16,1,1:6)=(/ 1405.500000,252.970000,60.785000,11.114000,3.252500,1.371900 /)
XPIZA_LKT(16,1,1:6)=(/ 0.950568,0.855086,0.639931,0.527052,0.164158,0.023570 /)
XCGA_LKT(16,1,1:6)=(/ 0.184997,0.070227,0.029860,0.011787,0.003737,0.000920 /)
XEXT_COEFF_550_LKT(16,1)=56.800000 !rg=0.0345903 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(16,2,1:6)=(/ 1641.400000,307.340000,71.091000,12.662000,3.395100,1.380700 /)
XPIZA_LKT(16,2,1:6)=(/ 0.956687,0.878568,0.689125,0.583166,0.198230,0.029520 /)
XCGA_LKT(16,2,1:6)=(/ 0.235637,0.089037,0.037873,0.014967,0.004747,0.001170 /)
XEXT_COEFF_550_LKT(16,2)=67.623000 !rg=0.0345903 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(16,3,1:6)=(/ 2094.000000,428.920000,95.848000,16.431000,3.742000,1.402200 /)
XPIZA_LKT(16,3,1:6)=(/ 0.964319,0.910038,0.764834,0.675955,0.270550,0.043732 /)
XCGA_LKT(16,3,1:6)=(/ 0.340120,0.134473,0.057470,0.022783,0.007243,0.001787 /)
XEXT_COEFF_550_LKT(16,3)=93.607000 !rg=0.0345903 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(16,4,1:6)=(/ 2793.400000,638.680000,145.550000,24.348000,4.473400,1.447400 /)
XPIZA_LKT(16,4,1:6)=(/ 0.971274,0.936588,0.840263,0.778085,0.386886,0.072458 /)
XCGA_LKT(16,4,1:6)=(/ 0.458097,0.220377,0.096963,0.038760,0.012370,0.003057 /)
XEXT_COEFF_550_LKT(16,4)=145.620000 !rg=0.0345903 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(16,5,1:6)=(/ 3673.300000,956.250000,233.050000,39.732000,5.914600,1.536300 /)
XPIZA_LKT(16,5,1:6)=(/ 0.976492,0.954789,0.895895,0.860930,0.532738,0.124424 /)
XCGA_LKT(16,5,1:6)=(/ 0.550893,0.334907,0.160813,0.064797,0.020773,0.005143 /)
XEXT_COEFF_550_LKT(16,5)=236.660000 !rg=0.0345903 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(16,6,1:6)=(/ 4596.100000,1384.000000,367.620000,66.923000,8.551600,1.699200 /)
XPIZA_LKT(16,6,1:6)=(/ 0.979922,0.966286,0.930303,0.914854,0.673179,0.206030 /)
XCGA_LKT(16,6,1:6)=(/ 0.617603,0.442847,0.248573,0.101187,0.032520,0.008077 /)
XEXT_COEFF_550_LKT(16,6)=376.000000 !rg=0.0345903 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(16,7,1:6)=(/ 5431.700000,1912.900000,564.800000,112.140000,13.269000,1.992400 /)
XPIZA_LKT(16,7,1:6)=(/ 0.981992,0.973582,0.951373,0.947114,0.785950,0.319912 /)
XCGA_LKT(16,7,1:6)=(/ 0.665050,0.532057,0.352730,0.151747,0.048730,0.012133 /)
XEXT_COEFF_550_LKT(16,7)=579.490000 !rg=0.0345903 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(16,8,1:6)=(/ 6012.700000,2488.200000,830.010000,181.140000,21.468000,2.510200 /)
XPIZA_LKT(16,8,1:6)=(/ 0.982909,0.978115,0.964291,0.965606,0.864747,0.456795 /)
XCGA_LKT(16,8,1:6)=(/ 0.696410,0.599527,0.446717,0.223043,0.071377,0.017817 /)
XEXT_COEFF_550_LKT(16,8)=850.770000 !rg=0.0345903 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(16,9,1:6)=(/ 6221.200000,3048.600000,1157.000000,281.750000,35.529000,3.430600 /)
XPIZA_LKT(16,9,1:6)=(/ 0.982783,0.980908,0.972258,0.976451,0.915821,0.598896 /)
XCGA_LKT(16,9,1:6)=(/ 0.714380,0.650500,0.528187,0.320517,0.103757,0.025923 /)
XEXT_COEFF_550_LKT(16,9)=1184.800000 !rg=0.0345903 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(16,10,1:6)=(/ 6019.900000,3496.300000,1524.500000,426.400000,58.470000,5.050200 /)
XPIZA_LKT(16,10,1:6)=(/ 0.981603,0.982376,0.977303,0.983187,0.946902,0.723983 /)
XCGA_LKT(16,10,1:6)=(/ 0.721193,0.685977,0.594760,0.423960,0.150607,0.037553 /)
XEXT_COEFF_550_LKT(16,10)=1556.600000 !rg=0.0345903 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(16,11,1:6)=(/ 5522.700000,3744.700000,1896.000000,612.450000,93.752000,7.917400 /)
XPIZA_LKT(16,11,1:6)=(/ 0.979476,0.982726,0.980438,0.987370,0.965336,0.820695 /)
XCGA_LKT(16,11,1:6)=(/ 0.721983,0.708000,0.646190,0.504423,0.219387,0.054390 /)
XEXT_COEFF_550_LKT(16,11)=1930.700000 !rg=0.0345903 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(16,12,1:6)=(/ 4922.700000,3730.500000,2209.400000,832.350000,144.110000,12.919000 /)
XPIZA_LKT(16,12,1:6)=(/ 0.976538,0.981975,0.982165,0.989940,0.976126,0.887331 /)
XCGA_LKT(16,12,1:6)=(/ 0.723173,0.717697,0.683297,0.577240,0.317160,0.078813 /)
XEXT_COEFF_550_LKT(16,12)=2243.000000 !rg=0.0345903 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(16,13,1:6)=(/ 4271.600000,3491.000000,2399.700000,1064.200000,216.500000,21.427000 /)
XPIZA_LKT(16,13,1:6)=(/ 0.972789,0.980105,0.982729,0.991552,0.982897,0.929792 /)
XCGA_LKT(16,13,1:6)=(/ 0.725450,0.718520,0.706837,0.632970,0.429297,0.114447 /)
XEXT_COEFF_550_LKT(16,13)=2428.300000 !rg=0.0345903 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(16,14,1:6)=(/ 3608.200000,3143.000000,2420.300000,1278.400000,312.090000,35.109000 /)
XPIZA_LKT(16,14,1:6)=(/ 0.968267,0.977394,0.982144,0.992507,0.987250,0.955368 /)
XCGA_LKT(16,14,1:6)=(/ 0.727190,0.717277,0.717317,0.675550,0.513813,0.166830 /)
XEXT_COEFF_550_LKT(16,14)=2441.300000 !rg=0.0345903 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(16,15,1:6)=(/ 3036.200000,2782.200000,2282.600000,1435.500000,422.340000,55.392000 /)
XPIZA_LKT(16,15,1:6)=(/ 0.962064,0.974177,0.980413,0.992948,0.989888,0.970302 /)
XCGA_LKT(16,15,1:6)=(/ 0.730563,0.721227,0.717567,0.703973,0.585737,0.244260 /)
XEXT_COEFF_550_LKT(16,15)=2296.300000 !rg=0.0345903 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(16,16,1:6)=(/ 2503.400000,2402.300000,2059.200000,1500.000000,542.800000,83.695000 /)
XPIZA_LKT(16,16,1:6)=(/ 0.955322,0.969600,0.977795,0.992920,0.991601,0.979091 /)
XCGA_LKT(16,16,1:6)=(/ 0.733420,0.724247,0.715290,0.719390,0.644320,0.352637 /)
XEXT_COEFF_550_LKT(16,16)=2069.400000 !rg=0.0345903 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(16,17,1:6)=(/ 2057.300000,2024.100000,1834.500000,1454.400000,657.150000,124.050000 /)
XPIZA_LKT(16,17,1:6)=(/ 0.948660,0.964982,0.974666,0.992412,0.992644,0.984753 /)
XCGA_LKT(16,17,1:6)=(/ 0.742273,0.727320,0.717790,0.722497,0.688540,0.464737 /)
XEXT_COEFF_550_LKT(16,17)=1840.500000 !rg=0.0345903 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(16,18,1:6)=(/ 1698.700000,1700.300000,1604.300000,1329.300000,746.620000,174.490000 /)
XPIZA_LKT(16,18,1:6)=(/ 0.938939,0.958970,0.970706,0.991458,0.993187,0.988388 /)
XCGA_LKT(16,18,1:6)=(/ 0.749343,0.731637,0.721833,0.718357,0.719183,0.541017 /)
XEXT_COEFF_550_LKT(16,18)=1607.100000 !rg=0.0345903 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(16,19,1:6)=(/ 1385.600000,1415.400000,1361.800000,1187.400000,792.800000,231.610000 /)
XPIZA_LKT(16,19,1:6)=(/ 0.928875,0.951163,0.965611,0.990159,0.993291,0.990615 /)
XCGA_LKT(16,19,1:6)=(/ 0.756060,0.738097,0.724123,0.715270,0.737330,0.610087 /)
XEXT_COEFF_550_LKT(16,19)=1362.700000 !rg=0.0345903 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(16,20,1:6)=(/ 1128.000000,1162.500000,1152.900000,1052.900000,782.860000,290.350000 /)
XPIZA_LKT(16,20,1:6)=(/ 0.918338,0.944377,0.960402,0.988680,0.992950,0.992041 /)
XCGA_LKT(16,20,1:6)=(/ 0.768093,0.745363,0.731920,0.718663,0.743463,0.662273 /)
XEXT_COEFF_550_LKT(16,20)=1153.500000 !rg=0.0345903 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(17,1,1:6)=(/ 1702.300000,311.270000,71.589000,12.730000,3.401600,1.381200 /)
XPIZA_LKT(17,1,1:6)=(/ 0.957880,0.879544,0.690674,0.585075,0.199570,0.029767 /)
XCGA_LKT(17,1,1:6)=(/ 0.219647,0.082373,0.035007,0.013827,0.004383,0.001080 /)
XEXT_COEFF_550_LKT(17,1)=68.140000 !rg=0.0374736 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(17,2,1:6)=(/ 1956.900000,378.440000,84.625000,14.697000,3.582600,1.392400 /)
XPIZA_LKT(17,2,1:6)=(/ 0.962434,0.898990,0.735415,0.638854,0.238850,0.037220 /)
XCGA_LKT(17,2,1:6)=(/ 0.279393,0.104473,0.044390,0.017553,0.005570,0.001373 /)
XEXT_COEFF_550_LKT(17,2)=81.832000 !rg=0.0374736 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(17,3,1:6)=(/ 2460.400000,524.430000,115.760000,19.479000,4.022900,1.419600 /)
XPIZA_LKT(17,3,1:6)=(/ 0.968397,0.924478,0.802234,0.724756,0.319960,0.054923 /)
XCGA_LKT(17,3,1:6)=(/ 0.390587,0.157577,0.067303,0.026707,0.008497,0.002097 /)
XEXT_COEFF_550_LKT(17,3)=114.500000 !rg=0.0374736 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(17,4,1:6)=(/ 3221.000000,767.810000,177.110000,29.499000,4.951000,1.476900 /)
XPIZA_LKT(17,4,1:6)=(/ 0.974064,0.945673,0.866235,0.815149,0.444356,0.090287 /)
XCGA_LKT(17,4,1:6)=(/ 0.499283,0.254203,0.113273,0.045387,0.014507,0.003587 /)
XEXT_COEFF_550_LKT(17,4)=178.590000 !rg=0.0374736 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(17,5,1:6)=(/ 4125.200000,1127.400000,281.940000,48.821000,6.779000,1.589700 /)
XPIZA_LKT(17,5,1:6)=(/ 0.978307,0.960344,0.911972,0.885436,0.590581,0.152885 /)
XCGA_LKT(17,5,1:6)=(/ 0.580860,0.372550,0.187140,0.075807,0.024347,0.006037 /)
XEXT_COEFF_550_LKT(17,5)=287.440000 !rg=0.0374736 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(17,6,1:6)=(/ 5019.600000,1598.200000,439.870000,82.442000,10.119000,1.796300 /)
XPIZA_LKT(17,6,1:6)=(/ 0.981049,0.969747,0.940102,0.929791,0.722137,0.247769 /)
XCGA_LKT(17,6,1:6)=(/ 0.639180,0.477420,0.285053,0.118363,0.038100,0.009473 /)
XEXT_COEFF_550_LKT(17,6)=450.780000 !rg=0.0374736 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(17,7,1:6)=(/ 5763.600000,2154.300000,666.270000,136.900000,16.074000,2.168400 /)
XPIZA_LKT(17,7,1:6)=(/ 0.982590,0.975742,0.957452,0.955830,0.821821,0.373660 /)
XCGA_LKT(17,7,1:6)=(/ 0.680063,0.559497,0.389223,0.177577,0.057070,0.014230 /)
XEXT_COEFF_550_LKT(17,7)=683.560000 !rg=0.0374736 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(17,8,1:6)=(/ 6193.700000,2734.400000,956.550000,217.760000,26.350000,2.825400 /)
XPIZA_LKT(17,8,1:6)=(/ 0.983065,0.979477,0.967986,0.970686,0.888559,0.515770 /)
XCGA_LKT(17,8,1:6)=(/ 0.705827,0.621177,0.478603,0.260170,0.083593,0.020887 /)
XEXT_COEFF_550_LKT(17,8)=980.190000 !rg=0.0374736 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(17,9,1:6)=(/ 6222.900000,3262.800000,1307.800000,335.540000,43.717000,3.992500 /)
XPIZA_LKT(17,9,1:6)=(/ 0.982514,0.981680,0.974631,0.979596,0.930586,0.653678 /)
XCGA_LKT(17,9,1:6)=(/ 0.719093,0.666173,0.556733,0.364420,0.121590,0.030383 /)
XEXT_COEFF_550_LKT(17,9)=1337.500000 !rg=0.0374736 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(17,10,1:6)=(/ 5875.500000,3636.700000,1679.800000,498.820000,71.302000,6.042700 /)
XPIZA_LKT(17,10,1:6)=(/ 0.980924,0.982674,0.978757,0.985166,0.955676,0.767746 /)
XCGA_LKT(17,10,1:6)=(/ 0.722637,0.696307,0.616300,0.458100,0.176700,0.044003 /)
XEXT_COEFF_550_LKT(17,10)=1714.100000 !rg=0.0374736 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(17,11,1:6)=(/ 5313.500000,3780.100000,2034.500000,697.620000,112.320000,9.658600 /)
XPIZA_LKT(17,11,1:6)=(/ 0.978456,0.982578,0.981276,0.988548,0.970431,0.851630 /)
XCGA_LKT(17,11,1:6)=(/ 0.723293,0.713240,0.662290,0.534357,0.257137,0.063723 /)
XEXT_COEFF_550_LKT(17,11)=2069.900000 !rg=0.0374736 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(17,12,1:6)=(/ 4679.100000,3665.000000,2305.900000,927.300000,170.580000,15.914000 /)
XPIZA_LKT(17,12,1:6)=(/ 0.975157,0.981374,0.982520,0.990699,0.979244,0.907375 /)
XCGA_LKT(17,12,1:6)=(/ 0.724547,0.718980,0.693937,0.601607,0.364900,0.092370 /)
XEXT_COEFF_550_LKT(17,12)=2338.700000 !rg=0.0374736 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(17,13,1:6)=(/ 4004.400000,3364.100000,2432.200000,1157.700000,253.940000,26.367000 /)
XPIZA_LKT(17,13,1:6)=(/ 0.971112,0.979121,0.982619,0.992014,0.984958,0.942032 /)
XCGA_LKT(17,13,1:6)=(/ 0.726613,0.717980,0.712247,0.651833,0.467903,0.134290 /)
XEXT_COEFF_550_LKT(17,13)=2458.400000 !rg=0.0374736 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(17,14,1:6)=(/ 3371.500000,3005.200000,2383.200000,1353.200000,355.150000,42.652000 /)
XPIZA_LKT(17,14,1:6)=(/ 0.965587,0.976220,0.981580,0.992740,0.988478,0.962556 /)
XCGA_LKT(17,14,1:6)=(/ 0.727933,0.719313,0.718343,0.688530,0.543160,0.196130 /)
XEXT_COEFF_550_LKT(17,14)=2400.700000 !rg=0.0374736 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(17,15,1:6)=(/ 2806.400000,2638.900000,2198.900000,1475.100000,472.190000,65.955000 /)
XPIZA_LKT(17,15,1:6)=(/ 0.959287,0.972547,0.979492,0.992991,0.990695,0.974467 /)
XCGA_LKT(17,15,1:6)=(/ 0.730907,0.722833,0.717010,0.711730,0.612460,0.286653 /)
XEXT_COEFF_550_LKT(17,15)=2210.600000 !rg=0.0374736 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(17,16,1:6)=(/ 2306.000000,2239.200000,1972.100000,1493.800000,591.960000,98.747000 /)
XPIZA_LKT(17,16,1:6)=(/ 0.953099,0.967865,0.976540,0.992771,0.992092,0.981712 /)
XCGA_LKT(17,16,1:6)=(/ 0.738810,0.725117,0.715703,0.721837,0.663870,0.402960 /)
XEXT_COEFF_550_LKT(17,16)=1977.300000 !rg=0.0374736 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(17,17,1:6)=(/ 1898.200000,1886.600000,1744.400000,1410.000000,698.420000,144.510000 /)
XPIZA_LKT(17,17,1:6)=(/ 0.944411,0.962838,0.973154,0.992071,0.992912,0.986507 /)
XCGA_LKT(17,17,1:6)=(/ 0.744427,0.729613,0.719633,0.721303,0.702703,0.499947 /)
XEXT_COEFF_550_LKT(17,17)=1751.600000 !rg=0.0374736 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(17,18,1:6)=(/ 1555.500000,1578.800000,1502.700000,1272.700000,772.150000,197.110000 /)
XPIZA_LKT(17,18,1:6)=(/ 0.934636,0.955544,0.968951,0.990925,0.993273,0.989425 /)
XCGA_LKT(17,18,1:6)=(/ 0.750323,0.734307,0.723380,0.716023,0.728200,0.570060 /)
XEXT_COEFF_550_LKT(17,18)=1503.800000 !rg=0.0374736 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(17,19,1:6)=(/ 1266.300000,1302.700000,1275.700000,1131.100000,795.250000,256.590000 /)
XPIZA_LKT(17,19,1:6)=(/ 0.924775,0.948544,0.962995,0.989628,0.993203,0.991305 /)
XCGA_LKT(17,19,1:6)=(/ 0.762930,0.740737,0.727837,0.716517,0.741223,0.634260 /)
XEXT_COEFF_550_LKT(17,19)=1272.100000 !rg=0.0374736 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(17,20,1:6)=(/ 1039.300000,1070.800000,1068.800000,996.650000,762.380000,314.330000 /)
XPIZA_LKT(17,20,1:6)=(/ 0.912221,0.940844,0.958095,0.987994,0.992664,0.992465 /)
XCGA_LKT(17,20,1:6)=(/ 0.772137,0.748047,0.733320,0.721123,0.742533,0.680560 /)
XEXT_COEFF_550_LKT(17,20)=1073.800000 !rg=0.0374736 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(18,1,1:6)=(/ 2036.000000,384.160000,85.275000,14.784000,3.590800,1.392900 /)
XPIZA_LKT(18,1,1:6)=(/ 0.963632,0.899950,0.736821,0.640642,0.240375,0.037528 /)
XCGA_LKT(18,1,1:6)=(/ 0.262000,0.096687,0.041033,0.016213,0.005143,0.001267 /)
XEXT_COEFF_550_LKT(18,1)=82.511000 !rg=0.0405973 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(18,2,1:6)=(/ 2307.600000,466.180000,101.740000,17.280000,3.820500,1.407100 /)
XPIZA_LKT(18,2,1:6)=(/ 0.966954,0.915841,0.776630,0.690841,0.284772,0.046829 /)
XCGA_LKT(18,2,1:6)=(/ 0.331223,0.122697,0.052023,0.020583,0.006537,0.001610 /)
XEXT_COEFF_550_LKT(18,2)=99.801000 !rg=0.0405973 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(18,3,1:6)=(/ 2871.500000,638.530000,140.740000,23.347000,4.379300,1.441700 /)
XPIZA_LKT(18,3,1:6)=(/ 0.971736,0.936221,0.834466,0.768494,0.373678,0.068765 /)
XCGA_LKT(18,3,1:6)=(/ 0.442410,0.184683,0.078817,0.031300,0.009967,0.002460 /)
XEXT_COEFF_550_LKT(18,3)=140.680000 !rg=0.0405973 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(18,4,1:6)=(/ 3682.100000,917.430000,215.780000,36.006000,5.557100,1.514300 /)
XPIZA_LKT(18,4,1:6)=(/ 0.976385,0.953070,0.887910,0.846957,0.503207,0.111961 /)
XCGA_LKT(18,4,1:6)=(/ 0.537680,0.291433,0.132237,0.053133,0.017007,0.004207 /)
XEXT_COEFF_550_LKT(18,4)=218.930000 !rg=0.0405973 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(18,5,1:6)=(/ 4587.000000,1320.100000,340.180000,60.161000,7.874700,1.657400 /)
XPIZA_LKT(18,5,1:6)=(/ 0.979811,0.964918,0.925208,0.905749,0.645807,0.186441 /)
XCGA_LKT(18,5,1:6)=(/ 0.608090,0.410800,0.216913,0.088653,0.028530,0.007080 /)
XEXT_COEFF_550_LKT(18,5)=347.860000 !rg=0.0405973 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(18,6,1:6)=(/ 5422.900000,1829.500000,524.220000,101.350000,12.100000,1.919600 /)
XPIZA_LKT(18,6,1:6)=(/ 0.981954,0.972620,0.948213,0.941894,0.766028,0.294762 /)
XCGA_LKT(18,6,1:6)=(/ 0.658440,0.509997,0.322187,0.138397,0.044623,0.011110 /)
XEXT_COEFF_550_LKT(18,6)=537.800000 !rg=0.0405973 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(18,7,1:6)=(/ 6045.900000,2405.700000,778.620000,166.170000,19.599000,2.391700 /)
XPIZA_LKT(18,7,1:6)=(/ 0.983015,0.977539,0.962405,0.962825,0.852479,0.430606 /)
XCGA_LKT(18,7,1:6)=(/ 0.692947,0.585110,0.423397,0.207547,0.066827,0.016687 /)
XEXT_COEFF_550_LKT(18,7)=798.480000 !rg=0.0405973 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(18,8,1:6)=(/ 6308.000000,2973.300000,1095.200000,260.680000,32.414000,3.225200 /)
XPIZA_LKT(18,8,1:6)=(/ 0.983049,0.980564,0.971064,0.974835,0.908258,0.574134 /)
XCGA_LKT(18,8,1:6)=(/ 0.713363,0.640447,0.510300,0.301070,0.097903,0.024483 /)
XEXT_COEFF_550_LKT(18,8)=1121.900000 !rg=0.0405973 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(18,9,1:6)=(/ 6164.000000,3451.400000,1461.100000,397.850000,53.653000,4.704300 /)
XPIZA_LKT(18,9,1:6)=(/ 0.982059,0.982252,0.976578,0.982223,0.942534,0.704436 /)
XCGA_LKT(18,9,1:6)=(/ 0.722323,0.679660,0.581650,0.404870,0.142520,0.035607 /)
XEXT_COEFF_550_LKT(18,9)=1492.700000 !rg=0.0405973 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(18,10,1:6)=(/ 5696.000000,3737.900000,1833.600000,575.870000,86.344000,7.296400 /)
XPIZA_LKT(18,10,1:6)=(/ 0.980125,0.982791,0.979974,0.986751,0.962687,0.806145 /)
XCGA_LKT(18,10,1:6)=(/ 0.724263,0.704440,0.636137,0.488770,0.207290,0.051557 /)
XEXT_COEFF_550_LKT(18,10)=1868.300000 !rg=0.0405973 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(18,11,1:6)=(/ 5084.600000,3773.300000,2162.200000,789.810000,133.680000,11.843000 /)
XPIZA_LKT(18,11,1:6)=(/ 0.977301,0.982260,0.981927,0.989542,0.974545,0.877703 /)
XCGA_LKT(18,11,1:6)=(/ 0.724257,0.716853,0.676303,0.563727,0.299963,0.074660 /)
XEXT_COEFF_550_LKT(18,11)=2196.700000 !rg=0.0405973 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(18,12,1:6)=(/ 4423.000000,3568.400000,2380.000000,1021.400000,201.640000,19.618000 /)
XPIZA_LKT(18,12,1:6)=(/ 0.973686,0.980622,0.982703,0.991307,0.981886,0.923811 /)
XCGA_LKT(18,12,1:6)=(/ 0.726290,0.719230,0.702700,0.622693,0.411113,0.108290 /)
XEXT_COEFF_550_LKT(18,12)=2410.000000 !rg=0.0405973 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(18,13,1:6)=(/ 3742.800000,3229.500000,2435.100000,1243.200000,294.060000,32.313000 /)
XPIZA_LKT(18,13,1:6)=(/ 0.968864,0.978104,0.982345,0.992374,0.986634,0.951878 /)
XCGA_LKT(18,13,1:6)=(/ 0.726590,0.719017,0.715910,0.668113,0.500117,0.157683 /)
XEXT_COEFF_550_LKT(18,13)=2457.500000 !rg=0.0405973 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(18,14,1:6)=(/ 3129.900000,2869.200000,2323.200000,1413.500000,401.510000,51.382000 /)
XPIZA_LKT(18,14,1:6)=(/ 0.963039,0.974872,0.980867,0.992902,0.989497,0.968273 /)
XCGA_LKT(18,14,1:6)=(/ 0.729057,0.720917,0.718460,0.699107,0.573093,0.230617 /)
XEXT_COEFF_550_LKT(18,14)=2337.400000 !rg=0.0405973 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(18,15,1:6)=(/ 2590.500000,2477.500000,2113.600000,1495.900000,521.210000,78.055000 /)
XPIZA_LKT(18,15,1:6)=(/ 0.956460,0.970831,0.978333,0.992955,0.991350,0.977846 /)
XCGA_LKT(18,15,1:6)=(/ 0.734773,0.724427,0.715877,0.716887,0.634823,0.334210 /)
XEXT_COEFF_550_LKT(18,15)=2120.600000 !rg=0.0405973 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(18,16,1:6)=(/ 2125.200000,2087.100000,1880.600000,1470.500000,638.760000,116.370000 /)
XPIZA_LKT(18,16,1:6)=(/ 0.949594,0.966195,0.975355,0.992542,0.992500,0.983959 /)
XCGA_LKT(18,16,1:6)=(/ 0.740793,0.727600,0.717467,0.722557,0.681437,0.448840 /)
XEXT_COEFF_550_LKT(18,16)=1891.500000 !rg=0.0405973 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(18,17,1:6)=(/ 1744.200000,1759.900000,1649.800000,1358.100000,734.070000,165.790000 /)
XPIZA_LKT(18,17,1:6)=(/ 0.940022,0.959750,0.971584,0.991628,0.993121,0.987912 /)
XCGA_LKT(18,17,1:6)=(/ 0.746107,0.731833,0.722103,0.718720,0.714477,0.529577 /)
XEXT_COEFF_550_LKT(18,17)=1652.200000 !rg=0.0405973 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(18,18,1:6)=(/ 1427.900000,1455.800000,1408.600000,1213.300000,788.020000,221.530000 /)
XPIZA_LKT(18,18,1:6)=(/ 0.929947,0.953030,0.966217,0.990421,0.993297,0.990301 /)
XCGA_LKT(18,18,1:6)=(/ 0.756850,0.735657,0.724673,0.715620,0.734680,0.599350 /)
XEXT_COEFF_550_LKT(18,18)=1403.800000 !rg=0.0405973 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(18,19,1:6)=(/ 1163.000000,1196.600000,1186.800000,1078.300000,787.810000,280.800000 /)
XPIZA_LKT(18,19,1:6)=(/ 0.919637,0.945882,0.961301,0.989006,0.993029,0.991853 /)
XCGA_LKT(18,19,1:6)=(/ 0.766757,0.744187,0.730067,0.719217,0.742700,0.654610 /)
XEXT_COEFF_550_LKT(18,19)=1192.200000 !rg=0.0405973 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(18,20,1:6)=(/ 952.390000,989.570000,988.140000,938.380000,735.720000,335.920000 /)
XPIZA_LKT(18,20,1:6)=(/ 0.905887,0.935907,0.955454,0.987051,0.992301,0.992796 /)
XCGA_LKT(18,20,1:6)=(/ 0.774550,0.752037,0.734687,0.722127,0.740680,0.696357 /)
XEXT_COEFF_550_LKT(18,20)=989.410000 !rg=0.0405973 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(19,1,1:6)=(/ 2398.900000,474.650000,102.600000,17.391000,3.830800,1.407800 /)
XPIZA_LKT(19,1,1:6)=(/ 0.968067,0.916815,0.777898,0.692473,0.286472,0.047212 /)
XCGA_LKT(19,1,1:6)=(/ 0.313953,0.113600,0.048097,0.019013,0.006037,0.001487 /)
XEXT_COEFF_550_LKT(19,1)=100.700000 !rg=0.0439813 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(19,2,1:6)=(/ 2696.400000,573.190000,123.340000,20.560000,4.122300,1.425800 /)
XPIZA_LKT(19,2,1:6)=(/ 0.970516,0.929625,0.812631,0.738194,0.335550,0.058761 /)
XCGA_LKT(19,2,1:6)=(/ 0.390697,0.144290,0.060970,0.024133,0.007670,0.001890 /)
XEXT_COEFF_550_LKT(19,2)=122.480000 !rg=0.0439813 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(19,3,1:6)=(/ 3329.500000,772.880000,171.910000,28.250000,4.831600,1.469600 /)
XPIZA_LKT(19,3,1:6)=(/ 0.974506,0.945721,0.861821,0.806891,0.430585,0.085771 /)
XCGA_LKT(19,3,1:6)=(/ 0.492647,0.216400,0.092297,0.036680,0.011690,0.002887 /)
XEXT_COEFF_550_LKT(19,3)=173.330000 !rg=0.0439813 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(19,4,1:6)=(/ 4169.400000,1089.000000,262.750000,44.205000,6.325900,1.561700 /)
XPIZA_LKT(19,4,1:6)=(/ 0.978324,0.959101,0.905835,0.873844,0.561795,0.138028 /)
XCGA_LKT(19,4,1:6)=(/ 0.572550,0.331630,0.154223,0.062180,0.019937,0.004937 /)
XEXT_COEFF_550_LKT(19,4)=267.840000 !rg=0.0439813 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(19,5,1:6)=(/ 5047.200000,1534.800000,409.080000,74.216000,9.262600,1.743200 /)
XPIZA_LKT(19,5,1:6)=(/ 0.981058,0.968701,0.936100,0.922419,0.697156,0.225358 /)
XCGA_LKT(19,5,1:6)=(/ 0.632573,0.448717,0.249683,0.103620,0.033423,0.008307 /)
XEXT_COEFF_550_LKT(19,5)=419.210000 !rg=0.0439813 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(19,6,1:6)=(/ 5795.200000,2075.200000,620.730000,124.140000,14.601000,2.076100 /)
XPIZA_LKT(19,6,1:6)=(/ 0.982664,0.974996,0.954876,0.951650,0.804556,0.346477 /)
XCGA_LKT(19,6,1:6)=(/ 0.675533,0.540327,0.358563,0.161700,0.052253,0.013030 /)
XEXT_COEFF_550_LKT(19,6)=637.000000 !rg=0.0439813 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(19,7,1:6)=(/ 6270.100000,2660.300000,902.120000,200.560000,24.013000,2.675100 /)
XPIZA_LKT(19,7,1:6)=(/ 0.983271,0.979036,0.966456,0.968462,0.878292,0.489307 /)
XCGA_LKT(19,7,1:6)=(/ 0.704007,0.608647,0.456723,0.241807,0.078243,0.019563 /)
XEXT_COEFF_550_LKT(19,7)=924.950000 !rg=0.0439813 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(19,8,1:6)=(/ 6350.900000,3201.800000,1243.800000,311.200000,39.882000,3.732100 /)
XPIZA_LKT(19,8,1:6)=(/ 0.982882,0.981443,0.973648,0.978275,0.924388,0.630297 /)
XCGA_LKT(19,8,1:6)=(/ 0.719443,0.657657,0.540033,0.343257,0.114670,0.028697 /)
XEXT_COEFF_550_LKT(19,8)=1273.000000 !rg=0.0439813 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(19,9,1:6)=(/ 6043.400000,3612.400000,1617.800000,467.170000,65.545000,5.605000 /)
XPIZA_LKT(19,9,1:6)=(/ 0.981500,0.982639,0.978170,0.984372,0.952141,0.750332 /)
XCGA_LKT(19,9,1:6)=(/ 0.724950,0.691320,0.604357,0.440103,0.167083,0.041723 /)
XEXT_COEFF_550_LKT(19,9)=1651.700000 !rg=0.0439813 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(19,10,1:6)=(/ 5494.100000,3797.900000,1978.000000,658.540000,103.760000,8.876200 /)
XPIZA_LKT(19,10,1:6)=(/ 0.979107,0.982740,0.980925,0.988035,0.968291,0.839223 /)
XCGA_LKT(19,10,1:6)=(/ 0.724980,0.710933,0.653423,0.519100,0.242867,0.060400 /)
XEXT_COEFF_550_LKT(19,10)=2013.600000 !rg=0.0439813 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(19,11,1:6)=(/ 4840.300000,3728.300000,2269.400000,884.530000,158.510000,14.571000 /)
XPIZA_LKT(19,11,1:6)=(/ 0.976025,0.981753,0.982385,0.990372,0.977935,0.899403 /)
XCGA_LKT(19,11,1:6)=(/ 0.726367,0.718923,0.688030,0.589483,0.346267,0.087487 /)
XEXT_COEFF_550_LKT(19,11)=2303.200000 !rg=0.0439813 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(19,12,1:6)=(/ 4157.200000,3450.900000,2426.300000,1116.300000,237.220000,24.155000 /)
XPIZA_LKT(19,12,1:6)=(/ 0.971773,0.979792,0.982696,0.991810,0.984114,0.937169 /)
XCGA_LKT(19,12,1:6)=(/ 0.726500,0.720017,0.709130,0.642397,0.451480,0.127017 /)
XEXT_COEFF_550_LKT(19,12)=2454.200000 !rg=0.0439813 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(19,13,1:6)=(/ 3478.700000,3097.300000,2410.800000,1322.200000,336.040000,39.353000 /)
XPIZA_LKT(19,13,1:6)=(/ 0.966555,0.976850,0.981894,0.992638,0.987974,0.959750 /)
XCGA_LKT(19,13,1:6)=(/ 0.727387,0.719673,0.717960,0.682047,0.529793,0.185280 /)
XEXT_COEFF_550_LKT(19,13)=2430.200000 !rg=0.0439813 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(19,14,1:6)=(/ 2905.300000,2715.600000,2249.800000,1460.900000,450.800000,61.364000 /)
XPIZA_LKT(19,14,1:6)=(/ 0.960166,0.973501,0.979930,0.992972,0.990366,0.972829 /)
XCGA_LKT(19,14,1:6)=(/ 0.732173,0.723627,0.717400,0.707827,0.601127,0.270787 /)
XEXT_COEFF_550_LKT(19,14)=2260.500000 !rg=0.0439813 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(19,15,1:6)=(/ 2383.400000,2316.300000,2021.300000,1498.600000,570.320000,92.157000 /)
XPIZA_LKT(19,15,1:6)=(/ 0.953827,0.969157,0.977233,0.992844,0.991881,0.980662 /)
XCGA_LKT(19,15,1:6)=(/ 0.736630,0.726197,0.716343,0.720393,0.654983,0.384207 /)
XEXT_COEFF_550_LKT(19,15)=2033.100000 !rg=0.0439813 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(19,16,1:6)=(/ 1950.400000,1952.500000,1792.100000,1433.000000,681.430000,136.140000 /)
XPIZA_LKT(19,16,1:6)=(/ 0.945774,0.963422,0.973964,0.992220,0.992802,0.985850 /)
XCGA_LKT(19,16,1:6)=(/ 0.742313,0.728993,0.720277,0.721520,0.696470,0.486567 /)
XEXT_COEFF_550_LKT(19,16)=1797.400000 !rg=0.0439813 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(19,17,1:6)=(/ 1603.700000,1625.100000,1552.300000,1298.800000,762.390000,187.850000 /)
XPIZA_LKT(19,17,1:6)=(/ 0.935561,0.956850,0.969569,0.991195,0.993238,0.989033 /)
XCGA_LKT(19,17,1:6)=(/ 0.751253,0.731910,0.723187,0.717430,0.724190,0.558403 /)
XEXT_COEFF_550_LKT(19,17)=1550.300000 !rg=0.0439813 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(19,18,1:6)=(/ 1305.200000,1341.500000,1311.800000,1159.200000,794.970000,246.600000 /)
XPIZA_LKT(19,18,1:6)=(/ 0.925633,0.950843,0.964299,0.989881,0.993236,0.991045 /)
XCGA_LKT(19,18,1:6)=(/ 0.760567,0.741167,0.725757,0.717003,0.739347,0.625093 /)
XEXT_COEFF_550_LKT(19,18)=1315.800000 !rg=0.0439813 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(19,19,1:6)=(/ 1066.700000,1107.100000,1101.000000,1024.600000,771.110000,305.050000 /)
XPIZA_LKT(19,19,1:6)=(/ 0.913614,0.941588,0.959349,0.988206,0.992779,0.992305 /)
XCGA_LKT(19,19,1:6)=(/ 0.769370,0.747243,0.732903,0.720183,0.742720,0.673493 /)
XEXT_COEFF_550_LKT(19,19)=1104.100000 !rg=0.0439813 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(19,20,1:6)=(/ 874.710000,907.350000,918.930000,876.120000,704.810000,355.530000 /)
XPIZA_LKT(19,20,1:6)=(/ 0.899320,0.931724,0.951856,0.986029,0.991853,0.993032 /)
XCGA_LKT(19,20,1:6)=(/ 0.781093,0.753317,0.737943,0.721953,0.737553,0.709740 /)
XEXT_COEFF_550_LKT(19,20)=914.540000 !rg=0.0439813 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(20,1,1:6)=(/ 2785.400000,585.850000,124.500000,20.701000,4.135400,1.426700 /)
XPIZA_LKT(20,1,1:6)=(/ 0.971407,0.930641,0.813779,0.739650,0.337400,0.059235 /)
XCGA_LKT(20,1,1:6)=(/ 0.377100,0.133677,0.056380,0.022297,0.007080,0.001747 /)
XEXT_COEFF_550_LKT(20,1)=123.700000 !rg=0.0476474 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(20,2,1:6)=(/ 3134.500000,701.720000,150.540000,24.723000,4.505400,1.449500 /)
XPIZA_LKT(20,2,1:6)=(/ 0.973380,0.940807,0.843558,0.780373,0.390337,0.073492 /)
XCGA_LKT(20,2,1:6)=(/ 0.454770,0.169980,0.071467,0.028290,0.008997,0.002220 /)
XEXT_COEFF_550_LKT(20,2)=151.020000 !rg=0.0476474 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(20,3,1:6)=(/ 3832.000000,928.790000,210.560000,34.461000,5.405600,1.505000 /)
XPIZA_LKT(20,3,1:6)=(/ 0.976832,0.953382,0.884740,0.839999,0.489248,0.106488 /)
XCGA_LKT(20,3,1:6)=(/ 0.538510,0.253233,0.108093,0.042973,0.013710,0.003387 /)
XEXT_COEFF_550_LKT(20,3)=213.760000 !rg=0.0476474 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(20,4,1:6)=(/ 4672.800000,1284.000000,319.240000,54.502000,7.301000,1.621900 /)
XPIZA_LKT(20,4,1:6)=(/ 0.979944,0.964033,0.920564,0.896283,0.618513,0.168983 /)
XCGA_LKT(20,4,1:6)=(/ 0.603630,0.373983,0.179557,0.072737,0.023367,0.005793 /)
XEXT_COEFF_550_LKT(20,4)=326.540000 !rg=0.0476474 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(20,5,1:6)=(/ 5492.800000,1769.900000,489.720000,91.498000,11.019000,1.852200 /)
XPIZA_LKT(20,5,1:6)=(/ 0.982073,0.971830,0.945045,0.935992,0.743747,0.269650 /)
XCGA_LKT(20,5,1:6)=(/ 0.654407,0.485243,0.284570,0.121033,0.039147,0.009743 /)
XEXT_COEFF_550_LKT(20,5)=502.500000 !rg=0.0476474 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(20,6,1:6)=(/ 6122.400000,2333.800000,729.160000,151.340000,17.749000,2.274500 /)
XPIZA_LKT(20,6,1:6)=(/ 0.983192,0.976977,0.960322,0.959502,0.837758,0.401988 /)
XCGA_LKT(20,6,1:6)=(/ 0.690433,0.568747,0.394113,0.188670,0.061177,0.015280 /)
XEXT_COEFF_550_LKT(20,6)=748.220000 !rg=0.0476474 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(20,7,1:6)=(/ 6427.800000,2911.800000,1038.200000,241.000000,29.507000,3.034500 /)
XPIZA_LKT(20,7,1:6)=(/ 0.983358,0.980251,0.969825,0.973055,0.899752,0.548138 /)
XCGA_LKT(20,7,1:6)=(/ 0.713170,0.629840,0.489723,0.279780,0.091597,0.022933 /)
XEXT_COEFF_550_LKT(20,7)=1064.200000 !rg=0.0476474 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(20,8,1:6)=(/ 6326.900000,3408.800000,1397.500000,369.910000,48.983000,4.374300 /)
XPIZA_LKT(20,8,1:6)=(/ 0.982538,0.982118,0.975777,0.981135,0.937485,0.682915 /)
XCGA_LKT(20,8,1:6)=(/ 0.723940,0.672690,0.566487,0.383273,0.134317,0.033630 /)
XEXT_COEFF_550_LKT(20,8)=1428.700000 !rg=0.0476474 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(20,9,1:6)=(/ 5883.900000,3736.400000,1775.300000,541.800000,79.571000,6.743200 /)
XPIZA_LKT(20,9,1:6)=(/ 0.980738,0.982854,0.979502,0.986099,0.959832,0.790937 /)
XCGA_LKT(20,9,1:6)=(/ 0.726780,0.700800,0.625410,0.471827,0.195847,0.048880 /)
XEXT_COEFF_550_LKT(20,9)=1810.200000 !rg=0.0476474 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(20,10,1:6)=(/ 5259.900000,3816.100000,2113.900000,748.460000,123.840000,10.861000 /)
XPIZA_LKT(20,10,1:6)=(/ 0.978103,0.982516,0.981665,0.989114,0.972806,0.867267 /)
XCGA_LKT(20,10,1:6)=(/ 0.726830,0.715693,0.668607,0.549173,0.283453,0.070760 /)
XEXT_COEFF_550_LKT(20,10)=2149.300000 !rg=0.0476474 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(20,11,1:6)=(/ 4586.400000,3647.900000,2357.000000,978.790000,187.630000,17.954000 /)
XPIZA_LKT(20,11,1:6)=(/ 0.974432,0.981138,0.982648,0.991039,0.980789,0.917265 /)
XCGA_LKT(20,11,1:6)=(/ 0.726977,0.720497,0.697953,0.611620,0.392407,0.102540 /)
XEXT_COEFF_550_LKT(20,11)=2388.500000 !rg=0.0476474 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(20,12,1:6)=(/ 3873.300000,3326.700000,2444.700000,1205.400000,275.980000,29.643000 /)
XPIZA_LKT(20,12,1:6)=(/ 0.969831,0.978734,0.982517,0.992217,0.985941,0.947945 /)
XCGA_LKT(20,12,1:6)=(/ 0.727130,0.720110,0.713927,0.659793,0.485280,0.149070 /)
XEXT_COEFF_550_LKT(20,12)=2469.200000 !rg=0.0476474 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(20,13,1:6)=(/ 3241.700000,2951.300000,2363.900000,1389.300000,380.960000,47.551000 /)
XPIZA_LKT(20,13,1:6)=(/ 0.963807,0.975698,0.981236,0.992838,0.989072,0.966019 /)
XCGA_LKT(20,13,1:6)=(/ 0.729983,0.722077,0.718503,0.693683,0.559843,0.217777 /)
XEXT_COEFF_550_LKT(20,13)=2379.400000 !rg=0.0476474 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(20,14,1:6)=(/ 2670.500000,2559.900000,2162.500000,1489.700000,500.020000,72.780000 /)
XPIZA_LKT(20,14,1:6)=(/ 0.957760,0.971906,0.978996,0.992978,0.991082,0.976505 /)
XCGA_LKT(20,14,1:6)=(/ 0.733200,0.725313,0.717207,0.714013,0.624837,0.316393 /)
XEXT_COEFF_550_LKT(20,14)=2176.100000 !rg=0.0476474 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(20,15,1:6)=(/ 2192.400000,2162.400000,1933.800000,1484.300000,618.510000,108.730000 /)
XPIZA_LKT(20,15,1:6)=(/ 0.950925,0.966891,0.976026,0.992646,0.992328,0.983066 /)
XCGA_LKT(20,15,1:6)=(/ 0.739137,0.726490,0.718137,0.721810,0.673490,0.431627 /)
XEXT_COEFF_550_LKT(20,15)=1940.400000 !rg=0.0476474 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(20,16,1:6)=(/ 1800.300000,1812.100000,1700.500000,1383.200000,719.970000,157.010000 /)
XPIZA_LKT(20,16,1:6)=(/ 0.941223,0.960881,0.972201,0.991859,0.993037,0.987381 /)
XCGA_LKT(20,16,1:6)=(/ 0.746597,0.730233,0.721403,0.720367,0.709220,0.517553 /)
XEXT_COEFF_550_LKT(20,16)=1702.200000 !rg=0.0476474 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(20,17,1:6)=(/ 1465.000000,1501.500000,1448.900000,1243.100000,782.400000,211.610000 /)
XPIZA_LKT(20,17,1:6)=(/ 0.931581,0.954074,0.967402,0.990679,0.993292,0.989966 /)
XCGA_LKT(20,17,1:6)=(/ 0.754030,0.736667,0.723840,0.716727,0.731767,0.588023 /)
XEXT_COEFF_550_LKT(20,17)=1450.900000 !rg=0.0476474 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(20,18,1:6)=(/ 1194.700000,1237.300000,1222.000000,1106.000000,791.280000,270.920000 /)
XPIZA_LKT(20,18,1:6)=(/ 0.921238,0.946992,0.962918,0.989287,0.993100,0.991643 /)
XCGA_LKT(20,18,1:6)=(/ 0.764420,0.742990,0.730600,0.718143,0.741727,0.646423 /)
XEXT_COEFF_550_LKT(20,18)=1226.500000 !rg=0.0476474 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(20,19,1:6)=(/ 980.470000,1017.000000,1022.400000,962.950000,746.950000,327.550000 /)
XPIZA_LKT(20,19,1:6)=(/ 0.907548,0.937179,0.956101,0.987536,0.992446,0.992674 /)
XCGA_LKT(20,19,1:6)=(/ 0.775210,0.748680,0.734273,0.722477,0.741000,0.690210 /)
XEXT_COEFF_550_LKT(20,19)=1019.000000 !rg=0.0476474 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(20,20,1:6)=(/ 798.030000,832.520000,848.410000,817.220000,671.730000,371.750000 /)
XPIZA_LKT(20,20,1:6)=(/ 0.892773,0.927031,0.948477,0.985379,0.991388,0.993211 /)
XCGA_LKT(20,20,1:6)=(/ 0.784923,0.759953,0.738920,0.726457,0.735663,0.720637 /)
XEXT_COEFF_550_LKT(20,20)=848.360000 !rg=0.0476474 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(21,1,1:6)=(/ 3204.500000,720.600000,152.130000,24.904000,4.522000,1.450600 /)
XPIZA_LKT(21,1,1:6)=(/ 0.973903,0.941885,0.844612,0.781648,0.392299,0.074074 /)
XCGA_LKT(21,1,1:6)=(/ 0.450850,0.157647,0.066100,0.026143,0.008307,0.002050 /)
XEXT_COEFF_550_LKT(21,1)=152.720000 !rg=0.0516191 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(21,2,1:6)=(/ 3638.000000,853.190000,184.660000,30.005000,4.991700,1.479600 /)
XPIZA_LKT(21,2,1:6)=(/ 0.975794,0.949807,0.869752,0.817212,0.447917,0.091546 /)
XCGA_LKT(21,2,1:6)=(/ 0.517100,0.200687,0.083797,0.033160,0.010553,0.002603 /)
XEXT_COEFF_550_LKT(21,2)=186.800000 !rg=0.0516191 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(21,3,1:6)=(/ 4370.500000,1107.500000,258.090000,42.318000,6.134000,1.550000 /)
XPIZA_LKT(21,3,1:6)=(/ 0.978805,0.959551,0.903731,0.868107,0.548045,0.131470 /)
XCGA_LKT(21,3,1:6)=(/ 0.578200,0.295383,0.126613,0.050337,0.016077,0.003973 /)
XEXT_COEFF_550_LKT(21,3)=263.370000 !rg=0.0516191 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(21,4,1:6)=(/ 5178.700000,1503.000000,386.490000,67.374000,8.537300,1.698200 /)
XPIZA_LKT(21,4,1:6)=(/ 0.981292,0.968081,0.932611,0.914811,0.671952,0.205190 /)
XCGA_LKT(21,4,1:6)=(/ 0.631010,0.417257,0.208440,0.085047,0.027383,0.006793 /)
XEXT_COEFF_550_LKT(21,4)=396.260000 !rg=0.0516191 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(21,5,1:6)=(/ 5910.000000,2023.600000,582.780000,112.560000,13.238000,1.990500 /)
XPIZA_LKT(21,5,1:6)=(/ 0.982885,0.974424,0.952368,0.946982,0.785089,0.318988 /)
XCGA_LKT(21,5,1:6)=(/ 0.673723,0.519900,0.320830,0.141253,0.045840,0.011427 /)
XEXT_COEFF_550_LKT(21,5)=598.370000 !rg=0.0516191 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(21,6,1:6)=(/ 6393.100000,2599.200000,850.100000,183.570000,21.697000,2.526400 /)
XPIZA_LKT(21,6,1:6)=(/ 0.983549,0.978626,0.964794,0.965835,0.865914,0.460002 /)
XCGA_LKT(21,6,1:6)=(/ 0.703337,0.594903,0.429437,0.219510,0.071603,0.017913 /)
XEXT_COEFF_550_LKT(21,6)=872.250000 !rg=0.0516191 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(21,7,1:6)=(/ 6511.900000,3155.800000,1185.400000,288.600000,36.297000,3.490300 /)
XPIZA_LKT(21,7,1:6)=(/ 0.983291,0.981242,0.972642,0.976838,0.917401,0.605459 /)
XCGA_LKT(21,7,1:6)=(/ 0.720703,0.648917,0.520937,0.319583,0.107217,0.026880 /)
XEXT_COEFF_550_LKT(21,7)=1213.900000 !rg=0.0516191 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(21,8,1:6)=(/ 6235.400000,3591.100000,1555.900000,435.870000,59.936000,5.187200 /)
XPIZA_LKT(21,8,1:6)=(/ 0.982070,0.982605,0.977527,0.983475,0.948045,0.730982 /)
XCGA_LKT(21,8,1:6)=(/ 0.727590,0.685860,0.590710,0.419217,0.157323,0.039407 /)
XEXT_COEFF_550_LKT(21,8)=1589.500000 !rg=0.0516191 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(21,9,1:6)=(/ 5681.600000,3820.900000,1925.900000,622.260000,95.907000,8.178500 /)
XPIZA_LKT(21,9,1:6)=(/ 0.979852,0.982901,0.980568,0.987500,0.965987,0.826168 /)
XCGA_LKT(21,9,1:6)=(/ 0.728033,0.708630,0.644000,0.502827,0.229297,0.057263 /)
XEXT_COEFF_550_LKT(21,9)=1961.700000 !rg=0.0516191 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(21,10,1:6)=(/ 5016.600000,3791.900000,2232.100000,842.440000,147.160000,13.342000 /)
XPIZA_LKT(21,10,1:6)=(/ 0.976860,0.982134,0.982224,0.990017,0.976505,0.890720 /)
XCGA_LKT(21,10,1:6)=(/ 0.728227,0.719017,0.681547,0.576213,0.327917,0.082900 /)
XEXT_COEFF_550_LKT(21,10)=2266.700000 !rg=0.0516191 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(21,11,1:6)=(/ 4300.200000,3547.600000,2417.600000,1074.500000,221.290000,22.112000 /)
XPIZA_LKT(21,11,1:6)=(/ 0.972810,0.980301,0.982742,0.991588,0.983198,0.931831 /)
XCGA_LKT(21,11,1:6)=(/ 0.728087,0.721083,0.705537,0.632173,0.434000,0.120227 /)
XEXT_COEFF_550_LKT(21,11)=2446.900000 !rg=0.0516191 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(21,12,1:6)=(/ 3611.100000,3187.400000,2435.800000,1288.100000,316.800000,36.181000 /)
XPIZA_LKT(21,12,1:6)=(/ 0.967326,0.977701,0.982153,0.992519,0.987405,0.956580 /)
XCGA_LKT(21,12,1:6)=(/ 0.728610,0.721503,0.716943,0.674717,0.515607,0.175063 /)
XEXT_COEFF_550_LKT(21,12)=2456.400000 !rg=0.0516191 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(21,13,1:6)=(/ 2992.300000,2805.300000,2295.600000,1443.300000,429.330000,56.966000 /)
XPIZA_LKT(21,13,1:6)=(/ 0.961702,0.974342,0.980486,0.992943,0.990005,0.971014 /)
XCGA_LKT(21,13,1:6)=(/ 0.731167,0.724183,0.718690,0.703297,0.588910,0.255750 /)
XEXT_COEFF_550_LKT(21,13)=2310.800000 !rg=0.0516191 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(21,14,1:6)=(/ 2455.000000,2395.900000,2076.700000,1501.800000,549.010000,86.018000 /)
XPIZA_LKT(21,14,1:6)=(/ 0.955384,0.969963,0.977884,0.992901,0.991658,0.979541 /)
XCGA_LKT(21,14,1:6)=(/ 0.735847,0.725723,0.717563,0.718583,0.645687,0.365487 /)
XEXT_COEFF_550_LKT(21,14)=2084.500000 !rg=0.0516191 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(21,15,1:6)=(/ 2021.400000,2013.400000,1844.800000,1452.900000,662.650000,127.640000 /)
XPIZA_LKT(21,15,1:6)=(/ 0.946647,0.964664,0.974674,0.992381,0.992671,0.985106 /)
XCGA_LKT(21,15,1:6)=(/ 0.741910,0.728060,0.719537,0.721970,0.689453,0.471827 /)
XEXT_COEFF_550_LKT(21,15)=1851.500000 !rg=0.0516191 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(21,16,1:6)=(/ 1645.900000,1681.700000,1597.700000,1330.600000,751.070000,178.590000 /)
XPIZA_LKT(21,16,1:6)=(/ 0.937365,0.957657,0.970619,0.991392,0.993192,0.988599 /)
XCGA_LKT(21,16,1:6)=(/ 0.747913,0.733200,0.723617,0.718410,0.719670,0.546407 /)
XEXT_COEFF_550_LKT(21,16)=1600.300000 !rg=0.0516191 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(21,17,1:6)=(/ 1337.800000,1384.200000,1353.000000,1187.000000,793.160000,236.590000 /)
XPIZA_LKT(21,17,1:6)=(/ 0.927954,0.951384,0.965504,0.990138,0.993262,0.990764 /)
XCGA_LKT(21,17,1:6)=(/ 0.759963,0.739130,0.727340,0.716710,0.737053,0.615240 /)
XEXT_COEFF_550_LKT(21,17)=1352.700000 !rg=0.0516191 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(21,18,1:6)=(/ 1099.300000,1138.200000,1138.300000,1049.900000,778.530000,295.220000 /)
XPIZA_LKT(21,18,1:6)=(/ 0.915219,0.942795,0.960294,0.988665,0.992884,0.992127 /)
XCGA_LKT(21,18,1:6)=(/ 0.769327,0.744633,0.731990,0.721387,0.742380,0.665823 /)
XEXT_COEFF_550_LKT(21,18)=1141.200000 !rg=0.0516191 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(21,19,1:6)=(/ 895.210000,935.590000,944.890000,904.240000,717.560000,347.850000 /)
XPIZA_LKT(21,19,1:6)=(/ 0.901599,0.932533,0.952818,0.986194,0.992056,0.992939 /)
XCGA_LKT(21,19,1:6)=(/ 0.777837,0.754090,0.735103,0.723250,0.739117,0.704397 /)
XEXT_COEFF_550_LKT(21,19)=943.650000 !rg=0.0516191 sigma=2.85 wvl=0.55
 
IF (LHOOK) CALL DR_HOOK('MODE_DUSTOPT:DUST_OPT_LKT_SET2',1,ZHOOK_HANDLE)
END SUBROUTINE DUST_OPT_LKT_SET2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE DUST_OPT_LKT_SET3()

  USE MODD_DUST_OPT_LKT
  
  IMPLICIT NONE

REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_DUSTOPT:DUST_OPT_LKT_SET3',0,ZHOOK_HANDLE)
XEXT_COEFF_WVL_LKT(21,20,1:6)=(/ 728.630000,762.500000,781.350000,762.600000,640.740000,384.220000 /)
XPIZA_LKT(21,20,1:6)=(/ 0.886500,0.922634,0.945140,0.984554,0.990860,0.993294 /)
XCGA_LKT(21,20,1:6)=(/ 0.791770,0.763860,0.744493,0.727810,0.734203,0.729583 /)
XEXT_COEFF_550_LKT(21,20)=780.930000 !rg=0.0516191 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(22,1,1:6)=(/ 3695.400000,880.730000,186.890000,30.239000,5.012700,1.480900 /)
XPIZA_LKT(22,1,1:6)=(/ 0.975903,0.950951,0.870743,0.818312,0.449941,0.092256 /)
XCGA_LKT(22,1,1:6)=(/ 0.528673,0.186480,0.077520,0.030647,0.009747,0.002403 /)
XEXT_COEFF_550_LKT(22,1)=189.200000 !rg=0.0559219 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(22,2,1:6)=(/ 4214.900000,1027.800000,227.220000,36.703000,5.608800,1.517700 /)
XPIZA_LKT(22,2,1:6)=(/ 0.977943,0.956989,0.891669,0.848839,0.506782,0.113478 /)
XCGA_LKT(22,2,1:6)=(/ 0.569943,0.237520,0.098300,0.038867,0.012380,0.003057 /)
XEXT_COEFF_550_LKT(22,2)=231.380000 !rg=0.0559219 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(22,3,1:6)=(/ 4929.100000,1310.200000,315.900000,52.237000,7.058400,1.606900 /)
XPIZA_LKT(22,3,1:6)=(/ 0.980479,0.964530,0.919328,0.891662,0.605349,0.161232 /)
XCGA_LKT(22,3,1:6)=(/ 0.611410,0.342427,0.148347,0.058957,0.018850,0.004663 /)
XEXT_COEFF_550_LKT(22,3)=323.570000 !rg=0.0559219 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(22,4,1:6)=(/ 5672.200000,1746.300000,465.620000,83.376000,10.104000,1.795100 /)
XPIZA_LKT(22,4,1:6)=(/ 0.982401,0.971418,0.942434,0.929980,0.721038,0.246803 /)
XCGA_LKT(22,4,1:6)=(/ 0.655040,0.460080,0.240920,0.099387,0.032083,0.007970 /)
XEXT_COEFF_550_LKT(22,4)=478.120000 !rg=0.0559219 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(22,5,1:6)=(/ 6284.200000,2293.600000,688.840000,137.990000,16.038000,2.166000 /)
XPIZA_LKT(22,5,1:6)=(/ 0.983508,0.976585,0.958359,0.955851,0.821056,0.372652 /)
XCGA_LKT(22,5,1:6)=(/ 0.690667,0.552363,0.358123,0.164637,0.053660,0.013400 /)
XEXT_COEFF_550_LKT(22,5)=707.440000 !rg=0.0559219 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(22,6,1:6)=(/ 6596.300000,2866.200000,984.380000,221.680000,26.627000,2.845900 /)
XPIZA_LKT(22,6,1:6)=(/ 0.983739,0.979980,0.968504,0.970980,0.889467,0.518965 /)
XCGA_LKT(22,6,1:6)=(/ 0.714263,0.618670,0.464443,0.253903,0.083783,0.021003 /)
XEXT_COEFF_550_LKT(22,6)=1009.700000 !rg=0.0559219 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(22,7,1:6)=(/ 6523.700000,3382.700000,1339.800000,344.010000,44.609000,4.067900 /)
XPIZA_LKT(22,7,1:6)=(/ 0.983045,0.982024,0.974969,0.979963,0.931786,0.659801 /)
XCGA_LKT(22,7,1:6)=(/ 0.726563,0.665800,0.549310,0.358480,0.125487,0.031500 /)
XEXT_COEFF_550_LKT(22,7)=1370.800000 !rg=0.0559219 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(22,8,1:6)=(/ 6094.000000,3739.300000,1717.000000,507.840000,72.938000,6.214900 /)
XPIZA_LKT(22,8,1:6)=(/ 0.981389,0.982920,0.978992,0.985364,0.956523,0.773901 /)
XCGA_LKT(22,8,1:6)=(/ 0.730027,0.696893,0.613230,0.452090,0.184210,0.046167 /)
XEXT_COEFF_550_LKT(22,8)=1752.100000 !rg=0.0559219 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(22,9,1:6)=(/ 5451.000000,3863.700000,2069.600000,710.090000,114.820000,9.983500 /)
XPIZA_LKT(22,9,1:6)=(/ 0.978842,0.982773,0.981406,0.988670,0.970941,0.856224 /)
XCGA_LKT(22,9,1:6)=(/ 0.729747,0.714640,0.660470,0.533623,0.267580,0.067077 /)
XEXT_COEFF_550_LKT(22,9)=2105.800000 !rg=0.0559219 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(22,10,1:6)=(/ 4743.900000,3733.700000,2332.300000,936.900000,174.480000,16.428000 /)
XPIZA_LKT(22,10,1:6)=(/ 0.975412,0.981582,0.982575,0.990750,0.979601,0.910110 /)
XCGA_LKT(22,10,1:6)=(/ 0.728967,0.721350,0.692647,0.599563,0.373290,0.097140 /)
XEXT_COEFF_550_LKT(22,10)=2365.400000 !rg=0.0559219 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(22,11,1:6)=(/ 4015.600000,3423.400000,2451.400000,1166.800000,258.500000,27.165000 /)
XPIZA_LKT(22,11,1:6)=(/ 0.970960,0.979420,0.982661,0.992042,0.985185,0.943615 /)
XCGA_LKT(22,11,1:6)=(/ 0.728983,0.721910,0.711537,0.650697,0.469330,0.141033 /)
XEXT_COEFF_550_LKT(22,11)=2477.800000 !rg=0.0559219 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(22,12,1:6)=(/ 3343.100000,3048.300000,2400.300000,1361.500000,360.320000,43.841000 /)
XPIZA_LKT(22,12,1:6)=(/ 0.965425,0.976500,0.981628,0.992754,0.988597,0.963469 /)
XCGA_LKT(22,12,1:6)=(/ 0.729443,0.723390,0.718577,0.687463,0.545793,0.205663 /)
XEXT_COEFF_550_LKT(22,12)=2419.100000 !rg=0.0559219 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(22,13,1:6)=(/ 2753.500000,2645.400000,2219.200000,1480.700000,478.560000,67.734000 /)
XPIZA_LKT(22,13,1:6)=(/ 0.959505,0.972694,0.979530,0.992986,0.990785,0.975026 /)
XCGA_LKT(22,13,1:6)=(/ 0.733413,0.725297,0.718493,0.710687,0.614047,0.299240 /)
XEXT_COEFF_550_LKT(22,13)=2227.800000 !rg=0.0559219 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(22,14,1:6)=(/ 2267.900000,2233.500000,1987.500000,1495.100000,598.110000,101.580000 /)
XPIZA_LKT(22,14,1:6)=(/ 0.951574,0.968034,0.976678,0.992745,0.992142,0.982115 /)
XCGA_LKT(22,14,1:6)=(/ 0.737963,0.726723,0.718167,0.720870,0.665050,0.413757 /)
XEXT_COEFF_550_LKT(22,14)=1997.500000 !rg=0.0559219 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(22,15,1:6)=(/ 1851.600000,1879.700000,1750.000000,1411.200000,703.530000,147.980000 /)
XPIZA_LKT(22,15,1:6)=(/ 0.942578,0.961668,0.973241,0.992022,0.992933,0.986775 /)
XCGA_LKT(22,15,1:6)=(/ 0.742860,0.730713,0.722533,0.720717,0.703173,0.504577 /)
XEXT_COEFF_550_LKT(22,15)=1754.700000 !rg=0.0559219 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(22,16,1:6)=(/ 1504.000000,1546.900000,1498.700000,1272.500000,775.150000,201.630000 /)
XPIZA_LKT(22,16,1:6)=(/ 0.933460,0.955292,0.968124,0.990921,0.993271,0.989598 /)
XCGA_LKT(22,16,1:6)=(/ 0.754717,0.734583,0.724787,0.717290,0.728307,0.576053 /)
XEXT_COEFF_550_LKT(22,16)=1495.600000 !rg=0.0559219 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(22,17,1:6)=(/ 1230.700000,1273.000000,1260.500000,1133.400000,793.850000,261.100000 /)
XPIZA_LKT(22,17,1:6)=(/ 0.922283,0.948332,0.963581,0.989551,0.993160,0.991416 /)
XCGA_LKT(22,17,1:6)=(/ 0.763550,0.741427,0.728940,0.718783,0.740447,0.637833 /)
XEXT_COEFF_550_LKT(22,17)=1266.900000 !rg=0.0559219 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(22,18,1:6)=(/ 1005.500000,1050.700000,1053.000000,993.330000,757.690000,318.560000 /)
XPIZA_LKT(22,18,1:6)=(/ 0.909242,0.938554,0.957271,0.987761,0.992593,0.992532 /)
XCGA_LKT(22,18,1:6)=(/ 0.771373,0.749483,0.732703,0.722330,0.741697,0.683460 /)
XEXT_COEFF_550_LKT(22,18)=1052.700000 !rg=0.0559219 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(22,19,1:6)=(/ 817.490000,856.020000,876.340000,840.270000,686.640000,365.600000 /)
XPIZA_LKT(22,19,1:6)=(/ 0.895165,0.928405,0.949190,0.985834,0.991575,0.993143 /)
XCGA_LKT(22,19,1:6)=(/ 0.785733,0.757390,0.739920,0.724753,0.736317,0.716207 /)
XEXT_COEFF_550_LKT(22,19)=872.020000 !rg=0.0559219 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(22,20,1:6)=(/ 670.540000,698.860000,718.290000,711.250000,609.750000,392.020000 /)
XPIZA_LKT(22,20,1:6)=(/ 0.879623,0.917837,0.942004,0.983534,0.990361,0.993317 /)
XCGA_LKT(22,20,1:6)=(/ 0.796407,0.767217,0.747447,0.729410,0.734550,0.736150 /)
XEXT_COEFF_550_LKT(22,20)=721.820000 !rg=0.0559219 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(23,1,1:6)=(/ 4322.200000,1066.200000,230.440000,37.007000,5.635500,1.519400 /)
XPIZA_LKT(23,1,1:6)=(/ 0.977847,0.958183,0.892630,0.849777,0.508810,0.114334 /)
XCGA_LKT(23,1,1:6)=(/ 0.595207,0.221453,0.090960,0.035927,0.011433,0.002820 /)
XEXT_COEFF_550_LKT(23,1)=234.860000 !rg=0.0605833 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(23,2,1:6)=(/ 4845.900000,1224.700000,279.930000,45.189000,6.392100,1.565900 /)
XPIZA_LKT(23,2,1:6)=(/ 0.979899,0.962672,0.909818,0.875598,0.565292,0.139835 /)
XCGA_LKT(23,2,1:6)=(/ 0.608610,0.281643,0.115400,0.045553,0.014520,0.003587 /)
XEXT_COEFF_550_LKT(23,2)=286.490000 !rg=0.0605833 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(23,3,1:6)=(/ 5487.400000,1539.000000,385.320000,64.726000,8.231100,1.679200 /)
XPIZA_LKT(23,3,1:6)=(/ 0.981886,0.968569,0.932041,0.911186,0.659693,0.196174 /)
XCGA_LKT(23,3,1:6)=(/ 0.639110,0.393017,0.173857,0.069043,0.022100,0.005470 /)
XEXT_COEFF_550_LKT(23,3)=395.670000 !rg=0.0605833 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(23,4,1:6)=(/ 6136.500000,2012.900000,557.690000,103.130000,12.086000,1.918000 /)
XPIZA_LKT(23,4,1:6)=(/ 0.983296,0.974183,0.950432,0.942312,0.765083,0.293668 /)
XCGA_LKT(23,4,1:6)=(/ 0.676133,0.501157,0.276887,0.116077,0.037577,0.009350 /)
XEXT_COEFF_550_LKT(23,4)=573.130000 !rg=0.0605833 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(23,5,1:6)=(/ 6600.700000,2574.700000,808.710000,168.460000,19.559000,2.388600 /)
XPIZA_LKT(23,5,1:6)=(/ 0.983955,0.978384,0.963280,0.963007,0.851813,0.429530 /)
XCGA_LKT(23,5,1:6)=(/ 0.705383,0.582200,0.396220,0.191453,0.062793,0.015713 /)
XEXT_COEFF_550_LKT(23,5)=830.550000 !rg=0.0605833 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(23,6,1:6)=(/ 6723.100000,3128.600000,1130.800000,266.590000,32.742000,3.251000 /)
XPIZA_LKT(23,6,1:6)=(/ 0.983763,0.981097,0.971591,0.975187,0.908940,0.577221 /)
XCGA_LKT(23,6,1:6)=(/ 0.723320,0.640137,0.497883,0.290613,0.098007,0.024617 /)
XEXT_COEFF_550_LKT(23,6)=1159.000000 !rg=0.0605833 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(23,7,1:6)=(/ 6459.700000,3587.500000,1500.900000,406.770000,54.670000,4.799200 /)
XPIZA_LKT(23,7,1:6)=(/ 0.982660,0.982615,0.976892,0.982514,0.943425,0.709990 /)
XCGA_LKT(23,7,1:6)=(/ 0.731183,0.680717,0.575563,0.394757,0.146827,0.036907 /)
XEXT_COEFF_550_LKT(23,7)=1534.500000 !rg=0.0605833 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(23,8,1:6)=(/ 5897.500000,3849.800000,1874.000000,586.110000,88.186000,7.512000 /)
XPIZA_LKT(23,8,1:6)=(/ 0.980595,0.983068,0.980180,0.986899,0.963319,0.811442 /)
XCGA_LKT(23,8,1:6)=(/ 0.731840,0.706217,0.633363,0.484100,0.215443,0.054080 /)
XEXT_COEFF_550_LKT(23,8)=1910.200000 !rg=0.0605833 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(23,9,1:6)=(/ 5200.500000,3862.600000,2198.600000,803.090000,136.790000,12.244000 /)
XPIZA_LKT(23,9,1:6)=(/ 0.977584,0.982497,0.982062,0.989651,0.974985,0.881492 /)
XCGA_LKT(23,9,1:6)=(/ 0.730287,0.719270,0.674747,0.561863,0.309910,0.078573 /)
XEXT_COEFF_550_LKT(23,9)=2234.100000 !rg=0.0605833 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(23,10,1:6)=(/ 4456.700000,3644.000000,2407.400000,1033.200000,206.230000,20.232000 /)
XPIZA_LKT(23,10,1:6)=(/ 0.973817,0.980875,0.982769,0.991351,0.982207,0.925977 /)
XCGA_LKT(23,10,1:6)=(/ 0.729723,0.722630,0.701513,0.621080,0.415450,0.113857 /)
XEXT_COEFF_550_LKT(23,10)=2438.400000 !rg=0.0605833 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(23,11,1:6)=(/ 3729.200000,3291.600000,2457.900000,1253.000000,298.100000,33.218000 /)
XPIZA_LKT(23,11,1:6)=(/ 0.968896,0.978422,0.982393,0.992386,0.986785,0.953082 /)
XCGA_LKT(23,11,1:6)=(/ 0.728637,0.723433,0.715557,0.666687,0.500563,0.165530 /)
XEXT_COEFF_550_LKT(23,11)=2481.100000 !rg=0.0605833 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(23,12,1:6)=(/ 3090.400000,2899.000000,2345.800000,1422.000000,407.450000,52.684000 /)
XPIZA_LKT(23,12,1:6)=(/ 0.963239,0.975077,0.980923,0.992898,0.989603,0.968957 /)
XCGA_LKT(23,12,1:6)=(/ 0.731363,0.724427,0.719350,0.698037,0.575623,0.241480 /)
XEXT_COEFF_550_LKT(23,12)=2359.000000 !rg=0.0605833 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(23,13,1:6)=(/ 2540.500000,2478.800000,2133.100000,1501.800000,527.470000,80.173000 /)
XPIZA_LKT(23,13,1:6)=(/ 0.955935,0.970949,0.978488,0.992941,0.991414,0.978313 /)
XCGA_LKT(23,13,1:6)=(/ 0.734260,0.726323,0.718300,0.716087,0.635760,0.346943 /)
XEXT_COEFF_550_LKT(23,13)=2144.100000 !rg=0.0605833 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(23,14,1:6)=(/ 2077.700000,2088.100000,1897.600000,1472.600000,643.820000,119.560000 /)
XPIZA_LKT(23,14,1:6)=(/ 0.947593,0.965260,0.975397,0.992511,0.992529,0.984307 /)
XCGA_LKT(23,14,1:6)=(/ 0.738817,0.728117,0.720553,0.721793,0.682003,0.456147 /)
XEXT_COEFF_550_LKT(23,14)=1903.900000 !rg=0.0605833 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(23,15,1:6)=(/ 1698.300000,1735.400000,1653.400000,1358.700000,737.620000,169.090000 /)
XPIZA_LKT(23,15,1:6)=(/ 0.938243,0.959313,0.971334,0.991635,0.993126,0.988105 /)
XCGA_LKT(23,15,1:6)=(/ 0.748907,0.731383,0.723640,0.719523,0.714487,0.533823 /)
XEXT_COEFF_550_LKT(23,15)=1653.200000 !rg=0.0605833 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(23,16,1:6)=(/ 1377.500000,1425.100000,1394.300000,1216.900000,789.660000,226.290000 /)
XPIZA_LKT(23,16,1:6)=(/ 0.928918,0.952964,0.966531,0.990433,0.993277,0.990452 /)
XCGA_LKT(23,16,1:6)=(/ 0.758243,0.738397,0.725827,0.718113,0.734370,0.604483 /)
XEXT_COEFF_550_LKT(23,16)=1401.000000 !rg=0.0605833 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(23,17,1:6)=(/ 1126.800000,1175.800000,1172.400000,1079.700000,785.100000,285.350000 /)
XPIZA_LKT(23,17,1:6)=(/ 0.916291,0.943869,0.961364,0.988937,0.992982,0.991937 /)
XCGA_LKT(23,17,1:6)=(/ 0.765860,0.744627,0.731120,0.720443,0.742007,0.657780 /)
XEXT_COEFF_550_LKT(23,17)=1174.300000 !rg=0.0605833 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(23,18,1:6)=(/ 921.100000,961.090000,977.600000,929.190000,731.070000,339.550000 /)
XPIZA_LKT(23,18,1:6)=(/ 0.902798,0.934581,0.953763,0.986997,0.992212,0.992833 /)
XCGA_LKT(23,18,1:6)=(/ 0.778717,0.750960,0.735760,0.723283,0.739483,0.698490 /)
XEXT_COEFF_550_LKT(23,18)=972.220000 !rg=0.0605833 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(23,19,1:6)=(/ 749.960000,782.430000,804.590000,784.700000,654.110000,379.370000 /)
XPIZA_LKT(23,19,1:6)=(/ 0.888167,0.924713,0.946357,0.984920,0.991094,0.993261 /)
XCGA_LKT(23,19,1:6)=(/ 0.790663,0.762520,0.742230,0.727240,0.734920,0.725783 /)
XEXT_COEFF_550_LKT(23,19)=807.100000 !rg=0.0605833 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(23,20,1:6)=(/ 614.740000,643.640000,660.830000,663.040000,581.550000,395.180000 /)
XPIZA_LKT(23,20,1:6)=(/ 0.872438,0.911646,0.938380,0.981909,0.989808,0.993251 /)
XCGA_LKT(23,20,1:6)=(/ 0.799347,0.772307,0.750250,0.731113,0.736353,0.740580 /)
XEXT_COEFF_550_LKT(23,20)=661.020000 !rg=0.0605833 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(24,1,1:6)=(/ 5099.200000,1274.300000,284.660000,45.588000,6.426000,1.568100 /)
XPIZA_LKT(24,1,1:6)=(/ 0.980017,0.963868,0.910783,0.876397,0.567270,0.140856 /)
XCGA_LKT(24,1,1:6)=(/ 0.633857,0.264210,0.106823,0.042113,0.013410,0.003310 /)
XEXT_COEFF_550_LKT(24,1)=291.630000 !rg=0.0656333 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(24,2,1:6)=(/ 5481.200000,1443.500000,344.530000,55.925000,7.386300,1.627200 /)
XPIZA_LKT(24,2,1:6)=(/ 0.981631,0.967142,0.924711,0.897965,0.621852,0.171104 /)
XCGA_LKT(24,2,1:6)=(/ 0.634597,0.333860,0.135633,0.053387,0.017027,0.004207 /)
XEXT_COEFF_550_LKT(24,2)=353.870000 !rg=0.0656333 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(24,3,1:6)=(/ 6024.100000,1795.600000,467.530000,80.388000,9.718500,1.770900 /)
XPIZA_LKT(24,3,1:6)=(/ 0.983044,0.971878,0.942344,0.927224,0.709915,0.236506 /)
XCGA_LKT(24,3,1:6)=(/ 0.662863,0.444830,0.203750,0.080847,0.025903,0.006420 /)
XEXT_COEFF_550_LKT(24,3)=480.780000 !rg=0.0656333 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(24,4,1:6)=(/ 6555.000000,2300.000000,663.690000,127.300000,14.593000,2.073900 /)
XPIZA_LKT(24,4,1:6)=(/ 0.983993,0.976485,0.956949,0.952288,0.803789,0.345271 /)
XCGA_LKT(24,4,1:6)=(/ 0.694670,0.539403,0.316027,0.135470,0.044000,0.010967 /)
XEXT_COEFF_550_LKT(24,4)=682.300000 !rg=0.0656333 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(24,5,1:6)=(/ 6846.600000,2861.500000,942.910000,204.730000,23.971000,2.671100 /)
XPIZA_LKT(24,5,1:6)=(/ 0.984228,0.979876,0.967346,0.968795,0.877727,0.488190 /)
XCGA_LKT(24,5,1:6)=(/ 0.717970,0.609300,0.434407,0.221690,0.073450,0.018420 /)
XEXT_COEFF_550_LKT(24,5)=968.060000 !rg=0.0656333 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(24,6,1:6)=(/ 6769.200000,3378.300000,1287.000000,318.920000,40.269000,3.764600 /)
XPIZA_LKT(24,6,1:6)=(/ 0.983619,0.981992,0.974144,0.978633,0.924886,0.633194 /)
XCGA_LKT(24,6,1:6)=(/ 0.730650,0.659293,0.529020,0.327780,0.114600,0.028847 /)
XEXT_COEFF_550_LKT(24,6)=1318.100000 !rg=0.0656333 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(24,7,1:6)=(/ 6332.600000,3761.200000,1666.400000,476.190000,66.698000,5.724300 /)
XPIZA_LKT(24,7,1:6)=(/ 0.982076,0.983031,0.978501,0.984581,0.952796,0.755252 /)
XCGA_LKT(24,7,1:6)=(/ 0.734520,0.693507,0.599977,0.428887,0.171700,0.043237 /)
XEXT_COEFF_550_LKT(24,7)=1701.800000 !rg=0.0656333 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(24,8,1:6)=(/ 5665.000000,3918.800000,2026.000000,671.920000,105.940000,9.144500 /)
XPIZA_LKT(24,8,1:6)=(/ 0.979626,0.983045,0.981130,0.988177,0.968790,0.843691 /)
XCGA_LKT(24,8,1:6)=(/ 0.733470,0.713687,0.651363,0.515810,0.251230,0.063340 /)
XEXT_COEFF_550_LKT(24,8)=2063.000000 !rg=0.0656333 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(24,9,1:6)=(/ 4912.000000,3824.600000,2311.100000,897.800000,162.520000,15.060000 /)
XPIZA_LKT(24,9,1:6)=(/ 0.976276,0.982034,0.982507,0.990452,0.978348,0.902473 /)
XCGA_LKT(24,9,1:6)=(/ 0.731593,0.722487,0.687127,0.586560,0.353923,0.092047 /)
XEXT_COEFF_550_LKT(24,9)=2345.600000 !rg=0.0656333 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(24,10,1:6)=(/ 4164.800000,3530.300000,2457.300000,1128.200000,241.770000,24.873000 /)
XPIZA_LKT(24,10,1:6)=(/ 0.971753,0.980086,0.982781,0.991852,0.984365,0.938853 /)
XCGA_LKT(24,10,1:6)=(/ 0.729050,0.724303,0.708753,0.640753,0.452040,0.133497 /)
XEXT_COEFF_550_LKT(24,10)=2485.200000 !rg=0.0656333 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(24,11,1:6)=(/ 3454.700000,3151.700000,2437.000000,1332.200000,340.290000,40.352000 /)
XPIZA_LKT(24,11,1:6)=(/ 0.966639,0.977206,0.981972,0.992656,0.988081,0.960649 /)
XCGA_LKT(24,11,1:6)=(/ 0.729487,0.724187,0.718380,0.680567,0.531043,0.194347 /)
XEXT_COEFF_550_LKT(24,11)=2456.600000 !rg=0.0656333 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(24,12,1:6)=(/ 2855.300000,2738.700000,2274.400000,1468.200000,456.390000,62.820000 /)
XPIZA_LKT(24,12,1:6)=(/ 0.960027,0.973564,0.980076,0.992976,0.990450,0.973356 /)
XCGA_LKT(24,12,1:6)=(/ 0.731887,0.726243,0.719310,0.706643,0.602150,0.282753 /)
XEXT_COEFF_550_LKT(24,12)=2288.300000 !rg=0.0656333 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(24,13,1:6)=(/ 2331.100000,2318.200000,2044.100000,1503.900000,577.010000,94.758000 /)
XPIZA_LKT(24,13,1:6)=(/ 0.952399,0.968580,0.977428,0.992829,0.991937,0.981080 /)
XCGA_LKT(24,13,1:6)=(/ 0.735287,0.726677,0.719833,0.719560,0.655913,0.395340 /)
XEXT_COEFF_550_LKT(24,13)=2052.400000 !rg=0.0656333 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(24,14,1:6)=(/ 1908.800000,1937.600000,1806.100000,1435.300000,686.490000,139.280000 /)
XPIZA_LKT(24,14,1:6)=(/ 0.943238,0.962838,0.973982,0.992217,0.992821,0.986121 /)
XCGA_LKT(24,14,1:6)=(/ 0.743580,0.728973,0.722190,0.721583,0.696657,0.490910 /)
XEXT_COEFF_550_LKT(24,14)=1809.400000 !rg=0.0656333 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(24,15,1:6)=(/ 1549.200000,1597.400000,1546.300000,1305.000000,765.350000,191.410000 /)
XPIZA_LKT(24,15,1:6)=(/ 0.934352,0.957280,0.969491,0.991180,0.993233,0.989184 /)
XCGA_LKT(24,15,1:6)=(/ 0.752087,0.735323,0.724673,0.718833,0.724047,0.563360 /)
XEXT_COEFF_550_LKT(24,15)=1550.900000 !rg=0.0656333 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(24,16,1:6)=(/ 1260.700000,1317.200000,1299.000000,1163.800000,794.760000,250.970000 /)
XPIZA_LKT(24,16,1:6)=(/ 0.923659,0.949043,0.964846,0.989844,0.993209,0.991161 /)
XCGA_LKT(24,16,1:6)=(/ 0.761000,0.740537,0.728997,0.718443,0.738780,0.628550 /)
XEXT_COEFF_550_LKT(24,16)=1303.700000 !rg=0.0656333 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(24,17,1:6)=(/ 1033.600000,1079.000000,1090.100000,1018.900000,767.840000,309.330000 /)
XPIZA_LKT(24,17,1:6)=(/ 0.910529,0.939936,0.958152,0.988295,0.992716,0.992376 /)
XCGA_LKT(24,17,1:6)=(/ 0.772437,0.745900,0.732970,0.723210,0.741703,0.676250 /)
XEXT_COEFF_550_LKT(24,17)=1084.800000 !rg=0.0656333 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(24,18,1:6)=(/ 841.200000,879.380000,901.400000,866.800000,699.960000,358.600000 /)
XPIZA_LKT(24,18,1:6)=(/ 0.896381,0.930947,0.951108,0.986385,0.991808,0.993061 /)
XCGA_LKT(24,18,1:6)=(/ 0.783383,0.757827,0.737393,0.726070,0.737747,0.711223 /)
XEXT_COEFF_550_LKT(24,18)=902.720000 !rg=0.0656333 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(24,19,1:6)=(/ 687.600000,720.940000,738.810000,735.140000,624.280000,389.180000 /)
XPIZA_LKT(24,19,1:6)=(/ 0.881044,0.918995,0.943672,0.983578,0.990579,0.993308 /)
XCGA_LKT(24,19,1:6)=(/ 0.793987,0.766537,0.746343,0.728897,0.734980,0.733397 /)
XEXT_COEFF_550_LKT(24,19)=740.980000 !rg=0.0656333 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(24,20,1:6)=(/ 564.260000,588.770000,610.230000,612.750000,553.860000,393.020000 /)
XPIZA_LKT(24,20,1:6)=(/ 0.866299,0.905690,0.933329,0.980655,0.989114,0.993111 /)
XCGA_LKT(24,20,1:6)=(/ 0.805837,0.774373,0.754227,0.730143,0.737357,0.742897 /)
XEXT_COEFF_550_LKT(24,20)=605.350000 !rg=0.0656333 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(25,1,1:6)=(/ 5876.000000,1500.400000,351.590000,56.457000,7.429400,1.629900 /)
XPIZA_LKT(25,1,1:6)=(/ 0.982226,0.968247,0.925708,0.898647,0.623726,0.172302 /)
XCGA_LKT(25,1,1:6)=(/ 0.644623,0.316663,0.125620,0.049363,0.015730,0.003883 /)
XEXT_COEFF_550_LKT(25,1)=361.550000 !rg=0.0711042 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(25,2,1:6)=(/ 6068.000000,1686.200000,422.610000,69.482000,8.647900,1.704900 /)
XPIZA_LKT(25,2,1:6)=(/ 0.983062,0.970664,0.936831,0.916472,0.675074,0.207637 /)
XCGA_LKT(25,2,1:6)=(/ 0.654257,0.393643,0.159667,0.062570,0.019967,0.004937 /)
XEXT_COEFF_550_LKT(25,2)=435.040000 !rg=0.0711042 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(25,3,1:6)=(/ 6520.500000,2081.300000,563.450000,99.930000,11.604000,1.887300 /)
XPIZA_LKT(25,3,1:6)=(/ 0.983968,0.974624,0.950660,0.940301,0.755238,0.282156 /)
XCGA_LKT(25,3,1:6)=(/ 0.684047,0.494920,0.238593,0.094670,0.030360,0.007530 /)
XEXT_COEFF_550_LKT(25,3)=579.790000 !rg=0.0711042 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(25,4,1:6)=(/ 6911.200000,2603.200000,784.540000,156.630000,17.757000,2.271800 /)
XPIZA_LKT(25,4,1:6)=(/ 0.984501,0.978407,0.962271,0.960333,0.837183,0.400690 /)
XCGA_LKT(25,4,1:6)=(/ 0.710883,0.574097,0.357687,0.157930,0.051507,0.012860 /)
XEXT_COEFF_550_LKT(25,4)=806.540000 !rg=0.0711042 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(25,5,1:6)=(/ 7010.000000,3146.900000,1090.800000,247.580000,29.472000,3.029400 /)
XPIZA_LKT(25,5,1:6)=(/ 0.984327,0.981110,0.970708,0.973491,0.899290,0.547005 /)
XCGA_LKT(25,5,1:6)=(/ 0.728497,0.633657,0.471550,0.254833,0.085877,0.021593 /)
XEXT_COEFF_550_LKT(25,5)=1119.200000 !rg=0.0711042 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(25,6,1:6)=(/ 6732.500000,3608.300000,1452.200000,378.680000,49.443000,4.415200 /)
XPIZA_LKT(25,6,1:6)=(/ 0.983290,0.982693,0.976265,0.981438,0.937844,0.685556 /)
XCGA_LKT(25,6,1:6)=(/ 0.736177,0.676283,0.558170,0.364153,0.133930,0.033800 /)
XEXT_COEFF_550_LKT(25,6)=1486.100000 !rg=0.0711042 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(25,7,1:6)=(/ 6141.700000,3898.500000,1830.700000,552.600000,80.914000,6.892400 /)
XPIZA_LKT(25,7,1:6)=(/ 0.981316,0.983279,0.979819,0.986266,0.960326,0.795192 /)
XCGA_LKT(25,7,1:6)=(/ 0.736567,0.704473,0.622027,0.462330,0.200533,0.050643 /)
XEXT_COEFF_550_LKT(25,7)=1867.600000 !rg=0.0711042 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(25,8,1:6)=(/ 5400.700000,3943.200000,2166.100000,763.710000,126.620000,11.192000 /)
XPIZA_LKT(25,8,1:6)=(/ 0.978405,0.982871,0.981888,0.989245,0.973245,0.870962 /)
XCGA_LKT(25,8,1:6)=(/ 0.733710,0.719690,0.667197,0.545293,0.291040,0.074183 /)
XEXT_COEFF_550_LKT(25,8)=2202.700000 !rg=0.0711042 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(25,9,1:6)=(/ 4615.200000,3748.300000,2400.600000,994.900000,192.500000,18.541000 /)
XPIZA_LKT(25,9,1:6)=(/ 0.974630,0.981461,0.982796,0.991109,0.981169,0.919709 /)
XCGA_LKT(25,9,1:6)=(/ 0.731717,0.724963,0.697387,0.609233,0.395943,0.107847 /)
XEXT_COEFF_550_LKT(25,9)=2433.300000 !rg=0.0711042 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(25,10,1:6)=(/ 3859.400000,3405.100000,2479.000000,1217.900000,280.040000,30.460000 /)
XPIZA_LKT(25,10,1:6)=(/ 0.969573,0.979065,0.982621,0.992241,0.986110,0.949225 /)
XCGA_LKT(25,10,1:6)=(/ 0.728993,0.725160,0.714003,0.657907,0.484330,0.156590 /)
XEXT_COEFF_550_LKT(25,10)=2504.600000 !rg=0.0711042 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(25,11,1:6)=(/ 3203.600000,2997.300000,2394.000000,1399.200000,386.090000,48.637000 /)
XPIZA_LKT(25,11,1:6)=(/ 0.963584,0.975978,0.981367,0.992842,0.989168,0.966684 /)
XCGA_LKT(25,11,1:6)=(/ 0.730313,0.726450,0.719763,0.692200,0.561467,0.228100 /)
XEXT_COEFF_550_LKT(25,11)=2411.200000 !rg=0.0711042 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(25,12,1:6)=(/ 2617.300000,2574.200000,2192.700000,1497.700000,505.230000,74.500000 /)
XPIZA_LKT(25,12,1:6)=(/ 0.956682,0.971610,0.979167,0.992972,0.991139,0.976938 /)
XCGA_LKT(25,12,1:6)=(/ 0.732037,0.726633,0.720213,0.712993,0.624913,0.328680 /)
XEXT_COEFF_550_LKT(25,12)=2203.700000 !rg=0.0711042 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(25,13,1:6)=(/ 2144.600000,2153.900000,1956.900000,1489.300000,624.270000,111.760000 /)
XPIZA_LKT(25,13,1:6)=(/ 0.948187,0.966235,0.976109,0.992640,0.992368,0.983436 /)
XCGA_LKT(25,13,1:6)=(/ 0.739297,0.727107,0.720823,0.721413,0.673890,0.439417 /)
XEXT_COEFF_550_LKT(25,13)=1962.000000 !rg=0.0711042 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(25,14,1:6)=(/ 1741.400000,1792.000000,1702.700000,1390.100000,723.680000,159.920000 /)
XPIZA_LKT(25,14,1:6)=(/ 0.939640,0.960835,0.972472,0.991835,0.993049,0.987574 /)
XCGA_LKT(25,14,1:6)=(/ 0.746423,0.733047,0.724563,0.720710,0.708953,0.520953 /)
XEXT_COEFF_550_LKT(25,14)=1710.000000 !rg=0.0711042 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(25,15,1:6)=(/ 1416.100000,1475.000000,1442.200000,1250.000000,783.860000,215.520000 /)
XPIZA_LKT(25,15,1:6)=(/ 0.930371,0.953925,0.967916,0.990664,0.993276,0.990098 /)
XCGA_LKT(25,15,1:6)=(/ 0.755960,0.736780,0.727103,0.717993,0.731083,0.592657 /)
XEXT_COEFF_550_LKT(25,15)=1445.300000 !rg=0.0711042 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(25,16,1:6)=(/ 1159.000000,1208.400000,1213.300000,1106.900000,790.480000,275.190000 /)
XPIZA_LKT(25,16,1:6)=(/ 0.917663,0.944946,0.961998,0.989300,0.993060,0.991728 /)
XCGA_LKT(25,16,1:6)=(/ 0.766560,0.741530,0.730487,0.721493,0.741020,0.649197 /)
XEXT_COEFF_550_LKT(25,16)=1211.400000 !rg=0.0711042 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(25,17,1:6)=(/ 942.040000,989.680000,1005.600000,958.250000,743.340000,331.110000 /)
XPIZA_LKT(25,17,1:6)=(/ 0.904824,0.935508,0.955029,0.987429,0.992395,0.992717 /)
XCGA_LKT(25,17,1:6)=(/ 0.775883,0.751870,0.733113,0.724603,0.740710,0.692193 /)
XEXT_COEFF_550_LKT(25,17)=1005.200000 !rg=0.0711042 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(25,18,1:6)=(/ 768.980000,807.180000,828.710000,810.340000,669.420000,373.760000 /)
XPIZA_LKT(25,18,1:6)=(/ 0.890167,0.925983,0.948694,0.985278,0.991327,0.993215 /)
XCGA_LKT(25,18,1:6)=(/ 0.788423,0.761103,0.742903,0.726203,0.736087,0.721523 /)
XEXT_COEFF_550_LKT(25,18)=830.930000 !rg=0.0711042 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(25,19,1:6)=(/ 631.880000,660.380000,683.370000,682.260000,595.300000,394.080000 /)
XPIZA_LKT(25,19,1:6)=(/ 0.873931,0.913067,0.939190,0.982327,0.990042,0.993275 /)
XCGA_LKT(25,19,1:6)=(/ 0.800207,0.768723,0.749433,0.728833,0.735697,0.738450 /)
XEXT_COEFF_550_LKT(25,19)=679.940000 !rg=0.0711042 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(25,20,1:6)=(/ 515.350000,539.360000,559.570000,567.010000,523.110000,386.150000 /)
XPIZA_LKT(25,20,1:6)=(/ 0.860698,0.899130,0.928922,0.979345,0.988534,0.992890 /)
XCGA_LKT(25,20,1:6)=(/ 0.809577,0.781777,0.755730,0.736033,0.739860,0.743393 /)
XEXT_COEFF_550_LKT(25,20)=558.740000 !rg=0.0711042 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(26,1,1:6)=(/ 6421.700000,1741.100000,433.160000,70.204000,8.702900,1.708300 /)
XPIZA_LKT(26,1,1:6)=(/ 0.983966,0.971541,0.937883,0.917063,0.676807,0.209022 /)
XCGA_LKT(26,1,1:6)=(/ 0.647320,0.380350,0.148010,0.057867,0.018450,0.004557 /)
XEXT_COEFF_550_LKT(26,1)=446.520000 !rg=0.0770312 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(26,2,1:6)=(/ 6587.100000,1959.900000,515.330000,86.553000,10.249000,1.803400 /)
XPIZA_LKT(26,2,1:6)=(/ 0.984145,0.973502,0.946615,0.931661,0.723914,0.249571 /)
XCGA_LKT(26,2,1:6)=(/ 0.674253,0.457803,0.188343,0.073343,0.023413,0.005793 /)
XEXT_COEFF_550_LKT(26,2)=531.040000 !rg=0.0770312 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(26,3,1:6)=(/ 6960.600000,2394.500000,673.860000,124.140000,13.993000,2.034900 /)
XPIZA_LKT(26,3,1:6)=(/ 0.984672,0.976933,0.957360,0.950900,0.795276,0.332700 /)
XCGA_LKT(26,3,1:6)=(/ 0.703377,0.540520,0.278730,0.110860,0.035577,0.008833 /)
XEXT_COEFF_550_LKT(26,3)=693.460000 !rg=0.0770312 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(26,4,1:6)=(/ 7189.300000,2916.000000,920.940000,191.860000,21.741000,2.522900 /)
XPIZA_LKT(26,4,1:6)=(/ 0.984826,0.980013,0.966633,0.966809,0.865543,0.458643 /)
XCGA_LKT(26,4,1:6)=(/ 0.724887,0.604997,0.400770,0.183777,0.060270,0.015080 /)
XEXT_COEFF_550_LKT(26,4)=946.490000 !rg=0.0770312 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(26,5,1:6)=(/ 7081.900000,3423.000000,1251.200000,297.650000,36.289000,3.483700 /)
XPIZA_LKT(26,5,1:6)=(/ 0.984248,0.982115,0.973492,0.977302,0.917049,0.604343 /)
XCGA_LKT(26,5,1:6)=(/ 0.737000,0.655370,0.506947,0.290007,0.100347,0.025307 /)
XEXT_COEFF_550_LKT(26,5)=1282.800000 !rg=0.0770312 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(26,6,1:6)=(/ 6612.700000,3810.000000,1623.600000,445.770000,60.506000,5.238600 /)
XPIZA_LKT(26,6,1:6)=(/ 0.982800,0.983212,0.978036,0.983717,0.948316,0.733325 /)
XCGA_LKT(26,6,1:6)=(/ 0.740373,0.691083,0.585243,0.399777,0.156393,0.039597 /)
XEXT_COEFF_550_LKT(26,6)=1659.600000 !rg=0.0770312 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(26,7,1:6)=(/ 5897.800000,3994.400000,1992.000000,636.830000,97.584000,8.364200 /)
XPIZA_LKT(26,7,1:6)=(/ 0.980393,0.983358,0.980888,0.987666,0.966395,0.829768 /)
XCGA_LKT(26,7,1:6)=(/ 0.738190,0.713567,0.641920,0.495347,0.233560,0.059310 /)
XEXT_COEFF_550_LKT(26,7)=2029.900000 !rg=0.0770312 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(26,8,1:6)=(/ 5097.200000,3926.700000,2291.600000,858.580000,150.840000,13.747000 /)
XPIZA_LKT(26,8,1:6)=(/ 0.977089,0.982516,0.982436,0.990122,0.976928,0.893723 /)
XCGA_LKT(26,8,1:6)=(/ 0.734553,0.724133,0.681087,0.571547,0.333050,0.086880 /)
XEXT_COEFF_550_LKT(26,8)=2327.600000 !rg=0.0770312 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(26,9,1:6)=(/ 4296.100000,3647.100000,2466.300000,1092.200000,226.380000,22.805000 /)
XPIZA_LKT(26,9,1:6)=(/ 0.972735,0.980711,0.982900,0.991658,0.983505,0.933741 /)
XCGA_LKT(26,9,1:6)=(/ 0.730890,0.726940,0.705907,0.630147,0.433377,0.126390 /)
XEXT_COEFF_550_LKT(26,9)=2496.200000 !rg=0.0770312 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(26,10,1:6)=(/ 3577.000000,3258.500000,2474.100000,1302.300000,320.960000,37.087000 /)
XPIZA_LKT(26,10,1:6)=(/ 0.967035,0.978060,0.982292,0.992549,0.987523,0.957534 /)
XCGA_LKT(26,10,1:6)=(/ 0.730053,0.726853,0.717893,0.672973,0.515340,0.183723 /)
XEXT_COEFF_550_LKT(26,10)=2495.600000 !rg=0.0770312 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(26,11,1:6)=(/ 2940.500000,2838.500000,2330.900000,1453.700000,434.440000,58.166000 /)
XPIZA_LKT(26,11,1:6)=(/ 0.960725,0.974446,0.980669,0.992954,0.990087,0.971514 /)
XCGA_LKT(26,11,1:6)=(/ 0.730197,0.727400,0.721070,0.702000,0.589250,0.267143 /)
XEXT_COEFF_550_LKT(26,11)=2346.100000 !rg=0.0770312 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(26,12,1:6)=(/ 2408.200000,2397.800000,2109.000000,1509.600000,554.870000,88.155000 /)
XPIZA_LKT(26,12,1:6)=(/ 0.952885,0.969553,0.978032,0.992898,0.991706,0.979929 /)
XCGA_LKT(26,12,1:6)=(/ 0.735637,0.726943,0.720577,0.717643,0.645853,0.376507 /)
XEXT_COEFF_550_LKT(26,12)=2116.000000 !rg=0.0770312 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(26,13,1:6)=(/ 1953.800000,2001.200000,1859.200000,1460.600000,668.440000,130.730000 /)
XPIZA_LKT(26,13,1:6)=(/ 0.944920,0.964221,0.974853,0.992373,0.992695,0.985401 /)
XCGA_LKT(26,13,1:6)=(/ 0.741337,0.730767,0.723380,0.721963,0.689487,0.476227 /)
XEXT_COEFF_550_LKT(26,13)=1868.200000 !rg=0.0770312 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(26,14,1:6)=(/ 1587.600000,1650.900000,1598.700000,1337.500000,754.480000,181.610000 /)
XPIZA_LKT(26,14,1:6)=(/ 0.936312,0.958079,0.970685,0.991403,0.993187,0.988743 /)
XCGA_LKT(26,14,1:6)=(/ 0.751053,0.733723,0.725993,0.719360,0.719333,0.550433 /)
XEXT_COEFF_550_LKT(26,14)=1600.500000 !rg=0.0770312 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(26,15,1:6)=(/ 1301.300000,1357.000000,1345.700000,1194.600000,793.750000,240.220000 /)
XPIZA_LKT(26,15,1:6)=(/ 0.924708,0.950353,0.965571,0.990164,0.993241,0.990869 /)
XCGA_LKT(26,15,1:6)=(/ 0.760553,0.738410,0.727717,0.719700,0.736470,0.618267 /)
XEXT_COEFF_550_LKT(26,15)=1349.400000 !rg=0.0770312 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(26,16,1:6)=(/ 1056.900000,1113.100000,1122.800000,1050.300000,776.730000,299.520000 /)
XPIZA_LKT(26,16,1:6)=(/ 0.912247,0.940713,0.959274,0.988641,0.992841,0.992200 /)
XCGA_LKT(26,16,1:6)=(/ 0.768693,0.746667,0.731303,0.723297,0.741787,0.668393 /)
XEXT_COEFF_550_LKT(26,16)=1121.700000 !rg=0.0770312 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(26,17,1:6)=(/ 858.790000,904.970000,929.860000,893.820000,715.520000,351.150000 /)
XPIZA_LKT(26,17,1:6)=(/ 0.899224,0.931633,0.951995,0.986695,0.991980,0.992969 /)
XCGA_LKT(26,17,1:6)=(/ 0.783290,0.755427,0.738327,0.725020,0.738580,0.705817 /)
XEXT_COEFF_550_LKT(26,17)=927.870000 !rg=0.0770312 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(26,18,1:6)=(/ 706.720000,740.490000,763.280000,755.900000,638.400000,385.470000 /)
XPIZA_LKT(26,18,1:6)=(/ 0.882832,0.920316,0.944792,0.984264,0.990826,0.993286 /)
XCGA_LKT(26,18,1:6)=(/ 0.794187,0.763723,0.744903,0.727453,0.735243,0.730077 /)
XEXT_COEFF_550_LKT(26,18)=765.510000 !rg=0.0770312 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(26,19,1:6)=(/ 577.420000,606.260000,626.870000,633.610000,566.610000,393.990000 /)
XPIZA_LKT(26,19,1:6)=(/ 0.868042,0.906913,0.934571,0.980541,0.989487,0.993168 /)
XCGA_LKT(26,19,1:6)=(/ 0.803210,0.775200,0.750787,0.731557,0.738290,0.741720 /)
XEXT_COEFF_550_LKT(26,19)=624.340000 !rg=0.0770312 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(26,20,1:6)=(/ 471.250000,493.350000,513.270000,522.760000,492.690000,375.400000 /)
XPIZA_LKT(26,20,1:6)=(/ 0.855020,0.892612,0.923897,0.978259,0.987433,0.992578 /)
XCGA_LKT(26,20,1:6)=(/ 0.816330,0.786527,0.762267,0.738617,0.740330,0.742127 /)
XEXT_COEFF_550_LKT(26,20)=512.930000 !rg=0.0770312 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(27,1,1:6)=(/ 6748.600000,2002.600000,530.900000,87.553000,10.319000,1.807700 /)
XPIZA_LKT(27,1,1:6)=(/ 0.984869,0.974005,0.947734,0.932185,0.725478,0.251140 /)
XCGA_LKT(27,1,1:6)=(/ 0.664543,0.454503,0.174863,0.067847,0.021633,0.005350 /)
XEXT_COEFF_550_LKT(27,1)=547.870000 !rg=0.0834522 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(27,2,1:6)=(/ 7053.200000,2274.600000,623.130000,107.960000,12.279000,1.928500 /)
XPIZA_LKT(27,2,1:6)=(/ 0.984915,0.975900,0.954449,0.944044,0.767708,0.296734 /)
XCGA_LKT(27,2,1:6)=(/ 0.697083,0.519850,0.222710,0.085997,0.027450,0.006797 /)
XEXT_COEFF_550_LKT(27,2)=642.140000 !rg=0.0834522 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(27,3,1:6)=(/ 7329.400000,2729.600000,799.590000,153.880000,17.015000,2.222300 /)
XPIZA_LKT(27,3,1:6)=(/ 0.985173,0.978890,0.962759,0.959447,0.829985,0.387311 /)
XCGA_LKT(27,3,1:6)=(/ 0.720877,0.579903,0.324013,0.129843,0.041680,0.010360 /)
XEXT_COEFF_550_LKT(27,3)=822.650000 !rg=0.0834522 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(27,4,1:6)=(/ 7375.200000,3230.000000,1073.100000,233.730000,26.739000,2.841500 /)
XPIZA_LKT(27,4,1:6)=(/ 0.984966,0.981349,0.970222,0.972019,0.889305,0.517579 /)
XCGA_LKT(27,4,1:6)=(/ 0.736663,0.632210,0.443920,0.213190,0.070493,0.017680 /)
XEXT_COEFF_550_LKT(27,4)=1102.300000 !rg=0.0834522 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(27,5,1:6)=(/ 7057.300000,3681.000000,1423.000000,355.350000,44.672000,4.059600 /)
XPIZA_LKT(27,5,1:6)=(/ 0.983983,0.982919,0.975808,0.980393,0.931552,0.658726 /)
XCGA_LKT(27,5,1:6)=(/ 0.743570,0.674570,0.540300,0.326513,0.117170,0.029653 /)
XEXT_COEFF_550_LKT(27,5)=1457.700000 !rg=0.0834522 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(27,6,1:6)=(/ 6421.900000,3976.300000,1797.000000,520.620000,73.707000,6.279200 /)
XPIZA_LKT(27,6,1:6)=(/ 0.982066,0.983562,0.979496,0.985580,0.956756,0.775923 /)
XCGA_LKT(27,6,1:6)=(/ 0.742787,0.703893,0.609897,0.435230,0.182373,0.046377 /)
XEXT_COEFF_550_LKT(27,6)=1834.900000 !rg=0.0834522 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(27,7,1:6)=(/ 5615.400000,4044.400000,2144.000000,727.690000,117.100000,10.212000 /)
XPIZA_LKT(27,7,1:6)=(/ 0.979180,0.983283,0.981753,0.988831,0.971327,0.859203 /)
XCGA_LKT(27,7,1:6)=(/ 0.738183,0.721010,0.659610,0.526367,0.270460,0.069447 /)
XEXT_COEFF_550_LKT(27,7)=2182.000000 !rg=0.0834522 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(27,8,1:6)=(/ 4781.700000,3867.900000,2396.400000,956.670000,179.090000,16.914000 /)
XPIZA_LKT(27,8,1:6)=(/ 0.975365,0.982040,0.982823,0.990844,0.979999,0.912499 /)
XCGA_LKT(27,8,1:6)=(/ 0.733830,0.727677,0.692903,0.595723,0.374173,0.101750 /)
XEXT_COEFF_550_LKT(27,8)=2430.900000 !rg=0.0834522 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(27,9,1:6)=(/ 3977.800000,3522.400000,2503.800000,1185.500000,263.280000,27.961000 /)
XPIZA_LKT(27,9,1:6)=(/ 0.970719,0.979789,0.982843,0.992095,0.985402,0.945075 /)
XCGA_LKT(27,9,1:6)=(/ 0.730593,0.728013,0.712523,0.648593,0.466737,0.148153 /)
XEXT_COEFF_550_LKT(27,9)=2531.300000 !rg=0.0834522 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(27,10,1:6)=(/ 3289.500000,3106.000000,2443.200000,1375.800000,365.460000,44.827000 /)
XPIZA_LKT(27,10,1:6)=(/ 0.964799,0.976816,0.981821,0.992777,0.988700,0.964170 /)
XCGA_LKT(27,10,1:6)=(/ 0.730147,0.728420,0.720550,0.685790,0.546313,0.215493 /)
XEXT_COEFF_550_LKT(27,10)=2463.100000 !rg=0.0834522 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(27,11,1:6)=(/ 2693.600000,2661.200000,2258.300000,1491.900000,483.200000,69.144000 /)
XPIZA_LKT(27,11,1:6)=(/ 0.958054,0.972755,0.979745,0.992991,0.990840,0.975430 /)
XCGA_LKT(27,11,1:6)=(/ 0.733447,0.728143,0.721557,0.709477,0.613213,0.311040 /)
XEXT_COEFF_550_LKT(27,11)=2268.400000 !rg=0.0834522 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(27,12,1:6)=(/ 2196.800000,2229.500000,2015.600000,1504.300000,603.440000,104.140000 /)
XPIZA_LKT(27,12,1:6)=(/ 0.949882,0.967826,0.976978,0.992751,0.992181,0.982467 /)
XCGA_LKT(27,12,1:6)=(/ 0.737623,0.729297,0.722503,0.720583,0.664877,0.421603 /)
XEXT_COEFF_550_LKT(27,12)=2026.500000 !rg=0.0834522 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(27,13,1:6)=(/ 1782.000000,1849.900000,1760.100000,1418.900000,708.540000,150.860000 /)
XPIZA_LKT(27,13,1:6)=(/ 0.941847,0.961962,0.973411,0.992044,0.992955,0.986988 /)
XCGA_LKT(27,13,1:6)=(/ 0.746260,0.731710,0.725250,0.721480,0.702830,0.507427 /)
XEXT_COEFF_550_LKT(27,13)=1763.800000 !rg=0.0834522 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(27,14,1:6)=(/ 1459.000000,1518.900000,1491.600000,1282.400000,777.160000,205.060000 /)
XPIZA_LKT(27,14,1:6)=(/ 0.931021,0.955137,0.968740,0.990962,0.993265,0.989722 /)
XCGA_LKT(27,14,1:6)=(/ 0.754720,0.735127,0.726300,0.719833,0.727497,0.580253 /)
XEXT_COEFF_550_LKT(27,14)=1498.500000 !rg=0.0834522 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(27,15,1:6)=(/ 1188.200000,1251.000000,1251.600000,1141.000000,793.600000,264.460000 /)
XPIZA_LKT(27,15,1:6)=(/ 0.919080,0.946100,0.963313,0.989601,0.993130,0.991489 /)
XCGA_LKT(27,15,1:6)=(/ 0.762380,0.741993,0.729433,0.721263,0.739583,0.639857 /)
XEXT_COEFF_550_LKT(27,15)=1252.400000 !rg=0.0834522 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(27,16,1:6)=(/ 964.560000,1016.600000,1040.400000,986.920000,756.160000,322.170000 /)
XPIZA_LKT(27,16,1:6)=(/ 0.906446,0.937053,0.955826,0.987883,0.992535,0.992583 /)
XCGA_LKT(27,16,1:6)=(/ 0.776967,0.749173,0.734217,0.724487,0.740917,0.685313 /)
XEXT_COEFF_550_LKT(27,16)=1034.900000 !rg=0.0834522 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(27,17,1:6)=(/ 789.630000,828.320000,855.200000,834.150000,684.230000,367.840000 /)
XPIZA_LKT(27,17,1:6)=(/ 0.891645,0.927414,0.949231,0.985798,0.991543,0.993158 /)
XCGA_LKT(27,17,1:6)=(/ 0.788090,0.759010,0.740827,0.726133,0.736947,0.716980 /)
XEXT_COEFF_550_LKT(27,17)=858.500000 !rg=0.0834522 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(27,18,1:6)=(/ 646.900000,680.430000,702.520000,706.610000,609.690000,392.190000 /)
XPIZA_LKT(27,18,1:6)=(/ 0.875539,0.914620,0.940542,0.982456,0.990288,0.993288 /)
XCGA_LKT(27,18,1:6)=(/ 0.796913,0.769400,0.746963,0.729827,0.736170,0.735927 /)
XEXT_COEFF_550_LKT(27,18)=701.410000 !rg=0.0834522 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(27,19,1:6)=(/ 528.150000,553.130000,576.950000,583.530000,538.050000,389.220000 /)
XPIZA_LKT(27,19,1:6)=(/ 0.861371,0.901021,0.929849,0.979533,0.988714,0.992977 /)
XCGA_LKT(27,19,1:6)=(/ 0.810800,0.779067,0.756537,0.733877,0.739373,0.742927 /)
XEXT_COEFF_550_LKT(27,19)=573.790000 !rg=0.0834522 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(27,20,1:6)=(/ 433.890000,451.680000,470.000000,481.330000,458.440000,361.080000 /)
XPIZA_LKT(27,20,1:6)=(/ 0.849070,0.886116,0.919288,0.976828,0.986920,0.992205 /)
XCGA_LKT(27,20,1:6)=(/ 0.820793,0.790910,0.766343,0.739873,0.741747,0.739897 /)
XEXT_COEFF_550_LKT(27,20)=471.620000 !rg=0.0834522 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(28,1,1:6)=(/ 7162.300000,2310.500000,645.280000,109.380000,12.369000,1.933900 /)
XPIZA_LKT(28,1,1:6)=(/ 0.985146,0.975993,0.955628,0.944522,0.769094,0.298472 /)
XCGA_LKT(28,1,1:6)=(/ 0.702147,0.532203,0.207323,0.079570,0.025367,0.006277 /)
XEXT_COEFF_550_LKT(28,1)=665.850000 !rg=0.0904084 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(28,2,1:6)=(/ 7473.600000,2634.800000,745.710000,134.650000,14.853000,2.087200 /)
XPIZA_LKT(28,2,1:6)=(/ 0.985458,0.978039,0.960668,0.954082,0.806182,0.348582 /)
XCGA_LKT(28,2,1:6)=(/ 0.720227,0.572093,0.263927,0.100887,0.032177,0.007973 /)
XEXT_COEFF_550_LKT(28,2)=767.890000 !rg=0.0904084 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(28,3,1:6)=(/ 7610.900000,3076.800000,941.680000,189.980000,20.836000,2.460100 /)
XPIZA_LKT(28,3,1:6)=(/ 0.985479,0.980551,0.967129,0.966313,0.859591,0.444788 /)
XCGA_LKT(28,3,1:6)=(/ 0.736170,0.612823,0.373473,0.152117,0.048827,0.012153 /)
XEXT_COEFF_550_LKT(28,3)=968.480000 !rg=0.0904084 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(28,4,1:6)=(/ 7456.700000,3535.800000,1240.700000,282.920000,32.984000,3.245600 /)
XPIZA_LKT(28,4,1:6)=(/ 0.984916,0.982448,0.973190,0.976213,0.908997,0.575856 /)
XCGA_LKT(28,4,1:6)=(/ 0.746160,0.656093,0.485807,0.246200,0.082410,0.020723 /)
XEXT_COEFF_550_LKT(28,4)=1273.600000 !rg=0.0904084 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(28,5,1:6)=(/ 6936.500000,3911.900000,1603.400000,421.040000,54.888000,4.788900 /)
XPIZA_LKT(28,5,1:6)=(/ 0.983496,0.983533,0.977739,0.982904,0.943323,0.708984 /)
XCGA_LKT(28,5,1:6)=(/ 0.748017,0.691410,0.571170,0.364043,0.136687,0.034737 /)
XEXT_COEFF_550_LKT(28,5)=1640.700000 !rg=0.0904084 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(28,6,1:6)=(/ 6156.400000,4100.300000,1969.500000,603.690000,89.331000,7.591500 /)
XPIZA_LKT(28,6,1:6)=(/ 0.981163,0.983744,0.980698,0.987121,0.963566,0.813136 /)
XCGA_LKT(28,6,1:6)=(/ 0.744003,0.714733,0.632250,0.470250,0.212113,0.054303 /)
XEXT_COEFF_550_LKT(28,6)=2008.600000 !rg=0.0904084 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(28,7,1:6)=(/ 5283.400000,4049.700000,2283.400000,823.060000,139.980000,12.523000 /)
XPIZA_LKT(28,7,1:6)=(/ 0.977843,0.983030,0.982408,0.989790,0.975383,0.883906 /)
XCGA_LKT(28,7,1:6)=(/ 0.738127,0.726827,0.675287,0.554573,0.309880,0.081307 /)
XEXT_COEFF_550_LKT(28,7)=2321.100000 !rg=0.0904084 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(28,8,1:6)=(/ 4437.300000,3778.100000,2478.500000,1056.200000,211.230000,20.806000 /)
XPIZA_LKT(28,8,1:6)=(/ 0.973488,0.981365,0.983028,0.991448,0.982538,0.927842 /)
XCGA_LKT(28,8,1:6)=(/ 0.732630,0.730217,0.702907,0.618107,0.411920,0.119177 /)
XEXT_COEFF_550_LKT(28,8)=2510.600000 !rg=0.0904084 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(28,9,1:6)=(/ 3678.600000,3377.000000,2514.300000,1274.800000,303.030000,34.111000 /)
XPIZA_LKT(28,9,1:6)=(/ 0.967989,0.978791,0.982616,0.992442,0.986940,0.954178 /)
XCGA_LKT(28,9,1:6)=(/ 0.729427,0.729777,0.717630,0.664950,0.498533,0.173690 /)
XEXT_COEFF_550_LKT(28,9)=2538.500000 !rg=0.0904084 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(28,10,1:6)=(/ 3017.000000,2938.800000,2393.100000,1438.300000,413.000000,53.772000 /)
XPIZA_LKT(28,10,1:6)=(/ 0.962265,0.975292,0.981164,0.992923,0.989696,0.969481 /)
XCGA_LKT(28,10,1:6)=(/ 0.731063,0.728737,0.722363,0.696803,0.575217,0.252320 /)
XEXT_COEFF_550_LKT(28,10)=2407.400000 !rg=0.0904084 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(28,11,1:6)=(/ 2462.900000,2482.800000,2172.400000,1513.400000,532.810000,81.942000 /)
XPIZA_LKT(28,11,1:6)=(/ 0.954889,0.971029,0.978799,0.992957,0.991457,0.978675 /)
XCGA_LKT(28,11,1:6)=(/ 0.734330,0.728980,0.722577,0.715317,0.635017,0.357710 /)
XEXT_COEFF_550_LKT(28,11)=2185.600000 !rg=0.0904084 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(28,12,1:6)=(/ 2002.000000,2068.400000,1922.100000,1483.200000,649.040000,122.250000 /)
XPIZA_LKT(28,12,1:6)=(/ 0.946809,0.965539,0.975698,0.992526,0.992549,0.984596 /)
XCGA_LKT(28,12,1:6)=(/ 0.740970,0.729690,0.724520,0.721887,0.681463,0.460257 /)
XEXT_COEFF_550_LKT(28,12)=1927.300000 !rg=0.0904084 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(28,13,1:6)=(/ 1635.400000,1703.400000,1652.600000,1369.900000,742.240000,171.960000 /)
XPIZA_LKT(28,13,1:6)=(/ 0.936758,0.959259,0.971623,0.991677,0.993130,0.988259 /)
XCGA_LKT(28,13,1:6)=(/ 0.749077,0.732680,0.726010,0.721557,0.714053,0.537057 /)
XEXT_COEFF_550_LKT(28,13)=1659.200000 !rg=0.0904084 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(28,14,1:6)=(/ 1331.600000,1403.900000,1389.100000,1229.700000,791.390000,229.630000 /)
XPIZA_LKT(28,14,1:6)=(/ 0.925457,0.951171,0.966696,0.990419,0.993260,0.990555 /)
XCGA_LKT(28,14,1:6)=(/ 0.756557,0.738107,0.727620,0.719960,0.733657,0.607310 /)
XEXT_COEFF_550_LKT(28,14)=1390.800000 !rg=0.0904084 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(28,15,1:6)=(/ 1087.100000,1144.200000,1165.600000,1081.000000,784.560000,288.880000 /)
XPIZA_LKT(28,15,1:6)=(/ 0.913075,0.942601,0.960038,0.988985,0.992944,0.991998 /)
XCGA_LKT(28,15,1:6)=(/ 0.769827,0.743420,0.731773,0.723537,0.741130,0.659700 /)
XEXT_COEFF_550_LKT(28,15)=1159.900000 !rg=0.0904084 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(28,16,1:6)=(/ 882.350000,928.270000,957.050000,922.990000,729.170000,342.980000 /)
XPIZA_LKT(28,16,1:6)=(/ 0.900405,0.933723,0.953408,0.987175,0.992184,0.992864 /)
XCGA_LKT(28,16,1:6)=(/ 0.781967,0.754607,0.736037,0.725743,0.739690,0.699807 /)
XEXT_COEFF_550_LKT(28,16)=959.660000 !rg=0.0904084 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(28,17,1:6)=(/ 723.330000,762.320000,785.110000,781.610000,653.960000,381.110000 /)
XPIZA_LKT(28,17,1:6)=(/ 0.883635,0.921705,0.946103,0.984478,0.991076,0.993256 /)
XCGA_LKT(28,17,1:6)=(/ 0.791140,0.763760,0.743220,0.727800,0.736570,0.726300 /)
XEXT_COEFF_550_LKT(28,17)=785.770000 !rg=0.0904084 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(28,18,1:6)=(/ 592.930000,620.990000,647.540000,651.580000,581.670000,394.290000 /)
XPIZA_LKT(28,18,1:6)=(/ 0.869019,0.909173,0.935961,0.981873,0.989720,0.993213 /)
XCGA_LKT(28,18,1:6)=(/ 0.804090,0.771860,0.751540,0.730290,0.737510,0.740173 /)
XEXT_COEFF_550_LKT(28,18)=642.880000 !rg=0.0904084 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(28,19,1:6)=(/ 484.570000,504.900000,527.240000,538.820000,505.310000,379.790000 /)
XPIZA_LKT(28,19,1:6)=(/ 0.855914,0.895174,0.925496,0.978626,0.988058,0.992707 /)
XCGA_LKT(28,19,1:6)=(/ 0.815727,0.785627,0.759627,0.737410,0.740840,0.742400 /)
XEXT_COEFF_550_LKT(28,19)=528.830000 !rg=0.0904084 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(28,20,1:6)=(/ 398.290000,415.230000,430.340000,446.390000,428.100000,345.570000 /)
XPIZA_LKT(28,20,1:6)=(/ 0.842885,0.878705,0.914493,0.974272,0.986104,0.991762 /)
XCGA_LKT(28,20,1:6)=(/ 0.823630,0.796563,0.769993,0.742730,0.744287,0.737573 /)
XEXT_COEFF_550_LKT(28,20)=430.250000 !rg=0.0904084 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(29,1,1:6)=(/ 7779.100000,2704.400000,775.280000,136.710000,14.970000,2.094100 /)
XPIZA_LKT(29,1,1:6)=(/ 0.985570,0.977944,0.961871,0.954536,0.807391,0.350463 /)
XCGA_LKT(29,1,1:6)=(/ 0.738977,0.597737,0.246903,0.093377,0.029743,0.007363 /)
XEXT_COEFF_550_LKT(29,1)=799.020000 !rg=0.0979445 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(29,2,1:6)=(/ 7822.200000,3027.500000,882.470000,167.670000,18.114000,2.288600 /)
XPIZA_LKT(29,2,1:6)=(/ 0.985836,0.979986,0.965567,0.962182,0.839380,0.404175 /)
XCGA_LKT(29,2,1:6)=(/ 0.739750,0.610080,0.313023,0.118447,0.037717,0.009353 /)
XEXT_COEFF_550_LKT(29,2)=907.770000 !rg=0.0979445 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(29,3,1:6)=(/ 7788.200000,3423.000000,1101.400000,233.270000,25.654000,2.761900 /)
XPIZA_LKT(29,3,1:6)=(/ 0.985594,0.981946,0.970693,0.971810,0.884498,0.503636 /)
XCGA_LKT(29,3,1:6)=(/ 0.748833,0.640300,0.425150,0.178247,0.057183,0.014253 /)
XEXT_COEFF_550_LKT(29,3)=1132.300000 !rg=0.0979445 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(29,4,1:6)=(/ 7425.100000,3822.900000,1422.600000,340.050000,40.743000,3.758200 /)
XPIZA_LKT(29,4,1:6)=(/ 0.984655,0.983333,0.975657,0.979593,0.925165,0.631893 /)
XCGA_LKT(29,4,1:6)=(/ 0.753260,0.677060,0.525240,0.282680,0.096283,0.024290 /)
XEXT_COEFF_550_LKT(29,4)=1458.900000 !rg=0.0979445 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(29,5,1:6)=(/ 6718.800000,4106.400000,1788.800000,495.240000,67.223000,5.711400 /)
XPIZA_LKT(29,5,1:6)=(/ 0.982806,0.983972,0.979340,0.984956,0.952839,0.754323 /)
XCGA_LKT(29,5,1:6)=(/ 0.750673,0.706020,0.599310,0.402330,0.159240,0.040683 /)
XEXT_COEFF_550_LKT(29,5)=1828.400000 !rg=0.0979445 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(29,6,1:6)=(/ 5838.800000,4176.800000,2135.300000,694.130000,107.740000,9.242300 /)
XPIZA_LKT(29,6,1:6)=(/ 0.979979,0.983759,0.981676,0.988398,0.969092,0.845074 /)
XCGA_LKT(29,6,1:6)=(/ 0.743973,0.723693,0.652257,0.503567,0.245493,0.063570 /)
XEXT_COEFF_550_LKT(29,6)=2175.100000 !rg=0.0979445 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(29,7,1:6)=(/ 4934.800000,4008.100000,2404.200000,922.740000,166.680000,15.396000 /)
XPIZA_LKT(29,7,1:6)=(/ 0.976159,0.982633,0.982894,0.990584,0.978741,0.904387 /)
XCGA_LKT(29,7,1:6)=(/ 0.736863,0.731317,0.688857,0.580763,0.349427,0.095183 /)
XEXT_COEFF_550_LKT(29,7)=2440.700000 !rg=0.0979445 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(29,8,1:6)=(/ 4095.700000,3656.100000,2533.200000,1153.400000,246.650000,25.535000 /)
XPIZA_LKT(29,8,1:6)=(/ 0.971469,0.980539,0.983075,0.991936,0.984605,0.940276 /)
XCGA_LKT(29,8,1:6)=(/ 0.731787,0.731820,0.711033,0.638090,0.446290,0.139597 /)
XEXT_COEFF_550_LKT(29,8)=2563.000000 !rg=0.0979445 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(29,9,1:6)=(/ 3373.100000,3221.000000,2497.600000,1354.800000,346.360000,41.339000 /)
XPIZA_LKT(29,9,1:6)=(/ 0.965208,0.977542,0.982244,0.992710,0.988216,0.961461 /)
XCGA_LKT(29,9,1:6)=(/ 0.728927,0.730543,0.721493,0.679087,0.530170,0.203557 /)
XEXT_COEFF_550_LKT(29,9)=2518.700000 !rg=0.0979445 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(29,10,1:6)=(/ 2768.200000,2758.700000,2325.000000,1485.300000,461.580000,64.093000 /)
XPIZA_LKT(29,10,1:6)=(/ 0.958656,0.973687,0.980369,0.993001,0.990517,0.973773 /)
XCGA_LKT(29,10,1:6)=(/ 0.730947,0.729597,0.723347,0.705527,0.600457,0.294017 /)
XEXT_COEFF_550_LKT(29,10)=2339.200000 !rg=0.0979445 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(29,11,1:6)=(/ 2249.300000,2310.300000,2083.300000,1518.000000,582.360000,96.950000 /)
XPIZA_LKT(29,11,1:6)=(/ 0.951420,0.968639,0.977712,0.992848,0.991977,0.981416 /)
XCGA_LKT(29,11,1:6)=(/ 0.735987,0.728373,0.724270,0.719263,0.655080,0.403103 /)
XEXT_COEFF_550_LKT(29,11)=2090.900000 !rg=0.0979445 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(29,12,1:6)=(/ 1837.000000,1912.000000,1820.500000,1448.600000,691.650000,141.750000 /)
XPIZA_LKT(29,12,1:6)=(/ 0.942035,0.963069,0.974154,0.992249,0.992842,0.986328 /)
XCGA_LKT(29,12,1:6)=(/ 0.743680,0.730840,0.725457,0.722707,0.695877,0.492840 /)
XEXT_COEFF_550_LKT(29,12)=1826.500000 !rg=0.0979445 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(29,13,1:6)=(/ 1491.900000,1572.900000,1543.400000,1318.800000,769.060000,194.700000 /)
XPIZA_LKT(29,13,1:6)=(/ 0.931499,0.955572,0.969801,0.991214,0.993238,0.989313 /)
XCGA_LKT(29,13,1:6)=(/ 0.750870,0.734397,0.726850,0.720827,0.723360,0.567160 /)
XEXT_COEFF_550_LKT(29,13)=1545.300000 !rg=0.0979445 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(29,14,1:6)=(/ 1219.800000,1285.400000,1297.600000,1171.900000,795.780000,253.940000 /)
XPIZA_LKT(29,14,1:6)=(/ 0.919707,0.947427,0.963887,0.989926,0.993187,0.991233 /)
XCGA_LKT(29,14,1:6)=(/ 0.763400,0.738707,0.729403,0.722127,0.737860,0.630080 /)
XEXT_COEFF_550_LKT(29,14)=1292.700000 !rg=0.0979445 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(29,15,1:6)=(/ 990.960000,1045.600000,1073.300000,1019.700000,767.030000,312.340000 /)
XPIZA_LKT(29,15,1:6)=(/ 0.907426,0.939351,0.957521,0.988345,0.992691,0.992422 /)
XCGA_LKT(29,15,1:6)=(/ 0.774507,0.749883,0.732160,0.725280,0.741337,0.677613 /)
XEXT_COEFF_550_LKT(29,15)=1074.500000 !rg=0.0979445 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(29,16,1:6)=(/ 807.480000,854.080000,880.460000,863.470000,700.470000,361.220000 /)
XPIZA_LKT(29,16,1:6)=(/ 0.893352,0.928379,0.951004,0.985970,0.991758,0.993085 /)
XCGA_LKT(29,16,1:6)=(/ 0.785670,0.758040,0.740343,0.725430,0.738237,0.711943 /)
XEXT_COEFF_550_LKT(29,16)=882.820000 !rg=0.0979445 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(29,17,1:6)=(/ 663.480000,696.640000,725.580000,726.680000,624.600000,389.780000 /)
XPIZA_LKT(29,17,1:6)=(/ 0.877097,0.916104,0.941552,0.983293,0.990554,0.993293 /)
XCGA_LKT(29,17,1:6)=(/ 0.798243,0.765537,0.746810,0.727747,0.736353,0.733173 /)
XEXT_COEFF_550_LKT(29,17)=719.980000 !rg=0.0979445 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(29,18,1:6)=(/ 542.260000,566.880000,592.550000,600.700000,551.160000,391.440000 /)
XPIZA_LKT(29,18,1:6)=(/ 0.862305,0.903488,0.931830,0.981030,0.989169,0.993054 /)
XCGA_LKT(29,18,1:6)=(/ 0.808873,0.779780,0.753207,0.735153,0.740307,0.742010 /)
XEXT_COEFF_550_LKT(29,18)=593.140000 !rg=0.0979445 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(29,19,1:6)=(/ 444.810000,464.560000,481.860000,498.500000,473.110000,367.330000 /)
XPIZA_LKT(29,19,1:6)=(/ 0.849837,0.887736,0.921659,0.976912,0.987215,0.992355 /)
XCGA_LKT(29,19,1:6)=(/ 0.819087,0.790440,0.765350,0.739000,0.742907,0.740717 /)
XEXT_COEFF_550_LKT(29,19)=483.240000 !rg=0.0979445 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(29,20,1:6)=(/ 365.970000,380.280000,396.180000,409.920000,400.090000,329.510000 /)
XPIZA_LKT(29,20,1:6)=(/ 0.837715,0.872187,0.908165,0.972828,0.985072,0.991253 /)
XCGA_LKT(29,20,1:6)=(/ 0.829337,0.799223,0.775127,0.742987,0.745653,0.734847 /)
XEXT_COEFF_550_LKT(29,20)=392.790000 !rg=0.0979445 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(30,1,1:6)=(/ 8176.600000,3189.800000,918.170000,170.710000,18.267000,2.297300 /)
XPIZA_LKT(30,1,1:6)=(/ 0.986259,0.980123,0.966721,0.962628,0.840424,0.406157 /)
XCGA_LKT(30,1,1:6)=(/ 0.754750,0.634890,0.295443,0.109677,0.034870,0.008637 /)
XEXT_COEFF_550_LKT(30,1)=944.410000 !rg=0.106109 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(30,2,1:6)=(/ 8061.900000,3421.200000,1033.900000,208.070000,22.242000,2.544300 /)
XPIZA_LKT(30,2,1:6)=(/ 0.986047,0.981706,0.969418,0.968689,0.867592,0.462216 /)
XCGA_LKT(30,2,1:6)=(/ 0.754210,0.635603,0.370133,0.139230,0.044210,0.010973 /)
XEXT_COEFF_550_LKT(30,2)=1062.700000 !rg=0.106109 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(30,3,1:6)=(/ 7845.400000,3755.400000,1280.000000,284.410000,31.715000,3.144900 /)
XPIZA_LKT(30,3,1:6)=(/ 0.985505,0.983092,0.973637,0.976203,0.905214,0.562220 /)
XCGA_LKT(30,3,1:6)=(/ 0.758603,0.663907,0.476253,0.208837,0.066963,0.016713 /)
XEXT_COEFF_550_LKT(30,3)=1315.100000 !rg=0.106109 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(30,4,1:6)=(/ 7276.200000,4081.200000,1616.100000,405.750000,50.318000,4.408000 /)
XPIZA_LKT(30,4,1:6)=(/ 0.984164,0.984021,0.977716,0.982325,0.938344,0.684361 /)
XCGA_LKT(30,4,1:6)=(/ 0.757850,0.695483,0.561340,0.322287,0.112423,0.028463 /)
XEXT_COEFF_550_LKT(30,4)=1655.600000 !rg=0.106109 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(30,5,1:6)=(/ 6419.900000,4256.900000,1975.400000,578.240000,81.991000,6.876800 /)
XPIZA_LKT(30,5,1:6)=(/ 0.981825,0.984236,0.980669,0.986645,0.960525,0.794355 /)
XCGA_LKT(30,5,1:6)=(/ 0.751253,0.718510,0.624713,0.440577,0.185123,0.047637 /)
XEXT_COEFF_550_LKT(30,5)=2016.500000 !rg=0.106109 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(30,6,1:6)=(/ 5474.500000,4202.800000,2290.500000,790.630000,129.390000,11.311000 /)
XPIZA_LKT(30,6,1:6)=(/ 0.978452,0.983610,0.982447,0.989452,0.973609,0.872064 /)
XCGA_LKT(30,6,1:6)=(/ 0.741920,0.730963,0.670087,0.534583,0.281653,0.074400 /)
XEXT_COEFF_550_LKT(30,6)=2330.200000 !rg=0.106109 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(30,7,1:6)=(/ 4562.100000,3926.500000,2503.300000,1025.000000,197.190000,18.938000 /)
XPIZA_LKT(30,7,1:6)=(/ 0.974138,0.982052,0.983203,0.991247,0.981506,0.921191 /)
XCGA_LKT(30,7,1:6)=(/ 0.734253,0.734700,0.700547,0.605033,0.387003,0.111413 /)
XEXT_COEFF_550_LKT(30,7)=2537.700000 !rg=0.106109 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(30,8,1:6)=(/ 3768.700000,3509.500000,2560.500000,1247.700000,285.190000,31.208000 /)
XPIZA_LKT(30,8,1:6)=(/ 0.968675,0.979573,0.982953,0.992328,0.986286,0.950290 /)
XCGA_LKT(30,8,1:6)=(/ 0.729340,0.733493,0.717550,0.655973,0.479143,0.163503 /)
XEXT_COEFF_550_LKT(30,8)=2587.100000 !rg=0.106109 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(30,9,1:6)=(/ 3090.000000,3041.400000,2458.800000,1424.900000,393.060000,49.739000 /)
XPIZA_LKT(30,9,1:6)=(/ 0.962263,0.976182,0.981673,0.992894,0.989293,0.967294 /)
XCGA_LKT(30,9,1:6)=(/ 0.730247,0.731467,0.724007,0.691383,0.560120,0.238197 /)
XEXT_COEFF_550_LKT(30,9)=2475.400000 !rg=0.106109 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(30,10,1:6)=(/ 2522.300000,2578.200000,2244.000000,1516.600000,511.170000,76.107000 /)
XPIZA_LKT(30,10,1:6)=(/ 0.954942,0.971566,0.979502,0.993007,0.991191,0.977310 /)
XCGA_LKT(30,10,1:6)=(/ 0.731050,0.729170,0.725040,0.712640,0.623277,0.339053 /)
XEXT_COEFF_550_LKT(30,10)=2255.400000 !rg=0.106109 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(30,11,1:6)=(/ 2061.200000,2138.200000,1988.400000,1505.200000,629.390000,114.140000 /)
XPIZA_LKT(30,11,1:6)=(/ 0.946682,0.966294,0.976432,0.992675,0.992389,0.983721 /)
XCGA_LKT(30,11,1:6)=(/ 0.738673,0.729077,0.725130,0.721870,0.672730,0.443170 /)
XEXT_COEFF_550_LKT(30,11)=1995.000000 !rg=0.106109 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(30,12,1:6)=(/ 1674.000000,1767.500000,1712.200000,1406.200000,728.280000,162.270000 /)
XPIZA_LKT(30,12,1:6)=(/ 0.937117,0.959614,0.972537,0.991876,0.993056,0.987717 /)
XCGA_LKT(30,12,1:6)=(/ 0.745123,0.731970,0.726797,0.722267,0.708040,0.522910 /)
XEXT_COEFF_550_LKT(30,12)=1714.300000 !rg=0.106109 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(30,13,1:6)=(/ 1365.900000,1442.900000,1441.800000,1262.600000,787.310000,218.950000 /)
XPIZA_LKT(30,13,1:6)=(/ 0.926129,0.952074,0.967232,0.990769,0.993269,0.990210 /)
XCGA_LKT(30,13,1:6)=(/ 0.757347,0.735050,0.727357,0.721820,0.730370,0.595457 /)
XEXT_COEFF_550_LKT(30,13)=1437.500000 !rg=0.106109 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(30,14,1:6)=(/ 1109.500000,1176.200000,1201.300000,1115.000000,790.980000,278.340000 /)
XPIZA_LKT(30,14,1:6)=(/ 0.914554,0.943897,0.961285,0.989360,0.993041,0.991784 /)
XCGA_LKT(30,14,1:6)=(/ 0.767330,0.744660,0.730030,0.724173,0.740347,0.650570 /)
XEXT_COEFF_550_LKT(30,14)=1202.200000 !rg=0.106109 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(30,15,1:6)=(/ 903.850000,957.980000,987.580000,955.760000,744.120000,333.840000 /)
XPIZA_LKT(30,15,1:6)=(/ 0.902109,0.934745,0.955369,0.987465,0.992352,0.992741 /)
XCGA_LKT(30,15,1:6)=(/ 0.779763,0.752673,0.736757,0.725207,0.740357,0.693007 /)
XEXT_COEFF_550_LKT(30,15)=989.570000 !rg=0.106109 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(30,16,1:6)=(/ 741.850000,781.270000,811.950000,804.230000,670.070000,375.890000 /)
XPIZA_LKT(30,16,1:6)=(/ 0.885371,0.922704,0.946879,0.984940,0.991290,0.993216 /)
XCGA_LKT(30,16,1:6)=(/ 0.792150,0.760093,0.742230,0.726097,0.736857,0.721983 /)
XEXT_COEFF_550_LKT(30,16)=809.500000 !rg=0.106109 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(30,17,1:6)=(/ 605.200000,638.060000,664.880000,673.410000,595.470000,394.080000 /)
XPIZA_LKT(30,17,1:6)=(/ 0.870728,0.910094,0.937317,0.981979,0.990031,0.993244 /)
XCGA_LKT(30,17,1:6)=(/ 0.802003,0.772917,0.747953,0.731907,0.738393,0.738157 /)
XEXT_COEFF_550_LKT(30,17)=663.320000 !rg=0.106109 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(30,18,1:6)=(/ 495.850000,519.890000,541.720000,556.380000,521.230000,384.070000 /)
XPIZA_LKT(30,18,1:6)=(/ 0.857329,0.896897,0.927855,0.979396,0.988257,0.992818 /)
XCGA_LKT(30,18,1:6)=(/ 0.814057,0.784507,0.760190,0.736477,0.741173,0.742367 /)
XEXT_COEFF_550_LKT(30,18)=543.170000 !rg=0.106109 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(30,19,1:6)=(/ 408.720000,425.570000,443.700000,458.800000,440.430000,352.320000 /)
XPIZA_LKT(30,19,1:6)=(/ 0.843641,0.880194,0.915670,0.974862,0.986485,0.991931 /)
XCGA_LKT(30,19,1:6)=(/ 0.824537,0.793360,0.769230,0.739027,0.743510,0.737980 /)
XEXT_COEFF_550_LKT(30,19)=441.850000 !rg=0.106109 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(30,20,1:6)=(/ 334.660000,348.430000,362.600000,376.490000,373.690000,313.680000 /)
XPIZA_LKT(30,20,1:6)=(/ 0.832496,0.866194,0.901783,0.970416,0.983760,0.990756 /)
XCGA_LKT(30,20,1:6)=(/ 0.832523,0.806213,0.777340,0.748673,0.745660,0.734367 /)
XEXT_COEFF_550_LKT(30,20)=361.780000 !rg=0.106109 sigma=2.95 wvl=0.55

IF (LHOOK) CALL DR_HOOK('MODE_DUSTOPT:DUST_OPT_LKT_SET3',1,ZHOOK_HANDLE)
END SUBROUTINE DUST_OPT_LKT_SET3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE DUST_OPT_LKT_SET4()

  USE MODD_DUST_OPT_LKT
  
  IMPLICIT NONE
 
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_DUSTOPT:DUST_OPT_LKT_SET4',0,ZHOOK_HANDLE)
XEXT_COEFF_WVL_LKT(31,1,1:6)=(/ 8217.700000,3668.300000,1071.000000,212.620000,22.443000,2.555300 /)
XPIZA_LKT(31,1,1:6)=(/ 0.986440,0.982322,0.970401,0.969142,0.868490,0.464248 /)
XCGA_LKT(31,1,1:6)=(/ 0.758317,0.644733,0.354757,0.129003,0.040880,0.010133 /)
XEXT_COEFF_550_LKT(31,1)=1099.400000 !rg=0.114954 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(31,2,1:6)=(/ 8171.600000,3783.600000,1203.400000,256.800000,27.459000,2.868800 /)
XPIZA_LKT(31,2,1:6)=(/ 0.986060,0.983121,0.972487,0.973894,0.891260,0.521150 /)
XCGA_LKT(31,2,1:6)=(/ 0.764690,0.655150,0.433227,0.163923,0.051813,0.012870 /)
XEXT_COEFF_550_LKT(31,2)=1236.800000 !rg=0.114954 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(31,3,1:6)=(/ 7769.400000,4062.200000,1477.000000,343.950000,39.310000,3.630900 /)
XPIZA_LKT(31,3,1:6)=(/ 0.985189,0.984006,0.976099,0.979710,0.922280,0.618930 /)
XCGA_LKT(31,3,1:6)=(/ 0.765350,0.684993,0.523843,0.244443,0.078407,0.019597 /)
XEXT_COEFF_550_LKT(31,3)=1516.400000 !rg=0.114954 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(31,4,1:6)=(/ 7011.600000,4300.300000,1817.700000,480.560000,62.036000,5.231100 /)
XPIZA_LKT(31,4,1:6)=(/ 0.983410,0.984521,0.979436,0.984544,0.949027,0.732263 /)
XCGA_LKT(31,4,1:6)=(/ 0.759940,0.711590,0.593687,0.364297,0.131163,0.033343 /)
XEXT_COEFF_550_LKT(31,4)=1859.700000 !rg=0.114954 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(31,5,1:6)=(/ 6043.800000,4355.700000,2157.700000,669.610000,99.543000,8.345800 /)
XPIZA_LKT(31,5,1:6)=(/ 0.980555,0.984328,0.981758,0.988036,0.966744,0.829028 /)
XCGA_LKT(31,5,1:6)=(/ 0.749740,0.728947,0.647410,0.477673,0.214430,0.055760 /)
XEXT_COEFF_550_LKT(31,5)=2199.800000 !rg=0.114954 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(31,6,1:6)=(/ 5071.600000,4177.900000,2429.100000,892.750000,154.660000,13.890000 /)
XPIZA_LKT(31,6,1:6)=(/ 0.976718,0.983268,0.983036,0.990327,0.977317,0.894573 /)
XCGA_LKT(31,6,1:6)=(/ 0.739773,0.736403,0.685687,0.563640,0.318990,0.087053 /)
XEXT_COEFF_550_LKT(31,6)=2467.900000 !rg=0.114954 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(31,7,1:6)=(/ 4186.800000,3807.700000,2575.800000,1126.600000,231.210000,23.262000 /)
XPIZA_LKT(31,7,1:6)=(/ 0.971890,0.981265,0.983348,0.991788,0.983759,0.934859 /)
XCGA_LKT(31,7,1:6)=(/ 0.732163,0.736600,0.710300,0.626933,0.422327,0.130390 /)
XEXT_COEFF_550_LKT(31,7)=2608.100000 !rg=0.114954 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(31,8,1:6)=(/ 3440.100000,3345.700000,2559.300000,1334.500000,327.410000,37.920000 /)
XPIZA_LKT(31,8,1:6)=(/ 0.965917,0.978339,0.982680,0.992639,0.987678,0.958321 /)
XCGA_LKT(31,8,1:6)=(/ 0.728437,0.733800,0.722680,0.671660,0.511697,0.191413 /)
XEXT_COEFF_550_LKT(31,8)=2582.400000 !rg=0.114954 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(31,9,1:6)=(/ 2809.300000,2855.000000,2396.800000,1480.800000,441.430000,59.464000 /)
XPIZA_LKT(31,9,1:6)=(/ 0.959285,0.974601,0.981021,0.993012,0.990185,0.972002 /)
XCGA_LKT(31,9,1:6)=(/ 0.729760,0.731793,0.726123,0.701477,0.586697,0.277593 /)
XEXT_COEFF_550_LKT(31,9)=2413.800000 !rg=0.114954 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(31,10,1:6)=(/ 2306.300000,2388.200000,2157.200000,1530.800000,561.430000,90.197000 /)
XPIZA_LKT(31,10,1:6)=(/ 0.950778,0.969356,0.978389,0.992940,0.991757,0.980279 /)
XCGA_LKT(31,10,1:6)=(/ 0.734697,0.728850,0.725933,0.717747,0.644427,0.384017 /)
XEXT_COEFF_550_LKT(31,10)=2164.600000 !rg=0.114954 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(31,11,1:6)=(/ 1874.700000,1981.900000,1883.700000,1479.500000,674.120000,132.930000 /)
XPIZA_LKT(31,11,1:6)=(/ 0.942355,0.963271,0.975084,0.992418,0.992716,0.985608 /)
XCGA_LKT(31,11,1:6)=(/ 0.739493,0.730570,0.727297,0.723120,0.688233,0.477187 /)
XEXT_COEFF_550_LKT(31,11)=1889.300000 !rg=0.114954 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(31,12,1:6)=(/ 1530.400000,1619.400000,1603.800000,1354.500000,759.030000,184.280000 /)
XPIZA_LKT(31,12,1:6)=(/ 0.931838,0.956472,0.970417,0.991499,0.993194,0.988857 /)
XCGA_LKT(31,12,1:6)=(/ 0.751217,0.731933,0.727000,0.722653,0.718450,0.553207 /)
XEXT_COEFF_550_LKT(31,12)=1601.100000 !rg=0.114954 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(31,13,1:6)=(/ 1241.300000,1322.500000,1338.200000,1208.000000,796.430000,243.310000 /)
XPIZA_LKT(31,13,1:6)=(/ 0.921329,0.948637,0.964803,0.990254,0.993232,0.990951 /)
XCGA_LKT(31,13,1:6)=(/ 0.760830,0.740293,0.728017,0.723360,0.735660,0.619563 /)
XEXT_COEFF_550_LKT(31,13)=1340.000000 !rg=0.114954 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(31,14,1:6)=(/ 1009.700000,1075.400000,1108.400000,1052.900000,778.010000,302.430000 /)
XPIZA_LKT(31,14,1:6)=(/ 0.909846,0.940240,0.958939,0.988643,0.992816,0.992247 /)
XCGA_LKT(31,14,1:6)=(/ 0.774173,0.747603,0.733873,0.725077,0.741223,0.669427 /)
XEXT_COEFF_550_LKT(31,14)=1107.600000 !rg=0.114954 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(31,15,1:6)=(/ 830.240000,876.930000,911.090000,891.430000,716.050000,353.480000 /)
XPIZA_LKT(31,15,1:6)=(/ 0.894798,0.929645,0.951938,0.986582,0.991970,0.992989 /)
XCGA_LKT(31,15,1:6)=(/ 0.785757,0.755177,0.738610,0.725663,0.738947,0.706147 /)
XEXT_COEFF_550_LKT(31,15)=913.680000 !rg=0.114954 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(31,16,1:6)=(/ 677.620000,716.670000,745.180000,751.840000,640.400000,386.650000 /)
XPIZA_LKT(31,16,1:6)=(/ 0.878532,0.917149,0.942573,0.983467,0.990810,0.993281 /)
XCGA_LKT(31,16,1:6)=(/ 0.795173,0.766587,0.743597,0.728910,0.737387,0.729973 /)
XEXT_COEFF_550_LKT(31,16)=742.230000 !rg=0.114954 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(31,17,1:6)=(/ 552.700000,582.400000,609.910000,619.910000,567.010000,393.160000 /)
XPIZA_LKT(31,17,1:6)=(/ 0.864466,0.904335,0.932660,0.981038,0.989412,0.993123 /)
XCGA_LKT(31,17,1:6)=(/ 0.809657,0.777483,0.754037,0.733283,0.739987,0.740957 /)
XEXT_COEFF_550_LKT(31,17)=608.510000 !rg=0.114954 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(31,18,1:6)=(/ 456.050000,476.460000,496.990000,513.170000,486.610000,372.770000 /)
XPIZA_LKT(31,18,1:6)=(/ 0.851195,0.889378,0.922886,0.977228,0.987723,0.992507 /)
XCGA_LKT(31,18,1:6)=(/ 0.819527,0.787840,0.763857,0.736833,0.742567,0.741107 /)
XEXT_COEFF_550_LKT(31,18)=498.630000 !rg=0.114954 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(31,19,1:6)=(/ 374.380000,390.530000,406.610000,423.570000,411.870000,336.420000 /)
XPIZA_LKT(31,19,1:6)=(/ 0.838298,0.873596,0.909511,0.972488,0.985532,0.991477 /)
XCGA_LKT(31,19,1:6)=(/ 0.827663,0.800110,0.771647,0.743723,0.744287,0.736273 /)
XEXT_COEFF_550_LKT(31,19)=404.710000 !rg=0.114954 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(31,20,1:6)=(/ 306.770000,318.970000,332.250000,344.540000,347.860000,299.380000 /)
XPIZA_LKT(31,20,1:6)=(/ 0.827072,0.859694,0.894837,0.968770,0.982169,0.990209 /)
XCGA_LKT(31,20,1:6)=(/ 0.838360,0.811013,0.784667,0.752720,0.747703,0.734800 /)
XEXT_COEFF_550_LKT(31,20)=331.900000 !rg=0.114954 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(32,1,1:6)=(/ 8362.400000,3998.600000,1234.700000,263.590000,27.726000,2.882800 /)
XPIZA_LKT(32,1,1:6)=(/ 0.986184,0.984028,0.973146,0.974365,0.892031,0.523168 /)
XCGA_LKT(32,1,1:6)=(/ 0.768243,0.647667,0.425307,0.152037,0.047923,0.011887 /)
XEXT_COEFF_550_LKT(32,1)=1266.300000 !rg=0.124536 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(32,2,1:6)=(/ 8139.200000,4104.000000,1397.100000,314.480000,34.039000,3.280600 /)
XPIZA_LKT(32,2,1:6)=(/ 0.985842,0.984188,0.975026,0.978038,0.910901,0.579323 /)
XCGA_LKT(32,2,1:6)=(/ 0.772363,0.675267,0.497000,0.193393,0.060730,0.015093 /)
XEXT_COEFF_550_LKT(32,2)=1436.400000 !rg=0.124536 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(32,3,1:6)=(/ 7553.000000,4333.600000,1690.100000,412.350000,48.778000,4.247300 /)
XPIZA_LKT(32,3,1:6)=(/ 0.984606,0.984700,0.978183,0.982510,0.936229,0.672359 /)
XCGA_LKT(32,3,1:6)=(/ 0.768977,0.704243,0.565700,0.285350,0.091797,0.022973 /)
XEXT_COEFF_550_LKT(32,3)=1733.200000 !rg=0.124536 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(32,4,1:6)=(/ 6642.700000,4470.500000,2022.200000,564.900000,76.248000,6.272800 /)
XPIZA_LKT(32,4,1:6)=(/ 0.982317,0.984837,0.980870,0.986354,0.957656,0.775021 /)
XCGA_LKT(32,4,1:6)=(/ 0.759137,0.725490,0.622273,0.407587,0.152863,0.039053 /)
XEXT_COEFF_550_LKT(32,4)=2066.100000 !rg=0.124536 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(32,5,1:6)=(/ 5612.600000,4397.700000,2330.600000,768.720000,120.260000,10.192000 /)
XPIZA_LKT(32,5,1:6)=(/ 0.978957,0.984238,0.982636,0.989186,0.971790,0.858563 /)
XCGA_LKT(32,5,1:6)=(/ 0.746873,0.737343,0.667560,0.512970,0.246850,0.065243 /)
XEXT_COEFF_550_LKT(32,5)=2372.800000 !rg=0.124536 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(32,6,1:6)=(/ 4663.400000,4101.000000,2546.500000,998.590000,183.680000,17.084000 /)
XPIZA_LKT(32,6,1:6)=(/ 0.974434,0.982771,0.983449,0.991056,0.980352,0.913134 /)
XCGA_LKT(32,6,1:6)=(/ 0.735620,0.740503,0.699270,0.590547,0.356030,0.101823 /)
XEXT_COEFF_550_LKT(32,6)=2583.600000 !rg=0.124536 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(32,7,1:6)=(/ 3830.800000,3653.700000,2620.200000,1226.500000,268.670000,28.481000 /)
XPIZA_LKT(32,7,1:6)=(/ 0.968953,0.980347,0.983328,0.992229,0.985597,0.945903 /)
XCGA_LKT(32,7,1:6)=(/ 0.728763,0.738207,0.718353,0.646687,0.456560,0.152550 /)
XEXT_COEFF_550_LKT(32,7)=2649.100000 !rg=0.124536 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(32,8,1:6)=(/ 3142.600000,3155.200000,2532.400000,1412.600000,373.220000,45.773000 /)
XPIZA_LKT(32,8,1:6)=(/ 0.962381,0.977021,0.982219,0.992865,0.988848,0.964763 /)
XCGA_LKT(32,8,1:6)=(/ 0.727957,0.734493,0.726317,0.685447,0.542767,0.223757 /)
XEXT_COEFF_550_LKT(32,8)=2551.800000 !rg=0.124536 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(32,9,1:6)=(/ 2553.200000,2659.900000,2321.600000,1521.800000,491.100000,70.784000 /)
XPIZA_LKT(32,9,1:6)=(/ 0.956040,0.972597,0.980185,0.993057,0.990919,0.975861 /)
XCGA_LKT(32,9,1:6)=(/ 0.730623,0.730610,0.727910,0.709900,0.610730,0.320620 /)
XEXT_COEFF_550_LKT(32,9)=2333.000000 !rg=0.124536 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(32,10,1:6)=(/ 2090.200000,2208.100000,2056.800000,1527.900000,609.870000,106.450000 /)
XPIZA_LKT(32,10,1:6)=(/ 0.947281,0.967452,0.977289,0.992806,0.992216,0.982774 /)
XCGA_LKT(32,10,1:6)=(/ 0.736510,0.730613,0.727687,0.721540,0.663223,0.424900 /)
XEXT_COEFF_550_LKT(32,10)=2067.100000 !rg=0.124536 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(32,11,1:6)=(/ 1706.300000,1817.100000,1776.600000,1440.300000,713.710000,152.840000 /)
XPIZA_LKT(32,11,1:6)=(/ 0.938086,0.960785,0.973297,0.992121,0.992970,0.987125 /)
XCGA_LKT(32,11,1:6)=(/ 0.746477,0.730847,0.727777,0.723857,0.701460,0.507973 /)
XEXT_COEFF_550_LKT(32,11)=1777.200000 !rg=0.124536 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(32,12,1:6)=(/ 1389.800000,1484.200000,1490.700000,1302.200000,781.420000,208.020000 /)
XPIZA_LKT(32,12,1:6)=(/ 0.927343,0.953684,0.968155,0.991057,0.993263,0.989825 /)
XCGA_LKT(32,12,1:6)=(/ 0.754997,0.736700,0.727143,0.723447,0.726470,0.582487 /)
XEXT_COEFF_550_LKT(32,12)=1493.200000 !rg=0.124536 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(32,13,1:6)=(/ 1128.400000,1206.700000,1241.200000,1149.700000,796.470000,267.640000 /)
XPIZA_LKT(32,13,1:6)=(/ 0.916919,0.945214,0.962579,0.989644,0.993121,0.991549 /)
XCGA_LKT(32,13,1:6)=(/ 0.768080,0.742507,0.731700,0.724233,0.739097,0.640807 /)
XEXT_COEFF_550_LKT(32,13)=1239.900000 !rg=0.124536 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(32,14,1:6)=(/ 927.640000,982.530000,1019.500000,987.820000,757.220000,324.660000 /)
XPIZA_LKT(32,14,1:6)=(/ 0.902927,0.935983,0.956179,0.987904,0.992533,0.992607 /)
XCGA_LKT(32,14,1:6)=(/ 0.779250,0.750207,0.734897,0.725870,0.741040,0.685770 /)
XEXT_COEFF_550_LKT(32,14)=1023.900000 !rg=0.124536 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(32,15,1:6)=(/ 759.180000,805.320000,836.430000,834.680000,686.810000,369.570000 /)
XPIZA_LKT(32,15,1:6)=(/ 0.887020,0.924052,0.948226,0.985159,0.991549,0.993159 /)
XCGA_LKT(32,15,1:6)=(/ 0.788530,0.760680,0.740113,0.726820,0.738450,0.716977 /)
XEXT_COEFF_550_LKT(32,15)=835.750000 !rg=0.124536 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(32,16,1:6)=(/ 619.540000,652.800000,685.810000,693.710000,611.790000,392.840000 /)
XPIZA_LKT(32,16,1:6)=(/ 0.871714,0.911930,0.938282,0.982750,0.990300,0.993266 /)
XCGA_LKT(32,16,1:6)=(/ 0.803460,0.770037,0.748840,0.730523,0.738367,0.735627 /)
XEXT_COEFF_550_LKT(32,16)=681.190000 !rg=0.124536 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(32,17,1:6)=(/ 508.500000,532.090000,557.460000,572.630000,534.750000,387.700000 /)
XPIZA_LKT(32,17,1:6)=(/ 0.858036,0.898637,0.928527,0.979915,0.988798,0.992923 /)
XCGA_LKT(32,17,1:6)=(/ 0.814540,0.782330,0.757763,0.734880,0.741827,0.742020 /)
XEXT_COEFF_550_LKT(32,17)=559.300000 !rg=0.124536 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(32,18,1:6)=(/ 418.280000,437.640000,455.620000,474.200000,454.780000,359.000000 /)
XPIZA_LKT(32,18,1:6)=(/ 0.844706,0.881887,0.917148,0.975322,0.986874,0.992127 /)
XCGA_LKT(32,18,1:6)=(/ 0.822337,0.794247,0.766537,0.740207,0.743503,0.739430 /)
XEXT_COEFF_550_LKT(32,18)=454.940000 !rg=0.124536 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(32,19,1:6)=(/ 343.310000,356.450000,373.240000,386.710000,385.470000,321.200000 /)
XPIZA_LKT(32,19,1:6)=(/ 0.832778,0.867452,0.903211,0.971086,0.984188,0.990951 /)
XCGA_LKT(32,19,1:6)=(/ 0.833840,0.804160,0.778170,0.746523,0.746453,0.734557 /)
XEXT_COEFF_550_LKT(32,19)=371.080000 !rg=0.124536 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(32,20,1:6)=(/ 282.790000,292.060000,303.830000,316.060000,319.990000,285.110000 /)
XPIZA_LKT(32,20,1:6)=(/ 0.822103,0.854898,0.887888,0.966801,0.981836,0.989618 /)
XCGA_LKT(32,20,1:6)=(/ 0.842177,0.815350,0.789413,0.754977,0.749923,0.736200 /)
XEXT_COEFF_550_LKT(32,20)=304.570000 !rg=0.124536 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(33,1,1:6)=(/ 8374.300000,4199.300000,1420.500000,324.480000,34.401000,3.298400 /)
XPIZA_LKT(33,1,1:6)=(/ 0.986178,0.984891,0.975265,0.978531,0.911574,0.581281 /)
XCGA_LKT(33,1,1:6)=(/ 0.781217,0.665920,0.503047,0.179680,0.056187,0.013943 /)
XEXT_COEFF_550_LKT(33,1)=1459.200000 !rg=0.134917 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(33,2,1:6)=(/ 7950.900000,4392.200000,1620.300000,381.350000,42.316000,3.803300 /)
XPIZA_LKT(33,2,1:6)=(/ 0.985354,0.984945,0.977249,0.981318,0.927063,0.635188 /)
XCGA_LKT(33,2,1:6)=(/ 0.777147,0.698197,0.553780,0.228703,0.071187,0.017703 /)
XEXT_COEFF_550_LKT(33,2)=1666.000000 !rg=0.134917 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(33,3,1:6)=(/ 7197.300000,4560.200000,1913.600000,490.120000,60.504000,5.029200 /)
XPIZA_LKT(33,3,1:6)=(/ 0.983701,0.985192,0.979952,0.984754,0.947558,0.721444 /)
XCGA_LKT(33,3,1:6)=(/ 0.769340,0.721650,0.601017,0.331330,0.107477,0.026927 /)
XEXT_COEFF_550_LKT(33,3)=1959.600000 !rg=0.134917 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(33,4,1:6)=(/ 6178.700000,4583.000000,2223.700000,658.920000,93.316000,7.589200 /)
XPIZA_LKT(33,4,1:6)=(/ 0.980894,0.984968,0.982056,0.987838,0.964611,0.812417 /)
XCGA_LKT(33,4,1:6)=(/ 0.755927,0.737163,0.647373,0.450777,0.177857,0.045727 /)
XEXT_COEFF_550_LKT(33,4)=2268.700000 !rg=0.134917 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(33,5,1:6)=(/ 5150.800000,4379.400000,2488.000000,874.820000,144.490000,12.503000 /)
XPIZA_LKT(33,5,1:6)=(/ 0.976854,0.983967,0.983321,0.990140,0.975895,0.883369 /)
XCGA_LKT(33,5,1:6)=(/ 0.741290,0.743830,0.685267,0.546173,0.281657,0.076310 /)
XEXT_COEFF_550_LKT(33,5)=2529.600000 !rg=0.134917 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(33,6,1:6)=(/ 4235.100000,3980.200000,2638.000000,1105.800000,216.390000,21.006000 /)
XPIZA_LKT(33,6,1:6)=(/ 0.971956,0.982036,0.983693,0.991656,0.982826,0.928303 /)
XCGA_LKT(33,6,1:6)=(/ 0.731430,0.742947,0.710813,0.615033,0.392377,0.119050 /)
XEXT_COEFF_550_LKT(33,6)=2672.800000 !rg=0.134917 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(33,7,1:6)=(/ 3470.800000,3477.300000,2634.200000,1320.600000,309.990000,34.703000 /)
XPIZA_LKT(33,7,1:6)=(/ 0.965923,0.979137,0.983154,0.992583,0.987120,0.954787 /)
XCGA_LKT(33,7,1:6)=(/ 0.726440,0.738353,0.724833,0.664220,0.490443,0.178353 /)
XEXT_COEFF_550_LKT(33,7)=2659.900000 !rg=0.134917 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(33,8,1:6)=(/ 2841.300000,2957.400000,2479.400000,1477.800000,421.270000,54.909000 /)
XPIZA_LKT(33,8,1:6)=(/ 0.959296,0.975372,0.981658,0.993023,0.989820,0.969960 /)
XCGA_LKT(33,8,1:6)=(/ 0.727027,0.734053,0.729337,0.697093,0.570823,0.260610 /)
XEXT_COEFF_550_LKT(33,8)=2497.700000 !rg=0.134917 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(33,9,1:6)=(/ 2330.800000,2460.900000,2232.800000,1545.700000,541.970000,84.053000 /)
XPIZA_LKT(33,9,1:6)=(/ 0.951602,0.970509,0.979168,0.993031,0.991535,0.979083 /)
XCGA_LKT(33,9,1:6)=(/ 0.731867,0.730290,0.728843,0.716340,0.633083,0.364497 /)
XEXT_COEFF_550_LKT(33,9)=2243.600000 !rg=0.134917 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(33,10,1:6)=(/ 1896.300000,2034.700000,1952.600000,1509.600000,656.490000,124.460000 /)
XPIZA_LKT(33,10,1:6)=(/ 0.943826,0.964882,0.975958,0.992602,0.992581,0.984826 /)
XCGA_LKT(33,10,1:6)=(/ 0.740150,0.730233,0.729037,0.723813,0.679877,0.460200 /)
XEXT_COEFF_550_LKT(33,10)=1956.600000 !rg=0.134917 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(33,11,1:6)=(/ 1549.800000,1662.200000,1658.400000,1394.600000,748.010000,174.170000 /)
XPIZA_LKT(33,11,1:6)=(/ 0.933793,0.958410,0.971499,0.991760,0.993141,0.988363 /)
XCGA_LKT(33,11,1:6)=(/ 0.749963,0.733513,0.728153,0.724633,0.712930,0.538537 /)
XEXT_COEFF_550_LKT(33,11)=1663.700000 !rg=0.134917 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(33,12,1:6)=(/ 1262.500000,1357.300000,1384.300000,1246.100000,795.340000,232.310000 /)
XPIZA_LKT(33,12,1:6)=(/ 0.923257,0.950299,0.966339,0.990561,0.993262,0.990632 /)
XCGA_LKT(33,12,1:6)=(/ 0.761190,0.738390,0.730063,0.724003,0.732883,0.607980 /)
XEXT_COEFF_550_LKT(33,12)=1384.500000 !rg=0.134917 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(33,13,1:6)=(/ 1035.700000,1102.800000,1143.500000,1088.100000,787.310000,292.130000 /)
XPIZA_LKT(33,13,1:6)=(/ 0.910317,0.941465,0.960152,0.989003,0.992939,0.992052 /)
XCGA_LKT(33,13,1:6)=(/ 0.772810,0.745317,0.732483,0.725863,0.740867,0.660537 /)
XEXT_COEFF_550_LKT(33,13)=1147.900000 !rg=0.134917 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(33,14,1:6)=(/ 848.020000,903.470000,937.810000,924.900000,732.260000,345.390000 /)
XPIZA_LKT(33,14,1:6)=(/ 0.895623,0.930349,0.953176,0.986775,0.992186,0.992883 /)
XCGA_LKT(33,14,1:6)=(/ 0.782157,0.754960,0.736827,0.725870,0.740553,0.699890 /)
XEXT_COEFF_550_LKT(33,14)=938.290000 !rg=0.134917 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(33,15,1:6)=(/ 695.350000,733.970000,770.860000,774.790000,657.880000,382.350000 /)
XPIZA_LKT(33,15,1:6)=(/ 0.879556,0.919060,0.943760,0.984414,0.991055,0.993252 /)
XCGA_LKT(33,15,1:6)=(/ 0.796490,0.762827,0.743897,0.727743,0.737760,0.726027 /)
XEXT_COEFF_550_LKT(33,15)=765.140000 !rg=0.134917 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(33,16,1:6)=(/ 567.240000,594.700000,625.760000,639.470000,581.400000,394.160000 /)
XPIZA_LKT(33,16,1:6)=(/ 0.865326,0.906939,0.934378,0.981777,0.989760,0.993179 /)
XCGA_LKT(33,16,1:6)=(/ 0.808947,0.777203,0.751197,0.733147,0.740607,0.739457 /)
XEXT_COEFF_550_LKT(33,16)=627.510000 !rg=0.134917 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(33,17,1:6)=(/ 466.750000,488.970000,509.650000,530.600000,502.240000,378.270000 /)
XPIZA_LKT(33,17,1:6)=(/ 0.851237,0.891131,0.924329,0.977905,0.988084,0.992646 /)
XCGA_LKT(33,17,1:6)=(/ 0.817587,0.788180,0.761677,0.737313,0.743353,0.741657 /)
XEXT_COEFF_550_LKT(33,17)=509.900000 !rg=0.134917 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(33,18,1:6)=(/ 383.840000,399.730000,418.850000,434.510000,425.410000,343.950000 /)
XPIZA_LKT(33,18,1:6)=(/ 0.838987,0.875673,0.911122,0.974272,0.985634,0.991664 /)
XCGA_LKT(33,18,1:6)=(/ 0.828397,0.797313,0.772290,0.741580,0.744213,0.736753 /)
XEXT_COEFF_550_LKT(33,18)=415.420000 !rg=0.134917 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(33,19,1:6)=(/ 315.110000,325.730000,340.430000,354.190000,357.390000,305.840000 /)
XPIZA_LKT(33,19,1:6)=(/ 0.827591,0.861956,0.896994,0.969643,0.983339,0.990469 /)
XCGA_LKT(33,19,1:6)=(/ 0.838187,0.810927,0.781970,0.751607,0.746953,0.735203 /)
XEXT_COEFF_550_LKT(33,19)=341.170000 !rg=0.134917 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(33,20,1:6)=(/ 259.760000,268.480000,278.020000,291.590000,296.190000,270.820000 /)
XPIZA_LKT(33,20,1:6)=(/ 0.816470,0.848637,0.881927,0.963651,0.980713,0.989084 /)
XCGA_LKT(33,20,1:6)=(/ 0.844800,0.820750,0.794083,0.758830,0.753560,0.738917 /)
XEXT_COEFF_550_LKT(33,20)=277.820000 !rg=0.134917 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(34,1,1:6)=(/ 8021.400000,4463.100000,1652.500000,395.470000,42.815000,3.825900 /)
XPIZA_LKT(34,1,1:6)=(/ 0.985488,0.985157,0.977180,0.981829,0.927658,0.637028 /)
XCGA_LKT(34,1,1:6)=(/ 0.784750,0.704137,0.575467,0.213120,0.065880,0.016353 /)
XEXT_COEFF_550_LKT(34,1)=1703.500000 !rg=0.146163 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(34,2,1:6)=(/ 7598.400000,4651.600000,1869.400000,457.100000,52.685000,4.466600 /)
XPIZA_LKT(34,2,1:6)=(/ 0.984533,0.985480,0.979268,0.983899,0.940266,0.687424 /)
XCGA_LKT(34,2,1:6)=(/ 0.778280,0.721250,0.597333,0.271010,0.083467,0.020757 /)
XEXT_COEFF_550_LKT(34,2)=1919.900000 !rg=0.146163 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(34,3,1:6)=(/ 6713.200000,4732.100000,2140.000000,577.920000,74.901000,6.020100 /)
XPIZA_LKT(34,3,1:6)=(/ 0.982385,0.985489,0.981445,0.986561,0.956713,0.765500 /)
XCGA_LKT(34,3,1:6)=(/ 0.766147,0.736827,0.630393,0.381293,0.125840,0.031557 /)
XEXT_COEFF_550_LKT(34,3)=2187.500000 !rg=0.146163 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(34,4,1:6)=(/ 5655.400000,4630.400000,2415.700000,762.430000,113.600000,9.249800 /)
XPIZA_LKT(34,4,1:6)=(/ 0.978964,0.984908,0.983020,0.989063,0.970215,0.844550 /)
XCGA_LKT(34,4,1:6)=(/ 0.749777,0.746547,0.669400,0.492537,0.206380,0.053520 /)
XEXT_COEFF_550_LKT(34,4)=2460.900000 !rg=0.146163 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(34,5,1:6)=(/ 4658.900000,4301.800000,2624.200000,986.140000,172.460000,15.380000 /)
XPIZA_LKT(34,5,1:6)=(/ 0.974471,0.983469,0.983825,0.990934,0.979232,0.903952 /)
XCGA_LKT(34,5,1:6)=(/ 0.735573,0.748200,0.700720,0.576830,0.318093,0.089207 /)
XEXT_COEFF_550_LKT(34,5)=2664.200000 !rg=0.146163 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(34,6,1:6)=(/ 3830.100000,3814.300000,2699.800000,1212.400000,252.910000,25.775000 /)
XPIZA_LKT(34,6,1:6)=(/ 0.969203,0.981105,0.983773,0.992150,0.984850,0.940610 /)
XCGA_LKT(34,6,1:6)=(/ 0.727710,0.743940,0.720513,0.637213,0.428447,0.139107 /)
XEXT_COEFF_550_LKT(34,6)=2731.400000 !rg=0.146163 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(34,7,1:6)=(/ 3142.500000,3270.700000,2619.500000,1407.200000,355.080000,42.037000 /)
XPIZA_LKT(34,7,1:6)=(/ 0.962720,0.977752,0.982791,0.992853,0.988393,0.961928 /)
XCGA_LKT(34,7,1:6)=(/ 0.725730,0.737943,0.729717,0.679763,0.522950,0.208203 /)
XEXT_COEFF_550_LKT(34,7)=2640.700000 !rg=0.146163 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(34,8,1:6)=(/ 2569.800000,2745.600000,2409.000000,1529.000000,471.110000,65.566000 /)
XPIZA_LKT(34,8,1:6)=(/ 0.956152,0.973424,0.980884,0.993109,0.990622,0.974206 /)
XCGA_LKT(34,8,1:6)=(/ 0.728563,0.732477,0.731450,0.706993,0.596390,0.301190 /)
XEXT_COEFF_550_LKT(34,8)=2421.400000 !rg=0.146163 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(34,9,1:6)=(/ 2110.700000,2276.200000,2132.000000,1552.700000,591.820000,99.429000 /)
XPIZA_LKT(34,9,1:6)=(/ 0.947042,0.967674,0.978050,0.992937,0.992041,0.981782 /)
XCGA_LKT(34,9,1:6)=(/ 0.732280,0.729603,0.730483,0.721320,0.653150,0.405533 /)
XEXT_COEFF_550_LKT(34,9)=2139.800000 !rg=0.146163 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(34,10,1:6)=(/ 1734.000000,1866.000000,1838.500000,1477.700000,698.970000,143.730000 /)
XPIZA_LKT(34,10,1:6)=(/ 0.938564,0.962023,0.974328,0.992347,0.992873,0.986482 /)
XCGA_LKT(34,10,1:6)=(/ 0.743493,0.730503,0.729210,0.725683,0.694290,0.491907 /)
XEXT_COEFF_550_LKT(34,10)=1843.600000 !rg=0.146163 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(34,11,1:6)=(/ 1410.800000,1525.200000,1542.400000,1343.500000,774.580000,197.310000 /)
XPIZA_LKT(34,11,1:6)=(/ 0.929013,0.954612,0.969783,0.991321,0.993247,0.989408 /)
XCGA_LKT(34,11,1:6)=(/ 0.753613,0.734110,0.729337,0.724650,0.722097,0.568593 /)
XEXT_COEFF_550_LKT(34,11)=1544.800000 !rg=0.146163 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(34,12,1:6)=(/ 1158.200000,1238.500000,1281.400000,1187.500000,800.210000,256.590000 /)
XPIZA_LKT(34,12,1:6)=(/ 0.916957,0.946520,0.964043,0.990010,0.993187,0.991286 /)
XCGA_LKT(34,12,1:6)=(/ 0.766127,0.740147,0.730693,0.725977,0.737237,0.630143 /)
XEXT_COEFF_550_LKT(34,12)=1286.200000 !rg=0.146163 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(34,13,1:6)=(/ 946.350000,1012.800000,1051.000000,1024.100000,770.930000,315.130000 /)
XPIZA_LKT(34,13,1:6)=(/ 0.903309,0.936170,0.957411,0.988155,0.992688,0.992456 /)
XCGA_LKT(34,13,1:6)=(/ 0.775703,0.749317,0.733683,0.726070,0.741600,0.677890 /)
XEXT_COEFF_550_LKT(34,13)=1051.400000 !rg=0.146163 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(34,14,1:6)=(/ 777.210000,823.330000,865.250000,860.220000,704.370000,362.970000 /)
XPIZA_LKT(34,14,1:6)=(/ 0.888103,0.925396,0.948996,0.985772,0.991771,0.993091 /)
XCGA_LKT(34,14,1:6)=(/ 0.790067,0.756430,0.739783,0.725597,0.739367,0.711620 /)
XEXT_COEFF_550_LKT(34,14)=859.430000 !rg=0.146163 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(34,15,1:6)=(/ 634.670000,668.850000,704.950000,716.890000,628.090000,390.490000 /)
XPIZA_LKT(34,15,1:6)=(/ 0.872445,0.913933,0.940099,0.983731,0.990581,0.993273 /)
XCGA_LKT(34,15,1:6)=(/ 0.801607,0.771033,0.745417,0.732203,0.739180,0.732530 /)
XEXT_COEFF_550_LKT(34,15)=705.480000 !rg=0.146163 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(34,16,1:6)=(/ 520.060000,546.650000,571.440000,592.410000,551.150000,390.830000 /)
XPIZA_LKT(34,16,1:6)=(/ 0.859198,0.899914,0.930907,0.980132,0.989126,0.993016 /)
XCGA_LKT(34,16,1:6)=(/ 0.812873,0.781893,0.757300,0.733870,0.742350,0.741427 /)
XEXT_COEFF_550_LKT(34,16)=573.350000 !rg=0.146163 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(34,17,1:6)=(/ 428.610000,447.330000,468.950000,487.220000,468.920000,365.580000 /)
XPIZA_LKT(34,17,1:6)=(/ 0.845541,0.883470,0.918526,0.976262,0.987217,0.992290 /)
XCGA_LKT(34,17,1:6)=(/ 0.823987,0.790893,0.766537,0.736680,0.743253,0.739863 /)
XEXT_COEFF_550_LKT(34,17)=464.970000 !rg=0.146163 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(34,18,1:6)=(/ 351.800000,365.330000,382.230000,397.560000,397.020000,328.130000 /)
XPIZA_LKT(34,18,1:6)=(/ 0.832817,0.868835,0.905410,0.972678,0.984674,0.991201 /)
XCGA_LKT(34,18,1:6)=(/ 0.832667,0.805253,0.774710,0.748023,0.744433,0.735677 /)
XEXT_COEFF_550_LKT(34,18)=382.190000 !rg=0.146163 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(34,19,1:6)=(/ 289.690000,299.840000,310.710000,326.610000,330.770000,291.800000 /)
XPIZA_LKT(34,19,1:6)=(/ 0.822287,0.855818,0.890792,0.967074,0.982086,0.989936 /)
XCGA_LKT(34,19,1:6)=(/ 0.841170,0.815557,0.788913,0.753860,0.750360,0.736917 /)
XEXT_COEFF_550_LKT(34,19)=311.790000 !rg=0.146163 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(34,20,1:6)=(/ 238.770000,246.350000,255.590000,266.940000,273.760000,255.730000 /)
XPIZA_LKT(34,20,1:6)=(/ 0.812061,0.842637,0.874955,0.961206,0.979131,0.988332 /)
XCGA_LKT(34,20,1:6)=(/ 0.849113,0.823490,0.799660,0.759410,0.755320,0.740200 /)
XEXT_COEFF_550_LKT(34,20)=253.540000 !rg=0.146163 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(35,1,1:6)=(/ 7695.500000,4847.200000,1948.700000,475.750000,53.387000,4.495500 /)
XPIZA_LKT(35,1,1:6)=(/ 0.984488,0.985601,0.979269,0.984413,0.940806,0.689118 /)
XCGA_LKT(35,1,1:6)=(/ 0.782563,0.740233,0.624567,0.253913,0.077270,0.019180 /)
XEXT_COEFF_550_LKT(35,1)=2011.500000 !rg=0.158347 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(35,2,1:6)=(/ 7084.000000,4865.300000,2127.900000,541.430000,65.603000,5.308100 /)
XPIZA_LKT(35,2,1:6)=(/ 0.983280,0.985850,0.981083,0.985919,0.950989,0.735080 /)
XCGA_LKT(35,2,1:6)=(/ 0.775093,0.740543,0.627007,0.321267,0.097910,0.024340 /)
XEXT_COEFF_550_LKT(35,2)=2179.500000 !rg=0.158347 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(35,3,1:6)=(/ 6121.900000,4839.000000,2360.400000,676.610000,92.391000,7.275200 /)
XPIZA_LKT(35,3,1:6)=(/ 0.980560,0.985594,0.982684,0.988033,0.964079,0.804230 /)
XCGA_LKT(35,3,1:6)=(/ 0.759010,0.749360,0.655287,0.433153,0.147363,0.036980 /)
XEXT_COEFF_550_LKT(35,3)=2408.400000 !rg=0.158347 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(35,4,1:6)=(/ 5086.000000,4607.300000,2591.600000,874.670000,137.430000,11.339000 /)
XPIZA_LKT(35,4,1:6)=(/ 0.976520,0.984638,0.983781,0.990080,0.974731,0.871745 /)
XCGA_LKT(35,4,1:6)=(/ 0.740763,0.753540,0.688757,0.531713,0.238520,0.062620 /)
XEXT_COEFF_550_LKT(35,4)=2636.000000 !rg=0.158347 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(35,5,1:6)=(/ 4188.300000,4164.000000,2733.500000,1100.700000,204.350000,18.940000 /)
XPIZA_LKT(35,5,1:6)=(/ 0.971454,0.982764,0.984156,0.991593,0.981948,0.920867 /)
XCGA_LKT(35,5,1:6)=(/ 0.728390,0.750710,0.714003,0.604743,0.355743,0.104223 /)
XEXT_COEFF_550_LKT(35,5)=2771.100000 !rg=0.158347 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(35,6,1:6)=(/ 3447.800000,3614.500000,2729.300000,1315.100000,293.530000,31.510000 /)
XPIZA_LKT(35,6,1:6)=(/ 0.965663,0.979917,0.983686,0.992552,0.986524,0.950546 /)
XCGA_LKT(35,6,1:6)=(/ 0.722173,0.743943,0.728367,0.657040,0.464323,0.162393 /)
XEXT_COEFF_550_LKT(35,6)=2757.700000 !rg=0.158347 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(35,7,1:6)=(/ 2829.200000,3052.100000,2575.000000,1482.300000,402.970000,50.630000 /)
XPIZA_LKT(35,7,1:6)=(/ 0.959131,0.976093,0.982295,0.993052,0.989450,0.967691 /)
XCGA_LKT(35,7,1:6)=(/ 0.723527,0.736833,0.733453,0.693160,0.552830,0.242237 /)
XEXT_COEFF_550_LKT(35,7)=2594.700000 !rg=0.158347 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(35,8,1:6)=(/ 2334.600000,2533.000000,2319.200000,1563.500000,522.610000,78.056000 /)
XPIZA_LKT(35,8,1:6)=(/ 0.951510,0.971305,0.979955,0.993126,0.991295,0.977730 /)
XCGA_LKT(35,8,1:6)=(/ 0.728617,0.731443,0.732713,0.714933,0.620213,0.343283 /)
XEXT_COEFF_550_LKT(35,8)=2331.700000 !rg=0.158347 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(35,9,1:6)=(/ 1917.900000,2085.700000,2024.300000,1542.700000,640.290000,116.640000 /)
XPIZA_LKT(35,9,1:6)=(/ 0.942296,0.965123,0.976695,0.992782,0.992447,0.984005 /)
XCGA_LKT(35,9,1:6)=(/ 0.737977,0.729027,0.730960,0.724903,0.671063,0.441807 /)
XEXT_COEFF_550_LKT(35,9)=2028.000000 !rg=0.158347 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(35,10,1:6)=(/ 1573.800000,1713.100000,1718.300000,1437.200000,736.600000,164.430000 /)
XPIZA_LKT(35,10,1:6)=(/ 0.933126,0.958210,0.972578,0.991998,0.993080,0.987829 /)
XCGA_LKT(35,10,1:6)=(/ 0.745093,0.731043,0.729637,0.726327,0.706867,0.522933 /)
XEXT_COEFF_550_LKT(35,10)=1719.600000 !rg=0.158347 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(35,11,1:6)=(/ 1292.900000,1392.300000,1433.000000,1287.100000,793.360000,221.410000 /)
XPIZA_LKT(35,11,1:6)=(/ 0.922646,0.950673,0.967226,0.990881,0.993281,0.990286 /)
XCGA_LKT(35,11,1:6)=(/ 0.759463,0.735193,0.728887,0.726253,0.729620,0.595410 /)
XEXT_COEFF_550_LKT(35,11)=1433.100000 !rg=0.158347 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(35,12,1:6)=(/ 1056.400000,1137.400000,1181.200000,1127.100000,795.860000,281.260000 /)
XPIZA_LKT(35,12,1:6)=(/ 0.910404,0.941466,0.961390,0.989320,0.993046,0.991831 /)
XCGA_LKT(35,12,1:6)=(/ 0.768727,0.743907,0.731610,0.726457,0.740170,0.650720 /)
XEXT_COEFF_550_LKT(35,12)=1181.300000 !rg=0.158347 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(35,13,1:6)=(/ 866.790000,922.940000,970.190000,955.270000,748.540000,336.720000 /)
XPIZA_LKT(35,13,1:6)=(/ 0.896567,0.931300,0.953617,0.987231,0.992364,0.992764 /)
XCGA_LKT(35,13,1:6)=(/ 0.783703,0.750600,0.735867,0.725930,0.741070,0.692970 /)
XEXT_COEFF_550_LKT(35,13)=964.520000 !rg=0.158347 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(35,14,1:6)=(/ 707.880000,751.600000,791.320000,799.940000,674.730000,377.490000 /)
XPIZA_LKT(35,14,1:6)=(/ 0.881127,0.920022,0.945052,0.984973,0.991349,0.993212 /)
XCGA_LKT(35,14,1:6)=(/ 0.794700,0.764420,0.740357,0.729943,0.739590,0.721590 /)
XEXT_COEFF_550_LKT(35,14)=790.460000 !rg=0.158347 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(35,15,1:6)=(/ 579.850000,612.050000,643.520000,662.520000,599.030000,394.090000 /)
XPIZA_LKT(35,15,1:6)=(/ 0.866902,0.908130,0.936446,0.982292,0.990076,0.993222 /)
XCGA_LKT(35,15,1:6)=(/ 0.807560,0.775740,0.751977,0.732047,0.740950,0.737410 /)
XEXT_COEFF_550_LKT(35,15)=644.860000 !rg=0.158347 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(35,16,1:6)=(/ 478.060000,500.410000,525.090000,545.680000,517.570000,383.120000 /)
XPIZA_LKT(35,16,1:6)=(/ 0.852000,0.892288,0.925351,0.978531,0.988439,0.992772 /)
XCGA_LKT(35,16,1:6)=(/ 0.818987,0.784837,0.760637,0.734287,0.743183,0.741577 /)
XEXT_COEFF_550_LKT(35,16)=523.760000 !rg=0.158347 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(35,17,1:6)=(/ 391.700000,409.780000,428.680000,448.370000,438.000000,350.810000 /)
XPIZA_LKT(35,17,1:6)=(/ 0.839961,0.876477,0.912591,0.974054,0.986216,0.991900 /)
XCGA_LKT(35,17,1:6)=(/ 0.827340,0.798440,0.768470,0.742173,0.743073,0.738500 /)
XEXT_COEFF_550_LKT(35,17)=427.190000 !rg=0.158347 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(35,18,1:6)=(/ 322.150000,335.180000,349.030000,365.200000,370.410000,313.500000 /)
XPIZA_LKT(35,18,1:6)=(/ 0.828300,0.862985,0.899267,0.970386,0.983047,0.990690 /)
XCGA_LKT(35,18,1:6)=(/ 0.837440,0.810370,0.782753,0.750337,0.746910,0.735530 /)
XEXT_COEFF_550_LKT(35,18)=349.890000 !rg=0.158347 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(35,19,1:6)=(/ 266.450000,275.090000,285.880000,299.680000,305.080000,277.660000 /)
XPIZA_LKT(35,19,1:6)=(/ 0.816859,0.849322,0.883536,0.964348,0.980750,0.989288 /)
XCGA_LKT(35,19,1:6)=(/ 0.845593,0.818510,0.793643,0.754863,0.752253,0.738123 /)
XEXT_COEFF_550_LKT(35,19)=285.110000 !rg=0.158347 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(35,20,1:6)=(/ 218.970000,225.910000,234.150000,244.810000,253.090000,239.760000 /)
XPIZA_LKT(35,20,1:6)=(/ 0.807614,0.837418,0.868340,0.957940,0.976631,0.987573 /)
XCGA_LKT(35,20,1:6)=(/ 0.851753,0.829300,0.802337,0.766007,0.754050,0.740803 /)
XEXT_COEFF_550_LKT(35,20)=233.450000 !rg=0.158347 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(36,1,1:6)=(/ 7035.200000,5081.100000,2273.300000,563.570000,66.610000,5.345200 /)
XPIZA_LKT(36,1,1:6)=(/ 0.983319,0.986285,0.981515,0.986402,0.951495,0.736600 /)
XCGA_LKT(36,1,1:6)=(/ 0.777573,0.755030,0.643170,0.303900,0.090677,0.022490 /)
XEXT_COEFF_550_LKT(36,1)=2333.600000 !rg=0.171546 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(36,2,1:6)=(/ 6427.000000,5010.600000,2372.900000,634.760000,81.572000,6.375300 /)
XPIZA_LKT(36,2,1:6)=(/ 0.981449,0.986052,0.982623,0.987502,0.959655,0.777587 /)
XCGA_LKT(36,2,1:6)=(/ 0.767047,0.754780,0.647943,0.379377,0.114923,0.028533 /)
XEXT_COEFF_550_LKT(36,2)=2422.600000 !rg=0.171546 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(36,3,1:6)=(/ 5456.900000,4870.700000,2566.900000,786.840000,113.370000,8.863200 /)
XPIZA_LKT(36,3,1:6)=(/ 0.978101,0.985496,0.983685,0.989247,0.969988,0.837666 /)
XCGA_LKT(36,3,1:6)=(/ 0.748080,0.758990,0.677207,0.484063,0.172590,0.043323 /)
XEXT_COEFF_550_LKT(36,3)=2614.600000 !rg=0.171546 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(36,4,1:6)=(/ 4511.900000,4511.600000,2744.400000,994.070000,165.110000,13.958000 /)
XPIZA_LKT(36,4,1:6)=(/ 0.973501,0.984133,0.984352,0.990926,0.978375,0.894465 /)
XCGA_LKT(36,4,1:6)=(/ 0.730613,0.757993,0.705723,0.567467,0.274220,0.073233 /)
XEXT_COEFF_550_LKT(36,4)=2787.200000 !rg=0.171546 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(36,5,1:6)=(/ 3725.600000,3975.500000,2810.900000,1215.900000,240.400000,23.307000 /)
XPIZA_LKT(36,5,1:6)=(/ 0.968020,0.981780,0.984314,0.992139,0.984170,0.934656 /)
XCGA_LKT(36,5,1:6)=(/ 0.720350,0.751290,0.725210,0.629900,0.394327,0.121670 /)
XEXT_COEFF_550_LKT(36,5)=2845.400000 !rg=0.171546 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(36,6,1:6)=(/ 3086.700000,3389.000000,2725.500000,1411.200000,338.140000,38.337000 /)
XPIZA_LKT(36,6,1:6)=(/ 0.961927,0.978356,0.983427,0.992868,0.987917,0.958554 /)
XCGA_LKT(36,6,1:6)=(/ 0.719233,0.741783,0.734603,0.674690,0.498987,0.189277 /)
XEXT_COEFF_550_LKT(36,6)=2749.400000 !rg=0.171546 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(36,7,1:6)=(/ 2547.300000,2823.700000,2506.700000,1543.900000,453.250000,60.692000 /)
XPIZA_LKT(36,7,1:6)=(/ 0.955233,0.973978,0.981606,0.993179,0.990327,0.972389 /)
XCGA_LKT(36,7,1:6)=(/ 0.723353,0.733897,0.736193,0.704710,0.580423,0.279923 /)
XEXT_COEFF_550_LKT(36,7)=2520.700000 !rg=0.171546 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(36,8,1:6)=(/ 2107.700000,2332.100000,2215.500000,1580.900000,573.870000,92.561000 /)
XPIZA_LKT(36,8,1:6)=(/ 0.946860,0.968493,0.978856,0.993075,0.991852,0.980666 /)
XCGA_LKT(36,8,1:6)=(/ 0.729603,0.729383,0.734047,0.721303,0.641793,0.383740 /)
XEXT_COEFF_550_LKT(36,8)=2223.100000 !rg=0.171546 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(36,9,1:6)=(/ 1734.000000,1905.800000,1902.100000,1518.700000,685.550000,135.280000 /)
XPIZA_LKT(36,9,1:6)=(/ 0.938171,0.962753,0.975194,0.992553,0.992777,0.985806 /)
XCGA_LKT(36,9,1:6)=(/ 0.741387,0.730683,0.731600,0.727473,0.686780,0.474543 /)
XEXT_COEFF_550_LKT(36,9)=1909.200000 !rg=0.171546 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(36,10,1:6)=(/ 1436.000000,1560.000000,1601.100000,1386.300000,767.380000,186.940000 /)
XPIZA_LKT(36,10,1:6)=(/ 0.927391,0.954746,0.970214,0.991637,0.993223,0.988959 /)
XCGA_LKT(36,10,1:6)=(/ 0.752603,0.731020,0.729077,0.727490,0.717290,0.553660 /)
XEXT_COEFF_550_LKT(36,10)=1597.000000 !rg=0.171546 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(36,11,1:6)=(/ 1176.000000,1274.900000,1323.800000,1229.500000,802.940000,245.640000 /)
XPIZA_LKT(36,11,1:6)=(/ 0.916931,0.946104,0.964801,0.990353,0.993246,0.991001 /)
XCGA_LKT(36,11,1:6)=(/ 0.761530,0.738813,0.729693,0.727310,0.735083,0.618660 /)
XEXT_COEFF_550_LKT(36,11)=1322.900000 !rg=0.171546 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(36,12,1:6)=(/ 966.800000,1035.100000,1089.700000,1058.300000,783.610000,304.970000 /)
XPIZA_LKT(36,12,1:6)=(/ 0.903737,0.937013,0.957778,0.988642,0.992831,0.992280 /)
XCGA_LKT(36,12,1:6)=(/ 0.776950,0.745000,0.732753,0.727593,0.741490,0.669127 /)
XEXT_COEFF_550_LKT(36,12)=1083.300000 !rg=0.171546 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(36,13,1:6)=(/ 788.950000,842.490000,888.930000,889.810000,721.810000,355.810000 /)
XPIZA_LKT(36,13,1:6)=(/ 0.889768,0.926061,0.949862,0.986404,0.992009,0.993007 /)
XCGA_LKT(36,13,1:6)=(/ 0.788037,0.758080,0.736517,0.728463,0.741027,0.705717 /)
XEXT_COEFF_550_LKT(36,13)=887.760000 !rg=0.171546 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(36,14,1:6)=(/ 645.420000,685.130000,724.570000,740.900000,646.110000,387.610000 /)
XPIZA_LKT(36,14,1:6)=(/ 0.874851,0.914917,0.941180,0.984097,0.990865,0.993270 /)
XCGA_LKT(36,14,1:6)=(/ 0.802477,0.769090,0.746657,0.730773,0.740150,0.729107 /)
XEXT_COEFF_550_LKT(36,14)=723.920000 !rg=0.171546 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(36,15,1:6)=(/ 533.130000,560.060000,589.310000,610.710000,567.600000,393.060000 /)
XPIZA_LKT(36,15,1:6)=(/ 0.860102,0.901355,0.931883,0.980632,0.989465,0.993091 /)
XCGA_LKT(36,15,1:6)=(/ 0.813533,0.779093,0.755253,0.732133,0.742437,0.740107 /)
XEXT_COEFF_550_LKT(36,15)=591.300000 !rg=0.171546 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(36,16,1:6)=(/ 437.590000,458.780000,480.690000,503.550000,484.300000,371.920000 /)
XPIZA_LKT(36,16,1:6)=(/ 0.846047,0.885059,0.919460,0.976296,0.987662,0.992465 /)
XCGA_LKT(36,16,1:6)=(/ 0.822080,0.792217,0.763017,0.737757,0.743567,0.740960 /)
XEXT_COEFF_550_LKT(36,16)=478.440000 !rg=0.171546 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(36,17,1:6)=(/ 358.850000,374.490000,392.710000,409.640000,410.030000,335.970000 /)
XPIZA_LKT(36,17,1:6)=(/ 0.833858,0.869545,0.906301,0.972535,0.984815,0.991427 /)
XCGA_LKT(36,17,1:6)=(/ 0.833967,0.803407,0.775757,0.745523,0.745107,0.736730 /)
XEXT_COEFF_550_LKT(36,17)=391.790000 !rg=0.171546 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(36,18,1:6)=(/ 296.560000,307.520000,319.890000,335.530000,340.650000,298.880000 /)
XPIZA_LKT(36,18,1:6)=(/ 0.823392,0.856638,0.892479,0.967557,0.982400,0.990166 /)
XCGA_LKT(36,18,1:6)=(/ 0.841807,0.813760,0.787713,0.751443,0.749370,0.736383 /)
XEXT_COEFF_550_LKT(36,18)=321.040000 !rg=0.171546 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(36,19,1:6)=(/ 244.390000,252.660000,262.000000,275.370000,282.390000,262.650000 /)
XPIZA_LKT(36,19,1:6)=(/ 0.812269,0.843567,0.876071,0.961115,0.979562,0.988692 /)
XCGA_LKT(36,19,1:6)=(/ 0.848397,0.824583,0.796617,0.760083,0.753273,0.740260 /)
XEXT_COEFF_550_LKT(36,19)=261.100000 !rg=0.171546 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(36,20,1:6)=(/ 201.060000,207.380000,214.590000,223.480000,234.010000,224.770000 /)
XPIZA_LKT(36,20,1:6)=(/ 0.802173,0.831316,0.862351,0.954974,0.975028,0.986187 /)
XCGA_LKT(36,20,1:6)=(/ 0.856147,0.833640,0.809093,0.770477,0.758077,0.741150 /)
XEXT_COEFF_550_LKT(36,20)=214.450000 !rg=0.171546 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(37,1,1:6)=(/ 6328.200000,5103.300000,2531.500000,657.200000,83.045000,6.423100 /)
XPIZA_LKT(37,1,1:6)=(/ 0.980862,0.986430,0.983468,0.987902,0.960145,0.778928 /)
XCGA_LKT(37,1,1:6)=(/ 0.768640,0.758550,0.645757,0.364793,0.106493,0.026373 /)
XEXT_COEFF_550_LKT(37,1)=2573.900000 !rg=0.185845 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(37,2,1:6)=(/ 5659.500000,5074.700000,2591.300000,739.400000,101.110000,7.728200 /)
XPIZA_LKT(37,2,1:6)=(/ 0.978824,0.986055,0.983822,0.988764,0.966628,0.814759 /)
XCGA_LKT(37,2,1:6)=(/ 0.753113,0.765107,0.667277,0.442953,0.135037,0.033450 /)
XEXT_COEFF_550_LKT(37,2)=2638.900000 !rg=0.185845 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(37,3,1:6)=(/ 4768.400000,4819.500000,2752.800000,908.450000,138.160000,10.869000 /)
XPIZA_LKT(37,3,1:6)=(/ 0.974728,0.985168,0.984461,0.990261,0.974717,0.866081 /)
XCGA_LKT(37,3,1:6)=(/ 0.731870,0.765590,0.697123,0.531117,0.202100,0.050747 /)
XEXT_COEFF_550_LKT(37,3)=2799.300000 !rg=0.185845 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(37,4,1:6)=(/ 3970.000000,4343.800000,2867.800000,1118.400000,196.930000,17.226000 /)
XPIZA_LKT(37,4,1:6)=(/ 0.969596,0.983368,0.984738,0.991633,0.981323,0.913245 /)
XCGA_LKT(37,4,1:6)=(/ 0.717043,0.759967,0.720460,0.599433,0.313213,0.085600 /)
XEXT_COEFF_550_LKT(37,4)=2908.000000 !rg=0.185845 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(37,5,1:6)=(/ 3306.200000,3741.400000,2852.700000,1328.600000,280.820000,28.618000 /)
XPIZA_LKT(37,5,1:6)=(/ 0.964148,0.980474,0.984295,0.992586,0.986000,0.945834 /)
XCGA_LKT(37,5,1:6)=(/ 0.714377,0.749547,0.734403,0.652357,0.433257,0.141893 /)
XEXT_COEFF_550_LKT(37,5)=2883.000000 !rg=0.185845 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(37,6,1:6)=(/ 2775.700000,3136.400000,2689.000000,1497.200000,386.140000,46.406000 /)
XPIZA_LKT(37,6,1:6)=(/ 0.957214,0.976628,0.982967,0.993110,0.989071,0.965021 /)
XCGA_LKT(37,6,1:6)=(/ 0.716670,0.739587,0.739013,0.690103,0.531510,0.219937 /)
XEXT_COEFF_550_LKT(37,6)=2708.700000 !rg=0.185845 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(37,7,1:6)=(/ 2307.300000,2588.600000,2416.000000,1589.400000,505.610000,72.490000 /)
XPIZA_LKT(37,7,1:6)=(/ 0.950211,0.971711,0.980695,0.993238,0.991062,0.976262 /)
XCGA_LKT(37,7,1:6)=(/ 0.724123,0.731890,0.737463,0.714283,0.606163,0.319613 /)
XEXT_COEFF_550_LKT(37,7)=2427.500000 !rg=0.185845 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(37,8,1:6)=(/ 1916.200000,2127.500000,2101.300000,1580.400000,624.340000,108.930000 /)
XPIZA_LKT(37,8,1:6)=(/ 0.941420,0.965620,0.977490,0.992963,0.992304,0.983083 /)
XCGA_LKT(37,8,1:6)=(/ 0.734467,0.728130,0.733930,0.726210,0.661207,0.420557 /)
XEXT_COEFF_550_LKT(37,8)=2105.000000 !rg=0.185845 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(37,9,1:6)=(/ 1570.900000,1737.800000,1777.200000,1481.800000,726.450000,155.390000 /)
XPIZA_LKT(37,9,1:6)=(/ 0.933951,0.959269,0.973458,0.992255,0.993021,0.987272 /)
XCGA_LKT(37,9,1:6)=(/ 0.745880,0.729663,0.731630,0.728863,0.700580,0.506277 /)
XEXT_COEFF_550_LKT(37,9)=1778.300000 !rg=0.185845 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(37,10,1:6)=(/ 1300.900000,1420.400000,1479.600000,1332.200000,790.860000,210.730000 /)
XPIZA_LKT(37,10,1:6)=(/ 0.922504,0.951682,0.967899,0.991201,0.993292,0.989912 /)
XCGA_LKT(37,10,1:6)=(/ 0.757257,0.735827,0.728817,0.728710,0.725930,0.581700 /)
XEXT_COEFF_550_LKT(37,10)=1481.900000 !rg=0.185845 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(37,11,1:6)=(/ 1071.600000,1157.600000,1224.000000,1164.200000,803.730000,270.380000 /)
XPIZA_LKT(37,11,1:6)=(/ 0.910827,0.942350,0.961518,0.989754,0.993141,0.991593 /)
XCGA_LKT(37,11,1:6)=(/ 0.771133,0.740577,0.731077,0.728467,0.739020,0.640113 /)
XEXT_COEFF_550_LKT(37,11)=1217.500000 !rg=0.185845 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(37,12,1:6)=(/ 879.280000,942.650000,998.030000,989.790000,764.190000,327.320000 /)
XPIZA_LKT(37,12,1:6)=(/ 0.897560,0.932402,0.954416,0.987843,0.992567,0.992626 /)
XCGA_LKT(37,12,1:6)=(/ 0.781833,0.752397,0.733003,0.728580,0.742297,0.685190 /)
XEXT_COEFF_550_LKT(37,12)=997.490000 !rg=0.185845 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(37,13,1:6)=(/ 718.570000,766.760000,813.910000,825.780000,694.180000,371.850000 /)
XPIZA_LKT(37,13,1:6)=(/ 0.883193,0.921271,0.946025,0.985569,0.991583,0.993161 /)
XCGA_LKT(37,13,1:6)=(/ 0.796610,0.762553,0.741730,0.728930,0.740717,0.716547 /)
XEXT_COEFF_550_LKT(37,13)=811.560000 !rg=0.185845 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(37,14,1:6)=(/ 593.780000,625.490000,661.290000,683.830000,616.140000,393.570000 /)
XPIZA_LKT(37,14,1:6)=(/ 0.867750,0.909383,0.937386,0.982882,0.990327,0.993250 /)
XCGA_LKT(37,14,1:6)=(/ 0.807947,0.773190,0.749687,0.731010,0.741150,0.734947 /)
XEXT_COEFF_550_LKT(37,14)=663.980000 !rg=0.185845 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(37,15,1:6)=(/ 488.730000,514.090000,539.030000,565.590000,535.090000,387.370000 /)
XPIZA_LKT(37,15,1:6)=(/ 0.853028,0.893968,0.926699,0.978805,0.988828,0.992890 /)
XCGA_LKT(37,15,1:6)=(/ 0.816577,0.785780,0.757903,0.734877,0.743690,0.741373 /)
XEXT_COEFF_550_LKT(37,15)=538.390000 !rg=0.185845 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(37,16,1:6)=(/ 401.240000,418.210000,440.670000,459.830000,452.810000,358.650000 /)
XPIZA_LKT(37,16,1:6)=(/ 0.840078,0.878407,0.913824,0.975235,0.986511,0.992069 /)
XCGA_LKT(37,16,1:6)=(/ 0.828943,0.796183,0.769513,0.740197,0.743930,0.739047 /)
XEXT_COEFF_550_LKT(37,16)=437.540000 !rg=0.185845 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(37,17,1:6)=(/ 330.350000,342.370000,358.380000,375.270000,379.820000,320.550000 /)
XPIZA_LKT(37,17,1:6)=(/ 0.828436,0.864327,0.900177,0.971001,0.984373,0.990953 /)
XCGA_LKT(37,17,1:6)=(/ 0.838203,0.808733,0.780357,0.748017,0.747340,0.736260 /)
XEXT_COEFF_550_LKT(37,17)=359.290000 !rg=0.185845 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(37,18,1:6)=(/ 272.480000,282.500000,293.130000,308.850000,315.280000,284.710000 /)
XPIZA_LKT(37,18,1:6)=(/ 0.817498,0.850142,0.885178,0.964926,0.981144,0.989626 /)
XCGA_LKT(37,18,1:6)=(/ 0.844407,0.819593,0.791183,0.755907,0.750743,0.738693 /)
XEXT_COEFF_550_LKT(37,18)=292.630000 !rg=0.185845 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(37,19,1:6)=(/ 224.640000,231.010000,240.400000,250.610000,261.170000,247.140000 /)
XPIZA_LKT(37,19,1:6)=(/ 0.807713,0.838363,0.869760,0.959222,0.977478,0.987796 /)
XCGA_LKT(37,19,1:6)=(/ 0.852793,0.828040,0.803293,0.763807,0.755300,0.741120 /)
XEXT_COEFF_550_LKT(37,19)=239.500000 !rg=0.185845 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(37,20,1:6)=(/ 185.320000,189.940000,196.550000,204.670000,214.080000,208.830000 /)
XPIZA_LKT(37,20,1:6)=(/ 0.797945,0.827278,0.855922,0.952060,0.974053,0.985797 /)
XCGA_LKT(37,20,1:6)=(/ 0.859117,0.837570,0.813753,0.773967,0.761157,0.743277 /)
XEXT_COEFF_550_LKT(37,20)=196.800000 !rg=0.185845 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(38,1,1:6)=(/ 5417.500000,5197.800000,2682.800000,757.760000,103.290000,7.790100 /)
XPIZA_LKT(38,1,1:6)=(/ 0.977965,0.986179,0.984664,0.989019,0.967117,0.815926 /)
XCGA_LKT(38,1,1:6)=(/ 0.747317,0.768917,0.656273,0.436650,0.125223,0.030920 /)
XEXT_COEFF_550_LKT(38,1)=2716.100000 !rg=0.201336 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(38,2,1:6)=(/ 4837.300000,5050.200000,2785.700000,859.310000,124.700000,9.441800 /)
XPIZA_LKT(38,2,1:6)=(/ 0.975075,0.985826,0.984688,0.989813,0.972217,0.846707 /)
XCGA_LKT(38,2,1:6)=(/ 0.732047,0.772660,0.689233,0.506293,0.158887,0.039207 /)
XEXT_COEFF_550_LKT(38,2)=2833.100000 !rg=0.201336 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(38,3,1:6)=(/ 4079.400000,4681.000000,2911.900000,1039.800000,167.050000,13.398000 /)
XPIZA_LKT(38,3,1:6)=(/ 0.970585,0.984571,0.985028,0.991120,0.978497,0.889917 /)
XCGA_LKT(38,3,1:6)=(/ 0.713590,0.769067,0.715263,0.572253,0.236470,0.059430 /)
XEXT_COEFF_550_LKT(38,3)=2956.500000 !rg=0.201336 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(38,4,1:6)=(/ 3452.800000,4112.100000,2955.900000,1244.500000,233.160000,21.280000 /)
XPIZA_LKT(38,4,1:6)=(/ 0.965331,0.982259,0.984941,0.992222,0.983719,0.928635 /)
XCGA_LKT(38,4,1:6)=(/ 0.706100,0.759073,0.732990,0.627643,0.354943,0.100000 /)
XEXT_COEFF_550_LKT(38,4)=2992.500000 !rg=0.201336 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(38,5,1:6)=(/ 2942.000000,3471.000000,2856.000000,1435.600000,325.570000,35.017000 /)
XPIZA_LKT(38,5,1:6)=(/ 0.959373,0.978873,0.984088,0.992946,0.987512,0.954867 /)
XCGA_LKT(38,5,1:6)=(/ 0.707527,0.746670,0.741550,0.672277,0.471467,0.165237 /)
XEXT_COEFF_550_LKT(38,5)=2882.000000 !rg=0.201336 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(38,6,1:6)=(/ 2476.400000,2882.000000,2619.700000,1570.200000,437.230000,55.911000 /)
XPIZA_LKT(38,6,1:6)=(/ 0.952871,0.974356,0.982350,0.993279,0.990032,0.970278 /)
XCGA_LKT(38,6,1:6)=(/ 0.715647,0.735750,0.742220,0.703503,0.561983,0.254110 /)
XEXT_COEFF_550_LKT(38,6)=2636.200000 !rg=0.201336 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(38,7,1:6)=(/ 2075.200000,2368.400000,2304.500000,1617.300000,558.530000,86.205000 /)
XPIZA_LKT(38,7,1:6)=(/ 0.945222,0.968723,0.979642,0.993229,0.991673,0.979467 /)
XCGA_LKT(38,7,1:6)=(/ 0.724897,0.728957,0.738560,0.722180,0.629640,0.358817 /)
XEXT_COEFF_550_LKT(38,7)=2313.000000 !rg=0.201336 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(38,8,1:6)=(/ 1727.600000,1938.700000,1971.200000,1564.400000,672.420000,126.870000 /)
XPIZA_LKT(38,8,1:6)=(/ 0.937216,0.962623,0.976020,0.992773,0.992674,0.985047 /)
XCGA_LKT(38,8,1:6)=(/ 0.737780,0.728653,0.734233,0.729790,0.678443,0.454363 /)
XEXT_COEFF_550_LKT(38,8)=1977.000000 !rg=0.201336 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(38,9,1:6)=(/ 1436.700000,1579.200000,1649.100000,1434.200000,761.370000,177.340000 /)
XPIZA_LKT(38,9,1:6)=(/ 0.927949,0.955610,0.971322,0.991922,0.993200,0.988496 /)
XCGA_LKT(38,9,1:6)=(/ 0.751593,0.729953,0.730237,0.730593,0.712337,0.537727 /)
XEXT_COEFF_550_LKT(38,9)=1651.800000 !rg=0.201336 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(38,10,1:6)=(/ 1180.600000,1291.100000,1365.400000,1272.000000,805.410000,234.910000 /)
XPIZA_LKT(38,10,1:6)=(/ 0.918199,0.947742,0.965811,0.990693,0.993296,0.990693 /)
XCGA_LKT(38,10,1:6)=(/ 0.764370,0.737243,0.730913,0.729180,0.732643,0.606173 /)
XEXT_COEFF_550_LKT(38,10)=1364.800000 !rg=0.201336 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(38,11,1:6)=(/ 975.700000,1050.800000,1119.300000,1096.500000,795.690000,294.690000 /)
XPIZA_LKT(38,11,1:6)=(/ 0.905188,0.938991,0.958855,0.989120,0.992976,0.992088 /)
XCGA_LKT(38,11,1:6)=(/ 0.777093,0.747013,0.731223,0.729380,0.741570,0.659603 /)
XEXT_COEFF_550_LKT(38,11)=1120.900000 !rg=0.201336 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(38,12,1:6)=(/ 800.300000,857.890000,914.120000,920.240000,740.870000,347.800000 /)
XPIZA_LKT(38,12,1:6)=(/ 0.891624,0.927925,0.951341,0.987017,0.992221,0.992903 /)
XCGA_LKT(38,12,1:6)=(/ 0.789910,0.756530,0.738137,0.728070,0.742180,0.699023 /)
XEXT_COEFF_550_LKT(38,12)=912.590000 !rg=0.201336 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(38,13,1:6)=(/ 660.860000,698.890000,742.290000,765.650000,664.150000,384.050000 /)
XPIZA_LKT(38,13,1:6)=(/ 0.875582,0.916324,0.942550,0.984524,0.991157,0.993253 /)
XCGA_LKT(38,13,1:6)=(/ 0.802260,0.766953,0.744463,0.729873,0.741377,0.725213 /)
XEXT_COEFF_550_LKT(38,13)=745.250000 !rg=0.201336 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(38,14,1:6)=(/ 544.520000,574.420000,604.000000,632.580000,585.390000,394.590000 /)
XPIZA_LKT(38,14,1:6)=(/ 0.860179,0.902277,0.933135,0.981100,0.989807,0.993158 /)
XCGA_LKT(38,14,1:6)=(/ 0.811257,0.779437,0.752750,0.732180,0.743087,0.738587 /)
XEXT_COEFF_550_LKT(38,14)=603.980000 !rg=0.201336 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(38,15,1:6)=(/ 448.360000,468.940000,495.210000,517.910000,501.410000,378.080000 /)
XPIZA_LKT(38,15,1:6)=(/ 0.846631,0.886994,0.920848,0.977727,0.987968,0.992607 /)
XCGA_LKT(38,15,1:6)=(/ 0.823367,0.788873,0.763510,0.735480,0.743767,0.740967 /)
XEXT_COEFF_550_LKT(38,15)=490.820000 !rg=0.201336 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(38,16,1:6)=(/ 367.960000,381.590000,401.480000,420.670000,421.950000,343.530000 /)
XPIZA_LKT(38,16,1:6)=(/ 0.834317,0.872002,0.908330,0.973944,0.985521,0.991655 /)
XCGA_LKT(38,16,1:6)=(/ 0.833780,0.804050,0.772897,0.745423,0.744053,0.737873 /)
XEXT_COEFF_550_LKT(38,16)=402.330000 !rg=0.201336 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(38,17,1:6)=(/ 303.780000,314.780000,327.290000,346.050000,351.600000,306.280000 /)
XPIZA_LKT(38,17,1:6)=(/ 0.822976,0.857890,0.894295,0.968324,0.983379,0.990435 /)
XCGA_LKT(38,17,1:6)=(/ 0.840997,0.814450,0.785680,0.751663,0.750223,0.737127 /)
XEXT_COEFF_550_LKT(38,17)=327.360000 !rg=0.201336 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(38,18,1:6)=(/ 250.650000,258.630000,269.370000,281.650000,291.970000,270.260000 /)
XPIZA_LKT(38,18,1:6)=(/ 0.812853,0.844604,0.877675,0.962986,0.979816,0.988903 /)
XCGA_LKT(38,18,1:6)=(/ 0.848997,0.822783,0.797477,0.757437,0.752633,0.739773 /)
XEXT_COEFF_550_LKT(38,18)=267.150000 !rg=0.201336 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(38,19,1:6)=(/ 206.430000,211.710000,219.480000,228.610000,240.180000,230.800000 /)
XPIZA_LKT(38,19,1:6)=(/ 0.802956,0.832888,0.863968,0.956719,0.976310,0.986923 /)
XCGA_LKT(38,19,1:6)=(/ 0.856237,0.834227,0.807113,0.770163,0.756160,0.740860 /)
XEXT_COEFF_550_LKT(38,19)=220.070000 !rg=0.201336 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(38,20,1:6)=(/ 170.410000,174.670000,179.840000,188.190000,196.350000,195.100000 /)
XPIZA_LKT(38,20,1:6)=(/ 0.793013,0.822004,0.851272,0.948015,0.972246,0.984962 /)
XCGA_LKT(38,20,1:6)=(/ 0.861180,0.842040,0.818360,0.778647,0.765717,0.746967 /)
XEXT_COEFF_550_LKT(38,20)=179.750000 !rg=0.201336 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(39,1,1:6)=(/ 4519.100000,5192.000000,2823.400000,873.110000,127.930000,9.522700 /)
XPIZA_LKT(39,1,1:6)=(/ 0.973137,0.986170,0.985078,0.989893,0.972716,0.847713 /)
XCGA_LKT(39,1,1:6)=(/ 0.721267,0.781630,0.688190,0.514497,0.147510,0.036250 /)
XEXT_COEFF_550_LKT(39,1)=2873.500000 !rg=0.218119 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(39,2,1:6)=(/ 4009.500000,4928.600000,2962.000000,997.490000,152.670000,11.610000 /)
XPIZA_LKT(39,2,1:6)=(/ 0.969949,0.985324,0.985298,0.990736,0.976676,0.873765 /)
XCGA_LKT(39,2,1:6)=(/ 0.703360,0.777287,0.712743,0.561683,0.187290,0.045953 /)
XEXT_COEFF_550_LKT(39,2)=3009.700000 !rg=0.218119 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(39,3,1:6)=(/ 3477.600000,4456.200000,3038.000000,1177.600000,200.220000,16.575000 /)
XPIZA_LKT(39,3,1:6)=(/ 0.965214,0.983649,0.985399,0.991847,0.981521,0.909692 /)
XCGA_LKT(39,3,1:6)=(/ 0.692357,0.769267,0.731343,0.606840,0.276077,0.069590 /)
XEXT_COEFF_550_LKT(39,3)=3079.700000 !rg=0.218119 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(39,4,1:6)=(/ 3024.100000,3821.800000,3003.300000,1368.900000,274.050000,26.274000 /)
XPIZA_LKT(39,4,1:6)=(/ 0.960283,0.980810,0.984954,0.992708,0.985677,0.941161 /)
XCGA_LKT(39,4,1:6)=(/ 0.696130,0.755637,0.743250,0.652413,0.398420,0.116740 /)
XEXT_COEFF_550_LKT(39,4)=3035.500000 !rg=0.218119 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(39,5,1:6)=(/ 2603.900000,3184.400000,2819.600000,1533.100000,374.380000,42.662000 /)
XPIZA_LKT(39,5,1:6)=(/ 0.954568,0.976752,0.983692,0.993226,0.988766,0.962165 /)
XCGA_LKT(39,5,1:6)=(/ 0.705393,0.741177,0.746823,0.689763,0.508117,0.191960 /)
XEXT_COEFF_550_LKT(39,5)=2840.700000 !rg=0.218119 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(39,6,1:6)=(/ 2219.000000,2617.500000,2524.100000,1627.300000,490.900000,67.075000 /)
XPIZA_LKT(39,6,1:6)=(/ 0.949053,0.971881,0.981475,0.993379,0.990837,0.974582 /)
XCGA_LKT(39,6,1:6)=(/ 0.719957,0.731150,0.743770,0.714867,0.590427,0.290720 /)
XEXT_COEFF_550_LKT(39,6)=2534.300000 !rg=0.218119 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(39,7,1:6)=(/ 1874.900000,2144.500000,2180.300000,1626.600000,611.320000,101.780000 /)
XPIZA_LKT(39,7,1:6)=(/ 0.940472,0.965847,0.978288,0.993156,0.992174,0.982099 /)
XCGA_LKT(39,7,1:6)=(/ 0.732393,0.726453,0.738093,0.728413,0.650920,0.395793 /)
XEXT_COEFF_550_LKT(39,7)=2183.500000 !rg=0.218119 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(39,8,1:6)=(/ 1560.300000,1755.300000,1838.500000,1532.500000,716.810000,146.400000 /)
XPIZA_LKT(39,8,1:6)=(/ 0.933378,0.959374,0.974159,0.992523,0.992960,0.986648 /)
XCGA_LKT(39,8,1:6)=(/ 0.745500,0.727413,0.733417,0.732183,0.693707,0.487113 /)
XEXT_COEFF_550_LKT(39,8)=1837.900000 !rg=0.218119 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(39,9,1:6)=(/ 1305.100000,1442.100000,1521.700000,1380.900000,789.370000,200.790000 /)
XPIZA_LKT(39,9,1:6)=(/ 0.921832,0.951007,0.969053,0.991485,0.993305,0.989527 /)
XCGA_LKT(39,9,1:6)=(/ 0.753873,0.732007,0.729940,0.731310,0.722197,0.566923 /)
XEXT_COEFF_550_LKT(39,9)=1520.800000 !rg=0.218119 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(39,10,1:6)=(/ 1082.900000,1173.100000,1254.400000,1206.700000,811.220000,259.670000 /)
XPIZA_LKT(39,10,1:6)=(/ 0.911527,0.943478,0.963158,0.990118,0.993233,0.991338 /)
XCGA_LKT(39,10,1:6)=(/ 0.770527,0.739267,0.730600,0.730530,0.737730,0.628617 /)
XEXT_COEFF_550_LKT(39,10)=1257.800000 !rg=0.218119 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(39,11,1:6)=(/ 890.140000,959.710000,1022.800000,1026.300000,780.880000,317.770000 /)
XPIZA_LKT(39,11,1:6)=(/ 0.899172,0.933650,0.956402,0.988245,0.992734,0.992475 /)
XCGA_LKT(39,11,1:6)=(/ 0.782223,0.750027,0.734590,0.728477,0.742930,0.676707 /)
XEXT_COEFF_550_LKT(39,11)=1024.400000 !rg=0.218119 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(39,12,1:6)=(/ 735.950000,781.680000,833.950000,854.670000,713.310000,365.320000 /)
XPIZA_LKT(39,12,1:6)=(/ 0.883769,0.922810,0.947867,0.986033,0.991845,0.993095 /)
XCGA_LKT(39,12,1:6)=(/ 0.795967,0.760213,0.739927,0.728590,0.742203,0.710767 /)
XEXT_COEFF_550_LKT(39,12)=837.220000 !rg=0.218119 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(39,13,1:6)=(/ 605.890000,641.590000,677.720000,709.420000,634.700000,392.160000 /)
XPIZA_LKT(39,13,1:6)=(/ 0.867541,0.909685,0.938739,0.982944,0.990665,0.993265 /)
XCGA_LKT(39,13,1:6)=(/ 0.805643,0.773033,0.747600,0.730513,0.742723,0.731870 /)
XEXT_COEFF_550_LKT(39,13)=677.730000 !rg=0.218119 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(39,14,1:6)=(/ 499.740000,524.290000,555.220000,581.520000,552.870000,391.210000 /)
XPIZA_LKT(39,14,1:6)=(/ 0.853661,0.895359,0.927346,0.979673,0.989096,0.992990 /)
XCGA_LKT(39,14,1:6)=(/ 0.818500,0.782030,0.757640,0.731820,0.743697,0.740637 /)
XEXT_COEFF_550_LKT(39,14)=549.900000 !rg=0.218119 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(39,15,1:6)=(/ 410.590000,428.090000,451.230000,473.470000,467.940000,365.590000 /)
XPIZA_LKT(39,15,1:6)=(/ 0.839989,0.879672,0.915588,0.976394,0.986989,0.992278 /)
XCGA_LKT(39,15,1:6)=(/ 0.827977,0.797673,0.765593,0.741607,0.743000,0.740267 /)
XEXT_COEFF_550_LKT(39,15)=451.080000 !rg=0.218119 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(39,16,1:6)=(/ 338.110000,350.980000,365.760000,387.290000,393.130000,328.850000 /)
XPIZA_LKT(39,16,1:6)=(/ 0.828986,0.865170,0.902981,0.971433,0.984479,0.991191 /)
XCGA_LKT(39,16,1:6)=(/ 0.837350,0.809100,0.780490,0.746867,0.747693,0.737383 /)
XEXT_COEFF_550_LKT(39,16)=367.240000 !rg=0.218119 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(39,17,1:6)=(/ 279.300000,288.550000,300.790000,316.530000,324.610000,292.210000 /)
XPIZA_LKT(39,17,1:6)=(/ 0.817719,0.851111,0.886779,0.965834,0.981961,0.989852 /)
XCGA_LKT(39,17,1:6)=(/ 0.846047,0.817340,0.791477,0.751700,0.750827,0.738050 /)
XEXT_COEFF_550_LKT(39,17)=298.460000 !rg=0.218119 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(39,18,1:6)=(/ 229.930000,237.000000,246.140000,257.240000,269.370000,254.170000 /)
XPIZA_LKT(39,18,1:6)=(/ 0.807232,0.838762,0.871503,0.960318,0.978135,0.988240 /)
XCGA_LKT(39,18,1:6)=(/ 0.852203,0.829410,0.800510,0.765157,0.752097,0.741130 /)
XEXT_COEFF_550_LKT(39,18)=245.800000 !rg=0.218119 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(39,19,1:6)=(/ 189.960000,194.970000,200.640000,210.800000,221.120000,215.510000 /)
XPIZA_LKT(39,19,1:6)=(/ 0.798421,0.827426,0.858157,0.952875,0.974386,0.986039 /)
XCGA_LKT(39,19,1:6)=(/ 0.858617,0.838363,0.814083,0.773300,0.760793,0.744197 /)
XEXT_COEFF_550_LKT(39,19)=201.510000 !rg=0.218119 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(39,20,1:6)=(/ 156.990000,160.620000,165.320000,172.440000,180.330000,181.720000 /)
XPIZA_LKT(39,20,1:6)=(/ 0.789013,0.816667,0.845377,0.944353,0.969958,0.983881 /)
XCGA_LKT(39,20,1:6)=(/ 0.864487,0.844490,0.823560,0.779707,0.768553,0.748413 /)
XEXT_COEFF_550_LKT(39,20)=164.470000 !rg=0.218119 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(40,1,1:6)=(/ 3571.600000,4970.700000,3055.000000,1018.600000,157.400000,11.717000 /)
XPIZA_LKT(40,1,1:6)=(/ 0.966341,0.985433,0.985377,0.990703,0.977191,0.874630 /)
XCGA_LKT(40,1,1:6)=(/ 0.677150,0.784697,0.728597,0.584770,0.174190,0.042497 /)
XEXT_COEFF_550_LKT(40,1)=3120.700000 !rg=0.236301 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(40,2,1:6)=(/ 3277.900000,4704.800000,3114.100000,1151.200000,185.160000,14.350000 /)
XPIZA_LKT(40,2,1:6)=(/ 0.963035,0.984485,0.985727,0.991573,0.980216,0.896401 /)
XCGA_LKT(40,2,1:6)=(/ 0.669207,0.778233,0.733803,0.603440,0.221230,0.053860 /)
XEXT_COEFF_550_LKT(40,2)=3159.800000 !rg=0.236301 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(40,3,1:6)=(/ 2945.300000,4152.000000,3124.300000,1317.000000,237.890000,20.553000 /)
XPIZA_LKT(40,3,1:6)=(/ 0.959099,0.982314,0.985577,0.992460,0.983946,0.925949 /)
XCGA_LKT(40,3,1:6)=(/ 0.672113,0.765927,0.744923,0.635633,0.320887,0.081483 /)
XEXT_COEFF_550_LKT(40,3)=3162.100000 !rg=0.236301 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(40,4,1:6)=(/ 2645.600000,3494.300000,3006.000000,1487.600000,319.730000,32.373000 /)
XPIZA_LKT(40,4,1:6)=(/ 0.954981,0.978876,0.984767,0.993103,0.987285,0.951306 /)
XCGA_LKT(40,4,1:6)=(/ 0.688960,0.749540,0.751163,0.674150,0.442340,0.136170 /)
XEXT_COEFF_550_LKT(40,4)=3033.100000 !rg=0.236301 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(40,5,1:6)=(/ 2331.600000,2878.500000,2745.600000,1617.700000,427.010000,51.734000 /)
XPIZA_LKT(40,5,1:6)=(/ 0.949306,0.974327,0.983049,0.993433,0.989809,0.968072 /)
XCGA_LKT(40,5,1:6)=(/ 0.707380,0.735097,0.749940,0.705003,0.542803,0.222073 /)
XEXT_COEFF_550_LKT(40,5)=2760.700000 !rg=0.236301 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(40,6,1:6)=(/ 2000.700000,2365.700000,2401.100000,1666.100000,545.960000,80.063000 /)
XPIZA_LKT(40,6,1:6)=(/ 0.943770,0.969026,0.980394,0.993412,0.991506,0.978115 /)
XCGA_LKT(40,6,1:6)=(/ 0.722107,0.727307,0.743917,0.724380,0.616493,0.328087 /)
XEXT_COEFF_550_LKT(40,6)=2410.600000 !rg=0.236301 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(40,7,1:6)=(/ 1691.500000,1938.500000,2038.300000,1617.800000,662.450000,119.060000 /)
XPIZA_LKT(40,7,1:6)=(/ 0.936031,0.962927,0.976739,0.993011,0.992588,0.984242 /)
XCGA_LKT(40,7,1:6)=(/ 0.737117,0.725817,0.737203,0.733243,0.669980,0.430727 /)
XEXT_COEFF_550_LKT(40,7)=2044.200000 !rg=0.236301 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(40,8,1:6)=(/ 1425.200000,1588.600000,1698.900000,1488.100000,755.990000,167.810000 /)
XPIZA_LKT(40,8,1:6)=(/ 0.927237,0.955810,0.972097,0.992217,0.993176,0.987983 /)
XCGA_LKT(40,8,1:6)=(/ 0.750703,0.727727,0.731617,0.734337,0.706950,0.519493 /)
XEXT_COEFF_550_LKT(40,8)=1702.800000 !rg=0.236301 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(40,9,1:6)=(/ 1190.700000,1304.400000,1404.500000,1317.500000,808.960000,224.910000 /)
XPIZA_LKT(40,9,1:6)=(/ 0.915584,0.946962,0.966056,0.991030,0.993347,0.990377 /)
XCGA_LKT(40,9,1:6)=(/ 0.763920,0.732653,0.729467,0.732427,0.730227,0.592727 /)
XEXT_COEFF_550_LKT(40,9)=1398.100000 !rg=0.236301 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(40,10,1:6)=(/ 988.070000,1073.300000,1147.500000,1138.700000,808.290000,284.470000 /)
XPIZA_LKT(40,10,1:6)=(/ 0.904598,0.937936,0.960135,0.989353,0.993104,0.991880 /)
XCGA_LKT(40,10,1:6)=(/ 0.773497,0.743707,0.730820,0.730087,0.741400,0.649233 /)
XEXT_COEFF_550_LKT(40,10)=1146.300000 !rg=0.236301 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(40,11,1:6)=(/ 817.830000,874.190000,937.620000,953.890000,759.470000,339.470000 /)
XPIZA_LKT(40,11,1:6)=(/ 0.891389,0.927856,0.952455,0.987350,0.992440,0.992786 /)
XCGA_LKT(40,11,1:6)=(/ 0.789747,0.752547,0.735573,0.728407,0.743450,0.691663 /)
XEXT_COEFF_550_LKT(40,11)=937.270000 !rg=0.236301 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(40,12,1:6)=(/ 673.900000,717.220000,760.800000,795.030000,684.920000,379.580000 /)
XPIZA_LKT(40,12,1:6)=(/ 0.875662,0.916253,0.944024,0.984587,0.991424,0.993219 /)
XCGA_LKT(40,12,1:6)=(/ 0.799327,0.766363,0.742150,0.729220,0.742960,0.720600 /)
XEXT_COEFF_550_LKT(40,12)=760.660000 !rg=0.236301 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(40,13,1:6)=(/ 555.810000,585.200000,622.040000,651.600000,604.300000,395.500000 /)
XPIZA_LKT(40,13,1:6)=(/ 0.861153,0.903082,0.933351,0.981542,0.990048,0.993210 /)
XCGA_LKT(40,13,1:6)=(/ 0.813370,0.775380,0.751740,0.729350,0.743163,0.736627 /)
XEXT_COEFF_550_LKT(40,13)=616.350000 !rg=0.236301 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(40,14,1:6)=(/ 456.490000,479.160000,506.020000,533.630000,518.210000,383.600000 /)
XPIZA_LKT(40,14,1:6)=(/ 0.847554,0.887891,0.922095,0.978131,0.988369,0.992759 /)
XCGA_LKT(40,14,1:6)=(/ 0.822620,0.790830,0.759207,0.737497,0.744057,0.741373 /)
XEXT_COEFF_550_LKT(40,14)=504.940000 !rg=0.236301 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(40,15,1:6)=(/ 375.710000,392.270000,411.520000,434.450000,437.480000,352.000000 /)
XPIZA_LKT(40,15,1:6)=(/ 0.834887,0.873212,0.910089,0.974397,0.985686,0.991867 /)
XCGA_LKT(40,15,1:6)=(/ 0.833477,0.803313,0.774047,0.743580,0.744657,0.738930 /)
XEXT_COEFF_550_LKT(40,15)=412.370000 !rg=0.236301 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(40,16,1:6)=(/ 310.540000,322.030000,336.120000,355.020000,362.890000,314.000000 /)
XPIZA_LKT(40,16,1:6)=(/ 0.822788,0.857813,0.895705,0.968720,0.983348,0.990702 /)
XCGA_LKT(40,16,1:6)=(/ 0.842180,0.812280,0.785270,0.747563,0.749190,0.737410 /)
XEXT_COEFF_550_LKT(40,16)=335.600000 !rg=0.236301 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(40,17,1:6)=(/ 255.830000,264.680000,275.240000,290.270000,300.720000,277.250000 /)
XPIZA_LKT(40,17,1:6)=(/ 0.813000,0.845064,0.879186,0.962815,0.980085,0.989302 /)
XCGA_LKT(40,17,1:6)=(/ 0.848773,0.823820,0.794263,0.758087,0.750270,0.740597 /)
XEXT_COEFF_550_LKT(40,17)=274.330000 !rg=0.236301 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(40,18,1:6)=(/ 210.940000,217.560000,225.080000,235.630000,248.700000,238.730000 /)
XPIZA_LKT(40,18,1:6)=(/ 0.802973,0.833395,0.864939,0.957354,0.975765,0.986966 /)
XCGA_LKT(40,18,1:6)=(/ 0.856143,0.834063,0.808253,0.769153,0.755377,0.741500 /)
XEXT_COEFF_550_LKT(40,18)=225.240000 !rg=0.236301 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(40,19,1:6)=(/ 174.810000,179.080000,184.700000,193.230000,201.890000,200.660000 /)
XPIZA_LKT(40,19,1:6)=(/ 0.793548,0.821861,0.852025,0.948766,0.972398,0.985330 /)
XCGA_LKT(40,19,1:6)=(/ 0.861863,0.840910,0.818673,0.774780,0.764360,0.746123 /)
XEXT_COEFF_550_LKT(40,19)=184.540000 !rg=0.236301 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(40,20,1:6)=(/ 144.110000,147.460000,151.770000,158.140000,166.240000,168.820000 /)
XPIZA_LKT(40,20,1:6)=(/ 0.784604,0.811720,0.839350,0.940455,0.966755,0.982084 /)
XCGA_LKT(40,20,1:6)=(/ 0.866417,0.848907,0.826240,0.786290,0.767627,0.747157 /)
XEXT_COEFF_550_LKT(40,20)=151.480000 !rg=0.236301 sigma=2.95 wvl=0.55
 
IF (LHOOK) CALL DR_HOOK('MODE_DUSTOPT:DUST_OPT_LKT_SET4',1,ZHOOK_HANDLE)
END SUBROUTINE DUST_OPT_LKT_SET4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE DUST_OPT_LKT_SET5()

  USE MODD_DUST_OPT_LKT
  
  IMPLICIT NONE

REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_DUSTOPT:DUST_OPT_LKT_SET5',0,ZHOOK_HANDLE)
XEXT_COEFF_WVL_LKT(41,1,1:6)=(/ 2811.300000,4764.400000,3267.400000,1203.300000,191.860000,14.492000 /)
XPIZA_LKT(41,1,1:6)=(/ 0.956493,0.984454,0.986046,0.991591,0.980746,0.897146 /)
XCGA_LKT(41,1,1:6)=(/ 0.622593,0.782433,0.751803,0.629970,0.206360,0.049820 /)
XEXT_COEFF_550_LKT(41,1)=3310.300000 !rg=0.255998 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(41,2,1:6)=(/ 2670.200000,4381.100000,3226.900000,1309.800000,222.060000,17.803000 /)
XPIZA_LKT(41,2,1:6)=(/ 0.954477,0.983208,0.985997,0.992320,0.983012,0.915145 /)
XCGA_LKT(41,2,1:6)=(/ 0.630387,0.774830,0.749910,0.631720,0.261817,0.063130 /)
XEXT_COEFF_550_LKT(41,2)=3268.200000 !rg=0.255998 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(41,3,1:6)=(/ 2542.400000,3782.900000,3164.200000,1452.800000,280.330000,25.505000 /)
XPIZA_LKT(41,3,1:6)=(/ 0.952794,0.980453,0.985558,0.992967,0.985902,0.939216 /)
XCGA_LKT(41,3,1:6)=(/ 0.662367,0.758537,0.755677,0.660133,0.370117,0.095397 /)
XEXT_COEFF_550_LKT(41,3)=3196.900000 !rg=0.255998 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(41,4,1:6)=(/ 2347.200000,3142.500000,2961.800000,1596.500000,370.220000,39.751000 /)
XPIZA_LKT(41,4,1:6)=(/ 0.949773,0.976365,0.984356,0.993415,0.988615,0.959493 /)
XCGA_LKT(41,4,1:6)=(/ 0.689660,0.740133,0.756627,0.693250,0.485330,0.158643 /)
XEXT_COEFF_550_LKT(41,4)=2982.800000 !rg=0.255998 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(41,5,1:6)=(/ 2084.000000,2584.700000,2633.900000,1685.800000,482.840000,62.422000 /)
XPIZA_LKT(41,5,1:6)=(/ 0.944643,0.971405,0.982193,0.993568,0.990681,0.972871 /)
XCGA_LKT(41,5,1:6)=(/ 0.710640,0.728313,0.751260,0.718080,0.575103,0.255150 /)
XEXT_COEFF_550_LKT(41,5)=2646.000000 !rg=0.255998 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(41,6,1:6)=(/ 1805.400000,2132.900000,2259.400000,1685.200000,601.610000,94.905000 /)
XPIZA_LKT(41,6,1:6)=(/ 0.938461,0.965351,0.979035,0.993373,0.992062,0.981003 /)
XCGA_LKT(41,6,1:6)=(/ 0.726397,0.721787,0.743137,0.732007,0.640230,0.364950 /)
XEXT_COEFF_550_LKT(41,6)=2262.400000 !rg=0.255998 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(41,7,1:6)=(/ 1532.100000,1750.600000,1891.700000,1591.600000,710.620000,138.090000 /)
XPIZA_LKT(41,7,1:6)=(/ 0.931269,0.958879,0.974900,0.992789,0.992915,0.985994 /)
XCGA_LKT(41,7,1:6)=(/ 0.742703,0.723097,0.735643,0.736477,0.686987,0.464867 /)
XEXT_COEFF_550_LKT(41,7)=1890.900000 !rg=0.255998 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(41,8,1:6)=(/ 1295.900000,1444.900000,1562.300000,1434.600000,788.770000,190.870000 /)
XPIZA_LKT(41,8,1:6)=(/ 0.920894,0.951102,0.969804,0.991803,0.993320,0.989103 /)
XCGA_LKT(41,8,1:6)=(/ 0.754197,0.728737,0.730610,0.735083,0.718247,0.549903 /)
XEXT_COEFF_550_LKT(41,8)=1561.000000 !rg=0.255998 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(41,9,1:6)=(/ 1081.200000,1179.800000,1283.600000,1250.500000,819.850000,249.750000 /)
XPIZA_LKT(41,9,1:6)=(/ 0.910143,0.943474,0.963282,0.990482,0.993323,0.991079 /)
XCGA_LKT(41,9,1:6)=(/ 0.770223,0.739120,0.728840,0.733170,0.736570,0.616357 /)
XEXT_COEFF_550_LKT(41,9)=1284.300000 !rg=0.255998 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(41,10,1:6)=(/ 904.590000,973.130000,1052.400000,1062.600000,797.340000,308.270000 /)
XPIZA_LKT(41,10,1:6)=(/ 0.897556,0.932998,0.956149,0.988603,0.992906,0.992313 /)
XCGA_LKT(41,10,1:6)=(/ 0.783013,0.745057,0.731747,0.730040,0.743597,0.667477 /)
XEXT_COEFF_550_LKT(41,10)=1045.100000 !rg=0.255998 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(41,11,1:6)=(/ 747.230000,799.550000,854.680000,888.030000,734.240000,358.460000 /)
XPIZA_LKT(41,11,1:6)=(/ 0.883796,0.922100,0.948264,0.986025,0.992102,0.993017 /)
XCGA_LKT(41,11,1:6)=(/ 0.792763,0.759577,0.736673,0.728560,0.744267,0.704417 /)
XEXT_COEFF_550_LKT(41,11)=852.100000 !rg=0.255998 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(41,12,1:6)=(/ 617.990000,653.020000,698.520000,731.910000,655.310000,389.760000 /)
XPIZA_LKT(41,12,1:6)=(/ 0.868452,0.910288,0.938791,0.983435,0.990950,0.993268 /)
XCGA_LKT(41,12,1:6)=(/ 0.807607,0.768553,0.746193,0.728530,0.743487,0.728197 /)
XEXT_COEFF_550_LKT(41,12)=692.100000 !rg=0.255998 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(41,13,1:6)=(/ 507.480000,534.930000,567.040000,599.400000,570.830000,394.320000 /)
XPIZA_LKT(41,13,1:6)=(/ 0.854553,0.896023,0.928204,0.979937,0.989484,0.993084 /)
XCGA_LKT(41,13,1:6)=(/ 0.817443,0.784180,0.753277,0.733840,0.744783,0.739677 /)
XEXT_COEFF_550_LKT(41,13)=565.150000 !rg=0.255998 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(41,14,1:6)=(/ 417.690000,437.610000,462.160000,487.230000,485.170000,373.100000 /)
XPIZA_LKT(41,14,1:6)=(/ 0.841351,0.880889,0.916406,0.976774,0.987133,0.992440 /)
XCGA_LKT(41,14,1:6)=(/ 0.829713,0.796570,0.767247,0.739463,0.744030,0.740783 /)
XEXT_COEFF_550_LKT(41,14)=461.600000 !rg=0.255998 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(41,15,1:6)=(/ 346.000000,359.530000,376.280000,398.080000,405.120000,337.100000 /)
XPIZA_LKT(41,15,1:6)=(/ 0.829480,0.866296,0.904187,0.971979,0.985053,0.991428 /)
XCGA_LKT(41,15,1:6)=(/ 0.838380,0.807103,0.778970,0.744393,0.747160,0.737980 /)
XEXT_COEFF_550_LKT(41,15)=377.870000 !rg=0.255998 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(41,16,1:6)=(/ 285.440000,295.630000,307.810000,325.670000,335.170000,299.660000 /)
XPIZA_LKT(41,16,1:6)=(/ 0.817676,0.852082,0.887672,0.966365,0.982454,0.990173 /)
XCGA_LKT(41,16,1:6)=(/ 0.845250,0.819057,0.788350,0.753110,0.749483,0.739140 /)
XEXT_COEFF_550_LKT(41,16)=306.510000 !rg=0.255998 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(41,17,1:6)=(/ 234.830000,242.480000,252.290000,264.450000,278.100000,262.160000 /)
XPIZA_LKT(41,17,1:6)=(/ 0.807202,0.838685,0.872321,0.960164,0.978534,0.988444 /)
XCGA_LKT(41,17,1:6)=(/ 0.853633,0.828313,0.801597,0.762340,0.753190,0.741350 /)
XEXT_COEFF_550_LKT(41,17)=251.900000 !rg=0.255998 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(41,18,1:6)=(/ 194.440000,199.920000,206.400000,216.180000,227.060000,221.560000 /)
XPIZA_LKT(41,18,1:6)=(/ 0.798813,0.827824,0.859045,0.953770,0.974890,0.986600 /)
XCGA_LKT(41,18,1:6)=(/ 0.859407,0.837193,0.813433,0.770903,0.759537,0.743460 /)
XEXT_COEFF_550_LKT(41,18)=206.870000 !rg=0.255998 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(41,19,1:6)=(/ 160.650000,164.690000,169.540000,177.320000,185.700000,187.300000 /)
XPIZA_LKT(41,19,1:6)=(/ 0.788909,0.817234,0.845650,0.945100,0.970629,0.984313 /)
XCGA_LKT(41,19,1:6)=(/ 0.864150,0.845670,0.821593,0.780850,0.765890,0.747043 /)
XEXT_COEFF_550_LKT(41,19)=169.040000 !rg=0.255998 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(41,20,1:6)=(/ 132.530000,135.680000,139.310000,144.370000,152.930000,156.590000 /)
XPIZA_LKT(41,20,1:6)=(/ 0.779260,0.806302,0.834012,0.937250,0.964464,0.980727 /)
XCGA_LKT(41,20,1:6)=(/ 0.869547,0.852437,0.831930,0.790570,0.772313,0.749530 /)
XEXT_COEFF_550_LKT(41,20)=139.270000 !rg=0.255998 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(42,1,1:6)=(/ 2180.300000,4346.400000,3315.000000,1401.000000,231.020000,17.996000 /)
XPIZA_LKT(42,1,1:6)=(/ 0.945101,0.983221,0.986467,0.992523,0.983543,0.915795 /)
XCGA_LKT(42,1,1:6)=(/ 0.572387,0.777217,0.757143,0.645573,0.245437,0.058410 /)
XEXT_COEFF_550_LKT(42,1)=3342.500000 !rg=0.277337 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(42,2,1:6)=(/ 2239.100000,3969.200000,3289.200000,1459.700000,263.130000,22.142000 /)
XPIZA_LKT(42,2,1:6)=(/ 0.946306,0.981346,0.986079,0.992947,0.985209,0.930537 /)
XCGA_LKT(42,2,1:6)=(/ 0.613353,0.766527,0.761540,0.652163,0.310050,0.074003 /)
XEXT_COEFF_550_LKT(42,2)=3325.300000 !rg=0.277337 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(42,3,1:6)=(/ 2261.700000,3367.400000,3151.900000,1580.200000,327.930000,31.626000 /)
XPIZA_LKT(42,3,1:6)=(/ 0.946426,0.977972,0.985322,0.993375,0.987494,0.949975 /)
XCGA_LKT(42,3,1:6)=(/ 0.658803,0.747450,0.763440,0.681813,0.422033,0.111690 /)
XEXT_COEFF_550_LKT(42,3)=3178.300000 !rg=0.277337 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(42,4,1:6)=(/ 2120.500000,2783.900000,2870.500000,1691.300000,425.330000,48.591000 /)
XPIZA_LKT(42,4,1:6)=(/ 0.944000,0.973340,0.983687,0.993648,0.989721,0.966092 /)
XCGA_LKT(42,4,1:6)=(/ 0.693583,0.730020,0.759467,0.709990,0.526157,0.184480 /)
XEXT_COEFF_550_LKT(42,4)=2885.500000 !rg=0.277337 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(42,5,1:6)=(/ 1879.900000,2302.500000,2492.900000,1734.400000,540.960000,74.884000 /)
XPIZA_LKT(42,5,1:6)=(/ 0.940170,0.967808,0.981024,0.993632,0.991409,0.976773 /)
XCGA_LKT(42,5,1:6)=(/ 0.718123,0.719663,0.750647,0.729073,0.604697,0.290413 /)
XEXT_COEFF_550_LKT(42,5)=2497.600000 !rg=0.277337 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(42,6,1:6)=(/ 1649.000000,1906.800000,2104.100000,1683.200000,656.340000,111.580000 /)
XPIZA_LKT(42,6,1:6)=(/ 0.932151,0.961572,0.977266,0.993269,0.992522,0.983357 /)
XCGA_LKT(42,6,1:6)=(/ 0.735410,0.718837,0.740147,0.738043,0.661583,0.401180 /)
XEXT_COEFF_550_LKT(42,6)=2103.500000 !rg=0.277337 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(42,7,1:6)=(/ 1403.500000,1575.000000,1741.500000,1548.800000,754.260000,159.070000 /)
XPIZA_LKT(42,7,1:6)=(/ 0.924832,0.954706,0.972585,0.992514,0.993170,0.987451 /)
XCGA_LKT(42,7,1:6)=(/ 0.751093,0.722717,0.732327,0.739087,0.701910,0.498567 /)
XEXT_COEFF_550_LKT(42,7)=1740.700000 !rg=0.277337 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(42,8,1:6)=(/ 1186.000000,1304.000000,1434.800000,1369.100000,813.660000,214.920000 /)
XPIZA_LKT(42,8,1:6)=(/ 0.914621,0.946360,0.966797,0.991370,0.993400,0.990029 /)
XCGA_LKT(42,8,1:6)=(/ 0.764483,0.729173,0.728673,0.736197,0.727717,0.577233 /)
XEXT_COEFF_550_LKT(42,8)=1428.000000 !rg=0.277337 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(42,9,1:6)=(/ 983.620000,1071.600000,1169.900000,1177.200000,821.950000,274.960000 /)
XPIZA_LKT(42,9,1:6)=(/ 0.905275,0.938627,0.960727,0.989779,0.993231,0.991669 /)
XCGA_LKT(42,9,1:6)=(/ 0.777440,0.742037,0.730777,0.732377,0.741347,0.638200 /)
XEXT_COEFF_550_LKT(42,9)=1168.800000 !rg=0.277337 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(42,10,1:6)=(/ 823.300000,883.650000,957.870000,987.920000,779.170000,331.050000 /)
XPIZA_LKT(42,10,1:6)=(/ 0.891051,0.928223,0.952439,0.987817,0.992665,0.992660 /)
XCGA_LKT(42,10,1:6)=(/ 0.788720,0.754050,0.731993,0.730673,0.745390,0.683607 /)
XEXT_COEFF_550_LKT(42,10)=957.030000 !rg=0.277337 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(42,11,1:6)=(/ 683.580000,726.040000,782.420000,818.640000,706.940000,374.650000 /)
XPIZA_LKT(42,11,1:6)=(/ 0.876098,0.916850,0.943538,0.985236,0.991684,0.993174 /)
XCGA_LKT(42,11,1:6)=(/ 0.802350,0.762763,0.740673,0.728877,0.744423,0.715397 /)
XEXT_COEFF_550_LKT(42,11)=775.760000 !rg=0.277337 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(42,12,1:6)=(/ 563.950000,595.650000,635.370000,671.860000,623.460000,395.530000 /)
XPIZA_LKT(42,12,1:6)=(/ 0.861744,0.903843,0.934046,0.982263,0.990450,0.993250 /)
XCGA_LKT(42,12,1:6)=(/ 0.812337,0.777960,0.747240,0.732197,0.745317,0.734080 /)
XEXT_COEFF_550_LKT(42,12)=634.160000 !rg=0.277337 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(42,13,1:6)=(/ 464.050000,487.790000,517.710000,548.000000,537.130000,389.020000 /)
XPIZA_LKT(42,13,1:6)=(/ 0.847927,0.889100,0.922637,0.978750,0.988562,0.992880 /)
XCGA_LKT(42,13,1:6)=(/ 0.825250,0.789857,0.760770,0.735930,0.744803,0.741107 /)
XEXT_COEFF_550_LKT(42,13)=516.220000 !rg=0.277337 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(42,14,1:6)=(/ 384.830000,400.200000,421.220000,446.390000,449.820000,359.790000 /)
XPIZA_LKT(42,14,1:6)=(/ 0.835073,0.874512,0.911074,0.975076,0.986658,0.992085 /)
XCGA_LKT(42,14,1:6)=(/ 0.834560,0.801513,0.772117,0.741083,0.745483,0.739937 /)
XEXT_COEFF_550_LKT(42,14)=422.710000 !rg=0.277337 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(42,15,1:6)=(/ 317.730000,330.360000,344.270000,366.490000,375.860000,322.620000 /)
XPIZA_LKT(42,15,1:6)=(/ 0.823086,0.859001,0.897311,0.969264,0.984039,0.990959 /)
XCGA_LKT(42,15,1:6)=(/ 0.841210,0.813640,0.782590,0.748447,0.748607,0.738460 /)
XEXT_COEFF_550_LKT(42,15)=343.870000 !rg=0.277337 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(42,16,1:6)=(/ 261.880000,270.030000,282.240000,296.280000,310.300000,285.450000 /)
XPIZA_LKT(42,16,1:6)=(/ 0.812596,0.846163,0.880756,0.964481,0.980814,0.989531 /)
XCGA_LKT(42,16,1:6)=(/ 0.850030,0.822637,0.795673,0.756123,0.751127,0.740213 /)
XEXT_COEFF_550_LKT(42,16)=280.610000 !rg=0.277337 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(42,17,1:6)=(/ 216.400000,221.980000,230.490000,241.410000,254.360000,245.110000 /)
XPIZA_LKT(42,17,1:6)=(/ 0.803328,0.833903,0.865577,0.958241,0.977742,0.987783 /)
XCGA_LKT(42,17,1:6)=(/ 0.856810,0.833303,0.806387,0.766413,0.755467,0.742243 /)
XEXT_COEFF_550_LKT(42,17)=230.960000 !rg=0.277337 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(42,18,1:6)=(/ 178.650000,183.830000,189.390000,198.670000,208.280000,207.130000 /)
XPIZA_LKT(42,18,1:6)=(/ 0.793711,0.822153,0.852891,0.949635,0.973320,0.985575 /)
XCGA_LKT(42,18,1:6)=(/ 0.861403,0.841883,0.816977,0.775773,0.762680,0.745100 /)
XEXT_COEFF_550_LKT(42,18)=188.910000 !rg=0.277337 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(42,19,1:6)=(/ 147.840000,150.870000,155.730000,161.420000,171.070000,174.230000 /)
XPIZA_LKT(42,19,1:6)=(/ 0.784899,0.812382,0.840300,0.942028,0.967905,0.982905 /)
XCGA_LKT(42,19,1:6)=(/ 0.867140,0.848483,0.827520,0.784760,0.768970,0.748543 /)
XEXT_COEFF_550_LKT(42,19)=155.360000 !rg=0.277337 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(42,20,1:6)=(/ 122.190000,124.380000,127.830000,132.190000,139.430000,144.230000 /)
XPIZA_LKT(42,20,1:6)=(/ 0.775612,0.802647,0.827911,0.934921,0.962821,0.980276 /)
XCGA_LKT(42,20,1:6)=(/ 0.871590,0.855690,0.836067,0.794753,0.776167,0.752537 /)
XEXT_COEFF_550_LKT(42,20)=127.860000 !rg=0.277337 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(43,1,1:6)=(/ 1877.200000,3906.600000,3347.400000,1553.700000,274.070000,22.411000 /)
XPIZA_LKT(43,1,1:6)=(/ 0.934644,0.980751,0.986251,0.993304,0.985715,0.931117 /)
XCGA_LKT(43,1,1:6)=(/ 0.549373,0.768077,0.763830,0.648107,0.293150,0.068490 /)
XEXT_COEFF_550_LKT(43,1)=3395.200000 !rg=0.300455 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(43,2,1:6)=(/ 2023.200000,3490.400000,3295.800000,1593.800000,308.370000,27.572000 /)
XPIZA_LKT(43,2,1:6)=(/ 0.939777,0.978675,0.985942,0.993431,0.986934,0.943092 /)
XCGA_LKT(43,2,1:6)=(/ 0.616407,0.752283,0.770097,0.671743,0.366207,0.086777 /)
XEXT_COEFF_550_LKT(43,2)=3325.700000 !rg=0.300455 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(43,3,1:6)=(/ 2045.500000,2940.900000,3083.700000,1695.300000,381.130000,39.123000 /)
XPIZA_LKT(43,3,1:6)=(/ 0.941701,0.974552,0.984834,0.993692,0.988806,0.958661 /)
XCGA_LKT(43,3,1:6)=(/ 0.674197,0.731173,0.768123,0.701550,0.474033,0.130780 /)
XEXT_COEFF_550_LKT(43,3)=3103.000000 !rg=0.300455 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(43,4,1:6)=(/ 1918.700000,2450.000000,2734.200000,1768.300000,484.520000,59.067000 /)
XPIZA_LKT(43,4,1:6)=(/ 0.938855,0.969409,0.982735,0.993806,0.990645,0.971408 /)
XCGA_LKT(43,4,1:6)=(/ 0.704350,0.716783,0.759827,0.724510,0.563823,0.213883 /)
XEXT_COEFF_550_LKT(43,4)=2743.000000 !rg=0.300455 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(43,5,1:6)=(/ 1719.700000,2041.000000,2325.800000,1760.900000,600.350000,89.226000 /)
XPIZA_LKT(43,5,1:6)=(/ 0.934314,0.963889,0.979501,0.993625,0.992017,0.979946 /)
XCGA_LKT(43,5,1:6)=(/ 0.727390,0.713840,0.747673,0.738073,0.631523,0.327177 /)
XEXT_COEFF_550_LKT(43,5)=2328.200000 !rg=0.300455 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(43,6,1:6)=(/ 1494.200000,1713.000000,1934.700000,1661.500000,708.890000,130.180000 /)
XPIZA_LKT(43,6,1:6)=(/ 0.927208,0.956982,0.975324,0.993080,0.992894,0.985287 /)
XCGA_LKT(43,6,1:6)=(/ 0.740757,0.717273,0.737363,0.742320,0.680720,0.437240 /)
XEXT_COEFF_550_LKT(43,6)=1934.700000 !rg=0.300455 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(43,7,1:6)=(/ 1276.900000,1426.800000,1590.500000,1494.000000,792.050000,181.830000 /)
XPIZA_LKT(43,7,1:6)=(/ 0.918745,0.949582,0.970090,0.992130,0.993353,0.988670 /)
XCGA_LKT(43,7,1:6)=(/ 0.754310,0.724430,0.730217,0.740190,0.714857,0.530493 /)
XEXT_COEFF_550_LKT(43,7)=1587.600000 !rg=0.300455 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(43,8,1:6)=(/ 1077.000000,1181.000000,1305.400000,1298.600000,830.000000,239.890000 /)
XPIZA_LKT(43,8,1:6)=(/ 0.909379,0.941645,0.963708,0.990824,0.993416,0.990797 /)
XCGA_LKT(43,8,1:6)=(/ 0.769933,0.735607,0.727473,0.736353,0.735457,0.602377 /)
XEXT_COEFF_550_LKT(43,8)=1304.100000 !rg=0.300455 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(43,9,1:6)=(/ 904.460000,972.300000,1064.900000,1098.800000,815.090000,299.510000 /)
XPIZA_LKT(43,9,1:6)=(/ 0.898160,0.933335,0.957115,0.989048,0.993080,0.992148 /)
XCGA_LKT(43,9,1:6)=(/ 0.785387,0.744910,0.730157,0.732130,0.744773,0.657707 /)
XEXT_COEFF_550_LKT(43,9)=1066.500000 !rg=0.300455 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(43,10,1:6)=(/ 749.930000,802.810000,871.250000,914.030000,756.750000,351.500000 /)
XPIZA_LKT(43,10,1:6)=(/ 0.885022,0.923279,0.949028,0.986868,0.992345,0.992931 /)
XCGA_LKT(43,10,1:6)=(/ 0.797257,0.758980,0.737143,0.729417,0.746310,0.697503 /)
XEXT_COEFF_550_LKT(43,10)=869.620000 !rg=0.300455 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(43,11,1:6)=(/ 624.960000,659.040000,710.460000,753.000000,676.490000,386.910000 /)
XPIZA_LKT(43,11,1:6)=(/ 0.868927,0.911966,0.939570,0.984369,0.991260,0.993261 /)
XCGA_LKT(43,11,1:6)=(/ 0.808653,0.772203,0.742547,0.731337,0.745693,0.724070 /)
XEXT_COEFF_550_LKT(43,11)=711.260000 !rg=0.300455 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(43,12,1:6)=(/ 515.150000,543.040000,579.140000,615.240000,591.290000,396.850000 /)
XPIZA_LKT(43,12,1:6)=(/ 0.855316,0.897774,0.929107,0.981030,0.989766,0.993160 /)
XCGA_LKT(43,12,1:6)=(/ 0.820173,0.783880,0.755007,0.732877,0.745810,0.738147 /)
XEXT_COEFF_550_LKT(43,12)=578.200000 !rg=0.300455 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(43,13,1:6)=(/ 427.410000,445.380000,471.220000,501.070000,500.000000,379.740000 /)
XPIZA_LKT(43,13,1:6)=(/ 0.841424,0.882668,0.917605,0.977306,0.987932,0.992621 /)
XCGA_LKT(43,13,1:6)=(/ 0.830347,0.795563,0.765450,0.737220,0.745250,0.741507 /)
XEXT_COEFF_550_LKT(43,13)=472.670000 !rg=0.300455 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(43,14,1:6)=(/ 353.650000,367.840000,384.380000,411.050000,418.620000,345.810000 /)
XPIZA_LKT(43,14,1:6)=(/ 0.828608,0.866813,0.905507,0.972173,0.985760,0.991679 /)
XCGA_LKT(43,14,1:6)=(/ 0.837510,0.808170,0.776580,0.744090,0.747793,0.739610 /)
XEXT_COEFF_550_LKT(43,14)=384.220000 !rg=0.300455 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(43,15,1:6)=(/ 292.280000,302.160000,316.150000,333.760000,347.400000,308.310000 /)
XPIZA_LKT(43,15,1:6)=(/ 0.817943,0.852900,0.889500,0.967604,0.982728,0.990454 /)
XCGA_LKT(43,15,1:6)=(/ 0.846277,0.816910,0.789260,0.749540,0.748993,0.739160 /)
XEXT_COEFF_550_LKT(43,15)=313.440000 !rg=0.300455 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(43,16,1:6)=(/ 240.770000,247.150000,257.400000,269.780000,286.110000,269.320000 /)
XPIZA_LKT(43,16,1:6)=(/ 0.807548,0.839931,0.874183,0.962391,0.979330,0.988947 /)
XCGA_LKT(43,16,1:6)=(/ 0.853800,0.829677,0.799463,0.763103,0.750923,0.742220 /)
XEXT_COEFF_550_LKT(43,16)=257.840000 !rg=0.300455 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(43,17,1:6)=(/ 199.260000,204.310000,210.700000,221.970000,233.600000,228.830000 /)
XPIZA_LKT(43,17,1:6)=(/ 0.798536,0.828769,0.860238,0.954778,0.976089,0.987047 /)
XCGA_LKT(43,17,1:6)=(/ 0.859190,0.838053,0.812107,0.771230,0.760057,0.744533 /)
XEXT_COEFF_550_LKT(43,17)=210.940000 !rg=0.300455 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(43,18,1:6)=(/ 164.800000,168.610000,174.010000,181.400000,191.580000,193.510000 /)
XPIZA_LKT(43,18,1:6)=(/ 0.789680,0.817227,0.846414,0.946504,0.971198,0.984420 /)
XCGA_LKT(43,18,1:6)=(/ 0.864697,0.844847,0.822580,0.778303,0.765053,0.746660 /)
XEXT_COEFF_550_LKT(43,18)=172.540000 !rg=0.300455 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(43,19,1:6)=(/ 136.050000,138.480000,142.390000,147.280000,156.980000,160.990000 /)
XPIZA_LKT(43,19,1:6)=(/ 0.780310,0.807156,0.834950,0.939180,0.965943,0.981747 /)
XCGA_LKT(43,19,1:6)=(/ 0.869697,0.853153,0.830807,0.791473,0.770070,0.748030 /)
XEXT_COEFF_550_LKT(43,19)=142.730000 !rg=0.300455 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(43,20,1:6)=(/ 112.500000,114.600000,117.070000,121.430000,127.780000,133.140000 /)
XPIZA_LKT(43,20,1:6)=(/ 0.770754,0.798080,0.824010,0.931539,0.959974,0.978975 /)
XCGA_LKT(43,20,1:6)=(/ 0.873343,0.859133,0.840187,0.799423,0.781677,0.756567 /)
XEXT_COEFF_550_LKT(43,20)=117.000000 !rg=0.300455 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(44,1,1:6)=(/ 1776.500000,3337.000000,3400.500000,1644.800000,320.050000,27.952000 /)
XPIZA_LKT(44,1,1:6)=(/ 0.930939,0.977772,0.986187,0.993760,0.987371,0.943623 /)
XCGA_LKT(44,1,1:6)=(/ 0.581197,0.746027,0.777597,0.661043,0.351203,0.080340 /)
XEXT_COEFF_550_LKT(44,1)=3428.900000 !rg=0.3255 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(44,2,1:6)=(/ 1910.700000,2978.600000,3241.400000,1714.300000,358.560000,34.326000 /)
XPIZA_LKT(44,2,1:6)=(/ 0.937285,0.974874,0.985550,0.993781,0.988307,0.953277 /)
XCGA_LKT(44,2,1:6)=(/ 0.653630,0.730980,0.775920,0.694053,0.428713,0.101803 /)
XEXT_COEFF_550_LKT(44,2)=3263.700000 !rg=0.3255 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(44,3,1:6)=(/ 1908.100000,2513.900000,2958.600000,1794.200000,440.110000,48.203000 /)
XPIZA_LKT(44,3,1:6)=(/ 0.937548,0.970343,0.984041,0.993924,0.989903,0.965643 /)
XCGA_LKT(44,3,1:6)=(/ 0.697603,0.712327,0.769583,0.719500,0.523173,0.153153 /)
XEXT_COEFF_550_LKT(44,3)=2970.100000 !rg=0.3255 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(44,4,1:6)=(/ 1763.600000,2129.900000,2560.700000,1823.600000,546.920000,71.342000 /)
XPIZA_LKT(44,4,1:6)=(/ 0.933673,0.965084,0.981391,0.993889,0.991421,0.975693 /)
XCGA_LKT(44,4,1:6)=(/ 0.722103,0.705270,0.757180,0.736820,0.597793,0.246907 /)
XEXT_COEFF_550_LKT(44,4)=2562.200000 !rg=0.3255 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(44,5,1:6)=(/ 1566.000000,1817.600000,2140.100000,1764.000000,659.570000,105.540000 /)
XPIZA_LKT(44,5,1:6)=(/ 0.928540,0.959027,0.977664,0.993539,0.992523,0.982531 /)
XCGA_LKT(44,5,1:6)=(/ 0.734573,0.707750,0.743787,0.744997,0.655583,0.365080 /)
XEXT_COEFF_550_LKT(44,5)=2138.100000 !rg=0.3255 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(44,6,1:6)=(/ 1354.100000,1528.300000,1768.800000,1619.300000,757.510000,150.840000 /)
XPIZA_LKT(44,6,1:6)=(/ 0.923620,0.952854,0.972766,0.992817,0.993191,0.986887 /)
XCGA_LKT(44,6,1:6)=(/ 0.756023,0.716407,0.733073,0.745047,0.697610,0.472940 /)
XEXT_COEFF_550_LKT(44,6)=1762.100000 !rg=0.3255 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(44,7,1:6)=(/ 1165.300000,1281.400000,1451.200000,1424.700000,822.370000,205.870000 /)
XPIZA_LKT(44,7,1:6)=(/ 0.912761,0.945205,0.966864,0.991693,0.993472,0.989679 /)
XCGA_LKT(44,7,1:6)=(/ 0.767817,0.725733,0.727400,0.740707,0.725897,0.559737 /)
XEXT_COEFF_550_LKT(44,7)=1443.000000 !rg=0.3255 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(44,8,1:6)=(/ 978.280000,1066.200000,1187.200000,1219.400000,837.380000,265.490000 /)
XPIZA_LKT(44,8,1:6)=(/ 0.904681,0.937369,0.960367,0.990156,0.993365,0.991442 /)
XCGA_LKT(44,8,1:6)=(/ 0.781053,0.739447,0.728223,0.735357,0.741540,0.625697 /)
XEXT_COEFF_550_LKT(44,8)=1182.000000 !rg=0.3255 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(44,9,1:6)=(/ 827.220000,889.310000,967.260000,1022.100000,800.940000,323.310000 /)
XPIZA_LKT(44,9,1:6)=(/ 0.890631,0.927123,0.953335,0.987953,0.992868,0.992535 /)
XCGA_LKT(44,9,1:6)=(/ 0.788320,0.751727,0.730927,0.730500,0.747390,0.675110 /)
XEXT_COEFF_550_LKT(44,9)=965.440000 !rg=0.3255 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(44,10,1:6)=(/ 690.550000,731.040000,790.950000,843.440000,729.420000,369.500000 /)
XPIZA_LKT(44,10,1:6)=(/ 0.877204,0.917735,0.945064,0.985806,0.991986,0.993122 /)
XCGA_LKT(44,10,1:6)=(/ 0.804023,0.763317,0.739040,0.729277,0.747093,0.709643 /)
XEXT_COEFF_550_LKT(44,10)=794.090000 !rg=0.3255 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(44,11,1:6)=(/ 572.070000,603.510000,644.390000,692.150000,645.250000,395.080000 /)
XPIZA_LKT(44,11,1:6)=(/ 0.862714,0.905413,0.935991,0.982822,0.990779,0.993280 /)
XCGA_LKT(44,11,1:6)=(/ 0.814053,0.777693,0.749517,0.730210,0.746880,0.731117 /)
XEXT_COEFF_550_LKT(44,11)=646.290000 !rg=0.3255 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(44,12,1:6)=(/ 474.660000,495.980000,526.600000,563.710000,554.680000,393.620000 /)
XPIZA_LKT(44,12,1:6)=(/ 0.848420,0.891150,0.924460,0.979559,0.989153,0.993001 /)
XCGA_LKT(44,12,1:6)=(/ 0.825637,0.788887,0.759377,0.733833,0.746470,0.740560 /)
XEXT_COEFF_550_LKT(44,12)=528.830000 !rg=0.3255 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(44,13,1:6)=(/ 392.880000,409.270000,429.240000,461.110000,465.520000,368.210000 /)
XPIZA_LKT(44,13,1:6)=(/ 0.834511,0.875060,0.912681,0.974764,0.987134,0.992293 /)
XCGA_LKT(44,13,1:6)=(/ 0.833463,0.802370,0.770320,0.739693,0.746907,0.741417 /)
XEXT_COEFF_550_LKT(44,13)=429.400000 !rg=0.3255 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(44,14,1:6)=(/ 325.040000,336.780000,353.020000,374.880000,388.300000,331.470000 /)
XPIZA_LKT(44,14,1:6)=(/ 0.823337,0.860002,0.898244,0.970395,0.984622,0.991206 /)
XCGA_LKT(44,14,1:6)=(/ 0.843107,0.811177,0.782920,0.744063,0.748347,0.739033 /)
XEXT_COEFF_550_LKT(44,14)=349.600000 !rg=0.3255 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(44,15,1:6)=(/ 267.920000,276.750000,288.570000,304.210000,320.330000,293.190000 /)
XPIZA_LKT(44,15,1:6)=(/ 0.812106,0.846465,0.882428,0.965216,0.981413,0.989935 /)
XCGA_LKT(44,15,1:6)=(/ 0.849727,0.824383,0.792353,0.757497,0.748263,0.741527 /)
XEXT_COEFF_550_LKT(44,15)=288.220000 !rg=0.3255 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(44,16,1:6)=(/ 221.280000,227.610000,234.950000,248.210000,263.290000,253.390000 /)
XPIZA_LKT(44,16,1:6)=(/ 0.803354,0.834180,0.867764,0.958925,0.977500,0.988027 /)
XCGA_LKT(44,16,1:6)=(/ 0.856257,0.834173,0.807430,0.765987,0.754680,0.743223 /)
XEXT_COEFF_550_LKT(44,16)=235.960000 !rg=0.3255 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(44,17,1:6)=(/ 183.300000,187.650000,193.810000,203.260000,214.110000,213.370000 /)
XPIZA_LKT(44,17,1:6)=(/ 0.793715,0.822585,0.853987,0.950786,0.973890,0.986127 /)
XCGA_LKT(44,17,1:6)=(/ 0.862800,0.840757,0.817607,0.772000,0.762310,0.745043 /)
XEXT_COEFF_550_LKT(44,17)=192.850000 !rg=0.3255 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(44,18,1:6)=(/ 151.260000,154.830000,159.390000,165.720000,176.320000,179.730000 /)
XPIZA_LKT(44,18,1:6)=(/ 0.784917,0.812512,0.840893,0.942581,0.968536,0.983381 /)
XCGA_LKT(44,18,1:6)=(/ 0.866753,0.849683,0.825730,0.785773,0.764800,0.746327 /)
XEXT_COEFF_550_LKT(44,18)=159.070000 !rg=0.3255 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(44,19,1:6)=(/ 125.200000,127.740000,130.490000,135.840000,143.960000,149.070000 /)
XPIZA_LKT(44,19,1:6)=(/ 0.775786,0.802763,0.829407,0.935353,0.963104,0.980332 /)
XCGA_LKT(44,19,1:6)=(/ 0.871350,0.856420,0.836883,0.794973,0.775160,0.752330 /)
XEXT_COEFF_550_LKT(44,19)=131.030000 !rg=0.3255 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(44,20,1:6)=(/ 103.610000,105.310000,107.750000,111.540000,116.890000,122.450000 /)
XPIZA_LKT(44,20,1:6)=(/ 0.765966,0.792976,0.818913,0.928326,0.956673,0.977228 /)
XCGA_LKT(44,20,1:6)=(/ 0.875813,0.860943,0.844453,0.800817,0.785310,0.758157 /)
XEXT_COEFF_550_LKT(44,20)=107.500000 !rg=0.3255 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(45,1,1:6)=(/ 1851.500000,2776.600000,3286.000000,1738.500000,368.870000,34.873000 /)
XPIZA_LKT(45,1,1:6)=(/ 0.934500,0.972892,0.985839,0.993917,0.988617,0.953778 /)
XCGA_LKT(45,1,1:6)=(/ 0.660900,0.719943,0.784697,0.695133,0.420183,0.094297 /)
XEXT_COEFF_550_LKT(45,1)=3299.800000 !rg=0.352632 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(45,2,1:6)=(/ 1884.100000,2467.700000,3121.200000,1824.400000,415.450000,42.656000 /)
XPIZA_LKT(45,2,1:6)=(/ 0.937434,0.969635,0.984849,0.994030,0.989437,0.961499 /)
XCGA_LKT(45,2,1:6)=(/ 0.702643,0.701507,0.778347,0.717463,0.493017,0.119523 /)
XEXT_COEFF_550_LKT(45,2)=3134.900000 !rg=0.352632 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(45,3,1:6)=(/ 1780.600000,2140.200000,2779.600000,1873.000000,504.500000,59.057000 /)
XPIZA_LKT(45,3,1:6)=(/ 0.934657,0.965005,0.982879,0.994075,0.990832,0.971240 /)
XCGA_LKT(45,3,1:6)=(/ 0.721913,0.691907,0.767627,0.735340,0.567030,0.179367 /)
XEXT_COEFF_550_LKT(45,3)=2783.700000 !rg=0.352632 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(45,4,1:6)=(/ 1613.200000,1862.800000,2353.800000,1853.900000,611.260000,85.561000 /)
XPIZA_LKT(45,4,1:6)=(/ 0.930213,0.960228,0.979678,0.993895,0.992071,0.979152 /)
XCGA_LKT(45,4,1:6)=(/ 0.737197,0.696213,0.752100,0.746860,0.627940,0.283440 /)
XEXT_COEFF_550_LKT(45,4)=2353.400000 !rg=0.352632 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(45,5,1:6)=(/ 1437.800000,1607.900000,1949.400000,1742.300000,717.150000,123.940000 /)
XPIZA_LKT(45,5,1:6)=(/ 0.922665,0.954273,0.975246,0.993375,0.992937,0.984649 /)
XCGA_LKT(45,5,1:6)=(/ 0.751103,0.705217,0.737273,0.749993,0.677043,0.403827 /)
XEXT_COEFF_550_LKT(45,5)=1941.800000 !rg=0.352632 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(45,6,1:6)=(/ 1240.000000,1370.200000,1598.900000,1560.000000,800.780000,173.400000 /)
XPIZA_LKT(45,6,1:6)=(/ 0.917961,0.948729,0.970076,0.992466,0.993415,0.988218 /)
XCGA_LKT(45,6,1:6)=(/ 0.764693,0.719180,0.728343,0.746593,0.712427,0.507127 /)
XEXT_COEFF_550_LKT(45,6)=1599.300000 !rg=0.352632 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(45,7,1:6)=(/ 1060.800000,1153.300000,1310.400000,1347.500000,844.290000,231.110000 /)
XPIZA_LKT(45,7,1:6)=(/ 0.907741,0.941552,0.963742,0.991145,0.993526,0.990518 /)
XCGA_LKT(45,7,1:6)=(/ 0.775970,0.733213,0.725050,0.740427,0.735130,0.586883 /)
XEXT_COEFF_550_LKT(45,7)=1310.600000 !rg=0.352632 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(45,8,1:6)=(/ 900.050000,965.420000,1073.000000,1136.200000,835.360000,290.820000 /)
XPIZA_LKT(45,8,1:6)=(/ 0.897537,0.932800,0.957005,0.989404,0.993255,0.991971 /)
XCGA_LKT(45,8,1:6)=(/ 0.788170,0.744160,0.727503,0.734353,0.746230,0.646717 /)
XEXT_COEFF_550_LKT(45,8)=1075.500000 !rg=0.352632 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(45,9,1:6)=(/ 758.110000,806.510000,882.700000,941.390000,780.580000,345.220000 /)
XPIZA_LKT(45,9,1:6)=(/ 0.882910,0.921829,0.948468,0.987027,0.992582,0.992844 /)
XCGA_LKT(45,9,1:6)=(/ 0.798550,0.754387,0.733353,0.729263,0.748780,0.690297 /)
XEXT_COEFF_550_LKT(45,9)=874.930000 !rg=0.352632 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(45,10,1:6)=(/ 633.120000,670.690000,718.860000,778.600000,700.100000,383.860000 /)
XPIZA_LKT(45,10,1:6)=(/ 0.869086,0.910726,0.940735,0.984204,0.991563,0.993246 /)
XCGA_LKT(45,10,1:6)=(/ 0.807377,0.770760,0.741627,0.728850,0.748113,0.719523 /)
XEXT_COEFF_550_LKT(45,10)=718.250000 !rg=0.352632 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(45,11,1:6)=(/ 526.220000,552.000000,588.490000,633.320000,610.930000,398.870000 /)
XPIZA_LKT(45,11,1:6)=(/ 0.855070,0.897853,0.930609,0.981034,0.990192,0.993226 /)
XCGA_LKT(45,11,1:6)=(/ 0.821040,0.781240,0.753023,0.729697,0.747453,0.736187 /)
XEXT_COEFF_550_LKT(45,11)=589.060000 !rg=0.352632 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(45,12,1:6)=(/ 436.020000,455.760000,479.900000,517.950000,518.300000,386.680000 /)
XPIZA_LKT(45,12,1:6)=(/ 0.841014,0.882954,0.919377,0.977246,0.988443,0.992775 /)
XCGA_LKT(45,12,1:6)=(/ 0.828910,0.796167,0.763537,0.735673,0.747217,0.741920 /)
XEXT_COEFF_550_LKT(45,12)=479.750000 !rg=0.352632 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(45,13,1:6)=(/ 360.930000,374.650000,394.060000,420.800000,432.640000,354.860000 /)
XPIZA_LKT(45,13,1:6)=(/ 0.828818,0.867477,0.905935,0.972844,0.986147,0.991901 /)
XCGA_LKT(45,13,1:6)=(/ 0.839570,0.805407,0.776480,0.739483,0.747120,0.740623 /)
XEXT_COEFF_550_LKT(45,13)=390.330000 !rg=0.352632 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(45,14,1:6)=(/ 297.570000,308.440000,322.420000,342.470000,358.930000,316.810000 /)
XPIZA_LKT(45,14,1:6)=(/ 0.817905,0.853443,0.890851,0.967677,0.982933,0.990733 /)
XCGA_LKT(45,14,1:6)=(/ 0.846273,0.818763,0.785627,0.751073,0.747077,0.740347 /)
XEXT_COEFF_550_LKT(45,14)=321.440000 !rg=0.352632 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(45,15,1:6)=(/ 245.770000,253.880000,263.640000,277.970000,296.450000,278.240000 /)
XPIZA_LKT(45,15,1:6)=(/ 0.807533,0.840504,0.875467,0.962616,0.979200,0.989229 /)
XCGA_LKT(45,15,1:6)=(/ 0.854210,0.829613,0.800983,0.761427,0.750897,0.742640 /)
XEXT_COEFF_550_LKT(45,15)=263.850000 !rg=0.352632 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(45,16,1:6)=(/ 203.420000,209.190000,216.070000,227.620000,240.250000,235.880000 /)
XPIZA_LKT(45,16,1:6)=(/ 0.798046,0.827996,0.861121,0.954847,0.975881,0.987469 /)
XCGA_LKT(45,16,1:6)=(/ 0.860120,0.837120,0.812483,0.767083,0.758127,0.744000 /)
XEXT_COEFF_550_LKT(45,16)=216.070000 !rg=0.352632 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(45,17,1:6)=(/ 168.210000,172.470000,177.710000,186.380000,197.040000,199.300000 /)
XPIZA_LKT(45,17,1:6)=(/ 0.789444,0.817547,0.847343,0.946519,0.971457,0.985000 /)
XCGA_LKT(45,17,1:6)=(/ 0.864817,0.845710,0.820490,0.778763,0.761873,0.745333 /)
XEXT_COEFF_550_LKT(45,17)=177.310000 !rg=0.352632 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(45,18,1:6)=(/ 139.070000,142.190000,146.100000,151.830000,161.990000,166.880000 /)
XPIZA_LKT(45,18,1:6)=(/ 0.780514,0.807531,0.835380,0.939677,0.965796,0.981342 /)
XCGA_LKT(45,18,1:6)=(/ 0.869873,0.853297,0.832067,0.790767,0.769203,0.748067 /)
XEXT_COEFF_550_LKT(45,18)=145.940000 !rg=0.352632 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(45,19,1:6)=(/ 115.220000,117.380000,120.250000,124.750000,131.050000,137.170000 /)
XPIZA_LKT(45,19,1:6)=(/ 0.771140,0.797561,0.824278,0.931668,0.960085,0.978851 /)
XCGA_LKT(45,19,1:6)=(/ 0.873780,0.858610,0.840793,0.796423,0.780000,0.755130 /)
XEXT_COEFF_550_LKT(45,19)=120.130000 !rg=0.352632 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(45,20,1:6)=(/ 95.276000,97.053000,99.023000,102.440000,107.610000,113.370000 /)
XPIZA_LKT(45,20,1:6)=(/ 0.761329,0.788766,0.813565,0.924857,0.952594,0.974675 /)
XCGA_LKT(45,20,1:6)=(/ 0.877270,0.864313,0.846840,0.806420,0.785487,0.757603 /)
XEXT_COEFF_550_LKT(45,20)=99.013000 !rg=0.352632 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(46,1,1:6)=(/ 1999.000000,2192.300000,3156.500000,1887.700000,423.140000,43.460000 /)
XPIZA_LKT(46,1,1:6)=(/ 0.938500,0.965919,0.984755,0.994067,0.989579,0.961988 /)
XCGA_LKT(46,1,1:6)=(/ 0.731800,0.674110,0.783413,0.734333,0.496907,0.110773 /)
XEXT_COEFF_550_LKT(46,1)=3176.700000 !rg=0.382026 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(46,2,1:6)=(/ 1855.100000,2011.900000,2935.200000,1919.500000,480.980000,52.819000 /)
XPIZA_LKT(46,2,1:6)=(/ 0.936632,0.962726,0.983755,0.994208,0.990416,0.968108 /)
XCGA_LKT(46,2,1:6)=(/ 0.742310,0.667930,0.776687,0.737990,0.551957,0.140480 /)
XEXT_COEFF_550_LKT(46,2)=2939.200000 !rg=0.382026 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(46,3,1:6)=(/ 1659.800000,1818.500000,2554.400000,1927.600000,573.110000,71.843000 /)
XPIZA_LKT(46,3,1:6)=(/ 0.932209,0.958725,0.981236,0.994149,0.991623,0.975718 /)
XCGA_LKT(46,3,1:6)=(/ 0.744037,0.670703,0.761930,0.748650,0.604447,0.210017 /)
XEXT_COEFF_550_LKT(46,3)=2550.400000 !rg=0.382026 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(46,4,1:6)=(/ 1479.900000,1635.600000,2130.700000,1856.800000,675.980000,101.870000 /)
XPIZA_LKT(46,4,1:6)=(/ 0.926872,0.954691,0.977417,0.993817,0.992615,0.981953 /)
XCGA_LKT(46,4,1:6)=(/ 0.751903,0.688190,0.744400,0.754527,0.654490,0.323173 /)
XEXT_COEFF_550_LKT(46,4)=2122.000000 !rg=0.382026 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(46,5,1:6)=(/ 1309.300000,1434.900000,1750.700000,1697.200000,771.300000,144.530000 /)
XPIZA_LKT(46,5,1:6)=(/ 0.918802,0.949930,0.972525,0.993116,0.993272,0.986396 /)
XCGA_LKT(46,5,1:6)=(/ 0.761887,0.708317,0.730657,0.753013,0.696000,0.442683 /)
XEXT_COEFF_550_LKT(46,5)=1748.500000 !rg=0.382026 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(46,6,1:6)=(/ 1134.600000,1240.100000,1440.300000,1486.000000,836.940000,197.590000 /)
XPIZA_LKT(46,6,1:6)=(/ 0.911714,0.943081,0.967014,0.991992,0.993572,0.989322 /)
XCGA_LKT(46,6,1:6)=(/ 0.770210,0.721810,0.724717,0.746080,0.725193,0.539117 /)
XEXT_COEFF_550_LKT(46,6)=1435.400000 !rg=0.382026 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(46,7,1:6)=(/ 968.220000,1046.800000,1182.000000,1262.100000,857.020000,257.220000 /)
XPIZA_LKT(46,7,1:6)=(/ 0.902619,0.936058,0.960673,0.990436,0.993514,0.991222 /)
XCGA_LKT(46,7,1:6)=(/ 0.783093,0.737350,0.725030,0.738237,0.742653,0.612107 /)
XEXT_COEFF_550_LKT(46,7)=1179.400000 !rg=0.382026 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(46,8,1:6)=(/ 825.230000,883.130000,970.680000,1053.100000,825.240000,315.690000 /)
XPIZA_LKT(46,8,1:6)=(/ 0.889608,0.926532,0.953345,0.988402,0.993079,0.992402 /)
XCGA_LKT(46,8,1:6)=(/ 0.791693,0.751000,0.728517,0.731777,0.749780,0.665623 /)
XEXT_COEFF_550_LKT(46,8)=969.180000 !rg=0.382026 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(46,9,1:6)=(/ 692.160000,731.650000,798.240000,865.460000,754.440000,364.930000 /)
XPIZA_LKT(46,9,1:6)=(/ 0.875553,0.916918,0.944157,0.986161,0.992261,0.993071 /)
XCGA_LKT(46,9,1:6)=(/ 0.804800,0.765607,0.734143,0.730423,0.750237,0.703670 /)
XEXT_COEFF_550_LKT(46,9)=797.420000 !rg=0.382026 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(46,10,1:6)=(/ 581.010000,610.760000,657.470000,710.930000,668.320000,394.460000 /)
XPIZA_LKT(46,10,1:6)=(/ 0.861663,0.904539,0.934959,0.982937,0.991061,0.993301 /)
XCGA_LKT(46,10,1:6)=(/ 0.815933,0.773423,0.746303,0.727577,0.748453,0.727770 /)
XEXT_COEFF_550_LKT(46,10)=650.750000 !rg=0.382026 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(46,11,1:6)=(/ 482.700000,506.080000,536.090000,581.900000,574.950000,398.060000 /)
XPIZA_LKT(46,11,1:6)=(/ 0.847797,0.890785,0.924621,0.979108,0.989563,0.993108 /)
XCGA_LKT(46,11,1:6)=(/ 0.824243,0.789820,0.755743,0.732003,0.747960,0.739813 /)
XEXT_COEFF_550_LKT(46,11)=534.170000 !rg=0.382026 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(46,12,1:6)=(/ 400.490000,416.660000,440.290000,472.100000,482.540000,376.590000 /)
XPIZA_LKT(46,12,1:6)=(/ 0.834664,0.875485,0.912675,0.975513,0.987425,0.992476 /)
XCGA_LKT(46,12,1:6)=(/ 0.835510,0.799260,0.769777,0.735210,0.746293,0.741987 /)
XEXT_COEFF_550_LKT(46,12)=435.650000 !rg=0.382026 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(46,13,1:6)=(/ 330.360000,343.320000,359.720000,384.390000,401.570000,340.470000 /)
XPIZA_LKT(46,13,1:6)=(/ 0.823177,0.860499,0.898962,0.970214,0.984819,0.991497 /)
XCGA_LKT(46,13,1:6)=(/ 0.842863,0.813420,0.779000,0.745553,0.746613,0.741037 /)
XEXT_COEFF_550_LKT(46,13)=358.280000 !rg=0.382026 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(46,14,1:6)=(/ 273.130000,282.600000,294.900000,311.510000,331.470000,302.240000 /)
XPIZA_LKT(46,14,1:6)=(/ 0.812215,0.846453,0.883281,0.965606,0.981397,0.990230 /)
XCGA_LKT(46,14,1:6)=(/ 0.851650,0.824137,0.794257,0.755547,0.748703,0.741847 /)
XEXT_COEFF_550_LKT(46,14)=294.600000 !rg=0.382026 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(46,15,1:6)=(/ 226.500000,233.210000,241.440000,254.520000,270.410000,261.150000 /)
XPIZA_LKT(46,15,1:6)=(/ 0.803129,0.834621,0.868877,0.959530,0.978422,0.988613 /)
XCGA_LKT(46,15,1:6)=(/ 0.857707,0.833177,0.806673,0.763163,0.754217,0.743760 /)
XEXT_COEFF_550_LKT(46,15)=242.130000 !rg=0.382026 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(46,16,1:6)=(/ 187.240000,192.220000,198.350000,208.500000,220.560000,220.280000 /)
XPIZA_LKT(46,16,1:6)=(/ 0.793488,0.822868,0.854061,0.951886,0.974684,0.986520 /)
XCGA_LKT(46,16,1:6)=(/ 0.862483,0.842487,0.815637,0.773797,0.760003,0.744477 /)
XEXT_COEFF_550_LKT(46,16)=197.700000 !rg=0.382026 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(46,17,1:6)=(/ 154.710000,158.250000,163.070000,169.730000,181.390000,185.820000 /)
XPIZA_LKT(46,17,1:6)=(/ 0.784357,0.811871,0.841406,0.942920,0.969010,0.983540 /)
XCGA_LKT(46,17,1:6)=(/ 0.868100,0.849143,0.826733,0.783270,0.765757,0.747253 /)
XEXT_COEFF_550_LKT(46,17)=163.110000 !rg=0.382026 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(46,18,1:6)=(/ 128.180000,130.950000,134.110000,139.410000,147.100000,153.240000 /)
XPIZA_LKT(46,18,1:6)=(/ 0.776444,0.803116,0.829825,0.936003,0.964224,0.980834 /)
XCGA_LKT(46,18,1:6)=(/ 0.872200,0.855853,0.836660,0.792947,0.774300,0.751243 /)
XEXT_COEFF_550_LKT(46,18)=134.200000 !rg=0.382026 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(46,19,1:6)=(/ 106.030000,108.060000,110.470000,114.440000,119.950000,126.340000 /)
XPIZA_LKT(46,19,1:6)=(/ 0.766491,0.793581,0.818487,0.929013,0.957839,0.977706 /)
XCGA_LKT(46,19,1:6)=(/ 0.875453,0.862057,0.843490,0.802273,0.782683,0.756013 /)
XEXT_COEFF_550_LKT(46,19)=110.170000 !rg=0.382026 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(46,20,1:6)=(/ 87.758000,89.230000,91.069000,93.586000,98.702000,104.620000 /)
XPIZA_LKT(46,20,1:6)=(/ 0.755969,0.783685,0.808810,0.921959,0.949164,0.972781 /)
XCGA_LKT(46,20,1:6)=(/ 0.879540,0.866760,0.851133,0.810007,0.790503,0.761263 /)
XEXT_COEFF_550_LKT(46,20)=91.099000 !rg=0.382026 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(47,1,1:6)=(/ 2010.700000,1727.500000,2940.100000,2012.400000,489.250000,54.010000 /)
XPIZA_LKT(47,1,1:6)=(/ 0.940632,0.955996,0.983884,0.994355,0.990421,0.968599 /)
XCGA_LKT(47,1,1:6)=(/ 0.784223,0.620833,0.779807,0.754760,0.570850,0.130303 /)
XEXT_COEFF_550_LKT(47,1)=2932.200000 !rg=0.41387 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(47,2,1:6)=(/ 1763.500000,1647.700000,2687.800000,1990.200000,555.280000,65.037000 /)
XPIZA_LKT(47,2,1:6)=(/ 0.934318,0.954009,0.982135,0.994320,0.991298,0.973400 /)
XCGA_LKT(47,2,1:6)=(/ 0.766057,0.628753,0.770303,0.753537,0.598913,0.165360 /)
XEXT_COEFF_550_LKT(47,2)=2683.200000 !rg=0.41387 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(47,3,1:6)=(/ 1548.700000,1567.500000,2294.600000,1953.700000,643.970000,86.682000 /)
XPIZA_LKT(47,3,1:6)=(/ 0.926990,0.952435,0.978980,0.994141,0.992294,0.979298 /)
XCGA_LKT(47,3,1:6)=(/ 0.763177,0.661773,0.751890,0.759127,0.635710,0.245640 /)
XEXT_COEFF_550_LKT(47,3)=2285.900000 !rg=0.41387 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(47,4,1:6)=(/ 1372.900000,1450.100000,1898.700000,1830.600000,739.250000,120.390000 /)
XPIZA_LKT(47,4,1:6)=(/ 0.921018,0.949371,0.974522,0.993651,0.993064,0.984233 /)
XCGA_LKT(47,4,1:6)=(/ 0.766847,0.689223,0.733467,0.759780,0.677837,0.365433 /)
XEXT_COEFF_550_LKT(47,4)=1891.000000 !rg=0.41387 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(47,5,1:6)=(/ 1194.200000,1288.400000,1565.000000,1629.500000,820.280000,167.230000 /)
XPIZA_LKT(47,5,1:6)=(/ 0.915187,0.944712,0.969310,0.992749,0.993532,0.987840 /)
XCGA_LKT(47,5,1:6)=(/ 0.772897,0.710837,0.723747,0.753797,0.712663,0.480570 /)
XEXT_COEFF_550_LKT(47,5)=1555.800000 !rg=0.41387 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(47,6,1:6)=(/ 1046.600000,1117.900000,1295.900000,1397.100000,864.740000,223.290000 /)
XPIZA_LKT(47,6,1:6)=(/ 0.904871,0.937567,0.963042,0.991446,0.993663,0.990241 /)
XCGA_LKT(47,6,1:6)=(/ 0.782370,0.725660,0.719490,0.745163,0.736050,0.569123 /)
XEXT_COEFF_550_LKT(47,6)=1287.600000 !rg=0.41387 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(47,7,1:6)=(/ 892.950000,948.980000,1066.300000,1169.800000,859.930000,283.470000 /)
XPIZA_LKT(47,7,1:6)=(/ 0.895355,0.930267,0.956438,0.989675,0.993441,0.991804 /)
XCGA_LKT(47,7,1:6)=(/ 0.792710,0.741513,0.722960,0.736547,0.748567,0.635023 /)
XEXT_COEFF_550_LKT(47,7)=1064.800000 !rg=0.41387 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(47,8,1:6)=(/ 757.370000,802.340000,882.550000,966.720000,807.400000,339.100000 /)
XPIZA_LKT(47,8,1:6)=(/ 0.882724,0.920438,0.948425,0.987334,0.992833,0.992752 /)
XCGA_LKT(47,8,1:6)=(/ 0.802067,0.753637,0.730093,0.729807,0.752003,0.682327 /)
XEXT_COEFF_550_LKT(47,8)=874.480000 !rg=0.41387 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(47,9,1:6)=(/ 632.140000,667.360000,722.430000,793.110000,725.290000,381.400000 /)
XPIZA_LKT(47,9,1:6)=(/ 0.869513,0.911140,0.940429,0.984900,0.991860,0.993232 /)
XCGA_LKT(47,9,1:6)=(/ 0.812180,0.771923,0.741507,0.728633,0.751137,0.714873 /)
XEXT_COEFF_550_LKT(47,9)=722.800000 !rg=0.41387 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(47,10,1:6)=(/ 530.860000,556.940000,596.610000,648.420000,632.910000,400.670000 /)
XPIZA_LKT(47,10,1:6)=(/ 0.854496,0.897863,0.929802,0.981683,0.990544,0.993286 /)
XCGA_LKT(47,10,1:6)=(/ 0.820850,0.784230,0.747883,0.731107,0.749643,0.734007 /)
XEXT_COEFF_550_LKT(47,10)=595.380000 !rg=0.41387 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(47,11,1:6)=(/ 442.550000,461.210000,490.720000,528.490000,537.860000,393.330000 /)
XPIZA_LKT(47,11,1:6)=(/ 0.840955,0.883756,0.918662,0.977954,0.988713,0.992917 /)
XCGA_LKT(47,11,1:6)=(/ 0.831423,0.793687,0.762993,0.732903,0.747317,0.741857 /)
XEXT_COEFF_550_LKT(47,11)=486.110000 !rg=0.41387 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(47,12,1:6)=(/ 366.480000,381.370000,401.210000,430.610000,448.160000,364.060000 /)
XPIZA_LKT(47,12,1:6)=(/ 0.828701,0.868015,0.906408,0.973303,0.986471,0.992144 /)
XCGA_LKT(47,12,1:6)=(/ 0.839210,0.808217,0.772120,0.741827,0.745893,0.742370 /)
XEXT_COEFF_550_LKT(47,12)=400.040000 !rg=0.41387 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(47,13,1:6)=(/ 303.130000,314.200000,329.050000,349.230000,371.620000,326.280000 /)
XPIZA_LKT(47,13,1:6)=(/ 0.816867,0.853292,0.891548,0.968270,0.983314,0.991023 /)
XCGA_LKT(47,13,1:6)=(/ 0.848700,0.818787,0.787680,0.749860,0.747617,0.741513 /)
XEXT_COEFF_550_LKT(47,13)=328.280000 !rg=0.41387 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(47,14,1:6)=(/ 251.890000,258.930000,269.320000,284.470000,303.330000,286.300000 /)
XPIZA_LKT(47,14,1:6)=(/ 0.807348,0.841265,0.876136,0.963395,0.980850,0.989642 /)
XCGA_LKT(47,14,1:6)=(/ 0.855247,0.828713,0.799900,0.758347,0.751393,0.743423 /)
XEXT_COEFF_550_LKT(47,14)=269.950000 !rg=0.41387 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(47,15,1:6)=(/ 208.280000,214.390000,221.310000,233.870000,247.850000,244.320000 /)
XPIZA_LKT(47,15,1:6)=(/ 0.798131,0.828133,0.862282,0.955766,0.976992,0.987899 /)
XCGA_LKT(47,15,1:6)=(/ 0.860007,0.838420,0.810637,0.768133,0.756877,0.744757 /)
XEXT_COEFF_550_LKT(47,15)=220.790000 !rg=0.41387 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(47,16,1:6)=(/ 172.120000,176.060000,182.090000,189.720000,202.910000,205.890000 /)
XPIZA_LKT(47,16,1:6)=(/ 0.789267,0.817887,0.848190,0.948565,0.972212,0.985256 /)
XCGA_LKT(47,16,1:6)=(/ 0.865687,0.845433,0.822257,0.777270,0.762387,0.745820 /)
XEXT_COEFF_550_LKT(47,16)=181.330000 !rg=0.41387 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(47,17,1:6)=(/ 142.580000,145.130000,149.400000,155.030000,165.200000,171.150000 /)
XPIZA_LKT(47,17,1:6)=(/ 0.780788,0.807551,0.835140,0.940517,0.967770,0.982970 /)
XCGA_LKT(47,17,1:6)=(/ 0.870303,0.853320,0.830973,0.788460,0.769060,0.748527 /)
XEXT_COEFF_550_LKT(47,17)=149.680000 !rg=0.41387 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(47,18,1:6)=(/ 117.890000,120.560000,123.230000,128.070000,134.500000,141.680000 /)
XPIZA_LKT(47,18,1:6)=(/ 0.771487,0.798017,0.824873,0.932340,0.961974,0.979603 /)
XCGA_LKT(47,18,1:6)=(/ 0.873700,0.859467,0.839760,0.797483,0.779180,0.753827 /)
XEXT_COEFF_550_LKT(47,18)=122.910000 !rg=0.41387 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(47,19,1:6)=(/ 97.618000,99.179000,101.550000,104.540000,110.460000,116.900000 /)
XPIZA_LKT(47,19,1:6)=(/ 0.762036,0.789221,0.814052,0.926163,0.954204,0.975424 /)
XCGA_LKT(47,19,1:6)=(/ 0.877557,0.864117,0.848007,0.805700,0.786883,0.758423 /)
XEXT_COEFF_550_LKT(47,19)=101.400000 !rg=0.41387 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(47,20,1:6)=(/ 80.878000,81.917000,83.654000,85.814000,90.137000,95.603000 /)
XPIZA_LKT(47,20,1:6)=(/ 0.751746,0.779518,0.803566,0.919906,0.946309,0.971633 /)
XCGA_LKT(47,20,1:6)=(/ 0.881060,0.869403,0.854290,0.814203,0.794590,0.764193 /)
XEXT_COEFF_550_LKT(47,20)=83.687000 !rg=0.41387 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(48,1,1:6)=(/ 1840.200000,1338.400000,2649.700000,2039.800000,573.870000,66.800000 /)
XPIZA_LKT(48,1,1:6)=(/ 0.938109,0.944449,0.981700,0.994502,0.991297,0.973903 /)
XCGA_LKT(48,1,1:6)=(/ 0.800133,0.568417,0.772177,0.759477,0.625370,0.153577 /)
XEXT_COEFF_550_LKT(48,1)=2647.500000 !rg=0.448369 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(48,2,1:6)=(/ 1609.500000,1378.600000,2391.300000,2030.200000,634.820000,79.459000 /)
XPIZA_LKT(48,2,1:6)=(/ 0.926989,0.945919,0.979800,0.994354,0.992094,0.977617 /)
XCGA_LKT(48,2,1:6)=(/ 0.782857,0.612347,0.758457,0.764830,0.631980,0.195017 /)
XEXT_COEFF_550_LKT(48,2)=2377.800000 !rg=0.448369 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(48,3,1:6)=(/ 1418.600000,1398.600000,2015.000000,1947.600000,714.720000,103.670000 /)
XPIZA_LKT(48,3,1:6)=(/ 0.920864,0.946185,0.976021,0.994046,0.992857,0.982162 /)
XCGA_LKT(48,3,1:6)=(/ 0.767833,0.660067,0.738337,0.766633,0.662113,0.286543 /)
XEXT_COEFF_550_LKT(48,3)=2001.500000 !rg=0.448369 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(48,4,1:6)=(/ 1260.400000,1312.000000,1670.300000,1775.500000,799.120000,141.250000 /)
XPIZA_LKT(48,4,1:6)=(/ 0.914530,0.943470,0.971087,0.993374,0.993430,0.986098 /)
XCGA_LKT(48,4,1:6)=(/ 0.769523,0.694587,0.722577,0.762370,0.698400,0.409143 /)
XEXT_COEFF_550_LKT(48,4)=1659.700000 !rg=0.448369 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(48,5,1:6)=(/ 1104.300000,1163.200000,1389.900000,1541.900000,862.290000,191.920000 /)
XPIZA_LKT(48,5,1:6)=(/ 0.908762,0.939664,0.965408,0.992277,0.993723,0.989038 /)
XCGA_LKT(48,5,1:6)=(/ 0.784080,0.717203,0.715350,0.753087,0.727130,0.516770 /)
XEXT_COEFF_550_LKT(48,5)=1385.500000 !rg=0.448369 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(48,6,1:6)=(/ 956.470000,1019.000000,1157.300000,1301.800000,882.910000,250.130000 /)
XPIZA_LKT(48,6,1:6)=(/ 0.899391,0.931398,0.959002,0.990720,0.993691,0.991012 /)
XCGA_LKT(48,6,1:6)=(/ 0.785667,0.736020,0.717160,0.742290,0.745083,0.597007 /)
XEXT_COEFF_550_LKT(48,6)=1152.600000 !rg=0.448369 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(48,7,1:6)=(/ 817.990000,868.750000,959.560000,1079.000000,853.460000,309.560000 /)
XPIZA_LKT(48,7,1:6)=(/ 0.888225,0.923908,0.952074,0.988591,0.993304,0.992283 /)
XCGA_LKT(48,7,1:6)=(/ 0.795043,0.751200,0.723273,0.733083,0.753213,0.655780 /)
XEXT_COEFF_550_LKT(48,7)=955.800000 !rg=0.448369 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(48,8,1:6)=(/ 691.490000,731.680000,796.590000,886.410000,782.730000,360.620000 /)
XPIZA_LKT(48,8,1:6)=(/ 0.875700,0.914281,0.943440,0.986131,0.992540,0.993018 /)
XCGA_LKT(48,8,1:6)=(/ 0.806670,0.765273,0.730893,0.729307,0.753927,0.697150 /)
XEXT_COEFF_550_LKT(48,8)=793.700000 !rg=0.448369 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(48,9,1:6)=(/ 582.410000,609.480000,655.730000,723.300000,691.880000,394.350000 /)
XPIZA_LKT(48,9,1:6)=(/ 0.861951,0.904613,0.935552,0.983373,0.991386,0.993322 /)
XCGA_LKT(48,9,1:6)=(/ 0.819230,0.776467,0.744923,0.727527,0.751527,0.724370 /)
XEXT_COEFF_550_LKT(48,9)=658.270000 !rg=0.448369 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(48,10,1:6)=(/ 485.460000,508.280000,542.150000,590.550000,596.710000,402.490000 /)
XPIZA_LKT(48,10,1:6)=(/ 0.848063,0.891563,0.924613,0.980208,0.989785,0.993206 /)
XCGA_LKT(48,10,1:6)=(/ 0.828393,0.791017,0.757140,0.731480,0.749277,0.738760 /)
XEXT_COEFF_550_LKT(48,10)=541.510000 !rg=0.448369 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(48,11,1:6)=(/ 406.300000,420.490000,445.750000,479.850000,499.930000,384.790000 /)
XPIZA_LKT(48,11,1:6)=(/ 0.833683,0.876687,0.913176,0.976615,0.987796,0.992678 /)
XCGA_LKT(48,11,1:6)=(/ 0.836597,0.803583,0.766043,0.738487,0.746273,0.743243 /)
XEXT_COEFF_550_LKT(48,11)=445.950000 !rg=0.448369 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(48,12,1:6)=(/ 336.120000,348.940000,366.290000,390.910000,416.490000,350.750000 /)
XPIZA_LKT(48,12,1:6)=(/ 0.822502,0.861186,0.899615,0.971651,0.984945,0.991741 /)
XCGA_LKT(48,12,1:6)=(/ 0.845457,0.814177,0.781460,0.745167,0.746777,0.742393 /)
XEXT_COEFF_550_LKT(48,12)=365.900000 !rg=0.448369 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(48,13,1:6)=(/ 279.360000,287.420000,300.080000,318.320000,340.100000,311.150000 /)
XPIZA_LKT(48,13,1:6)=(/ 0.811711,0.847876,0.884304,0.966396,0.982690,0.990504 /)
XCGA_LKT(48,13,1:6)=(/ 0.852590,0.824287,0.793287,0.752993,0.748800,0.742567 /)
XEXT_COEFF_550_LKT(48,13)=300.750000 !rg=0.448369 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(48,14,1:6)=(/ 231.590000,238.070000,246.250000,261.480000,278.440000,269.950000 /)
XPIZA_LKT(48,14,1:6)=(/ 0.801873,0.834679,0.870158,0.959751,0.979520,0.989001 /)
XCGA_LKT(48,14,1:6)=(/ 0.857727,0.834487,0.805140,0.763070,0.754763,0.744703 /)
XEXT_COEFF_550_LKT(48,14)=246.030000 !rg=0.448369 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(48,15,1:6)=(/ 191.870000,196.560000,203.380000,213.110000,227.830000,228.220000 /)
XPIZA_LKT(48,15,1:6)=(/ 0.793775,0.822968,0.854926,0.953107,0.974991,0.986716 /)
XCGA_LKT(48,15,1:6)=(/ 0.863363,0.841590,0.816963,0.770243,0.758890,0.744327 /)
XEXT_COEFF_550_LKT(48,15)=201.550000 !rg=0.448369 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(48,16,1:6)=(/ 158.540000,161.590000,166.420000,172.760000,186.210000,191.380000 /)
XPIZA_LKT(48,16,1:6)=(/ 0.784543,0.812021,0.842416,0.945324,0.969940,0.984216 /)
XCGA_LKT(48,16,1:6)=(/ 0.868317,0.850760,0.825703,0.785097,0.762483,0.745833 /)
XEXT_COEFF_550_LKT(48,16)=166.530000 !rg=0.448369 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(48,17,1:6)=(/ 131.330000,133.760000,136.760000,142.560000,151.270000,158.050000 /)
XPIZA_LKT(48,17,1:6)=(/ 0.776392,0.803257,0.830322,0.936810,0.965002,0.981811 /)
XCGA_LKT(48,17,1:6)=(/ 0.872097,0.856833,0.836233,0.793467,0.774623,0.751863 /)
XEXT_COEFF_550_LKT(48,17)=136.930000 !rg=0.448369 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(48,18,1:6)=(/ 108.760000,110.620000,113.410000,117.180000,123.660000,130.640000 /)
XPIZA_LKT(48,18,1:6)=(/ 0.767430,0.793328,0.819056,0.929580,0.958457,0.978016 /)
XCGA_LKT(48,18,1:6)=(/ 0.875907,0.861790,0.844420,0.800173,0.782110,0.755327 /)
XEXT_COEFF_550_LKT(48,18)=112.460000 !rg=0.448369 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(48,19,1:6)=(/ 89.994000,91.211000,93.095000,95.377000,101.370000,107.510000 /)
XPIZA_LKT(48,19,1:6)=(/ 0.757341,0.783988,0.809512,0.923711,0.950639,0.973897 /)
XCGA_LKT(48,19,1:6)=(/ 0.879457,0.867460,0.850693,0.811747,0.787987,0.758923 /)
XEXT_COEFF_550_LKT(48,19)=93.174000 !rg=0.448369 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(48,20,1:6)=(/ 74.502000,75.534000,76.717000,78.901000,82.652000,87.727000 /)
XPIZA_LKT(48,20,1:6)=(/ 0.746626,0.775532,0.799474,0.916772,0.942483,0.969495 /)
XCGA_LKT(48,20,1:6)=(/ 0.882513,0.871710,0.857820,0.818287,0.800283,0.769340 /)
XEXT_COEFF_550_LKT(48,20)=76.755000 !rg=0.448369 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(49,1,1:6)=(/ 1509.800000,1153.300000,2313.700000,2071.200000,672.190000,82.023000 /)
XPIZA_LKT(49,1,1:6)=(/ 0.928189,0.934502,0.979168,0.994418,0.992234,0.978136 /)
XCGA_LKT(49,1,1:6)=(/ 0.808493,0.552460,0.756477,0.767483,0.650747,0.181480 /)
XEXT_COEFF_550_LKT(49,1)=2289.900000 !rg=0.485743 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(49,2,1:6)=(/ 1424.300000,1249.100000,2062.600000,2036.200000,713.380000,96.120000 /)
XPIZA_LKT(49,2,1:6)=(/ 0.918021,0.939794,0.976494,0.994300,0.992780,0.980961 /)
XCGA_LKT(49,2,1:6)=(/ 0.775037,0.619503,0.740033,0.773207,0.655120,0.230473 /)
XEXT_COEFF_550_LKT(49,2)=2044.400000 !rg=0.485743 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(49,3,1:6)=(/ 1301.500000,1264.500000,1741.200000,1907.000000,782.980000,122.920000 /)
XPIZA_LKT(49,3,1:6)=(/ 0.913265,0.941394,0.971980,0.993847,0.993321,0.984462 /)
XCGA_LKT(49,3,1:6)=(/ 0.785007,0.674967,0.719377,0.771083,0.685190,0.332533 /)
XEXT_COEFF_550_LKT(49,3)=1721.900000 !rg=0.485743 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(49,4,1:6)=(/ 1158.700000,1186.600000,1464.800000,1692.000000,853.570000,164.500000 /)
XPIZA_LKT(49,4,1:6)=(/ 0.908110,0.938471,0.966700,0.992984,0.993717,0.987633 /)
XCGA_LKT(49,4,1:6)=(/ 0.786327,0.705013,0.708997,0.762407,0.716520,0.452960 /)
XEXT_COEFF_550_LKT(49,4)=1449.700000 !rg=0.485743 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(49,5,1:6)=(/ 1013.200000,1065.300000,1231.100000,1440.300000,895.550000,218.450000 /)
XPIZA_LKT(49,5,1:6)=(/ 0.901709,0.933395,0.961199,0.991634,0.993847,0.990037 /)
XCGA_LKT(49,5,1:6)=(/ 0.785140,0.727877,0.710363,0.749763,0.739533,0.550920 /)
XEXT_COEFF_550_LKT(49,5)=1223.800000 !rg=0.485743 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(49,6,1:6)=(/ 871.480000,921.330000,1040.400000,1197.200000,890.840000,277.550000 /)
XPIZA_LKT(49,6,1:6)=(/ 0.893888,0.927143,0.954075,0.989879,0.993650,0.991652 /)
XCGA_LKT(49,6,1:6)=(/ 0.799400,0.743763,0.716723,0.738293,0.752293,0.622497 /)
XEXT_COEFF_550_LKT(49,6)=1031.200000 !rg=0.485743 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(49,7,1:6)=(/ 749.650000,788.640000,869.920000,983.780000,838.250000,334.590000 /)
XPIZA_LKT(49,7,1:6)=(/ 0.880605,0.918811,0.946681,0.987585,0.993090,0.992675 /)
XCGA_LKT(49,7,1:6)=(/ 0.806140,0.755747,0.725880,0.729997,0.756280,0.674300 /)
XEXT_COEFF_550_LKT(49,7)=860.910000 !rg=0.485743 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(49,8,1:6)=(/ 631.160000,665.030000,722.570000,806.280000,753.500000,379.290000 /)
XPIZA_LKT(49,8,1:6)=(/ 0.868854,0.909021,0.938384,0.985010,0.992158,0.993218 /)
XCGA_LKT(49,8,1:6)=(/ 0.816337,0.772497,0.738073,0.727347,0.754813,0.709850 /)
XEXT_COEFF_550_LKT(49,8)=718.550000 !rg=0.485743 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(49,9,1:6)=(/ 534.720000,560.150000,595.330000,661.420000,655.760000,403.020000 /)
XPIZA_LKT(49,9,1:6)=(/ 0.853833,0.897241,0.930249,0.981365,0.990870,0.993347 /)
XCGA_LKT(49,9,1:6)=(/ 0.822420,0.785440,0.748293,0.727177,0.752130,0.731900 /)
XEXT_COEFF_550_LKT(49,9)=594.620000 !rg=0.485743 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(49,10,1:6)=(/ 447.590000,464.880000,492.570000,537.440000,556.230000,399.910000 /)
XPIZA_LKT(49,10,1:6)=(/ 0.840896,0.884606,0.919618,0.978517,0.989149,0.993063 /)
XCGA_LKT(49,10,1:6)=(/ 0.833907,0.796250,0.762530,0.732033,0.749067,0.742000 /)
XEXT_COEFF_550_LKT(49,10)=495.060000 !rg=0.485743 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(49,11,1:6)=(/ 372.480000,386.350000,404.960000,438.490000,464.640000,374.040000 /)
XPIZA_LKT(49,11,1:6)=(/ 0.828363,0.869540,0.907969,0.974333,0.986594,0.992361 /)
XCGA_LKT(49,11,1:6)=(/ 0.841140,0.809807,0.775947,0.739827,0.746983,0.743823 /)
XEXT_COEFF_550_LKT(49,11)=406.690000 !rg=0.485743 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(49,12,1:6)=(/ 309.870000,319.570000,333.930000,356.590000,382.040000,335.980000 /)
XPIZA_LKT(49,12,1:6)=(/ 0.816696,0.855085,0.893201,0.969567,0.984395,0.991328 /)
XCGA_LKT(49,12,1:6)=(/ 0.849723,0.819190,0.787430,0.747540,0.748110,0.743097 /)
XEXT_COEFF_550_LKT(49,12)=335.100000 !rg=0.485743 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(49,13,1:6)=(/ 257.150000,264.480000,273.910000,292.330000,312.500000,295.620000 /)
XPIZA_LKT(49,13,1:6)=(/ 0.806619,0.841293,0.878162,0.963261,0.981454,0.989990 /)
XCGA_LKT(49,13,1:6)=(/ 0.855147,0.830260,0.799403,0.757670,0.752000,0.744337 /)
XEXT_COEFF_550_LKT(49,13)=274.050000 !rg=0.485743 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(49,14,1:6)=(/ 213.090000,218.710000,226.360000,238.610000,254.920000,252.730000 /)
XPIZA_LKT(49,14,1:6)=(/ 0.797778,0.828444,0.862872,0.956916,0.977677,0.988203 /)
XCGA_LKT(49,14,1:6)=(/ 0.861403,0.837463,0.811670,0.763887,0.755927,0.744587 /)
XEXT_COEFF_550_LKT(49,14)=224.400000 !rg=0.485743 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(49,15,1:6)=(/ 176.090000,180.530000,186.120000,194.440000,209.090000,212.660000 /)
XPIZA_LKT(49,15,1:6)=(/ 0.789075,0.817700,0.848969,0.949106,0.972670,0.985814 /)
XCGA_LKT(49,15,1:6)=(/ 0.865590,0.846920,0.820173,0.778613,0.758177,0.744290 /)
XEXT_COEFF_550_LKT(49,15)=185.800000 !rg=0.485743 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(49,16,1:6)=(/ 145.660000,148.900000,152.400000,159.150000,170.690000,177.410000 /)
XPIZA_LKT(49,16,1:6)=(/ 0.780814,0.807711,0.836117,0.941151,0.967190,0.982809 /)
XCGA_LKT(49,16,1:6)=(/ 0.869910,0.854180,0.832597,0.788990,0.767287,0.748603 /)
XEXT_COEFF_550_LKT(49,16)=152.920000 !rg=0.485743 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(49,17,1:6)=(/ 120.980000,122.940000,125.930000,130.730000,138.140000,145.850000 /)
XPIZA_LKT(49,17,1:6)=(/ 0.771925,0.798129,0.825271,0.933297,0.962018,0.980228 /)
XCGA_LKT(49,17,1:6)=(/ 0.874597,0.859017,0.840693,0.794850,0.778273,0.753440 /)
XEXT_COEFF_550_LKT(49,17)=125.670000 !rg=0.485743 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(49,18,1:6)=(/ 99.938000,101.860000,103.980000,107.340000,113.760000,120.470000 /)
XPIZA_LKT(49,18,1:6)=(/ 0.762529,0.789371,0.814106,0.926064,0.954744,0.976012 /)
XCGA_LKT(49,18,1:6)=(/ 0.877383,0.865077,0.847247,0.806223,0.782373,0.754743 /)
XEXT_COEFF_550_LKT(49,18)=103.740000 !rg=0.485743 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(49,19,1:6)=(/ 82.730000,84.124000,85.467000,88.165000,93.050000,98.739000 /)
XPIZA_LKT(49,19,1:6)=(/ 0.752482,0.780061,0.804008,0.920044,0.946498,0.971809 /)
XCGA_LKT(49,19,1:6)=(/ 0.880510,0.869837,0.855227,0.815130,0.793347,0.763423 /)
XEXT_COEFF_550_LKT(49,19)=85.715000 !rg=0.485743 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(49,20,1:6)=(/ 68.721000,69.450000,70.699000,72.535000,75.518000,80.391000 /)
XPIZA_LKT(49,20,1:6)=(/ 0.741658,0.770126,0.795370,0.913972,0.939307,0.966936 /)
XCGA_LKT(49,20,1:6)=(/ 0.884463,0.873217,0.860847,0.819837,0.804457,0.772363 /)
XEXT_COEFF_550_LKT(49,20)=70.679000 !rg=0.485743 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(50,1,1:6)=(/ 1182.500000,1108.700000,1944.200000,2101.100000,760.820000,99.710000 /)
XPIZA_LKT(50,1,1:6)=(/ 0.915194,0.930989,0.974844,0.994417,0.993085,0.981492 /)
XCGA_LKT(50,1,1:6)=(/ 0.784587,0.585623,0.729213,0.781093,0.656810,0.215190 /)
XEXT_COEFF_550_LKT(50,1)=1925.400000 !rg=0.526233 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(50,2,1:6)=(/ 1241.300000,1184.500000,1730.600000,2004.500000,786.030000,114.940000 /)
XPIZA_LKT(50,2,1:6)=(/ 0.914318,0.937231,0.971795,0.994141,0.993331,0.983598 /)
XCGA_LKT(50,2,1:6)=(/ 0.795640,0.655350,0.713817,0.778847,0.675127,0.272847 /)
XEXT_COEFF_550_LKT(50,2)=1704.000000 !rg=0.526233 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(50,3,1:6)=(/ 1173.500000,1175.700000,1476.400000,1831.100000,846.690000,144.580000 /)
XPIZA_LKT(50,3,1:6)=(/ 0.910092,0.938242,0.967167,0.993525,0.993693,0.986320 /)
XCGA_LKT(50,3,1:6)=(/ 0.790547,0.699727,0.699380,0.772330,0.706053,0.382617 /)
XEXT_COEFF_550_LKT(50,3)=1466.800000 !rg=0.526233 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(50,4,1:6)=(/ 1053.900000,1083.500000,1271.900000,1585.000000,900.600000,190.120000 /)
XPIZA_LKT(50,4,1:6)=(/ 0.903849,0.935206,0.962011,0.992440,0.993931,0.988903 /)
XCGA_LKT(50,4,1:6)=(/ 0.793400,0.724243,0.698667,0.759727,0.732397,0.495533 /)
XEXT_COEFF_550_LKT(50,4)=1268.400000 !rg=0.526233 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(50,5,1:6)=(/ 930.240000,969.280000,1097.900000,1324.400000,918.620000,246.460000 /)
XPIZA_LKT(50,5,1:6)=(/ 0.894965,0.928244,0.956042,0.990875,0.993903,0.990871 /)
XCGA_LKT(50,5,1:6)=(/ 0.798227,0.735117,0.705733,0.745283,0.749907,0.582577 /)
XEXT_COEFF_550_LKT(50,5)=1086.200000 !rg=0.526233 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(50,6,1:6)=(/ 800.390000,836.190000,928.970000,1092.700000,887.770000,305.130000 /)
XPIZA_LKT(50,6,1:6)=(/ 0.887389,0.923380,0.949868,0.988902,0.993546,0.992183 /)
XCGA_LKT(50,6,1:6)=(/ 0.806683,0.755060,0.717143,0.734190,0.757920,0.645657 /)
XEXT_COEFF_550_LKT(50,6)=930.780000 !rg=0.526233 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(50,7,1:6)=(/ 685.840000,715.330000,782.080000,894.680000,814.360000,358.050000 /)
XPIZA_LKT(50,7,1:6)=(/ 0.873489,0.914606,0.942097,0.986518,0.992822,0.992983 /)
XCGA_LKT(50,7,1:6)=(/ 0.812833,0.769560,0.727687,0.728267,0.758570,0.690837 /)
XEXT_COEFF_550_LKT(50,7)=781.690000 !rg=0.526233 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(50,8,1:6)=(/ 581.900000,605.990000,652.370000,732.180000,718.600000,394.680000 /)
XPIZA_LKT(50,8,1:6)=(/ 0.861039,0.903953,0.933985,0.983640,0.991712,0.993346 /)
XCGA_LKT(50,8,1:6)=(/ 0.822553,0.779560,0.742490,0.726027,0.755260,0.720760 /)
XEXT_COEFF_550_LKT(50,8)=654.970000 !rg=0.526233 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(50,9,1:6)=(/ 491.070000,510.990000,545.020000,599.150000,617.120000,407.380000 /)
XPIZA_LKT(50,9,1:6)=(/ 0.846303,0.890685,0.923791,0.980005,0.990189,0.993305 /)
XCGA_LKT(50,9,1:6)=(/ 0.829990,0.788727,0.755463,0.726653,0.751320,0.737850 /)
XEXT_COEFF_550_LKT(50,9)=539.020000 !rg=0.526233 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(50,10,1:6)=(/ 411.460000,427.720000,448.570000,491.620000,516.920000,393.670000 /)
XPIZA_LKT(50,10,1:6)=(/ 0.833345,0.876257,0.914094,0.975902,0.988361,0.992854 /)
XCGA_LKT(50,10,1:6)=(/ 0.837187,0.804233,0.767233,0.734017,0.749050,0.744270 /)
XEXT_COEFF_550_LKT(50,10)=448.500000 !rg=0.526233 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(50,11,1:6)=(/ 342.540000,354.770000,370.700000,399.790000,427.920000,360.980000 /)
XPIZA_LKT(50,11,1:6)=(/ 0.821975,0.861649,0.901366,0.971164,0.985794,0.992002 /)
XCGA_LKT(50,11,1:6)=(/ 0.846483,0.813483,0.781670,0.740137,0.747833,0.744083 /)
XEXT_COEFF_550_LKT(50,11)=371.700000 !rg=0.526233 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(50,12,1:6)=(/ 284.930000,293.850000,305.000000,327.340000,350.920000,321.190000 /)
XPIZA_LKT(50,12,1:6)=(/ 0.810561,0.847712,0.886533,0.966144,0.983255,0.990847 /)
XCGA_LKT(50,12,1:6)=(/ 0.852443,0.825750,0.792860,0.751723,0.749763,0.744370 /)
XEXT_COEFF_550_LKT(50,12)=304.810000 !rg=0.526233 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(50,13,1:6)=(/ 236.440000,242.780000,251.770000,266.750000,286.520000,279.060000 /)
XPIZA_LKT(50,13,1:6)=(/ 0.801557,0.834415,0.870743,0.960230,0.979775,0.989299 /)
XCGA_LKT(50,13,1:6)=(/ 0.859300,0.833280,0.806043,0.758033,0.753057,0.744807 /)
XEXT_COEFF_550_LKT(50,13)=249.810000 !rg=0.526233 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(50,14,1:6)=(/ 195.520000,200.670000,207.390000,218.160000,234.630000,235.870000 /)
XPIZA_LKT(50,14,1:6)=(/ 0.793555,0.822993,0.855857,0.953136,0.974891,0.987237 /)
XCGA_LKT(50,14,1:6)=(/ 0.863700,0.843017,0.814687,0.771850,0.754940,0.744007 /)
XEXT_COEFF_550_LKT(50,14)=206.730000 !rg=0.526233 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(50,15,1:6)=(/ 162.020000,165.840000,170.600000,177.790000,192.050000,198.480000 /)
XPIZA_LKT(50,15,1:6)=(/ 0.784983,0.812494,0.842748,0.945726,0.970159,0.984390 /)
XCGA_LKT(50,15,1:6)=(/ 0.868867,0.850983,0.827350,0.784033,0.762547,0.746363 /)
XEXT_COEFF_550_LKT(50,15)=170.470000 !rg=0.526233 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(50,16,1:6)=(/ 134.090000,137.010000,140.340000,146.230000,154.910000,162.980000 /)
XPIZA_LKT(50,16,1:6)=(/ 0.776490,0.802750,0.830801,0.936507,0.965116,0.981629 /)
XCGA_LKT(50,16,1:6)=(/ 0.872723,0.856610,0.837013,0.790200,0.772633,0.750387 /)
XEXT_COEFF_550_LKT(50,16)=140.320000 !rg=0.526233 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(50,17,1:6)=(/ 111.120000,113.320000,115.670000,120.120000,126.790000,134.620000 /)
XPIZA_LKT(50,17,1:6)=(/ 0.767481,0.793952,0.819472,0.929345,0.959022,0.978481 /)
XCGA_LKT(50,17,1:6)=(/ 0.876110,0.862503,0.843343,0.800643,0.778837,0.753040 /)
XEXT_COEFF_550_LKT(50,17)=115.650000 !rg=0.526233 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(50,18,1:6)=(/ 92.122000,93.696000,95.555000,98.419000,104.200000,111.160000 /)
XPIZA_LKT(50,18,1:6)=(/ 0.758311,0.784863,0.809488,0.923483,0.951142,0.973829 /)
XCGA_LKT(50,18,1:6)=(/ 0.879657,0.867687,0.851803,0.811073,0.787610,0.758717 /)
XEXT_COEFF_550_LKT(50,18)=95.356000 !rg=0.526233 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(50,19,1:6)=(/ 76.255000,77.362000,78.782000,81.125000,84.529000,90.156000 /)
XPIZA_LKT(50,19,1:6)=(/ 0.747530,0.775205,0.799750,0.916567,0.943037,0.969581 /)
XCGA_LKT(50,19,1:6)=(/ 0.882697,0.871527,0.858307,0.816400,0.799170,0.767733 /)
XEXT_COEFF_550_LKT(50,19)=78.845000 !rg=0.526233 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(50,20,1:6)=(/ 63.190000,64.130000,64.982000,66.741000,69.529000,73.999000 /)
XPIZA_LKT(50,20,1:6)=(/ 0.736140,0.766301,0.790015,0.910733,0.935556,0.964046 /)
XCGA_LKT(50,20,1:6)=(/ 0.885757,0.875653,0.862747,0.824037,0.805483,0.772120 /)
XEXT_COEFF_550_LKT(50,20)=65.069000 !rg=0.526233 sigma=2.95 wvl=0.55
 
IF (LHOOK) CALL DR_HOOK('MODE_DUSTOPT:DUST_OPT_LKT_SET5',1,ZHOOK_HANDLE)
END SUBROUTINE DUST_OPT_LKT_SET5

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE DUST_OPT_LKT_SET6()

  USE MODD_DUST_OPT_LKT
  
  IMPLICIT NONE


REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_DUSTOPT:DUST_OPT_LKT_SET6',0,ZHOOK_HANDLE)
XEXT_COEFF_WVL_LKT(51,1,1:6)=(/ 1070.100000,1162.800000,1569.700000,2030.200000,820.670000,119.650000 /)
XPIZA_LKT(51,1,1:6)=(/ 0.902505,0.934400,0.968857,0.994241,0.993670,0.984126 /)
XCGA_LKT(51,1,1:6)=(/ 0.781037,0.663077,0.693750,0.787200,0.665610,0.256210 /)
XEXT_COEFF_550_LKT(51,1)=1536.100000 !rg=0.570098 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(51,2,1:6)=(/ 1128.700000,1170.700000,1412.800000,1932.100000,852.240000,135.810000 /)
XPIZA_LKT(51,2,1:6)=(/ 0.910202,0.937154,0.965524,0.993857,0.993748,0.985668 /)
XCGA_LKT(51,2,1:6)=(/ 0.797117,0.702913,0.679930,0.781053,0.696900,0.323037 /)
XEXT_COEFF_550_LKT(51,2)=1397.700000 !rg=0.570098 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(51,3,1:6)=(/ 1061.300000,1101.900000,1257.400000,1721.700000,903.980000,168.850000 /)
XPIZA_LKT(51,3,1:6)=(/ 0.906660,0.934863,0.961287,0.993047,0.993980,0.987837 /)
XCGA_LKT(51,3,1:6)=(/ 0.798900,0.723187,0.679757,0.770087,0.725117,0.434823 /)
XEXT_COEFF_550_LKT(51,3)=1238.800000 !rg=0.570098 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(51,4,1:6)=(/ 958.470000,997.150000,1116.000000,1458.600000,938.200000,217.980000 /)
XPIZA_LKT(51,4,1:6)=(/ 0.899939,0.930776,0.956962,0.991710,0.994074,0.989962 /)
XCGA_LKT(51,4,1:6)=(/ 0.800920,0.738720,0.692110,0.753837,0.746117,0.535663 /)
XEXT_COEFF_550_LKT(51,4)=1106.400000 !rg=0.570098 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(51,5,1:6)=(/ 849.080000,881.690000,973.840000,1205.600000,930.160000,275.500000 /)
XPIZA_LKT(51,5,1:6)=(/ 0.889071,0.924837,0.951079,0.989916,0.993892,0.991568 /)
XCGA_LKT(51,5,1:6)=(/ 0.804353,0.753737,0.704977,0.739217,0.758350,0.611497 /)
XEXT_COEFF_550_LKT(51,5)=972.680000 !rg=0.570098 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(51,6,1:6)=(/ 735.340000,769.920000,833.830000,991.660000,874.500000,332.070000 /)
XPIZA_LKT(51,6,1:6)=(/ 0.879661,0.916996,0.945775,0.987584,0.993366,0.992621 /)
XCGA_LKT(51,6,1:6)=(/ 0.810153,0.764553,0.721540,0.727800,0.761977,0.666447 /)
XEXT_COEFF_550_LKT(51,6)=833.940000 !rg=0.570098 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(51,7,1:6)=(/ 627.740000,655.620000,704.740000,811.000000,784.240000,379.020000 /)
XPIZA_LKT(51,7,1:6)=(/ 0.867216,0.908461,0.938392,0.985017,0.992463,0.993221 /)
XCGA_LKT(51,7,1:6)=(/ 0.818890,0.777273,0.736617,0.724027,0.759697,0.705227 /)
XEXT_COEFF_550_LKT(51,7)=706.110000 !rg=0.570098 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(51,8,1:6)=(/ 535.580000,557.940000,591.000000,666.150000,680.490000,406.000000 /)
XPIZA_LKT(51,8,1:6)=(/ 0.852387,0.896758,0.929523,0.981676,0.991202,0.993410 /)
XCGA_LKT(51,8,1:6)=(/ 0.825793,0.788443,0.748043,0.724567,0.755427,0.729743 /)
XEXT_COEFF_550_LKT(51,8)=591.600000 !rg=0.570098 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(51,9,1:6)=(/ 449.740000,466.300000,494.750000,542.470000,575.410000,407.220000 /)
XPIZA_LKT(51,9,1:6)=(/ 0.838452,0.883510,0.918288,0.978643,0.989463,0.993203 /)
XCGA_LKT(51,9,1:6)=(/ 0.834747,0.799957,0.758093,0.732000,0.750553,0.742343 /)
XEXT_COEFF_550_LKT(51,9)=494.030000 !rg=0.570098 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(51,10,1:6)=(/ 378.120000,391.360000,411.810000,446.130000,478.810000,384.230000 /)
XPIZA_LKT(51,10,1:6)=(/ 0.826776,0.868968,0.906973,0.974064,0.987277,0.992575 /)
XCGA_LKT(51,10,1:6)=(/ 0.843303,0.807433,0.774843,0.733883,0.747620,0.745307 /)
XEXT_COEFF_550_LKT(51,10)=407.200000 !rg=0.570098 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(51,11,1:6)=(/ 315.170000,325.820000,339.200000,365.780000,394.380000,346.950000 /)
XPIZA_LKT(51,11,1:6)=(/ 0.815484,0.854432,0.893694,0.968707,0.984745,0.991599 /)
XCGA_LKT(51,11,1:6)=(/ 0.849527,0.821137,0.785417,0.745707,0.747980,0.745050 /)
XEXT_COEFF_550_LKT(51,11)=338.270000 !rg=0.570098 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(51,12,1:6)=(/ 262.140000,269.680000,280.290000,298.000000,322.110000,305.770000 /)
XPIZA_LKT(51,12,1:6)=(/ 0.805992,0.840770,0.878549,0.963844,0.981686,0.990242 /)
XCGA_LKT(51,12,1:6)=(/ 0.856807,0.828953,0.800067,0.752213,0.749930,0.744673 /)
XEXT_COEFF_550_LKT(51,12)=277.590000 !rg=0.570098 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(51,13,1:6)=(/ 216.930000,222.930000,230.630000,243.950000,262.920000,261.540000 /)
XPIZA_LKT(51,13,1:6)=(/ 0.797268,0.828364,0.863190,0.956643,0.977478,0.988591 /)
XCGA_LKT(51,13,1:6)=(/ 0.861747,0.839297,0.809133,0.765657,0.751823,0.745167 /)
XEXT_COEFF_550_LKT(51,13)=229.810000 !rg=0.570098 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(51,14,1:6)=(/ 179.790000,184.470000,190.010000,198.620000,215.350000,220.280000 /)
XPIZA_LKT(51,14,1:6)=(/ 0.788311,0.816645,0.849345,0.949623,0.972851,0.986006 /)
XCGA_LKT(51,14,1:6)=(/ 0.867197,0.847223,0.821960,0.777370,0.758780,0.745130 /)
XEXT_COEFF_550_LKT(51,14)=190.000000 !rg=0.570098 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(51,15,1:6)=(/ 149.380000,152.610000,156.460000,163.100000,174.140000,182.720000 /)
XPIZA_LKT(51,15,1:6)=(/ 0.781087,0.807643,0.836856,0.941767,0.968747,0.983752 /)
XCGA_LKT(51,15,1:6)=(/ 0.871357,0.853860,0.832430,0.786583,0.767213,0.748607 /)
XEXT_COEFF_550_LKT(51,15)=156.640000 !rg=0.570098 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(51,16,1:6)=(/ 123.530000,126.030000,129.010000,133.970000,141.560000,150.520000 /)
XPIZA_LKT(51,16,1:6)=(/ 0.771998,0.797913,0.824871,0.934019,0.963178,0.980806 /)
XCGA_LKT(51,16,1:6)=(/ 0.874340,0.860470,0.839900,0.796833,0.775813,0.751627 /)
XEXT_COEFF_550_LKT(51,16)=128.670000 !rg=0.570098 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(51,17,1:6)=(/ 102.320000,104.010000,106.460000,109.610000,116.420000,124.250000 /)
XPIZA_LKT(51,17,1:6)=(/ 0.762519,0.789162,0.814473,0.926288,0.955711,0.976477 /)
XCGA_LKT(51,17,1:6)=(/ 0.878350,0.864703,0.848127,0.804640,0.783603,0.755680 /)
XEXT_COEFF_550_LKT(51,17)=106.450000 !rg=0.570098 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(51,18,1:6)=(/ 84.945000,86.266000,87.830000,90.440000,94.811000,101.230000 /)
XPIZA_LKT(51,18,1:6)=(/ 0.753617,0.780516,0.804685,0.920277,0.948101,0.972742 /)
XCGA_LKT(51,18,1:6)=(/ 0.881537,0.869487,0.855323,0.813383,0.793463,0.762633 /)
XEXT_COEFF_550_LKT(51,18)=87.738000 !rg=0.570098 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(51,19,1:6)=(/ 70.266000,71.314000,72.448000,74.430000,77.372000,82.611000 /)
XPIZA_LKT(51,19,1:6)=(/ 0.742395,0.771162,0.794822,0.914591,0.940306,0.968030 /)
XCGA_LKT(51,19,1:6)=(/ 0.884027,0.874027,0.860447,0.821107,0.802590,0.769847 /)
XEXT_COEFF_550_LKT(51,19)=72.376000 !rg=0.570098 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(51,20,1:6)=(/ 58.218000,58.929000,59.948000,61.107000,63.756000,68.052000 /)
XPIZA_LKT(51,20,1:6)=(/ 0.730420,0.760925,0.786242,0.908240,0.932877,0.961266 /)
XCGA_LKT(51,20,1:6)=(/ 0.887510,0.877307,0.865960,0.827137,0.810157,0.776567 /)
XEXT_COEFF_550_LKT(51,20)=59.895000 !rg=0.570098 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(52,1,1:6)=(/ 1089.000000,1245.700000,1224.500000,1958.700000,868.130000,141.400000 /)
XPIZA_LKT(52,1,1:6)=(/ 0.900207,0.938811,0.960342,0.993830,0.993947,0.986162 /)
XCGA_LKT(52,1,1:6)=(/ 0.782047,0.736940,0.644703,0.785907,0.693060,0.306300 /)
XEXT_COEFF_550_LKT(52,1)=1198.300000 !rg=0.617619 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(52,2,1:6)=(/ 1046.000000,1155.300000,1151.500000,1818.600000,913.730000,158.790000 /)
XPIZA_LKT(52,2,1:6)=(/ 0.905089,0.936083,0.957642,0.993412,0.994055,0.987298 /)
XCGA_LKT(52,2,1:6)=(/ 0.797837,0.741670,0.645890,0.779190,0.720530,0.380983 /)
XEXT_COEFF_550_LKT(52,2)=1126.000000 !rg=0.617619 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(52,3,1:6)=(/ 981.260000,1029.900000,1076.000000,1582.700000,952.950000,195.940000 /)
XPIZA_LKT(52,3,1:6)=(/ 0.899915,0.931634,0.954738,0.992380,0.994190,0.989091 /)
XCGA_LKT(52,3,1:6)=(/ 0.805807,0.742710,0.662303,0.764240,0.742233,0.486390 /)
XEXT_COEFF_550_LKT(52,3)=1070.300000 !rg=0.617619 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(52,4,1:6)=(/ 885.000000,917.120000,985.500000,1318.900000,964.520000,247.780000 /)
XPIZA_LKT(52,4,1:6)=(/ 0.892934,0.926293,0.951374,0.990795,0.994147,0.990847 /)
XCGA_LKT(52,4,1:6)=(/ 0.808793,0.750187,0.687590,0.746137,0.757667,0.572460 /)
XEXT_COEFF_550_LKT(52,4)=983.660000 !rg=0.617619 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(52,5,1:6)=(/ 774.350000,807.840000,872.000000,1084.800000,929.580000,305.010000 /)
XPIZA_LKT(52,5,1:6)=(/ 0.884180,0.920053,0.946727,0.988714,0.993804,0.992150 /)
XCGA_LKT(52,5,1:6)=(/ 0.811740,0.764327,0.710700,0.730723,0.764813,0.637640 /)
XEXT_COEFF_550_LKT(52,5)=869.120000 !rg=0.617619 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(52,6,1:6)=(/ 677.050000,704.370000,756.420000,891.140000,851.010000,357.750000 /)
XPIZA_LKT(52,6,1:6)=(/ 0.871802,0.910408,0.939975,0.986221,0.993106,0.992974 /)
XCGA_LKT(52,6,1:6)=(/ 0.818863,0.768490,0.725403,0.723473,0.764340,0.685050 /)
XEXT_COEFF_550_LKT(52,6)=750.310000 !rg=0.617619 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(52,7,1:6)=(/ 578.490000,600.520000,640.010000,731.280000,747.350000,396.940000 /)
XPIZA_LKT(52,7,1:6)=(/ 0.859125,0.901439,0.932899,0.983321,0.992020,0.993389 /)
XCGA_LKT(52,7,1:6)=(/ 0.826207,0.781780,0.741547,0.721660,0.759683,0.717743 /)
XEXT_COEFF_550_LKT(52,7)=641.810000 !rg=0.617619 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(52,8,1:6)=(/ 491.540000,510.660000,541.150000,601.830000,638.920000,413.060000 /)
XPIZA_LKT(52,8,1:6)=(/ 0.845453,0.889413,0.923211,0.979778,0.990540,0.993409 /)
XCGA_LKT(52,8,1:6)=(/ 0.833420,0.791177,0.755287,0.723040,0.754097,0.737063 /)
XEXT_COEFF_550_LKT(52,8)=535.080000 !rg=0.617619 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(52,9,1:6)=(/ 411.850000,427.270000,449.400000,491.820000,534.830000,403.130000 /)
XPIZA_LKT(52,9,1:6)=(/ 0.832187,0.876553,0.912961,0.976518,0.988385,0.993027 /)
XCGA_LKT(52,9,1:6)=(/ 0.840840,0.807050,0.769617,0.732707,0.749200,0.745447 /)
XEXT_COEFF_550_LKT(52,9)=450.310000 !rg=0.617619 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(52,10,1:6)=(/ 346.300000,358.460000,375.320000,404.860000,441.770000,372.110000 /)
XPIZA_LKT(52,10,1:6)=(/ 0.820603,0.861279,0.900472,0.971757,0.986236,0.992263 /)
XCGA_LKT(52,10,1:6)=(/ 0.846890,0.816583,0.777523,0.741600,0.746640,0.746543 /)
XEXT_COEFF_550_LKT(52,10)=374.320000 !rg=0.617619 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(52,11,1:6)=(/ 289.390000,297.940000,311.220000,331.450000,362.430000,332.160000 /)
XPIZA_LKT(52,11,1:6)=(/ 0.810074,0.847752,0.886228,0.966972,0.983250,0.991137 /)
XCGA_LKT(52,11,1:6)=(/ 0.854103,0.824820,0.793873,0.747810,0.747597,0.745633 /)
XEXT_COEFF_550_LKT(52,11)=308.670000 !rg=0.617619 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(52,12,1:6)=(/ 240.300000,247.270000,256.450000,271.760000,295.860000,288.710000 /)
XPIZA_LKT(52,12,1:6)=(/ 0.800997,0.834462,0.870988,0.960540,0.979548,0.989683 /)
XCGA_LKT(52,12,1:6)=(/ 0.859493,0.835530,0.803150,0.760650,0.748830,0.746097 /)
XEXT_COEFF_550_LKT(52,12)=255.580000 !rg=0.617619 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(52,13,1:6)=(/ 199.440000,204.580000,211.450000,221.730000,241.510000,244.690000 /)
XPIZA_LKT(52,13,1:6)=(/ 0.791716,0.821804,0.856258,0.953499,0.975517,0.987383 /)
XCGA_LKT(52,13,1:6)=(/ 0.865493,0.843577,0.816833,0.771043,0.755080,0.744793 /)
XEXT_COEFF_550_LKT(52,13)=211.290000 !rg=0.617619 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(52,14,1:6)=(/ 165.800000,169.150000,174.190000,181.570000,195.560000,203.920000 /)
XPIZA_LKT(52,14,1:6)=(/ 0.784558,0.812519,0.842414,0.946406,0.971712,0.985435 /)
XCGA_LKT(52,14,1:6)=(/ 0.869610,0.850897,0.827073,0.781453,0.762227,0.747010 /)
XEXT_COEFF_550_LKT(52,14)=174.360000 !rg=0.617619 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(52,15,1:6)=(/ 137.310000,140.330000,143.730000,149.880000,158.980000,168.630000 /)
XPIZA_LKT(52,15,1:6)=(/ 0.776359,0.802281,0.831199,0.937582,0.966785,0.982613 /)
XCGA_LKT(52,15,1:6)=(/ 0.872903,0.857620,0.836043,0.791593,0.772083,0.749980 /)
XEXT_COEFF_550_LKT(52,15)=143.320000 !rg=0.617619 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(52,16,1:6)=(/ 113.600000,115.650000,118.580000,122.270000,130.200000,139.060000 /)
XPIZA_LKT(52,16,1:6)=(/ 0.768309,0.793700,0.819435,0.930792,0.959945,0.978995 /)
XCGA_LKT(52,16,1:6)=(/ 0.876293,0.862673,0.845100,0.800133,0.779587,0.753523 /)
XEXT_COEFF_550_LKT(52,16)=118.180000 !rg=0.617619 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(52,17,1:6)=(/ 94.281000,95.536000,97.597000,100.310000,106.180000,113.680000 /)
XPIZA_LKT(52,17,1:6)=(/ 0.758691,0.784595,0.809308,0.923912,0.953104,0.975621 /)
XCGA_LKT(52,17,1:6)=(/ 0.879767,0.867973,0.851177,0.809940,0.787370,0.758323 /)
XEXT_COEFF_550_LKT(52,17)=97.791000 !rg=0.617619 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(52,18,1:6)=(/ 78.110000,79.338000,80.823000,83.085000,86.538000,92.770000 /)
XPIZA_LKT(52,18,1:6)=(/ 0.748263,0.775907,0.800215,0.917178,0.945202,0.971099 /)
XCGA_LKT(52,18,1:6)=(/ 0.882670,0.872000,0.857790,0.817363,0.799053,0.766703 /)
XEXT_COEFF_550_LKT(52,18)=80.558000 !rg=0.617619 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(52,19,1:6)=(/ 64.701000,65.534000,66.723000,68.184000,71.312000,76.137000 /)
XPIZA_LKT(52,19,1:6)=(/ 0.737520,0.766675,0.790814,0.911508,0.936757,0.965180 /)
XCGA_LKT(52,19,1:6)=(/ 0.885600,0.875533,0.863800,0.823907,0.806947,0.773240 /)
XEXT_COEFF_550_LKT(52,19)=66.637000 !rg=0.617619 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(52,20,1:6)=(/ 53.664000,54.200000,55.033000,56.093000,58.365000,62.086000 /)
XPIZA_LKT(52,20,1:6)=(/ 0.725363,0.756131,0.781008,0.906107,0.930363,0.959368 /)
XCGA_LKT(52,20,1:6)=(/ 0.888793,0.879407,0.868097,0.830727,0.813873,0.780297 /)
XEXT_COEFF_550_LKT(52,20)=55.167000 !rg=0.617619 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(53,1,1:6)=(/ 1095.900000,1249.100000,965.760000,1819.500000,936.570000,164.520000 /)
XPIZA_LKT(53,1,1:6)=(/ 0.907673,0.940059,0.948529,0.993469,0.994114,0.987708 /)
XCGA_LKT(53,1,1:6)=(/ 0.817427,0.783713,0.582280,0.782163,0.732503,0.367050 /)
XEXT_COEFF_550_LKT(53,1)=940.400000 !rg=0.669101 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(53,2,1:6)=(/ 982.090000,1091.800000,958.520000,1667.100000,969.820000,184.400000 /)
XPIZA_LKT(53,2,1:6)=(/ 0.899269,0.933780,0.948603,0.992752,0.994286,0.988601 /)
XCGA_LKT(53,2,1:6)=(/ 0.808130,0.767310,0.612577,0.772640,0.742607,0.444473 /)
XEXT_COEFF_550_LKT(53,2)=946.540000 !rg=0.669101 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(53,3,1:6)=(/ 902.240000,960.910000,942.210000,1423.300000,991.500000,225.880000 /)
XPIZA_LKT(53,3,1:6)=(/ 0.892244,0.925747,0.948643,0.991443,0.994327,0.990144 /)
XCGA_LKT(53,3,1:6)=(/ 0.802087,0.763390,0.660417,0.753843,0.757043,0.534427 /)
XEXT_COEFF_550_LKT(53,3)=934.800000 !rg=0.669101 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(53,4,1:6)=(/ 812.640000,851.650000,880.830000,1177.700000,977.890000,279.030000 /)
XPIZA_LKT(53,4,1:6)=(/ 0.885154,0.919626,0.946158,0.989563,0.994147,0.991590 /)
XCGA_LKT(53,4,1:6)=(/ 0.807910,0.766717,0.692133,0.734670,0.767003,0.605493 /)
XEXT_COEFF_550_LKT(53,4)=876.090000 !rg=0.669101 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(53,5,1:6)=(/ 714.840000,739.680000,785.500000,967.710000,916.210000,334.260000 /)
XPIZA_LKT(53,5,1:6)=(/ 0.876548,0.914704,0.941631,0.987346,0.993642,0.992632 /)
XCGA_LKT(53,5,1:6)=(/ 0.818800,0.771140,0.715387,0.723500,0.769337,0.661050 /)
XEXT_COEFF_550_LKT(53,5)=788.500000 !rg=0.669101 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(53,6,1:6)=(/ 619.160000,647.580000,683.650000,802.400000,818.440000,381.280000 /)
XPIZA_LKT(53,6,1:6)=(/ 0.864741,0.903794,0.934379,0.984394,0.992767,0.993252 /)
XCGA_LKT(53,6,1:6)=(/ 0.821817,0.781883,0.729967,0.719050,0.765640,0.701403 /)
XEXT_COEFF_550_LKT(53,6)=679.650000 !rg=0.669101 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(53,7,1:6)=(/ 531.080000,552.790000,581.150000,662.160000,705.770000,411.000000 /)
XPIZA_LKT(53,7,1:6)=(/ 0.850890,0.894526,0.927081,0.981140,0.991499,0.993491 /)
XCGA_LKT(53,7,1:6)=(/ 0.829267,0.792613,0.746170,0.720097,0.759277,0.728297 /)
XEXT_COEFF_550_LKT(53,7)=579.740000 !rg=0.669101 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(53,8,1:6)=(/ 449.770000,467.950000,492.270000,545.510000,594.580000,415.520000 /)
XPIZA_LKT(53,8,1:6)=(/ 0.837840,0.881754,0.916852,0.977544,0.989825,0.993346 /)
XCGA_LKT(53,8,1:6)=(/ 0.837237,0.802000,0.758057,0.727020,0.752907,0.742860 /)
XEXT_COEFF_550_LKT(53,8)=489.600000 !rg=0.669101 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(53,9,1:6)=(/ 379.550000,392.070000,409.880000,446.340000,491.840000,395.050000 /)
XPIZA_LKT(53,9,1:6)=(/ 0.825518,0.869041,0.907256,0.974327,0.987684,0.992800 /)
XCGA_LKT(53,9,1:6)=(/ 0.846013,0.811507,0.776270,0.733797,0.748480,0.747573 /)
XEXT_COEFF_550_LKT(53,9)=412.240000 !rg=0.669101 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(53,10,1:6)=(/ 317.810000,328.420000,342.650000,366.740000,407.350000,358.740000 /)
XPIZA_LKT(53,10,1:6)=(/ 0.814211,0.854316,0.893402,0.969903,0.984635,0.991875 /)
XCGA_LKT(53,10,1:6)=(/ 0.852593,0.822767,0.788167,0.745650,0.746500,0.747290 /)
XEXT_COEFF_550_LKT(53,10)=342.620000 !rg=0.669101 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(53,11,1:6)=(/ 266.310000,272.890000,283.830000,300.480000,331.750000,315.790000 /)
XPIZA_LKT(53,11,1:6)=(/ 0.803700,0.840374,0.879187,0.964932,0.981767,0.990645 /)
XCGA_LKT(53,11,1:6)=(/ 0.857773,0.832580,0.797600,0.757087,0.746147,0.747180 /)
XEXT_COEFF_550_LKT(53,11)=283.790000 !rg=0.669101 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(53,12,1:6)=(/ 220.870000,227.030000,234.860000,247.000000,271.330000,271.550000 /)
XPIZA_LKT(53,12,1:6)=(/ 0.795675,0.827407,0.863792,0.957897,0.977640,0.988767 /)
XCGA_LKT(53,12,1:6)=(/ 0.863553,0.840457,0.811710,0.766430,0.751287,0.745913 /)
XEXT_COEFF_550_LKT(53,12)=234.800000 !rg=0.669101 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(53,13,1:6)=(/ 183.810000,187.510000,193.490000,202.240000,219.460000,226.850000 /)
XPIZA_LKT(53,13,1:6)=(/ 0.788376,0.817066,0.849166,0.950961,0.974419,0.986819 /)
XCGA_LKT(53,13,1:6)=(/ 0.867927,0.848153,0.821923,0.776013,0.757860,0.745660 /)
XEXT_COEFF_550_LKT(53,13)=193.840000 !rg=0.669101 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(53,14,1:6)=(/ 152.550000,155.510000,159.500000,166.800000,178.460000,188.840000 /)
XPIZA_LKT(53,14,1:6)=(/ 0.780170,0.807264,0.837411,0.942287,0.969693,0.984499 /)
XCGA_LKT(53,14,1:6)=(/ 0.871427,0.855080,0.831883,0.787033,0.767763,0.749483 /)
XEXT_COEFF_550_LKT(53,14)=159.350000 !rg=0.669101 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(53,15,1:6)=(/ 126.590000,129.060000,132.120000,137.080000,145.760000,155.730000 /)
XPIZA_LKT(53,15,1:6)=(/ 0.772522,0.797923,0.824843,0.934580,0.963645,0.981080 /)
XCGA_LKT(53,15,1:6)=(/ 0.875047,0.860297,0.841107,0.794150,0.775057,0.750847 /)
XEXT_COEFF_550_LKT(53,15)=131.110000 !rg=0.669101 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(53,16,1:6)=(/ 104.860000,106.370000,108.630000,111.550000,119.440000,127.830000 /)
XPIZA_LKT(53,16,1:6)=(/ 0.763453,0.788479,0.814532,0.928051,0.956688,0.977162 /)
XCGA_LKT(53,16,1:6)=(/ 0.878297,0.866220,0.847930,0.807070,0.779997,0.752800 /)
XEXT_COEFF_550_LKT(53,16)=108.600000 !rg=0.669101 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(53,17,1:6)=(/ 86.838000,88.114000,89.526000,92.344000,97.269000,104.220000 /)
XPIZA_LKT(53,17,1:6)=(/ 0.753955,0.780859,0.804442,0.920613,0.948922,0.973554 /)
XCGA_LKT(53,17,1:6)=(/ 0.881237,0.870457,0.855537,0.814140,0.793507,0.762973 /)
XEXT_COEFF_550_LKT(53,17)=89.650000 !rg=0.669101 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(53,18,1:6)=(/ 72.045000,73.040000,74.271000,76.243000,79.671000,85.260000 /)
XPIZA_LKT(53,18,1:6)=(/ 0.743570,0.771395,0.795204,0.914653,0.940776,0.968494 /)
XCGA_LKT(53,18,1:6)=(/ 0.884303,0.873910,0.861123,0.819830,0.802327,0.769120 /)
XEXT_COEFF_550_LKT(53,18)=73.898000 !rg=0.669101 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(53,19,1:6)=(/ 59.680000,60.331000,61.249000,62.322000,65.522000,69.934000 /)
XPIZA_LKT(53,19,1:6)=(/ 0.731877,0.761365,0.786453,0.909294,0.933596,0.962473 /)
XCGA_LKT(53,19,1:6)=(/ 0.887277,0.877880,0.865783,0.828773,0.807673,0.773707 /)
XEXT_COEFF_550_LKT(53,19)=61.275000 !rg=0.669101 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(53,20,1:6)=(/ 49.459000,49.998000,50.563000,51.668000,53.617000,56.930000 /)
XPIZA_LKT(53,20,1:6)=(/ 0.719652,0.751685,0.776620,0.903078,0.926688,0.956004 /)
XCGA_LKT(53,20,1:6)=(/ 0.890210,0.881257,0.870823,0.833737,0.818403,0.786037 /)
XEXT_COEFF_550_LKT(53,20)=50.583000 !rg=0.669101 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(54,1,1:6)=(/ 986.360000,1138.300000,772.230000,1647.700000,1013.700000,189.210000 /)
XPIZA_LKT(54,1,1:6)=(/ 0.897858,0.937675,0.937647,0.992570,0.994391,0.988872 /)
XCGA_LKT(54,1,1:6)=(/ 0.816437,0.798877,0.549613,0.774817,0.760017,0.438437 /)
XEXT_COEFF_550_LKT(54,1)=778.990000 !rg=0.724875 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(54,2,1:6)=(/ 898.900000,996.230000,831.920000,1483.900000,1016.600000,213.600000 /)
XPIZA_LKT(54,2,1:6)=(/ 0.895142,0.925175,0.941687,0.991798,0.994455,0.989686 /)
XCGA_LKT(54,2,1:6)=(/ 0.807410,0.782427,0.613243,0.760597,0.760330,0.508183 /)
XEXT_COEFF_550_LKT(54,2)=830.390000 !rg=0.724875 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(54,3,1:6)=(/ 826.850000,880.990000,855.600000,1249.100000,1017.400000,258.400000 /)
XPIZA_LKT(54,3,1:6)=(/ 0.886022,0.919866,0.942975,0.990225,0.994391,0.991037 /)
XCGA_LKT(54,3,1:6)=(/ 0.811910,0.767137,0.667007,0.739637,0.769250,0.576760 /)
XEXT_COEFF_550_LKT(54,3)=846.640000 !rg=0.724875 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(54,4,1:6)=(/ 744.020000,781.060000,803.040000,1033.800000,977.050000,311.070000 /)
XPIZA_LKT(54,4,1:6)=(/ 0.878070,0.914390,0.940344,0.988136,0.994067,0.992213 /)
XCGA_LKT(54,4,1:6)=(/ 0.816943,0.769863,0.701183,0.722727,0.774017,0.634733 /)
XEXT_COEFF_550_LKT(54,4)=793.660000 !rg=0.724875 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(54,5,1:6)=(/ 656.740000,684.380000,711.710000,862.620000,890.720000,362.490000 /)
XPIZA_LKT(54,5,1:6)=(/ 0.868268,0.907608,0.936674,0.985552,0.993392,0.993026 /)
XCGA_LKT(54,5,1:6)=(/ 0.820723,0.783973,0.722983,0.714493,0.772053,0.681890 /)
XEXT_COEFF_550_LKT(54,5)=710.380000 !rg=0.724875 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(54,6,1:6)=(/ 567.320000,590.280000,624.650000,715.090000,778.570000,401.950000 /)
XPIZA_LKT(54,6,1:6)=(/ 0.855905,0.898716,0.928069,0.982858,0.992318,0.993460 /)
XCGA_LKT(54,6,1:6)=(/ 0.829563,0.788730,0.741213,0.715503,0.765157,0.715730 /)
XEXT_COEFF_550_LKT(54,6)=617.510000 !rg=0.724875 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(54,7,1:6)=(/ 487.360000,504.970000,532.910000,593.740000,660.670000,420.810000 /)
XPIZA_LKT(54,7,1:6)=(/ 0.842677,0.888338,0.920539,0.979630,0.990818,0.993528 /)
XCGA_LKT(54,7,1:6)=(/ 0.836167,0.796080,0.756257,0.719323,0.757153,0.737110 /)
XEXT_COEFF_550_LKT(54,7)=526.610000 !rg=0.724875 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(54,8,1:6)=(/ 412.360000,427.760000,449.610000,490.770000,550.850000,413.690000 /)
XPIZA_LKT(54,8,1:6)=(/ 0.829845,0.874507,0.910448,0.975851,0.988767,0.993212 /)
XCGA_LKT(54,8,1:6)=(/ 0.843963,0.808623,0.769670,0.729207,0.750473,0.747113 /)
XEXT_COEFF_550_LKT(54,8)=447.490000 !rg=0.724875 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(54,9,1:6)=(/ 348.910000,360.870000,374.470000,408.070000,452.670000,384.210000 /)
XPIZA_LKT(54,9,1:6)=(/ 0.818356,0.860780,0.900756,0.971244,0.986586,0.992515 /)
XCGA_LKT(54,9,1:6)=(/ 0.849287,0.819317,0.781057,0.737913,0.747883,0.749343 /)
XEXT_COEFF_550_LKT(54,9)=374.450000 !rg=0.724875 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(54,10,1:6)=(/ 293.070000,301.250000,312.770000,333.840000,370.750000,343.400000 /)
XPIZA_LKT(54,10,1:6)=(/ 0.808370,0.847874,0.886806,0.967544,0.983872,0.991451 /)
XCGA_LKT(54,10,1:6)=(/ 0.856627,0.827490,0.794973,0.748367,0.747050,0.748143 /)
XEXT_COEFF_550_LKT(54,10)=314.030000 !rg=0.724875 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(54,11,1:6)=(/ 244.340000,251.070000,259.010000,274.660000,304.590000,299.130000 /)
XPIZA_LKT(54,11,1:6)=(/ 0.799861,0.834299,0.871911,0.961895,0.979517,0.990014 /)
XCGA_LKT(54,11,1:6)=(/ 0.861047,0.837883,0.807390,0.761053,0.748370,0.747787 /)
XEXT_COEFF_550_LKT(54,11)=259.970000 !rg=0.724875 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(54,12,1:6)=(/ 203.790000,208.290000,214.930000,225.410000,246.080000,252.590000 /)
XPIZA_LKT(54,12,1:6)=(/ 0.791575,0.822604,0.856512,0.955252,0.976810,0.988155 /)
XCGA_LKT(54,12,1:6)=(/ 0.866427,0.844550,0.817500,0.770153,0.754243,0.746253 /)
XEXT_COEFF_550_LKT(54,12)=215.260000 !rg=0.724875 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(54,13,1:6)=(/ 169.330000,172.670000,177.020000,185.780000,200.120000,210.770000 /)
XPIZA_LKT(54,13,1:6)=(/ 0.784130,0.812151,0.843647,0.946898,0.972268,0.985952 /)
XCGA_LKT(54,13,1:6)=(/ 0.869940,0.852413,0.827747,0.781923,0.762987,0.748173 /)
XEXT_COEFF_550_LKT(54,13)=177.170000 !rg=0.724875 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(54,14,1:6)=(/ 140.610000,143.200000,146.660000,152.810000,162.950000,174.110000 /)
XPIZA_LKT(54,14,1:6)=(/ 0.776626,0.802127,0.831217,0.938594,0.966943,0.983163 /)
XCGA_LKT(54,14,1:6)=(/ 0.873850,0.857507,0.837360,0.788367,0.770957,0.749487 /)
XEXT_COEFF_550_LKT(54,14)=145.960000 !rg=0.724875 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(54,15,1:6)=(/ 116.300000,118.600000,121.330000,125.410000,133.990000,143.560000 /)
XPIZA_LKT(54,15,1:6)=(/ 0.768011,0.793582,0.819699,0.930815,0.960034,0.979492 /)
XCGA_LKT(54,15,1:6)=(/ 0.876577,0.863670,0.844217,0.801163,0.774713,0.750430 /)
XEXT_COEFF_550_LKT(54,15)=120.960000 !rg=0.724875 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(54,16,1:6)=(/ 96.372000,98.005000,99.752000,103.030000,109.590000,117.620000 /)
XPIZA_LKT(54,16,1:6)=(/ 0.759354,0.785002,0.808778,0.924444,0.952433,0.975076 /)
XCGA_LKT(54,16,1:6)=(/ 0.879560,0.868513,0.852987,0.811080,0.785313,0.756930 /)
XEXT_COEFF_550_LKT(54,16)=99.913000 !rg=0.724875 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(54,17,1:6)=(/ 80.131000,81.039000,82.528000,84.797000,88.699000,95.315000 /)
XPIZA_LKT(54,17,1:6)=(/ 0.749554,0.776011,0.800303,0.917896,0.945189,0.971300 /)
XCGA_LKT(54,17,1:6)=(/ 0.883093,0.872070,0.858850,0.815860,0.798107,0.765857 /)
XEXT_COEFF_550_LKT(54,17)=82.461000 !rg=0.724875 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(54,18,1:6)=(/ 66.260000,67.192000,68.336000,69.961000,73.433000,78.503000 /)
XPIZA_LKT(54,18,1:6)=(/ 0.738392,0.767039,0.790748,0.911426,0.936968,0.965455 /)
XCGA_LKT(54,18,1:6)=(/ 0.885477,0.875977,0.863483,0.824287,0.802763,0.768780 /)
XEXT_COEFF_550_LKT(54,18)=68.210000 !rg=0.724875 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(54,19,1:6)=(/ 54.916000,55.652000,56.328000,57.670000,60.180000,64.178000 /)
XPIZA_LKT(54,19,1:6)=(/ 0.726468,0.757465,0.781145,0.906402,0.930000,0.959225 /)
XCGA_LKT(54,19,1:6)=(/ 0.888337,0.879490,0.868993,0.831910,0.812510,0.778780 /)
XEXT_COEFF_550_LKT(54,19)=56.400000 !rg=0.724875 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(54,20,1:6)=(/ 45.638000,46.013000,46.630000,47.509000,49.025000,51.992000 /)
XPIZA_LKT(54,20,1:6)=(/ 0.714207,0.746182,0.772515,0.900771,0.924013,0.952367 /)
XCGA_LKT(54,20,1:6)=(/ 0.891900,0.882563,0.873103,0.835320,0.822767,0.789913 /)
XEXT_COEFF_550_LKT(54,20)=46.612000 !rg=0.724875 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(55,1,1:6)=(/ 809.490000,931.860000,710.300000,1435.100000,1051.700000,217.210000 /)
XPIZA_LKT(55,1,1:6)=(/ 0.886720,0.926547,0.932052,0.991553,0.994636,0.989789 /)
XCGA_LKT(55,1,1:6)=(/ 0.795933,0.805500,0.573180,0.757893,0.769503,0.515787 /)
XEXT_COEFF_550_LKT(55,1)=708.270000 !rg=0.785298 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(55,2,1:6)=(/ 819.120000,871.300000,783.780000,1280.900000,1050.200000,247.220000 /)
XPIZA_LKT(55,2,1:6)=(/ 0.886666,0.919742,0.937875,0.990435,0.994557,0.990639 /)
XCGA_LKT(55,2,1:6)=(/ 0.810413,0.783223,0.639860,0.741723,0.773643,0.564690 /)
XEXT_COEFF_550_LKT(55,2)=782.450000 !rg=0.785298 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(55,3,1:6)=(/ 757.430000,794.680000,785.340000,1078.400000,1028.400000,292.810000 /)
XPIZA_LKT(55,3,1:6)=(/ 0.878786,0.916628,0.938806,0.988587,0.994378,0.991797 /)
XCGA_LKT(55,3,1:6)=(/ 0.816573,0.789177,0.687767,0.721330,0.778713,0.612600 /)
XEXT_COEFF_550_LKT(55,3)=789.070000 !rg=0.785298 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(55,4,1:6)=(/ 681.490000,709.150000,730.850000,904.750000,961.120000,343.090000 /)
XPIZA_LKT(55,4,1:6)=(/ 0.870215,0.910941,0.935407,0.986398,0.993903,0.992733 /)
XCGA_LKT(55,4,1:6)=(/ 0.821677,0.789510,0.713110,0.710287,0.778713,0.660463 /)
XEXT_COEFF_550_LKT(55,4)=731.150000 !rg=0.785298 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(55,5,1:6)=(/ 601.680000,626.770000,653.860000,761.270000,853.950000,388.820000 /)
XPIZA_LKT(55,5,1:6)=(/ 0.860287,0.901987,0.930378,0.983772,0.993034,0.993342 /)
XCGA_LKT(55,5,1:6)=(/ 0.828317,0.785460,0.734930,0.708180,0.772557,0.700270 /)
XEXT_COEFF_550_LKT(55,5)=645.410000 !rg=0.785298 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(55,6,1:6)=(/ 521.460000,537.700000,565.730000,639.750000,731.120000,418.960000 /)
XPIZA_LKT(55,6,1:6)=(/ 0.848325,0.894075,0.923850,0.981144,0.991781,0.993600 /)
XCGA_LKT(55,6,1:6)=(/ 0.835817,0.799147,0.749150,0.714663,0.763793,0.728010 /)
XEXT_COEFF_550_LKT(55,6)=567.890000 !rg=0.785298 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(55,7,1:6)=(/ 447.380000,460.750000,484.500000,533.630000,611.720000,425.890000 /)
XPIZA_LKT(55,7,1:6)=(/ 0.833999,0.881662,0.915271,0.978141,0.990067,0.993504 /)
XCGA_LKT(55,7,1:6)=(/ 0.840983,0.807740,0.760493,0.724403,0.754983,0.744260 /)
XEXT_COEFF_550_LKT(55,7)=484.070000 !rg=0.785298 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(55,8,1:6)=(/ 380.220000,391.170000,408.930000,444.050000,504.980000,407.420000 /)
XPIZA_LKT(55,8,1:6)=(/ 0.822755,0.868409,0.904946,0.974102,0.987955,0.993026 /)
XCGA_LKT(55,8,1:6)=(/ 0.848883,0.815273,0.776880,0.731873,0.748690,0.750317 /)
XEXT_COEFF_550_LKT(55,8)=410.430000 !rg=0.785298 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(55,9,1:6)=(/ 321.030000,330.300000,344.410000,369.310000,415.150000,371.120000 /)
XPIZA_LKT(55,9,1:6)=(/ 0.811968,0.853828,0.893093,0.969369,0.985266,0.992148 /)
XCGA_LKT(55,9,1:6)=(/ 0.854063,0.822763,0.789910,0.739010,0.746267,0.750070 /)
XEXT_COEFF_550_LKT(55,9)=340.790000 !rg=0.785298 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(55,10,1:6)=(/ 269.600000,277.120000,286.060000,306.300000,338.580000,327.230000 /)
XPIZA_LKT(55,10,1:6)=(/ 0.802590,0.839978,0.879975,0.963841,0.982578,0.990970 /)
XCGA_LKT(55,10,1:6)=(/ 0.859267,0.834033,0.800507,0.753717,0.748283,0.749250 /)
XEXT_COEFF_550_LKT(55,10)=285.880000 !rg=0.785298 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(55,11,1:6)=(/ 224.930000,231.060000,237.870000,251.590000,275.840000,280.500000 /)
XPIZA_LKT(55,11,1:6)=(/ 0.795135,0.827571,0.864993,0.957441,0.978478,0.989373 /)
XCGA_LKT(55,11,1:6)=(/ 0.864573,0.841117,0.813640,0.762297,0.751297,0.747920 /)
XEXT_COEFF_550_LKT(55,11)=238.520000 !rg=0.785298 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(55,12,1:6)=(/ 187.410000,191.430000,196.790000,207.160000,224.380000,234.950000 /)
XPIZA_LKT(55,12,1:6)=(/ 0.786760,0.816520,0.850776,0.951034,0.975086,0.987339 /)
XCGA_LKT(55,12,1:6)=(/ 0.868533,0.849510,0.822753,0.776027,0.758843,0.747440 /)
XEXT_COEFF_550_LKT(55,12)=196.590000 !rg=0.785298 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(55,13,1:6)=(/ 155.900000,158.790000,162.960000,170.140000,182.440000,195.060000 /)
XPIZA_LKT(55,13,1:6)=(/ 0.779815,0.806230,0.837295,0.942696,0.969797,0.984843 /)
XCGA_LKT(55,13,1:6)=(/ 0.872653,0.855077,0.833410,0.782990,0.766110,0.748540 /)
XEXT_COEFF_550_LKT(55,13)=162.270000 !rg=0.785298 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(55,14,1:6)=(/ 129.150000,131.580000,134.800000,140.030000,149.690000,160.650000 /)
XPIZA_LKT(55,14,1:6)=(/ 0.772718,0.797498,0.824961,0.934538,0.963263,0.981288 /)
XCGA_LKT(55,14,1:6)=(/ 0.875357,0.861157,0.840240,0.795687,0.770107,0.748353 /)
XEXT_COEFF_550_LKT(55,14)=134.600000 !rg=0.785298 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(55,15,1:6)=(/ 107.220000,109.110000,111.360000,114.840000,122.570000,132.270000 /)
XPIZA_LKT(55,15,1:6)=(/ 0.764131,0.788996,0.814481,0.927910,0.957213,0.977406 /)
XCGA_LKT(55,15,1:6)=(/ 0.878797,0.866520,0.849220,0.806520,0.780487,0.753240 /)
XEXT_COEFF_550_LKT(55,15)=111.150000 !rg=0.785298 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(55,16,1:6)=(/ 88.832000,90.278000,91.840000,94.690000,99.286000,107.070000 /)
XPIZA_LKT(55,16,1:6)=(/ 0.755150,0.780608,0.804715,0.920099,0.949603,0.973674 /)
XCGA_LKT(55,16,1:6)=(/ 0.881607,0.870400,0.856223,0.812130,0.792183,0.761200 /)
XEXT_COEFF_550_LKT(55,16)=91.913000 !rg=0.785298 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(55,17,1:6)=(/ 73.686000,74.744000,75.868000,78.104000,81.450000,87.600000 /)
XPIZA_LKT(55,17,1:6)=(/ 0.744430,0.772056,0.795273,0.914119,0.941405,0.968974 /)
XCGA_LKT(55,17,1:6)=(/ 0.884307,0.874360,0.860973,0.820437,0.799733,0.766047 /)
XEXT_COEFF_550_LKT(55,17)=75.907000 !rg=0.785298 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(55,18,1:6)=(/ 61.159000,61.920000,62.803000,64.248000,67.171000,71.989000 /)
XPIZA_LKT(55,18,1:6)=(/ 0.733356,0.762494,0.786431,0.908747,0.933959,0.962940 /)
XCGA_LKT(55,18,1:6)=(/ 0.887323,0.877950,0.866453,0.828030,0.808017,0.773687 /)
XEXT_COEFF_550_LKT(55,18)=62.645000 !rg=0.785298 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(55,19,1:6)=(/ 50.678000,51.246000,51.934000,52.963000,54.805000,58.293000 /)
XPIZA_LKT(55,19,1:6)=(/ 0.721236,0.752399,0.777600,0.902995,0.927273,0.956439 /)
XCGA_LKT(55,19,1:6)=(/ 0.890217,0.881037,0.871077,0.832523,0.818143,0.784343 /)
XEXT_COEFF_550_LKT(55,19)=51.989000 !rg=0.785298 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(55,20,1:6)=(/ 42.020000,42.475000,42.892000,43.799000,45.161000,47.804000 /)
XPIZA_LKT(55,20,1:6)=(/ 0.708263,0.741513,0.767416,0.897646,0.920724,0.948292 /)
XCGA_LKT(55,20,1:6)=(/ 0.893177,0.884337,0.874593,0.838437,0.824443,0.790830 /)
XEXT_COEFF_550_LKT(55,20)=42.942000 !rg=0.785298 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(56,1,1:6)=(/ 754.920000,730.860000,729.130000,1211.400000,1070.800000,252.030000 /)
XPIZA_LKT(56,1,1:6)=(/ 0.879025,0.912839,0.933514,0.989775,0.994660,0.990627 /)
XCGA_LKT(56,1,1:6)=(/ 0.806477,0.772977,0.640473,0.732203,0.776533,0.586567 /)
XEXT_COEFF_550_LKT(56,1)=732.060000 !rg=0.850757 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(56,2,1:6)=(/ 749.450000,758.750000,764.080000,1072.600000,1068.300000,284.950000 /)
XPIZA_LKT(56,2,1:6)=(/ 0.878888,0.916407,0.936951,0.988530,0.994583,0.991501 /)
XCGA_LKT(56,2,1:6)=(/ 0.818133,0.795450,0.685207,0.715877,0.783873,0.608197 /)
XEXT_COEFF_550_LKT(56,2)=775.530000 !rg=0.850757 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(56,3,1:6)=(/ 692.510000,721.910000,734.200000,916.380000,1022.800000,328.060000 /)
XPIZA_LKT(56,3,1:6)=(/ 0.873115,0.912233,0.936100,0.986488,0.994282,0.992441 /)
XCGA_LKT(56,3,1:6)=(/ 0.822413,0.795293,0.715403,0.698407,0.785387,0.642537 /)
XEXT_COEFF_550_LKT(56,3)=737.110000 !rg=0.850757 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(56,4,1:6)=(/ 623.340000,649.180000,668.940000,788.150000,930.380000,374.150000 /)
XPIZA_LKT(56,4,1:6)=(/ 0.864127,0.905833,0.932280,0.984298,0.993635,0.993161 /)
XCGA_LKT(56,4,1:6)=(/ 0.827787,0.796870,0.734180,0.696430,0.780907,0.683093 /)
XEXT_COEFF_550_LKT(56,4)=670.830000 !rg=0.850757 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(56,5,1:6)=(/ 551.420000,571.380000,596.500000,675.700000,806.670000,412.350000 /)
XPIZA_LKT(56,5,1:6)=(/ 0.851461,0.896400,0.925050,0.981956,0.992579,0.993585 /)
XCGA_LKT(56,5,1:6)=(/ 0.832777,0.800243,0.741673,0.705687,0.771600,0.716390 /)
XEXT_COEFF_550_LKT(56,5)=594.970000 !rg=0.850757 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(56,6,1:6)=(/ 480.250000,496.630000,513.690000,576.380000,680.200000,431.620000 /)
XPIZA_LKT(56,6,1:6)=(/ 0.839755,0.886648,0.920140,0.978744,0.991100,0.993676 /)
XCGA_LKT(56,6,1:6)=(/ 0.839447,0.806633,0.761327,0.712860,0.761113,0.738427 /)
XEXT_COEFF_550_LKT(56,6)=517.480000 !rg=0.850757 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(56,7,1:6)=(/ 410.210000,423.490000,439.990000,482.230000,563.570000,426.330000 /)
XPIZA_LKT(56,7,1:6)=(/ 0.827341,0.874487,0.910730,0.975759,0.989045,0.993410 /)
XCGA_LKT(56,7,1:6)=(/ 0.846130,0.814693,0.774660,0.725277,0.751963,0.749813 /)
XEXT_COEFF_550_LKT(56,7)=442.460000 !rg=0.850757 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(56,8,1:6)=(/ 350.280000,360.440000,372.590000,405.270000,462.430000,397.820000 /)
XPIZA_LKT(56,8,1:6)=(/ 0.815607,0.860536,0.899973,0.971169,0.986918,0.992772 /)
XCGA_LKT(56,8,1:6)=(/ 0.851913,0.822650,0.784447,0.736110,0.747883,0.752690 /)
XEXT_COEFF_550_LKT(56,8)=373.500000 !rg=0.850757 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(56,9,1:6)=(/ 294.440000,302.860000,314.450000,334.450000,378.550000,355.750000 /)
XPIZA_LKT(56,9,1:6)=(/ 0.805672,0.846038,0.886232,0.967079,0.983904,0.991753 /)
XCGA_LKT(56,9,1:6)=(/ 0.857297,0.830743,0.793103,0.749707,0.744273,0.751310 /)
XEXT_COEFF_550_LKT(56,9)=314.070000 !rg=0.850757 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(56,10,1:6)=(/ 247.980000,254.450000,263.210000,278.590000,309.050000,309.870000 /)
XPIZA_LKT(56,10,1:6)=(/ 0.797602,0.833061,0.871756,0.961485,0.980823,0.990345 /)
XCGA_LKT(56,10,1:6)=(/ 0.862980,0.837387,0.808313,0.754553,0.748223,0.749073 /)
XEXT_COEFF_550_LKT(56,10)=260.560000 !rg=0.850757 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(56,11,1:6)=(/ 207.010000,212.320000,218.470000,230.500000,251.450000,261.850000 /)
XPIZA_LKT(56,11,1:6)=(/ 0.790183,0.821034,0.857355,0.954455,0.976905,0.988640 /)
XCGA_LKT(56,11,1:6)=(/ 0.866813,0.846630,0.817410,0.769777,0.752913,0.747947 /)
XEXT_COEFF_550_LKT(56,11)=217.900000 !rg=0.850757 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(56,12,1:6)=(/ 172.570000,176.280000,180.930000,189.390000,204.440000,218.100000 /)
XPIZA_LKT(56,12,1:6)=(/ 0.783302,0.810670,0.843693,0.947374,0.972730,0.986305 /)
XCGA_LKT(56,12,1:6)=(/ 0.871123,0.852373,0.829057,0.777283,0.760983,0.747287 /)
XEXT_COEFF_550_LKT(56,12)=179.630000 !rg=0.850757 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(56,13,1:6)=(/ 143.210000,146.140000,149.630000,156.010000,167.150000,180.080000 /)
XPIZA_LKT(56,13,1:6)=(/ 0.776079,0.801763,0.830599,0.938446,0.966885,0.983299 /)
XCGA_LKT(56,13,1:6)=(/ 0.874280,0.858983,0.836390,0.790410,0.765740,0.747320 /)
XEXT_COEFF_550_LKT(56,13)=149.410000 !rg=0.850757 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(56,14,1:6)=(/ 118.880000,121.200000,123.770000,127.830000,137.120000,148.320000 /)
XPIZA_LKT(56,14,1:6)=(/ 0.767772,0.792215,0.819468,0.931032,0.960593,0.979693 /)
XCGA_LKT(56,14,1:6)=(/ 0.877540,0.864183,0.845550,0.800577,0.775680,0.750887 /)
XEXT_COEFF_550_LKT(56,14)=123.840000 !rg=0.850757 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(56,15,1:6)=(/ 98.910000,100.440000,102.350000,105.570000,111.300000,120.450000 /)
XPIZA_LKT(56,15,1:6)=(/ 0.760019,0.784823,0.809242,0.924543,0.954641,0.976462 /)
XCGA_LKT(56,15,1:6)=(/ 0.880623,0.868423,0.853147,0.809187,0.786330,0.756750 /)
XEXT_COEFF_550_LKT(56,15)=102.270000 !rg=0.850757 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(56,16,1:6)=(/ 81.829000,83.162000,84.602000,86.907000,90.742000,98.074000 /)
XPIZA_LKT(56,16,1:6)=(/ 0.750228,0.776314,0.799917,0.917824,0.946200,0.972213 /)
XCGA_LKT(56,16,1:6)=(/ 0.882683,0.873063,0.858647,0.817503,0.796230,0.763353 /)
XEXT_COEFF_550_LKT(56,16)=84.404000 !rg=0.850757 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(56,17,1:6)=(/ 67.955000,68.662000,69.914000,71.449000,74.815000,80.573000 /)
XPIZA_LKT(56,17,1:6)=(/ 0.739541,0.767286,0.791032,0.911534,0.938155,0.966206 /)
XCGA_LKT(56,17,1:6)=(/ 0.885993,0.875870,0.864217,0.823697,0.804637,0.769833 /)
XEXT_COEFF_550_LKT(56,17)=69.939000 !rg=0.850757 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(56,18,1:6)=(/ 56.420000,56.974000,57.834000,59.138000,61.299000,65.431000 /)
XPIZA_LKT(56,18,1:6)=(/ 0.728147,0.758002,0.781972,0.906075,0.931370,0.960941 /)
XCGA_LKT(56,18,1:6)=(/ 0.888903,0.879080,0.869073,0.830220,0.813397,0.778750 /)
XEXT_COEFF_550_LKT(56,18)=57.756000 !rg=0.850757 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(56,19,1:6)=(/ 46.692000,47.239000,47.811000,48.768000,50.143000,53.327000 /)
XPIZA_LKT(56,19,1:6)=(/ 0.715348,0.747494,0.772826,0.900901,0.924483,0.953729 /)
XCGA_LKT(56,19,1:6)=(/ 0.891470,0.883000,0.872847,0.836293,0.821443,0.787643 /)
XEXT_COEFF_550_LKT(56,19)=47.772000 !rg=0.850757 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(56,20,1:6)=(/ 38.775000,39.066000,39.587000,40.173000,41.430000,43.896000 /)
XPIZA_LKT(56,20,1:6)=(/ 0.702564,0.735861,0.763169,0.894986,0.918181,0.944881 /)
XCGA_LKT(56,20,1:6)=(/ 0.894810,0.885623,0.876847,0.840497,0.828297,0.795647 /)
XEXT_COEFF_550_LKT(56,20)=39.569000 !rg=0.850757 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(57,1,1:6)=(/ 751.890000,676.490000,783.480000,973.640000,1101.100000,296.260000 /)
XPIZA_LKT(57,1,1:6)=(/ 0.880264,0.895121,0.938389,0.987331,0.994674,0.991523 /)
XCGA_LKT(57,1,1:6)=(/ 0.830903,0.777947,0.716433,0.695690,0.789690,0.634220 /)
XEXT_COEFF_550_LKT(57,1)=787.620000 !rg=0.921673 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(57,2,1:6)=(/ 690.750000,701.620000,754.170000,878.910000,1069.000000,324.580000 /)
XPIZA_LKT(57,2,1:6)=(/ 0.871952,0.908644,0.937199,0.985844,0.994524,0.992273 /)
XCGA_LKT(57,2,1:6)=(/ 0.824470,0.798797,0.729517,0.679870,0.791667,0.638293 /)
XEXT_COEFF_550_LKT(57,2)=762.360000 !rg=0.921673 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(57,3,1:6)=(/ 639.270000,656.460000,686.850000,775.560000,999.260000,362.960000 /)
XPIZA_LKT(57,3,1:6)=(/ 0.865038,0.906435,0.932924,0.984061,0.994086,0.992979 /)
XCGA_LKT(57,3,1:6)=(/ 0.828500,0.796310,0.738370,0.679237,0.789203,0.667997 /)
XEXT_COEFF_550_LKT(57,3)=694.470000 !rg=0.921673 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(57,4,1:6)=(/ 575.170000,593.410000,615.510000,689.000000,885.140000,403.290000 /)
XPIZA_LKT(57,4,1:6)=(/ 0.855666,0.899578,0.928107,0.982083,0.993252,0.993507 /)
XCGA_LKT(57,4,1:6)=(/ 0.833777,0.799263,0.748660,0.689077,0.780553,0.703037 /)
XEXT_COEFF_550_LKT(57,4)=622.130000 !rg=0.921673 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(57,5,1:6)=(/ 504.820000,523.340000,543.400000,601.180000,752.500000,432.200000 /)
XPIZA_LKT(57,5,1:6)=(/ 0.843922,0.890670,0.921338,0.979762,0.991967,0.993759 /)
XCGA_LKT(57,5,1:6)=(/ 0.839083,0.807347,0.760723,0.702347,0.768420,0.730347 /)
XEXT_COEFF_550_LKT(57,5)=544.740000 !rg=0.921673 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(57,6,1:6)=(/ 441.220000,456.310000,473.380000,516.870000,625.550000,439.430000 /)
XPIZA_LKT(57,6,1:6)=(/ 0.830845,0.878486,0.913644,0.976485,0.990276,0.993687 /)
XCGA_LKT(57,6,1:6)=(/ 0.845327,0.808723,0.770587,0.713633,0.756733,0.746977 /)
XEXT_COEFF_550_LKT(57,6)=470.510000 !rg=0.921673 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(57,7,1:6)=(/ 377.440000,389.480000,403.170000,436.380000,513.620000,421.900000 /)
XPIZA_LKT(57,7,1:6)=(/ 0.819571,0.866255,0.904688,0.972878,0.988126,0.993255 /)
XCGA_LKT(57,7,1:6)=(/ 0.851490,0.818273,0.782473,0.726920,0.748740,0.753937 /)
XEXT_COEFF_550_LKT(57,7)=405.520000 !rg=0.921673 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(57,8,1:6)=(/ 321.800000,330.940000,343.190000,367.630000,421.400000,385.070000 /)
XPIZA_LKT(57,8,1:6)=(/ 0.809017,0.852407,0.892716,0.968395,0.985635,0.992442 /)
XCGA_LKT(57,8,1:6)=(/ 0.856917,0.825707,0.793103,0.736360,0.745333,0.753880 /)
XEXT_COEFF_550_LKT(57,8)=339.760000 !rg=0.921673 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(57,9,1:6)=(/ 270.380000,278.050000,287.310000,304.030000,345.810000,339.180000 /)
XPIZA_LKT(57,9,1:6)=(/ 0.800079,0.838972,0.878721,0.964549,0.981852,0.991282 /)
XCGA_LKT(57,9,1:6)=(/ 0.861823,0.836623,0.803967,0.754890,0.744167,0.751927 /)
XEXT_COEFF_550_LKT(57,9)=287.960000 !rg=0.921673 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(57,10,1:6)=(/ 227.470000,233.510000,241.050000,253.890000,281.950000,290.850000 /)
XPIZA_LKT(57,10,1:6)=(/ 0.793184,0.826344,0.864370,0.957992,0.978551,0.989720 /)
XCGA_LKT(57,10,1:6)=(/ 0.865397,0.843463,0.811477,0.764707,0.746760,0.749373 /)
XEXT_COEFF_550_LKT(57,10)=240.300000 !rg=0.921673 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(57,11,1:6)=(/ 190.360000,194.700000,200.760000,209.550000,229.900000,243.700000 /)
XPIZA_LKT(57,11,1:6)=(/ 0.786089,0.815565,0.849979,0.951883,0.974709,0.987558 /)
XCGA_LKT(57,11,1:6)=(/ 0.869517,0.849920,0.824690,0.772837,0.755080,0.746953 /)
XEXT_COEFF_550_LKT(57,11)=199.430000 !rg=0.921673 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(57,12,1:6)=(/ 158.500000,161.830000,166.220000,173.290000,187.550000,202.000000 /)
XPIZA_LKT(57,12,1:6)=(/ 0.779669,0.805689,0.836726,0.942954,0.969331,0.985026 /)
XCGA_LKT(57,12,1:6)=(/ 0.872820,0.856647,0.832237,0.785723,0.760013,0.746633 /)
XEXT_COEFF_550_LKT(57,12)=165.780000 !rg=0.921673 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(57,13,1:6)=(/ 131.810000,134.210000,137.460000,142.210000,153.280000,166.030000 /)
XPIZA_LKT(57,13,1:6)=(/ 0.771424,0.796169,0.824685,0.934798,0.964020,0.981775 /)
XCGA_LKT(57,13,1:6)=(/ 0.876523,0.861850,0.842097,0.795330,0.770707,0.748607 /)
XEXT_COEFF_550_LKT(57,13)=137.620000 !rg=0.921673 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(57,14,1:6)=(/ 109.700000,111.210000,113.770000,117.120000,124.560000,135.320000 /)
XPIZA_LKT(57,14,1:6)=(/ 0.764517,0.788883,0.813193,0.928444,0.958638,0.978981 /)
XCGA_LKT(57,14,1:6)=(/ 0.879070,0.866780,0.849540,0.805103,0.780343,0.753323 /)
XEXT_COEFF_550_LKT(57,14)=113.770000 !rg=0.921673 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(57,15,1:6)=(/ 90.936000,92.384000,94.154000,96.967000,101.510000,110.310000 /)
XPIZA_LKT(57,15,1:6)=(/ 0.755195,0.780329,0.804617,0.920958,0.951652,0.974989 /)
XCGA_LKT(57,15,1:6)=(/ 0.881687,0.871113,0.855927,0.813627,0.792413,0.760523 /)
XEXT_COEFF_550_LKT(57,15)=93.852000 !rg=0.921673 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(57,16,1:6)=(/ 75.326000,76.354000,77.837000,79.539000,83.609000,90.214000 /)
XPIZA_LKT(57,16,1:6)=(/ 0.745544,0.771973,0.795153,0.914895,0.942430,0.969620 /)
XCGA_LKT(57,16,1:6)=(/ 0.884253,0.874660,0.862310,0.820480,0.800623,0.766277 /)
XEXT_COEFF_550_LKT(57,16)=77.670000 !rg=0.921673 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(57,17,1:6)=(/ 62.572000,63.272000,64.165000,65.463000,68.477000,73.475000 /)
XPIZA_LKT(57,17,1:6)=(/ 0.734440,0.763063,0.786338,0.909058,0.935494,0.964655 /)
XCGA_LKT(57,17,1:6)=(/ 0.887313,0.878247,0.866313,0.827993,0.808207,0.773070 /)
XEXT_COEFF_550_LKT(57,17)=64.313000 !rg=0.921673 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(57,18,1:6)=(/ 51.927000,52.477000,53.181000,54.333000,56.067000,59.757000 /)
XPIZA_LKT(57,18,1:6)=(/ 0.722423,0.753141,0.777881,0.903508,0.928869,0.958569 /)
XCGA_LKT(57,18,1:6)=(/ 0.889997,0.881240,0.870577,0.833510,0.818887,0.784187 /)
XEXT_COEFF_550_LKT(57,18)=53.111000 !rg=0.921673 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(57,19,1:6)=(/ 43.026000,43.445000,44.049000,44.694000,46.307000,49.136000 /)
XPIZA_LKT(57,19,1:6)=(/ 0.709887,0.742240,0.768225,0.898042,0.921638,0.949736 /)
XCGA_LKT(57,19,1:6)=(/ 0.892933,0.884200,0.875343,0.838343,0.825603,0.791920 /)
XEXT_COEFF_550_LKT(57,19)=43.974000 !rg=0.921673 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(57,20,1:6)=(/ 35.721000,36.012000,36.385000,36.922000,38.041000,40.114000 /)
XPIZA_LKT(57,20,1:6)=(/ 0.696585,0.730756,0.758045,0.892460,0.915622,0.942242 /)
XCGA_LKT(57,20,1:6)=(/ 0.896223,0.887393,0.878357,0.843283,0.831183,0.799377 /)
XEXT_COEFF_550_LKT(57,20)=36.434000 !rg=0.921673 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(58,1,1:6)=(/ 665.790000,675.900000,813.840000,757.730000,1092.300000,345.170000 /)
XPIZA_LKT(58,1,1:6)=(/ 0.867802,0.904200,0.941608,0.983725,0.994660,0.992457 /)
XCGA_LKT(58,1,1:6)=(/ 0.830813,0.787103,0.773773,0.643930,0.799963,0.653237 /)
XEXT_COEFF_550_LKT(58,1)=827.540000 !rg=0.9985 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(58,2,1:6)=(/ 638.410000,653.090000,728.790000,708.330000,1050.000000,362.990000 /)
XPIZA_LKT(58,2,1:6)=(/ 0.862864,0.902251,0.934850,0.982603,0.994366,0.992927 /)
XCGA_LKT(58,2,1:6)=(/ 0.831550,0.792557,0.760900,0.645350,0.796650,0.660037 /)
XEXT_COEFF_550_LKT(58,2)=726.120000 !rg=0.9985 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(58,3,1:6)=(/ 588.320000,607.680000,637.710000,668.140000,957.570000,396.360000 /)
XPIZA_LKT(58,3,1:6)=(/ 0.856362,0.898830,0.929466,0.981200,0.993771,0.993420 /)
XCGA_LKT(58,3,1:6)=(/ 0.833877,0.805570,0.754547,0.660340,0.790043,0.690457 /)
XEXT_COEFF_550_LKT(58,3)=639.010000 !rg=0.9985 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(58,4,1:6)=(/ 528.800000,547.950000,565.460000,612.260000,827.780000,429.500000 /)
XPIZA_LKT(58,4,1:6)=(/ 0.846840,0.891907,0.923610,0.979428,0.992728,0.993776 /)
XCGA_LKT(58,4,1:6)=(/ 0.838580,0.808697,0.758940,0.683900,0.777790,0.720593 /)
XEXT_COEFF_550_LKT(58,4)=566.900000 !rg=0.9985 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(58,5,1:6)=(/ 465.460000,479.950000,497.510000,538.340000,691.490000,447.540000 /)
XPIZA_LKT(58,5,1:6)=(/ 0.834977,0.883804,0.916775,0.977546,0.991221,0.993866 /)
XCGA_LKT(58,5,1:6)=(/ 0.844607,0.810843,0.771073,0.703837,0.763513,0.742240 /)
XEXT_COEFF_550_LKT(58,5)=502.560000 !rg=0.9985 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(58,6,1:6)=(/ 404.270000,418.500000,433.470000,469.520000,570.430000,442.000000 /)
XPIZA_LKT(58,6,1:6)=(/ 0.822852,0.871291,0.906866,0.973201,0.989312,0.993634 /)
XCGA_LKT(58,6,1:6)=(/ 0.849893,0.818293,0.774263,0.719163,0.752383,0.753880 /)
XEXT_COEFF_550_LKT(58,6)=430.380000 !rg=0.9985 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(58,7,1:6)=(/ 347.070000,358.190000,369.430000,398.520000,467.350000,413.240000 /)
XPIZA_LKT(58,7,1:6)=(/ 0.811813,0.858276,0.897741,0.969816,0.986931,0.993035 /)
XCGA_LKT(58,7,1:6)=(/ 0.855263,0.826253,0.787043,0.733457,0.745850,0.757063 /)
XEXT_COEFF_550_LKT(58,7)=369.280000 !rg=0.9985 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(58,8,1:6)=(/ 294.870000,303.770000,314.370000,335.020000,382.980000,369.620000 /)
XPIZA_LKT(58,8,1:6)=(/ 0.802931,0.844559,0.885016,0.965140,0.984038,0.992060 /)
XCGA_LKT(58,8,1:6)=(/ 0.860153,0.833053,0.795873,0.746137,0.742783,0.755090 /)
XEXT_COEFF_550_LKT(58,8)=312.590000 !rg=0.9985 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(58,9,1:6)=(/ 249.080000,255.760000,263.190000,277.400000,312.210000,320.520000 /)
XPIZA_LKT(58,9,1:6)=(/ 0.795334,0.832150,0.871780,0.961634,0.980885,0.990718 /)
XCGA_LKT(58,9,1:6)=(/ 0.865150,0.840427,0.811003,0.757800,0.745903,0.752003 /)
XEXT_COEFF_550_LKT(58,9)=264.370000 !rg=0.9985 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(58,10,1:6)=(/ 209.260000,214.600000,220.920000,230.790000,257.420000,271.980000 /)
XPIZA_LKT(58,10,1:6)=(/ 0.788141,0.819364,0.856795,0.955326,0.976372,0.988708 /)
XCGA_LKT(58,10,1:6)=(/ 0.868890,0.848117,0.820157,0.771487,0.749417,0.748350 /)
XEXT_COEFF_550_LKT(58,10)=220.980000 !rg=0.9985 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(58,11,1:6)=(/ 175.230000,178.890000,183.680000,190.690000,209.690000,225.770000 /)
XPIZA_LKT(58,11,1:6)=(/ 0.781463,0.808844,0.843595,0.948334,0.972422,0.986494 /)
XCGA_LKT(58,11,1:6)=(/ 0.871817,0.854963,0.828250,0.782837,0.754327,0.745997 /)
XEXT_COEFF_550_LKT(58,11)=183.490000 !rg=0.9985 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(58,12,1:6)=(/ 145.900000,149.000000,152.480000,158.040000,171.530000,186.650000 /)
XPIZA_LKT(58,12,1:6)=(/ 0.774872,0.799792,0.830528,0.939201,0.966965,0.983429 /)
XCGA_LKT(58,12,1:6)=(/ 0.875277,0.860070,0.838503,0.791617,0.765187,0.747327 /)
XEXT_COEFF_550_LKT(58,12)=152.590000 !rg=0.9985 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(58,13,1:6)=(/ 121.480000,123.230000,126.100000,129.960000,139.180000,151.910000 /)
XPIZA_LKT(58,13,1:6)=(/ 0.768306,0.792107,0.818343,0.932398,0.962382,0.980979 /)
XCGA_LKT(58,13,1:6)=(/ 0.877987,0.865240,0.846013,0.801000,0.774907,0.750337 /)
XEXT_COEFF_550_LKT(58,13)=126.310000 !rg=0.9985 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(58,14,1:6)=(/ 100.990000,102.490000,104.300000,107.580000,113.790000,123.820000 /)
XPIZA_LKT(58,14,1:6)=(/ 0.760335,0.784767,0.809294,0.924834,0.955422,0.977448 /)
XCGA_LKT(58,14,1:6)=(/ 0.880457,0.869627,0.853353,0.810330,0.787053,0.757270 /)
XEXT_COEFF_550_LKT(58,14)=104.220000 !rg=0.9985 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(58,15,1:6)=(/ 83.904000,85.026000,86.518000,88.915000,93.280000,101.110000 /)
XPIZA_LKT(58,15,1:6)=(/ 0.751171,0.776188,0.799331,0.918175,0.946861,0.972592 /)
XCGA_LKT(58,15,1:6)=(/ 0.883273,0.872983,0.859593,0.815943,0.796127,0.762350 /)
XEXT_COEFF_550_LKT(58,15)=86.072000 !rg=0.9985 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(58,16,1:6)=(/ 69.492000,70.294000,71.409000,72.870000,76.804000,82.782000 /)
XPIZA_LKT(58,16,1:6)=(/ 0.740368,0.767141,0.790765,0.912327,0.938943,0.966976 /)
XCGA_LKT(58,16,1:6)=(/ 0.885813,0.876977,0.864477,0.825663,0.801447,0.766070 /)
XEXT_COEFF_550_LKT(58,16)=71.385000 !rg=0.9985 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(58,17,1:6)=(/ 57.667000,58.355000,59.058000,60.309000,62.836000,67.278000 /)
XPIZA_LKT(58,17,1:6)=(/ 0.728971,0.758779,0.782169,0.906112,0.931708,0.961616 /)
XCGA_LKT(58,17,1:6)=(/ 0.888517,0.880163,0.869537,0.831103,0.813477,0.778830 /)
XEXT_COEFF_550_LKT(58,17)=59.103000 !rg=0.9985 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(58,18,1:6)=(/ 47.921000,48.349000,48.941000,49.937000,51.608000,54.940000 /)
XPIZA_LKT(58,18,1:6)=(/ 0.717174,0.748129,0.773110,0.900906,0.924604,0.954388 /)
XCGA_LKT(58,18,1:6)=(/ 0.891593,0.882683,0.873303,0.835230,0.821523,0.787263 /)
XEXT_COEFF_550_LKT(58,18)=48.855000 !rg=0.9985 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(58,19,1:6)=(/ 39.694000,40.038000,40.484000,41.052000,42.608000,45.122000 /)
XPIZA_LKT(58,19,1:6)=(/ 0.703790,0.736644,0.763582,0.895822,0.918871,0.946020 /)
XCGA_LKT(58,19,1:6)=(/ 0.894460,0.885977,0.876787,0.841810,0.826473,0.792633 /)
XEXT_COEFF_550_LKT(58,19)=40.487000 !rg=0.9985 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(58,20,1:6)=(/ 32.932000,33.220000,33.517000,33.991000,34.995000,36.815000 /)
XPIZA_LKT(58,20,1:6)=(/ 0.690510,0.725418,0.753423,0.889483,0.912473,0.938406 /)
XCGA_LKT(58,20,1:6)=(/ 0.897637,0.889103,0.880523,0.845450,0.834767,0.804827 /)
XEXT_COEFF_550_LKT(58,20)=33.533000 !rg=0.9985 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(59,1,1:6)=(/ 585.780000,677.650000,776.400000,596.740000,1065.800000,386.010000 /)
XPIZA_LKT(59,1,1:6)=(/ 0.851316,0.908914,0.937781,0.978919,0.994381,0.993257 /)
XCGA_LKT(59,1,1:6)=(/ 0.827483,0.815957,0.798303,0.587227,0.803117,0.657853 /)
XEXT_COEFF_550_LKT(59,1)=779.030000 !rg=1.08173 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(59,2,1:6)=(/ 580.150000,607.440000,669.650000,592.780000,1009.900000,398.200000 /)
XPIZA_LKT(59,2,1:6)=(/ 0.856673,0.898885,0.930817,0.978807,0.994087,0.993444 /)
XCGA_LKT(59,2,1:6)=(/ 0.838210,0.807000,0.775503,0.613377,0.798237,0.680183 /)
XEXT_COEFF_550_LKT(59,2)=663.160000 !rg=1.08173 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(59,3,1:6)=(/ 537.420000,557.500000,591.710000,579.980000,898.770000,427.270000 /)
XPIZA_LKT(59,3,1:6)=(/ 0.848344,0.893296,0.922984,0.978587,0.993305,0.993771 /)
XCGA_LKT(59,3,1:6)=(/ 0.839817,0.802830,0.772753,0.656310,0.787537,0.710867 /)
XEXT_COEFF_550_LKT(59,3)=581.880000 !rg=1.08173 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(59,4,1:6)=(/ 484.750000,501.310000,524.570000,544.030000,761.020000,451.770000 /)
XPIZA_LKT(59,4,1:6)=(/ 0.837683,0.886063,0.916780,0.977132,0.992017,0.993972 /)
XCGA_LKT(59,4,1:6)=(/ 0.844050,0.808787,0.774427,0.686587,0.771843,0.735943 /)
XEXT_COEFF_550_LKT(59,4)=516.910000 !rg=1.08173 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(59,5,1:6)=(/ 428.140000,442.670000,455.440000,488.970000,628.640000,457.630000 /)
XPIZA_LKT(59,5,1:6)=(/ 0.825856,0.875540,0.911471,0.974530,0.990304,0.993907 /)
XCGA_LKT(59,5,1:6)=(/ 0.849297,0.818983,0.777560,0.708637,0.757570,0.752160 /)
XEXT_COEFF_550_LKT(59,5)=456.940000 !rg=1.08173 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(59,6,1:6)=(/ 371.670000,382.690000,398.590000,421.510000,517.490000,439.370000 /)
XPIZA_LKT(59,6,1:6)=(/ 0.814152,0.863680,0.900330,0.971712,0.988052,0.993505 /)
XCGA_LKT(59,6,1:6)=(/ 0.854503,0.822680,0.787270,0.725663,0.746617,0.758987 /)
XEXT_COEFF_550_LKT(59,6)=394.450000 !rg=1.08173 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(59,7,1:6)=(/ 319.100000,327.720000,340.310000,360.190000,423.880000,400.830000 /)
XPIZA_LKT(59,7,1:6)=(/ 0.805003,0.850896,0.890479,0.968168,0.985445,0.992728 /)
XCGA_LKT(59,7,1:6)=(/ 0.859397,0.829773,0.797410,0.736783,0.742053,0.758667 /)
XEXT_COEFF_550_LKT(59,7)=336.890000 !rg=1.08173 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(59,8,1:6)=(/ 271.180000,278.570000,288.270000,302.980000,348.180000,352.320000 /)
XPIZA_LKT(59,8,1:6)=(/ 0.796126,0.836370,0.877230,0.962825,0.982165,0.991594 /)
XCGA_LKT(59,8,1:6)=(/ 0.864327,0.837870,0.806193,0.753297,0.741787,0.755310 /)
XEXT_COEFF_550_LKT(59,8)=287.670000 !rg=1.08173 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(59,9,1:6)=(/ 229.120000,235.140000,241.480000,255.130000,283.540000,300.940000 /)
XPIZA_LKT(59,9,1:6)=(/ 0.790255,0.824390,0.864568,0.957604,0.979282,0.990087 /)
XCGA_LKT(59,9,1:6)=(/ 0.867710,0.846003,0.815580,0.764813,0.747667,0.751843 /)
XEXT_COEFF_550_LKT(59,9)=241.120000 !rg=1.08173 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(59,10,1:6)=(/ 193.130000,197.080000,202.450000,210.900000,232.350000,251.540000 /)
XPIZA_LKT(59,10,1:6)=(/ 0.784297,0.814314,0.849464,0.952425,0.975373,0.988083 /)
XCGA_LKT(59,10,1:6)=(/ 0.871587,0.851720,0.826080,0.775493,0.753017,0.748120 /)
XEXT_COEFF_550_LKT(59,10)=202.760000 !rg=1.08173 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(59,11,1:6)=(/ 161.070000,164.610000,168.380000,175.010000,191.850000,209.170000 /)
XPIZA_LKT(59,11,1:6)=(/ 0.778285,0.804382,0.836044,0.944316,0.969104,0.985159 /)
XCGA_LKT(59,11,1:6)=(/ 0.874020,0.858763,0.835657,0.788377,0.759143,0.747033 /)
XEXT_COEFF_550_LKT(59,11)=168.650000 !rg=1.08173 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(59,12,1:6)=(/ 134.580000,136.710000,140.070000,144.640000,155.370000,170.450000 /)
XPIZA_LKT(59,12,1:6)=(/ 0.771705,0.796120,0.823584,0.936237,0.965492,0.982821 /)
XCGA_LKT(59,12,1:6)=(/ 0.877040,0.862957,0.843090,0.796157,0.769897,0.748513 /)
XEXT_COEFF_550_LKT(59,12)=140.080000 !rg=1.08173 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(59,13,1:6)=(/ 111.940000,113.580000,115.610000,119.550000,127.110000,139.230000 /)
XPIZA_LKT(59,13,1:6)=(/ 0.764372,0.788372,0.813531,0.928573,0.959128,0.979525 /)
XCGA_LKT(59,13,1:6)=(/ 0.879413,0.868017,0.850750,0.806410,0.781607,0.754177 /)
XEXT_COEFF_550_LKT(59,13)=115.730000 !rg=1.08173 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(59,14,1:6)=(/ 93.050000,94.275000,96.018000,98.932000,103.840000,113.330000 /)
XPIZA_LKT(59,14,1:6)=(/ 0.755913,0.780143,0.804300,0.921384,0.951533,0.975401 /)
XCGA_LKT(59,14,1:6)=(/ 0.882250,0.871233,0.857273,0.811930,0.791563,0.759300 /)
XEXT_COEFF_550_LKT(59,14)=95.849000 !rg=1.08173 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(59,15,1:6)=(/ 77.125000,78.220000,79.596000,81.600000,85.982000,92.952000 /)
XPIZA_LKT(59,15,1:6)=(/ 0.746425,0.772200,0.794923,0.914734,0.942222,0.969785 /)
XCGA_LKT(59,15,1:6)=(/ 0.884330,0.875127,0.862063,0.821213,0.796100,0.761813 /)
XEXT_COEFF_550_LKT(59,15)=79.471000 !rg=1.08173 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(59,16,1:6)=(/ 63.934000,64.762000,65.668000,67.272000,70.527000,75.935000 /)
XPIZA_LKT(59,16,1:6)=(/ 0.735521,0.763465,0.785760,0.909388,0.934699,0.963768 /)
XCGA_LKT(59,16,1:6)=(/ 0.887013,0.878440,0.867790,0.829097,0.806437,0.770987 /)
XEXT_COEFF_550_LKT(59,16)=65.680000 !rg=1.08173 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(59,17,1:6)=(/ 53.220000,53.659000,54.452000,55.473000,57.377000,61.330000 /)
XPIZA_LKT(59,17,1:6)=(/ 0.723883,0.753634,0.778128,0.903664,0.928379,0.958460 /)
XCGA_LKT(59,17,1:6)=(/ 0.890037,0.881323,0.871970,0.832943,0.818117,0.782933 /)
XEXT_COEFF_550_LKT(59,17)=54.397000 !rg=1.08173 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(59,18,1:6)=(/ 44.090000,44.538000,45.079000,45.901000,47.648000,50.613000 /)
XPIZA_LKT(59,18,1:6)=(/ 0.711187,0.743359,0.768501,0.898236,0.921108,0.949857 /)
XCGA_LKT(59,18,1:6)=(/ 0.892740,0.884353,0.875023,0.838713,0.822220,0.787157 /)
XEXT_COEFF_550_LKT(59,18)=45.057000 !rg=1.08173 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(59,19,1:6)=(/ 36.576000,36.921000,37.291000,37.960000,39.180000,41.442000 /)
XPIZA_LKT(59,19,1:6)=(/ 0.697842,0.731984,0.758337,0.893014,0.915266,0.941797 /)
XCGA_LKT(59,19,1:6)=(/ 0.895813,0.887333,0.879063,0.844033,0.830353,0.797617 /)
XEXT_COEFF_550_LKT(59,19)=37.281000 !rg=1.08173 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(59,20,1:6)=(/ 30.387000,30.581000,30.906000,31.300000,32.062000,33.621000 /)
XPIZA_LKT(59,20,1:6)=(/ 0.684673,0.719332,0.748615,0.886791,0.909689,0.935379 /)
XCGA_LKT(59,20,1:6)=(/ 0.899287,0.890317,0.882390,0.846943,0.838403,0.809217 /)
XEXT_COEFF_550_LKT(59,20)=30.874000 !rg=1.08173 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(60,1,1:6)=(/ 568.900000,602.110000,657.130000,479.260000,1022.700000,412.720000 /)
XPIZA_LKT(60,1,1:6)=(/ 0.853410,0.902301,0.932169,0.973830,0.994158,0.993763 /)
XCGA_LKT(60,1,1:6)=(/ 0.842340,0.816167,0.800420,0.539677,0.802000,0.670367 /)
XEXT_COEFF_550_LKT(60,1)=659.810000 !rg=1.1719 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(60,2,1:6)=(/ 533.440000,555.450000,601.630000,510.790000,948.600000,430.400000 /)
XPIZA_LKT(60,2,1:6)=(/ 0.844946,0.894055,0.919656,0.975653,0.993653,0.993831 /)
XCGA_LKT(60,2,1:6)=(/ 0.840243,0.808830,0.785333,0.606350,0.796027,0.702593 /)
XEXT_COEFF_550_LKT(60,2)=584.330000 !rg=1.1719 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(60,3,1:6)=(/ 493.250000,510.170000,540.810000,523.840000,824.770000,454.780000 /)
XPIZA_LKT(60,3,1:6)=(/ 0.838174,0.886236,0.916388,0.976428,0.992658,0.994038 /)
XCGA_LKT(60,3,1:6)=(/ 0.843820,0.812070,0.772007,0.665740,0.781653,0.729483 /)
XEXT_COEFF_550_LKT(60,3)=536.010000 !rg=1.1719 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(60,4,1:6)=(/ 444.760000,458.630000,480.980000,491.690000,687.170000,469.140000 /)
XPIZA_LKT(60,4,1:6)=(/ 0.827494,0.878194,0.911129,0.975335,0.991135,0.994099 /)
XCGA_LKT(60,4,1:6)=(/ 0.848280,0.817453,0.775437,0.699127,0.764127,0.749137 /)
XEXT_COEFF_550_LKT(60,4)=478.320000 !rg=1.1719 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(60,5,1:6)=(/ 393.510000,404.820000,421.310000,441.700000,565.940000,461.940000 /)
XPIZA_LKT(60,5,1:6)=(/ 0.816755,0.868799,0.904227,0.972357,0.989118,0.993876 /)
XCGA_LKT(60,5,1:6)=(/ 0.854140,0.821450,0.789630,0.713437,0.749160,0.760073 /)
XEXT_COEFF_550_LKT(60,5)=416.050000 !rg=1.1719 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(60,6,1:6)=(/ 342.200000,349.820000,363.520000,381.540000,464.910000,431.320000 /)
XPIZA_LKT(60,6,1:6)=(/ 0.806454,0.856492,0.895012,0.969863,0.986778,0.993312 /)
XCGA_LKT(60,6,1:6)=(/ 0.859413,0.830200,0.793300,0.735400,0.741633,0.762687 /)
XEXT_COEFF_550_LKT(60,6)=364.550000 !rg=1.1719 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(60,7,1:6)=(/ 293.230000,300.520000,311.030000,325.770000,382.150000,384.660000 /)
XPIZA_LKT(60,7,1:6)=(/ 0.797607,0.841991,0.884048,0.966357,0.983922,0.992368 /)
XCGA_LKT(60,7,1:6)=(/ 0.862777,0.836897,0.800550,0.750670,0.738730,0.759840 /)
XEXT_COEFF_550_LKT(60,7)=310.900000 !rg=1.1719 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(60,8,1:6)=(/ 249.950000,255.310000,263.570000,275.630000,313.660000,332.460000 /)
XPIZA_LKT(60,8,1:6)=(/ 0.791864,0.829898,0.869747,0.961165,0.981003,0.991050 /)
XCGA_LKT(60,8,1:6)=(/ 0.867387,0.843443,0.812390,0.759430,0.742020,0.755170 /)
XEXT_COEFF_550_LKT(60,8)=264.090000 !rg=1.1719 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(60,9,1:6)=(/ 211.110000,215.870000,222.320000,232.240000,257.560000,280.520000 /)
XPIZA_LKT(60,9,1:6)=(/ 0.786025,0.817970,0.856219,0.955344,0.977080,0.989271 /)
XCGA_LKT(60,9,1:6)=(/ 0.870247,0.849677,0.823083,0.766863,0.748390,0.750257 /)
XEXT_COEFF_550_LKT(60,9)=220.220000 !rg=1.1719 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(60,10,1:6)=(/ 177.490000,181.150000,185.630000,194.150000,210.970000,232.730000 /)
XPIZA_LKT(60,10,1:6)=(/ 0.780042,0.808016,0.843418,0.947795,0.973515,0.987174 /)
XCGA_LKT(60,10,1:6)=(/ 0.873330,0.856400,0.830970,0.782050,0.758017,0.748540 /)
XEXT_COEFF_550_LKT(60,10)=185.320000 !rg=1.1719 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(60,11,1:6)=(/ 148.490000,151.590000,154.810000,160.900000,172.960000,191.430000 /)
XPIZA_LKT(60,11,1:6)=(/ 0.775319,0.799230,0.830209,0.939488,0.967662,0.984319 /)
XCGA_LKT(60,11,1:6)=(/ 0.875987,0.861193,0.840880,0.790290,0.765217,0.748047 /)
XEXT_COEFF_550_LKT(60,11)=155.020000 !rg=1.1719 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(60,12,1:6)=(/ 123.800000,125.800000,128.390000,132.980000,141.580000,156.330000 /)
XPIZA_LKT(60,12,1:6)=(/ 0.767725,0.791596,0.818936,0.932382,0.962826,0.981563 /)
XCGA_LKT(60,12,1:6)=(/ 0.878477,0.866317,0.847347,0.801950,0.776573,0.751293 /)
XEXT_COEFF_550_LKT(60,12)=128.240000 !rg=1.1719 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(60,13,1:6)=(/ 103.190000,104.420000,106.510000,109.830000,115.790000,127.110000 /)
XPIZA_LKT(60,13,1:6)=(/ 0.760330,0.783675,0.808754,0.924968,0.955609,0.977632 /)
XCGA_LKT(60,13,1:6)=(/ 0.881220,0.869867,0.854757,0.808050,0.786320,0.755457 /)
XEXT_COEFF_550_LKT(60,13)=106.410000 !rg=1.1719 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(60,14,1:6)=(/ 85.625000,86.914000,88.329000,90.845000,95.527000,104.080000 /)
XPIZA_LKT(60,14,1:6)=(/ 0.751635,0.776483,0.799116,0.917699,0.946889,0.972471 /)
XCGA_LKT(60,14,1:6)=(/ 0.883340,0.873670,0.859583,0.817477,0.791853,0.758177 /)
XEXT_COEFF_550_LKT(60,14)=88.389000 !rg=1.1719 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(60,15,1:6)=(/ 71.163000,72.102000,73.146000,74.906000,78.563000,85.175000 /)
XPIZA_LKT(60,15,1:6)=(/ 0.741482,0.767988,0.790622,0.911773,0.939268,0.967445 /)
XCGA_LKT(60,15,1:6)=(/ 0.885990,0.877060,0.865277,0.825387,0.802400,0.766563 /)
XEXT_COEFF_550_LKT(60,15)=72.992000 !rg=1.1719 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(60,16,1:6)=(/ 59.010000,59.661000,60.489000,61.847000,64.107000,68.765000 /)
XPIZA_LKT(60,16,1:6)=(/ 0.730634,0.758978,0.782160,0.905905,0.931990,0.961918 /)
XCGA_LKT(60,16,1:6)=(/ 0.888590,0.879863,0.869990,0.830357,0.812887,0.776850 /)
XEXT_COEFF_550_LKT(60,16)=60.582000 !rg=1.1719 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(60,17,1:6)=(/ 48.976000,49.492000,50.039000,51.129000,52.725000,56.277000 /)
XPIZA_LKT(60,17,1:6)=(/ 0.717959,0.749157,0.773225,0.900906,0.925320,0.955209 /)
XCGA_LKT(60,17,1:6)=(/ 0.891363,0.882987,0.873403,0.835993,0.820127,0.784230 /)
XEXT_COEFF_550_LKT(60,17)=50.065000 !rg=1.1719 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(60,18,1:6)=(/ 40.660000,41.065000,41.503000,42.176000,43.612000,46.291000 /)
XPIZA_LKT(60,18,1:6)=(/ 0.705199,0.738093,0.764239,0.895351,0.918589,0.946599 /)
XCGA_LKT(60,18,1:6)=(/ 0.894270,0.885867,0.877170,0.841040,0.826707,0.792927 /)
XEXT_COEFF_550_LKT(60,18)=41.402000 !rg=1.1719 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(60,19,1:6)=(/ 33.747000,34.017000,34.368000,34.874000,35.786000,37.640000 /)
XPIZA_LKT(60,19,1:6)=(/ 0.692154,0.726397,0.754264,0.889895,0.912795,0.939164 /)
XCGA_LKT(60,19,1:6)=(/ 0.897540,0.888707,0.880587,0.845000,0.834947,0.804027 /)
XEXT_COEFF_550_LKT(60,19)=34.406000 !rg=1.1719 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(60,20,1:6)=(/ 27.992000,28.211000,28.432000,28.893000,29.540000,30.930000 /)
XPIZA_LKT(60,20,1:6)=(/ 0.678323,0.713943,0.742970,0.884074,0.906958,0.931889 /)
XCGA_LKT(60,20,1:6)=(/ 0.900757,0.891850,0.883590,0.849117,0.840013,0.810717 /)
XEXT_COEFF_550_LKT(60,20)=28.442000 !rg=1.1719 sigma=2.95 wvl=0.55
 
IF (LHOOK) CALL DR_HOOK('MODE_DUSTOPT:DUST_OPT_LKT_SET6',1,ZHOOK_HANDLE)
END SUBROUTINE DUST_OPT_LKT_SET6

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE DUST_OPT_LKT_SET7()

  USE MODD_DUST_OPT_LKT
  
  IMPLICIT NONE
 

REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_DUSTOPT:DUST_OPT_LKT_SET7',0,ZHOOK_HANDLE)
XEXT_COEFF_WVL_LKT(61,1,1:6)=(/ 511.090000,502.070000,520.650000,430.770000,946.590000,437.600000 /)
XPIZA_LKT(61,1,1:6)=(/ 0.842644,0.881515,0.918709,0.971421,0.993616,0.993986 /)
XCGA_LKT(61,1,1:6)=(/ 0.844147,0.792317,0.795810,0.562880,0.799427,0.702427 /)
XEXT_COEFF_550_LKT(61,1)=542.750000 !rg=1.26958 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(61,2,1:6)=(/ 490.760000,505.380000,525.210000,478.140000,867.890000,460.360000 /)
XPIZA_LKT(61,2,1:6)=(/ 0.836999,0.887094,0.914290,0.974009,0.993010,0.994117 /)
XCGA_LKT(61,2,1:6)=(/ 0.846640,0.811337,0.780853,0.633940,0.789457,0.726207 /)
XEXT_COEFF_550_LKT(61,2)=523.510000 !rg=1.26958 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(61,3,1:6)=(/ 452.030000,467.400000,488.290000,483.250000,740.690000,477.900000 /)
XPIZA_LKT(61,3,1:6)=(/ 0.829314,0.879829,0.912404,0.974485,0.991755,0.994230 /)
XCGA_LKT(61,3,1:6)=(/ 0.850230,0.818490,0.793013,0.681607,0.771577,0.746073 /)
XEXT_COEFF_550_LKT(61,3)=487.850000 !rg=1.26958 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(61,4,1:6)=(/ 407.810000,420.440000,437.130000,448.850000,612.920000,480.710000 /)
XPIZA_LKT(61,4,1:6)=(/ 0.818154,0.871427,0.906903,0.973244,0.989965,0.994153 /)
XCGA_LKT(61,4,1:6)=(/ 0.854460,0.823887,0.793887,0.709977,0.753403,0.760150 /)
XEXT_COEFF_550_LKT(61,4)=438.490000 !rg=1.26958 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(61,5,1:6)=(/ 360.690000,371.230000,385.950000,400.450000,504.360000,460.040000 /)
XPIZA_LKT(61,5,1:6)=(/ 0.807495,0.859678,0.898330,0.970687,0.987808,0.993775 /)
XCGA_LKT(61,5,1:6)=(/ 0.858093,0.828360,0.790413,0.730787,0.741707,0.766130 /)
XEXT_COEFF_550_LKT(61,5)=384.710000 !rg=1.26958 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(61,6,1:6)=(/ 315.190000,322.860000,331.170000,350.420000,417.650000,418.680000 /)
XPIZA_LKT(61,6,1:6)=(/ 0.799367,0.848000,0.890021,0.966859,0.985263,0.993033 /)
XCGA_LKT(61,6,1:6)=(/ 0.862693,0.835777,0.803820,0.742267,0.738150,0.764910 /)
XEXT_COEFF_550_LKT(61,6)=333.450000 !rg=1.26958 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(61,7,1:6)=(/ 269.340000,276.480000,284.120000,297.650000,345.770000,366.090000 /)
XPIZA_LKT(61,7,1:6)=(/ 0.792410,0.834693,0.876780,0.963348,0.981676,0.991909 /)
XCGA_LKT(61,7,1:6)=(/ 0.866567,0.842713,0.811660,0.756987,0.737463,0.759840 /)
XEXT_COEFF_550_LKT(61,7)=285.510000 !rg=1.26958 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(61,8,1:6)=(/ 230.360000,235.060000,240.960000,253.550000,283.710000,311.640000 /)
XPIZA_LKT(61,8,1:6)=(/ 0.787291,0.822939,0.863607,0.957769,0.979306,0.990405 /)
XCGA_LKT(61,8,1:6)=(/ 0.869737,0.848693,0.819363,0.767207,0.745300,0.754450 /)
XEXT_COEFF_550_LKT(61,8)=241.310000 !rg=1.26958 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(61,9,1:6)=(/ 193.840000,198.390000,203.620000,211.700000,234.330000,259.690000 /)
XPIZA_LKT(61,9,1:6)=(/ 0.781859,0.811609,0.849155,0.951591,0.974697,0.988351 /)
XCGA_LKT(61,9,1:6)=(/ 0.872237,0.854513,0.826497,0.778207,0.747643,0.748697 /)
XEXT_COEFF_550_LKT(61,9)=203.280000 !rg=1.26958 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(61,10,1:6)=(/ 163.570000,166.760000,170.770000,177.600000,191.890000,214.650000 /)
XPIZA_LKT(61,10,1:6)=(/ 0.777137,0.802310,0.835927,0.944335,0.970964,0.986026 /)
XCGA_LKT(61,10,1:6)=(/ 0.875257,0.859353,0.837253,0.783577,0.760923,0.747550 /)
XEXT_COEFF_550_LKT(61,10)=169.480000 !rg=1.26958 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(61,11,1:6)=(/ 136.600000,139.390000,142.450000,147.740000,157.550000,175.590000 /)
XPIZA_LKT(61,11,1:6)=(/ 0.770906,0.794134,0.823688,0.935825,0.965482,0.983044 /)
XCGA_LKT(61,11,1:6)=(/ 0.877373,0.864657,0.844070,0.797017,0.769610,0.748160 /)
XEXT_COEFF_550_LKT(61,11)=142.070000 !rg=1.26958 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(61,12,1:6)=(/ 114.160000,115.850000,118.090000,122.090000,129.130000,143.120000 /)
XPIZA_LKT(61,12,1:6)=(/ 0.764346,0.787050,0.813144,0.928653,0.959540,0.979865 /)
XCGA_LKT(61,12,1:6)=(/ 0.880197,0.868283,0.852033,0.803540,0.781077,0.752113 /)
XEXT_COEFF_550_LKT(61,12)=117.710000 !rg=1.26958 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(61,13,1:6)=(/ 94.836000,96.379000,97.928000,101.030000,106.240000,116.670000 /)
XPIZA_LKT(61,13,1:6)=(/ 0.756140,0.780185,0.803298,0.920847,0.951900,0.975426 /)
XCGA_LKT(61,13,1:6)=(/ 0.882357,0.872350,0.857313,0.813633,0.787353,0.755017 /)
XEXT_COEFF_550_LKT(61,13)=98.007000 !rg=1.26958 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(61,14,1:6)=(/ 78.893000,79.974000,81.301000,83.137000,87.544000,95.542000 /)
XPIZA_LKT(61,14,1:6)=(/ 0.746570,0.771842,0.794762,0.914578,0.943462,0.970273 /)
XCGA_LKT(61,14,1:6)=(/ 0.884890,0.875550,0.863090,0.820980,0.797897,0.762407 /)
XEXT_COEFF_550_LKT(61,14)=81.317000 !rg=1.26958 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(61,15,1:6)=(/ 65.637000,66.276000,67.367000,68.906000,71.644000,77.274000 /)
XPIZA_LKT(61,15,1:6)=(/ 0.736733,0.763865,0.786160,0.909114,0.936432,0.965806 /)
XCGA_LKT(61,15,1:6)=(/ 0.887533,0.878043,0.867973,0.827863,0.808073,0.771473 /)
XEXT_COEFF_550_LKT(61,15)=67.317000 !rg=1.26958 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(61,16,1:6)=(/ 54.335000,55.015000,55.705000,56.857000,58.599000,62.836000 /)
XPIZA_LKT(61,16,1:6)=(/ 0.724834,0.754485,0.777882,0.903513,0.928967,0.959486 /)
XCGA_LKT(61,16,1:6)=(/ 0.889660,0.881850,0.871713,0.833940,0.817157,0.780580 /)
XEXT_COEFF_550_LKT(61,16)=55.656000 !rg=1.26958 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(61,17,1:6)=(/ 45.183000,45.508000,46.138000,46.903000,48.475000,51.690000 /)
XPIZA_LKT(61,17,1:6)=(/ 0.712697,0.744135,0.769122,0.898202,0.922333,0.951492 /)
XCGA_LKT(61,17,1:6)=(/ 0.892773,0.884167,0.875573,0.838580,0.824747,0.788930 /)
XEXT_COEFF_550_LKT(61,17)=46.169000 !rg=1.26958 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(61,18,1:6)=(/ 37.529000,37.779000,38.261000,38.921000,39.944000,42.178000 /)
XPIZA_LKT(61,18,1:6)=(/ 0.699516,0.732918,0.759380,0.893097,0.916137,0.943660 /)
XCGA_LKT(61,18,1:6)=(/ 0.895900,0.886807,0.879050,0.843410,0.831283,0.798427 /)
XEXT_COEFF_550_LKT(61,18)=38.286000 !rg=1.26958 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(61,19,1:6)=(/ 31.097000,31.392000,31.659000,32.141000,32.806000,34.411000 /)
XPIZA_LKT(61,19,1:6)=(/ 0.686010,0.721112,0.749278,0.887303,0.910087,0.936204 /)
XCGA_LKT(61,19,1:6)=(/ 0.898770,0.890447,0.882037,0.847520,0.837900,0.807450 /)
XEXT_COEFF_550_LKT(61,19)=31.643000 !rg=1.26958 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(61,20,1:6)=(/ 25.831000,25.971000,26.225000,26.540000,27.147000,28.383000 /)
XPIZA_LKT(61,20,1:6)=(/ 0.672772,0.707799,0.738096,0.880636,0.904463,0.929209 /)
XCGA_LKT(61,20,1:6)=(/ 0.902323,0.893187,0.885253,0.850820,0.843280,0.815220 /)
XEXT_COEFF_550_LKT(61,20)=26.222000 !rg=1.26958 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(62,1,1:6)=(/ 468.690000,467.500000,428.760000,443.900000,857.720000,474.520000 /)
XPIZA_LKT(62,1,1:6)=(/ 0.832770,0.879949,0.906923,0.972118,0.992898,0.994170 /)
XCGA_LKT(62,1,1:6)=(/ 0.847993,0.807813,0.772673,0.633503,0.791613,0.741077 /)
XEXT_COEFF_550_LKT(62,1)=435.830000 !rg=1.37541 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(62,2,1:6)=(/ 452.580000,464.420000,464.220000,472.630000,771.270000,487.260000 /)
XPIZA_LKT(62,2,1:6)=(/ 0.827678,0.878765,0.912562,0.973315,0.992084,0.994333 /)
XCGA_LKT(62,2,1:6)=(/ 0.851170,0.817977,0.798207,0.673597,0.777680,0.747377 /)
XEXT_COEFF_550_LKT(62,2)=476.490000 !rg=1.37541 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(62,3,1:6)=(/ 416.780000,429.550000,443.690000,453.640000,649.490000,495.600000 /)
XPIZA_LKT(62,3,1:6)=(/ 0.818899,0.872931,0.908334,0.972984,0.990563,0.994349 /)
XCGA_LKT(62,3,1:6)=(/ 0.855573,0.822577,0.798267,0.703503,0.757377,0.760270 /)
XEXT_COEFF_550_LKT(62,3)=448.980000 !rg=1.37541 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(62,4,1:6)=(/ 375.600000,386.940000,399.610000,412.840000,538.290000,485.660000 /)
XPIZA_LKT(62,4,1:6)=(/ 0.808405,0.863705,0.901939,0.971226,0.988570,0.994135 /)
XCGA_LKT(62,4,1:6)=(/ 0.859640,0.827743,0.800713,0.721710,0.741003,0.768923 /)
XEXT_COEFF_550_LKT(62,4)=403.220000 !rg=1.37541 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(62,5,1:6)=(/ 331.170000,340.770000,352.340000,365.270000,449.300000,452.040000 /)
XPIZA_LKT(62,5,1:6)=(/ 0.798969,0.851901,0.892010,0.968644,0.986082,0.993587 /)
XCGA_LKT(62,5,1:6)=(/ 0.863440,0.834710,0.803997,0.740637,0.733560,0.770177 /)
XEXT_COEFF_550_LKT(62,5)=353.380000 !rg=1.37541 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(62,6,1:6)=(/ 289.500000,297.270000,305.800000,319.870000,372.930000,401.490000 /)
XPIZA_LKT(62,6,1:6)=(/ 0.792096,0.838992,0.882388,0.963974,0.983641,0.992665 /)
XCGA_LKT(62,6,1:6)=(/ 0.866600,0.839353,0.810777,0.745217,0.733900,0.765477 /)
XEXT_COEFF_550_LKT(62,6)=304.850000 !rg=1.37541 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(62,7,1:6)=(/ 247.830000,254.640000,260.970000,272.540000,309.350000,344.630000 /)
XPIZA_LKT(62,7,1:6)=(/ 0.787980,0.826624,0.869717,0.959736,0.980393,0.991347 /)
XCGA_LKT(62,7,1:6)=(/ 0.869710,0.846150,0.818573,0.759900,0.738297,0.758757 /)
XEXT_COEFF_550_LKT(62,7)=262.020000 !rg=1.37541 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(62,8,1:6)=(/ 211.830000,216.140000,222.070000,232.180000,256.170000,289.540000 /)
XPIZA_LKT(62,8,1:6)=(/ 0.782568,0.815290,0.855923,0.953956,0.977154,0.989619 /)
XCGA_LKT(62,8,1:6)=(/ 0.872550,0.852197,0.826273,0.767683,0.746293,0.752253 /)
XEXT_COEFF_550_LKT(62,8)=220.650000 !rg=1.37541 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(62,9,1:6)=(/ 178.550000,182.440000,186.820000,193.310000,213.710000,240.110000 /)
XPIZA_LKT(62,9,1:6)=(/ 0.778601,0.805512,0.841453,0.948530,0.971959,0.987106 /)
XCGA_LKT(62,9,1:6)=(/ 0.874960,0.858540,0.834007,0.785410,0.752123,0.747577 /)
XEXT_COEFF_550_LKT(62,9)=186.900000 !rg=1.37541 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(62,10,1:6)=(/ 150.250000,153.120000,156.900000,162.580000,175.750000,197.200000 /)
XPIZA_LKT(62,10,1:6)=(/ 0.774059,0.797739,0.828892,0.939865,0.967154,0.984581 /)
XCGA_LKT(62,10,1:6)=(/ 0.876657,0.862803,0.840513,0.792483,0.759943,0.746053 /)
XEXT_COEFF_550_LKT(62,10)=156.540000 !rg=1.37541 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(62,11,1:6)=(/ 126.070000,128.140000,130.940000,134.950000,144.220000,161.050000 /)
XPIZA_LKT(62,11,1:6)=(/ 0.768349,0.789872,0.817197,0.932657,0.962598,0.981462 /)
XCGA_LKT(62,11,1:6)=(/ 0.879037,0.867347,0.849227,0.800360,0.773813,0.748407 /)
XEXT_COEFF_550_LKT(62,11)=130.070000 !rg=1.37541 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(62,12,1:6)=(/ 104.970000,106.560000,108.680000,112.040000,118.720000,131.130000 /)
XPIZA_LKT(62,12,1:6)=(/ 0.760818,0.783116,0.807446,0.924596,0.955306,0.977400 /)
XCGA_LKT(62,12,1:6)=(/ 0.881197,0.870857,0.854717,0.810153,0.780717,0.750660 /)
XEXT_COEFF_550_LKT(62,12)=108.650000 !rg=1.37541 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(62,13,1:6)=(/ 87.382000,88.546000,90.207000,92.289000,97.482000,107.070000 /)
XPIZA_LKT(62,13,1:6)=(/ 0.751581,0.775903,0.798706,0.917672,0.948050,0.973036 /)
XCGA_LKT(62,13,1:6)=(/ 0.883917,0.873957,0.860997,0.817373,0.793100,0.758240 /)
XEXT_COEFF_550_LKT(62,13)=90.273000 !rg=1.37541 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(62,14,1:6)=(/ 72.724000,73.471000,74.749000,76.287000,79.863000,86.735000 /)
XPIZA_LKT(62,14,1:6)=(/ 0.742557,0.768070,0.789694,0.912363,0.940546,0.968983 /)
XCGA_LKT(62,14,1:6)=(/ 0.885963,0.877467,0.865713,0.825087,0.802663,0.766043 /)
XEXT_COEFF_550_LKT(62,14)=74.803000 !rg=1.37541 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(62,15,1:6)=(/ 60.388000,61.066000,61.914000,63.294000,65.492000,70.510000 /)
XPIZA_LKT(62,15,1:6)=(/ 0.731373,0.759482,0.782306,0.906479,0.933701,0.963651 /)
XCGA_LKT(62,15,1:6)=(/ 0.888457,0.880287,0.869540,0.831347,0.814083,0.776917 /)
XEXT_COEFF_550_LKT(62,15)=61.858000 !rg=1.37541 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(62,16,1:6)=(/ 50.091000,50.606000,51.307000,52.195000,54.158000,57.824000 /)
XPIZA_LKT(62,16,1:6)=(/ 0.719345,0.749893,0.773488,0.900912,0.925842,0.956036 /)
XCGA_LKT(62,16,1:6)=(/ 0.891230,0.883063,0.874253,0.836330,0.821010,0.784503 /)
XEXT_COEFF_550_LKT(62,16)=51.241000 !rg=1.37541 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(62,17,1:6)=(/ 41.622000,41.984000,42.396000,43.035000,44.523000,47.211000 /)
XPIZA_LKT(62,17,1:6)=(/ 0.706617,0.739286,0.764664,0.895731,0.919484,0.948476 /)
XCGA_LKT(62,17,1:6)=(/ 0.894287,0.885967,0.877033,0.841880,0.827240,0.792340 /)
XEXT_COEFF_550_LKT(62,17)=42.492000 !rg=1.37541 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(62,18,1:6)=(/ 34.547000,34.847000,35.171000,35.767000,36.656000,38.515000 /)
XPIZA_LKT(62,18,1:6)=(/ 0.693234,0.727734,0.754907,0.890584,0.913790,0.941042 /)
XCGA_LKT(62,18,1:6)=(/ 0.897180,0.888750,0.880103,0.845747,0.835687,0.804247 /)
XEXT_COEFF_550_LKT(62,18)=35.141000 !rg=1.37541 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(62,19,1:6)=(/ 28.669000,28.882000,29.200000,29.505000,30.314000,31.759000 /)
XPIZA_LKT(62,19,1:6)=(/ 0.679865,0.715114,0.744404,0.884149,0.907451,0.933012 /)
XCGA_LKT(62,19,1:6)=(/ 0.900417,0.891697,0.883973,0.849010,0.841180,0.811807 /)
XEXT_COEFF_550_LKT(62,19)=29.159000 !rg=1.37541 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(62,20,1:6)=(/ 23.811000,23.957000,24.129000,24.396000,24.983000,26.033000 /)
XPIZA_LKT(62,20,1:6)=(/ 0.666463,0.702123,0.732414,0.877538,0.901590,0.926620 /)
XCGA_LKT(62,20,1:6)=(/ 0.904007,0.894813,0.886607,0.852800,0.845117,0.818343 /)
XEXT_COEFF_550_LKT(62,20)=24.166000 !rg=1.37541 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(63,1,1:6)=(/ 437.190000,463.930000,441.400000,478.300000,744.750000,509.680000 /)
XPIZA_LKT(63,1,1:6)=(/ 0.823671,0.881232,0.897153,0.974358,0.991785,0.994466 /)
XCGA_LKT(63,1,1:6)=(/ 0.855760,0.833647,0.791067,0.712330,0.773503,0.763577 /)
XEXT_COEFF_550_LKT(63,1)=430.720000 !rg=1.49006 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(63,2,1:6)=(/ 413.650000,429.490000,435.750000,468.150000,664.730000,509.020000 /)
XPIZA_LKT(63,2,1:6)=(/ 0.817692,0.869756,0.906490,0.973296,0.990756,0.994487 /)
XCGA_LKT(63,2,1:6)=(/ 0.856187,0.824047,0.802363,0.718913,0.759183,0.763937 /)
XEXT_COEFF_550_LKT(63,2)=436.640000 !rg=1.49006 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(63,3,1:6)=(/ 383.650000,396.150000,404.280000,428.420000,559.890000,506.740000 /)
XPIZA_LKT(63,3,1:6)=(/ 0.808700,0.863980,0.902754,0.970615,0.989002,0.994395 /)
XCGA_LKT(63,3,1:6)=(/ 0.860893,0.829180,0.799823,0.727183,0.740007,0.771813 /)
XEXT_COEFF_550_LKT(63,3)=406.930000 !rg=1.49006 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(63,4,1:6)=(/ 345.510000,356.270000,365.360000,384.120000,470.630000,483.440000 /)
XPIZA_LKT(63,4,1:6)=(/ 0.798834,0.854644,0.895773,0.967845,0.986905,0.994036 /)
XCGA_LKT(63,4,1:6)=(/ 0.864583,0.834100,0.803387,0.737547,0.728953,0.775403 /)
XEXT_COEFF_550_LKT(63,4)=366.720000 !rg=1.49006 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(63,5,1:6)=(/ 305.010000,313.600000,322.460000,334.580000,396.590000,437.830000 /)
XPIZA_LKT(63,5,1:6)=(/ 0.791906,0.843604,0.886216,0.966320,0.984499,0.993315 /)
XCGA_LKT(63,5,1:6)=(/ 0.867570,0.838827,0.810930,0.747387,0.727023,0.772307 /)
XEXT_COEFF_550_LKT(63,5)=324.230000 !rg=1.49006 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(63,6,1:6)=(/ 266.020000,272.650000,281.170000,294.540000,334.240000,380.590000 /)
XPIZA_LKT(63,6,1:6)=(/ 0.786484,0.831053,0.873898,0.960051,0.981862,0.992206 /)
XCGA_LKT(63,6,1:6)=(/ 0.870247,0.845497,0.813047,0.757373,0.731310,0.765177 /)
XEXT_COEFF_550_LKT(63,6)=279.210000 !rg=1.49006 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(63,7,1:6)=(/ 228.070000,233.900000,239.980000,251.160000,279.230000,321.700000 /)
XPIZA_LKT(63,7,1:6)=(/ 0.782574,0.818491,0.861655,0.956201,0.978492,0.990692 /)
XCGA_LKT(63,7,1:6)=(/ 0.872400,0.851473,0.822300,0.769507,0.739383,0.757220 /)
XEXT_COEFF_550_LKT(63,7)=239.610000 !rg=1.49006 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(63,8,1:6)=(/ 194.520000,198.810000,204.000000,212.810000,232.830000,267.260000 /)
XPIZA_LKT(63,8,1:6)=(/ 0.779251,0.808902,0.847658,0.950073,0.974645,0.988708 /)
XCGA_LKT(63,8,1:6)=(/ 0.874577,0.856640,0.829450,0.777943,0.745907,0.750250 /)
XEXT_COEFF_550_LKT(63,8)=203.280000 !rg=1.49006 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(63,9,1:6)=(/ 164.510000,167.970000,171.540000,177.370000,192.320000,219.690000 /)
XPIZA_LKT(63,9,1:6)=(/ 0.775305,0.800315,0.834400,0.944669,0.970630,0.986287 /)
XCGA_LKT(63,9,1:6)=(/ 0.876913,0.861397,0.839600,0.788627,0.757870,0.747280 /)
XEXT_COEFF_550_LKT(63,9)=171.760000 !rg=1.49006 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(63,10,1:6)=(/ 138.390000,141.150000,143.990000,148.530000,160.610000,180.930000 /)
XPIZA_LKT(63,10,1:6)=(/ 0.769804,0.792029,0.822570,0.936096,0.964609,0.982879 /)
XCGA_LKT(63,10,1:6)=(/ 0.878610,0.865943,0.846157,0.798893,0.766410,0.746020 /)
XEXT_COEFF_550_LKT(63,10)=144.130000 !rg=1.49006 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(63,11,1:6)=(/ 115.920000,117.930000,120.160000,123.390000,132.190000,147.270000 /)
XPIZA_LKT(63,11,1:6)=(/ 0.763887,0.785099,0.811879,0.928671,0.959290,0.979835 /)
XCGA_LKT(63,11,1:6)=(/ 0.880367,0.870300,0.852503,0.807843,0.774247,0.747377 /)
XEXT_COEFF_550_LKT(63,11)=119.890000 !rg=1.49006 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(63,12,1:6)=(/ 96.678000,98.185000,99.868000,102.490000,108.690000,120.220000 /)
XPIZA_LKT(63,12,1:6)=(/ 0.756109,0.778637,0.802451,0.921044,0.952134,0.975428 /)
XCGA_LKT(63,12,1:6)=(/ 0.882787,0.873093,0.858730,0.814437,0.787337,0.754003 /)
XEXT_COEFF_550_LKT(63,12)=100.020000 !rg=1.49006 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(63,13,1:6)=(/ 80.498000,81.387000,82.816000,84.607000,88.910000,97.244000 /)
XPIZA_LKT(63,13,1:6)=(/ 0.747701,0.771674,0.793890,0.915135,0.945219,0.971814 /)
XCGA_LKT(63,13,1:6)=(/ 0.884853,0.876407,0.863480,0.822460,0.797610,0.761020 /)
XEXT_COEFF_550_LKT(63,13)=82.986000 !rg=1.49006 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(63,14,1:6)=(/ 66.984000,67.769000,68.622000,70.161000,73.227000,79.189000 /)
XPIZA_LKT(63,14,1:6)=(/ 0.737500,0.764608,0.785918,0.909068,0.936620,0.966447 /)
XCGA_LKT(63,14,1:6)=(/ 0.887200,0.879137,0.868543,0.829157,0.808840,0.771953 /)
XEXT_COEFF_550_LKT(63,14)=68.621000 !rg=1.49006 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(63,15,1:6)=(/ 55.707000,56.281000,56.988000,58.213000,60.212000,64.629000 /)
XPIZA_LKT(63,15,1:6)=(/ 0.726503,0.755043,0.777783,0.903610,0.929206,0.960190 /)
XCGA_LKT(63,15,1:6)=(/ 0.889873,0.881607,0.872460,0.833283,0.817457,0.780250 /)
XEXT_COEFF_550_LKT(63,15)=56.873000 !rg=1.49006 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(63,16,1:6)=(/ 46.176000,46.654000,47.169000,48.020000,49.770000,53.059000 /)
XPIZA_LKT(63,16,1:6)=(/ 0.713564,0.744913,0.769417,0.898259,0.922531,0.952497 /)
XCGA_LKT(63,16,1:6)=(/ 0.892613,0.884737,0.875800,0.839767,0.822063,0.785130 /)
XEXT_COEFF_550_LKT(63,16)=47.098000 !rg=1.49006 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(63,17,1:6)=(/ 38.357000,38.747000,39.099000,39.626000,40.905000,43.285000 /)
XPIZA_LKT(63,17,1:6)=(/ 0.700631,0.734424,0.760276,0.892877,0.916445,0.944416 /)
XCGA_LKT(63,17,1:6)=(/ 0.895427,0.887580,0.879373,0.843903,0.831587,0.798390 /)
XEXT_COEFF_550_LKT(63,17)=39.114000 !rg=1.49006 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(63,18,1:6)=(/ 31.885000,32.115000,32.416000,32.877000,33.692000,35.417000 /)
XPIZA_LKT(63,18,1:6)=(/ 0.687627,0.722043,0.750224,0.887468,0.909860,0.936585 /)
XCGA_LKT(63,18,1:6)=(/ 0.898660,0.889957,0.882270,0.846633,0.837930,0.807387 /)
XEXT_COEFF_550_LKT(63,18)=32.377000 !rg=1.49006 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(63,19,1:6)=(/ 26.450000,26.616000,26.846000,27.208000,27.896000,29.195000 /)
XPIZA_LKT(63,19,1:6)=(/ 0.673863,0.709056,0.739005,0.881520,0.904785,0.929902 /)
XCGA_LKT(63,19,1:6)=(/ 0.902017,0.893273,0.885240,0.851480,0.841873,0.812410 /)
XEXT_COEFF_550_LKT(63,19)=26.838000 !rg=1.49006 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(63,20,1:6)=(/ 21.950000,22.108000,22.256000,22.475000,22.958000,23.904000 /)
XPIZA_LKT(63,20,1:6)=(/ 0.660450,0.696347,0.727248,0.873998,0.898907,0.923204 /)
XCGA_LKT(63,20,1:6)=(/ 0.905437,0.896497,0.888387,0.854417,0.848050,0.822793 /)
XEXT_COEFF_550_LKT(63,20)=22.260000 !rg=1.49006 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(64,1,1:6)=(/ 398.630000,408.420000,449.150000,499.530000,624.760000,523.690000 /)
XPIZA_LKT(64,1,1:6)=(/ 0.811429,0.866921,0.910631,0.975302,0.990151,0.994661 /)
XCGA_LKT(64,1,1:6)=(/ 0.857763,0.831697,0.810147,0.765477,0.751050,0.770760 /)
XEXT_COEFF_550_LKT(64,1)=455.020000 !rg=1.61427 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(64,2,1:6)=(/ 380.260000,391.630000,406.400000,452.330000,555.770000,523.900000 /)
XPIZA_LKT(64,2,1:6)=(/ 0.806132,0.864511,0.899904,0.971621,0.988921,0.994571 /)
XCGA_LKT(64,2,1:6)=(/ 0.861950,0.832700,0.797993,0.751623,0.734453,0.776377 /)
XEXT_COEFF_550_LKT(64,2)=405.900000 !rg=1.61427 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(64,3,1:6)=(/ 353.840000,361.430000,374.220000,398.130000,476.660000,510.270000 /)
XPIZA_LKT(64,3,1:6)=(/ 0.797306,0.857290,0.895270,0.968502,0.986963,0.994363 /)
XCGA_LKT(64,3,1:6)=(/ 0.864873,0.835293,0.808413,0.739200,0.717717,0.780600 /)
XEXT_COEFF_550_LKT(64,3)=370.540000 !rg=1.61427 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(64,4,1:6)=(/ 318.610000,325.490000,337.400000,353.160000,411.060000,473.770000 /)
XPIZA_LKT(64,4,1:6)=(/ 0.789615,0.846933,0.887895,0.965858,0.984853,0.993845 /)
XCGA_LKT(64,4,1:6)=(/ 0.868210,0.839397,0.812083,0.742897,0.715503,0.779430 /)
XEXT_COEFF_550_LKT(64,4)=333.930000 !rg=1.61427 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(64,5,1:6)=(/ 280.520000,288.290000,295.730000,310.050000,352.260000,418.280000 /)
XPIZA_LKT(64,5,1:6)=(/ 0.784766,0.833994,0.879405,0.962465,0.982656,0.992939 /)
XCGA_LKT(64,5,1:6)=(/ 0.871223,0.844630,0.814733,0.759420,0.723280,0.772677 /)
XEXT_COEFF_550_LKT(64,5)=295.840000 !rg=1.61427 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(64,6,1:6)=(/ 244.920000,249.700000,257.780000,267.050000,301.020000,357.090000 /)
XPIZA_LKT(64,6,1:6)=(/ 0.782156,0.822509,0.866844,0.958814,0.979495,0.991598 /)
XCGA_LKT(64,6,1:6)=(/ 0.872970,0.849430,0.822180,0.765910,0.730570,0.762880 /)
XEXT_COEFF_550_LKT(64,6)=256.660000 !rg=1.61427 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(64,7,1:6)=(/ 210.220000,214.560000,220.780000,228.880000,252.650000,297.820000 /)
XPIZA_LKT(64,7,1:6)=(/ 0.779533,0.811484,0.853312,0.954245,0.976194,0.989845 /)
XCGA_LKT(64,7,1:6)=(/ 0.874400,0.855630,0.829743,0.772473,0.740883,0.753990 /)
XEXT_COEFF_550_LKT(64,7)=219.040000 !rg=1.61427 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(64,8,1:6)=(/ 179.070000,182.500000,187.310000,193.750000,212.250000,246.030000 /)
XPIZA_LKT(64,8,1:6)=(/ 0.774972,0.801708,0.839865,0.946520,0.971891,0.987486 /)
XCGA_LKT(64,8,1:6)=(/ 0.876760,0.860003,0.836383,0.784897,0.750423,0.747877 /)
XEXT_COEFF_550_LKT(64,8)=187.490000 !rg=1.61427 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(64,9,1:6)=(/ 151.230000,154.200000,157.720000,163.500000,174.670000,201.100000 /)
XPIZA_LKT(64,9,1:6)=(/ 0.771677,0.794349,0.827803,0.939729,0.968617,0.985153 /)
XCGA_LKT(64,9,1:6)=(/ 0.878210,0.864987,0.843477,0.794983,0.764053,0.747087 /)
XEXT_COEFF_550_LKT(64,9)=157.150000 !rg=1.61427 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(64,10,1:6)=(/ 127.750000,129.600000,132.430000,136.150000,145.200000,164.180000 /)
XPIZA_LKT(64,10,1:6)=(/ 0.767166,0.788574,0.815441,0.932935,0.963151,0.982115 /)
XCGA_LKT(64,10,1:6)=(/ 0.880050,0.868247,0.850540,0.803127,0.772347,0.747170 /)
XEXT_COEFF_550_LKT(64,10)=132.400000 !rg=1.61427 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(64,11,1:6)=(/ 106.600000,108.480000,110.440000,113.470000,121.090000,134.880000 /)
XPIZA_LKT(64,11,1:6)=(/ 0.760366,0.781804,0.805440,0.925352,0.955225,0.977211 /)
XCGA_LKT(64,11,1:6)=(/ 0.881827,0.872467,0.857253,0.813303,0.780557,0.749773 /)
XEXT_COEFF_550_LKT(64,11)=110.320000 !rg=1.61427 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(64,12,1:6)=(/ 89.184000,90.207000,91.890000,94.009000,98.894000,109.010000 /)
XPIZA_LKT(64,12,1:6)=(/ 0.752978,0.775371,0.796940,0.918490,0.949317,0.974343 /)
XCGA_LKT(64,12,1:6)=(/ 0.883803,0.875087,0.861903,0.818857,0.792653,0.757083 /)
XEXT_COEFF_550_LKT(64,12)=91.910000 !rg=1.61427 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(64,13,1:6)=(/ 74.152000,75.063000,76.046000,77.936000,81.495000,88.729000 /)
XPIZA_LKT(64,13,1:6)=(/ 0.743000,0.768385,0.789288,0.911741,0.940917,0.969405 /)
XCGA_LKT(64,13,1:6)=(/ 0.886093,0.878200,0.867030,0.826423,0.804217,0.766920 /)
XEXT_COEFF_550_LKT(64,13)=76.121000 !rg=1.61427 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(64,14,1:6)=(/ 61.836000,62.340000,63.270000,64.611000,66.884000,72.221000 /)
XPIZA_LKT(64,14,1:6)=(/ 0.732803,0.759647,0.782383,0.906319,0.933234,0.963536 /)
XCGA_LKT(64,14,1:6)=(/ 0.888753,0.880340,0.870990,0.830900,0.813767,0.775693 /)
XEXT_COEFF_550_LKT(64,14)=63.281000 !rg=1.61427 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(64,15,1:6)=(/ 51.311000,51.806000,52.502000,53.480000,55.630000,59.491000 /)
XPIZA_LKT(64,15,1:6)=(/ 0.720858,0.751085,0.773620,0.900961,0.925049,0.955861 /)
XCGA_LKT(64,15,1:6)=(/ 0.890980,0.882993,0.874130,0.837117,0.817817,0.779640 /)
XEXT_COEFF_550_LKT(64,15)=52.474000 !rg=1.61427 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(64,16,1:6)=(/ 42.538000,42.966000,43.467000,44.231000,45.764000,48.709000 /)
XPIZA_LKT(64,16,1:6)=(/ 0.707766,0.740435,0.764693,0.895981,0.918793,0.947671 /)
XCGA_LKT(64,16,1:6)=(/ 0.893970,0.885970,0.878013,0.842543,0.826197,0.790320 /)
XEXT_COEFF_550_LKT(64,16)=43.392000 !rg=1.61427 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(64,17,1:6)=(/ 35.409000,35.639000,36.064000,36.596000,37.453000,39.444000 /)
XPIZA_LKT(64,17,1:6)=(/ 0.694937,0.728576,0.755822,0.890497,0.913131,0.940791 /)
XCGA_LKT(64,17,1:6)=(/ 0.896997,0.888653,0.881237,0.845557,0.835563,0.803227 /)
XEXT_COEFF_550_LKT(64,17)=36.028000 !rg=1.61427 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(64,18,1:6)=(/ 29.380000,29.587000,29.862000,30.292000,31.120000,32.690000 /)
XPIZA_LKT(64,18,1:6)=(/ 0.681313,0.716820,0.745026,0.884876,0.906812,0.932579 /)
XCGA_LKT(64,18,1:6)=(/ 0.900100,0.891370,0.883580,0.849637,0.838797,0.807700 /)
XEXT_COEFF_550_LKT(64,18)=29.868000 !rg=1.61427 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(64,19,1:6)=(/ 24.376000,24.539000,24.733000,25.064000,25.667000,26.827000 /)
XPIZA_LKT(64,19,1:6)=(/ 0.667727,0.703516,0.733386,0.878447,0.901679,0.926109 /)
XCGA_LKT(64,19,1:6)=(/ 0.903537,0.894537,0.886967,0.853267,0.844713,0.817060 /)
XEXT_COEFF_550_LKT(64,19)=24.724000 !rg=1.61427 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(64,20,1:6)=(/ 20.264000,20.358000,20.535000,20.732000,21.115000,21.869000 /)
XPIZA_LKT(64,20,1:6)=(/ 0.654824,0.690072,0.721838,0.870609,0.896029,0.920287 /)
XCGA_LKT(64,20,1:6)=(/ 0.907200,0.897763,0.890103,0.855770,0.851060,0.826950 /)
XEXT_COEFF_550_LKT(64,20)=20.512000 !rg=1.61427 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(65,1,1:6)=(/ 365.740000,359.430000,419.120000,484.780000,506.130000,534.960000 /)
XPIZA_LKT(65,1,1:6)=(/ 0.802088,0.855063,0.906049,0.973370,0.987728,0.994653 /)
XCGA_LKT(65,1,1:6)=(/ 0.864407,0.831057,0.812623,0.786393,0.714127,0.779373 /)
XEXT_COEFF_550_LKT(65,1)=419.280000 !rg=1.74883 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(65,2,1:6)=(/ 350.250000,358.680000,373.860000,414.510000,455.940000,530.880000 /)
XPIZA_LKT(65,2,1:6)=(/ 0.795447,0.853948,0.895670,0.970241,0.986332,0.994577 /)
XCGA_LKT(65,2,1:6)=(/ 0.868217,0.836007,0.808687,0.768913,0.700077,0.786000 /)
XEXT_COEFF_550_LKT(65,2)=372.240000 !rg=1.74883 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(65,3,1:6)=(/ 323.700000,332.020000,343.360000,362.970000,403.230000,505.350000 /)
XPIZA_LKT(65,3,1:6)=(/ 0.786565,0.846746,0.889557,0.966968,0.984645,0.994244 /)
XCGA_LKT(65,3,1:6)=(/ 0.869570,0.838890,0.807127,0.765777,0.698737,0.786587 /)
XEXT_COEFF_550_LKT(65,3)=342.630000 !rg=1.74883 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(65,4,1:6)=(/ 291.850000,299.240000,308.730000,321.230000,359.190000,456.600000 /)
XPIZA_LKT(65,4,1:6)=(/ 0.781324,0.836500,0.881865,0.964353,0.982746,0.993556 /)
XCGA_LKT(65,4,1:6)=(/ 0.871880,0.843750,0.812817,0.766833,0.707673,0.781110 /)
XEXT_COEFF_550_LKT(65,4)=308.580000 !rg=1.74883 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(65,5,1:6)=(/ 258.660000,264.200000,272.750000,284.190000,314.450000,394.100000 /)
XPIZA_LKT(65,5,1:6)=(/ 0.779189,0.825412,0.870884,0.960399,0.980372,0.992429 /)
XCGA_LKT(65,5,1:6)=(/ 0.873703,0.849780,0.822340,0.760283,0.719850,0.770683 /)
XEXT_COEFF_550_LKT(65,5)=269.930000 !rg=1.74883 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(65,6,1:6)=(/ 225.400000,229.170000,235.770000,243.070000,270.230000,330.670000 /)
XPIZA_LKT(65,6,1:6)=(/ 0.777723,0.813957,0.858975,0.956742,0.977497,0.990910 /)
XCGA_LKT(65,6,1:6)=(/ 0.875287,0.855117,0.826240,0.777517,0.730430,0.760057 /)
XEXT_COEFF_550_LKT(65,6)=236.560000 !rg=1.74883 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(65,7,1:6)=(/ 193.200000,197.290000,202.080000,208.490000,228.940000,273.150000 /)
XPIZA_LKT(65,7,1:6)=(/ 0.775117,0.804274,0.845866,0.951077,0.973867,0.988900 /)
XCGA_LKT(65,7,1:6)=(/ 0.876297,0.859797,0.833393,0.785437,0.740890,0.750937 /)
XEXT_COEFF_550_LKT(65,7)=202.040000 !rg=1.74883 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(65,8,1:6)=(/ 164.940000,167.560000,171.740000,176.810000,191.550000,224.280000 /)
XPIZA_LKT(65,8,1:6)=(/ 0.773003,0.796497,0.831697,0.944382,0.970528,0.986547 /)
XCGA_LKT(65,8,1:6)=(/ 0.877980,0.864060,0.840890,0.791750,0.754730,0.746230 /)
XEXT_COEFF_550_LKT(65,8)=172.000000 !rg=1.74883 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(65,9,1:6)=(/ 139.490000,142.060000,144.960000,149.590000,159.400000,183.490000 /)
XPIZA_LKT(65,9,1:6)=(/ 0.768960,0.789886,0.820181,0.936738,0.965463,0.983613 /)
XCGA_LKT(65,9,1:6)=(/ 0.879530,0.867987,0.848817,0.797530,0.768440,0.745407 /)
XEXT_COEFF_550_LKT(65,9)=143.860000 !rg=1.74883 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(65,10,1:6)=(/ 117.390000,119.150000,121.520000,125.260000,132.250000,149.690000 /)
XPIZA_LKT(65,10,1:6)=(/ 0.763606,0.784464,0.810645,0.928788,0.960517,0.980709 /)
XCGA_LKT(65,10,1:6)=(/ 0.881187,0.871183,0.854263,0.808867,0.780170,0.749730 /)
XEXT_COEFF_550_LKT(65,10)=121.240000 !rg=1.74883 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(65,11,1:6)=(/ 98.388000,99.937000,101.580000,104.470000,109.520000,121.900000 /)
XPIZA_LKT(65,11,1:6)=(/ 0.757270,0.777883,0.801045,0.921298,0.952990,0.976109 /)
XCGA_LKT(65,11,1:6)=(/ 0.883147,0.874017,0.860813,0.815387,0.788343,0.753503 /)
XEXT_COEFF_550_LKT(65,11)=101.570000 !rg=1.74883 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(65,12,1:6)=(/ 82.133000,83.182000,84.352000,86.429000,90.469000,99.242000 /)
XPIZA_LKT(65,12,1:6)=(/ 0.748511,0.771935,0.793094,0.914930,0.945327,0.972302 /)
XCGA_LKT(65,12,1:6)=(/ 0.885097,0.877097,0.865013,0.823657,0.799810,0.762293 /)
XEXT_COEFF_550_LKT(65,12)=84.290000 !rg=1.74883 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(65,13,1:6)=(/ 68.468000,69.063000,70.106000,71.615000,74.339000,80.774000 /)
XPIZA_LKT(65,13,1:6)=(/ 0.738716,0.763975,0.785689,0.909158,0.937199,0.966575 /)
XCGA_LKT(65,13,1:6)=(/ 0.887553,0.879430,0.869597,0.828310,0.809453,0.770513 /)
XEXT_COEFF_550_LKT(65,13)=70.118000 !rg=1.74883 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(65,14,1:6)=(/ 56.878000,57.600000,58.189000,59.464000,61.622000,66.292000 /)
XPIZA_LKT(65,14,1:6)=(/ 0.727229,0.756363,0.777475,0.903033,0.929232,0.960085 /)
XCGA_LKT(65,14,1:6)=(/ 0.889797,0.882087,0.872653,0.834527,0.814987,0.775533 /)
XEXT_COEFF_550_LKT(65,14)=58.308000 !rg=1.74883 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(65,15,1:6)=(/ 47.273000,47.781000,48.274000,49.139000,50.883000,54.384000 /)
XPIZA_LKT(65,15,1:6)=(/ 0.714829,0.745807,0.769928,0.897987,0.922408,0.953005 /)
XCGA_LKT(65,15,1:6)=(/ 0.892367,0.884700,0.876087,0.839423,0.823010,0.785860 /)
XEXT_COEFF_550_LKT(65,15)=48.195000 !rg=1.74883 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(65,16,1:6)=(/ 39.266000,39.587000,40.021000,40.650000,41.771000,44.149000 /)
XPIZA_LKT(65,16,1:6)=(/ 0.702237,0.735155,0.760958,0.893051,0.916074,0.944860 /)
XCGA_LKT(65,16,1:6)=(/ 0.895453,0.887140,0.879570,0.843857,0.831437,0.797293 /)
XEXT_COEFF_550_LKT(65,16)=40.036000 !rg=1.74883 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(65,17,1:6)=(/ 32.596000,32.897000,33.157000,33.704000,34.447000,36.197000 /)
XPIZA_LKT(65,17,1:6)=(/ 0.688452,0.723723,0.750665,0.888096,0.910407,0.937466 /)
XCGA_LKT(65,17,1:6)=(/ 0.898463,0.890163,0.882330,0.847460,0.837337,0.805280 /)
XEXT_COEFF_550_LKT(65,17)=33.125000 !rg=1.74883 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(65,18,1:6)=(/ 27.091000,27.305000,27.487000,27.841000,28.515000,29.896000 /)
XPIZA_LKT(65,18,1:6)=(/ 0.675178,0.710836,0.740382,0.881704,0.904459,0.929924 /)
XCGA_LKT(65,18,1:6)=(/ 0.901573,0.892930,0.885027,0.851143,0.842133,0.812887 /)
XEXT_COEFF_550_LKT(65,18)=27.466000 !rg=1.74883 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(65,19,1:6)=(/ 22.497000,22.622000,22.797000,23.054000,23.544000,24.439000 /)
XPIZA_LKT(65,19,1:6)=(/ 0.661934,0.697476,0.728443,0.874925,0.899001,0.923628 /)
XCGA_LKT(65,19,1:6)=(/ 0.905233,0.895990,0.888213,0.854267,0.848103,0.822347 /)
XEXT_COEFF_550_LKT(65,19)=22.840000 !rg=1.74883 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(65,20,1:6)=(/ 18.669000,18.789000,18.897000,19.132000,19.425000,20.112000 /)
XPIZA_LKT(65,20,1:6)=(/ 0.648747,0.684448,0.715794,0.867373,0.893166,0.917434 /)
XCGA_LKT(65,20,1:6)=(/ 0.908817,0.899333,0.891243,0.857300,0.852103,0.828617 /)
XEXT_COEFF_550_LKT(65,20)=18.888000 !rg=1.74883 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(66,1,1:6)=(/ 336.570000,356.200000,349.900000,415.980000,391.390000,547.640000 /)
XPIZA_LKT(66,1,1:6)=(/ 0.789558,0.853688,0.888681,0.970478,0.984200,0.994689 /)
XCGA_LKT(66,1,1:6)=(/ 0.870337,0.842947,0.798413,0.778583,0.664070,0.792827 /)
XEXT_COEFF_550_LKT(66,1)=345.600000 !rg=1.89461 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(66,2,1:6)=(/ 323.030000,328.830000,340.270000,362.830000,367.900000,528.920000 /)
XPIZA_LKT(66,2,1:6)=(/ 0.782287,0.845822,0.890010,0.968696,0.983186,0.994495 /)
XCGA_LKT(66,2,1:6)=(/ 0.872387,0.841683,0.811320,0.782583,0.665647,0.793167 /)
XEXT_COEFF_550_LKT(66,2)=340.190000 !rg=1.89461 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(66,3,1:6)=(/ 297.650000,305.330000,314.570000,328.910000,347.350000,491.480000 /)
XPIZA_LKT(66,3,1:6)=(/ 0.775436,0.837889,0.881686,0.965554,0.981985,0.994021 /)
XCGA_LKT(66,3,1:6)=(/ 0.875070,0.845070,0.814700,0.776487,0.682457,0.789690 /)
XEXT_COEFF_550_LKT(66,3)=315.830000 !rg=1.89461 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(66,4,1:6)=(/ 268.550000,275.360000,282.800000,293.500000,319.160000,432.790000 /)
XPIZA_LKT(66,4,1:6)=(/ 0.774480,0.827309,0.873376,0.962439,0.980161,0.993135 /)
XCGA_LKT(66,4,1:6)=(/ 0.876007,0.849740,0.820347,0.776287,0.703503,0.780113 /)
XEXT_COEFF_550_LKT(66,4)=284.150000 !rg=1.89461 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(66,5,1:6)=(/ 237.270000,242.870000,249.520000,258.870000,281.830000,365.860000 /)
XPIZA_LKT(66,5,1:6)=(/ 0.774588,0.815878,0.863877,0.957647,0.978149,0.991811 /)
XCGA_LKT(66,5,1:6)=(/ 0.876123,0.854373,0.825237,0.778117,0.719753,0.767550 /)
XEXT_COEFF_550_LKT(66,5)=249.470000 !rg=1.89461 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(66,6,1:6)=(/ 207.650000,211.380000,215.660000,224.840000,244.270000,303.860000 /)
XPIZA_LKT(66,6,1:6)=(/ 0.774711,0.806244,0.851297,0.953087,0.975128,0.990039 /)
XCGA_LKT(66,6,1:6)=(/ 0.876973,0.859670,0.833593,0.783487,0.737437,0.756173 /)
XEXT_COEFF_550_LKT(66,6)=216.450000 !rg=1.89461 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(66,7,1:6)=(/ 177.600000,181.630000,185.450000,191.300000,208.900000,249.870000 /)
XPIZA_LKT(66,7,1:6)=(/ 0.772646,0.797906,0.837233,0.947418,0.970376,0.987633 /)
XCGA_LKT(66,7,1:6)=(/ 0.878217,0.863923,0.839920,0.792300,0.746997,0.747667 /)
XEXT_COEFF_550_LKT(66,7)=185.780000 !rg=1.89461 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(66,8,1:6)=(/ 152.100000,154.300000,157.380000,163.050000,174.030000,204.310000 /)
XPIZA_LKT(66,8,1:6)=(/ 0.770502,0.791992,0.825181,0.940301,0.968012,0.985354 /)
XCGA_LKT(66,8,1:6)=(/ 0.879420,0.867323,0.846547,0.798687,0.763673,0.745970 /)
XEXT_COEFF_550_LKT(66,8)=157.470000 !rg=1.89461 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(66,9,1:6)=(/ 128.190000,130.480000,133.240000,137.040000,146.280000,166.950000 /)
XPIZA_LKT(66,9,1:6)=(/ 0.766519,0.785795,0.813875,0.932075,0.961579,0.981963 /)
XCGA_LKT(66,9,1:6)=(/ 0.880377,0.870473,0.852410,0.805317,0.768107,0.743590 /)
XEXT_COEFF_550_LKT(66,9)=132.910000 !rg=1.89461 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(66,10,1:6)=(/ 108.270000,109.810000,111.710000,115.100000,120.720000,136.180000 /)
XPIZA_LKT(66,10,1:6)=(/ 0.760472,0.780387,0.804925,0.925043,0.956957,0.978797 /)
XCGA_LKT(66,10,1:6)=(/ 0.882517,0.872980,0.858613,0.810733,0.785393,0.750180 /)
XEXT_COEFF_550_LKT(66,10)=111.400000 !rg=1.89461 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(66,11,1:6)=(/ 90.484000,92.040000,93.560000,96.015000,99.915000,111.100000 /)
XPIZA_LKT(66,11,1:6)=(/ 0.752754,0.774224,0.796127,0.917552,0.949789,0.974367 /)
XCGA_LKT(66,11,1:6)=(/ 0.884010,0.876063,0.863183,0.820027,0.794447,0.756240 /)
XEXT_COEFF_550_LKT(66,11)=93.271000 !rg=1.89461 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(66,12,1:6)=(/ 75.743000,76.498000,77.713000,79.563000,82.635000,90.420000 /)
XPIZA_LKT(66,12,1:6)=(/ 0.744121,0.767679,0.788978,0.911679,0.941396,0.969786 /)
XCGA_LKT(66,12,1:6)=(/ 0.886533,0.878247,0.868010,0.825380,0.805080,0.765417 /)
XEXT_COEFF_550_LKT(66,12)=77.692000 !rg=1.89461 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(66,13,1:6)=(/ 62.984000,63.746000,64.485000,66.047000,68.350000,73.947000 /)
XPIZA_LKT(66,13,1:6)=(/ 0.733382,0.760470,0.781151,0.905405,0.933299,0.963913 /)
XCGA_LKT(66,13,1:6)=(/ 0.888530,0.880993,0.871383,0.832030,0.811293,0.771073 /)
XEXT_COEFF_550_LKT(66,13)=64.573000 !rg=1.89461 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(66,14,1:6)=(/ 52.419000,52.951000,53.704000,54.495000,56.501000,60.758000 /)
XPIZA_LKT(66,14,1:6)=(/ 0.721499,0.751305,0.774285,0.900569,0.926260,0.957002 /)
XCGA_LKT(66,14,1:6)=(/ 0.891223,0.883273,0.874927,0.836993,0.820080,0.781060 /)
XEXT_COEFF_550_LKT(66,14)=53.686000 !rg=1.89461 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(66,15,1:6)=(/ 43.663000,43.955000,44.522000,45.333000,46.597000,49.453000 /)
XPIZA_LKT(66,15,1:6)=(/ 0.709572,0.741169,0.765074,0.895866,0.919650,0.949995 /)
XCGA_LKT(66,15,1:6)=(/ 0.893907,0.885553,0.878130,0.841933,0.827717,0.791510 /)
XEXT_COEFF_550_LKT(66,15)=44.569000 !rg=1.89461 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(66,16,1:6)=(/ 36.154000,36.531000,36.864000,37.468000,38.216000,40.348000 /)
XPIZA_LKT(66,16,1:6)=(/ 0.695740,0.730192,0.756264,0.890681,0.913562,0.941633 /)
XCGA_LKT(66,16,1:6)=(/ 0.896640,0.888753,0.880843,0.845877,0.835627,0.801623 /)
XEXT_COEFF_550_LKT(66,16)=36.819000 !rg=1.89461 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(66,17,1:6)=(/ 30.080000,30.269000,30.593000,30.976000,31.703000,33.274000 /)
XPIZA_LKT(66,17,1:6)=(/ 0.682648,0.717932,0.746399,0.885010,0.907762,0.934123 /)
XCGA_LKT(66,17,1:6)=(/ 0.899957,0.891273,0.884000,0.849327,0.841360,0.810080 /)
XEXT_COEFF_550_LKT(66,17)=30.598000 !rg=1.89461 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(66,18,1:6)=(/ 25.002000,25.142000,25.370000,25.712000,26.231000,27.324000 /)
XPIZA_LKT(66,18,1:6)=(/ 0.669283,0.704936,0.734858,0.878983,0.901918,0.927315 /)
XCGA_LKT(66,18,1:6)=(/ 0.903330,0.894157,0.886737,0.852780,0.845457,0.817880 /)
XEXT_COEFF_550_LKT(66,18)=25.401000 !rg=1.89461 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(66,19,1:6)=(/ 20.734000,20.876000,21.010000,21.260000,21.566000,22.371000 /)
XPIZA_LKT(66,19,1:6)=(/ 0.655767,0.691749,0.722827,0.871451,0.896411,0.920748 /)
XCGA_LKT(66,19,1:6)=(/ 0.906737,0.897610,0.889597,0.855813,0.850873,0.825910 /)
XEXT_COEFF_550_LKT(66,19)=21.008000 !rg=1.89461 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(66,20,1:6)=(/ 17.226000,17.305000,17.438000,17.599000,17.888000,18.496000 /)
XPIZA_LKT(66,20,1:6)=(/ 0.643195,0.678131,0.710390,0.863051,0.890354,0.914806 /)
XCGA_LKT(66,20,1:6)=(/ 0.910483,0.900747,0.892763,0.858463,0.854810,0.832623 /)
XEXT_COEFF_550_LKT(66,20)=17.433000 !rg=1.89461 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(67,1,1:6)=(/ 309.090000,314.680000,306.510000,333.110000,307.140000,538.230000 /)
XPIZA_LKT(67,1,1:6)=(/ 0.776053,0.840280,0.878417,0.960778,0.979754,0.994609 /)
XCGA_LKT(67,1,1:6)=(/ 0.875490,0.843333,0.795063,0.769563,0.613143,0.801320 /)
XEXT_COEFF_550_LKT(67,1)=307.240000 !rg=2.05254 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(67,2,1:6)=(/ 297.340000,303.570000,310.730000,321.830000,307.410000,516.970000 /)
XPIZA_LKT(67,2,1:6)=(/ 0.769262,0.837215,0.882481,0.963633,0.979708,0.994310 /)
XCGA_LKT(67,2,1:6)=(/ 0.877640,0.846830,0.814920,0.783167,0.638057,0.797363 /)
XEXT_COEFF_550_LKT(67,2)=310.860000 !rg=2.05254 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(67,3,1:6)=(/ 274.010000,281.160000,288.040000,298.540000,302.990000,468.640000 /)
XPIZA_LKT(67,3,1:6)=(/ 0.768300,0.828829,0.874919,0.962973,0.979488,0.993673 /)
XCGA_LKT(67,3,1:6)=(/ 0.878987,0.849757,0.821453,0.777887,0.676253,0.789733 /)
XEXT_COEFF_550_LKT(67,3)=288.090000 !rg=2.05254 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(67,4,1:6)=(/ 247.100000,253.560000,259.220000,268.490000,283.570000,402.810000 /)
XPIZA_LKT(67,4,1:6)=(/ 0.770346,0.817963,0.866333,0.959797,0.978316,0.992573 /)
XCGA_LKT(67,4,1:6)=(/ 0.878677,0.854073,0.826890,0.778553,0.706130,0.776500 /)
XEXT_COEFF_550_LKT(67,4)=259.680000 !rg=2.05254 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(67,5,1:6)=(/ 218.590000,223.670000,228.920000,236.300000,255.900000,335.900000 /)
XPIZA_LKT(67,5,1:6)=(/ 0.771457,0.807053,0.854571,0.955484,0.975156,0.990985 /)
XCGA_LKT(67,5,1:6)=(/ 0.878710,0.859150,0.831590,0.786493,0.724547,0.762120 /)
XEXT_COEFF_550_LKT(67,5)=229.840000 !rg=2.05254 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(67,6,1:6)=(/ 190.770000,194.570000,199.010000,206.810000,219.690000,276.280000 /)
XPIZA_LKT(67,6,1:6)=(/ 0.771482,0.798449,0.842735,0.948875,0.972961,0.989041 /)
XCGA_LKT(67,6,1:6)=(/ 0.878627,0.863263,0.839047,0.783710,0.743770,0.750810 /)
XEXT_COEFF_550_LKT(67,6)=198.740000 !rg=2.05254 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(67,7,1:6)=(/ 163.730000,167.290000,170.550000,176.030000,187.570000,226.140000 /)
XPIZA_LKT(67,7,1:6)=(/ 0.771111,0.792109,0.829642,0.943255,0.969140,0.986590 /)
XCGA_LKT(67,7,1:6)=(/ 0.879573,0.866603,0.845640,0.794340,0.755617,0.744823 /)
XEXT_COEFF_550_LKT(67,7)=170.780000 !rg=2.05254 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(67,8,1:6)=(/ 140.080000,141.920000,144.970000,149.870000,158.020000,185.190000 /)
XPIZA_LKT(67,8,1:6)=(/ 0.767581,0.786874,0.818454,0.935883,0.964957,0.983838 /)
XCGA_LKT(67,8,1:6)=(/ 0.880823,0.869953,0.851443,0.799670,0.769527,0.743683 /)
XEXT_COEFF_550_LKT(67,8)=144.680000 !rg=2.05254 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(67,9,1:6)=(/ 118.310000,120.200000,122.350000,125.590000,133.820000,152.400000 /)
XPIZA_LKT(67,9,1:6)=(/ 0.763313,0.781641,0.807718,0.928473,0.958833,0.979793 /)
XCGA_LKT(67,9,1:6)=(/ 0.882100,0.872887,0.856853,0.811663,0.776307,0.744627 /)
XEXT_COEFF_550_LKT(67,9)=122.190000 !rg=2.05254 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(67,10,1:6)=(/ 99.600000,100.980000,102.880000,105.630000,111.260000,124.280000 /)
XPIZA_LKT(67,10,1:6)=(/ 0.757413,0.777051,0.799314,0.920918,0.952135,0.976054 /)
XCGA_LKT(67,10,1:6)=(/ 0.883257,0.874990,0.861363,0.816980,0.784857,0.748697 /)
XEXT_COEFF_550_LKT(67,10)=102.830000 !rg=2.05254 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(67,11,1:6)=(/ 83.655000,84.603000,86.057000,87.945000,91.971000,101.410000 /)
XPIZA_LKT(67,11,1:6)=(/ 0.749536,0.770650,0.791241,0.914642,0.945728,0.972099 /)
XCGA_LKT(67,11,1:6)=(/ 0.885390,0.877860,0.866367,0.823397,0.799267,0.758680 /)
XEXT_COEFF_550_LKT(67,11)=85.426000 !rg=2.05254 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(67,12,1:6)=(/ 69.700000,70.554000,71.477000,73.160000,76.124000,82.861000 /)
XPIZA_LKT(67,12,1:6)=(/ 0.739402,0.764625,0.784369,0.908108,0.936956,0.966418 /)
XCGA_LKT(67,12,1:6)=(/ 0.887393,0.879940,0.869847,0.829790,0.805803,0.764640 /)
XEXT_COEFF_550_LKT(67,12)=71.627000 !rg=2.05254 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(67,13,1:6)=(/ 58.114000,58.600000,59.481000,60.530000,62.787000,67.844000 /)
XPIZA_LKT(67,13,1:6)=(/ 0.728529,0.756005,0.777553,0.902924,0.930062,0.960675 /)
XCGA_LKT(67,13,1:6)=(/ 0.889917,0.882060,0.873640,0.834823,0.816550,0.775970 /)
XEXT_COEFF_550_LKT(67,13)=59.537000 !rg=2.05254 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(67,14,1:6)=(/ 48.323000,48.726000,49.342000,50.119000,51.782000,55.286000 /)
XPIZA_LKT(67,14,1:6)=(/ 0.716291,0.746721,0.769318,0.898465,0.923491,0.954653 /)
XCGA_LKT(67,14,1:6)=(/ 0.892270,0.884787,0.876557,0.840183,0.824053,0.785570 /)
XEXT_COEFF_550_LKT(67,14)=49.459000 !rg=2.05254 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(67,15,1:6)=(/ 40.173000,40.552000,40.938000,41.632000,42.739000,45.122000 /)
XPIZA_LKT(67,15,1:6)=(/ 0.703295,0.736561,0.761258,0.893503,0.917192,0.946927 /)
XCGA_LKT(67,15,1:6)=(/ 0.895060,0.887370,0.879180,0.844517,0.832753,0.798000 /)
XEXT_COEFF_550_LKT(67,15)=40.862000 !rg=2.05254 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(67,16,1:6)=(/ 33.395000,33.622000,33.996000,34.398000,35.381000,37.196000 /)
XPIZA_LKT(67,16,1:6)=(/ 0.690348,0.724429,0.751849,0.887819,0.910580,0.938250 /)
XCGA_LKT(67,16,1:6)=(/ 0.898223,0.889987,0.882673,0.847913,0.837833,0.805790 /)
XEXT_COEFF_550_LKT(67,16)=33.851000 !rg=2.05254 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(67,17,1:6)=(/ 27.734000,27.923000,28.127000,28.487000,29.172000,30.518000 /)
XPIZA_LKT(67,17,1:6)=(/ 0.676562,0.712108,0.741344,0.882405,0.905248,0.931329 /)
XCGA_LKT(67,17,1:6)=(/ 0.901497,0.892817,0.885070,0.851723,0.843153,0.812967 /)
XEXT_COEFF_550_LKT(67,17)=28.189000 !rg=2.05254 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(67,18,1:6)=(/ 23.040000,23.202000,23.355000,23.641000,24.117000,25.029000 /)
XPIZA_LKT(67,18,1:6)=(/ 0.663004,0.699412,0.729742,0.875744,0.899558,0.924917 /)
XCGA_LKT(67,18,1:6)=(/ 0.904853,0.895857,0.887853,0.854573,0.848553,0.823313 /)
XEXT_COEFF_550_LKT(67,18)=23.328000 !rg=2.05254 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(67,19,1:6)=(/ 19.138000,19.233000,19.386000,19.564000,19.973000,20.673000 /)
XPIZA_LKT(67,19,1:6)=(/ 0.650089,0.685582,0.717467,0.867710,0.893713,0.917942 /)
XCGA_LKT(67,19,1:6)=(/ 0.908543,0.898970,0.891250,0.857037,0.852360,0.829790 /)
XEXT_COEFF_550_LKT(67,19)=19.363000 !rg=2.05254 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(67,20,1:6)=(/ 15.885000,15.967000,16.051000,16.207000,16.479000,17.001000 /)
XPIZA_LKT(67,20,1:6)=(/ 0.637586,0.672132,0.704296,0.859183,0.887533,0.911987 /)
XCGA_LKT(67,20,1:6)=(/ 0.912197,0.902403,0.894070,0.860240,0.856233,0.834773 /)
XEXT_COEFF_550_LKT(67,20)=16.077000 !rg=2.05254 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(68,1,1:6)=(/ 287.130000,291.980000,305.100000,282.490000,252.930000,525.600000 /)
XPIZA_LKT(68,1,1:6)=(/ 0.758562,0.830638,0.881864,0.951206,0.975016,0.994318 /)
XCGA_LKT(68,1,1:6)=(/ 0.882543,0.848073,0.824857,0.708547,0.571327,0.803120 /)
XEXT_COEFF_550_LKT(68,1)=310.940000 !rg=2.22363 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(68,2,1:6)=(/ 272.610000,279.960000,286.380000,289.970000,268.220000,494.500000 /)
XPIZA_LKT(68,2,1:6)=(/ 0.755722,0.827545,0.874166,0.959169,0.976697,0.993997 /)
XCGA_LKT(68,2,1:6)=(/ 0.883043,0.852100,0.822203,0.767027,0.629413,0.798070 /)
XEXT_COEFF_550_LKT(68,2)=285.980000 !rg=2.22363 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(68,3,1:6)=(/ 251.850000,257.740000,264.770000,275.950000,273.820000,437.500000 /)
XPIZA_LKT(68,3,1:6)=(/ 0.761468,0.818410,0.868051,0.958441,0.977446,0.993168 /)
XCGA_LKT(68,3,1:6)=(/ 0.881853,0.855190,0.826107,0.786730,0.684110,0.786520 /)
XEXT_COEFF_550_LKT(68,3)=263.870000 !rg=2.22363 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(68,4,1:6)=(/ 227.170000,232.490000,238.530000,248.660000,257.060000,368.730000 /)
XPIZA_LKT(68,4,1:6)=(/ 0.766367,0.807468,0.858646,0.955330,0.976390,0.991834 /)
XCGA_LKT(68,4,1:6)=(/ 0.880707,0.859320,0.831123,0.787487,0.715077,0.770537 /)
XEXT_COEFF_550_LKT(68,4)=237.760000 !rg=2.22363 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(68,5,1:6)=(/ 201.300000,205.850000,210.210000,216.690000,229.610000,303.950000 /)
XPIZA_LKT(68,5,1:6)=(/ 0.769623,0.799307,0.846336,0.952541,0.973983,0.990034 /)
XCGA_LKT(68,5,1:6)=(/ 0.880250,0.862977,0.837683,0.789023,0.733797,0.755413 /)
XEXT_COEFF_550_LKT(68,5)=210.430000 !rg=2.22363 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(68,6,1:6)=(/ 175.950000,179.050000,183.330000,189.700000,200.100000,249.810000 /)
XPIZA_LKT(68,6,1:6)=(/ 0.770070,0.792683,0.833533,0.945600,0.971057,0.987871 /)
XCGA_LKT(68,6,1:6)=(/ 0.880217,0.866743,0.842697,0.794083,0.748270,0.745987 /)
XEXT_COEFF_550_LKT(68,6)=182.610000 !rg=2.22363 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(68,7,1:6)=(/ 150.580000,153.670000,157.040000,162.230000,170.820000,204.730000 /)
XPIZA_LKT(68,7,1:6)=(/ 0.768407,0.786350,0.821839,0.938778,0.966934,0.985233 /)
XCGA_LKT(68,7,1:6)=(/ 0.880537,0.869780,0.849137,0.801333,0.762607,0.742367 /)
XEXT_COEFF_550_LKT(68,7)=156.420000 !rg=2.22363 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(68,8,1:6)=(/ 128.750000,130.860000,133.290000,137.700000,144.800000,167.940000 /)
XPIZA_LKT(68,8,1:6)=(/ 0.765165,0.783047,0.811259,0.931268,0.962057,0.982122 /)
XCGA_LKT(68,8,1:6)=(/ 0.881747,0.872253,0.854883,0.806810,0.770857,0.741610 /)
XEXT_COEFF_550_LKT(68,8)=133.330000 !rg=2.22363 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(68,9,1:6)=(/ 109.100000,110.550000,112.590000,115.510000,121.130000,137.270000 /)
XPIZA_LKT(68,9,1:6)=(/ 0.760348,0.778688,0.801798,0.924760,0.956857,0.978738 /)
XCGA_LKT(68,9,1:6)=(/ 0.883343,0.874060,0.860607,0.814577,0.783880,0.747130 /)
XEXT_COEFF_550_LKT(68,9)=112.450000 !rg=2.22363 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(68,10,1:6)=(/ 91.749000,93.175000,94.513000,96.782000,101.890000,113.510000 /)
XPIZA_LKT(68,10,1:6)=(/ 0.752744,0.772860,0.794645,0.916906,0.948994,0.973818 /)
XCGA_LKT(68,10,1:6)=(/ 0.884503,0.876997,0.864667,0.821083,0.792613,0.752537 /)
XEXT_COEFF_550_LKT(68,10)=94.658000 !rg=2.22363 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(68,11,1:6)=(/ 76.888000,77.965000,79.085000,80.795000,84.664000,92.644000 /)
XPIZA_LKT(68,11,1:6)=(/ 0.744771,0.767911,0.787131,0.910544,0.941494,0.969515 /)
XCGA_LKT(68,11,1:6)=(/ 0.886133,0.879320,0.868987,0.827650,0.800193,0.758453 /)
XEXT_COEFF_550_LKT(68,11)=78.864000 !rg=2.22363 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(68,12,1:6)=(/ 64.234000,64.997000,65.876000,67.045000,69.815000,75.881000 /)
XPIZA_LKT(68,12,1:6)=(/ 0.734102,0.760200,0.780846,0.905225,0.933672,0.963754 /)
XCGA_LKT(68,12,1:6)=(/ 0.888690,0.881403,0.872313,0.832497,0.811823,0.770223 /)
XEXT_COEFF_550_LKT(68,12)=65.930000 !rg=2.22363 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(68,13,1:6)=(/ 53.522000,54.047000,54.635000,55.505000,57.556000,61.719000 /)
XPIZA_LKT(68,13,1:6)=(/ 0.723071,0.752219,0.773315,0.900384,0.927201,0.958669 /)
XCGA_LKT(68,13,1:6)=(/ 0.891077,0.883790,0.875160,0.838367,0.820307,0.779930 /)
XEXT_COEFF_550_LKT(68,13)=54.762000 !rg=2.22363 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(68,14,1:6)=(/ 44.537000,44.959000,45.371000,46.164000,47.610000,50.561000 /)
XPIZA_LKT(68,14,1:6)=(/ 0.710360,0.742490,0.765495,0.895502,0.919522,0.950719 /)
XCGA_LKT(68,14,1:6)=(/ 0.893587,0.886213,0.878490,0.842830,0.828430,0.792270 /)
XEXT_COEFF_550_LKT(68,14)=45.381000 !rg=2.22363 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(68,15,1:6)=(/ 37.078000,37.353000,37.727000,38.338000,39.256000,41.419000 /)
XPIZA_LKT(68,15,1:6)=(/ 0.697547,0.731291,0.756930,0.890697,0.913297,0.942343 /)
XCGA_LKT(68,15,1:6)=(/ 0.896580,0.888333,0.881270,0.845483,0.835763,0.801867 /)
XEXT_COEFF_550_LKT(68,15)=37.677000 !rg=2.22363 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(68,16,1:6)=(/ 30.768000,31.026000,31.271000,31.724000,32.538000,34.167000 /)
XPIZA_LKT(68,16,1:6)=(/ 0.683768,0.719345,0.746850,0.885139,0.907499,0.934652 /)
XCGA_LKT(68,16,1:6)=(/ 0.899707,0.891460,0.883983,0.849963,0.839267,0.806797 /)
XEXT_COEFF_550_LKT(68,16)=31.234000 !rg=2.22363 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(68,17,1:6)=(/ 25.571000,25.763000,25.946000,26.197000,26.828000,28.005000 /)
XPIZA_LKT(68,17,1:6)=(/ 0.670213,0.706533,0.736082,0.879000,0.902313,0.927906 /)
XCGA_LKT(68,17,1:6)=(/ 0.902990,0.894360,0.886780,0.853050,0.845917,0.818060 /)
XEXT_COEFF_550_LKT(68,17)=25.940000 !rg=2.22363 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(68,18,1:6)=(/ 21.265000,21.376000,21.545000,21.775000,22.195000,22.990000 /)
XPIZA_LKT(68,18,1:6)=(/ 0.657216,0.693183,0.724667,0.872299,0.896768,0.920735 /)
XCGA_LKT(68,18,1:6)=(/ 0.906570,0.897057,0.889633,0.855350,0.851057,0.826077 /)
XEXT_COEFF_550_LKT(68,18)=21.526000 !rg=2.22363 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(68,19,1:6)=(/ 17.645000,17.751000,17.850000,18.065000,18.354000,19.014000 /)
XPIZA_LKT(68,19,1:6)=(/ 0.644194,0.679791,0.711619,0.864108,0.890528,0.915160 /)
XCGA_LKT(68,19,1:6)=(/ 0.910173,0.900630,0.892507,0.858840,0.853083,0.830587 /)
XEXT_COEFF_550_LKT(68,19)=17.842000 !rg=2.22363 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(68,20,1:6)=(/ 14.655000,14.732000,14.809000,14.902000,15.157000,15.623000 /)
XPIZA_LKT(68,20,1:6)=(/ 0.632115,0.666311,0.698444,0.854307,0.884299,0.909295 /)
XCGA_LKT(68,20,1:6)=(/ 0.913867,0.904087,0.895710,0.861167,0.857993,0.838507 /)
XEXT_COEFF_550_LKT(68,20)=14.809000 !rg=2.22363 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(69,1,1:6)=(/ 262.030000,268.390000,282.790000,260.230000,230.620000,498.400000 /)
XPIZA_LKT(69,1,1:6)=(/ 0.737229,0.822873,0.874558,0.959473,0.972628,0.994075 /)
XCGA_LKT(69,1,1:6)=(/ 0.889493,0.856203,0.834623,0.791457,0.586847,0.801397 /)
XEXT_COEFF_550_LKT(69,1)=285.290000 !rg=2.40898 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(69,2,1:6)=(/ 252.090000,257.110000,265.150000,271.700000,251.850000,461.650000 /)
XPIZA_LKT(69,2,1:6)=(/ 0.749565,0.817612,0.865127,0.955986,0.975134,0.993518 /)
XCGA_LKT(69,2,1:6)=(/ 0.886030,0.856103,0.828473,0.785490,0.653370,0.794850 /)
XEXT_COEFF_550_LKT(69,2)=261.290000 !rg=2.40898 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(69,3,1:6)=(/ 232.220000,236.910000,244.260000,253.300000,253.990000,399.480000 /)
XPIZA_LKT(69,3,1:6)=(/ 0.759261,0.807204,0.858966,0.956078,0.975460,0.992456 /)
XCGA_LKT(69,3,1:6)=(/ 0.883093,0.860887,0.832407,0.777840,0.698710,0.779373 /)
XEXT_COEFF_550_LKT(69,3)=241.260000 !rg=2.40898 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(69,4,1:6)=(/ 209.590000,213.860000,219.620000,227.640000,235.590000,332.170000 /)
XPIZA_LKT(69,4,1:6)=(/ 0.765276,0.797591,0.849312,0.953404,0.973988,0.990853 /)
XCGA_LKT(69,4,1:6)=(/ 0.881497,0.864657,0.837203,0.784570,0.725180,0.761150 /)
XEXT_COEFF_550_LKT(69,4)=217.440000 !rg=2.40898 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(69,5,1:6)=(/ 185.140000,188.760000,193.350000,200.620000,209.240000,272.970000 /)
XPIZA_LKT(69,5,1:6)=(/ 0.768059,0.790705,0.838108,0.947573,0.972106,0.988878 /)
XCGA_LKT(69,5,1:6)=(/ 0.881173,0.867417,0.842053,0.796030,0.744937,0.748253 /)
XEXT_COEFF_550_LKT(69,5)=192.540000 !rg=2.40898 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(69,6,1:6)=(/ 161.710000,164.070000,168.190000,172.870000,183.470000,225.250000 /)
XPIZA_LKT(69,6,1:6)=(/ 0.768906,0.786842,0.826005,0.942674,0.967875,0.986349 /)
XCGA_LKT(69,6,1:6)=(/ 0.880937,0.869747,0.848550,0.798490,0.756543,0.740257 /)
XEXT_COEFF_550_LKT(69,6)=168.140000 !rg=2.40898 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(69,7,1:6)=(/ 139.100000,141.280000,144.260000,148.230000,156.330000,184.970000 /)
XPIZA_LKT(69,7,1:6)=(/ 0.767101,0.782464,0.813745,0.935871,0.964111,0.983510 /)
XCGA_LKT(69,7,1:6)=(/ 0.881473,0.872740,0.854133,0.804313,0.769240,0.739057 /)
XEXT_COEFF_550_LKT(69,7)=143.270000 !rg=2.40898 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(69,8,1:6)=(/ 118.610000,120.220000,122.680000,125.830000,133.080000,152.640000 /)
XPIZA_LKT(69,8,1:6)=(/ 0.762103,0.779124,0.804919,0.927102,0.958580,0.980022 /)
XCGA_LKT(69,8,1:6)=(/ 0.882900,0.874010,0.859037,0.811113,0.778443,0.741860 /)
XEXT_COEFF_550_LKT(69,8)=122.920000 !rg=2.40898 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(69,9,1:6)=(/ 100.260000,101.620000,103.540000,106.220000,110.440000,124.450000 /)
XPIZA_LKT(69,9,1:6)=(/ 0.756935,0.774600,0.797266,0.920247,0.954268,0.977055 /)
XCGA_LKT(69,9,1:6)=(/ 0.883880,0.876553,0.863133,0.819243,0.792027,0.750263 /)
XEXT_COEFF_550_LKT(69,9)=103.200000 !rg=2.40898 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(69,10,1:6)=(/ 84.661000,85.564000,87.105000,88.933000,92.774000,102.420000 /)
XPIZA_LKT(69,10,1:6)=(/ 0.749703,0.770377,0.789120,0.914237,0.946114,0.972670 /)
XCGA_LKT(69,10,1:6)=(/ 0.885433,0.878283,0.867600,0.824957,0.798597,0.756353 /)
XEXT_COEFF_550_LKT(69,10)=87.037000 !rg=2.40898 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(69,11,1:6)=(/ 70.810000,71.763000,72.795000,74.365000,77.610000,84.743000 /)
XPIZA_LKT(69,11,1:6)=(/ 0.740041,0.763871,0.783655,0.908187,0.937354,0.966103 /)
XCGA_LKT(69,11,1:6)=(/ 0.887533,0.880843,0.871430,0.832060,0.806410,0.763617 /)
XEXT_COEFF_550_LKT(69,11)=72.536000 !rg=2.40898 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(69,12,1:6)=(/ 59.243000,59.736000,60.626000,61.617000,63.863000,68.808000 /)
XPIZA_LKT(69,12,1:6)=(/ 0.729556,0.756317,0.776393,0.903057,0.930865,0.962083 /)
XCGA_LKT(69,12,1:6)=(/ 0.889707,0.882710,0.874297,0.835990,0.816563,0.774923 /)
XEXT_COEFF_550_LKT(69,12)=60.693000 !rg=2.40898 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(69,13,1:6)=(/ 49.328000,49.837000,50.332000,51.137000,52.890000,56.438000 /)
XPIZA_LKT(69,13,1:6)=(/ 0.717278,0.747809,0.769578,0.897573,0.923223,0.955043 /)
XCGA_LKT(69,13,1:6)=(/ 0.892210,0.885307,0.877510,0.841020,0.825213,0.786877 /)
XEXT_COEFF_550_LKT(69,13)=50.358000 !rg=2.40898 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(69,14,1:6)=(/ 41.130000,41.384000,41.854000,42.461000,43.577000,46.112000 /)
XPIZA_LKT(69,14,1:6)=(/ 0.704800,0.737124,0.761890,0.893438,0.916668,0.946705 /)
XCGA_LKT(69,14,1:6)=(/ 0.895210,0.887240,0.880207,0.844440,0.832863,0.797113 /)
XEXT_COEFF_550_LKT(69,14)=41.856000 !rg=2.40898 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(69,15,1:6)=(/ 34.162000,34.414000,34.761000,35.266000,36.279000,38.228000 /)
XPIZA_LKT(69,15,1:6)=(/ 0.691351,0.726058,0.752569,0.888234,0.909828,0.937602 /)
XCGA_LKT(69,15,1:6)=(/ 0.897863,0.889770,0.882407,0.848533,0.836457,0.801653 /)
XEXT_COEFF_550_LKT(69,15)=34.776000 !rg=2.40898 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(69,16,1:6)=(/ 28.343000,28.569000,28.845000,29.220000,29.934000,31.382000 /)
XPIZA_LKT(69,16,1:6)=(/ 0.677314,0.713416,0.742235,0.882839,0.904565,0.930530 /)
XCGA_LKT(69,16,1:6)=(/ 0.901223,0.892787,0.885623,0.852003,0.842227,0.811623 /)
XEXT_COEFF_550_LKT(69,16)=28.771000 !rg=2.40898 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(69,17,1:6)=(/ 23.598000,23.711000,23.937000,24.197000,24.659000,25.584000 /)
XPIZA_LKT(69,17,1:6)=(/ 0.664387,0.700283,0.731028,0.876076,0.899134,0.924326 /)
XCGA_LKT(69,17,1:6)=(/ 0.904660,0.895597,0.888440,0.854427,0.848643,0.822793 /)
XEXT_COEFF_550_LKT(69,17)=23.926000 !rg=2.40898 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(69,18,1:6)=(/ 19.600000,19.707000,19.848000,20.063000,20.443000,21.238000 /)
XPIZA_LKT(69,18,1:6)=(/ 0.651125,0.687155,0.718851,0.868857,0.893390,0.917264 /)
XCGA_LKT(69,18,1:6)=(/ 0.908113,0.898643,0.890727,0.857607,0.851853,0.826740 /)
XEXT_COEFF_550_LKT(69,18)=19.861000 !rg=2.40898 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(69,19,1:6)=(/ 16.266000,16.359000,16.472000,16.615000,16.901000,17.492000 /)
XPIZA_LKT(69,19,1:6)=(/ 0.638414,0.673587,0.706020,0.860212,0.887463,0.911773 /)
XCGA_LKT(69,19,1:6)=(/ 0.911850,0.902120,0.894157,0.860237,0.855383,0.833817 /)
XEXT_COEFF_550_LKT(69,19)=16.453000 !rg=2.40898 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(69,20,1:6)=(/ 13.523000,13.578000,13.663000,13.758000,13.972000,14.327000 /)
XPIZA_LKT(69,20,1:6)=(/ 0.626927,0.660248,0.692669,0.849752,0.880863,0.906244 /)
XCGA_LKT(69,20,1:6)=(/ 0.915647,0.905580,0.897343,0.862497,0.859850,0.841953 /)
XEXT_COEFF_550_LKT(69,20)=13.661000 !rg=2.40898 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(70,1,1:6)=(/ 240.290000,247.860000,243.880000,277.750000,236.250000,460.350000 /)
XPIZA_LKT(70,1,1:6)=(/ 0.732652,0.811541,0.853904,0.961644,0.973474,0.993435 /)
XCGA_LKT(70,1,1:6)=(/ 0.892207,0.861500,0.824557,0.810337,0.654830,0.798520 /)
XEXT_COEFF_550_LKT(70,1)=248.270000 !rg=2.60978 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(70,2,1:6)=(/ 231.970000,234.990000,241.580000,250.750000,248.700000,419.480000 /)
XPIZA_LKT(70,2,1:6)=(/ 0.750879,0.806686,0.859553,0.956961,0.974636,0.992814 /)
XCGA_LKT(70,2,1:6)=(/ 0.886673,0.861623,0.836547,0.789553,0.693553,0.787093 /)
XEXT_COEFF_550_LKT(70,2)=241.170000 !rg=2.60978 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(70,3,1:6)=(/ 213.310000,218.060000,222.990000,231.140000,238.240000,356.270000 /)
XPIZA_LKT(70,3,1:6)=(/ 0.761197,0.796780,0.851768,0.952585,0.973973,0.991508 /)
XCGA_LKT(70,3,1:6)=(/ 0.883367,0.866263,0.838660,0.791723,0.718397,0.768760 /)
XEXT_COEFF_550_LKT(70,3)=223.110000 !rg=2.60978 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(70,4,1:6)=(/ 192.620000,196.650000,200.820000,207.780000,217.320000,294.460000 /)
XPIZA_LKT(70,4,1:6)=(/ 0.766209,0.789095,0.841127,0.949378,0.971766,0.989674 /)
XCGA_LKT(70,4,1:6)=(/ 0.881783,0.869090,0.842773,0.795523,0.733657,0.750843 /)
XEXT_COEFF_550_LKT(70,4)=200.980000 !rg=2.60978 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(70,5,1:6)=(/ 170.710000,173.940000,177.690000,183.400000,191.600000,243.490000 /)
XPIZA_LKT(70,5,1:6)=(/ 0.767450,0.784010,0.828461,0.945476,0.969421,0.987409 /)
XCGA_LKT(70,5,1:6)=(/ 0.881730,0.871277,0.847740,0.796680,0.754133,0.739027 /)
XEXT_COEFF_550_LKT(70,5)=176.300000 !rg=2.60978 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(70,6,1:6)=(/ 148.960000,151.040000,154.000000,157.630000,167.890000,201.470000 /)
XPIZA_LKT(70,6,1:6)=(/ 0.767371,0.780908,0.817416,0.939476,0.965562,0.984833 /)
XCGA_LKT(70,6,1:6)=(/ 0.881860,0.873297,0.852490,0.806903,0.759237,0.735617 /)
XEXT_COEFF_550_LKT(70,6)=154.320000 !rg=2.60978 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(70,7,1:6)=(/ 127.860000,130.100000,132.420000,135.860000,143.640000,166.530000 /)
XPIZA_LKT(70,7,1:6)=(/ 0.764821,0.779157,0.806962,0.930823,0.960887,0.981802 /)
XCGA_LKT(70,7,1:6)=(/ 0.881977,0.874870,0.858277,0.810993,0.770027,0.736417 /)
XEXT_COEFF_550_LKT(70,7)=132.160000 !rg=2.60978 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(70,8,1:6)=(/ 109.220000,110.500000,112.580000,115.100000,121.100000,137.430000 /)
XPIZA_LKT(70,8,1:6)=(/ 0.760168,0.775320,0.798693,0.924070,0.956894,0.978686 /)
XCGA_LKT(70,8,1:6)=(/ 0.883353,0.876577,0.862033,0.817253,0.784100,0.742707 /)
XEXT_COEFF_550_LKT(70,8)=112.860000 !rg=2.60978 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(70,9,1:6)=(/ 92.524000,93.707000,95.092000,97.538000,101.370000,113.130000 /)
XPIZA_LKT(70,9,1:6)=(/ 0.753773,0.771529,0.791385,0.916829,0.949885,0.974622 /)
XCGA_LKT(70,9,1:6)=(/ 0.884967,0.878213,0.866860,0.822010,0.797330,0.752140 /)
XEXT_COEFF_550_LKT(70,9)=94.695000 !rg=2.60978 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(70,10,1:6)=(/ 77.888000,78.888000,79.950000,81.678000,84.895000,93.035000 /)
XPIZA_LKT(70,10,1:6)=(/ 0.745126,0.767351,0.785985,0.910598,0.942387,0.970607 /)
XCGA_LKT(70,10,1:6)=(/ 0.886510,0.879980,0.869957,0.829467,0.806157,0.762630 /)
XEXT_COEFF_550_LKT(70,10)=79.830000 !rg=2.60978 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(70,11,1:6)=(/ 65.346000,66.116000,66.986000,68.491000,70.621000,76.452000 /)
XPIZA_LKT(70,11,1:6)=(/ 0.735680,0.760579,0.779300,0.904717,0.934467,0.964594 /)
XCGA_LKT(70,11,1:6)=(/ 0.888777,0.881697,0.873733,0.833870,0.813270,0.770247 /)
XEXT_COEFF_550_LKT(70,11)=66.837000 !rg=2.60978 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(70,12,1:6)=(/ 54.572000,55.133000,55.683000,56.672000,58.702000,62.789000 /)
XPIZA_LKT(70,12,1:6)=(/ 0.724006,0.752888,0.772873,0.899939,0.926501,0.958928 /)
XCGA_LKT(70,12,1:6)=(/ 0.890880,0.884050,0.876267,0.839313,0.821883,0.781970 /)
XEXT_COEFF_550_LKT(70,12)=55.683000 !rg=2.60978 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(70,13,1:6)=(/ 45.537000,45.829000,46.408000,47.095000,48.363000,51.355000 /)
XPIZA_LKT(70,13,1:6)=(/ 0.711941,0.742688,0.765937,0.895369,0.919805,0.951258 /)
XCGA_LKT(70,13,1:6)=(/ 0.893600,0.886173,0.879323,0.842853,0.829930,0.792053 /)
XEXT_COEFF_550_LKT(70,13)=46.366000 !rg=2.60978 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(70,14,1:6)=(/ 37.863000,38.224000,38.505000,39.165000,40.181000,42.385000 /)
XPIZA_LKT(70,14,1:6)=(/ 0.698467,0.732643,0.757052,0.890486,0.913199,0.942568 /)
XCGA_LKT(70,14,1:6)=(/ 0.896320,0.888647,0.881363,0.846810,0.834747,0.798220 /)
XEXT_COEFF_550_LKT(70,14)=38.607000 !rg=2.60978 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(70,15,1:6)=(/ 31.495000,31.760000,31.997000,32.381000,33.225000,34.945000 /)
XPIZA_LKT(70,15,1:6)=(/ 0.684983,0.720417,0.748088,0.885279,0.907421,0.934834 /)
XCGA_LKT(70,15,1:6)=(/ 0.899283,0.891180,0.883943,0.849797,0.840080,0.807750 /)
XEXT_COEFF_550_LKT(70,15)=31.989000 !rg=2.60978 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(70,16,1:6)=(/ 26.171000,26.333000,26.548000,26.876000,27.441000,28.571000 /)
XPIZA_LKT(70,16,1:6)=(/ 0.671648,0.707487,0.737131,0.879669,0.902258,0.927800 /)
XCGA_LKT(70,16,1:6)=(/ 0.902817,0.893947,0.886907,0.853247,0.846050,0.817693 /)
XEXT_COEFF_550_LKT(70,16)=26.551000 !rg=2.60978 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(70,17,1:6)=(/ 21.743000,21.881000,22.018000,22.292000,22.646000,23.490000 /)
XPIZA_LKT(70,17,1:6)=(/ 0.657927,0.694620,0.725385,0.873212,0.896506,0.921463 /)
XCGA_LKT(70,17,1:6)=(/ 0.906243,0.897080,0.889477,0.856047,0.850250,0.824957 /)
XEXT_COEFF_550_LKT(70,17)=22.003000 !rg=2.60978 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(70,18,1:6)=(/ 18.074000,18.180000,18.287000,18.441000,18.770000,19.438000 /)
XPIZA_LKT(70,18,1:6)=(/ 0.645099,0.681049,0.713126,0.864857,0.890741,0.914929 /)
XCGA_LKT(70,18,1:6)=(/ 0.909770,0.900210,0.892297,0.858537,0.853827,0.830990 /)
XEXT_COEFF_550_LKT(70,18)=18.277000 !rg=2.60978 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(70,19,1:6)=(/ 15.013000,15.084000,15.173000,15.295000,15.544000,15.992000 /)
XPIZA_LKT(70,19,1:6)=(/ 0.633133,0.667465,0.699989,0.855658,0.884709,0.909142 /)
XCGA_LKT(70,19,1:6)=(/ 0.913573,0.903663,0.895507,0.861400,0.858017,0.838243 /)
XEXT_COEFF_550_LKT(70,19)=15.174000 !rg=2.60978 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(70,20,1:6)=(/ 12.469000,12.529000,12.586000,12.695000,12.842000,13.184000 /)
XPIZA_LKT(70,20,1:6)=(/ 0.621587,0.654547,0.686522,0.845352,0.877480,0.903590 /)
XCGA_LKT(70,20,1:6)=(/ 0.917323,0.907317,0.898657,0.863930,0.860810,0.843393 /)
XEXT_COEFF_550_LKT(70,20)=12.574000 !rg=2.60978 sigma=2.95 wvl=0.55
 
 IF (LHOOK) CALL DR_HOOK('MODE_DUSTOPT:DUST_OPT_LKT_SET7',1,ZHOOK_HANDLE)
 END SUBROUTINE DUST_OPT_LKT_SET7

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE DUST_OPT_LKT_SET8()

  USE MODD_DUST_OPT_LKT
  
  IMPLICIT NONE


REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_DUSTOPT:DUST_OPT_LKT_SET8',0,ZHOOK_HANDLE)
XEXT_COEFF_WVL_LKT(71,1,1:6)=(/ 220.270000,228.550000,230.480000,262.220000,252.450000,411.740000 /)
XPIZA_LKT(71,1,1:6)=(/ 0.751647,0.801675,0.854372,0.958853,0.975633,0.992690 /)
XCGA_LKT(71,1,1:6)=(/ 0.886540,0.862543,0.839427,0.805113,0.728693,0.787477 /)
XEXT_COEFF_550_LKT(71,1)=233.860000 !rg=2.82732 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(71,2,1:6)=(/ 213.230000,216.050000,221.200000,229.780000,245.280000,369.940000 /)
XPIZA_LKT(71,2,1:6)=(/ 0.755432,0.793764,0.848813,0.954824,0.974547,0.991804 /)
XCGA_LKT(71,2,1:6)=(/ 0.885607,0.867937,0.839433,0.794993,0.734320,0.773827 /)
XEXT_COEFF_550_LKT(71,2)=223.350000 !rg=2.82732 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(71,3,1:6)=(/ 196.900000,200.990000,204.750000,211.060000,225.330000,311.690000 /)
XPIZA_LKT(71,3,1:6)=(/ 0.762618,0.785511,0.841506,0.950406,0.971109,0.990203 /)
XCGA_LKT(71,3,1:6)=(/ 0.883903,0.870187,0.842573,0.798280,0.739073,0.753347 /)
XEXT_COEFF_550_LKT(71,3)=206.090000 !rg=2.82732 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(71,4,1:6)=(/ 177.740000,181.240000,184.580000,189.720000,202.010000,259.160000 /)
XPIZA_LKT(71,4,1:6)=(/ 0.766595,0.780815,0.830812,0.946662,0.968354,0.988142 /)
XCGA_LKT(71,4,1:6)=(/ 0.882733,0.872510,0.847230,0.802463,0.748383,0.737993 /)
XEXT_COEFF_550_LKT(71,4)=185.290000 !rg=2.82732 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(71,5,1:6)=(/ 156.920000,159.930000,163.060000,168.090000,176.760000,215.790000 /)
XPIZA_LKT(71,5,1:6)=(/ 0.767865,0.779082,0.819566,0.940427,0.966071,0.985827 /)
XCGA_LKT(71,5,1:6)=(/ 0.881767,0.874307,0.852913,0.803530,0.756397,0.731880 /)
XEXT_COEFF_550_LKT(71,5)=162.880000 !rg=2.82732 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(71,6,1:6)=(/ 137.310000,139.300000,141.440000,145.890000,153.790000,180.670000 /)
XPIZA_LKT(71,6,1:6)=(/ 0.766004,0.777321,0.808753,0.934652,0.962557,0.983026 /)
XCGA_LKT(71,6,1:6)=(/ 0.882157,0.875670,0.857877,0.810867,0.771793,0.733673 /)
XEXT_COEFF_550_LKT(71,6)=141.880000 !rg=2.82732 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(71,7,1:6)=(/ 117.820000,119.610000,121.780000,124.760000,131.970000,150.860000 /)
XPIZA_LKT(71,7,1:6)=(/ 0.762643,0.775372,0.800467,0.927249,0.957201,0.979212 /)
XCGA_LKT(71,7,1:6)=(/ 0.883263,0.876750,0.862133,0.817380,0.779480,0.736747 /)
XEXT_COEFF_550_LKT(71,7)=121.500000 !rg=2.82732 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(71,8,1:6)=(/ 100.580000,101.870000,103.370000,106.130000,110.820000,124.290000 /)
XPIZA_LKT(71,8,1:6)=(/ 0.756956,0.773213,0.793031,0.919778,0.953218,0.976754 /)
XCGA_LKT(71,8,1:6)=(/ 0.884200,0.878207,0.866173,0.822093,0.793283,0.747647 /)
XEXT_COEFF_550_LKT(71,8)=103.530000 !rg=2.82732 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(71,9,1:6)=(/ 85.088000,86.170000,87.624000,89.562000,93.765000,103.110000 /)
XPIZA_LKT(71,9,1:6)=(/ 0.749987,0.769082,0.787009,0.912831,0.944355,0.971621 /)
XCGA_LKT(71,9,1:6)=(/ 0.885580,0.879417,0.869540,0.826830,0.796580,0.751270 /)
XEXT_COEFF_550_LKT(71,9)=87.478000 !rg=2.82732 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(71,10,1:6)=(/ 71.887000,72.595000,73.632000,75.236000,77.699000,84.617000 /)
XPIZA_LKT(71,10,1:6)=(/ 0.740957,0.763561,0.782207,0.907301,0.938206,0.967817 /)
XCGA_LKT(71,10,1:6)=(/ 0.887810,0.880883,0.872623,0.831383,0.811617,0.766393 /)
XEXT_COEFF_550_LKT(71,10)=73.638000 !rg=2.82732 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(71,11,1:6)=(/ 60.163000,60.994000,61.711000,62.885000,64.529000,69.636000 /)
XPIZA_LKT(71,11,1:6)=(/ 0.730087,0.756903,0.776128,0.901840,0.931307,0.962204 /)
XCGA_LKT(71,11,1:6)=(/ 0.889643,0.883150,0.875107,0.836883,0.819517,0.775647 /)
XEXT_COEFF_550_LKT(71,11)=61.536000 !rg=2.82732 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(71,12,1:6)=(/ 50.388000,50.706000,51.354000,52.193000,53.673000,57.228000 /)
XPIZA_LKT(71,12,1:6)=(/ 0.718826,0.747980,0.769766,0.897553,0.923277,0.955390 /)
XCGA_LKT(71,12,1:6)=(/ 0.892390,0.884953,0.878060,0.841020,0.826717,0.786860 /)
XEXT_COEFF_550_LKT(71,12)=51.400000 !rg=2.82732 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(71,13,1:6)=(/ 41.911000,42.287000,42.663000,43.427000,44.492000,47.100000 /)
XPIZA_LKT(71,13,1:6)=(/ 0.705613,0.738230,0.761440,0.892963,0.916540,0.947391 /)
XCGA_LKT(71,13,1:6)=(/ 0.894820,0.887447,0.880397,0.845087,0.831853,0.793867 /)
XEXT_COEFF_550_LKT(71,13)=42.694000 !rg=2.82732 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(71,14,1:6)=(/ 34.944000,35.168000,35.562000,35.981000,36.863000,38.895000 /)
XPIZA_LKT(71,14,1:6)=(/ 0.692599,0.726922,0.753431,0.888063,0.910575,0.939137 /)
XCGA_LKT(71,14,1:6)=(/ 0.897847,0.889697,0.882960,0.848333,0.838460,0.803840 /)
XEXT_COEFF_550_LKT(71,14)=35.547000 !rg=2.82732 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(71,15,1:6)=(/ 29.069000,29.248000,29.530000,29.939000,30.560000,31.911000 /)
XPIZA_LKT(71,15,1:6)=(/ 0.679029,0.714758,0.743021,0.883018,0.904638,0.931980 /)
XCGA_LKT(71,15,1:6)=(/ 0.900933,0.892370,0.885443,0.851693,0.843490,0.813050 /)
XEXT_COEFF_550_LKT(71,15)=29.542000 !rg=2.82732 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(71,16,1:6)=(/ 24.113000,24.307000,24.467000,24.765000,25.125000,26.110000 /)
XPIZA_LKT(71,16,1:6)=(/ 0.665026,0.701915,0.731827,0.876801,0.899573,0.924999 /)
XCGA_LKT(71,16,1:6)=(/ 0.904320,0.895517,0.887997,0.854503,0.849497,0.822357 /)
XEXT_COEFF_550_LKT(71,16)=24.437000 !rg=2.82732 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(71,17,1:6)=(/ 20.061000,20.166000,20.318000,20.492000,20.884000,21.631000 /)
XPIZA_LKT(71,17,1:6)=(/ 0.652076,0.688348,0.720241,0.869541,0.893800,0.918488 /)
XCGA_LKT(71,17,1:6)=(/ 0.907927,0.898443,0.890980,0.857177,0.853287,0.829553 /)
XEXT_COEFF_550_LKT(71,17)=20.309000 !rg=2.82732 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(71,18,1:6)=(/ 16.681000,16.754000,16.870000,17.035000,17.317000,17.841000 /)
XPIZA_LKT(71,18,1:6)=(/ 0.639578,0.674887,0.707178,0.861191,0.887811,0.912187 /)
XCGA_LKT(71,18,1:6)=(/ 0.911547,0.901690,0.893793,0.859923,0.856037,0.834767 /)
XEXT_COEFF_550_LKT(71,18)=16.884000 !rg=2.82732 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(71,19,1:6)=(/ 13.843000,13.919000,13.988000,14.112000,14.263000,14.644000 /)
XPIZA_LKT(71,19,1:6)=(/ 0.627531,0.661641,0.693903,0.851437,0.881527,0.906733 /)
XCGA_LKT(71,19,1:6)=(/ 0.915263,0.905337,0.896893,0.862340,0.860187,0.841760 /)
XEXT_COEFF_550_LKT(71,19)=13.977000 !rg=2.82732 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(71,20,1:6)=(/ 11.505000,11.548000,11.616000,11.680000,11.840000,12.141000 /)
XPIZA_LKT(71,20,1:6)=(/ 0.616698,0.648563,0.680753,0.839954,0.874046,0.901068 /)
XCGA_LKT(71,20,1:6)=(/ 0.919070,0.908857,0.900340,0.864843,0.863040,0.846953 /)
XEXT_COEFF_550_LKT(71,20)=11.602000 !rg=2.82732 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(72,1,1:6)=(/ 203.440000,208.990000,217.190000,215.450000,262.210000,356.280000 /)
XPIZA_LKT(72,1,1:6)=(/ 0.766075,0.789502,0.846820,0.952100,0.976658,0.991421 /)
XCGA_LKT(72,1,1:6)=(/ 0.883193,0.869100,0.844530,0.788467,0.778473,0.769180 /)
XEXT_COEFF_550_LKT(72,1)=213.360000 !rg=3.06299 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(72,2,1:6)=(/ 195.970000,199.800000,203.130000,210.680000,237.030000,316.700000 /)
XPIZA_LKT(72,2,1:6)=(/ 0.765139,0.781069,0.839900,0.950760,0.973035,0.990361 /)
XCGA_LKT(72,2,1:6)=(/ 0.882343,0.873727,0.844610,0.797210,0.764030,0.754010 /)
XEXT_COEFF_550_LKT(72,2)=203.220000 !rg=3.06299 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(72,3,1:6)=(/ 181.370000,184.860000,188.520000,193.660000,206.240000,266.840000 /)
XPIZA_LKT(72,3,1:6)=(/ 0.766103,0.776167,0.831911,0.947800,0.970904,0.988545 /)
XCGA_LKT(72,3,1:6)=(/ 0.883320,0.874413,0.848057,0.800157,0.760100,0.734407 /)
XEXT_COEFF_550_LKT(72,3)=188.500000 !rg=3.06299 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(72,4,1:6)=(/ 163.830000,166.800000,169.960000,174.400000,182.900000,225.480000 /)
XPIZA_LKT(72,4,1:6)=(/ 0.767799,0.775238,0.821146,0.943341,0.968114,0.986420 /)
XCGA_LKT(72,4,1:6)=(/ 0.882810,0.875247,0.852673,0.804527,0.764923,0.725130 /)
XEXT_COEFF_550_LKT(72,4)=170.030000 !rg=3.06299 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(72,5,1:6)=(/ 144.960000,147.270000,149.890000,153.790000,163.100000,192.040000 /)
XPIZA_LKT(72,5,1:6)=(/ 0.767170,0.774217,0.810549,0.936790,0.962952,0.983741 /)
XCGA_LKT(72,5,1:6)=(/ 0.882623,0.876497,0.857437,0.810507,0.768917,0.725163 /)
XEXT_COEFF_550_LKT(72,5)=150.030000 !rg=3.06299 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(72,6,1:6)=(/ 126.210000,128.100000,130.470000,134.670000,139.380000,161.170000 /)
XPIZA_LKT(72,6,1:6)=(/ 0.764291,0.773509,0.801567,0.929207,0.959561,0.981239 /)
XCGA_LKT(72,6,1:6)=(/ 0.882747,0.877620,0.862223,0.812123,0.782543,0.731697 /)
XEXT_COEFF_550_LKT(72,6)=130.370000 !rg=3.06299 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(72,7,1:6)=(/ 108.570000,110.320000,112.020000,115.050000,119.420000,134.880000 /)
XPIZA_LKT(72,7,1:6)=(/ 0.760311,0.773171,0.794181,0.922504,0.955531,0.977934 /)
XCGA_LKT(72,7,1:6)=(/ 0.884077,0.878060,0.865850,0.819357,0.789093,0.739577 /)
XEXT_COEFF_550_LKT(72,7)=111.930000 !rg=3.06299 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(72,8,1:6)=(/ 92.882000,93.727000,95.250000,97.590000,101.030000,112.290000 /)
XPIZA_LKT(72,8,1:6)=(/ 0.754131,0.770142,0.788896,0.916065,0.949283,0.974395 /)
XCGA_LKT(72,8,1:6)=(/ 0.885363,0.879443,0.869163,0.824353,0.800213,0.750337 /)
XEXT_COEFF_550_LKT(72,8)=95.330000 !rg=3.06299 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(72,9,1:6)=(/ 78.437000,79.478000,80.521000,82.294000,85.776000,94.122000 /)
XPIZA_LKT(72,9,1:6)=(/ 0.745748,0.765682,0.783410,0.908992,0.941467,0.969000 /)
XCGA_LKT(72,9,1:6)=(/ 0.886543,0.880933,0.871827,0.830670,0.804930,0.757093 /)
XEXT_COEFF_550_LKT(72,9)=80.347000 !rg=3.06299 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(72,10,1:6)=(/ 66.189000,66.858000,67.760000,69.182000,71.749000,77.595000 /)
XPIZA_LKT(72,10,1:6)=(/ 0.736256,0.760230,0.778016,0.903817,0.933386,0.963922 /)
XCGA_LKT(72,10,1:6)=(/ 0.888543,0.882227,0.874273,0.835347,0.811917,0.765683 /)
XEXT_COEFF_550_LKT(72,10)=67.900000 !rg=3.06299 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(72,11,1:6)=(/ 55.565000,56.052000,56.871000,57.801000,59.625000,63.815000 /)
XPIZA_LKT(72,11,1:6)=(/ 0.725392,0.752496,0.772280,0.899157,0.926875,0.958998 /)
XCGA_LKT(72,11,1:6)=(/ 0.890963,0.884387,0.877150,0.839870,0.822850,0.780050 /)
XEXT_COEFF_550_LKT(72,11)=56.479000 !rg=3.06299 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(72,12,1:6)=(/ 46.354000,46.874000,47.222000,48.125000,49.494000,52.579000 /)
XPIZA_LKT(72,12,1:6)=(/ 0.712670,0.744343,0.765422,0.894560,0.919395,0.951309 /)
XCGA_LKT(72,12,1:6)=(/ 0.893380,0.886393,0.879253,0.843650,0.828407,0.787353 /)
XEXT_COEFF_550_LKT(72,12)=47.358000 !rg=3.06299 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(72,13,1:6)=(/ 38.675000,38.907000,39.362000,39.877000,40.922000,43.273000 /)
XPIZA_LKT(72,13,1:6)=(/ 0.699945,0.733065,0.757697,0.890317,0.913616,0.943468 /)
XCGA_LKT(72,13,1:6)=(/ 0.896210,0.888430,0.881887,0.847007,0.836310,0.799303 /)
XEXT_COEFF_550_LKT(72,13)=39.388000 !rg=3.06299 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(72,14,1:6)=(/ 32.201000,32.432000,32.698000,33.112000,33.894000,35.538000 /)
XPIZA_LKT(72,14,1:6)=(/ 0.686285,0.721802,0.748528,0.885671,0.907880,0.936470 /)
XCGA_LKT(72,14,1:6)=(/ 0.899220,0.891103,0.884027,0.850470,0.841217,0.808250 /)
XEXT_COEFF_550_LKT(72,14)=32.741000 !rg=3.06299 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(72,15,1:6)=(/ 26.786000,26.984000,27.171000,27.524000,28.101000,29.213000 /)
XPIZA_LKT(72,15,1:6)=(/ 0.672498,0.709327,0.738275,0.880206,0.902249,0.929236 /)
XCGA_LKT(72,15,1:6)=(/ 0.902453,0.893910,0.886540,0.853587,0.846707,0.818970 /)
XEXT_COEFF_550_LKT(72,15)=27.133000 !rg=3.06299 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(72,16,1:6)=(/ 22.267000,22.372000,22.583000,22.757000,23.281000,24.140000 /)
XPIZA_LKT(72,16,1:6)=(/ 0.659326,0.695530,0.726867,0.873487,0.896965,0.921837 /)
XCGA_LKT(72,16,1:6)=(/ 0.906040,0.896843,0.889563,0.856243,0.850847,0.825440 /)
XEXT_COEFF_550_LKT(72,16)=22.478000 !rg=3.06299 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(72,17,1:6)=(/ 18.501000,18.594000,18.714000,18.888000,19.215000,19.879000 /)
XPIZA_LKT(72,17,1:6)=(/ 0.646124,0.682179,0.714395,0.865914,0.891457,0.915873 /)
XCGA_LKT(72,17,1:6)=(/ 0.909623,0.899927,0.892157,0.859157,0.854497,0.831533 /)
XEXT_COEFF_550_LKT(72,17)=18.743000 !rg=3.06299 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(72,18,1:6)=(/ 15.380000,15.463000,15.536000,15.680000,15.940000,16.392000 /)
XPIZA_LKT(72,18,1:6)=(/ 0.633902,0.669039,0.701388,0.856942,0.884944,0.909886 /)
XCGA_LKT(72,18,1:6)=(/ 0.913207,0.903360,0.895090,0.861543,0.857787,0.838960 /)
XEXT_COEFF_550_LKT(72,18)=15.523000 !rg=3.06299 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(72,19,1:6)=(/ 12.779000,12.826000,12.907000,12.986000,13.201000,13.557000 /)
XPIZA_LKT(72,19,1:6)=(/ 0.622608,0.655418,0.688074,0.846429,0.878426,0.904006 /)
XCGA_LKT(72,19,1:6)=(/ 0.917033,0.906903,0.898520,0.863830,0.861303,0.843663 /)
XEXT_COEFF_550_LKT(72,19)=12.877000 !rg=3.06299 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(72,20,1:6)=(/ 10.612000,10.649000,10.702000,10.776000,10.905000,11.172000 /)
XPIZA_LKT(72,20,1:6)=(/ 0.611861,0.642685,0.674506,0.834662,0.870565,0.898457 /)
XCGA_LKT(72,20,1:6)=(/ 0.920783,0.910527,0.901740,0.866567,0.864207,0.848543 /)
XEXT_COEFF_550_LKT(72,20)=10.711000 !rg=3.06299 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(73,1,1:6)=(/ 188.340000,192.570000,192.720000,187.800000,255.510000,295.270000 /)
XPIZA_LKT(73,1,1:6)=(/ 0.773258,0.775407,0.832234,0.947109,0.973850,0.989641 /)
XCGA_LKT(73,1,1:6)=(/ 0.881267,0.875727,0.843343,0.782217,0.789337,0.742883 /)
XEXT_COEFF_550_LKT(73,1)=192.610000 !rg=3.31831 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(73,2,1:6)=(/ 180.080000,183.160000,187.260000,193.930000,214.060000,262.970000 /)
XPIZA_LKT(73,2,1:6)=(/ 0.769939,0.767528,0.831177,0.943764,0.972139,0.988340 /)
XCGA_LKT(73,2,1:6)=(/ 0.881997,0.878270,0.850237,0.796290,0.786103,0.725813 /)
XEXT_COEFF_550_LKT(73,2)=186.910000 !rg=3.31831 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(73,3,1:6)=(/ 166.760000,169.430000,173.330000,179.400000,189.150000,226.900000 /)
XPIZA_LKT(73,3,1:6)=(/ 0.769210,0.767176,0.822623,0.942272,0.969274,0.986478 /)
XCGA_LKT(73,3,1:6)=(/ 0.882553,0.878817,0.852930,0.804703,0.778307,0.714327 /)
XEXT_COEFF_550_LKT(73,3)=172.490000 !rg=3.31831 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(73,4,1:6)=(/ 150.580000,152.960000,156.360000,161.400000,167.760000,196.940000 /)
XPIZA_LKT(73,4,1:6)=(/ 0.768436,0.768803,0.811907,0.937587,0.966304,0.984467 /)
XCGA_LKT(73,4,1:6)=(/ 0.882537,0.878667,0.857130,0.808870,0.778787,0.714877 /)
XEXT_COEFF_550_LKT(73,4)=155.520000 !rg=3.31831 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(73,5,1:6)=(/ 133.580000,135.480000,138.020000,141.420000,147.550000,169.620000 /)
XPIZA_LKT(73,5,1:6)=(/ 0.766523,0.772094,0.801693,0.933101,0.962430,0.982075 /)
XCGA_LKT(73,5,1:6)=(/ 0.883263,0.877713,0.861687,0.813347,0.780780,0.721563 /)
XEXT_COEFF_550_LKT(73,5)=137.930000 !rg=3.31831 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(73,6,1:6)=(/ 116.290000,118.060000,120.110000,123.180000,127.730000,144.850000 /)
XPIZA_LKT(73,6,1:6)=(/ 0.762695,0.771784,0.794078,0.926014,0.958084,0.979410 /)
XCGA_LKT(73,6,1:6)=(/ 0.883400,0.878930,0.865597,0.819067,0.787527,0.731263 /)
XEXT_COEFF_550_LKT(73,6)=120.000000 !rg=3.31831 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(73,7,1:6)=(/ 99.906000,101.470000,103.240000,105.700000,109.080000,121.880000 /)
XPIZA_LKT(73,7,1:6)=(/ 0.757167,0.770398,0.789363,0.917700,0.952999,0.975961 /)
XCGA_LKT(73,7,1:6)=(/ 0.884500,0.879387,0.868400,0.823563,0.797130,0.742743 /)
XEXT_COEFF_550_LKT(73,7)=102.740000 !rg=3.31831 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(73,8,1:6)=(/ 85.422000,86.484000,87.621000,89.912000,92.918000,102.170000 /)
XPIZA_LKT(73,8,1:6)=(/ 0.750114,0.767945,0.784231,0.911097,0.945587,0.971992 /)
XCGA_LKT(73,8,1:6)=(/ 0.885963,0.880520,0.871607,0.828443,0.802033,0.750897 /)
XEXT_COEFF_550_LKT(73,8)=87.800000 !rg=3.31831 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(73,9,1:6)=(/ 72.470000,73.060000,74.207000,75.781000,78.232000,84.803000 /)
XPIZA_LKT(73,9,1:6)=(/ 0.741618,0.763023,0.779005,0.905875,0.938380,0.967630 /)
XCGA_LKT(73,9,1:6)=(/ 0.887990,0.881290,0.874153,0.833547,0.811210,0.763407 /)
XEXT_COEFF_550_LKT(73,9)=74.211000 !rg=3.31831 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(73,10,1:6)=(/ 61.014000,61.635000,62.363000,63.471000,65.775000,70.990000 /)
XPIZA_LKT(73,10,1:6)=(/ 0.730913,0.756500,0.774693,0.900933,0.930143,0.961206 /)
XCGA_LKT(73,10,1:6)=(/ 0.889757,0.883527,0.876157,0.837587,0.818217,0.772580 /)
XEXT_COEFF_550_LKT(73,10)=62.419000 !rg=3.31831 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(73,11,1:6)=(/ 51.123000,51.710000,52.239000,53.160000,54.979000,58.585000 /)
XPIZA_LKT(73,11,1:6)=(/ 0.719358,0.749124,0.768445,0.896191,0.922836,0.955292 /)
XCGA_LKT(73,11,1:6)=(/ 0.891897,0.885463,0.878827,0.841993,0.824017,0.780527 /)
XEXT_COEFF_550_LKT(73,11)=52.175000 !rg=3.31831 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(73,12,1:6)=(/ 42.723000,43.110000,43.631000,44.127000,45.397000,48.210000 /)
XPIZA_LKT(73,12,1:6)=(/ 0.706486,0.739262,0.762373,0.892370,0.916506,0.947590 /)
XCGA_LKT(73,12,1:6)=(/ 0.894813,0.887243,0.880890,0.845530,0.832817,0.793683 /)
XEXT_COEFF_550_LKT(73,12)=43.609000 !rg=3.31831 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(73,13,1:6)=(/ 35.638000,35.912000,36.191000,36.632000,37.642000,39.555000 /)
XPIZA_LKT(73,13,1:6)=(/ 0.693496,0.728131,0.753456,0.887936,0.910634,0.940583 /)
XCGA_LKT(73,13,1:6)=(/ 0.897733,0.889880,0.882940,0.849417,0.838687,0.803297 /)
XEXT_COEFF_550_LKT(73,13)=36.269000 !rg=3.31831 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(73,14,1:6)=(/ 29.689000,29.918000,30.132000,30.470000,31.203000,32.613000 /)
XPIZA_LKT(73,14,1:6)=(/ 0.679915,0.716339,0.744014,0.882855,0.904689,0.932279 /)
XCGA_LKT(73,14,1:6)=(/ 0.900673,0.892600,0.885660,0.852210,0.844183,0.813913 /)
XEXT_COEFF_550_LKT(73,14)=30.145000 !rg=3.31831 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(73,15,1:6)=(/ 24.724000,24.860000,25.052000,25.339000,25.834000,26.821000 /)
XPIZA_LKT(73,15,1:6)=(/ 0.666456,0.703187,0.733405,0.877273,0.899568,0.925179 /)
XCGA_LKT(73,15,1:6)=(/ 0.904167,0.895007,0.888133,0.854360,0.849507,0.822730 /)
XEXT_COEFF_550_LKT(73,15)=25.070000 !rg=3.31831 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(73,16,1:6)=(/ 20.516000,20.643000,20.772000,21.032000,21.382000,22.196000 /)
XPIZA_LKT(73,16,1:6)=(/ 0.652913,0.689832,0.721038,0.870253,0.894157,0.918407 /)
XCGA_LKT(73,16,1:6)=(/ 0.907600,0.898273,0.890827,0.857520,0.852033,0.826777 /)
XEXT_COEFF_550_LKT(73,16)=20.769000 !rg=3.31831 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(73,17,1:6)=(/ 17.058000,17.155000,17.244000,17.396000,17.700000,18.274000 /)
XPIZA_LKT(73,17,1:6)=(/ 0.640197,0.676222,0.708332,0.861964,0.888576,0.912809 /)
XCGA_LKT(73,17,1:6)=(/ 0.911223,0.901560,0.893570,0.860297,0.856590,0.835383 /)
XEXT_COEFF_550_LKT(73,17)=17.252000 !rg=3.31831 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(73,18,1:6)=(/ 14.192000,14.253000,14.339000,14.439000,14.668000,15.055000 /)
XPIZA_LKT(73,18,1:6)=(/ 0.628497,0.662726,0.695700,0.852468,0.882364,0.906946 /)
XCGA_LKT(73,18,1:6)=(/ 0.915017,0.904880,0.896690,0.862227,0.860493,0.841943 /)
XEXT_COEFF_550_LKT(73,18)=14.342000 !rg=3.31831 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(73,19,1:6)=(/ 11.783000,11.835000,11.889000,11.996000,12.137000,12.463000 /)
XPIZA_LKT(73,19,1:6)=(/ 0.617398,0.649876,0.681690,0.841573,0.874618,0.901079 /)
XCGA_LKT(73,19,1:6)=(/ 0.918720,0.908523,0.899963,0.865197,0.862083,0.844753 /)
XEXT_COEFF_550_LKT(73,19)=11.885000 !rg=3.31831 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(73,20,1:6)=(/ 9.790400,9.830300,9.865800,9.929900,10.049000,10.270000 /)
XPIZA_LKT(73,20,1:6)=(/ 0.607205,0.637415,0.668272,0.829048,0.866644,0.895450 /)
XCGA_LKT(73,20,1:6)=(/ 0.922480,0.912303,0.903327,0.867683,0.865567,0.850787 /)
XEXT_COEFF_550_LKT(73,20)=9.868300 !rg=3.31831 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(74,1,1:6)=(/ 174.250000,176.140000,181.950000,189.270000,218.110000,235.440000 /)
XPIZA_LKT(74,1,1:6)=(/ 0.777475,0.758477,0.827891,0.947301,0.971233,0.987017 /)
XCGA_LKT(74,1,1:6)=(/ 0.880113,0.882613,0.854793,0.807443,0.780947,0.704023 /)
XEXT_COEFF_550_LKT(74,1)=180.820000 !rg=3.59491 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(74,2,1:6)=(/ 166.530000,168.980000,172.860000,176.590000,190.510000,214.230000 /)
XPIZA_LKT(74,2,1:6)=(/ 0.774337,0.755587,0.821097,0.943806,0.969361,0.985681 /)
XCGA_LKT(74,2,1:6)=(/ 0.880817,0.883133,0.855177,0.807070,0.788377,0.693683 /)
XEXT_COEFF_550_LKT(74,2)=173.020000 !rg=3.59491 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(74,3,1:6)=(/ 153.890000,156.480000,158.810000,163.530000,172.120000,192.960000 /)
XPIZA_LKT(74,3,1:6)=(/ 0.771387,0.761340,0.811882,0.940649,0.966467,0.983954 /)
XCGA_LKT(74,3,1:6)=(/ 0.882253,0.881583,0.858577,0.810517,0.788453,0.692223 /)
XEXT_COEFF_550_LKT(74,3)=158.220000 !rg=3.59491 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(74,4,1:6)=(/ 138.840000,141.210000,143.330000,147.100000,153.890000,173.080000 /)
XPIZA_LKT(74,4,1:6)=(/ 0.768761,0.766169,0.801445,0.935445,0.963299,0.982153 /)
XCGA_LKT(74,4,1:6)=(/ 0.882657,0.880540,0.862430,0.813840,0.787877,0.705273 /)
XEXT_COEFF_550_LKT(74,4)=142.650000 !rg=3.59491 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(74,5,1:6)=(/ 122.660000,124.410000,126.940000,130.370000,134.910000,151.430000 /)
XPIZA_LKT(74,5,1:6)=(/ 0.765451,0.768596,0.794590,0.927420,0.960319,0.980140 /)
XCGA_LKT(74,5,1:6)=(/ 0.883203,0.880020,0.865370,0.817580,0.791547,0.721873 /)
XEXT_COEFF_550_LKT(74,5)=126.340000 !rg=3.59491 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(74,6,1:6)=(/ 107.020000,108.420000,110.510000,112.660000,117.940000,130.940000 /)
XPIZA_LKT(74,6,1:6)=(/ 0.760617,0.770517,0.788845,0.921266,0.954321,0.976817 /)
XCGA_LKT(74,6,1:6)=(/ 0.883867,0.879890,0.868837,0.821943,0.795007,0.733153 /)
XEXT_COEFF_550_LKT(74,6)=110.580000 !rg=3.59491 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(74,7,1:6)=(/ 92.251000,93.471000,94.815000,96.909000,100.490000,110.610000 /)
XPIZA_LKT(74,7,1:6)=(/ 0.754471,0.768316,0.783809,0.913943,0.948888,0.973456 /)
XCGA_LKT(74,7,1:6)=(/ 0.885347,0.880923,0.871350,0.827580,0.802803,0.745953 /)
XEXT_COEFF_550_LKT(74,7)=94.217000 !rg=3.59491 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(74,8,1:6)=(/ 78.830000,79.502000,80.769000,82.285000,85.646000,93.385000 /)
XPIZA_LKT(74,8,1:6)=(/ 0.746571,0.765164,0.780533,0.907602,0.941450,0.968906 /)
XCGA_LKT(74,8,1:6)=(/ 0.887073,0.881300,0.873650,0.831273,0.808373,0.756370 /)
XEXT_COEFF_550_LKT(74,8)=80.993000 !rg=3.59491 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(74,9,1:6)=(/ 66.574000,67.286000,68.214000,69.509000,71.545000,77.096000 /)
XPIZA_LKT(74,9,1:6)=(/ 0.736593,0.759979,0.776369,0.902897,0.935147,0.965483 /)
XCGA_LKT(74,9,1:6)=(/ 0.888613,0.882860,0.875277,0.836917,0.818283,0.770773 /)
XEXT_COEFF_550_LKT(74,9)=68.057000 !rg=3.59491 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(74,10,1:6)=(/ 56.271000,56.720000,57.452000,58.397000,60.242000,64.331000 /)
XPIZA_LKT(74,10,1:6)=(/ 0.726143,0.753121,0.770852,0.898704,0.926868,0.959508 /)
XCGA_LKT(74,10,1:6)=(/ 0.890733,0.884583,0.877923,0.840577,0.822980,0.778437 /)
XEXT_COEFF_550_LKT(74,10)=57.506000 !rg=3.59491 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(74,11,1:6)=(/ 47.168000,47.621000,48.165000,48.972000,50.413000,53.647000 /)
XPIZA_LKT(74,11,1:6)=(/ 0.713888,0.744728,0.765687,0.894154,0.919444,0.951419 /)
XCGA_LKT(74,11,1:6)=(/ 0.893283,0.886563,0.880203,0.844467,0.828883,0.787247 /)
XEXT_COEFF_550_LKT(74,11)=47.987000 !rg=3.59491 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(74,12,1:6)=(/ 39.412000,39.691000,40.115000,40.648000,41.698000,43.943000 /)
XPIZA_LKT(74,12,1:6)=(/ 0.700824,0.734177,0.758133,0.890236,0.913621,0.944834 /)
XCGA_LKT(74,12,1:6)=(/ 0.895913,0.888647,0.881877,0.848003,0.836260,0.798727 /)
XEXT_COEFF_550_LKT(74,12)=40.178000 !rg=3.59491 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(74,13,1:6)=(/ 32.850000,33.149000,33.396000,33.729000,34.631000,36.290000 /)
XPIZA_LKT(74,13,1:6)=(/ 0.687169,0.723030,0.749284,0.885240,0.907632,0.936470 /)
XCGA_LKT(74,13,1:6)=(/ 0.898930,0.891347,0.884637,0.851037,0.842197,0.809563 /)
XEXT_COEFF_550_LKT(74,13)=33.396000 !rg=3.59491 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(74,14,1:6)=(/ 27.399000,27.545000,27.800000,28.069000,28.631000,29.794000 /)
XPIZA_LKT(74,14,1:6)=(/ 0.673756,0.710063,0.739547,0.880351,0.902027,0.928930 /)
XCGA_LKT(74,14,1:6)=(/ 0.902337,0.893680,0.887100,0.853583,0.847490,0.818867 /)
XEXT_COEFF_550_LKT(74,14)=27.768000 !rg=3.59491 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(74,15,1:6)=(/ 22.788000,22.920000,23.085000,23.351000,23.817000,24.775000 /)
XPIZA_LKT(74,15,1:6)=(/ 0.660114,0.697115,0.728055,0.874242,0.896298,0.921061 /)
XCGA_LKT(74,15,1:6)=(/ 0.905637,0.896553,0.889093,0.856607,0.850520,0.823093 /)
XEXT_COEFF_550_LKT(74,15)=23.110000 !rg=3.59491 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(74,16,1:6)=(/ 18.919000,19.029000,19.165000,19.369000,19.698000,20.408000 /)
XPIZA_LKT(74,16,1:6)=(/ 0.646924,0.683563,0.715920,0.866731,0.891477,0.915192 /)
XCGA_LKT(74,16,1:6)=(/ 0.909257,0.899823,0.892173,0.858597,0.854300,0.830697 /)
XEXT_COEFF_550_LKT(74,16)=19.119000 !rg=3.59491 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(74,17,1:6)=(/ 15.741000,15.809000,15.914000,16.044000,16.302000,16.741000 /)
XPIZA_LKT(74,17,1:6)=(/ 0.634748,0.669854,0.702826,0.858001,0.885489,0.909371 /)
XCGA_LKT(74,17,1:6)=(/ 0.912970,0.903070,0.895097,0.861313,0.858453,0.838953 /)
XEXT_COEFF_550_LKT(74,17)=15.921000 !rg=3.59491 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(74,18,1:6)=(/ 13.089000,13.147000,13.212000,13.332000,13.507000,13.887000 /)
XPIZA_LKT(74,18,1:6)=(/ 0.623112,0.656757,0.689340,0.847968,0.878763,0.903499 /)
XCGA_LKT(74,18,1:6)=(/ 0.916703,0.906613,0.898013,0.864003,0.861300,0.842687 /)
XEXT_COEFF_550_LKT(74,18)=13.221000 !rg=3.59491 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(74,19,1:6)=(/ 10.868000,10.914000,10.972000,11.045000,11.176000,11.474000 /)
XPIZA_LKT(74,19,1:6)=(/ 0.612435,0.643920,0.676122,0.836309,0.871138,0.898199 /)
XCGA_LKT(74,19,1:6)=(/ 0.920440,0.910313,0.901510,0.866307,0.863773,0.847487 /)
XEXT_COEFF_550_LKT(74,19)=10.956000 !rg=3.59491 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(74,20,1:6)=(/ 9.033300,9.062200,9.106900,9.159400,9.271500,9.454200 /)
XPIZA_LKT(74,20,1:6)=(/ 0.602874,0.631847,0.662603,0.823124,0.862650,0.892242 /)
XCGA_LKT(74,20,1:6)=(/ 0.924163,0.913933,0.905060,0.868833,0.867103,0.853070 /)
XEXT_COEFF_550_LKT(74,20)=9.106000 !rg=3.59491 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(75,1,1:6)=(/ 159.730000,162.660000,163.900000,175.100000,179.350000,184.800000 /)
XPIZA_LKT(75,1,1:6)=(/ 0.779284,0.737236,0.814689,0.944541,0.957084,0.983139 /)
XCGA_LKT(75,1,1:6)=(/ 0.879927,0.889963,0.857783,0.820543,0.756777,0.650320 /)
XEXT_COEFF_550_LKT(75,1)=164.720000 !rg=3.89457 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(75,2,1:6)=(/ 152.600000,155.210000,158.610000,162.120000,168.580000,175.190000 /)
XPIZA_LKT(75,2,1:6)=(/ 0.776462,0.751138,0.810868,0.938961,0.963953,0.982254 /)
XCGA_LKT(75,2,1:6)=(/ 0.880117,0.886053,0.858977,0.809980,0.792470,0.655363 /)
XEXT_COEFF_550_LKT(75,2)=158.100000 !rg=3.89457 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(75,3,1:6)=(/ 141.450000,143.600000,146.230000,150.120000,158.710000,165.750000 /)
XPIZA_LKT(75,3,1:6)=(/ 0.772190,0.761177,0.800344,0.934915,0.961798,0.981405 /)
XCGA_LKT(75,3,1:6)=(/ 0.881460,0.882543,0.864093,0.811253,0.780913,0.680293 /)
XEXT_COEFF_550_LKT(75,3)=146.230000 !rg=3.89457 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(75,4,1:6)=(/ 127.700000,129.810000,132.010000,135.270000,142.530000,153.220000 /)
XPIZA_LKT(75,4,1:6)=(/ 0.768593,0.765671,0.791816,0.929538,0.958559,0.979927 /)
XCGA_LKT(75,4,1:6)=(/ 0.882317,0.881703,0.867400,0.815777,0.783713,0.703163 /)
XEXT_COEFF_550_LKT(75,4)=131.730000 !rg=3.89457 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(75,5,1:6)=(/ 113.120000,114.870000,116.430000,119.350000,123.780000,136.090000 /)
XPIZA_LKT(75,5,1:6)=(/ 0.763958,0.767944,0.786358,0.923593,0.956892,0.977753 /)
XCGA_LKT(75,5,1:6)=(/ 0.883577,0.880850,0.869810,0.822067,0.798523,0.722580 /)
XEXT_COEFF_550_LKT(75,5)=115.940000 !rg=3.89457 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(75,6,1:6)=(/ 98.807000,99.918000,101.310000,103.030000,108.670000,118.240000 /)
XPIZA_LKT(75,6,1:6)=(/ 0.758010,0.767807,0.783801,0.916978,0.951044,0.974557 /)
XCGA_LKT(75,6,1:6)=(/ 0.884823,0.881273,0.871673,0.827503,0.795730,0.734133 /)
XEXT_COEFF_550_LKT(75,6)=101.320000 !rg=3.89457 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(75,7,1:6)=(/ 84.841000,86.091000,87.297000,89.107000,92.956000,100.600000 /)
XPIZA_LKT(75,7,1:6)=(/ 0.750418,0.767019,0.779934,0.909139,0.944203,0.970888 /)
XCGA_LKT(75,7,1:6)=(/ 0.885873,0.881740,0.874193,0.830657,0.802470,0.746503 /)
XEXT_COEFF_550_LKT(75,7)=86.950000 !rg=3.89457 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(75,8,1:6)=(/ 72.542000,73.346000,74.186000,75.440000,78.377000,84.551000 /)
XPIZA_LKT(75,8,1:6)=(/ 0.742039,0.763090,0.776910,0.904307,0.938717,0.967435 /)
XCGA_LKT(75,8,1:6)=(/ 0.887853,0.882697,0.875307,0.835807,0.813087,0.760820 /)
XEXT_COEFF_550_LKT(75,8)=74.395000 !rg=3.89457 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(75,9,1:6)=(/ 61.459000,62.054000,62.738000,64.063000,65.821000,70.392000 /)
XPIZA_LKT(75,9,1:6)=(/ 0.731988,0.756570,0.772944,0.899633,0.930163,0.962206 /)
XCGA_LKT(75,9,1:6)=(/ 0.889767,0.883713,0.877510,0.839010,0.822527,0.775787 /)
XEXT_COEFF_550_LKT(75,9)=62.634000 !rg=3.89457 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(75,10,1:6)=(/ 51.828000,52.355000,52.832000,53.689000,55.402000,58.691000 /)
XPIZA_LKT(75,10,1:6)=(/ 0.720186,0.749764,0.768074,0.895683,0.922698,0.956361 /)
XCGA_LKT(75,10,1:6)=(/ 0.891957,0.885680,0.879437,0.843770,0.827997,0.786517 /)
XEXT_COEFF_550_LKT(75,10)=52.815000 !rg=3.89457 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(75,11,1:6)=(/ 43.525000,43.880000,44.339000,45.109000,46.116000,48.639000 /)
XPIZA_LKT(75,11,1:6)=(/ 0.707974,0.739975,0.761970,0.891537,0.916501,0.948672 /)
XCGA_LKT(75,11,1:6)=(/ 0.894827,0.887377,0.881520,0.846293,0.834117,0.794720 /)
XEXT_COEFF_550_LKT(75,11)=44.241000 !rg=3.89457 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(75,12,1:6)=(/ 36.325000,36.627000,36.909000,37.429000,38.391000,40.272000 /)
XPIZA_LKT(75,12,1:6)=(/ 0.694214,0.729485,0.753808,0.887574,0.909814,0.940522 /)
XCGA_LKT(75,12,1:6)=(/ 0.897430,0.889930,0.883550,0.849913,0.839653,0.805510 /)
XEXT_COEFF_550_LKT(75,12)=36.896000 !rg=3.89457 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(75,13,1:6)=(/ 30.339000,30.499000,30.808000,31.181000,31.752000,33.099000 /)
XPIZA_LKT(75,13,1:6)=(/ 0.681261,0.716976,0.744900,0.883112,0.904374,0.932734 /)
XCGA_LKT(75,13,1:6)=(/ 0.900590,0.892307,0.886127,0.852637,0.845603,0.814973 /)
XEXT_COEFF_550_LKT(75,13)=30.774000 !rg=3.89457 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(75,14,1:6)=(/ 25.246000,25.413000,25.579000,25.926000,26.400000,27.432000 /)
XPIZA_LKT(75,14,1:6)=(/ 0.667090,0.704448,0.733997,0.877635,0.899151,0.925192 /)
XCGA_LKT(75,14,1:6)=(/ 0.903843,0.895080,0.888050,0.855197,0.848763,0.820633 /)
XEXT_COEFF_550_LKT(75,14)=25.602000 !rg=3.89457 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(75,15,1:6)=(/ 21.012000,21.143000,21.277000,21.427000,21.856000,22.672000 /)
XPIZA_LKT(75,15,1:6)=(/ 0.653797,0.691007,0.722617,0.870887,0.893900,0.918448 /)
XCGA_LKT(75,15,1:6)=(/ 0.907297,0.898037,0.890597,0.857397,0.852800,0.827727 /)
XEXT_COEFF_550_LKT(75,15)=21.263000 !rg=3.89457 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(75,16,1:6)=(/ 17.458000,17.546000,17.659000,17.822000,18.117000,18.656000 /)
XPIZA_LKT(75,16,1:6)=(/ 0.641129,0.677323,0.709962,0.862723,0.888656,0.912768 /)
XCGA_LKT(75,16,1:6)=(/ 0.911060,0.901207,0.893630,0.860007,0.856533,0.835423 /)
XEXT_COEFF_550_LKT(75,16)=17.649000 !rg=3.89457 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(75,17,1:6)=(/ 14.515000,14.589000,14.656000,14.797000,14.960000,15.378000 /)
XPIZA_LKT(75,17,1:6)=(/ 0.628999,0.663927,0.696539,0.853861,0.882339,0.906855 /)
XCGA_LKT(75,17,1:6)=(/ 0.914700,0.904750,0.896437,0.862813,0.859623,0.841110 /)
XEXT_COEFF_550_LKT(75,17)=14.649000 !rg=3.89457 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(75,18,1:6)=(/ 12.074000,12.124000,12.184000,12.245000,12.405000,12.738000 /)
XPIZA_LKT(75,18,1:6)=(/ 0.617939,0.650765,0.683234,0.842869,0.875608,0.901056 /)
XCGA_LKT(75,18,1:6)=(/ 0.918457,0.908250,0.899680,0.864877,0.862923,0.845533 /)
XEXT_COEFF_550_LKT(75,18)=12.177000 !rg=3.89457 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(75,19,1:6)=(/ 10.029000,10.063000,10.113000,10.173000,10.309000,10.528000 /)
XPIZA_LKT(75,19,1:6)=(/ 0.607856,0.638184,0.669756,0.830477,0.867596,0.895754 /)
XCGA_LKT(75,19,1:6)=(/ 0.922193,0.911883,0.903180,0.867453,0.865217,0.850770 /)
XEXT_COEFF_550_LKT(75,19)=10.115000 !rg=3.89457 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(75,20,1:6)=(/ 8.333000,8.362700,8.392200,8.451600,8.514900,8.686900 /)
XPIZA_LKT(75,20,1:6)=(/ 0.598506,0.626662,0.656572,0.817272,0.858068,0.889278 /)
XCGA_LKT(75,20,1:6)=(/ 0.925863,0.915697,0.906557,0.870257,0.867817,0.854240 /)
XEXT_COEFF_550_LKT(75,20)=8.389100 !rg=3.89457 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(76,1,1:6)=(/ 147.140000,150.170000,154.370000,150.070000,139.460000,142.520000 /)
XPIZA_LKT(76,1,1:6)=(/ 0.776457,0.736825,0.806877,0.935689,0.959522,0.978438 /)
XCGA_LKT(76,1,1:6)=(/ 0.881400,0.890853,0.864457,0.808993,0.767063,0.595993 /)
XEXT_COEFF_550_LKT(76,1)=154.720000 !rg=4.21921 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(76,2,1:6)=(/ 140.960000,143.110000,145.010000,147.800000,154.520000,145.840000 /)
XPIZA_LKT(76,2,1:6)=(/ 0.774974,0.749119,0.799447,0.935394,0.959466,0.978861 /)
XCGA_LKT(76,2,1:6)=(/ 0.881660,0.887480,0.865060,0.816403,0.773850,0.633140 /)
XEXT_COEFF_550_LKT(76,2)=144.390000 !rg=4.21921 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(76,3,1:6)=(/ 130.430000,132.360000,134.640000,137.570000,145.550000,147.160000 /)
XPIZA_LKT(76,3,1:6)=(/ 0.771782,0.759047,0.790240,0.930486,0.958173,0.978777 /)
XCGA_LKT(76,3,1:6)=(/ 0.881673,0.883907,0.869180,0.818010,0.792587,0.676857 /)
XEXT_COEFF_550_LKT(76,3)=134.620000 !rg=4.21921 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(76,4,1:6)=(/ 117.960000,119.500000,121.430000,124.240000,130.880000,138.420000 /)
XPIZA_LKT(76,4,1:6)=(/ 0.767143,0.764796,0.784079,0.924641,0.955881,0.977147 /)
XCGA_LKT(76,4,1:6)=(/ 0.882873,0.882420,0.871613,0.821910,0.795223,0.706283 /)
XEXT_COEFF_550_LKT(76,4)=121.110000 !rg=4.21921 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(76,5,1:6)=(/ 104.230000,105.450000,107.360000,109.680000,114.920000,123.150000 /)
XPIZA_LKT(76,5,1:6)=(/ 0.761576,0.768630,0.780473,0.918173,0.951097,0.975250 /)
XCGA_LKT(76,5,1:6)=(/ 0.883853,0.881073,0.873117,0.825527,0.794983,0.724917 /)
XEXT_COEFF_550_LKT(76,5)=107.150000 !rg=4.21921 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(76,6,1:6)=(/ 90.780000,92.095000,93.288000,95.472000,99.710000,107.480000 /)
XPIZA_LKT(76,6,1:6)=(/ 0.754330,0.766965,0.778125,0.911629,0.947148,0.971908 /)
XCGA_LKT(76,6,1:6)=(/ 0.884770,0.881770,0.874863,0.831323,0.803760,0.743497 /)
XEXT_COEFF_550_LKT(76,6)=93.484000 !rg=4.21921 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(76,7,1:6)=(/ 78.320000,79.230000,80.362000,81.999000,85.137000,92.047000 /)
XPIZA_LKT(76,7,1:6)=(/ 0.746839,0.764757,0.777180,0.905805,0.940818,0.967476 /)
XCGA_LKT(76,7,1:6)=(/ 0.887000,0.882530,0.876120,0.835157,0.810583,0.754157 /)
XEXT_COEFF_550_LKT(76,7)=79.892000 !rg=4.21921 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(76,8,1:6)=(/ 66.838000,67.609000,68.332000,69.543000,72.042000,76.989000 /)
XPIZA_LKT(76,8,1:6)=(/ 0.736911,0.760205,0.773854,0.900678,0.934035,0.964602 /)
XCGA_LKT(76,8,1:6)=(/ 0.888693,0.883823,0.877637,0.839067,0.819020,0.770483 /)
XEXT_COEFF_550_LKT(76,8)=68.452000 !rg=4.21921 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(76,9,1:6)=(/ 56.606000,57.106000,57.868000,58.840000,60.977000,64.809000 /)
XPIZA_LKT(76,9,1:6)=(/ 0.726341,0.753706,0.769779,0.897006,0.924967,0.957533 /)
XCGA_LKT(76,9,1:6)=(/ 0.890747,0.884627,0.878920,0.842467,0.822373,0.775037 /)
XEXT_COEFF_550_LKT(76,9)=57.845000 !rg=4.21921 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(76,10,1:6)=(/ 47.860000,48.152000,48.722000,49.451000,50.706000,53.564000 /)
XPIZA_LKT(76,10,1:6)=(/ 0.714849,0.744833,0.765457,0.893479,0.919305,0.952561 /)
XCGA_LKT(76,10,1:6)=(/ 0.893440,0.886330,0.880793,0.845290,0.832670,0.792143 /)
XEXT_COEFF_550_LKT(76,10)=48.783000 !rg=4.21921 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(76,11,1:6)=(/ 40.067000,40.447000,40.860000,41.420000,42.234000,44.366000 /)
XPIZA_LKT(76,11,1:6)=(/ 0.701361,0.735443,0.758236,0.889630,0.913708,0.945452 /)
XCGA_LKT(76,11,1:6)=(/ 0.895903,0.888653,0.882400,0.848603,0.839100,0.801430 /)
XEXT_COEFF_550_LKT(76,11)=40.745000 !rg=4.21921 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(76,12,1:6)=(/ 33.539000,33.727000,34.050000,34.434000,35.189000,36.770000 /)
XPIZA_LKT(76,12,1:6)=(/ 0.688318,0.723540,0.750187,0.885655,0.907099,0.936751 /)
XCGA_LKT(76,12,1:6)=(/ 0.898970,0.890953,0.884880,0.851257,0.843363,0.810817 /)
XEXT_COEFF_550_LKT(76,12)=34.050000 !rg=4.21921 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(76,13,1:6)=(/ 27.930000,28.156000,28.325000,28.719000,29.224000,30.399000 /)
XPIZA_LKT(76,13,1:6)=(/ 0.674288,0.711809,0.739721,0.880934,0.901746,0.929270 /)
XCGA_LKT(76,13,1:6)=(/ 0.902077,0.893713,0.886957,0.853863,0.847290,0.817170 /)
XEXT_COEFF_550_LKT(76,13)=28.296000 !rg=4.21921 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(76,14,1:6)=(/ 23.297000,23.405000,23.602000,23.836000,24.242000,25.176000 /)
XPIZA_LKT(76,14,1:6)=(/ 0.661251,0.698082,0.729210,0.874537,0.896796,0.922260 /)
XCGA_LKT(76,14,1:6)=(/ 0.905410,0.896337,0.889373,0.856393,0.851943,0.825400 /)
XEXT_COEFF_550_LKT(76,14)=23.593000 !rg=4.21921 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(76,15,1:6)=(/ 19.391000,19.482000,19.628000,19.831000,20.154000,20.803000 /)
XPIZA_LKT(76,15,1:6)=(/ 0.647944,0.684767,0.716964,0.867783,0.891137,0.915493 /)
XCGA_LKT(76,15,1:6)=(/ 0.909007,0.899480,0.892027,0.858800,0.854993,0.831820 /)
XEXT_COEFF_550_LKT(76,15)=19.634000 !rg=4.21921 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(76,16,1:6)=(/ 16.099000,16.195000,16.276000,16.427000,16.624000,17.088000 /)
XPIZA_LKT(76,16,1:6)=(/ 0.635224,0.671390,0.704062,0.859192,0.886025,0.910197 /)
XCGA_LKT(76,16,1:6)=(/ 0.912713,0.902890,0.894767,0.861387,0.859313,0.839753 /)
XEXT_COEFF_550_LKT(76,16)=16.265000 !rg=4.21921 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(76,17,1:6)=(/ 13.393000,13.449000,13.526000,13.604000,13.812000,14.172000 /)
XPIZA_LKT(76,17,1:6)=(/ 0.623718,0.657819,0.690739,0.849056,0.879472,0.903972 /)
XCGA_LKT(76,17,1:6)=(/ 0.916447,0.906257,0.898017,0.863577,0.862257,0.844540 /)
XEXT_COEFF_550_LKT(76,17)=13.512000 !rg=4.21921 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(76,18,1:6)=(/ 11.140000,11.179000,11.239000,11.320000,11.454000,11.734000 /)
XPIZA_LKT(76,18,1:6)=(/ 0.613137,0.644877,0.677180,0.837792,0.871965,0.898274 /)
XCGA_LKT(76,18,1:6)=(/ 0.920130,0.909903,0.901240,0.866060,0.864090,0.848360 /)
XEXT_COEFF_550_LKT(76,18)=11.241000 !rg=4.21921 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(76,19,1:6)=(/ 9.251700,9.289900,9.322300,9.387000,9.466700,9.648700 /)
XPIZA_LKT(76,19,1:6)=(/ 0.603234,0.632918,0.663684,0.824991,0.863539,0.892805 /)
XCGA_LKT(76,19,1:6)=(/ 0.923917,0.913680,0.904583,0.868620,0.867220,0.853490 /)
XEXT_COEFF_550_LKT(76,19)=9.321100 !rg=4.21921 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(76,20,1:6)=(/ 7.688600,7.711900,7.744200,7.780700,7.860400,8.010900 /)
XPIZA_LKT(76,20,1:6)=(/ 0.594537,0.621430,0.650902,0.810690,0.853726,0.886311 /)
XCGA_LKT(76,20,1:6)=(/ 0.927480,0.917367,0.908267,0.871350,0.869793,0.856937 /)
XEXT_COEFF_550_LKT(76,20)=7.738700 !rg=4.21921 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(77,1,1:6)=(/ 135.230000,138.450000,139.750000,147.410000,139.520000,119.500000 /)
XPIZA_LKT(77,1,1:6)=(/ 0.776038,0.747038,0.794279,0.929890,0.956622,0.974345 /)
XCGA_LKT(77,1,1:6)=(/ 0.881143,0.888267,0.865740,0.818647,0.780017,0.574830 /)
XEXT_COEFF_550_LKT(77,1)=139.750000 !rg=4.57091 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(77,2,1:6)=(/ 129.560000,131.830000,133.360000,136.630000,143.020000,130.360000 /)
XPIZA_LKT(77,2,1:6)=(/ 0.773444,0.757283,0.785840,0.931109,0.956493,0.976297 /)
XCGA_LKT(77,2,1:6)=(/ 0.881087,0.885530,0.871330,0.820173,0.793257,0.635933 /)
XEXT_COEFF_550_LKT(77,2)=133.660000 !rg=4.57091 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(77,3,1:6)=(/ 120.490000,121.670000,123.960000,126.690000,131.280000,132.900000 /)
XPIZA_LKT(77,3,1:6)=(/ 0.770434,0.763560,0.778473,0.926425,0.958879,0.976832 /)
XCGA_LKT(77,3,1:6)=(/ 0.881923,0.882630,0.873020,0.821960,0.800497,0.686270 /)
XEXT_COEFF_550_LKT(77,3)=123.960000 !rg=4.57091 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(77,4,1:6)=(/ 108.850000,109.900000,111.810000,114.330000,118.530000,124.580000 /)
XPIZA_LKT(77,4,1:6)=(/ 0.764952,0.767033,0.776425,0.920269,0.955152,0.975706 /)
XCGA_LKT(77,4,1:6)=(/ 0.883777,0.881553,0.874543,0.825290,0.802890,0.716160 /)
XEXT_COEFF_550_LKT(77,4)=111.660000 !rg=4.57091 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(77,5,1:6)=(/ 95.971000,97.307000,98.611000,100.930000,105.370000,112.680000 /)
XPIZA_LKT(77,5,1:6)=(/ 0.758251,0.766665,0.777047,0.912700,0.948936,0.972140 /)
XCGA_LKT(77,5,1:6)=(/ 0.884397,0.882410,0.875417,0.829527,0.805140,0.733360 /)
XEXT_COEFF_550_LKT(77,5)=98.371000 !rg=4.57091 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(77,6,1:6)=(/ 83.747000,84.743000,86.010000,88.001000,90.542000,96.826000 /)
XPIZA_LKT(77,6,1:6)=(/ 0.751004,0.765625,0.775033,0.906575,0.943450,0.969678 /)
XCGA_LKT(77,6,1:6)=(/ 0.886060,0.882667,0.876897,0.833453,0.812183,0.752397 /)
XEXT_COEFF_550_LKT(77,6)=85.975000 !rg=4.57091 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(77,7,1:6)=(/ 72.261000,73.038000,73.945000,75.590000,77.625000,82.848000 /)
XPIZA_LKT(77,7,1:6)=(/ 0.742294,0.762452,0.773934,0.901704,0.937562,0.966233 /)
XCGA_LKT(77,7,1:6)=(/ 0.888273,0.882950,0.877743,0.837320,0.816877,0.763587 /)
XEXT_COEFF_550_LKT(77,7)=73.713000 !rg=4.57091 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(77,8,1:6)=(/ 61.750000,62.140000,62.980000,63.991000,65.804000,69.920000 /)
XPIZA_LKT(77,8,1:6)=(/ 0.732352,0.756646,0.771481,0.898239,0.929686,0.961484 /)
XCGA_LKT(77,8,1:6)=(/ 0.889943,0.884280,0.879240,0.841520,0.824653,0.777337 /)
XEXT_COEFF_550_LKT(77,8)=62.999000 !rg=4.57091 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(77,9,1:6)=(/ 52.175000,52.759000,53.191000,54.077000,55.742000,59.314000 /)
XPIZA_LKT(77,9,1:6)=(/ 0.720542,0.749657,0.767603,0.894008,0.921943,0.954844 /)
XCGA_LKT(77,9,1:6)=(/ 0.891803,0.885867,0.880037,0.844227,0.828280,0.783623 /)
XEXT_COEFF_550_LKT(77,9)=53.131000 !rg=4.57091 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(77,10,1:6)=(/ 44.062000,44.517000,44.819000,45.583000,46.811000,49.338000 /)
XPIZA_LKT(77,10,1:6)=(/ 0.708355,0.741478,0.761280,0.890826,0.915022,0.947915 /)
XCGA_LKT(77,10,1:6)=(/ 0.894523,0.887693,0.881703,0.847680,0.834093,0.792267 /)
XEXT_COEFF_550_LKT(77,10)=44.936000 !rg=4.57091 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(77,11,1:6)=(/ 36.973000,37.243000,37.618000,38.132000,39.048000,40.874000 /)
XPIZA_LKT(77,11,1:6)=(/ 0.695260,0.730047,0.754424,0.887316,0.909255,0.940917 /)
XCGA_LKT(77,11,1:6)=(/ 0.897380,0.889790,0.883767,0.850407,0.840670,0.805810 /)
XEXT_COEFF_550_LKT(77,11)=37.495000 !rg=4.57091 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(77,12,1:6)=(/ 30.894000,31.151000,31.338000,31.823000,32.464000,33.851000 /)
XPIZA_LKT(77,12,1:6)=(/ 0.681565,0.718625,0.745231,0.883059,0.903724,0.932653 /)
XCGA_LKT(77,12,1:6)=(/ 0.900313,0.892307,0.885827,0.853003,0.845023,0.812263 /)
XEXT_COEFF_550_LKT(77,12)=31.385000 !rg=4.57091 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(77,13,1:6)=(/ 25.774000,25.911000,26.150000,26.420000,26.905000,27.969000 /)
XPIZA_LKT(77,13,1:6)=(/ 0.668080,0.705493,0.735480,0.877822,0.899262,0.925851 /)
XCGA_LKT(77,13,1:6)=(/ 0.903663,0.894813,0.888287,0.855197,0.850833,0.822270 /)
XEXT_COEFF_550_LKT(77,13)=26.142000 !rg=4.57091 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(77,14,1:6)=(/ 21.480000,21.594000,21.725000,21.932000,22.350000,23.114000 /)
XPIZA_LKT(77,14,1:6)=(/ 0.654696,0.692066,0.723493,0.871531,0.894118,0.919491 /)
XCGA_LKT(77,14,1:6)=(/ 0.907167,0.897843,0.890493,0.857790,0.853517,0.828733 /)
XEXT_COEFF_550_LKT(77,14)=21.758000 !rg=4.57091 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(77,15,1:6)=(/ 17.879000,17.985000,18.072000,18.258000,18.550000,19.106000 /)
XPIZA_LKT(77,15,1:6)=(/ 0.641828,0.678803,0.711295,0.864017,0.888542,0.912915 /)
XCGA_LKT(77,15,1:6)=(/ 0.910727,0.901103,0.893287,0.860410,0.856590,0.836257 /)
XEXT_COEFF_550_LKT(77,15)=18.055000 !rg=4.57091 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(77,16,1:6)=(/ 14.857000,14.920000,15.023000,15.117000,15.369000,15.799000 /)
XPIZA_LKT(77,16,1:6)=(/ 0.629784,0.664790,0.698339,0.854683,0.882921,0.907097 /)
XCGA_LKT(77,16,1:6)=(/ 0.914487,0.904463,0.896353,0.862800,0.859850,0.841497 /)
XEXT_COEFF_550_LKT(77,16)=14.983000 !rg=4.57091 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(77,17,1:6)=(/ 12.356000,12.401000,12.465000,12.547000,12.729000,13.041000 /)
XPIZA_LKT(77,17,1:6)=(/ 0.618564,0.651651,0.684648,0.844343,0.876472,0.901845 /)
XCGA_LKT(77,17,1:6)=(/ 0.918177,0.907873,0.899340,0.865280,0.863133,0.846153 /)
XEXT_COEFF_550_LKT(77,17)=12.469000 !rg=4.57091 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(77,18,1:6)=(/ 10.277000,10.318000,10.357000,10.437000,10.551000,10.791000 /)
XPIZA_LKT(77,18,1:6)=(/ 0.608270,0.639342,0.671026,0.832358,0.868191,0.895612 /)
XCGA_LKT(77,18,1:6)=(/ 0.921927,0.911660,0.902753,0.867423,0.865057,0.850647 /)
XEXT_COEFF_550_LKT(77,18)=10.349000 !rg=4.57091 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(77,19,1:6)=(/ 8.537900,8.563300,8.606800,8.647500,8.751500,8.942000 /)
XPIZA_LKT(77,19,1:6)=(/ 0.599045,0.627325,0.657996,0.818734,0.859491,0.889970 /)
XCGA_LKT(77,19,1:6)=(/ 0.925607,0.915367,0.906337,0.870100,0.868150,0.854633 /)
XEXT_COEFF_550_LKT(77,19)=8.588000 !rg=4.57091 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(77,20,1:6)=(/ 7.094000,7.113700,7.139100,7.179900,7.246100,7.371300 /)
XPIZA_LKT(77,20,1:6)=(/ 0.590719,0.616428,0.645056,0.804195,0.849167,0.883702 /)
XCGA_LKT(77,20,1:6)=(/ 0.929080,0.919050,0.909860,0.872837,0.870847,0.858463 /)
XEXT_COEFF_550_LKT(77,20)=7.144300 !rg=4.57091 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(78,1,1:6)=(/ 124.720000,126.860000,128.540000,132.510000,145.200000,113.320000 /)
XPIZA_LKT(78,1,1:6)=(/ 0.772805,0.764428,0.781769,0.929533,0.961530,0.972725 /)
XCGA_LKT(78,1,1:6)=(/ 0.880820,0.883137,0.873403,0.825473,0.820463,0.603210 /)
XEXT_COEFF_550_LKT(78,1)=129.030000 !rg=4.95192 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(78,2,1:6)=(/ 119.390000,121.350000,123.440000,126.990000,132.180000,123.820000 /)
XPIZA_LKT(78,2,1:6)=(/ 0.772081,0.764434,0.772708,0.924754,0.954759,0.974977 /)
XCGA_LKT(78,2,1:6)=(/ 0.881860,0.883757,0.876797,0.823090,0.803070,0.662880 /)
XEXT_COEFF_550_LKT(78,2)=123.220000 !rg=4.95192 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(78,3,1:6)=(/ 110.620000,111.770000,114.060000,116.430000,120.090000,123.490000 /)
XPIZA_LKT(78,3,1:6)=(/ 0.767428,0.766516,0.770749,0.921066,0.956391,0.975402 /)
XCGA_LKT(78,3,1:6)=(/ 0.883103,0.882470,0.876907,0.826380,0.806890,0.706450 /)
XEXT_COEFF_550_LKT(78,3)=113.610000 !rg=4.95192 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(78,4,1:6)=(/ 99.850000,101.050000,102.870000,104.870000,108.180000,114.110000 /)
XPIZA_LKT(78,4,1:6)=(/ 0.761688,0.767564,0.771974,0.914509,0.952603,0.973992 /)
XCGA_LKT(78,4,1:6)=(/ 0.884107,0.882190,0.876933,0.829443,0.809740,0.730597 /)
XEXT_COEFF_550_LKT(78,4)=102.430000 !rg=4.95192 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(78,5,1:6)=(/ 88.703000,89.438000,90.891000,92.812000,95.910000,101.590000 /)
XPIZA_LKT(78,5,1:6)=(/ 0.755116,0.767210,0.771554,0.908149,0.946971,0.971318 /)
XCGA_LKT(78,5,1:6)=(/ 0.885507,0.881827,0.877850,0.833063,0.810930,0.744900 /)
XEXT_COEFF_550_LKT(78,5)=90.829000 !rg=4.95192 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(78,6,1:6)=(/ 77.165000,78.158000,79.143000,80.749000,82.832000,88.440000 /)
XPIZA_LKT(78,6,1:6)=(/ 0.746315,0.764446,0.772308,0.903557,0.940787,0.968127 /)
XCGA_LKT(78,6,1:6)=(/ 0.886840,0.883247,0.878547,0.837170,0.816647,0.758153 /)
XEXT_COEFF_550_LKT(78,6)=79.116000 !rg=4.95192 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(78,7,1:6)=(/ 66.420000,67.175000,68.165000,69.261000,70.906000,75.523000 /)
XPIZA_LKT(78,7,1:6)=(/ 0.736637,0.760022,0.771666,0.898348,0.934286,0.964081 /)
XCGA_LKT(78,7,1:6)=(/ 0.888807,0.883697,0.878820,0.840317,0.823647,0.771907 /)
XEXT_COEFF_550_LKT(78,7)=67.883000 !rg=4.95192 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(78,8,1:6)=(/ 56.753000,57.354000,57.894000,59.013000,60.505000,64.110000 /)
XPIZA_LKT(78,8,1:6)=(/ 0.726188,0.753717,0.768391,0.895152,0.925636,0.958659 /)
XCGA_LKT(78,8,1:6)=(/ 0.890863,0.885197,0.880257,0.843503,0.826947,0.779427 /)
XEXT_COEFF_550_LKT(78,8)=57.954000 !rg=4.95192 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(78,9,1:6)=(/ 48.190000,48.495000,49.139000,49.941000,51.143000,53.820000 /)
XPIZA_LKT(78,9,1:6)=(/ 0.715078,0.745505,0.764046,0.892051,0.918332,0.952763 /)
XCGA_LKT(78,9,1:6)=(/ 0.893353,0.886417,0.881550,0.846403,0.832750,0.790843 /)
XEXT_COEFF_550_LKT(78,9)=49.137000 !rg=4.95192 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(78,10,1:6)=(/ 40.607000,40.958000,41.407000,41.818000,42.908000,45.259000 /)
XPIZA_LKT(78,10,1:6)=(/ 0.701871,0.735939,0.758961,0.888859,0.912103,0.944377 /)
XCGA_LKT(78,10,1:6)=(/ 0.895947,0.888683,0.883043,0.849147,0.838277,0.799533 /)
XEXT_COEFF_550_LKT(78,10)=41.378000 !rg=4.95192 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(78,11,1:6)=(/ 34.059000,34.328000,34.647000,35.123000,36.030000,37.695000 /)
XPIZA_LKT(78,11,1:6)=(/ 0.688564,0.724851,0.750252,0.885067,0.906074,0.936566 /)
XCGA_LKT(78,11,1:6)=(/ 0.898683,0.890890,0.885033,0.851777,0.842150,0.806457 /)
XEXT_COEFF_550_LKT(78,11)=34.619000 !rg=4.95192 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(78,12,1:6)=(/ 28.507000,28.662000,28.941000,29.232000,29.776000,31.075000 /)
XPIZA_LKT(78,12,1:6)=(/ 0.675373,0.712413,0.741283,0.880633,0.901490,0.929444 /)
XCGA_LKT(78,12,1:6)=(/ 0.901893,0.893277,0.887087,0.854280,0.848447,0.817763 /)
XEXT_COEFF_550_LKT(78,12)=28.918000 !rg=4.95192 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(78,13,1:6)=(/ 23.769000,23.906000,24.055000,24.320000,24.781000,25.686000 /)
XPIZA_LKT(78,13,1:6)=(/ 0.661779,0.699329,0.730205,0.875216,0.896743,0.922834 /)
XCGA_LKT(78,13,1:6)=(/ 0.905280,0.896253,0.889183,0.856933,0.852280,0.825190 /)
XEXT_COEFF_550_LKT(78,13)=24.103000 !rg=4.95192 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(78,14,1:6)=(/ 19.805000,19.930000,20.044000,20.198000,20.549000,21.253000 /)
XPIZA_LKT(78,14,1:6)=(/ 0.648442,0.686046,0.718203,0.868036,0.891566,0.915817 /)
XCGA_LKT(78,14,1:6)=(/ 0.908747,0.899490,0.891947,0.859107,0.855700,0.832910 /)
XEXT_COEFF_550_LKT(78,14)=20.042000 !rg=4.95192 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(78,15,1:6)=(/ 16.496000,16.567000,16.680000,16.798000,17.071000,17.538000 /)
XPIZA_LKT(78,15,1:6)=(/ 0.636151,0.672243,0.705614,0.859983,0.886343,0.909897 /)
XCGA_LKT(78,15,1:6)=(/ 0.912467,0.902413,0.894763,0.861067,0.859383,0.839713 /)
XEXT_COEFF_550_LKT(78,15)=16.684000 !rg=4.95192 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(78,16,1:6)=(/ 13.698000,13.764000,13.831000,13.960000,14.141000,14.529000 /)
XPIZA_LKT(78,16,1:6)=(/ 0.624181,0.658840,0.691756,0.850476,0.879953,0.904130 /)
XCGA_LKT(78,16,1:6)=(/ 0.916147,0.906047,0.897770,0.863763,0.861380,0.842913 /)
XEXT_COEFF_550_LKT(78,16)=13.835000 !rg=4.95192 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(78,17,1:6)=(/ 11.394000,11.444000,11.491000,11.583000,11.715000,12.009000 /)
XPIZA_LKT(78,17,1:6)=(/ 0.613369,0.645930,0.678258,0.839226,0.872733,0.898927 /)
XCGA_LKT(78,17,1:6)=(/ 0.919900,0.909640,0.900850,0.866460,0.864390,0.848807 /)
XEXT_COEFF_550_LKT(78,17)=11.491000 !rg=4.95192 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(78,18,1:6)=(/ 9.483300,9.511900,9.560000,9.608700,9.717500,9.921400 /)
XPIZA_LKT(78,18,1:6)=(/ 0.603795,0.633536,0.665238,0.826492,0.864823,0.893348 /)
XCGA_LKT(78,18,1:6)=(/ 0.923660,0.913297,0.904417,0.868320,0.867430,0.853813 /)
XEXT_COEFF_550_LKT(78,18)=9.558600 !rg=4.95192 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(78,19,1:6)=(/ 7.875300,7.901300,7.928900,7.987200,8.047300,8.212900 /)
XPIZA_LKT(78,19,1:6)=(/ 0.594829,0.622229,0.651798,0.812781,0.854880,0.886675 /)
XCGA_LKT(78,19,1:6)=(/ 0.927260,0.917087,0.907943,0.871170,0.869100,0.855733 /)
XEXT_COEFF_550_LKT(78,19)=7.930700 !rg=4.95192 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(78,20,1:6)=(/ 6.544600,6.564500,6.584300,6.620800,6.676200,6.790100 /)
XPIZA_LKT(78,20,1:6)=(/ 0.586978,0.611709,0.639375,0.797295,0.844031,0.880416 /)
XCGA_LKT(78,20,1:6)=(/ 0.930700,0.920823,0.911503,0.874010,0.872183,0.860250 /)
XEXT_COEFF_550_LKT(78,20)=6.584200 !rg=4.95192 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(79,1,1:6)=(/ 114.540000,116.760000,119.110000,119.380000,136.400000,118.620000 /)
XPIZA_LKT(79,1,1:6)=(/ 0.769857,0.775021,0.765133,0.920185,0.957914,0.974015 /)
XCGA_LKT(79,1,1:6)=(/ 0.882340,0.880670,0.879650,0.821100,0.818470,0.672673 /)
XEXT_COEFF_550_LKT(79,1)=118.920000 !rg=5.36469 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(79,2,1:6)=(/ 110.170000,111.930000,113.130000,115.820000,121.070000,122.730000 /)
XPIZA_LKT(79,2,1:6)=(/ 0.768919,0.770263,0.760765,0.921297,0.955013,0.974789 /)
XCGA_LKT(79,2,1:6)=(/ 0.882463,0.882653,0.881200,0.827697,0.804790,0.703710 /)
XEXT_COEFF_550_LKT(79,2)=112.810000 !rg=5.36469 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(79,3,1:6)=(/ 101.870000,103.510000,104.600000,106.860000,110.340000,115.990000 /)
XPIZA_LKT(79,3,1:6)=(/ 0.763729,0.769537,0.763322,0.915464,0.952594,0.973704 /)
XCGA_LKT(79,3,1:6)=(/ 0.883643,0.881710,0.880750,0.830953,0.810980,0.727013 /)
XEXT_COEFF_550_LKT(79,3)=104.300000 !rg=5.36469 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(79,4,1:6)=(/ 92.079000,93.302000,94.407000,96.447000,99.339000,105.090000 /)
XPIZA_LKT(79,4,1:6)=(/ 0.758153,0.768361,0.767471,0.908823,0.948541,0.971664 /)
XCGA_LKT(79,4,1:6)=(/ 0.884790,0.882193,0.879770,0.833810,0.814297,0.743317 /)
XEXT_COEFF_550_LKT(79,4)=94.055000 !rg=5.36469 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(79,5,1:6)=(/ 81.445000,82.356000,83.573000,84.973000,87.642000,92.835000 /)
XPIZA_LKT(79,5,1:6)=(/ 0.750342,0.766520,0.770514,0.903536,0.943652,0.969473 /)
XCGA_LKT(79,5,1:6)=(/ 0.886130,0.882680,0.878763,0.837173,0.817553,0.758053 /)
XEXT_COEFF_550_LKT(79,5)=83.310000 !rg=5.36469 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(79,6,1:6)=(/ 71.126000,71.901000,72.951000,74.026000,76.633000,81.308000 /)
XPIZA_LKT(79,6,1:6)=(/ 0.741792,0.763107,0.770463,0.898869,0.936367,0.964846 /)
XCGA_LKT(79,6,1:6)=(/ 0.887863,0.883383,0.879867,0.840467,0.821273,0.766873 /)
XEXT_COEFF_550_LKT(79,6)=72.955000 !rg=5.36469 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(79,7,1:6)=(/ 61.358000,61.932000,62.589000,63.800000,65.532000,69.244000 /)
XPIZA_LKT(79,7,1:6)=(/ 0.731955,0.757265,0.768931,0.895295,0.928797,0.960776 /)
XCGA_LKT(79,7,1:6)=(/ 0.890077,0.884663,0.880090,0.843197,0.827073,0.778550 /)
XEXT_COEFF_550_LKT(79,7)=62.371000 !rg=5.36469 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(79,8,1:6)=(/ 52.423000,52.748000,53.410000,54.161000,55.720000,59.028000 /)
XPIZA_LKT(79,8,1:6)=(/ 0.721027,0.749894,0.766140,0.892644,0.921661,0.954869 /)
XCGA_LKT(79,8,1:6)=(/ 0.892160,0.885810,0.881253,0.845870,0.831707,0.786623 /)
XEXT_COEFF_550_LKT(79,8)=53.529000 !rg=5.36469 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(79,9,1:6)=(/ 44.333000,44.756000,45.135000,45.813000,46.936000,49.129000 /)
XPIZA_LKT(79,9,1:6)=(/ 0.708298,0.741612,0.761349,0.890031,0.915197,0.949681 /)
XCGA_LKT(79,9,1:6)=(/ 0.894453,0.887823,0.882097,0.848787,0.837667,0.799063 /)
XEXT_COEFF_550_LKT(79,9)=45.040000 !rg=5.36469 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(79,10,1:6)=(/ 37.449000,37.691000,38.085000,38.530000,39.449000,41.284000 /)
XPIZA_LKT(79,10,1:6)=(/ 0.695642,0.730699,0.754451,0.887017,0.909114,0.941528 /)
XCGA_LKT(79,10,1:6)=(/ 0.897217,0.889780,0.884080,0.850917,0.841550,0.805077 /)
XEXT_COEFF_550_LKT(79,10)=38.158000 !rg=5.36469 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(79,11,1:6)=(/ 31.428000,31.657000,31.937000,32.334000,33.023000,34.512000 /)
XPIZA_LKT(79,11,1:6)=(/ 0.682183,0.719293,0.746280,0.882517,0.903302,0.933034 /)
XCGA_LKT(79,11,1:6)=(/ 0.900230,0.892210,0.886150,0.852693,0.845367,0.812813 /)
XEXT_COEFF_550_LKT(79,11)=31.826000 !rg=5.36469 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(79,12,1:6)=(/ 26.284000,26.441000,26.615000,26.912000,27.445000,28.466000 /)
XPIZA_LKT(79,12,1:6)=(/ 0.668710,0.706664,0.736058,0.878141,0.898826,0.926537 /)
XCGA_LKT(79,12,1:6)=(/ 0.903510,0.894737,0.887913,0.855783,0.850327,0.821963 /)
XEXT_COEFF_550_LKT(79,12)=26.650000 !rg=5.36469 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(79,13,1:6)=(/ 21.917000,22.058000,22.196000,22.354000,22.808000,23.605000 /)
XPIZA_LKT(79,13,1:6)=(/ 0.655245,0.693374,0.724799,0.871757,0.894100,0.919331 /)
XCGA_LKT(79,13,1:6)=(/ 0.906897,0.897790,0.890580,0.857870,0.854460,0.829857 /)
XEXT_COEFF_550_LKT(79,13)=22.181000 !rg=5.36469 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(79,14,1:6)=(/ 18.283000,18.357000,18.497000,18.639000,18.922000,19.460000 /)
XPIZA_LKT(79,14,1:6)=(/ 0.642600,0.679406,0.712661,0.864657,0.888900,0.912763 /)
XCGA_LKT(79,14,1:6)=(/ 0.910547,0.900810,0.893483,0.860357,0.858023,0.837000 /)
XEXT_COEFF_550_LKT(79,14)=18.478000 !rg=5.36469 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(79,15,1:6)=(/ 15.213000,15.298000,15.358000,15.513000,15.740000,16.190000 /)
XPIZA_LKT(79,15,1:6)=(/ 0.630240,0.666313,0.699297,0.856046,0.883000,0.906384 /)
XCGA_LKT(79,15,1:6)=(/ 0.914183,0.904177,0.895853,0.862773,0.860253,0.840743 /)
XEXT_COEFF_550_LKT(79,15)=15.377000 !rg=5.36469 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(79,16,1:6)=(/ 12.637000,12.693000,12.760000,12.852000,13.027000,13.367000 /)
XPIZA_LKT(79,16,1:6)=(/ 0.618863,0.652737,0.685924,0.845320,0.876850,0.901445 /)
XCGA_LKT(79,16,1:6)=(/ 0.917947,0.907717,0.899213,0.864360,0.862640,0.845780 /)
XEXT_COEFF_550_LKT(79,16)=12.747000 !rg=5.36469 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(79,17,1:6)=(/ 10.515000,10.551000,10.606000,10.663000,10.799000,11.038000 /)
XPIZA_LKT(79,17,1:6)=(/ 0.608676,0.639997,0.672430,0.833671,0.869367,0.895883 /)
XCGA_LKT(79,17,1:6)=(/ 0.921700,0.911310,0.902510,0.867170,0.866077,0.851377 /)
XEXT_COEFF_550_LKT(79,17)=10.611000 !rg=5.36469 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(79,18,1:6)=(/ 8.747000,8.780500,8.808900,8.876600,8.959300,9.144200 /)
XPIZA_LKT(79,18,1:6)=(/ 0.599327,0.628436,0.658883,0.820706,0.860401,0.889887 /)
XCGA_LKT(79,18,1:6)=(/ 0.925333,0.915030,0.905903,0.869937,0.868200,0.854430 /)
XEXT_COEFF_550_LKT(79,18)=8.812100 !rg=5.36469 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(79,19,1:6)=(/ 7.265800,7.288000,7.315800,7.358400,7.427400,7.562100 /)
XPIZA_LKT(79,19,1:6)=(/ 0.590948,0.617144,0.646190,0.805904,0.850588,0.883885 /)
XCGA_LKT(79,19,1:6)=(/ 0.928870,0.918830,0.909583,0.872173,0.870403,0.857857 /)
XEXT_COEFF_550_LKT(79,19)=7.309600 !rg=5.36469 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(79,20,1:6)=(/ 6.039400,6.055200,6.076100,6.100600,6.160300,6.257900 /)
XPIZA_LKT(79,20,1:6)=(/ 0.583578,0.607000,0.633991,0.790092,0.839071,0.876945 /)
XCGA_LKT(79,20,1:6)=(/ 0.932270,0.922550,0.913303,0.875013,0.873377,0.862033 /)
XEXT_COEFF_550_LKT(79,20)=6.078200 !rg=5.36469 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(80,1,1:6)=(/ 105.630000,106.560000,108.740000,112.850000,112.340000,126.820000 /)
XPIZA_LKT(80,1,1:6)=(/ 0.766318,0.776465,0.745792,0.919315,0.950307,0.975977 /)
XCGA_LKT(80,1,1:6)=(/ 0.882440,0.880480,0.887053,0.834320,0.799720,0.742607 /)
XEXT_COEFF_550_LKT(80,1)=108.540000 !rg=5.81187 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(80,2,1:6)=(/ 101.270000,102.630000,104.350000,105.830000,110.710000,121.450000 /)
XPIZA_LKT(80,2,1:6)=(/ 0.763553,0.774276,0.750736,0.915189,0.952265,0.974198 /)
XCGA_LKT(80,2,1:6)=(/ 0.883563,0.881140,0.885483,0.830853,0.807447,0.741927 /)
XEXT_COEFF_550_LKT(80,2)=104.540000 !rg=5.81187 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(80,3,1:6)=(/ 94.087000,94.866000,96.522000,98.435000,102.550000,108.960000 /)
XPIZA_LKT(80,3,1:6)=(/ 0.759977,0.771937,0.759582,0.909995,0.945801,0.971958 /)
XCGA_LKT(80,3,1:6)=(/ 0.884473,0.881307,0.882640,0.835560,0.807283,0.742480 /)
XEXT_COEFF_550_LKT(80,3)=96.434000 !rg=5.81187 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(80,4,1:6)=(/ 84.948000,85.704000,87.159000,88.687000,92.391000,97.394000 /)
XPIZA_LKT(80,4,1:6)=(/ 0.753475,0.769398,0.765359,0.903509,0.941832,0.968887 /)
XCGA_LKT(80,4,1:6)=(/ 0.885677,0.881880,0.881153,0.838287,0.811217,0.749237 /)
XEXT_COEFF_550_LKT(80,4)=86.956000 !rg=5.81187 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(80,5,1:6)=(/ 75.184000,76.048000,76.816000,78.439000,80.548000,85.082000 /)
XPIZA_LKT(80,5,1:6)=(/ 0.745775,0.765712,0.768409,0.898475,0.938860,0.966784 /)
XCGA_LKT(80,5,1:6)=(/ 0.887370,0.882597,0.880630,0.840067,0.822183,0.767667 /)
XEXT_COEFF_550_LKT(80,5)=76.656000 !rg=5.81187 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(80,6,1:6)=(/ 65.640000,66.219000,66.996000,67.801000,70.697000,74.725000 /)
XPIZA_LKT(80,6,1:6)=(/ 0.736362,0.759879,0.769418,0.895862,0.932050,0.961956 /)
XCGA_LKT(80,6,1:6)=(/ 0.889157,0.884127,0.880630,0.843713,0.822793,0.768297 /)
XEXT_COEFF_550_LKT(80,6)=66.933000 !rg=5.81187 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(80,7,1:6)=(/ 56.442000,56.979000,57.757000,58.661000,60.698000,63.836000 /)
XPIZA_LKT(80,7,1:6)=(/ 0.725750,0.754115,0.767245,0.892398,0.923650,0.956961 /)
XCGA_LKT(80,7,1:6)=(/ 0.890963,0.885190,0.881467,0.845453,0.827317,0.778633 /)
XEXT_COEFF_550_LKT(80,7)=57.602000 !rg=5.81187 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(80,8,1:6)=(/ 48.294000,48.710000,49.093000,49.720000,51.218000,53.885000 /)
XPIZA_LKT(80,8,1:6)=(/ 0.714594,0.746195,0.763723,0.890190,0.918131,0.952787 /)
XCGA_LKT(80,8,1:6)=(/ 0.893577,0.886977,0.881917,0.848323,0.835093,0.791427 /)
XEXT_COEFF_550_LKT(80,8)=49.218000 !rg=5.81187 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(80,9,1:6)=(/ 40.944000,41.206000,41.594000,42.251000,43.160000,45.071000 /)
XPIZA_LKT(80,9,1:6)=(/ 0.702251,0.736723,0.758473,0.887897,0.911033,0.945119 /)
XCGA_LKT(80,9,1:6)=(/ 0.896073,0.888313,0.883537,0.849953,0.840930,0.804503 /)
XEXT_COEFF_550_LKT(80,9)=41.634000 !rg=5.81187 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(80,10,1:6)=(/ 34.535000,34.811000,35.032000,35.469000,36.360000,37.864000 /)
XPIZA_LKT(80,10,1:6)=(/ 0.688970,0.725938,0.750593,0.884469,0.905147,0.937301 /)
XCGA_LKT(80,10,1:6)=(/ 0.898700,0.891013,0.885213,0.852753,0.844403,0.811907 /)
XEXT_COEFF_550_LKT(80,10)=35.028000 !rg=5.81187 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(80,11,1:6)=(/ 29.008000,29.149000,29.417000,29.817000,30.317000,31.478000 /)
XPIZA_LKT(80,11,1:6)=(/ 0.675784,0.713238,0.741784,0.880258,0.900860,0.929911 /)
XCGA_LKT(80,11,1:6)=(/ 0.901887,0.893127,0.887223,0.854380,0.848930,0.818927 /)
XEXT_COEFF_550_LKT(80,11)=29.377000 !rg=5.81187 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(80,12,1:6)=(/ 24.239000,24.398000,24.554000,24.765000,25.252000,26.179000 /)
XPIZA_LKT(80,12,1:6)=(/ 0.662107,0.700667,0.731305,0.875204,0.895923,0.922291 /)
XCGA_LKT(80,12,1:6)=(/ 0.905123,0.896300,0.889347,0.857063,0.852473,0.826707 /)
XEXT_COEFF_550_LKT(80,12)=24.546000 !rg=5.81187 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(80,13,1:6)=(/ 20.224000,20.311000,20.480000,20.658000,20.992000,21.604000 /)
XPIZA_LKT(80,13,1:6)=(/ 0.649196,0.686779,0.719499,0.868800,0.891141,0.915601 /)
XCGA_LKT(80,13,1:6)=(/ 0.908633,0.899103,0.892090,0.859190,0.856377,0.834220 /)
XEXT_COEFF_550_LKT(80,13)=20.471000 !rg=5.81187 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(80,14,1:6)=(/ 16.851000,16.942000,17.023000,17.203000,17.424000,17.929000 /)
XPIZA_LKT(80,14,1:6)=(/ 0.636409,0.673465,0.706291,0.861434,0.886226,0.909650 /)
XCGA_LKT(80,14,1:6)=(/ 0.912240,0.902377,0.894577,0.861493,0.858897,0.838650 /)
XEXT_COEFF_550_LKT(80,14)=17.014000 !rg=5.81187 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(80,15,1:6)=(/ 14.034000,14.099000,14.183000,14.247000,14.460000,14.825000 /)
XPIZA_LKT(80,15,1:6)=(/ 0.624723,0.659830,0.693557,0.851516,0.880541,0.904137 /)
XCGA_LKT(80,15,1:6)=(/ 0.915927,0.905793,0.897533,0.863663,0.862047,0.843993 /)
XEXT_COEFF_550_LKT(80,15)=14.161000 !rg=5.81187 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(80,16,1:6)=(/ 11.662000,11.706000,11.764000,11.842000,11.989000,12.264000 /)
XPIZA_LKT(80,16,1:6)=(/ 0.613856,0.646648,0.679741,0.840277,0.873547,0.898984 /)
XCGA_LKT(80,16,1:6)=(/ 0.919777,0.909380,0.900743,0.865770,0.864670,0.849147 /)
XEXT_COEFF_550_LKT(80,16)=11.762000 !rg=5.81187 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(80,17,1:6)=(/ 9.699900,9.736900,9.774200,9.853500,9.914000,10.119000 /)
XPIZA_LKT(80,17,1:6)=(/ 0.603936,0.634405,0.666078,0.828331,0.864914,0.893131 /)
XCGA_LKT(80,17,1:6)=(/ 0.923433,0.913113,0.904067,0.868653,0.866283,0.852900 /)
XEXT_COEFF_550_LKT(80,17)=9.773200 !rg=5.81187 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(80,18,1:6)=(/ 8.070000,8.097300,8.132900,8.161000,8.251800,8.394500 /)
XPIZA_LKT(80,18,1:6)=(/ 0.595133,0.622937,0.653408,0.814106,0.856192,0.887431 /)
XCGA_LKT(80,18,1:6)=(/ 0.927017,0.916797,0.907587,0.870930,0.869213,0.856637 /)
XEXT_COEFF_550_LKT(80,18)=8.125500 !rg=5.81187 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(80,19,1:6)=(/ 6.704600,6.724000,6.746300,6.779700,6.846000,6.963800 /)
XPIZA_LKT(80,19,1:6)=(/ 0.587301,0.612232,0.640431,0.798981,0.845529,0.880901 /)
XCGA_LKT(80,19,1:6)=(/ 0.930483,0.920560,0.911307,0.873557,0.871813,0.859893 /)
XEXT_COEFF_550_LKT(80,19)=6.746800 !rg=5.81187 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(80,20,1:6)=(/ 5.572300,5.587500,5.603800,5.631800,5.664000,5.743900 /)
XPIZA_LKT(80,20,1:6)=(/ 0.580339,0.602616,0.628518,0.782986,0.832845,0.873046 /)
XCGA_LKT(80,20,1:6)=(/ 0.933757,0.924280,0.915040,0.876503,0.873933,0.862717 /)
XEXT_COEFF_550_LKT(80,20)=5.601900 !rg=5.81187 sigma=2.95 wvl=0.55
 
IF (LHOOK) CALL DR_HOOK('MODE_DUSTOPT:DUST_OPT_LKT_SET8',1,ZHOOK_HANDLE)
END SUBROUTINE DUST_OPT_LKT_SET8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE DUST_OPT_LKT_SET9()

  USE MODD_DUST_OPT_LKT
  
  IMPLICIT NONE

REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_DUSTOPT:DUST_OPT_LKT_SET9',0,ZHOOK_HANDLE)
XEXT_COEFF_WVL_LKT(81,1,1:6)=(/ 98.528000,99.670000,100.560000,101.780000,98.949000,132.390000 /)
XPIZA_LKT(81,1,1:6)=(/ 0.764207,0.777747,0.730800,0.911237,0.947466,0.975518 /)
XCGA_LKT(81,1,1:6)=(/ 0.884603,0.880620,0.892387,0.830657,0.794577,0.774727 /)
XEXT_COEFF_550_LKT(81,1)=99.359000 !rg=6.29632 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(81,2,1:6)=(/ 93.675000,94.724000,95.962000,97.386000,101.660000,115.450000 /)
XPIZA_LKT(81,2,1:6)=(/ 0.760202,0.774840,0.750235,0.908731,0.948202,0.973100 /)
XCGA_LKT(81,2,1:6)=(/ 0.884617,0.881323,0.886860,0.835417,0.809150,0.764460 /)
XEXT_COEFF_550_LKT(81,2)=95.589000 !rg=6.29632 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(81,3,1:6)=(/ 86.445000,87.830000,88.548000,90.763000,93.996000,102.030000 /)
XPIZA_LKT(81,3,1:6)=(/ 0.754182,0.771970,0.760822,0.902827,0.944738,0.968311 /)
XCGA_LKT(81,3,1:6)=(/ 0.885210,0.882100,0.882970,0.837000,0.815820,0.759407 /)
XEXT_COEFF_550_LKT(81,3)=88.606000 !rg=6.29632 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(81,4,1:6)=(/ 78.173000,79.179000,79.944000,81.761000,84.546000,90.360000 /)
XPIZA_LKT(81,4,1:6)=(/ 0.748076,0.767979,0.766451,0.897426,0.939990,0.965383 /)
XCGA_LKT(81,4,1:6)=(/ 0.886507,0.882970,0.881397,0.840130,0.819470,0.763560 /)
XEXT_COEFF_550_LKT(81,4)=79.795000 !rg=6.29632 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(81,5,1:6)=(/ 69.276000,69.865000,70.895000,72.075000,74.746000,78.768000 /)
XPIZA_LKT(81,5,1:6)=(/ 0.740220,0.763735,0.767556,0.894608,0.932588,0.962548 /)
XCGA_LKT(81,5,1:6)=(/ 0.888303,0.883277,0.881307,0.844133,0.820893,0.767830 /)
XEXT_COEFF_550_LKT(81,5)=70.867000 !rg=6.29632 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(81,6,1:6)=(/ 60.392000,61.068000,61.731000,62.904000,64.980000,68.577000 /)
XPIZA_LKT(81,6,1:6)=(/ 0.730280,0.758181,0.766845,0.892321,0.926847,0.958602 /)
XCGA_LKT(81,6,1:6)=(/ 0.889927,0.884393,0.881823,0.846633,0.826763,0.779843 /)
XEXT_COEFF_550_LKT(81,6)=61.779000 !rg=6.29632 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(81,7,1:6)=(/ 52.123000,52.563000,53.067000,54.022000,55.409000,58.586000 /)
XPIZA_LKT(81,7,1:6)=(/ 0.719914,0.750558,0.765504,0.890067,0.920154,0.953873 /)
XCGA_LKT(81,7,1:6)=(/ 0.892400,0.886230,0.881997,0.847563,0.832873,0.788400 /)
XEXT_COEFF_550_LKT(81,7)=52.867000 !rg=6.29632 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(81,8,1:6)=(/ 44.488000,44.960000,45.320000,45.806000,47.161000,49.380000 /)
XPIZA_LKT(81,8,1:6)=(/ 0.708028,0.742007,0.761165,0.887975,0.913885,0.948751 /)
XCGA_LKT(81,8,1:6)=(/ 0.894537,0.888200,0.883303,0.850340,0.838683,0.799913 /)
XEXT_COEFF_550_LKT(81,8)=45.335000 !rg=6.29632 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(81,9,1:6)=(/ 37.713000,37.974000,38.327000,38.840000,39.897000,41.719000 /)
XPIZA_LKT(81,9,1:6)=(/ 0.695481,0.731460,0.754993,0.885946,0.906843,0.939505 /)
XCGA_LKT(81,9,1:6)=(/ 0.897243,0.889670,0.884013,0.852383,0.842050,0.803710 /)
XEXT_COEFF_550_LKT(81,9)=38.382000 !rg=6.29632 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(81,10,1:6)=(/ 31.885000,32.031000,32.344000,32.665000,33.322000,34.623000 /)
XPIZA_LKT(81,10,1:6)=(/ 0.682626,0.719649,0.747088,0.882930,0.902779,0.933419 /)
XCGA_LKT(81,10,1:6)=(/ 0.900400,0.891947,0.886350,0.854017,0.847897,0.817390 /)
XEXT_COEFF_550_LKT(81,10)=32.358000 !rg=6.29632 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(81,11,1:6)=(/ 26.736000,26.895000,27.075000,27.424000,27.854000,28.783000 /)
XPIZA_LKT(81,11,1:6)=(/ 0.669060,0.707635,0.736931,0.878439,0.898677,0.927004 /)
XCGA_LKT(81,11,1:6)=(/ 0.903367,0.894640,0.887813,0.856267,0.852343,0.825140 /)
XEXT_COEFF_550_LKT(81,11)=27.077000 !rg=6.29632 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(81,12,1:6)=(/ 22.368000,22.468000,22.653000,22.826000,23.226000,23.948000 /)
XPIZA_LKT(81,12,1:6)=(/ 0.655846,0.693985,0.726179,0.872416,0.893541,0.919035 /)
XCGA_LKT(81,12,1:6)=(/ 0.906827,0.897487,0.890783,0.858173,0.855233,0.831403 /)
XEXT_COEFF_550_LKT(81,12)=22.620000 !rg=6.29632 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(81,13,1:6)=(/ 18.642000,18.744000,18.839000,19.044000,19.289000,19.852000 /)
XPIZA_LKT(81,13,1:6)=(/ 0.642717,0.680773,0.713507,0.865801,0.888528,0.912658 /)
XCGA_LKT(81,13,1:6)=(/ 0.910340,0.900600,0.893067,0.860477,0.857453,0.836447 /)
XEXT_COEFF_550_LKT(81,13)=18.829000 !rg=6.29632 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(81,14,1:6)=(/ 15.547000,15.611000,15.712000,15.826000,16.042000,16.472000 /)
XPIZA_LKT(81,14,1:6)=(/ 0.630723,0.666824,0.700776,0.856830,0.883679,0.907156 /)
XCGA_LKT(81,14,1:6)=(/ 0.913983,0.903873,0.895960,0.862227,0.861357,0.842570 /)
XEXT_COEFF_550_LKT(81,14)=15.701000 !rg=6.29632 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(81,15,1:6)=(/ 12.949000,12.994000,13.070000,13.158000,13.337000,13.658000 /)
XPIZA_LKT(81,15,1:6)=(/ 0.619428,0.653553,0.687140,0.846969,0.877411,0.901203 /)
XCGA_LKT(81,15,1:6)=(/ 0.917710,0.907383,0.898973,0.864690,0.863343,0.846680 /)
XEXT_COEFF_550_LKT(81,15)=13.075000 !rg=6.29632 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(81,16,1:6)=(/ 10.757000,10.800000,10.848000,10.926000,11.016000,11.246000 /)
XPIZA_LKT(81,16,1:6)=(/ 0.608900,0.640960,0.673426,0.835357,0.870025,0.896406 /)
XCGA_LKT(81,16,1:6)=(/ 0.921507,0.911087,0.902220,0.867393,0.866080,0.852380 /)
XEXT_COEFF_550_LKT(81,16)=10.836000 !rg=6.29632 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(81,17,1:6)=(/ 8.949700,8.978900,9.019900,9.070300,9.165200,9.344200 /)
XPIZA_LKT(81,17,1:6)=(/ 0.599554,0.628750,0.660150,0.822086,0.861518,0.890456 /)
XCGA_LKT(81,17,1:6)=(/ 0.925140,0.914783,0.905780,0.869697,0.868750,0.855890 /)
XEXT_COEFF_550_LKT(81,17)=9.014900 !rg=6.29632 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(81,18,1:6)=(/ 7.446500,7.466600,7.497100,7.536200,7.606800,7.742800 /)
XPIZA_LKT(81,18,1:6)=(/ 0.591216,0.617627,0.647230,0.807761,0.851619,0.884240 /)
XCGA_LKT(81,18,1:6)=(/ 0.928677,0.918547,0.909273,0.872163,0.870313,0.858383 /)
XEXT_COEFF_550_LKT(81,18)=7.501100 !rg=6.29632 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(81,19,1:6)=(/ 6.186200,6.204600,6.223500,6.255100,6.292300,6.388500 /)
XPIZA_LKT(81,19,1:6)=(/ 0.583754,0.607635,0.634828,0.792212,0.840185,0.877675 /)
XCGA_LKT(81,19,1:6)=(/ 0.932080,0.922300,0.912980,0.874937,0.873123,0.862277 /)
XEXT_COEFF_550_LKT(81,19)=6.220200 !rg=6.29632 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(81,20,1:6)=(/ 5.141900,5.154100,5.171000,5.188800,5.226300,5.302200 /)
XPIZA_LKT(81,20,1:6)=(/ 0.577256,0.598343,0.623415,0.775396,0.827446,0.869504 /)
XCGA_LKT(81,20,1:6)=(/ 0.935263,0.925903,0.916803,0.877827,0.875827,0.864947 /)
XEXT_COEFF_550_LKT(81,20)=5.168000 !rg=6.29632 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(82,1,1:6)=(/ 90.449000,91.209000,92.448000,93.310000,98.930000,121.320000 /)
XPIZA_LKT(82,1,1:6)=(/ 0.758046,0.776991,0.745618,0.907750,0.948667,0.975361 /)
XCGA_LKT(82,1,1:6)=(/ 0.886483,0.879573,0.887877,0.838320,0.818307,0.798257 /)
XEXT_COEFF_550_LKT(82,1)=91.446000 !rg=6.82116 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(82,2,1:6)=(/ 86.027000,86.962000,88.493000,90.260000,93.283000,106.470000 /)
XPIZA_LKT(82,2,1:6)=(/ 0.754138,0.775134,0.752191,0.902273,0.944008,0.968843 /)
XCGA_LKT(82,2,1:6)=(/ 0.885360,0.880197,0.886950,0.840813,0.814480,0.777823 /)
XEXT_COEFF_550_LKT(82,2)=88.159000 !rg=6.82116 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(82,3,1:6)=(/ 79.950000,80.533000,81.819000,83.339000,85.981000,92.024000 /)
XPIZA_LKT(82,3,1:6)=(/ 0.750181,0.772071,0.760812,0.896775,0.941541,0.968313 /)
XCGA_LKT(82,3,1:6)=(/ 0.886737,0.881703,0.883723,0.841367,0.819220,0.774963 /)
XEXT_COEFF_550_LKT(82,3)=81.684000 !rg=6.82116 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(82,4,1:6)=(/ 72.219000,72.774000,73.776000,75.213000,77.341000,81.641000 /)
XPIZA_LKT(82,4,1:6)=(/ 0.743332,0.767078,0.765065,0.892533,0.936397,0.965441 /)
XCGA_LKT(82,4,1:6)=(/ 0.888133,0.882490,0.882683,0.843663,0.822873,0.777520 /)
XEXT_COEFF_550_LKT(82,4)=73.737000 !rg=6.82116 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(82,5,1:6)=(/ 63.816000,64.588000,65.088000,66.301000,68.368000,72.597000 /)
XPIZA_LKT(82,5,1:6)=(/ 0.734093,0.761248,0.767997,0.890498,0.929389,0.959702 /)
XCGA_LKT(82,5,1:6)=(/ 0.889303,0.884253,0.881503,0.845290,0.827353,0.780067 /)
XEXT_COEFF_550_LKT(82,5)=65.088000 !rg=6.82116 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(82,6,1:6)=(/ 55.819000,56.255000,56.867000,57.757000,59.154000,62.086000 /)
XPIZA_LKT(82,6,1:6)=(/ 0.725016,0.754640,0.766756,0.889290,0.922997,0.955990 /)
XCGA_LKT(82,6,1:6)=(/ 0.891640,0.885467,0.881907,0.847863,0.832853,0.790857 /)
XEXT_COEFF_550_LKT(82,6)=56.892000 !rg=6.82116 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(82,7,1:6)=(/ 48.092000,48.334000,48.932000,49.831000,50.840000,53.166000 /)
XPIZA_LKT(82,7,1:6)=(/ 0.713740,0.746097,0.763270,0.887958,0.916116,0.951752 /)
XCGA_LKT(82,7,1:6)=(/ 0.893843,0.886440,0.883017,0.849113,0.837263,0.796900 /)
XEXT_COEFF_550_LKT(82,7)=48.795000 !rg=6.82116 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(82,8,1:6)=(/ 41.116000,41.326000,41.793000,42.338000,43.184000,44.990000 /)
XPIZA_LKT(82,8,1:6)=(/ 0.702062,0.736553,0.758238,0.886748,0.909598,0.944640 /)
XCGA_LKT(82,8,1:6)=(/ 0.896200,0.888777,0.884430,0.851977,0.842470,0.807033 /)
XEXT_COEFF_550_LKT(82,8)=41.738000 !rg=6.82116 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(82,9,1:6)=(/ 34.764000,35.045000,35.291000,35.672000,36.500000,38.188000 /)
XPIZA_LKT(82,9,1:6)=(/ 0.688438,0.726168,0.751243,0.883790,0.903971,0.936527 /)
XCGA_LKT(82,9,1:6)=(/ 0.898700,0.890913,0.885210,0.853150,0.845357,0.811450 /)
XEXT_COEFF_550_LKT(82,9)=35.264000 !rg=6.82116 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(82,10,1:6)=(/ 29.363000,29.616000,29.745000,30.167000,30.774000,31.934000 /)
XPIZA_LKT(82,10,1:6)=(/ 0.675516,0.714676,0.742056,0.880401,0.899423,0.928995 /)
XCGA_LKT(82,10,1:6)=(/ 0.901780,0.893437,0.887060,0.855510,0.849750,0.818460 /)
XEXT_COEFF_550_LKT(82,10)=29.823000 !rg=6.82116 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(82,11,1:6)=(/ 24.673000,24.786000,24.965000,25.219000,25.706000,26.581000 /)
XPIZA_LKT(82,11,1:6)=(/ 0.662725,0.701183,0.732225,0.875648,0.894954,0.922125 /)
XCGA_LKT(82,11,1:6)=(/ 0.905083,0.895893,0.889280,0.857037,0.853103,0.828040 /)
XEXT_COEFF_550_LKT(82,11)=24.972000 !rg=6.82116 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(82,12,1:6)=(/ 20.617000,20.732000,20.846000,21.087000,21.392000,22.068000 /)
XPIZA_LKT(82,12,1:6)=(/ 0.649363,0.687851,0.720202,0.869904,0.891137,0.915256 /)
XCGA_LKT(82,12,1:6)=(/ 0.908423,0.898990,0.891697,0.859457,0.856427,0.833227 /)
XEXT_COEFF_550_LKT(82,12)=20.854000 !rg=6.82116 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(82,13,1:6)=(/ 17.201000,17.281000,17.388000,17.512000,17.788000,18.286000 /)
XPIZA_LKT(82,13,1:6)=(/ 0.636762,0.674126,0.707967,0.861880,0.885944,0.909764 /)
XCGA_LKT(82,13,1:6)=(/ 0.912120,0.902110,0.894497,0.861303,0.859927,0.840727 /)
XEXT_COEFF_550_LKT(82,13)=17.376000 !rg=6.82116 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(82,14,1:6)=(/ 14.341000,14.402000,14.468000,14.594000,14.786000,15.164000 /)
XPIZA_LKT(82,14,1:6)=(/ 0.625075,0.660596,0.694330,0.852954,0.880914,0.904265 /)
XCGA_LKT(82,14,1:6)=(/ 0.915773,0.905513,0.897210,0.863743,0.862473,0.844653 /)
XEXT_COEFF_550_LKT(82,14)=14.487000 !rg=6.82116 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(82,15,1:6)=(/ 11.944000,11.994000,12.039000,12.130000,12.280000,12.581000 /)
XPIZA_LKT(82,15,1:6)=(/ 0.614176,0.647761,0.680806,0.841994,0.873731,0.898510 /)
XCGA_LKT(82,15,1:6)=(/ 0.919497,0.909073,0.900420,0.865820,0.864003,0.849070 /)
XEXT_COEFF_550_LKT(82,15)=12.036000 !rg=6.82116 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(82,16,1:6)=(/ 9.925000,9.957400,10.008000,10.066000,10.185000,10.404000 /)
XPIZA_LKT(82,16,1:6)=(/ 0.604190,0.634920,0.667493,0.829571,0.866182,0.893360 /)
XCGA_LKT(82,16,1:6)=(/ 0.923260,0.912830,0.903853,0.868520,0.866567,0.853080 /)
XEXT_COEFF_550_LKT(82,16)=9.995700 !rg=6.82116 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(82,17,1:6)=(/ 8.256800,8.281600,8.314300,8.359200,8.441000,8.607400 /)
XPIZA_LKT(82,17,1:6)=(/ 0.595299,0.623370,0.653955,0.816085,0.857427,0.888102 /)
XCGA_LKT(82,17,1:6)=(/ 0.926823,0.916527,0.907323,0.871133,0.869757,0.857007 /)
XEXT_COEFF_550_LKT(82,17)=8.316300 !rg=6.82116 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(82,18,1:6)=(/ 6.870700,6.891900,6.910400,6.948700,7.007600,7.130000 /)
XPIZA_LKT(82,18,1:6)=(/ 0.587465,0.612963,0.641220,0.800979,0.846555,0.881037 /)
XCGA_LKT(82,18,1:6)=(/ 0.930303,0.920277,0.910967,0.873323,0.871593,0.859700 /)
XEXT_COEFF_550_LKT(82,18)=6.910700 !rg=6.82116 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(82,19,1:6)=(/ 5.708300,5.722200,5.743100,5.765300,5.814400,5.907500 /)
XPIZA_LKT(82,19,1:6)=(/ 0.580479,0.603006,0.629537,0.784816,0.834867,0.874225 /)
XCGA_LKT(82,19,1:6)=(/ 0.933610,0.923970,0.914753,0.876370,0.873967,0.863287 /)
XEXT_COEFF_550_LKT(82,19)=5.736900 !rg=6.82116 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(82,20,1:6)=(/ 0.000000,4.755000,4.768700,4.788300,4.821900,4.889100 /)
XPIZA_LKT(82,20,1:6)=(/ 0.900,0.594268,0.618308,0.768095,0.821452,0.865808 /)
XCGA_LKT(82,20,1:6)=(/ 0.900,0.927553,0.918447,0.879327,0.876680,0.865860 /)
XEXT_COEFF_550_LKT(82,20)=4.768600 !rg=6.82116 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(83,1,1:6)=(/ 82.572000,83.550000,85.468000,86.000000,90.604000,98.295000 /)
XPIZA_LKT(83,1,1:6)=(/ 0.751040,0.776566,0.759142,0.898377,0.945383,0.971116 /)
XCGA_LKT(83,1,1:6)=(/ 0.885427,0.879603,0.885200,0.840620,0.828737,0.812767 /)
XEXT_COEFF_550_LKT(83,1)=83.806000 !rg=7.38975 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(83,2,1:6)=(/ 79.770000,80.556000,81.402000,83.085000,84.996000,90.864000 /)
XPIZA_LKT(83,2,1:6)=(/ 0.750563,0.774234,0.761586,0.894232,0.941109,0.967763 /)
XCGA_LKT(83,2,1:6)=(/ 0.887100,0.881953,0.884430,0.842873,0.821383,0.792993 /)
XEXT_COEFF_550_LKT(83,2)=81.524000 !rg=7.38975 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(83,3,1:6)=(/ 73.402000,74.226000,75.056000,76.084000,78.627000,83.574000 /)
XPIZA_LKT(83,3,1:6)=(/ 0.743567,0.770111,0.765590,0.890611,0.936973,0.966346 /)
XCGA_LKT(83,3,1:6)=(/ 0.887337,0.882293,0.882450,0.846033,0.824003,0.790710 /)
XEXT_COEFF_550_LKT(83,3)=74.999000 !rg=7.38975 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(83,4,1:6)=(/ 66.351000,67.060000,67.826000,68.723000,70.779000,74.633000 /)
XPIZA_LKT(83,4,1:6)=(/ 0.736814,0.765021,0.767828,0.888921,0.932302,0.963443 /)
XCGA_LKT(83,4,1:6)=(/ 0.888687,0.883323,0.881730,0.847190,0.828450,0.790147 /)
XEXT_COEFF_550_LKT(83,4)=67.649000 !rg=7.38975 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(83,5,1:6)=(/ 58.965000,59.327000,60.179000,61.165000,62.698000,65.768000 /)
XPIZA_LKT(83,5,1:6)=(/ 0.728755,0.758227,0.766764,0.887962,0.924886,0.959246 /)
XCGA_LKT(83,5,1:6)=(/ 0.890890,0.884497,0.882523,0.848007,0.831120,0.789867 /)
XEXT_COEFF_550_LKT(83,5)=60.055000 !rg=7.38975 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(83,6,1:6)=(/ 51.351000,51.892000,52.430000,53.172000,54.137000,56.847000 /)
XPIZA_LKT(83,6,1:6)=(/ 0.718318,0.750829,0.765479,0.887853,0.918331,0.954344 /)
XCGA_LKT(83,6,1:6)=(/ 0.892420,0.886540,0.882927,0.849633,0.837190,0.796080 /)
XEXT_COEFF_550_LKT(83,6)=52.363000 !rg=7.38975 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(83,7,1:6)=(/ 44.258000,44.604000,45.000000,45.685000,46.552000,48.555000 /)
XPIZA_LKT(83,7,1:6)=(/ 0.706836,0.741973,0.761056,0.886738,0.912722,0.949016 /)
XCGA_LKT(83,7,1:6)=(/ 0.894867,0.887880,0.882847,0.851397,0.842300,0.805013 /)
XEXT_COEFF_550_LKT(83,7)=45.004000 !rg=7.38975 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(83,8,1:6)=(/ 37.819000,38.171000,38.412000,38.974000,39.700000,41.343000 /)
XPIZA_LKT(83,8,1:6)=(/ 0.694519,0.732198,0.754602,0.885365,0.906102,0.940971 /)
XCGA_LKT(83,8,1:6)=(/ 0.897500,0.890047,0.884853,0.852727,0.844803,0.809487 /)
XEXT_COEFF_550_LKT(83,8)=38.359000 !rg=7.38975 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(83,9,1:6)=(/ 32.095000,32.264000,32.579000,32.999000,33.625000,34.902000 /)
XPIZA_LKT(83,9,1:6)=(/ 0.681973,0.720324,0.747131,0.882478,0.900635,0.933298 /)
XCGA_LKT(83,9,1:6)=(/ 0.900417,0.891970,0.886433,0.854453,0.848260,0.816977 /)
XEXT_COEFF_550_LKT(83,9)=32.599000 !rg=7.38975 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(83,10,1:6)=(/ 27.094000,27.245000,27.509000,27.695000,28.192000,29.306000 /)
XPIZA_LKT(83,10,1:6)=(/ 0.668965,0.708015,0.738313,0.878264,0.897199,0.925652 /)
XCGA_LKT(83,10,1:6)=(/ 0.903413,0.894510,0.888410,0.856453,0.852500,0.824037 /)
XEXT_COEFF_550_LKT(83,10)=27.462000 !rg=7.38975 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(83,11,1:6)=(/ 22.735000,22.865000,23.008000,23.247000,23.725000,24.564000 /)
XPIZA_LKT(83,11,1:6)=(/ 0.655859,0.694989,0.726747,0.873114,0.892786,0.918296 /)
XCGA_LKT(83,11,1:6)=(/ 0.906673,0.897403,0.890373,0.858447,0.854917,0.828990 /)
XEXT_COEFF_550_LKT(83,11)=23.011000 !rg=7.38975 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(83,12,1:6)=(/ 19.022000,19.100000,19.237000,19.402000,19.679000,20.260000 /)
XPIZA_LKT(83,12,1:6)=(/ 0.643198,0.681174,0.714817,0.866117,0.888867,0.912568 /)
XCGA_LKT(83,12,1:6)=(/ 0.910153,0.900363,0.893033,0.860357,0.858967,0.837517 /)
XEXT_COEFF_550_LKT(83,12)=19.222000 !rg=7.38975 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(83,13,1:6)=(/ 15.866000,15.936000,16.022000,16.153000,16.374000,16.827000 /)
XPIZA_LKT(83,13,1:6)=(/ 0.630835,0.667716,0.701736,0.857990,0.883771,0.907067 /)
XCGA_LKT(83,13,1:6)=(/ 0.913913,0.903640,0.895663,0.862927,0.861117,0.842533 /)
XEXT_COEFF_550_LKT(83,13)=16.045000 !rg=7.38975 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(83,14,1:6)=(/ 13.230000,13.288000,13.349000,13.416000,13.594000,13.942000 /)
XPIZA_LKT(83,14,1:6)=(/ 0.619604,0.654505,0.688225,0.847848,0.877816,0.901514 /)
XCGA_LKT(83,14,1:6)=(/ 0.917543,0.907227,0.898750,0.864437,0.863757,0.847563 /)
XEXT_COEFF_550_LKT(83,14)=13.345000 !rg=7.38975 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(83,15,1:6)=(/ 11.023000,11.057000,11.115000,11.176000,11.305000,11.547000 /)
XPIZA_LKT(83,15,1:6)=(/ 0.609365,0.641594,0.674958,0.836725,0.870987,0.896533 /)
XCGA_LKT(83,15,1:6)=(/ 0.921273,0.910730,0.901963,0.866873,0.866257,0.852673 /)
XEXT_COEFF_550_LKT(83,15)=11.118000 !rg=7.38975 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(83,16,1:6)=(/ 9.155500,9.187200,9.222300,9.288600,9.376300,9.579800 /)
XPIZA_LKT(83,16,1:6)=(/ 0.599642,0.629423,0.660965,0.823896,0.862332,0.890772 /)
XCGA_LKT(83,16,1:6)=(/ 0.924973,0.914547,0.905463,0.869587,0.868043,0.854773 /)
XEXT_COEFF_550_LKT(83,16)=9.220800 !rg=7.38975 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(83,17,1:6)=(/ 7.617500,7.642600,7.666500,7.711600,7.779900,7.925800 /)
XPIZA_LKT(83,17,1:6)=(/ 0.591269,0.618306,0.647916,0.809488,0.852658,0.884667 /)
XCGA_LKT(83,17,1:6)=(/ 0.928477,0.918273,0.908990,0.872080,0.870493,0.858630 /)
XEXT_COEFF_550_LKT(83,17)=7.667900 !rg=7.38975 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(83,18,1:6)=(/ 6.340200,6.356300,6.379400,6.405300,6.462000,6.558100 /)
XPIZA_LKT(83,18,1:6)=(/ 0.583955,0.608033,0.635960,0.793945,0.842025,0.878400 /)
XCGA_LKT(83,18,1:6)=(/ 0.931903,0.921997,0.912683,0.874660,0.873443,0.862157 /)
XEXT_COEFF_550_LKT(83,18)=6.377900 !rg=7.38975 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(83,19,1:6)=(/ 5.266700,5.279900,5.294900,5.321900,5.358500,5.430900 /)
XPIZA_LKT(83,19,1:6)=(/ 0.577370,0.598766,0.624050,0.777632,0.829103,0.870251 /)
XCGA_LKT(83,19,1:6)=(/ 0.935070,0.925667,0.916440,0.877600,0.875123,0.864000 /)
XEXT_COEFF_550_LKT(83,19)=5.295100 !rg=7.38975 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(83,20,1:6)=(/ 0.000000,4.388000,4.398700,4.417700,4.445200,4.502600 /)
XPIZA_LKT(83,20,1:6)=(/ 0.900,0.590511,0.613345,0.760255,0.814943,0.861454 /)
XCGA_LKT(83,20,1:6)=(/ 0.900,0.929187,0.920183,0.880657,0.877540,0.867343 /)
XEXT_COEFF_550_LKT(83,20)=4.399000 !rg=7.38975 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(84,1,1:6)=(/ 76.850000,77.201000,78.552000,79.438000,79.873000,80.824000 /)
XPIZA_LKT(84,1,1:6)=(/ 0.746676,0.772759,0.772437,0.891968,0.932713,0.962032 /)
XCGA_LKT(84,1,1:6)=(/ 0.887767,0.882443,0.881247,0.846247,0.813707,0.753587 /)
XEXT_COEFF_550_LKT(84,1)=78.024000 !rg=8.00573 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(84,2,1:6)=(/ 72.996000,73.963000,74.988000,76.132000,77.291000,80.850000 /)
XPIZA_LKT(84,2,1:6)=(/ 0.743130,0.770901,0.768019,0.886565,0.937875,0.964120 /)
XCGA_LKT(84,2,1:6)=(/ 0.887363,0.882500,0.882893,0.846003,0.829330,0.790543 /)
XEXT_COEFF_550_LKT(84,2)=74.884000 !rg=8.00573 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(84,3,1:6)=(/ 67.951000,68.484000,69.033000,70.324000,72.170000,75.851000 /)
XPIZA_LKT(84,3,1:6)=(/ 0.738572,0.767610,0.768400,0.884845,0.932178,0.963379 /)
XCGA_LKT(84,3,1:6)=(/ 0.889343,0.882450,0.882257,0.848030,0.829590,0.797920 /)
XEXT_COEFF_550_LKT(84,3)=69.301000 !rg=8.00573 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(84,4,1:6)=(/ 61.351000,61.851000,62.346000,63.634000,65.160000,68.311000 /)
XPIZA_LKT(84,4,1:6)=(/ 0.731097,0.761908,0.768464,0.885050,0.926744,0.960296 /)
XCGA_LKT(84,4,1:6)=(/ 0.890490,0.883623,0.882113,0.848643,0.832837,0.797480 /)
XEXT_COEFF_550_LKT(84,4)=62.387000 !rg=8.00573 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(84,5,1:6)=(/ 54.201000,54.786000,55.222000,55.954000,57.562000,60.099000 /)
XPIZA_LKT(84,5,1:6)=(/ 0.721455,0.755029,0.767104,0.886161,0.920366,0.956727 /)
XCGA_LKT(84,5,1:6)=(/ 0.891807,0.885580,0.882247,0.850550,0.835920,0.799783 /)
XEXT_COEFF_550_LKT(84,5)=55.113000 !rg=8.00573 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(84,6,1:6)=(/ 47.378000,47.771000,48.296000,48.856000,50.174000,52.518000 /)
XPIZA_LKT(84,6,1:6)=(/ 0.711662,0.746369,0.763297,0.885870,0.914030,0.950456 /)
XCGA_LKT(84,6,1:6)=(/ 0.894083,0.887043,0.883687,0.851917,0.840990,0.802533 /)
XEXT_COEFF_550_LKT(84,6)=48.305000 !rg=8.00573 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(84,7,1:6)=(/ 40.872000,41.094000,41.450000,42.091000,43.035000,44.734000 /)
XPIZA_LKT(84,7,1:6)=(/ 0.700440,0.736612,0.758301,0.884993,0.907102,0.944302 /)
XCGA_LKT(84,7,1:6)=(/ 0.896540,0.888817,0.884150,0.852430,0.844030,0.810160 /)
XEXT_COEFF_550_LKT(84,7)=41.482000 !rg=8.00573 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(84,8,1:6)=(/ 34.903000,35.108000,35.467000,35.855000,36.548000,38.139000 /)
XPIZA_LKT(84,8,1:6)=(/ 0.687870,0.726224,0.751824,0.883482,0.903023,0.936642 /)
XCGA_LKT(84,8,1:6)=(/ 0.899030,0.890863,0.885847,0.854213,0.848930,0.815100 /)
XEXT_COEFF_550_LKT(84,8)=35.470000 !rg=8.00573 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(84,9,1:6)=(/ 29.571000,29.797000,29.960000,30.325000,30.967000,31.977000 /)
XPIZA_LKT(84,9,1:6)=(/ 0.674887,0.714739,0.743010,0.880369,0.898202,0.929802 /)
XCGA_LKT(84,9,1:6)=(/ 0.901857,0.893403,0.887177,0.856300,0.850870,0.823377 /)
XEXT_COEFF_550_LKT(84,9)=29.910000 !rg=8.00573 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(84,10,1:6)=(/ 24.989000,25.125000,25.290000,25.572000,25.985000,26.870000 /)
XPIZA_LKT(84,10,1:6)=(/ 0.662358,0.701939,0.732871,0.876043,0.894638,0.922594 /)
XCGA_LKT(84,10,1:6)=(/ 0.904967,0.895893,0.889163,0.857930,0.854297,0.828273 /)
XEXT_COEFF_550_LKT(84,10)=25.307000 !rg=8.00573 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(84,11,1:6)=(/ 20.970000,21.089000,21.216000,21.406000,21.766000,22.472000 /)
XPIZA_LKT(84,11,1:6)=(/ 0.649256,0.688454,0.721476,0.869803,0.890377,0.915122 /)
XCGA_LKT(84,11,1:6)=(/ 0.908387,0.898840,0.891670,0.859237,0.856457,0.833673 /)
XEXT_COEFF_550_LKT(84,11)=21.179000 !rg=8.00573 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(84,12,1:6)=(/ 17.543000,17.627000,17.711000,17.866000,18.140000,18.642000 /)
XPIZA_LKT(84,12,1:6)=(/ 0.636916,0.674869,0.708585,0.862843,0.886021,0.909648 /)
XCGA_LKT(84,12,1:6)=(/ 0.911997,0.901950,0.894170,0.861707,0.859767,0.840267 /)
XEXT_COEFF_550_LKT(84,12)=17.735000 !rg=8.00573 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(84,13,1:6)=(/ 14.632000,14.701000,14.768000,14.873000,15.084000,15.483000 /)
XPIZA_LKT(84,13,1:6)=(/ 0.625047,0.661416,0.695410,0.853701,0.881056,0.904108 /)
XCGA_LKT(84,13,1:6)=(/ 0.915637,0.905367,0.897007,0.863810,0.862737,0.845687 /)
XEXT_COEFF_550_LKT(84,13)=14.767000 !rg=8.00573 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(84,14,1:6)=(/ 12.209000,12.252000,12.319000,12.387000,12.556000,12.811000 /)
XPIZA_LKT(84,14,1:6)=(/ 0.614403,0.648231,0.682176,0.843008,0.874588,0.898626 /)
XCGA_LKT(84,14,1:6)=(/ 0.919387,0.908870,0.900333,0.865717,0.864870,0.850547 /)
XEXT_COEFF_550_LKT(84,14)=12.315000 !rg=8.00573 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(84,15,1:6)=(/ 10.166000,10.208000,10.241000,10.326000,10.424000,10.655000 /)
XPIZA_LKT(84,15,1:6)=(/ 0.604429,0.635980,0.668365,0.831329,0.866984,0.893101 /)
XCGA_LKT(84,15,1:6)=(/ 0.923030,0.912553,0.903437,0.868407,0.867163,0.853443 /)
XEXT_COEFF_550_LKT(84,15)=10.243000 !rg=8.00573 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(84,16,1:6)=(/ 8.447500,8.473500,8.510400,8.548000,8.640100,8.810500 /)
XPIZA_LKT(84,16,1:6)=(/ 0.595409,0.623887,0.655005,0.817187,0.858150,0.887872 /)
XCGA_LKT(84,16,1:6)=(/ 0.926663,0.916303,0.907143,0.870217,0.868953,0.856373 /)
XEXT_COEFF_550_LKT(84,16)=8.502600 !rg=8.00573 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(84,17,1:6)=(/ 7.030100,7.049800,7.075600,7.101500,7.177100,7.293700 /)
XPIZA_LKT(84,17,1:6)=(/ 0.587557,0.613221,0.642330,0.802578,0.848346,0.881862 /)
XCGA_LKT(84,17,1:6)=(/ 0.930153,0.920067,0.910677,0.873087,0.872087,0.860863 /)
XEXT_COEFF_550_LKT(84,17)=7.076400 !rg=8.00573 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(84,18,1:6)=(/ 5.849300,5.866100,5.881800,5.916800,5.954400,6.048400 /)
XPIZA_LKT(84,18,1:6)=(/ 0.580570,0.603607,0.630187,0.786932,0.836175,0.874665 /)
XCGA_LKT(84,18,1:6)=(/ 0.933430,0.923727,0.914380,0.876123,0.873750,0.863320 /)
XEXT_COEFF_550_LKT(84,18)=5.881600 !rg=8.00573 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(84,19,1:6)=(/ 4.859800,4.871200,4.885600,4.905000,4.940700,5.010700 /)
XPIZA_LKT(84,19,1:6)=(/ 0.574439,0.594647,0.619039,0.769762,0.823077,0.866728 /)
XCGA_LKT(84,19,1:6)=(/ 0.936520,0.927343,0.918190,0.878703,0.876193,0.865357 /)
XEXT_COEFF_550_LKT(84,19)=4.884100 !rg=8.00573 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(84,20,1:6)=(/ 0.000000,4.048400,4.059300,4.070500,4.101600,4.151900 /)
XPIZA_LKT(84,20,1:6)=(/ 0.900,0.586850,0.608801,0.752260,0.808681,0.857356 /)
XCGA_LKT(84,20,1:6)=(/ 0.900,0.930790,0.921897,0.881993,0.878940,0.868793 /)
XEXT_COEFF_550_LKT(84,20)=4.059300 !rg=8.00573 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(85,1,1:6)=(/ 70.769000,71.365000,71.586000,72.455000,76.049000,71.081000 /)
XPIZA_LKT(85,1,1:6)=(/ 0.741324,0.770455,0.776832,0.882355,0.935913,0.954159 /)
XCGA_LKT(85,1,1:6)=(/ 0.889457,0.882953,0.879337,0.849100,0.833303,0.769713 /)
XEXT_COEFF_550_LKT(85,1)=71.346000 !rg=8.67306 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(85,2,1:6)=(/ 67.670000,68.236000,69.146000,70.237000,71.949000,74.118000 /)
XPIZA_LKT(85,2,1:6)=(/ 0.738310,0.767802,0.772989,0.879289,0.932318,0.959630 /)
XCGA_LKT(85,2,1:6)=(/ 0.889580,0.883007,0.881827,0.853070,0.830723,0.796017 /)
XEXT_COEFF_550_LKT(85,2)=69.220000 !rg=8.67306 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(85,3,1:6)=(/ 62.504000,63.031000,63.782000,64.750000,66.791000,70.238000 /)
XPIZA_LKT(85,3,1:6)=(/ 0.731538,0.763067,0.771152,0.880313,0.926900,0.957764 /)
XCGA_LKT(85,3,1:6)=(/ 0.890137,0.883573,0.881177,0.852083,0.830770,0.788633 /)
XEXT_COEFF_550_LKT(85,3)=63.898000 !rg=8.67306 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(85,4,1:6)=(/ 56.482000,56.908000,57.674000,58.430000,60.319000,63.381000 /)
XPIZA_LKT(85,4,1:6)=(/ 0.724195,0.757331,0.769045,0.882968,0.920753,0.954444 /)
XCGA_LKT(85,4,1:6)=(/ 0.891550,0.884913,0.881883,0.852080,0.833700,0.791860 /)
XEXT_COEFF_550_LKT(85,4)=57.665000 !rg=8.67306 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(85,5,1:6)=(/ 50.116000,50.408000,50.872000,51.735000,52.929000,55.087000 /)
XPIZA_LKT(85,5,1:6)=(/ 0.715362,0.750346,0.766144,0.884505,0.914849,0.953019 /)
XCGA_LKT(85,5,1:6)=(/ 0.893637,0.885900,0.882870,0.851543,0.839837,0.806000 /)
XEXT_COEFF_550_LKT(85,5)=51.001000 !rg=8.67306 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(85,6,1:6)=(/ 43.698000,44.017000,44.437000,44.877000,46.211000,48.460000 /)
XPIZA_LKT(85,6,1:6)=(/ 0.704692,0.741150,0.761518,0.885175,0.909701,0.946683 /)
XCGA_LKT(85,6,1:6)=(/ 0.895437,0.888157,0.883857,0.853360,0.842937,0.802800 /)
XEXT_COEFF_550_LKT(85,6)=44.428000 !rg=8.67306 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(85,7,1:6)=(/ 37.634000,37.888000,38.218000,38.747000,39.792000,41.465000 /)
XPIZA_LKT(85,7,1:6)=(/ 0.693059,0.731465,0.755111,0.884052,0.903265,0.938878 /)
XCGA_LKT(85,7,1:6)=(/ 0.897860,0.889970,0.884947,0.854253,0.845753,0.809140 /)
XEXT_COEFF_550_LKT(85,7)=38.238000 !rg=8.67306 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(85,8,1:6)=(/ 32.188000,32.399000,32.608000,33.000000,33.663000,34.981000 /)
XPIZA_LKT(85,8,1:6)=(/ 0.681031,0.720344,0.747822,0.881805,0.899932,0.933685 /)
XCGA_LKT(85,8,1:6)=(/ 0.900593,0.892147,0.886327,0.855623,0.850673,0.819010 /)
XEXT_COEFF_550_LKT(85,8)=32.689000 !rg=8.67306 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(85,9,1:6)=(/ 27.295000,27.433000,27.667000,27.928000,28.430000,29.364000 /)
XPIZA_LKT(85,9,1:6)=(/ 0.668272,0.708065,0.738799,0.878330,0.895885,0.925304 /)
XCGA_LKT(85,9,1:6)=(/ 0.903613,0.894450,0.888400,0.856743,0.853960,0.827820 /)
XEXT_COEFF_550_LKT(85,9)=27.685000 !rg=8.67306 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(85,10,1:6)=(/ 23.042000,23.174000,23.320000,23.534000,23.942000,24.738000 /)
XPIZA_LKT(85,10,1:6)=(/ 0.655484,0.695511,0.727880,0.873194,0.891849,0.918079 /)
XCGA_LKT(85,10,1:6)=(/ 0.906730,0.897437,0.890463,0.858537,0.855847,0.832537 /)
XEXT_COEFF_550_LKT(85,10)=23.298000 !rg=8.67306 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(85,11,1:6)=(/ 19.350000,19.421000,19.570000,19.785000,20.038000,20.606000 /)
XPIZA_LKT(85,11,1:6)=(/ 0.643068,0.681659,0.715690,0.867203,0.888490,0.912321 /)
XCGA_LKT(85,11,1:6)=(/ 0.910167,0.900117,0.892937,0.860843,0.859300,0.838287 /)
XEXT_COEFF_550_LKT(85,11)=19.594000 !rg=8.67306 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(85,12,1:6)=(/ 16.180000,16.270000,16.346000,16.437000,16.681000,17.157000 /)
XPIZA_LKT(85,12,1:6)=(/ 0.630906,0.668559,0.702818,0.858677,0.883673,0.906337 /)
XCGA_LKT(85,12,1:6)=(/ 0.913690,0.903687,0.895597,0.862503,0.861557,0.843443 /)
XEXT_COEFF_550_LKT(85,12)=16.340000 !rg=8.67306 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(85,13,1:6)=(/ 13.503000,13.553000,13.628000,13.718000,13.913000,14.212000 /)
XPIZA_LKT(85,13,1:6)=(/ 0.619685,0.654898,0.689438,0.849304,0.878143,0.900805 /)
XCGA_LKT(85,13,1:6)=(/ 0.917440,0.906980,0.898573,0.864753,0.863920,0.848477 /)
XEXT_COEFF_550_LKT(85,13)=13.634000 !rg=8.67306 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(85,14,1:6)=(/ 11.260000,11.306000,11.351000,11.436000,11.550000,11.793000 /)
XPIZA_LKT(85,14,1:6)=(/ 0.609222,0.642307,0.675681,0.838321,0.871218,0.896055 /)
XCGA_LKT(85,14,1:6)=(/ 0.921157,0.910647,0.901727,0.867027,0.865580,0.851877 /)
XEXT_COEFF_550_LKT(85,14)=11.340000 !rg=8.67306 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(85,15,1:6)=(/ 9.378400,9.411100,9.458200,9.497600,9.598000,9.783400 /)
XPIZA_LKT(85,15,1:6)=(/ 0.599815,0.630115,0.662556,0.825352,0.863598,0.891047 /)
XCGA_LKT(85,15,1:6)=(/ 0.924763,0.914250,0.905190,0.869387,0.868293,0.855823 /)
XEXT_COEFF_550_LKT(85,15)=9.447600 !rg=8.67306 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(85,16,1:6)=(/ 7.795000,7.815900,7.846700,7.887700,7.960200,8.102100 /)
XPIZA_LKT(85,16,1:6)=(/ 0.591346,0.618560,0.648875,0.810953,0.854237,0.885201 /)
XCGA_LKT(85,16,1:6)=(/ 0.928383,0.918043,0.908790,0.871617,0.870793,0.859073 /)
XEXT_COEFF_550_LKT(85,16)=7.845700 !rg=8.67306 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(85,17,1:6)=(/ 6.485800,6.504800,6.524600,6.562400,6.590400,6.691600 /)
XPIZA_LKT(85,17,1:6)=(/ 0.583904,0.608503,0.636487,0.795978,0.842725,0.878253 /)
XCGA_LKT(85,17,1:6)=(/ 0.931763,0.921793,0.912420,0.874617,0.872507,0.861480 /)
XEXT_COEFF_550_LKT(85,17)=6.524900 !rg=8.67306 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(85,18,1:6)=(/ 5.397300,5.410600,5.429400,5.445600,5.489700,5.565800 /)
XPIZA_LKT(85,18,1:6)=(/ 0.577427,0.599130,0.625075,0.779341,0.830806,0.871268 /)
XCGA_LKT(85,18,1:6)=(/ 0.934940,0.925420,0.916163,0.877227,0.875250,0.864390 /)
XEXT_COEFF_550_LKT(85,18)=5.425300 !rg=8.67306 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(85,19,1:6)=(/ 0.000000,4.494100,4.506700,4.522600,4.555800,4.615800 /)
XPIZA_LKT(85,19,1:6)=(/ 0.900,0.590750,0.614062,0.761997,0.816790,0.862615 /)
XCGA_LKT(85,19,1:6)=(/ 0.900,0.928953,0.919927,0.880217,0.877490,0.867103 /)
XEXT_COEFF_550_LKT(85,19)=4.506000 !rg=8.67306 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(85,20,1:6)=(/ 0.000000,3.735600,3.744300,3.759700,3.772700,3.814300 /)
XPIZA_LKT(85,20,1:6)=(/ 0.900,0.583472,0.604242,0.744724,0.801508,0.852207 /)
XCGA_LKT(85,20,1:6)=(/ 0.900,0.932327,0.923607,0.883657,0.879767,0.869407 /)
XEXT_COEFF_550_LKT(85,20)=3.743800 !rg=8.67306 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(86,1,1:6)=(/ 64.909000,65.343000,66.106000,67.692000,70.067000,69.986000 /)
XPIZA_LKT(86,1,1:6)=(/ 0.734050,0.765440,0.776817,0.868619,0.929981,0.959728 /)
XCGA_LKT(86,1,1:6)=(/ 0.890627,0.882803,0.880603,0.856283,0.832260,0.798617 /)
XEXT_COEFF_550_LKT(86,1)=66.115000 !rg=9.39601 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(86,2,1:6)=(/ 62.212000,62.792000,63.425000,64.256000,66.520000,69.457000 /)
XPIZA_LKT(86,2,1:6)=(/ 0.731013,0.764310,0.775297,0.874004,0.928078,0.957436 /)
XCGA_LKT(86,2,1:6)=(/ 0.890393,0.883860,0.881090,0.855110,0.833837,0.789170 /)
XEXT_COEFF_550_LKT(86,2)=63.181000 !rg=9.39601 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(86,3,1:6)=(/ 57.648000,58.235000,58.648000,59.363000,61.163000,64.404000 /)
XPIZA_LKT(86,3,1:6)=(/ 0.725199,0.759505,0.772201,0.878627,0.922075,0.955616 /)
XCGA_LKT(86,3,1:6)=(/ 0.891560,0.884777,0.881223,0.852690,0.834947,0.800233 /)
XEXT_COEFF_550_LKT(86,3)=58.621000 !rg=9.39601 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(86,4,1:6)=(/ 52.079000,52.575000,52.934000,53.683000,55.060000,58.119000 /)
XPIZA_LKT(86,4,1:6)=(/ 0.717187,0.753242,0.768751,0.881308,0.916024,0.952404 /)
XCGA_LKT(86,4,1:6)=(/ 0.892860,0.885960,0.882167,0.852507,0.838317,0.802587 /)
XEXT_COEFF_550_LKT(86,4)=52.942000 !rg=9.39601 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(86,5,1:6)=(/ 46.143000,46.488000,46.922000,47.537000,48.908000,51.137000 /)
XPIZA_LKT(86,5,1:6)=(/ 0.707979,0.745243,0.764687,0.883473,0.909307,0.947045 /)
XCGA_LKT(86,5,1:6)=(/ 0.894750,0.887357,0.882767,0.853773,0.841670,0.802570 /)
XEXT_COEFF_550_LKT(86,5)=47.037000 !rg=9.39601 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(86,6,1:6)=(/ 40.303000,40.616000,40.952000,41.632000,42.561000,44.494000 /)
XPIZA_LKT(86,6,1:6)=(/ 0.697479,0.736358,0.758331,0.884182,0.904551,0.942110 /)
XCGA_LKT(86,6,1:6)=(/ 0.896923,0.889143,0.884630,0.854690,0.845310,0.809547 /)
XEXT_COEFF_550_LKT(86,6)=40.928000 !rg=9.39601 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(86,7,1:6)=(/ 34.707000,34.954000,35.193000,35.617000,36.346000,37.906000 /)
XPIZA_LKT(86,7,1:6)=(/ 0.685755,0.725686,0.751944,0.882264,0.900163,0.936000 /)
XCGA_LKT(86,7,1:6)=(/ 0.899403,0.891180,0.885663,0.854663,0.848683,0.817127 /)
XEXT_COEFF_550_LKT(86,7)=35.120000 !rg=9.39601 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(86,8,1:6)=(/ 29.673000,29.900000,30.099000,30.327000,31.008000,32.168000 /)
XPIZA_LKT(86,8,1:6)=(/ 0.673745,0.714443,0.743485,0.879642,0.896964,0.929164 /)
XCGA_LKT(86,8,1:6)=(/ 0.902147,0.893643,0.887410,0.856610,0.853110,0.824173 /)
XEXT_COEFF_550_LKT(86,8)=30.067000 !rg=9.39601 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(86,9,1:6)=(/ 25.161000,25.306000,25.460000,25.754000,26.249000,27.169000 /)
XPIZA_LKT(86,9,1:6)=(/ 0.661148,0.701675,0.733336,0.876110,0.892749,0.920319 /)
XCGA_LKT(86,9,1:6)=(/ 0.905220,0.896063,0.889073,0.858413,0.855180,0.828160 /)
XEXT_COEFF_550_LKT(86,9)=25.508000 !rg=9.39601 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(86,10,1:6)=(/ 21.265000,21.355000,21.505000,21.666000,22.015000,22.644000 /)
XPIZA_LKT(86,10,1:6)=(/ 0.649092,0.688584,0.722412,0.870467,0.890026,0.914874 /)
XCGA_LKT(86,10,1:6)=(/ 0.908477,0.898763,0.891793,0.859737,0.858510,0.836937 /)
XEXT_COEFF_550_LKT(86,10)=21.489000 !rg=9.39601 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(86,11,1:6)=(/ 17.843000,17.935000,18.007000,18.217000,18.443000,18.893000 /)
XPIZA_LKT(86,11,1:6)=(/ 0.636675,0.675451,0.709434,0.863901,0.886313,0.909683 /)
XCGA_LKT(86,11,1:6)=(/ 0.911957,0.901967,0.893847,0.862010,0.860903,0.842810 /)
XEXT_COEFF_550_LKT(86,11)=18.017000 !rg=9.39601 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(86,12,1:6)=(/ 14.933000,14.990000,15.088000,15.187000,15.384000,15.734000 /)
XPIZA_LKT(86,12,1:6)=(/ 0.625256,0.661859,0.696804,0.854721,0.881042,0.903361 /)
XCGA_LKT(86,12,1:6)=(/ 0.915533,0.905170,0.897173,0.863800,0.863367,0.846763 /)
XEXT_COEFF_550_LKT(86,12)=15.073000 !rg=9.39601 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(86,13,1:6)=(/ 12.454000,12.507000,12.556000,12.659000,12.769000,13.059000 /)
XPIZA_LKT(86,13,1:6)=(/ 0.614184,0.648770,0.682833,0.844690,0.874813,0.898226 /)
XCGA_LKT(86,13,1:6)=(/ 0.919253,0.908740,0.899970,0.866123,0.864537,0.850153 /)
XEXT_COEFF_550_LKT(86,13)=12.549000 !rg=9.39601 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(86,14,1:6)=(/ 10.388000,10.423000,10.476000,10.524000,10.646000,10.857000 /)
XPIZA_LKT(86,14,1:6)=(/ 0.604395,0.636182,0.669614,0.832530,0.867819,0.893697 /)
XCGA_LKT(86,14,1:6)=(/ 0.922920,0.912317,0.903397,0.867757,0.867393,0.855083 /)
XEXT_COEFF_550_LKT(86,14)=10.461000 !rg=9.39601 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(86,15,1:6)=(/ 8.653200,8.679300,8.715000,8.763900,8.847400,9.013900 /)
XPIZA_LKT(86,15,1:6)=(/ 0.595459,0.624458,0.656094,0.819330,0.859431,0.888156 /)
XCGA_LKT(86,15,1:6)=(/ 0.926513,0.916007,0.906737,0.870533,0.869377,0.857147 /)
XEXT_COEFF_550_LKT(86,15)=8.720100 !rg=9.39601 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(86,16,1:6)=(/ 7.192100,7.212300,7.234100,7.277400,7.328300,7.438500 /)
XPIZA_LKT(86,16,1:6)=(/ 0.587530,0.613602,0.642853,0.804547,0.849210,0.882233 /)
XCGA_LKT(86,16,1:6)=(/ 0.930027,0.919817,0.910367,0.873077,0.871667,0.860793 /)
XEXT_COEFF_550_LKT(86,16)=7.233600 !rg=9.39601 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(86,17,1:6)=(/ 5.984700,5.999700,6.020300,6.045500,6.095000,6.184800 /)
XPIZA_LKT(86,17,1:6)=(/ 0.580561,0.603773,0.631034,0.788499,0.837970,0.875684 /)
XCGA_LKT(86,17,1:6)=(/ 0.933317,0.923487,0.914183,0.875787,0.874493,0.864297 /)
XEXT_COEFF_550_LKT(86,17)=6.015600 !rg=9.39601 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(86,18,1:6)=(/ 4.980400,4.991600,5.006500,5.025800,5.060900,5.128200 /)
XPIZA_LKT(86,18,1:6)=(/ 0.574527,0.594961,0.619718,0.771759,0.824758,0.867549 /)
XCGA_LKT(86,18,1:6)=(/ 0.936370,0.927090,0.917873,0.878650,0.876360,0.865690 /)
XEXT_COEFF_550_LKT(86,18)=5.005100 !rg=9.39601 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(86,19,1:6)=(/ 0.000000,4.147100,4.157000,4.174200,4.194300,4.240700 /)
XPIZA_LKT(86,19,1:6)=(/ 0.900,0.587139,0.609283,0.754472,0.810063,0.858274 /)
XCGA_LKT(86,19,1:6)=(/ 0.900,0.930553,0.921597,0.881837,0.878440,0.868643 /)
XEXT_COEFF_550_LKT(86,19)=4.156100 !rg=9.39601 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(86,20,1:6)=(/ 0.000000,3.446600,3.455000,3.465600,3.485400,3.523300 /)
XPIZA_LKT(86,20,1:6)=(/ 0.900,0.580162,0.600041,0.736500,0.794875,0.847942 /)
XCGA_LKT(86,20,1:6)=(/ 0.900,0.933837,0.925277,0.885173,0.881457,0.871540 /)
XEXT_COEFF_550_LKT(86,20)=3.453500 !rg=9.39601 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(87,1,1:6)=(/ 59.702000,60.722000,60.959000,62.295000,62.369000,71.927000 /)
XPIZA_LKT(87,1,1:6)=(/ 0.726686,0.763454,0.778354,0.859716,0.924199,0.961031 /)
XCGA_LKT(87,1,1:6)=(/ 0.889877,0.883893,0.878713,0.857193,0.830993,0.815000 /)
XEXT_COEFF_550_LKT(87,1)=61.203000 !rg=10.1792 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(87,2,1:6)=(/ 57.393000,57.822000,58.458000,59.678000,61.125000,64.686000 /)
XPIZA_LKT(87,2,1:6)=(/ 0.723969,0.759534,0.774606,0.871878,0.922478,0.953199 /)
XCGA_LKT(87,2,1:6)=(/ 0.891703,0.884680,0.881367,0.856120,0.834770,0.800000 /)
XEXT_COEFF_550_LKT(87,2)=58.697000 !rg=10.1792 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(87,3,1:6)=(/ 53.241000,53.545000,54.234000,54.866000,56.149000,58.569000 /)
XPIZA_LKT(87,3,1:6)=(/ 0.718612,0.754884,0.771659,0.878748,0.916935,0.955609 /)
XCGA_LKT(87,3,1:6)=(/ 0.893040,0.885773,0.882240,0.854123,0.839147,0.805800 /)
XEXT_COEFF_550_LKT(87,3)=54.076000 !rg=10.1792 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(87,4,1:6)=(/ 48.107000,48.338000,48.967000,49.644000,50.684000,52.831000 /)
XPIZA_LKT(87,4,1:6)=(/ 0.710657,0.747751,0.767360,0.881472,0.910264,0.951333 /)
XCGA_LKT(87,4,1:6)=(/ 0.894643,0.886677,0.883080,0.853893,0.842073,0.808343 /)
XEXT_COEFF_550_LKT(87,4)=48.861000 !rg=10.1792 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(87,5,1:6)=(/ 42.539000,42.892000,43.190000,43.582000,44.626000,46.896000 /)
XPIZA_LKT(87,5,1:6)=(/ 0.700808,0.740139,0.762166,0.883143,0.904894,0.944743 /)
XCGA_LKT(87,5,1:6)=(/ 0.896213,0.888487,0.883580,0.854210,0.845567,0.811523 /)
XEXT_COEFF_550_LKT(87,5)=43.182000 !rg=10.1792 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(87,6,1:6)=(/ 37.196000,37.376000,37.724000,38.068000,38.938000,40.380000 /)
XPIZA_LKT(87,6,1:6)=(/ 0.690653,0.730136,0.755936,0.883166,0.901086,0.939043 /)
XCGA_LKT(87,6,1:6)=(/ 0.898593,0.890153,0.884943,0.855550,0.849453,0.817540 /)
XEXT_COEFF_550_LKT(87,6)=37.793000 !rg=10.1792 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(87,7,1:6)=(/ 32.033000,32.168000,32.480000,32.948000,33.474000,34.648000 /)
XPIZA_LKT(87,7,1:6)=(/ 0.678957,0.719169,0.747975,0.881335,0.897036,0.932366 /)
XCGA_LKT(87,7,1:6)=(/ 0.901040,0.892073,0.886743,0.856337,0.852157,0.822463 /)
XEXT_COEFF_550_LKT(87,7)=32.557000 !rg=10.1792 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(87,8,1:6)=(/ 27.387000,27.506000,27.777000,28.037000,28.490000,29.385000 /)
XPIZA_LKT(87,8,1:6)=(/ 0.667152,0.707678,0.739043,0.878106,0.893890,0.924440 /)
XCGA_LKT(87,8,1:6)=(/ 0.903863,0.894673,0.888713,0.857753,0.855150,0.829520 /)
XEXT_COEFF_550_LKT(87,8)=27.750000 !rg=10.1792 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(87,9,1:6)=(/ 23.204000,23.331000,23.479000,23.630000,24.029000,24.864000 /)
XPIZA_LKT(87,9,1:6)=(/ 0.654401,0.695051,0.728033,0.873770,0.891003,0.917134 /)
XCGA_LKT(87,9,1:6)=(/ 0.906903,0.897437,0.890500,0.859183,0.857347,0.833223 /)
XEXT_COEFF_550_LKT(87,9)=23.457000 !rg=10.1792 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(87,10,1:6)=(/ 19.603000,19.710000,19.806000,20.035000,20.304000,20.890000 /)
XPIZA_LKT(87,10,1:6)=(/ 0.642368,0.682264,0.716109,0.867841,0.887244,0.910922 /)
XCGA_LKT(87,10,1:6)=(/ 0.910193,0.900343,0.892780,0.861077,0.859327,0.838843 /)
XEXT_COEFF_550_LKT(87,10)=19.834000 !rg=10.1792 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(87,11,1:6)=(/ 16.465000,16.529000,16.628000,16.712000,16.995000,17.432000 /)
XPIZA_LKT(87,11,1:6)=(/ 0.630819,0.668604,0.703774,0.859592,0.883108,0.905270 /)
XCGA_LKT(87,11,1:6)=(/ 0.913713,0.903420,0.895563,0.862313,0.861600,0.844347 /)
XEXT_COEFF_550_LKT(87,11)=16.620000 !rg=10.1792 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(87,12,1:6)=(/ 13.770000,13.834000,13.890000,14.008000,14.165000,14.493000 /)
XPIZA_LKT(87,12,1:6)=(/ 0.619472,0.655680,0.690106,0.850777,0.878351,0.900702 /)
XCGA_LKT(87,12,1:6)=(/ 0.917337,0.906890,0.898370,0.864770,0.864023,0.848417 /)
XEXT_COEFF_550_LKT(87,12)=13.877000 !rg=10.1792 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(87,13,1:6)=(/ 11.492000,11.532000,11.589000,11.642000,11.789000,12.047000 /)
XPIZA_LKT(87,13,1:6)=(/ 0.609144,0.642546,0.676659,0.839339,0.871930,0.895593 /)
XCGA_LKT(87,13,1:6)=(/ 0.921047,0.910383,0.901620,0.866783,0.866827,0.853240 /)
XEXT_COEFF_550_LKT(87,13)=11.576000 !rg=10.1792 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(87,14,1:6)=(/ 9.584300,9.613800,9.652700,9.714300,9.807800,10.001000 /)
XPIZA_LKT(87,14,1:6)=(/ 0.599758,0.630247,0.663044,0.826871,0.864306,0.891139 /)
XCGA_LKT(87,14,1:6)=(/ 0.924687,0.914083,0.904883,0.869333,0.868493,0.856353 /)
XEXT_COEFF_550_LKT(87,14)=9.661500 !rg=10.1792 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(87,15,1:6)=(/ 7.984900,8.009600,8.034600,8.084100,8.146300,8.303600 /)
XPIZA_LKT(87,15,1:6)=(/ 0.591374,0.619217,0.649836,0.812843,0.854813,0.885005 /)
XCGA_LKT(87,15,1:6)=(/ 0.928243,0.917833,0.908393,0.871653,0.869940,0.858737 /)
XEXT_COEFF_550_LKT(87,15)=8.035200 !rg=10.1792 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(87,16,1:6)=(/ 6.636300,6.653200,6.675300,6.703400,6.765800,6.878200 /)
XPIZA_LKT(87,16,1:6)=(/ 0.583894,0.608670,0.637264,0.797357,0.844255,0.878993 /)
XCGA_LKT(87,16,1:6)=(/ 0.931667,0.921550,0.912127,0.874200,0.872663,0.861717 /)
XEXT_COEFF_550_LKT(87,16)=6.674600 !rg=10.1792 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(87,17,1:6)=(/ 5.522200,5.535300,5.551100,5.574400,5.613300,5.695600 /)
XPIZA_LKT(87,17,1:6)=(/ 0.577379,0.599312,0.625443,0.781144,0.832338,0.872264 /)
XCGA_LKT(87,17,1:6)=(/ 0.934827,0.925213,0.915843,0.877107,0.875227,0.865043 /)
XEXT_COEFF_550_LKT(87,17)=5.551200 !rg=10.1792 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(87,18,1:6)=(/ 0.000000,4.606100,4.617700,4.634800,4.661900,4.723300 /)
XPIZA_LKT(87,18,1:6)=(/ 0.900,0.591044,0.614680,0.764070,0.818400,0.863259 /)
XCGA_LKT(87,18,1:6)=(/ 0.900,0.928773,0.919590,0.879923,0.876943,0.866687 /)
XEXT_COEFF_550_LKT(87,18)=4.616900 !rg=10.1792 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(87,19,1:6)=(/ 0.000000,3.826200,3.835800,3.849300,3.873500,3.919400 /)
XPIZA_LKT(87,19,1:6)=(/ 0.900,0.583592,0.604838,0.746391,0.803515,0.853874 /)
XCGA_LKT(87,19,1:6)=(/ 0.900,0.932140,0.923300,0.883310,0.879487,0.869477 /)
XEXT_COEFF_550_LKT(87,19)=3.834900 !rg=10.1792 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(87,20,1:6)=(/ 0.000000,3.180300,3.187100,3.197300,3.213100,3.247100 /)
XPIZA_LKT(87,20,1:6)=(/ 0.900,0.577131,0.595806,0.728612,0.787548,0.842974 /)
XCGA_LKT(87,20,1:6)=(/ 0.900,0.935297,0.926933,0.886720,0.882423,0.872650 /)
XEXT_COEFF_550_LKT(87,20)=3.187200 !rg=10.1792 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(88,1,1:6)=(/ 55.512000,55.760000,56.132000,57.432000,59.899000,64.820000 /)
XPIZA_LKT(88,1,1:6)=(/ 0.720871,0.757781,0.777807,0.867217,0.920708,0.959773 /)
XCGA_LKT(88,1,1:6)=(/ 0.892477,0.886097,0.880157,0.855663,0.842850,0.816537 /)
XEXT_COEFF_550_LKT(88,1)=56.564000 !rg=11.0277 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(88,2,1:6)=(/ 53.054000,53.267000,53.679000,54.664000,55.461000,58.982000 /)
XPIZA_LKT(88,2,1:6)=(/ 0.717444,0.754211,0.773865,0.873844,0.918098,0.952099 /)
XCGA_LKT(88,2,1:6)=(/ 0.893630,0.885700,0.880580,0.857163,0.840557,0.805727 /)
XEXT_COEFF_550_LKT(88,2)=53.894000 !rg=11.0277 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(88,3,1:6)=(/ 48.972000,49.514000,49.777000,50.311000,51.897000,53.626000 /)
XPIZA_LKT(88,3,1:6)=(/ 0.710452,0.750144,0.771103,0.878967,0.911288,0.952407 /)
XCGA_LKT(88,3,1:6)=(/ 0.894263,0.886987,0.882287,0.855123,0.843040,0.810517 /)
XEXT_COEFF_550_LKT(88,3)=49.609000 !rg=11.0277 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(88,4,1:6)=(/ 44.264000,44.714000,44.935000,45.497000,46.655000,48.268000 /)
XPIZA_LKT(88,4,1:6)=(/ 0.702817,0.743322,0.765695,0.882190,0.905515,0.948084 /)
XCGA_LKT(88,4,1:6)=(/ 0.895723,0.888200,0.883060,0.855063,0.846423,0.814323 /)
XEXT_COEFF_550_LKT(88,4)=44.839000 !rg=11.0277 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(88,5,1:6)=(/ 39.268000,39.485000,39.905000,40.381000,41.179000,42.792000 /)
XPIZA_LKT(88,5,1:6)=(/ 0.693744,0.734140,0.759346,0.883275,0.899720,0.942125 /)
XCGA_LKT(88,5,1:6)=(/ 0.897890,0.889607,0.884670,0.855080,0.848613,0.816303 /)
XEXT_COEFF_550_LKT(88,5)=39.869000 !rg=11.0277 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(88,6,1:6)=(/ 34.273000,34.552000,34.762000,35.154000,35.686000,36.931000 /)
XPIZA_LKT(88,6,1:6)=(/ 0.683179,0.724744,0.752188,0.882176,0.896942,0.935577 /)
XCGA_LKT(88,6,1:6)=(/ 0.900023,0.891810,0.885863,0.856087,0.852300,0.822393 /)
XEXT_COEFF_550_LKT(88,6)=34.794000 !rg=11.0277 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(88,7,1:6)=(/ 29.518000,29.714000,29.865000,30.268000,30.787000,31.678000 /)
XPIZA_LKT(88,7,1:6)=(/ 0.671639,0.713291,0.743173,0.880340,0.894963,0.928961 /)
XCGA_LKT(88,7,1:6)=(/ 0.902570,0.893870,0.887073,0.857470,0.855160,0.828790 /)
XEXT_COEFF_550_LKT(88,7)=29.868000 !rg=11.0277 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(88,8,1:6)=(/ 25.237000,25.389000,25.519000,25.817000,26.196000,26.985000 /)
XPIZA_LKT(88,8,1:6)=(/ 0.659694,0.701481,0.733643,0.876537,0.891546,0.920618 /)
XCGA_LKT(88,8,1:6)=(/ 0.905527,0.896143,0.889343,0.858707,0.856493,0.832300 /)
XEXT_COEFF_550_LKT(88,8)=25.499000 !rg=11.0277 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(88,9,1:6)=(/ 21.407000,21.502000,21.648000,21.846000,22.189000,22.847000 /)
XPIZA_LKT(88,9,1:6)=(/ 0.647857,0.688207,0.722461,0.871182,0.888391,0.913423 /)
XCGA_LKT(88,9,1:6)=(/ 0.908617,0.898867,0.891693,0.859977,0.858830,0.837037 /)
XEXT_COEFF_550_LKT(88,9)=21.672000 !rg=11.0277 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(88,10,1:6)=(/ 18.087000,18.160000,18.291000,18.426000,18.639000,19.155000 /)
XPIZA_LKT(88,10,1:6)=(/ 0.636227,0.675225,0.710705,0.864157,0.885514,0.908142 /)
XCGA_LKT(88,10,1:6)=(/ 0.911957,0.901750,0.894213,0.861883,0.861553,0.842667 /)
XEXT_COEFF_550_LKT(88,10)=18.261000 !rg=11.0277 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(88,11,1:6)=(/ 15.184000,15.242000,15.319000,15.450000,15.649000,16.112000 /)
XPIZA_LKT(88,11,1:6)=(/ 0.624696,0.662136,0.697200,0.855985,0.880729,0.902346 /)
XCGA_LKT(88,11,1:6)=(/ 0.915560,0.905010,0.896757,0.864003,0.862733,0.845810 /)
XEXT_COEFF_550_LKT(88,11)=15.326000 !rg=11.0277 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(88,12,1:6)=(/ 12.705000,12.752000,12.822000,12.902000,13.048000,13.328000 /)
XPIZA_LKT(88,12,1:6)=(/ 0.614059,0.649068,0.684122,0.845597,0.875715,0.898266 /)
XCGA_LKT(88,12,1:6)=(/ 0.919177,0.908530,0.899953,0.865513,0.866113,0.851677 /)
XEXT_COEFF_550_LKT(88,12)=12.809000 !rg=11.0277 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(88,13,1:6)=(/ 10.602000,10.635000,10.682000,10.750000,10.870000,11.081000 /)
XPIZA_LKT(88,13,1:6)=(/ 0.604248,0.636412,0.670210,0.834136,0.868783,0.893656 /)
XCGA_LKT(88,13,1:6)=(/ 0.922827,0.912093,0.903053,0.868387,0.867567,0.854707 /)
XEXT_COEFF_550_LKT(88,13)=10.685000 !rg=11.0277 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(88,14,1:6)=(/ 8.842900,8.871900,8.902400,8.947700,9.037900,9.199800 /)
XPIZA_LKT(88,14,1:6)=(/ 0.595350,0.624839,0.656568,0.820689,0.860367,0.888327 /)
XCGA_LKT(88,14,1:6)=(/ 0.926450,0.915877,0.906547,0.870243,0.869580,0.857923 /)
XEXT_COEFF_550_LKT(88,14)=8.901700 !rg=11.0277 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(88,15,1:6)=(/ 7.367900,7.388700,7.414800,7.446700,7.519600,7.635400 /)
XPIZA_LKT(88,15,1:6)=(/ 0.587507,0.613989,0.644001,0.806114,0.850950,0.882826 /)
XCGA_LKT(88,15,1:6)=(/ 0.929907,0.919593,0.910173,0.872877,0.872117,0.861100 /)
XEXT_COEFF_550_LKT(88,15)=7.413000 !rg=11.0277 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(88,16,1:6)=(/ 6.122700,6.138900,6.156700,6.188300,6.231100,6.334600 /)
XPIZA_LKT(88,16,1:6)=(/ 0.580456,0.603986,0.631530,0.790383,0.839185,0.876047 /)
XCGA_LKT(88,16,1:6)=(/ 0.933227,0.923333,0.913827,0.875597,0.873893,0.863257 /)
XEXT_COEFF_550_LKT(88,16)=6.156000 !rg=11.0277 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(88,17,1:6)=(/ 5.095600,5.107700,5.121300,5.142500,5.176600,5.247500 /)
XPIZA_LKT(88,17,1:6)=(/ 0.574385,0.595141,0.620155,0.773560,0.826330,0.868093 /)
XCGA_LKT(88,17,1:6)=(/ 0.936320,0.926930,0.917620,0.878367,0.876190,0.865673 /)
XEXT_COEFF_550_LKT(88,17)=5.120500 !rg=11.0277 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(88,18,1:6)=(/ 0.000000,4.249700,4.261000,4.275600,4.304900,4.355200 /)
XPIZA_LKT(88,18,1:6)=(/ 0.900,0.587285,0.609903,0.756180,0.812313,0.859799 /)
XCGA_LKT(88,18,1:6)=(/ 0.900,0.930373,0.921370,0.881480,0.878717,0.868900 /)
XEXT_COEFF_550_LKT(88,18)=4.260100 !rg=11.0277 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(88,19,1:6)=(/ 0.000000,3.530700,3.538500,3.552300,3.569100,3.611300 /)
XPIZA_LKT(88,19,1:6)=(/ 0.900,0.580343,0.600342,0.738588,0.796466,0.849181 /)
XCGA_LKT(88,19,1:6)=(/ 0.900,0.933657,0.925037,0.884827,0.880917,0.870790 /)
XEXT_COEFF_550_LKT(88,19)=3.537700 !rg=11.0277 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(88,20,1:6)=(/ 0.000000,0.000000,2.940600,2.949000,2.963200,2.992200 /)
XPIZA_LKT(88,20,1:6)=(/ 0.900,0.900,0.591926,0.720609,0.780093,0.837368 /)
XCGA_LKT(88,20,1:6)=(/ 0.900,0.900,0.928547,0.888270,0.883763,0.873367 /)
XEXT_COEFF_550_LKT(88,20)=2.940200 !rg=11.0277 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(89,1,1:6)=(/ 51.195000,50.761000,51.720000,52.883000,53.924000,54.191000 /)
XPIZA_LKT(89,1,1:6)=(/ 0.714388,0.750069,0.775068,0.878802,0.911335,0.941271 /)
XCGA_LKT(89,1,1:6)=(/ 0.895020,0.884810,0.880437,0.854543,0.836173,0.785943 /)
XEXT_COEFF_550_LKT(89,1)=51.749000 !rg=11.9469 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(89,2,1:6)=(/ 48.757000,49.291000,49.549000,50.466000,50.897000,53.463000 /)
XPIZA_LKT(89,2,1:6)=(/ 0.709401,0.749712,0.772543,0.879358,0.911656,0.950959 /)
XCGA_LKT(89,2,1:6)=(/ 0.894563,0.887233,0.881720,0.856143,0.846417,0.807563 /)
XEXT_COEFF_550_LKT(89,2)=49.699000 !rg=11.9469 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(89,3,1:6)=(/ 45.213000,45.427000,45.910000,46.516000,47.603000,49.190000 /)
XPIZA_LKT(89,3,1:6)=(/ 0.703828,0.743405,0.768504,0.881553,0.905422,0.948003 /)
XCGA_LKT(89,3,1:6)=(/ 0.895903,0.887460,0.882787,0.855533,0.846740,0.814897 /)
XEXT_COEFF_550_LKT(89,3)=45.972000 !rg=11.9469 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(89,4,1:6)=(/ 40.880000,41.066000,41.465000,42.060000,42.872000,44.302000 /)
XPIZA_LKT(89,4,1:6)=(/ 0.695723,0.736589,0.763119,0.882906,0.899815,0.943662 /)
XCGA_LKT(89,4,1:6)=(/ 0.897497,0.888777,0.883880,0.855350,0.849467,0.819043 /)
XEXT_COEFF_550_LKT(89,4)=41.548000 !rg=11.9469 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(89,5,1:6)=(/ 36.177000,36.503000,36.671000,37.116000,37.973000,39.192000 /)
XPIZA_LKT(89,5,1:6)=(/ 0.686003,0.728518,0.756229,0.882652,0.896064,0.937998 /)
XCGA_LKT(89,5,1:6)=(/ 0.899373,0.891143,0.885233,0.856393,0.851810,0.821830 /)
XEXT_COEFF_550_LKT(89,5)=36.591000 !rg=11.9469 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(89,6,1:6)=(/ 31.606000,31.764000,32.095000,32.349000,33.081000,34.198000 /)
XPIZA_LKT(89,6,1:6)=(/ 0.675789,0.717688,0.748306,0.880847,0.894243,0.931094 /)
XCGA_LKT(89,6,1:6)=(/ 0.901687,0.892697,0.887097,0.856883,0.855577,0.826447 /)
XEXT_COEFF_550_LKT(89,6)=32.063000 !rg=11.9469 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(89,7,1:6)=(/ 27.246000,27.367000,27.566000,27.837000,28.338000,29.241000 /)
XPIZA_LKT(89,7,1:6)=(/ 0.664756,0.706280,0.738836,0.878296,0.890961,0.923322 /)
XCGA_LKT(89,7,1:6)=(/ 0.904300,0.894967,0.888580,0.857593,0.855943,0.831993 /)
XEXT_COEFF_550_LKT(89,7)=27.572000 !rg=11.9469 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(89,8,1:6)=(/ 23.286000,23.405000,23.556000,23.752000,24.144000,24.882000 /)
XPIZA_LKT(89,8,1:6)=(/ 0.652839,0.694408,0.728637,0.874047,0.889350,0.916647 /)
XCGA_LKT(89,8,1:6)=(/ 0.907310,0.897537,0.890580,0.859330,0.859313,0.836880 /)
XEXT_COEFF_550_LKT(89,8)=23.537000 !rg=11.9469 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(89,9,1:6)=(/ 19.746000,19.855000,19.944000,20.154000,20.461000,21.023000 /)
XPIZA_LKT(89,9,1:6)=(/ 0.641153,0.681697,0.716516,0.868032,0.886254,0.910103 /)
XCGA_LKT(89,9,1:6)=(/ 0.910493,0.900537,0.892890,0.861293,0.859757,0.840993 /)
XEXT_COEFF_550_LKT(89,9)=19.921000 !rg=11.9469 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(89,10,1:6)=(/ 16.683000,16.753000,16.834000,16.970000,17.212000,17.626000 /)
XPIZA_LKT(89,10,1:6)=(/ 0.629931,0.668596,0.704015,0.860692,0.883071,0.905197 /)
XCGA_LKT(89,10,1:6)=(/ 0.913823,0.903393,0.895280,0.862783,0.862327,0.845390 /)
XEXT_COEFF_550_LKT(89,10)=16.856000 !rg=11.9469 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(89,11,1:6)=(/ 14.007000,14.064000,14.123000,14.220000,14.399000,14.750000 /)
XPIZA_LKT(89,11,1:6)=(/ 0.619083,0.655505,0.690933,0.851455,0.878309,0.899818 /)
XCGA_LKT(89,11,1:6)=(/ 0.917280,0.906757,0.898143,0.864803,0.863983,0.848580 /)
XEXT_COEFF_550_LKT(89,11)=14.114000 !rg=11.9469 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(89,12,1:6)=(/ 11.722000,11.763000,11.812000,11.899000,12.027000,12.280000 /)
XPIZA_LKT(89,12,1:6)=(/ 0.608948,0.642710,0.677398,0.840808,0.872655,0.895453 /)
XCGA_LKT(89,12,1:6)=(/ 0.920970,0.910270,0.901320,0.866940,0.866813,0.853110 /)
XEXT_COEFF_550_LKT(89,12)=11.824000 !rg=11.9469 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(89,13,1:6)=(/ 9.779100,9.814400,9.848800,9.916600,10.006000,10.209000 /)
XPIZA_LKT(89,13,1:6)=(/ 0.599450,0.630648,0.663642,0.828367,0.864937,0.890883 /)
XCGA_LKT(89,13,1:6)=(/ 0.924617,0.913950,0.904613,0.869457,0.868540,0.856700 /)
XEXT_COEFF_550_LKT(89,13)=9.849000 !rg=11.9469 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(89,14,1:6)=(/ 8.159200,8.182800,8.215900,8.254800,8.345300,8.481300 /)
XPIZA_LKT(89,14,1:6)=(/ 0.591221,0.619369,0.650675,0.814353,0.856206,0.885380 /)
XCGA_LKT(89,14,1:6)=(/ 0.928157,0.917613,0.908273,0.871480,0.870537,0.859563 /)
XEXT_COEFF_550_LKT(89,14)=8.215500 !rg=11.9469 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(89,15,1:6)=(/ 6.797400,6.817900,6.837900,6.878200,6.920100,7.042200 /)
XPIZA_LKT(89,15,1:6)=(/ 0.583797,0.609092,0.637962,0.799517,0.845666,0.879652 /)
XCGA_LKT(89,15,1:6)=(/ 0.931527,0.921367,0.911863,0.874123,0.872550,0.862480 /)
XEXT_COEFF_550_LKT(89,15)=6.837400 !rg=11.9469 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(89,16,1:6)=(/ 5.649400,5.663500,5.680700,5.701200,5.745700,5.829200 /)
XPIZA_LKT(89,16,1:6)=(/ 0.577196,0.599458,0.626048,0.782778,0.833599,0.872462 /)
XCGA_LKT(89,16,1:6)=(/ 0.934777,0.925050,0.915640,0.876720,0.874743,0.864303 /)
XEXT_COEFF_550_LKT(89,16)=5.678200 !rg=11.9469 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(89,17,1:6)=(/ 0.000000,4.712200,4.725700,4.740300,4.778500,4.838300 /)
XPIZA_LKT(89,17,1:6)=(/ 0.900,0.591093,0.615170,0.765783,0.820325,0.864739 /)
XCGA_LKT(89,17,1:6)=(/ 0.900,0.928600,0.919380,0.879717,0.877390,0.867603 /)
XEXT_COEFF_550_LKT(89,17)=4.725100 !rg=11.9469 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(89,18,1:6)=(/ 0.000000,3.921600,3.930500,3.946800,3.963800,4.011300 /)
XPIZA_LKT(89,18,1:6)=(/ 0.900,0.583804,0.605208,0.748514,0.805223,0.854828 /)
XCGA_LKT(89,18,1:6)=(/ 0.900,0.931957,0.923050,0.882953,0.879303,0.869210 /)
XEXT_COEFF_550_LKT(89,18)=3.930700 !rg=11.9469 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(89,19,1:6)=(/ 0.000000,3.257900,3.264900,3.273800,3.292400,3.328200 /)
XPIZA_LKT(89,19,1:6)=(/ 0.900,0.577284,0.596232,0.730297,0.789318,0.844226 /)
XCGA_LKT(89,19,1:6)=(/ 0.900,0.935123,0.926683,0.886260,0.881903,0.871950 /)
XEXT_COEFF_550_LKT(89,19)=3.264100 !rg=11.9469 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(89,20,1:6)=(/ 0.000000,0.000000,2.713400,2.720100,2.736200,2.760600 /)
XPIZA_LKT(89,20,1:6)=(/ 0.900,0.900,0.588246,0.712654,0.772831,0.832031 /)
XCGA_LKT(89,20,1:6)=(/ 0.900,0.900,0.930157,0.889963,0.885193,0.874817 /)
XEXT_COEFF_550_LKT(89,20)=2.713200 !rg=11.9469 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(90,1,1:6)=(/ 46.898000,47.745000,47.831000,48.518000,49.875000,48.610000 /)
XPIZA_LKT(90,1,1:6)=(/ 0.704683,0.747466,0.771654,0.884624,0.907904,0.947709 /)
XCGA_LKT(90,1,1:6)=(/ 0.896217,0.888050,0.883020,0.854570,0.845580,0.798470 /)
XEXT_COEFF_550_LKT(90,1)=47.929000 !rg=12.9427 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(90,2,1:6)=(/ 44.980000,45.181000,45.666000,46.293000,47.482000,49.120000 /)
XPIZA_LKT(90,2,1:6)=(/ 0.701931,0.742768,0.768571,0.883565,0.905158,0.947852 /)
XCGA_LKT(90,2,1:6)=(/ 0.896203,0.887640,0.882937,0.855350,0.848657,0.810993 /)
XEXT_COEFF_550_LKT(90,2)=45.605000 !rg=12.9427 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(90,3,1:6)=(/ 41.652000,42.125000,42.308000,42.758000,43.877000,45.588000 /)
XPIZA_LKT(90,3,1:6)=(/ 0.695711,0.738909,0.765201,0.883122,0.898171,0.942340 /)
XCGA_LKT(90,3,1:6)=(/ 0.897317,0.889427,0.883080,0.855873,0.849877,0.814103 /)
XEXT_COEFF_550_LKT(90,3)=42.434000 !rg=12.9427 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(90,4,1:6)=(/ 37.672000,37.931000,38.192000,38.622000,39.603000,41.125000 /)
XPIZA_LKT(90,4,1:6)=(/ 0.687812,0.730591,0.759405,0.883329,0.894298,0.937504 /)
XCGA_LKT(90,4,1:6)=(/ 0.898963,0.890573,0.884130,0.856407,0.852223,0.817603 /)
XEXT_COEFF_550_LKT(90,4)=38.306000 !rg=12.9427 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(90,5,1:6)=(/ 33.387000,33.533000,33.881000,34.219000,34.790000,35.977000 /)
XPIZA_LKT(90,5,1:6)=(/ 0.678992,0.721197,0.752383,0.882513,0.893056,0.933108 /)
XCGA_LKT(90,5,1:6)=(/ 0.901030,0.891830,0.886250,0.856640,0.854803,0.826703 /)
XEXT_COEFF_550_LKT(90,5)=33.900000 !rg=12.9427 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(90,6,1:6)=(/ 29.144000,29.285000,29.499000,29.912000,30.389000,31.551000 /)
XPIZA_LKT(90,6,1:6)=(/ 0.668461,0.710883,0.743206,0.880414,0.891661,0.926624 /)
XCGA_LKT(90,6,1:6)=(/ 0.903317,0.894240,0.887603,0.857990,0.856503,0.828157 /)
XEXT_COEFF_550_LKT(90,6)=29.533000 !rg=12.9427 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(90,7,1:6)=(/ 25.115000,25.243000,25.396000,25.665000,26.161000,27.094000 /)
XPIZA_LKT(90,7,1:6)=(/ 0.657330,0.699878,0.733251,0.876697,0.889298,0.917962 /)
XCGA_LKT(90,7,1:6)=(/ 0.906040,0.896427,0.889377,0.859360,0.857987,0.832257 /)
XEXT_COEFF_550_LKT(90,7)=25.414000 !rg=12.9427 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(90,8,1:6)=(/ 21.479000,21.582000,21.713000,21.902000,22.215000,22.891000 /)
XPIZA_LKT(90,8,1:6)=(/ 0.646007,0.687492,0.722639,0.871426,0.887814,0.913223 /)
XCGA_LKT(90,8,1:6)=(/ 0.909103,0.898990,0.891620,0.860813,0.860297,0.839453 /)
XEXT_COEFF_550_LKT(90,8)=21.752000 !rg=12.9427 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(90,9,1:6)=(/ 18.220000,18.287000,18.418000,18.525000,18.787000,19.280000 /)
XPIZA_LKT(90,9,1:6)=(/ 0.634884,0.674504,0.710811,0.864663,0.885110,0.906854 /)
XCGA_LKT(90,9,1:6)=(/ 0.912310,0.901907,0.894297,0.861693,0.862060,0.844740 /)
XEXT_COEFF_550_LKT(90,9)=18.423000 !rg=12.9427 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(90,10,1:6)=(/ 15.389000,15.461000,15.530000,15.624000,15.820000,16.245000 /)
XPIZA_LKT(90,10,1:6)=(/ 0.623943,0.661994,0.697831,0.856533,0.880868,0.901716 /)
XCGA_LKT(90,10,1:6)=(/ 0.915637,0.905227,0.896753,0.863790,0.863633,0.847883 /)
XEXT_COEFF_550_LKT(90,10)=15.528000 !rg=12.9427 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(90,11,1:6)=(/ 12.924000,12.971000,13.033000,13.131000,13.290000,13.570000 /)
XPIZA_LKT(90,11,1:6)=(/ 0.613595,0.648981,0.684422,0.846923,0.875963,0.897481 /)
XCGA_LKT(90,11,1:6)=(/ 0.919220,0.908447,0.899723,0.865887,0.865847,0.851923 /)
XEXT_COEFF_550_LKT(90,11)=13.048000 !rg=12.9427 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(90,12,1:6)=(/ 10.815000,10.853000,10.894000,10.944000,11.064000,11.294000 /)
XPIZA_LKT(90,12,1:6)=(/ 0.603893,0.636731,0.670712,0.834985,0.869193,0.893090 /)
XCGA_LKT(90,12,1:6)=(/ 0.922817,0.912050,0.902923,0.867653,0.867557,0.855430 /)
XEXT_COEFF_550_LKT(90,12)=10.891000 !rg=12.9427 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(90,13,1:6)=(/ 9.024500,9.052000,9.089800,9.125400,9.233900,9.402700 /)
XPIZA_LKT(90,13,1:6)=(/ 0.595088,0.624871,0.657492,0.822082,0.861477,0.888035 /)
XCGA_LKT(90,13,1:6)=(/ 0.926393,0.915723,0.906380,0.870150,0.869950,0.858307 /)
XEXT_COEFF_550_LKT(90,13)=9.091900 !rg=12.9427 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(90,14,1:6)=(/ 7.527800,7.550600,7.573600,7.616900,7.666700,7.796900 /)
XPIZA_LKT(90,14,1:6)=(/ 0.587219,0.614171,0.644481,0.807922,0.851539,0.882573 /)
XCGA_LKT(90,14,1:6)=(/ 0.929867,0.919427,0.909863,0.872860,0.871223,0.860550 /)
XEXT_COEFF_550_LKT(90,14)=7.570900 !rg=12.9427 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(90,15,1:6)=(/ 6.272500,6.288100,6.310800,6.333300,6.383300,6.475000 /)
XPIZA_LKT(90,15,1:6)=(/ 0.580383,0.604176,0.632386,0.792070,0.840784,0.876859 /)
XCGA_LKT(90,15,1:6)=(/ 0.933120,0.923120,0.913637,0.875237,0.874030,0.863870 /)
XEXT_COEFF_550_LKT(90,15)=6.305100 !rg=12.9427 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(90,16,1:6)=(/ 5.213300,5.224500,5.240400,5.262800,5.298300,5.369300 /)
XPIZA_LKT(90,16,1:6)=(/ 0.574247,0.595168,0.620636,0.775356,0.828166,0.869511 /)
XCGA_LKT(90,16,1:6)=(/ 0.936263,0.926733,0.917403,0.878237,0.876570,0.866383 /)
XEXT_COEFF_550_LKT(90,16)=5.240900 !rg=12.9427 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(90,17,1:6)=(/ 0.000000,4.348600,4.358800,4.376500,4.397700,4.441700 /)
XPIZA_LKT(90,17,1:6)=(/ 0.900,0.587396,0.610180,0.758078,0.813665,0.860258 /)
XCGA_LKT(90,17,1:6)=(/ 0.900,0.930233,0.921120,0.881283,0.878227,0.868173 /)
XEXT_COEFF_550_LKT(90,17)=4.358400 !rg=12.9427 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(90,18,1:6)=(/ 0.000000,3.618000,3.626900,3.638400,3.657300,3.696400 /)
XPIZA_LKT(90,18,1:6)=(/ 0.900,0.580439,0.600847,0.740262,0.798539,0.850577 /)
XCGA_LKT(90,18,1:6)=(/ 0.900,0.933477,0.924773,0.884447,0.880830,0.871013 /)
XEXT_COEFF_550_LKT(90,18)=3.625400 !rg=12.9427 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(90,19,1:6)=(/ 0.000000,3.006000,3.012400,3.020500,3.036400,3.067800 /)
XPIZA_LKT(90,19,1:6)=(/ 0.900,0.574304,0.592284,0.722400,0.782244,0.839253 /)
XCGA_LKT(90,19,1:6)=(/ 0.900,0.936600,0.928307,0.887863,0.883580,0.873883 /)
XEXT_COEFF_550_LKT(90,19)=3.012100 !rg=12.9427 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(90,20,1:6)=(/ 0.000000,0.000000,2.503500,2.510500,2.520000,2.537000 /)
XPIZA_LKT(90,20,1:6)=(/ 0.900,0.900,0.584701,0.704877,0.764919,0.825679 /)
XCGA_LKT(90,20,1:6)=(/ 0.900,0.900,0.931747,0.891703,0.886403,0.875340 /)
XEXT_COEFF_550_LKT(90,20)=2.503200 !rg=12.9427 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(91,1,1:6)=(/ 43.083000,43.620000,43.935000,44.660000,45.419000,48.941000 /)
XPIZA_LKT(91,1,1:6)=(/ 0.696820,0.740044,0.767753,0.889189,0.902857,0.947601 /)
XCGA_LKT(91,1,1:6)=(/ 0.896167,0.888927,0.882840,0.854867,0.849970,0.825273 /)
XEXT_COEFF_550_LKT(91,1)=44.046000 !rg=14.0216 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(91,2,1:6)=(/ 41.493000,41.861000,42.159000,42.613000,43.420000,45.702000 /)
XPIZA_LKT(91,2,1:6)=(/ 0.694209,0.738008,0.765279,0.886242,0.899016,0.942680 /)
XCGA_LKT(91,2,1:6)=(/ 0.897797,0.889540,0.883540,0.855633,0.853320,0.814310 /)
XEXT_COEFF_550_LKT(91,2)=42.027000 !rg=14.0216 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(91,3,1:6)=(/ 38.433000,38.720000,39.121000,39.327000,39.925000,41.990000 /)
XPIZA_LKT(91,3,1:6)=(/ 0.687990,0.731588,0.761480,0.885151,0.892995,0.939669 /)
XCGA_LKT(91,3,1:6)=(/ 0.898990,0.890497,0.884413,0.856250,0.854630,0.820143 /)
XEXT_COEFF_550_LKT(91,3)=38.965000 !rg=14.0216 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(91,4,1:6)=(/ 34.734000,34.968000,35.191000,35.489000,36.081000,37.706000 /)
XPIZA_LKT(91,4,1:6)=(/ 0.680323,0.724102,0.754582,0.884094,0.890482,0.934701 /)
XCGA_LKT(91,4,1:6)=(/ 0.900580,0.891760,0.885467,0.856913,0.855437,0.824150 /)
XEXT_COEFF_550_LKT(91,4)=35.178000 !rg=14.0216 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(91,5,1:6)=(/ 30.768000,31.025000,31.140000,31.495000,32.169000,33.289000 /)
XPIZA_LKT(91,5,1:6)=(/ 0.671002,0.715263,0.747063,0.881644,0.889013,0.927481 /)
XCGA_LKT(91,5,1:6)=(/ 0.902617,0.893667,0.886590,0.857690,0.856737,0.826893 /)
XEXT_COEFF_550_LKT(91,5)=31.235000 !rg=14.0216 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(91,6,1:6)=(/ 26.879000,27.032000,27.184000,27.527000,27.996000,29.002000 /)
XPIZA_LKT(91,6,1:6)=(/ 0.661072,0.704495,0.737926,0.879087,0.889260,0.921227 /)
XCGA_LKT(91,6,1:6)=(/ 0.905040,0.895447,0.888873,0.859013,0.858023,0.831600 /)
XEXT_COEFF_550_LKT(91,6)=27.178000 !rg=14.0216 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(91,7,1:6)=(/ 23.162000,23.283000,23.413000,23.596000,23.980000,24.725000 /)
XPIZA_LKT(91,7,1:6)=(/ 0.650448,0.692656,0.728095,0.874266,0.887609,0.914681 /)
XCGA_LKT(91,7,1:6)=(/ 0.907640,0.897943,0.890627,0.859950,0.859093,0.837333 /)
XEXT_COEFF_550_LKT(91,7)=23.375000 !rg=14.0216 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(91,8,1:6)=(/ 19.806000,19.908000,20.012000,20.157000,20.472000,21.087000 /)
XPIZA_LKT(91,8,1:6)=(/ 0.639294,0.680671,0.716499,0.868435,0.886010,0.909178 /)
XCGA_LKT(91,8,1:6)=(/ 0.910850,0.900677,0.892800,0.861493,0.861963,0.842680 /)
XEXT_COEFF_550_LKT(91,8)=19.995000 !rg=14.0216 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(91,9,1:6)=(/ 16.802000,16.893000,16.948000,17.122000,17.345000,17.811000 /)
XPIZA_LKT(91,9,1:6)=(/ 0.628501,0.668216,0.703908,0.861339,0.882632,0.902823 /)
XCGA_LKT(91,9,1:6)=(/ 0.914110,0.903687,0.895293,0.863440,0.863623,0.845997 /)
XEXT_COEFF_550_LKT(91,9)=16.963000 !rg=14.0216 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(91,10,1:6)=(/ 14.201000,14.253000,14.333000,14.418000,14.598000,14.904000 /)
XPIZA_LKT(91,10,1:6)=(/ 0.618277,0.655151,0.691423,0.852257,0.878713,0.899172 /)
XCGA_LKT(91,10,1:6)=(/ 0.917507,0.906840,0.898420,0.864907,0.865453,0.850977 /)
XEXT_COEFF_550_LKT(91,10)=14.315000 !rg=14.0216 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(91,11,1:6)=(/ 11.924000,11.971000,12.015000,12.099000,12.239000,12.484000 /)
XPIZA_LKT(91,11,1:6)=(/ 0.608306,0.642830,0.677755,0.841837,0.873084,0.895392 /)
XCGA_LKT(91,11,1:6)=(/ 0.921067,0.910280,0.901160,0.866940,0.866717,0.854563 /)
XEXT_COEFF_550_LKT(91,11)=12.011000 !rg=14.0216 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(91,12,1:6)=(/ 9.979100,10.009000,10.054000,10.100000,10.226000,10.398000 /)
XPIZA_LKT(91,12,1:6)=(/ 0.599194,0.630622,0.664493,0.829403,0.865667,0.890389 /)
XCGA_LKT(91,12,1:6)=(/ 0.924647,0.913873,0.904617,0.869007,0.868377,0.857587 /)
XEXT_COEFF_550_LKT(91,12)=10.053000 !rg=14.0216 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(91,13,1:6)=(/ 8.325800,8.352100,8.379200,8.437600,8.478500,8.625400 /)
XPIZA_LKT(91,13,1:6)=(/ 0.590822,0.619366,0.651043,0.816112,0.856613,0.885274 /)
XCGA_LKT(91,13,1:6)=(/ 0.928137,0.917560,0.908007,0.871660,0.869917,0.859383 /)
XEXT_COEFF_550_LKT(91,13)=8.378800 !rg=14.0216 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(91,14,1:6)=(/ 6.946100,6.964400,6.988900,7.014300,7.077800,7.189200 /)
XPIZA_LKT(91,14,1:6)=(/ 0.583572,0.609033,0.638609,0.800742,0.846995,0.879803 /)
XCGA_LKT(91,14,1:6)=(/ 0.931487,0.921200,0.911647,0.873870,0.872933,0.862873 /)
XEXT_COEFF_550_LKT(91,14)=6.983600 !rg=14.0216 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(91,15,1:6)=(/ 5.787600,5.801100,5.818200,5.844600,5.884200,5.967000 /)
XPIZA_LKT(91,15,1:6)=(/ 0.577125,0.599664,0.626522,0.784825,0.835199,0.873376 /)
XCGA_LKT(91,15,1:6)=(/ 0.934677,0.924820,0.915377,0.876607,0.874987,0.864920 /)
XEXT_COEFF_550_LKT(91,15)=5.817700 !rg=14.0216 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(91,16,1:6)=(/ 0.000000,4.821400,4.832600,4.853700,4.880100,4.939600 /)
XPIZA_LKT(91,16,1:6)=(/ 0.900,0.591163,0.615419,0.767632,0.821745,0.865410 /)
XCGA_LKT(91,16,1:6)=(/ 0.900,0.928473,0.919087,0.879557,0.877047,0.867467 /)
XEXT_COEFF_550_LKT(91,16)=4.832700 !rg=14.0216 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(91,17,1:6)=(/ 0.000000,4.011700,4.022600,4.034600,4.059900,4.108300 /)
XPIZA_LKT(91,17,1:6)=(/ 0.900,0.583739,0.605614,0.749975,0.807126,0.856391 /)
XCGA_LKT(91,17,1:6)=(/ 0.900,0.931837,0.922877,0.882590,0.879347,0.870313 /)
XEXT_COEFF_550_LKT(91,17)=4.020400 !rg=14.0216 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(91,18,1:6)=(/ 0.000000,3.338300,3.345500,3.356700,3.373800,3.407000 /)
XPIZA_LKT(91,18,1:6)=(/ 0.900,0.577271,0.596519,0.732272,0.791320,0.845607 /)
XCGA_LKT(91,18,1:6)=(/ 0.900,0.935000,0.926430,0.886047,0.881943,0.872207 /)
XEXT_COEFF_550_LKT(91,18)=3.345500 !rg=14.0216 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(91,19,1:6)=(/ 0.000000,0.000000,2.779100,2.787100,2.798800,2.821400 /)
XPIZA_LKT(91,19,1:6)=(/ 0.900,0.900,0.588411,0.714432,0.774547,0.833311 /)
XCGA_LKT(91,19,1:6)=(/ 0.900,0.900,0.929967,0.889603,0.884670,0.874490 /)
XEXT_COEFF_550_LKT(91,19)=2.778600 !rg=14.0216 sigma=2.85 wvl=0.55
 
IF (LHOOK) CALL DR_HOOK('MODE_DUSTOPT:DUST_OPT_LKT_SET9',1,ZHOOK_HANDLE)
END SUBROUTINE DUST_OPT_LKT_SET9

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE DUST_OPT_LKT_SET10()

  USE MODD_DUST_OPT_LKT
  
  IMPLICIT NONE 

REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_DUSTOPT:DUST_OPT_LKT_SET10',0,ZHOOK_HANDLE)
XEXT_COEFF_WVL_LKT(91,20,1:6)=(/ 0.000000,0.000000,2.310300,2.315500,2.325900,2.345100 /)
XPIZA_LKT(91,20,1:6)=(/ 0.900,0.900,0.581451,0.696978,0.757133,0.819978 /)
XCGA_LKT(91,20,1:6)=(/ 0.900,0.900,0.933277,0.893380,0.887830,0.877120 /)
XEXT_COEFF_550_LKT(91,20)=2.309500 !rg=14.0216 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(92,1,1:6)=(/ 39.896000,40.228000,40.329000,40.793000,41.609000,43.289000 /)
XPIZA_LKT(92,1,1:6)=(/ 0.689281,0.734307,0.762122,0.889002,0.897673,0.944111 /)
XCGA_LKT(92,1,1:6)=(/ 0.898053,0.890873,0.883123,0.854543,0.858030,0.829960 /)
XEXT_COEFF_550_LKT(92,1)=40.360000 !rg=15.1904 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(92,2,1:6)=(/ 38.253000,38.533000,38.877000,39.304000,40.085000,41.706000 /)
XPIZA_LKT(92,2,1:6)=(/ 0.686222,0.730793,0.762288,0.888554,0.889335,0.939827 /)
XCGA_LKT(92,2,1:6)=(/ 0.899353,0.890703,0.884700,0.855893,0.854430,0.822037 /)
XEXT_COEFF_550_LKT(92,2)=38.823000 !rg=15.1904 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(92,3,1:6)=(/ 35.469000,35.627000,35.966000,36.295000,36.950000,38.360000 /)
XPIZA_LKT(92,3,1:6)=(/ 0.680671,0.724477,0.756925,0.886029,0.886681,0.936280 /)
XCGA_LKT(92,3,1:6)=(/ 0.900600,0.891540,0.885267,0.856120,0.856197,0.824533 /)
XEXT_COEFF_550_LKT(92,3)=36.026000 !rg=15.1904 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(92,4,1:6)=(/ 32.046000,32.215000,32.492000,32.862000,33.385000,34.530000 /)
XPIZA_LKT(92,4,1:6)=(/ 0.672747,0.716830,0.750229,0.883909,0.886070,0.930499 /)
XCGA_LKT(92,4,1:6)=(/ 0.902243,0.893073,0.886583,0.857123,0.857403,0.827827 /)
XEXT_COEFF_550_LKT(92,4)=32.513000 !rg=15.1904 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(92,5,1:6)=(/ 28.383000,28.565000,28.824000,28.932000,29.387000,30.482000 /)
XPIZA_LKT(92,5,1:6)=(/ 0.663601,0.707850,0.742668,0.880997,0.887652,0.923455 /)
XCGA_LKT(92,5,1:6)=(/ 0.904320,0.895023,0.888227,0.858657,0.859117,0.832103 /)
XEXT_COEFF_550_LKT(92,5)=28.704000 !rg=15.1904 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(92,6,1:6)=(/ 24.804000,24.907000,25.074000,25.259000,25.746000,26.437000 /)
XPIZA_LKT(92,6,1:6)=(/ 0.653979,0.696880,0.732702,0.877131,0.887100,0.917362 /)
XCGA_LKT(92,6,1:6)=(/ 0.906850,0.897040,0.889780,0.859643,0.859893,0.837260 /)
XEXT_COEFF_550_LKT(92,6)=25.135000 !rg=15.1904 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(92,7,1:6)=(/ 21.372000,21.467000,21.612000,21.831000,22.130000,22.728000 /)
XPIZA_LKT(92,7,1:6)=(/ 0.643485,0.685446,0.722010,0.872086,0.885883,0.910571 /)
XCGA_LKT(92,7,1:6)=(/ 0.909643,0.899447,0.891917,0.860773,0.861233,0.841277 /)
XEXT_COEFF_550_LKT(92,7)=21.650000 !rg=15.1904 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(92,8,1:6)=(/ 18.275000,18.349000,18.464000,18.603000,18.884000,19.316000 /)
XPIZA_LKT(92,8,1:6)=(/ 0.632983,0.673406,0.710436,0.865406,0.883914,0.905038 /)
XCGA_LKT(92,8,1:6)=(/ 0.912663,0.902190,0.894307,0.862357,0.862693,0.846107 /)
XEXT_COEFF_550_LKT(92,8)=18.478000 !rg=15.1904 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(92,9,1:6)=(/ 15.498000,15.568000,15.668000,15.727000,15.943000,16.306000 /)
XPIZA_LKT(92,9,1:6)=(/ 0.622366,0.661008,0.698242,0.857158,0.880768,0.900506 /)
XCGA_LKT(92,9,1:6)=(/ 0.915970,0.905327,0.896987,0.864093,0.864183,0.849197 /)
XEXT_COEFF_550_LKT(92,9)=15.643000 !rg=15.1904 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(92,10,1:6)=(/ 13.098000,13.147000,13.202000,13.305000,13.447000,13.751000 /)
XPIZA_LKT(92,10,1:6)=(/ 0.612624,0.648617,0.684465,0.848050,0.876171,0.896264 /)
XCGA_LKT(92,10,1:6)=(/ 0.919353,0.908557,0.899713,0.865993,0.865980,0.852310 /)
XEXT_COEFF_550_LKT(92,10)=13.204000 !rg=15.1904 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(92,11,1:6)=(/ 11.002000,11.037000,11.091000,11.135000,11.285000,11.494000 /)
XPIZA_LKT(92,11,1:6)=(/ 0.603266,0.636471,0.671506,0.836321,0.870090,0.891839 /)
XCGA_LKT(92,11,1:6)=(/ 0.922903,0.911983,0.902910,0.867590,0.868123,0.855757 /)
XEXT_COEFF_550_LKT(92,11)=11.084000 !rg=15.1904 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(92,12,1:6)=(/ 9.204300,9.236000,9.267100,9.329800,9.402600,9.570500 /)
XPIZA_LKT(92,12,1:6)=(/ 0.594579,0.624952,0.657851,0.823772,0.861823,0.887941 /)
XCGA_LKT(92,12,1:6)=(/ 0.926407,0.915673,0.906217,0.870293,0.869003,0.858580 /)
XEXT_COEFF_550_LKT(92,12)=9.262900 !rg=15.1904 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(92,13,1:6)=(/ 7.682400,7.703800,7.732600,7.768800,7.838100,7.965100 /)
XPIZA_LKT(92,13,1:6)=(/ 0.586884,0.613936,0.644912,0.809206,0.852946,0.882840 /)
XCGA_LKT(92,13,1:6)=(/ 0.929830,0.919340,0.909827,0.872697,0.872200,0.862103 /)
XEXT_COEFF_550_LKT(92,13)=7.728500 !rg=15.1904 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(92,14,1:6)=(/ 6.408900,6.424500,6.444400,6.476800,6.525200,6.616500 /)
XPIZA_LKT(92,14,1:6)=(/ 0.580059,0.604185,0.632601,0.793842,0.842099,0.877205 /)
XCGA_LKT(92,14,1:6)=(/ 0.933110,0.922927,0.913380,0.875403,0.873923,0.864093 /)
XEXT_COEFF_550_LKT(92,14)=6.446900 !rg=15.1904 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(92,15,1:6)=(/ 5.340600,5.352900,5.367000,5.386900,5.418200,5.493000 /)
XPIZA_LKT(92,15,1:6)=(/ 0.574124,0.595331,0.621173,0.777058,0.829389,0.869568 /)
XCGA_LKT(92,15,1:6)=(/ 0.936177,0.926590,0.917080,0.877917,0.875657,0.865563 /)
XEXT_COEFF_550_LKT(92,15)=5.365800 !rg=15.1904 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(92,16,1:6)=(/ 0.000000,4.448000,4.460200,4.474000,4.505100,4.560900 /)
XPIZA_LKT(92,16,1:6)=(/ 0.900,0.587292,0.610541,0.759589,0.815524,0.861180 /)
XCGA_LKT(92,16,1:6)=(/ 0.900,0.930100,0.920957,0.880933,0.878250,0.868197 /)
XEXT_COEFF_550_LKT(92,16)=4.458400 !rg=15.1904 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(92,17,1:6)=(/ 0.000000,3.701400,3.709900,3.724500,3.742100,3.782800 /)
XPIZA_LKT(92,17,1:6)=(/ 0.900,0.580341,0.600977,0.742099,0.800270,0.851864 /)
XCGA_LKT(92,17,1:6)=(/ 0.900,0.933397,0.924573,0.884207,0.880347,0.870947 /)
XEXT_COEFF_550_LKT(92,17)=3.709800 !rg=15.1904 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(92,18,1:6)=(/ 0.000000,3.080700,3.086600,3.095000,3.109300,3.138400 /)
XPIZA_LKT(92,18,1:6)=(/ 0.900,0.574380,0.592482,0.724137,0.783977,0.840253 /)
XCGA_LKT(92,18,1:6)=(/ 0.900,0.936470,0.928097,0.887630,0.883237,0.872823 /)
XEXT_COEFF_550_LKT(92,18)=3.086300 !rg=15.1904 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(92,19,1:6)=(/ 0.000000,0.000000,2.564500,2.571400,2.583800,2.607200 /)
XPIZA_LKT(92,19,1:6)=(/ 0.900,0.900,0.584903,0.706459,0.766921,0.827508 /)
XCGA_LKT(92,19,1:6)=(/ 0.900,0.900,0.931557,0.891307,0.886127,0.875200 /)
XEXT_COEFF_550_LKT(92,19)=2.564000 !rg=15.1904 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(92,20,1:6)=(/ 0.000000,0.000000,2.131500,2.137200,2.144900,2.160900 /)
XPIZA_LKT(92,20,1:6)=(/ 0.900,0.900,0.578275,0.689389,0.749311,0.813468 /)
XCGA_LKT(92,20,1:6)=(/ 0.900,0.900,0.934757,0.895220,0.889167,0.877860 /)
XEXT_COEFF_550_LKT(92,20)=2.131200 !rg=15.1904 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(93,1,1:6)=(/ 36.990000,36.816000,37.184000,37.812000,38.475000,37.410000 /)
XPIZA_LKT(93,1,1:6)=(/ 0.682672,0.724868,0.758878,0.890291,0.883452,0.936677 /)
XCGA_LKT(93,1,1:6)=(/ 0.900260,0.890493,0.884283,0.856220,0.854717,0.824933 /)
XEXT_COEFF_550_LKT(93,1)=37.725000 !rg=16.4566 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(93,2,1:6)=(/ 35.338000,35.497000,35.761000,35.937000,37.054000,37.972000 /)
XPIZA_LKT(93,2,1:6)=(/ 0.678898,0.723577,0.756582,0.888717,0.882396,0.936110 /)
XCGA_LKT(93,2,1:6)=(/ 0.901043,0.891690,0.885537,0.857160,0.860723,0.826177 /)
XEXT_COEFF_550_LKT(93,2)=35.601000 !rg=16.4566 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(93,3,1:6)=(/ 32.704000,32.937000,33.080000,33.453000,34.115000,35.242000 /)
XPIZA_LKT(93,3,1:6)=(/ 0.672695,0.718526,0.751555,0.885636,0.881763,0.930512 /)
XCGA_LKT(93,3,1:6)=(/ 0.902310,0.893113,0.886267,0.856583,0.859240,0.827983 /)
XEXT_COEFF_550_LKT(93,3)=33.046000 !rg=16.4566 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(93,4,1:6)=(/ 29.567000,29.775000,29.892000,30.265000,30.806000,31.724000 /)
XPIZA_LKT(93,4,1:6)=(/ 0.665176,0.710354,0.744755,0.882890,0.884076,0.925502 /)
XCGA_LKT(93,4,1:6)=(/ 0.904093,0.894733,0.887627,0.858197,0.859377,0.832247 /)
XEXT_COEFF_550_LKT(93,4)=29.828000 !rg=16.4566 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(93,5,1:6)=(/ 26.184000,26.302000,26.512000,26.733000,27.153000,27.991000 /)
XPIZA_LKT(93,5,1:6)=(/ 0.656297,0.700349,0.736502,0.879487,0.885735,0.918676 /)
XCGA_LKT(93,5,1:6)=(/ 0.906090,0.896343,0.889253,0.858790,0.860190,0.835997 /)
XEXT_COEFF_550_LKT(93,5)=26.547000 !rg=16.4566 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(93,6,1:6)=(/ 22.872000,23.007000,23.089000,23.309000,23.600000,24.194000 /)
XPIZA_LKT(93,6,1:6)=(/ 0.646820,0.690346,0.726158,0.874834,0.885913,0.912543 /)
XCGA_LKT(93,6,1:6)=(/ 0.908587,0.898643,0.890927,0.859930,0.861870,0.841663 /)
XEXT_COEFF_550_LKT(93,6)=23.108000 !rg=16.4566 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(93,7,1:6)=(/ 19.715000,19.822000,19.899000,20.104000,20.395000,20.861000 /)
XPIZA_LKT(93,7,1:6)=(/ 0.636676,0.678774,0.715391,0.869092,0.884860,0.907382 /)
XCGA_LKT(93,7,1:6)=(/ 0.911513,0.901203,0.893030,0.861853,0.862533,0.845640 /)
XEXT_COEFF_550_LKT(93,7)=19.885000 !rg=16.4566 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(93,8,1:6)=(/ 16.856000,16.936000,17.001000,17.167000,17.333000,17.738000 /)
XPIZA_LKT(93,8,1:6)=(/ 0.626510,0.666545,0.703535,0.861951,0.881938,0.901902 /)
XCGA_LKT(93,8,1:6)=(/ 0.914547,0.903970,0.895473,0.863740,0.863177,0.848427 /)
XEXT_COEFF_550_LKT(93,8)=16.992000 !rg=16.4566 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(93,9,1:6)=(/ 14.299000,14.349000,14.429000,14.524000,14.695000,15.026000 /)
XPIZA_LKT(93,9,1:6)=(/ 0.616582,0.653946,0.691027,0.852901,0.878449,0.897348 /)
XCGA_LKT(93,9,1:6)=(/ 0.917847,0.907027,0.898367,0.864993,0.865133,0.851567 /)
XEXT_COEFF_550_LKT(93,9)=14.446000 !rg=16.4566 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(93,10,1:6)=(/ 12.084000,12.124000,12.180000,12.259000,12.382000,12.611000 /)
XPIZA_LKT(93,10,1:6)=(/ 0.607298,0.642030,0.677933,0.842461,0.873521,0.894109 /)
XCGA_LKT(93,10,1:6)=(/ 0.921217,0.910283,0.901290,0.866547,0.867607,0.855357 /)
XEXT_COEFF_550_LKT(93,10)=12.167000 !rg=16.4566 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(93,11,1:6)=(/ 10.149000,10.181000,10.221000,10.287000,10.372000,10.599000 /)
XPIZA_LKT(93,11,1:6)=(/ 0.598449,0.630405,0.664694,0.830945,0.866323,0.889579 /)
XCGA_LKT(93,11,1:6)=(/ 0.924707,0.913807,0.904367,0.869230,0.868760,0.857073 /)
XEXT_COEFF_550_LKT(93,11)=10.222000 !rg=16.4566 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(93,12,1:6)=(/ 8.493400,8.518200,8.551900,8.586600,8.673100,8.815000 /)
XPIZA_LKT(93,12,1:6)=(/ 0.590382,0.619201,0.651649,0.817154,0.857902,0.885692 /)
XCGA_LKT(93,12,1:6)=(/ 0.928160,0.917457,0.907953,0.871150,0.870723,0.861230 /)
XEXT_COEFF_550_LKT(93,12)=8.543800 !rg=16.4566 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(93,13,1:6)=(/ 7.088200,7.106300,7.129700,7.163800,7.219600,7.337900 /)
XPIZA_LKT(93,13,1:6)=(/ 0.583126,0.608808,0.638663,0.802487,0.848424,0.880612 /)
XCGA_LKT(93,13,1:6)=(/ 0.931497,0.921137,0.911493,0.874053,0.873257,0.863220 /)
XEXT_COEFF_550_LKT(93,13)=7.130300 !rg=16.4566 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(93,14,1:6)=(/ 5.913700,5.928200,5.944200,5.971400,6.012500,6.097400 /)
XPIZA_LKT(93,14,1:6)=(/ 0.576790,0.599582,0.626923,0.786412,0.836644,0.874081 /)
XCGA_LKT(93,14,1:6)=(/ 0.934683,0.924747,0.915063,0.876543,0.875063,0.865507 /)
XEXT_COEFF_550_LKT(93,14)=5.943700 !rg=16.4566 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(93,15,1:6)=(/ 4.928000,4.939000,4.952400,4.969300,5.006500,5.067200 /)
XPIZA_LKT(93,15,1:6)=(/ 0.571281,0.591187,0.615979,0.769276,0.823624,0.866546 /)
XCGA_LKT(93,15,1:6)=(/ 0.937647,0.928317,0.918913,0.879227,0.877280,0.867777 /)
XEXT_COEFF_550_LKT(93,15)=4.951100 !rg=16.4566 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(93,16,1:6)=(/ 0.000000,4.104200,4.113900,4.129800,4.148200,4.200900 /)
XPIZA_LKT(93,16,1:6)=(/ 0.900,0.583709,0.605697,0.751790,0.808617,0.857191 /)
XCGA_LKT(93,16,1:6)=(/ 0.900,0.931700,0.922653,0.882487,0.878917,0.869600 /)
XEXT_COEFF_550_LKT(93,16)=4.113300 !rg=16.4566 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(93,17,1:6)=(/ 0.000000,3.415700,3.422500,3.433500,3.452300,3.487300 /)
XPIZA_LKT(93,17,1:6)=(/ 0.900,0.577254,0.596625,0.733879,0.793250,0.846897 /)
XCGA_LKT(93,17,1:6)=(/ 0.900,0.934883,0.926277,0.885803,0.881793,0.872100 /)
XEXT_COEFF_550_LKT(93,17)=3.422900 !rg=16.4566 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(93,18,1:6)=(/ 0.000000,0.000000,2.848400,2.855500,2.871600,2.897300 /)
XPIZA_LKT(93,18,1:6)=(/ 0.900,0.900,0.588706,0.716084,0.776646,0.835169 /)
XCGA_LKT(93,18,1:6)=(/ 0.900,0.900,0.929773,0.889260,0.884640,0.874760 /)
XEXT_COEFF_550_LKT(93,18)=2.847800 !rg=16.4566 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(93,19,1:6)=(/ 0.000000,0.000000,2.366200,2.372500,2.381700,2.400900 /)
XPIZA_LKT(93,19,1:6)=(/ 0.900,0.900,0.581558,0.698580,0.759082,0.821498 /)
XCGA_LKT(93,19,1:6)=(/ 0.900,0.900,0.933073,0.893113,0.887473,0.876490 /)
XEXT_COEFF_550_LKT(93,19)=2.365800 !rg=16.4566 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(93,20,1:6)=(/ 0.000000,0.000000,1.966700,1.971200,1.978900,1.993200 /)
XPIZA_LKT(93,20,1:6)=(/ 0.900,0.900,0.575305,0.681725,0.741253,0.806839 /)
XCGA_LKT(93,20,1:6)=(/ 0.900,0.900,0.936200,0.897077,0.890600,0.879137 /)
XEXT_COEFF_550_LKT(93,20)=0.000000 !rg=16.4566 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(94,1,1:6)=(/ 34.106000,34.325000,34.362000,34.877000,35.458000,37.678000 /)
XPIZA_LKT(94,1,1:6)=(/ 0.675023,0.720843,0.754795,0.890264,0.871234,0.934095 /)
XCGA_LKT(94,1,1:6)=(/ 0.902353,0.892447,0.885840,0.858617,0.861243,0.835743 /)
XEXT_COEFF_550_LKT(94,1)=34.364000 !rg=17.8284 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(94,2,1:6)=(/ 32.585000,32.794000,32.908000,33.284000,33.875000,34.584000 /)
XPIZA_LKT(94,2,1:6)=(/ 0.671221,0.716895,0.750815,0.887905,0.875266,0.932553 /)
XCGA_LKT(94,2,1:6)=(/ 0.902867,0.893573,0.886547,0.857240,0.862053,0.834193 /)
XEXT_COEFF_550_LKT(94,2)=32.910000 !rg=17.8284 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(94,3,1:6)=(/ 30.182000,30.276000,30.586000,30.855000,31.186000,32.334000 /)
XPIZA_LKT(94,3,1:6)=(/ 0.665277,0.710199,0.746586,0.884864,0.880967,0.926068 /)
XCGA_LKT(94,3,1:6)=(/ 0.904023,0.894220,0.887613,0.857343,0.861897,0.833723 /)
XEXT_COEFF_550_LKT(94,3)=30.601000 !rg=17.8284 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(94,4,1:6)=(/ 27.271000,27.378000,27.649000,27.838000,28.201000,29.143000 /)
XPIZA_LKT(94,4,1:6)=(/ 0.657789,0.702418,0.739658,0.881610,0.883986,0.920186 /)
XCGA_LKT(94,4,1:6)=(/ 0.905720,0.895833,0.889033,0.858313,0.861973,0.837003 /)
XEXT_COEFF_550_LKT(94,4)=27.621000 !rg=17.8284 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(94,5,1:6)=(/ 24.155000,24.292000,24.396000,24.652000,25.082000,25.790000 /)
XPIZA_LKT(94,5,1:6)=(/ 0.649011,0.693622,0.730261,0.877337,0.883759,0.913527 /)
XCGA_LKT(94,5,1:6)=(/ 0.907923,0.897893,0.890377,0.859700,0.860773,0.839783 /)
XEXT_COEFF_550_LKT(94,5)=24.373000 !rg=17.8284 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(94,6,1:6)=(/ 21.103000,21.190000,21.347000,21.468000,21.878000,22.450000 /)
XPIZA_LKT(94,6,1:6)=(/ 0.639906,0.682677,0.720735,0.871903,0.884924,0.908316 /)
XCGA_LKT(94,6,1:6)=(/ 0.910447,0.900067,0.892490,0.860657,0.862170,0.844960 /)
XEXT_COEFF_550_LKT(94,6)=21.308000 !rg=17.8284 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(94,7,1:6)=(/ 18.190000,18.266000,18.386000,18.488000,18.764000,19.217000 /)
XPIZA_LKT(94,7,1:6)=(/ 0.630197,0.671501,0.709517,0.865592,0.883608,0.902428 /)
XCGA_LKT(94,7,1:6)=(/ 0.913363,0.902663,0.894727,0.862160,0.863943,0.847763 /)
XEXT_COEFF_550_LKT(94,7)=18.368000 !rg=17.8284 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(94,8,1:6)=(/ 15.553000,15.612000,15.697000,15.767000,15.995000,16.359000 /)
XPIZA_LKT(94,8,1:6)=(/ 0.620510,0.659517,0.697115,0.857605,0.880223,0.898703 /)
XCGA_LKT(94,8,1:6)=(/ 0.916417,0.905547,0.897117,0.863957,0.865537,0.851953 /)
XEXT_COEFF_550_LKT(94,8)=15.676000 !rg=17.8284 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(94,9,1:6)=(/ 13.194000,13.247000,13.290000,13.395000,13.530000,13.861000 /)
XPIZA_LKT(94,9,1:6)=(/ 0.611087,0.647667,0.683923,0.848216,0.875740,0.894484 /)
XCGA_LKT(94,9,1:6)=(/ 0.919730,0.908790,0.899820,0.865977,0.865880,0.853270 /)
XEXT_COEFF_550_LKT(94,9)=13.294000 !rg=17.8284 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(94,10,1:6)=(/ 11.148000,11.185000,11.226000,11.300000,11.419000,11.639000 /)
XPIZA_LKT(94,10,1:6)=(/ 0.602224,0.635768,0.670975,0.837496,0.870095,0.891551 /)
XCGA_LKT(94,10,1:6)=(/ 0.923063,0.912083,0.902817,0.867963,0.867723,0.856687 /)
XEXT_COEFF_550_LKT(94,10)=11.232000 !rg=17.8284 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(94,11,1:6)=(/ 9.362900,9.391900,9.426400,9.469800,9.558200,9.731500 /)
XPIZA_LKT(94,11,1:6)=(/ 0.593891,0.624499,0.658061,0.824681,0.862815,0.887469 /)
XCGA_LKT(94,11,1:6)=(/ 0.926487,0.915593,0.906123,0.870033,0.869540,0.858763 /)
XEXT_COEFF_550_LKT(94,11)=9.422600 !rg=17.8284 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(94,12,1:6)=(/ 7.836000,7.857300,7.885000,7.928800,7.990100,8.123700 /)
XPIZA_LKT(94,12,1:6)=(/ 0.586331,0.613751,0.645161,0.810744,0.853947,0.883156 /)
XCGA_LKT(94,12,1:6)=(/ 0.929890,0.919237,0.909600,0.872670,0.871837,0.862113 /)
XEXT_COEFF_550_LKT(94,12)=7.889000 !rg=17.8284 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(94,13,1:6)=(/ 6.540200,6.557600,6.575500,6.605500,6.653900,6.756900 /)
XPIZA_LKT(94,13,1:6)=(/ 0.579620,0.603992,0.632709,0.795174,0.843084,0.877255 /)
XCGA_LKT(94,13,1:6)=(/ 0.933123,0.922930,0.913247,0.875030,0.873640,0.864317 /)
XEXT_COEFF_550_LKT(94,13)=6.576200 !rg=17.8284 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(94,14,1:6)=(/ 5.456900,5.469900,5.485400,5.503400,5.551400,5.627400 /)
XPIZA_LKT(94,14,1:6)=(/ 0.573764,0.595174,0.621466,0.778703,0.831317,0.870693 /)
XCGA_LKT(94,14,1:6)=(/ 0.936200,0.926523,0.916943,0.877643,0.876113,0.866520 /)
XEXT_COEFF_550_LKT(94,14)=5.486500 !rg=17.8284 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(94,15,1:6)=(/ 0.000000,4.557100,4.568700,4.586000,4.610200,4.665100 /)
XPIZA_LKT(94,15,1:6)=(/ 0.900,0.587295,0.610851,0.761691,0.817184,0.862188 /)
XCGA_LKT(94,15,1:6)=(/ 0.900,0.929997,0.920697,0.880720,0.878017,0.868040 /)
XEXT_COEFF_550_LKT(94,15)=4.567600 !rg=17.8284 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(94,16,1:6)=(/ 0.000000,3.787100,3.795400,3.808100,3.827700,3.870800 /)
XPIZA_LKT(94,16,1:6)=(/ 0.900,0.580263,0.601182,0.743611,0.801910,0.852637 /)
XCGA_LKT(94,16,1:6)=(/ 0.900,0.933317,0.924353,0.883937,0.880420,0.870477 /)
XEXT_COEFF_550_LKT(94,16)=3.794700 !rg=17.8284 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(94,17,1:6)=(/ 0.000000,3.152000,3.158200,3.166500,3.185600,3.217200 /)
XPIZA_LKT(94,17,1:6)=(/ 0.900,0.574269,0.592599,0.725669,0.786007,0.841994 /)
XCGA_LKT(94,17,1:6)=(/ 0.900,0.936400,0.927950,0.887333,0.883220,0.873530 /)
XEXT_COEFF_550_LKT(94,17)=3.158600 !rg=17.8284 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(94,18,1:6)=(/ 0.000000,0.000000,2.628100,2.635000,2.646700,2.667000 /)
XPIZA_LKT(94,18,1:6)=(/ 0.900,0.900,0.585091,0.708165,0.768939,0.828997 /)
XCGA_LKT(94,18,1:6)=(/ 0.900,0.900,0.931363,0.890967,0.885933,0.875117 /)
XEXT_COEFF_550_LKT(94,18)=2.627500 !rg=17.8284 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(94,19,1:6)=(/ 0.000000,0.000000,2.183400,2.187900,2.197500,2.214400 /)
XPIZA_LKT(94,19,1:6)=(/ 0.900,0.900,0.578370,0.690717,0.751146,0.815071 /)
XCGA_LKT(94,19,1:6)=(/ 0.900,0.900,0.934610,0.894807,0.888833,0.877463 /)
XEXT_COEFF_550_LKT(94,19)=2.182900 !rg=17.8284 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(94,20,1:6)=(/ 0.000000,0.000000,0.000000,1.818400,1.826000,1.839800 /)
XPIZA_LKT(94,20,1:6)=(/ 0.900,0.900,0.900,0.674164,0.733394,0.800154 /)
XCGA_LKT(94,20,1:6)=(/ 0.900,0.900,0.900,0.898920,0.892233,0.880450 /)
XEXT_COEFF_550_LKT(94,20)=0.000000 !rg=17.8284 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(95,1,1:6)=(/ 31.303000,31.603000,31.656000,31.684000,32.339000,33.635000 /)
XPIZA_LKT(95,1,1:6)=(/ 0.666339,0.713413,0.747991,0.887778,0.862005,0.926199 /)
XCGA_LKT(95,1,1:6)=(/ 0.903757,0.894873,0.886507,0.857707,0.865080,0.829387 /)
XEXT_COEFF_550_LKT(95,1)=31.342000 !rg=19.3145 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(95,2,1:6)=(/ 30.016000,30.174000,30.497000,30.619000,31.335000,32.139000 /)
XPIZA_LKT(95,2,1:6)=(/ 0.662946,0.709479,0.746376,0.886564,0.873934,0.927228 /)
XCGA_LKT(95,2,1:6)=(/ 0.904463,0.894653,0.888113,0.858617,0.864033,0.835853 /)
XEXT_COEFF_550_LKT(95,2)=30.321000 !rg=19.3145 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(95,3,1:6)=(/ 27.817000,28.004000,28.066000,28.398000,28.906000,29.789000 /)
XPIZA_LKT(95,3,1:6)=(/ 0.657330,0.703430,0.739426,0.883521,0.878845,0.920839 /)
XCGA_LKT(95,3,1:6)=(/ 0.905730,0.896140,0.888270,0.858410,0.863213,0.835703 /)
XEXT_COEFF_550_LKT(95,3)=28.161000 !rg=19.3145 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(95,4,1:6)=(/ 25.142000,25.296000,25.373000,25.661000,26.113000,26.892000 /)
XPIZA_LKT(95,4,1:6)=(/ 0.650087,0.695585,0.732450,0.879639,0.882169,0.914394 /)
XCGA_LKT(95,4,1:6)=(/ 0.907533,0.897580,0.889597,0.859567,0.862717,0.838950 /)
XEXT_COEFF_550_LKT(95,4)=25.442000 !rg=19.3145 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(95,5,1:6)=(/ 22.290000,22.374000,22.544000,22.690000,22.975000,23.639000 /)
XPIZA_LKT(95,5,1:6)=(/ 0.642067,0.685780,0.724565,0.874823,0.884581,0.908888 /)
XCGA_LKT(95,5,1:6)=(/ 0.909750,0.899357,0.891723,0.860337,0.863117,0.843880 /)
XEXT_COEFF_550_LKT(95,5)=22.554000 !rg=19.3145 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(95,6,1:6)=(/ 19.460000,19.548000,19.640000,19.886000,20.074000,20.652000 /)
XPIZA_LKT(95,6,1:6)=(/ 0.632943,0.675552,0.713669,0.869326,0.883511,0.904264 /)
XCGA_LKT(95,6,1:6)=(/ 0.912327,0.901680,0.893480,0.862110,0.862893,0.846957 /)
XEXT_COEFF_550_LKT(95,6)=19.663000 !rg=19.3145 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(95,7,1:6)=(/ 16.774000,16.842000,16.923000,17.065000,17.269000,17.766000 /)
XPIZA_LKT(95,7,1:6)=(/ 0.623782,0.664192,0.702528,0.861980,0.881652,0.898831 /)
XCGA_LKT(95,7,1:6)=(/ 0.915183,0.904377,0.895720,0.863857,0.865060,0.849417 /)
XEXT_COEFF_550_LKT(95,7)=16.934000 !rg=19.3145 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(95,8,1:6)=(/ 14.347000,14.397000,14.465000,14.570000,14.749000,15.050000 /)
XPIZA_LKT(95,8,1:6)=(/ 0.614684,0.652512,0.690163,0.853531,0.878690,0.896435 /)
XCGA_LKT(95,8,1:6)=(/ 0.918287,0.907240,0.898430,0.865657,0.866250,0.853440 /)
XEXT_COEFF_550_LKT(95,8)=14.474000 !rg=19.3145 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(95,9,1:6)=(/ 12.173000,12.213000,12.274000,12.335000,12.474000,12.706000 /)
XPIZA_LKT(95,9,1:6)=(/ 0.605759,0.640815,0.677625,0.843126,0.873944,0.893004 /)
XCGA_LKT(95,9,1:6)=(/ 0.921610,0.910577,0.901437,0.866983,0.867877,0.856160 /)
XEXT_COEFF_550_LKT(95,9)=12.269000 !rg=19.3145 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(95,10,1:6)=(/ 10.286000,10.321000,10.360000,10.398000,10.496000,10.705000 /)
XPIZA_LKT(95,10,1:6)=(/ 0.597442,0.629662,0.664455,0.831395,0.867067,0.889325 /)
XCGA_LKT(95,10,1:6)=(/ 0.924897,0.913930,0.904467,0.868680,0.869013,0.858513 /)
XEXT_COEFF_550_LKT(95,10)=10.357000 !rg=19.3145 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(95,11,1:6)=(/ 8.639500,8.663200,8.695000,8.742400,8.836200,8.978200 /)
XPIZA_LKT(95,11,1:6)=(/ 0.589638,0.618748,0.651534,0.818511,0.859028,0.885187 /)
XCGA_LKT(95,11,1:6)=(/ 0.928277,0.917413,0.907790,0.871250,0.870793,0.860817 /)
XEXT_COEFF_550_LKT(95,11)=8.699700 !rg=19.3145 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(95,12,1:6)=(/ 7.230100,7.250200,7.270600,7.305000,7.362800,7.476700 /)
XPIZA_LKT(95,12,1:6)=(/ 0.582589,0.608665,0.638779,0.803801,0.849366,0.880484 /)
XCGA_LKT(95,12,1:6)=(/ 0.931567,0.921083,0.911327,0.873690,0.872920,0.863413 /)
XEXT_COEFF_550_LKT(95,12)=7.271600 !rg=19.3145 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(95,13,1:6)=(/ 6.035400,6.050100,6.068100,6.087600,6.143300,6.227200 /)
XPIZA_LKT(95,13,1:6)=(/ 0.576381,0.599291,0.627079,0.787639,0.838242,0.874476 /)
XCGA_LKT(95,13,1:6)=(/ 0.934727,0.924753,0.915037,0.876203,0.875103,0.865857 /)
XEXT_COEFF_550_LKT(95,13)=6.068300 !rg=19.3145 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(95,14,1:6)=(/ 5.035300,5.046800,5.060200,5.081000,5.105500,5.166600 /)
XPIZA_LKT(95,14,1:6)=(/ 0.570956,0.591008,0.616066,0.771061,0.824662,0.866786 /)
XCGA_LKT(95,14,1:6)=(/ 0.937653,0.928250,0.918760,0.879157,0.876487,0.866963 /)
XEXT_COEFF_550_LKT(95,14)=5.058100 !rg=19.3145 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(95,15,1:6)=(/ 0.000000,4.204400,4.215100,4.228400,4.251600,4.299500 /)
XPIZA_LKT(95,15,1:6)=(/ 0.900,0.583631,0.606008,0.753375,0.810550,0.858436 /)
XCGA_LKT(95,15,1:6)=(/ 0.900,0.931603,0.922477,0.882183,0.878973,0.869823 /)
XEXT_COEFF_550_LKT(95,15)=4.213300 !rg=19.3145 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(95,16,1:6)=(/ 0.000000,3.494500,3.501900,3.512800,3.533400,3.568900 /)
XPIZA_LKT(95,16,1:6)=(/ 0.900,0.577066,0.596773,0.735417,0.795061,0.848437 /)
XCGA_LKT(95,16,1:6)=(/ 0.900,0.934853,0.926120,0.885510,0.881727,0.872497 /)
XEXT_COEFF_550_LKT(95,16)=3.502200 !rg=19.3145 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(95,17,1:6)=(/ 0.000000,0.000000,2.914100,2.922300,2.934600,2.958800 /)
XPIZA_LKT(95,17,1:6)=(/ 0.900,0.900,0.588739,0.717577,0.778354,0.836321 /)
XCGA_LKT(95,17,1:6)=(/ 0.900,0.900,0.929630,0.889000,0.884317,0.874313 /)
XEXT_COEFF_550_LKT(95,17)=2.913500 !rg=19.3145 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(95,18,1:6)=(/ 0.000000,0.000000,2.425000,2.430600,2.440500,2.459600 /)
XPIZA_LKT(95,18,1:6)=(/ 0.900,0.900,0.581682,0.700061,0.760979,0.823260 /)
XCGA_LKT(95,18,1:6)=(/ 0.900,0.900,0.932923,0.892673,0.887207,0.876537 /)
XEXT_COEFF_550_LKT(95,18)=2.424200 !rg=19.3145 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(95,19,1:6)=(/ 0.000000,0.000000,2.014600,2.019000,2.027400,2.042000 /)
XPIZA_LKT(95,19,1:6)=(/ 0.900,0.900,0.575408,0.683069,0.743324,0.808939 /)
XCGA_LKT(95,19,1:6)=(/ 0.900,0.900,0.936043,0.896693,0.890533,0.879180 /)
XEXT_COEFF_550_LKT(95,19)=2.014600 !rg=19.3145 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(95,20,1:6)=(/ 0.000000,0.000000,0.000000,1.678300,1.684000,1.693900 /)
XPIZA_LKT(95,20,1:6)=(/ 0.900,0.900,0.900,0.667020,0.725198,0.792918 /)
XCGA_LKT(95,20,1:6)=(/ 0.900,0.900,0.900,0.900807,0.893517,0.881467 /)
XEXT_COEFF_550_LKT(95,20)=0.000000 !rg=19.3145 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(96,1,1:6)=(/ 28.781000,28.957000,29.243000,29.824000,30.094000,30.629000 /)
XPIZA_LKT(96,1,1:6)=(/ 0.657201,0.704540,0.742826,0.886452,0.865611,0.922810 /)
XCGA_LKT(96,1,1:6)=(/ 0.905170,0.895800,0.888413,0.858210,0.867473,0.832000 /)
XEXT_COEFF_550_LKT(96,1)=29.679000 !rg=20.9245 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(96,2,1:6)=(/ 27.646000,27.793000,27.904000,28.301000,28.686000,29.738000 /)
XPIZA_LKT(96,2,1:6)=(/ 0.654698,0.701696,0.738663,0.884713,0.874103,0.921754 /)
XCGA_LKT(96,2,1:6)=(/ 0.906087,0.895990,0.888333,0.858043,0.864420,0.838807 /)
XEXT_COEFF_550_LKT(96,2)=28.070000 !rg=20.9245 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(96,3,1:6)=(/ 25.639000,25.786000,26.040000,26.139000,26.445000,27.217000 /)
XPIZA_LKT(96,3,1:6)=(/ 0.649369,0.695640,0.734847,0.881583,0.880914,0.915586 /)
XCGA_LKT(96,3,1:6)=(/ 0.907427,0.897413,0.890333,0.859977,0.864340,0.840247 /)
XEXT_COEFF_550_LKT(96,3)=25.917000 !rg=20.9245 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(96,4,1:6)=(/ 23.188000,23.317000,23.499000,23.567000,23.893000,24.547000 /)
XPIZA_LKT(96,4,1:6)=(/ 0.642696,0.687819,0.727198,0.877433,0.883243,0.909260 /)
XCGA_LKT(96,4,1:6)=(/ 0.909267,0.899140,0.891490,0.860807,0.863450,0.843233 /)
XEXT_COEFF_550_LKT(96,4)=23.402000 !rg=20.9245 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(96,5,1:6)=(/ 20.554000,20.664000,20.725000,20.980000,21.277000,21.821000 /)
XPIZA_LKT(96,5,1:6)=(/ 0.634937,0.678621,0.716993,0.872216,0.883669,0.903260 /)
XCGA_LKT(96,5,1:6)=(/ 0.911610,0.901167,0.892720,0.861787,0.864337,0.846187 /)
XEXT_COEFF_550_LKT(96,5)=20.752000 !rg=20.9245 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(96,6,1:6)=(/ 17.948000,18.030000,18.130000,18.251000,18.529000,19.009000 /)
XPIZA_LKT(96,6,1:6)=(/ 0.626348,0.668049,0.707365,0.865739,0.882667,0.899866 /)
XCGA_LKT(96,6,1:6)=(/ 0.914167,0.903447,0.894990,0.862993,0.864060,0.848903 /)
XEXT_COEFF_550_LKT(96,6)=18.103000 !rg=20.9245 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(96,7,1:6)=(/ 15.473000,15.534000,15.605000,15.687000,15.883000,16.239000 /)
XPIZA_LKT(96,7,1:6)=(/ 0.617563,0.657066,0.695474,0.857678,0.880101,0.896101 /)
XCGA_LKT(96,7,1:6)=(/ 0.917100,0.906093,0.897357,0.864340,0.865280,0.851850 /)
XEXT_COEFF_550_LKT(96,7)=15.590000 !rg=20.9245 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(96,8,1:6)=(/ 13.234000,13.286000,13.335000,13.437000,13.568000,13.879000 /)
XPIZA_LKT(96,8,1:6)=(/ 0.608982,0.645804,0.683074,0.848755,0.876149,0.893684 /)
XCGA_LKT(96,8,1:6)=(/ 0.920203,0.909130,0.899890,0.866527,0.867293,0.855527 /)
XEXT_COEFF_550_LKT(96,8)=13.334000 !rg=20.9245 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(96,9,1:6)=(/ 11.229000,11.271000,11.308000,11.408000,11.497000,11.736000 /)
XPIZA_LKT(96,9,1:6)=(/ 0.600615,0.634588,0.670267,0.838084,0.870576,0.890138 /)
XCGA_LKT(96,9,1:6)=(/ 0.923470,0.912397,0.902970,0.868410,0.868017,0.857970 /)
XEXT_COEFF_550_LKT(96,9)=11.310000 !rg=20.9245 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(96,10,1:6)=(/ 9.490500,9.520700,9.558800,9.605000,9.712500,9.858200 /)
XPIZA_LKT(96,10,1:6)=(/ 0.592887,0.623748,0.657815,0.825422,0.863654,0.887170 /)
XCGA_LKT(96,10,1:6)=(/ 0.926720,0.915753,0.906213,0.870173,0.869757,0.860577 /)
XEXT_COEFF_550_LKT(96,10)=9.556800 !rg=20.9245 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(96,11,1:6)=(/ 7.972100,7.993900,8.018100,8.060500,8.135000,8.268000 /)
XPIZA_LKT(96,11,1:6)=(/ 0.585688,0.613345,0.645054,0.811753,0.854694,0.883011 /)
XCGA_LKT(96,11,1:6)=(/ 0.929993,0.919273,0.909490,0.872543,0.871353,0.862230 /)
XEXT_COEFF_550_LKT(96,11)=8.016400 !rg=20.9245 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(96,12,1:6)=(/ 6.671700,6.689100,6.709400,6.735700,6.800600,6.898300 /)
XPIZA_LKT(96,12,1:6)=(/ 0.579095,0.603669,0.632968,0.796449,0.844721,0.877486 /)
XCGA_LKT(96,12,1:6)=(/ 0.933217,0.922910,0.913150,0.874860,0.873863,0.864330 /)
XEXT_COEFF_550_LKT(96,12)=6.710600 !rg=20.9245 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(96,13,1:6)=(/ 5.568800,5.582200,5.597400,5.624400,5.644100,5.715200 /)
XPIZA_LKT(96,13,1:6)=(/ 0.573297,0.594909,0.621415,0.780399,0.832101,0.870647 /)
XCGA_LKT(96,13,1:6)=(/ 0.936273,0.926487,0.916880,0.877723,0.875447,0.865927 /)
XEXT_COEFF_550_LKT(96,13)=5.596900 !rg=20.9245 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(96,14,1:6)=(/ 0.000000,4.656100,4.668700,4.682400,4.710600,4.769600 /)
XPIZA_LKT(96,14,1:6)=(/ 0.900,0.587081,0.610961,0.762876,0.818835,0.863303 /)
XCGA_LKT(96,14,1:6)=(/ 0.900,0.929897,0.920567,0.880570,0.878293,0.868997 /)
XEXT_COEFF_550_LKT(96,14)=4.666200 !rg=20.9245 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(96,15,1:6)=(/ 0.000000,3.879400,3.888100,3.901700,3.922800,3.965600 /)
XPIZA_LKT(96,15,1:6)=(/ 0.900,0.580170,0.601331,0.745397,0.803817,0.854000 /)
XCGA_LKT(96,15,1:6)=(/ 0.900,0.933200,0.924183,0.883657,0.879877,0.870827 /)
XEXT_COEFF_550_LKT(96,15)=3.888100 !rg=20.9245 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(96,16,1:6)=(/ 0.000000,3.224700,3.231000,3.240100,3.256800,3.287400 /)
XPIZA_LKT(96,16,1:6)=(/ 0.900,0.574123,0.592622,0.727194,0.787580,0.843127 /)
XCGA_LKT(96,16,1:6)=(/ 0.900,0.936330,0.927813,0.887037,0.882733,0.872997 /)
XEXT_COEFF_550_LKT(96,16)=3.230400 !rg=20.9245 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(96,17,1:6)=(/ 0.000000,0.000000,2.689000,2.695200,2.708700,2.732800 /)
XPIZA_LKT(96,17,1:6)=(/ 0.900,0.900,0.585140,0.709413,0.770578,0.830792 /)
XCGA_LKT(96,17,1:6)=(/ 0.900,0.900,0.931230,0.890633,0.885510,0.875370 /)
XEXT_COEFF_550_LKT(96,17)=2.688100 !rg=20.9245 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(96,18,1:6)=(/ 0.000000,0.000000,2.237400,2.243100,2.251900,2.269800 /)
XPIZA_LKT(96,18,1:6)=(/ 0.900,0.900,0.578450,0.692342,0.753124,0.816926 /)
XCGA_LKT(96,18,1:6)=(/ 0.900,0.900,0.934440,0.894500,0.888463,0.877417 /)
XEXT_COEFF_550_LKT(96,18)=2.237000 !rg=20.9245 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(96,19,1:6)=(/ 0.000000,0.000000,0.000000,1.862900,1.869500,1.881600 /)
XPIZA_LKT(96,19,1:6)=(/ 0.900,0.900,0.900,0.675591,0.735051,0.801841 /)
XCGA_LKT(96,19,1:6)=(/ 0.900,0.900,0.900,0.898507,0.891747,0.879957 /)
XEXT_COEFF_550_LKT(96,19)=0.000000 !rg=20.9245 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(96,20,1:6)=(/ 0.000000,0.000000,0.000000,1.548400,1.553600,1.563500 /)
XPIZA_LKT(96,20,1:6)=(/ 0.900,0.900,0.900,0.659858,0.717186,0.785783 /)
XCGA_LKT(96,20,1:6)=(/ 0.900,0.900,0.900,0.902737,0.895133,0.882873 /)
XEXT_COEFF_550_LKT(96,20)=0.000000 !rg=20.9245 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(97,1,1:6)=(/ 26.588000,26.632000,26.860000,27.202000,27.567000,28.750000 /)
XPIZA_LKT(97,1,1:6)=(/ 0.650027,0.695955,0.735528,0.883836,0.878720,0.920984 /)
XCGA_LKT(97,1,1:6)=(/ 0.906587,0.896570,0.889700,0.860387,0.865810,0.842993 /)
XEXT_COEFF_550_LKT(97,1)=26.812000 !rg=22.6687 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(97,2,1:6)=(/ 25.507000,25.638000,25.894000,26.130000,26.525000,27.275000 /)
XPIZA_LKT(97,2,1:6)=(/ 0.646896,0.693682,0.734075,0.881628,0.880093,0.915774 /)
XCGA_LKT(97,2,1:6)=(/ 0.907873,0.897733,0.890330,0.859663,0.864847,0.838803 /)
XEXT_COEFF_550_LKT(97,2)=25.694000 !rg=22.6687 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(97,3,1:6)=(/ 23.662000,23.785000,23.927000,24.123000,24.438000,25.067000 /)
XPIZA_LKT(97,3,1:6)=(/ 0.642060,0.688049,0.727364,0.878851,0.882350,0.909977 /)
XCGA_LKT(97,3,1:6)=(/ 0.909297,0.898943,0.891187,0.860260,0.864523,0.844173 /)
XEXT_COEFF_550_LKT(97,3)=23.988000 !rg=22.6687 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(97,4,1:6)=(/ 21.396000,21.478000,21.628000,21.784000,22.071000,22.635000 /)
XPIZA_LKT(97,4,1:6)=(/ 0.635531,0.679714,0.719926,0.874452,0.883504,0.903369 /)
XCGA_LKT(97,4,1:6)=(/ 0.911270,0.900723,0.892707,0.861107,0.863910,0.846623 /)
XEXT_COEFF_550_LKT(97,4)=21.668000 !rg=22.6687 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(97,5,1:6)=(/ 18.955000,19.039000,19.181000,19.254000,19.481000,19.924000 /)
XPIZA_LKT(97,5,1:6)=(/ 0.627997,0.670921,0.711037,0.868599,0.883598,0.899567 /)
XCGA_LKT(97,5,1:6)=(/ 0.913463,0.902703,0.894610,0.862523,0.864683,0.850050 /)
XEXT_COEFF_550_LKT(97,5)=19.132000 !rg=22.6687 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(97,6,1:6)=(/ 16.563000,16.630000,16.713000,16.813000,17.040000,17.431000 /)
XPIZA_LKT(97,6,1:6)=(/ 0.620133,0.660728,0.700044,0.861561,0.881678,0.896312 /)
XCGA_LKT(97,6,1:6)=(/ 0.916100,0.905193,0.896510,0.863867,0.865110,0.852667 /)
XEXT_COEFF_550_LKT(97,6)=16.723000 !rg=22.6687 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(97,7,1:6)=(/ 14.277000,14.327000,14.391000,14.499000,14.678000,14.981000 /)
XPIZA_LKT(97,7,1:6)=(/ 0.611770,0.650025,0.688433,0.853408,0.878236,0.893066 /)
XCGA_LKT(97,7,1:6)=(/ 0.919047,0.907823,0.898833,0.865437,0.866430,0.854627 /)
XEXT_COEFF_550_LKT(97,7)=14.412000 !rg=22.6687 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(97,8,1:6)=(/ 12.212000,12.253000,12.309000,12.357000,12.523000,12.765000 /)
XPIZA_LKT(97,8,1:6)=(/ 0.603721,0.639132,0.676316,0.843311,0.873989,0.890947 /)
XCGA_LKT(97,8,1:6)=(/ 0.922097,0.910917,0.901670,0.866977,0.868167,0.857283 /)
XEXT_COEFF_550_LKT(97,8)=12.312000 !rg=22.6687 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(97,9,1:6)=(/ 10.360000,10.393000,10.442000,10.481000,10.587000,10.764000 /)
XPIZA_LKT(97,9,1:6)=(/ 0.595803,0.628188,0.663825,0.831911,0.867768,0.888465 /)
XCGA_LKT(97,9,1:6)=(/ 0.925317,0.914223,0.904747,0.869047,0.869303,0.859440 /)
XEXT_COEFF_550_LKT(97,9)=10.430000 !rg=22.6687 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(97,10,1:6)=(/ 8.756200,8.783900,8.814000,8.864500,8.931000,9.082600 /)
XPIZA_LKT(97,10,1:6)=(/ 0.588589,0.618004,0.651226,0.819419,0.859661,0.884930 /)
XCGA_LKT(97,10,1:6)=(/ 0.928507,0.917623,0.907863,0.871420,0.870493,0.861250 /)
XEXT_COEFF_550_LKT(97,10)=8.806900 !rg=22.6687 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(97,11,1:6)=(/ 7.355400,7.375300,7.399200,7.428300,7.503600,7.616500 /)
XPIZA_LKT(97,11,1:6)=(/ 0.581866,0.608020,0.638937,0.804731,0.850669,0.880737 /)
XCGA_LKT(97,11,1:6)=(/ 0.931720,0.921163,0.911293,0.873560,0.873217,0.863873 /)
XEXT_COEFF_550_LKT(97,11)=7.398300 !rg=22.6687 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(97,12,1:6)=(/ 6.156200,6.171800,6.188800,6.217200,6.250600,6.340700 /)
XPIZA_LKT(97,12,1:6)=(/ 0.575814,0.598943,0.627048,0.789280,0.838978,0.874495 /)
XCGA_LKT(97,12,1:6)=(/ 0.934803,0.924720,0.914927,0.876280,0.874437,0.865123 /)
XEXT_COEFF_550_LKT(97,12)=6.186600 !rg=22.6687 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(97,13,1:6)=(/ 5.138800,5.149700,5.164400,5.184200,5.217700,5.282700 /)
XPIZA_LKT(97,13,1:6)=(/ 0.570486,0.590651,0.616090,0.772336,0.826598,0.867982 /)
XCGA_LKT(97,13,1:6)=(/ 0.937747,0.928200,0.918683,0.879060,0.877203,0.868427 /)
XEXT_COEFF_550_LKT(97,13)=5.161500 !rg=22.6687 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(97,14,1:6)=(/ 0.000000,4.296000,4.306400,4.323100,4.347100,4.398100 /)
XPIZA_LKT(97,14,1:6)=(/ 0.900,0.583345,0.606008,0.755106,0.812336,0.859505 /)
XCGA_LKT(97,14,1:6)=(/ 0.900,0.931560,0.922290,0.882107,0.879190,0.869893 /)
XEXT_COEFF_550_LKT(97,14)=4.305800 !rg=22.6687 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(97,15,1:6)=(/ 0.000000,3.580000,3.587200,3.598000,3.614200,3.648200 /)
XPIZA_LKT(97,15,1:6)=(/ 0.900,0.576975,0.596876,0.737067,0.796626,0.849258 /)
XCGA_LKT(97,15,1:6)=(/ 0.900,0.934770,0.925930,0.885177,0.881403,0.871627 /)
XEXT_COEFF_550_LKT(97,15)=3.587000 !rg=22.6687 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(97,16,1:6)=(/ 0.000000,0.000000,2.981500,2.989200,3.006400,3.033300 /)
XPIZA_LKT(97,16,1:6)=(/ 0.900,0.900,0.588778,0.719001,0.780270,0.837948 /)
XCGA_LKT(97,16,1:6)=(/ 0.900,0.900,0.929483,0.888717,0.884233,0.874307 /)
XEXT_COEFF_550_LKT(97,16)=2.981000 !rg=22.6687 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(97,17,1:6)=(/ 0.000000,0.000000,2.480800,2.488000,2.497300,2.517500 /)
XPIZA_LKT(97,17,1:6)=(/ 0.900,0.900,0.581614,0.701568,0.762665,0.824847 /)
XCGA_LKT(97,17,1:6)=(/ 0.900,0.900,0.932807,0.892437,0.886807,0.876133 /)
XEXT_COEFF_550_LKT(97,17)=2.480500 !rg=22.6687 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(97,18,1:6)=(/ 0.000000,0.000000,2.064600,2.068800,2.076800,2.090000 /)
XPIZA_LKT(97,18,1:6)=(/ 0.900,0.900,0.575470,0.684479,0.745103,0.810367 /)
XCGA_LKT(97,18,1:6)=(/ 0.900,0.900,0.935897,0.896257,0.889963,0.878673 /)
XEXT_COEFF_550_LKT(97,18)=2.064300 !rg=22.6687 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(97,19,1:6)=(/ 0.000000,0.000000,0.000000,1.718900,1.725500,1.736800 /)
XPIZA_LKT(97,19,1:6)=(/ 0.900,0.900,0.900,0.668227,0.727157,0.795033 /)
XCGA_LKT(97,19,1:6)=(/ 0.900,0.900,0.900,0.900440,0.893397,0.881387 /)
XEXT_COEFF_550_LKT(97,19)=0.000000 !rg=22.6687 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(97,20,1:6)=(/ 0.000000,0.000000,0.000000,1.428900,1.433000,1.441200 /)
XPIZA_LKT(97,20,1:6)=(/ 0.900,0.900,0.900,0.653004,0.709204,0.778351 /)
XCGA_LKT(97,20,1:6)=(/ 0.900,0.900,0.900,0.904690,0.896753,0.883963 /)
XEXT_COEFF_550_LKT(97,20)=0.000000 !rg=22.6687 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(98,1,1:6)=(/ 24.558000,24.732000,24.876000,24.877000,25.296000,26.332000 /)
XPIZA_LKT(98,1,1:6)=(/ 0.642210,0.689813,0.730140,0.880403,0.886410,0.912463 /)
XCGA_LKT(98,1,1:6)=(/ 0.908480,0.898147,0.891077,0.858863,0.863940,0.844343 /)
XEXT_COEFF_550_LKT(98,1)=24.732000 !rg=24.5583 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(98,2,1:6)=(/ 23.552000,23.686000,23.850000,23.979000,24.301000,24.803000 /)
XPIZA_LKT(98,2,1:6)=(/ 0.639673,0.686505,0.727052,0.878564,0.883944,0.910771 /)
XCGA_LKT(98,2,1:6)=(/ 0.909720,0.899267,0.891763,0.860230,0.862403,0.846020 /)
XEXT_COEFF_550_LKT(98,2)=23.880000 !rg=24.5583 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(98,3,1:6)=(/ 21.843000,21.938000,22.018000,22.249000,22.479000,23.208000 /)
XPIZA_LKT(98,3,1:6)=(/ 0.635003,0.680368,0.720287,0.875655,0.883483,0.903506 /)
XCGA_LKT(98,3,1:6)=(/ 0.911320,0.900657,0.892333,0.861270,0.862977,0.847943 /)
XEXT_COEFF_550_LKT(98,3)=21.998000 !rg=24.5583 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(98,4,1:6)=(/ 19.738000,19.831000,19.904000,20.102000,20.387000,20.899000 /)
XPIZA_LKT(98,4,1:6)=(/ 0.628596,0.672743,0.712365,0.870996,0.883121,0.898323 /)
XCGA_LKT(98,4,1:6)=(/ 0.913210,0.902343,0.893947,0.862180,0.863533,0.850223 /)
XEXT_COEFF_550_LKT(98,4)=19.889000 !rg=24.5583 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(98,5,1:6)=(/ 17.489000,17.554000,17.650000,17.788000,17.986000,18.387000 /)
XPIZA_LKT(98,5,1:6)=(/ 0.621440,0.663226,0.703442,0.864884,0.882689,0.895324 /)
XCGA_LKT(98,5,1:6)=(/ 0.915457,0.904400,0.895777,0.863567,0.865020,0.852383 /)
XEXT_COEFF_550_LKT(98,5)=17.677000 !rg=24.5583 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(98,6,1:6)=(/ 15.282000,15.342000,15.405000,15.523000,15.658000,15.959000 /)
XPIZA_LKT(98,6,1:6)=(/ 0.614145,0.653641,0.692709,0.857619,0.880314,0.893308 /)
XCGA_LKT(98,6,1:6)=(/ 0.918057,0.906900,0.897937,0.864503,0.866527,0.855663 /)
XEXT_COEFF_550_LKT(98,6)=15.394000 !rg=24.5583 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(98,7,1:6)=(/ 13.174000,13.220000,13.268000,13.367000,13.512000,13.792000 /)
XPIZA_LKT(98,7,1:6)=(/ 0.606324,0.643183,0.681134,0.848486,0.876042,0.891416 /)
XCGA_LKT(98,7,1:6)=(/ 0.920950,0.909687,0.900380,0.866653,0.866450,0.857030 /)
XEXT_COEFF_550_LKT(98,7)=13.261000 !rg=24.5583 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(98,8,1:6)=(/ 11.266000,11.305000,11.344000,11.437000,11.494000,11.713000 /)
XPIZA_LKT(98,8,1:6)=(/ 0.598612,0.632599,0.669110,0.838247,0.870561,0.888849 /)
XCGA_LKT(98,8,1:6)=(/ 0.923993,0.912800,0.903183,0.868613,0.867863,0.858523 /)
XEXT_COEFF_550_LKT(98,8)=11.343000 !rg=24.5583 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(98,9,1:6)=(/ 9.559300,9.586700,9.621700,9.679300,9.762100,9.924900 /)
XPIZA_LKT(98,9,1:6)=(/ 0.591346,0.622176,0.656664,0.825902,0.864142,0.886265 /)
XCGA_LKT(98,9,1:6)=(/ 0.927103,0.916020,0.906370,0.870377,0.870183,0.860633 /)
XEXT_COEFF_550_LKT(98,9)=9.622300 !rg=24.5583 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(98,10,1:6)=(/ 8.079300,8.101100,8.132700,8.159400,8.237200,8.358200 /)
XPIZA_LKT(98,10,1:6)=(/ 0.584576,0.612360,0.644794,0.812422,0.855630,0.882927 /)
XCGA_LKT(98,10,1:6)=(/ 0.930257,0.919450,0.909680,0.872077,0.872110,0.863617 /)
XEXT_COEFF_550_LKT(98,10)=8.124200 !rg=24.5583 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(98,11,1:6)=(/ 6.786800,6.804100,6.824700,6.862300,6.895900,7.006500 /)
XPIZA_LKT(98,11,1:6)=(/ 0.578334,0.603012,0.632666,0.797871,0.845431,0.877670 /)
XCGA_LKT(98,11,1:6)=(/ 0.933373,0.922987,0.913113,0.875033,0.874087,0.864887 /)
XEXT_COEFF_550_LKT(98,11)=6.823700 !rg=24.5583 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(98,12,1:6)=(/ 5.680900,5.693700,5.710600,5.727900,5.769600,5.848100 /)
XPIZA_LKT(98,12,1:6)=(/ 0.572792,0.594447,0.621403,0.781500,0.833698,0.871430 /)
XCGA_LKT(98,12,1:6)=(/ 0.936350,0.926473,0.916777,0.877453,0.875977,0.866963 /)
XEXT_COEFF_550_LKT(98,12)=5.706600 !rg=24.5583 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(98,13,1:6)=(/ 0.000000,4.751500,4.763100,4.780600,4.807200,4.866000 /)
XPIZA_LKT(98,13,1:6)=(/ 0.900,0.586671,0.610768,0.764318,0.820330,0.864402 /)
XCGA_LKT(98,13,1:6)=(/ 0.900,0.929913,0.920437,0.880420,0.878030,0.869240 /)
XEXT_COEFF_550_LKT(98,13)=4.762800 !rg=24.5583 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(98,14,1:6)=(/ 0.000000,3.964200,3.972800,3.986900,4.007600,4.050400 /)
XPIZA_LKT(98,14,1:6)=(/ 0.900,0.579892,0.601227,0.746754,0.805407,0.855074 /)
XCGA_LKT(98,14,1:6)=(/ 0.900,0.933190,0.924083,0.883457,0.880067,0.871097 /)
XEXT_COEFF_550_LKT(98,14)=3.973000 !rg=24.5583 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(98,15,1:6)=(/ 0.000000,3.303300,3.310300,3.319300,3.339200,3.371900 /)
XPIZA_LKT(98,15,1:6)=(/ 0.900,0.573956,0.592746,0.728811,0.789585,0.844694 /)
XCGA_LKT(98,15,1:6)=(/ 0.900,0.936267,0.927670,0.886830,0.882783,0.873440 /)
XEXT_COEFF_550_LKT(98,15)=3.309600 !rg=24.5583 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(98,16,1:6)=(/ 0.000000,0.000000,2.750800,2.758400,2.768600,2.791100 /)
XPIZA_LKT(98,16,1:6)=(/ 0.900,0.900,0.585039,0.710930,0.772226,0.832019 /)
XCGA_LKT(98,16,1:6)=(/ 0.900,0.900,0.931117,0.890433,0.885263,0.874877 /)
XEXT_COEFF_550_LKT(98,16)=2.750300 !rg=24.5583 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(98,17,1:6)=(/ 0.000000,0.000000,2.289300,2.294400,2.304300,2.321500 /)
XPIZA_LKT(98,17,1:6)=(/ 0.900,0.900,0.578430,0.693490,0.754860,0.818636 /)
XCGA_LKT(98,17,1:6)=(/ 0.900,0.900,0.934337,0.894183,0.888247,0.877313 /)
XEXT_COEFF_550_LKT(98,17)=2.289000 !rg=24.5583 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(98,18,1:6)=(/ 0.000000,0.000000,0.000000,1.909000,1.917300,1.930800 /)
XPIZA_LKT(98,18,1:6)=(/ 0.900,0.900,0.900,0.676848,0.737132,0.803860 /)
XCGA_LKT(98,18,1:6)=(/ 0.900,0.900,0.900,0.898183,0.891590,0.880043 /)
XEXT_COEFF_550_LKT(98,18)=0.000000 !rg=24.5583 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(98,19,1:6)=(/ 0.000000,0.000000,0.000000,1.586200,1.590500,1.600300 /)
XPIZA_LKT(98,19,1:6)=(/ 0.900,0.900,0.900,0.661093,0.718803,0.787535 /)
XCGA_LKT(98,19,1:6)=(/ 0.900,0.900,0.900,0.902370,0.894723,0.882267 /)
XEXT_COEFF_550_LKT(98,19)=0.000000 !rg=24.5583 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(98,20,1:6)=(/ 0.000000,0.000000,0.000000,1.318300,1.322300,1.329900 /)
XPIZA_LKT(98,20,1:6)=(/ 0.900,0.900,0.900,0.646262,0.701320,0.770825 /)
XCGA_LKT(98,20,1:6)=(/ 0.900,0.900,0.900,0.906643,0.898437,0.885433 /)
XEXT_COEFF_550_LKT(98,20)=0.000000 !rg=24.5583 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(99,1,1:6)=(/ 22.735000,22.875000,22.971000,23.140000,23.567000,24.274000 /)
XPIZA_LKT(99,1,1:6)=(/ 0.635815,0.681996,0.723315,0.877620,0.889178,0.906424 /)
XCGA_LKT(99,1,1:6)=(/ 0.910807,0.900413,0.891407,0.860570,0.866137,0.846480 /)
XEXT_COEFF_550_LKT(99,1)=22.961000 !rg=26.6054 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(99,2,1:6)=(/ 21.753000,21.863000,21.925000,22.110000,22.417000,22.765000 /)
XPIZA_LKT(99,2,1:6)=(/ 0.632818,0.678912,0.719304,0.875994,0.887731,0.903601 /)
XCGA_LKT(99,2,1:6)=(/ 0.911760,0.901160,0.892477,0.860800,0.865370,0.850573 /)
XEXT_COEFF_550_LKT(99,2)=21.910000 !rg=26.6054 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(99,3,1:6)=(/ 20.147000,20.243000,20.347000,20.457000,20.717000,21.213000 /)
XPIZA_LKT(99,3,1:6)=(/ 0.627974,0.672550,0.713814,0.871880,0.885796,0.897915 /)
XCGA_LKT(99,3,1:6)=(/ 0.913190,0.902520,0.894057,0.861903,0.865013,0.851070 /)
XEXT_COEFF_550_LKT(99,3)=20.340000 !rg=26.6054 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(99,4,1:6)=(/ 18.209000,18.285000,18.390000,18.484000,18.704000,19.128000 /)
XPIZA_LKT(99,4,1:6)=(/ 0.621975,0.664876,0.706049,0.867083,0.884516,0.894222 /)
XCGA_LKT(99,4,1:6)=(/ 0.915107,0.904087,0.895467,0.863017,0.865420,0.853417 /)
XEXT_COEFF_550_LKT(99,4)=18.390000 !rg=26.6054 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(99,5,1:6)=(/ 16.141000,16.197000,16.261000,16.395000,16.550000,16.991000 /)
XPIZA_LKT(99,5,1:6)=(/ 0.615363,0.656049,0.695949,0.860514,0.881747,0.891120 /)
XCGA_LKT(99,5,1:6)=(/ 0.917457,0.906177,0.897180,0.864263,0.865190,0.854697 /)
XEXT_COEFF_550_LKT(99,5)=16.266000 !rg=26.6054 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(99,6,1:6)=(/ 14.099000,14.147000,14.220000,14.268000,14.509000,14.816000 /)
XPIZA_LKT(99,6,1:6)=(/ 0.608292,0.646552,0.685662,0.852454,0.878672,0.891009 /)
XCGA_LKT(99,6,1:6)=(/ 0.919997,0.908670,0.899623,0.865560,0.866440,0.856830 /)
XEXT_COEFF_550_LKT(99,6)=14.197000 !rg=26.6054 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(99,7,1:6)=(/ 12.154000,12.195000,12.248000,12.296000,12.459000,12.682000 /)
XPIZA_LKT(99,7,1:6)=(/ 0.600998,0.636388,0.674233,0.842908,0.874546,0.889325 /)
XCGA_LKT(99,7,1:6)=(/ 0.922857,0.911580,0.902093,0.867190,0.868820,0.858803 /)
XEXT_COEFF_550_LKT(99,7)=12.244000 !rg=26.6054 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(99,8,1:6)=(/ 10.395000,10.426000,10.470000,10.520000,10.625000,10.811000 /)
XPIZA_LKT(99,8,1:6)=(/ 0.593866,0.626170,0.662196,0.831986,0.868213,0.886909 /)
XCGA_LKT(99,8,1:6)=(/ 0.925857,0.914647,0.905033,0.869433,0.870160,0.861360 /)
XEXT_COEFF_550_LKT(99,8)=10.462000 !rg=26.6054 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(99,9,1:6)=(/ 8.820700,8.844900,8.873500,8.917100,8.981100,9.144400 /)
XPIZA_LKT(99,9,1:6)=(/ 0.587053,0.616340,0.649991,0.819436,0.860150,0.883721 /)
XCGA_LKT(99,9,1:6)=(/ 0.928970,0.917957,0.908033,0.871273,0.870433,0.861587 /)
XEXT_COEFF_550_LKT(99,9)=8.873700 !rg=26.6054 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(99,10,1:6)=(/ 7.454900,7.473300,7.496700,7.534100,7.594700,7.706900 /)
XPIZA_LKT(99,10,1:6)=(/ 0.580825,0.606926,0.638203,0.805637,0.851286,0.880503 /)
XCGA_LKT(99,10,1:6)=(/ 0.931977,0.921327,0.911400,0.873647,0.873027,0.864370 /)
XEXT_COEFF_550_LKT(99,10)=7.499700 !rg=26.6054 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(99,11,1:6)=(/ 6.262200,6.277500,6.295700,6.313100,6.356300,6.445500 /)
XPIZA_LKT(99,11,1:6)=(/ 0.575043,0.598224,0.626655,0.789984,0.840177,0.875007 /)
XCGA_LKT(99,11,1:6)=(/ 0.934993,0.924803,0.914977,0.876077,0.874643,0.865660 /)
XEXT_COEFF_550_LKT(99,11)=6.293800 !rg=26.6054 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(99,12,1:6)=(/ 5.242000,5.252800,5.267000,5.289300,5.322100,5.385300 /)
XPIZA_LKT(99,12,1:6)=(/ 0.569974,0.590168,0.615887,0.773668,0.827859,0.868537 /)
XCGA_LKT(99,12,1:6)=(/ 0.937830,0.928210,0.918547,0.878993,0.876900,0.867987 /)
XEXT_COEFF_550_LKT(99,12)=5.266900 !rg=26.6054 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(99,13,1:6)=(/ 0.000000,4.384500,4.394600,4.409500,4.434000,4.483500 /)
XPIZA_LKT(99,13,1:6)=(/ 0.900,0.582982,0.605781,0.756146,0.813656,0.860055 /)
XCGA_LKT(99,13,1:6)=(/ 0.900,0.931583,0.922243,0.881777,0.878940,0.869653 /)
XEXT_COEFF_550_LKT(99,13)=4.394700 !rg=26.6054 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(99,14,1:6)=(/ 0.000000,3.658000,3.665900,3.675600,3.699400,3.739000 /)
XPIZA_LKT(99,14,1:6)=(/ 0.900,0.576671,0.596801,0.738364,0.798611,0.850886 /)
XCGA_LKT(99,14,1:6)=(/ 0.900,0.934753,0.925853,0.884953,0.881490,0.872163 /)
XEXT_COEFF_550_LKT(99,14)=3.665600 !rg=26.6054 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(99,15,1:6)=(/ 0.000000,3.048200,3.054100,3.063100,3.076000,3.102300 /)
XPIZA_LKT(99,15,1:6)=(/ 0.900,0.571179,0.588727,0.720684,0.781904,0.839367 /)
XCGA_LKT(99,15,1:6)=(/ 0.900,0.937700,0.929370,0.888443,0.883767,0.874083 /)
XEXT_COEFF_550_LKT(99,15)=3.053300 !rg=26.6054 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(99,16,1:6)=(/ 0.000000,0.000000,2.538200,2.543900,2.555200,2.575500 /)
XPIZA_LKT(99,16,1:6)=(/ 0.900,0.900,0.581568,0.702691,0.764530,0.826387 /)
XCGA_LKT(99,16,1:6)=(/ 0.900,0.900,0.932707,0.892093,0.886723,0.876287 /)
XEXT_COEFF_550_LKT(99,16)=2.537700 !rg=26.6054 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(99,17,1:6)=(/ 0.000000,0.000000,2.112300,2.116900,2.126400,2.142600 /)
XPIZA_LKT(99,17,1:6)=(/ 0.900,0.900,0.575364,0.685629,0.746921,0.812285 /)
XCGA_LKT(99,17,1:6)=(/ 0.900,0.900,0.935833,0.896007,0.889900,0.878773 /)
XEXT_COEFF_550_LKT(99,17)=2.112200 !rg=26.6054 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(99,18,1:6)=(/ 0.000000,0.000000,0.000000,1.761600,1.767300,1.778800 /)
XPIZA_LKT(99,18,1:6)=(/ 0.900,0.900,0.900,0.669452,0.728816,0.796838 /)
XCGA_LKT(99,18,1:6)=(/ 0.900,0.900,0.900,0.900083,0.892987,0.881007 /)
XEXT_COEFF_550_LKT(99,18)=0.000000 !rg=26.6054 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(99,19,1:6)=(/ 0.000000,0.000000,0.000000,1.463500,1.468000,1.476900 /)
XPIZA_LKT(99,19,1:6)=(/ 0.900,0.900,0.900,0.654006,0.710905,0.780193 /)
XCGA_LKT(99,19,1:6)=(/ 0.900,0.900,0.900,0.904310,0.896457,0.883687 /)
XEXT_COEFF_550_LKT(99,19)=0.000000 !rg=26.6054 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(99,20,1:6)=(/ 0.000000,0.000000,0.000000,0.000000,1.220100,1.226900 /)
XPIZA_LKT(99,20,1:6)=(/ 0.900,0.900,0.900,0.900,0.693473,0.763204 /)
XCGA_LKT(99,20,1:6)=(/ 0.900,0.900,0.900,0.900,0.900167,0.886890 /)
XEXT_COEFF_550_LKT(99,20)=0.000000 !rg=26.6054 sigma=2.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(100,1,1:6)=(/ 20.982000,21.090000,21.183000,21.503000,21.559000,22.000000 /)
XPIZA_LKT(100,1,1:6)=(/ 0.628567,0.674456,0.715980,0.875687,0.890660,0.898705 /)
XCGA_LKT(100,1,1:6)=(/ 0.913093,0.902423,0.893757,0.864203,0.861800,0.848947 /)
XEXT_COEFF_550_LKT(100,1)=21.258000 !rg=28.8231 sigma=1.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(100,2,1:6)=(/ 20.075000,20.151000,20.254000,20.353000,20.653000,21.275000 /)
XPIZA_LKT(100,2,1:6)=(/ 0.625903,0.670913,0.712071,0.871715,0.888637,0.895532 /)
XCGA_LKT(100,2,1:6)=(/ 0.913840,0.902910,0.894330,0.862043,0.862753,0.852503 /)
XEXT_COEFF_550_LKT(100,2)=20.224000 !rg=28.8231 sigma=1.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(100,3,1:6)=(/ 18.577000,18.660000,18.731000,18.953000,19.163000,19.591000 /)
XPIZA_LKT(100,3,1:6)=(/ 0.620957,0.664710,0.705692,0.868986,0.886689,0.890430 /)
XCGA_LKT(100,3,1:6)=(/ 0.915167,0.904067,0.895297,0.863143,0.865127,0.854350 /)
XEXT_COEFF_550_LKT(100,3)=18.750000 !rg=28.8231 sigma=1.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(100,4,1:6)=(/ 16.797000,16.868000,16.925000,17.132000,17.295000,17.680000 /)
XPIZA_LKT(100,4,1:6)=(/ 0.615520,0.657266,0.697855,0.863263,0.884524,0.888880 /)
XCGA_LKT(100,4,1:6)=(/ 0.917080,0.905847,0.896800,0.864633,0.866537,0.856007 /)
XEXT_COEFF_550_LKT(100,4)=16.931000 !rg=28.8231 sigma=1.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(100,5,1:6)=(/ 14.886000,14.947000,15.011000,15.088000,15.264000,15.546000 /)
XPIZA_LKT(100,5,1:6)=(/ 0.609268,0.648801,0.688935,0.855729,0.881061,0.889654 /)
XCGA_LKT(100,5,1:6)=(/ 0.919370,0.908093,0.898893,0.865333,0.867023,0.857547 /)
XEXT_COEFF_550_LKT(100,5)=15.005000 !rg=28.8231 sigma=1.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(100,6,1:6)=(/ 13.007000,13.049000,13.099000,13.209000,13.316000,13.588000 /)
XPIZA_LKT(100,6,1:6)=(/ 0.602794,0.639501,0.678241,0.848071,0.876330,0.888564 /)
XCGA_LKT(100,6,1:6)=(/ 0.921933,0.910470,0.901033,0.867110,0.866917,0.858390 /)
XEXT_COEFF_550_LKT(100,6)=13.102000 !rg=28.8231 sigma=1.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(100,7,1:6)=(/ 11.212000,11.249000,11.292000,11.377000,11.455000,11.695000 /)
XPIZA_LKT(100,7,1:6)=(/ 0.595902,0.629824,0.666843,0.837764,0.871358,0.887061 /)
XCGA_LKT(100,7,1:6)=(/ 0.924787,0.913460,0.903757,0.868803,0.869207,0.860453 /)
XEXT_COEFF_550_LKT(100,7)=11.291000 !rg=28.8231 sigma=1.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(100,8,1:6)=(/ 9.590800,9.616800,9.652000,9.705400,9.785100,9.955100 /)
XPIZA_LKT(100,8,1:6)=(/ 0.589363,0.620109,0.655028,0.826054,0.864980,0.885791 /)
XCGA_LKT(100,8,1:6)=(/ 0.927687,0.916500,0.906670,0.870750,0.871073,0.862310 /)
XEXT_COEFF_550_LKT(100,8)=9.652000 !rg=28.8231 sigma=1.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(100,9,1:6)=(/ 8.138300,8.160600,8.188200,8.221600,8.300300,8.415300 /)
XPIZA_LKT(100,9,1:6)=(/ 0.583082,0.610815,0.643422,0.812517,0.856505,0.882574 /)
XCGA_LKT(100,9,1:6)=(/ 0.930700,0.919807,0.909920,0.872517,0.872507,0.863810 /)
XEXT_COEFF_550_LKT(100,9)=8.186600 !rg=28.8231 sigma=1.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(100,10,1:6)=(/ 6.879000,6.896200,6.915100,6.942400,6.996000,7.083500 /)
XPIZA_LKT(100,10,1:6)=(/ 0.577330,0.602031,0.631807,0.798260,0.846362,0.877976 /)
XCGA_LKT(100,10,1:6)=(/ 0.933660,0.923143,0.913227,0.874723,0.873727,0.865413 /)
XEXT_COEFF_550_LKT(100,10)=6.915800 !rg=28.8231 sigma=1.95 wvl=0.55
 
XEXT_COEFF_WVL_LKT(100,11,1:6)=(/ 5.778600,5.791300,5.807600,5.829900,5.874100,5.956500 /)
XPIZA_LKT(100,11,1:6)=(/ 0.572042,0.593693,0.620903,0.782410,0.834804,0.872015 /)
XCGA_LKT(100,11,1:6)=(/ 0.936537,0.926570,0.916773,0.877370,0.875550,0.866947 /)
XEXT_COEFF_550_LKT(100,11)=5.807600 !rg=28.8231 sigma=2.05 wvl=0.55
 
XEXT_COEFF_WVL_LKT(100,12,1:6)=(/ 4.837100,4.847100,4.858500,4.877000,4.904700,4.961700 /)
XPIZA_LKT(100,12,1:6)=(/ 0.567323,0.586216,0.610544,0.765554,0.821489,0.864977 /)
XCGA_LKT(100,12,1:6)=(/ 0.939293,0.929940,0.920347,0.880310,0.877807,0.869093 /)
XEXT_COEFF_550_LKT(100,12)=4.858400 !rg=28.8231 sigma=2.15 wvl=0.55
 
XEXT_COEFF_WVL_LKT(100,13,1:6)=(/ 0.000000,4.045600,4.055000,4.065900,4.093800,4.138200 /)
XPIZA_LKT(100,13,1:6)=(/ 0.900,0.579464,0.601067,0.747849,0.807032,0.856428 /)
XCGA_LKT(100,13,1:6)=(/ 0.900,0.933227,0.924053,0.883227,0.880237,0.871270 /)
XEXT_COEFF_550_LKT(100,13)=4.054500 !rg=28.8231 sigma=2.25 wvl=0.55
 
XEXT_COEFF_WVL_LKT(100,14,1:6)=(/ 0.000000,3.375400,3.382200,3.393800,3.403600,3.436200 /)
XPIZA_LKT(100,14,1:6)=(/ 0.900,0.573672,0.592563,0.730328,0.790882,0.845374 /)
XCGA_LKT(100,14,1:6)=(/ 0.900,0.936253,0.927587,0.886657,0.882310,0.872583 /)
XEXT_COEFF_550_LKT(100,14)=3.381600 !rg=28.8231 sigma=2.35 wvl=0.55
 
XEXT_COEFF_WVL_LKT(100,15,1:6)=(/ 0.000000,0.000000,2.818000,2.825200,2.837300,2.861800 /)
XPIZA_LKT(100,15,1:6)=(/ 0.900,0.900,0.585054,0.712345,0.774175,0.833906 /)
XCGA_LKT(100,15,1:6)=(/ 0.900,0.900,0.930977,0.890030,0.885067,0.875200 /)
XEXT_COEFF_550_LKT(100,15)=2.817100 !rg=28.8231 sigma=2.45 wvl=0.55
 
XEXT_COEFF_WVL_LKT(100,16,1:6)=(/ 0.000000,0.000000,2.342000,2.347200,2.358900,2.377500 /)
XPIZA_LKT(100,16,1:6)=(/ 0.900,0.900,0.578301,0.694652,0.756646,0.820422 /)
XCGA_LKT(100,16,1:6)=(/ 0.900,0.900,0.934253,0.893887,0.888227,0.877450 /)
XEXT_COEFF_550_LKT(100,16)=2.341700 !rg=28.8231 sigma=2.55 wvl=0.55
 
XEXT_COEFF_WVL_LKT(100,17,1:6)=(/ 0.000000,0.000000,0.000000,1.953300,1.960200,1.972500 /)
XPIZA_LKT(100,17,1:6)=(/ 0.900,0.900,0.900,0.677978,0.738648,0.805453 /)
XCGA_LKT(100,17,1:6)=(/ 0.900,0.900,0.900,0.897820,0.891203,0.879753 /)
XEXT_COEFF_550_LKT(100,17)=0.000000 !rg=28.8231 sigma=2.65 wvl=0.55
 
XEXT_COEFF_WVL_LKT(100,18,1:6)=(/ 0.000000,0.000000,0.000000,1.625400,1.630400,1.640400 /)
XPIZA_LKT(100,18,1:6)=(/ 0.900,0.900,0.900,0.662174,0.720632,0.789541 /)
XCGA_LKT(100,18,1:6)=(/ 0.900,0.900,0.900,0.901963,0.894483,0.882260 /)
XEXT_COEFF_550_LKT(100,18)=0.000000 !rg=28.8231 sigma=2.75 wvl=0.55
 
XEXT_COEFF_WVL_LKT(100,19,1:6)=(/ 0.000000,0.000000,0.000000,1.350400,1.354700,1.362800 /)
XPIZA_LKT(100,19,1:6)=(/ 0.900,0.900,0.900,0.647238,0.702931,0.772925 /)
XCGA_LKT(100,19,1:6)=(/ 0.900,0.900,0.900,0.906290,0.898083,0.885340 /)
XEXT_COEFF_550_LKT(100,19)=0.000000 !rg=28.8231 sigma=2.85 wvl=0.55
 
XEXT_COEFF_WVL_LKT(100,20,1:6)=(/ 0.000000,0.000000,0.000000,0.000000,1.125500,1.130900 /)
XPIZA_LKT(100,20,1:6)=(/ 0.900,0.900,0.900,0.900,0.685665,0.755238 /)
XCGA_LKT(100,20,1:6)=(/ 0.900,0.900,0.900,0.900,0.901823,0.888070 /)
XEXT_COEFF_550_LKT(100,20)=0.000000 !rg=28.8231 sigma=2.95 wvl=0.55
 
IF (LHOOK) CALL DR_HOOK('MODE_DUSTOPT:DUST_OPT_LKT_SET10',1,ZHOOK_HANDLE)
END SUBROUTINE DUST_OPT_LKT_SET10


END MODULE MODE_DUSTOPT
