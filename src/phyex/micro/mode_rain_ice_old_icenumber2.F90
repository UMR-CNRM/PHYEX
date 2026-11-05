MODULE MODE_RAIN_ICE_OLD_ICENUMBER2

  IMPLICIT NONE

  CONTAINS

  PURE FUNCTION ICENUMBER2(Q_ICE, T3D) RESULT(ICENUMBER)

    IMPLICIT NONE

    REAL, INTENT(IN) :: Q_ICE
    REAL, INTENT(IN) :: T3D

    REAL, PARAMETER:: ICE_DENSITY = 890.0
    REAL, PARAMETER:: PI = 4.0*ATAN(1.)
    INTEGER IDX_REI
    REAL CORR, REICE, DEICE
    DOUBLE PRECISION LAMBDA

    REAL :: ICENUMBER

!+---+-----------------------------------------------------------------+
!..Table of lookup values of radiative effective radius of ice crystals
!.. as a function of Temperature from -94C to 0C.  Taken from WRF RRTMG
!.. radiation code where it is attributed to Jon Egill Kristjansson
!.. and coauthors.
!+---+-----------------------------------------------------------------+

    REAL, DIMENSION(95), PARAMETER :: RETAB = (/5.92779, 6.26422, 6.61973, 6.99539, 7.39234,          &
                                                7.81177, 8.25496, 8.72323, 9.21800, 9.74075, 10.2930, &
                                                10.8765, 11.4929, 12.1440, 12.8317, 13.5581, 14.2319, &
                                                15.0351, 15.8799, 16.7674, 17.6986, 18.6744, 19.6955, &
                                                20.7623, 21.8757, 23.0364, 24.2452, 25.5034, 26.8125, &
                                                27.7895, 28.6450, 29.4167, 30.1088, 30.7306, 31.2943, &
                                                31.8151, 32.3077, 32.7870, 33.2657, 33.7540, 34.2601, &
                                                34.7892, 35.3442, 35.9255, 36.5316, 37.1602, 37.8078, &
                                                38.4720, 39.1508, 39.8442, 40.5552, 41.2912, 42.0635, &
                                                42.8876, 43.7863, 44.7853, 45.9170, 47.2165, 48.7221, &
                                                50.4710, 52.4980, 54.8315, 57.4898, 60.4785, 63.7898, &
                                                65.5604, 71.2885, 75.4113, 79.7368, 84.2351, 88.8833, &
                                                93.6658, 98.5739, 103.603, 108.752, 114.025, 119.424, &
                                                124.954, 130.630, 136.457, 142.446, 148.608, 154.956, &
                                                161.503, 168.262, 175.248, 182.473, 189.952, 197.699, &
                                                205.728, 214.055, 222.694, 231.661, 240.971, 250.639/)

!+---+-----------------------------------------------------------------+
!..From the model 3D temperature field, subtract 179K for which
!.. index value of retab as a start.  Value of corr is for
!.. interpolating between neighboring values in the table.
!+---+-----------------------------------------------------------------+

    IDX_REI = INT(T3D-179.)
    IDX_REI = MIN(MAX(IDX_REI,1),95)
    CORR = T3D - INT(T3D)
    REICE = RETAB(IDX_REI)*(1.-CORR) + RETAB(MIN(95,IDX_REI+1))*CORR
    DEICE = 2.*REICE * 1.E-6

!+---+-----------------------------------------------------------------+
!..Now we have the final radiative effective size of ice (as function
!.. of temperature only).  This size represents 3rd moment divided by
!.. second moment of the ice size distribution, so we can compute a
!.. number concentration from the mean size and mass mixing ratio.
!.. The mean (radiative effective) diameter is 3./Slope for an inverse
!.. exponential size distribution.  So, starting with slope, work
!.. backwords to get number concentration.
!+---+-----------------------------------------------------------------+

    LAMBDA = 3.0 / DEICE
    ICENUMBER = Q_ICE * LAMBDA*LAMBDA*LAMBDA / (PI*ICE_DENSITY)

!+---+-----------------------------------------------------------------+
!..Example1: Common ice size coming from Thompson scheme is about 30 microns.
!.. An example ice mixing ratio could be 0.001 g/kg for a temperature of -50C.
!.. Remember to convert both into MKS units.  This gives N_ice=357652 per kg.
!..Example2: Lower in atmosphere at T=-10C matching ~162 microns in retab,
!.. and assuming we have 0.1 g/kg mixing ratio, then N_ice=28122 per kg,
!.. which is 28 crystals per liter of air if the air density is 1.0.
!+---+-----------------------------------------------------------------+

  END FUNCTION ICENUMBER2

END MODULE MODE_RAIN_ICE_OLD_ICENUMBER2
