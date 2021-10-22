!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source: /home/tulet/MESONH/prepsource/test/isocom.f,v $ $Revision: 1.1 $ $Date: 2005/06/02 10:27:31 $
!-----------------------------------------------------------------
!=======================================================================
!
! *** ISORROPIA CODE
! *** SUBROUTINE ISOROPIA
! *** THIS SUBROUTINE IS THE MASTER ROUTINE FOR THE ISORROPIA
!     THERMODYNAMIC EQUILIBRIUM AEROSOL MODEL (VERSION 1.1 and above)
!
! ======================== ARGUMENTS / USAGE ===========================
!
!  INPUT:
!  1. [WI]
!     DOUBLE PRECISION array of length [5].
!     Concentrations, expressed in moles/m3. Depending on the type of
!     problem solved (specified in CNTRL(1)), WI contains either
!     GAS+AEROSOL or AEROSOL only concentratios.
!     WI(1) - sodium
!     WI(2) - sulfate
!     WI(3) - ammonium
!     WI(4) - nitrate
!     WI(5) - chloride
!
!  2. [RHI]
!     DOUBLE PRECISION variable.
!     Ambient relative humidity expressed on a (0,1) scale.
!
!  3. [TEMPI]
!     DOUBLE PRECISION variable.
!     Ambient temperature expressed in Kelvins.
!
!  4. [CNTRL]
!     DOUBLE PRECISION array of length [2].
!     Parameters that control the type of problem solved.
!
!     CNTRL(1): Defines the type of problem solved.
!     0 - Forward problem is solved. In this case, array WI contains
!         GAS and AEROSOL concentrations together.
!     1 - Reverse problem is solved. In this case, array WI contains
!         AEROSOL concentrations only.
!
!     CNTRL(2): Defines the state of the aerosol
!     0 - The aerosol can have both solid+liquid phases (deliquescent)
!     1 - The aerosol is in only liquid state (metastable aerosol)
!
!  OUTPUT:
!  1. [WT]
!     DOUBLE PRECISION array of length [5].
!     Total concentrations (GAS+AEROSOL) of species, expressed in moles/m3.
!     If the foreward probelm is solved (CNTRL(1)=0), array WT is
!     identical to array WI.
!     WT(1) - total sodium
!     WT(2) - total sulfate
!     WT(3) - total ammonium
!     WT(4) - total nitrate
!     WT(5) - total chloride
!
!  2. [GAS]
!     DOUBLE PRECISION array of length [03].
!     Gaseous species concentrations, expressed in moles/m3.
!     GAS(1) - NH3
!     GAS(2) - HNO3
!     GAS(3) - HCl
!
!  3. [AERLIQ]
!     DOUBLE PRECISION array of length [11].
!     Liquid aerosol species concentrations, expressed in moles/m3.
!     AERLIQ(01) - H+(aq)
!     AERLIQ(02) - Na+(aq)
!     AERLIQ(03) - NH4+(aq)
!     AERLIQ(04) - Cl-(aq)
!     AERLIQ(05) - SO4--(aq)
!     AERLIQ(06) - HSO4-(aq)
!     AERLIQ(07) - NO3-(aq)
!     AERLIQ(08) - H2O
!     AERLIQ(09) - NH3(aq) (undissociated)
!     AERLIQ(10) - HNCl(aq) (undissociated)
!     AERLIQ(11) - HNO3(aq) (undissociated)
!     AERLIQ(12) - OH-(aq)
!
!  4. [AERSLD]
!     DOUBLE PRECISION array of length [09].
!     Solid aerosol species concentrations, expressed in moles/m3.
!     AERSLD(01) - NaNO3(s)
!     AERSLD(02) - NH4NO3(s)
!     AERSLD(03) - NaCl(s)
!     AERSLD(04) - NH4Cl(s)
!     AERSLD(05) - Na2SO4(s)
!     AERSLD(06) - (NH4)2SO4(s)
!     AERSLD(07) - NaHSO4(s)
!     AERSLD(08) - NH4HSO4(s)
!     AERSLD(09) - (NH4)4H(SO4)2(s)
!
!  5. [SCASI]
!     CHARACTER*15 variable.
!     Returns the subcase which the input corresponds to.
!
!  6. [OTHER]
!     DOUBLE PRECISION array of length [6].
!     Returns solution information.
!
!     OTHER(1): Shows if aerosol water exists.
!     0 - Aerosol is WET
!     1 - Aerosol is DRY
!
!     OTHER(2): Aerosol Sulfate ratio, defined as (in moles/m3) :
!               (total ammonia + total Na) / (total sulfate)
!
!     OTHER(3): Sulfate ratio based on aerosol properties that defines
!               a sulfate poor system:
!               (aerosol ammonia + aerosol Na) / (aerosol sulfate)
!
!     OTHER(4): Aerosol sodium ratio, defined as (in moles/m3) :
!               (total Na) / (total sulfate)
!
!     OTHER(5): Ionic strength of the aqueous aerosol (if it exists).
!
!     OTHER(6): Total number of calls to the activity coefficient
!               calculation subroutine.
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!=======================================================================
!
      SUBROUTINE ISOROPIA (WI, RHI, TEMPI,  CNTRL,                      &
     &                     WT, GAS, AERLIQ, AERSLD, SCASI, OTHER)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
      IMPLICIT NONE
      INCLUDE 'isrpia.inc'
      REAL(KIND=8) :: WI,RHI,TEMPI,CNTRL,WT,GAS,AERLIQ,AERSLD,OTHER
      INTEGER :: NCTRL,NOTHER,I
      PARAMETER (NCTRL=2,NOTHER=6)
      CHARACTER SCASI*15
      DIMENSION WI(NCOMP), WT(NCOMP),   GAS(NGASAQ),  AERSLD(NSLDS),    &
     &          AERLIQ(NIONS+NGASAQ+2), CNTRL(NCTRL), OTHER(NOTHER)
!
! *** PROBLEM TYPE (0=FOREWARD, 1=REVERSE) ******************************
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('ISOROPIA',0,ZHOOK_HANDLE)
      IPROB   = NINT(CNTRL(1))
!
! *** AEROSOL STATE (0=SOLID+LIQUID, 1=METASTABLE) **********************
!
      METSTBL = NINT(CNTRL(2))
!
! *** SOLVE FOREWARD PROBLEM ********************************************
!
50    IF (IPROB.EQ.0) THEN
         IF (WI(1)+WI(2)+WI(3)+WI(4)+WI(5) .LE. TINY) THEN ! Everything=0
            CALL INIT1 (WI, RHI, TEMPI)
         ELSE IF (WI(1)+WI(4)+WI(5) .LE. TINY) THEN        ! Na,Cl,NO3=0
            CALL ISRP1F (WI, RHI, TEMPI)
         ELSE IF (WI(1)+WI(5) .LE. TINY) THEN              ! Na,Cl=0
            CALL ISRP2F (WI, RHI, TEMPI)
         ELSE
            CALL ISRP3F (WI, RHI, TEMPI)
         ENDIF
!
! *** SOLVE REVERSE PROBLEM *********************************************
!
      ELSE
         IF (WI(1)+WI(2)+WI(3)+WI(4)+WI(5) .LE. TINY) THEN ! Everything=0
            CALL INIT1 (WI, RHI, TEMPI)
         ELSE IF (WI(1)+WI(4)+WI(5) .LE. TINY) THEN        ! Na,Cl,NO3=0
            CALL ISRP1R (WI, RHI, TEMPI)
         ELSE IF (WI(1)+WI(5) .LE. TINY) THEN              ! Na,Cl=0
            CALL ISRP2R (WI, RHI, TEMPI)
         ELSE
            CALL ISRP3R (WI, RHI, TEMPI)
         ENDIF
      ENDIF
!cC
!cC *** IF METASTABLE AND NO WATER - RESOLVE AS NORMAL ********************
!cC
!c      IF (WATER.LE.TINY .AND. METSTBL.EQ.1) THEN
!c         METSTBL = 0
!c         GOTO 50
!c      ENDIF
!
! *** SAVE RESULTS TO ARRAYS (units = mole/m3) ****************************
!
      GAS(1) = GNH3                ! Gaseous aerosol species
      GAS(2) = GHNO3
      GAS(3) = GHCL
!
      DO 10 I=1,NIONS              ! Liquid aerosol species
         AERLIQ(I) = MOLAL(I)
  10  CONTINUE
      DO 20 I=1,NGASAQ
         AERLIQ(NIONS+1+I) = GASAQ(I)
  20  CONTINUE
      AERLIQ(NIONS+1)        = WATER*1.0D3/18.0D0
      AERLIQ(NIONS+NGASAQ+2) = COH
!
      AERSLD(1) = CNANO3           ! Solid aerosol species
      AERSLD(2) = CNH4NO3
      AERSLD(3) = CNACL
      AERSLD(4) = CNH4CL
      AERSLD(5) = CNA2SO4
      AERSLD(6) = CNH42S4
      AERSLD(7) = CNAHSO4
      AERSLD(8) = CNH4HS4
      AERSLD(9) = CLC
!
      IF(WATER.LE.TINY) THEN       ! Dry flag
        OTHER(1) = 1.d0
      ELSE
        OTHER(1) = 0.d0
      ENDIF
!
      OTHER(2) = SULRAT            ! Other stuff
      OTHER(3) = SULRATW
      OTHER(4) = SODRAT
      OTHER(5) = IONIC
      OTHER(6) = ICLACT
!
      SCASI = SCASE
!
      WT(1) = WI(1)                ! Total gas+aerosol phase
      WT(2) = WI(2)
      WT(3) = WI(3)
      WT(4) = WI(4)
      WT(5) = WI(5)
      IF (IPROB.GT.0 .AND. WATER.GT.TINY) THEN
         WT(3) = WT(3) + GNH3
         WT(4) = WT(4) + GHNO3
         WT(5) = WT(5) + GHCL
      ENDIF
!
      IF (LHOOK) CALL DR_HOOK('ISOROPIA',1,ZHOOK_HANDLE)
!
! *** END OF SUBROUTINE ISOROPIA ******************************************
!
      ENDSUBROUTINE ISOROPIA
!=======================================================================
!
! *** ISORROPIA CODE
! *** SUBROUTINE SETPARM
! *** THIS SUBROUTINE REDEFINES THE SOLUTION PARAMETERS OF ISORROPIA
!
! ======================== ARGUMENTS / USAGE ===========================
!
! *** NOTE: IF NEGATIVE VALUES ARE GIVEN FOR A PARAMETER, IT IS
!     IGNORED AND THE CURRENT VALUE IS USED INSTEAD.
!
!  INPUT:
!  1. [WFTYPI]
!     INTEGER variable.
!     Defines the type of weighting algorithm for the solution in Mutual
!     Deliquescence Regions (MDR's):
!     0 - MDR's are assumed dry. This is equivalent to the approach
!         used by SEQUILIB.
!     1 - The solution is assumed "half" dry and "half" wet throughout
!         the MDR.
!     2 - The solution is a relative-humidity weighted mean of the
!         dry and wet solutions (as defined in Nenes et al., 1998)
!
!  2. [IACALCI]
!     INTEGER variable.
!     Method of activity coefficient calculation:
!     0 - Calculate coefficients during runtime
!     1 - Use precalculated tables
!
!  3. [EPSI]
!     DOUBLE PRECITION variable.
!     Defines the convergence criterion for all iterative processes
!     in ISORROPIA, except those for activity coefficient calculations
!     (EPSACTI controls that).
!
!  4. [MAXITI]
!     INTEGER variable.
!     Defines the maximum number of iterations for all iterative
!     processes in ISORROPIA, except for activity coefficient calculations
!     (NSWEEPI controls that).
!
!  5. [NSWEEPI]
!     INTEGER variable.
!     Defines the maximum number of iterations for activity coefficient
!     calculations.
!
!  6. [EPSACTI]
!     DOUBLE PRECISION variable.
!     Defines the convergence criterion for activity coefficient
!     calculations.
!
!  7. [NDIV]
!     INTEGER variable.
!     Defines the number of subdivisions needed for the initial root
!     tracking for the bisection method. Usually this parameter should
!     not be altered, but is included for completeness.
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!=======================================================================
!
!=======================================================================
!
! *** ISORROPIA CODE
! *** BLOCK DATA BLKISO
! *** THIS SUBROUTINE PROVIDES INITIAL (DEFAULT) VALUES TO PROGRAM
!     PARAMETERS VIA DATA STATEMENTS
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!=======================================================================
!
      BLOCK DATA BLKISO
      IMPLICIT NONE
      INCLUDE 'isrpia.inc'
!
! *** DEFAULT VALUES *************************************************
!
      DATA TEMP/298.0/, R/82.0567D-6/, RH/0.9D0/, EPS/1D-6/, MAXIT/100/,&
     &     TINY/1D-20/, GREAT/1D10/, ZERO/0.0D0/, ONE/1.0D0/,NSWEEP/4/, &
     &     TINY2/1D-11/,NDIV/5/
!
      DATA MOLAL/NIONS*0.0D0/, MOLALR/NPAIR*0.0D0/, GAMA/NPAIR*0.1D0/,  &
     &     GAMOU/NPAIR*1D10/,  GAMIN/NPAIR*1D10/,   CALAIN/.TRUE./,     &
     &     CALAOU/.TRUE./,     EPSACT/5D-2/,        ICLACT/0/,          &
     &     IACALC/1/,          WFTYP/2/
!
      DATA ERRSTK/NERRMX*0/,   ERRMSG/NERRMX*' '/,  NOFER/0/,           &
     &     STKOFL/.FALSE./
!
      DATA IPROB/0/, METSTBL/0/
!
      DATA VERSION /'1.5 (12/03/03)'/
!
! *** OTHER PARAMETERS ***********************************************
!
      DATA SMW/58.5,142.,85.0,132.,80.0,53.5,98.0,98.0,115.,63.0,       &
     &         36.5,120.,247./                                          &
     &     IMW/ 1.0,23.0,18.0,35.5,96.0,97.0,63.0/,                     &
     &     WMW/23.0,98.0,17.0,63.0,36.5/
!
      DATA ZZ/1,2,1,2,1,1,2,1,1,1,1,1,2/, Z /1,1,1,1,2,1,1/
!
! *** ZSR RELATIONSHIP PARAMETERS **************************************
!
! awas= ammonium sulfate
!
      DATA AWAS/33*100.,30,30,30,29.54,28.25,27.06,25.94,               &
     & 24.89,23.90,22.97,22.10,21.27,20.48,19.73,19.02,18.34,17.69,     &
     & 17.07,16.48,15.91,15.37,14.85,14.34,13.86,13.39,12.94,12.50,     &
     & 12.08,11.67,11.27,10.88,10.51,10.14, 9.79, 9.44, 9.10, 8.78,     &
     &  8.45, 8.14, 7.83, 7.53, 7.23, 6.94, 6.65, 6.36, 6.08, 5.81,     &
     &  5.53, 5.26, 4.99, 4.72, 4.46, 4.19, 3.92, 3.65, 3.38, 3.11,     &
     &  2.83, 2.54, 2.25, 1.95, 1.63, 1.31, 0.97, 0.63, 0.30, 0.001/
!
! awsn= sodium nitrate
!
      DATA AWSN/ 9*1.e5,685.59,                                         &
     & 451.00,336.46,268.48,223.41,191.28,                              &
     & 167.20,148.46,133.44,121.12,110.83,                              &
     & 102.09,94.57,88.03,82.29,77.20,72.65,68.56,64.87,61.51,58.44,    &
     & 55.62,53.03,50.63,48.40,46.32,44.39,42.57,40.87,39.27,37.76,     &
     & 36.33,34.98,33.70,32.48,31.32,30.21,29.16,28.14,27.18,26.25,     &
     & 25.35,24.50,23.67,22.87,22.11,21.36,20.65,19.95,19.28,18.62,     &
     & 17.99,17.37,16.77,16.18,15.61,15.05,14.51,13.98,13.45,12.94,     &
     & 12.44,11.94,11.46,10.98,10.51,10.04, 9.58, 9.12, 8.67, 8.22,     &
     &  7.77, 7.32, 6.88, 6.43, 5.98, 5.53, 5.07, 4.61, 4.15, 3.69,     &
     &  3.22, 2.76, 2.31, 1.87, 1.47, 1.10, 0.77, 0.48, 0.23, 0.001/
!
! awsc= sodium chloride
!
      DATA AWSC/                                                        &
     &  100., 100., 100., 100., 100., 100., 100., 100., 100., 100.,     &
     &  100., 100., 100., 100., 100., 100., 100., 100., 100.,16.34,     &
     & 16.28,16.22,16.15,16.09,16.02,15.95,15.88,15.80,15.72,15.64,     &
     & 15.55,15.45,15.36,15.25,15.14,15.02,14.89,14.75,14.60,14.43,     &
     & 14.25,14.04,13.81,13.55,13.25,12.92,12.56,12.19,11.82,11.47,     &
     & 11.13,10.82,10.53,10.26,10.00, 9.76, 9.53, 9.30, 9.09, 8.88,     &
     &  8.67, 8.48, 8.28, 8.09, 7.90, 7.72, 7.54, 7.36, 7.17, 6.99,     &
     &  6.81, 6.63, 6.45, 6.27, 6.09, 5.91, 5.72, 5.53, 5.34, 5.14,     &
     &  4.94, 4.74, 4.53, 4.31, 4.09, 3.86, 3.62, 3.37, 3.12, 2.85,     &
     &  2.58, 2.30, 2.01, 1.72, 1.44, 1.16, 0.89, 0.64, 0.40, 0.18/
!
! awac= ammonium chloride
!
      DATA AWAC/                                                        &
     &  100., 100., 100., 100., 100., 100., 100., 100., 100., 100.,     &
     &  100., 100., 100., 100., 100., 100., 100., 100., 100.,31.45,     &
     & 31.30,31.14,30.98,30.82,30.65,30.48,30.30,30.11,29.92,29.71,     &
     & 29.50,29.29,29.06,28.82,28.57,28.30,28.03,27.78,27.78,27.77,     &
     & 27.77,27.43,27.07,26.67,26.21,25.73,25.18,24.56,23.84,23.01,     &
     & 22.05,20.97,19.85,18.77,17.78,16.89,16.10,15.39,14.74,14.14,     &
     & 13.59,13.06,12.56,12.09,11.65,11.22,10.81,10.42,10.03, 9.66,     &
     &  9.30, 8.94, 8.59, 8.25, 7.92, 7.59, 7.27, 6.95, 6.63, 6.32,     &
     &  6.01, 5.70, 5.39, 5.08, 4.78, 4.47, 4.17, 3.86, 3.56, 3.25,     &
     &  2.94, 2.62, 2.30, 1.98, 1.65, 1.32, 0.97, 0.62, 0.26, 0.13/
!
! awss= sodium sulfate
!
      DATA AWSS/34*1.e5,23*14.30,14.21,12.53,11.47,                     &
     & 10.66,10.01, 9.46, 8.99, 8.57, 8.19, 7.85, 7.54, 7.25, 6.98,     &
     &  6.74, 6.50, 6.29, 6.08, 5.88, 5.70, 5.52, 5.36, 5.20, 5.04,     &
     &  4.90, 4.75, 4.54, 4.34, 4.14, 3.93, 3.71, 3.49, 3.26, 3.02,     &
     &  2.76, 2.49, 2.20, 1.89, 1.55, 1.18, 0.82, 0.49, 0.22, 0.001/
!
! awab= ammonium bisulfate
!
      DATA AWAB/356.45,296.51,253.21,220.47,194.85,                     &
     & 174.24,157.31,143.16,131.15,120.82,                              &
     & 111.86,103.99,97.04,90.86,85.31,80.31,75.78,71.66,67.90,64.44,   &
     &  61.25,58.31,55.58,53.04,50.68,48.47,46.40,44.46,42.63,40.91,    &
     &  39.29,37.75,36.30,34.92,33.61,32.36,31.18,30.04,28.96,27.93,    &
     &  26.94,25.99,25.08,24.21,23.37,22.57,21.79,21.05,20.32,19.63,    &
     &  18.96,18.31,17.68,17.07,16.49,15.92,15.36,14.83,14.31,13.80,    &
     &  13.31,12.83,12.36,11.91,11.46,11.03,10.61,10.20, 9.80, 9.41,    &
     &   9.02, 8.64, 8.28, 7.91, 7.56, 7.21, 6.87, 6.54, 6.21, 5.88,    &
     &   5.56, 5.25, 4.94, 4.63, 4.33, 4.03, 3.73, 3.44, 3.14, 2.85,    &
     &   2.57, 2.28, 1.99, 1.71, 1.42, 1.14, 0.86, 0.57, 0.29, 0.001/
!
! awsa= sulfuric acid
!
      DATA AWSA/                                                        &
     & 34.0,33.56,29.22,26.55,24.61,23.11,21.89,20.87,19.99,            &
     & 19.21,18.51,17.87,17.29,16.76,16.26,15.8,15.37,14.95,14.56,      &
     & 14.20,13.85,13.53,13.22,12.93,12.66,12.40,12.14,11.90,11.67,     &
     & 11.44,11.22,11.01,10.8,10.60,10.4,10.2,10.01,9.83,9.65,9.47,     &
     & 9.3,9.13,8.96,8.81,8.64,8.48,8.33,8.17,8.02,7.87,7.72,7.58,      &
     & 7.44,7.30,7.16,7.02,6.88,6.75,6.61,6.48,6.35,6.21,6.08,5.95,     &
     & 5.82,5.69,5.56,5.44,5.31,5.18,5.05,4.92,4.79,4.66,4.53,4.40,     &
     & 4.27,4.14,4.,3.87,3.73,3.6,3.46,3.31,3.17,3.02,2.87,2.72,        &
     & 2.56,2.4,2.23,2.05,1.87,1.68,1.48,1.27,1.05,0.807,0.552,0.281/
!
! awlc= (NH4)3H(SO4)2
!
      DATA AWLC/34*1.e5,17.0,16.5,15.94,15.31,14.71,14.14,              &
     & 13.60,13.08,12.59,12.12,11.68,11.25,10.84,10.44,10.07, 9.71,     &
     &  9.36, 9.02, 8.70, 8.39, 8.09, 7.80, 7.52, 7.25, 6.99, 6.73,     &
     &  6.49, 6.25, 6.02, 5.79, 5.57, 5.36, 5.15, 4.95, 4.76, 4.56,     &
     &  4.38, 4.20, 4.02, 3.84, 3.67, 3.51, 3.34, 3.18, 3.02, 2.87,     &
     &  2.72, 2.57, 2.42, 2.28, 2.13, 1.99, 1.85, 1.71, 1.57, 1.43,     &
     &  1.30, 1.16, 1.02, 0.89, 0.75, 0.61, 0.46, 0.32, 0.16, 0.001/
!
! awan= ammonium nitrate
!
      DATA AWAN/31*1.e5,                                                &
     &       97.17,92.28,87.66,83.15,78.87,74.84,70.98,67.46,64.11,     &
     & 60.98,58.07,55.37,52.85,50.43,48.24,46.19,44.26,42.40,40.70,     &
     & 39.10,37.54,36.10,34.69,33.35,32.11,30.89,29.71,28.58,27.46,     &
     & 26.42,25.37,24.33,23.89,22.42,21.48,20.56,19.65,18.76,17.91,     &
     & 17.05,16.23,15.40,14.61,13.82,13.03,12.30,11.55,10.83,10.14,     &
     &  9.44, 8.79, 8.13, 7.51, 6.91, 6.32, 5.75, 5.18, 4.65, 4.14,     &
     &  3.65, 3.16, 2.71, 2.26, 1.83, 1.42, 1.03, 0.66, 0.30, 0.001/
!
! awsb= sodium bisulfate
!
      DATA AWSB/173.72,156.88,142.80,130.85,120.57,                     &
     & 111.64,103.80,96.88,90.71,85.18,                                 &
     & 80.20,75.69,71.58,67.82,64.37,61.19,58.26,55.53,53.00,50.64,     &
     & 48.44,46.37,44.44,42.61,40.90,39.27,37.74,36.29,34.91,33.61,     &
     & 32.36,31.18,30.05,28.97,27.94,26.95,26.00,25.10,24.23,23.39,     &
     & 22.59,21.81,21.07,20.35,19.65,18.98,18.34,17.71,17.11,16.52,     &
     & 15.95,15.40,14.87,14.35,13.85,13.36,12.88,12.42,11.97,11.53,     &
     & 11.10,10.69,10.28, 9.88, 9.49, 9.12, 8.75, 8.38, 8.03, 7.68,     &
     &  7.34, 7.01, 6.69, 6.37, 6.06, 5.75, 5.45, 5.15, 4.86, 4.58,     &
     &  4.30, 4.02, 3.76, 3.49, 3.23, 2.98, 2.73, 2.48, 2.24, 2.01,     &
     &  1.78, 1.56, 1.34, 1.13, 0.92, 0.73, 0.53, 0.35, 0.17, 0.001/
!
! *** END OF BLOCK DATA SUBPROGRAM *************************************
!
      ENDBLOCKDATA BLKISO
!=======================================================================
!
! *** ISORROPIA CODE
! *** SUBROUTINE INIT1
! *** THIS SUBROUTINE INITIALIZES ALL GLOBAL VARIABLES FOR AMMONIUM
!     SULFATE AEROSOL SYSTEMS (SUBROUTINE ISRP1)
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!=======================================================================
!
      SUBROUTINE INIT1 (WI, RHI, TEMPI)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
      IMPLICIT NONE
      INCLUDE 'isrpia.inc'
      REAL(KIND=8) :: WI,RHI,TEMPI,COEF,G130,G13I,T0,T0T,TCF
      INTEGER :: I,IRH
      DIMENSION WI(NCOMP)
      REAL      IC,GII,GI0,XX,LN10
      PARAMETER (LN10=2.3025851)
!
! *** SAVE INPUT VARIABLES IN COMMON BLOCK ******************************
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('INIT1',0,ZHOOK_HANDLE)
      IF (IPROB.EQ.0) THEN                 ! FORWARD CALCULATION
         DO 10 I=1,NCOMP
            W(I) = MAX(WI(I), TINY)
10       CONTINUE
      ELSE
         DO 15 I=1,NCOMP                   ! REVERSE CALCULATION
            WAER(I) = MAX(WI(I), TINY)
            W(I)    = ZERO
15       CONTINUE
      ENDIF
      RH      = RHI
      TEMP    = TEMPI
!
! *** CALCULATE EQUILIBRIUM CONSTANTS ***********************************
!
      XK1  = 1.015e-2  ! HSO4(aq)         <==> H(aq)     + SO4(aq)
      XK21 = 57.639    ! NH3(g)           <==> NH3(aq)
      XK22 = 1.805e-5  ! NH3(aq)          <==> NH4(aq)   + OH(aq)
      XK7  = 1.817     ! (NH4)2SO4(s)     <==> 2*NH4(aq) + SO4(aq)
      XK12 = 1.382e2   ! NH4HSO4(s)       <==> NH4(aq)   + HSO4(aq)
      XK13 = 29.268    ! (NH4)3H(SO4)2(s) <==> 3*NH4(aq) + HSO4(aq) + SO4(aq)
      XKW  = 1.010e-14 ! H2O              <==> H(aq)     + OH(aq)
!
      IF (INT(TEMP) .NE. 298) THEN   ! FOR T != 298K or 298.15K
         T0  = 298.15
         T0T = T0/TEMP
         COEF= 1.0+LOG(T0T)-T0T
         XK1 = XK1 *EXP(  8.85*(T0T-1.0) + 25.140*COEF)
         XK21= XK21*EXP( 13.79*(T0T-1.0) -  5.393*COEF)
         XK22= XK22*EXP( -1.50*(T0T-1.0) + 26.920*COEF)
         XK7 = XK7 *EXP( -2.65*(T0T-1.0) + 38.570*COEF)
         XK12= XK12*EXP( -2.87*(T0T-1.0) + 15.830*COEF)
         XK13= XK13*EXP( -5.19*(T0T-1.0) + 54.400*COEF)
         XKW = XKW *EXP(-22.52*(T0T-1.0) + 26.920*COEF)
      ENDIF
      XK2 = XK21*XK22
!
! *** CALCULATE DELIQUESCENCE RELATIVE HUMIDITIES (UNICOMPONENT) ********
!
      DRH2SO4  = 0.0000D0
      DRNH42S4 = 0.7997D0
      DRNH4HS4 = 0.4000D0
      DRLC     = 0.6900D0
      IF (INT(TEMP) .NE. 298) THEN
         T0       = 298.15d0
         TCF      = 1.0/TEMP - 1.0/T0
         DRNH42S4 = DRNH42S4*EXP( 80.*TCF)
         DRNH4HS4 = DRNH4HS4*EXP(384.*TCF)
         DRLC     = DRLC    *EXP(186.*TCF)
      ENDIF
!
! *** CALCULATE MUTUAL DELIQUESCENCE RELATIVE HUMIDITIES ****************
!
      DRMLCAB = 0.3780D0              ! (NH4)3H(SO4)2 & NH4HSO4
      DRMLCAS = 0.6900D0              ! (NH4)3H(SO4)2 & (NH4)2SO4
!CC      IF (INT(TEMP) .NE. 298) THEN      ! For the time being.
!CC         T0       = 298.15d0
!CC         TCF      = 1.0/TEMP - 1.0/T0
!CC         DRMLCAB  = DRMLCAB*EXP(507.506*TCF)
!CC         DRMLCAS  = DRMLCAS*EXP(133.865*TCF)
!CC      ENDIF
!
! *** LIQUID PHASE ******************************************************
!
      CHNO3  = ZERO
      CHCL   = ZERO
      CH2SO4 = ZERO
      COH    = ZERO
      WATER  = TINY
!
      DO 20 I=1,NPAIR
         MOLALR(I)=ZERO
         GAMA(I)  =0.1
         GAMIN(I) =GREAT
         GAMOU(I) =GREAT
         M0(I)    =1d5
 20   CONTINUE
!
      DO 30 I=1,NPAIR
         GAMA(I) = 0.1d0
 30   CONTINUE
!
      DO 40 I=1,NIONS
         MOLAL(I)=ZERO
40    CONTINUE
      COH = ZERO
!
      DO 50 I=1,NGASAQ
         GASAQ(I)=ZERO
50    CONTINUE
!
! *** SOLID PHASE *******************************************************
!
      CNH42S4= ZERO
      CNH4HS4= ZERO
      CNACL  = ZERO
      CNA2SO4= ZERO
      CNANO3 = ZERO
      CNH4NO3= ZERO
      CNH4CL = ZERO
      CNAHSO4= ZERO
      CLC    = ZERO
!
! *** GAS PHASE *********************************************************
!
      GNH3   = ZERO
      GHNO3  = ZERO
      GHCL   = ZERO
!
! *** CALCULATE ZSR PARAMETERS ******************************************
!
      IRH    = MIN (INT(RH*NZSR+0.5),NZSR)  ! Position in ZSR arrays
      IRH    = MAX (IRH, 1)
!
      M0(01) = AWSC(IRH)      ! NACl
      IF (M0(01) .LT. 100.0) THEN
         IC = M0(01)
         CALL KMTAB(IC,298.0,     GI0,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX)
         CALL KMTAB(IC,SNGL(TEMP),GII,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX)
         M0(01) = M0(01)*EXP(LN10*(GI0-GII))
      ENDIF
!
      M0(02) = AWSS(IRH)      ! (NA)2SO4
      IF (M0(02) .LT. 100.0) THEN
         IC = 3.0*M0(02)
         CALL KMTAB(IC,298.0,     XX,GI0,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX)
         CALL KMTAB(IC,SNGL(TEMP),XX,GII,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX)
         M0(02) = M0(02)*EXP(LN10*(GI0-GII))
      ENDIF
!
      M0(03) = AWSN(IRH)      ! NANO3
      IF (M0(03) .LT. 100.0) THEN
         IC = M0(03)
         CALL KMTAB(IC,298.0,     XX,XX,GI0,XX,XX,XX,XX,XX,XX,XX,XX,XX)
         CALL KMTAB(IC,SNGL(TEMP),XX,XX,GII,XX,XX,XX,XX,XX,XX,XX,XX,XX)
         M0(03) = M0(03)*EXP(LN10*(GI0-GII))
      ENDIF
!
      M0(04) = AWAS(IRH)      ! (NH4)2SO4
      IF (M0(04) .LT. 100.0) THEN
         IC = 3.0*M0(04)
         CALL KMTAB(IC,298.0,     XX,XX,XX,GI0,XX,XX,XX,XX,XX,XX,XX,XX)
         CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,GII,XX,XX,XX,XX,XX,XX,XX,XX)
         M0(04) = M0(04)*EXP(LN10*(GI0-GII))
      ENDIF
!
      M0(05) = AWAN(IRH)      ! NH4NO3
      IF (M0(05) .LT. 100.0) THEN
         IC     = M0(05)
         CALL KMTAB(IC,298.0,     XX,XX,XX,XX,GI0,XX,XX,XX,XX,XX,XX,XX)
         CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,XX,GII,XX,XX,XX,XX,XX,XX,XX)
         M0(05) = M0(05)*EXP(LN10*(GI0-GII))
      ENDIF
!
      M0(06) = AWAC(IRH)      ! NH4CL
      IF (M0(06) .LT. 100.0) THEN
         IC = M0(06)
         CALL KMTAB(IC,298.0,     XX,XX,XX,XX,XX,GI0,XX,XX,XX,XX,XX,XX)
         CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,XX,XX,GII,XX,XX,XX,XX,XX,XX)
         M0(06) = M0(06)*EXP(LN10*(GI0-GII))
      ENDIF
!
      M0(07) = AWSA(IRH)      ! 2H-SO4
      IF (M0(07) .LT. 100.0) THEN
         IC = 3.0*M0(07)
         CALL KMTAB(IC,298.0,     XX,XX,XX,XX,XX,XX,GI0,XX,XX,XX,XX,XX)
         CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,XX,XX,XX,GII,XX,XX,XX,XX,XX)
         M0(07) = M0(07)*EXP(LN10*(GI0-GII))
      ENDIF
!
      M0(08) = AWSA(IRH)      ! H-HSO4
!CC      IF (M0(08) .LT. 100.0) THEN     ! These are redundant, because M0(8) is not used
!CC         IC = M0(08)
!CC         CALL KMTAB(IC,298.0,     XX,XX,XX,XX,XX,XX,XX,GI0,XX,XX,XX,XX)
!CC         CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,XX,XX,XX,XX,GI0,XX,XX,XX,XX)
!CCCCC         CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,XX,XX,XX,XX,GII,XX,XX,XX,XX)
!CC         M0(08) = M0(08)*EXP(LN10*(GI0-GII))
!CC      ENDIF
!
      M0(09) = AWAB(IRH)      ! NH4HSO4
      IF (M0(09) .LT. 100.0) THEN
         IC = M0(09)
         CALL KMTAB(IC,298.0,     XX,XX,XX,XX,XX,XX,XX,XX,GI0,XX,XX,XX)
         CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,XX,XX,XX,XX,XX,GII,XX,XX,XX)
         M0(09) = M0(09)*EXP(LN10*(GI0-GII))
      ENDIF
!
      M0(12) = AWSB(IRH)      ! NAHSO4
      IF (M0(12) .LT. 100.0) THEN
         IC = M0(12)
         CALL KMTAB(IC,298.0,     XX,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX,GI0)
         CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX,GII)
         M0(12) = M0(12)*EXP(LN10*(GI0-GII))
      ENDIF
!
      M0(13) = AWLC(IRH)      ! (NH4)3H(SO4)2
      IF (M0(13) .LT. 100.0) THEN
         IC     = 4.0*M0(13)
         CALL KMTAB(IC,298.0,     XX,XX,XX,GI0,XX,XX,XX,XX,GII,XX,XX,XX)
         G130   = 0.2*(3.0*GI0+2.0*GII)
         CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,GI0,XX,XX,XX,XX,GII,XX,XX,XX)
         G13I   = 0.2*(3.0*GI0+2.0*GII)
         M0(13) = M0(13)*EXP(LN10*SNGL(G130-G13I))
      ENDIF
!
! *** OTHER INITIALIZATIONS *********************************************
!
      ICLACT  = 0
      CALAOU  = .TRUE.
      CALAIN  = .TRUE.
      FRST    = .TRUE.
      SCASE   = '?? '
      SULRATW = 2.D0
      SODRAT  = ZERO
      NOFER   = 0
      STKOFL  =.FALSE.
      DO 60 I=1,NERRMX
         ERRSTK(I) =-999
         ERRMSG(I) = 'MESSAGE N/A'
   60 CONTINUE
!
! *** END OF SUBROUTINE INIT1 *******************************************
!
      IF (LHOOK) CALL DR_HOOK('INIT1',1,ZHOOK_HANDLE)
      ENDSUBROUTINE INIT1



!=======================================================================
!
! *** ISORROPIA CODE
! *** SUBROUTINE INIT2
! *** THIS SUBROUTINE INITIALIZES ALL GLOBAL VARIABLES FOR AMMONIUM,
!     NITRATE, SULFATE AEROSOL SYSTEMS (SUBROUTINE ISRP2)
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!=======================================================================
!
      SUBROUTINE INIT2 (WI, RHI, TEMPI)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
      IMPLICIT NONE
      INCLUDE 'isrpia.inc'
      REAL(KIND=8) :: WI,RHI,TEMPI
      DIMENSION WI(NCOMP)
      REAL :: IC,GII,GI0,XX,LN10
      REAL(KIND=8) :: COEF,G130,G13I,T0,T0T,TCF
      PARAMETER (LN10=2.3025851)
      INTEGER :: I,IRH
!
! *** SAVE INPUT VARIABLES IN COMMON BLOCK ******************************
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('INIT2',0,ZHOOK_HANDLE)
      IF (IPROB.EQ.0) THEN                 ! FORWARD CALCULATION
         DO 10 I=1,NCOMP
            W(I) = MAX(WI(I), TINY)
10       CONTINUE
      ELSE
         DO 15 I=1,NCOMP                   ! REVERSE CALCULATION
            WAER(I) = MAX(WI(I), TINY)
            W(I)    = ZERO
15       CONTINUE
      ENDIF
      RH      = RHI
      TEMP    = TEMPI
!
! *** CALCULATE EQUILIBRIUM CONSTANTS ***********************************
!
      XK1  = 1.015e-2  ! HSO4(aq)         <==> H(aq)     + SO4(aq)
      XK21 = 57.639    ! NH3(g)           <==> NH3(aq)
      XK22 = 1.805e-5  ! NH3(aq)          <==> NH4(aq)   + OH(aq)
      XK4  = 2.511e6   ! HNO3(g)          <==> H(aq)     + NO3(aq) ! ISORR
!CC      XK4  = 3.638e6   ! HNO3(g)          <==> H(aq)     + NO3(aq) ! SEQUIL
      XK41 = 2.100e5   ! HNO3(g)          <==> HNO3(aq)
      XK7  = 1.817     ! (NH4)2SO4(s)     <==> 2*NH4(aq) + SO4(aq)
      XK10 = 5.746e-17 ! NH4NO3(s)        <==> NH3(g)    + HNO3(g) ! ISORR
!CC      XK10 = 2.985e-17 ! NH4NO3(s)        <==> NH3(g)    + HNO3(g) ! SEQUIL
      XK12 = 1.382e2   ! NH4HSO4(s)       <==> NH4(aq)   + HSO4(aq)
      XK13 = 29.268    ! (NH4)3H(SO4)2(s) <==> 3*NH4(aq) + HSO4(aq) + SO4(aq)
      XKW  = 1.010e-14 ! H2O              <==> H(aq)     + OH(aq)
!
      IF (INT(TEMP) .NE. 298) THEN   ! FOR T != 298K or 298.15K
         T0  = 298.15D0
         T0T = T0/TEMP
         COEF= 1.0+LOG(T0T)-T0T
         XK1 = XK1 *EXP(  8.85*(T0T-1.0) + 25.140*COEF)
         XK21= XK21*EXP( 13.79*(T0T-1.0) -  5.393*COEF)
         XK22= XK22*EXP( -1.50*(T0T-1.0) + 26.920*COEF)
         XK4 = XK4 *EXP( 29.17*(T0T-1.0) + 16.830*COEF) !ISORR
!CC         XK4 = XK4 *EXP( 29.47*(T0T-1.0) + 16.840*COEF) ! SEQUIL
         XK41= XK41*EXP( 29.17*(T0T-1.0) + 16.830*COEF)
         XK7 = XK7 *EXP( -2.65*(T0T-1.0) + 38.570*COEF)
         XK10= XK10*EXP(-74.38*(T0T-1.0) +  6.120*COEF) ! ISORR
!CC         XK10= XK10*EXP(-75.11*(T0T-1.0) + 13.460*COEF) ! SEQUIL
         XK12= XK12*EXP( -2.87*(T0T-1.0) + 15.830*COEF)
         XK13= XK13*EXP( -5.19*(T0T-1.0) + 54.400*COEF)
         XKW = XKW *EXP(-22.52*(T0T-1.0) + 26.920*COEF)
      ENDIF
      XK2  = XK21*XK22
      XK42 = XK4/XK41
!
! *** CALCULATE DELIQUESCENCE RELATIVE HUMIDITIES (UNICOMPONENT) ********
!
      DRH2SO4  = ZERO
      DRNH42S4 = 0.7997D0
      DRNH4HS4 = 0.4000D0
      DRNH4NO3 = 0.6183D0
      DRLC     = 0.6900D0
      IF (INT(TEMP) .NE. 298) THEN
         T0       = 298.15D0
         TCF      = 1.0/TEMP - 1.0/T0
         DRNH4NO3 = DRNH4NO3*EXP(852.*TCF)
         DRNH42S4 = DRNH42S4*EXP( 80.*TCF)
         DRNH4HS4 = DRNH4HS4*EXP(384.*TCF)
         DRLC     = DRLC    *EXP(186.*TCF)
      ENDIF
!
! *** CALCULATE MUTUAL DELIQUESCENCE RELATIVE HUMIDITIES ****************
!
      DRMLCAB = 0.3780D0              ! (NH4)3H(SO4)2 & NH4HSO4
      DRMLCAS = 0.6900D0              ! (NH4)3H(SO4)2 & (NH4)2SO4
      DRMASAN = 0.6000D0              ! (NH4)2SO4     & NH4NO3
!CC      IF (INT(TEMP) .NE. 298) THEN    ! For the time being
!CC         T0       = 298.15d0
!CC         TCF      = 1.0/TEMP - 1.0/T0
!CC         DRMLCAB  = DRMLCAB*EXP( 507.506*TCF)
!CC         DRMLCAS  = DRMLCAS*EXP( 133.865*TCF)
!CC         DRMASAN  = DRMASAN*EXP(1269.068*TCF)
!CC      ENDIF
!
! *** LIQUID PHASE ******************************************************
!
      CHNO3  = ZERO
      CHCL   = ZERO
      CH2SO4 = ZERO
      COH    = ZERO
      WATER  = TINY
!
      DO 20 I=1,NPAIR
         MOLALR(I)=ZERO
         GAMA(I)  =0.1
         GAMIN(I) =GREAT
         GAMOU(I) =GREAT
         M0(I)    =1d5
 20   CONTINUE
!
      DO 30 I=1,NPAIR
         GAMA(I) = 0.1d0
 30   CONTINUE
!
      DO 40 I=1,NIONS
         MOLAL(I)=ZERO
40    CONTINUE
      COH = ZERO
!
      DO 50 I=1,NGASAQ
         GASAQ(I)=ZERO
50    CONTINUE
!
! *** SOLID PHASE *******************************************************
!
      CNH42S4= ZERO
      CNH4HS4= ZERO
      CNACL  = ZERO
      CNA2SO4= ZERO
      CNANO3 = ZERO
      CNH4NO3= ZERO
      CNH4CL = ZERO
      CNAHSO4= ZERO
      CLC    = ZERO
!
! *** GAS PHASE *********************************************************
!
      GNH3   = ZERO
      GHNO3  = ZERO
      GHCL   = ZERO
!
! *** CALCULATE ZSR PARAMETERS ******************************************
!
      IRH    = MIN (INT(RH*NZSR+0.5),NZSR)  ! Position in ZSR arrays
      IRH    = MAX (IRH, 1)
!
      M0(01) = AWSC(IRH)      ! NACl
      IF (M0(01) .LT. 100.0) THEN
         IC = M0(01)
         CALL KMTAB(IC,298.0,     GI0,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX)
         CALL KMTAB(IC,SNGL(TEMP),GII,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX)
         M0(01) = M0(01)*EXP(LN10*(GI0-GII))
      ENDIF
!
      M0(02) = AWSS(IRH)      ! (NA)2SO4
      IF (M0(02) .LT. 100.0) THEN
         IC = 3.0*M0(02)
         CALL KMTAB(IC,298.0,     XX,GI0,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX)
         CALL KMTAB(IC,SNGL(TEMP),XX,GII,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX)
         M0(02) = M0(02)*EXP(LN10*(GI0-GII))
      ENDIF
!
      M0(03) = AWSN(IRH)      ! NANO3
      IF (M0(03) .LT. 100.0) THEN
         IC = M0(03)
         CALL KMTAB(IC,298.0,     XX,XX,GI0,XX,XX,XX,XX,XX,XX,XX,XX,XX)
         CALL KMTAB(IC,SNGL(TEMP),XX,XX,GII,XX,XX,XX,XX,XX,XX,XX,XX,XX)
         M0(03) = M0(03)*EXP(LN10*(GI0-GII))
      ENDIF
!
      M0(04) = AWAS(IRH)      ! (NH4)2SO4
      IF (M0(04) .LT. 100.0) THEN
         IC = 3.0*M0(04)
         CALL KMTAB(IC,298.0,     XX,XX,XX,GI0,XX,XX,XX,XX,XX,XX,XX,XX)
         CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,GII,XX,XX,XX,XX,XX,XX,XX,XX)
         M0(04) = M0(04)*EXP(LN10*(GI0-GII))
      ENDIF
!
      M0(05) = AWAN(IRH)      ! NH4NO3
      IF (M0(05) .LT. 100.0) THEN
         IC     = M0(05)
         CALL KMTAB(IC,298.0,     XX,XX,XX,XX,GI0,XX,XX,XX,XX,XX,XX,XX)
         CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,XX,GII,XX,XX,XX,XX,XX,XX,XX)
         M0(05) = M0(05)*EXP(LN10*(GI0-GII))
      ENDIF
!
      M0(06) = AWAC(IRH)      ! NH4CL
      IF (M0(06) .LT. 100.0) THEN
         IC = M0(06)
         CALL KMTAB(IC,298.0,     XX,XX,XX,XX,XX,GI0,XX,XX,XX,XX,XX,XX)
         CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,XX,XX,GII,XX,XX,XX,XX,XX,XX)
         M0(06) = M0(06)*EXP(LN10*(GI0-GII))
      ENDIF
!
      M0(07) = AWSA(IRH)      ! 2H-SO4
      IF (M0(07) .LT. 100.0) THEN
         IC = 3.0*M0(07)
         CALL KMTAB(IC,298.0,     XX,XX,XX,XX,XX,XX,GI0,XX,XX,XX,XX,XX)
         CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,XX,XX,XX,GII,XX,XX,XX,XX,XX)
         M0(07) = M0(07)*EXP(LN10*(GI0-GII))
      ENDIF
!
      M0(08) = AWSA(IRH)      ! H-HSO4
!CC      IF (M0(08) .LT. 100.0) THEN     ! These are redundant, because M0(8) is not used
!CC         IC = M0(08)
!CC         CALL KMTAB(IC,298.0,     XX,XX,XX,XX,XX,XX,XX,GI0,XX,XX,XX,XX)
!CC         CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,XX,XX,XX,XX,GI0,XX,XX,XX,XX)
!CCCCC         CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,XX,XX,XX,XX,GII,XX,XX,XX,XX)
!CC         M0(08) = M0(08)*EXP(LN10*(GI0-GII))
!CC      ENDIF
!
      M0(09) = AWAB(IRH)      ! NH4HSO4
      IF (M0(09) .LT. 100.0) THEN
         IC = M0(09)
         CALL KMTAB(IC,298.0,     XX,XX,XX,XX,XX,XX,XX,XX,GI0,XX,XX,XX)
         CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,XX,XX,XX,XX,XX,GII,XX,XX,XX)
         M0(09) = M0(09)*EXP(LN10*(GI0-GII))
      ENDIF
!
      M0(12) = AWSB(IRH)      ! NAHSO4
      IF (M0(12) .LT. 100.0) THEN
         IC = M0(12)
         CALL KMTAB(IC,298.0,     XX,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX,GI0)
         CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX,GII)
         M0(12) = M0(12)*EXP(LN10*(GI0-GII))
      ENDIF
!
      M0(13) = AWLC(IRH)      ! (NH4)3H(SO4)2
      IF (M0(13) .LT. 100.0) THEN
         IC     = 4.0*M0(13)
         CALL KMTAB(IC,298.0,     XX,XX,XX,GI0,XX,XX,XX,XX,GII,XX,XX,XX)
         G130   = 0.2*(3.0*GI0+2.0*GII)
         CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,GI0,XX,XX,XX,XX,GII,XX,XX,XX)
         G13I   = 0.2*(3.0*GI0+2.0*GII)
         M0(13) = M0(13)*EXP(LN10*SNGL(G130-G13I))
      ENDIF
!
! *** OTHER INITIALIZATIONS *********************************************
!
      ICLACT  = 0
      CALAOU  = .TRUE.
      CALAIN  = .TRUE.
      FRST    = .TRUE.
      SCASE   = '?? '
      SULRATW = 2.D0
      SODRAT  = ZERO
      NOFER   = 0
      STKOFL  =.FALSE.
      DO 60 I=1,NERRMX
         ERRSTK(I) =-999
         ERRMSG(I) = 'MESSAGE N/A'
   60 CONTINUE
!
! *** END OF SUBROUTINE INIT2 *******************************************
!
      IF (LHOOK) CALL DR_HOOK('INIT2',1,ZHOOK_HANDLE)
      ENDSUBROUTINE INIT2





!=======================================================================
!
! *** ISORROPIA CODE
!=======================================================================
!
      SUBROUTINE INIT3 (WI, RHI, TEMPI)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
      IMPLICIT NONE
      INCLUDE 'isrpia.inc'
      REAL(KIND=8) :: WI,RHI,TEMPI
      DIMENSION WI(NCOMP)
      REAL :: IC,GII,GI0,XX,LN10
      REAL(KIND=8) :: COEF,G130,G13I,T0,T0T,TCF
      PARAMETER (LN10=2.3025851)
      INTEGER :: I,IRH
!
! *** SAVE INPUT VARIABLES IN COMMON BLOCK ******************************
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('INIT3',0,ZHOOK_HANDLE)
      IF (IPROB.EQ.0) THEN                 ! FORWARD CALCULATION
         DO 10 I=1,NCOMP
            W(I) = MAX(WI(I), TINY)
10       CONTINUE
      ELSE
         DO 15 I=1,NCOMP                   ! REVERSE CALCULATION
            WAER(I) = MAX(WI(I), TINY)
            W(I)    = ZERO
15       CONTINUE
      ENDIF
      RH      = RHI
      TEMP    = TEMPI
!
! *** CALCULATE EQUILIBRIUM CONSTANTS ***********************************
!
      XK1  = 1.015D-2  ! HSO4(aq)         <==> H(aq)     + SO4(aq)
      XK21 = 57.639D0  ! NH3(g)           <==> NH3(aq)
      XK22 = 1.805D-5  ! NH3(aq)          <==> NH4(aq)   + OH(aq)
      XK3  = 1.971D6   ! HCL(g)           <==> H(aq)     + CL(aq)
      XK31 = 2.500e3   ! HCL(g)           <==> HCL(aq)
      XK4  = 2.511e6   ! HNO3(g)          <==> H(aq)     + NO3(aq) ! ISORR
!CC      XK4  = 3.638e6   ! HNO3(g)          <==> H(aq)     + NO3(aq) ! SEQUIL
      XK41 = 2.100e5   ! HNO3(g)          <==> HNO3(aq)
      XK5  = 0.4799D0  ! NA2SO4(s)        <==> 2*NA(aq)  + SO4(aq)
      XK6  = 1.086D-16 ! NH4CL(s)         <==> NH3(g)    + HCL(g)
      XK7  = 1.817D0   ! (NH4)2SO4(s)     <==> 2*NH4(aq) + SO4(aq)
      XK8  = 37.661D0  ! NACL(s)          <==> NA(aq)    + CL(aq)
      XK10 = 5.746D-17 ! NH4NO3(s)        <==> NH3(g)    + HNO3(g) ! ISORR
!CC      XK10 = 2.985e-17 ! NH4NO3(s)        <==> NH3(g)    + HNO3(g) ! SEQUIL
      XK11 = 2.413D4   ! NAHSO4(s)        <==> NA(aq)    + HSO4(aq)
      XK12 = 1.382D2   ! NH4HSO4(s)       <==> NH4(aq)   + HSO4(aq)
      XK13 = 29.268D0  ! (NH4)3H(SO4)2(s) <==> 3*NH4(aq) + HSO4(aq) + SO4(aq)
      XK14 = 22.05D0   ! NH4CL(s)         <==> NH4(aq)   + CL(aq)
      XKW  = 1.010D-14 ! H2O              <==> H(aq)     + OH(aq)
      XK9  = 11.977D0  ! NANO3(s)         <==> NA(aq)    + NO3(aq)
!
      IF (INT(TEMP) .NE. 298) THEN   ! FOR T != 298K or 298.15K
         T0  = 298.15D0
         T0T = T0/TEMP
         COEF= 1.0+LOG(T0T)-T0T
         XK1 = XK1 *EXP(  8.85*(T0T-1.0) + 25.140*COEF)
         XK21= XK21*EXP( 13.79*(T0T-1.0) -  5.393*COEF)
         XK22= XK22*EXP( -1.50*(T0T-1.0) + 26.920*COEF)
         XK3 = XK3 *EXP( 30.20*(T0T-1.0) + 19.910*COEF)
         XK31= XK31*EXP( 30.20*(T0T-1.0) + 19.910*COEF)
         XK4 = XK4 *EXP( 29.17*(T0T-1.0) + 16.830*COEF) !ISORR
!CC         XK4 = XK4 *EXP( 29.47*(T0T-1.0) + 16.840*COEF) ! SEQUIL
         XK41= XK41*EXP( 29.17*(T0T-1.0) + 16.830*COEF)
         XK5 = XK5 *EXP(  0.98*(T0T-1.0) + 39.500*COEF)
         XK6 = XK6 *EXP(-71.00*(T0T-1.0) +  2.400*COEF)
         XK7 = XK7 *EXP( -2.65*(T0T-1.0) + 38.570*COEF)
         XK8 = XK8 *EXP( -1.56*(T0T-1.0) + 16.900*COEF)
         XK9 = XK9 *EXP( -8.22*(T0T-1.0) + 16.010*COEF)
         XK10= XK10*EXP(-74.38*(T0T-1.0) +  6.120*COEF) ! ISORR
!CC         XK10= XK10*EXP(-75.11*(T0T-1.0) + 13.460*COEF) ! SEQUIL
         XK11= XK11*EXP(  0.79*(T0T-1.0) + 14.746*COEF)
         XK12= XK12*EXP( -2.87*(T0T-1.0) + 15.830*COEF)
         XK13= XK13*EXP( -5.19*(T0T-1.0) + 54.400*COEF)
         XK14= XK14*EXP( 24.55*(T0T-1.0) + 16.900*COEF)
         XKW = XKW *EXP(-22.52*(T0T-1.0) + 26.920*COEF)
      ENDIF
      XK2  = XK21*XK22
      XK42 = XK4/XK41
      XK32 = XK3/XK31
!
! *** CALCULATE DELIQUESCENCE RELATIVE HUMIDITIES (UNICOMPONENT) ********
!
      DRH2SO4  = ZERO
      DRNH42S4 = 0.7997D0
      DRNH4HS4 = 0.4000D0
      DRLC     = 0.6900D0
      DRNACL   = 0.7528D0
      DRNANO3  = 0.7379D0
      DRNH4CL  = 0.7710D0
      DRNH4NO3 = 0.6183D0
      DRNA2SO4 = 0.9300D0
      DRNAHSO4 = 0.5200D0
      IF (INT(TEMP) .NE. 298) THEN
         T0       = 298.15D0
         TCF      = 1.0/TEMP - 1.0/T0
         DRNACL   = DRNACL  *EXP( 25.*TCF)
         DRNANO3  = DRNANO3 *EXP(304.*TCF)
         DRNA2SO4 = DRNA2SO4*EXP( 80.*TCF)
         DRNH4NO3 = DRNH4NO3*EXP(852.*TCF)
         DRNH42S4 = DRNH42S4*EXP( 80.*TCF)
         DRNH4HS4 = DRNH4HS4*EXP(384.*TCF)
         DRLC     = DRLC    *EXP(186.*TCF)
         DRNH4CL  = DRNH4Cl *EXP(239.*TCF)
         DRNAHSO4 = DRNAHSO4*EXP(-45.*TCF)
      ENDIF
!
! *** CALCULATE MUTUAL DELIQUESCENCE RELATIVE HUMIDITIES ****************
!
      DRMLCAB = 0.378D0    ! (NH4)3H(SO4)2 & NH4HSO4
      DRMLCAS = 0.690D0    ! (NH4)3H(SO4)2 & (NH4)2SO4
      DRMASAN = 0.600D0    ! (NH4)2SO4     & NH4NO3
      DRMG1   = 0.460D0    ! (NH4)2SO4, NH4NO3, NA2SO4, NH4CL
      DRMG2   = 0.691D0    ! (NH4)2SO4, NA2SO4, NH4CL
      DRMG3   = 0.697D0    ! (NH4)2SO4, NA2SO4
      DRMH1   = 0.240D0    ! NA2SO4, NANO3, NACL, NH4NO3, NH4CL
      DRMH2   = 0.596D0    ! NA2SO4, NANO3, NACL, NH4CL
      DRMI1   = 0.240D0    ! LC, NAHSO4, NH4HSO4, NA2SO4, (NH4)2SO4
      DRMI2   = 0.363D0    ! LC, NAHSO4, NA2SO4, (NH4)2SO4  - NO DATA -
      DRMI3   = 0.610D0    ! LC, NA2SO4, (NH4)2SO4
      DRMQ1   = 0.494D0    ! (NH4)2SO4, NH4NO3, NA2SO4
      DRMR1   = 0.663D0    ! NA2SO4, NANO3, NACL
      DRMR2   = 0.735D0    ! NA2SO4, NACL
      DRMR3   = 0.673D0    ! NANO3, NACL
      DRMR4   = 0.694D0    ! NA2SO4, NACL, NH4CL
      DRMR5   = 0.731D0    ! NA2SO4, NH4CL
      DRMR6   = 0.596D0    ! NA2SO4, NANO3, NH4CL
      DRMR7   = 0.380D0    ! NA2SO4, NANO3, NACL, NH4NO3
      DRMR8   = 0.380D0    ! NA2SO4, NACL, NH4NO3
      DRMR9   = 0.494D0    ! NA2SO4, NH4NO3
      DRMR10  = 0.476D0    ! NA2SO4, NANO3, NH4NO3
      DRMR11  = 0.340D0    ! NA2SO4, NACL, NH4NO3, NH4CL
      DRMR12  = 0.460D0    ! NA2SO4, NH4NO3, NH4CL
      DRMR13  = 0.438D0    ! NA2SO4, NANO3, NH4NO3, NH4CL
!CC      IF (INT(TEMP) .NE. 298) THEN
!CC         T0       = 298.15d0
!CC         TCF      = 1.0/TEMP - 1.0/T0
!CC         DRMLCAB  = DRMLCAB*EXP( 507.506*TCF)
!CC         DRMLCAS  = DRMLCAS*EXP( 133.865*TCF)
!CC         DRMASAN  = DRMASAN*EXP(1269.068*TCF)
!CC         DRMG1    = DRMG1  *EXP( 572.207*TCF)
!CC         DRMG2    = DRMG2  *EXP(  58.166*TCF)
!CC         DRMG3    = DRMG3  *EXP(  22.253*TCF)
!CC         DRMH1    = DRMH1  *EXP(2116.542*TCF)
!CC         DRMH2    = DRMH2  *EXP( 650.549*TCF)
!CC         DRMI1    = DRMI1  *EXP( 565.743*TCF)
!CC         DRMI2    = DRMI2  *EXP(  91.745*TCF)
!CC         DRMI3    = DRMI3  *EXP( 161.272*TCF)
!CC         DRMQ1    = DRMQ1  *EXP(1616.621*TCF)
!CC         DRMR1    = DRMR1  *EXP( 292.564*TCF)
!CC         DRMR2    = DRMR2  *EXP(  14.587*TCF)
!CC         DRMR3    = DRMR3  *EXP( 307.907*TCF)
!CC         DRMR4    = DRMR4  *EXP(  97.605*TCF)
!CC         DRMR5    = DRMR5  *EXP(  98.523*TCF)
!CC         DRMR6    = DRMR6  *EXP( 465.500*TCF)
!CC         DRMR7    = DRMR7  *EXP( 324.425*TCF)
!CC         DRMR8    = DRMR8  *EXP(2660.184*TCF)
!CC         DRMR9    = DRMR9  *EXP(1617.178*TCF)
!CC         DRMR10   = DRMR10 *EXP(1745.226*TCF)
!CC         DRMR11   = DRMR11 *EXP(3691.328*TCF)
!CC         DRMR12   = DRMR12 *EXP(1836.842*TCF)
!CC         DRMR13   = DRMR13 *EXP(1967.938*TCF)
!CC      ENDIF
!
! *** LIQUID PHASE ******************************************************
!
      CHNO3  = ZERO
      CHCL   = ZERO
      CH2SO4 = ZERO
      COH    = ZERO
      WATER  = TINY
!
      DO 20 I=1,NPAIR
         MOLALR(I)=ZERO
         GAMA(I)  =0.1
         GAMIN(I) =GREAT
         GAMOU(I) =GREAT
         M0(I)    =1d5
 20   CONTINUE
!
      DO 30 I=1,NPAIR
         GAMA(I) = 0.1d0
 30   CONTINUE
!
      DO 40 I=1,NIONS
         MOLAL(I)=ZERO
40    CONTINUE
      COH = ZERO
!
      DO 50 I=1,NGASAQ
         GASAQ(I)=ZERO
50    CONTINUE
!
! *** SOLID PHASE *******************************************************
!
      CNH42S4= ZERO
      CNH4HS4= ZERO
      CNACL  = ZERO
      CNA2SO4= ZERO
      CNANO3 = ZERO
      CNH4NO3= ZERO
      CNH4CL = ZERO
      CNAHSO4= ZERO
      CLC    = ZERO
!
! *** GAS PHASE *********************************************************
!
      GNH3   = ZERO
      GHNO3  = ZERO
      GHCL   = ZERO
!
! *** CALCULATE ZSR PARAMETERS ******************************************
!
      IRH    = MIN (INT(RH*NZSR+0.5),NZSR)  ! Position in ZSR arrays
      IRH    = MAX (IRH, 1)
!
      M0(01) = AWSC(IRH)      ! NACl
      IF (M0(01) .LT. 100.0) THEN
         IC = M0(01)
         CALL KMTAB(IC,298.0,     GI0,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX)
         CALL KMTAB(IC,SNGL(TEMP),GII,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX)
         M0(01) = M0(01)*EXP(LN10*(GI0-GII))
      ENDIF
!
      M0(02) = AWSS(IRH)      ! (NA)2SO4
      IF (M0(02) .LT. 100.0) THEN
         IC = 3.0*M0(02)
         CALL KMTAB(IC,298.0,     XX,GI0,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX)
         CALL KMTAB(IC,SNGL(TEMP),XX,GII,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX)
         M0(02) = M0(02)*EXP(LN10*(GI0-GII))
      ENDIF
!
      M0(03) = AWSN(IRH)      ! NANO3
      IF (M0(03) .LT. 100.0) THEN
         IC = M0(03)
         CALL KMTAB(IC,298.0,     XX,XX,GI0,XX,XX,XX,XX,XX,XX,XX,XX,XX)
         CALL KMTAB(IC,SNGL(TEMP),XX,XX,GII,XX,XX,XX,XX,XX,XX,XX,XX,XX)
         M0(03) = M0(03)*EXP(LN10*(GI0-GII))
      ENDIF
!
      M0(04) = AWAS(IRH)      ! (NH4)2SO4
      IF (M0(04) .LT. 100.0) THEN
         IC = 3.0*M0(04)
         CALL KMTAB(IC,298.0,     XX,XX,XX,GI0,XX,XX,XX,XX,XX,XX,XX,XX)
         CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,GII,XX,XX,XX,XX,XX,XX,XX,XX)
         M0(04) = M0(04)*EXP(LN10*(GI0-GII))
      ENDIF
!
      M0(05) = AWAN(IRH)      ! NH4NO3
      IF (M0(05) .LT. 100.0) THEN
         IC     = M0(05)
         CALL KMTAB(IC,298.0,     XX,XX,XX,XX,GI0,XX,XX,XX,XX,XX,XX,XX)
         CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,XX,GII,XX,XX,XX,XX,XX,XX,XX)
         M0(05) = M0(05)*EXP(LN10*(GI0-GII))
      ENDIF
!
      M0(06) = AWAC(IRH)      ! NH4CL
      IF (M0(06) .LT. 100.0) THEN
         IC = M0(06)
         CALL KMTAB(IC,298.0,     XX,XX,XX,XX,XX,GI0,XX,XX,XX,XX,XX,XX)
         CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,XX,XX,GII,XX,XX,XX,XX,XX,XX)
         M0(06) = M0(06)*EXP(LN10*(GI0-GII))
      ENDIF
!
      M0(07) = AWSA(IRH)      ! 2H-SO4
      IF (M0(07) .LT. 100.0) THEN
         IC = 3.0*M0(07)
         CALL KMTAB(IC,298.0,     XX,XX,XX,XX,XX,XX,GI0,XX,XX,XX,XX,XX)
         CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,XX,XX,XX,GII,XX,XX,XX,XX,XX)
         M0(07) = M0(07)*EXP(LN10*(GI0-GII))
      ENDIF
!
      M0(08) = AWSA(IRH)      ! H-HSO4
!CC      IF (M0(08) .LT. 100.0) THEN     ! These are redundant, because M0(8) is not used
!CC         IC = M0(08)
!CC         CALL KMTAB(IC,298.0,     XX,XX,XX,XX,XX,XX,XX,GI0,XX,XX,XX,XX)
!CC         CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,XX,XX,XX,XX,GI0,XX,XX,XX,XX)
!CCCCC         CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,XX,XX,XX,XX,GII,XX,XX,XX,XX)
!CC         M0(08) = M0(08)*EXP(LN10*(GI0-GII))
!CC      ENDIF
!
      M0(09) = AWAB(IRH)      ! NH4HSO4
      IF (M0(09) .LT. 100.0) THEN
         IC = M0(09)
         CALL KMTAB(IC,298.0,     XX,XX,XX,XX,XX,XX,XX,XX,GI0,XX,XX,XX)
         CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,XX,XX,XX,XX,XX,GII,XX,XX,XX)
         M0(09) = M0(09)*EXP(LN10*(GI0-GII))
      ENDIF
!
      M0(12) = AWSB(IRH)      ! NAHSO4
      IF (M0(12) .LT. 100.0) THEN
         IC = M0(12)
         CALL KMTAB(IC,298.0,     XX,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX,GI0)
         CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX,GII)
         M0(12) = M0(12)*EXP(LN10*(GI0-GII))
      ENDIF
!
      M0(13) = AWLC(IRH)      ! (NH4)3H(SO4)2
      IF (M0(13) .LT. 100.0) THEN
         IC     = 4.0*M0(13)
         CALL KMTAB(IC,298.0,     XX,XX,XX,GI0,XX,XX,XX,XX,GII,XX,XX,XX)
         G130   = 0.2*(3.0*GI0+2.0*GII)
         CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,GI0,XX,XX,XX,XX,GII,XX,XX,XX)
         G13I   = 0.2*(3.0*GI0+2.0*GII)
         M0(13) = M0(13)*EXP(LN10*SNGL(G130-G13I))
      ENDIF
!
! *** OTHER INITIALIZATIONS *********************************************
!
      ICLACT  = 0
      CALAOU  = .TRUE.
      CALAIN  = .TRUE.
      FRST    = .TRUE.
      SCASE   = '?? '
      SULRATW = 2.D0
      NOFER   = 0
      STKOFL  =.FALSE.
      DO 60 I=1,NERRMX
         ERRSTK(I) =-999
         ERRMSG(I) = 'MESSAGE N/A'
   60 CONTINUE
!
! *** END OF SUBROUTINE INIT3 *******************************************
!
      IF (LHOOK) CALL DR_HOOK('INIT3',1,ZHOOK_HANDLE)
      ENDSUBROUTINE INIT3
!=======================================================================
!
! *** ISORROPIA CODE
! *** FUNCTION GETASR
! *** CALCULATES THE LIMITING NH4+/SO4 RATIO OF A SULFATE POOR SYSTEM
!     (i.e. SULFATE RATIO = 2.0) FOR GIVEN SO4 LEVEL AND RH
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!=======================================================================
!
!     DOUBLE PRECISION FUNCTION GETASR (SO4I, RHI)
                       FUNCTION GETASR (SO4I, RHI)
                       USE PARKIND1, ONLY : JPRB
                       USE YOMHOOK , ONLY : LHOOK, DR_HOOK
      IMPLICIT NONE
      INTEGER :: NSO4S,NRHS,NASRD,IA1,INDS,INDR,INDSL,INDSH,IPOSL,IPOSH
      REAL :: A1,RAT,WF
      REAL :: ASRAT,ASSO4
      PARAMETER (NSO4S=14, NRHS=20, NASRD=NSO4S*NRHS)
      COMMON /ASRC/ ASRAT(NASRD), ASSO4(NSO4S)
      REAL(KIND=8) :: SO4I, RHI
      REAL(KIND=8) :: GETASR
!CC
!CC *** SOLVE USING FULL COMPUTATIONS, NOT LOOK-UP TABLES **************
!CC
!CC         W(2) = WAER(2)
!CC         W(3) = WAER(2)*2.0001D0
!CC         CALL CALCA2
!CC         SULRATW = MOLAL(3)/WAER(2)
!CC         CALL INIT1 (WI, RHI, TEMPI)   ! Re-initialize COMMON BLOCK
!
! *** CALCULATE INDICES ************************************************
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('GETASR',0,ZHOOK_HANDLE)
      RAT    = SO4I/1.E-9
      A1     = INT(LOG10(RAT))                    ! Magnitude of RAT
      IA1    = INT(RAT/2.5/10.0**A1)
!
      INDS   = 4.0*A1 + MIN(IA1,4)
      INDS   = MIN(MAX(0, INDS), NSO4S-1) + 1     ! SO4 component of IPOS
!
      INDR   = INT(99.0-RHI*100.0) + 1
      INDR   = MIN(MAX(1, INDR), NRHS)            ! RH component of IPOS
!
! *** GET VALUE AND RETURN *********************************************
!
      INDSL  = INDS
      INDSH  = MIN(INDSL+1, NSO4S)
      IPOSL  = (INDSL-1)*NRHS + INDR              ! Low position in array
      IPOSH  = (INDSH-1)*NRHS + INDR              ! High position in array
!
      WF     = (SO4I-ASSO4(INDSL))/(ASSO4(INDSH)-ASSO4(INDSL) + 1e-7)
      WF     = MIN(MAX(WF, 0.0), 1.0)
!
      GETASR = WF*ASRAT(IPOSH) + (1.0-WF)*ASRAT(IPOSL)
!
! *** END OF FUNCTION GETASR *******************************************
!
      IF (LHOOK) CALL DR_HOOK('GETASR',1,ZHOOK_HANDLE)
      ENDFUNCTION GETASR



!=======================================================================
!
! *** ISORROPIA CODE
! *** BLOCK DATA AERSR
! *** CONTAINS DATA FOR AEROSOL SULFATE RATIO ARRAY NEEDED IN FUNCTION
!     GETASR
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!=======================================================================
!
      BLOCK DATA AERSR
      IMPLICIT NONE
      INTEGER :: NSO4S,NRHS,NASRD,I
      REAL :: ASRAT,ASSO4
      PARAMETER (NSO4S=14, NRHS=20, NASRD=NSO4S*NRHS)
      COMMON /ASRC/ ASRAT(NASRD), ASSO4(NSO4S)
!
      DATA ASSO4/1.0E-9, 2.5E-9, 5.0E-9, 7.5E-9, 1.0E-8,                &
     &           2.5E-8, 5.0E-8, 7.5E-8, 1.0E-7, 2.5E-7,                &
     &           5.0E-7, 7.5E-7, 1.0E-6, 5.0E-6/
!
      DATA (ASRAT(I), I=1,100)/                                         &
     & 1.020464, 0.9998130, 0.9960167, 0.9984423, 1.004004,             &
     & 1.010885,  1.018356,  1.026726,  1.034268, 1.043846,             &
     & 1.052933,  1.062230,  1.062213,  1.080050, 1.088350,             &
     & 1.096603,  1.104289,  1.111745,  1.094662, 1.121594,             &
     & 1.268909,  1.242444,  1.233815,  1.232088, 1.234020,             &
     & 1.238068,  1.243455,  1.250636,  1.258734, 1.267543,             &
     & 1.276948,  1.286642,  1.293337,  1.305592, 1.314726,             &
     & 1.323463,  1.333258,  1.343604,  1.344793, 1.355571,             &
     & 1.431463,  1.405204,  1.395791,  1.393190, 1.394403,             &
     & 1.398107,  1.403811,  1.411744,  1.420560, 1.429990,             &
     & 1.439742,  1.449507,  1.458986,  1.468403, 1.477394,             &
     & 1.487373,  1.495385,  1.503854,  1.512281, 1.520394,             &
     & 1.514464,  1.489699,  1.480686,  1.478187, 1.479446,             &
     & 1.483310,  1.489316,  1.497517,  1.506501, 1.515816,             &
     & 1.524724,  1.533950,  1.542758,  1.551730, 1.559587,             &
     & 1.568343,  1.575610,  1.583140,  1.590440, 1.596481,             &
     & 1.567743,  1.544426,  1.535928,  1.533645, 1.535016,             &
     & 1.539003,  1.545124,  1.553283,  1.561886, 1.570530,             &
     & 1.579234,  1.587813,  1.595956,  1.603901, 1.611349,             &
     & 1.618833,  1.625819,  1.632543,  1.639032, 1.645276/

      DATA (ASRAT(I), I=101,200)/                                       &
     & 1.707390,  1.689553,  1.683198,  1.681810, 1.683490,             &
     & 1.687477,  1.693148,  1.700084,  1.706917, 1.713507,             &
     & 1.719952,  1.726190,  1.731985,  1.737544, 1.742673,             &
     & 1.747756,  1.752431,  1.756890,  1.761141, 1.765190,             &
     & 1.785657,  1.771851,  1.767063,  1.766229, 1.767901,             &
     & 1.771455,  1.776223,  1.781769,  1.787065, 1.792081,             &
     & 1.796922,  1.801561,  1.805832,  1.809896, 1.813622,             &
     & 1.817292,  1.820651,  1.823841,  1.826871, 1.829745,             &
     & 1.822215,  1.810497,  1.806496,  1.805898, 1.807480,             &
     & 1.810684,  1.814860,  1.819613,  1.824093, 1.828306,             &
     & 1.832352,  1.836209,  1.839748,  1.843105, 1.846175,             &
     & 1.849192,  1.851948,  1.854574,  1.857038, 1.859387,             &
     & 1.844588,  1.834208,  1.830701,  1.830233, 1.831727,             &
     & 1.834665,  1.838429,  1.842658,  1.846615, 1.850321,             &
     & 1.853869,  1.857243,  1.860332,  1.863257, 1.865928,             &
     & 1.868550,  1.870942,  1.873208,  1.875355, 1.877389,             &
     & 1.899556,  1.892637,  1.890367,  1.890165, 1.891317,             &
     & 1.893436,  1.896036,  1.898872,  1.901485, 1.903908,             &
     & 1.906212,  1.908391,  1.910375,  1.912248, 1.913952,             &
     & 1.915621,  1.917140,  1.918576,  1.919934, 1.921220/

      DATA (ASRAT(I), I=201,280)/                                       &
     & 1.928264,  1.923245,  1.921625,  1.921523, 1.922421,             &
     & 1.924016,  1.925931,  1.927991,  1.929875, 1.931614,             &
     & 1.933262,  1.934816,  1.936229,  1.937560, 1.938769,             &
     & 1.939951,  1.941026,  1.942042,  1.943003, 1.943911,             &
     & 1.941205,  1.937060,  1.935734,  1.935666, 1.936430,             &
     & 1.937769,  1.939359,  1.941061,  1.942612, 1.944041,             &
     & 1.945393,  1.946666,  1.947823,  1.948911, 1.949900,             &
     & 1.950866,  1.951744,  1.952574,  1.953358, 1.954099,             &
     & 1.948985,  1.945372,  1.944221,  1.944171, 1.944850,             &
     & 1.946027,  1.947419,  1.948902,  1.950251, 1.951494,             &
     & 1.952668,  1.953773,  1.954776,  1.955719, 1.956576,             &
     & 1.957413,  1.958174,  1.958892,  1.959571, 1.960213,             &
     & 1.977193,  1.975540,  1.975023,  1.975015, 1.975346,             &
     & 1.975903,  1.976547,  1.977225,  1.977838, 1.978401,             &
     & 1.978930,  1.979428,  1.979879,  1.980302, 1.980686,             &
     & 1.981060,  1.981401,  1.981722,  1.982025, 1.982312/
!
! *** END OF BLOCK DATA AERSR ******************************************
!
       ENDBLOCKDATA AERSR


!=======================================================================
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCHA
! *** CALCULATES CHLORIDES SPECIATION
!
!     HYDROCHLORIC ACID IN THE LIQUID PHASE IS ASSUMED A MINOR SPECIES,
!     AND DOES NOT SIGNIFICANTLY PERTURB THE HSO4-SO4 EQUILIBRIUM. THE
!     HYDROCHLORIC ACID DISSOLVED IS CALCULATED FROM THE
!     HCL(G) <-> (H+) + (CL-)
!     EQUILIBRIUM, USING THE (H+) FROM THE SULFATES.
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!=======================================================================
!
      SUBROUTINE CALCHA
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
      IMPLICIT NONE
      INCLUDE 'isrpia.inc'
      REAL(KIND=8) :: ALFA,DET,DELT,DIAK,X
      REAL(KIND=8) :: KAPA
!C      CHARACTER ERRINF*40
!
! *** CALCULATE HCL DISSOLUTION *****************************************
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('CALCHA',0,ZHOOK_HANDLE)
      X    = W(5)
      DELT = 0.0d0
      IF (WATER.GT.TINY) THEN
         KAPA = MOLAL(1)
         ALFA = XK3*R*TEMP*(WATER/GAMA(11))**2.0
         DIAK = SQRT( (KAPA+ALFA)**2.0 + 4.0*ALFA*X)
         DELT = 0.5*(-(KAPA+ALFA) + DIAK)
!C         IF (DELT/KAPA.GT.0.1d0) THEN
!C            WRITE (ERRINF,'(1PE10.3)') DELT/KAPA*100.0
!C            CALL PUSHERR (0033, ERRINF)
!C         ENDIF
      ENDIF
!
! *** CALCULATE HCL SPECIATION IN THE GAS PHASE *************************
!
      GHCL     = MAX(X-DELT, 0.0d0)  ! GAS HCL
!
! *** CALCULATE HCL SPECIATION IN THE LIQUID PHASE **********************
!
      MOLAL(4) = DELT                ! CL-
      MOLAL(1) = MOLAL(1) + DELT     ! H+
!
      IF (LHOOK) CALL DR_HOOK('CALCHA',1,ZHOOK_HANDLE)
!
! *** END OF SUBROUTINE CALCHA ******************************************
!
      ENDSUBROUTINE CALCHA





!=======================================================================
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCHAP
! *** CALCULATES CHLORIDES SPECIATION
!
!     HYDROCHLORIC ACID IN THE LIQUID PHASE IS ASSUMED A MINOR SPECIES,
!     THAT DOES NOT SIGNIFICANTLY PERTURB THE HSO4-SO4 EQUILIBRIUM.
!     THE HYDROCHLORIC ACID DISSOLVED IS CALCULATED FROM THE
!     HCL(G) -> HCL(AQ)   AND  HCL(AQ) ->  (H+) + (CL-)
!     EQUILIBRIA, USING (H+) FROM THE SULFATES.
!
!     THIS IS THE VERSION USED BY THE INVERSE PROBLEM SOVER
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!=======================================================================
!
      SUBROUTINE CALCHAP
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
      IMPLICIT NONE
      INCLUDE 'isrpia.inc'
      REAL(KIND=8) :: ALFA,DELT,DIAK,X
!
! *** IS THERE A LIQUID PHASE? ******************************************
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('CALCHAP',0,ZHOOK_HANDLE)
      IF (WATER.LE.TINY)  THEN
        IF (LHOOK) CALL DR_HOOK('CALCHAP',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
!
! *** CALCULATE HCL SPECIATION IN THE GAS PHASE *************************
!
      CALL CALCCLAQ (MOLAL(4), MOLAL(1), DELT)
      ALFA     = XK3*R*TEMP*(WATER/GAMA(11))**2.0
      GASAQ(3) = DELT
      MOLAL(1) = MOLAL(1) - DELT
      MOLAL(4) = MOLAL(4) - DELT
      GHCL     = MOLAL(1)*MOLAL(4)/ALFA
!
      IF (LHOOK) CALL DR_HOOK('CALCHAP',1,ZHOOK_HANDLE)
!
! *** END OF SUBROUTINE CALCHAP *****************************************
!
      ENDSUBROUTINE CALCHAP


!=======================================================================
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCNA
! *** CALCULATES NITRATES SPECIATION
!
!     NITRIC ACID IN THE LIQUID PHASE IS ASSUMED A MINOR SPECIES, THAT
!     DOES NOT SIGNIFICANTLY PERTURB THE HSO4-SO4 EQUILIBRIUM. THE NITRIC
!     ACID DISSOLVED IS CALCULATED FROM THE HNO3(G) -> (H+) + (NO3-)
!     EQUILIBRIUM, USING THE (H+) FROM THE SULFATES.
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!=======================================================================
!
      SUBROUTINE CALCNA
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
      IMPLICIT NONE
      INCLUDE 'isrpia.inc'
      REAL(KIND=8) :: ALFA,DELT,DIAK,X
      REAL(KIND=8) :: KAPA
!C      CHARACTER ERRINF*40
!
! *** CALCULATE HNO3 DISSOLUTION ****************************************
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('CALCNA',0,ZHOOK_HANDLE)
      X    = W(4)
      DELT = 0.0d0
      IF (WATER.GT.TINY) THEN
         KAPA = MOLAL(1)
         ALFA = XK4*R*TEMP*(WATER/GAMA(10))**2.0
         DIAK = SQRT( (KAPA+ALFA)**2.0 + 4.0*ALFA*X)
         DELT = 0.5*(-(KAPA+ALFA) + DIAK)
!C         IF (DELT/KAPA.GT.0.1d0) THEN
!C            WRITE (ERRINF,'(1PE10.3)') DELT/KAPA*100.0
!C            CALL PUSHERR (0019, ERRINF)    ! WARNING ERROR: NO SOLUTION
!C         ENDIF
      ENDIF
!
! *** CALCULATE HNO3 SPECIATION IN THE GAS PHASE ************************
!
      GHNO3    = MAX(X-DELT, 0.0d0)  ! GAS HNO3
!
! *** CALCULATE HNO3 SPECIATION IN THE LIQUID PHASE *********************
!
      MOLAL(7) = DELT                ! NO3-
      MOLAL(1) = MOLAL(1) + DELT     ! H+
!
      IF (LHOOK) CALL DR_HOOK('CALCNA',1,ZHOOK_HANDLE)
!
! *** END OF SUBROUTINE CALCNA ******************************************
!
      ENDSUBROUTINE CALCNA



!=======================================================================
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCNAP
! *** CALCULATES NITRATES SPECIATION
!
!     NITRIC ACID IN THE LIQUID PHASE IS ASSUMED A MINOR SPECIES, THAT
!     DOES NOT SIGNIFICANTLY PERTURB THE HSO4-SO4 EQUILIBRIUM. THE NITRIC
!     ACID DISSOLVED IS CALCULATED FROM THE HNO3(G) -> HNO3(AQ) AND
!     HNO3(AQ) -> (H+) + (CL-) EQUILIBRIA, USING (H+) FROM THE SULFATES.
!
!     THIS IS THE VERSION USED BY THE INVERSE PROBLEM SOVER
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!=======================================================================
!
      SUBROUTINE CALCNAP
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
      IMPLICIT NONE
      INCLUDE 'isrpia.inc'
      REAL(KIND=8) :: ALFA,DELT
!
! *** IS THERE A LIQUID PHASE? ******************************************
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('CALCNAP',0,ZHOOK_HANDLE)
      IF (WATER.LE.TINY)  THEN
        IF (LHOOK) CALL DR_HOOK('CALCNAP',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
!
! *** CALCULATE HNO3 SPECIATION IN THE GAS PHASE ************************
!
      CALL CALCNIAQ (MOLAL(7), MOLAL(1), DELT)
      ALFA     = XK4*R*TEMP*(WATER/GAMA(10))**2.0
      GASAQ(3) = DELT
      MOLAL(1) = MOLAL(1) - DELT
      MOLAL(7) = MOLAL(7) - DELT
      GHNO3    = MOLAL(1)*MOLAL(7)/ALFA
!
      IF (LHOOK) CALL DR_HOOK('CALCNAP',1,ZHOOK_HANDLE)
!
! *** END OF SUBROUTINE CALCNAP *****************************************
!
      ENDSUBROUTINE CALCNAP


!=======================================================================
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCNH3
! *** CALCULATES AMMONIA IN GAS PHASE
!
!     AMMONIA IN THE GAS PHASE IS ASSUMED A MINOR SPECIES, THAT
!     DOES NOT SIGNIFICANTLY PERTURB THE AEROSOL EQUILIBRIUM.
!     AMMONIA GAS IS CALCULATED FROM THE NH3(g) + (H+)(l) <==> (NH4+)(l)
!     EQUILIBRIUM, USING (H+), (NH4+) FROM THE AEROSOL SOLUTION.
!
!     THIS IS THE VERSION USED BY THE DIRECT PROBLEM
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!=======================================================================
!
      SUBROUTINE CALCNH3
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
      IMPLICIT NONE
      INCLUDE 'isrpia.inc'
      REAL(KIND=8) :: A1,CHI1,CHI2,BB,CC,DIAK,PSI
!
! *** IS THERE A LIQUID PHASE? ******************************************
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('CALCNH3',0,ZHOOK_HANDLE)
      IF (WATER.LE.TINY)  THEN
        IF (LHOOK) CALL DR_HOOK('CALCNH3',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
!
! *** CALCULATE NH3 SUBLIMATION *****************************************
!
      A1   = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
      CHI1 = MOLAL(3)
      CHI2 = MOLAL(1)
!
      BB   =(CHI2 + ONE/A1)          ! a=1; b!=1; c!=1
      CC   =-CHI1/A1
      DIAK = SQRT(BB*BB - 4.D0*CC)   ! Always > 0
      PSI  = 0.5*(-BB + DIAK)        ! One positive root
      PSI  = MAX(TINY, MIN(PSI,CHI1))! Constrict in acceptible range
!
! *** CALCULATE NH3 SPECIATION IN THE GAS PHASE *************************
!
      GNH3     = PSI                 ! GAS HNO3
!
! *** CALCULATE NH3 AFFECT IN THE LIQUID PHASE **************************
!
      MOLAL(3) = CHI1 - PSI          ! NH4+
      MOLAL(1) = CHI2 + PSI          ! H+
!
      IF (LHOOK) CALL DR_HOOK('CALCNH3',1,ZHOOK_HANDLE)
!
! *** END OF SUBROUTINE CALCNH3 *****************************************
!
      ENDSUBROUTINE CALCNH3



!=======================================================================
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCNH3P
! *** CALCULATES AMMONIA IN GAS PHASE
!
!     AMMONIA GAS IS CALCULATED FROM THE NH3(g) + (H+)(l) <==> (NH4+)(l)
!     EQUILIBRIUM, USING (H+), (NH4+) FROM THE AEROSOL SOLUTION.
!
!     THIS IS THE VERSION USED BY THE INVERSE PROBLEM SOLVER
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!=======================================================================
!
      SUBROUTINE CALCNH3P
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
      IMPLICIT NONE
      INCLUDE 'isrpia.inc'
      REAL(KIND=8) :: A1
!
! *** IS THERE A LIQUID PHASE? ******************************************
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('CALCNH3P',0,ZHOOK_HANDLE)
      IF (WATER.LE.TINY)  THEN
        IF (LHOOK) CALL DR_HOOK('CALCNH3P',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
!
! *** CALCULATE NH3 GAS PHASE CONCENTRATION *****************************
!
      A1   = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
      GNH3 = MOLAL(3)/MOLAL(1)/A1
!
      IF (LHOOK) CALL DR_HOOK('CALCNH3P',1,ZHOOK_HANDLE)
!
! *** END OF SUBROUTINE CALCNH3P ****************************************
!
      ENDSUBROUTINE CALCNH3P


!=======================================================================
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCNHA
!
!     THIS SUBROUTINE CALCULATES THE DISSOLUTION OF HCL, HNO3 AT
!     THE PRESENCE OF (H,SO4). HCL, HNO3 ARE CONSIDERED MINOR SPECIES,
!     THAT DO NOT SIGNIFICANTLY AFFECT THE EQUILIBRIUM POINT.
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!=======================================================================
!
      SUBROUTINE CALCNHA
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
      IMPLICIT NONE
      INCLUDE 'isrpia.inc'
      REAL(KIND=8) :: A3,A4,C1,C2,C3,CHI3,CHI4,DELCL,DELNO,OMEGA
      INTEGER :: ISLV
      REAL(KIND=8) :: M1, M2, M3
      CHARACTER ERRINF*40
!
! *** SPECIAL CASE; WATER=ZERO ******************************************
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('CALCNHA',0,ZHOOK_HANDLE)
      IF (WATER.LE.TINY) THEN
         GOTO 55
!
! *** SPECIAL CASE; HCL=HNO3=ZERO ***************************************
!
      ELSEIF (W(5).LE.TINY .AND. W(4).LE.TINY) THEN
         GOTO 60
!
! *** SPECIAL CASE; HCL=ZERO ********************************************
!
      ELSE IF (W(5).LE.TINY) THEN
         CALL CALCNA              ! CALL HNO3 DISSOLUTION ROUTINE
         GOTO 60
!
! *** SPECIAL CASE; HNO3=ZERO *******************************************
!
      ELSE IF (W(4).LE.TINY) THEN
         CALL CALCHA              ! CALL HCL DISSOLUTION ROUTINE
         GOTO 60
      ENDIF
!
! *** CALCULATE EQUILIBRIUM CONSTANTS ***********************************
!
      A3 = XK4*R*TEMP*(WATER/GAMA(10))**2.0   ! HNO3
      A4 = XK3*R*TEMP*(WATER/GAMA(11))**2.0   ! HCL
!
! *** CALCULATE CUBIC EQUATION COEFFICIENTS *****************************
!
      DELCL = ZERO
      DELNO = ZERO
!
      OMEGA = MOLAL(1)       ! H+
      CHI3  = W(4)           ! HNO3
      CHI4  = W(5)           ! HCL
!
      C1    = A3*CHI3
      C2    = A4*CHI4
      C3    = A3 - A4
!
      M1    = (C1 + C2 + (OMEGA+A4)*C3)/C3
      M2    = ((OMEGA+A4)*C2 - A4*C3*CHI4)/C3
      M3    =-A4*C2*CHI4/C3
!
! *** CALCULATE ROOTS ***************************************************
!
      CALL POLY3 (M1, M2, M3, DELCL, ISLV) ! HCL DISSOLUTION
      IF (ISLV.NE.0) THEN
         DELCL = TINY       ! TINY AMOUNTS OF HCL ASSUMED WHEN NO ROOT
         WRITE (ERRINF,'(1PE7.1)') TINY
         CALL PUSHERR (0022, ERRINF)    ! WARNING ERROR: NO SOLUTION
      ENDIF
      DELCL = MIN(DELCL, CHI4)
!
      DELNO = C1*DELCL/(C2 + C3*DELCL)
      DELNO = MIN(DELNO, CHI3)
!
      IF (DELCL.LT.ZERO .OR. DELNO.LT.ZERO .OR.                         &
     &   DELCL.GT.CHI4 .OR. DELNO.GT.CHI3       ) THEN
         DELCL = TINY  ! TINY AMOUNTS OF HCL ASSUMED WHEN NO ROOT
         DELNO = TINY
         WRITE (ERRINF,'(1PE7.1)') TINY
         CALL PUSHERR (0022, ERRINF)    ! WARNING ERROR: NO SOLUTION
      ENDIF
!CC
!CC *** COMPARE DELTA TO TOTAL H+ ; ESTIMATE EFFECT TO HSO4 ***************
!CC
!C      IF ((DELCL+DELNO)/MOLAL(1).GT.0.1d0) THEN
!C         WRITE (ERRINF,'(1PE10.3)') (DELCL+DELNO)/MOLAL(1)*100.0
!C         CALL PUSHERR (0021, ERRINF)
!C      ENDIF
!
! *** EFFECT ON LIQUID PHASE ********************************************
!
50    MOLAL(1) = MOLAL(1) + (DELNO+DELCL)  ! H+   CHANGE
      MOLAL(4) = MOLAL(4) + DELCL          ! CL-  CHANGE
      MOLAL(7) = MOLAL(7) + DELNO          ! NO3- CHANGE
!
! *** EFFECT ON GAS PHASE ***********************************************
!
55    GHCL     = MAX(W(5) - MOLAL(4), TINY)
      GHNO3    = MAX(W(4) - MOLAL(7), TINY)
!
60    CONTINUE
      IF (LHOOK) CALL DR_HOOK('CALCNHA',1,ZHOOK_HANDLE)
!
! *** END OF SUBROUTINE CALCNHA *****************************************
!
      ENDSUBROUTINE CALCNHA



!=======================================================================
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCNHP
!
!     THIS SUBROUTINE CALCULATES THE GAS PHASE NITRIC AND HYDROCHLORIC
!     ACID. CONCENTRATIONS ARE CALCULATED FROM THE DISSOLUTION
!     EQUILIBRIA, USING (H+), (Cl-), (NO3-) IN THE AEROSOL PHASE.
!
!     THIS IS THE VERSION USED BY THE INVERSE PROBLEM SOLVER
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!=======================================================================
!
      SUBROUTINE CALCNHP
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
      IMPLICIT NONE
      INCLUDE 'isrpia.inc'
      REAL(KIND=8) :: A3,A4,DELT
!
! *** IS THERE A LIQUID PHASE? ******************************************
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('CALCNHP',0,ZHOOK_HANDLE)
      IF (WATER.LE.TINY)  THEN
        IF (LHOOK) CALL DR_HOOK('CALCNHP',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
!
! *** CALCULATE EQUILIBRIUM CONSTANTS ***********************************
!
      A3       = XK3*R*TEMP*(WATER/GAMA(11))**2.0
      A4       = XK4*R*TEMP*(WATER/GAMA(10))**2.0
      MOLAL(1) = MOLAL(1) + WAER(4) + WAER(5)
!
! *** CALCULATE CONCENTRATIONS ******************************************
! *** ASSUME THAT 'DELT' FROM HNO3 >> 'DELT' FROM HCL
!
      CALL CALCNIAQ (WAER(4), MOLAL(1)+MOLAL(7)+MOLAL(4), DELT)
      MOLAL(1) = MOLAL(1) - DELT
      MOLAL(7) = WAER(4)  - DELT
      GASAQ(3) = DELT
!
      CALL CALCCLAQ (WAER(5), MOLAL(1)+MOLAL(7)+MOLAL(4), DELT)
      MOLAL(1) = MOLAL(1) - DELT
      MOLAL(4) = WAER(5)  - DELT
      GASAQ(2) = DELT
!
      GHNO3    = MOLAL(1)*MOLAL(7)/A4
      GHCL     = MOLAL(1)*MOLAL(4)/A3
!
      IF (LHOOK) CALL DR_HOOK('CALCNHP',1,ZHOOK_HANDLE)
!
! *** END OF SUBROUTINE CALCNHP *****************************************
!
      ENDSUBROUTINE CALCNHP

!=======================================================================
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCAMAQ
! *** THIS SUBROUTINE CALCULATES THE NH3(aq) GENERATED FROM (H,NH4+).
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!=======================================================================
!
      SUBROUTINE CALCAMAQ (NH4I, OHI, DELT)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
      IMPLICIT NONE
      INCLUDE 'isrpia.inc'
      REAL(KIND=8) :: A22,AKW,BB,CC,DD,OM1,OM2,DEL1,DEL2
      REAL(KIND=8) :: OHI,DELT
      REAL(KIND=8) :: NH4I
!C      CHARACTER ERRINF*40
!
! *** EQUILIBRIUM CONSTANTS
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('CALCAMAQ',0,ZHOOK_HANDLE)
      A22  = XK22/XKW/WATER*(GAMA(8)/GAMA(9))**2. ! GAMA(NH3) ASSUMED 1
      AKW  = XKW *RH*WATER*WATER
!
! *** FIND ROOT
!
      OM1  = NH4I
      OM2  = OHI
      BB   =-(OM1+OM2+A22*AKW)
      CC   = OM1*OM2
      DD   = SQRT(BB*BB-4.D0*CC)

      DEL1 = 0.5D0*(-BB - DD)
      DEL2 = 0.5D0*(-BB + DD)
!
! *** GET APPROPRIATE ROOT.
!
      IF (DEL1.LT.ZERO) THEN
         IF (DEL2.GT.NH4I .OR. DEL2.GT.OHI) THEN
            DELT = ZERO
         ELSE
            DELT = DEL2
         ENDIF
      ELSE
         DELT = DEL1
      ENDIF
!C
!C *** COMPARE DELTA TO TOTAL NH4+ ; ESTIMATE EFFECT *********************
!C
!C      IF (DELTA/HYD.GT.0.1d0) THEN
!C         WRITE (ERRINF,'(1PE10.3)') DELTA/HYD*100.0
!C         CALL PUSHERR (0020, ERRINF)
!C      ENDIF
!
      IF (LHOOK) CALL DR_HOOK('CALCAMAQ',1,ZHOOK_HANDLE)
!
! *** END OF SUBROUTINE CALCAMAQ ****************************************
!
      ENDSUBROUTINE CALCAMAQ



!=======================================================================
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCAMAQ2
!
!     THIS SUBROUTINE CALCULATES THE NH3(aq) GENERATED FROM (H,NH4+).
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!=======================================================================
!
      SUBROUTINE CALCAMAQ2 (GGNH3, NH4I, OHI, NH3AQ)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
      IMPLICIT NONE
      INCLUDE 'isrpia.inc'
      REAL(KIND=8) :: A22,ALF1,ALF2,AKW,BB,CC,DEL
      REAL(KIND=8) :: GGNH3,OHI
      REAL(KIND=8) :: NH4I, NH3AQ
!
! *** EQUILIBRIUM CONSTANTS
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('CALCAMAQ2',0,ZHOOK_HANDLE)
      A22  = XK22/XKW/WATER*(GAMA(8)/GAMA(9))**2. ! GAMA(NH3) ASSUMED 1
      AKW  = XKW *RH*WATER*WATER
!
! *** FIND ROOT
!
      ALF1 = NH4I - GGNH3
      ALF2 = GGNH3
      BB   = ALF1 + A22*AKW
      CC   =-A22*AKW*ALF2
      DEL  = 0.5D0*(-BB + SQRT(BB*BB-4.D0*CC))
!
! *** ADJUST CONCENTRATIONS
!
      NH4I  = ALF1 + DEL
      OHI   = DEL
      IF (OHI.LE.TINY) OHI = SQRT(AKW)   ! If solution is neutral.
      NH3AQ = ALF2 - DEL
!
      IF (LHOOK) CALL DR_HOOK('CALCAMAQ2',1,ZHOOK_HANDLE)
!
! *** END OF SUBROUTINE CALCAMAQ2 ****************************************
!
      ENDSUBROUTINE CALCAMAQ2



!=======================================================================
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCCLAQ
!
!     THIS SUBROUTINE CALCULATES THE HCL(aq) GENERATED FROM (H+,CL-).
!
! *** COPYRIGHT 1996-98, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!=======================================================================
!
      SUBROUTINE CALCCLAQ (CLI, HI, DELT)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
      IMPLICIT NONE
      INCLUDE 'isrpia.inc'
      REAL(KIND=8) :: HI,DELT
      REAL(KIND=8) :: A32,BB,CC,DD,DEL1,DEL2,OM1,OM2,X
      REAL(KIND=8) :: CLI
!
! *** EQUILIBRIUM CONSTANTS
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('CALCCLAQ',0,ZHOOK_HANDLE)
      A32  = XK32*WATER/(GAMA(11))**2. ! GAMA(HCL) ASSUMED 1
!
! *** FIND ROOT
!
      OM1  = CLI
      OM2  = HI
      BB   =-(OM1+OM2+A32)
      CC   = OM1*OM2
      DD   = SQRT(BB*BB-4.D0*CC)

      DEL1 = 0.5D0*(-BB - DD)
      DEL2 = 0.5D0*(-BB + DD)
!
! *** GET APPROPRIATE ROOT.
!
      IF (DEL1.LT.ZERO) THEN
         IF (DEL2.LT.ZERO .OR. DEL2.GT.CLI .OR. DEL2.GT.HI) THEN
            DELT = ZERO
         ELSE
            DELT = DEL2
         ENDIF
      ELSE
         DELT = DEL1
      ENDIF
!
      IF (LHOOK) CALL DR_HOOK('CALCCLAQ',1,ZHOOK_HANDLE)
!
! *** END OF SUBROUTINE CALCCLAQ ****************************************
!
      ENDSUBROUTINE CALCCLAQ



!=======================================================================
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCCLAQ2
!
!     THIS SUBROUTINE CALCULATES THE HCL(aq) GENERATED FROM (H+,CL-).
!
! *** COPYRIGHT 1996-98, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!=======================================================================
!
      SUBROUTINE CALCCLAQ2 (GGCL, CLI, HI, CLAQ)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
      IMPLICIT NONE
      INCLUDE 'isrpia.inc'
      REAL(KIND=8) :: A32,AKW,ALF1,ALF2,COEF,DEL1
      REAL(KIND=8) :: GGCL,HI,CLAQ
      REAL(KIND=8) :: CLI
!
! *** EQUILIBRIUM CONSTANTS
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('CALCCLAQ2',0,ZHOOK_HANDLE)
      A32  = XK32*WATER/(GAMA(11))**2. ! GAMA(HCL) ASSUMED 1
      AKW  = XKW *RH*WATER*WATER
!
! *** FIND ROOT
!
      ALF1  = CLI - GGCL
      ALF2  = GGCL
      COEF  = (ALF1+A32)
      DEL1  = 0.5*(-COEF + SQRT(COEF*COEF+4.D0*A32*ALF2))
!
! *** CORRECT CONCENTRATIONS
!
      CLI  = ALF1 + DEL1
      HI   = DEL1
      IF (HI.LE.TINY) HI = SQRT(AKW)   ! If solution is neutral.
      CLAQ = ALF2 - DEL1
!
      IF (LHOOK) CALL DR_HOOK('CALCCLAQ2',1,ZHOOK_HANDLE)
!
! *** END OF SUBROUTINE CALCCLAQ2 ****************************************
!
      ENDSUBROUTINE CALCCLAQ2



!=======================================================================
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCNIAQ
!
!     THIS SUBROUTINE CALCULATES THE HNO3(aq) GENERATED FROM (H,NO3-).
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!=======================================================================
!
      SUBROUTINE CALCNIAQ (NO3I, HI, DELT)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
      IMPLICIT NONE
      INCLUDE 'isrpia.inc'
      REAL(KIND=8) :: HI,DELT
      REAL(KIND=8) :: A42,BB,CC,DD,OM1,OM2,DEL1,DEL2
      REAL(KIND=8) :: NO3I
!
! *** EQUILIBRIUM CONSTANTS
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('CALCNIAQ',0,ZHOOK_HANDLE)
      A42  = XK42*WATER/(GAMA(10))**2. ! GAMA(HNO3) ASSUMED 1
!
! *** FIND ROOT
!
      OM1  = NO3I
      OM2  = HI
      BB   =-(OM1+OM2+A42)
      CC   = OM1*OM2
      DD   = SQRT(BB*BB-4.D0*CC)

      DEL1 = 0.5D0*(-BB - DD)
      DEL2 = 0.5D0*(-BB + DD)
!
! *** GET APPROPRIATE ROOT.
!
      IF (DEL1.LT.ZERO .OR. DEL1.GT.HI .OR. DEL1.GT.NO3I) THEN
         DELT = ZERO
      ELSE
         DELT = DEL1
         IF (LHOOK) CALL DR_HOOK('CALCNIAQ',1,ZHOOK_HANDLE)
         RETURN
      ENDIF
!
      IF (DEL2.LT.ZERO .OR. DEL2.GT.NO3I .OR. DEL2.GT.HI) THEN
         DELT = ZERO
      ELSE
         DELT = DEL2
      ENDIF
!
      IF (LHOOK) CALL DR_HOOK('CALCNIAQ',1,ZHOOK_HANDLE)
!
! *** END OF SUBROUTINE CALCNIAQ ****************************************
!
      ENDSUBROUTINE CALCNIAQ



!=======================================================================
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCNIAQ2
!
!     THIS SUBROUTINE CALCULATES THE UNDISSOCIATED HNO3(aq)
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!=======================================================================
!
      SUBROUTINE CALCNIAQ2 (GGNO3, NO3I, HI, NO3AQ)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
      IMPLICIT NONE
      INCLUDE 'isrpia.inc'
      REAL(KIND=8) :: A42,AKW,ALF1,ALF2,ALF3,BB,CC,DEL1
      REAL(KIND=8) :: GGNO3,HI
      REAL(KIND=8) :: NO3I, NO3AQ
!
! *** EQUILIBRIUM CONSTANTS
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('CALCNIAQ2',0,ZHOOK_HANDLE)
      A42  = XK42*WATER/(GAMA(10))**2. ! GAMA(HNO3) ASSUMED 1
      AKW  = XKW *RH*WATER*WATER
!
! *** FIND ROOT
!
      ALF1  = NO3I - GGNO3
      ALF2  = GGNO3
      ALF3  = HI
!
      BB    = ALF3 + ALF1 + A42
      CC    = ALF3*ALF1 - A42*ALF2
      DEL1  = 0.5*(-BB + SQRT(BB*BB-4.D0*CC))
!
! *** CORRECT CONCENTRATIONS
!
      NO3I  = ALF1 + DEL1
      HI    = ALF3 + DEL1
      IF (HI.LE.TINY) HI = SQRT(AKW)   ! If solution is neutral.
      NO3AQ = ALF2 - DEL1
!
      IF (LHOOK) CALL DR_HOOK('CALCNIAQ2',1,ZHOOK_HANDLE)
!
! *** END OF SUBROUTINE CALCNIAQ2 ****************************************
!
      ENDSUBROUTINE CALCNIAQ2


!=======================================================================
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCMR
! *** THIS SUBROUTINE CALCULATES:
!     1. ION PAIR CONCENTRATIONS (FROM [MOLAR] ARRAY)
!     2. WATER CONTENT OF LIQUID AEROSOL PHASE (FROM ZSR CORRELATION)
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!=======================================================================
!
      SUBROUTINE CALCMR
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
      IMPLICIT NONE
      INCLUDE 'isrpia.inc'
      REAL(KIND=8) :: AML5,FRCL,FRNH4,FRNO3,HSO4I,SO4I,TOTS4
      INTEGER :: I
      REAL(KIND=8) :: CHI1,CHI2,CHI3,CHI4,CHI5,CHI6,CHI7,CHI8
      REAL(KIND=8) :: PSI1,PSI2,PSI3,PSI4,PSI5,PSI6,PSI7,PSI8
      REAL(KIND=8) :: A1,A2,A3,A4,A5,A6,A7,A8
      COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,    &
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,    &
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8
      CHARACTER SC*1
!
! *** CALCULATE ION PAIR CONCENTRATIONS ACCORDING TO SPECIFIC CASE ****
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('CALCMR',0,ZHOOK_HANDLE)
      SC =SCASE(1:1)                   ! SULRAT & SODRAT case
!
! *** NH4-SO4 SYSTEM ; SULFATE POOR CASE
!
      IF (SC.EQ.'A') THEN
         MOLALR(4) = MOLAL(5)+MOLAL(6) ! (NH4)2SO4 - CORRECT FOR SO4 TO HSO4
!
! *** NH4-SO4 SYSTEM ; SULFATE RICH CASE ; NO FREE ACID
!
      ELSE IF (SC.EQ.'B') THEN
         SO4I  = MOLAL(5)-MOLAL(1)     ! CORRECT FOR HSO4 DISSOCIATION
         HSO4I = MOLAL(6)+MOLAL(1)
         IF (SO4I.LT.HSO4I) THEN
            MOLALR(13) = SO4I                   ! [LC] = [SO4]
            MOLALR(9)  = MAX(HSO4I-SO4I, ZERO)  ! NH4HSO4
         ELSE
            MOLALR(13) = HSO4I                  ! [LC] = [HSO4]
            MOLALR(4)  = MAX(SO4I-HSO4I, ZERO)  ! (NH4)2SO4
         ENDIF
!
! *** NH4-SO4 SYSTEM ; SULFATE RICH CASE ; FREE ACID
!
      ELSE IF (SC.EQ.'C') THEN
         MOLALR(4) = MOLAL(3)                     ! NH4HSO4
         MOLALR(7) = MAX(W(2)-W(3), ZERO)         ! H2SO4
!
! *** NH4-SO4-NO3 SYSTEM ; SULFATE POOR CASE
!
      ELSE IF (SC.EQ.'D') THEN
         MOLALR(4) = MOLAL(5) + MOLAL(6)          ! (NH4)2SO4
         AML5      = MOLAL(3)-2.D0*MOLALR(4)      ! "free" NH4
         MOLALR(5) = MAX(MIN(AML5,MOLAL(7)), ZERO)! NH4NO3 = MIN("free", NO3)
!
! *** NH4-SO4-NO3 SYSTEM ; SULFATE RICH CASE ; NO FREE ACID
!
      ELSE IF (SC.EQ.'E') THEN
         SO4I  = MAX(MOLAL(5)-MOLAL(1),ZERO)      ! FROM HSO4 DISSOCIATION
         HSO4I = MOLAL(6)+MOLAL(1)
         IF (SO4I.LT.HSO4I) THEN
            MOLALR(13) = SO4I                     ! [LC] = [SO4]
            MOLALR(9)  = MAX(HSO4I-SO4I, ZERO)    ! NH4HSO4
         ELSE
            MOLALR(13) = HSO4I                    ! [LC] = [HSO4]
            MOLALR(4)  = MAX(SO4I-HSO4I, ZERO)    ! (NH4)2SO4
         ENDIF
!
! *** NH4-SO4-NO3 SYSTEM ; SULFATE RICH CASE ; FREE ACID
!
      ELSE IF (SC.EQ.'F') THEN
         MOLALR(4) = MOLAL(3)                              ! NH4HSO4
         MOLALR(7) = MAX(MOLAL(5)+MOLAL(6)-MOLAL(3),ZERO)  ! H2SO4
!
! *** NA-NH4-SO4-NO3-CL SYSTEM ; SULFATE POOR ; SODIUM POOR CASE
!
      ELSE IF (SC.EQ.'G') THEN
         MOLALR(2) = 0.5*MOLAL(2)                          ! NA2SO4
         TOTS4     = MOLAL(5)+MOLAL(6)                     ! Total SO4
         MOLALR(4) = MAX(TOTS4 - MOLALR(2), ZERO)          ! (NH4)2SO4
         FRNH4     = MAX(MOLAL(3) - 2.D0*MOLALR(4), ZERO)
         MOLALR(5) = MIN(MOLAL(7),FRNH4)                   ! NH4NO3
         FRNH4     = MAX(FRNH4 - MOLALR(5), ZERO)
         MOLALR(6) = MIN(MOLAL(4), FRNH4)                  ! NH4CL
!
! *** NA-NH4-SO4-NO3-CL SYSTEM ; SULFATE POOR ; SODIUM RICH CASE
! *** RETREIVE DISSOLVED SALTS DIRECTLY FROM COMMON BLOCK /SOLUT/
!
      ELSE IF (SC.EQ.'H') THEN
         MOLALR(1) = PSI7                                  ! NACL
         MOLALR(2) = PSI1                                  ! NA2SO4
         MOLALR(3) = PSI8                                  ! NANO3
         MOLALR(4) = ZERO                                  ! (NH4)2SO4
         FRNO3     = MAX(MOLAL(7) - MOLALR(3), ZERO)       ! "FREE" NO3
         FRCL      = MAX(MOLAL(4) - MOLALR(1), ZERO)       ! "FREE" CL
         MOLALR(5) = MIN(MOLAL(3),FRNO3)                   ! NH4NO3
         FRNH4     = MAX(MOLAL(3) - MOLALR(5), ZERO)       ! "FREE" NH3
         MOLALR(6) = MIN(FRCL, FRNH4)                      ! NH4CL
!
! *** NA-NH4-SO4-NO3-CL SYSTEM ; SULFATE RICH CASE ; NO FREE ACID
! *** RETREIVE DISSOLVED SALTS DIRECTLY FROM COMMON BLOCK /SOLUT/
!
      ELSE IF (SC.EQ.'I') THEN
         MOLALR(04) = PSI5                                 ! (NH4)2SO4
         MOLALR(02) = PSI4                                 ! NA2SO4
         MOLALR(09) = PSI1                                 ! NH4HSO4
         MOLALR(12) = PSI3                                 ! NAHSO4
         MOLALR(13) = PSI2                                 ! LC
!
! *** NA-NH4-SO4-NO3-CL SYSTEM ; SULFATE RICH CASE ; FREE ACID
!
      ELSE IF (SC.EQ.'J') THEN
         MOLALR(09) = MOLAL(3)                             ! NH4HSO4
         MOLALR(12) = MOLAL(2)                             ! NAHSO4
         MOLALR(07) = MOLAL(5)+MOLAL(6)-MOLAL(3)-MOLAL(2)  ! H2SO4
         MOLALR(07) = MAX(MOLALR(07),ZERO)
!
! ======= REVERSE PROBLEMS ===========================================
!
! *** NH4-SO4-NO3 SYSTEM ; SULFATE POOR CASE
!
      ELSE IF (SC.EQ.'N') THEN
         MOLALR(4) = MOLAL(5) + MOLAL(6)          ! (NH4)2SO4
         AML5      = WAER(3)-2.D0*MOLALR(4)       ! "free" NH4
         MOLALR(5) = MAX(MIN(AML5,WAER(4)), ZERO) ! NH4NO3 = MIN("free", NO3)
!
! *** NH4-SO4-NO3-NA-CL SYSTEM ; SULFATE POOR, SODIUM POOR CASE
!
      ELSE IF (SC.EQ.'Q') THEN
         MOLALR(2) = PSI1                                  ! NA2SO4
         MOLALR(4) = PSI6                                  ! (NH4)2SO4
         MOLALR(5) = PSI5                                  ! NH4NO3
         MOLALR(6) = PSI4                                  ! NH4CL
!
! *** NH4-SO4-NO3-NA-CL SYSTEM ; SULFATE POOR, SODIUM RICH CASE
!
      ELSE IF (SC.EQ.'R') THEN
         MOLALR(1) = PSI3                                  ! NACL
         MOLALR(2) = PSI1                                  ! NA2SO4
         MOLALR(3) = PSI2                                  ! NANO3
         MOLALR(4) = ZERO                                  ! (NH4)2SO4
         MOLALR(5) = PSI5                                  ! NH4NO3
         MOLALR(6) = PSI4                                  ! NH4CL
!
! *** UNKNOWN CASE
!
      ELSE
         CALL PUSHERR (1001, ' ') ! FATAL ERROR: CASE NOT SUPPORTED
      ENDIF
!
! *** CALCULATE WATER CONTENT ; ZSR CORRELATION ***********************
!
      WATER = ZERO
      DO 10 I=1,NPAIR
         WATER = WATER + MOLALR(I)/M0(I)
10    CONTINUE
      WATER = MAX(WATER, TINY)
!
      IF (LHOOK) CALL DR_HOOK('CALCMR',1,ZHOOK_HANDLE)
!
! *** END OF SUBROUTINE CALCMR ******************************************
!
      ENDSUBROUTINE CALCMR
!=======================================================================
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCMDRH
!
!     THIS IS THE CASE WHERE THE RELATIVE HUMIDITY IS IN THE MUTUAL
!     DRH REGION. THE SOLUTION IS ASSUMED TO BE THE SUM OF TWO WEIGHTED
!     SOLUTIONS ; THE 'DRY' SOLUTION (SUBROUTINE DRYCASE) AND THE
!     'SATURATED LIQUID' SOLUTION (SUBROUTINE LIQCASE).
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!=======================================================================
!
      SUBROUTINE CALCMDRH (RHI, RHDRY, RHLIQ, DRYCASE, LIQCASE)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
      IMPLICIT NONE
      INCLUDE 'isrpia.inc'
      REAL(KIND=8) :: CLCO,CNA2SO,CNACLO,CNAHSO,CNANO
      REAL(KIND=8) :: CNH42SO,CNH4N3O,CNH4CLO,CNH4HSO
      REAL(KIND=8) :: DAMBIS,DAMCHL,DAMG,DAMNIT,DAMSUL,DHAG,DLC,DNAG
      REAL(KIND=8) :: DSOBIS,DSONIT,DSOSUL,DSOCHL
      REAL(KIND=8) :: GHCLO,GHO3O,GNH3O,GHNO3O,ONEMWF,WF
      INTEGER :: I
      REAL(KIND=8) :: RHI,RHDRY,RHLIQ
      EXTERNAL DRYCASE, LIQCASE
!
! *** FIND WEIGHT FACTOR **********************************************
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('CALCMDRH',0,ZHOOK_HANDLE)
      IF (WFTYP.EQ.0) THEN
         WF = ONE
      ELSEIF (WFTYP.EQ.1) THEN
         WF = 0.5D0
      ELSE
         WF = (RHLIQ-RHI)/(RHLIQ-RHDRY)
      ENDIF
      ONEMWF  = ONE - WF
!
! *** FIND FIRST SECTION ; DRY ONE ************************************
!
      CALL DRYCASE
      IF (ABS(ONEMWF).LE.1D-5) GOTO 200  ! DRY AEROSOL
!
      CNH42SO = CNH42S4                  ! FIRST (DRY) SOLUTION
      CNH4HSO = CNH4HS4
      CLCO    = CLC
      CNH4N3O = CNH4NO3
      CNH4CLO = CNH4CL
      CNA2SO  = CNA2SO4
      CNAHSO  = CNAHSO4
      CNANO   = CNANO3
      CNACLO  = CNACL
      GNH3O   = GNH3
      GHNO3O  = GHNO3
      GHCLO   = GHCL
!
! *** FIND SECOND SECTION ; DRY & LIQUID ******************************
!
      CNH42S4 = ZERO
      CNH4HS4 = ZERO
      CLC     = ZERO
      CNH4NO3 = ZERO
      CNH4CL  = ZERO
      CNA2SO4 = ZERO
      CNAHSO4 = ZERO
      CNANO3  = ZERO
      CNACL   = ZERO
      GNH3    = ZERO
      GHNO3   = ZERO
      GHCL    = ZERO
      CALL LIQCASE                   ! SECOND (LIQUID) SOLUTION
!
! *** ADJUST THINGS FOR THE CASE THAT THE LIQUID SUB PREDICTS DRY AEROSOL
!
      IF (WATER.LE.TINY) THEN
         DO 100 I=1,NIONS
            MOLAL(I)= ZERO           ! Aqueous phase
  100    CONTINUE
         WATER   = ZERO
!
         CNH42S4 = CNH42SO           ! Solid phase
         CNA2SO4 = CNA2SO
         CNAHSO4 = CNAHSO
         CNH4HS4 = CNH4HSO
         CLC     = CLCO
         CNH4NO3 = CNH4N3O
         CNANO3  = CNANO
         CNACL   = CNACLO
         CNH4CL  = CNH4CLO
!
         GNH3    = GNH3O             ! Gas phase
         GHNO3   = GHNO3O
         GHCL    = GHCLO
!
         GOTO 200
      ENDIF
!
! *** FIND SALT DISSOLUTIONS BETWEEN DRY & LIQUID SOLUTIONS.
!
      DAMSUL  = CNH42SO - CNH42S4
      DSOSUL  = CNA2SO  - CNA2SO4
      DAMBIS  = CNH4HSO - CNH4HS4
      DSOBIS  = CNAHSO  - CNAHSO4
      DLC     = CLCO    - CLC
      DAMNIT  = CNH4N3O - CNH4NO3
      DAMCHL  = CNH4CLO - CNH4CL
      DSONIT  = CNANO   - CNANO3
      DSOCHL  = CNACLO  - CNACL
!
! *** FIND GAS DISSOLUTIONS BETWEEN DRY & LIQUID SOLUTIONS.
!
      DAMG    = GNH3O   - GNH3
      DHAG    = GHCLO   - GHCL
      DNAG    = GHNO3O  - GHNO3
!
! *** FIND SOLUTION AT MDRH BY WEIGHTING DRY & LIQUID SOLUTIONS.
!
!     LIQUID
!
      MOLAL(1)= ONEMWF*MOLAL(1)                                 ! H+
      MOLAL(2)= ONEMWF*(2.D0*DSOSUL + DSOBIS + DSONIT + DSOCHL) ! NA+
      MOLAL(3)= ONEMWF*(2.D0*DAMSUL + DAMG   + DAMBIS + DAMCHL +        &
     &                  3.D0*DLC    + DAMNIT )                  ! NH4+
      MOLAL(4)= ONEMWF*(     DAMCHL + DSOCHL + DHAG)            ! CL-
      MOLAL(5)= ONEMWF*(     DAMSUL + DSOSUL + DLC - MOLAL(6))  ! SO4-- !VB 17 Sept 2001
      MOLAL(6)= ONEMWF*(   MOLAL(6) + DSOBIS + DAMBIS + DLC)    ! HSO4-
      MOLAL(7)= ONEMWF*(     DAMNIT + DSONIT + DNAG)            ! NO3-
      WATER   = ONEMWF*WATER
!
!     SOLID
!
      CNH42S4 = WF*CNH42SO + ONEMWF*CNH42S4
      CNA2SO4 = WF*CNA2SO  + ONEMWF*CNA2SO4
      CNAHSO4 = WF*CNAHSO  + ONEMWF*CNAHSO4
      CNH4HS4 = WF*CNH4HSO + ONEMWF*CNH4HS4
      CLC     = WF*CLCO    + ONEMWF*CLC
      CNH4NO3 = WF*CNH4N3O + ONEMWF*CNH4NO3
      CNANO3  = WF*CNANO   + ONEMWF*CNANO3
      CNACL   = WF*CNACLO  + ONEMWF*CNACL
      CNH4CL  = WF*CNH4CLO + ONEMWF*CNH4CL
!
!     GAS
!
      GNH3    = WF*GNH3O   + ONEMWF*GNH3
      GHNO3   = WF*GHNO3O  + ONEMWF*GHNO3
      GHCL    = WF*GHCLO   + ONEMWF*GHCL
!
! *** RETURN POINT
!
200   CONTINUE
      IF (LHOOK) CALL DR_HOOK('CALCMDRH',1,ZHOOK_HANDLE)
!
! *** END OF SUBROUTINE CALCMDRH ****************************************
!
      ENDSUBROUTINE CALCMDRH






!=======================================================================
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCMDRP
!
!     THIS IS THE CASE WHERE THE RELATIVE HUMIDITY IS IN THE MUTUAL
!     DRH REGION. THE SOLUTION IS ASSUMED TO BE THE SUM OF TWO WEIGHTED
!     SOLUTIONS ; THE 'DRY' SOLUTION (SUBROUTINE DRYCASE) AND THE
!     'SATURATED LIQUID' SOLUTION (SUBROUTINE LIQCASE).
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!=======================================================================
!
      SUBROUTINE CALCMDRP (RHI, RHDRY, RHLIQ, DRYCASE, LIQCASE)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
      INCLUDE 'isrpia.inc'
      REAL(KIND=8) :: A2,A3,A4,A8,CLC0,CNA2SO,CNACLO,CNAHSO,CNANO
      REAL(KIND=8) :: CNH4HSO,CNH4N30,CNH42SO,CNH4CLO,CNH4SO,CNH4N3O
      REAL(KIND=8) :: CLCO
      REAL(KIND=8) :: DAMBIS,DLC,DSOBIS,HIEN,HIEQ
      INTEGER :: I
      REAL(KIND=8) :: WF,ONEMWF
      REAL(KIND=8) :: RHI,RHDRY,RHLIQ
      EXTERNAL DRYCASE, LIQCASE
!
! *** FIND WEIGHT FACTOR **********************************************
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('CALCMDRP',0,ZHOOK_HANDLE)
      IF (WFTYP.EQ.0) THEN
         WF = ONE
      ELSEIF (WFTYP.EQ.1) THEN
         WF = 0.5D0
      ELSE
         WF = (RHLIQ-RHI)/(RHLIQ-RHDRY)
      ENDIF
      ONEMWF  = ONE - WF
!
! *** FIND FIRST SECTION ; DRY ONE ************************************
!
      CALL DRYCASE
      IF (ABS(ONEMWF).LE.1D-5) GOTO 200  ! DRY AEROSOL
!
      CNH42SO = CNH42S4              ! FIRST (DRY) SOLUTION
      CNH4HSO = CNH4HS4
      CLCO    = CLC
      CNH4N3O = CNH4NO3
      CNH4CLO = CNH4CL
      CNA2SO  = CNA2SO4
      CNAHSO  = CNAHSO4
      CNANO   = CNANO3
      CNACLO  = CNACL
!
! *** FIND SECOND SECTION ; DRY & LIQUID ******************************
!
      CNH42S4 = ZERO
      CNH4HS4 = ZERO
      CLC     = ZERO
      CNH4NO3 = ZERO
      CNH4CL  = ZERO
      CNA2SO4 = ZERO
      CNAHSO4 = ZERO
      CNANO3  = ZERO
      CNACL   = ZERO
      GNH3    = ZERO
      GHNO3   = ZERO
      GHCL    = ZERO
      CALL LIQCASE                   ! SECOND (LIQUID) SOLUTION
!
! *** ADJUST THINGS FOR THE CASE THAT THE LIQUID SUB PREDICTS DRY AEROSOL
!
      IF (WATER.LE.TINY) THEN
         WATER = ZERO
         DO 100 I=1,NIONS
            MOLAL(I)= ZERO
 100     CONTINUE
         CALL DRYCASE
         GOTO 200
      ENDIF
!
! *** FIND SALT DISSOLUTIONS BETWEEN DRY & LIQUID SOLUTIONS.
!
      DAMBIS  = CNH4HSO - CNH4HS4
      DSOBIS  = CNAHSO  - CNAHSO4
      DLC     = CLCO    - CLC
!
! *** FIND SOLUTION AT MDRH BY WEIGHTING DRY & LIQUID SOLUTIONS.
!
! *** SOLID
!
      CNH42S4 = WF*CNH42SO + ONEMWF*CNH42S4
      CNA2SO4 = WF*CNA2SO  + ONEMWF*CNA2SO4
      CNAHSO4 = WF*CNAHSO  + ONEMWF*CNAHSO4
      CNH4HS4 = WF*CNH4HSO + ONEMWF*CNH4HS4
      CLC     = WF*CLCO    + ONEMWF*CLC
      CNH4NO3 = WF*CNH4N3O + ONEMWF*CNH4NO3
      CNANO3  = WF*CNANO   + ONEMWF*CNANO3
      CNACL   = WF*CNACLO  + ONEMWF*CNACL
      CNH4CL  = WF*CNH4CLO + ONEMWF*CNH4CL
!
! *** LIQUID
!
      WATER   = ONEMWF*WATER
!
      MOLAL(2)= WAER(1) - 2.D0*CNA2SO4 - CNAHSO4 - CNANO3 -             &
     &                         CNACL                            ! NA+
      MOLAL(3)= WAER(3) - 2.D0*CNH42S4 - CNH4HS4 - CNH4CL -             &
     &                    3.D0*CLC     - CNH4NO3                ! NH4+
      MOLAL(4)= WAER(5) - CNACL - CNH4CL                        ! CL-
      MOLAL(7)= WAER(4) - CNANO3 - CNH4NO3                      ! NO3-
      MOLAL(6)= ONEMWF*(MOLAL(6) + DSOBIS + DAMBIS + DLC)       ! HSO4-
      MOLAL(5)= WAER(2) - MOLAL(6) - CLC - CNH42S4 - CNA2SO4    ! SO4--
!
      A8      = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.
      IF (MOLAL(5).LE.TINY) THEN
         HIEQ = SQRT(XKW *RH*WATER*WATER)  ! Neutral solution
      ELSE
         HIEQ = A8*MOLAL(6)/MOLAL(5)
      ENDIF
      HIEN    = MOLAL(4) + MOLAL(7) + MOLAL(6) + 2.D0*MOLAL(5) -        &
     &          MOLAL(2) - MOLAL(3)
      MOLAL(1)= MAX (HIEQ, HIEN)                                ! H+
!
! *** GAS (ACTIVITY COEFS FROM LIQUID SOLUTION)
!
      A2      = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2. ! NH3  <==> NH4+
      A3      = XK4 *R*TEMP*(WATER/GAMA(10))**2.        ! HNO3 <==> NO3-
      A4      = XK3 *R*TEMP*(WATER/GAMA(11))**2.        ! HCL  <==> CL-
!
      GNH3    = MOLAL(3)/MAX(MOLAL(1),TINY)/A2
      GHNO3   = MOLAL(1)*MOLAL(7)/A3
      GHCL    = MOLAL(1)*MOLAL(4)/A4
!
200   CONTINUE
      IF (LHOOK) CALL DR_HOOK('CALCMDRP',1,ZHOOK_HANDLE)
!
! *** END OF SUBROUTINE CALCMDRP ****************************************
!
      ENDSUBROUTINE CALCMDRP
!=======================================================================
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCHS4
! *** THIS SUBROUTINE CALCULATES THE HSO4 GENERATED FROM (H,SO4).
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!=======================================================================
!
      SUBROUTINE CALCHS4 (HI, SO4I, HSO4I, DELTA)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
      IMPLICIT NONE
      INCLUDE 'isrpia.inc'
      REAL(KIND=8) :: A8,BB,CC,DD,DELTA1,DELTA2,SQDD
      REAL(KIND=8) :: HI,SO4I,HSO4I,DELTA
!C      CHARACTER ERRINF*40
!
! *** IF TOO LITTLE WATER, DONT SOLVE
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('CALCHS4',0,ZHOOK_HANDLE)
      IF (WATER.LE.1d1*TINY) THEN
         DELTA = ZERO
         IF (LHOOK) CALL DR_HOOK('CALCHS4',1,ZHOOK_HANDLE)
         RETURN
      ENDIF
!
! *** CALCULATE HSO4 SPECIATION *****************************************
!
      A8 = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.
!
      BB =-(HI + SO4I + A8)
      CC = HI*SO4I - HSO4I*A8
      DD = BB*BB - 4.D0*CC
!
      IF (DD.GE.ZERO) THEN
         SQDD   = SQRT(DD)
         DELTA1 = 0.5*(-BB + SQDD)
         DELTA2 = 0.5*(-BB - SQDD)
         IF (HSO4I.LE.TINY) THEN
            DELTA = DELTA2
         ELSEIF( HI*SO4I .GE. A8*HSO4I ) THEN
            DELTA = DELTA2
         ELSEIF( HI*SO4I .LT. A8*HSO4I ) THEN
            DELTA = DELTA1
         ELSE
            DELTA = ZERO
         ENDIF
      ELSE
         DELTA  = ZERO
      ENDIF
!CC
!CC *** COMPARE DELTA TO TOTAL H+ ; ESTIMATE EFFECT OF HSO4 ***************
!CC
!C      HYD = MAX(HI, MOLAL(1))
!C      IF (HYD.GT.TINY) THEN
!C         IF (DELTA/HYD.GT.0.1d0) THEN
!C            WRITE (ERRINF,'(1PE10.3)') DELTA/HYD*100.0
!C            CALL PUSHERR (0020, ERRINF)
!C         ENDIF
!C      ENDIF
!
      IF (LHOOK) CALL DR_HOOK('CALCHS4',1,ZHOOK_HANDLE)
!
! *** END OF SUBROUTINE CALCHS4 *****************************************
!
      ENDSUBROUTINE CALCHS4


!=======================================================================
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCPH
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!=======================================================================
!
      SUBROUTINE CALCPH (GG, HI, OHI)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
      IMPLICIT NONE
      INCLUDE 'isrpia.inc'
      REAL(KIND=8) :: AKW,BB,CC,DD,CN
      REAL(KIND=8) :: GG,HI,OHI
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('CALCPH',0,ZHOOK_HANDLE)
      AKW  = XKW *RH*WATER*WATER
      CN   = SQRT(AKW)
!
! *** GG = (negative charge) - (positive charge)
!
      IF (GG.GT.TINY) THEN                        ! H+ in excess
         BB =-GG
         CC =-AKW
         DD = BB*BB - 4.D0*CC
         HI = MAX(0.5D0*(-BB + SQRT(DD)),CN)
         OHI= AKW/HI
      ELSE                                        ! OH- in excess
         BB = GG
         CC =-AKW
         DD = BB*BB - 4.D0*CC
         OHI= MAX(0.5D0*(-BB + SQRT(DD)),CN)
         HI = AKW/OHI
      ENDIF
!
      IF (LHOOK) CALL DR_HOOK('CALCPH',1,ZHOOK_HANDLE)
!
! *** END OF SUBROUTINE CALCPH ******************************************
!
      ENDSUBROUTINE CALCPH

!=======================================================================
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCACT
! *** CALCULATES MULTI-COMPONENET ACTIVITY COEFFICIENTS FROM BROMLEYS
!     METHOD. THE BINARY ACTIVITY COEFFICIENTS ARE CALCULATED BY
!     KUSIK-MEISNER RELATION (SUBROUTINE KMTAB or SUBROUTINE KMFUL).
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!=======================================================================
!
      SUBROUTINE CALCACT
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
      IMPLICIT NONE
      INCLUDE 'isrpia.inc'
      INTEGER :: I,J
      REAL(KIND=8) :: ERRIN,ERROU,G
!
      REAL EX10, URF
      REAL G0(3,4),ZPL,ZMI,AGAMA,SION,H,CH,F1(3),F2(4)
      REAL(KIND=8) :: MPL, XIJ, YJI
      PARAMETER (URF=0.5)
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      G(I,J)= (F1(I)/Z(I) + F2(J)/Z(J+3)) / (Z(I)+Z(J+3)) - H
      IF (LHOOK) CALL DR_HOOK('CALCACT',0,ZHOOK_HANDLE)
!
! *** SAVE ACTIVITIES IN OLD ARRAY *************************************
!
      IF (FRST) THEN               ! Outer loop
         DO 10 I=1,NPAIR
            GAMOU(I) = GAMA(I)
10       CONTINUE
      ENDIF
!
      DO 20 I=1,NPAIR              ! Inner loop
         GAMIN(I) = GAMA(I)
20    CONTINUE
!
! *** CALCULATE IONIC ACTIVITY OF SOLUTION *****************************
!
      IONIC=0.0
      DO 30 I=1,NIONS
         IONIC=IONIC + MOLAL(I)*Z(I)*Z(I)
30    CONTINUE
      IONIC = MAX(MIN(0.5*IONIC/WATER,20.d0), TINY)
!
! *** CALCULATE BINARY ACTIVITY COEFFICIENTS ***************************
!
!  G0(1,1)=G11;G0(1,2)=G07;G0(1,3)=G08;G0(1,4)=G10;G0(2,1)=G01;G0(2,2)=G02
!  G0(2,3)=G12;G0(2,4)=G03;G0(3,1)=G06;G0(3,2)=G04;G0(3,3)=G09;G0(3,4)=G05
!
      IF (IACALC.EQ.0) THEN              ! K.M.; FULL
         CALL KMFUL (IONIC, SNGL(TEMP),G0(2,1),G0(2,2),G0(2,4),         &
     &               G0(3,2),G0(3,4),G0(3,1),G0(1,2),G0(1,3),G0(3,3),   &
     &               G0(1,4),G0(1,1),G0(2,3))
      ELSE                               ! K.M.; TABULATED
         CALL KMTAB (IONIC, SNGL(TEMP),G0(2,1),G0(2,2),G0(2,4),         &
     &               G0(3,2),G0(3,4),G0(3,1),G0(1,2),G0(1,3),G0(3,3),   &
     &               G0(1,4),G0(1,1),G0(2,3))
      ENDIF
!
! *** CALCULATE MULTICOMPONENT ACTIVITY COEFFICIENTS *******************
!
      AGAMA = 0.511*(298.0/TEMP)**1.5    ! Debye Huckel const. at T
      SION  = SQRT(IONIC)
      H     = AGAMA*SION/(1+SION)
!
      DO 100 I=1,3
         F1(I)=0.0
         F2(I)=0.0
100   CONTINUE
      F2(4)=0.0
!
      DO 110 I=1,3
         ZPL = Z(I)
         MPL = MOLAL(I)/WATER
         DO 110 J=1,4
            ZMI   = Z(J+3)
            CH    = 0.25*(ZPL+ZMI)*(ZPL+ZMI)/IONIC
            XIJ   = CH*MPL
            YJI   = CH*MOLAL(J+3)/WATER
            F1(I) = F1(I) + SNGL(YJI*(G0(I,J) + ZPL*ZMI*H))
            F2(J) = F2(J) + SNGL(XIJ*(G0(I,J) + ZPL*ZMI*H))
110   CONTINUE
!
! *** LOG10 OF ACTIVITY COEFFICIENTS ***********************************
!
      GAMA(01) = G(2,1)*ZZ(01)                     ! NACL
      GAMA(02) = G(2,2)*ZZ(02)                     ! NA2SO4
      GAMA(03) = G(2,4)*ZZ(03)                     ! NANO3
      GAMA(04) = G(3,2)*ZZ(04)                     ! (NH4)2SO4
      GAMA(05) = G(3,4)*ZZ(05)                     ! NH4NO3
      GAMA(06) = G(3,1)*ZZ(06)                     ! NH4CL
      GAMA(07) = G(1,2)*ZZ(07)                     ! 2H-SO4
      GAMA(08) = G(1,3)*ZZ(08)                     ! H-HSO4
      GAMA(09) = G(3,3)*ZZ(09)                     ! NH4HSO4
      GAMA(10) = G(1,4)*ZZ(10)                     ! HNO3
      GAMA(11) = G(1,1)*ZZ(11)                     ! HCL
      GAMA(12) = G(2,3)*ZZ(12)                     ! NAHSO4
      GAMA(13) = 0.20*(3.0*GAMA(04)+2.0*GAMA(09))  ! LC ; SCAPE
!C      GAMA(13) = 0.50*(GAMA(04)+GAMA(09))          ! LC ; SEQUILIB
!C      GAMA(13) = 0.25*(3.0*GAMA(04)+GAMA(07))      ! LC ; AIM
!
! *** CONVERT LOG (GAMA) COEFFICIENTS TO GAMA **************************
!
      DO 200 I=1,NPAIR
!C         GAMA(I)=MAX(-5.0d0, MIN(GAMA(I),5.0d0) ) ! F77 LIBRARY ROUTINE
!C         GAMA(I)=10.0**GAMA(I)
         GAMA(I)=EX10(SNGL(GAMA(I)), 5.0)    ! CUTOFF SET TO [-5,5]
         GAMA(I) = GAMIN(I)*(1.0-URF) + URF*GAMA(I)  ! Under-relax GAMA's
  200 CONTINUE
!
! *** SETUP ACTIVITY CALCULATION FLAGS *********************************
!
! OUTER CALCULATION LOOP ; ONLY IF FRST=.TRUE.
!
      IF (FRST) THEN
         ERROU = ZERO                    ! CONVERGENCE CRITERION
         DO 210 I=1,NPAIR
            ERROU=MAX(ERROU, ABS((GAMOU(I)-GAMA(I))/GAMOU(I)))
210      CONTINUE
         CALAOU = ERROU .GE. EPSACT      ! SETUP FLAGS
         FRST   =.FALSE.
      ENDIF
!
! INNER CALCULATION LOOP ; ALWAYS
!
      ERRIN = ZERO                       ! CONVERGENCE CRITERION
      DO 220 I=1,NPAIR
         ERRIN = MAX (ERRIN, ABS((GAMIN(I)-GAMA(I))/GAMIN(I)))
220   CONTINUE
      CALAIN = ERRIN .GE. EPSACT
!
      ICLACT = ICLACT + 1                ! Increment ACTIVITY call counter
!
! *** END OF SUBROUTINE ACTIVITY ****************************************
!
      IF (LHOOK) CALL DR_HOOK('CALCACT',1,ZHOOK_HANDLE)
      ENDSUBROUTINE CALCACT


!=======================================================================
!
! *** ISORROPIA CODE
! *** SUBROUTINE RSTGAM
! *** RESETS ACTIVITY COEFFICIENT ARRAYS TO DEFAULT VALUE OF 0.1
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!=======================================================================
!
      SUBROUTINE RSTGAM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
      IMPLICIT NONE
      INCLUDE 'isrpia.inc'
      INTEGER :: I
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('RSTGAM',0,ZHOOK_HANDLE)
      DO 10 I=1, NPAIR
         GAMA(I) = 0.1
10    CONTINUE
!
! *** END OF SUBROUTINE RSTGAM ******************************************
!
      IF (LHOOK) CALL DR_HOOK('RSTGAM',1,ZHOOK_HANDLE)
      ENDSUBROUTINE RSTGAM
!=======================================================================
!
! *** ISORROPIA CODE
! *** SUBROUTINE KMFUL
! *** CALCULATES BINARY ACTIVITY COEFFICIENTS BY KUSIK-MEISSNER METHOD.
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!=======================================================================
!
      SUBROUTINE KMFUL (IONIC,TEMP,G01,G02,G03,G04,G05,G06,G07,G08,G09, &
     &                  G10,G11,G12)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
      IMPLICIT NONE
      REAL :: CF1,CF2,SION,TI
      REAL Ionic, TEMP
      REAL :: G01,G02,G03,G04,G05,G06,G07,G08,G09,G10,G11,G12
      REAL :: Z01,Z02,Z03,Z04,Z05,Z06,Z07,Z08,Z10,Z11
      DATA Z01,Z02,Z03,Z04,Z05,Z06,Z07,Z08,Z10,Z11                      &
     &    /1,  2,  1,  2,  1,  1,  2,  1,  1,  1/
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('KMFUL',0,ZHOOK_HANDLE)
      SION = SQRT(IONIC)
!
! *** Coefficients at 25 oC
!
      CALL MKBI(2.230, IONIC, SION, Z01, G01)
      CALL MKBI(-0.19, IONIC, SION, Z02, G02)
      CALL MKBI(-0.39, IONIC, SION, Z03, G03)
      CALL MKBI(-0.25, IONIC, SION, Z04, G04)
      CALL MKBI(-1.15, IONIC, SION, Z05, G05)
      CALL MKBI(0.820, IONIC, SION, Z06, G06)
      CALL MKBI(-.100, IONIC, SION, Z07, G07)
      CALL MKBI(8.000, IONIC, SION, Z08, G08)
      CALL MKBI(2.600, IONIC, SION, Z10, G10)
      CALL MKBI(6.000, IONIC, SION, Z11, G11)
!
! *** Correct for T other than 298 K
!
      TI  = TEMP-273.0
      IF (ABS(TI) .GT. 1.0) THEN
         CF1 = 1.125-0.005*TI
         CF2 = (0.125-0.005*TI)*(0.039*IONIC**0.92-0.41*SION/(1.+SION))
         G01 = CF1*G01 - CF2*Z01
         G02 = CF1*G02 - CF2*Z02
         G03 = CF1*G03 - CF2*Z03
         G04 = CF1*G04 - CF2*Z04
         G05 = CF1*G05 - CF2*Z05
         G06 = CF1*G06 - CF2*Z06
         G07 = CF1*G07 - CF2*Z07
         G08 = CF1*G08 - CF2*Z08
         G10 = CF1*G10 - CF2*Z10
         G11 = CF1*G11 - CF2*Z11
      ENDIF
!
      G09 = G06 + G08 - G11
      G12 = G01 + G08 - G11
!
! *** Return point ; End of subroutine
!
      IF (LHOOK) CALL DR_HOOK('KMFUL',1,ZHOOK_HANDLE)
      ENDSUBROUTINE KMFUL


!=======================================================================
!
! *** ISORROPIA CODE
! *** SUBROUTINE MKBI
! *** CALCULATES BINARY ACTIVITY COEFFICIENTS BY KUSIK-MEISSNER METHOD.
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!=======================================================================
!
      SUBROUTINE MKBI(Q,IONIC,SION,ZIP,BI)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
      IMPLICIT NONE
!
      REAL :: Q,SION,ZIP,BI
      REAL :: B,C,XX
      REAL IONIC
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('MKBI',0,ZHOOK_HANDLE)
      B=.75-.065*Q
      C= 1.0
      IF (IONIC.LT.6.0) C=1.+.055*Q*EXP(-.023*IONIC*IONIC*IONIC)
      XX=-0.5107*SION/(1.+C*SION)
      BI=(1.+B*(1.+.1*IONIC)**Q-B)
      BI=ZIP*LOG10(BI) + ZIP*XX
!
      IF (LHOOK) CALL DR_HOOK('MKBI',1,ZHOOK_HANDLE)
      ENDSUBROUTINE MKBI
!=======================================================================
!
! *** ISORROPIA CODE
! *** SUBROUTINE KMTAB
! *** CALCULATES BINARY ACTIVITY COEFFICIENTS BY KUSIK-MEISSNER METHOD.
!     THE COMPUTATIONS HAVE BEEN PERFORMED AND THE RESULTS ARE STORED IN
!     LOOKUP TABLES. THE IONIC ACTIVITY 'IONIC' IS INPUT, AND THE ARRAY
!     'BINARR' IS RETURNED WITH THE BINARY COEFFICIENTS.
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!=======================================================================
!
      SUBROUTINE KMTAB (IN,TEMP,G01,G02,G03,G04,G05,G06,G07,G08,G09,G10,&
     &                          G11,G12)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
      IMPLICIT NONE
      INTEGER :: IND
      REAL :: IN,TEMP,G01,G02,G03,G04,G05,G06,G07,G08,G09,G10,G11,G12
!
! *** Find temperature range
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('KMTAB',0,ZHOOK_HANDLE)
      IND = NINT((TEMP-198.0)/25.0) + 1
      IND = MIN(MAX(IND,1),6)
!
! *** Call appropriate routine
!
      IF (IND.EQ.1) THEN
         CALL KM198 (IN,G01,G02,G03,G04,G05,G06,G07,G08,G09,G10,G11,G12)
      ELSEIF (IND.EQ.2) THEN
         CALL KM223 (IN,G01,G02,G03,G04,G05,G06,G07,G08,G09,G10,G11,G12)
      ELSEIF (IND.EQ.3) THEN
         CALL KM248 (IN,G01,G02,G03,G04,G05,G06,G07,G08,G09,G10,G11,G12)
      ELSEIF (IND.EQ.4) THEN
         CALL KM273 (IN,G01,G02,G03,G04,G05,G06,G07,G08,G09,G10,G11,G12)
      ELSEIF (IND.EQ.5) THEN
         CALL KM298 (IN,G01,G02,G03,G04,G05,G06,G07,G08,G09,G10,G11,G12)
      ELSE
         CALL KM323 (IN,G01,G02,G03,G04,G05,G06,G07,G08,G09,G10,G11,G12)
      ENDIF
!
! *** Return point; End of subroutine
!
      IF (LHOOK) CALL DR_HOOK('KMTAB',1,ZHOOK_HANDLE)
      ENDSUBROUTINE KMTAB


!=======================================================================
!
! *** ISORROPIA CODE
! *** SUBROUTINE KM198
! *** CALCULATES BINARY ACTIVITY COEFFICIENTS BY KUSIK-MEISSNER METHOD.
!     THE COMPUTATIONS HAVE BEEN PERFORMED AND THE RESULTS ARE STORED IN
!     LOOKUP TABLES. THE IONIC ACTIVITY 'IN' IS INPUT, AND THE ARRAY
!     'BINARR' IS RETURNED WITH THE BINARY COEFFICIENTS.
!
!     TEMPERATURE IS 198K
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!=======================================================================
!
      SUBROUTINE KM198 (IN,G01,G02,G03,G04,G05,G06,G07,G08,G09,G10,     &
     &                     G11,G12)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
      IMPLICIT NONE
      INTEGER :: IPOS
      REAL :: IN,G01,G02,G03,G04,G05,G06,G07,G08,G09,G10,G11,G12
      REAL :: BNC01M,BNC02M,BNC03M,BNC04M,BNC05M,BNC06M,BNC07M
      REAL :: BNC08M,BNC09M,BNC10M,BNC11M,BNC12M,BNC13M
!
! *** Common block definition
!
      COMMON /KMC198/                                                   &
     &BNC01M(  741),BNC02M(  741),BNC03M(  741),BNC04M(  741),          &
     &BNC05M(  741),BNC06M(  741),BNC07M(  741),BNC08M(  741),          &
     &BNC09M(  741),BNC10M(  741),BNC11M(  741),BNC12M(  741),          &
     &BNC13M(  741)
!      REAL IN
!
! *** Find position in arrays for bincoef
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('KM198',0,ZHOOK_HANDLE)
      IF (IN.LE. 0.300000E+02) THEN
         ipos = NINT( 0.200000E+02*IN) + 1
      ELSE
         ipos =   600+NINT( 0.200000E+01*IN- 0.300000E+02)
      ENDIF
      ipos = min(ipos,  741)
!
! *** Assign values to return array
!
      G01 = BNC01M(ipos)
      G02 = BNC02M(ipos)
      G03 = BNC03M(ipos)
      G04 = BNC04M(ipos)
      G05 = BNC05M(ipos)
      G06 = BNC06M(ipos)
      G07 = BNC07M(ipos)
      G08 = BNC08M(ipos)
      G09 = BNC09M(ipos)
      G10 = BNC10M(ipos)
      G11 = BNC11M(ipos)
      G12 = BNC12M(ipos)
!
! *** Return point ; End of subroutine
!
      IF (LHOOK) CALL DR_HOOK('KM198',1,ZHOOK_HANDLE)
      ENDSUBROUTINE KM198


!=======================================================================
!
! *** ISORROPIA CODE
! *** SUBROUTINE KM223
! *** CALCULATES BINARY ACTIVITY COEFFICIENTS BY KUSIK-MEISSNER METHOD.
!     THE COMPUTATIONS HAVE BEEN PERFORMED AND THE RESULTS ARE STORED IN
!     LOOKUP TABLES. THE IONIC ACTIVITY 'IN' IS INPUT, AND THE ARRAY
!     'BINARR' IS RETURNED WITH THE BINARY COEFFICIENTS.
!
!     TEMPERATURE IS 223K
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!=======================================================================
!
      SUBROUTINE KM223 (IN,G01,G02,G03,G04,G05,G06,G07,G08,G09,G10,     &
     &                     G11,G12)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
      IMPLICIT NONE
!
      REAL :: IN,G01,G02,G03,G04,G05,G06,G07,G08,G09,G10,G11,G12
!
! *** Common block definition
!
      INTEGER :: IPOS
      REAL :: BNC01M,BNC02M,BNC03M,BNC04M,BNC05M,BNC06M,BNC07M
      REAL :: BNC08M,BNC09M,BNC10M,BNC11M,BNC12M,BNC13M
      COMMON /KMC223/                                                   &
     &BNC01M(  741),BNC02M(  741),BNC03M(  741),BNC04M(  741),          &
     &BNC05M(  741),BNC06M(  741),BNC07M(  741),BNC08M(  741),          &
     &BNC09M(  741),BNC10M(  741),BNC11M(  741),BNC12M(  741),          &
     &BNC13M(  741)
!      REAL IN
!
! *** Find position in arrays for bincoef
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('KM223',0,ZHOOK_HANDLE)
      IF (IN.LE. 0.300000E+02) THEN
         ipos = NINT( 0.200000E+02*IN) + 1
      ELSE
         ipos =   600+NINT( 0.200000E+01*IN- 0.300000E+02)
      ENDIF
      ipos = min(ipos,  741)
!
! *** Assign values to return array
!
      G01 = BNC01M(ipos)
      G02 = BNC02M(ipos)
      G03 = BNC03M(ipos)
      G04 = BNC04M(ipos)
      G05 = BNC05M(ipos)
      G06 = BNC06M(ipos)
      G07 = BNC07M(ipos)
      G08 = BNC08M(ipos)
      G09 = BNC09M(ipos)
      G10 = BNC10M(ipos)
      G11 = BNC11M(ipos)
      G12 = BNC12M(ipos)
!
! *** Return point ; End of subroutine
!
      IF (LHOOK) CALL DR_HOOK('KM223',1,ZHOOK_HANDLE)
      ENDSUBROUTINE KM223



!=======================================================================
!
! *** ISORROPIA CODE
! *** SUBROUTINE KM248
! *** CALCULATES BINARY ACTIVITY COEFFICIENTS BY KUSIK-MEISSNER METHOD.
!     THE COMPUTATIONS HAVE BEEN PERFORMED AND THE RESULTS ARE STORED IN
!     LOOKUP TABLES. THE IONIC ACTIVITY 'IN' IS INPUT, AND THE ARRAY
!     'BINARR' IS RETURNED WITH THE BINARY COEFFICIENTS.
!
!     TEMPERATURE IS 248K
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!=======================================================================
!
      SUBROUTINE KM248 (IN,G01,G02,G03,G04,G05,G06,G07,G08,G09,G10,     &
     &                     G11,G12)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
      IMPLICIT NONE
!
      INTEGER :: IPOS
      REAL :: IN,G01,G02,G03,G04,G05,G06,G07,G08,G09,G10,G11,G12
!
! *** Common block definition
!
      REAL :: BNC01M,BNC02M,BNC03M,BNC04M,BNC05M,BNC06M,BNC07M
      REAL :: BNC08M,BNC09M,BNC10M,BNC11M,BNC12M,BNC13M
      COMMON /KMC248/                                                   &
     &BNC01M(  741),BNC02M(  741),BNC03M(  741),BNC04M(  741),          &
     &BNC05M(  741),BNC06M(  741),BNC07M(  741),BNC08M(  741),          &
     &BNC09M(  741),BNC10M(  741),BNC11M(  741),BNC12M(  741),          &
     &BNC13M(  741)
!      REAL IN
!
! *** Find position in arrays for bincoef
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('KM248',0,ZHOOK_HANDLE)
      IF (IN.LE. 0.300000E+02) THEN
         ipos = NINT( 0.200000E+02*IN) + 1
      ELSE
         ipos =   600+NINT( 0.200000E+01*IN- 0.300000E+02)
      ENDIF
      ipos = min(ipos,  741)
!
! *** Assign values to return array
!
      G01 = BNC01M(ipos)
      G02 = BNC02M(ipos)
      G03 = BNC03M(ipos)
      G04 = BNC04M(ipos)
      G05 = BNC05M(ipos)
      G06 = BNC06M(ipos)
      G07 = BNC07M(ipos)
      G08 = BNC08M(ipos)
      G09 = BNC09M(ipos)
      G10 = BNC10M(ipos)
      G11 = BNC11M(ipos)
      G12 = BNC12M(ipos)
!
! *** Return point ; End of subroutine
!
      IF (LHOOK) CALL DR_HOOK('KM248',1,ZHOOK_HANDLE)
      ENDSUBROUTINE KM248


!=======================================================================
!
! *** ISORROPIA CODE
! *** SUBROUTINE KMTAB273
! *** CALCULATES BINARY ACTIVITY COEFFICIENTS BY KUSIK-MEISSNER METHOD.
!     THE COMPUTATIONS HAVE BEEN PERFORMED AND THE RESULTS ARE STORED IN
!     LOOKUP TABLES. THE IONIC ACTIVITY 'IN' IS INPUT, AND THE ARRAY
!     'BINARR' IS RETURNED WITH THE BINARY COEFFICIENTS.
!
!     TEMPERATURE IS 273K
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!=======================================================================
!
      SUBROUTINE KM273 (IN,G01,G02,G03,G04,G05,G06,G07,G08,G09,G10,     &
     &                     G11,G12)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
      IMPLICIT NONE
!
      INTEGER :: IPOS
      REAL :: IN,G01,G02,G03,G04,G05,G06,G07,G08,G09,G10,G11,G12
!
! *** Common block definition
!
      REAL :: BNC01M,BNC02M,BNC03M,BNC04M,BNC05M,BNC06M,BNC07M
      REAL :: BNC08M,BNC09M,BNC10M,BNC11M,BNC12M,BNC13M
      COMMON /KMC273/                                                   &
     &BNC01M(  741),BNC02M(  741),BNC03M(  741),BNC04M(  741),          &
     &BNC05M(  741),BNC06M(  741),BNC07M(  741),BNC08M(  741),          &
     &BNC09M(  741),BNC10M(  741),BNC11M(  741),BNC12M(  741),          &
     &BNC13M(  741)
!      REAL IN
!
! *** Find position in arrays for bincoef
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('KM273',0,ZHOOK_HANDLE)
      IF (IN.LE. 0.300000E+02) THEN
         ipos = NINT( 0.200000E+02*IN) + 1
      ELSE
         ipos =   600+NINT( 0.200000E+01*IN- 0.300000E+02)
      ENDIF
      ipos = min(ipos,  741)
!
! *** Assign values to return array
!
      G01 = BNC01M(ipos)
      G02 = BNC02M(ipos)
      G03 = BNC03M(ipos)
      G04 = BNC04M(ipos)
      G05 = BNC05M(ipos)
      G06 = BNC06M(ipos)
      G07 = BNC07M(ipos)
      G08 = BNC08M(ipos)
      G09 = BNC09M(ipos)
      G10 = BNC10M(ipos)
      G11 = BNC11M(ipos)
      G12 = BNC12M(ipos)
!
! *** Return point ; End of subroutine
!
      IF (LHOOK) CALL DR_HOOK('KM273',1,ZHOOK_HANDLE)
      ENDSUBROUTINE KM273


!=======================================================================
!
! *** ISORROPIA CODE
! *** SUBROUTINE KM298
! *** CALCULATES BINARY ACTIVITY COEFFICIENTS BY KUSIK-MEISSNER METHOD.
!     THE COMPUTATIONS HAVE BEEN PERFORMED AND THE RESULTS ARE STORED IN
!     LOOKUP TABLES. THE IONIC ACTIVITY 'IN' IS INPUT, AND THE ARRAY
!     'BINARR' IS RETURNED WITH THE BINARY COEFFICIENTS.
!
!     TEMPERATURE IS 298K
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!=======================================================================
!
      SUBROUTINE KM298 (IN,G01,G02,G03,G04,G05,G06,G07,G08,G09,G10,     &
     &                     G11,G12)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
      IMPLICIT NONE
!
      INTEGER :: IPOS
      REAL :: IN,G01,G02,G03,G04,G05,G06,G07,G08,G09,G10,G11,G12
!
! *** Common block definition
!
      REAL :: BNC01M,BNC02M,BNC03M,BNC04M,BNC05M,BNC06M,BNC07M
      REAL :: BNC08M,BNC09M,BNC10M,BNC11M,BNC12M,BNC13M
      COMMON /KMC298/                                                   &
     &BNC01M(  741),BNC02M(  741),BNC03M(  741),BNC04M(  741),          &
     &BNC05M(  741),BNC06M(  741),BNC07M(  741),BNC08M(  741),          &
     &BNC09M(  741),BNC10M(  741),BNC11M(  741),BNC12M(  741),          &
     &BNC13M(  741)
!      REAL IN
!
! *** Find position in arrays for bincoef
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('KM298',0,ZHOOK_HANDLE)
      IF (IN.LE. 0.300000E+02) THEN
         ipos = NINT( 0.200000E+02*IN) + 1
      ELSE
         ipos =   600+NINT( 0.200000E+01*IN- 0.300000E+02)
      ENDIF
      ipos = min(ipos,  741)
!
! *** Assign values to return array
!
      G01 = BNC01M(ipos)
      G02 = BNC02M(ipos)
      G03 = BNC03M(ipos)
      G04 = BNC04M(ipos)
      G05 = BNC05M(ipos)
      G06 = BNC06M(ipos)
      G07 = BNC07M(ipos)
      G08 = BNC08M(ipos)
      G09 = BNC09M(ipos)
      G10 = BNC10M(ipos)
      G11 = BNC11M(ipos)
      G12 = BNC12M(ipos)
!
! *** Return point ; End of subroutine
!
      IF (LHOOK) CALL DR_HOOK('KM298',1,ZHOOK_HANDLE)
      ENDSUBROUTINE KM298


!=======================================================================
!
! *** ISORROPIA CODE
! *** SUBROUTINE KM323
! *** CALCULATES BINARY ACTIVITY COEFFICIENTS BY KUSIK-MEISSNER METHOD.
!     THE COMPUTATIONS HAVE BEEN PERFORMED AND THE RESULTS ARE STORED IN
!     LOOKUP TABLES. THE IONIC ACTIVITY 'IN' IS INPUT, AND THE ARRAY
!     'BINARR' IS RETURNED WITH THE BINARY COEFFICIENTS.
!
!     TEMPERATURE IS 323K
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!=======================================================================
!
      SUBROUTINE KM323 (IN,G01,G02,G03,G04,G05,G06,G07,G08,G09,G10,     &
     &                     G11,G12)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
      IMPLICIT NONE
!
      INTEGER :: IPOS
      REAL :: IN,G01,G02,G03,G04,G05,G06,G07,G08,G09,G10,G11,G12
!
! *** Common block definition
!
      REAL :: BNC01M,BNC02M,BNC03M,BNC04M,BNC05M,BNC06M,BNC07M
      REAL :: BNC08M,BNC09M,BNC10M,BNC11M,BNC12M,BNC13M
      COMMON /KMC323/                                                   &
     &BNC01M(  741),BNC02M(  741),BNC03M(  741),BNC04M(  741),          &
     &BNC05M(  741),BNC06M(  741),BNC07M(  741),BNC08M(  741),          &
     &BNC09M(  741),BNC10M(  741),BNC11M(  741),BNC12M(  741),          &
     &BNC13M(  741)
!      REAL IN
!
! *** Find position in arrays for bincoef
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('KM323',0,ZHOOK_HANDLE)
      IF (IN.LE. 0.300000E+02) THEN
         ipos = NINT( 0.200000E+02*IN) + 1
      ELSE
         ipos =   600+NINT( 0.200000E+01*IN- 0.300000E+02)
      ENDIF
      ipos = min(ipos,  741)
!
! *** Assign values to return array
!
      G01 = BNC01M(ipos)
      G02 = BNC02M(ipos)
      G03 = BNC03M(ipos)
      G04 = BNC04M(ipos)
      G05 = BNC05M(ipos)
      G06 = BNC06M(ipos)
      G07 = BNC07M(ipos)
      G08 = BNC08M(ipos)
      G09 = BNC09M(ipos)
      G10 = BNC10M(ipos)
      G11 = BNC11M(ipos)
      G12 = BNC12M(ipos)
!
! *** Return point ; End of subroutine
!
      IF (LHOOK) CALL DR_HOOK('KM323',1,ZHOOK_HANDLE)
      ENDSUBROUTINE KM323

!C*************************************************************************
!C
!C  TOOLBOX LIBRARY v.1.0 (May 1995)
!C
!C  Program unit   : SUBROUTINE CHRBLN
!C  Purpose        : Position of last non-blank character in a string
!C  Author         : Athanasios Nenes
!C
!C  ======================= ARGUMENTS / USAGE =============================
!C
!C  STR        is the CHARACTER variable containing the string examined
!C  IBLK       is a INTEGER variable containing the position of last non
!C             blank character. If string is all spaces (ie '   '), then
!C             the value returned is 1.
!C
!C  EXAMPLE:
!C             STR = 'TEST1.DAT     '
!C             CALL CHRBLN (STR, IBLK)
!C
!C  after execution of this code segment, "IBLK" has the value "9", which
!C  is the position of the last non-blank character of "STR".
!C
!C***********************************************************************
!C
      SUBROUTINE CHRBLN (STR, IBLK)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
      IMPLICIT NONE
!C
!C***********************************************************************
      INTEGER :: IBLK
      INTEGER ::I,ILEN
      CHARACTER*(*) STR
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('CHRBLN',0,ZHOOK_HANDLE)
      IBLK = 1                       ! Substring pointer (default=1)
      ILEN = LEN(STR)                ! Length of string
      DO 10 i=ILEN,1,-1
         IF (STR(i:i).NE.' ' .AND. STR(i:i).NE.CHAR(0)) THEN
            IBLK = i
            IF (LHOOK) CALL DR_HOOK('CHRBLN',1,ZHOOK_HANDLE)
            RETURN
         ENDIF
10    CONTINUE
      IF (LHOOK) CALL DR_HOOK('CHRBLN',1,ZHOOK_HANDLE)
!
      ENDSUBROUTINE CHRBLN




!C*************************************************************************
!C
!C  TOOLBOX LIBRARY v.1.0 (May 1995)
!C
!C  Program unit   : SUBROUTINE SHFTRGHT
!C  Purpose        : RIGHT-JUSTIFICATION FUNCTION ON A STRING
!C  Author         : Athanasios Nenes
!C
!C  ======================= ARGUMENTS / USAGE =============================
!C
!C  STRING     is the CHARACTER variable with the string to be justified
!C
!C  EXAMPLE:
!C             STRING    = 'AAAA    '
!C             CALL SHFTRGHT (STRING)
!C
!C  after execution of this code segment, STRING contains the value
!C  '    AAAA'.
!C
!C*************************************************************************
!C
      SUBROUTINE SHFTRGHT (CHR)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
      IMPLICIT NONE
!C
!C***********************************************************************
      INTEGER :: I,I1,I2
      CHARACTER CHR*(*)
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('SHFTRGHT',0,ZHOOK_HANDLE)
      I1  = LEN(CHR)             ! Total length of string
      CALL CHRBLN(CHR,I2)        ! Position of last non-blank character
      IF (I2.EQ.I1)  THEN
        IF (LHOOK) CALL DR_HOOK('SHFTRGHT',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
!
      DO 10 I=I2,1,-1            ! Shift characters
         CHR(I1+I-I2:I1+I-I2) = CHR(I:I)
         CHR(I:I) = ' '
10    CONTINUE
      IF (LHOOK) CALL DR_HOOK('SHFTRGHT',1,ZHOOK_HANDLE)
!
      ENDSUBROUTINE SHFTRGHT




!C*************************************************************************
!C
!C  TOOLBOX LIBRARY v.1.0 (May 1995)
!C
!C  Program unit   : SUBROUTINE RPLSTR
!C  Purpose        : REPLACE CHARACTERS OCCURING IN A STRING
!C  Author         : Athanasios Nenes
!C
!C  ======================= ARGUMENTS / USAGE =============================
!C
!C  STRING     is the CHARACTER variable with the string to be edited
!C  OLD        is the old character which is to be replaced
!C  NEW        is the new character which OLD is to be replaced with
!C  IERR       is 0 if everything went well, is 1 if 'NEW' contains 'OLD'.
!C             In this case, this is invalid, and no change is done.
!C
!C  EXAMPLE:
!C             STRING    = 'AAAA'
!C             OLD       = 'A'
!C             NEW       = 'B'
!C             CALL RPLSTR (STRING, OLD, NEW)
!C
!C  after execution of this code segment, STRING contains the value
!C  'BBBB'.
!C
!C*************************************************************************
!C
      SUBROUTINE RPLSTR (STRING, OLD, NEW, IERR)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
      IMPLICIT NONE
!C
!C***********************************************************************
      INTEGER :: IERR
      INTEGER :: IP,ILO
      CHARACTER STRING*(*), OLD*(*), NEW*(*)
!
! *** INITIALIZE ********************************************************
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('RPLSTR',0,ZHOOK_HANDLE)
      ILO = LEN(OLD)
!
! *** CHECK AND SEE IF 'NEW' CONTAINS 'OLD', WHICH CANNOT ***************
!
      IP = INDEX(NEW,OLD)
      IF (IP.NE.0) THEN
         IERR = 1
         IF (LHOOK) CALL DR_HOOK('RPLSTR',1,ZHOOK_HANDLE)
         RETURN
      ELSE
         IERR = 0
      ENDIF
!
! *** PROCEED WITH REPLACING *******************************************
!
10    IP = INDEX(STRING,OLD)      ! SEE IF 'OLD' EXISTS IN 'STRING'
      IF (IP.EQ.0)  THEN          ! 'OLD' DOES NOT EXIST ; RETURN
        IF (LHOOK) CALL DR_HOOK('RPLSTR',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
      STRING(IP:IP+ILO-1) = NEW   ! REPLACE SUBSTRING 'OLD' WITH 'NEW'
      GOTO 10                     ! GO FOR NEW OCCURANCE OF 'OLD'
!
      IF (LHOOK) CALL DR_HOOK('RPLSTR',1,ZHOOK_HANDLE)
      ENDSUBROUTINE RPLSTR


!C*************************************************************************
!C
!C  TOOLBOX LIBRARY v.1.0 (May 1995)
!C
!C  Program unit   : SUBROUTINE INPTD
!C  Purpose        : Prompts user for a value (DOUBLE). A default value
!C                   is provided, so if user presses <Enter>, the default
!C                   is used.
!C  Author         : Athanasios Nenes
!C
!C  ======================= ARGUMENTS / USAGE =============================
!C
!C  VAR        is the DOUBLE PRECISION variable which value is to be saved
!C  DEF        is a DOUBLE PRECISION variable, with the default value of VAR.
!C  PROMPT     is a CHARACTER varible containing the prompt string.
!C  PRFMT      is a CHARACTER variable containing the FORMAT specifier
!C             for the default value DEF.
!C  IERR       is an INTEGER error flag, and has the values:
!C             0 - No error detected.
!C             1 - Invalid FORMAT and/or Invalid default value.
!C             2 - Bad value specified by user
!C
!C  EXAMPLE:
!C             CALL INPTD (VAR, 1.0D0, 'Give value for A ', '*', Ierr)
!C
!C  after execution of this code segment, the user is prompted for the
!C  value of variable VAR. If <Enter> is pressed (ie no value is specified)
!C  then 1.0 is assigned to VAR. The default value is displayed in free-
!C  format. The error status is specified by variable Ierr
!C
!C***********************************************************************
!C
      SUBROUTINE INPTD (VAR, DEF, PROMPT, PRFMT, IERR)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
      IMPLICIT NONE
!C
!C***********************************************************************
      INTEGER :: IEND
      CHARACTER PROMPT*(*), PRFMT*(*), BUFFER*128
      REAL(KIND=8) :: DEF, VAR
      INTEGER IERR
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('INPTD',0,ZHOOK_HANDLE)
      IERR = 0
!
! *** WRITE DEFAULT VALUE TO WORK BUFFER *******************************
!
      WRITE (BUFFER, FMT=PRFMT, ERR=10) DEF
      CALL CHRBLN (BUFFER, IEND)
!
! *** PROMPT USER FOR INPUT AND READ IT ********************************
!
      WRITE (*,*) PROMPT,' [',BUFFER(1:IEND),']: '
      READ  (*, '(A)', ERR=20, END=20) BUFFER
      CALL CHRBLN (BUFFER,IEND)
!
! *** READ DATA OR SET DEFAULT ? ****************************************
!
      IF (IEND.EQ.1 .AND. BUFFER(1:1).EQ.' ') THEN
         VAR = DEF
      ELSE
         READ (BUFFER, *, ERR=20, END=20) VAR
      ENDIF
!
! *** RETURN POINT ******************************************************
!
30    CONTINUE
      IF (LHOOK) CALL DR_HOOK('INPTD',1,ZHOOK_HANDLE)
      RETURN
!
! *** ERROR HANDLER *****************************************************
!
10    IERR = 1       ! Bad FORMAT and/or bad default value
      GOTO 30
!
20    IERR = 2       ! Bad number given by user
      GOTO 30
!
      IF (LHOOK) CALL DR_HOOK('INPTD',1,ZHOOK_HANDLE)
      ENDSUBROUTINE INPTD


!C*************************************************************************
!C
!C  TOOLBOX LIBRARY v.1.0 (May 1995)
!C
!C  Program unit   : SUBROUTINE Pushend
!C  Purpose        : Positions the pointer of a sequential file at its end
!C                   Simulates the ACCESS='APPEND' clause of a F77L OPEN
!C                   statement with Standard Fortran commands.
!C
!C  ======================= ARGUMENTS / USAGE =============================
!C
!C  Iunit      is a INTEGER variable, the file unit which the file is
!C             connected to.
!C
!C  EXAMPLE:
!C             CALL PUSHEND (10)
!C
!C  after execution of this code segment, the pointer of unit 10 is
!C  pushed to its end.
!C
!C***********************************************************************
!C
      SUBROUTINE Pushend (Iunit)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
      IMPLICIT NONE
!C
!C***********************************************************************
!
      INTEGER :: IUNIT
      LOGICAL OPNED
!
! *** INQUIRE IF Iunit CONNECTED TO FILE ********************************
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('PUSHEND',0,ZHOOK_HANDLE)
      INQUIRE (UNIT=Iunit, OPENED=OPNED)
      IF (.NOT.OPNED) GOTO 25
!
! *** Iunit CONNECTED, PUSH POINTER TO END ******************************
!
10    READ (Iunit,'()', ERR=20, END=20)
      GOTO 10
!
! *** RETURN POINT ******************************************************
!
20    BACKSPACE (Iunit)
25    CONTINUE
      IF (LHOOK) CALL DR_HOOK('PUSHEND',1,ZHOOK_HANDLE)
      ENDSUBROUTINE PUSHEND



!C*************************************************************************
!C
!C  TOOLBOX LIBRARY v.1.0 (May 1995)
!C
!C  Program unit   : SUBROUTINE APPENDEXT
!C  Purpose        : Fix extension in file name string
!C
!C  ======================= ARGUMENTS / USAGE =============================
!C
!C  Filename   is the CHARACTER variable with the file name
!C  Defext     is the CHARACTER variable with extension (including '.',
!C             ex. '.DAT')
!C  Overwrite  is a LOGICAL value, .TRUE. overwrites any existing extension
!C             in "Filename" with "Defext", .FALSE. puts "Defext" only if
!C             there is no extension in "Filename".
!C
!C  EXAMPLE:
!C             FILENAME1 = 'TEST.DAT'
!C             FILENAME2 = 'TEST.DAT'
!C             CALL APPENDEXT (FILENAME1, '.TXT', .FALSE.)
!C             CALL APPENDEXT (FILENAME2, '.TXT', .TRUE. )
!C
!C  after execution of this code segment, "FILENAME1" has the value
!C  'TEST.DAT', while "FILENAME2" has the value 'TEST.TXT'
!C
!C***********************************************************************
!C
      SUBROUTINE Appendext (Filename, Defext, Overwrite)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
      IMPLICIT NONE
!C
!C***********************************************************************
      INTEGER :: IDOT,IEND
      CHARACTER*(*) Filename, Defext
      LOGICAL       Overwrite
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('APPENDEXT',0,ZHOOK_HANDLE)
      CALL CHRBLN (Filename, Iend)
      IF (Filename(1:1).EQ.' ' .AND. Iend.EQ.1)  THEN   ! Filename empty
        IF (LHOOK) CALL DR_HOOK('APPENDEXT',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
      Idot = INDEX (Filename, '.')                      ! Append extension ?
      IF (Idot.EQ.0) Filename = Filename(1:Iend)//Defext
      IF (Overwrite .AND. Idot.NE.0)                                    &
     &              Filename = Filename(:Idot-1)//Defext
      IF (LHOOK) CALL DR_HOOK('APPENDEXT',1,ZHOOK_HANDLE)
      ENDSUBROUTINE APPENDEXT





!=======================================================================
!
! *** ISORROPIA CODE
! *** SUBROUTINE POLY3
! *** FINDS THE REAL ROOTS OF THE THIRD ORDER ALGEBRAIC EQUATION:
!     X**3 + A1*X**2 + A2*X + A3 = 0.0
!     THE EQUATION IS SOLVED ANALYTICALLY.
!
!     PARAMETERS A1, A2, A3 ARE SPECIFIED BY THE USER. THE MINIMUM
!     NONEGATIVE ROOT IS RETURNED IN VARIABLE 'ROOT'. IF NO ROOT IS
!     FOUND (WHICH IS GREATER THAN ZERO), ROOT HAS THE VALUE 1D30.
!     AND THE FLAG ISLV HAS A VALUE GREATER THAN ZERO.
!
!     SOLUTION FORMULA IS FOUND IN PAGE 32 OF:
!     MATHEMATICAL HANDBOOK OF FORMULAS AND TABLES
!     SCHAUM'S OUTLINE SERIES
!     MURRAY SPIEGER, McGRAW-HILL, NEW YORK, 1968
!     (GREEK TRANSLATION: BY SOTIRIOS PERSIDES, ESPI, ATHENS, 1976)
!
!     A SPECIAL CASE IS CONSIDERED SEPERATELY ; WHEN A3 = 0, THEN
!     ONE ROOT IS X=0.0, AND THE OTHER TWO FROM THE SOLUTION OF THE
!     QUADRATIC EQUATION X**2 + A1*X + A2 = 0.0
!     THIS SPECIAL CASE IS CONSIDERED BECAUSE THE ANALYTICAL FORMULA
!     DOES NOT YIELD ACCURATE RESULTS (DUE TO NUMERICAL ROUNDOFF ERRORS)
!
! *** COPYRIGHT 1996-98, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!=======================================================================
!
      SUBROUTINE POLY3 (A1, A2, A3, ROOT, ISLV)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
      IMPLICIT NONE
!
!      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      REAL(KIND=8) :: A1,A2,A3,ROOT,COEF,D,SQD,Q,R,S,SSIG,T,THET,TSIG
      INTEGER :: ISLV,IX,I
      REAL(KIND=8) :: EXPON,ZERO,THET1,THET2,PI,EPS
      PARAMETER (EXPON=1.D0/3.D0,     ZERO=0.D0, THET1=120.D0/180.D0,   &
     &           THET2=240.D0/180.D0, PI=3.14159265358932, EPS=1D-50)
      REAL(KIND=8) :: X(3)
!
! *** SPECIAL CASE : QUADRATIC*X EQUATION *****************************
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('POLY3',0,ZHOOK_HANDLE)
      IF (ABS(A3).LE.EPS) THEN
         ISLV = 1
         IX   = 1
         X(1) = ZERO
         D    = A1*A1-4.D0*A2
         IF (D.GE.ZERO) THEN
            IX   = 3
            SQD  = SQRT(D)
            X(2) = 0.5*(-A1+SQD)
            X(3) = 0.5*(-A1-SQD)
         ENDIF
      ELSE
!
! *** NORMAL CASE : CUBIC EQUATION ************************************
!
! DEFINE PARAMETERS Q, R, S, T, D
!
         ISLV= 1
         Q   = (3.D0*A2 - A1*A1)/9.D0
         R   = (9.D0*A1*A2 - 27.D0*A3 - 2.D0*A1*A1*A1)/54.D0
         D   = Q*Q*Q + R*R
!
! *** CALCULATE ROOTS *************************************************
!
!  D < 0, THREE REAL ROOTS
!
         IF (D.LT.-EPS) THEN        ! D < -EPS  : D < ZERO
            IX   = 3
            THET = EXPON*ACOS(R/SQRT(-Q*Q*Q))
            COEF = 2.D0*SQRT(-Q)
            X(1) = COEF*COS(THET)            - EXPON*A1
            X(2) = COEF*COS(THET + THET1*PI) - EXPON*A1
            X(3) = COEF*COS(THET + THET2*PI) - EXPON*A1
!
!  D = 0, THREE REAL (ONE DOUBLE) ROOTS
!
         ELSE IF (D.LE.EPS) THEN    ! -EPS <= D <= EPS  : D = ZERO
            IX   = 2
            SSIG = SIGN (1._8, R)
            S    = SSIG*(ABS(R))**EXPON
            X(1) = 2.D0*S  - EXPON*A1
            X(2) =     -S  - EXPON*A1
!
!  D > 0, ONE REAL ROOT
!
         ELSE                       ! D > EPS  : D > ZERO
            IX   = 1
            SQD  = SQRT(D)
            SSIG = SIGN (1._8, R+SQD)       ! TRANSFER SIGN TO SSIG
            TSIG = SIGN (1._8, R-SQD)
            S    = SSIG*(ABS(R+SQD))**EXPON ! EXPONENTIATE ABS()
            T    = TSIG*(ABS(R-SQD))**EXPON
            X(1) = S + T - EXPON*A1
         ENDIF
      ENDIF
!
! *** SELECT APPROPRIATE ROOT *****************************************
!
      ROOT = 1.D30
      DO 10 I=1,IX
         IF (X(I).GT.ZERO) THEN
            ROOT = MIN (ROOT, X(I))
            ISLV = 0
         ENDIF
10    CONTINUE
!
! *** END OF SUBROUTINE POLY3 *****************************************
!
      IF (LHOOK) CALL DR_HOOK('POLY3',1,ZHOOK_HANDLE)
      ENDSUBROUTINE POLY3




!=======================================================================
!
! *** ISORROPIA CODE
! *** SUBROUTINE POLY3B
! *** FINDS A REAL ROOT OF THE THIRD ORDER ALGEBRAIC EQUATION:
!     X**3 + A1*X**2 + A2*X + A3 = 0.0
!     THE EQUATION IS SOLVED NUMERICALLY (BISECTION).
!
!     PARAMETERS A1, A2, A3 ARE SPECIFIED BY THE USER. THE MINIMUM
!     NONEGATIVE ROOT IS RETURNED IN VARIABLE 'ROOT'. IF NO ROOT IS
!     FOUND (WHICH IS GREATER THAN ZERO), ROOT HAS THE VALUE 1D30.
!     AND THE FLAG ISLV HAS A VALUE GREATER THAN ZERO.
!
!     RTLW, RTHI DEFINE THE INTERVAL WHICH THE ROOT IS LOOKED FOR.
!
!=======================================================================
!
      SUBROUTINE POLY3B (A1, A2, A3, RTLW, RTHI, ROOT, ISLV)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
      IMPLICIT NONE
!
!      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      REAL(KIND=8) :: A1,A2,A3,RTLW,RTHI,ROOT
      INTEGER :: ISLV
      REAL(KIND=8) :: DX,X1,X2,X3,Y1,Y2,Y3,X,FUNC
      INTEGER :: I
      REAL(KIND=8) :: ZERO,EPS
      INTEGER :: MAXIT,NDIV
      PARAMETER (ZERO=0.D0, EPS=1D-15, MAXIT=100, NDIV=5)
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      FUNC(X) = X**3.d0 + A1*X**2.0 + A2*X + A3
      IF (LHOOK) CALL DR_HOOK('POLY3B',0,ZHOOK_HANDLE)
!
! *** INITIAL VALUES FOR BISECTION *************************************
!
      X1   = RTLW
      Y1   = FUNC(X1)
      IF (ABS(Y1).LE.EPS) THEN     ! Is low a root?
         ROOT = RTLW
         GOTO 50
      ENDIF
!
! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO ***********************
!
      DX = (RTHI-RTLW)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1+DX
         Y2 = FUNC (X2)
         IF (SIGN(1._8,Y1)*SIGN(1._8,Y2) .LT. ZERO) GOTO 20 ! (Y1*Y2.LT.ZERO)
         X1 = X2
         Y1 = Y2
10    CONTINUE
!
! *** NO SUBDIVISION WITH SOLUTION FOUND
!
      IF (ABS(Y2) .LT. EPS) THEN   ! X2 is a root
         ROOT = X2
      ELSE
         ROOT = 1.d30
         ISLV = 1
      ENDIF
      GOTO 50
!
! *** BISECTION *******************************************************
!
20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNC (X3)
         IF (SIGN(1._8,Y1)*SIGN(1._8,Y3) .LE. ZERO) THEN  ! (Y1*Y3 .LE. ZERO)
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
!
! *** CONVERGED ; RETURN ***********************************************
!
40    X3   = 0.5*(X1+X2)
      Y3   = FUNC (X3)
      ROOT = X3
      ISLV = 0
!
50    CONTINUE
      IF (LHOOK) CALL DR_HOOK('POLY3B',1,ZHOOK_HANDLE)
!
! *** END OF SUBROUTINE POLY3B *****************************************
!
      ENDSUBROUTINE POLY3B



!cc      PROGRAM DRIVER
!cc      DOUBLE PRECISION ROOT
!ccC
!cc      CALL POLY3 (-1.d0, 1.d0, -1.d0, ROOT, ISLV)
!cc      IF (ISLV.NE.0) STOP 'Error in POLY3'
!cc      WRITE (*,*) 'Root=', ROOT
!ccC
!cc      CALL POLY3B (-1.d0, 1.d0, -1.d0, -10.d0, 10.d0, ROOT, ISLV)
!cc      IF (ISLV.NE.0) STOP 'Error in POLY3B'
!cc      WRITE (*,*) 'Root=', ROOT
!ccC
!cc      END
!=======================================================================
!
! *** ISORROPIA CODE
! *** FUNCTION EX10
! *** 10^X FUNCTION ; ALTERNATE OF LIBRARY ROUTINE ; USED BECAUSE IT IS
!     MUCH FASTER BUT WITHOUT GREAT LOSS IN ACCURACY. ,
!     MAXIMUM ERROR IS 2%, EXECUTION TIME IS 42% OF THE LIBRARY ROUTINE
!     (ON A 80286/80287 MACHINE, using Lahey FORTRAN 77 v.3.0).
!
!     EXPONENT RANGE IS BETWEEN -K AND K (K IS THE REAL ARGUMENT 'K')
!     MAX VALUE FOR K: 9.999
!     IF X < -K, X IS SET TO -K, IF X > K, X IS SET TO K
!
!     THE EXPONENT IS CALCULATED BY THE PRODUCT ADEC*AINT, WHERE ADEC
!     IS THE MANTISSA AND AINT IS THE MAGNITUDE (EXPONENT). BOTH
!     MANTISSA AND MAGNITUDE ARE PRE-CALCULATED AND STORED IN LOOKUP
!     TABLES ; THIS LEADS TO THE INCREASED SPEED.
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!=======================================================================
!
      FUNCTION EX10(X,K)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
      IMPLICIT NONE
      REAL    X, EX10, Y, AINT10, ADEC10, K
      INTEGER K1, K2
      COMMON /EXPNC/ AINT10(20), ADEC10(200)
!
! *** LIMIT X TO [-K, K] RANGE *****************************************
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('EX10',0,ZHOOK_HANDLE)
      Y    = MAX(-K, MIN(X,K))   ! MIN: -9.999, MAX: 9.999
!
! *** GET INTEGER AND DECIMAL PART *************************************
!
      K1   = INT(Y)
      K2   = INT(100*(Y-K1))
!
! *** CALCULATE EXP FUNCTION *******************************************
!
      EX10 = AINT10(K1+10)*ADEC10(K2+100)
!
! *** END OF EXP FUNCTION **********************************************
!
      IF (LHOOK) CALL DR_HOOK('EX10',1,ZHOOK_HANDLE)
      ENDFUNCTION EX10


!=======================================================================
!
! *** ISORROPIA CODE
! *** BLOCK DATA EXPON
! *** CONTAINS DATA FOR EXPONENT ARRAYS NEEDED IN FUNCTION EXP10
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!=======================================================================
!
      BLOCK DATA EXPON
      IMPLICIT NONE
!
! *** Common block definition
!
      INTEGER :: I
      REAL AINT10, ADEC10
      COMMON /EXPNC/ AINT10(20), ADEC10(200)
!
! *** Integer part
!
      DATA AINT10/                                                      &
     & 0.1000E-08, 0.1000E-07, 0.1000E-06, 0.1000E-05, 0.1000E-04,      &
     & 0.1000E-03, 0.1000E-02, 0.1000E-01, 0.1000E+00, 0.1000E+01,      &
     & 0.1000E+02, 0.1000E+03, 0.1000E+04, 0.1000E+05, 0.1000E+06,      &
     & 0.1000E+07, 0.1000E+08, 0.1000E+09, 0.1000E+10, 0.1000E+11       &
     & /
!
! *** decimal part
!
      DATA (ADEC10(I),I=1,100)/                                         &
     & 0.1023E+00, 0.1047E+00, 0.1072E+00, 0.1096E+00, 0.1122E+00,      &
     & 0.1148E+00, 0.1175E+00, 0.1202E+00, 0.1230E+00, 0.1259E+00,      &
     & 0.1288E+00, 0.1318E+00, 0.1349E+00, 0.1380E+00, 0.1413E+00,      &
     & 0.1445E+00, 0.1479E+00, 0.1514E+00, 0.1549E+00, 0.1585E+00,      &
     & 0.1622E+00, 0.1660E+00, 0.1698E+00, 0.1738E+00, 0.1778E+00,      &
     & 0.1820E+00, 0.1862E+00, 0.1905E+00, 0.1950E+00, 0.1995E+00,      &
     & 0.2042E+00, 0.2089E+00, 0.2138E+00, 0.2188E+00, 0.2239E+00,      &
     & 0.2291E+00, 0.2344E+00, 0.2399E+00, 0.2455E+00, 0.2512E+00,      &
     & 0.2570E+00, 0.2630E+00, 0.2692E+00, 0.2754E+00, 0.2818E+00,      &
     & 0.2884E+00, 0.2951E+00, 0.3020E+00, 0.3090E+00, 0.3162E+00,      &
     & 0.3236E+00, 0.3311E+00, 0.3388E+00, 0.3467E+00, 0.3548E+00,      &
     & 0.3631E+00, 0.3715E+00, 0.3802E+00, 0.3890E+00, 0.3981E+00,      &
     & 0.4074E+00, 0.4169E+00, 0.4266E+00, 0.4365E+00, 0.4467E+00,      &
     & 0.4571E+00, 0.4677E+00, 0.4786E+00, 0.4898E+00, 0.5012E+00,      &
     & 0.5129E+00, 0.5248E+00, 0.5370E+00, 0.5495E+00, 0.5623E+00,      &
     & 0.5754E+00, 0.5888E+00, 0.6026E+00, 0.6166E+00, 0.6310E+00,      &
     & 0.6457E+00, 0.6607E+00, 0.6761E+00, 0.6918E+00, 0.7079E+00,      &
     & 0.7244E+00, 0.7413E+00, 0.7586E+00, 0.7762E+00, 0.7943E+00,      &
     & 0.8128E+00, 0.8318E+00, 0.8511E+00, 0.8710E+00, 0.8913E+00,      &
     & 0.9120E+00, 0.9333E+00, 0.9550E+00, 0.9772E+00, 0.1000E+01/

      DATA (ADEC10(I),I=101,200)/                                       &
     & 0.1023E+01, 0.1047E+01, 0.1072E+01, 0.1096E+01, 0.1122E+01,      &
     & 0.1148E+01, 0.1175E+01, 0.1202E+01, 0.1230E+01, 0.1259E+01,      &
     & 0.1288E+01, 0.1318E+01, 0.1349E+01, 0.1380E+01, 0.1413E+01,      &
     & 0.1445E+01, 0.1479E+01, 0.1514E+01, 0.1549E+01, 0.1585E+01,      &
     & 0.1622E+01, 0.1660E+01, 0.1698E+01, 0.1738E+01, 0.1778E+01,      &
     & 0.1820E+01, 0.1862E+01, 0.1905E+01, 0.1950E+01, 0.1995E+01,      &
     & 0.2042E+01, 0.2089E+01, 0.2138E+01, 0.2188E+01, 0.2239E+01,      &
     & 0.2291E+01, 0.2344E+01, 0.2399E+01, 0.2455E+01, 0.2512E+01,      &
     & 0.2570E+01, 0.2630E+01, 0.2692E+01, 0.2754E+01, 0.2818E+01,      &
     & 0.2884E+01, 0.2951E+01, 0.3020E+01, 0.3090E+01, 0.3162E+01,      &
     & 0.3236E+01, 0.3311E+01, 0.3388E+01, 0.3467E+01, 0.3548E+01,      &
     & 0.3631E+01, 0.3715E+01, 0.3802E+01, 0.3890E+01, 0.3981E+01,      &
     & 0.4074E+01, 0.4169E+01, 0.4266E+01, 0.4365E+01, 0.4467E+01,      &
     & 0.4571E+01, 0.4677E+01, 0.4786E+01, 0.4898E+01, 0.5012E+01,      &
     & 0.5129E+01, 0.5248E+01, 0.5370E+01, 0.5495E+01, 0.5623E+01,      &
     & 0.5754E+01, 0.5888E+01, 0.6026E+01, 0.6166E+01, 0.6310E+01,      &
     & 0.6457E+01, 0.6607E+01, 0.6761E+01, 0.6918E+01, 0.7079E+01,      &
     & 0.7244E+01, 0.7413E+01, 0.7586E+01, 0.7762E+01, 0.7943E+01,      &
     & 0.8128E+01, 0.8318E+01, 0.8511E+01, 0.8710E+01, 0.8913E+01,      &
     & 0.9120E+01, 0.9333E+01, 0.9550E+01, 0.9772E+01, 0.1000E+02       &
     & /
!
! *** END OF BLOCK DATA EXPON ******************************************
!
      ENDBLOCKDATA EXPON
!=======================================================================
!
! *** ISORROPIA CODE
! *** SUBROUTINE ISOPLUS
! *** THIS SUBROUTINE IS THE MASTER ROUTINE FOR THE ISORROPIA-PLUS
!     THERMODYNAMIC EQUILIBRIUM AEROSOL MODEL (VERSION 1.0)
!
! *** NOTE: THIS SUBROUTINE IS INCLUDED FOR BACKWARD COMPATABILITY ONLY.
!     A PROGRAMMER SHOULD USE THE MORE COMPLETE SUBROUTINE ISOROPIA INSTEAD
!
! ======================== ARGUMENTS / USAGE ===========================
!
!  INPUT:
!  1. [WI]
!     DOUBLE PRECISION array of length [5].
!     Concentrations, expressed in moles/m3. Depending on the type of
!     problem solved, WI contains either GAS+AEROSOL or AEROSOL only
!     concentratios.
!     WI(1) - sodium
!     WI(2) - sulfate
!     WI(3) - ammonium
!     WI(4) - nitrate
!     WI(5) - chloride
!
!  2. [RHI]
!     DOUBLE PRECISION variable.
!     Ambient relative humidity expressed in a (0,1) scale.
!
!  3. [TEMPI]
!     DOUBLE PRECISION variable.
!     Ambient temperature expressed in Kelvins.
!
!  4. [IPROB]
!     INTEGER variable.
!     The type of problem solved.
!     IPROB = 0  - Forward problem is solved. In this case, array WI
!                  contains GAS and AEROSOL concentrations together.
!     IPROB = 1  - Reverse problem is solved. In this case, array WI
!                  contains AEROSOL concentrations only.
!
!  OUTPUT:
!  1. [GAS]
!     DOUBLE PRECISION array of length [03].
!     Gaseous species concentrations, expressed in moles/m3.
!     GAS(1) - NH3
!     GAS(2) - HNO3
!     GAS(3) - HCl
!
!  2. [AERLIQ]
!     DOUBLE PRECISION array of length [11].
!     Liquid aerosol species concentrations, expressed in moles/m3.
!     AERLIQ(01) - H+(aq)
!     AERLIQ(02) - Na+(aq)
!     AERLIQ(03) - NH4+(aq)
!     AERLIQ(04) - Cl-(aq)
!     AERLIQ(05) - SO4--(aq)
!     AERLIQ(06) - HSO4-(aq)
!     AERLIQ(07) - NO3-(aq)
!     AERLIQ(08) - H2O
!     AERLIQ(09) - NH3(aq) (undissociated)
!     AERLIQ(10) - HNCl(aq) (undissociated)
!     AERLIQ(11) - HNO3(aq) (undissociated)
!
!  3. [AERSLD]
!     DOUBLE PRECISION array of length [09].
!     Solid aerosol species concentrations, expressed in moles/m3.
!     AERSLD(01) - NaNO3(s)
!     AERSLD(02) - NH4NO3(s)
!     AERSLD(03) - NaCl(s)
!     AERSLD(04) - NH4Cl(s)
!     AERSLD(05) - Na2SO4(s)
!     AERSLD(06) - (NH4)2SO4(s)
!     AERSLD(07) - NaHSO4(s)
!     AERSLD(08) - NH4HSO4(s)
!     AERSLD(09) - (NH4)4H(SO4)2(s)
!
!  4. [DRY]
!     LOGICAL variable.
!     Contains information about the physical state of the system.
!     .TRUE. - There is no aqueous phase present
!     .FALSE.- There is an aqueous phase present
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!=======================================================================
!
      SUBROUTINE ISOPLUS (WI,  RHI,    TEMPI,  IPROBI,                  &
     &                    GAS, AERLIQ, AERSLD, DRYI   )
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
      IMPLICIT NONE
      INCLUDE 'isrpia.inc'
      REAL(KIND=8) :: WI,RHI,TEMPI,GAS,AERLIQ,AERSLD,DRVI
      INTEGER :: IPROBI
      INTEGER :: I
      DIMENSION WI(NCOMP), GAS(NGASAQ), AERLIQ(NIONS+NGASAQ+1),         &
     &          AERSLD(NSLDS)
      LOGICAL   DRYI
!
! *** PROBLEM TYPE (0=FOREWARD, 1=REVERSE) ******************************
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('ISOPLUS',0,ZHOOK_HANDLE)
      IPROB = IPROBI
!
! *** SOLVE FOREWARD PROBLEM ********************************************
!
      IF (IPROB.EQ.0) THEN
         IF (WI(1)+WI(2)+WI(3)+WI(4)+WI(5) .LE. TINY) THEN ! Everything=0
            CALL INIT1 (WI, RHI, TEMPI)
         ELSE IF (WI(1)+WI(4)+WI(5) .LE. TINY) THEN        ! Na,Cl,NO3=0
            CALL ISRP1F (WI, RHI, TEMPI)
         ELSE IF (WI(1)+WI(5) .LE. TINY) THEN              ! Na,Cl=0
            CALL ISRP2F (WI, RHI, TEMPI)
         ELSE
            CALL ISRP3F (WI, RHI, TEMPI)
         ENDIF
!
! *** SOLVE REVERSE PROBLEM *********************************************
!
      ELSE
         IF (WI(1)+WI(2)+WI(3)+WI(4)+WI(5) .LE. TINY) THEN ! Everything=0
            CALL INIT1 (WI, RHI, TEMPI)
         ELSE IF (WI(1)+WI(4)+WI(5) .LE. TINY) THEN        ! Na,Cl,NO3=0
            CALL ISRP1R (WI, RHI, TEMPI)
         ELSE IF (WI(1)+WI(5) .LE. TINY) THEN              ! Na,Cl=0
            CALL ISRP2R (WI, RHI, TEMPI)
         ELSE
            CALL ISRP3R (WI, RHI, TEMPI)
         ENDIF
      ENDIF
!
! *** SAVE RESULTS TO ARRAYS (units = mole/m3, kg/m3 for water) *********
!
      GAS(1) = GNH3
      GAS(2) = GHNO3
      GAS(3) = GHCL
!
      DO 10 I=1,NIONS
         AERLIQ(I) = MOLAL(I)
  10  CONTINUE
      AERLIQ(NIONS+1) = WATER*1.0D3/18.0D0
      DO 20 I=1,NGASAQ
         AERLIQ(NIONS+1+I) = GASAQ(I)
  20  CONTINUE
!
      AERSLD(1) = CNANO3
      AERSLD(2) = CNH4NO3
      AERSLD(3) = CNACL
      AERSLD(4) = CNH4CL
      AERSLD(5) = CNA2SO4
      AERSLD(6) = CNH42S4
      AERSLD(7) = CNAHSO4
      AERSLD(8) = CNH4HS4
      AERSLD(9) = CLC
!
      DRYI = WATER.LE.TINY
!
      IF (LHOOK) CALL DR_HOOK('ISOPLUS',1,ZHOOK_HANDLE)
!
! *** END OF SUBROUTINE ISOPLUS ******************************************
!
      ENDSUBROUTINE ISOPLUS




!=======================================================================
!
! *** ISORROPIA CODE
! *** SUBROUTINE ISRPIA
! *** THIS SUBROUTINE IS THE MASTER ROUTINE FOR THE ISORROPIA-PLUS
!     THERMODYNAMIC EQUILIBRIUM AEROSOL MODEL (VERSIONS 0.x)
!
! *** NOTE: THIS SUBROUTINE IS INCLUDED FOR BACKWARD COMPATABILITY ONLY.
!     A PROGRAMMER SHOULD USE THE MORE COMPLETE SUBROUTINE ISOROPIA INSTEAD
!
!
!     DEPENDING ON THE INPUT VALUES PROVIDED, THE FOLLOWING MODEL
!     SUBVERSIONS ARE CALLED:
!
!     FOREWARD PROBLEM (IPROB=0):
!     Na      SO4      NH4       NO3      CL       SUBROUTINE CALLED
!     ----    ----     ----      ----     ----     -----------------
!     0.0     >0.0     >0.0       0.0      0.0     SUBROUTINE ISRP1F
!     0.0     >0.0     >0.0      >0.0      0.0     SUBROUTINE ISRP2F
!     >0.0    >0.0     >0.0      >0.0     >0.0     SUBROUTINE ISRP3F
!
!     REVERSE PROBLEM (IPROB=1):
!     Na      SO4      NH4       NO3      CL       SUBROUTINE CALLED
!     ----    ----     ----      ----     ----     -----------------
!     0.0     >0.0     >0.0       0.0      0.0     SUBROUTINE ISRP1R
!     0.0     >0.0     >0.0      >0.0      0.0     SUBROUTINE ISRP2R
!     >0.0    >0.0     >0.0      >0.0     >0.0     SUBROUTINE ISRP3R
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!=======================================================================
!
      SUBROUTINE ISRPIA (WI, RHI, TEMPI, IPROBI)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
      IMPLICIT NONE
      INCLUDE 'isrpia.inc'
      REAL(KIND=8) :: WI,RHI,TEMPI
      INTEGER :: IPROBI
      DIMENSION WI(NCOMP)
!
! *** PROBLEM TYPE (0=FOREWARD, 1=REVERSE) ******************************
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('ISRPIA',0,ZHOOK_HANDLE)
      IPROB = IPROBI
!
! *** SOLVE FOREWARD PROBLEM ********************************************
!
      IF (IPROB.EQ.0) THEN
         IF (WI(1)+WI(2)+WI(3)+WI(4)+WI(5) .LE. TINY) THEN ! Everything=0
            CALL INIT1 (WI, RHI, TEMPI)
         ELSE IF (WI(1)+WI(4)+WI(5) .LE. TINY) THEN        ! Na,Cl,NO3=0
            CALL ISRP1F (WI, RHI, TEMPI)
         ELSE IF (WI(1)+WI(5) .LE. TINY) THEN              ! Na,Cl=0
            CALL ISRP2F (WI, RHI, TEMPI)
         ELSE
            CALL ISRP3F (WI, RHI, TEMPI)
         ENDIF
!
! *** SOLVE REVERSE PROBLEM *********************************************
!
      ELSE
         IF (WI(1)+WI(2)+WI(3)+WI(4)+WI(5) .LE. TINY) THEN ! Everything=0
            CALL INIT1 (WI, RHI, TEMPI)
         ELSE IF (WI(1)+WI(4)+WI(5) .LE. TINY) THEN        ! Na,Cl,NO3=0
            CALL ISRP1R (WI, RHI, TEMPI)
         ELSE IF (WI(1)+WI(5) .LE. TINY) THEN              ! Na,Cl=0
            CALL ISRP2R (WI, RHI, TEMPI)
         ELSE
            CALL ISRP3R (WI, RHI, TEMPI)
         ENDIF
      ENDIF
!
! *** SETUP 'DRY' FLAG ***************************************************
!
      DRYF = WATER.LE.TINY
!
! *** FIND TOTALS *******************************************************
!
      IF (IPROB.EQ.0) THEN
         CONTINUE
      ELSE
         W(1) = WAER(1)
         W(2) = WAER(2)
         W(3) = WAER(3)
         W(4) = WAER(4)
         W(5) = WAER(5)
!
         IF (.NOT.DRYF) THEN
            W(3) = W(3) + GNH3
            W(4) = W(4) + GHNO3
            W(5) = W(5) + GHCL
         ENDIF
      ENDIF
!
      IF (LHOOK) CALL DR_HOOK('ISRPIA',1,ZHOOK_HANDLE)
!
! *** END OF SUBROUTINE ISRPIA *******************************************
!
      ENDSUBROUTINE ISRPIA
!=======================================================================
!
! *** ISORROPIA CODE
! *** SUBROUTINE PUSHERR
! *** THIS SUBROUTINE SAVES AN ERROR MESSAGE IN THE ERROR STACK
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!=======================================================================
!
      SUBROUTINE PUSHERR (IERR,ERRINF)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
      IMPLICIT NONE
      INCLUDE 'isrpia.inc'
      INTEGER :: IERR
      CHARACTER ERRINF*(*)
!
! *** SAVE ERROR CODE IF THERE IS ANY SPACE ***************************
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('PUSHERR',0,ZHOOK_HANDLE)
      IF (NOFER.LT.NERRMX) THEN
         NOFER         = NOFER + 1
         ERRSTK(NOFER) = IERR
         ERRMSG(NOFER) = ERRINF
         STKOFL        =.FALSE.
      ELSE
         STKOFL        =.TRUE.      ! STACK OVERFLOW
      ENDIF
!
! *** END OF SUBROUTINE PUSHERR ****************************************
!
      IF (LHOOK) CALL DR_HOOK('PUSHERR',1,ZHOOK_HANDLE)
      ENDSUBROUTINE PUSHERR



!=======================================================================
!
! *** ISORROPIA CODE
! *** SUBROUTINE ISERRINF
! *** THIS SUBROUTINE OBTAINS A COPY OF THE ERROR STACK (& MESSAGES)
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!=======================================================================
!
      SUBROUTINE ISERRINF (ERRSTKI, ERRMSGI, NOFERI, STKOFLI)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
      IMPLICIT NONE
      INCLUDE 'isrpia.inc'
      INTEGER :: I,NOFERI
      CHARACTER ERRMSGI*40
      INTEGER   ERRSTKI
      LOGICAL   STKOFLI
      DIMENSION ERRMSGI(NERRMX), ERRSTKI(NERRMX)
!
! *** OBTAIN WHOLE ERROR STACK ****************************************
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('ISERRINF',0,ZHOOK_HANDLE)
      DO 10 I=1,NOFER              ! Error messages & codes
        ERRSTKI(I) = ERRSTK(I)
        ERRMSGI(I) = ERRMSG(I)
  10  CONTINUE
!
      STKOFLI = STKOFL
      NOFERI  = NOFER
!
      IF (LHOOK) CALL DR_HOOK('ISERRINF',1,ZHOOK_HANDLE)
!
! *** END OF SUBROUTINE ISERRINF ***************************************
!
      ENDSUBROUTINE ISERRINF



!=======================================================================
!
! *** ISORROPIA CODE
! *** SUBROUTINE ISORINF
! *** THIS SUBROUTINE PROVIDES INFORMATION ABOUT ISORROPIA
!
! ======================== ARGUMENTS / USAGE ===========================
!
!  OUTPUT:
!  1. [VERSI]
!     CHARACTER*14 variable.
!     Contains version-date information of ISORROPIA
!
!  2. [NCMP]
!     INTEGER variable.
!     The number of components needed in input array WI
!     (or, the number of major species accounted for by ISORROPIA)
!
!  3. [NION]
!     INTEGER variable
!     The number of ions considered in the aqueous phase
!
!  4. [NAQGAS]
!     INTEGER variable
!     The number of undissociated species found in aqueous aerosol
!     phase
!
!  5. [NSOL]
!     INTEGER variable
!     The number of solids considered in the solid aerosol phase
!
!  6. [NERR]
!     INTEGER variable
!     The size of the error stack (maximum number of errors that can
!     be stored before the stack exhausts).
!
!  7. [TIN]
!     DOUBLE PRECISION variable
!     The value used for a very small number.
!
!  8. [GRT]
!     DOUBLE PRECISION variable
!     The value used for a very large number.
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!=======================================================================
!
      SUBROUTINE ISORINF (VERSI, NCMP, NION, NAQGAS, NSOL, NERR, TIN,   &
     &                    GRT)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
      IMPLICIT NONE
      INCLUDE 'isrpia.inc'
      REAL(KIND=8) :: TIN,GRT
      INTEGER :: NCMP,NION,NAQGAS,NSOL,NERR
      CHARACTER*14 VERSI
!
! *** ASSIGN INFO *******************************************************
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('ISORINF',0,ZHOOK_HANDLE)
      VERSI  = VERSION
      NCMP   = NCOMP
      NION   = NIONS
      NAQGAS = NGASAQ
      NSOL   = NSLDS
      NERR   = NERRMX
      TIN    = TINY
      GRT    = GREAT
!
      IF (LHOOK) CALL DR_HOOK('ISORINF',1,ZHOOK_HANDLE)
!
! *** END OF SUBROUTINE ISORINF *******************************************
!
      ENDSUBROUTINE ISORINF
