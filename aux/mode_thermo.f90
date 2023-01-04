!MNH_LIC Copyright 1994-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     #######################
      MODULE MODE_THERMO
!     #######################
!
!!****  *MODE_THERMO_MONO* -  module for routines SM_FOES,SM_PMR_HU
!!
!!    PURPOSE
!!    -------
!       The purpose of this executive module  is to package
!     the routine SM_FOES, SM_PMR_HU without use of comlib parallel routine
!
!
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!       NONE
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!      V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    28/08/94
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!  J.Escobar : 5/10/2018 : add FLUSH , for better logging in case of PB
!  P. Wautelet 10/04/2019: replace ABORT and STOP calls by Print_msg
!--------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!-------------------------------------------------------------------------------
USE MODE_MSG
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
IMPLICIT NONE

PRIVATE

PUBLIC :: DQSAT, DQSATI, QSAT, QSATI, SM_FOES, SM_PMR_HU

INTERFACE SM_FOES
  MODULE PROCEDURE SM_FOES_0D
  MODULE PROCEDURE SM_FOES_1D
  MODULE PROCEDURE SM_FOES_2D
  MODULE PROCEDURE SM_FOES_2D_MASK
  MODULE PROCEDURE SM_FOES_3D
END INTERFACE
INTERFACE QSAT
  MODULE PROCEDURE QSATW_3D
  MODULE PROCEDURE QSATW_2D
  MODULE PROCEDURE QSATW_2D_MASK
  MODULE PROCEDURE QSATW_1D
  MODULE PROCEDURE QSATW_0D
END INTERFACE
INTERFACE DQSAT
  MODULE PROCEDURE DQSATW_O_DT_2D_MASK
  MODULE PROCEDURE DQSATW_O_DT_1D
  MODULE PROCEDURE DQSATW_O_DT_3D
END INTERFACE
INTERFACE QSATI
  MODULE PROCEDURE QSATI_3D
  MODULE PROCEDURE QSATI_2D
  MODULE PROCEDURE QSATI_2D_MASK
  MODULE PROCEDURE QSATI_1D
  MODULE PROCEDURE QSATI_0D
END INTERFACE
INTERFACE DQSATI
  MODULE PROCEDURE DQSATI_O_DT_2D_MASK
  MODULE PROCEDURE DQSATI_O_DT_1D
  MODULE PROCEDURE DQSATI_O_DT_3D
END INTERFACE
INTERFACE SM_PMR_HU
  MODULE PROCEDURE SM_PMR_HU_1D
  MODULE PROCEDURE SM_PMR_HU_3D
END INTERFACE
CONTAINS
!-------------------------------------------------------------------------------
!     ####################################
      FUNCTION SM_FOES_3D(PT) RESULT(PFOES)
!     ####################################
!
!!****  *SM_FOES_3D * - function to compute saturation vapor pressure from
!!                    temperature
!!
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the saturation vapor
!     pressure from temperature
!
!
!!**  METHOD
!!    ------
!!       Given temperature T (PT), the saturation vapor pressure es(T)
!!    (FOES(PT)) is computed by integration of the Clapeyron equation
!!    from the triple point temperature Tt (XTT) and the saturation vapor
!!    pressure of the triple point es(Tt) (XESTT), i.e
!!
!!         es(T)= EXP( alphaw - betaw /T - gammaw Log(T) )
!!
!!     with :
!!       alphaw (XALPW) = LOG(es(Tt))+ betaw/Tt + gammaw Log(Tt)
!!       betaw (XBETAW) = Lv(Tt)/Rv + gammaw Tt
!!       gammaw (XGAMW) = (Cl -Cpv) /Rv
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST : comtains physical constants
!!        XALPW   : Constant for saturation vapor pressure function
!!        XBETAW  : Constant for saturation vapor pressure function
!!        XGAMW   : Constant for saturation vapor pressure function
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH
!!
!!
!!    AUTHOR
!!    ------
!!      V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    28/08/94
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments and results
!
!
REAL, DIMENSION(:,:,:), INTENT(IN)                :: PT     ! Temperature
                                                            ! (Kelvin)
REAL, DIMENSION(SIZE(PT,1),SIZE(PT,2),SIZE(PT,3)) :: PFOES  ! saturation vapor
                                                            ! pressure
                                                            ! (Pascal)
!
!*       0.2   Declarations of local variables
!
!-------------------------------------------------------------------------------
!
!*       1.    COMPUTE SATURATION VAPOR PRESSURE
!              ---------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_THERMO:SM_FOES_3D',0,ZHOOK_HANDLE)
PFOES(:,:,:) = EXP( XALPW - XBETAW/PT(:,:,:) - XGAMW*LOG(PT(:,:,:))  )
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_THERMO:SM_FOES_3D',1,ZHOOK_HANDLE)
END FUNCTION SM_FOES_3D
!     ####################################
      FUNCTION SM_FOES_1D(PT) RESULT(PFOES)
!     ####################################
!
!!****  *SM_FOES_1D * - function to compute saturation vapor pressure from
!!                    temperature
!!
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the saturation vapor
!     pressure from temperature
!
!
!!**  METHOD
!!    ------
!!       Given temperature T (PT), the saturation vapor pressure es(T)
!!    (FOES(PT)) is computed by integration of the Clapeyron equation
!!    from the triple point temperature Tt (XTT) and the saturation vapor
!!    pressure of the triple point es(Tt) (XESTT), i.e
!!
!!         es(T)= EXP( alphaw - betaw /T - gammaw Log(T) )
!!
!!     with :
!!       alphaw (XALPW) = LOG(es(Tt))+ betaw/Tt + gammaw Log(Tt)
!!       betaw (XBETAW) = Lv(Tt)/Rv + gammaw Tt
!!       gammaw (XGAMW) = (Cl -Cpv) /Rv
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST : comtains physical constants
!!        XALPW   : Constant for saturation vapor pressure function
!!        XBETAW  : Constant for saturation vapor pressure function
!!        XGAMW   : Constant for saturation vapor pressure function
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH
!!
!!
!!    AUTHOR
!!    ------
!!      V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    28/08/94
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments and results
!
!
REAL, DIMENSION(:), INTENT(IN)                :: PT     ! Temperature
                                                            ! (Kelvin)
REAL, DIMENSION(SIZE(PT)) :: PFOES  ! saturation vapor
                                                            ! pressure
                                                            ! (Pascal)
!
!*       0.2   Declarations of local variables
!
!-------------------------------------------------------------------------------
!
!*       1.    COMPUTE SATURATION VAPOR PRESSURE
!              ---------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_THERMO:SM_FOES_1D',0,ZHOOK_HANDLE)
PFOES(:) = EXP( XALPW - XBETAW/PT(:) - XGAMW*LOG(PT(:))  )
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_THERMO:SM_FOES_1D',1,ZHOOK_HANDLE)
END FUNCTION SM_FOES_1D
!-------------------------------------------------------------------------------
!     ####################################################
      FUNCTION SM_PMR_HU_3D(PP,PTV,PHU,PR,KITERMAX) RESULT(PMR)
!     ####################################################
!
!!****  *SM_PMR_HU_3D * - function to compute vapor mixing ratio
!!
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the vapor mixing ratio
!     from pressure, virtual  temperature and relative humidity
!
!
!!**  METHOD
!!    ------
!!      Given Pressure (PP), Virtual temperature (PTV) and Relative
!!    humidity (PHU), the vapor mixing ratio is computed by iterating
!!    the following procedure :
!!      T          ----> es(T)
!!      es(T) ,HU  ----> es(Td)
!!      es(Td), P  ----> r
!!      r , Tv     ----> T
!!
!!     at the beginning T=Tv
!!
!!    EXTERNAL
!!    --------
!!      SM_FOES   : to compute saturation vapor pressure
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST : comtains physical constants
!!        XRV : gas constant for vapor
!!        XRD : gas constant for dry air
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH
!!
!!
!!    AUTHOR
!!    ------
!!      V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       30/08/94
!!      Modification   16/03/95    remove the EPSILON function
!!      Modification   15/09/97    (V. Masson) add solid and liquid water phases
!!                                 in thetav computation
!!      Modification    22/01/2019 (P. Wautelet) use standard FLUSH statement
!!                                 instead of non standard intrinsics!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
USE MODD_LUNIT_n, ONLY: TLUOUT
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments and results
!
!
REAL, DIMENSION(:,:,:), INTENT(IN)               :: PP     ! Pressure
                                                           ! (Pa)
REAL, DIMENSION(:,:,:), INTENT(IN)               :: PTV    ! Virtual Temperature
                                                           ! (Kelvin)
REAL, DIMENSION(:,:,:), INTENT(IN)               :: PHU    ! Relative humidity
                                                           ! (percent)
REAL, DIMENSION(:,:,:,:), INTENT(IN)             :: PR     ! vapor, liquid and
                                                           ! solid water mixing
                                                           ! ratio

INTEGER,                  INTENT(IN), OPTIONAL   :: KITERMAX ! maximum number
                                                             ! of iterations
                                                             ! (default 10)
!
REAL, DIMENSION(SIZE(PP,1),SIZE(PP,2),SIZE(PP,3)) :: PMR   ! vapor mixing ratio
!
!*       0.2   Declarations of local variables
!
REAL, DIMENSION(SIZE(PP,1),SIZE(PP,2),SIZE(PP,3)) :: ZT      ! temperature
REAL, DIMENSION(SIZE(PP,1),SIZE(PP,2),SIZE(PP,3)) :: ZDT     ! increment of
                              ! temperature between two iterations
REAL                                              :: ZRDSRV  ! Rd/Rv
REAL, DIMENSION(SIZE(PP,1),SIZE(PP,2),SIZE(PP,3)) :: ZESTD   !  es(Td)
REAL, DIMENSION(SIZE(PP,1),SIZE(PP,2),SIZE(PP,3)) :: ZRSLW   ! total solid and liquid water mixing ratio
INTEGER                                           :: ITERMAX ! Maximum number
                                                             ! of iteration
INTEGER                                           :: ITER    ! iteration number of
REAL                                              :: ZEPS    ! a small number
INTEGER, DIMENSION(3)                             :: IMAXLOC ! localisation of
                                                             ! a maximum
INTEGER                                           :: ILUOUT
                                                             ! logical unit for
                                                             ! output-listing
                                                             ! and error code
INTEGER                                           :: JRR     ! loop counter
!-------------------------------------------------------------------------------
!
!*       1.    COMPUTE VAPOR MIXING RATIO
!              --------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_THERMO:SM_PMR_HU_3D',0,ZHOOK_HANDLE)
ITERMAX = 10
IF (PRESENT(KITERMAX)) ITERMAX=KITERMAX
ZRDSRV = XRD /XRV
ZEPS = XEPS_DT
!
ZRSLW(:,:,:)=0.
DO JRR=2,SIZE(PR,4)
  ZRSLW(:,:,:)=ZRSLW(:,:,:)+PR(:,:,:,JRR)
END DO
!
ZT(:,:,:) = PTV(:,:,:)
DO ITER=1,ITERMAX
  ZESTD(:,:,:) = PHU(:,:,:) * SM_FOES(ZT(:,:,:)) * 0.01
  PMR (:,:,:)  = ZRDSRV * ZESTD(:,:,:) /(PP(:,:,:) - ZESTD(:,:,:))
  ZDT(:,:,:)   = ZT(:,:,:)
  ZT(:,:,:)    = PTV(:,:,:) * (1.+PMR(:,:,:)+ZRSLW(:,:,:)) / (1.+ PMR(:,:,:)/ZRDSRV)
  ZDT(:,:,:)   = ABS(ZDT(:,:,:) - ZT(:,:,:))
END DO
!-------------------------------------------------------------------------------
!
!*       2.    NO CONVERGENCE
!              --------------
!
IF ( ANY(ZDT > ZEPS) ) THEN
  ILUOUT = TLUOUT%NLU
  WRITE(ILUOUT,*) 'ERROR IN FUNCTION SM_PMR_HU (module MODE_THERMO)'
  WRITE(ILUOUT,*) 'FUNCTION FAILS TO CONVERGE AFTER ', ITERMAX,' ITERATIONS'
  WRITE(ILUOUT,*) 'EPS = ' , ZEPS
  IMAXLOC(:) = MAXLOC(ZDT)
  WRITE(ILUOUT,*) 'MAXIMUM RESIDUAL DT :', MAXVAL(ZDT)
!  WRITE(ILUOUT,*) 'LOCATION OF THIS MAXIMUM I=',IMAXLOC(1),' J=',IMAXLOC(2), &
!                  ' K=',IMAXLOC(3)
  WRITE(ILUOUT,*) 'MR AT THIS MAXIMUM : ', PMR(IMAXLOC(1),IMAXLOC(2),IMAXLOC(3))
  WRITE(ILUOUT,*) 'T AT THIS MAXIMUM : ', ZT(IMAXLOC(1),IMAXLOC(2),IMAXLOC(3))
  WRITE(ILUOUT,*) 'JOB ABORTED '
  FLUSH(unit=ILUOUT)
  CALL PRINT_MSG( NVERB_FATAL, 'GEN', 'SM_PMR_HU_3D', 'failed to converge' )
END IF
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_THERMO:SM_PMR_HU_3D',1,ZHOOK_HANDLE)
END FUNCTION SM_PMR_HU_3D
!     ################################################################
      FUNCTION SM_PMR_HU_1D(PP,PTV,PHU,PR,KITERMAX) RESULT(PMR)
!     ################################################################
!
!!****  *SM_PMR_HU_1D * - function to compute vapor mixing ratio
!!
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the vapor mixing ratio
!     from pressure, virtual  temperature and relative humidity
!
!
!!**  METHOD
!!    ------
!!      Given Pressure (PP), Virtual temperature (PTV) and Relative
!!    humidity (PHU), the vapor mixing ratio is computed by iterating
!!    the following procedure :
!!      T          ----> es(T)
!!      es(T) ,HU  ----> es(Td)
!!      es(Td), P  ----> r
!!      r , Tv     ----> T
!!
!!     at the beginning T=Tv
!!
!!    EXTERNAL
!!    --------
!!      SM_FOES   : to compute saturation vapor pressure
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST : comtains physical constants
!!        XRV : gas constant for vapor
!!        XRD : gas constant for dry air
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH
!!
!!
!!    AUTHOR
!!    ------
!!      V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       30/08/94
!!      Modification   16/03/95    remove the EPSILON function
!!      Modification   15/09/97    (V. Masson) add solid and liquid water phases
!!                                 in thetav computation
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
USE MODD_LUNIT_n, ONLY: TLUOUT
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments and results
!
!
REAL, DIMENSION(:), INTENT(IN)               :: PP     ! Pressure
                                                       ! (Pa)
REAL, DIMENSION(:), INTENT(IN)               :: PTV    ! Virtual Temperature
                                                       ! (Kelvin)
REAL, DIMENSION(:), INTENT(IN)               :: PHU    ! Relative humidity
                                                       ! (percent)
REAL, DIMENSION(:,:), INTENT(IN)             :: PR     ! vapor, liquid and solid
                                                       ! water mixing ratio
INTEGER,              INTENT(IN), OPTIONAL   :: KITERMAX ! maximum number
                                                         ! of iterations
                                                         ! (default 10)
!
REAL, DIMENSION(SIZE(PP)) :: PMR   ! vapor mixing ratio
!
!*       0.2   Declarations of local variables
!
REAL, DIMENSION(SIZE(PP)) :: ZT      ! temperature
REAL, DIMENSION(SIZE(PP)) :: ZDT     ! increment of
                                     ! temperature between two iterations
REAL                                              :: ZRDSRV  ! Rd/Rv
REAL, DIMENSION(SIZE(PP)) :: ZESTD   ! es(Td)
REAL, DIMENSION(SIZE(PP)) :: ZRSLW   ! total solid and liquid water mixing ratio
INTEGER                                           :: ITERMAX ! Maximum number
                                                             ! of iteration
INTEGER                                           :: ITER    ! iteration number of
REAL                                              :: ZEPS    ! a small number
INTEGER,DIMENSION(1)                              :: IMAXLOC ! localisation of
                                                             ! a maximum
INTEGER                                           :: ILUOUT,IRESP
                                                             ! logical unit for
                                                             ! output-listing
                                                             ! and error code
!-------------------------------------------------------------------------------
!
!*       1.    COMPUTE VAPOR MIXING RATIO
!              --------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_THERMO:SM_PMR_HU_1D',0,ZHOOK_HANDLE)
ITERMAX = 10
IF (PRESENT(KITERMAX)) ITERMAX=KITERMAX
ZRDSRV = XRD /XRV
ZEPS = 1.E-5
!
IF (SIZE(PR,2)>1) THEN
  ZRSLW(:)=SUM(PR(:,2:),DIM=2)
ELSE
  ZRSLW(:)=0.
END IF
!
ZT(:) = PTV(:)
DO ITER=1,ITERMAX
  ZESTD(:) = PHU(:) * SM_FOES(ZT(:)) * 0.01
  PMR (:)  = ZRDSRV * ZESTD(:) /(PP(:) - ZESTD(:))
  ZDT(:)   = ZT(:)
  ZT(:)    = PTV(:) * (1.+PMR(:)+ZRSLW(:)) / (1.+ PMR(:)/ZRDSRV)
  ZDT(:)   = ABS(ZDT(:) - ZT(:))
END DO
!-------------------------------------------------------------------------------
!
!*       2.    NO CONVERGENCE
!              --------------
!
IF (ANY(ZDT>ZEPS)) THEN
  ILUOUT = TLUOUT%NLU
  WRITE(ILUOUT,*) 'ERROR IN FUNCTION SM_PMR_HU (module MODE_THERMO)'
  WRITE(ILUOUT,*) 'FUNCTION FAILS TO CONVERGE AFTER ', ITERMAX,' ITERATIONS'
  WRITE(ILUOUT,*) 'EPS = ' , ZEPS
  IMAXLOC = MAXLOC(ZDT)
  WRITE(ILUOUT,*) 'MAXIMUM RESIDUAL DT :', MAXVAL(ZDT)
  WRITE(ILUOUT,*) 'MR AT THIS MAXIMUM : ', PMR(IMAXLOC)
  WRITE(ILUOUT,*) 'T AT THIS MAXIMUM : ', ZT(IMAXLOC)
  WRITE(ILUOUT,*) 'JOB ABORTED '
  CALL PRINT_MSG( NVERB_FATAL, 'GEN', 'SM_PMR_HU_1D', 'failed to converge' )
END IF
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_THERMO:SM_PMR_HU_1D',1,ZHOOK_HANDLE)
END FUNCTION SM_PMR_HU_1D
!     ####################################
      FUNCTION SM_FOES_0D(PT) RESULT(PFOES)
!     ####################################
!
!!****  *SM_FOES_0D * - function to compute saturation vapor pressure from
!!                    temperature
!!
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the saturation vapor
!     pressure from temperature
!
!
!!**  METHOD
!!    ------
!!       Given temperature T (PT), the saturation vapor pressure es(T)
!!    (FOES(PT)) is computed by integration of the Clapeyron equation
!!    from the triple point temperature Tt (XTT) and the saturation vapor
!!    pressure of the triple point es(Tt) (XESTT), i.e
!!
!!         es(T)= EXP( alphaw - betaw /T - gammaw Log(T) )
!!
!!     with :
!!       alphaw (XALPW) = LOG(es(Tt))+ betaw/Tt + gammaw Log(Tt)
!!       betaw (XBETAW) = Lv(Tt)/Rv + gammaw Tt
!!       gammaw (XGAMW) = (Cl -Cpv) /Rv
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST : comtains physical constants
!!        XALPW   : Constant for saturation vapor pressure function
!!        XBETAW  : Constant for saturation vapor pressure function
!!        XGAMW   : Constant for saturation vapor pressure function
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH
!!
!!
!!    AUTHOR
!!    ------
!!      V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    28/08/94
!!                  24/12/97 (V. Masson) version for 0D arrays
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments and results
!
!
REAL,                 INTENT(IN)       :: PT     ! Temperature
                                                 ! (Kelvin)
REAL                                   :: PFOES  ! saturation vapor
                                                 ! pressure
                                                 ! (Pascal)
!
!*       0.2   Declarations of local variables
!
!-------------------------------------------------------------------------------
!
!*       1.    COMPUTE SATURATION VAPOR PRESSURE
!              ---------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_THERMO:SM_FOES_0D',0,ZHOOK_HANDLE)
PFOES = EXP( XALPW - XBETAW/PT - XGAMW*LOG(PT)  )
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_THERMO:SM_FOES_0D',1,ZHOOK_HANDLE)
END FUNCTION SM_FOES_0D
!
!-------------------------------------------------------------------------------
!     ####################################
      FUNCTION SM_FOES_2D(PT) RESULT(PFOES)
!     ####################################
!
!!****  *SM_FOES_2D * - function to compute saturation vapor pressure from
!!                    temperature
!!
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the saturation vapor
!     pressure from temperature
!
!
!!**  METHOD
!!    ------
!!       Given temperature T (PT), the saturation vapor pressure es(T)
!!    (FOES(PT)) is computed by integration of the Clapeyron equation
!!    from the triple point temperature Tt (XTT) and the saturation vapor
!!    pressure of the triple point es(Tt) (XESTT), i.e
!!
!!         es(T)= EXP( alphaw - betaw /T - gammaw Log(T) )
!!
!!     with :
!!       alphaw (XALPW) = LOG(es(Tt))+ betaw/Tt + gammaw Log(Tt)
!!       betaw (XBETAW) = Lv(Tt)/Rv + gammaw Tt
!!       gammaw (XGAMW) = (Cl -Cpv) /Rv
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST : comtains physical constants
!!        XALPW   : Constant for saturation vapor pressure function
!!        XBETAW  : Constant for saturation vapor pressure function
!!        XGAMW   : Constant for saturation vapor pressure function
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH
!!
!!
!!    AUTHOR
!!    ------
!!      V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    28/08/94
!!                  24/12/97 (V. Masson) version for 2D arrays
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments and results
!
!
REAL, DIMENSION(:,:), INTENT(IN)       :: PT     ! Temperature
                                                 ! (Kelvin)
REAL, DIMENSION(SIZE(PT,1),SIZE(PT,2)) :: PFOES  ! saturation vapor
                                                 ! pressure
                                                 ! (Pascal)
!
!*       0.2   Declarations of local variables
!
!-------------------------------------------------------------------------------
!
!*       1.    COMPUTE SATURATION VAPOR PRESSURE
!              ---------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_THERMO:SM_FOES_2D',0,ZHOOK_HANDLE)
PFOES(:,:) = EXP( XALPW - XBETAW/PT(:,:) - XGAMW*LOG(PT(:,:))  )
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_THERMO:SM_FOES_2D',1,ZHOOK_HANDLE)
END FUNCTION SM_FOES_2D
!
!-------------------------------------------------------------------------------
!
!     ################################################
      FUNCTION SM_FOES_2D_MASK(OMASK,PT) RESULT(PFOES)
!     ################################################
!
!!****  *SM_FOES_2D * - function to compute saturation vapor pressure from
!!                    temperature
!!
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the saturation vapor
!     pressure from temperature
!
!
!!**  METHOD
!!    ------
!!       Given temperature T (PT), the saturation vapor pressure es(T)
!!    (FOES(PT)) is computed by integration of the Clapeyron equation
!!    from the triple point temperature Tt (XTT) and the saturation vapor
!!    pressure of the triple point es(Tt) (XESTT), i.e
!!
!!         es(T)= EXP( alphaw - betaw /T - gammaw Log(T) )
!!
!!     with :
!!       alphaw (XALPW) = LOG(es(Tt))+ betaw/Tt + gammaw Log(Tt)
!!       betaw (XBETAW) = Lv(Tt)/Rv + gammaw Tt
!!       gammaw (XGAMW) = (Cl -Cpv) /Rv
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST : comtains physical constants
!!        XALPW   : Constant for saturation vapor pressure function
!!        XBETAW  : Constant for saturation vapor pressure function
!!        XGAMW   : Constant for saturation vapor pressure function
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH
!!
!!
!!    AUTHOR
!!    ------
!!      V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    28/08/94
!!                  24/12/97 (V. Masson) version for 2D arrays
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments and results
!
!
LOGICAL, DIMENSION(:,:), INTENT(IN)    :: OMASK  ! Localization mask
REAL, DIMENSION(:,:), INTENT(IN)       :: PT     ! Temperature
                                                 ! (Kelvin)
REAL, DIMENSION(SIZE(PT,1),SIZE(PT,2)) :: PFOES  ! saturation vapor
                                                 ! pressure
                                                 ! (Pascal)
!
!*       0.2   Declarations of local variables
!
!-------------------------------------------------------------------------------
!
!*       1.    COMPUTE SATURATION VAPOR PRESSURE
!              ---------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_THERMO:SM_FOES_2D_MASK',0,ZHOOK_HANDLE)
WHERE (OMASK(:,:))
  PFOES(:,:) = EXP( XALPW - XBETAW/PT(:,:) - XGAMW*LOG(PT(:,:))  )
END WHERE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_THERMO:SM_FOES_2D_MASK',1,ZHOOK_HANDLE)
END FUNCTION SM_FOES_2D_MASK
!
!-------------------------------------------------------------------------------
!
!     ######################################
      FUNCTION QSATW_3D(PT,PP) RESULT(PQSAT)
!     ######################################
!
!!****  *QSATW * - function to compute saturation vapor humidity from
!!                 temperature
!!
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the saturation vapor
!     pressure from temperature
!
!
!!**  METHOD
!!    ------
!!       Given temperature T (PT), the saturation vapor pressure es(T)
!!    (FOES(PT)) is computed by integration of the Clapeyron equation
!!    from the triple point temperature Tt (XTT) and the saturation vapor
!!    pressure of the triple point es(Tt) (XESTT), i.e
!!
!!         es(T)= EXP( alphaw - betaw /T - gammaw Log(T) )
!!
!!     with :
!!       alphaw (XALPW) = LOG(es(Tt))+ betaw/Tt + gammaw Log(Tt)
!!       betaw (XBETAW) = Lv(Tt)/Rv + gammaw Tt
!!       gammaw (XGAMW) = (Cl -Cpv) /Rv
!!
!!      Then, the specific humidity at saturation is deduced.
!!
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST : comtains physical constants
!!        XALPW   : Constant for saturation vapor pressure function
!!        XBETAW  : Constant for saturation vapor pressure function
!!        XGAMW   : Constant for saturation vapor pressure function
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH
!!
!!
!!    AUTHOR
!!    ------
!!      V. Masson       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    21/09/98
!!      S. Riette   april 2011 : protection in high statosphere where ZFOES > PP
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments and results
!
!
REAL, DIMENSION(:,:,:), INTENT(IN)                :: PT     ! Temperature
                                                            ! (Kelvin)
REAL, DIMENSION(:,:,:), INTENT(IN)                :: PP     ! Pressure
                                                            ! (Pa)
REAL, DIMENSION(SIZE(PT,1),SIZE(PT,2),SIZE(PT,3)) :: PQSAT  ! saturation vapor
                                                            ! specific humidity
                                                            ! with respect to
                                                            ! water (kg/kg)
!
!*       0.2   Declarations of local variables
!
REAL, DIMENSION(SIZE(PT,1),SIZE(PT,2),SIZE(PT,3)) :: ZFOES  ! saturation vapor
                                                            ! pressure
                                                            ! (Pascal)
!
!-------------------------------------------------------------------------------
!
!*       1.    COMPUTE SATURATION VAPOR PRESSURE
!              ---------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_THERMO:QSATW_3D',0,ZHOOK_HANDLE)
ZFOES(:,:,:) = MIN(EXP( XALPW - XBETAW/PT(:,:,:) - XGAMW*LOG(PT(:,:,:))  ), 0.99*PP(:,:,:))
!
!*       2.    COMPUTE SATURATION HUMIDITY
!              ---------------------------
!
PQSAT(:,:,:) = XRD/XRV*ZFOES(:,:,:)/PP(:,:,:)   &
                   / (1.+(XRD/XRV-1.)*ZFOES(:,:,:)/PP(:,:,:))
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_THERMO:QSATW_3D',1,ZHOOK_HANDLE)
END FUNCTION QSATW_3D
!
!-------------------------------------------------------------------------------
!
!     ######################################
      FUNCTION QSATW_2D(PT,PP) RESULT(PQSAT)
!     ######################################
!
!!****  *QSATW * - function to compute saturation vapor humidity from
!!                 temperature
!!
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the saturation vapor
!     pressure from temperature
!
!
!!**  METHOD
!!    ------
!!       Given temperature T (PT), the saturation vapor pressure es(T)
!!    (FOES(PT)) is computed by integration of the Clapeyron equation
!!    from the triple point temperature Tt (XTT) and the saturation vapor
!!    pressure of the triple point es(Tt) (XESTT), i.e
!!
!!         es(T)= EXP( alphaw - betaw /T - gammaw Log(T) )
!!
!!     with :
!!       alphaw (XALPW) = LOG(es(Tt))+ betaw/Tt + gammaw Log(Tt)
!!       betaw (XBETAW) = Lv(Tt)/Rv + gammaw Tt
!!       gammaw (XGAMW) = (Cl -Cpv) /Rv
!!
!!      Then, the specific humidity at saturation is deduced.
!!
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST : comtains physical constants
!!        XALPW   : Constant for saturation vapor pressure function
!!        XBETAW  : Constant for saturation vapor pressure function
!!        XGAMW   : Constant for saturation vapor pressure function
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH
!!
!!
!!    AUTHOR
!!    ------
!!      V. Masson       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    21/09/98
!!      S. Riette   april 2011 : protection in high statosphere where ZFOES > PP
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments and results
!
!
REAL, DIMENSION(:,:), INTENT(IN)                :: PT     ! Temperature
                                                          ! (Kelvin)
REAL, DIMENSION(:,:), INTENT(IN)                :: PP     ! Pressure
                                                          ! (Pa)
REAL, DIMENSION(SIZE(PT,1),SIZE(PT,2))          :: PQSAT  ! saturation vapor
                                                          ! specific humidity
                                                          ! with respect to
                                                          ! water (kg/kg)
!
!*       0.2   Declarations of local variables
!
REAL, DIMENSION(SIZE(PT,1),SIZE(PT,2))          :: ZFOES  ! saturation vapor
                                                          ! pressure
                                                          ! (Pascal)
!
!-------------------------------------------------------------------------------
!
!*       1.    COMPUTE SATURATION VAPOR PRESSURE
!              ---------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_THERMO:QSATW_2D',0,ZHOOK_HANDLE)
ZFOES(:,:) = MIN(EXP( XALPW - XBETAW/PT(:,:) - XGAMW*LOG(PT(:,:))  ), 0.99*PP(:,:))
!
!*       2.    COMPUTE SATURATION HUMIDITY
!              ---------------------------
!
PQSAT(:,:) = XRD/XRV*ZFOES(:,:)/PP(:,:)   &
                   / (1.+(XRD/XRV-1.)*ZFOES(:,:)/PP(:,:))
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_THERMO:QSATW_2D',1,ZHOOK_HANDLE)
END FUNCTION QSATW_2D
!
!-------------------------------------------------------------------------------
!
!     #################################################
      FUNCTION QSATW_2D_MASK(OMASK,PT,PP) RESULT(PQSAT)
!     #################################################
!
!!****  *QSATW * - function to compute saturation vapor humidity from
!!                 temperature
!!
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the saturation vapor
!     pressure from temperature
!
!
!!**  METHOD
!!    ------
!!       Given temperature T (PT), the saturation vapor pressure es(T)
!!    (FOES(PT)) is computed by integration of the Clapeyron equation
!!    from the triple point temperature Tt (XTT) and the saturation vapor
!!    pressure of the triple point es(Tt) (XESTT), i.e
!!
!!         es(T)= EXP( alphaw - betaw /T - gammaw Log(T) )
!!
!!     with :
!!       alphaw (XALPW) = LOG(es(Tt))+ betaw/Tt + gammaw Log(Tt)
!!       betaw (XBETAW) = Lv(Tt)/Rv + gammaw Tt
!!       gammaw (XGAMW) = (Cl -Cpv) /Rv
!!
!!      Then, the specific humidity at saturation is deduced.
!!
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST : comtains physical constants
!!        XALPW   : Constant for saturation vapor pressure function
!!        XBETAW  : Constant for saturation vapor pressure function
!!        XGAMW   : Constant for saturation vapor pressure function
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH
!!
!!
!!    AUTHOR
!!    ------
!!      V. Masson       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    21/09/98
!!      S. Riette   april 2011 : protection in high statosphere where ZFOES > PP
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments and results
!
!
LOGICAL, DIMENSION(:,:), INTENT(IN)             :: OMASK  ! Localization mask
REAL, DIMENSION(:,:), INTENT(IN)                :: PT     ! Temperature
                                                          ! (Kelvin)
REAL, DIMENSION(:,:), INTENT(IN)                :: PP     ! Pressure
                                                          ! (Pa)
REAL, DIMENSION(SIZE(PT,1),SIZE(PT,2))          :: PQSAT  ! saturation vapor
                                                          ! specific humidity
                                                          ! with respect to
                                                          ! water (kg/kg)
!
!*       0.2   Declarations of local variables
!
REAL, DIMENSION(SIZE(PT,1),SIZE(PT,2))          :: ZFOES  ! saturation vapor
                                                          ! pressure
                                                          ! (Pascal)
!
!-------------------------------------------------------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_THERMO:QSATW_2D_MASK',0,ZHOOK_HANDLE)
WHERE (OMASK(:,:))
!
!*       1.    COMPUTE SATURATION VAPOR PRESSURE
!              ---------------------------------
!
  ZFOES(:,:) = MIN(EXP( XALPW - XBETAW/PT(:,:) - XGAMW*LOG(PT(:,:))  ), 0.99*PP(:,:))
!
!*       2.    COMPUTE SATURATION HUMIDITY
!              ---------------------------
!
  PQSAT(:,:) = XRD/XRV*ZFOES(:,:)/PP(:,:)   &
                     / (1.+(XRD/XRV-1.)*ZFOES(:,:)/PP(:,:))
ELSEWHERE
!
!*       3.    BOGUS VALUE
!              -----------
!
  PQSAT(:,:) = 0.
END WHERE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_THERMO:QSATW_2D_MASK',1,ZHOOK_HANDLE)
END FUNCTION QSATW_2D_MASK
!
!-------------------------------------------------------------------------------
!
!     ######################################
      FUNCTION QSATW_1D(PT,PP) RESULT(PQSAT)
!     ######################################
!
!!****  *QSATW * - function to compute saturation vapor humidity from
!!                 temperature
!!
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the saturation vapor
!     pressure from temperature
!
!
!!**  METHOD
!!    ------
!!       Given temperature T (PT), the saturation vapor pressure es(T)
!!    (FOES(PT)) is computed by integration of the Clapeyron equation
!!    from the triple point temperature Tt (XTT) and the saturation vapor
!!    pressure of the triple point es(Tt) (XESTT), i.e
!!
!!         es(T)= EXP( alphaw - betaw /T - gammaw Log(T) )
!!
!!     with :
!!       alphaw (XALPW) = LOG(es(Tt))+ betaw/Tt + gammaw Log(Tt)
!!       betaw (XBETAW) = Lv(Tt)/Rv + gammaw Tt
!!       gammaw (XGAMW) = (Cl -Cpv) /Rv
!!
!!      Then, the specific humidity at saturation is deduced.
!!
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST : comtains physical constants
!!        XALPW   : Constant for saturation vapor pressure function
!!        XBETAW  : Constant for saturation vapor pressure function
!!        XGAMW   : Constant for saturation vapor pressure function
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH
!!
!!
!!    AUTHOR
!!    ------
!!      V. Masson       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    21/09/98
!!      S. Riette   april 2011 : protection in high statosphere where ZFOES > PP
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments and results
!
!
REAL, DIMENSION(:), INTENT(IN)                :: PT     ! Temperature
                                                        ! (Kelvin)
REAL, DIMENSION(:), INTENT(IN)                :: PP     ! Pressure
                                                        ! (Pa)
REAL, DIMENSION(SIZE(PT,1))                   :: PQSAT  ! saturation vapor
                                                        ! specific humidity
                                                        ! with respect to
                                                        ! water (kg/kg)
!
!*       0.2   Declarations of local variables
!
REAL, DIMENSION(SIZE(PT,1))                   :: ZFOES  ! saturation vapor
                                                        ! pressure
                                                        ! (Pascal)
!
!-------------------------------------------------------------------------------
!
!*       1.    COMPUTE SATURATION VAPOR PRESSURE
!              ---------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_THERMO:QSATW_1D',0,ZHOOK_HANDLE)
ZFOES(:) = MIN(EXP( XALPW - XBETAW/PT(:) - XGAMW*LOG(PT(:))  ), 0.99*PP(:))
!
!*       2.    COMPUTE SATURATION HUMIDITY
!              ---------------------------
!
PQSAT(:) = XRD/XRV*ZFOES(:)/PP(:)   &
                   / (1.+(XRD/XRV-1.)*ZFOES(:)/PP(:))
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_THERMO:QSATW_1D',1,ZHOOK_HANDLE)
END FUNCTION QSATW_1D
!
!-------------------------------------------------------------------------------
!
!     ######################################
      FUNCTION QSATW_0D(PT,PP) RESULT(PQSAT)
!     ######################################
!
!!****  *QSATW * - function to compute saturation vapor humidity from
!!                 temperature
!!
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the saturation vapor
!     pressure from temperature
!
!
!!**  METHOD
!!    ------
!!       Given temperature T (PT), the saturation vapor pressure es(T)
!!    (FOES(PT)) is computed by integration of the Clapeyron equation
!!    from the triple point temperature Tt (XTT) and the saturation vapor
!!    pressure of the triple point es(Tt) (XESTT), i.e
!!
!!         es(T)= EXP( alphaw - betaw /T - gammaw Log(T) )
!!
!!     with :
!!       alphaw (XALPW) = LOG(es(Tt))+ betaw/Tt + gammaw Log(Tt)
!!       betaw (XBETAW) = Lv(Tt)/Rv + gammaw Tt
!!       gammaw (XGAMW) = (Cl -Cpv) /Rv
!!
!!      Then, the specific humidity at saturation is deduced.
!!
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST : comtains physical constants
!!        XALPW   : Constant for saturation vapor pressure function
!!        XBETAW  : Constant for saturation vapor pressure function
!!        XGAMW   : Constant for saturation vapor pressure function
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH
!!
!!
!!    AUTHOR
!!    ------
!!      V. Masson       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    21/09/98
!!      S. Riette   april 2011 : protection in high statosphere where ZFOES > PP
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments and results
!
!
REAL, INTENT(IN)                :: PT     ! Temperature
                                          ! (Kelvin)
REAL, INTENT(IN)                :: PP     ! Pressure
                                          ! (Pa)
REAL                            :: PQSAT  ! saturation vapor
                                          ! specific humidity
                                          ! with respect to
                                          ! water (kg/kg)
!
!*       0.2   Declarations of local variables
!
REAL                            :: ZFOES  ! saturation vapor
                                          ! pressure
                                          ! (Pascal)
!
!-------------------------------------------------------------------------------
!
!*       1.    COMPUTE SATURATION VAPOR PRESSURE
!              ---------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_THERMO:QSATW_0D',0,ZHOOK_HANDLE)
ZFOES = MIN(EXP( XALPW - XBETAW/PT - XGAMW*LOG(PT)  ), 0.99*PP)
!
!*       2.    COMPUTE SATURATION HUMIDITY
!              ---------------------------
!
PQSAT = XRD/XRV*ZFOES/PP / (1.+(XRD/XRV-1.)*ZFOES/PP)
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_THERMO:QSATW_0D',1,ZHOOK_HANDLE)
END FUNCTION QSATW_0D
!
!-------------------------------------------------------------------------------
!
!     ##############################################################
      FUNCTION DQSATW_O_DT_2D_MASK(OMASK,PT,PP,PQSAT) RESULT(PDQSAT)
!     ##############################################################
!
!!****  *QSATW * - function to compute saturation vapor humidity from
!!                 temperature
!!
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the saturation vapor
!     pressure from temperature
!
!
!!**  METHOD
!!    ------
!!       Given temperature T (PT), the saturation vapor pressure es(T)
!!    (FOES(PT)) is computed by integration of the Clapeyron equation
!!    from the triple point temperature Tt (XTT) and the saturation vapor
!!    pressure of the triple point es(Tt) (XESTT), i.e
!!
!!         es(T)= EXP( alphaw - betaw /T - gammaw Log(T) )
!!
!!     with :
!!       alphaw (XALPW) = LOG(es(Tt))+ betaw/Tt + gammaw Log(Tt)
!!       betaw (XBETAW) = Lv(Tt)/Rv + gammaw Tt
!!       gammaw (XGAMW) = (Cl -Cpv) /Rv
!!
!!      Then, the specific humidity at saturation is deduced.
!!
!!      Finally, dqsat / dT  (T) is computed.
!!
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST : comtains physical constants
!!        XALPW   : Constant for saturation vapor pressure function
!!        XBETAW  : Constant for saturation vapor pressure function
!!        XGAMW   : Constant for saturation vapor pressure function
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH
!!
!!
!!    AUTHOR
!!    ------
!!      V. Masson       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    21/09/98
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments and results
!
!
LOGICAL, DIMENSION(:,:), INTENT(IN)             :: OMASK  ! Localization mask
REAL,    DIMENSION(:,:), INTENT(IN)             :: PT     ! Temperature
                                                          ! (Kelvin)
REAL,    DIMENSION(:,:), INTENT(IN)             :: PP     ! Pressure
                                                          ! (Pa)
REAL,    DIMENSION(:,:), INTENT(IN)             :: PQSAT  ! saturation vapor
                                                          ! specific humidity
                                                          ! with respect to
                                                          ! water (kg/kg))
REAL,    DIMENSION(SIZE(PT,1),SIZE(PT,2))       :: PDQSAT ! derivative according
                                                          ! to temperature of
                                                          ! saturation vapor
                                                          ! specific humidity
                                                          ! with respect to
                                                          ! water (kg/kg))
!
!*       0.2   Declarations of local variables
!
REAL, DIMENSION(SIZE(PT,1),SIZE(PT,2))          :: ZFOES  ! saturation vapor
                                                          ! pressure
                                                          ! (Pascal)
!
!-------------------------------------------------------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_THERMO:DQSATW_O_DT_2D_MASK',0,ZHOOK_HANDLE)
WHERE (OMASK(:,:))
!
!*       1.    COMPUTE SATURATION VAPOR PRESSURE
!              ---------------------------------
!
  ZFOES(:,:) = PP(:,:) / (1.+XRD/XRV*(1./PQSAT(:,:)-1.))
!
!*       2.    DERIVATION ACCORDING TO TEMPERATURE
!              -----------------------------------
!
  PDQSAT(:,:) = PQSAT(:,:) / (1.+(XRD/XRV-1.)*ZFOES(:,:)/PP(:,:) ) &
                   * (XBETAW/PT(:,:)**2 - XGAMW/PT(:,:))
ELSEWHERE
!
!*       3.    BOGUS VALUE
!              -----------
!
  PDQSAT(:,:) = 0.
END WHERE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_THERMO:DQSATW_O_DT_2D_MASK',1,ZHOOK_HANDLE)
END FUNCTION DQSATW_O_DT_2D_MASK
!
!-------------------------------------------------------------------------------
!     ##############################################################
      FUNCTION DQSATW_O_DT_1D(PT,PP,PQSAT) RESULT(PDQSAT)
!     ##############################################################
!
!!****  *QSATW * - function to compute saturation vapor humidity from
!!                 temperature
!!
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the saturation vapor
!     pressure from temperature
!
!
!!**  METHOD
!!    ------
!!       Given temperature T (PT), the saturation vapor pressure es(T)
!!    (FOES(PT)) is computed by integration of the Clapeyron equation
!!    from the triple point temperature Tt (XTT) and the saturation vapor
!!    pressure of the triple point es(Tt) (XESTT), i.e
!!
!!         es(T)= EXP( alphaw - betaw /T - gammaw Log(T) )
!!
!!     with :
!!       alphaw (XALPW) = LOG(es(Tt))+ betaw/Tt + gammaw Log(Tt)
!!       betaw (XBETAW) = Lv(Tt)/Rv + gammaw Tt
!!       gammaw (XGAMW) = (Cl -Cpv) /Rv
!!
!!      Then, the specific humidity at saturation is deduced.
!!
!!      Finally, dqsat / dT  (T) is computed.
!!
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST : comtains physical constants
!!        XALPW   : Constant for saturation vapor pressure function
!!        XBETAW  : Constant for saturation vapor pressure function
!!        XGAMW   : Constant for saturation vapor pressure function
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH
!!
!!
!!    AUTHOR
!!    ------
!!      V. Masson       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    21/09/98
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments and results
!
!
REAL,    DIMENSION(:), INTENT(IN)             :: PT     ! Temperature
                                                          ! (Kelvin)
REAL,    DIMENSION(:), INTENT(IN)               :: PP     ! Pressure
                                                          ! (Pa)
REAL,    DIMENSION(:), INTENT(IN)               :: PQSAT  ! saturation vapor
                                                          ! specific humidity
                                                          ! with respect to
                                                          ! water (kg/kg))
REAL,    DIMENSION(SIZE(PT))                    :: PDQSAT ! derivative according
                                                          ! to temperature of
                                                          ! saturation vapor
                                                          ! specific humidity
                                                          ! with respect to
                                                          ! water (kg/kg))
!
!*       0.2   Declarations of local variables
!
REAL, DIMENSION(SIZE(PT))                       :: ZFOES  ! saturation vapor
                                                          ! pressure
                                                          ! (Pascal)
!
!-------------------------------------------------------------------------------
!
!
!*       1.    COMPUTE SATURATION VAPOR PRESSURE
!              ---------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_THERMO:DQSATW_O_DT_1D',0,ZHOOK_HANDLE)
ZFOES(:) = PP(:) / (1.+XRD/XRV*(1./PQSAT(:)-1.))
!
!*       2.    DERIVATION ACCORDING TO TEMPERATURE
!              -----------------------------------
!
PDQSAT(:) = PQSAT(:) / (1.+(XRD/XRV-1.)*ZFOES(:)/PP(:) ) &
                   * (XBETAW/PT(:)**2 - XGAMW/PT(:))
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_THERMO:DQSATW_O_DT_1D',1,ZHOOK_HANDLE)
END FUNCTION DQSATW_O_DT_1D
!
!-------------------------------------------------------------------------------
!     ##############################################################
      FUNCTION DQSATW_O_DT_3D(PT,PP,PQSAT) RESULT(PDQSAT)
!     ##############################################################
!
!!****  *QSATW * - function to compute saturation vapor humidity from
!!                 temperature
!!
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the saturation vapor
!     pressure from temperature
!
!
!!**  METHOD
!!    ------
!!       Given temperature T (PT), the saturation vapor pressure es(T)
!!    (FOES(PT)) is computed by integration of the Clapeyron equation
!!    from the triple point temperature Tt (XTT) and the saturation vapor
!!    pressure of the triple point es(Tt) (XESTT), i.e
!!
!!         es(T)= EXP( alphaw - betaw /T - gammaw Log(T) )
!!
!!     with :
!!       alphaw (XALPW) = LOG(es(Tt))+ betaw/Tt + gammaw Log(Tt)
!!       betaw (XBETAW) = Lv(Tt)/Rv + gammaw Tt
!!       gammaw (XGAMW) = (Cl -Cpv) /Rv
!!
!!      Then, the specific humidity at saturation is deduced.
!!
!!      Finally, dqsat / dT  (T) is computed.
!!
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST : comtains physical constants
!!        XALPW   : Constant for saturation vapor pressure function
!!        XBETAW  : Constant for saturation vapor pressure function
!!        XGAMW   : Constant for saturation vapor pressure function
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH
!!
!!
!!    AUTHOR
!!    ------
!!      V. Masson       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    21/09/98
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments and results
!
!
REAL,    DIMENSION(:,:,:), INTENT(IN)             :: PT     ! Temperature
                                                          ! (Kelvin)
REAL,    DIMENSION(:,:,:), INTENT(IN)               :: PP     ! Pressure
                                                          ! (Pa)
REAL,    DIMENSION(:,:,:), INTENT(IN)               :: PQSAT  ! saturation vapor
                                                          ! specific humidity
                                                          ! with respect to
                                                          ! water (kg/kg))
REAL,    DIMENSION(SIZE(PT,1),SIZE(PT,2),SIZE(PT,3))   :: PDQSAT ! derivative according
                                                          ! to temperature of
                                                          ! saturation vapor
                                                          ! specific humidity
                                                          ! with respect to
                                                          ! water (kg/kg))
!
!*       0.2   Declarations of local variables
!
REAL, DIMENSION(SIZE(PT,1),SIZE(PT,2),SIZE(PT,3))   :: ZFOES  ! saturation vapor
                                                          ! pressure
                                                          ! (Pascal)
!
!-------------------------------------------------------------------------------
!
!
!*       1.    COMPUTE SATURATION VAPOR PRESSURE
!              ---------------------------------
!
ZFOES(:,:,:) = PP(:,:,:) / (1.+XRD/XRV*(1./PQSAT(:,:,:)-1.))
!
!*       2.    DERIVATION ACCORDING TO TEMPERATURE
!              -----------------------------------
!
PDQSAT(:,:,:) = PQSAT(:,:,:) / (1.+(XRD/XRV-1.)*ZFOES(:,:,:)/PP(:,:,:) ) &
                   * (XBETAW/PT(:,:,:)**2 - XGAMW/PT(:,:,:))
!
!-------------------------------------------------------------------------------
!
END FUNCTION DQSATW_O_DT_3D
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
!     ##############################################################
      FUNCTION DQSATI_O_DT_2D_MASK(OMASK,PT,PP,PQSAT) RESULT(PDQSAT)
!     ##############################################################
!
!!****  *QSATW * - function to compute saturation vapor humidity from
!!                 temperature (with respect to ice)
!!
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the saturation vapor
!     pressure from temperature
!
!
!!**  METHOD
!!    ------
!!       Given temperature T (PT), the saturation vapor pressure es(T)
!!    (FOES(PT)) is computed by integration of the Clapeyron equation
!!    from the triple point temperature Tt (XTT) and the saturation vapor
!!    pressure of the triple point es(Tt) (XESTT), i.e
!!
!!         es(T)= EXP( alphaw - betaw /T - gammaw Log(T) )
!!
!!     with :
!!       alphaw (XALPW) = LOG(es(Tt))+ betaw/Tt + gammaw Log(Tt)
!!       betaw (XBETAW) = Lv(Tt)/Rv + gammaw Tt
!!       gammaw (XGAMW) = (Cl -Cpv) /Rv
!!
!!      Then, the specific humidity at saturation is deduced.
!!
!!      Finally, dqsat / dT  (T) is computed.
!!
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST : comtains physical constants
!!        XALPW   : Constant for saturation vapor pressure function
!!        XBETAW  : Constant for saturation vapor pressure function
!!        XGAMW   : Constant for saturation vapor pressure function
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH
!!
!!
!!    AUTHOR
!!    ------
!!      V. Masson       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    21/09/98
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments and results
!
!
LOGICAL, DIMENSION(:,:), INTENT(IN)             :: OMASK  ! Localization mask
REAL,    DIMENSION(:,:), INTENT(IN)             :: PT     ! Temperature
                                                          ! (Kelvin)
REAL,    DIMENSION(:,:), INTENT(IN)             :: PP     ! Pressure
                                                          ! (Pa)
REAL,    DIMENSION(:,:), INTENT(IN)             :: PQSAT  ! saturation vapor
                                                          ! specific humidity
                                                          ! with respect to
                                                          ! water (kg/kg))
REAL,    DIMENSION(SIZE(PT,1),SIZE(PT,2))       :: PDQSAT ! derivative according
                                                          ! to temperature of
                                                          ! saturation vapor
                                                          ! specific humidity
                                                          ! with respect to
                                                          ! water (kg/kg))
!
!*       0.2   Declarations of local variables
!
REAL, DIMENSION(SIZE(PT,1),SIZE(PT,2))          :: ZFOES  ! saturation vapor
                                                          ! pressure
                                                          ! (Pascal)
!
!-------------------------------------------------------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_THERMO:DQSATI_O_DT_2D_MASK',0,ZHOOK_HANDLE)
WHERE (OMASK(:,:))
!
!*       1.    COMPUTE SATURATION VAPOR PRESSURE
!              ---------------------------------
!
  ZFOES(:,:) = PP(:,:) / (1.+XRD/XRV*(1./PQSAT(:,:)-1.))
!
!*       3.    DERIVATION ACCORDING TO TEMPERATURE
!              -----------------------------------
!
  PDQSAT(:,:) = PQSAT(:,:) / (1.+(XRD/XRV-1.)*ZFOES(:,:)/PP(:,:) ) &
                   * (XBETAI/PT(:,:)**2 - XGAMI/PT(:,:))
ELSEWHERE
!
!*       3.    BOGUS VALUE
!              -----------
!
  PDQSAT(:,:) = 0.
END WHERE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_THERMO:DQSATI_O_DT_2D_MASK',1,ZHOOK_HANDLE)
END FUNCTION DQSATI_O_DT_2D_MASK
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!     ##############################################################
      FUNCTION DQSATI_O_DT_1D(PT,PP,PQSAT) RESULT(PDQSAT)
!     ##############################################################
!
!!****  *QSATW * - function to compute saturation vapor humidity from
!!                 temperature (with respect to ice)
!!
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the saturation vapor
!     pressure from temperature
!
!
!!**  METHOD
!!    ------
!!       Given temperature T (PT), the saturation vapor pressure es(T)
!!    (FOES(PT)) is computed by integration of the Clapeyron equation
!!    from the triple point temperature Tt (XTT) and the saturation vapor
!!    pressure of the triple point es(Tt) (XESTT), i.e
!!
!!         es(T)= EXP( alphaw - betaw /T - gammaw Log(T) )
!!
!!     with :
!!       alphaw (XALPW) = LOG(es(Tt))+ betaw/Tt + gammaw Log(Tt)
!!       betaw (XBETAW) = Lv(Tt)/Rv + gammaw Tt
!!       gammaw (XGAMW) = (Cl -Cpv) /Rv
!!
!!      Then, the specific humidity at saturation is deduced.
!!
!!      Finally, dqsat / dT  (T) is computed.
!!
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST : comtains physical constants
!!        XALPW   : Constant for saturation vapor pressure function
!!        XBETAW  : Constant for saturation vapor pressure function
!!        XGAMW   : Constant for saturation vapor pressure function
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH
!!
!!
!!    AUTHOR
!!    ------
!!      V. Masson       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    21/09/98
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments and results
!
!
REAL,    DIMENSION(:), INTENT(IN)               :: PT     ! Temperature
                                                          ! (Kelvin)
REAL,    DIMENSION(:), INTENT(IN)               :: PP     ! Pressure
                                                          ! (Pa)
REAL,    DIMENSION(:), INTENT(IN)               :: PQSAT  ! saturation vapor
                                                          ! specific humidity
                                                          ! with respect to
                                                          ! water (kg/kg))
REAL,    DIMENSION(SIZE(PT))                    :: PDQSAT ! derivative according
                                                          ! to temperature of
                                                          ! saturation vapor
                                                          ! specific humidity
                                                          ! with respect to
                                                          ! water (kg/kg))
!
!*       0.2   Declarations of local variables
!
REAL, DIMENSION(SIZE(PT))                       :: ZFOES  ! saturation vapor
                                                          ! pressure
                                                          ! (Pascal)
!
!-------------------------------------------------------------------------------
!
!
!*       1.    COMPUTE SATURATION VAPOR PRESSURE
!              ---------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_THERMO:DQSATI_O_DT_1D',0,ZHOOK_HANDLE)
ZFOES(:) = PP(:) / (1.+XRD/XRV*(1./PQSAT(:)-1.))
!
!*       3.    DERIVATION ACCORDING TO TEMPERATURE
!              -----------------------------------
!
PDQSAT(:) = PQSAT(:) / (1.+(XRD/XRV-1.)*ZFOES(:)/PP(:) ) &
                   * (XBETAI/PT(:)**2 - XGAMI/PT(:))
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_THERMO:DQSATI_O_DT_1D',1,ZHOOK_HANDLE)
END FUNCTION DQSATI_O_DT_1D
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!     ##############################################################
      FUNCTION DQSATI_O_DT_3D(PT,PP,PQSAT) RESULT(PDQSAT)
!     ##############################################################
!
!!****  *QSATW * - function to compute saturation vapor humidity from
!!                 temperature (with respect to ice)
!!
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the saturation vapor
!     pressure from temperature
!
!
!!**  METHOD
!!    ------
!!       Given temperature T (PT), the saturation vapor pressure es(T)
!!    (FOES(PT)) is computed by integration of the Clapeyron equation
!!    from the triple point temperature Tt (XTT) and the saturation vapor
!!    pressure of the triple point es(Tt) (XESTT), i.e
!!
!!         es(T)= EXP( alphaw - betaw /T - gammaw Log(T) )
!!
!!     with :
!!       alphaw (XALPW) = LOG(es(Tt))+ betaw/Tt + gammaw Log(Tt)
!!       betaw (XBETAW) = Lv(Tt)/Rv + gammaw Tt
!!       gammaw (XGAMW) = (Cl -Cpv) /Rv
!!
!!      Then, the specific humidity at saturation is deduced.
!!
!!      Finally, dqsat / dT  (T) is computed.
!!
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST : comtains physical constants
!!        XALPW   : Constant for saturation vapor pressure function
!!        XBETAW  : Constant for saturation vapor pressure function
!!        XGAMW   : Constant for saturation vapor pressure function
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH
!!
!!
!!    AUTHOR
!!    ------
!!      V. Masson       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    21/09/98
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments and results
!
!
REAL,    DIMENSION(:,:,:), INTENT(IN)               :: PT     ! Temperature
                                                          ! (Kelvin)
REAL,    DIMENSION(:,:,:), INTENT(IN)               :: PP     ! Pressure
                                                          ! (Pa)
REAL,    DIMENSION(:,:,:), INTENT(IN)               :: PQSAT  ! saturation vapor
                                                          ! specific humidity
                                                          ! with respect to
                                                          ! water (kg/kg))
REAL,    DIMENSION(SIZE(PT,1),SIZE(PT,2),SIZE(PT,3))  :: PDQSAT ! derivative according
                                                          ! to temperature of
                                                          ! saturation vapor
                                                          ! specific humidity
                                                          ! with respect to
                                                          ! water (kg/kg))
!
!*       0.2   Declarations of local variables
!
REAL, DIMENSION(SIZE(PT,1),SIZE(PT,2),SIZE(PT,3))   :: ZFOES  ! saturation vapor
                                                          ! pressure
                                                          ! (Pascal)
!
!-------------------------------------------------------------------------------
!
!
!*       1.    COMPUTE SATURATION VAPOR PRESSURE
!              ---------------------------------
!
ZFOES(:,:,:) = PP(:,:,:) / (1.+XRD/XRV*(1./PQSAT(:,:,:)-1.))
!
!*       3.    DERIVATION ACCORDING TO TEMPERATURE
!              -----------------------------------
!
PDQSAT(:,:,:) = PQSAT(:,:,:) / (1.+(XRD/XRV-1.)*ZFOES(:,:,:)/PP(:,:,:) ) &
                   * (XBETAI/PT(:,:,:)**2 - XGAMI/PT(:,:,:))
!
!-------------------------------------------------------------------------------
!
END FUNCTION DQSATI_O_DT_3D
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
!     ######################################
      FUNCTION QSATI_3D(PT,PP) RESULT(PQSAT)
!     ######################################
!
!!****  *QSATI * - function to compute saturation vapor humidity from
!!                 temperature
!!
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the saturation vapor
!     pressure from temperature
!
!
!!**  METHOD
!!    ------
!!       Given temperature T (PT), the saturation vapor pressure es(T)
!!    (FOES(PT)) is computed by integration of the Clapeyron equation
!!    from the triple point temperature Tt (XTT) and the saturation vapor
!!    pressure of the triple point es(Tt) (XESTT), i.e
!!
!!         es(T)= EXP( alphaw - betaw /T - gammaw Log(T) )
!!
!!     with :
!!       alphaw (XALPI) = LOG(es(Tt))+ betaw/Tt + gammaw Log(Tt)
!!       betaw (XBETAI) = Lv(Tt)/Rv + gammaw Tt
!!       gammaw (XGAMI) = (Cl -Cpv) /Rv
!!
!!      Then, the specific humidity at saturation is deduced.
!!
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST : comtains physical constants
!!        XALPI   : Constant for saturation vapor pressure function
!!        XBETAI  : Constant for saturation vapor pressure function
!!        XGAMI   : Constant for saturation vapor pressure function
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH
!!
!!
!!    AUTHOR
!!    ------
!!      V. Masson       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    21/09/98
!!      S. Riette   april 2011 : protection in high statosphere where ZFOES > PP
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments and results
!
!
REAL, DIMENSION(:,:,:), INTENT(IN)                :: PT     ! Temperature
                                                            ! (Kelvin)
REAL, DIMENSION(:,:,:), INTENT(IN)                :: PP     ! Pressure
                                                            ! (Pa)
REAL, DIMENSION(SIZE(PT,1),SIZE(PT,2),SIZE(PT,3)) :: PQSAT  ! saturation vapor
                                                            ! specific humidity
                                                            ! with respect to
                                                            ! water (kg/kg)
!
!*       0.2   Declarations of local variables
!
REAL, DIMENSION(SIZE(PT,1),SIZE(PT,2),SIZE(PT,3)) :: ZFOES  ! saturation vapor
                                                            ! pressure
                                                            ! (Pascal)
!
!-------------------------------------------------------------------------------
!
!*       1.    COMPUTE SATURATION VAPOR PRESSURE
!              ---------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_THERMO:QSATI_3D',0,ZHOOK_HANDLE)
ZFOES(:,:,:) = MIN(EXP( XALPI - XBETAI/PT(:,:,:) - XGAMI*LOG(PT(:,:,:))  ), 0.99*PP(:,:,:))
!
!*       2.    COMPUTE SATURATION HUMIDITY
!              ---------------------------
!
PQSAT(:,:,:) = XRD/XRV*ZFOES(:,:,:)/PP(:,:,:)   &
                   / (1.+(XRD/XRV-1.)*ZFOES(:,:,:)/PP(:,:,:))
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_THERMO:QSATI_3D',1,ZHOOK_HANDLE)
END FUNCTION QSATI_3D
!
!-------------------------------------------------------------------------------
!
!     ######################################
      FUNCTION QSATI_2D(PT,PP) RESULT(PQSAT)
!     ######################################
!
!!****  *QSATI * - function to compute saturation vapor humidity from
!!                 temperature
!!
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the saturation vapor
!     pressure from temperature
!
!
!!**  METHOD
!!    ------
!!       Given temperature T (PT), the saturation vapor pressure es(T)
!!    (FOES(PT)) is computed by integration of the Clapeyron equation
!!    from the triple point temperature Tt (XTT) and the saturation vapor
!!    pressure of the triple point es(Tt) (XESTT), i.e
!!
!!         es(T)= EXP( alphaw - betaw /T - gammaw Log(T) )
!!
!!     with :
!!       alphaw (XALPI) = LOG(es(Tt))+ betaw/Tt + gammaw Log(Tt)
!!       betaw (XBETAI) = Lv(Tt)/Rv + gammaw Tt
!!       gammaw (XGAMI) = (Cl -Cpv) /Rv
!!
!!      Then, the specific humidity at saturation is deduced.
!!
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST : comtains physical constants
!!        XALPI   : Constant for saturation vapor pressure function
!!        XBETAI  : Constant for saturation vapor pressure function
!!        XGAMI   : Constant for saturation vapor pressure function
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH
!!
!!
!!    AUTHOR
!!    ------
!!      V. Masson       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    21/09/98
!!      S. Riette   april 2011 : protection in high statosphere where ZFOES > PP
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments and results
!
!
REAL, DIMENSION(:,:), INTENT(IN)                :: PT     ! Temperature
                                                          ! (Kelvin)
REAL, DIMENSION(:,:), INTENT(IN)                :: PP     ! Pressure
                                                          ! (Pa)
REAL, DIMENSION(SIZE(PT,1),SIZE(PT,2))          :: PQSAT  ! saturation vapor
                                                          ! specific humidity
                                                          ! with respect to
                                                          ! water (kg/kg)
!
!*       0.2   Declarations of local variables
!
REAL, DIMENSION(SIZE(PT,1),SIZE(PT,2))          :: ZFOES  ! saturation vapor
                                                          ! pressure
                                                          ! (Pascal)
!
!-------------------------------------------------------------------------------
!
!*       1.    COMPUTE SATURATION VAPOR PRESSURE
!              ---------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_THERMO:QSATI_2D',0,ZHOOK_HANDLE)
ZFOES(:,:) = MIN(EXP( XALPI - XBETAI/PT(:,:) - XGAMI*LOG(PT(:,:))  ), 0.99*PP(:,:))
!
!*       2.    COMPUTE SATURATION HUMIDITY
!              ---------------------------
!
PQSAT(:,:) = XRD/XRV*ZFOES(:,:)/PP(:,:)   &
                   / (1.+(XRD/XRV-1.)*ZFOES(:,:)/PP(:,:))
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_THERMO:QSATI_2D',1,ZHOOK_HANDLE)
END FUNCTION QSATI_2D
!
!-------------------------------------------------------------------------------
!
!     #################################################
      FUNCTION QSATI_2D_MASK(OMASK,PT,PP) RESULT(PQSAT)
!     #################################################
!
!!****  *QSATI * - function to compute saturation vapor humidity from
!!                 temperature
!!
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the saturation vapor
!     pressure from temperature
!
!
!!**  METHOD
!!    ------
!!       Given temperature T (PT), the saturation vapor pressure es(T)
!!    (FOES(PT)) is computed by integration of the Clapeyron equation
!!    from the triple point temperature Tt (XTT) and the saturation vapor
!!    pressure of the triple point es(Tt) (XESTT), i.e
!!
!!         es(T)= EXP( alphaw - betaw /T - gammaw Log(T) )
!!
!!     with :
!!       alphaw (XALPI) = LOG(es(Tt))+ betaw/Tt + gammaw Log(Tt)
!!       betaw (XBETAI) = Lv(Tt)/Rv + gammaw Tt
!!       gammaw (XGAMI) = (Cl -Cpv) /Rv
!!
!!      Then, the specific humidity at saturation is deduced.
!!
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST : comtains physical constants
!!        XALPI   : Constant for saturation vapor pressure function
!!        XBETAI  : Constant for saturation vapor pressure function
!!        XGAMI   : Constant for saturation vapor pressure function
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH
!!
!!
!!    AUTHOR
!!    ------
!!      V. Masson       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    21/09/98
!!      S. Riette   april 2011 : protection in high statosphere where ZFOES > PP
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments and results
!
!
LOGICAL, DIMENSION(:,:), INTENT(IN)             :: OMASK  ! Localization mask
REAL, DIMENSION(:,:), INTENT(IN)                :: PT     ! Temperature
                                                          ! (Kelvin)
REAL, DIMENSION(:,:), INTENT(IN)                :: PP     ! Pressure
                                                          ! (Pa)
REAL, DIMENSION(SIZE(PT,1),SIZE(PT,2))          :: PQSAT  ! saturation vapor
                                                          ! specific humidity
                                                          ! with respect to
                                                          ! water (kg/kg)
!
!*       0.2   Declarations of local variables
!
REAL, DIMENSION(SIZE(PT,1),SIZE(PT,2))          :: ZFOES  ! saturation vapor
                                                          ! pressure
                                                          ! (Pascal)
!
!-------------------------------------------------------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_THERMO:QSATI_2D_MASK',0,ZHOOK_HANDLE)
WHERE (OMASK(:,:))
!
!*       1.    COMPUTE SATURATION VAPOR PRESSURE
!              ---------------------------------
!
  ZFOES(:,:) = MIN(EXP( XALPI - XBETAI/PT(:,:) - XGAMI*LOG(PT(:,:))  ), 0.99*PP(:,:))
!
!*       2.    COMPUTE SATURATION HUMIDITY
!              ---------------------------
!
  PQSAT(:,:) = XRD/XRV*ZFOES(:,:)/PP(:,:)   &
                     / (1.+(XRD/XRV-1.)*ZFOES(:,:)/PP(:,:))
ELSEWHERE
!
!*       3.    BOGUS VALUE
!              -----------
!
  PQSAT(:,:) = 0.
END WHERE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_THERMO:QSATI_2D_MASK',1,ZHOOK_HANDLE)
END FUNCTION QSATI_2D_MASK
!
!-------------------------------------------------------------------------------
!
!     ######################################
      FUNCTION QSATI_1D(PT,PP) RESULT(PQSAT)
!     ######################################
!
!!****  *QSATI * - function to compute saturation vapor humidity from
!!                 temperature
!!
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the saturation vapor
!     pressure from temperature
!
!
!!**  METHOD
!!    ------
!!       Given temperature T (PT), the saturation vapor pressure es(T)
!!    (FOES(PT)) is computed by integration of the Clapeyron equation
!!    from the triple point temperature Tt (XTT) and the saturation vapor
!!    pressure of the triple point es(Tt) (XESTT), i.e
!!
!!         es(T)= EXP( alphaw - betaw /T - gammaw Log(T) )
!!
!!     with :
!!       alphaw (XALPI) = LOG(es(Tt))+ betaw/Tt + gammaw Log(Tt)
!!       betaw (XBETAI) = Lv(Tt)/Rv + gammaw Tt
!!       gammaw (XGAMI) = (Cl -Cpv) /Rv
!!
!!      Then, the specific humidity at saturation is deduced.
!!
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST : comtains physical constants
!!        XALPI   : Constant for saturation vapor pressure function
!!        XBETAI  : Constant for saturation vapor pressure function
!!        XGAMI   : Constant for saturation vapor pressure function
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH
!!
!!
!!    AUTHOR
!!    ------
!!      V. Masson       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    21/09/98
!!      S. Riette   april 2011 : protection in high statosphere where ZFOES > PP
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments and results
!
!
REAL, DIMENSION(:), INTENT(IN)                :: PT     ! Temperature
                                                        ! (Kelvin)
REAL, DIMENSION(:), INTENT(IN)                :: PP     ! Pressure
                                                        ! (Pa)
REAL, DIMENSION(SIZE(PT,1))                   :: PQSAT  ! saturation vapor
                                                        ! specific humidity
                                                        ! with respect to
                                                        ! water (kg/kg)
!
!*       0.2   Declarations of local variables
!
REAL, DIMENSION(SIZE(PT,1))                   :: ZFOES  ! saturation vapor
                                                        ! pressure
                                                        ! (Pascal)
!
!-------------------------------------------------------------------------------
!
!*       1.    COMPUTE SATURATION VAPOR PRESSURE
!              ---------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_THERMO:QSATI_1D',0,ZHOOK_HANDLE)
ZFOES(:) = MIN(EXP( XALPI - XBETAI/PT(:) - XGAMI*LOG(PT(:))  ), 0.99*PP(:))
!
!*       2.    COMPUTE SATURATION HUMIDITY
!              ---------------------------
!
PQSAT(:) = XRD/XRV*ZFOES(:)/PP(:)   &
                   / (1.+(XRD/XRV-1.)*ZFOES(:)/PP(:))
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_THERMO:QSATI_1D',1,ZHOOK_HANDLE)
END FUNCTION QSATI_1D
!
!-------------------------------------------------------------------------------
!
!     ######################################
      FUNCTION QSATI_0D(PT,PP) RESULT(PQSAT)
!     ######################################
!
!!****  *QSATI * - function to compute saturation vapor humidity from
!!                 temperature
!!
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the saturation vapor
!     pressure from temperature
!
!
!!**  METHOD
!!    ------
!!       Given temperature T (PT), the saturation vapor pressure es(T)
!!    (FOES(PT)) is computed by integration of the Clapeyron equation
!!    from the triple point temperature Tt (XTT) and the saturation vapor
!!    pressure of the triple point es(Tt) (XESTT), i.e
!!
!!         es(T)= EXP( alphaw - betaw /T - gammaw Log(T) )
!!
!!     with :
!!       alphaw (XALPI) = LOG(es(Tt))+ betaw/Tt + gammaw Log(Tt)
!!       betaw (XBETAI) = Lv(Tt)/Rv + gammaw Tt
!!       gammaw (XGAMI) = (Cl -Cpv) /Rv
!!
!!      Then, the specific humidity at saturation is deduced.
!!
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST : comtains physical constants
!!        XALPI   : Constant for saturation vapor pressure function
!!        XBETAI  : Constant for saturation vapor pressure function
!!        XGAMI   : Constant for saturation vapor pressure function
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH
!!
!!
!!    AUTHOR
!!    ------
!!      V. Masson       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    21/09/98
!!      S. Riette   april 2011 : protection in high statosphere where ZFOES > PP
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments and results
!
!
REAL, INTENT(IN)                :: PT     ! Temperature
                                          ! (Kelvin)
REAL, INTENT(IN)                :: PP     ! Pressure
                                          ! (Pa)
REAL                            :: PQSAT  ! saturation vapor
                                          ! specific humidity
                                          ! with respect to
                                          ! water (kg/kg)
!
!*       0.2   Declarations of local variables
!
REAL                            :: ZFOES  ! saturation vapor
                                          ! pressure
                                          ! (Pascal)
!
!-------------------------------------------------------------------------------
!
!*       1.    COMPUTE SATURATION VAPOR PRESSURE
!              ---------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_THERMO:QSATI_0D',0,ZHOOK_HANDLE)
ZFOES = MIN(EXP( XALPI - XBETAI/PT - XGAMI*LOG(PT)  ), 0.99*PP)
!
!*       2.    COMPUTE SATURATION HUMIDITY
!              ---------------------------
!
PQSAT = XRD/XRV*ZFOES/PP / (1.+(XRD/XRV-1.)*ZFOES/PP)
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_THERMO:QSATI_0D',1,ZHOOK_HANDLE)
END FUNCTION QSATI_0D
!
!-------------------------------------------------------------------------------
END MODULE MODE_THERMO
