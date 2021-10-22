!     ######spl
      SUBROUTINE CH_GET_RATES(PRATE,KMI,KVECNPT,KREAC)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!!    ############################################
!! This code has been created automatically by preprocessor m10,
!! version: 9.7, copyright 1995-1999 by Meteo France/Universite Paul Sabatier.
!! Please report all bugs to K. Suhre (Lab. d'Aerologie UPS/CNRS).
!!
!!*** *MODD_CH_GETRATES*
!!
!!    PURPOSE
!!    -------
!       retrieve reaction rates from TPK in an array
!!
!!**  METHOD
!!    ------
!!      simple
!!
!!    REFERENCE
!!    ---------
!!    none
!!
!!    AUTHOR
!!    ------
!!    Karsten Suhre (LA)
!!
!!    MODIFICATIONS
!!    -------------
!!    Original 26/07/96
!!    Modified 05/05/98: Vectorization (Vincent Crassier & KS)
!!    Modified 31/10/03: New interface for better MesoNH compilation (D. Gazen)
!!
!!----------------------------------------------------------------------
!!
!!    EXTERNAL
!!    --------
!!    none
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
USE MODD_CH_M9_SCHEME
!!
!!    EXPLICIT ARGUMENTS
!!    ------------------
IMPLICIT NONE
INTEGER, INTENT(IN)                         :: KVECNPT
INTEGER, INTENT(IN)                         :: KREAC
REAL, INTENT(OUT), DIMENSION(KVECNPT,KREAC) :: PRATE
INTEGER, INTENT(IN)                         :: KMI
!!
!!    LOCAL VARIABLES
!!    ---------------
TYPE(CCSTYPE), POINTER                      :: TPK
!!----------------------------------------------------------------------
!!
!!    EXECUTABLE STATEMENTS
!!    ---------------------
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('CH_GET_RATES',0,ZHOOK_HANDLE)
TPK=>TACCS(KMI)
!!
PRATE(:,1) = TPK%K001(:)
PRATE(:,2) = TPK%K002(:)
PRATE(:,3) = TPK%K003(:)
PRATE(:,4) = TPK%K004(:)
PRATE(:,5) = TPK%K005(:)
PRATE(:,6) = TPK%K006(:)
PRATE(:,7) = TPK%K007(:)
PRATE(:,8) = TPK%K008(:)
PRATE(:,9) = TPK%K009(:)
PRATE(:,10) = TPK%K010(:)
PRATE(:,11) = TPK%K011(:)
PRATE(:,12) = TPK%K012(:)
PRATE(:,13) = TPK%K013(:)
PRATE(:,14) = TPK%K014(:)
PRATE(:,15) = TPK%K015(:)
PRATE(:,16) = TPK%K016(:)
PRATE(:,17) = TPK%K017(:)
PRATE(:,18) = TPK%K018(:)
PRATE(:,19) = TPK%K019(:)
PRATE(:,20) = TPK%K020(:)
PRATE(:,21) = TPK%K021(:)
PRATE(:,22) = TPK%K022(:)
PRATE(:,23) = TPK%K023(:)
PRATE(:,24) = TPK%K024(:)
PRATE(:,25) = TPK%K025(:)
PRATE(:,26) = TPK%K026(:)
PRATE(:,27) = TPK%K027(:)
PRATE(:,28) = TPK%K028(:)
PRATE(:,29) = TPK%K029(:)
PRATE(:,30) = TPK%K030(:)
PRATE(:,31) = TPK%K031(:)
PRATE(:,32) = TPK%K032(:)
PRATE(:,33) = TPK%K033(:)
PRATE(:,34) = TPK%K034(:)
PRATE(:,35) = TPK%K035(:)
PRATE(:,36) = TPK%K036(:)
PRATE(:,37) = TPK%K037(:)
PRATE(:,38) = TPK%K038(:)
PRATE(:,39) = TPK%K039(:)
PRATE(:,40) = TPK%K040(:)
PRATE(:,41) = TPK%K041(:)
PRATE(:,42) = TPK%K042(:)
PRATE(:,43) = TPK%K043(:)
PRATE(:,44) = TPK%K044(:)
PRATE(:,45) = TPK%K045(:)
PRATE(:,46) = TPK%K046(:)
PRATE(:,47) = TPK%K047(:)
PRATE(:,48) = TPK%K048(:)
PRATE(:,49) = TPK%K049(:)
PRATE(:,50) = TPK%K050(:)
PRATE(:,51) = TPK%K051(:)
PRATE(:,52) = TPK%K052(:)
PRATE(:,53) = TPK%K053(:)
PRATE(:,54) = TPK%K054(:)
PRATE(:,55) = TPK%K055(:)
PRATE(:,56) = TPK%K056(:)
PRATE(:,57) = TPK%K057(:)
PRATE(:,58) = TPK%K058(:)
PRATE(:,59) = TPK%K059(:)
PRATE(:,60) = TPK%K060(:)
PRATE(:,61) = TPK%K061(:)
PRATE(:,62) = TPK%K062(:)
PRATE(:,63) = TPK%K063(:)
PRATE(:,64) = TPK%K064(:)
PRATE(:,65) = TPK%K065(:)
PRATE(:,66) = TPK%K066(:)
PRATE(:,67) = TPK%K067(:)
PRATE(:,68) = TPK%K068(:)
PRATE(:,69) = TPK%K069(:)
PRATE(:,70) = TPK%K070(:)
PRATE(:,71) = TPK%K071(:)
PRATE(:,72) = TPK%K072(:)
PRATE(:,73) = TPK%K073(:)
PRATE(:,74) = TPK%K074(:)
PRATE(:,75) = TPK%K075(:)
PRATE(:,76) = TPK%K076(:)
PRATE(:,77) = TPK%K077(:)
PRATE(:,78) = TPK%K078(:)
PRATE(:,79) = TPK%K079(:)
PRATE(:,80) = TPK%K080(:)
PRATE(:,81) = TPK%K081(:)
PRATE(:,82) = TPK%K082(:)
PRATE(:,83) = TPK%K083(:)
PRATE(:,84) = TPK%K084(:)
PRATE(:,85) = TPK%K085(:)
PRATE(:,86) = TPK%K086(:)
PRATE(:,87) = TPK%K087(:)
PRATE(:,88) = TPK%K088(:)
PRATE(:,89) = TPK%K089(:)
PRATE(:,90) = TPK%K090(:)
PRATE(:,91) = TPK%K091(:)
PRATE(:,92) = TPK%K092(:)
PRATE(:,93) = TPK%K093(:)
PRATE(:,94) = TPK%K094(:)
PRATE(:,95) = TPK%K095(:)
PRATE(:,96) = TPK%K096(:)
PRATE(:,97) = TPK%K097(:)
PRATE(:,98) = TPK%K098(:)
PRATE(:,99) = TPK%K099(:)
PRATE(:,100) = TPK%K100(:)
PRATE(:,101) = TPK%K101(:)
PRATE(:,102) = TPK%K102(:)
PRATE(:,103) = TPK%K103(:)
PRATE(:,104) = TPK%K104(:)
PRATE(:,105) = TPK%K105(:)
PRATE(:,106) = TPK%K106(:)
PRATE(:,107) = TPK%K107(:)
PRATE(:,108) = TPK%K108(:)
PRATE(:,109) = TPK%K109(:)
PRATE(:,110) = TPK%K110(:)
PRATE(:,111) = TPK%K111(:)
PRATE(:,112) = TPK%K112(:)
PRATE(:,113) = TPK%K113(:)
PRATE(:,114) = TPK%K114(:)
PRATE(:,115) = TPK%K115(:)
PRATE(:,116) = TPK%K116(:)
PRATE(:,117) = TPK%K117(:)
PRATE(:,118) = TPK%K118(:)
PRATE(:,119) = TPK%K119(:)
PRATE(:,120) = TPK%K120(:)
PRATE(:,121) = TPK%K121(:)
PRATE(:,122) = TPK%K122(:)
PRATE(:,123) = TPK%K123(:)
PRATE(:,124) = TPK%K124(:)
PRATE(:,125) = TPK%K125(:)
PRATE(:,126) = TPK%K126(:)
PRATE(:,127) = TPK%K127(:)
PRATE(:,128) = TPK%K128(:)
IF (LHOOK) CALL DR_HOOK('CH_GET_RATES',1,ZHOOK_HANDLE)
RETURN
END SUBROUTINE CH_GET_RATES
