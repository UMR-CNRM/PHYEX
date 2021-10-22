!     ######spl
       FUNCTION TROE_EQUIL(PKO, PNEXP, PKINF, PMEXP, PAFACT, PB, PM, PT,KVECNPT)
       USE PARKIND1, ONLY : JPRB
       USE YOMHOOK , ONLY : LHOOK, DR_HOOK
 !!    #########################################################################
 !!
 !!*** *TROE_EQUIL*
 !!
 !!    PURPOSE
 !!    -------
 !     this function implements the TROE_EQUIL reaction rate for RACM
 !!
 !!    REFERENCE
 !!    ---------
 !!    Stockwell et al., JGR, 1997
 !!
 !!    AUTHOR
 !!    ------
 !!    Karsten Suhre (LA)
 !!    
 !!    MODIFICATIONS
 !!    -------------
 !!    Original 27/01/98
 !!
 !!------------------------------------------------------------------------------
 !!
 !!    EXTERNAL
 !!    --------
 !!    none
 !!
 !!    IMPLICIT ARGUMENTS
 !!    ------------------
 !!    none
 !!
 !!    EXPLICIT ARGUMENTS
 !!    ------------------
 IMPLICIT NONE
 INTEGER,                  INTENT(IN) :: KVECNPT ! no. of points in vector mask
 REAL,DIMENSION(KVECNPT)              :: TROE_EQUIL 
 REAL,                     INTENT(IN) :: PKO, PNEXP, PKINF, PMEXP, PAFACT, PB
 REAL, DIMENSION(KVECNPT), INTENT(IN) :: PM, PT
 !!
 !!    LOCAL VARIABLES
 !!    ---------------
 REAL, DIMENSION(KVECNPT) :: ZKOTM, ZKINFT, ZFACT
 !!
 !------------------------------------------------------------------------------
 !!
 !!    EXECUTABLE STATEMENTS
 !!    ---------------------
 !
 !*        1. THE EXPRESSION
 !         -----------------
 !
 REAL(KIND=JPRB) :: ZHOOK_HANDLE
 IF (LHOOK) CALL DR_HOOK('TROE_EQUIL',0,ZHOOK_HANDLE)
 ZKOTM(:)      = PM(:) * PKO * ( (PT(:)/300.)**(-PNEXP) )
 ZKINFT(:)     = PKINF * ( (PT(:)/300.)**(-PMEXP) )
 ZFACT(:)      = 0.6**(1./(1.+ALOG10(ZKOTM(:)/ZKINFT(:))**2 ))
 TROE_EQUIL(:) = PAFACT*exp(-PB/PT)*(ZKOTM/(1.+ZKOTM(:)/ZKINFT(:)))*ZFACT(:)
 !
 IF (LHOOK) CALL DR_HOOK('TROE_EQUIL',1,ZHOOK_HANDLE)
 END FUNCTION TROE_EQUIL
