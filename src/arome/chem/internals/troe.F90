!     ######spl
       FUNCTION TROE(PCOEF,PKO, PNEXP, PKINF, PMEXP, PM, PT,KVECNPT)
       USE PARKIND1, ONLY : JPRB
       USE YOMHOOK , ONLY : LHOOK, DR_HOOK
 !!    #############################################################
 !!
 !!*** *TROE*
 !!
 !!    PURPOSE
 !!    -------
 !     this function implements the TROE reaction rate for RACM
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
 REAL, DIMENSION(KVECNPT)             :: TROE 
 REAL,                     INTENT(IN) :: PCOEF,PKO, PNEXP, PKINF, PMEXP
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
 IF (LHOOK) CALL DR_HOOK('TROE',0,ZHOOK_HANDLE)
 ZKOTM(:)  = PM(:) * PKO * ( (PT(:)/300.)**(-PNEXP) )
 ZKINFT(:) = PKINF * ( (PT(:)/300.)**(-PMEXP) )
 ZFACT(:)  = 0.6**(1./(1.+ALOG10(ZKOTM(:)/ZKINFT(:))**2 ))
 TROE(:)   = PCOEF*(ZKOTM(:)/(1.+ZKOTM(:)/ZKINFT(:)))*ZFACT(:)
 !
 IF (LHOOK) CALL DR_HOOK('TROE',1,ZHOOK_HANDLE)
 END FUNCTION TROE
