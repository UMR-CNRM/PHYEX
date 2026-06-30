MODULE MODI_SHALLOW_CONVECTION_PART1

IMPLICIT NONE

INTERFACE
SUBROUTINE SHALLOW_CONVECTION_PART1(&
                              CVPEXT, CVP_SHAL, CST, D, NSV, CONVPAR, KBDIA, KTDIA, &
                              KICE, OSETTADJ, PTADJS, PPABST, PZZ, &
                              PTKECLS, PTT, PRVT, PRCT, PRIT, PWT, &
                              PTTEN, PRVTEN, PRCTEN, PRITEN,       &
                              KCLTOP, KCLBAS, PUMF, OCH1CONV, KCH1,&
                              PCH1, PCH1TEN, PTHT, PSTHV, PSTHES,  &
                              KSDPL, KSPBL, KSLCL, PSTHLCL, PSTLCL,&
                              PSRVLCL, PSWLCL, PSZLCL, PSTHVELCL, OTRIG1)

!$ACDC singlecolumn

USE MODD_CST, ONLY : CST_T
USE MODD_CONVPAREXT, ONLY : CONVPAREXT
USE MODD_CONVPAR, ONLY: CONVPAR_T
USE MODD_CONVPAR_SHAL, ONLY: CONVPAR_SHAL
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_T
USE MODD_NSV, ONLY: NSV_T
IMPLICIT NONE
TYPE(CONVPAREXT),             INTENT(IN) :: CVPEXT
TYPE(CONVPAR_SHAL),           INTENT(IN) :: CVP_SHAL
TYPE(CST_T),                  INTENT(IN) :: CST
TYPE(DIMPHYEX_T),             INTENT(IN) :: D
TYPE(NSV_T),                  INTENT(IN) :: NSV
TYPE(CONVPAR_T),              INTENT(IN) :: CONVPAR
INTEGER,                      INTENT(IN) :: KBDIA    
INTEGER,                      INTENT(IN) :: KTDIA    
INTEGER,                      INTENT(IN) :: KICE     
LOGICAL,                      INTENT(IN) :: OSETTADJ 
REAL,                         INTENT(IN) :: PTADJS   
REAL,                         INTENT(IN) :: PTT(D%NIJT,D%NKT)      
REAL,                         INTENT(IN) :: PRVT(D%NIJT,D%NKT)     
REAL,                         INTENT(IN) :: PRCT(D%NIJT,D%NKT)     
REAL,                         INTENT(IN) :: PRIT(D%NIJT,D%NKT)     
REAL,                         INTENT(IN) :: PWT(D%NIJT,D%NKT)      
REAL,                         INTENT(IN) :: PPABST(D%NIJT,D%NKT)   
REAL,                         INTENT(IN) :: PZZ(D%NIJT,D%NKT)      
REAL,                         INTENT(IN) :: PTKECLS(D%NIJT)  
REAL,                         INTENT(INOUT):: PTTEN(D%NIJT,D%NKT)  
REAL,                         INTENT(INOUT):: PRVTEN(D%NIJT,D%NKT) 
REAL,                         INTENT(INOUT):: PRCTEN(D%NIJT,D%NKT) 
REAL,                         INTENT(INOUT):: PRITEN(D%NIJT,D%NKT) 
INTEGER,                      INTENT(INOUT):: KCLTOP(D%NIJT) 
INTEGER,                      INTENT(INOUT):: KCLBAS(D%NIJT) 
REAL,                         INTENT(INOUT):: PUMF(D%NIJT,D%NKT)   
LOGICAL,                      INTENT(IN) :: OCH1CONV 
INTEGER,                      INTENT(IN) :: KCH1     
REAL,                         INTENT(IN) :: PCH1(D%NIJT,D%NKT,KCH1)
REAL,                         INTENT(INOUT):: PCH1TEN(D%NIJT,D%NKT,KCH1)
REAL,                         INTENT(OUT)  :: PTHT(D%NIJT,D%NKT) 
REAL,                         INTENT(OUT)  :: PSTHV(D%NIJT,D%NKT) 
REAL,                         INTENT(OUT)  :: PSTHES(D%NIJT,D%NKT)  
INTEGER,                      INTENT(OUT)  :: KSDPL(D%NIJT)   
INTEGER,                      INTENT(OUT)  :: KSPBL(D%NIJT)   
INTEGER,                      INTENT(OUT)  :: KSLCL(D%NIJT)   
REAL,                         INTENT(OUT)  :: PSTHLCL(D%NIJT) 
REAL,                         INTENT(OUT)  :: PSTLCL(D%NIJT)  
REAL,                         INTENT(OUT)  :: PSRVLCL(D%NIJT) 
REAL,                         INTENT(OUT)  :: PSWLCL(D%NIJT)  
REAL,                         INTENT(OUT)  :: PSZLCL(D%NIJT)  
REAL,                         INTENT(OUT)  :: PSTHVELCL(D%NIJT)
LOGICAL,                      INTENT(OUT)  :: OTRIG1(D%NIJT)  

END SUBROUTINE
END INTERFACE

END MODULE
