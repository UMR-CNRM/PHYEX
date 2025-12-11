MODULE MODI_SHALLOW_CONVECTION_PART1

INTERFACE
SUBROUTINE SHALLOW_CONVECTION_PART1&
                             (CVPEXT, CVP_SHAL, CST, D, NSV, CONVPAR, KBDIA, KTDIA, &
                              KICE, OSETTADJ, PTADJS, PPABST, PZZ, &
                              PTKECLS, PTT, PRVT, PRCT, PRIT, PWT, &
                              PTTEN, PRVTEN, PRCTEN, PRITEN,       &
                              KCLTOP, KCLBAS, PUMF, OCH1CONV, KCH1,&
                              PCH1, PCH1TEN, PTHT, PSTHV, PSTHES,  &
                              KSDPL, KSPBL, KSLCL, PSTHLCL, PSTLCL,&
                              PSRVLCL, PSWLCL, PSZLCL, PSTHVELCL, OTRIG1)
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
REAL,                         INTENT(IN) :: PTT(D%NIT,D%NKT)      
REAL,                         INTENT(IN) :: PRVT(D%NIT,D%NKT)     
REAL,                         INTENT(IN) :: PRCT(D%NIT,D%NKT)     
REAL,                         INTENT(IN) :: PRIT(D%NIT,D%NKT)     
REAL,                         INTENT(IN) :: PWT(D%NIT,D%NKT)      
REAL,                         INTENT(IN) :: PPABST(D%NIT,D%NKT)   
REAL,                         INTENT(IN) :: PZZ(D%NIT,D%NKT)      
REAL,                         INTENT(IN) :: PTKECLS(D%NIT)  
REAL,                         INTENT(INOUT):: PTTEN(D%NIT,D%NKT)  
REAL,                         INTENT(INOUT):: PRVTEN(D%NIT,D%NKT) 
REAL,                         INTENT(INOUT):: PRCTEN(D%NIT,D%NKT) 
REAL,                         INTENT(INOUT):: PRITEN(D%NIT,D%NKT) 
INTEGER,                      INTENT(INOUT):: KCLTOP(D%NIT) 
INTEGER,                      INTENT(INOUT):: KCLBAS(D%NIT) 
REAL,                         INTENT(INOUT):: PUMF(D%NIT,D%NKT)   
LOGICAL,                      INTENT(IN) :: OCH1CONV 
INTEGER,                      INTENT(IN) :: KCH1     
REAL,                         INTENT(IN) :: PCH1(D%NIT,D%NKT,KCH1)
REAL,                         INTENT(INOUT):: PCH1TEN(D%NIT,D%NKT,KCH1)
REAL,                         INTENT(OUT)  :: PTHT(D%NIT,D%NKT) 
REAL,                         INTENT(OUT)  :: PSTHV(D%NIT,D%NKT) 
REAL,                         INTENT(OUT)  :: PSTHES(D%NIT,D%NKT)  
INTEGER,                      INTENT(OUT)  :: KSDPL(D%NIT)   
INTEGER,                      INTENT(OUT)  :: KSPBL(D%NIT)   
INTEGER,                      INTENT(OUT)  :: KSLCL(D%NIT)   
REAL,                         INTENT(OUT)  :: PSTHLCL(D%NIT) 
REAL,                         INTENT(OUT)  :: PSTLCL(D%NIT)  
REAL,                         INTENT(OUT)  :: PSRVLCL(D%NIT) 
REAL,                         INTENT(OUT)  :: PSWLCL(D%NIT)  
REAL,                         INTENT(OUT)  :: PSZLCL(D%NIT)  
REAL,                         INTENT(OUT)  :: PSTHVELCL(D%NIT)
LOGICAL,                      INTENT(OUT)  :: OTRIG1(D%NIT)  

END SUBROUTINE
END INTERFACE

END MODULE
