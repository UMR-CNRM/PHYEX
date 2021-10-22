!     ######spl
     SUBROUTINE CH_AER_MOD_INIT
     USE PARKIND1, ONLY : JPRB
     USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!!   ####################################
!!
!!    PURPOSE
!!    -------
!!     initialize the aerosol module (to be called only once)
!!
!!    METHOD
!!    ------
!!
!!      allocate all arrays and initialize the basic variables (i.e. densities
!!    and molar weights)
!!
!!    REFERENCE
!!    ---------
!!    none
!!
!!    AUTHOR
!!    ------
!!    Vincent Crassier (LA)
!!
!!    MODIFICATIONS
!!    -------------
!!    20/03/03   P . Tulet (CNRM/GMEI)   add  initialization tabulation
!!
!!    EXTERNAL
!!    --------
!!    none
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
USE MODD_CH_AEROSOL
!!

!*       0.     DECLARATIONS
!               ------------
!
IMPLICIT NONE
!
!* arguments
!
!
!* local variables 
INTEGER, PARAMETER :: nc=22, nh=16, nt=11 ! inorganic interpolation
INTEGER             :: JI ! loop counter
INTEGER :: i,j,k,l,m, JJ !loop
INTEGER :: NRESP     ! return code in FM routines
INTEGER :: NLUOUT    ! unit for output listing count
REAL :: ZPI
!
!---------------------------------------------------------------------------
!
!
!
!        1.1    initialisation of tables

! Initialize the mineral tablution 
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('CH_AER_MOD_INIT',0,ZHOOK_HANDLE)
NLUOUT = 77 
IF (CMINERAL == 'NARES') THEN
!       .. the file ares.w contains the weights of the model
        OPEN(NLUOUT,FILE="ares1A.w",STATUS="OLD") 
        READ(NLUOUT,*) I1IA,J1JA,K1KA 
        DO I=1,I1IA 
          READ(NLUOUT,*) X1MAXA(1,I),X1MINA(1,I),X1MODA(1,I) 
        ENDDO
        DO I=1,K1KA 
          READ(NLUOUT,*) X1MAXA(2,I),X1MINA(2,I),X1MODA(2,I) 
        ENDDO
        DO I=1,I1IA+1 
          READ(NLUOUT,*) (W1IJA(I,J),J=1,J1JA) 
        ENDDO
        DO J=1,J1JA+1 
          READ(NLUOUT,*) (W1JKA(J,K),K=1,K1KA) 
        ENDDO
        CLOSE(NLUOUT)

        OPEN(NLUOUT,FILE="ares1C.w",STATUS="OLD") 
        READ(NLUOUT,*) I1IC,J1JC,K1KC 
        DO I=1,I1IC 
          READ(NLUOUT,*) X1MAXC(1,I),X1MINC(1,I),X1MODC(1,I) 
        ENDDO
        DO I=1,K1KC 
          READ(NLUOUT,*) X1MAXC(2,I),X1MINC(2,I),X1MODC(2,I) 
        ENDDO
        DO I=1,I1IC+1 
          READ(NLUOUT,*) (W1IJC(I,J),J=1,J1JC) 
        ENDDO
        DO J=1,J1JC+1 
          READ(NLUOUT,*) (W1JKC(J,K),K=1,K1KC) 
        ENDDO
        CLOSE(NLUOUT)

        OPEN(NLUOUT,FILE="ares2A.w",STATUS="OLD") 
        READ(NLUOUT,*) I2IA,J2JA,K2KA 
        DO I=1,I2IA 
          READ(NLUOUT,*) X2MAXA(1,I),X2MINA(1,I),X2MODA(1,I) 
        ENDDO
        DO I=1,K2KA 
          READ(NLUOUT,*) X2MAXA(2,I),X2MINA(2,I),X2MODA(2,I) 
        ENDDO
        DO I=1,I2IA+1 
          READ(NLUOUT,*) (W2IJA(I,J),J=1,J2JA) 
        ENDDO
        DO J=1,J2JA+1 
          READ(NLUOUT,*) (W2JKA(J,K),K=1,K2KA) 
        ENDDO
        CLOSE(NLUOUT)

        OPEN(NLUOUT,FILE="ares2B.w",STATUS="OLD") 
        READ(NLUOUT,*) I2IB,J2JB,K2KB 
        DO I=1,I2IB 
          READ(NLUOUT,*) X2MAXB(1,I),X2MINB(1,I),X2MODB(1,I) 
        ENDDO
        DO I=1,K2KB 
          READ(NLUOUT,*) X2MAXB(2,I),X2MINB(2,I),X2MODB(2,I) 
        ENDDO
        DO I=1,I2IB+1 
          READ(NLUOUT,*) (W2IJB(I,J),J=1,J2JB) 
        ENDDO
        DO J=1,J2JB+1 
          READ(NLUOUT,*) (W2JKB(J,K),K=1,K2KB) 
        ENDDO
        CLOSE(NLUOUT)

        OPEN(NLUOUT,FILE="ares2C.w",STATUS="OLD") 
        READ(NLUOUT,*) I2IC,J2JC,K2KC 
        DO I=1,I2IC 
          READ(NLUOUT,*) X2MAXC(1,I),X2MINC(1,I),X2MODC(1,I) 
        ENDDO
        DO I=1,K2KC 
          READ(NLUOUT,*) X2MAXC(2,I),X2MINC(2,I),X2MODC(2,I) 
        ENDDO
        DO I=1,I2IC+1 
          READ(NLUOUT,*) (W2IJC(I,J),J=1,J2JC) 
        ENDDO
        DO J=1,J2JC+1 
          READ(NLUOUT,*) (W2JKC(J,K),K=1,K2KC) 
        ENDDO
        CLOSE(NLUOUT)


END IF

IF (CMINERAL == 'TABUL') THEN
IF(.NOT.ALLOCATED(rhi)) ALLOCATE(rhi(16))
IF(.NOT.ALLOCATED(tempi)) ALLOCATE(tempi(11))
IF(.NOT.ALLOCATED(zsu)) ALLOCATE(zsu(22))
IF(.NOT.ALLOCATED(znh)) ALLOCATE(znh(22))
IF(.NOT.ALLOCATED(zni)) ALLOCATE(zni(22))
IF(.NOT.ALLOCATED(zf)) ALLOCATE(zf(16,11,22,22,22,3))
OPEN(NLUOUT,FILE="AEROMIN_NEW",STATUS="OLD") 

WRITE(*,*) 'LOADING MINERAL AEROSOL DATA ...'
DO i=1,nh
READ(NLUOUT,*) rhi(i)
ENDDO
DO i=1,nt
READ(NLUOUT,*) tempi(i)
ENDDO
DO i=1,nc
READ(NLUOUT,*) zsu(i)
ENDDO
DO i=1,nc
READ(NLUOUT,*) znh(i)
ENDDO
DO i=1,nc
READ(NLUOUT,*) zni(i)
ENDDO
DO i=1,nh
DO j=1,nt
DO k=1,nc
DO l=1,nc
DO m=1,nc
  READ (NLUOUT,*) zf(i,j,k,l,m,1:3)
ENDDO
ENDDO
ENDDO
ENDDO
ENDDO
CLOSE(NLUOUT)
WRITE(*,*) 'END LOADING'
ENDIF



IF (LHOOK) CALL DR_HOOK('CH_AER_MOD_INIT',1,ZHOOK_HANDLE)
END SUBROUTINE CH_AER_MOD_INIT
