!     ######spl
      SUBROUTINE CH_SET_PHOTO_RATES(PTIME,PCONC,KL,TPM,KMI,KOUT,KVERB,KVECNPT,KVECMASK,KEQ,PJVALUES)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!!    #############################################################
!! This code has been created automatically by preprocessor m10,
!! version: 9.7, copyright 1995-1999 by Meteo France/Universite Paul Sabatier.
!! Please report all bugs to K. Suhre (Lab. d'Aerologie UPS/CNRS).
!!*** *MODD_CH_SET_PHOTO_RATES*
!!
!!    PURPOSE
!!    -------
!       set or calculate photolysis rates
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
!!    Modified 29/03/01: Vectorization + nesting (C. Mari)
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
USE MODD_CH_M9,        ONLY : METEOTRANSTYPE
USE MODI_CH_ALLOCATE_TACCS
!     USER DEFINED FUNCTIONS
USE MODI_TROE
USE MODI_TROE_EQUIL
!!
!!    EXPLICIT ARGUMENTS
!!    ------------------
IMPLICIT NONE
REAL,    INTENT(IN)                              :: PTIME
INTEGER, INTENT(IN)                              :: KVECNPT,KL,KEQ,KMI
INTEGER, DIMENSION(:,:), INTENT(IN)              :: KVECMASK
REAL,    INTENT(IN),  DIMENSION(KVECNPT,KEQ)     :: PCONC
TYPE(METEOTRANSTYPE), DIMENSION(KVECNPT), INTENT(IN) :: TPM
INTEGER, INTENT(IN)                              :: KOUT,KVERB
REAL,DIMENSION(:,:,:,:), INTENT(IN) :: PJVALUES    ! Tuv coefficient
!!
! /BEGIN_SET_PHOTO_RATES/
 ! parameter for use by subroutine JVALUES,
 ! contains the actual photolysis rates
 REAL, DIMENSION(KVECNPT,21) :: ZRATESIO ! TUV photolysis rates at one level
 REAL, DIMENSION(KVECNPT,17) :: ZRATES   ! photolysis rates of RACM (vector)
 INTEGER                     :: JITPK    ! loop counter for J-Value transfer
 INTEGER                     :: IDTI,IDTJ
 INTEGER                     :: JITPKPLUS
 INTEGER, DIMENSION(KVECNPT) :: ITABI, ITABJ
 INTEGER, DIMENSION(KVECNPT) :: IMODELLEVEL
 TYPE(CCSTYPE), POINTER      :: TPK
 !
 ! Normally allocated in CH_SET_RATES but who knows ?
 REAL(KIND=JPRB) :: ZHOOK_HANDLE
 IF (LHOOK) CALL DR_HOOK('CH_SET_PHOTO_RATES',0,ZHOOK_HANDLE)
 IF (.NOT. ASSOCIATED(TACCS(KMI)%NVERB)) THEN
   CALL CH_ALLOCATE_TACCS(KMI, KVECNPT)
 END IF
 !
 ! TPK is set for current model 
 TPK=>TACCS(KMI)
 
 ! calculation of photolysis rates and transfer into local variables
 !
 IDTI=KVECMASK(2,KL)-KVECMASK(1,KL)+1
 IDTJ=KVECMASK(4,KL)-KVECMASK(3,KL)+1
 DO JITPK = 0, KVECNPT-1
 !
   JITPKPLUS=JITPK+1
   ITABI(JITPKPLUS)=JITPK-IDTI*(JITPK/IDTI)+KVECMASK(1,KL)
   ITABJ(JITPKPLUS)=JITPK/IDTI-IDTJ*(JITPK/(IDTI*IDTJ))+KVECMASK(3,KL)
! modification for VPP optimization
 ZRATESIO(JITPK+1,:) = PJVALUES(ITABI(JITPKPLUS),ITABJ(JITPKPLUS),TPK%MODELLEVEL(JITPK+1),:)
 ENDDO 
 !
 DO JITPK = 0, KVECNPT-1
 !
 ! associate TUV J-Values to ReLACS J-Values
 !  
   ZRATES(JITPK+1, 1) = ZRATESIO(JITPK+1,2) 
   ZRATES(JITPK+1, 2) = ZRATESIO(JITPK+1,3) 
   ZRATES(JITPK+1, 3) = ZRATESIO(JITPK+1,4) 
   ZRATES(JITPK+1, 4) = ZRATESIO(JITPK+1,9) 
   ZRATES(JITPK+1, 5) = ZRATESIO(JITPK+1,10) 
   ZRATES(JITPK+1, 6) = ZRATESIO(JITPK+1,11) 
   ZRATES(JITPK+1, 7) = ZRATESIO(JITPK+1,5) 
   ZRATES(JITPK+1, 8) = ZRATESIO(JITPK+1,6) 
   ZRATES(JITPK+1, 9) = ZRATESIO(JITPK+1,12) 
   ZRATES(JITPK+1, 10) = ZRATESIO(JITPK+1,14) 
   ZRATES(JITPK+1, 11) = ZRATESIO(JITPK+1,13) 
   ZRATES(JITPK+1, 12) = ZRATESIO(JITPK+1,17) 
   ZRATES(JITPK+1, 13) = ZRATESIO(JITPK+1,15)  
   ZRATES(JITPK+1, 14) =  0.962055 *ZRATESIO(JITPK+1,15)+&
                     &  1.06247E-02 *ZRATESIO(JITPK+1,12)  
   ZRATES(JITPK+1, 15) =  ZRATESIO(JITPK+1,20)  
   ZRATES(JITPK+1, 16) =  3.16657 *ZRATESIO(JITPK+1,15)&
                     &+ 0.372446 *ZRATESIO(JITPK+1,15)&
                     &+ 8.42257 *ZRATESIO(JITPK+1,15)&
                     &+ 207.5913 *ZRATESIO(JITPK+1,20)&
                     &+ 0.0 *ZRATESIO(JITPK+1,20)&
                     &+ 8.44837E-04 *ZRATESIO(JITPK+1,20)   
   ZRATES(JITPK+1, 17) = ZRATESIO(JITPK+1,16)
 !
 END DO
 !
! /END_SET_PHOTO_RATES/
 TPK%K001=ZRATES(:,001)
 TPK%K002=ZRATES(:,002)
 TPK%K003=ZRATES(:,003)
 TPK%K004=ZRATES(:,004)
 TPK%K005=ZRATES(:,005)
 TPK%K006=ZRATES(:,006)
 TPK%K007=ZRATES(:,007)
 TPK%K008=ZRATES(:,008)
 TPK%K009=ZRATES(:,009)
 TPK%K010=ZRATES(:,010)
 TPK%K011=ZRATES(:,011)
 TPK%K012=ZRATES(:,012)
 TPK%K013=ZRATES(:,013)
 TPK%K014=ZRATES(:,014)
 TPK%K015=ZRATES(:,015)
 TPK%K016=ZRATES(:,016)
 TPK%K017=ZRATES(:,017)
TPK%NOUT  = KOUT
TPK%NVERB = KVERB
IF (LHOOK) CALL DR_HOOK('CH_SET_PHOTO_RATES',1,ZHOOK_HANDLE)
RETURN
END SUBROUTINE CH_SET_PHOTO_RATES
