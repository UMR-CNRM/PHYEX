!     ######spl
SUBROUTINE CH_EXQSSA(PTSIMUL, PDTACT, PCONC, PNEWCONC, KEQ, KVECNPT, KMI, &
                       PSLOW, PFAST, PDTMAX)
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!!    #########################################################################
!!
!!*** *CH_EXQSSA*
!!
!!    PURPOSE
!!    -------
!    solve one time-step of the chemical differential equation d/dt C = P - L*C
!!
!!**  METHOD
!!    ------
!!    Quasy steady state approximation (EXQSSA):
!!      C_n+1 = Q(C_n,h) = P_n / L_n + ( C_n - P_n / L_n ) exp( - L_n Dt)
!!
!!    Extrapolation
!!    C_n+1(h) = Q(C_n,h)
!!    C_n+0.5(h/2)=Q(C_n,h/2)
!!    C_n+1(h/2)=Q(C_n+0.5,h/2)
!!    C_n+1 = 2 C_n+1(h/2) - C_n+1(h)
!!
!!    REFERENCE
!!    ---------
!!    O. Hertel et al ,
!!    Test of two Numerical schemes for use in atmospheric transport 
!!    chemistry models,
!!    Atmospheric environment, 16, 2591--2611, 1993.
!!
!!    A. Sandu et al
!!    Benchmarking Stiff Ode Solvers For Atmospheric Chemistry Problems I:
!!    Implicit versus Explicit
!!    submitted to Atmospheric Environment 26/01/1996
!!
!!    AUTHOR
!!    ------
!!    A. FASSI FIHRI    *Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!    Original 13/02/96
!!    10/04/96 (K. Suhre) take out performance control statements
!!    31/07/96 (K. Suhre) restructured
!!    05/03/98 (V. Crassier) vectorized with fixed timestep
!!    31/03/99 (K. Suhre) change minor details (i.e. numerical pbs with
!!                        LOSS=0 species under Absoft and HP f90 compilers
!!                        and control of correct time stepping if user input
!!                        is uncorrect) and restructure
!!    01/08/01 (C. Mari)   add PSLOW, PFAST, PDTMAX arguments
!!    01/12/03 (D. Gazen)  change Chemical scheme interface
!!
!!    EXTERNAL
!!    --------
USE MODI_CH_PRODLOSS
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!        -----------------
IMPLICIT NONE
!
!*       0.1  declaration of arguments
!
REAL,    INTENT(IN) :: PTSIMUL  ! time of simulation
REAL,    INTENT(IN) :: PDTACT   ! actual time-step
INTEGER, INTENT(IN) :: KEQ      ! dimension of the problem to solve
INTEGER, INTENT(IN) :: KVECNPT
REAL, INTENT(IN),  DIMENSION(KVECNPT,KEQ) :: PCONC    
                              ! concentration vector at PTSIMUL
REAL, INTENT(OUT), DIMENSION(KVECNPT,KEQ) :: PNEWCONC 
                              ! solution at PTSIMUL + PDTACT
INTEGER, INTENT(IN) :: KMI    ! model number
                              ! reac. rates, auxiliary variables
REAL,    INTENT(IN) :: PSLOW, PFAST, PDTMAX
!
!*       0.2  declaration of local variables
!
REAL, DIMENSION(KVECNPT,KEQ)     :: ZCONC, ZPROD, ZLOSS
REAL, DIMENSION(KVECNPT,KEQ)     :: ZCONC1, ZCONC2
REAL, DIMENSION(KVECNPT,KEQ)     :: ZCONC1A, ZCONC2A
REAL, DIMENSION(KVECNPT,KEQ)     :: ZCONC1B, ZCONC2B
REAL, DIMENSION(KVECNPT,KEQ)     :: ZCONC1C, ZCONC2C
REAL, DIMENSION(KVECNPT,KEQ)     :: ZTMP, ZRATIO
REAL                             :: ZTIME, ZSTEP
INTEGER                          :: JI
REAL, DIMENSION(KVECNPT,KEQ)     :: ZA, ZB, ZC
REAL, PARAMETER                  :: ZEPS = 1E-14 ! to avoid division by zero
INTEGER                          :: IITER  ! number of internal timesteps
!
!------------------------------------------------------------------------------
!
!*       1.   INITIALIZE SOME VARIABLES
!        -------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('CH_EXQSSA',0,ZHOOK_HANDLE)
ZTIME  = PTSIMUL 
ZCONC  = PCONC
IITER  = MAX(1,INT(1.E-3+(PDTACT/PDTMAX)))
ZSTEP  = PDTACT / IITER
!
!*       2.   INTERATION LOOP
!        --------------------
!
DO JI=1, IITER
!
!*       2.1  calculate production and loss terms
!
  ZCONC(:,:) = MAX(ZCONC(:,:), 1E-10)
  CALL CH_PRODLOSS(ZTIME,ZCONC,ZPROD,ZLOSS,KMI,KVECNPT,KEQ)
!
! elimiate negative loss rates that may be inherited from the transport scheme
!
  ZLOSS(:,:) = MAX(0.,ZLOSS(:,:))
  ZA(:,:) = 0.
  ZB(:,:) = 0.
  WHERE ( ZLOSS .LE. 1./(PSLOW*ZSTEP) ) ZA(:,:) = 1.
  WHERE ( ZLOSS .GT. 1./(PFAST*ZSTEP) ) ZB(:,:) = 1.
  ZC(:,:) = 1. - ZA(:,:) - ZB(:,:)
!        
!*       2.2  first guess
!
! especes lentes ZLOSS<1./(PSLOW*ZSTEP)
! calcul des constantes multiplicatives
!***************************************

  ZTMP(:,:)    = ZSTEP*(ZPROD(:,:) - ZLOSS(:,:)*ZCONC(:,:)) 
  ZCONC1A(:,:) = ZCONC(:,:) + ZTMP(:,:) 
  ZCONC2A(:,:) = ZCONC(:,:) + 0.5 * ZTMP(:,:)
  
! especes rapides ZLOSS>1./(PFAST*ZSTEP)
! calcul des constantes multiplicatives
!***************************************
  
  ZCONC1B(:,:) = ZPROD(:,:) / MAX(ZEPS, ZLOSS(:,:))
  ZCONC2B(:,:) = ZCONC1B(:,:)
  
! especes intermediaires
! calcul des constantes multiplicatives
!***************************************

  ZRATIO(:,:)  = ZPROD(:,:) / MAX(ZEPS, ZLOSS(:,:))
  ZTMP(:,:)    = EXP(-ZLOSS(:,:)*ZSTEP*0.5)
  ZCONC1C(:,:) = ZTMP(:,:) * ZTMP(:,:) * ( ZCONC(:,:) - ZRATIO(:,:) ) &
                 + ZRATIO(:,:)
  ZCONC2C(:,:) = ZTMP(:,:) * ( ZCONC(:,:) - ZRATIO(:,:) ) + ZRATIO(:,:)

! calcul des concentrations intermediaires
!*****************************************

  ZCONC1(:,:)= ZA(:,:)*ZCONC1A(:,:) &
             + ZB(:,:)*ZCONC1B(:,:) &
             + ZC(:,:)*ZCONC1C(:,:)
  ZCONC2(:,:)= ZA(:,:)*ZCONC2A(:,:) &
             + ZB(:,:)*ZCONC2B(:,:) &
             + ZC(:,:)*ZCONC2C(:,:)

!
!*      2.3  calculate new production and loss terms
!
  ZCONC2(:,:) = MAX(ZCONC2(:,:), 1E-10)
  CALL CH_PRODLOSS(ZTIME+0.5*ZSTEP,ZCONC2,ZPROD,ZLOSS,KMI,KVECNPT,KEQ)
!
! elimiate negative loss rates that may be inherited from the transport scheme
!
  ZLOSS(:,:) = MAX(0.,ZLOSS(:,:))
  ZA(:,:) = 0.
  ZB(:,:) = 0.
  WHERE ( ZLOSS .LE. 1./(PSLOW*ZSTEP) ) ZA(:,:) = 1.
  WHERE ( ZLOSS .GT. 1./(PFAST*ZSTEP) ) ZB(:,:) = 1.
  ZC(:,:) = 1. - ZA(:,:) - ZB(:,:)
!
! especes lentes ZLOSS<1./(PSLOW*ZSTEP)
! calcul des constantes multiplicatives
!***************************************
  
  ZTMP(:,:) = ZSTEP*(ZPROD(:,:) - ZLOSS(:,:)*ZCONC(:,:)) 
  ZCONC2A(:,:) = ZCONC2(:,:) + 0.5 * ZTMP(:,:)
!
! especes rapides ZLOSS>1./(PFAST*ZSTEP)
! calcul des constantes multiplicatives
!***************************************

  ZCONC2B(:,:) = ZPROD(:,:) / MAX(ZEPS, ZLOSS(:,:))
  
!
! especes intermediaires
! calcul des constantes multiplicatives
!***************************************
  
  ZRATIO(:,:)  = ZPROD(:,:) /  MAX(ZEPS, ZLOSS(:,:))
  ZTMP(:,:)    = EXP(-ZLOSS(:,:)*ZSTEP*0.5)
  ZCONC2C(:,:) = ZTMP(:,:) * ( ZCONC2(:,:) - ZRATIO(:,:) ) + ZRATIO(:,:) 
  
! calcul des concentrations intermediaires
!*****************************************

  ZCONC2(:,:) = ZA(:,:)*ZCONC2A(:,:) &
              + ZB(:,:)*ZCONC2B(:,:) &
              + ZC(:,:)*ZCONC2C(:,:)
  
! IF ZSTEP*ZLOSS < 1. --> A=1 else B=1
!*************************************

  ZA(:,:) = 0.
  WHERE ( (ZSTEP*ZLOSS) .LT. 1. ) ZA(:,:) = 1.
  ZB(:,:) = 1. - ZA(:,:)

  ZCONC(:,:) = ZA(:,:)*MAX(0.,2.*ZCONC2(:,:)-ZCONC1(:,:)) &
             + ZB(:,:)*MAX(0.,ZCONC2(:,:))

  ZTIME=ZTIME+ZSTEP

ENDDO
!
PNEWCONC = ZCONC
!     
IF (LHOOK) CALL DR_HOOK('CH_EXQSSA',1,ZHOOK_HANDLE)
END SUBROUTINE CH_EXQSSA
