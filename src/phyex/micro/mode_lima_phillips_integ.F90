MODULE MODE_LIMA_PHILLIPS_INTEG
  IMPLICIT NONE
CONTAINS
!     ######################################################################
  SUBROUTINE LIMA_PHILLIPS_INTEG (LIMAP, CST, ISIZE, PZT, PSI, PSI0, PSW, PZY, P_FRAC_ACT)
!     ######################################################################
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the fraction of each aerosol 
!!    species (DM1, DM2, BC, O) that may be activated, following Phillips (2008)
!!
!!
!!    REFERENCE
!!    ---------
!!
!!      Phillips et al., 2008: An empirical parameterization of heterogeneous
!!        ice nucleation for multiple chemical species of aerosols, J. Atmos. Sci. 
!!
!!    AUTHOR
!!    ------
!!      J.-M. Cohard     * Laboratoire d'Aerologie*
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!      S.    Berthet    * Laboratoire d'Aerologie*
!!      B.    Vié        * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original             ??/??/13 
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST,            ONLY: CST_T
USE MODE_LIMA_FUNCTIONS,  ONLY : DELTA, DELTA_VEC
USE MODD_PARAM_LIMA, ONLY:PARAM_LIMA_T
USE YOMHOOK, ONLY:LHOOK, DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(PARAM_LIMA_T),INTENT(IN)::LIMAP
TYPE(CST_T),                    INTENT(IN)    :: CST
INTEGER,                        INTENT(IN)    :: ISIZE
REAL, DIMENSION(ISIZE),         INTENT(IN)    :: PZT
REAL, DIMENSION(ISIZE),         INTENT(IN)    :: PSI
REAL, DIMENSION(ISIZE,LIMAP%NSPECIE), INTENT(IN)    :: PSI0
REAL, DIMENSION(ISIZE),         INTENT(IN)    :: PSW
REAL, DIMENSION(ISIZE),         INTENT(IN)    :: PZY
REAL, DIMENSION(ISIZE,LIMAP%NSPECIE), INTENT(OUT)   :: P_FRAC_ACT
!
!*       0.2   Declarations of local variables :
!
INTEGER :: ISPECIE, IL, IL2
REAL :: ZB
!
REAL, DIMENSION(:), ALLOCATABLE :: ZZX,      & ! Work array
                                   ZFACTOR,  &
                                   ZSUBSAT,  &
                                   ZEMBRYO
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
LOGICAL, DIMENSION(:),   ALLOCATABLE :: GINTEG ! Mask to integrate over the
                                               ! AP size spectrum
!
!
!-------------------------------------------------------------------------------
!
P_FRAC_ACT(:,:)=0.
!
IF (LHOOK) CALL DR_HOOK('LIMA_PHILLIPS_INTEG', 0, ZHOOK_HANDLE)
DO ISPECIE = 1, LIMAP%NSPECIE        ! = 4 = {DM1, DM2, BC, O} respectively  
!
   ALLOCATE(ZZX     (SIZE(PZT)) ) ; ZZX(:) = 0.0
   ALLOCATE(ZFACTOR (SIZE(PZT)) )
   ALLOCATE(ZSUBSAT (SIZE(PZT)) )    
   ALLOCATE(ZEMBRYO (SIZE(PZT)) )
   ALLOCATE(GINTEG  (SIZE(PZT)) )

! Compute log in advance for efficiency
   ZB = LOG(0.1E-6/LIMAP%XMDIAM_IFN(ISPECIE))/(SQRT(2.0)*LOG(LIMAP%XSIGMA_IFN(ISPECIE)))
! ZFACTOR = f_c
   ZFACTOR(:) = DELTA(1.,LIMAP%XH(ISPECIE),PZT(:),LIMAP%XT0(ISPECIE),LIMAP%XT0(ISPECIE)+LIMAP%XDT0(ISPECIE))         &
        * DELTA_VEC(0.,1.,PSI(:),PSI0(:,ISPECIE),PSI0(:,ISPECIE)+LIMAP%XDSI0(ISPECIE)) / LIMAP%XGAMMA
! ZSUBSAT = H_X
   ZSUBSAT(:) = MIN(ZFACTOR(:)+(1.0-ZFACTOR(:))*DELTA(0.,1.,PSW(:),LIMAP%XSW0,1.) , 1.0)
! ZEMBRYO = µ_X/(pi*(D_X)**2) = A
   ZEMBRYO(:) = ZSUBSAT(:)*DELTA(1.,0.,PZT(:),LIMAP%XTX1(ISPECIE),LIMAP%XTX2(ISPECIE))          &
        * LIMAP%XFRAC_REF(ISPECIE)*PZY(:)/LIMAP%XAREA1(ISPECIE) 
!
! For T warmer than -35°C, the integration is approximated with µ_X << 1
! Error function : GAMMA_INC(1/2, x**2) = ERF(x) !!! for x>=0 !!!
!   
!   WHERE (PZT(:)>(CST%XTT-35.) .AND. ZEMBRYO(:)>1.0E-8)   
!      ZZX(:) = ZZX(:) + ZEMBRYO(:) * CST%XPI * (LIMAP%XMDIAM_IFN(ISPECIE))**2 / 2.0           &
!           * EXP(2*(LOG(LIMAP%XSIGMA_IFN(ISPECIE)))**2)                                   &
!           * (1.0+GAMMA_INC(0.5,(SQRT(2.0)*LOG(LIMAP%XSIGMA_IFN(ISPECIE))-ZB)**2))     
!   END WHERE

   DO IL = 1, SIZE(PZT)
      IF (PZT(IL)>(CST%XTT-35.) .AND. ZEMBRYO(IL)>1.0E-8) THEN
         ZZX(IL) = ZZX(IL) + ZEMBRYO(IL) * CST%XPI * (LIMAP%XMDIAM_IFN(ISPECIE))**2 / 2.0        &
              * EXP(2*(LOG(LIMAP%XSIGMA_IFN(ISPECIE)))**2)                                   &
              * (1.0+SIGN(1.,SQRT(2.0)*LOG(LIMAP%XSIGMA_IFN(ISPECIE))-ZB)*LIMAP%XGINC_IFN(ISPECIE))
      END IF
   ENDDO

!
! For other T, integration between 0 and infinity is made with a Gauss-Hermite
! quadrature method and integration between 0 and 0.1 uses e(x) ~ 1+x+O(x**2)
! Beware : here, weights are normalized : LIMAP%XWEIGHT = wi/sqrt(pi)
!
   GINTEG(:) = PZT(:)<=(CST%XTT-35.) .AND. PSI(:)>1.0 .AND. ZEMBRYO(:)>1.0E-8
!
   DO IL = 1, LIMAP%NDIAM
      DO IL2 = 1, SIZE(GINTEG)
         IF (GINTEG(IL2)) THEN
            ZZX(IL2) = ZZX(IL2) - LIMAP%XWEIGHT(IL)*EXP(-ZEMBRYO(IL2)*CST%XPI*(LIMAP%XMDIAM_IFN(ISPECIE))**2 & 
                 * EXP(2.0*SQRT(2.0)*LOG(LIMAP%XSIGMA_IFN(ISPECIE)) * LIMAP%XABSCISS(IL)) ) 
         END IF
      ENDDO
   ENDDO
!
!   DO IL2 = 1, SIZE(GINTEG)
!      IF (GINTEG(IL2)) THEN
!         ZZX(IL2) = ZZX(IL2) + 0.5* CST%XPI*ZEMBRYO(IL2)*(LIMAP%XMDIAM_IFN(ISPECIE))**2                &
!              * (1.0-( 1.0-GAMMA_INC(0.5,(SQRT(2.0)*LOG(LIMAP%XSIGMA_IFN(ISPECIE))-ZB)**2))  &
!              * EXP( 2.0*(LOG(LIMAP%XSIGMA_IFN(ISPECIE)))**2)  )
!      END IF
!   ENDDO
   DO IL2 = 1, SIZE(GINTEG)
      IF (GINTEG(IL2)) THEN
         ZZX(IL2) = 1 + ZZX(IL2)  &
              - ( 0.5* CST%XPI*ZEMBRYO(IL2)*(LIMAP%XMDIAM_IFN(ISPECIE))**2 * EXP( 2.0*(LOG(LIMAP%XSIGMA_IFN(ISPECIE)))**2)   &
              * ( 1.0-SIGN(1.,SQRT(2.0)*LOG(LIMAP%XSIGMA_IFN(ISPECIE))-ZB)*LIMAP%XGINC_IFN(ISPECIE)) )
      END IF
   ENDDO
! 
   P_FRAC_ACT(:,ISPECIE)=ZZX(:)
!
   DEALLOCATE(ZZX)
   DEALLOCATE(ZFACTOR)
   DEALLOCATE(ZSUBSAT)    
   DEALLOCATE(ZEMBRYO)
   DEALLOCATE(GINTEG)
!
ENDDO
!
IF (LHOOK) CALL DR_HOOK('LIMA_PHILLIPS_INTEG', 1, ZHOOK_HANDLE)
END SUBROUTINE LIMA_PHILLIPS_INTEG
END MODULE MODE_LIMA_PHILLIPS_INTEG
