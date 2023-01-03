!      ###############################
       MODULE MODI_LIMA_PHILLIPS_INTEG
!      ###############################
!
INTERFACE
      SUBROUTINE LIMA_PHILLIPS_INTEG (ZZT, ZSI, ZSI0, ZSW, ZZY, Z_FRAC_ACT)
!
REAL, DIMENSION(:),   INTENT(IN)    :: ZZT
REAL, DIMENSION(:),   INTENT(IN)    :: ZSI
REAL, DIMENSION(:,:), INTENT(IN)    :: ZSI0
REAL, DIMENSION(:),   INTENT(IN)    :: ZSW
REAL, DIMENSION(:),   INTENT(IN)    :: ZZY
REAL, DIMENSION(:,:), INTENT(INOUT) :: Z_FRAC_ACT
!
END SUBROUTINE LIMA_PHILLIPS_INTEG
END INTERFACE
END MODULE MODI_LIMA_PHILLIPS_INTEG
!
!     ######################################################################
      SUBROUTINE LIMA_PHILLIPS_INTEG (ZZT, ZSI, ZSI0, ZSW, ZZY, Z_FRAC_ACT)
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
USE MODD_CST,             ONLY : XTT, XPI
USE MODD_PARAM_LIMA,      ONLY : XMDIAM_IFN, XSIGMA_IFN, NSPECIE, XFRAC_REF, &
                                 XH, XAREA1, XGAMMA, XABSCISS, XWEIGHT, NDIAM,     &
                                 XT0, XDT0, XDSI0, XSW0, XTX1, XTX2
USE MODI_LIMA_FUNCTIONS,  ONLY : DELTA, DELTA_VEC
USE MODI_GAMMA_INC
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
REAL, DIMENSION(:),   INTENT(IN)    :: ZZT
REAL, DIMENSION(:),   INTENT(IN)    :: ZSI
REAL, DIMENSION(:,:), INTENT(IN)    :: ZSI0
REAL, DIMENSION(:),   INTENT(IN)    :: ZSW
REAL, DIMENSION(:),   INTENT(IN)    :: ZZY
REAL, DIMENSION(:,:), INTENT(INOUT) :: Z_FRAC_ACT
!
!*       0.2   Declarations of local variables :
!
INTEGER :: JSPECIE, JL, JL2
REAL :: XB
!
REAL, DIMENSION(:), ALLOCATABLE :: ZZX,      & ! Work array
                                   ZFACTOR,  &
                                   ZSUBSAT,  &
                                   ZEMBRYO
!
LOGICAL, DIMENSION(:),   ALLOCATABLE :: GINTEG ! Mask to integrate over the
                                               ! AP size spectrum
!
!
!-------------------------------------------------------------------------------
!
!
DO JSPECIE = 1, NSPECIE        ! = 4 = {DM1, DM2, BC, O} respectively  
!
   ALLOCATE(ZZX     (SIZE(ZZT)) ) ; ZZX(:) = 0.0
   ALLOCATE(ZFACTOR (SIZE(ZZT)) )
   ALLOCATE(ZSUBSAT (SIZE(ZZT)) )    
   ALLOCATE(ZEMBRYO (SIZE(ZZT)) )
   ALLOCATE(GINTEG  (SIZE(ZZT)) )

! Compute log in advance for efficiency
   XB = LOG(0.1E-6/XMDIAM_IFN(JSPECIE))/(SQRT(2.0)*LOG(XSIGMA_IFN(JSPECIE)))
! ZFACTOR = f_c
   ZFACTOR(:) = DELTA(1.,XH(JSPECIE),ZZT(:),XT0(JSPECIE),XT0(JSPECIE)+XDT0(JSPECIE))         &
        * DELTA_VEC(0.,1.,ZSI(:),ZSI0(:,JSPECIE),ZSI0(:,JSPECIE)+XDSI0(JSPECIE)) / XGAMMA
! ZSUBSAT = H_X
   ZSUBSAT(:) = MIN(ZFACTOR(:)+(1.0-ZFACTOR(:))*DELTA(0.,1.,ZSW(:),XSW0,1.) , 1.0)
! ZEMBRYO = µ_X/(pi*(D_X)**2) = A
   ZEMBRYO(:) = ZSUBSAT(:)*DELTA(1.,0.,ZZT(:),XTX1(JSPECIE),XTX2(JSPECIE))          &
        * XFRAC_REF(JSPECIE)*ZZY(:)/XAREA1(JSPECIE) 
!
! For T warmer than -35°C, the integration is approximated with µ_X << 1
! Error function : GAMMA_INC(1/2, x**2) = ERF(x) !!! for x>=0 !!!
!   
!   WHERE (ZZT(:)>(XTT-35.) .AND. ZEMBRYO(:)>1.0E-8)   
!      ZZX(:) = ZZX(:) + ZEMBRYO(:) * XPI * (XMDIAM_IFN(JSPECIE))**2 / 2.0           &
!           * EXP(2*(LOG(XSIGMA_IFN(JSPECIE)))**2)                                   &
!           * (1.0+GAMMA_INC(0.5,(SQRT(2.0)*LOG(XSIGMA_IFN(JSPECIE))-XB)**2))     
!   END WHERE

   DO JL = 1, SIZE(ZZT)
      IF (ZZT(JL)>(XTT-35.) .AND. ZEMBRYO(JL)>1.0E-8) THEN
         ZZX(JL) = ZZX(JL) + ZEMBRYO(JL) * XPI * (XMDIAM_IFN(JSPECIE))**2 / 2.0        &
              * EXP(2*(LOG(XSIGMA_IFN(JSPECIE)))**2)                                   &
              * (1.0+SIGN(1.,SQRT(2.0)*LOG(XSIGMA_IFN(JSPECIE))-XB)*GAMMA_INC(0.5,(SQRT(2.0)*LOG(XSIGMA_IFN(JSPECIE))-XB)**2))     
      END IF
   ENDDO

!
! For other T, integration between 0 and infinity is made with a Gauss-Hermite
! quadrature method and integration between 0 and 0.1 uses e(x) ~ 1+x+O(x**2)
! Beware : here, weights are normalized : XWEIGHT = wi/sqrt(pi)
!
   GINTEG(:) = ZZT(:)<=(XTT-35.) .AND. ZSI(:)>1.0 .AND. ZEMBRYO(:)>1.0E-8
!
   DO JL = 1, NDIAM
      DO JL2 = 1, SIZE(GINTEG)
         IF (GINTEG(JL2)) THEN
            ZZX(JL2) = ZZX(JL2) - XWEIGHT(JL)*EXP(-ZEMBRYO(JL2)*XPI*(XMDIAM_IFN(JSPECIE))**2 & 
                 * EXP(2.0*SQRT(2.0)*LOG(XSIGMA_IFN(JSPECIE)) * XABSCISS(JL)) ) 
         END IF
      ENDDO
   ENDDO
!
!   DO JL2 = 1, SIZE(GINTEG)
!      IF (GINTEG(JL2)) THEN
!         ZZX(JL2) = ZZX(JL2) + 0.5* XPI*ZEMBRYO(JL2)*(XMDIAM_IFN(JSPECIE))**2                &
!              * (1.0-( 1.0-GAMMA_INC(0.5,(SQRT(2.0)*LOG(XSIGMA_IFN(JSPECIE))-XB)**2))  &
!              * EXP( 2.0*(LOG(XSIGMA_IFN(JSPECIE)))**2)  )
!      END IF
!   ENDDO
   DO JL2 = 1, SIZE(GINTEG)
      IF (GINTEG(JL2)) THEN
         ZZX(JL2) = 1 + ZZX(JL2)  &
              - ( 0.5* XPI*ZEMBRYO(JL2)*(XMDIAM_IFN(JSPECIE))**2 * EXP( 2.0*(LOG(XSIGMA_IFN(JSPECIE)))**2)   &
              * ( 1.0-SIGN(1.,SQRT(2.0)*LOG(XSIGMA_IFN(JSPECIE))-XB)*GAMMA_INC(0.5,(SQRT(2.0)*LOG(XSIGMA_IFN(JSPECIE))-XB)**2)) )
      END IF
   ENDDO
! 
   Z_FRAC_ACT(:,JSPECIE)=ZZX(:)
!
   DEALLOCATE(ZZX)
   DEALLOCATE(ZFACTOR)
   DEALLOCATE(ZSUBSAT)    
   DEALLOCATE(ZEMBRYO)
   DEALLOCATE(GINTEG)
!
ENDDO
!
END SUBROUTINE LIMA_PHILLIPS_INTEG
