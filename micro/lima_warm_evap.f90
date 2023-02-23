!MNH_LIC Copyright 2013-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!      ##########################
       MODULE MODI_LIMA_WARM_EVAP
!      ##########################
!
INTERFACE
      SUBROUTINE LIMA_WARM_EVAP (PTSTEP, KMI,                                &
                                 PRHODREF, PEXNREF, PPABST, ZT,              &
                                 ZWLBDC3, ZWLBDC, ZWLBDR3, ZWLBDR,           &
                                 PRVT, PRCT, PRRT, PCRT,                     &
                                 PRVS, PRCS, PRRS, PCCS, PCRS, PTHS,         &
                                 PEVAP3D)
!
REAL,                     INTENT(IN)    :: PTSTEP     ! Double Time step
                                                      ! (single if cold start)
INTEGER,                  INTENT(IN)    :: KMI        ! Model index 
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF   ! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PEXNREF    ! Reference Exner function
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABST     ! abs. pressure at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: ZT         ! Temperature
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: ZWLBDC3    ! Lambda(cloud) **3
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: ZWLBDC     ! Lambda(cloud)
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: ZWLBDR3    ! Lambda(rain) **3
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: ZWLBDR     ! Lambda(rain)
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRVT       ! Water vapor m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRCT       ! Cloud water m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRRT       ! Rain water m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCRT       ! Rain water C. at t
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRVS       ! Water vapor m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRCS       ! Cloud water m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRRS       ! Rain water m.r. source
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCCS       ! Cloud water C. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCRS       ! Rain water C. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PTHS       ! Theta source
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PEVAP3D    ! Rain evap profile
!
      END SUBROUTINE LIMA_WARM_EVAP
END INTERFACE
END MODULE MODI_LIMA_WARM_EVAP
!     #############################################################################
      SUBROUTINE LIMA_WARM_EVAP (PTSTEP, KMI,                                &
                                 PRHODREF, PEXNREF, PPABST, ZT,              &
                                 ZWLBDC3, ZWLBDC, ZWLBDR3, ZWLBDR,           &
                                 PRVT, PRCT, PRRT, PCRT,                     &
                                 PRVS, PRCS, PRRS, PCCS, PCRS, PTHS,         &
                                 PEVAP3D)
!     #############################################################################
!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the raindrop evaporation
!!
!!
!!    AUTHOR
!!    ------
!!      J.-M. Cohard     * Laboratoire d'Aerologie*
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!      S.    Berthet    * Laboratoire d'Aerologie*
!!      B.    Vié        * Laboratoire d'Aerologie*
!!
!!
!!    MODIFICATIONS
!!    -------------
!!      Original             ??/??/13 
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!  P. Wautelet 28/05/2019: move COUNTJV function to tools.f90
! Delbeke/Vie     03/2022 : KHKO option
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
USE MODD_PARAMETERS,      ONLY: JPHEXT, JPVEXT
USE MODD_PARAM_LIMA
USE MODD_PARAM_LIMA_WARM
!
use mode_tools,           only: Countjv
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
REAL,                     INTENT(IN)    :: PTSTEP     ! Double Time step
                                                      ! (single if cold start)
INTEGER,                  INTENT(IN)    :: KMI        ! Model index 
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF   ! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PEXNREF    ! Reference Exner function
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABST     ! abs. pressure at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: ZT         ! Temperature
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: ZWLBDC3    ! Lambda(cloud) **3
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: ZWLBDC     ! Lambda(cloud)
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: ZWLBDR3    ! Lambda(rain) **3
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: ZWLBDR     ! Lambda(rain)
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRVT       ! Water vapor m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRCT       ! Cloud water m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRRT       ! Rain water m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCRT       ! Rain water C. at t
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRVS       ! Water vapor m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRCS       ! Cloud water m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRRS       ! Rain water m.r. source
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCCS       ! Cloud water C. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCRS       ! Rain water C. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PTHS       ! Theta source
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PEVAP3D    ! Rain evap profile
!
!*       0.1   Declarations of local variables :
!
! Packing variables
LOGICAL, DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3)) :: GEVAP, GMICRO 
INTEGER :: IEVAP, IMICRO
INTEGER , DIMENSION(SIZE(GEVAP))   :: I1,I2,I3 ! Used to replace the COUNT
INTEGER                            :: JL       ! and PACK intrinsics 
!
! Packed micophysical variables
REAL, DIMENSION(:)  , ALLOCATABLE :: ZRVT     ! Water vapor m.r. at t 
REAL, DIMENSION(:)  , ALLOCATABLE :: ZRCT     ! Cloud water m.r. at t 
REAL, DIMENSION(:)  , ALLOCATABLE :: ZRRT     ! Rain water m.r. at t 
REAL, DIMENSION(:)  , ALLOCATABLE :: ZCRT     ! rain conc. at t
!
REAL, DIMENSION(:)  , ALLOCATABLE :: ZRVS     ! Water vapor m.r. source
REAL, DIMENSION(:)  , ALLOCATABLE :: ZRRS     ! Rain water m.r. source
REAL, DIMENSION(:)  , ALLOCATABLE :: ZCRS     ! Rain water m.r. source
REAL, DIMENSION(:)  , ALLOCATABLE :: ZTHS     ! Theta source
!
! Other packed variables
REAL, DIMENSION(:)  , ALLOCATABLE :: ZRHODREF ! RHO Dry REFerence
REAL, DIMENSION(:)  , ALLOCATABLE :: ZEXNREF  ! EXNer Pressure REFerence
REAL, DIMENSION(:)  , ALLOCATABLE :: ZZT      ! Temperature
REAL, DIMENSION(:)  , ALLOCATABLE :: ZLBDR    ! Lambda(rain)
!
! Work arrays
REAL, DIMENSION(:), ALLOCATABLE   :: ZZW1, ZZW2, ZZW3, &
                                     ZRTMIN, ZCTMIN, &
                                     ZZLV     ! Latent heat of vaporization at T
!
REAL,    DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3))   &
                                  :: ZW, ZW2, ZRVSAT, ZDR, ZLV
!
!
REAL    :: ZEPS, ZFACT
INTEGER :: IIB, IIE, IJB, IJE, IKB, IKE        ! Physical domain
!
!-------------------------------------------------------------------------------
!
!
!*       1.     PREPARE COMPUTATIONS - PACK
!   	        ---------------------------
!
!
IIB=1+JPHEXT
IIE=SIZE(PRHODREF,1) - JPHEXT
IJB=1+JPHEXT
IJE=SIZE(PRHODREF,2) - JPHEXT
IKB=1+JPVEXT
IKE=SIZE(PRHODREF,3) - JPVEXT
!
ALLOCATE(ZRTMIN(SIZE(XRTMIN)))
ALLOCATE(ZCTMIN(SIZE(XCTMIN)))
ZRTMIN(:) = XRTMIN(:) / PTSTEP
ZCTMIN(:) = XCTMIN(:) / PTSTEP
!
ZEPS= XMV / XMD
ZRVSAT(:,:,:) = ZEPS / (PPABST(:,:,:) * &
                   EXP(-XALPW+XBETAW/ZT(:,:,:)+XGAMW*ALOG(ZT(:,:,:))) - 1.0)
ZLV(:,:,:) = XLVTT + (XCPV-XCL)*(ZT(:,:,:)-XTT)
!
GEVAP(:,:,:) = .FALSE.
GEVAP(IIB:IIE,IJB:IJE,IKB:IKE) =                              &
     PRRS(IIB:IIE,IJB:IJE,IKB:IKE)>ZRTMIN(3) .AND.            &
     PCRS(IIB:IIE,IJB:IJE,IKB:IKE)>ZCTMIN(3) .AND.            &
     PRRT(IIB:IIE,IJB:IJE,IKB:IKE)>ZRTMIN(3) .AND.            &
     PCRT(IIB:IIE,IJB:IJE,IKB:IKE)>ZCTMIN(3) .AND.            &
     PRVT(IIB:IIE,IJB:IJE,IKB:IKE)<ZRVSAT(IIB:IIE,IJB:IJE,IKB:IKE)
!
IEVAP = COUNTJV( GEVAP(:,:,:),I1(:),I2(:),I3(:))
!
IF( IEVAP >= 0 ) THEN
   ALLOCATE(ZRVT(IEVAP))
   ALLOCATE(ZRCT(IEVAP))
   ALLOCATE(ZRRT(IEVAP))
   ALLOCATE(ZCRT(IEVAP))
!
   ALLOCATE(ZRVS(IEVAP))
   ALLOCATE(ZRRS(IEVAP))
   ALLOCATE(ZCRS(IEVAP))
   ALLOCATE(ZTHS(IEVAP))
!
   ALLOCATE(ZLBDR(IEVAP))
!
   ALLOCATE(ZRHODREF(IEVAP))
   ALLOCATE(ZEXNREF(IEVAP))
!
   ALLOCATE(ZZT(IEVAP))  
   ALLOCATE(ZZLV(IEVAP)) 
   ALLOCATE(ZZW1(IEVAP)) 
   DO JL=1,IEVAP
      ZRVT(JL) = PRVT(I1(JL),I2(JL),I3(JL))
      ZRCT(JL) = PRCT(I1(JL),I2(JL),I3(JL))
      ZRRT(JL) = PRRT(I1(JL),I2(JL),I3(JL))
      ZCRT(JL) = PCRT(I1(JL),I2(JL),I3(JL))
      ZRRS(JL) = PRRS(I1(JL),I2(JL),I3(JL))
      ZCRS(JL) = PCRS(I1(JL),I2(JL),I3(JL))
      ZRVS(JL) = PRVS(I1(JL),I2(JL),I3(JL))
      ZTHS(JL) = PTHS(I1(JL),I2(JL),I3(JL))
      ZZT(JL) = ZT(I1(JL),I2(JL),I3(JL))
      ZZW1(JL) = ZRVSAT(I1(JL),I2(JL),I3(JL))
      ZLBDR(JL) = ZWLBDR(I1(JL),I2(JL),I3(JL))
      ZRHODREF(JL) = PRHODREF(I1(JL),I2(JL),I3(JL))
      ZEXNREF(JL)  = PEXNREF(I1(JL),I2(JL),I3(JL))
      ZZLV(JL) = ZLV(I1(JL),I2(JL),I3(JL))
   END DO
!
   ALLOCATE(ZZW2(IEVAP)) 
   ALLOCATE(ZZW3(IEVAP))
!
!
!-------------------------------------------------------------------------------
!
!
!*       2. compute the evaporation of rain drops    
!   	 ----------------------------------------
!
!
   ZZW3(:) = MAX((1.0 - ZRVT(:)/ZZW1(:)),0.0)  ! Subsaturation
!
! Compute the function G(T)
!
   ZZW2(:) = 1. / ( XRHOLW*((((ZZLV(:)/ZZT(:))**2)/(XTHCO*XRV)) +          & ! G
          (XRV*ZZT(:))/(XDIVA*EXP(XALPW-XBETAW/ZZT(:)-XGAMW*ALOG(ZZT(:))))))
!
! Compute the evaporation tendency
!
   IF (LKHKO) THEN
    ZZW2(:) = 3.0 * XCEVAP * ZZW2(:) * (4.*XPI*XRHOLW/(3.*ZRHODREF(:)))**(2./3.) *    &
                               (ZRRT(:))**(1./3.) * (ZCRT(:))**(2./3.) * ZZW3(:)
    ZZW2(:) = MIN(ZZW2(:),ZRRS(:))        
   ELSE
      ZZW2(:) = MIN( ZZW2(:) * ZZW3(:) * ZRRT(:) *        &
                (X0EVAR*ZLBDR(:)**XEX0EVAR + X1EVAR*ZRHODREF(:)**XEX2EVAR* &
                 ZLBDR(:)**XEX1EVAR),ZRRS(:) )
      ZZW2(:) = MAX(ZZW2(:),0.0)
   END IF
!
! Adjust sources
!
   ZRVS(:) = ZRVS(:) + ZZW2(:)
   ZRRS(:) = ZRRS(:) - ZZW2(:)
   ZTHS(:) = ZTHS(:) - ZZW2(:)*ZZLV(:) /                                        &
                    ( ZEXNREF(:)*(XCPD + XCPV*ZRVT(:) + XCL*(ZRCT(:) + ZRRT(:)) ) )
!
!
!-------------------------------------------------------------------------------
!
!
!*       3. Unpack and clean
!   	 -------------------
!
!
   ZW(:,:,:) = PRVS(:,:,:)
   PRVS(:,:,:) = UNPACK( ZRVS(:),MASK=GEVAP(:,:,:),FIELD=ZW(:,:,:) )
   ZW(:,:,:) = PRRS(:,:,:)
   PRRS(:,:,:) = UNPACK( ZRRS(:),MASK=GEVAP(:,:,:),FIELD=ZW(:,:,:) )
   ZW(:,:,:) = PTHS(:,:,:)
   PTHS(:,:,:) = UNPACK( ZTHS(:),MASK=GEVAP(:,:,:),FIELD=ZW(:,:,:) )
   ZW(:,:,:)= PEVAP3D(:,:,:)
   PEVAP3D(:,:,:) = UNPACK( ZZW2(:),MASK=GEVAP(:,:,:),FIELD=ZW(:,:,:) )
!
   IF (LKHKO) THEN
      ZZW2(:) = MIN(ZZW2(:) * ZCRT(:)/ZRRT(:),ZCRS(:))
      ZCRS(:) = ZCRS(:) - ZZW2(:)
      ZW(:,:,:) = PCRS(:,:,:)
      PCRS(:,:,:) = UNPACK( ZCRS(:),MASK=GEVAP(:,:,:),FIELD=ZW(:,:,:) )
   ENDIF

   DEALLOCATE(ZRCT)
   DEALLOCATE(ZRRT)
   DEALLOCATE(ZRVT)
   DEALLOCATE(ZCRT)
   DEALLOCATE(ZRVS)
   DEALLOCATE(ZRRS)
   DEALLOCATE(ZCRS)
   DEALLOCATE(ZTHS)
   DEALLOCATE(ZZLV)
   DEALLOCATE(ZZT)
   DEALLOCATE(ZRHODREF)
   DEALLOCATE(ZEXNREF)
   DEALLOCATE(ZZW1)
   DEALLOCATE(ZZW2)
   DEALLOCATE(ZZW3)
   DEALLOCATE(ZLBDR)
!
!
!-----------------------------------------------------------------------------
!
!
!*       4. Update Nr if:  80 microns < Dr < D_h
!   	 ---------------------------------------
!
!
   IF (LKHKO) THEN
!*             correct negative values for rain
!   	      --------------------------------
!
      WHERE (PRRS(:,:,:)<0.) 
         PRCS(:,:,:) = PRCS(:,:,:)+PRRS(:,:,:)
         PRRS(:,:,:) = 0.
         PCRS(:,:,:) = 0.
      END WHERE
!
!*           REMOVES NON-PHYSICAL LOW VALUES
      GEVAP(:,:,:) = PRRS(:,:,:)<ZRTMIN(3) .AND. PCRS(:,:,:)< ZCTMIN(3)
      WHERE (GEVAP(:,:,:))
         PRVS(:,:,:) = PRVS(:,:,:) + PRRS(:,:,:)
         PTHS(:,:,:) = PTHS(:,:,:) - PRRS(:,:,:) * ZLV(:,:,:) /                         &
              ( PEXNREF(:,:,:)*(XCPD + XCPV*PRVT(:,:,:) + XCL*(PRCT(:,:,:) + PRRT(:,:,:)) ) )
         PCRS(:,:,:) = 0.0
         PRRS(:,:,:) = 0.0
      END WHERE
   ELSE
      GEVAP(:,:,:) = PRRS(:,:,:)>ZRTMIN(3) .AND. PCRS(:,:,:)>ZCTMIN(3)
      ZDR(:,:,:) = 9999.
      WHERE (GEVAP(:,:,:))
         ZDR(:,:,:)=(6.*PRRS(:,:,:)/XPI/XRHOLW/PCRS(:,:,:))**0.33
         ZWLBDR3(:,:,:) = XLBR * PCRS(:,:,:) / PRRS(:,:,:)
         ZWLBDR(:,:,:)  = ZWLBDR3(:,:,:)**XLBEXR
      END WHERE
   !
      WHERE (GEVAP(:,:,:) .AND. ZDR(:,:,:).LT.82.E-6)
         PCCS(:,:,:)   = PCCS(:,:,:)+PCRS(:,:,:)
         PCRS(:,:,:)   = 0.
         PRCS(:,:,:)   = PRCS(:,:,:)+PRRS(:,:,:)
         PRRS(:,:,:)   = 0.
      END WHERE
!!$   GMICRO(:,:,:) = GEVAP(:,:,:) .AND. ZWLBDR(:,:,:)/XACCR1>ZWLBDC3(:,:,:)
!!$                          ! the raindrops are too small, that is lower than D_h
!!$   ZFACT = 1.2E4*XACCR1
!!$   WHERE (GMICRO(:,:,:))
!!$      ZWLBDC(:,:,:) = XLBR / MIN( ZFACT,ZWLBDC3(:,:,:) )**3
!!$      ZW(:,:,:) = MIN( MAX(                                                      &
!!$           (PRHODREF(:,:,:)*PRRS(:,:,:) - ZWLBDC(:,:,:)*PCRS(:,:,:)) / &
!!$           (PRHODREF(:,:,:)*PRCS(:,:,:)/PCCS(:,:,:) - ZWLBDC(:,:,:)) , &
!!$                    0.0 ),PCRS(:,:,:),                                         &
!!$                          PCCS(:,:,:)*PRRS(:,:,:)/(PRCS(:,:,:)))
!!$!
!!$! Compute the percent (=1 if (ZWLBDR/XACCR1) >= 1.2E4
!!$! of transfer with    (=0 if (ZWLBDR/XACCR1) <= (XACCR4/ZWLBDC-XACCR5)/XACCR3
!!$!
!!$      ZW(:,:,:) = ZW(:,:,:)*( (MIN(ZWLBDR(:,:,:),1.2E4*XACCR1)-ZWLBDC3(:,:,:)) / &
!!$                            (                  1.2E4*XACCR1 -ZWLBDC3(:,:,:))   )
!!$!
!!$      ZW2(:,:,:) = PCCS(:,:,:)      !temporary storage
!!$      PCCS(:,:,:)   = PCCS(:,:,:)+ZW(:,:,:)
!!$      PCRS(:,:,:)   = PCRS(:,:,:)-ZW(:,:,:)
!!$      ZW(:,:,:)     = ZW(:,:,:) * (PRHODREF(:,:,:)*PRCS(:,:,:)/ZW2(:,:,:))
!!$      PRCS(:,:,:)   = PRCS(:,:,:)+ZW(:,:,:)
!!$      PRRS(:,:,:)   = PRRS(:,:,:)-ZW(:,:,:)
!!$   END WHERE
!!$!
!!$   GEVAP(:,:,:) = PRRS(:,:,:)<ZRTMIN(3) .OR. PCRS(:,:,:)<ZCTMIN(3)
!!$   WHERE (GEVAP(:,:,:))
!!$      PCRS(:,:,:) = 0.0
!!$      PRRS(:,:,:) = 0.0
!!$   END WHERE
!
   END IF ! LKHKO
   !
END IF ! IEVAP
!
!++cb++
DEALLOCATE(ZRTMIN)
DEALLOCATE(ZCTMIN)
!--cb--
!
!-----------------------------------------------------------------------------
!
END SUBROUTINE LIMA_WARM_EVAP
