!MNH_LIC Copyright 2013-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-------------------------------------------------------------------------------
!      #################################
       MODULE MODI_LIMA_CCN_HOM_FREEZING
!      #################################
!
INTERFACE
   SUBROUTINE LIMA_CCN_HOM_FREEZING (PRHODREF, PEXNREF, PPABST, PW_NU,         &
                                     PTHT, PRVT, PRCT, PRRT, PRIT, PRST, PRGT, &
                                     PCCT, PCRT, PCIT, PNFT, PNHT,             &
                                     PICEFR                                    )
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PEXNREF ! Reference Exner function
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABST  ! abs. pressure at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PW_NU   ! updraft velocity used for
                                                   ! the nucleation param.
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PTHT    ! Theta at time t
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRVT    ! Water vapor m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRCT    ! Cloud water m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRRT    ! Rain water m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRIT    ! Cloud ice m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRST    ! Snow/aggregate m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRGT    ! Graupel m.r. at t 
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCCT    ! Cloud water C. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCRT    ! Rain water C. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCIT    ! Ice crystal C. source
!
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PNFT    ! Free CCN conc. 
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PNHT    ! haze homogeneous freezing
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PICEFR  ! Ice fraction
!
END SUBROUTINE LIMA_CCN_HOM_FREEZING
END INTERFACE
END MODULE MODI_LIMA_CCN_HOM_FREEZING
!
!     ##########################################################################
   SUBROUTINE LIMA_CCN_HOM_FREEZING (PRHODREF, PEXNREF, PPABST, PW_NU,         &
                                     PTHT, PRVT, PRCT, PRRT, PRIT, PRST, PRGT, &
                                     PCCT, PCRT, PCIT, PNFT, PNHT ,            &
                                     PICEFR                                    )
!     ##########################################################################
!
!!    PURPOSE
!!    -------
!!      Compute the homogeneous freezing of CCN where T<-35°C
!!
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
!!      Original             15/03/2018 
!  P. Wautelet 28/05/2019: move COUNTJV function to tools.f90
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST,             ONLY: XP00, XRD, XRV, XMV, XMD, XCPD, XCPV, XCL, XCI,   &
                                XTT, XLSTT, XLVTT, XALPI, XBETAI, XGAMI,          &
                                XG
USE MODD_NSV
USE MODD_PARAMETERS,      ONLY: JPHEXT, JPVEXT
USE MODD_PARAM_LIMA,      ONLY: NMOD_CCN, NMOD_IMM, XRTMIN, XCTMIN, XNUC
USE MODD_PARAM_LIMA_COLD, ONLY: XRCOEF_HONH, XCEXP_DIFVAP_HONH, XCOEF_DIFVAP_HONH,&
                                XCRITSAT1_HONH, XCRITSAT2_HONH, XTMAX_HONH,       &
                                XTMIN_HONH, XC1_HONH, XC2_HONH, XC3_HONH,         &
                                XDLNJODT1_HONH, XDLNJODT2_HONH, XRHOI_HONH,       &
                                XC_HONC, XTEXP1_HONC, XTEXP2_HONC, XTEXP3_HONC,   &
                                XTEXP4_HONC, XTEXP5_HONC
USE MODD_PARAM_LIMA_WARM, ONLY: XLBC
!
use mode_tools,           only: Countjv
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PEXNREF ! Reference Exner function
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABST  ! abs. pressure at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PW_NU   ! updraft velocity used for
                                                   ! the nucleation param.
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PTHT    ! Theta at time t
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRVT    ! Water vapor m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRCT    ! Cloud water m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRRT    ! Rain water m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRIT    ! Cloud ice m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRST    ! Snow/aggregate m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRGT    ! Graupel m.r. at t 
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCCT    ! Cloud water C. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCRT    ! Rain water C. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCIT    ! Ice crystal C. source
!
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PNFT    ! Free CCN conc. 
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PNHT    ! haze homogeneous freezing
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PICEFR  ! Ice fraction
!
!*       0.2   Declarations of local variables :
!
REAL, DIMENSION(:), ALLOCATABLE :: ZRVT    ! Water vapor m.r. at t
REAL, DIMENSION(:), ALLOCATABLE :: ZRCT    ! Cloud water m.r. at t
REAL, DIMENSION(:), ALLOCATABLE :: ZRRT    ! Rain water m.r. at t
REAL, DIMENSION(:), ALLOCATABLE :: ZRIT    ! Pristine ice m.r. at t
REAL, DIMENSION(:), ALLOCATABLE :: ZRST    ! Snow/aggregate m.r. at t
REAL, DIMENSION(:), ALLOCATABLE :: ZRGT    ! Graupel/hail m.r. at t
!
REAL, DIMENSION(:), ALLOCATABLE :: ZTHT    ! Theta source
!
REAL, DIMENSION(:),   ALLOCATABLE :: ZCCT    ! Cloud water conc. source
REAL, DIMENSION(:),   ALLOCATABLE :: ZCRT    ! Rain water conc. source
REAL, DIMENSION(:,:), ALLOCATABLE :: ZNFT    ! available nucleus conc. source
REAL, DIMENSION(:),   ALLOCATABLE :: ZCIT    ! Pristine ice conc. source
REAL, DIMENSION(:),   ALLOCATABLE :: ZZNHT   ! Nucleated Ice nuclei conc. source
                                             !by Homogeneous freezing
!
REAL, DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3))   &
                                  :: ZNHT  ! Nucleated Ice nuclei conc. source
                                           ! by Homogeneous freezing of haze
REAL, DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3))   &
                                  :: ZW, ZT ! work arrays
!
REAL, DIMENSION(:), ALLOCATABLE &
                           :: ZRHODREF, & ! RHO Dry REFerence
                              ZRHODJ,   & ! RHO times Jacobian
                              ZZT,      & ! Temperature
                              ZPRES,    & ! Pressure
                              ZEXNREF,  & ! EXNer Pressure REFerence
                              ZZW,      & ! Work array
                              ZZX,      & ! Work array
                              ZZY,      & ! Work array
                              ZLSFACT,  & ! L_s/(Pi_ref*C_ph)
                              ZLVFACT,  & ! L_v/(Pi_ref*C_ph)
                              ZLBDAC,   & ! Slope parameter of the cloud droplet distr.
                              ZSI,      & ! Saturation over ice
                              ZTCELSIUS,&
                              ZLS,      &
                              ZPSI1,    &
                              ZPSI2,    &
                              ZTAU,     &
                              ZBFACT,   &
                              ZW_NU,    &
                              ZFREECCN, &
                              ZCCNFROZEN
!
INTEGER :: IIB, IIE, IJB, IJE, IKB, IKE   ! Physical domain
INTEGER :: JL, JMOD_CCN, JMOD_IMM         ! Loop index
!
INTEGER :: INEGT                          ! Case number of hom. nucleation
LOGICAL, DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3)) &
			  :: GNEGT        ! Test where to compute the hom. nucleation
INTEGER , DIMENSION(SIZE(GNEGT)) :: I1,I2,I3 ! Used to replace the COUNT
!
REAL    :: ZEPS                           ! molar mass ratio
!
!-------------------------------------------------------------------------------
!
!
!*       1.     Preliminary computations and packing
!	        ------------------------------------
!
!
! Physical domain
IIB=1+JPHEXT
IIE=SIZE(PTHT,1) - JPHEXT
IJB=1+JPHEXT
IJE=SIZE(PTHT,2) - JPHEXT
IKB=1+JPVEXT
IKE=SIZE(PTHT,3) - JPVEXT
!
! Temperature
ZT(:,:,:) = PTHT(:,:,:) * ( PPABST(:,:,:)/XP00 ) ** (XRD/XCPD)
!
ZNHT(:,:,:) = PNHT(:,:,:)
!
! Computations only where the temperature is below -35°C
! PACK variables
!
GNEGT(:,:,:) = .FALSE.
GNEGT(IIB:IIE,IJB:IJE,IKB:IKE) = ZT(IIB:IIE,IJB:IJE,IKB:IKE)<XTT-35.0
INEGT = COUNTJV( GNEGT(:,:,:),I1(:),I2(:),I3(:))
!
IF (INEGT.GT.0) THEN

   ALLOCATE(ZRVT(INEGT)) 
   ALLOCATE(ZRCT(INEGT)) 
   ALLOCATE(ZRRT(INEGT)) 
   ALLOCATE(ZRIT(INEGT)) 
   ALLOCATE(ZRST(INEGT)) 
   ALLOCATE(ZRGT(INEGT)) 
   !
   ALLOCATE(ZTHT(INEGT))
   !
   ALLOCATE(ZCCT(INEGT))
   ALLOCATE(ZCRT(INEGT))
   ALLOCATE(ZCIT(INEGT))
   !
   ALLOCATE(ZNFT(INEGT,NMOD_CCN))
   ALLOCATE(ZZNHT(INEGT))
   !
   ALLOCATE(ZRHODREF(INEGT)) 
   ALLOCATE(ZZT(INEGT)) 
   ALLOCATE(ZPRES(INEGT)) 
   ALLOCATE(ZEXNREF(INEGT))
   !
   DO JL=1,INEGT
      ZRVT(JL) = PRVT(I1(JL),I2(JL),I3(JL))
      ZRCT(JL) = PRCT(I1(JL),I2(JL),I3(JL))
      ZRRT(JL) = PRRT(I1(JL),I2(JL),I3(JL))
      ZRIT(JL) = PRIT(I1(JL),I2(JL),I3(JL))
      ZRST(JL) = PRST(I1(JL),I2(JL),I3(JL))
      ZRGT(JL) = PRGT(I1(JL),I2(JL),I3(JL))
      !
      ZTHT(JL) = PTHT(I1(JL),I2(JL),I3(JL))
      !
      ZCCT(JL) = PCCT(I1(JL),I2(JL),I3(JL))
      ZCRT(JL) = PCRT(I1(JL),I2(JL),I3(JL))
      ZCIT(JL) = PCIT(I1(JL),I2(JL),I3(JL))
      !
      DO JMOD_CCN = 1, NMOD_CCN
         ZNFT(JL,JMOD_CCN) = PNFT(I1(JL),I2(JL),I3(JL),JMOD_CCN)
      ENDDO
      ZZNHT(JL) = ZNHT(I1(JL),I2(JL),I3(JL))
      ZRHODREF(JL) = PRHODREF(I1(JL),I2(JL),I3(JL))
      ZZT(JL)      = ZT(I1(JL),I2(JL),I3(JL))
      ZPRES(JL)    = PPABST(I1(JL),I2(JL),I3(JL))
      ZEXNREF(JL)  = PEXNREF(I1(JL),I2(JL),I3(JL))
   ENDDO
!
! PACK : done
! Prepare computations
!
   ALLOCATE( ZLSFACT    (INEGT) )
   ALLOCATE( ZLVFACT    (INEGT) )
   ALLOCATE( ZSI        (INEGT) )
   ALLOCATE( ZTCELSIUS  (INEGT) )
   ALLOCATE( ZLBDAC     (INEGT) )
!
   ALLOCATE( ZZW (INEGT) ) ; ZZW(:) = 0.0
   ALLOCATE( ZZX (INEGT) ) ; ZZX(:) = 0.0
   ALLOCATE( ZZY (INEGT) ) ; ZZY(:) = 0.0
!
   ZTCELSIUS(:) = ZZT(:)-XTT                                    ! T [°C]
   ZZW(:)  = ZEXNREF(:)*( XCPD+XCPV*ZRVT(:)+XCL*(ZRCT(:)+ZRRT(:)) &
             +XCI*(ZRIT(:)+ZRST(:)+ZRGT(:)) )
   ZLSFACT(:) = (XLSTT+(XCPV-XCI)*ZTCELSIUS(:))/ZZW(:)          ! L_s/(Pi_ref*C_ph)
   ZLVFACT(:) = (XLVTT+(XCPV-XCL)*ZTCELSIUS(:))/ZZW(:)          ! L_v/(Pi_ref*C_ph)
!
   ZZW(:)  = EXP( XALPI - XBETAI/ZZT(:) - XGAMI*ALOG(ZZT(:) ) ) ! es_i
   ZSI(:)  = ZRVT(:)*(ZPRES(:)-ZZW(:))/((XMV/XMD)*ZZW(:))       ! Saturation over ice
!
!
!-------------------------------------------------------------------------------
!
!
!*       2.     Haze homogeneous freezing
!	        ------------------------
!
!
!  Compute the haze homogeneous nucleation source: RHHONI
!
   IF( NMOD_CCN.GT.0 ) THEN

! Sum of the available CCN
      ALLOCATE( ZFREECCN(INEGT) )
      ALLOCATE( ZCCNFROZEN(INEGT) )
      ZFREECCN(:)=0.
      ZCCNFROZEN(:)=0.
      DO JMOD_CCN = 1, NMOD_CCN
         ZFREECCN(:) = ZFREECCN(:) + ZNFT(:,JMOD_CCN)
      END DO
!
      ALLOCATE(ZW_NU(INEGT))
      DO JL=1,INEGT
         ZW_NU(JL) = PW_NU(I1(JL),I2(JL),I3(JL))
      END DO
!
      ZZW(:)  = 0.0
      ZZX(:)  = 0.0
      ZEPS    = XMV / XMD
      ZZY(:)  = XCRITSAT1_HONH -                              &  ! Critical Sat.
              (MIN( XTMAX_HONH,MAX( XTMIN_HONH,ZZT(:) ) )/XCRITSAT2_HONH)
!
      ALLOCATE(ZLS(INEGT))
      ALLOCATE(ZPSI1(INEGT))
      ALLOCATE(ZPSI2(INEGT))
      ALLOCATE(ZTAU(INEGT))
      ALLOCATE(ZBFACT(INEGT))
!
      WHERE( (ZZT(:)<XTT-35.0) .AND. (ZSI(:)>ZZY(:)) )
            ZLS(:)   = XLSTT+(XCPV-XCI)*ZTCELSIUS(:)          ! Ls
!
            ZPSI1(:) = ZZY(:) * (XG/(XRD*ZZT(:)))*(ZEPS*ZLS(:)/(XCPD*ZZT(:))-1.)
!                                                         ! Psi1 (a1*Scr in KL01)
! BV correction PSI2 enlever 1/ZEPS ?
!            ZPSI2(:) = ZSI(:) * (1.0/ZEPS+1.0/ZRVT(:)) +                           &
            ZPSI2(:) = ZSI(:) * (1.0/ZRVT(:)) +                           &
                 ZZY(:) * ((ZLS(:)/ZZT(:))**2)/(XCPD*XRV) 
!                                                         ! Psi2 (a2+a3*Scr in KL01)
            ZTAU(:) = 1.0 / ( MAX( XC1_HONH,XC1_HONH*(XC2_HONH-XC3_HONH*ZZT(:)) ) *&
                 ABS( (XDLNJODT1_HONH - XDLNJODT2_HONH*ZZT(:))       *             &
                 ((ZPRES(:)/XP00)**(XRD/XCPD))*ZTHT(:) ) )
!
            ZBFACT(:) = (XRHOI_HONH/ZRHODREF(:)) * (ZSI(:)/(ZZY(:)-1.0))           &
! BV correction ZBFACT enlever 1/ZEPS ?
!                 * (1.0/ZRVT(:)+1.0/ZEPS)                                          &
                 * (1.0/ZRVT(:))                                          &
                 / (XCOEF_DIFVAP_HONH*(ZZT(:)**XCEXP_DIFVAP_HONH /ZPRES(:)))
!
! BV correction ZZX rho_i{-1} ?
!            ZZX(:) = MAX( MIN( XRHOI_HONH*ZBFACT(:)**1.5 * (ZPSI1(:)/ZPSI2(:))     &
            ZZX(:) = MAX( MIN( (1/XRHOI_HONH)*ZBFACT(:)**1.5 * (ZPSI1(:)/ZPSI2(:))     &
                 * (ZW_NU(:)/SQRT(ZTAU(:))) , ZFREECCN(:) ) , 0.)
!
            ZZW(:) = MIN( XRCOEF_HONH*ZZX(:)*(ZTAU(:)/ZBFACT(:))**1.5 , ZRVT(:) )
      END WHERE
!
! Apply the changes 
      DO JMOD_CCN = 1, NMOD_CCN
         WHERE(ZFREECCN(:)>1.)
            ZCCNFROZEN(:) = ZZX(:) * ZNFT(:,JMOD_CCN)/ZFREECCN(:)
         END WHERE
         PNFT(:,:,:,JMOD_CCN) = PNFT(:,:,:,JMOD_CCN) - UNPACK( ZCCNFROZEN(:), MASK=GNEGT(:,:,:),FIELD=0.)
      END DO
!
      PTHT(:,:,:) = PTHT(:,:,:) + UNPACK( ZZW(:)*(ZLSFACT(:)-ZLVFACT(:)), MASK=GNEGT(:,:,:),FIELD=0.)
      PRVT(:,:,:) = PRVT(:,:,:) - UNPACK( ZZW(:), MASK=GNEGT(:,:,:),FIELD=0.)
      PRIT(:,:,:) = PRIT(:,:,:) + UNPACK( ZZW(:), MASK=GNEGT(:,:,:),FIELD=0.)
      PCIT(:,:,:) = PCIT(:,:,:) + UNPACK( ZZX(:), MASK=GNEGT(:,:,:),FIELD=0.)
      PNHT(:,:,:) = PNHT(:,:,:) + UNPACK( ZZX(:), MASK=GNEGT(:,:,:),FIELD=0.)

      DEALLOCATE(ZFREECCN)
      DEALLOCATE(ZCCNFROZEN)
      DEALLOCATE(ZLS)
      DEALLOCATE(ZPSI1)
      DEALLOCATE(ZPSI2)
      DEALLOCATE(ZTAU)
      DEALLOCATE(ZBFACT)
      DEALLOCATE(ZW_NU)
!
   END IF
!
!
   DEALLOCATE(ZRVT) 
   DEALLOCATE(ZRCT) 
   DEALLOCATE(ZRRT) 
   DEALLOCATE(ZRIT) 
   DEALLOCATE(ZRST) 
   DEALLOCATE(ZRGT) 
!
   DEALLOCATE(ZTHT)
!
   DEALLOCATE(ZCCT)
   DEALLOCATE(ZCRT)
   DEALLOCATE(ZCIT)
!
   DEALLOCATE(ZNFT)
   DEALLOCATE(ZZNHT)
!
   DEALLOCATE(ZRHODREF) 
   DEALLOCATE(ZZT) 
   DEALLOCATE(ZPRES) 
   DEALLOCATE(ZEXNREF)
!
   DEALLOCATE(ZLSFACT)
   DEALLOCATE(ZLVFACT)
   DEALLOCATE(ZSI)
   DEALLOCATE(ZTCELSIUS)
   DEALLOCATE(ZLBDAC)
!
   DEALLOCATE(ZZW) 
   DEALLOCATE(ZZX)
   DEALLOCATE(ZZY)
!
!
END IF ! INEGT>0
!
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_CCN_HOM_FREEZING
