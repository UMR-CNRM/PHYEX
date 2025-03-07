!MNH_LIC Copyright 2013-2024 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-------------------------------------------------------------------------------
MODULE MODE_LIMA_CCN_HOM_FREEZING
  IMPLICIT NONE
CONTAINS
!     ##########################################################################
  SUBROUTINE LIMA_CCN_HOM_FREEZING (LIMAP, LIMAC, D, CST, PRHODREF, PEXNREF, PPABST, PW_NU, &
                                    PTHT, PRVT, PRCT, PRRT, PRIT, PRST, PRGT, &
                                    PCCT, PCRT, PCIT, PCIT_SHAPE, PNFT, PNHT ,            &
                                    PICEFR, PTOT_RV_HONH                      )
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
!  C. Barthe   07/06/2022: save mixing ratio change for cloud electrification
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_T
USE MODD_CST,            ONLY: CST_T
!
USE MODE_TOOLS,           only: COUNTJV
USE MODD_PARAM_LIMA_COLD, ONLY:PARAM_LIMA_COLD_T
USE MODD_PARAM_LIMA, ONLY:PARAM_LIMA_T
USE YOMHOOK, ONLY:LHOOK, DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(PARAM_LIMA_COLD_T),INTENT(IN)::LIMAC
TYPE(PARAM_LIMA_T),INTENT(IN)::LIMAP
TYPE(DIMPHYEX_T),         INTENT(IN)    :: D
TYPE(CST_T),              INTENT(IN)    :: CST
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PEXNREF ! Reference Exner function
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PPABST  ! abs. pressure at time t
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PW_NU   ! updraft velocity used for
                                                   ! the nucleation param.
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PTHT    ! Theta at time t
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PRVT    ! Water vapor m.r. at t 
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PRCT    ! Cloud water m.r. at t 
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PRRT    ! Rain water m.r. at t 
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PRIT    ! Cloud ice m.r. at t 
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PRST    ! Snow/aggregate m.r. at t 
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PRGT    ! Graupel m.r. at t 
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PCCT    ! Cloud water C. at t
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PCRT    ! Rain water C. source
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PCIT    ! Ice crystal C. source
REAL, DIMENSION(D%NIJT,D%NKT,LIMAP%NNB_CRYSTAL_SHAPE), INTENT(INOUT) :: PCIT_SHAPE ! Ice crystal conc. at t for each shape
!
REAL, DIMENSION(D%NIJT,D%NKT,LIMAP%NMOD_CCN), INTENT(INOUT) :: PNFT    ! Free CCN conc. 
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PNHT    ! haze homogeneous freezing
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PICEFR  ! Ice fraction
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PTOT_RV_HONH ! Mixing ratio change due to HONH
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
REAL, DIMENSION(:,:), ALLOCATABLE :: ZCIT_SHAPE ! Ice crystal conc. at t for each shape
REAL, DIMENSION(:),   ALLOCATABLE :: ZZNHT   ! Nucleated Ice nuclei conc. source
                                             !by Homogeneous freezing
!
REAL, DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2))   &
                                  :: ZNHT  ! Nucleated Ice nuclei conc. source
                                           ! by Homogeneous freezing of haze
REAL, DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2))   &
                                  :: ZT ! work arrays
!
REAL, DIMENSION(:), ALLOCATABLE &
                           :: ZRHODREF, & ! RHO Dry REFerence
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
INTEGER :: IIJB, IIJE, IKB, IKE   ! Physical domain
INTEGER :: IL, IMOD_CCN, ISH      ! Loop index
!
INTEGER :: INEGT                          ! Case number of hom. nucleation
LOGICAL, DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2)) &
        :: GNEGT        ! Test where to compute the hom. nucleation
INTEGER , DIMENSION(SIZE(GNEGT)) :: I1,I3 ! Used to replace the COUNT
!
REAL    :: ZEPS                           ! molar mass ratio
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!
!*       1.     Preliminary computations and packing
!               ------------------------------------
!
IF (LHOOK) CALL DR_HOOK('LIMA_CCN_HOM_FREEZING', 0, ZHOOK_HANDLE)
!
! Temperature
ZT(:,:) = PTHT(:,:) * ( PPABST(:,:)/CST%XP00 ) ** (CST%XRD/CST%XCPD)
!
ZNHT(:,:) = PNHT(:,:)
!
! Computations only where the temperature is below -35°C
! PACK variables
!
GNEGT(:,:) = .FALSE.
GNEGT(D%NIJB:D%NIJE,D%NKB:D%NKE) = ZT(D%NIJB:D%NIJE,D%NKB:D%NKE)<CST%XTT-35.0
INEGT = COUNTJV( GNEGT(:,:),I1(:),I3(:))
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
   IF (LIMAP%LCRYSTAL_SHAPE) ALLOCATE(ZCIT_SHAPE(INEGT,LIMAP%NNB_CRYSTAL_SHAPE)) !++cb--
   !
   ALLOCATE(ZNFT(INEGT,LIMAP%NMOD_CCN))
   ALLOCATE(ZZNHT(INEGT))
   !
   ALLOCATE(ZRHODREF(INEGT)) 
   ALLOCATE(ZZT(INEGT)) 
   ALLOCATE(ZPRES(INEGT)) 
   ALLOCATE(ZEXNREF(INEGT))
   !
   DO IL=1,INEGT
      ZRVT(IL) = PRVT(I1(IL),I3(IL))
      ZRCT(IL) = PRCT(I1(IL),I3(IL))
      ZRRT(IL) = PRRT(I1(IL),I3(IL))
      ZRIT(IL) = PRIT(I1(IL),I3(IL))
      ZRST(IL) = PRST(I1(IL),I3(IL))
      ZRGT(IL) = PRGT(I1(IL),I3(IL))
      !
      ZTHT(IL) = PTHT(I1(IL),I3(IL))
      !
      ZCCT(IL) = PCCT(I1(IL),I3(IL))
      ZCRT(IL) = PCRT(I1(IL),I3(IL))
      ZCIT(IL) = PCIT(I1(IL),I3(IL))
      !
      IF (LIMAP%LCRYSTAL_SHAPE) THEN
         DO ISH = 1, LIMAP%NNB_CRYSTAL_SHAPE
            ZCIT_SHAPE(IL,ISH) = PCIT_SHAPE(I1(IL),I3(IL),ISH)
         END DO
      END IF
      DO IMOD_CCN = 1, LIMAP%NMOD_CCN
         ZNFT(IL,IMOD_CCN) = PNFT(I1(IL),I3(IL),IMOD_CCN)
      ENDDO
      ZZNHT(IL) = ZNHT(I1(IL),I3(IL))
      ZRHODREF(IL) = PRHODREF(I1(IL),I3(IL))
      ZZT(IL)      = ZT(I1(IL),I3(IL))
      ZPRES(IL)    = PPABST(I1(IL),I3(IL))
      ZEXNREF(IL)  = PEXNREF(I1(IL),I3(IL))
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
   ZTCELSIUS(:) = ZZT(:)-CST%XTT                                    ! T [°C]
   ZZW(:)  = ZEXNREF(:)*( CST%XCPD+CST%XCPV*ZRVT(:)+CST%XCL*(ZRCT(:)+ZRRT(:)) &
             +CST%XCI*(ZRIT(:)+ZRST(:)+ZRGT(:)) )
   ZLSFACT(:) = (CST%XLSTT+(CST%XCPV-CST%XCI)*ZTCELSIUS(:))/ZZW(:)          ! L_s/(Pi_ref*C_ph)
   ZLVFACT(:) = (CST%XLVTT+(CST%XCPV-CST%XCL)*ZTCELSIUS(:))/ZZW(:)          ! L_v/(Pi_ref*C_ph)
!
   ZZW(:)  = EXP( CST%XALPI - CST%XBETAI/ZZT(:) - CST%XGAMI*ALOG(ZZT(:) ) ) ! es_i
   ZSI(:)  = ZRVT(:)*(ZPRES(:)-ZZW(:))/((CST%XMV/CST%XMD)*ZZW(:))       ! Saturation over ice
!
!
!-------------------------------------------------------------------------------
!
!
!*       2.     Haze homogeneous freezing
!               ------------------------
!
!
!  Compute the haze homogeneous nucleation source: RHHONI
!
   IF( LIMAP%NMOD_CCN.GT.0 ) THEN

! Sum of the available CCN
      ALLOCATE( ZFREECCN(INEGT) )
      ALLOCATE( ZCCNFROZEN(INEGT) )
      ZFREECCN(:)=0.
      ZCCNFROZEN(:)=0.
      DO IMOD_CCN = 1, LIMAP%NMOD_CCN
         ZFREECCN(:) = ZFREECCN(:) + ZNFT(:,IMOD_CCN)
      END DO
!
      ALLOCATE(ZW_NU(INEGT))
      DO IL=1,INEGT
         ZW_NU(IL) = PW_NU(I1(IL),I3(IL))
      END DO
!
      ZZW(:)  = 0.0
      ZZX(:)  = 0.0
      ZEPS    = CST%XMV / CST%XMD
      ZZY(:)  = LIMAC%XCRITSAT1_HONH -                              &  ! Critical Sat.
              (MIN( LIMAC%XTMAX_HONH,MAX( LIMAC%XTMIN_HONH,ZZT(:) ) )/LIMAC%XCRITSAT2_HONH)
!
      ALLOCATE(ZLS(INEGT))
      ALLOCATE(ZPSI1(INEGT))
      ALLOCATE(ZPSI2(INEGT))
      ALLOCATE(ZTAU(INEGT))
      ALLOCATE(ZBFACT(INEGT))
!
      WHERE( (ZZT(:)<CST%XTT-35.0) .AND. (ZSI(:)>ZZY(:)) )
            ZLS(:)   = CST%XLSTT+(CST%XCPV-CST%XCI)*ZTCELSIUS(:)          ! Ls
!
            ZPSI1(:) = ZZY(:) * (CST%XG/(CST%XRD*ZZT(:)))*(ZEPS*ZLS(:)/(CST%XCPD*ZZT(:))-1.)
!                                                         ! Psi1 (a1*Scr in KL01)
! BV correction PSI2 enlever 1/ZEPS ?
!            ZPSI2(:) = ZSI(:) * (1.0/ZEPS+1.0/ZRVT(:)) +                           &
            ZPSI2(:) = ZSI(:) * (1.0/ZRVT(:)) +                           &
                 ZZY(:) * ((ZLS(:)/ZZT(:))**2)/(CST%XCPD*CST%XRV) 
!                                                         ! Psi2 (a2+a3*Scr in KL01)
            ZTAU(:) = 1.0 / ( MAX( LIMAC%XC1_HONH,LIMAC%XC1_HONH*(LIMAC%XC2_HONH-LIMAC%XC3_HONH*ZZT(:)) ) *&
                 ABS( (LIMAC%XDLNJODT1_HONH - LIMAC%XDLNJODT2_HONH*ZZT(:))       *             &
                 ((ZPRES(:)/CST%XP00)**(CST%XRD/CST%XCPD))*ZTHT(:) ) )
!
            ZBFACT(:) = (LIMAC%XRHOI_HONH/ZRHODREF(:)) * (ZSI(:)/(ZZY(:)-1.0))           &
! BV correction ZBFACT enlever 1/ZEPS ?
!                 * (1.0/ZRVT(:)+1.0/ZEPS)                                          &
                 * (1.0/ZRVT(:))                                          &
                 / (LIMAC%XCOEF_DIFVAP_HONH*(ZZT(:)**LIMAC%XCEXP_DIFVAP_HONH /ZPRES(:)))
!
! BV correction ZZX rho_i{-1} ?
!            ZZX(:) = MAX( MIN( LIMAC%XRHOI_HONH*ZBFACT(:)**1.5 * (ZPSI1(:)/ZPSI2(:))     &
            ZZX(:) = MAX( MIN( (1/LIMAC%XRHOI_HONH)*ZBFACT(:)**1.5 * (ZPSI1(:)/ZPSI2(:))     &
                 * (ZW_NU(:)/SQRT(ZTAU(:))) , ZFREECCN(:) ) , 0.)
!
            ZZW(:) = MIN( LIMAC%XRCOEF_HONH*ZZX(:)*(ZTAU(:)/ZBFACT(:))**1.5 , ZRVT(:) )
      END WHERE
!
! Apply the changes 
      DO IMOD_CCN = 1, LIMAP%NMOD_CCN
         WHERE(ZFREECCN(:)>1.)
            ZCCNFROZEN(:) = ZZX(:) * ZNFT(:,IMOD_CCN)/ZFREECCN(:)
         END WHERE
         PNFT(:,:,IMOD_CCN) = PNFT(:,:,IMOD_CCN) - UNPACK( ZCCNFROZEN(:), MASK=GNEGT(:,:),FIELD=0.)
      END DO
!
      PTOT_RV_HONH(:,:) = UNPACK( ZZW(:), MASK=GNEGT(:,:),FIELD=0.)
!
      PTHT(:,:) = PTHT(:,:) + UNPACK( ZZW(:)*(ZLSFACT(:)-ZLVFACT(:)), MASK=GNEGT(:,:),FIELD=0.)
      PRVT(:,:) = PRVT(:,:) - UNPACK( ZZW(:), MASK=GNEGT(:,:),FIELD=0.)
      PRIT(:,:) = PRIT(:,:) + UNPACK( ZZW(:), MASK=GNEGT(:,:),FIELD=0.)
      IF (.NOT. LIMAP%LCRYSTAL_SHAPE) THEN
        PCIT(:,:) = PCIT(:,:) + UNPACK( ZZX(:), MASK=GNEGT(:,:),FIELD=0.)
      ELSE
        !--> Hyp : etant donnee la gamme de temperatures, formation de 20% de colonnes et
        ! 80% de polycristaux
        !PCIT_SHAPE(:,:,:,2) = PCIT_SHAPE(:,:,:,2) + &
        !                      0.2 * UNPACK( ZZX(:), MASK=GNEGT(:,:,:),FIELD=0.)
        !PCIT_SHAPE(:,:,:,3) = PCIT_SHAPE(:,:,:,3) + &
        !                      0.8 * UNPACK( ZZX(:), MASK=GNEGT(:,:,:),FIELD=0.)
        !--> Hyp : la congelation des CCN et des gouttelettes produit des droxtals
        PCIT_SHAPE(:,:,4) = PCIT_SHAPE(:,:,4) + UNPACK( ZZX(:), MASK=GNEGT(:,:),FIELD=0.)
        PCIT(:,:) = SUM(PCIT_SHAPE, DIM=3)
      END IF
      PNHT(:,:) = PNHT(:,:) + UNPACK( ZZX(:), MASK=GNEGT(:,:),FIELD=0.)

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
   IF (ALLOCATED(ZCIT_SHAPE)) DEALLOCATE(ZCIT_SHAPE) !++cb-- 19/02/24
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
IF (LHOOK) CALL DR_HOOK('LIMA_CCN_HOM_FREEZING', 1, ZHOOK_HANDLE)
END SUBROUTINE LIMA_CCN_HOM_FREEZING
END MODULE MODE_LIMA_CCN_HOM_FREEZING
