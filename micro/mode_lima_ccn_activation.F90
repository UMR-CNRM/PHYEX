!MNH_LIC Copyright 2013-2024 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_LIMA_CCN_ACTIVATION
  IMPLICIT NONE
CONTAINS
!     ##############################################################################
    SUBROUTINE LIMA_CCN_ACTIVATION (LIMAP, LIMAW, TNSV, D, CST, NEBN,              &
                                    KCARB, KSOA, KSP, ODUST, OSALT, OORILAM,       &
                                    PRHODREF, PEXNREF, PPABST, PT, PDTHRAD, PW_NU, &
                                    PAERO,PSOLORG, PMI,  HACTCCN,                  &
                                    PTHT, PRVT, PRCT, PCCT, PRRT, PNFT, PNAT,      &
                                    PCLDFR, PTOT_RV_HENU                           )
!     ##############################################################################
!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the activation of CCN 
!!    according to Cohard and Pinty, QJRMS, 2000
!!
!!
!!**  METHOD
!!    ------
!!      The activation of CCN is checked for quasi-saturated air parcels 
!!    to update the cloud droplet number concentration.
!!
!!    Computation steps :
!!      1- Check where computations are necessary
!!      2- and 3- Compute the maximum of supersaturation using the iterative 
!!                Ridder algorithm
!!      4- Compute the nucleation source
!!      5- Deallocate local variables
!! 
!!    Contains :
!!      6- Functions : Ridder algorithm
!!
!!
!!    REFERENCE
!!    ---------
!!
!!      Cohard, J.-M. and J.-P. Pinty, 2000: A comprehensive two-moment warm 
!!      microphysical bulk scheme. 
!!        Part I: Description and tests
!!        Part II: 2D experiments with a non-hydrostatic model
!!      Accepted for publication in Quart. J. Roy. Meteor. Soc. 
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
!  B. Vie      03/03/2020: use DTHRAD instead of dT/dt in Smax diagnostic computation
!  P. Wautelet 10/04/2019: replace ABORT and STOP calls by Print_msg
!  P. Wautelet 26/04/2019: replace non-standard FLOAT function by REAL function
!  P. Wautelet 28/05/2019: move COUNTJV function to tools.f90
!  C. Barthe      06/2022: save mixing ratio change for cloud electrification
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_T
USE MODD_CST,            ONLY: CST_T
!use modd_field,           only: TFIELDDATA, TYPEREAL
!USE MODD_IO,              ONLY: TFILEDATA
!USE MODD_LUNIT_n,         ONLY: TLUOUT
USE MODD_NEB_N,           ONLY: NEB_T
USE MODD_NEB_n,           ONLY: LSUBG_COND
USE MODD_NSV,        ONLY : NSV_T
USE MODI_CH_AER_ACTIVATION

!USE MODE_IO_FIELD_WRITE,  only: IO_Field_write
USE MODE_TOOLS,           only: COUNTJV

USE MODD_PARAM_LIMA_WARM, ONLY:PARAM_LIMA_WARM_T
USE MODD_PARAM_LIMA, ONLY:PARAM_LIMA_T
USE YOMHOOK, ONLY:LHOOK, DR_HOOK, JPHOOK

IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(PARAM_LIMA_WARM_T),INTENT(IN)::LIMAW
TYPE(PARAM_LIMA_T),INTENT(IN)::LIMAP
TYPE(NSV_T),              INTENT(IN)    :: TNSV
TYPE(DIMPHYEX_T),         INTENT(IN)    :: D
TYPE(CST_T),              INTENT(IN)    :: CST
TYPE(NEB_T),              INTENT(IN)    :: NEBN
!TYPE(TFILEDATA),          INTENT(IN)    :: TPFILE     ! Output file
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PRHODREF   ! Reference density
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PEXNREF    ! Reference Exner function
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PPABST     ! abs. pressure at time t
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PT         ! Temperature
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PDTHRAD    ! Radiative temperature tendency
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PW_NU      ! updraft velocity used for
                                                      ! the nucleation param.
REAL, DIMENSION(D%NIJT, D%NKT ,TNSV%NSV), INTENT(INOUT) :: PAERO   ! Aerosol concentration
REAL, DIMENSION(D%NIJT, D%NKT, 10), INTENT(IN)    :: PSOLORG ![%] solubility fraction of soa
REAL, DIMENSION(D%NIJT, D%NKT, KSP+KCARB+KSOA), INTENT(IN)    :: PMI
CHARACTER(LEN=4),         INTENT(IN)    :: HACTCCN  ! kind of CCN activation
INTEGER,                  INTENT(IN)    :: KCARB, KSOA, KSP ! for array size declarations
LOGICAL,                  INTENT(IN)    :: ODUST, OSALT, OORILAM

!   
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PTHT       ! Theta at t 
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PRVT       ! Water vapor m.r. at t 
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PRCT       ! Cloud water m.r. at t 
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PCCT       ! Cloud water m.r. at t 
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PRRT       ! Cloud water m.r. at t 
REAL, DIMENSION(D%NIJT,D%NKT,LIMAP%NMOD_CCN), INTENT(INOUT) :: PNFT       ! CCN C. available at t
REAL, DIMENSION(D%NIJT,D%NKT,LIMAP%NMOD_CCN), INTENT(INOUT) :: PNAT       ! CCN C. activated at t
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PCLDFR     ! Precipitation fraction
!
REAL, DIMENSION(D%NIJT,D%NKT), OPTIONAL, INTENT(INOUT) :: PTOT_RV_HENU  ! Mixing ratio change due to HENU
!
!*       0.1   Declarations of local variables :
!
! Packing variables
LOGICAL, DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2)) :: GNUCT 
INTEGER :: INUCT
INTEGER , DIMENSION(SIZE(GNUCT))   :: I1,I3 ! Used to replace the COUNT
INTEGER                            :: IL       ! and PACK intrinsics 
!
! Packed micophysical variables
REAL, DIMENSION(:)  , ALLOCATABLE  :: ZRCT     ! cloud mr
REAL, DIMENSION(:)  , ALLOCATABLE  :: ZCCT     ! cloud conc.
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZNFT     ! available nucleus conc.
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZNAT     ! activated nucleus conc.
!
! Other packed variables
REAL, DIMENSION(:)  , ALLOCATABLE  :: ZRHODREF ! RHO Dry REFerence
REAL, DIMENSION(:)  , ALLOCATABLE  :: ZEXNREF  ! EXNer Pressure REFerence
REAL, DIMENSION(:)  , ALLOCATABLE  :: ZZT      ! Temperature
!
! Work arrays
REAL, DIMENSION(:), ALLOCATABLE    :: ZZW1, ZZW2, ZZW3, ZZW4, ZZW5, ZZW6, &
                                      ZZTDT,          & ! dT/dt
                                      ZSW,            & ! real supersaturation                                      
                                      ZSMAX,          & ! Maximum supersaturation
                                      ZVEC1
!
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZTMP, ZCHEN_MULTI
!
REAL, DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2)) :: ZTDT, ZRVSAT, ZW, ZW2, ZCLDFR  
REAL, DIMENSION(SIZE(PNFT,1),SIZE(PNFT,2)) :: ZCONC_TOT         ! total CCN C. available
!
INTEGER, DIMENSION(:), ALLOCATABLE :: IVEC1             ! Vectors of indices for
                                                        ! interpolations
REAL, DIMENSION(:)  , ALLOCATABLE  :: ZPABST, ZMCN, ZNATOLD
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZAERO, ZSOLORG, ZMI
!
! 
REAL    :: ZEPS                                ! molar mass ratio
REAL    :: ZS1, ZS2, ZXACC 
INTEGER :: IMOD
INTEGER :: IIJB, IIJE, IKB, IKE        ! Physical domain
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
!!$INTEGER                  :: ILUOUT     ! Logical unit of output listing 
!!$TYPE(TFIELDMETADATA) :: TZFIELD
!-------------------------------------------------------------------------------
!
!ILUOUT = TLUOUT%NLU
!
!*       1.     PREPARE COMPUTATIONS - PACK
!               ---------------------------
!
IF (LHOOK) CALL DR_HOOK('LIMA_CCN_ACTIVATION', 0, ZHOOK_HANDLE)
!
!  Saturation vapor mixing ratio and radiative tendency                    
!
ZEPS= CST%XMV / CST%XMD
ZRVSAT(:,:) = ZEPS / (PPABST(:,:)*EXP(-CST%XALPW+CST%XBETAW/PT(:,:)+CST%XGAMW*ALOG(PT(:,:))) - 1.0)
ZTDT(:,:)   = 0.
IF (LIMAP%LACTIT .AND. SIZE(PDTHRAD).GT.0) ZTDT(:,:)   = PDTHRAD(:,:) * PEXNREF(:,:)
!
!  find locations where CCN are available
!
ZCONC_TOT(:,:) = 0.0
DO IMOD = 1, LIMAP%NMOD_CCN 
   ZCONC_TOT(:,:) = ZCONC_TOT(:,:) + PNFT(:,:,IMOD) ! sum over the free CCN
ENDDO
!
!  optimization by looking for locations where
!  the updraft velocity is positive!!!
!
GNUCT(:,:) = .FALSE.
!
IF (LIMAP%LADJ) THEN
   GNUCT(D%NIJB:D%NIJE,D%NKTB:D%NKTE) =      PW_NU(D%NIJB:D%NIJE,D%NKTB:D%NKTE)>LIMAW%XWMIN                          &
                                    .OR. PRVT(D%NIJB:D%NIJE,D%NKTB:D%NKTE)>ZRVSAT(D%NIJB:D%NIJE,D%NKTB:D%NKTE)
   IF (LIMAP%LACTIT) GNUCT(D%NIJB:D%NIJE,D%NKTB:D%NKTE) =      GNUCT(D%NIJB:D%NIJE,D%NKTB:D%NKTE)      &
                                                .OR. ZTDT(D%NIJB:D%NIJE,D%NKTB:D%NKTE)<LIMAW%XTMIN
!
   GNUCT(D%NIJB:D%NIJE,D%NKTB:D%NKTE) =       GNUCT(D%NIJB:D%NIJE,D%NKTB:D%NKTE)                       &
                                    .AND. PT(D%NIJB:D%NIJE,D%NKTB:D%NKTE)>(CST%XTT-22.)                &
                                    .AND. ZCONC_TOT(D%NIJB:D%NIJE,D%NKTB:D%NKTE)>LIMAP%XCTMIN(2)
!
   IF (NEBN%LSUBG_COND) GNUCT(D%NIJB:D%NIJE,D%NKTB:D%NKTE) = GNUCT(D%NIJB:D%NIJE,D%NKTB:D%NKTE)       &
                                              .AND. PCLDFR(D%NIJB:D%NIJE,D%NKTB:D%NKTE)>0.01
   IF (.NOT. NEBN%LSUBG_COND) GNUCT(D%NIJB:D%NIJE,D%NKTB:D%NKTE) = GNUCT(D%NIJB:D%NIJE,D%NKTB:D%NKTE)       &
                                                    .AND. PRVT(D%NIJB:D%NIJE,D%NKTB:D%NKTE).GE.ZRVSAT(D%NIJB:D%NIJE,D%NKTB:D%NKTE)
ELSE
   GNUCT(D%NIJB:D%NIJE,D%NKTB:D%NKTE) =       PRVT(D%NIJB:D%NIJE,D%NKTB:D%NKTE).GE.ZRVSAT(D%NIJB:D%NIJE,D%NKTB:D%NKTE) &
                                    .AND. PT(D%NIJB:D%NIJE,D%NKTB:D%NKTE)>(CST%XTT-22.)                            &
                                    .AND. ZCONC_TOT(D%NIJB:D%NIJE,D%NKTB:D%NKTE)>LIMAP%XCTMIN(2)
END IF
!
IF (.NOT. NEBN%LSUBG_COND) THEN
   ZCLDFR(:,:) = 1.
ELSE
   ZCLDFR(:,:) = PCLDFR(:,:)
END IF
!
INUCT = COUNTJV( GNUCT(:,:),I1(:),I3(:))
!
IF( INUCT >= 1 ) THEN
!
   ALLOCATE(ZNFT(INUCT,LIMAP%NMOD_CCN))
   ALLOCATE(ZNAT(INUCT,LIMAP%NMOD_CCN))
   ALLOCATE(ZTMP(INUCT,LIMAP%NMOD_CCN))
   ALLOCATE(ZRCT(INUCT))
   ALLOCATE(ZCCT(INUCT))
   ALLOCATE(ZZT(INUCT)) 
   ALLOCATE(ZZTDT(INUCT)) 
   ALLOCATE(ZSW(INUCT))    
   ALLOCATE(ZZW1(INUCT))
   ALLOCATE(ZZW2(INUCT))
   ALLOCATE(ZZW3(INUCT))
   ALLOCATE(ZZW4(INUCT))
   ALLOCATE(ZZW5(INUCT))
   ALLOCATE(ZZW6(INUCT))
   ALLOCATE(ZCHEN_MULTI(INUCT,LIMAP%NMOD_CCN))
   ALLOCATE(ZVEC1(INUCT))
   ALLOCATE(IVEC1(INUCT))
   ALLOCATE(ZRHODREF(INUCT)) 
   ALLOCATE(ZEXNREF(INUCT)) 
   ALLOCATE(ZPABST(INUCT))
   ALLOCATE(ZAERO(INUCT,SIZE(PAERO,3)))
   ALLOCATE(ZSOLORG(INUCT,SIZE(PSOLORG,3)))
   ALLOCATE(ZMI(INUCT,SIZE(PMI,3)))
   DO IL=1,INUCT
      ZRCT(IL) = PRCT(I1(IL),I3(IL))/ZCLDFR(I1(IL),I3(IL))
      ZCCT(IL) = PCCT(I1(IL),I3(IL))/ZCLDFR(I1(IL),I3(IL))
      ZZT(IL)  = PT(I1(IL),I3(IL))
      ZZW1(IL) = ZRVSAT(I1(IL),I3(IL))
      ZZW2(IL) = PW_NU(I1(IL),I3(IL))
      ZZTDT(IL)  = ZTDT(I1(IL),I3(IL))
      ZSW(IL)  = PRVT(I1(IL),I3(IL))/ZRVSAT(I1(IL),I3(IL)) - 1.
      ZRHODREF(IL) = PRHODREF(I1(IL),I3(IL))
      ZEXNREF(IL)  = PEXNREF(I1(IL),I3(IL))
      ZPABST(IL)  = PPABST(I1(IL),I3(IL))
      IF ((OORILAM).OR.(ODUST).OR.(OSALT)) THEN
         ZAERO(IL,:)  = PAERO(I1(IL),I3(IL),:)
      ELSE
         ZAERO(IL,:) = 0.
      END IF

      IF (OORILAM) THEN
         ZSOLORG(IL,:) = PSOLORG(I1(IL),I3(IL),:)
         ZMI(IL,:) = PMI(I1(IL),I3(IL),:)
      ELSE
         ZSOLORG(IL,:) = 0.
         ZMI(IL,:) = 0.
      END IF
      
      DO IMOD = 1,LIMAP%NMOD_CCN
         ZNFT(IL,IMOD)        = PNFT(I1(IL),I3(IL),IMOD)
         ZNAT(IL,IMOD)        = PNAT(I1(IL),I3(IL),IMOD)
         ZCHEN_MULTI(IL,IMOD) = (ZNFT(IL,IMOD)+ZNAT(IL,IMOD))*ZRHODREF(IL) &
                                                             / LIMAP%XLIMIT_FACTOR(IMOD)
      ENDDO
   ENDDO
!
IF ((HACTCCN == 'ABRK').AND.((OORILAM).OR.(ODUST).OR.(OSALT))) THEN  ! CCN activation from Abdul-Razack (only if prognostic aerosols)
!
     ALLOCATE(ZMCN(INUCT))
     ALLOCATE(ZNATOLD(INUCT))
     ALLOCATE(ZSMAX(INUCT))

     ZZW1(:) = 0.
     ZZW2(:) = 0.
     ZZW3(:) = 0.
     ZSMAX(:) = 0.

     ZMCN(:) = 0. !masse activée (non utilisée!!)

     ZNATOLD(:) = ZNAT(:,1)
   !
     !ZZW2 veetical activation velocity

     CALL CH_AER_ACTIVATION(ZAERO, ZZT, ZZW2, ZZTDT, ZRHODREF, ZPABST,&
                             ZNAT(:,1), ZMCN, ZSOLORG, ZMI, ZSMAX)

     ZZW1(:) = MAX(ZNATOLD(:)- ZNAT(:,1) , 0.0 )
!
     ZW(:,:) = UNPACK( ZZW1(:),MASK=GNUCT(:,:),FIELD=0.0 )
     PNAT(:,:,1) = PNAT(:,:,1) + ZW(:,:)
     ! Je sais pas ce que c'est:
     ZZW4(:) = 1.
     ZZW5(:) = 1.
   !
   !* prepare to update the cloud water concentration
   !
      ZZW6(:) = ZZW1(:) 
     DEALLOCATE(ZMCN)
     DEALLOCATE(ZNATOLD)
!
ELSE ! CCN activation from Cohard-Pinty
   ALLOCATE(ZSMAX(INUCT))
   IF (LIMAP%LADJ) THEN
      ZZW1(:) = 1.0/ZEPS + 1.0/ZZW1(:)                                   &
                + (((CST%XLVTT+(CST%XCPV-CST%XCL)*(ZZT(:)-CST%XTT))/ZZT(:))**2)/(CST%XCPD*CST%XRV) ! Psi2
!
!
!-------------------------------------------------------------------------------
!
!
!*       2. compute the constant term (ZZW3) relative to smax    
!        ----------------------------------------------------
!
!  Remark : in LIMA's nucleation parameterization, Smax=0.01 for a supersaturation of 1% !
!
!
      ZVEC1(:) = MAX( 1.0001, MIN( REAL(LIMAW%NAHEN)-0.0001, LIMAW%XAHENINTP1 * ZZT(:) + LIMAW%XAHENINTP2 ) )
      IVEC1(:) = INT( ZVEC1(:) )
      ZVEC1(:) = ZVEC1(:) - REAL( IVEC1(:) )
!
!
      IF (LIMAP%LACTIT) THEN ! including a cooling rate
!
!       Compute the tabulation of function of ZZW3 :
!
!                                                  (Psi1*w+Psi3*DT/Dt)**1.5 
!       ZZW3 = LIMAW%XAHENG*(Psi1*w + Psi3*DT/Dt)**1.5 = ------------------------ 
!                                                   2*pi*rho_l*G**(3/2)     
!
!
         ZZW4(:)=LIMAW%XPSI1( IVEC1(:)+1)*ZZW2(:)+LIMAW%XPSI3(IVEC1(:)+1)*ZZTDT(:)
         ZZW5(:)=LIMAW%XPSI1( IVEC1(:)  )*ZZW2(:)+LIMAW%XPSI3(IVEC1(:)  )*ZZTDT(:)
         WHERE (ZZW4(:) < 0. .OR. ZZW5(:) < 0.)
            ZZW4(:) = 0.
            ZZW5(:) = 0.
         END WHERE
         ZZW3(:) =   LIMAW%XAHENG( IVEC1(:)+1)*(ZZW4(:)**1.5)* ZVEC1(:)      &
                   - LIMAW%XAHENG( IVEC1(:)  )*(ZZW5(:)**1.5)*(ZVEC1(:) - 1.0)
                       ! Cste*((Psi1*w+Psi3*dT/dt)/(G))**1.5
         ZZW6(:) =   LIMAW%XAHENG2( IVEC1(:)+1)*(ZZW4(:)**0.5)* ZVEC1(:)      &
                   - LIMAW%XAHENG2( IVEC1(:)  )*(ZZW5(:)**0.5)*(ZVEC1(:) - 1.0)
!
!
      ELSE ! LIMAP%LACTIT , for clouds
!
!
!       Compute the tabulation of function of ZZW3 :
!
!                                             (Psi1 * w)**1.5       
!       ZZW3 = LIMAW%XAHENG * (Psi1 * w)**1.5  = -------------------------
!                                            2 pi rho_l * G**(3/2)  
!
!
         ZZW2(:)=MAX(ZZW2(:),0.)
         ZZW3(:)=LIMAW%XAHENG(IVEC1(:)+1)*((LIMAW%XPSI1(IVEC1(:)+1)*ZZW2(:))**1.5)* ZVEC1(:)    &
                -LIMAW%XAHENG(IVEC1(:)  )*((LIMAW%XPSI1(IVEC1(:)  )*ZZW2(:))**1.5)*(ZVEC1(:)-1.0)
!
         ZZW6(:)=LIMAW%XAHENG2(IVEC1(:)+1)*((LIMAW%XPSI1(IVEC1(:)+1)*ZZW2(:))**0.5)* ZVEC1(:)    &
                -LIMAW%XAHENG2(IVEC1(:)  )*((LIMAW%XPSI1(IVEC1(:)  )*ZZW2(:))**0.5)*(ZVEC1(:)-1.0)
!
      END IF ! LIMAP%LACTIT
!
!
!              (Psi1*w+Psi3*DT/Dt)**1.5   rho_air
!       ZZW3 = ------------------------ * -------
!                 2*pi*rho_l*G**(3/2)       Psi2
!
      ZZW5(:) = 1.
      ZZW3(:) = (ZZW3(:)/ZZW1(:))*ZRHODREF(:) ! R.H.S. of Eq 9 of CPB 98 but
      ! for multiple aerosol modes
      WHERE (ZRCT(:) > LIMAP%XRTMIN(2) .AND. ZCCT(:) > LIMAP%XCTMIN(2))
         ZZW6(:) = ZZW6(:) * ZRHODREF(:) * ZCCT(:) / (LIMAW%XLBC*ZCCT(:)/ZRCT(:))**LIMAW%XLBEXC
      ELSEWHERE
         ZZW6(:)=0.
      END WHERE

      WHERE (ZZW3(:) == 0. .AND. .NOT.(ZSW>0.))
         ZZW5(:) = -1.
      END WHERE
!
!-------------------------------------------------------------------------------
!
!
!*       3. Compute the maximum of supersaturation
!        -----------------------------------------
!
!
! estimate S_max for the CPB98 parameterization with SEVERAL aerosols mode
! Reminder : Smax=0.01 for a 1% supersaturation
!
! Interval bounds to tabulate sursaturation Smax
! Check with values used for tabulation in ini_lima_warm.f90
      ZS1 = 1.0E-5                   ! corresponds to  0.001% supersaturation
      ZS2 = 5.0E-2                   ! corresponds to 5.0% supersaturation 
      ZXACC = 1.0E-10                ! Accuracy needed for the search in [NO UNITS]
!
      ZSMAX(:) = ZRIDDR(ZS1,ZS2,ZXACC,ZZW3(:),ZZW6(:),INUCT)    ! ZSMAX(:) is in [NO UNITS]
      ZSMAX(:) = MIN(MAX(ZSMAX(:), ZSW(:)),ZS2)
      !
   ELSE
      ZSMAX(:) = ZSW(:)
      ZZW5(:) = 1.
   END IF
!
!-------------------------------------------------------------------------------
!
!
!*       4. Compute the nucleus source
!        -----------------------------
!
!
! Again : Smax=0.01 for a 1% supersaturation
! Modified values for Beta and C (see in init_aerosol_properties) account for that
!
   WHERE (ZZW5(:) > 0. .AND. ZSMAX(:) > 0.)
      ZVEC1(:) = MAX( 1.0001, MIN( REAL(LIMAW%NHYP)-0.0001, LIMAW%XHYPINTP1*LOG(ZSMAX(:))+LIMAW%XHYPINTP2 ) )
      IVEC1(:) = INT( ZVEC1(:) )
      ZVEC1(:) = ZVEC1(:) - REAL( IVEC1(:) )
   END WHERE
   ZZW6(:)  = 0. ! initialize the change of cloud droplet concentration
!
   ZTMP(:,:)=0.0
!
! Compute the concentration of activable aerosols for each mode
! based on the max of supersaturation ( -> ZTMP )
!
   DO IMOD = 1, LIMAP%NMOD_CCN                     ! iteration on mode number
      ZZW1(:) = 0.
      ZZW2(:) = 0.
      ZZW3(:) = 0.
   !
      WHERE( ZZW5(:) > 0. .AND. ZSMAX(:)>0.0 )
         ZZW2(:) =  LIMAW%XHYPF12( IVEC1(:)+1,IMOD )* ZVEC1(:)      & ! hypergeo function
                  - LIMAW%XHYPF12( IVEC1(:)  ,IMOD )*(ZVEC1(:) - 1.0) ! XHYPF12 is tabulated
   !
         ZTMP(:,IMOD) = ZCHEN_MULTI(:,IMOD)/ZRHODREF(:)*ZSMAX(:)**LIMAP%XKHEN_MULTI(IMOD)*ZZW2(:)
      ENDWHERE
   ENDDO
!
! Compute the concentration of aerosols activated at this time step
! as the difference between ZTMP and the aerosols already activated at t-dt (ZZW1)
!
   DO IMOD = 1, LIMAP%NMOD_CCN                     ! iteration on mode number
      ZZW1(:) = 0.
      ZZW2(:) = 0.
      ZZW3(:) = 0.
   !
      WHERE( SUM(ZTMP(:,:),DIM=2) .GT. 0.01E6/ZRHODREF(:) ) 
         ZZW1(:) = MIN( ZNFT(:,IMOD),MAX( ZTMP(:,IMOD)- ZNAT(:,IMOD) , 0.0 ) )
      ENDWHERE
   !
   !* update the concentration of activated CCN = Na
   !
      PNAT(:,:,IMOD) = PNAT(:,:,IMOD) + ZCLDFR(:,:) * UNPACK( ZZW1(:), MASK=GNUCT(:,:), FIELD=0.0 )
   !
   !* update the concentration of free CCN = Nf
   !
      PNFT(:,:,IMOD) = PNFT(:,:,IMOD) - ZCLDFR(:,:) * UNPACK( ZZW1(:), MASK=GNUCT(:,:), FIELD=0.0 )
   !
   !* prepare to update the cloud water concentration 
   !
      ZZW6(:) = ZZW6(:) + ZZW1(:)
   ENDDO
END IF ! AER_ACTIVATION
!
! Output tendencies
!
   ZZW1(:)=0.
   WHERE (ZZW5(:)>0.0 .AND. ZSMAX(:)>0.0) ! ZZW1 is computed with ZSMAX [NO UNIT]
      ZZW1(:) = MIN(LIMAW%XCSTDCRIT*ZZW6(:)/(((ZZT(:)*ZSMAX(:))**3)*ZRHODREF(:)),1.E-5)
   END WHERE
!
   IF(PRESENT(PTOT_RV_HENU)) PTOT_RV_HENU(:,:) = 0.
   IF (.NOT.NEBN%LSUBG_COND) THEN
      ZW(:,:) = MIN( UNPACK( ZZW1(:),MASK=GNUCT(:,:),FIELD=0.0 ),PRVT(:,:) )
      IF(PRESENT(PTOT_RV_HENU)) PTOT_RV_HENU(:,:) = ZW(:,:)
      PTHT(:,:) = PTHT(:,:) + ZW(:,:) * (CST%XLVTT+(CST%XCPV-CST%XCL)*(PT(:,:)-CST%XTT))/                &
            (PEXNREF(:,:)*(CST%XCPD+CST%XCPV*PRVT(:,:)+CST%XCL*(PRCT(:,:)+PRRT(:,:))))
      PRVT(:,:) = PRVT(:,:) - ZW(:,:) 
      PRCT(:,:) = PRCT(:,:) + ZW(:,:) 
      PCCT(:,:) = PCCT(:,:) + UNPACK( ZZW6(:),MASK=GNUCT(:,:),FIELD=0. ) 
   ELSE
      ZW(:,:) = MIN( ZCLDFR(:,:) * UNPACK( ZZW1(:),MASK=GNUCT(:,:),FIELD=0.0 ),PRVT(:,:) )
      PCCT(:,:) = PCCT(:,:) + ZCLDFR(:,:) * UNPACK( ZZW6(:),MASK=GNUCT(:,:),FIELD=0. ) 
   END IF
!
!++cb-- A quoi servent ces 2 dernieres lignes ? variables locales, non sauvees, et ne servent pas 
! a calculer quoi que ce soit (fin de la routine)
   ZW(:,:)   = UNPACK( 100.0*ZSMAX(:),MASK=GNUCT(:,:),FIELD=0.0 )
   ZW2(:,:)  = ZCLDFR(:,:) * UNPACK( ZZW6(:),MASK=GNUCT(:,:),FIELD=0.0 )
!
!
!-------------------------------------------------------------------------------
!
!
!*       5. Cleaning
!        -----------
!
!
   DEALLOCATE(IVEC1)
   DEALLOCATE(ZVEC1)
   DEALLOCATE(ZNFT)
   DEALLOCATE(ZNAT)
   DEALLOCATE(ZCCT)
   DEALLOCATE(ZRCT)
   DEALLOCATE(ZZT)
   DEALLOCATE(ZSMAX)
   DEALLOCATE(ZZW1)
   DEALLOCATE(ZZW2)
   DEALLOCATE(ZZW3)
   DEALLOCATE(ZZW4)
   DEALLOCATE(ZZW5)
   DEALLOCATE(ZZW6)
   DEALLOCATE(ZZTDT)
   DEALLOCATE(ZSW)
   DEALLOCATE(ZRHODREF)
   DEALLOCATE(ZCHEN_MULTI)
   DEALLOCATE(ZEXNREF)
   DEALLOCATE(ZPABST)
   DEALLOCATE(ZAERO)
   DEALLOCATE(ZSOLORG)
   DEALLOCATE(ZMI)
!
END IF ! INUCT
!
!!$IF ( tpfile%lopened ) THEN
!!$  IF ( INUCT == 0 ) THEN
!!$    ZW (:,:) = 0.
!!$    ZW2(:,:) = 0.
!!$  END IF
!!$
!!$  TZFIELD%CMNHNAME   ='SMAX'
!!$  TZFIELD%CSTDNAME   = ''
!!$  TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
!!$  TZFIELD%CUNITS     = ''
!!$  TZFIELD%CDIR       = 'XY'
!!$  TZFIELD%CCOMMENT   = 'X_Y_Z_SMAX'
!!$  TZFIELD%NGRID      = 1
!!$  TZFIELD%NTYPE      = TYPEREAL
!!$  TZFIELD%NDIMS      = 3
!!$  TZFIELD%LTIMEDEP   = .TRUE.
!!$  CALL IO_Field_write(TPFILE,TZFIELD,ZW)
!!$  !
!!$  TZFIELD%CMNHNAME   ='NACT'
!!$  TZFIELD%CSTDNAME   = ''
!!$  TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
!!$  TZFIELD%CUNITS     = 'kg-1'
!!$  TZFIELD%CDIR       = 'XY'
!!$  TZFIELD%CCOMMENT   = 'X_Y_Z_NACT'
!!$  TZFIELD%NGRID      = 1
!!$  TZFIELD%NTYPE      = TYPEREAL
!!$  TZFIELD%NDIMS      = 3
!!$  TZFIELD%LTIMEDEP   = .TRUE.
!!$  CALL IO_Field_write(TPFILE,TZFIELD,ZW2)
!!$END IF
!
!
!-------------------------------------------------------------------------------
!
!
!*       6. Functions used to compute the maximum of supersaturation
!        -----------------------------------------------------------
!
!
IF (LHOOK) CALL DR_HOOK('LIMA_CCN_ACTIVATION', 1, ZHOOK_HANDLE)
CONTAINS
!------------------------------------------------------------------------------
!
  FUNCTION ZRIDDR(PX1,PX2INIT,PXACC,PZZW3,PZZW6,KPTS)  RESULT(PZRIDDR)
!
!
!!****  *ZRIDDR* - iterative algorithm to find root of a function
!!
!!
!!    PURPOSE
!!    -------
!!       The purpose of this function is to find the root of a given function
!!     the arguments are the brackets bounds (the interval where to find the root)
!!     the accuracy needed and the input parameters of the given function.
!!     Using Ridders' method, return the root of a function known to lie between 
!!     PX1 and PX2. The root, returned as PZRIDDR, will be refined to an approximate
!!     accuracy PXACC.
!! 
!!**  METHOD
!!    ------
!!       Ridders' method
!!
!!    EXTERNAL
!!    --------
!!       FUNCSMAX  
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!      NUMERICAL RECIPES IN FORTRAN 77: THE ART OF SCIENTIFIC COMPUTING 
!!     (ISBN 0-521-43064-X)
!!      Copyright (C) 1986-1992 by Cambridge University Press.
!!      Programs Copyright (C) 1986-1992 by Numerical Recipes Software.
!!
!!    AUTHOR
!!    ------
!!      Frederick Chosson *CERFACS*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original     12/07/07
!!      S.BERTHET        2008 vectorization 
!------------------------------------------------------------------------------
!
!*       0. DECLARATIONS
!
!
USE MODE_MSG
USE YOMHOOK, ONLY:LHOOK, DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
!*       0.1 declarations of arguments and result
!
INTEGER,            INTENT(IN)     :: KPTS
REAL, DIMENSION(KPTS), INTENT(IN)     :: PZZW3
REAL, DIMENSION(KPTS), INTENT(IN)     :: PZZW6
REAL,               INTENT(IN)     :: PX1, PX2INIT, PXACC
REAL, DIMENSION(KPTS)    :: PZRIDDR
!
!*       0.2 declarations of local variables
!
!
INTEGER, PARAMETER                 :: JPMAXIT=60
REAL,    PARAMETER                 :: PPUNUSED=0.0 !-1.11e30
REAL,    DIMENSION(:), ALLOCATABLE :: ZFH,ZFL, ZFM,ZFNEW
REAL                               :: ZS,ZXH,ZXL,ZXM,ZXNEW
REAL                               :: ZX2
INTEGER                            :: IJ, IL
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('ZRIDDR', 0, ZHOOK_HANDLE)
ALLOCATE(  ZFH(KPTS))
ALLOCATE(  ZFL(KPTS))
ALLOCATE(  ZFM(KPTS))
ALLOCATE(ZFNEW(KPTS))
!
PZRIDDR(:)= PPUNUSED
ZX2       = PX2INIT 
ZFL(:)     = FUNCSMAX(PX1,PZZW3(:),PZZW6(:),KPTS)
ZFH(:)     = FUNCSMAX(ZX2,PZZW3(:),PZZW6(:),KPTS)
!
DO IL = 1, KPTS
   ZX2 = PX2INIT
100 IF ((ZFL(IL) > 0.0 .AND. ZFH(IL) < 0.0) .OR. (ZFL(IL) < 0.0 .AND. ZFH(IL) > 0.0)) then
      ZXL         = PX1
      ZXH         = ZX2
      DO IJ=1,JPMAXIT
         ZXM     = 0.5*(ZXL+ZXH)
         ZFM(IL) = SINGL_FUNCSMAX(ZXM,PZZW3(IL),PZZW6(IL),IL)
         ZS      = SQRT(ZFM(IL)**2-ZFL(IL)*ZFH(IL))
         IF (ZS == 0.0) then
            GO TO 101
         ENDIF
         ZXNEW  = ZXM+(ZXM-ZXL)*(SIGN(1.0,ZFL(IL)-ZFH(IL))*ZFM(IL)/ZS)
         IF (ABS(ZXNEW - PZRIDDR(IL)) <= PXACC) then
            GO TO 101 
         ENDIF
         PZRIDDR(IL) = ZXNEW
         ZFNEW(IL)  = SINGL_FUNCSMAX(PZRIDDR(IL),PZZW3(IL),PZZW6(IL),IL)
         IF (ZFNEW(IL) == 0.0) then
            GO TO 101
         ENDIF
         IF (SIGN(ZFM(IL),ZFNEW(IL)) /= ZFM(IL)) then
            ZXL    =ZXM
            ZFL(IL)=ZFM(IL)
            ZXH    =PZRIDDR(IL)
            ZFH(IL)=ZFNEW(IL)
         ELSE IF (SIGN(ZFL(IL),ZFNEW(IL)) /= ZFL(IL)) then
            ZXH    =PZRIDDR(IL)
            ZFH(IL)=ZFNEW(IL)
         ELSE IF (SIGN(ZFH(IL),ZFNEW(IL)) /= ZFH(IL)) then
            ZXL    =PZRIDDR(IL)
            ZFL(IL)=ZFNEW(IL)
         ELSE IF (ZX2 .LT. 0.05) then
            ZX2 = ZX2 + 1.0E-2
!            PRINT*, 'ZX2 ALWAYS too small, we put a greater one : ZX2 =',ZX2
            ZFH(IL)   = SINGL_FUNCSMAX(ZX2,PZZW3(IL),PZZW6(IL),IL)
            GO TO 100
         END IF
         IF (ABS(ZXH-ZXL) <= PXACC) then
            GO TO 101 
         ENDIF
!!SB
!!$      if (ij == JPMAXIT .and. (abs(zxh-zxl) > PXACC) ) then
!!$        PZRIDDR(IL)=0.0
!!$        go to 101
!!$      endif   
!!SB
      END DO
      CALL PRINT_MSG( NVERB_FATAL, 'GEN', 'ZRIDDR', 'EXCEEDED MAXIMUM ITERATIONS' )
   ELSE IF (ZFL(IL) == 0.0) then
      PZRIDDR(IL)=PX1
   ELSE IF (ZFH(IL) == 0.0) then
      PZRIDDR(IL)=ZX2
   ELSE IF (ZX2 .LT. 0.05) then
      ZX2 = ZX2 + 1.0E-2
!      PRINT*, 'ZX2 too small, we put a greater one : ZX2 =',ZX2
      ZFH(IL)   = SINGL_FUNCSMAX(ZX2,PZZW3(IL),PZZW6(IL),IL)
      GO TO 100
   ELSE
!!$      print*, 'PZRIDDR: root must be bracketed'
!!$      print*,'npts ',KPTS,'jl',IL
!!$      print*, 'PX1,ZX2,zfl,zfh',PX1,ZX2,zfl(IL),zfh(IL)
!!$      print*, 'ZX2 = 30 % of supersaturation, there is no solution for Smax'
!!$      print*, 'try to put greater ZX2 (upper bound for Smax research)'
!!$      STOP
      PZRIDDR(IL)=0.0
      GO TO 101
   END IF
101 ENDDO
!
DEALLOCATE(  ZFH)
DEALLOCATE(  ZFL)
DEALLOCATE(  ZFM)
DEALLOCATE(ZFNEW)
!
IF (LHOOK) CALL DR_HOOK('ZRIDDR', 1, ZHOOK_HANDLE)
END FUNCTION ZRIDDR
!
!------------------------------------------------------------------------------
!
  FUNCTION FUNCSMAX(PPZSMAX,PPZZW3,PPZZW6,KPTS)  RESULT(PFUNCSMAX)
!
USE YOMHOOK, ONLY:LHOOK, DR_HOOK, JPHOOK
!
!!****  *FUNCSMAX* - function describing SMAX function that you want to find the root
!!
!!
!!    PURPOSE
!!    -------
!!       This function describe the equilibrium between Smax and two aerosol mode
!!     acting as CCN. This function is derive from eq. (9) of CPB98 but for two
!!     aerosols mode described by their respective parameters C, k, Mu, Beta.
!!     the arguments are the supersaturation in "no unit" and the r.h.s. of this eq.
!!     and the ratio of concentration of injected aerosols on maximum concentration
!!     of injected aerosols ever.
!!**  METHOD
!!    ------
!!       This function is called by zriddr.f90
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!    Module MODD_PARAM_LIMA_WARM
!!        LIMAW%XHYPF32
!!
!!        LIMAW%XHYPINTP1
!!        LIMAW%XHYPINTP2
!!
!!    Module MODD_PARAM_C2R2
!!        LIMAP%XKHEN_MULTI()
!!        LIMAP%NMOD_CCN
!!       
!!    REFERENCE
!!    ---------
!!    Cohard, J.M., J.P.Pinty, K.Suhre, 2000:"On the parameterization of activation
!!             spectra from cloud condensation nuclei microphysical properties",
!!             J. Geophys. Res., Vol.105, N0.D9, pp. 11753-11766
!!
!!    AUTHOR
!!    ------
!!      Frederick Chosson *CERFACS*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original     12/07/07
!!      S.Berthet    19/03/08 Extension a une population multimodale d aerosols
!
!------------------------------------------------------------------------------
!
!*       0. DECLARATIONS
!
IMPLICIT NONE
!
!*       0.1 declarations of arguments and result
!
INTEGER,            INTENT(IN)  :: KPTS
REAL,               INTENT(IN)  :: PPZSMAX   ! supersaturation is already in no units
REAL, DIMENSION(KPTS), INTENT(IN)  :: PPZZW3    ! 
REAL, DIMENSION(KPTS), INTENT(IN)  :: PPZZW6    ! 
REAL, DIMENSION(KPTS) :: PFUNCSMAX ! 
!
!*       0.2 declarations of local variables
!
REAL                           :: ZHYPF
!
REAL                           :: ZVEC1
INTEGER                        :: IVEC1
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('FUNCSMAX', 0, ZHOOK_HANDLE)
PFUNCSMAX(:) = 0.
ZVEC1 = MAX( ( 1.0 + 10.0 * CST%XMNH_EPSILON ) ,MIN( REAL(LIMAW%NHYP)*( 1.0 - 10.0 * CST%XMNH_EPSILON ) ,               &
                           LIMAW%XHYPINTP1*LOG(PPZSMAX)+LIMAW%XHYPINTP2 ) )
IVEC1 = INT( ZVEC1 )
ZVEC1 = ZVEC1 - REAL( IVEC1 )
DO IMOD = 1, LIMAP%NMOD_CCN
   ZHYPF        = 0.          ! LIMAW%XHYPF32 is tabulated with ZSMAX in [NO UNITS]
   ZHYPF        =   LIMAW%XHYPF32( IVEC1+1,IMOD ) * ZVEC1              &
                  - LIMAW%XHYPF32( IVEC1  ,IMOD ) *(ZVEC1 - 1.0)
                             ! sum of s**(ki+2) * F32 * Ci * ki * beta(ki/2,3/2)
   PFUNCSMAX(:) =  PFUNCSMAX(:) + (PPZSMAX)**(LIMAP%XKHEN_MULTI(IMOD) + 2) &
                 * ZHYPF* LIMAP%XKHEN_MULTI(IMOD) * ZCHEN_MULTI(:,IMOD)    &
                 * LIMAP%XGMULTI(IMOD)
ENDDO
! function l.h.s. minus r.h.s. of eq. (9) of CPB98 but for LIMAP%NMOD_CCN aerosol mode
PFUNCSMAX(:) = PFUNCSMAX(:) + PPZZW6(:)*PPZSMAX - PPZZW3(:)
!
IF (LHOOK) CALL DR_HOOK('FUNCSMAX', 1, ZHOOK_HANDLE)
END FUNCTION FUNCSMAX
!
!------------------------------------------------------------------------------
!
  FUNCTION SINGL_FUNCSMAX(PPZSMAX,PPZZW3,PPZZW6,KINDEX)  RESULT(PSINGL_FUNCSMAX)
USE YOMHOOK, ONLY:LHOOK, DR_HOOK, JPHOOK
!
!
!!****  *SINGL_FUNCSMAX* - same function as FUNCSMAX
!!
!!
!!    PURPOSE
!!    -------
!        As for FUNCSMAX but for a scalar
!!
!!**  METHOD
!!    ------
!!       This function is called by zriddr.f90
!!
!------------------------------------------------------------------------------
!
!*       0. DECLARATIONS
!
IMPLICIT NONE
!
!*       0.1 declarations of arguments and result
!
INTEGER,            INTENT(IN)  :: KINDEX
REAL,               INTENT(IN)  :: PPZSMAX   ! supersaturation is "no unit"
REAL,               INTENT(IN)  :: PPZZW3    ! 
REAL,               INTENT(IN)  :: PPZZW6    ! 
REAL                            :: PSINGL_FUNCSMAX ! 
!
!*       0.2 declarations of local variables
!
REAL                           :: ZHYPF
!
REAL                           :: ZVEC1
INTEGER                        :: IVEC1
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('SINGL_FUNCSMAX', 0, ZHOOK_HANDLE)
PSINGL_FUNCSMAX = 0.
ZVEC1    = MAX( 1.0001,MIN( REAL(LIMAW%NHYP)-0.0001,               &
                              LIMAW%XHYPINTP1*LOG(PPZSMAX)+LIMAW%XHYPINTP2 ) )
IVEC1 = INT( ZVEC1 )
ZVEC1 = ZVEC1 - REAL( IVEC1 )
DO IMOD = 1, LIMAP%NMOD_CCN
   ZHYPF        = 0.          ! LIMAW%XHYPF32 is tabulated with ZSMAX in [NO UNITS]
   ZHYPF        =   LIMAW%XHYPF32( IVEC1+1,IMOD ) * ZVEC1              &
                  - LIMAW%XHYPF32( IVEC1  ,IMOD ) *(ZVEC1 - 1.0)
                             ! sum of s**(ki+2) * F32 * Ci * ki * bêta(ki/2,3/2)
   PSINGL_FUNCSMAX = PSINGL_FUNCSMAX + (PPZSMAX)**(LIMAP%XKHEN_MULTI(IMOD) + 2)   &
                   * ZHYPF* LIMAP%XKHEN_MULTI(IMOD) * ZCHEN_MULTI(KINDEX,IMOD) &
                   * LIMAP%XGMULTI(IMOD)
ENDDO
! function l.h.s. minus r.h.s. of eq. (9) of CPB98 but for LIMAP%NMOD_CCN aerosol mode
PSINGL_FUNCSMAX = PSINGL_FUNCSMAX + PPZZW6*PPZSMAX - PPZZW3
!
IF (LHOOK) CALL DR_HOOK('SINGL_FUNCSMAX', 1, ZHOOK_HANDLE)
END FUNCTION SINGL_FUNCSMAX
!
!-----------------------------------------------------------------------------
!
END SUBROUTINE LIMA_CCN_ACTIVATION
END MODULE MODE_LIMA_CCN_ACTIVATION
