!MNH_LIC Copyright 2018-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-------------------------------------------------------------------------------
!      ########################################
       MODULE MODI_LIMA_COLLISIONAL_ICE_BREAKUP
!      ########################################
!
INTERFACE
   SUBROUTINE LIMA_COLLISIONAL_ICE_BREAKUP (LDCOMPUTE,              &
                                            PRHODREF,               &
                                            PRIT, PRST, PRGT, PCIT, PCST, PCGT, &
                                            PLBDS, PLBDG,           &
                                            P_RI_CIBU, P_CI_CIBU    )
!
LOGICAL, DIMENSION(:),INTENT(IN)    :: LDCOMPUTE
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRHODREF
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRIT
REAL, DIMENSION(:),   INTENT(IN)    :: PRST
REAL, DIMENSION(:),   INTENT(IN)    :: PRGT
REAL, DIMENSION(:),   INTENT(IN)    :: PCIT
REAL, DIMENSION(:),   INTENT(IN)    :: PCST
REAL, DIMENSION(:),   INTENT(IN)    :: PCGT
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDS 
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDG 
!
REAL, DIMENSION(:),   INTENT(OUT)   :: P_RI_CIBU
REAL, DIMENSION(:),   INTENT(OUT)   :: P_CI_CIBU
!
END SUBROUTINE LIMA_COLLISIONAL_ICE_BREAKUP
END INTERFACE
END MODULE MODI_LIMA_COLLISIONAL_ICE_BREAKUP
!
!     #######################################################################
      SUBROUTINE LIMA_COLLISIONAL_ICE_BREAKUP (LDCOMPUTE,              &
                                               PRHODREF,               &
                                               PRIT, PRST, PRGT, PCIT, PCST, PCGT, &
                                               PLBDS, PLBDG,           &
                                               P_RI_CIBU, P_CI_CIBU    )
!     #######################################################################
!
!!    PURPOSE
!!    -------
!!      Compute the collisional ice break-up (secondary ice production process)
!!
!!    AUTHOR
!!    ------
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!      T.    Hoarau     * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/2022  duplicate from original for LIMA_SPLIT
!       B.    Vié   04/2022  Adapt to the new snow characteristics        
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAM_LIMA,       ONLY : LCIBU, XRTMIN, XCTMIN, XCEXVT, XALPHAS, XNUS, XNDEBRIS_CIBU

USE MODD_PARAM_LIMA_COLD,  ONLY : XBS, XCS, XDS, XFVELOS, XMNU0
USE MODD_PARAM_LIMA_MIXED, ONLY : XCG, XDG, XCXG,                              &
                                  XCIBUINTP_S, XCIBUINTP1_S, XCIBUINTP2_S,     &
                                  XCIBUINTP_G, XCIBUINTP1_G,                   &
                                  XFACTOR_CIBU_NI, XFACTOR_CIBU_RI,            &
                                  XMOMGG_CIBU_1, XMOMGG_CIBU_2,                &
                                  XMOMGS_CIBU_1, XMOMGS_CIBU_2, XMOMGS_CIBU_3, &
                                  NGAMINC, XGAMINC_CIBU_S, XGAMINC_CIBU_G
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
LOGICAL, DIMENSION(:),INTENT(IN)    :: LDCOMPUTE
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRHODREF
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRIT
REAL, DIMENSION(:),   INTENT(IN)    :: PRST
REAL, DIMENSION(:),   INTENT(IN)    :: PRGT
REAL, DIMENSION(:),   INTENT(IN)    :: PCIT
REAL, DIMENSION(:),   INTENT(IN)    :: PCST
REAL, DIMENSION(:),   INTENT(IN)    :: PCGT
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDS 
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDG 
!
REAL, DIMENSION(:),   INTENT(OUT)   :: P_RI_CIBU
REAL, DIMENSION(:),   INTENT(OUT)   :: P_CI_CIBU
!
!
!*       0.2   Declarations of local variables :
!
LOGICAL, DIMENSION(SIZE(PRST)) :: GCIBU ! Test where to compute collision process
LOGICAL, SAVE                 :: GFIRSTCALL = .TRUE. ! control switch for the first call
!
INTEGER                            :: ICIBU
INTEGER, DIMENSION(:), ALLOCATABLE :: IVEC2_S1,IVEC2_S2         ! Snow indice vector
INTEGER, DIMENSION(:), ALLOCATABLE :: IVEC2_G                   ! Graupel indice vector
INTEGER, PARAMETER                 :: I_SEED_PARAM = 26032012
INTEGER, DIMENSION(:), ALLOCATABLE :: I_SEED
INTEGER                            :: NI_SEED
!
REAL,    DIMENSION(:), ALLOCATABLE :: ZVEC1_S, ZVEC1_SW, ZVEC1_S1, ZVEC1_S2,  & ! Work vectors
                                      ZVEC1_S3, ZVEC1_S4,           &
                                      ZVEC1_S11, ZVEC1_S12,         & ! for snow
                                      ZVEC1_S21, ZVEC1_S22,         &
                                      ZVEC1_S31, ZVEC1_S32,         &
                                      ZVEC1_S41, ZVEC1_S42,         &
                                      ZVEC2_S1, ZVEC2_S2
REAL,    DIMENSION(:), ALLOCATABLE :: ZVEC1_G, ZVEC1_G1, ZVEC1_G2, & ! Work vectors
                                      ZVEC2_G                        ! for graupel
REAL,    DIMENSION(:), ALLOCATABLE :: ZFRAGMENTS, ZHARVEST
REAL,    DIMENSION(SIZE(PRST))     :: ZINTG_SNOW_1, & ! incomplete gamma function
                                      ZINTG_SNOW_2, & ! for snow
                                      ZINTG_SNOW_3, &
                                      ZINTG_SNOW_4
REAL,    DIMENSION(SIZE(PRST))     :: ZINTG_GRAUPEL_1, &  ! incomplete gamma
                                      ZINTG_GRAUPEL_2     ! function for graupel
REAL,    DIMENSION(SIZE(PRST))     :: ZNI_CIBU, ZRI_CIBU  ! CIBU rates
REAL,    DIMENSION(SIZE(PRST))     :: ZFRAG_CIBU
REAL                               :: ZFACT1_XNDEBRIS, ZFACT2_XNDEBRIS
!
!-------------------------------------------------------------------------------

GCIBU(:) = LCIBU        .AND. PRST(:)>XRTMIN(5) .AND. PRGT(:)>XRTMIN(6) .AND. &
           LDCOMPUTE(:) .AND. PCST(:)>XCTMIN(5) .AND. PCGT(:)>XCTMIN(6)
ICIBU    = COUNT( GCIBU(:) )
!
P_RI_CIBU(:)=0.
P_CI_CIBU(:)=0.
!
IF (ICIBU > 0) THEN
!
!       1.3.0 randomization of XNDEBRIS_CIBU values
!
  IF (GFIRSTCALL) THEN
    CALL RANDOM_SEED(SIZE=NI_SEED) ! get size of seed
    ALLOCATE(I_SEED(NI_SEED))
    I_SEED(:) = I_SEED_PARAM !
    CALL RANDOM_SEED(PUT=I_SEED)
    GFIRSTCALL = .FALSE.
  END IF
!
  ALLOCATE(ZFRAGMENTS(ICIBU))
!
  IF (XNDEBRIS_CIBU >= 0.0) THEN
    ZFRAGMENTS(:) = XNDEBRIS_CIBU
  ELSE
!
! Mantissa gives the mean value (randomization around 10**MANTISSA)
! First digit after the comma provides the full range around 10**MANTISSA
!
    ALLOCATE(ZHARVEST(ICIBU))
!
    ZFACT1_XNDEBRIS = AINT(XNDEBRIS_CIBU)
    ZFACT2_XNDEBRIS = ABS(ANINT(10.0*(XNDEBRIS_CIBU - ZFACT1_XNDEBRIS)))
!
    CALL RANDOM_NUMBER(ZHARVEST(:))
!
    ZFRAGMENTS(:) = 10.0**(ZFACT2_XNDEBRIS*ZHARVEST(:) + ZFACT1_XNDEBRIS)
!
    DEALLOCATE(ZHARVEST)
!
! ZFRAGMENTS is a random variable containing the number of fragments per collision
! For XNDEBRIS_CIBU=-1.2345  => ZFRAGMENTS(:) = 10.0**(2.0*RANDOM_NUMBER(ZHARVEST(:)) - 1.0)
! and ZFRAGMENTS=[0.1, 10.0] centered around 1.0
!
  END IF
!
!
!       1.3.1 To compute the partial integration of snow gamma function
!
!       1.3.1.0 allocations
!
  ALLOCATE(ZVEC1_S(ICIBU))
  ALLOCATE(ZVEC1_SW(ICIBU))
  ALLOCATE(ZVEC1_S1(ICIBU))
  ALLOCATE(ZVEC1_S2(ICIBU))
  ALLOCATE(ZVEC1_S3(ICIBU))
  ALLOCATE(ZVEC1_S4(ICIBU))
  ALLOCATE(ZVEC1_S11(ICIBU))
  ALLOCATE(ZVEC1_S12(ICIBU))
  ALLOCATE(ZVEC1_S21(ICIBU))
  ALLOCATE(ZVEC1_S22(ICIBU))
  ALLOCATE(ZVEC1_S31(ICIBU))
  ALLOCATE(ZVEC1_S32(ICIBU))
  ALLOCATE(ZVEC1_S41(ICIBU))
  ALLOCATE(ZVEC1_S42(ICIBU))
  ALLOCATE(ZVEC2_S1(ICIBU))
  ALLOCATE(IVEC2_S1(ICIBU))
  ALLOCATE(ZVEC2_S2(ICIBU))
  ALLOCATE(IVEC2_S2(ICIBU))
!
!
!       1.3.1.1 select the PLBDS
!
  ZVEC1_S(:) = PACK( PLBDS(:),MASK=GCIBU(:) )
  ZVEC1_SW(:)= ( XFVELOS**XALPHAS + ZVEC1_S(:)**XALPHAS ) ** (1./XALPHAS) ! modified equivalent lambda 
!
!
!       1.3.1.2 find the next lower indice for the PLBDS in the
!               geometrical set of Lbda_s used to tabulate some moments of the
!               incomplete gamma function, for boundary 1 (0.2 mm)
!
  ZVEC2_S1(1:ICIBU) = MAX( 1.0001, MIN( FLOAT(NGAMINC)-0.0001,XCIBUINTP_S  &
                      * LOG( ZVEC1_S(1:ICIBU) ) + XCIBUINTP1_S  ) )
  IVEC2_S1(1:ICIBU) = INT( ZVEC2_S1(1:ICIBU) )
  ZVEC2_S1(1:ICIBU) = ZVEC2_S1(1:ICIBU) - FLOAT( IVEC2_S1(1:ICIBU) )
!
!
!       1.3.1.3 find the next lower indice for the PLBDS in the
!               geometrical set of Lbda_s used to tabulate some moments of the
!               incomplete gamma function, for boundary 2 (1 mm)
!
  ZVEC2_S2(1:ICIBU) = MAX( 1.0001, MIN( FLOAT(NGAMINC)-0.0001,XCIBUINTP_S  &
                      * LOG( ZVEC1_S(1:ICIBU) ) + XCIBUINTP2_S  ) )
  IVEC2_S2(1:ICIBU) = INT( ZVEC2_S2(1:ICIBU) )
  ZVEC2_S2(1:ICIBU) = ZVEC2_S2(1:ICIBU) - FLOAT( IVEC2_S2(1:ICIBU) )
!
!
!       1.3.1.4 perform the linear interpolation of the
!               normalized "0"-moment of the incomplete gamma function
!
! For lower boundary (0.2 mm)
  ZVEC1_S11(1:ICIBU) = XGAMINC_CIBU_S(1,IVEC2_S1(1:ICIBU)+1) *  ZVEC2_S1(1:ICIBU) &
                     - XGAMINC_CIBU_S(1,IVEC2_S1(1:ICIBU))   * (ZVEC2_S1(1:ICIBU)-1.0)
!
! For upper boundary (1 mm)
  ZVEC1_S12(1:ICIBU) = XGAMINC_CIBU_S(1,IVEC2_S2(1:ICIBU)+1) *  ZVEC2_S2(1:ICIBU) &
                     - XGAMINC_CIBU_S(1,IVEC2_S2(1:ICIBU))   * (ZVEC2_S2(1:ICIBU)-1.0)
!
! Computation of spectrum from 0.2 mm to 1 mm
  ZVEC1_S1(1:ICIBU) = ZVEC1_S12(1:ICIBU) - ZVEC1_S11(1:ICIBU)
!
!
!       1.3.1.5 perform the linear interpolation of the
!               normalized "XBS"-moment of the incomplete gamma function
!
! For lower boundary (0.2 mm)
  ZVEC1_S31(1:ICIBU) = XGAMINC_CIBU_S(3,IVEC2_S1(1:ICIBU)+1) *  ZVEC2_S1(1:ICIBU) &
                     - XGAMINC_CIBU_S(3,IVEC2_S1(1:ICIBU))   * (ZVEC2_S1(1:ICIBU)-1.0)
!
! For upper boundary (1 mm)
  ZVEC1_S32(1:ICIBU) = XGAMINC_CIBU_S(3,IVEC2_S2(1:ICIBU)+1) *  ZVEC2_S2(1:ICIBU) &
                     - XGAMINC_CIBU_S(3,IVEC2_S2(1:ICIBU))   * (ZVEC2_S2(1:ICIBU)-1.0)
!
! From 0.2 mm to 1 mm we need
  ZVEC1_S3(1:ICIBU) = XMOMGS_CIBU_2 * (ZVEC1_S32(1:ICIBU) - ZVEC1_S31(1:ICIBU))
!
!
!       1.3.1.2 find the next lower indice for the PLBDS in the
!               geometrical set of Lbda_s used to tabulate some moments of the
!               incomplete gamma function, for boundary 1 (0.2 mm) for modified lambda (Wurtz snow fall speed)
!
  ZVEC2_S1(1:ICIBU) = MAX( 1.0001, MIN( FLOAT(NGAMINC)-0.0001,XCIBUINTP_S  &
                      * LOG( ZVEC1_SW(1:ICIBU) ) + XCIBUINTP1_S  ) )
  IVEC2_S1(1:ICIBU) = INT( ZVEC2_S1(1:ICIBU) )
  ZVEC2_S1(1:ICIBU) = ZVEC2_S1(1:ICIBU) - FLOAT( IVEC2_S1(1:ICIBU) )
!
!
!       1.3.1.3 find the next lower indice for the PLBDS in the
!               geometrical set of Lbda_s used to tabulate some moments of the
!               incomplete gamma function, for boundary 2 (1 mm) for modified lambda (Wurtz snow fall speed)
!
  ZVEC2_S2(1:ICIBU) = MAX( 1.0001, MIN( FLOAT(NGAMINC)-0.0001,XCIBUINTP_S  &
                      * LOG( ZVEC1_SW(1:ICIBU) ) + XCIBUINTP2_S  ) )
  IVEC2_S2(1:ICIBU) = INT( ZVEC2_S2(1:ICIBU) )
  ZVEC2_S2(1:ICIBU) = ZVEC2_S2(1:ICIBU) - FLOAT( IVEC2_S2(1:ICIBU) )
!
!
!       1.3.1.5 perform the linear interpolation of the
!               normalized "XDS"-moment of the incomplete gamma function
!
! For lower boundary (0.2 mm)
  ZVEC1_S21(1:ICIBU) = XGAMINC_CIBU_S(2,IVEC2_S1(1:ICIBU)+1) *  ZVEC2_S1(1:ICIBU) &
                     - XGAMINC_CIBU_S(2,IVEC2_S1(1:ICIBU))   * (ZVEC2_S1(1:ICIBU)-1.0)
!
! For upper boundary (1 mm)
  ZVEC1_S22(1:ICIBU) = XGAMINC_CIBU_S(2,IVEC2_S2(1:ICIBU)+1) *  ZVEC2_S2(1:ICIBU) &
                     - XGAMINC_CIBU_S(2,IVEC2_S2(1:ICIBU))   * (ZVEC2_S2(1:ICIBU)-1.0)
!
! From 0.2 mm to 1 mm we need
  ZVEC1_S2(1:ICIBU) = XMOMGS_CIBU_1 * (ZVEC1_S22(1:ICIBU) - ZVEC1_S21(1:ICIBU))
!
!
!       1.3.1.6 perform the linear interpolation of the
!               normalized "XBS+XDS"-moment of the incomplete gamma function
!
! For lower boundary (0.2 mm)
  ZVEC1_S41(1:ICIBU) = XGAMINC_CIBU_S(4,IVEC2_S1(1:ICIBU)+1) *  ZVEC2_S1(1:ICIBU) &
                     - XGAMINC_CIBU_S(4,IVEC2_S1(1:ICIBU))   * (ZVEC2_S1(1:ICIBU)-1.0)
!
! For upper boundary (1 mm)
  ZVEC1_S42(1:ICIBU) = XGAMINC_CIBU_S(4,IVEC2_S2(1:ICIBU)+1) *  ZVEC2_S2(1:ICIBU) &
                     - XGAMINC_CIBU_S(4,IVEC2_S2(1:ICIBU))   * (ZVEC2_S2(1:ICIBU)-1.0)
!
! From 0.2 mm to 1 mm we need
  ZVEC1_S4(1:ICIBU) = XMOMGS_CIBU_3 * (ZVEC1_S42(1:ICIBU) - ZVEC1_S41(1:ICIBU))
!
  ZINTG_SNOW_1(:) = UNPACK ( VECTOR=ZVEC1_S1(:),MASK=GCIBU,FIELD=0.0 )
  ZINTG_SNOW_2(:) = UNPACK ( VECTOR=ZVEC1_S2(:),MASK=GCIBU,FIELD=0.0 )
  ZINTG_SNOW_3(:) = UNPACK ( VECTOR=ZVEC1_S3(:),MASK=GCIBU,FIELD=0.0 )
  ZINTG_SNOW_4(:) = UNPACK ( VECTOR=ZVEC1_S4(:),MASK=GCIBU,FIELD=0.0 )
!
!
!       1.3.2 Compute the partial integration of graupel gamma function
!
!       1.3.2.0 allocations
!
  ALLOCATE(ZVEC1_G(ICIBU))
  ALLOCATE(ZVEC1_G1(ICIBU))
  ALLOCATE(ZVEC1_G2(ICIBU))
  ALLOCATE(ZVEC2_G(ICIBU))
  ALLOCATE(IVEC2_G(ICIBU))
!
!
!       1.3.2.1 select the PLBDG
!
  ZVEC1_G(:) = PACK( PLBDG(:),MASK=GCIBU(:) )
!
!
!       1.3.2.2 find the next lower indice for the PLBDG in the
!               geometrical set of Lbda_g used to tabulate some moments of the
!               incomplete gamma function, for the "2mm" boundary
!
  ZVEC2_G(1:ICIBU) = MAX( 1.0001, MIN( FLOAT(NGAMINC)-0.0001,XCIBUINTP_G  &
                     * LOG( ZVEC1_G(1:ICIBU) ) + XCIBUINTP1_G  ) )
  IVEC2_G(1:ICIBU) = INT( ZVEC2_G(1:ICIBU) )
  ZVEC2_G(1:ICIBU) = ZVEC2_G(1:ICIBU) - FLOAT( IVEC2_G(1:ICIBU) )
!
!
!       1.3.2.3 perform the linear interpolation of the
!               normalized "2+XDG"-moment of the incomplete gamma function
!
  ZVEC1_G1(1:ICIBU) = XGAMINC_CIBU_G(1,IVEC2_G(1:ICIBU)+1) *  ZVEC2_G(1:ICIBU)    &
                    - XGAMINC_CIBU_G(1,IVEC2_G(1:ICIBU))   * (ZVEC2_G(1:ICIBU)-1.0)
!
! From 2 mm to infinity we need
  ZVEC1_G1(1:ICIBU) = XMOMGG_CIBU_1 * (1.0 - ZVEC1_G1(1:ICIBU))
!
!
!       1.3.2.4 perform the linear interpolation of the
!               normalized "2.0"-moment of the incomplete gamma function
!
  ZVEC1_G2(1:ICIBU) = XGAMINC_CIBU_G(2,IVEC2_G(1:ICIBU)+1) *  ZVEC2_G(1:ICIBU)    &
                    - XGAMINC_CIBU_G(2,IVEC2_G(1:ICIBU))   * (ZVEC2_G(1:ICIBU)-1.0)
!
! From 2 mm to infinity we need
  ZVEC1_G2(1:ICIBU) = XMOMGG_CIBU_2 * (1.0 - ZVEC1_G2(1:ICIBU))
!
!
  ZINTG_GRAUPEL_1(:) = UNPACK ( VECTOR=ZVEC1_G1(:),MASK=GCIBU,FIELD=0.0 )
  ZINTG_GRAUPEL_2(:) = UNPACK ( VECTOR=ZVEC1_G2(:),MASK=GCIBU,FIELD=0.0 )
!
!
!        1.3.3 To compute final "CIBU" contributions
!
  ZFRAG_CIBU(:) = UNPACK ( VECTOR=ZFRAGMENTS(:),MASK=GCIBU,FIELD=0.0 )
  ZNI_CIBU(:) = ZFRAG_CIBU(:) * (XFACTOR_CIBU_NI * PCST(:) * PCGT(:) / (PRHODREF(:)**XCEXVT)) * &
                (XCG * ZINTG_GRAUPEL_1(:) * ZINTG_SNOW_1(:) *                                               &
                 PLBDG(:)**(-(XDG+2.0))                                             &
               - XCS * ZINTG_GRAUPEL_2(:) * ZINTG_SNOW_2(:) *                                               &
                 PLBDS(:)**(-XDS) * PLBDG(:)**(-2.0) *                                            &
                 (1+(XFVELOS/PLBDS(:))**XALPHAS)**(-XNUS-XDS/XALPHAS) )

  P_CI_CIBU(:) = MAX(ZNI_CIBU(:),0.)
!
  DEALLOCATE(ZFRAGMENTS)
!
! Max value of rs removed by CIBU
  ZRI_CIBU(:) = (XFACTOR_CIBU_RI * PCST(:) * PCGT(:) / (PRHODREF(:)**XCEXVT)) * &
                 (XCG * ZINTG_GRAUPEL_1(:) * ZINTG_SNOW_3(:) *                              &
                  PLBDS(:)**(-XBS) * PLBDG(:)**(-(XDG+2.0))                                               &
                - XCS * ZINTG_GRAUPEL_2(:) * ZINTG_SNOW_4(:) *                              &
                  PLBDS(:)**(-XBS-XDS) * PLBDG(:)**(-2.0) *                               &
                  (1+(XFVELOS/PLBDS(:))**XALPHAS)**(-XNUS-(XBS+XDS)/XALPHAS) )
!
! The value of rs removed by CIBU is determined by the mean mass of pristine ice
  WHERE( PRIT(:)>XRTMIN(4) .AND. PCIT(:)>XCTMIN(4) )
    ZRI_CIBU(:) = MIN( ZRI_CIBU(:), ZNI_CIBU(:)*PRIT(:)/PCIT(:) )
  ELSE WHERE
    ZRI_CIBU(:) = MIN( ZRI_CIBU(:), MAX( ZNI_CIBU(:)*XMNU0,XRTMIN(4) ) )
  END WHERE
!
  P_RI_CIBU(:) = MAX(ZRI_CIBU(:), 0.)
!
  DEALLOCATE(ZVEC1_S)
  DEALLOCATE(ZVEC1_SW)
  DEALLOCATE(ZVEC1_S1)
  DEALLOCATE(ZVEC1_S2)
  DEALLOCATE(ZVEC1_S3)
  DEALLOCATE(ZVEC1_S4)
  DEALLOCATE(ZVEC1_S11)
  DEALLOCATE(ZVEC1_S12)
  DEALLOCATE(ZVEC1_S21)
  DEALLOCATE(ZVEC1_S22)
  DEALLOCATE(ZVEC1_S31)
  DEALLOCATE(ZVEC1_S32)
  DEALLOCATE(ZVEC1_S41)
  DEALLOCATE(ZVEC1_S42)
  DEALLOCATE(ZVEC2_S1)
  DEALLOCATE(IVEC2_S1)
  DEALLOCATE(ZVEC2_S2)
  DEALLOCATE(IVEC2_S2)
  DEALLOCATE(ZVEC1_G)
  DEALLOCATE(ZVEC1_G1)
  DEALLOCATE(ZVEC1_G2)
  DEALLOCATE(ZVEC2_G)
  DEALLOCATE(IVEC2_G)
END IF
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_COLLISIONAL_ICE_BREAKUP
