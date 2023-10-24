!MNH_LIC Copyright 2010-2020 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ###########################
      MODULE MODI_ION_ATTACH_ELEC
!     ###########################
!
INTERFACE
        SUBROUTINE ION_ATTACH_ELEC(KTCOUNT, KRR, HCLOUD, PTSTEP, PRHODREF, &
                                   PRHODJ, PSVS, PRS, PTHT, PCIT, PPABST,  &
                                   PEFIELDU, PEFIELDV, PEFIELDW, GATTACH,  &
                                   PTOWN, PSEA,                            &
                                   PCCS, PCRS, PCSS, PCGS, PCHS            )
!
INTEGER,                             INTENT(IN)    :: KTCOUNT  ! Temporal loop counter
INTEGER,                             INTENT(IN)    :: KRR      ! Number of moist variables
CHARACTER(LEN=4),                    INTENT(IN)    :: HCLOUD   ! kind of cloud paramerization
REAL,                                INTENT(IN)    :: PTSTEP   ! Time step          
!
REAL,    DIMENSION(:,:,:),           INTENT(IN)    :: PRHODREF ! Reference dry air density 
REAL,    DIMENSION(:,:,:),           INTENT(IN)    :: PRHODJ   ! dry air density* Jacobian
REAL,    DIMENSION(:,:,:,:),         INTENT(INOUT) :: PSVS     ! Scalar variable vol. source
REAL,    DIMENSION(:,:,:,:),         INTENT(INOUT) :: PRS      ! Moist variable vol. source
REAL,    DIMENSION(:,:,:),           INTENT(IN)    :: PTHT     ! Theta (K) at time t
REAL,    DIMENSION(:,:,:),           INTENT(IN)    :: PCIT     ! Pristine ice n.c.
                                                               ! - at time t (for ICE schemes)
                                                               ! - source (for LIMA)
REAL,    DIMENSION(:,:,:),           INTENT(IN)    :: PPABST   ! Absolute pressure at t
REAL,    DIMENSION(:,:,:),           INTENT(IN)    :: PEFIELDU, PEFIELDV, PEFIELDW
                                                               ! Electric field components
LOGICAL, DIMENSION(:,:,:),           INTENT(IN)    :: GATTACH  ! Recombination and
                                                               ! Attachment if true
REAL,    DIMENSION(:,:),   OPTIONAL, INTENT(IN)    :: PTOWN    ! Town fraction
REAL,    DIMENSION(:,:),   OPTIONAL, INTENT(IN)    :: PSEA     ! Land-sea mask
!
REAL,    DIMENSION(:,:,:), OPTIONAL, INTENT(IN)    :: PCCS   ! Cld droplets nb conc source
REAL,    DIMENSION(:,:,:), OPTIONAL, INTENT(IN)    :: PCRS   ! Rain nb conc source
REAL,    DIMENSION(:,:,:), OPTIONAL, INTENT(IN)    :: PCSS   ! Snow nb conc source
REAL,    DIMENSION(:,:,:), OPTIONAL, INTENT(IN)    :: PCGS   ! Graupel nb conc source
REAL,    DIMENSION(:,:,:), OPTIONAL, INTENT(IN)    :: PCHS   ! Hail nb conc source
!
END SUBROUTINE ION_ATTACH_ELEC
END INTERFACE
END MODULE MODI_ION_ATTACH_ELEC
!
!
!       ####################################################################
        SUBROUTINE ION_ATTACH_ELEC(KTCOUNT, KRR, HCLOUD, PTSTEP, PRHODREF, &
                                   PRHODJ, PSVS, PRS, PTHT, PCIT, PPABST,  &
                                   PEFIELDU, PEFIELDV, PEFIELDW, GATTACH,  &
                                   PTOWN, PSEA,                            &
                                   PCCS, PCRS, PCSS, PCGS, PCHS            )
!       ####################################################################
!!
!!    PURPOSE
!!    -------
!!      This routine computes the ion capture by (or attachment to) hydrometeors
!!    providing a source of charge for hydrometeors and a sink for positive
!!    negative ion mixing ratio. It is assumed as resulting from both ionic
!!    diffusion and conduction (electrical attraction).
!!
!!
!!    METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      M. Chong     *Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    2010
!!      Modifications:
!!      J.Escobar : 18/12/2015 : Correction of bug in bound in // for NHALO <>1 
!  P. Wautelet    03/2020: use the new data structures and subroutines for budgets
!  C. Barthe      09/2022: enable the use of LIMA for cloud electrification
!  C. Barthe      09/2023: enable the use of LIMA2 for cloud electrification
!
!-------------------------------------------------------------------------------
!
!*	0.	DECLARATIONS
!		------------
!
use modd_budget,          only : lbudget_sv, NBUDGET_SV1, tbudgets
USE MODD_CONF,            ONLY : CCONF
USE MODD_CST
USE MODD_ELEC_DESCR
USE MODD_ELEC_n
USE MODD_ELEC_PARAM
USE MODD_NSV,             ONLY : NSV_ELECBEG, NSV_ELEC
USE MODD_PARAMETERS,      ONLY : JPVEXT
USE MODD_PARAM_LIMA,      ONLY : XALPHAC_L=>XALPHAC, XNUC_L=>XNUC, &
                                 XALPHAI_L=>XALPHAI, XNUI_L=>XNUI, &
                                 XCEXVT_L=>XCEXVT,                 &
                                 NMOM_S, NMOM_G, NMOM_H
USE MODD_PARAM_LIMA_COLD, ONLY : XDI, XLBI, XLBEXI, XFSEDRI,                        &
                                 XDS, XCCS_L=>XCCS, XCXS_L=>XCXS, XLBS_L=>XLBS,     &
                                 XEXSEDS_L=>XEXSEDS, XLBEXS_L=>XLBEXS, XFSEDS_L=>XFSEDS
USE MODD_PARAM_LIMA_MIXED,ONLY : XDG, XCCG_L=>XCCG, XCXG_L=>XCXG, XLBG_L=>XLBG,          &
                                 XEXSEDG_L=>XEXSEDG, XLBEXG_L=>XLBEXG, XFSEDG_L=>XFSEDG, &
                                 XDH, XALPHAH_L=>XALPHAH, XNUH_L=>XNUH,                  &
                                 XCCH_L=>XCCH, XCXH_L=>XCXH, XLBH_L=>XLBH,               &
                                 XEXSEDH_L=>XEXSEDH, XLBEXH_L=>XLBEXH, XFSEDH_L=>XFSEDH
USE MODD_PARAM_LIMA_WARM, ONLY : XCC_L=>XCC, XDC_L=>XDC, XLBC_L=>XLBC, XLBEXC_L=>XLBEXC, &
                                 XFSEDC_L=>XFSEDRC,                                      &
                                 XLBR_L=>XLBR, XLBEXR_L=>XLBEXR,                         &
                                 XFSEDR_L=>XFSEDRR, XBR, XDR
USE MODD_RAIN_ICE_DESCR_n,ONLY : XALPHAC_I=>XALPHAC, XNUC_I=>XNUC,                           &
                                 XALPHAI_I=>XALPHAI, XNUI_I=>XNUI,                           &
                                 XALPHAH_I=>XALPHAH, XNUH_I=>XNUH,                           &
                                 XCC_I=>XCC, XDC_I=>XDC, XLBC_I=>XLBC, XLBEXC_I=>XLBEXC,     &
                                 XCONC_SEA, XCONC_LAND, XCONC_URBAN, XALPHAC2, XNUC2,        &
                                 XCCR, XLBR_I=>XLBR, XLBEXR_I=>XLBEXR,                       &
                                 XCCS_I=>XCCS, XCXS_I=>XCXS, XLBS_I=>XLBS, XLBEXS_I=>XLBEXS, &
                                 XCCG_I=>XCCG, XCXG_I=>XCXG, XLBG_I=>XLBG, XLBEXG_I=>XLBEXG, &
                                 XCCH_I=>XCCH, XCXH_I=>XCXH, XLBH_I=>XLBH, XLBEXH_I=>XLBEXH, &
                                 XCEXVT_I=>XCEXVT
USE MODD_RAIN_ICE_PARAM_n,ONLY : XFSEDC_I=>XFSEDC,                     &
                                 XFSEDR_I=>XFSEDR, XEXSEDR,            &
                                 XFSEDS_I=>XFSEDS, XEXSEDS_I=>XEXSEDS, &
                                 XFSEDG_I=>XFSEDG, XEXSEDG_I=>XEXSEDG, &
                                 XFSEDH_I=>XFSEDH, XEXSEDH_I=>XEXSEDH
USE MODD_REF,             ONLY : XTHVREFZ

use mode_budget,          only : Budget_store_init, Budget_store_end
use mode_tools_ll,        only : GET_INDICE_ll

USE MODI_MOMG

IMPLICIT NONE
!
!	0.1	Declaration of arguments
!
INTEGER,                             INTENT(IN)    :: KTCOUNT  ! Temporal loop counter
INTEGER,                             INTENT(IN)    :: KRR      ! Number of moist variables
CHARACTER(LEN=4),                    INTENT(IN)    :: HCLOUD   ! kind of cloud paramerization
REAL,                                INTENT(IN)    :: PTSTEP   ! Time step          
!
REAL,    DIMENSION(:,:,:),           INTENT(IN)    :: PRHODREF ! Reference dry air density 
REAL,    DIMENSION(:,:,:),           INTENT(IN)    :: PRHODJ   ! dry air density* Jacobian
REAL,    DIMENSION(:,:,:,:),         INTENT(INOUT) :: PSVS     ! Scalar variable vol. source
REAL,    DIMENSION(:,:,:,:),         INTENT(INOUT) :: PRS      ! Moist variable vol. source
REAL,    DIMENSION(:,:,:),           INTENT(IN)    :: PTHT     ! Theta (K) at time t
REAL,    DIMENSION(:,:,:),           INTENT(IN)    :: PCIT     ! Pristine ice n.c.
                                                               ! - at time t (for ICE schemes)
                                                               ! - source (for LIMA)
REAL,    DIMENSION(:,:,:),           INTENT(IN)    :: PPABST   ! Absolute pressure at t
REAL,    DIMENSION(:,:,:),           INTENT(IN)    :: PEFIELDU, PEFIELDV, PEFIELDW
                                                               ! Electric field components
LOGICAL, DIMENSION(:,:,:),           INTENT(IN)    :: GATTACH  ! Recombination and
                                                               ! Attachment if true
REAL,    DIMENSION(:,:),   OPTIONAL, INTENT(IN)    :: PTOWN    ! Town fraction
REAL,    DIMENSION(:,:),   OPTIONAL, INTENT(IN)    :: PSEA     ! Land-sea mask
!
REAL,    DIMENSION(:,:,:), OPTIONAL, INTENT(IN)    :: PCCS   ! Cld droplets nb conc source
REAL,    DIMENSION(:,:,:), OPTIONAL, INTENT(IN)    :: PCRS   ! Rain nb conc source
REAL,    DIMENSION(:,:,:), OPTIONAL, INTENT(IN)    :: PCSS   ! Snow nb conc source
REAL,    DIMENSION(:,:,:), OPTIONAL, INTENT(IN)    :: PCGS   ! Graupel nb conc source
REAL,    DIMENSION(:,:,:), OPTIONAL, INTENT(IN)    :: PCHS   ! Hail nb conc source
!
!
!	0.2	Declaration of local variables
!
REAL,    DIMENSION(:), ALLOCATABLE :: ZT        ! Temperature (K)
REAL,    DIMENSION(:), ALLOCATABLE :: ZCONC, &  ! Number concentration
                                      ZVIT,  &  ! Fallspeed
                                      ZRADIUS   ! Radius
REAL                           :: ZCQD, ZCDIF   ! Computed coefficients
INTEGER, DIMENSION(SIZE(PTHT)) :: IGI, IGJ, IGK ! Valid grid index
INTEGER :: IVALID                               ! Nb of valid grid
INTEGER :: IIB           !  Beginning (B) and end (E) grid points
INTEGER :: IIE           !  along i axis,
INTEGER :: IJB           !        j axis,
INTEGER :: IJE           !
INTEGER :: IKB           !    and k axis
INTEGER :: IKE           !

INTEGER :: II, IJ, IK, JRR, JSV      ! Loop index for variable
INTEGER :: ITYPE         ! Hydrometeor category (2: cloud, 3: rain,
                         ! 4: ice crystal, 5: snow, 6: graupel, 7: hail)
REAL    :: ZCOMB         ! Recombination
!
! variables used to select between common parameters between ICEx and LIMA
REAL :: ZALPHAC, ZNUC, ZCC, ZDC,                   &
        ZFSEDC1, ZFSEDC2, ZLBC1, ZLBC2, ZLBEXC,    & 
        ZALPHAI, ZNUI,                             &
        ZLBR, ZLBEXR, ZFSEDR,                      &
        ZCCS, ZCXS, ZLBS, ZLBEXS, ZFSEDS, ZEXSEDS, &
        ZCCG, ZCXG, ZLBG, ZLBEXG, ZFSEDG, ZEXSEDG, &
        ZCCH, ZCXH, ZLBH, ZLBEXH, ZFSEDH, ZEXSEDH, &
        ZALPHAH, ZNUH,                             &
        ZCEXVT
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZCCT, & ! nb conc at t for cld droplets
                                       ZCRT, & !                  rain
                                       ZCIT, & !                  ice crystal
                                       ZCST, & !                  snow
                                       ZCGT, & !                  graupel
                                       ZCHT    !                  hail
!
!-------------------------------------------------------------------------------
!
if ( lbudget_sv ) then
  do jrr = 1, nsv_elec
    call Budget_store_init( tbudgets( NBUDGET_SV1 - 1 + nsv_elecbeg - 1 + jrr), 'NEUT', psvs(:, :, :, jrr) )
  end do
end if
!
!*       1.     PRELIMINARIES
!               -------------
!
! select parameters between ICEx and LIMA
!
IF (HCLOUD(1:3) == 'ICE') THEN
  !
  ZALPHAC = XALPHAC_I
  ZNUC    = XNUC_I
  ZCC     = XCC_I
  ZDC     = XDC_I
  ZFSEDC1 = XFSEDC_I(1)
  ZFSEDC2 = XFSEDC_I(2)
  ZLBC1   = XLBC_I(1)
  ZLBC2   = XLBC_I(2)
  ZLBEXC  = XLBEXC_I
  !
  ZALPHAI = XALPHAI_I
  ZNUI    = XNUI_I
  !
  ZLBR    = XLBR_I
  ZLBEXR  = XLBEXR_I
  ZFSEDR  = XFSEDR_I
  !
  ZCCS    = XCCS_I
  ZCXS    = XCXS_I
  ZLBS    = XLBS_I
  ZLBEXS  = XLBEXS_I
  ZFSEDS  = XFSEDS_I
  ZEXSEDS = XEXSEDS_I
  !
  ZCCG    = XCCG_I
  ZCXG    = XCXG_I
  ZLBG    = XLBG_I
  ZLBEXG  = XLBEXG_I
  ZFSEDG  = XFSEDG_I
  ZEXSEDG = XEXSEDG_I
  !
  ZALPHAH = XALPHAH_I
  ZNUH    = XNUH_I
  ZCCH    = XCCH_I
  ZCXH    = XCXH_I
  ZLBH    = XLBH_I
  ZLBEXH  = XLBEXH_I
  ZFSEDH  = XFSEDH_I
  ZEXSEDH = XEXSEDH_I
  !
  ZCEXVT = XCEXVT_I
  !
ELSE IF (HCLOUD == 'LIMA') THEN
  !
  ZALPHAC = XALPHAC_L
  ZNUC    = XNUC_L
  ZCC     = XCC_L
  ZDC     = XDC_L
  ZFSEDC1 = XFSEDC_L
  ZLBC1   = XLBC_L
  ZLBEXC  = XLBEXC_L
  !
  ZALPHAI = XALPHAI_L
  ZNUI    = XNUI_L
  !
  ZLBR    = XLBR_L
  ZLBEXR  = XLBEXR_L
  ZFSEDR  = XFSEDR_L
  !
  ZCCS    = XCCS_L
  ZCXS    = XCXS_L
  ZLBS    = XLBS_L
  ZLBEXS  = XLBEXS_L
  ZFSEDS  = XFSEDS_L
  ZEXSEDS = XEXSEDS_L
  !
  ZCCG    = XCCG_L
  ZCXG    = XCXG_L
  ZLBG    = XLBG_L
  ZLBEXG  = XLBEXG_L
  ZFSEDG  = XFSEDG_L
  ZEXSEDG = XEXSEDG_L
  !
  ZALPHAH = XALPHAH_L
  ZNUH    = XNUH_L
  ZCCH    = XCCH_L
  ZCXH    = XCXH_L
  ZLBH    = XLBH_L
  ZLBEXH  = XLBEXH_L
  ZFSEDH  = XFSEDH_L
  ZEXSEDH = XEXSEDH_L
  !
  ZCEXVT = XCEXVT_L
END IF
!
!-------------------------------------------------------------------------------
!
!*       2.     COMPUTE THE ION RECOMBINATION and TEMPERATURE
!               ---------------------------------------------
!
ZCQD = 4 * XPI * XEPSILON * XBOLTZ / XECHARGE
ZCDIF = XBOLTZ / XECHARGE
!
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
IKB = 1 + JPVEXT
IKE = SIZE(PTHT,3) - JPVEXT
!
!
!*       2.1    Add Ion Recombination source (PSVS in 1/(m3.s))
!               and count and localize valid grid points for ion source terms
!
IVALID = 0
DO IK = IKB, IKE
  DO IJ = IJB, IJE
    DO II = IIB, IIE
      IF (GATTACH(II,IJ,IK)) THEN
! Recombination rate
        ZCOMB = XIONCOMB * (PSVS(II,IJ,IK,1) * PTSTEP) *        &
                           (PSVS(II,IJ,IK,NSV_ELEC) * PTSTEP) * &
                           PRHODREF(II,IJ,IK) / PRHODJ(II,IJ,IK)
        ZCOMB = MIN(ZCOMB, PSVS(II,IJ,IK,1), PSVS(II,IJ,IK,NSV_ELEC))
        !
! Update the sources
        PSVS(II,IJ,IK,1)        = PSVS(II,IJ,IK,1)        - ZCOMB
        PSVS(II,IJ,IK,NSV_ELEC) = PSVS(II,IJ,IK,NSV_ELEC) - ZCOMB
        !
! Counting
        IVALID = IVALID + 1
        IGI(IVALID) = II
        IGJ(IVALID) = IJ
        IGK(IVALID) = IK
      END IF
    ENDDO
  ENDDO
ENDDO
!
!
!*       2.2    Compute the temperature
!
IF( IVALID /= 0 ) THEN
  ALLOCATE (ZT(IVALID))
  DO II = 1, IVALID
    ZT(II) = PTHT(IGI(II),IGJ(II),IGK(II)) * &
            (PPABST(IGI(II),IGJ(II),IGK(II)) / XP00) ** (XRD / XCPD)
  ENDDO
END IF
!
!-------------------------------------------------------------------------------
!
!*	 3.     TRANSFORM VOLUM. SOURCE TERMS INTO MIXING RATIO
!               FOR WATER SPECIES, AND VOLUMIC CONTENT FOR ELECTRIC VARIABLES
!               -------------------------------------------------------------
!
DO JRR = 1, KRR
  PRS(:,:,:,JRR) = PRS(:,:,:,JRR) * PTSTEP / PRHODJ(:,:,:)
ENDDO
!
ALLOCATE(ZCIT(SIZE(PTHT,1),SIZE(PTHT,2),SIZE(PTHT,3)))
! ICEx : pcit is really pcit
IF (HCLOUD(1:3) == 'ICE') ZCIT(:,:,:) = PCIT(:,:,:)
!
IF (HCLOUD == 'LIMA') THEN
  ALLOCATE(ZCCT(SIZE(PTHT,1),SIZE(PTHT,2),SIZE(PTHT,3)))
  ALLOCATE(ZCRT(SIZE(PTHT,1),SIZE(PTHT,2),SIZE(PTHT,3)))
  !
  ZCCT(:,:,:) = PCCS(:,:,:) * PTSTEP / PRHODJ(:,:,:)
  ZCRT(:,:,:) = PCRS(:,:,:) * PTSTEP / PRHODJ(:,:,:)
! LIMA : pcit is pcis !
  ZCIT(:,:,:) = PCIT(:,:,:) * PTSTEP / PRHODJ(:,:,:)
  !
  IF (PRESENT(PCSS)) THEN
    ALLOCATE(ZCST(SIZE(PTHT,1),SIZE(PTHT,2),SIZE(PTHT,3)))
    ZCST(:,:,:) = PCSS(:,:,:) * PTSTEP / PRHODJ(:,:,:)
  END IF
  IF (PRESENT(PCGS)) THEN
    ALLOCATE(ZCGT(SIZE(PTHT,1),SIZE(PTHT,2),SIZE(PTHT,3)))
    ZCGT(:,:,:) = PCGS(:,:,:) * PTSTEP / PRHODJ(:,:,:)
  END IF
  IF (PRESENT(PCHS)) THEN
    ALLOCATE(ZCHT(SIZE(PTHT,1),SIZE(PTHT,2),SIZE(PTHT,3)))
    ZCHT(:,:,:) = PCHS(:,:,:) * PTSTEP / PRHODJ(:,:,:)
  END IF
END IF
!
DO JSV = 1, NSV_ELEC
  PSVS(:,:,:,JSV) = PSVS(:,:,:,JSV) * PTSTEP * PRHODREF(:,:,:) / PRHODJ(:,:,:)
ENDDO
!
!-------------------------------------------------------------------------------
!
!*	 4.	COMPUTE ATTACHMENT DUE TO ION DIFFUSION AND CONDUCTION
!               ------------------------------------------------------
! 
!  Attachment to cloud droplets, rain, cloud ice, snow, graupel, 
!                and hail (optional)
!
!
IF( IVALID /= 0 ) THEN
!
!
!*	4.1	Attachment to cloud droplets
!
  ALLOCATE (ZCONC(IVALID))
  ALLOCATE (ZVIT (IVALID))
  ALLOCATE (ZRADIUS(IVALID))

  ITYPE = 2
  IF (PRESENT(PSEA)) THEN
    CALL HYDROPARAM (IGI, IGJ, IGK, ZCONC, ZVIT, ZRADIUS, ITYPE, PSEA, PTOWN)
  ELSE
    CALL HYDROPARAM (IGI, IGJ, IGK, ZCONC, ZVIT, ZRADIUS, ITYPE)
  ENDIF
!
  CALL DIFF_COND (IGI, IGJ, IGK, PSVS(:,:,:,1), PSVS(:,:,:,NSV_ELEC),  &
                  PSVS(:,:,:,ITYPE))
!
!
!*	4.2	Attachment to raindrops, ice crystals, snow, graupel, 
!                             and hail (if activated)
!
  DO ITYPE = 3, KRR
    CALL HYDROPARAM (IGI, IGJ, IGK, ZCONC, ZVIT, ZRADIUS, ITYPE)
!
    CALL DIFF_COND (IGI, IGJ, IGK, PSVS(:,:,:,1), PSVS(:,:,:,NSV_ELEC),  &
                    PSVS(:,:,:,ITYPE))
  END DO
!
  DEALLOCATE (ZCONC, ZVIT, ZRADIUS)
  DEALLOCATE (ZT)
ENDIF
!
IF (ALLOCATED(ZCCT)) DEALLOCATE(ZCCT)
IF (ALLOCATED(ZCRT)) DEALLOCATE(ZCRT)
IF (ALLOCATED(ZCIT)) DEALLOCATE(ZCIT)
IF (ALLOCATED(ZCST)) DEALLOCATE(ZCST)
IF (ALLOCATED(ZCGT)) DEALLOCATE(ZCGT)
IF (ALLOCATED(ZCHT)) DEALLOCATE(ZCHT)
!
!-------------------------------------------------------------------------------
!
!*	5.	RETURN TO VOLUMETRIC SOURCE (Prognostic units)
!               ---------------------------
!
DO JRR = 1, KRR
  PRS(:,:,:,JRR) = PRS(:,:,:,JRR) * PRHODJ(:,:,:) / PTSTEP
ENDDO
!
DO JSV = 1, NSV_ELEC
  PSVS(:,:,:,JSV) = PSVS(:,:,:,JSV) * PRHODJ(:,:,:) / (PTSTEP * PRHODREF(:,:,:))
ENDDO
!
!-------------------------------------------------------------------------------
!
!*	6.	BUDGET
!               ------
!
if ( lbudget_sv ) then
  do jrr = 1, nsv_elec
    call Budget_store_end( tbudgets( NBUDGET_SV1 - 1 + nsv_elecbeg - 1 + jrr), 'NEUT', psvs(:, :, :, jrr) )
  end do
end if
!
!------------------------------------------------------------------------------
!
CONTAINS
!
!------------------------------------------------------------------------------
!
 SUBROUTINE HYDROPARAM (IGRIDX, IGRIDY, IGRIDZ, ZCONC, &
                        ZVIT, ZRADIUS, ITYPE, PSEA, PTOWN)
!
!       Purpose : Compute in regions of valid grid points (IGRIDX, IGRIDY, IGRIDZ)
!                 the hydrometeor parameters: concentration (ZCONC),
!                                             fallspeed (ZVIT),
!                                         and mean radius (ZRADIUS)
!                 involved in the evaluation of ion attachment
!
!
!*      0.      DECLARATIONS
!               ------------
IMPLICIT NONE
!
!*      0.1     declaration of dummy arguments
!
INTEGER, DIMENSION(:), INTENT(IN)    :: IGRIDX, &  !  Index of
                                        IGRIDY, &  !   valid
                                        IGRIDZ     ! gridpoints
INTEGER,               INTENT(IN)    :: ITYPE      ! Hydrometeor category
              ! ITYPE= 2: cloud, 3: rain, 4: ice, 5: snow, 6: graupel, 7: hail
REAL,    DIMENSION(:), INTENT(INOUT) :: ZCONC, &   ! Number concentration
                                        ZVIT,  &   ! Fallspeed
                                        ZRADIUS    ! Radius
REAL,    DIMENSION(:,:), OPTIONAL, INTENT(IN) :: PTOWN ! Town fraction
REAL,    DIMENSION(:,:), OPTIONAL, INTENT(IN) :: PSEA  ! Land-sea mask
!
!*      0.2     declaration of local variables
!
REAL    :: ZCONC1, ZCONC2 ! for cloud
REAL    :: ZLBC
REAL    :: ZFSEDC
REAL    :: ZRAY 
REAL    :: ZEXP1, ZEXP2, ZMOM1, ZMOM2
REAL    :: ZVCOEF, ZRHO00, ZLBI
REAL    :: ZLAMBDA
REAL    :: ZCOR    ! correction factor for cloud droplet terminal fall speed
INTEGER :: JI, JJ, JK, IV
!
!
ZCONC(:) = 0.
ZVIT (:) = 0.
ZRADIUS(:) = 0.
!
SELECT CASE (ITYPE)
!
!*	1.	PARAMETERS FOR CLOUD
!               --------------------
  CASE (2)
!
    IF (HCLOUD(1:3) == 'ICE') THEN
      IF (PRESENT(PSEA)) THEN
        ZMOM1 = 0.5 * MOMG(ZALPHAC,ZNUC,1.)
        ZMOM2 = 0.5 * MOMG(XALPHAC2,XNUC2,1.)      
        DO IV = 1, IVALID
          JI = IGRIDX(IV)
          JJ = IGRIDY(IV)
          JK = IGRIDZ(IV)
          IF( PRS(JI,JJ,JK,2)/PRHODREF(JI,JJ,JK) > XRTMIN_ELEC(2) .AND. &
              PSVS(JI,JJ,JK,2) /= 0. ) THEN
            ZCONC1 = PSEA(JI,JJ) * XCONC_SEA + (1. - PSEA(JI,JJ)) * XCONC_LAND
            ZLBC   = PSEA(JI,JJ) * ZLBC2     + (1. - PSEA(JI,JJ)) * ZLBC1
            ZFSEDC = PSEA(JI,JJ) * ZFSEDC2   + (1. - PSEA(JI,JJ)) * ZFSEDC1
            ZFSEDC = MAX(MIN(ZFSEDC1,ZFSEDC2), ZFSEDC)
            ZCONC2 = (1. - PTOWN(JI,JJ)) * ZCONC1 + PTOWN(JI,JJ) * XCONC_URBAN
            ZRAY   = (1. - PSEA(JI,JJ))  * ZMOM1  + PSEA(JI,JJ) * ZMOM2
            ZCONC(IV) = ZCONC2                 ! Number concentration
            ZLAMBDA   = (ZLBC * ZCONC2 / (PRHODREF(JI,JJ,JK) * PRS(JI,JJ,JK,2)))**ZLBEXC
            ZRADIUS(IV) = ZRAY / ZLAMBDA
            ZVIT(IV)    = ZCC * ZFSEDC * ZLAMBDA**(-ZDC) * PRHODREF(JI,JJ,JK)**(-ZCEXVT)
          END IF
        ENDDO
      ELSE
        ZRAY = 0.5 * MOMG(ZALPHAC,ZNUC,1.)
        ZLBC = ZLBC1 * XCONC_LAND
        DO IV = 1, IVALID
          JI = IGRIDX(IV)
          JJ = IGRIDY(IV)
          JK = IGRIDZ(IV)
          IF( PRS(JI,JJ,JK,2)/PRHODREF(JI,JJ,JK) > XRTMIN_ELEC(2) .AND. &
              PSVS(JI,JJ,JK,2) /= 0. ) THEN
            ZCONC(IV) = XCONC_LAND             ! Number concentration
            ZLAMBDA   = (ZLBC / (PRHODREF(JI,JJ,JK) * PRS(JI,JJ,JK,2)))**ZLBEXC
            ZRADIUS(IV) = ZRAY / ZLAMBDA
            ZVIT(IV)    = ZCC * ZFSEDC1 * ZLAMBDA**(-ZDC) * PRHODREF(JI,JJ,JK)**(-ZCEXVT)       
          END IF
        ENDDO
      END IF
    ELSE IF (HCLOUD == 'LIMA') THEN
      ZRAY = 0.5 * MOMG(ZALPHAC,ZNUC,1.)
      DO IV = 1, IVALID
        JI = IGRIDX(IV)
        JJ = IGRIDY(IV)
        JK = IGRIDZ(IV)
        IF( PRS(JI,JJ,JK,2)/PRHODREF(JI,JJ,JK) > XRTMIN_ELEC(2) .AND. &
            ZCCT(JI,JJ,JK) > 0. .AND. PSVS(JI,JJ,JK,2) /= 0. ) THEN
          ZCONC(IV)   = ZCCT(JI,JJ,JK) * PRHODREF(JI,JJ,JK)  ! Number concentration (m-3)
          ZLAMBDA     = (ZLBC1 * ZCCT(JI,JJ,JK) / PRS(JI,JJ,JK,2))**ZLBEXC
          ZRADIUS(IV) = ZRAY / ZLAMBDA
! correction factor for cloud droplet terminal fall speed
          ZCOR        = 1. + 1.26 * 6.6E-8 * (101325. / PPABST(JI,JJ,JK)) * (ZT(IV) / 293.15) / ZRADIUS(IV)
          ZVIT(IV)    = ZCOR * ZFSEDC1 * ZLAMBDA**(-ZDC) * PRHODREF(JI,JJ,JK)**(-ZCEXVT)
        END IF
      ENDDO
    END IF
!
!
!*	2.	PARAMETERS FOR RAIN
!               -------------------
  CASE (3)
!
    IF (HCLOUD(1:3) == 'ICE') THEN
      ZEXP1 = XEXSEDR - 1.
      ZEXP2 = ZEXP1 - ZCEXVT
!
      DO IV = 1, IVALID
        JI = IGRIDX(IV)
        JJ = IGRIDY(IV)
        JK = IGRIDZ(IV)
        IF( PRS(JI,JJ,JK,3)/PRHODREF(JI,JJ,JK) > XRTMIN_ELEC(3) .AND. &
            PSVS(JI,JJ,JK,3) /= 0. ) THEN
          ZLAMBDA = ZLBR * (PRHODREF(JI,JJ,JK) * PRS(JI,JJ,JK,3))**ZLBEXR
! dans ice3, alpha_r = 1 et nu_r = 1
          ZRADIUS(IV) = 0.5 / ZLAMBDA
          ZCONC(IV)   = XCCR / ZLAMBDA
          ZVIT(IV)    = ZFSEDR * PRHODREF(JI,JJ,JK)**ZEXP2 &
                                  * PRS(JI,JJ,JK,3)**ZEXP1
        END IF
      ENDDO
    ELSE IF (HCLOUD == 'LIMA') THEN
      DO IV = 1, IVALID
        JI = IGRIDX(IV)
        JJ = IGRIDY(IV)
        JK = IGRIDZ(IV)
        IF( PRS(JI,JJ,JK,3)/PRHODREF(JI,JJ,JK) > XRTMIN_ELEC(3) .AND. &
            ZCRT(JI,JJ,JK) > 0. .AND. PSVS(JI,JJ,JK,3) /= 0. ) THEN
          ZLAMBDA     = (ZLBR * ZCRT(JI,JJ,JK) / PRS(JI,JJ,JK,3))**ZLBEXR
! dans lima, alpha_r = 1 et nu_r = 2
          ZRADIUS(IV) = 1. / ZLAMBDA
          ZCONC(IV)   = ZCRT(JI,JJ,JK) * PRHODREF(JI,JJ,JK)  ! Number concentration (m-3)
          ! zvit = zwsedr / (r * rho_dref)
          ZVIT(IV)    = ZFSEDR * ZLAMBDA**(-XDR) * PRHODREF(JI,JJ,JK)**(-ZCEXVT)
        END IF
      ENDDO
    END IF
!
!
!*	3.	PARAMETERS FOR ICE
!               ------------------
!
  CASE (4)
!
    ZRAY = 0.5 * MOMG(ZALPHAI,ZNUI,1.)
    !
    IF (HCLOUD(1:3) == 'ICE') THEN
      ZRHO00 = XP00 / (XRD * XTHVREFZ(IKB))
! Computations for Columns (see ini_rain_ice_elec.f90)
      ZVCOEF = 2.1E5 * MOMG(ZALPHAI,ZNUI, 3.285) / MOMG(ZALPHAI,ZNUI, 1.7)   &
                     * ZRHO00**ZCEXVT
      ZLBI = (2.14E-3 * MOMG(ZALPHAI,ZNUI,1.7)) **0.588235
      !
      DO IV = 1, IVALID
        JI = IGRIDX(IV)
        JJ = IGRIDY(IV)
        JK = IGRIDZ(IV)
        IF( PRS(JI, JJ, JK, 4)/PRHODREF(JI, JJ, JK) > XRTMIN_ELEC(4) .AND. &
            PSVS(JI, JJ, JK, 4) /=0.) THEN
          ZCONC(IV)   = XFCI * PRHODREF(JI,JJ,JK) * PRS(JI,JJ,JK,4) * &
                        MAX(0.05E6, -0.15319E6 - 0.021454E6 *         &
                            ALOG(PRHODREF(JI,JJ,JK) * PRS(JI,JJ,JK,4)))**3
          ZLAMBDA     = ZLBI * (ZCONC(IV) / (PRHODREF(JI,JJ,JK) * &
                                PRS(JI,JJ,JK,4)))**0.588235
          ZRADIUS(IV) = ZRAY / ZLAMBDA
          ZVIT(IV)    = ZVCOEF * ZLAMBDA**(-1.585) * &      !(-XDI) * &
                        PRHODREF(JI,JJ,JK)**(-ZCEXVT) 
        END IF
      ENDDO
    ELSE IF (HCLOUD == 'LIMA') THEN
      DO IV = 1, IVALID
        JI = IGRIDX(IV)
        JJ = IGRIDY(IV)
        JK = IGRIDZ(IV)
        IF( PRS(JI,JJ,JK,4)/PRHODREF(JI,JJ,JK) > XRTMIN_ELEC(4) .AND. &
            ZCIT(JI,JJ,JK) > 0. .AND. PSVS(JI,JJ,JK,4) /= 0.) THEN
          ZCONC(IV)   = ZCIT(JI,JJ,JK) * PRHODREF(JI,JJ,JK)  ! Number concentration (m-3)
          ZLAMBDA     = (XLBI * ZCIT(JI,JJ,JK) / PRS(JI,JJ,JK,4))**XLBEXI
          ZRADIUS(IV) = ZRAY / ZLAMBDA
          ! zvit = zwsedr / (r * rho_dref)
          ZVIT(IV)    = XFSEDRI * ZLAMBDA**(-XDI) * PRHODREF(JI,JJ,JK)**(-ZCEXVT) 
        END IF
      ENDDO
    END IF
!
!
!*	4.	PARAMETERS FOR SNOW 
!               -------------------
!
  CASE (5)
!
    IF ((HCLOUD(1:3) == 'ICE') .OR. (HCLOUD == 'LIMA' .AND. NMOM_S == 1)) THEN
      ZEXP1 = ZEXSEDS - 1.
      ZEXP2 = ZEXP1 - ZCEXVT
!
      DO IV = 1, IVALID
        JI = IGRIDX(IV)
        JJ = IGRIDY(IV)
        JK = IGRIDZ(IV)
        IF( PRS(JI,JJ,JK,5)/PRHODREF(JI,JJ,JK) > XRTMIN_ELEC(5) .AND. &
            PSVS(JI,JJ,JK,5) /= 0. ) THEN
          ZLAMBDA = ZLBS * (PRHODREF(JI,JJ,JK) * PRS(JI,JJ,JK,5))**ZLBEXS
          ZRADIUS(IV) = 0.5 / ZLAMBDA
          ZCONC(IV)   = ZCCS * ZLAMBDA**ZCXS
          ZVIT(IV)    = ZFSEDS * PRHODREF(JI,JJ,JK)**ZEXP2 &
                                  * PRS(JI,JJ,JK,5)**ZEXP1
        END IF
      ENDDO
    ELSE IF (HCLOUD == 'LIMA' .AND. NMOM_S == 2) THEN
      DO IV = 1, IVALID
        JI = IGRIDX(IV)
        JJ = IGRIDY(IV)
        JK = IGRIDZ(IV)
        IF( PRS(JI,JJ,JK,5)/PRHODREF(JI,JJ,JK) > XRTMIN_ELEC(5) .AND. &
            ZCST(JI,JJ,JK) > 0. .AND. PSVS(JI,JJ,JK,5) /= 0. ) THEN
          ZLAMBDA     = (ZLBS * ZCST(JI,JJ,JK) / PRS(JI,JJ,JK,5))**ZLBEXS
          ZRADIUS(IV) = 1. / ZLAMBDA
          ZCONC(IV)   = ZCST(JI,JJ,JK) * PRHODREF(JI,JJ,JK)  ! Number concentration (m-3)
          ZVIT(IV)    = ZFSEDS * ZLAMBDA**(-XDS) * PRHODREF(JI,JJ,JK)**(-ZCEXVT)
        END IF
      ENDDO
    END IF
!
!
!*	5.	PARAMETERS FOR GRAUPEL
!               ----------------------
!
  CASE (6)
!
    IF ((HCLOUD(1:3) == 'ICE') .OR. (HCLOUD == 'LIMA' .AND. NMOM_G == 1)) THEN
      ZEXP1 = ZEXSEDG - 1.
      ZEXP2 = ZEXP1 - ZCEXVT
!
      DO IV = 1, IVALID
        JI = IGRIDX(IV)
        JJ = IGRIDY(IV)
        JK = IGRIDZ(IV)
        IF( PRS(JI,JJ,JK,6)/PRHODREF(JI,JJ,JK) > XRTMIN_ELEC(6) .AND. &
            PSVS(JI,JJ,JK,6) /= 0. ) THEN
          ZLAMBDA = ZLBG * (PRHODREF(JI,JJ,JK) * PRS(JI,JJ,JK,6))**ZLBEXG
          ZRADIUS(IV) = 0.5 / ZLAMBDA
          ZCONC(IV)   = ZCCG * ZLAMBDA**ZCXG
          ZVIT(IV)    = ZFSEDG * PRHODREF(JI,JJ,JK)**ZEXP2 &
                                  * PRS(JI,JJ,JK,6)**ZEXP1
        END IF
      ENDDO
    ELSE IF (HCLOUD == 'LIMA' .AND. NMOM_G == 2) THEN
      DO IV = 1, IVALID
        JI = IGRIDX(IV)
        JJ = IGRIDY(IV)
        JK = IGRIDZ(IV)
        IF( PRS(JI,JJ,JK,6)/PRHODREF(JI,JJ,JK) > XRTMIN_ELEC(6) .AND. &
            ZCGT(JI,JJ,JK) > 0. .AND. PSVS(JI,JJ,JK,6) /= 0. ) THEN
          ZLAMBDA     = (ZLBG * ZCGT(JI,JJ,JK) / PRS(JI,JJ,JK,6))**ZLBEXG
          ZRADIUS(IV) = 1. / ZLAMBDA
          ZCONC(IV)   = ZCGT(JI,JJ,JK) * PRHODREF(JI,JJ,JK)  ! Number concentration (m-3)
          ZVIT(IV)    = ZFSEDG * ZLAMBDA**(-XDG) * PRHODREF(JI,JJ,JK)**(-ZCEXVT)
        END IF
      ENDDO
    END IF
!
!
!*	6.	PARAMETERS FOR HAIL
!               -------------------
!
  CASE (7)
!
    IF ((HCLOUD(1:3) == 'ICE') .OR. (HCLOUD == 'LIMA' .AND. NMOM_H == 1)) THEN
      ZEXP1 = ZEXSEDH - 1.
      ZEXP2 = ZEXP1 - ZCEXVT
!
      DO IV = 1, IVALID
        JI = IGRIDX(IV)
        JJ = IGRIDY(IV)
        JK = IGRIDZ(IV)
        IF( PRS(JI,JJ,JK,7)/PRHODREF(JI,JJ,JK) > XRTMIN_ELEC(7) .AND. &
            PSVS(JI,JJ,JK,7) /= 0. ) THEN
          ZLAMBDA = ZLBH * (PRHODREF(JI,JJ,JK) * PRS(JI,JJ,JK,7))**ZLBEXH
          ZRADIUS(IV) = 0.5 / ZLAMBDA
          ZCONC(IV)   = ZCCH * ZLAMBDA**ZCXH
          ZVIT(IV)    = ZFSEDH * PRHODREF(JI,JJ,JK)**ZEXP2 &
                                  * PRS(JI,JJ,JK,7)**ZEXP1
        END IF
      ENDDO
    ELSE IF (HCLOUD == 'LIMA' .AND. NMOM_H == 2) THEN
      DO IV = 1, IVALID
        JI = IGRIDX(IV)
        JJ = IGRIDY(IV)
        JK = IGRIDZ(IV)
        IF( PRS(JI,JJ,JK,7)/PRHODREF(JI,JJ,JK) > XRTMIN_ELEC(7) .AND. &
            ZCHT(JI,JJ,JK) > 0. .AND. PSVS(JI,JJ,JK,7) /= 0. ) THEN
          ZLAMBDA     = (ZLBH * ZCHT(JI,JJ,JK) / PRS(JI,JJ,JK,7))**ZLBEXH
          ZRADIUS(IV) = 1. / ZLAMBDA
          ZCONC(IV)   = ZCHT(JI,JJ,JK) * PRHODREF(JI,JJ,JK)  ! Number concentration (m-3)
          ZVIT(IV)    = ZFSEDH * ZLAMBDA**(-XDH) * PRHODREF(JI,JJ,JK)**(-ZCEXVT)
        END IF
      ENDDO
    END IF
!
END SELECT
!
END SUBROUTINE HYDROPARAM
!
!------------------------------------------------------------------------------
!
 SUBROUTINE DIFF_COND (IGRIDX, IGRIDY, IGRIDZ, PQPIS, PQNIS, PQVS)
!
!       Purpose : Compute in regions of valid grid points (IGRIDX, IGRIDY, IGRIDZ)
!                 the attachment of positive (sink for PQPIS) and negative 
!                 (sink for PQNIS) ions to the hydrometeor variable (charge
!                 source for PQVS)
!
!
!*      0.	DECLARATIONS
!        	------------
IMPLICIT NONE
!
!*      0.1	declaration of dummy arguments
!
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PQPIS     ! Positive ion concentration
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PQNIS     ! Negative ion concentration
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PQVS      !Hydrom volumetric charge
INTEGER, DIMENSION(:), INTENT(IN)     :: IGRIDX, IGRIDY, IGRIDZ ! Index of 
                                                   ! valid gridpoints

!
!*      0.2	declaration of local variables
!
INTEGER :: JI, JJ, JK, IV
REAL    :: ZNC, ZRADI, ZVT ! Nb conc., radius, fallspeed of the hydrometeor category
REAL    :: ZQ              ! net particule charge
REAL    :: ZX, ZFXP, ZFXN  ! Limiting diffusion function ZFX = +/- ZX /(exp(+/-ZX) -1)
REAL    :: ZDIFP, ZDPIDT_D ! Diffusion of positive ions
REAL    :: ZDIFM, ZDNIDT_D ! Diffusion of negative ions
REAL    :: ZDPIDT_C        ! Conduction of positive ions
REAL    :: ZDNIDT_C        ! Conduction of negative ions
REAL    :: ZDELPI, ZDELNI  ! Total attachment of pos/neg ions
REAL    :: ZEFIELD         ! Electric field magnitude
REAL    :: ZQBOUND         ! Limit charge for conduction
!
!
!*	1.	COMPUTE ION ATTACHMENT
!               ----------------------
!
DO IV = 1, IVALID
  IF (ZCONC(IV) .NE. 0.) THEN
    JI = IGRIDX(IV)
    JJ = IGRIDY(IV)
    JK = IGRIDZ(IV)
!
    ZNC   = ZCONC(IV)
    ZRADI = ZRADIUS(IV)
    ZVT   = ZVIT(IV)
!
!*	1.0	Ion diffusion to a particle
!
    ZDPIDT_D = 0.
    ZDNIDT_D = 0.
!
    ZQ = PQVS(JI,JJ,JK) / ZNC
    ZX = ZQ / (ZCQD * ZRADI * ZT(IV))
!
    IF (ZX /= 0. .AND. ABS(ZX) <= 20.0) THEN
      IF (ABS(ZX) < 1.0E-15) THEN
        ZFXP = 1.
        ZFXN = 1.
      ELSE
        ZFXP =  ZX / (EXP(ZX) - 1.)
        ZFXN = -ZX / (EXP(-ZX) -1.)
      ENDIF
!  
      ZDIFP = 4. * XPI * XMOBIL_POS(JI,JJ,JK) * ZCDIF * ZT(IV)
      ZDPIDT_D = ZRADI * ZDIFP * PQPIS(JI,JJ,JK) * ZFXP * &
                 (1. + (2. * ZRADI * ZVT / ZDIFP)**0.5)
!    
      ZDIFM = 4. * XPI * XMOBIL_NEG(JI,JJ,JK) * ZCDIF * ZT(IV)
      ZDNIDT_D = ZRADI * ZDIFM * PQNIS(JI,JJ,JK) * ZFXN * &
                 (1. + (2. * ZRADI * ZVT / ZDIFM)**0.5)
!
      ZDELPI = MIN(ZDPIDT_D*PTSTEP*ZNC, PQPIS(JI,JJ,JK))
      ZDELNI = MIN(ZDNIDT_D*PTSTEP*ZNC, PQNIS(JI,JJ,JK))
!
      PQPIS(JI,JJ,JK) = PQPIS(JI,JJ,JK) - ZDELPI
      PQNIS(JI,JJ,JK) = PQNIS(JI,JJ,JK) - ZDELNI
      PQVS(JI,JJ,JK)  = PQVS(JI,JJ,JK)  + XECHARGE * (ZDELPI - ZDELNI)
    ENDIF
!
!
!*	1.1	Ion conduction to a particle
!
     ZDPIDT_C = 0.
     ZDNIDT_C = 0.
     ZEFIELD = SQRT(PEFIELDU(JI,JJ,JK)**2+PEFIELDV(JI,JJ,JK)**2+ &
                    PEFIELDW(JI,JJ,JK)**2)
     ZQBOUND = 12. * XPI * XEPSILON * ZEFIELD * ZRADI**2
     ZQ = PQVS(JI,JJ,JK) / ZNC
!
     IF (ABS(ZQ) < ZQBOUND) THEN
       IF (PEFIELDW(JI,JJ,JK) > 0.) THEN  ! opposite to fall velocity direction
         ZDPIDT_C = 3. * XPI * ZRADI**2 * ZEFIELD * PQPIS(JI,JJ,JK) * &
                    XMOBIL_POS(JI,JJ,JK) * (1. - ZQ / ZQBOUND)**2
         IF (ZVT < XMOBIL_NEG(JI,JJ,JK)*ZEFIELD) THEN
           ZDNIDT_C = 3. * XPI * ZRADI**2 * ZEFIELD * PQNIS(JI,JJ,JK) * &
                      XMOBIL_NEG(JI,JJ,JK) * (1. + ZQ / ZQBOUND)**2
         ELSE IF (ZQ > 0.) THEN
           ZDNIDT_C = PQNIS(JI,JJ,JK) * XMOBIL_NEG(JI,JJ,JK) * ZQ / XEPSILON
         ENDIF
       ELSE IF (PEFIELDW(JI,JJ,JK) < 0.) THEN  ! in the direction of fall veloc.
         IF( ZVT < XMOBIL_POS(JI,JJ,JK)*ZEFIELD) THEN
           ZDPIDT_C = 3. * XPI * ZRADI**2 * ZEFIELD * PQPIS(JI,JJ,JK) * &
                      XMOBIL_POS(JI,JJ,JK) * (1. - ZQ / ZQBOUND)**2
         ELSE IF (ZQ < 0.) THEN
           ZDPIDT_C = -PQPIS(JI,JJ,JK) * XMOBIL_POS(JI,JJ,JK) * ZQ / XEPSILON
         ENDIF
         ZDNIDT_C = 3. * XPI * ZRADI**2 * ZEFIELD * PQNIS(JI,JJ,JK) * &
                    XMOBIL_NEG(JI,JJ,JK) * (1. + ZQ / ZQBOUND)**2
       ENDIF
     ELSE IF (ZQ >= ZQBOUND) THEN
       ZDPIDT_C = 0.
       ZDNIDT_C = PQNIS(JI,JJ,JK) * XMOBIL_NEG(JI,JJ,JK) * ZQ / XEPSILON
     ELSE IF (ZQ <= -ZQBOUND) THEN
       ZDPIDT_C = -PQPIS(JI,JJ,JK) * XMOBIL_POS(JI,JJ,JK) * ZQ / XEPSILON
       ZDNIDT_C = 0.
     ENDIF
!
     ZDELPI = MIN(ZDPIDT_C*PTSTEP*ZNC, PQPIS(JI,JJ,JK))
     ZDELNI = MIN(ZDNIDT_C*PTSTEP*ZNC, PQNIS(JI,JJ,JK))
!
     PQPIS(JI,JJ,JK) = PQPIS(JI,JJ,JK) - ZDELPI
     PQNIS(JI,JJ,JK) = PQNIS(JI,JJ,JK) - ZDELNI
     PQVS(JI,JJ,JK) = PQVS(JI,JJ,JK) + XECHARGE * (ZDELPI - ZDELNI)
  END IF
ENDDO
!
END SUBROUTINE DIFF_COND
!
!-----------------------------------------------------------------------------
!
END SUBROUTINE ION_ATTACH_ELEC
