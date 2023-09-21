!MNH_LIC Copyright 2010-2020 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ############################
      MODULE MODI_ION_ATTACH_ELEC
!     ############################
!
INTERFACE
        SUBROUTINE ION_ATTACH_ELEC(KTCOUNT, KRR, PTSTEP, PRHODREF,           &
                            PRHODJ,PSVS, PRS, PTHT, PCIT, PPABST, PEFIELDU,  &
                            PEFIELDV, PEFIELDW, GATTACH, PTOWN, PSEA      )


INTEGER,                  INTENT(IN)    :: KTCOUNT  ! Temporal loop counter
INTEGER,                  INTENT(IN)    :: KRR      ! Number of moist variables
REAL,                     INTENT(IN)    :: PTSTEP   ! Time step
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF ! Reference dry air density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ   ! dry air density* Jacobian
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PSVS   ! Scalar variable vol. source
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRS    ! Moist variable vol. source
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT     ! Theta (K) at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCIT    ! Pristine ice n.c. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABST   ! Absolute pressure at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PEFIELDU, PEFIELDV, PEFIELDW
                                                    ! Electric field components
LOGICAL, DIMENSION(:,:,:),    INTENT(IN)     :: GATTACH !Recombination and
                                                    !Attachment if true
REAL, DIMENSION(:,:), OPTIONAL, INTENT(IN)   :: PTOWN ! town fraction
REAL, DIMENSION(:,:), OPTIONAL, INTENT(IN)   :: PSEA  ! Land-sea mask

       END SUBROUTINE ION_ATTACH_ELEC
END INTERFACE
END MODULE MODI_ION_ATTACH_ELEC



!       ######################################################################
        SUBROUTINE ION_ATTACH_ELEC(KTCOUNT, KRR, PTSTEP, PRHODREF,           &
                            PRHODJ,PSVS, PRS, PTHT, PCIT, PPABST, PEFIELDU,  &
                            PEFIELDV, PEFIELDW, GATTACH, PTOWN, PSEA      )
!       ######################################################################


!
!!****  * -
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
!-------------------------------------------------------------------------------
!
!*	0.	DECLARATIONS
!		------------
!
use modd_budget,          only : lbudget_sv, NBUDGET_SV1, tbudgets
USE MODD_CONF,            ONLY: CCONF
USE MODD_CST
USE MODD_ELEC_DESCR
USE MODD_ELEC_n
USE MODD_ELEC_PARAM
USE MODD_NSV,             ONLY: NSV_ELECBEG, NSV_ELEC
USE MODD_PARAMETERS,      ONLY: JPHEXT, JPVEXT
USE MODD_RAIN_ICE_DESCR_n
USE MODD_RAIN_ICE_PARAM_n
USE MODD_REF,             ONLY: XTHVREFZ

use mode_budget,          only: Budget_store_init, Budget_store_end
use mode_tools_ll,        only: GET_INDICE_ll

USE MODI_MOMG

IMPLICIT NONE
!
!	0.1	Declaration of arguments
!
INTEGER,                  INTENT(IN)    :: KTCOUNT  ! Temporal loop counter
INTEGER,                  INTENT(IN)    :: KRR      ! Number of moist variables
REAL,                     INTENT(IN)    :: PTSTEP   ! Time step          
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF ! Reference dry air density 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ   ! dry air density* Jacobian
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PSVS   ! Scalar variable vol. source
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRS    ! Moist variable vol. source
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT     ! Theta (K) at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCIT    ! Pristine ice n.c. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABST   ! Absolute pressure at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PEFIELDU, PEFIELDV, PEFIELDW
                                                    ! Electric field components
LOGICAL, DIMENSION(:,:,:),    INTENT(IN)     :: GATTACH !Recombination and
                                                    !Attachment if true
REAL, DIMENSION(:,:), OPTIONAL, INTENT(IN)   :: PTOWN ! town fraction
REAL, DIMENSION(:,:), OPTIONAL, INTENT(IN)   :: PSEA  ! Land-sea mask

!
!
!	0.2	Declaration of local variables
!
REAL,    DIMENSION(:), ALLOCATABLE :: ZT            ! Temperature (K)
REAL,    DIMENSION(:), ALLOCATABLE :: ZCONC, ZVIT, ZRADIUS ! Number concentration
                                                       !fallspeed, radius
REAL                           :: ZCQD, ZCDIF   ! computed coefficients
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
!
!-------------------------------------------------------------------------------
if ( lbudget_sv ) then
  do jrr = 1, nsv_elec
    call Budget_store_init( tbudgets( NBUDGET_SV1 - 1 + nsv_elecbeg - 1 + jrr), 'NEUT', psvs(:, :, :, jrr) )
  end do
end if
!
!*       1.     COMPUTE THE ION RECOMBINATION and TEMPERATURE
!               ---------------------------------------------
!
!
ZCQD = 4 * XPI * XEPSILON * XBOLTZ / XECHARGE
ZCDIF = XBOLTZ /XECHARGE
!
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
IKB = 1 + JPVEXT
IKE = SIZE(PTHT,3) - JPVEXT
!
!*      1.1     Add Ion Recombination source (PSVS in 1/(m3.s))
!               and count and localize valid grid points for ion source terms
!
IVALID = 0
DO IK = IKB, IKE
  DO IJ = IJB, IJE
    DO II = IIB, IIE
      IF (GATTACH(II,IJ,IK)) THEN
! Recombination
        ZCOMB = XIONCOMB * (PSVS(II,IJ,IK,1)*PTSTEP) *    &
                           (PSVS(II,IJ,IK,NSV_ELEC)*PTSTEP) * &
                           PRHODREF(II,IJ,IK) / PRHODJ(II,IJ,IK)
        ZCOMB = MIN(ZCOMB, PSVS(II,IJ,IK,1), PSVS(II,IJ,IK,NSV_ELEC))
        PSVS(II,IJ,IK,1) = PSVS(II,IJ,IK,1) - ZCOMB
        PSVS(II,IJ,IK,NSV_ELEC) = PSVS(II,IJ,IK,NSV_ELEC) - ZCOMB
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
!*      1.2     Compute the temperature
!
IF( IVALID /= 0 ) THEN
  ALLOCATE (ZT(IVALID))
  DO II = 1, IVALID
    ZT(II) = PTHT(IGI(II),IGJ(II),IGK(II)) * &
            (PPABST(IGI(II),IGJ(II),IGK(II)) / XP00) ** (XRD / XCPD)
  ENDDO
END IF
!
!
!*	2.	TRANSFORM VOLUM. SOURCE TERMS INTO MIXING RATIO
!               FOR WATER SPECIES, AND VOLUMIC CONTENT FOR ELECTRIC VARIABLES
!               -------------------------------------------------------------
!
DO JRR = 1, KRR
  PRS(:,:,:,JRR) = PRS(:,:,:,JRR) *PTSTEP / PRHODJ(:,:,:)
ENDDO
!
DO JSV = 1, NSV_ELEC
  PSVS(:,:,:,JSV) = PSVS(:,:,:,JSV) *PTSTEP *PRHODREF(:,:,:) / PRHODJ(:,:,:)
ENDDO
!
!
!*	3.	COMPUTE ATTACHMENT DUE TO ION DIFFUSION AND CONDUCTION
!               ------------------------------------------------------
! 
!  Attachment to cloud droplets, rain, cloud ice, snow, graupel, 
!                and hail (optional)
!
!
IF( IVALID /= 0 ) THEN
!
!*	3.1	Attachment to cloud droplets
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
!*	3.2	Attachment to raindrops, ice crystals, snow, graupel, 
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
!
!*	4.	RETURN TO VOLUMETRIC SOURCE (Prognostic units)
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
!
!*	5.	BUDGET
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
 SUBROUTINE HYDROPARAM (IGRIDX, IGRIDY, IGRIDZ, ZCONC,       &
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
INTEGER, DIMENSION(:), INTENT(IN)     :: IGRIDX, IGRIDY, IGRIDZ ! Index of
                                                   ! valid gridpoints
INTEGER,               INTENT(IN)     :: ITYPE     ! Hydrometeor category
              ! ITYPE= 2: cloud, 3: rain, 4: ice, 5: snow, 6: graupel, 7: hail
REAL, DIMENSION(:), INTENT(INOUT) :: ZCONC, ZVIT, ZRADIUS
!                                      Number concentration, fallspeed, radius
REAL, DIMENSION(:,:), OPTIONAL, INTENT(IN)   :: PTOWN ! town fraction
REAL, DIMENSION(:,:), OPTIONAL, INTENT(IN)   :: PSEA  ! Land-sea mask
!
!*      0.2     declaration of local variables
!
REAL                           :: ZCONC1, ZCONC2 ! for cloud
REAL                           :: ZLBC
REAL                           :: ZFSEDC
REAL                           :: ZRAY 
REAL                           :: ZEXP1, ZEXP2, ZMOM1, ZMOM2
REAL                           :: ZVCOEF, ZRHO00, ZLBI
REAL                           :: ZLAMBDA
INTEGER                        :: JI, JJ, JK, IV
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
    IF (PRESENT(PSEA)) THEN

      ZMOM1 = 0.5*MOMG(XALPHAC,XNUC,1.)
      ZMOM2 = 0.5*MOMG(XALPHAC2,XNUC2,1.)      
      DO IV = 1, IVALID
        JI = IGRIDX(IV)
        JJ = IGRIDY(IV)
        JK = IGRIDZ(IV)
        IF( PRS(JI, JJ, JK, 2)/PRHODREF(JI, JJ, JK) >XRTMIN_ELEC(2) .AND. &
            PSVS(JI, JJ, JK, 2) /=0. ) THEN
          ZCONC1 = PSEA(JI,JJ) * XCONC_SEA + (1. - PSEA(JI,JJ)) * XCONC_LAND
          ZLBC   = PSEA(JI,JJ) * XLBC(2)   + (1. - PSEA(JI,JJ)) * XLBC(1)
          ZFSEDC = PSEA(JI,JJ) * XFSEDC(2) + (1. - PSEA(JI,JJ)) * XFSEDC(1)
          ZFSEDC = MAX(MIN(XFSEDC(1),XFSEDC(2)), ZFSEDC)
          ZCONC2 = (1. - PTOWN(JI,JJ)) * ZCONC1 + PTOWN(JI,JJ) * XCONC_URBAN
          ZRAY   = (1. - PSEA(JI,JJ)) * ZMOM1 + PSEA(JI,JJ) * ZMOM2
          ZCONC (IV) = ZCONC2                 ! Number concentration
          ZLAMBDA    = (ZLBC *ZCONC2 / (PRHODREF(JI,JJ,JK) * PRS(JI,JJ,JK,2)))**XLBEXC
          ZRADIUS (IV) = ZRAY / ZLAMBDA
          ZVIT (IV)    = XCC * ZFSEDC * ZLAMBDA**(-XDC) * &
                                PRHODREF(JI,JJ,JK)**(-XCEXVT)
        END IF
      ENDDO
    ELSE
      ZRAY = 0.5*MOMG(XALPHAC,XNUC,1.)
      ZLBC = XLBC(1) * XCONC_LAND
      DO IV = 1, IVALID
        JI = IGRIDX(IV)
        JJ = IGRIDY(IV)
        JK = IGRIDZ(IV)
        IF( PRS(JI, JJ, JK, 2)/PRHODREF(JI, JJ, JK) >XRTMIN_ELEC(2) .AND. &
            PSVS(JI, JJ, JK, 2) /=0. ) THEN
          ZCONC (IV) = XCONC_LAND             ! Number concentration
          ZLAMBDA    = (ZLBC / (PRHODREF(JI,JJ,JK) * PRS(JI,JJ,JK,2)))**XLBEXC
          ZRADIUS (IV) = ZRAY / ZLAMBDA
          ZVIT (IV)    = XCC * XFSEDC(1) * ZLAMBDA**(-XDC) *    &
                               PRHODREF(JI,JJ,JK)**(-XCEXVT)       
        END IF
      ENDDO
    END IF
!
!
!*	2.	PARAMETERS FOR RAIN
!               -------------------
  CASE (3)
    ZEXP1 = XEXSEDR - 1.
    ZEXP2 = ZEXP1 - XCEXVT
!
    DO IV = 1, IVALID
      JI = IGRIDX(IV)
      JJ = IGRIDY(IV)
      JK = IGRIDZ(IV)
      IF( PRS(JI, JJ, JK, 3)/PRHODREF(JI, JJ, JK) >XRTMIN_ELEC(3) .AND. &
          PSVS(JI, JJ, JK, 3) /=0. ) THEN
         ZLAMBDA = XLBR * (PRHODREF(JI,JJ,JK) * PRS(JI,JJ,JK,3))**XLBEXR
         ZRADIUS (IV) = 0.5 / ZLAMBDA
         ZCONC (IV) = XCCR / ZLAMBDA
         ZVIT (IV) = XFSEDR * PRHODREF(JI,JJ,JK)**ZEXP2 &
                                * PRS(JI,JJ,JK,3)**ZEXP1
      END IF
    ENDDO
!
!
!*	3.	PARAMETERS FOR ICE
!               ------------------
!
  CASE (4)
!
    ZRAY = 0.5*MOMG(XALPHAI,XNUI,1.)
    ZRHO00 = XP00 / (XRD * XTHVREFZ(IKB))
!   ZVCOEF= XC_I * (GAMMA(XNUI+(XBI+XDI)/XALPHAI) / GAMMA(XNUI+XBI/XALPHAI)) &
!                * ZRHO00**XCEXVT
! Computations for Columns (see ini_rain_ice_elec.f90)
    ZVCOEF = 2.1E5 * MOMG(XALPHAI,XNUI, 3.285) / MOMG(XALPHAI,XNUI, 1.7)   &
                   * ZRHO00**XCEXVT
    ZLBI = (2.14E-3 * MOMG(XALPHAI,XNUI,1.7)) **0.588235

    DO IV = 1, IVALID
      JI = IGRIDX(IV)
      JJ = IGRIDY(IV)
      JK = IGRIDZ(IV)
      IF( PRS(JI, JJ, JK, 4)/PRHODREF(JI, JJ, JK) > XRTMIN_ELEC(4) .AND. &
          PSVS(JI, JJ, JK, 4) /=0.) THEN
        ZCONC (IV) = XFCI * PRHODREF(JI,JJ,JK) * PRS(JI,JJ,JK,4) *  &
              MAX(0.05E6, -0.15319E6 - 0.021454E6 *                &
                  ALOG(PRHODREF(JI,JJ,JK) * PRS(JI,JJ,JK,4)))**3
        ZLAMBDA = ZLBI * (ZCONC(IV) / (PRHODREF(JI,JJ,JK) *   &
                                       PRS(JI,JJ,JK,4)))**0.588235
        ZRADIUS (IV) = ZRAY / ZLAMBDA
        ZVIT (IV) = ZVCOEF * ZLAMBDA**(-1.585) *   &      !(-XDI) * &
                             PRHODREF(JI,JJ,JK)**(-XCEXVT) 
      END IF
    ENDDO
!
!
!*	4.	PARAMETERS FOR SNOW 
!               -------------------
!
  CASE (5)
!
    ZEXP1 = XEXSEDS - 1. 
    ZEXP2 = ZEXP1 - XCEXVT
!
    DO IV = 1, IVALID
      JI = IGRIDX(IV)
      JJ = IGRIDY(IV)
      JK = IGRIDZ(IV)
      IF( PRS(JI, JJ, JK, 5)/PRHODREF(JI, JJ, JK) >XRTMIN_ELEC(5) .AND. &
          PSVS(JI, JJ, JK, 5) /=0. ) THEN
        ZLAMBDA = XLBS * (PRHODREF(JI,JJ,JK) * PRS(JI,JJ,JK,5))**XLBEXS
        ZRADIUS (IV) = 0.5 / ZLAMBDA
        ZCONC (IV) = XCCS * ZLAMBDA**XCXS
        ZVIT (IV) = XFSEDS * PRHODREF(JI,JJ,JK)**ZEXP2 &
                               * PRS(JI,JJ,JK,5)**ZEXP1
      END IF
    ENDDO
!
!
!*	5.	PARAMETERS FOR GRAUPEL
!               ----------------------
!
  CASE (6)
!
    ZEXP1 = XEXSEDG - 1.
    ZEXP2 = ZEXP1 - XCEXVT
!
    DO IV = 1, IVALID
      JI = IGRIDX(IV)
      JJ = IGRIDY(IV)
      JK = IGRIDZ(IV)
      IF( PRS(JI, JJ, JK, 6)/PRHODREF(JI, JJ, JK) >XRTMIN_ELEC(6) .AND. &
          PSVS(JI, JJ, JK, 6) /=0. ) THEN
        ZLAMBDA = XLBG * (PRHODREF(JI,JJ,JK) * PRS(JI,JJ,JK,6))**XLBEXG
        ZRADIUS (IV) = 0.5 / ZLAMBDA
        ZCONC (IV) = XCCG * ZLAMBDA**XCXG
        ZVIT (IV) = XFSEDG * PRHODREF(JI,JJ,JK)**ZEXP2 &
                               * PRS(JI,JJ,JK,6)**ZEXP1
      END IF
    ENDDO
!
!
!*	6.	PARAMETERS FOR HAIL
!               -------------------
!
  CASE (7)
!
    ZEXP1 = XEXSEDH - 1.
    ZEXP2 = ZEXP1-XCEXVT
    ZRAY = 0.5*MOMG(XALPHAH, XNUH, 1.)
!
    DO IV = 1, IVALID
      JI = IGRIDX(IV)
      JJ = IGRIDY(IV)
      JK = IGRIDZ(IV)
      IF( PRS(JI, JJ, JK, 7)/PRHODREF(JI, JJ, JK) >XRTMIN_ELEC(7) .AND. &
          PSVS(JI, JJ, JK, 7) /=0. ) THEN
        ZLAMBDA = XLBH * (PRHODREF(JI,JJ,JK) * PRS(JI,JJ,JK,7))**XLBEXH
        ZRADIUS (IV) = ZRAY / ZLAMBDA
        ZCONC (IV) = XCCG * ZLAMBDA**XCXG
        ZVIT (IV) = XFSEDH * PRHODREF(JI,JJ,JK)**ZEXP2 &
                               * PRS(JI,JJ,JK,7)**ZEXP1
      END IF
    ENDDO
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
    IF(ZX /= 0. .AND. ABS(ZX) <= 20.0) THEN
      IF( ABS(ZX) < 1.0E-15) THEN
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
      PQVS(JI,JJ,JK) = PQVS(JI,JJ,JK) + XECHARGE * (ZDELPI - ZDELNI)
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
     PQVS(JI,JJ,JK) = PQVS(JI,JJ,JK) + XECHARGE *(ZDELPI - ZDELNI)
  END IF
ENDDO
!
END SUBROUTINE DIFF_COND
!
!-----------------------------------------------------------------------------
!
END SUBROUTINE ION_ATTACH_ELEC
