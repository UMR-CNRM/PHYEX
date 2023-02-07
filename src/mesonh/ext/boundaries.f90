!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!#####################
MODULE MODI_BOUNDARIES
!#####################
!
INTERFACE
!
      SUBROUTINE BOUNDARIES (                                       &
            PTSTEP,HLBCX,HLBCY,KRR,KSV,KTCOUNT,                     &
            PLBXUM,PLBXVM,PLBXWM,PLBXTHM,PLBXTKEM,PLBXRM,PLBXSVM,   &
            PLBYUM,PLBYVM,PLBYWM,PLBYTHM,PLBYTKEM,PLBYRM,PLBYSVM,   &
            PLBXUS,PLBXVS,PLBXWS,PLBXTHS,PLBXTKES,PLBXRS,PLBXSVS,   &
            PLBYUS,PLBYVS,PLBYWS,PLBYTHS,PLBYTKES,PLBYRS,PLBYSVS,   &
            PRHODJ,PRHODREF,                                        &
            PUT,PVT,PWT,PTHT,PTKET,PRT,PSVT,PSRCT                   )
!
REAL,                  INTENT(IN) :: PTSTEP        ! time step dt
CHARACTER(LEN=4), DIMENSION(2), INTENT(IN) :: HLBCX,HLBCY   ! X and Y-direc. LBC type
!
INTEGER,               INTENT(IN) :: KRR           ! Number of moist  variables
INTEGER,               INTENT(IN) :: KSV           ! Number of Scalar Variables
INTEGER,               INTENT(IN) :: KTCOUNT       ! Temporal loop COUNTer
                                                   ! (=1 at the segment beginning)
!
! Lateral Boundary fields at time t
REAL, DIMENSION(:,:,:),          INTENT(IN) :: PLBXUM,PLBXVM,PLBXWM ! Wind
REAL, DIMENSION(:,:,:),          INTENT(IN) :: PLBXTHM              ! Mass
REAL, DIMENSION(:,:,:),          INTENT(IN) :: PLBYUM,PLBYVM,PLBYWM ! Wind
REAL, DIMENSION(:,:,:),          INTENT(IN) :: PLBYTHM              ! Mass
REAL, DIMENSION(:,:,:),          INTENT(IN) :: PLBXTKEM          ! TKE
REAL, DIMENSION(:,:,:),          INTENT(IN) :: PLBYTKEM
REAL, DIMENSION(:,:,:,:),        INTENT(IN) :: PLBXRM  ,PLBXSVM  ! Moisture and SV
REAL, DIMENSION(:,:,:,:),        INTENT(IN) :: PLBYRM  ,PLBYSVM  ! in x and y-dir.
! temporal derivative of the Lateral Boundary fields
REAL, DIMENSION(:,:,:),          INTENT(IN) :: PLBXUS,PLBXVS,PLBXWS ! Wind
REAL, DIMENSION(:,:,:),          INTENT(IN) :: PLBXTHS              ! Mass
REAL, DIMENSION(:,:,:),          INTENT(IN) :: PLBYUS,PLBYVS,PLBYWS ! Wind
REAL, DIMENSION(:,:,:),          INTENT(IN) :: PLBYTHS              ! Mass
REAL, DIMENSION(:,:,:),          INTENT(IN) :: PLBXTKES          ! TKE
REAL, DIMENSION(:,:,:),          INTENT(IN) :: PLBYTKES
REAL, DIMENSION(:,:,:,:),        INTENT(IN) :: PLBXRS  ,PLBXSVS  ! Moisture and SV
REAL, DIMENSION(:,:,:,:),        INTENT(IN) :: PLBYRS  ,PLBYSVS  ! in x and y-dir.
!
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PRHODJ    ! Jacobian * dry density of
                                                  !  the reference state
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PRHODREF
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PUT,PVT,PWT,PTHT,PTKET,PSRCT
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRT,PSVT
                                                      ! Variables at t
!
END SUBROUTINE BOUNDARIES
!
END INTERFACE
!

END MODULE MODI_BOUNDARIES
!
!
!     ####################################################################
      SUBROUTINE BOUNDARIES (                                       &
            PTSTEP,HLBCX,HLBCY,KRR,KSV,KTCOUNT,                     &
            PLBXUM,PLBXVM,PLBXWM,PLBXTHM,PLBXTKEM,PLBXRM,PLBXSVM,   &
            PLBYUM,PLBYVM,PLBYWM,PLBYTHM,PLBYTKEM,PLBYRM,PLBYSVM,   &
            PLBXUS,PLBXVS,PLBXWS,PLBXTHS,PLBXTKES,PLBXRS,PLBXSVS,   &
            PLBYUS,PLBYVS,PLBYWS,PLBYTHS,PLBYTKES,PLBYRS,PLBYSVS,   &
            PRHODJ,PRHODREF,                                        &
            PUT,PVT,PWT,PTHT,PTKET,PRT,PSVT,PSRCT                   )
!     ####################################################################
!
!!****  *BOUNDARIES* - routine to prepare the Lateral Boundary Conditions for
!!                 all variables at a scalar localization relative to the 
!!                 considered boundary.
!!
!!    PURPOSE
!!    -------
!       Fill up the left and right lateral EXTernal zones, for all prognostic
!       variables, at time t and t-dt, to avoid particular cases close to
!       the Lateral Boundaries in routines computing the evolution terms, in
!       particular in the advection routines.
!
!!**  METHOD
!!    ------
!!      3 different options are proposed: 'WALL' 'CYCL' 'OPEN'
!!                    to define the Boundary Condition type,
!!      though the variables HLBCX and HLBCY (for the X and Y-directions
!!      respectively).
!!       For the 'OPEN' type of LBC, the treatment depends
!!      on the flow configuration: i.e. INFLOW or OUTFLOW conditions.
!!   
!!    EXTERNAL 
!!    --------  
!!    GET_INDICE_ll  : get physical sub-domain bounds
!!    LWEAST_ll,LEAST_ll,LNORTH_ll,LSOUTH_ll : position functions
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------  
!!      Module MODD_PARAMETERS : 
!!        JPHEXT ,JPVEXT 
!!
!!      Module MODD_CONF :
!!        CCONF
!!
!!      Module MODE_UPDATE_NSV :
!!        NSV_CHEM, NSV_CHEMBEG, NSV_CHEMEND
!!
!!      Module MODD_CTURB :
!!        XTKEMIN 
!!
!!    REFERENCE
!!    ---------
!!      Book1 and book2 of documentation (routine BOUNDARIES)
!!
!!    AUTHOR
!!    ------
!!	J.-P. Lafore J. Stein     * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original        17/10/94 
!!      Modification    02/11/94  (J.Stein) copy for t-dt at the external points
!!                                          + change the copy formulation
!!      Modification    18/11/94  (J.Stein) bug correction in the normal velocity
!!                                          prescription in the WALL cases
!!      Modification    13/02/95  (Lafore)  to account for the OPEN case and
!!                                          for the LS fields introduction
!!      Modification    03/03/95  (Mallet)  corrections in variables names in 
!!                                          the Y-OPEN case
!!                      16/03/95  (J.Stein) remove R from the historical variables
!!      Modification    31/05/95  (Lafore)  MASTER_DEV2.1 preparation after the
!!                                          LBC tests performed by I. Mallet
!!      Modification    15/03/96  (Richard) bug correction for OPEN CASE: (TOP Y-LBC)
!!                                          Rv case 
!!      Modification    15/03/96  (Shure)   bug correction for SV variable in
!!                                          open x  right case
!!      Modification    24/10/96  (Masson)  initialization of outer points in
!!                                          wall cases for spawning interpolations
!!      Modification    13/03/97  (Lafore)  "surfacic" LS-fields introduction
!!      Modification    10/04/97  (Lafore)  proper treatment of minima for TKE and EPS
!!      Modification    01/09/97  (Masson)  minimum value for water and passive
!!                                          scalars set to zero at instants M,T
!!      Modification    20/10/97  (Lafore)  introduction of DAVI type of lbc
!!                                           suppression of NEST type
!!      Modification    12/11/97  ( Stein ) use the lB fields
!!      Modification    02/06/98  (Lafore)  declaration of local variables (PLBXUM
!!                                          and PLBXWM do'nt have the same size)
!!      Modification    24/08/98   (Jabouille) parallelize the code 
!!      Modification    20/04/99  ( Stein ) use the same conditions for times t
!!                                          and t-dt
!!      Modification    11/04/00  (Mari)    special conditions for chemical variables
!!      Modification    10/01/01  (Tulet)   update for MOCAGE boundary conditions
!!      Modification    22/01/01  (Gazen)   use NSV_CHEM,NSV_CHEMBEG,NSV_CHEMEND variables
!!      Modification    22/06/01(Jabouille) use XSVMIN
!!      Modification    20/11/01(Gazen & Escobar) rewrite GCHBOUNDARY for portability
!!      Modification    14/03/05 (Tulet)    bug : in case of CYCL do not call ch_boundaries
!!      Modification    14/05/05 (Tulet)    add aerosols / dust
!!      Modification    05/06               Suppression of DAVI type of lbc
!!      Modification    05/06               Remove EPS
!!      Modification    12/2010  (Chong)    Add boundary condition for ions
!!                                          (fair weather profiles)
!!      Modification    07/2013  (Bosseur & Filippi) adds Forefire
!!      Modification    04/2013  (C.Lac)    Remove instant M               
!!      Modification    01/2015  (JL Redelsperger) Introduction of ponderation
!!                                 for non normal velocity and potential temp
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1
!!      Redelsperger & Pianezze : 08/2015 : add XPOND coefficient
!!      Modification    01/2016  (JP Pinty) Add LIMA that is LBC for CCN and IFN
!!      Modification    18/07/17 (Vionnet)  Add blowing snow variables 
!!      Modification    01/2018  (JL Redelsperger) Correction for TKE treatment
!!      Modification    03/02/2020 (B. Vié)  Correction for SV with LIMA
!  P. Wautelet 04/06/2020: correct call to Set_conc_lima
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!         
USE MODD_BLOWSNOW,  ONLY : LBLOWSNOW,NBLOWSNOW_2D
USE MODD_BLOWSNOW_n
USE MODD_CH_AEROSOL , ONLY : LORILAM
USE MODD_CH_MNHC_n,   ONLY : LUSECHEM, LUSECHIC
USE MODD_CONDSAMP,    ONLY : LCONDSAMP
USE MODD_CONF
USE MODD_CTURB
USE MODD_DUST
USE MODD_GRID_n,    ONLY : XZZ
USE MODD_ELEC_DESCR
USE MODD_ELEC_n
#ifdef MNH_FOREFIRE
USE MODD_FOREFIRE,   ONLY : LFOREFIRE
#endif
USE MODD_LBC_n,      ONLY : XPOND
USE MODE_ll
USE MODD_NESTING,      ONLY : NDAD
USE MODD_NSV
USE MODD_PARAMETERS
USE MODD_PARAM_LIMA, ONLY : NMOD_CCN, NMOD_IFN, LBOUND
USE MODD_PARAM_n,    ONLY : CELEC,CCLOUD
USE MODD_PASPOL,      ONLY : LPASPOL
USE MODD_PRECISION,   ONLY: MNHREAL32
USE MODD_REF_n
USE MODD_SALT,        ONLY : LSALT

USE MODE_MODELN_HANDLER
USE MODE_SET_CONC_LIMA

USE MODI_CH_BOUNDARIES
USE MODI_INIT_AEROSOL_CONCENTRATION
USE MODI_ION_BOUNDARIES

IMPLICIT NONE
!
!
!*       0.1   declarations of arguments
!
!
!
!
REAL,                  INTENT(IN) :: PTSTEP        ! time step dt
CHARACTER(LEN=4), DIMENSION(2), INTENT(IN) :: HLBCX,HLBCY   ! X and Y-direc. LBC type
!
INTEGER,               INTENT(IN) :: KRR           ! Number of moist  variables
INTEGER,               INTENT(IN) :: KSV           ! Number of Scalar Variables
INTEGER,               INTENT(IN) :: KTCOUNT       ! Temporal loop COUNTer
                                                   ! (=1 at the segment beginning)
!
! Lateral Boundary fields at time t
REAL, DIMENSION(:,:,:),          INTENT(IN) :: PLBXUM,PLBXVM,PLBXWM ! Wind
REAL, DIMENSION(:,:,:),          INTENT(IN) :: PLBXTHM              ! Mass
REAL, DIMENSION(:,:,:),          INTENT(IN) :: PLBYUM,PLBYVM,PLBYWM ! Wind
REAL, DIMENSION(:,:,:),          INTENT(IN) :: PLBYTHM              ! Mass
REAL, DIMENSION(:,:,:),          INTENT(IN) :: PLBXTKEM          ! TKE
REAL, DIMENSION(:,:,:),          INTENT(IN) :: PLBYTKEM
REAL, DIMENSION(:,:,:,:),        INTENT(IN) :: PLBXRM  ,PLBXSVM  ! Moisture and SV
REAL, DIMENSION(:,:,:,:),        INTENT(IN) :: PLBYRM  ,PLBYSVM  ! in x and y-dir.
! temporal derivative of the Lateral Boundary fields
REAL, DIMENSION(:,:,:),          INTENT(IN) :: PLBXUS,PLBXVS,PLBXWS ! Wind
REAL, DIMENSION(:,:,:),          INTENT(IN) :: PLBXTHS              ! Mass
REAL, DIMENSION(:,:,:),          INTENT(IN) :: PLBYUS,PLBYVS,PLBYWS ! Wind
REAL, DIMENSION(:,:,:),          INTENT(IN) :: PLBYTHS              ! Mass
REAL, DIMENSION(:,:,:),          INTENT(IN) :: PLBXTKES          ! TKE
REAL, DIMENSION(:,:,:),          INTENT(IN) :: PLBYTKES
REAL, DIMENSION(:,:,:,:),        INTENT(IN) :: PLBXRS  ,PLBXSVS  ! Moisture and SV
REAL, DIMENSION(:,:,:,:),        INTENT(IN) :: PLBYRS  ,PLBYSVS  ! in x and y-dir.
!
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PRHODJ    ! Jacobian * dry density of
                                                  !  the reference state
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PRHODREF
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PUT,PVT,PWT,PTHT,PTKET,PSRCT
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRT,PSVT
                                                      ! Variables at t
!
!*       0.2   declarations of local variables
!
INTEGER             :: IIB       ! indice I Beginning in x direction
INTEGER             :: IJB       ! indice J Beginning in y direction
INTEGER             :: IKB       ! indice K Beginning in z direction
INTEGER             :: IIE       ! indice I End       in x direction 
INTEGER             :: IJE       ! indice J End       in y direction 
INTEGER             :: IKE       ! indice K End       in z direction 
INTEGER             :: JEXT      ! Loop index for EXTernal points
INTEGER             :: JRR       ! Loop index for RR variables (water)
INTEGER             :: JSV       ! Loop index for Scalar Variables
INTEGER             :: IMI       ! Model Index
REAL                :: ZTSTEP    ! effective time step
REAL                :: ZPOND      !  Coeff PONDERATION LS 
INTEGER             :: ILBX,ILBY ! size of LB fields' arrays
LOGICAL, SAVE, DIMENSION(:), ALLOCATABLE :: GCHBOUNDARY, GAERBOUNDARY,&
                    GDSTBOUNDARY, GSLTBOUNDARY, GPPBOUNDARY,          &
                    GCSBOUNDARY, GICBOUNDARY, GLIMABOUNDARY,GSNWBOUNDARY 
LOGICAL, SAVE        :: GFIRSTCALL1 = .TRUE.
LOGICAL, SAVE        :: GFIRSTCALL2 = .TRUE.
LOGICAL, SAVE        :: GFIRSTCALL3 = .TRUE.
LOGICAL, SAVE        :: GFIRSTCALL5 = .TRUE. 
LOGICAL, SAVE        :: GFIRSTCALLPP = .TRUE.                         
LOGICAL, SAVE        :: GFIRSTCALLCS = .TRUE.                         
LOGICAL, SAVE        :: GFIRSTCALLIC = .TRUE.                 
LOGICAL, SAVE        :: GFIRSTCALLLIMA = .TRUE.                 
!
REAL, DIMENSION(SIZE(PLBXWM,1),SIZE(PLBXWM,2),SIZE(PLBXWM,3)) ::  &
                       ZLBXVT,ZLBXWT,ZLBXTHT
REAL, DIMENSION(SIZE(PLBYWM,1),SIZE(PLBYWM,2),SIZE(PLBYWM,3)) ::  &
                       ZLBYUT,ZLBYWT,ZLBYTHT
REAL, DIMENSION(SIZE(PLBXTKEM,1),SIZE(PLBXTKEM,2),SIZE(PLBXTKEM,3)) ::  &
                       ZLBXTKET
REAL, DIMENSION(SIZE(PLBYTKEM,1),SIZE(PLBYTKEM,2),SIZE(PLBYTKEM,3)) ::  &
                       ZLBYTKET
REAL, DIMENSION(SIZE(PLBXRM,1),SIZE(PLBXRM,2),SIZE(PLBXRM,3),SIZE(PLBXRM,4)) :: &
                       ZLBXRT
REAL, DIMENSION(SIZE(PLBYRM,1),SIZE(PLBYRM,2),SIZE(PLBYRM,3),SIZE(PLBYRM,4)) :: &
                       ZLBYRT
REAL, DIMENSION(SIZE(PLBXSVM,1),SIZE(PLBXSVM,2),SIZE(PLBXSVM,3),SIZE(PLBXSVM,4)) :: &
                       ZLBXSVT
REAL, DIMENSION(SIZE(PLBYSVM,1),SIZE(PLBYSVM,2),SIZE(PLBYSVM,3),SIZE(PLBYSVM,4)) :: &
                       ZLBYSVT
LOGICAL              :: GCHTMP
LOGICAL              :: GPPTMP
LOGICAL              :: GCSTMP
!
LOGICAL, SAVE        :: GFIRSTCALL4 = .TRUE.
!
#ifdef MNH_FOREFIRE
LOGICAL, SAVE, DIMENSION(:), ALLOCATABLE ::  GFFBOUNDARY
LOGICAL, SAVE        :: GFIRSTCALLFF = .TRUE.                         
LOGICAL              :: GFFTMP
#endif
!
INTEGER              :: JI,JJ
!
REAL, DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3),SIZE(PSVT,4)) :: ZSVT
REAL, DIMENSION(SIZE(PRT,1),SIZE(PRT,2),SIZE(PRT,3),SIZE(PRT,4)) :: ZRT
!
!-------------------------------------------------------------------------------
!
!*       1.    COMPUTE DIMENSIONS OF ARRAYS AND OTHER INDICES:
!              ----------------------------------------------
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
IKB = 1 + JPVEXT
IKE = SIZE(PUT,3) - JPVEXT
IMI = GET_CURRENT_MODEL_INDEX()
!
!-------------------------------------------------------------------------------
!
!*       2.    UPPER AND LOWER BC FILLING:   
!              ---------------------------
!
!*       2.1    COMPUTE THE FIELD EXTRAPOLATIONS AT THE GROUND
!

!
!          at the instant t
!
IF(SIZE(PUT) /= 0) PUT  (:,:,IKB-1)   = PUT  (:,:,IKB) 
IF(SIZE(PVT) /= 0) PVT  (:,:,IKB-1)   = PVT  (:,:,IKB) 
IF(SIZE(PWT) /= 0) PWT  (:,:,IKB-1)   = PWT  (:,:,IKB)  
IF(SIZE(PTHT) /= 0) PTHT (:,:,IKB-1)  = PTHT (:,:,IKB) 
IF(SIZE(PTKET) /= 0) PTKET(:,:,IKB-1) = PTKET(:,:,IKB)
IF(SIZE(PRT) /= 0)  PRT  (:,:,IKB-1,:)= PRT  (:,:,IKB,:)
IF(SIZE(PSVT)/= 0)  PSVT (:,:,IKB-1,:)= PSVT (:,:,IKB,:)
IF(SIZE(PSRCT) /= 0) PSRCT(:,:,IKB-1)   = PSRCT(:,:,IKB)
!
!
!*       2.2    COMPUTE THE FIELD EXTRAPOLATIONS AT THE TOP
!
!          at the instant t
!
IF(SIZE(PWT) /= 0) PWT  (:,:,IKE+1)   = 0.           
IF(SIZE(PUT) /= 0) PUT  (:,:,IKE+1)   = PUT  (:,:,IKE) 
IF(SIZE(PVT) /= 0) PVT  (:,:,IKE+1)   = PVT  (:,:,IKE)
IF(SIZE(PTHT) /= 0) PTHT (:,:,IKE+1)  = PTHT (:,:,IKE)
IF(SIZE(PTKET) /= 0) PTKET(:,:,IKE+1) = PTKET(:,:,IKE)
IF(SIZE(PRT) /= 0) PRT  (:,:,IKE+1,:) = PRT  (:,:,IKE,:)
IF(SIZE(PSVT)/= 0) PSVT (:,:,IKE+1,:) = PSVT (:,:,IKE,:)
IF(SIZE(PSRCT) /= 0) PSRCT(:,:,IKE+1)   = PSRCT(:,:,IKE)

! specific for positive and negative ions mixing ratios (1/kg)

IF (NSV_ELEC .NE. 0) THEN
!
   IF (SIZE(PWT) /= 0) THEN
     WHERE ( PWT(:,:,IKE+1)  .GE. 0.)    ! Outflow
         PSVT (:,:,IKE+1,NSV_ELECBEG) = 2.*PSVT (:,:,IKE,NSV_ELECBEG) -  &
                                           PSVT (:,:,IKE-1,NSV_ELECBEG)
         PSVT (:,:,IKE+1,NSV_ELECEND) = 2.*PSVT (:,:,IKE,NSV_ELECEND) -  &
                                           PSVT (:,:,IKE-1,NSV_ELECEND)
     ELSE WHERE                         ! Inflow from the top
         PSVT (:,:,IKE+1,NSV_ELECBEG) = XCION_POS_FW(:,:,IKE+1)
         PSVT (:,:,IKE+1,NSV_ELECEND) = XCION_NEG_FW(:,:,IKE+1)
     END WHERE
   ENDIF
!
END IF

!
!
!-------------------------------------------------------------------------------
!
!*       3.    COMPUTE LB FIELDS AT TIME T
!              ---------------------------
!
!
IF ( KTCOUNT == 1) THEN
  ZTSTEP = 0.
ELSE
  ZTSTEP = PTSTEP
END IF
!
!
IF ( SIZE(PLBXTHS,1) /= 0 .AND.                      &
   (     HLBCX(1)=='OPEN' .OR. HLBCX(2)=='OPEN')     ) THEN            
  ZLBXVT(:,:,:) = PLBXVM(:,:,:) + ZTSTEP * PLBXVS(:,:,:)
  ZLBXWT(:,:,:) = PLBXWM(:,:,:) + ZTSTEP * PLBXWS(:,:,:)
  ZLBXTHT(:,:,:) = PLBXTHM(:,:,:) + ZTSTEP * PLBXTHS(:,:,:)
  IF ( SIZE(PTKET,1) /= 0 ) THEN
    ZLBXTKET(:,:,:) = PLBXTKEM(:,:,:) + ZTSTEP * PLBXTKES(:,:,:)
  END IF
  IF ( KRR > 0) THEN
    ZLBXRT(:,:,:,:) = PLBXRM(:,:,:,:) + ZTSTEP * PLBXRS(:,:,:,:)
  END IF
  IF ( KSV > 0) THEN
    ZLBXSVT(:,:,:,:) = PLBXSVM(:,:,:,:) + ZTSTEP * PLBXSVS(:,:,:,:)
  END IF
!
ELSE
!
  ZLBXVT(:,:,:) = PLBXVM(:,:,:)
  ZLBXWT(:,:,:) = PLBXWM(:,:,:)
  ZLBXTHT(:,:,:) = PLBXTHM(:,:,:)
  IF ( SIZE(PTKET,1) /= 0 ) THEN
    ZLBXTKET(:,:,:) = PLBXTKEM(:,:,:) 
  END IF
  IF ( KRR > 0) THEN
    ZLBXRT(:,:,:,:) = PLBXRM(:,:,:,:) 
  END IF
  IF ( KSV > 0) THEN
    ZLBXSVT(:,:,:,:) = PLBXSVM(:,:,:,:) 
  END IF
!
END IF
!
!     ============================================================    
!
!  Reproductibility for RSTART -> truncate ZLB to real(knd=4) to have reproductible result
!
ZLBXVT(:,:,:)  = real(ZLBXVT(:,:,:),kind=MNHREAL32)
ZLBXWT(:,:,:)  = real(ZLBXWT(:,:,:),kind=MNHREAL32)
ZLBXTHT(:,:,:) = real(ZLBXTHT(:,:,:),kind=MNHREAL32)
IF ( SIZE(PTKET,1) /= 0 ) THEN
   ZLBXTKET(:,:,:) = real(ZLBXTKET(:,:,:),kind=MNHREAL32)
END IF
IF ( KRR > 0) THEN
   ZLBXRT(:,:,:,:) = real(ZLBXRT(:,:,:,:),kind=MNHREAL32)
END IF
IF ( KSV > 0) THEN
   ZLBXSVT(:,:,:,:) = real(ZLBXSVT(:,:,:,:),kind=MNHREAL32)
END IF
!     ============================================================ 
!
IF ( SIZE(PLBYTHS,1) /= 0 .AND.                       &
     (    HLBCY(1)=='OPEN' .OR. HLBCY(2)=='OPEN'    ))  THEN          
  ZLBYUT(:,:,:) = PLBYUM(:,:,:) + ZTSTEP * PLBYUS(:,:,:)
  ZLBYWT(:,:,:) = PLBYWM(:,:,:) + ZTSTEP * PLBYWS(:,:,:)
  ZLBYTHT(:,:,:) = PLBYTHM(:,:,:) + ZTSTEP * PLBYTHS(:,:,:)
  IF ( SIZE(PTKET,1) /= 0 ) THEN
    ZLBYTKET(:,:,:) = PLBYTKEM(:,:,:) + ZTSTEP * PLBYTKES(:,:,:)
  END IF
  IF ( KRR > 0) THEN
    ZLBYRT(:,:,:,:) = PLBYRM(:,:,:,:) + ZTSTEP * PLBYRS(:,:,:,:)
  END IF
  IF ( KSV > 0) THEN
    ZLBYSVT(:,:,:,:) = PLBYSVM(:,:,:,:) + ZTSTEP * PLBYSVS(:,:,:,:)
  END IF
!
ELSE
!
  ZLBYUT(:,:,:) = PLBYUM(:,:,:)
  ZLBYWT(:,:,:) = PLBYWM(:,:,:)
  ZLBYTHT(:,:,:) = PLBYTHM(:,:,:)
  IF ( SIZE(PTKET,1) /= 0 ) THEN
    ZLBYTKET(:,:,:) = PLBYTKEM(:,:,:) 
  END IF
  IF ( KRR > 0) THEN
    ZLBYRT(:,:,:,:) = PLBYRM(:,:,:,:) 
  END IF
  IF ( KSV > 0) THEN
    ZLBYSVT(:,:,:,:) = PLBYSVM(:,:,:,:) 
  END IF
!
END IF
!
!
!     ============================================================    
!
!  Reproductibility for RSTART -> truncate ZLB to real(knd=4) to have reproductible result
!
ZLBYUT(:,:,:)  = real(ZLBYUT(:,:,:),kind=MNHREAL32)
ZLBYWT(:,:,:)  = real(ZLBYWT(:,:,:),kind=MNHREAL32)
ZLBYTHT(:,:,:) = real(ZLBYTHT(:,:,:),kind=MNHREAL32)
IF ( SIZE(PTKET,1) /= 0 ) THEN
   ZLBYTKET(:,:,:) = real(ZLBYTKET(:,:,:),kind=MNHREAL32)
END IF
IF ( KRR > 0) THEN
   ZLBYRT(:,:,:,:) = real(ZLBYRT(:,:,:,:),kind=MNHREAL32)
END IF
IF ( KSV > 0) THEN
   ZLBYSVT(:,:,:,:) = real(ZLBYSVT(:,:,:,:),kind=MNHREAL32)
END IF
!     ============================================================ 
!
!-------------------------------------------------------------------------------
! PONDERATION COEFF for  Non-Normal velocities and pot temperature
!
ZPOND = XPOND
!
!*       4.    LBC FILLING IN THE X DIRECTION (LEFT WEST SIDE):   
!              ------------------------------------------------
IF (LWEST_ll( )) THEN
!
!
SELECT CASE ( HLBCX(1) )  
!
!*       4.1  WALL CASE:  
!             ========= 
!
  CASE ('WALL')
!
    DO JEXT=1,JPHEXT
       IF(SIZE(PUT) /= 0) PUT  (IIB-JEXT,:,:)   = PUT  (IIB       ,:,:)   ! never used during run
       IF(SIZE(PVT) /= 0) PVT  (IIB-JEXT,:,:)   = PVT  (IIB-1+JEXT,:,:)
       IF(SIZE(PWT) /= 0) PWT  (IIB-JEXT,:,:)   = PWT  (IIB-1+JEXT,:,:)
       IF(SIZE(PTHT) /= 0) PTHT(IIB-JEXT,:,:)   = PTHT (IIB-1+JEXT,:,:)
       IF(SIZE(PTKET)/= 0) PTKET(IIB-JEXT,:,:)  = PTKET(IIB-1+JEXT,:,:)
       IF(SIZE(PRT) /= 0) PRT  (IIB-JEXT,:,:,:) = PRT  (IIB-1+JEXT,:,:,:)
       IF(SIZE(PSVT) /= 0) PSVT(IIB-JEXT,:,:,:) = PSVT (IIB-1+JEXT,:,:,:)
       IF(SIZE(PSRCT) /= 0) PSRCT (IIB-JEXT,:,:)   = PSRCT (IIB-1+JEXT,:,:)
       IF(LBLOWSNOW) XSNWCANO(IIB-JEXT,:,:)     = XSNWCANO(IIB-1+JEXT,:,:)            
!
    END DO
!
    IF(SIZE(PUT) /= 0) PUT(IIB     ,:,:)   = 0.    ! set the normal velocity
!
!
!*       4.2  OPEN CASE:
!             =========
!
  CASE ('OPEN')
!
    IF(SIZE(PUT) /= 0) THEN
       DO JI=JPHEXT,1,-1
          PUT(JI,:,:)=0.
          WHERE ( PUT(IIB,:,:) <= 0. )               !  OUTFLOW condition
             PVT  (JI,:,:) = 2.*PVT  (JI+1,:,:)  -PVT  (JI+2,:,:)
             PWT  (JI,:,:) = 2.*PWT  (JI+1,:,:)  -PWT  (JI+2,:,:)
             PTHT (JI,:,:) = 2.*PTHT (JI+1,:,:)  -PTHT (JI+2,:,:)
             !
          ELSEWHERE                                   !  INFLOW  condition
             PVT  (JI,:,:) = ZPOND*ZLBXVT   (JI,:,:) + (1.-ZPOND)* PVT(JI+1,:,:) ! 1 
             PWT  (JI,:,:) = ZPOND*ZLBXWT   (JI,:,:) + (1.-ZPOND)* PWT(JI+1,:,:) ! 1
             PTHT (JI,:,:) = ZPOND*ZLBXTHT  (JI,:,:) + (1.-ZPOND)* PTHT(JI+1,:,:)! 1
          ENDWHERE
       ENDDO
    ENDIF
!
!
  IF(SIZE(PTKET) /= 0) THEN
     DO JI=JPHEXT,1,-1
        WHERE ( PUT(IIB,:,:) <= 0. )               !  OUTFLOW condition
           PTKET(JI,:,:) = MAX(XTKEMIN, 2.*PTKET(JI+1,:,:)-PTKET(JI+2,:,:))  
        ELSEWHERE                                  !  INFLOW  condition
           PTKET(JI,:,:) = MAX(XTKEMIN, ZPOND*ZLBXTKET(JI,:,:) + (1.-ZPOND)*PTKET(JI+1,:,:))
        ENDWHERE
     ENDDO
  END IF
    !
!                      Case with KRR moist variables 
! 
! 
!
  DO JRR =1 ,KRR
     IF(SIZE(PUT) /= 0) THEN
        DO JI=JPHEXT,1,-1
           WHERE ( PUT(IIB,:,:) <= 0. )         !  OUTFLOW condition
              PRT(JI,:,:,JRR) = MAX(0.,2.*PRT(JI+1,:,:,JRR) -PRT(JI+2,:,:,JRR))
           ELSEWHERE                            !  INFLOW  condition
              PRT(JI,:,:,JRR) = MAX(0.,ZLBXRT(JI,:,:,JRR)) ! 1
           END WHERE
        END DO
     END IF
     !
  END DO
!
  IF(SIZE(PSRCT) /= 0) THEN
     DO JI=JPHEXT,1,-1
        PSRCT (JI,:,:)   = PSRCT (JI+1,:,:)
     END DO
  END IF
!
!                       Case with KSV scalar variables  
  DO JSV=1 ,KSV
    IF(SIZE(PUT) /= 0) THEN
       DO JI=JPHEXT,1,-1
          WHERE ( PUT(IIB,:,:) <= 0. )         !  OUTFLOW condition
             PSVT(JI,:,:,JSV) = MAX(XSVMIN(JSV),2.*PSVT(JI+1,:,:,JSV) - &
                  PSVT(JI+2,:,:,JSV))
          ELSEWHERE                            !  INFLOW  condition
             PSVT(JI,:,:,JSV) = MAX(XSVMIN(JSV),ZLBXSVT(JI,:,:,JSV)) ! 1
          END WHERE
       END DO
    END IF
    !
  END DO
  !
    IF(LBLOWSNOW) THEN
    DO JSV=1 ,NBLOWSNOW_2D
      WHERE ( PUT(IIB,:,IKB) <= 0. )         !  OUTFLOW condition
        XSNWCANO(IIB-1,:,JSV) = MAX(0.,2.*XSNWCANO(IIB,:,JSV) - &
                                             XSNWCANO(IIB+1,:,JSV))
      ELSEWHERE                            !  INFLOW  condition
        XSNWCANO(IIB-1,:,JSV) = 0.         ! Assume no snow enter throug
                                           ! boundaries
      END WHERE
    END DO
   DO JSV=NSV_SNWBEG ,NSV_SNWEND
    IF(SIZE(PUT) /= 0) THEN
      WHERE ( PUT(IIB,:,:) <= 0. )         !  OUTFLOW condition
        PSVT(IIB-1,:,:,JSV) = MAX(0.,2.*PSVT(IIB,:,:,JSV) - &
                                               PSVT(IIB+1,:,:,JSV))
      ELSEWHERE                            !  INFLOW  condition
        PSVT(IIB-1,:,:,JSV) = 0.           ! Assume no snow enter throug
                                           ! boundaries
      END WHERE
    END IF
    !
   END DO
  ENDIF  
!
!
END SELECT
!
END IF
!-------------------------------------------------------------------------------
!
!*       5    LBC FILLING IN THE X DIRECTION (RIGHT EAST SIDE): 
!              ===============--------------------------------
!
IF (LEAST_ll( )) THEN
!
SELECT CASE ( HLBCX(2) ) 
!
!*       5.1  WALL CASE:
!             =========
!
  CASE ('WALL')
!
    DO JEXT=1,JPHEXT
       IF(SIZE(PUT) /= 0) PUT  (IIE+JEXT,:,:)   = PUT  (IIE       ,:,:)   ! never used during run
       IF(SIZE(PVT) /= 0) PVT  (IIE+JEXT,:,:)   = PVT  (IIE+1-JEXT,:,:)
       IF(SIZE(PWT) /= 0) PWT  (IIE+JEXT,:,:)   = PWT  (IIE+1-JEXT,:,:)
       IF(SIZE(PTHT) /= 0) PTHT (IIE+JEXT,:,:)  = PTHT (IIE+1-JEXT,:,:)
       IF(SIZE(PTKET) /= 0) PTKET(IIE+JEXT,:,:) = PTKET(IIE+1-JEXT,:,:)
       IF(SIZE(PRT) /= 0) PRT  (IIE+JEXT,:,:,:) = PRT  (IIE+1-JEXT,:,:,:)
       IF(SIZE(PSVT) /= 0) PSVT(IIE+JEXT,:,:,:) = PSVT (IIE+1-JEXT,:,:,:)
       IF(SIZE(PSRCT) /= 0) PSRCT (IIE+JEXT,:,:)= PSRCT (IIE+1-JEXT,:,:)
       IF(LBLOWSNOW) XSNWCANO(IIE+JEXT,:,:) = XSNWCANO(IIE+1-JEXT,:,:)           
!
    END DO
!
    IF(SIZE(PUT) /= 0) PUT(IIE+1   ,:,:)   = 0.     ! set the normal velocity
!
!*       5.2  OPEN CASE:
!             =========
!
  CASE ('OPEN')
!
     ILBX = SIZE(PLBXVM,1)
     IF(SIZE(PUT) /= 0) THEN
        DO JI=1,JPHEXT
           WHERE ( PUT(IIE+1,:,:) >= 0. )               !  OUTFLOW condition
              PVT  (IIE+JI,:,:) = 2.*PVT  (IIE+JI-1,:,:)   -PVT  (IIE+JI-2,:,:)
              PWT  (IIE+JI,:,:) = 2.*PWT  (IIE+JI-1,:,:)   -PWT  (IIE+JI-2,:,:)
              PTHT (IIE+JI,:,:) = 2.*PTHT (IIE+JI-1,:,:)   -PTHT (IIE+JI-2,:,:)
              !
           ELSEWHERE                                   !  INFLOW  condition
              PVT  (IIE+JI,:,:) = ZPOND*ZLBXVT   (ILBX-JPHEXT+JI,:,:) + (1.-ZPOND)* PVT(IIE+JI-1,:,:)
              PWT  (IIE+JI,:,:) = ZPOND*ZLBXWT   (ILBX-JPHEXT+JI,:,:) + (1.-ZPOND)* PWT(IIE+JI-1,:,:)
              PTHT (IIE+JI,:,:) = ZPOND*ZLBXTHT  (ILBX-JPHEXT+JI,:,:) + (1.-ZPOND)* PTHT(IIE+JI-1,:,:)
           ENDWHERE
        END DO
     ENDIF
     !
     IF(SIZE(PTKET) /= 0) THEN
        ILBX = SIZE(PLBXTKEM,1)
        DO JI=1,JPHEXT
           WHERE ( PUT(IIE+1,:,:) >= 0. )             !  OUTFLOW condition
              PTKET(IIE+JI,:,:) = MAX(XTKEMIN, 2.*PTKET(IIE+JI-1,:,:)-PTKET(IIE+JI-2,:,:))  
           ELSEWHERE                                  !  INFLOW  condition
              PTKET(IIE+JI,:,:) =  MAX(XTKEMIN, ZPOND*ZLBXTKET(ILBX-JPHEXT+JI,:,:) + &
                                   (1.-ZPOND)*PTKET(IIE+JI-1,:,:))
           ENDWHERE
        END DO
     END IF
    !
!
!                      Case with KRR moist variables 
! 
! 
  DO JRR =1 ,KRR
    ILBX=SIZE(PLBXRM,1)
    !
    IF(SIZE(PUT) /= 0) THEN
       DO JI=1,JPHEXT  
          WHERE ( PUT(IIE+1,:,:) >= 0. )       !  OUTFLOW condition
             PRT(IIE+JI,:,:,JRR) = MAX(0.,2.*PRT(IIE+JI-1,:,:,JRR) -PRT(IIE+JI-2,:,:,JRR))
          ELSEWHERE                            !  INFLOW  condition
             PRT(IIE+JI,:,:,JRR) = MAX(0.,ZLBXRT(ILBX-JPHEXT+JI,:,:,JRR))
          END WHERE
       END DO
    END IF
    !
  END DO
!
  IF(SIZE(PSRCT) /= 0) THEN
     DO JI=1,JPHEXT 
        PSRCT (IIE+JI,:,:)   = PSRCT (IIE+JI-1,:,:)
     END DO
  END IF
!                       Case with KSV scalar variables  
  DO JSV=1 ,KSV
    ILBX=SIZE(PLBXSVM,1)
    IF(SIZE(PUT) /= 0) THEN
       DO JI=1,JPHEXT 
          WHERE ( PUT(IIE+1,:,:) >= 0. )       !  OUTFLOW condition
             PSVT(IIE+JI,:,:,JSV) = MAX(XSVMIN(JSV),2.*PSVT(IIE+JI-1,:,:,JSV) - &
                  PSVT(IIE+JI-2,:,:,JSV))
          ELSEWHERE                            !  INFLOW  condition
             PSVT(IIE+JI,:,:,JSV) = MAX(XSVMIN(JSV),ZLBXSVT(ILBX-JPHEXT+JI,:,:,JSV))
          END WHERE
       END DO
    END IF
    !
  END DO
!
  IF(LBLOWSNOW) THEN
    DO JSV=1 ,3
      WHERE ( PUT(IIE+1,:,IKB) >= 0. )       !  OUTFLOW condition
        XSNWCANO(IIE+1,:,JSV) = MAX(0.,2.*XSNWCANO(IIE,:,JSV) - &
                                              XSNWCANO(IIE-1,:,JSV))
      ELSEWHERE                            !  INFLOW  condition
        XSNWCANO(IIE+1,:,JSV) = 0.         ! Assume no snow enter throug
                                           ! boundaries
      END WHERE
    END DO
    DO JSV=NSV_SNWBEG ,NSV_SNWEND
    IF(SIZE(PUT) /= 0) THEN
      WHERE ( PUT(IIE+1,:,:) >= 0. )       !  OUTFLOW condition
        PSVT(IIE+1,:,:,JSV) = MAX(0.,2.*PSVT(IIE,:,:,JSV) - &
                                               PSVT(IIE-1,:,:,JSV))
      ELSEWHERE                            !  INFLOW  condition
        PSVT(IIE+1,:,:,JSV) = 0.           ! Assume no snow enter throug
                                           ! boundaries
      END WHERE
    END IF
    !
  END DO
  END IF
!
END SELECT
!
END IF
!-------------------------------------------------------------------------------
!
!*       6.    LBC FILLING IN THE Y DIRECTION (BOTTOM SOUTH SIDE): 
!              ------------------------------
IF (LSOUTH_ll( )) THEN
!
SELECT CASE ( HLBCY(1) )              
!
!*       6.1  WALL CASE:
!             ========= 
!
  CASE ('WALL')
!
    DO JEXT=1,JPHEXT
       IF(SIZE(PUT) /= 0) PUT  (:,IJB-JEXT,:)   = PUT  (:,IJB-1+JEXT,:)
       IF(SIZE(PVT) /= 0) PVT  (:,IJB-JEXT,:)   = PVT  (:,IJB       ,:)   ! never used during run
       IF(SIZE(PWT) /= 0) PWT  (:,IJB-JEXT,:)   = PWT  (:,IJB-1+JEXT,:)
       IF(SIZE(PTHT) /= 0) PTHT (:,IJB-JEXT,:)  = PTHT (:,IJB-1+JEXT,:)
       IF(SIZE(PTKET) /= 0) PTKET(:,IJB-JEXT,:) = PTKET(:,IJB-1+JEXT,:)
       IF(SIZE(PRT) /= 0) PRT  (:,IJB-JEXT,:,:) = PRT  (:,IJB-1+JEXT,:,:)
       IF(SIZE(PSVT) /= 0) PSVT (:,IJB-JEXT,:,:)= PSVT (:,IJB-1+JEXT,:,:)
       IF(SIZE(PSRCT) /= 0) PSRCT(:,IJB-JEXT,:) = PSRCT(:,IJB-1+JEXT,:)
       IF(LBLOWSNOW) XSNWCANO(:,IJB-JEXT,:)     = XSNWCANO(:,IJB-1+JEXT,:)     
!
    END DO
!
    IF(SIZE(PVT) /= 0) PVT(:,IJB     ,:)   = 0.       ! set the normal velocity
!
!*       6.2  OPEN CASE:
!             =========
!
  CASE ('OPEN')
!
    IF(SIZE(PVT) /= 0) THEN
       DO JJ=JPHEXT,1,-1
          PVT(:,JJ,:)=0.
          WHERE ( PVT(:,IJB,:) <= 0. )               !  OUTFLOW condition
             PUT  (:,JJ,:) = 2.*PUT  (:,JJ+1,:)   -PUT  (:,JJ+2,:)
             PWT  (:,JJ,:) = 2.*PWT  (:,JJ+1,:)   -PWT  (:,JJ+2,:)
             PTHT (:,JJ,:) = 2.*PTHT (:,JJ+1,:)   -PTHT (:,JJ+2,:)
          ELSEWHERE                                   !  INFLOW  condition
             PUT  (:,JJ,:) = ZPOND*ZLBYUT   (:,JJ,:) + (1.-ZPOND)* PUT(:,JJ+1,:)
             PWT  (:,JJ,:) = ZPOND*ZLBYWT   (:,JJ,:) + (1.-ZPOND)* PWT(:,JJ+1,:)
             PTHT (:,JJ,:) = ZPOND*ZLBYTHT  (:,JJ,:) + (1.-ZPOND)* PTHT(:,JJ+1,:)
          ENDWHERE
       END DO
    ENDIF
!
  IF(SIZE(PTKET) /= 0) THEN
     DO JJ=JPHEXT,1,-1
        WHERE ( PVT(:,IJB,:) <= 0. )             !  OUTFLOW condition
           PTKET(:,JJ,:) = MAX(XTKEMIN, 2.*PTKET(:,JJ+1,:)-PTKET(:,JJ+2,:))  
        ELSEWHERE                                !  INFLOW  condition
           PTKET(:,JJ,:) =  MAX(XTKEMIN,ZPOND*ZLBYTKET(:,JJ,:) +  &
                            (1.-ZPOND)*PTKET(:,JJ+1,:))
        ENDWHERE
     END DO
  END IF
    !
!
!                      Case with KRR moist variables 
! 
! 
  DO JRR =1 ,KRR
    IF(SIZE(PVT) /= 0) THEN
       DO JJ=JPHEXT,1,-1
          WHERE ( PVT(:,IJB,:) <= 0. )         !  OUTFLOW condition
             PRT(:,JJ,:,JRR) = MAX(0.,2.*PRT(:,JJ+1,:,JRR) -PRT(:,JJ+2,:,JRR))
          ELSEWHERE                            !  INFLOW  condition
             PRT(:,JJ,:,JRR) = MAX(0.,ZLBYRT(:,JJ,:,JRR))
          END WHERE
       END DO
    END IF
    !
  END DO
!
  IF(SIZE(PSRCT) /= 0) THEN
     DO JJ=JPHEXT,1,-1
        PSRCT(:,JJ,:)   = PSRCT(:,JJ+1,:)
     END DO
  END IF
!
!                       Case with KSV scalar variables  
!
  DO JSV=1 ,KSV
    IF(SIZE(PVT) /= 0) THEN
       DO JJ=JPHEXT,1,-1
          WHERE ( PVT(:,IJB,:) <= 0. )         !  OUTFLOW condition
             PSVT(:,JJ,:,JSV) = MAX(XSVMIN(JSV),2.*PSVT(:,JJ+1,:,JSV) - &
                  PSVT(:,JJ+2,:,JSV))
          ELSEWHERE                            !  INFLOW  condition
             PSVT(:,JJ,:,JSV) = MAX(XSVMIN(JSV),ZLBYSVT(:,JJ,:,JSV))
          END WHERE
       END DO
    END IF
    !
  END DO
!
   IF(LBLOWSNOW) THEN
    DO JSV=1 ,3
      WHERE ( PVT(:,IJB,IKB) <= 0. )         !  OUTFLOW condition
        XSNWCANO(:,IJB-1,JSV) = MAX(0.,2.*XSNWCANO(:,IJB,JSV) - &
                                             XSNWCANO(:,IJB+1,JSV))
      ELSEWHERE                            !  INFLOW  condition
        XSNWCANO(:,IJB-1,JSV) = 0.         ! Assume no snow enter throug
                                           ! boundaries
      END WHERE
    END DO
   DO JSV=NSV_SNWBEG ,NSV_SNWEND
    IF(SIZE(PVT) /= 0) THEN
      WHERE ( PVT(:,IJB,:) <= 0. )         !  OUTFLOW condition
        PSVT(:,IJB-1,:,JSV) = MAX(0.,2.*PSVT(:,IJB,:,JSV) - &
                                               PSVT(:,IJB+1,:,JSV))
      ELSEWHERE                            !  INFLOW  condition
        PSVT(:,IJB-1,:,JSV) = 0.           ! Assume no snow enter throug
                                           ! boundaries
      END WHERE
    END IF
    !
  END DO
  END IF 
!
!
END SELECT
!
END IF
!-------------------------------------------------------------------------------
!
!*       7.    LBC FILLING IN THE Y DIRECTION (TOP NORTH SIDE): 
!              ===============
!
IF (LNORTH_ll( )) THEN
!
SELECT CASE ( HLBCY(2) ) 
!
!*       4.3.1  WALL CASE:
!               ========= 
!
  CASE ('WALL')
!
    DO JEXT=1,JPHEXT
       IF(SIZE(PUT) /= 0) PUT  (:,IJE+JEXT,:)   = PUT  (:,IJE+1-JEXT,:)
       IF(SIZE(PVT) /= 0) PVT  (:,IJE+JEXT,:)   = PVT  (:,IJE       ,:)   ! never used during run
       IF(SIZE(PWT) /= 0) PWT  (:,IJE+JEXT,:)   = PWT  (:,IJE+1-JEXT,:)
       IF(SIZE(PTHT) /= 0) PTHT (:,IJE+JEXT,:)  = PTHT (:,IJE+1-JEXT,:)
       IF(SIZE(PTKET) /= 0) PTKET(:,IJE+JEXT,:) = PTKET(:,IJE+1-JEXT,:)
       IF(SIZE(PRT) /= 0) PRT  (:,IJE+JEXT,:,:) = PRT  (:,IJE+1-JEXT,:,:)
       IF(SIZE(PSVT) /= 0) PSVT (:,IJE+JEXT,:,:)= PSVT (:,IJE+1-JEXT,:,:)
       IF(SIZE(PSRCT) /= 0) PSRCT(:,IJE+JEXT,:) = PSRCT(:,IJE+1-JEXT,:)
       IF(LBLOWSNOW) XSNWCANO(:,IJE+JEXT,:)     = XSNWCANO(:,IJE+1-JEXT,:)   
!
    END DO
!
    IF(SIZE(PVT) /= 0) PVT(:,IJE+1   ,:)   = 0.    ! set the normal velocity
!
!*       4.3.2  OPEN CASE:
!               ========= 
!
  CASE ('OPEN')
!
!
    ILBY=SIZE(PLBYUM,2)
    IF(SIZE(PVT) /= 0) THEN
       DO JJ=1,JPHEXT
          WHERE ( PVT(:,IJE+1,:) >= 0. )               !  OUTFLOW condition
             PUT  (:,IJE+JJ,:) = 2.*PUT  (:,IJE+JJ-1,:)   -PUT  (:,IJE+JJ-2,:)
             PWT  (:,IJE+JJ,:) = 2.*PWT  (:,IJE+JJ-1,:)   -PWT  (:,IJE+JJ-2,:)
             PTHT (:,IJE+JJ,:) = 2.*PTHT (:,IJE+JJ-1,:)   -PTHT (:,IJE+JJ-2,:)
          ELSEWHERE                                   !  INFLOW  condition
             PUT  (:,IJE+JJ,:) = ZPOND*ZLBYUT   (:,ILBY-JPHEXT+JJ,:) + (1.-ZPOND)* PUT(:,IJE+JJ-1,:)
             PWT  (:,IJE+JJ,:) = ZPOND*ZLBYWT   (:,ILBY-JPHEXT+JJ,:) + (1.-ZPOND)* PWT(:,IJE+JJ-1,:)
             PTHT (:,IJE+JJ,:) = ZPOND*ZLBYTHT  (:,ILBY-JPHEXT+JJ,:) + (1.-ZPOND)* PTHT(:,IJE+JJ-1,:)
          ENDWHERE
       END DO
    ENDIF
!
  IF(SIZE(PTKET) /= 0) THEN
    ILBY=SIZE(PLBYTKEM,2)
    DO JJ=1,JPHEXT
       WHERE ( PVT(:,IJE+1,:) >= 0. )             !  OUTFLOW condition
          PTKET(:,IJE+JJ,:) = MAX(XTKEMIN, 2.*PTKET(:,IJE+JJ-1,:)-PTKET(:,IJE+JJ-2,:))  
       ELSEWHERE                                  !  INFLOW  condition
          PTKET(:,IJE+JJ,:) = MAX(XTKEMIN,ZPOND*ZLBYTKET(:,ILBY-JPHEXT+JJ,:) + &
                              (1.-ZPOND)*PTKET(:,IJE+JJ-1,:))
       ENDWHERE
    END DO
  ENDIF
    !
!                      Case with KRR moist variables 
! 
! 
  DO JRR =1 ,KRR
    ILBY=SIZE(PLBYRM,2)
    !
    IF(SIZE(PVT) /= 0) THEN
       DO JJ=1,JPHEXT
          WHERE ( PVT(:,IJE+1,:) >= 0. )         !  OUTFLOW condition
             PRT(:,IJE+JJ,:,JRR) = MAX(0.,2.*PRT(:,IJE+JJ-1,:,JRR) -PRT(:,IJE+JJ-2,:,JRR))
          ELSEWHERE                            !  INFLOW  condition
             PRT(:,IJE+JJ,:,JRR) = MAX(0.,ZLBYRT(:,ILBY-JPHEXT+JJ,:,JRR))
          END WHERE
       END DO
    END IF
    !
  END DO
!
  IF(SIZE(PSRCT) /= 0) THEN
      DO JJ=1,JPHEXT
         PSRCT(:,IJE+JJ,:)   = PSRCT(:,IJE+JJ-1,:)
      END DO
  END IF
!
!                       Case with KSV scalar variables  
  DO JSV=1 ,KSV
    ILBY=SIZE(PLBYSVM,2)
    !
    IF(SIZE(PVT) /= 0) THEN
       DO JJ=1,JPHEXT
          WHERE ( PVT(:,IJE+1,:) >= 0. )         !  OUTFLOW condition
             PSVT(:,IJE+JJ,:,JSV) = MAX(XSVMIN(JSV),2.*PSVT(:,IJE+JJ-1,:,JSV) - &
                  PSVT(:,IJE+JJ-2,:,JSV))
          ELSEWHERE                            !  INFLOW  condition
             PSVT(:,IJE+JJ,:,JSV) = MAX(XSVMIN(JSV),ZLBYSVT(:,ILBY-JPHEXT+JJ,:,JSV))
          END WHERE
       END DO
    END IF
    !
  END DO
!
    IF(LBLOWSNOW) THEN
  DO JSV=1 ,3
    WHERE ( PVT(:,IJE+1,IKB) >= 0. )         !  OUTFLOW condition
      XSNWCANO(:,IJE+1,JSV) = MAX(0.,2.*XSNWCANO(:,IJE,JSV) - &
                                            XSNWCANO(:,IJE-1,JSV))
    ELSEWHERE                            !  INFLOW  condition
      XSNWCANO(:,IJE+1,JSV) = 0.         ! Assume no snow enter throug
                                         ! boundaries
    END WHERE
  END DO
  DO JSV=NSV_SNWBEG ,NSV_SNWEND
    !
    IF(SIZE(PVT) /= 0) THEN
      WHERE ( PVT(:,IJE+1,:) >= 0. )         !  OUTFLOW condition
        PSVT(:,IJE+1,:,JSV) = MAX(0.,2.*PSVT(:,IJE,:,JSV) - &
                                               PSVT(:,IJE-1,:,JSV))
      ELSEWHERE                            !  INFLOW  condition
        PSVT(:,IJE+1,:,JSV) = 0.           ! Assume no snow enter throug
                                           ! boundaries
      END WHERE
    END IF
    !
  END DO
  ENDIF
!
END SELECT
END IF
!
!
IF (CCLOUD == 'LIMA' .AND. IMI == 1 .AND. CPROGRAM=='MESONH') THEN

   ZSVT=PSVT
   ZRT=PRT

  IF (GFIRSTCALLLIMA) THEN
    ALLOCATE(GLIMABOUNDARY(NSV_LIMA))
    GFIRSTCALLLIMA = .FALSE.
    DO JSV=NSV_LIMA_BEG,NSV_LIMA_END
       GCHTMP = .FALSE.
       IF (LWEST_ll().AND.HLBCX(1)=='OPEN')  GCHTMP = GCHTMP .OR. ALL(PLBXSVM(JPHEXT,:,:,JSV)==0)
       IF (LEAST_ll().AND.HLBCX(2)=='OPEN')  GCHTMP = GCHTMP .OR. ALL(PLBXSVM(ILBX-JPHEXT+1,:,:,JSV)==0)
       IF (LSOUTH_ll().AND.HLBCY(1)=='OPEN') GCHTMP = GCHTMP .OR. ALL(PLBYSVM(:,JPHEXT,:,JSV)==0)
       IF (LNORTH_ll().AND.HLBCY(2)=='OPEN') GCHTMP = GCHTMP .OR. ALL(PLBYSVM(:,ILBY-JPHEXT+1,:,JSV)==0)
       GLIMABOUNDARY(JSV-NSV_LIMA_BEG+1) = GCHTMP
    ENDDO
  ENDIF
  CALL INIT_AEROSOL_CONCENTRATION(PRHODREF,ZSVT,XZZ)
  DO JSV=NSV_LIMA_CCN_FREE,NSV_LIMA_CCN_FREE+NMOD_CCN-1 ! LBC for CCN
     IF (GLIMABOUNDARY(JSV-NSV_LIMA_BEG+1)) THEN
        PSVT(IIB-1,:,:,JSV)=ZSVT(IIB-1,:,:,JSV)
        PSVT(IIE+1,:,:,JSV)=ZSVT(IIE+1,:,:,JSV)
        PSVT(:,IJB-1,:,JSV)=ZSVT(:,IJB-1,:,JSV)
        PSVT(:,IJE+1,:,JSV)=ZSVT(:,IJE+1,:,JSV)
     ENDIF
  END DO
  DO JSV=NSV_LIMA_IFN_FREE,NSV_LIMA_IFN_FREE+NMOD_IFN-1 ! LBC for IFN
     IF (GLIMABOUNDARY(JSV-NSV_LIMA_BEG+1)) THEN
        PSVT(IIB-1,:,:,JSV)=ZSVT(IIB-1,:,:,JSV)
        PSVT(IIE+1,:,:,JSV)=ZSVT(IIE+1,:,:,JSV)
        PSVT(:,IJB-1,:,JSV)=ZSVT(:,IJB-1,:,JSV)
        PSVT(:,IJE+1,:,JSV)=ZSVT(:,IJE+1,:,JSV)
     ENDIF
  END DO

  CALL SET_CONC_LIMA( IMI, 'NONE', PRHODREF, ZRT(:, :, :, :), ZSVT(:, :, :, NSV_LIMA_BEG:NSV_LIMA_END) )
  IF (NSV_LIMA_NC.GE.1) THEN
     IF (GLIMABOUNDARY(NSV_LIMA_NC-NSV_LIMA_BEG+1)) THEN
        PSVT(IIB-1,:,:,NSV_LIMA_NC)=ZSVT(IIB-1,:,:,NSV_LIMA_NC) ! cloud
        PSVT(IIE+1,:,:,NSV_LIMA_NC)=ZSVT(IIE+1,:,:,NSV_LIMA_NC)
        PSVT(:,IJB-1,:,NSV_LIMA_NC)=ZSVT(:,IJB-1,:,NSV_LIMA_NC)
        PSVT(:,IJE+1,:,NSV_LIMA_NC)=ZSVT(:,IJE+1,:,NSV_LIMA_NC)
     ENDIF
  ENDIF
  IF (NSV_LIMA_NR.GE.1) THEN
     IF (GLIMABOUNDARY(NSV_LIMA_NR-NSV_LIMA_BEG+1)) THEN
        PSVT(IIB-1,:,:,NSV_LIMA_NR)=ZSVT(IIB-1,:,:,NSV_LIMA_NR) ! rain
        PSVT(IIE+1,:,:,NSV_LIMA_NR)=ZSVT(IIE+1,:,:,NSV_LIMA_NR)
        PSVT(:,IJB-1,:,NSV_LIMA_NR)=ZSVT(:,IJB-1,:,NSV_LIMA_NR)
        PSVT(:,IJE+1,:,NSV_LIMA_NR)=ZSVT(:,IJE+1,:,NSV_LIMA_NR)
     ENDIF
  ENDIF
  IF (NSV_LIMA_NI.GE.1) THEN
     IF (GLIMABOUNDARY(NSV_LIMA_NI-NSV_LIMA_BEG+1)) THEN
        PSVT(IIB-1,:,:,NSV_LIMA_NI)=ZSVT(IIB-1,:,:,NSV_LIMA_NI) ! ice
        PSVT(IIE+1,:,:,NSV_LIMA_NI)=ZSVT(IIE+1,:,:,NSV_LIMA_NI)
        PSVT(:,IJB-1,:,NSV_LIMA_NI)=ZSVT(:,IJB-1,:,NSV_LIMA_NI)
        PSVT(:,IJE+1,:,NSV_LIMA_NI)=ZSVT(:,IJE+1,:,NSV_LIMA_NI)
     ENDIF
  END IF
END IF
!
!
IF (LUSECHEM .AND. IMI == 1) THEN
  IF (GFIRSTCALL1) THEN
    ALLOCATE(GCHBOUNDARY(NSV_CHEM))
    GFIRSTCALL1 = .FALSE.
    DO JSV=NSV_CHEMBEG,NSV_CHEMEND
       GCHTMP = .FALSE.
       IF (LWEST_ll().AND.HLBCX(1)=='OPEN')  GCHTMP = GCHTMP .OR. ALL(PLBXSVM(JPHEXT,:,:,JSV)==0)
       IF (LEAST_ll().AND.HLBCX(2)=='OPEN')  GCHTMP = GCHTMP .OR. ALL(PLBXSVM(ILBX-JPHEXT+1,:,:,JSV)==0)
       IF (LSOUTH_ll().AND.HLBCY(1)=='OPEN') GCHTMP = GCHTMP .OR. ALL(PLBYSVM(:,JPHEXT,:,JSV)==0)
       IF (LNORTH_ll().AND.HLBCY(2)=='OPEN') GCHTMP = GCHTMP .OR. ALL(PLBYSVM(:,ILBY-JPHEXT+1,:,JSV)==0)
       GCHBOUNDARY(JSV-NSV_CHEMBEG+1) = GCHTMP
    ENDDO
  ENDIF

  DO JSV=NSV_CHEMBEG,NSV_CHEMEND
    IF (GCHBOUNDARY(JSV-NSV_CHEMBEG+1))  THEN
      IF (SIZE(PSVT)>0) THEN
        CALL CH_BOUNDARIES (HLBCX,HLBCY,PUT,PVT,PSVT(:,:,:,JSV),XSVMIN(JSV))
      ENDIF
    ENDIF
  ENDDO
ENDIF
!
IF (LUSECHIC .AND. IMI == 1) THEN
  IF (GFIRSTCALLIC) THEN
    ALLOCATE(GICBOUNDARY(NSV_CHIC))
    GFIRSTCALLIC = .FALSE.
    DO JSV=NSV_CHICBEG,NSV_CHICEND
       GCHTMP = .FALSE.
       IF (LWEST_ll().AND.HLBCX(1)=='OPEN')  GCHTMP = GCHTMP .OR. ALL(PLBXSVM(JPHEXT,:,:,JSV)==0)
       IF (LEAST_ll().AND.HLBCX(2)=='OPEN')  GCHTMP = GCHTMP .OR. ALL(PLBXSVM(ILBX-JPHEXT+1,:,:,JSV)==0)
       IF (LSOUTH_ll().AND.HLBCY(1)=='OPEN') GCHTMP = GCHTMP .OR. ALL(PLBYSVM(:,JPHEXT,:,JSV)==0)
       IF (LNORTH_ll().AND.HLBCY(2)=='OPEN') GCHTMP = GCHTMP .OR. ALL(PLBYSVM(:,ILBY-JPHEXT+1,:,JSV)==0)
       GICBOUNDARY(JSV-NSV_CHICBEG+1) = GCHTMP
    ENDDO
  ENDIF

  DO JSV=NSV_CHICBEG,NSV_CHICEND
    IF (GICBOUNDARY(JSV-NSV_CHICBEG+1))  THEN
      IF (SIZE(PSVT)>0) THEN
        CALL CH_BOUNDARIES (HLBCX,HLBCY,PUT,PVT,PSVT(:,:,:,JSV),XSVMIN(JSV))
      ENDIF
    ENDIF
  ENDDO
ENDIF
IF (LORILAM .AND. IMI == 1) THEN
  IF (GFIRSTCALL2) THEN
    ALLOCATE(GAERBOUNDARY(NSV_AER))
    GFIRSTCALL2 = .FALSE.
    DO JSV=NSV_AERBEG,NSV_AEREND
       GCHTMP = .FALSE.
       IF (LWEST_ll().AND.HLBCX(1)=='OPEN')  GCHTMP = GCHTMP .OR. ALL(PLBXSVM(JPHEXT,:,:,JSV)==0)
       IF (LEAST_ll().AND.HLBCX(2)=='OPEN')  GCHTMP = GCHTMP .OR. ALL(PLBXSVM(ILBX-JPHEXT+1,:,:,JSV)==0)
       IF (LSOUTH_ll().AND.HLBCY(1)=='OPEN') GCHTMP = GCHTMP .OR. ALL(PLBYSVM(:,JPHEXT,:,JSV)==0)
       IF (LNORTH_ll().AND.HLBCY(2)=='OPEN') GCHTMP = GCHTMP .OR. ALL(PLBYSVM(:,ILBY-JPHEXT+1,:,JSV)==0)
       GAERBOUNDARY(JSV-NSV_AERBEG+1) = GCHTMP
    ENDDO
  ENDIF

  DO JSV=NSV_AERBEG,NSV_AEREND
    IF (GAERBOUNDARY(JSV-NSV_AERBEG+1)) THEN
      IF (SIZE(PSVT)>0) THEN
        CALL CH_BOUNDARIES (HLBCX,HLBCY,PUT,PVT,PSVT(:,:,:,JSV),XSVMIN(JSV))
      ENDIF
    ENDIF
  ENDDO
ENDIF
!
IF (LDUST .AND. IMI == 1) THEN
  IF (GFIRSTCALL3) THEN
    ALLOCATE(GDSTBOUNDARY(NSV_DST))
    GFIRSTCALL3 = .FALSE.
    DO JSV=NSV_DSTBEG,NSV_DSTEND
       GCHTMP = .FALSE.
       IF (LWEST_ll().AND.HLBCX(1)=='OPEN')  GCHTMP = GCHTMP .OR. ALL(PLBXSVM(JPHEXT,:,:,JSV)==0)
       IF (LEAST_ll().AND.HLBCX(2)=='OPEN')  GCHTMP = GCHTMP .OR. ALL(PLBXSVM(ILBX-JPHEXT+1,:,:,JSV)==0)
       IF (LSOUTH_ll().AND.HLBCY(1)=='OPEN') GCHTMP = GCHTMP .OR. ALL(PLBYSVM(:,JPHEXT,:,JSV)==0)
       IF (LNORTH_ll().AND.HLBCY(2)=='OPEN') GCHTMP = GCHTMP .OR. ALL(PLBYSVM(:,ILBY-JPHEXT+1,:,JSV)==0)
       GDSTBOUNDARY(JSV-NSV_DSTBEG+1) = GCHTMP
    ENDDO
    ENDIF

  DO JSV=NSV_DSTBEG,NSV_DSTEND
    IF (GDSTBOUNDARY(JSV-NSV_DSTBEG+1)) THEN
      IF (SIZE(PSVT)>0) THEN
        CALL CH_BOUNDARIES (HLBCX,HLBCY,PUT,PVT,PSVT(:,:,:,JSV),XSVMIN(JSV))
      ENDIF
    ENDIF
  ENDDO 
ENDIF
!
IF (LSALT .AND. IMI == 1) THEN
  IF (GFIRSTCALL5) THEN
    ALLOCATE(GSLTBOUNDARY(NSV_SLT))
    GFIRSTCALL5 = .FALSE.
    DO JSV=NSV_SLTBEG,NSV_SLTEND
       GCHTMP = .FALSE.
       IF (LWEST_ll().AND.HLBCX(1)=='OPEN')  GCHTMP = GCHTMP .OR. ALL(PLBXSVM(JPHEXT,:,:,JSV)==0)
       IF (LEAST_ll().AND.HLBCX(2)=='OPEN')  GCHTMP = GCHTMP .OR. ALL(PLBXSVM(ILBX-JPHEXT+1,:,:,JSV)==0)
       IF (LSOUTH_ll().AND.HLBCY(1)=='OPEN') GCHTMP = GCHTMP .OR. ALL(PLBYSVM(:,JPHEXT,:,JSV)==0)
       IF (LNORTH_ll().AND.HLBCY(2)=='OPEN') GCHTMP = GCHTMP .OR. ALL(PLBYSVM(:,ILBY-JPHEXT+1,:,JSV)==0)
       GSLTBOUNDARY(JSV-NSV_SLTBEG+1) = GCHTMP
    ENDDO
  ENDIF

  DO JSV=NSV_SLTBEG,NSV_SLTEND
    IF (GSLTBOUNDARY(JSV-NSV_SLTBEG+1)) THEN
      IF (SIZE(PSVT)>0) THEN
        CALL CH_BOUNDARIES (HLBCX,HLBCY,PUT,PVT,PSVT(:,:,:,JSV),XSVMIN(JSV))
      ENDIF
    ENDIF
  ENDDO
ENDIF
!
IF ( LPASPOL .AND. IMI == 1) THEN
  IF (GFIRSTCALLPP) THEN
    ALLOCATE(GPPBOUNDARY(NSV_PP))
    GFIRSTCALLPP = .FALSE.
    DO JSV=NSV_PPBEG,NSV_PPEND
      GPPTMP = .FALSE.
      IF (LWEST_ll().AND.HLBCX(1)=='OPEN') GPPTMP = GPPTMP .OR. ALL(PLBXSVM(JPHEXT,:,:,JSV)==0)
      IF (LEAST_ll().AND.HLBCX(2)=='OPEN') GPPTMP = GPPTMP .OR. ALL(PLBXSVM(ILBX-JPHEXT+1,:,:,JSV)==0)
      IF (LSOUTH_ll().AND.HLBCY(1)=='OPEN') GPPTMP = GPPTMP .OR. ALL(PLBYSVM(:,JPHEXT,:,JSV)==0)
      IF (LNORTH_ll().AND.HLBCY(2)=='OPEN') GPPTMP = GPPTMP .OR. ALL(PLBYSVM(:,ILBY-JPHEXT+1,:,JSV)==0)
      GPPBOUNDARY(JSV-NSV_PPBEG+1) = GPPTMP
    ENDDO
  ENDIF

  DO JSV=NSV_PPBEG,NSV_PPEND
    IF (GPPBOUNDARY(JSV-NSV_PPBEG+1)) THEN
      IF (SIZE(PSVT)>0) THEN
        CALL CH_BOUNDARIES (HLBCX,HLBCY,PUT,PVT,PSVT(:,:,:,JSV),XSVMIN(JSV))
      ENDIF
    ENDIF
  ENDDO
ENDIF
!
IF ( LCONDSAMP .AND. IMI == 1) THEN
  IF (GFIRSTCALLCS) THEN
    ALLOCATE(GCSBOUNDARY(NSV_CS))
    GFIRSTCALLCS = .FALSE.
    DO JSV=NSV_CSBEG,NSV_CSEND
      GCSTMP = .FALSE.
      IF (LWEST_ll().AND.HLBCX(1)=='OPEN') GCSTMP = GCSTMP .OR. ALL(PLBXSVM(JPHEXT,:,:,JSV)==0)
      IF (LEAST_ll().AND.HLBCX(2)=='OPEN') GCSTMP = GCSTMP .OR. ALL(PLBXSVM(ILBX-JPHEXT+1,:,:,JSV)==0)
      IF (LSOUTH_ll().AND.HLBCY(1)=='OPEN') GCSTMP = GCSTMP .OR. ALL(PLBYSVM(:,JPHEXT,:,JSV)==0)
      IF (LNORTH_ll().AND.HLBCY(2)=='OPEN') GCSTMP = GCSTMP .OR. ALL(PLBYSVM(:,ILBY-JPHEXT+1,:,JSV)==0)
      GCSBOUNDARY(JSV-NSV_CSBEG+1) = GCSTMP
    ENDDO
  ENDIF

  DO JSV=NSV_CSBEG,NSV_CSEND
    IF (GCSBOUNDARY(JSV-NSV_CSBEG+1)) THEN
      IF (SIZE(PSVT)>0) THEN
        CALL CH_BOUNDARIES (HLBCX,HLBCY,PUT,PVT,PSVT(:,:,:,JSV),XSVMIN(JSV))
      ENDIF
    ENDIF
  ENDDO
ENDIF

IF (LBLOWSNOW .AND. IMI == 1) THEN
  IF (GFIRSTCALL3) THEN
    ALLOCATE(GSNWBOUNDARY(NSV_SNW))
    GFIRSTCALL3 = .FALSE.
    DO JSV=NSV_SNWBEG,NSV_SNWEND
       GCHTMP = .FALSE.
       IF (LWEST_ll().AND.HLBCX(1)=='OPEN')  GCHTMP = GCHTMP .OR. ALL(PLBXSVM(1,:,:,JSV)==0)
       IF (LEAST_ll().AND.HLBCX(2)=='OPEN')  GCHTMP = GCHTMP .OR. ALL(PLBXSVM(ILBX,:,:,JSV)==0)
       IF (LSOUTH_ll().AND.HLBCY(1)=='OPEN') GCHTMP = GCHTMP .OR. ALL(PLBYSVM(:,1,:,JSV)==0)
       IF (LNORTH_ll().AND.HLBCY(2)=='OPEN') GCHTMP = GCHTMP .OR. ALL(PLBYSVM(:,ILBY,:,JSV)==0)
       GSNWBOUNDARY(JSV-NSV_SNWBEG+1) = GCHTMP
    ENDDO
  ENDIF
ENDIF

#ifdef MNH_FOREFIRE
!ForeFire 
IF ( LFOREFIRE .AND. IMI == 1) THEN
  IF (GFIRSTCALLFF) THEN
    ALLOCATE(GFFBOUNDARY(NSV_FF))
    GFIRSTCALLFF = .FALSE.
    DO JSV=NSV_FFBEG,NSV_FFEND
      GFFTMP = .FALSE.
      IF (LWEST_ll().AND.HLBCX(1)=='OPEN') GFFTMP = GFFTMP .OR. ALL(PLBXSVM(JPHEXT,:,:,JSV)==0)
      IF (LEAST_ll().AND.HLBCX(2)=='OPEN') GFFTMP = GFFTMP .OR. ALL(PLBXSVM(ILBX-JPHEXT+1,:,:,JSV)==0)
      IF (LSOUTH_ll().AND.HLBCY(1)=='OPEN') GFFTMP = GFFTMP .OR. ALL(PLBYSVM(:,JPHEXT,:,JSV)==0)
      IF (LNORTH_ll().AND.HLBCY(2)=='OPEN') GFFTMP = GFFTMP .OR. ALL(PLBYSVM(:,ILBY-JPHEXT+1,:,JSV)==0)
      GFFBOUNDARY(JSV-NSV_FFBEG+1) = GFFTMP
    ENDDO
  ENDIF

  DO JSV=NSV_FFBEG,NSV_FFEND
    IF (GFFBOUNDARY(JSV-NSV_FFBEG+1)) THEN
      IF (SIZE(PSVT)>0) THEN
        CALL CH_BOUNDARIES (HLBCX,HLBCY,PUT,PVT,PSVT(:,:,:,JSV),XSVMIN(JSV))
      ENDIF
    ENDIF
  ENDDO
ENDIF
#endif
!
IF ( CELEC /= 'NONE' .AND. (NSV_ELEC_A(NDAD(IMI)) == 0 .OR. IMI == 1)) THEN
  CALL ION_BOUNDARIES (HLBCX,HLBCY,PUT,PVT,PSVT)
ENDIF
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE BOUNDARIES
