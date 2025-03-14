!MNH_LIC Copyright 1994-2024 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!    #################### 
     MODULE MODE_PRANDTL
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
!    #################### 
!
!* modification 08/2010  V. Masson  smoothing of the discontinuity in functions 
!                                   used for implicitation of exchange coefficients
!               05/2020   V. Masson and C. Lac : bug in D_PHI3DTDZ2_O_DDTDZ
!               06/2023   S. Riette add the LSMOOTH_PRANDTL key
!
USE MODD_CTURB,      ONLY : CSTURB_t
USE MODD_TURB_n,     ONLY : TURB_t
USE MODD_DIMPHYEX,   ONLY : DIMPHYEX_t
USE MODD_PARAMETERS, ONLY : JPVEXT_TURB
!
USE MODE_SHUMAN_PHY, ONLY: MZM_PHY,MZF_PHY
USE MODE_GRADIENT_M_PHY
IMPLICIT NONE
!----------------------------------------------------------------------------
CONTAINS
!----------------------------------------------------------------------------
      SUBROUTINE PRANDTL(D,CST,CSTURB,TURBN,KRR,KSV,KRRI,OTURB_DIAG,&
                         HTURBDIM,OOCEAN,OHARAT,O2D,OCOMPUTE_SRC,&
                         TPFILE, OFLAT,                        &
                         PDXX,PDYY,PDZZ,PDZX,PDZY,             &
                         PTHVREF,PLOCPEXNM,PATHETA,PAMOIST,    &
                         PLM,PLEPS,PTKEM,PTHLM,PRM,PSVM,PSRCM, &
                         PREDTH1,PREDR1,                       &
                         PRED2TH3, PRED2R3, PRED2THR3,         &
                         PREDS1,PRED2THS3, PRED2RS3,           &
                         PBLL_O_E,                             &
                         PETHETA, PEMOIST                      )
!     ###########################################################
!
!
!!****  *PRANDTL* - routine to compute the Prandtl turbulent numbers
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to compute the Redelsperger 
!     numbers and then get the turbulent Prandtl and Schmidt numbers:
!       * for the heat fluxes     - PHI3 = 1/ Prandtl
!       * for the moisture fluxes - PSI3 = 1/ Schmidt
!
!!**  METHOD
!!    ------
!!    The following steps are performed:
!!
!!     1 - default values of 1 are taken for phi3 and psi3 and different masks
!!       are defined depending on the presence of turbulence, stratification and
!!       humidity. The 1D Redelsperger numbers are computed  
!!         * ZREDTH1 : (g / THVREF ) (LT**2 / TKE ) ETHETA (D Theta / Dz)   
!!         * ZREDR1  : (g / THVREF ) (LT**2 / TKE ) EMOIST (D TW    / Dz)  
!!     2 - 3D Redelsperger numbers are computed only for turbulent 
!!       grid points where  ZREDTH1 or ZREDR1 are > 0.
!!     3 - PHI3 is  computed only for turbulent grid points where ZREDTH1 > 0
!!      (turbulent thermally stratified points)
!!     4 - PSI3 is computed only for turbulent grid points where ZREDR1 > 0
!!      (turbulent moist points)
!!    
!!
!!    EXTERNAL
!!    --------
!!      FUNCTIONs ETHETA and EMOIST  :  
!!            allows to compute the coefficients 
!!            for the turbulent correlation between any variable
!!            and the virtual potential temperature, of its correlations 
!!            with the conservative potential temperature and the humidity 
!!            conservative variable:
!!            -------              -------              -------
!!            A' Thv'  =  ETHETA   A' Thl'  +  EMOIST   A' Rnp'  
!!
!!      GX_M_M, GY_M_M, GZ_M_M :  Cartesian gradient operators
!!      MZM : Shuman function (mean operator in the z direction)
!!      Module MODI_ETHETA    : interface module for ETHETA
!!      Module MODI_EMOIST    : interface module for EMOIST
!!      Module MODI_SHUMAN    : interface module for Shuman operators
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST : contains physical constants
!!           CST%XG         : gravity constant
!!
!!      Module MODD_CTURB: contains the set of constants for
!!                        the turbulence scheme
!!       TURBN%XCTV,XCPR2    : constants for the turbulent prandtl numbers
!!       XTKEMIN        : minimum value allowed for the TKE
!!
!!      Module MODD_PARAMETERS 
!!       JPVEXT_TURB         : number of vertical marginal points
!!
!!    REFERENCE
!!    ---------
!!      Book 2 of documentation (routine PRANDTL)
!!      Book 1 of documentation (Chapter: Turbulence)
!!
!!    AUTHOR
!!    ------
!!      Joan Cuxart             * INM and Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original         18/10/94
!!      Modifications: Feb 14, 1995 (J.Cuxart and J.Stein) 
!!                                  Doctorization and Optimization
!!      Modifications: March 21, 1995 (J.M. Carriere) 
!!                                  Introduction of cloud water
!!      Modifications: March 21, 1995 (J. Cuxart and J.Stein)
!!                                  Phi3 and Psi3 at w point + cleaning
!!      Modifications: July 2, 1995 (J.Cuxart and Ph.Bougeault)
!!                         change the value of Phi3 and Psi3 if negative
!!      Modifications: Sept 20, 1995 (J. Stein, J. Cuxart, J.L. Redelsperger)
!!                         remove the Where + use REDTH1+REDR1 for the tests
!!      Modifications: October 10, 1995 (J. Cuxart and J.Stein)
!!                         Psi3 for tPREDS1he scalar variables
!!      Modifications: February 27, 1996 (J.Stein) optimization
!!      Modifications: June 15, 1996 (P.Jabouille) return to the previous
!!                          computation of Phi3 and Psi3
!!      Modifications: October 10, 1996 (J. Stein) change the temporal
!!                          discretization
!!      Modifications: May 23, 1997 (J. Stein) bug in 3D Redels number at ground
!!                                             with orography
!!      Modifications: Feb 20, 1998 (J. Stein) bug in all the 3D cases due to
!!                                             the use of ZW1 instead of ZW2
!!                     Feb 20, 2003 (JP Pinty) Add PFRAC_ICE
!!                     July    2005 (Tomas, Masson) implicitation of PHI3 and PSI3
!!                     October 2009 (G. Tanguy) add ILENCH=LEN(YCOMMENT) after
!!                                              change of YCOMMENT
!!                     2012-02 Y. Seity,  add possibility to run with reversed 
!!                                               vertical levels
!!      Modifications: July 2015    (Wim de Rooy) OHARAT (Racmo turbulence) switch
!!                     2017-09 J.Escobar, use epsilon XMNH_TINY_12 for R*4 
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!! JL Redelsperger 03/2021 : adding Ocean case for temperature only 
!! --------------------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
!
!
USE MODD_CST,            ONLY: CST_t
USE MODD_CTURB,          ONLY: CSTURB_t
USE MODD_DIMPHYEX,       ONLY: DIMPHYEX_t
USE MODD_FIELD,          ONLY: TFIELDMETADATA, TYPEREAL
USE MODD_TURB_n,         ONLY: TURB_t
USE MODD_IO,             ONLY: TFILEDATA
!
USE MODE_EMOIST,         ONLY: EMOIST
USE MODE_ETHETA,         ONLY: ETHETA
USE MODE_GRADIENT_M_PHY, ONLY: GX_M_M_PHY, GY_M_M_PHY
USE MODE_IO_FIELD_WRITE_PHY, ONLY: IO_FIELD_WRITE_PHY
!
! These macro are handled by pft_tool.py --craybyPassDOCONCURRENT applied on Cray Rules
#ifdef MNH_COMPILER_CCE
!$mnh_undef(LOOP)
!$mnh_undef(OPENACC)
#endif
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
TYPE(DIMPHYEX_t),       INTENT(IN)   :: D
TYPE(CST_t),            INTENT(IN)   :: CST
TYPE(CSTURB_t),         INTENT(IN)   :: CSTURB
TYPE(TURB_t),           INTENT(IN)   :: TURBN
INTEGER,                INTENT(IN)   :: KSV           ! number of scalar variables

INTEGER,                INTENT(IN)   :: KRR           ! number of moist var.
INTEGER,                INTENT(IN)   :: KRRI          ! number of ice var.
!
LOGICAL,                INTENT(IN)   ::  OTURB_DIAG   ! switch to write some
                                 ! diagnostic fields in the syncronous FM-file
LOGICAL,                INTENT(IN)   ::  OOCEAN       ! switch for Ocean model version
LOGICAL,                INTENT(IN)   ::  OHARAT
LOGICAL,                INTENT(IN)   ::  OCOMPUTE_SRC ! flag to define dimensions of SIGS and
LOGICAL,                INTENT(IN)   :: O2D           ! Logical for 2D model version (modd_conf)
LOGICAL,                INTENT(IN)   ::  OFLAT        ! Logical for zero ororography
CHARACTER(LEN=4),       INTENT(IN)   ::  HTURBDIM     ! Kind of turbulence param.
TYPE(TFILEDATA),        INTENT(INOUT)   ::  TPFILE       ! Output file
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PDXX,PDYY,PDZZ,PDZX,PDZY
                                                  ! metric coefficients
!
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PTHVREF  ! Virtual Potential Temp.
                                                  ! of the reference state
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PLOCPEXNM ! Lv(T)/Cp/Exner at t-1 
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PATHETA      ! coefficients between 
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PAMOIST      ! s and Thetal and Rnp
!
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PLM      ! Turbulent Mixing length
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PLEPS    ! Dissipative length
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PTHLM,PTKEM! Conservative Potential 
                                                  ! Temperature and TKE at t-1
REAL, DIMENSION(D%NIJT,D%NKT,KRR), INTENT(IN) ::  PRM      ! Mixing ratios at  t-1
                                                  ! with PRM(:,:,:,1) = cons.
                                                  ! mixing ratio
REAL, DIMENSION(D%NIJT,D%NKT,KSV), INTENT(IN) ::  PSVM     ! Scalars at t-1      
REAL, DIMENSION(MERGE(D%NIT,0,OCOMPUTE_SRC),&
                MERGE(D%NJT,0,OCOMPUTE_SRC),&
                MERGE(D%NKT,0,OCOMPUTE_SRC)),  INTENT(IN)   ::  PSRCM
                                  ! s'r'c/2Sigma_s2 at t-1 multiplied by Lambda_3
!
!
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)  ::  PREDTH1 ! Redelsperger number R_theta
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)  ::  PREDR1  ! Redelsperger number R_q
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)  ::  PRED2TH3 ! Redelsperger number R*2_theta
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)  ::  PRED2R3  ! Redelsperger number R*2_q
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)  ::  PRED2THR3! Redelsperger number R*2_thq
REAL, DIMENSION(D%NIJT,D%NKT,KSV), INTENT(OUT)::  PREDS1   ! Redelsperger number R_s
REAL, DIMENSION(D%NIJT,D%NKT,KSV), INTENT(OUT)::  PRED2THS3! Redelsperger number R*2_thsv
REAL, DIMENSION(D%NIJT,D%NKT,KSV), INTENT(OUT)::  PRED2RS3 ! Redelsperger number R*2_qsv
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)  ::  PBLL_O_E! beta*Lk*Leps/tke
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)  ::  PETHETA ! coefficient E_theta
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)  ::  PEMOIST ! coefficient E_moist
!
!
!       0.2  declaration of local variables
!
REAL, DIMENSION(D%NIJT,D%NKT) ::  &
                  ZW1, ZW2,  &
!                                                 working variables
                  ZWORK1,ZWORK2,ZWORK3,ZWORK4, &
                  ZGXMM_PTH,ZGYMM_PTH,ZGXMM_PRM,ZGYMM_PRM, ZGXMM_PSV,ZGYMM_PSV
!                                                 working variables for explicit array
!                                                     
INTEGER :: IKB      ! vertical index value for the first inner mass point
INTEGER :: IKE      ! vertical index value for the last inner mass point
INTEGER::  JSV,JIJ,JK ! loop index
INTEGER :: IIJB,IIJE,IKT,IKA,IKL

REAL    :: ZMINVAL
TYPE(TFIELDMETADATA) :: TZFIELD
! ---------------------------------------------------------------------------
!
!*      1.  DEFAULT VALUES,  1D REDELSPERGER NUMBERS 
!           ----------------------------------------
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('PRANDTL',0,ZHOOK_HANDLE)

IF (OHARAT) THEN
PREDTH1(:,:)=0.
PREDR1(:,:)=0.
PRED2TH3(:,:)=0.
PRED2R3(:,:)=0.
PRED2THR3(:,:)=0.
PREDS1(:,:,:)=0.
PRED2THS3(:,:,:)=0.
PRED2RS3(:,:,:)=0.
PBLL_O_E(:,:)=0.
ENDIF
!
IKB=D%NKTB
IKE=D%NKTE 
IIJE=D%NIJE
IIJB=D%NIJB
IKT=D%NKT
IKA=D%NKA
IKL=D%NKL
!
CALL ETHETA(D,CST,KRR,KRRI,PTHLM,PRM,PLOCPEXNM,PATHETA,PSRCM,OOCEAN,OCOMPUTE_SRC,ZWORK1)
CALL EMOIST(D,CST,KRR,KRRI,PTHLM,PRM,PLOCPEXNM,PAMOIST,PSRCM,OOCEAN,ZWORK2)
CALL MZM_PHY(D,ZWORK1,PETHETA)
CALL MZM_PHY(D,ZWORK2,PEMOIST)

DO JIJ=IIJB, IIJE
  PETHETA(JIJ, IKA) = 2.*PETHETA(JIJ, IKB) - PETHETA(JIJ, IKB+IKL)
  PEMOIST(JIJ, IKA) = 2.*PEMOIST(JIJ, IKB) - PEMOIST(JIJ, IKB+IKL)
END DO

!
!---------------------------------------------------------------------------
IF (.NOT. OHARAT) THEN
!
!          1.3 1D Redelsperger numbers
!
IF (OOCEAN) THEN

  DO JK=1, IKT
    DO JIJ=IIJB, IIJE
      ZWORK1(JIJ, JK) = CST%XG * CST%XALPHAOC * PLM(JIJ, JK) & 
                                        * PLEPS(JIJ, JK) / PTKEM(JIJ, JK)
    END DO
  END DO

ELSE

  DO JK=1, IKT
    DO JIJ=IIJB, IIJE
      ZWORK1(JIJ, JK) = CST%XG / PTHVREF(JIJ, JK) * PLM(JIJ, JK) & 
                                        * PLEPS(JIJ, JK) / PTKEM(JIJ, JK)
    END DO
  END DO

END IF
!
CALL MZM_PHY(D,ZWORK1,PBLL_O_E)
CALL GZ_M_W_PHY(D,PTHLM,PDZZ,ZWORK1)
!
IF (OOCEAN) THEN

  DO JK=1, IKT
    DO JIJ=IIJB, IIJE
      PREDTH1(JIJ, JK)= TURBN%XCTV*PBLL_O_E(JIJ, JK)*ZWORK1(JIJ, JK)
    END DO
  END DO
  PREDR1(IIJB:IIJE,1:IKT) = 0.

ELSE
  IF (KRR /= 0) THEN                ! moist case
    CALL GZ_M_W_PHY(D,PRM(:,:,1),PDZZ,ZWORK2)

    DO JK=1, IKT
      DO JIJ=IIJB, IIJE
        PREDTH1(JIJ, JK)= TURBN%XCTV*PBLL_O_E(JIJ, JK) * PETHETA(JIJ, JK) &
                                          * ZWORK1(JIJ, JK)
        PREDR1(JIJ, JK) = TURBN%XCTV*PBLL_O_E(JIJ, JK) * PEMOIST(JIJ, JK) &
                                          * ZWORK2(JIJ, JK)
      END DO
    END DO

  ELSE    
                          ! dry case

    DO JK=1, IKT
      DO JIJ=IIJB, IIJE
        PREDTH1(JIJ, JK)= TURBN%XCTV*PBLL_O_E(JIJ, JK)  * ZWORK1(JIJ, JK)
        PREDR1(JIJ, JK) = 0.
      END DO
    END DO

  END IF
END IF
!
!       3. Limits on 1D Redelperger numbers
!          --------------------------------
!
ZMINVAL = (1.-1./CSTURB%XPHI_LIM)
!

DO JK=1,IKT 
  DO JIJ=IIJB,IIJE 
   ZW1(JIJ,JK) = 1.
   ZW2(JIJ,JK) = 1.
   !
   IF (PREDTH1(JIJ,JK)+PREDR1(JIJ,JK)<-ZMINVAL)THEN
    ZW1(JIJ,JK) = (-ZMINVAL) / (PREDTH1(JIJ,JK)+PREDR1(JIJ,JK))
   ENDIF
   !
   IF (PREDTH1(JIJ,JK)<-ZMINVAL)THEN
    ZW2(JIJ,JK) = (-ZMINVAL) / (PREDTH1(JIJ,JK))
   ENDIF
   ZW2(JIJ,JK) = MIN(ZW1(JIJ,JK),ZW2(JIJ,JK))
   !
   ZW1(JIJ,JK) = 1.
   IF (PREDR1(JIJ,JK)<-ZMINVAL)THEN
    ZW1(JIJ,JK) = (-ZMINVAL) / (PREDR1(JIJ,JK))
   ENDIF
   ZW1(JIJ,JK) = MIN(ZW2(JIJ,JK),ZW1(JIJ,JK))
   !
   !
   !       3. Modification of Mixing length and dissipative length
   !          ----------------------------------------------------
   !
   PBLL_O_E(JIJ,JK) = PBLL_O_E(JIJ,JK) * ZW1(JIJ,JK)
   PREDTH1(JIJ,JK)  = PREDTH1(JIJ,JK)  * ZW1(JIJ,JK)
   PREDR1(JIJ,JK)   = PREDR1(JIJ,JK)   * ZW1(JIJ,JK)
   !
   !       4. Threshold for very small (in absolute value) Redelperger numbers
   !          ----------------------------------------------------------------
   !
   IF(PREDTH1(JIJ,JK) < 0.) THEN
    ZW2(JIJ,JK)=-1.
   ELSE
    ZW2(JIJ,JK)=1.
   END IF
   PREDTH1(JIJ,JK)= ZW2(JIJ,JK) * MAX(CST%XMNH_TINY_12, ZW2(JIJ,JK)*PREDTH1(JIJ,JK))
   !
   IF (KRR /= 0) THEN                ! moist case
    IF(PREDR1(JIJ,JK) < 0.) THEN
     ZW2(JIJ,JK)=-1.
    ELSE
     ZW2(JIJ,JK)=1.
    END IF
    PREDR1(JIJ,JK)= ZW2(JIJ,JK) * MAX(CST%XMNH_TINY_12, ZW2(JIJ,JK)*PREDR1(JIJ,JK))
   END IF
  ENDDO
ENDDO

!
!
!---------------------------------------------------------------------------
!
!          For the scalar variables
DO JSV=1,KSV
  CALL GZ_M_W_PHY(D,PSVM(:,:,JSV),PDZZ,ZWORK1)

  DO JK=1, IKT
    DO JIJ=IIJB, IIJE
      PREDS1(JIJ, JK, JSV)=TURBN%XCTV*PBLL_O_E(JIJ, JK)*ZWORK1(JIJ, JK)
    END DO
  END DO

END DO
!

DO JSV=1,KSV
 DO JK=1,IKT 
   DO JIJ=IIJB,IIJE 
    IF(PREDS1(JIJ,JK,JSV) < 0.) THEN
     ZW2(JIJ,JK)=-1.
    ELSE
     ZW2(JIJ,JK)=1.
    END IF
    PREDS1(JIJ,JK,JSV)= ZW2(JIJ,JK) * MAX(CST%XMNH_TINY_12, ZW2(JIJ,JK)*PREDS1(JIJ,JK,JSV))
   ENDDO
 ENDDO
ENDDO

!
!---------------------------------------------------------------------------
!
!*      2.  3D REDELSPERGER NUMBERS
!           ------------------------
!
IF(HTURBDIM=='1DIM') THEN        ! 1D case
!
!

  DO JK=1, IKT
    DO JIJ=IIJB, IIJE
      PRED2TH3(JIJ, JK)  = PREDTH1(JIJ, JK)**2
    !
      PRED2R3(JIJ, JK)   = PREDR1(JIJ, JK) **2
    !
      PRED2THR3(JIJ, JK) = PREDTH1(JIJ, JK) * PREDR1(JIJ, JK)
    END DO
  END DO

!
ELSE IF (O2D) THEN                      ! 3D case in a 2D model
!
    CALL GX_M_M_PHY(D,OFLAT,PTHLM,PDXX,PDZZ,PDZX,ZGXMM_PTH)

    DO JK=1, IKT
      DO JIJ=IIJB, IIJE
        ZWORK1(JIJ, JK) = ZGXMM_PTH(JIJ, JK)**2
      END DO
    END DO

    CALL MZM_PHY(D,ZWORK1,ZWORK2)
    !
  IF (KRR /= 0) THEN                 ! moist 3D case
    CALL GX_M_M_PHY(D,OFLAT,PRM(:,:,1),PDXX,PDZZ,PDZX,ZGXMM_PRM)

    DO JK=1, IKT
      DO JIJ=IIJB, IIJE
        ZWORK1(JIJ, JK) = ZGXMM_PRM(JIJ, JK)**2
      END DO
    END DO

    CALL MZM_PHY(D,ZWORK1,ZWORK3)
    !

    DO JK=1, IKT
      DO JIJ=IIJB, IIJE
        ZWORK1(JIJ, JK) = ZGXMM_PTH(JIJ, JK) * ZGXMM_PRM(JIJ, JK)
      END DO
    END DO

    CALL MZM_PHY(D,ZWORK1,ZWORK4)
    !

    DO JK=1, IKT
      DO JIJ=IIJB, IIJE
        PRED2TH3(JIJ, JK)= PREDTH1(JIJ, JK)**2+(TURBN%XCTV*PBLL_O_E(JIJ, JK) &
                                           *PETHETA(JIJ, JK) )**2 * ZWORK2(JIJ, JK)
    !
        PRED2R3(JIJ, JK)= PREDR1(JIJ, JK)**2 + (TURBN%XCTV*PBLL_O_E(JIJ, JK) &
                                         * PEMOIST(JIJ, JK))**2 * ZWORK3(JIJ, JK)
    !
        PRED2THR3(JIJ, JK)= PREDR1(JIJ, JK) * PREDTH1(JIJ, JK) +  TURBN%XCTV**2 &
                                            * PBLL_O_E(JIJ, JK)**2   &
                                            * PEMOIST(JIJ, JK) * PETHETA(JIJ, JK) &
                                            * ZWORK4(JIJ, JK)
      END DO
    END DO
!
    PRED2TH3(IIJB:IIJE,IKB)=PRED2TH3(IIJB:IIJE,IKB+IKL)
    PRED2R3(IIJB:IIJE,IKB)=PRED2R3(IIJB:IIJE,IKB+IKL) 
    PRED2THR3(IIJB:IIJE,IKB)=PRED2THR3(IIJB:IIJE,IKB+IKL)

!
  ELSE                 ! dry 3D case in a 2D model

    DO JK=1, IKT
      DO JIJ=IIJB, IIJE
        PRED2TH3(JIJ, JK) = PREDTH1(JIJ, JK)**2 +  TURBN%XCTV**2 & 
                                           * PBLL_O_E(JIJ, JK)**2 * ZWORK2(JIJ, JK)      
      END DO
    END DO
    PRED2TH3(IIJB:IIJE,IKB)=PRED2TH3(IIJB:IIJE,IKB+IKL)
!
    PRED2R3(IIJB:IIJE,1:IKT) = 0.
!
    PRED2THR3(IIJB:IIJE,1:IKT) = 0.

!
  END IF
!
ELSE                                 ! 3D case in a 3D model
!
  CALL GX_M_M_PHY(D,OFLAT,PTHLM,PDXX,PDZZ,PDZX,ZGXMM_PTH)
  CALL GY_M_M_PHY(D,OFLAT,PTHLM,PDYY,PDZZ,PDZY,ZGYMM_PTH)

  DO JK=1, IKT
    DO JIJ=IIJB, IIJE
      ZWORK1(JIJ, JK) = ZGXMM_PTH(JIJ, JK)**2 + ZGYMM_PTH(JIJ, JK)**2
    END DO
  END DO

  CALL MZM_PHY(D,ZWORK1,ZWORK2)
  !
  IF (KRR /= 0) THEN                 ! moist 3D case
    CALL GX_M_M_PHY(D,OFLAT,PRM(:,:,1),PDXX,PDZZ,PDZX,ZGXMM_PRM)
    CALL GY_M_M_PHY(D,OFLAT,PRM(:,:,1),PDYY,PDZZ,PDZY,ZGYMM_PRM)

    DO JK=1, IKT
      DO JIJ=IIJB, IIJE
        ZWORK1(JIJ, JK) = ZGXMM_PRM(JIJ, JK)**2 + ZGYMM_PRM(JIJ, JK)**2
      END DO
    END DO

    CALL MZM_PHY(D,ZWORK1,ZWORK3)
    !

    DO JK=1, IKT
      DO JIJ=IIJB, IIJE
        ZWORK1(JIJ, JK) = ZGXMM_PRM(JIJ, JK) * ZGXMM_PTH(JIJ, JK) &
                                        + ZGYMM_PRM(JIJ, JK) * ZGYMM_PTH(JIJ, JK)
      END DO
    END DO

    CALL MZM_PHY(D,ZWORK1,ZWORK4)
    !

    DO JK=1, IKT
      DO JIJ=IIJB, IIJE
        PRED2TH3(JIJ, JK)= PREDTH1(JIJ, JK)**2 +  ( TURBN%XCTV*PBLL_O_E(JIJ, JK) &
                                          * PETHETA(JIJ, JK) )**2 * ZWORK2(JIJ, JK)      
    !
        PRED2R3(JIJ, JK)= PREDR1(JIJ, JK)**2 + (TURBN%XCTV*PBLL_O_E(JIJ, JK) &
                                         * PEMOIST(JIJ, JK))**2 * ZWORK3(JIJ, JK)
    !
        PRED2THR3(JIJ, JK)= PREDR1(JIJ, JK) * PREDTH1(JIJ, JK) + TURBN%XCTV**2 &
                                           * PBLL_O_E(JIJ, JK)**2 *   &
                                PEMOIST(JIJ, JK) * PETHETA(JIJ, JK) * ZWORK4(JIJ, JK)

      END DO
    END DO
!
    PRED2TH3(IIJB:IIJE,IKB)=PRED2TH3(IIJB:IIJE,IKB+IKL)
    PRED2R3(IIJB:IIJE,IKB)=PRED2R3(IIJB:IIJE,IKB+IKL)
    PRED2THR3(IIJB:IIJE,IKB)=PRED2THR3(IIJB:IIJE,IKB+IKL)

!
  ELSE                 ! dry 3D case in a 3D model

    DO JK=1, IKT
      DO JIJ=IIJB, IIJE
        PRED2TH3(JIJ, JK) = PREDTH1(JIJ, JK)**2 + TURBN%XCTV**2 &
                                            * PBLL_O_E(JIJ, JK)**2 * ZWORK2(JIJ, JK)
      END DO
    END DO
! 
    PRED2TH3(IIJB:IIJE,IKB)=PRED2TH3(IIJB:IIJE,IKB+IKL)
!
    PRED2R3(IIJB:IIJE,1:IKT) = 0.
!
    PRED2THR3(IIJB:IIJE,1:IKT) = 0.
!

  END IF
!
END IF   ! end of the if structure on the turbulence dimensionnality
!
!
!---------------------------------------------------------------------------
!
!           5. Prandtl numbers for scalars
!              ---------------------------
DO JSV=1,KSV
!
  IF(HTURBDIM=='1DIM') THEN
!        1D case

    DO JK=1, IKT
      DO JIJ=IIJB, IIJE
        PRED2THS3(JIJ, JK, JSV)  = PREDS1(JIJ, JK, JSV) * PREDTH1(JIJ, JK)
        IF (KRR /= 0) THEN
          PRED2RS3(JIJ, JK, JSV)   = PREDR1(JIJ, JK) *PREDS1(JIJ, JK, JSV)
        ELSE
          PRED2RS3(JIJ, JK, JSV)   = 0.
        END IF
      END DO
    END DO

!
  ELSE  IF (O2D) THEN ! 3D case in a 2D model
!
    IF (OOCEAN) THEN
!      !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
!      ZWORK1(IIJB:IIJE,1:IKT) = (CST%XG *CST%XALPHAOC * PLM(IIJB:IIJE,1:IKT) * PLEPS(IIJB:IIJE,1:IKT) &
!                                       / PTKEM(IIJB:IIJE,1:IKT))**2
!      !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
!      CALL MZM_PHY(D,ZWORK1,ZWORK2)
!      IF (KRR /= 0) THEN
!        !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
!        ZW1(IIJB:IIJE,1:IKT) = ZWORK2(IIJB:IIJE,1:IKT) * PETHETA(IIJB:IIJE,1:IKT)
!        !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
!      ELSE
!        ZW1 = ZWORK2
!      END IF
    ELSE
      !Compute only once and reuse in next JSV iterations
      IF ( JSV == 1 ) THEN

        DO JK=1, IKT
          DO JIJ=IIJB, IIJE
            ZWORK1(JIJ, JK) = (CST%XG / PTHVREF(JIJ, JK) * PLM(JIJ, JK) &
                                          * PLEPS(JIJ, JK) / PTKEM(JIJ, JK))**2
          END DO
        END DO

        CALL MZM_PHY(D,ZWORK1,ZW1)
      END IF
      !
      CALL GX_M_M_PHY(D,OFLAT,PSVM(:,:,JSV),PDXX,PDZZ,PDZX,ZGXMM_PSV)
      !Already computed CALL GX_M_M_PHY(D,OFLAT,PTHLM,PDXX,PDZZ,PDZX,ZGXMM_PTH)
      !

      DO JK=1, IKT
        DO JIJ=IIJB, IIJE
          ZWORK1(JIJ, JK) = ZGXMM_PSV(JIJ, JK) * ZGXMM_PTH(JIJ, JK)
        END DO
      END DO

      CALL MZM_PHY(D,ZWORK1,ZWORK2)
      !

      DO JK=1, IKT
        DO JIJ=IIJB, IIJE
          IF (KRR /= 0) THEN
            ZWORK1(JIJ, JK) = ZW1(JIJ, JK)*PETHETA(JIJ, JK)
          ELSE
            ZWORK1(JIJ, JK) = ZW1(JIJ, JK)
          END IF
          PRED2THS3(JIJ, JK, JSV) = PREDTH1(JIJ, JK) * PREDS1(JIJ, JK, JSV)   +        &
                             ZWORK1(JIJ, JK) * ZWORK2(JIJ, JK)
        END DO
      END DO

      !
      IF (KRR /= 0) THEN
        !Already computed CALL GX_M_M_PHY(D,OFLAT,PRM(:,:,1),PDXX,PDZZ,PDZX,ZGXMM_PRM)
        !

        DO JK=1, IKT
          DO JIJ=IIJB, IIJE
            ZWORK1(JIJ, JK) = ZGXMM_PSV(JIJ, JK) * ZGXMM_PRM(JIJ, JK)
          END DO
        END DO

        CALL MZM_PHY(D,ZWORK1,ZWORK3)
        !

        DO JK=1, IKT
          DO JIJ=IIJB, IIJE
            PRED2RS3(JIJ, JK, JSV) = PREDR1(JIJ, JK) * PREDS1(JIJ, JK, JSV)   +        &
                             ZW1(JIJ, JK) * PEMOIST(JIJ, JK) * ZWORK3(JIJ, JK)
          END DO
        END DO

      ELSE

        DO JK=1, IKT
          DO JIJ=IIJB, IIJE
            PRED2RS3(JIJ, JK, JSV) = 0.
          END DO
        END DO

      END IF
    END IF
!
  ELSE ! 3D case in a 3D model
!
    IF (OOCEAN) THEN
!      !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
!      ZWORK1(IIJB:IIJE,1:IKT) = (CST%XG *CST%XALPHAOC * PLM(IIJB:IIJE,1:IKT) * PLEPS(IIJB:IIJE,1:IKT) &
!                                       / PTKEM(IIJB:IIJE,1:IKT))**2
!      !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
!      CALL MZM_PHY(D,ZWORK1,ZWORK2)
!      IF (KRR /= 0) THEN
!        !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
!        ZW1(IIJB:IIJE,1:IKT) = ZWORK2(IIJB:IIJE,1:IKT) * PETHETA(IIJB:IIJE,1:IKT)
!        !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
!      ELSE
!        ZW1 = ZWORK2
!      END IF
    ELSE
      !Compute only once and reuse in next JSV iterations
      IF ( JSV == 1 ) THEN

        DO JK=1, IKT
          DO JIJ=IIJB, IIJE
            ZWORK1(JIJ, JK) = (CST%XG / PTHVREF(JIJ, JK) * PLM(JIJ, JK) &
                                          * PLEPS(JIJ, JK) / PTKEM(JIJ, JK))**2
          END DO
        END DO

        CALL MZM_PHY(D,ZWORK1,ZW1)
      END IF
      !
      CALL GX_M_M_PHY(D,OFLAT,PSVM(:,:,JSV),PDXX,PDZZ,PDZX,ZGXMM_PSV)
      !Already computed CALL GX_M_M_PHY(D,OFLAT,PTHLM,PDXX,PDZZ,PDZX,ZGXMM_PTH)
      CALL GY_M_M_PHY(D,OFLAT,PSVM(:,:,JSV),PDYY,PDZZ,PDZY,ZGYMM_PSV)
      !Already computed CALL GY_M_M_PHY(D,OFLAT,PTHLM,PDYY,PDZZ,PDZY,ZGYMM_PTH)
      !

      DO JK=1, IKT
        DO JIJ=IIJB, IIJE
          ZWORK1(JIJ, JK) = ZGXMM_PSV(JIJ, JK) * ZGXMM_PTH(JIJ, JK) &
                                          + ZGYMM_PSV(JIJ, JK) * ZGYMM_PTH(JIJ, JK)
        END DO
      END DO

      CALL MZM_PHY(D,ZWORK1,ZWORK2)
      !
      IF (KRR /= 0) THEN
        !Already computed CALL GX_M_M_PHY(D,OFLAT,PRM(:,:,1),PDXX,PDZZ,PDZX,ZGXMM_PRM)
        !Already computed CALL GY_M_M_PHY(D,OFLAT,PRM(:,:,1),PDYY,PDZZ,PDZY,ZGYMM_PRM)
        !

        DO JK=1, IKT
          DO JIJ=IIJB, IIJE
            ZWORK1(JIJ, JK) = ZGXMM_PSV(JIJ, JK) * ZGXMM_PRM(JIJ, JK) &
                                          + ZGYMM_PSV(JIJ, JK) * ZGYMM_PRM(JIJ, JK)
          END DO
        END DO

        CALL MZM_PHY(D,ZWORK1,ZWORK3)
        !

        DO JK=1, IKT
          DO JIJ=IIJB, IIJE
            ZWORK1(JIJ, JK) = ZW1(JIJ, JK)*PETHETA(JIJ, JK)
          END DO
        END DO

      ELSE

        DO JK=1, IKT
          DO JIJ=IIJB, IIJE
            ZWORK1(JIJ, JK) = ZW1(JIJ, JK)
          END DO
        END DO

      END IF

      DO JK=1, IKT
        DO JIJ=IIJB, IIJE
          PRED2THS3(JIJ, JK, JSV) = PREDTH1(JIJ, JK) * PREDS1(JIJ, JK, JSV)   +        &
                             ZWORK1(JIJ, JK)*ZWORK2(JIJ, JK)
        END DO
      END DO

      IF (KRR /= 0) THEN

        DO JK=1, IKT
          DO JIJ=IIJB, IIJE
            PRED2RS3(JIJ, JK, JSV) = PREDR1(JIJ, JK) * PREDS1(JIJ, JK, JSV)   +        &
                             ZW1(JIJ, JK) * PEMOIST(JIJ, JK) * ZWORK3(JIJ, JK)
          END DO
        END DO

      ELSE

        DO JK=1, IKT
          DO JIJ=IIJB, IIJE
            PRED2RS3(JIJ, JK, JSV) = 0.
          END DO
        END DO

      END IF
  END IF 
!
  END IF ! end of HTURBDIM if-block
!
END DO
!
!---------------------------------------------------------------------------
!
!*          6. SAVES THE REDELSPERGER NUMBERS
!              ------------------------------
!
IF ( OTURB_DIAG .AND. TPFILE%LOPENED ) THEN
  !
  ! stores the RED_TH1
  TZFIELD = TFIELDMETADATA(       &
    CMNHNAME   = 'RED_TH1',       &
    CSTDNAME   = '',              &
    CLONGNAME  = 'RED_TH1',       &
    CUNITS     = '1',             &
    CDIR       = 'XY',            &
    CCOMMENT   = 'X_Y_Z_RED_TH1', &
    NGRID      = 4,               &
    NTYPE      = TYPEREAL,        &
    NDIMS      = 3,               &
    LTIMEDEP   = .TRUE.           )
  CALL IO_FIELD_WRITE_PHY(D,TPFILE,TZFIELD,PREDTH1)
  !
  ! stores the RED_R1
  TZFIELD = TFIELDMETADATA(      &
    CMNHNAME   = 'RED_R1',       &
    CSTDNAME   = '',             &
    CLONGNAME  = 'RED_R1',       &
    CUNITS     = '1',            &
    CDIR       = 'XY',           &
    CCOMMENT   = 'X_Y_Z_RED_R1', &
    NGRID      = 4,              &
    NTYPE      = TYPEREAL,       &
    NDIMS      = 3,              &
    LTIMEDEP   = .TRUE.          )
  CALL IO_FIELD_WRITE_PHY(D,TPFILE,TZFIELD,PREDR1)
  !
  ! stores the RED2_TH3
  TZFIELD = TFIELDMETADATA(        &
    CMNHNAME   = 'RED2_TH3',       &
    CSTDNAME   = '',               &
    CLONGNAME  = 'RED2_TH3',       &
    CUNITS     = '1',              &
    CDIR       = 'XY',             &
    CCOMMENT   = 'X_Y_Z_RED2_TH3', &
    NGRID      = 4,                &
    NTYPE      = TYPEREAL,         &
    NDIMS      = 3,                &
    LTIMEDEP   = .TRUE.            )
  CALL IO_FIELD_WRITE_PHY(D,TPFILE,TZFIELD,PRED2TH3)
  !
  ! stores the RED2_R3
  TZFIELD = TFIELDMETADATA(       &
    CMNHNAME   = 'RED2_R3',       &
    CSTDNAME   = '',              &
    CLONGNAME  = 'RED2_R3',       &
    CUNITS     = '1',             &
    CDIR       = 'XY',            &
    CCOMMENT   = 'X_Y_Z_RED2_R3', &
    NGRID      = 4,               &
    NTYPE      = TYPEREAL,        &
    NDIMS      = 3,               &
    LTIMEDEP   = .TRUE.           )
  CALL IO_FIELD_WRITE_PHY(D,TPFILE,TZFIELD,PRED2R3)
  !
  ! stores the RED2_THR3
  TZFIELD = TFIELDMETADATA(         &
    CMNHNAME   = 'RED2_THR3',       &
    CSTDNAME   = '',                &
    CLONGNAME  = 'RED2_THR3',       &
    CUNITS     = '1',               &
    CDIR       = 'XY',              &
    CCOMMENT   = 'X_Y_Z_RED2_THR3', &
    NGRID      = 4,                 &
    NTYPE      = TYPEREAL,          &
    NDIMS      = 3,                 &
    LTIMEDEP   = .TRUE.             )
  CALL IO_FIELD_WRITE_PHY(D,TPFILE,TZFIELD,PRED2THR3)
  !
END IF
!
!---------------------------------------------------------------------------
ENDIF ! (Done only if OHARAT is FALSE)
!
IF (LHOOK) CALL DR_HOOK('PRANDTL',1,ZHOOK_HANDLE)
END SUBROUTINE PRANDTL
!
SUBROUTINE SMOOTH_TURB_FUNCT(D,CSTURB,TURBN,PPHI3,PF_LIM,PF)
!
TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
TYPE(TURB_t),                    INTENT(IN)    :: TURBN
TYPE(DIMPHYEX_t),       INTENT(IN)   :: D
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PPHI3   ! Phi3
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PF_LIM  ! Value of F when Phi3 is
!                                                ! larger than Phi_lim
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(INOUT) :: PF      ! function F to smooth
!
REAL, DIMENSION(D%NIJT,D%NKT) :: ZCOEF
INTEGER :: JIJ,JK, IIJB,IIJE,IKT
!
!* adds a artificial correction to smooth the function near the discontinuity
!  point at Phi3 = Phi_lim
!  This smoothing is applied between 0.9*phi_lim (=2.7) and Phi_lim (=3)
!   Note that in the Boundary layer, phi is usually between 0.8 and 1
!
IIJE=D%NIJE
IIJB=D%NIJB
IKT=D%NKT
!
IF(TURBN%LSMOOTH_PRANDTL) THEN

  DO JK=1, IKT
    DO JIJ=IIJB, IIJE
      ZCOEF(JIJ, JK) = MAX(MIN((  10.*(1.-PPHI3(JIJ, JK)/CSTURB%XPHI_LIM)) ,1.), 0.) 
      !
      PF(JIJ, JK) =     ZCOEF(JIJ, JK)   * PF(JIJ, JK)    &
                + (1.-ZCOEF(JIJ, JK))  * PF_LIM(JIJ, JK)
    END DO
  END DO

ENDIF
!
END SUBROUTINE SMOOTH_TURB_FUNCT
!----------------------------------------------------------------------------
SUBROUTINE PHI3(D,CSTURB,TURBN,PREDTH1,PREDR1,PRED2TH3,PRED2R3,PRED2THR3,HTURBDIM,OUSERV,PPHI3)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  TYPE(TURB_t),                    INTENT(IN)    :: TURBN
  TYPE(DIMPHYEX_t),                   INTENT(IN)   :: D
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PRED2TH3
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PRED2R3
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PRED2THR3
  CHARACTER(LEN=4),       INTENT(IN) :: HTURBDIM  ! 1DIM or 3DIM turb. scheme
  LOGICAL,                INTENT(IN) :: OUSERV    ! flag to use vapor
  REAL, DIMENSION(D%NIJT,D%NKT),INTENT(OUT) :: PPHI3
!
  REAL, DIMENSION(D%NIJT,D%NKT) :: ZW1, ZW2
  INTEGER :: IKB, IKE, JIJ,JK, IIJB,IIJE, IKT
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:PHI3',0,ZHOOK_HANDLE)
IKB=D%NKTB
IKE=D%NKTE
IIJE=D%NIJE
IIJB=D%NIJB
IKT=D%NKT
!

IF (HTURBDIM=='3DIM') THEN
        !* 3DIM case
  DO JK=1, IKT
    DO JIJ=IIJB, IIJE
      IF (OUSERV) THEN
        ZW1(JIJ, JK) = 1. + 1.5* (PREDTH1(JIJ, JK)+PREDR1(JIJ, JK)) +      &
                       ( 0.5 * (PREDTH1(JIJ, JK)**2+PREDR1(JIJ, JK)**2)  &
                             + PREDTH1(JIJ, JK) * PREDR1(JIJ, JK)        &
                       )

        ZW2(JIJ, JK) = 0.5 * (PRED2TH3(JIJ, JK)-PRED2R3(JIJ, JK))

        PPHI3(JIJ, JK)= 1. -                                          &
        ( ( (1.+PREDR1(JIJ, JK)) *                                   &
            (PRED2THR3(JIJ, JK) + PRED2TH3(JIJ, JK)) / PREDTH1(JIJ, JK)  &
          ) + ZW2(JIJ, JK)                                           &
        ) / ZW1(JIJ, JK)
      ELSE
        ZW1(JIJ, JK) = 1. + 1.5* PREDTH1(JIJ, JK) + &
                     0.5* PREDTH1(JIJ, JK)**2

        ZW2(JIJ, JK) = 0.5* PRED2TH3(JIJ, JK)

        PPHI3(JIJ, JK)= 1. -                                       &
                (PRED2TH3(JIJ, JK) / PREDTH1(JIJ, JK) + ZW2(JIJ, JK)) &
                / ZW1(JIJ, JK)
      END IF
    END DO
  END DO

  DO JK=1, IKT
    DO JIJ=IIJB, IIJE
      IF ( PPHI3(JIJ, JK) <= 0. .OR. PPHI3(JIJ, JK) > CSTURB%XPHI_LIM ) THEN
        PPHI3(JIJ, JK) = CSTURB%XPHI_LIM
      END IF
    END DO
  END DO
ELSE
        !* 1DIM case
  DO JK=1, IKT
    DO JIJ=IIJB, IIJE
      IF (OUSERV) THEN
        PPHI3(JIJ, JK)= 1./(1.+PREDTH1(JIJ, JK)+PREDR1(JIJ, JK))
      ELSE
        PPHI3(JIJ, JK)= 1./(1.+PREDTH1(JIJ, JK))
      END IF
    END DO
  END DO
END IF
!
PPHI3(IIJB:IIJE,IKB-1)=PPHI3(IIJB:IIJE,IKB)
PPHI3(IIJB:IIJE,IKE+1)=PPHI3(IIJB:IIJE,IKE)

!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:PHI3',1,ZHOOK_HANDLE)
END SUBROUTINE PHI3
!----------------------------------------------------------------------------
SUBROUTINE PSI_SV(D,CSTURB,TURBN,KSV,PREDTH1,PREDR1,PREDS1,PRED2THS,PRED2RS,PPHI3,PPSI3,PPSI_SV)
  TYPE(CSTURB_t),                  INTENT(IN)      :: CSTURB
  TYPE(TURB_t),                    INTENT(IN)    :: TURBN
  TYPE(DIMPHYEX_t),                INTENT(IN)      :: D
  INTEGER,                         INTENT(IN)      :: KSV
  REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIJT,D%NKT,KSV), INTENT(IN) :: PREDS1
  REAL, DIMENSION(D%NIJT,D%NKT,KSV), INTENT(IN) :: PRED2THS
  REAL, DIMENSION(D%NIJT,D%NKT,KSV), INTENT(IN) :: PRED2RS
  REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN) :: PPHI3
  REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN) :: PPSI3
  REAL, DIMENSION(D%NIJT,D%NKT,KSV),INTENT(OUT) :: PPSI_SV
!
  INTEGER :: IKB, IKE, IIJB,IIJE, IKT
  INTEGER :: JSV,JIJ,JK
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:PSI_SV',0,ZHOOK_HANDLE)
IKB=D%NKTB
IKE=D%NKTE
IIJE=D%NIJE
IIJB=D%NIJB
IKT=D%NKT
!

DO JSV=1,KSV
  DO JK=1, IKT
    DO JIJ=IIJB, IIJE
      PPSI_SV(JIJ, JK, JSV) = ( 1.                                             &
        - (CSTURB%XCPR3+CSTURB%XCPR5) * &
        (PRED2THS(JIJ, JK, JSV)/PREDS1(JIJ, JK, JSV)-PREDTH1(JIJ, JK)) &
        - (CSTURB%XCPR4+CSTURB%XCPR5) * &
        (PRED2RS(JIJ, JK, JSV)/PREDS1(JIJ, JK, JSV)-PREDR1(JIJ, JK)) &
        - CSTURB%XCPR3 * &
        PREDTH1(JIJ, JK) * PPHI3(JIJ, JK) &
        - CSTURB%XCPR4 * PREDR1(JIJ, JK) * PPSI3(JIJ, JK)                 &
                    ) / (  1. + CSTURB%XCPR5 * ( PREDTH1(JIJ, JK) + PREDR1(JIJ, JK) ) )           
  
    END DO
  END DO
!        control of the PSI_SV positivity
  DO JK=1, IKT
    DO JIJ=IIJB, IIJE
      IF ( (PPSI_SV(JIJ, JK, JSV) <=0.).AND. (PREDTH1(JIJ, JK)+PREDR1(JIJ, JK))<=0.) THEN
        PPSI_SV(JIJ, JK, JSV)=CSTURB%XPHI_LIM
      END IF
    END DO
  END DO

  DO JK=1, IKT
    DO JIJ=IIJB, IIJE
      PPSI_SV(JIJ, JK, JSV) = MAX( 1.E-4, MIN(CSTURB%XPHI_LIM,PPSI_SV(JIJ, JK, JSV)) )
    END DO
  END DO
!
  PPSI_SV(IIJB:IIJE,IKB-1,JSV)=PPSI_SV(IIJB:IIJE,IKB,JSV)
  PPSI_SV(IIJB:IIJE,IKE+1,JSV)=PPSI_SV(IIJB:IIJE,IKE,JSV)
END DO

!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:PSI_SV',1,ZHOOK_HANDLE)
END SUBROUTINE PSI_SV
!----------------------------------------------------------------------------
SUBROUTINE D_PHI3DTDZ_O_DDTDZ(D,CSTURB,TURBN,PPHI3,PREDTH1,PREDR1,PRED2TH3,PRED2THR3,HTURBDIM,OUSERV,PD_PHI3DTDZ_O_DDTDZ)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  TYPE(TURB_t),                    INTENT(IN)    :: TURBN
  TYPE(DIMPHYEX_t),                   INTENT(IN)   :: D
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PPHI3
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PRED2TH3
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PRED2THR3
  CHARACTER(LEN=4),       INTENT(IN) :: HTURBDIM  ! 1DIM or 3DIM turb. scheme
  LOGICAL,                INTENT(IN) :: OUSERV    ! flag to use vapor
  REAL, DIMENSION(D%NIJT,D%NKT),INTENT(OUT) :: PD_PHI3DTDZ_O_DDTDZ
  INTEGER :: IKB, IKE,JIJ,JK, IIJB,IIJE,IKT
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_PHI3DTDZ_O_DDTDZ',0,ZHOOK_HANDLE)
IKB=D%NKTB
IKE=D%NKTE
IIJE=D%NIJE
IIJB=D%NIJB
IKT=D%NKT
!

IF (HTURBDIM=='3DIM') THEN
        !* 3DIM case
  IF (OUSERV) THEN
   DO JK=1, IKT
      DO JIJ=IIJB, IIJE
        IF (PPHI3(JIJ, JK)<=CSTURB%XPHI_LIM) THEN
          PD_PHI3DTDZ_O_DDTDZ(JIJ, JK) = PPHI3(JIJ, JK)                       &
              * (1. - PREDTH1(JIJ, JK) * (3./2.+PREDTH1(JIJ, JK)+PREDR1(JIJ, JK))   &
                   /((1.+PREDTH1(JIJ, JK)+PREDR1(JIJ, JK)) &
                   *(1.+1./2.*(PREDTH1(JIJ, JK)+PREDR1(JIJ, JK))))) &
              + (1.+PREDR1(JIJ, JK))*(PRED2THR3(JIJ, JK)+PRED2TH3(JIJ, JK))         &
                   / (PREDTH1(JIJ, JK)*(1.+PREDTH1(JIJ, JK)+PREDR1(JIJ, JK))* &
                     (1.+1./2.*(PREDTH1(JIJ, JK)+PREDR1(JIJ, JK)))) &
              - (1./2.*PREDTH1(JIJ, JK)+PREDR1(JIJ, JK) &
              * (1.+PREDTH1(JIJ, JK)+PREDR1(JIJ, JK)))           &
                   / ((1.+PREDTH1(JIJ, JK)+PREDR1(JIJ, JK))&
                   *(1.+1./2.*(PREDTH1(JIJ, JK)+PREDR1(JIJ, JK))))
        ELSE
          PD_PHI3DTDZ_O_DDTDZ(JIJ, JK) = PPHI3(JIJ, JK)
        END IF
     END DO
    END DO

!
  ELSE
   DO JK=1, IKT
      DO JIJ=IIJB, IIJE
        IF (PPHI3(JIJ, JK)<=CSTURB%XPHI_LIM) THEN
        PD_PHI3DTDZ_O_DDTDZ(JIJ, JK) = PPHI3(JIJ, JK)             &
              * (1. - PREDTH1(JIJ, JK) * (3./2.+PREDTH1(JIJ, JK))      &
                   /((1.+PREDTH1(JIJ, JK))*(1.+1./2.*PREDTH1(JIJ, JK))))        &
              + PRED2TH3(JIJ, JK) &
              / (PREDTH1(JIJ, JK)*(1.+PREDTH1(JIJ, JK))*(1.+1./2.*PREDTH1(JIJ, JK))) &
              - 1./2.*PREDTH1(JIJ, JK) &
              / ((1.+PREDTH1(JIJ, JK))*(1.+1./2.*PREDTH1(JIJ, JK)))
        ELSE
          PD_PHI3DTDZ_O_DDTDZ(JIJ, JK) = PPHI3(JIJ, JK)
        END IF
     END DO
    END DO
!
  END IF
ELSE
        !* 1DIM case
DO JK=1,IKT
    DO JIJ=IIJB,IIJE
      IF ( ABS(PPHI3(JIJ,JK)-CSTURB%XPHI_LIM) < 1.E-12 ) THEN
         PD_PHI3DTDZ_O_DDTDZ(JIJ,JK)=PPHI3(JIJ,JK)*&
&       (1. - PREDTH1(JIJ,JK)*PPHI3(JIJ,JK))
      ELSE
         PD_PHI3DTDZ_O_DDTDZ(JIJ,JK)=PPHI3(JIJ,JK)
      ENDIF
  ENDDO
ENDDO
END IF

!
!* smoothing
CALL SMOOTH_TURB_FUNCT(D,CSTURB,TURBN,PPHI3,PPHI3,PD_PHI3DTDZ_O_DDTDZ)
!

PD_PHI3DTDZ_O_DDTDZ(IIJB:IIJE,IKB-1)=PD_PHI3DTDZ_O_DDTDZ(IIJB:IIJE,IKB)
PD_PHI3DTDZ_O_DDTDZ(IIJB:IIJE,IKE+1)=PD_PHI3DTDZ_O_DDTDZ(IIJB:IIJE,IKE)

!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_PHI3DTDZ_O_DDTDZ',1,ZHOOK_HANDLE)
END SUBROUTINE D_PHI3DTDZ_O_DDTDZ
!----------------------------------------------------------------------------
SUBROUTINE D_PHI3DRDZ_O_DDRDZ(D,CSTURB,TURBN,PPHI3,PREDTH1,PREDR1,PRED2TH3,PRED2THR3,HTURBDIM,OUSERV,PD_PHI3DRDZ_O_DDRDZ)
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  TYPE(TURB_t),                    INTENT(IN)    :: TURBN
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PPHI3
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PRED2TH3
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PRED2THR3
  CHARACTER(LEN=4),       INTENT(IN) :: HTURBDIM  ! 1DIM or 3DIM turb. scheme
  LOGICAL,                INTENT(IN) :: OUSERV    ! flag to use vapor
  REAL, DIMENSION(D%NIJT,D%NKT),INTENT(OUT) :: PD_PHI3DRDZ_O_DDRDZ
  INTEGER :: IKB, IKE, JIJ,JK, IIJB,IIJE,IKT
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_PHI3DRDZ_O_DDRDZ',0,ZHOOK_HANDLE)
IKB=D%NKTB
IKE=D%NKTE
IIJE=D%NIJE
IIJB=D%NIJB
IKT=D%NKT
!
!

IF (HTURBDIM=='3DIM') THEN
        !* 3DIM case
  IF (OUSERV) THEN
   DO JK=1, IKT
      DO JIJ=IIJB, IIJE
        IF (PPHI3(JIJ, JK)<=CSTURB%XPHI_LIM) THEN
          PD_PHI3DRDZ_O_DDRDZ(JIJ, JK) = PPHI3(JIJ, JK) &
          * (1.-PREDR1(JIJ, JK)*(3./2.+PREDTH1(JIJ, JK)+PREDR1(JIJ, JK)) &
          / ((1.+PREDTH1(JIJ, JK)+PREDR1(JIJ, JK)) & 
          *(1.+1./2.*(PREDTH1(JIJ, JK)+PREDR1(JIJ, JK)))))  &
          - PREDR1(JIJ, JK) &
          * (PRED2THR3(JIJ, JK)+PRED2TH3(JIJ, JK)) / (PREDTH1(JIJ, JK)         &
          * (1.+PREDTH1(JIJ, JK)+PREDR1(JIJ, JK))*&
          (1.+1./2.*(PREDTH1(JIJ, JK)+PREDR1(JIJ, JK))))    &
          + PREDR1(JIJ, JK) * (1./2.+PREDTH1(JIJ, JK)+PREDR1(JIJ, JK))         &
          / ((1.+PREDTH1(JIJ, JK)+PREDR1(JIJ, JK))&
          *(1.+1./2.*(PREDTH1(JIJ, JK)+PREDR1(JIJ, JK))))
        ELSE
          PD_PHI3DRDZ_O_DDRDZ(JIJ, JK) = PPHI3(JIJ, JK)
        END IF
     END DO
    END DO
  ELSE
    PD_PHI3DRDZ_O_DDRDZ(IIJB:IIJE,1:IKT) = PPHI3(IIJB:IIJE,1:IKT)
  END IF
ELSE
        !* 1DIM case
    DO JK=1, IKT
      DO JIJ=IIJB, IIJE
        IF (PPHI3(JIJ, JK)<=CSTURB%XPHI_LIM) THEN
        PD_PHI3DRDZ_O_DDRDZ(JIJ, JK) = PPHI3(JIJ, JK)                           &
              * (1. - PREDR1(JIJ, JK)*PPHI3(JIJ, JK))
      ELSE
        PD_PHI3DRDZ_O_DDRDZ(JIJ, JK) = PPHI3(JIJ, JK)
      END IF
    END DO
    END DO
END IF

!
!* smoothing
CALL SMOOTH_TURB_FUNCT(D,CSTURB,TURBN,PPHI3,PPHI3,PD_PHI3DRDZ_O_DDRDZ)
!

PD_PHI3DRDZ_O_DDRDZ(IIJB:IIJE,IKB-1)=PD_PHI3DRDZ_O_DDRDZ(IIJB:IIJE,IKB)
PD_PHI3DRDZ_O_DDRDZ(IIJB:IIJE,IKE+1)=PD_PHI3DRDZ_O_DDRDZ(IIJB:IIJE,IKE)

!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_PHI3DRDZ_O_DDRDZ',1,ZHOOK_HANDLE)
END SUBROUTINE D_PHI3DRDZ_O_DDRDZ
!----------------------------------------------------------------------------
SUBROUTINE D_PHI3DTDZ2_O_DDTDZ(D,CSTURB,TURBN,PPHI3,PREDTH1,PREDR1,PRED2TH3,PRED2THR3,PDTDZ,HTURBDIM,OUSERV,PD_PHI3DTDZ2_O_DDTDZ)
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  TYPE(TURB_t),                    INTENT(IN)    :: TURBN
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PPHI3
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PRED2TH3
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PRED2THR3
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PDTDZ
  CHARACTER(LEN=4),       INTENT(IN) :: HTURBDIM  ! 1DIM or 3DIM turb. scheme
  LOGICAL,                INTENT(IN) :: OUSERV    ! flag to use vapor
  REAL, DIMENSION(D%NIJT,D%NKT),INTENT(OUT) :: PD_PHI3DTDZ2_O_DDTDZ
  REAL, DIMENSION(D%NIJT,D%NKT) :: ZWORK1 ! working array
  INTEGER :: IKB, IKE, JIJ,JK, IIJB,IIJE,IKT
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_PHI3DTDZ2_O_DDTDZ',0,ZHOOK_HANDLE)
IKB=D%NKTB
IKE=D%NKTE
IIJE=D%NIJE
IIJB=D%NIJB
IKT=D%NKT
!
!
IF (HTURBDIM=='3DIM') THEN
   ! by derivation of (phi3 dtdz) * dtdz according to dtdz we obtain:
   CALL D_PHI3DTDZ_O_DDTDZ(D,CSTURB,TURBN,PPHI3,PREDTH1,PREDR1,PRED2TH3,PRED2THR3,HTURBDIM,OUSERV,ZWORK1)

   DO JK=1, IKT
     DO JIJ=IIJB, IIJE
       PD_PHI3DTDZ2_O_DDTDZ(JIJ, JK) = PDTDZ(JIJ, JK) &
       * (PPHI3(JIJ, JK) +  ZWORK1(JIJ, JK))
     END DO
   END DO

ELSE
        !* 1DIM case

    DO JK=1, IKT
      DO JIJ=IIJB, IIJE
        IF (PPHI3(JIJ, JK)<=CSTURB%XPHI_LIM) THEN
          PD_PHI3DTDZ2_O_DDTDZ(JIJ, JK) = PPHI3(JIJ, JK)*PDTDZ(JIJ, JK)             &
              * (2. - PREDTH1(JIJ, JK)*PPHI3(JIJ, JK))
        ELSE
          PD_PHI3DTDZ2_O_DDTDZ(JIJ, JK) = PPHI3(JIJ, JK) * 2. * PDTDZ(JIJ, JK)
        END IF
      END DO
    END DO

END IF
!
!* smoothing
DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZWORK1(JIJ, JK)=PPHI3(JIJ, JK)*2.*PDTDZ(JIJ, JK)
  END DO
END DO
CALL SMOOTH_TURB_FUNCT(D,CSTURB,TURBN,PPHI3,ZWORK1,PD_PHI3DTDZ2_O_DDTDZ)
!
!

PD_PHI3DTDZ2_O_DDTDZ(IIJB:IIJE,IKB-1)=PD_PHI3DTDZ2_O_DDTDZ(IIJB:IIJE,IKB)
PD_PHI3DTDZ2_O_DDTDZ(IIJB:IIJE,IKE+1)=PD_PHI3DTDZ2_O_DDTDZ(IIJB:IIJE,IKE)

!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_PHI3DTDZ2_O_DDTDZ',1,ZHOOK_HANDLE)
END SUBROUTINE D_PHI3DTDZ2_O_DDTDZ
!----------------------------------------------------------------------------
SUBROUTINE M3_WTH_WTH2(D,CSTURB,TURBN,PREDTH1,PREDR1,PD,PBLL_O_E,PETHETA,PM3_WTH_WTH2)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  TYPE(TURB_t),                    INTENT(IN)    :: TURBN
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PETHETA
  REAL, DIMENSION(D%NIJT,D%NKT),INTENT(OUT) :: PM3_WTH_WTH2
  INTEGER :: IKB, IKE, JIJ,JK, IIJB,IIJE,IKT
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_WTH_WTH2',0,ZHOOK_HANDLE)
IKB=D%NKTB
IKE=D%NKTE
IIJE=D%NIJE
IIJB=D%NIJB
IKT=D%NKT
!

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    PM3_WTH_WTH2(JIJ, JK) = TURBN%XCSHF*PBLL_O_E(JIJ, JK)&
                       * PETHETA(JIJ, JK)*0.5/CSTURB%XCTD        &
                       * (1.+0.5*PREDTH1(JIJ, JK)+PREDR1(JIJ, JK)) / PD(JIJ, JK)
  END DO
END DO
PM3_WTH_WTH2(IIJB:IIJE,IKB-1)=PM3_WTH_WTH2(IIJB:IIJE,IKB)
PM3_WTH_WTH2(IIJB:IIJE,IKE+1)=PM3_WTH_WTH2(IIJB:IIJE,IKE)

!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_WTH_WTH2',1,ZHOOK_HANDLE)
END SUBROUTINE M3_WTH_WTH2
!----------------------------------------------------------------------------
SUBROUTINE D_M3_WTH_WTH2_O_DDTDZ(D,CSTURB,TURBN,PM3_WTH_WTH2,PREDTH1,PREDR1,PD,PBLL_O_E,PETHETA,PD_M3_WTH_WTH2_O_DDTDZ)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  TYPE(TURB_t),                    INTENT(IN)    :: TURBN
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PM3_WTH_WTH2
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PETHETA
  REAL, DIMENSION(D%NIJT,D%NKT),INTENT(OUT) :: PD_M3_WTH_WTH2_O_DDTDZ
  INTEGER :: IKB, IKE, JIJ,JK, IIJB,IIJE,IKT
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_WTH_WTH2_O_DDTDZ',0,ZHOOK_HANDLE)
IKB=D%NKTB
IKE=D%NKTE
IIJE=D%NIJE
IIJB=D%NIJB
IKT=D%NKT
!

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    PD_M3_WTH_WTH2_O_DDTDZ(JIJ, JK) = &
    (0.5*TURBN%XCSHF*PBLL_O_E(JIJ, JK)*PETHETA(JIJ, JK)*0.5/CSTURB%XCTD/PD(JIJ, JK) &
    - PM3_WTH_WTH2(JIJ, JK)/PD(JIJ, JK)&
    *(1.5+PREDTH1(JIJ, JK)+PREDR1(JIJ, JK))  )&
    * PBLL_O_E(JIJ, JK) * PETHETA(JIJ, JK) * TURBN%XCTV
  END DO
END DO
!
PD_M3_WTH_WTH2_O_DDTDZ(IIJB:IIJE,IKB-1)=PD_M3_WTH_WTH2_O_DDTDZ(IIJB:IIJE,IKB)
PD_M3_WTH_WTH2_O_DDTDZ(IIJB:IIJE,IKE+1)=PD_M3_WTH_WTH2_O_DDTDZ(IIJB:IIJE,IKE)

!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_WTH_WTH2_O_DDTDZ',1,ZHOOK_HANDLE)
END SUBROUTINE D_M3_WTH_WTH2_O_DDTDZ
!----------------------------------------------------------------------------
SUBROUTINE M3_WTH_W2TH(D,CSTURB,TURBN,PREDTH1,PREDR1,PD,PKEFF,PTKE,PM3_WTH_W2TH)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  TYPE(TURB_t),                    INTENT(IN)    :: TURBN
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PKEFF
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PTKE
  REAL, DIMENSION(D%NIJT,D%NKT),INTENT(OUT) :: PM3_WTH_W2TH
  REAL, DIMENSION(D%NIJT,D%NKT) :: ZWORK1 ! working array
  INTEGER :: IKB, IKE, JIJ,JK, IIJB,IIJE,IKT
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_WTH_W2TH',0,ZHOOK_HANDLE)
IKB=D%NKTB
IKE=D%NKTE
IIJE=D%NIJE
IIJB=D%NIJB
IKT=D%NKT
!
CALL MZM_PHY(D,PTKE,ZWORK1)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    PM3_WTH_W2TH(JIJ, JK) = TURBN%XCSHF*PKEFF(JIJ, JK)*1.5/ZWORK1(JIJ, JK) &
      * (1. - 0.5*PREDR1(JIJ, JK)*(1.+PREDR1(JIJ, JK))/PD(JIJ, JK) ) &
      / (1.+PREDTH1(JIJ, JK))
  END DO
END DO
!
PM3_WTH_W2TH(IIJB:IIJE,IKB-1)=PM3_WTH_W2TH(IIJB:IIJE,IKB)
PM3_WTH_W2TH(IIJB:IIJE,IKE+1)=PM3_WTH_W2TH(IIJB:IIJE,IKE)

!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_WTH_W2TH',1,ZHOOK_HANDLE)
END SUBROUTINE M3_WTH_W2TH
!----------------------------------------------------------------------------
SUBROUTINE D_M3_WTH_W2TH_O_DDTDZ(D,CSTURB,TURBN,PREDTH1,PREDR1,PD,PBLL_O_E,PETHETA,PKEFF,PTKE,PD_M3_WTH_W2TH_O_DDTDZ)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  TYPE(TURB_t),                    INTENT(IN)    :: TURBN
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PETHETA
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PKEFF
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PTKE
  REAL, DIMENSION(D%NIJT,D%NKT),INTENT(OUT) :: PD_M3_WTH_W2TH_O_DDTDZ
  REAL, DIMENSION(D%NIJT,D%NKT) :: ZWORK1 ! working array
  INTEGER :: IKB, IKE, JIJ,JK, IIJB,IIJE,IKT
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_WTH_W2TH_O_DDTDZ',0,ZHOOK_HANDLE)
IKB=D%NKTB
IKE=D%NKTE
IIJE=D%NIJE
IIJB=D%NIJB
IKT=D%NKT
!
CALL MZM_PHY(D,PTKE,ZWORK1)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    PD_M3_WTH_W2TH_O_DDTDZ(JIJ, JK) = &
     - TURBN%XCSHF*PKEFF(JIJ, JK)*1.5/ZWORK1(JIJ, JK)/(1.+PREDTH1(JIJ, JK))**2 &
     * TURBN%XCTV*PBLL_O_E(JIJ, JK)*PETHETA(JIJ, JK)  &
     * (1. - 0.5*PREDR1(JIJ, JK)*(1.+PREDR1(JIJ, JK))/PD(JIJ, JK)* &
       ( 1.+(1.+PREDTH1(JIJ, JK))*(1.5+PREDR1(JIJ, JK)+PREDTH1(JIJ, JK))&
       /PD(JIJ, JK)) )
  END DO
END DO
!
PD_M3_WTH_W2TH_O_DDTDZ(IIJB:IIJE,IKB-1)=PD_M3_WTH_W2TH_O_DDTDZ(IIJB:IIJE,IKB)
PD_M3_WTH_W2TH_O_DDTDZ(IIJB:IIJE,IKE+1)=PD_M3_WTH_W2TH_O_DDTDZ(IIJB:IIJE,IKE)

!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_WTH_W2TH_O_DDTDZ',1,ZHOOK_HANDLE)
END SUBROUTINE D_M3_WTH_W2TH_O_DDTDZ
!----------------------------------------------------------------------------
SUBROUTINE M3_WTH_W2R(D,CSTURB,TURBN,PD,PKEFF,PTKE,PBLL_O_E,PEMOIST,PDTDZ,PM3_WTH_W2R)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  TYPE(TURB_t),                    INTENT(IN)    :: TURBN
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PKEFF
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PTKE
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PEMOIST
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PDTDZ
  REAL, DIMENSION(D%NIJT,D%NKT),INTENT(OUT) :: PM3_WTH_W2R
  REAL, DIMENSION(D%NIJT,D%NKT) :: ZWORK1 ! working array
  INTEGER :: IKB, IKE, JIJ,JK,IIJB,IIJE,IKT
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_WTH_W2R',0,ZHOOK_HANDLE)
IKB=D%NKTB
IKE=D%NKTE
IIJE=D%NIJE
IIJB=D%NIJB
IKT=D%NKT
!
CALL MZM_PHY(D,PTKE,ZWORK1)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    PM3_WTH_W2R(JIJ, JK) = &
      - TURBN%XCSHF*PKEFF(JIJ, JK)*0.75*TURBN%XCTV*PBLL_O_E(JIJ, JK) &
      /ZWORK1(JIJ, JK)*PEMOIST(JIJ, JK)*PDTDZ(JIJ, JK)/PD(JIJ, JK)
  END DO
END DO
!
PM3_WTH_W2R(IIJB:IIJE,IKB-1)=PM3_WTH_W2R(IIJB:IIJE,IKB)
PM3_WTH_W2R(IIJB:IIJE,IKE+1)=PM3_WTH_W2R(IIJB:IIJE,IKE)

!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_WTH_W2R',1,ZHOOK_HANDLE)
END SUBROUTINE M3_WTH_W2R
!----------------------------------------------------------------------------
SUBROUTINE D_M3_WTH_W2R_O_DDTDZ(D,CSTURB,TURBN,PREDTH1,PREDR1,PD,PKEFF,PTKE,PBLL_O_E,PEMOIST,PD_M3_WTH_W2R_O_DDTDZ)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  TYPE(TURB_t),                    INTENT(IN)    :: TURBN
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PKEFF
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PTKE
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PEMOIST
  REAL, DIMENSION(D%NIJT,D%NKT),INTENT(OUT) :: PD_M3_WTH_W2R_O_DDTDZ
  REAL, DIMENSION(D%NIJT,D%NKT) :: ZWORK1 ! working array
  INTEGER :: IKB, IKE, JIJ,JK,IIJB,IIJE,IKT
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_WTH_W2R_O_DDTDZ',0,ZHOOK_HANDLE)
IKB=D%NKTB
IKE=D%NKTE
IIJE=D%NIJE
IIJB=D%NIJB
IKT=D%NKT
!
CALL MZM_PHY(D,PTKE,ZWORK1)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    PD_M3_WTH_W2R_O_DDTDZ(JIJ, JK) = &
    - TURBN%XCSHF*PKEFF(JIJ, JK)*0.75*TURBN%XCTV*PBLL_O_E(JIJ, JK) &
                                   /ZWORK1(JIJ, JK)*PEMOIST(JIJ, JK)/PD(JIJ, JK) &
                                         * (1. -  PREDTH1(JIJ, JK)*(1.5+PREDTH1(JIJ, JK)& 
                                         +PREDR1(JIJ, JK))/PD(JIJ, JK))
  END DO
END DO
!
PD_M3_WTH_W2R_O_DDTDZ(IIJB:IIJE,IKB-1)=PD_M3_WTH_W2R_O_DDTDZ(IIJB:IIJE,IKB)
PD_M3_WTH_W2R_O_DDTDZ(IIJB:IIJE,IKE+1)=PD_M3_WTH_W2R_O_DDTDZ(IIJB:IIJE,IKE)

!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_WTH_W2R_O_DDTDZ',1,ZHOOK_HANDLE)
END SUBROUTINE D_M3_WTH_W2R_O_DDTDZ
!----------------------------------------------------------------------------
SUBROUTINE M3_WTH_WR2(D,CSTURB,TURBN,PD,PKEFF,PTKE,PSQRT_TKE,PBLL_O_E,PBETA,PLEPS,PEMOIST,PDTDZ,PM3_WTH_WR2)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  TYPE(TURB_t),                    INTENT(IN)    :: TURBN
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PKEFF
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PTKE
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PBETA
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PEMOIST
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PDTDZ
  REAL, DIMENSION(D%NIJT,D%NKT),INTENT(OUT) :: PM3_WTH_WR2
  REAL, DIMENSION(D%NIJT,D%NKT) :: ZWORK1,ZWORK2 ! working array
  INTEGER :: IKB, IKE, JIJ,JK,IIJB,IIJE,IKT
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_WTH_WR2',0,ZHOOK_HANDLE)
IKB=D%NKTB
IKE=D%NKTE
IIJE=D%NIJE
IIJB=D%NIJB
IKT=D%NKT
!

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZWORK1(JIJ, JK) = PBETA(JIJ, JK)*PLEPS(JIJ, JK) &
                                     /(PSQRT_TKE(JIJ, JK)*PTKE(JIJ, JK))
  END DO
END DO

CALL MZM_PHY(D,ZWORK1,ZWORK2)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    PM3_WTH_WR2(JIJ, JK) = - TURBN%XCSHF*PKEFF(JIJ, JK)& 
                               *0.25*PBLL_O_E(JIJ, JK)*TURBN%XCTV*PEMOIST(JIJ, JK)**2 &
                               *ZWORK2(JIJ, JK)/CSTURB%XCTD*PDTDZ(JIJ, JK)/PD(JIJ, JK)
  END DO
END DO
!
PM3_WTH_WR2(IIJB:IIJE,IKB-1)=PM3_WTH_WR2(IIJB:IIJE,IKB)
PM3_WTH_WR2(IIJB:IIJE,IKE+1)=PM3_WTH_WR2(IIJB:IIJE,IKE)

!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_WTH_WR2',1,ZHOOK_HANDLE)
END SUBROUTINE M3_WTH_WR2
!----------------------------------------------------------------------------
SUBROUTINE D_M3_WTH_WR2_O_DDTDZ(D,CSTURB,TURBN,PREDTH1,PREDR1,PD,PKEFF,PTKE,PSQRT_TKE,PBLL_O_E,&
                               &PBETA,PLEPS,PEMOIST,PD_M3_WTH_WR2_O_DDTDZ)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  TYPE(TURB_t),                    INTENT(IN)    :: TURBN
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PKEFF
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PTKE
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PBETA
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PEMOIST
  REAL, DIMENSION(D%NIJT,D%NKT),INTENT(OUT) :: PD_M3_WTH_WR2_O_DDTDZ
  REAL, DIMENSION(D%NIJT,D%NKT) :: ZWORK1,ZWORK2 ! working array
  INTEGER :: IKB, IKE, JIJ,JK,IIJB,IIJE,IKT
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_WTH_WR2_O_DDTDZ',0,ZHOOK_HANDLE)
IKB=D%NKTB
IKE=D%NKTE
IIJE=D%NIJE
IIJB=D%NIJB
IKT=D%NKT
!

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZWORK1(JIJ, JK) = PBETA(JIJ, JK)*PLEPS(JIJ, JK)&
                                      /(PSQRT_TKE(JIJ, JK)*PTKE(JIJ, JK))
  END DO
END DO

CALL MZM_PHY(D,ZWORK1,ZWORK2)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    PD_M3_WTH_WR2_O_DDTDZ(JIJ, JK) = - TURBN%XCSHF*PKEFF(JIJ, JK)& 
                               *0.25*PBLL_O_E(JIJ, JK)*TURBN%XCTV*PEMOIST(JIJ, JK)**2 &
                               *ZWORK2(JIJ, JK)/CSTURB%XCTD/PD(JIJ, JK)     &
                               * (1. -  PREDTH1(JIJ, JK)* &
                               (1.5+PREDTH1(JIJ, JK)+PREDR1(JIJ, JK))/PD(JIJ, JK))
  END DO
END DO
!
PD_M3_WTH_WR2_O_DDTDZ(IIJB:IIJE,IKB-1)=PD_M3_WTH_WR2_O_DDTDZ(IIJB:IIJE,IKB)
PD_M3_WTH_WR2_O_DDTDZ(IIJB:IIJE,IKE+1)=PD_M3_WTH_WR2_O_DDTDZ(IIJB:IIJE,IKE)

!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_WTH_WR2_O_DDTDZ',1,ZHOOK_HANDLE)
END SUBROUTINE D_M3_WTH_WR2_O_DDTDZ
!----------------------------------------------------------------------------
SUBROUTINE M3_WTH_WTHR(D,CSTURB,TURBN,PREDR1,PD,PKEFF,PTKE,PSQRT_TKE,PBETA,PLEPS,PEMOIST,PM3_WTH_WTHR)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  TYPE(TURB_t),                    INTENT(IN)    :: TURBN
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PKEFF
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PTKE
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PBETA
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PEMOIST
  REAL, DIMENSION(D%NIJT,D%NKT),INTENT(OUT) :: PM3_WTH_WTHR
  REAL, DIMENSION(D%NIJT,D%NKT) :: ZWORK1,ZWORK2 ! working array
  INTEGER :: IKB, IKE, JIJ,JK,IIJB,IIJE,IKT
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_WTH_WTHR',0,ZHOOK_HANDLE)
IKB=D%NKTB
IKE=D%NKTE
IIJE=D%NIJE
IIJB=D%NIJB
IKT=D%NKT
!

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZWORK1(JIJ, JK) = PBETA(JIJ, JK)*PLEPS(JIJ, JK)&
                                      /(PSQRT_TKE(JIJ, JK)*PTKE(JIJ, JK))
  END DO
END DO

CALL MZM_PHY(D,ZWORK1,ZWORK2)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    PM3_WTH_WTHR(JIJ, JK) = &
                       TURBN%XCSHF*PKEFF(JIJ, JK)*PEMOIST(JIJ, JK)*ZWORK2(JIJ, JK) &
                       *0.5*PLEPS(JIJ, JK)/CSTURB%XCTD*(1+PREDR1(JIJ, JK))/PD(JIJ, JK)
  END DO
END DO
!
PM3_WTH_WTHR(IIJB:IIJE,IKB-1)=PM3_WTH_WTHR(IIJB:IIJE,IKB)
PM3_WTH_WTHR(IIJB:IIJE,IKE+1)=PM3_WTH_WTHR(IIJB:IIJE,IKE)

!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_WTH_WTHR',1,ZHOOK_HANDLE)
END SUBROUTINE M3_WTH_WTHR
!----------------------------------------------------------------------------
SUBROUTINE D_M3_WTH_WTHR_O_DDTDZ(D,CSTURB,TURBN,PM3_WTH_WTHR,PREDTH1,PREDR1,PD,PBLL_O_E,PETHETA,PD_M3_WTH_WTHR_O_DDTDZ)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  TYPE(TURB_t),                    INTENT(IN)    :: TURBN
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PM3_WTH_WTHR
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PETHETA
  REAL, DIMENSION(D%NIJT,D%NKT),INTENT(OUT) :: PD_M3_WTH_WTHR_O_DDTDZ
  INTEGER :: IKB, IKE, JIJ,JK,IIJB,IIJE,IKT
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_WTH_WTHR_O_DDTDZ',0,ZHOOK_HANDLE)
IKB=D%NKTB
IKE=D%NKTE
IIJE=D%NIJE
IIJB=D%NIJB
IKT=D%NKT
!

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    PD_M3_WTH_WTHR_O_DDTDZ(JIJ, JK) = &
                    - PM3_WTH_WTHR(JIJ, JK) * (1.5+PREDTH1(JIJ, JK)+PREDR1(JIJ, JK))&
                    /PD(JIJ, JK)*TURBN%XCTV*PBLL_O_E(JIJ, JK)*PETHETA(JIJ, JK)
  END DO
END DO
!
PD_M3_WTH_WTHR_O_DDTDZ(IIJB:IIJE,IKB-1)=PD_M3_WTH_WTHR_O_DDTDZ(IIJB:IIJE,IKB)
PD_M3_WTH_WTHR_O_DDTDZ(IIJB:IIJE,IKE+1)=PD_M3_WTH_WTHR_O_DDTDZ(IIJB:IIJE,IKE)

!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_WTH_WTHR_O_DDTDZ',1,ZHOOK_HANDLE)
END SUBROUTINE D_M3_WTH_WTHR_O_DDTDZ
!----------------------------------------------------------------------------
SUBROUTINE M3_TH2_W2TH(D,CSTURB,TURBN,PREDTH1,PREDR1,PD,PDTDZ,PLM,PLEPS,PTKE,PM3_TH2_W2TH)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  TYPE(TURB_t),                    INTENT(IN)    :: TURBN
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PDTDZ
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PLM
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PTKE
  REAL, DIMENSION(D%NIJT,D%NKT),INTENT(OUT) :: PM3_TH2_W2TH
  REAL, DIMENSION(D%NIJT,D%NKT) :: ZWORK1,ZWORK2 ! working array
  INTEGER :: IKB, IKE, JIJ,JK,IIJB,IIJE,IKT
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_TH2_W2TH',0,ZHOOK_HANDLE)
IKB=D%NKTB
IKE=D%NKTE
IIJE=D%NIJE
IIJB=D%NIJB
IKT=D%NKT
!

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZWORK1(JIJ, JK) = (1.-0.5*PREDR1(JIJ, JK)*(1.+PREDR1(JIJ, JK))& 
                                    /PD(JIJ, JK))/(1.+PREDTH1(JIJ, JK))*PDTDZ(JIJ, JK)
  END DO
END DO

CALL MZF_PHY(D,ZWORK1,ZWORK2)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    PM3_TH2_W2TH(JIJ, JK) = - ZWORK2(JIJ, JK) &
                           * 1.5*PLM(JIJ, JK)*PLEPS(JIJ, JK)/PTKE(JIJ, JK)*TURBN%XCTV
  END DO
END DO
!
PM3_TH2_W2TH(IIJB:IIJE,IKB-1)=PM3_TH2_W2TH(IIJB:IIJE,IKB)
PM3_TH2_W2TH(IIJB:IIJE,IKE+1)=PM3_TH2_W2TH(IIJB:IIJE,IKE)

!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_TH2_W2TH',1,ZHOOK_HANDLE)
END SUBROUTINE M3_TH2_W2TH
!----------------------------------------------------------------------------
SUBROUTINE D_M3_TH2_W2TH_O_DDTDZ(D,CSTURB,TURBN,PREDTH1,PREDR1,PD,PLM,PLEPS,PTKE,OUSERV,PD_M3_TH2_W2TH_O_DDTDZ)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  TYPE(TURB_t),                    INTENT(IN)    :: TURBN
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PLM
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PTKE
  LOGICAL,                INTENT(IN) :: OUSERV
  REAL, DIMENSION(D%NIJT,D%NKT),INTENT(OUT) :: PD_M3_TH2_W2TH_O_DDTDZ
  REAL, DIMENSION(D%NIJT,D%NKT) :: ZWORK1,ZWORK2 ! working array
  INTEGER :: IKB, IKE, JIJ,JK,IIJB,IIJE,IKT
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_TH2_W2TH_O_DDTDZ',0,ZHOOK_HANDLE)
IKB=D%NKTB
IKE=D%NKTE
IIJE=D%NIJE
IIJB=D%NIJB
IKT=D%NKT
!
IF (OUSERV) THEN
!  D_M3_TH2_W2TH_O_DDTDZ(IIJB:IIJE,1:IKT) = - 1.5*PLM*PLEPS/PTKE*TURBN%XCTV * MZF(                    &
!          (1.-0.5*PREDR1*(1.+PREDR1)/PD)*(1.-(1.5+PREDTH1+PREDR1)*(1.+PREDTH1)/PD )  &
!        / (1.+PREDTH1)**2, IKA, IKU, IKL)

  DO JK=1, IKT
    DO JIJ=IIJB, IIJE
      ZWORK1(JIJ, JK) = (1.-0.5*PREDR1(JIJ, JK)*(1.+PREDR1(JIJ, JK))&
                 / PD(JIJ, JK))*(1.-(1.5+PREDTH1(JIJ, JK)+PREDR1(JIJ, JK))   &
                 * PREDTH1(JIJ, JK)*(1.+PREDTH1(JIJ, JK))/PD(JIJ, JK) ) &
                 / (1.+PREDTH1(JIJ, JK))**2
    END DO
  END DO

  CALL MZF_PHY(D,ZWORK1,ZWORK2)

  DO JK=1, IKT
    DO JIJ=IIJB, IIJE
      PD_M3_TH2_W2TH_O_DDTDZ(JIJ, JK) = - 1.5*PLM(JIJ, JK)*PLEPS(JIJ, JK) &
                                                       /PTKE(JIJ, JK)*TURBN%XCTV * ZWORK2(JIJ, JK)
    END DO
  END DO

ELSE

  DO JK=1, IKT
    DO JIJ=IIJB, IIJE
      ZWORK1(JIJ, JK) = 1./(1.+PREDTH1(JIJ, JK))**2
    END DO
  END DO

  CALL MZF_PHY(D,ZWORK1,ZWORK2)

  DO JK=1, IKT
    DO JIJ=IIJB, IIJE
      PD_M3_TH2_W2TH_O_DDTDZ(JIJ, JK) = - 1.5*PLM(JIJ, JK)*PLEPS(JIJ, JK) & 
                                                       /PTKE(JIJ, JK)*TURBN%XCTV * ZWORK2(JIJ, JK)
    END DO
  END DO

END IF
!

PD_M3_TH2_W2TH_O_DDTDZ(IIJB:IIJE,IKB-1)=PD_M3_TH2_W2TH_O_DDTDZ(IIJB:IIJE,IKB)
PD_M3_TH2_W2TH_O_DDTDZ(IIJB:IIJE,IKE+1)=PD_M3_TH2_W2TH_O_DDTDZ(IIJB:IIJE,IKE)

!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_TH2_W2TH_O_DDTDZ',1,ZHOOK_HANDLE)
END SUBROUTINE D_M3_TH2_W2TH_O_DDTDZ
!----------------------------------------------------------------------------
SUBROUTINE M3_TH2_WTH2(D,CSTURB,TURBN,PREDTH1,PREDR1,PD,PLEPS,PSQRT_TKE,PM3_TH2_WTH2)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  TYPE(TURB_t),                    INTENT(IN)    :: TURBN
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(D%NIJT,D%NKT),INTENT(OUT) :: PM3_TH2_WTH2
  REAL, DIMENSION(D%NIJT,D%NKT) :: ZWORK1,ZWORK2 ! working array
  INTEGER :: IKB, IKE, JIJ,JK,IIJB,IIJE,IKT
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_TH2_WTH2',0,ZHOOK_HANDLE)
IKB=D%NKTB
IKE=D%NKTE
IIJE=D%NIJE
IIJB=D%NIJB
IKT=D%NKT
!

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZWORK1(JIJ, JK) = (1.+0.5*PREDTH1(JIJ, JK) &
                             +1.5*PREDR1(JIJ, JK)+0.5*PREDR1(JIJ, JK)**2)/PD(JIJ, JK)
  END DO
END DO

CALL MZF_PHY(D,ZWORK1,ZWORK2)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    PM3_TH2_WTH2(JIJ, JK) = PLEPS(JIJ, JK)*0.5/CSTURB%XCTD/PSQRT_TKE(JIJ, JK) &
                         * ZWORK2(JIJ, JK)
  END DO
END DO
!
PM3_TH2_WTH2(IIJB:IIJE,IKB-1)=PM3_TH2_WTH2(IIJB:IIJE,IKB)
PM3_TH2_WTH2(IIJB:IIJE,IKE+1)=PM3_TH2_WTH2(IIJB:IIJE,IKE)

!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_TH2_WTH2',1,ZHOOK_HANDLE)
END SUBROUTINE M3_TH2_WTH2
!----------------------------------------------------------------------------
SUBROUTINE D_M3_TH2_WTH2_O_DDTDZ(D,CSTURB,TURBN,PREDTH1,PREDR1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PETHETA,PD_M3_TH2_WTH2_O_DDTDZ)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  TYPE(TURB_t),                    INTENT(IN)    :: TURBN
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PETHETA
  REAL, DIMENSION(D%NIJT,D%NKT),INTENT(OUT) :: PD_M3_TH2_WTH2_O_DDTDZ
  REAL, DIMENSION(D%NIJT,D%NKT) :: ZWORK1,ZWORK2 ! working array
  INTEGER :: IKB, IKE, JIJ,JK,IIJB,IIJE,IKT
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_TH2_WTH2_O_DDTDZ',0,ZHOOK_HANDLE)
IKB=D%NKTB
IKE=D%NKTE
IIJE=D%NIJE
IIJB=D%NIJB
IKT=D%NKT
!

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZWORK1(JIJ, JK) = PBLL_O_E(JIJ, JK)*PETHETA(JIJ, JK) &
                 * (0.5/PD(JIJ, JK) - (1.5+PREDTH1(JIJ, JK)+PREDR1(JIJ, JK))& 
                 *(1.+0.5*PREDTH1(JIJ, JK)+1.5*PREDR1(JIJ, JK)& 
                 +0.5*PREDR1(JIJ, JK)**2)/PD(JIJ, JK)**2)
   END DO
END DO

CALL MZF_PHY(D,ZWORK1,ZWORK2)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    PD_M3_TH2_WTH2_O_DDTDZ(JIJ, JK) = PLEPS(JIJ, JK) & 
                                     *0.5/CSTURB%XCTD/PSQRT_TKE(JIJ, JK)*TURBN%XCTV * ZWORK2(JIJ, JK)
  END DO
END DO
!
PD_M3_TH2_WTH2_O_DDTDZ(IIJB:IIJE,IKB-1)=PD_M3_TH2_WTH2_O_DDTDZ(IIJB:IIJE,IKB)
PD_M3_TH2_WTH2_O_DDTDZ(IIJB:IIJE,IKE+1)=PD_M3_TH2_WTH2_O_DDTDZ(IIJB:IIJE,IKE)

!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_TH2_WTH2_O_DDTDZ',1,ZHOOK_HANDLE)
END SUBROUTINE D_M3_TH2_WTH2_O_DDTDZ
!----------------------------------------------------------------------------
SUBROUTINE M3_TH2_W2R(D,CSTURB,TURBN,PD,PLM,PLEPS,PTKE,PBLL_O_E,PEMOIST,PDTDZ,PM3_TH2_W2R)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  TYPE(TURB_t),                    INTENT(IN)    :: TURBN
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PLM
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PTKE
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PEMOIST
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PDTDZ
  REAL, DIMENSION(D%NIJT,D%NKT),INTENT(OUT) :: PM3_TH2_W2R
  REAL, DIMENSION(D%NIJT,D%NKT) :: ZWORK1,ZWORK2 ! working array
  INTEGER :: IKB, IKE, JIJ,JK,IIJB,IIJE,IKT
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_TH2_W2R',0,ZHOOK_HANDLE)
IKB=D%NKTB
IKE=D%NKTE
IIJE=D%NIJE
IIJB=D%NIJB
IKT=D%NKT
!

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZWORK1(JIJ, JK) = PBLL_O_E(JIJ, JK)*PEMOIST(JIJ, JK) & 
                                     /PD(JIJ, JK)*PDTDZ(JIJ, JK)**2
  END DO
END DO

CALL MZF_PHY(D,ZWORK1,ZWORK2)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    PM3_TH2_W2R(JIJ, JK) = 0.75*TURBN%XCTV**2*ZWORK2(JIJ, JK) &
                        *PLM(JIJ, JK)*PLEPS(JIJ, JK)/PTKE(JIJ, JK)
  END DO
END DO
!
PM3_TH2_W2R(IIJB:IIJE,IKB-1)=PM3_TH2_W2R(IIJB:IIJE,IKB)
PM3_TH2_W2R(IIJB:IIJE,IKE+1)=PM3_TH2_W2R(IIJB:IIJE,IKE)

!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_TH2_W2R',1,ZHOOK_HANDLE)
END SUBROUTINE M3_TH2_W2R
!----------------------------------------------------------------------------
SUBROUTINE D_M3_TH2_W2R_O_DDTDZ(D,CSTURB,TURBN,PREDTH1,PREDR1,PD,PLM,PLEPS,PTKE,PBLL_O_E,PEMOIST,PDTDZ,PD_M3_TH2_W2R_O_DDTDZ)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  TYPE(TURB_t),                    INTENT(IN)    :: TURBN
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PLM
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PTKE
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PEMOIST
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PDTDZ
  REAL, DIMENSION(D%NIJT,D%NKT),INTENT(OUT) :: PD_M3_TH2_W2R_O_DDTDZ
  REAL, DIMENSION(D%NIJT,D%NKT) :: ZWORK1,ZWORK2 ! working array
  INTEGER :: IKB, IKE, JIJ,JK,IIJB,IIJE,IKT
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_TH2_W2R_O_DDTDZ',0,ZHOOK_HANDLE)
IKB=D%NKTB
IKE=D%NKTE
IIJE=D%NIJE
IIJB=D%NIJB
IKT=D%NKT
!

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZWORK1(JIJ, JK) =  PBLL_O_E(JIJ, JK)*PEMOIST(JIJ, JK)& 
     /PD(JIJ, JK)*PDTDZ(JIJ, JK)*(2.-PREDTH1(JIJ, JK)* & 
     (1.5+PREDTH1(JIJ, JK)+PREDR1(JIJ, JK))/PD(JIJ, JK))
  END DO
END DO

CALL MZF_PHY(D,ZWORK1,ZWORK2)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    PD_M3_TH2_W2R_O_DDTDZ(JIJ, JK) = 0.75*TURBN%XCTV**2*PLM(JIJ, JK) *PLEPS(JIJ, JK) &
                                                    /PTKE(JIJ, JK) * ZWORK2(JIJ, JK)
  END DO
END DO
!
PD_M3_TH2_W2R_O_DDTDZ(IIJB:IIJE,IKB-1)=PD_M3_TH2_W2R_O_DDTDZ(IIJB:IIJE,IKB)
PD_M3_TH2_W2R_O_DDTDZ(IIJB:IIJE,IKE+1)=PD_M3_TH2_W2R_O_DDTDZ(IIJB:IIJE,IKE)

!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_TH2_W2R_O_DDTDZ',1,ZHOOK_HANDLE)
END SUBROUTINE D_M3_TH2_W2R_O_DDTDZ
!----------------------------------------------------------------------------
SUBROUTINE M3_TH2_WR2(D,CSTURB,TURBN,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PEMOIST,PDTDZ,PM3_TH2_WR2)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  TYPE(TURB_t),                    INTENT(IN)    :: TURBN
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PEMOIST
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PDTDZ
  REAL, DIMENSION(D%NIJT,D%NKT),INTENT(OUT) :: PM3_TH2_WR2
  REAL, DIMENSION(D%NIJT,D%NKT) :: ZWORK1,ZWORK2 ! working array
  INTEGER :: IKB, IKE, JIJ,JK,IIJB,IIJE,IKT
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_TH2_WR2',0,ZHOOK_HANDLE)
IKB=D%NKTB
IKE=D%NKTE
IIJE=D%NIJE
IIJB=D%NIJB
IKT=D%NKT
!

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZWORK1(JIJ, JK) = (PBLL_O_E(JIJ, JK)*PEMOIST(JIJ, JK)& 
                                      *PDTDZ(JIJ, JK))**2/PD(JIJ, JK)
  END DO
END DO

CALL MZF_PHY(D,ZWORK1,ZWORK2)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    PM3_TH2_WR2(JIJ, JK) = 0.25*TURBN%XCTV**2*ZWORK2(JIJ, JK)&
                        *PLEPS(JIJ, JK)/PSQRT_TKE(JIJ, JK)/CSTURB%XCTD
  END DO
END DO
!
PM3_TH2_WR2(IIJB:IIJE,IKB-1)=PM3_TH2_WR2(IIJB:IIJE,IKB)
PM3_TH2_WR2(IIJB:IIJE,IKE+1)=PM3_TH2_WR2(IIJB:IIJE,IKE)

!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_TH2_WR2',1,ZHOOK_HANDLE)
END SUBROUTINE M3_TH2_WR2
!----------------------------------------------------------------------------
SUBROUTINE D_M3_TH2_WR2_O_DDTDZ(D,CSTURB,TURBN,PREDTH1,PREDR1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PEMOIST,PDTDZ,PD_M3_TH2_WR2_O_DDTDZ)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  TYPE(TURB_t),                    INTENT(IN)    :: TURBN
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PEMOIST
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PDTDZ
  REAL, DIMENSION(D%NIJT,D%NKT),INTENT(OUT) :: PD_M3_TH2_WR2_O_DDTDZ
  REAL, DIMENSION(D%NIJT,D%NKT) :: ZWORK1,ZWORK2 ! working array
  INTEGER :: IKB, IKE, JIJ,JK,IIJB,IIJE,IKT
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_TH2_WR2_O_DDTDZ',0,ZHOOK_HANDLE)
IKB=D%NKTB
IKE=D%NKTE
IIJE=D%NIJE
IIJB=D%NIJB
IKT=D%NKT
!

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZWORK1(JIJ, JK) = (PBLL_O_E(JIJ, JK)*PEMOIST(JIJ, JK))**2 & 
    *PDTDZ(JIJ, JK)/PD(JIJ, JK)*(2.-PREDTH1(JIJ, JK) & 
    *(1.5+PREDTH1(JIJ, JK)+PREDR1(JIJ, JK))/PD(JIJ, JK))
  END DO
END DO

CALL MZF_PHY(D,ZWORK1,ZWORK2)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    PD_M3_TH2_WR2_O_DDTDZ(JIJ, JK) = 0.25*TURBN%XCTV**2*PLEPS(JIJ, JK) & 
                                                   / PSQRT_TKE(JIJ, JK)/CSTURB%XCTD * ZWORK2(JIJ, JK)
  END DO
END DO
!
PD_M3_TH2_WR2_O_DDTDZ(IIJB:IIJE,IKB-1)=PD_M3_TH2_WR2_O_DDTDZ(IIJB:IIJE,IKB)
PD_M3_TH2_WR2_O_DDTDZ(IIJB:IIJE,IKE+1)=PD_M3_TH2_WR2_O_DDTDZ(IIJB:IIJE,IKE)

!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_TH2_WR2_O_DDTDZ',1,ZHOOK_HANDLE)
END SUBROUTINE D_M3_TH2_WR2_O_DDTDZ
!----------------------------------------------------------------------------
SUBROUTINE M3_TH2_WTHR(D,CSTURB,TURBN,PREDR1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PEMOIST,PDTDZ,PM3_TH2_WTHR)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  TYPE(TURB_t),                    INTENT(IN)    :: TURBN
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PEMOIST
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PDTDZ
  REAL, DIMENSION(D%NIJT,D%NKT),INTENT(OUT) :: PM3_TH2_WTHR
  REAL, DIMENSION(D%NIJT,D%NKT) :: ZWORK1,ZWORK2 ! working array
  INTEGER :: IKB, IKE, JIJ,JK,IIJB,IIJE,IKT
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_TH2_WTHR',0,ZHOOK_HANDLE)
IKB=D%NKTB
IKE=D%NKTE
IIJE=D%NIJE
IIJB=D%NIJB
IKT=D%NKT
!

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZWORK1(JIJ, JK) = PBLL_O_E(JIJ, JK)*PEMOIST(JIJ, JK) & 
                                    * PDTDZ(JIJ, JK)*(1.+PREDR1(JIJ, JK))/PD(JIJ, JK)
  END DO
END DO

CALL MZF_PHY(D,ZWORK1,ZWORK2)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    PM3_TH2_WTHR(JIJ, JK) = - 0.5*TURBN%XCTV*PLEPS(JIJ, JK) & 
                                           / PSQRT_TKE(JIJ, JK)/CSTURB%XCTD * ZWORK2(JIJ, JK)
  END DO
END DO
!
PM3_TH2_WTHR(IIJB:IIJE,IKB-1)=PM3_TH2_WTHR(IIJB:IIJE,IKB)
PM3_TH2_WTHR(IIJB:IIJE,IKE+1)=PM3_TH2_WTHR(IIJB:IIJE,IKE)

!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_TH2_WTHR',1,ZHOOK_HANDLE)
END SUBROUTINE M3_TH2_WTHR
!----------------------------------------------------------------------------
SUBROUTINE D_M3_TH2_WTHR_O_DDTDZ(D,CSTURB,TURBN,PREDTH1,PREDR1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PEMOIST,PDTDZ,PD_M3_TH2_WTHR_O_DDTDZ)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  TYPE(TURB_t),                    INTENT(IN)    :: TURBN
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PEMOIST
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PDTDZ
  REAL, DIMENSION(D%NIJT,D%NKT),INTENT(OUT) :: PD_M3_TH2_WTHR_O_DDTDZ
  REAL, DIMENSION(D%NIJT,D%NKT) :: ZWORK1,ZWORK2 ! working array
  INTEGER :: IKB, IKE, JIJ,JK,IIJB,IIJE,IKT
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_TH2_WTHR_O_DDTDZ',0,ZHOOK_HANDLE)
IKB=D%NKTB
IKE=D%NKTE
IIJE=D%NIJE
IIJB=D%NIJB
IKT=D%NKT
!

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZWORK1(JIJ, JK) = PBLL_O_E(JIJ, JK)*PEMOIST(JIJ, JK)* & 
                     (1.+PREDR1(JIJ, JK))/PD(JIJ, JK) * (1. -PREDTH1(JIJ, JK)*& 
                     (1.5+PREDTH1(JIJ, JK)+PREDR1(JIJ, JK))/PD(JIJ, JK))
  END DO
END DO

CALL MZF_PHY(D,ZWORK1,ZWORK2)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    PD_M3_TH2_WTHR_O_DDTDZ(JIJ, JK) = - 0.5*TURBN%XCTV*PLEPS(JIJ, JK) & 
                                                    / PSQRT_TKE(JIJ, JK)/CSTURB%XCTD * ZWORK2(JIJ, JK) 
  END DO
END DO
!
PD_M3_TH2_WTHR_O_DDTDZ(IIJB:IIJE,IKB-1)=PD_M3_TH2_WTHR_O_DDTDZ(IIJB:IIJE,IKB)
PD_M3_TH2_WTHR_O_DDTDZ(IIJB:IIJE,IKE+1)=PD_M3_TH2_WTHR_O_DDTDZ(IIJB:IIJE,IKE)

!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_TH2_WTHR_O_DDTDZ',1,ZHOOK_HANDLE)
END SUBROUTINE D_M3_TH2_WTHR_O_DDTDZ
!----------------------------------------------------------------------------
SUBROUTINE M3_THR_WTHR(D,CSTURB,TURBN,PREDTH1,PREDR1,PD,PLEPS,PSQRT_TKE,PM3_THR_WTHR)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  TYPE(TURB_t),                    INTENT(IN)    :: TURBN
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(D%NIJT,D%NKT),INTENT(OUT) :: PM3_THR_WTHR
  REAL, DIMENSION(D%NIJT,D%NKT) :: ZWORK1,ZWORK2 ! working array
  INTEGER :: IKB, IKE, JIJ,JK,IIJB,IIJE,IKT
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_THR_WTHR',0,ZHOOK_HANDLE)
IKB=D%NKTB
IKE=D%NKTE
IIJE=D%NIJE
IIJB=D%NIJB
IKT=D%NKT
!

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZWORK1(JIJ, JK) =  (1.+PREDTH1(JIJ, JK))* & 
                                       (1.+PREDR1(JIJ, JK))/PD(JIJ, JK)
  END DO
END DO

CALL MZF_PHY(D,ZWORK1,ZWORK2)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    PM3_THR_WTHR(JIJ, JK) = 0.5*PLEPS(JIJ, JK)/PSQRT_TKE(JIJ, JK)/CSTURB%XCTD &
                         * ZWORK2(JIJ, JK)
  END DO
END DO
!
PM3_THR_WTHR(IIJB:IIJE,IKB-1)=PM3_THR_WTHR(IIJB:IIJE,IKB)
PM3_THR_WTHR(IIJB:IIJE,IKE+1)=PM3_THR_WTHR(IIJB:IIJE,IKE)

!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_THR_WTHR',1,ZHOOK_HANDLE)
END SUBROUTINE M3_THR_WTHR
!----------------------------------------------------------------------------
SUBROUTINE D_M3_THR_WTHR_O_DDTDZ(D,CSTURB,TURBN,PREDTH1,PREDR1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PETHETA,PD_M3_THR_WTHR_O_DDTDZ)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  TYPE(TURB_t),                    INTENT(IN)    :: TURBN
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PETHETA
  REAL, DIMENSION(D%NIJT,D%NKT),INTENT(OUT) :: PD_M3_THR_WTHR_O_DDTDZ
  REAL, DIMENSION(D%NIJT,D%NKT) :: ZWORK1,ZWORK2 ! working array
  INTEGER :: IKB, IKE, JIJ,JK,IIJB,IIJE,IKT
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_THR_WTHR_O_DDTDZ',0,ZHOOK_HANDLE)
IKB=D%NKTB
IKE=D%NKTE
IIJE=D%NIJE
IIJB=D%NIJB
IKT=D%NKT
!

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZWORK1(JIJ, JK) = PETHETA(JIJ, JK)*PBLL_O_E(JIJ, JK)/PD(JIJ, JK) & 
                                 *(1.+PREDR1(JIJ, JK))*(1.-(1.+PREDTH1(JIJ, JK)) & 
                                 *(1.5+PREDTH1(JIJ, JK)+PREDR1(JIJ, JK))/PD(JIJ, JK))
  END DO
END DO

CALL MZF_PHY(D,ZWORK1,ZWORK2)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    PD_M3_THR_WTHR_O_DDTDZ(JIJ, JK) = 0.5*PLEPS(JIJ, JK)/PSQRT_TKE(JIJ, JK) & 
                                                    / CSTURB%XCTD * TURBN%XCTV * ZWORK2(JIJ, JK)
  END DO
END DO
!
PD_M3_THR_WTHR_O_DDTDZ(IIJB:IIJE,IKB-1)=PD_M3_THR_WTHR_O_DDTDZ(IIJB:IIJE,IKB)
PD_M3_THR_WTHR_O_DDTDZ(IIJB:IIJE,IKE+1)=PD_M3_THR_WTHR_O_DDTDZ(IIJB:IIJE,IKE)

!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_THR_WTHR_O_DDTDZ',1,ZHOOK_HANDLE)
END SUBROUTINE D_M3_THR_WTHR_O_DDTDZ
!----------------------------------------------------------------------------
SUBROUTINE M3_THR_WTH2(D,CSTURB,TURBN,PREDR1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PETHETA,PDRDZ,PM3_THR_WTH2)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  TYPE(TURB_t),                    INTENT(IN)    :: TURBN
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PETHETA
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PDRDZ
  REAL, DIMENSION(D%NIJT,D%NKT),INTENT(OUT) :: PM3_THR_WTH2
  REAL, DIMENSION(D%NIJT,D%NKT) :: ZWORK1,ZWORK2 ! working array
  INTEGER :: IKB, IKE, JIJ,JK,IIJB,IIJE,IKT
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_THR_WTH2',0,ZHOOK_HANDLE)
IKB=D%NKTB
IKE=D%NKTE
IIJE=D%NIJE
IIJB=D%NIJB
IKT=D%NKT
!

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZWORK1(JIJ, JK) = (1.+PREDR1(JIJ, JK))*PBLL_O_E(JIJ, JK)* & 
                                      PETHETA(JIJ, JK)*PDRDZ(JIJ, JK)/PD(JIJ, JK)
  END DO
END DO

CALL MZF_PHY(D,ZWORK1,ZWORK2)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    PM3_THR_WTH2(JIJ, JK) = - 0.25*PLEPS(JIJ, JK) & 
                                        / PSQRT_TKE(JIJ, JK)/CSTURB%XCTD*TURBN%XCTV * ZWORK2(JIJ, JK)
  END DO
END DO
!
PM3_THR_WTH2(IIJB:IIJE,IKB-1)=PM3_THR_WTH2(IIJB:IIJE,IKB)
PM3_THR_WTH2(IIJB:IIJE,IKE+1)=PM3_THR_WTH2(IIJB:IIJE,IKE)

!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_THR_WTH2',1,ZHOOK_HANDLE)
END SUBROUTINE M3_THR_WTH2
!----------------------------------------------------------------------------
SUBROUTINE D_M3_THR_WTH2_O_DDTDZ(D,CSTURB,TURBN,PREDTH1,PREDR1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PETHETA,PDRDZ,PD_M3_THR_WTH2_O_DDTDZ)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  TYPE(TURB_t),                    INTENT(IN)    :: TURBN
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PETHETA
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PDRDZ
  REAL, DIMENSION(D%NIJT,D%NKT),INTENT(OUT) :: PD_M3_THR_WTH2_O_DDTDZ
  REAL, DIMENSION(D%NIJT,D%NKT) :: ZWORK1,ZWORK2 ! working array
  INTEGER :: IKB, IKE, JIJ,JK,IIJB,IIJE,IKT
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_THR_WTH2_O_DDTDZ',0,ZHOOK_HANDLE)
IKB=D%NKTB
IKE=D%NKTE
IIJE=D%NIJE
IIJB=D%NIJB
IKT=D%NKT
!

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZWORK1(JIJ, JK) = -(1.+PREDR1(JIJ, JK))*(PBLL_O_E(JIJ, JK) & 
                                     *PETHETA(JIJ, JK)/PD(JIJ, JK))**2&
                                     *PDRDZ(JIJ, JK)&
                                     *(1.5+PREDTH1(JIJ, JK)+PREDR1(JIJ, JK))
  END DO
END DO

CALL MZF_PHY(D,ZWORK1,ZWORK2)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    PD_M3_THR_WTH2_O_DDTDZ(JIJ, JK) = - 0.25*PLEPS(JIJ, JK) &
                                     /PSQRT_TKE(JIJ, JK)/CSTURB%XCTD*TURBN%XCTV**2 * ZWORK2(JIJ, JK) 
  END DO
END DO
!
PD_M3_THR_WTH2_O_DDTDZ(IIJB:IIJE,IKB-1)=PD_M3_THR_WTH2_O_DDTDZ(IIJB:IIJE,IKB)
PD_M3_THR_WTH2_O_DDTDZ(IIJB:IIJE,IKE+1)=PD_M3_THR_WTH2_O_DDTDZ(IIJB:IIJE,IKE)

!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_THR_WTH2_O_DDTDZ',1,ZHOOK_HANDLE)
END SUBROUTINE D_M3_THR_WTH2_O_DDTDZ
!----------------------------------------------------------------------------
SUBROUTINE D_M3_THR_WTH2_O_DDRDZ(D,CSTURB,TURBN,PREDTH1,PREDR1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PETHETA,PD_M3_THR_WTH2_O_DDRDZ)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  TYPE(TURB_t),                    INTENT(IN)    :: TURBN
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PETHETA
  REAL, DIMENSION(D%NIJT,D%NKT),INTENT(OUT) :: PD_M3_THR_WTH2_O_DDRDZ
  REAL, DIMENSION(D%NIJT,D%NKT) :: ZWORK1,ZWORK2 ! working array
  INTEGER :: IKB, IKE, JIJ,JK,IIJB,IIJE,IKT
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_THR_WTH2_O_DDRDZ',0,ZHOOK_HANDLE)
IKB=D%NKTB
IKE=D%NKTE
IIJE=D%NIJE
IIJB=D%NIJB
IKT=D%NKT
!

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZWORK1(JIJ, JK) = PBLL_O_E(JIJ, JK)*PETHETA(JIJ, JK)/PD(JIJ, JK)&
           *(-(1.+PREDR1(JIJ, JK))*PREDR1(JIJ, JK)/PD(JIJ, JK)&
           *(1.5+PREDTH1(JIJ, JK)+PREDR1(JIJ, JK))+(1.+2.*PREDR1(JIJ, JK)))
  END DO
END DO

CALL MZF_PHY(D,ZWORK1,ZWORK2)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    PD_M3_THR_WTH2_O_DDRDZ(JIJ, JK) = - 0.25*PLEPS(JIJ, JK)/PSQRT_TKE(JIJ, JK)& 
                                                    / CSTURB%XCTD*TURBN%XCTV * ZWORK2(JIJ, JK)
  END DO
END DO
!
PD_M3_THR_WTH2_O_DDRDZ(IIJB:IIJE,IKB-1)=PD_M3_THR_WTH2_O_DDRDZ(IIJB:IIJE,IKB)
PD_M3_THR_WTH2_O_DDRDZ(IIJB:IIJE,IKE+1)=PD_M3_THR_WTH2_O_DDRDZ(IIJB:IIJE,IKE)

!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_THR_WTH2_O_DDRDZ',1,ZHOOK_HANDLE)
END SUBROUTINE D_M3_THR_WTH2_O_DDRDZ
!----------------------------------------------------------------------------
SUBROUTINE M3_THR_W2TH(D,CSTURB,TURBN,PREDR1,PD,PLM,PLEPS,PTKE,PDRDZ,PM3_THR_W2TH)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  TYPE(TURB_t),                    INTENT(IN)    :: TURBN
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PLM
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PTKE
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PDRDZ
  REAL, DIMENSION(D%NIJT,D%NKT),INTENT(OUT) :: PM3_THR_W2TH
  REAL, DIMENSION(D%NIJT,D%NKT) :: ZWORK1,ZWORK2 ! working array
  INTEGER :: IKB, IKE, JIJ,JK,IIJB,IIJE,IKT
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_THR_W2TH',0,ZHOOK_HANDLE)
IKB=D%NKTB
IKE=D%NKTE
IIJE=D%NIJE
IIJB=D%NIJB
IKT=D%NKT
!

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZWORK1(JIJ, JK) = (1.+PREDR1(JIJ, JK))*PDRDZ(JIJ, JK)/PD(JIJ, JK)
  END DO
END DO

CALL MZF_PHY(D,ZWORK1,ZWORK2)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    PM3_THR_W2TH(JIJ, JK) = - 0.75*PLM(JIJ, JK)*PLEPS(JIJ, JK)& 
                                          / PTKE(JIJ, JK) * TURBN%XCTV * ZWORK2(JIJ, JK)
  END DO
END DO
!
PM3_THR_W2TH(IIJB:IIJE,IKB-1)=PM3_THR_W2TH(IIJB:IIJE,IKB)
PM3_THR_W2TH(IIJB:IIJE,IKE+1)=PM3_THR_W2TH(IIJB:IIJE,IKE)

!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_THR_W2TH',1,ZHOOK_HANDLE)
END SUBROUTINE M3_THR_W2TH
!----------------------------------------------------------------------------
SUBROUTINE D_M3_THR_W2TH_O_DDTDZ(D,CSTURB,TURBN,PREDTH1,PREDR1,PD,PLM,PLEPS,PTKE,PBLL_O_E,PDRDZ,PETHETA,PD_M3_THR_W2TH_O_DDTDZ)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  TYPE(TURB_t),                    INTENT(IN)    :: TURBN
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PLM
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PTKE
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PDRDZ
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PETHETA
  REAL, DIMENSION(D%NIJT,D%NKT),INTENT(OUT) :: PD_M3_THR_W2TH_O_DDTDZ
  REAL, DIMENSION(D%NIJT,D%NKT) :: ZWORK1,ZWORK2 ! working array
  INTEGER :: IKB, IKE, JIJ,JK,IIJB,IIJE,IKT
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_THR_W2TH_O_DDTDZ',0,ZHOOK_HANDLE)
IKB=D%NKTB
IKE=D%NKTE
IIJE=D%NIJE
IIJB=D%NIJB
IKT=D%NKT
!

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZWORK1(JIJ, JK) =  -PETHETA(JIJ, JK)*PBLL_O_E(JIJ, JK)*& 
    (1.+PREDR1(JIJ, JK))*PDRDZ(JIJ, JK)& 
    *(1.5+PREDTH1(JIJ, JK)+PREDR1(JIJ, JK))/PD(JIJ, JK)**2
  END DO
END DO

CALL MZF_PHY(D,ZWORK1,ZWORK2)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    PD_M3_THR_W2TH_O_DDTDZ(JIJ, JK) = - 0.75*PLM(JIJ, JK)*PLEPS(JIJ, JK)&
                                                    / PTKE(JIJ, JK) * TURBN%XCTV**2 * ZWORK1(JIJ, JK)
  END DO
END DO
!
PD_M3_THR_W2TH_O_DDTDZ(IIJB:IIJE,IKB-1)=PD_M3_THR_W2TH_O_DDTDZ(IIJB:IIJE,IKB)
PD_M3_THR_W2TH_O_DDTDZ(IIJB:IIJE,IKE+1)=PD_M3_THR_W2TH_O_DDTDZ(IIJB:IIJE,IKE)

!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_THR_W2TH_O_DDTDZ',1,ZHOOK_HANDLE)
END SUBROUTINE D_M3_THR_W2TH_O_DDTDZ
!----------------------------------------------------------------------------
SUBROUTINE D_M3_THR_W2TH_O_DDRDZ(D,CSTURB,TURBN,PREDTH1,PREDR1,PD,PLM,PLEPS,PTKE,PD_M3_THR_W2TH_O_DDRDZ)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  TYPE(TURB_t),                    INTENT(IN)    :: TURBN
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PLM
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PTKE
  REAL, DIMENSION(D%NIJT,D%NKT),INTENT(OUT) :: PD_M3_THR_W2TH_O_DDRDZ
  REAL, DIMENSION(D%NIJT,D%NKT) :: ZWORK1,ZWORK2 ! working array
  INTEGER :: IKB, IKE, JIJ,JK,IIJB,IIJE,IKT
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_THR_W2TH_O_DDRDZ',0,ZHOOK_HANDLE)
IKB=D%NKTB
IKE=D%NKTE
IIJE=D%NIJE
IIJB=D%NIJB
IKT=D%NKT
!

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZWORK1(JIJ, JK) =  -(1.+PREDR1(JIJ, JK))*PREDR1(JIJ, JK)&
    * (1.5+PREDTH1(JIJ, JK)+PREDR1(JIJ, JK))/PD(JIJ, JK)**2          &
            +(1.+2.*PREDR1(JIJ, JK))/PD(JIJ, JK)
  END DO
END DO

CALL MZF_PHY(D,ZWORK1,ZWORK2)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    PD_M3_THR_W2TH_O_DDRDZ(JIJ, JK) = - 0.75*PLM(JIJ, JK)*PLEPS(JIJ, JK)&
                                                    / PTKE(JIJ, JK) * TURBN%XCTV * ZWORK2(JIJ, JK)
  END DO
END DO
!
PD_M3_THR_W2TH_O_DDRDZ(IIJB:IIJE,IKB-1)=PD_M3_THR_W2TH_O_DDRDZ(IIJB:IIJE,IKB)
PD_M3_THR_W2TH_O_DDRDZ(IIJB:IIJE,IKE+1)=PD_M3_THR_W2TH_O_DDRDZ(IIJB:IIJE,IKE)

!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_THR_W2TH_O_DDRDZ',1,ZHOOK_HANDLE)
END SUBROUTINE D_M3_THR_W2TH_O_DDRDZ
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!
SUBROUTINE PSI3(D,CSTURB,TURBN,PREDR1,PREDTH1,PRED2R3,PRED2TH3,PRED2THR3,HTURBDIM,OUSERV,PPSI3)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  TYPE(TURB_t),                    INTENT(IN)    :: TURBN
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PRED2TH3
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PRED2R3
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PRED2THR3
  CHARACTER(LEN=4),       INTENT(IN) :: HTURBDIM  ! 1DIM or 3DIM turb. scheme
  LOGICAL,                INTENT(IN) :: OUSERV    ! flag to use vapor
  REAL, DIMENSION(D%NIJT,D%NKT),INTENT(OUT) :: PPSI3
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:PSI3',0,ZHOOK_HANDLE)
CALL PHI3(D,CSTURB,TURBN,PREDR1,PREDTH1,PRED2R3,PRED2TH3,PRED2THR3,HTURBDIM,OUSERV,PPSI3)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:PSI3',1,ZHOOK_HANDLE)
END SUBROUTINE PSI3
!----------------------------------------------------------------------------
SUBROUTINE D_PSI3DRDZ_O_DDRDZ(D,CSTURB,TURBN,PPSI3,PREDR1,PREDTH1,PRED2R3,PRED2THR3,HTURBDIM,OUSERV,PD_PSI3DRDZ_O_DDRDZ)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  TYPE(TURB_t),                    INTENT(IN)    :: TURBN
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PPSI3
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PRED2R3
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PRED2THR3
  CHARACTER(LEN=4),       INTENT(IN) :: HTURBDIM  ! 1DIM or 3DIM turb. scheme
  LOGICAL,                INTENT(IN) :: OUSERV    ! flag to use vapor
  REAL, DIMENSION(D%NIJT,D%NKT),INTENT(OUT) :: PD_PSI3DRDZ_O_DDRDZ

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_PSI3DRDZ_O_DDRDZ',0,ZHOOK_HANDLE)
CALL D_PHI3DTDZ_O_DDTDZ(D,CSTURB,TURBN,PPSI3,PREDR1,PREDTH1,PRED2R3,PRED2THR3,HTURBDIM,OUSERV,PD_PSI3DRDZ_O_DDRDZ)
!
!C'est ok?!
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_PSI3DRDZ_O_DDRDZ',1,ZHOOK_HANDLE)
END SUBROUTINE D_PSI3DRDZ_O_DDRDZ
!----------------------------------------------------------------------------
SUBROUTINE D_PSI3DTDZ_O_DDTDZ(D,CSTURB,TURBN,PPSI3,PREDR1,PREDTH1,PRED2R3,PRED2THR3,HTURBDIM,OUSERV,PD_PSI3DTDZ_O_DDTDZ)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  TYPE(TURB_t),                    INTENT(IN)    :: TURBN
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PPSI3
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PRED2R3
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PRED2THR3
  CHARACTER(LEN=4),       INTENT(IN) :: HTURBDIM  ! 1DIM or 3DIM turb. scheme
  LOGICAL,                INTENT(IN) :: OUSERV    ! flag to use vapor
  REAL, DIMENSION(D%NIJT,D%NKT),INTENT(OUT) :: PD_PSI3DTDZ_O_DDTDZ
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_PSI3DTDZ_O_DDTDZ',0,ZHOOK_HANDLE)
CALL D_PHI3DRDZ_O_DDRDZ(D,CSTURB,TURBN,PPSI3,PREDR1,PREDTH1,PRED2R3,PRED2THR3,HTURBDIM,OUSERV,PD_PSI3DTDZ_O_DDTDZ)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_PSI3DTDZ_O_DDTDZ',1,ZHOOK_HANDLE)
END SUBROUTINE D_PSI3DTDZ_O_DDTDZ
!----------------------------------------------------------------------------
SUBROUTINE D_PSI3DRDZ2_O_DDRDZ(D,CSTURB,TURBN,PPSI3,PREDR1,PREDTH1,PRED2R3,PRED2THR3,PDRDZ,HTURBDIM,OUSERV,PD_PSI3DRDZ2_O_DDRDZ)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  TYPE(TURB_t),                    INTENT(IN)    :: TURBN
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PPSI3
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PRED2R3
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PRED2THR3
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PDRDZ
  CHARACTER(LEN=4),       INTENT(IN) :: HTURBDIM  ! 1DIM or 3DIM turb. scheme
  LOGICAL,                INTENT(IN) :: OUSERV    ! flag to use vapor
  REAL, DIMENSION(D%NIJT,D%NKT),INTENT(OUT) :: PD_PSI3DRDZ2_O_DDRDZ
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_PSI3DRDZ2_O_DDRDZ',0,ZHOOK_HANDLE)
CALL D_PHI3DTDZ2_O_DDTDZ(D,CSTURB,TURBN,PPSI3,PREDR1,PREDTH1,PRED2R3,PRED2THR3,PDRDZ,HTURBDIM,OUSERV,PD_PSI3DRDZ2_O_DDRDZ)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_PSI3DRDZ2_O_DDRDZ',1,ZHOOK_HANDLE)
END SUBROUTINE D_PSI3DRDZ2_O_DDRDZ
!----------------------------------------------------------------------------
SUBROUTINE M3_WR_WR2(D,CSTURB,TURBN,PREDR1,PREDTH1,PD,PBLL_O_E,PEMOIST,PM3_WR_WR2)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  TYPE(TURB_t),                    INTENT(IN)    :: TURBN
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PEMOIST
  REAL, DIMENSION(D%NIJT,D%NKT),INTENT(OUT) :: PM3_WR_WR2
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_WR_WR2',0,ZHOOK_HANDLE)
CALL M3_WTH_WTH2(D,CSTURB,TURBN,PREDR1,PREDTH1,PD,PBLL_O_E,PEMOIST,PM3_WR_WR2)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_WR_WR2',1,ZHOOK_HANDLE)
END SUBROUTINE M3_WR_WR2
!----------------------------------------------------------------------------
SUBROUTINE D_M3_WR_WR2_O_DDRDZ(D,CSTURB,TURBN,PM3_WR_WR2,PREDR1,PREDTH1,PD,PBLL_O_E,PEMOIST,PD_M3_WR_WR2_O_DDRDZ)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  TYPE(TURB_t),                    INTENT(IN)    :: TURBN
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PM3_WR_WR2
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PEMOIST
  REAL, DIMENSION(D%NIJT,D%NKT),INTENT(OUT) :: PD_M3_WR_WR2_O_DDRDZ
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_WR_WR2_O_DDRDZ',0,ZHOOK_HANDLE)
CALL D_M3_WTH_WTH2_O_DDTDZ(D,CSTURB,TURBN,PM3_WR_WR2,PREDR1,PREDTH1,PD,PBLL_O_E,PEMOIST,PD_M3_WR_WR2_O_DDRDZ)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_WR_WR2_O_DDRDZ',1,ZHOOK_HANDLE)
END SUBROUTINE D_M3_WR_WR2_O_DDRDZ
!----------------------------------------------------------------------------
SUBROUTINE M3_WR_W2R(D,CSTURB,TURBN,PREDR1,PREDTH1,PD,PKEFF,PTKE,PM3_WR_W2R)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  TYPE(TURB_t),                    INTENT(IN)    :: TURBN
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PKEFF
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PTKE
  REAL, DIMENSION(D%NIJT,D%NKT),INTENT(OUT) :: PM3_WR_W2R
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_WR_W2R',0,ZHOOK_HANDLE)
CALL M3_WTH_W2TH(D,CSTURB,TURBN,PREDR1,PREDTH1,PD,PKEFF,PTKE,PM3_WR_W2R)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_WR_W2R',1,ZHOOK_HANDLE)
END SUBROUTINE M3_WR_W2R
!----------------------------------------------------------------------------
SUBROUTINE D_M3_WR_W2R_O_DDRDZ(D,CSTURB,TURBN,PREDR1,PREDTH1,PD,PBLL_O_E,PEMOIST,PKEFF,PTKE,PD_M3_WR_W2R_O_DDRDZ)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  TYPE(TURB_t),                    INTENT(IN)    :: TURBN
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PEMOIST
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PKEFF
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PTKE
  REAL, DIMENSION(D%NIJT,D%NKT),INTENT(OUT) :: PD_M3_WR_W2R_O_DDRDZ
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_WR_W2R_O_DDRDZ',0,ZHOOK_HANDLE)
CALL D_M3_WTH_W2TH_O_DDTDZ(D,CSTURB,TURBN,PREDR1,PREDTH1,PD,PBLL_O_E,PEMOIST,PKEFF,PTKE,PD_M3_WR_W2R_O_DDRDZ)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_WR_W2R_O_DDRDZ',1,ZHOOK_HANDLE)
END SUBROUTINE D_M3_WR_W2R_O_DDRDZ
!----------------------------------------------------------------------------
SUBROUTINE M3_WR_W2TH(D,CSTURB,TURBN,PD,PKEFF,PTKE,PBLL_O_E,PETHETA,PDRDZ,PM3_WR_W2TH)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  TYPE(TURB_t),                    INTENT(IN)    :: TURBN
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PKEFF
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PTKE
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PETHETA
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PDRDZ
  REAL, DIMENSION(D%NIJT,D%NKT),INTENT(OUT) :: PM3_WR_W2TH
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_WR_W2TH',0,ZHOOK_HANDLE)
CALL M3_WTH_W2R(D,CSTURB,TURBN,PD,PKEFF,PTKE,PBLL_O_E,PETHETA,PDRDZ,PM3_WR_W2TH)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_WR_W2TH',1,ZHOOK_HANDLE)
END SUBROUTINE M3_WR_W2TH
!----------------------------------------------------------------------------
SUBROUTINE D_M3_WR_W2TH_O_DDRDZ(D,CSTURB,TURBN,PREDR1,PREDTH1,PD,PKEFF,PTKE,PBLL_O_E,PETHETA,PD_M3_WR_W2TH_O_DDRDZ)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  TYPE(TURB_t),                    INTENT(IN)    :: TURBN
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PKEFF
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PTKE
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PETHETA
  REAL, DIMENSION(D%NIJT,D%NKT),INTENT(OUT) :: PD_M3_WR_W2TH_O_DDRDZ
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_WR_W2TH_O_DDRDZ',0,ZHOOK_HANDLE)
CALL D_M3_WTH_W2R_O_DDTDZ(D,CSTURB,TURBN,PREDR1,PREDTH1,PD,PKEFF,PTKE,PBLL_O_E,PETHETA,PD_M3_WR_W2TH_O_DDRDZ)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_WR_W2TH_O_DDRDZ',1,ZHOOK_HANDLE)
END SUBROUTINE D_M3_WR_W2TH_O_DDRDZ
!----------------------------------------------------------------------------
SUBROUTINE M3_WR_WTH2(D,CSTURB,TURBN,PD,PKEFF,PTKE,PSQRT_TKE,PBLL_O_E,PBETA,PLEPS,PETHETA,PDRDZ,PM3_WR_WTH2)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  TYPE(TURB_t),                    INTENT(IN)    :: TURBN
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PKEFF
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PTKE
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PBETA
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PETHETA
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PDRDZ
  REAL, DIMENSION(D%NIJT,D%NKT),INTENT(OUT) :: PM3_WR_WTH2
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_WR_WTH2',0,ZHOOK_HANDLE)
CALL M3_WTH_WR2(D,CSTURB,TURBN,PD,PKEFF,PTKE,PSQRT_TKE,PBLL_O_E,PBETA,PLEPS,PETHETA,PDRDZ,PM3_WR_WTH2)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_WR_WTH2',1,ZHOOK_HANDLE)
END SUBROUTINE M3_WR_WTH2
!----------------------------------------------------------------------------
SUBROUTINE D_M3_WR_WTH2_O_DDRDZ(D,CSTURB,TURBN,PREDR1,PREDTH1,PD,PKEFF,PTKE,PSQRT_TKE,PBLL_O_E,&
                               &PBETA,PLEPS,PETHETA,PD_M3_WR_WTH2_O_DDRDZ)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  TYPE(TURB_t),                    INTENT(IN)    :: TURBN
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PKEFF
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PTKE
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PBETA
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PETHETA
  REAL, DIMENSION(D%NIJT,D%NKT),INTENT(OUT) :: PD_M3_WR_WTH2_O_DDRDZ
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_WR_WTH2_O_DDRDZ',0,ZHOOK_HANDLE)
CALL D_M3_WTH_WR2_O_DDTDZ(D,CSTURB,TURBN,PREDR1,PREDTH1,PD,PKEFF,PTKE,PSQRT_TKE,PBLL_O_E,PBETA,PLEPS,PETHETA,PD_M3_WR_WTH2_O_DDRDZ)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_WR_WTH2_O_DDRDZ',1,ZHOOK_HANDLE)
END SUBROUTINE D_M3_WR_WTH2_O_DDRDZ
!----------------------------------------------------------------------------
SUBROUTINE M3_WR_WTHR(D,CSTURB,TURBN,PREDTH1,PD,PKEFF,PTKE,PSQRT_TKE,PBETA,PLEPS,PETHETA,PM3_WR_WTHR)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  TYPE(TURB_t),                    INTENT(IN)    :: TURBN
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PKEFF
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PTKE
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PBETA
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PETHETA
  REAL, DIMENSION(D%NIJT,D%NKT),INTENT(OUT) :: PM3_WR_WTHR
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_WR_WTHR',0,ZHOOK_HANDLE)
CALL M3_WTH_WTHR(D,CSTURB,TURBN,PREDTH1,PD,PKEFF,PTKE,PSQRT_TKE,PBETA,PLEPS,PETHETA,PM3_WR_WTHR)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_WR_WTHR',1,ZHOOK_HANDLE)
END SUBROUTINE M3_WR_WTHR
!----------------------------------------------------------------------------
SUBROUTINE D_M3_WR_WTHR_O_DDRDZ(D,CSTURB,TURBN,PM3_WR_WTHR,PREDR1,PREDTH1,PD,PBLL_O_E,PEMOIST,PD_M3_WR_WTHR_O_DDRDZ)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  TYPE(TURB_t),                    INTENT(IN)    :: TURBN
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PM3_WR_WTHR
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PEMOIST
  REAL, DIMENSION(D%NIJT,D%NKT),INTENT(OUT) :: PD_M3_WR_WTHR_O_DDRDZ
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_WR_WTHR_O_DDRDZ',0,ZHOOK_HANDLE)
CALL D_M3_WTH_WTHR_O_DDTDZ(D,CSTURB,TURBN,PM3_WR_WTHR,PREDR1,PREDTH1,PD,PBLL_O_E,PEMOIST,PD_M3_WR_WTHR_O_DDRDZ)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_WR_WTHR_O_DDRDZ',1,ZHOOK_HANDLE)
END SUBROUTINE D_M3_WR_WTHR_O_DDRDZ
!----------------------------------------------------------------------------
SUBROUTINE M3_R2_W2R(D,CSTURB,TURBN,PREDR1,PREDTH1,PD,PDRDZ,PLM,PLEPS,PTKE,PM3_R2_W2R)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  TYPE(TURB_t),                    INTENT(IN)    :: TURBN
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PDRDZ
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PLM
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PTKE
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT) :: PM3_R2_W2R
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_R2_W2R',0,ZHOOK_HANDLE)
CALL M3_TH2_W2TH(D,CSTURB,TURBN,PREDR1,PREDTH1,PD,PDRDZ,PLM,PLEPS,PTKE,PM3_R2_W2R)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_R2_W2R',1,ZHOOK_HANDLE)
END SUBROUTINE M3_R2_W2R
!----------------------------------------------------------------------------
SUBROUTINE D_M3_R2_W2R_O_DDRDZ(D,CSTURB,TURBN,PREDR1,PREDTH1,PD,PLM,PLEPS,PTKE,OUSERV,PD_M3_R2_W2R_O_DDRDZ)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  TYPE(TURB_t),                    INTENT(IN)    :: TURBN
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PLM
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PTKE
  LOGICAL,                INTENT(IN) :: OUSERV
  REAL, DIMENSION(D%NIJT,D%NKT),INTENT(OUT) :: PD_M3_R2_W2R_O_DDRDZ
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_R2_W2R_O_DDRDZ',0,ZHOOK_HANDLE)
CALL D_M3_TH2_W2TH_O_DDTDZ(D,CSTURB,TURBN,PREDR1,PREDTH1,PD,PLM,PLEPS,PTKE,OUSERV,PD_M3_R2_W2R_O_DDRDZ)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_R2_W2R_O_DDRDZ',1,ZHOOK_HANDLE)
END SUBROUTINE D_M3_R2_W2R_O_DDRDZ
!----------------------------------------------------------------------------
SUBROUTINE M3_R2_WR2(D,CSTURB,TURBN,PREDR1,PREDTH1,PD,PLEPS,PSQRT_TKE,PM3_R2_WR2)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  TYPE(TURB_t),                    INTENT(IN)    :: TURBN
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(D%NIJT,D%NKT),INTENT(OUT) :: PM3_R2_WR2
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_R2_WR2',0,ZHOOK_HANDLE)
CALL M3_TH2_WTH2(D,CSTURB,TURBN,PREDR1,PREDTH1,PD,PLEPS,PSQRT_TKE,PM3_R2_WR2)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_R2_WR2',1,ZHOOK_HANDLE)
END SUBROUTINE M3_R2_WR2
!----------------------------------------------------------------------------
SUBROUTINE D_M3_R2_WR2_O_DDRDZ(D,CSTURB,TURBN,PREDR1,PREDTH1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PEMOIST,PD_M3_R2_WR2_O_DDRDZ)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  TYPE(TURB_t),                    INTENT(IN)    :: TURBN
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PEMOIST
  REAL, DIMENSION(D%NIJT,D%NKT),INTENT(OUT) :: PD_M3_R2_WR2_O_DDRDZ
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_R2_WR2_O_DDRDZ',0,ZHOOK_HANDLE)
CALL D_M3_TH2_WTH2_O_DDTDZ(D,CSTURB,TURBN,PREDR1,PREDTH1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PEMOIST,PD_M3_R2_WR2_O_DDRDZ)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_R2_WR2_O_DDRDZ',1,ZHOOK_HANDLE)
END SUBROUTINE D_M3_R2_WR2_O_DDRDZ
!----------------------------------------------------------------------------
SUBROUTINE M3_R2_W2TH(D,CSTURB,TURBN,PD,PLM,PLEPS,PTKE,PBLL_O_E,PETHETA,PDRDZ,PM3_R2_W2TH)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  TYPE(TURB_t),                    INTENT(IN)    :: TURBN
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PLM
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PTKE
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PETHETA
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PDRDZ
  REAL, DIMENSION(D%NIJT,D%NKT),INTENT(OUT) :: PM3_R2_W2TH
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_R2_W2TH',0,ZHOOK_HANDLE)
CALL M3_TH2_W2R(D,CSTURB,TURBN,PD,PLM,PLEPS,PTKE,PBLL_O_E,PETHETA,PDRDZ,PM3_R2_W2TH)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_R2_W2TH',1,ZHOOK_HANDLE)
END SUBROUTINE M3_R2_W2TH
!----------------------------------------------------------------------------
SUBROUTINE D_M3_R2_W2TH_O_DDRDZ(D,CSTURB,TURBN,PREDR1,PREDTH1,PD,PLM,PLEPS,PTKE,PBLL_O_E,PETHETA,PDRDZ,PD_M3_R2_W2TH_O_DDRDZ)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  TYPE(TURB_t),                    INTENT(IN)    :: TURBN
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PLM
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PTKE
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PETHETA
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PDRDZ
  REAL, DIMENSION(D%NIJT,D%NKT),INTENT(OUT) :: PD_M3_R2_W2TH_O_DDRDZ
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_R2_W2TH_O_DDRDZ',0,ZHOOK_HANDLE)
CALL D_M3_TH2_W2R_O_DDTDZ(D,CSTURB,TURBN,PREDR1,PREDTH1,PD,PLM,PLEPS,PTKE,PBLL_O_E,PETHETA,PDRDZ,PD_M3_R2_W2TH_O_DDRDZ)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_R2_W2TH_O_DDRDZ',1,ZHOOK_HANDLE)
END SUBROUTINE D_M3_R2_W2TH_O_DDRDZ
!----------------------------------------------------------------------------
SUBROUTINE M3_R2_WTH2(D,CSTURB,TURBN,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PETHETA,PDRDZ,PM3_R2_WTH2)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  TYPE(TURB_t),                    INTENT(IN)    :: TURBN
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PETHETA
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PDRDZ
  REAL, DIMENSION(D%NIJT,D%NKT),INTENT(OUT) :: PM3_R2_WTH2
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_R2_WTH2',0,ZHOOK_HANDLE)
CALL M3_TH2_WR2(D,CSTURB,TURBN,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PETHETA,PDRDZ,PM3_R2_WTH2)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_R2_WTH2',1,ZHOOK_HANDLE)
END SUBROUTINE M3_R2_WTH2
!----------------------------------------------------------------------------
SUBROUTINE D_M3_R2_WTH2_O_DDRDZ(D,CSTURB,TURBN,PREDR1,PREDTH1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PETHETA,PDRDZ,PD_M3_R2_WTH2_O_DDRDZ)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  TYPE(TURB_t),                    INTENT(IN)    :: TURBN
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PETHETA
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PDRDZ
  REAL, DIMENSION(D%NIJT,D%NKT),INTENT(OUT) :: PD_M3_R2_WTH2_O_DDRDZ
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_R2_WTH2_O_DDRDZ',0,ZHOOK_HANDLE)
CALL D_M3_TH2_WR2_O_DDTDZ(D,CSTURB,TURBN,PREDR1,PREDTH1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PETHETA,PDRDZ,PD_M3_R2_WTH2_O_DDRDZ)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_R2_WTH2_O_DDRDZ',1,ZHOOK_HANDLE)
END SUBROUTINE D_M3_R2_WTH2_O_DDRDZ
!----------------------------------------------------------------------------
SUBROUTINE M3_R2_WTHR(D,CSTURB,TURBN,PREDTH1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PETHETA,PDRDZ,PM3_R2_WTHR)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  TYPE(TURB_t),                    INTENT(IN)    :: TURBN
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PETHETA
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PDRDZ
  REAL, DIMENSION(D%NIJT,D%NKT),INTENT(OUT) :: PM3_R2_WTHR
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_R2_WTHR',0,ZHOOK_HANDLE)
CALL M3_TH2_WTHR(D,CSTURB,TURBN,PREDTH1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PETHETA,PDRDZ,PM3_R2_WTHR)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_R2_WTHR',1,ZHOOK_HANDLE)
END SUBROUTINE M3_R2_WTHR
!----------------------------------------------------------------------------
SUBROUTINE D_M3_R2_WTHR_O_DDRDZ(D,CSTURB,TURBN,PREDR1,PREDTH1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PETHETA,PDRDZ,PD_M3_R2_WTHR_O_DDRDZ)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  TYPE(TURB_t),                    INTENT(IN)    :: TURBN
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PETHETA
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PDRDZ
  REAL, DIMENSION(D%NIJT,D%NKT),INTENT(OUT) :: PD_M3_R2_WTHR_O_DDRDZ
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_R2_WTHR_O_DDRDZ',0,ZHOOK_HANDLE)
CALL D_M3_TH2_WTHR_O_DDTDZ(D,CSTURB,TURBN,PREDR1,PREDTH1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PETHETA,PDRDZ,PD_M3_R2_WTHR_O_DDRDZ)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_R2_WTHR_O_DDRDZ',1,ZHOOK_HANDLE)
END SUBROUTINE D_M3_R2_WTHR_O_DDRDZ
!----------------------------------------------------------------------------
SUBROUTINE D_M3_THR_WTHR_O_DDRDZ(D,CSTURB,TURBN,PREDR1,PREDTH1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PEMOIST,PD_M3_THR_WTHR_O_DDRDZ)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  TYPE(TURB_t),                    INTENT(IN)    :: TURBN
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PEMOIST
  REAL, DIMENSION(D%NIJT,D%NKT),INTENT(OUT) :: PD_M3_THR_WTHR_O_DDRDZ
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_THR_WTHR_O_DDRDZ',0,ZHOOK_HANDLE)
CALL D_M3_THR_WTHR_O_DDTDZ(D,CSTURB,TURBN,PREDR1,PREDTH1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PEMOIST,PD_M3_THR_WTHR_O_DDRDZ)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_THR_WTHR_O_DDRDZ',1,ZHOOK_HANDLE)
END SUBROUTINE D_M3_THR_WTHR_O_DDRDZ
!----------------------------------------------------------------------------
SUBROUTINE M3_THR_WR2(D,CSTURB,TURBN,PREDTH1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PEMOIST,PDTDZ,PM3_THR_WR2)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  TYPE(TURB_t),                    INTENT(IN)    :: TURBN
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PEMOIST
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PDTDZ
  REAL, DIMENSION(D%NIJT,D%NKT),INTENT(OUT) :: PM3_THR_WR2
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_THR_WR2',0,ZHOOK_HANDLE)
CALL M3_THR_WTH2(D,CSTURB,TURBN,PREDTH1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PEMOIST,PDTDZ,PM3_THR_WR2)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_THR_WR2',1,ZHOOK_HANDLE)
END SUBROUTINE M3_THR_WR2
!----------------------------------------------------------------------------
SUBROUTINE D_M3_THR_WR2_O_DDRDZ(D,CSTURB,TURBN,PREDR1,PREDTH1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PEMOIST,PDTDZ,PD_M3_THR_WR2_O_DDRDZ)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  TYPE(TURB_t),                    INTENT(IN)    :: TURBN
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PEMOIST
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PDTDZ
  REAL, DIMENSION(D%NIJT,D%NKT),INTENT(OUT) :: PD_M3_THR_WR2_O_DDRDZ
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_THR_WR2_O_DDRDZ',0,ZHOOK_HANDLE)
CALL D_M3_THR_WTH2_O_DDTDZ(D,CSTURB,TURBN,PREDR1,PREDTH1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PEMOIST,PDTDZ,PD_M3_THR_WR2_O_DDRDZ)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_THR_WR2_O_DDRDZ',1,ZHOOK_HANDLE)
END SUBROUTINE D_M3_THR_WR2_O_DDRDZ
!----------------------------------------------------------------------------
SUBROUTINE D_M3_THR_WR2_O_DDTDZ(D,CSTURB,TURBN,PREDR1,PREDTH1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PEMOIST,PD_M3_THR_WR2_O_DDTDZ)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  TYPE(TURB_t),                    INTENT(IN)    :: TURBN
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PEMOIST
  REAL, DIMENSION(D%NIJT,D%NKT),INTENT(OUT) :: PD_M3_THR_WR2_O_DDTDZ
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_THR_WR2_O_DDTDZ',0,ZHOOK_HANDLE)
CALL D_M3_THR_WTH2_O_DDRDZ(D,CSTURB,TURBN,PREDR1,PREDTH1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PEMOIST,PD_M3_THR_WR2_O_DDTDZ)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_THR_WR2_O_DDTDZ',1,ZHOOK_HANDLE)
END SUBROUTINE D_M3_THR_WR2_O_DDTDZ
!----------------------------------------------------------------------------
SUBROUTINE M3_THR_W2R(D,CSTURB,TURBN,PREDTH1,PD,PLM,PLEPS,PTKE,PDTDZ,PM3_THR_W2R)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  TYPE(TURB_t),                    INTENT(IN)    :: TURBN
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PLM
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PTKE
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PDTDZ
  REAL, DIMENSION(D%NIJT,D%NKT),INTENT(OUT) :: PM3_THR_W2R
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_THR_W2R',0,ZHOOK_HANDLE)
CALL M3_THR_W2TH(D,CSTURB,TURBN,PREDTH1,PD,PLM,PLEPS,PTKE,PDTDZ,PM3_THR_W2R)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_THR_W2R',1,ZHOOK_HANDLE)
END SUBROUTINE M3_THR_W2R
!----------------------------------------------------------------------------
SUBROUTINE D_M3_THR_W2R_O_DDRDZ(D,CSTURB,TURBN,PREDR1,PREDTH1,PD,PLM,PLEPS,PTKE,PBLL_O_E,PDTDZ,PEMOIST,PD_M3_THR_W2R_O_DDRDZ)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  TYPE(TURB_t),                    INTENT(IN)    :: TURBN
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PLM
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PTKE
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PDTDZ
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PEMOIST
  REAL, DIMENSION(D%NIJT,D%NKT),INTENT(OUT) :: PD_M3_THR_W2R_O_DDRDZ
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_THR_W2R_O_DDRDZ',0,ZHOOK_HANDLE)
CALL D_M3_THR_W2TH_O_DDTDZ(D,CSTURB,TURBN,PREDR1,PREDTH1,PD,PLM,PLEPS,PTKE,PBLL_O_E,PDTDZ,PEMOIST,PD_M3_THR_W2R_O_DDRDZ)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_THR_W2R_O_DDRDZ',1,ZHOOK_HANDLE)
END SUBROUTINE D_M3_THR_W2R_O_DDRDZ
!----------------------------------------------------------------------------
SUBROUTINE D_M3_THR_W2R_O_DDTDZ(D,CSTURB,TURBN,PREDR1,PREDTH1,PD,PLM,PLEPS,PTKE,PD_M3_THR_W2R_O_DDTDZ)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  TYPE(TURB_t),                    INTENT(IN)    :: TURBN
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PLM
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PTKE
  REAL, DIMENSION(D%NIJT,D%NKT),INTENT(OUT) :: PD_M3_THR_W2R_O_DDTDZ
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_THR_W2R_O_DDTDZ',0,ZHOOK_HANDLE)
CALL D_M3_THR_W2TH_O_DDRDZ(D,CSTURB,TURBN,PREDR1,PREDTH1,PD,PLM,PLEPS,PTKE,PD_M3_THR_W2R_O_DDTDZ)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_THR_W2R_O_DDTDZ',1,ZHOOK_HANDLE)
END SUBROUTINE D_M3_THR_W2R_O_DDTDZ
!----------------------------------------------------------------------------
!
END MODULE MODE_PRANDTL

