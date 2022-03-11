!MNH_LIC Copyright 1994-2020 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!    #################### 
     MODULE MODE_PRANDTL
     USE PARKIND1, ONLY : JPRB
     USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!    #################### 
!
!* modification 08/2010  V. Masson  smoothing of the discontinuity in functions 
!                                   used for implicitation of exchange coefficients
!               05/2020   V. Masson and C. Lac : bug in D_PHI3DTDZ2_O_DDTDZ
!
USE MODD_CTURB,      ONLY : CSTURB_t
USE MODD_PARAMETERS, ONLY : JPVEXT_TURB
!
USE MODI_SHUMAN, ONLY: MZM, MZF
IMPLICIT NONE
!----------------------------------------------------------------------------
CONTAINS
!----------------------------------------------------------------------------
      SUBROUTINE PRANDTL(CST,CSTURB,KKA,KKU,KKL,KRR,KRRI,OTURB_DIAG,&
                         HTURBDIM,OOCEAN,OHARAT,O2D,           &
                         TPFILE,                               &
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
!!       CSTURB%XCTV,XCPR2    : constants for the turbulent prandtl numbers
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
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!
USE MODD_CST, ONLY : CST_t
USE MODD_CTURB, ONLY : CSTURB_t
USE MODD_FIELD,          ONLY: TFIELDDATA, TYPEREAL
USE MODD_IO,             ONLY: TFILEDATA
USE MODD_PARAMETERS
!
USE MODI_GRADIENT_M
USE MODE_EMOIST
USE MODE_ETHETA
USE MODI_SHUMAN, ONLY: MZM
USE MODE_IO_FIELD_WRITE, ONLY: IO_FIELD_WRITE
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
TYPE(CST_t),                  INTENT(IN)    :: CST
TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
INTEGER,                INTENT(IN)   :: KKA           !near ground array index  
INTEGER,                INTENT(IN)   :: KKU           !uppest atmosphere array index
INTEGER,                INTENT(IN)   :: KKL           !vert. levels type 1=MNH -1=ARO

INTEGER,                INTENT(IN)   :: KRR           ! number of moist var.
INTEGER,                INTENT(IN)   :: KRRI          ! number of ice var.
!
LOGICAL,                INTENT(IN)   ::  OTURB_DIAG   ! switch to write some
                                 ! diagnostic fields in the syncronous FM-file
LOGICAL,                INTENT(IN)   ::  OOCEAN       ! switch for Ocean model version
LOGICAL,                INTENT(IN)   ::  OHARAT
LOGICAL, INTENT(IN) :: O2D               ! Logical for 2D model version (modd_conf)
CHARACTER(LEN=4),       INTENT(IN)   ::  HTURBDIM     ! Kind of turbulence param.
TYPE(TFILEDATA),        INTENT(IN)   ::  TPFILE       ! Output file
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PDXX,PDYY,PDZZ,PDZX,PDZY
                                                  ! metric coefficients
!
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PTHVREF  ! Virtual Potential Temp.
                                                  ! of the reference state
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PLOCPEXNM ! Lv(T)/Cp/Exner at t-1 
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PATHETA      ! coefficients between 
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PAMOIST      ! s and Thetal and Rnp
!
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PLM      ! Turbulent Mixing length
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PLEPS    ! Dissipative length
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PTHLM,PTKEM! Conservative Potential 
                                                  ! Temperature and TKE at t-1
REAL, DIMENSION(:,:,:,:), INTENT(IN) ::  PRM      ! Mixing ratios at  t-1
                                                  ! with PRM(:,:,:,1) = cons.
                                                  ! mixing ratio
REAL, DIMENSION(:,:,:,:), INTENT(IN) ::  PSVM     ! Scalars at t-1      
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PSRCM
                                  ! s'r'c/2Sigma_s2 at t-1 multiplied by Lambda_3
!
!
REAL, DIMENSION(:,:,:), INTENT(OUT)  ::  PREDTH1 ! Redelsperger number R_theta
REAL, DIMENSION(:,:,:), INTENT(OUT)  ::  PREDR1  ! Redelsperger number R_q
REAL, DIMENSION(:,:,:), INTENT(OUT)  ::  PRED2TH3 ! Redelsperger number R*2_theta
REAL, DIMENSION(:,:,:), INTENT(OUT)  ::  PRED2R3  ! Redelsperger number R*2_q
REAL, DIMENSION(:,:,:), INTENT(OUT)  ::  PRED2THR3! Redelsperger number R*2_thq
REAL, DIMENSION(:,:,:,:), INTENT(OUT)::  PREDS1   ! Redelsperger number R_s
REAL, DIMENSION(:,:,:,:), INTENT(OUT)::  PRED2THS3! Redelsperger number R*2_thsv
REAL, DIMENSION(:,:,:,:), INTENT(OUT)::  PRED2RS3 ! Redelsperger number R*2_qsv
REAL, DIMENSION(:,:,:), INTENT(OUT)  ::  PBLL_O_E! beta*Lk*Leps/tke
REAL, DIMENSION(:,:,:), INTENT(OUT)  ::  PETHETA ! coefficient E_theta
REAL, DIMENSION(:,:,:), INTENT(OUT)  ::  PEMOIST ! coefficient E_moist
!
!
!       0.2  declaration of local variables
!
REAL, DIMENSION(SIZE(PTHLM,1),SIZE(PTHLM,2),SIZE(PTHLM,3)) ::  &
                  ZW1, ZW2, ZW3
!                                                 working variables
!                                                     
INTEGER :: IKB      ! vertical index value for the first inner mass point
INTEGER :: IKE      ! vertical index value for the last inner mass point
INTEGER             :: IRESP        ! Return code of FM routines
INTEGER             :: ILENG        ! Length of the data field in LFIFM file
INTEGER             :: IGRID        ! C-grid indicator in LFIFM file
INTEGER             :: ILENCH       ! Length of comment string in LFIFM file
CHARACTER (LEN=100) :: YCOMMENT     ! comment string in LFIFM file
CHARACTER (LEN=16)  :: YRECFM       ! Name of the desired field in LFIFM file
INTEGER::  ISV                      ! number of scalar variables       
INTEGER::  JSV                      ! loop index for the scalar variables  

INTEGER :: JLOOP
REAL    :: ZMINVAL
TYPE(TFIELDDATA)  :: TZFIELD
! ---------------------------------------------------------------------------
!
!*      1.  DEFAULT VALUES,  1D REDELSPERGER NUMBERS 
!           ----------------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('PRANDTL',0,ZHOOK_HANDLE)

IF (OHARAT) THEN
PREDTH1(:,:,:)=0.
PREDR1(:,:,:)=0.
PRED2TH3(:,:,:)=0.
PRED2R3(:,:,:)=0.
PRED2THR3(:,:,:)=0.
PREDS1(:,:,:,:)=0.
PRED2THS3(:,:,:,:)=0.
PRED2RS3(:,:,:,:)=0.
PBLL_O_E(:,:,:)=0.
ENDIF
!
IKB = KKA+JPVEXT_TURB*KKL
IKE = KKU-JPVEXT_TURB*KKL 
ILENG=SIZE(PTHLM,1)*SIZE(PTHLM,2)*SIZE(PTHLM,3)
ISV  =SIZE(PSVM,4)
!
PETHETA(:,:,:) = MZM(ETHETA(KRR,KRRI,PTHLM,PRM,PLOCPEXNM,PATHETA,PSRCM,OOCEAN), KKA, KKU, KKL)
PEMOIST(:,:,:) = MZM(EMOIST(KRR,KRRI,PTHLM,PRM,PLOCPEXNM,PAMOIST,PSRCM,OOCEAN), KKA, KKU, KKL)
PETHETA(:,:,KKA) = 2.*PETHETA(:,:,IKB) - PETHETA(:,:,IKB+KKL)
PEMOIST(:,:,KKA) = 2.*PEMOIST(:,:,IKB) - PEMOIST(:,:,IKB+KKL)
!
!---------------------------------------------------------------------------
IF (.NOT. OHARAT) THEN
!
!          1.3 1D Redelsperger numbers
!
PBLL_O_E(:,:,:) = MZM(CST%XG / PTHVREF(:,:,:) * PLM(:,:,:) * PLEPS(:,:,:) / PTKEM(:,:,:), KKA, KKU, KKL)
IF (KRR /= 0) THEN                ! moist case
  PREDTH1(:,:,:)= CSTURB%XCTV*PBLL_O_E(:,:,:) * PETHETA(:,:,:) * &
                   & GZ_M_W(KKA, KKU, KKL,PTHLM,PDZZ)
  PREDR1(:,:,:) = CSTURB%XCTV*PBLL_O_E(:,:,:) * PEMOIST(:,:,:) * &
                   & GZ_M_W(KKA, KKU, KKL,PRM(:,:,:,1),PDZZ)
ELSE                              ! dry case
  PREDTH1(:,:,:)= CSTURB%XCTV*PBLL_O_E(:,:,:)  * GZ_M_W(KKA, KKU, KKL,PTHLM,PDZZ)
  PREDR1(:,:,:) = 0.
END IF
!
!       3. Limits on 1D Redelperger numbers
!          --------------------------------
!
ZMINVAL = (1.-1./CSTURB%XPHI_LIM)
!
ZW1 = 1.
ZW2 = 1.
!
WHERE (PREDTH1+PREDR1<-ZMINVAL)
  ZW1 = (-ZMINVAL) / (PREDTH1+PREDR1)
END WHERE
!
WHERE (PREDTH1<-ZMINVAL)
  ZW2 = (-ZMINVAL) / (PREDTH1)
END WHERE
ZW2 = MIN(ZW1,ZW2)
!
ZW1 = 1.
WHERE (PREDR1<-ZMINVAL)
  ZW1 = (-ZMINVAL) / (PREDR1)
END WHERE
ZW1 = MIN(ZW2,ZW1)
!
!
!       3. Modification of Mixing length and dissipative length
!          ----------------------------------------------------
!
PBLL_O_E(:,:,:) = PBLL_O_E(:,:,:) * ZW1(:,:,:)
PREDTH1 (:,:,:) = PREDTH1 (:,:,:) * ZW1(:,:,:)
PREDR1  (:,:,:) = PREDR1  (:,:,:) * ZW1(:,:,:)
!
!       4. Threshold for very small (in absolute value) Redelperger numbers
!          ----------------------------------------------------------------
!
ZW2=SIGN(1.,PREDTH1(:,:,:))
PREDTH1(:,:,:)= ZW2(:,:,:) * MAX(1.E-30, ZW2(:,:,:)*PREDTH1(:,:,:))
!
IF (KRR /= 0) THEN                ! dry case
  ZW2=SIGN(1.,PREDR1(:,:,:))
  PREDR1(:,:,:)= ZW2(:,:,:) * MAX(1.E-30, ZW2(:,:,:)*PREDR1(:,:,:))
END IF
!
!
!---------------------------------------------------------------------------
!
!          For the scalar variables
DO JSV=1,ISV
  PREDS1(:,:,:,JSV)=CSTURB%XCTV*PBLL_O_E(:,:,:)*GZ_M_W(KKA, KKU, KKL,PSVM(:,:,:,JSV),PDZZ)
END DO
!
DO JSV=1,ISV
  ZW2=SIGN(1.,PREDS1(:,:,:,JSV))
  PREDS1(:,:,:,JSV)= ZW2(:,:,:) * MAX(1.E-30, ZW2(:,:,:)*PREDS1(:,:,:,JSV))
END DO
!
!---------------------------------------------------------------------------
!
!*      2.  3D REDELSPERGER NUMBERS
!           ------------------------
!
IF(HTURBDIM=='1DIM') THEN        ! 1D case
!
!
  PRED2TH3(:,:,:)  = PREDTH1(:,:,:)**2
!
  PRED2R3(:,:,:)   = PREDR1(:,:,:) **2
!
  PRED2THR3(:,:,:) = PREDTH1(:,:,:) * PREDR1(:,:,:)
!
ELSE IF (O2D) THEN                      ! 3D case in a 2D model
!
  IF (KRR /= 0) THEN                 ! moist 3D case
    PRED2TH3(:,:,:)= PREDTH1(:,:,:)**2+(CSTURB%XCTV*PBLL_O_E(:,:,:)*PETHETA(:,:,:) )**2 * &
      MZM(GX_M_M(PTHLM,PDXX,PDZZ,PDZX, KKA, KKU, KKL)**2, KKA, KKU, KKL)
    PRED2TH3(:,:,IKB)=PRED2TH3(:,:,IKB+KKL)
!
    PRED2R3(:,:,:)= PREDR1(:,:,:)**2 + (CSTURB%XCTV*PBLL_O_E(:,:,:)*PEMOIST(:,:,:))**2 * &
        MZM(GX_M_M(PRM(:,:,:,1),PDXX,PDZZ,PDZX, KKA, KKU, KKL)**2, KKA, KKU, KKL)
    PRED2R3(:,:,IKB)=PRED2R3(:,:,IKB+KKL)
!
    PRED2THR3(:,:,:)= PREDR1(:,:,:) * PREDTH1(:,:,:) +  CSTURB%XCTV**2*PBLL_O_E(:,:,:)**2 *   &
                  PEMOIST(:,:,:) * PETHETA(:,:,:) *                         &
      MZM(GX_M_M(PRM(:,:,:,1),PDXX,PDZZ,PDZX, KKA, KKU, KKL)*     &
                     GX_M_M(PTHLM,PDXX,PDZZ,PDZX, KKA, KKU, KKL), KKA, KKU, KKL)
    PRED2THR3(:,:,IKB)=PRED2THR3(:,:,IKB+KKL)
!
  ELSE                 ! dry 3D case in a 2D model
    PRED2TH3(:,:,:) = PREDTH1(:,:,:)**2 +  CSTURB%XCTV**2*PBLL_O_E(:,:,:)**2 *     &
      MZM(GX_M_M(PTHLM,PDXX,PDZZ,PDZX, KKA, KKU, KKL)**2, KKA, KKU, KKL)
    PRED2TH3(:,:,IKB)=PRED2TH3(:,:,IKB+KKL)
!
    PRED2R3(:,:,:) = 0.
!
    PRED2THR3(:,:,:) = 0.
!
  END IF
!
ELSE                                 ! 3D case in a 3D model
!
  IF (KRR /= 0) THEN                 ! moist 3D case
    PRED2TH3(:,:,:)= PREDTH1(:,:,:)**2 +  ( CSTURB%XCTV*PBLL_O_E(:,:,:)*PETHETA(:,:,:) )**2 * &
      MZM(GX_M_M(PTHLM,PDXX,PDZZ,PDZX, KKA, KKU, KKL)**2 &
      + GY_M_M(PTHLM,PDYY,PDZZ,PDZY, KKA, KKU, KKL)**2, KKA, KKU, KKL)
    PRED2TH3(:,:,IKB)=PRED2TH3(:,:,IKB+KKL)
!
    PRED2R3(:,:,:)= PREDR1(:,:,:)**2 + (CSTURB%XCTV*PBLL_O_E(:,:,:)*PEMOIST(:,:,:))**2 *      &
        MZM(GX_M_M(PRM(:,:,:,1),PDXX,PDZZ,PDZX, KKA, KKU, KKL)**2 + &
        GY_M_M(PRM(:,:,:,1),PDYY,PDZZ,PDZY, KKA, KKU, KKL)**2, KKA, KKU, KKL)
    PRED2R3(:,:,IKB)=PRED2R3(:,:,IKB+KKL)
!
    PRED2THR3(:,:,:)= PREDR1(:,:,:) * PREDTH1(:,:,:) +  CSTURB%XCTV**2*PBLL_O_E(:,:,:)**2 *   &
         PEMOIST(:,:,:) * PETHETA(:,:,:) *                            &
         MZM(GX_M_M(PRM(:,:,:,1),PDXX,PDZZ,PDZX, KKA, KKU, KKL)*   &
         GX_M_M(PTHLM,PDXX,PDZZ,PDZX, KKA, KKU, KKL)+                           &
         GY_M_M(PRM(:,:,:,1),PDYY,PDZZ,PDZY, KKA, KKU, KKL)*                    &
         GY_M_M(PTHLM,PDYY,PDZZ,PDZY, KKA, KKU, KKL), KKA, KKU, KKL)
    PRED2THR3(:,:,IKB)=PRED2THR3(:,:,IKB+KKL)
!
  ELSE                 ! dry 3D case in a 3D model
    PRED2TH3(:,:,:) = PREDTH1(:,:,:)**2 +  CSTURB%XCTV**2*PBLL_O_E(:,:,:)**2 *                &
      MZM(GX_M_M(PTHLM,PDXX,PDZZ,PDZX, KKA, KKU, KKL)**2 &
      + GY_M_M(PTHLM,PDYY,PDZZ,PDZY, KKA, KKU, KKL)**2, KKA, KKU, KKL)
    PRED2TH3(:,:,IKB)=PRED2TH3(:,:,IKB+KKL)
!
    PRED2R3(:,:,:) = 0.
!
    PRED2THR3(:,:,:) = 0.
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
DO JSV=1,ISV
!
  IF(HTURBDIM=='1DIM') THEN
!        1D case
    PRED2THS3(:,:,:,JSV)  = PREDS1(:,:,:,JSV) * PREDTH1(:,:,:)
    IF (KRR /= 0) THEN
      PRED2RS3(:,:,:,JSV)   = PREDR1(:,:,:) *PREDS1(:,:,:,JSV)
    ELSE
      PRED2RS3(:,:,:,JSV)   = 0.
    END IF
!
  ELSE  IF (O2D) THEN ! 3D case in a 2D model
!
    IF (KRR /= 0) THEN
      ZW1 = MZM((CST%XG / PTHVREF * PLM * PLEPS / PTKEM)**2, KKA, KKU, KKL) *PETHETA
    ELSE
      ZW1 = MZM((CST%XG / PTHVREF * PLM * PLEPS / PTKEM)**2, KKA, KKU, KKL)
    END IF
    PRED2THS3(:,:,:,JSV) = PREDTH1(:,:,:) * PREDS1(:,:,:,JSV)   +        &
                       ZW1*                                              &
                       MZM(GX_M_M(PSVM(:,:,:,JSV),PDXX,PDZZ,PDZX, KKA, KKU, KKL)*       &
                           GX_M_M(PTHLM,PDXX,PDZZ,PDZX, KKA, KKU, KKL),                 &
                           KKA, KKU, KKL)
!
    IF (KRR /= 0) THEN
      PRED2RS3(:,:,:,JSV) = PREDR1(:,:,:) * PREDS1(:,:,:,JSV)   +        &
                       ZW1 * PEMOIST *                                   &
                       MZM(GX_M_M(PSVM(:,:,:,JSV),PDXX,PDZZ,PDZX, KKA, KKU, KKL)*       &
                           GX_M_M(PRM(:,:,:,1),PDXX,PDZZ,PDZX, KKA, KKU, KKL),          &
                           KKA, KKU, KKL)
    ELSE
      PRED2RS3(:,:,:,JSV) = 0.
    END IF
!
  ELSE ! 3D case in a 3D model
!
    IF (KRR /= 0) THEN
      ZW1 = MZM((CST%XG / PTHVREF * PLM * PLEPS / PTKEM)**2, KKA, KKU, KKL) *PETHETA
    ELSE
      ZW1 = MZM((CST%XG / PTHVREF * PLM * PLEPS / PTKEM)**2, KKA, KKU, KKL)
    END IF
    PRED2THS3(:,:,:,JSV) = PREDTH1(:,:,:) * PREDS1(:,:,:,JSV)   +        &
                       ZW1*                                              &
                       MZM(GX_M_M(PSVM(:,:,:,JSV),PDXX,PDZZ,PDZX, KKA, KKU, KKL)*       &
                           GX_M_M(PTHLM,PDXX,PDZZ,PDZX, KKA, KKU, KKL)                  &
                          +GY_M_M(PSVM(:,:,:,JSV),PDYY,PDZZ,PDZY, KKA, KKU, KKL)*       &
                           GY_M_M(PTHLM,PDYY,PDZZ,PDZY, KKA, KKU, KKL),                 &
                           KKA, KKU, KKL)
!
    IF (KRR /= 0) THEN
      PRED2RS3(:,:,:,JSV) = PREDR1(:,:,:) * PREDS1(:,:,:,JSV)   +        &
                       ZW1 * PEMOIST *                                   &
                       MZM(GX_M_M(PSVM(:,:,:,JSV),PDXX,PDZZ,PDZX, KKA, KKU, KKL)*       &
                           GX_M_M(PRM(:,:,:,1),PDXX,PDZZ,PDZX, KKA, KKU, KKL)           &
                          +GY_M_M(PSVM(:,:,:,JSV),PDYY,PDZZ,PDZY, KKA, KKU, KKL)*       &
                           GY_M_M(PRM(:,:,:,1),PDYY,PDZZ,PDZY, KKA, KKU, KKL),          &
                           KKA, KKU, KKL)
    ELSE
      PRED2RS3(:,:,:,JSV) = 0.
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
  TZFIELD%CMNHNAME   = 'RED_TH1'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'RED_TH1'
  TZFIELD%CUNITS     = '1'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Z_RED_TH1'
  TZFIELD%NGRID      = 4
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,PREDTH1)
  !
  ! stores the RED_R1
  TZFIELD%CMNHNAME   = 'RED_R1'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'RED_R1'
  TZFIELD%CUNITS     = '1'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Z_RED_R1'
  TZFIELD%NGRID      = 4
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,PREDR1)
  !
  ! stores the RED2_TH3
  TZFIELD%CMNHNAME   = 'RED2_TH3'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'RED2_TH3'
  TZFIELD%CUNITS     = '1'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Z_RED2_TH3'
  TZFIELD%NGRID      = 4
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,PRED2TH3)
  !
  ! stores the RED2_R3
  TZFIELD%CMNHNAME   = 'RED2_R3'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'RED2_R3'
  TZFIELD%CUNITS     = '1'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Z_RED2_R3'
  TZFIELD%NGRID      = 4
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,PRED2R3)
  !
  ! stores the RED2_THR3
  TZFIELD%CMNHNAME   = 'RED2_THR3'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'RED2_THR3'
  TZFIELD%CUNITS     = '1'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Z_RED2_THR3'
  TZFIELD%NGRID      = 4
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,PRED2THR3)
  !
END IF
!
!---------------------------------------------------------------------------
ENDIF ! (Done only if OHARAT is FALSE)
!
IF (LHOOK) CALL DR_HOOK('PRANDTL',1,ZHOOK_HANDLE)
END SUBROUTINE PRANDTL
!
SUBROUTINE SMOOTH_TURB_FUNCT(CSTURB,PPHI3,PF_LIM,PF)
!
TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PPHI3   ! Phi3
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PF_LIM  ! Value of F when Phi3 is
!                                                ! larger than Phi_lim
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PF      ! function F to smooth
!
REAL, DIMENSION(SIZE(PF,1),SIZE(PF,2),SIZE(PF,3)) :: ZCOEF
!
!* adds a artificial correction to smooth the function near the discontinuity
!  point at Phi3 = Phi_lim
!  This smoothing is applied between 0.9*phi_lim (=2.7) and Phi_lim (=3)
!   Note that in the Boundary layer, phi is usually between 0.8 and 1
!
!
ZCOEF = MAX(MIN((  10.*(1.-PPHI3/CSTURB%XPHI_LIM)) ,1.), 0.) 
!
PF(:,:,:) =     ZCOEF(:,:,:)   * PF    &
          + (1.-ZCOEF(:,:,:))  * PF_LIM
!
END SUBROUTINE SMOOTH_TURB_FUNCT
!----------------------------------------------------------------------------
FUNCTION PHI3(CSTURB,PREDTH1,PREDR1,PRED2TH3,PRED2R3,PRED2THR3,HTURBDIM,OUSERV)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDR1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PRED2TH3
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PRED2R3
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PRED2THR3
  CHARACTER(len=4),       INTENT(IN) :: HTURBDIM  ! 1DIM or 3DIM turb. scheme
  LOGICAL,                INTENT(IN) :: OUSERV    ! flag to use vapor
  REAL, DIMENSION(SIZE(PREDTH1,1),SIZE(PREDTH1,2),SIZE(PREDTH1,3)) :: PHI3
!
  REAL, DIMENSION(SIZE(PREDTH1,1),SIZE(PREDTH1,2),SIZE(PREDTH1,3)) :: ZW1, ZW2
  INTEGER :: IKB, IKE
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:PHI3',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PREDTH1,3)-JPVEXT_TURB
!
IF (HTURBDIM=='3DIM') THEN
        !* 3DIM case
  IF (OUSERV) THEN
    ZW1(:,:,:) = 1. + 1.5* (PREDTH1(:,:,:)+PREDR1(:,:,:)) +      &
                   ( 0.5 * (PREDTH1(:,:,:)**2+PREDR1(:,:,:)**2)  &
                         + PREDTH1(:,:,:) * PREDR1(:,:,:)        &
                   )

    ZW2(:,:,:) = 0.5 * (PRED2TH3(:,:,:)-PRED2R3(:,:,:))

    PHI3(:,:,:)= 1. -                                          &
    ( ( (1.+PREDR1(:,:,:)) *                                   &
        (PRED2THR3(:,:,:) + PRED2TH3(:,:,:)) / PREDTH1(:,:,:)  &
      ) + ZW2(:,:,:)                                           &
    ) / ZW1(:,:,:)
  ELSE
    ZW1(:,:,:) = 1. + 1.5* PREDTH1(:,:,:) + &
                 0.5* PREDTH1(:,:,:)**2

    ZW2(:,:,:) = 0.5* PRED2TH3(:,:,:)

    PHI3(:,:,:)= 1. -                                       &
            (PRED2TH3(:,:,:) / PREDTH1(:,:,:) + ZW2(:,:,:)) / ZW1(:,:,:)
  END IF
  WHERE( PHI3 <= 0. .OR. PHI3 > CSTURB%XPHI_LIM )
    PHI3 = CSTURB%XPHI_LIM
  END WHERE

ELSE
        !* 1DIM case
  IF (OUSERV) THEN
    PHI3(:,:,:)= 1./(1.+PREDTH1(:,:,:)+PREDR1(:,:,:))
  ELSE
    PHI3(:,:,:)= 1./(1.+PREDTH1(:,:,:))
  END IF
END IF
!
PHI3(:,:,IKB-1)=PHI3(:,:,IKB)
PHI3(:,:,IKE+1)=PHI3(:,:,IKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:PHI3',1,ZHOOK_HANDLE)
END FUNCTION PHI3
!----------------------------------------------------------------------------
FUNCTION PSI_SV(CSTURB,PREDTH1,PREDR1,PREDS1,PRED2THS,PRED2RS,PPHI3,PPSI3)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  REAL, DIMENSION(:,:,:),   INTENT(IN) :: PREDTH1
  REAL, DIMENSION(:,:,:),   INTENT(IN) :: PREDR1
  REAL, DIMENSION(:,:,:,:), INTENT(IN) :: PREDS1
  REAL, DIMENSION(:,:,:,:), INTENT(IN) :: PRED2THS
  REAL, DIMENSION(:,:,:,:), INTENT(IN) :: PRED2RS
  REAL, DIMENSION(:,:,:),   INTENT(IN) :: PPHI3
  REAL, DIMENSION(:,:,:),   INTENT(IN) :: PPSI3
  REAL, DIMENSION(SIZE(PRED2THS,1),SIZE(PRED2THS,2),SIZE(PRED2THS,3),SIZE(PRED2THS,4)) :: PSI_SV
!
  INTEGER :: IKB, IKE
  INTEGER :: JSV
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:PSI_SV',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PREDTH1,3)-JPVEXT_TURB
!
DO JSV=1,SIZE(PSI_SV,4)
  PSI_SV(:,:,:,JSV) = ( 1.                                             &
    - (CSTURB%XCPR3+CSTURB%XCPR5) * (PRED2THS(:,:,:,JSV)/PREDS1(:,:,:,JSV)-PREDTH1) &
    - (CSTURB%XCPR4+CSTURB%XCPR5) * (PRED2RS (:,:,:,JSV)/PREDS1(:,:,:,JSV)-PREDR1 ) &
    - CSTURB%XCPR3 * PREDTH1 * PPHI3 - CSTURB%XCPR4 * PREDR1 * PPSI3                 &
                ) / (  1. + CSTURB%XCPR5 * ( PREDTH1 + PREDR1 ) )           
  
!        control of the PSI_SV positivity
  WHERE ( (PSI_SV(:,:,:,JSV) <=0.).AND. (PREDTH1+PREDR1) <= 0. )
    PSI_SV(:,:,:,JSV)=CSTURB%XPHI_LIM
  END WHERE
  PSI_SV(:,:,:,JSV) = MAX( 1.E-4, MIN(CSTURB%XPHI_LIM,PSI_SV(:,:,:,JSV)) )
!
  PSI_SV(:,:,IKB-1,JSV)=PSI_SV(:,:,IKB,JSV)
  PSI_SV(:,:,IKE+1,JSV)=PSI_SV(:,:,IKE,JSV)
END DO
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:PSI_SV',1,ZHOOK_HANDLE)
END FUNCTION PSI_SV
!----------------------------------------------------------------------------
FUNCTION D_PHI3DTDZ_O_DDTDZ(CSTURB,PPHI3,PREDTH1,PREDR1,PRED2TH3,PRED2THR3,HTURBDIM,OUSERV)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PPHI3
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDR1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PRED2TH3
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PRED2THR3
  CHARACTER(len=4),       INTENT(IN) :: HTURBDIM  ! 1DIM or 3DIM turb. scheme
  LOGICAL,                INTENT(IN) :: OUSERV    ! flag to use vapor
  REAL, DIMENSION(SIZE(PREDTH1,1),SIZE(PREDTH1,2),SIZE(PREDTH1,3)) :: D_PHI3DTDZ_O_DDTDZ
  INTEGER :: IKB, IKE,JL,JK,JJ
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_PHI3DTDZ_O_DDTDZ',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PREDTH1,3)-JPVEXT_TURB
!
IF (HTURBDIM=='3DIM') THEN
        !* 3DIM case
  IF (OUSERV) THEN
#ifdef REPRO48
    WHERE (PPHI3(:,:,:)/=CSTURB%XPHI_LIM)
#else
    WHERE (PPHI3(:,:,:)<=CSTURB%XPHI_LIM)
#endif
    D_PHI3DTDZ_O_DDTDZ(:,:,:) = PPHI3(:,:,:)                       &
          * (1. - PREDTH1(:,:,:) * (3./2.+PREDTH1+PREDR1)          &
               /((1.+PREDTH1+PREDR1)*(1.+1./2.*(PREDTH1+PREDR1)))) &
          + (1.+PREDR1)*(PRED2THR3+PRED2TH3)                       &
               / (PREDTH1*(1.+PREDTH1+PREDR1)*(1.+1./2.*(PREDTH1+PREDR1))) &
          - (1./2.*PREDTH1+PREDR1 * (1.+PREDTH1+PREDR1))           &
               / ((1.+PREDTH1+PREDR1)*(1.+1./2.*(PREDTH1+PREDR1)))
    ELSEWHERE
      D_PHI3DTDZ_O_DDTDZ(:,:,:) = PPHI3(:,:,:)
    ENDWHERE

!
  ELSE
#ifdef REPRO48
    WHERE (PPHI3(:,:,:)/=CSTURB%XPHI_LIM)
#else
    WHERE (PPHI3(:,:,:)<=CSTURB%XPHI_LIM)
#endif
    D_PHI3DTDZ_O_DDTDZ(:,:,:) = PPHI3(:,:,:)             &
          * (1. - PREDTH1(:,:,:) * (3./2.+PREDTH1)      &
               /((1.+PREDTH1)*(1.+1./2.*PREDTH1)))        &
          + PRED2TH3 / (PREDTH1*(1.+PREDTH1)*(1.+1./2.*PREDTH1)) &
          - 1./2.*PREDTH1 / ((1.+PREDTH1)*(1.+1./2.*PREDTH1))
    ELSEWHERE
      D_PHI3DTDZ_O_DDTDZ(:,:,:) = PPHI3(:,:,:)
    ENDWHERE
!
  END IF
ELSE
        !* 1DIM case
DO JJ=1,SIZE(PPHI3,2)
  DO JL=1,SIZE(PPHI3,1)
    DO JK=1,SIZE(PPHI3,3)
      IF ( ABS(PPHI3(JL,JJ,JK)-CSTURB%XPHI_LIM) < 1.E-12 ) THEN
         D_PHI3DTDZ_O_DDTDZ(JL,JJ,JK)=PPHI3(JL,JJ,JK)*&
&       (1. - PREDTH1(JL,JJ,JK)*PPHI3(JL,JJ,JK))
      ELSE
         D_PHI3DTDZ_O_DDTDZ(JL,JJ,JK)=PPHI3(JL,JJ,JK)
      ENDIF
    ENDDO
  ENDDO
ENDDO
END IF
!
#ifdef REPRO48
#else
!* smoothing
CALL SMOOTH_TURB_FUNCT(CSTURB,PPHI3,PPHI3,D_PHI3DTDZ_O_DDTDZ)
#endif
!
D_PHI3DTDZ_O_DDTDZ(:,:,IKB-1)=D_PHI3DTDZ_O_DDTDZ(:,:,IKB)
D_PHI3DTDZ_O_DDTDZ(:,:,IKE+1)=D_PHI3DTDZ_O_DDTDZ(:,:,IKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_PHI3DTDZ_O_DDTDZ',1,ZHOOK_HANDLE)
END FUNCTION D_PHI3DTDZ_O_DDTDZ
!----------------------------------------------------------------------------
FUNCTION D_PHI3DRDZ_O_DDRDZ(CSTURB,PPHI3,PREDTH1,PREDR1,PRED2TH3,PRED2THR3,HTURBDIM,OUSERV)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PPHI3
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDR1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PRED2TH3
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PRED2THR3
  CHARACTER(len=4),       INTENT(IN) :: HTURBDIM  ! 1DIM or 3DIM turb. scheme
  LOGICAL,                INTENT(IN) :: OUSERV    ! flag to use vapor
   REAL, DIMENSION(SIZE(PREDTH1,1),SIZE(PREDTH1,2),SIZE(PREDTH1,3)) :: D_PHI3DRDZ_O_DDRDZ
  INTEGER :: IKB, IKE
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_PHI3DRDZ_O_DDRDZ',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PREDTH1,3)-JPVEXT_TURB
!
!
IF (HTURBDIM=='3DIM') THEN
        !* 3DIM case
  IF (OUSERV) THEN
#ifdef REPRO48
    WHERE (PPHI3(:,:,:)/=CSTURB%XPHI_LIM)
#else
    WHERE (PPHI3(:,:,:)<=CSTURB%XPHI_LIM)
#endif
      D_PHI3DRDZ_O_DDRDZ(:,:,:) =          &
                PPHI3(:,:,:) * (1.-PREDR1(:,:,:)*(3./2.+PREDTH1+PREDR1) &
                  / ((1.+PREDTH1+PREDR1)*(1.+1./2.*(PREDTH1+PREDR1))))  &
              - PREDR1(:,:,:) * (PRED2THR3+PRED2TH3) / (PREDTH1         &
                  * (1.+PREDTH1+PREDR1)*(1.+1./2.*(PREDTH1+PREDR1)))    &
              + PREDR1(:,:,:) * (1./2.+PREDTH1+PREDR1)                  &
                  / ((1.+PREDTH1+PREDR1)*(1.+1./2.*(PREDTH1+PREDR1)))
    ELSEWHERE
      D_PHI3DRDZ_O_DDRDZ(:,:,:) = PPHI3(:,:,:)
    END WHERE
  ELSE
    D_PHI3DRDZ_O_DDRDZ(:,:,:) = PPHI3(:,:,:)
  END IF
ELSE
        !* 1DIM case
#ifdef REPRO48
    WHERE (PPHI3(:,:,:)/=CSTURB%XPHI_LIM)
#else
    WHERE (PPHI3(:,:,:)<=CSTURB%XPHI_LIM)
#endif
    D_PHI3DRDZ_O_DDRDZ(:,:,:) = PPHI3(:,:,:)                           &
          * (1. - PREDR1(:,:,:)*PPHI3(:,:,:))
  ELSEWHERE
    D_PHI3DRDZ_O_DDRDZ(:,:,:) = PPHI3(:,:,:)
  END WHERE
END IF
!
#ifdef REPRO48
#else
!* smoothing
CALL SMOOTH_TURB_FUNCT(CSTURB,PPHI3,PPHI3,D_PHI3DRDZ_O_DDRDZ)
#endif
!
D_PHI3DRDZ_O_DDRDZ(:,:,IKB-1)=D_PHI3DRDZ_O_DDRDZ(:,:,IKB)
D_PHI3DRDZ_O_DDRDZ(:,:,IKE+1)=D_PHI3DRDZ_O_DDRDZ(:,:,IKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_PHI3DRDZ_O_DDRDZ',1,ZHOOK_HANDLE)
END FUNCTION D_PHI3DRDZ_O_DDRDZ
!----------------------------------------------------------------------------
FUNCTION D_PHI3DTDZ2_O_DDTDZ(CSTURB,PPHI3,PREDTH1,PREDR1,PRED2TH3,PRED2THR3,PDTDZ,HTURBDIM,OUSERV)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PPHI3
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDR1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PRED2TH3
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PRED2THR3
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PDTDZ
  CHARACTER(len=4),       INTENT(IN) :: HTURBDIM  ! 1DIM or 3DIM turb. scheme
  LOGICAL,                INTENT(IN) :: OUSERV    ! flag to use vapor
  REAL, DIMENSION(SIZE(PREDTH1,1),SIZE(PREDTH1,2),SIZE(PREDTH1,3)) :: D_PHI3DTDZ2_O_DDTDZ
  INTEGER :: IKB, IKE
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_PHI3DTDZ2_O_DDTDZ',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PREDTH1,3)-JPVEXT_TURB
!
!
IF (HTURBDIM=='3DIM') THEN
   ! by derivation of (phi3 dtdz) * dtdz according to dtdz we obtain:
   D_PHI3DTDZ2_O_DDTDZ(:,:,:) = PDTDZ * (PPHI3 +  &
           D_PHI3DTDZ_O_DDTDZ(CSTURB,PPHI3,PREDTH1,PREDR1,PRED2TH3,PRED2THR3,HTURBDIM,OUSERV) )

ELSE
        !* 1DIM case
#ifdef REPRO48
    WHERE (PPHI3(:,:,:)/=CSTURB%XPHI_LIM)
#else
    WHERE (PPHI3(:,:,:)<=CSTURB%XPHI_LIM)
#endif
      D_PHI3DTDZ2_O_DDTDZ(:,:,:) = PPHI3(:,:,:)*PDTDZ(:,:,:)             &
          * (2. - PREDTH1(:,:,:)*PPHI3(:,:,:))
    ELSEWHERE
      D_PHI3DTDZ2_O_DDTDZ(:,:,:) = PPHI3(:,:,:) * 2. * PDTDZ(:,:,:)
    END WHERE
END IF
!
#ifdef REPRO48
#else
!* smoothing
CALL SMOOTH_TURB_FUNCT(CSTURB,PPHI3,PPHI3*2.*PDTDZ,D_PHI3DTDZ2_O_DDTDZ)
#endif
!
!
D_PHI3DTDZ2_O_DDTDZ(:,:,IKB-1)=D_PHI3DTDZ2_O_DDTDZ(:,:,IKB)
D_PHI3DTDZ2_O_DDTDZ(:,:,IKE+1)=D_PHI3DTDZ2_O_DDTDZ(:,:,IKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_PHI3DTDZ2_O_DDTDZ',1,ZHOOK_HANDLE)
END FUNCTION D_PHI3DTDZ2_O_DDTDZ
!----------------------------------------------------------------------------
FUNCTION M3_WTH_WTH2(CSTURB,PREDTH1,PREDR1,PD,PBLL_O_E,PETHETA)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDR1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PD
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PETHETA
  REAL, DIMENSION(SIZE(PD,1),SIZE(PD,2),SIZE(PD,3)) :: M3_WTH_WTH2
  INTEGER :: IKB, IKE
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_WTH_WTH2',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PD,3)-JPVEXT_TURB

M3_WTH_WTH2(:,:,:) = CSTURB%XCSHF*PBLL_O_E*PETHETA*0.5/CSTURB%XCTD        &
                   * (1.+0.5*PREDTH1+PREDR1) / PD
M3_WTH_WTH2(:,:,IKB-1)=M3_WTH_WTH2(:,:,IKB)
M3_WTH_WTH2(:,:,IKE+1)=M3_WTH_WTH2(:,:,IKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_WTH_WTH2',1,ZHOOK_HANDLE)
END FUNCTION M3_WTH_WTH2
!----------------------------------------------------------------------------
FUNCTION D_M3_WTH_WTH2_O_DDTDZ(CSTURB,PM3_WTH_WTH2,PREDTH1,PREDR1,PD,PBLL_O_E,PETHETA)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PM3_WTH_WTH2
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDR1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PD
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PETHETA
  REAL, DIMENSION(SIZE(PD,1),SIZE(PD,2),SIZE(PD,3)) :: D_M3_WTH_WTH2_O_DDTDZ
  INTEGER :: IKB, IKE
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_WTH_WTH2_O_DDTDZ',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PD,3)-JPVEXT_TURB

D_M3_WTH_WTH2_O_DDTDZ(:,:,:) = (  0.5*CSTURB%XCSHF*PBLL_O_E*PETHETA*0.5/CSTURB%XCTD/PD &
                                - PM3_WTH_WTH2/PD*(1.5+PREDTH1+PREDR1)  )&
                             * PBLL_O_E * PETHETA * CSTURB%XCTV
!
D_M3_WTH_WTH2_O_DDTDZ(:,:,IKB-1)=D_M3_WTH_WTH2_O_DDTDZ(:,:,IKB)
D_M3_WTH_WTH2_O_DDTDZ(:,:,IKE+1)=D_M3_WTH_WTH2_O_DDTDZ(:,:,IKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_WTH_WTH2_O_DDTDZ',1,ZHOOK_HANDLE)
END FUNCTION D_M3_WTH_WTH2_O_DDTDZ
!----------------------------------------------------------------------------
FUNCTION M3_WTH_W2TH(CSTURB,KKA,KKU,KKL,PREDTH1,PREDR1,PD,PKEFF,PTKE)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  INTEGER,                INTENT(IN) :: KKA 
  INTEGER,                INTENT(IN) :: KKU  
  INTEGER,                INTENT(IN) :: KKL
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDR1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PD
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PKEFF
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PTKE
  REAL, DIMENSION(SIZE(PD,1),SIZE(PD,2),SIZE(PD,3)) :: M3_WTH_W2TH
  INTEGER :: IKB, IKE
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_WTH_W2TH',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PD,3)-JPVEXT_TURB

M3_WTH_W2TH(:,:,:) = CSTURB%XCSHF*PKEFF*1.5/MZM(PTKE, KKA, KKU, KKL)              &
  * (1. - 0.5*PREDR1*(1.+PREDR1)/PD ) / (1.+PREDTH1)
!
M3_WTH_W2TH(:,:,IKB-1)=M3_WTH_W2TH(:,:,IKB)
M3_WTH_W2TH(:,:,IKE+1)=M3_WTH_W2TH(:,:,IKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_WTH_W2TH',1,ZHOOK_HANDLE)
END FUNCTION M3_WTH_W2TH
!----------------------------------------------------------------------------
FUNCTION D_M3_WTH_W2TH_O_DDTDZ(CSTURB,KKA,KKU,KKL,PREDTH1,PREDR1,PD,PBLL_O_E,PETHETA,PKEFF,PTKE)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  INTEGER,                INTENT(IN) :: KKA 
  INTEGER,                INTENT(IN) :: KKU  
  INTEGER,                INTENT(IN) :: KKL
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDR1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PD
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PETHETA
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PKEFF
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PTKE
  REAL, DIMENSION(SIZE(PD,1),SIZE(PD,2),SIZE(PD,3)) :: D_M3_WTH_W2TH_O_DDTDZ
  INTEGER :: IKB, IKE
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_WTH_W2TH_O_DDTDZ',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PD,3)-JPVEXT_TURB

D_M3_WTH_W2TH_O_DDTDZ(:,:,:) = &
 - CSTURB%XCSHF*PKEFF*1.5/MZM(PTKE, KKA, KKU, KKL)/(1.+PREDTH1)**2*CSTURB%XCTV*PBLL_O_E*PETHETA  &
 * (1. - 0.5*PREDR1*(1.+PREDR1)/PD*( 1.+(1.+PREDTH1)*(1.5+PREDR1+PREDTH1)/PD) )
!
D_M3_WTH_W2TH_O_DDTDZ(:,:,IKB-1)=D_M3_WTH_W2TH_O_DDTDZ(:,:,IKB)
D_M3_WTH_W2TH_O_DDTDZ(:,:,IKE+1)=D_M3_WTH_W2TH_O_DDTDZ(:,:,IKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_WTH_W2TH_O_DDTDZ',1,ZHOOK_HANDLE)
END FUNCTION D_M3_WTH_W2TH_O_DDTDZ
!----------------------------------------------------------------------------
FUNCTION M3_WTH_W2R(CSTURB,KKA,KKU,KKL,PD,PKEFF,PTKE,PBLL_O_E,PEMOIST,PDTDZ)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  INTEGER,                INTENT(IN) :: KKA 
  INTEGER,                INTENT(IN) :: KKU  
  INTEGER,                INTENT(IN) :: KKL
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PD
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PKEFF
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PTKE
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PEMOIST
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PDTDZ
  REAL, DIMENSION(SIZE(PD,1),SIZE(PD,2),SIZE(PD,3)) :: M3_WTH_W2R
  INTEGER :: IKB, IKE
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_WTH_W2R',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PD,3)-JPVEXT_TURB

M3_WTH_W2R(:,:,:) = - CSTURB%XCSHF*PKEFF*0.75*CSTURB%XCTV*PBLL_O_E/MZM(PTKE, KKA, KKU, KKL)*PEMOIST*PDTDZ/PD
!
M3_WTH_W2R(:,:,IKB-1)=M3_WTH_W2R(:,:,IKB)
M3_WTH_W2R(:,:,IKE+1)=M3_WTH_W2R(:,:,IKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_WTH_W2R',1,ZHOOK_HANDLE)
END FUNCTION M3_WTH_W2R
!----------------------------------------------------------------------------
FUNCTION D_M3_WTH_W2R_O_DDTDZ(CSTURB,KKA,KKU,KKL,PREDTH1,PREDR1,PD,PKEFF,PTKE,PBLL_O_E,PEMOIST)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  INTEGER,                INTENT(IN) :: KKA 
  INTEGER,                INTENT(IN) :: KKU  
  INTEGER,                INTENT(IN) :: KKL
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDR1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PD
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PKEFF
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PTKE
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PEMOIST
  REAL, DIMENSION(SIZE(PD,1),SIZE(PD,2),SIZE(PD,3)) :: D_M3_WTH_W2R_O_DDTDZ
  INTEGER :: IKB, IKE
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_WTH_W2R_O_DDTDZ',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PD,3)-JPVEXT_TURB

D_M3_WTH_W2R_O_DDTDZ(:,:,:) = - CSTURB%XCSHF*PKEFF*0.75*CSTURB%XCTV*PBLL_O_E/MZM(PTKE, KKA, KKU, KKL)*PEMOIST/PD &
                                     * (1. -  PREDTH1*(1.5+PREDTH1+PREDR1)/PD)
!
D_M3_WTH_W2R_O_DDTDZ(:,:,IKB-1)=D_M3_WTH_W2R_O_DDTDZ(:,:,IKB)
D_M3_WTH_W2R_O_DDTDZ(:,:,IKE+1)=D_M3_WTH_W2R_O_DDTDZ(:,:,IKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_WTH_W2R_O_DDTDZ',1,ZHOOK_HANDLE)
END FUNCTION D_M3_WTH_W2R_O_DDTDZ
!----------------------------------------------------------------------------
FUNCTION M3_WTH_WR2(CSTURB,KKA,KKU,KKL,PD,PKEFF,PTKE,PSQRT_TKE,PBLL_O_E,PBETA,PLEPS,PEMOIST,PDTDZ)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  INTEGER,                INTENT(IN) :: KKA 
  INTEGER,                INTENT(IN) :: KKU  
  INTEGER,                INTENT(IN) :: KKL
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PD
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PKEFF
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PTKE
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PBETA
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PLEPS
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PEMOIST
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PDTDZ
  REAL, DIMENSION(SIZE(PD,1),SIZE(PD,2),SIZE(PD,3)) :: M3_WTH_WR2
  INTEGER :: IKB, IKE
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_WTH_WR2',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PD,3)-JPVEXT_TURB

M3_WTH_WR2(:,:,:) = - CSTURB%XCSHF*PKEFF*0.25*PBLL_O_E*CSTURB%XCTV*PEMOIST**2       &
                           *MZM(PBETA*PLEPS/(PSQRT_TKE*PTKE), KKA, KKU, KKL)/CSTURB%XCTD*PDTDZ/PD
!
M3_WTH_WR2(:,:,IKB-1)=M3_WTH_WR2(:,:,IKB)
M3_WTH_WR2(:,:,IKE+1)=M3_WTH_WR2(:,:,IKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_WTH_WR2',1,ZHOOK_HANDLE)
END FUNCTION M3_WTH_WR2
!----------------------------------------------------------------------------
FUNCTION D_M3_WTH_WR2_O_DDTDZ(CSTURB,KKA,KKU,KKL,PREDTH1,PREDR1,PD,PKEFF,PTKE,PSQRT_TKE,PBLL_O_E,PBETA,PLEPS,PEMOIST)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  INTEGER,                INTENT(IN) :: KKA 
  INTEGER,                INTENT(IN) :: KKU  
  INTEGER,                INTENT(IN) :: KKL
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDR1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PD
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PKEFF
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PTKE
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PBETA
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PLEPS
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PEMOIST
  REAL, DIMENSION(SIZE(PD,1),SIZE(PD,2),SIZE(PD,3)) :: D_M3_WTH_WR2_O_DDTDZ
  INTEGER :: IKB, IKE
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_WTH_WR2_O_DDTDZ',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PD,3)-JPVEXT_TURB

D_M3_WTH_WR2_O_DDTDZ(:,:,:) = - CSTURB%XCSHF*PKEFF*0.25*PBLL_O_E*CSTURB%XCTV*PEMOIST**2 &
                           *MZM(PBETA*PLEPS/(PSQRT_TKE*PTKE), KKA, KKU, KKL)/CSTURB%XCTD/PD     &
                           * (1. -  PREDTH1*(1.5+PREDTH1+PREDR1)/PD)
!
D_M3_WTH_WR2_O_DDTDZ(:,:,IKB-1)=D_M3_WTH_WR2_O_DDTDZ(:,:,IKB)
D_M3_WTH_WR2_O_DDTDZ(:,:,IKE+1)=D_M3_WTH_WR2_O_DDTDZ(:,:,IKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_WTH_WR2_O_DDTDZ',1,ZHOOK_HANDLE)
END FUNCTION D_M3_WTH_WR2_O_DDTDZ
!----------------------------------------------------------------------------
FUNCTION M3_WTH_WTHR(CSTURB,KKA,KKU,KKL,PREDR1,PD,PKEFF,PTKE,PSQRT_TKE,PBETA,PLEPS,PEMOIST)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  INTEGER,                INTENT(IN) :: KKA 
  INTEGER,                INTENT(IN) :: KKU  
  INTEGER,                INTENT(IN) :: KKL
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDR1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PD
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PKEFF
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PTKE
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PBETA
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PLEPS
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PEMOIST
  REAL, DIMENSION(SIZE(PREDR1,1),SIZE(PREDR1,2),SIZE(PREDR1,3)) :: M3_WTH_WTHR
  INTEGER :: IKB, IKE
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_WTH_WTHR',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PD,3)-JPVEXT_TURB

!M3_WTH_WTHR(:,:,:) = CSTURB%XCSHF*PKEFF*PEMOIST/MZM(PBETA*PTKE*PSQRT_TKE, KKA, KKU, KKL) &
!                         *0.5*PLEPS/CSTURB%XCTD*(1+PREDR1)/PD
M3_WTH_WTHR(:,:,:) = CSTURB%XCSHF*PKEFF*PEMOIST*MZM(PBETA/PTKE*PSQRT_TKE, KKA, KKU, KKL) &
                         *0.5*PLEPS/CSTURB%XCTD*(1+PREDR1)/PD
!
M3_WTH_WTHR(:,:,IKB-1)=M3_WTH_WTHR(:,:,IKB)
M3_WTH_WTHR(:,:,IKE+1)=M3_WTH_WTHR(:,:,IKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_WTH_WTHR',1,ZHOOK_HANDLE)
END FUNCTION M3_WTH_WTHR
!----------------------------------------------------------------------------
FUNCTION D_M3_WTH_WTHR_O_DDTDZ(CSTURB,PM3_WTH_WTHR,PREDTH1,PREDR1,PD,PBLL_O_E,PETHETA)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PM3_WTH_WTHR
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDR1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PD
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PETHETA
  REAL, DIMENSION(SIZE(PD,1),SIZE(PD,2),SIZE(PD,3)) :: D_M3_WTH_WTHR_O_DDTDZ
  INTEGER :: IKB, IKE
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_WTH_WTHR_O_DDTDZ',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PD,3)-JPVEXT_TURB

D_M3_WTH_WTHR_O_DDTDZ(:,:,:) = - PM3_WTH_WTHR * (1.5+PREDTH1+PREDR1)/PD*CSTURB%XCTV*PBLL_O_E*PETHETA
!
D_M3_WTH_WTHR_O_DDTDZ(:,:,IKB-1)=D_M3_WTH_WTHR_O_DDTDZ(:,:,IKB)
D_M3_WTH_WTHR_O_DDTDZ(:,:,IKE+1)=D_M3_WTH_WTHR_O_DDTDZ(:,:,IKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_WTH_WTHR_O_DDTDZ',1,ZHOOK_HANDLE)
END FUNCTION D_M3_WTH_WTHR_O_DDTDZ
!----------------------------------------------------------------------------
FUNCTION M3_TH2_W2TH(CSTURB,KKA,KKU,KKL,PREDTH1,PREDR1,PD,PDTDZ,PLM,PLEPS,PTKE)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  INTEGER,                INTENT(IN) :: KKA 
  INTEGER,                INTENT(IN) :: KKU  
  INTEGER,                INTENT(IN) :: KKL
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDR1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PD
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PDTDZ
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PLM
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PLEPS
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PTKE
  REAL, DIMENSION(SIZE(PD,1),SIZE(PD,2),SIZE(PD,3)) :: M3_TH2_W2TH
  INTEGER :: IKB, IKE
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_TH2_W2TH',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PD,3)-JPVEXT_TURB

M3_TH2_W2TH(:,:,:) = - MZF((1.-0.5*PREDR1*(1.+PREDR1)/PD)/(1.+PREDTH1)*PDTDZ, KKA, KKU, KKL) &
                       * 1.5*PLM*PLEPS/PTKE*CSTURB%XCTV
!
M3_TH2_W2TH(:,:,IKB-1)=M3_TH2_W2TH(:,:,IKB)
M3_TH2_W2TH(:,:,IKE+1)=M3_TH2_W2TH(:,:,IKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_TH2_W2TH',1,ZHOOK_HANDLE)
END FUNCTION M3_TH2_W2TH
!----------------------------------------------------------------------------
FUNCTION D_M3_TH2_W2TH_O_DDTDZ(CSTURB,KKA,KKU,KKL,PREDTH1,PREDR1,PD,PLM,PLEPS,PTKE,OUSERV)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  INTEGER,                INTENT(IN) :: KKA 
  INTEGER,                INTENT(IN) :: KKU  
  INTEGER,                INTENT(IN) :: KKL
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDR1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PD
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PLM
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PLEPS
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PTKE
  LOGICAL,                INTENT(IN) :: OUSERV
  REAL, DIMENSION(SIZE(PD,1),SIZE(PD,2),SIZE(PD,3)) :: D_M3_TH2_W2TH_O_DDTDZ
  INTEGER :: IKB, IKE
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_TH2_W2TH_O_DDTDZ',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PD,3)-JPVEXT_TURB

IF (OUSERV) THEN
!  D_M3_TH2_W2TH_O_DDTDZ(:,:,:) = - 1.5*PLM*PLEPS/PTKE*CSTURB%XCTV * MZF(                    &
!          (1.-0.5*PREDR1*(1.+PREDR1)/PD)*(1.-(1.5+PREDTH1+PREDR1)*(1.+PREDTH1)/PD )  &
!        / (1.+PREDTH1)**2, KKA, KKU, KKL)
  D_M3_TH2_W2TH_O_DDTDZ(:,:,:) = - 1.5*PLM*PLEPS/PTKE*CSTURB%XCTV * MZF( &
          (1.-0.5*PREDR1*(1.+PREDR1)/PD)*(1.-(1.5+PREDTH1+PREDR1)*   &
             PREDTH1*(1.+PREDTH1)/PD ) / (1.+PREDTH1)**2, KKA, KKU, KKL)

ELSE
  D_M3_TH2_W2TH_O_DDTDZ(:,:,:) = - 1.5*PLM*PLEPS/PTKE*CSTURB%XCTV * MZF(1./(1.+PREDTH1)**2, KKA, KKU, KKL)
END IF
!
D_M3_TH2_W2TH_O_DDTDZ(:,:,IKB-1)=D_M3_TH2_W2TH_O_DDTDZ(:,:,IKB)
D_M3_TH2_W2TH_O_DDTDZ(:,:,IKE+1)=D_M3_TH2_W2TH_O_DDTDZ(:,:,IKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_TH2_W2TH_O_DDTDZ',1,ZHOOK_HANDLE)
END FUNCTION D_M3_TH2_W2TH_O_DDTDZ
!----------------------------------------------------------------------------
FUNCTION M3_TH2_WTH2(CSTURB,KKA,KKU,KKL,PREDTH1,PREDR1,PD,PLEPS,PSQRT_TKE)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  INTEGER,                INTENT(IN) :: KKA 
  INTEGER,                INTENT(IN) :: KKU  
  INTEGER,                INTENT(IN) :: KKL
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDR1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PD
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PLEPS
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(SIZE(PD,1),SIZE(PD,2),SIZE(PD,3)) :: M3_TH2_WTH2
  INTEGER :: IKB, IKE
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_TH2_WTH2',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PD,3)-JPVEXT_TURB

M3_TH2_WTH2(:,:,:) = PLEPS*0.5/CSTURB%XCTD/PSQRT_TKE          &
  * MZF((1.+0.5*PREDTH1+1.5*PREDR1+0.5*PREDR1**2)/PD, KKA, KKU, KKL)
!
M3_TH2_WTH2(:,:,IKB-1)=M3_TH2_WTH2(:,:,IKB)
M3_TH2_WTH2(:,:,IKE+1)=M3_TH2_WTH2(:,:,IKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_TH2_WTH2',1,ZHOOK_HANDLE)
END FUNCTION M3_TH2_WTH2
!----------------------------------------------------------------------------
FUNCTION D_M3_TH2_WTH2_O_DDTDZ(CSTURB,KKA,KKU,KKL,PREDTH1,PREDR1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PETHETA)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  INTEGER,                INTENT(IN) :: KKA 
  INTEGER,                INTENT(IN) :: KKU  
  INTEGER,                INTENT(IN) :: KKL
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDR1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PD
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PLEPS
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PETHETA
  REAL, DIMENSION(SIZE(PD,1),SIZE(PD,2),SIZE(PD,3)) :: D_M3_TH2_WTH2_O_DDTDZ
  INTEGER :: IKB, IKE
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_TH2_WTH2_O_DDTDZ',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PD,3)-JPVEXT_TURB

D_M3_TH2_WTH2_O_DDTDZ(:,:,:) = PLEPS*0.5/CSTURB%XCTD/PSQRT_TKE*CSTURB%XCTV                        &
 * MZF(PBLL_O_E*PETHETA* (0.5/PD                                                   &
             - (1.5+PREDTH1+PREDR1)*(1.+0.5*PREDTH1+1.5*PREDR1+0.5*PREDR1**2)/PD**2 &
                           ), KKA, KKU, KKL)
!
D_M3_TH2_WTH2_O_DDTDZ(:,:,IKB-1)=D_M3_TH2_WTH2_O_DDTDZ(:,:,IKB)
D_M3_TH2_WTH2_O_DDTDZ(:,:,IKE+1)=D_M3_TH2_WTH2_O_DDTDZ(:,:,IKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_TH2_WTH2_O_DDTDZ',1,ZHOOK_HANDLE)
END FUNCTION D_M3_TH2_WTH2_O_DDTDZ
!----------------------------------------------------------------------------
FUNCTION M3_TH2_W2R(CSTURB,KKA,KKU,KKL,PD,PLM,PLEPS,PTKE,PBLL_O_E,PEMOIST,PDTDZ)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  INTEGER,                INTENT(IN) :: KKA 
  INTEGER,                INTENT(IN) :: KKU  
  INTEGER,                INTENT(IN) :: KKL
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PD
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PLM
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PLEPS
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PTKE
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PEMOIST
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PDTDZ
  REAL, DIMENSION(SIZE(PD,1),SIZE(PD,2),SIZE(PD,3)) :: M3_TH2_W2R
  INTEGER :: IKB, IKE
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_TH2_W2R',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PD,3)-JPVEXT_TURB

M3_TH2_W2R(:,:,:) = 0.75*CSTURB%XCTV**2*MZF(PBLL_O_E*PEMOIST/PD*PDTDZ**2, KKA, KKU, KKL)*PLM*PLEPS/PTKE
!
M3_TH2_W2R(:,:,IKB-1)=M3_TH2_W2R(:,:,IKB)
M3_TH2_W2R(:,:,IKE+1)=M3_TH2_W2R(:,:,IKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_TH2_W2R',1,ZHOOK_HANDLE)
END FUNCTION M3_TH2_W2R
!----------------------------------------------------------------------------
FUNCTION D_M3_TH2_W2R_O_DDTDZ(CSTURB,KKA,KKU,KKL,PREDTH1,PREDR1,PD,PLM,PLEPS,PTKE,PBLL_O_E,PEMOIST,PDTDZ)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  INTEGER,                INTENT(IN) :: KKA 
  INTEGER,                INTENT(IN) :: KKU  
  INTEGER,                INTENT(IN) :: KKL
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDR1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PD
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PLM
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PLEPS
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PTKE
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PEMOIST
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PDTDZ
  REAL, DIMENSION(SIZE(PD,1),SIZE(PD,2),SIZE(PD,3)) :: D_M3_TH2_W2R_O_DDTDZ
  INTEGER :: IKB, IKE
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_TH2_W2R_O_DDTDZ',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PD,3)-JPVEXT_TURB

D_M3_TH2_W2R_O_DDTDZ(:,:,:) = 0.75*CSTURB%XCTV**2*PLM*PLEPS/PTKE &
 * MZF(PBLL_O_E*PEMOIST/PD*PDTDZ*(2.-PREDTH1*(1.5+PREDTH1+PREDR1)/PD), KKA, KKU, KKL)
!
D_M3_TH2_W2R_O_DDTDZ(:,:,IKB-1)=D_M3_TH2_W2R_O_DDTDZ(:,:,IKB)
D_M3_TH2_W2R_O_DDTDZ(:,:,IKE+1)=D_M3_TH2_W2R_O_DDTDZ(:,:,IKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_TH2_W2R_O_DDTDZ',1,ZHOOK_HANDLE)
END FUNCTION D_M3_TH2_W2R_O_DDTDZ
!----------------------------------------------------------------------------
FUNCTION M3_TH2_WR2(CSTURB,KKA,KKU,KKL,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PEMOIST,PDTDZ)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  INTEGER,                INTENT(IN) :: KKA 
  INTEGER,                INTENT(IN) :: KKU  
  INTEGER,                INTENT(IN) :: KKL
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PD
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PLEPS
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PEMOIST
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PDTDZ
  REAL, DIMENSION(SIZE(PD,1),SIZE(PD,2),SIZE(PD,3)) :: M3_TH2_WR2
  INTEGER :: IKB, IKE
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_TH2_WR2',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PD,3)-JPVEXT_TURB

M3_TH2_WR2(:,:,:) = 0.25*CSTURB%XCTV**2*MZF((PBLL_O_E*PEMOIST*PDTDZ)**2/PD, KKA, KKU, KKL)*PLEPS/PSQRT_TKE/CSTURB%XCTD
!
M3_TH2_WR2(:,:,IKB-1)=M3_TH2_WR2(:,:,IKB)
M3_TH2_WR2(:,:,IKE+1)=M3_TH2_WR2(:,:,IKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_TH2_WR2',1,ZHOOK_HANDLE)
END FUNCTION M3_TH2_WR2
!----------------------------------------------------------------------------
FUNCTION D_M3_TH2_WR2_O_DDTDZ(CSTURB,KKA,KKU,KKL,PREDTH1,PREDR1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PEMOIST,PDTDZ)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  INTEGER,                INTENT(IN) :: KKA 
  INTEGER,                INTENT(IN) :: KKU  
  INTEGER,                INTENT(IN) :: KKL
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDR1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PD
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PLEPS
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PEMOIST
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PDTDZ
  REAL, DIMENSION(SIZE(PD,1),SIZE(PD,2),SIZE(PD,3)) :: D_M3_TH2_WR2_O_DDTDZ
  INTEGER :: IKB, IKE
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_TH2_WR2_O_DDTDZ',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PD,3)-JPVEXT_TURB

D_M3_TH2_WR2_O_DDTDZ(:,:,:) = 0.25*CSTURB%XCTV**2*PLEPS/PSQRT_TKE/CSTURB%XCTD &
  *  MZF((PBLL_O_E*PEMOIST)**2*PDTDZ/PD*(2.-PREDTH1*(1.5+PREDTH1+PREDR1)/PD), KKA, KKU, KKL)
!
D_M3_TH2_WR2_O_DDTDZ(:,:,IKB-1)=D_M3_TH2_WR2_O_DDTDZ(:,:,IKB)
D_M3_TH2_WR2_O_DDTDZ(:,:,IKE+1)=D_M3_TH2_WR2_O_DDTDZ(:,:,IKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_TH2_WR2_O_DDTDZ',1,ZHOOK_HANDLE)
END FUNCTION D_M3_TH2_WR2_O_DDTDZ
!----------------------------------------------------------------------------
FUNCTION M3_TH2_WTHR(CSTURB,KKA,KKU,KKL,PREDR1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PEMOIST,PDTDZ)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  INTEGER,                INTENT(IN) :: KKA 
  INTEGER,                INTENT(IN) :: KKU  
  INTEGER,                INTENT(IN) :: KKL
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDR1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PD
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PLEPS
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PEMOIST
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PDTDZ
  REAL, DIMENSION(SIZE(PD,1),SIZE(PD,2),SIZE(PD,3)) :: M3_TH2_WTHR
  INTEGER :: IKB, IKE
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_TH2_WTHR',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PD,3)-JPVEXT_TURB

M3_TH2_WTHR(:,:,:) = - 0.5*CSTURB%XCTV*PLEPS/PSQRT_TKE/CSTURB%XCTD &
 * MZF(PBLL_O_E*PEMOIST*PDTDZ*(1.+PREDR1)/PD, KKA, KKU, KKL)
!
M3_TH2_WTHR(:,:,IKB-1)=M3_TH2_WTHR(:,:,IKB)
M3_TH2_WTHR(:,:,IKE+1)=M3_TH2_WTHR(:,:,IKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_TH2_WTHR',1,ZHOOK_HANDLE)
END FUNCTION M3_TH2_WTHR
!----------------------------------------------------------------------------
FUNCTION D_M3_TH2_WTHR_O_DDTDZ(CSTURB,KKA,KKU,KKL,PREDTH1,PREDR1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PEMOIST,PDTDZ)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  INTEGER,                INTENT(IN) :: KKA 
  INTEGER,                INTENT(IN) :: KKU  
  INTEGER,                INTENT(IN) :: KKL
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDR1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PD
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PLEPS
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PEMOIST
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PDTDZ
  REAL, DIMENSION(SIZE(PD,1),SIZE(PD,2),SIZE(PD,3)) :: D_M3_TH2_WTHR_O_DDTDZ
  INTEGER :: IKB, IKE
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_TH2_WTHR_O_DDTDZ',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PD,3)-JPVEXT_TURB

D_M3_TH2_WTHR_O_DDTDZ(:,:,:) = - 0.5*CSTURB%XCTV*PLEPS/PSQRT_TKE/CSTURB%XCTD &
 * MZF(PBLL_O_E*PEMOIST*(1.+PREDR1)/PD * (1. -PREDTH1*(1.5+PREDTH1+PREDR1)/PD), KKA, KKU, KKL)
!
D_M3_TH2_WTHR_O_DDTDZ(:,:,IKB-1)=D_M3_TH2_WTHR_O_DDTDZ(:,:,IKB)
D_M3_TH2_WTHR_O_DDTDZ(:,:,IKE+1)=D_M3_TH2_WTHR_O_DDTDZ(:,:,IKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_TH2_WTHR_O_DDTDZ',1,ZHOOK_HANDLE)
END FUNCTION D_M3_TH2_WTHR_O_DDTDZ
!----------------------------------------------------------------------------
FUNCTION M3_THR_WTHR(CSTURB,KKA,KKU,KKL,PREDTH1,PREDR1,PD,PLEPS,PSQRT_TKE)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  INTEGER,                INTENT(IN) :: KKA 
  INTEGER,                INTENT(IN) :: KKU  
  INTEGER,                INTENT(IN) :: KKL
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDR1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PD
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PLEPS
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(SIZE(PD,1),SIZE(PD,2),SIZE(PD,3)) :: M3_THR_WTHR
  INTEGER :: IKB, IKE
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_THR_WTHR',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PD,3)-JPVEXT_TURB

M3_THR_WTHR(:,:,:) = 0.5*PLEPS/PSQRT_TKE/CSTURB%XCTD &
 * MZF((1.+PREDTH1)*(1.+PREDR1)/PD, KKA, KKU, KKL)
!
M3_THR_WTHR(:,:,IKB-1)=M3_THR_WTHR(:,:,IKB)
M3_THR_WTHR(:,:,IKE+1)=M3_THR_WTHR(:,:,IKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_THR_WTHR',1,ZHOOK_HANDLE)
END FUNCTION M3_THR_WTHR
!----------------------------------------------------------------------------
FUNCTION D_M3_THR_WTHR_O_DDTDZ(CSTURB,KKA,KKU,KKL,PREDTH1,PREDR1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PETHETA)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  INTEGER,                INTENT(IN) :: KKA 
  INTEGER,                INTENT(IN) :: KKU  
  INTEGER,                INTENT(IN) :: KKL
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDR1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PD
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PLEPS
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PETHETA
  REAL, DIMENSION(SIZE(PD,1),SIZE(PD,2),SIZE(PD,3)) :: D_M3_THR_WTHR_O_DDTDZ
  INTEGER :: IKB, IKE
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_THR_WTHR_O_DDTDZ',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PD,3)-JPVEXT_TURB

D_M3_THR_WTHR_O_DDTDZ(:,:,:) = 0.5*PLEPS/PSQRT_TKE/CSTURB%XCTD * CSTURB%XCTV &
 * MZF(PETHETA*PBLL_O_E/PD*(1.+PREDR1)*(1.-(1.+PREDTH1)*(1.5+PREDTH1+PREDR1)/PD), KKA, KKU, KKL)
!
D_M3_THR_WTHR_O_DDTDZ(:,:,IKB-1)=D_M3_THR_WTHR_O_DDTDZ(:,:,IKB)
D_M3_THR_WTHR_O_DDTDZ(:,:,IKE+1)=D_M3_THR_WTHR_O_DDTDZ(:,:,IKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_THR_WTHR_O_DDTDZ',1,ZHOOK_HANDLE)
END FUNCTION D_M3_THR_WTHR_O_DDTDZ
!----------------------------------------------------------------------------
FUNCTION M3_THR_WTH2(CSTURB,KKA,KKU,KKL,PREDR1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PETHETA,PDRDZ)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  INTEGER,                INTENT(IN) :: KKA 
  INTEGER,                INTENT(IN) :: KKU  
  INTEGER,                INTENT(IN) :: KKL
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDR1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PD
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PLEPS
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PETHETA
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PDRDZ
  REAL, DIMENSION(SIZE(PD,1),SIZE(PD,2),SIZE(PD,3)) :: M3_THR_WTH2
  INTEGER :: IKB, IKE
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_THR_WTH2',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PD,3)-JPVEXT_TURB

M3_THR_WTH2(:,:,:) = - 0.25*PLEPS/PSQRT_TKE/CSTURB%XCTD*CSTURB%XCTV &
 * MZF((1.+PREDR1)*PBLL_O_E*PETHETA*PDRDZ/PD, KKA, KKU, KKL)
!
M3_THR_WTH2(:,:,IKB-1)=M3_THR_WTH2(:,:,IKB)
M3_THR_WTH2(:,:,IKE+1)=M3_THR_WTH2(:,:,IKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_THR_WTH2',1,ZHOOK_HANDLE)
END FUNCTION M3_THR_WTH2
!----------------------------------------------------------------------------
FUNCTION D_M3_THR_WTH2_O_DDTDZ(CSTURB,KKA,KKU,KKL,PREDTH1,PREDR1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PETHETA,PDRDZ)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  INTEGER,                INTENT(IN) :: KKA 
  INTEGER,                INTENT(IN) :: KKU  
  INTEGER,                INTENT(IN) :: KKL
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDR1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PD
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PLEPS
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PETHETA
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PDRDZ
  REAL, DIMENSION(SIZE(PD,1),SIZE(PD,2),SIZE(PD,3)) :: D_M3_THR_WTH2_O_DDTDZ
  INTEGER :: IKB, IKE
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_THR_WTH2_O_DDTDZ',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PD,3)-JPVEXT_TURB

D_M3_THR_WTH2_O_DDTDZ(:,:,:) = - 0.25*PLEPS/PSQRT_TKE/CSTURB%XCTD*CSTURB%XCTV**2 &
 * MZF(-(1.+PREDR1)*(PBLL_O_E*PETHETA/PD)**2*PDRDZ*(1.5+PREDTH1+PREDR1), KKA, KKU, KKL)
!
D_M3_THR_WTH2_O_DDTDZ(:,:,IKB-1)=D_M3_THR_WTH2_O_DDTDZ(:,:,IKB)
D_M3_THR_WTH2_O_DDTDZ(:,:,IKE+1)=D_M3_THR_WTH2_O_DDTDZ(:,:,IKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_THR_WTH2_O_DDTDZ',1,ZHOOK_HANDLE)
END FUNCTION D_M3_THR_WTH2_O_DDTDZ
!----------------------------------------------------------------------------
FUNCTION D_M3_THR_WTH2_O_DDRDZ(CSTURB,KKA,KKU,KKL,PREDTH1,PREDR1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PETHETA)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  INTEGER,                INTENT(IN) :: KKA 
  INTEGER,                INTENT(IN) :: KKU  
  INTEGER,                INTENT(IN) :: KKL
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDR1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PD
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PLEPS
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PETHETA
  REAL, DIMENSION(SIZE(PD,1),SIZE(PD,2),SIZE(PD,3)) :: D_M3_THR_WTH2_O_DDRDZ
  INTEGER :: IKB, IKE
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_THR_WTH2_O_DDRDZ',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PD,3)-JPVEXT_TURB

D_M3_THR_WTH2_O_DDRDZ(:,:,:) = - 0.25*PLEPS/PSQRT_TKE/CSTURB%XCTD*CSTURB%XCTV          &
 * MZF(PBLL_O_E*PETHETA/PD                                              &
       *(-(1.+PREDR1)*PREDR1/PD*(1.5+PREDTH1+PREDR1)+(1.+2.*PREDR1)),     &
       KKA, KKU, KKL)
!
D_M3_THR_WTH2_O_DDRDZ(:,:,IKB-1)=D_M3_THR_WTH2_O_DDRDZ(:,:,IKB)
D_M3_THR_WTH2_O_DDRDZ(:,:,IKE+1)=D_M3_THR_WTH2_O_DDRDZ(:,:,IKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_THR_WTH2_O_DDRDZ',1,ZHOOK_HANDLE)
END FUNCTION D_M3_THR_WTH2_O_DDRDZ
!----------------------------------------------------------------------------
FUNCTION M3_THR_W2TH(CSTURB,KKA,KKU,KKL,PREDR1,PD,PLM,PLEPS,PTKE,PDRDZ)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  INTEGER,                INTENT(IN) :: KKA 
  INTEGER,                INTENT(IN) :: KKU  
  INTEGER,                INTENT(IN) :: KKL
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDR1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PD
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PLM
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PLEPS
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PTKE
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PDRDZ
  REAL, DIMENSION(SIZE(PD,1),SIZE(PD,2),SIZE(PD,3)) :: M3_THR_W2TH
  INTEGER :: IKB, IKE
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_THR_W2TH',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PD,3)-JPVEXT_TURB

M3_THR_W2TH(:,:,:) = - 0.75*PLM*PLEPS/PTKE * CSTURB%XCTV      &
 * MZF((1.+PREDR1)*PDRDZ/PD, KKA, KKU, KKL)
!
M3_THR_W2TH(:,:,IKB-1)=M3_THR_W2TH(:,:,IKB)
M3_THR_W2TH(:,:,IKE+1)=M3_THR_W2TH(:,:,IKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_THR_W2TH',1,ZHOOK_HANDLE)
END FUNCTION M3_THR_W2TH
!----------------------------------------------------------------------------
FUNCTION D_M3_THR_W2TH_O_DDTDZ(CSTURB,KKA,KKU,KKL,PREDTH1,PREDR1,PD,PLM,PLEPS,PTKE,PBLL_O_E,PDRDZ,PETHETA)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  INTEGER,                INTENT(IN) :: KKA 
  INTEGER,                INTENT(IN) :: KKU  
  INTEGER,                INTENT(IN) :: KKL
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDR1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PD
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PLM
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PLEPS
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PTKE
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PDRDZ
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PETHETA
  REAL, DIMENSION(SIZE(PD,1),SIZE(PD,2),SIZE(PD,3)) :: D_M3_THR_W2TH_O_DDTDZ
  INTEGER :: IKB, IKE
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_THR_W2TH_O_DDTDZ',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PD,3)-JPVEXT_TURB

D_M3_THR_W2TH_O_DDTDZ(:,:,:) = - 0.75*PLM*PLEPS/PTKE * CSTURB%XCTV**2    &
 * MZF(-PETHETA*PBLL_O_E*(1.+PREDR1)*PDRDZ*(1.5+PREDTH1+PREDR1)/PD**2, KKA, KKU, KKL)

!
D_M3_THR_W2TH_O_DDTDZ(:,:,IKB-1)=D_M3_THR_W2TH_O_DDTDZ(:,:,IKB)
D_M3_THR_W2TH_O_DDTDZ(:,:,IKE+1)=D_M3_THR_W2TH_O_DDTDZ(:,:,IKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_THR_W2TH_O_DDTDZ',1,ZHOOK_HANDLE)
END FUNCTION D_M3_THR_W2TH_O_DDTDZ
!----------------------------------------------------------------------------
FUNCTION D_M3_THR_W2TH_O_DDRDZ(CSTURB,KKA,KKU,KKL,PREDTH1,PREDR1,PD,PLM,PLEPS,PTKE)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  INTEGER,                INTENT(IN) :: KKA 
  INTEGER,                INTENT(IN) :: KKU  
  INTEGER,                INTENT(IN) :: KKL
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDR1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PD
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PLM
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PLEPS
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PTKE
  REAL, DIMENSION(SIZE(PD,1),SIZE(PD,2),SIZE(PD,3)) :: D_M3_THR_W2TH_O_DDRDZ
  INTEGER :: IKB, IKE
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_THR_W2TH_O_DDRDZ',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PD,3)-JPVEXT_TURB

D_M3_THR_W2TH_O_DDRDZ(:,:,:) = - 0.75*PLM*PLEPS/PTKE * CSTURB%XCTV     &
 * MZF(-(1.+PREDR1)*PREDR1*(1.5+PREDTH1+PREDR1)/PD**2          &
        +(1.+2.*PREDR1)/PD,                                    &
       KKA, KKU, KKL)

!
D_M3_THR_W2TH_O_DDRDZ(:,:,IKB-1)=D_M3_THR_W2TH_O_DDRDZ(:,:,IKB)
D_M3_THR_W2TH_O_DDRDZ(:,:,IKE+1)=D_M3_THR_W2TH_O_DDRDZ(:,:,IKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_THR_W2TH_O_DDRDZ',1,ZHOOK_HANDLE)
END FUNCTION D_M3_THR_W2TH_O_DDRDZ
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!
FUNCTION PSI3(CSTURB,PREDR1,PREDTH1,PRED2R3,PRED2TH3,PRED2THR3,HTURBDIM,OUSERV)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDR1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PRED2TH3
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PRED2R3
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PRED2THR3
  CHARACTER(len=4),       INTENT(IN) :: HTURBDIM  ! 1DIM or 3DIM turb. scheme
  LOGICAL,                INTENT(IN) :: OUSERV    ! flag to use vapor
  REAL, DIMENSION(SIZE(PREDTH1,1),SIZE(PREDTH1,2),SIZE(PREDTH1,3)) :: PSI3
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:PSI3',0,ZHOOK_HANDLE)
PSI3 = PHI3(CSTURB,PREDR1,PREDTH1,PRED2R3,PRED2TH3,PRED2THR3,HTURBDIM,OUSERV)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:PSI3',1,ZHOOK_HANDLE)
END FUNCTION PSI3
!----------------------------------------------------------------------------
FUNCTION D_PSI3DRDZ_O_DDRDZ(CSTURB,PPSI3,PREDR1,PREDTH1,PRED2R3,PRED2THR3,HTURBDIM,OUSERV)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PPSI3
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDR1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PRED2R3
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PRED2THR3
  CHARACTER(len=4),       INTENT(IN) :: HTURBDIM  ! 1DIM or 3DIM turb. scheme
  LOGICAL,                INTENT(IN) :: OUSERV    ! flag to use vapor
  REAL, DIMENSION(SIZE(PREDTH1,1),SIZE(PREDTH1,2),SIZE(PREDTH1,3)) :: D_PSI3DRDZ_O_DDRDZ

REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_PSI3DRDZ_O_DDRDZ',0,ZHOOK_HANDLE)
D_PSI3DRDZ_O_DDRDZ = D_PHI3DTDZ_O_DDTDZ(CSTURB,PPSI3,PREDR1,PREDTH1,PRED2R3,PRED2THR3,HTURBDIM,OUSERV)
!
!C'est ok?!
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_PSI3DRDZ_O_DDRDZ',1,ZHOOK_HANDLE)
END FUNCTION D_PSI3DRDZ_O_DDRDZ
!----------------------------------------------------------------------------
FUNCTION D_PSI3DTDZ_O_DDTDZ(CSTURB,PPSI3,PREDR1,PREDTH1,PRED2R3,PRED2THR3,HTURBDIM,OUSERV)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PPSI3
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDR1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PRED2R3
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PRED2THR3
  CHARACTER(len=4),       INTENT(IN) :: HTURBDIM  ! 1DIM or 3DIM turb. scheme
  LOGICAL,                INTENT(IN) :: OUSERV    ! flag to use vapor
  REAL, DIMENSION(SIZE(PREDTH1,1),SIZE(PREDTH1,2),SIZE(PREDTH1,3)) :: D_PSI3DTDZ_O_DDTDZ
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_PSI3DTDZ_O_DDTDZ',0,ZHOOK_HANDLE)
D_PSI3DTDZ_O_DDTDZ = D_PHI3DRDZ_O_DDRDZ(CSTURB,PPSI3,PREDR1,PREDTH1,PRED2R3,PRED2THR3,HTURBDIM,OUSERV)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_PSI3DTDZ_O_DDTDZ',1,ZHOOK_HANDLE)
END FUNCTION D_PSI3DTDZ_O_DDTDZ
!----------------------------------------------------------------------------
FUNCTION D_PSI3DRDZ2_O_DDRDZ(CSTURB,PPSI3,PREDR1,PREDTH1,PRED2R3,PRED2THR3,PDRDZ,HTURBDIM,OUSERV)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PPSI3
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDR1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PRED2R3
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PRED2THR3
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PDRDZ
  CHARACTER(len=4),       INTENT(IN) :: HTURBDIM  ! 1DIM or 3DIM turb. scheme
  LOGICAL,                INTENT(IN) :: OUSERV    ! flag to use vapor
  REAL, DIMENSION(SIZE(PREDTH1,1),SIZE(PREDTH1,2),SIZE(PREDTH1,3)) :: D_PSI3DRDZ2_O_DDRDZ
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_PSI3DRDZ2_O_DDRDZ',0,ZHOOK_HANDLE)
D_PSI3DRDZ2_O_DDRDZ = D_PHI3DTDZ2_O_DDTDZ(CSTURB,PPSI3,PREDR1,PREDTH1,PRED2R3,PRED2THR3,PDRDZ,HTURBDIM,OUSERV)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_PSI3DRDZ2_O_DDRDZ',1,ZHOOK_HANDLE)
END FUNCTION D_PSI3DRDZ2_O_DDRDZ
!----------------------------------------------------------------------------
FUNCTION M3_WR_WR2(CSTURB,PREDR1,PREDTH1,PD,PBLL_O_E,PEMOIST)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDR1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PD
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PEMOIST
  REAL, DIMENSION(SIZE(PD,1),SIZE(PD,2),SIZE(PD,3)) :: M3_WR_WR2
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_WR_WR2',0,ZHOOK_HANDLE)
M3_WR_WR2 = M3_WTH_WTH2(CSTURB,PREDR1,PREDTH1,PD,PBLL_O_E,PEMOIST)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_WR_WR2',1,ZHOOK_HANDLE)
END FUNCTION M3_WR_WR2
!----------------------------------------------------------------------------
FUNCTION D_M3_WR_WR2_O_DDRDZ(CSTURB,PM3_WR_WR2,PREDR1,PREDTH1,PD,PBLL_O_E,PEMOIST)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PM3_WR_WR2
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDR1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PD
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PEMOIST
  REAL, DIMENSION(SIZE(PD,1),SIZE(PD,2),SIZE(PD,3)) :: D_M3_WR_WR2_O_DDRDZ
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_WR_WR2_O_DDRDZ',0,ZHOOK_HANDLE)
D_M3_WR_WR2_O_DDRDZ = D_M3_WTH_WTH2_O_DDTDZ(CSTURB,PM3_WR_WR2,PREDR1,PREDTH1,PD,PBLL_O_E,PEMOIST)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_WR_WR2_O_DDRDZ',1,ZHOOK_HANDLE)
END FUNCTION D_M3_WR_WR2_O_DDRDZ
!----------------------------------------------------------------------------
FUNCTION M3_WR_W2R(CSTURB,KKA,KKU,KKL,PREDR1,PREDTH1,PD,PKEFF,PTKE)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  INTEGER,                INTENT(IN) :: KKA 
  INTEGER,                INTENT(IN) :: KKU  
  INTEGER,                INTENT(IN) :: KKL
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDR1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PD
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PKEFF
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PTKE
  REAL, DIMENSION(SIZE(PD,1),SIZE(PD,2),SIZE(PD,3)) :: M3_WR_W2R
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_WR_W2R',0,ZHOOK_HANDLE)
M3_WR_W2R = M3_WTH_W2TH(CSTURB,KKA,KKU,KKL,PREDR1,PREDTH1,PD,PKEFF,PTKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_WR_W2R',1,ZHOOK_HANDLE)
END FUNCTION M3_WR_W2R
!----------------------------------------------------------------------------
FUNCTION D_M3_WR_W2R_O_DDRDZ(CSTURB,KKA,KKU,KKL,PREDR1,PREDTH1,PD,PBLL_O_E,PEMOIST,PKEFF,PTKE)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  INTEGER,                INTENT(IN) :: KKA 
  INTEGER,                INTENT(IN) :: KKU  
  INTEGER,                INTENT(IN) :: KKL
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDR1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PD
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PEMOIST
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PKEFF
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PTKE
  REAL, DIMENSION(SIZE(PD,1),SIZE(PD,2),SIZE(PD,3)) :: D_M3_WR_W2R_O_DDRDZ
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_WR_W2R_O_DDRDZ',0,ZHOOK_HANDLE)
D_M3_WR_W2R_O_DDRDZ = D_M3_WTH_W2TH_O_DDTDZ(CSTURB,KKA,KKU,KKL,PREDR1,PREDTH1,PD,PBLL_O_E,PEMOIST,PKEFF,PTKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_WR_W2R_O_DDRDZ',1,ZHOOK_HANDLE)
END FUNCTION D_M3_WR_W2R_O_DDRDZ
!----------------------------------------------------------------------------
FUNCTION M3_WR_W2TH(CSTURB,KKA,KKU,KKL,PD,PKEFF,PTKE,PBLL_O_E,PETHETA,PDRDZ)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  INTEGER,                INTENT(IN) :: KKA 
  INTEGER,                INTENT(IN) :: KKU  
  INTEGER,                INTENT(IN) :: KKL
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PD
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PKEFF
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PTKE
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PETHETA
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PDRDZ
  REAL, DIMENSION(SIZE(PD,1),SIZE(PD,2),SIZE(PD,3)) :: M3_WR_W2TH
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_WR_W2TH',0,ZHOOK_HANDLE)
M3_WR_W2TH = M3_WTH_W2R(CSTURB,KKA,KKU,KKL,PD,PKEFF,PTKE,PBLL_O_E,PETHETA,PDRDZ)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_WR_W2TH',1,ZHOOK_HANDLE)
END FUNCTION M3_WR_W2TH
!----------------------------------------------------------------------------
FUNCTION D_M3_WR_W2TH_O_DDRDZ(CSTURB,KKA,KKU,KKL,PREDR1,PREDTH1,PD,PKEFF,PTKE,PBLL_O_E,PETHETA)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  INTEGER,                INTENT(IN) :: KKA 
  INTEGER,                INTENT(IN) :: KKU  
  INTEGER,                INTENT(IN) :: KKL
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDR1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PD
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PKEFF
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PTKE
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PETHETA
  REAL, DIMENSION(SIZE(PD,1),SIZE(PD,2),SIZE(PD,3)) :: D_M3_WR_W2TH_O_DDRDZ
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_WR_W2TH_O_DDRDZ',0,ZHOOK_HANDLE)
D_M3_WR_W2TH_O_DDRDZ = D_M3_WTH_W2R_O_DDTDZ(CSTURB,KKA,KKU,KKL,PREDR1,PREDTH1,PD,PKEFF,PTKE,PBLL_O_E,PETHETA)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_WR_W2TH_O_DDRDZ',1,ZHOOK_HANDLE)
END FUNCTION D_M3_WR_W2TH_O_DDRDZ
!----------------------------------------------------------------------------
FUNCTION M3_WR_WTH2(CSTURB,KKA,KKU,KKL,PD,PKEFF,PTKE,PSQRT_TKE,PBLL_O_E,PBETA,PLEPS,PETHETA,PDRDZ)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  INTEGER,                INTENT(IN) :: KKA 
  INTEGER,                INTENT(IN) :: KKU  
  INTEGER,                INTENT(IN) :: KKL
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PD
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PKEFF
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PTKE
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PBETA
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PLEPS
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PETHETA
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PDRDZ
  REAL, DIMENSION(SIZE(PD,1),SIZE(PD,2),SIZE(PD,3)) :: M3_WR_WTH2
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_WR_WTH2',0,ZHOOK_HANDLE)
M3_WR_WTH2 = M3_WTH_WR2(CSTURB,KKA,KKU,KKL,PD,PKEFF,PTKE,PSQRT_TKE,PBLL_O_E,PBETA,PLEPS,PETHETA,PDRDZ)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_WR_WTH2',1,ZHOOK_HANDLE)
END FUNCTION M3_WR_WTH2
!----------------------------------------------------------------------------
FUNCTION D_M3_WR_WTH2_O_DDRDZ(CSTURB,KKA,KKU,KKL,PREDR1,PREDTH1,PD,PKEFF,PTKE,PSQRT_TKE,PBLL_O_E,PBETA,PLEPS,PETHETA)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  INTEGER,                INTENT(IN) :: KKA 
  INTEGER,                INTENT(IN) :: KKU  
  INTEGER,                INTENT(IN) :: KKL
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDR1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PD
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PKEFF
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PTKE
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PBETA
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PLEPS
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PETHETA
  REAL, DIMENSION(SIZE(PD,1),SIZE(PD,2),SIZE(PD,3)) :: D_M3_WR_WTH2_O_DDRDZ
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_WR_WTH2_O_DDRDZ',0,ZHOOK_HANDLE)
D_M3_WR_WTH2_O_DDRDZ = D_M3_WTH_WR2_O_DDTDZ(CSTURB,KKA,KKU,KKL,PREDR1,PREDTH1,PD,PKEFF,PTKE,PSQRT_TKE,PBLL_O_E,PBETA,PLEPS,PETHETA)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_WR_WTH2_O_DDRDZ',1,ZHOOK_HANDLE)
END FUNCTION D_M3_WR_WTH2_O_DDRDZ
!----------------------------------------------------------------------------
FUNCTION M3_WR_WTHR(CSTURB,KKA,KKU,KKL,PREDTH1,PD,PKEFF,PTKE,PSQRT_TKE,PBETA,PLEPS,PETHETA)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  INTEGER,                INTENT(IN) :: KKA 
  INTEGER,                INTENT(IN) :: KKU  
  INTEGER,                INTENT(IN) :: KKL
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PD
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PKEFF
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PTKE
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PBETA
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PLEPS
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PETHETA
  REAL, DIMENSION(SIZE(PD,1),SIZE(PD,2),SIZE(PD,3)) :: M3_WR_WTHR
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_WR_WTHR',0,ZHOOK_HANDLE)
M3_WR_WTHR = M3_WTH_WTHR(CSTURB,KKA,KKU,KKL,PREDTH1,PD,PKEFF,PTKE,PSQRT_TKE,PBETA,PLEPS,PETHETA)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_WR_WTHR',1,ZHOOK_HANDLE)
END FUNCTION M3_WR_WTHR
!----------------------------------------------------------------------------
FUNCTION D_M3_WR_WTHR_O_DDRDZ(CSTURB,KKA,KKU,KKL,PM3_WR_WTHR,PREDR1,PREDTH1,PD,PBLL_O_E,PEMOIST)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  INTEGER,                INTENT(IN) :: KKA 
  INTEGER,                INTENT(IN) :: KKU  
  INTEGER,                INTENT(IN) :: KKL
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PM3_WR_WTHR
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDR1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PD
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PEMOIST
  REAL, DIMENSION(SIZE(PD,1),SIZE(PD,2),SIZE(PD,3)) :: D_M3_WR_WTHR_O_DDRDZ
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_WR_WTHR_O_DDRDZ',0,ZHOOK_HANDLE)
D_M3_WR_WTHR_O_DDRDZ = D_M3_WTH_WTHR_O_DDTDZ(CSTURB,PM3_WR_WTHR,PREDR1,PREDTH1,PD,PBLL_O_E,PEMOIST)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_WR_WTHR_O_DDRDZ',1,ZHOOK_HANDLE)
END FUNCTION D_M3_WR_WTHR_O_DDRDZ
!----------------------------------------------------------------------------
FUNCTION M3_R2_W2R(CSTURB,KKA,KKU,KKL,PREDR1,PREDTH1,PD,PDRDZ,PLM,PLEPS,PTKE)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  INTEGER,                INTENT(IN) :: KKA 
  INTEGER,                INTENT(IN) :: KKU  
  INTEGER,                INTENT(IN) :: KKL
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDR1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PD
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PDRDZ
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PLM
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PLEPS
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PTKE
  REAL, DIMENSION(SIZE(PD,1),SIZE(PD,2),SIZE(PD,3)) :: M3_R2_W2R
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_R2_W2R',0,ZHOOK_HANDLE)
M3_R2_W2R = M3_TH2_W2TH(CSTURB,KKA,KKU,KKL,PREDR1,PREDTH1,PD,PDRDZ,PLM,PLEPS,PTKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_R2_W2R',1,ZHOOK_HANDLE)
END FUNCTION M3_R2_W2R
!----------------------------------------------------------------------------
FUNCTION D_M3_R2_W2R_O_DDRDZ(CSTURB,KKA,KKU,KKL,PREDR1,PREDTH1,PD,PLM,PLEPS,PTKE,OUSERV)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  INTEGER,                INTENT(IN) :: KKA 
  INTEGER,                INTENT(IN) :: KKU  
  INTEGER,                INTENT(IN) :: KKL
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDR1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PD
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PLM
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PLEPS
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PTKE
  LOGICAL,                INTENT(IN) :: OUSERV
  REAL, DIMENSION(SIZE(PD,1),SIZE(PD,2),SIZE(PD,3)) :: D_M3_R2_W2R_O_DDRDZ
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_R2_W2R_O_DDRDZ',0,ZHOOK_HANDLE)
D_M3_R2_W2R_O_DDRDZ = D_M3_TH2_W2TH_O_DDTDZ(CSTURB,KKA,KKU,KKL,PREDR1,PREDTH1,PD,PLM,PLEPS,PTKE,OUSERV)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_R2_W2R_O_DDRDZ',1,ZHOOK_HANDLE)
END FUNCTION D_M3_R2_W2R_O_DDRDZ
!----------------------------------------------------------------------------
FUNCTION M3_R2_WR2(CSTURB,KKA,KKU,KKL,PREDR1,PREDTH1,PD,PLEPS,PSQRT_TKE)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  INTEGER,                INTENT(IN) :: KKA 
  INTEGER,                INTENT(IN) :: KKU  
  INTEGER,                INTENT(IN) :: KKL
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDR1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PD
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PLEPS
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(SIZE(PD,1),SIZE(PD,2),SIZE(PD,3)) :: M3_R2_WR2
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_R2_WR2',0,ZHOOK_HANDLE)
M3_R2_WR2 = M3_TH2_WTH2(CSTURB,KKA,KKU,KKL,PREDR1,PREDTH1,PD,PLEPS,PSQRT_TKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_R2_WR2',1,ZHOOK_HANDLE)
END FUNCTION M3_R2_WR2
!----------------------------------------------------------------------------
FUNCTION D_M3_R2_WR2_O_DDRDZ(CSTURB,KKA,KKU,KKL,PREDR1,PREDTH1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PEMOIST)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  INTEGER,                INTENT(IN) :: KKA 
  INTEGER,                INTENT(IN) :: KKU  
  INTEGER,                INTENT(IN) :: KKL
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDR1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PD
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PLEPS
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PEMOIST
  REAL, DIMENSION(SIZE(PD,1),SIZE(PD,2),SIZE(PD,3)) :: D_M3_R2_WR2_O_DDRDZ
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_R2_WR2_O_DDRDZ',0,ZHOOK_HANDLE)
D_M3_R2_WR2_O_DDRDZ = D_M3_TH2_WTH2_O_DDTDZ(CSTURB,KKA,KKU,KKL,PREDR1,PREDTH1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PEMOIST)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_R2_WR2_O_DDRDZ',1,ZHOOK_HANDLE)
END FUNCTION D_M3_R2_WR2_O_DDRDZ
!----------------------------------------------------------------------------
FUNCTION M3_R2_W2TH(CSTURB,KKA,KKU,KKL,PD,PLM,PLEPS,PTKE,PBLL_O_E,PETHETA,PDRDZ)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  INTEGER,                INTENT(IN) :: KKA 
  INTEGER,                INTENT(IN) :: KKU  
  INTEGER,                INTENT(IN) :: KKL
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PD
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PLM
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PLEPS
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PTKE
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PETHETA
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PDRDZ
  REAL, DIMENSION(SIZE(PD,1),SIZE(PD,2),SIZE(PD,3)) :: M3_R2_W2TH
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_R2_W2TH',0,ZHOOK_HANDLE)
M3_R2_W2TH = M3_TH2_W2R(CSTURB,KKA,KKU,KKL,PD,PLM,PLEPS,PTKE,PBLL_O_E,PETHETA,PDRDZ)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_R2_W2TH',1,ZHOOK_HANDLE)
END FUNCTION M3_R2_W2TH
!----------------------------------------------------------------------------
FUNCTION D_M3_R2_W2TH_O_DDRDZ(CSTURB,KKA,KKU,KKL,PREDR1,PREDTH1,PD,PLM,PLEPS,PTKE,PBLL_O_E,PETHETA,PDRDZ)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  INTEGER,                INTENT(IN) :: KKA 
  INTEGER,                INTENT(IN) :: KKU  
  INTEGER,                INTENT(IN) :: KKL
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDR1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PD
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PLM
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PLEPS
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PTKE
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PETHETA
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PDRDZ
  REAL, DIMENSION(SIZE(PD,1),SIZE(PD,2),SIZE(PD,3)) :: D_M3_R2_W2TH_O_DDRDZ
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_R2_W2TH_O_DDRDZ',0,ZHOOK_HANDLE)
D_M3_R2_W2TH_O_DDRDZ = D_M3_TH2_W2R_O_DDTDZ(CSTURB,KKA,KKU,KKL,PREDR1,PREDTH1,PD,PLM,PLEPS,PTKE,PBLL_O_E,PETHETA,PDRDZ)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_R2_W2TH_O_DDRDZ',1,ZHOOK_HANDLE)
END FUNCTION D_M3_R2_W2TH_O_DDRDZ
!----------------------------------------------------------------------------
FUNCTION M3_R2_WTH2(CSTURB,KKA,KKU,KKL,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PETHETA,PDRDZ)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  INTEGER,                INTENT(IN) :: KKA 
  INTEGER,                INTENT(IN) :: KKU  
  INTEGER,                INTENT(IN) :: KKL
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PD
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PLEPS
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PETHETA
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PDRDZ
  REAL, DIMENSION(SIZE(PD,1),SIZE(PD,2),SIZE(PD,3)) :: M3_R2_WTH2
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_R2_WTH2',0,ZHOOK_HANDLE)
M3_R2_WTH2 = M3_TH2_WR2(CSTURB,KKA,KKU,KKL,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PETHETA,PDRDZ)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_R2_WTH2',1,ZHOOK_HANDLE)
END FUNCTION M3_R2_WTH2
!----------------------------------------------------------------------------
FUNCTION D_M3_R2_WTH2_O_DDRDZ(CSTURB,KKA,KKU,KKL,PREDR1,PREDTH1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PETHETA,PDRDZ)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  INTEGER,                INTENT(IN) :: KKA 
  INTEGER,                INTENT(IN) :: KKU  
  INTEGER,                INTENT(IN) :: KKL
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDR1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PD
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PLEPS
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PETHETA
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PDRDZ
  REAL, DIMENSION(SIZE(PD,1),SIZE(PD,2),SIZE(PD,3)) :: D_M3_R2_WTH2_O_DDRDZ
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_R2_WTH2_O_DDRDZ',0,ZHOOK_HANDLE)
D_M3_R2_WTH2_O_DDRDZ = D_M3_TH2_WR2_O_DDTDZ(CSTURB,KKA,KKU,KKL,PREDR1,PREDTH1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PETHETA,PDRDZ)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_R2_WTH2_O_DDRDZ',1,ZHOOK_HANDLE)
END FUNCTION D_M3_R2_WTH2_O_DDRDZ
!----------------------------------------------------------------------------
FUNCTION M3_R2_WTHR(CSTURB,KKA,KKU,KKL,PREDTH1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PETHETA,PDRDZ)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  INTEGER,                INTENT(IN) :: KKA 
  INTEGER,                INTENT(IN) :: KKU  
  INTEGER,                INTENT(IN) :: KKL
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PD
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PLEPS
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PETHETA
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PDRDZ
  REAL, DIMENSION(SIZE(PD,1),SIZE(PD,2),SIZE(PD,3)) :: M3_R2_WTHR
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_R2_WTHR',0,ZHOOK_HANDLE)
M3_R2_WTHR = M3_TH2_WTHR(CSTURB,KKA,KKU,KKL,PREDTH1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PETHETA,PDRDZ)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_R2_WTHR',1,ZHOOK_HANDLE)
END FUNCTION M3_R2_WTHR
!----------------------------------------------------------------------------
FUNCTION D_M3_R2_WTHR_O_DDRDZ(CSTURB,KKA,KKU,KKL,PREDR1,PREDTH1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PETHETA,PDRDZ)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  INTEGER,                INTENT(IN) :: KKA 
  INTEGER,                INTENT(IN) :: KKU  
  INTEGER,                INTENT(IN) :: KKL
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDR1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PD
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PLEPS
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PETHETA
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PDRDZ
  REAL, DIMENSION(SIZE(PD,1),SIZE(PD,2),SIZE(PD,3)) :: D_M3_R2_WTHR_O_DDRDZ
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_R2_WTHR_O_DDRDZ',0,ZHOOK_HANDLE)
D_M3_R2_WTHR_O_DDRDZ = D_M3_TH2_WTHR_O_DDTDZ(CSTURB,KKA,KKU,KKL,PREDR1,PREDTH1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PETHETA,PDRDZ)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_R2_WTHR_O_DDRDZ',1,ZHOOK_HANDLE)
END FUNCTION D_M3_R2_WTHR_O_DDRDZ
!----------------------------------------------------------------------------
FUNCTION D_M3_THR_WTHR_O_DDRDZ(CSTURB,KKA,KKU,KKL,PREDR1,PREDTH1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PEMOIST)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  INTEGER,                INTENT(IN) :: KKA 
  INTEGER,                INTENT(IN) :: KKU  
  INTEGER,                INTENT(IN) :: KKL
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDR1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PD
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PLEPS
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PEMOIST
  REAL, DIMENSION(SIZE(PD,1),SIZE(PD,2),SIZE(PD,3)) :: D_M3_THR_WTHR_O_DDRDZ
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_THR_WTHR_O_DDRDZ',0,ZHOOK_HANDLE)
D_M3_THR_WTHR_O_DDRDZ = D_M3_THR_WTHR_O_DDTDZ(CSTURB,KKA,KKU,KKL,PREDR1,PREDTH1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PEMOIST)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_THR_WTHR_O_DDRDZ',1,ZHOOK_HANDLE)
END FUNCTION D_M3_THR_WTHR_O_DDRDZ
!----------------------------------------------------------------------------
FUNCTION M3_THR_WR2(CSTURB,KKA,KKU,KKL,PREDTH1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PEMOIST,PDTDZ)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  INTEGER,                INTENT(IN) :: KKA 
  INTEGER,                INTENT(IN) :: KKU  
  INTEGER,                INTENT(IN) :: KKL
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PD
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PLEPS
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PEMOIST
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PDTDZ
  REAL, DIMENSION(SIZE(PD,1),SIZE(PD,2),SIZE(PD,3)) :: M3_THR_WR2
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_THR_WR2',0,ZHOOK_HANDLE)
M3_THR_WR2 = M3_THR_WTH2(CSTURB,KKA,KKU,KKL,PREDTH1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PEMOIST,PDTDZ)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_THR_WR2',1,ZHOOK_HANDLE)
END FUNCTION M3_THR_WR2
!----------------------------------------------------------------------------
FUNCTION D_M3_THR_WR2_O_DDRDZ(CSTURB,KKA,KKU,KKL,PREDR1,PREDTH1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PEMOIST,PDTDZ)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  INTEGER,                INTENT(IN) :: KKA 
  INTEGER,                INTENT(IN) :: KKU  
  INTEGER,                INTENT(IN) :: KKL
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDR1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PD
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PLEPS
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PEMOIST
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PDTDZ
  REAL, DIMENSION(SIZE(PD,1),SIZE(PD,2),SIZE(PD,3)) :: D_M3_THR_WR2_O_DDRDZ
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_THR_WR2_O_DDRDZ',0,ZHOOK_HANDLE)
D_M3_THR_WR2_O_DDRDZ = D_M3_THR_WTH2_O_DDTDZ(CSTURB,KKA,KKU,KKL,PREDR1,PREDTH1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PEMOIST,PDTDZ)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_THR_WR2_O_DDRDZ',1,ZHOOK_HANDLE)
END FUNCTION D_M3_THR_WR2_O_DDRDZ
!----------------------------------------------------------------------------
FUNCTION D_M3_THR_WR2_O_DDTDZ(CSTURB,KKA,KKU,KKL,PREDR1,PREDTH1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PEMOIST)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  INTEGER,                INTENT(IN) :: KKA 
  INTEGER,                INTENT(IN) :: KKU  
  INTEGER,                INTENT(IN) :: KKL
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDR1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PD
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PLEPS
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PEMOIST
  REAL, DIMENSION(SIZE(PD,1),SIZE(PD,2),SIZE(PD,3)) :: D_M3_THR_WR2_O_DDTDZ
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_THR_WR2_O_DDTDZ',0,ZHOOK_HANDLE)
D_M3_THR_WR2_O_DDTDZ = D_M3_THR_WTH2_O_DDRDZ(CSTURB,KKA,KKU,KKL,PREDR1,PREDTH1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PEMOIST)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_THR_WR2_O_DDTDZ',1,ZHOOK_HANDLE)
END FUNCTION D_M3_THR_WR2_O_DDTDZ
!----------------------------------------------------------------------------
FUNCTION M3_THR_W2R(CSTURB,KKA,KKU,KKL,PREDTH1,PD,PLM,PLEPS,PTKE,PDTDZ)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  INTEGER,                INTENT(IN) :: KKA 
  INTEGER,                INTENT(IN) :: KKU  
  INTEGER,                INTENT(IN) :: KKL
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PD
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PLM
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PLEPS
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PTKE
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PDTDZ
  REAL, DIMENSION(SIZE(PD,1),SIZE(PD,2),SIZE(PD,3)) :: M3_THR_W2R
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_THR_W2R',0,ZHOOK_HANDLE)
M3_THR_W2R = M3_THR_W2TH(CSTURB,KKA,KKU,KKL,PREDTH1,PD,PLM,PLEPS,PTKE,PDTDZ)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_THR_W2R',1,ZHOOK_HANDLE)
END FUNCTION M3_THR_W2R
!----------------------------------------------------------------------------
FUNCTION D_M3_THR_W2R_O_DDRDZ(CSTURB,KKA,KKU,KKL,PREDR1,PREDTH1,PD,PLM,PLEPS,PTKE,PBLL_O_E,PDTDZ,PEMOIST)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  INTEGER,                INTENT(IN) :: KKA 
  INTEGER,                INTENT(IN) :: KKU  
  INTEGER,                INTENT(IN) :: KKL
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDR1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PD
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PLM
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PLEPS
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PTKE
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PDTDZ
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PEMOIST
  REAL, DIMENSION(SIZE(PD,1),SIZE(PD,2),SIZE(PD,3)) :: D_M3_THR_W2R_O_DDRDZ
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_THR_W2R_O_DDRDZ',0,ZHOOK_HANDLE)
D_M3_THR_W2R_O_DDRDZ = D_M3_THR_W2TH_O_DDTDZ(CSTURB,KKA,KKU,KKL,PREDR1,PREDTH1,PD,PLM,PLEPS,PTKE,PBLL_O_E,PDTDZ,PEMOIST)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_THR_W2R_O_DDRDZ',1,ZHOOK_HANDLE)
END FUNCTION D_M3_THR_W2R_O_DDRDZ
!----------------------------------------------------------------------------
FUNCTION D_M3_THR_W2R_O_DDTDZ(CSTURB,KKA,KKU,KKL,PREDR1,PREDTH1,PD,PLM,PLEPS,PTKE)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  INTEGER,                INTENT(IN) :: KKA 
  INTEGER,                INTENT(IN) :: KKU  
  INTEGER,                INTENT(IN) :: KKL
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDR1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PD
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PLM
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PLEPS
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PTKE
  REAL, DIMENSION(SIZE(PD,1),SIZE(PD,2),SIZE(PD,3)) :: D_M3_THR_W2R_O_DDTDZ
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_THR_W2R_O_DDTDZ',0,ZHOOK_HANDLE)
D_M3_THR_W2R_O_DDTDZ = D_M3_THR_W2TH_O_DDRDZ(CSTURB,KKA,KKU,KKL,PREDR1,PREDTH1,PD,PLM,PLEPS,PTKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_THR_W2R_O_DDTDZ',1,ZHOOK_HANDLE)
END FUNCTION D_M3_THR_W2R_O_DDTDZ
!----------------------------------------------------------------------------
!
END MODULE MODE_PRANDTL

