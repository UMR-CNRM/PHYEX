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
USE MODD_DIMPHYEX,   ONLY : DIMPHYEX_t
USE MODD_PARAMETERS, ONLY : JPVEXT_TURB
!
USE MODI_SHUMAN, ONLY: MZM, MZF
IMPLICIT NONE
!----------------------------------------------------------------------------
CONTAINS
!----------------------------------------------------------------------------
      SUBROUTINE PRANDTL(D,CST,CSTURB,KRR,KSV,KRRI,OTURB_DIAG,&
                         HTURBDIM,OOCEAN,OHARAT,O2D,OCOMPUTE_SRC,&
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
USE MODD_DIMPHYEX,   ONLY: DIMPHYEX_t
USE MODD_FIELD,          ONLY: TFIELDDATA, TYPEREAL
USE MODD_IO,             ONLY: TFILEDATA
USE MODD_PARAMETERS, ONLY: JPVEXT_TURB
!
USE MODI_GRADIENT_M
USE MODE_EMOIST, ONLY : EMOIST
USE MODE_ETHETA, ONLY : ETHETA
USE MODI_SHUMAN, ONLY: MZM
USE MODE_IO_FIELD_WRITE, ONLY: IO_FIELD_WRITE
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
TYPE(DIMPHYEX_t),       INTENT(IN)   :: D
TYPE(CST_t),            INTENT(IN)   :: CST
TYPE(CSTURB_t),         INTENT(IN)   :: CSTURB
INTEGER,                INTENT(IN)   :: KSV           ! number of scalar variables

INTEGER,                INTENT(IN)   :: KRR           ! number of moist var.
INTEGER,                INTENT(IN)   :: KRRI          ! number of ice var.
!
LOGICAL,                INTENT(IN)   ::  OTURB_DIAG   ! switch to write some
                                 ! diagnostic fields in the syncronous FM-file
LOGICAL,                INTENT(IN)   ::  OOCEAN       ! switch for Ocean model version
LOGICAL,                INTENT(IN)   ::  OHARAT
LOGICAL,                INTENT(IN)   ::  OCOMPUTE_SRC ! flag to define dimensions of SIGS and
LOGICAL, INTENT(IN) :: O2D               ! Logical for 2D model version (modd_conf)
CHARACTER(LEN=4),       INTENT(IN)   ::  HTURBDIM     ! Kind of turbulence param.
TYPE(TFILEDATA),        INTENT(IN)   ::  TPFILE       ! Output file
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PDXX,PDYY,PDZZ,PDZX,PDZY
                                                  ! metric coefficients
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PTHVREF  ! Virtual Potential Temp.
                                                  ! of the reference state
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PLOCPEXNM ! Lv(T)/Cp/Exner at t-1 
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PATHETA      ! coefficients between 
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PAMOIST      ! s and Thetal and Rnp
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PLM      ! Turbulent Mixing length
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PLEPS    ! Dissipative length
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PTHLM,PTKEM! Conservative Potential 
                                                  ! Temperature and TKE at t-1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT,KRR), INTENT(IN) ::  PRM      ! Mixing ratios at  t-1
                                                  ! with PRM(:,:,:,1) = cons.
                                                  ! mixing ratio
REAL, DIMENSION(D%NIT,D%NJT,D%NKT,KSV), INTENT(IN) ::  PSVM     ! Scalars at t-1      
REAL, DIMENSION(MERGE(D%NIT,0,OCOMPUTE_SRC),&
                MERGE(D%NJT,0,OCOMPUTE_SRC),&
                MERGE(D%NKT,0,OCOMPUTE_SRC)),  INTENT(IN)   ::  PSRCM
                                  ! s'r'c/2Sigma_s2 at t-1 multiplied by Lambda_3
!
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT)  ::  PREDTH1 ! Redelsperger number R_theta
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT)  ::  PREDR1  ! Redelsperger number R_q
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT)  ::  PRED2TH3 ! Redelsperger number R*2_theta
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT)  ::  PRED2R3  ! Redelsperger number R*2_q
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT)  ::  PRED2THR3! Redelsperger number R*2_thq
REAL, DIMENSION(D%NIT,D%NJT,D%NKT,KSV), INTENT(OUT)::  PREDS1   ! Redelsperger number R_s
REAL, DIMENSION(D%NIT,D%NJT,D%NKT,KSV), INTENT(OUT)::  PRED2THS3! Redelsperger number R*2_thsv
REAL, DIMENSION(D%NIT,D%NJT,D%NKT,KSV), INTENT(OUT)::  PRED2RS3 ! Redelsperger number R*2_qsv
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT)  ::  PBLL_O_E! beta*Lk*Leps/tke
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT)  ::  PETHETA ! coefficient E_theta
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT)  ::  PEMOIST ! coefficient E_moist
!
!
!       0.2  declaration of local variables
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::  &
                  ZW1, ZW2, ZW3, &
!                                                 working variables
                  ZWORK1,ZWORK2,ZWORK3 
!                                                 working variables for explicit array
!                                                     
INTEGER :: IKB      ! vertical index value for the first inner mass point
INTEGER :: IKE      ! vertical index value for the last inner mass point
INTEGER::  JSV,JI,JJ,JK ! loop index

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
IKB = D%NKB
IKE = D%NKE 
!
PETHETA(:,:,:) = MZM(ETHETA(D,CST,KRR,KRRI,PTHLM,PRM,PLOCPEXNM,PATHETA,PSRCM,OOCEAN,OCOMPUTE_SRC), D%NKA, D%NKU, D%NKL)
PEMOIST(:,:,:) = MZM(EMOIST(D,CST,KRR,KRRI,PTHLM,PRM,PLOCPEXNM,PAMOIST,PSRCM,OOCEAN), D%NKA, D%NKU, D%NKL)
!$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT)
PETHETA(:,:,D%NKA) = 2.*PETHETA(:,:,IKB) - PETHETA(:,:,IKB+D%NKL)
PEMOIST(:,:,D%NKA) = 2.*PEMOIST(:,:,IKB) - PEMOIST(:,:,IKB+D%NKL)
!$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT)
!
!---------------------------------------------------------------------------
IF (.NOT. OHARAT) THEN
!
!          1.3 1D Redelsperger numbers
!
PBLL_O_E(:,:,:) = MZM(CST%XG / PTHVREF(:,:,:) * PLM(:,:,:) * PLEPS(:,:,:) / PTKEM(:,:,:), D%NKA, D%NKU, D%NKL)
ZWORK1 =  GZ_M_W(D%NKA, D%NKU, D%NKL,PTHLM,PDZZ)
IF (KRR /= 0) THEN                ! moist case
  ZWORK2 =  GZ_M_W(D%NKA, D%NKU, D%NKL,PRM(:,:,:,1),PDZZ)
  !$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
  PREDTH1(:,:,:)= CSTURB%XCTV*PBLL_O_E(:,:,:) * PETHETA(:,:,:) * ZWORK1(:,:,:)
  PREDR1(:,:,:) = CSTURB%XCTV*PBLL_O_E(:,:,:) * PEMOIST(:,:,:) * ZWORK2(:,:,:)
  !$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
ELSE                              ! dry case
  !$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
  PREDTH1(:,:,:)= CSTURB%XCTV*PBLL_O_E(:,:,:)  * ZWORK1(:,:,:)
  PREDR1(:,:,:) = 0.
  !$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
END IF
!
!       3. Limits on 1D Redelperger numbers
!          --------------------------------
!
ZMINVAL = (1.-1./CSTURB%XPHI_LIM)
!
DO JK=1,D%NKT 
 DO JJ=1,D%NJT 
  DO JI=1,D%NIT 
   ZW1(JI,JJ,JK) = 1.
   ZW2(JI,JJ,JK) = 1.
   !
   IF (PREDTH1(JI,JJ,JK)+PREDR1(JI,JJ,JK)<-ZMINVAL)THEN
    ZW1(JI,JJ,JK) = (-ZMINVAL) / (PREDTH1(JI,JJ,JK)+PREDR1(JI,JJ,JK))
   ENDIF
   !
   IF (PREDTH1(JI,JJ,JK)<-ZMINVAL)THEN
    ZW2(JI,JJ,JK) = (-ZMINVAL) / (PREDTH1(JI,JJ,JK))
   ENDIF
   ZW2(JI,JJ,JK) = MIN(ZW1(JI,JJ,JK),ZW2(JI,JJ,JK))
   !
   ZW1(JI,JJ,JK) = 1.
   IF (PREDR1(JI,JJ,JK)<-ZMINVAL)THEN
    ZW1(JI,JJ,JK) = (-ZMINVAL) / (PREDR1(JI,JJ,JK))
   ENDIF
   ZW1(JI,JJ,JK) = MIN(ZW2(JI,JJ,JK),ZW1(JI,JJ,JK))
   !
   !
   !       3. Modification of Mixing length and dissipative length
   !          ----------------------------------------------------
   !
   PBLL_O_E(JI,JJ,JK) = PBLL_O_E(JI,JJ,JK) * ZW1(JI,JJ,JK)
   PREDTH1(JI,JJ,JK)  = PREDTH1(JI,JJ,JK)  * ZW1(JI,JJ,JK)
   PREDR1(JI,JJ,JK)   = PREDR1(JI,JJ,JK)   * ZW1(JI,JJ,JK)
   !
   !       4. Threshold for very small (in absolute value) Redelperger numbers
   !          ----------------------------------------------------------------
   !
   IF(PREDTH1(JI,JJ,JK) < 0.) THEN
    ZW2(JI,JJ,JK)=-1.
   ELSE
    ZW2(JI,JJ,JK)=1.
   END IF
   PREDTH1(JI,JJ,JK)= ZW2(JI,JJ,JK) * MAX(1.E-30, ZW2(JI,JJ,JK)*PREDTH1(JI,JJ,JK))
   !
   IF (KRR /= 0) THEN                ! moist case
    IF(PREDR1(JI,JJ,JK) < 0.) THEN
     ZW2(JI,JJ,JK)=-1.
    ELSE
     ZW2(JI,JJ,JK)=1.
    END IF
    PREDR1(JI,JJ,JK)= ZW2(JI,JJ,JK) * MAX(1.E-30, ZW2(JI,JJ,JK)*PREDR1(JI,JJ,JK))
   END IF
  ENDDO
 ENDDO
ENDDO
!
!
!---------------------------------------------------------------------------
!
!          For the scalar variables
DO JSV=1,KSV
  ZWORK1 = GZ_M_W(D%NKA, D%NKU, D%NKL,PSVM(:,:,:,JSV),PDZZ)
  !$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
  PREDS1(:,:,:,JSV)=CSTURB%XCTV*PBLL_O_E(:,:,:)*ZWORK1(:,:,:)
  !$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
END DO
!
DO JSV=1,KSV
 DO JK=1,D%NKT 
  DO JJ=1,D%NJT 
   DO JI=1,D%NIT 
    IF(PREDS1(JI,JJ,JK,JSV) < 0.) THEN
     ZW2(JI,JJ,JK)=-1.
    ELSE
     ZW2(JI,JJ,JK)=1.
    END IF
    PREDS1(JI,JJ,JK,JSV)= ZW2(JI,JJ,JK) * MAX(1.E-30, ZW2(JI,JJ,JK)*PREDS1(JI,JJ,JK,JSV))
   ENDDO
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
  !$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
  PRED2TH3(:,:,:)  = PREDTH1(:,:,:)**2
!
  PRED2R3(:,:,:)   = PREDR1(:,:,:) **2
!
  PRED2THR3(:,:,:) = PREDTH1(:,:,:) * PREDR1(:,:,:)
  !$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
!
ELSE IF (O2D) THEN                      ! 3D case in a 2D model
!
    ZWORK1 = MZM(GX_M_M(PTHLM,PDXX,PDZZ,PDZX, D%NKA, D%NKU, D%NKL)**2, D%NKA, D%NKU, D%NKL)
  IF (KRR /= 0) THEN                 ! moist 3D case
    ZWORK2 = MZM(GX_M_M(PRM(:,:,:,1),PDXX,PDZZ,PDZX, D%NKA, D%NKU, D%NKL)**2, D%NKA, D%NKU, D%NKL)
    ZWORK3 = MZM(GX_M_M(PRM(:,:,:,1),PDXX,PDZZ,PDZX, D%NKA, D%NKU, D%NKL)*     &
                 GX_M_M(PTHLM,PDXX,PDZZ,PDZX, D%NKA, D%NKU, D%NKL), D%NKA, D%NKU, D%NKL)
!
    !$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
    PRED2TH3(:,:,:)= PREDTH1(:,:,:)**2+(CSTURB%XCTV*PBLL_O_E(:,:,:)*PETHETA(:,:,:) )**2 * &
                     ZWORK1(:,:,:)
!
    PRED2R3(:,:,:)= PREDR1(:,:,:)**2 + (CSTURB%XCTV*PBLL_O_E(:,:,:)*PEMOIST(:,:,:))**2 * &
                    ZWORK2(:,:,:)
!
    PRED2THR3(:,:,:)= PREDR1(:,:,:) * PREDTH1(:,:,:) +  CSTURB%XCTV**2*PBLL_O_E(:,:,:)**2 *   &
                      PEMOIST(:,:,:) * PETHETA(:,:,:) *                         &
                      ZWORK3(:,:,:)
    !$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
!
    PRED2TH3(:,:,IKB)=PRED2TH3(:,:,IKB+D%NKL)
    PRED2R3(:,:,IKB)=PRED2R3(:,:,IKB+D%NKL) 
    PRED2THR3(:,:,IKB)=PRED2THR3(:,:,IKB+D%NKL)
!
  ELSE                 ! dry 3D case in a 2D model
    !$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
    PRED2TH3(:,:,:) = PREDTH1(:,:,:)**2 +  CSTURB%XCTV**2*PBLL_O_E(:,:,:)**2 *     &
                      ZWORK1(:,:,:)      
    !$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
    PRED2TH3(:,:,IKB)=PRED2TH3(:,:,IKB+D%NKL)
!
    PRED2R3(:,:,:) = 0.
!
    PRED2THR3(:,:,:) = 0.
!
  END IF
!
ELSE                                 ! 3D case in a 3D model
!
  ZWORK1 = MZM(GX_M_M(PTHLM,PDXX,PDZZ,PDZX, D%NKA, D%NKU, D%NKL)**2 &
      + GY_M_M(PTHLM,PDYY,PDZZ,PDZY, D%NKA, D%NKU, D%NKL)**2, D%NKA, D%NKU, D%NKL)

  IF (KRR /= 0) THEN                 ! moist 3D case
    ZWORK2 = MZM(GX_M_M(PRM(:,:,:,1),PDXX,PDZZ,PDZX, D%NKA, D%NKU, D%NKL)**2 + &
        GY_M_M(PRM(:,:,:,1),PDYY,PDZZ,PDZY, D%NKA, D%NKU, D%NKL)**2, D%NKA, D%NKU, D%NKL)
    ZWORK3 = MZM(GX_M_M(PRM(:,:,:,1),PDXX,PDZZ,PDZX, D%NKA, D%NKU, D%NKL)*   &
         GX_M_M(PTHLM,PDXX,PDZZ,PDZX, D%NKA, D%NKU, D%NKL)+                           &
         GY_M_M(PRM(:,:,:,1),PDYY,PDZZ,PDZY, D%NKA, D%NKU, D%NKL)*                    &
         GY_M_M(PTHLM,PDYY,PDZZ,PDZY, D%NKA, D%NKU, D%NKL), D%NKA, D%NKU, D%NKL)

    !$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
    PRED2TH3(:,:,:)= PREDTH1(:,:,:)**2 +  ( CSTURB%XCTV*PBLL_O_E(:,:,:)*PETHETA(:,:,:) )**2 * &
                     ZWORK1(:,:,:)      
!
    PRED2R3(:,:,:)= PREDR1(:,:,:)**2 + (CSTURB%XCTV*PBLL_O_E(:,:,:)*PEMOIST(:,:,:))**2 *      &
                    ZWORK2(:,:,:)
!
    PRED2THR3(:,:,:)= PREDR1(:,:,:) * PREDTH1(:,:,:) +  CSTURB%XCTV**2*PBLL_O_E(:,:,:)**2 *   &
         PEMOIST(:,:,:) * PETHETA(:,:,:) * ZWORK3(:,:,:)

    !$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
!
    PRED2TH3(:,:,IKB)=PRED2TH3(:,:,IKB+D%NKL)
    PRED2R3(:,:,IKB)=PRED2R3(:,:,IKB+D%NKL)
    PRED2THR3(:,:,IKB)=PRED2THR3(:,:,IKB+D%NKL)
!
  ELSE                 ! dry 3D case in a 3D model
    !$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
    PRED2TH3(:,:,:) = PREDTH1(:,:,:)**2 +  CSTURB%XCTV**2*PBLL_O_E(:,:,:)**2 *                &
                      ZWORK1(:,:,:)
    !$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
! 
    PRED2TH3(:,:,IKB)=PRED2TH3(:,:,IKB+D%NKL)
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
DO JSV=1,KSV
!
  IF(HTURBDIM=='1DIM') THEN
!        1D case
    !$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
    PRED2THS3(:,:,:,JSV)  = PREDS1(:,:,:,JSV) * PREDTH1(:,:,:)
    IF (KRR /= 0) THEN
      PRED2RS3(:,:,:,JSV)   = PREDR1(:,:,:) *PREDS1(:,:,:,JSV)
    ELSE
      PRED2RS3(:,:,:,JSV)   = 0.
    END IF
    !$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
!
  ELSE  IF (O2D) THEN ! 3D case in a 2D model
!
    ZW1 = MZM((CST%XG / PTHVREF * PLM * PLEPS / PTKEM)**2, D%NKA, D%NKU, D%NKL)
    ZWORK2 = MZM(GX_M_M(PSVM(:,:,:,JSV),PDXX,PDZZ,PDZX, D%NKA, D%NKU, D%NKL)*       &
                           GX_M_M(PTHLM,PDXX,PDZZ,PDZX, D%NKA, D%NKU, D%NKL),                 &
                           D%NKA, D%NKU, D%NKL)
    ZWORK3 =  MZM(GX_M_M(PSVM(:,:,:,JSV),PDXX,PDZZ,PDZX, D%NKA, D%NKU, D%NKL)*       &
                           GX_M_M(PRM(:,:,:,1),PDXX,PDZZ,PDZX, D%NKA, D%NKU, D%NKL),          &
                           D%NKA, D%NKU, D%NKL)
!
    !$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)    
    IF (KRR /= 0) THEN
      ZWORK1(:,:,:) = ZW1(:,:,:)*PETHETA(:,:,:)
    ELSE
      ZWORK1(:,:,:) = ZW1(:,:,:)
    END IF
    PRED2THS3(:,:,:,JSV) = PREDTH1(:,:,:) * PREDS1(:,:,:,JSV)   +        &
                       ZWORK1(:,:,:) * ZWORK2(:,:,:)
                       !
    IF (KRR /= 0) THEN
      PRED2RS3(:,:,:,JSV) = PREDR1(:,:,:) * PREDS1(:,:,:,JSV)   +        &
                       ZW1(:,:,:) * PEMOIST(:,:,:) * ZWORK3(:,:,:)
    ELSE
      PRED2RS3(:,:,:,JSV) = 0.
    END IF
    !$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
!
  ELSE ! 3D case in a 3D model
!
    ZW1 = MZM((CST%XG / PTHVREF * PLM * PLEPS / PTKEM)**2, D%NKA, D%NKU, D%NKL)
    ZWORK2 =  MZM(GX_M_M(PSVM(:,:,:,JSV),PDXX,PDZZ,PDZX, D%NKA, D%NKU, D%NKL)*       &
                           GX_M_M(PTHLM,PDXX,PDZZ,PDZX, D%NKA, D%NKU, D%NKL)                  &
                          +GY_M_M(PSVM(:,:,:,JSV),PDYY,PDZZ,PDZY, D%NKA, D%NKU, D%NKL)*       &
                           GY_M_M(PTHLM,PDYY,PDZZ,PDZY, D%NKA, D%NKU, D%NKL),                 &
                           D%NKA, D%NKU, D%NKL)
    ZWORK3 =  MZM(GX_M_M(PSVM(:,:,:,JSV),PDXX,PDZZ,PDZX, D%NKA, D%NKU, D%NKL)*       &
                           GX_M_M(PRM(:,:,:,1),PDXX,PDZZ,PDZX, D%NKA, D%NKU, D%NKL)           &
                          +GY_M_M(PSVM(:,:,:,JSV),PDYY,PDZZ,PDZY, D%NKA, D%NKU, D%NKL)*       &
                           GY_M_M(PRM(:,:,:,1),PDYY,PDZZ,PDZY, D%NKA, D%NKU, D%NKL),          &
                           D%NKA, D%NKU, D%NKL)

    !$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)    
    IF (KRR /= 0) THEN
      ZWORK1(:,:,:) = ZW1(:,:,:)*PETHETA(:,:,:)
    ELSE
      ZWORK1(:,:,:) = ZW1(:,:,:)
    END IF
    PRED2THS3(:,:,:,JSV) = PREDTH1(:,:,:) * PREDS1(:,:,:,JSV)   +        &
                       ZWORK1(:,:,:)*ZWORK2(:,:,:)
                      !
    IF (KRR /= 0) THEN
      PRED2RS3(:,:,:,JSV) = PREDR1(:,:,:) * PREDS1(:,:,:,JSV)   +        &
                       ZW1(:,:,:) * PEMOIST(:,:,:) * ZWORK3(:,:,:)
    ELSE
      PRED2RS3(:,:,:,JSV) = 0.
    END IF
    !$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)    
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
SUBROUTINE SMOOTH_TURB_FUNCT(D,CSTURB,PPHI3,PF_LIM,PF)
!
TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
TYPE(DIMPHYEX_t),       INTENT(IN)   :: D
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)    :: PPHI3   ! Phi3
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)    :: PF_LIM  ! Value of F when Phi3 is
!                                                ! larger than Phi_lim
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(INOUT) :: PF      ! function F to smooth
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: ZCOEF
INTEGER :: JI,JJ,JK
!
!* adds a artificial correction to smooth the function near the discontinuity
!  point at Phi3 = Phi_lim
!  This smoothing is applied between 0.9*phi_lim (=2.7) and Phi_lim (=3)
!   Note that in the Boundary layer, phi is usually between 0.8 and 1
!
!
!$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)    
ZCOEF(:,:,:) = MAX(MIN((  10.*(1.-PPHI3(:,:,:)/CSTURB%XPHI_LIM)) ,1.), 0.) 
!
PF(:,:,:) =     ZCOEF(:,:,:)   * PF(:,:,:)    &
          + (1.-ZCOEF(:,:,:))  * PF_LIM(:,:,:)
!$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)    
!
END SUBROUTINE SMOOTH_TURB_FUNCT
!----------------------------------------------------------------------------
FUNCTION PHI3(D,CSTURB,PREDTH1,PREDR1,PRED2TH3,PRED2R3,PRED2THR3,HTURBDIM,OUSERV)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  TYPE(DIMPHYEX_t),                   INTENT(IN)   :: D
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PRED2TH3
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PRED2R3
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PRED2THR3
  CHARACTER(LEN=4),       INTENT(IN) :: HTURBDIM  ! 1DIM or 3DIM turb. scheme
  LOGICAL,                INTENT(IN) :: OUSERV    ! flag to use vapor
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: PHI3
!
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: ZW1, ZW2
  INTEGER :: IKB, IKE, JI,JJ,JK
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:PHI3',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PREDTH1,3)-JPVEXT_TURB
!
IF (HTURBDIM=='3DIM') THEN
        !* 3DIM case
  !$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)    
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
  !$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)    

  !$mnh_expand_where(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)    
  WHERE( PHI3(:,:,:) <= 0. .OR. PHI3(:,:,:) > CSTURB%XPHI_LIM )
    PHI3(:,:,:) = CSTURB%XPHI_LIM
  END WHERE
  !$mnh_end_expand_where(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)    
ELSE
        !* 1DIM case
  !$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)    
  IF (OUSERV) THEN
    PHI3(:,:,:)= 1./(1.+PREDTH1(:,:,:)+PREDR1(:,:,:))
  ELSE
    PHI3(:,:,:)= 1./(1.+PREDTH1(:,:,:))
  END IF
  !$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)    
END IF
!
PHI3(:,:,IKB-1)=PHI3(:,:,IKB)
PHI3(:,:,IKE+1)=PHI3(:,:,IKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:PHI3',1,ZHOOK_HANDLE)
END FUNCTION PHI3
!----------------------------------------------------------------------------
FUNCTION PSI_SV(D,CSTURB,KSV,PREDTH1,PREDR1,PREDS1,PRED2THS,PRED2RS,PPHI3,PPSI3)
  TYPE(CSTURB_t),                  INTENT(IN)      :: CSTURB
  TYPE(DIMPHYEX_t),                INTENT(IN)      :: D
  INTEGER,                         INTENT(IN)      :: KSV
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT,KSV), INTENT(IN) :: PREDS1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT,KSV), INTENT(IN) :: PRED2THS
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT,KSV), INTENT(IN) :: PRED2RS
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN) :: PPHI3
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN) :: PPSI3
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT,KSV) :: PSI_SV
!
  INTEGER :: IKB, IKE
  INTEGER :: JSV,JI,JJ,JK
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:PSI_SV',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PREDTH1,3)-JPVEXT_TURB
!
DO JSV=1,SIZE(PSI_SV,4)
  !$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)    
  PSI_SV(:,:,:,JSV) = ( 1.                                             &
    - (CSTURB%XCPR3+CSTURB%XCPR5) * (PRED2THS(:,:,:,JSV)/PREDS1(:,:,:,JSV)-PREDTH1(:,:,:)) &
    - (CSTURB%XCPR4+CSTURB%XCPR5) * (PRED2RS(:,:,:,JSV)/PREDS1(:,:,:,JSV)-PREDR1(:,:,:)) &
    - CSTURB%XCPR3 * PREDTH1(:,:,:) * PPHI3(:,:,:) - CSTURB%XCPR4 * PREDR1(:,:,:) * PPSI3(:,:,:)                 &
                ) / (  1. + CSTURB%XCPR5 * ( PREDTH1(:,:,:) + PREDR1(:,:,:) ) )           
  
  !$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)    
!        control of the PSI_SV positivity
  !$mnh_expand_where(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)    
  WHERE ( (PSI_SV(:,:,:,JSV) <=0.).AND. (PREDTH1(:,:,:)+PREDR1(:,:,:)) <= 0. )
    PSI_SV(:,:,:,JSV)=CSTURB%XPHI_LIM
  END WHERE
  !$mnh_end_expand_where(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)    

  !$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)    
  PSI_SV(:,:,:,JSV) = MAX( 1.E-4, MIN(CSTURB%XPHI_LIM,PSI_SV(:,:,:,JSV)) )
  !$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)    
!
  PSI_SV(:,:,IKB-1,JSV)=PSI_SV(:,:,IKB,JSV)
  PSI_SV(:,:,IKE+1,JSV)=PSI_SV(:,:,IKE,JSV)
END DO
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:PSI_SV',1,ZHOOK_HANDLE)
END FUNCTION PSI_SV
!----------------------------------------------------------------------------
FUNCTION D_PHI3DTDZ_O_DDTDZ(D,CSTURB,PPHI3,PREDTH1,PREDR1,PRED2TH3,PRED2THR3,HTURBDIM,OUSERV)
  TYPE(CSTURB_t),                  INTENT(IN)    :: CSTURB
  TYPE(DIMPHYEX_t),                   INTENT(IN)   :: D
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PPHI3
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PRED2TH3
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PRED2THR3
  CHARACTER(LEN=4),       INTENT(IN) :: HTURBDIM  ! 1DIM or 3DIM turb. scheme
  LOGICAL,                INTENT(IN) :: OUSERV    ! flag to use vapor
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: D_PHI3DTDZ_O_DDTDZ
  INTEGER :: IKB, IKE,JL,JK,JJ,JI
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_PHI3DTDZ_O_DDTDZ',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PREDTH1,3)-JPVEXT_TURB
!
IF (HTURBDIM=='3DIM') THEN
        !* 3DIM case
  IF (OUSERV) THEN
   !$mnh_expand_where(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)    
#ifdef REPRO48
    WHERE (PPHI3(:,:,:)/=CSTURB%XPHI_LIM)
#else
    WHERE (PPHI3(:,:,:)<=CSTURB%XPHI_LIM)
#endif
      D_PHI3DTDZ_O_DDTDZ(:,:,:) = PPHI3(:,:,:)                       &
          * (1. - PREDTH1(:,:,:) * (3./2.+PREDTH1(:,:,:)+PREDR1(:,:,:))          &
               /((1.+PREDTH1(:,:,:)+PREDR1(:,:,:))*(1.+1./2.*(PREDTH1(:,:,:)+PREDR1(:,:,:))))) &
          + (1.+PREDR1(:,:,:))*(PRED2THR3(:,:,:)+PRED2TH3(:,:,:))                       &
               / (PREDTH1(:,:,:)*(1.+PREDTH1(:,:,:)+PREDR1(:,:,:))* &
                 (1.+1./2.*(PREDTH1(:,:,:)+PREDR1(:,:,:)))) &
          - (1./2.*PREDTH1(:,:,:)+PREDR1(:,:,:) * (1.+PREDTH1(:,:,:)+PREDR1(:,:,:)))           &
               / ((1.+PREDTH1(:,:,:)+PREDR1(:,:,:))*(1.+1./2.*(PREDTH1(:,:,:)+PREDR1(:,:,:))))
    ELSEWHERE
      D_PHI3DTDZ_O_DDTDZ(:,:,:) = PPHI3(:,:,:)
    ENDWHERE
   !$mnh_end_expand_where(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)    

!
  ELSE
   !$mnh_expand_where(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)    
#ifdef REPRO48
    WHERE (PPHI3(:,:,:)/=CSTURB%XPHI_LIM)
#else
    WHERE (PPHI3(:,:,:)<=CSTURB%XPHI_LIM)
#endif
    D_PHI3DTDZ_O_DDTDZ(:,:,:) = PPHI3(:,:,:)             &
          * (1. - PREDTH1(:,:,:) * (3./2.+PREDTH1(:,:,:))      &
               /((1.+PREDTH1(:,:,:))*(1.+1./2.*PREDTH1(:,:,:))))        &
          + PRED2TH3(:,:,:) / (PREDTH1(:,:,:)*(1.+PREDTH1(:,:,:))*(1.+1./2.*PREDTH1(:,:,:))) &
          - 1./2.*PREDTH1(:,:,:) / ((1.+PREDTH1(:,:,:))*(1.+1./2.*PREDTH1(:,:,:)))
    ELSEWHERE
      D_PHI3DTDZ_O_DDTDZ(:,:,:) = PPHI3(:,:,:)
    ENDWHERE
   !$mnh_end_expand_where(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)    
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
CALL SMOOTH_TURB_FUNCT(D,CSTURB,PPHI3,PPHI3,D_PHI3DTDZ_O_DDTDZ)
#endif
!
D_PHI3DTDZ_O_DDTDZ(:,:,IKB-1)=D_PHI3DTDZ_O_DDTDZ(:,:,IKB)
D_PHI3DTDZ_O_DDTDZ(:,:,IKE+1)=D_PHI3DTDZ_O_DDTDZ(:,:,IKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_PHI3DTDZ_O_DDTDZ',1,ZHOOK_HANDLE)
END FUNCTION D_PHI3DTDZ_O_DDTDZ
!----------------------------------------------------------------------------
FUNCTION D_PHI3DRDZ_O_DDRDZ(D,CSTURB,PPHI3,PREDTH1,PREDR1,PRED2TH3,PRED2THR3,HTURBDIM,OUSERV)
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PPHI3
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PRED2TH3
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PRED2THR3
  CHARACTER(LEN=4),       INTENT(IN) :: HTURBDIM  ! 1DIM or 3DIM turb. scheme
  LOGICAL,                INTENT(IN) :: OUSERV    ! flag to use vapor
   REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: D_PHI3DRDZ_O_DDRDZ
  INTEGER :: IKB, IKE, JI,JJ,JK
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
   !$mnh_expand_where(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)    
#ifdef REPRO48
    WHERE (PPHI3(:,:,:)/=CSTURB%XPHI_LIM)
#else
    WHERE (PPHI3(:,:,:)<=CSTURB%XPHI_LIM)
#endif
      D_PHI3DRDZ_O_DDRDZ(:,:,:) =          &
                PPHI3(:,:,:) * (1.-PREDR1(:,:,:)*(3./2.+PREDTH1(:,:,:)+PREDR1(:,:,:)) &
                  / ((1.+PREDTH1(:,:,:)+PREDR1(:,:,:))*(1.+1./2.*(PREDTH1(:,:,:)+PREDR1(:,:,:)))))  &
              - PREDR1(:,:,:) * (PRED2THR3(:,:,:)+PRED2TH3(:,:,:)) / (PREDTH1(:,:,:)         &
                  * (1.+PREDTH1(:,:,:)+PREDR1(:,:,:))*(1.+1./2.*(PREDTH1(:,:,:)+PREDR1(:,:,:))))    &
              + PREDR1(:,:,:) * (1./2.+PREDTH1(:,:,:)+PREDR1(:,:,:))                  &
                  / ((1.+PREDTH1(:,:,:)+PREDR1(:,:,:))*(1.+1./2.*(PREDTH1(:,:,:)+PREDR1(:,:,:))))
    ELSEWHERE
      D_PHI3DRDZ_O_DDRDZ(:,:,:) = PPHI3(:,:,:)
    END WHERE
   !$mnh_end_expand_where(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)    
  ELSE
    D_PHI3DRDZ_O_DDRDZ(:,:,:) = PPHI3(:,:,:)
  END IF
ELSE
        !* 1DIM case
    !$mnh_expand_where(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)    
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
  !$mnh_end_expand_where(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)    
END IF
!
#ifdef REPRO48
#else
!* smoothing
CALL SMOOTH_TURB_FUNCT(D,CSTURB,PPHI3,PPHI3,D_PHI3DRDZ_O_DDRDZ)
#endif
!
D_PHI3DRDZ_O_DDRDZ(:,:,IKB-1)=D_PHI3DRDZ_O_DDRDZ(:,:,IKB)
D_PHI3DRDZ_O_DDRDZ(:,:,IKE+1)=D_PHI3DRDZ_O_DDRDZ(:,:,IKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_PHI3DRDZ_O_DDRDZ',1,ZHOOK_HANDLE)
END FUNCTION D_PHI3DRDZ_O_DDRDZ
!----------------------------------------------------------------------------
FUNCTION D_PHI3DTDZ2_O_DDTDZ(D,CSTURB,PPHI3,PREDTH1,PREDR1,PRED2TH3,PRED2THR3,PDTDZ,HTURBDIM,OUSERV)
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PPHI3
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PRED2TH3
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PRED2THR3
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PDTDZ
  CHARACTER(LEN=4),       INTENT(IN) :: HTURBDIM  ! 1DIM or 3DIM turb. scheme
  LOGICAL,                INTENT(IN) :: OUSERV    ! flag to use vapor
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: D_PHI3DTDZ2_O_DDTDZ
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: ZWORK1 ! working array
  INTEGER :: IKB, IKE, JI,JJ,JK
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_PHI3DTDZ2_O_DDTDZ',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PREDTH1,3)-JPVEXT_TURB
!
!
IF (HTURBDIM=='3DIM') THEN
   ! by derivation of (phi3 dtdz) * dtdz according to dtdz we obtain:
   ZWORK1 = D_PHI3DTDZ_O_DDTDZ(D,CSTURB,PPHI3,PREDTH1,PREDR1,PRED2TH3,PRED2THR3,HTURBDIM,OUSERV) 
   !$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
   D_PHI3DTDZ2_O_DDTDZ(:,:,:) = PDTDZ(:,:,:) * (PPHI3(:,:,:) +  ZWORK1(:,:,:))
   !$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
ELSE
        !* 1DIM case
    !$mnh_expand_where(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)    
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
    !$mnh_end_expand_where(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)    
END IF
!
#ifdef REPRO48
#else
!* smoothing
CALL SMOOTH_TURB_FUNCT(D,CSTURB,PPHI3,PPHI3*2.*PDTDZ,D_PHI3DTDZ2_O_DDTDZ)
#endif
!
!
D_PHI3DTDZ2_O_DDTDZ(:,:,IKB-1)=D_PHI3DTDZ2_O_DDTDZ(:,:,IKB)
D_PHI3DTDZ2_O_DDTDZ(:,:,IKE+1)=D_PHI3DTDZ2_O_DDTDZ(:,:,IKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_PHI3DTDZ2_O_DDTDZ',1,ZHOOK_HANDLE)
END FUNCTION D_PHI3DTDZ2_O_DDTDZ
!----------------------------------------------------------------------------
FUNCTION M3_WTH_WTH2(D,CSTURB,PREDTH1,PREDR1,PD,PBLL_O_E,PETHETA)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PETHETA
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: M3_WTH_WTH2
  INTEGER :: IKB, IKE, JI,JJ,JK
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_WTH_WTH2',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PD,3)-JPVEXT_TURB

!$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
M3_WTH_WTH2(:,:,:) = CSTURB%XCSHF*PBLL_O_E(:,:,:)*PETHETA(:,:,:)*0.5/CSTURB%XCTD        &
                   * (1.+0.5*PREDTH1(:,:,:)+PREDR1(:,:,:)) / PD(:,:,:)
!$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
M3_WTH_WTH2(:,:,IKB-1)=M3_WTH_WTH2(:,:,IKB)
M3_WTH_WTH2(:,:,IKE+1)=M3_WTH_WTH2(:,:,IKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_WTH_WTH2',1,ZHOOK_HANDLE)
END FUNCTION M3_WTH_WTH2
!----------------------------------------------------------------------------
FUNCTION D_M3_WTH_WTH2_O_DDTDZ(D,CSTURB,PM3_WTH_WTH2,PREDTH1,PREDR1,PD,PBLL_O_E,PETHETA)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PM3_WTH_WTH2
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PETHETA
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: D_M3_WTH_WTH2_O_DDTDZ
  INTEGER :: IKB, IKE, JI,JJ,JK
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_WTH_WTH2_O_DDTDZ',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PD,3)-JPVEXT_TURB

!$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
D_M3_WTH_WTH2_O_DDTDZ(:,:,:) = (  0.5*CSTURB%XCSHF*PBLL_O_E(:,:,:)*PETHETA(:,:,:)*0.5/CSTURB%XCTD/PD(:,:,:) &
                                - PM3_WTH_WTH2(:,:,:)/PD(:,:,:)*(1.5+PREDTH1(:,:,:)+PREDR1(:,:,:))  )&
                             * PBLL_O_E(:,:,:) * PETHETA(:,:,:) * CSTURB%XCTV
!$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
!
D_M3_WTH_WTH2_O_DDTDZ(:,:,IKB-1)=D_M3_WTH_WTH2_O_DDTDZ(:,:,IKB)
D_M3_WTH_WTH2_O_DDTDZ(:,:,IKE+1)=D_M3_WTH_WTH2_O_DDTDZ(:,:,IKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_WTH_WTH2_O_DDTDZ',1,ZHOOK_HANDLE)
END FUNCTION D_M3_WTH_WTH2_O_DDTDZ
!----------------------------------------------------------------------------
FUNCTION M3_WTH_W2TH(D,CSTURB,PREDTH1,PREDR1,PD,PKEFF,PTKE)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PKEFF
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PTKE
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: M3_WTH_W2TH
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: ZWORK1 ! working array
  INTEGER :: IKB, IKE, JI,JJ,JK
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_WTH_W2TH',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PD,3)-JPVEXT_TURB
ZWORK1 = MZM(PTKE, D%NKA, D%NKU, D%NKL)
!$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
M3_WTH_W2TH(:,:,:) = CSTURB%XCSHF*PKEFF(:,:,:)*1.5/ZWORK1(:,:,:)              &
  * (1. - 0.5*PREDR1(:,:,:)*(1.+PREDR1(:,:,:))/PD(:,:,:) ) / (1.+PREDTH1(:,:,:))
!$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
!
M3_WTH_W2TH(:,:,IKB-1)=M3_WTH_W2TH(:,:,IKB)
M3_WTH_W2TH(:,:,IKE+1)=M3_WTH_W2TH(:,:,IKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_WTH_W2TH',1,ZHOOK_HANDLE)
END FUNCTION M3_WTH_W2TH
!----------------------------------------------------------------------------
FUNCTION D_M3_WTH_W2TH_O_DDTDZ(D,CSTURB,PREDTH1,PREDR1,PD,PBLL_O_E,PETHETA,PKEFF,PTKE)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PETHETA
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PKEFF
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PTKE
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: D_M3_WTH_W2TH_O_DDTDZ
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: ZWORK1 ! working array
  INTEGER :: IKB, IKE, JI,JJ,JK
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_WTH_W2TH_O_DDTDZ',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PD,3)-JPVEXT_TURB

ZWORK1 = MZM(PTKE, D%NKA, D%NKU, D%NKL)
!$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
D_M3_WTH_W2TH_O_DDTDZ(:,:,:) = &
 - CSTURB%XCSHF*PKEFF(:,:,:)*1.5/ZWORK1(:,:,:)/(1.+PREDTH1(:,:,:))**2 &
 * CSTURB%XCTV*PBLL_O_E(:,:,:)*PETHETA(:,:,:)  &
 * (1. - 0.5*PREDR1(:,:,:)*(1.+PREDR1(:,:,:))/PD(:,:,:)* &
   ( 1.+(1.+PREDTH1(:,:,:))*(1.5+PREDR1(:,:,:)+PREDTH1(:,:,:))/PD(:,:,:)) )
!$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
!
D_M3_WTH_W2TH_O_DDTDZ(:,:,IKB-1)=D_M3_WTH_W2TH_O_DDTDZ(:,:,IKB)
D_M3_WTH_W2TH_O_DDTDZ(:,:,IKE+1)=D_M3_WTH_W2TH_O_DDTDZ(:,:,IKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_WTH_W2TH_O_DDTDZ',1,ZHOOK_HANDLE)
END FUNCTION D_M3_WTH_W2TH_O_DDTDZ
!----------------------------------------------------------------------------
FUNCTION M3_WTH_W2R(D,CSTURB,PD,PKEFF,PTKE,PBLL_O_E,PEMOIST,PDTDZ)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PKEFF
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PTKE
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PEMOIST
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PDTDZ
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: M3_WTH_W2R
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: ZWORK1 ! working array
  INTEGER :: IKB, IKE, JI,JJ,JK
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_WTH_W2R',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PD,3)-JPVEXT_TURB

ZWORK1 = MZM(PTKE, D%NKA, D%NKU, D%NKL)
!$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
M3_WTH_W2R(:,:,:) = - CSTURB%XCSHF*PKEFF(:,:,:)*0.75*CSTURB%XCTV*PBLL_O_E(:,:,:) &
                    /ZWORK1(:,:,:)*PEMOIST(:,:,:)*PDTDZ(:,:,:)/PD(:,:,:)
!$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
!
M3_WTH_W2R(:,:,IKB-1)=M3_WTH_W2R(:,:,IKB)
M3_WTH_W2R(:,:,IKE+1)=M3_WTH_W2R(:,:,IKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_WTH_W2R',1,ZHOOK_HANDLE)
END FUNCTION M3_WTH_W2R
!----------------------------------------------------------------------------
FUNCTION D_M3_WTH_W2R_O_DDTDZ(D,CSTURB,PREDTH1,PREDR1,PD,PKEFF,PTKE,PBLL_O_E,PEMOIST)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PKEFF
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PTKE
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PEMOIST
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: D_M3_WTH_W2R_O_DDTDZ
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: ZWORK1 ! working array
  INTEGER :: IKB, IKE, JI,JJ,JK
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_WTH_W2R_O_DDTDZ',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PD,3)-JPVEXT_TURB

ZWORK1 = MZM(PTKE, D%NKA, D%NKU, D%NKL)
!$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
D_M3_WTH_W2R_O_DDTDZ(:,:,:) = - CSTURB%XCSHF*PKEFF(:,:,:)*0.75*CSTURB%XCTV*PBLL_O_E(:,:,:) &
                               /ZWORK1(:,:,:)*PEMOIST(:,:,:)/PD(:,:,:) &
                                     * (1. -  PREDTH1(:,:,:)*(1.5+PREDTH1(:,:,:)+PREDR1(:,:,:))/PD(:,:,:))
!$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
!
D_M3_WTH_W2R_O_DDTDZ(:,:,IKB-1)=D_M3_WTH_W2R_O_DDTDZ(:,:,IKB)
D_M3_WTH_W2R_O_DDTDZ(:,:,IKE+1)=D_M3_WTH_W2R_O_DDTDZ(:,:,IKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_WTH_W2R_O_DDTDZ',1,ZHOOK_HANDLE)
END FUNCTION D_M3_WTH_W2R_O_DDTDZ
!----------------------------------------------------------------------------
FUNCTION M3_WTH_WR2(D,CSTURB,PD,PKEFF,PTKE,PSQRT_TKE,PBLL_O_E,PBETA,PLEPS,PEMOIST,PDTDZ)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PKEFF
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PTKE
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PBETA
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PEMOIST
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PDTDZ
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: M3_WTH_WR2
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: ZWORK1 ! working array
  INTEGER :: IKB, IKE, JI,JJ,JK
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_WTH_WR2',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PD,3)-JPVEXT_TURB

ZWORK1 = MZM(PBETA*PLEPS/(PSQRT_TKE*PTKE), D%NKA, D%NKU, D%NKL)
!$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
M3_WTH_WR2(:,:,:) = - CSTURB%XCSHF*PKEFF(:,:,:)*0.25*PBLL_O_E(:,:,:)*CSTURB%XCTV*PEMOIST(:,:,:)**2       &
                           *ZWORK1(:,:,:)/CSTURB%XCTD*PDTDZ(:,:,:)/PD(:,:,:)
!$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
!
M3_WTH_WR2(:,:,IKB-1)=M3_WTH_WR2(:,:,IKB)
M3_WTH_WR2(:,:,IKE+1)=M3_WTH_WR2(:,:,IKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_WTH_WR2',1,ZHOOK_HANDLE)
END FUNCTION M3_WTH_WR2
!----------------------------------------------------------------------------
FUNCTION D_M3_WTH_WR2_O_DDTDZ(D,CSTURB,PREDTH1,PREDR1,PD,PKEFF,PTKE,PSQRT_TKE,PBLL_O_E,PBETA,PLEPS,PEMOIST)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PKEFF
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PTKE
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PBETA
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PEMOIST
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: D_M3_WTH_WR2_O_DDTDZ
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: ZWORK1 ! working array
  INTEGER :: IKB, IKE, JI,JJ,JK
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_WTH_WR2_O_DDTDZ',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PD,3)-JPVEXT_TURB

ZWORK1 = MZM(PBETA*PLEPS/(PSQRT_TKE*PTKE), D%NKA, D%NKU, D%NKL)
!$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
D_M3_WTH_WR2_O_DDTDZ(:,:,:) = - CSTURB%XCSHF*PKEFF(:,:,:)*0.25*PBLL_O_E(:,:,:)*CSTURB%XCTV*PEMOIST(:,:,:)**2 &
                           *ZWORK1(:,:,:)/CSTURB%XCTD/PD(:,:,:)     &
                           * (1. -  PREDTH1(:,:,:)*(1.5+PREDTH1(:,:,:)+PREDR1(:,:,:))/PD(:,:,:))
!$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
!
D_M3_WTH_WR2_O_DDTDZ(:,:,IKB-1)=D_M3_WTH_WR2_O_DDTDZ(:,:,IKB)
D_M3_WTH_WR2_O_DDTDZ(:,:,IKE+1)=D_M3_WTH_WR2_O_DDTDZ(:,:,IKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_WTH_WR2_O_DDTDZ',1,ZHOOK_HANDLE)
END FUNCTION D_M3_WTH_WR2_O_DDTDZ
!----------------------------------------------------------------------------
FUNCTION M3_WTH_WTHR(D,CSTURB,PREDR1,PD,PKEFF,PTKE,PSQRT_TKE,PBETA,PLEPS,PEMOIST)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PKEFF
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PTKE
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PBETA
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PEMOIST
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: M3_WTH_WTHR
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: ZWORK1 ! working array
  INTEGER :: IKB, IKE, JI,JJ,JK
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_WTH_WTHR',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PD,3)-JPVEXT_TURB

ZWORK1 = MZM(PBETA/PTKE*PSQRT_TKE, D%NKA, D%NKU, D%NKL)
!$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
M3_WTH_WTHR(:,:,:) = CSTURB%XCSHF*PKEFF(:,:,:)*PEMOIST(:,:,:)*ZWORK1(:,:,:) &
                         *0.5*PLEPS(:,:,:)/CSTURB%XCTD*(1+PREDR1(:,:,:))/PD(:,:,:)
!$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
!
M3_WTH_WTHR(:,:,IKB-1)=M3_WTH_WTHR(:,:,IKB)
M3_WTH_WTHR(:,:,IKE+1)=M3_WTH_WTHR(:,:,IKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_WTH_WTHR',1,ZHOOK_HANDLE)
END FUNCTION M3_WTH_WTHR
!----------------------------------------------------------------------------
FUNCTION D_M3_WTH_WTHR_O_DDTDZ(D,CSTURB,PM3_WTH_WTHR,PREDTH1,PREDR1,PD,PBLL_O_E,PETHETA)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PM3_WTH_WTHR
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PETHETA
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: D_M3_WTH_WTHR_O_DDTDZ
  INTEGER :: IKB, IKE, JI,JJ,JK
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_WTH_WTHR_O_DDTDZ',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PD,3)-JPVEXT_TURB

!$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
D_M3_WTH_WTHR_O_DDTDZ(:,:,:) = - PM3_WTH_WTHR(:,:,:) * (1.5+PREDTH1(:,:,:)+PREDR1(:,:,:))&
                               /PD(:,:,:)*CSTURB%XCTV*PBLL_O_E(:,:,:)*PETHETA(:,:,:)
!$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
!
D_M3_WTH_WTHR_O_DDTDZ(:,:,IKB-1)=D_M3_WTH_WTHR_O_DDTDZ(:,:,IKB)
D_M3_WTH_WTHR_O_DDTDZ(:,:,IKE+1)=D_M3_WTH_WTHR_O_DDTDZ(:,:,IKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_WTH_WTHR_O_DDTDZ',1,ZHOOK_HANDLE)
END FUNCTION D_M3_WTH_WTHR_O_DDTDZ
!----------------------------------------------------------------------------
FUNCTION M3_TH2_W2TH(D,CSTURB,PREDTH1,PREDR1,PD,PDTDZ,PLM,PLEPS,PTKE)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PDTDZ
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PLM
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PTKE
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: M3_TH2_W2TH
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: ZWORK1 ! working array
  INTEGER :: IKB, IKE, JI,JJ,JK
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_TH2_W2TH',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PD,3)-JPVEXT_TURB

ZWORK1 = MZF((1.-0.5*PREDR1*(1.+PREDR1)/PD)/(1.+PREDTH1)*PDTDZ, D%NKA, D%NKU, D%NKL)
!$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
M3_TH2_W2TH(:,:,:) = - ZWORK1(:,:,:) &
                       * 1.5*PLM(:,:,:)*PLEPS(:,:,:)/PTKE(:,:,:)*CSTURB%XCTV
!$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
!
M3_TH2_W2TH(:,:,IKB-1)=M3_TH2_W2TH(:,:,IKB)
M3_TH2_W2TH(:,:,IKE+1)=M3_TH2_W2TH(:,:,IKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_TH2_W2TH',1,ZHOOK_HANDLE)
END FUNCTION M3_TH2_W2TH
!----------------------------------------------------------------------------
FUNCTION D_M3_TH2_W2TH_O_DDTDZ(D,CSTURB,PREDTH1,PREDR1,PD,PLM,PLEPS,PTKE,OUSERV)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PLM
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PTKE
  LOGICAL,                INTENT(IN) :: OUSERV
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: D_M3_TH2_W2TH_O_DDTDZ
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: ZWORK1 ! working array
  INTEGER :: IKB, IKE, JI,JJ,JK
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_TH2_W2TH_O_DDTDZ',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PD,3)-JPVEXT_TURB

IF (OUSERV) THEN
!  D_M3_TH2_W2TH_O_DDTDZ(:,:,:) = - 1.5*PLM*PLEPS/PTKE*CSTURB%XCTV * MZF(                    &
!          (1.-0.5*PREDR1*(1.+PREDR1)/PD)*(1.-(1.5+PREDTH1+PREDR1)*(1.+PREDTH1)/PD )  &
!        / (1.+PREDTH1)**2, D%NKA, D%NKU, D%NKL)

  ZWORK1 = MZF((1.-0.5*PREDR1*(1.+PREDR1)/PD)*(1.-(1.5+PREDTH1+PREDR1)*   &
             PREDTH1*(1.+PREDTH1)/PD ) / (1.+PREDTH1)**2, D%NKA, D%NKU, D%NKL)

  !$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
  D_M3_TH2_W2TH_O_DDTDZ(:,:,:) = - 1.5*PLM(:,:,:)*PLEPS(:,:,:)/PTKE(:,:,:)*CSTURB%XCTV * &
                                 ZWORK1(:,:,:)
  !$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
ELSE
  ZWORK1 = MZF(1./(1.+PREDTH1)**2, D%NKA, D%NKU, D%NKL)
  !$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
  D_M3_TH2_W2TH_O_DDTDZ(:,:,:) = - 1.5*PLM(:,:,:)*PLEPS(:,:,:)/PTKE(:,:,:)*CSTURB%XCTV &
                                 * ZWORK1(:,:,:)
  !$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
END IF
!
D_M3_TH2_W2TH_O_DDTDZ(:,:,IKB-1)=D_M3_TH2_W2TH_O_DDTDZ(:,:,IKB)
D_M3_TH2_W2TH_O_DDTDZ(:,:,IKE+1)=D_M3_TH2_W2TH_O_DDTDZ(:,:,IKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_TH2_W2TH_O_DDTDZ',1,ZHOOK_HANDLE)
END FUNCTION D_M3_TH2_W2TH_O_DDTDZ
!----------------------------------------------------------------------------
FUNCTION M3_TH2_WTH2(D,CSTURB,PREDTH1,PREDR1,PD,PLEPS,PSQRT_TKE)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: M3_TH2_WTH2
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: ZWORK1 ! working array
  INTEGER :: IKB, IKE, JI,JJ,JK
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_TH2_WTH2',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PD,3)-JPVEXT_TURB

ZWORK1 = MZF((1.+0.5*PREDTH1+1.5*PREDR1+0.5*PREDR1**2)/PD, D%NKA, D%NKU, D%NKL)
!$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
M3_TH2_WTH2(:,:,:) = PLEPS(:,:,:)*0.5/CSTURB%XCTD/PSQRT_TKE(:,:,:)          &
                     * ZWORK1(:,:,:)
!$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
!
M3_TH2_WTH2(:,:,IKB-1)=M3_TH2_WTH2(:,:,IKB)
M3_TH2_WTH2(:,:,IKE+1)=M3_TH2_WTH2(:,:,IKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_TH2_WTH2',1,ZHOOK_HANDLE)
END FUNCTION M3_TH2_WTH2
!----------------------------------------------------------------------------
FUNCTION D_M3_TH2_WTH2_O_DDTDZ(D,CSTURB,PREDTH1,PREDR1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PETHETA)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PETHETA
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: D_M3_TH2_WTH2_O_DDTDZ
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: ZWORK1 ! working array
  INTEGER :: IKB, IKE, JI,JJ,JK
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_TH2_WTH2_O_DDTDZ',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PD,3)-JPVEXT_TURB

ZWORK1 = MZF(PBLL_O_E*PETHETA* (0.5/PD                                                   &
             - (1.5+PREDTH1+PREDR1)*(1.+0.5*PREDTH1+1.5*PREDR1+0.5*PREDR1**2)/PD**2 &
                           ), D%NKA, D%NKU, D%NKL)
 
!$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
D_M3_TH2_WTH2_O_DDTDZ(:,:,:) = PLEPS(:,:,:)*0.5/CSTURB%XCTD/PSQRT_TKE(:,:,:)*CSTURB%XCTV &
                               * ZWORK1(:,:,:)
!$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
!
D_M3_TH2_WTH2_O_DDTDZ(:,:,IKB-1)=D_M3_TH2_WTH2_O_DDTDZ(:,:,IKB)
D_M3_TH2_WTH2_O_DDTDZ(:,:,IKE+1)=D_M3_TH2_WTH2_O_DDTDZ(:,:,IKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_TH2_WTH2_O_DDTDZ',1,ZHOOK_HANDLE)
END FUNCTION D_M3_TH2_WTH2_O_DDTDZ
!----------------------------------------------------------------------------
FUNCTION M3_TH2_W2R(D,CSTURB,PD,PLM,PLEPS,PTKE,PBLL_O_E,PEMOIST,PDTDZ)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PLM
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PTKE
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PEMOIST
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PDTDZ
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: M3_TH2_W2R
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: ZWORK1 ! working array
  INTEGER :: IKB, IKE, JI,JJ,JK
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_TH2_W2R',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PD,3)-JPVEXT_TURB

ZWORK1 = MZF(PBLL_O_E*PEMOIST/PD*PDTDZ**2, D%NKA, D%NKU, D%NKL)
!$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
M3_TH2_W2R(:,:,:) = 0.75*CSTURB%XCTV**2*ZWORK1(:,:,:) &
                    *PLM(:,:,:)*PLEPS(:,:,:)/PTKE(:,:,:)
!$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
!
M3_TH2_W2R(:,:,IKB-1)=M3_TH2_W2R(:,:,IKB)
M3_TH2_W2R(:,:,IKE+1)=M3_TH2_W2R(:,:,IKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_TH2_W2R',1,ZHOOK_HANDLE)
END FUNCTION M3_TH2_W2R
!----------------------------------------------------------------------------
FUNCTION D_M3_TH2_W2R_O_DDTDZ(D,CSTURB,PREDTH1,PREDR1,PD,PLM,PLEPS,PTKE,PBLL_O_E,PEMOIST,PDTDZ)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PLM
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PTKE
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PEMOIST
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PDTDZ
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: D_M3_TH2_W2R_O_DDTDZ
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: ZWORK1 ! working array
  INTEGER :: IKB, IKE, JI,JJ,JK
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_TH2_W2R_O_DDTDZ',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PD,3)-JPVEXT_TURB

ZWORK1 =  MZF(PBLL_O_E*PEMOIST/PD*PDTDZ*(2.-PREDTH1*(1.5+PREDTH1+PREDR1)/PD), D%NKA, D%NKU, D%NKL)
!$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
D_M3_TH2_W2R_O_DDTDZ(:,:,:) = 0.75*CSTURB%XCTV**2*PLM(:,:,:)*PLEPS(:,:,:)/PTKE(:,:,:) &
                              * ZWORK1(:,:,:)
!$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
!
D_M3_TH2_W2R_O_DDTDZ(:,:,IKB-1)=D_M3_TH2_W2R_O_DDTDZ(:,:,IKB)
D_M3_TH2_W2R_O_DDTDZ(:,:,IKE+1)=D_M3_TH2_W2R_O_DDTDZ(:,:,IKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_TH2_W2R_O_DDTDZ',1,ZHOOK_HANDLE)
END FUNCTION D_M3_TH2_W2R_O_DDTDZ
!----------------------------------------------------------------------------
FUNCTION M3_TH2_WR2(D,CSTURB,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PEMOIST,PDTDZ)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PEMOIST
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PDTDZ
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: M3_TH2_WR2
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: ZWORK1 ! working array
  INTEGER :: IKB, IKE, JI,JJ,JK
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_TH2_WR2',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PD,3)-JPVEXT_TURB

ZWORK1 = MZF((PBLL_O_E*PEMOIST*PDTDZ)**2/PD, D%NKA, D%NKU, D%NKL)
!$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
M3_TH2_WR2(:,:,:) = 0.25*CSTURB%XCTV**2*ZWORK1(:,:,:)&
                    *PLEPS(:,:,:)/PSQRT_TKE(:,:,:)/CSTURB%XCTD
!$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
!
M3_TH2_WR2(:,:,IKB-1)=M3_TH2_WR2(:,:,IKB)
M3_TH2_WR2(:,:,IKE+1)=M3_TH2_WR2(:,:,IKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_TH2_WR2',1,ZHOOK_HANDLE)
END FUNCTION M3_TH2_WR2
!----------------------------------------------------------------------------
FUNCTION D_M3_TH2_WR2_O_DDTDZ(D,CSTURB,PREDTH1,PREDR1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PEMOIST,PDTDZ)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PEMOIST
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PDTDZ
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: D_M3_TH2_WR2_O_DDTDZ
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: ZWORK1 ! working array
  INTEGER :: IKB, IKE, JI,JJ,JK
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_TH2_WR2_O_DDTDZ',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PD,3)-JPVEXT_TURB

ZWORK1 = MZF((PBLL_O_E*PEMOIST)**2*PDTDZ/PD*(2.-PREDTH1*(1.5+PREDTH1+PREDR1)/PD), D%NKA, D%NKU, D%NKL)
!$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
D_M3_TH2_WR2_O_DDTDZ(:,:,:) = 0.25*CSTURB%XCTV**2*PLEPS(:,:,:)/PSQRT_TKE(:,:,:)/CSTURB%XCTD &
                              *  ZWORK1(:,:,:)
!$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
!
D_M3_TH2_WR2_O_DDTDZ(:,:,IKB-1)=D_M3_TH2_WR2_O_DDTDZ(:,:,IKB)
D_M3_TH2_WR2_O_DDTDZ(:,:,IKE+1)=D_M3_TH2_WR2_O_DDTDZ(:,:,IKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_TH2_WR2_O_DDTDZ',1,ZHOOK_HANDLE)
END FUNCTION D_M3_TH2_WR2_O_DDTDZ
!----------------------------------------------------------------------------
FUNCTION M3_TH2_WTHR(D,CSTURB,PREDR1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PEMOIST,PDTDZ)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PEMOIST
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PDTDZ
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: M3_TH2_WTHR
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: ZWORK1 ! working array
  INTEGER :: IKB, IKE, JI,JJ,JK
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_TH2_WTHR',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PD,3)-JPVEXT_TURB

ZWORK1 = MZF(PBLL_O_E*PEMOIST*PDTDZ*(1.+PREDR1)/PD, D%NKA, D%NKU, D%NKL)
!$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
M3_TH2_WTHR(:,:,:) = - 0.5*CSTURB%XCTV*PLEPS(:,:,:)/PSQRT_TKE(:,:,:)/CSTURB%XCTD &
                       *ZWORK1(:,:,:)
!$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
!
M3_TH2_WTHR(:,:,IKB-1)=M3_TH2_WTHR(:,:,IKB)
M3_TH2_WTHR(:,:,IKE+1)=M3_TH2_WTHR(:,:,IKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_TH2_WTHR',1,ZHOOK_HANDLE)
END FUNCTION M3_TH2_WTHR
!----------------------------------------------------------------------------
FUNCTION D_M3_TH2_WTHR_O_DDTDZ(D,CSTURB,PREDTH1,PREDR1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PEMOIST,PDTDZ)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PEMOIST
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PDTDZ
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: D_M3_TH2_WTHR_O_DDTDZ
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: ZWORK1 ! working array
  INTEGER :: IKB, IKE, JI,JJ,JK
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_TH2_WTHR_O_DDTDZ',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PD,3)-JPVEXT_TURB

ZWORK1 = MZF(PBLL_O_E*PEMOIST*(1.+PREDR1)/PD * (1. -PREDTH1*(1.5+PREDTH1+PREDR1)/PD), D%NKA, D%NKU, D%NKL)
!$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
D_M3_TH2_WTHR_O_DDTDZ(:,:,:) = - 0.5*CSTURB%XCTV*PLEPS(:,:,:)/PSQRT_TKE(:,:,:)/CSTURB%XCTD &
                               * ZWORK1(:,:,:) 
!$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
!
D_M3_TH2_WTHR_O_DDTDZ(:,:,IKB-1)=D_M3_TH2_WTHR_O_DDTDZ(:,:,IKB)
D_M3_TH2_WTHR_O_DDTDZ(:,:,IKE+1)=D_M3_TH2_WTHR_O_DDTDZ(:,:,IKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_TH2_WTHR_O_DDTDZ',1,ZHOOK_HANDLE)
END FUNCTION D_M3_TH2_WTHR_O_DDTDZ
!----------------------------------------------------------------------------
FUNCTION M3_THR_WTHR(D,CSTURB,PREDTH1,PREDR1,PD,PLEPS,PSQRT_TKE)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: M3_THR_WTHR
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: ZWORK1 ! working array
  INTEGER :: IKB, IKE, JI,JJ,JK
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_THR_WTHR',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PD,3)-JPVEXT_TURB

ZWORK1 =  MZF((1.+PREDTH1)*(1.+PREDR1)/PD, D%NKA, D%NKU, D%NKL)
!$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
M3_THR_WTHR(:,:,:) = 0.5*PLEPS(:,:,:)/PSQRT_TKE(:,:,:)/CSTURB%XCTD &
                     * ZWORK1(:,:,:)
!$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
!
M3_THR_WTHR(:,:,IKB-1)=M3_THR_WTHR(:,:,IKB)
M3_THR_WTHR(:,:,IKE+1)=M3_THR_WTHR(:,:,IKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_THR_WTHR',1,ZHOOK_HANDLE)
END FUNCTION M3_THR_WTHR
!----------------------------------------------------------------------------
FUNCTION D_M3_THR_WTHR_O_DDTDZ(D,CSTURB,PREDTH1,PREDR1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PETHETA)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PETHETA
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: D_M3_THR_WTHR_O_DDTDZ
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: ZWORK1 ! working array
  INTEGER :: IKB, IKE, JI,JJ,JK
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_THR_WTHR_O_DDTDZ',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PD,3)-JPVEXT_TURB

ZWORK1 = MZF(PETHETA*PBLL_O_E/PD*(1.+PREDR1)*(1.-(1.+PREDTH1)*(1.5+PREDTH1+PREDR1)/PD), D%NKA, D%NKU, D%NKL)
!$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
D_M3_THR_WTHR_O_DDTDZ(:,:,:) = 0.5*PLEPS(:,:,:)/PSQRT_TKE(:,:,:)/CSTURB%XCTD * CSTURB%XCTV &
                               * ZWORK1(:,:,:)
!$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
!
D_M3_THR_WTHR_O_DDTDZ(:,:,IKB-1)=D_M3_THR_WTHR_O_DDTDZ(:,:,IKB)
D_M3_THR_WTHR_O_DDTDZ(:,:,IKE+1)=D_M3_THR_WTHR_O_DDTDZ(:,:,IKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_THR_WTHR_O_DDTDZ',1,ZHOOK_HANDLE)
END FUNCTION D_M3_THR_WTHR_O_DDTDZ
!----------------------------------------------------------------------------
FUNCTION M3_THR_WTH2(D,CSTURB,PREDR1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PETHETA,PDRDZ)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PETHETA
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PDRDZ
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: M3_THR_WTH2
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: ZWORK1 ! working array
  INTEGER :: IKB, IKE, JI,JJ,JK
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_THR_WTH2',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PD,3)-JPVEXT_TURB

ZWORK1 = MZF((1.+PREDR1)*PBLL_O_E*PETHETA*PDRDZ/PD, D%NKA, D%NKU, D%NKL)
!$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
M3_THR_WTH2(:,:,:) = - 0.25*PLEPS(:,:,:)/PSQRT_TKE(:,:,:)/CSTURB%XCTD*CSTURB%XCTV &
                     * ZWORK1(:,:,:)
!$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
!
M3_THR_WTH2(:,:,IKB-1)=M3_THR_WTH2(:,:,IKB)
M3_THR_WTH2(:,:,IKE+1)=M3_THR_WTH2(:,:,IKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_THR_WTH2',1,ZHOOK_HANDLE)
END FUNCTION M3_THR_WTH2
!----------------------------------------------------------------------------
FUNCTION D_M3_THR_WTH2_O_DDTDZ(D,CSTURB,PREDTH1,PREDR1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PETHETA,PDRDZ)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PETHETA
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PDRDZ
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: D_M3_THR_WTH2_O_DDTDZ
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: ZWORK1 ! working array
  INTEGER :: IKB, IKE, JI,JJ,JK
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_THR_WTH2_O_DDTDZ',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PD,3)-JPVEXT_TURB

ZWORK1 = MZF(-(1.+PREDR1)*(PBLL_O_E*PETHETA/PD)**2*PDRDZ*(1.5+PREDTH1+PREDR1), D%NKA, D%NKU, D%NKL)
!$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
D_M3_THR_WTH2_O_DDTDZ(:,:,:) = - 0.25*PLEPS(:,:,:)/PSQRT_TKE(:,:,:)/CSTURB%XCTD*CSTURB%XCTV**2 &
                               * ZWORK1(:,:,:) 
!$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
!
D_M3_THR_WTH2_O_DDTDZ(:,:,IKB-1)=D_M3_THR_WTH2_O_DDTDZ(:,:,IKB)
D_M3_THR_WTH2_O_DDTDZ(:,:,IKE+1)=D_M3_THR_WTH2_O_DDTDZ(:,:,IKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_THR_WTH2_O_DDTDZ',1,ZHOOK_HANDLE)
END FUNCTION D_M3_THR_WTH2_O_DDTDZ
!----------------------------------------------------------------------------
FUNCTION D_M3_THR_WTH2_O_DDRDZ(D,CSTURB,PREDTH1,PREDR1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PETHETA)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PETHETA
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: D_M3_THR_WTH2_O_DDRDZ
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: ZWORK1 ! working array
  INTEGER :: IKB, IKE, JI,JJ,JK
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_THR_WTH2_O_DDRDZ',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PD,3)-JPVEXT_TURB

ZWORK1 = MZF(PBLL_O_E*PETHETA/PD                                              &
       *(-(1.+PREDR1)*PREDR1/PD*(1.5+PREDTH1+PREDR1)+(1.+2.*PREDR1)),     &
       D%NKA, D%NKU, D%NKL)
!$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
D_M3_THR_WTH2_O_DDRDZ(:,:,:) = - 0.25*PLEPS(:,:,:)/PSQRT_TKE(:,:,:)/CSTURB%XCTD*CSTURB%XCTV          &
                               * ZWORK1(:,:,:)
!$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
!
D_M3_THR_WTH2_O_DDRDZ(:,:,IKB-1)=D_M3_THR_WTH2_O_DDRDZ(:,:,IKB)
D_M3_THR_WTH2_O_DDRDZ(:,:,IKE+1)=D_M3_THR_WTH2_O_DDRDZ(:,:,IKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_THR_WTH2_O_DDRDZ',1,ZHOOK_HANDLE)
END FUNCTION D_M3_THR_WTH2_O_DDRDZ
!----------------------------------------------------------------------------
FUNCTION M3_THR_W2TH(D,CSTURB,PREDR1,PD,PLM,PLEPS,PTKE,PDRDZ)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PLM
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PTKE
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PDRDZ
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: M3_THR_W2TH
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: ZWORK1 ! working array
  INTEGER :: IKB, IKE, JI,JJ,JK
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_THR_W2TH',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PD,3)-JPVEXT_TURB

ZWORK1 = MZF((1.+PREDR1)*PDRDZ/PD, D%NKA, D%NKU, D%NKL)
!$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
M3_THR_W2TH(:,:,:) = - 0.75*PLM(:,:,:)*PLEPS(:,:,:)/PTKE(:,:,:) * CSTURB%XCTV      &
                     * ZWORK1(:,:,:)
!$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
!
M3_THR_W2TH(:,:,IKB-1)=M3_THR_W2TH(:,:,IKB)
M3_THR_W2TH(:,:,IKE+1)=M3_THR_W2TH(:,:,IKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_THR_W2TH',1,ZHOOK_HANDLE)
END FUNCTION M3_THR_W2TH
!----------------------------------------------------------------------------
FUNCTION D_M3_THR_W2TH_O_DDTDZ(D,CSTURB,PREDTH1,PREDR1,PD,PLM,PLEPS,PTKE,PBLL_O_E,PDRDZ,PETHETA)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PLM
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PTKE
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PDRDZ
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PETHETA
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: D_M3_THR_W2TH_O_DDTDZ
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: ZWORK1 ! working array
  INTEGER :: IKB, IKE, JI,JJ,JK
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_THR_W2TH_O_DDTDZ',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PD,3)-JPVEXT_TURB

ZWORK1 =  MZF(-PETHETA*PBLL_O_E*(1.+PREDR1)*PDRDZ*(1.5+PREDTH1+PREDR1)/PD**2, D%NKA, D%NKU, D%NKL)
!$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
D_M3_THR_W2TH_O_DDTDZ(:,:,:) = - 0.75*PLM(:,:,:)*PLEPS(:,:,:)/PTKE(:,:,:) * CSTURB%XCTV**2    &
                               * ZWORK1(:,:,:)
!$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
!
D_M3_THR_W2TH_O_DDTDZ(:,:,IKB-1)=D_M3_THR_W2TH_O_DDTDZ(:,:,IKB)
D_M3_THR_W2TH_O_DDTDZ(:,:,IKE+1)=D_M3_THR_W2TH_O_DDTDZ(:,:,IKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_THR_W2TH_O_DDTDZ',1,ZHOOK_HANDLE)
END FUNCTION D_M3_THR_W2TH_O_DDTDZ
!----------------------------------------------------------------------------
FUNCTION D_M3_THR_W2TH_O_DDRDZ(D,CSTURB,PREDTH1,PREDR1,PD,PLM,PLEPS,PTKE)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PLM
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PTKE
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: D_M3_THR_W2TH_O_DDRDZ
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: ZWORK1 ! working array
  INTEGER :: IKB, IKE, JI,JJ,JK
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_THR_W2TH_O_DDRDZ',0,ZHOOK_HANDLE)
IKB = 1+JPVEXT_TURB
IKE = SIZE(PD,3)-JPVEXT_TURB

ZWORK1 =  MZF(-(1.+PREDR1)*PREDR1*(1.5+PREDTH1+PREDR1)/PD**2          &
        +(1.+2.*PREDR1)/PD,                                    &
       D%NKA, D%NKU, D%NKL)
!$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
D_M3_THR_W2TH_O_DDRDZ(:,:,:) = - 0.75*PLM(:,:,:)*PLEPS(:,:,:)/PTKE(:,:,:) * CSTURB%XCTV     &
                               * ZWORK1(:,:,:)
!$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
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
FUNCTION PSI3(D,CSTURB,PREDR1,PREDTH1,PRED2R3,PRED2TH3,PRED2THR3,HTURBDIM,OUSERV)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PRED2TH3
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PRED2R3
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PRED2THR3
  CHARACTER(LEN=4),       INTENT(IN) :: HTURBDIM  ! 1DIM or 3DIM turb. scheme
  LOGICAL,                INTENT(IN) :: OUSERV    ! flag to use vapor
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: PSI3
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:PSI3',0,ZHOOK_HANDLE)
PSI3 = PHI3(D,CSTURB,PREDR1,PREDTH1,PRED2R3,PRED2TH3,PRED2THR3,HTURBDIM,OUSERV)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:PSI3',1,ZHOOK_HANDLE)
END FUNCTION PSI3
!----------------------------------------------------------------------------
FUNCTION D_PSI3DRDZ_O_DDRDZ(D,CSTURB,PPSI3,PREDR1,PREDTH1,PRED2R3,PRED2THR3,HTURBDIM,OUSERV)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PPSI3
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PRED2R3
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PRED2THR3
  CHARACTER(LEN=4),       INTENT(IN) :: HTURBDIM  ! 1DIM or 3DIM turb. scheme
  LOGICAL,                INTENT(IN) :: OUSERV    ! flag to use vapor
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: D_PSI3DRDZ_O_DDRDZ

REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_PSI3DRDZ_O_DDRDZ',0,ZHOOK_HANDLE)
D_PSI3DRDZ_O_DDRDZ = D_PHI3DTDZ_O_DDTDZ(D,CSTURB,PPSI3,PREDR1,PREDTH1,PRED2R3,PRED2THR3,HTURBDIM,OUSERV)
!
!C'est ok?!
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_PSI3DRDZ_O_DDRDZ',1,ZHOOK_HANDLE)
END FUNCTION D_PSI3DRDZ_O_DDRDZ
!----------------------------------------------------------------------------
FUNCTION D_PSI3DTDZ_O_DDTDZ(D,CSTURB,PPSI3,PREDR1,PREDTH1,PRED2R3,PRED2THR3,HTURBDIM,OUSERV)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PPSI3
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PRED2R3
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PRED2THR3
  CHARACTER(LEN=4),       INTENT(IN) :: HTURBDIM  ! 1DIM or 3DIM turb. scheme
  LOGICAL,                INTENT(IN) :: OUSERV    ! flag to use vapor
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: D_PSI3DTDZ_O_DDTDZ
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_PSI3DTDZ_O_DDTDZ',0,ZHOOK_HANDLE)
D_PSI3DTDZ_O_DDTDZ = D_PHI3DRDZ_O_DDRDZ(D,CSTURB,PPSI3,PREDR1,PREDTH1,PRED2R3,PRED2THR3,HTURBDIM,OUSERV)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_PSI3DTDZ_O_DDTDZ',1,ZHOOK_HANDLE)
END FUNCTION D_PSI3DTDZ_O_DDTDZ
!----------------------------------------------------------------------------
FUNCTION D_PSI3DRDZ2_O_DDRDZ(D,CSTURB,PPSI3,PREDR1,PREDTH1,PRED2R3,PRED2THR3,PDRDZ,HTURBDIM,OUSERV)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PPSI3
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PRED2R3
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PRED2THR3
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PDRDZ
  CHARACTER(LEN=4),       INTENT(IN) :: HTURBDIM  ! 1DIM or 3DIM turb. scheme
  LOGICAL,                INTENT(IN) :: OUSERV    ! flag to use vapor
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: D_PSI3DRDZ2_O_DDRDZ
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_PSI3DRDZ2_O_DDRDZ',0,ZHOOK_HANDLE)
D_PSI3DRDZ2_O_DDRDZ = D_PHI3DTDZ2_O_DDTDZ(D,CSTURB,PPSI3,PREDR1,PREDTH1,PRED2R3,PRED2THR3,PDRDZ,HTURBDIM,OUSERV)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_PSI3DRDZ2_O_DDRDZ',1,ZHOOK_HANDLE)
END FUNCTION D_PSI3DRDZ2_O_DDRDZ
!----------------------------------------------------------------------------
FUNCTION M3_WR_WR2(D,CSTURB,PREDR1,PREDTH1,PD,PBLL_O_E,PEMOIST)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PEMOIST
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: M3_WR_WR2
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_WR_WR2',0,ZHOOK_HANDLE)
M3_WR_WR2 = M3_WTH_WTH2(D,CSTURB,PREDR1,PREDTH1,PD,PBLL_O_E,PEMOIST)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_WR_WR2',1,ZHOOK_HANDLE)
END FUNCTION M3_WR_WR2
!----------------------------------------------------------------------------
FUNCTION D_M3_WR_WR2_O_DDRDZ(D,CSTURB,PM3_WR_WR2,PREDR1,PREDTH1,PD,PBLL_O_E,PEMOIST)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PM3_WR_WR2
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PEMOIST
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: D_M3_WR_WR2_O_DDRDZ
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_WR_WR2_O_DDRDZ',0,ZHOOK_HANDLE)
D_M3_WR_WR2_O_DDRDZ = D_M3_WTH_WTH2_O_DDTDZ(D,CSTURB,PM3_WR_WR2,PREDR1,PREDTH1,PD,PBLL_O_E,PEMOIST)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_WR_WR2_O_DDRDZ',1,ZHOOK_HANDLE)
END FUNCTION D_M3_WR_WR2_O_DDRDZ
!----------------------------------------------------------------------------
FUNCTION M3_WR_W2R(D,CSTURB,PREDR1,PREDTH1,PD,PKEFF,PTKE)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PKEFF
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PTKE
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: M3_WR_W2R
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_WR_W2R',0,ZHOOK_HANDLE)
M3_WR_W2R = M3_WTH_W2TH(D,CSTURB,PREDR1,PREDTH1,PD,PKEFF,PTKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_WR_W2R',1,ZHOOK_HANDLE)
END FUNCTION M3_WR_W2R
!----------------------------------------------------------------------------
FUNCTION D_M3_WR_W2R_O_DDRDZ(D,CSTURB,PREDR1,PREDTH1,PD,PBLL_O_E,PEMOIST,PKEFF,PTKE)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PEMOIST
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PKEFF
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PTKE
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: D_M3_WR_W2R_O_DDRDZ
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_WR_W2R_O_DDRDZ',0,ZHOOK_HANDLE)
D_M3_WR_W2R_O_DDRDZ = D_M3_WTH_W2TH_O_DDTDZ(D,CSTURB,PREDR1,PREDTH1,PD,PBLL_O_E,PEMOIST,PKEFF,PTKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_WR_W2R_O_DDRDZ',1,ZHOOK_HANDLE)
END FUNCTION D_M3_WR_W2R_O_DDRDZ
!----------------------------------------------------------------------------
FUNCTION M3_WR_W2TH(D,CSTURB,PD,PKEFF,PTKE,PBLL_O_E,PETHETA,PDRDZ)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PKEFF
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PTKE
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PETHETA
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PDRDZ
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: M3_WR_W2TH
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_WR_W2TH',0,ZHOOK_HANDLE)
M3_WR_W2TH = M3_WTH_W2R(D,CSTURB,PD,PKEFF,PTKE,PBLL_O_E,PETHETA,PDRDZ)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_WR_W2TH',1,ZHOOK_HANDLE)
END FUNCTION M3_WR_W2TH
!----------------------------------------------------------------------------
FUNCTION D_M3_WR_W2TH_O_DDRDZ(D,CSTURB,PREDR1,PREDTH1,PD,PKEFF,PTKE,PBLL_O_E,PETHETA)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PKEFF
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PTKE
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PETHETA
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: D_M3_WR_W2TH_O_DDRDZ
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_WR_W2TH_O_DDRDZ',0,ZHOOK_HANDLE)
D_M3_WR_W2TH_O_DDRDZ = D_M3_WTH_W2R_O_DDTDZ(D,CSTURB,PREDR1,PREDTH1,PD,PKEFF,PTKE,PBLL_O_E,PETHETA)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_WR_W2TH_O_DDRDZ',1,ZHOOK_HANDLE)
END FUNCTION D_M3_WR_W2TH_O_DDRDZ
!----------------------------------------------------------------------------
FUNCTION M3_WR_WTH2(D,CSTURB,PD,PKEFF,PTKE,PSQRT_TKE,PBLL_O_E,PBETA,PLEPS,PETHETA,PDRDZ)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PKEFF
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PTKE
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PBETA
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PETHETA
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PDRDZ
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: M3_WR_WTH2
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_WR_WTH2',0,ZHOOK_HANDLE)
M3_WR_WTH2 = M3_WTH_WR2(D,CSTURB,PD,PKEFF,PTKE,PSQRT_TKE,PBLL_O_E,PBETA,PLEPS,PETHETA,PDRDZ)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_WR_WTH2',1,ZHOOK_HANDLE)
END FUNCTION M3_WR_WTH2
!----------------------------------------------------------------------------
FUNCTION D_M3_WR_WTH2_O_DDRDZ(D,CSTURB,PREDR1,PREDTH1,PD,PKEFF,PTKE,PSQRT_TKE,PBLL_O_E,PBETA,PLEPS,PETHETA)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PKEFF
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PTKE
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PBETA
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PETHETA
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: D_M3_WR_WTH2_O_DDRDZ
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_WR_WTH2_O_DDRDZ',0,ZHOOK_HANDLE)
D_M3_WR_WTH2_O_DDRDZ=D_M3_WTH_WR2_O_DDTDZ(D,CSTURB,PREDR1,PREDTH1,PD,PKEFF,PTKE,PSQRT_TKE,PBLL_O_E,PBETA,PLEPS,PETHETA)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_WR_WTH2_O_DDRDZ',1,ZHOOK_HANDLE)
END FUNCTION D_M3_WR_WTH2_O_DDRDZ
!----------------------------------------------------------------------------
FUNCTION M3_WR_WTHR(D,CSTURB,PREDTH1,PD,PKEFF,PTKE,PSQRT_TKE,PBETA,PLEPS,PETHETA)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PKEFF
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PTKE
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PBETA
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PETHETA
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: M3_WR_WTHR
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_WR_WTHR',0,ZHOOK_HANDLE)
M3_WR_WTHR = M3_WTH_WTHR(D,CSTURB,PREDTH1,PD,PKEFF,PTKE,PSQRT_TKE,PBETA,PLEPS,PETHETA)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_WR_WTHR',1,ZHOOK_HANDLE)
END FUNCTION M3_WR_WTHR
!----------------------------------------------------------------------------
FUNCTION D_M3_WR_WTHR_O_DDRDZ(D,CSTURB,PM3_WR_WTHR,PREDR1,PREDTH1,PD,PBLL_O_E,PEMOIST)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PM3_WR_WTHR
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PEMOIST
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: D_M3_WR_WTHR_O_DDRDZ
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_WR_WTHR_O_DDRDZ',0,ZHOOK_HANDLE)
D_M3_WR_WTHR_O_DDRDZ = D_M3_WTH_WTHR_O_DDTDZ(D,CSTURB,PM3_WR_WTHR,PREDR1,PREDTH1,PD,PBLL_O_E,PEMOIST)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_WR_WTHR_O_DDRDZ',1,ZHOOK_HANDLE)
END FUNCTION D_M3_WR_WTHR_O_DDRDZ
!----------------------------------------------------------------------------
FUNCTION M3_R2_W2R(D,CSTURB,PREDR1,PREDTH1,PD,PDRDZ,PLM,PLEPS,PTKE)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PDRDZ
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PLM
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PTKE
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: M3_R2_W2R
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_R2_W2R',0,ZHOOK_HANDLE)
M3_R2_W2R = M3_TH2_W2TH(D,CSTURB,PREDR1,PREDTH1,PD,PDRDZ,PLM,PLEPS,PTKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_R2_W2R',1,ZHOOK_HANDLE)
END FUNCTION M3_R2_W2R
!----------------------------------------------------------------------------
FUNCTION D_M3_R2_W2R_O_DDRDZ(D,CSTURB,PREDR1,PREDTH1,PD,PLM,PLEPS,PTKE,OUSERV)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PLM
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PTKE
  LOGICAL,                INTENT(IN) :: OUSERV
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: D_M3_R2_W2R_O_DDRDZ
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_R2_W2R_O_DDRDZ',0,ZHOOK_HANDLE)
D_M3_R2_W2R_O_DDRDZ = D_M3_TH2_W2TH_O_DDTDZ(D,CSTURB,PREDR1,PREDTH1,PD,PLM,PLEPS,PTKE,OUSERV)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_R2_W2R_O_DDRDZ',1,ZHOOK_HANDLE)
END FUNCTION D_M3_R2_W2R_O_DDRDZ
!----------------------------------------------------------------------------
FUNCTION M3_R2_WR2(D,CSTURB,PREDR1,PREDTH1,PD,PLEPS,PSQRT_TKE)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: M3_R2_WR2
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_R2_WR2',0,ZHOOK_HANDLE)
M3_R2_WR2 = M3_TH2_WTH2(D,CSTURB,PREDR1,PREDTH1,PD,PLEPS,PSQRT_TKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_R2_WR2',1,ZHOOK_HANDLE)
END FUNCTION M3_R2_WR2
!----------------------------------------------------------------------------
FUNCTION D_M3_R2_WR2_O_DDRDZ(D,CSTURB,PREDR1,PREDTH1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PEMOIST)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PEMOIST
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: D_M3_R2_WR2_O_DDRDZ
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_R2_WR2_O_DDRDZ',0,ZHOOK_HANDLE)
D_M3_R2_WR2_O_DDRDZ = D_M3_TH2_WTH2_O_DDTDZ(D,CSTURB,PREDR1,PREDTH1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PEMOIST)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_R2_WR2_O_DDRDZ',1,ZHOOK_HANDLE)
END FUNCTION D_M3_R2_WR2_O_DDRDZ
!----------------------------------------------------------------------------
FUNCTION M3_R2_W2TH(D,CSTURB,PD,PLM,PLEPS,PTKE,PBLL_O_E,PETHETA,PDRDZ)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PLM
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PTKE
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PETHETA
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PDRDZ
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: M3_R2_W2TH
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_R2_W2TH',0,ZHOOK_HANDLE)
M3_R2_W2TH = M3_TH2_W2R(D,CSTURB,PD,PLM,PLEPS,PTKE,PBLL_O_E,PETHETA,PDRDZ)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_R2_W2TH',1,ZHOOK_HANDLE)
END FUNCTION M3_R2_W2TH
!----------------------------------------------------------------------------
FUNCTION D_M3_R2_W2TH_O_DDRDZ(D,CSTURB,PREDR1,PREDTH1,PD,PLM,PLEPS,PTKE,PBLL_O_E,PETHETA,PDRDZ)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PLM
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PTKE
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PETHETA
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PDRDZ
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: D_M3_R2_W2TH_O_DDRDZ
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_R2_W2TH_O_DDRDZ',0,ZHOOK_HANDLE)
D_M3_R2_W2TH_O_DDRDZ = D_M3_TH2_W2R_O_DDTDZ(D,CSTURB,PREDR1,PREDTH1,PD,PLM,PLEPS,PTKE,PBLL_O_E,PETHETA,PDRDZ)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_R2_W2TH_O_DDRDZ',1,ZHOOK_HANDLE)
END FUNCTION D_M3_R2_W2TH_O_DDRDZ
!----------------------------------------------------------------------------
FUNCTION M3_R2_WTH2(D,CSTURB,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PETHETA,PDRDZ)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PETHETA
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PDRDZ
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: M3_R2_WTH2
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_R2_WTH2',0,ZHOOK_HANDLE)
M3_R2_WTH2 = M3_TH2_WR2(D,CSTURB,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PETHETA,PDRDZ)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_R2_WTH2',1,ZHOOK_HANDLE)
END FUNCTION M3_R2_WTH2
!----------------------------------------------------------------------------
FUNCTION D_M3_R2_WTH2_O_DDRDZ(D,CSTURB,PREDR1,PREDTH1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PETHETA,PDRDZ)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PETHETA
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PDRDZ
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: D_M3_R2_WTH2_O_DDRDZ
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_R2_WTH2_O_DDRDZ',0,ZHOOK_HANDLE)
D_M3_R2_WTH2_O_DDRDZ = D_M3_TH2_WR2_O_DDTDZ(D,CSTURB,PREDR1,PREDTH1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PETHETA,PDRDZ)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_R2_WTH2_O_DDRDZ',1,ZHOOK_HANDLE)
END FUNCTION D_M3_R2_WTH2_O_DDRDZ
!----------------------------------------------------------------------------
FUNCTION M3_R2_WTHR(D,CSTURB,PREDTH1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PETHETA,PDRDZ)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PETHETA
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PDRDZ
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: M3_R2_WTHR
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_R2_WTHR',0,ZHOOK_HANDLE)
M3_R2_WTHR = M3_TH2_WTHR(D,CSTURB,PREDTH1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PETHETA,PDRDZ)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_R2_WTHR',1,ZHOOK_HANDLE)
END FUNCTION M3_R2_WTHR
!----------------------------------------------------------------------------
FUNCTION D_M3_R2_WTHR_O_DDRDZ(D,CSTURB,PREDR1,PREDTH1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PETHETA,PDRDZ)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PETHETA
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PDRDZ
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: D_M3_R2_WTHR_O_DDRDZ
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_R2_WTHR_O_DDRDZ',0,ZHOOK_HANDLE)
D_M3_R2_WTHR_O_DDRDZ = D_M3_TH2_WTHR_O_DDTDZ(D,CSTURB,PREDR1,PREDTH1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PETHETA,PDRDZ)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_R2_WTHR_O_DDRDZ',1,ZHOOK_HANDLE)
END FUNCTION D_M3_R2_WTHR_O_DDRDZ
!----------------------------------------------------------------------------
FUNCTION D_M3_THR_WTHR_O_DDRDZ(D,CSTURB,PREDR1,PREDTH1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PEMOIST)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PEMOIST
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: D_M3_THR_WTHR_O_DDRDZ
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_THR_WTHR_O_DDRDZ',0,ZHOOK_HANDLE)
D_M3_THR_WTHR_O_DDRDZ = D_M3_THR_WTHR_O_DDTDZ(D,CSTURB,PREDR1,PREDTH1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PEMOIST)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_THR_WTHR_O_DDRDZ',1,ZHOOK_HANDLE)
END FUNCTION D_M3_THR_WTHR_O_DDRDZ
!----------------------------------------------------------------------------
FUNCTION M3_THR_WR2(D,CSTURB,PREDTH1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PEMOIST,PDTDZ)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PEMOIST
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PDTDZ
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: M3_THR_WR2
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_THR_WR2',0,ZHOOK_HANDLE)
M3_THR_WR2 = M3_THR_WTH2(D,CSTURB,PREDTH1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PEMOIST,PDTDZ)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_THR_WR2',1,ZHOOK_HANDLE)
END FUNCTION M3_THR_WR2
!----------------------------------------------------------------------------
FUNCTION D_M3_THR_WR2_O_DDRDZ(D,CSTURB,PREDR1,PREDTH1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PEMOIST,PDTDZ)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PEMOIST
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PDTDZ
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: D_M3_THR_WR2_O_DDRDZ
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_THR_WR2_O_DDRDZ',0,ZHOOK_HANDLE)
D_M3_THR_WR2_O_DDRDZ = D_M3_THR_WTH2_O_DDTDZ(D,CSTURB,PREDR1,PREDTH1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PEMOIST,PDTDZ)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_THR_WR2_O_DDRDZ',1,ZHOOK_HANDLE)
END FUNCTION D_M3_THR_WR2_O_DDRDZ
!----------------------------------------------------------------------------
FUNCTION D_M3_THR_WR2_O_DDTDZ(D,CSTURB,PREDR1,PREDTH1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PEMOIST)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PSQRT_TKE
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PEMOIST
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: D_M3_THR_WR2_O_DDTDZ
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_THR_WR2_O_DDTDZ',0,ZHOOK_HANDLE)
D_M3_THR_WR2_O_DDTDZ = D_M3_THR_WTH2_O_DDRDZ(D,CSTURB,PREDR1,PREDTH1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PEMOIST)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_THR_WR2_O_DDTDZ',1,ZHOOK_HANDLE)
END FUNCTION D_M3_THR_WR2_O_DDTDZ
!----------------------------------------------------------------------------
FUNCTION M3_THR_W2R(D,CSTURB,PREDTH1,PD,PLM,PLEPS,PTKE,PDTDZ)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PLM
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PTKE
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PDTDZ
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: M3_THR_W2R
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_THR_W2R',0,ZHOOK_HANDLE)
M3_THR_W2R = M3_THR_W2TH(D,CSTURB,PREDTH1,PD,PLM,PLEPS,PTKE,PDTDZ)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:M3_THR_W2R',1,ZHOOK_HANDLE)
END FUNCTION M3_THR_W2R
!----------------------------------------------------------------------------
FUNCTION D_M3_THR_W2R_O_DDRDZ(D,CSTURB,PREDR1,PREDTH1,PD,PLM,PLEPS,PTKE,PBLL_O_E,PDTDZ,PEMOIST)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PLM
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PTKE
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PBLL_O_E
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PDTDZ
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PEMOIST
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: D_M3_THR_W2R_O_DDRDZ
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_THR_W2R_O_DDRDZ',0,ZHOOK_HANDLE)
D_M3_THR_W2R_O_DDRDZ = D_M3_THR_W2TH_O_DDTDZ(D,CSTURB,PREDR1,PREDTH1,PD,PLM,PLEPS,PTKE,PBLL_O_E,PDTDZ,PEMOIST)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_THR_W2R_O_DDRDZ',1,ZHOOK_HANDLE)
END FUNCTION D_M3_THR_W2R_O_DDRDZ
!----------------------------------------------------------------------------
FUNCTION D_M3_THR_W2R_O_DDTDZ(D,CSTURB,PREDR1,PREDTH1,PD,PLM,PLEPS,PTKE)
  TYPE(DIMPHYEX_t),                   INTENT(IN) :: D
  TYPE(CSTURB_t),                     INTENT(IN) :: CSTURB
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDR1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PREDTH1
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PD
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PLM
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PLEPS
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PTKE
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: D_M3_THR_W2R_O_DDTDZ
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_THR_W2R_O_DDTDZ',0,ZHOOK_HANDLE)
D_M3_THR_W2R_O_DDTDZ = D_M3_THR_W2TH_O_DDRDZ(D,CSTURB,PREDR1,PREDTH1,PD,PLM,PLEPS,PTKE)
!
IF (LHOOK) CALL DR_HOOK('MODE_PRANDTL:D_M3_THR_W2R_O_DDTDZ',1,ZHOOK_HANDLE)
END FUNCTION D_M3_THR_W2R_O_DDTDZ
!----------------------------------------------------------------------------
!
END MODULE MODE_PRANDTL

