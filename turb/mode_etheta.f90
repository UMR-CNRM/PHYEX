!MNH_LIC Copyright 1995-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
MODULE MODE_ETHETA
IMPLICIT NONE
CONTAINS
SUBROUTINE ETHETA(D,CST,KRR,KRRI,PTHLM,PRM,PLOCPEXNM,PATHETA,PSRCM,OOCEAN,OCOMPUTE_SRC,PETHETA)
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!   ############################################################################
!
!      PURPOSE
!!     -------
!      ETHETA computes the coefficient Etheta in the flottability turbulent
!      flux. This coefficient relates the vertical flux of the virtual potential
!      temperature ( <Thv' W'> ) to the vertical flux of the conservative potential
!      temperature ( <Thl' W'> ).
!
!!**   METHOD
!!     ------
!!
!!     The virtual potential temperature perturbation is linearized in function
!!     of Thl' and Rnp'. The result is
!!        Thv'= ( ZA + ZC * Atheta * 2 * SRC ) Thl' 
!!             +( ZB + ZC * Amoist * 2 * SRC ) Rnp'
!!     From this relation, we can compute the vertical turbulent fluxes.
!!
!!     EXTERNAL
!!     --------
!!
!!        NONE
!!
!!     IMPLICIT ARGUMENTS
!!     ------------------
!!       Module MODD_CST : contains physical constants.
!!         XRV, XRD  : R for water vapor and dry air
!!            
!!     REFERENCE
!!     ---------
!!
!!
!!     AUTHOR
!!     ------
!!       Jean-Marie Carriere      * Meteo-France *
!!
!!     MODIFICATIONS
!!     -------------
!!       Original       20/03/95
!!     
!!       J. Stein       Feb  28, 1996   optimization + Doctorization
!!       J. Stein       Sept 15, 1996   Atheta previously computed
!!       J.-P. Pinty    May  20, 2003   Improve ETHETA expression
!!       J.L Redelsperger    03, 2021   Ocean Model Case 
!! ----------------------------------------------------------------------
!
!*       0. DECLARATIONS
!           ------------
USE MODD_CST, ONLY : CST_t
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
!
IMPLICIT NONE
!
!*       0.1 declarations of arguments and result
!
!
TYPE(DIMPHYEX_t),       INTENT(IN)   :: D
TYPE(CST_t),            INTENT(IN)   :: CST
INTEGER                              :: KRR          ! number of moist var.
INTEGER                              :: KRRI         ! number of ice var.
LOGICAL,                INTENT(IN)   ::  OOCEAN       ! switch for Ocean model version
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),  INTENT(IN)  ::   PTHLM     ! Conservative pot. temperature
REAL, DIMENSION(D%NIT,D%NJT,D%NKT,KRR), INTENT(IN) ::   PRM       ! Mixing ratios, where
!                                      PRM(:,:,:,1) = conservative mixing ratio
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),  INTENT(IN)  ::   PLOCPEXNM ! Lv(T)/Cp/Exner at time t-1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),  INTENT(IN)  ::   PATHETA   ! Atheta
!                                                    
LOGICAL,                INTENT(IN)   ::  OCOMPUTE_SRC ! flag to define dimensions of SIGS and
REAL, DIMENSION(MERGE(D%NIT,0,OCOMPUTE_SRC),&
                MERGE(D%NJT,0,OCOMPUTE_SRC),&
                MERGE(D%NKT,0,OCOMPUTE_SRC)),   INTENT(IN)  ::   PSRCM     ! Normalized 2dn_order
                                                    ! moment s'r'c/2Sigma_s2
!
REAL,DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT) :: PETHETA ! result
!
!
!
!*       0.2 declarations of local variables
!
REAL,DIMENSION(D%NIT,D%NJT,D%NKT) ::       &
                                        ZA, ZRW
!                ZA = coeft A, ZRW = total mixing ratio rw
REAL                                  :: ZDELTA  ! = Rv/Rd - 1
INTEGER                               :: JRR     ! moist loop counter
INTEGER                               :: JI,JJ,JK ! loop counter
INTEGER                               :: IIB,IJB,IIE,IJE,IKT
!
!---------------------------------------------------------------------------
!
!
!*       1. COMPUTE ETHETA
!           --------------
!
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('ETHETA',0,ZHOOK_HANDLE)
!
IIE=D%NIEC
IIB=D%NIBC
IJE=D%NJEC
IJB=D%NJBC
IKT=D%NKT
!
IF (OOCEAN) THEN                                    ! ocean case
  !$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:IKT)
  PETHETA(IIB:IIE,IJB:IJE,1:IKT) =  1.
  !$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:IKT)
ELSE   
 IF ( KRR == 0) THEN                                ! dry case
 !$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:IKT)
  PETHETA(IIB:IIE,IJB:IJE,1:IKT) = 1.
 !$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:IKT)
 ELSE IF ( KRR == 1 ) THEN                           ! only vapor
  ZDELTA = (CST%XRV/CST%XRD) - 1.
  !$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:IKT)
  PETHETA(IIB:IIE,IJB:IJE,1:IKT) = 1. + ZDELTA*PRM(IIB:IIE,IJB:IJE,1:IKT,1)
  !$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:IKT)
 ELSE                                                ! liquid water & ice present
  ZDELTA = (CST%XRV/CST%XRD) - 1.
  !$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:IKT)
  ZRW(IIB:IIE,IJB:IJE,1:IKT) = PRM(IIB:IIE,IJB:IJE,1:IKT,1)
  !$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:IKT)
!
  IF ( KRRI>0 ) THEN  ! rc and ri case
    !$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:IKT)
    ZRW(IIB:IIE,IJB:IJE,1:IKT) = ZRW(IIB:IIE,IJB:IJE,1:IKT) + PRM(IIB:IIE,IJB:IJE,1:IKT,3)
    !$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:IKT)
    DO JRR=5,KRR
      !$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:IKT)
      ZRW(IIB:IIE,IJB:IJE,1:IKT) = ZRW(IIB:IIE,IJB:IJE,1:IKT) + PRM(IIB:IIE,IJB:IJE,1:IKT,JRR)
      !$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:IKT)
    ENDDO
    !$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:IKT)
    ZA(IIB:IIE,IJB:IJE,1:IKT) = 1. + (                                    &  ! Compute A
              (1.+ZDELTA) * (PRM(IIB:IIE,IJB:IJE,1:IKT,1) - PRM(IIB:IIE,IJB:IJE,1:IKT,2) - PRM(IIB:IIE,IJB:IJE,1:IKT,4)) &
              -ZRW(IIB:IIE,IJB:IJE,1:IKT)                                                &
                     )  /  (1. + ZRW(IIB:IIE,IJB:IJE,1:IKT))
  !
  !   Etheta = ZA + ZC * Atheta  
  !   ZC is computed from line 2 to line 5
  !   - Atheta * 2. * SRC is computed at line 6 
  !
    PETHETA(IIB:IIE,IJB:IJE,1:IKT) = ZA(IIB:IIE,IJB:IJE,1:IKT)                                                 &
        +( PLOCPEXNM(IIB:IIE,IJB:IJE,1:IKT) * ZA(IIB:IIE,IJB:IJE,1:IKT)                                        &
               -(1.+ZDELTA) * (PTHLM(IIB:IIE,IJB:IJE,1:IKT) + PLOCPEXNM(IIB:IIE,IJB:IJE,1:IKT)*(               &
                                                    PRM(IIB:IIE,IJB:IJE,1:IKT,2)+PRM(IIB:IIE,IJB:IJE,1:IKT,4)))&
                            / (1. + ZRW(IIB:IIE,IJB:IJE,1:IKT))                                &
         ) * PATHETA(IIB:IIE,IJB:IJE,1:IKT) * 2. * PSRCM(IIB:IIE,IJB:IJE,1:IKT)
    !$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:IKT)
  ELSE
    DO JRR=3,KRR
      !$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:IKT)
      ZRW(IIB:IIE,IJB:IJE,1:IKT) = ZRW(IIB:IIE,IJB:IJE,1:IKT) + PRM(IIB:IIE,IJB:IJE,1:IKT,JRR)
      !$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:IKT)
    ENDDO
    !$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:IKT)
    ZA(IIB:IIE,IJB:IJE,1:IKT) = 1. + (                                    &  ! Compute A
              (1.+ZDELTA) * (PRM(IIB:IIE,IJB:IJE,1:IKT,1) - PRM(IIB:IIE,IJB:IJE,1:IKT,2)) &
              -ZRW(IIB:IIE,IJB:IJE,1:IKT)                                 &
                     )  /  (1. + ZRW(IIB:IIE,IJB:IJE,1:IKT))
  !
  !   Etheta = ZA + ZC * Atheta  
  !   ZC is computed from line 2 to line 5
  !   - Atheta * 2. * SRC is computed at line 6 
  !
    PETHETA(IIB:IIE,IJB:IJE,1:IKT) = ZA(IIB:IIE,IJB:IJE,1:IKT)                                                 &
        +( PLOCPEXNM(IIB:IIE,IJB:IJE,1:IKT) * ZA(IIB:IIE,IJB:IJE,1:IKT) -(1.+ZDELTA) * (PTHLM(IIB:IIE,IJB:IJE,1:IKT) &
        + PLOCPEXNM(IIB:IIE,IJB:IJE,1:IKT)*PRM(IIB:IIE,IJB:IJE,1:IKT,2))   &
         / (1. + ZRW(IIB:IIE,IJB:IJE,1:IKT))                                 &
         ) * PATHETA(IIB:IIE,IJB:IJE,1:IKT) * 2. * PSRCM(IIB:IIE,IJB:IJE,1:IKT)
    !$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:IKT)
  END IF
 END IF
!
END IF
!---------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('ETHETA',1,ZHOOK_HANDLE)
END SUBROUTINE ETHETA
END MODULE MODE_ETHETA
