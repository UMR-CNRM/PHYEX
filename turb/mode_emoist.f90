!MNH_LIC Copyright 1995-2022 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
MODULE MODE_EMOIST
IMPLICIT NONE
CONTAINS
SUBROUTINE EMOIST(D,CST,KRR,KRRI,PTHLM,PRM,PLOCPEXNM,PAMOIST,PSRCM,OOCEAN,PEMOIST)
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!   ############################################################################
!
!      PURPOSE
!!     -------
!      EMOIST computes the coefficient Emoist in the flottability turbulent
!      flux. This coefficient relates the vertical flux of the virtual potential
!      temperature ( <Thv' W'> ) to the vertical flux of the conservative mixing
!      ratio ( <Rnp' W'> ). 
!
!!**   METHOD
!!     ------
!!     The virtual potential temperature perturbation is linearized in function
!!     of Thl' and Rnp'. The result is
!!        Thv'= ( ZA + ZC * Atheta * 2 * SRC ) Thl' 
!!             +( ZB + ZC * Amoist * 2 * SRC ) Rnp'
!!     From this relation, we can compute the verical turbulent fluxes.
!! 
!!     EXTERNAL
!!     --------
!!       NONE
!!
!!     IMPLICIT ARGUMENTS
!!     ------------------
!!       Module MODD_CST : contains physical constants.
!!         XRV, XRD  : R for water vapor and dry air
!!   
!!     REFERENCE
!!     ---------
!! 
!!       NONE
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
!!       J. Stein       Spet 15, 1996   Amoist previously computed
!!       J.-P. Pinty    May  20, 2003   Improve EMOIST expression
!!                  03/2021 (JL Redelsperger) Ocean model case 
!! 
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
INTEGER                              :: KRR        ! number of moist var.
INTEGER                              :: KRRI       ! number of ice var.
LOGICAL,                INTENT(IN)   ::  OOCEAN       ! switch for Ocean model version
!
REAL, DIMENSION(D%NIJT,D%NKT),  INTENT(IN)  ::   PTHLM    ! Conservative pot. temperature
REAL, DIMENSION(D%NIJT,D%NKT,KRR), INTENT(IN) ::   PRM      ! Mixing ratios, where
!                                    PRM(:,:,:,1) = conservative mixing ratio
REAL, DIMENSION(D%NIJT,D%NKT),  INTENT(IN)  ::   PLOCPEXNM ! Lv(T)/Cp/Exner at time t-1
REAL, DIMENSION(D%NIJT,D%NKT),  INTENT(IN)  ::   PAMOIST   ! Amoist
REAL, DIMENSION(D%NIJT,D%NKT),  INTENT(IN)  ::   PSRCM     ! Normalized 2dn_order
                                                    ! moment s'r'c/2Sigma_s2
!
REAL,DIMENSION(D%NIJT,D%NKT), INTENT(OUT) :: PEMOIST ! result
!
!*       0.2 declarations of local variables
!
REAL,DIMENSION(D%NIJT,D%NKT) ::       &
                                        ZA, ZRW
!                ZA = coeft A, ZRW = total mixing ratio rw 
REAL                                  :: ZDELTA  ! = Rv/Rd - 1
INTEGER                               :: JRR     ! moist loop counter
INTEGER                               :: JIJ,JK ! loop counter
INTEGER                               :: IIJB,IIJE,IKT
!
!---------------------------------------------------------------------------
!
!
!*       1. COMPUTE EMOIST
!           --------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('EMOIST',0,ZHOOK_HANDLE)
!
IIJB=D%NIJB
IIJE=D%NIJE
IKT=D%NKT
!
IF (OOCEAN) THEN
 IF ( KRR == 0 ) THEN                                ! Unsalted
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
   PEMOIST(IIJB:IIJE,:) = 0.
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
 ELSE
   !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
   PEMOIST(IIJB:IIJE,:) = 1.                              ! Salted case
   !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
 END IF
!
ELSE
!
 IF ( KRR == 0 ) THEN                                ! dry case
   PEMOIST(IIJB:IIJE,:) = 0.
 ELSE IF ( KRR == 1 ) THEN                           ! only vapor
  ZDELTA = (CST%XRV/CST%XRD) - 1.
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  PEMOIST(IIJB:IIJE,:) = ZDELTA*PTHLM(IIJB:IIJE,:)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
 ELSE                                                ! liquid water & ice present
  ZDELTA = (CST%XRV/CST%XRD) - 1.
  ZRW(IIJB:IIJE,:) = PRM(IIJB:IIJE,:,1)
!
  IF ( KRRI>0) THEN  ! rc and ri case
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ZRW(IIJB:IIJE,:) = ZRW(IIJB:IIJE,:) + PRM(IIJB:IIJE,:,3)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    DO JRR=5,KRR
      !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
      ZRW(IIJB:IIJE,:) = ZRW(IIJB:IIJE,:) + PRM(IIJB:IIJE,:,JRR)
      !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ENDDO
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ZA(IIJB:IIJE,:) = 1. + (                                    &  ! Compute A
              (1.+ZDELTA) * (PRM(IIJB:IIJE,:,1) - PRM(IIJB:IIJE,:,2) - PRM(IIJB:IIJE,:,4)) &
              -ZRW(IIJB:IIJE,:)                                                &
                     )  /  (1. + ZRW(IIJB:IIJE,:)) 
  !
  !   Emoist = ZB + ZC * Amoist
  !   ZB is computed from line 1 to line 2
  !   ZC is computed from line 3 to line 5
  !   Amoist* 2 * SRC is computed at line 6
  !
    PEMOIST(IIJB:IIJE,:) = ZDELTA * (PTHLM(IIJB:IIJE,:) + PLOCPEXNM(IIJB:IIJE,:)*(           &
                                                    PRM(IIJB:IIJE,:,2)+PRM(IIJB:IIJE,:,4)))&
                            / (1. + ZRW(IIJB:IIJE,:))                                &
        +( PLOCPEXNM(IIJB:IIJE,:) * ZA(IIJB:IIJE,:)                                        &
               -(1.+ZDELTA) * (PTHLM(IIJB:IIJE,:) + PLOCPEXNM(IIJB:IIJE,:)*(               &
                                                    PRM(IIJB:IIJE,:,2)+PRM(IIJB:IIJE,:,4)))&
                            / (1. + ZRW(IIJB:IIJE,:))                                &
         ) * PAMOIST(IIJB:IIJE,:) * 2. * PSRCM(IIJB:IIJE,:)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  ELSE
    DO JRR=3,KRR
      !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
      ZRW(IIJB:IIJE,:) = ZRW(IIJB:IIJE,:) + PRM(IIJB:IIJE,:,JRR)
      !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ENDDO
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ZA(IIJB:IIJE,:) = 1. + (                                    &  ! Compute ZA
              (1.+ZDELTA) * (PRM(IIJB:IIJE,:,1) - PRM(IIJB:IIJE,:,2)) &
              -ZRW(IIJB:IIJE,:)                                 &
                     )  /  (1. + ZRW(IIJB:IIJE,:)) 
  !
  !   Emoist = ZB + ZC * Amoist
  !   ZB is computed from line 1 to line 2
  !   ZC is computed from line 3 to line 5
  !   Amoist* 2 * SRC is computed at line 6
  !
    PEMOIST(IIJB:IIJE,:) = ZDELTA * (PTHLM(IIJB:IIJE,:) + PLOCPEXNM(IIJB:IIJE,:)* &
                                       PRM(IIJB:IIJE,:,2)) / (1. + ZRW(IIJB:IIJE,:))          &
        +( PLOCPEXNM(IIJB:IIJE,:) * ZA(IIJB:IIJE,:)                                        &
               -(1.+ZDELTA) * (PTHLM(IIJB:IIJE,:) + PLOCPEXNM(IIJB:IIJE,:)* &
               PRM(IIJB:IIJE,:,2)) / (1. + ZRW(IIJB:IIJE,:))                                &
         ) * PAMOIST(IIJB:IIJE,:) * 2. * PSRCM(IIJB:IIJE,:)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  END IF
 END IF
!
END IF
!---------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('EMOIST',1,ZHOOK_HANDLE)
END SUBROUTINE EMOIST
END MODULE MODE_EMOIST
