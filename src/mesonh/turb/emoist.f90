!MNH_LIC Copyright 1995-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!#################
MODULE MODI_EMOIST
!#################
!
INTERFACE
!
FUNCTION EMOIST(KRR,KRRI,PTHLM,PRM,PLOCPEXNM,PAMOIST,PSRCM) RESULT(PEMOIST)
!
INTEGER                              :: KRR        ! number of moist var.
INTEGER                              :: KRRI       ! number of ice var.
!
REAL, DIMENSION(:,:,:),  INTENT(IN)  ::   PTHLM    ! Conservative pot. temperature
REAL, DIMENSION(:,:,:,:), INTENT(IN) ::   PRM      ! Mixing ratios, where
!                                    PRM(:,:,:,1) = conservative mixing ratio
REAL, DIMENSION(:,:,:),  INTENT(IN)  ::   PLOCPEXNM ! Lv(T)/Cp/Exner at time t-1
REAL, DIMENSION(:,:,:),  INTENT(IN)  ::   PAMOIST   ! Amoist
REAL, DIMENSION(:,:,:),  INTENT(IN)  ::   PSRCM     ! Normalized 2dn_order
                                                    ! moment s'r'c/2Sigma_s2
!
REAL,DIMENSION(SIZE(PTHLM,1),SIZE(PTHLM,2),SIZE(PTHLM,3)):: PEMOIST ! result
!
END FUNCTION EMOIST
!
END INTERFACE
!
END MODULE MODI_EMOIST
!
!   ############################################################################
FUNCTION EMOIST(KRR,KRRI,PTHLM,PRM,PLOCPEXNM,PAMOIST,PSRCM) RESULT(PEMOIST)
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
USE MODD_CST
USE MODD_DYN_n, ONLY : LOCEAN
!
IMPLICIT NONE
!
!*       0.1 declarations of arguments and result
!
!
INTEGER                              :: KRR        ! number of moist var.
INTEGER                              :: KRRI       ! number of ice var.
!
REAL, DIMENSION(:,:,:),  INTENT(IN)  ::   PTHLM    ! Conservative pot. temperature
REAL, DIMENSION(:,:,:,:), INTENT(IN) ::   PRM      ! Mixing ratios, where
!                                    PRM(:,:,:,1) = conservative mixing ratio
REAL, DIMENSION(:,:,:),  INTENT(IN)  ::   PLOCPEXNM ! Lv(T)/Cp/Exner at time t-1
REAL, DIMENSION(:,:,:),  INTENT(IN)  ::   PAMOIST   ! Amoist
REAL, DIMENSION(:,:,:),  INTENT(IN)  ::   PSRCM     ! Normalized 2dn_order
                                                    ! moment s'r'c/2Sigma_s2
!
REAL,DIMENSION(SIZE(PTHLM,1),SIZE(PTHLM,2),SIZE(PTHLM,3)):: PEMOIST ! result
!
!*       0.2 declarations of local variables
!
REAL,DIMENSION(SIZE(PTHLM,1),SIZE(PTHLM,2),SIZE(PTHLM,3)) ::       &
                                        ZA, ZRW
!                ZA = coeft A, ZRW = total mixing ratio rw 
REAL                                  :: ZDELTA  ! = Rv/Rd - 1
INTEGER                               :: JRR     ! moist loop counter
!
!---------------------------------------------------------------------------
!
!
!*       1. COMPUTE EMOIST
!           --------------
IF (LOCEAN) THEN
 IF ( KRR == 0 ) THEN                                ! Unsalted
   PEMOIST(:,:,:) = 0.
 ELSE
   PEMOIST(:,:,:) = 1.                              ! Salted case
 END IF
!
ELSE
!
 IF ( KRR == 0 ) THEN                                ! dry case
   PEMOIST(:,:,:) = 0.
 ELSE IF ( KRR == 1 ) THEN                           ! only vapor
  ZDELTA = (XRV/XRD) - 1.
  PEMOIST(:,:,:) = ZDELTA*PTHLM(:,:,:)
 ELSE                                                ! liquid water & ice present
  ZDELTA = (XRV/XRD) - 1.
  ZRW(:,:,:) = PRM(:,:,:,1)
!
  IF ( KRRI>0) THEN  ! rc and ri case  
    ZRW(:,:,:) = ZRW(:,:,:) + PRM(:,:,:,3)
    DO JRR=5,KRR
      ZRW(:,:,:) = ZRW(:,:,:) + PRM(:,:,:,JRR)
    ENDDO
    ZA(:,:,:) = 1. + (                                    &  ! Compute A
              (1.+ZDELTA) * (PRM(:,:,:,1) - PRM(:,:,:,2) - PRM(:,:,:,4)) &
              -ZRW(:,:,:)                                                &
                     )  /  (1. + ZRW(:,:,:)) 
  !
  !   Emoist = ZB + ZC * Amoist
  !   ZB is computed from line 1 to line 2
  !   ZC is computed from line 3 to line 5
  !   Amoist* 2 * SRC is computed at line 6
  !
    PEMOIST(:,:,:) = ZDELTA * (PTHLM(:,:,:) + PLOCPEXNM(:,:,:)*(               &
                                                    PRM(:,:,:,2)+PRM(:,:,:,4)))&
                            / (1. + ZRW(:,:,:))                                &
        +( PLOCPEXNM(:,:,:) * ZA(:,:,:)                                        &
               -(1.+ZDELTA) * (PTHLM(:,:,:) + PLOCPEXNM(:,:,:)*(               &
                                                    PRM(:,:,:,2)+PRM(:,:,:,4)))&
                            / (1. + ZRW(:,:,:))                                &
         ) * PAMOIST(:,:,:) * 2. * PSRCM(:,:,:)
  ELSE
    DO JRR=3,KRR
      ZRW(:,:,:) = ZRW(:,:,:) + PRM(:,:,:,JRR)
    ENDDO
    ZA(:,:,:) = 1. + (                                    &  ! Compute ZA
              (1.+ZDELTA) * (PRM(:,:,:,1) - PRM(:,:,:,2)) &
              -ZRW(:,:,:)                                 &
                     )  /  (1. + ZRW(:,:,:)) 
  !
  !   Emoist = ZB + ZC * Amoist
  !   ZB is computed from line 1 to line 2
  !   ZC is computed from line 3 to line 5
  !   Amoist* 2 * SRC is computed at line 6
  !
    PEMOIST(:,:,:) = ZDELTA * (PTHLM(:,:,:) + PLOCPEXNM(:,:,:)*PRM(:,:,:,2))   &
                            / (1. + ZRW(:,:,:))                                &
        +( PLOCPEXNM(:,:,:) * ZA(:,:,:)                                        &
               -(1.+ZDELTA) * (PTHLM(:,:,:) + PLOCPEXNM(:,:,:)*PRM(:,:,:,2))   &
                            / (1. + ZRW(:,:,:))                                &
         ) * PAMOIST(:,:,:) * 2. * PSRCM(:,:,:)
  END IF
 END IF
!
END IF
!---------------------------------------------------------------------------
!
END FUNCTION EMOIST
