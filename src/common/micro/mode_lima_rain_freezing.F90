!MNH_LIC Copyright 2018-2024 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_LIMA_RAIN_FREEZING
  IMPLICIT NONE
CONTAINS
!     #######################################################################################
  SUBROUTINE LIMA_RAIN_FREEZING (CST, LIMAP, LIMAM, KSIZE, ODCOMPUTE,                                      &
                                 PRHODREF, PT, PLVFACT, PLSFACT,                        &
                                 PRRT, PCRT, PRIT, PCIT, PLBDR,                         &
                                 P_TH_CFRZ, P_RR_CFRZ, P_CR_CFRZ, P_RI_CFRZ, P_CI_CFRZ  )
!     #######################################################################################
!
!!    PURPOSE
!!    -------
!!      Compute the rain freezing by contact with an ice crystal
!!
!!
!!    AUTHOR
!!    ------
!!      J.-M. Cohard     * Laboratoire d'Aerologie*
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!      S.    Berthet    * Laboratoire d'Aerologie*
!!      B.    Vié        * CNRM *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original             15/03/2018
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAM_LIMA_MIXED, ONLY:PARAM_LIMA_MIXED_T
USE MODD_PARAM_LIMA, ONLY:PARAM_LIMA_T
USE MODD_CST, ONLY:CST_T
USE YOMHOOK, ONLY:LHOOK, DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(PARAM_LIMA_MIXED_T),INTENT(IN)::LIMAM
TYPE(PARAM_LIMA_T),INTENT(IN)::LIMAP
TYPE(CST_T),INTENT(IN)::CST
INTEGER, INTENT(IN) :: KSIZE
LOGICAL, DIMENSION(KSIZE),INTENT(IN)    :: ODCOMPUTE
!
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRHODREF ! Reference Exner function
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PT       !
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLVFACT  !
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLSFACT  !
!
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRRT    ! 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PCRT    ! 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRIT    ! 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PCIT    ! 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLBDR   ! 
!
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_TH_CFRZ
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_RR_CFRZ
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_CR_CFRZ
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_RI_CFRZ
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_CI_CFRZ
!
!*       0.2   Declarations of local variables :
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
REAL, DIMENSION(SIZE(PRRT)) :: ZW0, ZW1, ZW2, ZW3 ! work arrays
!
!-------------------------------------------------------------------------------
!
!
!*       1.     PRELIMINARY COMPUTATIONS
!               ------------------------
!
!
IF (LHOOK) CALL DR_HOOK('LIMA_RAIN_FREEZING', 0, ZHOOK_HANDLE)
P_TH_CFRZ(:)=0.
P_RR_CFRZ(:)=0.
P_CR_CFRZ(:)=0.
P_RI_CFRZ(:)=0.
P_CI_CFRZ(:)=0.
!
ZW0(:)=0.
ZW1(:)=0.
ZW2(:)=0.
ZW3(:)=0.
!
WHERE( PRIT(:)>LIMAP%XRTMIN(4) .AND. PRRT(:)>LIMAP%XRTMIN(3) .AND. PT(:)<CST%XTT .AND. &
       PCIT(:)>LIMAP%XCTMIN(4) .AND. PCRT(:)>LIMAP%XCTMIN(3) .AND. ODCOMPUTE(:) )
!
   ZW1(:) = LIMAM%XICFRR * PRIT(:) * PCRT(:)                    & ! RICFRRG
                                     * PLBDR(:)**LIMAM%XEXICFRR         &
                                     * PRHODREF(:)**(-LIMAP%XCEXVT-1.0)
!
   ZW2(:) = LIMAM%XRCFRI * PCIT(:) * PCRT(:)                    & ! RRCFRIG
                                     * PLBDR(:)**LIMAM%XEXRCFRI         &
                                     * PRHODREF(:)**(-LIMAP%XCEXVT-1.0)
!
! Comparison between heat to be released (to freeze rain) and heat sink (rain and ice temperature change)
! ZW0 is the proportion of process that can take place
   ZW0(:) = MAX(0., MIN(1., (ZW1(:)*CST%XCI+ZW2(:)*CST%XCL)*(CST%XTT-PT(:)) / &
        MAX(1.E-20, CST%XLVTT*ZW2(:))))
   ZW2(:) = ZW0(:) * ZW2(:)       ! Part of rain that can be frozen
   ZW3(:) = (1.- ZW0(:)) * ZW1(:) ! Part of collected pristine ice converted to rain
   ZW1(:) = ZW0(:) * ZW1(:)       ! Part of collected pristine ice that leads to graupel
!   
   P_RR_CFRZ(:) = - ZW2(:) + ZW3(:)
   P_CR_CFRZ(:) =   P_RR_CFRZ(:) * (PCRT(:)/PRRT(:))
   P_RI_CFRZ(:) = - ZW1(:) - ZW3(:)
   P_CI_CFRZ(:) =   P_RI_CFRZ(:) * (PCIT(:)/PRIT(:))
   P_TH_CFRZ(:) = - P_RR_CFRZ(:) * (PLSFACT(:)-PLVFACT(:))
!
END WHERE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('LIMA_RAIN_FREEZING', 1, ZHOOK_HANDLE)
END SUBROUTINE LIMA_RAIN_FREEZING
END MODULE MODE_LIMA_RAIN_FREEZING
