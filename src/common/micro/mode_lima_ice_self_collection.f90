!MNH_LIC Copyright 2018-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-------------------------------------------------------------------------------
MODULE MODE_LIMA_ICE_SELF_COLLECTION
  IMPLICIT NONE
CONTAINS
! ####################################################################################
  SUBROUTINE LIMA_ICE_SELF_COLLECTION (LDCOMPUTE,                                    &
                                       PRHODREF, PT,                                 &
                                       PRIT, PCIT, PLBDI,                            &
                                       PRIT_SHAPE, PCIT_SHAPE, PLBDI_SHAPE,          &
                                       P_CI_ISC, P_SHCI_ISC, P_SHRI_ISCS, P_SHCI_ISCS)
! ####################################################################################
!
!!    PURPOSE
!!    -------
!!      Compute the self-collection of ice crystals
!!
!!
!!    AUTHOR
!!    ------
!!      M.    Taufour    * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original             15/03/2024
!!      C. Barthe   12/04/24   change in ice to snow conversion, simplify the code
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST,             ONLY : XTT
USE MODD_PARAM_LIMA,      ONLY : XRTMIN, XCTMIN, XCEXVT, NNB_CRYSTAL_SHAPE, LCRYSTAL_SHAPE
USE MODD_PARAM_LIMA_COLD, ONLY : NSCLBDAI, XSCINTP1I, XSCINTP2I, XKER_N_ISCI, XFNISCI, XCOLEXII, &
                                 XFISCS, XCOLII, XLBISCS1, XLBISCS2, XLBISCS3, &
                                 XLBNISCI1, XLBNISCI2, XKER_N_ISCI_SHAPE, XKER_R_ISCI_SHAPE, &
                                 XBI_SHAPE
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
LOGICAL, DIMENSION(:),INTENT(IN)    :: LDCOMPUTE
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRHODREF ! Reference Exner function
REAL, DIMENSION(:),   INTENT(IN)    :: PT       ! Temperature
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRIT    ! Ice mr at t
REAL, DIMENSION(:),   INTENT(IN)    :: PCIT    ! Ice C. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDI   !
!
REAL, DIMENSION(:,:), INTENT(IN)    :: PRIT_SHAPE    ! Ice mr at t
REAL, DIMENSION(:,:), INTENT(IN)    :: PCIT_SHAPE    ! Ice C. at t
REAL, DIMENSION(:,:), INTENT(IN)    :: PLBDI_SHAPE   !
!
REAL, DIMENSION(:),   INTENT(OUT)   :: P_CI_ISC 
REAL, DIMENSION(:,:), INTENT(OUT)   :: P_SHCI_ISC 
REAL, DIMENSION(:,:), INTENT(OUT)   :: P_SHCI_ISCS 
REAL, DIMENSION(:),   INTENT(OUT)   :: P_SHRI_ISCS
!
!*       0.2   Declarations of local variables :
!
LOGICAL, DIMENSION(SIZE(PT)) :: GISC
INTEGER :: IGISC, JJ, JSH, JSH2, JDUO
INTEGER, DIMENSION(:), ALLOCATABLE :: IVEC1, IVEC2        ! Vectors of indices
REAL,    DIMENSION(:), ALLOCATABLE :: ZVEC1, ZVEC2, ZVEC3 ! Work vectors
REAL,    DIMENSION(SIZE(PT)) :: ZONE
REAL,    DIMENSION(SIZE(PT)) :: ZW1, ZW2, ZW3, ZW4  ! work arrays
!
!-------------------------------------------------------------------------------
!
!*       1.     ICE CRYSTALS SELF-COLLECTION
!	        ----------------------------
!
ZONE(:) = 1.
P_CI_ISC(:)      = 0.
P_SHCI_ISC(:,:)  = 0.
P_SHRI_ISCS(:)   = 0.
P_SHCI_ISCS(:,:) = 0.
!
IF (LCRYSTAL_SHAPE) THEN
!
!
!*       1.1    Plate + Plate / Column + Column / Droxtal + Droxtal --> Bullet-rosette
!               ----------------------------------------------------------------------
!
  DO JSH = 1, NNB_CRYSTAL_SHAPE
    IF (JSH .NE. 3) THEN
      ZW1(:) = 0.
      ZW2(:) = 0.
      !
      GISC(:) = PCIT_SHAPE(:,JSH) > XCTMIN(4) .AND. PRIT_SHAPE(:,JSH) > XRTMIN(4)
      IGISC   = COUNT(GISC(:))
      !
      IF (IGISC > 0) THEN
!
!*       1.1.0  allocations
!
        ALLOCATE(ZVEC1(IGISC))
        ALLOCATE(IVEC1(IGISC))
        ALLOCATE(ZVEC3(IGISC))
!
!*       1.1.1  select the (ZLBDAI,ZLBDAI) couplet
!
        ZVEC1(:) = PACK( PLBDI_SHAPE(:,JSH),MASK=GISC(:) )
!
!*       1.1.2  find the next lower indice for the ZLBDAI and for the ZLBDAI
!               in the geometrical set of (Lbda_I,Lbda_I) couplet use to
!               tabulate the SACCS-kernel
!
        ZVEC1(1:IGISC) = MAX( 1.0001, MIN( FLOAT(NSCLBDAI)-0.0001,           &
                                           XSCINTP1I * LOG( ZVEC1(1:IGISC) ) + XSCINTP2I ) )
        IVEC1(1:IGISC) = INT( ZVEC1(1:IGISC) )
        ZVEC1(1:IGISC) = ZVEC1(1:IGISC) - FLOAT( IVEC1(1:IGISC) )
!
!*       1.1.3  perform the bilinear interpolation of the normalized ISCI-kernel
!
        DO JJ = 1, IGISC
          ZVEC3(JJ) =  (  XKER_N_ISCI_SHAPE(JSH,JSH,IVEC1(JJ)+1,IVEC1(JJ)+1) *  ZVEC1(JJ)          &
                        - XKER_N_ISCI_SHAPE(JSH,JSH,IVEC1(JJ)+1,IVEC1(JJ)  ) * (ZVEC1(JJ) - 1.0) ) &
                                                                             *  ZVEC1(JJ) &
                     - (  XKER_N_ISCI_SHAPE(JSH,JSH,IVEC1(JJ)  ,IVEC1(JJ)+1) *  ZVEC1(JJ)          &
                        - XKER_N_ISCI_SHAPE(JSH,JSH,IVEC1(JJ)  ,IVEC1(JJ)  ) * (ZVEC1(JJ) - 1.0) ) &
                                                                             * (ZVEC1(JJ) - 1.0)
        END DO
        ZW1(:) = UNPACK( VECTOR=ZVEC3(:),MASK=GISC(:),FIELD=0.0 ) !! 
        !
        WHERE( GISC(:) )
          ZW2(:) = ZW1(:) * XFNISCI * MIN(XCOLII * EXP( XCOLEXII*(PT(:)-XTT) ),ZONE(:)) * PCIT_SHAPE(:,JSH)**2 &  
                          * PRHODREF(:)**(-XCEXVT-1.) * (XLBNISCI1 + XLBNISCI2) / PLBDI_SHAPE(:,JSH)**2
          P_SHCI_ISC(:,JSH) = P_SHCI_ISC(:,JSH) - 2. * ZW2(:)
          P_SHCI_ISC(:,3)   = P_SHCI_ISC(:,3)   +      ZW2(:)       ! add to BURO
        END WHERE
        !
        DEALLOCATE(IVEC1)
        DEALLOCATE(ZVEC1)
        DEALLOCATE(ZVEC3)
      END IF
    END IF
  END DO  
!
!
!*       1.2    Plate + Column / Plate + Droxtal / Column + Droxtal --> Bullet-rosette
!               ----------------------------------------------------------------------
!
!*       1.2.1  select the two interacting pristine ice shapes
! 
  DO JDUO = 1, 3
    IF (JDUO == 1) THEN       ! plate + column --> bullet_rosette
      JSH  = 1 ! plates
      JSH2 = 2 ! columns
    ELSE IF (JDUO == 2) THEN  ! plate + droxtal --> bullet_rosette
      JSH  = 1 ! plates
      JSH2 = 4 ! droxtals
    ELSE                      ! column + droxtal --> bullet_rosette
      JSH  = 2 ! columns
      JSH2 = 4 ! droxtals
    END IF
    ! 
    ZW1(:) = 0.
    ZW2(:) = 0.
    !
    GISC(:) = PCIT_SHAPE(:,JSH)  > XCTMIN(4) .AND. PRIT_SHAPE(:,JSH)  > XRTMIN(4) .AND. & 
              PCIT_SHAPE(:,JSH2) > XCTMIN(4) .AND. PRIT_SHAPE(:,JSH2) > XRTMIN(4)
    IGISC = COUNT(GISC(:))
    !
    IF (IGISC > 0) THEN
!
!*       1.2.2  allocations
!
      ALLOCATE(ZVEC1(IGISC))
      ALLOCATE(IVEC1(IGISC))
      ALLOCATE(ZVEC2(IGISC))
      ALLOCATE(IVEC2(IGISC))       
      ALLOCATE(ZVEC3(IGISC))
!
!*       1.2.3  select the (ZLBDAI,ZLBDAI) couplet
!
      ZVEC1(:) = PACK( PLBDI_SHAPE(:,JSH),MASK=GISC(:) )
      ZVEC2(:) = PACK( PLBDI_SHAPE(:,JSH2),MASK=GISC(:) )
!
!*       1.2.4  find the next lower indice for the ZLBDAI and for the ZLBDAI
!               in the geometrical set of (Lbda_I,Lbda_I) couplet use to
!               tabulate the SACCS-kernel
!
      ZVEC1(1:IGISC) = MAX( 1.0001, MIN( FLOAT(NSCLBDAI)-0.0001,           &
                                         XSCINTP1I * LOG( ZVEC1(1:IGISC) ) + XSCINTP2I ) )
      IVEC1(1:IGISC) = INT( ZVEC1(1:IGISC) )
      ZVEC1(1:IGISC) = ZVEC1(1:IGISC) - FLOAT( IVEC1(1:IGISC) )
      !
      ZVEC2(1:IGISC) = MAX( 1.0001, MIN( FLOAT(NSCLBDAI)-0.0001,           &
                                         XSCINTP1I * LOG( ZVEC2(1:IGISC) ) + XSCINTP2I ) )
      IVEC2(1:IGISC) = INT( ZVEC2(1:IGISC) )
      ZVEC2(1:IGISC) = ZVEC2(1:IGISC) - FLOAT( IVEC2(1:IGISC) )
!
!*       1.2.5  perform the bilinear interpolation of the normalized ISCI-kernel
!
      DO JJ = 1, IGISC
        ZVEC3(JJ) =  (  XKER_N_ISCI_SHAPE(JSH,JSH2,IVEC1(JJ)+1,IVEC2(JJ)+1) *  ZVEC2(JJ)          &
                      - XKER_N_ISCI_SHAPE(JSH,JSH2,IVEC1(JJ)+1,IVEC2(JJ)  ) * (ZVEC2(JJ) - 1.0) ) &
                                                                            *  ZVEC1(JJ) &
                   - (  XKER_N_ISCI_SHAPE(JSH,JSH2,IVEC1(JJ)  ,IVEC2(JJ)+1) *  ZVEC2(JJ)          &
                      - XKER_N_ISCI_SHAPE(JSH,JSH2,IVEC1(JJ)  ,IVEC2(JJ)  ) * (ZVEC2(JJ) - 1.0) ) &
                                                                            * (ZVEC1(JJ) - 1.0)
      END DO
      ZW1(:) = UNPACK( VECTOR=ZVEC3(:),MASK=GISC(:),FIELD=0.0 )
      !
      WHERE( GISC(:) )
        ZW2(:) = ZW1(:) * XFNISCI * MIN(XCOLII * EXP( XCOLEXII*(PT(:)-XTT) ),ZONE(:)) * &
                          PCIT_SHAPE(:,JSH) * PCIT_SHAPE(:,JSH2) *                      &  
                          PRHODREF(:)**(-XCEXVT-1.) *                                   &
                         (XLBNISCI1 / (2. * PLBDI_SHAPE(:,JSH)**2) +                    &
                          XLBNISCI2 / (PLBDI_SHAPE(:,JSH) * PLBDI_SHAPE(:,JSH2) )  +    &
                          XLBNISCI1 / (2. * PLBDI_SHAPE(:,JSH2)**2)                     ) !++cb-- JSH --> JSH2
        !
        P_SHCI_ISC(:,JSH)  = P_SHCI_ISC(:,JSH)  - ZW2(:)
        P_SHCI_ISC(:,JSH2) = P_SHCI_ISC(:,JSH2) - ZW2(:) 
        P_SHCI_ISC(:,3)    = P_SHCI_ISC(:,3)    + ZW2(:)     ! add to BURO
      END WHERE
      !
      DEALLOCATE(IVEC1)
      DEALLOCATE(ZVEC1)
      DEALLOCATE(ZVEC2)
      DEALLOCATE(IVEC2)
      DEALLOCATE(ZVEC3)
    END IF     
  END DO
!
!
!*       2.     ICE CRYSTALS SELF-COLLECTION AND CONVERSION TO AGGREGATES
!	        ---------------------------------------------------------
!  
  JSH2 = 3
  DO JSH = 1, NNB_CRYSTAL_SHAPE
    ZW1(:) = 0.
    ZW2(:) = 0.
    !
    GISC(:) = PCIT_SHAPE(:,JSH)  > XCTMIN(4) .AND. PRIT_SHAPE(:,JSH)  > XRTMIN(4) .AND. & 
              PCIT_SHAPE(:,JSH2) > XCTMIN(4) .AND. PRIT_SHAPE(:,JSH2) > XRTMIN(4)
    IGISC   = COUNT(GISC(:))
    !
    IF (IGISC > 0) THEN
!
!*       2.0    allocations
!
      ALLOCATE(ZVEC1(IGISC))
      ALLOCATE(IVEC1(IGISC))
      ALLOCATE(ZVEC2(IGISC))
      ALLOCATE(IVEC2(IGISC))
      ALLOCATE(ZVEC3(IGISC))
!
!*       2.1    select the (ZLBDAI,ZLBDAI) couplet
!
      ZVEC1(:) = PACK( PLBDI_SHAPE(:,JSH),MASK=GISC(:) )
      ZVEC2(:) = PACK( PLBDI_SHAPE(:,JSH2),MASK=GISC(:) )
!
!*       2.2    find the next lower indice for the ZLBDAI and for the ZLBDAI
!               in the geometrical set of (Lbda_I,Lbda_I) couplet use to
!               tabulate the SACCS-kernel
!
      ZVEC1(1:IGISC) = MAX( 1.0001, MIN( FLOAT(NSCLBDAI)-0.0001,           &
                                         XSCINTP1I * LOG( ZVEC1(1:IGISC) ) + XSCINTP2I ) )
      IVEC1(1:IGISC) = INT( ZVEC1(1:IGISC) )
      ZVEC1(1:IGISC) = ZVEC1(1:IGISC) - FLOAT( IVEC1(1:IGISC) )
      !
      ZVEC2(1:IGISC) = MAX( 1.0001, MIN( FLOAT(NSCLBDAI)-0.0001,           &
                                         XSCINTP1I * LOG( ZVEC2(1:IGISC) ) + XSCINTP2I ) )
      IVEC2(1:IGISC) = INT( ZVEC2(1:IGISC) )
      ZVEC2(1:IGISC) = ZVEC2(1:IGISC) - FLOAT( IVEC2(1:IGISC) )
!
!*       2.3    perform the bilinear interpolation of the normalized ISCI-kernel
!
      DO JJ = 1, IGISC
        ZVEC3(JJ) = (  XKER_N_ISCI_SHAPE(JSH,JSH2,IVEC1(JJ)+1,IVEC2(JJ)+1) *  ZVEC2(JJ)          &
                     - XKER_N_ISCI_SHAPE(JSH,JSH2,IVEC1(JJ)+1,IVEC2(JJ)  ) * (ZVEC2(JJ) - 1.0) ) &
                                                                           *  ZVEC1(JJ)          &
                  - (  XKER_N_ISCI_SHAPE(JSH,JSH2,IVEC1(JJ)  ,IVEC2(JJ)+1) *  ZVEC2(JJ)          &
                     - XKER_N_ISCI_SHAPE(JSH,JSH2,IVEC1(JJ)  ,IVEC2(JJ)  ) * (ZVEC2(JJ) - 1.0) ) &
                                                                           * (ZVEC1(JJ) - 1.0)
      END DO
      ZW1(:) = UNPACK( VECTOR=ZVEC3(:),MASK=GISC(:),FIELD=0.0 ) !! 
      !
      DO JJ = 1, IGISC
        ZVEC3(JJ) = (  XKER_R_ISCI_SHAPE(JSH,IVEC1(JJ)+1,IVEC2(JJ)+1) *  ZVEC2(JJ)          &
                     - XKER_R_ISCI_SHAPE(JSH,IVEC1(JJ)+1,IVEC2(JJ)  ) * (ZVEC2(JJ) - 1.0) ) &
                                                                      *  ZVEC1(JJ)          &
                  - (  XKER_R_ISCI_SHAPE(JSH,IVEC1(JJ)  ,IVEC2(JJ)+1) *  ZVEC2(JJ)          &
                     - XKER_R_ISCI_SHAPE(JSH,IVEC1(JJ)  ,IVEC2(JJ)  ) * (ZVEC2(JJ) - 1.0) ) &
                                                                      * (ZVEC1(JJ) - 1.0)
      END DO
      ZW3(:) = UNPACK( VECTOR=ZVEC3(:),MASK=GISC(:),FIELD=0.0 ) !! 
      !
      WHERE( GISC(:) )
        ZW2(:) = ZW1(:) * XFNISCI * MIN( XCOLII*EXP( XCOLEXII*(PT(:)-XTT) ),ZONE(:) ) * &
                          PCIT_SHAPE(:,JSH) * PCIT_SHAPE(:,JSH2) *                      &  
                          PRHODREF(:)**(-XCEXVT-1.) *                                   &
                         (XLBNISCI1 / (2. * PLBDI_SHAPE(:,JSH)**2) +                    &
                          XLBNISCI2 / (PLBDI_SHAPE(:,JSH) * PLBDI_SHAPE(:,JSH2))  +     &
                          XLBNISCI1 / (2. * PLBDI_SHAPE(:,JSH2)**2))                       !++cb-- JSH --> JSH2
        ZW4(:) = ZW3(:) * XFISCS(JSH) * MIN( XCOLII*EXP( XCOLEXII*(PT(:)-XTT) ),ZONE(:)) * &
                          PCIT_SHAPE(:,JSH) * PCIT_SHAPE(:,JSH2) *                         &
                          PRHODREF(:)**(-XCEXVT-1.) *                                      &
!++cb++ 15/04/23 Y=BUR (JSH2 = 3), X=PLA/COL/DRO/BUR (JSH)
!                         (XLBISCS1(JSH) / (PLBDI_SHAPE(:,JSH)**2) +                        &
!                          XLBISCS2(JSH) / (PLBDI_SHAPE(:,JSH) * PLBDI_SHAPE(:,JSH2))  +    &
!                          XLBISCS3(JSH) / (PLBDI_SHAPE(:,JSH2)**2))
                         (XLBISCS1(JSH2) / (PLBDI_SHAPE(:,JSH)**2 * PLBDI_SHAPE(:,JSH2)**XBI_SHAPE(JSH2)) + &
                          XLBISCS2(JSH2) / (PLBDI_SHAPE(:,JSH)    * PLBDI_SHAPE(:,JSH2)**(1.+XBI_SHAPE(JSH2))) + &
                          XLBISCS3(JSH2) / (                        PLBDI_SHAPE(:,JSH2)**(2.+XBI_SHAPE(JSH2))))
!--cb--
        !
        P_SHCI_ISCS(:,JSH)  = P_SHCI_ISCS(:,JSH)  - ZW2(:)
        P_SHCI_ISCS(:,JSH2) = P_SHCI_ISCS(:,JSH2) - ZW2(:) 
        P_SHRI_ISCS(:)      = P_SHRI_ISCS(:)      - ZW4(:)     ! add to SNOW
      END WHERE
      !
      DEALLOCATE(IVEC1)
      DEALLOCATE(ZVEC1)
      DEALLOCATE(ZVEC2)
      DEALLOCATE(IVEC2) 
      DEALLOCATE(ZVEC3)
    END IF
  END DO
!
ELSE
!
!
!*       3.     ICE CRYSTALS SELF-COLLECTION WHEN ICE SHAPES NOT CONSIDERED
!	        -----------------------------------------------------------
!
  ZW1(:) = 0.
  ZW2(:) = 0.
  !
  GISC(:) = PCIT(:) > XCTMIN(4) .AND. PRIT(:) > XRTMIN(4)
  IGISC   = COUNT(GISC(:))
  !
  IF (IGISC > 0) THEN
!
!*       3.0    allocations
!
    ALLOCATE(ZVEC1(IGISC))
    ALLOCATE(IVEC1(IGISC))
    ALLOCATE(ZVEC3(IGISC))
!
!*       3.1    select the (ZLBDAI,ZLBDAI) couplet
!
    ZVEC1(:) = PACK( PLBDI(:),MASK=GISC(:) )
!
!*       3.2    find the next lower indice for the ZLBDAI and for the ZLBDAI
!               in the geometrical set of (Lbda_I,Lbda_I) couplet use to
!               tabulate the SACCS-kernel
!
    ZVEC1(1:IGISC) = MAX( 1.0001, MIN( FLOAT(NSCLBDAI)-0.0001,           &
                                       XSCINTP1I * LOG( ZVEC1(1:IGISC) ) + XSCINTP2I ) )
    IVEC1(1:IGISC) = INT( ZVEC1(1:IGISC) )
    ZVEC1(1:IGISC) = ZVEC1(1:IGISC) - FLOAT( IVEC1(1:IGISC) )
!
!*       3.3    perform the bilinear interpolation of the normalized ISCI-kernel
!
    DO JJ = 1, IGISC
      ZVEC3(JJ) = (  XKER_N_ISCI(IVEC1(JJ)+1,IVEC1(JJ)+1) *  ZVEC1(JJ)          &
                   - XKER_N_ISCI(IVEC1(JJ)+1,IVEC1(JJ)  ) * (ZVEC1(JJ) - 1.0) ) &
                                                          * ZVEC1(JJ)           &
                - (  XKER_N_ISCI(IVEC1(JJ)  ,IVEC1(JJ)+1) *  ZVEC1(JJ)          &
                   - XKER_N_ISCI(IVEC1(JJ)  ,IVEC1(JJ)  ) * (ZVEC1(JJ) - 1.0) ) &
                                                          * (ZVEC1(JJ) - 1.0)
    END DO
    ZW1(:) = UNPACK( VECTOR=ZVEC3(:),MASK=GISC(:),FIELD=0.0 ) !! 
    !
    WHERE( GISC(:) )
      P_CI_ISC(:) = - ZW1(:) * XFNISCI * MIN(XCOLII * EXP( XCOLEXII*(PT(:)-XTT) ),ZONE(:)) * PCIT(:)**2 &  
                             * PRHODREF(:)**(-XCEXVT-1.) * (XLBNISCI1 + XLBNISCI2) / PLBDI(:)**2
    END WHERE
    DEALLOCATE(IVEC1)
    DEALLOCATE(ZVEC1)
    DEALLOCATE(ZVEC3)
  END IF
END IF
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_ICE_SELF_COLLECTION
END MODULE MODE_LIMA_ICE_SELF_COLLECTION
