!MNH_LIC Copyright 2018-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-------------------------------------------------------------------------------
MODULE MODE_LIMA_ICE_SELF_COLLECTION
  IMPLICIT NONE
CONTAINS
! ####################################################################################
  SUBROUTINE LIMA_ICE_SELF_COLLECTION (CST, LIMAP, LIMAC, KSIZE, ODCOMPUTE,          &
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
USE MODD_PARAM_LIMA_COLD, ONLY:PARAM_LIMA_COLD_T
USE MODD_PARAM_LIMA, ONLY:PARAM_LIMA_T
USE MODD_CST, ONLY:CST_T
USE YOMHOOK, ONLY:LHOOK, DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(PARAM_LIMA_COLD_T),INTENT(IN)::LIMAC
TYPE(PARAM_LIMA_T),INTENT(IN)::LIMAP
TYPE(CST_T),INTENT(IN)::CST
INTEGER, INTENT(IN) :: KSIZE
LOGICAL, DIMENSION(KSIZE),INTENT(IN)    :: ODCOMPUTE
!
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRHODREF ! Reference Exner function
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PT       ! Temperature
!
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRIT    ! Ice mr at t
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PCIT    ! Ice C. at t
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLBDI   !
!
REAL, DIMENSION(KSIZE,LIMAP%NNB_CRYSTAL_SHAPE), INTENT(IN)    :: PRIT_SHAPE    ! Ice mr at t
REAL, DIMENSION(KSIZE,LIMAP%NNB_CRYSTAL_SHAPE), INTENT(IN)    :: PCIT_SHAPE    ! Ice C. at t
REAL, DIMENSION(KSIZE,LIMAP%NNB_CRYSTAL_SHAPE), INTENT(IN)    :: PLBDI_SHAPE   !
!
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_CI_ISC 
REAL, DIMENSION(KSIZE,LIMAP%NNB_CRYSTAL_SHAPE), INTENT(OUT)   :: P_SHCI_ISC 
REAL, DIMENSION(KSIZE,LIMAP%NNB_CRYSTAL_SHAPE), INTENT(OUT)   :: P_SHCI_ISCS 
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_SHRI_ISCS
!
!*       0.2   Declarations of local variables :
!
LOGICAL, DIMENSION(SIZE(PT)) :: GISC
INTEGER :: IGISC, IJ, ISH, ISH2, IDUO
INTEGER, DIMENSION(:), ALLOCATABLE :: IVEC1, IVEC2        ! Vectors of indices
REAL,    DIMENSION(:), ALLOCATABLE :: ZVEC1, ZVEC2, ZVEC3 ! Work vectors
REAL,    DIMENSION(SIZE(PT)) :: ZONE
REAL,    DIMENSION(SIZE(PT)) :: ZW1, ZW2, ZW3, ZW4  ! work arrays
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!*       1.     ICE CRYSTALS SELF-COLLECTION
!	        ----------------------------
!
IF (LHOOK) CALL DR_HOOK('LIMA_ICE_SELF_COLLECTION', 0, ZHOOK_HANDLE)
ZONE(:) = 1.
P_CI_ISC(:)      = 0.
P_SHCI_ISC(:,:)  = 0.
P_SHRI_ISCS(:)   = 0.
P_SHCI_ISCS(:,:) = 0.
!
IF (LIMAP%LCRYSTAL_SHAPE) THEN
!
!
!*       1.1    Plate + Plate / Column + Column / Droxtal + Droxtal --> Bullet-rosette
!               ----------------------------------------------------------------------
!
  DO ISH = 1, LIMAP%NNB_CRYSTAL_SHAPE
    IF (ISH .NE. 3) THEN
      ZW1(:) = 0.
      ZW2(:) = 0.
      !
      GISC(:) = PCIT_SHAPE(:,ISH) > LIMAP%XCTMIN(4) .AND. PRIT_SHAPE(:,ISH) > LIMAP%XRTMIN(4)
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
        ZVEC1(:) = PACK( PLBDI_SHAPE(:,ISH),MASK=GISC(:) )
!
!*       1.1.2  find the next lower indice for the ZLBDAI and for the ZLBDAI
!               in the geometrical set of (Lbda_I,Lbda_I) couplet use to
!               tabulate the SACCS-kernel
!
        ZVEC1(1:IGISC) = MAX( 1.0001, MIN( FLOAT(LIMAC%NSCLBDAI)-0.0001,           &
                                           LIMAC%XSCINTP1I * LOG( ZVEC1(1:IGISC) ) + LIMAC%XSCINTP2I ) )
        IVEC1(1:IGISC) = INT( ZVEC1(1:IGISC) )
        ZVEC1(1:IGISC) = ZVEC1(1:IGISC) - FLOAT( IVEC1(1:IGISC) )
!
!*       1.1.3  perform the bilinear interpolation of the normalized ISCI-kernel
!
        DO IJ = 1, IGISC
          ZVEC3(IJ) =  (  LIMAC%XKER_N_ISCI_SHAPE(ISH,ISH,IVEC1(IJ)+1,IVEC1(IJ)+1) *  ZVEC1(IJ)          &
                        - LIMAC%XKER_N_ISCI_SHAPE(ISH,ISH,IVEC1(IJ)+1,IVEC1(IJ)  ) * (ZVEC1(IJ) - 1.0) ) &
                                                                             *  ZVEC1(IJ) &
                     - (  LIMAC%XKER_N_ISCI_SHAPE(ISH,ISH,IVEC1(IJ)  ,IVEC1(IJ)+1) *  ZVEC1(IJ)          &
                        - LIMAC%XKER_N_ISCI_SHAPE(ISH,ISH,IVEC1(IJ)  ,IVEC1(IJ)  ) * (ZVEC1(IJ) - 1.0) ) &
                                                                             * (ZVEC1(IJ) - 1.0)
        END DO
        ZW1(:) = UNPACK( VECTOR=ZVEC3(:),MASK=GISC(:),FIELD=0.0 ) !! 
        !
        WHERE( GISC(:) )
           ZW2(:) = ZW1(:) * LIMAC%XFNISCI * &
                MIN(LIMAC%XCOLII * EXP( LIMAC%XCOLEXII*(PT(:)-CST%XTT) ), ZONE(:) ) * &
                PCIT_SHAPE(:,ISH)**2 * PRHODREF(:)**(-LIMAP%XCEXVT-1.) * &
                (LIMAC%XLBNISCI1 + LIMAC%XLBNISCI2) / PLBDI_SHAPE(:,ISH)**2
          P_SHCI_ISC(:,ISH) = P_SHCI_ISC(:,ISH) - 2. * ZW2(:)
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
  DO IDUO = 1, 3
    IF (IDUO == 1) THEN       ! plate + column --> bullet_rosette
      ISH  = 1 ! plates
      ISH2 = 2 ! columns
    ELSE IF (IDUO == 2) THEN  ! plate + droxtal --> bullet_rosette
      ISH  = 1 ! plates
      ISH2 = 4 ! droxtals
    ELSE                      ! column + droxtal --> bullet_rosette
      ISH  = 2 ! columns
      ISH2 = 4 ! droxtals
    END IF
    ! 
    ZW1(:) = 0.
    ZW2(:) = 0.
    !
    GISC(:) = PCIT_SHAPE(:,ISH)  > LIMAP%XCTMIN(4) .AND. PRIT_SHAPE(:,ISH)  > LIMAP%XRTMIN(4) .AND. & 
              PCIT_SHAPE(:,ISH2) > LIMAP%XCTMIN(4) .AND. PRIT_SHAPE(:,ISH2) > LIMAP%XRTMIN(4)
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
      ZVEC1(:) = PACK( PLBDI_SHAPE(:,ISH),MASK=GISC(:) )
      ZVEC2(:) = PACK( PLBDI_SHAPE(:,ISH2),MASK=GISC(:) )
!
!*       1.2.4  find the next lower indice for the ZLBDAI and for the ZLBDAI
!               in the geometrical set of (Lbda_I,Lbda_I) couplet use to
!               tabulate the SACCS-kernel
!
      ZVEC1(1:IGISC) = MAX( 1.0001, MIN( FLOAT(LIMAC%NSCLBDAI)-0.0001,           &
                                         LIMAC%XSCINTP1I * LOG( ZVEC1(1:IGISC) ) + LIMAC%XSCINTP2I ) )
      IVEC1(1:IGISC) = INT( ZVEC1(1:IGISC) )
      ZVEC1(1:IGISC) = ZVEC1(1:IGISC) - FLOAT( IVEC1(1:IGISC) )
      !
      ZVEC2(1:IGISC) = MAX( 1.0001, MIN( FLOAT(LIMAC%NSCLBDAI)-0.0001,           &
                                         LIMAC%XSCINTP1I * LOG( ZVEC2(1:IGISC) ) + LIMAC%XSCINTP2I ) )
      IVEC2(1:IGISC) = INT( ZVEC2(1:IGISC) )
      ZVEC2(1:IGISC) = ZVEC2(1:IGISC) - FLOAT( IVEC2(1:IGISC) )
!
!*       1.2.5  perform the bilinear interpolation of the normalized ISCI-kernel
!
      DO IJ = 1, IGISC
        ZVEC3(IJ) =  (  LIMAC%XKER_N_ISCI_SHAPE(ISH,ISH2,IVEC1(IJ)+1,IVEC2(IJ)+1) *  ZVEC2(IJ)          &
                      - LIMAC%XKER_N_ISCI_SHAPE(ISH,ISH2,IVEC1(IJ)+1,IVEC2(IJ)  ) * (ZVEC2(IJ) - 1.0) ) &
                                                                            *  ZVEC1(IJ) &
                   - (  LIMAC%XKER_N_ISCI_SHAPE(ISH,ISH2,IVEC1(IJ)  ,IVEC2(IJ)+1) *  ZVEC2(IJ)          &
                      - LIMAC%XKER_N_ISCI_SHAPE(ISH,ISH2,IVEC1(IJ)  ,IVEC2(IJ)  ) * (ZVEC2(IJ) - 1.0) ) &
                                                                            * (ZVEC1(IJ) - 1.0)
      END DO
      ZW1(:) = UNPACK( VECTOR=ZVEC3(:),MASK=GISC(:),FIELD=0.0 )
      !
      WHERE( GISC(:) )
         ZW2(:) = ZW1(:) * LIMAC%XFNISCI * &
              MIN(LIMAC%XCOLII * EXP( LIMAC%XCOLEXII*(PT(:)-CST%XTT) ),ZONE(:)) * &
              PCIT_SHAPE(:,ISH) * PCIT_SHAPE(:,ISH2) *                            &  
              PRHODREF(:)**(-LIMAP%XCEXVT-1.) *                                   &
              (LIMAC%XLBNISCI1 / (2. * PLBDI_SHAPE(:,ISH)**2) +                   &
               LIMAC%XLBNISCI2 / (PLBDI_SHAPE(:,ISH) * PLBDI_SHAPE(:,ISH2) ) +    &
               LIMAC%XLBNISCI1 / (2. * PLBDI_SHAPE(:,ISH2)**2)                    ) !++cb-- ISH --> ISH2
        !
        P_SHCI_ISC(:,ISH)  = P_SHCI_ISC(:,ISH)  - ZW2(:)
        P_SHCI_ISC(:,ISH2) = P_SHCI_ISC(:,ISH2) - ZW2(:) 
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
  ISH2 = 3
  DO ISH = 1, LIMAP%NNB_CRYSTAL_SHAPE
    ZW1(:) = 0.
    ZW2(:) = 0.
    !
    GISC(:) = PCIT_SHAPE(:,ISH)  > LIMAP%XCTMIN(4) .AND. PRIT_SHAPE(:,ISH)  > LIMAP%XRTMIN(4) .AND. & 
              PCIT_SHAPE(:,ISH2) > LIMAP%XCTMIN(4) .AND. PRIT_SHAPE(:,ISH2) > LIMAP%XRTMIN(4)
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
      ZVEC1(:) = PACK( PLBDI_SHAPE(:,ISH),MASK=GISC(:) )
      ZVEC2(:) = PACK( PLBDI_SHAPE(:,ISH2),MASK=GISC(:) )
!
!*       2.2    find the next lower indice for the ZLBDAI and for the ZLBDAI
!               in the geometrical set of (Lbda_I,Lbda_I) couplet use to
!               tabulate the SACCS-kernel
!
      ZVEC1(1:IGISC) = MAX( 1.0001, MIN( FLOAT(LIMAC%NSCLBDAI)-0.0001,           &
                                         LIMAC%XSCINTP1I * LOG( ZVEC1(1:IGISC) ) + LIMAC%XSCINTP2I ) )
      IVEC1(1:IGISC) = INT( ZVEC1(1:IGISC) )
      ZVEC1(1:IGISC) = ZVEC1(1:IGISC) - FLOAT( IVEC1(1:IGISC) )
      !
      ZVEC2(1:IGISC) = MAX( 1.0001, MIN( FLOAT(LIMAC%NSCLBDAI)-0.0001,           &
                                         LIMAC%XSCINTP1I * LOG( ZVEC2(1:IGISC) ) + LIMAC%XSCINTP2I ) )
      IVEC2(1:IGISC) = INT( ZVEC2(1:IGISC) )
      ZVEC2(1:IGISC) = ZVEC2(1:IGISC) - FLOAT( IVEC2(1:IGISC) )
!
!*       2.3    perform the bilinear interpolation of the normalized ISCI-kernel
!
      DO IJ = 1, IGISC
        ZVEC3(IJ) = (  LIMAC%XKER_N_ISCI_SHAPE(ISH,ISH2,IVEC1(IJ)+1,IVEC2(IJ)+1) *  ZVEC2(IJ)          &
                     - LIMAC%XKER_N_ISCI_SHAPE(ISH,ISH2,IVEC1(IJ)+1,IVEC2(IJ)  ) * (ZVEC2(IJ) - 1.0) ) &
                                                                           *  ZVEC1(IJ)          &
                  - (  LIMAC%XKER_N_ISCI_SHAPE(ISH,ISH2,IVEC1(IJ)  ,IVEC2(IJ)+1) *  ZVEC2(IJ)          &
                     - LIMAC%XKER_N_ISCI_SHAPE(ISH,ISH2,IVEC1(IJ)  ,IVEC2(IJ)  ) * (ZVEC2(IJ) - 1.0) ) &
                                                                           * (ZVEC1(IJ) - 1.0)
      END DO
      ZW1(:) = UNPACK( VECTOR=ZVEC3(:),MASK=GISC(:),FIELD=0.0 ) !! 
      !
      DO IJ = 1, IGISC
        ZVEC3(IJ) = (  LIMAC%XKER_R_ISCI_SHAPE(ISH,IVEC1(IJ)+1,IVEC2(IJ)+1) *  ZVEC2(IJ)          &
                     - LIMAC%XKER_R_ISCI_SHAPE(ISH,IVEC1(IJ)+1,IVEC2(IJ)  ) * (ZVEC2(IJ) - 1.0) ) &
                                                                      *  ZVEC1(IJ)          &
                  - (  LIMAC%XKER_R_ISCI_SHAPE(ISH,IVEC1(IJ)  ,IVEC2(IJ)+1) *  ZVEC2(IJ)          &
                     - LIMAC%XKER_R_ISCI_SHAPE(ISH,IVEC1(IJ)  ,IVEC2(IJ)  ) * (ZVEC2(IJ) - 1.0) ) &
                                                                      * (ZVEC1(IJ) - 1.0)
      END DO
      ZW3(:) = UNPACK( VECTOR=ZVEC3(:),MASK=GISC(:),FIELD=0.0 ) !! 
      !
      WHERE( GISC(:) )
         ZW2(:) = ZW1(:) * LIMAC%XFNISCI * &
              MIN( LIMAC%XCOLII*EXP( LIMAC%XCOLEXII*(PT(:)-CST%XTT) ),ZONE(:) ) * &
                          PCIT_SHAPE(:,ISH) * PCIT_SHAPE(:,ISH2) *                      &  
                          PRHODREF(:)**(-LIMAP%XCEXVT-1.) *                                   &
                         (LIMAC%XLBNISCI1 / (2. * PLBDI_SHAPE(:,ISH)**2) +                    &
                          LIMAC%XLBNISCI2 / (PLBDI_SHAPE(:,ISH) * PLBDI_SHAPE(:,ISH2))  +     &
                          LIMAC%XLBNISCI1 / (2. * PLBDI_SHAPE(:,ISH2)**2))                       !++cb-- ISH --> ISH2
         ZW4(:) = ZW3(:) * LIMAC%XFISCS(ISH) * &
              MIN( LIMAC%XCOLII*EXP( LIMAC%XCOLEXII*(PT(:)-CST%XTT) ),ZONE(:)) * &
                          PCIT_SHAPE(:,ISH) * PCIT_SHAPE(:,ISH2) *                         &
                          PRHODREF(:)**(-LIMAP%XCEXVT-1.) *                                      &
!++cb++ 15/04/23 Y=BUR (ISH2 = 3), X=PLA/COL/DRO/BUR (ISH)
!                         (XLBISCS1(ISH) / (PLBDI_SHAPE(:,ISH)**2) +                        &
!                          XLBISCS2(ISH) / (PLBDI_SHAPE(:,ISH) * PLBDI_SHAPE(:,ISH2))  +    &
!                          XLBISCS3(ISH) / (PLBDI_SHAPE(:,ISH2)**2))
                         (LIMAC%XLBISCS1(ISH2) / (PLBDI_SHAPE(:,ISH)**2 * PLBDI_SHAPE(:,ISH2)**LIMAC%XBI_SHAPE(ISH2)) + &
                          LIMAC%XLBISCS2(ISH2) / (PLBDI_SHAPE(:,ISH)    * PLBDI_SHAPE(:,ISH2)**(1.+LIMAC%XBI_SHAPE(ISH2))) + &
                          LIMAC%XLBISCS3(ISH2) / (                        PLBDI_SHAPE(:,ISH2)**(2.+LIMAC%XBI_SHAPE(ISH2))))
!--cb--
        !
        P_SHCI_ISCS(:,ISH)  = P_SHCI_ISCS(:,ISH)  - ZW2(:)
        P_SHCI_ISCS(:,ISH2) = P_SHCI_ISCS(:,ISH2) - ZW2(:) 
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
  GISC(:) = PCIT(:) > LIMAP%XCTMIN(4) .AND. PRIT(:) > LIMAP%XRTMIN(4)
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
    ZVEC1(1:IGISC) = MAX( 1.0001, MIN( FLOAT(LIMAC%NSCLBDAI)-0.0001,           &
                                       LIMAC%XSCINTP1I * LOG( ZVEC1(1:IGISC) ) + LIMAC%XSCINTP2I ) )
    IVEC1(1:IGISC) = INT( ZVEC1(1:IGISC) )
    ZVEC1(1:IGISC) = ZVEC1(1:IGISC) - FLOAT( IVEC1(1:IGISC) )
!
!*       3.3    perform the bilinear interpolation of the normalized ISCI-kernel
!
    DO IJ = 1, IGISC
      ZVEC3(IJ) = (  LIMAC%XKER_N_ISCI(IVEC1(IJ)+1,IVEC1(IJ)+1) *  ZVEC1(IJ)          &
                   - LIMAC%XKER_N_ISCI(IVEC1(IJ)+1,IVEC1(IJ)  ) * (ZVEC1(IJ) - 1.0) ) &
                                                          * ZVEC1(IJ)           &
                - (  LIMAC%XKER_N_ISCI(IVEC1(IJ)  ,IVEC1(IJ)+1) *  ZVEC1(IJ)          &
                   - LIMAC%XKER_N_ISCI(IVEC1(IJ)  ,IVEC1(IJ)  ) * (ZVEC1(IJ) - 1.0) ) &
                                                          * (ZVEC1(IJ) - 1.0)
    END DO
    ZW1(:) = UNPACK( VECTOR=ZVEC3(:),MASK=GISC(:),FIELD=0.0 ) !! 
    !
    WHERE( GISC(:) )
       P_CI_ISC(:) = - ZW1(:) * LIMAC%XFNISCI * &
            MIN(LIMAC%XCOLII * EXP( LIMAC%XCOLEXII*(PT(:)-CST%XTT) ),ZONE(:)) * PCIT(:)**2 &  
                             * PRHODREF(:)**(-LIMAP%XCEXVT-1.) * (LIMAC%XLBNISCI1 + LIMAC%XLBNISCI2) / PLBDI(:)**2
    END WHERE
    DEALLOCATE(IVEC1)
    DEALLOCATE(ZVEC1)
    DEALLOCATE(ZVEC3)
  END IF
END IF
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('LIMA_ICE_SELF_COLLECTION', 1, ZHOOK_HANDLE)
END SUBROUTINE LIMA_ICE_SELF_COLLECTION
END MODULE MODE_LIMA_ICE_SELF_COLLECTION
