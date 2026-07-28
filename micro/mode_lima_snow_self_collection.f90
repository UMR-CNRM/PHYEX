!MNH_LIC Copyright 2018-2024 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-------------------------------------------------------------------------------
MODULE MODE_LIMA_SNOW_SELF_COLLECTION
  IMPLICIT NONE
CONTAINS
!     #############################################################
  SUBROUTINE LIMA_SNOW_SELF_COLLECTION (CST, LIMAP, LIMAC, KSIZE, ODCOMPUTE,   &
                                        PRHODREF, PT,       &
                                        PRST, PCST, PLBDS,  &
                                        P_CS_SSC            )
!     #############################################################
!
!!    PURPOSE
!!    -------
!!      Compute the self-collection and physical break-up of snow
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
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRST    ! Snow mr at t
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PCST    ! Snow C. at t
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLBDS   ! 
!
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_CS_SSC 
!
!*       0.2   Declarations of local variables :
!
REAL, DIMENSION(SIZE(PCST)) :: &
                                           ZW1, & ! work arrays
                                           ZW2
LOGICAL, DIMENSION(SIZE(PCST)) :: GSSC
INTEGER :: IGSSC, IL
INTEGER, DIMENSION(:), ALLOCATABLE :: IVEC1        ! Vectors of indices
REAL,    DIMENSION(:), ALLOCATABLE :: ZVEC1, ZVEC3 ! Work vectors
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!
!*       1.     Snow self-collection and break-up
!               ---------------------------------
!
!
IF (LHOOK) CALL DR_HOOK('LIMA_SNOW_SELF_COLLECTION', 0, ZHOOK_HANDLE)
P_CS_SSC(:)=0.
!
ZW1(:) =0.
ZW2(:) =0.
!
GSSC(:) = PCST(:)>LIMAP%XCTMIN(5) .AND. PRST(:)>LIMAP%XRTMIN(5)
IGSSC = COUNT(GSSC(:))
!
IF( IGSSC>0 ) THEN
!
!        1.3N.0  allocations
!
   ALLOCATE(ZVEC1(IGSSC))
   ALLOCATE(IVEC1(IGSSC))
!
!        1.3N.1  select the (ZLBDAS,ZLBDAS) couplet
!
   ZVEC1(:) = PACK( PLBDS(:),MASK=GSSC(:) )
!
!        1.3N.2  find the next lower indice for the ZLBDAS and for the ZLBDAS
!               in the geometrical set of (Lbda_s,Lbda_s) couplet use to
!               tabulate the SACCS-kernel
!
   ZVEC1(1:IGSSC) = MAX( 1.0001, MIN( FLOAT(LIMAC%NSCLBDAS)-0.0001,           &
        LIMAC%XSCINTP1S * LOG( ZVEC1(1:IGSSC) ) + LIMAC%XSCINTP2S ) )
   IVEC1(1:IGSSC) = INT( ZVEC1(1:IGSSC) )
   ZVEC1(1:IGSSC) = ZVEC1(1:IGSSC) - FLOAT( IVEC1(1:IGSSC) )
!
!        1.3N.3 perform the bilinear interpolation of the normalized
!               SSCS-kernel
!
   ALLOCATE(ZVEC3(IGSSC))
   DO IL = 1,IGSSC
      ZVEC3(IL) =  (   LIMAC%XKER_N_SSCS(IVEC1(IL)+1,IVEC1(IL)+1)* ZVEC1(IL)          &
                    -  LIMAC%XKER_N_SSCS(IVEC1(IL)+1,IVEC1(IL)  )*(ZVEC1(IL) - 1.0) ) &
                                                                         * ZVEC1(IL) &
                 - (   LIMAC%XKER_N_SSCS(IVEC1(IL)  ,IVEC1(IL)+1)* ZVEC1(IL)          &
                    -  LIMAC%XKER_N_SSCS(IVEC1(IL)  ,IVEC1(IL)  )*(ZVEC1(IL) - 1.0) ) &
                                                           * (ZVEC1(IL) - 1.0)
   END DO
   ZW1(:) = UNPACK( VECTOR=ZVEC3(:),MASK=GSSC(:),FIELD=0.0 ) !! NSACCS
   DEALLOCATE(ZVEC3)
!
   WHERE( GSSC(:) )
      P_CS_SSC(:) = - LIMAC%XFNSSCS * ZW1(:) * EXP( LIMAC%XCOLEXSS*(PT(:)-CST%XTT) ) * PCST(:)**2 &  
                    * PRHODREF(:)**(-LIMAP%XCEXVT-1.) * (LIMAC%XLBNSSCS1+LIMAC%XLBNSSCS2) / PLBDS(:)**2
   END WHERE
   DEALLOCATE(IVEC1)
   DEALLOCATE(ZVEC1)
END IF
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('LIMA_SNOW_SELF_COLLECTION', 1, ZHOOK_HANDLE)
END SUBROUTINE LIMA_SNOW_SELF_COLLECTION
END MODULE MODE_LIMA_SNOW_SELF_COLLECTION
