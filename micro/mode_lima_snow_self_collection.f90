!MNH_LIC Copyright 2018-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-------------------------------------------------------------------------------
MODULE MODE_LIMA_SNOW_SELF_COLLECTION
  IMPLICIT NONE
CONTAINS
!     #############################################################
  SUBROUTINE LIMA_SNOW_SELF_COLLECTION (LDCOMPUTE,          &
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
USE MODD_CST,             ONLY : XTT
USE MODD_PARAM_LIMA,      ONLY : XRTMIN, XCTMIN, XCEXVT
USE MODD_PARAM_LIMA_COLD, ONLY : NSCLBDAS, XSCINTP1S, XSCINTP2S, XKER_N_SSCS, XFNSSCS, XCOLEXSS, &
                                 XLBNSSCS1, XLBNSSCS2
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
REAL, DIMENSION(:),   INTENT(IN)    :: PRST    ! Snow mr at t
REAL, DIMENSION(:),   INTENT(IN)    :: PCST    ! Snow C. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDS   ! 
!
REAL, DIMENSION(:),   INTENT(OUT)   :: P_CS_SSC 
!
!*       0.2   Declarations of local variables :
!
REAL, DIMENSION(SIZE(PCST)) :: &
                                           ZW1, & ! work arrays
                                           ZW2
LOGICAL, DIMENSION(SIZE(PCST)) :: GSSC
INTEGER :: IGSSC, JJ
INTEGER, DIMENSION(:), ALLOCATABLE :: IVEC1        ! Vectors of indices
REAL,    DIMENSION(:), ALLOCATABLE :: ZVEC1, ZVEC3 ! Work vectors
!
!-------------------------------------------------------------------------------
!
!
!*       1.     Snow self-collection and break-up
!	        ---------------------------------
!
!
P_CS_SSC(:)=0.
!
ZW1(:) =0.
ZW2(:) =0.
!
GSSC(:) = PCST(:)>XCTMIN(5) .AND. PRST(:)>XRTMIN(5)
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
   ZVEC1(1:IGSSC) = MAX( 1.0001, MIN( FLOAT(NSCLBDAS)-0.0001,           &
        XSCINTP1S * LOG( ZVEC1(1:IGSSC) ) + XSCINTP2S ) )
   IVEC1(1:IGSSC) = INT( ZVEC1(1:IGSSC) )
   ZVEC1(1:IGSSC) = ZVEC1(1:IGSSC) - FLOAT( IVEC1(1:IGSSC) )
!
!        1.3N.3 perform the bilinear interpolation of the normalized
!               SSCS-kernel
!
   ALLOCATE(ZVEC3(IGSSC))
   DO JJ = 1,IGSSC
      ZVEC3(JJ) =  (   XKER_N_SSCS(IVEC1(JJ)+1,IVEC1(JJ)+1)* ZVEC1(JJ)          &
                    -  XKER_N_SSCS(IVEC1(JJ)+1,IVEC1(JJ)  )*(ZVEC1(JJ) - 1.0) ) &
                                                                         * ZVEC1(JJ) &
                 - (   XKER_N_SSCS(IVEC1(JJ)  ,IVEC1(JJ)+1)* ZVEC1(JJ)          &
                    -  XKER_N_SSCS(IVEC1(JJ)  ,IVEC1(JJ)  )*(ZVEC1(JJ) - 1.0) ) &
                                                           * (ZVEC1(JJ) - 1.0)
   END DO
   ZW1(:) = UNPACK( VECTOR=ZVEC3(:),MASK=GSSC(:),FIELD=0.0 ) !! NSACCS
   DEALLOCATE(ZVEC3)
!
   WHERE( GSSC(:) )
      P_CS_SSC(:) = - XFNSSCS * ZW1(:) * EXP( XCOLEXSS*(PT(:)-XTT) ) * PCST(:)**2 &  
                    * PRHODREF(:)**(-XCEXVT-1.) * (XLBNSSCS1+XLBNSSCS2) / PLBDS(:)**2
   END WHERE
   DEALLOCATE(IVEC1)
   DEALLOCATE(ZVEC1)
END IF
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_SNOW_SELF_COLLECTION
END MODULE MODE_LIMA_SNOW_SELF_COLLECTION
