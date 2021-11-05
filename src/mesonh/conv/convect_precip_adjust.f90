!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 conv 2006/05/18 13:07:25
!-----------------------------------------------------------------
!     #################
      MODULE MODI_CONVECT_PRECIP_ADJUST
!     #################
!
INTERFACE
!
       SUBROUTINE CONVECT_PRECIP_ADJUST( KLON, KLEV,                        &
                                        PPRES, PUMF, PUER, PUDR,           &
                                        PUPR, PUTPR, PURW,                 &
                                        PDMF, PDER, PDDR, PDTHL, PDRW,     &
                                        PPREF, PTPR, PMIXF, PDTEVR,        &
                                        KLFS, KDBL, KLCL, KCTL, KETL,      &
                                        PDTEVRF )

!
INTEGER,                    INTENT(IN) :: KLON  ! horizontal dimension
INTEGER,                    INTENT(IN) :: KLEV  ! vertical dimension
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PPRES ! pressure (Pa)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PURW  ! updraft total water (kg/kg) 
REAL, DIMENSION(KLON),      INTENT(IN) :: PUTPR ! updraft  total precipit. (kg/s
REAL, DIMENSION(KLON),      INTENT(IN) :: PPREF ! precipitation efficiency
REAL, DIMENSION(KLON),      INTENT(IN) :: PMIXF ! critical mixed fraction at LCL
INTEGER, DIMENSION(KLON),   INTENT(IN) :: KLCL  ! contains vert. index of LCL
INTEGER, DIMENSION(KLON),   INTENT(IN) :: KCTL  ! contains vert. index of CTL
INTEGER, DIMENSION(KLON),   INTENT(IN) :: KETL  ! contains vert. index of equilibrium 
						! (zero buoyancy) level 
INTEGER, DIMENSION(KLON),  INTENT(INOUT) :: KLFS ! contains vert. index of LFS
INTEGER, DIMENSION(KLON),  INTENT(INOUT) :: KDBL ! contains vert. index of DBL
!
REAL, DIMENSION(KLON),      INTENT(INOUT) :: PDTEVR ! total downdraft evaporation
                                                    ! rate at LFS   
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PDTEVRF! downdraft evaporation rate
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PUMF   ! updraft mass flux (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PUER   ! updraft entrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PUDR   ! updraft detrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PUPR   ! updraft  precipit. (kg/s)     
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PDMF   ! downdraft mass flux (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PDER   ! downdraft entrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PDDR   ! downdraft detrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PDTHL  ! downdraft enthalpy (J/kg)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PDRW   ! downdraft total water (kg/kg)
!
REAL, DIMENSION(KLON),     INTENT(OUT)   :: PTPR    ! total precipitation (kg/s) 
                                                 ! = downdraft precipitation
!
END SUBROUTINE CONVECT_PRECIP_ADJUST
!
END INTERFACE
!
END MODULE MODI_CONVECT_PRECIP_ADJUST
!     ######################################################################
      SUBROUTINE CONVECT_PRECIP_ADJUST( KLON, KLEV,                        &
                                        PPRES, PUMF, PUER, PUDR,           &
                                        PUPR, PUTPR, PURW,                 &
                                        PDMF, PDER, PDDR, PDTHL, PDRW,     &
                                        PPREF, PTPR, PMIXF, PDTEVR,        &
                                        KLFS, KDBL, KLCL, KCTL, KETL,      &
                                        PDTEVRF )
!     ######################################################################
!
!!**** Adjust up- and downdraft mass fluxes to be consistent with the
!!     mass transport at the LFS given by the precipitation efficiency
!!     relation. 
!!
!!
!!    PURPOSE                                                       
!!    -------
!!      The purpose of this routine is to adjust up- and downdraft mass
!!      fluxes below the LFS to be consistent with the precipitation
!!      efficiency relation
!!
!!
!!
!!**  METHOD
!!    ------
!!      
!!
!!    EXTERNAL
!!    --------
!!     None
!!     
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!     Module MODD_CONVPAREXT
!!          JCVEXB, JCVEXT     ! extra levels on the vertical boundaries
!!
!!     Module MODD_CONVPAR
!!        XUSRDPTH             ! pressure depth to compute updraft humidity
!!                             ! supply rate for downdraft
!!
!!    REFERENCE
!!    ---------
!!
!!      Book1,2 of documentation ( routine CONVECT_PRECIP_ADJUST)
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    07/11/95 
!!   Last modified  04/10/97
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CONVPAREXT
USE MODD_CONVPAR
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
INTEGER,                    INTENT(IN) :: KLON  ! horizontal dimension
INTEGER,                    INTENT(IN) :: KLEV  ! vertical dimension
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PPRES ! pressure (Pa)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PURW  ! updraft total water (kg/kg) 
REAL, DIMENSION(KLON),      INTENT(IN) :: PUTPR ! updraft  total precipit. (kg/s
REAL, DIMENSION(KLON),      INTENT(IN) :: PPREF ! precipitation efficiency
REAL, DIMENSION(KLON),      INTENT(IN) :: PMIXF ! critical mixed fraction at LCL
INTEGER, DIMENSION(KLON),   INTENT(IN) :: KLCL  ! contains vert. index of LCL
INTEGER, DIMENSION(KLON),   INTENT(IN) :: KCTL  ! contains vert. index of CTL
INTEGER, DIMENSION(KLON),   INTENT(IN) :: KETL  ! contains vert. index of equilibrium 
						! (zero buoyancy) level 
INTEGER, DIMENSION(KLON),  INTENT(INOUT) :: KLFS ! contains vert. index of LFS
INTEGER, DIMENSION(KLON),  INTENT(INOUT) :: KDBL ! contains vert. index of DBL
!
REAL, DIMENSION(KLON),      INTENT(INOUT) :: PDTEVR ! total downdraft evaporation
                                                    ! rate at LFS   
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PDTEVRF! downdraft evaporation rate
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PUMF   ! updraft mass flux (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PUER   ! updraft entrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PUDR   ! updraft detrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PUPR   ! updraft  precipit. (kg/s)     
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PDMF   ! downdraft mass flux (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PDER   ! downdraft entrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PDDR   ! downdraft detrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PDTHL  ! downdraft enthalpy (J/kg)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PDRW   ! downdraft total water (kg/kg)
!
REAL, DIMENSION(KLON),     INTENT(OUT)   :: PTPR    ! total precipitation (kg/s) 
                                                 ! = downdraft precipitation
!
!*       0.2   Declarations of local variables :
!
INTEGER :: IIE, IKB, IKE        ! horizontal + vertical loop bounds
INTEGER :: JK, JKT1, JKT2, JKT3 ! vertical loop index
INTEGER :: JI                   ! horizontal loop index
!
INTEGER, DIMENSION(KLON) :: IPRL
REAL, DIMENSION(KLON)    :: ZWORK1, ZWORK2, ZWORK3,     &
				    ZWORK4, ZWORK5, ZWORK6 ! work arrays
!
!
!-------------------------------------------------------------------------------
!
!        0.3   Set loop bounds
!              ---------------
!
IKB  = 1 + JCVEXB 
IKE  = KLEV - JCVEXT 
IIE  = KLON
JKT1 = MAXVAL( KLFS(:) )
JKT2 = MAXVAL( KCTL(:) )
JKT3 = MINVAL( KLCL(:) )
!
!
!        1.    Set some output variables for columns where no downdraft 
!              exists. Exit if there is no downdraft at all.
!              --------------------------------------------------------
!
IPRL(:) = IKB
PTPR(:) = 0.
!
WHERE ( PDTEVR(:) == 0. )
     PTPR(:)    = PUTPR(:)  ! no downdraft evaporation => no downdraft, all
			    ! precipitation occurs in updraft
END WHERE
IF ( COUNT( PDTEVR(:) > 0. ) == 0 )  THEN  ! exit routine if no downdraft exists
  RETURN
ENDIF
!
!*       2.     The total mass transported from the updraft to the down-  
!               draft at the LFS must be consistent with the three water
!               budget terms :
!               ---------------------------------------------------------
!
!*       2.1    Downdraft evaporation rate at the DBL. The evaporation
!               rate in downdraft must be consistent with precipitation
!               efficiency relation.
!               --------------------------------------------------------
!
!
DO JI = 1, IIE
     JK = KLFS(JI)
     ZWORK1(JI) = PDTEVR(JI) / MIN( -1.E-1, PDMF(JI,JK) )
     ZWORK6(JI) = PDMF(JI,JK)
END DO
!
!*       2.2    Some preliminar computations for downdraft = total 
!               precipitation rate. The precipitation is evaluated in 
!               a layer thickness DP=XUSRDPTH=165 hPa above the LCL.
!               The difference between updraft precipitation and downdraft
!               precipitation (updraft supply rate) is used to drive the
!               downdraft through evaporational cooling.
!               --------------------------------------------------------
!
DO JI = 1, IIE
     JK = KLCL(JI)
     ZWORK5(JI) = PPRES(JI,JK)
END DO
!
PTPR(:) = 0.
DO JK = JKT3, JKT2
    WHERE ( JK >= KLCL(:) .AND. PPRES(:,JK) >= ZWORK5(:) - XUSRDPTH )
	PTPR(:) = PTPR(:) + PUPR(:,JK)
	IPRL(:) = JK
    END WHERE
END DO
IPRL(:) = MIN( KETL(:), IPRL(:) )
!
DO JI = 1, IIE
     JK = IPRL(JI)
     PTPR(JI) = PUMF(JI,JK+1) * PURW(JI,JK+1) + PTPR(JI) 
END DO
!
PTPR(:) = PPREF(:) * MIN( PUTPR(:), PTPR(:) )
ZWORK4(:) = PUTPR(:) - PTPR(:) 
!
!
!*       2.3    Total amount of precipitation that falls out of the up-
!               draft between the LCL and the LFS.
!               Condensate transfer from up to downdraft at LFS
!               ---------------------------------------------------------
!
ZWORK5(:) = 0.
DO JK = JKT3, JKT1
     WHERE ( JK >= KLCL(:) .AND. JK <= KLFS(:) )
	   ZWORK5(:) = ZWORK5(:) +  PUPR(:,JK)
     END WHERE
END DO
!
DO JI = 1, IIE
     JK = KLFS(JI)
     ZWORK2(JI) = ( 1. - PPREF(JI) ) * ZWORK5(JI) *                     &
                  ( 1. - PMIXF(JI) ) / MAX( 1.E-1, PUMF(JI,JK) )
END DO
!
!
!*       2.4    Increase the first guess downdraft mass flux to satisfy
!               precipitation efficiency relation.
!               If downdraft does not evaporate any water at the DBL for  
!               the specified relative humidity, or if the corrected mass 
!               flux at the LFS is positive no downdraft is allowed
!               ---------------------------------------------------------
!    
!
!ZWORK1(:) = ZWORK4(:) / ( ZWORK1(:) + ZWORK2(:) + 1.E-8 ) 
ZWORK1(:) = -ZWORK4(:) / ( -ZWORK1(:) + ZWORK2(:) + 1.E-8 )
ZWORK2(:) = ZWORK1(:) / MIN( -1.E-1, ZWORK6(:) ) ! ratio of budget consistent to actual DMF
!
ZWORK3(:) = 1.
ZWORK6(:) = 1.
WHERE ( ZWORK1(:) > 0. .OR. PDTEVR(:) < 1. ) 
   KDBL(:)   = IKB
   KLFS(:)   = IKB
   PDTEVR(:) = 0. 
   ZWORK2(:) = 0.
   ZWORK3(:) = 0.
   ZWORK6(:) = 0.
END WHERE
!
DO JK = IKB, JKT1   
     PDMF(:,JK)  = PDMF(:,JK)  * ZWORK2(:)
     PDER(:,JK)  = PDER(:,JK)  * ZWORK2(:)  
     PDDR(:,JK)  = PDDR(:,JK)  * ZWORK2(:)  
   PDTEVRF(:,JK) = PDTEVRF(:,JK)* ZWORK2(:)  
     PDRW(:,JK)  = PDRW(:,JK)  * ZWORK3(:)  
     PDTHL(:,JK) = PDTHL(:,JK) * ZWORK3(:)  
END DO     
ZWORK4(:) = ZWORK2(:)
!
!
!*       3.     Increase updraft mass flux, mass detrainment rate, and water  
!               substance detrainment rates to be consistent with the transfer
!               of the estimated mass from the up- to the downdraft at the LFS
!               --------------------------------------------------------------
!
DO JI = 1, IIE
    JK = KLFS(JI)
    ZWORK2(JI) = ( 1. - ZWORK6(JI) ) + ZWORK6(JI) *                   &
		  ( PUMF(JI,JK) - ( 1. - PMIXF(JI) ) * ZWORK1(JI) ) / &
		  MAX( 1.E-1, PUMF(JI,JK) )
END DO
!
!
JKT1  = MAXVAL( KLFS(:) )  ! value of KLFS might have been reset to IKB above
DO JK = IKB, JKT1
    DO JI = 1, IIE
      IF ( JK <= KLFS(JI) ) THEN
	PUMF(JI,JK)  = PUMF(JI,JK)  * ZWORK2(JI) 
	PUER(JI,JK)  = PUER(JI,JK)  * ZWORK2(JI)
	PUDR(JI,JK)  = PUDR(JI,JK)  * ZWORK2(JI)
	PUPR(JI,JK)  = PUPR(JI,JK)  * ZWORK2(JI)
      END IF
    END DO
END DO
!
!
!*       4.     Increase total = downdraft precipitation and evaporation rate
!               -------------------------------------------------------------
!
WHERE ( PDTEVR(:) > 0. )
    PTPR(:)    = PTPR(:) + PPREF(:) * ZWORK5(:) * ( ZWORK2(:) - 1. )
    PDTEVR(:)  = PUTPR(:) - PTPR(:)
    PDTEVRF(:,IKB+1)  = PDTEVR(:)
ELSEWHERE
    PTPR(:)    = PUTPR(:)
END WHERE
!
!
END SUBROUTINE CONVECT_PRECIP_ADJUST
