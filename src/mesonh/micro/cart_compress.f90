!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! masdev4_7 BUG2 2007/06/29 16:52:14
!-----------------------------------------------------------------
!#########################
 MODULE MODI_CART_COMPRESS
!#########################
!
INTERFACE
!
FUNCTION CART_COMPRESS(PVARS) RESULT(PCOMPRESS)
!
USE MODD_BUDGET
!
REAL, DIMENSION(:,:,:), INTENT(IN)       :: PVARS     ! Source
REAL, DIMENSION(NBUIMAX,NBUJMAX,NBUKMAX) :: PCOMPRESS ! result
!
END FUNCTION CART_COMPRESS
!
END INTERFACE
!
END MODULE MODI_CART_COMPRESS
!     ###############################################
      FUNCTION CART_COMPRESS(PVARS) RESULT(PCOMPRESS) 
!     ###############################################
!
!!****  *CART_COMPRESS* - function to compress the Source in CART case. 
!!                           
!!
!!    PURPOSE
!!    -------
!       This function compresses or not the Source XVARS of the VARiable
!     VAR whose budget is analysed. This compression is controlled by 3 
!     logical switches for the budget in I,J and K directions (LBU_ICP,
!     LBU_JCP, LBU_KCP), in the budget box described by the lowest and
!     highest values of the I,J and K indices.  
!
!!**  METHOD
!!    ------
!!      The source PVARS is first transfered in a local array whose 
!!    dimensions correspond to the budget box. Then compressions
!!    are or aren't achieved depending on the logical switches.
!!
!!    EXTERNAL
!!    --------
!!       NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!       Module MODD_BUDGET
!!           LBU_ICP   : switch for compression in I direction
!!           LBU_JCP   : switch for compression in J direction
!!           LBU_KCP   : switch for compression in K direction
!!           NBUIL     : lowest I indice value of the budget box
!!           NBUJL     : lowest J indice value of the budget box
!!           NBUKL     : lowest K indice value of the budget box
!!           NBUIH     : highest I indice value of the budget box
!!           NBUJH     : highest J indice value of the budget box
!!           NBUKH     : highest K indice value of the budget box
!!           NBUIMAX   : dimension along I of the budget tabular
!!           NBUJMAX   : dimension along J of the budget tabular
!!           NBUKMAX   : dimension along K of the budget tabular
!!          
!!
!!
!!    REFERENCE
!!    ---------
!!      Book2 of MESO-NH documentation (function CART_COMPRESS)
!!
!!
!!    AUTHOR
!!    ------
!!	J. Nicolau       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original             27/02/95
!!      JP Pinty & J Escobar 12/10/98 Enable vectorization and remove 
!!                                     SUM functions
!!      V. Ducrocq           4/06/99  //
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_BUDGET
USE MODD_PARAMETERS , ONLY : JPVEXT
!
!
IMPLICIT NONE
!  
!  
!*       0.1   Declarations of arguments and result :
!
REAL, DIMENSION(:,:,:), INTENT(IN)       :: PVARS     ! Source 
REAL, DIMENSION(NBUIMAX,NBUJMAX,NBUKMAX) :: PCOMPRESS ! result
!
!*       0.2   Declarations of local variables :
! 
! 
REAL, DIMENSION (NBUSIH-NBUSIL+1,NBUSJH-NBUSJL+1,NBUKH-NBUKL+1) :: ZVARS ! 3D Work 
                                                                     ! array
REAL, DIMENSION (NBUSIH-NBUSIL+1,NBUKH-NBUKL+1) :: ZWORKIK ! 2D Work array
REAL, DIMENSION (NBUSJH-NBUSJL+1,NBUKH-NBUKL+1) :: ZWORKJK ! 2D Work array
REAL, DIMENSION (NBUSIH-NBUSIL+1,NBUSJH-NBUSJL+1) :: ZWORKIJ ! 2D Work array
! 
INTEGER                                          :: JJ,JK   ! loop indexes 
! 
!
!-------------------------------------------------------------------------------
!
!*	 1.     SOURCE TRANSFERT IN A LOCAL ARRAY 
!	        ---------------------------------
!JUAN
IF (SIZE (PCOMPRESS) .EQ. 0 ) RETURN
!JUAN
!
ZVARS(1:NBUSIH-NBUSIL+1,1:NBUSJH-NBUSJL+1,1:NBUKH-NBUKL+1) = &
            PVARS(NBUSIL:NBUSIH,NBUSJL:NBUSJH,NBUKL+JPVEXT:NBUKH+JPVEXT)
!
!-------------------------------------------------------------------------------
!
!*       2.     COMPRESSIONS IN I,J AND K DIRECTIONS
!	        ------------------------------------
!                                 
!
IF (LBU_ICP.AND.LBU_JCP.AND.LBU_KCP) THEN
  PCOMPRESS(1,1,1)=SUM(ZVARS)
!
ELSE IF (LBU_ICP.AND.LBU_JCP.AND..NOT.LBU_KCP) THEN
  ZWORKJK(:,:)    =SUM(ZVARS,1)
  PCOMPRESS(1,1,:)=SUM(ZWORKJK,1)
!
ELSE IF (LBU_ICP.AND..NOT.LBU_JCP.AND.LBU_KCP) THEN
  ZWORKIJ(:,:)=0.0
  DO JK = 1,NBUKH-NBUKL+1
    ZWORKIJ(:,:) = ZWORKIJ(:,:) + ZVARS(:,:,JK)
  END DO
  PCOMPRESS(1,:,1)=SUM(ZWORKIJ,1)
!
ELSE IF (.NOT.LBU_ICP.AND.LBU_JCP.AND.LBU_KCP) THEN
  ZWORKIK(:,:)=0.0
  DO JJ = 1,NBUSJH-NBUSJL+1
    ZWORKIK(:,:) = ZWORKIK(:,:) + ZVARS(:,JJ,:)
  END DO
  PCOMPRESS(:,1,1)=SUM(ZWORKIK,2)
!
ELSE IF (LBU_ICP.AND..NOT.LBU_JCP.AND..NOT.LBU_KCP) THEN 
  PCOMPRESS(1,:,:)=SUM(ZVARS,1)
!
ELSE IF (.NOT.LBU_ICP.AND.LBU_JCP.AND..NOT.LBU_KCP) THEN
  PCOMPRESS(:,1,:)=SUM(ZVARS,2)
!
ELSE IF (.NOT.LBU_ICP.AND..NOT.LBU_JCP.AND.LBU_KCP) THEN
  PCOMPRESS(:,:,1)=SUM(ZVARS,3)
!
ELSE  
  PCOMPRESS=ZVARS
!
END IF
!
!
END FUNCTION CART_COMPRESS                           
