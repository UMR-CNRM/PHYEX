!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source: /home/cvsroot/MESONH_RCS/CODEMNH/hypser.f90,v $ $Revision: 1.7 $
! MASDEV4_7 operators 2006/05/18 13:07:25
!-----------------------------------------------------------------
!####################
MODULE MODI_HYPSER
!####################
!
INTERFACE
!
SUBROUTINE HYPSER(PA,PB,PC,PX,PHYP)  
REAL, INTENT(IN)                                  :: PA,PB,PC
REAL, INTENT(IN)                                  :: PX
REAL, INTENT(INOUT)                               :: PHYP
END SUBROUTINE HYPSER
!
END INTERFACE
!
END MODULE MODI_HYPSER
!     #############################################
      SUBROUTINE HYPSER(PA,PB,PC,PX,PHYP)
!     #############################################
!
!
!!****  *HYPSER* -  hypergeometric  function
!!
!!
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the hypergeometric
!!   function of its argument.
!!
!!
!!                          A*B        (A+1)A*(B+1)B   X^2
!!    HYPSER(A,B,C,X)= 1 + ----- * X + ------------- * --- + ... +
!!                           C            (C+1)C        2
!!
!!                          (A+n)...A*(B+n)...B     X^n
!!                         --------------------- * ----- + ... ...
!!                               (C+n)...C           n!
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!    HYPSER
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!      Press, Teukolsky, Vetterling and Flannery: Numerical Recipes, 272
!!
!!
!!    AUTHOR
!!    ------
!!         Jean-Martial Cohard *LA/OMP*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original     31/12/96
!
!------------------------------------------------------------------------------
!
!*       0. DECLARATIONS
!
!
IMPLICIT NONE
!
!*       0.1 declarations of arguments and result
!
REAL, INTENT(IN)                                  :: PA,PB,PC
REAL, INTENT(IN)                                  :: PX
REAL, INTENT(INOUT)                               :: PHYP
!
!
!
!*       0.2 declarations of local variables
!
INTEGER                                           :: JN,JFLAG
REAL                                              :: ZXH,ZZA,ZZB,ZZC,ZFAC,ZTEMP
REAL                                              :: ZPREC
!
!------------------------------------------------------------------------------
!
ZPREC = 1.0E-04
ZXH = PX
ZFAC = 1.0
ZTEMP = ZFAC
ZZA = PA
ZZB = PB
ZZC = PC
JFLAG = 0
SERIE: DO JN = 1,5000
         ZFAC = ZFAC * ZZA * ZZB / ZZC
         ZFAC = ZFAC * ZXH / FLOAT(JN)
         PHYP = ZTEMP + ZFAC
         IF (ABS(PHYP-ZTEMP).LE.ZPREC) THEN
	   JFLAG = 1
	   EXIT SERIE
         END IF
  ZTEMP = PHYP
  ZZA = ZZA + 1.
  ZZB = ZZB + 1.
  ZZC = ZZC + 1.
END DO SERIE
IF (JFLAG == 0) THEN
  PRINT *,'CONVERGENCE FAILURE IN HYPSER'
!callabortstop
CALL ABORT
  STOP
END IF
!
END
