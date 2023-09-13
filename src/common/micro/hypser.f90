!MNH_LIC Copyright 1996-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!####################
MODULE MODI_HYPSER
!####################
!
IMPLICIT NONE
INTERFACE
!
SUBROUTINE HYPSER(PA,PB,PC,PX,PHYP)  
IMPLICIT NONE
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
!  P. Wautelet 10/04/2019: replace ABORT and STOP calls by Print_msg
!  P. Wautelet 26/04/2019: replace non-standard FLOAT function by REAL function
!
!------------------------------------------------------------------------------
!
!*       0. DECLARATIONS
!
!
use mode_msg
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
         ZFAC = ZFAC * ZXH / REAL(JN)
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
IF (JFLAG == 0) call Print_msg(NVERB_FATAL,'GEN','HYPSER','convergence failure')
!
END SUBROUTINE HYPSER
