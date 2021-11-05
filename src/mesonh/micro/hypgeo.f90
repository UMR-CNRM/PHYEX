!MNH_LIC Copyright 1996-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!####################
MODULE MODI_HYPGEO
!####################
!
INTERFACE
!
FUNCTION HYPGEO(PA,PB,PC,PF,PX)  RESULT(PHYPGEO)
REAL, INTENT(IN)                                  :: PA,PB,PC,PF
REAL, INTENT(IN)                                  :: PX
REAL                                              :: PHYPGEO
END FUNCTION HYPGEO
!
END INTERFACE
!
END MODULE MODI_HYPGEO
!     #############################################
      FUNCTION HYPGEO(PA,PB,PC,PF,PX)  RESULT(PHYPGEO)
!     #############################################
!
!
!!****  *HYPGEO* -  hypergeometric  function
!!
!!
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the hypergeometric
!!   function of its argument.
!!
!!                             
!!                          A*B        (A+1)A*(B+1)B   X^2           
!!    HYPGEO(A,B,C,X)= 1 + ----- * X + ------------- * --- + ... +
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
USE MODI_GAMMA
USE MODI_HYPSER
!
IMPLICIT NONE
!
!*       0.1 declarations of arguments and result
!
REAL, INTENT(IN)                     :: PA,PB,PC,PF
REAL, INTENT(IN)                     :: PX
REAL                                 :: PHYPGEO
!
!*       0.2 declarations of local variables
!
!
INTEGER                              :: JN
INTEGER                              :: ITMAX=100
REAL                                 :: ZEPS,ZTEMP
REAL                                 :: ZFPMIN=1.E-30
REAL                                 :: ZXH
REAL                                 :: ZX0,ZX1,ZZA,ZZB,ZZC,ZZD,Y(2)
!
!------------------------------------------------------------------------------
!
!
ZEPS = 4.E-2
ZXH = PF * PX**2.0
IF (ZXH.LT.(1-ZEPS)) THEN
  CALL HYPSER(PA,PB,PC,-ZXH,PHYPGEO)
ELSE IF (ZXH.GT.(1.+ZEPS)) THEN
  ZXH = 1./ZXH
  CALL HYPSER(PA,PA-PC+1.,PA-PB+1.,-ZXH,PHYPGEO)
  PHYPGEO = PHYPGEO*ZXH**(PA)*                         &
           (GAMMA(PC)*GAMMA(PB-PA)/(GAMMA(PB)*GAMMA(PC-PA)))
  CALL HYPSER(PB,PB-PC+1.,PB-PA+1.,-ZXH,ZTEMP)
  PHYPGEO = PHYPGEO+ZTEMP*ZXH**(PB)*                         &
           (GAMMA(PC)*GAMMA(PA-PB)/(GAMMA(PA)*GAMMA(PC-PB)))
ELSE
  ZX0 = (1.-ZEPS)
  ZX1 = 1./(1.+ZEPS)
  CALL HYPSER(PA,PA-PC+1.,PA-PB+1.,-ZX1,PHYPGEO)
  PHYPGEO = PHYPGEO*ZX1**(PA)*                         &
           (GAMMA(PC)*GAMMA(PB-PA)/(GAMMA(PB)*GAMMA(PC-PA)))
  CALL HYPSER(PB,PB-PC+1.,PB-PA+1.,-ZX1,ZTEMP)
  PHYPGEO = PHYPGEO+ZTEMP*ZX1**(PB)*                         &
           (GAMMA(PC)*GAMMA(PA-PB)/(GAMMA(PA)*GAMMA(PC-PB)))
  CALL HYPSER(PA,PB,PC,-ZX0,ZTEMP)
  PHYPGEO = ZTEMP + (ZXH-ZX0)*(PHYPGEO-ZTEMP)/(2.*ZEPS)
ENDIF
END
