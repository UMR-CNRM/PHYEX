!MNH_LIC Copyright 2007-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-------------------------------------------------------------------------------
!      ####################
       MODULE MODI_LIMA_INIT_CCN_ACTIVATION_SPECTRUM
INTERFACE
   SUBROUTINE LIMA_INIT_CCN_ACTIVATION_SPECTRUM (CTYPE_CCN,XD,XSIGMA,XLIMIT_FACTOR,XK,XMU,XBETA,XKAPPA)
     !
     CHARACTER(LEN=*), INTENT(IN)  :: CTYPE_CCN         ! Aerosol type
     REAL,             INTENT(IN)  :: XD            ! Aerosol PSD modal diameter          
     REAL,             INTENT(IN)  :: XSIGMA        ! Aerosol PSD width
     REAL,             INTENT(OUT) :: XLIMIT_FACTOR ! C/Naer
     REAL,             INTENT(OUT) :: XK            ! k
     REAL,             INTENT(OUT) :: XMU           ! mu
     REAL,             INTENT(OUT) :: XBETA         ! beta
     REAL,             INTENT(OUT) :: XKAPPA        ! kappa
!
   END SUBROUTINE LIMA_INIT_CCN_ACTIVATION_SPECTRUM
END INTERFACE
END MODULE MODI_LIMA_INIT_CCN_ACTIVATION_SPECTRUM
!      ####################
!
!     #############################################################
      SUBROUTINE LIMA_INIT_CCN_ACTIVATION_SPECTRUM (CTYPE_CCN,XD,XSIGMA,XLIMIT_FACTOR,XK,XMU,XBETA,XKAPPA)
!     #############################################################

!!
!!
!!      PURPOSE
!!      -------
!!      
!!      Compute mu, k and beta parameters of the activation spectrum based on CCN
!!      characteristics (type and PSD)        
!!
!!
!!    AUTHOR
!!    ------
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!      S.    Berthet    * Laboratoire d'Aerologie*
!!      B.    Vié        * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original             ??/??/13 
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST, ONLY : XMV, XAVOGADRO, XBOLTZ, XRHOLW
!
USE MODI_GAMMA_INC
USE MODI_HYPGEO
USE MODI_HYPSER
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments : 
!
CHARACTER(LEN=*), INTENT(IN)  :: CTYPE_CCN          ! Aerosol type
REAL,             INTENT(IN)  :: XD             ! Aerosol PSD modal diameter          
REAL,             INTENT(IN)  :: XSIGMA         ! Aerosol PSD width
REAL,             INTENT(OUT) :: XLIMIT_FACTOR  ! C/Naer
REAL,             INTENT(OUT) :: XK             ! k
REAL,             INTENT(OUT) :: XMU            ! mu
REAL,             INTENT(OUT) :: XBETA          ! beta
REAL,             INTENT(OUT) :: XKAPPA         ! kappa
!
!*       0.2   Declarations of local variables :
!
INTEGER, PARAMETER            :: M = 1000        ! Number of points (S,Nccn) used to fit the spectra
INTEGER, PARAMETER            :: N = 3          ! Number of parameters to adjust
REAL, DIMENSION(N)            :: PARAMS         ! Parameters to adjust by the LM algorithm (k, mu, beta)
REAL, DIMENSION(M)            :: FVEC           ! Array to store the distance between theoretical and fitted spectra
INTEGER                       :: IFLAG          ! 
INTEGER                       :: INFO           ! 
REAL                          :: TOL = 1.E-16   ! Fit precision required
!
INTEGER                       :: II, IJ         ! Loop indices
!
REAL                          :: XW             ! 
REAL                          :: XDDRY = 0.1E-6 ! Dry diameter for which to compute Scrit
REAL                          :: XSCRIT         ! Scrit for dry diameter XDDRY
REAL                          :: XMIN  = 0.1E-6 ! minimum diameter for root search (m)
REAL                          :: XMAX  = 10.E-6 ! maximum diameter for root search (m)
REAL                          :: XPREC = 1.E-8  ! precision wanted for root (m)
!
!REAL                          :: XKAPPA         ! kappa coefficient
REAL, DIMENSION(M)            :: XS             ! saturation ratio (S=1.01 for a 1% supersaturation)
REAL, DIMENSION(M)            :: XDCRIT         ! critical diameters (m) for the chosen S values
REAL, DIMENSION(M)            :: XNCCN          ! fraction of the aerosols larger than XDCRIT (ie activable)
REAL, DIMENSION(1)            :: XT             ! temperature
!
!
!-------------------------------------------------------------------------------
!
!*       1.     Select kappa value based on CTYPE_CCN
!	        ---------------------------------
!
! Kappa values are from Petters and Kreidenweis (2007), table 1.
!
SELECT CASE (CTYPE_CCN)
CASE('NH42SO4','C') ! Ammonium sulfate
   XKAPPA = 0.61
CASE('NH4NO3')      ! Ammonium nitrate
   XKAPPA = 0.67
CASE('NaCl','M')    ! Sea Salt
   XKAPPA = 1.28
CASE('H2SO4')       ! Sulfuric acid
   XKAPPA = 0.90
CASE('NaNO3')       ! Sodium nitrate
   XKAPPA = 0.88
CASE('NaHSO4')      ! Sodium bisulfate
   XKAPPA = 0.91
CASE('Na2SO4')      ! Sodium sulfate
   XKAPPA = 0.80
CASE('NH43HSO42')   ! Letovicite (rare ammonium sulfate mineral)
   XKAPPA = 0.65
CASE('SOA')         ! Secondary organic aerosol (alpha-pinene, beta-pinene)
   XKAPPA = 0.1
CASE DEFAULT
   XKAPPA = 1.
END SELECT
!
!XT = (/ 270., 271., 272., 273., 274., 275., 276., 277., 278., 279., 280., 281., 282., 283., 284., 285., 286., 287., 288., 289. /)
XT = (/ 280. /)

!
! Initialize supersaturation values (in %)
!
DO II=1, SIZE(XS)
   XS(II)=EXP( LOG(10.**(-3.)) + REAL(II) / REAL(SIZE(XS)) * (LOG(10.**2.)-LOG(10.**(-3.))) )
END DO

DO IJ=1, SIZE(XT)
!
!*       2.     Compute Nccn(s) for several supersaturation values
!	        --------------------------------------------------
!
! Get the value of Scrit at Ddry=0.1 micron
!
   XDDRY  = XD
   XMIN   = XD
   XMAX   = XD*10.
   XPREC  = XD/100.
   XW     = 4 * 0.072 * XMV / XAVOGADRO / XBOLTZ / XT(IJ) / XRHOLW
   XSCRIT = ZRIDDR(XMIN,XMAX,XPREC,XDDRY,XKAPPA,XT(IJ))                             ! wet diameter at Scrit
   XSCRIT = (XSCRIT**3-XDDRY**3) * EXP(XW/XSCRIT) / (XSCRIT**3-(1-XKAPPA)*XDDRY**3) ! Saturation ratio at Scrit
   XSCRIT = (XSCRIT - 1.) * 100.                                                    ! Scrit (in %)
!
! Get the XDCRIT values for XS using the approx.
! ln(100*(Sw))~Dcrit^(-3/2) where Sw is in % (Sw=1 for a 1% supersaturation)
!
   XW = XDDRY * XSCRIT**0.66             ! "a" factor in Ddry_crit = a*S**-0.66
   XDCRIT(:) = XW * XS(:)**(-0.66)       ! Ddry_crit for each value of S
!
! Compute Nccn(S) as the incomplete integral of n(D) from 0 to Ddry_crit(S)
!
   DO II=1, SIZE(XS)
      XNCCN(II) = 1- ( 0.5 + SIGN(0.5,XDCRIT(II)-XD) * GAMMA_INC(0.5,(LOG(XDCRIT(II)/XD)/SQRT(2.)/LOG(XSIGMA))**2) )
   END DO
!
!-------------------------------------------------------------------------------
!
!*       3.     Compute C, k, mu, beta, using the Levenberg-Marquardt algorithm
!	        ---------------------------------------------------------------
!
   PARAMS(1:3) = (/ 1., 1., 1000. /)
   IFLAG = 1
   call lmdif1 ( DISTANCE, M, N, PARAMS, FVEC, TOL, INFO )
!
   XLIMIT_FACTOR = gamma(PARAMS(2))*PARAMS(3)**(PARAMS(1)/2)/gamma(1+PARAMS(1)/2)/gamma(PARAMS(2)-PARAMS(1)/2)
   XK            = PARAMS(1)
   XMU           = PARAMS(2)
   XBETA         = PARAMS(3)
!
END DO ! loop on temperatures
!
!-------------------------------------------------------------------------------
!
!*       6.     Functions used to compute Scrit at Ddry=0.1 micron
!   	        --------------------------------------------------
!
CONTAINS
!
!------------------------------------------------------------------------------
!
  FUNCTION ZRIDDR(PX1,PX2,PXACC,XDDRY,XKAPPA,XT)  RESULT(PZRIDDR)
!
!
!!****  *ZRIDDR* - iterative algorithm to find root of a function
!!
!!
!!    PURPOSE
!!    -------
!!       The purpose of this function is to find the root of a given function
!!     the arguments are the brackets bounds (the interval where to find the root)
!!     the accuracy needed and the input parameters of the given function.
!!     Using Ridders' method, return the root of a function known to lie between 
!!     PX1 and PX2. The root, returned as PZRIDDR, will be refined to an approximate
!!     accuracy PXACC.
!! 
!!**  METHOD
!!    ------
!!       Ridders' method
!!
!!    EXTERNAL
!!    --------
!!       FUNCSMAX  
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!      NUMERICAL RECIPES IN FORTRAN 77: THE ART OF SCIENTIFIC COMPUTING 
!!     (ISBN 0-521-43064-X)
!!      Copyright (C) 1986-1992 by Cambridge University Press.
!!      Programs Copyright (C) 1986-1992 by Numerical Recipes Software.
!!
!!    AUTHOR
!!    ------
!!      Frederick Chosson *CERFACS*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original     12/07/07
!!      S.BERTHET        2008 vectorization 
!------------------------------------------------------------------------------
!
!*       0. DECLARATIONS
!
!
IMPLICIT NONE
!
!*       0.1 declarations of arguments and result
!
REAL,    INTENT(INOUT)  :: PX1, PX2, PXACC
REAL,    INTENT(IN)     :: XDDRY, XKAPPA, XT
REAL                    :: PZRIDDR
!
!*       0.2 declarations of local variables
!
!
INTEGER, PARAMETER      :: MAXIT=60
REAL,    PARAMETER      :: UNUSED=0.0 !-1.11e30
REAL                    :: fh,fl, fm,fnew
REAL                    :: s,xh,xl,xm,xnew
INTEGER                 :: j, JL
!
PZRIDDR= 999999.
fl     = DSDD(PX1,XDDRY,XKAPPA,XT)
fh     = DSDD(PX2,XDDRY,XKAPPA,XT)
!
100 if ((fl > 0.0 .and. fh < 0.0) .or. (fl < 0.0 .and. fh > 0.0)) then
      xl         = PX1
      xh         = PX2
      do j=1,MAXIT
         xm     = 0.5*(xl+xh)
         fm = DSDD(xm,XDDRY,XKAPPA,XT)
         s      = sqrt(fm**2-fl*fh)
         if (s == 0.0) then
            GO TO 101
         endif
         xnew  = xm+(xm-xl)*(sign(1.0,fl-fh)*fm/s)
         if (abs(xnew - PZRIDDR) <= PXACC) then
            GO TO 101 
         endif
         PZRIDDR = xnew
         fnew  = DSDD(PZRIDDR,XDDRY,XKAPPA,XT)
         if (fnew == 0.0) then
            GO TO 101
         endif
         if (sign(fm,fnew) /= fm) then
            xl    =xm
            fl=fm
            xh    =PZRIDDR
            fh=fnew
         else if (sign(fl,fnew) /= fl) then
            xh    =PZRIDDR
            fh=fnew
         else if (sign(fh,fnew) /= fh) then
            xl    =PZRIDDR
            fl=fnew
         else if (PX2 .lt. 0.05) then
            PX2 = PX2 + 1.0E-2
!            PRINT*, 'PX2 ALWAYS too small, we put a greater one : PX2 =',PX2
            fh   = DSDD(PX2,XDDRY,XKAPPA,XT)
            go to 100
            STOP
         end if
         if (abs(xh-xl) <= PXACC) then
            GO TO 101 
         endif
      end do
      STOP
   else if (fl == 0.0) then
      PZRIDDR=PX1
   else if (fh == 0.0) then
      PZRIDDR=PX2
   else if (PX2 .lt. 0.05) then
      PX2 = PX2 + 1.0E-2
!      PRINT*, 'PX2 too small, we put a greater one : PX2 =',PX2
      fh   = DSDD(PX2,XDDRY,XKAPPA,XT)
      go to 100
   else
      PZRIDDR=0.0
      go to 101
   end if
!
101 END FUNCTION ZRIDDR
!
!------------------------------------------------------------------------------
!
  FUNCTION DSDD(XD,XDDRY,XKAPPA, XT)  RESULT(DS)
!!
!!    PURPOSE
!!    -------
!!       Derivative of S(D) from Petters and Kreidenweis 2007 (eq. 6) to get Dcrit and Scrit
!!    
!!**  METHOD
!!    ------
!!       This function is called by zriddr
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!       
!!    REFERENCE
!!    ---------
!!    Petters and Kreidenweis, 2007: "A single parameter representation of hygroscopic 
!!             growth and cloud condensation nucleus activity",
!!             ACP, 7, 1961-1971
!!
!!    AUTHOR
!!    ------
!!      Benoit Vie *CNRM*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original     13/11/17
!!
!------------------------------------------------------------------------------
!
!*       0. DECLARATIONS
!
    USE MODD_CST, ONLY : XMV, XAVOGADRO, XBOLTZ, XRHOLW
!
    IMPLICIT NONE
!
!*       0.1 declarations of arguments and result
!
    REAL, INTENT(IN)  :: XD     ! supersaturation is already in no units
    REAL, INTENT(IN)  :: XDDRY  ! supersaturation is already in no units
    REAL, INTENT(IN)  :: XKAPPA ! supersaturation is already in no units
    REAL, INTENT(IN)  :: XT     ! supersaturation is already in no units
!
    REAL              :: DS     ! result
!
!*       0.2 declarations of local variables
!
    REAL              :: XA     ! factor inside the exponential
!    
    XA = 4 * 0.072 * XMV / XAVOGADRO / XBOLTZ / XT / XRHOLW
    DS = (XD**3-XDDRY**3) * (XD**3-(1-XKAPPA)*XDDRY**3) * XA - 3. * XKAPPA * XD**4 * XDDRY**3
    DS = DS * EXP(XA/XD) / (XD**3-(1-XKAPPA)*XDDRY**3)**2
!
END FUNCTION DSDD
!
!-------------------------------------------------------------------------------
!
!*       7.     Functions used to fit the CCN activation spectra with C s**k F()
!   	        ----------------------------------------------------------------
!
  SUBROUTINE DISTANCE(M,N,X,FVEC,IFLAG)
!!
!!    PURPOSE
!!    -------
!!       Derivative of S(D) from Petters and Kreidenweis 2007 (eq. 6) to get Dcrit and Scrit
!!    
!!**  METHOD
!!    ------
!!       This function is called by zriddr
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!       
!!    REFERENCE
!!    ---------
!!    Petters and Kreidenweis, 2007: "A single parameter representation of hygroscopic 
!!             growth and cloud condensation nucleus activity",
!!             ACP, 7, 1961-1971
!!
!!    AUTHOR
!!    ------
!!      Benoit Vie *CNRM*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original     13/11/17
!!
!------------------------------------------------------------------------------
!
!*       0. DECLARATIONS
!
!*       0.1 declarations of arguments and result
!
    integer M
    integer N
    real    X(N)
    real    FVEC(M)
    integer IFLAG
!
!*       0.2 declarations of local variables
!
    integer I
    real    C
    real    ZW, ZW2
!    
    ! print *, "X = ", X
    IF ( ANY(X .LT.0.) .OR. X(1).gt.2*X(2)) THEN
       FVEC(:) = 999999.
    ELSE
       C=gamma(X(2))*X(3)**(X(1)/2)/gamma(1+X(1)/2)/gamma(X(2)-X(1)/2)
       DO I=1, M
          ! XS in "no units", ie XS=0.01 for a 1% suersaturation
          !       ZW= C * (XS(I)/100)**X(1) * HYPGEO(X(2),X(1)/2,X(1)/2+1,X(3),XS(I)/100)
          ZW= C * (XS(I))**X(1) * HYPGEO(X(2),X(1)/2,X(1)/2+1,X(3),XS(I))
!!$       IF (X(3)*(XS(I)/100)**2 .LT. 0.98) THEN
!!$          CALL HYPSER(X(2),X(1)/2,X(1)/2+1,-X(3)*(XS(I)/100)**2,ZW2)
!!$          print *, "args= ", X(2), X(1)/2, X(1)/2+1, -X(3)*(XS(I)/100)**2, " hypser = ", ZW2
!!$          CALL HYPSER(27.288,0.82/2,0.82/2+1,-38726*(0.5/100)**2,ZW2)
!!$          print *, "args= ", 27.288, 0.82/2, 0.82/2+1, -38726*(0.5/100)**2, " hypser = ", ZW2
!!$       END IF
          !       print *, I, XS(I), C, ZW, XNCCN(I)
          IF ( ZW.GT.0. .AND. XNCCN(I).GT.0.) THEN
             FVEC(I) = LOG(ZW) - LOG(XNCCN(I))
          ELSE
             FVEC(I) = 0.
          END IF
          !FVEC(I) = LOG(MAX(ZW,1.E-24)) - LOG(MAX(XNCCN(I),1.E-24))
          !FVEC(I) = ZW - XNCCN(I)
       END DO
    END IF
!    print *, "distance : ", SUM(FVEC*FVEC)
!
  END SUBROUTINE DISTANCE
!
!------------------------------------------------------------------------------
END SUBROUTINE LIMA_INIT_CCN_ACTIVATION_SPECTRUM
