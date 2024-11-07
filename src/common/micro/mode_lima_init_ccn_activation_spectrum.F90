!MNH_LIC Copyright 2007-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-------------------------------------------------------------------------------
MODULE MODE_LIMA_INIT_CCN_ACTIVATION_SPECTRUM
  IMPLICIT NONE
CONTAINS
!     #############################################################
  SUBROUTINE LIMA_INIT_CCN_ACTIVATION_SPECTRUM (CST, HTYPE_CCN,PD,PSIGMA,PLIMIT_FACTOR,PK,PMU,PBETA,PKAPPA)
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
USE MODI_GAMMA_INC
USE MODI_HYPGEO
USE MODI_HYPSER
USE MODD_CST, ONLY:CST_t
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments : 
!
TYPE(CST_t),      INTENT(IN)  :: CST
CHARACTER(LEN=*), INTENT(IN)  :: HTYPE_CCN          ! Aerosol type
REAL,             INTENT(IN)  :: PD             ! Aerosol PSD modal diameter          
REAL,             INTENT(IN)  :: PSIGMA         ! Aerosol PSD width
REAL,             INTENT(OUT) :: PLIMIT_FACTOR  ! C/Naer
REAL,             INTENT(OUT) :: PK             ! k
REAL,             INTENT(OUT) :: PMU            ! mu
REAL,             INTENT(OUT) :: PBETA          ! beta
REAL,             INTENT(OUT) :: PKAPPA         ! kappa
!
!*       0.2   Declarations of local variables :
!
INTEGER, PARAMETER            :: IM = 1000        ! Number of points (S,Nccn) used to fit the spectra
INTEGER, PARAMETER            :: IN = 3          ! Number of parameters to adjust
REAL, DIMENSION(IN)           :: ZPARAMS         ! Parameters to adjust by the LM algorithm (k, mu, beta)
REAL, DIMENSION(IM)           :: ZFVEC           ! Array to store the distance between theoretical and fitted spectra
INTEGER                       :: IFLAG          ! 
INTEGER                       :: IINFO           ! 
REAL                          :: ZTOL = 1.E-16   ! Fit precision required
!
INTEGER                       :: II, IJ         ! Loop indices
!
REAL                          :: ZW             ! 
REAL                          :: ZDDRY = 0.1E-6 ! Dry diameter for which to compute Scrit
REAL                          :: ZSCRIT         ! Scrit for dry diameter ZDDRY
REAL                          :: ZMIN  = 0.1E-6 ! minimum diameter for root search (m)
REAL                          :: ZMAX  = 10.E-6 ! maximum diameter for root search (m)
REAL                          :: ZPREC = 1.E-8  ! precision wanted for root (m)
!
REAL, DIMENSION(IM)           :: ZS             ! saturation ratio (S=1.01 for a 1% supersaturation)
REAL, DIMENSION(IM)           :: ZDCRIT         ! critical diameters (m) for the chosen S values
REAL, DIMENSION(IM)           :: ZNCCN          ! fraction of the aerosols larger than ZDCRIT (ie activable)
REAL, DIMENSION(1)            :: ZT             ! temperature
!
!
!-------------------------------------------------------------------------------
!
!*       1.     Select kappa value based on HTYPE_CCN
!	        ---------------------------------
!
! Kappa values are from Petters and Kreidenweis (2007), table 1.
!
SELECT CASE (HTYPE_CCN)
CASE('NH42SO4','C') ! Ammonium sulfate
   PKAPPA = 0.61
CASE('NH4NO3')      ! Ammonium nitrate
   PKAPPA = 0.67
CASE('NaCl','M')    ! Sea Salt
   PKAPPA = 1.28
CASE('H2SO4')       ! Sulfuric acid
   PKAPPA = 0.90
CASE('NaNO3')       ! Sodium nitrate
   PKAPPA = 0.88
CASE('NaHSO4')      ! Sodium bisulfate
   PKAPPA = 0.91
CASE('Na2SO4')      ! Sodium sulfate
   PKAPPA = 0.80
CASE('NH43HSO42')   ! Letovicite (rare ammonium sulfate mineral)
   PKAPPA = 0.65
CASE('SOA')         ! Secondary organic aerosol (alpha-pinene, beta-pinene)
   PKAPPA = 0.1
CASE DEFAULT
   PKAPPA = 1.
END SELECT
!
!ZT = (/ 270., 271., 272., 273., 274., 275., 276., 277., 278., 279., 280., 281., 282., 283., 284., 285., 286., 287., 288., 289. /)
ZT = (/ 280. /)

!
! Initialize supersaturation values (in %)
!
DO II=1, SIZE(ZS)
   ZS(II)=EXP( LOG(10.**(-3.)) + REAL(II) / REAL(SIZE(ZS)) * (LOG(10.**2.)-LOG(10.**(-3.))) )
END DO

DO IJ=1, SIZE(ZT)
!
!*       2.     Compute Nccn(s) for several supersaturation values
!	        --------------------------------------------------
!
! Get the value of Scrit at Ddry=0.1 micron
!
   ZDDRY  = PD
   ZMIN   = PD
   ZMAX   = PD*10.
   ZPREC  = PD/100.
   ZW     = 4 * 0.072 * CST%XMV / CST%XAVOGADRO / CST%XBOLTZ / ZT(IJ) / CST%XRHOLW
   ZSCRIT = ZRIDDR(ZMIN,ZMAX,ZPREC,ZDDRY,PKAPPA,ZT(IJ))                             ! wet diameter at Scrit
   ZSCRIT = (ZSCRIT**3-ZDDRY**3) * EXP(ZW/ZSCRIT) / (ZSCRIT**3-(1-PKAPPA)*ZDDRY**3) ! Saturation ratio at Scrit
   ZSCRIT = (ZSCRIT - 1.) * 100.                                                    ! Scrit (in %)
!
! Get the ZDCRIT values for ZS using the approx.
! ln(100*(Sw))~Dcrit^(-3/2) where Sw is in % (Sw=1 for a 1% supersaturation)
!
   ZW = ZDDRY * ZSCRIT**0.66             ! "a" factor in Ddry_crit = a*S**-0.66
   ZDCRIT(:) = ZW * ZS(:)**(-0.66)       ! Ddry_crit for each value of S
!
! Compute Nccn(S) as the incomplete integral of n(D) from 0 to Ddry_crit(S)
!
   DO II=1, SIZE(ZS)
      ZNCCN(II) = 1- ( 0.5 + SIGN(0.5,ZDCRIT(II)-PD) * GAMMA_INC(0.5,(LOG(ZDCRIT(II)/PD)/SQRT(2.)/LOG(PSIGMA))**2) )
   END DO
!
!-------------------------------------------------------------------------------
!
!*       3.     Compute C, k, mu, beta, using the Levenberg-Marquardt algorithm
!	        ---------------------------------------------------------------
!
   ZPARAMS(1:3) = (/ 1., 1., 1000. /)
   IFLAG = 1
   call lmdif1 ( DISTANCE, IM, IN, ZPARAMS, ZFVEC, ZTOL, IINFO )
!
   PLIMIT_FACTOR = gamma(ZPARAMS(2))*ZPARAMS(3)**(ZPARAMS(1)/2)/gamma(1+ZPARAMS(1)/2)/gamma(ZPARAMS(2)-ZPARAMS(1)/2)
   PK            = ZPARAMS(1)
   PMU           = ZPARAMS(2)
   PBETA         = ZPARAMS(3)
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
  FUNCTION ZRIDDR(PX1,PX2,PXACC,PDDRY,PKAPPA,PT)  RESULT(PZRIDDR)
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
REAL,    INTENT(IN)     :: PDDRY, PKAPPA, PT
REAL                    :: PZRIDDR
!
!*       0.2 declarations of local variables
!
!
INTEGER, PARAMETER      :: IMAXIT=60
REAL                    :: zfh,zfl, zfm,zfnew
REAL                    :: zs,zxh,zxl,zxm,zxnew
INTEGER                 :: ij
!
PZRIDDR= 999999.
zfl     = DSDD(PX1,PDDRY,PKAPPA,PT)
zfh     = DSDD(PX2,PDDRY,PKAPPA,PT)
!
100 if ((zfl > 0.0 .and. zfh < 0.0) .or. (zfl < 0.0 .and. zfh > 0.0)) then
      zxl         = PX1
      zxh         = PX2
      do ij=1,IMAXIT
         zxm     = 0.5*(zxl+zxh)
         zfm = DSDD(zxm,PDDRY,PKAPPA,PT)
         zs      = sqrt(zfm**2-zfl*zfh)
         if (zs == 0.0) then
            GO TO 101
         endif
         zxnew  = zxm+(zxm-zxl)*(sign(1.0,zfl-zfh)*zfm/zs)
         if (abs(zxnew - PZRIDDR) <= PXACC) then
            GO TO 101 
         endif
         PZRIDDR = zxnew
         zfnew  = DSDD(PZRIDDR,PDDRY,PKAPPA,PT)
         if (zfnew == 0.0) then
            GO TO 101
         endif
         if (sign(zfm,zfnew) /= zfm) then
            zxl    =zxm
            zfl=zfm
            zxh    =PZRIDDR
            zfh=zfnew
         else if (sign(zfl,zfnew) /= zfl) then
            zxh    =PZRIDDR
            zfh=zfnew
         else if (sign(zfh,zfnew) /= zfh) then
            zxl    =PZRIDDR
            zfl=zfnew
         else if (PX2 .lt. 0.05) then
            PX2 = PX2 + 1.0E-2
!            PRINT*, 'PX2 ALWAYS too small, we put a greater one : PX2 =',PX2
            zfh   = DSDD(PX2,PDDRY,PKAPPA,PT)
            go to 100
            STOP
         end if
         if (abs(zxh-zxl) <= PXACC) then
            GO TO 101 
         endif
      end do
      STOP
   else if (zfl == 0.0) then
      PZRIDDR=PX1
   else if (zfh == 0.0) then
      PZRIDDR=PX2
   else if (PX2 .lt. 0.05) then
      PX2 = PX2 + 1.0E-2
!      PRINT*, 'PX2 too small, we put a greater one : PX2 =',PX2
      zfh   = DSDD(PX2,PDDRY,PKAPPA,PT)
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
  FUNCTION DSDD(PD,PDDRY,PKAPPA, PT)  RESULT(ZDS)
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
!
    IMPLICIT NONE
!
!*       0.1 declarations of arguments and result
!
    REAL, INTENT(IN)  :: PD     ! supersaturation is already in no units
    REAL, INTENT(IN)  :: PDDRY  ! supersaturation is already in no units
    REAL, INTENT(IN)  :: PKAPPA ! supersaturation is already in no units
    REAL, INTENT(IN)  :: PT     ! supersaturation is already in no units
!
    REAL              :: ZDS     ! result
!
!*       0.2 declarations of local variables
!
    REAL              :: ZA     ! factor inside the exponential
!    
    ZA = 4 * 0.072 * CST%XMV / CST%XAVOGADRO / CST%XBOLTZ / PT / CST%XRHOLW
    ZDS = (PD**3-PDDRY**3) * (PD**3-(1-PKAPPA)*PDDRY**3) * ZA - 3. * PKAPPA * PD**4 * PDDRY**3
    ZDS = ZDS * EXP(ZA/PD) / (PD**3-(1-PKAPPA)*PDDRY**3)**2
!
END FUNCTION DSDD
!
!-------------------------------------------------------------------------------
!
!*       7.     Functions used to fit the CCN activation spectra with C s**k F()
!   	        ----------------------------------------------------------------
!
  SUBROUTINE DISTANCE(KM,KN,PX,PFVEC,IFLAG)
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
    integer, intent(in) :: KM
    integer, intent(in) :: KN
    real, intent(in) ::    PX(KN)
    real, intent(out) ::    PFVEC(KM)
    integer, intent(inout) :: IFLAG
!
!*       0.2 declarations of local variables
!
    integer II
    real    ZC
    real    ZW
!    
    ! print *, "X = ", X
    IF ( ANY(PX .LT.0.) .OR. PX(1).gt.2*PX(2)) THEN
       PFVEC(:) = 999999.
    ELSE
       ZC=gamma(PX(2))*PX(3)**(PX(1)/2)/gamma(1+PX(1)/2)/gamma(PX(2)-PX(1)/2)
       DO II=1, KM
          ! ZS in "no units", ie ZS=0.01 for a 1% suersaturation
          !       ZW= C * (ZS(I)/100)**X(1) * HYPGEO(X(2),X(1)/2,X(1)/2+1,X(3),ZS(I)/100)
          ZW= ZC * (ZS(II))**PX(1) * HYPGEO(PX(2),PX(1)/2,PX(1)/2+1,PX(3),ZS(II))
!!$       IF (X(3)*(ZS(I)/100)**2 .LT. 0.98) THEN
!!$          CALL HYPSER(X(2),X(1)/2,X(1)/2+1,-X(3)*(ZS(I)/100)**2,ZW2)
!!$          print *, "args= ", X(2), X(1)/2, X(1)/2+1, -X(3)*(ZS(I)/100)**2, " hypser = ", ZW2
!!$          CALL HYPSER(27.288,0.82/2,0.82/2+1,-38726*(0.5/100)**2,ZW2)
!!$          print *, "args= ", 27.288, 0.82/2, 0.82/2+1, -38726*(0.5/100)**2, " hypser = ", ZW2
!!$       END IF
          !       print *, I, ZS(I), C, ZW, ZNCCN(I)
          IF ( ZW.GT.0. .AND. ZNCCN(II).GT.0.) THEN
             PFVEC(II) = LOG(ZW) - LOG(ZNCCN(II))
          ELSE
             PFVEC(II) = 0.
          END IF
          !FVEC(I) = LOG(MAX(ZW,1.E-24)) - LOG(MAX(ZNCCN(I),1.E-24))
          !FVEC(I) = ZW - ZNCCN(I)
       END DO
    END IF
!    print *, "distance : ", SUM(FVEC*FVEC)
!
  END SUBROUTINE DISTANCE
!
!------------------------------------------------------------------------------
END SUBROUTINE LIMA_INIT_CCN_ACTIVATION_SPECTRUM
END MODULE MODE_LIMA_INIT_CCN_ACTIVATION_SPECTRUM
