!MNH_LIC Copyright 1995-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ######spl
 MODULE MODI_SET_MSK
!####################
!
INTERFACE
!
SUBROUTINE SET_MSK(PRT,PRHODREF,OBU_MSK) 
!
REAL , DIMENSION (:,:,:,:),INTENT(IN) :: PRT
REAL , DIMENSION (:,:,:),INTENT(IN)   :: PRHODREF
LOGICAL , DIMENSION (:,:,:),INTENT(OUT)  :: OBU_MSK
!
END SUBROUTINE SET_MSK  
!
END INTERFACE
!
END MODULE MODI_SET_MSK
!
!     ######spl
      SUBROUTINE SET_MSK(PRT,PRHODREF,OBU_MSK) 
!     ###############################
!
!!****SET_MSK** -routine to define the mask based on SET_MASK
!!                           
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to test the occurence or not of the
!     different criteria, used to compute the budgets. It also updates the 
!     number of occurence of the different criteria.
!
!!**  METHOD
!!    ------
!!      According to each criterion associated to one zone, the mask is
!!    set to TRUE at each point where the criterion is confirmed, at each 
!!    time step of the model.
!!
!!
!!    EXTERNAL
!!    --------
!!       NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!      Book2 of MESO-NH documentation (routine BUDGET)
!!
!!
!!    AUTHOR
!!    ------
!!	J. Nicolau       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    27/02/95
!!      T.Montmerle 15/07/96  Computation of masks for convective and stratiform parts
!!      Biju Thomas 29/03/99  Identified nonprecipitating convective cells and only
!!                            precipitating anvils as stratiform part
!!      O. Caumont  09/04/08  Use in RADAR_SIMULATOR
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_FIELD_n
USE MODD_RAIN_ICE_PARAM_n , ONLY : XFSEDR,XEXSEDR
USE MODD_RAIN_ICE_DESCR_n , ONLY : XCEXVT
USE MODD_CST  , ONLY : XRHOLW
USE MODD_PARAMETERS
USE MODD_CONF
USE MODE_ll
USE MODD_LUNIT, ONLY : TLUOUT0
USE MODD_ARGSLIST_ll, ONLY : LIST_ll
!
!
IMPLICIT NONE
!  
!  
!*       0.1   Declarations of arguments :
!
REAL , DIMENSION (:,:,:,:),INTENT(IN) :: PRT
REAL , DIMENSION (:,:,:),INTENT(IN)   :: PRHODREF
LOGICAL , DIMENSION (:,:,:),INTENT(OUT)  :: OBU_MSK
!
!*       0.2   Declarations of local variables :
!
INTEGER                    :: IIB,IJB       ! Lower bounds and Upper bounds
INTEGER                    :: IIE,IJE       ! of the physical sub-domain 
INTEGER                    :: IKB,IKE       ! in x, y and z directions
INTEGER                    :: IIU,IJU!,IKU
!
REAL,DIMENSION(:,:,:), ALLOCATABLE  :: ZMASK     ! signature de l'insertion 
                                                 ! dans un masque (0 ou 1.)
REAL,DIMENSION(:,:),   ALLOCATABLE  :: ZCONVECT  ! signature du domaine convectif
REAL,DIMENSION(:,:),   ALLOCATABLE  :: ZSURFPP   ! precipitation au sol 
REAL,DIMENSION(:,:),   ALLOCATABLE  :: ZMAXWATER ! teneur maximale en eau 
                                                 ! recensee sur la verticale
REAL,DIMENSION(:,:),   ALLOCATABLE  :: ZMIMX,ZMIPX ! I,I+1 and I,I-1 precipitation sums
REAL,DIMENSION(:,:),   ALLOCATABLE  :: ZMEANX_MY,ZMEANX_PY ! J,J+1 and J,J-1 precipitation sums
REAL,DIMENSION(:,:),   ALLOCATABLE  :: ZMEANX, ZMEANXY
REAL       ::   ZAVER_PR,ZREPSILON,ZTOTWATER,ZREPSILON1
REAL       ::   ZCRS,ZCEXRS,ZCEXVT,ZREPSILON2,ZREPSILON3            
INTEGER    ::   I,J,JILOOP,JJLOOP,JKLOOP
INTEGER :: ILUOUT0
INTEGER :: IRESP
INTEGER :: IBUIL,IBUJL,IBUIH,IBUJH
!INTEGER :: IBUSIL,IBUSJL,IBUSIH,IBUSJH
!INTEGER                :: IINFO_ll       ! return code of parallel routine
!TYPE(LIST_ll), POINTER :: TZFIELDS_ll    ! list of fields to exchange
!-------------------------------------------------------------------------------
!
ILUOUT0 = TLUOUT0%NLU
!
!*       1.    COMPUTES THE PHYSICAL SUBDOMAIN BOUNDS
!              ---------------------------------------
!
IKB = 1 + JPVEXT
IKE = SIZE(PRT,3) - JPVEXT
IIU = SIZE(PRT,1)
IJU = SIZE(PRT,2)
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
!
!               ----------------------
ALLOCATE( ZMASK(IIU,IJU,4) )
ALLOCATE( ZSURFPP(IIU,IJU) )
ALLOCATE(ZMIMX(IIU,IJU),ZMIPX(IIU,IJU),ZMEANX(IIU,IJU))
ALLOCATE(ZMEANX_MY(IIU,IJU),ZMEANX_PY(IIU,IJU),ZMEANXY(IIU,IJU))
ALLOCATE( ZCONVECT(IIU,IJU) )
ALLOCATE( ZMAXWATER(IIU,IJU) )
!
!*	 2.     DEFINITION OF THE MASK
!           ----------------------
!  initialization to FALSE on the extended subdomain
OBU_MSK(:,:,:)=.FALSE.
ZMASK(:,:,:)=0.
ZSURFPP(:,:)=0.
ZCONVECT(:,:)=0.
ZMAXWATER(:,:)=0.
ZREPSILON=5.E-6
ZREPSILON1=5.E-4
ZREPSILON2=5.0
ZREPSILON3=5.E-6
ZAVER_PR=0.

!**********************************************************************
!                      CAUTION: Definition of parameters
!              depends on the model configuration WARM or COLD      
!              -----------------------------------------------

!**********************************************************************
!partie a activer pour le cas chaud, en activant USE MODD_CLOUDPAR et en 
!desactivant USE MODD_RAIN_ICE_PARAM et USE MODD_RAIN_ICE_DESCR qui servent
!au cas froid. En activant tout, XCEXVT est defini deux fois, donc une fois
!de trop.
!**********************************************************************
!IF (CCLOUD == 'REVE' .OR. CCLOUD == 'KESS' .OR. CCLOUD == 'KES2') THEN
!  ZCRS=XCRS
!  ZCEXRS=XCEXRS
!  ZCEXVT=XCEXVT
!ELSE IF (CCLOUD == 'ICE3') THEN
!**********************************************************************

  ZCRS=XFSEDR
  ZCEXRS=XEXSEDR
  ZCEXVT=XCEXVT
!END IF

!            Total solid and liquid water (qr+qc+qs+qi+qg) (= cloudy area)
!            -------------------------------------------------------------
 
DO JKLOOP=IKB,IKE
  DO JJLOOP=IJB,IJE
     DO JILOOP=IIB,IIE
       ZTOTWATER = PRT(JILOOP,JJLOOP,JKLOOP,2) &
                  +PRT(JILOOP,JJLOOP,JKLOOP,3) &
                  +PRT(JILOOP,JJLOOP,JKLOOP,4) &
                  +PRT(JILOOP,JJLOOP,JKLOOP,5) &
                  +PRT(JILOOP,JJLOOP,JKLOOP,6)  
       ZMAXWATER(JILOOP,JJLOOP)=MAX(ZMAXWATER(JILOOP,JJLOOP),ZTOTWATER)
     END DO
  END DO
END DO

!           Computation of ground precipitation
!           -----------------------------------
 
!      Precipitation (mm/h)
ZSURFPP(IIB:IIE,IJB:IJE)=ZCRS*PRT(IIB:IIE,IJB:IJE,IKB,3)**ZCEXRS  &
     *PRHODREF(IIB:IIE,IJB:IJE,IKB)**(ZCEXRS-ZCEXVT)*3.6E6/XRHOLW

!    Lateral Boundaries for Precipitation
!  (cyclic case in Y-direction, OPEN in X-direction)
    ZSURFPP(1,IJB:IJE)=ZSURFPP(IIB,IJB:IJE)
    ZSURFPP(IIU,IJB:IJE)=ZSURFPP(IIE,IJB:IJE)
    ZSURFPP(1:IIU,1)=ZSURFPP(1:IIU,IJB)
    ZSURFPP(1:IIU,IJU)=ZSURFPP(1:IIU,IJE)

!
!       Predefinition of the Convective region criteria
!       ------------------------------------------------
ZMIPX(:,:)=0.
ZMIMX(:,:)=0.
ZMEANX(:,:)=0.
!
ZMIPX(1:IIU-1,:)=ZSURFPP(1:IIU-1,:)+ZSURFPP(2:IIU,:)
ZMIMX(2:IIU,:)=ZSURFPP(2:IIU,:)+ZSURFPP(1:IIU-1,:)

DO J=IJB+1,IJE-1
   DO I=3,IIE-1
     ZAVER_PR=(SUM(ZSURFPP(I-2:I+2,J-2:J+2))-ZSURFPP(I,J))/24.

!                        threshold at  4 mm/h
     IF(ZSURFPP(I,J) >= MAX(4.,2.*ZAVER_PR)  &
     .AND.(ZMAXWATER(I,J) >= ZREPSILON))  ZCONVECT(I-1:I+1,J-1:J+1)=1.
     IF(ZSURFPP(I,J) >= 20.) ZCONVECT(I,J)=1.
     IF(ZMAXWATER(I,J) >= ZREPSILON)THEN
       DO JKLOOP=2,IKE
           IF(PRT(I,J,JKLOOP,2) >= ZREPSILON1) ZCONVECT(I,J)=1.
           IF(XWT(I,J,JKLOOP) >= ZREPSILON2) ZCONVECT(I,J)=1.
       END DO
     END IF
  END DO
END DO 

!------------------------------------------
!*  MASK Definition
!------------------------------------------
IBUIL=IIB+1
IBUIH = IIE-1
IBUJL = IJB+1
IBUJH = IJE-1
DO JILOOP=IBUIL,IBUIH
  DO JJLOOP=IBUJL,IBUJH
!------------------------------------------
!*   Zone 1: Convective Zone
!------------------------------------------
     ZMASK(JILOOP,JJLOOP,1)=ZCONVECT(JILOOP,JJLOOP) 
!------------------------------------------
!*   Zone 2: Stratiforme Zone
!------------------------------------------
     IF (ZMAXWATER(JILOOP,JJLOOP) >= 10.*ZREPSILON.AND.ZMASK(JILOOP,JJLOOP,1)/=1.) THEN 
         DO JKLOOP=IKB,IKE
            IF(PRT(JILOOP,JJLOOP,JKLOOP,3) >= ZREPSILON3)  ZMASK(JILOOP,JJLOOP,2)=1.
         END DO
     END IF
!------------------------------------------
!*   Zone 3: Clear air Zone
!------------------------------------------
     IF (ZMASK(JILOOP,JJLOOP,1)/=1. .AND. ZMASK(JILOOP,JJLOOP,2)/=1.)  ZMASK(JILOOP,JJLOOP,3)=1.
!------------------------------------------
!*   Zone 4: Total Domain
!------------------------------------------
     ZMASK(JILOOP,JJLOOP,4)=1.
 
  END DO
END DO
!
!-----------------------------------------------------------------------
!

OBU_MSK(IIB:IIE,IJB:IJE,:)=ZMASK(IIB:IIE,IJB:IJE,:)>0.8


!
!*	 2.     INCREASE IN SURFACE ARRAY
!               -------------------------
!
DEALLOCATE( ZMASK )
DEALLOCATE( ZCONVECT )
DEALLOCATE( ZSURFPP )
DEALLOCATE( ZMAXWATER )
DEALLOCATE(ZMIMX,ZMIPX,ZMEANX)
DEALLOCATE(ZMEANX_MY,ZMEANX_PY,ZMEANXY)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SET_MSK
