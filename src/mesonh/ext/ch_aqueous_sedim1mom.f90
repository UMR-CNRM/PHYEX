!MNH_LIC Copyright 2007-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-------------------------------------------------------------------------------
!      ################################
       MODULE MODI_CH_AQUEOUS_SEDIM1MOM
!      ################################
!
INTERFACE
      SUBROUTINE CH_AQUEOUS_SEDIM1MOM (KSPLITR, HCLOUD, OUSECHIC, PTSTEP,  &
                                       PZZ, PRHODREF, PRHODJ, PRRS,        &
                                       PRSS, PRGS, PRRSVS, PSGRSVS, PINPRR )
!
CHARACTER (LEN=4),        INTENT(IN)    :: HCLOUD  ! Cloud parameterization
INTEGER,                  INTENT(IN)    :: KSPLITR ! Current time
REAL,                     INTENT(IN)    :: PTSTEP  ! Time step          
LOGICAL,                  INTENT(IN)    :: OUSECHIC ! flag for ice chem.
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PZZ     ! Height (z)
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRRS    ! Rainwater m.r. source
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRSS    ! Snow m.r. source
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRGS    ! Graupel m.r. source
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRRSVS  ! Rainwater aq. species source
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PSGRSVS ! Precip. ice species source
REAL, DIMENSION(:,:),     INTENT(OUT)   :: PINPRR  ! instantaneaous precip.
!
END SUBROUTINE CH_AQUEOUS_SEDIM1MOM
END INTERFACE
END MODULE MODI_CH_AQUEOUS_SEDIM1MOM
!
!     ######################################################################
      SUBROUTINE CH_AQUEOUS_SEDIM1MOM (KSPLITR, HCLOUD, OUSECHIC, PTSTEP,  &
                                       PZZ, PRHODREF, PRHODJ, PRRS,        &
                                       PRSS, PRGS, PRRSVS, PSGRSVS, PINPRR )
!     ######################################################################
!
!!****  * -  compute the explicit microphysical sources 
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the sedimentation
!!   of chemical species in the raindrops for the Kessler, ICE2, ICE3 and
!!   ICE4 cloud microphysical scheme
!!      The sedimentation rates are computed with a time spliting technique: 
!!    an upstream scheme, written as a difference of non-advective fluxes. 
!!    This source term is added to the next coming time step (split-implicit 
!!    process). see rain_ice.f90
!!
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!      None
!!     
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_PARAMETERS
!!          JPHEXT       : Horizontal external points number
!!          JPVEXT       : Vertical external points number
!!      Module MODD_CONF :
!!          CCONF configuration of the model for the first time step
!!
!!    REFERENCE
!!    ---------
!!      Book1 of the documentation ( routine CH_AQUEOUS_SEDIM1MOM )
!!
!!    AUTHOR
!!    ------
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!    Original    22/07/07
!!    04/11/08 (M Leriche) add ICE3    
!!    17/09/10 (M Leriche) add LUSECHIC flag
!!    J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!!    16/12/15 (M Leriche) compute instantaneous rain at the surface
!  P. Wautelet 12/02/2019: bugfix: ZRR_SEDIM was not initialized everywhere
!  P. Wautelet 26/04/2019: replace non-standard FLOAT function by REAL function
!  P. Wautelet 28/05/2019: move COUNTJV function to tools.f90
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAMETERS,      ONLY : JPHEXT, JPVEXT
USE MODD_CONF
USE MODD_CST,             ONLY : XRHOLW
USE MODD_CLOUDPAR,        ONLY : VCEXVT=>XCEXVT, XCRS, XCEXRS
USE MODD_RAIN_ICE_DESCR,  ONLY : WCEXVT=>XCEXVT, WRTMIN=>XRTMIN
USE MODD_RAIN_ICE_PARAM,  ONLY : XFSEDR, XEXSEDR, &
                                 XFSEDS, XEXSEDS, &
                                 XFSEDG, XEXSEDG

use mode_tools,           only: Countjv
use mode_tools_ll,        only: GET_INDICE_ll

IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
CHARACTER (LEN=4),        INTENT(IN)    :: HCLOUD  ! Cloud parameterization
INTEGER,                  INTENT(IN)    :: KSPLITR ! Current time
REAL,                     INTENT(IN)    :: PTSTEP  ! Time step          
LOGICAL,                  INTENT(IN)    :: OUSECHIC ! flag for ice chem.
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PZZ     ! Height (z)
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRRS    ! Rainwater m.r. source
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRSS    ! Snow m.r. source
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRGS    ! Graupel m.r. source
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRRSVS  ! Rainwater aq. species source
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PSGRSVS ! Precip. ice species source
REAL, DIMENSION(:,:),     INTENT(OUT)   :: PINPRR  ! instantaneaous precip.
!
!*       0.2   Declarations of local variables :
!
INTEGER :: JK,JI,JJ            ! Vertical loop index for the rain sedimentation 
INTEGER :: JN            ! Temporal loop index for the rain sedimentation
INTEGER :: IIB           !  Define the domain where is 
INTEGER :: IIE           !  the microphysical sources have to be computed
INTEGER :: IJB           ! 
INTEGER :: IJE           !
INTEGER :: IKB           ! 
INTEGER :: IKE           !
!
REAL    :: ZTSPLITR      ! Small time step for rain sedimentation
!
INTEGER :: ISEDIMR, ISEDIMS, ISEDIMG ! Case number of sedimentation
LOGICAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)) &
                                :: GSEDIMR ! where to compute the SED processes
LOGICAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)) &
                                :: GSEDIMS ! where to compute the SED processes
LOGICAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)) &
                                :: GSEDIMG ! where to compute the SED processes
REAL,    DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3))   &
                                :: ZRRS     ! rainwater m.r.source phys.tendency
REAL,    DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3))   &
                                :: ZRSS     ! snow m.r.source phys.tendency
REAL,    DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3))   &
                                :: ZRGS     ! graupel m.r.source phys.tendency
REAL,    DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3))   &
                                :: ZW     ! work array
REAL,    DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3))   &
                                :: ZWSED  ! sedimentation fluxes
REAL,    DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3))   &
                                :: ZZRRS  ! Rainwater m.r. source phys.tendency *dt
REAL,    DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3))   &
                                :: ZZRSS  ! Snow m.r. source phys.tendency *dt
REAL,    DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3))   &
                                :: ZZRGS  ! Graupel m.r. source phys.tendency *dt
REAL,    DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3))   &
                                :: ZRR_SEDIM       ! Drain/Dt sur ZTSPLIT
REAL,    DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3))   &
                                :: ZSV_SEDIM_FACTR  ! Cumul des Dsv/DT
REAL,    DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3))   &
                                :: ZSV_SEDIM_FACTS  ! Cumul des Dsv/DT
REAL,    DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3))   &
                                :: ZSV_SEDIM_FACTG  ! Cumul des Dsv/DT
REAL, DIMENSION(:), ALLOCATABLE :: ZZZRRS  ! Rainwater m.r. source
REAL, DIMENSION(:), ALLOCATABLE :: ZZZRSS  ! Snow m.r. source
REAL, DIMENSION(:), ALLOCATABLE :: ZZZRGS  ! Graupel m.r. source
!
REAL, DIMENSION(:), ALLOCATABLE :: ZRHODREF, & ! RHO Dry REFerence
                                   ZZW         ! Work array
REAL, DIMENSION(7), SAVE        :: Z_XRTMIN
!
REAL                            :: ZVTRMAX, ZT
LOGICAL, SAVE                   :: GSFIRSTCALL = .TRUE.
REAL,    SAVE                   :: ZFSEDR, ZEXSEDR, ZCEXVT
!
INTEGER , DIMENSION(SIZE(GSEDIMR)) :: IR1,IR2,IR3 ! Used to replace the COUNT
INTEGER , DIMENSION(SIZE(GSEDIMS)) :: IS1,IS2,IS3 ! Used to replace the COUNT
INTEGER , DIMENSION(SIZE(GSEDIMG)) :: IG1,IG2,IG3 ! Used to replace the COUNT
INTEGER                           :: JL       ! and PACK intrinsics
!  
!-------------------------------------------------------------------------------
!
!*       1.     COMPUTE THE LOOP BOUNDS
!   	        -----------------------
!
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
IKB=1+JPVEXT
IKE=SIZE(PZZ,3) - JPVEXT
PINPRR(:,:) = 0. ! initialize instantaneous precip.
!
!-------------------------------------------------------------------------------
!
!!*       2.     TRANSFORMATION INTO PHYSICAL TENDENCIES
!               ---------------------------------------
!
ZRRS(:,:,:)  = PRRS(:,:,:) / PRHODJ(:,:,:)
IF (HCLOUD(1:3) == 'ICE') THEN
  ZRSS(:,:,:)  = PRSS(:,:,:) / PRHODJ(:,:,:)
  ZRGS(:,:,:)  = PRGS(:,:,:) / PRHODJ(:,:,:)
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       3.     COMPUTE THE SEDIMENTATION (RS) SOURCE
!	        -------------------------------------
!
!*       3.1    Initialize some constants
!  
firstcall : IF (GSFIRSTCALL) THEN
  GSFIRSTCALL = .FALSE.
  SELECT CASE ( HCLOUD)
    CASE('KESS')
      ZVTRMAX = 20.                          
    CASE('ICE3')
      ZVTRMAX = 10.
    CASE('ICE4')
      ZVTRMAX = 40.
  END SELECT
!
  SELECT CASE ( HCLOUD )  ! constants for rain sedimentation
    CASE('KESS')
      Z_XRTMIN(2:3) = 1.0E-20     ! Default values
      ZFSEDR  = XCRS
      ZEXSEDR = XCEXRS
      ZCEXVT  = VCEXVT
    CASE('ICE3','ICE4')
      Z_XRTMIN(1:SIZE(WRTMIN)) = WRTMIN ! Values given in ICEx schemes
      ZFSEDR  = XFSEDR
      ZEXSEDR = XEXSEDR
      ZCEXVT  = WCEXVT
  END SELECT
END IF firstcall
!
!*       3.2    time splitting loop initialization
!
ZTSPLITR = PTSTEP / REAL(KSPLITR)       ! Small time step
!
!*       3.3    compute the fluxes
!
ZSV_SEDIM_FACTR(:,:,:) = 1.0
ZZRRS(:,:,:) = ZRRS(:,:,:) * PTSTEP
IF (HCLOUD(1:3) == 'ICE') THEN
  ZZRSS(:,:,:) = ZRSS(:,:,:) * PTSTEP
  ZZRGS(:,:,:) = ZRGS(:,:,:) * PTSTEP
  ZSV_SEDIM_FACTS(:,:,:) = 1.0
  ZSV_SEDIM_FACTG(:,:,:) = 1.0
ENDIF
DO JN = 1 , KSPLITR
  IF( JN==1 ) THEN
    ZW(:,:,:) = 0.0
    DO JK = IKB , IKE-1
      ZW(:,:,JK) =ZTSPLITR*2./(PRHODREF(:,:,JK)*(PZZ(:,:,JK+2)-PZZ(:,:,JK)))
    END DO
    ZW(:,:,IKE)  =ZTSPLITR/(PRHODREF(:,:,IKE)*(PZZ(:,:,IKE+1)-PZZ(:,:,IKE)))
  END IF
!
!*       3.3.1   for rain
!
  GSEDIMR(:,:,:) = .FALSE.
  GSEDIMR(IIB:IIE,IJB:IJE,IKB:IKE) = ZZRRS(IIB:IIE,IJB:IJE,IKB:IKE) > Z_XRTMIN(3)
  ISEDIMR = COUNTJV( GSEDIMR(:,:,:),IR1(:),IR2(:),IR3(:))
!
  IF ( ISEDIMR >= 1 ) THEN
    ALLOCATE(ZZZRRS(ISEDIMR)) 
    ALLOCATE(ZRHODREF(ISEDIMR))
    DO JL=1,ISEDIMR
      ZZZRRS(JL) = ZZRRS(IR1(JL),IR2(JL),IR3(JL))
      ZRHODREF(JL) =  PRHODREF(IR1(JL),IR2(JL),IR3(JL))
    ENDDO
    ALLOCATE(ZZW(ISEDIMR)) ; ZZW(:) = 0.0
!
    ZZW(:) = ZFSEDR * ZZZRRS(:)**(ZEXSEDR) * ZRHODREF(:)**(ZEXSEDR-ZCEXVT)
    ZWSED(:,:,:) = UNPACK( ZZW(:),MASK=GSEDIMR(:,:,:),FIELD=0.0 )
    ZRR_SEDIM(:,:,:) = 0.0
    DO JK = IKB , IKE
      ZRR_SEDIM(:,:,JK) = ZW(:,:,JK)*(ZWSED(:,:,JK+1)-ZWSED(:,:,JK))
    END DO
    ZZRRS(:,:,:) = ZZRRS(:,:,:) + ZRR_SEDIM(:,:,:)
    PINPRR(:,:) = PINPRR(:,:) + ZWSED(:,:,IKB)/XRHOLW/KSPLITR
!
    ZZW(:) = ZFSEDR * ZZZRRS(:)**(ZEXSEDR-1.0) * ZRHODREF(:)**(ZEXSEDR-ZCEXVT)
    ZWSED(:,:,:) = UNPACK( ZZW(:),MASK=GSEDIMR(:,:,:),FIELD=0.0 )
    ZRR_SEDIM(:,:,:) = 0.0
    DO JK = IKB , IKE
      ZRR_SEDIM(:,:,JK) = ZW(:,:,JK)*(ZWSED(:,:,JK+1)-ZWSED(:,:,JK))
    END DO
    DEALLOCATE(ZRHODREF)
    DEALLOCATE(ZZZRRS)
    DEALLOCATE(ZZW)
    ZSV_SEDIM_FACTR(:,:,:) =   ZSV_SEDIM_FACTR(:,:,:) * (1.0 + ZRR_SEDIM(:,:,:))
!!                       (1.0 + ZRR_SEDIM(:,:,:)/MAX(ZZRRS(:,:,:),XRTMIN_AQ))
  END IF      
  IF (HCLOUD == 'KESS') EXIT
!
!*       3.3.1   for iced precip.hydrometeors
!
  GSEDIMS(:,:,:) = .FALSE.
  GSEDIMG(:,:,:) = .FALSE.
  GSEDIMS(IIB:IIE,IJB:IJE,IKB:IKE) = ZZRSS(IIB:IIE,IJB:IJE,IKB:IKE) > Z_XRTMIN(5)
  GSEDIMG(IIB:IIE,IJB:IJE,IKB:IKE) = ZZRGS(IIB:IIE,IJB:IJE,IKB:IKE) > Z_XRTMIN(6)
  ISEDIMS = COUNTJV( GSEDIMS(:,:,:),IS1(:),IS2(:),IS3(:))
  ISEDIMG = COUNTJV( GSEDIMG(:,:,:),IG1(:),IG2(:),IG3(:))
! for snow
  IF ( ISEDIMS >= 1) THEN
    ALLOCATE(ZZZRSS(ISEDIMS))
    ALLOCATE(ZRHODREF(ISEDIMS))
    DO JL=1,ISEDIMS
      ZZZRSS(JL) = ZZRSS(IS1(JL),IS2(JL),IS3(JL))
      ZRHODREF(JL) =  PRHODREF(IS1(JL),IS2(JL),IS3(JL))
    ENDDO
    ALLOCATE(ZZW(ISEDIMS)) ; ZZW(:) = 0.0
!
    ZZW(:) = XFSEDS * ZZZRSS(:)**(XEXSEDS) * ZRHODREF(:)**(XEXSEDS-ZCEXVT)
    ZWSED(:,:,:) = UNPACK( ZZW(:),MASK=GSEDIMS(:,:,:),FIELD=0.0 )
    ZRR_SEDIM(:,:,:) = 0.0
    DO JK = IKB , IKE
      ZRR_SEDIM(:,:,JK) = ZW(:,:,JK)*(ZWSED(:,:,JK+1)-ZWSED(:,:,JK))
    END DO
    ZZRSS(:,:,:) = ZZRSS(:,:,:) + ZRR_SEDIM(:,:,:)
!
    ZZW(:) = XFSEDS * ZZZRSS(:)**(XEXSEDS-1.0) * ZRHODREF(:)**(XEXSEDS-ZCEXVT)
    ZWSED(:,:,:) = UNPACK( ZZW(:),MASK=GSEDIMS(:,:,:),FIELD=0.0 )
    ZRR_SEDIM(:,:,:) = 0.0
    DO JK = IKB , IKE
      ZRR_SEDIM(:,:,JK) = ZW(:,:,JK)*(ZWSED(:,:,JK+1)-ZWSED(:,:,JK))
    END DO
    DEALLOCATE(ZRHODREF)
    DEALLOCATE(ZZZRSS)
    DEALLOCATE(ZZW)
    ZSV_SEDIM_FACTS(:,:,:) =   ZSV_SEDIM_FACTS(:,:,:) * (1.0 + ZRR_SEDIM(:,:,:))
  ENDIF
! for graupel
  IF ( ISEDIMG >= 1) THEN
    ALLOCATE(ZZZRGS(ISEDIMG))
    ALLOCATE(ZRHODREF(ISEDIMG))
    DO JL=1,ISEDIMG
      ZZZRGS(JL) = ZZRGS(IG1(JL),IG2(JL),IG3(JL))
      ZRHODREF(JL) =  PRHODREF(IG1(JL),IG2(JL),IG3(JL))
    ENDDO
    ALLOCATE(ZZW(ISEDIMG)) ; ZZW(:) = 0.0
!
    ZZW(:) = XFSEDG * ZZZRGS(:)**(XEXSEDG) * ZRHODREF(:)**(XEXSEDG-ZCEXVT)
    ZWSED(:,:,:) = UNPACK( ZZW(:),MASK=GSEDIMG(:,:,:),FIELD=0.0 )
    ZRR_SEDIM(:,:,:) = 0.0
    DO JK = IKB , IKE
      ZRR_SEDIM(:,:,JK) = ZW(:,:,JK)*(ZWSED(:,:,JK+1)-ZWSED(:,:,JK))
    END DO
    ZZRGS(:,:,:) = ZZRGS(:,:,:) + ZRR_SEDIM(:,:,:)
!
    ZZW(:) = XFSEDG * ZZZRGS(:)**(XEXSEDG-1.0) * ZRHODREF(:)**(XEXSEDG-ZCEXVT)
    ZWSED(:,:,:) = UNPACK( ZZW(:),MASK=GSEDIMG(:,:,:),FIELD=0.0 )
    ZRR_SEDIM(:,:,:) = 0.0
    DO JK = IKB , IKE
      ZRR_SEDIM(:,:,JK) = ZW(:,:,JK)*(ZWSED(:,:,JK+1)-ZWSED(:,:,JK))
    END DO
    DEALLOCATE(ZRHODREF)
    DEALLOCATE(ZZZRGS)
    DEALLOCATE(ZZW)
    ZSV_SEDIM_FACTG(:,:,:) =   ZSV_SEDIM_FACTG(:,:,:) * (1.0 + ZRR_SEDIM(:,:,:))
  ENDIF
END DO
!
! Apply the rain sedimentation rate to the WR_xxx aqueous species
DO JL= 1, SIZE(PRRSVS,4)
  PRRSVS(:,:,:,JL) = MAX( 0.0,ZSV_SEDIM_FACTR(:,:,:)*PRRSVS(:,:,:,JL) )
ENDDO
!ice phase
IF (OUSECHIC) THEN
  DO JL= 1, SIZE(PSGRSVS,4)
    PSGRSVS(:,:,:,JL) = MAX( 0.0, &
                      ((ZSV_SEDIM_FACTS(:,:,:)+ZSV_SEDIM_FACTG(:,:,:))/2.) &
                                        *PSGRSVS(:,:,:,JL) )
  ENDDO
ENDIF
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE CH_AQUEOUS_SEDIM1MOM
