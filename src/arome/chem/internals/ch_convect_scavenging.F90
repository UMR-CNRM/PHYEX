!     ######spl
      SUBROUTINE CH_CONVECT_SCAVENGING( KLON, KLEV, KCH, PCH1, PCH1C,        &
                                         KDPL, KPBL, KLCL, KCTL, KLFS, KDBL, &
                                         PUMF, PUER, PUDR, PDMF, PDER, PDDR, &
                                         PTIMEC, PDXDY, PMIXF, PLMASS, PWSUB,&
                                         KFTSTEPS,                           &
                                         PURC, PURR, PURI, PURS, PUTT, PPRESS,&
                                         PRHODREF                            )
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!     #######################################################################
!
!!**** Compute modified soluble  chemical tracer values due to convective
!!     precipitations + transport
!!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to determine the final adjusted
!!      environmental values of the chemical tracers
!!      The final convective tendencies can then be evaluated in the main
!!      routine DEEP_CONVECT by (PCH1C-PCH1)/PTIMEC
!!
!!
!!**  METHOD
!!    ------
!!      Identical to the computation of the conservative variables in the
!!      main deep convection code
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST
!!          XG                 ! gravity constant
!!
!!     Module MODD_CONVPAREXT
!!          JCVEXB, JCVEXT     ! extra levels on the vertical boundaries
!!
!!    AUTHOR
!!    ------
!!      C. MARI       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      Original    10/04/00
!!      P. Tulet    25/04/05  Aerosols/ Dust scavenging
!!
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
USE MODD_PARAMETERS
USE MODD_CONVPAREXT
USE MODD_CH_M9
USE MODD_CH_CONST, ONLY : XSREALHENRYVAL
USE MODD_NSV,     ONLY : NSV_CHEMBEG, NSV_CHEMEND, &
                         NSV_AERBEG, NSV_AEREND, &
                         NSV_DSTBEG, NSV_DSTEND, &
                         NSV_SLTBEG, NSV_SLTEND
!USE MODD_SALT,    ONLY : LSALT, NMODE_SLT, JPSALTORDER
!USE MODD_DUST,    ONLY : LDUST, NMODE_DST, JPDUSTORDER
!USE MODD_CH_AEROSOL, ONLY : JPMODE, LORILAM, NSP, NSOA, NCARB
USE MODD_SALT
USE MODD_DUST
USE MODD_CH_AEROSOL
USE MODD_CSTS_SALT
USE MODD_CSTS_DUST
USE MODE_DUST_PSD
USE MODE_SALT_PSD
USE MODE_AERO_PSD
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER,                INTENT(IN) :: KLON     ! horizontal dimension
INTEGER,                INTENT(IN) :: KLEV     ! vertical dimension
INTEGER,                INTENT(IN) :: KCH      ! number of passive tracers
!
REAL,DIMENSION(KLON,KLEV,KCH),INTENT(IN) :: PCH1 ! grid scale tracer concentr.
REAL,DIMENSION(KLON,KLEV,KCH),INTENT(OUT):: PCH1C! conv adjusted tracer concntr.
!
INTEGER, DIMENSION(KLON), INTENT(IN) :: KDPL   ! index for departure level
INTEGER, DIMENSION(KLON), INTENT(IN) :: KPBL   ! index for top of source layer
INTEGER, DIMENSION(KLON), INTENT(IN) :: KLCL   ! index lifting condens. level
INTEGER, DIMENSION(KLON), INTENT(IN) :: KCTL   ! index for cloud top level
INTEGER, DIMENSION(KLON), INTENT(IN) :: KLFS   ! index for level of free sink
INTEGER, DIMENSION(KLON), INTENT(IN) :: KDBL   ! index for downdraft base level
!
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PUMF ! updraft mass flux (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PUER ! updraft entrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PUDR ! updraft detrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PDMF ! downdraft mass flux (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PDER ! downdraft entrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PDDR ! downdraft detrainment (kg/s)
!
REAL, DIMENSION(KLON),     INTENT(IN) :: PTIMEC! convection time step
REAL, DIMENSION(KLON),     INTENT(IN) :: PDXDY ! grid area (m^2)
REAL, DIMENSION(KLON),     INTENT(IN) :: PMIXF ! mixed fraction at LFS
REAL, DIMENSION(KLON,KLEV),INTENT(IN) :: PLMASS! mass of model layer (kg)
REAL, DIMENSION(KLON,KLEV),INTENT(IN) :: PWSUB ! envir. compensating subsidence(Pa/s)
INTEGER,                INTENT(IN) :: KFTSTEPS  ! maximum fractional time steps
REAL, DIMENSION(KLON,KLEV),   INTENT(IN) :: PURC      ! microphysical
REAL, DIMENSION(KLON,KLEV),   INTENT(IN) :: PURR      ! reservoirs in the
REAL, DIMENSION(KLON,KLEV),   INTENT(IN) :: PURI      ! updraft
REAL, DIMENSION(KLON,KLEV),   INTENT(IN) :: PURS      !
REAL, DIMENSION(KLON,KLEV),   INTENT(IN) :: PUTT      ! updraft temperature (K)
REAL, DIMENSION(KLON,KLEV),   INTENT(IN) :: PPRESS    ! Pa
REAL, DIMENSION(KLON,KLEV),   INTENT(IN) :: PRHODREF
!
!*       0.2   Declarations of local variables :
!
INTEGER :: INCH1          ! number of chemical tracers
INTEGER :: IKB,IKE
INTEGER :: IKS            ! vertical dimension
INTEGER :: JI             ! horizontal loop index
INTEGER :: JK, JKP        ! vertical loop index
INTEGER :: JN             ! chemical tracer loop index
INTEGER :: IMODEIDX       ! mode order index
INTEGER :: JSTEP          ! fractional time loop index
INTEGER :: JKLC, JKLD, JKLP, JKMAX ! loop index for levels
INTEGER :: JKAQ
REAL    :: ZCONV
!
REAL, DIMENSION(KLON,KLEV)     :: ZOMG ! compensat. subsidence (Pa/s)
REAL, DIMENSION(KLON,KLEV,KCH) :: ZUCH1, ZDCH1 ! updraft/downdraft values
REAL, DIMENSION(KLON)          :: ZTIMEC  ! fractional convective time step
REAL, DIMENSION(KLON,KLEV)     :: ZTIMC! 2D work array for ZTIMEC
REAL, DIMENSION(KLON,KLEV,KCH) :: ZCH1MFIN, ZCH1MFOUT
                                   ! work arrays for environm. compensat. mass
REAL, DIMENSION(KLON,KCH)      :: ZWORK1, ZWORK2, ZWORK3
!
! scavenging in updraft
!
REAL, DIMENSION(KLON,KLEV)     :: ZT, ZLWCC, ZLWCI, ZKS1, ZKS2
REAL, DIMENSION(KLON,KLEV)     :: ZKHDUST, ZKHCDUST
REAL, DIMENSION(KLON,KLEV)     :: ZKHSALT, ZKHCSALT
REAL, DIMENSION(KLON,KLEV,KCH) :: ZKH, ZKHC, ZKHI
REAL, DIMENSION(KLON,KLEV,KCH) :: ZFACTOR
REAL, DIMENSION(KLON,KLEV)     :: ZPROS,ZPROR
REAL, DIMENSION(KLON,KLEV,KCH) :: ZPARTSCAV, ZRIM
REAL, DIMENSION(KLON,KLEV)     :: ZPARTAERO
REAL                           :: ZHP
REAL, DIMENSION(KLON,1,KLEV,JPMODE) :: ZRGAER,ZSIGAER, ZNAER, ZBCMIN
REAL, DIMENSION(KLON,1,KLEV,NMODE_DST) :: ZRGDST,ZSIGDST,ZNDST, ZMINMASS_DST,ZRGDSTMIN
REAL, DIMENSION(KLON,1,KLEV,NMODE_SLT) :: ZRGSLT,ZSIGSLT,ZNSLT, ZMINMASS_SLT,ZRGSLTMIN
REAL, DIMENSION(NMODE_DST)    :: ZMIN_DST
REAL, DIMENSION(NMODE_SLT)    :: ZMIN_SLT
REAL, DIMENSION(KLON,1,KLEV,KCH)    :: ZSV, ZSVC
REAL, DIMENSION(KLON,1,KLEV)        :: ZRHODREF
REAL, DIMENSION(NMODE_DST)   :: ZINIRADIUS_DST
REAL, DIMENSION(NMODE_SLT)   :: ZINIRADIUS_SLT
REAL :: ZRHOP, ZFAC, ZMI

LOGICAL, SAVE :: GCHFIRSTCALL = .TRUE.
!
!-------------------------------------------------------------------------------
!
!*       0.1   Compute loop bounds
!              -------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('CH_CONVECT_SCAVENGING',0,ZHOOK_HANDLE)
INCH1  = KCH
IKB    = 1 + JCVEXB
IKS    = KLEV
IKE    = KLEV - JCVEXT
JKMAX  = MAXVAL( KCTL(:) )
!
!
!*       0.2   Initialisation
!              --------------
ZKH(:,:,:) = 0.
ZKHDUST(:,:) = 0.
ZKHSALT(:,:) = 0.
ZKHI(:,:,:) = 0.
ZKHC(:,:,:) = 0.
ZKHCDUST(:,:) = 0.
ZKHCSALT(:,:) = 0.
ZPARTSCAV(:,:,:) = 0.
ZFACTOR(:,:,:) = 0.
ZPROS(:,:) = 0.
ZPROR(:,:) = 0.
ZKS2(:,:) = 0.
ZKS1(:,:) = 0.
!
!*  1.      Fraction of tracer present in the liquid/ice cloud condensate
!           --------------------------------------------------------------
!
!   1.1   Total liquid water content m3/m3
!         --------------------------------
!
ZLWCC(:,:)   = 0.
ZLWCI(:,:)   = 0.
WHERE ((PURC(:,:)+PURR(:,:)).NE.0.) &
       ZLWCC(:,:)   = (PURC(:,:)+PURR(:,:))*PRHODREF(:,:)/XRHOLW
! bulk density for ice varies widely from 0.3 to 0.92 g/cm3
WHERE ((PURI(:,:)+PURS(:,:)).NE.0.) &
       ZLWCI(:,:)   = (PURI(:,:)+PURS(:,:))*PRHODREF(:,:)/910.
!
ZT(:,:)     = PUTT(:,:)
ZHP         = 1.E-5
!
!*    1.2   Fraction of tracer present in the liquid cloud condensate
!           ---------------------------------------------------------
DO JKAQ = NSV_CHEMBEG, NSV_CHEMEND
 WHERE (ZT(:,:).GT.0.)
  ZKS1(:,:) = 1.7E-2 * EXP(-2090.*((1./ZT(:,:))-(1./298.)))
  ZKS2(:,:) = 6.3E-8 * EXP(-1495.*((1./ZT(:,:))-(1./298.)))
  ZKH(:,:,JKAQ) = XSREALHENRYVAL(JKAQ-NSV_CHEMBEG+1,1)* &
         EXP(-XSREALHENRYVAL(JKAQ-NSV_CHEMBEG+1,2)*((1./ZT(:,:))-(1./298.)))
 ENDWHERE
ENDDO
!
DO JKAQ = NSV_CHEMBEG, NSV_CHEMEND
 IF (CNAMES(JKAQ-NSV_CHEMBEG+1)=='N2O5') THEN
    WHERE (ZT(:,:).GT.0.)
       ZKH(:,:,JKAQ)  = 2.00E+0*EXP(3400.*((1./ZT(:,:))-(1./298.)))*1.E20
    END WHERE
  IF (GCHFIRSTCALL) THEN
    WRITE(*,*)'In CH_CONVECT_SCAVENGING: special treatment for',CNAMES(JKAQ-NSV_CHEMBEG+1)
  ENDIF
!
 ELSE IF (CNAMES(JKAQ-NSV_CHEMBEG+1)=='HONO' .OR. CNAMES(JKAQ-NSV_CHEMBEG+1)=='HNO2') THEN
    WHERE (ZT(:,:).GT.0.)
     ZKH(:,:,JKAQ)  = 5.00E+1*EXP(4900.*((1./ZT(:,:))-(1./298.)))* &
                         5.10E-4*EXP(1250.*((1./ZT(:,:))-(1./298.)))/ZHP
    END WHERE
  IF (GCHFIRSTCALL) THEN
    WRITE(*,*)'In CH_CONVECT_SCAVENGING: special treatment for',CNAMES(JKAQ-NSV_CHEMBEG+1)
  ENDIF
 ELSE IF (CNAMES(JKAQ-NSV_CHEMBEG+1)=='HNO3') THEN
    WHERE (ZT(:,:).GT.0.)
      ZKH(:,:,JKAQ)  = 2.10E+5*EXP(8700.*((1./ZT(:,:))-(1./298.)))*15.4/ZHP
    END WHERE
  IF (GCHFIRSTCALL) THEN
    WRITE(*,*)'In CH_CONVECT_SCAVENGING: special treatment for',CNAMES(JKAQ-NSV_CHEMBEG+1)
  ENDIF
 ELSE IF (CNAMES(JKAQ-NSV_CHEMBEG+1)=='SO2') THEN
    WHERE (ZT(:,:).GT.0.)
      ZKH(:,:,JKAQ)   = 1.4*EXP(2900.*((1./ZT(:,:))-(1./298.))) * &
                        ( 1.+(ZKS1(:,:)/ZHP)+(ZKS1(:,:)*ZKS2(:,:)/ZHP/ZHP))
    END WHERE
  IF (GCHFIRSTCALL) THEN
    WRITE(*,*)'In CH_CONVECT_SCAVENGING: special treatment for',CNAMES(JKAQ-NSV_CHEMBEG+1)
  ENDIF
 ELSE IF (CNAMES(JKAQ-NSV_CHEMBEG+1)=='HO2') THEN
    WHERE (ZT(:,:).GT.0.)
      ZKH(:,:,JKAQ)   = 4.00E+3*EXP(5900.*((1./ZT(:,:))-(1./298.)))*2.5E-5/ZHP
    END WHERE
  IF (GCHFIRSTCALL) THEN
    WRITE(*,*)'In CH_CONVECT_SCAVENGING: special treatment for',CNAMES(JKAQ-NSV_CHEMBEG+1)
  ENDIF
 ELSE IF (CNAMES(JKAQ-NSV_CHEMBEG+1)=='ORA2') THEN
    WHERE (ZT(:,:).GT.0.)
      ZKH(:,:,JKAQ)  = 8.80E+3*EXP(6400.*((1./ZT(:,:))-(1./298.)))* &
                          1.70E-5*EXP(50.*((1./ZT(:,:))-(1./298.)))/ZHP
    END WHERE
  IF (GCHFIRSTCALL) THEN
    WRITE(*,*)'In CH_CONVECT_SCAVENGING: special treatment for',CNAMES(JKAQ-NSV_CHEMBEG+1)
  ENDIF
 ENDIF
!
ENDDO
!
IF (LDUST) THEN
  WHERE (ZT(:,:).GT.0.)
    ZKHDUST(:,:)  = 2.10E+5*EXP(8700.*((1./ZT(:,:))-(1./298.)))*15.4/ZHP
  END WHERE
END IF

IF (LSALT) THEN
  WHERE (ZT(:,:).GT.0.)
    ZKHSALT(:,:)  = 2.10E+5*EXP(8700.*((1./ZT(:,:))-(1./298.)))*15.4/ZHP
  END WHERE
END IF
!
GCHFIRSTCALL = .FALSE.
!
!               Convert KH from mol/l/atm in ppp/ppp
!               ------------------------------------
DO JKAQ = NSV_CHEMBEG, NSV_CHEMEND
 ZKHC(:,:,JKAQ) = ZKH(:,:,JKAQ)*0.08205*ZT(:,:)*ZLWCC(:,:)
ENDDO

IF (LDUST) ZKHCDUST(:,:) = ZKHDUST(:,:)*0.08205*ZT(:,:)*ZLWCC(:,:)
IF (LSALT) ZKHCSALT(:,:) = ZKHSALT(:,:)*0.08205*ZT(:,:)*ZLWCC(:,:)

!
!*    1.3   Fraction of tracer present in the ice cloud condensate
!           ---------------------------------------------------------
!
DO JKAQ = NSV_CHEMBEG, NSV_CHEMEND
!
 IF (CNAMES(JKAQ-NSV_CHEMBEG+1)=='HNO3') THEN
    WHERE (ZT(:,:).GT.0.)
      ZKHI(:,:,JKAQ) = 1.E9
    ENDWHERE
 ENDIF
 IF (CNAMES(JKAQ-NSV_CHEMBEG+1)=='H2O2') THEN
! Co-Condensation
    WHERE (ZT(:,:).GT.0.)
      ZKHI(:,:,JKAQ) = (0.2/1.0)*SQRT(18./34.)*PPRESS(:,:)/           &
                     (611.0 * EXP(6138.*((1./273.)-(1./ZT(:,:))))) &
                     * 34./18. / (34./29.) * PURI(:,:)
    ENDWHERE
 ENDIF
ENDDO
!
!
!          Convert KHI from cm3(air)/cm3(ice) in ppp/ppp
!          ---------------------------------------------
DO JKAQ = NSV_CHEMBEG, NSV_CHEMEND
 IF (CNAMES(JKAQ-NSV_CHEMBEG+1)=='HNO3') THEN
  WHERE (ZLWCI(:,:).NE.0.)
    ZKHI(:,:,JKAQ)  = ZKHI(:,:,JKAQ) * ZLWCI(:,:)
  ELSEWHERE
    ZKHI(:,:,JKAQ)  = 0.
  ENDWHERE
 ENDIF
ENDDO
!
!*     1.4  Retention efficiency (R<1 => volatilization during riming
!            in mixed clouds
!           ----------------------------------------------------------
ZRIM(:,:,:)=1.
!
DO JKAQ = NSV_CHEMBEG, NSV_CHEMEND
!
 IF (CNAMES(JKAQ-NSV_CHEMBEG+1)=='OP1' .OR. CNAMES(JKAQ-NSV_CHEMBEG+1)=='CH3OOH') THEN
   WHERE ((PURI(:,:) /=0.).AND.(PURC(:,:) /=0.))
    ZRIM(:,:,JKAQ) =0.02
   ENDWHERE
 ENDIF
 IF (CNAMES(JKAQ-NSV_CHEMBEG+1)=='HCHO' .OR. CNAMES(JKAQ-NSV_CHEMBEG+1)=='CH2O') THEN
   WHERE ((PURI(:,:) /=0.).AND.(PURC(:,:) /=0.))
    ZRIM(:,:,JKAQ)=0.02
   ENDWHERE
 ENDIF
 IF (CNAMES(JKAQ-NSV_CHEMBEG+1)=='H2O2') THEN
   WHERE ((PURI(:,:) /=0.).AND.(PURC(:,:) /=0.))
    ZRIM(:,:,JKAQ)=0.05
   ENDWHERE
 ENDIF
 IF (CNAMES(JKAQ-NSV_CHEMBEG+1)=='HNO3') THEN
   WHERE ((PURI(:,:) /=0.).AND.(PURC(:,:) /=0.))
    ZRIM(:,:,JKAQ)=1.00
   ENDWHERE
 ENDIF
!
ENDDO
!
!*     1.5  Calculate partitioning factor
!           -----------------------------
!
WHERE (PURC(:,:)+PURR(:,:) /= 0.)
  ZPROR(:,:) = PURR(:,:) / ( PURC(:,:)+PURR(:,:) )
ENDWHERE
!
WHERE (PURI(:,:)+PURS(:,:) /= 0.)
  ZPROS(:,:) = PURS(:,:) / ( PURI(:,:)+PURS(:,:) )
ENDWHERE
!
DO JKAQ = NSV_CHEMBEG, NSV_CHEMEND
!
  ZFACTOR(:,:,JKAQ)  = ( 1. + ZKHC(:,:,JKAQ) + ZKHI(:,:,JKAQ) )
  ZPARTSCAV(:,:,JKAQ) =   &
             (ZPROR(:,:) * ZKHC(:,:,JKAQ) * ZRIM(:,:,JKAQ) + &
              ZPROS(:,:) * ZKHI(:,:,JKAQ)) / &
              ZFACTOR(:,:,JKAQ)
  IF (CNAMES(JKAQ-NSV_CHEMBEG+1)=='HNO3') THEN
    ZPARTAERO(:,:) = (ZPROR(:,:) * ZKHC(:,:,JKAQ) ) / ( 1. + ZKHC(:,:,JKAQ) )
  END IF
!
ENDDO
!
IF (LORILAM) THEN
  DO JKAQ = NSV_AERBEG, NSV_AERBEG+(NSP+NSOA+NCARB)*JPMODE-1
    ZPARTSCAV(:,:,JKAQ) = ZPARTAERO(:,:)
  ENDDO
END IF
!
IF (LDUST) THEN
  DO JN=1,NMODE_DST
    IF (LVARSIG) THEN
      ZPARTSCAV(:,:,NSV_DSTBEG-1+2+(JN-1)*3) = &
                           (ZPROR(:,:) * ZKHCDUST(:,:) ) / ( 1. + ZKHCDUST(:,:) )
    ELSE
      ZPARTSCAV(:,:,NSV_DSTBEG-1+2+(JN-1)*2) = &
                           (ZPROR(:,:) * ZKHCDUST(:,:) ) / ( 1. + ZKHCDUST(:,:) )
    END IF
  ENDDO
END IF
!
IF (LSALT) THEN
  DO JN=1,NMODE_SLT
    IF (LVARSIG_SLT) THEN
      ZPARTSCAV(:,:,NSV_SLTBEG-1+2+(JN-1)*3) = &
                           (ZPROR(:,:) * ZKHCSALT(:,:) ) / ( 1. + ZKHCSALT(:,:) )
    ELSE
      ZPARTSCAV(:,:,NSV_SLTBEG-1+2+(JN-1)*2) = &
                           (ZPROR(:,:) * ZKHCSALT(:,:) ) / ( 1. + ZKHCSALT(:,:) )
    END IF
  ENDDO
END IF
!
!
!*      2.      Updraft computations
!               --------------------
!
ZUCH1(:,:,:) = 0.
!
!*      2.1     Initialization  at LCL
!               ----------------------------------
!
DO JI = 1, KLON
    JKLC = KLCL(JI)
    JKLD = KDPL(JI)
    JKLP = KPBL(JI)
    ZWORK1(JI,:) = .5 * ( PCH1(JI,JKLD,:) + PCH1(JI,JKLP,:) )
END DO
!
!*      2.2     Final updraft loop
!               ------------------
!
DO JK = MINVAL( KDPL(:) ), JKMAX
JKP = JK + 1
!
    DO JN = 1, INCH1
     DO JI = 1, KLON
       IF ( KDPL(JI) <= JK .AND. KLCL(JI) > JK )                             &
            ZUCH1(JI,JK,JN) = ZWORK1(JI,JN)
!
       IF ( KLCL(JI) - 1 <= JK .AND. KCTL(JI) > JK ) THEN
           ZUCH1(JI,JKP,JN) = ( PUMF(JI,JK) * ZUCH1(JI,JK,JN) +              &
                                PUER(JI,JKP) * PCH1(JI,JK,JN) )  /           &
           ((1. + ZPARTSCAV(JI,JKP,JN)) * PUMF(JI,JKP) + PUDR(JI,JKP) )
       END IF
     END DO
   END DO
!
END DO
!
!*      3.      Downdraft computations
!               ----------------------
!
ZDCH1(:,:,:) = 0.
!
!*      3.1     Initialization at the LFS
!               -------------------------
!
ZWORK1(:,:) = SPREAD( PMIXF(:), DIM=2, NCOPIES=INCH1 )
DO JI = 1, KLON
     JK = KLFS(JI)
     ZDCH1(JI,JK,:) = ZWORK1(JI,:) * PCH1(JI,JK,:) +                          &
                                       ( 1. - ZWORK1(JI,:) ) * ZUCH1(JI,JK,:)
END DO
!
!*      3.2     Final downdraft loop
!               --------------------
!
DO JK = MAXVAL( KLFS(:) ), IKB + 1, -1
JKP = JK - 1
    DO JN = 1, INCH1
    DO JI = 1, KLON
      IF ( JK <= KLFS(JI) .AND. JKP >= KDBL(JI) ) THEN
       ZDCH1(JI,JKP,JN) = ( ZDCH1(JI,JK,JN) * PDMF(JI,JK) -              &
                            PCH1(JI,JK,JN) *  PDER(JI,JKP) ) /           &
                          ( PDMF(JI,JKP) - PDDR(JI,JKP) - 1.E-7 )
      END IF
    END DO
    END DO
END DO
!
!
!*      4.      Final closure (environmental) computations
!               ------------------------------------------
!
PCH1C(:,:,:) = PCH1(:,:,:) ! initialize adjusted envir. values
!
DO JK = IKB, IKE
   ZOMG(:,JK) = PWSUB(:,JK) * PDXDY(:) / XG ! environmental subsidence
END DO
!
ZTIMEC(:) = PTIMEC(:) / REAL( KFTSTEPS ) ! adjust  fractional time step
                                         ! to be an integer multiple of PTIMEC
WHERE ( PTIMEC(:) < 1. ) ZTIMEC(:) = 0.
ZTIMC(:,:)= SPREAD( ZTIMEC(:), DIM=2, NCOPIES=IKS )
!
ZCH1MFIN(:,:,:)   = 0.
ZCH1MFOUT(:,:,:)  = 0.
!
DO JSTEP = 1, KFTSTEPS ! Enter the fractional time step loop
!
      DO JK = IKB + 1, JKMAX
          JKP = MAX( IKB + 1, JK - 1 )
          ZWORK3(:,:) = SPREAD( ZOMG(:,JK), DIM=2, NCOPIES=INCH1 )
          ZWORK1(:,:) = SIGN( 1., ZWORK3(:,:) )
          ZWORK2(:,:) = 0.5 * ( 1. + ZWORK1(:,:) )
          ZWORK1(:,:) = 0.5 * ( 1. - ZWORK1(:,:) )
          ZCH1MFIN(:,JK,:)  = - ZWORK3(:,:) * PCH1C(:,JKP,:) * ZWORK1(:,:)
          ZCH1MFOUT(:,JK,:) =   ZWORK3(:,:) * PCH1C(:,JK,:)  * ZWORK2(:,:)
         ZCH1MFIN(:,JKP,:) = ZCH1MFIN(:,JKP,:) + ZCH1MFOUT(:,JK,:) * ZWORK2(:,:)
         ZCH1MFOUT(:,JKP,:)= ZCH1MFOUT(:,JKP,:) + ZCH1MFIN(:,JK,:) * ZWORK1(:,:)
      END DO
!
      DO JN = 1, INCH1
       DO JK = IKB + 1, JKMAX
         PCH1C(:,JK,JN) = PCH1C(:,JK,JN) + ZTIMC(:,JK) / PLMASS(:,JK) *  (    &
                      ZCH1MFIN(:,JK,JN) + PUDR(:,JK) * ZUCH1(:,JK,JN) +       &
                      PDDR(:,JK) * ZDCH1(:,JK,JN) - ZCH1MFOUT(:,JK,JN) -      &
                      ( PUER(:,JK) + PDER(:,JK) ) * PCH1(:,JK,JN)    )
         PCH1C(:,JK,JN) = MAX( 0., PCH1C(:,JK,JN) )
       END DO
      END DO
!
END DO ! Exit the fractional time step loop

IF (LORILAM) THEN ! ORILAM chemical aerosol scavenging
ZSV(:,1,:,NSV_AERBEG:NSV_AEREND) = MAX(PCH1(:,:,NSV_AERBEG:NSV_AEREND),1.E-80)
ZSVC(:,1,:,NSV_AERBEG:NSV_AEREND) = MAX(PCH1(:,:,NSV_AERBEG:NSV_AEREND),1.E-80)
DO JN = NSV_AERBEG, NSV_AERBEG+(NSP+NSOA+NCARB)*JPMODE-1
     ZSVC(:,1,IKB:IKE,JN) = MIN(ZSVC(:,IKB:IKE,1,JN),ZSV(:,IKB:IKE,1,JN))
ENDDO

ZRHODREF(:,1,IKB:IKE) = PRHODREF(:,IKB:IKE)

! Compute RG and SIGMA with old concentration PCH1
CALL PPP2AERO(ZSV(:,:,IKB:IKE,NSV_AERBEG:NSV_AEREND),&
              ZRHODREF(:,:,IKB:IKE), PSIG3D=ZSIGAER(:,:,IKB:IKE,:),&
              PRG3D=ZRGAER(:,:,IKB:IKE,:),PN3D=ZNAER(:,:,IKB:IKE,:))

DO JN=NSV_AERBEG,NSV_AEREND,2
  WHERE (ZNAER(:,:,IKB:IKE,1) .LE. XN0IMIN) ! minimum number of particles
                                          ! => no scavenging
    ZSVC(:,:,IKB:IKE,JN) = ZSV(:,:,IKB:IKE,JN)
  END WHERE
ENDDO

DO JN=NSV_AERBEG+1,NSV_AEREND,2
  WHERE (ZNAER(:,:,IKB:IKE,2) .LE. XN0JMIN) ! minimum number of particles
                                          ! => no scavenging
    ZSVC(:,:,IKB:IKE,JN) = ZSV(:,:,IKB:IKE,JN)
  END WHERE
ENDDO
! Limit low values to minimum particles number; minimum mass has been
! put into BC
ZBCMIN(:,:,IKB:IKE,1) = XN0IMIN * (ZRGAER(:,:,IKB:IKE,1)**3)*EXP(4.5 * (LOG(ZSIGAER(:,:,IKB:IKE,1))**2)) *&
                  (4./3.)*3.14292654*1.8e3*1.e-9 * &
                   6.0221367E+11/(1E-6 * XAVOGADRO / XMD*12.*ZRHODREF(:,:,IKB:IKE))

ZBCMIN(:,:,IKB:IKE,2) = XN0JMIN * (ZRGAER(:,:,IKB:IKE,2)**3)*&
                   EXP(4.5 * (LOG(ZSIGAER(:,:,IKB:IKE,2))**2)) *&
                  (4./3.)*3.14292654*1.8e3*1.e-9 * &
                   6.0221367E+11/(1E-6 * XAVOGADRO / XMD*12.*ZRHODREF(:,:,IKB:IKE))

  ZSVC(:,:,IKB:IKE,NSV_AERBEG+JP_CH_BCi-1)  = MAX(ZBCMIN(:,:,IKB:IKE,1), &
                                      ZSVC(:,:,IKB:IKE,NSV_AERBEG+JP_CH_BCi-1))
  ZSVC(:,:,IKB:IKE,NSV_AERBEG+JP_CH_BCj-1)  = MAX(ZBCMIN(:,:,IKB:IKE,2), &
                                      ZSVC(:,:,IKB:IKE,NSV_AERBEG+JP_CH_BCj-1))

END IF
!
! Dust scavenging
IF (LDUST) THEN
  ! Initialization
  ZRHOP = XDENSITY_DUST
  ZFAC  = (4./3.)*XPI*ZRHOP*1.e-9
  ZMI   = XMOLARWEIGHT_DUST ! molecular mass in kg/mol
  ZCONV = ZMI * XM3TOUM3  / (XMD  * ZRHOP * XPI * 4./3.)
  !
  DO JN=1,NMODE_DST
    IMODEIDX = JPDUSTORDER(JN)
    ZMIN_DST(JN) = XN0MIN(IMODEIDX)
    IF (CRGUNITD=="MASS") THEN
      ZINIRADIUS_DST(JN) = XINIRADIUS(IMODEIDX) * EXP(-3.*(LOG(XINISIG(IMODEIDX)))**2)
    ELSE
      ZINIRADIUS_DST(JN) = XINIRADIUS(IMODEIDX)
    END IF
  END DO
  !
  ZRHODREF(:,1,:) = PRHODREF(:,:)
  !
  ZSV(:,1,:,NSV_DSTBEG:NSV_DSTEND)  = PCH1(:,:,NSV_DSTBEG:NSV_DSTEND)
  !
  ZSVC(:,1,:,NSV_DSTBEG:NSV_DSTEND) = ZSV(:,1,:,NSV_DSTBEG:NSV_DSTEND)
  !
  ! Compute RG and SIGMA with old concentration PCH1
  CALL PPP2DUST(ZSV(:,:,IKB:IKE,NSV_DSTBEG:NSV_DSTEND), ZRHODREF(:,:,IKB:IKE),&
                PSIG3D=ZSIGDST(:,:,IKB:IKE,:), PRG3D=ZRGDST(:,:,IKB:IKE,:),   &
                PN3D=ZNDST(:,:,IKB:IKE,:))
  !
  IF (LVARSIG) THEN ! Case with variable standard devisation
  !
    DO JN=1,NMODE_DST
      IMODEIDX = JPDUSTORDER(JN)
      !New mass concentration do not be less than minimum value possible
      ZRGDSTMIN(:,:,:,JN) = ZINIRADIUS_DST(JN)
      ZRGDSTMIN(:,:,IKB:IKE,JN) = MIN(ZRGDSTMIN(:,:,IKB:IKE,JN), ZRGDST(:,:,IKB:IKE,JN))
      ZMINMASS_DST(:,:,IKB:IKE,JN) = XN0MIN(IMODEIDX) * ZRGDSTMIN(:,:,IKB:IKE,JN)**3 * &
                             EXP(4.5 * LOG(ZSIGDST(:,:,IKB:IKE,JN))**2)   * &
                             (ZCONV * ZRHODREF(:,:,IKB:IKE))

      ZSVC(:,1,IKB:IKE,NSV_DSTBEG-1+2+(JN-1)*3) = &
                    PCH1C(:,IKB:IKE,NSV_DSTBEG-1+2+(JN-1)*3)
    ! ZSVC(:,1,IKB:IKE,NSV_DSTBEG-1+2+(JN-1)*3) = &
         !      MIN(PCH1C(:,IKB:IKE,NSV_DSTBEG-1+2+(JN-1)*3),PCH1(:,IKB:IKE,NSV_DSTBEG-1+2+(JN-1)*3))
      WHERE (ZNDST(:,1,IKB:IKE,JN) .GT. ZMIN_DST(JN))
        PCH1C(:,IKB:IKE,NSV_DSTBEG-1+2+(JN-1)*3) = &
                        MAX(ZSVC(:,1,IKB:IKE,NSV_DSTBEG-1+2+(JN-1)*3),&
                        ZMINMASS_DST(:,1,IKB:IKE,JN))
      END WHERE
    ENDDO
    !
  ELSE  ! Case with fixed standard deviation
    DO JN=1,NMODE_DST
      IMODEIDX = JPDUSTORDER(JN)
      !New mass concentration should not be less than minimum value possible
      ZRGDSTMIN(:,:,:,JN) = ZINIRADIUS_DST(JN)
      ZRGDSTMIN(:,:,IKB:IKE,JN) = MIN(ZRGDSTMIN(:,:,IKB:IKE,JN), ZRGDST(:,:,IKB:IKE,JN))
      ZMINMASS_DST(:,:,IKB:IKE,JN) = XN0MIN(IMODEIDX) * ZRGDSTMIN(:,:,IKB:IKE,JN)**3 * &
                             EXP(4.5 * LOG(ZSIGDST(:,:,IKB:IKE,JN))**2)   / &
                             (ZCONV * ZRHODREF(:,:,IKB:IKE))
      ZSVC(:,1,IKB:IKE,NSV_DSTBEG-1+2+(JN-1)*2) = &
                     PCH1C(:,IKB:IKE,NSV_DSTBEG-1+2+(JN-1)*2)
    ! ZSVC(:,1,IKB:IKE,NSV_DSTBEG-1+2+(JN-1)*2) = &
    !          MIN(PCH1C(:,IKB:IKE,NSV_DSTBEG-1+2+(JN-1)*2),PCH1(:,IKB:IKE,NSV_DSTBEG-1+2+(JN-1)*2))
      WHERE (ZNDST(:,1,IKB:IKE,JN) .GT. ZMIN_DST(JN))
         PCH1C(:,IKB:IKE,NSV_DSTBEG-1+2+(JN-1)*2) = MAX(ZSVC(:,1,IKB:IKE,NSV_DSTBEG-1+2+(JN-1)*2),&
                                                  ZMINMASS_DST(:,1,IKB:IKE,JN))
      END WHERE
    ENDDO
  ENDIF
END IF
!
! Sea Salt scavenging
IF (LSALT) THEN
  ! Initialization
  ZRHOP = XDENSITY_SALT
  ZFAC  = (4./3.)*XPI*ZRHOP*1.e-9
  ZMI   = XMOLARWEIGHT_SALT ! molecular mass in kg/mol
  ZCONV = ZMI * XM3TOUM3_SALT  / (XMD  * ZRHOP * XPI * 4./3.)
  !
  DO JN=1,NMODE_SLT
    IMODEIDX = JPSALTORDER(JN)
    ZMIN_SLT(JN) = XN0MIN_SLT(IMODEIDX)
    IF (CRGUNITD=="MASS") THEN
      ZINIRADIUS_SLT(JN) = XINIRADIUS_SLT(IMODEIDX) * EXP(-3.*(LOG(XINISIG_SLT(IMODEIDX)))**2)
    ELSE
      ZINIRADIUS_SLT(JN) = XINIRADIUS_SLT(IMODEIDX)
    END IF
  END DO
  !
  ZRHODREF(:,1,:) = PRHODREF(:,:)
  !
  ZSV(:,1,:,NSV_SLTBEG:NSV_SLTEND)  = PCH1(:,:,NSV_SLTBEG:NSV_SLTEND)
  !
  ZSVC(:,1,:,NSV_SLTBEG:NSV_SLTEND) = ZSV(:,1,:,NSV_SLTBEG:NSV_SLTEND)
  !
  ! Compute RG and SIGMA with old concentration PCH1
  CALL PPP2SALT(ZSV(:,:,IKB:IKE,NSV_SLTBEG:NSV_SLTEND), ZRHODREF(:,:,IKB:IKE),&
                PSIG3D=ZSIGSLT(:,:,IKB:IKE,:), PRG3D=ZRGSLT(:,:,IKB:IKE,:),   &
                PN3D=ZNSLT(:,:,IKB:IKE,:))
  !
  IF (LVARSIG_SLT) THEN ! Case with variable standard devisation
  !
    DO JN=1,NMODE_SLT
      IMODEIDX = JPSALTORDER(JN)
      !New mass concentration do not be less than minimum value possible
      ZRGSLTMIN(:,:,:,JN) = ZINIRADIUS_SLT(JN)
      ZRGSLTMIN(:,:,IKB:IKE,JN) = MIN(ZRGSLTMIN(:,:,IKB:IKE,JN), ZRGSLT(:,:,IKB:IKE,JN))
      ZMINMASS_SLT(:,:,IKB:IKE,JN) = XN0MIN(IMODEIDX) * ZRGSLTMIN(:,:,IKB:IKE,JN)**3 * &
                             EXP(4.5 * LOG(ZSIGSLT(:,:,IKB:IKE,JN))**2)   * &
                             (ZCONV * ZRHODREF(:,:,IKB:IKE))

      ZSVC(:,1,IKB:IKE,NSV_SLTBEG-1+2+(JN-1)*3) = &
                    PCH1C(:,IKB:IKE,NSV_SLTBEG-1+2+(JN-1)*3)
    ! ZSVC(:,1,IKB:IKE,NSV_SLTBEG-1+2+(JN-1)*3) = &
         !      MIN(PCH1C(:,IKB:IKE,NSV_SLTBEG-1+2+(JN-1)*3),PCH1(:,IKB:IKE,NSV_SLTBEG-1+2+(JN-1)*3))
      WHERE (ZNSLT(:,1,IKB:IKE,JN) .GT. ZMIN_SLT(JN))
        PCH1C(:,IKB:IKE,NSV_SLTBEG-1+2+(JN-1)*3) = &
                        MAX(ZSVC(:,1,IKB:IKE,NSV_SLTBEG-1+2+(JN-1)*3),&
                        ZMINMASS_SLT(:,1,IKB:IKE,JN))
      END WHERE
    ENDDO
    !
  ELSE  ! Case with fixed standard deviation
    DO JN=1,NMODE_SLT
      IMODEIDX = JPSALTORDER(JN)
      !New mass concentration should not be less than minimum value possible
      ZRGSLTMIN(:,:,:,JN) = ZINIRADIUS_SLT(JN)
      ZRGSLTMIN(:,:,IKB:IKE,JN) = MIN(ZRGSLTMIN(:,:,IKB:IKE,JN), ZRGSLT(:,:,IKB:IKE,JN))
      ZMINMASS_SLT(:,:,IKB:IKE,JN) = XN0MIN(IMODEIDX) *&
                        ZRGSLTMIN(:,:,IKB:IKE,JN)**3 * &
                        EXP(4.5 * LOG(ZSIGSLT(:,:,IKB:IKE,JN))**2)   / &
                        (ZCONV * ZRHODREF(:,:,IKB:IKE))
      ZSVC(:,1,IKB:IKE,NSV_SLTBEG-1+2+(JN-1)*2) = &
                     PCH1C(:,IKB:IKE,NSV_SLTBEG-1+2+(JN-1)*2)
      WHERE (ZNSLT(:,1,IKB:IKE,JN) .GT. ZMIN_SLT(JN))
         PCH1C(:,IKB:IKE,NSV_SLTBEG-1+2+(JN-1)*2) = &
                MAX(ZSVC(:,1,IKB:IKE,NSV_SLTBEG-1+2+(JN-1)*2),&
                ZMINMASS_SLT(:,1,IKB:IKE,JN))
      END WHERE
    ENDDO
  ENDIF
END IF
!
!
IF (LHOOK) CALL DR_HOOK('CH_CONVECT_SCAVENGING',1,ZHOOK_HANDLE)
END SUBROUTINE CH_CONVECT_SCAVENGING
