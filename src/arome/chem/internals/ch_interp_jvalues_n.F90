!     ######spl
       SUBROUTINE CH_INTERP_JVALUES_n(PJVALUES,PZZ, PZS, PALBUV, PZENITH,&
                                      KLON, KLAT, KLEV,                  &
                                      NIB,NIE,NJB,NJE,NIU,NJU, KVERB, KLUOUT)
       USE PARKIND1, ONLY : JPRB
       USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!      ######################################
! 
!!
!!*** *CH_INTERP_JVALUES_n* interpolate photolysis rates in time and space
!!
!!    PURPOSE
!!    -------
!!      Interpolates J-Values in time and space 
!!
!!**  METHOD
!!       
!!    REFERENCE
!!    ---------
!!    MesoNH documentation
!!
!!    AUTHOR
!!    ------
!!    C. Mari (LA)
!!    
!!    MODIFICATIONS
!!    -------------
!!    Original 21/03/01
!!    P. Tulet 01/11/03  externalisation surface/ UV albedos from
!                        radiations
!!    P. Tulet 01/06/05  updates for arome
!!
!!------------------------------------------------------------------------------
!!
!!    EXTERNAL
!!    --------
!
USE MODD_PARAMETERS
USE MODD_CST
USE MODD_CONF
USE MODD_CH_INIT_JVALUES
!
!!    EXPLICIT ARGUMENTS
!!    ------------------
!
IMPLICIT NONE
!
INTEGER,                  INTENT(IN) :: KLUOUT  
INTEGER,                  INTENT(IN)   :: KLON     ! dimension I
INTEGER,                  INTENT(IN)   :: KLAT     ! dimension J
INTEGER,                  INTENT(IN)   :: KLEV     ! Number of vertical levels
REAL, DIMENSION(KLON,KLAT),       INTENT(IN) :: PZENITH, PALBUV, PZS
REAL, DIMENSION(KLON,KLAT,KLEV),  INTENT(IN) :: PZZ
REAL, DIMENSION(KLON,KLAT,KLEV,JPJVMAX), INTENT(INOUT) :: PJVALUES
INTEGER,                  INTENT(IN)    :: NIB,NIE,NJB,NJE,NIU,NJU   !  domain dim
INTEGER,                   INTENT(IN)    :: KVERB      ! verbosity level
!
!!    LOCAL VARIABLES
!!    ---------------
!
!
INTEGER       :: JI, JJ, JH     ! loop counters
INTEGER       :: IIU, IJU, IKU , IKB, IKE
INTEGER       :: JALB, JJVAL
!
REAL, DIMENSION(:,:,:,:,:), ALLOCATABLE :: ZJDATAALB
!
REAL, DIMENSION(:,:), ALLOCATABLE :: ZSZA
REAL, DIMENSION(:,:), ALLOCATABLE :: ZCOSZEN
                                    ! cosine of zenithal angle
REAL, DIMENSION(:,:), ALLOCATABLE :: ZSINZEN
                                    ! sine of zenithal angle
REAL, DIMENSION(:,:), ALLOCATABLE :: ZAZIMSOL
                                    ! azimuthal solar angle (not used)
REAL, DIMENSION(:,:,:), ALLOCATABLE  :: ZWT
REAL, DIMENSION(:,:),   ALLOCATABLE  :: ZALB1, ZALB2
INTEGER            :: JSZA, JNHT
INTEGER, DIMENSION(:,:), ALLOCATABLE :: JSZAN, JKHTA
REAL,    DIMENSION(:,:), ALLOCATABLE :: ZFSZA, ZOMSZA, ZOMFHTA, ZFHTA
INTEGER :: IIB  ! physical boundary
INTEGER :: IIE  ! physical boundary
INTEGER :: IJB  ! physical boundary
INTEGER :: IJE  ! physical boundary

!
!---------------------------------------------------------------------------
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('CH_INTERP_JVALUES_N',0,ZHOOK_HANDLE)
IF (CPROGRAM == "MOD0D") THEN

IIU=1
IJU=1
IKB=1
IKE=1
IKU=1
IJB=1
IJE=1
IIB=1
IIE=1

ELSE

IIU=NIU
IJU=NJU
IKU=SIZE(PZZ,3)

IIB = NIB 
IIE = NIE
IJB = NJB
IJE = NJE
IKB = 1 + JPVEXT
IKE = IKU - JPVEXT

END IF
!
IF (.NOT.ALLOCATED(ZCOSZEN)) ALLOCATE(ZCOSZEN(IIU,IJU))
IF (.NOT.ALLOCATED(ZSINZEN)) ALLOCATE(ZSINZEN(IIU,IJU))
IF (.NOT.ALLOCATED(ZAZIMSOL)) ALLOCATE(ZAZIMSOL(IIU,IJU))
IF (.NOT.ALLOCATED(ZSZA)) ALLOCATE(ZSZA(IIU,IJU))
IF (.NOT.ALLOCATED(JSZAN)) ALLOCATE(JSZAN(IIU,IJU))
IF (.NOT.ALLOCATED(ZFSZA)) ALLOCATE(ZFSZA(IIU,IJU))
IF (.NOT.ALLOCATED(ZOMSZA)) ALLOCATE(ZOMSZA(IIU,IJU))
IF (.NOT.ALLOCATED(ZOMFHTA)) ALLOCATE(ZOMFHTA(IIU,IJU))
IF (.NOT.ALLOCATED(ZFHTA)) ALLOCATE(ZFHTA(IIU,IJU))
IF (.NOT.ALLOCATED(JKHTA)) ALLOCATE(JKHTA(IIU,IJU))
IF (.NOT.ALLOCATED(ZWT)) ALLOCATE(ZWT(IIU,IJU,4))
IF (.NOT.ALLOCATED(ZALB1)) ALLOCATE(ZALB1(IIU,IJU))
IF (.NOT.ALLOCATED(ZALB2)) ALLOCATE(ZALB2(IIU,IJU))
IF (.NOT.ALLOCATED(ZJDATAALB)) &
                 ALLOCATE(ZJDATAALB(IIU,IJU,IKU,JPJVMAX,NBALB))

 ZJDATAALB(:,:,:,:,:)=0.
!
ZCOSZEN(:,:) = COS(PZENITH(:,:))
ZSZA(:,:) = ACOS(ZCOSZEN(:,:)) * 180. / (2.*ASIN(1.))
!
   !* find index for solar zenith angle
   !  ---------------------------------
   !
JSZAN(:,:)  = 2
DO JSZA = 2, NSZA_INCR-1
  WHERE(ZSZA(:,:).GT.XSZA_JVAL(JSZA))
    JSZAN(:,:) = JSZA+1
  END WHERE
END DO

ZFSZA(:,:)  = 0.
DO JJ=IJB,IJE
  DO JI=IIB,IIE
    ZFSZA(JI,JJ)  = (XSZA_JVAL(JSZAN(JI,JJ)) - ZSZA(JI,JJ))/ &
         (XSZA_JVAL(JSZAN(JI,JJ)) - XSZA_JVAL(JSZAN(JI,JJ)-1))
  END DO
END DO
ZOMSZA(:,:) = 1. - ZFSZA(:,:)

   !
   ! * Vertical interpolation of J values
   !   ----------------------------------
   !
!
ZFHTA(:,:) = 0.
DO JH = 1, IKE-IKB+1
  !
  JKHTA(:,:) = 2
  DO JNHT = 2, NZZ_JVAL - 1
    WHERE((PZZ(:,:,JH)-PZS(:,:)).GT.XZZ_JVAL(JNHT)) JKHTA(:,:) = JNHT + 1
  END DO
  DO JJ=IJB,IJE
    DO JI=IIB,IIE
      ZFHTA(JI,JJ) = (XZZ_JVAL(JKHTA(JI,JJ)) - (PZZ(JI,JJ,JH)-PZS(JI,JJ))) &
           /(XZZ_JVAL(JKHTA(JI,JJ))-XZZ_JVAL(JKHTA(JI,JJ)-1))
    END DO
  END DO
  ZOMFHTA(:,:) = 1. - ZFHTA(:,:)
  !
  ZWT(:,:,1) = ZOMSZA(:,:) * ZFHTA(:,:)
  ZWT(:,:,2) = ZFSZA(:,:)  * ZFHTA(:,:)
  ZWT(:,:,3) = ZFSZA(:,:)  * ZOMFHTA(:,:)
  ZWT(:,:,4) = ZOMSZA(:,:) * ZOMFHTA(:,:)
  !

  !
  DO JALB = 1, NBALB
    DO JJVAL = 1, JPJVMAX
      DO JJ=IJB,IJE
        DO JI=IIB,IIE
          ZJDATAALB(JI,JJ,JH,JJVAL,JALB) =                              &
           ZWT(JI,JJ,1) * XJDATA(JSZAN(JI,JJ)-1,JKHTA(JI,JJ)-1,JJVAL,JALB)  &
         + ZWT(JI,JJ,2) * XJDATA(JSZAN(JI,JJ),JKHTA(JI,JJ)-1,JJVAL,JALB)    &
         + ZWT(JI,JJ,3) * XJDATA(JSZAN(JI,JJ)-1,JKHTA(JI,JJ),JJVAL,JALB)    &
         + ZWT(JI,JJ,4) * XJDATA(JSZAN(JI,JJ),JKHTA(JI,JJ),JJVAL,JALB)
        END DO
      END DO
    ENDDO
  ENDDO
ENDDO

ZJDATAALB(:,:,IKE,:,:) = ZJDATAALB(:,:,IKE-IKB+1,:,:)
ZJDATAALB(:,:,IKU,:,:) = ZJDATAALB(:,:,IKE-IKB+1,:,:)


!
DO JALB = 1, NBALB
 DO JJVAL = 1, JPJVMAX
  DO JH = IKB, IKE
      WHERE(ZSZA(:,:).GT.100.) ZJDATAALB(:,:,JH,JJVAL,JALB)=0.
  ENDDO
 ENDDO
ENDDO
!
!
PJVALUES(:,:,:,:) =  ZJDATAALB(:,:,:,:,1)
!

DO JALB=1,NBALB-1
  ZALB1(:,:) = 0.02+0.20*FLOAT(JALB-1)/FLOAT(NBALB-1)
  ZALB2(:,:) = 0.02+0.20*FLOAT(JALB  )/FLOAT(NBALB-1)
  
  DO JJVAL = 1, JPJVMAX
    DO JH = IKB, IKE
      WHERE(PALBUV(:,:)>ZALB1(:,:) .AND. PALBUV(:,:)<=ZALB2(:,:))
        PJVALUES(:,:,JH,JJVAL) = &
             ( PALBUV(:,:) - ZALB1(:,:) ) / ( ZALB2(:,:) - ZALB1(:,:) ) * ZJDATAALB(:,:,JH,JJVAL,JALB  ) &
             + ( ZALB2(:,:) - PALBUV(:,:) ) / ( ZALB2(:,:) - ZALB1(:,:)) * ZJDATAALB(:,:,JH,JJVAL,JALB+1)
      END WHERE
      WHERE(PALBUV(:,:) > ZALB2(:,:)) PJVALUES(:,:,JH,JJVAL) = ZJDATAALB(:,:,JH,JJVAL,JALB+1)
    END DO
  END DO
END DO
!
!
PJVALUES(:,:,:,:) = MAX(0.,PJVALUES(:,:,:,:))
!
IF (KVERB >= 6) THEN
WRITE(KLUOUT,*) "photolysis rates interpolated in time and space" 
END IF

DEALLOCATE(ZCOSZEN)
DEALLOCATE(ZSINZEN)
DEALLOCATE(ZAZIMSOL)
DEALLOCATE(ZSZA)
DEALLOCATE(ZJDATAALB)
DEALLOCATE(JSZAN)
DEALLOCATE(ZFSZA)
DEALLOCATE(ZOMSZA)
DEALLOCATE(ZOMFHTA)
DEALLOCATE(ZFHTA)
DEALLOCATE(ZWT)
DEALLOCATE(ZALB1)
DEALLOCATE(ZALB2)

!
IF (LHOOK) CALL DR_HOOK('CH_INTERP_JVALUES_N',1,ZHOOK_HANDLE)
END SUBROUTINE CH_INTERP_JVALUES_n
