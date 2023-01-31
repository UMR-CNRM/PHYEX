!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ##########################
      MODULE MODI_LIMA_NOTADJUST
!     ##########################
!
INTERFACE
!
      SUBROUTINE LIMA_NOTADJUST(KMI, TPFILE, HRAD,                                       &
                                PTSTEP, PRHODJ, PPABSTT,  PPABST, PRHODREF, PEXNREF, PZZ, &
                                PTHT,PRT, PSVT, PTHS, PRS,PSVS, PCLDFR, PICEFR, PRAINFR, PSRCS )
!
USE MODD_IO, ONLY: TFILEDATA
USE MODD_NSV,   only: NSV_LIMA_BEG
!
INTEGER,                  INTENT(IN)    :: KMI        ! Model index
TYPE(TFILEDATA),          INTENT(IN)    :: TPFILE   ! Output file
CHARACTER(len=4),         INTENT(IN)    :: HRAD     ! Radiation scheme name
REAL,                     INTENT(IN)    :: PTSTEP   ! Double Time step
                                                    ! (single if cold start)
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PPABSTT  ! Absolute Pressure at t+dt
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PPABST  ! Absolute Pressure at t     
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF ! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PEXNREF ! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PZZ      ! Reference density
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT     ! Theta
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRT      ! m.r. at t
REAL, DIMENSION(:,:,:,NSV_LIMA_BEG:), INTENT(IN)    :: PSVT     ! Concentrations at t
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PTHS     ! Theta source
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRS      ! m.r. source
REAL, DIMENSION(:,:,:,NSV_LIMA_BEG:), INTENT(INOUT) :: PSVS     ! Concentrations source
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PSRCS   ! Second-order flux
                                                   ! s'rc'/2Sigma_s2 at time t+1
                                                   ! multiplied by Lambda_3
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PCLDFR  ! Cloud fraction
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PICEFR    ! Cloud fraction          
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRAINFR   ! Cloud fraction          
!
!
END SUBROUTINE LIMA_NOTADJUST
!
END INTERFACE
!
END MODULE MODI_LIMA_NOTADJUST
!
!     ####################################################################################
      SUBROUTINE LIMA_NOTADJUST(KMI, TPFILE, HRAD,                                       &
                                PTSTEP, PRHODJ, PPABSTT,  PPABST, PRHODREF, PEXNREF, PZZ, &
                                PTHT,PRT, PSVT, PTHS, PRS,PSVS, PCLDFR, PICEFR, PRAINFR, PSRCS )
!     ####################################################################################
!
!!****  * -  compute pseudo-prognostic of supersaturation according to Thouron
!                                                                     et al. 2012
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!
!!    REFERENCE
!!    ---------
!!
!!      Thouron, O., J.-L. Brenguier, and F. Burnet, Supersaturation calculation
!!      in large eddy simulation models for prediction of the droplet number
!!      concentration, Geosci. Model Dev., 5, 761-772, 2012.
!!
!!    AUTHOR
!!    ------
!!      B.Vie forked from lima_adjust.f90
!!
!!    MODIFICATIONS
!!    -------------
!
!*       0.    DECLARATIONS
!
use modd_budget,           only: lbu_enable, nbumod,                                          &
                                 lbudget_th, lbudget_rv, lbudget_rc, lbudget_ri, lbudget_sv,  &
                                 NBUDGET_TH, NBUDGET_RV, NBUDGET_RC, NBUDGET_RI, NBUDGET_SV1, &
                                 tbudgets
USE MODD_CONF
USE MODD_CST
USE MODD_FIELD,            ONLY: TFIELDMETADATA, TYPEREAL
USE MODD_IO,               ONLY: TFILEDATA
USE MODD_LUNIT_n,          ONLY: TLUOUT
USE MODD_NSV
USE MODD_PARAMETERS
USE MODD_PARAM_LIMA
USE MODD_PARAM_LIMA_COLD

!
use mode_budget,           only: Budget_store_init, Budget_store_end
USE MODE_IO_FIELD_WRITE,   only: IO_Field_write
USE MODE_MSG
use mode_tools,            only: Countjv
use mode_tools_ll,         only: GET_INDICE_ll
!
USE MODI_PROGNOS_LIMA
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
INTEGER,                  INTENT(IN)    :: KMI        ! Model index
TYPE(TFILEDATA),          INTENT(IN)    :: TPFILE   ! Output file
CHARACTER(len=4),              INTENT(IN)    :: HRAD     ! Radiation scheme name
REAL,                     INTENT(IN)    :: PTSTEP   ! Double Time step
                                                    ! (single if cold start)
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PPABSTT  ! Absolute Pressure at t+dt
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PPABST  ! Absolute Pressure at t     
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF ! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PEXNREF ! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PZZ      ! Reference density
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT     ! Theta
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRT      ! m.r. at t
REAL, DIMENSION(:,:,:,NSV_LIMA_BEG:), INTENT(IN)    :: PSVT     ! Concentrations at t
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PTHS     ! Theta source
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRS      ! m.r. source
REAL, DIMENSION(:,:,:,NSV_LIMA_BEG:), INTENT(INOUT) :: PSVS     ! Concentrations source
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PSRCS   ! Second-order flux
                                                   ! s'rc'/2Sigma_s2 at time t+1
                                                   ! multiplied by Lambda_3
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PCLDFR  ! Cloud fraction
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PICEFR    ! Cloud fraction          
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRAINFR   ! Cloud fraction          
!
!
!*       0.2   Declarations of local variables :
!
!
!
INTEGER             :: IRESP      ! Return code of FM routines
INTEGER             :: ILUOUT     ! Logical unit of output listing 

! For Activation :                       
LOGICAL, DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3)) &
            :: GNUCT, GMICRO  ! Test where to compute the HEN process
INTEGER , DIMENSION(SIZE(GNUCT))  :: I1,I2,I3 ! Used to replace the COUNT
INTEGER                           :: JL, JMOD       ! and PACK intrinsics
REAL, DIMENSION(:), ALLOCATABLE ::ZPRES,ZRHOD,ZRR,ZTT,ZRV,ZRC,ZS0,ZCCL, &
                                  ZZDZ, ZZLV, ZZCPH,                    &
                                  ZRVT, ZRIT, ZCIT, ZRVS, ZRIS, ZCIS,   &
                                  ZTHS, ZRHODREF, ZZT, ZEXNREF, ZZW,    &
                                  ZLSFACT, ZRVSATI, ZRVSATI_PRIME,      &
                                  ZDELTI, ZAI, ZKA, ZDV, ZITI, ZAII, ZDEP, &
                                  ZCJ 
!
INTEGER :: INUCT
INTEGER :: IMICRO
INTEGER :: IIB           !  Define the domain where 
INTEGER :: IIE           !  the microphysical sources have to be computed
INTEGER :: IJB           ! 
INTEGER :: IJE           !
INTEGER :: IKB           ! 
INTEGER :: IKE           !

REAL, DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3)) ::&
                       ZEXNT,ZEXNS,ZT,ZRVSAT,ZWORK,ZLV,ZLS,ZCPH, ZW1,        &
                       ZDZ, ZW
REAL, DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3)) ::&
                       ZSAT,ZCCS
INTEGER           :: JK            ! For loop
integer :: idx
TYPE(TFIELDMETADATA) :: TZFIELD
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZNFS     ! CCN C. available source
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZNAS     ! Cloud  C. nuclei C. source
REAL, DIMENSION(:,:), ALLOCATABLE :: ZZNFS     ! CCN C. available source
REAL, DIMENSION(:,:), ALLOCATABLE :: ZZNAS     ! Cloud  C. nuclei C. source
REAL :: ZEPS

!-------------------------------------------------------------------------------
!
!*       1.     PRELIMINARIES
!               -------------
!
ILUOUT = TLUOUT%NLU
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
IKB=1+JPVEXT
IKE=SIZE(PZZ,3) - JPVEXT
!
!-------------------------------------------------------------------------------
if ( nbumod == kmi .and. lbu_enable ) then
  if ( lbudget_th ) call Budget_store_init( tbudgets(NBUDGET_TH), 'CEDS', pths(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rv ) call Budget_store_init( tbudgets(NBUDGET_RV), 'CEDS', prs(:, :, :, 1) * prhodj(:, :, :) )
  if ( lbudget_rc ) call Budget_store_init( tbudgets(NBUDGET_RC), 'CEDS', prs(:, :, :, 2) * prhodj(:, :, :) )
  if ( lbudget_ri ) call Budget_store_init( tbudgets(NBUDGET_RI), 'CEDS', prs(:, :, :, 4) * prhodj(:, :, :) )
  if ( lbudget_sv ) then
    if ( nmom_c.ge.2) then
      call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nc), 'CEDS', psvs(:, :, :, nsv_lima_nc) * prhodj(:, :, :) )
      do jl = nsv_lima_ccn_free, nsv_lima_ccn_free + nmod_ccn - 1
        idx = NBUDGET_SV1 - 1 + jl
        call Budget_store_init( tbudgets(idx), 'CEDS', psvs(:, :, :, jl) * prhodj(:, :, :) )
      end do
      do jl = nsv_lima_ccn_acti, nsv_lima_ccn_acti + nmod_ccn - 1
        idx = NBUDGET_SV1 - 1 + jl
        call Budget_store_init( tbudgets(idx), 'CEDS', psvs(:, :, :, jl) * prhodj(:, :, :) )
      end do
    end if
!     if ( lscav .and. laero_mass ) &
!       call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_scavmass), 'CEDS', psvs(:, :, :, nsv_lima_scavmass) &
!                                                                                     * prhodj(:, :, :) )
!     if ( nmom_i.ge.2 ) then
!       call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ni), 'CEDS', psvs(:, :, :, nsv_lima_ni) * prhodj(:, :, :) )
!       do jl = 1, nsv_lima_ifn_free, nsv_lima_ifn_free + nmod_ifn - 1
!         idx = NBUDGET_SV1 - 1 + jl
!         call Budget_store_init( tbudgets(idx), 'CEDS', psvs(:, :, :, jl) * prhodj(:, :, :) )
!       end do
!       do jl = 1, nsv_lima_ifn_nucl, nsv_lima_ifn_nucl + nmod_ifn - 1
!         idx = NBUDGET_SV1 - 1 + jl
!         call Budget_store_init( tbudgets(idx), 'CEDS', psvs(:, :, :, jl) * prhodj(:, :, :) )
!       end do
!       do jl = 1, nsv_lima_imm_nucl, nsv_lima_imm_nucl + nmod_ifn - 1
!         idx = NBUDGET_SV1 - 1 + jl
!         call Budget_store_init( tbudgets(idx), 'CEDS', psvs(:, :, :, jl) * prhodj(:, :, :) )
!       end do
!     end if
      call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_spro), 'CEDS', psvs(:, :, :, nsv_lima_spro) * prhodj(:, :, :) )
  end if
end if
!
!*       2.     COMPUTE QUANTITIES WITH THE GUESS OF THE FUTURE INSTANT
!               -------------------------------------------------------
!
!*       2.1    remove negative non-precipitating negative water
!               ------------------------------------------------
!
IF (ANY(PRS(:,:,:,2) < 0. .OR. PSVS(:,:,:,NSV_LIMA_NC) < 0.)) THEN
  WRITE(ILUOUT,*) 'LIMA_NOTADJUST beginning:  negative values of PRCS or PCCS'
  WRITE(ILUOUT,*) '  location of minimum of PRCS:', MINLOC(PRS(:,:,:,2))
  WRITE(ILUOUT,*) ' value of minimum   :', MINVAL(PRS(:,:,:,2))
  WRITE(ILUOUT,*) '  location of minimum of PCCS:', MINLOC(PSVS(:,:,:,NSV_LIMA_NC))
  WRITE(ILUOUT,*) ' value of minimum   :', MINVAL(PSVS(:,:,:,NSV_LIMA_NC))
END IF
!
IF (ANY(PRS(:,:,:,2)+PRS(:,:,:,1) < 0.) .AND. NVERB>5) THEN
  WRITE(ILUOUT,*) 'LIMA_NOT_ADJUST:  negative values of total water (reset to zero)'
  WRITE(ILUOUT,*) '  location of minimum:', MINLOC(PRS(:,:,:,2)+PRS(:,:,:,1))
  WRITE(ILUOUT,*) '  value of minimum   :', MINVAL(PRS(:,:,:,2)+PRS(:,:,:,1))
!callabortstop
  CALL PRINT_MSG(NVERB_FATAL,'GEN','LIMA_NOTADJUST','')
END IF
!
!
!*       2.2    estimate the Exner function at t+1 and t respectively
!
ZEXNS(:,:,:)=(PPABSTT(:,:,:)/XP00 )**(XRD/XCPD)
ZEXNT(:,:,:)=(PPABST(:,:,:)/XP00 )**(XRD/XCPD)
!sources terms *dt
PRS(:,:,:,:)   = PRS(:,:,:,:) * PTSTEP
PSVS(:,:,:,:)  = PSVS(:,:,:,:) * PTSTEP
ZSAT(:,:,:) = PSVS(:,:,:,NSV_LIMA_SPRO)-1.0
ZCCS(:,:,:) = PSVS(:,:,:,NSV_LIMA_NC)
IF ( NMOD_CCN .GE. 1 ) THEN
   ALLOCATE( ZNFS(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3),NMOD_CCN) )
   ALLOCATE( ZNAS(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3),NMOD_CCN) )
   ZNFS(:,:,:,:) = PSVS(:,:,:,NSV_LIMA_CCN_FREE:NSV_LIMA_CCN_FREE+NMOD_CCN-1)
   ZNAS(:,:,:,:) = PSVS(:,:,:,NSV_LIMA_CCN_ACTI:NSV_LIMA_CCN_ACTI+NMOD_CCN-1)
ELSE
   ALLOCATE( ZNFS(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3),1) )
   ALLOCATE( ZNAS(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3),1) )
   ZNFS(:,:,:,:) = 0.
   ZNAS(:,:,:,:) = 0.
END IF
ZW(:,:,:)=SUM(ZNAS,4)
!
!state temperature at t+dt
PTHS(:,:,:)   = PTHS(:,:,:) * PTSTEP * ZEXNS(:,:,:)

!state temperature at t
ZT(:,:,:)=PTHT(:,:,:)*ZEXNT(:,:,:)
!Lv and Cph at t
ZLV(:,:,:) = XLVTT+(XCPV-XCL)*(ZT(:,:,:)-XTT)                
ZLS(:,:,:) = XLSTT + ( XCPV - XCI ) * ( ZT(:,:,:) -XTT )
ZCPH(:,:,:)= XCPD+XCPV*PRT(:,:,:,1)+XCL*(PRT(:,:,:,2)+PRT(:,:,:,3)) &
                                   +XCI*(PRT(:,:,:,4)+PRT(:,:,:,5)+PRT(:,:,:,6))
!dz
DO JK=1,IKE                 
  ZDZ(:,:,JK)=PZZ(:,:,JK+1)-PZZ(:,:,JK)
END DO
!
!*       2.3    compute the latent heat of vaporization Lv(T*) at t+1
!
!Removed negligible values
!
WHERE ( ((PRS(:,:,:,2).LT.XRTMIN(2)) .AND. (ZSAT(:,:,:).LT.0.0)) .OR. &
        ((PRS(:,:,:,2).GT.0.0)       .AND. (ZCCS(:,:,:).LE.0.0)) )
 PTHS(:,:,:)  = PTHS(:,:,:)-(ZLV(:,:,:)/ZCPH(:,:,:))*PRS(:,:,:,2)
 PRS(:,:,:,1)  = PRS(:,:,:,1)+PRS(:,:,:,2)
 PRS(:,:,:,2)  = 0.0
!ZSAT(:,:,:)   = 0.0
 ZCCS(:,:,:)   = 0.0
!ZNFS(:,:,:,1:NMOD_CCN) = ZNFS(:,:,:,1:NMOD_CCN) + ZNAS(:,:,:,1:NMOD_CCN)
!ZNAS(:,:,:,1:NMOD_CCN) = 0.
END WHERE
!


!
! Ice deposition/sublimation
! 
ZEPS= XMV / XMD
GMICRO(:,:,:)=.FALSE.
GMICRO(IIB:IIE,IJB:IJE,IKB:IKE) = (PRS(IIB:IIE,IJB:IJE,IKB:IKE,4)>XRTMIN(4)/PTSTEP .AND.      &
                                   PSVS(IIB:IIE,IJB:IJE,IKB:IKE,NSV_LIMA_NI)>XCTMIN(4)/PTSTEP )
IMICRO = COUNTJV( GMICRO(:,:,:),I1(:),I2(:),I3(:))
IF( IMICRO >= 1 .AND. .NOT.LPTSPLIT) THEN
   ALLOCATE(ZRVT(IMICRO))
   ALLOCATE(ZRIT(IMICRO))
   ALLOCATE(ZCIT(IMICRO))
!
   ALLOCATE(ZRVS(IMICRO))
   ALLOCATE(ZRIS(IMICRO))
   ALLOCATE(ZCIS(IMICRO))  !!!BVIE!!!
   ALLOCATE(ZTHS(IMICRO))
!
   ALLOCATE(ZRHODREF(IMICRO))
   ALLOCATE(ZZT(IMICRO))
   ALLOCATE(ZPRES(IMICRO))
   ALLOCATE(ZEXNREF(IMICRO))
   ALLOCATE(ZZCPH(IMICRO))
   DO JL=1,IMICRO
      ZRVT(JL) = PRT(I1(JL),I2(JL),I3(JL),1)
      ZRIT(JL) = PRT(I1(JL),I2(JL),I3(JL),4)
      ZCIT(JL) = PSVT(I1(JL),I2(JL),I3(JL),NSV_LIMA_NI)
!
      ZRVS(JL) = PRS(I1(JL),I2(JL),I3(JL),1)
      ZRIS(JL) = PRS(I1(JL),I2(JL),I3(JL),4)
      ZCIS(JL) = PSVS(I1(JL),I2(JL),I3(JL),NSV_LIMA_NI) !!!BVIE!!!
      ZTHS(JL) = PTHS(I1(JL),I2(JL),I3(JL))
!
      ZRHODREF(JL) = PRHODREF(I1(JL),I2(JL),I3(JL))
      ZZT(JL) = ZT(I1(JL),I2(JL),I3(JL))
      ZPRES(JL) = PPABSTT(I1(JL),I2(JL),I3(JL))
      ZEXNREF(JL) = PEXNREF(I1(JL),I2(JL),I3(JL))
      ZZCPH(JL) = ZCPH(I1(JL),I2(JL),I3(JL))
   ENDDO
   ALLOCATE(ZZW(IMICRO))
   ALLOCATE(ZLSFACT(IMICRO))
   ZLSFACT(:) = (XLSTT+(XCPV-XCI)*(ZZT(:)-XTT))/ZZCPH(:) ! L_s/C_ph
   ALLOCATE(ZRVSATI(IMICRO))
   ALLOCATE(ZRVSATI_PRIME(IMICRO))
   ALLOCATE(ZDELTI(IMICRO))
   ALLOCATE(ZAI(IMICRO))
   ALLOCATE(ZCJ(IMICRO))
   ALLOCATE(ZKA(IMICRO))
   ALLOCATE(ZDV(IMICRO))
   ALLOCATE(ZITI(IMICRO))
!
   ZKA(:) = 2.38E-2 + 0.0071E-2 * ( ZZT(:) - XTT )          ! k_a
   ZDV(:) = 0.211E-4 * (ZZT(:)/XTT)**1.94 * (XP00/ZPRES(:)) ! D_v
   ZCJ(:) = XSCFAC * ZRHODREF(:)**0.3 / SQRT( 1.718E-5+0.0049E-5*(ZZT(:)-XTT) )
!
   ZZW(:) = EXP( XALPI - XBETAI/ZZT(:) - XGAMI*ALOG(ZZT(:) ) ) ! es_i
   ZRVSATI(:) = ZEPS*ZZW(:) / ( ZPRES(:)-ZZW(:) )              ! r_si
   ZRVSATI_PRIME(:) = (( XBETAI/ZZT(:) - XGAMI ) / ZZT(:))  &  ! r'_si
                       * ZRVSATI(:) * ( 1. + ZRVSATI(:)/ZEPS )
!
   ZDELTI(:) = ZRVS(:)*PTSTEP - ZRVSATI(:)
   ZAI(:) = ( XLSTT + (XCPV-XCI)*(ZZT(:)-XTT) )**2 / (ZKA(:)*XRV*ZZT(:)**2) &
                                  + ( XRV*ZZT(:) ) / (ZDV(:)*ZZW(:))
   ZZW(:) = MIN(1.E8,( XLBI* MAX(ZCIT(:),XCTMIN(4))                       &
                           /(MAX(ZRIT(:),XRTMIN(4))) )**XLBEXI)
                                                                  ! Lbda_I
   ZITI(:) = ZCIT(:) * (X0DEPI/ZZW(:) + X2DEPI*ZCJ(:)*ZCJ(:)/ZZW(:)**(XDI+2.0)) &
                     / (ZRVSATI(:)*ZAI(:))
!
   ALLOCATE(ZAII(IMICRO))
   ALLOCATE(ZDEP(IMICRO))
!
   ZAII(:) = 1.0 + ZRVSATI_PRIME(:)*ZLSFACT(:)
   ZDEP(:) = 0.0
!
   ZZW(:)  = ZAII(:)*ZITI(:)*PTSTEP ! R*delta_T
   WHERE( ZZW(:)<1.0E-2 )
      ZDEP(:) = ZITI(:)*ZDELTI(:)*(1.0 - (ZZW(:)/2.0)*(1.0-ZZW(:)/3.0))
   ELSEWHERE          
      ZDEP(:) = ZITI(:)*ZDELTI(:)*(1.0 - EXP(-ZZW(:)))/ZZW(:)
   END WHERE
!
! Integration
!
   WHERE( ZDEP(:) < 0.0 )
      ZDEP(:) = MAX ( ZDEP(:), -ZRIS(:) )
   ELSEWHERE
      ZDEP(:) = MIN ( ZDEP(:),  ZRVS(:) )
!      ZDEP(:) = MIN ( ZDEP(:),  ZCIS(:)*5.E-10 ) !!!BVIE!!!
   END WHERE
   WHERE( ZRIS(:) < XRTMIN(4)/PTSTEP )
      ZDEP(:) = 0.0
   END WHERE
   ZRVS(:) = ZRVS(:) - ZDEP(:)
   ZRIS(:) = ZRIS(:) + ZDEP(:)
   ZTHS(:) = ZTHS(:) + ZDEP(:) * ZLSFACT(:) / ZEXNREF(:)
!
!  Implicit ice crystal sublimation if ice saturated conditions are not met
!
   ZZT(:) = ( ZTHS(:) * PTSTEP ) * ( ZPRES(:) / XP00 ) ** (XRD/XCPD)
   ZZW(:) = EXP( XALPI - XBETAI/ZZT(:) - XGAMI*ALOG(ZZT(:) ) ) ! es_i
   ZRVSATI(:) = ZEPS*ZZW(:) / ( ZPRES(:)-ZZW(:) )              ! r_si
   WHERE( ZRVS(:)*PTSTEP<ZRVSATI(:) )
      ZZW(:)  = ZRVS(:) + ZRIS(:)
      ZRVS(:) = MIN( ZZW(:),ZRVSATI(:)/PTSTEP )
      ZTHS(:) = ZTHS(:) + ( MAX( 0.0,ZZW(:)-ZRVS(:) )-ZRIS(:) ) &
                          * ZLSFACT(:) / ZEXNREF(:)
      ZRIS(:) = MAX( 0.0,ZZW(:)-ZRVS(:) )
   END WHERE
!
!
   ZW(:,:,:) = PRS(:,:,:,1)
   PRS(:,:,:,1) = UNPACK( ZRVS(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
   ZW(:,:,:) = PRS(:,:,:,4)
   PRS(:,:,:,4) = UNPACK( ZRIS(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
   ZW(:,:,:) = PTHS(:,:,:)
   PTHS(:,:,:) = UNPACK( ZTHS(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
!
   DEALLOCATE(ZRVT)
   DEALLOCATE(ZRIT)
   DEALLOCATE(ZCIT)
   DEALLOCATE(ZRVS)
   DEALLOCATE(ZRIS)
   DEALLOCATE(ZCIS) !!!BVIE!!!
   DEALLOCATE(ZTHS)
   DEALLOCATE(ZRHODREF)
   DEALLOCATE(ZZT)
   DEALLOCATE(ZPRES)
   DEALLOCATE(ZEXNREF)
   DEALLOCATE(ZZCPH)
   DEALLOCATE(ZZW)
   DEALLOCATE(ZLSFACT)
   DEALLOCATE(ZRVSATI)
   DEALLOCATE(ZRVSATI_PRIME)
   DEALLOCATE(ZDELTI)
   DEALLOCATE(ZAI)
   DEALLOCATE(ZCJ)
   DEALLOCATE(ZKA)
   DEALLOCATE(ZDV)
   DEALLOCATE(ZITI)
   DEALLOCATE(ZAII)
   DEALLOCATE(ZDEP)
END IF ! IMICRO
!
!selection of mesh where condensation/evaportion/activation is performed
GNUCT(:,:,:) = .FALSE.
!GNUCT(IIB:IIE,IJB:IJE,IKB:IKE) = ZSAT(IIB:IIE,IJB:IJE,IKB:IKE)>0.0 .OR.   &
!                                 ZCCS(IIB:IIE,IJB:IJE,IKB:IKE)>0.0
!GNUCT(IIB:IIE,IJB:IJE,IKB:IKE) = ZSAT(IIB:IIE,IJB:IJE,IKB:IKE)>0.0 
GNUCT(IIB:IIE,IJB:IJE,IKB:IKE) = ZSAT(IIB:IIE,IJB:IJE,IKB:IKE)>-1.0 .AND.  &
                                 ( ZSAT(IIB:IIE,IJB:IJE,IKB:IKE)>0.0 .OR.  &
!                                ZCCS(IIB:IIE,IJB:IJE,IKB:IKE)>1.E+05
                                 ZCCS(IIB:IIE,IJB:IJE,IKB:IKE)>XCTMIN(2) )
INUCT = COUNTJV( GNUCT(:,:,:),I1(:),I2(:),I3(:))
!3D array to 1D array
!
IF( INUCT >= 1 ) THEN
  ALLOCATE(ZZNFS(INUCT,NMOD_CCN))
  ALLOCATE(ZZNAS(INUCT,NMOD_CCN))
  ALLOCATE(ZPRES(INUCT))
  ALLOCATE(ZRHOD(INUCT))
  ALLOCATE(ZRR(INUCT))
  ALLOCATE(ZTT(INUCT))
  ALLOCATE(ZRV(INUCT))
  ALLOCATE(ZRC(INUCT))
  ALLOCATE(ZS0(INUCT))
  ALLOCATE(ZCCL(INUCT))
  ALLOCATE(ZZDZ(INUCT))
  ALLOCATE(ZZLV(INUCT))
  ALLOCATE(ZZCPH(INUCT))
  DO JL=1,INUCT
   ZPRES(JL) = PPABSTT(I1(JL),I2(JL),I3(JL))
   ZRHOD(JL) = PRHODREF(I1(JL),I2(JL),I3(JL))
   ZRR(JL)   = PRS(I1(JL),I2(JL),I3(JL),3)
   ZTT(JL)   = PTHS(I1(JL),I2(JL),I3(JL))
   ZRV(JL)   = PRS(I1(JL),I2(JL),I3(JL),1)
   ZRC(JL)   = PRS(I1(JL),I2(JL),I3(JL),2)
   ZS0(JL)   = ZSAT(I1(JL),I2(JL),I3(JL))
   DO JMOD = 1,NMOD_CCN
         ZZNFS(JL,JMOD)        = ZNFS(I1(JL),I2(JL),I3(JL),JMOD)
         ZZNAS(JL,JMOD)        = ZNAS(I1(JL),I2(JL),I3(JL),JMOD)
   ENDDO
   ZCCL(JL)  = ZCCS(I1(JL),I2(JL),I3(JL))
   ZZDZ(JL)=ZDZ(I1(JL),I2(JL),I3(JL))
   ZZLV(JL)=ZLV(I1(JL),I2(JL),I3(JL))
   ZZCPH(JL)=ZCPH(I1(JL),I2(JL),I3(JL))
  ENDDO
  !
  !Evaporation/Condensation/activation
   CALL PROGNOS_LIMA(PTSTEP,ZZDZ,ZZLV,ZZCPH,ZPRES,ZRHOD,  &
                ZRR,ZTT,ZRV,ZRC,ZS0,ZZNAS,ZCCL,ZZNFS)
  !
!1D array to 3D array
  DO JMOD = 1, NMOD_CCN 
   ZWORK(:,:,:)  = ZNAS(:,:,:,JMOD)
   ZNAS(:,:,:,JMOD)  = UNPACK( ZZNAS(:,JMOD) ,MASK=GNUCT(:,:,:),FIELD=ZWORK(:,:,:)  )
   ZWORK(:,:,:)  = ZNFS(:,:,:,JMOD)
   ZNFS(:,:,:,JMOD) = UNPACK( ZZNFS(:,JMOD) ,MASK=GNUCT(:,:,:),FIELD=ZWORK(:,:,:)  )
  END DO
  PSVS(:,:,:,NSV_LIMA_CCN_FREE:NSV_LIMA_CCN_FREE+NMOD_CCN-1) = ZNFS(:,:,:,:)
  PSVS(:,:,:,NSV_LIMA_CCN_ACTI:NSV_LIMA_CCN_ACTI+NMOD_CCN-1) = ZNAS(:,:,:,:)
  !
  ZWORK(:,:,:)  = ZCCS(:,:,:)
  ZCCS(:,:,:)   = UNPACK( ZCCL(:),MASK=GNUCT(:,:,:),FIELD=ZWORK(:,:,:) )
  PSVS(:,:,:,NSV_LIMA_NC) = ZCCS(:,:,:)
  !
  ZWORK(:,:,:)  = PTHS(:,:,:)
  PTHS(:,:,:)   = UNPACK( ZTT(:),MASK=GNUCT(:,:,:),FIELD=ZWORK(:,:,:) )
  ZWORK(:,:,:)  = PRS(:,:,:,1)
  PRS(:,:,:,1)   = UNPACK( ZRV(:),MASK=GNUCT(:,:,:),FIELD=ZWORK(:,:,:) )
  ZWORK(:,:,:)  = PRS(:,:,:,2)
  PRS(:,:,:,2)   = UNPACK( ZRC(:),MASK=GNUCT(:,:,:),FIELD=ZWORK(:,:,:) )
  ZWORK(:,:,:)  = ZSAT(:,:,:)
  ZSAT(:,:,:)   = UNPACK( ZS0(:),MASK=GNUCT(:,:,:),FIELD=ZWORK(:,:,:) )
  !
  DEALLOCATE(ZPRES)
  DEALLOCATE(ZRHOD)
  DEALLOCATE(ZRR)
  DEALLOCATE(ZTT)
  DEALLOCATE(ZRV)
  DEALLOCATE(ZRC)
  DEALLOCATE(ZS0)
  DEALLOCATE(ZZNFS)
  DEALLOCATE(ZZNAS)
  DEALLOCATE(ZCCL)
  DEALLOCATE(ZZDZ)
!
ENDIF
!
!Computation of saturation in the meshes where there is no
!condensation/evaporation/activation
WHERE(.NOT.GNUCT(:,:,:) )
 ZRVSAT(:,:,:) = EXP(XALPW-XBETAW/PTHS(:,:,:)-XGAMW*ALOG(PTHS(:,:,:)))
 !rvsat
 ZRVSAT(:,:,:) = (XMV / XMD)*ZRVSAT(:,:,:)/(PPABSTT(:,:,:)-ZRVSAT(:,:,:))
 ZSAT(:,:,:)   = (PRS(:,:,:,1)/ZRVSAT(:,:,:))-1D0
ENDWHERE
!
!source terms /dt
PRS(:,:,:,:)   = PRS(:,:,:,:)/PTSTEP
PTHS(:,:,:)   = PTHS(:,:,:)/PTSTEP/ZEXNS(:,:,:)
ZSAT(:,:,:)   = ZSAT(:,:,:)+1.0
PSVS(:,:,:,NSV_LIMA_SPRO) = ZSAT(:,:,:)
PSVS(:,:,:,:) = PSVS(:,:,:,:)/PTSTEP
!
IF (ANY(PRS(:,:,:,2)+PRS(:,:,:,1) < 0.) .AND. NVERB>5) THEN
  WRITE(*,*) 'LIMA_NOTADJUST:  negative values of total water (reset to zero)'
  WRITE(*,*) '  location of minimum:', MINLOC(PRS(:,:,:,2)+PRS(:,:,:,1))
  WRITE(*,*) '  value of minimum   :', MINVAL(PRS(:,:,:,2)+PRS(:,:,:,1))
  CALL PRINT_MSG(NVERB_FATAL,'GEN','LIMA_NOTADJUST','')
END IF
!
!*              compute the cloud fraction PCLDFR
!
WHERE (PRS(:,:,:,2) > 0. )
    ZW1(:,:,:)  = 1.
ELSEWHERE
    ZW1(:,:,:)  = 0. 
ENDWHERE 
IF ( SIZE(PSRCS,3) /= 0 ) THEN
    PSRCS(:,:,:) = ZW1(:,:,:) 
END IF
!
IF ( HRAD /= 'NONE' ) THEN
     PCLDFR(:,:,:) = ZW1(:,:,:)
END IF
!
ZW1(:,:,:)=0.
IF (SIZE(PRS,4)>3) ZW1(:,:,:)=ZW1(:,:,:) + PRS(:,:,:,4)
WHERE (ZW1(:,:,:) > 1.E-15)
   PICEFR(:,:,:)  = 1.
ELSEWHERE
   PICEFR(:,:,:)  = 0.
ENDWHERE
ZW1(:,:,:)=0.
IF (SIZE(PRS,4)>2) ZW1(:,:,:)=ZW1(:,:,:) + PRS(:,:,:,3)
IF (SIZE(PRS,4)>4) ZW1(:,:,:)=ZW1(:,:,:) + PRS(:,:,:,5)
IF (SIZE(PRS,4)>5) ZW1(:,:,:)=ZW1(:,:,:) + PRS(:,:,:,6)
WHERE (ZW1(:,:,:) > 1.E-15)
   PRAINFR(:,:,:)  = 1.
ELSEWHERE
   PRAINFR(:,:,:)  = 0.
ENDWHERE
!
IF ( tpfile%lopened ) THEN
  ZW(:,:,:)=SUM(ZNAS,4)-ZW(:,:,:)
  TZFIELD = TFIELDMETADATA(    &
    CMNHNAME   = 'NACT',       &
    CSTDNAME   = '',           &
    CLONGNAME  = 'NACT',       &
    CUNITS     = 'kg-1',       &
    CDIR       = 'XY',         &
    CCOMMENT   = 'X_Y_Z_NACT', &
    NGRID      = 1,            &
    NTYPE      = TYPEREAL,     &
    NDIMS      = 3,            &
    LTIMEDEP   = .TRUE.        )
  CALL IO_Field_write(TPFILE,TZFIELD,ZW)
END IF
!
!*       7.  STORE THE BUDGET TERMS
!            ----------------------
!
if ( nbumod == kmi .and. lbu_enable ) then
  if ( lbudget_th ) call Budget_store_end( tbudgets(NBUDGET_TH), 'CEDS', pths(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rv ) call Budget_store_end( tbudgets(NBUDGET_RV), 'CEDS', prs(:, :, :, 1) * prhodj(:, :, :) )
  if ( lbudget_rc ) call Budget_store_end( tbudgets(NBUDGET_RC), 'CEDS', prs(:, :, :, 2) * prhodj(:, :, :) )
  if ( lbudget_ri ) call Budget_store_end( tbudgets(NBUDGET_RI), 'CEDS', prs(:, :, :, 4) * prhodj(:, :, :) )
  if ( lbudget_sv ) then
    if (nmom_c.ge.2) then
      call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nc), 'CEDS', psvs(:, :, :, nsv_lima_nc) * prhodj(:, :, :) )
      do jl = nsv_lima_ccn_free, nsv_lima_ccn_free + nmod_ccn - 1
        idx = NBUDGET_SV1 - 1 + jl
        call Budget_store_end( tbudgets(idx), 'CEDS', psvs(:, :, :, jl) * prhodj(:, :, :) )
      end do
      do jl = nsv_lima_ccn_acti, nsv_lima_ccn_acti + nmod_ccn - 1
        idx = NBUDGET_SV1 - 1 + jl
        call Budget_store_end( tbudgets(idx), 'CEDS', psvs(:, :, :, jl) * prhodj(:, :, :) )
      end do
    end if
!     if ( lscav .and. laero_mass ) &
!       call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_scavmass), 'CEDS', psvs(:, :, :, nsv_lima_scavmass) &
!                                                                                     * prhodj(:, :, :) )
!     if ( nmom_i.ge.2 ) then
!       call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ni), 'CEDS', psvs(:, :, :, nsv_lima_ni) * prhodj(:, :, :) )
!       do jl = 1, nsv_lima_ifn_free, nsv_lima_ifn_free + nmod_ifn - 1
!         idx = NBUDGET_SV1 - 1 + jl
!         call Budget_store_end( tbudgets(idx), 'CEDS', psvs(:, :, :, jl) * prhodj(:, :, :) )
!       end do
!       do jl = 1, nsv_lima_ifn_nucl, nsv_lima_ifn_nucl + nmod_ifn - 1
!         idx = NBUDGET_SV1 - 1 + jl
!         call Budget_store_end( tbudgets(idx), 'CEDS', psvs(:, :, :, jl) * prhodj(:, :, :) )
!       end do
!       do jl = 1, nsv_lima_imm_nucl, nsv_lima_imm_nucl + nmod_ifn - 1
!         idx = NBUDGET_SV1 - 1 + jl
!         call Budget_store_end( tbudgets(idx), 'CEDS', psvs(:, :, :, jl) * prhodj(:, :, :) )
!       end do
!     end if
      call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_spro), 'CEDS', psvs(:, :, :, nsv_lima_spro) * prhodj(:, :, :) )
  end if
end if
!
END SUBROUTINE LIMA_NOTADJUST
