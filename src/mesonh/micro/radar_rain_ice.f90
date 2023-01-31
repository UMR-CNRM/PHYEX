!MNH_LIC Copyright 1996-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.

!      ##########################
       MODULE MODI_RADAR_RAIN_ICE 
!      ##########################
!
INTERFACE
      SUBROUTINE RADAR_RAIN_ICE(PRT,PCIT,PRHODREF,PTEMP,PRARE,PVDOP,PRZDR,PRKDP,&
                             PCRT,PCST,PCGT,PCHT)
!
REAL,  DIMENSION(:,:,:,:), INTENT(IN)  :: PRT  ! microphysical  mix. ratios at t
REAL,  DIMENSION(:,:,:),   INTENT(IN)  :: PCIT ! pristine ice concentration at t
REAL,  DIMENSION(:,:,:),   INTENT(IN)  :: PRHODREF ! density of the ref. state
REAL,  DIMENSION(:,:,:),   INTENT(IN)  :: PTEMP    ! air temperature
!
REAL,  DIMENSION(:,:,:), INTENT(OUT) :: PRARE! radar reflectivity in dBZ
REAL,  DIMENSION(:,:,:), INTENT(OUT) :: PVDOP! radar Doppler fall speed
REAL,  DIMENSION(:,:,:), INTENT(OUT) :: PRZDR! radar differential reflectivity
                                             ! H-V in dBZ
REAL,  DIMENSION(:,:,:), INTENT(OUT) :: PRKDP! radar differential phase shift
                                             ! H-V in degree/km
REAL,  DIMENSION(:,:,:),   INTENT(IN),OPTIONAL  :: PCRT ! rain concentration at t
REAL,  DIMENSION(:,:,:),   INTENT(IN),OPTIONAL  :: PCST ! snow concentration at t !
REAL,  DIMENSION(:,:,:),   INTENT(IN),OPTIONAL  :: PCGT ! graupel concentration at t ! 
REAL,  DIMENSION(:,:,:),   INTENT(IN),OPTIONAL  :: PCHT ! hail concentration at t !                                             
!
END SUBROUTINE RADAR_RAIN_ICE
!
END INTERFACE
!
END MODULE MODI_RADAR_RAIN_ICE
!     #########################################################################
      SUBROUTINE RADAR_RAIN_ICE(PRT,PCIT,PRHODREF,PTEMP,PRARE,PVDOP,PRZDR,PRKDP,PCRT,PCST,PCGT,PCHT)
!     #########################################################################
!
!!****  *RADAR_RAIN_ICE * - computes some pertinent radar parameters
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the equivalent reflectivity,
!!    the Doppler reflectivity and the H and V polarized reflectivities of a 
!!    mixed phase cloud.
!!
!!**  METHOD
!!    ------
!!      The reflectivities are computed using the n(D) * D**6 formula. The 
!!    equivalent reflectiviy is the sum of the reflectivity produced by the
!!    the raindrops and the equivalent reflectivities of the ice crystals.
!!    The latter are computed using the melted diameter. The Doppler 
!!    reflectivity is the 'fall speed'-moment of individual particle
!!    reflectivity. Ice crystal are assumed to have no preferred orientation.
!!    the Z_VV formula is taken from Brandes et al. (MWR, 1995).
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST
!!        XPI                  !
!!        XRHOLW               ! Liquid water density
!!      Module MODD_RAIN_ICE_DESCR
!!      Module MODD_RAIN_ICE_PARAM
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation ( routine RADAR_RAIN_ICE )
!!      Smith P.L., 1984: Equivalent Radar Reflectivity Factors for Snow and
!!                        Ice Particles, JCAM, 23, 1258-1260.
!!      Andsager K., K. V. Beard, and N. F. Laird, 1999: Laboratory Measurements
!!                        of Axis Ratio for Large Raindrops, JAS, 56, 2673-2683.
!!
!!    AUTHOR
!!    ------
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/05/96
!!                  03/12/96 change the arg. list
!!                  15/12/00 change the reflectivity factor
!!                  01/07/02 (E.Richard) bug in reflectivity formula
!!                          for graupeln when warmer than 0°C
!!                  19/12/00 (JP Pinty) change the hailstone reflectivity
!!                  19/05/04 (JP Pinty) add ZDR and KDP for raindops at 10.71 cm
!! J.-P. Chaboureau 17/06/10 bug correction in reflectivity calculation of icy hydrometeors
!! J.-P. Chaboureau 03/02/12 set undef values for radar reflectivities
!!       O. Caumont 09/04/14 correction of ZDR calculation
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
USE MODD_REF
USE MODD_PARAM_ICE,      ONLY: LSNOW_T_I=>LSNOW_T
USE MODD_RAIN_ICE_DESCR, ONLY: XALPHAR_I=>XALPHAR,XNUR_I=>XNUR,XLBEXR_I=>XLBEXR,&
                               XLBR_I=>XLBR,XCCR_I=>XCCR,XBR_I=>XBR,XAR_I=>XAR,&
                               XALPHAC_I=>XALPHAC,XNUC_I=>XNUC,&
                               XLBC_I=>XLBC,XBC_I=>XBC,XAC_I=>XAC,&
                               XALPHAC2_I=>XALPHAC2,XNUC2_I=>XNUC2,&
                               XALPHAS_I=>XALPHAS,XNUS_I=>XNUS,XLBEXS_I=>XLBEXS,&
                               XLBS_I=>XLBS,XCCS_I=>XCCS,XNS_I=>XNS,XAS_I=>XAS,XBS_I=>XBS,XCXS_I=>XCXS,&
                               XALPHAG_I=>XALPHAG,XNUG_I=>XNUG,XDG_I=>XDG,XLBEXG_I=>XLBEXG,&
                               XLBG_I=>XLBG,XCCG_I=>XCCG,XAG_I=>XAG,XBG_I=>XBG,XCXG_I=>XCXG,XCG_I=>XCG,&
                               XALPHAI_I=>XALPHAI,XNUI_I=>XNUI,XDI_I=>XDI,XLBEXI_I=>XLBEXI,&
                               XLBI_I=>XLBI,XAI_I=>XAI,XBI_I=>XBI,XC_I_I=>XC_I,XCS_I=>XCS,XDS_I=>XDS,&
                               XRTMIN_I=>XRTMIN,XCONC_LAND,XCONC_SEA,XCR_I=>XCR,XDR_I=>XDR,&
                               XAH_I=>XAH,XLBH_I=>XLBH,XLBEXH_I=>XLBEXH,XCCH_I=>XCCH,&
                               XALPHAH_I=>XALPHAH,XNUH_I=>XNUH,XCXH_I=>XCXH,XDH_I=>XDH,XCH_I=>XCH,XBH_I=>XBH,        &
                               XLBDAS_MAX_I=>XLBDAS_MAX,XLBDAS_MIN_I=>XLBDAS_MIN,XTRANS_MP_GAMMAS_I=>XTRANS_MP_GAMMAS
USE MODD_PARAM_LIMA_WARM, ONLY: XLBEXR_L=>XLBEXR,XLBR_L=>XLBR,XBR_L=>XBR,XAR_L=>XAR,&
                                XBC_L=>XBC,XAC_L=>XAC,XCR_L=>XCR,XDR_L=>XDR
USE MODD_PARAM_LIMA_COLD, ONLY: XDI_L=>XDI,XLBEXI_L=>XLBEXI,XLBI_L=>XLBI,XAI_L=>XAI,XBI_L=>XBI,XC_I_L=>XC_I,&
                                XLBEXS_L=>XLBEXS,XLBS_L=>XLBS,XCCS_L=>XCCS,XNS_L=>XNS,&
                                XAS_L=>XAS,XBS_L=>XBS,XCXS_L=>XCXS,XCS_L=>XCS,XDS_L=>XDS,&
                                XLBDAS_MAX_L=>XLBDAS_MAX,XLBDAS_MIN_L=>XLBDAS_MIN,XTRANS_MP_GAMMAS_L=>XTRANS_MP_GAMMAS

USE MODD_PARAM_LIMA_MIXED, ONLY:XDG_L=>XDG,XLBEXG_L=>XLBEXG,XLBG_L=>XLBG,XCCG_L=>XCCG,&
                                XAG_L=>XAG,XBG_L=>XBG,XCXG_L=>XCXG,XCG_L=>XCG,&
                                XAH_L=>XAH,XLBH_L=>XLBH,XLBEXH_L=>XLBEXH,XCCH_L=>XCCH,&
                                XCXH_L=>XCXH,XDH_L=>XDH,XCH_L=>XCH,XALPHAH_L=>XALPHAH,XNUH_L=>XNUH,XBH_L=>XBH

USE MODD_PARAM_LIMA, ONLY: XALPHAR_L=>XALPHAR,XNUR_L=>XNUR,XALPHAS_L=>XALPHAS,XNUS_L=>XNUS,&
                           XALPHAG_L=>XALPHAG,XNUG_L=>XNUG, XALPHAI_L=>XALPHAI,XNUI_L=>XNUI,&
                           XRTMIN_L=>XRTMIN,XALPHAC_L=>XALPHAC,XNUC_L=>XNUC,LSNOW_T_L=>LSNOW_T,NMOM_S,NMOM_G,NMOM_H                      
USE MODD_PARAMETERS
USE MODD_PARAM_n, ONLY : CCLOUD
USE MODD_LUNIT
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
REAL,  DIMENSION(:,:,:,:), INTENT(IN)  :: PRT  ! microphysical  mix. ratios at t
REAL,  DIMENSION(:,:,:),   INTENT(IN)  :: PCIT ! pristine ice concentration at t
REAL,  DIMENSION(:,:,:),   INTENT(IN)  :: PRHODREF ! density of the ref. state
REAL,  DIMENSION(:,:,:),   INTENT(IN)  :: PTEMP    ! air temperature
!
REAL,  DIMENSION(:,:,:), INTENT(OUT) :: PRARE! radar reflectivity in dBZ
REAL,  DIMENSION(:,:,:), INTENT(OUT) :: PVDOP! radar Doppler fall speed
REAL,  DIMENSION(:,:,:), INTENT(OUT) :: PRZDR! radar differential reflectivity
                                             ! H-V in dBZ
REAL,  DIMENSION(:,:,:), INTENT(OUT) :: PRKDP! radar differential phase shift
                                             ! H-V in degree/km
REAL,  DIMENSION(:,:,:),   INTENT(IN),OPTIONAL  :: PCRT ! rain concentration at t
REAL,  DIMENSION(:,:,:),   INTENT(IN),OPTIONAL  :: PCST ! snow concentration at t     
REAL,  DIMENSION(:,:,:),   INTENT(IN),OPTIONAL  :: PCGT ! graupel concentration at t  
REAL,  DIMENSION(:,:,:),   INTENT(IN),OPTIONAL  :: PCHT ! hail concentration at t   
!
!*       0.2   Declarations of local variables :
!
INTEGER :: IKB           ! Coordinates of the first physical points along z
INTEGER :: IND           ! Number of interval to integrate the kernels
REAL :: ZALPHA, ZNU, ZP  ! Parameters to compute the value of the p_moment
       			 ! of the generalized Gamma function
REAL :: ZDINFTY          ! Factor used to define the "infinite" diameter
!
REAL :: ZCXR=-1.0                     ! for rain N ~ 1/N_0 
                                      ! (in Kessler parameterization)
REAL :: ZSLOPE, ZINTERCEPT, ZEXPONENT ! parameters defining the mean axis ratio
       				      ! functionnal
REAL :: ZDMELT_FACT                   ! factor used to compute the equivalent
	                	      ! melted diameter
REAL :: ZEQICE                        ! factor used to convert the ice crystals
               			      ! reflectivity into an equivalent  liquid
		                      ! water reflectivity (from Smith, JCAM 84)
REAL :: ZEXP                          ! anciliary parameter
REAL :: ZRHO00                        ! Surface reference air density
!
LOGICAL, DIMENSION(SIZE(PTEMP,1),SIZE(PTEMP,2),SIZE(PTEMP,3)) :: GRAIN
REAL,    DIMENSION(SIZE(PTEMP,1),SIZE(PTEMP,2),SIZE(PTEMP,3)) :: ZLBDA 
                                      ! slope distribution parameter
REAL,    DIMENSION(SIZE(PTEMP,1),SIZE(PTEMP,2),SIZE(PTEMP,3)) :: ZN 
                                      ! number concentration
REAL,    DIMENSION(SIZE(PTEMP,1),SIZE(PTEMP,2),SIZE(PTEMP,3)) :: ZW
REAL,    DIMENSION(SIZE(PTEMP,1),SIZE(PTEMP,2),SIZE(PTEMP,3)) :: ZREFL_MELT_CONV
INTEGER                                                       :: JLBDA
REAL                                                          :: ZFRAC_WATER
!
LOGICAL  :: GFLAG   ! Logical flag for printing the constatnts on the output
                    ! listing
!
REAL     ::   ZR0, ZR1, ZR2 ! r(D) parameters
!REAL     ::   ZREXP, ZSCALE ! parameters to compute Zhh from Zvv
REAL     ::   Z1, Z2, Z3    ! expansion coefficients
! 
INTEGER  :: II, IJ, IK 
!
!REAL                           :: ZA,ZB,ZCX,ZALPHA,ZNU,ZLB,ZLBEX,ZRHOHYD   ! generic microphysical parameters
!REAL,DIMENSION(:),ALLOCATABLE  :: ZRTMIN ! local values for XRTMIN

REAL,    DIMENSION(SIZE(PTEMP,1),SIZE(PTEMP,2),SIZE(PTEMP,3)) :: ZLB_L
REAL :: ZLB,ZLBEX ,  ZCC,ZCX,ZC,ZD

!-------------------------------------------------------------------------------
!
!
!*       1.     FUNCTION STATEMENTS
!   	        -------------------
!
!
!*       1.1    p_moment of the Generalized GAMMA function
!
! Recall that MOMG(ZALPHA,ZNU,ZP)=GAMMA(ZNU+ZP/ZALPHA)/GAMMA(ZNU)
!
!
!-------------------------------------------------------------------------------
!
!
!        2.     INTIALIZE OUTPUT LISTING AND OTHER ARRAYS
!               -----------------------------------------
!
!
PRARE(:,:,:) = 0.0    ! radar reflectivity
PVDOP(:,:,:) = 0.0    ! radar Doppler fall speed
PRZDR(:,:,:) = 0.0    ! radar differential reflectivity
PRKDP(:,:,:) = 0.0    ! radar differential phase shift
!
!-------------------------------------------------------------------------------
!
!
!*       3.     RAINDROPS
!               ---------
!
IF (SIZE(PRT,4) >= 3) THEN
  IND        = 50
  ZSLOPE     = 0.62     ! the mean axis ratio function writes as r**ZEXPONENT
  ZINTERCEPT = 1.03     !                  with
  ZEXPONENT  = 7.0/3.0  ! r = ZSLOPE*D+ZINTERCEPT where D is the drop diameter
  ZDINFTY    = 20.0
!
! The raindrop aspect ratio is given by Andsager et al. (1999)
!   r(D) = ZR0 + ZR1*D + ZR2*D**2
!
  ZR0 = 1.012
  ZR1 = -0.144E2
  ZR2 = -1.03E4
!
!  ZREXP  = 7.0/3.0
!  ZSCALE = ZR0**ZREXP
!  Z1     = ZREXP*(ZR1/ZR0)
!  Z2     = ZREXP*(ZR2/ZR0)+ZREXP*(ZREXP-1.0)*0.5*(ZR1/ZR0)**2
  Z1=.97
  Z2=.64
  Z3=7.8
!
  ZLBDA(:,:,:) = 0.0
  IF (CCLOUD == 'LIMA') THEN
    GRAIN(:,:,:) =( (PRT(:,:,:,3).GT.XRTMIN_L(3)).AND.  PCRT(:,:,:).GT.0.0)     
    ZLBEX=1.0/(-XBR_L)
    ZLB_L(:,:,:)=( XAR_L*PCRT(:,:,:)*PRHODREF(:,:,:)*MOMG(XALPHAR_L,XNUR_L,XBR_L) )**(-ZLBEX)
    WHERE( GRAIN(:,:,:) )
      ZLBDA(:,:,:) =ZLB_L(:,:,:) *( PRHODREF(:,:,:)*PRT(:,:,:,3) )**(ZLBEX)
      PRARE(:,:,:) = 1.E18*PCRT(:,:,:)*PRHODREF(:,:,:)*(ZLBDA(:,:,:)**(-6.0))*MOMG(XALPHAR_L,XNUR_L,6.0)
      PVDOP(:,:,:) = 1.E18*PCRT(:,:,:)*PRHODREF(:,:,:)*XCR_L*(ZLBDA(:,:,:)**(-6.0-XDR_L))              &
                  *MOMG(XALPHAR_L,XNUR_L,6.0+XDR_L)
      PRZDR(:,:,:) = Z1+Z2*(PRHODREF(:,:,:)*PRT(:,:,:,3))**(-ZLBEX)+Z3*(PRHODREF(:,:,:)*PRT(:,:,:,3))**(-2.*ZLBEX)
      PRZDR(:,:,:) = 10.0*LOG10( PRZDR(:,:,:) ) ! now in dBZ
      PRKDP(:,:,:) = 6.7E3*( PRHODREF(:,:,:)*PRT(:,:,:,3) )*                     &
                  (-ZR1*(MOMG(XALPHAR_L,XNUR_L,4.0)/MOMG(XALPHAR_L,XNUR_L,3.0))*(1.0/ZLBDA(:,:,:))   &
                   -ZR2*(MOMG(XALPHAR_L,XNUR_L,5.0)/MOMG(XALPHAR_L,XNUR_L,3.0))*(1.0/ZLBDA(:,:,:)**2))
                                              ! in degree/km
    END WHERE    
  ELSE
    IF( SIZE(PRT,4) == 3 ) THEN
      GRAIN(:,:,:) = PRT(:,:,:,3).GT.1.0E-15
      ZCC = 1.E7; ZLBEX = -0.25          ! Marshall-Palmer law
      ZALPHA = 1.0; ZNU = 1.0            ! Marshall-Palmer law
      ZC = 842.; ZD = 0.8                ! Raindrop fall-speed
      ZLB = (XPI*XRHOLW*ZCC)**(-XLBEXR_I)
    ELSE
      ZLB=XLBR_I
      ZLBEX=XLBEXR_I
      ZCC=XCCR_I
      ZALPHA=XALPHAR_I
      ZNU=XNUR_I
      ZC=XCR_I
      ZD=XDR_I
      GRAIN(:,:,:) = PRT(:,:,:,3).GT.XRTMIN_I(3)
    END IF
    WHERE( GRAIN(:,:,:) )
      ZLBDA(:,:,:) = ZLB*( PRHODREF(:,:,:)*PRT(:,:,:,3) )**ZLBEX
      PRARE(:,:,:) = 1.E18*ZCC*(ZLBDA(:,:,:)**(ZCXR-6.0))*MOMG(ZALPHA,ZNU,6.0)
      PVDOP(:,:,:) = 1.E18*ZCC*ZC*(ZLBDA(:,:,:)**(ZCXR-6.0-ZD))              &
                     *MOMG(ZALPHA,ZNU,6.0+ZD)
      PRZDR(:,:,:) = Z1+Z2*(PRHODREF(:,:,:)*PRT(:,:,:,3))**(-ZLBEX)+Z3*(PRHODREF(:,:,:)*PRT(:,:,:,3))**(-2.*ZLBEX)
      PRZDR(:,:,:) = 10.0*LOG10( PRZDR(:,:,:) ) ! now in dBZ
      PRKDP(:,:,:) = 6.7E3*( PRHODREF(:,:,:)*PRT(:,:,:,3) )*                     &
                    (-ZR1*(MOMG(ZALPHA,ZNU,4.0)/MOMG(ZALPHA,ZNU,3.0))*(1.0/ZLBDA(:,:,:))   &
                     -ZR2*(MOMG(ZALPHA,ZNU,5.0)/MOMG(ZALPHA,ZNU,3.0))*(1.0/ZLBDA(:,:,:)**2))
                                              ! in degree/km
    END WHERE
  ENDIF

END IF


!
!*       4.     PRISTINE ICE
!               ------------
!
IF (SIZE(PRT,4) >= 4) THEN
  ZEQICE = 0.224
  IF (CCLOUD == 'LIMA') THEN
    ZDMELT_FACT = ( (6.0*XAI_L)/(XPI*XRHOLW) )**(2.0)
    ZEXP = 2.0*XBI_L
    WHERE( PRT(:,:,:,4).GT.XRTMIN_L(4) .AND. PCIT(:,:,:).GT.0.0 )
      ZLBDA(:,:,:) = XLBI_L**(XLBEXI_L)* (PRT(:,:,:,4)/PCIT(:,:,:))**(-XLBEXI_L)
      ZW(:,:,:) = ZEQICE*ZDMELT_FACT *1.E18*PRHODREF(:,:,:)*PCIT(:,:,:)*(ZLBDA(:,:,:)**(-ZEXP))*MOMG(XALPHAI_L,XNUI_L,ZEXP)
      PVDOP(:,:,:) = PVDOP(:,:,:)+ZEQICE*ZDMELT_FACT*MOMG(XALPHAI_L,XNUI_L,ZEXP+XDI_L) &
                     *1.E18*PRHODREF(:,:,:)*PCIT(:,:,:)*XC_I_L*(ZLBDA(:,:,:)**(-ZEXP-XDI_L))
      PRARE(:,:,:) = PRARE(:,:,:) + ZW(:,:,:)
    END WHERE
  ELSE
    ZDMELT_FACT = ( (6.0*XAI_I)/(XPI*XRHOLW) )**(2.0)
    ZEXP = 2.0*XBI_I
    WHERE( PRT(:,:,:,4).GT.XRTMIN_I(4) .AND. PCIT(:,:,:).GT.0.0 )
      ZLBDA(:,:,:) = XLBI_I*( PRHODREF(:,:,:)*PRT(:,:,:,4)/PCIT(:,:,:) )**XLBEXI_I
      ZW(:,:,:) = ZEQICE*ZDMELT_FACT*1.E18*PCIT(:,:,:)*(ZLBDA(:,:,:)**(-ZEXP))*MOMG(XALPHAI_I,XNUI_I,ZEXP)
      PVDOP(:,:,:) = PVDOP(:,:,:)+ZEQICE*ZDMELT_FACT*MOMG(XALPHAI_I,XNUI_I,ZEXP+XDI_I) &
                     *1.E18*PCIT(:,:,:)*XC_I_I*(ZLBDA(:,:,:)**(-ZEXP-XDI_I))
      PRARE(:,:,:) = PRARE(:,:,:) + ZW(:,:,:)
    END WHERE
  END IF
END IF
!
!*       5.     SNOW/AGGREGATES
!               ---------------
!
IF (SIZE(PRT,4) >= 5) THEN
   IF ( (CCLOUD=='LIMA' .AND. LSNOW_T_L) ) THEN
      ZDMELT_FACT = ( (6.0*XAS_L)/(XPI*XRHOLW) )**(2.0)
      ZEXP = 2.0*XBS_L
      WHERE(PTEMP(:,:,:)>263.15 .AND. PRT(:,:,:,5).GT.XRTMIN_L(5))
         ZLBDA(:,:,:) = MAX(MIN(XLBDAS_MAX_L, 10**(14.554-0.0423*PTEMP(:,:,:))),XLBDAS_MIN_L)*XTRANS_MP_GAMMAS_L
      END WHERE
      WHERE(PTEMP(:,:,:)<=263.15 .AND. PRT(:,:,:,5).GT.XRTMIN_L(5))
         ZLBDA(:,:,:) = MAX(MIN(XLBDAS_MAX_L, 10**(6.226 -0.0106*PTEMP(:,:,:))),XLBDAS_MIN_L)*XTRANS_MP_GAMMAS_L
      END WHERE
      IF (NMOM_S.GE.2) THEN
         ZN(:,:,:)=PCST(:,:,:)
      ELSE
         WHERE( PRT(:,:,:,5).GT.XRTMIN_L(5) )
            ZN(:,:,:)=XNS_L*PRHODREF(:,:,:)*PRT(:,:,:,5)*ZLBDA(:,:,:)**XBS_L
         END WHERE
      END IF
      WHERE( PRT(:,:,:,5).GT.XRTMIN_L(5) )
         ZW(:,:,:) = ZEQICE*ZDMELT_FACT                                             &
              *1.E18*PRHODREF(:,:,:)*ZN(:,:,:)*(ZLBDA(:,:,:)**(-ZEXP))*MOMG(XALPHAS_L,XNUS_L,ZEXP)
         PVDOP(:,:,:) = PVDOP(:,:,:)+ZEQICE*ZDMELT_FACT*MOMG(XALPHAS_L,XNUS_L,ZEXP+XDS_L) &
              *1.E18*PRHODREF(:,:,:)*ZN(:,:,:)*XCS_L*(ZLBDA(:,:,:)**(-ZEXP-XDS_L))
         PRARE(:,:,:) = PRARE(:,:,:) + ZW(:,:,:)
      END WHERE
   ELSEIF ( (CCLOUD=='ICE3' .AND. LSNOW_T_I) ) THEN
    ZDMELT_FACT = ( (6.0*XAS_I)/(XPI*XRHOLW) )**(2.0)
    ZEXP = 2.0*XBS_I
    WHERE(PTEMP(:,:,:)>263.15 .AND. PRT(:,:,:,5).GT.XRTMIN_I(5))
       ZLBDA(:,:,:) = MAX(MIN(XLBDAS_MAX_I, 10**(14.554-0.0423*PTEMP(:,:,:))),XLBDAS_MIN_I)*XTRANS_MP_GAMMAS_I
    END WHERE
    WHERE(PTEMP(:,:,:)<=263.15 .AND. PRT(:,:,:,5).GT.XRTMIN_I(5))
       ZLBDA(:,:,:) = MAX(MIN(XLBDAS_MAX_I, 10**(6.226- 0.0106*PTEMP(:,:,:))),XLBDAS_MIN_I)*XTRANS_MP_GAMMAS_I
    END WHERE
    ZN(:,:,:)=XNS_I*PRHODREF(:,:,:)*PRT(:,:,:,5)*ZLBDA(:,:,:)**XBS_I
    WHERE( PRT(:,:,:,5).GT.XRTMIN_I(5) )
      ZW(:,:,:) = ZEQICE*ZDMELT_FACT                                             &
                  *1.E18*PRHODREF(:,:,:)*ZN(:,:,:)*(ZLBDA(:,:,:)**(-ZEXP))*MOMG(XALPHAS_I,XNUS_I,ZEXP)
      PVDOP(:,:,:) = PVDOP(:,:,:)+ZEQICE*ZDMELT_FACT*MOMG(XALPHAS_I,XNUS_I,ZEXP+XDS_I) &
                     *1.E18*PRHODREF(:,:,:)*ZN(:,:,:)*XCS_I*(ZLBDA(:,:,:)**(-ZEXP-XDS_I))
      PRARE(:,:,:) = PRARE(:,:,:) + ZW(:,:,:)
    END WHERE
  ELSEIF (CCLOUD=='LIMA') THEN
    ZDMELT_FACT = ( (6.0*XAS_L)/(XPI*XRHOLW) )**(2.0)
    ZEXP = 2.0*XBS_L
    if (NMOM_S.GE.2) then
      WHERE( PRT(:,:,:,5).GT.XRTMIN_L(5) .AND. PCST(:,:,:).GT.0.0)      
        ZLBDA(:,:,:) = XLBS_L**(XLBEXS_L)*(PRT(:,:,:,5)/PCST(:,:,:))**(-XLBEXS_L)  
        ZW(:,:,:) = ZEQICE*ZDMELT_FACT                                             &
                    *1.E18*PRHODREF(:,:,:)*PCST(:,:,:)*(ZLBDA(:,:,:)**(-ZEXP))*MOMG(XALPHAS_L,XNUS_L,ZEXP)
        PVDOP(:,:,:) = PVDOP(:,:,:)+ZEQICE*ZDMELT_FACT*MOMG(XALPHAS_L,XNUS_L,ZEXP+XDS_L) &
                       *1.E18*PRHODREF(:,:,:)*PCST(:,:,:)*XCS_L*(ZLBDA(:,:,:)**(-ZEXP-XDS_L))
        PRARE(:,:,:) = PRARE(:,:,:) + ZW(:,:,:)
      END WHERE             
    else
      WHERE( PRT(:,:,:,5).GT.XRTMIN_L(5) )
        ZLBDA(:,:,:) = XLBS_L*( PRHODREF(:,:,:)*PRT(:,:,:,5) )**XLBEXS_L
        ZW(:,:,:) = ZEQICE*ZDMELT_FACT                                             &
                  *1.E18*XCCS_L*(ZLBDA(:,:,:)**(XCXS_L-ZEXP))*MOMG(XALPHAS_L,XNUS_L,ZEXP)
        PVDOP(:,:,:) = PVDOP(:,:,:)+ZEQICE*ZDMELT_FACT*MOMG(XALPHAS_L,XNUS_L,ZEXP+XDS_L) &
                     *1.E18*XCCS_L*XCS_L*(ZLBDA(:,:,:)**(XCXS_L-ZEXP-XDS_L))
        PRARE(:,:,:) = PRARE(:,:,:) + ZW(:,:,:)
      END WHERE
    end if
  ELSE
    ZDMELT_FACT = ( (6.0*XAS_I)/(XPI*XRHOLW) )**(2.0)
    ZEXP = 2.0*XBS_I
    WHERE( PRT(:,:,:,5).GT.XRTMIN_I(5) )
      ZLBDA(:,:,:) = XLBS_I*( PRHODREF(:,:,:)*PRT(:,:,:,5) )**XLBEXS_I
      ZW(:,:,:) = ZEQICE*ZDMELT_FACT                                             &
                  *1.E18*XCCS_I*(ZLBDA(:,:,:)**(XCXS_I-ZEXP))*MOMG(XALPHAS_I,XNUS_I,ZEXP)
      PVDOP(:,:,:) = PVDOP(:,:,:)+ZEQICE*ZDMELT_FACT*MOMG(XALPHAS_I,XNUS_I,ZEXP+XDS_I) &
                     *1.E18*XCCS_I*XCS_I*(ZLBDA(:,:,:)**(XCXS_I-ZEXP-XDS_I))
      PRARE(:,:,:) = PRARE(:,:,:) + ZW(:,:,:)
    END WHERE
  ENDIF
END IF
!
!*       6.     GRAUPELN
!               --------
!
IF (SIZE(PRT,4) >= 6) THEN
  IF (CCLOUD=='LIMA') THEN
    ZFRAC_WATER = 0.14
    ZDMELT_FACT = ( (6.0*XAG_L)/(XPI*XRHOLW) )**(2.0)
    WHERE( PTEMP(:,:,:).GT.XTT )
      ZREFL_MELT_CONV(:,:,:) = ((1.0-ZFRAC_WATER)*ZEQICE+ZFRAC_WATER)*ZDMELT_FACT
    ELSEWHERE
      ZREFL_MELT_CONV(:,:,:) = ZEQICE*ZDMELT_FACT
    END WHERE
!
    ZEXP = 2.0*XBG_L
    if(NMOM_G.GE.2) then
      WHERE( PRT(:,:,:,6).GT.XRTMIN_L(6) .AND. PCGT(:,:,:).GT.1.0E-3 )      
        ZLBDA(:,:,:) = XLBG_L**(XLBEXG_L)*(PRT(:,:,:,6)/PCGT(:,:,:))**(-XLBEXG_L)  
        ZW(:,:,:)    = ZREFL_MELT_CONV(:,:,:)          &
                  *1.E18*PRHODREF(:,:,:)*PCGT(:,:,:)*(ZLBDA(:,:,:)**(-ZEXP))*MOMG(XALPHAG_L,XNUG_L,ZEXP)   
        PVDOP(:,:,:) = PVDOP(:,:,:) +                                            &
                       ZREFL_MELT_CONV(:,:,:)*1.E18                      &
                       *1.E18*PRHODREF(:,:,:)*PCGT(:,:,:)*XCG_L*(ZLBDA(:,:,:)**(-ZEXP-XDG_L))
        PRARE(:,:,:) = PRARE(:,:,:) + ZW(:,:,:)
      END WHERE
    else
      WHERE( PRT(:,:,:,6).GT.XRTMIN_L(6) )
        ZLBDA(:,:,:) = XLBG_L*( PRHODREF(:,:,:)*PRT(:,:,:,6) )**XLBEXG_L
        ZW(:,:,:)    = ZREFL_MELT_CONV(:,:,:)*1.E18*XCCG_L*                &
                      (ZLBDA(:,:,:)**(XCXG_L-ZEXP))*MOMG(XALPHAG_L,XNUG_L,ZEXP)
        PVDOP(:,:,:) = PVDOP(:,:,:) +                                            &
                       ZREFL_MELT_CONV(:,:,:)*1.E18*XCCG_L*XCG_L*                    &
                       (ZLBDA(:,:,:)**(XCXG_L-ZEXP-XDG_L))*MOMG(XALPHAG_L,XNUG_L,ZEXP+XDG_L)
        PRARE(:,:,:) = PRARE(:,:,:) + ZW(:,:,:)
      END WHERE
    end if
  ELSE
    ZFRAC_WATER = 0.14
    ZDMELT_FACT = ( (6.0*XAG_I)/(XPI*XRHOLW) )**(2.0)
    WHERE( PTEMP(:,:,:).GT.XTT )
      ZREFL_MELT_CONV(:,:,:) = ((1.0-ZFRAC_WATER)*ZEQICE+ZFRAC_WATER)*ZDMELT_FACT
    ELSEWHERE
      ZREFL_MELT_CONV(:,:,:) = ZEQICE*ZDMELT_FACT
    END WHERE
!
    ZEXP = 2.0*XBG_I
    WHERE( PRT(:,:,:,6).GT.XRTMIN_I(6) )
      ZLBDA(:,:,:) = XLBG_I*( PRHODREF(:,:,:)*PRT(:,:,:,6) )**XLBEXG_I
      ZW(:,:,:)    = ZREFL_MELT_CONV(:,:,:)*1.E18*XCCG_I*                &
                    (ZLBDA(:,:,:)**(XCXG_I-ZEXP))*MOMG(XALPHAG_I,XNUG_I,ZEXP)
      PVDOP(:,:,:) = PVDOP(:,:,:) +                                            &
                     ZREFL_MELT_CONV(:,:,:)*1.E18*XCCG_I*XCG_I*                    &
                     (ZLBDA(:,:,:)**(XCXG_I-ZEXP-XDG_I))*MOMG(XALPHAG_I,XNUG_I,ZEXP+XDG_I)
      PRARE(:,:,:) = PRARE(:,:,:) + ZW(:,:,:)
    END WHERE
  ENDIF
END IF
!
!*       7.     HAILSTONES
!               ----------
!²
IF (SIZE(PRT,4) >= 7) THEN
  IF (CCLOUD=='LIMA') THEN
    ZFRAC_WATER = 1.
    ZDMELT_FACT = ( (6.0*XAH_L)/(XPI*XRHOLW) )**(2.0)
    ZREFL_MELT_CONV(:,:,:) = ((1.0-ZFRAC_WATER)*ZEQICE+ZFRAC_WATER)*ZDMELT_FACT
!
    ZEXP = 2.0*XBH_L
    if (NMOM_H.GE.2) then
      WHERE( PRT(:,:,:,7).GT.XRTMIN_L(7) .AND. PCHT(:,:,:).GT.1.0E-3 ) 
        ZLBDA(:,:,:) =  XLBH_L**(XLBEXH_L)*(PRT(:,:,:,7)/PCHT(:,:,:))**(-XLBEXH_L)  
        ZW(:,:,:)    = ZREFL_MELT_CONV(:,:,:)     &
                    *1.E18*PRHODREF(:,:,:)*PCHT(:,:,:)*(ZLBDA(:,:,:)**(-ZEXP))*MOMG(XALPHAH_L,XNUH_L,ZEXP)   
        PVDOP(:,:,:) = PVDOP(:,:,:) +                                            &
                     ZREFL_MELT_CONV(:,:,:)            &
                       *1.E18*PRHODREF(:,:,:)*PCHT(:,:,:)*XCH_L*(ZLBDA(:,:,:)**(-ZEXP-XDH_L)) 
        PRARE(:,:,:) = PRARE(:,:,:) + ZW(:,:,:)
      END WHERE            
    else
      WHERE( PRT(:,:,:,7).GT.XRTMIN_L(7) )
        ZLBDA(:,:,:) = XLBH_L*( PRHODREF(:,:,:)*PRT(:,:,:,7) )**XLBEXH_L
        ZW(:,:,:)    = ZREFL_MELT_CONV(:,:,:)*1.E18*XCCH_L*                &
                     (ZLBDA(:,:,:)**(XCXH_L-ZEXP))*MOMG(XALPHAH_L,XNUH_L,ZEXP)
        PVDOP(:,:,:) = PVDOP(:,:,:) +                                            &
                     ZREFL_MELT_CONV(:,:,:)*1.E18*XCCH_L*XCH_L*                    &
                     (ZLBDA(:,:,:)**(XCXH_L-ZEXP-XDH_L))*MOMG(XALPHAH_L,XNUH_L,ZEXP+XDH_L)
        PRARE(:,:,:) = PRARE(:,:,:) + ZW(:,:,:)
      END WHERE
    end if
  ELSE
    ZFRAC_WATER = 1.
    ZDMELT_FACT = ( (6.0*XAH_I)/(XPI*XRHOLW) )**(2.0)
    ZREFL_MELT_CONV(:,:,:) = ((1.0-ZFRAC_WATER)*ZEQICE+ZFRAC_WATER)*ZDMELT_FACT
!
    ZEXP = 2.0*XBH_I
    WHERE( PRT(:,:,:,7).GT.XRTMIN_I(7) )
      ZLBDA(:,:,:) = XLBH_I*( PRHODREF(:,:,:)*PRT(:,:,:,7) )**XLBEXH_I
      ZW(:,:,:)    = ZREFL_MELT_CONV(:,:,:)*1.E18*XCCH_I*                &
                   (ZLBDA(:,:,:)**(XCXH_I-ZEXP))*MOMG(XALPHAH_I,XNUH_I,ZEXP)
      PVDOP(:,:,:) = PVDOP(:,:,:) +                                            &
                   ZREFL_MELT_CONV(:,:,:)*1.E18*XCCH_I*XCH_I*                    &
                   (ZLBDA(:,:,:)**(XCXH_I-ZEXP-XDH_I))*MOMG(XALPHAH_I,XNUH_I,ZEXP+XDH_I)
      PRARE(:,:,:) = PRARE(:,:,:) + ZW(:,:,:)
    END WHERE    
  END IF
END IF
!
!*       8.     UNIT CONVERSION
!               ---------------
!
IKB = 1 + JPVEXT
ZRHO00 = XP00/(XRD*XTHVREFZ(IKB))
WHERE( PRARE(:,:,:) >= 1.0 )
  PVDOP(:,:,:) = PVDOP(:,:,:)/PRARE(:,:,:)  ! Doppler speed normalization in m/s
  PVDOP(:,:,:) = PVDOP(:,:,:)*(ZRHO00/PRHODREF(:,:,:))**0.4
                                            ! air density correction
ELSEWHERE
  PVDOP(:,:,:) = 0.0
END WHERE
!
! MODIF FP FEB 2012
!WHERE( PRARE(:,:,:) > 0.0 )
WHERE( PRARE(:,:,:) > 1.E-3 )
! END MODIF
  PRARE(:,:,:) = 10.0*LOG10( PRARE(:,:,:) ) ! Z_equiv in dBZ
ELSEWHERE
  PRARE(:,:,:) = XUNDEF
END WHERE


!
!-------------------------------------------------------------------------------
!
CONTAINS
!
  FUNCTION MOMG (PALPHA,PNU,PP) RESULT (PMOMG)
!
! auxiliary routine used to compute the Pth moment order of the generalized
! gamma law
!
  USE MODI_GAMMA
!
  IMPLICIT NONE
!
  REAL     :: PALPHA ! first shape parameter of the dimensionnal distribution
  REAL     :: PNU    ! second shape parameter of the dimensionnal distribution
  REAL     :: PP     ! order of the moment
  REAL     :: PMOMG  ! result: moment of order ZP
!
!------------------------------------------------------------------------------
!
!
  PMOMG = GAMMA(PNU+PP/PALPHA)/GAMMA(PNU)
!
  END FUNCTION MOMG
!
!------------------------------------------------------------------------------
!
END SUBROUTINE RADAR_RAIN_ICE
