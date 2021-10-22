!     ######spl
     MODULE MODE_DUST_PSD
     USE PARKIND1, ONLY : JPRB
     USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!!   ########################
!!
!!    PURPOSE
!!    -------
!! MODULE DUST PSD (Particle Size Distribution)
!! Purpose: Contains subroutines to convert from transported variables (ppp)
!! to understandable aerosol variables, e.g. #/m3, kg/m3, sigma, R_{n}
!!
!!    AUTHOR
!!    ------
!!      Alf Grini (CNRM/GMEI)
!!
!!    MODIFICATIONS
!!    -------------
!!
!-------------------------------------------------------------------------------
!
USE MODD_CSTS_DUST         !Constants which are important for dust calculations
USE MODD_DUST              !Dust module which contains even more constants
!
IMPLICIT NONE
!
CONTAINS
!
!!   ############################################################
  SUBROUTINE PPP2DUST(             &
       PSVT                         & !I [ppp] input scalar variables (moment of distribution)
       , PRHODREF                   & !I [kg/m3] density of air       
       , PSIG3D                     & !O [-] standard deviation of aerosol distribution
       , PRG3D                      & !O [um] number median diameter of aerosol distribution
       , PN3D                       & !O [#/m3] number concentration of aerosols
       , PMASS3D                    & !O [kg/m3] mass concentration of aerosol
       , PM3D                       & !O aerosols moments 0, 3 and 6
       )
!!   ############################################################
!
!!
!!    PURPOSE
!!    -------
!!    Translate the three moments M0, M3 and M6 given in ppp into
!!    Values which can be understood more easily (R, sigma, N, M)
!! 
!!    CALLING STRUCTURE NOTE: OPTIONAL VARIABLES
!!    -------
!!    CALL PPP2AEROS(PSVT, PRHODREF, PSIG3D=SIGVAR,  &
!!       PRG3D=RVAR, PN3D=NVAR, PM3D=MASSVAR)
!!
!!    REFERENCE
!!    ---------
!!    none
!!
!!    AUTHOR
!!    ------
!!    Pierre TULET (LA)
!!
!!    MODIFICATIONS
!!    -------------
!!    2005 Alf Grini (CNRM)
!!    2006 Jean-Pierre Chaboureau (LA)
!!
!!    EXTERNAL
!!    --------
!!    None
!!
    IMPLICIT NONE
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
!*      0.1    declarations of arguments
!
REAL,       DIMENSION(:,:,:,:),  INTENT(INOUT)  :: PSVT      !I [ppp] first moment
REAL,       DIMENSION(:,:,:),    INTENT(IN)     :: PRHODREF !I [kg/m3] density of air

REAL,       DIMENSION(:,:,:,:),  OPTIONAL, INTENT(OUT)     :: PSIG3D   !O [-] standard deviation
REAL,       DIMENSION(:,:,:,:),  OPTIONAL, INTENT(OUT)     :: PRG3D    !O [um] number median diameter
REAL,       DIMENSION(:,:,:,:),  OPTIONAL, INTENT(OUT)     :: PN3D     !O [#/m3] number concentration
REAL,       DIMENSION(:,:,:,:),  OPTIONAL, INTENT(OUT)     :: PMASS3D  !O [kg_{aer}/m3] mass concentration
REAL,       DIMENSION(:,:,:,:),  OPTIONAL, INTENT(OUT)     :: PM3D     !O aerosols moments 
!
!
!*      0.2    declarations local variables
!
REAL                                  :: ZRHOI               ! [kg/m3] density of aerosol
REAL                                  :: ZMI                 ! [kg/mol] molar weight of aerosol
REAL                                  :: ZRGMIN              ! [um] minimum radius accepted
REAL                                  :: ZSIGMIN             ! minimum standard deviation accepted
REAL,DIMENSION(:,:,:,:), ALLOCATABLE  :: ZM                  ! [aerosol units] local array which goes to output later
REAL,DIMENSION(:,:,:,:), ALLOCATABLE  :: ZSV                 ! [dusts moment concentration]
REAL,DIMENSION(:,:,:),   ALLOCATABLE  :: ZSIGMA              ! [-] standard deviation
REAL,DIMENSION(:,:,:),   ALLOCATABLE  :: ZRG                 ! [um] number median diameter
REAL,DIMENSION(:),       ALLOCATABLE  :: ZMMIN               ! [aerosol units] minimum values for N, sigma, M
INTEGER,DIMENSION(:),    ALLOCATABLE  :: NM0                 ! [idx] index for Mode 0 in passed variables
INTEGER,DIMENSION(:),    ALLOCATABLE  :: NM3                 ! [idx] indexes for Mode 3 in passed variables
INTEGER,DIMENSION(:),    ALLOCATABLE  :: NM6                 ! [idx] indexes for Mode 6 in passed variables
REAL,DIMENSION(:),       ALLOCATABLE  :: ZINIRADIUS          ! initial mean radius
INTEGER                               :: JN,IMODEIDX,JJ      ! [idx] loop counters
!
!-------------------------------------------------------------------------------
!
!        1.1    initialisation 
!
!Calculations here are for one mode only
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_DUST_PSD:PPP2DUST',0,ZHOOK_HANDLE)
ALLOCATE (NM0(NMODE_DST))
ALLOCATE (NM3(NMODE_DST))
ALLOCATE (NM6(NMODE_DST))
ALLOCATE (ZM(SIZE(PSVT,1), SIZE(PSVT,2), SIZE(PSVT,3), NMODE_DST*3))
ALLOCATE (ZMMIN(NMODE_DST*3))
ALLOCATE (ZSIGMA(SIZE(PSVT,1), SIZE(PSVT,2), SIZE(PSVT,3)))
ALLOCATE (ZRG(SIZE(PSVT,1), SIZE(PSVT,2), SIZE(PSVT,3)))
ALLOCATE (ZSV(SIZE(PSVT,1), SIZE(PSVT,2), SIZE(PSVT,3), SIZE(PSVT,4)))
ALLOCATE (ZINIRADIUS(NMODE_DST))
    
ZSV(:,:,:,:) = MAX(PSVT(:,:,:,:), 1.E-80)

DO JN=1,NMODE_DST
  IMODEIDX = JPDUSTORDER(JN)
  !Calculations here are for one mode only
  IF (CRGUNITD=="MASS") THEN
    ZINIRADIUS(JN) = XINIRADIUS(IMODEIDX) * EXP(-3.*(LOG(XINISIG(IMODEIDX)))**2)
  ELSE
    ZINIRADIUS(JN) = XINIRADIUS(IMODEIDX)
  END IF

  !Set counter for number, M3 and M6
  NM0(JN) = 1+(JN-1)*3
  NM3(JN) = 2+(JN-1)*3
  NM6(JN) = 3+(JN-1)*3
  !Get minimum values possible
  ZMMIN(NM0(JN)) = XN0MIN(IMODEIDX)
  ZRGMIN         = ZINIRADIUS(JN)
  IF (LVARSIG) THEN
    ZSIGMIN = XSIGMIN
  ELSE
    ZSIGMIN = XINISIG(IMODEIDX)
  ENDIF
  ZMMIN(NM3(JN)) = ZMMIN(NM0(JN)) * (ZRGMIN**3)*EXP(4.5 * LOG(ZSIGMIN)**2) 
  ZMMIN(NM6(JN)) = ZMMIN(NM0(JN)) * (ZRGMIN**6)*EXP(18. * LOG(ZSIGMIN)**2)
END DO
!
!Set density of aerosol, here dust (kg/m3)
ZRHOI = XDENSITY_DUST
!Set molecular weight of dust !NOTE THAT THIS IS NOW IN KG
ZMI   = XMOLARWEIGHT_DUST
!
!
DO JN=1,NMODE_DST
  !
  IF (LVARSIG) THEN ! give M6 (case of variable standard deviation)
  !
  !Get number concentration (#/molec_{air}==>#/m3)
    ZM(:,:,:,NM0(JN))=                         &
         ZSV(:,:,:,1+(JN-1)*3)                 & !#/molec_{air}
         * XAVOGADRO                           & !==>#/mole
         / XMD                                 & !==>#/kg_{air}
         * PRHODREF(:,:,:)                       !==>#/m3
  ! 
  !calculate moment 3 from total aerosol mass (molec_{aer}/molec_{air} ==>
    ZM(:,:,:,NM3(JN)) =                        &
         ZSV(:,:,:,2+(JN-1)*3)                 & !molec_{aer}/molec_{aer}
         * (ZMI/XMD)                           & !==>kg_{aer}/kg_{aer}
         * PRHODREF(:,:,:)                     & !==>kg_{aer}/m3_{air}
         * (1.d0/ZRHOI)                        & !==>m3_{aer}/m3_{air}
         * XM3TOUM3                            & !==>um3_{aer}/m3_{air}
         / (XPI * 4./3.)                         !==>um3_{aer}/m3_{air} (volume ==> 3rd moment)
     !Limit mass concentration to minimum value
     ZM(:,:,:,NM3(JN)) = MAX(ZM(:,:,:,NM3(JN)), ZMMIN(NM3(JN)))
  ! 
    ZM(:,:,:,NM6(JN)) = ZSV(:,:,:,3+(JN-1)*3)  & !um6/molec_{air}*(cm3/m3)
         * 1.d-6                               & !==> um6/molec_{air}
         * XAVOGADRO                           & !==> um6/mole_{air}
         / XMD                                 & !==> um6/kg_{air}
         * PRHODREF(:,:,:)                       !==> um6/m3_{air}
     !Limit m6 concentration to minimum value
     ZM(:,:,:,NM6(JN)) =  MAX(ZM(:,:,:,NM6(JN)), ZMMIN(NM6(JN)))
  !
  !Get sigma (only if sigma is allowed to vary)
    !Get intermediate values for sigma M3^2/(M0*M6) (ORILAM paper, eqn 8)
    ZSIGMA(:,:,:)=ZM(:,:,:,NM3(JN))**2/(ZM(:,:,:,NM0(JN))*ZM(:,:,:,NM6(JN)))
    !Limit the intermediate value, can not be larger than 1
    ZSIGMA(:,:,:)=MIN(1-1E-10,ZSIGMA(:,:,:))
    !Limit the value for intermediate, can not be smaller than 0
    ZSIGMA(:,:,:)=MAX(1E-10,ZSIGMA(:,:,:))
    !Calculate log(sigma)
    ZSIGMA(:,:,:)= LOG(ZSIGMA(:,:,:))
    !Finally get the real sigma the negative sign is because of 
    !The way the equation is written (M3^2/(M0*M6)) instead of (M0*M6)/M3^3
    ZSIGMA(:,:,:)= EXP(1./3.*SQRT(-ZSIGMA(:,:,:)))
    !Limit the value to reasonable ones
    ZSIGMA(:,:,:) =  MAX( XSIGMIN, MIN( XSIGMAX, ZSIGMA(:,:,:) ) )

  !
    !Put back M6 so that it fits the sigma which is possibly modified above
    !The following makes M6 consistent with N, R, SIGMA
    ZM(:,:,:,NM6(JN)) = ZM(:,:,:,NM0(JN)) &
         * ( (ZM(:,:,:,NM3(JN))/ZM(:,:,:,NM0(JN)))**(1./3.) &
         * exp(-(3./2.)*log(ZSIGMA(:,:,:))**2))**6 &
         * exp(18.*log(ZSIGMA(:,:,:))**2)

  ELSE ! compute M6 from M0, M3 and SIGMA
    ! 
    ZSIGMA(:,:,:) = XINISIG(JPDUSTORDER(JN))
 
    IF (LRGFIX_DST) THEN
    !calculate moment 3 from total aerosol mass (molec_{aer}/molec_{air} ==>
    ZM(:,:,:,NM3(JN)) =                        &
         ZSV(:,:,:,JN)                         & !molec_{aer}/molec_{aer}
         * (ZMI/XMD)                           & !==>kg_{aer}/kg_{aer}
         * PRHODREF(:,:,:)                     & !==>kg_{aer}/m3_{air}
         * (1.d0/ZRHOI)                        & !==>m3_{aer}/m3_{air}
         * XM3TOUM3                            & !==>um3_{aer}/m3_{air}
         / (XPI * 4./3.)                         !==>um3_{aer}/m3_{air} (volume ==> 3rd moment)

     ZM(:,:,:,NM3(JN)) = MAX(ZM(:,:,:,NM3(JN)), ZMMIN(NM3(JN)))

     ZM(:,:,:,NM0(JN))=  ZM(:,:,:,NM3(JN))/&
       ((ZINIRADIUS(JN)**3)*EXP(4.5 * LOG(XINISIG(JPDUSTORDER(JN)))**2))

    ELSE
    !calculate moment 3 from total aerosol mass (molec_{aer}/molec_{air} ==>
    ZM(:,:,:,NM3(JN)) =                        &
         ZSV(:,:,:,2+(JN-1)*2)                 & !molec_{aer}/molec_{aer}
         * (ZMI/XMD)                           & !==>kg_{aer}/kg_{aer}
         * PRHODREF(:,:,:)                     & !==>kg_{aer}/m3_{air}
         * (1.d0/ZRHOI)                        & !==>m3_{aer}/m3_{air}
         * XM3TOUM3                            & !==>um3_{aer}/m3_{air}
         / (XPI * 4./3.)                         !==>um3_{aer}/m3_{air} (volume ==> 3rd moment)

    !Get number concentration (#/molec_{air}==>#/m3)
     ZM(:,:,:,NM0(JN))=                         &
         ZSV(:,:,:,1+(JN-1)*2)                 & !#/molec_{air}
         * XAVOGADRO                           & !==>#/mole
         / XMD                                 & !==>#/kg_{air}
         * PRHODREF(:,:,:)                       !==>#/m3

    ! Limit concentration to minimum values
    WHERE ((ZM(:,:,:,NM0(JN)) < ZMMIN(NM0(JN)) ).OR. &
           (ZM(:,:,:,NM3(JN)) < ZMMIN(NM3(JN)) )) 
       ZM(:,:,:,NM0(JN)) = ZMMIN(NM0(JN))
       ZM(:,:,:,NM3(JN)) = ZMMIN(NM3(JN))
       PSVT(:,:,:,1+(JN-1)*2) = ZM(:,:,:,NM0(JN)) * XMD / &
       (XAVOGADRO * PRHODREF(:,:,:) )
       PSVT(:,:,:,2+(JN-1)*2) = ZM(:,:,:,NM3(JN)) * XMD * XPI * 4./3. * ZRHOI  / &
                              (ZMI*PRHODREF(:,:,:)*XM3TOUM3)
    ENDWHERE
    END IF

    ! 

    ZM(:,:,:,NM6(JN)) = ZM(:,:,:,NM0(JN))                   &
         * ( (ZM(:,:,:,NM3(JN))/ZM(:,:,:,NM0(JN)))**(1./3.) &
         * exp(-(3./2.)*log(ZSIGMA(:,:,:))**2))**6          &
         * exp(18.*log(ZSIGMA(:,:,:))**2)
    
  !
  END IF
  !   
  !Get number median radius (eqn. 7 in Orilam manuscript)
  ZRG(:,:,:)=    &
          (      &
          ZM(:,:,:,NM3(JN))*ZM(:,:,:,NM3(JN))*ZM(:,:,:,NM3(JN))*ZM(:,:,:,NM3(JN))    &
          /(ZM(:,:,:,NM6(JN))*ZM(:,:,:,NM0(JN))*ZM(:,:,:,NM0(JN))*ZM(:,:,:,NM0(JN))) &
          )                                                                          &
          ** XSIXTH 
  !ZRG(:,:,:)=MIN(ZRG(:,:,:),ZINIRADIUS(JN))
  !Give the sigma-values to the passed array
  IF(PRESENT(PSIG3D)) PSIG3D(:,:,:,JN) = ZSIGMA(:,:,:)
  !
  !Set the number concentrations in the passed array
  IF(PRESENT(PN3D)) PN3D(:,:,:,JN) = ZM(:,:,:,NM0(JN))
  !
  !Get the number median radius
  IF(PRESENT(PRG3D)) PRG3D(:,:,:,JN)= ZRG(:,:,:)
  !
  IF(PRESENT(PMASS3D))THEN
       PMASS3D(:,:,:,JN)=     &
            ZM(:,:,:,NM0(JN)) &    !#/m^3_{air}
            * XPI*4./3.       &    
            * ZRHOI           &    !==>kg/m^3_{aeros}/m^3_{air}
            * ZRG(:,:,:) * ZRG(:,:,:) * ZRG(:,:,:) &
            * XUM3TOM3        &    !==>kg/m^3_{air}
            * exp(4.5*log(ZSIGMA(:,:,:))*log(ZSIGMA(:,:,:)))
  ENDIF
!
END DO  !Loop on modes
!
IF(PRESENT(PM3D)) PM3D(:,:,:,:) = ZM(:,:,:,:)
!
DEALLOCATE(ZINIRADIUS)
DEALLOCATE(ZSV)
DEALLOCATE(ZRG)
DEALLOCATE(ZSIGMA)
DEALLOCATE(ZMMIN)
DEALLOCATE(ZM)
DEALLOCATE(NM6)
DEALLOCATE(NM3)
DEALLOCATE(NM0)
!
!
IF (LHOOK) CALL DR_HOOK('MODE_DUST_PSD:PPP2DUST',1,ZHOOK_HANDLE)
END SUBROUTINE PPP2DUST

!!   ############################################################
  SUBROUTINE DUST2PPP(             &
       PSVT                         & !IO [ppp] input scalar variables (moment of distribution)
       , PRHODREF                   & !I [kg/m3] density of air       
       , PSIG3D                     & !I [-] standard deviation of aerosol distribution
       , PRG3D                      & !I [um] number median diameter of aerosol distribution
       )
!!   ############################################################
!
!!
!!    PURPOSE
!!    -------
!!    Translate the dust Mass, RG and SIGMA in the  three moments M0, M3 and M6 given in ppp 
!! 
!!    CALLING STRUCTURE NOTE: OPTIONAL VARIABLES
!!    -------
!!    CALL PPP2AEROS(PSVT, PRHODREF, PSIG3D=SIGVAR,  &
!!       PRG3D=RVAR, PN3D=NVAR)
!!
!!    REFERENCE
!!    ---------
!!    none
!!
!!    AUTHOR
!!    ------
!!    Pierre TULET (LA)
!!
!!    MODIFICATIONS
!!    -------------
!!    Alf Grini (CNRM)
!!
!!    EXTERNAL
!!    --------
!!    None
!!
    IMPLICIT NONE
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
!*      0.1    declarations of arguments
!
    !INPUT
    REAL,       DIMENSION(:,:,:),    INTENT(IN)     :: PRHODREF !I [kg/m3] density of air
    REAL,       DIMENSION(:,:,:,:),  INTENT(IN)     :: PSIG3D   !O [-] standard deviation
    REAL,       DIMENSION(:,:,:,:),  INTENT(IN)     :: PRG3D    !O [um] number median diameter

    !OUTPUT
    REAL,       DIMENSION(:,:,:,:),  INTENT(INOUT)  :: PSVT  !IO [#/molec_{air}] first moment
                                                                !IO [molec_{aer}/molec_{air} 3rd moment
                                                                !IO [um6/molec_{air}*(cm3/m3)] 6th moment
!
!
!*      0.2    declarations local variables
!
    REAL                                  :: ZRHOI               ! [kg/m3] density of aerosol
    REAL                                  :: ZMI                 ! [kg/mol] molar weight of aerosol
    REAL                                  :: ZRGMIN              ! [um] minimum radius accepted
    REAL                                  :: ZSIGMIN             ! minimum standard deviation accepted
    REAL,DIMENSION(:,:,:,:), ALLOCATABLE  :: ZM                  ! [aerosol units] local array which goes to output later
    REAL,DIMENSION(:,:,:),   ALLOCATABLE  :: ZSIGMA              ! aersol standard deviation
    REAL,DIMENSION(:),       ALLOCATABLE  :: ZMMIN               ! [aerosol units] minimum values for N, sigma, M
    REAL,DIMENSION(:),       ALLOCATABLE  :: ZINIRADIUS          ! initial mean radius
    INTEGER,DIMENSION(:),    ALLOCATABLE  :: NM0                 ! [idx] index for Mode 0 in passed variables
    INTEGER,DIMENSION(:),    ALLOCATABLE  :: NM3                 ! [idx] indexes for Mode 3 in passed variables
    INTEGER,DIMENSION(:),    ALLOCATABLE  :: NM6                 ! [idx] indexes for Mode 6 in passed variables
    INTEGER                               :: JJ, JN              ! [idx] loop counters
    INTEGER                               :: IMODEIDX
!
!-------------------------------------------------------------------------------
!
!        1.1    initialisation 


    REAL(KIND=JPRB) :: ZHOOK_HANDLE
    IF (LHOOK) CALL DR_HOOK('MODE_DUST_PSD:DUST2PPP',0,ZHOOK_HANDLE)
    ALLOCATE (NM0(NMODE_DST))
    ALLOCATE (NM3(NMODE_DST))
    ALLOCATE (NM6(NMODE_DST))
    ALLOCATE (ZINIRADIUS(NMODE_DST))
    ALLOCATE (ZMMIN(NMODE_DST*3))
    ALLOCATE (ZM(SIZE(PSVT,1), SIZE(PSVT,2), SIZE(PSVT,3), NMODE_DST*3))
    ALLOCATE (ZSIGMA(SIZE(PSVT,1), SIZE(PSVT,2), SIZE(PSVT,3)))

    !Set density of aerosol, here dust (kg/m3)
    ZRHOI = XDENSITY_DUST
    !Set molecular weight of dust !NOTE THAT THIS IS NOW IN KG
    ZMI   = XMOLARWEIGHT_DUST
!

    ! PSVT need to be positive
    PSVT(:,:,:,:) = MAX(PSVT(:,:,:,:), 1E-80)
    
    DO JN=1,NMODE_DST
    IMODEIDX = JPDUSTORDER(JN)
    !Calculations here are for one mode only
    IF (CRGUNITD=="MASS") THEN
      ZINIRADIUS(JN) = XINIRADIUS(IMODEIDX) * EXP(-3.*(LOG(XINISIG(IMODEIDX)))**2)
    ELSE
      ZINIRADIUS(JN) = XINIRADIUS(IMODEIDX)
    END IF

    !Set counter for number, M3 and M6
       NM0(JN) = 1+(JN-1)*3
       NM3(JN) = 2+(JN-1)*3
       NM6(JN) = 3+(JN-1)*3

    !Get minimum values possible
    ZMMIN(NM0(JN)) = XN0MIN(IMODEIDX)
    ZRGMIN     =  ZINIRADIUS(JN)
    IF (LVARSIG) THEN
      ZSIGMIN = XSIGMIN
    ELSE
      ZSIGMIN = XINISIG(IMODEIDX)
    ENDIF
    ZMMIN(NM3(JN)) = ZMMIN(NM0(JN)) * (ZRGMIN**3)*EXP(4.5 * LOG(ZSIGMIN)**2) 
    ZMMIN(NM6(JN)) = ZMMIN(NM0(JN)) * (ZRGMIN**6)*EXP(18. * LOG(ZSIGMIN)**2)
    END DO

    !Set density of aerosol, here dust (kg/m3)
    ZRHOI = XDENSITY_DUST
    !Set molecular weight of dust !NOTE THAT THIS IS NOW IN KG
    ZMI   = XMOLARWEIGHT_DUST
!
    DO JN=1,NMODE_DST
     !calculate moment 3 from total aerosol mass (molec_{aer}/molec_{air} ==>
     IF (LVARSIG) THEN
     ZM(:,:,:,NM3(JN)) =                        &
          PSVT(:,:,:,2+(JN-1)*3)                & !molec_{aer}/molec_{aer}
          * (ZMI/XMD)                           & !==>kg_{aer}/kg_{aer}
          * PRHODREF(:,:,:)                     & !==>kg_{aer}/m3_{air}
          * (1.d0/ZRHOI)                        & !==>m3_{aer}/m3_{air}
          * XM3TOUM3                            & !==>um3_{aer}/m3_{air}
          / (XPI * 4./3.)                         !==>um3_{aer}/m3_{air} (volume ==> 3rd moment)
    ELSE
       IF ((LRGFIX_DST)) THEN
          ZM(:,:,:,NM3(JN)) =                   &
          PSVT(:,:,:,JN)                        & !molec_{aer}/molec_{aer}
          * (ZMI/XMD)                           & !==>kg_{aer}/kg_{aer}
          * PRHODREF(:,:,:)                     & !==>kg_{aer}/m3_{air}
          * (1.d0/ZRHOI)                        & !==>m3_{aer}/m3_{air}
          * XM3TOUM3                            & !==>um3_{aer}/m3_{air}
          / (XPI * 4./3.)                         !==>um3_{aer}/m3_{air} (volume ==> 3rd moment)
         ZM(:,:,:,NM3(JN)) = MAX(ZM(:,:,:,NM3(JN)), ZMMIN(NM3(JN)))
      ELSE
          ZM(:,:,:,NM3(JN)) =                   &
          PSVT(:,:,:,2+(JN-1)*2)                & !molec_{aer}/molec_{aer}
          * (ZMI/XMD)                           & !==>kg_{aer}/kg_{aer}
          * PRHODREF(:,:,:)                     & !==>kg_{aer}/m3_{air}
          * (1.d0/ZRHOI)                        & !==>m3_{aer}/m3_{air}
          * XM3TOUM3                            & !==>um3_{aer}/m3_{air}
          / (XPI * 4./3.)                         !==>um3_{aer}/m3_{air} (volume ==> 3rd moment)
      END IF
    END IF
! calculate moment 0 from dispersion and mean radius
     ZM(:,:,:,NM0(JN))=  ZM(:,:,:,NM3(JN))/&
       ((PRG3D(:,:,:,JN)**3)*EXP(4.5 * LOG(PSIG3D(:,:,:,JN))**2))


! calculate moment 6 from dispersion and mean radius
     ZM(:,:,:,NM6(JN)) = ZM(:,:,:,NM0(JN)) * (PRG3D(:,:,:,JN)**6) * &
               EXP(18 *(LOG(PSIG3D(:,:,:,JN)))**2)

     IF (LVARSIG) THEN
     WHERE ((ZM(:,:,:,NM0(JN)) .LT. ZMMIN(NM0(JN))).OR.&
            (ZM(:,:,:,NM3(JN)) .LT. ZMMIN(NM3(JN))).OR.&
            (ZM(:,:,:,NM6(JN)) .LT. ZMMIN(NM6(JN))))
     ZM(:,:,:,NM0(JN)) = ZMMIN(NM0(JN))
     ZM(:,:,:,NM3(JN)) = ZMMIN(NM3(JN))
     ZM(:,:,:,NM6(JN)) = ZMMIN(NM6(JN))
     END WHERE

     ELSE IF (.NOT.(LRGFIX_DST)) THEN
     
     WHERE ((ZM(:,:,:,NM0(JN)) .LT. ZMMIN(NM0(JN))).OR.&
            (ZM(:,:,:,NM3(JN)) .LT. ZMMIN(NM3(JN))))
     ZM(:,:,:,NM0(JN)) = ZMMIN(NM0(JN))
     ZM(:,:,:,NM3(JN)) = ZMMIN(NM3(JN))
     END WHERE
     ENDIF

     
     ! return to concentration #/m3 =>  (#/molec_{air}
     IF (LVARSIG) THEN
     PSVT(:,:,:,1+(JN-1)*3) = ZM(:,:,:,NM0(JN)) * XMD / &
                              (XAVOGADRO*PRHODREF(:,:,:))

     PSVT(:,:,:,2+(JN-1)*3) = ZM(:,:,:,NM3(JN)) * XMD  * XPI * 4./3 * ZRHOI / &
                              (ZMI*PRHODREF(:,:,:)*XM3TOUM3)

     PSVT(:,:,:,3+(JN-1)*3) = ZM(:,:,:,NM6(JN)) * XMD  / &
                              ( XAVOGADRO*PRHODREF(:,:,:) * 1.d-6) 
     ELSE IF (LRGFIX_DST) THEN
      PSVT(:,:,:,JN)        = ZM(:,:,:,NM3(JN)) * XMD * XPI * 4./3. * ZRHOI  / &
                              (ZMI*PRHODREF(:,:,:)*XM3TOUM3)
     ELSE
     PSVT(:,:,:,1+(JN-1)*2) = ZM(:,:,:,NM0(JN)) * XMD / &
                              (XAVOGADRO*PRHODREF(:,:,:))

     PSVT(:,:,:,2+(JN-1)*2) = ZM(:,:,:,NM3(JN)) * XMD * XPI * 4./3. * ZRHOI  / &
                              (ZMI*PRHODREF(:,:,:)*XM3TOUM3)
     END IF

!
 END DO  !Loop on modes

DEALLOCATE(ZINIRADIUS)
DEALLOCATE(ZMMIN)
DEALLOCATE(ZSIGMA)
DEALLOCATE(ZM)
DEALLOCATE(NM6)
DEALLOCATE(NM3)
DEALLOCATE(NM0)
!
IF (LHOOK) CALL DR_HOOK('MODE_DUST_PSD:DUST2PPP',1,ZHOOK_HANDLE)
END SUBROUTINE DUST2PPP
!
END MODULE MODE_DUST_PSD
