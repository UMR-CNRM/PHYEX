!     ######spl
      SUBROUTINE DUSTLFI_n(PSV, PRHODREF)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!     ############################################################
!
!!    PURPOSE
!!    -------
!!    Realise l'équilibre des moments à partir du sigma et du diametre moyen
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
!!    none
!!
!!    EXTERNAL
!!    --------
!!    None
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_DUST
USE MODD_NSV
USE MODD_CSTS_DUST
USE MODE_DUST_PSD
! 
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL,   DIMENSION(:,:,:,:),    INTENT(INOUT) :: PSV
REAL,   DIMENSION(:,:,:),      INTENT(IN) :: PRHODREF
!
!
!*      0.2    declarations local variables
!
REAL   :: ZDEN2MOL, ZRHOI, ZMI, ZFAC, ZRGMIN
REAL,DIMENSION(:,:,:,:), ALLOCATABLE  :: ZCTOTA
REAL,DIMENSION(:,:,:,:), ALLOCATABLE  :: ZM
REAL,DIMENSION(:,:,:),   ALLOCATABLE  :: ZSIGMA
REAL,DIMENSION(:,:,:,:),   ALLOCATABLE  :: ZSIG, ZRG, ZN
INTEGER,DIMENSION(:),    ALLOCATABLE  :: IM0, IM3, IM6
REAL,DIMENSION(:),       ALLOCATABLE  :: ZMMIN
REAL,DIMENSION(:),       ALLOCATABLE  :: ZINIRADIUS, ZINISIGMA
INTEGER :: IKU
INTEGER :: JJ, JN, JK  ! loop counter
INTEGER :: IMODEIDX  ! index mode
REAL, PARAMETER  :: ZMMR_DUST=5.d-11!kg_{dust}/kg_{air}
REAL    :: ZMMR_DUSTN
!
!-------------------------------------------------------------------------------
!
!*       1.     TRANSFER FROM GAS TO AEROSOL MODULE
!               -----------------------------------
!
!        1.1    initialisation 
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('DUSTLFI_N',0,ZHOOK_HANDLE)
IKU=SIZE(PSV,3)
!
ALLOCATE (IM0(NMODE_DST))
ALLOCATE (IM3(NMODE_DST))
ALLOCATE (IM6(NMODE_DST))
ALLOCATE (ZCTOTA(SIZE(PSV,1), SIZE(PSV,2), SIZE(PSV,3), NMODE_DST))
ALLOCATE (ZM(SIZE(PSV,1), SIZE(PSV,2), SIZE(PSV,3), NMODE_DST*3))
ALLOCATE (ZSIGMA(SIZE(PSV,1), SIZE(PSV,2), SIZE(PSV,3)))
ALLOCATE (ZSIG(SIZE(PSV,1), SIZE(PSV,2), SIZE(PSV,3),NMODE_DST))
ALLOCATE (ZRG(SIZE(PSV,1), SIZE(PSV,2), SIZE(PSV,3),NMODE_DST))
ALLOCATE (ZN(SIZE(PSV,1), SIZE(PSV,2), SIZE(PSV,3),NMODE_DST))
ALLOCATE (ZINIRADIUS(NMODE_DST))
ALLOCATE (ZINISIGMA(NMODE_DST))
ALLOCATE (ZMMIN(NMODE_DST*3))
!
!
DO JN = 1, NMODE_DST
  IM0(JN) = 1+(JN-1)*3
  IM3(JN) = 2+(JN-1)*3
  IM6(JN) = 3+(JN-1)*3
  !
  !Get the dust mode we are talking about, MODE 2 is treated first, then mode 3, then 1
  !This index is only needed to get the right radius out of the XINIRADIUS array and the
  !right XINISIG out of the XINISIG-array
  IMODEIDX = JPDUSTORDER(JN)
  !
  !Convert initial mass median radius to number median radius
  IF (CRGUNITD=="MASS") THEN
    ZINIRADIUS(JN) = XINIRADIUS(IMODEIDX) * EXP(-3.*(LOG(XINISIG(IMODEIDX)))**2)
  ELSE
    ZINIRADIUS(JN) = XINIRADIUS(IMODEIDX)
  END IF
  ZINISIGMA(JN)  = XINISIG(IMODEIDX)
  !
  ZMMIN(IM0(JN)) = XN0MIN(IMODEIDX)
! ZRGMIN   = XCOEFRADMIN * ZINIRADIUS(JN)
  ZRGMIN   = ZINIRADIUS(JN)
  ZMMIN(IM3(JN)) = XN0MIN(IMODEIDX) * (ZRGMIN**3)*EXP(4.5 * LOG(ZINISIGMA(JN))**2) 
  ZMMIN(IM6(JN)) = XN0MIN(IMODEIDX) * (ZRGMIN**6)*EXP(18. * LOG(ZINISIGMA(JN))**2)
ENDDO
!
!
ZRHOI = XDENSITY_DUST !1.8e3 !++changed alfgr
!ZMI   = XMOLARWEIGHT_DUST*1.D3 !100.  !++changed alfgr
ZMI   = XMOLARWEIGHT_DUST
ZDEN2MOL = 1E-6 * XAVOGADRO / XMD
ZFAC=(4./3.)*XPI*ZRHOI*1.e-9

! conversion into mol.cm-3
DO JJ=1,SIZE(PSV,4)
  PSV(:,:,:,JJ) =  PSV(:,:,:,JJ) * ZDEN2MOL * PRHODREF(:,:,:)
END DO
!
DO JN=1,NMODE_DST

!*       1.1    calculate moment 0 from XN0MIN
!
!  IF (JN == 1) ZMMR_DUSTN = 0.09 *  ZMMR_DUST
!  IF (JN == 2) ZMMR_DUSTN = 0.45 *  ZMMR_DUST
!  IF (JN == 3) ZMMR_DUSTN = 0.46 *  ZMMR_DUST
  IF (JN == 1) ZMMR_DUSTN = 0.40 *  ZMMR_DUST
  IF (JN == 2) ZMMR_DUSTN = 0.45 *  ZMMR_DUST
  IF (JN == 3) ZMMR_DUSTN = 0.15 *  ZMMR_DUST
  DO JK=1,IKU
    ZM(:,:,JK,IM0(JN)) =                   &
                ZMMR_DUSTN                 &![kg_{dust}/kg_{air}  
                *PRHODREF(:,:,JK)          &![kg_{air}/m3_{air}]==> kg/m3
                /XDENSITY_DUST             &![kg__{dust}/m3_{dust}==>m3_{dust}/m3{air}
                *(6.d0/XPI)                 &
                /(2.d0*ZINIRADIUS(JN)*1.d-6)**3 &![particle/m_dust^{-3}]==> particle/m3
                *EXP(-4.5*(LOG(XINISIG(JPDUSTORDER(JN))))**2) !Take into account distribution
  END DO
  ZM(:,:,:,IM0(JN)) = MAX(ZMMIN(IM0(JN)), ZM(:,:,:,IM0(JN)))
!
!*       1.2    calculate moment 3 from m0,  RG and SIG 
!
  ZM(:,:,:,IM3(JN)) = ZM(:,:,:,IM0(JN)) * &
            (ZINIRADIUS(JN)**3)*EXP(4.5 * LOG(ZINISIGMA(JN))**2) 
  ZM(:,:,:,IM3(JN)) = MAX(ZMMIN(IM3(JN)), ZM(:,:,:,IM3(JN)))
!
!*       1.3    calculate moment 6 from m0,  RG and SIG 
!
  ZM(:,:,:,IM6(JN))= ZM(:,:,:,IM0(JN)) * ((ZINIRADIUS(JN)**6)*&
                        EXP(18. * (LOG(ZINISIGMA(JN)))**2))
  ZM(:,:,:,IM6(JN)) = MAX(ZMMIN(IM6(JN)), ZM(:,:,:,IM6(JN)))
!
!*       1.4    output concentration
!
  PSV(:,:,:,1+(JN-1)*3) = ZM(:,:,:,IM0(JN)) * XMD / (XAVOGADRO*PRHODREF(:,:,:))
  PSV(:,:,:,2+(JN-1)*3) = ZM(:,:,:,IM3(JN)) * XMD*XPI * 4./3. * ZRHOI  / &
                           (ZMI*PRHODREF(:,:,:)*XM3TOUM3)

  PSV(:,:,:,3+(JN-1)*3) = ZM(:,:,:,IM6(JN)) *  XMD / (XAVOGADRO*PRHODREF(:,:,:)*1.d-6)
!
END DO
CALL PPP2DUST(PSV, PRHODREF, PSIG3D=ZSIG, PRG3D=ZRG, PN3D=ZN) 
!
DEALLOCATE(ZMMIN)
DEALLOCATE(ZINISIGMA)
DEALLOCATE(ZINIRADIUS)
DEALLOCATE(ZSIGMA)
DEALLOCATE(ZM)
DEALLOCATE(ZCTOTA)
DEALLOCATE(IM6)
DEALLOCATE(IM3)
DEALLOCATE(IM0)
!
!
IF (LHOOK) CALL DR_HOOK('DUSTLFI_N',1,ZHOOK_HANDLE)
END SUBROUTINE DUSTLFI_n
