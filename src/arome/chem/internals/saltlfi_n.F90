!     ######spl
      SUBROUTINE SALTLFI_n(PSV, PRHODREF)
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
USE MODD_SALT
USE MODD_NSV
USE MODD_CSTS_SALT
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
INTEGER,DIMENSION(:),    ALLOCATABLE  :: IM0, IM3, IM6
REAL,DIMENSION(:),       ALLOCATABLE  :: ZMMIN
REAL,DIMENSION(:),       ALLOCATABLE  :: ZINIRADIUS, ZINISIGMA
INTEGER :: IKU
INTEGER :: JJ, JN, JK  ! loop counter
INTEGER :: IMODEIDX  ! index mode
REAL, PARAMETER  :: ZMMR_SALT=1.d-7 !kg_{sea salt}/kg_{air}
REAL    :: ZMMR_SALTN
!
!-------------------------------------------------------------------------------
!
!*       1.     TRANSFER FROM GAS TO AEROSOL MODULE
!               -----------------------------------
!
!        1.1    initialisation 
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('SALTLFI_N',0,ZHOOK_HANDLE)
IKU=SIZE(PSV,3)
!
ALLOCATE (IM0(NMODE_SLT))
ALLOCATE (IM3(NMODE_SLT))
ALLOCATE (IM6(NMODE_SLT))
ALLOCATE (ZCTOTA(SIZE(PSV,1), SIZE(PSV,2), SIZE(PSV,3), NMODE_SLT))
ALLOCATE (ZM(SIZE(PSV,1), SIZE(PSV,2), SIZE(PSV,3), NMODE_SLT*3))
ALLOCATE (ZSIGMA(SIZE(PSV,1), SIZE(PSV,2), SIZE(PSV,3)))
ALLOCATE (ZINIRADIUS(NMODE_SLT))
ALLOCATE (ZINISIGMA(NMODE_SLT))
ALLOCATE (ZMMIN(NMODE_SLT*3))
!
!
DO JN = 1, NMODE_SLT
  IM0(JN) = 1+(JN-1)*3
  IM3(JN) = 2+(JN-1)*3
  IM6(JN) = 3+(JN-1)*3
  !
  !Get the sea salt mode we are talking about, MODE 2 is treated first, then mode 3, then 1
  !This index is only needed to get the right radius out of the XINIRADIUS_SLT array and the
  !right XINISIG_SLT out of the XINISIG_SLT-array
  IMODEIDX = JPSALTORDER(JN)
  !
  !Convert initial mass median radius to number median radius
  IF (CRGUNITS=="MASS") THEN
    ZINIRADIUS(JN) = XINIRADIUS_SLT(IMODEIDX) * EXP(-3.*(LOG(XINISIG_SLT(IMODEIDX)))**2)
  ELSE
    ZINIRADIUS(JN) = XINIRADIUS_SLT(IMODEIDX)
  END IF
  ZINISIGMA(JN)  = XINISIG_SLT(IMODEIDX)
  !
  ZMMIN(IM0(JN)) = XN0MIN_SLT(IMODEIDX)
! ZRGMIN   = XCOEFRADMIN * ZINIRADIUS(JN)
  ZRGMIN   = ZINIRADIUS(JN)
  ZMMIN(IM3(JN)) = XN0MIN_SLT(IMODEIDX) * (ZRGMIN**3)*EXP(4.5 * LOG(ZINISIGMA(JN))**2) 
  ZMMIN(IM6(JN)) = XN0MIN_SLT(IMODEIDX) * (ZRGMIN**6)*EXP(18. * LOG(ZINISIGMA(JN))**2)
ENDDO
!
!
ZRHOI = XDENSITY_SALT !1.8e3 !++changed alfgr
!ZMI   = XMOLARWEIGHT_SALT*1.D3 !100.  !++changed alfgr
ZMI   = XMOLARWEIGHT_SALT 
ZDEN2MOL = 1E-6 * XAVOGADRO / XMD
ZFAC=(4./3.)*XPI*ZRHOI*1.e-9

! conversion into mol.cm-3
DO JJ=1,SIZE(PSV,4)
  PSV(:,:,:,JJ) =  PSV(:,:,:,JJ) * ZDEN2MOL * PRHODREF(:,:,:)
END DO
!
DO JN=1,NMODE_SLT

!*       1.1    calculate moment 0 from XN0MIN_SLT
!
! change initial proportion 
  IF (JN == 1) ZMMR_SALTN = 0.80 *  ZMMR_SALT
  IF (JN == 2) ZMMR_SALTN = 0.10 *  ZMMR_SALT
  IF (JN == 3) ZMMR_SALTN = 0.10 *  ZMMR_SALT
  DO JK=1,IKU
        ZM(:,:,JK,IM0(JN)) =               &
                ZMMR_SALTN                 &![kg_{sea salt}/kg_{air}  
                *PRHODREF(:,:,JK)          &![kg_{air}/m3_{air}]==> kg/m3
                /XDENSITY_SALT             &![kg__{sea salt}/m3_{sea salt}==>m3_{sea salt}/m3{air}
                *(6.d0/XPI)                 &
                /(2.d0*ZINIRADIUS(JN)*1.d-6)**3 &![particle/m_sea salt^{-3}]==> particle/m3
                *EXP(-4.5*(LOG(XINISIG_SLT(JPSALTORDER(JN))))**2) !Take into account distribution
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
  PSV(:,:,:,2+(JN-1)*3) = ZM(:,:,:,IM3(JN)) * XMD*XPI * 4./3.  / &
                           (ZMI*PRHODREF(:,:,:)*(1.d0/ZRHOI)*XM3TOUM3_SALT)

  PSV(:,:,:,3+(JN-1)*3) = ZM(:,:,:,IM6(JN)) *  XMD / (XAVOGADRO*PRHODREF(:,:,:)*1.d-6)
!
END DO
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
IF (LHOOK) CALL DR_HOOK('SALTLFI_N',1,ZHOOK_HANDLE)
END SUBROUTINE SALTLFI_n
