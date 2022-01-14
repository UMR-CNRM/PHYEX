!     ######spl
      MODULE MODI_TKE_EPS_SOURCES
!     ###########################
INTERFACE
!
      SUBROUTINE TKE_EPS_SOURCES(KKA,KKU,KKL,KMI,PTKEM,PLM,PLEPS,PDP,  &
                    & PTRH,PRHODJ,PDZZ,PDXX,PDYY,PDZX,PDZY,PZZ,        &
                    & PTSTEP,PIMPL,PEXPL,                              &
                    & HTURBLEN,HTURBDIM,                               &
                    & TPFILE,OTURB_DIAG,                               &
                    & PTP,PRTKES,PRTHLS,PCOEF_DISS,PTDIFF,PTDISS,&
                    & PEDR, YDDDH, YDLDDH, YDMDDH,TBUDGETS, KBUDGETS,  &
                    & PTR,PDISS,PRTKESM                                )
!
USE DDH_MIX, ONLY : TYP_DDH
USE YOMLDDH, ONLY : TLDDH
USE YOMMDDH, ONLY : TMDDH
USE MODD_IO, ONLY: TFILEDATA
USE MODD_BUDGET, ONLY : TBUDGETDATA
!
INTEGER,                 INTENT(IN)   :: KKA           !near ground array index  
INTEGER,                 INTENT(IN)   :: KKU           !uppest atmosphere array index
INTEGER,                 INTENT(IN)   :: KKL           !vert. levels type 1=MNH -1=ARO

INTEGER,                 INTENT(IN)   ::  KMI          ! model index number  
REAL, DIMENSION(:,:,:),  INTENT(IN)   ::  PTKEM        ! TKE at t-deltat
REAL, DIMENSION(:,:,:),  INTENT(IN)   ::  PLM          ! mixing length         
REAL, DIMENSION(:,:,:),  INTENT(IN)   ::  PLEPS        ! dissipative length
REAL, DIMENSION(:,:,:),  INTENT(IN)   ::  PRHODJ       ! density * grid volume
REAL, DIMENSION(:,:,:),  INTENT(IN)   ::  PDXX,PDYY,PDZZ,PDZX,PDZY
                                                       ! metric coefficients
REAL, DIMENSION(:,:,:),  INTENT(IN)   ::  PZZ          ! physical height w-pt
REAL,                    INTENT(IN)   ::  PTSTEP       ! Time step 
REAL,                    INTENT(IN)   ::  PEXPL, PIMPL ! Coef. temporal. disc.
CHARACTER(LEN=4),        INTENT(IN)   ::  HTURBDIM     ! dimensionality of the
                                                       ! turbulence scheme
CHARACTER(LEN=4),        INTENT(IN)   ::  HTURBLEN     ! kind of mixing length
TYPE(TFILEDATA),         INTENT(IN)   ::  TPFILE       ! Output file
LOGICAL,                 INTENT(IN)   ::  OTURB_DIAG   ! switch to write some
                                  ! diagnostic fields in the syncronous FM-file
REAL, DIMENSION(:,:,:),  INTENT(INOUT)::  PDP          ! Dyn. prod. of TKE
REAL, DIMENSION(:,:,:),  INTENT(IN)   ::  PTRH
REAL, DIMENSION(:,:,:),  INTENT(IN)   ::  PTP          ! Ther. prod. of TKE
REAL, DIMENSION(:,:,:),  INTENT(INOUT)::  PRTKES       ! RHOD * Jacobian *
                                                       ! TKE at t+deltat
REAL, DIMENSION(:,:,:),  INTENT(INOUT)::  PRTHLS       ! Source of Theta_l
REAL, DIMENSION(:,:,:),  INTENT(IN)   ::  PCOEF_DISS   ! 1/(Cph*Exner)
REAL, DIMENSION(:,:,:),  INTENT(OUT)  ::  PTDIFF       ! Diffusion TKE term
REAL, DIMENSION(:,:,:),  INTENT(OUT)  ::  PTDISS       ! Dissipation TKE term
TYPE(TYP_DDH), INTENT(INOUT) :: YDDDH
TYPE(TLDDH), INTENT(IN) :: YDLDDH
TYPE(TMDDH), INTENT(IN) :: YDMDDH
TYPE(TBUDGETDATA), DIMENSION(KBUDGETS), INTENT(INOUT) :: TBUDGETS
INTEGER, INTENT(IN) :: KBUDGETS
REAL, DIMENSION(:,:,:),  INTENT(OUT), OPTIONAL  ::  PTR          ! Transport prod. of TKE
REAL, DIMENSION(:,:,:),  INTENT(OUT), OPTIONAL  ::  PDISS        ! Dissipation of TKE
REAL, DIMENSION(:,:,:),  INTENT(OUT), OPTIONAL  ::  PEDR         ! EDR 
REAL, DIMENSION(:,:,:),  INTENT(IN),  OPTIONAL  ::  PRTKESM      ! Advection source 
!
!
!
END SUBROUTINE TKE_EPS_SOURCES
!
END INTERFACE
!
END MODULE MODI_TKE_EPS_SOURCES
