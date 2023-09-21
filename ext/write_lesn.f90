!MNH_LIC Copyright 2000-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!######################
module mode_write_les_n
!######################

use modd_field, only: tfieldmetadata_base

implicit none

private

public :: Write_les_n


character(len=:), allocatable :: cgroup
character(len=:), allocatable :: cgroupcomment

logical :: ldoavg    ! Compute and store time average
logical :: ldonorm   ! Compute and store normalized field

type(tfieldmetadata_base) :: tfield
type(tfieldmetadata_base) :: tfieldx
type(tfieldmetadata_base) :: tfieldy

interface Les_diachro_write
  module procedure Les_diachro_write_1D, Les_diachro_write_2D, Les_diachro_write_3D, Les_diachro_write_4D
end interface

contains

!###################################
subroutine  Write_les_n( tpdiafile )
!###################################
!
!
!!****  *WRITE_LES_n* writes the LES final diagnostics for model _n
!!
!!
!!    PURPOSE
!!    -------
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      V. Masson
!!
!!    MODIFICATIONS
!!    -------------
!!      Original         07/02/00
!!                       01/02/01 (D. Gazen) add module MODD_NSV for NSV variable
!!                       06/11/02 (V. Masson) some minor bugs
!!                       01/04/03 (V. Masson) idem
!!                       10/10/09 (P. Aumond) Add user multimaskS
!!                          11/15 (C.Lac) Add production terms of TKE
!!                    10/2016 (C.Lac) Add droplet deposition
!  P. Wautelet 05/2016-04/2018: new data structures and calls for I/O
!  C. Lac         02/2019: add rain fraction as a LES diagnostic
!  P. Wautelet 13/09/2019: budget: simplify and modernize date/time management
!  P. Wautelet 12/10/2020: remove HLES_AVG dummy argument and group all 4 calls
!  P. Wautelet 13/10/2020: bugfix: correct some names for LES_DIACHRO_2PT diagnostics (Ri)
!  P. Wautelet 26/10/2020: bugfix: correct some comments and conditions + add missing RES_RTPZ
!  P. Wautelet 26/10/2020: restructure subroutines to use tfieldmetadata_base type
! --------------------------------------------------------------------------
!
!*      0. DECLARATIONS
!          ------------
!
use modd_conf_n,     only: luserv, luserc, luserr, luseri, lusers, luserg, luserh
use modd_io,         only: tfiledata
use modd_field,      only: NMNHDIM_BUDGET_LES_TIME, NMNHDIM_BUDGET_LES_LEVEL, NMNHDIM_BUDGET_LES_SV, NMNHDIM_BUDGET_LES_MASK, &
                           NMNHDIM_BUDGET_LES_PDF,                                                                            &
                           NMNHDIM_SPECTRA_2PTS_NI, NMNHDIM_SPECTRA_2PTS_NJ,  NMNHDIM_SPECTRA_LEVEL, NMNHDIM_UNUSED,          &
                           TYPEREAL
use modd_grid_n,     only: xdxhat, xdyhat
use modd_nsv,        only: nsv
use modd_les
use modd_les_n
use modd_param_n,    only: ccloud
use modd_param_c2r2, only: ldepoc
USE MODD_PARAM_ICE_n,  only: ldeposc
use modd_parameters, only: XUNDEF

use mode_les_spec_n,            only: Les_spec_n
use mode_modeln_handler,        only: Get_current_model_index
use mode_write_les_budget_n,    only: Write_les_budget_n
use mode_write_les_rt_budget_n, only: Write_les_rt_budget_n
use mode_write_les_sv_budget_n, only: Write_les_sv_budget_n

IMPLICIT NONE
!
!*      0.1  declarations of arguments
!
TYPE(TFILEDATA),  INTENT(IN) :: TPDIAFILE! file to write
!
!
!*      0.2  declaration of local variables
!
INTEGER :: IMASK
!
INTEGER :: JSV       ! scalar loop counter
INTEGER :: JI        ! loop counter
!
character(len=3)                        :: ynum
CHARACTER(len=5)                        :: YGROUP
character(len=7), dimension(nles_masks) :: ymasks
!
logical :: gdoavg    ! Compute and store time average
logical :: gdonorm   ! Compute and store normalized field
REAL, DIMENSION(:,:,:), ALLOCATABLE     :: ZAVG_PTS_ll
REAL, DIMENSION(:,:,:), ALLOCATABLE     :: ZUND_PTS_ll
REAL                                    :: ZCART_PTS_ll
INTEGER                                 :: IMI ! Current model inde
!
!-------------------------------------------------------------------------------
!
IF (.NOT. LLES) RETURN
!
!
!*      1.   Initializations
!            ---------------
!
IMI = GET_CURRENT_MODEL_INDEX()
!
!
!*      1.1  Normalization variables
!            -----------------------
!
IF (CLES_NORM_TYPE/='NONE' ) THEN
  CALL LES_ALLOCATE('XLES_NORM_M',  (/NLES_TIMES/))
  CALL LES_ALLOCATE('XLES_NORM_S',  (/NLES_TIMES/))
  CALL LES_ALLOCATE('XLES_NORM_K',  (/NLES_TIMES/))
  CALL LES_ALLOCATE('XLES_NORM_RHO',(/NLES_TIMES/))
  CALL LES_ALLOCATE('XLES_NORM_RV', (/NLES_TIMES/))
  CALL LES_ALLOCATE('XLES_NORM_SV', (/NLES_TIMES,NSV/))
  CALL LES_ALLOCATE('XLES_NORM_P',  (/NLES_TIMES/))
  !
  IF (CLES_NORM_TYPE=='CONV') THEN
    WHERE (XLES_WSTAR(:)>0.)
      XLES_NORM_M(:)   = XLES_BL_HEIGHT(:)
      XLES_NORM_S(:)   = XLES_NORM_M(:) / XLES_WSTAR(:)
      XLES_NORM_K(:)   = XLES_Q0(:) / XLES_WSTAR(:)
      XLES_NORM_RHO(:) = XLES_MEAN_RHO(1,:,1)
      XLES_NORM_RV(:)  = XLES_E0(:) / XLES_WSTAR(:)
      XLES_NORM_P(:)   = XLES_MEAN_RHO(1,:,1) * XLES_WSTAR(:)**2
    ELSEWHERE
      XLES_NORM_M(:)   = 0.
      XLES_NORM_S(:)   = 0.
      XLES_NORM_K(:)   = 0.
      XLES_NORM_RHO(:) = 0.
      XLES_NORM_RV(:)  = 0.
      XLES_NORM_P(:)   = 0.
    END WHERE
    DO JSV=1,NSV
      WHERE (XLES_WSTAR(:)>0.)
        XLES_NORM_SV(:,JSV)= XLES_SV0(:,JSV) / XLES_WSTAR(:)
      ELSEWHERE
        XLES_NORM_SV(:,JSV)= 0.
      END WHERE
    END DO
  ELSE IF (CLES_NORM_TYPE=='EKMA') THEN
    WHERE (XLES_USTAR(:)>0.)
      XLES_NORM_M(:)   = XLES_BL_HEIGHT(:)
      XLES_NORM_S(:)   = XLES_NORM_M(:) / XLES_USTAR(:)
      XLES_NORM_K(:)   = XLES_Q0(:) / XLES_USTAR(:)
      XLES_NORM_RHO(:) = XLES_MEAN_RHO(1,:,1)
      XLES_NORM_RV(:)  = XLES_E0(:) / XLES_USTAR(:)
      XLES_NORM_P(:)   = XLES_MEAN_RHO(1,:,1) * XLES_USTAR(:)**2
    ELSEWHERE
      XLES_NORM_M(:)   = 0.
      XLES_NORM_S(:)   = 0.
      XLES_NORM_K(:)   = 0.
      XLES_NORM_RHO(:) = 0.
      XLES_NORM_RV(:)  = 0.
      XLES_NORM_P(:)   = 0.
    END WHERE
    DO JSV=1,NSV
      WHERE (XLES_USTAR(:)>0.)
        XLES_NORM_SV(:,JSV)= XLES_SV0(:,JSV) / XLES_USTAR(:)
      ELSEWHERE
        XLES_NORM_SV(:,JSV)= 0.
      END WHERE
    END DO
  ELSE IF (CLES_NORM_TYPE=='MOBU') THEN
    XLES_NORM_M(:) = XLES_MO_LENGTH(:)
    WHERE (XLES_USTAR(:)>0.)
      XLES_NORM_S(:)   = XLES_NORM_M(:) / XLES_USTAR(:)
      XLES_NORM_K(:)   = XLES_Q0(:) / XLES_USTAR(:)
      XLES_NORM_RHO(:) = XLES_MEAN_RHO(1,:,1)
      XLES_NORM_RV(:)  = XLES_E0(:) / XLES_USTAR(:)
      XLES_NORM_P(:)   = XLES_MEAN_RHO(1,:,1) * XLES_USTAR(:)**2
    ELSEWHERE
      XLES_NORM_S(:)   = 0.
      XLES_NORM_K(:)   = 0.
      XLES_NORM_RHO(:) = 0.
      XLES_NORM_RV(:)  = 0.
      XLES_NORM_P(:)   = 0.
    END WHERE
    DO JSV=1,NSV
      WHERE (XLES_USTAR(:)>0.)
        XLES_NORM_SV(:,JSV)= XLES_SV0(:,JSV) / XLES_USTAR(:)
      ELSEWHERE
        XLES_NORM_SV(:,JSV)= 0.
      END WHERE
    END DO
  END IF
END IF
!
!*      1.2  Initializations for WRITE_DIACHRO
!            ---------------------------------
!
NLES_CURRENT_TIMES=NLES_TIMES
!
CALL LES_ALLOCATE('XLES_CURRENT_Z',(/NLES_K/))

XLES_CURRENT_Z(:) = XLES_Z(:)
!
XLES_CURRENT_ZS = XLES_ZS
!
NLES_CURRENT_IINF=NLESn_IINF(IMI)
NLES_CURRENT_ISUP=NLESn_ISUP(IMI)
NLES_CURRENT_JINF=NLESn_JINF(IMI)
NLES_CURRENT_JSUP=NLESn_JSUP(IMI)
!
XLES_CURRENT_DOMEGAX=XDXHAT(1)
XLES_CURRENT_DOMEGAY=XDYHAT(1)

tfield%ngrid = 0 !Not on the Arakawa grid
tfield%ntype = TYPEREAL
!
!*      2.   (z,t) profiles (all masks)
!            --------------
IMASK = 1
ymasks(imask) = 'cart'
IF (LLES_NEB_MASK) THEN
  IMASK=IMASK+1
  ymasks(imask) = 'neb'
  IMASK=IMASK+1
  ymasks(imask) = 'clear'
END IF
IF (LLES_CORE_MASK) THEN
  IMASK=IMASK+1
  ymasks(imask) = 'core'
  IMASK=IMASK+1
  ymasks(imask) = 'env'
END IF
IF (LLES_MY_MASK) THEN
   DO JI=1,NLES_MASKS_USER
        IMASK=IMASK+1
        Write( ynum, '( i3.3 )' ) ji
        ymasks(imask) = 'user' // ynum
   END DO
END IF
IF (LLES_CS_MASK) THEN
  IMASK=IMASK+1
  ymasks(imask) = 'cs1'
  IMASK=IMASK+1
  ymasks(imask) = 'cs2'
  IMASK=IMASK+1
  ymasks(imask) = 'cs3'
END IF
!
!*      2.0  averaging diagnostics
!            ---------------------
!
ALLOCATE(ZAVG_PTS_ll (NLES_K,NLES_TIMES,NLES_MASKS))
ALLOCATE(ZUND_PTS_ll (NLES_K,NLES_TIMES,NLES_MASKS))

ZAVG_PTS_ll(:,:,:) = NLES_AVG_PTS_ll(:,:,:)
ZUND_PTS_ll(:,:,:) = NLES_UND_PTS_ll(:,:,:)
ZCART_PTS_ll       = (NLESn_ISUP(IMI)-NLESn_IINF(IMI)+1) * (NLESn_JSUP(IMI)-NLESn_JINF(IMI)+1)

tfield%ndims = 3
tfield%ndimlist(1)  = NMNHDIM_BUDGET_LES_LEVEL
tfield%ndimlist(2)  = NMNHDIM_BUDGET_LES_TIME
tfield%ndimlist(3)  = NMNHDIM_BUDGET_LES_MASK
tfield%ndimlist(4:) = NMNHDIM_UNUSED

ldoavg  = xles_temp_mean_start /= XUNDEF .and. xles_temp_mean_end /= XUNDEF
ldonorm = .false.

cgroup        = 'Miscellaneous'
cgroupcomment = 'Miscellaneous terms (geometry, various unclassified averaged terms...)'

call Les_diachro_write( tpdiafile, zavg_pts_ll,                'AVG_PTS',  'number of points used for averaging',   '1', ymasks )
call Les_diachro_write( tpdiafile, zavg_pts_ll / zcart_pts_ll, 'AVG_PTSF', 'fraction of points used for averaging', '1', ymasks )
call Les_diachro_write( tpdiafile, zund_pts_ll,                'UND_PTS',  'number of points below orography',      '1', ymasks )
call Les_diachro_write( tpdiafile, zund_pts_ll / zcart_pts_ll, 'UND_PTSF', 'fraction of points below orography',    '1', ymasks )

DEALLOCATE(ZAVG_PTS_ll)
DEALLOCATE(ZUND_PTS_ll)
!
!*      2.1  mean quantities
!            ---------------
!
cgroup = 'Mean'
cgroupcomment = 'Mean vertical profiles of the model variables'

tfield%ndims = 3
tfield%ndimlist(1)  = NMNHDIM_BUDGET_LES_LEVEL
tfield%ndimlist(2)  = NMNHDIM_BUDGET_LES_TIME
tfield%ndimlist(3)  = NMNHDIM_BUDGET_LES_MASK
tfield%ndimlist(4:) = NMNHDIM_UNUSED

ldoavg  = xles_temp_mean_start /= XUNDEF .and. xles_temp_mean_end /= XUNDEF
ldonorm = trim(cles_norm_type) /= 'NONE'

call Les_diachro_write( tpdiafile, XLES_MEAN_U,      'MEAN_U',      'Mean U Profile',                        'm s-1',  ymasks )
call Les_diachro_write( tpdiafile, XLES_MEAN_V,      'MEAN_V',      'Mean V Profile',                        'm s-1',  ymasks )
call Les_diachro_write( tpdiafile, XLES_MEAN_W,      'MEAN_W',      'Mean W Profile',                        'm s-1',  ymasks )
call Les_diachro_write( tpdiafile, XLES_MEAN_P,      'MEAN_PRE',    'Mean pressure Profile',                 'Pa',     ymasks )
call Les_diachro_write( tpdiafile, XLES_MEAN_DP,     'MEAN_DP',     'Mean Dyn production TKE Profile',       'm2 s-3', ymasks )
call Les_diachro_write( tpdiafile, XLES_MEAN_TP,     'MEAN_TP',     'Mean Thermal production TKE Profile',   'm2 s-3', ymasks )
call Les_diachro_write( tpdiafile, XLES_MEAN_TR,     'MEAN_TR',     'Mean transport production TKE Profile', 'm2 s-3', ymasks )
call Les_diachro_write( tpdiafile, XLES_MEAN_DISS,   'MEAN_DISS',   'Mean Dissipation TKE Profile',          'm2 s-3', ymasks )
call Les_diachro_write( tpdiafile, XLES_MEAN_LM,     'MEAN_LM',     'Mean mixing length Profile',            'm',      ymasks )
call Les_diachro_write( tpdiafile, XLES_MEAN_RHO,    'MEAN_RHO',    'Mean density Profile',                  'kg m-3', ymasks )
call Les_diachro_write( tpdiafile, XLES_MEAN_Th,     'MEAN_TH',     'Mean potential temperature Profile',    'K',      ymasks )
call Les_diachro_write( tpdiafile, XLES_MEAN_Mf,     'MEAN_MF',     'Mass-flux Profile',                     'm s-1',  ymasks )
if ( luserc ) &
call Les_diachro_write( tpdiafile, XLES_MEAN_Thl,    'MEAN_THL',    'Mean liquid potential temperature Profile',  'K', ymasks )
if ( luserv ) &
call Les_diachro_write( tpdiafile, XLES_MEAN_Thv,    'MEAN_THV',    'Mean virtual potential temperature Profile', 'K', ymasks )
if ( luserc ) &
call Les_diachro_write( tpdiafile, XLES_MEAN_Rt,     'MEAN_RT',     'Mean Rt Profile', 'kg kg-1', ymasks )
if ( luserv ) &
call Les_diachro_write( tpdiafile, XLES_MEAN_Rv,     'MEAN_RV',     'Mean Rv Profile', 'kg kg-1', ymasks )
if ( luserv ) &
call Les_diachro_write( tpdiafile, XLES_MEAN_Rehu,   'MEAN_REHU',   'Mean Rh Profile', 'percent', ymasks )
if ( luserv ) &
call Les_diachro_write( tpdiafile, XLES_MEAN_Qs,     'MEAN_QS',     'Mean Qs Profile', 'kg kg-1', ymasks )
if ( luserc ) &
call Les_diachro_write( tpdiafile, XLES_MEAN_KHt,    'MEAN_KHT',    'Eddy-diffusivity (temperature) Profile', 'm2 s-1', ymasks )
if ( luserc ) &
call Les_diachro_write( tpdiafile, XLES_MEAN_KHr,    'MEAN_KHR',    'Eddy-diffusivity (vapor) Profile',      'm2 s-1', ymasks )
if ( luserc ) &
call Les_diachro_write( tpdiafile, XLES_MEAN_Rc,     'MEAN_RC',     'Mean Rc Profile',              'kg kg-1', ymasks )
if ( luserc ) &
call Les_diachro_write( tpdiafile, XLES_MEAN_Cf,     'MEAN_CF',     'Mean Cf Profile',              '1',       ymasks )
if ( luserc ) &
call Les_diachro_write( tpdiafile, XLES_MEAN_INDCf,  'MEAN_INDCF',  'Mean Cf>1-6 Profile (0 or 1)', '1',       ymasks )
if ( luserc ) &
call Les_diachro_write( tpdiafile, XLES_MEAN_INDCf2, 'MEAN_INDCF2', 'Mean Cf>1-5 Profile (0 or 1)', '1',       ymasks )
if ( luserr ) &
call Les_diachro_write( tpdiafile, XLES_MEAN_Rr,     'MEAN_RR',     'Mean Rr Profile',              'kg kg-1', ymasks )
if ( luserr ) &
call Les_diachro_write( tpdiafile, XLES_MEAN_RF,     'MEAN_RF',     'Mean RF Profile',              '1',       ymasks )
if ( luseri ) &
call Les_diachro_write( tpdiafile, XLES_MEAN_Ri,     'MEAN_RI',     'Mean Ri Profile',              'kg kg-1', ymasks )
if ( luseri ) &
call Les_diachro_write( tpdiafile, XLES_MEAN_If,     'MEAN_IF',     'Mean If Profile',              '1',       ymasks )
if ( lusers ) &
call Les_diachro_write( tpdiafile, XLES_MEAN_Rs,     'MEAN_RS',     'Mean Rs Profile',              'kg kg-1', ymasks )
if ( luserg ) &
call Les_diachro_write( tpdiafile, XLES_MEAN_Rg,     'MEAN_RG',     'Mean Rg Profile',              'kg kg-1', ymasks )
if ( luserh ) &
call Les_diachro_write( tpdiafile, XLES_MEAN_Rh,     'MEAN_RH',     'Mean Rh Profile',              'kg kg-1', ymasks )

if ( nsv > 0 ) then
  tfield%ndims = 4
  tfield%ndimlist(1)  = NMNHDIM_BUDGET_LES_LEVEL
  tfield%ndimlist(2)  = NMNHDIM_BUDGET_LES_TIME
  tfield%ndimlist(3)  = NMNHDIM_BUDGET_LES_MASK
  tfield%ndimlist(4)  = NMNHDIM_BUDGET_LES_SV
  tfield%ndimlist(5:) = NMNHDIM_UNUSED

  call Les_diachro_write( tpdiafile, XLES_MEAN_Sv, 'MEAN_SV', 'Mean Sv Profiles', 'kg kg-1', ymasks )

  tfield%ndims = 3
  !tfield%ndimlist(1)  = NMNHDIM_BUDGET_LES_LEVEL
  !tfield%ndimlist(2)  = NMNHDIM_BUDGET_LES_TIME
  !tfield%ndimlist(3)  = NMNHDIM_BUDGET_LES_MASK
  tfield%ndimlist(4)  = NMNHDIM_UNUSED
  !tfield%ndimlist(5:) = NMNHDIM_UNUSED
end if

call Les_diachro_write( tpdiafile, XLES_MEAN_WIND, 'MEANWIND',       'Profile of Mean Modulus of Wind', 'm s-1',      ymasks )
call Les_diachro_write( tpdiafile, XLES_RESOLVED_MASSFX, 'MEANMSFX', 'Total updraft mass flux',         'kg m-2 s-1', ymasks )

if ( lles_pdf ) then
  cgroup = 'PDF'
  cgroupcomment = ''

  tfield%ndims = 4
  !tfield%ndimlist(1)  = NMNHDIM_BUDGET_LES_LEVEL
  !tfield%ndimlist(2)  = NMNHDIM_BUDGET_LES_TIME
  !tfield%ndimlist(3)  = NMNHDIM_BUDGET_LES_MASK
  tfield%ndimlist(4)  = NMNHDIM_BUDGET_LES_PDF
  tfield%ndimlist(5:) = NMNHDIM_UNUSED

  call Les_diachro_write( tpdiafile,   XLES_PDF_TH,  'PDF_TH',  'Pdf potential temperature Profiles', '1', ymasks )
  call Les_diachro_write( tpdiafile,   XLES_PDF_W,   'PDF_W',   'Pdf vertical velocity Profiles',     '1', ymasks )
  call Les_diachro_write( tpdiafile,   XLES_PDF_THV, 'PDF_THV', 'Pdf virtual pot. temp. Profiles',    '1', ymasks )
  if ( luserv ) &
  call Les_diachro_write( tpdiafile,   XLES_PDF_RV,  'PDF_RV',  'Pdf Rv Profiles',                    '1', ymasks )
  if ( luserc ) then
    call Les_diachro_write( tpdiafile, XLES_PDF_RC,  'PDF_RC',  'Pdf Rc Profiles',                    '1', ymasks )
    call Les_diachro_write( tpdiafile, XLES_PDF_RT,  'PDF_RT',  'Pdf Rt Profiles',                    '1', ymasks )
    call Les_diachro_write( tpdiafile, XLES_PDF_THL, 'PDF_THL', 'Pdf Thl Profiles',                   '1', ymasks )
  end if
  if ( luserr ) &
  call Les_diachro_write( tpdiafile,   XLES_PDF_RR,  'PDF_RR',  'Pdf Rr Profiles',                    '1', ymasks )
  if ( luseri ) &
  call Les_diachro_write( tpdiafile,   XLES_PDF_RI,  'PDF_RI',  'Pdf Ri Profiles',                    '1', ymasks )
  if ( lusers ) &
  call Les_diachro_write( tpdiafile,   XLES_PDF_RS,  'PDF_RS',  'Pdf Rs Profiles',                    '1', ymasks )
  if ( luserg ) &
  call Les_diachro_write( tpdiafile,   XLES_PDF_RG,  'PDF_RG',  'Pdf Rg Profiles',                    '1', ymasks )
end if
!
!*      2.2  resolved quantities
!            -------------------
!
if ( lles_resolved ) then
  !Prepare metadata (used in Les_diachro_write calls)
  ldoavg  = xles_temp_mean_start /= XUNDEF .and. xles_temp_mean_end /= XUNDEF
  ldonorm = trim(cles_norm_type) /= 'NONE'

  cgroup = 'Resolved'
  cgroupcomment = 'Mean vertical profiles of the resolved fluxes, variances and covariances'

  tfield%ndims = 3
  tfield%ndimlist(1)  = NMNHDIM_BUDGET_LES_LEVEL
  tfield%ndimlist(2)  = NMNHDIM_BUDGET_LES_TIME
  tfield%ndimlist(3)  = NMNHDIM_BUDGET_LES_MASK
  tfield%ndimlist(4:) = NMNHDIM_UNUSED

  call Les_diachro_write( tpdiafile, XLES_RESOLVED_U2, 'RES_U2',  'Resolved <u2> variance',        'm2 s-2', ymasks )
  call Les_diachro_write( tpdiafile, XLES_RESOLVED_V2, 'RES_V2',  'Resolved <v2> variance',        'm2 s-2', ymasks )
  call Les_diachro_write( tpdiafile, XLES_RESOLVED_W2, 'RES_W2',  'Resolved <w2> variance',        'm2 s-2', ymasks )
  call Les_diachro_write( tpdiafile, XLES_RESOLVED_UV, 'RES_UV',  'Resolved <uv> Flux',            'm2 s-2', ymasks )
  call Les_diachro_write( tpdiafile, XLES_RESOLVED_WU, 'RES_WU',  'Resolved <wu> Flux',            'm2 s-2', ymasks )
  call Les_diachro_write( tpdiafile, XLES_RESOLVED_WV, 'RES_WV',  'Resolved <wv> Flux',            'm2 s-2', ymasks )
  call Les_diachro_write( tpdiafile, XLES_RESOLVED_Ke, 'RES_KE',  'Resolved TKE Profile',          'm2 s-2', ymasks )
  call Les_diachro_write( tpdiafile, XLES_RESOLVED_P2, 'RES_P2',  'Resolved pressure variance',    'Pa2',    ymasks )
  call Les_diachro_write( tpdiafile, XLES_RESOLVED_UP, 'RES_UPZ', 'Resolved <up> horizontal Flux', 'Pa s-1', ymasks )
  call Les_diachro_write( tpdiafile, XLES_RESOLVED_VP, 'RES_VPZ', 'Resolved <vp> horizontal Flux', 'Pa s-1', ymasks )
  call Les_diachro_write( tpdiafile, XLES_RESOLVED_WP, 'RES_WPZ', 'Resolved <wp> vertical Flux',   'Pa s-1', ymasks )

  if ( luserv ) &
  call Les_diachro_write( tpdiafile, XLES_RESOLVED_ThThv, 'RES_THTV', &
                          'Resolved potential temperature - virtual potential temperature covariance',        'K2', ymasks )
  if ( luserc ) &
  call Les_diachro_write( tpdiafile, XLES_RESOLVED_ThlThv, 'RES_TLTV', &
                          'Resolved liquid potential temperature - virtual potential temperature covariance', 'K2', ymasks )
  call Les_diachro_write( tpdiafile, XLES_RESOLVED_Th2, 'RES_TH2', 'Resolved potential temperature variance', 'K2', ymasks )
  if ( luserc ) &
  call Les_diachro_write( tpdiafile, XLES_RESOLVED_Thl2, 'RES_THL2', 'Resolved liquid potential temperature variance', 'K2',&
                          ymasks )
  call Les_diachro_write( tpdiafile, XLES_RESOLVED_UTh, 'RES_UTH', 'Resolved <uth> horizontal Flux', 'm K s-1', ymasks )
  call Les_diachro_write( tpdiafile, XLES_RESOLVED_VTh, 'RES_VTH', 'Resolved <vth> horizontal Flux', 'm K s-1', ymasks )
  call Les_diachro_write( tpdiafile, XLES_RESOLVED_WTh, 'RES_WTH', 'Resolved <wth> vertical Flux',   'm K s-1', ymasks )

  if ( luserc ) then
    call Les_diachro_write( tpdiafile, XLES_RESOLVED_UThl, 'RES_UTHL', 'Resolved <uthl> horizontal Flux', 'm K s-1',       ymasks )
    call Les_diachro_write( tpdiafile, XLES_RESOLVED_VThl, 'RES_VTHL', 'Resolved <vthl> horizontal Flux', 'm K s-1',       ymasks )
    call Les_diachro_write( tpdiafile, XLES_RESOLVED_WThl, 'RES_WTHL', 'Resolved <wthl> vertical Flux',   'm K s-1',       ymasks )
    call Les_diachro_write( tpdiafile, XLES_RESOLVED_Rt2,  'RES_RT2',  'Resolved total water variance',   'kg2 kg-2',      ymasks )
    call Les_diachro_write( tpdiafile, XLES_RESOLVED_WRt,  'RES_WRT',  'Resolved <wrt> vertical Flux',    'm kg kg-1 s-1', ymasks )
  end if

  if ( luserv ) then
    call Les_diachro_write( tpdiafile, XLES_RESOLVED_UThv,  'RES_UTHV', 'Resolved <uthv> horizontal Flux', 'm K s-1',       ymasks )
    call Les_diachro_write( tpdiafile, XLES_RESOLVED_VThv,  'RES_VTHV', 'Resolved <vthv> horizontal Flux', 'm K s-1',       ymasks )
    call Les_diachro_write( tpdiafile, XLES_RESOLVED_WThv,  'RES_WTHV', 'Resolved <wthv> vertical Flux',   'm K s-1',       ymasks )
    call Les_diachro_write( tpdiafile, XLES_RESOLVED_Rv2,   'RES_RV2',  'Resolved water vapor variance',   'kg2 kg-2',      ymasks )
    call Les_diachro_write( tpdiafile, XLES_RESOLVED_ThRv,  'RES_THRV', 'Resolved <thrv> covariance',      'K kg kg-1',     ymasks )
    if ( luserc ) &
    call Les_diachro_write( tpdiafile, XLES_RESOLVED_ThlRv, 'RES_TLRV', 'Resolved <thlrv> covariance',     'K kg kg-1',     ymasks )
    call Les_diachro_write( tpdiafile, XLES_RESOLVED_ThvRv, 'RES_TVRV', 'Resolved <thvrv> covariance',     'K kg kg-1',     ymasks )
    call Les_diachro_write( tpdiafile, XLES_RESOLVED_URv,   'RES_URV',  'Resolved <urv> horizontal flux',  'm kg kg-1 s-1', ymasks )
    call Les_diachro_write( tpdiafile, XLES_RESOLVED_VRv,   'RES_VRV',  'Resolved <vrv> horizontal flux',  'm kg kg-1 s-1', ymasks )
    call Les_diachro_write( tpdiafile, XLES_RESOLVED_WRv,   'RES_WRV',  'Resolved <wrv> vertical flux',    'm kg kg-1 s-1', ymasks )
  end if

  if ( luserc ) then
    call Les_diachro_write( tpdiafile, XLES_RESOLVED_Rc2,   'RES_RC2',  'Resolved cloud water variance',  'kg2 kg-2',      ymasks )
    call Les_diachro_write( tpdiafile, XLES_RESOLVED_ThRc,  'RES_THRC', 'Resolved <thrc> covariance',     'K kg kg-1',     ymasks )
    call Les_diachro_write( tpdiafile, XLES_RESOLVED_ThlRc, 'RES_TLRC', 'Resolved <thlrc> covariance',    'K kg kg-1',     ymasks )
    call Les_diachro_write( tpdiafile, XLES_RESOLVED_ThvRc, 'RES_TVRC', 'Resolved <thvrc> covariance',    'K kg kg-1',     ymasks )
    call Les_diachro_write( tpdiafile, XLES_RESOLVED_URc,   'RES_URC',  'Resolved <urc> horizontal flux', 'm kg kg-1 s-1', ymasks )
    call Les_diachro_write( tpdiafile, XLES_RESOLVED_VRc,   'RES_VRC',  'Resolved <vrc> horizontal flux', 'm kg kg-1 s-1', ymasks )
    call Les_diachro_write( tpdiafile, XLES_RESOLVED_WRc,   'RES_WRC',  'Resolved <wrc> vertical flux',   'm kg kg-1 s-1', ymasks )
  end if

  if ( luseri ) then
    call Les_diachro_write( tpdiafile, XLES_RESOLVED_Ri2,   'RES_RI2',  'Resolved cloud ice variance',    'kg2 kg-2',      ymasks )
    call Les_diachro_write( tpdiafile, XLES_RESOLVED_ThRi,  'RES_THRI', 'Resolved <thri> covariance',     'K kg kg-1',     ymasks )
    call Les_diachro_write( tpdiafile, XLES_RESOLVED_ThlRi, 'RES_TLRI', 'Resolved <thlri> covariance',    'K kg kg-1',     ymasks )
    call Les_diachro_write( tpdiafile, XLES_RESOLVED_ThvRi, 'RES_TVRI', 'Resolved <thvri> covariance',    'K kg kg-1',     ymasks )
    call Les_diachro_write( tpdiafile, XLES_RESOLVED_URi,   'RES_URI',  'Resolved <uri> horizontal flux', 'm kg kg-1 s-1', ymasks )
    call Les_diachro_write( tpdiafile, XLES_RESOLVED_VRi,   'RES_VRI',  'Resolved <vri> horizontal flux', 'm kg kg-1 s-1', ymasks )
    call Les_diachro_write( tpdiafile, XLES_RESOLVED_WRi,   'RES_WRI',  'Resolved <wri> vertical flux',   'm kg kg-1 s-1', ymasks )
  end if

  if ( luserr ) then
    call Les_diachro_write( tpdiafile, XLES_RESOLVED_WRr,   'RES_WRR',   'Resolved <wrr> vertical flux', 'm kg kg-1 s-1', ymasks )
    call Les_diachro_write( tpdiafile, XLES_INPRR3D,        'INPRR3D',   'Precipitation flux',           'm s-1',         ymasks )
    call Les_diachro_write( tpdiafile, XLES_MAX_INPRR3D,    'MAXINPR3D', 'Max Precip flux',              'm s-1',         ymasks )
    call Les_diachro_write( tpdiafile, XLES_EVAP3D,         'EVAP3D',    'Evaporation profile',          'kg kg-1 s-1',   ymasks )
  end if

  if ( nsv > 0 ) then
    tfield%ndims = 4
    tfield%ndimlist(1)  = NMNHDIM_BUDGET_LES_LEVEL
    tfield%ndimlist(2)  = NMNHDIM_BUDGET_LES_TIME
    tfield%ndimlist(3)  = NMNHDIM_BUDGET_LES_MASK
    tfield%ndimlist(4)  = NMNHDIM_BUDGET_LES_SV
    tfield%ndimlist(5:) = NMNHDIM_UNUSED

    call Les_diachro_write( tpdiafile, XLES_RESOLVED_Sv2,   'RES_SV2',  'Resolved scalar variables variances', 'kg2 kg-2', ymasks )
    call Les_diachro_write( tpdiafile, XLES_RESOLVED_ThSv,  'RES_THSV', 'Resolved <ThSv> variance',  'K kg kg-1',          ymasks )
    if ( luserc ) &
    call Les_diachro_write( tpdiafile, XLES_RESOLVED_ThlSv, 'RES_TLSV', 'Resolved <ThlSv> variance', 'K kg kg-1',          ymasks )
    if ( luserv ) &
    call Les_diachro_write( tpdiafile, XLES_RESOLVED_ThvSv, 'RES_TVSV', 'Resolved <ThvSv> variance', 'K kg kg-1',          ymasks )
    call Les_diachro_write( tpdiafile, XLES_RESOLVED_USv,   'RES_USV',  'Resolved <uSv> horizontal flux', 'm kg kg-1 s-1', ymasks )
    call Les_diachro_write( tpdiafile, XLES_RESOLVED_VSv,   'RES_VSV',  'Resolved <vSv> horizontal flux', 'm kg kg-1 s-1', ymasks )
    call Les_diachro_write( tpdiafile, XLES_RESOLVED_WSv,   'RES_WSV',  'Resolved <wSv> vertical flux', 'm kg kg-1 s-1',   ymasks )

    tfield%ndims = 3
    !tfield%ndimlist(1)  = NMNHDIM_BUDGET_LES_LEVEL
    !tfield%ndimlist(2)  = NMNHDIM_BUDGET_LES_TIME
    !tfield%ndimlist(3)  = NMNHDIM_BUDGET_LES_MASK
    tfield%ndimlist(4)  = NMNHDIM_UNUSED
    !tfield%ndimlist(5:) = NMNHDIM_UNUSED
  end if

  call Les_diachro_write( tpdiafile, XLES_RESOLVED_U3, 'RES_U3', 'Resolved <u3>', 'm3 s-3', ymasks )
  call Les_diachro_write( tpdiafile, XLES_RESOLVED_V3, 'RES_V3', 'Resolved <v3>', 'm3 s-3', ymasks )
  call Les_diachro_write( tpdiafile, XLES_RESOLVED_W3, 'RES_W3', 'Resolved <w3>', 'm3 s-3', ymasks )
  call Les_diachro_write( tpdiafile, XLES_RESOLVED_U4, 'RES_U4', 'Resolved <u4>', 'm4 s-4', ymasks )
  call Les_diachro_write( tpdiafile, XLES_RESOLVED_V4, 'RES_V4', 'Resolved <v4>', 'm4 s-4', ymasks )
  call Les_diachro_write( tpdiafile, XLES_RESOLVED_W4, 'RES_W4', 'Resolved <w4>', 'm4 s-4', ymasks )

  call Les_diachro_write( tpdiafile, XLES_RESOLVED_WThl2, 'RES_WTL2', 'Resolved <wThl2>', 'm K2 s-1', ymasks )
  call Les_diachro_write( tpdiafile, XLES_RESOLVED_W2Thl, 'RES_W2TL', 'Resolved <w2Thl>', 'm2 K s-2', ymasks )

  if ( luserv ) then
    call Les_diachro_write( tpdiafile, XLES_RESOLVED_WRv2,   'RES_WRV2', 'Resolved <wRv2>',   'm kg2 kg-2 s-1',  ymasks )
    call Les_diachro_write( tpdiafile, XLES_RESOLVED_W2Rv,   'RES_W2RV', 'Resolved <w2Rv>',   'm2 kg kg-1 s-2',  ymasks )
    call Les_diachro_write( tpdiafile, XLES_RESOLVED_WRt2,   'RES_WRT2', 'Resolved <wRt2>',   'm kg2 kg-2 s-1',  ymasks )
    call Les_diachro_write( tpdiafile, XLES_RESOLVED_W2Rt,   'RES_W2RT', 'Resolved <w2Rt>',   'm2 kg kg-1 s-2',  ymasks )
    call Les_diachro_write( tpdiafile, XLES_RESOLVED_WThlRv, 'RE_WTLRV', 'Resolved <wThlRv>', 'm K kg kg-1 s-1', ymasks )
    call Les_diachro_write( tpdiafile, XLES_RESOLVED_WThlRt, 'RE_WTLRT', 'Resolved <wThlRt>', 'm K kg kg-1 s-1', ymasks )
  end if

  if ( luserc ) then
    call Les_diachro_write( tpdiafile, XLES_RESOLVED_WRc2,   'RES_WRC2', 'Resolved <wRc2>',   'm kg2 kg-2 s-1',  ymasks )
    call Les_diachro_write( tpdiafile, XLES_RESOLVED_W2Rc,   'RES_W2RC', 'Resolved <w2Rc>',   'm2 kg kg-1 s-2',  ymasks )
    call Les_diachro_write( tpdiafile, XLES_RESOLVED_WThlRc, 'RE_WTLRC', 'Resolved <wThlRc>', 'm K kg kg-1 s-1', ymasks )
    call Les_diachro_write( tpdiafile, XLES_RESOLVED_WRvRc,  'RE_WRVRC', 'Resolved <wRvRc>',  'm kg2 kg-2 s-1',  ymasks )
  end if

  if ( luseri ) then
    call Les_diachro_write( tpdiafile, XLES_RESOLVED_WRi2,   'RES_WRI2', 'Resolved <wRi2>',   'm kg2 kg-2 s-1',  ymasks )
    call Les_diachro_write( tpdiafile, XLES_RESOLVED_W2Ri,   'RES_W2RI', 'Resolved <w2Ri>',   'm2 kg kg-1 s-2',  ymasks )
    call Les_diachro_write( tpdiafile, XLES_RESOLVED_WThlRi, 'RE_WTLRI', 'Resolved <wThlRi>', 'm K kg kg-1 s-1', ymasks )
    call Les_diachro_write( tpdiafile, XLES_RESOLVED_WRvRi,  'RE_WRVRI', 'Resolved <wRvRi>',  'm kg2 kg-2 s-1',  ymasks )
  end if

  if ( nsv > 0 ) then
    tfield%ndims = 4
    tfield%ndimlist(1)  = NMNHDIM_BUDGET_LES_LEVEL
    tfield%ndimlist(2)  = NMNHDIM_BUDGET_LES_TIME
    tfield%ndimlist(3)  = NMNHDIM_BUDGET_LES_MASK
    tfield%ndimlist(4)  = NMNHDIM_BUDGET_LES_SV
    tfield%ndimlist(5:) = NMNHDIM_UNUSED

    call Les_diachro_write( tpdiafile, XLES_RESOLVED_WSv2,   'RES_WSV2', 'Resolved <wSv2>',   'm kg2 kg-2 s-1',  ymasks )
    call Les_diachro_write( tpdiafile, XLES_RESOLVED_W2Sv,   'RES_W2SV', 'Resolved <w2Sv>',   'm2 kg kg-1 s-2',  ymasks )
    call Les_diachro_write( tpdiafile, XLES_RESOLVED_WThlSv, 'RE_WTLSV', 'Resolved <wThlSv>', 'm K kg kg-1 s-1', ymasks )
    if ( luserv ) &
    call Les_diachro_write( tpdiafile, XLES_RESOLVED_WRvSv,  'RE_WRVSV', 'Resolved <wRvSv>',  'm kg2 kg-2 s-1',  ymasks )

    tfield%ndims = 3
    !tfield%ndimlist(1)  = NMNHDIM_BUDGET_LES_LEVEL
    !tfield%ndimlist(2)  = NMNHDIM_BUDGET_LES_TIME
    !tfield%ndimlist(3)  = NMNHDIM_BUDGET_LES_MASK
    tfield%ndimlist(4)  = NMNHDIM_UNUSED
    !tfield%ndimlist(5:) = NMNHDIM_UNUSED
  end if

  call Les_diachro_write( tpdiafile, XLES_RESOLVED_ThlPz, 'RES_TLPZ', 'Resolved <Thldp/dz>', 'K Pa m-1',        ymasks )
  if ( luserc ) &
  call Les_diachro_write( tpdiafile, XLES_RESOLVED_RtPz,  'RES_RTPZ', 'Resolved <Rtdp/dz>',  'kg2 kg-2 Pa m-1', ymasks )
  if ( luserv ) &
  call Les_diachro_write( tpdiafile, XLES_RESOLVED_RvPz,  'RES_RVPZ', 'Resolved <Rvdp/dz>',  'kg2 kg-2 Pa m-1', ymasks )
  if ( luserc ) &
  call Les_diachro_write( tpdiafile, XLES_RESOLVED_RcPz,  'RES_RCPZ', 'Resolved <Rcdp/dz>',  'kg2 kg-2 Pa m-1', ymasks )
  if ( luseri ) &
  call Les_diachro_write( tpdiafile, XLES_RESOLVED_RiPz,  'RES_RIPZ', 'Resolved <Ridp/dz>',  'kg2 kg-2 Pa m-1', ymasks )

  if ( nsv > 0 ) then
    tfield%ndims = 4
    tfield%ndimlist(1)  = NMNHDIM_BUDGET_LES_LEVEL
    tfield%ndimlist(2)  = NMNHDIM_BUDGET_LES_TIME
    tfield%ndimlist(3)  = NMNHDIM_BUDGET_LES_MASK
    tfield%ndimlist(4)  = NMNHDIM_BUDGET_LES_SV
    tfield%ndimlist(5:) = NMNHDIM_UNUSED

    call Les_diachro_write( tpdiafile, XLES_RESOLVED_SvPz, 'RES_SVPZ', 'Resolved <Svdp/dz>', 'kg2 kg-2 Pa m-1', ymasks )

    tfield%ndims = 3
    !tfield%ndimlist(1)  = NMNHDIM_BUDGET_LES_LEVEL
    !tfield%ndimlist(2)  = NMNHDIM_BUDGET_LES_TIME
    !tfield%ndimlist(3)  = NMNHDIM_BUDGET_LES_MASK
    tfield%ndimlist(4)  = NMNHDIM_UNUSED
    !tfield%ndimlist(5:) = NMNHDIM_UNUSED
  end if

  call Les_diachro_write( tpdiafile, XLES_RESOLVED_UKe, 'RES_UKE', 'Resolved flux of resolved kinetic energy', 'm3 s-3', ymasks )
  call Les_diachro_write( tpdiafile, XLES_RESOLVED_VKe, 'RES_VKE', 'Resolved flux of resolved kinetic energy', 'm3 s-3', ymasks )
  call Les_diachro_write( tpdiafile, XLES_RESOLVED_WKe, 'RES_WKE', 'Resolved flux of resolved kinetic energy', 'm3 s-3', ymasks )
end if
!
!
!*      2.3  subgrid quantities
!            ------------------
!
if ( lles_subgrid ) then
  !Prepare metadata (used in Les_diachro_write calls)
  ldoavg  = xles_temp_mean_start /= XUNDEF .and. xles_temp_mean_end /= XUNDEF
  ldonorm = trim(cles_norm_type) /= 'NONE'

  cgroup = 'Subgrid'
  cgroupcomment = 'Mean vertical profiles of the subgrid fluxes, variances and covariances'

  tfield%ndims = 3
  tfield%ndimlist(1)  = NMNHDIM_BUDGET_LES_LEVEL
  tfield%ndimlist(2)  = NMNHDIM_BUDGET_LES_TIME
  tfield%ndimlist(3)  = NMNHDIM_BUDGET_LES_MASK
  tfield%ndimlist(4:) = NMNHDIM_UNUSED

  call Les_diachro_write( tpdiafile, XLES_SUBGRID_Tke,  'SBG_TKE',  'Subgrid TKE',           'm2 s-2', ymasks )
  call Les_diachro_write( tpdiafile, XLES_SUBGRID_U2,   'SBG_U2',   'Subgrid <u2> variance', 'm2 s-2', ymasks )
  call Les_diachro_write( tpdiafile, XLES_SUBGRID_V2,   'SBG_V2',   'Subgrid <v2> variance', 'm2 s-2', ymasks )
  call Les_diachro_write( tpdiafile, XLES_SUBGRID_W2,   'SBG_W2',   'Subgrid <w2> variance', 'm2 s-2', ymasks )
  call Les_diachro_write( tpdiafile, XLES_SUBGRID_UV,   'SBG_UV',   'Subgrid <uv> flux',     'm2 s-2', ymasks )
  call Les_diachro_write( tpdiafile, XLES_SUBGRID_WU,   'SBG_WU',   'Subgrid <wu> flux',     'm2 s-2', ymasks )
  call Les_diachro_write( tpdiafile, XLES_SUBGRID_WV,   'SBG_WV',   'Subgrid <wv> flux',     'm2 s-2', ymasks )
  call Les_diachro_write( tpdiafile, XLES_SUBGRID_Thl2, 'SBG_THL2', 'Subgrid liquid potential temperature variance', &
                          'K2', ymasks )
  call Les_diachro_write( tpdiafile, XLES_SUBGRID_UThl, 'SBG_UTHL', 'Subgrid horizontal flux of liquid potential temperature',  &
                          'm K s-1', ymasks )
  call Les_diachro_write( tpdiafile, XLES_SUBGRID_VThl, 'SBG_VTHL', 'Subgrid horizontal flux of liquid potential temperature',  &
                          'm K s-1', ymasks )
  call Les_diachro_write( tpdiafile, XLES_SUBGRID_WThl, 'SBG_WTHL', 'Subgrid vertical flux of liquid potential temperature', &
                          'm K s-1', ymasks )
  call Les_diachro_write( tpdiafile, XLES_SUBGRID_WP,   'SBG_WP',   'Subgrid <wp> vertical Flux', 'm Pa s-1', ymasks )

  call Les_diachro_write( tpdiafile, XLES_SUBGRID_THLUP_MF, 'THLUP_MF', 'Subgrid <thl> of updraft',    'K',          ymasks )
  call Les_diachro_write( tpdiafile, XLES_SUBGRID_RTUP_MF,  'RTUP_MF',  'Subgrid <rt> of updraft',     'kg kg-1',    ymasks )
  call Les_diachro_write( tpdiafile, XLES_SUBGRID_RVUP_MF,  'RVUP_MF',  'Subgrid <rv> of updraft',     'kg kg-1',    ymasks )
  call Les_diachro_write( tpdiafile, XLES_SUBGRID_RCUP_MF,  'RCUP_MF',  'Subgrid <rc> of updraft',     'kg kg-1',    ymasks )
  call Les_diachro_write( tpdiafile, XLES_SUBGRID_RIUP_MF,  'RIUP_MF',  'Subgrid <ri> of updraft',     'kg kg-1',    ymasks )
  call Les_diachro_write( tpdiafile, XLES_SUBGRID_WUP_MF,   'WUP_MF',   'Subgrid <w> of updraft',      'm s-1',      ymasks )
  call Les_diachro_write( tpdiafile, XLES_SUBGRID_MASSFLUX, 'MAFLX_MF', 'Subgrid <MF> of updraft',     'kg m-2 s-1', ymasks )
  call Les_diachro_write( tpdiafile, XLES_SUBGRID_DETR,     'DETR_MF',  'Subgrid <detr> of updraft',   'kg m-3 s-1', ymasks )
  call Les_diachro_write( tpdiafile, XLES_SUBGRID_ENTR,     'ENTR_MF',  'Subgrid <entr> of updraft',   'kg m-3 s-1', ymasks )
  call Les_diachro_write( tpdiafile, XLES_SUBGRID_FRACUP,   'FRCUP_MF', 'Subgrid <FracUp> of updraft', '1',          ymasks )
  call Les_diachro_write( tpdiafile, XLES_SUBGRID_THVUP_MF, 'THVUP_MF', 'Subgrid <thv> of updraft',    'K',          ymasks )
  call Les_diachro_write( tpdiafile, XLES_SUBGRID_WTHLMF,   'WTHL_MF',  'Subgrid <wthl> of mass flux convection scheme', &
                          'm K s-1', ymasks )
  call Les_diachro_write( tpdiafile, XLES_SUBGRID_WRTMF,    'WRT_MF',   'Subgrid <wrt> of mass flux convection scheme',  &
                          'm kg kg-1 s-1', ymasks )
  call Les_diachro_write( tpdiafile, XLES_SUBGRID_WTHVMF,   'WTHV_MF',  'Subgrid <wthv> of mass flux convection scheme', &
                          'm K s-1', ymasks )
  call Les_diachro_write( tpdiafile, XLES_SUBGRID_WUMF,     'WU_MF',    'Subgrid <wu> of mass flux convection scheme',   &
                          'm2 s-2', ymasks )
  call Les_diachro_write( tpdiafile, XLES_SUBGRID_WVMF,     'WV_MF',    'Subgrid <wv> of mass flux convection scheme',   &
                          'm2 s-2', ymasks )

  call Les_diachro_write( tpdiafile, XLES_SUBGRID_PHI3,  'SBG_PHI3', 'Subgrid Phi3 function',         '1',      ymasks )
  call Les_diachro_write( tpdiafile, XLES_SUBGRID_LMix,  'SBG_LMIX', 'Subgrid Mixing Length',         '1',      ymasks )
  call Les_diachro_write( tpdiafile, XLES_SUBGRID_LDiss, 'SBG_LDIS', 'Subgrid Dissipation Length',    '1',      ymasks )
  call Les_diachro_write( tpdiafile, XLES_SUBGRID_Km,    'SBG_KM',   'Eddy diffusivity for momentum', 'm2 s-1', ymasks )
  call Les_diachro_write( tpdiafile, XLES_SUBGRID_Kh,    'SBG_KH',   'Eddy diffusivity for heat',     'm2 s-1', ymasks )

  if ( luserv ) then
    call Les_diachro_write( tpdiafile, XLES_SUBGRID_WThv,  'SBG_WTHV', 'Subgrid vertical flux of liquid potential temperature', &
                            'm K s-1', ymasks )
    call Les_diachro_write( tpdiafile, XLES_SUBGRID_Rt2,   'SBG_RT2',  'Subgrid total water variance', 'kg2 kg-2', ymasks )
    call Les_diachro_write( tpdiafile, XLES_SUBGRID_ThlRt, 'SBG_TLRT', 'Subgrid <thlrt> covariance',  'K kg kg-1', ymasks )
    call Les_diachro_write( tpdiafile, XLES_SUBGRID_URt,   'SBG_URT',  'Subgrid total water horizontal flux', &
                            'm kg kg-1 s-1', ymasks )
    call Les_diachro_write( tpdiafile, XLES_SUBGRID_VRt,   'SBG_VRT',  'Subgrid total water horizontal flux', &
                            'm kg kg-1 s-1', ymasks )
    call Les_diachro_write( tpdiafile, XLES_SUBGRID_WRt,   'SBG_WRT',  'Subgrid total water vertical flux',   &
                            'm kg kg-1 s-1', ymasks )
    call Les_diachro_write( tpdiafile, XLES_SUBGRID_PSI3,  'SBG_PSI3', 'Subgrid Psi3 function', '1', ymasks )
  end if

  if ( luserc ) then
    call Les_diachro_write( tpdiafile, XLES_SUBGRID_Rc2, 'SBG_RC2', 'Subgrid cloud water variance',        'kg2 kg-2', ymasks )
    call Les_diachro_write( tpdiafile, XLES_SUBGRID_URc, 'SBG_URC', 'Subgrid cloud water horizontal flux', 'm kg kg-1 s-1', &
                            ymasks )
    call Les_diachro_write( tpdiafile, XLES_SUBGRID_VRc, 'SBG_VRC', 'Subgrid cloud water horizontal flux', 'm kg kg-1 s-1', &
                            ymasks )
    call Les_diachro_write( tpdiafile, XLES_SUBGRID_WRc, 'SBG_WRC', 'Subgrid cloud water vertical flux',   'm kg kg-1 s-1', &
                            ymasks )
  end if

  if ( nsv > 0 ) then
    tfield%ndims = 4
    tfield%ndimlist(1)  = NMNHDIM_BUDGET_LES_LEVEL
    tfield%ndimlist(2)  = NMNHDIM_BUDGET_LES_TIME
    tfield%ndimlist(3)  = NMNHDIM_BUDGET_LES_MASK
    tfield%ndimlist(4)  = NMNHDIM_BUDGET_LES_SV
    tfield%ndimlist(5:) = NMNHDIM_UNUSED

    call Les_diachro_write( tpdiafile, XLES_SUBGRID_USv, 'SBG_USV', 'Subgrid <uSv> horizontal flux', 'm kg kg-1 s-1', ymasks )
    call Les_diachro_write( tpdiafile, XLES_SUBGRID_VSv, 'SBG_VSV', 'Subgrid <vSv> horizontal flux', 'm kg kg-1 s-1', ymasks )
    call Les_diachro_write( tpdiafile, XLES_SUBGRID_WSv, 'SBG_WSV', 'Subgrid <wSv> vertical flux',   'm kg kg-1 s-1', ymasks )

    tfield%ndims = 3
    !tfield%ndimlist(1)  = NMNHDIM_BUDGET_LES_LEVEL
    !tfield%ndimlist(2)  = NMNHDIM_BUDGET_LES_TIME
    !tfield%ndimlist(3)  = NMNHDIM_BUDGET_LES_MASK
    tfield%ndimlist(4)  = NMNHDIM_UNUSED
    !tfield%ndimlist(5:) = NMNHDIM_UNUSED


  end if

  call Les_diachro_write( tpdiafile, XLES_SUBGRID_UTke,  'SBG_UTKE', 'Subgrid flux of subgrid kinetic energy', 'm3 s-3',   ymasks )
  call Les_diachro_write( tpdiafile, XLES_SUBGRID_VTke,  'SBG_VTKE', 'Subgrid flux of subgrid kinetic energy', 'm3 s-3',   ymasks )
  call Les_diachro_write( tpdiafile, XLES_SUBGRID_WTke,  'SBG_WTKE', 'Subgrid flux of subgrid kinetic energy', 'm3 s-3',   ymasks )
  call Les_diachro_write( tpdiafile, XLES_SUBGRID_W2Thl, 'SBG_W2TL', 'Subgrid flux of subgrid kinetic energy', 'm2 K s-2', ymasks )
  call Les_diachro_write( tpdiafile, XLES_SUBGRID_WThl2, 'SBG_WTL2', 'Subgrid flux of subgrid kinetic energy', 'm K2 s-1', ymasks )
end if


!Prepare metadata (used in Les_diachro_write calls)
tfield%ndims = 2
tfield%ndimlist(1)  = NMNHDIM_BUDGET_LES_LEVEL
tfield%ndimlist(2)  = NMNHDIM_BUDGET_LES_TIME
tfield%ndimlist(3:) = NMNHDIM_UNUSED

ldoavg  = xles_temp_mean_start /= XUNDEF .and. xles_temp_mean_end /= XUNDEF
ldonorm = trim(cles_norm_type) /= 'NONE'
!
!*      2.4  Updraft quantities
!            ------------------
!
if ( lles_updraft ) then
  cgroup = 'Updraft'
  cgroupcomment = 'Updraft vertical profiles of some resolved and subgrid fluxes, variances and covariances'

  call Les_diachro_write( tpdiafile, XLES_UPDRAFT,     'UP_FRAC', 'Updraft fraction',                                 '1' )
  call Les_diachro_write( tpdiafile, XLES_UPDRAFT_W,   'UP_W',    'Updraft W mean value',                             'm s-1' )
  call Les_diachro_write( tpdiafile, XLES_UPDRAFT_Th,  'UP_TH',   'Updraft potential temperature mean value',         'K' )
  if ( luserc ) &
  call Les_diachro_write( tpdiafile, XLES_UPDRAFT_Thl, 'UP_THL',  'Updraft liquid potential temperature mean value',  'K' )
  if ( luserv ) &
  call Les_diachro_write( tpdiafile, XLES_UPDRAFT_Thv, 'UP_THV',  'Updraft virtual potential temperature mean value', 'K' )
  call Les_diachro_write( tpdiafile, XLES_UPDRAFT_Ke,  'UP_KE',   'Updraft resolved TKE mean value',                   'm2 s-2' )
  call Les_diachro_write( tpdiafile, XLES_UPDRAFT_Tke, 'UP_TKE',  'Updraft subgrid TKE mean value',                    'm2 s-2' )
  if ( luserv ) &
  call Les_diachro_write( tpdiafile, XLES_UPDRAFT_Rv,  'UP_RV',   'Updraft water vapor mean value',                    'kg kg-1' )
  if ( luserc ) &
  call Les_diachro_write( tpdiafile, XLES_UPDRAFT_Rc,  'UP_RC',   'Updraft cloud water mean value',                    'kg kg-1' )
  if ( luserr ) &
  call Les_diachro_write( tpdiafile, XLES_UPDRAFT_Rr,  'UP_RR',   'Updraft rain mean value',                           'kg kg-1' )
  if ( luseri ) &
  call Les_diachro_write( tpdiafile, XLES_UPDRAFT_Ri,  'UP_RI',   'Updraft ice mean value',                            'kg kg-1' )
  if ( lusers ) &
  call Les_diachro_write( tpdiafile, XLES_UPDRAFT_Rs,  'UP_RS',   'Updraft snow mean value',                           'kg kg-1' )
  if ( luserg ) &
  call Les_diachro_write( tpdiafile, XLES_UPDRAFT_Rg,  'UP_RG',   'Updraft graupel mean value',                        'kg kg-1' )
  if ( luserh ) &
  call Les_diachro_write( tpdiafile, XLES_UPDRAFT_Rh,  'UP_RH',   'Updraft hail mean value',                           'kg kg-1' )

  if ( nsv > 0 ) then
    tfield%ndims = 3
    tfield%ndimlist(1)  = NMNHDIM_BUDGET_LES_LEVEL
    tfield%ndimlist(2)  = NMNHDIM_BUDGET_LES_TIME
    tfield%ndimlist(3)  = NMNHDIM_BUDGET_LES_SV
    tfield%ndimlist(4:) = NMNHDIM_UNUSED

    call Les_diachro_write( tpdiafile, XLES_UPDRAFT_Sv, 'UP_SV', 'Updraft scalar variables mean values', 'kg kg-1' )

    tfield%ndims = 2
    !tfield%ndimlist(1)  = NMNHDIM_BUDGET_LES_LEVEL
    !tfield%ndimlist(2)  = NMNHDIM_BUDGET_LES_TIME
    tfield%ndimlist(3) = NMNHDIM_UNUSED
    !tfield%ndimlist(4:) = NMNHDIM_UNUSED
  end if

  call Les_diachro_write( tpdiafile, XLES_UPDRAFT_Th2,    'UP_TH2',  'Updraft resolved Theta variance',             'K2' )
  if ( luserc ) &
  call Les_diachro_write( tpdiafile, XLES_UPDRAFT_Thl2,   'UP_THL2', 'Updraft resolved Theta_l variance',           'K2' )
  if ( luserv ) &
  call Les_diachro_write( tpdiafile, XLES_UPDRAFT_ThThv,  'UP_THTV', 'Updraft resolved Theta Theta_v covariance',   'K2' )
  if ( luserc ) &
  call Les_diachro_write( tpdiafile, XLES_UPDRAFT_ThlThv, 'UP_TLTV', 'Updraft resolved Theta_l Theta_v covariance', 'K2' )
  call Les_diachro_write( tpdiafile, XLES_UPDRAFT_WTh,    'UP_WTH',  'Updraft resolved WTh flux',                    'm K s-1' )
  if ( luserc ) &
  call Les_diachro_write( tpdiafile, XLES_UPDRAFT_WThl,   'UP_WTHL', 'Updraft resolved WThl flux',                   'm K s-1' )
  if ( luserv ) &
  call Les_diachro_write( tpdiafile, XLES_UPDRAFT_WThv,   'UP_WTHV', 'Updraft resolved WThv flux',                   'm K s-1' )

  if ( luserv ) then
    call Les_diachro_write( tpdiafile, XLES_UPDRAFT_Rv2,   'UP_RV2',   'Updraft resolved water vapor variance', 'kg2 kg-2' )
    call Les_diachro_write( tpdiafile, XLES_UPDRAFT_ThRv,  'UP_THRV',  'Updraft resolved <thrv> covariance',    'K kg kg-1' )
    if ( luserc ) &
    call Les_diachro_write( tpdiafile, XLES_UPDRAFT_ThlRv, 'UP_THLRV', 'Updraft resolved <thlrv> covariance',   'K kg kg-1' )
    call Les_diachro_write( tpdiafile, XLES_UPDRAFT_ThvRv, 'UP_THVRV', 'Updraft resolved <thvrv> covariance',   'K kg kg-1' )
    call Les_diachro_write( tpdiafile, XLES_UPDRAFT_WRv,   'UP_WRV',   'Updraft resolved <wrv> vertical flux',  'm kg kg-1 s-1' )
  end if

  if ( luserc ) then
    call Les_diachro_write( tpdiafile, XLES_UPDRAFT_Rc2,   'UP_RC2',   'Updraft resolved cloud water variance', 'kg2 kg-2' )
    call Les_diachro_write( tpdiafile, XLES_UPDRAFT_ThRc,  'UP_THRC',  'Updraft resolved <thrc> covariance',    'K kg kg-1' )
    call Les_diachro_write( tpdiafile, XLES_UPDRAFT_ThlRc, 'UP_THLRC', 'Updraft resolved <thlrc> covariance',   'K kg kg-1' )
    call Les_diachro_write( tpdiafile, XLES_UPDRAFT_ThvRc, 'UP_THVRC', 'Updraft resolved <thvrc> covariance',   'K kg kg-1' )
    call Les_diachro_write( tpdiafile, XLES_UPDRAFT_WRc,   'UP_WRC',   'Updraft resolved <wrc> vertical flux',  'm kg kg-1 s-1' )
  end if

  if ( luseri ) then
    call Les_diachro_write( tpdiafile, XLES_UPDRAFT_Ri2,   'UP_RI2',   'Updraft resolved cloud ice variance',   'kg2 kg-2' )
    call Les_diachro_write( tpdiafile, XLES_UPDRAFT_ThRi,  'UP_THRI',  'Updraft resolved <thri> covariance',    'K kg kg-1' )
    call Les_diachro_write( tpdiafile, XLES_UPDRAFT_ThlRi, 'UP_THLRI', 'Updraft resolved <thlri> covariance',   'K kg kg-1' )
    call Les_diachro_write( tpdiafile, XLES_UPDRAFT_ThvRi, 'UP_THVRI', 'Updraft resolved <thvri> covariance',   'K kg kg-1' )
    call Les_diachro_write( tpdiafile, XLES_UPDRAFT_WRi,   'UP_WRI',   'Updraft resolved <wri> vertical flux',  'm kg kg-1 s-1' )
  end if


  if ( nsv > 0 ) then
    tfield%ndims = 3
    tfield%ndimlist(1)  = NMNHDIM_BUDGET_LES_LEVEL
    tfield%ndimlist(2)  = NMNHDIM_BUDGET_LES_TIME
    tfield%ndimlist(3)  = NMNHDIM_BUDGET_LES_SV
    tfield%ndimlist(4:) = NMNHDIM_UNUSED

    call Les_diachro_write( tpdiafile, XLES_UPDRAFT_Sv2,   'UP_SV2',   'Updraft resolved scalar variables variances', 'kg2 kg-2' )
    call Les_diachro_write( tpdiafile, XLES_UPDRAFT_ThSv,  'UP_THSV',  'Updraft resolved <ThSv> variance',            'K kg kg-1' )
    if ( luserc ) &
    call Les_diachro_write( tpdiafile, XLES_UPDRAFT_ThlSv, 'UP_THLSV', 'Updraft resolved <ThlSv> variance',           'K kg kg-1' )
    if ( luserv ) &
    call Les_diachro_write( tpdiafile, XLES_UPDRAFT_ThvSv, 'UP_THVSV', 'Updraft resolved <ThvSv> variance',           'K kg kg-1' )
    call Les_diachro_write( tpdiafile, XLES_UPDRAFT_WSv,   'UP_WSV',   'Updraft resolved <wSv> vertical flux', 'm kg kg-1 s-1' )

    tfield%ndims = 2
    !tfield%ndimlist(1)  = NMNHDIM_BUDGET_LES_LEVEL
    !tfield%ndimlist(2)  = NMNHDIM_BUDGET_LES_TIME
    tfield%ndimlist(3) = NMNHDIM_UNUSED
    !tfield%ndimlist(4:) = NMNHDIM_UNUSED
  end if
end if
!
!
!*      2.5  Downdraft quantities
!            --------------------
!
if ( lles_downdraft ) then
  cgroup = 'Downdraft'
  cgroupcomment = 'Downdraft vertical profiles of some resolved and subgrid fluxes, variances and covariances'

  call Les_diachro_write( tpdiafile, XLES_DOWNDRAFT,     'DW_FRAC', 'Downdraft fraction',                                 '1' )
  call Les_diachro_write( tpdiafile, XLES_DOWNDRAFT_W,   'DW_W',    'Downdraft W mean value',                             'm s-1' )
  call Les_diachro_write( tpdiafile, XLES_DOWNDRAFT_Th,  'DW_TH',   'Downdraft potential temperature mean value',         'K' )
  if ( luserc ) &
  call Les_diachro_write( tpdiafile, XLES_DOWNDRAFT_Thl, 'DW_THL',  'Downdraft liquid potential temperature mean value',  'K' )
  if ( luserv ) &
  call Les_diachro_write( tpdiafile, XLES_DOWNDRAFT_Thv, 'DW_THV',  'Downdraft virtual potential temperature mean value', 'K' )
  call Les_diachro_write( tpdiafile, XLES_DOWNDRAFT_Ke,  'DW_KE',   'Downdraft resolved TKE mean value',                 'm2 s-2' )
  call Les_diachro_write( tpdiafile, XLES_DOWNDRAFT_Tke, 'DW_TKE',  'Downdraft subgrid TKE mean value',                  'm2 s-2' )
  if ( luserv ) &
  call Les_diachro_write( tpdiafile, XLES_DOWNDRAFT_Rv,  'DW_RV',   'Downdraft water vapor mean value',                 'kg kg-1' )
  if ( luserc ) &
  call Les_diachro_write( tpdiafile, XLES_DOWNDRAFT_Rc,  'DW_RC',   'Downdraft cloud water mean value',                 'kg kg-1' )
  if ( luserr ) &
  call Les_diachro_write( tpdiafile, XLES_DOWNDRAFT_Rr,  'DW_RR',   'Downdraft rain mean value',                        'kg kg-1' )
  if ( luseri ) &
  call Les_diachro_write( tpdiafile, XLES_DOWNDRAFT_Ri,  'DW_RI',   'Downdraft ice mean value',                         'kg kg-1' )
  if ( lusers ) &
  call Les_diachro_write( tpdiafile, XLES_DOWNDRAFT_Rs,  'DW_RS',   'Downdraft snow mean value',                        'kg kg-1' )
  if ( luserg ) &
  call Les_diachro_write( tpdiafile, XLES_DOWNDRAFT_Rg,  'DW_RG',   'Downdraft graupel mean value',                     'kg kg-1' )
  if ( luserh ) &
  call Les_diachro_write( tpdiafile, XLES_DOWNDRAFT_Rh,  'DW_RH',   'Downdraft hail mean value',                        'kg kg-1' )

  if ( nsv > 0 ) then
    tfield%ndims = 3
    tfield%ndimlist(1)  = NMNHDIM_BUDGET_LES_LEVEL
    tfield%ndimlist(2)  = NMNHDIM_BUDGET_LES_TIME
    tfield%ndimlist(3)  = NMNHDIM_BUDGET_LES_SV
    tfield%ndimlist(4:) = NMNHDIM_UNUSED

    call Les_diachro_write( tpdiafile, XLES_DOWNDRAFT_Sv, 'DW_SV', 'Downdraft scalar variables mean values', 'kg kg-1' )

    tfield%ndims = 2
    !tfield%ndimlist(1)  = NMNHDIM_BUDGET_LES_LEVEL
    !tfield%ndimlist(2)  = NMNHDIM_BUDGET_LES_TIME
    tfield%ndimlist(3) = NMNHDIM_UNUSED
    !tfield%ndimlist(4:) = NMNHDIM_UNUSED
  end if

  call Les_diachro_write( tpdiafile, XLES_DOWNDRAFT_Th2,    'DW_TH2',  'Downdraft resolved Theta variance',             'K2' )
  if ( luserc ) &
  call Les_diachro_write( tpdiafile, XLES_DOWNDRAFT_Thl2,   'DW_THL2', 'Downdraft resolved Theta_l variance',           'K2' )
  if ( luserv ) &
  call Les_diachro_write( tpdiafile, XLES_DOWNDRAFT_ThThv,  'DW_THTV', 'Downdraft resolved Theta Theta_v covariance',   'K2' )
  if ( luserc ) &
  call Les_diachro_write( tpdiafile, XLES_DOWNDRAFT_ThlThv, 'DW_TLTV', 'Downdraft resolved Theta_l Theta_v covariance', 'K2' )
  call Les_diachro_write( tpdiafile, XLES_DOWNDRAFT_WTh,    'DW_WTH',  'Downdraft resolved WTh flux',                   'm K s-1' )
  if ( luserc ) &
  call Les_diachro_write( tpdiafile, XLES_DOWNDRAFT_WThl,   'DW_WTHL', 'Downdraft resolved WThl flux',                  'm K s-1' )
  if ( luserv ) &
  call Les_diachro_write( tpdiafile, XLES_DOWNDRAFT_WThv,   'DW_WTHV', 'Downdraft resolved WThv flux',                  'm K s-1' )

  if ( luserv ) then
    call Les_diachro_write( tpdiafile, XLES_DOWNDRAFT_Rv2,   'DW_RV2',   'Downdraft resolved water vapor variance', 'kg2 kg-2' )
    call Les_diachro_write( tpdiafile, XLES_DOWNDRAFT_ThRv,  'DW_THRV',  'Downdraft resolved <thrv> covariance',    'K kg kg-1' )
    if ( luserc ) &
    call Les_diachro_write( tpdiafile, XLES_DOWNDRAFT_ThlRv, 'DW_THLRV', 'Downdraft resolved <thlrv> covariance',   'K kg kg-1' )
    call Les_diachro_write( tpdiafile, XLES_DOWNDRAFT_ThvRv, 'DW_THVRV', 'Downdraft resolved <thvrv> covariance',   'K kg kg-1' )
    call Les_diachro_write( tpdiafile, XLES_DOWNDRAFT_WRv,   'DW_WRV',   'Downdraft resolved <wrv> vertical flux',  &
                            'm kg kg-1 s-1' )
  end if

  if ( luserc ) then
    call Les_diachro_write( tpdiafile, XLES_DOWNDRAFT_Rc2,   'DW_RC2',   'Downdraft resolved cloud water variance', 'kg2 kg-2' )
    call Les_diachro_write( tpdiafile, XLES_DOWNDRAFT_ThRc,  'DW_THRC',  'Downdraft resolved <thrc> covariance',    'K kg kg-1' )
    call Les_diachro_write( tpdiafile, XLES_DOWNDRAFT_ThlRc, 'DW_THLRC', 'Downdraft resolved <thlrc> covariance',   'K kg kg-1' )
    call Les_diachro_write( tpdiafile, XLES_DOWNDRAFT_ThvRc, 'DW_THVRC', 'Downdraft resolved <thvrc> covariance',   'K kg kg-1' )
    call Les_diachro_write( tpdiafile, XLES_DOWNDRAFT_WRc,   'DW_WRC',   'Downdraft resolved <wrc> vertical flux',  &
                            'm kg kg-1 s-1' )
  end if

  if ( luseri ) then
    call Les_diachro_write( tpdiafile, XLES_DOWNDRAFT_Ri2,   'DW_RI2',   'Downdraft resolved cloud ice variance',   'kg2 kg-2' )
    call Les_diachro_write( tpdiafile, XLES_DOWNDRAFT_ThRi,  'DW_THRI',  'Downdraft resolved <thri> covariance',    'K kg kg-1' )
    call Les_diachro_write( tpdiafile, XLES_DOWNDRAFT_ThlRi, 'DW_THLRI', 'Downdraft resolved <thlri> covariance',   'K kg kg-1' )
    call Les_diachro_write( tpdiafile, XLES_DOWNDRAFT_ThvRi, 'DW_THVRI', 'Downdraft resolved <thvri> covariance',   'K kg kg-1' )
    call Les_diachro_write( tpdiafile, XLES_DOWNDRAFT_WRi,   'DW_WRI',   'Downdraft resolved <wri> vertical flux',  &
                            'm kg kg-1 s-1' )
  end if


  if ( nsv > 0 ) then
    tfield%ndims = 3
    tfield%ndimlist(1)  = NMNHDIM_BUDGET_LES_LEVEL
    tfield%ndimlist(2)  = NMNHDIM_BUDGET_LES_TIME
    tfield%ndimlist(3)  = NMNHDIM_BUDGET_LES_SV
    tfield%ndimlist(4:) = NMNHDIM_UNUSED

    call Les_diachro_write( tpdiafile, XLES_DOWNDRAFT_Sv2,   'DW_SV2',   'Downdraft resolved scalar variables variances', &
                            'kg2 kg-2' )
    call Les_diachro_write( tpdiafile, XLES_DOWNDRAFT_ThSv,  'DW_THSV',  'Downdraft resolved <ThSv> variance',            &
                            'K kg kg-1' )
    if ( luserc ) &
    call Les_diachro_write( tpdiafile, XLES_DOWNDRAFT_ThlSv, 'DW_THLSV', 'Downdraft resolved <ThlSv> variance',           &
                            'K kg kg-1' )
    if ( luserv ) &
    call Les_diachro_write( tpdiafile, XLES_DOWNDRAFT_ThvSv, 'DW_THVSV', 'Downdraft resolved <ThvSv> variance',           &
                            'K kg kg-1' )
    call Les_diachro_write( tpdiafile, XLES_DOWNDRAFT_WSv,   'DW_WSV',   'Downdraft resolved <wSv> vertical flux',        &
                            'm kg kg-1 s-1' )

    tfield%ndims = 2
    !tfield%ndimlist(1)  = NMNHDIM_BUDGET_LES_LEVEL
    !tfield%ndimlist(2)  = NMNHDIM_BUDGET_LES_TIME
    tfield%ndimlist(3) = NMNHDIM_UNUSED
    !tfield%ndimlist(4:) = NMNHDIM_UNUSED
  end if
end if
!
!-------------------------------------------------------------------------------
!
!*      3.   surface normalization parameters
!            --------------------------------
!
cgroup = 'Radiation'
cgroupcomment = 'Radiative terms'

!Prepare metadata (used in Les_diachro_write calls)
tfield%ndims = 2
tfield%ndimlist(1)  = NMNHDIM_BUDGET_LES_LEVEL
tfield%ndimlist(2)  = NMNHDIM_BUDGET_LES_TIME
tfield%ndimlist(3:) = NMNHDIM_UNUSED

ldoavg  = xles_temp_mean_start /= XUNDEF .and. xles_temp_mean_end /= XUNDEF
ldonorm = .false.

call Les_diachro_write( tpdiafile, XLES_SWU,      'SWU',      'SW upward radiative flux',          'W m-2' )
call Les_diachro_write( tpdiafile, XLES_SWD,      'SWD',      'SW downward radiative flux',        'W m-2' )
call Les_diachro_write( tpdiafile, XLES_LWU,      'LWU',      'LW upward radiative flux',          'W m-2' )
call Les_diachro_write( tpdiafile, XLES_LWD,      'LWD',      'LW downward radiative flux',        'W m-2' )
call Les_diachro_write( tpdiafile, XLES_DTHRADSW, 'DTHRADSW', 'SW radiative temperature tendency', 'K s-1' )
call Les_diachro_write( tpdiafile, XLES_DTHRADLW, 'DTHRADLW', 'LW radiative temperature tendency', 'K s-1' )
!writes mean_effective radius at all levels
call Les_diachro_write( tpdiafile, XLES_RADEFF,   'RADEFF',   'Mean effective radius',             'micron' )


cgroup = 'Surface'
cgroupcomment = 'Averaged surface fields'

! !Prepare metadate (used in Les_diachro_write calls)
tfield%ndims = 1
tfield%ndimlist(1)  = NMNHDIM_BUDGET_LES_TIME
tfield%ndimlist(2:) = NMNHDIM_UNUSED

call Les_diachro_write( tpdiafile, XLES_Q0, 'Q0', 'Sensible heat flux at the surface', 'm K s-1' )
if ( luserv ) &
call Les_diachro_write( tpdiafile, XLES_E0, 'E0', 'Latent heat flux at the surface',   'kg kg-1 m s-1' )

if ( nsv > 0 ) then
  tfield%ndims = 2
  tfield%ndimlist(1)  = NMNHDIM_BUDGET_LES_TIME
  tfield%ndimlist(2)  = NMNHDIM_BUDGET_LES_SV
  tfield%ndimlist(3:) = NMNHDIM_UNUSED

  call Les_diachro_write( tpdiafile, XLES_SV0, 'SV0', 'Scalar variable fluxes at the surface', 'kg kg-1 m s-1' )

  tfield%ndims = 1
  !tfield%ndimlist(1)  = NMNHDIM_BUDGET_LES_TIME
  tfield%ndimlist(2)  = NMNHDIM_UNUSED
  !tfield%ndimlist(3:) = NMNHDIM_UNUSED
end if

call Les_diachro_write( tpdiafile, XLES_USTAR,      'Ustar',      'Friction velocity',                   'm s-1' )
call Les_diachro_write( tpdiafile, XLES_WSTAR,      'Wstar',      'Convective velocity',                 'm s-1' )
call Les_diachro_write( tpdiafile, XLES_MO_LENGTH,  'L_MO',       'Monin-Obukhov length',                'm' )
if ( luserr ) &
call Les_diachro_write( tpdiafile, XLES_PRECFR,     'PREC_FRAC',  'Fraction of columns where rain at surface',  '1' )
if ( luserr ) &
call Les_diachro_write( tpdiafile, XLES_INPRR,      'INST_PREC',  'Instantaneous precipitation rate',       'mm day-1' )
if ( luserc ) &
call Les_diachro_write( tpdiafile, XLES_INPRC,      'INST_SEDIM', 'Instantaneous cloud precipitation rate', 'mm day-1' )
if ( luserc .and. ( ldeposc .or. ldepoc ) ) &
call Les_diachro_write( tpdiafile, XLES_INDEP,      'INST_DEPOS', 'Instantaneous cloud deposition rate',    'mm day-1' )
if ( luserr ) &
call Les_diachro_write( tpdiafile, XLES_RAIN_INPRR, 'RAIN_PREC',  'Instantaneous precipitation rate over rainy grid cells', &
                        'mm day-1' )
if ( luserr ) &
call Les_diachro_write( tpdiafile, XLES_ACPRR,      'ACCU_PREC',  'Accumulated precipitation rate',             'mm' )


cgroup        = 'Miscellaneous'
cgroupcomment = 'Miscellaneous terms (geometry, various unclassified averaged terms...)'

call Les_diachro_write( tpdiafile, XLES_BL_HEIGHT,  'BL_H',       'Boundary Layer Height',               'm' )
call Les_diachro_write( tpdiafile, XLES_INT_TKE,    'INT_TKE',    'Vertical integrated TKE',             'm2 s-2' )
if ( luserc ) &
call Les_diachro_write( tpdiafile, XLES_ZCB,        'ZCB',        'Cloud base Height',                   'm' )
if ( luserc ) &
call Les_diachro_write( tpdiafile, XLES_CFtot,      'ZCFTOT',     'Total cloud cover (rc>1e-6)',         '1' )
if ( luserc ) &
call Les_diachro_write( tpdiafile, XLES_CF2tot,     'ZCF2TOT',    'Total cloud cover (rc>1e-5)',         '1' )
if ( luserc ) &
call Les_diachro_write( tpdiafile, XLES_LWP,        'LWP',        'Liquid Water path',                   'kg m-2' )
if ( luserc ) &
call Les_diachro_write( tpdiafile, XLES_LWPVAR,     'LWPVAR',     'Liquid Water path variance',          'kg m-4' )
if ( luserr ) &
call Les_diachro_write( tpdiafile, XLES_RWP,        'RWP',        'Rain Water path',                     'kg m-2' )
if ( luseri ) &
call Les_diachro_write( tpdiafile, XLES_IWP,        'IWP',        'Ice Water path',                      'kg m-2' )
if ( lusers ) &
call Les_diachro_write( tpdiafile, XLES_SWP,        'SWP',        'Snow Water path',                     'kg m-2' )
if ( luserg ) &
call Les_diachro_write( tpdiafile, XLES_GWP,        'GWP',        'Graupel Water path',                  'kg m-2' )
if ( luserh ) &
call Les_diachro_write( tpdiafile, XLES_HWP,        'HWP',        'Hail Water path',                     'kg m-2' )
if ( luserc ) &
call Les_diachro_write( tpdiafile, XLES_ZMAXCF,     'ZMAXCF',     'Height of Cloud fraction maximum (rc>1e-6)', 'm' )
if ( luserc ) &
call Les_diachro_write( tpdiafile, XLES_ZMAXCF2,    'ZMAXCF2',    'Height of Cloud fraction maximum (rc>1e-5)', 'm' )

!-------------------------------------------------------------------------------
!
!*      4.   LES budgets
!            -----------
!
call Write_les_budget_n( tpdiafile )

if ( luserv )  call Write_les_rt_budget_n( tpdiafile )

if ( nsv > 0 ) call Write_les_sv_budget_n( tpdiafile )
!
!-------------------------------------------------------------------------------
!
!*      5.   (ni,z,t) and (nj,z,t) 2points correlations
!            ------------------------------------------
!
if ( nspectra_k > 0 ) then
  tfieldx%cstdname = ''
  tfieldx%ngrid    = 0 !Not on the Arakawa grid
  tfieldx%ntype    = TYPEREAL
  tfieldx%ndims    = 3
  tfieldx%ndimlist(1)  = NMNHDIM_SPECTRA_2PTS_NI
  tfieldx%ndimlist(2)  = NMNHDIM_SPECTRA_LEVEL
  tfieldx%ndimlist(3)  = NMNHDIM_BUDGET_LES_TIME
  tfieldx%ndimlist(4:) = NMNHDIM_UNUSED

  tfieldy%cstdname = ''
  tfieldy%ngrid    = 0 !Not on the Arakawa grid
  tfieldy%ntype    = TYPEREAL
  tfieldy%ndims    = 3
  tfieldy%ndimlist(1)  = NMNHDIM_SPECTRA_2PTS_NJ
  tfieldy%ndimlist(2)  = NMNHDIM_SPECTRA_LEVEL
  tfieldy%ndimlist(3)  = NMNHDIM_BUDGET_LES_TIME
  tfieldy%ndimlist(4:) = NMNHDIM_UNUSED

  call Les_diachro_2pt_write( tpdiafile, XCORRi_UU, XCORRj_UU, 'UU', 'U*U     2 points correlations', 'm2 s-2' )
  call Les_diachro_2pt_write( tpdiafile, XCORRi_VV, XCORRj_VV, 'VV', 'V*V     2 points correlations', 'm2 s-2' )
  call Les_diachro_2pt_write( tpdiafile, XCORRi_WW, XCORRj_WW, 'WW', 'W*W     2 points correlations', 'm2 s-2' )
  call Les_diachro_2pt_write( tpdiafile, XCORRi_UV, XCORRj_UV, 'UV', 'U*V     2 points correlations', 'm2 s-2' )
  call Les_diachro_2pt_write( tpdiafile, XCORRi_WU, XCORRj_WU, 'WU', 'W*U     2 points correlations', 'm2 s-2' )
  call Les_diachro_2pt_write( tpdiafile, XCORRi_WV, XCORRj_WV, 'WV', 'W*V     2 points correlations', 'm2 s-2' )

  call Les_diachro_2pt_write( tpdiafile, XCORRi_ThTh, XCORRj_ThTh, 'THTH', 'Th*Th   2 points correlations', 'K2' )
  if ( luserc ) &
  call Les_diachro_2pt_write( tpdiafile, XCORRi_ThlThl, XCORRj_ThlThl, 'TLTL', 'Thl*Thl 2 points correlations', 'K2' )
  call Les_diachro_2pt_write( tpdiafile, XCORRi_WTh,    XCORRj_WTh,    'WTH',  'W*Th    2 points correlations', 'm K s-1' )
  if ( luserc ) &
  call Les_diachro_2pt_write( tpdiafile, XCORRi_WThl,   XCORRj_WThl,   'WTHL', 'W*Thl   2 points correlations', 'm K s-1' )

  if ( luserv ) then
    call Les_diachro_2pt_write( tpdiafile, XCORRi_RvRv,  XCORRj_RvRv,  'RVRV', 'rv*rv   2 points correlations', 'kg2 kg-2' )
    call Les_diachro_2pt_write( tpdiafile, XCORRi_ThRv,  XCORRj_ThRv,  'THRV', 'TH*RV   2 points correlations', 'K kg kg-1' )
    if ( luserc ) &
    call Les_diachro_2pt_write( tpdiafile, XCORRi_ThlRv, XCORRj_ThlRv, 'TLRV', 'thl*rv  2 points correlations', 'K kg kg-1' )
    call Les_diachro_2pt_write( tpdiafile, XCORRi_WRv,   XCORRj_WRv,   'WRV',  'W*rv    2 points correlations', 'm kg s-1 kg-1' )
  end if

  if ( luserc ) then
    call Les_diachro_2pt_write( tpdiafile, XCORRi_RcRc,  XCORRj_RcRc,  'RCRC', 'rc*rc   2 points correlations', 'kg2 kg-2' )
    call Les_diachro_2pt_write( tpdiafile, XCORRi_ThRc,  XCORRj_ThRc,  'THRC', 'th*rc   2 points correlations', 'K kg kg-1' )
    call Les_diachro_2pt_write( tpdiafile, XCORRi_ThlRc, XCORRj_ThlRc, 'TLRC', 'thl*rc  2 points correlations', 'K kg kg-1' )
    call Les_diachro_2pt_write( tpdiafile, XCORRi_WRc,   XCORRj_WRc,   'WRC',  'W*rc    2 points correlations', 'm kg s-1 kg-1' )
  end if

  if ( luseri ) then
    call Les_diachro_2pt_write( tpdiafile, XCORRi_RiRi,  XCORRj_RiRi,  'RIRI', 'ri*ri   2 points correlations', 'kg2 kg-2' )
    call Les_diachro_2pt_write( tpdiafile, XCORRi_ThRi,  XCORRj_ThRi,  'THRI', 'th*ri   2 points correlations', 'K kg kg-1' )
    call Les_diachro_2pt_write( tpdiafile, XCORRi_ThlRi, XCORRj_ThlRi, 'TLRI', 'thl*ri  2 points correlations', 'K kg kg-1' )
    call Les_diachro_2pt_write( tpdiafile, XCORRi_WRi,   XCORRj_WRi,   'WRI',  'W*ri    2 points correlations', 'm kg s-1 kg-1' )
  end if

!PW: TODO: ameliorer le ygroup (tenir compte de ce qu'est la variable scalaire et pas juste son jsv!)
  do jsv = 1, nsv
    Write( ygroup, fmt = "( a2, i3.3 )" ) "SS", jsv
    call Les_diachro_2pt_write( tpdiafile, XCORRi_SvSv(:,:,:,JSV), XCORRj_SvSv(:,:,:,JSV), ygroup, &
                                'Sv*Sv   2 points correlations','kg2 kg-2' )
  end do

!PW: TODO: ameliorer le ygroup (tenir compte de ce qu'est la variable scalaire et pas juste son jsv!)
  do jsv = 1, nsv
    Write( ygroup, fmt = "( a2, i3.3 )" ) "WS", jsv
    call Les_diachro_2pt_write( tpdiafile, XCORRi_WSv(:,:,:,JSV), XCORRj_WSv(:,:,:,JSV), ygroup, &
                                'W*Sv    2 points correlations','m kg s-1 kg-1' )
  end do
end if
!
!-------------------------------------------------------------------------------
!
!*      6.   spectra and time-averaged profiles (if first call to WRITE_LES_n)
!            ----------------------------------
!
call Les_spec_n( tpdiafile )
!
!-------------------------------------------------------------------------------
!
!*      7.   deallocations
!            -------------
!
CALL LES_DEALLOCATE('XLES_CURRENT_Z')

IF (CLES_NORM_TYPE/='NONE' ) THEN
  CALL LES_DEALLOCATE('XLES_NORM_M')
  CALL LES_DEALLOCATE('XLES_NORM_S')
  CALL LES_DEALLOCATE('XLES_NORM_K')
  CALL LES_DEALLOCATE('XLES_NORM_RHO')
  CALL LES_DEALLOCATE('XLES_NORM_RV')
  CALL LES_DEALLOCATE('XLES_NORM_SV')
  CALL LES_DEALLOCATE('XLES_NORM_P')
END IF

end subroutine Write_les_n

!------------------------------------------------------------------------------

subroutine Les_diachro_write_1D( tpdiafile, pdata, hmnhname, hcomment, hunits )

use modd_io,          only: tfiledata

use mode_les_diachro, only: Les_diachro

type(tfiledata),    intent(in) :: tpdiafile ! file to write
real, dimension(:), intent(in) :: pdata
character(len=*),   intent(in) :: hmnhname
character(len=*),   intent(in) :: hcomment
character(len=*),   intent(in) :: hunits

tfield%cmnhname  = hmnhname
tfield%clongname = hmnhname
tfield%ccomment  = hcomment
tfield%cunits    = hunits

call Les_diachro( tpdiafile, tfield, cgroup, cgroupcomment, ldoavg, ldonorm, pdata )

end subroutine Les_diachro_write_1D

!------------------------------------------------------------------------------

subroutine Les_diachro_write_2D( tpdiafile, pdata, hmnhname, hcomment, hunits )

use modd_io,          only: tfiledata

use mode_les_diachro, only: Les_diachro

type(tfiledata),      intent(in) :: tpdiafile ! file to write
real, dimension(:,:), intent(in) :: pdata
character(len=*),     intent(in) :: hmnhname
character(len=*),     intent(in) :: hcomment
character(len=*),     intent(in) :: hunits

tfield%cmnhname  = hmnhname
tfield%clongname = hmnhname
tfield%ccomment  = hcomment
tfield%cunits    = hunits

call Les_diachro( tpdiafile, tfield, cgroup, cgroupcomment, ldoavg, ldonorm, pdata )

end subroutine Les_diachro_write_2D

!------------------------------------------------------------------------------

subroutine Les_diachro_write_3D( tpdiafile, pdata, hmnhname, hcomment, hunits, hmasks )

use modd_io,          only: tfiledata

use mode_les_diachro, only: Les_diachro

type(tfiledata),                          intent(in) :: tpdiafile ! file to write
real,             dimension(:,:,:),       intent(in) :: pdata
character(len=*),                         intent(in) :: hmnhname
character(len=*),                         intent(in) :: hcomment
character(len=*),                         intent(in) :: hunits
character(len=*), dimension(:), optional, intent(in) :: hmasks

tfield%cmnhname  = hmnhname
tfield%clongname = hmnhname
tfield%ccomment  = hcomment
tfield%cunits    = hunits

call Les_diachro( tpdiafile, tfield, cgroup, cgroupcomment, ldoavg, ldonorm, pdata, hmasks = hmasks )

end subroutine Les_diachro_write_3D

!------------------------------------------------------------------------------

subroutine Les_diachro_write_4D( tpdiafile, pdata, hmnhname, hcomment, hunits, hmasks )

use modd_io,          only: tfiledata

use mode_les_diachro, only: Les_diachro

type(tfiledata),                            intent(in) :: tpdiafile ! file to write
real,             dimension(:,:,:,:),       intent(in) :: pdata
character(len=*),                           intent(in) :: hmnhname
character(len=*),                           intent(in) :: hcomment
character(len=*),                           intent(in) :: hunits
character(len=*), dimension(:),   optional, intent(in) :: hmasks

tfield%cmnhname  = hmnhname
tfield%clongname = hmnhname
tfield%ccomment  = hcomment
tfield%cunits    = hunits

call Les_diachro( tpdiafile, tfield, cgroup, cgroupcomment, ldoavg, ldonorm, pdata, hmasks = hmasks )

end subroutine Les_diachro_write_4D

!------------------------------------------------------------------------------

subroutine Les_diachro_2pt_write( tpdiafile, zcorri, zcorrj, hmnhname, hcomment, hunits )

use modd_io,          only: tfiledata

use mode_les_diachro, only: Les_diachro_2pt

type(tfiledata),          intent(in) :: tpdiafile ! file to write
real, dimension(:,:,:),   intent(in) :: zcorri    ! 2 pts correlation data
real, dimension(:,:,:),   intent(in) :: zcorrj    ! 2 pts correlation data
character(len=*),         intent(in) :: hmnhname
character(len=*),         intent(in) :: hcomment
character(len=*),         intent(in) :: hunits

tfieldx%cmnhname  = hmnhname
tfieldx%clongname = hmnhname
tfieldx%ccomment  = hcomment
tfieldx%cunits    = hunits

tfieldy%cmnhname  = hmnhname
tfieldy%clongname = hmnhname
tfieldy%ccomment  = hcomment
tfieldy%cunits    = hunits

call Les_diachro_2pt( tpdiafile, tfieldx, tfieldy, zcorri, zcorrj )

end subroutine Les_diachro_2pt_write

!------------------------------------------------------------------------------

end module mode_write_les_n
