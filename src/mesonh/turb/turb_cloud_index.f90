!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!    ################
     MODULE MODI_TURB_CLOUD_INDEX
!    ################
!
INTERFACE
!
      SUBROUTINE TURB_CLOUD_INDEX(PTSTEP,TPFILE,                            &
                                  OTURB_DIAG,KRRI,                          &
                                  PRRS,PRM,PRHODJ,PDXX,PDYY,PDZZ,PDZX,PDZY, &
                                  PCEI                                      )
!
USE MODD_IO, ONLY: TFILEDATA
!
REAL,                   INTENT(IN)   ::  PTSTEP       ! Double Time step
TYPE(TFILEDATA),        INTENT(IN)   ::  TPFILE       ! Output file
LOGICAL,                INTENT(IN)   ::  OTURB_DIAG   ! switch to write some
                                 ! diagnostic fields in the syncronous FM-file
INTEGER,                INTENT(IN)   :: KRRI          ! number of ice water var.
REAL, DIMENSION(:,:,:,:), INTENT(IN) :: PRRS ! Sources term of RR
REAL, DIMENSION(:,:,:,:), INTENT(IN) :: PRM  ! Variable at t-dt
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PRHODJ  ! Jacobian * dry density of
                                                  !  the reference state
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PDXX,PDYY,PDZZ,PDZX,PDZY
                                                  !  metric coefficients
REAL, DIMENSION(:,:,:),   INTENT(OUT):: PCEI ! Cloud Entrainment instability
                                             ! index to emphasize locally
                                             ! turbulent fluxes
!
END SUBROUTINE TURB_CLOUD_INDEX
!
END INTERFACE
!
END MODULE MODI_TURB_CLOUD_INDEX
!
!     #######################
      SUBROUTINE TURB_CLOUD_INDEX(PTSTEP,TPFILE,                            &
                                  OTURB_DIAG,KRRI,                          &
                                  PRRS,PRM,PRHODJ,PDXX,PDYY,PDZZ,PDZX,PDZY, &
                                  PCEI                                      )
!     #######################
      !
!!    PURPOSE
!!    -------
!!    CEI (cloud Entrainment Instability) index calculation
!!    It permits to localize cloudy points where a different mixing length
!!       from the one in clear sky can be applicated
!!    It permits to quantify also, at those cloudy points, an instability
!!       that can emphasize sub-grid turbulence.
!!    If such an instability exists, mixing length is increased proportionnaly
!!       to that CEI criterium
!!
!!**  METHOD
!!    ------
!!
!!    Criteria: For a cloudy point or a point adjacent to a cloudy point,
!!              G   = NORM( dVAR/dx_j ) > threshold
!!              Q_j = DG_j/Dt   of the same sign as G_j
!!              where VAR=rv+rc+ri and j=x or y
!!              then CEI= NORM(Q)
!!
!!    EXTERNAL
!!    --------
!!      GX_M_M, GY_M_M : Cartesian gradient operators
!!      FMWRIT         : FM-routine to write a record
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODI_GRADIENT_M : GX_M_M, GY_M_M
!!
!!    AUTHOR
!!    ------
!!      M. Tomasini                  * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    15/09/94
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!
!-------------------------------------------------------------------------------
!
use modd_field,          only: tfielddata, TYPEREAL
USE MODD_IO,             ONLY: TFILEDATA
USE MODD_PARAMETERS,     ONLY: JPVEXT
!
USE MODE_IO_FIELD_WRITE, only: IO_Field_write
use mode_tools_ll,       only: GET_INDICE_ll
!
USE MODI_GRADIENT_M
!
IMPLICIT NONE
!
!*       0.     DECLARATIONS
!               ------------
!
!*       0.1   declarations of arguments
!
REAL,                   INTENT(IN)   ::  PTSTEP       ! Double Time step
TYPE(TFILEDATA),        INTENT(IN)   ::  TPFILE       ! Output file
LOGICAL,                INTENT(IN)   ::  OTURB_DIAG   ! switch to write some
                                 ! diagnostic fields in the syncronous FM-file
INTEGER,                INTENT(IN)   :: KRRI          ! number of ice water var.
REAL, DIMENSION(:,:,:,:), INTENT(IN) :: PRRS ! Sources term of RR
REAL, DIMENSION(:,:,:,:), INTENT(IN) :: PRM  ! Variable at t-dt
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PRHODJ  ! Jacobian * dry density of
                                                  !  the reference state
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PDXX,PDYY,PDZZ,PDZX,PDZY
                                                  !  metric coefficients
REAL, DIMENSION(:,:,:),   INTENT(OUT):: PCEI ! Cloud Entrainment instability
                                             ! index to emphasize locally
                                             ! turbulent fluxes
!
!*       0.2   declarations of local variables
!
REAL, DIMENSION(SIZE(PRM,1),SIZE(PRM,2),SIZE(PRM,3)) :: ZWORK,ZRVCI0 ! Work arrays
REAL, DIMENSION(SIZE(PRM,1),SIZE(PRM,2),SIZE(PRM,3)) :: ZCLOUD
                             ! rc+ri at time after ADVECTION routine
                             ! for the CEI criterium
REAL, DIMENSION(SIZE(PRM,1),SIZE(PRM,2),SIZE(PRM,3)) :: ZRVCI,ZGNORM_RVCI,ZQNORM_RVCI
                                  ! rv+rc+ri at time after ADVECTION routine
                                  ! horizontal norm of the vector PG_RVCI
                                  ! horizontal norm of the vector PQ_RVCI
REAL, DIMENSION(SIZE(PRM,1),SIZE(PRM,2),SIZE(PRM,3),2) :: ZG_RVCI,ZQ_RVCI
                                  ! x and y gradient of rv+rc+ri
                                  ! x and y gradient of the advection of rv+rc+ri
!
INTEGER             :: JI,JJ,JK     ! loop counters
INTEGER             :: IIB,IJB,IKB  ! Begin of physical dimensions
INTEGER             :: IIE,IJE,IKE  ! End   of physical dimensions
INTEGER, DIMENSION(SIZE(PRM,1),SIZE(PRM,2),SIZE(PRM,3)) :: IMASK_CLOUD
                             ! 0 except cloudy points or adjacent points (1)
TYPE(TFIELDDATA)    :: TZFIELD
!
!-------------------------------------------------------------------------------
!
!*       1.     INITIALISATION
!               --------------
!
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
IKB = 1   + JPVEXT
IKE = SIZE(PRM,3) - JPVEXT
!
IMASK_CLOUD(:,:,:) = 0
PCEI(:,:,:) = 0.
!
!-------------------------------------------------------------------------------
!
!*       2.     CALCULATION
!               -----------
!*       2.1    Gradients calculation of the variable :
!               VAR at time (t+1)=VAR at time (t-1) + 2*dt*ADV at time t
!               VAR is a source term (i.e. x by RHODJ)
!
! To avoid negative mixing ratios at external points
!                                 but also in the physical domain !
ZRVCI0(:,:,:) = MAX ( PRRS(:,:,:,1) , 0. ) + MAX ( PRRS(:,:,:,2) , 0. )
IF (KRRI>=1) ZRVCI0(:,:,:) = ZRVCI0(:,:,:) + MAX ( PRRS(:,:,:,4) , 0. )
!
ZRVCI(:,:,:)= PTSTEP *ZRVCI0(:,:,:) /PRHODJ(:,:,:)
ZG_RVCI(:,:,:,1) = GX_M_M(ZRVCI,PDXX,PDZZ,PDZX)
ZG_RVCI(:,:,:,2) = GY_M_M(ZRVCI,PDYY,PDZZ,PDZY)
!
ZGNORM_RVCI(:,:,:) = SQRT( ZG_RVCI(:,:,:,1)*ZG_RVCI(:,:,:,1) +        &
                           ZG_RVCI(:,:,:,2)*ZG_RVCI(:,:,:,2) )
!
!
!*       2.2    Frontogenetic terms calculation
!               (gradient of the advection)
!               Q_j=DG_j/Dt=d(DVAR/Dt)dx_j - d(u_k*dVAR/dx_k)/dx_j
!               As DVAR/Dt=0 if the VAR is conserved during the movement,
!               Q_j = dADV/dx_j
!               VAR=rv+rc+ri
!
ZWORK(:,:,:) = ZRVCI0 / PRHODJ(:,:,:) -   &
               ( PRM(:,:,:,1)+ PRM(:,:,:,2) ) / PTSTEP
IF (KRRI>=1) ZWORK(:,:,:) = ZWORK(:,:,:) - PRM(:,:,:,4) / PTSTEP
!
ZQ_RVCI(:,:,:,1) = GX_M_M(ZWORK,PDXX,PDZZ,PDZX)
ZQ_RVCI(:,:,:,2) = GY_M_M(ZWORK,PDYY,PDZZ,PDZY)
!
ZQNORM_RVCI(:,:,:) = SQRT( ZQ_RVCI(:,:,:,1)*ZQ_RVCI(:,:,:,1) +         &
                           ZQ_RVCI(:,:,:,2)*ZQ_RVCI(:,:,:,2) )
!
!
!*       2.3    Cloud mask
!
ZCLOUD(:,:,:)= MAX ( PRRS(:,:,:,2) , 0. )
IF (KRRI>=1) ZCLOUD(:,:,:) = ZCLOUD(:,:,:) + MAX ( PRRS(:,:,:,4) , 0. )
ZCLOUD(:,:,:) = PTSTEP * ZCLOUD / PRHODJ(:,:,:)
!
DO JK=IKB,IKE
DO JJ=IJB,IJE
DO JI=IIB,IIE
   ! rc+ri threshold to avoid white noise and calculations
   IF ( ZCLOUD(JI,JJ,JK) > 1.E-6 ) THEN
      IMASK_CLOUD(JI-1,JJ  ,JK  ) = 1
      IMASK_CLOUD(JI  ,JJ  ,JK  ) = 1
      IMASK_CLOUD(JI+1,JJ  ,JK  ) = 1
      IMASK_CLOUD(JI  ,JJ-1,JK  ) = 1
      IMASK_CLOUD(JI  ,JJ+1,JK  ) = 1
      IMASK_CLOUD(JI  ,JJ  ,JK-1) = 1
      IMASK_CLOUD(JI  ,JJ  ,JK+1) = 1
      ! The cloudy points where the criteria will not be satisfied
      ! will have the cloudy mixing length not amplified
      ! We put in the CEI index a negative number to mark those points
      ! in turb.f90
      PCEI(JI,JJ,JK) = -1.
   ENDIF
ENDDO
ENDDO
ENDDO
!
!*       2.4    Cloud Entrainment Instability index
!
! CEI(:,:,:)=NORM_Q
!
!   if the considered point is cloudy or surrounded by at least one cloudy point
!
!   and if the characteristic time >0 in at least one direction that is to say
!          |grad(rv+rc+ri)| increasing with time that is to say
!           grad(rv+rc+ri) has the same sign as Q_RVCI
!
!   and if NORM_G_RVCI >= 0.1 g/kg/km
!
DO JK=IKB,IKE
DO JJ=IJB,IJE
DO JI=IIB,IIE
   IF ( IMASK_CLOUD(JI,JJ,JK) == 1 ) THEN
      IF ( ZGNORM_RVCI(JI,JJ,JK) >= 1.E-07 ) THEN
         IF ( SIGN(1.0,ZG_RVCI(JI,JJ,JK,1))==SIGN(1.0,ZQ_RVCI(JI,JJ,JK,1)) .OR.   &
              SIGN(1.0,ZG_RVCI(JI,JJ,JK,2))==SIGN(1.0,ZQ_RVCI(JI,JJ,JK,2)) ) THEN
            PCEI(JI,JJ,JK) = ZQNORM_RVCI(JI,JJ,JK)
         ENDIF
      ENDIF
   ENDIF
ENDDO
ENDDO
ENDDO
!
!*       2.5    Writing
!
IF ( OTURB_DIAG .AND. tpfile%lopened ) THEN
  TZFIELD%CMNHNAME   = 'RVCI'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'RVCI'
  TZFIELD%CUNITS     = 'kg kg-1'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Z_RVCI'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZRVCI)
  !
  TZFIELD%CMNHNAME   = 'GX_RVCI'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'GX_RVCI'
  TZFIELD%CUNITS     = 'kg kg-1 m-1'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Z_GX_RVCI'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZG_RVCI(:,:,:,1))
  !
  TZFIELD%CMNHNAME   = 'GY_RVCI'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'GY_RVCI'
  TZFIELD%CUNITS     = 'kg kg-1 m-1'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Z_GY_RVCI'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZG_RVCI(:,:,:,2))
  !
  TZFIELD%CMNHNAME   = 'GNORM_RVCI'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'GNORM_RVCI'
  TZFIELD%CUNITS     = 'kg kg-1 m-1'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Z_NORM G'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZGNORM_RVCI)
  !
  TZFIELD%CMNHNAME   = 'QX_RVCI'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'QX_RVCI'
  TZFIELD%CUNITS     = 'kg kg-1 m-1'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Z_QX_RVCI'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZQ_RVCI(:,:,:,1))
  !
  TZFIELD%CMNHNAME   = 'QY_RVCI'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'QY_RVCI'
  TZFIELD%CUNITS     = 'kg kg-1 m-1'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Z_QY_RVCI'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZQ_RVCI(:,:,:,2))
  !
  TZFIELD%CMNHNAME   = 'QNORM_RVCI'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'QNORM_RVCI'
  TZFIELD%CUNITS     = 'kg kg-1 m-1'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Z_QNORM_RVCI'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZQNORM_RVCI)
  !
  TZFIELD%CMNHNAME   = 'CEI'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'CEI'
  TZFIELD%CUNITS     = 'kg kg-1 m-1 s-1'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Z_CEI'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,PCEI)
END IF
!
END SUBROUTINE TURB_CLOUD_INDEX
