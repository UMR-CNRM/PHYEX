!     ######spl
      MODULE MODI_CH_CONVECT_LINOX
!     ############################
!
INTERFACE
!
      SUBROUTINE CH_CONVECT_LINOX( KLON, KLEV, PCH1, PCH1C,            &
                                   KDPL, KPBL, KLCL, KCTL, KLFS, KDBL, &
                                   PUMF, PUER, PUDR, PDMF, PDER, PDDR, &
                                   PTIMEC, PDXDY, PMIXF, PLMASS, PWSUB,&
                                   KFTSTEPS, PUTT, PRHODREF,           &
                                   OUSECHEM, PZZ, PIC_RATE, PCG_RATE   )
!
INTEGER,                INTENT(IN) :: KLON     ! horizontal dimension
INTEGER,                INTENT(IN) :: KLEV     ! vertical dimension
!
REAL,DIMENSION(KLON,KLEV),INTENT(IN)   :: PCH1 ! grid scale tracer concentr.
REAL,DIMENSION(KLON,KLEV),INTENT(OUT)  :: PCH1C! conv adjusted tracer concntr.
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
REAL, DIMENSION(KLON,KLEV),   INTENT(IN) :: PUTT      ! updraft temperature (K)
REAL, DIMENSION(KLON,KLEV),   INTENT(IN) :: PRHODREF  
!
LOGICAL,                      INTENT(IN) :: OUSECHEM ! to indicate if chemistry is used
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PZZ      ! height of model layer (m)
REAL, DIMENSION(KLON),     INTENT(INOUT) :: PIC_RATE ! IC lightning frequency
REAL, DIMENSION(KLON),     INTENT(INOUT) :: PCG_RATE ! CG lightning frequency
!
END SUBROUTINE CH_CONVECT_LINOX
!
END INTERFACE
!
END MODULE MODI_CH_CONVECT_LINOX
