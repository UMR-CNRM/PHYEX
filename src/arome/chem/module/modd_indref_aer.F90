!     ######spl
     MODULE MODD_INDREF_AER
!    ########################
!
!!****  *MODD_INDREF_AER* - declaration of indices of refraction for the
!                           different aerosol types
!!
!!    PURPOSE
!!    -------
!!       1D Arrays to store the real and the imaginary parts of the indices
!        of refraction for each of the 6 SW bands for every type of aerosol
!        that can enter into the cloud droplet composition
!!
!-------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
!------------------------------
IMPLICIT NONE
!
! the wavelenghts corresponding to the 6 SW bands
!lambda = 0.217, 0.345, 0.565, 0.94, 1.78, 3.19 micron
!
REAL, DIMENSION(6):: XLAMBDA=(/0.217,0.345,0.565,0.94,1.78,3.19/)
!
!real and imaginary parts of the refractive indice of water(Hale and Querry
!1973)for the 6 SW bands:
REAL,DIMENSION(6)::XN_W=(/1.382, 1.3436, 1.333, 1.3274,&
                                   1.3125, 1.479/)
REAL,DIMENSION(6)::XK_W=(/7.E-08, 7.3E-09, 2.944E-09, 2.182E-06,&
                                   1.1205E-04, 0.10092/) 
!                                   
!!real and imaginary parts of the refractive indice for the 6 SW bands for
!differents types of aerosols inclusions into the droplet:
! first indice for the SW band, second for the aerosol type
!aer=1: dust like (d'almeida,1992, table 4.3) 
!aer=2: waso = water soluble (d'almeida,1992, table 4.3) 
!aer=3: soot (d'almeida,1992, table 4.3) 
!aer=4: sesa = sea salt (OPAC) 
!aer=5: sulf = sulfates (d'almeida,1992, table 4.3) 
!
!
REAL,DIMENSION(6)::XN_DUST=(/1.53, 1.53, 1.53, 1.52, 1.340, 1.22/)
REAL,DIMENSION(6)::XK_DUST=(/8.E-03, 8.E-03, 8.E-03, 8.E-03, 8.E-03, 1.E-02/)
!
REAL,DIMENSION(6)::XN_WATSOL=(/1.53, 1.53, 1.53, 1.52, 1.42, 1.43/)
REAL,DIMENSION(6)::XK_WATSOL=(/8.E-03, 5.E-03, 6.E-03, 1.31E-02, 2.5E-02, 8.E-03/) 
!
REAL,DIMENSION(6)::XN_SOOT=(/1.62, 1.75, 1.75, 1.754, 1.8, 1.859/)
REAL,DIMENSION(6)::XK_SOOT=(/0.45, 0.465, 0.44, 0.4354, 0.481 , 0.54 /) 
!
REAL,DIMENSION(6)::XN_SEASALT=(/1.51, 1.51, 1.5, 1.475, 1.45, 1.49/)
REAL,DIMENSION(6)::XK_SEASALT=(/5.E-06, 3.24E-07, 1.E-08, 9.E-05, 7.62E-04, 3.E-03/) 
!
REAL,DIMENSION(6)::XN_SULFATES=(/1.484, 1.452, 1.43, 1.4235, 1.394, 1.311/)
REAL,DIMENSION(6)::XK_SULFATES=(/1.E-08, 1.E-08, 1.15E-08, 7.08E-07, 6.35E-04, 1.35E-01/) 
!
END MODULE MODD_INDREF_AER
