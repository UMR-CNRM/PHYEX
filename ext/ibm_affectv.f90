!MNH_LIC Copyright 2019-2022 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!
!    #######################
MODULE MODI_IBM_AFFECTV  
  !    ####################### 
  !
  INTERFACE
     !
     SUBROUTINE IBM_AFFECTV(PVAR,PVAR2,PVAR3,HVAR,KIBM_LAYER,HIBM_MODE_INTE3,&
          HIBM_FORC_BOUNR,PRADIUS,PPOWERS,&
          HIBM_MODE_INT1N,HIBM_TYPE_BOUNN,HIBM_MODE_BOUNN,HIBM_FORC_BOUNN,PIBM_FORC_BOUNN,&
          HIBM_MODE_INT1T,HIBM_TYPE_BOUNT,HIBM_MODE_BOUNT,HIBM_FORC_BOUNT,PIBM_FORC_BOUNT,&
          HIBM_MODE_INT1C,HIBM_TYPE_BOUNC,HIBM_MODE_BOUNC,HIBM_FORC_BOUNC,PIBM_FORC_BOUNC,PXMU,PDIV)
       !
       REAL, DIMENSION(:,:,:)   ,INTENT(INOUT) :: PVAR
       REAL, DIMENSION(:,:,:,:) ,INTENT(IN)    :: PVAR2,PVAR3
       CHARACTER(LEN=1)         ,INTENT(IN)    :: HVAR
       INTEGER                  ,INTENT(IN)    :: KIBM_LAYER
       REAL                     ,INTENT(IN)    :: PRADIUS
       REAL                     ,INTENT(IN)    :: PPOWERS
       CHARACTER(LEN=3)         ,INTENT(IN)    :: HIBM_FORC_BOUNR
       CHARACTER(LEN=3)         ,INTENT(IN)    :: HIBM_MODE_INTE3
       CHARACTER(LEN=3)         ,INTENT(IN)    :: HIBM_MODE_INT1N
       CHARACTER(LEN=3)         ,INTENT(IN)    :: HIBM_TYPE_BOUNN
       CHARACTER(LEN=3)         ,INTENT(IN)    :: HIBM_MODE_BOUNN 
       CHARACTER(LEN=3)         ,INTENT(IN)    :: HIBM_FORC_BOUNN
       REAL                     ,INTENT(IN)    :: PIBM_FORC_BOUNN
       CHARACTER(LEN=3)         ,INTENT(IN)    :: HIBM_MODE_INT1T
       CHARACTER(LEN=3)         ,INTENT(IN)    :: HIBM_TYPE_BOUNT
       CHARACTER(LEN=3)         ,INTENT(IN)    :: HIBM_MODE_BOUNT
       CHARACTER(LEN=3)         ,INTENT(IN)    :: HIBM_FORC_BOUNT
       REAL                     ,INTENT(IN)    :: PIBM_FORC_BOUNT
       CHARACTER(LEN=3)         ,INTENT(IN)    :: HIBM_MODE_INT1C
       CHARACTER(LEN=3)         ,INTENT(IN)    :: HIBM_TYPE_BOUNC
       CHARACTER(LEN=3)         ,INTENT(IN)    :: HIBM_MODE_BOUNC
       CHARACTER(LEN=3)         ,INTENT(IN)    :: HIBM_FORC_BOUNC
       REAL                     ,INTENT(IN)    :: PIBM_FORC_BOUNC
       REAL, DIMENSION(:,:,:)   ,INTENT(IN)    :: PXMU
       REAL, DIMENSION(:,:,:)   ,INTENT(IN)    :: PDIV
       !
     END SUBROUTINE IBM_AFFECTV
     !
  END INTERFACE
  !
END MODULE MODI_IBM_AFFECTV
!
!     ########################################################
SUBROUTINE IBM_AFFECTV(PVAR,PVAR2,PVAR3,HVAR,KIBM_LAYER,HIBM_MODE_INTE3,&
     HIBM_FORC_BOUNR,PRADIUS,PPOWERS,&
     HIBM_MODE_INT1N,HIBM_TYPE_BOUNN,HIBM_MODE_BOUNN,HIBM_FORC_BOUNN,PIBM_FORC_BOUNN,&
     HIBM_MODE_INT1T,HIBM_TYPE_BOUNT,HIBM_MODE_BOUNT,HIBM_FORC_BOUNT,PIBM_FORC_BOUNT,&
     HIBM_MODE_INT1C,HIBM_TYPE_BOUNC,HIBM_MODE_BOUNC,HIBM_FORC_BOUNC,PIBM_FORC_BOUNC,PXMU,PDIV)
  !     ########################################################
  !
  !
  !****  IBM_AFFECTV computes the variable PVAR on desired ghost points :
  !                 - the V type of the ghost/image 
  !                 - the 3D interpolation mode (HIBM_MODE_INTE3)
  !                 - the 1D interpolation mode (HIBM_MODE_INTE1)
  !                 - the boundary condition    (HIBM_TYPE_BOUND)  
  !                 - the symmetry character    (HIBM_MODE_BOUND)
  !                 - the forcing type          (HIBM_FORC_BOUND)    
  !                 - the forcing term          (HIBM_FORC_BOUND)
  !       Choice of forcing type is depending on 
  !                 the normal, binormal, tangent vectors (N,C,T)
  !                   
  !                
  !    PURPOSE
  !    -------
  !****  Ghosts (resp. Images) locations are stored in KIBM_STOR_GHOST (resp. KIBM_STOR_IMAGE).
  !      Solutions are computed in regard of the symmetry character of the solution:
  !                                      HIBM_MODE_BOUND = 'SYM' (Symmetrical)
  !                                      HIBM_MODE_BOUND = 'ASY' (Anti-symmetrical)
  !      The ghost value is depending on the variable value at the interface:
  !                                      HIBM_TYPE_BOUND = "CST" (constant value)
  !                                      HIBM_TYPE_BOUND = "LAW" (wall models)
  !                                      HIBM_TYPE_BOUND = "LIN" (linear evolution, only IMAGE2 type)
  !                                      HIBM_TYPE_BOUND = "LOG" (logarithmic evol, only IMAGE2 type)
  !      Three 3D interpolations exists  HIBM_MODE_INTE3 = "IDW" (Inverse  Distance Weighting)
  !                                      HIBM_MODE_INTE3 = "MDW" (Modified Distance Weighting)
  !                                      HIBM_MODE_INTE3 = "LAG" (Trilinear Lagrange interp. )
  !      Three 1D interpolations exists  HIBM_MODE_INTE1 = "CL0" (Lagrange Polynomials - 1 points - MIRROR)
  !                                      HIBM_MODE_INTE1 = "CL1" (Lagrange Polynomials - 2 points - IMAGE1)
  !                                      HIBM_MODE_INTE1 = "CL2" (Lagrange Polynomials - 3 points - IMAGE2)
  !    METHOD
  !    ------
  !      - loop on ghosts
  !      - functions storage
  !      - computations of the location of the corners cell containing MIRROR/IMAGE1/IMAGE2
  !      - 3D interpolation (IDW, MDW, CLI)  to obtain the MIRROR/IMAGE1/IMAGE2 values
  !      - computation of the value at the interface 
  !      - 1D interpolation (CLI1,CLI2,CLI3) to obtain the GHOSTS values
  !      - Affectation
  !
  !    EXTERNAL
  !    --------
  !      SUBROUTINE ?
  !
  !    IMPLICIT ARGUMENTS
  !    ------------------
  !       MODD_?   
  !
  !    REFERENCE
  !    ---------
  !
  !    AUTHOR
  !    ------
  !      Franck Auguste (CERFACS-AE)
  !
  !    MODIFICATIONS
  !    -------------
  !      Original         01/01/2019
  !
  !------------------------------------------------------------------------------
  !       
  !**** 0. DECLARATIONS
  !     ---------------
  ! module
  USE MODE_POS
  USE MODE_ll
  USE MODE_IO
  USE MODD_ARGSLIST_ll, ONLY : LIST_ll
  !
  ! declaration
  USE MODD_IBM_PARAM_n
  USE MODD_FIELD_n
  USE MODD_PARAM_n, ONLY: CTURB
  USE MODD_GRID_n,  ONLY: XDXHAT, XDYHAT
  USE MODD_VAR_ll,  ONLY: IP
  USE MODD_LBC_n
  USE MODD_REF_n, ONLY: XRHODJ,XRHODREF
  !
  ! interface
  USE MODI_IBM_VALUECORN
  USE MODI_IBM_LOCATCORN
  USE MODI_IBM_3DINT
  USE MODI_IBM_1DINT
  USE MODI_IBM_0DINT
  USE MODI_IBM_VALUEMAT1
  USE MODI_IBM_VALUEMAT2
  USE MODI_SHUMAN
  USE MODD_DYN_n
  USE MODD_FIELD_n
  USE MODD_CST
  USE MODD_CTURB
  USE MODD_RADIATIONS_n
  !
  IMPLICIT NONE
  !
  !------------------------------------------------------------------------------
  !
  !       0.1  declarations of arguments
  !
  REAL, DIMENSION(:,:,:)   ,INTENT(INOUT) :: PVAR
  REAL, DIMENSION(:,:,:,:) ,INTENT(IN)    :: PVAR2,PVAR3
  CHARACTER(LEN=1)         ,INTENT(IN)    :: HVAR
  INTEGER                  ,INTENT(IN)    :: KIBM_LAYER
  REAL                     ,INTENT(IN)    :: PRADIUS
  REAL                     ,INTENT(IN)    :: PPOWERS
  CHARACTER(LEN=3)         ,INTENT(IN)    :: HIBM_FORC_BOUNR
  CHARACTER(LEN=3)         ,INTENT(IN)    :: HIBM_MODE_INTE3
  CHARACTER(LEN=3)         ,INTENT(IN)    :: HIBM_MODE_INT1N
  CHARACTER(LEN=3)         ,INTENT(IN)    :: HIBM_TYPE_BOUNN
  CHARACTER(LEN=3)         ,INTENT(IN)    :: HIBM_MODE_BOUNN
  CHARACTER(LEN=3)         ,INTENT(IN)    :: HIBM_FORC_BOUNN
  REAL                     ,INTENT(IN)    :: PIBM_FORC_BOUNN
  CHARACTER(LEN=3)         ,INTENT(IN)    :: HIBM_MODE_INT1T
  CHARACTER(LEN=3)         ,INTENT(IN)    :: HIBM_TYPE_BOUNT
  CHARACTER(LEN=3)         ,INTENT(IN)    :: HIBM_MODE_BOUNT
  CHARACTER(LEN=3)         ,INTENT(IN)    :: HIBM_FORC_BOUNT
  REAL                     ,INTENT(IN)    :: PIBM_FORC_BOUNT
  CHARACTER(LEN=3)         ,INTENT(IN)    :: HIBM_MODE_INT1C
  CHARACTER(LEN=3)         ,INTENT(IN)    :: HIBM_TYPE_BOUNC
  CHARACTER(LEN=3)         ,INTENT(IN)    :: HIBM_MODE_BOUNC
  CHARACTER(LEN=3)         ,INTENT(IN)    :: HIBM_FORC_BOUNC
  REAL                     ,INTENT(IN)    :: PIBM_FORC_BOUNC
  REAL, DIMENSION(:,:,:)   ,INTENT(IN)    :: PXMU
  REAL, DIMENSION(:,:,:)   ,INTENT(IN)    :: PDIV
  !
  !------------------------------------------------------------------------------
  !
  !       0.2  declaration of local variables
  !
  INTEGER                                    :: JI,JJ,JK,JL,JM,JMM,JN,JNN,JH,JLL         ! loop index
  INTEGER, DIMENSION(:)  , ALLOCATABLE       :: I_INDEX_CORN                             ! reference corner index 
  INTEGER                                    :: I_GHOST_NUMB                             ! ghost number per layer
  REAL   , DIMENSION(:,:), ALLOCATABLE       :: Z_LOCAT_CORN,Z_LOCAT_IMAG                ! corners coordinates
  REAL   , DIMENSION(:)  , ALLOCATABLE       :: Z_TESTS_CORN                             ! interface distance dependence 
  REAL   , DIMENSION(:)  , ALLOCATABLE       :: Z_VALUE_CORN                             ! value variables at corners
  REAL   , DIMENSION(:,:), ALLOCATABLE       :: Z_VALUE_IMAG,Z_VALUE_TEMP,Z_VALUE_ZLKE   ! value at mirror/image1/image2 
  REAL   , DIMENSION(:)  , ALLOCATABLE       :: Z_LOCAT_BOUN,Z_LOCAT_GHOS,Z_TEMP_ZLKE    ! location of bound and ghost
  REAL                                       :: Z_DELTA_IMAG,ZIBM_VISC,ZIBM_DIVK
  CHARACTER(LEN=3),DIMENSION(:), ALLOCATABLE :: Y_TYPE_BOUND,Y_FORC_BOUND,Y_MODE_BOUND,Y_MODE_INTE1
  REAL   , DIMENSION(:)  , ALLOCATABLE       :: Z_FORC_BOUND,Z_VALUE_GHOS
  REAL   , DIMENSION(:,:), ALLOCATABLE       :: Z_VALUE_MAT1,Z_VALUE_MAT2
  REAL                                       :: ZIBM_HALO
  !
  !------------------------------------------------------------------------------
  !
  !       0.3  Allocation
  !
  ALLOCATE(I_INDEX_CORN(3))
  ALLOCATE(Z_LOCAT_CORN(8,3))
  ALLOCATE(Z_VALUE_CORN(8))
  ALLOCATE(Z_TESTS_CORN(8))
  ALLOCATE(Z_LOCAT_IMAG(3,3))
  ALLOCATE(Z_VALUE_IMAG(4,3))
  ALLOCATE(Z_VALUE_TEMP(4,3))
  ALLOCATE(Z_LOCAT_BOUN(3))
  ALLOCATE(Z_LOCAT_GHOS(3))
  ALLOCATE(Z_VALUE_GHOS(3))
  ALLOCATE(Y_TYPE_BOUND(3),Y_FORC_BOUND(3))
  ALLOCATE(Y_MODE_BOUND(3),Y_MODE_INTE1(3))
  ALLOCATE(Z_FORC_BOUND(3))
  ALLOCATE(Z_VALUE_MAT1(3,3))
  ALLOCATE(Z_VALUE_MAT2(3,3))
  !
  !------------------------------------------------------------------------------
  !
  !**** 1. PRELIMINARIES
  !     ----------------
  I_INDEX_CORN(:)    = 0
  Z_LOCAT_CORN(:,:)  = 0.
  Z_VALUE_CORN(:)    = 0.
  Z_TESTS_CORN(:)    = 0.
  Z_LOCAT_IMAG(:,:)  = 0.
  Z_VALUE_IMAG(:,:)  = 0.
  Z_VALUE_TEMP(:,:)  = 0.
  Z_LOCAT_GHOS(:)    = 0.
  Z_LOCAT_BOUN(:)    = 0.
  Z_VALUE_GHOS(:)    = 0.
  Z_VALUE_MAT1(:,:)  = 0.
  Z_VALUE_MAT2(:,:)  = 0.
  IF (HVAR=='U') JH = 1
  IF (HVAR=='V') JH = 2
  IF (HVAR=='W') JH = 3
  Y_TYPE_BOUND(1) = HIBM_TYPE_BOUNN
  Y_TYPE_BOUND(2) = HIBM_TYPE_BOUNT
  Y_TYPE_BOUND(3) = HIBM_TYPE_BOUNC
  Y_FORC_BOUND(1) = HIBM_FORC_BOUNN
  Y_FORC_BOUND(2) = HIBM_FORC_BOUNT
  Y_FORC_BOUND(3) = HIBM_FORC_BOUNC
  Y_MODE_BOUND(1) = HIBM_MODE_BOUNN
  Y_MODE_BOUND(2) = HIBM_MODE_BOUNT
  Y_MODE_BOUND(3) = HIBM_MODE_BOUNC
  Y_MODE_INTE1(1) = HIBM_MODE_INT1N
  Y_MODE_INTE1(2) = HIBM_MODE_INT1T
  Y_MODE_INTE1(3) = HIBM_MODE_INT1C
  Z_FORC_BOUND(1) = PIBM_FORC_BOUNN
  Z_FORC_BOUND(2) = PIBM_FORC_BOUNT
  Z_FORC_BOUND(3) = PIBM_FORC_BOUNC
  !
  ALLOCATE(Z_VALUE_ZLKE(4,3))
  ALLOCATE(Z_TEMP_ZLKE(3))
  Z_VALUE_ZLKE(:,:) = 0.
  Z_TEMP_ZLKE(:)  = 0.
  !
  DO JMM=1,KIBM_LAYER
     !
     ! searching number of ghosts 
     JM = size(NIBM_GHOST_V,1)
     JI = 0
     JJ = 0
     JK = 0
     DO WHILE ((JI==0.and.JJ==0.and.JK==0).and.JM>0)
        JI = NIBM_GHOST_V(JM,JMM,JH,1)
        JJ = NIBM_GHOST_V(JM,JMM,JH,2)
        JK = NIBM_GHOST_V(JM,JMM,JH,3)
        IF (JI==0.and.JJ==0.and.JK==0) JM = JM - 1
     ENDDO
     I_GHOST_NUMB = JM
     !
     ! Loop on each P Ghosts 
     IF (I_GHOST_NUMB<=0) GO TO 666
     DO JM = 1,I_GHOST_NUMB
        !
        ! ghost index/ls
        JI = NIBM_GHOST_V(JM,JMM,JH,1)
        JJ = NIBM_GHOST_V(JM,JMM,JH,2)
        JK = NIBM_GHOST_V(JM,JMM,JH,3)
        IF (JI==0.or.JJ==0.or.JK==0) GO TO 777 
        Z_LOCAT_GHOS(:) =  XIBM_GHOST_V(JM,JMM,JH,:)
        Z_LOCAT_BOUN(:) = 2.0*XIBM_IMAGE_V(JM,JMM,JH,1,:)-1.0*XIBM_IMAGE_V(JM,JMM,JH,2,:) 
        ZIBM_HALO = 1.
        !
        DO JN = 1,3
           !
           Z_LOCAT_IMAG(JN,:)= XIBM_IMAGE_V(JM,JMM,JH  ,JN,:)
           Z_DELTA_IMAG      = ( XDXHAT(JI) * XDYHAT(JJ) ) ** 0.5
           !
           DO JLL=1,3
              I_INDEX_CORN(:)       = NIBM_IMAGE_V(JM,JMM,JH,JLL,JN,:)
              IF (I_INDEX_CORN(1)==0.AND.JN==2) ZIBM_HALO=0.
              IF (I_INDEX_CORN(2)==0.AND.JN==2) ZIBM_HALO=0.
              Z_LOCAT_CORN(:,:)     = IBM_LOCATCORN(I_INDEX_CORN,JLL+1)
              Z_TESTS_CORN(:)       = XIBM_TESTI_V(JM,JMM,JH,JLL,JN,:)
              Z_VALUE_CORN(:)       = IBM_VALUECORN(PVAR2(:,:,:,JLL),I_INDEX_CORN)          
              Z_VALUE_IMAG(JN,JLL)  = IBM_3DINT(JN,Z_VALUE_IMAG(:,JLL),Z_LOCAT_BOUN,Z_TESTS_CORN,&
                   Z_LOCAT_CORN,Z_VALUE_CORN,Z_LOCAT_IMAG(JN,:),&
                   HIBM_MODE_INTE3,PRADIUS,PPOWERS)
           ENDDO
           !
        ENDDO
        ZIBM_VISC = PXMU(JI,JJ,JK)
        ZIBM_DIVK = PDIV(JI,JJ,JK)
        !
        ! projection step
        Z_VALUE_MAT1(:,:) = IBM_VALUEMAT1(Z_LOCAT_IMAG(1,:),Z_LOCAT_BOUN,Z_VALUE_IMAG,HIBM_FORC_BOUNR)
        DO JN=1,3 
           Z_VALUE_TEMP(JN,:)= Z_VALUE_MAT1(:,1)*Z_VALUE_IMAG(JN,1) +&
                Z_VALUE_MAT1(:,2)*Z_VALUE_IMAG(JN,2) +&
                Z_VALUE_MAT1(:,3)*Z_VALUE_IMAG(JN,3) 
        ENDDO
        !
        ! === BOUND computation ===
        !
        JN=4
        DO JLL=1,3
           Z_VALUE_TEMP(JN,JLL)  = IBM_0DINT(Z_DELTA_IMAG,Z_VALUE_TEMP(:,JLL),Y_TYPE_BOUND(JLL),Y_FORC_BOUND(JLL), & 
                Z_FORC_BOUND(JLL),ZIBM_VISC,ZIBM_DIVK)
        ENDDO
        !
        ! inverse projection step
        Z_VALUE_MAT2(:,:) = IBM_VALUEMAT2(Z_VALUE_MAT1)
        Z_VALUE_IMAG(JN,:)= Z_VALUE_MAT2(:,1)*Z_VALUE_TEMP(JN,1) +&
             Z_VALUE_MAT2(:,2)*Z_VALUE_TEMP(JN,2) +&
             Z_VALUE_MAT2(:,3)*Z_VALUE_TEMP(JN,3)
        !
        ! === GHOST computation ===
        !
        ! functions storage
        Z_LOCAT_IMAG(1,3) = ((XIBM_GHOST_V(JM,JMM,JH,1)-Z_LOCAT_BOUN(1))**2.+&
             (XIBM_GHOST_V(JM,JMM,JH,2)-Z_LOCAT_BOUN(2))**2.+&
             (XIBM_GHOST_V(JM,JMM,JH,3)-Z_LOCAT_BOUN(3))**2.)**0.5
        IF (Z_LOCAT_IMAG(1,3)>Z_DELTA_IMAG.AND.ZIBM_HALO>0.5) THEN
           Z_LOCAT_IMAG(1,1) = ((XIBM_IMAGE_V(JM,JMM,JH,1,1)-Z_LOCAT_BOUN(1))**2.+&
                (XIBM_IMAGE_V(JM,JMM,JH,1,2)-Z_LOCAT_BOUN(2))**2.+&
                (XIBM_IMAGE_V(JM,JMM,JH,1,3)-Z_LOCAT_BOUN(3))**2.)**0.5
           Z_LOCAT_IMAG(1,2) = ((XIBM_IMAGE_V(JM,JMM,JH,2,1)-Z_LOCAT_BOUN(1))**2.+&
                (XIBM_IMAGE_V(JM,JMM,JH,2,2)-Z_LOCAT_BOUN(2))**2.+&
                (XIBM_IMAGE_V(JM,JMM,JH,2,3)-Z_LOCAT_BOUN(3))**2.)**0.5
        ELSE
           Z_LOCAT_IMAG(1,1) = ((XIBM_IMAGE_V(JM,JMM,JH,3,1)-Z_LOCAT_BOUN(1))**2.+&
                (XIBM_IMAGE_V(JM,JMM,JH,3,2)-Z_LOCAT_BOUN(2))**2.+&
                (XIBM_IMAGE_V(JM,JMM,JH,3,3)-Z_LOCAT_BOUN(3))**2.)**0.5
           Z_LOCAT_IMAG(1,2) = ((XIBM_IMAGE_V(JM,JMM,JH,1,1)-Z_LOCAT_BOUN(1))**2.+&
                (XIBM_IMAGE_V(JM,JMM,JH,1,2)-Z_LOCAT_BOUN(2))**2.+&
                (XIBM_IMAGE_V(JM,JMM,JH,1,3)-Z_LOCAT_BOUN(3))**2.)**0.5
           Z_VALUE_TEMP(2,:) = Z_VALUE_TEMP(1,:)
           Z_VALUE_TEMP(1,:) = Z_VALUE_TEMP(3,:)
        ENDIF
        ! 
        DO JLL=1,3
           Z_VALUE_GHOS(JLL) = IBM_1DINT(Z_LOCAT_IMAG(1,:),Z_VALUE_TEMP(:,JLL),Y_MODE_INTE1(JLL))
           IF (Y_MODE_BOUND(JLL)=='SYM') Z_VALUE_GHOS(JLL) = +Z_VALUE_GHOS(JLL)
           IF (Y_MODE_BOUND(JLL)=='ASY') Z_VALUE_GHOS(JLL) = -Z_VALUE_GHOS(JLL) + 2.*Z_VALUE_TEMP(4,JLL)
           IF (Y_MODE_BOUND(JLL)=='CST') Z_VALUE_GHOS(JLL) = Z_VALUE_TEMP(4,JLL) 
        ENDDO
        !
        PVAR(JI,JJ,JK) =  Z_VALUE_MAT2(JH,1)*Z_VALUE_GHOS(1) +&
             Z_VALUE_MAT2(JH,2)*Z_VALUE_GHOS(2) +&
             Z_VALUE_MAT2(JH,3)*Z_VALUE_GHOS(3)
        !
        IF ((JH==3).AND.(JK==2)) THEN
           PVAR(JI,JJ,JK) = 0.
        ENDIF
        !
777     CONTINUE
        !
     ENDDO
  ENDDO
  !
666 CONTINUE
  !
  !**** X. DEALLOCATIONS/CLOSES
  !     -----------------------
  !
  DEALLOCATE(I_INDEX_CORN)
  DEALLOCATE(Z_LOCAT_CORN)
  DEALLOCATE(Z_VALUE_CORN)
  DEALLOCATE(Z_LOCAT_IMAG)
  DEALLOCATE(Z_VALUE_IMAG)
  DEALLOCATE(Z_VALUE_TEMP)
  DEALLOCATE(Z_LOCAT_BOUN)
  DEALLOCATE(Z_LOCAT_GHOS)
  DEALLOCATE(Z_VALUE_GHOS)
  DEALLOCATE(Z_TESTS_CORN)
  DEALLOCATE(Y_TYPE_BOUND,Y_FORC_BOUND)
  DEALLOCATE(Y_MODE_BOUND,Y_MODE_INTE1)
  DEALLOCATE(Z_FORC_BOUND)
  DEALLOCATE(Z_VALUE_MAT1)
  DEALLOCATE(Z_VALUE_MAT2)
  DEALLOCATE(Z_VALUE_ZLKE)
  DEALLOCATE(Z_TEMP_ZLKE)
  ! 
  RETURN
  !
END SUBROUTINE IBM_AFFECTV
