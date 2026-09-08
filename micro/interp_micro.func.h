!These routines are intented to be included in the contains part of other subroutines.
!To allow the transformation for GPU, no local array must be declared.
!If a temporary local array is needed, it must be added as a buffer in the interface (IBUF?, ZBUF?)

SUBROUTINE INTERP_MICRO_1D(D, PIN, KNUM, P1, P2, &
                           LDPACK, LDMASK, KBUF1, KBUF2, PBUF1, PBUF2, &
                           KLEN, &
                           PLT1, POUT1, PLT2, POUT2, PLT3, POUT3)

USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t

IMPLICIT NONE

TYPE(DIMPHYEX_t),           INTENT(IN)  :: D
REAL,    DIMENSION(D%NIJT), INTENT(IN)  :: PIN          !Input array
INTEGER,                    INTENT(IN)  :: KNUM         !Number of points in the look-up table
REAL,                       INTENT(IN)  :: P1           !Scaling factor
REAL,                       INTENT(IN)  :: P2           !Scaling factor
LOGICAL,                    INTENT(IN)  :: LDPACK       !.TRUE. to perform packing
LOGICAL, DIMENSION(D%NIJT), INTENT(IN)  :: LDMASK       !Computation mask
INTEGER, DIMENSION(D%NIJT), INTENT(OUT) :: KBUF1, KBUF2 !Buffer arrays
REAL,    DIMENSION(D%NIJT), INTENT(OUT) :: PBUF1, PBUF2 !Buffer arrays
INTEGER,                    INTENT(OUT) :: KLEN         !Number of active points
REAL,    DIMENSION(KNUM),   INTENT(IN)            :: PLT1  !Look-up table
REAL,    DIMENSION(D%NIJT), INTENT(OUT)           :: POUT1 !Interpolated values
REAL,    DIMENSION(KNUM),   INTENT(IN) , OPTIONAL :: PLT2
REAL,    DIMENSION(D%NIJT), INTENT(OUT), OPTIONAL :: POUT2
REAL,    DIMENSION(KNUM),   INTENT(IN) , OPTIONAL :: PLT3
REAL,    DIMENSION(D%NIJT), INTENT(OUT), OPTIONAL :: POUT3

INTEGER :: JIJ
INTEGER :: IINDEX
REAL :: ZINDEX

IF (LDPACK) THEN

  !Pack input array
  KLEN=0
  DO JIJ=D%NIJB, D%NIJE
    IF (LDMASK(JIJ)) THEN
      KLEN=KLEN+1
      PBUF1(KLEN)=PIN(JIJ)
      KBUF1(KLEN)=JIJ
    ENDIF
  ENDDO

  IF (KLEN>0) THEN
    !Index computation
    DO JIJ=1, KLEN
      PBUF1(JIJ) = MAX(1.00001, MIN(REAL(KNUM)-0.00001, P1 * LOG(PBUF1(JIJ)) + P2))
      KBUF2(JIJ) = INT(PBUF1(JIJ))
      PBUF1(JIJ) = PBUF1(JIJ) - REAL(KBUF2(JIJ))
    END DO

    !Interpolation and unpack
    DO JIJ=1, KLEN
      PBUF2(JIJ) = PLT1(KBUF2(JIJ)+1) *  PBUF1(JIJ)       &
                    &-PLT1(KBUF2(JIJ)  ) * (PBUF1(JIJ) - 1.0)
    END DO
    POUT1(:)=0.
    DO JIJ=1, KLEN
      POUT1(KBUF1(JIJ))=PBUF2(JIJ)
    ENDDO

    !Interpolation and unpack 2
    IF(PRESENT(PLT2)) THEN
      DO JIJ=1, KLEN
        PBUF2(JIJ) = PLT2(KBUF2(JIJ)+1) *  PBUF1(JIJ)       &
                      &-PLT2(KBUF2(JIJ)  ) * (PBUF1(JIJ) - 1.0)
      END DO
      POUT2(:)=0.
      DO JIJ=1, KLEN
        POUT2(KBUF1(JIJ))=PBUF2(JIJ)
      ENDDO
    ENDIF

    !Interpolation and unpack 3
    IF(PRESENT(PLT3)) THEN
      DO JIJ=1, KLEN
        PBUF2(JIJ) = PLT3(KBUF2(JIJ)+1) *  PBUF1(JIJ)       &
                      &-PLT3(KBUF2(JIJ)  ) * (PBUF1(JIJ) - 1.0)
      END DO
      POUT3(:)=0.
      DO JIJ=1, KLEN
        POUT3(KBUF1(JIJ))=PBUF2(JIJ)
      ENDDO
    ENDIF

  ENDIF

ELSE

  KLEN=0
  DO JIJ=D%NIJB, D%NIJE
    IF (LDMASK(JIJ)) THEN
      KLEN=KLEN+1

      !Index computation
      ZINDEX = MAX(1.00001, MIN(REAL(KNUM)-0.00001, P1 * LOG(PIN(JIJ)) + P2))
      IINDEX = INT(ZINDEX)
      ZINDEX = ZINDEX - REAL(IINDEX)

      !Interpolations
      POUT1(JIJ) = PLT1(IINDEX+1) *  ZINDEX       &
                &-PLT1(IINDEX  ) * (ZINDEX - 1.0)

      IF(PRESENT(PLT2)) THEN
        POUT2(JIJ) = PLT2(IINDEX+1) *  ZINDEX       &
                  &-PLT2(IINDEX  ) * (ZINDEX - 1.0)
      ENDIF

      IF(PRESENT(PLT3)) THEN
        POUT3(JIJ) = PLT3(IINDEX+1) *  ZINDEX       &
                  &-PLT3(IINDEX  ) * (ZINDEX - 1.0)
      ENDIF

    ELSE
      POUT1(JIJ) = 0.
      IF(PRESENT(PLT2)) POUT2(JIJ) = 0.
      IF(PRESENT(PLT3)) POUT3(JIJ) = 0.
    ENDIF
  ENDDO

ENDIF
END SUBROUTINE INTERP_MICRO_1D

SUBROUTINE INTERP_MICRO_2D(D, PIN1, PIN2, KNUM1, KNUM2, P11, P12, P21, P22,&
                           LDPACK, LDMASK, KBUF1, KBUF2, KBUF3, PBUF1, PBUF2, PBUF3, &
                           KLEN, &
                           PLT1, POUT1, PLT2, POUT2, PLT3, POUT3)

USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t

IMPLICIT NONE

TYPE(DIMPHYEX_t),           INTENT(IN)  :: D
REAL,    DIMENSION(D%NIJT), INTENT(IN)  :: PIN1                !Input array
REAL,    DIMENSION(D%NIJT), INTENT(IN)  :: PIN2                !Input array
INTEGER,                    INTENT(IN)  :: KNUM1               !First dimension of the look-up table
INTEGER,                    INTENT(IN)  :: KNUM2               !Second dimension of the look-up table
REAL,                       INTENT(IN)  :: P11                 !Scaling factor
REAL,                       INTENT(IN)  :: P12                 !Scaling factor
REAL,                       INTENT(IN)  :: P21                 !Scaling factor
REAL,                       INTENT(IN)  :: P22                 !Scaling factor
LOGICAL,                    INTENT(IN)  :: LDPACK              !.TRUE. to perform packing
LOGICAL, DIMENSION(D%NIJT), INTENT(IN)  :: LDMASK              !Computation mask
INTEGER, DIMENSION(D%NIJT), INTENT(OUT) :: KBUF1, KBUF2, KBUF3 !Buffer arrays
REAL,    DIMENSION(D%NIJT), INTENT(OUT) :: PBUF1, PBUF2, PBUF3 !Buffer arrays
INTEGER,                    INTENT(OUT) :: KLEN                !Number of active points
REAL,    DIMENSION(KNUM1, KNUM2),   INTENT(IN)            :: PLT1  !Look-up table
REAL,    DIMENSION(D%NIJT),         INTENT(OUT)           :: POUT1 !Interpolated values from the first look-up table
REAL,    DIMENSION(KNUM1, KNUM2),   INTENT(IN) , OPTIONAL :: PLT2  !Other look-up table
REAL,    DIMENSION(D%NIJT),         INTENT(OUT), OPTIONAL :: POUT2 !Interpolated values from the second look-up table
REAL,    DIMENSION(KNUM2, KNUM1),   INTENT(IN) , OPTIONAL :: PLT3  !Another look-up table **CAUTION, TABLE IS REVERSED**
REAL,    DIMENSION(D%NIJT),         INTENT(OUT), OPTIONAL :: POUT3 !Interpolated values from the third look-up table

INTEGER :: JIJ
INTEGER :: IINDEX1, IINDEX2
REAL :: ZINDEX1, ZINDEX2

IF (LDPACK) THEN

  !Pack input array
  KLEN=0
  DO JIJ=D%NIJB, D%NIJE
    IF (LDMASK(JIJ)) THEN
      KLEN=KLEN+1
      PBUF1(KLEN)=PIN1(JIJ)
      PBUF2(KLEN)=PIN2(JIJ)
      KBUF3(KLEN)=JIJ
    ENDIF
  ENDDO

  IF (KLEN>0) THEN
    !Index computation
    DO JIJ=1, KLEN
      PBUF1(JIJ) = MAX(1.00001, MIN(REAL(KNUM1)-0.00001, P11 * LOG(PBUF1(JIJ)) + P12))
      KBUF1(JIJ) = INT(PBUF1(JIJ))
      PBUF1(JIJ) = PBUF1(JIJ) - REAL(KBUF1(JIJ))

      PBUF2(JIJ) = MAX(1.00001, MIN(REAL(KNUM2)-0.00001, P21 * LOG(PBUF2(JIJ)) + P22))
      KBUF2(JIJ) = INT(PBUF2(JIJ))
      PBUF2(JIJ) = PBUF2(JIJ) - REAL(KBUF2(JIJ))
    END DO

    !Interpolation and unpack 1
    DO JIJ=1, KLEN
      PBUF3(JIJ) = ( PLT1(KBUF1(JIJ)+1,KBUF2(JIJ)+1)* PBUF2(JIJ)         &
                   -PLT1(KBUF1(JIJ)+1,KBUF2(JIJ)  )*(PBUF2(JIJ) - 1.0)) *  PBUF1(JIJ) &
                 -( PLT1(KBUF1(JIJ)  ,KBUF2(JIJ)+1)* PBUF2(JIJ)         &
                   -PLT1(KBUF1(JIJ)  ,KBUF2(JIJ)  )*(PBUF2(JIJ) - 1.0)) * (PBUF1(JIJ) - 1.0)
    ENDDO
    POUT1(:)=0.
    DO JIJ=1, KLEN
      POUT1(KBUF3(JIJ))=PBUF3(JIJ)
    ENDDO

    !Interpolation and unpack 2
    IF(PRESENT(PLT2)) THEN
      DO JIJ=1, KLEN
        PBUF3(JIJ) = ( PLT2(KBUF1(JIJ)+1,KBUF2(JIJ)+1)* PBUF2(JIJ)         &
                     -PLT2(KBUF1(JIJ)+1,KBUF2(JIJ)  )*(PBUF2(JIJ) - 1.0)) *  PBUF1(JIJ) &
                   -( PLT2(KBUF1(JIJ)  ,KBUF2(JIJ)+1)* PBUF2(JIJ)         &
                     -PLT2(KBUF1(JIJ)  ,KBUF2(JIJ)  )*(PBUF2(JIJ) - 1.0)) * (PBUF1(JIJ) - 1.0)
      ENDDO
      POUT2(:)=0.
      DO JIJ=1, KLEN
        POUT2(KBUF3(JIJ))=PBUF3(JIJ)
      ENDDO
    ENDIF

    !Interpolation and unpack 3
    IF(PRESENT(PLT3)) THEN
      DO JIJ=1, KLEN
        PBUF3(JIJ) = ( PLT3(KBUF2(JIJ)+1,KBUF1(JIJ)+1)* PBUF1(JIJ)         &
                     -PLT3(KBUF2(JIJ)+1,KBUF1(JIJ)  )*(PBUF1(JIJ) - 1.0)) *  PBUF2(JIJ) &
                   -( PLT3(KBUF2(JIJ)  ,KBUF1(JIJ)+1)* PBUF1(JIJ)         &
                     -PLT3(KBUF2(JIJ)  ,KBUF1(JIJ)  )*(PBUF1(JIJ) - 1.0)) * (PBUF2(JIJ) - 1.0)
      ENDDO
      POUT3(:)=0.
      DO JIJ=1, KLEN
        POUT3(KBUF3(JIJ))=PBUF3(JIJ)
      ENDDO
    ENDIF
  ENDIF

ELSE

  KLEN=0
  DO JIJ=D%NIJB, D%NIJE
    IF (LDMASK(JIJ)) THEN
      KLEN=KLEN+1

      !Indexes computation
      ZINDEX1 = MAX(1.00001, MIN(REAL(KNUM1)-0.00001, P11 * LOG(PIN1(JIJ)) + P12))
      IINDEX1 = INT(ZINDEX1)
      ZINDEX1 = ZINDEX1 - REAL(IINDEX1)
  
      ZINDEX2 = MAX(1.00001, MIN(REAL(KNUM1)-0.00001, P21 * LOG(PIN2(JIJ)) + P22))
      IINDEX2 = INT(ZINDEX2)
      ZINDEX2 = ZINDEX2 - REAL(IINDEX2)
  
      !Interpolations
      POUT1(JIJ) = ( PLT1(IINDEX1+1,IINDEX2+1)* ZINDEX2         &
                   -PLT1(IINDEX1+1,IINDEX2  )*(ZINDEX2 - 1.0)) *  ZINDEX1 &
                 -( PLT1(IINDEX1  ,IINDEX2+1)* ZINDEX2         &
                   -PLT1(IINDEX1  ,IINDEX2  )*(ZINDEX2 - 1.0)) * (ZINDEX1 - 1.0)

      IF(PRESENT(PLT2)) THEN
        POUT2(JIJ) = ( PLT2(IINDEX1+1,IINDEX2+1)* ZINDEX2         &
                     -PLT2(IINDEX1+1,IINDEX2  )*(ZINDEX2 - 1.0)) *  ZINDEX1 &
                   -( PLT2(IINDEX1  ,IINDEX2+1)* ZINDEX2         &
                     -PLT2(IINDEX1  ,IINDEX2  )*(ZINDEX2 - 1.0)) * (ZINDEX1 - 1.0)
      ENDIF

      IF(PRESENT(PLT3)) THEN
        POUT3(JIJ) = ( PLT3(IINDEX2+1,IINDEX1+1)* ZINDEX1         &
                     -PLT3(IINDEX2+1,IINDEX1  )*(ZINDEX1 - 1.0)) *  ZINDEX2 &
                   -( PLT3(IINDEX2  ,IINDEX1+1)* ZINDEX1         &
                     -PLT3(IINDEX2  ,IINDEX1  )*(ZINDEX1 - 1.0)) * (ZINDEX2 - 1.0)
      ENDIF

    ELSE
      POUT1(JIJ)=0.
      IF(PRESENT(PLT2)) POUT2(JIJ)=0.
      IF(PRESENT(PLT3)) POUT3(JIJ)=0.
    ENDIF
  ENDDO

ENDIF
END SUBROUTINE INTERP_MICRO_2D
