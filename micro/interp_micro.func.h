!These routines are intented to be included in the contains part of other subroutines.
!To allow the transformation for GPU, no local array must be declared.
!If a temporary local array is needed, it must be added as a buffer in the interface (IBUF?, ZBUF?)

SUBROUTINE INTERP_MICRO_1D(KPROMA, KSIZE, PIN, KNUM, P1, P2, &
                           LDPACK, LDMASK, KBUF1, KBUF2, PBUF1, PBUF2, &
                           KLEN, &
                           PLT1, POUT1, PLT2, POUT2, PLT3, POUT3)

IMPLICIT NONE

INTEGER,                    INTENT(IN)  :: KPROMA       !Array size
INTEGER,                    INTENT(IN)  :: KSIZE        !Last usefull array index
REAL,    DIMENSION(KPROMA), INTENT(IN)  :: PIN          !Input array
INTEGER,                    INTENT(IN)  :: KNUM         !Number of points in the look-up table
REAL,                       INTENT(IN)  :: P1           !Scaling factor
REAL,                       INTENT(IN)  :: P2           !Scaling factor
LOGICAL,                    INTENT(IN)  :: LDPACK       !.TRUE. to perform packing
LOGICAL, DIMENSION(KPROMA), INTENT(IN)  :: LDMASK       !Computation mask
INTEGER, DIMENSION(KPROMA), INTENT(OUT) :: KBUF1, KBUF2 !Buffer arrays
REAL,    DIMENSION(KPROMA), INTENT(OUT) :: PBUF1, PBUF2 !Buffer arrays
INTEGER,                    INTENT(OUT) :: KLEN         !Number of active points
REAL,    DIMENSION(KNUM),   INTENT(IN)            :: PLT1  !Look-up table
REAL,    DIMENSION(KPROMA), INTENT(OUT)           :: POUT1 !Interpolated values
REAL,    DIMENSION(KNUM),   INTENT(IN) , OPTIONAL :: PLT2
REAL,    DIMENSION(KPROMA), INTENT(OUT), OPTIONAL :: POUT2
REAL,    DIMENSION(KNUM),   INTENT(IN) , OPTIONAL :: PLT3
REAL,    DIMENSION(KPROMA), INTENT(OUT), OPTIONAL :: POUT3

INTEGER :: JL
INTEGER :: IINDEX
REAL :: ZINDEX

IF (LDPACK) THEN

  !Pack input array
  KLEN=0
  DO JL=1, KSIZE
    IF (LDMASK(JL)) THEN
      KLEN=KLEN+1
      PBUF1(KLEN)=PIN(JL)
      KBUF1(KLEN)=JL
    ENDIF
  ENDDO

  IF (KLEN>0) THEN
    !Index computation
    !$mnh_expand_array(JL=1:KLEN)
    PBUF1(1:KLEN) = MAX(1.00001, MIN(REAL(KNUM)-0.00001, P1 * LOG(PBUF1(1:KLEN)) + P2))
    KBUF2(1:KLEN) = INT(PBUF1(1:KLEN))
    PBUF1(1:KLEN) = PBUF1(1:KLEN) - REAL(KBUF2(1:KLEN))
    !$mnh_end_expand_array(JL=1:KLEN)

    !Interpolation and unpack
    !$mnh_expand_array(JL=1:KLEN)
    PBUF2(1:KLEN) = PLT1(KBUF2(1:KLEN)+1) *  PBUF1(1:KLEN)       &
                  &-PLT1(KBUF2(1:KLEN)  ) * (PBUF1(1:KLEN) - 1.0)
    !$mnh_end_expand_array(JL=1:KLEN)
    POUT1(:)=0.
    DO JL=1, KLEN
      POUT1(KBUF1(JL))=PBUF2(JL)
    ENDDO

    !Interpolation and unpack 2
    IF(PRESENT(PLT2)) THEN
      !$mnh_expand_array(JL=1:KLEN)
      PBUF2(1:KLEN) = PLT2(KBUF2(1:KLEN)+1) *  PBUF1(1:KLEN)       &
                    &-PLT2(KBUF2(1:KLEN)  ) * (PBUF1(1:KLEN) - 1.0)
      !$mnh_end_expand_array(JL=1:KLEN)
      POUT2(:)=0.
      DO JL=1, KLEN
        POUT2(KBUF1(JL))=PBUF2(JL)
      ENDDO
    ENDIF

    !Interpolation and unpack 3
    IF(PRESENT(PLT3)) THEN
      !$mnh_expand_array(JL=1:KLEN)
      PBUF2(1:KLEN) = PLT3(KBUF2(1:KLEN)+1) *  PBUF1(1:KLEN)       &
                    &-PLT3(KBUF2(1:KLEN)  ) * (PBUF1(1:KLEN) - 1.0)
      !$mnh_end_expand_array(JL=1:KLEN)
      POUT3(:)=0.
      DO JL=1, KLEN
        POUT3(KBUF1(JL))=PBUF2(JL)
      ENDDO
    ENDIF

  ENDIF

ELSE

  KLEN=0
  DO JL=1, KSIZE
    IF (LDMASK(JL)) THEN
      KLEN=KLEN+1

      !Index computation
      ZINDEX = MAX(1.00001, MIN(REAL(KNUM)-0.00001, P1 * LOG(PIN(JL)) + P2))
      IINDEX = INT(ZINDEX)
      ZINDEX = ZINDEX - REAL(IINDEX)

      !Interpolations
      POUT1(JL) = PLT1(IINDEX+1) *  ZINDEX       &
                &-PLT1(IINDEX  ) * (ZINDEX - 1.0)

      IF(PRESENT(PLT2)) THEN
        POUT2(JL) = PLT2(IINDEX+1) *  ZINDEX       &
                  &-PLT2(IINDEX  ) * (ZINDEX - 1.0)
      ENDIF

      IF(PRESENT(PLT3)) THEN
        POUT3(JL) = PLT3(IINDEX+1) *  ZINDEX       &
                  &-PLT3(IINDEX  ) * (ZINDEX - 1.0)
      ENDIF

    ELSE
      POUT1(JL) = 0.
      IF(PRESENT(PLT2)) POUT2(JL) = 0.
      IF(PRESENT(PLT3)) POUT3(JL) = 0.
    ENDIF
  ENDDO

ENDIF
END SUBROUTINE INTERP_MICRO_1D

SUBROUTINE INTERP_MICRO_2D(KPROMA, KSIZE, PIN1, PIN2, KNUM1, KNUM2, P11, P12, P21, P22,&
                           LDPACK, LDMASK, KBUF1, KBUF2, KBUF3, PBUF1, PBUF2, PBUF3, &
                           KLEN, &
                           PLT1, POUT1, PLT2, POUT2, PLT3, POUT3)

IMPLICIT NONE

INTEGER,                    INTENT(IN)  :: KPROMA       !Array size
INTEGER,                    INTENT(IN)  :: KSIZE        !Last usefull array index
REAL,    DIMENSION(KPROMA), INTENT(IN)  :: PIN1                !Input array
REAL,    DIMENSION(KPROMA), INTENT(IN)  :: PIN2                !Input array
INTEGER,                    INTENT(IN)  :: KNUM1               !First dimension of the look-up table
INTEGER,                    INTENT(IN)  :: KNUM2               !Second dimension of the look-up table
REAL,                       INTENT(IN)  :: P11                 !Scaling factor
REAL,                       INTENT(IN)  :: P12                 !Scaling factor
REAL,                       INTENT(IN)  :: P21                 !Scaling factor
REAL,                       INTENT(IN)  :: P22                 !Scaling factor
LOGICAL,                    INTENT(IN)  :: LDPACK              !.TRUE. to perform packing
LOGICAL, DIMENSION(KPROMA), INTENT(IN)  :: LDMASK              !Computation mask
INTEGER, DIMENSION(KPROMA), INTENT(OUT) :: KBUF1, KBUF2, KBUF3 !Buffer arrays
REAL,    DIMENSION(KPROMA), INTENT(OUT) :: PBUF1, PBUF2, PBUF3 !Buffer arrays
INTEGER,                    INTENT(OUT) :: KLEN                !Number of active points
REAL,    DIMENSION(KNUM1, KNUM2),   INTENT(IN)            :: PLT1  !Look-up table
REAL,    DIMENSION(KPROMA),         INTENT(OUT)           :: POUT1 !Interpolated values from the first look-up table
REAL,    DIMENSION(KNUM1, KNUM2),   INTENT(IN) , OPTIONAL :: PLT2  !Other look-up table
REAL,    DIMENSION(KPROMA),         INTENT(OUT), OPTIONAL :: POUT2 !Interpolated values from the second look-up table
REAL,    DIMENSION(KNUM2, KNUM1),   INTENT(IN) , OPTIONAL :: PLT3  !Another look-up table **CAUTION, TABLE IS REVERSED**
REAL,    DIMENSION(KPROMA),         INTENT(OUT), OPTIONAL :: POUT3 !Interpolated values from the third look-up table

INTEGER :: JL
INTEGER :: IINDEX1, IINDEX2
REAL :: ZINDEX1, ZINDEX2

IF (LDPACK) THEN

  !Pack input array
  KLEN=0
  DO JL=1, KSIZE
    IF (LDMASK(JL)) THEN
      KLEN=KLEN+1
      PBUF1(KLEN)=PIN1(JL)
      PBUF2(KLEN)=PIN2(JL)
      KBUF3(KLEN)=JL
    ENDIF
  ENDDO

  IF (KLEN>0) THEN
    !Index computation
    !$mnh_expand_array(JL=1:KLEN)
    PBUF1(1:KLEN) = MAX(1.00001, MIN(REAL(KNUM1)-0.00001, P11 * LOG(PBUF1(1:KLEN)) + P12))
    KBUF1(1:KLEN) = INT(PBUF1(1:KLEN))
    PBUF1(1:KLEN) = PBUF1(1:KLEN) - REAL(KBUF1(1:KLEN))

    PBUF2(1:KLEN) = MAX(1.00001, MIN(REAL(KNUM2)-0.00001, P21 * LOG(PBUF2(1:KLEN)) + P22))
    KBUF2(1:KLEN) = INT(PBUF2(1:KLEN))
    PBUF2(1:KLEN) = PBUF2(1:KLEN) - REAL(KBUF2(1:KLEN))
    !$mnh_end_expand_array(JL=1:KLEN)

    !Interpolation and unpack 1
    DO JL=1, KLEN
      PBUF3(JL) = ( PLT1(KBUF1(JL)+1,KBUF2(JL)+1)* PBUF2(JL)         &
                   -PLT1(KBUF1(JL)+1,KBUF2(JL)  )*(PBUF2(JL) - 1.0)) *  PBUF1(JL) &
                 -( PLT1(KBUF1(JL)  ,KBUF2(JL)+1)* PBUF2(JL)         &
                   -PLT1(KBUF1(JL)  ,KBUF2(JL)  )*(PBUF2(JL) - 1.0)) * (PBUF1(JL) - 1.0)
    ENDDO
    POUT1(:)=0.
    DO JL=1, KLEN
      POUT1(KBUF3(JL))=PBUF3(JL)
    ENDDO

    !Interpolation and unpack 2
    IF(PRESENT(PLT2)) THEN
      DO JL=1, KLEN
        PBUF3(JL) = ( PLT2(KBUF1(JL)+1,KBUF2(JL)+1)* PBUF2(JL)         &
                     -PLT2(KBUF1(JL)+1,KBUF2(JL)  )*(PBUF2(JL) - 1.0)) *  PBUF1(JL) &
                   -( PLT2(KBUF1(JL)  ,KBUF2(JL)+1)* PBUF2(JL)         &
                     -PLT2(KBUF1(JL)  ,KBUF2(JL)  )*(PBUF2(JL) - 1.0)) * (PBUF1(JL) - 1.0)
      ENDDO
      POUT2(:)=0.
      DO JL=1, KLEN
        POUT2(KBUF3(JL))=PBUF3(JL)
      ENDDO
    ENDIF

    !Interpolation and unpack 3
    IF(PRESENT(PLT3)) THEN
      DO JL=1, KLEN
        PBUF3(JL) = ( PLT3(KBUF2(JL)+1,KBUF1(JL)+1)* PBUF1(JL)         &
                     -PLT3(KBUF2(JL)+1,KBUF1(JL)  )*(PBUF1(JL) - 1.0)) *  PBUF2(JL) &
                   -( PLT3(KBUF2(JL)  ,KBUF1(JL)+1)* PBUF1(JL)         &
                     -PLT3(KBUF2(JL)  ,KBUF1(JL)  )*(PBUF1(JL) - 1.0)) * (PBUF2(JL) - 1.0)
      ENDDO
      POUT3(:)=0.
      DO JL=1, KLEN
        POUT3(KBUF3(JL))=PBUF3(JL)
      ENDDO
    ENDIF
  ENDIF

ELSE

  KLEN=0
  DO JL=1, KSIZE
    IF (LDMASK(JL)) THEN
      KLEN=KLEN+1

      !Indexes computation
      ZINDEX1 = MAX(1.00001, MIN(REAL(KNUM1)-0.00001, P11 * LOG(PIN1(JL)) + P12))
      IINDEX1 = INT(ZINDEX1)
      ZINDEX1 = ZINDEX1 - REAL(IINDEX1)
  
      ZINDEX2 = MAX(1.00001, MIN(REAL(KNUM1)-0.00001, P21 * LOG(PIN2(JL)) + P22))
      IINDEX2 = INT(ZINDEX2)
      ZINDEX2 = ZINDEX2 - REAL(IINDEX2)
  
      !Interpolations
      POUT1(JL) = ( PLT1(IINDEX1+1,IINDEX2+1)* ZINDEX2         &
                   -PLT1(IINDEX1+1,IINDEX2  )*(ZINDEX2 - 1.0)) *  ZINDEX1 &
                 -( PLT1(IINDEX1  ,IINDEX2+1)* ZINDEX2         &
                   -PLT1(IINDEX1  ,IINDEX2  )*(ZINDEX2 - 1.0)) * (ZINDEX1 - 1.0)

      IF(PRESENT(PLT2)) THEN
        POUT2(JL) = ( PLT2(IINDEX1+1,IINDEX2+1)* ZINDEX2         &
                     -PLT2(IINDEX1+1,IINDEX2  )*(ZINDEX2 - 1.0)) *  ZINDEX1 &
                   -( PLT2(IINDEX1  ,IINDEX2+1)* ZINDEX2         &
                     -PLT2(IINDEX1  ,IINDEX2  )*(ZINDEX2 - 1.0)) * (ZINDEX1 - 1.0)
      ENDIF

      IF(PRESENT(PLT3)) THEN
        POUT3(JL) = ( PLT3(IINDEX2+1,IINDEX1+1)* ZINDEX1         &
                     -PLT3(IINDEX2+1,IINDEX1  )*(ZINDEX1 - 1.0)) *  ZINDEX2 &
                   -( PLT3(IINDEX2  ,IINDEX1+1)* ZINDEX1         &
                     -PLT3(IINDEX2  ,IINDEX1  )*(ZINDEX1 - 1.0)) * (ZINDEX2 - 1.0)
      ENDIF

    ELSE
      POUT1(JL)=0.
      IF(PRESENT(PLT2)) POUT2(JL)=0.
      IF(PRESENT(PLT3)) POUT3(JL)=0.
    ENDIF
  ENDDO

ENDIF
END SUBROUTINE INTERP_MICRO_2D
