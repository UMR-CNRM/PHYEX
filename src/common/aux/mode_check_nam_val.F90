MODULE MODE_CHECK_NAM_VAL
!> @file
!!      *MODE_CHECK_NAM_VAL" - Module containing the routines to control the different kind of variables
!!                             read from namelist
!!
USE MODE_MSG, ONLY: PRINT_MSG, NVERB_FATAL
IMPLICIT NONE
CONTAINS
SUBROUTINE CHECK_NAM_VAL_CHAR(KLUOUT, HNAME, HVAR, HVALUE1, HVALUE2, HVALUE3, HVALUE4, HVALUE5, &
                             &HVALUE6, HVALUE7, HVALUE8, HVALUE9, HVALUE10, HVALUE11, HVALUE12)
!!
!!      *CHECK_NAM_VAL* - Control of CHARACTER variables
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to control the validity of CHARACTER variables
!!
!!
!!    AUTHOR
!!    ------
!!     S. Riette
!!
!!    MODIFICATIONS
!!    -------------
!!
!!    -  Original      Feb 2023, from Méso-NH code
!!
!-------------------------------------------------------------------------------
!
!**       DECLARATIONS
!
IMPLICIT NONE
INTEGER,          INTENT(IN)           :: KLUOUT   !< output listing logical unit
CHARACTER(LEN=*), INTENT(IN)           :: HNAME    !< name of the variable to test
CHARACTER(LEN=*), INTENT(IN)           :: HVAR     !< variable to test
CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: HVALUE1  !< Authorised value 1
CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: HVALUE2  !< Authorised value 2
CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: HVALUE3  !< Authorised value 3
CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: HVALUE4  !< Authorised value 4
CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: HVALUE5  !< Authorised value 5
CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: HVALUE6  !< Authorised value 6
CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: HVALUE7  !< Authorised value 7
CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: HVALUE8  !< Authorised value 8
CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: HVALUE9  !< Authorised value 9
CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: HVALUE10 !< Authorised value 10
CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: HVALUE11 !< Authorised value 11
CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: HVALUE12 !< Authorised value 12
!
!**       CONTROLS
!
IF ( PRESENT (HVALUE1) ) THEN
  IF ( HVAR==HVALUE1 ) RETURN
END IF
!
IF ( PRESENT (HVALUE2) ) THEN
  IF ( HVAR==HVALUE2 ) RETURN
END IF
!
IF ( PRESENT (HVALUE3) ) THEN
  IF ( HVAR==HVALUE3 ) RETURN
END IF
!
IF ( PRESENT (HVALUE4) ) THEN
  IF ( HVAR==HVALUE4 ) RETURN
END IF
!
IF ( PRESENT (HVALUE5) ) THEN
  IF ( HVAR==HVALUE5 ) RETURN
END IF
!
IF ( PRESENT (HVALUE6) ) THEN
  IF ( HVAR==HVALUE6 ) RETURN
END IF
!
IF ( PRESENT (HVALUE7) ) THEN
  IF ( HVAR==HVALUE7 ) RETURN
END IF
!
IF ( PRESENT (HVALUE8) ) THEN
  IF ( HVAR==HVALUE8 ) RETURN
END IF
!
IF ( PRESENT (HVALUE9) ) THEN
  IF ( HVAR==HVALUE9 ) RETURN
END IF
!
IF ( PRESENT (HVALUE10) ) THEN
  IF ( HVAR==HVALUE10 ) RETURN
END IF
!
IF ( PRESENT (HVALUE11) ) THEN
  IF ( HVAR==HVALUE11 ) RETURN
END IF
!
IF ( PRESENT (HVALUE12) ) THEN
  IF ( HVAR==HVALUE12 ) RETURN
END IF
!
!** PRINTS AND ABORT
!
WRITE (KLUOUT,*) ' '
WRITE (KLUOUT,*) 'FATAL ERROR:'
WRITE (KLUOUT,*) '-----------'
WRITE (KLUOUT,*) ' '
WRITE (KLUOUT,*) 'Value "', HVAR, '" is not allowed for variable ', HNAME
WRITE (KLUOUT,*) ' '
WRITE (KLUOUT,*) 'Possible values are:'
IF ( PRESENT (HVALUE1) ) WRITE (KLUOUT,*) '"',HVALUE1,'"'
IF ( PRESENT (HVALUE2) ) WRITE (KLUOUT,*) '"',HVALUE2,'"'
IF ( PRESENT (HVALUE3) ) WRITE (KLUOUT,*) '"',HVALUE3,'"'
IF ( PRESENT (HVALUE4) ) WRITE (KLUOUT,*) '"',HVALUE4,'"'
IF ( PRESENT (HVALUE5) ) WRITE (KLUOUT,*) '"',HVALUE5,'"'
IF ( PRESENT (HVALUE6) ) WRITE (KLUOUT,*) '"',HVALUE6,'"'
IF ( PRESENT (HVALUE7) ) WRITE (KLUOUT,*) '"',HVALUE7,'"'
IF ( PRESENT (HVALUE8) ) WRITE (KLUOUT,*) '"',HVALUE8,'"'
IF ( PRESENT (HVALUE9) ) WRITE (KLUOUT,*) '"',HVALUE9,'"'
IF ( PRESENT (HVALUE10) ) WRITE (KLUOUT,*) '"',HVALUE10,'"'
IF ( PRESENT (HVALUE11) ) WRITE (KLUOUT,*) '"',HVALUE11,'"'
IF ( PRESENT (HVALUE12) ) WRITE (KLUOUT,*) '"',HVALUE12,'"'
FLUSH(UNIT=KLUOUT)
!
CALL PRINT_MSG(NVERB_FATAL, 'GEN', 'CHECK_NAM_VAL_CHAR', TRIM(HVAR) // ' is not allowed for variable ' // TRIM(HNAME))
!
END SUBROUTINE CHECK_NAM_VAL_CHAR

SUBROUTINE CHECK_NAM_VAL_REAL(KLUOUT, HNAME, PVALUE, CDSIGN1, PVAL1, CDSIGN2, PVAL2)
!!
!!      *CHECK_NAM_VAL* - Control of CHARACTER variables
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to control the validity of REAL variables
!!
!!
!!    AUTHOR
!!    ------
!!     S. Riette
!!
!!    MODIFICATIONS
!!    -------------
!!
!!    -  Original      Feb 2023
!!
!-------------------------------------------------------------------------------
!
!**       DECLARATIONS
!
IMPLICIT NONE
INTEGER,          INTENT(IN)           :: KLUOUT   !< output listing logical unit
CHARACTER(LEN=*), INTENT(IN)           :: HNAME    !< name of the variable to test
REAL,             INTENT(IN)           :: PVALUE   !< variable to test
CHARACTER(LEN=*), INTENT(IN)           :: CDSIGN1  !< sign for the first verification
REAL,             INTENT(IN)           :: PVAL1    !< bound for the first verification
CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: CDSIGN2  !< sign for the second verification
REAL,             INTENT(IN), OPTIONAL :: PVAL2    !< bound for the second verification

INTEGER :: II, INUM
REAL :: ZVAL
CHARACTER(LEN=2) :: CSIGN
LOGICAL :: LOK
CHARACTER(LEN=10) :: CHAR_VAL
!
!**       CONTROLS
!
LOK=.TRUE.
INUM=1
IF(PRESENT(CDSIGN2)) INUM=2
DO II=1, INUM
  IF(II==1) THEN
    ZVAL=PVAL1
    CSIGN=CDSIGN1(1:MIN(2, LEN(CDSIGN1)))
  ELSE
    ZVAL=PVAL2
    CSIGN=CDSIGN2(1:MIN(2, LEN(CDSIGN2)))
  ENDIF
  SELECT CASE (CSIGN)
    CASE ('<') 
      LOK=LOK .AND. PVALUE < ZVAL
    CASE ('<=')
      LOK=LOK .AND. PVALUE <= ZVAL
    CASE ('>') 
      LOK=LOK .AND. PVALUE > ZVAL
    CASE ('>=')
      LOK=LOK .AND. PVALUE >= ZVAL
    CASE ('==')
      LOK=LOK .AND. PVALUE == ZVAL
    CASE DEFAULT
       CALL PRINT_MSG(NVERB_FATAL, 'GEN', 'CHECK_NAM_VAL_REAL', TRIM(CSIGN) // ' is not allowed as comparator')
  END SELECT
ENDDO
!
!** PRINTS AND ABORT
!
IF(.NOT. LOK) THEN
  WRITE(KLUOUT,*) ' '
  WRITE(KLUOUT,*) 'FATAL ERROR:'
  WRITE(KLUOUT,*) '-----------'
  WRITE(KLUOUT,*) ' '
  WRITE(KLUOUT,*) 'Value "', PVALUE, '" is not allowed for variable ', HNAME
  WRITE(KLUOUT,*) ' '
  WRITE(KLUOUT,*) 'Possible values are such as:'
  WRITE(KLUOUT,*) CDSIGN1, PVAL1
  IF (PRESENT(CDSIGN2)) WRITE(KLUOUT, *) CDSIGN2, PVAL2
  FLUSH(UNIT=KLUOUT)
  !
  WRITE(UNIT=CHAR_VAL, FMT=*) PVALUE
  CALL PRINT_MSG(NVERB_FATAL, 'GEN', 'CHECK_NAM_VAL_REAL', TRIM(CHAR_VAL) // ' is not allowed for variable ' // TRIM(HNAME))
ENDIF
!
END SUBROUTINE CHECK_NAM_VAL_REAL

SUBROUTINE CHECK_NAM_VAL_INT(KLUOUT, HNAME, KVALUE, CDSIGN1, KVAL1, CDSIGN2, KVAL2)
!!
!!      *CHECK_NAM_VAL* - Control of CHARACTER variables
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to control the validity of REAL variables
!!
!!
!!    AUTHOR
!!    ------
!!     S. Riette
!!
!!    MODIFICATIONS
!!    -------------
!!
!!    -  Original      Feb 2023
!!
!-------------------------------------------------------------------------------
!
!**       DECLARATIONS
!
IMPLICIT NONE
INTEGER,          INTENT(IN)           :: KLUOUT   !< output listing logical unit
CHARACTER(LEN=*), INTENT(IN)           :: HNAME    !< name of the variable to test
INTEGER,          INTENT(IN)           :: KVALUE   !< variable to test
CHARACTER(LEN=*), INTENT(IN)           :: CDSIGN1  !< sign for the first verification
INTEGER,          INTENT(IN)           :: KVAL1    !< bound for the first verification
CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: CDSIGN2  !< sign for the second verification
INTEGER,          INTENT(IN), OPTIONAL :: KVAL2    !< bound for the second verification

INTEGER :: II, INUM
INTEGER :: IVAL
CHARACTER(LEN=2) :: CSIGN
LOGICAL :: LOK
CHARACTER(LEN=10) :: CHAR_VAL
!
!**       CONTROLS
!
LOK=.TRUE.
INUM=1
IF(PRESENT(CDSIGN2)) INUM=2
DO II=1, INUM
  IF(II==1) THEN
    IVAL=KVAL1
    CSIGN=CDSIGN1(1:MIN(2, LEN(CDSIGN1)))
  ELSE
    IVAL=KVAL2
    CSIGN=CDSIGN2(1:MIN(2, LEN(CDSIGN2)))
  ENDIF
  SELECT CASE (CSIGN)
    CASE ('<') 
      LOK=LOK .AND. KVALUE < IVAL
    CASE ('<=')
      LOK=LOK .AND. KVALUE <= IVAL
    CASE ('>') 
      LOK=LOK .AND. KVALUE > IVAL
    CASE ('>=')
      LOK=LOK .AND. KVALUE >= IVAL
    CASE ('==')
      LOK=LOK .AND. KVALUE == IVAL
    CASE DEFAULT
       CALL PRINT_MSG(NVERB_FATAL, 'GEN', 'CHECK_NAM_VAL_REAL', TRIM(CSIGN) // ' is not allowed as comparator')
  END SELECT
ENDDO
!
!** PRINTS AND ABORT
!
IF(.NOT. LOK) THEN
  WRITE(KLUOUT,*) ' '
  WRITE(KLUOUT,*) 'FATAL ERROR:'
  WRITE(KLUOUT,*) '-----------'
  WRITE(KLUOUT,*) ' '
  WRITE(KLUOUT,*) 'Value "', KVALUE, '" is not allowed for variable ', HNAME
  WRITE(KLUOUT,*) ' '
  WRITE(KLUOUT,*) 'Possible values are such as:'
  WRITE(KLUOUT,*) CDSIGN1, KVAL1
  IF (PRESENT(CDSIGN2)) WRITE(KLUOUT, *) CDSIGN2, KVAL2
  FLUSH(UNIT=KLUOUT)
  !
  WRITE(UNIT=CHAR_VAL, FMT=*) KVALUE
  CALL PRINT_MSG(NVERB_FATAL, 'GEN', 'CHECK_NAM_VAL_REAL', TRIM(CHAR_VAL) // ' is not allowed for variable ' // TRIM(HNAME))
ENDIF
!
END SUBROUTINE CHECK_NAM_VAL_INT
!
END MODULE MODE_CHECK_NAM_VAL
