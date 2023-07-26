MODULE MODE_BUDGET_PHY
USE MODD_BUDGET, ONLY : TBUDGETDATA
IMPLICIT NONE
CONTAINS

SUBROUTINE BUDGET_STORE_INIT(TPBUDGET, HSOURCE, PVARS)
  TYPE(TBUDGETDATA),      INTENT(INOUT) :: TPBUDGET ! Budget datastructure
  CHARACTER(LEN=*),       INTENT(IN)    :: HSOURCE  ! Name of the source term
  REAL, DIMENSION(:,:,:), INTENT(IN)    :: PVARS    ! Current value to be stored
END SUBROUTINE BUDGET_STORE_INIT
!
SUBROUTINE BUDGET_STORE_INIT_PHY(D,TPBUDGET, HSOURCE, PVARS)
  USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
  TYPE(DIMPHYEX_t),       INTENT(IN)    :: D
  TYPE(TBUDGETDATA),      INTENT(INOUT) :: TPBUDGET ! Budget datastructure
  CHARACTER(LEN=*),       INTENT(IN)    :: HSOURCE  ! Name of the source term
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)    :: PVARS    ! Current value to be stored
  CALL BUDGET_STORE_INIT(TPBUDGET, HSOURCE, PVARS)
END SUBROUTINE BUDGET_STORE_INIT_PHY
!
SUBROUTINE BUDGET_STORE_END(TPBUDGET, HSOURCE, PVARS)
  TYPE(TBUDGETDATA),      INTENT(INOUT) :: TPBUDGET ! Budget datastructure
  CHARACTER(LEN=*),       INTENT(IN)    :: HSOURCE  ! Name of the source term
  REAL, DIMENSION(:,:,:), INTENT(IN)    :: PVARS    ! Current value to be stored
  CALL BUDGET_DDH(PVARS, TPBUDGET%NBUDGET, HSOURCE, TPBUDGET%YDDDH, TPBUDGET%YDLDDH, TPBUDGET%YDMDDH)
END SUBROUTINE BUDGET_STORE_END
!
SUBROUTINE BUDGET_STORE_END_PHY(D,TPBUDGET, HSOURCE, PVARS)
  USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
  TYPE(DIMPHYEX_t),       INTENT(IN)    :: D
  TYPE(TBUDGETDATA),      INTENT(INOUT) :: TPBUDGET ! Budget datastructure
  CHARACTER(LEN=*),       INTENT(IN)    :: HSOURCE  ! Name of the source term
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)    :: PVARS    ! Current value to be stored
  CALL BUDGET_STORE_END(TPBUDGET, HSOURCE, PVARS)
END SUBROUTINE BUDGET_STORE_END_PHY
!
SUBROUTINE BUDGET_STORE_ADD_PHY(D,TPBUDGET, HSOURCE, PVARS)
  USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
  TYPE(DIMPHYEX_t),       INTENT(IN)    :: D
  TYPE(TBUDGETDATA),      INTENT(INOUT) :: TPBUDGET ! Budget datastructure
  CHARACTER(LEN=*),       INTENT(IN)    :: HSOURCE  ! Name of the source term
  REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)    :: PVARS    ! Current value to be stored
  CALL BUDGET_STORE_ADD(TPBUDGET, HSOURCE, PVARS)
END SUBROUTINE BUDGET_STORE_ADD_PHY
!
SUBROUTINE BUDGET_STORE_ADD(TPBUDGET, HSOURCE, PVARS)
  TYPE(TBUDGETDATA),      INTENT(INOUT) :: TPBUDGET ! Budget datastructure
  CHARACTER(LEN=*),       INTENT(IN)    :: HSOURCE  ! Name of the source term
  REAL, DIMENSION(:,:,:), INTENT(IN)    :: PVARS    ! Current value to be stored
  CALL BUDGET_DDH(PVARS, TPBUDGET%NBUDGET, HSOURCE, TPBUDGET%YDDDH, TPBUDGET%YDLDDH, TPBUDGET%YDMDDH, &
                 &LDISDIFF=.TRUE.)
END SUBROUTINE BUDGET_STORE_ADD

      SUBROUTINE BUDGET_DDH(PVARS,KBUDN,HBUVAR,YDDDH, YDLDDH, YDMDDH, LDISDIFF)
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
!     #####################################
!
!!   BUDGET - routine to call the BUDGET routine for AROME.
!!
!!   BEWARE THIS ROUTINE iS DIFFERENT FROM THE MNH ROUTINE BUDGET
!!   EVEN IF IT WEARS THE SAME NAME !!!
!!
!!    PURPOSE
!!    -------
!        This routine is an interface for the add_ddh subroutine.
!        It converts the selected field into klev reversed vertical
!        levels and attributes to the selected field are created.
!
!!**  METHOD
!!    ------
!!
!!      1st step: substract previous step (sequential approach in MNH)
!!      2nd step: reverse levels
!!      3rd step: multiplication by conversion factor for r-> q
!!                or Theta-> T
!!
!!      4rd step: call to add_ddh now that the field has an Aladin shape
!!
!!
!!    EXTERNAL
!!    --------
!!      ADD_FIELD_3D
!!      INVERT_VLEV
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_INTBUDGET:
!!
!!      PVARS_M(nlon,1,nlev,13) !13 different budgets
!!      VARMULT(:,:,13) ! allows to convert variables
!!
!!    REFERENCE
!!    ---------
!!      "New data flow for diagnostics in Arome/Arpege"
!!
!!    AUTHOR
!!    ------
!!      O.Riviere   17/07/08     * Meteo France *
!!      
!!
!!    MODIFICATIONS
!!    -------------
!!      F.Voitus    16/05/17  : Introduction of new DDH superstructure for budget
!!      S.Riette    Jan 2022  : LDISDIFF case
!!
!-------------------------------------------------------------------------------

!

USE MODDB_INTBUDGET,ONLY:TAB_VARMULT,TVARSM,CVARNAME,NLON
USE DDH_MIX, ONLY:NFLEVGDDH,NPROMADDH,ADD_FIELD_3D,  &
                 & TYP_DDH, NEW_ADD_FIELD_3D ! reference is Arpege
USE OML_MOD, ONLY : OML_MY_THREAD
USE YOMLDDH, ONLY : TLDDH
USE YOMMDDH, ONLY : TMDDH


IMPLICIT NONE
!
!
!*       0.1   Declarations of arguments :
!
REAL, DIMENSION(:,:,:), INTENT(IN) :: PVARS    ! source of the variable
INTEGER               , INTENT(IN) :: KBUDN    ! variable number

CHARACTER (LEN=*)     , INTENT(IN) :: HBUVAR   ! Identifier of the Budget
TYPE(TYP_DDH)         , INTENT(INOUT) :: YDDDH
TYPE(TLDDH)           , INTENT(IN)    :: YDLDDH
TYPE(TMDDH)           , INTENT(IN)    :: YDMDDH

LOGICAL, OPTIONAL     , INTENT(IN)    :: LDISDIFF ! PVARS contains the increment (default is .FALSE.)

!*       0.2   Declaration of local variables:
REAL,DIMENSION(NPROMADDH,NFLEVGDDH):: ZVARS
LOGICAL :: LINST,LDDH
INTEGER::IINCR,JLON,JLEV,IFDIA,IOFF
CHARACTER (LEN=4) :: CLPROC
CHARACTER (LEN=11) :: CLDDH
LOGICAL :: LISDIFF

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('BUDGET_DDH',0,ZHOOK_HANDLE)

IF (PRESENT(LDISDIFF)) THEN
  LISDIFF=LDISDIFF
ELSE
  LISDIFF=.FALSE.
ENDIF

IFDIA=SIZE(PVARS,1)
ZVARS(:,:)=0.
IF (SIZE(PVARS,3)==NFLEVGDDH+2) THEN
  IOFF=1
ELSE
  IOFF=0
ENDIF
!if length is less than 4, fill with budget old names
IF(LEN(HBUVAR)==1) THEN
  CLPROC=HBUVAR(1:MIN(4, LEN(HBUVAR)))//'_BU'
ELSE IF(LEN(HBUVAR)==2) THEN
  CLPROC=HBUVAR(1:MIN(4, LEN(HBUVAR)))//'_B'
ELSE IF(LEN(HBUVAR)==3) THEN
  CLPROC=HBUVAR(1:MIN(4, LEN(HBUVAR)))//'_'
ELSE
  CLPROC=HBUVAR(1:MIN(4, LEN(HBUVAR)))
END IF
!
IF (YDLDDH%LDDH_OMP) THEN
  CLDDH='T'//YDDDH%YVARMULT(KBUDN)%CNAME//CLPROC
ELSE
  CLDDH='T'//CVARNAME(KBUDN)//CLPROC
ENDIF

! depi not stored through call to budget but add_field
IF ((CLPROC=='DEPI').OR.(CLPROC=='CEDS')) THEN
  IF (LHOOK) CALL DR_HOOK('BUDGET_DDH',1,ZHOOK_HANDLE)
  RETURN
ENDIF

!1. Substraction of value at previous process and updates PVARSM

IF (YDLDDH%LDDH_OMP) THEN
  IF (CLPROC=='INIF') THEN
    DO JLEV=1,NFLEVGDDH
      DO JLON=1,IFDIA
        YDDDH%RVARSM(JLON,1,JLEV,KBUDN)=PVARS(JLON,1,JLEV+IOFF)
        ZVARS(JLON,JLEV)=PVARS(JLON,1,JLEV+IOFF)
      ENDDO
    ENDDO
  ELSE
    IF (LISDIFF) THEN
      DO JLEV=1,NFLEVGDDH
        DO JLON=1,IFDIA
          ZVARS(JLON,JLEV)=PVARS(JLON,1,JLEV+IOFF)
          YDDDH%RVARSM(JLON,1,JLEV,KBUDN)=YDDDH%RVARSM(JLON,1,JLEV,KBUDN)+PVARS(JLON,1,JLEV+IOFF)
        ENDDO
      ENDDO
    ELSE
      DO JLEV=1,NFLEVGDDH
        DO JLON=1,IFDIA
          ZVARS(JLON,JLEV)=PVARS(JLON,1,JLEV+IOFF)-YDDDH%RVARSM(JLON,1,JLEV,KBUDN)
          YDDDH%RVARSM(JLON,1,JLEV,KBUDN)=PVARS(JLON,1,JLEV+IOFF)
        ENDDO
      ENDDO
    ENDIF
  ENDIF
ELSE
  IF (CLPROC=='INIF') THEN
    DO JLEV=1,NFLEVGDDH
       DO JLON=1,IFDIA
         TVARSM(JLON,1,JLEV,KBUDN)=PVARS(JLON,1,JLEV+IOFF)
         ZVARS(JLON,JLEV)=PVARS(JLON,1,JLEV+IOFF)
       ENDDO
    ENDDO
  ELSE
    IF (LISDIFF) THEN
      DO JLEV=1,NFLEVGDDH
        DO JLON=1,IFDIA
          ZVARS(JLON,JLEV)=PVARS(JLON,1,JLEV+IOFF)
          TVARSM(JLON,1,JLEV,KBUDN)=TVARSM(JLON,1,JLEV,KBUDN)+PVARS(JLON,1,JLEV+IOFF)
        ENDDO
      ENDDO
    ELSE
      DO JLEV=1,NFLEVGDDH
        DO JLON=1,IFDIA
          ZVARS(JLON,JLEV)=PVARS(JLON,1,JLEV+IOFF)-TVARSM(JLON,1,JLEV,KBUDN)
          TVARSM(JLON,1,JLEV,KBUDN)=PVARS(JLON,1,JLEV+IOFF)
        ENDDO
      ENDDO
    ENDIF
  ENDIF
ENDIF



!2. Reverse levels MNH-> ALD
!IINCR=-1
!CALL INVERT_VLEV(1,NLON,NFLEVGDDH,IINCR,ZVARS,PVARS2)

!3. CONVERSION
! converting to desired budget variables

IF (YDLDDH%LDDH_OMP) THEN
  ZVARS(:,:)=ZVARS(:,:)*YDDDH%YVARMULT(KBUDN)%RVAL(:,:)
ELSE
  ZVARS(:,:)=ZVARS(:,:)*TAB_VARMULT(KBUDN)%VARMULT(:,:)
ENDIF


!4. CALL TO ADD_FIELD


LDDH=.TRUE.
LINST=.TRUE.
! saves ZVARS with NAME HBUVAR as a Tendency from AROME
! and it is an INSTantaneous field
IF (CLPROC/='INIF') THEN
  IF (YDLDDH%LDDH_OMP) THEN
    CALL NEW_ADD_FIELD_3D(YDMDDH,ZVARS,CLDDH,YDDDH)
  ELSE
    CALL ADD_FIELD_3D(YDLDDH,ZVARS,CLDDH,'T','AROME',LINST,LDDH)
  ENDIF
ENDIF

IF (LHOOK) CALL DR_HOOK('BUDGET_DDH',1,ZHOOK_HANDLE)
END SUBROUTINE BUDGET_DDH
END MODULE MODE_BUDGET_PHY
