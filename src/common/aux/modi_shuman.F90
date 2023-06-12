!     ##################
      MODULE MODI_SHUMAN
!     ##################
!
IMPLICIT NONE
INTERFACE
!
FUNCTION DXF(PA)  RESULT(PDXF)
IMPLICIT NONE
REAL, DIMENSION(:,:,:), INTENT(IN)                :: PA     ! variable at flux
                                                            !  side
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PDXF   ! result at mass
                                                            ! localization 
END FUNCTION DXF
!
FUNCTION DXM(PA)  RESULT(PDXM)
IMPLICIT NONE
REAL, DIMENSION(:,:,:), INTENT(IN)                :: PA     ! variable at mass
                                                            ! localization
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PDXM   ! result at flux
                                                            ! side
END FUNCTION DXM
!
FUNCTION DYF(PA)  RESULT(PDYF)
IMPLICIT NONE
REAL, DIMENSION(:,:,:), INTENT(IN)                :: PA     ! variable at flux
                                                            !  side
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PDYF   ! result at mass
                                                            ! localization 
END FUNCTION DYF
!
FUNCTION DYM(PA)  RESULT(PDYM)
IMPLICIT NONE
REAL, DIMENSION(:,:,:), INTENT(IN)                :: PA     ! variable at mass
                                                            ! localization
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PDYM   ! result at flux
                                                            ! side
END FUNCTION DYM
!
FUNCTION DZF(PA,KKA,KKU,KL)  RESULT(PDZF)
IMPLICIT NONE
REAL, DIMENSION(:,:,:), INTENT(IN)                :: PA     ! variable at flux
                                                            !  side
INTEGER,              INTENT(IN),OPTIONAL         :: KKA, KKU ! near ground and uppest atmosphere array indexes
INTEGER,              INTENT(IN),OPTIONAL         :: KL     ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PDZF   ! result at mass localization 
END FUNCTION DZF
!
FUNCTION DZM(PA,KKA,KKU,KL)  RESULT(PDZM)
IMPLICIT NONE
REAL, DIMENSION(:,:,:), INTENT(IN)                :: PA     ! variable at mass
                                                            ! localization
INTEGER,              INTENT(IN),OPTIONAL         :: KKA, KKU ! near ground and uppest atmosphere array indexes
INTEGER,              INTENT(IN),OPTIONAL         :: KL     ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PDZM   ! result at flux side
END FUNCTION DZM
!
FUNCTION MXF(PA)  RESULT(PMXF)
IMPLICIT NONE
REAL, DIMENSION(:,:,:), INTENT(IN)                :: PA     ! variable at flux
                                                            !  side
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PMXF   ! result at mass
                                                            ! localization 
END FUNCTION MXF
!
FUNCTION MXM(PA)  RESULT(PMXM)
IMPLICIT NONE
REAL, DIMENSION(:,:,:), INTENT(IN)                :: PA     ! variable at mass localization
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PMXM   ! result at flux localization 
END FUNCTION MXM
!
FUNCTION MYF(PA)  RESULT(PMYF)
IMPLICIT NONE
REAL, DIMENSION(:,:,:), INTENT(IN)                :: PA     ! variable at flux
                                                            !   side
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PMYF   ! result at mass 
                                                            ! localization 
END FUNCTION MYF
!
FUNCTION MYM(PA)  RESULT(PMYM)
IMPLICIT NONE
REAL, DIMENSION(:,:,:), INTENT(IN)                :: PA     ! variable at mass localization
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PMYM   ! result at flux localization 
END  FUNCTION MYM
!
FUNCTION MZF(PA,KKA,KKU,KL)  RESULT(PMZF)
IMPLICIT NONE
REAL, DIMENSION(:,:,:), INTENT(IN)                :: PA     ! variable at flux side
INTEGER,              INTENT(IN),OPTIONAL         :: KKA, KKU ! near ground and uppest atmosphere array indexes
INTEGER,              INTENT(IN),OPTIONAL         :: KL     ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PMZF   ! result at mass
                                                            ! localization 
END FUNCTION MZF
!
FUNCTION MZM(PA,KKA,KKU,KL)  RESULT(PMZM)
IMPLICIT NONE
REAL, DIMENSION(:,:,:), INTENT(IN)                :: PA     ! variable at mass localization
INTEGER,              INTENT(IN),OPTIONAL         :: KKA, KKU ! near ground and uppest atmosphere array indexes
INTEGER,              INTENT(IN),OPTIONAL         :: KL     ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PMZM   ! result at flux localization 
END FUNCTION MZM
!
END INTERFACE
!
END MODULE MODI_SHUMAN
