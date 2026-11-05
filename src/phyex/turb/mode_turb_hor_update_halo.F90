MODULE MODE_TURB_HOR_UPDATE_HALO
IMPLICIT NONE
CONTAINS
SUBROUTINE TURB_HOR_UPDATE_HALO(D, HLBCX, HLBCY, ODOUPDATE, TZFIELDS_ll, KRR, KSV, &
                                PUM, PVM, PWM, PTHLM, PTKEM, PRM, PSVM)
!       
!*      0. DECLARATIONS
!          ------------
!
USE MODE_ll, ONLY: LEAST_LL, LNORTH_LL, LSOUTH_LL, LWEST_LL, UPDATE_HALO_LL
USE MODD_ARGSLIST_ll, ONLY : LIST_ll
USE MODD_DIMPHYEX, ONLY : DIMPHYEX_t
!
IMPLICIT NONE
!
!
!*       0.1  declaration of arguments
!
!
TYPE(DIMPHYEX_t),       INTENT(IN)   :: D
CHARACTER (LEN=*), DIMENSION(:), INTENT(IN)       ::  HLBCX,HLBCY
LOGICAL, INTENT(IN) :: ODOUPDATE
TYPE(LIST_ll), POINTER, INTENT(INOUT) :: TZFIELDS_ll
INTEGER, INTENT(IN) :: KRR, KSV
REAL, DIMENSION(D%NIT, D%NJT, D%NKT), INTENT(INOUT) :: PUM, PVM, PWM, PTHLM, PTKEM
REAL, DIMENSION(D%NIT, D%NJT, D%NKT, KRR), INTENT(INOUT) :: PRM
REAL, DIMENSION(D%NIT, D%NJT, D%NKT, KSV), INTENT(INOUT) :: PSVM
!
!
!
!*       0.2  declaration of local variables
!
INTEGER :: IINFO_ll, IIB, IIE, IJB, IJE, IKB, IKE
!
! ---------------------------------------------------------------------------
!
#ifndef MNH_OPENACC
    IF (ODOUPDATE) CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
#else
    !
    ! Pour le update_halo_ll, on fait des update self, puis le update_halo_ll sur CPU
    ! on copie ensuite les champs avec le halo a jour sur le device
    !
    IF (ODOUPDATE) THEN
      !$update self(PUM,PVM,PWM,PTHLM,PTKEM)
      IF (KSV>0) THEN
        !$update self(PSVM(:,:,:,1:KSV))
      END IF
      IF (KRR>0) THEN
        !$update self(PRM(:,:,:,1:KRR))
      END IF
      CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
      IF (KRR>0) THEN
        !$update device(PRM(:,:,:,1:KRR))
      END IF
      IF (KSV>0) THEN
        !$update device(PSVM(:,:,:,1:KSV))
      END IF
      !$update device(PUM,PVM,PWM,PTHLM,PTKEM)
    ENDIF
#endif
    !
    IF ( HLBCX(1) /= "CYCL" .AND. LWEST_ll()) THEN
!$acc kernels present_cr(PSVM,PRM,PUM,PVM,PWM,PTHLM,PTKEM)
       PUM(IIB  ,:,:)=PUM(IIB  ,:,:)
       PVM(IIB-1,:,:)=PVM(IIB-1,:,:)
       PWM(IIB-1,:,:)=PWM(IIB-1,:,:)
       PTHLM(IIB-1,:,:)=PTHLM(IIB-1,:,:)
       PTKEM(IIB-1,:,:)=PTKEM(IIB-1,:,:)
       IF (KSV>0) PSVM(IIB-1,:,:,:)=PSVM(IIB-1,:,:,:)
       IF (KRR>0) PRM (IIB-1,:,:,:)=PRM (IIB-1,:,:,:)
!$acc end kernels
     ENDIF
     !
     IF ( HLBCX(2) /= "CYCL" .AND. LEAST_ll()) THEN
!$acc kernels present_cr(PSVM,PRM,PUM,PVM,PWM,PTHLM,PTKEM)
       PUM(IIE+1,:,:)=PUM(IIE+1,:,:)
       PVM(IIE+1,:,:)=PVM(IIE+1,:,:)
       PWM(IIE+1,:,:)=PWM(IIE+1,:,:)
       PTHLM(IIE+1,:,:)=PTHLM(IIE+1,:,:)
       PTKEM(IIE+1,:,:)=PTKEM(IIE+1,:,:)
       IF (KSV>0) PSVM(IIE+1,:,:,:)=PSVM(IIE+1,:,:,:)
       IF (KRR>0) PRM (IIE+1,:,:,:)=PRM(IIE+1,:,:,:)
!$acc end kernels
     ENDIF
     !
     IF ( HLBCY(1) /= "CYCL" .AND. LSOUTH_ll()) THEN
!$acc kernels present_cr(PSVM,PRM,PUM,PVM,PWM,PTHLM,PTKEM)
       PUM(:,IJB-1,:)=PUM(:,IJB-1,:)
       PVM(:,IJB  ,:)=PVM(:,IJB  ,:)
       PWM(:,IJB-1,:)=PWM(:,IJB-1,:)
       PTHLM(:,IJB-1,:)=PTHLM(:,IJB-1,:)
       PTKEM(:,IJB-1,:)=PTKEM(:,IJB-1,:)
       IF (KSV>0) PSVM(:,IJB-1,:,:)=PSVM(:,IJB-1,:,:)
       IF (KRR>0) PRM (:,IJB-1,:,:)=PRM (:,IJB-1,:,:)
!$acc end kernels
     ENDIF
     !
     IF ( HLBCY(2) /= "CYCL" .AND. LNORTH_ll()) THEN
!$acc kernels present_cr(PSVM,PRM,PUM,PVM,PWM,PTHLM,PTKEM)
       PUM(:,IJE+1,:)=PUM(:,IJE+1,:)
       PVM(:,IJE+1,:)=PVM(:,IJE+1,:)
       PWM(:,IJE+1,:)=PWM(:,IJE+1,:)
       PTHLM(:,IJE+1,:)=PTHLM(:,IJE+1,:)
       PTKEM(:,IJE+1,:)=PTKEM(:,IJE+1,:)
       IF (KSV>0) PSVM(:,IJE+1,:,:)=PSVM(:,IJE+1,:,:)
       IF (KRR>0) PRM (:,IJE+1,:,:)=PRM(:,IJE+1,:,:)
!$acc end kernels
     ENDIF
     !
!$acc kernels present_cr(PSVM,PRM)
     PUM(:,:,IKB-1)=PUM(:,:,IKB)
     PVM(:,:,IKB-1)=PVM(:,:,IKB)
     PWM(:,:,IKB-1)=PWM(:,:,IKB)
     PTHLM(:,:,IKB-1)=PTHLM(:,:,IKB)
     PTKEM(:,:,IKB-1)=PTKEM(:,:,IKB)
     IF (KSV>0) PSVM(:,:,IKB-1,:)=PSVM(:,:,IKB,:)
     IF (KRR>0) PRM (:,:,IKB-1,:)=PRM (:,:,IKB,:)
     !
     PUM(:,:,IKE+1)=PUM(:,:,IKE)
     PVM(:,:,IKE+1)=PVM(:,:,IKE)
     PWM(:,:,IKE+1)=PWM(:,:,IKE)
     PTHLM(:,:,IKE+1)=PTHLM(:,:,IKE)
     PTKEM(:,:,IKE+1)=PTKEM(:,:,IKE)
     IF (KSV>0) PSVM(:,:,IKE+1,:)=PSVM(:,:,IKE,:)
     IF (KRR>0) PRM (:,:,IKE+1,:)=PRM (:,:,IKE,:)
!$acc end kernels
     !
END SUBROUTINE TURB_HOR_UPDATE_HALO
END MODULE MODE_TURB_HOR_UPDATE_HALO

