!     ######spl
     SUBROUTINE CH_AER_COAG(PM, PSIG0, PRG0, PN0,PDMINTRA,PDMINTER,PTGAS,PMU,&
                              PLAMBDA,PRHOP0)
     USE PARKIND1, ONLY : JPRB
     USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!!   #############################################
!!
!!   PURPOSE
!!   -------
!!
!!   compute the terms due to Brownian, turbulent and Gravitational
!!   coagulation:
!!   a set of arrays are used to evaluate the double integral
!!    REFERENCE
!!    ---------
!!    none
!!
!!    AUTHOR
!!    ------
!!    Vincent Crassier (LA)
!!
!!    MODIFICATIONS
!!    -------------
!*****************************************************************
! * Arrays of numerical evaluation of coagulation terms
!   in the free-molecule regime (computed from the ESMAP code)
!
! ZINTRA     - Intamodal coagulation, mode i,j 0th and 6th Moment
!
! ZINTER0I   - Intermodal coagulation, mode i, 0th Moment
! ZINTER3I   - Intermodal coagulation, mode i, 3rd Moment
! ZINTER6I   - Intermodal coagulation, mode i, 6th Moment
! ZINTER6J   - Intermodal coagulation, mode j, 6th Moment
!
! * Variables used during the coefficients evaluation
! ZXI(i)     - Variables values at the array nodes
! ZXINT(i)   - Variables values where the interpolation
!             is to be made
!
! intramodal coagulation
!
! ZXINTRAMIN     - Minimal value of ln(sigma)
! ZXINTRAMAX     - Maximal value of ln(sigma)
! ZDXINTRA       - Step of ln(sigma) in the array
!
! intermodal coagulation
!
! ZXINTERMIN(i)  - Minimal value of the variable i
! ZXINTERMAX(i)  - Maximal value of the variable i
! ZDXINTER(i)    - Step of the variable i in the arrays
!
! i=1           - ln(sigmaj)
! i=2           - ln(sigmai)
! i=3           - ln((ZR=Rgj/Rgi)**2)
!
!***************************************************************
!!
!!   EXTERNAL
!!   -------
!!
!!   IMPLICIT ARGUMENTS
!!   ------------------
!!
USE MODD_CH_AEROSOL
!!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PM,PRHOP0 
REAL, DIMENSION(:),   INTENT(INOUT) :: PLAMBDA, PMU
REAL, DIMENSION(:,:), INTENT(INOUT) :: PDMINTRA
REAL, DIMENSION(:,:), INTENT(INOUT) :: PDMINTER
REAL, DIMENSION(:),INTENT(IN) :: PTGAS
REAL,   DIMENSION(:,:), INTENT(IN) :: PSIG0, PRG0, PN0
!
!*       0.2   Declarations of local variables
!
INTEGER :: JI,JJ
!
REAL :: ZTURBDS ! Rate of dissipation of kinetic energy per unit mass (m2/s3)
!
REAL, DIMENSION(SIZE(PM,1)) :: ZKFM,ZKNC
!REAL, DIMENSION(SIZE(PM,1)) :: ZKTURB,ZKGRAV,ZR3,ZRM4
REAL, DIMENSION(SIZE(PM,1)) :: ZR,ZR2,ZR4
REAL, DIMENSION(SIZE(PM,1)) :: ZRM,ZRM2,ZRM3
REAL, DIMENSION(SIZE(PM,1)) :: ZKNG
REAL, DIMENSION(SIZE(PM,1)) :: ZAI,ZKNGI,ZAJ,ZKNGJ
REAL, DIMENSION(SIZE(PM,1)) :: ZINTRA0NC,ZINTRA0FM,ZINTRA0
REAL, DIMENSION(SIZE(PM,1)) :: ZINTRA3NC,ZINTRA3FM,ZINTRA3
REAL, DIMENSION(SIZE(PM,1)) :: ZINTRA6NC,ZINTRA6FM,ZINTRA6
REAL, DIMENSION(SIZE(PM,1)) :: ZINTERNC,ZINTERFM,ZINTER
REAL, DIMENSION(SIZE(PM,1)) :: ZAPPROX
!
REAL, DIMENSION(SIZE(PM,1)) :: ZA,ZB,ZC,ZD
REAL, DIMENSION(SIZE(PM,1)) :: ZRGJ, ZRGI, ZRG
!
REAL, DIMENSION(SIZE(PM,1)) :: ZERF0,ZPHI0,ZXi,ZSOL
REAL, DIMENSION(SIZE(PM,1)) :: ZERF3,ZPHI3
REAL, DIMENSION(SIZE(PM,1)) :: ZERF6,ZPHI6
REAL, DIMENSION(SIZE(PM,1),JPMODE) :: ZINVSIG,ZLNDG
!
REAL, DIMENSION(SIZE(PM,1),JPMODE) :: ZESG01,ZESG04,ZESG05,ZESG08,ZESG09
REAL, DIMENSION(SIZE(PM,1),JPMODE) :: ZESG12,ZESG16
REAL, DIMENSION(SIZE(PM,1),JPMODE) :: ZESG20,ZESG24,ZESG25,ZESG28
REAL, DIMENSION(SIZE(PM,1),JPMODE) :: ZESG32,ZESG36
REAL, DIMENSION(SIZE(PM,1),JPMODE) :: ZESG49
REAL, DIMENSION(SIZE(PM,1),JPMODE) :: ZESG52
REAL, DIMENSION(SIZE(PM,1),JPMODE) :: ZESG64
REAL, DIMENSION(SIZE(PM,1),JPMODE) :: ZESG81,ZESG85
REAL, DIMENSION(SIZE(PM,1),JPMODE) :: ZESG100,ZESG121,ZESG144,ZESG169,ZESG196
REAL, DIMENSION(SIZE(PM,1),JPMODE) :: ZESG256
REAL, DIMENSION(SIZE(PM,1)) :: ZRB0,ZRB6
REAL, DIMENSION(SIZE(PM,1)) :: ZRES
!-------------------------------------------------------------------------------!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('CH_AER_COAG',0,ZHOOK_HANDLE)
ZTURBDS=0.001
ZKNC(:)=2.*XBOLTZ*PTGAS(:)/(3.*PMU(:))
!
PDMINTRA(:,:)=0.
PDMINTER(:,:)=0.
!
!****************************************************************
! Initialisation des variables utilisees dans le calcul des
! coefficients de coagulation
!****************************************************************

ZESG01(:,:) = exp(0.125*PSIG0(:,1:JPMODE)**2)            
ZESG04(:,:)  = ZESG01(:,:) ** 4            
ZESG05(:,:)  = ZESG04(:,:) * ZESG01(:,:)                        
ZESG08(:,:)  = ZESG04(:,:) * ZESG04(:,:)                        
ZESG09(:,:)  = ZESG04(:,:) * ZESG05(:,:)            
ZESG12(:,:)  = ZESG04(:,:) * ZESG04(:,:) * ZESG04(:,:)
ZESG16(:,:)  = ZESG08(:,:) * ZESG08(:,:)
ZESG20(:,:)  = ZESG16(:,:) * ZESG04(:,:)
ZESG24(:,:)  = ZESG12(:,:) * ZESG12(:,:)
ZESG25(:,:)  = ZESG16(:,:) * ZESG09(:,:)
ZESG28(:,:)  = ZESG20(:,:) * ZESG08(:,:)
ZESG32(:,:)  = ZESG16(:,:) * ZESG16(:,:)
ZESG36(:,:)  = ZESG16(:,:) * ZESG20(:,:)
ZESG49(:,:)  = ZESG25(:,:) * ZESG20(:,:) * ZESG04(:,:)
ZESG52(:,:)  = ZESG16(:,:) * ZESG36(:,:)
ZESG64(:,:)  = ZESG32(:,:) * ZESG32(:,:)
ZESG81(:,:)  = ZESG49(:,:) * ZESG32(:,:)
ZESG85(:,:)  = ZESG64(:,:) * ZESG20(:,:) * ZESG01(:,:)
ZESG100(:,:) = ZESG36(:,:) * ZESG64(:,:) 
ZESG121(:,:) = ZESG85(:,:) * ZESG36(:,:)
ZESG144(:,:) = ZESG100(:,:) * ZESG36(:,:) * ZESG08(:,:)
ZESG169(:,:) = ZESG144(:,:) * ZESG25(:,:)
ZESG196(:,:) = ZESG144(:,:) * ZESG52(:,:)
ZESG256(:,:) = ZESG144(:,:) * ZESG100(:,:) * ZESG12(:,:)

!***************************************************************
! Transfert de moments entre les modes i et j
!***************************************************************

ZINVSIG(:,:)=1./PSIG0(:,1:JPMODE)**2
ZLNDG(:,:)=log(2.*PRG0(:,1:JPMODE))

ZA(:)=0.5*(ZINVSIG(:,1)-ZINVSIG(:,2))
ZD(:) = 0.
ZXi(:)= 0.

WHERE (ABS(ZA(:)) > 1E-4)
  ZB(:)=ZINVSIG(:,2)*ZLNDG(:,2)-ZINVSIG(:,1)*ZLNDG(:,1)
  ZC(:)=0.5*(ZINVSIG(:,1)*ZLNDG(:,1)**2-ZINVSIG(:,2)*ZLNDG(:,2)**2) - &
       &log((PN0(:,1)*PSIG0(:,2))/(PN0(:,2)*PSIG0(:,1)))

  ZD(:)=ZB(:)**2-4.*ZA(:)*ZC(:)

  ZSOL(:)=(-ZB(:)+sqrt(ABS(ZD(:))))/(2.*ZA(:))
  WHERE (ZSOL(:) < 5.E+2)
    ZSOL(:)=exp(ZSOL(:))/2.
    ZXi(:)=log(ZSOL(:)/PRG0(:,1))/(sqrt(2.)*PSIG0(:,1))
  ENDWHERE
ENDWHERE

!*********************************************************************
!      calculate the intramodal moment coefficients (log-normal model)
!*********************************************************************
       
do JI=1,JPMODE

  ZKFM(:)=sqrt(3.*XBOLTZ*PTGAS(:)/PRHOP0(:,JI))*1.e-3
  !ZKTURB(:)=sqrt(XPI*ZTURBDS*PMU(:)/(120.*PRHOP0(:,JI)))*1.e-18
  !ZKGRAV(:)=1.5/4.*0.544*XPI*PRHOP0(:,JI)/PMU(:)*1.e-24
!*************************************************************
!      calculate ZVG,ln2(sigma) and sigma
!      (log-normal model)
!*************************************************************

  ZRG(:)=PRG0(:,JI)
  ZKNG(:)=PLAMBDA(:)/ZRG(:)
  ZAI(:)=1.392*ZKNG(:)**0.0783

!***********************
! Brownian Coagulation  
!***********************
       
  ZRB0(:)=0.8
  ZRB6(:)=ZRB0
       
  ZINTRA0FM(:)=ZKFM(:)*ZRB0(:)*sqrt(2.*ZRG(:))*(ZESG01(:,JI)+ZESG25(:,JI)+2.*ZESG05(:,JI))
  ZINTRA3FM(:)=ZKFM(:)*ZRB0(:)*sqrt(ZRG(:))**7*sqrt(2.)*(ZESG49(:,JI)+ZESG36(:,JI)*ZESG01(:,JI)+&
              &2.*ZESG25(:,JI)*ZESG04(:,JI)+ZESG09(:,JI)*ZESG16(:,JI)+ZESG100(:,JI)*ZESG09(:,JI)+&
              &2.*ZESG64(:,JI)*ZESG01(:,JI))
  ZINTRA6FM(:)=ZKFM(:)*ZRB6(:)*sqrt(ZRG(:))**13*sqrt(2.)*ZESG85(:,JI)*&
              (1.+2.*ZESG04(:,JI)+ZESG24(:,JI))
  ZINTRA0NC(:)=ZKNC(:)*(1.+ZESG08(:,JI)+ZAI(:)*ZKNG(:)*(ZESG20(:,JI)+ZESG04(:,JI)))
  ZINTRA3NC(:)=ZKNC(:)*ZRG(:)**3*(2.*ZESG36(:,JI)+ZAI(:)*ZKNG(:)*(ZESG16(:,JI)+ZESG04(:,JI)*ZESG04(:,JI)+&
              &ZESG36(:,JI)*ZESG04(:,JI)+ZESG64(:,JI)*ZESG16(:,JI))+ZESG16(:,JI)*ZESG04(:,JI)+&
              &ZESG64(:,JI)*ZESG04(:,JI))
  ZINTRA6NC(:)=2.*ZKNC(:)*(ZRG(:))**6*ZESG52(:,JI)*(ZESG20(:,JI)+ZESG28(:,JI)+ZAI(:)*ZKNG(:)*(1.+ZESG16(:,JI)))
  ZINTRA0(:)=ZINTRA0FM(:)*(ZINTRA0NC(:)/(ZINTRA0FM(:)+ZINTRA0NC(:)))
  ZINTRA3(:)=ZINTRA3FM(:)*(ZINTRA3NC(:)/(ZINTRA3FM(:)+ZINTRA3NC(:)))
  ZINTRA6(:)=ZINTRA6FM(:)*(ZINTRA6NC(:)/(ZINTRA6FM(:)+ZINTRA6NC(:)))
   
  PDMINTRA(:,NM0(JI))=ZINTRA0(:)
  PDMINTRA(:,NM3(JI))=ZINTRA3(:)
  PDMINTRA(:,NM6(JI))=ZINTRA6(:)
  !print*,'PDMINTRA(:,NM0(',JI,') =',MINVAL(PDMINTRA(:,NM0(JI))), MAXVAL(PDMINTRA(:,NM0(JI)))
  !print*,'PDMINTRA(:,NM3(',JI,') =',MINVAL(PDMINTRA(:,NM3(JI))), MAXVAL(PDMINTRA(:,NM3(JI)))
  !print*,'PDMINTRA(:,NM6(',JI,') =',MINVAL(PDMINTRA(:,NM6(JI))), MAXVAL(PDMINTRA(:,NM6(JI)))

enddo
!print*,'=============================='
!print*,'=============================='

WHERE (ZD(:) > 0. .AND. ZXi(:) > (6.*PSIG0(:,1)/sqrt(2.)))

! transfert du moment d'ordre 0 (nombre)
!**************************************

  ZERF0(:)=sqrt(1.-exp(-4.*(ZXi(:))**2/XPI))
  ZPHI0(:)=0.5*(1.+ZERF0(:))

! transfert du moment d'ordre 3 (masse)
!**************************************

  ZERF3(:)=sqrt(1.-exp(-4.*(ZXi(:)-3.*PSIG0(:,1)/sqrt(2.))**2/XPI))
  ZPHI3(:)=0.5*(1.+ZERF3(:))
  
! transfert du moment d'ordre 6 (dispersion)
!**************************************

  ZERF6(:)=sqrt(1.-exp(-4.*(ZXi(:)-6.*PSIG0(:,1)/sqrt(2.))**2/XPI))
  ZPHI6(:)=0.5*(1.+ZERF6(:))
  
  PDMINTRA(:,NM0(2))=PDMINTRA(:,NM0(2))-(1.-ZPHI0(:)**2)*PDMINTRA(:,NM0(1))*(PM(:,NM0(1))/PM(:,NM0(2)))**2
  PDMINTRA(:,NM0(1))=(2.-ZPHI0(:)**2)*PDMINTRA(:,NM0(1))

  PDMINTRA(:,NM3(2))=PDMINTRA(:,NM3(1))*(1.-ZPHI0(:)*ZPHI3(:))*PM(:,NM0(1))**2
  PDMINTRA(:,NM3(1))=PDMINTRA(:,NM3(1))*(ZPHI0(:)*ZPHI3(:)-1.)*PM(:,NM0(1))**2
  
  ZKFM(:)=sqrt(3.*XBOLTZ*PTGAS(:)/PRHOP0(:,1))*1.e-3
  ZRG(:)=PRG0(:,1)
  ZKNG(:)=PLAMBDA(:)/ZRG(:)
  ZAI(:)=1.392*ZKNG(:)**0.0783
  
  ZINTRA6FM(:)=ZKFM(:)*sqrt(2.)*sqrt(ZRG(:))**13*(ZESG169(:,1)+ZESG144(:,1)*ZESG01(:,1)+&
               2.*ZESG121(:,1)*ZESG04(:,1)+ZESG81(:,1)*ZESG16(:,1)+&
               ZESG256(:,1)*ZESG09(:,1)+ZESG196(:,1)*ZESG01(:,1))

  ZINTRA6NC(:)=ZKNC(:)*(ZRG(:))**6*(2.*ZESG144(:,1)+ZAI(:)*ZKNG(:)*(ZESG100(:,1)+&
             ZESG64(:,1)*ZESG04(:,1))+ZAI(:)*ZKNG(:)*(ZESG144(:,1)*ZESG04(:,1)+&
             ZESG196(:,1)*ZESG16(:,1))+ZESG100(:,1)*ZESG04(:,1)+&
             ZESG196(:,1)*ZESG04(:,1))
             
  ZINTRA6(:)=ZINTRA6FM(:)*(ZINTRA6NC(:)/(ZINTRA6FM(:)+ZINTRA6NC(:)))

  PDMINTRA(:,NM6(2))=PDMINTRA(:,NM6(2))+(PDMINTRA(:,NM6(1))*(1.-ZPHI3(:)**2)+ZINTRA6(:)*(1.-ZPHI0(:)*ZPHI6(:)))&
                    &*(PM(:,NM0(1))/PM(:,NM0(2)))**2
                   
  PDMINTRA(:,NM6(1))=PDMINTRA(:,NM6(1))*(ZPHI3(:)**2)+ZINTRA6(:)*(ZPHI0(:)*ZPHI6(:)-1.)
  
ELSEWHERE

  PDMINTRA(:,NM3(1))=0.
  PDMINTRA(:,NM3(2))=0.


ENDWHERE

do JI=1,JPMODE
!print*,'2.-ZPHI0(:)**2 =',MINVAL(2.-ZPHI0(:)**2), MAXVAL(2.-ZPHI0(:)**2)
!  print*,'apres corr PDMINTRA(:,NM0(',JI,') =',MINVAL(PDMINTRA(:,NM0(JI))), MAXVAL(PDMINTRA(:,NM0(JI)))
!  print*,'apres corr PDMINTRA(:,NM3(',JI,') =',MINVAL(PDMINTRA(:,NM3(JI))), MAXVAL(PDMINTRA(:,NM3(JI)))
!  print*,'apres corr PDMINTRA(:,NM6(',JI,') =',MINVAL(PDMINTRA(:,NM6(JI))), MAXVAL(PDMINTRA(:,NM6(JI)))
 enddo

!*********************************************************************
!   calculate the intermodal moment coefficients (log-normal model)
!*********************************************************************

do JI=1,(JPMODE-1)
  do JJ=(JI+1),JPMODE

    ZRGI(:)=PRG0(:,JI)
    ZKNGI(:)=PLAMBDA(:)/ZRGI(:)
    ZAI(:)=1.392*ZKNGI(:)**0.0783

    ZRGJ(:)=PRG0(:,JJ)
    ZKNGJ(:)=PLAMBDA(:)/ZRGJ(:)
    ZAJ(:)=1.392*ZKNGJ(:)**0.0783
          
    ZR(:)=sqrt(ZRGJ(:)/ZRGI(:))
    ZR2(:)=ZR(:)*ZR(:)
    !ZR3(:)=ZR(:)*ZR2(:)
    ZR4(:)=ZR2(:)*ZR2(:)
    ZRM(:)=1./ZR(:)
    ZRM2(:)=ZRM(:)*ZRM(:)
    ZRM3(:)=ZRM(:)*ZRM2(:)
    !ZRM4(:)=ZRM2(:)*ZRM2(:)

!**********************
! Brownian Coagulation
!**********************

      ZRES(:)=0.9

      ZAPPROX(:)=sqrt(2.*ZRGI(:))*(ZESG01(:,JI)+ZR(:)*ZESG01(:,JJ)+2.*ZR2(:)*ZESG01(:,JI)*ZESG04(:,JJ)&
                 +ZR4(:)*ZESG09(:,JI)*ZESG16(:,JJ)+ZRM3(:)*ZESG16(:,JI)*ZESG09(:,JJ)+&
                 2.*ZRM(:)*ZESG04(:,JI)*ZESG01(:,JJ))
       
      ZINTERFM(:)=ZKFM(:)*ZRES(:)*ZAPPROX(:)

      ZAPPROX(:)=2.+ZAI(:)*ZKNGI(:)*(ZESG04(:,JI)+ZR2(:)*ZESG16(:,JI)*ZESG04(:,JJ))+&
                ZAJ(:)*ZKNGJ(:)*(ZESG04(:,JJ)+ZRM2(:)*ZESG16(:,JJ)*ZESG04(:,JI))+&
                (ZR2(:)+ZRM2(:))*(ZESG04(:,JI)*ZESG04(:,JJ))

      ZINTERNC(:)=ZKNC(:)*ZAPPROX(:)

      ZINTER(:)=ZINTERNC(:)*(ZINTERFM(:)/(ZINTERNC(:)+ZINTERFM(:)))

      PDMINTER(:,NM0(JI))=PM(:,NM0(JJ))*ZINTER(:)
      PDMINTER(:,NM0(JJ))=-PM(:,NM0(JJ))*ZINTER(:)

      ZAPPROX(:)=sqrt(2.)*sqrt(ZRGI(:))**7*(ZESG49(:,JI)+ZR(:)*ZESG36(:,JI)*ZESG01(:,JJ)+2.*ZR2(:)*&
                 ZESG25(:,JI)*ZESG04(:,JJ)+ZR4(:)*ZESG09(:,JI)*ZESG16(:,JJ)+ZRM3(:)*&
                 ZESG100(:,JI)*ZESG09(:,JJ)+2.*ZRM(:)*ZESG64(:,JI)*ZESG01(:,JJ))

      ZINTERFM(:)=ZKFM(:)*ZRES(:)*ZAPPROX(:)

      ZAPPROX(:)=(2.*ZESG36(:,JI)+ZAI(:)*ZKNGI(:)*(ZESG16(:,JI)+ZR2(:)*ZESG04(:,JI)*ZESG04(:,JJ))+&
      ZAJ(:)*ZKNGJ(:)*(ZESG36(:,JI)*ZESG04(:,JJ)+ZRM2(:)*ZESG16(:,JJ)*ZESG64(:,JI))+&
      ZR2(:)*ZESG16(:,JI)*ZESG04(:,JJ)+ZRM2(:)*ZESG64(:,JI)*ZESG04(:,JJ))*(ZRGI(:))**3      

      ZINTERNC(:)=ZKNC(:)*ZAPPROX(:)

      ZINTER(:)=ZINTERNC(:)*(ZINTERFM(:)/(ZINTERNC(:)+ZINTERFM(:)))
       
      PDMINTER(:,NM3(JI))=-PM(:,NM0(JI))*PM(:,NM0(JJ))*ZINTER(:)
      PDMINTER(:,NM3(JJ))=PM(:,NM0(JI))*PM(:,NM0(JJ))*ZINTER(:)

      ZAPPROX(:)=sqrt(2.)*sqrt(ZRGI(:))**13*(ZESG169(:,JI)+ZR(:)*ZESG144(:,JI)*ZESG01(:,JJ)+&
             2.*ZR2(:)*ZESG121(:,JI)*ZESG04(:,JJ)+ZR4(:)*ZESG81(:,JI)*ZESG16(:,JJ)+&
             ZRM3(:)*ZESG256(:,JI)*ZESG09(:,JJ)+2*ZRM(:)*ZESG196(:,JI)*ZESG01(:,JJ))
      
      ZINTERFM(:)=ZKFM(:)*ZRES(:)*ZAPPROX(:)

      ZAPPROX(:)=(ZRGI(:))**6*(2.*ZESG144(:,JI)+ZAI(:)*ZKNGI(:)*(ZESG100(:,JI)+&
       ZR2(:)*ZESG64(:,JI)*ZESG04(:,JJ))+ZAJ(:)*ZKNGJ(:)*(ZESG144(:,JI)*ZESG04(:,JJ)+&
       ZRM2(:)*ZESG196(:,JI)*ZESG16(:,JJ))+ZR2(:)*ZESG100(:,JI)*ZESG04(:,JJ)+&
       ZRM2(:)*ZESG196(:,JI)*ZESG04(:,JJ))

      ZINTERNC(:)=ZKNC(:)*ZAPPROX(:)

      ZINTER(:)=ZINTERNC(:)*(ZINTERFM(:)/(ZINTERNC(:)+ZINTERFM(:)))

      PDMINTER(:,NM6(JI))=-PM(:,NM0(JI))*PM(:,NM0(JJ))*ZINTER(:)
      PDMINTER(:,NM6(JJ))=PM(:,NM0(JI))*PM(:,NM0(JJ))*ZINTER(:)

      ZAPPROX(:)=sqrt(2.)*sqrt(ZRGI(:))**7*sqrt(ZRGJ(:))**6*(ZESG49(:,JI)*&
           ZESG36(:,JJ)+ZR(:)*ZESG36(:,JI)*ZESG49(:,JJ)+2.*ZR2(:)*ZESG25(:,JI)*&
           ZESG64(:,JJ)+ZR4(:)*ZESG09(:,JI)*ZESG100(:,JJ)+ZRM3(:)*ZESG100(:,JI)*&
           ZESG09(:,JJ)+2.*ZRM(:)*ZESG64(:,JI)*ZESG25(:,JJ))
       
      ZINTERFM(:)=ZKFM(:)*ZRES(:)*ZAPPROX(:)

      ZAPPROX(:)=(ZRGI(:))**3*(ZRGJ(:))**3*(2.*ZESG36(:,JI)*ZESG36(:,JJ)+&
       ZAI(:)*ZKNGI(:)*(ZESG16(:,JI)*ZESG16(:,JJ)+ZR2(:)*ZESG04(:,JI)*ZESG64(:,JJ))+&
       ZAJ(:)*ZKNGJ(:)*(ZESG36(:,JI)*ZESG16(:,JJ)+ZRM2(:)*ZESG64(:,JI)*ZESG04(:,JJ))+&
       ZR2(:)*ZESG16(:,JI)*ZESG64(:,JJ)+ZRM2(:)*ZESG64(:,JI)*ZESG16(:,JJ))

      ZINTERNC(:)=ZKNC(:)*ZAPPROX(:)

      ZINTER(:)=ZINTERNC(:)*(ZINTERFM(:)/(ZINTERNC(:)+ZINTERFM(:)))
       
      PDMINTER(:,NM6(JJ))=PDMINTER(:,NM6(JJ))+2.*PM(:,NM0(JI))*PM(:,NM0(JJ))*ZINTER(:)
 
  enddo
enddo

IF (LHOOK) CALL DR_HOOK('CH_AER_COAG',1,ZHOOK_HANDLE)
END SUBROUTINE CH_AER_COAG 
