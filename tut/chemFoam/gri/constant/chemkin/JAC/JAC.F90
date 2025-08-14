MODULE example1
CONTAINS
SUBROUTINE FEX(NEQ, T, Y, YDOT,RWORK,IWORK)
!SUBROUTINE FEX(NEQ, T, Y, YDOT)
!     IMPLICIT NONE
!     INTEGER NEQ
!     DOUBLE PRECISION T, Y, YDOT
!     DIMENSION Y(NEQ), YDOT(NEQ)
        IMPLICIT NONE
        INTEGER, INTENT (IN) :: NEQ
        DOUBLE PRECISION, INTENT (IN) :: T
        DOUBLE PRECISION, INTENT (IN) :: Y(NEQ)
        DOUBLE PRECISION, INTENT (OUT) :: YDOT(NEQ)
         DOUBLE PRECISION, INTENT (INOUT) :: RWORK,IWORK
        DIMENSION RWORK(*),IWORK(*)
     YDOT(1) = -.04D0*Y(1) + 1.D4*Y(2)*Y(3)
     YDOT(3) = 3.E7*Y(2)*Y(2)
     YDOT(2) = -YDOT(1) - YDOT(3)
     RETURN
     END SUBROUTINE FEX
     
SUBROUTINE JEX(NEQ, T, Y, ML, MU, PD, NRPD)
     IMPLICIT NONE
     INTEGER , INTENT (IN) ::NEQ,ML,MU,NRPD
!     DOUBLE PRECISION PD
!     DOUBLE PRECISION, INTENT (IN) :: T, Y(NEQ)
!     DIMENSION PD(NRPD,NEQ)
        DOUBLE PRECISION, INTENT (IN) :: T
        DOUBLE PRECISION, INTENT (IN) :: Y(NEQ)
        DOUBLE PRECISION, INTENT (OUT) :: PD(NRPD,NEQ)
     PD(1,1) = -.04D0
     PD(1,2) = 1.D4*Y(3)
     PD(1,3) = 1.D4*Y(2)
     PD(2,1) = .04D0
     PD(2,3) = -PD(1,3)
     PD(3,2) = 6.E7*Y(2)
     PD(2,2) = -PD(1,2) - PD(3,2)
     RETURN
     END SUBROUTINE JEX
     
 END MODULE example1
! ******************************************************************

     PROGRAM runexample1
     USE example1
!     IMPLICIT NONE
    IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER(I-N)
     parameter (LENIWK=500000, LENRWK=500000, LENCWK=50000, NK=5, NLMAX=55)
     parameter (LIN=5, LOUT=6,LINKCK=25,NDIM=10,SMALL=1.D-200)
     parameter (KMAX=500,IMAX=5000, ITOL=1, IOPT=0, ATOL=1.0E-10, RTOL=1.0E-6)
     INTEGER NEQ, ITASK, ISTATE, ISTATS, IOUT, IERROR, I
     DIMENSION Y(3), RSTATS(22), ISTATS(31),YDOT(3),KI(NDIM),NU(NDIM)
     DIMENSION RWORK(LENRWK),IWORK(LENIWK),X(KMAX), Z(KMAX),ITHB(IMAX), IREV(IMAX)
     Dimension NUKI(KMAX,IMAX),PD(KMAX,IMAX),prod(KMAX,IMAX),dest(KMAX,IMAX),NSPinREV(IMAX*NDIM),&
     NSPinREV2(IMAX*NDIM),IFOP(IMAX), KFAL(IMAX), FPAR(NDIM,IMAX),AKI(KMAX,IMAX),&
     KNReacSpec(NDIM),KNProdSpec(NDIM),KNUProdSpec(KMAX),KNUReacSpec(KMAX),NSPinTH(NDIM*IMAX),NSPinTH2(KMAX),&
     IcRVIcRV(IMAX),IcLTIcLT(IMAX),IcRVNRV(IMAX),IcLTNLT(IMAX)
     CHARACTER*16 CWORK(LENCWK),KSYM(KMAX)
     CHARACTER*50 REAC(IMAX),ISTR
     LOGICAL KERR, IERR
    double precision, save :: atol_save = 1.0E-15, rtol_save = 1.0E-6
    double precision :: gas_pressure,gas_temperature,endtime,dtime,MOL_H2,MOL_O2,time,tout,MOL_N2
    character(*), parameter :: file_out = 'out.f90'
 !   DIMENSION Y(3), RSTATS(22), ISTATS(31),YDOT(3)
!     OPEN(UNIT=6, FILE = 'example1.dat')
COMMON /ICONS/ KK, NP, NWT, NH, NWDOT
OPEN (unit=LINKCK, file = 'chem.bin',form ='UNFORMATTED')
!write(*,'(200A20)')'subroutine','CKJAC','ABC'
open(21,file=file_out)
WRITE (21, 7100)
      
      CALL CKLEN  (LINKCK, LOUT, LENI, LENR, LENC)
      CALL CKINIT (LENIWK, LENRWK, LENCWK, LINKCK, LOUT, IWORK, RWORK, CWORK)
      CALL CKINDX (IWORK, RWORK, MM, KK, II, NFIT)

      NEQ   = KK + 1           ! 方程总数
      LRW   = 22 + 9*NEQ + 2*NEQ**2  !                                            应该是存放热力学数据
      NVODE = LENR + 1              ! LENR是实数数组RWORK所需要的最少长度
      NP    = NVODE + LRW         ! 存放压力
      NWT   = NP + 1              ! 存放组分的摩尔质量
      NH    = NWT  + KK       ! 存放H  焓值
      NWDOT = NH   + KK       ! 存放wt，组分k的产生速率
      NTOT  = NWDOT+ KK - 1
!
      LIW   = 30 + NEQ
      IVODE = LENI + 1
      ITOT  = IVODE + LIW - 1
      
      CALL CKSYMS (CWORK, LOUT, KSYM, IERR)
      
      IF (KK.GT.KMAX .OR. LENRWK.LT.NTOT .OR. LENIWK.LT.ITOT) THEN
         IF (KK .GT. KMAX)  WRITE (LOUT, *)' Error...KMAX too small...must be at least ', KK
         IF (LENRWK .LT. NTOT) WRITE (LOUT, *)' Error...LENRWK too small...must be at least', NTOT
         IF (LENIWK .LT. ITOT) WRITE (LOUT, *)' Error...LENIWK too small...must be at least', ITOT
         STOP
      ENDIF
      IF (IERR) KERR = .TRUE.
      
      CALL CKWT   (IWORK, RWORK, RWORK(NWT))
      CALL CKRP   (IWORK, RWORK, RU, RUC, PATM)
      
      CALL CKNU (KK, IWORK, RWORK, NUKI(1:KK,1:II))

write(*,*)'ksym = ',ksym
!write(*,*)'NUKI = ',NUKI(1:KK,1:II)
call CKcommon (IWORK, RWORK,NNREV,IIcRV,NNcTT,NNcAA,NNCP2,NNCP2T,NNPAR,NNcCO,NNFAR,NNcFL,NNcLT,&
                IIcLT,NNcRL,IIcRL,NNLAN,NNLAR,NNRLT,NNcRV)

NRV = 1
NLT = 1
do i=1,II
    IcRVIcRV(i) = 0
    IcLTIcLT(i) = 0
    IcRVNRV(i) = 0
    IcLTNLT(i) = 0
    if (i == IWORK(IIcRV+NRV-1) .AND. NRV <= NNREV)then
!   if reaction i with explicit reverse Arrhenius coefficients
        IcRVIcRV(i)=i
        IcRVNRV(i) = NRV
        NRV = NRV+1
    endif
    if (i == IWORK(IIcLT+NLT-1) .AND. NLT <= NNLAN)then
!   if reaction i is Landau-Teller reaction
        IcLTIcLT(i)=i
        IcLTNLT(i)=NLT
        NLT = NLT+1
    endif
enddo

! SUBROUTINE pTY2rhoC (P, T, Y, RHO, C)
write(21,100)
k=1
do k=1,KK
write(21,10000)ksym(k)
write(21,10001)k,k,RWORK(K+MM)
enddo

write(21,10002)KK,KK
write(21,101)
! SUBROUTINE pTY2rhoC (P, T, Y, RHO, C)


! SUBROUTINE pTC2rhoY (p, T, C, RHO, Y)
write(21,119)KK
write(21,20002)KK,KK
write(21,101)
! SUBROUTINE pTC2rhoY (p, T, C, RHO, Y)


! SUBROUTINE pTY2rhoC_OF (P, T, Y, RHO, C)
write(21,122)
k=1
do k=1,KK
write(21,10000)ksym(k)
write(21,10001)k,k,RWORK(K+MM)
enddo

write(21,30002)KK,KK
write(21,101)
! SUBROUTINE pTY2rhoC_OF (P, T, Y, RHO, C)


! SUBROUTINE pTC2rhoY_OF (p, T, C, RHO, Y)
write(21,123)KK
write(21,30003)KK,KK
write(21,101)
! SUBROUTINE pTC2rhoY_OF (p, T, C, RHO, Y)


!
!! SUBROUTINE rhoTY2C (RHO, T, Y, C)
!write(21,102)
!k=1
!do K=1,KK
!write(21,10000)ksym(k)
!write(21,10001)k,k,RWORK(K+MM)
!enddo
!write(21,10003)KK
!write(21,101)
!! SUBROUTINE rhoTY2C (RHO, T, Y, C)
!
!

!
!! SUBROUTINE SKSMH  (T, SMH)
!! NDIM = 10 == MAXSP in cklib.f
!! IIcRV 反应中 给定逆反应速率的反应标号
!! NNcTT 
!write(21,103)
!
call CKITR (IWORK, RWORK, ITHB(1:II), IREV(1:II))
! ITHB: third-body indices for reactions
! IREV: reversibility indices and species count
!
!NSP=0
!do i=1,II
!    ! IF reversible reaction
!    if (IREV(i) > 0)then
!        CALL CKINU (i, NDIM, IWORK, RWORK, NSPEC, KI, NU)
!! KI(N) is the index of the Nth species in reaction I
!! NU(N) is the stoichiometric coefficient of the Nth species in reaction I
!        DO K=1,NSPEC
!            NSP = NSP+1
!            NSPinREV(NSP)=KI(K)
!        enddo
!    endif
!enddo
!!    write(*,*)NSPinREV(1:NSP)
!    call delSameNumbANDsort(NSP,NSPinREV,NNSP,NSPinREV2)
!!    write(*,*)NSPinREV2(1:NNSP)
!! NNSP is the total species involving in reversible reaction
!DO K=1,NNSP
!!    write(21,10000)ksym(k)
!    write(21,10000)ksym(NSPinREV2(k))
!    TEMP = RWORK(NNcTT + (NSPinREV2(k)-1)*3 + 1)
!    
!    L=2
!    NA1 = NNcAA + (L-1)*NNCP2 + (NSPinREV2(k)-1)*NNCP2T
!    WRITE(21,10004)TEMP,NSPinREV2(K),RWORK(NA1+6),&
!    -RWORK(NA1+5),RWORK(NA1),RWORK(NA1+1)/2.D0,RWORK(NA1+2)/6.D0,&
!    RWORK(NA1+3)/12.D0,RWORK(NA1+4)/20.D0
!    
!    L=1
!    NA1 = NNcAA + (L-1)*NNCP2 + (NSPinREV2(k)-1)*NNCP2T
!    WRITE(21,10005)NSPinREV2(K),RWORK(NA1+6),&
!    -RWORK(NA1+5),RWORK(NA1),RWORK(NA1+1)/2.D0,RWORK(NA1+2)/6.D0,&
!    RWORK(NA1+3)/12.D0,RWORK(NA1+4)/20.D0
!enddo
!
!write(21,101)
!! SUBROUTINE SKSMH  (T, SMH)
!


!
!! SUBROUTINE SKRATT (T, RF, RB, RKLOW)
!write(21,104)KK,KK
!DO K=1,NNSP
!    write(21,10006)K,K
!ENDDO
!write(21,10007)
!
!do i=1,II
!    NRFArru = NNcCO + (i-1)*(NNPAR+1)
!    CALL CKINU (i, NDIM, IWORK, RWORK, NSPEC, KI, NU)
!    call CKSYMR (i, LOUT, IWORK, RWORK, CWORK, LT, ISTR,KERR)
!    write(21,10017)i,ISTR
!! RF
!!    IF (i == IWORK(IIcLT+LT-1) .AND. LT <= NNLAN)THEN ! IF Landau-Teller reactions
!    IF (IcLTIcLT(i) > 0)THEN ! IF Landau-Teller reactions
!        NRFLT = NNcLT + (IcLTNLT(i)-1)*NNLAR ! NNLAR = 4
!        if(RWORK(NRFArru) > SMALL)then
!            WRITE(21,10037)i,log(RWORK(NRFArru)),RWORK(NRFArru+1),-RWORK(NRFArru+2),& ! RF
!                        RWORK(NRFLT),RWORK(NRFLT+1)
!        else
!            write(21,10133)i
!        endif
!    ELSE
!        if(RWORK(NRFArru) > SMALL)then
!            WRITE(21,10008)i,log(RWORK(NRFArru)),RWORK(NRFArru+1),-RWORK(NRFArru+2)! RF ! RWORK(NRFArru+2) = Ea/1.987 (unit K) 
!        else
!            write(21,10133)i
!        endif
!    ENDIF
!!    if ((RWORK(NRFArru+1) .NE. 0.D0) .AND. (RWORK(NRFArru+2) .NE. 0.D0))then
!!      WRITE(21,10008)i,log(RWORK(NRFArru)),RWORK(NRFArru+1),-RWORK(NRFArru+2)
!!    endif
!!    if ((RWORK(NRFArru+1) == 0.D0) .AND. (RWORK(NRFArru+2) .NE. 0.D0))then
!!      WRITE(21,10010)i,log(RWORK(NRFArru)),-RWORK(NRFArru+2)
!!    endif
!!    if ((RWORK(NRFArru+1) .NE. 0.D0) .AND. (RWORK(NRFArru+2) == 0.D0))then
!!      WRITE(21,10012)i,log(RWORK(NRFArru)),RWORK(NRFArru+1)
!!    endif
!!    if ((RWORK(NRFArru+1) == 0.D0) .AND. (RWORK(NRFArru+2) == 0.D0))then
!!      WRITE(21,10012)i,log(RWORK(NRFArru))
!!    endif
!! NNLAR== 0
!! RB
!! IF reversible reaction
!if (IREV(i) > 0)then
!    if (NNRLT>0) then
!        write(21,*)'**ERROR--need to supply reverse Landau-Teller reactions'
!    endif
!!    IF (i == IWORK(IIcRV+LIcRV-1) .AND. LIcRV <= NNREV) THEN ! if reaction i with explicit reverse Arrhenius coefficients
!    IF (IcRVIcRV(i) > 0) THEN ! if reaction i with explicit reverse Arrhenius coefficients
!        NRBArru = NNcRV + (IcRVNRV(i)-1)*(NNPAR+1) ! NNPAR = 3
!        if(RWORK(NRBArru) > SMALL)then
!            WRITE(21,10009)i,log(RWORK(NRBArru)),RWORK(NRBArru+1),-RWORK(NRBArru+2)
!        else
!            WRITE(21,10132)i
!        endif
!    ELSE ! if reaction i without explicit reverse Arrhenius coefficients
!        write(21,10021,advance='no') ! EQK
!        NUSUM = 0
!        DO K=1,NSPEC
!            WRITE(21,10022,advance='no')NUKI(KI(K),I),KI(K)
!            NUSUM = NUSUM + NU(K)
!        ENDDO
!        WRITE(21,10124)NUSUM
!      WRITE(21,10016)i,i ! RB(X) = RF(X) / MAX(EQK, SMALL)
!    ENDIF
!else ! IF not reversible reaction
!    WRITE(21,10049)i ! RB(X) = 0.D0
!endif
!
!enddo
!
!WRITE(21,*)
!
call CKFAL  (NDIM, IWORK, RWORK, IFOP, KFAL, FPAR)
! FPAR(1,I), FPAR(2,I), FPAR(3,I) are always the parameters entered on the LOW auxiliary keyword line
!do i=1,II
!    IF (ITHB(i) >= 0 .AND. IFOP(i) > 0) THEN ! if three body reaction and pressure depend reaction
!!        NthreeBody = NNcFL + (i-1)*NNFAR
!!        WRITE(21,10023)i,log(RWORK(NthreeBody)),RWORK(NthreeBody+1),-RWORK(NthreeBody+2)
!        WRITE(21,10023)i,log(FPAR(1,i)),FPAR(2,i),-FPAR(3,i)
!    ENDIF
!enddo
!
!write(21,101)
!! SUBROUTINE SKRATT (T, RF, RB, RKLOW)
!


!
!! SUBROUTINE SKRATX (T, C, RF, RB, RKLOW)
!write(21,105)KK
call CKTHB (KK, IWORK, RWORK, AKI(1:KK,1:II))
! AKI(K,I) is the enhanced efficiency of species K in reaction I
!do i=1,II
!    CALL CKINU (i, NDIM, IWORK, RWORK, NSPEC, KI, NU)
!! KI(N) is the index of the Nth species in reaction I
!! NU(N) is the stoichiometric coefficient of the Nth species in reaction I
!    call CKSYMR (i, LOUT, IWORK, RWORK, CWORK, LT, ISTR,KERR)
!    write(21,10017)i,ISTR
!    if (ITHB(i) >= 0 ) then! if three body reaction 
!        if ( IFOP(i) > 0) then! if pressure depend reaction
!            write(21,10024,advance='no')
!            KKTHBinI = 0
!            DO K=1,KK
!                IF (AKI(k,i) .NE. 1.D0)THEN
!                    KKTHBinI = KKTHBinI + 1
!                    IF (MOD(KKTHBinI,5) == 0)THEN
!                        write(21,10127)
!                    ENDIF
!                    write(21,10025,advance='no')AKI(k,i)-1.D0,K ! CM=CTOT+[M]
!                ENDIF
!            ENDDO
!            WRITE(21,10031)i,i
!            if (IFOP(i) == 3) then
!                WRITE(21,10032)1-FPAR(4,i),FPAR(5,i),FPAR(4,i),FPAR(6,i)
!                write(21,10034)
!            elseif (IFOP(i) == 4) then
!                WRITE(21,10033)1-FPAR(4,i),FPAR(5,i),FPAR(4,i),FPAR(6,i),FPAR(7,i)
!                write(21,10034)
!            elseif (IFOP(i) == 2) then
!                write(21,*)'ERROR--need to supply SRI fomat output'
!            elseif (IFOP(i) == 1) then
!                write(21,10125)
!            endif
!            write(21,10035)i,i
!            write(21,10036)i,i
!            WRITE(21,10029,advance='no')I,I
!            DO K=1,NSPEC
!                IF(NU(K)<0)THEN
!                    WRITE(21,10027,advance='no')KI(K),ABS(NU(K))
!                ENDIF
!            ENDDO
!            WRITE(21,*)
!            WRITE(21,10030,advance='no')I,I
!            DO K=1,NSPEC
!                IF(NU(K)>0)THEN
!                    WRITE(21,10027,advance='no')KI(K),NU(K)
!                ENDIF
!            ENDDO
!            WRITE(21,*)
!        elseif (IFOP(i) == 0) then ! if not pressure depend reaction
!            write(21,10024,advance='no')
!            KKKTHBinI = 0
!            DO K=1,KK
!                IF (AKI(k,i) .NE. 1.D0)THEN
!                    KKKTHBinI = KKKTHBinI + 1
!                    IF (MOD(KKKTHBinI,5) == 0)THEN
!                        write(21,10127)
!                    ENDIF
!                    write(21,10025,advance='no')AKI(k,i)-1.D0,K ! CM=CTOT+[M]
!                ENDIF
!            ENDDO
!            WRITE(21,*)
!            WRITE(21,10026,advance='no')I,I
!            DO K=1,NSPEC
!                IF(NU(K)<0)THEN
!                    WRITE(21,10027,advance='no')KI(K),ABS(NU(K))
!                ENDIF
!            ENDDO
!            WRITE(21,*)
!            WRITE(21,10028,advance='no')I,I
!            DO K=1,NSPEC
!                IF(NU(K)>0)THEN
!                    WRITE(21,10027,advance='no')KI(K),NU(K)
!                ENDIF
!            ENDDO
!            WRITE(21,*)
!        endif
!    elseif (ITHB(i) == -1) then ! if not three body reaction 
!            WRITE(21,10029,advance='no')I,I
!            DO K=1,NSPEC
!                IF(NU(K)<0)THEN
!                    WRITE(21,10027,advance='no')KI(K),ABS(NU(K))
!                ENDIF
!            ENDDO
!            WRITE(21,*)
!            WRITE(21,10030,advance='no')I,I
!            DO K=1,NSPEC
!                IF(NU(K)>0)THEN
!                    WRITE(21,10027,advance='no')KI(K),NU(K)
!                ENDIF
!            ENDDO
!            WRITE(21,*)
!    endif
!    WRITE(21,*)
!enddo
!write(21,101)
!! SUBROUTINE SKRATX (T, C, RF, RB, RKLOW)
!

!
!! SUBROUTINE SKWDOT(RF, RB, WDOT)
!write(21,106)KK
!DO i=1,II
!    write(21,10038)i,i
!    CALL CKINU (i, NDIM, IWORK, RWORK, NSPEC, KI, NU)
!    DO K=1,NSPEC
!        IF(NU(K)>0)THEN
!            WRITE(21,10039)KI(K),KI(K),NU(K)
!        ELSE
!            WRITE(21,10040)KI(K),KI(K),ABS(NU(K))
!        ENDIF
!    ENDDO
!ENDDO
!write(21,101)
!! SUBROUTINE SKWDOT(RF, RB, WDOT)
!

! SUBROUTINE SKHML  (T, HML)
write(21,107)
DO K=1,KK
    write(21,10000)ksym(k)
    RU=8.31451D7
    TEMP = RWORK(NNcTT + (K-1)*3 + 1)
    
    L=2
    NA1 = NNcAA + (L-1)*NNCP2 + (K-1)*NNCP2T
    WRITE(21,10041)TEMP,K,RWORK(NA1+5)*RU,&
    RWORK(NA1)*RU,RWORK(NA1+1)/2.D0*RU,RWORK(NA1+2)/3.D0*RU,&
    RWORK(NA1+3)/4.D0*RU,RWORK(NA1+4)/5.D0*RU
    
    L=1
    NA1 = NNcAA + (L-1)*NNCP2 + (K-1)*NNCP2T
    WRITE(21,10042)K,RWORK(NA1+5)*RU,&
    RWORK(NA1)*RU,RWORK(NA1+1)/2.D0*RU,RWORK(NA1+2)/3.D0*RU,&
    RWORK(NA1+3)/4.D0*RU,RWORK(NA1+4)/5.D0*RU
ENDDO
write(21,101)
! SUBROUTINE SKHML  (T, HML)

! SUBROUTINE SKCPML  (T, CPML)
write(21,108)
DO K=1,KK
    write(21,10000)ksym(k)
    RU=8.31451D7
    TEMP = RWORK(NNcTT + (K-1)*3 + 1)
    
    L=2
    NA1 = NNcAA + (L-1)*NNCP2 + (K-1)*NNCP2T
    WRITE(21,10043)TEMP,K,RWORK(NA1)*RU,&
    RWORK(NA1+1)*RU,RWORK(NA1+2)*RU,&
    RWORK(NA1+3)*RU,RWORK(NA1+4)*RU
    
    L=1
    NA1 = NNcAA + (L-1)*NNCP2 + (K-1)*NNCP2T
    WRITE(21,10044)K,RWORK(NA1)*RU,&
    RWORK(NA1+1)*RU,RWORK(NA1+2)*RU,&
    RWORK(NA1+3)*RU,RWORK(NA1+4)*RU
ENDDO
write(21,101)
! SUBROUTINE SKCPML  (T, CPML)

! SUBROUTINE GDG1(T, G0, DG)
write(21,109)
DO K=1,KK
    write(21,10000)ksym(k)
    RU=8.31451D7
    TEMP = RWORK(NNcTT + (K-1)*3 + 1)
    
    L=2
    NA1 = NNcAA + (L-1)*NNCP2 + (K-1)*NNCP2T
    WRITE(21,10045)TEMP,K,-RWORK(NA1),-RWORK(NA1+1)/2.D0,-RWORK(NA1+2)/6.D0,&
    -RWORK(NA1+3)/12,-RWORK(NA1+4)/20,RWORK(NA1+5),&
    RWORK(NA1)-RWORK(NA1+6) ! (a1-a7)
    WRITE(21,10047)K,-RWORK(NA1),-RWORK(NA1+1)/2.D0,-RWORK(NA1+2)/3.D0,&
    -RWORK(NA1+3)/4,-RWORK(NA1+4)/5,-RWORK(NA1+5)
    
    L=1
    NA1 = NNcAA + (L-1)*NNCP2 + (K-1)*NNCP2T
    WRITE(21,10046)K,-RWORK(NA1),-RWORK(NA1+1)/2.D0,-RWORK(NA1+2)/6.D0,&
    -RWORK(NA1+3)/12,-RWORK(NA1+4)/20,RWORK(NA1+5),&
    RWORK(NA1)-RWORK(NA1+6)
    WRITE(21,10047)K,-RWORK(NA1),-RWORK(NA1+1)/2.D0,-RWORK(NA1+2)/3.D0,&
    -RWORK(NA1+3)/4,-RWORK(NA1+4)/5,-RWORK(NA1+5)
    WRITE(21,10094)
ENDDO
write(21,101)
! SUBROUTINE GDG1(T, G0, DG)


! SUBROUTINE SKG0(T, G0)
write(21,117)
DO K=1,KK
    write(21,10000)ksym(k)
    RU=8.31451D7
    TEMP = RWORK(NNcTT + (K-1)*3 + 1)
    
    L=2
    NA1 = NNcAA + (L-1)*NNCP2 + (K-1)*NNCP2T
    WRITE(21,10045)TEMP,K,-RWORK(NA1),-RWORK(NA1+1)/2.D0,-RWORK(NA1+2)/6.D0,&
    -RWORK(NA1+3)/12,-RWORK(NA1+4)/20,RWORK(NA1+5),&
    RWORK(NA1)-RWORK(NA1+6) ! (a1-a7)
    
    L=1
    NA1 = NNcAA + (L-1)*NNCP2 + (K-1)*NNCP2T
    WRITE(21,10046)K,-RWORK(NA1),-RWORK(NA1+1)/2.D0,-RWORK(NA1+2)/6.D0,&
    -RWORK(NA1+3)/12,-RWORK(NA1+4)/20,RWORK(NA1+5),&
    RWORK(NA1)-RWORK(NA1+6)
    WRITE(21,10094)
ENDDO
write(21,101)
! SUBROUTINE SKG0(T, G0)


! SUBROUTINE SKWT (WT)
write(21,110)
DO K=1,KK
    write(21,10000)ksym(k)
    write(21,10048)k,RWORK(K+MM)
ENDDO
write(21,101)
! SUBROUTINE SKWT (WT)


! SUBROUTINE AJ1(T, C, A, KK, G0, DG, DM, CTOT, P, PL)
write(21,111)
write(21,10104)KK

!LLIcRV = 1 ! starts an array of reaction indices for those with explicit reverse Arrhenius coefficients
!LIcRV = 1  ! starts an array of reaction indices for those with explicit reverse Arrhenius coefficients
!LLLIcRV = 1  ! starts an array of reaction indices for those with explicit reverse Arrhenius coefficients
!LT=1 ! Landau-Teller reactions index

DO i=1,II
    call CKSYMR (i, LOUT, IWORK, RWORK, CWORK, LT, ISTR,KERR)
    write(21,10017)i,ISTR
    CALL CKINU (i, NDIM, IWORK, RWORK, NSPEC, KI, NU)
    
    NReacSpec = 0
    NProdSpec = 0
    DO K=1,NSPEC
        IF(NU(k)<0)THEN
            NReacSpec = NReacSpec+1
            KNReacSpec(NReacSpec) = KI(k) ! reactant indice in reaction i
            KNUReacSpec(KNReacSpec(NReacSpec)) = NU(k) ! stoichiometric coefficient of reactant in reaction i
        ENDIF
    ENDDO
    DO K=1,NSPEC
        IF(NU(k)>0)THEN
            NProdSpec = NProdSpec+1
            KNProdSpec(NProdSpec) = KI(k)! product indice in reaction i
            KNUProdSpec(KNProdSpec(NProdSpec)) = NU(k)! stoichiometric coefficient of product in reaction i 
        ENDIF
    ENDDO
! KI(N) is the index of the Nth species in reaction I
! NU(N) is the stoichiometric coefficient of the Nth species in reaction I
! IREV(i) > 0 -- reversible reaction
    IF (IREV(i) > 0 .AND. IcRVIcRV(i) == 0) THEN ! if reversible but without explicit reverse Arrhenius coefficients
        write(21,10050,advance='no') ! G0_SUM
        DO k=1,IREV(i)
            write(21,10051,advance='no')NU(k),KI(k)
        ENDDO
        write(21,*)
        write(21,10052,advance='no') ! DGDT
        NUSUM = 0
        DO k=1,IREV(i)
            write(21,10053,advance='no')NU(k),KI(k) ! DG()
            NUSUM = NUSUM + NU(K)
        ENDDO
        write(21,*)
        write(21,10055)NUSUM,NUSUM ! EQINV, DGDT
    ENDIF
    
    NRFArru = NNcCO + (i-1)*(NNPAR+1)! NNPAR = 3
    IF(IcLTIcLT(i) == 0)THEN
        if(RWORK(NRFArru) > SMALL)then
            WRITE(21,10054)log(RWORK(NRFArru)),RWORK(NRFArru+1),-RWORK(NRFArru+2)! RWORK(NRFArru+2) = Ea/1.987 (unit K) ! RF
            write(21,10056)RWORK(NRFArru+1),RWORK(NRFArru+2) ! RFLGDT 
        else
            WRITE(21,11133)
            write(21,10134)
        endif
    ENDIF
!    write(21,10056)FPAR(2,i),FPAR(3,i) ! RFLGDT 
    
    
    if (ITHB(i) >= 0) then! if three body reaction 
        write(21,10059,advance='no')
        KTHBinI = 0
        DO K=1,KK
            IF (AKI(k,i) .NE. 1.D0)THEN
                KTHBinI = KTHBinI + 1
                IF (MOD(KTHBinI,5) == 0)THEN
                    write(21,10127)
                ENDIF
                write(21,10060,advance='no')AKI(k,i)-1.D0,K ! CM=CTOT+[M]
            ENDIF
        ENDDO
        WRITE(21,*)
        if ( IFOP(i) > 0) then ! if pressure depend reaction
            IF (IREV(i) > 0)THEN ! if reversible
                write(21,10058)log(FPAR(1,i)),FPAR(2,i),-FPAR(3,i) ! RL
                write(21,10061)! RB = RF/EQ; PR = RL*CM/RF
                write(21,10062)FPAR(2,i)-RWORK(NRFArru+1),FPAR(3,i)-RWORK(NRFArru+2) ! PRLGDT = D[log(PR)]/DT
                if (IFOP(i) == 3) then
                    write(21,10063) ! PRLGDM 'PC = PR/(1+PR)' 'PCLGDT = PRLGDT/TEMP1' 'PCLGDM = PRLGDM/TEMP1' 'PRLG = LOG(MAX(PR, SMALL))')
                    
                    IF(ABS(FPAR(7,i)) > 1.D+90)THEN
                        write(21,11064)1-FPAR(4,i),FPAR(5,i),FPAR(4,i),FPAR(6,i),-FPAR(7,i)/1.D10 ! TEMP1= ;TEMP2= ;TEMP3= ;
                        write(21,11165)-1.D0/FPAR(5,i),-1.D0/FPAR(6,i),FPAR(7,i)/1.D10 ! FCENT= ,FCNTDT= ,FTLGDT= ,FTLG= 
                    ELSE
                        write(21,10064)1-FPAR(4,i),FPAR(5,i),FPAR(4,i),FPAR(6,i),-FPAR(7,i) ! TEMP1= ;TEMP2= ;TEMP3= ;
                        write(21,10065)-1.D0/FPAR(5,i),-1.D0/FPAR(6,i),FPAR(7,i) ! FCENT= ,FCNTDT= ,FTLGDT= ,FTLG= 
                    ENDIF
                    
                    write(21,10066)
                elseif (IFOP(i) == 4) then
                    write(21,10063) ! PRLGDM 'PC = PR/(1+PR)' 'PCLGDT = PRLGDT/TEMP1' 'PCLGDM = PRLGDM/TEMP1' 'PRLG = LOG(MAX(PR, SMALL))')
                    IF(ABS(FPAR(7,i)) > 1.D+90)THEN
                        write(21,11064)1-FPAR(4,i),FPAR(5,i),FPAR(4,i),FPAR(6,i),-FPAR(7,i)/1.D10 ! TEMP1= ;TEMP2= ;TEMP3= ;
                        write(21,11165)-1.D0/FPAR(5,i),-1.D0/FPAR(6,i),FPAR(7,i)/1.D10 ! FCENT= ,FCNTDT= ,FTLGDT= ,FTLG= 
                    ELSE
                        write(21,10064)1-FPAR(4,i),FPAR(5,i),FPAR(4,i),FPAR(6,i),-FPAR(7,i) ! TEMP1= ;TEMP2= ;TEMP3= ;
                        write(21,10065)-1.D0/FPAR(5,i),-1.D0/FPAR(6,i),FPAR(7,i) ! FCENT= ,FCNTDT= ,FTLGDT= ,FTLG= 
                    ENDIF
                    write(21,10066)
                elseif (IFOP(i) == 2) then
                    write(21,*)'ERROR--need to supply SRI fomat output'
                elseif (IFOP(i) == 1) then
                    write(21,10126)
                endif
                
                WRITE(21,10075,advance='no') ! WF = RF * C(X)**',I4
                DO K=1,NSPEC
                    IF(NU(k)<0)THEN
                        WRITE(21,10076,advance='no')KI(k),ABS(NU(k)) ! *C(X)**
                    ENDIF
                ENDDO
                WRITE(21,*)
                WRITE(21,10077) ! WFDT ; WFDM
                WRITE(21,10078,advance='no') ! WB = RB
                DO K=1,NSPEC
                    IF(NU(k)>0)THEN
                        WRITE(21,10076,advance='no')KI(k),ABS(NU(k))
                    ENDIF
                ENDDO
                WRITE(21,*)
                WRITE(21,10079) ! WBDT ; WBDM ; DWDT

                ! DTDT(K) = DTDT(K) +/-DWDT
                DO K=1,NSPEC
                   WRITE(21,10105)KI(k),KI(k),NU(k)
                   WRITE(21,10110)KI(k),KI(k),NU(k)
                ENDDO
                
                ! A(K, KK+1) = A(K, KK+1) +/-DWDT
                DO K=1,NSPEC 
                   WRITE(21,10084)KI(k)+1,1,KI(k)+1,1,NU(k)
                ENDDO
                
    ! KNReacSpec(NReacSpec) =  reactant indice in reaction i
    ! KNUReacSpec(KNReacSpec(NReacSpec)) = stoichiometric coefficient of reactant in reaction i  
                DO K=1,NReacSpec
                    WRITE(21,10085,advance='no')! WFDC
                    WRITE(21,10087,advance='no')ABS(KNUReacSpec(KNReacSpec(k))),& 
                                            KNReacSpec(K),ABS(KNUReacSpec(KNReacSpec(k)))-1 ! I4,'*C(',I4,')**',I4
                    DO K2=1,NReacSpec
                        IF(K2 .NE. K)THEN
                            WRITE(21,10076,advance='no')KNReacSpec(K2),ABS(KNUReacSpec(KNReacSpec(k2)))
                        ENDIF
                    ENDDO
                    WRITE(21,*)
                    
                    ! DTDC(K) = DTDC(K) +/-WFDC
                    WRITE(21,10113)KNReacSpec(K),KNReacSpec(K),KNUReacSpec(KNReacSpec(k))
                    
                    ! A(K3, K) = A(K3, K) +/-WFDC
                    DO K3=1,NSPEC 
                        WRITE(21,10089)KI(k3)+1,KNReacSpec(k)+1,&
                                    KI(k3)+1,KNReacSpec(k)+1,NUKI(KI(k3),I)
                    ENDDO
                ENDDO  

                DO K=1,NProdSpec
                    WRITE(21,10086,advance='no')! WBDC
                    WRITE(21,10087,advance='no')ABS(KNUProdSpec(KNProdSpec(k))),&
                                            KNProdSpec(K),ABS(KNUProdSpec(KNProdSpec(k)))-1
                    DO K2=1,NProdSpec
                        IF(K2 .NE. K)THEN
                            WRITE(21,10076,advance='no')KNProdSpec(K2),ABS(KNUProdSpec(KNProdSpec(k2)))
                        ENDIF
                    ENDDO
                    WRITE(21,*)

                    ! DTDC(K) = DTDC(K) +/-WBDC
                    WRITE(21,10114)KNProdSpec(K),KNProdSpec(K),KNUProdSpec(KNProdSpec(k))
                    
                    ! A(K3, K) = A(K3, K) +/-WBDC
                    DO K3=1,NSPEC 
                        WRITE(21,10090)KI(k3)+1,KNProdSpec(K)+1,&
                                    KI(k3)+1,KNProdSpec(K)+1,-NUKI(KI(k3),I)
                    ENDDO
                ENDDO
                
                WRITE(21,10091) ! DWDM
                DO k1=1,NSPEC
                    DO K=1,KK
                        IF (AKI(k,i) .NE. 1.D0)THEN
                            write(21,10092)KI(k1)+1,K+1,KI(k1)+1,K+1,AKI(k,i)-1.D0,NU(k1)
                        ENDIF
                    ENDDO
                    write(21,10093)KI(k1),KI(k1),NU(k1)
                ENDDO
    !            DO K=1,KK ! A(KK+1,K)
    !                IF (AKI(k,i) .NE. 1.D0)THEN
    !                    write(21,10121)1,K+1,1,K+1,AKI(k,i)-1.D0,k,k
    !                ENDIF
    !            ENDDO 
            ELSEIF(IREV(i) < 0)THEN ! if irreversible
                write(21,10058)log(FPAR(1,i)),FPAR(2,i),-FPAR(3,i) ! RL
                WRITE(21,10069)! RB = 0
                write(21,10070)! RBLGDT = 0
                write(21,11070)! CM = MAX(CM, SMALL), PR = RL*CM/RF'
                write(21,10062)FPAR(2,i)-RWORK(NRFArru+1),FPAR(3,i)-RWORK(NRFArru+2) ! PRLGDT = D[log(PR)]/DT
                if (IFOP(i) == 3) then
                    write(21,10063) ! PRLGDM 'PC = PR/(1+PR)' 'PCLGDT = PRLGDT/TEMP1' 'PCLGDM = PRLGDM/TEMP1' 'PRLG = LOG(MAX(PR, SMALL))')
                    IF(ABS(FPAR(7,i)) > 1.D+90)THEN
                        write(21,11064)1-FPAR(4,i),FPAR(5,i),FPAR(4,i),FPAR(6,i),-FPAR(7,i)/1.D10 ! TEMP1= ;TEMP2= ;TEMP3= ;
                        write(21,11165)-1.D0/FPAR(5,i),-1.D0/FPAR(6,i),FPAR(7,i)/1.D10 ! FCENT= ,FCNTDT= ,FTLGDT= ,FTLG= 
                    ELSE
                        write(21,10064)1-FPAR(4,i),FPAR(5,i),FPAR(4,i),FPAR(6,i),-FPAR(7,i) ! TEMP1= ;TEMP2= ;TEMP3= ;
                        write(21,10065)-1.D0/FPAR(5,i),-1.D0/FPAR(6,i),FPAR(7,i) ! FCENT= ,FCNTDT= ,FTLGDT= ,FTLG= 
                    ENDIF
                    write(21,10066)
                elseif (IFOP(i) == 4) then
                    write(21,10063) ! PRLGDM 'PC = PR/(1+PR)' 'PCLGDT = PRLGDT/TEMP1' 'PCLGDM = PRLGDM/TEMP1' 'PRLG = LOG(MAX(PR, SMALL))')
                    IF(ABS(FPAR(7,i)) > 1.D+90)THEN
                        write(21,11064)1-FPAR(4,i),FPAR(5,i),FPAR(4,i),FPAR(6,i),-FPAR(7,i)/1.D10 ! TEMP1= ;TEMP2= ;TEMP3= ;
                        write(21,11165)-1.D0/FPAR(5,i),-1.D0/FPAR(6,i),FPAR(7,i)/1.D10 ! FCENT= ,FCNTDT= ,FTLGDT= ,FTLG= 
                    ELSE
                        write(21,10064)1-FPAR(4,i),FPAR(5,i),FPAR(4,i),FPAR(6,i),-FPAR(7,i) ! TEMP1= ;TEMP2= ;TEMP3= ;
                        write(21,10065)-1.D0/FPAR(5,i),-1.D0/FPAR(6,i),FPAR(7,i) ! FCENT= ,FCNTDT= ,FTLGDT= ,FTLG= 
                    ENDIF
                    write(21,10066)
                elseif (IFOP(i) == 2) then
                    write(21,*)'ERROR--need to supply SRI fomat output'
                elseif (IFOP(i) == 1) then
                    write(21,10126)
                endif
                
                WRITE(21,10075,advance='no') ! WF = RF * C(X)**',I4
                DO K=1,NSPEC
                    IF(NU(k)<0)THEN
                        WRITE(21,10076,advance='no')KI(k),ABS(NU(k)) ! *C(X)**
                    ENDIF
                ENDDO
                WRITE(21,*)
                WRITE(21,10077) ! WFDT ; WFDM
                WRITE(21,10078,advance='no') ! WB = RB
                DO K=1,NSPEC
                    IF(NU(k)>0)THEN
                        WRITE(21,10076,advance='no')KI(k),ABS(NU(k))
                    ENDIF
                ENDDO
                WRITE(21,*)
                WRITE(21,10079) ! WBDT ; WBDM ; DWDT

                ! DTDT(K) = DTDT(K) +/-DWDT
                DO K=1,NSPEC
                   WRITE(21,10105)KI(k),KI(k),NU(k)
                   WRITE(21,10110)KI(k),KI(k),NU(k)
                ENDDO
                
                ! A(K, KK+1) = A(K, KK+1) +/-DWDT
                DO K=1,NSPEC 
                   WRITE(21,10084)KI(k)+1,1,KI(k)+1,1,NU(k)
                ENDDO
                
    ! KNReacSpec(NReacSpec) =  reactant indice in reaction i
    ! KNUReacSpec(KNReacSpec(NReacSpec)) = stoichiometric coefficient of reactant in reaction i  
                DO K=1,NReacSpec
                    WRITE(21,10085,advance='no')! WFDC
                    WRITE(21,10087,advance='no')ABS(KNUReacSpec(KNReacSpec(k))),& 
                                            KNReacSpec(K),ABS(KNUReacSpec(KNReacSpec(k)))-1 ! I4,'*C(',I4,')**',I4
                    DO K2=1,NReacSpec
                        IF(K2 .NE. K)THEN
                            WRITE(21,10076,advance='no')KNReacSpec(K2),ABS(KNUReacSpec(KNReacSpec(k2)))
                        ENDIF
                    ENDDO
                    WRITE(21,*)
                    
                    ! DTDC(K) = DTDC(K) +/-WFDC
                    WRITE(21,10113)KNReacSpec(K),KNReacSpec(K),KNUReacSpec(KNReacSpec(k))
                    
                    ! A(K3, K) = A(K3, K) +/-WFDC
                    DO K3=1,NSPEC 
                        WRITE(21,10089)KI(k3)+1,KNReacSpec(k)+1,&
                                    KI(k3)+1,KNReacSpec(k)+1,NUKI(KI(k3),I)
                    ENDDO
                ENDDO  

                DO K=1,NProdSpec
                    WRITE(21,10086,advance='no')! WBDC
                    WRITE(21,10087,advance='no')ABS(KNUProdSpec(KNProdSpec(k))),&
                                            KNProdSpec(K),ABS(KNUProdSpec(KNProdSpec(k)))-1
                    DO K2=1,NProdSpec
                        IF(K2 .NE. K)THEN
                            WRITE(21,10076,advance='no')KNProdSpec(K2),ABS(KNUProdSpec(KNProdSpec(k2)))
                        ENDIF
                    ENDDO
                    WRITE(21,*)

                    ! DTDC(K) = DTDC(K) +/-WBDC
                    WRITE(21,10114)KNProdSpec(K),KNProdSpec(K),KNUProdSpec(KNProdSpec(k))
                    
                    ! A(K3, K) = A(K3, K) +/-WBDC
                    DO K3=1,NSPEC 
                        WRITE(21,10090)KI(k3)+1,KNProdSpec(K)+1,&
                                    KI(k3)+1,KNProdSpec(K)+1,-NUKI(KI(k3),I)
                    ENDDO
                ENDDO
                
                WRITE(21,10091) ! DWDM
                DO k1=1,NSPEC
                    DO K=1,KK
                        IF (AKI(k,i) .NE. 1.D0)THEN
                            write(21,10092)KI(k1)+1,K+1,KI(k1)+1,K+1,AKI(k,i)-1.D0,NU(k1)
                        ENDIF
                    ENDDO
                    write(21,10093)KI(k1),KI(k1),NU(k1)
                ENDDO
            ENDIF
        else ! if three body reaction but not pressure depend reaction
            IF (IREV(i) > 0) THEN  ! if reverse reaction
!                IF (i == IWORK(IIcRV+LIcRV-1) .AND. LIcRV <= NNREV) THEN ! if reaction i with explicit reverse Arrhenius coefficients
                if (IcLTIcLT(i) == 0)then ! if not Landau-Teller reaction
                    IF (IcRVIcRV(i)>0) THEN ! if reaction i with explicit reverse Arrhenius coefficients
                        NRBArru = NNcRV + (IcRVNRV(i)-1)*(NNPAR+1) ! NNPAR = 3
                        if(RWORK(NRBArru)>SMALL)then
                            WRITE(21,10067)log(RWORK(NRBArru)),RWORK(NRBArru+1),-RWORK(NRBArru+2)! RB
                            write(21,10068)RWORK(NRBArru+1),RWORK(NRBArru+2)! RBLGDT
                        else
                            write(21,10069)
                            write(21,10070)
                        endif
                        WRITE(21,10136)! PC = CM, PCLGDM = 1.D0/CM, RF = RF*PC, RB = RB*PC
                        
                        WRITE(21,10075,advance='no')
                        DO K=1,NSPEC
                            IF(NU(k)<0)THEN
                                WRITE(21,10076,advance='no')KI(k),ABS(NU(k))
                            ENDIF
                        ENDDO
                        WRITE(21,*)
                        WRITE(21,10080)
                        WRITE(21,10078,advance='no')
                        DO K=1,NSPEC
                            IF(NU(k)>0)THEN
                                WRITE(21,10076,advance='no')KI(k),ABS(NU(k))
                            ENDIF
                        ENDDO
                        WRITE(21,*)
                        WRITE(21,10081)! WBDT ; WBDM ; DWDT
                        
                        ! DTDT(K) = DTDT(K) +/-DWDT
                        DO K=1,NSPEC
                           WRITE(21,10105)KI(k),KI(k),NU(k)
                           WRITE(21,10110)KI(k),KI(k),NU(k)
                        ENDDO
                        
                        ! A(K, KK+1) = A(K, KK+1) +/-DWDT
                        DO K=1,NSPEC 
                           WRITE(21,10084)KI(k)+1,1,KI(k)+1,1,NU(k)
                        ENDDO
                        
                        DO K=1,NReacSpec
                            WRITE(21,10085,advance='no')! WFDC
                            WRITE(21,10087,advance='no')ABS(KNUReacSpec(KNReacSpec(k))),&
                                                    KNReacSpec(K),ABS(KNUReacSpec(KNReacSpec(k)))-1
                            DO K2=1,NReacSpec
                                IF(K2 .NE. K)THEN
                                    WRITE(21,10076,advance='no')KNReacSpec(K2),ABS(KNUReacSpec(KNReacSpec(k2)))
                                ENDIF
                            ENDDO
                            WRITE(21,*)
                            
                            ! DTDC(K) = DTDC(K) +/-WFDC
                            WRITE(21,10113)KNReacSpec(K),KNReacSpec(K),KNUReacSpec(KNReacSpec(k))
                            
                            ! A(K3, K) = A(K3, K) +/-WFDC
                            DO K3=1,NSPEC 
                                WRITE(21,10089)KI(k3)+1,KNReacSpec(K)+1,&
                                    KI(k3)+1,KNReacSpec(K)+1,NUKI(KI(k3),I)
                            ENDDO
                        ENDDO  
                        
                        DO K=1,NProdSpec
                            WRITE(21,10086,advance='no')! WBDC
                            WRITE(21,10087,advance='no')ABS(KNUProdSpec(KNProdSpec(k))),&
                                                    KNProdSpec(K),ABS(KNUProdSpec(KNProdSpec(k)))-1
                            DO K2=1,NProdSpec
                                IF(K2 .NE. K)THEN
                                    WRITE(21,10076,advance='no')KNProdSpec(K2),KNUProdSpec(KNProdSpec(k2))
                                ENDIF
                            ENDDO
                            WRITE(21,*)

                            ! DTDC(K) = DTDC(K) +/-WBDC
                            WRITE(21,10114)KNProdSpec(K),KNProdSpec(K),KNUProdSpec(KNProdSpec(k))
                            
                            ! A(K3, K) = A(K3, K) +/-WBDC
                            DO K3=1,NSPEC 
                                WRITE(21,10090)KI(k3)+1,KNProdSpec(K)+1,&
                                    KI(k3)+1,KNProdSpec(K)+1,-NUKI(KI(k3),I)
                            ENDDO
                        ENDDO
                        
                        WRITE(21,10091) ! DWDM
                        DO k1=1,NSPEC
                            DO K=1,KK
                                IF (AKI(k,i) .NE. 1.D0)THEN
                                    write(21,10092)KI(k1)+1,K+1,KI(k1)+1,K+1,AKI(k,i)-1.D0,NU(k1)
                                ENDIF
                            ENDDO
                            write(21,10093)KI(k1),KI(k1),NU(k1)
                        ENDDO
    !                    DO K=1,KK ! A(KK+1,K)
    !                        IF (AKI(k,i) .NE. 1.D0)THEN
    !                            write(21,10121)1,K+1,1,K+1,AKI(k,i)-1.D0,k,k
    !                        ENDIF
    !                    ENDDO
                        
    !                    WRITE(21,10091) ! DWDM
    !                    DO k1=1,NSPEC
    !                        DO K=1,KK
    !                            IF (AKI(k,i) .NE. 1.D0)THEN
    !                                write(21,10092)KI(k1),K,KI(k1),K,AKI(k,i)-1.D0
    !                            ENDIF
    !                        ENDDO
    !                        write(21,10093)KI(k1),KI(k1),NU(k1)
    !                    ENDDO
                        
                    ELSE ! if reaction i without explicit reverse Arrhenius coefficients
                        WRITE(21,10071)
                        WRITE(21,10075,advance='no')
                        DO K=1,NSPEC
                            IF(NU(k)<0)THEN
                                WRITE(21,10076,advance='no')KI(k),ABS(NU(k))
                            ENDIF
                        ENDDO
                        WRITE(21,*)
                        WRITE(21,10080)
                        WRITE(21,10078,advance='no')
                        DO K=1,NSPEC
                            IF(NU(k)>0)THEN
                                WRITE(21,10076,advance='no')KI(k),ABS(NU(k))
                            ENDIF
                        ENDDO
                        WRITE(21,*)
                        WRITE(21,10081)! WBDT ; WBDM ; DWDT
                        
                        ! DTDT(K) = DTDT(K) +/-DWDT
                        DO K=1,NSPEC
                           WRITE(21,10105)KI(k),KI(k),NU(k)
                           WRITE(21,10110)KI(k),KI(k),NU(k)
                        ENDDO
                        
                        ! A(K, KK+1) = A(K, KK+1) +/-DWDT
                        DO K=1,NSPEC 
                           WRITE(21,10084)KI(k)+1,1,KI(k)+1,1,NU(k)
                        ENDDO
                        
                        DO K=1,NReacSpec
                            WRITE(21,10085,advance='no')! WFDC
                            WRITE(21,10087,advance='no')ABS(KNUReacSpec(KNReacSpec(k))),&
                                                    KNReacSpec(K),ABS(KNUReacSpec(KNReacSpec(k)))-1
                            DO K2=1,NReacSpec
                                IF(K2 .NE. K)THEN
                                    WRITE(21,10076,advance='no')KNReacSpec(K2),ABS(KNUReacSpec(KNReacSpec(k2)))
                                ENDIF
                            ENDDO
                            WRITE(21,*)
                            
                            ! DTDC(K) = DTDC(K) +/-WFDC
                            WRITE(21,10113)KNReacSpec(K),KNReacSpec(K),KNUReacSpec(KNReacSpec(k))
                            
                            ! A(K3, K) = A(K3, K) +/-WFDC
                            DO K3=1,NSPEC 
                                WRITE(21,10089)KI(k3)+1,KNReacSpec(K)+1,&
                                    KI(k3)+1,KNReacSpec(K)+1,NUKI(KI(k3),I)
                            ENDDO
                        ENDDO  
                        
                        DO K=1,NProdSpec
                            WRITE(21,10086,advance='no')! WBDC
                            WRITE(21,10087,advance='no')ABS(KNUProdSpec(KNProdSpec(k))),&
                                                    KNProdSpec(K),ABS(KNUProdSpec(KNProdSpec(k)))-1
                            DO K2=1,NProdSpec
                                IF(K2 .NE. K)THEN
                                    WRITE(21,10076,advance='no')KNProdSpec(K2),KNUProdSpec(KNProdSpec(k2))
                                ENDIF
                            ENDDO
                            WRITE(21,*)
                            
                            ! DTDC(K) = DTDC(K) +/-WBDC
                            WRITE(21,10114)KNProdSpec(K),KNProdSpec(K),KNUProdSpec(KNProdSpec(k))
                            
                            ! A(K3, K) = A(K3, K) +/-WBDC
                            DO K3=1,NSPEC 
                                WRITE(21,10090)KI(k3)+1,KNProdSpec(K)+1,&
                                    KI(k3)+1,KNProdSpec(K)+1,-NUKI(KI(k3),I)
                            ENDDO
                        ENDDO
                        
                        WRITE(21,10091) ! DWDM
                        DO k1=1,NSPEC
                            DO K=1,KK
                                IF (AKI(k,i) .NE. 1.D0)THEN
                                    write(21,10092)KI(k1)+1,K+1,KI(k1)+1,K+1,AKI(k,i)-1.D0,NU(k1)
                                ENDIF
                            ENDDO
                            write(21,10093)KI(k1),KI(k1),NU(k1)
                        ENDDO
    !                    DO K=1,KK ! A(KK+1,K)
    !                        IF (AKI(k,i) .NE. 1.D0)THEN
    !                            write(21,10121)1,K+1,1,K+1,AKI(k,i)-1.D0,k,k
    !                        ENDIF
    !                    ENDDO
                        
                    ENDIF
                elseif (IcLTIcLT(i) > 0)then ! if Landau-Teller
                    write(*,*)'ERROR--need to supply Landau-Teller and pressure depend and reverse reaction'
                endif
            ELSE ! if not reverse reaction
                WRITE(21,10069)! RB = 0
                write(21,10070)! RBLGDT = 0
                WRITE(21,10072)
                
                WRITE(21,10075,advance='no')
                DO K=1,NSPEC
                    IF(NU(k)<0)THEN
                        WRITE(21,10076,advance='no')KI(k),ABS(NU(k))
                    ENDIF
                ENDDO
                WRITE(21,*)
                WRITE(21,10080)
                WRITE(21,10078,advance='no')
                DO K=1,NSPEC
                    IF(NU(k)>0)THEN
                        WRITE(21,10076,advance='no')KI(k),ABS(NU(k))
                    ENDIF
                ENDDO
                WRITE(21,*)
                WRITE(21,10081)! WBDT ; WBDM ; DWDT
                
                ! DTDT(K) = DTDT(K) +/-DWDT
                DO K=1,NSPEC
                   WRITE(21,10105)KI(k),KI(k),NU(k)
                   WRITE(21,10110)KI(k),KI(k),NU(k)
                ENDDO
                
                ! A(K, KK+1) = A(K, KK+1) +/-DWDT
                DO K=1,NSPEC 
                   WRITE(21,10084)KI(k)+1,1,KI(k)+1,1,NU(k)
                ENDDO
                
                DO K=1,NReacSpec
                    WRITE(21,10085,advance='no')! WFDC
                    WRITE(21,10087,advance='no')ABS(KNUReacSpec(KNReacSpec(k))),&
                                            KNReacSpec(K),ABS(KNUReacSpec(KNReacSpec(k)))-1
                    DO K2=1,NReacSpec
                        IF(K2 .NE. K)THEN
                            WRITE(21,10076,advance='no')KNReacSpec(K2),ABS(KNUReacSpec(KNReacSpec(k2)))
                        ENDIF
                    ENDDO
                    WRITE(21,*)

                    ! DTDC(K) = DTDC(K) +/-WFDC
                    WRITE(21,10113)KNReacSpec(K),KNReacSpec(K),KNUReacSpec(KNReacSpec(k))
                    
                    ! A(K3, K) = A(K3, K) +/-WFDC
                    DO K3=1,NSPEC 
                        WRITE(21,10089)KI(k3)+1,KNReacSpec(K)+1,&
                                KI(k3)+1,KNReacSpec(K)+1,NUKI(KI(k3),I)
                    ENDDO
                ENDDO  
                
                DO K=1,NProdSpec
                    WRITE(21,10086,advance='no')! WBDC
                    WRITE(21,10087,advance='no')ABS(KNUProdSpec(KNProdSpec(k))),&
                                            KNProdSpec(K),ABS(KNUProdSpec(KNProdSpec(k)))-1
                    DO K2=1,NProdSpec
                        IF(K2 .NE. K)THEN
                            WRITE(21,10076,advance='no')KNProdSpec(K2),KNUProdSpec(KNProdSpec(k2))
                        ENDIF
                    ENDDO
                    WRITE(21,*)

                    ! DTDC(K) = DTDC(K) +/-WBDC
                    WRITE(21,10114)KNProdSpec(K),KNProdSpec(K),KNUProdSpec(KNProdSpec(k))
                    
                    ! A(K3, K) = A(K3, K) +/-WBDC
                    DO K3=1,NSPEC 
                        WRITE(21,10090)KI(k3)+1,KNProdSpec(K)+1,&
                                KI(k3)+1,KNProdSpec(K)+1,-NUKI(KI(k3),I)
                    ENDDO
                ENDDO
                
                WRITE(21,10091) ! DWDM
                DO k1=1,NSPEC
                    DO K=1,KK
                        IF (AKI(k,i) .NE. 1.D0)THEN
                            write(21,10092)KI(k1)+1,K+1,KI(k1)+1,K+1,AKI(k,i)-1.D0,NU(k1)
                        ENDIF
                    ENDDO
                    write(21,10093)KI(k1),KI(k1),NU(k1)
                ENDDO
!                DO K=1,KK ! A(KK+1,K)
!                    IF (AKI(k,i) .NE. 1.D0)THEN
!                        write(21,10121)1,K+1,1,K+1,AKI(k,i)-1.D0,k,k
!                    ENDIF
!                ENDDO
                
            ENDIF
        endif
    else ! if not three body reaction 
!        IF (i == IWORK(IIcLT+LT-1) .AND. LT <= NNLAN)THEN ! IF not three body reactions AND Landau-Teller reactions
        IF (IcLTIcLT(i)>0)THEN ! IF not three body reactions AND Landau-Teller reactions
            NRFLT = NNcLT + (IcLTNLT(i)-1)*NNLAR ! NNLAR = 4
            if(RWORK(NRFArru) > SMALL)then
                WRITE(21,10074)log(RWORK(NRFArru)),RWORK(NRFArru+1),-RWORK(NRFArru+2),& ! RF
                            RWORK(NRFLT),RWORK(NRFLT+1)
                write(21,10057)RWORK(NRFArru+1),RWORK(NRFArru+2),&
                            -RWORK(NRFLT)*(1.D0/3.D0),-RWORK(NRFLT+1)*(2.D0/3.D0)! RFLGDT
            else
                write(21,11133)
                write(21,10134)
            endif
            
            IF (IREV(i) > 0) THEN  ! if reverse reaction
                WRITE(21,10073) ! RB = RF*EQINV; RBLGDT = RFLGDT +DGDT
                
                WRITE(21,10075,advance='no')
                DO K=1,NSPEC
                    IF(NU(k)<0)THEN
                        WRITE(21,10076,advance='no')KI(k),ABS(NU(k))
                    ENDIF
                ENDDO
                WRITE(21,*)
                WRITE(21,10082)
                WRITE(21,10078,advance='no')
                DO K=1,NSPEC
                    IF(NU(k)>0)THEN
                        WRITE(21,10076,advance='no')KI(k),ABS(NU(k))
                    ENDIF
                ENDDO
                WRITE(21,*)
                WRITE(21,10083) ! WBDT ; DWDT
                
                ! DTDT(K) = DTDT(K) +/-DWDT
                DO K=1,NSPEC
                   WRITE(21,10105)KI(k),KI(k),NU(k)
                   WRITE(21,10110)KI(k),KI(k),NU(k)
                ENDDO
                
                ! A(K, KK+1) = A(K, KK+1) +/-DWDT
                DO K=1,NSPEC 
                   WRITE(21,10084)KI(k)+1,1,KI(k)+1,1,NU(k)
                ENDDO
                
                DO K=1,NReacSpec
                    WRITE(21,10085,advance='no')! WFDC
                    WRITE(21,10087,advance='no')ABS(KNUReacSpec(KNReacSpec(k))),&
                                            KNReacSpec(K),ABS(KNUReacSpec(KNReacSpec(k)))-1
                    DO K2=1,NReacSpec
                        IF(K2 .NE. K)THEN
                            WRITE(21,10076,advance='no')KNReacSpec(K2),ABS(KNUReacSpec(KNReacSpec(k2)))
                        ENDIF
                    ENDDO
                    WRITE(21,*)
                    
                    ! DTDC(K) = DTDC(K) +/-WFDC
                    WRITE(21,10113)KNReacSpec(K),KNReacSpec(K),KNUReacSpec(KNReacSpec(k))

                    ! A(K3, K) = A(K3, K) +/-WFDC
                    DO K3=1,NSPEC 
                        WRITE(21,10089)KI(k3)+1,KNReacSpec(K)+1,&
                                KI(k3)+1,KNReacSpec(K)+1,NUKI(KI(k3),I)
                    ENDDO
                ENDDO  
                
                DO K=1,NProdSpec
                    WRITE(21,10086,advance='no')! WBDC
                    WRITE(21,10087,advance='no')ABS(KNUProdSpec(KNProdSpec(k))),&
                                            KNProdSpec(K),ABS(KNUProdSpec(KNProdSpec(k)))-1
                    DO K2=1,NProdSpec
                        IF(K2 .NE. K)THEN
                            WRITE(21,10076,advance='no')KNProdSpec(K2),KNUProdSpec(KNProdSpec(k2))
                        ENDIF
                    ENDDO
                    WRITE(21,*)

                    ! DTDC(K) = DTDC(K) +/-WBDC
                    WRITE(21,10114)KNProdSpec(K),KNProdSpec(K),KNUProdSpec(KNProdSpec(k))
                    
                    ! A(K3, K) = A(K3, K) +/-WBDC
                    DO K3=1,NSPEC 
                        WRITE(21,10090)KI(k3)+1,KNProdSpec(K)+1,&
                                KI(k3)+1,KNProdSpec(K)+1,-NUKI(KI(k3),I)
                    ENDDO
                ENDDO
                
            ELSE! if Landau-Teller and irreverse reaction
                WRITE(21,10069)! RB = 0
                write(21,10070)! RBLGDT = 0
                
                WRITE(21,10075,advance='no')
                DO K=1,NSPEC
                    IF(NU(k)<0)THEN
                        WRITE(21,10076,advance='no')KI(k),ABS(NU(k))
                    ENDIF
                ENDDO
                WRITE(21,*)
                WRITE(21,10082)
                WRITE(21,10078,advance='no')
                DO K=1,NSPEC
                    IF(NU(k)>0)THEN
                        WRITE(21,10076,advance='no')KI(k),ABS(NU(k))
                    ENDIF
                ENDDO
                WRITE(21,*)
                WRITE(21,10083) ! WBDT ; DWDT
                
                ! DTDT(K) = DTDT(K) +/-DWDT
                DO K=1,NSPEC
                   WRITE(21,10105)KI(k),KI(k),NU(k)
                   WRITE(21,10110)KI(k),KI(k),NU(k)
                ENDDO
                
                ! A(K, KK+1) = A(K, KK+1) +/-DWDT
                DO K=1,NSPEC 
                   WRITE(21,10084)KI(k)+1,1,KI(k)+1,1,NU(k)
                ENDDO
                
                DO K=1,NReacSpec
                    WRITE(21,10085,advance='no')! WFDC
                    WRITE(21,10087,advance='no')ABS(KNUReacSpec(KNReacSpec(k))),&
                                            KNReacSpec(K),ABS(KNUReacSpec(KNReacSpec(k)))-1
                    DO K2=1,NReacSpec
                        IF(K2 .NE. K)THEN
                            WRITE(21,10076,advance='no')KNReacSpec(K2),ABS(KNUReacSpec(KNReacSpec(k2)))
                        ENDIF
                    ENDDO
                    WRITE(21,*)
                    
                    ! DTDC(K) = DTDC(K) +/-WFDC
                    WRITE(21,10113)KNReacSpec(K),KNReacSpec(K),KNUReacSpec(KNReacSpec(k))
                    
                    ! A(K3, K) = A(K3, K) +/-WFDC
                    DO K3=1,NSPEC 
                        WRITE(21,10089)KI(k3)+1,KNReacSpec(K)+1,&
                                KI(k3)+1,KNReacSpec(K)+1,NUKI(KI(k3),I)
                    ENDDO
                ENDDO  
                
                DO K=1,NProdSpec
                    WRITE(21,10086,advance='no')! WBDC
                    WRITE(21,10087,advance='no')ABS(KNUProdSpec(KNProdSpec(k))),&
                                            KNProdSpec(K),ABS(KNUProdSpec(KNProdSpec(k)))-1
                    DO K2=1,NProdSpec
                        IF(K2 .NE. K)THEN
                            WRITE(21,10076,advance='no')KNProdSpec(K2),KNUProdSpec(KNProdSpec(k2))
                        ENDIF
                    ENDDO
                    WRITE(21,*)

                    ! DTDC(K) = DTDC(K) +/-WBDC
                    WRITE(21,10114)KNProdSpec(K),KNProdSpec(K),KNUProdSpec(KNProdSpec(k))

                    ! A(K3, K) = A(K3, K) +/-WBDC
                    DO K3=1,NSPEC 
                        WRITE(21,10090)KI(k3)+1,KNProdSpec(K)+1,&
                                KI(k3)+1,KNProdSpec(K)+1,-NUKI(KI(k3),I)
                    ENDDO
                ENDDO
                
            ENDIF
!        IF (i == IWORK(IIcLT+LT-1) .AND. LT <= NNLAN) THEN ! IF Landau-Teller reactions
        ELSE ! IF not Landau-Teller reactions
            IF (IREV(i) > 0) THEN  ! if reverse reaction
                IF (IcRVIcRV(i)>0) THEN ! if reaction i with explicit reverse Arrhenius coefficients
!                IF (i == IWORK(IIcRV+LLLIcRV-1) .AND. LLLIcRV <= NNREV) THEN ! if reaction i with explicit reverse Arrhenius coefficients
                    NRBArru = NNcRV + (IcRVNRV(i)-1)*(NNPAR+1) ! NNPAR = 3
                    if(RWORK(NRBArru)>SMALL)then
                        WRITE(21,10067)log(RWORK(NRBArru)),RWORK(NRBArru+1),-RWORK(NRBArru+2)! RB
                        write(21,10068)RWORK(NRBArru+1),RWORK(NRBArru+2)! RBLGDT
                    else
                        WRITE(21,10069)
                        WRITE(21,10070)
                    endif
!                    LLLIcRV = LLLIcRV + 1
                    
                    WRITE(21,10075,advance='no')
                    DO K=1,NSPEC
                        IF(NU(k)<0)THEN
                            WRITE(21,10076,advance='no')KI(k),ABS(NU(k))
                        ENDIF
                    ENDDO
                    WRITE(21,*)
                    WRITE(21,10082)
                    WRITE(21,10078,advance='no')
                    DO K=1,NSPEC
                        IF(NU(k)>0)THEN
                            WRITE(21,10076,advance='no')KI(k),ABS(NU(k))
                        ENDIF
                    ENDDO
                    WRITE(21,*)
                    WRITE(21,10083)
                    
                    ! DTDT(K) = DTDT(K) +/-DWDT
                    DO K=1,NSPEC
                       WRITE(21,10105)KI(k),KI(k),NU(k)
                       WRITE(21,10110)KI(k),KI(k),NU(k)
                    ENDDO
                    
                    ! A(K, KK+1) = A(K, KK+1) +/-DWDT
                    DO K=1,NSPEC 
                       WRITE(21,10084)KI(k)+1,1,KI(k)+1,1,NU(k)
                    ENDDO
                    
                    DO K=1,NReacSpec
                        WRITE(21,10085,advance='no')! WFDC
                        WRITE(21,10087,advance='no')ABS(KNUReacSpec(KNReacSpec(k))),&
                                                KNReacSpec(K),ABS(KNUReacSpec(KNReacSpec(k)))-1
                        DO K2=1,NReacSpec
                            IF(K2 .NE. K)THEN
                                WRITE(21,10076,advance='no')KNReacSpec(K2),ABS(KNUReacSpec(KNReacSpec(k2)))
                            ENDIF
                        ENDDO
                        WRITE(21,*)
                        
                        ! DTDC(K) = DTDC(K) +/-WFDC
                        WRITE(21,10113)KNReacSpec(K),KNReacSpec(K),KNUReacSpec(KNReacSpec(k))
                        
                        ! A(K3, K) = A(K3, K) +/-WFDC
                        DO K3=1,NSPEC 
                            WRITE(21,10089)KI(k3)+1,KNReacSpec(K)+1,&
                                KI(k3)+1,KNReacSpec(K)+1,NUKI(KI(k3),I)
                        ENDDO
                    ENDDO  
                    
                    DO K=1,NProdSpec
                        WRITE(21,10086,advance='no')! WBDC
                        WRITE(21,10087,advance='no')ABS(KNUProdSpec(KNProdSpec(k))),&
                                                KNProdSpec(K),ABS(KNUProdSpec(KNProdSpec(k)))-1
                        DO K2=1,NProdSpec
                            IF(K2 .NE. K)THEN
                                WRITE(21,10076,advance='no')KNProdSpec(K2),KNUProdSpec(KNProdSpec(k2))
                            ENDIF
                        ENDDO
                        WRITE(21,*)

                        ! DTDC(K) = DTDC(K) +/-WBDC
                        WRITE(21,10114)KNProdSpec(K),KNProdSpec(K),KNUProdSpec(KNProdSpec(k))
                        
                        ! A(K3, K) = A(K3, K) +/-WBDC
                        DO K3=1,NSPEC 
                            WRITE(21,10090)KI(k3)+1,KNProdSpec(K)+1,&
                                KI(k3)+1,KNProdSpec(K)+1,-NUKI(KI(k3),I)
                        ENDDO
                    ENDDO
                    
                ELSE ! if reaction i without explicit reverse Arrhenius coefficients
                    WRITE(21,10073)! RB = RF*EQINV; RBLGDT = RFLGDT +DGDT
                    
                    WRITE(21,10075,advance='no')
                    DO K=1,NSPEC
                        IF(NU(k)<0)THEN
                            WRITE(21,10076,advance='no')KI(k),ABS(NU(k))
                        ENDIF
                    ENDDO
                    WRITE(21,*)
                    WRITE(21,10082)
                    WRITE(21,10078,advance='no')
                    DO K=1,NSPEC
                        IF(NU(k)>0)THEN
                            WRITE(21,10076,advance='no')KI(k),ABS(NU(k))
                        ENDIF
                    ENDDO
                    WRITE(21,*)
                    WRITE(21,10083)
                    
                    ! DTDT(K) = DTDT(K) +/-DWDT
                    DO K=1,NSPEC
                       WRITE(21,10105)KI(k),KI(k),NU(k)
                       WRITE(21,10110)KI(k),KI(k),NU(k)
                    ENDDO
                    
                    ! A(K, KK+1) = A(K, KK+1) +/-DWDT
                    DO K=1,NSPEC 
                       WRITE(21,10084)KI(k)+1,1,KI(k)+1,1,NU(k)
                    ENDDO
                    
!                    write(*,*)'i = ',i
                    DO K=1,NReacSpec
                        WRITE(21,10085,advance='no')! WFDC
                        WRITE(21,10087,advance='no')ABS(KNUReacSpec(KNReacSpec(k))),&
                                                KNReacSpec(K),ABS(KNUReacSpec(KNReacSpec(k)))-1
                        DO K2=1,NReacSpec
                            IF(K2 .NE. K)THEN
                                WRITE(21,10076,advance='no')KNReacSpec(K2),ABS(KNUReacSpec(KNReacSpec(k2)))
                            ENDIF
                        ENDDO
                        WRITE(21,*)

                        ! DTDC(K) = DTDC(K) +/-WFDC
                        WRITE(21,10113)KNReacSpec(K),KNReacSpec(K),KNUReacSpec(KNReacSpec(k))

                        ! A(K3, K) = A(K3, K) +/-WFDC
                        DO K3=1,NSPEC 
                            WRITE(21,10089)KI(k3)+1,KNReacSpec(K)+1,&
                                KI(k3)+1,KNReacSpec(K)+1,NUKI(KI(k3),I)
                        ENDDO
                    ENDDO  
                    
                    DO K=1,NProdSpec
                        WRITE(21,10086,advance='no')! WBDC
                        WRITE(21,10087,advance='no')ABS(KNUProdSpec(KNProdSpec(k))),&
                                                KNProdSpec(K),ABS(KNUProdSpec(KNProdSpec(k)))-1
                        DO K2=1,NProdSpec
                            IF(K2 .NE. K)THEN
                                WRITE(21,10076,advance='no')KNProdSpec(K2),KNUProdSpec(KNProdSpec(k2))
                            ENDIF
                        ENDDO
                        WRITE(21,*)

                        ! DTDC(K) = DTDC(K) +/-WBDC
                        WRITE(21,10114)KNProdSpec(K),KNProdSpec(K),KNUProdSpec(KNProdSpec(k))
                        
                        ! A(K3, K) = A(K3, K) +/-WBDC
                        DO K3=1,NSPEC 
                            WRITE(21,10090)KI(k3)+1,KNProdSpec(K)+1,&
                                KI(k3)+1,KNProdSpec(K)+1,-NUKI(KI(k3),I)
                        ENDDO
                    ENDDO
                    
                ENDIF
            ELSE ! if not reverse reaction
                WRITE(21,10069)! RB = 0
                write(21,10070)! RBLGDT = 0
                
                WRITE(21,10075,advance='no')
                DO K=1,NSPEC
                    IF(NU(k)<0)THEN
                        WRITE(21,10076,advance='no')KI(k),ABS(NU(k))
                    ENDIF
                ENDDO
                WRITE(21,*)
                WRITE(21,10082)
                WRITE(21,10078,advance='no')
                DO K=1,NSPEC
                    IF(NU(k)>0)THEN
                        WRITE(21,10076,advance='no')KI(k),ABS(NU(k))
                    ENDIF
                ENDDO
                WRITE(21,*)
                WRITE(21,10083)
                
                ! DTDT(K) = DTDT(K) +/-DWDT
                DO K=1,NSPEC
                   WRITE(21,10105)KI(k),KI(k),NU(k)
                   WRITE(21,10110)KI(k),KI(k),NU(k)
                ENDDO
                
                ! A(K, KK+1) = A(K, KK+1) +/-DWDT
                DO K=1,NSPEC 
                   WRITE(21,10084)KI(k)+1,1,KI(k)+1,1,NU(k)
                ENDDO
                
                DO K=1,NReacSpec
                    WRITE(21,10085,advance='no')! WFDC
                    WRITE(21,10087,advance='no')ABS(KNUReacSpec(KNReacSpec(k))),&
                                            KNReacSpec(K),ABS(KNUReacSpec(KNReacSpec(k)))-1
                    DO K2=1,NReacSpec
                        IF(K2 .NE. K)THEN
                            WRITE(21,10076,advance='no')KNReacSpec(K2),ABS(KNUReacSpec(KNReacSpec(k2)))
                        ENDIF
                    ENDDO
                    WRITE(21,*)
                    
                    ! DTDC(K) = DTDC(K) +/-WFDC
                    WRITE(21,10113)KNReacSpec(K),KNReacSpec(K),KNUReacSpec(KNReacSpec(k))
                    
                    ! A(K3, K) = A(K3, K) +/-WFDC
                    DO K3=1,NSPEC 
                        WRITE(21,10089)KI(k3)+1,KNReacSpec(K)+1,&
                                KI(k3)+1,KNReacSpec(K)+1,NUKI(KI(k3),I)
                    ENDDO
                ENDDO  
                
                DO K=1,NProdSpec
                    WRITE(21,10086,advance='no')! WBDC
                    WRITE(21,10087,advance='no')ABS(KNUProdSpec(KNProdSpec(k))),&
                                            KNProdSpec(K),ABS(KNUProdSpec(KNProdSpec(k)))-1
                    DO K2=1,NProdSpec
                        IF(K2 .NE. K)THEN
                            WRITE(21,10076,advance='no')KNProdSpec(K2),KNUProdSpec(KNProdSpec(k2))
                        ENDIF
                    ENDDO
                    WRITE(21,*)

                    ! DTDC(K) = DTDC(K) +/-WBDC
                    WRITE(21,10114)KNProdSpec(K),KNProdSpec(K),KNUProdSpec(KNProdSpec(k))

                    ! A(K3, K) = A(K3, K) +/-WBDC
                    DO K3=1,NSPEC 
                        WRITE(21,10090)KI(k3)+1,KNProdSpec(K)+1,&
                                KI(k3)+1,KNProdSpec(K)+1,-NUKI(KI(k3),I)
                    ENDDO
                ENDDO
                
            ENDIF
        ENDIF
    endif
    
    



ENDDO
write(21,*)
!write(21,10108)
! A(1,1)
write(21,10106)1,1,1,1
DO K=1,KK
    write(21,10107)K,K,K,K
ENDDO
write(21,10111)

!DO K=1,KK
!WRITE(21,10115)1,K+1,1,K+1,K,K
!ENDDO

!DO K=1,KK+1
!write(21,10116)KK+1,K,KK+1,K
!ENDDO
write(21,*)
!write(21,10118)KK+1,1,1
write(21,10118)1,1,1,1
!write(21,*)
!write(21,10117)KK+1,KK
!write(21,*)
!write(21,10123)KK

write(21,101)
! SUBROUTINE AJ1(T, C, A, KK, G0, DG, DM, CTOT, P, PL)

! SUBROUTINE DMDC(A, KK, DM, HML, WT, RHO, CPB)
write(21,112)

write(21,*)
write(21,10131)KK

write(21,*)
write(21,10095)KK+1

NTHSP=0
do i=1,II
    ! IF third-body reaction
    if (ITHB(i) > 0)then
        CALL CKINU (i, NDIM, IWORK, RWORK, NSPEC, KI, NU)
! KI(N) is the index of the Nth species in reaction I
! NU(N) is the stoichiometric coefficient of the Nth species in reaction I
        DO K=1,NSPEC
            NTHSP = NTHSP+1
            NSPinTH(NTHSP)=KI(K)
        enddo
    endif
enddo
!write(*,*)'NSPinTH = ',NSPinTH
call delSameNumbANDsort(NTHSP,NSPinTH,NNTHSP,NSPinTH2)
!write(*,*)'NSPinTH2 = ',NSPinTH2
DO K=1,NNTHSP
!    write(21,10096)NSPinTH2(K)+1,NSPinTH2(K)+1,NSPinTH2(K),NSPinTH2(K)
    write(21,10096)NSPinTH2(K)+1,NSPinTH2(K)+1,NSPinTH2(K)
ENDDO
WRITE(21,10102)
WRITE(21,*)

! A(1,k+1)
WRITE(21,10129)KK,KK
write(21,*)

! DFDC ----> DFDY
WRITE(21,10128)KK,KK
write(21,*)

!write(21,10122)KK,1,1
!write(21,10122)KK,1,1

write(21,101)
! SUBROUTINE DMDC(A, KK, DM, HML, WT, RHO, CPB)



! SUBROUTINE DMDC_OF_C(A, KK, DM, HML, RHO, CPB)
write(21,120)

write(21,*)
write(21,10095)KK+1

NTHSP=0
do i=1,II
    ! IF third-body reaction
    if (ITHB(i) > 0)then
        CALL CKINU (i, NDIM, IWORK, RWORK, NSPEC, KI, NU)
! KI(N) is the index of the Nth species in reaction I
! NU(N) is the stoichiometric coefficient of the Nth species in reaction I
        DO K=1,NSPEC
            NTHSP = NTHSP+1
            NSPinTH(NTHSP)=KI(K)
        enddo
    endif
enddo
!write(*,*)'NSPinTH = ',NSPinTH
call delSameNumbANDsort(NTHSP,NSPinTH,NNTHSP,NSPinTH2)
!write(*,*)'NSPinTH2 = ',NSPinTH2
DO K=1,NNTHSP
!    write(21,10096)NSPinTH2(K)+1,NSPinTH2(K)+1,NSPinTH2(K),NSPinTH2(K)
    write(21,10096)NSPinTH2(K)+1,NSPinTH2(K)+1,NSPinTH2(K)
ENDDO
WRITE(21,10102)
WRITE(21,*)

! A(1,k+1)
WRITE(21,11129)KK,KK
write(21,*)

write(21,101)
! SUBROUTINE DMDC_OF_C(A, KK, DM, HML, RHO, CPB)


!! SUBROUTINE DWDCT(T, C, A, KK)
!write(21,113)
!write(21,10101)KK,KK,KK
!write(21,10097)KK
!write(21,10098)KK+1,KK+1 ! A(KK+1,KK+1)
!write(21,10099)KK
!write(21,10100)KK
!write(21,101)
!! SUBROUTINE DWDCT(T, C, A, KK)


! SUBROUTINE JAC (NEQ, T, Y, ML, MU, PD, RPAR, IPAR)
write(21,114)
WRITE(21,10120)KK

write(21,101)
! SUBROUTINE JAC (NEQ, T, Y, ML, MU, PD, RPAR, IPAR)


! SUBROUTINE JAC_OF_Y 
write(21,115)
WRITE(21,10137)
WRITE(21,20041)

write(21,101)
! SUBROUTINE JAC_OF_Y 


! SUBROUTINE JAC_OF_C 
write(21,118)
WRITE(21,20137)
WRITE(21,20041)

write(21,101)
! SUBROUTINE JAC_OF_C 



! SUBROUTINE omega_OF_Y(T, C, A, KK, G0, DG, DM, CTOT, P, PL)
write(21,116)
write(21,11104)KK

!LLIcRV = 1 ! starts an array of reaction indices for those with explicit reverse Arrhenius coefficients
!LIcRV = 1  ! starts an array of reaction indices for those with explicit reverse Arrhenius coefficients
!LLLIcRV = 1  ! starts an array of reaction indices for those with explicit reverse Arrhenius coefficients
!LT=1 ! Landau-Teller reactions index

DO i=1,II
    call CKSYMR (i, LOUT, IWORK, RWORK, CWORK, LT, ISTR,KERR)
    write(21,10017)i,ISTR
    CALL CKINU (i, NDIM, IWORK, RWORK, NSPEC, KI, NU)
    
    NReacSpec = 0
    NProdSpec = 0
    DO K=1,NSPEC
        IF(NU(k)<0)THEN
            NReacSpec = NReacSpec+1
            KNReacSpec(NReacSpec) = KI(k) ! reactant indice in reaction i
            KNUReacSpec(KNReacSpec(NReacSpec)) = NU(k) ! stoichiometric coefficient of reactant in reaction i
        ENDIF
    ENDDO
    DO K=1,NSPEC
        IF(NU(k)>0)THEN
            NProdSpec = NProdSpec+1
            KNProdSpec(NProdSpec) = KI(k)! product indice in reaction i
            KNUProdSpec(KNProdSpec(NProdSpec)) = NU(k)! stoichiometric coefficient of product in reaction i 
        ENDIF
    ENDDO
! KI(N) is the index of the Nth species in reaction I
! NU(N) is the stoichiometric coefficient of the Nth species in reaction I
! IREV(i) > 0 -- reversible reaction
    IF (IREV(i) > 0 .AND. IcRVIcRV(i) == 0) THEN ! if reversible but without explicit reverse Arrhenius coefficients
        write(21,10050,advance='no') ! G0_SUM
        DO k=1,IREV(i)
            write(21,10051,advance='no')NU(k),KI(k)
        ENDDO

        NUSUM = 0
        DO k=1,IREV(i)
            NUSUM = NUSUM + NU(K)
        ENDDO
        write(21,*)
        write(21,10155)NUSUM! DGDT
    ENDIF
    
    NRFArru = NNcCO + (i-1)*(NNPAR+1)! NNPAR = 3
    IF(IcLTIcLT(i) == 0)THEN
        if(RWORK(NRFArru) > SMALL)then
            WRITE(21,10054)log(RWORK(NRFArru)),RWORK(NRFArru+1),-RWORK(NRFArru+2)! RWORK(NRFArru+2) = Ea/1.987 (unit K) ! RF
        else
            WRITE(21,11133)
        endif
    ENDIF
    
    
    if (ITHB(i) >= 0) then! if three body reaction 
        write(21,10059,advance='no')
        KTHBinI = 0
        DO K=1,KK
            IF (AKI(k,i) .NE. 1.D0)THEN
                KTHBinI = KTHBinI + 1
                IF (MOD(KTHBinI,5) == 0)THEN
                    write(21,10127)
                ENDIF
                write(21,10060,advance='no')AKI(k,i)-1.D0,K ! CM=CTOT+[M]
            ENDIF
        ENDDO
        WRITE(21,*)
        if ( IFOP(i) > 0) then ! if pressure depend reaction
            IF (IREV(i) > 0)THEN ! if reversible
                write(21,10058)log(FPAR(1,i)),FPAR(2,i),-FPAR(3,i) ! RL
                write(21,10161)! RB = RF/EQ; PR = RL*CM/RF

                if (IFOP(i) == 3) then
                    write(21,11063) ! PRLGDM 'PC = PR/(1+PR)' 'PCLGDT = PRLGDT/TEMP1' 'PCLGDM = PRLGDM/TEMP1' 'PRLG = LOG(MAX(PR, SMALL))')
                    IF(ABS(FPAR(7,i)) > 1.D+90)THEN
                        write(21,11064)1-FPAR(4,i),FPAR(5,i),FPAR(4,i),FPAR(6,i),-FPAR(7,i)/1.D10 ! TEMP1= ;TEMP2= ;TEMP3= ;
                    ELSE
                        write(21,10064)1-FPAR(4,i),FPAR(5,i),FPAR(4,i),FPAR(6,i),-FPAR(7,i) ! TEMP1= ;TEMP2= ;TEMP3= ;
                    ENDIF

                    write(21,11065)
                    write(21,11066)
                elseif (IFOP(i) == 4) then
                    write(21,11063) ! PRLGDM 'PC = PR/(1+PR)' 'PCLGDT = PRLGDT/TEMP1' 'PCLGDM = PRLGDM/TEMP1' 'PRLG = LOG(MAX(PR, SMALL))')
                    IF(ABS(FPAR(7,i)) > 1.D+90)THEN
                        write(21,11064)1-FPAR(4,i),FPAR(5,i),FPAR(4,i),FPAR(6,i),-FPAR(7,i)/1.D10 ! TEMP1= ;TEMP2= ;TEMP3= ;
                    ELSE
                        write(21,10064)1-FPAR(4,i),FPAR(5,i),FPAR(4,i),FPAR(6,i),-FPAR(7,i) ! TEMP1= ;TEMP2= ;TEMP3= ;
                    ENDIF
                    write(21,11065)
                    write(21,11066)
                elseif (IFOP(i) == 2) then
                    write(21,*)'ERROR--need to supply SRI fomat output'
                elseif (IFOP(i) == 1) then
                    write(21,11126)
                endif
                
                WRITE(21,10075,advance='no') ! WF = RF * C(X)**',I4
                DO K=1,NSPEC
                    IF(NU(k)<0)THEN
                        WRITE(21,10076,advance='no')KI(k),ABS(NU(k)) ! *C(X)**
                    ENDIF
                ENDDO

                write(21,*)
				WRITE(21,10078,advance='no') ! WB = RB
                DO K=1,NSPEC
                    IF(NU(k)>0)THEN
                        WRITE(21,10076,advance='no')KI(k),ABS(NU(k))
                    ENDIF
                ENDDO
				
				write(21,*)  
                DO K=1,NSPEC
                    IF(NU(k)>0)THEN
                        WRITE(21,20039)KI(k),KI(k),NU(k) ! WDOT(k) = WDOT(k) + (WF-WB)
					ELSE
						WRITE(21,20040)KI(k),KI(k),ABS(NU(k)) ! WDOT(k) = WDOT(k) - (WF-WB)
                    ENDIF
                ENDDO
				write(21,*)
                
    ! KNReacSpec(NReacSpec) =  reactant indice in reaction i
    ! KNUReacSpec(KNReacSpec(NReacSpec)) = stoichiometric coefficient of reactant in reaction i  

            ELSEIF(IREV(i) < 0)THEN
                write(21,10058)log(FPAR(1,i)),FPAR(2,i),-FPAR(3,i) ! RL
                WRITE(21,10069)! RB = 0
                write(21,11070)! CM = MAX(CM, SMALL), PR = RL*CM/RF'

                if (IFOP(i) == 3) then
                    write(21,11063) ! PRLGDM 'PC = PR/(1+PR)' 'PCLGDT = PRLGDT/TEMP1' 'PCLGDM = PRLGDM/TEMP1' 'PRLG = LOG(MAX(PR, SMALL))')
                    IF(ABS(FPAR(7,i)) > 1.D+90)THEN
                        write(21,11064)1-FPAR(4,i),FPAR(5,i),FPAR(4,i),FPAR(6,i),-FPAR(7,i)/1.D10 ! TEMP1= ;TEMP2= ;TEMP3= ;
                    ELSE
                        write(21,10064)1-FPAR(4,i),FPAR(5,i),FPAR(4,i),FPAR(6,i),-FPAR(7,i) ! TEMP1= ;TEMP2= ;TEMP3= ;
                    ENDIF
                    write(21,11065)
                    write(21,11066)
                elseif (IFOP(i) == 4) then
                    write(21,11063) ! PRLGDM 'PC = PR/(1+PR)' 'PCLGDT = PRLGDT/TEMP1' 'PCLGDM = PRLGDM/TEMP1' 'PRLG = LOG(MAX(PR, SMALL))')
                    IF(ABS(FPAR(7,i)) > 1.D+90)THEN
                        write(21,11064)1-FPAR(4,i),FPAR(5,i),FPAR(4,i),FPAR(6,i),-FPAR(7,i)/1.D10 ! TEMP1= ;TEMP2= ;TEMP3= ;
                    ELSE
                        write(21,10064)1-FPAR(4,i),FPAR(5,i),FPAR(4,i),FPAR(6,i),-FPAR(7,i) ! TEMP1= ;TEMP2= ;TEMP3= ;
                    ENDIF
                    write(21,11065)
                    write(21,11066)
                elseif (IFOP(i) == 2) then
                    write(21,*)'ERROR--need to supply SRI fomat output'
                elseif (IFOP(i) == 1) then
                    write(21,11126)
                endif
                
                WRITE(21,10075,advance='no') ! WF = RF * C(X)**',I4
                DO K=1,NSPEC
                    IF(NU(k)<0)THEN
                        WRITE(21,10076,advance='no')KI(k),ABS(NU(k)) ! *C(X)**
                    ENDIF
                ENDDO

                write(21,*)
				WRITE(21,10078,advance='no') ! WB = RB
                DO K=1,NSPEC
                    IF(NU(k)>0)THEN
                        WRITE(21,10076,advance='no')KI(k),ABS(NU(k))
                    ENDIF
                ENDDO
				write(21,*)
				DO K=1,NSPEC
                    IF(NU(k)>0)THEN
                        WRITE(21,20039)KI(k),KI(k),NU(k) ! WDOT(k) = WDOT(k) + (WF-WB)
					ELSE
						WRITE(21,20040)KI(k),KI(k),ABS(NU(k)) ! WDOT(k) = WDOT(k) - (WF-WB)
                    ENDIF
                ENDDO
				write(21,*)
				
                
    ! KNReacSpec(NReacSpec) =  reactant indice in reaction i
    ! KNUReacSpec(KNReacSpec(NReacSpec)) = stoichiometric coefficient of reactant in reaction i  


!                write(*,*)'**ERROR--need to supply pressure depend and irreversible reaction'
!                write(*,*)'i = ',i
            ENDIF
        else ! if three body reaction but not pressure depend reaction
            IF (IREV(i) > 0) THEN  ! if reverse reaction
!                IF (i == IWORK(IIcRV+LIcRV-1) .AND. LIcRV <= NNREV) THEN ! if reaction i with explicit reverse Arrhenius coefficients
                if (IcLTIcLT(i) == 0)then ! if not Landau-Teller reaction
                    IF (IcRVIcRV(i)>0) THEN ! if reaction i with explicit reverse Arrhenius coefficients
                        NRBArru = NNcRV + (IcRVNRV(i)-1)*(NNPAR+1) ! NNPAR = 3
                        if(RWORK(NRBArru)>SMALL)then
                            WRITE(21,10067)log(RWORK(NRBArru)),RWORK(NRBArru+1),-RWORK(NRBArru+2)! RB
                        else
                            write(21,10069)

                        endif
                        WRITE(21,11136)! PC = CM, PCLGDM = 1.D0/CM, RF = RF*PC, RB = RB*PC
                        
                        WRITE(21,10075,advance='no')
                        DO K=1,NSPEC
                            IF(NU(k)<0)THEN
                                WRITE(21,10076,advance='no')KI(k),ABS(NU(k))
                            ENDIF
                        ENDDO

                        write(21,*)
						WRITE(21,10078,advance='no')
                        DO K=1,NSPEC
                            IF(NU(k)>0)THEN
                                WRITE(21,10076,advance='no')KI(k),ABS(NU(k))
                            ENDIF
                        ENDDO
						write(21,*)
						DO K=1,NSPEC
							IF(NU(k)>0)THEN
								WRITE(21,20039)KI(k),KI(k),NU(k) ! WDOT(k) = WDOT(k) + (WF-WB)
							ELSE
								WRITE(21,20040)KI(k),KI(k),ABS(NU(k)) ! WDOT(k) = WDOT(k) - (WF-WB)
							ENDIF
						ENDDO
						write(21,*)

                    ELSE ! if reaction i without explicit reverse Arrhenius coefficients
                        WRITE(21,11071)
                        WRITE(21,10075,advance='no')
                        DO K=1,NSPEC
                            IF(NU(k)<0)THEN
                                WRITE(21,10076,advance='no')KI(k),ABS(NU(k))
                            ENDIF
                        ENDDO

                        write(21,*)
						WRITE(21,10078,advance='no')
                        DO K=1,NSPEC
                            IF(NU(k)>0)THEN
                                WRITE(21,10076,advance='no')KI(k),ABS(NU(k))
                            ENDIF
                        ENDDO
						write(21,*)
						DO K=1,NSPEC
							IF(NU(k)>0)THEN
								WRITE(21,20039)KI(k),KI(k),NU(k) ! WDOT(k) = WDOT(k) + (WF-WB)
							ELSE
								WRITE(21,20040)KI(k),KI(k),ABS(NU(k)) ! WDOT(k) = WDOT(k) - (WF-WB)
							ENDIF
						ENDDO
						write(21,*)
                        
                    ENDIF
                elseif (IcLTIcLT(i) > 0)then ! if Landau-Teller
                    write(*,*)'ERROR--need to supply Landau-Teller and pressure depend and reverse reaction'
                endif
            ELSE ! if not reverse reaction
                WRITE(21,10069)! RB = 0

                WRITE(21,11072)
                
                WRITE(21,10075,advance='no')
                DO K=1,NSPEC
                    IF(NU(k)<0)THEN
                        WRITE(21,10076,advance='no')KI(k),ABS(NU(k))
                    ENDIF
                ENDDO

                write(21,*)
				WRITE(21,10078,advance='no')
                DO K=1,NSPEC
                    IF(NU(k)>0)THEN
                        WRITE(21,10076,advance='no')KI(k),ABS(NU(k))
                    ENDIF
                ENDDO
				write(21,*)
				DO K=1,NSPEC
                    IF(NU(k)>0)THEN
                        WRITE(21,20039)KI(k),KI(k),NU(k) ! WDOT(k) = WDOT(k) + (WF-WB)
					ELSE
						WRITE(21,20040)KI(k),KI(k),ABS(NU(k)) ! WDOT(k) = WDOT(k) - (WF-WB)
                    ENDIF
                ENDDO
				write(21,*)
            ENDIF
        endif
    else ! if not three body reaction 
!        IF (i == IWORK(IIcLT+LT-1) .AND. LT <= NNLAN)THEN ! IF not three body reactions AND Landau-Teller reactions
        IF (IcLTIcLT(i)>0)THEN ! IF not three body reactions AND Landau-Teller reactions
            NRFLT = NNcLT + (IcLTNLT(i)-1)*NNLAR ! NNLAR = 4
            if(RWORK(NRFArru) > SMALL)then
                WRITE(21,10074)log(RWORK(NRFArru)),RWORK(NRFArru+1),-RWORK(NRFArru+2),& ! RF
                            RWORK(NRFLT),RWORK(NRFLT+1)
            else
                write(21,11133)
            endif
            
            IF (IREV(i) > 0) THEN  ! if reverse reaction
                WRITE(21,11073) ! RB = RF*EQINV; RBLGDT = RFLGDT +DGDT
                
                WRITE(21,10075,advance='no')
                DO K=1,NSPEC
                    IF(NU(k)<0)THEN
                        WRITE(21,10076,advance='no')KI(k),ABS(NU(k))
                    ENDIF
                ENDDO

                write(21,*)
				WRITE(21,10078,advance='no')
                DO K=1,NSPEC
                    IF(NU(k)>0)THEN
                        WRITE(21,10076,advance='no')KI(k),ABS(NU(k))
                    ENDIF
                ENDDO
				write(21,*)
				DO K=1,NSPEC
                    IF(NU(k)>0)THEN
                        WRITE(21,20039)KI(k),KI(k),NU(k) ! WDOT(k) = WDOT(k) + (WF-WB)
					ELSE
						WRITE(21,20040)KI(k),KI(k),ABS(NU(k)) ! WDOT(k) = WDOT(k) - (WF-WB)
                    ENDIF
                ENDDO
				write(21,*)

 
            ELSE! if Landau-Teller and irreverse reaction
                WRITE(21,10069)! RB = 0
                
                WRITE(21,10075,advance='no')
                DO K=1,NSPEC
                    IF(NU(k)<0)THEN
                        WRITE(21,10076,advance='no')KI(k),ABS(NU(k))
                    ENDIF
                ENDDO

                write(21,*)
				WRITE(21,10078,advance='no')
                DO K=1,NSPEC
                    IF(NU(k)>0)THEN
                        WRITE(21,10076,advance='no')KI(k),ABS(NU(k))
                    ENDIF
                ENDDO
				write(21,*)
				DO K=1,NSPEC
                    IF(NU(k)>0)THEN
                        WRITE(21,20039)KI(k),KI(k),NU(k) ! WDOT(k) = WDOT(k) + (WF-WB)
					ELSE
						WRITE(21,20040)KI(k),KI(k),ABS(NU(k)) ! WDOT(k) = WDOT(k) - (WF-WB)
                    ENDIF
                ENDDO
				write(21,*)


            ENDIF
!        IF (i == IWORK(IIcLT+LT-1) .AND. LT <= NNLAN) THEN ! IF Landau-Teller reactions
        ELSE ! IF not Landau-Teller reactions
            IF (IREV(i) > 0) THEN  ! if reverse reaction
                IF (IcRVIcRV(i)>0) THEN ! if reaction i with explicit reverse Arrhenius coefficients
!                IF (i == IWORK(IIcRV+LLLIcRV-1) .AND. LLLIcRV <= NNREV) THEN ! if reaction i with explicit reverse Arrhenius coefficients
                    NRBArru = NNcRV + (IcRVNRV(i)-1)*(NNPAR+1) ! NNPAR = 3
                    if(RWORK(NRBArru)>SMALL)then
                        WRITE(21,10067)log(RWORK(NRBArru)),RWORK(NRBArru+1),-RWORK(NRBArru+2)! RB
                    else
                        WRITE(21,10069)
                    endif
!                    LLLIcRV = LLLIcRV + 1
                    
                    WRITE(21,10075,advance='no')
                    DO K=1,NSPEC
                        IF(NU(k)<0)THEN
                            WRITE(21,10076,advance='no')KI(k),ABS(NU(k))
                        ENDIF
                    ENDDO

                    write(21,*)
					WRITE(21,10078,advance='no')
                    DO K=1,NSPEC
                        IF(NU(k)>0)THEN
                            WRITE(21,10076,advance='no')KI(k),ABS(NU(k))
                        ENDIF
                    ENDDO
					write(21,*)
					DO K=1,NSPEC
						IF(NU(k)>0)THEN
							WRITE(21,20039)KI(k),KI(k),NU(k) ! WDOT(k) = WDOT(k) + (WF-WB)
						ELSE
							WRITE(21,20040)KI(k),KI(k),ABS(NU(k)) ! WDOT(k) = WDOT(k) - (WF-WB)
						ENDIF
					ENDDO
					write(21,*)

                ELSE ! if reaction i without explicit reverse Arrhenius coefficients
                    WRITE(21,11073)! RB = RF*EQINV; RBLGDT = RFLGDT +DGDT
                    
                    WRITE(21,10075,advance='no')
                    DO K=1,NSPEC
                        IF(NU(k)<0)THEN
                            WRITE(21,10076,advance='no')KI(k),ABS(NU(k))
                        ENDIF
                    ENDDO

                    write(21,*)
					WRITE(21,10078,advance='no')
                    DO K=1,NSPEC
                        IF(NU(k)>0)THEN
                            WRITE(21,10076,advance='no')KI(k),ABS(NU(k))
                        ENDIF
                    ENDDO
					write(21,*)
					DO K=1,NSPEC
						IF(NU(k)>0)THEN
							WRITE(21,20039)KI(k),KI(k),NU(k) ! WDOT(k) = WDOT(k) + (WF-WB)
						ELSE
							WRITE(21,20040)KI(k),KI(k),ABS(NU(k)) ! WDOT(k) = WDOT(k) - (WF-WB)
						ENDIF
					ENDDO
					write(21,*)
                    
                ENDIF
            ELSE ! if not reverse reaction
                WRITE(21,10069)! RB = 0
                
                WRITE(21,10075,advance='no')
                DO K=1,NSPEC
                    IF(NU(k)<0)THEN
                        WRITE(21,10076,advance='no')KI(k),ABS(NU(k))
                    ENDIF
                ENDDO

                write(21,*)
				WRITE(21,10078,advance='no')
                DO K=1,NSPEC
                    IF(NU(k)>0)THEN
                        WRITE(21,10076,advance='no')KI(k),ABS(NU(k))
                    ENDIF
                ENDDO
				write(21,*)
				 DO K=1,NSPEC
                    IF(NU(k)>0)THEN
                        WRITE(21,20039)KI(k),KI(k),NU(k) ! WDOT(k) = WDOT(k) + (WF-WB)
					ELSE
						WRITE(21,20040)KI(k),KI(k),ABS(NU(k)) ! WDOT(k) = WDOT(k) - (WF-WB)
                    ENDIF
                ENDDO
				write(21,*)

            ENDIF
        ENDIF
    endif
ENDDO

write(21,20042)KK

write(21,101)
! SUBROUTINE omega_OF_Y(T, C, A, KK, G0, DG, DM, CTOT, P, PL)



! SUBROUTINE omega_OF_C(T, C, A, KK, G0, DG, DM, CTOT, P, PL)
write(21,121)
write(21,11114)KK

!LLIcRV = 1 ! starts an array of reaction indices for those with explicit reverse Arrhenius coefficients
!LIcRV = 1  ! starts an array of reaction indices for those with explicit reverse Arrhenius coefficients
!LLLIcRV = 1  ! starts an array of reaction indices for those with explicit reverse Arrhenius coefficients
!LT=1 ! Landau-Teller reactions index

DO i=1,II
    call CKSYMR (i, LOUT, IWORK, RWORK, CWORK, LT, ISTR,KERR)
    write(21,10017)i,ISTR
    CALL CKINU (i, NDIM, IWORK, RWORK, NSPEC, KI, NU)
    
    NReacSpec = 0
    NProdSpec = 0
    DO K=1,NSPEC
        IF(NU(k)<0)THEN
            NReacSpec = NReacSpec+1
            KNReacSpec(NReacSpec) = KI(k) ! reactant indice in reaction i
            KNUReacSpec(KNReacSpec(NReacSpec)) = NU(k) ! stoichiometric coefficient of reactant in reaction i
        ENDIF
    ENDDO
    DO K=1,NSPEC
        IF(NU(k)>0)THEN
            NProdSpec = NProdSpec+1
            KNProdSpec(NProdSpec) = KI(k)! product indice in reaction i
            KNUProdSpec(KNProdSpec(NProdSpec)) = NU(k)! stoichiometric coefficient of product in reaction i 
        ENDIF
    ENDDO
! KI(N) is the index of the Nth species in reaction I
! NU(N) is the stoichiometric coefficient of the Nth species in reaction I
! IREV(i) > 0 -- reversible reaction
    IF (IREV(i) > 0 .AND. IcRVIcRV(i) == 0) THEN ! if reversible but without explicit reverse Arrhenius coefficients
        write(21,10050,advance='no') ! G0_SUM
        DO k=1,IREV(i)
            write(21,10051,advance='no')NU(k),KI(k)
        ENDDO

        NUSUM = 0
        DO k=1,IREV(i)
            NUSUM = NUSUM + NU(K)
        ENDDO
        write(21,*)
        write(21,10155)NUSUM! DGDT
    ENDIF
    
    NRFArru = NNcCO + (i-1)*(NNPAR+1)! NNPAR = 3
    IF(IcLTIcLT(i) == 0)THEN
        if(RWORK(NRFArru) > SMALL)then
            WRITE(21,10054)log(RWORK(NRFArru)),RWORK(NRFArru+1),-RWORK(NRFArru+2)! RWORK(NRFArru+2) = Ea/1.987 (unit K) ! RF
        else
            WRITE(21,11133)
        endif
    ENDIF
    
    
    if (ITHB(i) >= 0) then! if three body reaction 
        write(21,10059,advance='no')
        KTHBinI = 0
        DO K=1,KK
            IF (AKI(k,i) .NE. 1.D0)THEN
                KTHBinI = KTHBinI + 1
                IF (MOD(KTHBinI,5) == 0)THEN
                    write(21,10127)
                ENDIF
                write(21,10060,advance='no')AKI(k,i)-1.D0,K ! CM=CTOT+[M]
            ENDIF
        ENDDO
        WRITE(21,*)
        if ( IFOP(i) > 0) then ! if pressure depend reaction
            IF (IREV(i) > 0)THEN ! if reversible
                write(21,10058)log(FPAR(1,i)),FPAR(2,i),-FPAR(3,i) ! RL
                write(21,10161)! RB = RF/EQ; PR = RL*CM/RF

                if (IFOP(i) == 3) then
                    write(21,11063) ! PRLGDM 'PC = PR/(1+PR)' 'PCLGDT = PRLGDT/TEMP1' 'PCLGDM = PRLGDM/TEMP1' 'PRLG = LOG(MAX(PR, SMALL))')
                    IF(ABS(FPAR(7,i)) > 1.D+90)THEN
                        write(21,11064)1-FPAR(4,i),FPAR(5,i),FPAR(4,i),FPAR(6,i),-FPAR(7,i)/1.D10 ! TEMP1= ;TEMP2= ;TEMP3= ;
                    ELSE
                        write(21,10064)1-FPAR(4,i),FPAR(5,i),FPAR(4,i),FPAR(6,i),-FPAR(7,i) ! TEMP1= ;TEMP2= ;TEMP3= ;
                    ENDIF
                    write(21,11065)
                    write(21,11066)
                elseif (IFOP(i) == 4) then
                    write(21,11063) ! PRLGDM 'PC = PR/(1+PR)' 'PCLGDT = PRLGDT/TEMP1' 'PCLGDM = PRLGDM/TEMP1' 'PRLG = LOG(MAX(PR, SMALL))')
                    IF(ABS(FPAR(7,i)) > 1.D+90)THEN
                        write(21,11064)1-FPAR(4,i),FPAR(5,i),FPAR(4,i),FPAR(6,i),-FPAR(7,i)/1.D10 ! TEMP1= ;TEMP2= ;TEMP3= ;
                    ELSE
                        write(21,10064)1-FPAR(4,i),FPAR(5,i),FPAR(4,i),FPAR(6,i),-FPAR(7,i) ! TEMP1= ;TEMP2= ;TEMP3= ;
                    ENDIF
                    write(21,11065)
                    write(21,11066)
                elseif (IFOP(i) == 2) then
                    write(21,*)'ERROR--need to supply SRI fomat output'
                elseif (IFOP(i) == 1) then
                    write(21,11126)
                endif
                
                WRITE(21,10075,advance='no') ! WF = RF * C(X)**',I4
                DO K=1,NSPEC
                    IF(NU(k)<0)THEN
                        WRITE(21,10076,advance='no')KI(k),ABS(NU(k)) ! *C(X)**
                    ENDIF
                ENDDO

                write(21,*)
				WRITE(21,10078,advance='no') ! WB = RB
                DO K=1,NSPEC
                    IF(NU(k)>0)THEN
                        WRITE(21,10076,advance='no')KI(k),ABS(NU(k))
                    ENDIF
                ENDDO
				
				write(21,*)  
                DO K=1,NSPEC
                    IF(NU(k)>0)THEN
                        WRITE(21,20039)KI(k),KI(k),NU(k) ! WDOT(k) = WDOT(k) + (WF-WB)
					ELSE
						WRITE(21,20040)KI(k),KI(k),ABS(NU(k)) ! WDOT(k) = WDOT(k) - (WF-WB)
                    ENDIF
                ENDDO
				write(21,*)
                
    ! KNReacSpec(NReacSpec) =  reactant indice in reaction i
    ! KNUReacSpec(KNReacSpec(NReacSpec)) = stoichiometric coefficient of reactant in reaction i  

            ELSEIF(IREV(i) < 0)THEN
                write(21,10058)log(FPAR(1,i)),FPAR(2,i),-FPAR(3,i) ! RL
                WRITE(21,10069)! RB = 0
                write(21,11070)! CM = MAX(CM, SMALL), PR = RL*CM/RF'

                if (IFOP(i) == 3) then
                    write(21,11063) ! PRLGDM 'PC = PR/(1+PR)' 'PCLGDT = PRLGDT/TEMP1' 'PCLGDM = PRLGDM/TEMP1' 'PRLG = LOG(MAX(PR, SMALL))')
                    IF(ABS(FPAR(7,i)) > 1.D+90)THEN
                        write(21,11064)1-FPAR(4,i),FPAR(5,i),FPAR(4,i),FPAR(6,i),-FPAR(7,i)/1.D10 ! TEMP1= ;TEMP2= ;TEMP3= ;
                    ELSE
                        write(21,10064)1-FPAR(4,i),FPAR(5,i),FPAR(4,i),FPAR(6,i),-FPAR(7,i) ! TEMP1= ;TEMP2= ;TEMP3= ;
                    ENDIF
                    write(21,11065)
                    write(21,11066)
                elseif (IFOP(i) == 4) then
                    write(21,11063) ! PRLGDM 'PC = PR/(1+PR)' 'PCLGDT = PRLGDT/TEMP1' 'PCLGDM = PRLGDM/TEMP1' 'PRLG = LOG(MAX(PR, SMALL))')
                    IF(ABS(FPAR(7,i)) > 1.D+90)THEN
                        write(21,11064)1-FPAR(4,i),FPAR(5,i),FPAR(4,i),FPAR(6,i),-FPAR(7,i)/1.D10 ! TEMP1= ;TEMP2= ;TEMP3= ;
                    ELSE
                        write(21,10064)1-FPAR(4,i),FPAR(5,i),FPAR(4,i),FPAR(6,i),-FPAR(7,i) ! TEMP1= ;TEMP2= ;TEMP3= ;
                    ENDIF
                    write(21,11065)
                    write(21,11066)
                elseif (IFOP(i) == 2) then
                    write(21,*)'ERROR--need to supply SRI fomat output'
                elseif (IFOP(i) == 1) then
                    write(21,11126)
                endif
                
                WRITE(21,10075,advance='no') ! WF = RF * C(X)**',I4
                DO K=1,NSPEC
                    IF(NU(k)<0)THEN
                        WRITE(21,10076,advance='no')KI(k),ABS(NU(k)) ! *C(X)**
                    ENDIF
                ENDDO

                write(21,*)
				WRITE(21,10078,advance='no') ! WB = RB
                DO K=1,NSPEC
                    IF(NU(k)>0)THEN
                        WRITE(21,10076,advance='no')KI(k),ABS(NU(k))
                    ENDIF
                ENDDO
				write(21,*)
				DO K=1,NSPEC
                    IF(NU(k)>0)THEN
                        WRITE(21,20039)KI(k),KI(k),NU(k) ! WDOT(k) = WDOT(k) + (WF-WB)
					ELSE
						WRITE(21,20040)KI(k),KI(k),ABS(NU(k)) ! WDOT(k) = WDOT(k) - (WF-WB)
                    ENDIF
                ENDDO
				write(21,*)
				
                
    ! KNReacSpec(NReacSpec) =  reactant indice in reaction i
    ! KNUReacSpec(KNReacSpec(NReacSpec)) = stoichiometric coefficient of reactant in reaction i  


!                write(*,*)'**ERROR--need to supply pressure depend and irreversible reaction'
!                write(*,*)'i = ',i
            ENDIF
        else ! if three body reaction but not pressure depend reaction
            IF (IREV(i) > 0) THEN  ! if reverse reaction
!                IF (i == IWORK(IIcRV+LIcRV-1) .AND. LIcRV <= NNREV) THEN ! if reaction i with explicit reverse Arrhenius coefficients
                if (IcLTIcLT(i) == 0)then ! if not Landau-Teller reaction
                    IF (IcRVIcRV(i)>0) THEN ! if reaction i with explicit reverse Arrhenius coefficients
                        NRBArru = NNcRV + (IcRVNRV(i)-1)*(NNPAR+1) ! NNPAR = 3
                        if(RWORK(NRBArru)>SMALL)then
                            WRITE(21,10067)log(RWORK(NRBArru)),RWORK(NRBArru+1),-RWORK(NRBArru+2)! RB
                        else
                            write(21,10069)

                        endif
                        WRITE(21,11136)! PC = CM, PCLGDM = 1.D0/CM, RF = RF*PC, RB = RB*PC
                        
                        WRITE(21,10075,advance='no')
                        DO K=1,NSPEC
                            IF(NU(k)<0)THEN
                                WRITE(21,10076,advance='no')KI(k),ABS(NU(k))
                            ENDIF
                        ENDDO

                        write(21,*)
						WRITE(21,10078,advance='no')
                        DO K=1,NSPEC
                            IF(NU(k)>0)THEN
                                WRITE(21,10076,advance='no')KI(k),ABS(NU(k))
                            ENDIF
                        ENDDO
						write(21,*)
						DO K=1,NSPEC
							IF(NU(k)>0)THEN
								WRITE(21,20039)KI(k),KI(k),NU(k) ! WDOT(k) = WDOT(k) + (WF-WB)
							ELSE
								WRITE(21,20040)KI(k),KI(k),ABS(NU(k)) ! WDOT(k) = WDOT(k) - (WF-WB)
							ENDIF
						ENDDO
						write(21,*)

                    ELSE ! if reaction i without explicit reverse Arrhenius coefficients
                        WRITE(21,11071)
                        WRITE(21,10075,advance='no')
                        DO K=1,NSPEC
                            IF(NU(k)<0)THEN
                                WRITE(21,10076,advance='no')KI(k),ABS(NU(k))
                            ENDIF
                        ENDDO

                        write(21,*)
						WRITE(21,10078,advance='no')
                        DO K=1,NSPEC
                            IF(NU(k)>0)THEN
                                WRITE(21,10076,advance='no')KI(k),ABS(NU(k))
                            ENDIF
                        ENDDO
						write(21,*)
						DO K=1,NSPEC
							IF(NU(k)>0)THEN
								WRITE(21,20039)KI(k),KI(k),NU(k) ! WDOT(k) = WDOT(k) + (WF-WB)
							ELSE
								WRITE(21,20040)KI(k),KI(k),ABS(NU(k)) ! WDOT(k) = WDOT(k) - (WF-WB)
							ENDIF
						ENDDO
						write(21,*)
                        
                    ENDIF
                elseif (IcLTIcLT(i) > 0)then ! if Landau-Teller
                    write(*,*)'ERROR--need to supply Landau-Teller and pressure depend and reverse reaction'
                endif
            ELSE ! if not reverse reaction
                WRITE(21,10069)! RB = 0

                WRITE(21,11072)
                
                WRITE(21,10075,advance='no')
                DO K=1,NSPEC
                    IF(NU(k)<0)THEN
                        WRITE(21,10076,advance='no')KI(k),ABS(NU(k))
                    ENDIF
                ENDDO

                write(21,*)
				WRITE(21,10078,advance='no')
                DO K=1,NSPEC
                    IF(NU(k)>0)THEN
                        WRITE(21,10076,advance='no')KI(k),ABS(NU(k))
                    ENDIF
                ENDDO
				write(21,*)
				DO K=1,NSPEC
                    IF(NU(k)>0)THEN
                        WRITE(21,20039)KI(k),KI(k),NU(k) ! WDOT(k) = WDOT(k) + (WF-WB)
					ELSE
						WRITE(21,20040)KI(k),KI(k),ABS(NU(k)) ! WDOT(k) = WDOT(k) - (WF-WB)
                    ENDIF
                ENDDO
				write(21,*)
            ENDIF
        endif
    else ! if not three body reaction 
!        IF (i == IWORK(IIcLT+LT-1) .AND. LT <= NNLAN)THEN ! IF not three body reactions AND Landau-Teller reactions
        IF (IcLTIcLT(i)>0)THEN ! IF not three body reactions AND Landau-Teller reactions
            NRFLT = NNcLT + (IcLTNLT(i)-1)*NNLAR ! NNLAR = 4
            if(RWORK(NRFArru) > SMALL)then
                WRITE(21,10074)log(RWORK(NRFArru)),RWORK(NRFArru+1),-RWORK(NRFArru+2),& ! RF
                            RWORK(NRFLT),RWORK(NRFLT+1)
            else
                write(21,11133)
            endif
            
            IF (IREV(i) > 0) THEN  ! if reverse reaction
                WRITE(21,11073) ! RB = RF*EQINV; RBLGDT = RFLGDT +DGDT
                
                WRITE(21,10075,advance='no')
                DO K=1,NSPEC
                    IF(NU(k)<0)THEN
                        WRITE(21,10076,advance='no')KI(k),ABS(NU(k))
                    ENDIF
                ENDDO

                write(21,*)
				WRITE(21,10078,advance='no')
                DO K=1,NSPEC
                    IF(NU(k)>0)THEN
                        WRITE(21,10076,advance='no')KI(k),ABS(NU(k))
                    ENDIF
                ENDDO
				write(21,*)
				DO K=1,NSPEC
                    IF(NU(k)>0)THEN
                        WRITE(21,20039)KI(k),KI(k),NU(k) ! WDOT(k) = WDOT(k) + (WF-WB)
					ELSE
						WRITE(21,20040)KI(k),KI(k),ABS(NU(k)) ! WDOT(k) = WDOT(k) - (WF-WB)
                    ENDIF
                ENDDO
				write(21,*)

 
            ELSE! if Landau-Teller and irreverse reaction
                WRITE(21,10069)! RB = 0
                
                WRITE(21,10075,advance='no')
                DO K=1,NSPEC
                    IF(NU(k)<0)THEN
                        WRITE(21,10076,advance='no')KI(k),ABS(NU(k))
                    ENDIF
                ENDDO

                write(21,*)
				WRITE(21,10078,advance='no')
                DO K=1,NSPEC
                    IF(NU(k)>0)THEN
                        WRITE(21,10076,advance='no')KI(k),ABS(NU(k))
                    ENDIF
                ENDDO
				write(21,*)
				DO K=1,NSPEC
                    IF(NU(k)>0)THEN
                        WRITE(21,20039)KI(k),KI(k),NU(k) ! WDOT(k) = WDOT(k) + (WF-WB)
					ELSE
						WRITE(21,20040)KI(k),KI(k),ABS(NU(k)) ! WDOT(k) = WDOT(k) - (WF-WB)
                    ENDIF
                ENDDO
				write(21,*)


            ENDIF
!        IF (i == IWORK(IIcLT+LT-1) .AND. LT <= NNLAN) THEN ! IF Landau-Teller reactions
        ELSE ! IF not Landau-Teller reactions
            IF (IREV(i) > 0) THEN  ! if reverse reaction
                IF (IcRVIcRV(i)>0) THEN ! if reaction i with explicit reverse Arrhenius coefficients
!                IF (i == IWORK(IIcRV+LLLIcRV-1) .AND. LLLIcRV <= NNREV) THEN ! if reaction i with explicit reverse Arrhenius coefficients
                    NRBArru = NNcRV + (IcRVNRV(i)-1)*(NNPAR+1) ! NNPAR = 3
                    if(RWORK(NRBArru)>SMALL)then
                        WRITE(21,10067)log(RWORK(NRBArru)),RWORK(NRBArru+1),-RWORK(NRBArru+2)! RB
                    else
                        WRITE(21,10069)
                    endif
!                    LLLIcRV = LLLIcRV + 1
                    
                    WRITE(21,10075,advance='no')
                    DO K=1,NSPEC
                        IF(NU(k)<0)THEN
                            WRITE(21,10076,advance='no')KI(k),ABS(NU(k))
                        ENDIF
                    ENDDO

                    write(21,*)
					WRITE(21,10078,advance='no')
                    DO K=1,NSPEC
                        IF(NU(k)>0)THEN
                            WRITE(21,10076,advance='no')KI(k),ABS(NU(k))
                        ENDIF
                    ENDDO
					write(21,*)
					DO K=1,NSPEC
						IF(NU(k)>0)THEN
							WRITE(21,20039)KI(k),KI(k),NU(k) ! WDOT(k) = WDOT(k) + (WF-WB)
						ELSE
							WRITE(21,20040)KI(k),KI(k),ABS(NU(k)) ! WDOT(k) = WDOT(k) - (WF-WB)
						ENDIF
					ENDDO
					write(21,*)

                ELSE ! if reaction i without explicit reverse Arrhenius coefficients
                    WRITE(21,11073)! RB = RF*EQINV; RBLGDT = RFLGDT +DGDT
                    
                    WRITE(21,10075,advance='no')
                    DO K=1,NSPEC
                        IF(NU(k)<0)THEN
                            WRITE(21,10076,advance='no')KI(k),ABS(NU(k))
                        ENDIF
                    ENDDO

                    write(21,*)
					WRITE(21,10078,advance='no')
                    DO K=1,NSPEC
                        IF(NU(k)>0)THEN
                            WRITE(21,10076,advance='no')KI(k),ABS(NU(k))
                        ENDIF
                    ENDDO
					write(21,*)
					DO K=1,NSPEC
						IF(NU(k)>0)THEN
							WRITE(21,20039)KI(k),KI(k),NU(k) ! WDOT(k) = WDOT(k) + (WF-WB)
						ELSE
							WRITE(21,20040)KI(k),KI(k),ABS(NU(k)) ! WDOT(k) = WDOT(k) - (WF-WB)
						ENDIF
					ENDDO
					write(21,*)
                    
                ENDIF
            ELSE ! if not reverse reaction
                WRITE(21,10069)! RB = 0
                
                WRITE(21,10075,advance='no')
                DO K=1,NSPEC
                    IF(NU(k)<0)THEN
                        WRITE(21,10076,advance='no')KI(k),ABS(NU(k))
                    ENDIF
                ENDDO

                write(21,*)
				WRITE(21,10078,advance='no')
                DO K=1,NSPEC
                    IF(NU(k)>0)THEN
                        WRITE(21,10076,advance='no')KI(k),ABS(NU(k))
                    ENDIF
                ENDDO
				write(21,*)
				 DO K=1,NSPEC
                    IF(NU(k)>0)THEN
                        WRITE(21,20039)KI(k),KI(k),NU(k) ! WDOT(k) = WDOT(k) + (WF-WB)
					ELSE
						WRITE(21,20040)KI(k),KI(k),ABS(NU(k)) ! WDOT(k) = WDOT(k) - (WF-WB)
                    ENDIF
                ENDDO
				write(21,*)

            ENDIF
        ENDIF
    endif
ENDDO

write(21,20042)KK

write(21,101)
! SUBROUTINE omega_OF_C(T, C, A, KK, G0, DG, DM, CTOT, P, PL)

     write(*,*)'****************************************************'
     IERROR = 0
     NEQ = 3
     Y(1) = 1.0D0
     Y(2) = 0.0D0
     Y(3) = 0.0D0
     T = 0.0D0
     TOUT = 0.4D0
     ITASK = 1
     ISTATE = 1
    
60  FORMAT(/'  No. steps =',I4,'   No. f-s =',I4,        &
         '  No. J-s =',I4,'   No. LU-s =',I4/         &
         '  No. nonlinear iterations =',I4/           &
         '  No. nonlinear convergence failures =',I4/ &
         '  No. error test failures =',I4/)
61  FORMAT(/' An error occurred.')
62  FORMAT(/' No errors occurred.')
63  FORMAT(' At t =',D12.4,'   y =',3D14.6)
64  FORMAT(///' Error halt: ISTATE =',I3)
7100 FORMAT (5X, 100(1X,A10))
7105 FORMAT (5X, 100(1X,I5))

100 FORMAT (/'SUBROUTINE pTY2rhoC (P, T, Y, RHO, C)'/  &
             'IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)'/ &
             'DIMENSION Y(*), C(*)'/)
101 FORMAT (/'END'/)
102 FORMAT (/'SUBROUTINE rhoTY2C (RHO, T, Y, C)'/ &
            'IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)'/ &
            'DIMENSION Y(*), C(*)' / &
            'DATA SMALL/1D-50/'/)
103 FORMAT ('SUBROUTINE SKSMH  (T, SMH)'/ &
      'IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)'/ &
      'DIMENSION SMH(*), TN(5)'/ &
      'TLOG = LOG(T)'/ &
      'TI = 1D0/T'/ &
      'TN(1) = TLOG - 1D0'/ &
      'TN(2) = T'/ &
      'TN(3) = TN(2)*T'/ &
      'TN(4) = TN(3)*T'/ &
      'TN(5) = TN(4)*T'/)

104 FORMAT ('SUBROUTINE SKRATT (T, RF, RB, RKLOW)'/&
      'IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)'/&
      'PARAMETER (RU=8.31451D7, SMALL=1.D-200, PATM=1.01325D6)'/&
      'DIMENSION RF(*), RB(*), RKLOW(*)'/&
      'DIMENSION SMH(',I4,'), EG(',I4,')'/&
      '!'/&
      'ALOGT = LOG(T)'/&
      'TI = 1D0/T'/&
      'TI2 = TI*TI'/&
      '!'/&
      'CALL SKSMH (T, SMH)'/)
105 FORMAT ('SUBROUTINE SKRATX (T, C, RF, RB, RKLOW)'/&
      'IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)'/&
      'PARAMETER (SMALL = 1D-200)'/&
      'DIMENSION C(*), RF(*), RB(*), RKLOW(*)'/&
      '!'/&
      'ALOGT = LOG(T)'/&
      'CTOT = 0.0'/&
      'DO K = 1, ',I4/&
      '  CTOT = CTOT + C(K)'/&
      'ENDDO'/)
106 FORMAT ('SUBROUTINE SKWDOT(RF, RB, WDOT)'/&
      'IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)'/&
      'DIMENSION RF(*), RB(*), WDOT(*)'/&
      '!'/&
      'DO K = 1, ',I4/&
      '  WDOT(K) = 0D0'/&
      'ENDDO'/)
107 FORMAT ('SUBROUTINE SKHML  (T, HML) ! ergs/mole'/&
      'IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)'/&
      'DIMENSION HML(*), TN(5)'/&
      '!'/&
      'TN(1) = T'/&
      'TN(2) = TN(1)*T'/&
      'TN(3) = TN(2)*T'/&
      'TN(4) = TN(3)*T'/&
      'TN(5) = TN(4)*T'/)
108 FORMAT ('SUBROUTINE SKCPML  (T, CPML) ! ergs/(mole*K)'/&
      'IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)'/&
      'DIMENSION CPML(*), TN(4)'/&
      '!'/&
      'TN(1) = T'/&
      'TN(2) = TN(1)*T'/&
      'TN(3) = TN(2)*T'/&
      'TN(4) = TN(3)*T')
109 FORMAT ('SUBROUTINE GDG1(T, G0, DG)'/&
      'IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)'/&
      'DIMENSION G0(*), DG(*)'/&
      '!'/&
      'T2 = T*T'/&
      'T3 = T2*T'/&
      'T4 = T3*T'/&
      'TLOG = LOG(T)')
110 FORMAT ('SUBROUTINE SKWT (WT)'/&
      'IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)'/&
      'DIMENSION WT(*)'/&
      '!') 
111 FORMAT ('SUBROUTINE AJ1(T, C, A, KK, G0, DG, DM, CTOT, WT, RHO, CPB, CPML, HML)'/&
      'IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)'/&
      'PARAMETER (EXPARG=690.776, SMALL=1.D-300)'/&
      'DIMENSION C(*), A(KK+1,*), G0(*), DG(*), DM(*),CPML(*),HML(*),WT(*)'/&
      'DIMENSION DTDT(KK+1),DTDC(KK+1),WK(KK)'/&
      '!'/ &
      'T2 = T*T'/&
      'TLOG = LOG(T)'/&
      'TI = 1.D0/T'/&
      'TI2 = TI/T'/&
      'PFAC = 1.218652692702276D-2/T'/&
      'PFAC2 = PFAC*PFAC')
112 FORMAT ('SUBROUTINE DMDC(A, KK, DM, HML, WT, RHO, CPB)'/&
      'IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)' /&
      'DIMENSION A(KK+1,*), DM(*), HML(*), WT(*)')
113 FORMAT ('SUBROUTINE DWDCT(T, C, A, KK)' /&
'!' /&
'!     Analytic Jacobian for species net production rate w(T, C)' /&
'!     T: input, Temperature, K, scalar' /&
'!     C: input, species mole concentration, mole/cm^3, vector, size at least KK' /&
'!     A: output, Jacobian matrix [dg/dc dg/dT], size at least KK*(KK+1)' /&
'!     KK+1: input, leading dimension of matrix A, scalar, value at least KK' /&
'!' /&
'IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)' /&
'PARAMETER (EXPARG=690.776, SMALL=1.D-300)')
114 FORMAT ('SUBROUTINE JAC (NEQ, T, Z, ML, MU, PD, NRPD, RPAR, IPAR)'/&
      'IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)'/&
      'DOUBLE PRECISION PD, RPAR, T, Z'/&
      'COMMON /ICONS/ KK, NP, NWT, NH, NWDOT'/&
      'DIMENSION Z(NEQ), PD(NEQ,NEQ),CPML(KK),HML(KK),WT(KK),&'/&
      'G0(KK),DG(KK),C(KK),RPAR(*),IPAR(*),DM(KK),DTDT(KK+1),DTDC(KK+1),WK(KK)')

115 FORMAT ('SUBROUTINE JAC_OF_Y (KK, p, T, Y, PD_OF)'/&
      'IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)'/&
      'DIMENSION Y(KK), PD(KK+1,KK+1),CPML(KK),HML(KK),WT(KK),&'/&
        'G0(KK),DG(KK),C(KK),DM(KK),PD_OF(0:(KK+1)*(KK+1))')

116 FORMAT ('SUBROUTINE omega_OF_Y(T, Y, P, KK, WDOT, DTDt)'/&
      'IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)'/&
      'PARAMETER (EXPARG=690.776, SMALL=1.D-300)'/&
      'DIMENSION Y(*), C(KK), G0(KK), WDOT(*), CPML(KK), WT(KK), HML(Kk)'/&
      'DTdt = 0.d0' /&
      'RHO = 0.D0' /&
      '!'/ &
      'T2 = T*T'/&
      'TLOG = LOG(T)'/&
      'TI = 1.D0/T'/&
      'TI2 = TI/T'/&
      'PFAC = 1.218652692702276D-2/T'/ &
      '! G0 is used to calculate equilibrium constants for computing the reverse rate constants'/&
      'CALL SKG0(T, G0)')
117 FORMAT ('SUBROUTINE SKG0(T, G0)'/&
      'IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)'/&
      'DIMENSION G0(*)'/&
      '!'/&
      'T2 = T*T'/&
      'T3 = T2*T'/&
      'T4 = T3*T'/&
      'TLOG = LOG(T)')
118 FORMAT ('SUBROUTINE JAC_OF_C (KK, p, T, C, PD_OF)'/&
      'IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)'/&
      'DIMENSION Y(KK),PD(KK+1,KK+1),CPML(KK),HML(KK),WT(KK),&'/&
        'G0(KK),DG(KK),C(KK),DM(KK),PD_OF(0:(KK+1)*(KK+1))')
        
119 FORMAT (/'SUBROUTINE pTC2rhoY (p, T, C, RHO, Y)'/ &
            'IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)'/ &
            'DIMENSION Y(*), C(*), WT(',I4,')' / &
            'DATA SMALL/1D-50/'/)
120 FORMAT ('SUBROUTINE DMDC_OF_C(A, KK, DM, HML, RHO, CPB)'/&
      'IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)' /&
      'DIMENSION A(KK+1,*), DM(*), HML(*)')

121 FORMAT ('SUBROUTINE omega_OF_C(T, C, P, KK, WDOT, DTDt)'/&
      'IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)'/&
      'PARAMETER (EXPARG=690.776, SMALL=1.D-300)'/&
      'DIMENSION Y(KK), C(KK), G0(KK), WDOT(*), CPML(KK), WT(KK), HML(KK)'/&
      'DTdt = 0.d0' /&
      'RHO = 0.D0' /&
      '!'/ &
      'T2 = T*T'/&
      'TLOG = LOG(T)'/&
      'TI = 1.D0/T'/&
      'TI2 = TI/T'/&
      'PFAC = 1.218652692702276D-2/T'/ &
      '! G0 is used to calculate equilibrium constants for computing the reverse rate constants'/&
      'CALL SKG0(T, G0)')

122 FORMAT (/'SUBROUTINE pTY2rhoC_OF (P, T, Y, RHO, C)'/  &
             'IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)'/ &
             'DIMENSION Y(*), C(*)'/)
             
123 FORMAT (/'SUBROUTINE pTC2rhoY_OF (p, T, C, RHO, Y)'/ &
            'IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)'/ &
            'DIMENSION Y(*), C(*), WT(',I4,')' / &
            'DATA SMALL/1D-50/'/)


10000 FORMAT ('! ',A50)
10001 FORMAT (5X,'C(',I4,') = ','Y(',I4,')/',D20.14)
10002 FORMAT (  /'SUM = 0.D0'/ &
            'DO K = 1, ',I4/  &
            '  SUM = SUM + C(K)'/ &
            'ENDDO'/ &
            'RHO = P/(SUM*T*8.314510D7)'/ &
            / &
            'DO K = 1, ',I4/  &
            '  C(K) = C(K) * RHO'/ &
            'ENDDO'/ )
10003 FORMAT ('DO K = 1, ',I4/  &
            '  C(K) = RHO * C(K)'/ &
            'ENDDO'/)
10004 FORMAT ('IF (T .GT. ',D12.6,') THEN'/ &
     '  SMH(',I4,') = ',SP,D20.12,'       ',SP,D20.12,'*TI &'/ &
     '            ',SP,D20.12,'*TN(1) ',SP,D20.12,'*TN(2) &'/ &
     '            ',SP,D20.12,'*TN(3) ',SP,D20.12,'*TN(4) &'/ & 
     '            ',SP,D20.12,'*TN(5) ')
10005 FORMAT ('ELSE'/ &
     '  SMH(',I4,') = ',SP,D20.12,'       ',SP,D20.12,'*TI &'/ &
     '            ',SP,D20.12,'*TN(1) ',SP,D20.12,'*TN(2) &'/ &
     '            ',SP,D20.12,'*TN(3) ',SP,D20.12,'*TN(4) &'/ & 
     '            ',SP,D20.12,'*TN(5) '/ &
      'ENDIF'/)
10006 FORMAT ('EG(',I4,') = EXP(SMH(',I4,'))')
10007 FORMAT ('PFAC1 = PATM / (RU*T)'/ &
            'PFAC2 = PFAC1*PFAC1'/ &
            'PFAC3 = PFAC2*PFAC1'/)
10008 FORMAT ('RF(',I4,') = EXP(',SP,D20.12,' ',SP,D20.12,'*ALOGT &'/&
                SP,D20.12,'*TI)')
10009 FORMAT ('RB(',I4,') = EXP(',SP,D20.12,' ',SP,D20.12,'*ALOGT &'/&
                SP,D20.12,'*TI)')
10010 FORMAT ('RF(',I4,') = EXP(',SP,D20.12,' ',SP,D20.12,'*TI)'/)
10011 FORMAT ('RB(',I4,') = EXP(',SP,D20.12,' ',SP,D20.12,'*TI)'/)
10012 FORMAT ('RF(',I4,') = EXP(',SP,D20.12,' ',SP,D20.12,'*ALOGT )'/)
10013 FORMAT ('RB(',I4,') = EXP(',SP,D20.12,' ',SP,D20.12,'*ALOGT )'/)
10014 FORMAT ('RF(',I4,') = EXP(',SP,D20.12,' )'/)
10015 FORMAT ('RB(',I4,') = EXP(',SP,D20.12,' )'/)
10016 FORMAT ('RB(',I4,') = RF(',I4,') / MAX(EQK, SMALL)')
10017 FORMAT ('! R',I4,': ',A50)
10018 FORMAT ('EQK = EG(',I4,')/EG(',I4,')')
10019 FORMAT ('EQK = EG(',I4,')/EG(',I4,')/EG(',I4,')*PFAC1**',I4)
10020 FORMAT ('EQK = EG(',I4,')/EG(',I4,')/EG(',I4,')*PFAC1**',I4)
10021 FORMAT ('EQK = EXP(')
10022 FORMAT (SP,I4,'*','SMH(',SS,I4,')')
10023 FORMAT ('RKLOW(',I4,') = EXP(',SP,D20.12,' ',SP,D20.12,'*ALOGT &'/&
                SP,D20.12,'*TI)')
10024 FORMAT ('CTB = CTOT')
10025 FORMAT (SP,D20.12,'*C(',SS,I4,')')
10026 FORMAT ('RF(',I4,') = RF(',I4,')*CTB')
10027 FORMAT ('*C(',I4,')**',I4)
10028 FORMAT ('RB(',I4,') = RB(',I4,')*CTB')
10029 FORMAT ('RF(',I4,') = RF(',I4,')')
10030 FORMAT ('RB(',I4,') = RB(',I4,')')
10031 FORMAT (/'PR = RKLOW(',I4,') * CTB / RF(',I4,')'/&
                'PCOR = PR / (1.0 + PR)'/&
                'PRLOG = LOG10(MAX(PR,SMALL))')
10032 FORMAT ('FCENT = ',D20.12,'*EXP(-T/',D20.12,')+',D20.12,'*EXP(T/',D20.12,')')
10033 FORMAT ('FCENT = ',D20.12,'*EXP(-T/',D20.12,')+',D20.12,'*EXP(T/',D20.12,')+ EXP(-',D20.12,'/T)')
10034 FORMAT ('FCLOG = LOG10(MAX(FCENT,SMALL))'/&
      'XN    = 0.75 - 1.27*FCLOG'/&
      'CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)'/&
      'FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)'/&
      'FC = 10.0**FLOG'/&
      'PCOR = FC * PCOR')
10035 FORMAT ('RF(',I4,') = RF(',I4,') * PCOR')
10036 FORMAT ('RB(',I4,') = RB(',I4,') * PCOR')
10037 FORMAT ('RF(',I4,') = EXP(',SP,D20.12,' ',SP,D20.12,'*ALOGT &'/&
                SP,D20.12,'*TI',SP,D20.12,'/T**(1.0/3.0)',SP,D20.12,'/T**(2.0/3.0))')
10038 FORMAT ('ROP = RF(',I4,')-RB(',I4,')')
10039 FORMAT ('WDOT(',I4,') = WDOT(',I4,') +',SS,I2,'*ROP')
10040 FORMAT ('WDOT(',I4,') = WDOT(',I4,') -',SS,I2,'*ROP')
10041 FORMAT ('IF (T .GT. ',D12.6,') THEN'/ &
     'HML(',I4,') = ',SP,D20.12,'&'/&
      SP,D20.12,'*TN(1) ',SP,D20.12,'*TN(2)&'/ &
      SP,D20.12,'*TN(3) ',SP,D20.12,'*TN(4)&'/ &
      SP,D20.12,'*TN(5)' )
10042 FORMAT ('ELSE'/ &
     'HML(',I4,') = ',SP,D20.12,'&'/&
     SP,D20.12,'*TN(1) ',SP,D20.12,'*TN(2)&'/ &
     SP,D20.12,'*TN(3) ',SP,D20.12,'*TN(4)&'/ &
     SP,D20.12,'*TN(5)'/ &
     'ENDIF')
10043 FORMAT ('IF (T .GT. ',D12.6,') THEN'/ &
     'CPML(',I4,') = ',SP,D20.12,'&'/&
      SP,D20.12,'*TN(1) ',SP,D20.12,'*TN(2)&'/ &
      SP,D20.12,'*TN(3) ',SP,D20.12,'*TN(4)' )
10044 FORMAT ('ELSE'/ &
     'CPML(',I4,') = ',SP,D20.12,'&'/&
     SP,D20.12,'*TN(1) ',SP,D20.12,'*TN(2)&'/ &
     SP,D20.12,'*TN(3) ',SP,D20.12,'*TN(4)' / &
     'ENDIF')
10045 FORMAT ('IF (T .GT. ',D12.6,') THEN'/ &
     'G0(',I4,') = ',SP,D20.12,'*TLOG ',SP,D20.12,'*T ',SP,D20.12,'*T2&'/ &
       SP,D20.12,'*T3 ',SP,D20.12,'*T4 ',SP,D20.12,'/T &'/ &
       SP,D20.12 )
10046 FORMAT ('ELSE'/ &
     'G0(',I4,') = ',SP,D20.12,'*TLOG ',SP,D20.12,'*T ',SP,D20.12,'*T2&'/ &
       SP,D20.12,'*T3 ',SP,D20.12,'*T4 ',SP,D20.12,'/T &'/ &
       SP,D20.12 )
10047 FORMAT ('DG(',I4,') = ',SP,D20.12,'/T ',SP,D20.12,SP,D20.12,'*T&'/ &
       SP,D20.12,'*T2 ',SP,D20.12,'*T3 ',SP,D20.12,'/T2' )
10048 FORMAT ('WT(',I4,') = ',D20.12)
10049 FORMAT ('RB(',I4,') = 0.D0')
10050 FORMAT ('G0_SUM = ')
10051 FORMAT (SP,I2,'*G0(',SS,I4,')')
10052 FORMAT ('DGDT = ')
10053 FORMAT (SP,I2,'*DG(',SS,I4,')')
10054 FORMAT ('RF = EXP(',SP,D20.12,' ',SP,D20.12,'*TLOG &'/&
                SP,D20.12,'*TI)')
10055 FORMAT ('EQINV = EXP(MIN(G0_SUM, EXPARG))'/ &
      'EQINV = EQINV / PFAC ** (',SP,I4,')'/ &
      'DGDT = DGDT',SP,I4,' * TI')
10056 FORMAT ('RFLGDT = ',SP,D20.12,'/T',SP,D20.12,'/T2')
10057 FORMAT ('RFLGDT = ',SP,D20.12,'/T',SP,D20.12,'/T2 &'/ &
                            SP,D20.12,'*T**(-4/3)',SP,D20.12,'*T**(-5/3)')
10058 FORMAT ('RL = EXP(',SP,D20.12,' ',SP,D20.12,'*TLOG &'/&
                SP,D20.12,'*TI)')
10059 FORMAT ('CM = CTOT ')
10060 FORMAT (SP,D20.12,'*C(',SS,I4,')')
10061 FORMAT ('RB = RF*EQINV'/ &
              'RBLGDT = RFLGDT +DGDT'/ &
              'CM = MAX(CM, SMALL)'/ &
              'PR = RL*CM/RF')
10062 FORMAT ('PRLGDT = ',SP,D20.12,'/T ',SP,D20.12,'/T2')
11070 FORMAT ('CM = MAX(CM, SMALL)'/&
'PR = RL*CM/RF')
10063 FORMAT ('PRLGDM = 1.D0/CM'/ &
      'TEMP1 = 1.D0 +PR'/ &
      'PC = PR/TEMP1'/ &
      'PCLGDT = PRLGDT/TEMP1'/ &
      'PCLGDM = PRLGDM/TEMP1'/ &
      'PRLG = LOG(MAX(PR, SMALL))') 
10064 FORMAT ('TEMP1 = ',SP,D20.12,'*EXP(-T/ (',SP,D20.12,') )'/ &
                'TEMP2 = ',SP,D20.12,'*EXP(-T/ (',SP,D20.12,') )'/ &
                'TEMP3 = EXP(',SP,D20.12,'/T)')
                
11064 FORMAT ('TEMP1 = ',SP,D20.12,'*EXP(-T/ (',SP,D20.12,') )'/ &
                'TEMP2 = ',SP,D20.12,'*EXP(-T/ (',SP,D20.12,') )'/ &
                'TEMP3 = EXP(',SP,D20.12,'*1.D10/T)')
                
10065 FORMAT ('FCENT = TEMP1 +TEMP2 +TEMP3'/ &
              'FCNTDT = ',SP,D20.12,'*TEMP1 ',SP,D20.12,'*TEMP2',SP,D20.12,'/T2*TEMP3'/ &
              'FTLGDT = FCNTDT/FCENT'/ &
              'FTLG = LOG(MAX(FCENT, SMALL))')
              
11165 FORMAT ('FCENT = TEMP1 +TEMP2 +TEMP3'/ &
              'FCNTDT = ',SP,D20.12,'*TEMP1 ',SP,D20.12,'*TEMP2',SP,D20.12,'*1.D10/T2*TEMP3'/ &
              'FTLGDT = FCNTDT/FCENT'/ &
              'FTLG = LOG(MAX(FCENT, SMALL))')
              
!10066 FORMAT ('TEMP1 = -0.4 +PRLG -6.7D-1*FTLG '/ &
!                'TEMP2 = -0.19 -1.4D-1*PRLG -1.1762D0*FTLG')
10066 FORMAT ('TEMP1 = -9.210340371976184D-1 +PRLG -6.7D-1*FTLG '/ &
                'TEMP2 = 1.855883584953201D0 -1.4D-1*PRLG -1.1762D0*FTLG'/ &
                'XP = TEMP1/TEMP2'/ &
                'TEMP3 = 1.D0 +1.4D-1*XP'/ &
                'XPDT = (TEMP3*PRLGDT +(1.1762D0*XP -6.7D-1)*FTLGDT)/TEMP2'/ &
                'XPDM = TEMP3*PRLGDM/TEMP2'/ &
                'TEMP1 = 1.D0 +XP*XP'/ &
                'FCLG = FTLG/TEMP1'/ &
                'FC = EXP(FCLG)'/ &
                'TEMP2 = FCLG*(XP +XP)'/ &
                'FCLGDT = (FTLGDT -TEMP2*XPDT)/TEMP1'/ &
                'FCLGDM = -TEMP2*XPDM/TEMP1'/ &
                'PC = PC*FC'/ &
                'PCLGDT = PCLGDT +FCLGDT'/ &
                'PCLGDM = PCLGDM +FCLGDM'/ &
                'RF = RF*PC'/ &
                'RB = RB*PC') 
10067 FORMAT ('RB = EXP(',SP,D20.12,' ',SP,D20.12,'*TLOG &'/&
                SP,D20.12,'*TI)')
10068 FORMAT ('RBLGDT = ',SP,D20.12,'/T',SP,D20.12,'/T2')
10069 FORMAT ('RB = 0.D0')
10070 FORMAT ('RBLGDT = 0.D0')
10071 FORMAT ('RB = RF*EQINV'/ &
              'RBLGDT = RFLGDT +DGDT'/ &
              'CM = MAX(CM, SMALL)'/ &
                 'PC = CM'/ &
                 'PCLGDM = 1.D0/CM'/ &
                 'RF = RF*PC'/ &
                 'RB = RB*PC')
10072 FORMAT ('CM = MAX(CM, SMALL)'/ &
                 'PC = CM'/ &
                 'PCLGDM = 1.D0/CM'/ &
                 'RF = RF*PC')
10073 FORMAT ('RB = RF*EQINV'/ &
              'RBLGDT = RFLGDT +DGDT')
10074 FORMAT ('RF = EXP(',SP,D20.12,' ',SP,D20.12,'*ALOGT &'/&
                SP,D20.12,'*TI',SP,D20.12,'/T**(1.0/3.0)',SP,D20.12,'/T**(2.0/3.0))')
10075 FORMAT ('WF = RF')
10076 FORMAT ('*C(',I4,')**',I4)
10077 FORMAT ('WFDT = WF*(RFLGDT +PCLGDT)'/&
                'WFDM = WF*PCLGDM')
10078 FORMAT ('WB = RB')
10079 FORMAT ('WBDT = WB*(RBLGDT +PCLGDT)'/&
            'WBDM = WB*PCLGDM'/&
            'DWDT = WFDT -WBDT')
10080 FORMAT ('WFDT = WF*RFLGDT'/&
             'WFDM = WF*PCLGDM')
10081 FORMAT ('WBDT = WB*RBLGDT'/&
             'WBDM = WB*PCLGDM'/&
             'DWDT = WFDT -WBDT')
10082 FORMAT ('WFDT = WF*RFLGDT')
10083 FORMAT ('WBDT = WB*RBLGDT'/&
             'DWDT = WFDT -WBDT')
10084 FORMAT ('A(',I4,',',I4,') = A(',I4',',I4,')',SP,I4,'*DWDT')
10085 FORMAT ('WFDC = RF *')
10086 FORMAT ('WBDC = RB *')
10087 FORMAT (I4,'*C(',I4,')**',I4)
10089 FORMAT ('A(',I4,',',I4,') = A(',I4',',I4,')',SP,I4,'* WFDC')
10090 FORMAT ('A(',I4,',',I4,') = A(',I4',',I4,')',SP,I4,'* WBDC')
10091 FORMAT ('DWDM = WFDM -WBDM')
10092 FORMAT ('A(',I4,',',I4,') = A(',I4',',I4,')',SP,D20.12,' * (',SP,I4,') * DWDM')
10093 FORMAT ('DM(',I4,') = DM(',I4,') ',SP,I4,'* DWDM')
10094 FORMAT ('ENDIF')
10095 FORMAT ('DO N=2,',I4)
10096 FORMAT ('A(',I4,', N) = A(',I4,', N) + DM(',I4,')')
!10096 FORMAT ('A(',I4,', N) = A(',I4,', N) + DM(',I4,') * WT(',I4,') / WT(N-1)')
!10096 FORMAT ('A(',I4,', N) = A(',I4,', N) + DM(',I4,')/ WT(',I4,') * RHO')
10097 FORMAT ('IF (KK+1 .LT. ',I4,') THEN'/&
'WRITE(*,*) "DA TOO SMALL"' /&
'STOP'/&
'ENDIF')
10098 FORMAT ('DO M=1,',I4/&
'DO N=1,',I4/&
'  A(N,M) = 0.D0'/&
'ENDDO'/&
'ENDDO'/&
'!')
10099 FORMAT ('DO N=1,',I4/&
'  G0(N) = 0.D0'/&
'  DG(N) = 0.D0'/&
'  DM(N) = 0.D0'/&
'ENDDO'/&
'!')
10100 FORMAT ('CTOT = 0.D0'/&
'DO N=1,',I4/&
'  CTOT = CTOT +C(N)'/&
'ENDDO'/&
'!'/&
'CALL GDG1(T, G0, DG)'/&
'!'/&
'CALL CKRHOY (RPAR(NP), Z(1), Z(2), IPAR, RPAR, RHO)'/&
'!'/&
'CALL CKCPBS (Z(1), Z(2), IPAR, RPAR, CPB)'/&
'!'/&
'CALL SKWT (WT)'/&
'!'/&
'CALL AJ1(T, C, A, KK, G0, DG, DM, CTOT, WT, RHO, CPB, CPML, HML)'/&
'!'/&
'CALL DMDC(A, KK, DM)')
10101 FORMAT ('DIMENSION G0(',I4,'), DG(',I4,'), DM(',I4,'), C(*), A(KK+1,*)')
10102 FORMAT ('ENDDO')
10104 FORMAT ('!'/&
'DO K=1,',I4/&
'  DTDT(K) = 0.D0'/&
'  DTDC(K) = 0.D0' /&
'  WK(K) = 0.D0' /&
'  DM(K) = 0.D0' /&
'ENDDO' /&
'!')
10105 FORMAT ('DTDT(',SS,I4,') = DTDT(',SS,I4,') ',SP,I4,' *DWDT')
10106 FORMAT ('A(',I4,',',I4,') = A(',I4',',I4,') &')
10107 FORMAT (' + DTDT(',I4,') * HML(',I4,') ',' + CPML(',I4,') * WK(',I4,') &')
10108 FORMAT ('CALL SKCPML  (T, CPML)'/&
            'CALL SKHML  (T, HML)')
10110 FORMAT ('WK(',SS,I4,') = WK(',SS,I4,') ',SP,I4,' *(WF - WB)')
10111 FORMAT (' + 0.D0')
10113 FORMAT ('DTDC(',SS,I4,') = DTDC(',SS,I4,') ',SP,I4,' *WFDC')
10114 FORMAT ('DTDC(',SS,I4,') = DTDC(',SS,I4,') ',SP,I4,' *(-1)*WBDC')
!10115 FORMAT (' + DTDC(',I4,') * HML(',I4,') * WT(',I4,')')
10115 FORMAT ('A(',I4,',',I4,') = A(',I4',',I4,')  + DTDC(',I4,') * HML(',I4,') / CPB')
!10115 FORMAT ('A(',I4,',',I4,') = A(',I4',',I4,')  + DTDC(',I4,') * HML(',I4,') * WT(',I4,') / CPB')
10116 FORMAT ('A(',I4,',',I4,') = A(',I4,',',I4,') /(RHO * CPB)')
10117 FORMAT ('DO K=2,',I4 /&
                '  DO K1=1,',I4 /&
                '    A(K,K1+1) = A(K,K1+1) / WT(K1) * RHO'/&
                '  ENDDO'/&
                'ENDDO')
!10118 FORMAT ('DO K=1,',I4 /&
!                'A(',I4,',K) = A(',I4,',K) / (RHO * CPB)'/&
!                'ENDDO')
10118 FORMAT ('A(',I4,',',I4,') = A(',I4,',',I4,') / (RHO * CPB)')
10120 FORMAT ('CALL pTY2rhoC (RPAR(NP), Z(1), Z(2), RHO, C)' / &
'CTOT = 0.D0'/&
'DO N=1,',I4/&
'  CTOT = CTOT +C(N)'/&
'ENDDO'/&
'!'/&
'CALL GDG1(Z(1), G0, DG)'/&
'!'/&
'CALL CKRHOY (RPAR(NP), Z(1), Z(2), IPAR, RPAR, RHO)'/&
'!'/&
'CALL CKCPBS (Z(1), Z(2), IPAR, RPAR, CPB)'/&
'!'/&
'CALL SKWT (WT)'/&
'!'/&
'CALL SKCPML  (Z(1), CPML)'/&
'!'/&
'CALL SKHML  (Z(1), HML)'/&
'!'/&
'CALL AJ1(Z(1), C, PD, KK, G0, DG, DM, CTOT, WT, RHO, CPB, CPML, HML)'/&
'!'/&
'CALL DMDC(PD, KK, DM, HML, WT, RHO, CPB)')
!10121 FORMAT ('A(',I4,',',I4,') = A(',I4',',I4,')',SP,D20.12,'* DWDM * HML(',SS,I4,') * WT(',SS,I4,')')
!10122 FORMAT ('DO K=1,',I4 /&
!                'A(',I4,',K+1) = A(',I4,',K+1) + DM(K) * HML(K) * WT(K) / WT(K) * RHO / (RHO * CPB)'/&
!                'ENDDO')
10122 FORMAT ('DO K=1,',I4 /&
                '  A(',I4,',K+1) = A(',I4,',K+1) + DM(K) * HML(K) / CPB'/&
                'ENDDO')
10123 FORMAT ('DO K=1,',I4 /&
                '  A(1,K+1) = A(1,K+1) / WT(K) * RHO'/&
                'ENDDO')
10124 FORMAT (') / PFAC1 ** (',SP,I4,')')
10125 FORMAT ('FC = 1.D0'/&
      'PCOR = FC * PCOR')
10126 FORMAT ('PRLGDM = 1.D0/CM'/ &
      'TEMP1 = 1.D0 +PR'/ &
      'PC = PR/TEMP1'/ &
      'PCLGDT = PRLGDT/TEMP1'/ &
      'PCLGDM = PRLGDM/TEMP1'/&
      'FC = 1.D0'/ &
      'FCLGDT = 0.D0'/ &
      'FCLGDM = 0.D0'/ &
      'PC = PC*FC'/ &
      'PCLGDT = PCLGDT +FCLGDT'/ &
      'PCLGDM = PCLGDM +FCLGDM'/ &
      'RF = RF*PC'/ &
      'RB = RB*PC') 
10127 FORMAT (' &')
10128 FORMAT ('DO K=1,',I4 /&
                'DO K1=1,',I4 /&
                '  A(K+1,K1+1) = A(K+1,K1+1) * WT(K) / WT(K1)'/&
                'ENDDO'/&
                'ENDDO')
10129 FORMAT ('DO K=1,',I4 /&
                '  DO K1=1,',I4 /&
!                'A(1,K+1) = A(1,K+1) + A(K1+1,K+1) * HML(K1) * WT(K1) / (CPB * WT(K))'/&
                '    A(1,K+1) = A(1,K+1) + A(K1+1,K+1) * HML(K1) / (CPB * WT(K))'/&
                '  ENDDO'/&
                'ENDDO')
10130 FORMAT ('DO K=1,',I4 /&
                '  A(1,K+1) = A(1,K+1) / CPB / WT(K)'/&
                'ENDDO')
10131 FORMAT ('DO K=1,',I4 /&
                '  A(k+1,1) = A(k+1,1) * WT(K) / RHO'/&
                'ENDDO')
10132 FORMAT ('RB(',I4,') = 0.d0')
10133 FORMAT ('RF(',I4,') = 0.d0')
10134 FORMAT ('RFLGDT = 0.d0')
!10135 FORMAT ('RBLGDT = 0.d0')
10136 FORMAT ('CM = MAX(CM, SMALL)'/ &
                 'PC = CM'/ &
                 'PCLGDM = 1.D0/CM'/ &
                 'RF = RF*PC'/ &
                 'RB = RB*PC')


10137 FORMAT ('!'/&
'DO k1 = 0, KK'/&
'  DO k2 = 0, KK'/&
'    PD(k1+1,k2+1) = 0.d0'/&
'  ENDDO'/&
'ENDDO'/&
'DO N=1, KK'/&
'  CPML(N) = 0.D0'/&
'  HML(N) = 0.D0'/&
'  WT(N) = 0.D0'/&
'  G0(N) = 0.D0'/&
'  DG(N) = 0.D0'/&
'  C(N) = 0.D0'/&
'  DM(N) = 0.D0'/&
'ENDDO'/&
'!'/&
'CALL pTY2rhoC_OF (p, T, Y(1), RHO, C)' / &
'CTOT = 0.D0'/&
'DO N=1, KK'/&
'  CTOT = CTOT +C(N)'/&
'ENDDO'/&
'! WT (g/mol)'/&
'CALL SKWT (WT)'/&
'! standard-state Gibbs free energy G0: ergs/mole, DG: dG0/dT'/&
'! G0 is used to calculate equilibrium constants for computing the reverse rate constants'/&
'CALL GDG1(T, G0, DG)'/&
'! CPML: ergs/(mole*K)'/&
'CALL SKCPML  (T, CPML)'/&
'! HML: ergs/mole'/&
'CALL SKHML  (T, HML)'/&
'! CPB: ergs/(g*K)'/&
'CPB = 0.D0' /&
'DO K=1, KK' /&
'  CPB = CPB + CPML(K) / WT(K) * Y(K)'/&
'ENDDO' /&
'!'/&
'CALL AJ1(T, C, PD, KK, G0, DG, DM, CTOT, WT, RHO, CPB, CPML, HML)'/&
'! DM is related to three body reactions'/&
'CALL DMDC(PD, KK, DM, HML, WT, RHO, CPB)'/)

10155 FORMAT ('EQINV = EXP(MIN(G0_SUM, EXPARG))'/ &
      'EQINV = EQINV / PFAC ** (',SP,I4,')')
10161 FORMAT ('RB = RF*EQINV'/ &
              'CM = MAX(CM, SMALL)'/ &
              'PR = RL*CM/RF')
11126 FORMAT ('TEMP1 = 1.D0 +PR'/ &
      'PC = PR/TEMP1'/ &
      'FC = 1.D0'/ &
      'PC = PC*FC'/ &
      'RF = RF*PC'/ &
      'RB = RB*PC') 
11063 FORMAT ('TEMP1 = 1.D0 +PR'/ &
      'PC = PR/TEMP1'/ &
      'PRLG = LOG(MAX(PR, SMALL))') 
11136 FORMAT ('CM = MAX(CM, SMALL)'/ &
                 'PC = CM'/ &
                 'RF = RF*PC'/ &
                 'RB = RB*PC')
11073 FORMAT ('RB = RF*EQINV')
11072 FORMAT ('CM = MAX(CM, SMALL)'/ &
                 'PC = CM'/ &
                 'RF = RF*PC')
11071 FORMAT ('RB = RF*EQINV'/ &
              'CM = MAX(CM, SMALL)'/ &
                 'PC = CM'/ &
                 'RF = RF*PC'/ &
                 'RB = RB*PC')
11133 FORMAT ('RF = 0.D0')

11065 FORMAT ('FCENT = TEMP1 +TEMP2 +TEMP3'/ &
              'FTLG = LOG(MAX(FCENT, SMALL))')
11066 FORMAT ('TEMP1 = -9.210340371976184D-1 +PRLG -6.7D-1*FTLG '/ &
                'TEMP2 = 1.855883584953201D0 -1.4D-1*PRLG -1.1762D0*FTLG'/ &
                'XP = TEMP1/TEMP2'/ &
                'TEMP1 = 1.D0 +XP*XP'/ &
                'FCLG = FTLG/TEMP1'/ &
                'FC = EXP(FCLG)'/ &
                'PC = PC*FC'/ &
                'RF = RF*PC'/ &
                'RB = RB*PC') 
11104 FORMAT ('CALL pTY2rhoC_OF (p, T, Y(1), RHO, C)'/&
'CTOT = 0.D0'/&
'DO K=1,',I4/&
'  WDOT(K) = 0.D0 ! ergs/(cm3*s)' /&
'  CTOT = CTOT +C(K)'/&
'ENDDO' /&
'!')
11129 FORMAT ('DO K=1,',I4 /&
                '  DO K1=1,',I4 /&
                '    A(1,K+1) = A(1,K+1) + A(K1+1,K+1) * HML(K1) / (CPB * RHO)'/&
                '  ENDDO'/&
                'ENDDO')
11114 FORMAT ('CALL pTC2rhoY_OF (p, T, C, RHO, Y)'/&
'CTOT = 0.D0'/&
'DO K=1,',I4/&
'  WDOT(K) = 0.D0 ! ergs/(cm3*s)' /&
'  CTOT = CTOT +C(K)'/&
'ENDDO' /&
'!')


20002 FORMAT ('CALL SKWT (WT)' /&
'WM1 = 0.D0' /&
'WM2 = 0.D0' /&
'DO K=1,',I4 /&
'  WM1 = WM1 + C(K) '/ &
'  WM2 = WM2 + C(K)*WT(K) '/ &
'ENDDO'/&
'WM = WM2/max(WM1,1.D-50) '/ &
'RHO = P*WM/(T*8.314510D7)' / &
'!' / &
'DO K=1,',I4 /&
'  Y(K) = WT(K) * C(K) / RHO'/ &
'ENDDO'/&
)

20039 FORMAT ('WDOT(',I4,') = WDOT(',I4,') +',SS,I2,'*(WF-WB)')
20040 FORMAT ('WDOT(',I4,') = WDOT(',I4,') -',SS,I2,'*(WF-WB)')
20041 FORMAT ('DO k1 = 0, KK '/ &
                '  DO k2 = 0, KK'/ &
                '    PD_OF(k1*(KK+1)+k2) = PD(k1+1,k2+1)'/ &
                '  ENDDO'/ &
                'ENDDO') 
20042 FORMAT ('! WT (g/mol)'/&
'CALL SKWT (WT)'/&
'! CPML: ergs/(mole*K)'/&
'CALL SKCPML  (T, CPML)'/&
'! HML: ergs/mole'/&
'CALL SKHML  (T, HML)'/ &
'! CPB: ergs/(g*K)'/&
'CPB = 0.D0' /&
'DTDt = 0.d0' /&
'DO K=1,',I4 /&
'  CPB = CPB + CPML(K) / WT(K) * Y(K)'/&
'  DTDt = DTDt + HML(k)*WDOT(k)'/&
'ENDDO' /&
'DTDt = -DTDt / (RHO*CPB)')       

20137 FORMAT ('!' /&
'DO k1 = 0, KK'/&
'  DO k2 = 0, KK'/&
'    PD(k1+1,k2+1) = 0.d0'/&
'  ENDDO'/&
'ENDDO'/&
'!' /&
'DO N=1, KK'/&
'  CPML(N) = 0.D0'/&
'  HML(N) = 0.D0'/&
'  WT(N) = 0.D0'/&
'  G0(N) = 0.D0'/&
'  DG(N) = 0.D0'/&
'  Y(N) = 0.D0'/&
'  DM(N) = 0.D0'/&
'ENDDO'/&
'!'/&
'CALL pTC2rhoY_OF (p, T, C, RHO, Y)' / &
'CTOT = 0.D0'/&
'DO N=1, KK'/&
'  CTOT = CTOT +C(N)'/&
'ENDDO'/&
'! WT (g/mol)'/&
'CALL SKWT (WT)'/&
'! standard-state Gibbs free energy G0: ergs/mole, DG: dG0/dT'/&
'! G0 is used to calculate equilibrium constants for computing the reverse rate constants'/&
'CALL GDG1(T, G0, DG)'/&
'! CPML: ergs/(mole*K)'/&
'CALL SKCPML  (T, CPML)'/&
'! HML: ergs/mole'/&
'CALL SKHML  (T, HML)'/&
'! CPB: ergs/(g*K)'/&
'CPB = 0.D0' /&
'DO K=1, KK'/&
'  CPB = CPB + CPML(K) / WT(K) * Y(K)'/&
'ENDDO' /&
'!'/&
'CALL AJ1(T, C, PD, KK, G0, DG, DM, CTOT, WT, RHO, CPB, CPML, HML)'/&
'! DM is related to three body reactions'/&
'CALL DMDC_OF_C(PD, KK, DM, HML, RHO, CPB)'/)
                
30002 FORMAT (  /'SUM = 0.D0'/ &
            'DO K = 1, ',I4/  &
            '  SUM = SUM + C(K)'/ &
            'ENDDO'/ &
            '! 1 N/m2  = 1 Pa (OF_P) = 10 dynes/cm2 (CK_P)  ' / &
            'RHO = P*10/(SUM*T*8.314510D7)'/ &
            / &
            'DO K = 1, ',I4/  &
            '  C(K) = C(K) * RHO'/ &
            'ENDDO'/ )           
30003 FORMAT ('CALL SKWT (WT)' /&
'WM1 = 0.D0' /&
'WM2 = 0.D0' /&
'DO K=1,',I4 /&
'  WM1 = WM1 + C(K)  '/ &
'  WM2 = WM2 + C(K)*WT(K) '/ &
'ENDDO'/&
'WM = WM2/max(WM1,1.D-50) '/ &
'! 1 N/m2  = 1 Pa (OF_P) = 10 dynes/cm2 (CK_P)  ' / &
'RHO = P*10*WM/(T*8.314510D7)' / &
'!' / &
'DO K=1,',I4 /&
'  Y(K) = WT(K) * C(K) / RHO'/ &
'ENDDO'/&
)        


     STOP
END PROGRAM runexample1



SUBROUTINE delSameNumbANDsort(NSP,NSPinREV,NNSP,NNSPinREV)
     IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER(I-N)
     INTEGER , INTENT (INOUT) ::NSP,NSPinREV
     INTEGER , INTENT (OUT) ::NNSP,NNSPinREV
     DIMENSION NSPinREV(NSP),NNSPinREV(NSP)
     
     do N=1,NSP-1
       DO NN=N+1,NSP
        IF (NSPinREV(N) == NSPinREV(NN) .AND. NSPinREV(N)>0)THEN
          NSPinREV(NN) = -1
        ENDIF
       ENDDO
     enddo
     
     NNSP = 0
     DO N=1,NSP
        IF(NSPinREV(N) .NE. -1)THEN
            NNSP = NNSP+1
            NNSPinREV(NNSP) = NSPinREV(N)
        ENDIF
     ENDDO
     
    do N=1,NNSP-1
      Nmin = NNSPinREV(N)
       DO NN=N+1,NNSP
        IF (Nmin > NNSPinREV(NN))THEN
          Nmin = NNSPinREV(NN)
          Ntemp = NNSPinREV(NN)
          NNSPinREV(NN) = NNSPinREV(N)
          NNSPinREV(N) = Ntemp
        ENDIF
       ENDDO
     enddo

end
