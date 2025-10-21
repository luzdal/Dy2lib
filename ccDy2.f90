    !-----------------------------------------------------------------------------------------
    SUBROUTINE DRIVE_CC_SINGLE(E,scatlen)
    USE SYSTEM_PARAMETER,ONLY: RVDW
    IMPLICIT NONE
    REAL*8,PARAMETER :: R_MIN = 0.01D0
    REAL*8,PARAMETER :: R_MID =  1.D0
    REAL*8,PARAMETER :: R_MIDD = 1.D1
    REAL*8,PARAMETER :: R_MAX = 1.D2
    REAL*8,PARAMETER :: R_MAXX =1.D3
    INTEGER,PARAMETER :: N_STEP_S = 40000
    INTEGER,PARAMETER :: N_STEP_L = 40000
    INTEGER,PARAMETER :: N_STEP_LL =40000
    INTEGER,PARAMETER :: N_STEP_LLL=40000
    INTEGER,PARAMETER :: N_STEP=N_STEP_S+N_STEP_L+N_STEP_LL+N_STEP_LLL
    REAL*8,PARAMETER :: DR_S = (R_MID-R_MIN)/N_STEP_S
    REAL*8,PARAMETER :: DR_L = (R_MIDD-R_MID)/N_STEP_L
    REAL*8,PARAMETER :: DR_LL = (R_MAX-R_MIDD)/N_STEP_LL
    REAL*8,PARAMETER :: DR_LLL=(R_MAXX-R_MAX)/N_STEP_LLL
    REAL*8 Y14,Y23,Y,PP,PPH,R,H,W_REF
    REAL*8,DIMENSION(1):: WVEC
    REAL*8,DIMENSION(1,1) :: KM,YY
    INTEGER,DIMENSION(1) :: L
    REAL*8 E,PHASE,scatlen
    INTEGER I,J
    REAL*8,EXTERNAL :: POTENTIAL
    !POTENTIAL_ID=1
    WVEC=DSQRT(DABS(E))
    L=0
    Y=0.D0
    Y14=0.D0
    Y23=0.D0
    R=R_MIN
    !W_REF=POTENTIAL_4HE(R)-E
    W_REF=POTENTIAL(R)-E
    Y=DSQRT(DABS(W_REF))
    
    DO I=1,N_STEP
    !=====================================
        IF(I.LE.N_STEP_S) THEN
            H=DR_S
        ELSE IF(I.LE.(N_STEP_S+N_STEP_L))THEN
            H=DR_L
        ELSE IF(I.LE.(N_STEP_S+N_STEP_L+N_STEP_LL)) THEN
            H=DR_LL
        ELSE
            H=DR_LLL
        END IF
    !=====================================
        R=R+H
        !W_REF=POTENTIAL_4HE(R)-E
        W_REF=POTENTIAL(R)-E
        PP=DSQRT(DABS(W_REF))
        PPH=PP*H
        IF (W_REF.LT.0.D0) THEN
               Y14=PP/DTAN(PPH)
               Y23=PP/DSIN(PPH)
           ELSE IF (W_REF.GT.0.D0)THEN
               
               Y14=PP/DTANH(PPH)
               Y23=PP/DSINH(PPH)
           ELSE
               Y14 = 1.D0/H
               Y23 = 1.D0/H
           END IF
        Y=Y14-Y23*Y23/(Y+Y14)
    END DO
    YY=Y
    CALL YTOK(WVEC,L,1,1,YY,KM,R)
    PHASE=DATAN(KM(1,1))
    scatlen=-KM(1,1)/DSQRT(E)*rvdw
    PRINT*,'scatering length in a0=',scatlen
    END
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    !Some Arithmetic
    !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    !*************************************************************

      SUBROUTINE YTOK(WVEC,L,N,NOPEN,Y,YK,RUP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
!     ROUTINE TO OBTAIN THE K MATRIX FROM 
!     THE LOG DERIVATIVE MATRIX
!     ON ENTRY,  Y HOLDS THE LOG DERIVATIVE MATRIX (N,N)
!     ON EXIT,   YK HOLDS THE AUGMENTED   K MATRIX (NOPEN,NOPEN)
!
!     THESE ARE INPUTS
      DIMENSION WVEC(N), L(N), IS(N), JS(N)     
      DIMENSION Y(N,N),YK(NOPEN,NOPEN),YKK(N,N)
!
!     WVEC(I)=DQSRT(ABS(ETOT-EINT(I))) 
!     L(I)= ORBITAL QUANTUM NUMBER OF CHANNEL I
!
!     THESE ARE INTERNAL VARIABLES
      DIMENSION ABES(N,N), BBES(N,N),  &
          SJ(N), SJP(N), SN(N), SNP(N)

!     NOTE: Y IS A TWO DIMENSIONAL ARRAY IN DAPROP
      IF(NOPEN.EQ.0) RETURN
      DO I = 1,NOPEN
         DW = WVEC(I)
         DARG = DW*RUP
         CALL BESOPEN(DARG, L(I), RJ, RN, RJP, RNP)
         ROOTDW = SQRT(DW)
         SJ(I) = RJ/ROOTDW
         SJP(I) = RJP*ROOTDW
         SN(I) = RN/ROOTDW
         SNP(I) = RNP*ROOTDW
      ENDDO 
! 
!
      NCLOSE = N - NOPEN
      IF(NCLOSE.GT.0) THEN
       DO I = (NOPEN+1),N
         DW = WVEC(I)
         DARG = DW*RUP
         CALL BESCLOSED(DARG, L(I), RJ, RN, RJP, RNP)
         ROOTDW = SQRT(DW)
         SJ(I) = RJ/ROOTDW
         SJP(I) = RJP*ROOTDW
         SN(I) = RN/ROOTDW
         SNP(I) = RNP*ROOTDW
       ENDDO 
      ENDIF
!
!     FILL BBES=J'-Y*J
!
      DO J=1,N
        DO I=1,N
          BBES(I,J)=-Y(I,J)*SJ(J)
        ENDDO
        BBES(J,J)=BBES(J,J)+SJP(J)
      ENDDO  
!
!     FILL ABES=N'-Y*N
!
      DO J=1,N
        DO I=1,N
          ABES(I,J)=-Y(I,J)*SN(J)
        ENDDO
        ABES(J,J)=ABES(J,J)+SNP(J)
      ENDDO 
!
!     SOLVE ABES*K=BBES TO GET THE AUGMENTED K-MATRIX
!     ON EXIT BBES CONTAIN THE SOLUTION MATRIX
!
      !CALL DGESV(N,N,ABES,N,SNP,BBES,N,IERR)
		CALL BRINV(ABES,N)
		YKK=MATMUL(ABES,BBES)	
      DO J=1,NOPEN
       DO I=1,NOPEN
         YK(I,J)=YKK(I,J) !BBES(I,J)
       ENDDO
      ENDDO
      CALL KSYM(YK,NOPEN)
      RETURN    
 901  FORMAT('0***** ERROR IN LINEAR EQUATION SOLVER IN YTOK.',  &
           '  IER =',I4,'.  RUN HALTED.')
      END  

	       Subroutine besopen(x,n,rj,ry,rjp,ryp)     
!
!-------------------------------------------------------------
!     
!     THIS SUBROUTINES COMPUTE ODINARY AND MODIFIED 
!     RICCATI-BESSEL AND RICCATI-NEUMANN 
!     FUNCTIONS AND THEIR DERIVATIVES.
! 
!     - OPEN CHANNELS:
!       RICCATI-BESSEL (RJ) AND RICCATI-NEUMANN (RY)
!       ARE DEFINED ACCORDINGLY TO ABRAMOWITZ&STEGUN
!            F_N(Z)=Z*FSPH_N(Z) 
!       WHERE 
!            FSPH_N(Z)=DSQRT(PI/2/Z)*F_(N+0.5)(Z)
!       ARE THE USUAL SPHERICAL FUNCTIONS. 
!       THE WRONSKIAN W{RJ,RY} IS 1.
!
!     - CLOSED CHANNELS:
!       SCALED RICCATI-BESSEL (RI) AND RICCATI-NEUMANN (RK)
!       ARE DEFINED AS
!       RI_N(Z) =  DSQRT(Z)*I_(N+0.5) (Z) * EXP(-Z)
!       RK_N(Z) = -DSQRT(Z)*K_(N+0.5) (Z) * EXP(Z)
!       RIP_N(Z)= (RI_N(Z) * EXP(Z))' * EXP(-Z)
!       RKP_N(Z)= (RK_N(Z) * EXP(-Z))' * EXP(Z)
!       THE WRONSKIAN W{RI*EXP(Z),RK*EXP(-Z)} 
!       (NON-SCALED QUANTITIES) IS 1.
!       THE SCALING AVOIDS NUMERICAL UNDERFLOWS AND OVERFLOWS
!
!-----------------------------------------------------------------  
      implicit double precision(a-h,o-z)
      PARAMETER (XMIN=2.,PI=3.141592653589793d0)
!     
!     PI is used for consistency with bessjy
!    
      xnu=dble(n)+0.50d0
      dsx=dsqrt(x)
      factor=dsqrt(pi/2.0d0)
      call bessjy(x,xnu,rj,ry,rjp,ryp)
      rjp=factor*(1.0d0/2.0d0/dsx*rj+dsx*rjp)
      ryp=factor*(1.0d0/2.0d0/dsx*ry+dsx*ryp)
      rj=factor*dsx*rj
      ry=factor*dsx*ry 
      return
      end 

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      Subroutine besclosed(x,n,ri,rk,rip,rkp)
      implicit double precision(a-h,o-z)
      PARAMETER (XMIN=2.,PI=3.141592653589793d0)
!
!     change XMIN accordingly to bessik
!     PI is used for consistency with bessik
!
      xnu=dble(n)+0.50d0
      dsx=dsqrt(x)
      call bessik(x,xnu,ri,rk,rip,rkp)
      rip=(1.0d0/2.0d0/dsx*ri+dsx*rip)
      rkp=-(1.0d0/2.0d0/dsx*rk+dsx*rkp)
      ri=dsx*ri
      rk=-dsx*rk
!
!     bessik computes scaled quantities only for x.ge.XMIN
!
      if(x.lt.XMIN) then
       scale=dexp(x)
       ri=ri/scale
       rip=rip/scale
       rk=rk*scale
       rkp=rkp*scale
      endif
      return
      end


	  SUBROUTINE bessik(x,xnu,ri,rk,rip,rkp) 
      INTEGER MAXIT 
      DOUBLE PRECISION ri,rip,rk,rkp,x,xnu,XMIN 
      DOUBLE PRECISION EPS,FPMIN,PI 
      PARAMETER (EPS=1.e-16,FPMIN=1.e-30,MAXIT=1000000,XMIN=2., &
           PI=3.141592653589793d0) 
      
!     USES beschb
!     Returns the modified Bessel functions ri=I_v, rk=K_v and their derivatives
!     rip=I'_v, rkp=K'_v for positive x and for xnu=v.ge.0. The relative accuracy
!     is whithin one or two significant digits of EPS. FPMIN is a number close to 
!     the machine's smallest floating point number. All internal arithmetic is in 
!     double precision. To convert the entire routine to double precision, change
!     REAL declaration above and decrease EPS to 10**-16.
!     Also convert the subroutine beschb. 
      
      INTEGER i,l,nl 
      DOUBLE PRECISION a,a1,b,c,d,del,del1,delh,dels,e,f,fact,  &
           fact2,ff,gam1,gam2,gammi,gampl,h,p,pimu,q,q1,q2,     &
           qnew,ril,ril1,rimu,rip1,ripl,ritemp,rk1,rkmu,rkmup,  &
           rktemp,s,sum,sum1,x2,xi,xi2,xmu,xmu2
      
      if(x.le.0..or.xnu.lt.0.) pause 'bad arguments in bessik'
      nl=int(xnu+.5d0)
      xmu=xnu-nl 
      xmu2=xmu*xmu 
      xi=1.d0/x 
      xi2=2.d0*xi 
      h=xnu*xi 
      if(h.lt.FPMIN)h=FPMIN 
      b=xi2*xnu 
      d=0.d0 
      c=h 
      do i=1,MAXIT
         b=b+xi2 
         d=1.d0/(b+d) 
         c=b+1.d0/c 
         del=c*d 
         h=del*h
         if(abs(del-1.d0).lt.EPS)goto 1 
      enddo  

      pause 'x too large in bessik;try asymptotic expansion' 
 1    continue
      ril=FPMIN 
      ripl=h*ril
      ril1=ril
      rip1=ripl 
      fact=xnu*xi 
      do l=nl,1,-1
         ritemp=fact*ril+ripl 
         fact=fact-xi 
         ripl=fact*ritemp+ril 
         ril=ritemp
      enddo  
      f=ripl/ril  
      if(x.lt.XMIN)then 
         x2=.5d0*x 
         pimu=PI*xmu 
         if(abs(pimu).lt.EPS)then
            fact=1.d0 
         else
            fact=pimu/sin(pimu) 
         endif 
         d=-log(x2) 
         e=xmu*d 
         if(abs(e).lt.EPS)then
            fact2=1.d0 
         else
            fact2=sinh(e)/e 
         endif 
         call beschb(xmu,gam1,gam2,gampl,gammi)  
         ff=fact*(gam1*cosh(e)+gam2*fact2*d) 
         sum=ff 
         e=exp(e)
         p=0.5d0*e/gampl 
         q=0.5d0/(e*gammi)  
         c=1.d0 
         d=x2*x2 
         sum1=p 
         do i=1,MAXIT
            ff=(i*ff+p+q)/(i*i-xmu2) 
            c=c*d/i 
            p=p/(i-xmu) 
            q=q/(i+xmu) 
            del=c*ff
            sum=sum+del 
            del1=c*(p-i*ff)
            sum1=sum1+del1 
            if(abs(del).lt.abs(sum)*EPS)goto 2 
         enddo  
         pause 'bessk series failed to converge' 
 2       continue
         
         rkmu=sum 
         rk1=sum1*xi2 
      else  
         b=2.d0*(1.d0+x)
         d=1.d0/b 
         delh=d 
         h=delh 
         q1=0.d0  
         q2=1.d0 
         a1=.25d0-xmu2 
         c=a1 
         q=c 
         a=-a1 
         s=1.d0+q*delh 
         do i=2,MAXIT
            a=a-2*(i-1) 
            c=-a*c/i 
            qnew=(q1-b*q2)/a 
            q1=q2 
            q2=qnew 
            q=q+c*qnew
            b=b+2.d0 
            d=1.d0/(b+a*d) 
            delh=(b*d-1.d0)*delh 
            h=h+delh 
            dels=q*delh
            s=s+dels 
            if(abs(dels/s).lt.EPS)goto 3 
         enddo  
         pause 'bessik: failure to converge in cf2' 
 3       continue
         h=a1*h 
!
!     To compute scaled quantities
!
!         rkmu=sqrt(PI/(2.d0*x))*exp(-x)/s 
!
         rkmu=sqrt(PI/(2.d0*x))/s 
         rk1=rkmu*(xmu+x+.5d0-h)*xi
      endif 
      rkmup=xmu*xi*rkmu-rk1 
      rimu=xi/(f*rkmu-rkmup)  
      ri=(rimu*ril1)/ril  
      rip=(rimu*rip1)/ril 
      do i=1,nl 
         rktemp=(xmu+i)*xi2*rk1+rkmu 
         rkmu=rk1 
         rk1=rktemp 
      enddo  
      rk=rkmu
      rkp=xnu*xi*rkmu-rk1 
      return 
    END
    SUBROUTINE GEINV(N,U)
    !SUBROUTINE FOR INVERT A SYMMETRICAL MATRIX USE LAPACK ROUTINE SYTRF AND SYTRI
    !INPUT IS U(N,N)
    !OUTPUT IS U CANTAIN THE INVERT AND COVER THE ORGIN MATRIX
    IMPLICIT NONE
    INTEGER N,INFO
    INTEGER,DIMENSION(:) :: IPIV(N)
    REAL*8,DIMENSION(:,:)::U(N,N) 
    REAL*8,DIMENSION(:) :: V(N)
    INTEGER NN
    NN=N
    !CALL DSYTRF('U',NN,U,NN,IPIV,V,2*NN,INFO)
    
    !PRINT*,'INFO=',INFO
    !IF(INFO.NE.0) STOP 'THE INFO IN DSYTRF IS WRONG'
    !CALL DSYTRI('U',NN,U,NN,IPIV,V,INFO)
    !IF(INFO.NE.0) STOP 'THE INFO IN DSYTRI IS WRONG '
    !PRINT*,'U=',U
    CALL DGETRF(NN,NN,U,NN,IPIV,INFO)
    CALL DGETRI(N,U,NN,IPIV,V,N,INFO)
    !STOP
    END

    SUBROUTINE SYMINV(N,U)
    !SUBROUTINE FOR INVERT A SYMMETRICAL MATRIX USE LAPACK ROUTINE SYTRF AND SYTRI
    !INPUT IS U(N,N)
    !OUTPUT IS U CANTAIN THE INVERT AND COVER THE ORGIN MATRIX
    IMPLICIT NONE
    INTEGER N,INFO
    INTEGER,DIMENSION(:) :: IPIV(N)
    REAL*8,DIMENSION(:,:)::U(N,N)
    REAL*8,DIMENSION(:) :: V(2*N)
    INTEGER NN
    NN=N
    CALL DSYTRF('U',NN,U,NN,IPIV,V,2*NN,INFO)

    !PRINT*,'INFO=',INFO
    !IF(INFO.NE.0) STOP 'THE INFO IN DSYTRF IS WRONG'
    CALL DSYTRI('U',NN,U,NN,IPIV,V,INFO)
    !IF(INFO.NE.0) STOP 'THE INFO IN DSYTRI IS WRONG '
    !PRINT*,'U=',U
!    CALL DGETRF(NN,NN,U,NN,IPIV,INFO)
!    CALL DGETRI(N,U,NN,IPIV,V,N,INFO)
    !STOP
    END

    
      SUBROUTINE KSYM(AK,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION AK(N,N)
      DO 10 I=1,N
      DO 10 J=1,I
      TMP=0.5D0*(AK(I,J)+AK(J,I))
      AK(I,J)=TMP
      AK(J,I)=TMP
   10 CONTINUE
      RETURN
    END
	SUBROUTINE KTOS(XK, PR, PI, N)

	IMPLICIT DOUBLE PRECISION (A-H, O-Z)

	DIMENSION XK(N,N), XK2(N,N)
	DIMENSION UMASK2(N,N), UMENOSK2(N,N)
	DIMENSION PR(N,N), PI(N,N)
	DIMENSION IS(N), JS(N)

	XK2 = MATMUL(XK,XK)

	DO I = 1, N
	DO J = 1, N
		IF(I==J)THEN
			UMASK2(I,J) = 1.d0+XK2(I,J)
		ELSE
			UMASK2(I,J) = XK2(I,J)
		END IF
	END DO
	END DO

	CALL GEINV(N,UMASK2)

	DO I = 1, N
	DO J = 1, N
		IF(I==J)THEN
			UMENOSK2(I,J) = 1.d0-XK2(I,J)
		ELSE
			UMENOSK2(I,J) = -XK2(I,J)
		END IF
	END DO
	END DO

	PR = MATMUL(UMENOSK2,UMASK2)
	PI = MATMUL(XK,UMASK2)
	PI = 2.d0*PI

	RETURN

    END SUBROUTINE
    
	SUBROUTINE BRINV(A,N)
	DIMENSION A(N,N),IS(N),JS(N)
	DOUBLE PRECISION A,T,D
	L=1
	DO 100 K=1,N
	  D=0.0
	  DO 10 I=K,N
	  DO 10 J=K,N
	    IF (ABS(A(I,J)).GT.D) THEN
	      D=ABS(A(I,J))
	      IS(K)=I
	      JS(K)=J
	    END IF
10	  CONTINUE
	  IF (D+1.0.EQ.1.0) THEN
	    L=0
	    WRITE(*,20)
	    RETURN
	  END IF
20	  FORMAT(1X,'ERR**NOT INV')
	  DO 30 J=1,N
	    T=A(K,J)
	    A(K,J)=A(IS(K),J)
	    A(IS(K),J)=T
30	  CONTINUE
	  DO 40 I=1,N
	    T=A(I,K)
	    A(I,K)=A(I,JS(K))
	    A(I,JS(K))=T
40	  CONTINUE
	  A(K,K)=1/A(K,K)
	  DO 50 J=1,N
	    IF (J.NE.K) THEN
	      A(K,J)=A(K,J)*A(K,K)
	    END IF
50	  CONTINUE
	  DO 70 I=1,N
	    IF (I.NE.K) THEN
	      DO 60 J=1,N
	        IF (J.NE.K) THEN
	          A(I,J)=A(I,J)-A(I,K)*A(K,J)
	        END IF
60	      CONTINUE
	    END IF
70	  CONTINUE
	  DO 80 I=1,N
	    IF (I.NE.K) THEN
	      A(I,K)=-A(I,K)*A(K,K)
	    END IF
80	  CONTINUE
100	CONTINUE
	DO 130 K=N,1,-1
	  DO 110 J=1,N
	    T=A(K,J)
	    A(K,J)=A(JS(K),J)
	    A(JS(K),J)=T
110	  CONTINUE
	  DO 120 I=1,N
	    T=A(I,K)
	    A(I,K)=A(I,IS(K))
	    A(I,IS(K))=T
120	  CONTINUE
130	CONTINUE
	RETURN
	END
