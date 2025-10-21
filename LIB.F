
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE bessjy(x,xnu,rj,ry,rjp,ryp) 
      INTEGER MAXIT 
      DOUBLE PRECISION rj,rjp,ry,ryp,x,xnu,XMIN 
      DOUBLE PRECISION EPS,FPMIN,PI 
      PARAMETER(EPS=1.e-16,FPMIN=1.e-30,MAXIT=1000000,XMIN=2., 
     *     PI=3.141592653589793d0) 

C     USES beschb
C     Return the Bessel functions rj=J_v, ry=Y_v and their derivatives
C     rjp=J'_v, ryp=Y'_v for positive x and for xvu=v.ge.0.
C     The relative accuracy is whithin one or two digits of EPS, except
C     near a zero of one of the functions, where EPS controls its absolute 
C     accuracy. FPMIN is a number close to the machine'smallest floating number.
C     All internal arithmetic is in souble precision. To convert 
C     the entire routine
C     to double precision change REAL declaration above and 
C     decrease EPS to 10**-16.
C     Also convert beschb.

      INTEGER i,isign,l,nl 
      DOUBLE PRECISIONa,b,br,bi,c,cr,ci,d,del,del1,den,di,dlr,dli, 
     *     dr,e,f,fact,fact2,fact3,ff,gam,gam1,gam2,gammi,gampl,h,
     *     p,pimu,pimu2,q,r,rjl,rjl1,rjmu,rjp1,rjpl,rjtemp,ry1, 
     *     rymu,rymup,rytemp,sum,sum1,temp,w,x2,xi,xi2,xmu,xmu2
      
      if(x.le.0..or.xnu.lt.0.) pause 'bad arguments in bessjy'
      
      if(x.lt.XMIN)then 
         nl=int(xnu+.5d0) 
      else
         nl=max(0,int(xnu-x+1.5d0)) 
      endif 
         xmu=xnu-nl 
         xmu2=xmu*xmu 
         xi=1.d0/x
         xi2=2.d0*xi
         w=xi2/PI
         isign=1 
         h=xnu*xi
         if(h.lt.FPMIN)h=FPMIN 
         b=xi2*xnu 
         d=0.d0 
         c=h 
         do i=1,MAXIT
            b=b+xi2 
            d=b-d 
            if(abs(d).lt.FPMIN)d=FPMIN 
            c=b-1.d0/c
            if(abs(c).lt.FPMIN)c=FPMIN 
            d=1.d0/d 
            del=c*d 
            h=del*h
            if(d.lt.0.d0)isign=-isign 
            if(abs(del-1.d0).lt.EPS)goto 1 
         enddo 

c         pause 'x too large in bessjy; try asymptotic expansion' 
         write(6,*)'x too large in bessjy; trying asymptotic expansion' 
         wcos=dcos(x-pi*.5d0*xnu-pi*.25d0)
         wsin=dsin(x-pi*.5d0*xnu-pi*.25d0)
         wk=sqrt(2.d0/(x*pi))
         wkk=x*sqrt(2.d0*pi*x)
         rj=wk*wcos
         ry=wk*wsin
         rjp=(-1.d0/wkk)*wcos-ry
         ryp=(-1.d0/wkk)*wsin+rj
         goto 4
1        continue
         
         rjl=isign*FPMIN 
         rjpl=h*rjl 
         rjl1=rjl 
         rjp1=rjpl 
         fact=xnu*xi 
         do l=nl,1,-1
            rjtemp=fact*rjl+rjpl 
            fact=fact-xi 
            rjpl=fact*rjtemp-rjl 
            rjl=rjtemp
         enddo  
         if(rjl.eq.0.d0)rjl=EPS 
         f=rjpl/rjl
         if(x.lt.XMIN) then
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
            enDif 
            call beschb(xmu,gam1,gam2,gampl,gammi)  
            ff=2.d0/PI*fact*(gam1*cosh(e)+gam2*fact2*d) 
            e=exp(e)
            p=e/(gampl*PI) 
            q=1.d0/(e*PI*gammi) 
            pimu2=0.5d0*pimu
            if(abs(pimu2).lt.EPS)then
               fact3=1.d0 
            else
               fact3=sin(pimu2)/pimu2 
            endif 
            r=PI*pimu2*fact3*fact3 
            c=1.d0 
            d=-x2*x2
            sum=ff+r*q 
            sum1=p
            do  i=1,MAXIT
               ff=(i*ff+p+q)/(i*i-xmu2) 
               c=c*d/i 
               p=p/(i-xmu) 
               q=q/(i+xmu)
               del=c*(ff+r*q) 
               sum=sum+del 
               del1=c*p-i*del 
               sum1=sum1+del1
               if(abs(del).lt.(1.d0+abs(sum))*EPS)goto 2 
            enddo  
            pause 'bessy series failed to converge'
2           continue
            rymu=-sum 
            ry1=-sum1*xi2 
            rymup=xmu*xi*rymu-ry1 
            rjmu=w/(rymup-f*rymu)
         else 
            a=.25d0-xmu2 
            p=-.5d0*xi 
            q=1.d0 
            br=2.d0*x 
            bi=2.d0
            fact=a*xi/(p*p+q*q) 
            cr=br+q*fact 
            ci=bi+p*fact 
            den=br*br+bi*bi
            dr=br/den 
            di=-bi/den 
            dlr=cr*dr-ci*di 
            dli=cr*di+ci*dr 
            temp=p*dlr-q*dli
            q=p*dli+q*dlr 
            p=temp 
            do i=2,MAXIT
               a=a+2*(i-1) 
               bi=bi+2.d0 
               dr=a*dr+br 
               di=a*di+bi
               if(abs(dr)+abs(di).lt.FPMIN)dr=FPMIN 
               fact=a/(cr*cr+ci*ci)
               cr=br+cr*fact 
               ci=bi-ci*fact 
               if(abs(cr)+abs(ci).lt.FPMIN)cr=FPMIN
               den=dr*dr+di*di 
               dr=dr/den 
               di=-di/den 
               dlr=cr*dr-ci*di 
               dli=cr*di+ci*dr
               temp=p*dlr-q*dli 
               q=p*dli+q*dlr 
               p=temp
               if(abs(dlr-1.d0)+abs(dli).lt.EPS)goto 3 
            enddo  
            pause 'cf2 failed in bessjy' 
3            continue
            gam=(p-f)/q 
            rjmu=sqrt(w/((p-f)*gam+q))
            rjmu=sign(rjmu,rjl) 
            rymu=rjmu*gam 
            rymup=rymu*(p+q/gam)
            ry1=xmu*xi*rymu-rymup 
         endif 
         fact=rjmu/rjl
         rj=rjl1*fact  
         rjp=rjp1*fact 
         do  i=1,nl
            rytemp=(xmu+i)*xi2*ry1-rymu 
            rymu=ry1 
            ry1=rytemp 
         enddo  
         ry=rymu
         ryp=xnu*xi*rymu-ry1 
4        return 
         END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE beschb(x,gam1,gam2,gampl,gammi) 
      INTEGER NUSE1,NUSE2 
      DOUBLE PRECISION gam1,gam2,gammi,gampl,x 
      PARAMETER (NUSE1=7,NUSE2=8) 

C     USES chebev
C     Evluata Gamma1 and Gamma2 by Chebischev expansion for abs(x).le.0.5.
C     Also returns 1/Gamma(1+x) and 1/Gamma(1-x).
C     If converting to double precision, set NUSE1=7, NUSE2=8.

      DOUBLE PRECISION xx,c1(7),c2(8),chebev 
      SAVE c1,c2 
      DATA c1/-1.142022680371168d0,6.5165112670737d-3, 
     *     3.087090173086d-4,-3.4706269649d-6,6.9437664d-9, 
     *     3.67795d-11,-1.356d-13/
      
      DATA c2/1.843740587300905d0,-7.68528408447867d-2, 
     *     1.2719271366546d-3,-4.9717367042d-6,-3.31261198d-8, 
     *     2.423096d-10,-1.702d-13,-1.49d-15/
      
      xx=8.d0*x*x-1.d0 
      gam1=chebev(-1.0d0,1.0d0,c1,NUSE1,xx)
      gam2=chebev(-1.0d0,1.0d0,c2,NUSE2,xx) 
      gampl=gam2-x*gam1 
      gammi=gam2+x*gam1
      return 
      END
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function chebev(a,b,c,m,x)
      implicit double precision (a-h,o-z)
      dimension c(m)
      if ((x-a)*(x-b).gt.0.d0) pause 'x not in range.'
      d=0.d0
      dd=0.d0
      y=(2.d0*x-a-b)/(b-a)
      y2=2.d0*y
      do 11 j=m,2,-1
         sv=d
         d=y2*d-dd+c(j)
         dd=sv
 11   continue
      chebev=y*d-dd+0.5d0*c(1)
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	DOUBLE PRECISION FUNCTION DGAMMA(X)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	COMMON/A/ A(11)
	DATA A/1.0000000000D0,0.4227843370D0,
     *         0.4118402518D0,0.0815782188D0,
     *        0.0742379076D0,-0.0002109075D0,
     *        0.0109736958D0,-0.0024667480D0,
     *        0.0015397681D0,-0.0003442342D0,
     *        0.0000677106D0/	
	IF(X.GT.1.D0) GOTO 10
	T=X
	GOTO 50
10	IF(X.GT.2.D0) GOTO 20
	T=X-1.D0
	GOTO 50
20	IF(X.GT.3.D0) RETURN
	T=X-2.D0
50	CONTINUE
	P=A(11)
	DO 1 K=10,1,-1
1	P=T*P+A(K)
	IF(X.GT.1.D0) GOTO 110
	DGAMMA=P/X/(X+1.D0)
	GOTO 550
110	IF(X.GT.2.D0) GOTO 220
	DGAMMA=P/X
	GOTO 550
220	IF(X.GT.3.D0) RETURN
	DGAMMA=P
550	CONTINUE
	RETURN
	END

	  !------Spline expansion:Cubic spline interpolation----------
      SUBROUTINE SPLINE(N,X,Y,A,B,C,D)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(N),Y(N),A(N),B(N),C(N),D(N),U(N)
      NM1=N-1
      NM2=N-2
      DO 1 I=1,N
 1    A(I)=Y(I)
      DO 2 I=2,NM1
      IM1=I-1  
      IP1=I+1
      HI=X(IP1)-X(I)
      HIM1=X(I)-X(IM1)
      B(IM1)=HIM1 
      U(IM1)=HI
      D(IM1)=2*(HI+HIM1)
 2    Y(IM1)=3*((A(IP1)-A(I))/HI-(A(I)-A(IM1))/HIM1)
      CALL TRIMAT(NM2,N,Y,B,U,D,C)
      DO 3 I=1,NM2
      J=N-I
      K=J-1  
 3    C(J)=C(K)
      C(1)=0.
      C(N)=0.
      DO 4 I=1,NM1
      IP1=I+1
      HI=X(IP1)-X(I)
      B(I)=(A(IP1)-A(I))/HI-.333333332*HI*(C(IP1)+C(I)+C(I))
 4    D(I)=.333333332*(C(IP1)-C(I))/HI
      RETURN
      END SUBROUTINE


       SUBROUTINE TRIMAT(N,NQ,A,B,C,D,X)
      IMPLICIT REAL*8 (A-H,O-Z)   
      DIMENSION D(nq),B(nq),C(nq),X(nq),A(nq),GAMMA(nq),G(nq)
      NM1=N-1
      ALPHA=D(1)
      GAMMA(1)=C(1)/ALPHA
      G(1)=A(1)/ALPHA
      DO 1 I=2,N
      ALPHA   =D(I)-B(I)*GAMMA(I-1)
      G(I)=(A(I)-B(I)*G(I-1))/ALPHA
      IF(I.EQ.N)GOTO 1   
      GAMMA(I)=C(I)/ALPHA
 1    CONTINUE  
      X(N)=G(N)
      DO 2 I=1,NM1
      J=N-I
2     X(J)=G(J)-GAMMA(J)*X(J+1)
      RETURN  
      END SUBROUTINE 
c-------------------------------------------------------------------------
      SUBROUTINE GAUSSQ(KIND, N, ALPHA, BETA, KPTS, ENDPTS, B, T, W)    
      IMPLICIT REAL*8 (A-H,O-Z) 
C                                                                       
C           THIS SET OF ROUTINES COMPUTES THE NODES X(I) AND WEIGHTS    
C        C(I) FOR GAUSSIAN-TYPE QUADRATURE RULES WITH PRE-ASSIGNED      
C        NODES.  THESE ARE USED WHEN ONE WISHES TO APPROXIMATE          
C                                                                       
C                 INTEGRAL (FROM A TO B)  F(X) W(X) DX                  
C                                                                       
C                              N                                         
C        BY                   SUM C  F(X )                               
C                             I=1  I    I                                
C                                                                        
C        HERE W(X) IS ONE OF SIX POSSIBLE NON-NEGATIVE WEIGHT            
C        FUNCTIONS (LISTED BELOW), AND F(X) IS THE                       
C        FUNCTION TO BE INTEGRATED.  GAUSSIAN QUADRATURE IS PARTICULARLY 
C        USEFUL ON INFINITE INTERVALS (WITH APPROPRIATE WEIGHT           
C        FUNCTIONS), SINCE THEN OTHER TECHNIQUES OFTEN FAIL.             
C                                                                        
C           ASSOCIATED WITH EACH WEIGHT FUNCTION W(X) IS A SET OF        
C        ORTHOGONAL POLYNOMIALS.  THE NODES X(I) ARE JUST THE ZEROES     
C        OF THE PROPER N-TH DEGREE POLYNOMIAL.                           
C                                                                        
C     INPUT PARAMETERS                                                   
C                                                                        
C        KIND     AN INTEGER BETWEEN 1 AND 6 GIVING THE TYPE OF          
C                 QUADRATURE RULE                                        
C                                                                        
C        KIND = 1=  LEGENDRE QUADRATURE, W(X) = 1 ON (-1, 1)             
C        KIND = 2=  CHEBYSHEV QUADRATURE OF THE FIRST KIND               
C                   W(X) = 1/SQRT(1 - X*X) ON (-1, +1)                   
C        KIND = 3=  CHEBYSHEV QUADRATURE OF THE SECOND KIND              
C                   W(X) = SQRT(1 - X*X) ON (-1, 1)                      
C        KIND = 4=  HERMITE QUADRATURE, W(X) = EXP(-X*X) ON              
C                   (-INFINITY, +INFINITY)                               
C        KIND = 5=  JACOBI QUADRATURE, W(X) = (1-X)**ALPHA * (1+X)**     
C                   BETA ON (-1, 1), ALPHA, BETA .GT. -1.                
C                   NOTE= KIND=2 AND 3 ARE A SPECIAL CASE OF THIS.       
C        KIND = 6=  GENERALIZED LAGUERRE QUADRATURE, W(X) = EXP(-X)*     
C                   X**ALPHA ON (0, +INFINITY), ALPHA .GT. -1            
C                                                                        
C        N        THE NUMBER OF POINTS USED FOR THE QUADRATURE RULE      
C        ALPHA    REAL PARAMETER USED ONLY FOR GAUSS-JACOBI AND GAUSS-   
C                 LAGUERRE QUADRATURE (OTHERWISE USE 0.).                
C        BETA     REAL PARAMETER USED ONLY FOR GAUSS-JACOBI QUADRATURE-- 
C                 (OTHERWISE USE 0.).                                    
C        KPTS     (INTEGER) NORMALLY 0, UNLESS THE LEFT OR RIGHT END-    
C                 POINT (OR BOTH) OF THE INTERVAL IS REQUIRED TO BE A    
C                 NODE (THIS IS CALLED GAUSS-RADAU OR GAUSS-LOBATTO      
C                 QUADRATURE).  THEN KPTS IS THE NUMBER OF FIXED         
C                 ENDPOINTS (1 OR 2).                                    
C        ENDPTS   REAL ARRAY OF LENGTH 2.  CONTAINS THE VALUES OF        
C                 ANY FIXED ENDPOINTS, IF KPTS = 1 OR 2.                 
C        B        REAL SCRATCH ARRAY OF LENGTH N                         
C                                                                        
C     OUTPUT PARAMETERS (BOTH ARRAYS OF LENGTH N)                        
C                                                                        
C        T        WILL CONTAIN THE DESIRED NODES X(1),,,X(N)             
C        W        WILL CONTAIN THE DESIRED WEIGHTS C(1),,,C(N)           
C                                                                        
C     SUBROUTINES REQUIRED                                               
C                                                                        
C        GBSLVE, CLASS, AND GBTQL2 ARE PROVIDED. UNDERFLOW MAY SOMETIMES 
C        OCCUR, BUT IT IS HARMLESS IF THE UNDERFLOW INTERRUPTS ARE       
C        TURNED OFF AS THEY ARE ON THIS MACHINE.                         
C                                                                        
C     ACCURACY                                                           
C                                                                        
C        THE ROUTINE WAS TESTED UP TO N = 512 FOR LEGENDRE QUADRATURE,   
C        UP TO N = 136 FOR HERMITE, UP TO N = 68 FOR LAGUERRE, AND UP    
C        TO N = 10 OR 20 IN OTHER CASES.  IN ALL BUT TWO INSTANCES,      
C        COMPARISON WITH TABLES IN REF. 3 SHOWED 12 OR MORE SIGNIFICANT  
C        DIGITS OF ACCURACY.  THE TWO EXCEPTIONS WERE THE WEIGHTS FOR    
C        HERMITE AND LAGUERRE QUADRATURE, WHERE UNDERFLOW CAUSED SOME    
C        VERY SMALL WEIGHTS TO BE SET TO ZERO.  THIS IS, OF COURSE,      
C        COMPLETELY HARMLESS.                                            
C                                                                        
C     METHOD                                                             
C                                                                        
C           THE COEFFICIENTS OF THE THREE-TERM RECURRENCE RELATION       
C        FOR THE CORRESPONDING SET OF ORTHOGONAL POLYNOMIALS ARE         
C        USED TO FORM A SYMMETRIC TRIDIAGONAL MATRIX, WHOSE              
C        EIGENVALUES (DETERMINED BY THE IMPLICIT QL-METHOD WITH          
C        SHIFTS) ARE JUST THE DESIRED NODES.  THE FIRST COMPONENTS OF    
C        THE ORTHONORMALIZED EIGENVECTORS, WHEN PROPERLY SCALED,         
C        YIELD THE WEIGHTS.  THIS TECHNIQUE IS MUCH FASTER THAN USING A  
C        ROOT-FINDER TO LOCATE THE ZEROES OF THE ORTHOGONAL POLYNOMIAL.  
C        FOR FURTHER DETAILS, SEE REF. 1.  REF. 2 CONTAINS DETAILS OF    
C        GAUSS-RADAU AND GAUSS-LOBATTO QUADRATURE ONLY.                  
C                                                                        
C     REFERENCES                                                         
C                                                                        
C        1.  GOLUB, G. H., AND WELSCH, J. H.,  CALCULATION OF GAUSSIAN   
C            QUADRATURE RULES,  MATHEMATICS OF COMPUTATION 23 (APRIL,    
C            1969), PP. 221-230.                                         
C        2.  GOLUB, G. H.,  SOME MODIFIED MATRIX EIGENVALUE PROBLEMS,    
C            SIAM REVIEW 15 (APRIL, 1973), PP. 318-334 (SECTION 7).      
C        3.  STROUD AND SECREST, GAUSSIAN QUADRATURE FORMULAS, PRENTICE- 
C            HALL, ENGLEWOOD CLIFFS, N.J., 1966.                         
C                                                                       
C     ..................................................................
C                                                                       
      REAL*8  MUZERO                                                    
      DIMENSION  B(N),T(N),W(N),ENDPTS(2)                               
C                                                                       
      CALL CLASS (KIND, N, ALPHA, BETA, B, T, MUZERO)                   
C                                                                       
C           THE MATRIX OF COEFFICIENTS IS ASSUMED TO BE SYMMETRIC.      
C           THE ARRAY T CONTAINS THE DIAGONAL ELEMENTS, THE ARRAY       
C           B THE OFF-DIAGONAL ELEMENTS.                                
C           MAKE APPROPRIATE CHANGES IN THE LOWER RIGHT 2 BY 2          
C           SUBMATRIX.                                                  
C                                                                       
      IF (KPTS.EQ.0)  GO TO 100                                         
      IF (KPTS.EQ.2)  GO TO  50                                         
C                                                                       
C           IF KPTS=1, ONLY T(N) MUST BE CHANGED                        
C                                                                       
      T(N) =GBSLVE(ENDPTS(1), N, T, B)*B(N-1)**2 + ENDPTS(1)            
      GO TO 100                                                         
C                                                                       
C           IF KPTS=2, T(N) AND B(N-1) MUST BE RECOMPUTED               
C                                                                       
   50 GAM =GBSLVE(ENDPTS(1), N, T, B)                                   
      T1 = ((ENDPTS(1) - ENDPTS(2))/(GBSLVE(ENDPTS(2), N, T, B) - GAM)) 
      B(N-1) =  SQRT(T1)                                                
      T(N) = ENDPTS(1) + GAM*T1                                         
C                                                                       
C           NOTE THAT THE INDICES OF THE ELEMENTS OF B RUN FROM 1 TO N-1
C           AND THUS THE VALUE OF B(N) IS ARBITRARY.                    
C           NOW COMPUTE THE EIGENVALUES OF THE SYMMETRIC TRIDIAGONAL    
C           MATRIX, WHICH HAS BEEN MODIFIED AS NECESSARY.               
C           THE METHOD USED IS A QL-TYPE METHOD WITH ORIGIN SHIFTING    
C                                                                       
  100 W(1) = 1.0D0                                                      
      DO 105 I = 2, N                                                   
  105    W(I) = 0.0D0                                                   
C                                                                       
      CALL GBTQL2 (N, T, B, W, IERR)                                    
      DO 110 I = 1, N                                                   
  110    W(I) = MUZERO * W(I) * W(I)                                    
C                                                                       
      RETURN                                                            
      END                                                               
C                                                                       
C                                                                       
C                                                                       
      DOUBLE PRECISION FUNCTION GBSLVE(SHIFT, N, A, B)                                   
      IMPLICIT REAL*8 (A-H,O-Z) 
C                                                                       
C       THIS PROCEDURE PERFORMS ELIMINATION TO SOLVE FOR THE            
C       N-TH COMPONENT OF THE SOLUTION DELTA TO THE EQUATION            
C                                                                       
C             (JN - SHIFT*IDENTITY) * DELTA  = EN,                      
C                                                                       
C       WHERE EN IS THE VECTOR OF ALL ZEROES EXCEPT FOR 1 IN            
C       THE N-TH POSITION.                                              
C                                                                       
C       THE MATRIX JN IS SYMMETRIC TRIDIAGONAL, WITH DIAGONAL           
C       ELEMENTS A(I), OFF-DIAGONAL ELEMENTS B(I).  THIS EQUATION       
C       MUST BE SOLVED TO OBTAIN THE APPROPRIATE CHANGES IN THE LOWER   
C       2 BY 2 SUBMATRIX OF COEFFICIENTS FOR ORTHOGONAL POLYNOMIALS.    
C                                                                       
C                                                                       
      DIMENSION  A(N),B(N)                                              
C                                                                       
      ALPHA = A(1) - SHIFT                                              
      NM1 = N - 1                                                       
      DO 10 I = 2, NM1                                                  
   10    ALPHA = A(I) - SHIFT - B(I-1)**2/ALPHA                         
      GBSLVE = 1.0D0  /ALPHA                                            
      RETURN                                                            
      END                                                               
C                                                                       
C                                                                       
C                                                                       
      SUBROUTINE CLASS(KIND, N, ALPHA, BETA, B, A, MUZERO)              
      IMPLICIT REAL*8 (A-H,O-Z) 
C                                                                       
C           THIS PROCEDURE SUPPLIES THE COEFFICIENTS A(J), B(J) OF THE  
C        RECURRENCE RELATION                                            
C                                                                       
C             B P (X) = (X - A ) P   (X) - B   P   (X)                  
C              J J            J   J-1       J-1 J-2                     
C                                                                       
C        FOR THE VARIOUS CLASSICAL (NORMALIZED) ORTHOGONAL POLYNOMIALS, 
C        AND THE ZERO-TH MOMENT                                         
C                                                                       
C             MUZERO = INTEGRAL W(X) DX                                 
C                                                                       
C        OF THE GIVEN POLYNOMIAL   WEIGHT FUNCTION W(X).  SINCE THE     
C        POLYNOMIALS ARE ORTHONORMALIZED, THE TRIDIAGONAL MATRIX IS     
C        GUARANTEED TO BE SYMMETRIC.                                    
C                                                                       
C           THE INPUT PARAMETER ALPHA IS USED ONLY FOR LAGUERRE AND     
C        JACOBI POLYNOMIALS, AND THE PARAMETER BETA IS USED ONLY FOR    
C        JACOBI POLYNOMIALS.  THE LAGUERRE AND JACOBI POLYNOMIALS       
C        REQUIRE THE GAMMA FUNCTION.                                    
C                                                                       
C     ..................................................................
C                                                                       
      DIMENSION  A(N),B(N)                                              
      REAL*8  MUZERO                                                    
      DATA PI / 3.141592653589793D0  /                                  
C                                                                       
      NM1 = N - 1                                                       
      GO TO (10, 20, 30, 40, 50, 60), KIND                              
C                                                                       
C              KIND = 1=  LEGENDRE POLYNOMIALS P(X)                     
C              ON (-1, +1), W(X) = 1.                                   
C                                                                       
   10 MUZERO = 2.0D0                                                    
      DO 11 I = 1, NM1                                                  
         A(I) = 0.0D0                                                   
         ABI = I                                                        
   11    B(I) = ABI/ SQRT(4*ABI*ABI - 1.0D0  )                          
      A(N) = 0.0D0                                                      
      RETURN                                                            
C                                                                       
C              KIND = 2=  CHEBYSHEV POLYNOMIALS OF THE FIRST KIND T(X)  
C              ON (-1, +1), W(X) = 1 / SQRT(1 - X*X)                    
C                                                                       
   20 MUZERO = PI                                                       
      DO 21 I = 1, NM1                                                  
         A(I) = 0.0D0                                                   
   21    B(I) = 0.5D0                                                   
      B(1) =  SQRT(0.5D0  )                                             
      A(N) = 0.0D0                                                      
      RETURN                                                            
C                                                                       
C              KIND = 3=  CHEBYSHEV POLYNOMIALS OF THE SECOND KIND U(X) 
C              ON (-1, +1), W(X) = SQRT(1 - X*X)                        
C                                                                       
   30 MUZERO = PI/2.0D0                                                 
      DO 31 I = 1, NM1                                                  
         A(I) = 0.0D0                                                   
   31    B(I) = 0.5D0                                                   
      A(N) = 0.0D0                                                      
      RETURN                                                            
C                                                                       
C              KIND = 4=  HERMITE POLYNOMIALS H(X) ON (-INFINITY,       
C              +INFINITY), W(X) = EXP(-X**2)                            
C                                                                       
   40 MUZERO =  SQRT(PI)                                                
      DO 41 I = 1, NM1                                                  
         A(I) = 0.0D0                                                   
   41    B(I) =  SQRT(I/2.0D0  )                                        
      A(N) = 0.0D0                                                      
      RETURN                                                            
C                                                                       
C              KIND = 5=  JACOBI POLYNOMIALS P(ALPHA, BETA)(X) ON       
C              (-1, +1), W(X) = (1-X)**ALPHA + (1+X)**BETA, ALPHA AND   
C              BETA GREATER THAN -1                                     
C                                                                       
   50 AB = ALPHA + BETA                                                 
      ABI = 2.0D0   + AB                                                
      MUZERO = 2.0D0   ** (AB + 1.0D0  ) * DGAMMA(ALPHA + 1.0D0  ) * DGA
     VMMA(
     X BETA + 1.0D0  ) / DGAMMA(ABI)                                    
      A(1) = (BETA - ALPHA)/ABI                                         
      B(1) =  SQRT(4.0D0  *(1.0D0   + ALPHA)*(1.0D0   + BETA)/((ABI + 1.
     V0D0  )* 
     1  ABI*ABI))                                                       
      A2B2 = BETA*BETA - ALPHA*ALPHA                                    
      DO 51 I = 2, NM1                                                  
         ABI = 2.0D0  *I + AB                                           
         A(I) = A2B2/((ABI - 2.0D0  )*ABI)                              
   51    B(I) =  SQRT (4.0D0  *I*(I + ALPHA)*(I + BETA)*(I + AB)/       
     1   ((ABI*ABI - 1)*ABI*ABI))                                       
      ABI = 2.0D0  *N + AB                                              
      A(N) = A2B2/((ABI - 2.0D0  )*ABI)                                 
      RETURN                                                            
C                                                                       
C              KIND = 6=  LAGUERRE POLYNOMIALS L(ALPHA)(X) ON           
C              (0, +INFINITY), W(X) = EXP(-X) * X**ALPHA, ALPHA GREATER 
C              THAN -1.                                                 
C                                                                       
   60 MUZERO = DGAMMA(ALPHA + 1.0D0  )                                  
      DO 61 I = 1, NM1                                                  
         A(I) = 2.0D0  *I - 1.0D0   + ALPHA                             
   61    B(I) =  SQRT(I*(I + ALPHA))                                    
      A(N) = 2.0D0  *N - 1 + ALPHA                                      
      RETURN                                                            
      END                                                               
C     ------------------------------------------------------------------
C                                                                       
      SUBROUTINE GBTQL2(N, D, E, Z, IERR)                               
      IMPLICIT REAL*8 (A-H,O-Z) 
C                                                                       
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE IMTQL2,   
C     NUM. MATH. 12, 377-383(1968) BY MARTIN AND WILKINSON,             
C     AS MODIFIED IN NUM. MATH. 15, 450(1970) BY DUBRULLE.              
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 241-248(1971).   
C                                                                       
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND FIRST COMPONENTS OF THE 
C     EIGENVECTORS OF A SYMMETRIC TRIDIAGONAL MATRIX BY THE IMPLICIT QL 
C     METHOD, AND IS ADAPTED FROM THE EISPAK ROUTINE IMTQL2             
C                                                                       
C     ON INPUT=                                                         
C                                                                       
C        N IS THE ORDER OF THE MATRIX;                                  
C                                                                       
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX;          
C                                                                       
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX        
C          IN ITS FIRST N-1 POSITIONS.  E(N) IS ARBITRARY;              
C                                                                       
C        Z CONTAINS THE FIRST ROW OF THE IDENTITY MATRIX.               
C                                                                       
C      ON OUTPUT=                                                       
C                                                                       
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN          
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT BUT          
C          UNORDERED FOR INDICES 1, 2, ..., IERR-1;                     
C                                                                       
C        E HAS BEEN DESTROYED;                                          
C                                                                       
C        Z CONTAINS THE FIRST COMPONENTS OF THE ORTHONORMAL EIGENVECTORS
C          OF THE SYMMETRIC TRIDIAGONAL MATRIX.  IF AN ERROR EXIT IS    
C          MADE, Z CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE STORED 
C          EIGENVALUES;                                                 
C                                                                       
C        IERR IS SET TO                                                 
C                                                                       
C        IERR IS SET TO                                                 
C          ZERO       FOR NORMAL RETURN,                                
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN               
C                     DETERMINED AFTER 30 ITERATIONS.                   
C                                                                       
C     ------------------------------------------------------------------
C                                                                       
      INTEGER I, J, K, L, M, N, II, MML, IERR                           
      DIMENSION  D(N),E(N),Z(N)                                         
      REAL*8  MACHEP                                                    
C                                                                       
C     ========== MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING     
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.   
C                MACHEP = 16.0D0**(-13) FOR LONG FORM ARITHMETIC        
C                ON S360 ==========                                     
       MACHEP=1.0D-14                                                   
C                                                                       
      IERR = 0                                                          
      IF (N .EQ. 1) GO TO 1001                                          
C                                                                       
      E(N) = 0.0D0                                                      
      DO 240 L = 1, N                                                   
         J = 0                                                          
C     ========== LOOK FOR SMALL SUB-DIAGONAL ELEMENT ==========         
  105    DO 110 M = L, N                                                
            IF (M .EQ. N) GO TO 120                                     
            IF ( ABS(E(M)) .LE. MACHEP * ( ABS(D(M)) +  ABS(D(M+1))))   
     X         GO TO 120                                                
  110    CONTINUE                                                       
C                                                                       
  120    P = D(L)                                                       
         IF (M .EQ. L) GO TO 240                                        
         IF (J .EQ. 30) GO TO 1000                                      
         J = J + 1                                                      
C     ========== FORM SHIFT ==========                                  
         G = (D(L+1) - P) / (2.0D0   * E(L))                            
         R =  SQRT(G*G+1.0D0  )                                         
         G = D(M) - P + E(L) / (G +  SIGN(R, G))                        
         S = 1.0D0                                                      
         C = 1.0D0                                                      
         P = 0.0D0                                                      
         MML = M - L                                                    
C     ========== FOR I=M-1 STEP -1 UNTIL L DO -- ==========             
         DO 200 II = 1, MML                                             
            I = M - II                                                  
            F = S * E(I)                                                
            B = C * E(I)                                                
            IF ( ABS(F) .LT.  ABS(G)) GO TO 150                         
            C = G / F                                                   
            R =  SQRT(C*C+1.0D0  )                                      
            E(I+1) = F * R                                              
            S = 1.0D0   / R                                             
            C = C * S                                                   
            GO TO 160                                                   
  150       S = F / G                                                   
            R =  SQRT(S*S+1.0D0  )                                      
            E(I+1) = G * R                                              
            C = 1.0D0   / R                                             
            S = S * C                                                   
  160       G = D(I+1) - P                                              
            R = (D(I) - G) * S + 2.0D0   * C * B                        
            P = S * R                                                   
            D(I+1) = G + P                                              
            G = C * R - B                                               
C     ========== FORM FIRST COMPONENT OF VECTOR ==========              
            F = Z(I+1)                                                  
            Z(I+1) = S * Z(I) + C * F                                   
            Z(I) = C * Z(I) - S * F                                     
C                                                                       
  200    CONTINUE                                                       
C                                                                       
         D(L) = D(L) - P                                                
         E(L) = G                                                       
         E(M) = 0.0D0                                                   
         GO TO 105                                                      
  240 CONTINUE                                                          
C     ========== ORDER EIGENVALUES AND EIGENVECTORS ==========          
      DO 300 II = 2, N                                                  
         I = II - 1                                                     
         K = I                                                          
         P = D(I)                                                       
C                                                                       
         DO 260 J = II, N                                               
            IF (D(J) .GE. P) GO TO 260                                  
            K = J                                                       
            P = D(J)                                                    
  260    CONTINUE                                                       
C                                                                       
         IF (K .EQ. I) GO TO 300                                        
         D(K) = D(I)                                                    
         D(I) = P                                                       
C                                                                       
         P = Z(I)                                                       
         Z(I) = Z(K)                                                    
         Z(K) = P                                                       
C                                                                       
  300 CONTINUE                                                          
C                                                                       
      GO TO 1001                                                        
C     ========== SET ERROR -- NO CONVERGENCE TO AN                      
C                EIGENVALUE AFTER 30 ITERATIONS ==========              
 1000 IERR = L                                                          
 1001 RETURN                                                       
      END
      
      SUBROUTINE GET_INTERP(LENGTH,X,Y,NUM_A,A,B)
	IMPLICIT NONE
	INTEGER LENGTH,NUM_A,ISP,I,J
      REAL*8 DD,EPS0
	REAL*8,DIMENSION(NUM_A) :: A,B
	REAL*8,DIMENSION(LENGTH):: X,Y,ASP,BSP,CSP,DSP
	EPS0=1.D-10

	CALL SPLINE(LENGTH,X,Y,ASP,BSP,CSP,DSP)
	DO I=1,NUM_A
	IF(A(I)<X(1).OR.A(I)>(X(LENGTH)+EPS0)) THEN
	PRINT*,'A(',I,')=',A(I),'IS OUT OF RANGE !'
	STOP
	ENDIF
	DO J=1,LENGTH
	IF(J.NE.LENGTH) THEN
	IF(X(J)>A(I)) THEN
	ISP=J-1
	DD=A(I)-X(ISP)
	B(I)=ASP(ISP)+BSP(ISP)*DD+CSP(ISP)*DD*DD+DSP(ISP)*DD*DD*DD
	EXIT
	ENDIF
	ELSEIF(J.EQ.LENGTH) THEN
	IF(X(J)>A(I)) THEN
	ISP=J-1
	DD=A(I)-X(ISP)
	B(I)=ASP(ISP)+BSP(ISP)*DD+CSP(ISP)*DD*DD+DSP(ISP)*DD*DD*DD
	EXIT
	ELSE
	B(I)=Y(J)
	ENDIF
	ENDIF
	ENDDO
	ENDDO
	END 
