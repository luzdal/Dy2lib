  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
   !The module for getinterpolation
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
   MODULE INTEP
   CONTAINS
	  !----------------------------------------------------------
	  FUNCTION GETINTERPOLATION(LENGTH,X,Y,NUM_A,A)
	     IMPLICIT NONE
		 INTEGER         :: LENGTH,NUM_A,ISP,I,J
		 DOUBLE PRECISION::DD,EPS0=1.0D-10
		 DOUBLE PRECISION::GETINTERPOLATION(NUM_A),A(NUM_A),X(LENGTH),Y(LENGTH),X0(LENGTH),Y0(LENGTH)
		 DOUBLE PRECISION::ASP(LENGTH),BSP(LENGTH),CSP(LENGTH),DSP(LENGTH)
		 X0=X
		 Y0=Y
		 CALL SPLINE(LENGTH,X0,Y0,ASP,BSP,CSP,DSP)
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
				     GETINTERPOLATION(I)=ASP(ISP)+BSP(ISP)*DD+CSP(ISP)*DD*DD+DSP(ISP)*DD*DD*DD
				     EXIT
			      ENDIF
				ELSEIF(J.EQ.LENGTH) THEN
				  IF(X(J)>A(I)) THEN
				     ISP=J-1
					 DD=A(I)-X(ISP)
					 GETINTERPOLATION(I)=ASP(ISP)+BSP(ISP)*DD+CSP(ISP)*DD*DD+DSP(ISP)*DD*DD*DD
					 EXIT
					ELSE
					GETINTERPOLATION(I)=Y(J)
				   ENDIF
			 ENDIF
		  ENDDO
		ENDDO
	 END FUNCTION GETINTERPOLATION
	 END
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! parameter for mfg
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
MODULE MOD_SETPOT
!USE SP_PARAMETER,ONLY: RANK_S,RANK_T,RANK_B
IMPLICIT NONE
!REAL*8,PARAMETER :: R_LARGE_S = 10.D0
!REAL*8,PARAMETER :: R_LARGE_T = 10.D0
REAL*8,PARAMETER :: R_END=50.D0
REAL*8,PARAMETER :: R_MID=2.D0
REAL*8,PARAMETER :: NR_S=4000
REAL*8,PARAMETER :: NR_L=4000
REAL*8,PARAMETER :: EMAX=10.D0
REAL*8,PARAMETER :: FMASS=0.5D0 ! IN THE SCALED UNIT
REAL*8,PARAMETER :: DEGR=0.D0   !?
REAL*8,PARAMETER :: BETA=1.D0
INTEGER,PARAMETER :: NRMAX=NR_S+NR_L
INTEGER,PARAMETER :: NXMAX=1024
REAL*8,PARAMETER :: POT_TE=0.D0
REAL*8,PARAMETER :: POT_ATE=0.D0
REAL*8,PARAMETER :: ENERGY_DIFF=0.D0
REAL*8,ALLOCATABLE,DIMENSION(:) :: KE_MFG,DELTA_MFG,R_LARGE_MFG
INTEGER,ALLOCATABLE,DIMENSION (:) :: N_LARGE_MFG(:)
REAL*8 DR_S,DR_L,DX,XLENGTH,HAMI_V_MAX,HAMI_V_MIN,PMAX_X,DETX
REAL*8,DIMENSION(NRMAX) :: POT_R,JX_XTOR,WAV_MFG_B,WAV_MFG_ZERO,WAV_MFG_UP,WAV_MFG_B1,WAV_MFG_B2
!----------------------------------------------------------------------------------------
REAL*8,DIMENSION(NRMAX,3) :: POT_R_MUL   ! THE 3 COLUMN IS 1.D0/R/R
REAL*8,DIMENSION(0:NXMAX,3) :: POT_X_MUL
!----------------------------------------------------------------------------------------
REAL*8,DIMENSION(NRMAX,2) :: RTOX
REAL*8,DIMENSION(0:NXMAX) :: X_GRID,POT_X,JX,Trap_X
REAL*8,DIMENSION(1:NXMAX-1) :: R_MFG
REAL*8,DIMENSION (:,:) :: MATRIXT(NXMAX-1,NXMAX-1)
INTEGER,PARAMETER  :: MATZ=1
REAL*8,ALLOCATABLE:: EIGENVAL(:),EIGENVEC(:,:) 
INTEGER N_ATE
REAL*8 R_BEG
REAL*8 c6MFG,lamMFG,nu
real*8 potminr,potmin
END
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!SUBROUTINE AND FUNCTION FOR MFG
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!----------------------------------------------------------------
	    SUBROUTINE SETPOT()
        USE MOD_SETPOT
		IMPLICIT NONE
		REAL*8,PARAMETER :: PI=DACOS(-1.D0)
		INTEGER          :: GRID_MIN
		INTEGER          :: I,J,K
		DOUBLE PRECISION :: V_MIN
		DOUBLE PRECISION :: POT_ENV(NRMAX)
		REAL*8,EXTERNAL :: POTENTIAL
		REAL*8 R_TEMP
    
        DETX=BETA
		PMAX_X=PI/DETX   !IF IT IS THAT IN SCALED UNITS
        
        DR_S=(R_MID-R_BEG)/(NR_S-1)
		DR_L=(R_END-R_MID)/NR_L
		DO I=1,NRMAX
		IF(I.LE.NR_S) THEN
		R_TEMP=R_BEG+(I-1)*DR_S
		RTOX(I,1)=R_TEMP
		POT_R(I)=POTENTIAL(R_TEMP)
		ELSE
		R_TEMP=R_MID+(I-NR_S)*DR_L
		RTOX(I,1)=R_TEMP
		POT_R(I)=POTENTIAL(R_TEMP)
		END IF
		END DO
		CALL GET_POT_ENVELOPE(POT_ENV,V_MIN,GRID_MIN)
        
        
		CALL GET_RTOX(GRID_MIN,POT_ENV)
        PRINT*,'-------RTOX IS OBTAIN---------'
		DO I=1,NRMAX
			JX_XTOR(I)=PMAX_X/DSQRT(2.0D0*FMASS*(EMAX-POT_ENV(I)))

		ENDDO
		CALL GET_POT_X(V_MIN,GRID_MIN,POT_ENV)
		WRITE(*,*) 'XLENGTH IS ',XLENGTH,' .au.'
		!WRITE(*,*) 'THE NUMBER OF GRIDS IN X COORDINATE IS ',INT(XLENGTH)+1

		RETURN
	    END SUBROUTINE SETPOT


	    SUBROUTINE GET_POT_ENVELOPE(POT_ENV,VMIN,GRID_MIN)
        USE MOD_SETPOT
		IMPLICIT NONE
		DOUBLE PRECISION :: POT_ENV(NRMAX)
		INTEGER          :: GRID_MIN
		DOUBLE PRECISION :: R_ENV(NRMAX)
		DOUBLE PRECISION :: POT(NRMAX)
		DOUBLE PRECISION :: VMIN
		INTEGER          :: NVMIN
		INTEGER          :: I,J
		INTEGER          :: I_STATE
		
		R_ENV=RTOX(:,1)

		POT(1:NRMAX)=POT_R(1:NRMAX)
		POT_ENV(1:NRMAX)=POT(1:NRMAX)

		VMIN=MINVAL(POT(:))
		NVMIN=MINLOC(POT(:),DIM=1)
		DO I=1,NVMIN
			POT_ENV(I)=VMIN    !SO THAT THE JX AND THE RTOX
		ENDDO
		GRID_MIN=NVMIN
		RETURN
        END 

	    SUBROUTINE GET_RTOX(GRID_MIN,POT_ENV)
	    USE MOD_SETPOT

		IMPLICIT NONE
		INTEGER          :: GRID_MIN
		DOUBLE PRECISION :: POT_ENV(NRMAX)
		DOUBLE PRECISION :: SUM1,SUM2,SUM3,SUM4,INTERGRATION
		DOUBLE PRECISION :: S1,S2,SUM
		DOUBLE PRECISION :: L1,L2
		INTEGER          :: N_SPLINE
		INTEGER          :: I
		DOUBLE PRECISION :: EPS
		DOUBLE PRECISION,ALLOCATABLE :: A(:),B(:),C(:),D(:),X(:),Y(:)

		IF(NR_S<GRID_MIN) THEN
			PRINT*,'ERROR,SHORT RANGE NOT LONG ENOUGH !'
		ENDIF
		N_SPLINE=NR_S-GRID_MIN+1
                PRINT*,'N_SPLINE IS',N_SPLINE
		ALLOCATE(A(N_SPLINE))
		ALLOCATE(B(N_SPLINE))
		ALLOCATE(C(N_SPLINE))
		ALLOCATE(D(N_SPLINE))
		ALLOCATE(X(N_SPLINE))
		ALLOCATE(Y(N_SPLINE))
		X(:)=RTOX(GRID_MIN:NR_S,1)
		Y(:)=POT_ENV(GRID_MIN:NR_S)
		!PRINT*,'THE MIN POINT IS ',X(1)
		!CALL SPLINE(N_SPLINE,X,Y,A,B,C,D)
		
		EPS = 1.0D-13
		S1  = 0.0D0
		DO I=1,GRID_MIN
			IF(I==1) THEN
				S1=0.0D0
				RTOX(I,2)=S1
			ELSE
				S1=S1+DSQRT(2.0D0*FMASS*(EMAX-POT_ENV(I)))/PMAX_X*DR_S
				RTOX(I,2)=S1
			ENDIF
               ENDDO
               PRINT*,'BEFORE THE GRID_MIN IS DONE'
		!S1=0.D0
		S2   = 0.0D0
		SUM1 = 0.0D0
		SUM2 = 0.0D0
		SUM3 = 0.0D0
		SUM4 = 0.0D0
		INTERGRATION=0.0D0
		DO I=GRID_MIN+1,NRMAX
                        !PRINT*,'I=',I
			L1=RTOX(I-1,1)
			L2=RTOX(I,1)
			CALL FSIMP_CI(L1,L2,EPS,INTERGRATION)
			IF(POT_ENV(I)<-1.0D-5)      THEN
				SUM1=SUM1+INTERGRATION
			ELSEIF(POT_ENV(I)<-1.0D-10) THEN
				SUM2=SUM2+INTERGRATION
			ELSEIF(POT_ENV(I)<-1.0D-15) THEN
				SUM3=SUM3+INTERGRATION
			ELSE
				SUM4=SUM4+INTERGRATION
			ENDIF
			S2=SUM4+SUM3+SUM2+SUM1
			RTOX(I,2)=S1+S2			
        ENDDO
        PRINT*,'S1=',S1,'S2=',S2
		XLENGTH=S1+S2


		RETURN
	END SUBROUTINE GET_RTOX

	!--------------------------------------------------------------------------------
	  FUNCTION FR_SPLINE(Xab,A_r,B_r,nj,rj,asp0,bsp0,csp0,dsp0) 
	  USE MOD_SETPOT
	   implicit none	 
	   integer          :: nj,j,isp  
	   double precision :: eps0=1.0d-15
	   double precision :: FR_SPLINE,Xab,dd,A_r,B_r
       double precision :: asp0(nj),bsp0(nj),csp0(nj),dsp0(nj),rj(nj)
       if(Xab<A_r.or.Xab>B_r+eps0) then
			   print*,'Xab=',Xab,'is out of the range <FR>'
			   stop
	   endif
	   do j=1,nj     
		    if(rj(j)>Xab) then
		       isp=j-1
			   dd = Xab - rj(isp)
!			   print*,Emax-(asp0(isp) + bsp0(isp)*dd + csp0(isp)*dd**2 + dsp0(isp)*dd**3)
		       FR_SPLINE = DSQRT(2.0d0*fmass*(Emax-(asp0(isp) + bsp0(isp)*dd + csp0(isp)*dd**2 + dsp0(isp)*dd**3)))/pmax_x   ! cubic interpolation           ֵ
			   exit
			 ELSEIF(RJ(NJ)-XAB<1.0D-20) THEN
			   ISP=NJ-1
			   dd = Xab - rj(isp)
		       FR_SPLINE = DSQRT(2.0d0*fmass*(Emax-(asp0(isp) + bsp0(isp)*dd + csp0(isp)*dd**2 + dsp0(isp)*dd**3)))/pmax_x
		    endif
	   enddo
	END FUNCTION FR_SPLINE

	SUBROUTINE FSIMP_SPLINE(A,B,EPS,T,nj,rj,asp0,bsp0,csp0,dsp0)
	USE MOD_SETPOT
	IMPLICIT NONE
	INTEGER          :: N,K,nj
	DOUBLE PRECISION :: A,B,T,H,T1,S1,P,X,T2,S2,EPS
	DOUBLE PRECISION :: rj(nj),asp0(nj),bsp0(nj),csp0(nj),dsp0(nj)
	REAL*8,EXTERNAL :: FR_SPLINE
	N=1
	H=B-A
	T1=H*(FR_SPLINE(A,A,B,nj,rj,asp0,bsp0,csp0,dsp0)+FR_SPLINE(B,A,B,nj,rj,asp0,bsp0,csp0,dsp0))/2.0
	S1=T1
10	P=0.0
	DO 20 K=0,N-1
	  X=A+(K+0.5)*H
	  P=P+FR_SPLINE(X,A,B,nj,rj,asp0,bsp0,csp0,dsp0)
20	CONTINUE
	T2=(T1+H*P)/2.0
	S2=(4*T2-T1)/3.0
	IF(ABS(S2-S1).GE.EPS) THEN
	  T1=T2
	  N=N+N
	  H=H/2.0
	  S1=S2
	  GOTO 10
	END IF
	T=S2
	RETURN
	END SUBROUTINE FSIMP_SPLINE
!--------------------------------------------------------------------------------
	FUNCTION FR_CI(R) 
	USE MOD_SETPOT
!	   USE MOD_SYSPARAMETER,ONLY:EMAX,C3_1A,C6_1A,C8_1A,C3_3B,C6_3B,C8_3B,EXCHEX,EXCHA, &
!								 VEX_1A,VEX_3B,C6_1X,C8_1X,TE1X
	   IMPLICIT NONE	
	   DOUBLE PRECISION :: R
	   INTEGER          :: NS
	   DOUBLE PRECISION :: FR_CI
	   DOUBLE PRECISION :: YR
	   INTEGER          :: I
	   REAL*8,EXTERNAL:: POTENTIAL
	   
           YR=POTENTIAL(R)
	   FR_CI = DSQRT(2.0D0*FMASS*(EMAX-YR))/PMAX_X
	   RETURN
	END FUNCTION FR_CI
	SUBROUTINE FSIMP_CI(A,B,EPS,T)
	IMPLICIT NONE
	INTEGER          :: NS
	DOUBLE PRECISION :: A,B,EPS,T
	DOUBLE PRECISION :: H,T1,S1,P,X,T2,S2
	INTEGER          :: K,N
	REAL*8,EXTERNAL :: FR_CI
	N=1
	H=B-A
	T1=H*(FR_CI(A)+FR_CI(B))/2.0
	S1=T1
10	P=0.0
	DO 20 K=0,N-1
	  X=A+(K+0.5)*H
	  P=P+FR_CI(X)
20	CONTINUE
	T2=(T1+H*P)/2.0
	S2=(4*T2-T1)/3.0
	IF(ABS(S2-S1).GE.EPS) THEN
	  T1=T2
	  N=N+N
	  H=H/2.0
	  S1=S2
	  GOTO 10
	END IF
	T=S2
	RETURN
	END SUBROUTINE FSIMP_CI


!--------------------------------------------------------------------------------
	SUBROUTINE GET_POT_X(V_MIN,GRID_MIN,POT_ENV)
	USE MOD_SETPOT
    USE INTEP
		IMPLICIT NONE
		INTEGER          :: I
		INTEGER          :: GRID_MIN
		INTEGER          :: GRID_MIN_X
		DOUBLE PRECISION :: V_MIN
		DOUBLE PRECISION :: POT_ENV(NRMAX)
		DOUBLE PRECISION :: POT_ENV_X(0:NXMAX)
		DX = XLENGTH/NXMAX
		DO I=0,NXMAX
			X_GRID(I) = DBLE(I)*DX
		ENDDO
		POT_X(0)=POT_R(1)
		POT_X(NXMAX)=POT_R(NRMAX)
		
			POT_X(1:NXMAX-1)=GETINTERPOLATION(NRMAX,RTOX(:,2),POT_R(:),NXMAX-1,X_GRID(1:NXMAX-1))
		DO I=0,NXMAX
			IF(X_GRID(I)>RTOX(GRID_MIN,2)) THEN
				GRID_MIN_X=I-1
				EXIT
			ENDIF
		ENDDO
		DO I=0,GRID_MIN_X
			JX(I)=PMAX_X/DSQRT(2.0D0*FMASS*(EMAX-V_MIN))
		ENDDO
	
		JX(GRID_MIN_X+1:NXMAX)=GETINTERPOLATION(NRMAX-GRID_MIN+1,RTOX(GRID_MIN:NRMAX,2),JX_XTOR(GRID_MIN:NRMAX),NXMAX-GRID_MIN_X,X_GRID(GRID_MIN_X+1:NXMAX))
		RETURN
	    END SUBROUTINE GET_POT_X


		SUBROUTINE  GET_TIJ
		USE MOD_SETPOT
		IMPLICIT NONE
		INTEGER          :: I,J,K,N
		DOUBLE PRECISION :: AI(0:NXMAX)
		DOUBLE PRECISION :: MATRIXD(0:NXMAX,1:NXMAX-1)
		DOUBLE PRECISION :: PI
		DOUBLE PRECISION :: CON

		PI=DATAN(1.0D0)*4.0D0
		AI=1.0D0
		AI(0)=1.0D0/DSQRT(2.0D0)
		AI(NXMAX)=1.0D0/DSQRT(2.0D0)
		DO I=0,NXMAX
		   DO J=1,NXMAX-1
		      IF(I==J) THEN
			     MATRIXD(I,J)=-AI(I)*0.5D0*DCOTAN(PI*DBLE(I)/DBLE(NXMAX))
				ELSE
				 MATRIXD(I,J)=-AI(I)*0.5D0*(-1)**(I+J)*(DCOTAN(PI*DBLE(I+J)/(2.0D0*DBLE(NXMAX)))-DCOTAN(PI*DBLE(I-J)/(2.0D0*DBLE(NXMAX))))
			  ENDIF
			ENDDO
		ENDDO
   
		N=NXMAX-1
		CON=PI*PI/2.0D0/FMASS/XLENGTH/XLENGTH
		MATRIXT=0.0D0
		DO I=1,N
		   DO J=1,N
		      DO K=0,NXMAX
			     MATRIXT(I,J)=MATRIXT(I,J)+CON/DSQRT(JX(I))*MATRIXD(K,I)/JX(K)*MATRIXD(K,J)/DSQRT(JX(J))
			  ENDDO
		   ENDDO
		ENDDO
		RETURN
	    END

	SUBROUTINE GET_EIGEN(state)
        use system_parameter, only: rvdw,EvdW,AU2CMinv,MHz2AU
	USE MOD_SETPOT
        USE INTEP
!    USE SYSTEM_PARAMETER,ONLY: N_MFG,WAVE_MFG,E_VDW,MFG_E2B,POTENTIAL_ID,R_VDW
!    USE SP_PARAMETER,ONLY:RANK_S,RANK_T,RANK_B
!	USE UNIT_TRANSFORM,ONLY: MHZ2AU,AU2CM,AU2A
		IMPLICIT NONE
        REAL*8,PARAMETER :: PI=DACOS(-1.D0)
        character(len=*), intent(in) :: state
        INTEGER,PARAMETER :: NB_MAX = 5
        INTEGER,PARAMETER :: NC_MAX = 5
		REAL*8 R_LARGE
        INTEGER :: RANK_MFG
		INTEGER          :: I,J,K,M
		INTEGER          :: N1,COUNT_MFG,LWORK,LIWORK,INFO
        INTEGER,ALLOCATABLE :: EV_ID(:)
		DOUBLE PRECISION,ALLOCATABLE :: HAMILTON(:,:)
        DOUBLE PRECISION,ALLOCATABLE :: FV1(:),FV2(:)
		REAL*8,ALLOCATABLE :: WAV_X_R1(:),WAV_X_R2(:),WAV_X_R3(:),WAV_X_R(:,:),WAVe_R(:,:)
		REAL*8,ALLOCATABLE :: X_WAV(:),NORM_MFG(:), WORK(:)
        INTEGER,ALLOCATABLE ::IWORK(:)
		REAL*8 DXX,SJX,SUM_WAVE,SUM_WAVE1,SUM_WAVE2,SUM_WAVE3,PSI_PRE,PSI,R,SUM,T1,T2
        INTEGER IERR_RS
	INTEGER COUNT,NB
	INTEGER I_MFG_BOUND

      print*,'state=',state

      R_BEG=0.06D0 ! in vdW unit  
        
        
        CALL CPU_TIME(T1)
		CALL SETPOT		
        
		CALL GET_TIJ
		M = NXMAX-1
		N1= M


		ALLOCATE(EIGENVAL(N1))
		ALLOCATE(EIGENVEC(N1,N1))
		ALLOCATE(HAMILTON(N1,N1))
        LWORK=2*N1*N1+6*N1+1
        LIWORK=5*N1+3
        ALLOCATE(WORK(LWORK))
        ALLOCATE(IWORK(LIWORK))
        
        
		HAMILTON=MATRIXT
		DO I=1,N1
			HAMILTON(I,I)=HAMILTON(I,I)+POT_X(I)
        ENDDO
        CALL CPU_TIME(T2)
        
        PRINT*,'HAMILTONIAN IS OBTAIN'

        CALL DSYEVD('V','U',N1,HAMILTON,N1,EIGENVAL,WORK,LWORK,IWORK,LIWORK,INFO)
        EIGENVEC=HAMILTON
        !CALL RS(N1,N1,HAMILTON,EIGENVAL,1,EIGENVEC,FV1,FV2,IERR_RS)
         count=0
         DO I=1,N1
         IF (EIGENVAL(I).LT.0.D0) THEN
          !       PRINT*,'EB in Evdw =',EIGENVAL(I)
          !       PRINT*,'EB in cm^-1 =',EIGENVAL(I)*EvdW*AU2CMinv
                  print*,'EB in MHz=',EIGENVAL(I)*EvdW/MHz2AU
          count=count+1
         END IF
         END DO
         NB=COUNT
         print*,'the number of bound state is',count
         ALLOCATE(WAVE_R(N1,NB+2))
         WAVE_R(:,1)=GETINTERPOLATION(NRMAX,RTOX(:,2),RTOX(:,1),N1,X_GRID(1))
         DO J=1,NB+1
         DO I=1,N1
         WAVE_R(I,J+1)=EIGENVEC(I,J)/DSQRT(JX(I))/DSQRT(DX)
         END DO
         END DO

         OPEN(888,FILE=trim(state)//'Eb.txt')
         write(888,*) 'Eb(Evdw) ','Eb(MHz)      ','Eb(cm^-1)    '        
         DO I=1,NB+1
         WRITE(888,*) EIGENVAL(I),EIGENVAL(I)*EvdW/MHz2AU, EIGENVAL(I)*EvdW*AU2CMinv
         END DO
         CLOSE(888)
         open(999,FILE=trim(state)//'wave.txt')
         do i=1,n1
         write(999,'(100g20.10)') Wave_r(i,:)
         end do
         close(999)
!stop


         DEALLOCATE(WAVE_R)


		
		DEALLOCATE(HAMILTON)
        DEALLOCATE(WORK)
        DEALLOCATE(IWORK)
        
        
        DEALLOCATE(EIGENVAL)
		DEALLOCATE(EIGENVEC)
		RETURN
    END SUBROUTINE GET_EIGEN
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	SUBROUTINE get_eigen_trap(state,nu_kHz)
        use system_parameter, only: rvdw,EvdW,AU2CMinv,MHz2AU
	USE MOD_SETPOT
        USE INTEP
!    USE SYSTEM_PARAMETER,ONLY: N_MFG,WAVE_MFG,E_VDW,MFG_E2B,POTENTIAL_ID,R_VDW
!    USE SP_PARAMETER,ONLY:RANK_S,RANK_T,RANK_B
!	USE UNIT_TRANSFORM,ONLY: MHZ2AU,AU2CM,AU2A
		IMPLICIT NONE
        REAL*8,PARAMETER :: PI=DACOS(-1.D0)
        character(len=*), intent(in) :: state
        INTEGER,PARAMETER :: NB_MAX = 5
        INTEGER,PARAMETER :: NC_MAX = 5
		REAL*8 R_LARGE
        INTEGER :: RANK_MFG
		INTEGER          :: I,J,K,M
		INTEGER          :: N1,COUNT_MFG,LWORK,LIWORK,INFO
        INTEGER,ALLOCATABLE :: EV_ID(:)
		DOUBLE PRECISION,ALLOCATABLE :: HAMILTON(:,:)
        DOUBLE PRECISION,ALLOCATABLE :: FV1(:),FV2(:)
		REAL*8,ALLOCATABLE :: WAV_X_R1(:),WAV_X_R2(:),WAV_X_R3(:),WAV_X_R(:,:),WAVe_R(:,:)
		REAL*8,ALLOCATABLE :: X_WAV(:),NORM_MFG(:), WORK(:)
        INTEGER,ALLOCATABLE ::IWORK(:)
		REAL*8 DXX,SJX,SUM_WAVE,SUM_WAVE1,SUM_WAVE2,SUM_WAVE3,PSI_PRE,PSI,R,SUM,T1,T2
        INTEGER IERR_RS
        real*8 nu_khz
	INTEGER COUNT,NB
	INTEGER I_MFG_BOUND

      print*,'state=',state
      nu=nu_khz/1.d3*MHz2AU/EvdW
      print*,'\omega*\hbar in vdW=',nu 


      R_BEG=0.06D0 ! in vdW unit  
        
        
        CALL CPU_TIME(T1)
	CALL SETPOT		
	CALL GET_TIJ
	M = NXMAX-1
	N1= M

        R_MFG=GETINTERPOLATION(NRMAX,RTOX(:,2),RTOX(:,1),N1,X_GRID(1))
        do i=1,n1
        Trap_X(i)=0.5d0*0.5d0*nu**2.d0*R_MFG(i)**2.d0
        end do
        !print*,'nu=',nu
        !stop

	ALLOCATE(EIGENVAL(N1))
	ALLOCATE(EIGENVEC(N1,N1))
	ALLOCATE(HAMILTON(N1,N1))
        LWORK=2*N1*N1+6*N1+1
        LIWORK=5*N1+3
        ALLOCATE(WORK(LWORK))
        ALLOCATE(IWORK(LIWORK))
        
        
		HAMILTON=MATRIXT
		DO I=1,N1
			HAMILTON(I,I)=HAMILTON(I,I)+POT_X(I)+Trap_X(I)
        ENDDO
        CALL CPU_TIME(T2)
        
        PRINT*,'HAMILTONIAN IS OBTAIN'

        CALL DSYEVD('V','U',N1,HAMILTON,N1,EIGENVAL,WORK,LWORK,IWORK,LIWORK,INFO)
        EIGENVEC=HAMILTON
        !CALL RS(N1,N1,HAMILTON,EIGENVAL,1,EIGENVEC,FV1,FV2,IERR_RS)
         count=0
         DO I=1,N1
         IF (EIGENVAL(I).LT.0.D0) THEN
          !       PRINT*,'EB in Evdw =',EIGENVAL(I)
          !       PRINT*,'EB in cm^-1 =',EIGENVAL(I)*EvdW*AU2CMinv
                  print*,'EB in MHz=',EIGENVAL(I)*EvdW/MHz2AU
          count=count+1
         END IF
         ! print*,'EV in MHz=',EIGENVAL(I)*EvdW/MHz2AU
         END DO
         NB=COUNT
         print*,'the number of bound state is',count
         ALLOCATE(WAVE_R(N1,NB+2))
         WAVE_R(:,1)=GETINTERPOLATION(NRMAX,RTOX(:,2),RTOX(:,1),N1,X_GRID(1))
         DO J=1,NB+1
         DO I=1,N1
         WAVE_R(I,J+1)=EIGENVEC(I,J)/DSQRT(JX(I))/DSQRT(DX)
         END DO
         END DO

         OPEN(888,FILE=trim(state)//'Eb.txt')
         write(888,*) 'Eb(Evdw) ','Eb(MHz)      ','Eb(cm^-1)    '        
         DO I=1,NB+1
         WRITE(888,*) EIGENVAL(I),EIGENVAL(I)*EvdW/MHz2AU, EIGENVAL(I)*EvdW*AU2CMinv
         END DO
         CLOSE(888)
         open(999,FILE=trim(state)//'wave.txt')
         do i=1,n1
         write(999,'(100g20.10)') Wave_r(i,:)
         end do
         close(999)
!stop


         DEALLOCATE(WAVE_R)


		
		DEALLOCATE(HAMILTON)
        DEALLOCATE(WORK)
        DEALLOCATE(IWORK)
        
        
        DEALLOCATE(EIGENVAL)
		DEALLOCATE(EIGENVEC)
		RETURN
    END SUBROUTINE GET_EIGEN_TRAP
 !------------------------------------------------------------------------------------------------------- 
   function potential(r)
          use mod_setpot, only: lammfg
          implicit none
         real*8 potential,r
!--------------------------------------------------------------
! in an abitary unit
         potential=(lammfg**6.d0/r**6.d0-1.d0)*16.d0/r**6.d0

!------------------------------------------------------------

        return
        end
   subroutine getMFGParameter(c6,c12)
      use MOD_setPOT, only: c6mfg,lammfg,potminr,potmin
      use system_parameter,only: mass, rvdw,Evdw,AU2CMinv
      implicit none
      real*8 c6,lam,c12
      c6mfg=c6
      !lammfg=lam
      rvdw=(mass*c6mfg)**0.25d0/2
      EvdW=1.D0/mass/rvdw/rvdw
      lam=(c12/c6)**(1.d0/6.d0)/rvdw
      lammfg=lam
      print*,'rvdw=',rvdW,'EvdW=',EvdW,'lam=',lammfg
      potminr=lammfg*(2.d0)**(1.d0/6.d0)
      potmin=-4.d0/lammfg**6.d0
      print*,'vminr in a0=', potminr*rvdw,'vmin in cm^-1=',potmin*EvdW*AU2CMinv
      end     
!subroutine readSystemParameter
!       implicit none
!       open(10,file='MFG.inp',status='old')
!       read(10,*)
!       read(10,) c6g, lambdag
!       read(10,*)
!       read(10,*)
!       read(10,) c6e,lambdae
!      return
!      end 
