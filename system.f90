module system_parameter
        implicit none
        real*8,parameter :: AMU2AU = 1822.888486D0
        real*8,parameter :: AU2A    =    0.52917721092d0
        !real*8,parameter :: mass = 163.9291712d0*amu2au
        real*8,parameter :: MHz2AU  = 1.5198298460570D-10
        real*8,parameter :: AU2T = 3.157853405D5
        real*8,parameter :: AU2CMinv = 219474.631371
        real*8 rvdw,Evdw
        real*8 c6g,lambdag,c6e,lambdae,c12g,c12e
        real*8 mass,massamu
        !----------------------------------------------------------
        ! short-range boundary condition according to New J Phys. 17, 035015(2015)
        real*8,parameter :: PSR = 1.d0  ! [0,1] the short-range loss probability, 1 mean full loss in universal regime
        real*8,parameter :: PhaseSR = dcos(-1.d0)/2.d0 ! [0,pi] phase shift accumulated at short-range
        !----------------------------------------------------------
end module system_parameter

subroutine readSystemParameter
       use system_parameter
       implicit none
       open(10,file='MFG.inp',status='old')
       read(10,*)
       read(10,*) massamu
       read(10,*) 
       read(10,*)
       read(10,*) c6g, c12g
       read(10,*)
       read(10,*)
       read(10,*) c6e,c12e
       print*,'c6g=',c6g,'c6e=',c6e,'c12g=',c12g,'c12e=',c12e,'massamu=',massamu
       mass=massamu*amu2au
      !return
      end
