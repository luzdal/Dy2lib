program MFG
        use system_parameter
implicit none
real*8 E,scatlen
call readSystemParameter
call getMFGParameter(c6g,c12g)
E=1.d-4 ! in EvdW
call  DRIVE_CC_SINGLE(E,scatlen)
! scatlen in a0
!print*,'scatlength is',scatlen
!stop
call get_eigen_trap('ground',2.d2)
!the second argument of get_eigen_trap is trap frequency in KHz
call getMFGParameter(c6e,c12e)
call get_eigen_trap('excited',2.d2)
end
