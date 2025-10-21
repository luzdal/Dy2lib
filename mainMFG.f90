program MFG
        use system_parameter
implicit none
call readSystemParameter
call getMFGParameter(c6g,c12g)
call get_eigen_trap('ground',2.d2)
!the second argument of get_eigen_trap is trap frequency in KHz
call getMFGParameter(c6e,c12e)
call get_eigen_trap('excited',2.d2)
end
