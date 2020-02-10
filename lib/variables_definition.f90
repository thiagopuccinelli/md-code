module variables_definition  
    
    implicit none 
    integer, parameter :: Nsteps = 15000
    integer, parameter :: Nsteps_equil = 5000
    real, parameter :: dt = 0.001d0 
    real, parameter :: t = dt * Nsteps
    integer :: step  
    integer, parameter :: npart = 1000 
    real, parameter :: mass = 1.0d0
    real, parameter :: eps = 1.0d0
    real, parameter :: sigma = 1.0d0
    real, parameter :: rho = 0.75d0
    character(10), parameter  :: res_rho = "_0.75_"
    real, parameter :: BOX = (npart/rho) ** (1./3.)
    real, parameter :: Rcutoff = 2.5d0 * sigma  ! cutoff distance 2.5*sigma
  	real, parameter :: phicutoff =  1./(Rcutoff**12) - 1.d0/(Rcutoff**6)    ! potential at cutoff
    real, parameter :: temprqs = 1.d0
    real, parameter :: gamma = 1.0d0
    real :: temperature  
    integer, parameter :: dim = 3  
    real :: r(npart,dim), f(npart,dim), v(npart,dim)
    real :: en_pot(npart), en_kin 
end module variables_definition