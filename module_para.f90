module module_para

    implicit none

    real*8,  parameter       ::    omg0 = 7.292d-5	    ! angular velocity of the earth rotation
    real*8,  parameter       ::    a    = 6371000d0	    ! radius of the earth
     													       
    integer, parameter       ::    p  = 720		    ! to define zonal resolution
    integer, parameter       ::    kn = p/2		    ! to define zonal resolution
    integer, parameter       ::    np = p+1		    ! zonal grid nmuber 
    integer, parameter       ::    n1 = np+1			   
     													       
    integer, parameter       ::    q0 = 180  		    ! to define meridional resolution
    integer, parameter       ::    nq = q0+1			   
    integer, parameter       ::    q  = q0*2			   
    integer, parameter       ::    n  = q+1		    ! meridional grid number
      														   
    integer, parameter       ::    nm = (np-1)*n	    ! total grid number

    integer, parameter       ::    tlo = 200		    ! time stepsize
    integer, parameter       ::    t0 = 0		    ! initial time
    integer, parameter       ::    t1 = 8640000		    ! final time
							    	
    integer, parameter       ::    n0 = 0		    ! 0: using rh waves as initial conditon (ic); 1: read ic from files
    integer, parameter       ::    nyn = 0		    ! screen output 0: at each step; 1: once 12 hours
    integer, parameter       ::    thalf = 432000	    ! the time for saving the intermediate result

     														   
    character*7, parameter   ::    fu='uui.dat'		    ! initial field of zonal wind for self-defined experiment
    character*7, parameter   ::    fv='vvi.dat'	   	    ! initial field of meridional wind for self-defined experiment
    character*7, parameter   ::    fh='hhi.dat'	 	    ! initial field of geopotential height for self-defined experiment
     														   
    character*6, parameter   ::    ffu='ui.dat'		    ! result of zonal wind for self-defined experiment
    character*6, parameter   ::    ffv='vi.dat'		    ! result of meridional wind for self-defined experiment
    character*6, parameter   ::    ffh='hi.dat'		    ! result of geopotential height for self-defined experiment
     														   
    real*8                   ::    pi					   
    real*8                   ::    deta,deta1,deta2,detb   

    real*8,  dimension(1:n)  ::    c1, s1		           ! c1(j)=cos(theta(j)), s1(j)=sin(theta(j)), where theta is latitude
    real*8,  dimension(1:n)  ::    c2, c11, c12, c13, c14  ! 
    real*8,  dimension(1:n)  ::    f1, f2		           ! 
end module module_para
