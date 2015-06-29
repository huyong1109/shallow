module module_para

    implicit none

    real*8,  parameter       ::    omg0 = 7.292d-5	    ! angular velocity of the earth rotation
    real*8,  parameter       ::    a    = 6371000d0	    ! radius of the earth
     													       
    integer, parameter       ::    nx  = 180		    ! to define zonal resolution
    !integer, parameter       ::    nx  = 320 ! to define zonal resolution
    integer, parameter       ::    p = nx/2		    ! to define zonal resolution
    integer, parameter       ::    n1 = nx+2			   
     													       
    integer, parameter       ::    ny = 89  		    ! to define meridional resolution
    !integer, parameter       ::    ny = 384     ! to define meridional resolution
    !integer, parameter       ::    q  = ny/2			   
    integer, parameter       ::    n  = ny+2		    ! meridional grid number
      														   
    integer, parameter       ::    nm = nx*ny	    ! total grid number

    integer, parameter       ::    tlo = 800		    ! time stepsize
    integer, parameter       ::    t0 = 0		    ! initial time
    integer, parameter       ::    t1 = 8640000		    ! final time
							    	
    integer, parameter       ::    runcase = 0      ! 0: using longitude-latitude grid ; 1: using general curvilinear grid 
    character (len = *), parameter :: grid_file = "pop_grid.nc" ! if  runcase == 1


    integer, parameter       ::    n0 = 0    ! 0: using rh waves as initial conditon (ic); 1: read ic from files
    integer, parameter       ::    nyn = 0  ! screen output 0: at each step; 1: once 12 hours
    integer, parameter       ::    thalf = 432000    ! the time for saving the intermediate result

     														   
    character*7, parameter   ::    fu='uui.dat'		    ! initial field of zonal wind for self-defined experiment
    character*7, parameter   ::    fv='vvi.dat'	   	    ! initial field of meridional wind for self-defined experiment
    character*7, parameter   ::    fh='hhi.dat'	 	    ! initial field of geopotential height for self-defined experiment
     														   
    character*6, parameter   ::    ffu='ui.dat'		    ! result of zonal wind for self-defined experiment
    character*6, parameter   ::    ffv='vi.dat'		    ! result of meridional wind for self-defined experiment
    character*6, parameter   ::    ffh='hi.dat'		    ! result of geopotential height for self-defined experiment
     														   
    real*8                   ::    pi					   
    real*8                   ::    deta,deta1,detap,detaq   

    ! grid variabls 
    real*8,  dimension(1:n1,1:n)  ::    lon,lat ! lon-lat degree
    real*8,  dimension(1:n1,1:n)  ::    dxu,dyu ! dx, dy between grid point
    real*8,  dimension(1:n1,1:n)  ::    dxr,dyr ! dx, dy between grid point
    real*8,  dimension(1:n1,1:n)  ::    hte,htn ! reverease of distance along grid box
    real*8,  dimension(1:n1,1:n)  ::    c1, s1       ! c1(j)=cos(theta(j)), s1(j)=sin(theta(j)), where theta is latitude
    real*8,  dimension(1:n1,1:n)  ::    f1, f2           ! 

end module module_para
