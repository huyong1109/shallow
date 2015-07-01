subroutine llgrid 
    ! for longitude-latitude coordinates
    use module_para

    implicit none

    real*8  :: ai, aj    ! working variables
    integer :: i, j         ! working variable
    
    pi=datan(1d0)*4
    deta1=2*pi/nx
    deta=pi/90d0
    !deta=pi/(ny+1)
    detap=2*pi/nx*a
    !detaq=pi/(ny+1)*a
    detaq=pi/90d0*a

    do i=1,n1
        if ((i-2)*deta1 <= pi ) then 
            lon(i,:)    = (i-2)*deta1
        else 
            lon(i,:)    = (i-2)*deta1-2*pi
        end if 
    end do

    !do j=2,n-1
    !        ai            =  (j-1)*deta-pi/2.0d0
    !        lat(:,j)      = ai
    !        c1(:,j)       = dcos(ai)
    !        s1(:,j)       = dsin(ai)
    !end do

    !ai =0.5*deta-pi/2.0d0
    !lat(:,1) = ai
    !c1(:,1)=dcos(ai)/4.
    !    
    !ai=(ny+0.5)*deta-pi/2.0d0
    !c1(:,n)=dcos(ai)/4.
    !  
    !s1(:,1)=-1.0d0
    !s1(:,n)=1.0d0
    
    do j=1,n
            ai            =  (j-1)*deta-pi/2.0d0 + 5.0/90d0*pi
            lat(:,j)      = ai
            c1(:,j)       = dcos(ai)
            s1(:,j)       = dsin(ai)
    end do
    write(*,*) "lat(1,:)"
    write(*,*) lat(1,:)

    
    do j=1,n
        ai=s1(1,j)
        aj=1/c1(2,j)/a
        f1(:,j)=2.0d0*omg0*ai !coriolis 
        f2(:,j)=s1(:,j)*aj
    end do


    dxu(:,:)=detap*c1(:,:)
    dyu(:,:)=detaq*c1(:,:)
    dxr(:,:)=1.0d0/dxu(:,:)
    dyr(:,:)=1.0d0/dyu(:,:)
    hte(:,:)=dxr(:,:)
    htn(:,:)=dyr(:,:)

    write(*,*) "n =" ,n ,"n1= ", n1
    write(*,*) "c1(10,:)"
    write(*,*) c1(10,:)
    write(*,*) "s1(1,:)"
    write(*,*) s1(10,:)
    write(*,*) "f1(1,:)"
    write(*,*) f1(10,:)
    write(*,*) "f2(1,:)"
    write(*,*) f2(10,:)
    write(*,*) "dxr"
    write(*,*) dxr(10,:)
    write(*,*) "dyr"
    write(*,*) dyr(10,:)

    return
end subroutine llgrid

subroutine  gcgrid
    ! for general curvilinear grid 
    ! read in from POP output file
    use module_para
    use module_array
    use netcdf 

    implicit none

    real*8  :: ai, aj    ! working variables
    integer :: i, j         ! working variable
    ! for POP general grid 
    integer :: ncid, varid,stat

    stat = nf90_open(grid_file,nf90_write,ncid)
    if (stat /= nf90_noerr ) then
        print *,trim(nf90_strerror(stat))
        stop 
    end if 

    call readnc(ncid, "ULONG", lon(2:nx+1,2:ny+1),nx,ny)
    call readnc(ncid, "ULAT", lat(2:nx+1,2:ny+1),nx,ny)
    call readnc(ncid, "DXU", dxu(2:nx+1,2:ny+1),nx,ny)
    call readnc(ncid, "DYU", dyu(2:nx+1,2:ny+1),nx,ny)
    call readnc(ncid, "HTE", hte(2:nx+1,2:ny+1),nx,ny)
    call readnc(ncid, "HTN", htn(2:nx+1,2:ny+1),nx,ny)

    c1(2:nx+1,2:ny+1) = dcos(lat(2:nx+1,2:ny+1))/4.0
    s1(2:nx+1,2:ny+1) = dcos(lat(2:nx+1,2:ny+1))/4.0
    f1(2:nx+1,2:ny+1) = 2.0d0*omg0*s1(2:nx+1,2:ny+1)
    f2(2:nx+1,2:ny+1) = s1(2:nx+1,2:ny+1)/c1(2:nx+1,2:ny+1)/a


    stat = nf90_close(ncid)
    if (stat /= nf90_noerr ) then
        print *,trim(nf90_strerror(stat))
        stop
    end if 
    !write(*,*) dxu(2:nx+1,2:ny+1)
    !write(*,*) dyu(2:nx+1,2:ny+1)
    !write(*,*) hte(2:nx+1,2:ny+1)
    !write(*,*) htn(2:nx+1,2:ny+1)

    return

end subroutine gcgrid

subroutine readnc(ncid, field, indata,nx,ny)
    use netcdf 
    implicit none
    integer, intent(in) :: ncid
    integer, intent(in) :: nx,ny
    character (len=10), intent(in) :: field
    real*8, dimension(nx,ny), intent(inout) :: indata
    !local 
    integer :: varid,stat

    stat = nf90_inq_varid(ncid,field,varid)
    if (stat /= nf90_noerr ) then
        print *,trim(nf90_strerror(stat))
        stop "Stopped in read "
    end if 
    stat = nf90_get_var(ncid,varid,indata)
    if (stat /= nf90_noerr ) then
        print *,trim(nf90_strerror(stat))
        stop "Stopped in read "
    end if 
end subroutine readnc
