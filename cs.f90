subroutine cs
    use module_para
    use module_array

    implicit none

    real*8  :: ai, aj    ! working variables
    integer :: i, j         ! working variable

    pi=datan(1d0)*4
    deta1=2*pi/p
    deta=pi/q
    detap=2*pi/p*a
    detaq=pi/q*a

    ! for longitude-latitude coordinates
    do i=1,n1
        if ((i-2)*deta1 <= pi ) then 
            lon(i,:)    = (i-2)*deta1
        else 
            lon(i,:)    = (i-2)*deta1-2*pi
        end if 
    end do

    do j=2,n-1
            ai            =  (j-1)*deta-pi/2
            lat(:,j)      = ai
            c1(:,j)       = dcos(ai)
            s1(:,j)       = dsin(ai)
    end do

    ai =0.5*deta-pi/2
    lat(:,1) = ai
    c1(:,1)=dcos(ai)/4.
        
    ai=(n-1.5)*deta-pi/2
    c1(:,n)=dcos(ai)/4.
      
    s1(:,1)=-1.0d0
    s1(:,n)=1.0d0

    do j=1,n
        ai=s1(1,j)
        aj=c1(1,j)*a
        f1(:,j)=2.0d0*omg0*ai !coriolis 
        f2(:,j)=s1(:,j)/aj
    end do

    dxu(:,:)=detap*c1(:,:)
    dyu(:,:)=detaq*c1(:,:)
    dxr(:,:)=1/dxu(:,:)
    dyr(:,:)=1/dyu(:,:)


    ! for general curvilinear grid 
    ! read in from POP output file
     
    return
end subroutine cs
