subroutine difuh(wu,wv,du,dh,h)
    use module_para
    implicit none
    
    real*8,dimension(1:n1,1:n)  :: dh,h,hy
    real*8,dimension(1:n1,1:n)  :: wu,du,u
    real*8,dimension(1:n1,1:n)  :: wv,v,vv
    real*8                      :: hyn,hys
    real*8                      :: ff
    integer                     :: i,j
    
    do j=2,n-1
	do i=2,np
	    u (i,j)=wu(i,j)/h(i,j)
	    v (i,j)=wv(i,j)/h(i,j)
	    vv(i,j)=v (i,j)*c1(i,j)
	enddo
    enddo
    
    call advct(u,vv,wu,du)
       
    do j=2,n-1
	do i=2,np
	    ff=f1(i,j)+u(i,j)*f2(i,j)
	    du(i,j)=du(i,j)-ff*wv(i,j)
	    
	    hy(i,j)=wv(i,j)*h(i,j)*c1(i,j)
	end do
    end do
    
    do i=2,np
	    du(i,1)=0.0
	    du(i,n)=0.0
	    
	    hy(i,1)=0.0
	    hy(i,n)=0.0
    enddo
       
    do j=2,n-1
	do i=2,np
	    dh(i,j)=(hy(i,j+1)-hy(i,j-1))*0.5*dyr(i,j)
	end do
    end do
    
    hys=0.0
    hyn=0.0
    do i=2,np
	hys=hys+hy(i,2)
	hyn=hyn-hy(i,n-1)
    enddo
    hys=hys*0.5*dyr(1,1)/(np-1)
    hyn=hyn*0.5*dyr(1,n)/(np-1)
    
    do i=2,np
	dh(i,1)=hys
	dh(i,n)=hyn
    enddo
    
    return
end subroutine difuh
    
subroutine difv(wu,wv,wh,dv,h)
    use module_para
    implicit none
    
    real*8,dimension(1:n1,1:n)  :: wh,h,hh
    real*8,dimension(1:n1,1:n)  :: wu,u
    real*8,dimension(1:n1,1:n)  :: wv,dv,v,vv
    real*8                      :: ff
    integer                     :: i,j
    
    do j=2,n-1
	do i=2,np
	    hh(i,j)=h(i,j)*c1(i,j)
	    u(i,j)=wu(i,j)/h(i,j)
	    v(i,j)=wv(i,j)/h(i,j)
	    vv(i,j)=v(i,j)*c1(i,j)
	enddo
    enddo
    
    call advct(u,vv,wv,dv)
       
    do j=2,n-1
	do i=2,np
	    ff=f1(i,j)+u(i,j)*f2(i,j)
	    dv(i,j)=dv(i,j)+hh(i,j)*(wh(i,j+1)-wh(i,j-1))*0.5*dyr(i,j)+ff*wu(i,j)
	end do
    end do
    
    do i=2,np
	dv(i,1)=0.0
	dv(i,n)=0.0
    enddo
    
    return
end subroutine difv
    
subroutine advct(u,v,f,df)
    use module_para
    implicit none
    
    real*8,dimension(1:n1,1:n) :: f,df
    real*8,dimension(1:n1,1:n) :: u,v
    real*8,dimension(1:n1,1:n) :: su,sv
    real*8                     :: dx,dy
    integer                    :: i,j
    
    do j=2,n-1
	do i=2,np
	    su(i,j)=f(i,j)*u(i,j)
	    sv(i,j)=f(i,j)*v(i,j)
	enddo
	su(1,j)=su(np,j)
	su(n1,j)=su(2,j)
	
	f(1,j)=f(np,j)
	f(n1,j)=f(2,j)
    enddo
    
    do i=2,np
	sv(i,1)=0
	sv(i,n)=0
	df(i,1)=0.0
	df(i,n)=0.0
    enddo
    
    do j=2,n-1
	do i=2,np
	    dx=u(i,j)*(f(i+1,j)-f(i-1,j))+su(i+1,j)-su(i-1,j)
	    dy=v(i,j)*(f(i,j+1)-f(i,j-1))+sv(i,j+1)-sv(i,j-1)
	    df(i,j)=0.25*(dx*dxr(i,j)+dy*dyr(i,j))
	enddo
    enddo
    return
end subroutine advct
