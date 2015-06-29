subroutine euler(dt,iter)
    use module_para
    use module_array
    implicit none

    integer i,i1,i2,j,k,iter
     
    real*8,dimension(1:n1,1:n)  :: tu,tv,th,h 
    real*8,dimension(1:n1,1:n)  :: du,dv,dh 
    real*8,dimension(1:n1)      :: a1,a2,rh,ru 
    real*8,dimension(1:p)       :: fm,fp,f0,gm,gp,g0,rf,rg 
    real*8                      :: ai,aj,dt,dt2,en,en0,den
    real*8,external             :: inner

    dt2=dt*0.5d0

    do j=1,n
	do i=2,nx+1
	    tu(i,j)=wu(i,j)
	    tv(i,j)=wv(i,j)
	    th(i,j)=wh(i,j)
	end do
    end do

    en0=inner(wu,wv,wh,wu,wv,wh)

    do k=1,1000

	do j=1,ny+2
	    do i=2,nx+1
		h(i,j)=dsqrt(th(i,j))
	    end do
	end do

	call difuh(tu,tv,du,dh,h)

	do j=1,ny+2
	    do i=2,nx+1
		tu(i,j)=wu(i,j)-dt2*du(i,j)
		th(i,j)=wh(i,j)-dt2*dh(i,j)
	    end do
	end do

	do j=2,ny+1

	    do i=2,nx+1
	    ai=dt2*0.5*dxr(i,j)
		aj=ai*h(i,j)
		a1(i)=aj
		ru(i)=tu(i,j)*aj

		a2(i)=aj*aj
	    end do

	    ru(1)=ru(nx+1)
	    ru(n1)=ru(2)

	    do i=2,nx+1
		rh(i)=th(i,j)-ru(i+1)+ru(i-1)
	    enddo

	    do i=1,p
		i1=i*2
		i2=i1+1

		fp(i)=-a2(i2)
		rf(i)=rh(i1)

		gm(i)=-a2(i1)
		rg(i)=rh(i2)
	    enddo

	    do i=2,p
		fm(i)=fp(i-1)
	    enddo
	    fm(1)=fp(p)

	    do i=1,p-1
		gp(i)=gm(i+1)
	    enddo
	    gp(p)=gm(1)

	    do i=1,p
		f0(i)=1.0-fm(i)-fp(i)
		g0(i)=1.0-gm(i)-gp(i)
	    enddo

	    call lu0(fm,f0,fp,rf,p)
	    call lu0(gm,g0,gp,rg,p)

	    do i=1,p
		i1=i*2
		i2=i1+1

		th(i1,j)=rf(i)
		th(i2,j)=rg(i)
	    enddo

	    th(1,j)=th(nx+1,j)
	    th(n1,j)=th(2,j)

	    do i=2,nx+1
		tu(i,j)=tu(i,j)-a1(i)*(th(i+1,j)-th(i-1,j))
	    enddo
	end do
	call difv(tu,tv,th,dv,h)

	do j=1,ny+2
	    do i=2,nx+1
		tv(i,j)=wv(i,j)-dt2*dv(i,j)
	    end do
	end do

	en=inner(tu,tv,th,tu,tv,th)

	den=dabs(en-en0)*2.0/(en+en0)
	en0=en
	if (den.lt.1.0d-15) goto 10

    end do

    10   continue
    iter=k

    do j=1,ny+2
	do i=2,nx+1
	    wu(i,j)=tu(i,j)*2.0d0-wu(i,j)
	    wv(i,j)=tv(i,j)*2.0d0-wv(i,j)
	    wh(i,j)=th(i,j)*2.0d0-wh(i,j)
	end do
    end do

    return
end subroutine euler

subroutine lu0(a,b,c,r,n)
    implicit none

    integer               :: i,n
    real*8                :: ai,sn,rn
    real*8,dimension(1:n) :: a,b,c,r
    real*8,dimension(1:n) :: s,t

    s(1)=a(1)
    t(1)=c(n)

    sn=0.0
    rn=0.0

    do i=2,n-1
	ai=a(i)/b(i-1)
	b(i)=b(i)-ai*c(i-1)
	r(i)=r(i)-ai*r(i-1)
	s(i)=-ai*s(i-1)

	ai=t(i-1)/b(i-1)
	t(i)=-ai*c(i-1)
	sn  =sn-ai*s(i-1)
	rn  =rn-ai*r(i-1)
    enddo

    a(n)=a(n)+t(n-1)
    b(n)=b(n)+sn
    c(n-1)=c(n-1)+s(n-1)
    r(n)=r(n)+rn

    ai=a(n)/b(n-1)
    b(n)=b(n)-ai*c(n-1)
    r(n)=r(n)-ai*r(n-1)

    r(n)=r(n)/b(n)
    r(n-1)=(r(n-1)-c(n-1)*r(n))/b(n-1)

    do i=n-2,1,-1
	ai=r(i)-s(i)*r(n)-c(i)*r(i+1)
	r(i)=ai/b(i)
    enddo

    return
end subroutine 
