!*****************************************************************************!
!			The Barotropic Model on Arakawa A-grid		      ! 		  
!                           by Bin Wang					      !
!	        with the implicit scheme of energy conservation               !
!*****************************************************************************!	
program main
    use module_para
    use module_array

    implicit none

    real*8                      :: tener,tener0		          ! total energy at tn and t0, respectively
    real*8                      :: tmass,tmass0			  ! total mass at tn and t0, respectively
    real*4, dimension(1:n1,1:n) :: pu,pv,ph                       ! for grads saving
    real*8                      :: ai,dt                          ! working variables
    integer                     :: tlp,iter,irecu,irecv,irech     ! working variables
    integer                     :: i,j,iws,nt,nw,iwr,tt           ! working variables
    
    real*8, external            :: inner		    ! a external function to calculate inner product

    tt = t0
    iwr=(t1-tt)/tlo
    tlp=t1-tt-tlo*iwr
     
    if (tlp.gt.0) then
	iwr=iwr+1
    else
	tlp=tlo
    end if

    print *,'initial time is',tt
    print *,'final time is',t1
    print *,'time stepsize is',tlo

    if (n0.eq.0) then
	print *,'this is a rh-wvae experiment,'
	print *,'i.e., the ic is rh-wave'
    else
	print *,'this is self-defined experiment,'
	print *,'i.e., the ic is read from files'
    end if

    if (nyn.ne.0) then
	print *,'the energy... will be shown once 12 hours'
    else
	print *,'the energy... will be shown at each time step'
    end if

    !   parameter setting
    call cs

    !	initial condition
    if (n0.eq.0) then
	!      using rossby-haurwitz waves as initial condition
	call haurwitz
    else
	!      read initial conditon from your files
	open(10,file=fu)
	open(11,file=fv)
	open(12,file=fh)
	 
	read(10,130)(( u(i,j),j=1,n),i=1,n1)
	read(11,130)(( v(i,j),j=1,n),i=1,n1)
	read(12,130)((wh(i,j),j=1,n),i=1,n1)

	close(10)
	close(11)
	close(12)
    end if
     
    call opf(np,n,'uu.dat',12)
    call opf(np,n,'vv.dat',13)
    call opf(np,n,'hh.dat',14)
     
    irecu=0
    irecv=0
    irech=0
     
    do j=1,n
	do i=2,np
	    pu(i,j)=u(i,j)
	    pv(i,j)=v(i,j)
	    ph(i,j)=wh(i,j)
	    ai=dsqrt(wh(i,j))
	    wu(i,j)=u(i,j)*ai
	    wv(i,j)=v(i,j)*ai
	end do
	pu(1,j)=pu(np,j)
	pv(1,j)=pv(np,j)
	ph(1,j)=ph(np,j)

	wu(1,j)=wu(np,j)
	wu(n1,j)=wu(2,j)

	wv(1,j)=wv(np,j)
	wv(n1,j)=wv(2,j)

	wh(1,j)=wh(np,j)
	wh(n1,j)=wh(2,j)
    end do

    call wr(pu,np,n,12,irecu)
    call wr(pv,np,n,13,irecv)
    call wr(ph,np,n,14,irech)

    tener0=inner(wu,wv,wh,wu,wv,wh)

    tmass0 = 0
    do j=1,n
	do i=2,np
	    tmass0=tmass0+wh(i,j)*c1(i,j)
	end do
    end do

    print *,'the total energy is ',tener0
    print *,'the total mass is ',tmass0

    nt=23
    dt=tlo
    nw=0
    !
    print *,'the main part of this program has started'

    print *,'number of integration steps is',iwr

    do iws=1,iwr

	if (iws.eq.iwr) then
	    dt=tlp
	    tt=tt+tlp
	else
	    tt=tt+tlo
	end if

	nw=nw+1

	if (nt.eq.23) then
	    print *,'-------------------------------------------------------------------------------'
	    print *,'       the energy           the total-mass        the iteration number'
	    nt=1
	end if
	!------------------------------------------------
	!       the time integration
	!------------------------------------------------
	call euler(dt,iter)

	do j=1,n
	    do i=2,np
		ai=dsqrt(wh(i,j))
		u(i,j)=wu(i,j)/ai
		v(i,j)=wv(i,j)/ai
	    end do
	end do
	!
	if ((tt-thalf.ge.0).and.(tt-thalf.lt.dt)) then
	    !
	    do j=1,n
		do i=2,np
		    ai=dsqrt(wh(i,j))
		    pu(i,j)=u(i,j)
		    pv(i,j)=v(i,j)
		    ph(i,j)=wh(i,j)
		end do
		pu(1,j)=pu(np,j)
		pv(1,j)=pv(np,j)
		ph(1,j)=ph(np,j)
	    end do
	    !
	    call wr(pu,np,n,12,irecu)
	    call wr(pv,np,n,12,irecu)
	    call wr(ph,np,n,14,irech)
	    !
	end if
	!
	if ((nyn.eq.1).or.(int(tt/43200)*43200.eq.tt)) then
	    tener=inner(wu,wv,wh,wu,wv,wh)

	    tmass = 0
	    do j=1,n
		do i=2,np
		    tmass=tmass+wh(i,j)*c1(i,j)
		end do
	    end do
	    ! 
	    if (nyn.eq.0) then
		print *,'       (the integral time is      ',tt,')'
		nt=nt+1
	    endif
	    !
	    print *,tener,tmass,iter
	    nt=nt+1
	end if
    end do
    !------------------------- main part of the program end------------
    print *,'the main part of this program has ended'
    !-------------------------- out put the u,v,h ---------------------
    print *,'now,the working is to output result'

    open(1,file=ffu)
    open(2,file=ffv)
    open(3,file=ffh)
    !
    130	format(200f16.8)

    write(1,130) (( u(i,j),j=1,n),i=1,n1)
    write(2,130) (( v(i,j),j=1,n),i=1,n1)
    write(3,130) ((wh(i,j),j=1,n),i=1,n1)

    close(1)
    close(2)
    close(3)
    !
    do j=1,n
	do i=2,np
	    pu(i,j)=u(i,j)
	    pv(i,j)=v(i,j)
	    ph(i,j)=wh(i,j)
	end do
	pu(1,j)=pu(np,j)
	pv(1,j)=pv(np,j)
	ph(1,j)=ph(np,j)
    end do
    !
    call wr(pu,np,n,12,irecu)
    call wr(pv,np,n,12,irecu)
    call wr(ph,np,n,14,irech)
    !
    close(12)
    close(14)
    !
    stop
end
     
subroutine opf(nx,ny,ffn,nffn)
    character*6 :: ffn
    open(unit=nffn,file=ffn,form='unformatted',access='direct',recl=nx*ny*4)
    return
end subroutine opf
     
subroutine wr(rdata,nx,ny,nffn,irec)
    dimension rdata(nx+1,ny)
    irec = irec + 1
    write(unit=nffn,rec=irec)((rdata(i,j),i=1,nx),j=1,ny)
    return
end subroutine wr
     
function inner(u1,v1,h1,u2,v2,h2)
    use module_para
    implicit none

    real*8, dimension(1:n1,1:n) :: u1,v1,h1
    real*8, dimension(1:n1,1:n) :: u2,v2,h2
    real*8                      :: inner
    integer                     :: i,j
    
    inner = 0.0d0
    do j=1,n
	do i=2,np
	    inner=inner+(u1(i,j)*u2(i,j)+v1(i,j)*v2(i,j)+h1(i,j)*h2(i,j))*c1(i,j)
	end do
    end do
    return
end function inner
