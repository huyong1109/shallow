subroutine cs
    use module_para
    use module_array

    implicit none

    real*8  :: ai, aj	    ! working variables
    integer :: j	    ! working variable

    pi=datan(1d0)*4
    deta1=4*pi/p*a
    deta=pi/q
    deta2=deta*a*2
    detb=deta*nq

    do j=2,n-1
	ai    = j*deta-detb
	c1(j) = dcos(ai)
	s1(j) = dsin(ai)
    end do

    ai=1.5*deta-detb
    c1(1)=dcos(ai)/4
        
    ai=(n-0.5)*deta-detb
    c1(n)=dcos(ai)/4.
      
    s1(1)=-1.0d0
    s1(n)=1.0d0

    do j=1,n
	c2(j)=c1(j)
    end do

    do j=1,n
	ai=s1(j)
	aj=c1(j)*a
	f1(j)=2.0d0*omg0*ai
	f2(j)=ai/aj
    end do

    do j=1,n
	ai=deta1*c1(j)
	aj=deta2*c1(j)
	c11(j)=1/ai
	c12(j)=1/aj
	c13(j)=c11(j)*0.5
	c14(j)=c12(j)*0.5
    end do
     
    return
end subroutine cs
