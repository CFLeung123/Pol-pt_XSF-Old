

	MODULE	TreeLevel_lVt
	!Module with the subroutines necessary for the calculation of correlation
	!functions at tree level in perturbation theory.
	!
	! Contains the subroutines:
	!    SF:
	!	-fX0(l,S_prop,fX_0)
	!
	USE GammaMatrices
	USE Parameters
	USE Propagators

	implicit none

	CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	XSF: Tree level diagrams
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE lVt0(l,f1,f2,lVt_0)
	!Subroutine to cocmpute the gX correlation functions at tree level
	!in perturbation theory.
	!
	!	- X= A,P,V,S
	!
	!	g(x_0)^(0)f1f2 = 1/2 Tr{Q5f2f1 S(1,x_0)f1f1 gX S(x_0,1)f2f2}
	!
	!
	implicit none
	integer,intent(in) :: l, f1, f2		!size, gamma, flavour1, flavour2
	complex(kind=8),dimension(c0:c1),intent(out) :: lVt_0

	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop
	complex(kind=8),dimension(c0:c1,3) :: trace_1, trace_2
	integer,dimension(3) :: nvec
	integer :: x0,i
	integer :: f,k
	real(kind=8),dimension(3) :: p

	nvec(:) = 0

	lVt_0(:) = dcmplx(0.0d0,0.0d0)	
	trace_1(:,:) = dcmplx(0.0d0,0.0d0)	
	trace_2(:,:) = dcmplx(0.0d0,0.0d0)	
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	f = f1f2(f2,f1)
	


	!Prepare momenta for the exponential fator:
	 p(1)=(2.0d0*pi*nvec(1))/(l*1.0d0)	!momenta, plus theta			
	 p(2)=(2.0d0*pi*nvec(2))/(l*1.0d0)
	 p(3)=(2.0d0*pi*nvec(3))/(l*1.0d0)


	call Fermion_Propagator(l,nvec,S_Prop)

	do x0=c0,c1
	do k=1,3
	   trace_1(x0,k) = Trace4(Qk(:,:,k,f),S_prop(:,:,1,x0,f1),Pkplus(:,:,k),S_prop(:,:,x0,1,f2))!*exp(-imag*p(k))
	   trace_2(x0,k) = Trace4(Qk(:,:,k,f),S_prop(:,:,1,x0,f1),Pkminus(:,:,k),S_prop(:,:,x0,1,f2))!*exp(imag*p(k))
	enddo
	enddo


	lVt_0(:) = trace_1(:,1)+trace_1(:,2)+trace_1(:,3) - trace_2(:,1)-trace_2(:,2)-trace_2(:,3)

	lVt_0(:) = (1.0d0/6.0d0)*dimrep*lVt_0(:)

	END SUBROUTINE lVt0


	END MODULE TreeLevel_lVt

