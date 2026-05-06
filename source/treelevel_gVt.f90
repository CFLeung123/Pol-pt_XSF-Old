

	MODULE	TreeLevel_gVt
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


	SUBROUTINE gVt0(l,f1,f2,gVt_0)
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
	complex(kind=8),dimension(c0:c1),intent(out) :: gVt_0

	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop
	complex(kind=8),dimension(c0:c1) :: trace_1, trace_2
	integer,dimension(3) :: nvec
	integer :: x0,i
	integer :: f

	nvec(:) = 0

	gVt_0(:) = dcmplx(0.0d0,0.0d0)	
	trace_1(:) = dcmplx(0.0d0,0.0d0)	
	trace_2(:) = dcmplx(0.0d0,0.0d0)	
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	f = f1f2(f2,f1)
	
	call Fermion_Propagator(l,nvec,S_Prop)

	do x0=c0,c1

	   trace_1(x0) = Trace4(Q5(:,:,f),S_prop(:,:,1,x0+1,f1),Pplus,S_prop(:,:,x0,1,f2))
	   trace_2(x0) = Trace4(Q5(:,:,f),S_prop(:,:,1,x0,f1),Pminus,S_prop(:,:,x0+1,1,f2))

	enddo

 	gVt_0(:) = 0.5d0*dimrep*(trace_1(:) - trace_2(:))

	END SUBROUTINE gVt0


	END MODULE TreeLevel_gVt

