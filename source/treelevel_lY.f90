

	MODULE	TreeLevel_lY

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


	SUBROUTINE lY0(l,g_i,f1,f2,lY_0)
	!Subroutine to cocmpute the gX correlation functions at tree level
	!in perturbation theory.
	!
	!	- X= A,P,V,S
	!
	!	g(x_0)^(0)f1f2 = 1/2 Tr{Q5f2f1 S(1,x_0)f1f1 gX S(x_0,1)f2f2}
	!
	!
	implicit none
	integer,intent(in) :: l, g_i, f1, f2		!size, gamma, flavour1, flavour2
	complex(kind=8),dimension(c0:c1),intent(out) :: lY_0

	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop
	complex(kind=8),dimension(4,4) :: mat1, mat2, mat3
	complex(kind=8),dimension(c0:c1,3) :: l0
	integer,dimension(3) :: nvec
	integer :: x0,i,k
	integer :: f

	nvec(:) = 0

	lY_0(:) = dcmplx(0.0d0,0.0d0)	
	l0(:,:) = dcmplx(0.0d0,0.0d0)	
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	f = f1f2(f2,f1)
	
	call Fermion_Propagator(l,nvec,S_Prop)

	do k=1,3
	   do x0=c0,c1

	   mat1 = matmul(Qk(:,:,k,f),S_prop(:,:,1,x0,f1))
	   mat2 = matmul(mat1,gY(:,:,g_i,k))
	   mat3 = matmul(mat2,S_prop(:,:,x0,1,f2))	

	   l0(x0,k) = mat3(1,1)+mat3(2,2)+mat3(3,3)+mat3(4,4)

	   enddo

	enddo

 	lY_0(:) = (1.0d0/6.0d0)*dimrep*(l0(:,1)+l0(:,2)+l0(:,3))

	END SUBROUTINE lY0


	END MODULE TreeLevel_lY

