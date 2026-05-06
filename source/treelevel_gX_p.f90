

	MODULE	TreeLevel_gX_p
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


	SUBROUTINE gX0_p(l,g_i,f1,f2,gX_0_p)
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
	complex(kind=8),dimension(c0:c1),intent(out) :: gX_0_p

	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop
	complex(kind=8),dimension(4,4) :: mat1, mat2, mat3
	integer,dimension(3) :: nvec
	integer :: x0,i
	integer :: f

	nvec(:) = 0

	gX_0_p(:) = dcmplx(0.0d0,0.0d0)	
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	f = f1f2(f2,f1)
	
	call Fermion_Propagator(l,nvec,S_Prop)

	   do x0=c0,c1

	   mat1 = matmul(Q5p(:,:,f),S_prop(:,:,l-1,x0,f1))
	   mat2 = matmul(mat1,gX(:,:,g_i))
	   mat3 = matmul(mat2,S_prop(:,:,x0,l-1,f2))	

	   gX_0_p(x0) = mat3(1,1)+mat3(2,2)+mat3(3,3)+mat3(4,4)

	   enddo

 	gX_0_p(:) = 0.5d0*dimrep*gX_0_p(:)

	END SUBROUTINE gX0_p


	END MODULE TreeLevel_gX_p

