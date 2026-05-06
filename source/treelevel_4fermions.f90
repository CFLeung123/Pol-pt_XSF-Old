

	MODULE	TreeLevel_4fermions
	
	USE GammaMatrices
	USE Parameters
	USE Propagators

	implicit none

	CONTAINS


	SUBROUTINE G1_4f_0(l,g_i1,g_i2,flav,G14f0)

	implicit none
	integer,intent(in) :: l, g_i1, g_i2		!size, gamma1, gamma2 
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G14f0

	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop
	complex(kind=8) :: part1, part2
	integer,dimension(3) :: nvec
	integer :: x0,i
	integer :: f

	nvec(:) = 0

	G14f0(:) = dcmplx(0.0d0,0.0d0)	
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

	
	call Fermion_Propagator(l,nvec,S_Prop)

	do x0=c0,c1

	   f = f1f2(flav(2),flav(1))
	   part1 = Trace4(Q5p(:,:,f),S_prop(:,:,l-1,x0,flav(1)),g4f(:,:,g_i1),S_prop(:,:,x0,l-1,flav(2)))

	   f = f1f2(flav(4),flav(3))
	   part2 = Trace4(Q5(:,:,f),S_prop(:,:,1,x0,flav(3)),g4f(:,:,g_i2),S_prop(:,:,x0,1,flav(4)))

	   G14f0(x0) = part1*part2*dimrep**2

	enddo

	END SUBROUTINE G1_4f_0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G2_4f_0(l,g_i1,g_i2,flav,G24f0)

	implicit none
	integer,intent(in) :: l, g_i1, g_i2		!size, gamma1, gamma2 
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G24f0

	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop
	complex(kind=8),dimension(4,4) :: mat1, mat2, mat3 ,mat4
	integer,dimension(3) :: nvec
	integer :: x0,i
	integer :: f1,f2

	nvec(:) = 0

	G24f0(:) = dcmplx(0.0d0,0.0d0)	
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)


	
	call Fermion_Propagator(l,nvec,S_Prop)

	f1 = f1f2(flav(2),flav(1))
	f2 = f1f2(flav(4),flav(3))

	do x0=c0,c1

	   mat1 = matmul(Q5p(:,:,f1),S_prop(:,:,l-1,x0,flav(1)))
	   mat2 = matmul(g4f(:,:,g_i1),S_prop(:,:,x0,1,flav(4)))
	   mat3 = matmul(Q5(:,:,f2),S_prop(:,:,1,x0,flav(3)))
	   mat4 = matmul(g4f(:,:,g_i2),S_prop(:,:,x0,l-1,flav(2)))

	   G24f0(x0) = -dimrep*Trace4(mat1,mat2,mat3,mat4)

	enddo

	END SUBROUTINE G2_4f_0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE Gk1_4f_0(l,g_i1,g_i2,flav,Gk14f0)

	implicit none
	integer,intent(in) :: l, g_i1, g_i2		!size, gamma1, gamma2 
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: Gk14f0

	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop_1, S_prop_2
	complex(kind=8) :: part1, part2
	integer,dimension(3) :: nvec_1, nvec_2
	integer :: x0,i
	integer :: f,k

	nvec_1(:) = 0
	nvec_2(:) = 0

	Gk14f0(:) = dcmplx(0.0d0,0.0d0)	
	S_prop_1(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop_2(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

	
	call Fermion_Propagator(l,nvec_1,S_Prop_1)
	call Fermion_Propagator(l,nvec_2,S_Prop_2)

!	do k=1,3
	do x0=c0,c1

	   f = f1f2(flav(2),flav(1))
	   part1 = Trace4(Qkp(:,:,kdir,f),S_prop_1(:,:,l-1,x0,flav(1)),g4f(:,:,g_i1),S_prop_1(:,:,x0,l-1,flav(2)))

	   f = f1f2(flav(4),flav(3))
	   part2 = Trace4(Qk(:,:,kdir,f),S_prop_2(:,:,1,x0,flav(3)),g4f(:,:,g_i2),S_prop_2(:,:,x0,1,flav(4)))

	   Gk14f0(x0) = Gk14f0(x0) + part1*part2*dimrep**2

	enddo
!	enddo

        Gk14f0(:) = Gk14f0(:)/3.0d0

	END SUBROUTINE Gk1_4f_0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE Gk2_4f_0(l,g_i1,g_i2,flav,Gk24f0)

	implicit none
	integer,intent(in) :: l, g_i1, g_i2		!size, gamma1, gamma2 
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: Gk24f0

	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop_1, S_prop_2
	complex(kind=8),dimension(4,4) :: mat1, mat2, mat3 ,mat4
	integer,dimension(3) :: nvec_1, nvec_2
	integer :: x0,i,k
	integer :: f1,f2

	nvec_1(:) = 0
	nvec_2(:) = 0

	Gk24f0(:) = dcmplx(0.0d0,0.0d0)	
	S_prop_1(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop_2(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

	
	call Fermion_Propagator(l,nvec_1,S_Prop_1)
	call Fermion_Propagator(l,nvec_2,S_Prop_2)

	f1 = f1f2(flav(4),flav(3))
	f2 = f1f2(flav(2),flav(1))

!	do k=1,3
	do x0=c0,c1

	   mat1 = matmul(Qk(:,:,kdir,f1),S_prop_1(:,:,1,x0,flav(3)))
	   mat2 = matmul(g4f(:,:,g_i2),S_prop_2(:,:,x0,l-1,flav(2)))
	   mat3 = matmul(Qkp(:,:,kdir,f2),S_prop_2(:,:,l-1,x0,flav(1)))
	   mat4 = matmul(g4f(:,:,g_i1),S_prop_1(:,:,x0,1,flav(4)))

	   Gk24f0(x0) = Gk24f0(x0) -dimrep*Trace4(mat1,mat2,mat3,mat4)

	enddo
!	enddo

	Gk24f0(:) = Gk24f0(:)/3.0d0 

	END SUBROUTINE Gk2_4f_0


	END MODULE TreeLevel_4fermions
