
	MODULE Compute

	implicit none

	CONTAINS



	SUBROUTINE Compute_stuff(l)

	implicit none

	integer,intent(in) :: l
!	integer :: g_i, f1, f2
	integer,dimension(4) :: flav

	flav(1) = 1
	flav(2) = 1
	flav(3) = 1
	flav(4) = 2

	call compute_G1_4f(l,1,3,flav)
	call compute_G2_4f(l,1,3,flav)


	flav(1) = 1
	flav(2) = 2
	flav(3) = 1
	flav(4) = 1

	call compute_G1_4f(l,1,3,flav)
	call compute_G2_4f(l,1,3,flav)


	flav(1) = 1
	flav(2) = 1
	flav(3) = 1
	flav(4) = 1

	call compute_G1_4f(l,1,1,flav)
	call compute_G2_4f(l,1,1,flav)




	END SUBROUTINE Compute_stuff



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE compute_G1_4f(l,g_i1,g_i2,flav)

	USE Parameters
	USE TreeLevel_4fermions
	USE OneLoop_4fermions

	implicit none

	integer,intent(in) :: l, g_i1, g_i2		!size, gamma1, gamma2 
	integer,dimension(4),intent(in) :: flav

	complex(kind=8),dimension(c0:c1) :: G14f0
	complex(kind=8),dimension(c0:c1) :: G14f1a,G14f1b,G14f1c,G14f1d
	complex(kind=8),dimension(c0:c1) :: G14f2a,G14f2b
	complex(kind=8),dimension(c0:c1) :: G14f3a,G14f3b,G14f3c,G14f3d
	complex(kind=8),dimension(c0:c1) :: G14f4a,G14f4b,G14f4c,G14f4d
	complex(kind=8),dimension(c0:c1) :: G14f5a,G14f5b,G14f5c,G14f5d
	complex(kind=8),dimension(c0:c1) :: G14f6a,G14f6b,G14f6c,G14f6d
	complex(kind=8),dimension(c0:c1) :: G14f7a,G14f7b
	complex(kind=8),dimension(c0:c1) :: G14f1loop

	G14f0(:) = dcmplx(0.0d0,0.0d0)
	G14f1a(:) = dcmplx(0.0d0,0.0d0)
	G14f1b(:) = dcmplx(0.0d0,0.0d0)
	G14f1c(:) = dcmplx(0.0d0,0.0d0)
	G14f1d(:) = dcmplx(0.0d0,0.0d0)	
	G14f2a(:) = dcmplx(0.0d0,0.0d0)
	G14f2b(:) = dcmplx(0.0d0,0.0d0)
	G14f3a(:) = dcmplx(0.0d0,0.0d0)
	G14f3b(:) = dcmplx(0.0d0,0.0d0)
	G14f3c(:) = dcmplx(0.0d0,0.0d0)
	G14f3d(:) = dcmplx(0.0d0,0.0d0)
	G14f4a(:) = dcmplx(0.0d0,0.0d0)
	G14f4b(:) = dcmplx(0.0d0,0.0d0)
	G14f4c(:) = dcmplx(0.0d0,0.0d0)
	G14f4d(:) = dcmplx(0.0d0,0.0d0)
	G14f5a(:) = dcmplx(0.0d0,0.0d0)
	G14f5b(:) = dcmplx(0.0d0,0.0d0)
	G14f5c(:) = dcmplx(0.0d0,0.0d0)
	G14f5d(:) = dcmplx(0.0d0,0.0d0)
	G14f6a(:) = dcmplx(0.0d0,0.0d0)
	G14f6b(:) = dcmplx(0.0d0,0.0d0)
	G14f6c(:) = dcmplx(0.0d0,0.0d0)
	G14f6d(:) = dcmplx(0.0d0,0.0d0)
	G14f7a(:) = dcmplx(0.0d0,0.0d0)
	G14f7b(:) = dcmplx(0.0d0,0.0d0)
	G14f1loop(:) = dcmplx(0.0d0,0.0d0)
	
	call G1_4f_0(l,g_i1,g_i2,flav,G14f0)

	print*, ""
	print*, "G1_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_0=", G14f0(c0)
	print*, ""

	call G1_4f_1a(l,g_i1,g_i2,flav,G14f1a)
	print*, "G1_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_1a=", G14f1a(c0)

	G14f1b = G14f1a 
	print*, "G1_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_1b=", G14f1b(c0)

	call G1_4f_1c(l,g_i1,g_i2,flav,G14f1c)
	print*, "G1_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_1c=", G14f1c(c0)

	G14f1d = G14f1c 
	print*, "G1_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_1d=", G14f1d(c0)

	call G1_4f_2a(l,g_i1,g_i2,flav,G14f2a)
	print*, "G1_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_2a=", G14f2a(c0)
	call G1_4f_2b(l,g_i1,g_i2,flav,G14f2b)
	print*, "G1_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_2b=", G14f2b(c0)
	call G1_4f_3a(l,g_i1,g_i2,flav,G14f3a)
	print*, "G1_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_3a=", G14f3a(c0)
	call G1_4f_3b(l,g_i1,g_i2,flav,G14f3b)
	print*, "G1_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_3b=", G14f3b(c0)
	call G1_4f_3c(l,g_i1,g_i2,flav,G14f3c)
	print*, "G1_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_3c=", G14f3c(c0)
	call G1_4f_3d(l,g_i1,g_i2,flav,G14f3d)
	print*, "G1_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_3d=", G14f3d(c0)
	call G1_4f_4a(l,g_i1,g_i2,flav,G14f4a)
	print*, "G1_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_4a=", G14f4a(c0)
	call G1_4f_4b(l,g_i1,g_i2,flav,G14f4b)
	print*, "G1_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_4b=", G14f4b(c0)
	call G1_4f_4c(l,g_i1,g_i2,flav,G14f4c)
	print*, "G1_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_4c=", G14f4c(c0)
	call G1_4f_4d(l,g_i1,g_i2,flav,G14f4d)
	print*, "G1_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_4d=", G14f4d(c0)
	call G1_4f_5a(l,g_i1,g_i2,flav,G14f5a)
	print*, "G1_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_5a=", G14f5a(c0)
	call G1_4f_5b(l,g_i1,g_i2,flav,G14f5b)
	print*, "G1_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_5b=", G14f5b(c0)
	call G1_4f_5c(l,g_i1,g_i2,flav,G14f5c)
	print*, "G1_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_5c=", G14f5c(c0)
!	call G1_4f_5d(l,g_i1,g_i2,flav,G14f5d)
	G14f5d(:) = G14f5c(:) 
	print*, "G1_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_5d=", G14f5d(c0)
	call G1_4f_6a(l,g_i1,g_i2,flav,G14f6a)
	print*, "G1_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_6a=", G14f6a(c0)
	call G1_4f_6b(l,g_i1,g_i2,flav,G14f6b)
	print*, "G1_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_6b=", G14f6b(c0)
	call G1_4f_6c(l,g_i1,g_i2,flav,G14f6c)
	print*, "G1_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_6c=", G14f6c(c0)
	call G1_4f_6d(l,g_i1,g_i2,flav,G14f6d)
	print*, "G1_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_6d=", G14f6d(c0)
	call G1_4f_7a(l,g_i1,g_i2,flav,G14f7a)
	print*, "G1_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_7a=", G14f7a(c0)
	call G1_4f_7b(l,g_i1,g_i2,flav,G14f7b)
	print*, "G1_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_7b=", G14f7b(c0)

	G14f1loop(:) = G14f1a(:) + G14f1b(:) + G14f1c(:) + G14f1d(:)
	G14f1loop(:) = G14f1loop(:) + G14f2a(:) + G14f2b(:)
	G14f1loop(:) = G14f1loop(:) + G14f3a(:) + G14f3b(:) + G14f3c(:) + G14f3d(:)
	G14f1loop(:) = G14f1loop(:) + G14f4a(:) + G14f4b(:) + G14f4c(:) + G14f4d(:)
	G14f1loop(:) = G14f1loop(:) + G14f5a(:) + G14f5b(:) + G14f5c(:) + G14f5d(:)
	G14f1loop(:) = G14f1loop(:) + G14f6a(:) + G14f6b(:) + G14f6c(:) + G14f6d(:)
	G14f1loop(:) = G14f1loop(:) + G14f7a(:) + G14f7b(:)

	print*, ""
	print*, "G1_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_1loop=", G14f1loop(c0)
	print*, ""

	END SUBROUTINE compute_G1_4f

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE compute_G2_4f(l,g_i1,g_i2,flav)

	USE Parameters
	USE TreeLevel_4fermions
	USE OneLoop_4fermions

	implicit none

	integer,intent(in) :: l, g_i1, g_i2		!size, gamma1, gamma2 
	integer,dimension(4),intent(in) :: flav

	complex(kind=8),dimension(c0:c1) :: G24f0
	complex(kind=8),dimension(c0:c1) :: G24f1a,G24f1b,G24f1c,G24f1d
	complex(kind=8),dimension(c0:c1) :: G24f2a,G24f2b
	complex(kind=8),dimension(c0:c1) :: G24f3a,G24f3b,G24f3c,G24f3d
	complex(kind=8),dimension(c0:c1) :: G24f4a,G24f4b,G24f4c,G24f4d
	complex(kind=8),dimension(c0:c1) :: G24f5a,G24f5b,G24f5c,G24f5d
	complex(kind=8),dimension(c0:c1) :: G24f6a,G24f6b,G24f6c,G24f6d
	complex(kind=8),dimension(c0:c1) :: G24f7a,G24f7b
	complex(kind=8),dimension(c0:c1) :: G24f8a,G24f8b
	complex(kind=8),dimension(c0:c1) :: G24f9a,G24f9b
	complex(kind=8),dimension(c0:c1) :: G24f10a,G24f10b,G24f10c,G24f10d
	complex(kind=8),dimension(c0:c1) :: G24f11a,G24f11b,G24f11c,G24f11d
	complex(kind=8),dimension(c0:c1) :: G24f12a,G24f12b
	complex(kind=8),dimension(c0:c1) :: G24f13a,G24f13b
	complex(kind=8),dimension(c0:c1) :: G24f1loop

	G24f0(:) = dcmplx(0.0d0,0.0d0)
	G24f1a(:) = dcmplx(0.0d0,0.0d0)
	G24f1b(:) = dcmplx(0.0d0,0.0d0)
	G24f1c(:) = dcmplx(0.0d0,0.0d0)
	G24f1d(:) = dcmplx(0.0d0,0.0d0)	
	G24f2a(:) = dcmplx(0.0d0,0.0d0)
	G24f2b(:) = dcmplx(0.0d0,0.0d0)
	G24f3a(:) = dcmplx(0.0d0,0.0d0)
	G24f3b(:) = dcmplx(0.0d0,0.0d0)
	G24f3c(:) = dcmplx(0.0d0,0.0d0)
	G24f3d(:) = dcmplx(0.0d0,0.0d0)
	G24f4a(:) = dcmplx(0.0d0,0.0d0)
	G24f4b(:) = dcmplx(0.0d0,0.0d0)
	G24f4c(:) = dcmplx(0.0d0,0.0d0)
	G24f4d(:) = dcmplx(0.0d0,0.0d0)
	G24f5a(:) = dcmplx(0.0d0,0.0d0)
	G24f5b(:) = dcmplx(0.0d0,0.0d0)
	G24f5c(:) = dcmplx(0.0d0,0.0d0)
	G24f5d(:) = dcmplx(0.0d0,0.0d0)
	G24f6a(:) = dcmplx(0.0d0,0.0d0)
	G24f6b(:) = dcmplx(0.0d0,0.0d0)
	G24f6c(:) = dcmplx(0.0d0,0.0d0)
	G24f6d(:) = dcmplx(0.0d0,0.0d0)
	G24f7a(:) = dcmplx(0.0d0,0.0d0)
	G24f7b(:) = dcmplx(0.0d0,0.0d0)
	G24f8a(:) = dcmplx(0.0d0,0.0d0)
	G24f8b(:) = dcmplx(0.0d0,0.0d0)
	G24f9a(:) = dcmplx(0.0d0,0.0d0)
	G24f9b(:) = dcmplx(0.0d0,0.0d0)
	G24f10a(:) = dcmplx(0.0d0,0.0d0)
	G24f10b(:) = dcmplx(0.0d0,0.0d0)
	G24f10c(:) = dcmplx(0.0d0,0.0d0)
	G24f10d(:) = dcmplx(0.0d0,0.0d0)
	G24f11a(:) = dcmplx(0.0d0,0.0d0)
	G24f11b(:) = dcmplx(0.0d0,0.0d0)
	G24f11c(:) = dcmplx(0.0d0,0.0d0)
	G24f11d(:) = dcmplx(0.0d0,0.0d0)
	G24f12a(:) = dcmplx(0.0d0,0.0d0)
	G24f12b(:) = dcmplx(0.0d0,0.0d0)
	G24f13a(:) = dcmplx(0.0d0,0.0d0)
	G24f13b(:) = dcmplx(0.0d0,0.0d0)
	G24f1loop(:) = dcmplx(0.0d0,0.0d0)
	
	call G2_4f_0(l,g_i1,g_i2,flav,G24f0)

	print*, ""
	print*, "G2_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_0=", G24f0(c0)
	print*, ""

	call G2_4f_1a(l,g_i1,g_i2,flav,G24f1a)
	print*, "G2_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_1a=", G24f1a(c0)

	G24f1b = G24f1a 
	print*, "G2_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_1b=", G24f1b(c0)

	call G2_4f_1c(l,g_i1,g_i2,flav,G24f1c)
	print*, "G2_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_1c=", G24f1c(c0)

	G24f1d = G24f1c 
	print*, "G2_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_1d=", G24f1d(c0)

	call G2_4f_2a(l,g_i1,g_i2,flav,G24f2a)
	print*, "G2_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_2a=", G24f2a(c0)
	call G2_4f_2b(l,g_i1,g_i2,flav,G24f2b)
	print*, "G2_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_2b=", G24f2b(c0)
	call G2_4f_3a(l,g_i1,g_i2,flav,G24f3a)
	print*, "G2_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_3a=", G24f3a(c0)
	call G2_4f_3b(l,g_i1,g_i2,flav,G24f3b)
	print*, "G2_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_3b=", G24f3b(c0)
	call G2_4f_3c(l,g_i1,g_i2,flav,G24f3c)
	print*, "G2_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_3c=", G24f3c(c0)
	call G2_4f_3d(l,g_i1,g_i2,flav,G24f3d)
	print*, "G2_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_3d=", G24f3d(c0)
	call G2_4f_4a(l,g_i1,g_i2,flav,G24f4a)
	print*, "G2_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_4a=", G24f4a(c0)
	call G2_4f_4b(l,g_i1,g_i2,flav,G24f4b)
	print*, "G2_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_4b=", G24f4b(c0)
	call G2_4f_4c(l,g_i1,g_i2,flav,G24f4c)
	print*, "G2_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_4c=", G24f4c(c0)
	call G2_4f_4d(l,g_i1,g_i2,flav,G24f4d)
	print*, "G2_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_4d=", G24f4d(c0)
	call G2_4f_5a(l,g_i1,g_i2,flav,G24f5a)
	print*, "G2_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_5a=", G24f5a(c0)
	call G2_4f_5b(l,g_i1,g_i2,flav,G24f5b)
	print*, "G2_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_5b=", G24f5b(c0)
	call G2_4f_5c(l,g_i1,g_i2,flav,G24f5c)
	print*, "G2_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_5c=", G24f5c(c0)
!	call G2_4f_5d(l,g_i1,g_i2,flav,G24f5d)
	G24f5d(:)=G24f5c(:)
	print*, "G2_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_5d=", G24f5d(c0)
	call G2_4f_6a(l,g_i1,g_i2,flav,G24f6a)
	print*, "G2_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_6a=", G24f6a(c0)
	call G2_4f_6b(l,g_i1,g_i2,flav,G24f6b)
	print*, "G2_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_6b=", G24f6b(c0)
	call G2_4f_6c(l,g_i1,g_i2,flav,G24f6c)
	print*, "G2_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_6c=", G24f6c(c0)
!	call G2_4f_6d(l,g_i1,g_i2,flav,G24f6d)
	G24f6d(:)=G24f6c(:)
	print*, "G2_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_6d=", G24f6d(c0)
	call G2_4f_7a(l,g_i1,g_i2,flav,G24f7a)
	print*, "G2_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_7a=", G24f7a(c0)
	call G2_4f_7b(l,g_i1,g_i2,flav,G24f7b)
	print*, "G2_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_7b=", G24f7b(c0)
	call G2_4f_8a(l,g_i1,g_i2,flav,G24f8a)
	print*, "G2_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_8a=", G24f8a(c0)
	call G2_4f_8b(l,g_i1,g_i2,flav,G24f8b)
	print*, "G2_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_8b=", G24f8b(c0)
	call G2_4f_9a(l,g_i1,g_i2,flav,G24f9a)
	print*, "G2_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_9a=", G24f9a(c0)
	call G2_4f_9b(l,g_i1,g_i2,flav,G24f9b)
	print*, "G2_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_9b=", G24f9b(c0)
	call G2_4f_10a(l,g_i1,g_i2,flav,G24f10a)
	print*, "G2_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_10a=", G24f10a(c0)
	call G2_4f_10b(l,g_i1,g_i2,flav,G24f10b)
	print*, "G2_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_10b=", G24f10b(c0)
	call G2_4f_10c(l,g_i1,g_i2,flav,G24f10c)
	print*, "G2_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_10c=", G24f10c(c0)
	call G2_4f_10d(l,g_i1,g_i2,flav,G24f10d)
	print*, "G2_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_10d=", G24f10d(c0)
	call G2_4f_11a(l,g_i1,g_i2,flav,G24f11a)
	print*, "G2_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_11a=", G24f11a(c0)
	call G2_4f_11b(l,g_i1,g_i2,flav,G24f11b)
	print*, "G2_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_11b=", G24f11b(c0)
	call G2_4f_11c(l,g_i1,g_i2,flav,G24f11c)
	print*, "G2_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_11c=", G24f11c(c0)
	call G2_4f_11d(l,g_i1,g_i2,flav,G24f11d)
	print*, "G2_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_11d=", G24f11d(c0)
	call G2_4f_12a(l,g_i1,g_i2,flav,G24f12a)
	print*, "G2_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_12a=", G24f12a(c0)
	call G2_4f_12b(l,g_i1,g_i2,flav,G24f12b)
	print*, "G2_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_12b=", G24f12b(c0)
	call G2_4f_13a(l,g_i1,g_i2,flav,G24f13a)
	print*, "G2_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_13a=", G24f13a(c0)
	call G2_4f_13b(l,g_i1,g_i2,flav,G24f13b)
	print*, "G2_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_13b=", G24f13b(c0)


	G24f1loop(:) = G24f1a(:) + G24f1b(:) + G24f1c(:) + G24f1d(:)
	G24f1loop(:) = G24f1loop(:) + G24f2a(:) + G24f2b(:)
	G24f1loop(:) = G24f1loop(:) + G24f3a(:) + G24f3b(:) + G24f3c(:) + G24f3d(:)
	G24f1loop(:) = G24f1loop(:) + G24f4a(:) + G24f4b(:) + G24f4c(:) + G24f4d(:)
	G24f1loop(:) = G24f1loop(:) + G24f5a(:) + G24f5b(:) + G24f5c(:) + G24f5d(:)
	G24f1loop(:) = G24f1loop(:) + G24f6a(:) + G24f6b(:) + G24f6c(:) + G24f6d(:)
	G24f1loop(:) = G24f1loop(:) + G24f7a(:) + G24f7b(:)
	G24f1loop(:) = G24f1loop(:) + G24f8a(:) + G24f8b(:)
	G24f1loop(:) = G24f1loop(:) + G24f9a(:) + G24f9b(:)
	G24f1loop(:) = G24f1loop(:) + G24f10a(:) + G24f10b(:) + G24f10c(:) + G24f10d(:)
	G24f1loop(:) = G24f1loop(:) + G24f11a(:) + G24f11b(:) + G24f11c(:) + G24f11d(:)
	G24f1loop(:) = G24f1loop(:) + G24f12a(:) + G24f12b(:)
	G24f1loop(:) = G24f1loop(:) + G24f13a(:) + G24f13b(:)
	print*, ""
	print*, "G2_", corrname_gXff(g_i1),corrname_gXff(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_1loop=", G24f1loop(c0)
	print*, ""



	END SUBROUTINE compute_G2_4f

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	END MODULE Compute
