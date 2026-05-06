
	MODULE Compute
	!
	!  Module containing:
	!
	! 	- A control subroutine "Compute_stuff" that gets information 
	!	  about what to compute and calls the subroutines which compute
	!	  the requested observables.
	!	
	!	- A set of subroutines "compute_X", where X is a specific observable.
	!	  These are control subroutines that call the functions to calculate
	!	  every diagram.	 
	!	  Then the gauge invariant sums over all diagrams for a specific
	!	  observable are performed.
	!	  Finally, results are printed on screen
	!
	implicit none
	CONTAINS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE Compute_stuff(l,g_i,g_i2,g4f_i,g4f_i2,f1,f2,f3,f4)
	!
	!  Get information about what to compute, and call the specific 
	!  routines to compute observables.
	!
	USE Input
	implicit none

	integer,intent(in) :: l
	integer,intent(in) :: g_i,g_i2,g4f_i,g4f_i2, f1, f2, f3 , f4
	integer,dimension(4) :: flav, flavp

	flav(1) = f1
	flav(2) = f2
	flav(3) = f3
	flav(4) = f4

	if(comp_gX.eq.1)then
	  call compute_gX(l,g_i,f1,f2)
	endif

	if(comp_lY.eq.1)then
	  call compute_lY(l,g_i,f1,f2)
	endif

	if(comp_g1.eq.1)then
	  call compute_g1(l,f1,f2)
	endif

	if(comp_l1.eq.1)then
	  call compute_l1(l,f1,f2)
	endif

	if(comp_gVt.eq.1)then
	  call compute_gVt(l,f1,f2)
	endif

	if(comp_lVt.eq.1)then
	  call compute_lVt(l,f1,f2)
	endif

	if(comp_G1_4f.eq.1)then
	  call compute_G1_4f(l,g4f_i,g4f_i2,flav)
	endif

	if(comp_G2_4f.eq.1)then
	  call compute_G2_4f(l,g4f_i,g4f_i2,flav)
	endif

	if(comp_Gk1_4f.eq.1)then
	  call compute_Gk1_4f(l,g4f_i,g4f_i2,flav)
	endif

	if(comp_Gk2_4f.eq.1)then
	  call compute_Gk2_4f(l,g4f_i,g4f_i2,flav)
	endif


	END SUBROUTINE Compute_stuff



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE compute_gX(l,g_i,f1,f2)
	!
	! Control subroutine to compute all diagrams contributing to gX
	!
	USE Parameters
	USE Input
	USE TreeLevel_gX
	USE OneLoop_gX


	implicit none

	integer,intent(in) :: l,g_i,f1,f2
	complex(kind=8),dimension(c0:c1) :: gX_0
	complex(kind=8),dimension(c0:c1) :: gX_1a, gX_1b
	complex(kind=8),dimension(c0:c1) :: gX_2
	complex(kind=8),dimension(c0:c1) :: gX_3a, gX_3b
	complex(kind=8),dimension(c0:c1) :: gX_4a, gX_4b
	complex(kind=8),dimension(c0:c1) :: gX_5a, gX_5b
	complex(kind=8),dimension(c0:c1) :: gX_6a, gX_6b
	complex(kind=8),dimension(c0:c1) :: gX_7
	complex(kind=8),dimension(c0:c1) :: gX_1loop
	complex(kind=8),dimension(c0:c1) :: gX_dm, gX_zf, gX_ds

	integer :: i

	gX_1loop(:) = dcmplx(0.0d0,0.0d0)

	if(comp_tree.eq.1)then	
	call gX0(l,g_i,f1,f2,gX_0)
	endif

	if(comp_1loop.eq.1)then	
	call gX1a(l,g_i,f1,f2,gX_1a)
	gX_1b(:) = gX_1a(:)
	call gX2(l,g_i,f1,f2,gX_2)
	call gX3a(l,g_i,f1,f2,gX_3a)
	call gX3b(l,g_i,f1,f2,gX_3b)
	call gX4a(l,g_i,f1,f2,gX_4a)
	call gX4b(l,g_i,f1,f2,gX_4b)
	call gX5a(l,g_i,f1,f2,gX_5a)
	call gX5b(l,g_i,f1,f2,gX_5b)
	call gX6a(l,g_i,f1,f2,gX_6a)
	call gX6b(l,g_i,f1,f2,gX_6b)
	call gX7(l,g_i,f1,f2,gX_7)

	gX_1loop(:) = 2.0d0*gX_1a(:)
	gX_1loop(:) = gX_1loop(:) + gX_2(:)
	gX_1loop(:) = gX_1loop(:) + gX_3a(:)
	gX_1loop(:) = gX_1loop(:) + gX_3b(:)
	gX_1loop(:) = gX_1loop(:) + gX_4a(:)
	gX_1loop(:) = gX_1loop(:) + gX_4b(:)
	gX_1loop(:) = gX_1loop(:) + gX_5a(:)
	gX_1loop(:) = gX_1loop(:) + gX_5b(:)
	gX_1loop(:) = gX_1loop(:) + gX_6a(:)
	gX_1loop(:) = gX_1loop(:) + gX_6b(:)
	gX_1loop(:) = gX_1loop(:) + gX_7(:)


	endif

	if(comp_count.eq.1)then	
	call gXdm(l,g_i,f1,f2,gX_dm)
	call gXzf(l,g_i,f1,f2,gX_zf)
	call gXds(l,g_i,f1,f2,gX_ds)
	endif
	

	do i=c0,c1

	print*, ""
	
	if(comp_tree.eq.1)then	
	print*, corrname_gXff(g_i),flavname(f1f2(f1,f2)),"_0:",i, real(gX_0(i)), aimag(gX_0(i))
	endif

	if(comp_1loop.eq.1)then	
	print*, corrname_gXff(g_i),flavname(f1f2(f1,f2)),"_1a:",i, real(gX_1a(i)), aimag(gX_1a(i))
	print*, corrname_gXff(g_i),flavname(f1f2(f1,f2)),"_1a:",i, real(gX_1b(i)), aimag(gX_1b(i))
	print*, corrname_gXff(g_i),flavname(f1f2(f1,f2)),"_2:",i, real(gX_2(i)), aimag(gX_2(i))
	print*, corrname_gXff(g_i),flavname(f1f2(f1,f2)),"_3a:",i, real(gX_3a(i)), aimag(gX_3a(i))
	print*, corrname_gXff(g_i),flavname(f1f2(f1,f2)),"_3b:",i, real(gX_3b(i)), aimag(gX_3b(i))
	print*, corrname_gXff(g_i),flavname(f1f2(f1,f2)),"_4a:",i, real(gX_4a(i)), aimag(gX_4a(i))
	print*, corrname_gXff(g_i),flavname(f1f2(f1,f2)),"_4b:",i, real(gX_4b(i)), aimag(gX_4b(i))
	print*, corrname_gXff(g_i),flavname(f1f2(f1,f2)),"_5a:",i, real(gX_5a(i)), aimag(gX_5a(i))
	print*, corrname_gXff(g_i),flavname(f1f2(f1,f2)),"_5b:",i, real(gX_5b(i)), aimag(gX_5b(i))
	print*, corrname_gXff(g_i),flavname(f1f2(f1,f2)),"_6a:",i, real(gX_6a(i)), aimag(gX_6a(i))
	print*, corrname_gXff(g_i),flavname(f1f2(f1,f2)),"_6b:",i, real(gX_6b(i)), aimag(gX_6b(i))
	print*, corrname_gXff(g_i),flavname(f1f2(f1,f2)),"_7:",i, real(gX_7(i)), aimag(gX_7(i))
	print*, corrname_gXff(g_i),flavname(f1f2(f1,f2)),"_1loop:",i, real(gX_1loop(i)), aimag(gX_1loop(i))
	endif

	if(comp_count.eq.1)then	
	print*, corrname_gXff(g_i),flavname(f1f2(f1,f2)),"_dm:",i, real(gX_dm(i)), aimag(gX_dm(i))
	print*, corrname_gXff(g_i),flavname(f1f2(f1,f2)),"_zf:",i, real(gX_zf(i)), aimag(gX_zf(i))
	print*, corrname_gXff(g_i),flavname(f1f2(f1,f2)),"_ds:",i, real(gX_ds(i)), aimag(gX_ds(i))
	endif

	print*, ""
	print*, ""
	enddo

	END SUBROUTINE compute_gX




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!









!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE compute_g1(l,f1,f2)
	!
	! Control subroutine to compute all diagrams contributing to g1
	!
	USE Parameters
	USE Input
	USE TreeLevel_g1
	USE OneLoop_g1

	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8) :: g1_0
	complex(kind=8) :: g1_1a, g1_1b
	complex(kind=8) :: g1_2, g1_7
	complex(kind=8) :: g1_4a, g1_4b
	complex(kind=8) :: g1_3a, g1_3b
	complex(kind=8) :: g1_5a, g1_5b
	complex(kind=8) :: g1_6a, g1_6b
	complex(kind=8) :: g1_8a, g1_8b
	complex(kind=8) :: g1_9
	complex(kind=8) :: g1_10a, g1_10b
	complex(kind=8) :: g1_11a, g1_11b
	complex(kind=8) :: g1_12a, g1_12b
	complex(kind=8) :: g1_13a, g1_13b
	complex(kind=8) :: g1_dm, g1_zf, g1_ds
	complex(kind=8) :: g1_1loop

	g1_1loop = dcmplx(0.0d0,0.0d0) 	

	if(comp_tree.eq.1)then	
	call g10(l,f1,f2,g1_0)	
	endif

	if(comp_1loop.eq.1)then	
	call g11a(l,f1,f2,g1_1a)
	g1_1b = g1_1a	
	call g12(l,f1,f2,g1_2)	
	call g13a(l,f1,f2,g1_3a)
	call g13b(l,f1,f2,g1_3b)
	call g14a(l,f1,f2,g1_4a)
	call g14b(l,f1,f2,g1_4b)
	call g15a(l,f1,f2,g1_5a)
	call g15b(l,f1,f2,g1_5b)
	call g16a(l,f1,f2,g1_6a)							
	call g16b(l,f1,f2,g1_6b)	
	call g17(l,f1,f2,g1_7)	
	call g18a(l,f1,f2,g1_8a)	
	g1_8b = g1_8a
	call g19(l,f1,f2,g1_9)	
	call g110a(l,f1,f2,g1_10a)	
	call g110b(l,f1,f2,g1_10b)	
	call g111a(l,f1,f2,g1_11a)	
	call g111b(l,f1,f2,g1_11b)	
	call g112a(l,f1,f2,g1_12a)
	call g112b(l,f1,f2,g1_12b)	
	call g113a(l,f1,f2,g1_13a)	
	g1_13b = g1_13a

	g1_1loop = g1_1loop + g1_1a
	g1_1loop = g1_1loop + g1_1b
	g1_1loop = g1_1loop + g1_2
	g1_1loop = g1_1loop + g1_3a
	g1_1loop = g1_1loop + g1_3b
	g1_1loop = g1_1loop + g1_4a
	g1_1loop = g1_1loop + g1_4b
	g1_1loop = g1_1loop + g1_5a
	g1_1loop = g1_1loop + g1_5b
	g1_1loop = g1_1loop + g1_6a
	g1_1loop = g1_1loop + g1_6b
	g1_1loop = g1_1loop + g1_7
	g1_1loop = g1_1loop + g1_8a
	g1_1loop = g1_1loop + g1_8b
	g1_1loop = g1_1loop + g1_9
	g1_1loop = g1_1loop + g1_10a
	g1_1loop = g1_1loop + g1_10b
	g1_1loop = g1_1loop + g1_11a
	g1_1loop = g1_1loop + g1_11b
	g1_1loop = g1_1loop + g1_12a
	g1_1loop = g1_1loop + g1_12b
	g1_1loop = g1_1loop + g1_13a
	g1_1loop = g1_1loop + g1_13b
	
	endif


	if(comp_count.eq.1)then	
	call g1dm(l,f1,f2,g1_dm)	
	call g1zf(l,f1,f2,g1_zf)	
	call g1ds(l,f1,f2,g1_ds)	
	endif







	print*, ""
	if(comp_tree.eq.1)then	
	print*, "g1",flavname(f1f2(f1,f2)),"_0:", real(g1_0), aimag(g1_0)
	endif
	if(comp_1loop.eq.1)then	
	print*, "g1",flavname(f1f2(f1,f2)),"_1a:", real(g1_1a), aimag(g1_1a)
	print*, "g1",flavname(f1f2(f1,f2)),"_1b:", real(g1_1b), aimag(g1_1b)
	print*, "g1",flavname(f1f2(f1,f2)),"_2:", real(g1_2), aimag(g1_2)
	print*, "g1",flavname(f1f2(f1,f2)),"_3a:", real(g1_3a), aimag(g1_3a)
	print*, "g1",flavname(f1f2(f1,f2)),"_3b:", real(g1_3b), aimag(g1_3b)
	print*, "g1",flavname(f1f2(f1,f2)),"_4a:", real(g1_4a), aimag(g1_4a)
	print*, "g1",flavname(f1f2(f1,f2)),"_4b:", real(g1_4b), aimag(g1_4b)
	print*, "g1",flavname(f1f2(f1,f2)),"_5a:", real(g1_5a), aimag(g1_5a)
	print*, "g1",flavname(f1f2(f1,f2)),"_5b:", real(g1_5b), aimag(g1_5b)
	print*, "g1",flavname(f1f2(f1,f2)),"_6a:", real(g1_6a), aimag(g1_6a)
	print*, "g1",flavname(f1f2(f1,f2)),"_6b:", real(g1_6b), aimag(g1_6b)
	print*, "g1",flavname(f1f2(f1,f2)),"_7:", real(g1_7), aimag(g1_7)
	print*, "g1",flavname(f1f2(f1,f2)),"_8a:", real(g1_8a), aimag(g1_8a)
	print*, "g1",flavname(f1f2(f1,f2)),"_8b:", real(g1_8b), aimag(g1_8b)
	print*, "g1",flavname(f1f2(f1,f2)),"_9:", real(g1_9), aimag(g1_9)
	print*, "g1",flavname(f1f2(f1,f2)),"_10a:", real(g1_10a), aimag(g1_10a)
	print*, "g1",flavname(f1f2(f1,f2)),"_10b:", real(g1_10b), aimag(g1_10b)
	print*, "g1",flavname(f1f2(f1,f2)),"_11a:", real(g1_11a), aimag(g1_11a)
	print*, "g1",flavname(f1f2(f1,f2)),"_11b:", real(g1_11b), aimag(g1_11b)
	print*, "g1",flavname(f1f2(f1,f2)),"_12a:", real(g1_12a), aimag(g1_12a)
	print*, "g1",flavname(f1f2(f1,f2)),"_12b:", real(g1_12b), aimag(g1_12b)
	print*, "g1",flavname(f1f2(f1,f2)),"_13a:", real(g1_13a), aimag(g1_13a)
	print*, "g1",flavname(f1f2(f1,f2)),"_13b:", real(g1_13b), aimag(g1_13b)
	print*, "g1",flavname(f1f2(f1,f2)),"_1loop:", real(g1_1loop), aimag(g1_1loop)
	endif

	if(comp_count.eq.1)then	
	print*, "g1",flavname(f1f2(f1,f2)),"_dm:", real(g1_dm), aimag(g1_dm)
	print*, "g1",flavname(f1f2(f1,f2)),"_zf:", real(g1_zf), aimag(g1_zf)
	print*, "g1",flavname(f1f2(f1,f2)),"_ds:", real(g1_ds), aimag(g1_ds)
	endif
	print*, ""
	print*, ""

	END SUBROUTINE compute_g1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE compute_gVt(l,f1,f2)
	!
	! Control subroutine to compute all diagrams contributing to gVt
	!
	USE Parameters
	USE Input
	USE TreeLevel_gVt
	USE OneLoop_gVt

	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8),dimension(c0:c1) :: gVt_0
	complex(kind=8),dimension(c0:c1) :: gVt_1a, gVt_1b
	complex(kind=8),dimension(c0:c1) :: gVt_2
	complex(kind=8),dimension(c0:c1) :: gVt_3a, gVt_3b
	complex(kind=8),dimension(c0:c1) :: gVt_4a, gVt_4b
	complex(kind=8),dimension(c0:c1) :: gVt_5a, gVt_5b
	complex(kind=8),dimension(c0:c1) :: gVt_6a, gVt_6b
	complex(kind=8),dimension(c0:c1) :: gVt_7
	complex(kind=8),dimension(c0:c1) :: gVt_8
	complex(kind=8),dimension(c0:c1) :: gVt_9a, gVt_9b
	complex(kind=8),dimension(c0:c1) :: gVt_10a, gVt_10b
	complex(kind=8),dimension(c0:c1) :: gVt_dm, gVt_zf, gVt_ds

	complex(kind=8),dimension(c0:c1) :: gVt_1loop
	integer :: i

	gVt_1loop(:) = dcmplx(0.0d0,0.0d0)
	
	if(comp_tree.eq.1)then	
	call gVt0(l,f1,f2,gVt_0)		
	endif

	if(comp_1loop.eq.1)then	
	call gVt1a(l,f1,f2,gVt_1a)	
	gVt_1b(:) = gVt_1a(:)
	call gVt2(l,f1,f2,gVt_2)	
	call gVt3a(l,f1,f2,gVt_3a)	
	call gVt3b(l,f1,f2,gVt_3b)	
	call gVt4a(l,f1,f2,gVt_4a)	
	call gVt4b(l,f1,f2,gVt_4b)	
	call gVt5a(l,f1,f2,gVt_5a)
	call gVt5b(l,f1,f2,gVt_5b)
	call gVt6a(l,f1,f2,gVt_6a)
	call gVt6b(l,f1,f2,gVt_6b)
	call gVt7(l,f1,f2,gVt_7)
	call gVt8(l,f1,f2,gVt_8)
	call gVt9a(l,f1,f2,gVt_9a)							
	call gVt9b(l,f1,f2,gVt_9b)
	call gVt10a(l,f1,f2,gVt_10a)	
	call gVt10b(l,f1,f2,gVt_10b)		
	gVt_1loop(:) = gVt_1loop(:) + gVt_1a(:)
	gVt_1loop(:) = gVt_1loop(:) + gVt_1b(:)
	gVt_1loop(:) = gVt_1loop(:) + gVt_2(:)
	gVt_1loop(:) = gVt_1loop(:) + gVt_3a(:)
	gVt_1loop(:) = gVt_1loop(:) + gVt_3b(:)
	gVt_1loop(:) = gVt_1loop(:) + gVt_4a(:)
	gVt_1loop(:) = gVt_1loop(:) + gVt_4b(:)
	gVt_1loop(:) = gVt_1loop(:) + gVt_5a(:)
	gVt_1loop(:) = gVt_1loop(:) + gVt_5b(:)
	gVt_1loop(:) = gVt_1loop(:) + gVt_6a(:)
	gVt_1loop(:) = gVt_1loop(:) + gVt_6b(:)
	gVt_1loop(:) = gVt_1loop(:) + gVt_7(:)
	gVt_1loop(:) = gVt_1loop(:) + gVt_8(:)
	gVt_1loop(:) = gVt_1loop(:) + gVt_9a(:)
	gVt_1loop(:) = gVt_1loop(:) + gVt_9b(:)
	gVt_1loop(:) = gVt_1loop(:) + gVt_10a(:)
	gVt_1loop(:) = gVt_1loop(:) + gVt_10b(:)
	endif


	if(comp_count.eq.1)then	
	call gVtdm(l,f1,f2,gVt_dm)
	call gVtzf(l,f1,f2,gVt_zf)
	call gVtds(l,f1,f2,gVt_ds)
	endif


	
	do i=c0,c1

	print*, ""
	if(comp_tree.eq.1)then	
	print*, "gVt",flavname(f1f2(f1,f2)),"_0:",i, real(gVt_0(i)), aimag(gVt_0(i))
	endif
	if(comp_1loop.eq.1)then	
	print*, "gVt",flavname(f1f2(f1,f2)),"_1a:",i, real(gVt_1a(i)), aimag(gVt_1a(i)) 
	print*, "gVt",flavname(f1f2(f1,f2)),"_1b:",i, real(gVt_1b(i)), aimag(gVt_1b(i))
	print*, "gVt",flavname(f1f2(f1,f2)),"_2:",i, real(gVt_2(i)), aimag(gVt_2(i))
	print*, "gVt",flavname(f1f2(f1,f2)),"_3a:",i, real(gVt_3a(i)), aimag(gVt_3a(i))
	print*, "gVt",flavname(f1f2(f1,f2)),"_3b:",i, real(gVt_3b(i)), aimag(gVt_3b(i))
	print*, "gVt",flavname(f1f2(f1,f2)),"_4a:",i, real(gVt_4a(i)), aimag(gVt_4a(i))
	print*, "gVt",flavname(f1f2(f1,f2)),"_4b:",i, real(gVt_4b(i)), aimag(gVt_4b(i))
	print*, "gVt",flavname(f1f2(f1,f2)),"_5a:",i, real(gVt_5a(i)), aimag(gVt_5a(i))
	print*, "gVt",flavname(f1f2(f1,f2)),"_5b:",i, real(gVt_5b(i)), aimag(gVt_5b(i))
	print*, "gVt",flavname(f1f2(f1,f2)),"_6a:",i, real(gVt_6a(i)), aimag(gVt_6a(i))
	print*, "gVt",flavname(f1f2(f1,f2)),"_6b:",i, real(gVt_6b(i)), aimag(gVt_6b(i))
	print*, "gVt",flavname(f1f2(f1,f2)),"_7:",i, real(gVt_7(i)), aimag(gVt_7(i))
	print*, "gVt",flavname(f1f2(f1,f2)),"_8:",i, real(gVt_8(i)), aimag(gVt_8(i))
	print*, "gVt",flavname(f1f2(f1,f2)),"_9a:",i, real(gVt_9a(i)), aimag(gVt_9a(i))
	print*, "gVt",flavname(f1f2(f1,f2)),"_9b:",i, real(gVt_9b(i)), aimag(gVt_9b(i))
	print*, "gVt",flavname(f1f2(f1,f2)),"_10a:",i, real(gVt_10a(i)), aimag(gVt_10a(i))
	print*, "gVt",flavname(f1f2(f1,f2)),"_10b:",i, real(gVt_10b(i)), aimag(gVt_10b(i))
	print*, "gVt",flavname(f1f2(f1,f2)),"_1loop:",i, real(gVt_1loop(i)), aimag(gVt_1loop(i))
	endif
	if(comp_count.eq.1)then	
	print*, "gVt",flavname(f1f2(f1,f2)),"_dm:",i, real(gVt_dm(i)), aimag(gVt_dm(i))
	print*, "gVt",flavname(f1f2(f1,f2)),"_zf:",i, real(gVt_zf(i)), aimag(gVt_zf(i))
	print*, "gVt",flavname(f1f2(f1,f2)),"_ds:",i, real(gVt_ds(i)), aimag(gVt_ds(i))
	endif
	print*, ""	
	print*, ""

	enddo

	END SUBROUTINE compute_gVt




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE compute_lY(l,g_i,f1,f2)
	!
	! Control subroutine to compute all diagrams contributing to lY
	!
	USE Parameters
	USE Input
	USE TreeLevel_lY
	USE OneLoop_lY

	implicit none

	integer,intent(in) :: l,g_i,f1,f2
	complex(kind=8),dimension(c0:c1) :: lY_0
	complex(kind=8),dimension(c0:c1) :: lY_1a, lY_1b
	complex(kind=8),dimension(c0:c1) :: lY_2
	complex(kind=8),dimension(c0:c1) :: lY_3a, lY_3b
	complex(kind=8),dimension(c0:c1) :: lY_4a, lY_4b
	complex(kind=8),dimension(c0:c1) :: lY_5a, lY_5b
	complex(kind=8),dimension(c0:c1) :: lY_6a, lY_6b
	complex(kind=8),dimension(c0:c1) :: lY_7
	complex(kind=8),dimension(c0:c1) :: lY_1loop
	complex(kind=8),dimension(c0:c1) :: lY_dm, lY_zf, lY_ds

	integer::i

	lY_1loop(:) = dcmplx(0.0d0,0.0d0)
	

	if(comp_tree.eq.1)then		
	call lY0(l,g_i,f1,f2,lY_0)
	endif
	if(comp_1loop.eq.1)then	
	call lY1a(l,g_i,f1,f2,lY_1a)
	lY_1b(:) = lY_1a(:)
	call lY2(l,g_i,f1,f2,lY_2)
	call lY3a(l,g_i,f1,f2,lY_3a)
	call lY3b(l,g_i,f1,f2,lY_3b)
	call lY4a(l,g_i,f1,f2,lY_4a)
	call lY4b(l,g_i,f1,f2,lY_4b)
	call lY5a(l,g_i,f1,f2,lY_5a)
	call lY5b(l,g_i,f1,f2,lY_5b)
	call lY6a(l,g_i,f1,f2,lY_6a)
	call lY6b(l,g_i,f1,f2,lY_6b)
	call lY7(l,g_i,f1,f2,lY_7)

	lY_1loop(:) = 2.0d0*lY_1a(:)
	lY_1loop(:) = lY_1loop(:) + lY_2(:)
	lY_1loop(:) = lY_1loop(:) + lY_3a(:)
	lY_1loop(:) = lY_1loop(:) + lY_3b(:)
	lY_1loop(:) = lY_1loop(:) + lY_4a(:)
	lY_1loop(:) = lY_1loop(:) + lY_4b(:)
	lY_1loop(:) = lY_1loop(:) + lY_5a(:)
	lY_1loop(:) = lY_1loop(:) + lY_5b(:)
	lY_1loop(:) = lY_1loop(:) + lY_6a(:)
	lY_1loop(:) = lY_1loop(:) + lY_6b(:)
	lY_1loop(:) = lY_1loop(:) + lY_7(:)
	endif


	if(comp_count.eq.1)then	
	call lYdm(l,g_i,f1,f2,lY_dm)
	call lYzf(l,g_i,f1,f2,lY_zf)
	call lYds(l,g_i,f1,f2,lY_ds)
	endif
	

	do i=c0,c1

	print*, ""
	if(comp_tree.eq.1)then	
	print*, corrname_lYff(g_i),flavname(f1f2(f1,f2)),"_0:",i, real(lY_0(i)), aimag(lY_0(i))
	endif
	if(comp_1loop.eq.1)then	
	print*, corrname_lYff(g_i),flavname(f1f2(f1,f2)),"_1a:",i, real(lY_1a(i)), aimag(lY_1a(i))
	print*, corrname_lYff(g_i),flavname(f1f2(f1,f2)),"_1a:",i, real(lY_1b(i)), aimag(lY_1b(i))
	print*, corrname_lYff(g_i),flavname(f1f2(f1,f2)),"_2:",i, real(lY_2(i)), aimag(lY_2(i))
	print*, corrname_lYff(g_i),flavname(f1f2(f1,f2)),"_3a:",i, real(lY_3a(i)), aimag(lY_3a(i))
	print*, corrname_lYff(g_i),flavname(f1f2(f1,f2)),"_3b:",i, real(lY_3b(i)), aimag(lY_3b(i))
	print*, corrname_lYff(g_i),flavname(f1f2(f1,f2)),"_4a:",i, real(lY_4a(i)), aimag(lY_4a(i))
	print*, corrname_lYff(g_i),flavname(f1f2(f1,f2)),"_4b:",i, real(lY_4b(i)), aimag(lY_4b(i))
	print*, corrname_lYff(g_i),flavname(f1f2(f1,f2)),"_5a:",i, real(lY_5a(i)), aimag(lY_5a(i))
	print*, corrname_lYff(g_i),flavname(f1f2(f1,f2)),"_5b:",i, real(lY_5b(i)), aimag(lY_5b(i))
	print*, corrname_lYff(g_i),flavname(f1f2(f1,f2)),"_6a:",i, real(lY_6a(i)), aimag(lY_6a(i))
	print*, corrname_lYff(g_i),flavname(f1f2(f1,f2)),"_6b:",i, real(lY_6b(i)), aimag(lY_6b(i))
	print*, corrname_lYff(g_i),flavname(f1f2(f1,f2)),"_7:",i, real(lY_7(i)), aimag(lY_7(i))
	print*, corrname_lYff(g_i),flavname(f1f2(f1,f2)),"_1loop:",i, real(lY_1loop(i)), aimag(lY_1loop(i))
	endif
	if(comp_count.eq.1)then	
	print*, corrname_lYff(g_i),flavname(f1f2(f1,f2)),"_dm:",i, real(lY_dm(i)), aimag(lY_dm(i))
	print*, corrname_lYff(g_i),flavname(f1f2(f1,f2)),"_zf:",i, real(lY_zf(i)), aimag(lY_zf(i))
	print*, corrname_lYff(g_i),flavname(f1f2(f1,f2)),"_ds:",i, real(lY_ds(i)), aimag(lY_ds(i))
	endif
	print*, ""
	print*, ""
	enddo


	END SUBROUTINE compute_lY




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE compute_l1(l,f1,f2)
	!
	! Control subroutine to compute all diagrams contributing to l1
	!
	USE Parameters
	USE Input
	USE TreeLevel_l1
	USE OneLoop_l1

	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8) :: l1_0
	complex(kind=8) :: l1_1a, l1_1b
	complex(kind=8) :: l1_2, l1_7
	complex(kind=8) :: l1_4a, l1_4b
	complex(kind=8) :: l1_3a, l1_3b
	complex(kind=8) :: l1_5a, l1_5b
	complex(kind=8) :: l1_6a, l1_6b
	complex(kind=8) :: l1_8a, l1_8b
	complex(kind=8) :: l1_9
	complex(kind=8) :: l1_10a, l1_10b
	complex(kind=8) :: l1_11a, l1_11b
	complex(kind=8) :: l1_12a, l1_12b
	complex(kind=8) :: l1_13a, l1_13b
	complex(kind=8) :: l1_dm, l1_zf, l1_ds
	complex(kind=8) :: l1_1loop

	l1_1loop = dcmplx(0.0d0,0.0d0) 	

	if(comp_tree.eq.1)then	
	call l10(l,f1,f2,l1_0)	
	endif

	if(comp_1loop.eq.1)then	
	call l11a(l,f1,f2,l1_1a)
	l1_1b = l1_1a	
	call l12(l,f1,f2,l1_2)	
	call l13a(l,f1,f2,l1_3a)
	call l13b(l,f1,f2,l1_3b)
	call l14a(l,f1,f2,l1_4a)
	call l14b(l,f1,f2,l1_4b)
	call l15a(l,f1,f2,l1_5a)
	call l15b(l,f1,f2,l1_5b)
	call l16a(l,f1,f2,l1_6a)							
	call l16b(l,f1,f2,l1_6b)
	call l17(l,f1,f2,l1_7)	
	call l18a(l,f1,f2,l1_8a)	
	l1_8b = l1_8a
	call l19(l,f1,f2,l1_9)	
	call l110a(l,f1,f2,l1_10a)	
	call l110b(l,f1,f2,l1_10b)
	call l111a(l,f1,f2,l1_11a)	
	call l111b(l,f1,f2,l1_11b)
	call l112a(l,f1,f2,l1_12a)	
	call l112b(l,f1,f2,l1_12b)
	call l113a(l,f1,f2,l1_13a)	
	l1_13b = l1_13a
	l1_1loop = l1_1loop + l1_1a
	l1_1loop = l1_1loop + l1_1b
	l1_1loop = l1_1loop + l1_2
	l1_1loop = l1_1loop + l1_3a
	l1_1loop = l1_1loop + l1_3b
	l1_1loop = l1_1loop + l1_4a
	l1_1loop = l1_1loop + l1_4b
	l1_1loop = l1_1loop + l1_5a
	l1_1loop = l1_1loop + l1_5b
	l1_1loop = l1_1loop + l1_6a
	l1_1loop = l1_1loop + l1_6b
	l1_1loop = l1_1loop + l1_7
	l1_1loop = l1_1loop + l1_8a
	l1_1loop = l1_1loop + l1_8b
	l1_1loop = l1_1loop + l1_9
	l1_1loop = l1_1loop + l1_10a
	l1_1loop = l1_1loop + l1_10b
	l1_1loop = l1_1loop + l1_11a
	l1_1loop = l1_1loop + l1_11b
	l1_1loop = l1_1loop + l1_12a
	l1_1loop = l1_1loop + l1_12b
	l1_1loop = l1_1loop + l1_13a
	l1_1loop = l1_1loop + l1_13b

	endif


	if(comp_count.eq.1)then	
	call l1dm(l,f1,f2,l1_dm)	
	call l1zf(l,f1,f2,l1_zf)	
	call l1ds(l,f1,f2,l1_ds)	
	endif


	print*, ""
	if(comp_tree.eq.1)then	
	print*, "l1",flavname(f1f2(f1,f2)),"_0:", real(l1_0), aimag(l1_0)
	endif
	if(comp_1loop.eq.1)then	
	print*, "l1",flavname(f1f2(f1,f2)),"_1a:", real(l1_1a), aimag(l1_1a)
	print*, "l1",flavname(f1f2(f1,f2)),"_1b:", real(l1_1b), aimag(l1_1b)
	print*, "l1",flavname(f1f2(f1,f2)),"_2:", real(l1_2), aimag(l1_2)
	print*, "l1",flavname(f1f2(f1,f2)),"_3a:", real(l1_3a), aimag(l1_3a)
	print*, "l1",flavname(f1f2(f1,f2)),"_3b:", real(l1_3b), aimag(l1_3b)
	print*, "l1",flavname(f1f2(f1,f2)),"_4a:", real(l1_4a), aimag(l1_4a)
	print*, "l1",flavname(f1f2(f1,f2)),"_4b:", real(l1_4b), aimag(l1_4b)
	print*, "l1",flavname(f1f2(f1,f2)),"_5a:", real(l1_5a), aimag(l1_5a)
	print*, "l1",flavname(f1f2(f1,f2)),"_5b:", real(l1_5b), aimag(l1_5b)
	print*, "l1",flavname(f1f2(f1,f2)),"_6a:", real(l1_6a), aimag(l1_6a)
	print*, "l1",flavname(f1f2(f1,f2)),"_6b:", real(l1_6b), aimag(l1_6b)
	print*, "l1",flavname(f1f2(f1,f2)),"_7:", real(l1_7), aimag(l1_7)
	print*, "l1",flavname(f1f2(f1,f2)),"_8a:", real(l1_8a), aimag(l1_8a)
	print*, "l1",flavname(f1f2(f1,f2)),"_8b:", real(l1_8b), aimag(l1_8b)
	print*, "l1",flavname(f1f2(f1,f2)),"_9:", real(l1_9), aimag(l1_9)
	print*, "l1",flavname(f1f2(f1,f2)),"_10a:", real(l1_10a), aimag(l1_10a)
	print*, "l1",flavname(f1f2(f1,f2)),"_10b:", real(l1_10b), aimag(l1_10b)
	print*, "l1",flavname(f1f2(f1,f2)),"_11a:", real(l1_11a), aimag(l1_11a)
	print*, "l1",flavname(f1f2(f1,f2)),"_11b:", real(l1_11b), aimag(l1_11b)
	print*, "l1",flavname(f1f2(f1,f2)),"_12a:", real(l1_12a), aimag(l1_12a)
	print*, "l1",flavname(f1f2(f1,f2)),"_12b:", real(l1_12b), aimag(l1_12b)
	print*, "l1",flavname(f1f2(f1,f2)),"_13a:", real(l1_13a), aimag(l1_13a)
	print*, "l1",flavname(f1f2(f1,f2)),"_13b:", real(l1_13b), aimag(l1_13b)
	print*, "l1",flavname(f1f2(f1,f2)),"_1loop:", real(l1_1loop), aimag(l1_1loop)
	endif
	if(comp_count.eq.1)then	
	print*, "l1",flavname(f1f2(f1,f2)),"_dm:", real(l1_dm), aimag(l1_dm)
	print*, "l1",flavname(f1f2(f1,f2)),"_zf:", real(l1_zf), aimag(l1_zf)
	print*, "l1",flavname(f1f2(f1,f2)),"_ds:", real(l1_ds), aimag(l1_ds)
	endif
	print*, ""
	print*, ""

	END SUBROUTINE compute_l1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE compute_lVt(l,f1,f2)
	!
	! Control subroutine to compute all diagrams contributing to lVt
	!
	USE Parameters
	USE Input
	USE TreeLevel_lVt
	USE OneLoop_lVt

	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8),dimension(c0:c1) ::lVt_0
	complex(kind=8),dimension(c0:c1) :: lVt_1a,lVt_1b
	complex(kind=8),dimension(c0:c1) :: lVt_2
	complex(kind=8),dimension(c0:c1) :: lVt_3a, lVt_3b
	complex(kind=8),dimension(c0:c1) :: lVt_4a, lVt_4b
	complex(kind=8),dimension(c0:c1) :: lVt_5a, lVt_5b
	complex(kind=8),dimension(c0:c1) :: lVt_6a, lVt_6b
	complex(kind=8),dimension(c0:c1) :: lVt_7
	complex(kind=8),dimension(c0:c1) :: lVt_8
	complex(kind=8),dimension(c0:c1) :: lVt_9a, lVt_9b
	complex(kind=8),dimension(c0:c1) :: lVt_10a, lVt_10b
	complex(kind=8),dimension(c0:c1) :: lVt_dm, lVt_zf, lVt_ds
	complex(kind=8),dimension(c0:c1) :: lVt_1loop
	integer :: i

	lVt_1loop(:) = dcmplx(0.0d0,0.0d0)

	if(comp_tree.eq.1)then	
	call lVt0(l,f1,f2,lVt_0)
	endif

	if(comp_1loop.eq.1)then			
	call lVt1a(l,f1,f2,lVt_1a)	
	lVt_1b(:) = lVt_1a(:)
	call lVt2(l,f1,f2,lVt_2)	
	call lVt3a(l,f1,f2,lVt_3a)	
	call lVt3b(l,f1,f2,lVt_3b)
	call lVt4a(l,f1,f2,lVt_4a)	
	call lVt4b(l,f1,f2,lVt_4b)
	call lVt5a(l,f1,f2,lVt_5a)
	call lVt5b(l,f1,f2,lVt_5b)
	call lVt6a(l,f1,f2,lVt_6a)
	call lVt6b(l,f1,f2,lVt_6b)
	call lVt7(l,f1,f2,lVt_7)
	call lVt8(l,f1,f2,lVt_8)
	call lVt9a(l,f1,f2,lVt_9a)							
	call lVt9b(l,f1,f2,lVt_9b)
	call lVt10a(l,f1,f2,lVt_10a)	
	call lVt10b(l,f1,f2,lVt_10b)


	lVt_1loop(:) = lVt_1loop(:) + lVt_1a(:)
	lVt_1loop(:) = lVt_1loop(:) + lVt_1b(:)
	lVt_1loop(:) = lVt_1loop(:) + lVt_2(:)
	lVt_1loop(:) = lVt_1loop(:) + lVt_3a(:)
	lVt_1loop(:) = lVt_1loop(:) + lVt_3b(:)
	lVt_1loop(:) = lVt_1loop(:) + lVt_4a(:)
	lVt_1loop(:) = lVt_1loop(:) + lVt_4b(:)
	lVt_1loop(:) = lVt_1loop(:) + lVt_5a(:)
	lVt_1loop(:) = lVt_1loop(:) + lVt_5b(:)
	lVt_1loop(:) = lVt_1loop(:) + lVt_6a(:)
	lVt_1loop(:) = lVt_1loop(:) + lVt_6b(:)
	lVt_1loop(:) = lVt_1loop(:) + lVt_7(:)
	lVt_1loop(:) = lVt_1loop(:) + lVt_8(:)
	lVt_1loop(:) = lVt_1loop(:) + lVt_9a(:)
	lVt_1loop(:) = lVt_1loop(:) + lVt_9b(:)
	lVt_1loop(:) = lVt_1loop(:) + lVt_10a(:)
	lVt_1loop(:) = lVt_1loop(:) + lVt_10b(:)

	endif


	if(comp_count.eq.1)then	
	call lVtdm(l,f1,f2,lVt_dm)
	call lVtzf(l,f1,f2,lVt_zf)
	call lVtds(l,f1,f2,lVt_ds)
	endif

	
	do i=c0,c1

	print*, ""
	if(comp_tree.eq.1)then	
	print*, "lVt",flavname(f1f2(f1,f2)),"_0:",i, real(lVt_0(i)), aimag(lVt_0(i))
	endif
	if(comp_1loop.eq.1)then	
	print*, "lVt",flavname(f1f2(f1,f2)),"_1a:",i, real(lVt_1a(i)), aimag(lVt_1a(i)) 
	print*, "lVt",flavname(f1f2(f1,f2)),"_1b:",i, real(lVt_1b(i)), aimag(lVt_1b(i))
	print*, "lVt",flavname(f1f2(f1,f2)),"_2:",i, real(lVt_2(i)), aimag(lVt_2(i))
	print*, "lVt",flavname(f1f2(f1,f2)),"_3a:",i, real(lVt_3a(i)), aimag(lVt_3a(i))
	print*, "lVt",flavname(f1f2(f1,f2)),"_3b:",i, real(lVt_3b(i)), aimag(lVt_3b(i))
	print*, "lVt",flavname(f1f2(f1,f2)),"_4a:",i, real(lVt_4a(i)), aimag(lVt_4a(i))
	print*, "lVt",flavname(f1f2(f1,f2)),"_4b:",i, real(lVt_4b(i)), aimag(lVt_4b(i))
	print*, "lVt",flavname(f1f2(f1,f2)),"_5a:",i, real(lVt_5a(i)), aimag(lVt_5a(i))
	print*, "lVt",flavname(f1f2(f1,f2)),"_5b:",i, real(lVt_5b(i)), aimag(lVt_5b(i))
	print*, "lVt",flavname(f1f2(f1,f2)),"_6a:",i, real(lVt_6a(i)), aimag(lVt_6a(i))
	print*, "lVt",flavname(f1f2(f1,f2)),"_6b:",i, real(lVt_6b(i)), aimag(lVt_6b(i))
	print*, "lVt",flavname(f1f2(f1,f2)),"_7:",i, real(lVt_7(i)), aimag(lVt_7(i))
	print*, "lVt",flavname(f1f2(f1,f2)),"_8:",i, real(lVt_8(i)), aimag(lVt_8(i))
	print*, "lVt",flavname(f1f2(f1,f2)),"_9a:",i, real(lVt_9a(i)), aimag(lVt_9a(i))
	print*, "lVt",flavname(f1f2(f1,f2)),"_9b:",i, real(lVt_9b(i)), aimag(lVt_9b(i))
	print*, "lVt",flavname(f1f2(f1,f2)),"_10a:",i, real(lVt_10a(i)), aimag(lVt_10a(i))
	print*, "lVt",flavname(f1f2(f1,f2)),"_10b:",i, real(lVt_10b(i)), aimag(lVt_10b(i))
	print*, "lVt",flavname(f1f2(f1,f2)),"_1loop:",i, real(lVt_1loop(i)), aimag(lVt_1loop(i))
	endif

	if(comp_count.eq.1)then	
	print*, "lVt",flavname(f1f2(f1,f2)),"_dm:",i, real(lVt_dm(i)), aimag(lVt_dm(i))
	print*, "lVt",flavname(f1f2(f1,f2)),"_zf:",i, real(lVt_zf(i)), aimag(lVt_zf(i))
	print*, "lVt",flavname(f1f2(f1,f2)),"_ds:",i, real(lVt_ds(i)), aimag(lVt_ds(i))
	endif

	print*, ""	
	print*, ""

	enddo



	END SUBROUTINE compute_lVt


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



	SUBROUTINE compute_G1_4f(l,g_i1,g_i2,flav)
	!
	! Control subroutine to compute all diagrams contributing to G1_4f
	!
	USE Parameters
	USE Input
	USE TreeLevel_4fermions
	USE OneLoop_4f_G1

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
	complex(kind=8),dimension(c0:c1) :: G14fdm,G14fzf,G14fds
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
	G14fdm(:) = dcmplx(0.0d0,0.0d0)
	G14fzf(:) = dcmplx(0.0d0,0.0d0)
	G14fds(:) = dcmplx(0.0d0,0.0d0)	



	if(comp_tree.eq.1)then	
	call G1_4f_0(l,g_i1,g_i2,flav,G14f0)
	endif

	if(comp_1loop.eq.1)then	
	call G1_4f_1a(l,g_i1,g_i2,flav,G14f1a)
	G14f1b = G14f1a 
	call G1_4f_1c(l,g_i1,g_i2,flav,G14f1c)
	G14f1d = G14f1c 
	call G1_4f_2a(l,g_i1,g_i2,flav,G14f2a)
	call G1_4f_2b(l,g_i1,g_i2,flav,G14f2b)
	call G1_4f_3a(l,g_i1,g_i2,flav,G14f3a)
	call G1_4f_3b(l,g_i1,g_i2,flav,G14f3b)
	call G1_4f_3c(l,g_i1,g_i2,flav,G14f3c)
	call G1_4f_3d(l,g_i1,g_i2,flav,G14f3d)
	call G1_4f_4a(l,g_i1,g_i2,flav,G14f4a)
	call G1_4f_4b(l,g_i1,g_i2,flav,G14f4b)
	call G1_4f_4c(l,g_i1,g_i2,flav,G14f4c)
	call G1_4f_4d(l,g_i1,g_i2,flav,G14f4d)
	call G1_4f_5a(l,g_i1,g_i2,flav,G14f5a)
	call G1_4f_5b(l,g_i1,g_i2,flav,G14f5b)
	call G1_4f_5c(l,g_i1,g_i2,flav,G14f5c)
	call G1_4f_5d(l,g_i1,g_i2,flav,G14f5d)
	call G1_4f_6a(l,g_i1,g_i2,flav,G14f6a)
	call G1_4f_6b(l,g_i1,g_i2,flav,G14f6b)
	call G1_4f_6c(l,g_i1,g_i2,flav,G14f6c)
	call G1_4f_6d(l,g_i1,g_i2,flav,G14f6d)
	call G1_4f_7a(l,g_i1,g_i2,flav,G14f7a)
	call G1_4f_7b(l,g_i1,g_i2,flav,G14f7b)
	endif

	if(comp_count.eq.1)then	
	call G1_4f_dm(l,g_i1,g_i2,flav,G14fdm)
	call G1_4f_zf(l,g_i1,g_i2,flav,G14fzf)
	call G1_4f_ds(l,g_i1,g_i2,flav,G14fds)
	endif

	print*, ""
	if(comp_tree.eq.1)then	
	print*, "G1_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_0=", real(G14f0(c0)), aimag(G14f0(c0))
	endif
	print*, ""
	if(comp_1loop.eq.1)then	
	print*, "G1_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_1a=", real(G14f1a(c0)), aimag(G14f1a(c0))
	print*, "G1_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_1b=", real(G14f1b(c0)), aimag(G14f1b(c0))
	print*, "G1_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_1c=", real(G14f1c(c0)), aimag(G14f1c(c0))
	print*, "G1_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_1d=", real(G14f1d(c0)), aimag(G14f1d(c0))
	print*, "G1_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_2a=", real(G14f2a(c0)), aimag(G14f2a(c0))
	print*, "G1_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_2b=", real(G14f2b(c0)), aimag(G14f2b(c0))
	print*, "G1_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_3a=", real(G14f3a(c0)), aimag(G14f3a(c0))
	print*, "G1_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_3b=", real(G14f3b(c0)), aimag(G14f3b(c0))
	print*, "G1_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_3c=", real(G14f3c(c0)), aimag(G14f3c(c0))
	print*, "G1_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_3d=", real(G14f3d(c0)), aimag(G14f3d(c0))
	print*, "G1_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_4a=", real(G14f4a(c0)), aimag(G14f4a(c0))
	print*, "G1_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_4b=", real(G14f4b(c0)), aimag(G14f4b(c0))
	print*, "G1_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_4c=", real(G14f4c(c0)), aimag(G14f4c(c0))
	print*, "G1_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_4d=", real(G14f4d(c0)), aimag(G14f4d(c0))
	print*, "G1_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_5a=", real(G14f5a(c0)), aimag(G14f5a(c0))
	print*, "G1_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_5b=", real(G14f5b(c0)), aimag(G14f5b(c0))
	print*, "G1_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_5c=", real(G14f5c(c0)), aimag(G14f5c(c0))
	print*, "G1_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_5d=", real(G14f5d(c0)), aimag(G14f5d(c0))
	print*, "G1_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_6a=", real(G14f6a(c0)), aimag(G14f6a(c0))
	print*, "G1_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_6b=", real(G14f6b(c0)), aimag(G14f6b(c0))
	print*, "G1_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_6c=", real(G14f6c(c0)), aimag(G14f6c(c0))
	print*, "G1_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_6d=", real(G14f6d(c0)), aimag(G14f6d(c0))
	print*, "G1_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_7a=", real(G14f7a(c0)), aimag(G14f7a(c0))
	print*, "G1_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_7b=", real(G14f7b(c0)), aimag(G14f7b(c0))

	G14f1loop(:) = G14f1a(:) + G14f1b(:) + G14f1c(:) + G14f1d(:)
	G14f1loop(:) = G14f1loop(:) + G14f2a(:) + G14f2b(:)
	G14f1loop(:) = G14f1loop(:) + G14f3a(:) + G14f3b(:) + G14f3c(:) + G14f3d(:)
	G14f1loop(:) = G14f1loop(:) + G14f4a(:) + G14f4b(:) + G14f4c(:) + G14f4d(:)
	G14f1loop(:) = G14f1loop(:) + G14f5a(:) + G14f5b(:) + G14f5c(:) + G14f5d(:)
	G14f1loop(:) = G14f1loop(:) + G14f6a(:) + G14f6b(:) + G14f6c(:) + G14f6d(:)
	G14f1loop(:) = G14f1loop(:) + G14f7a(:) + G14f7b(:)

	print*, ""
print*, "G1_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_1loop=",real(G14f1loop(c0)),aimag(G14f1loop(c0)) 
	print*, ""


	endif

	if(comp_count.eq.1)then	
	print*, ""
	print*, "G1_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_dm=", real(G14fdm(c0)), aimag(G14fdm(c0))
	print*, "G1_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_zf=", real(G14fzf(c0)), aimag(G14fzf(c0))
	print*, "G1_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_ds=", real(G14fds(c0)), aimag(G14fds(c0))
	endif


	END SUBROUTINE compute_G1_4f





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




	SUBROUTINE compute_G2_4f(l,g_i1,g_i2,flav)
	!
	! Control subroutine to compute all diagrams contributing to G2_4f
	!
	USE Parameters
	USE Input
	USE TreeLevel_4fermions
	USE OneLoop_4f_G2

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
	complex(kind=8),dimension(c0:c1) :: G24f7a,G24f7b,G24f8a,G24f8b
	complex(kind=8),dimension(c0:c1) :: G24f9a,G24f9b
	complex(kind=8),dimension(c0:c1) :: G24f10a,G24f10b,G24f10c,G24f10d
	complex(kind=8),dimension(c0:c1) :: G24f11a,G24f11b,G24f11c,G24f11d
	complex(kind=8),dimension(c0:c1) :: G24f12a,G24f12b,G24f13a,G24f13b
	complex(kind=8),dimension(c0:c1) :: G24fdm,G24fzf,G24fds

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
	G24fdm(:) = dcmplx(0.0d0,0.0d0)
	G24fzf(:) = dcmplx(0.0d0,0.0d0)
	G24fds(:) = dcmplx(0.0d0,0.0d0)	


	
	if(comp_tree.eq.1)then	
	call G2_4f_0(l,g_i1,g_i2,flav,G24f0)
	endif
	if(comp_1loop.eq.1)then	
	call G2_4f_1a(l,g_i1,g_i2,flav,G24f1a)
	G24f1b = G24f1a 
	call G2_4f_1c(l,g_i1,g_i2,flav,G24f1c)
	G24f1d = G24f1c 
	call G2_4f_2a(l,g_i1,g_i2,flav,G24f2a)
	call G2_4f_2b(l,g_i1,g_i2,flav,G24f2b)
	call G2_4f_3a(l,g_i1,g_i2,flav,G24f3a)
	call G2_4f_3b(l,g_i1,g_i2,flav,G24f3b)
	call G2_4f_3c(l,g_i1,g_i2,flav,G24f3c)
	call G2_4f_3d(l,g_i1,g_i2,flav,G24f3d)
	call G2_4f_4a(l,g_i1,g_i2,flav,G24f4a)
	call G2_4f_4b(l,g_i1,g_i2,flav,G24f4b)
	call G2_4f_4c(l,g_i1,g_i2,flav,G24f4c)
	call G2_4f_4d(l,g_i1,g_i2,flav,G24f4d)
	call G2_4f_5a(l,g_i1,g_i2,flav,G24f5a)
	call G2_4f_5b(l,g_i1,g_i2,flav,G24f5b)
	call G2_4f_5c(l,g_i1,g_i2,flav,G24f5c)
	call G2_4f_5d(l,g_i1,g_i2,flav,G24f5d)
	call G2_4f_6a(l,g_i1,g_i2,flav,G24f6a)
	call G2_4f_6b(l,g_i1,g_i2,flav,G24f6b)
	call G2_4f_6c(l,g_i1,g_i2,flav,G24f6c)
	call G2_4f_6d(l,g_i1,g_i2,flav,G24f6d)
	call G2_4f_7a(l,g_i1,g_i2,flav,G24f7a)
	call G2_4f_7b(l,g_i1,g_i2,flav,G24f7b)
	call G2_4f_8a(l,g_i1,g_i2,flav,G24f8a)
	call G2_4f_8b(l,g_i1,g_i2,flav,G24f8b)
	call G2_4f_9a(l,g_i1,g_i2,flav,G24f9a)
	call G2_4f_9b(l,g_i1,g_i2,flav,G24f9b)
	call G2_4f_10a(l,g_i1,g_i2,flav,G24f10a)
	call G2_4f_10b(l,g_i1,g_i2,flav,G24f10b)
	call G2_4f_10c(l,g_i1,g_i2,flav,G24f10c)
	call G2_4f_10d(l,g_i1,g_i2,flav,G24f10d)
	call G2_4f_11a(l,g_i1,g_i2,flav,G24f11a)
	call G2_4f_11b(l,g_i1,g_i2,flav,G24f11b)
	call G2_4f_11c(l,g_i1,g_i2,flav,G24f11c)
	call G2_4f_11d(l,g_i1,g_i2,flav,G24f11d)
	call G2_4f_12a(l,g_i1,g_i2,flav,G24f12a)
	call G2_4f_12b(l,g_i1,g_i2,flav,G24f12b)
	call G2_4f_13a(l,g_i1,g_i2,flav,G24f13a)
	call G2_4f_13b(l,g_i1,g_i2,flav,G24f13b)
	endif

	if(comp_count.eq.1)then	
	call G2_4f_dm(l,g_i1,g_i2,flav,G24fdm)
	call G2_4f_zf(l,g_i1,g_i2,flav,G24fzf)
	call G2_4f_ds(l,g_i1,g_i2,flav,G24fds)
	endif

	print*, ""
	if(comp_tree.eq.1)then	
	print*,  "G2_",corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_0=", real(G24f0(c0)), aimag(G24f0(c0))
	endif
	print*, ""
	if(comp_1loop.eq.1)then	
	print*, "G2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_1a=", real(G24f1a(c0)), aimag(G24f1a(c0))
	print*, "G2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_1b=", real(G24f1b(c0)), aimag(G24f1b(c0))
	print*, "G2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_1c=", real(G24f1c(c0)), aimag(G24f1c(c0))
	print*, "G2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_1d=", real(G24f1d(c0)), aimag(G24f1d(c0))
	print*, "G2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_2a=", real(G24f2a(c0)), aimag(G24f2a(c0))
	print*, "G2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_2b=", real(G24f2b(c0)), aimag(G24f2b(c0))
	print*, "G2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_3a=", real(G24f3a(c0)), aimag(G24f3a(c0))
	print*, "G2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_3b=", real(G24f3b(c0)), aimag(G24f3b(c0))
	print*, "G2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_3c=", real(G24f3c(c0)), aimag(G24f3c(c0))
	print*, "G2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_3d=", real(G24f3d(c0)), aimag(G24f3d(c0))
	print*, "G2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_4a=", real(G24f4a(c0)), aimag(G24f4a(c0))
	print*, "G2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_4b=", real(G24f4b(c0)), aimag(G24f4b(c0))
	print*, "G2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_4c=", real(G24f4c(c0)), aimag(G24f4c(c0))
	print*, "G2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_4d=", real(G24f4d(c0)), aimag(G24f4d(c0))
	print*, "G2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_5a=", real(G24f5a(c0)), aimag(G24f5a(c0))
	print*, "G2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_5b=", real(G24f5b(c0)), aimag(G24f5b(c0))
	print*, "G2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_5c=", real(G24f5c(c0)), aimag(G24f5c(c0))
	print*, "G2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_5d=", real(G24f5d(c0)), aimag(G24f5d(c0))
	print*, "G2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_6a=", real(G24f6a(c0)), aimag(G24f6a(c0))
	print*, "G2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_6b=", real(G24f6b(c0)), aimag(G24f6b(c0))
	print*, "G2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_6c=", real(G24f6c(c0)), aimag(G24f6c(c0))
	print*, "G2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_6d=", real(G24f6d(c0)), aimag(G24f6d(c0))
	print*, "G2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_7a=", real(G24f7a(c0)), aimag(G24f7a(c0))
	print*, "G2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_7b=", real(G24f7b(c0)), aimag(G24f7b(c0))
	print*, "G2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_8a=", real(G24f8a(c0)), aimag(G24f8a(c0))
	print*, "G2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_8b=", real(G24f8b(c0)), aimag(G24f8b(c0))
	print*, "G2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_9a=", real(G24f9a(c0)), aimag(G24f9a(c0))
	print*, "G2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_9b=", real(G24f9b(c0)), aimag(G24f9b(c0))
	print*, "G2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_10a=", real(G24f10a(c0)), aimag(G24f10a(c0))
	print*, "G2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_10b=", real(G24f10b(c0)), aimag(G24f10b(c0))
	print*, "G2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_10c=", real(G24f10c(c0)), aimag(G24f10c(c0))
	print*, "G2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_10d=", real(G24f10d(c0)), aimag(G24f10d(c0))
	print*, "G2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_11a=", real(G24f11a(c0)), aimag(G24f11a(c0))
	print*, "G2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_11b=", real(G24f11b(c0)), aimag(G24f11b(c0))
	print*, "G2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_11c=", real(G24f11c(c0)), aimag(G24f11c(c0))
	print*, "G2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_11d=", real(G24f11d(c0)), aimag(G24f11d(c0))
	print*, "G2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_12a=", real(G24f12a(c0)), aimag(G24f12a(c0))
	print*, "G2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_12b=", real(G24f12b(c0)), aimag(G24f12b(c0))
	print*, "G2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_13a=", real(G24f13a(c0)), aimag(G24f13a(c0))
	print*, "G2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_13b=", real(G24f13b(c0)), aimag(G24f13b(c0))



	G24f1loop(:) = G24f1a(:) + G24f1b(:) + G24f1c(:) + G24f1d(:)
	G24f1loop(:) = G24f1loop(:) + G24f2a(:) + G24f2b(:)
	G24f1loop(:) = G24f1loop(:) + G24f3a(:) + G24f3b(:) + G24f3c(:) + G24f3d(:)
	G24f1loop(:) = G24f1loop(:) + G24f4a(:) + G24f4b(:) + G24f4c(:) + G24f4d(:)
	G24f1loop(:) = G24f1loop(:) + G24f5a(:) + G24f5b(:) + G24f5c(:) + G24f5d(:)
	G24f1loop(:) = G24f1loop(:) + G24f6a(:) + G24f6b(:) + G24f6c(:) + G24f6d(:)
	G24f1loop(:) = G24f1loop(:) + G24f7a(:) + G24f7b(:) + G24f8a(:) + G24f8b(:)
	G24f1loop(:) = G24f1loop(:) + G24f9a(:) + G24f9b(:)
	G24f1loop(:) = G24f1loop(:) + G24f10a(:) + G24f10b(:) + G24f10c(:) + G24f10d(:)
	G24f1loop(:) = G24f1loop(:) + G24f11a(:) + G24f11b(:) + G24f11c(:) + G24f11d(:)
	G24f1loop(:) = G24f1loop(:) + G24f12a(:) + G24f12b(:) + G24f13a(:) + G24f13b(:)

	print*, ""
print*, "G2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_1loop=",real(G24f1loop(c0)),aimag(G24f1loop(c0)) 
	print*, ""

	endif


	if(comp_count.eq.1)then	

	print*, ""
	print*, "G2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_dm=", real(G24fdm(c0)), aimag(G24fdm(c0))
	print*, "G2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_zf=", real(G24fzf(c0)), aimag(G24fzf(c0))
	print*, "G2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_ds=", real(G24fds(c0)), aimag(G24fds(c0))
	endif



	END SUBROUTINE compute_G2_4f




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



	SUBROUTINE compute_Gk1_4f(l,g_i1,g_i2,flav)
	!
	! Control subroutine to compute all diagrams contributing to Gk1_4f
	!
	USE Parameters
	USE Input
	USE TreeLevel_4fermions
	USE OneLoop_4f_Gk1

	implicit none

	integer,intent(in) :: l, g_i1, g_i2		!size, gamma1, gamma2 
	integer,dimension(4),intent(in) :: flav

	complex(kind=8),dimension(c0:c1) :: Gk14f0
	complex(kind=8),dimension(c0:c1) :: Gk14f1a,Gk14f1b,Gk14f1c,Gk14f1d
	complex(kind=8),dimension(c0:c1) :: Gk14f2a,Gk14f2b
	complex(kind=8),dimension(c0:c1) :: Gk14f3a,Gk14f3b,Gk14f3c,Gk14f3d
	complex(kind=8),dimension(c0:c1) :: Gk14f4a,Gk14f4b,Gk14f4c,Gk14f4d
	complex(kind=8),dimension(c0:c1) :: Gk14f5a,Gk14f5b,Gk14f5c,Gk14f5d
	complex(kind=8),dimension(c0:c1) :: Gk14f6a,Gk14f6b,Gk14f6c,Gk14f6d
	complex(kind=8),dimension(c0:c1) :: Gk14f7a,Gk14f7b
	complex(kind=8),dimension(c0:c1) :: Gk14fdm,Gk14fzf,Gk14fds
	complex(kind=8),dimension(c0:c1) :: Gk14f1loop


	Gk14f0(:) = dcmplx(0.0d0,0.0d0)
	Gk14f1a(:) = dcmplx(0.0d0,0.0d0)
	Gk14f1b(:) = dcmplx(0.0d0,0.0d0)
	Gk14f1c(:) = dcmplx(0.0d0,0.0d0)
	Gk14f1d(:) = dcmplx(0.0d0,0.0d0)	
	Gk14f2a(:) = dcmplx(0.0d0,0.0d0)
	Gk14f2b(:) = dcmplx(0.0d0,0.0d0)
	Gk14f3a(:) = dcmplx(0.0d0,0.0d0)
	Gk14f3b(:) = dcmplx(0.0d0,0.0d0)
	Gk14f3c(:) = dcmplx(0.0d0,0.0d0)
	Gk14f3d(:) = dcmplx(0.0d0,0.0d0)
	Gk14f4a(:) = dcmplx(0.0d0,0.0d0)
	Gk14f4b(:) = dcmplx(0.0d0,0.0d0)
	Gk14f4c(:) = dcmplx(0.0d0,0.0d0)
	Gk14f4d(:) = dcmplx(0.0d0,0.0d0)
	Gk14f5a(:) = dcmplx(0.0d0,0.0d0)
	Gk14f5b(:) = dcmplx(0.0d0,0.0d0)
	Gk14f5c(:) = dcmplx(0.0d0,0.0d0)
	Gk14f5d(:) = dcmplx(0.0d0,0.0d0)
	Gk14f6a(:) = dcmplx(0.0d0,0.0d0)
	Gk14f6b(:) = dcmplx(0.0d0,0.0d0)
	Gk14f6c(:) = dcmplx(0.0d0,0.0d0)
	Gk14f6d(:) = dcmplx(0.0d0,0.0d0)
	Gk14f7a(:) = dcmplx(0.0d0,0.0d0)
	Gk14f7b(:) = dcmplx(0.0d0,0.0d0)
	Gk14f1loop(:) = dcmplx(0.0d0,0.0d0)
	Gk14fdm(:) = dcmplx(0.0d0,0.0d0)
	Gk14fzf(:) = dcmplx(0.0d0,0.0d0)
	Gk14fds(:) = dcmplx(0.0d0,0.0d0)	





	if(comp_tree.eq.1)then	
	call Gk1_4f_0(l,g_i1,g_i2,flav,Gk14f0)
	endif

	if(comp_1loop.eq.1)then	
	call Gk1_4f_1a(l,g_i1,g_i2,flav,Gk14f1a)
	Gk14f1b = Gk14f1a 
	call Gk1_4f_1c(l,g_i1,g_i2,flav,Gk14f1c)
	Gk14f1d = Gk14f1c 
	call Gk1_4f_2a(l,g_i1,g_i2,flav,Gk14f2a)
	call Gk1_4f_2b(l,g_i1,g_i2,flav,Gk14f2b)
	call Gk1_4f_3a(l,g_i1,g_i2,flav,Gk14f3a)
	call Gk1_4f_3b(l,g_i1,g_i2,flav,Gk14f3b)
	call Gk1_4f_3c(l,g_i1,g_i2,flav,Gk14f3c)
	call Gk1_4f_3d(l,g_i1,g_i2,flav,Gk14f3d)
	call Gk1_4f_4a(l,g_i1,g_i2,flav,Gk14f4a)
	call Gk1_4f_4b(l,g_i1,g_i2,flav,Gk14f4b)
	call Gk1_4f_4c(l,g_i1,g_i2,flav,Gk14f4c)
	call Gk1_4f_4d(l,g_i1,g_i2,flav,Gk14f4d)
	call Gk1_4f_5a(l,g_i1,g_i2,flav,Gk14f5a)
	call Gk1_4f_5b(l,g_i1,g_i2,flav,Gk14f5b)
	call Gk1_4f_5c(l,g_i1,g_i2,flav,Gk14f5c)
	call Gk1_4f_5d(l,g_i1,g_i2,flav,Gk14f5d)
	call Gk1_4f_6a(l,g_i1,g_i2,flav,Gk14f6a)
	call Gk1_4f_6b(l,g_i1,g_i2,flav,Gk14f6b)
	call Gk1_4f_6c(l,g_i1,g_i2,flav,Gk14f6c)
	call Gk1_4f_6d(l,g_i1,g_i2,flav,Gk14f6d)
	call Gk1_4f_7a(l,g_i1,g_i2,flav,Gk14f7a)
	call Gk1_4f_7b(l,g_i1,g_i2,flav,Gk14f7b)
	endif

	if(comp_count.eq.1)then	
	call Gk1_4f_dm(l,g_i1,g_i2,flav,Gk14fdm)
	call Gk1_4f_zf(l,g_i1,g_i2,flav,Gk14fzf)
	call Gk1_4f_ds(l,g_i1,g_i2,flav,Gk14fds)
	endif



	print*, ""
	if(comp_tree.eq.1)then	
	print*, "Gk1_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_0=", real(Gk14f0(c0)), aimag(Gk14f0(c0))
	endif
	print*, ""

	if(comp_1loop.eq.1)then	
	print*, "Gk1_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_1a=", real(Gk14f1a(c0)), aimag(Gk14f1a(c0))
	print*, "Gk1_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_1b=", real(Gk14f1b(c0)), aimag(Gk14f1b(c0))
	print*, "Gk1_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_1c=", real(Gk14f1c(c0)), aimag(Gk14f1c(c0))
	print*, "Gk1_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_1d=", real(Gk14f1d(c0)), aimag(Gk14f1d(c0))
	print*, "Gk1_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_2a=", real(Gk14f2a(c0)), aimag(Gk14f2a(c0))
	print*, "Gk1_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_2b=", real(Gk14f2b(c0)), aimag(Gk14f2b(c0))
	print*, "Gk1_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_3a=", real(Gk14f3a(c0)), aimag(Gk14f3a(c0))
	print*, "Gk1_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_3b=", real(Gk14f3b(c0)), aimag(Gk14f3b(c0))
	print*, "Gk1_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_3c=", real(Gk14f3c(c0)), aimag(Gk14f3c(c0))
	print*, "Gk1_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_3d=", real(Gk14f3d(c0)), aimag(Gk14f3d(c0))
	print*, "Gk1_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_4a=", real(Gk14f4a(c0)), aimag(Gk14f4a(c0))
	print*, "Gk1_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_4b=", real(Gk14f4b(c0)), aimag(Gk14f4b(c0))
	print*, "Gk1_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_4c=", real(Gk14f4c(c0)), aimag(Gk14f4c(c0))
	print*, "Gk1_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_4d=", real(Gk14f4d(c0)), aimag(Gk14f4d(c0))
	print*, "Gk1_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_5a=", real(Gk14f5a(c0)), aimag(Gk14f5a(c0))
	print*, "Gk1_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_5b=", real(Gk14f5b(c0)), aimag(Gk14f5b(c0))
	print*, "Gk1_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_5c=", real(Gk14f5c(c0)), aimag(Gk14f5c(c0))
	print*, "Gk1_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_5d=", real(Gk14f5d(c0)), aimag(Gk14f5d(c0))
	print*, "Gk1_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_6a=", real(Gk14f6a(c0)), aimag(Gk14f6a(c0))
	print*, "Gk1_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_6b=", real(Gk14f6b(c0)), aimag(Gk14f6b(c0))
	print*, "Gk1_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_6c=", real(Gk14f6c(c0)), aimag(Gk14f6c(c0))
	print*, "Gk1_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_6d=", real(Gk14f6d(c0)), aimag(Gk14f6d(c0))
	print*, "Gk1_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_7a=", real(Gk14f7a(c0)), aimag(Gk14f7a(c0))
	print*, "Gk1_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_7b=", real(Gk14f7b(c0)), aimag(Gk14f7b(c0))



	Gk14f1loop(:) = Gk14f1a(:) + Gk14f1b(:) + Gk14f1c(:) + Gk14f1d(:)
	Gk14f1loop(:) = Gk14f1loop(:) + Gk14f2a(:) + Gk14f2b(:)
	Gk14f1loop(:) = Gk14f1loop(:) + Gk14f3a(:) + Gk14f3b(:) + Gk14f3c(:) + Gk14f3d(:)
	Gk14f1loop(:) = Gk14f1loop(:) + Gk14f4a(:) + Gk14f4b(:) + Gk14f4c(:) + Gk14f4d(:)
	Gk14f1loop(:) = Gk14f1loop(:) + Gk14f5a(:) + Gk14f5b(:) + Gk14f5c(:) + Gk14f5d(:)
	Gk14f1loop(:) = Gk14f1loop(:) + Gk14f6a(:) + Gk14f6b(:) + Gk14f6c(:) + Gk14f6d(:)
	Gk14f1loop(:) = Gk14f1loop(:) + Gk14f7a(:) + Gk14f7b(:)

	print*, ""
print*, "Gk1_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_1loop=",real(Gk14f1loop(c0)),aimag(Gk14f1loop(c0)) 
	print*, ""

	endif

	if(comp_count.eq.1)then	

	print*, ""
	print*, "Gk1_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_dm=", real(Gk14fdm(c0)), aimag(Gk14fdm(c0))
	print*, "Gk1_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_zf=", real(Gk14fzf(c0)), aimag(Gk14fzf(c0))
	print*, "Gk1_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_ds=", real(Gk14fds(c0)), aimag(Gk14fds(c0))

	endif




	END SUBROUTINE compute_Gk1_4f





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



	SUBROUTINE compute_Gk2_4f(l,g_i1,g_i2,flav)
	!
	! Control subroutine to compute all diagrams contributing to Gk2_4f
	!
	USE Parameters
	USE Input
	USE TreeLevel_4fermions
	USE OneLoop_4f_Gk2


	implicit none

	integer,intent(in) :: l, g_i1, g_i2		!size, gamma1, gamma2 
	integer,dimension(4),intent(in) :: flav

	complex(kind=8),dimension(c0:c1) :: Gk24f0
	complex(kind=8),dimension(c0:c1) :: Gk24f1a,Gk24f1b,Gk24f1c,Gk24f1d
	complex(kind=8),dimension(c0:c1) :: Gk24f2a,Gk24f2b
	complex(kind=8),dimension(c0:c1) :: Gk24f3a,Gk24f3b,Gk24f3c,Gk24f3d
	complex(kind=8),dimension(c0:c1) :: Gk24f4a,Gk24f4b,Gk24f4c,Gk24f4d
	complex(kind=8),dimension(c0:c1) :: Gk24f5a,Gk24f5b,Gk24f5c,Gk24f5d
	complex(kind=8),dimension(c0:c1) :: Gk24f6a,Gk24f6b,Gk24f6c,Gk24f6d
	complex(kind=8),dimension(c0:c1) :: Gk24f7a,Gk24f7b,Gk24f8a,Gk24f8b
	complex(kind=8),dimension(c0:c1) :: Gk24f9a,Gk24f9b
	complex(kind=8),dimension(c0:c1) :: Gk24f10a,Gk24f10b,Gk24f10c,Gk24f10d
	complex(kind=8),dimension(c0:c1) :: Gk24f11a,Gk24f11b,Gk24f11c,Gk24f11d
	complex(kind=8),dimension(c0:c1) :: Gk24f12a,Gk24f12b,Gk24f13a,Gk24f13b
	complex(kind=8),dimension(c0:c1) :: Gk24fdm,Gk24fzf,Gk24fds

	complex(kind=8),dimension(c0:c1) :: Gk24f1loop



	Gk24f0(:) = dcmplx(0.0d0,0.0d0)
	Gk24f1a(:) = dcmplx(0.0d0,0.0d0)
	Gk24f1b(:) = dcmplx(0.0d0,0.0d0)
	Gk24f1c(:) = dcmplx(0.0d0,0.0d0)
	Gk24f1d(:) = dcmplx(0.0d0,0.0d0)	
	Gk24f2a(:) = dcmplx(0.0d0,0.0d0)
	Gk24f2b(:) = dcmplx(0.0d0,0.0d0)
	Gk24f3a(:) = dcmplx(0.0d0,0.0d0)
	Gk24f3b(:) = dcmplx(0.0d0,0.0d0)
	Gk24f3c(:) = dcmplx(0.0d0,0.0d0)
	Gk24f3d(:) = dcmplx(0.0d0,0.0d0)
	Gk24f4a(:) = dcmplx(0.0d0,0.0d0)
	Gk24f4b(:) = dcmplx(0.0d0,0.0d0)
	Gk24f4c(:) = dcmplx(0.0d0,0.0d0)
	Gk24f4d(:) = dcmplx(0.0d0,0.0d0)
	Gk24f5a(:) = dcmplx(0.0d0,0.0d0)
	Gk24f5b(:) = dcmplx(0.0d0,0.0d0)
	Gk24f5c(:) = dcmplx(0.0d0,0.0d0)
	Gk24f5d(:) = dcmplx(0.0d0,0.0d0)
	Gk24f6a(:) = dcmplx(0.0d0,0.0d0)
	Gk24f6b(:) = dcmplx(0.0d0,0.0d0)
	Gk24f6c(:) = dcmplx(0.0d0,0.0d0)
	Gk24f6d(:) = dcmplx(0.0d0,0.0d0)
	Gk24f7a(:) = dcmplx(0.0d0,0.0d0)
	Gk24f7b(:) = dcmplx(0.0d0,0.0d0)
	Gk24f8a(:) = dcmplx(0.0d0,0.0d0)
	Gk24f8b(:) = dcmplx(0.0d0,0.0d0)
	Gk24f9a(:) = dcmplx(0.0d0,0.0d0)
	Gk24f9b(:) = dcmplx(0.0d0,0.0d0)
	Gk24f10a(:) = dcmplx(0.0d0,0.0d0)
	Gk24f10b(:) = dcmplx(0.0d0,0.0d0)
	Gk24f10c(:) = dcmplx(0.0d0,0.0d0)
	Gk24f10d(:) = dcmplx(0.0d0,0.0d0)
	Gk24f11a(:) = dcmplx(0.0d0,0.0d0)
	Gk24f11b(:) = dcmplx(0.0d0,0.0d0)
	Gk24f11c(:) = dcmplx(0.0d0,0.0d0)
	Gk24f11d(:) = dcmplx(0.0d0,0.0d0)
	Gk24f12a(:) = dcmplx(0.0d0,0.0d0)
	Gk24f12b(:) = dcmplx(0.0d0,0.0d0)
	Gk24f13a(:) = dcmplx(0.0d0,0.0d0)
	Gk24f13b(:) = dcmplx(0.0d0,0.0d0)
	Gk24f1loop(:) = dcmplx(0.0d0,0.0d0)
	Gk24fds(:) = dcmplx(0.0d0,0.0d0)
	Gk24fzf(:) = dcmplx(0.0d0,0.0d0)
	Gk24fds(:) = dcmplx(0.0d0,0.0d0)

	
	if(comp_tree.eq.1)then	
	call Gk2_4f_0(l,g_i1,g_i2,flav,Gk24f0)
	endif

	if(comp_1loop.eq.1)then	
	call Gk2_4f_1a(l,g_i1,g_i2,flav,Gk24f1a)
	Gk24f1b = Gk24f1a 
	call Gk2_4f_1c(l,g_i1,g_i2,flav,Gk24f1c)
	Gk24f1d = Gk24f1c 
	call Gk2_4f_2a(l,g_i1,g_i2,flav,Gk24f2a)
	call Gk2_4f_2b(l,g_i1,g_i2,flav,Gk24f2b)
	call Gk2_4f_3a(l,g_i1,g_i2,flav,Gk24f3a)
	call Gk2_4f_3b(l,g_i1,g_i2,flav,Gk24f3b)
	call Gk2_4f_3c(l,g_i1,g_i2,flav,Gk24f3c)
	call Gk2_4f_3d(l,g_i1,g_i2,flav,Gk24f3d)
	call Gk2_4f_4a(l,g_i1,g_i2,flav,Gk24f4a)
	call Gk2_4f_4b(l,g_i1,g_i2,flav,Gk24f4b)
	call Gk2_4f_4c(l,g_i1,g_i2,flav,Gk24f4c)
	call Gk2_4f_4d(l,g_i1,g_i2,flav,Gk24f4d)
	call Gk2_4f_5a(l,g_i1,g_i2,flav,Gk24f5a)
	call Gk2_4f_5b(l,g_i1,g_i2,flav,Gk24f5b)
	call Gk2_4f_5c(l,g_i1,g_i2,flav,Gk24f5c)
	call Gk2_4f_5d(l,g_i1,g_i2,flav,Gk24f5d)
	call Gk2_4f_6a(l,g_i1,g_i2,flav,Gk24f6a)
	call Gk2_4f_6b(l,g_i1,g_i2,flav,Gk24f6b)
	call Gk2_4f_6c(l,g_i1,g_i2,flav,Gk24f6c)
	call Gk2_4f_6d(l,g_i1,g_i2,flav,Gk24f6d)
	call Gk2_4f_7a(l,g_i1,g_i2,flav,Gk24f7a)
	call Gk2_4f_7b(l,g_i1,g_i2,flav,Gk24f7b)
	call Gk2_4f_8a(l,g_i1,g_i2,flav,Gk24f8a)
	call Gk2_4f_8b(l,g_i1,g_i2,flav,Gk24f8b)
	call Gk2_4f_9a(l,g_i1,g_i2,flav,Gk24f9a)
	call Gk2_4f_9b(l,g_i1,g_i2,flav,Gk24f9b)
	call Gk2_4f_10a(l,g_i1,g_i2,flav,Gk24f10a)
	call Gk2_4f_10b(l,g_i1,g_i2,flav,Gk24f10b)
	call Gk2_4f_10c(l,g_i1,g_i2,flav,Gk24f10c)
	call Gk2_4f_10d(l,g_i1,g_i2,flav,Gk24f10d)
	call Gk2_4f_11a(l,g_i1,g_i2,flav,Gk24f11a)
	call Gk2_4f_11b(l,g_i1,g_i2,flav,Gk24f11b)
	call Gk2_4f_11c(l,g_i1,g_i2,flav,Gk24f11c)
	call Gk2_4f_11d(l,g_i1,g_i2,flav,Gk24f11d)
	call Gk2_4f_12a(l,g_i1,g_i2,flav,Gk24f12a)
	call Gk2_4f_12b(l,g_i1,g_i2,flav,Gk24f12b)
	call Gk2_4f_13a(l,g_i1,g_i2,flav,Gk24f13a)
	call Gk2_4f_13b(l,g_i1,g_i2,flav,Gk24f13b)
	endif

	if(comp_count.eq.1)then	
	call Gk2_4f_dm(l,g_i1,g_i2,flav,Gk24fdm)
	call Gk2_4f_zf(l,g_i1,g_i2,flav,Gk24fzf)
	call Gk2_4f_ds(l,g_i1,g_i2,flav,Gk24fds)
	endif




	print*, ""
	if(comp_tree.eq.1)then	
	print*,  "Gk2_",corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_0=", real(Gk24f0(c0)), aimag(Gk24f0(c0))
	endif
	print*, ""

	if(comp_1loop.eq.1)then	
	print*, "Gk2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_1a=", real(Gk24f1a(c0)), aimag(Gk24f1a(c0))
	print*, "Gk2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_1b=", real(Gk24f1b(c0)), aimag(Gk24f1b(c0))
	print*, "Gk2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_1c=", real(Gk24f1c(c0)), aimag(Gk24f1c(c0))
	print*, "Gk2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_1d=", real(Gk24f1d(c0)), aimag(Gk24f1d(c0))
	print*, "Gk2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_2a=", real(Gk24f2a(c0)), aimag(Gk24f2a(c0))
	print*, "Gk2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_2b=", real(Gk24f2b(c0)), aimag(Gk24f2b(c0))
	print*, "Gk2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_3a=", real(Gk24f3a(c0)), aimag(Gk24f3a(c0))
	print*, "Gk2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_3b=", real(Gk24f3b(c0)), aimag(Gk24f3b(c0))
	print*, "Gk2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_3c=", real(Gk24f3c(c0)), aimag(Gk24f3c(c0))
	print*, "Gk2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_3d=", real(Gk24f3d(c0)), aimag(Gk24f3d(c0))
	print*, "Gk2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_4a=", real(Gk24f4a(c0)), aimag(Gk24f4a(c0))
	print*, "Gk2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_4b=", real(Gk24f4b(c0)), aimag(Gk24f4b(c0))
	print*, "Gk2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_4c=", real(Gk24f4c(c0)), aimag(Gk24f4c(c0))
	print*, "Gk2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_4d=", real(Gk24f4d(c0)), aimag(Gk24f4d(c0))
	print*, "Gk2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_5a=", real(Gk24f5a(c0)), aimag(Gk24f5a(c0))
	print*, "Gk2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_5b=", real(Gk24f5b(c0)), aimag(Gk24f5b(c0))
	print*, "Gk2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_5c=", real(Gk24f5c(c0)), aimag(Gk24f5c(c0))
	print*, "Gk2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_5d=", real(Gk24f5d(c0)), aimag(Gk24f5d(c0))
	print*, "Gk2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_6a=", real(Gk24f6a(c0)), aimag(Gk24f6a(c0))
	print*, "Gk2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_6b=", real(Gk24f6b(c0)), aimag(Gk24f6b(c0))
	print*, "Gk2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_6c=", real(Gk24f6c(c0)), aimag(Gk24f6c(c0))
	print*, "Gk2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_6d=", real(Gk24f6d(c0)), aimag(Gk24f6d(c0))
	print*, "Gk2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_7a=", real(Gk24f7a(c0)), aimag(Gk24f7a(c0))
	print*, "Gk2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_7b=", real(Gk24f7b(c0)), aimag(Gk24f7b(c0))
	print*, "Gk2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_8a=", real(Gk24f8a(c0)), aimag(Gk24f8a(c0))
	print*, "Gk2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_8b=", real(Gk24f8b(c0)), aimag(Gk24f8b(c0))
	print*, "Gk2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_9a=", real(Gk24f9a(c0)), aimag(Gk24f9a(c0))
	print*, "Gk2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_9b=", real(Gk24f9b(c0)), aimag(Gk24f9b(c0))
	print*, "Gk2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_10a=", real(Gk24f10a(c0)), aimag(Gk24f10a(c0))
	print*, "Gk2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_10b=", real(Gk24f10b(c0)), aimag(Gk24f10b(c0))
	print*, "Gk2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_10c=", real(Gk24f10c(c0)), aimag(Gk24f10c(c0))
	print*, "Gk2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_10d=", real(Gk24f10d(c0)), aimag(Gk24f10d(c0))
	print*, "Gk2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_11a=", real(Gk24f11a(c0)), aimag(Gk24f11a(c0))
	print*, "Gk2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_11b=", real(Gk24f11b(c0)), aimag(Gk24f11b(c0))
	print*, "Gk2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_11c=", real(Gk24f11c(c0)), aimag(Gk24f11c(c0))
	print*, "Gk2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_11d=", real(Gk24f11d(c0)), aimag(Gk24f11d(c0))
	print*, "Gk2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_12a=", real(Gk24f12a(c0)), aimag(Gk24f12a(c0))
	print*, "Gk2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_12b=", real(Gk24f12b(c0)), aimag(Gk24f12b(c0))
	print*, "Gk2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_13a=", real(Gk24f13a(c0)), aimag(Gk24f13a(c0))
	print*, "Gk2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_13b=", real(Gk24f13b(c0)), aimag(Gk24f13b(c0))


	Gk24f1loop(:) = Gk24f1a(:) + Gk24f1b(:) + Gk24f1c(:) + Gk24f1d(:)
	Gk24f1loop(:) = Gk24f1loop(:) + Gk24f2a(:) + Gk24f2b(:)
	Gk24f1loop(:) = Gk24f1loop(:) + Gk24f3a(:) + Gk24f3b(:) + Gk24f3c(:) + Gk24f3d(:)
	Gk24f1loop(:) = Gk24f1loop(:) + Gk24f4a(:) + Gk24f4b(:) + Gk24f4c(:) + Gk24f4d(:)
	Gk24f1loop(:) = Gk24f1loop(:) + Gk24f5a(:) + Gk24f5b(:) + Gk24f5c(:) + Gk24f5d(:)
	Gk24f1loop(:) = Gk24f1loop(:) + Gk24f6a(:) + Gk24f6b(:) + Gk24f6c(:) + Gk24f6d(:)
	Gk24f1loop(:) = Gk24f1loop(:) + Gk24f7a(:) + Gk24f7b(:) + Gk24f8a(:) + Gk24f8b(:)
	Gk24f1loop(:) = Gk24f1loop(:) + Gk24f9a(:) + Gk24f9b(:)
	Gk24f1loop(:) = Gk24f1loop(:) + Gk24f10a(:) + Gk24f10b(:) + Gk24f10c(:) + Gk24f10d(:)
	Gk24f1loop(:) = Gk24f1loop(:) + Gk24f11a(:) + Gk24f11b(:) + Gk24f11c(:) + Gk24f11d(:)
	Gk24f1loop(:) = Gk24f1loop(:) + Gk24f12a(:) + Gk24f12b(:) + Gk24f13a(:) + Gk24f13b(:)

	print*, ""
print*, "Gk2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_1loop=",real(Gk24f1loop(c0)),aimag(Gk24f1loop(c0)) 
	print*, ""


	endif

	if(comp_count.eq.1)then	
	print*, ""
	print*, "Gk2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_dm=", real(Gk24fdm(c0)), aimag(Gk24fdm(c0))
	print*, "Gk2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_zf=", real(Gk24fzf(c0)), aimag(Gk24fzf(c0))
	print*, "Gk2_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_ds=", real(Gk24fds(c0)), aimag(Gk24fds(c0))
	endif





	END SUBROUTINE compute_Gk2_4f





	END MODULE Compute
