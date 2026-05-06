
	MODULE Compute

	implicit none

	CONTAINS



	SUBROUTINE Compute_stuff(l)

	implicit none

	integer,intent(in) :: l
!	integer :: g_i, f1, f2



	!call compute_gX(l,1,1,1)
!	call compute_gX(l,1,1,2)

!	call compute_lY(l,2,1,1)
!	call compute_lY(l,2,1,2)


!	call compute_g1(l,1,1)
!	call compute_g1(l,1,2)
!	call compute_g1(l,2,1)
!	call compute_g1(l,2,2)
!	call compute_l1(l,1,1)
!	call compute_l1(l,1,2)


!	call compute_gVt(l,1,1)
!	call compute_gVt(l,1,2)

	call compute_lVt(l,1,1)
!	call compute_lVt(l,1,2)

	END SUBROUTINE Compute_stuff


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE compute_lVt(l,f1,f2)

	USE Parameters
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

	call lVt0(l,f1,f2,lVt_0)		
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
	call lVtdm(l,f1,f2,lVt_dm)
	call lVtzf(l,f1,f2,lVt_zf)
	call lVtds(l,f1,f2,lVt_ds)

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
	
	do i=c0,c1

	print*, ""
	print*, "lVt",flavname(f1f2(f1,f2)),"_0:",i, real(lVt_0(i)), aimag(lVt_0(i))
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
	print*, "lVt",flavname(f1f2(f1,f2)),"_dm:",i, real(lVt_dm(i)), aimag(lVt_dm(i))
	print*, "lVt",flavname(f1f2(f1,f2)),"_zf:",i, real(lVt_zf(i)), aimag(lVt_zf(i))
	print*, "lVt",flavname(f1f2(f1,f2)),"_ds:",i, real(lVt_ds(i)), aimag(lVt_ds(i))
	print*, "lVt",flavname(f1f2(f1,f2)),"_1loop:",i, real(lVt_1loop(i)), aimag(lVt_1a(i))
	print*, ""	
	print*, ""

	enddo


	END SUBROUTINE compute_lVt


	END MODULE Compute
