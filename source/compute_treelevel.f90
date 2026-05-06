
	MODULE Compute

	implicit none

	CONTAINS



	SUBROUTINE Compute_stuff(l,g_i,g_i2,g4f_i,g4f_i2,f1,f2,f3,f4)

	USE Input

	implicit none

	integer,intent(in) :: l
	integer,intent(in) :: g_i,g_i2,g4f_i,g4f_i2, f1, f2, f3 , f4
	integer,dimension(4) :: flav, flavp


	if(comp_gX.eq.1)then
	  call compute_gX(l,g_i,f1,f2)
!	  call compute_gX_p(l,g_i,f1,f2)
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
	flav(1) = f1
	flav(2) = f2
	flav(3) = f3
	flav(4) = f4
!	flavp(1)= f1
!	flavp(2)= f4
!	flavp(3)= f3
!	flavp(4)= f2
	  call compute_G1_4f(l,g4f_i,g4f_i2,flav)
	  call compute_G2_4f(l,g4f_i,g4f_i2,flav)
	  call compute_Gk1_4f(l,g4f_i,g4f_i2,flav)
	  call compute_Gk2_4f(l,g4f_i,g4f_i2,flav)
	endif



	END SUBROUTINE Compute_stuff



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE compute_gX(l,g_i,f1,f2)

	USE Parameters
	USE TreeLevel_gX

	implicit none

	integer,intent(in) :: l,g_i,f1,f2
	complex(kind=8),dimension(c0:c1) :: gX_0

	integer :: i

	
	call gX0(l,g_i,f1,f2,gX_0)

	do i=c0,c1

	print*, ""
	print*, corrname_gXff(g_i),flavname(f1f2(f1,f2)),"_0:",i, real(gX_0(i)), aimag(gX_0(i))
	print*, ""
	print*, ""
	enddo



	END SUBROUTINE compute_gX




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!	SUBROUTINE compute_gX_p(l,g_i,f1,f2)

!	USE Parameters
!	USE TreeLevel_gX_p

!	implicit none

!	integer,intent(in) :: l,g_i,f1,f2
!	complex(kind=8),dimension(c0:c1) :: gX_0_p

!	integer :: i

	
!	call gX0_p(l,g_i,f1,f2,gX_0_p)

!	do i=c0,c1

!	print*, ""
!	print*, corrname_gXff(g_i),flavname(f1f2(f1,f2)),"_0_p:",i, real(gX_0_p(i)), aimag(gX_0_p(i))
!	print*, ""
!	print*, ""
!	enddo



!	END SUBROUTINE compute_gX_p




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!








!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE compute_g1(l,f1,f2)

	USE Parameters
	USE TreeLevel_g1

	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8) :: g1_0

	call g10(l,f1,f2,g1_0)	

	print*, ""
	print*, "g1",flavname(f1f2(f1,f2)),"_0:", real(g1_0), aimag(g1_0)
	print*, ""
	print*, ""

	END SUBROUTINE compute_g1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE compute_gVt(l,f1,f2)

	USE Parameters
	USE TreeLevel_gVt

	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8),dimension(c0:c1) :: gVt_0
	integer :: i

	
	call gVt0(l,f1,f2,gVt_0)		
	
	do i=c0,c1

	print*, ""
	print*, "gVt",flavname(f1f2(f1,f2)),"_0:",i, real(gVt_0(i)), aimag(gVt_0(i))
	print*, ""	
	print*, ""

	enddo

	END SUBROUTINE compute_gVt




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE compute_lY(l,g_i,f1,f2)

	USE Parameters
	USE TreeLevel_lY

	implicit none

	integer,intent(in) :: l,g_i,f1,f2
	complex(kind=8),dimension(c0:c1) :: lY_0

	integer::i

	
	call lY0(l,g_i,f1,f2,lY_0)

	do i=c0,c1

	print*, ""
	print*, corrname_lYff(g_i),flavname(f1f2(f1,f2)),"_0:",i, real(lY_0(i)), aimag(lY_0(i))
	print*, ""
	print*, ""
	enddo


	END SUBROUTINE compute_lY




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE compute_l1(l,f1,f2)

	USE Parameters
	USE TreeLevel_l1

	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8) :: l1_0

	call l10(l,f1,f2,l1_0)	

	print*, ""
	print*, "l1",flavname(f1f2(f1,f2)),"_0:", real(l1_0), aimag(l1_0)
	print*, ""
	print*, ""

	END SUBROUTINE compute_l1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE compute_lVt(l,f1,f2)

	USE Parameters
	USE TreeLevel_lVt


	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8),dimension(c0:c1) ::lVt_0
	integer :: i


	call lVt0(l,f1,f2,lVt_0)		
	
	do i=c0,c1

	print*, ""
	print*, "lVt",flavname(f1f2(f1,f2)),"_0:",i, real(lVt_0(i)), aimag(lVt_0(i))
	print*, ""	
	print*, ""

	enddo


	END SUBROUTINE compute_lVt


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



	SUBROUTINE compute_G1_4f(l,g_i1,g_i2,flav)

	USE Parameters
	USE TreeLevel_4fermions

	implicit none

	integer,intent(in) :: l, g_i1, g_i2		!size, gamma1, gamma2 
	integer,dimension(4),intent(in) :: flav

	complex(kind=8),dimension(c0:c1) :: G14f0

	G14f0(:) = dcmplx(0.0d0,0.0d0)
	
	call G1_4f_0(l,g_i1,g_i2,flav,G14f0)

	print*, ""
	print*, "G1_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_0=", real(G14f0(c0)), aimag(G14f0(c0))

	END SUBROUTINE compute_G1_4f





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




	SUBROUTINE compute_G2_4f(l,g_i1,g_i2,flav)

	USE Parameters
	USE TreeLevel_4fermions

	implicit none

	integer,intent(in) :: l, g_i1, g_i2		!size, gamma1, gamma2 
	integer,dimension(4),intent(in) :: flav

	complex(kind=8),dimension(c0:c1) :: G24f0

	G24f0(:) = dcmplx(0.0d0,0.0d0)
	
	call G2_4f_0(l,g_i1,g_i2,flav,G24f0)

	print*, ""
	print*,  "G2_",corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_0=", real(G24f0(c0)), aimag(G24f0(c0))
	print*, ""


	END SUBROUTINE compute_G2_4f




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



	SUBROUTINE compute_Gk1_4f(l,g_i1,g_i2,flav)

	USE Parameters
	USE TreeLevel_4fermions

	implicit none

	integer,intent(in) :: l, g_i1, g_i2		!size, gamma1, gamma2 
	integer,dimension(4),intent(in) :: flav

	complex(kind=8),dimension(c0:c1) :: Gk14f0

	Gk14f0(:) = dcmplx(0.0d0,0.0d0)
	call Gk1_4f_0(l,g_i1,g_i2,flav,Gk14f0)
	print*, ""
	print*, "Gk1_", corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_0=", real(Gk14f0(c0)), aimag(Gk14f0(c0))

	END SUBROUTINE compute_Gk1_4f





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




	SUBROUTINE compute_Gk2_4f(l,g_i1,g_i2,flav)

	USE Parameters
	USE TreeLevel_4fermions

	implicit none

	integer,intent(in) :: l, g_i1, g_i2		!size, gamma1, gamma2 
	integer,dimension(4),intent(in) :: flav

	complex(kind=8),dimension(c0:c1) :: Gk24f0

	Gk24f0(:) = dcmplx(0.0d0,0.0d0)

	call Gk2_4f_0(l,g_i1,g_i2,flav,Gk24f0)
	print*, ""
	print*,  "Gk2_",corrname_4f(g_i1),corrname_4f(g_i2),flavname(f1f2(flav(1),flav(2))),flavname(f1f2(flav(3),flav(4))),"_0=", real(Gk24f0(c0)), aimag(Gk24f0(c0))
	print*, ""


	END SUBROUTINE compute_Gk2_4f



	END MODULE Compute
