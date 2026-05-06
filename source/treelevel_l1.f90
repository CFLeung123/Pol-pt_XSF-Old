	MODULE	TreeLevel_l1


	USE GammaMatrices
	USE Parameters
	USE Propagators

	implicit none

	CONTAINS


	SUBROUTINE l10(l,f1,f2,l1_0)


	implicit none
	integer,intent(in) :: l, f1, f2		!size
	complex(kind=8),intent(out) :: l1_0
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop_0
	complex(kind=8),dimension(3) :: l_temp
	integer,dimension(3) :: nvec
	integer :: f,fp,k

	nvec(:) = 0
	f = f1f2(f1,f2)
	fp = f1f2(f2,f1)

	l1_0 = dcmplx(0.0d0,0.0d0)	
	l_temp(:) = dcmplx(0.0d0,0.0d0)	
	S_prop_0(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	call Fermion_Propagator(l,nvec,S_Prop_0)

	do k = 1,3
	l_temp(k) = 0.5*dimrep*Trace4(Qk(:,:,k,f),S_prop_0(:,:,1,l-1,f2),Qkp(:,:,k,fp),S_prop_0(:,:,l-1,1,f1))
	enddo

	l1_0 = (l_temp(1)+l_temp(2)+l_temp(3))/3.0d0

	END SUBROUTINE l10

	END MODULE TreeLevel_l1
