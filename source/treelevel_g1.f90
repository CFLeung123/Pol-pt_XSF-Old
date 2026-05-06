	MODULE	TreeLevel_g1


	USE GammaMatrices
	USE Parameters
	USE Propagators

	implicit none

	CONTAINS


	SUBROUTINE g10(l,f1,f2,g1_0)


	implicit none
	integer,intent(in) :: l, f1, f2		!size
	complex(kind=8),intent(out) :: g1_0
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop_0
	integer,dimension(3) :: nvec
	integer :: f,fp

	nvec(:) = 0
	f = f1f2(f2,f1)
	fp = f1f2(f1,f2)

	g1_0 = dcmplx(0.0d0,0.0d0)	
	S_prop_0(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	call Fermion_Propagator(l,nvec,S_Prop_0)

	g1_0 = 0.5*dimrep*Trace4(Q5(:,:,f),S_prop_0(:,:,1,l-1,f1),Q5p(:,:,fp),S_prop_0(:,:,l-1,1,f2))


	END SUBROUTINE g10

	END MODULE TreeLevel_g1
