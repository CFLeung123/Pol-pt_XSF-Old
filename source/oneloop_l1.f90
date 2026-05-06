


	MODULE OneLoop_l1
	!
	! Contains the set of routines for the calculaton of the 1-loop
	! diagrams contributing to the l1 correlation functions (including
	! the contribution of counterterms.
	!
	USE GammaMatrices
	USE Input
	USE Parameters
	USE Vertices
	USE Propagators
	USE TreeLevel_l1
	USE omp_lib

	implicit none

	CONTAINS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE l11a(l,f1,f2,l1_1a)
	!
	! Subroutine to compute diagram l1_1a 
	!
	implicit none

	integer,intent(in) :: l, f1, f2
	complex(kind=8),intent(out) :: l1_1a



	integer :: n1,n2,n3
	integer,dimension(3) :: qloop
	complex(kind=8) :: loop
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop!,K_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop
	complex(kind=8) :: l1_0
	complex(kind=8),dimension(3) :: l_temp
	real(kind=8) :: multi
	integer,dimension(3) :: nvec
	integer :: f,fp,k

	f = f1f2(f1,f2)
	fp = f1f2(f2,f1)
	nvec(:) = 0

	!initialization
	l1_1a = dcmplx(0.0d0,0.0d0)
	loop = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	l1_0 = dcmplx(0.0d0,0.0d0)
	l_temp(:) = dcmplx(0.0d0,0.0d0)



	! Construct tree-level part:
	call Fermion_Propagator(l,nvec,S_Prop)
	do k = 1,3
	l_temp(k) = 0.5*dimrep*Trace4(Qk(:,:,k,f),S_prop(:,:,1,l-1,f2),Qkp(:,:,k,fp),S_prop(:,:,l-1,1,f1))
	enddo
	l1_0 = (l_temp(1)+l_temp(2)+l_temp(3))/3.0d0

	!Construct the loop
	do n1=0,l-1			
	do n2=n_limit(mom_deg,n1),l-1
	do n3=n_limit(mom_deg,n2),l-1


	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3

	   multi = 1.0d0*multiplicity(n1,n2,n3)

  	   D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	   call Gluon_propagator(l,qloop,D_prop)	
	   loop = loop + D_prop(0,0,0,0)*multi

	enddo
	enddo
	enddo

	l1_1a = -0.5d0*C2_R*loop*l1_0/((l*1.0d0)**3)

	END SUBROUTINE l11a

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE l12(l,f1,f2,l1_2)
	!
	! Subroutine to compute diagram l1_2
	!
	implicit none

	integer,intent(in) :: l, f1, f2
	complex(kind=8),intent(out) :: l1_2

	integer,dimension(3) :: qloop
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop!,K_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop
	complex(kind=8) :: l1_0
	complex(kind=8),dimension(3) :: l_temp
	integer,dimension(3) :: nvec
	integer :: f, fp, k

	f = f1f2(f1,f2)
	fp = f1f2(f2,f1)
	nvec(:) = 0	

	!initialization
	l1_2 = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	l1_0 = dcmplx(0.0d0,0.0d0)
	l_temp(:) = dcmplx(0.0d0,0.0d0)

	! Construct tree-level part:
	call Fermion_Propagator(l,nvec,S_Prop)
	do k = 1,3
	l_temp(k) = 0.5*dimrep*Trace4(Qk(:,:,k,f),S_prop(:,:,1,l-1,f2),Qkp(:,:,k,fp),S_prop(:,:,l-1,1,f1))
	enddo
	l1_0 = (l_temp(1)+l_temp(2)+l_temp(3))/3.0d0

	!Construct the loop
	qloop(1)=0
	qloop(2)=0
	qloop(3)=0

  	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	call Gluon_propagator(l,qloop,D_prop)	

	l1_2 = C2_R*l1_0*D_prop(0,0,0,0)/((l*1.0d0)**3)

	END SUBROUTINE l12

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE l13a(l,f1,f2,l1_3a)
	!
	! Subroutine to compute diagram l1_3a 
	!
	implicit none

	integer,intent(in) :: l, f1, f2
	complex(kind=8),intent(out) :: l1_3a



	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,0:l,3) :: Qk_S_Qkp_S
	complex(kind=8),dimension(4,4) :: V1mu	
	complex(kind=8),dimension(4,4) :: mat1	
	integer,dimension(3) :: qloop, p_ext
	complex(kind=8),dimension(3) :: l_loop
	integer ::  t0, u0, s0, mu, i
	integer :: n1, n2, n3
	integer :: t0_a, t0_b, u0_a, u0_b
	real(kind=8) :: multi
	complex(kind=8) :: trace
	integer :: f, fp,k

	f = f1f2(f1,f2)
	fp = f1f2(f2,f1)

	!initialization
	l1_3a = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	Qk_S_Qkp_S(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	p_ext(:) = 0
	trace = dcmplx(0.0d0,0.0d0)
	l_loop(:) = dcmplx(0.0d0,0.0d0)


	!part outside the loop
	call Fermion_Propagator(l,p_ext,S_prop_0)
	do k=1,3
	mat1 = matmul(matmul(Qk(:,:,k,f),S_prop_0(:,:,1,l-1,f2)),Qkp(:,:,k,fp))
	do s0=0,l
	  Qk_S_Qkp_S(:,:,s0,k) = matmul(mat1,S_prop_0(:,:,l-1,s0,f1))
	enddo
	enddo

	!Construct  loop
	do n1=0,l-1			
	do n2=n_limit(mom_deg,n1),l-1
	do n3=n_limit(mom_deg,n2),l-1


	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3

	   multi = 1.0d0*multiplicity(n1,n2,n3)
  	   D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	   S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	   call Gluon_propagator(l,qloop,D_prop)	
	   call Fermion_Propagator(l,p_ext+qloop,S_prop)

	do mu=mu_low(0),mu_up(0)

	do s0=0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0=t0_a,t0_b
	do u0=u0_a,u0_b

	call V1(l,p_ext,-p_ext-qloop,qloop,s0,t0,u0,mu,V1mu)

	do k=1,3
	trace = Trace3(Qk_S_Qkp_S(:,:,s0,k),V1mu,S_prop(:,:,t0,1,f1))
	l_loop(k) = l_loop(k) + D_prop(0,mu,0,u0)*trace*multi
	enddo

	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu


	enddo	!q3
	enddo	!q2
	enddo	!q1

	l1_3a = -C2_R*dimrep*(l_loop(1)+l_loop(2)+l_loop(3))/(6.0d0*l**3)


	END SUBROUTINE l13a




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE l13b(l,f1,f2,l1_3b)
	!
	! Subroutine to compute diagram l1_3b 
	!
	implicit none

	integer,intent(in) :: l, f1, f2
	complex(kind=8),intent(out) :: l1_3b



	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,0:l,3) :: S_Qkp_S_Qk
	complex(kind=8),dimension(4,4) :: V1mu	
	complex(kind=8),dimension(4,4) :: mat1	
	integer,dimension(3) :: qloop, p_ext
	complex(kind=8),dimension(3) :: l_loop
	integer ::  t0, u0, s0, mu, i
	integer :: n1, n2, n3
	integer :: t0_a, t0_b, u0_a, u0_b
	real(kind=8) :: multi
	complex(kind=8) :: trace
	integer :: f, fp, k

	f = f1f2(f1,f2)
	fp = f1f2(f2,f1)

	!initialization
	l1_3b = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_Qkp_S_Qk(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	p_ext(:) = 0
	trace = dcmplx(0.0d0,0.0d0)
	l_loop(:) = dcmplx(0.0d0,0.0d0)


	!part outside the loop
	call Fermion_Propagator(l,p_ext,S_prop_0)
	do k=1,3
	mat1 = matmul(matmul(Qkp(:,:,k,fp),S_prop_0(:,:,l-1,1,f1)),Qk(:,:,k,f))
	do t0=0,l
	  S_Qkp_S_Qk(:,:,t0,k) = matmul(S_prop_0(:,:,t0,l-1,f2),mat1)
	enddo
	enddo

	!Construct  loop
	do n1=0,l-1			
	do n2=n_limit(mom_deg,n1),l-1
	do n3=n_limit(mom_deg,n2),l-1


	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3

	   multi = 1.0d0*multiplicity(n1,n2,n3)
  	   D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	   S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	   call Gluon_propagator(l,qloop,D_prop)	
	   call Fermion_Propagator(l,p_ext+qloop,S_prop)

	do mu=mu_low(0),mu_up(0)

	do s0=0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0=t0_a,t0_b
	do u0=u0_a,u0_b

	call V1(l,p_ext+qloop,-p_ext,-qloop,s0,t0,u0,mu,V1mu)

	do k=1,3
	trace = Trace3(S_Qkp_S_Qk(:,:,t0,k),S_prop(:,:,1,s0,f2),V1mu)
	l_loop(k) = l_loop(k) + D_prop(mu,0,u0,0)*trace*multi
	enddo


	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu


	enddo	!q3
	enddo	!q2
	enddo	!q1

	l1_3b = C2_R*dimrep*(l_loop(1)+l_loop(2)+l_loop(3))/(6.0d0*l**3)


	END SUBROUTINE l13b
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE l14a(l,f1,f2,l1_4a)
	!
	! Subroutine to compute diagram l1_4a 
	!
	implicit none

	integer,intent(in) :: l, f1, f2
	complex(kind=8),intent(out) :: l1_4a



	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop_0
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop, p_ext
	complex(kind=8),dimension(3) :: l_loop
	integer ::  t0, u0, s0, mu, i
	integer :: n1, n2, n3
	integer :: t0_a, t0_b, u0_a, u0_b
	complex(kind=8) :: trace
	integer :: f, fp, k

	f = f1f2(f1,f2)
	fp = f1f2(f2,f1)

	!initialization
	l1_4a = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	p_ext(:) = 0
	qloop(:) = 0
	trace = dcmplx(0.0d0,0.0d0)
	l_loop(:) = dcmplx(0.0d0,0.0d0)


	!part outside the loop
	call Fermion_Propagator(l,p_ext,S_prop_0)
	call Gluon_propagator(l,qloop,D_prop)	

	do mu=mu_low(0),mu_up(0)

	do s0=0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0=t0_a,t0_b
	do u0=u0_a,u0_b

	call V1(l,p_ext,-p_ext-qloop,qloop,s0,t0,u0,mu,V1mu)

	do k=1,3
	trace = Trace6(Qk(:,:,k,f),S_prop_0(:,:,1,s0,f2),V1mu,S_prop_0(:,:,t0,l-1,f2),Qkp(:,:,k,fp),S_prop_0(:,:,l-1,1,f1))
	l_loop(k) = l_loop(k) + D_prop(0,mu,0,u0)*trace
	enddo

	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu

	l1_4a = -C2_R*dimrep*(l_loop(1)+l_loop(2)+l_loop(3))/(6.0d0*l**3)


	END SUBROUTINE l14a




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE l14b(l,f1,f2,l1_4b)
	!
	! Subroutine to compute diagram l1_4b 
	!
	implicit none

	integer,intent(in) :: l, f1, f2
	complex(kind=8),intent(out) :: l1_4b



	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) ::  S_prop_0
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop, p_ext
	complex(kind=8),dimension(3) :: l_loop
	integer ::  t0, u0, s0, mu, i
	integer :: n1, n2, n3
	integer :: t0_a, t0_b, u0_a, u0_b
	complex(kind=8) :: trace
	integer :: f, fp, k

	f = f1f2(f1,f2)
	fp = f1f2(f2,f1)

	!initialization
	l1_4b = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	p_ext(:) = 0
	qloop(:) = 0
	trace = dcmplx(0.0d0,0.0d0)
	l_loop(:) = dcmplx(0.0d0,0.0d0)


	!part outside the loop
	call Fermion_Propagator(l,p_ext,S_prop_0)
  	call Gluon_propagator(l,qloop,D_prop)	
	
	do mu=mu_low(0),mu_up(0)

	do s0=0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0=t0_a,t0_b
	do u0=u0_a,u0_b

	call V1(l,p_ext+qloop,-p_ext,-qloop,s0,t0,u0,mu,V1mu)

	do k=1,3
	trace = Trace6(Qk(:,:,k,f),S_prop_0(:,:,1,l-1,f2),Qkp(:,:,k,fp),S_prop_0(:,:,l-1,s0,f1),V1mu,S_prop_0(:,:,t0,1,f1))
	l_loop(k) = l_loop(k) + D_prop(mu,0,u0,0)*trace
	enddo

	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu

	l1_4b = C2_R*dimrep*(l_loop(1)+l_loop(2)+l_loop(3))/(6.0d0*l**3)


	END SUBROUTINE l14b
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE l15a(l,f1,f2,l1_5a)
	!
	! Subroutine to compute diagram l1_5a 
	!
	implicit none

	integer,intent(in) :: l, f1, f2
	complex(kind=8),intent(out) :: l1_5a

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,0:l,0:l,3) :: SQSQS
	complex(kind=8),dimension(4,4) :: V1mu, V1nu	
	complex(kind=8),dimension(4,4) :: mat1	
	integer,dimension(3) :: qloop, p_ext
	complex(kind=8),dimension(3) :: l_loop
	integer :: mu, nu, i
	integer :: n1, n2, n3
	integer :: s0,s0p,t0,t0p,u0,u0p
	integer :: t0_a,t0_b,t0p_a,t0p_b,u0_a,u0_b,u0p_a,u0p_b
	real(kind=8) :: multi
	complex(kind=8) :: trace
	integer :: f, fp, k

	!Variables for the openMP paralelization
	integer :: chunk
	chunk = 1

	f = f1f2(f1,f2)
	fp = f1f2(f2,f1)

	!initialization
	l1_5a = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	SQSQS(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	p_ext(:) = 0
	trace = dcmplx(0.0d0,0.0d0)
	l_loop(:) = dcmplx(0.0d0,0.0d0)


	!part outside the loop
	call Fermion_Propagator(l,p_ext,S_prop_0)
	do k=1,3
	mat1 = matmul(matmul(Qk(:,:,k,f),S_prop_0(:,:,1,l-1,f2)),Qkp(:,:,k,fp))
	do t0=0,l
	do s0=0,l
	  SQSQS(:,:,t0,s0,k) = matmul(matmul(S_prop_0(:,:,t0,1,f1),mat1),S_prop_0(:,:,l-1,s0,f1))
	enddo
	enddo
	enddo

	!Construct  loop

!$OMP PARALLEL DO &
!$OMP DEFAULT(PRIVATE) SHARED(l,p_ext,chunk,f1,f2,SQSQS,multiplicity,mom_deg) &
!$OMP REDUCTION(+:l_loop) &
!$OMP SCHEDULE(DYNAMIC,chunk) &
!$OMP COLLAPSE(1)
	do n1=0,l-1			
	do n2=n_limit(mom_deg,n1),l-1
	do n3=n_limit(mom_deg,n2),l-1


	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3

	multi = 1.0d0*multiplicity(n1,n2,n3)
  	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	call Gluon_propagator(l,qloop,D_prop)	
	call Fermion_Propagator(l,p_ext+qloop,S_prop)

	trace = dcmplx(0.0d0,0.0d0)

	do s0=0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do s0p=0,l
	call loop_limits(l,s0p,t0p_a,t0p_b,u0p_a,u0p_b)
	do t0=t0_a,t0_b
	do t0p=t0p_a,t0p_b		
	do u0=u0_a,u0_b		
	do u0p=u0p_a,u0p_b

	do mu=0,3	
	do nu=mu_low(mu),mu_up(mu)

	call V1(l,p_ext,-p_ext-qloop,qloop,s0,t0,u0,mu,V1mu)
	call V1(l,p_ext+qloop,-p_ext,-qloop,s0p,t0p,u0p,nu,V1nu)

	do k=1,3
	trace = Trace4(SQSQS(:,:,t0p,s0,k),V1mu(:,:),S_prop(:,:,t0,s0p,f1),V1nu(:,:))
	l_loop(k) = l_loop(k) + D_prop(nu,mu,u0p,u0)*trace*multi
	enddo
	

	enddo	!u0
	enddo	!t0
	enddo	!s0
	enddo 
	enddo
	enddo

	enddo	!nu
	enddo	!mu


	enddo	!q3
	enddo	!q2
	enddo	!q1

!$OMP END PARALLEL DO

	l1_5a = -C2_R*dimrep*(l_loop(1)+l_loop(2)+l_loop(3))/(6.0d0*l**3)


	END SUBROUTINE l15a




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE l15b(l,f1,f2,l1_5b)
	!
	! Subroutine to compute diagram l1_5b 
	!
	implicit none

	integer,intent(in) :: l, f1, f2
	complex(kind=8),intent(out) :: l1_5b

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,0:l,0:l,3) :: SQSQS
	complex(kind=8),dimension(4,4) :: V1mu, V1nu	
	complex(kind=8),dimension(4,4) :: mat1	
	integer,dimension(3) :: qloop, p_ext
	complex(kind=8),dimension(3) :: l_loop
	integer :: mu, nu, i
	integer :: n1, n2, n3
	integer :: s0,s0p,t0,t0p,u0,u0p
	integer :: t0_a,t0_b,t0p_a,t0p_b,u0_a,u0_b,u0p_a,u0p_b
	real(kind=8) :: multi
	complex(kind=8) :: trace
	integer :: f, fp,k 
	!Variables for the openMP paralelization
	integer :: chunk
	chunk = 1


	f = f1f2(f1,f2)
	fp = f1f2(f2,f1)

	!initialization
	l1_5b = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	SQSQS(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	p_ext(:) = 0
	trace = dcmplx(0.0d0,0.0d0)
	l_loop(:) = dcmplx(0.0d0,0.0d0)


	!part outside the loop
	call Fermion_Propagator(l,p_ext,S_prop_0)
	do k=1,3
	mat1 = matmul(matmul(Qkp(:,:,k,fp),S_prop_0(:,:,l-1,1,f1)),Qk(:,:,k,f))
	do t0=0,l
	do s0=0,l
	  SQSQS(:,:,t0,s0,k) = matmul(matmul(S_prop_0(:,:,t0,l-1,f2),mat1),S_prop_0(:,:,1,s0,f2))
	enddo
	enddo
	enddo

	!Construct  loop
!$OMP PARALLEL DO &
!$OMP DEFAULT(PRIVATE) SHARED(l,p_ext,chunk,f1,f2,SQSQS,multiplicity,mom_deg) &
!$OMP REDUCTION(+:l_loop) &
!$OMP SCHEDULE(DYNAMIC,chunk) &
!$OMP COLLAPSE(1)

	do n1=0,l-1			
	do n2=n_limit(mom_deg,n1),l-1
	do n3=n_limit(mom_deg,n2),l-1


	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3

	multi = 1.0d0*multiplicity(n1,n2,n3)
  	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	call Gluon_propagator(l,qloop,D_prop)	
	call Fermion_Propagator(l,p_ext+qloop,S_prop)


	trace = dcmplx(0.0d0,0.0d0)

	do s0=0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do s0p=0,l
	call loop_limits(l,s0p,t0p_a,t0p_b,u0p_a,u0p_b)
	do t0=t0_a,t0_b
	do t0p=t0p_a,t0p_b		
	do u0=u0_a,u0_b		
	do u0p=u0p_a,u0p_b

	do mu=0,3	
	do nu=mu_low(mu),mu_up(mu)

	call V1(l,p_ext,-p_ext-qloop,qloop,s0,t0,u0,mu,V1mu)
	call V1(l,p_ext+qloop,-p_ext,-qloop,s0p,t0p,u0p,nu,V1nu)

	do k=1,3
	trace = Trace4(SQSQS(:,:,t0p,s0,k),V1mu(:,:),S_prop(:,:,t0,s0p,f2),V1nu(:,:))
	l_loop(k) = l_loop(k) + D_prop(nu,mu,u0p,u0)*trace*multi
	enddo

	enddo	!u0
	enddo	!t0
	enddo	!s0
	enddo 
	enddo
	enddo

	enddo	!nu
	enddo	!mu


	enddo	!q3
	enddo	!q2
	enddo	!q1
!$OMP END PARALLEL DO


	l1_5b = -C2_R*dimrep*(l_loop(1)+l_loop(2)+l_loop(3))/(6.0d0*l**3)


	END SUBROUTINE l15b
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE l16a(l,f1,f2,l1_6a)
	!
	! Subroutine to compute diagram l1_6a 
	!
	implicit none

	integer,intent(in) :: l, f1, f2
	complex(kind=8),intent(out) :: l1_6a

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop_0
	complex(kind=8),dimension(4,4,0:l,0:l,3) :: SQSQS
	complex(kind=8),dimension(4,4) :: V2mu	
	complex(kind=8),dimension(4,4) :: mat1	
	integer,dimension(3) :: qloop, p_ext
	complex(kind=8),dimension(3) :: l_loop
	integer :: mu, i
	integer :: n1, n2, n3
	integer :: s0,t0,u0
	integer :: t0_a,t0_b,u0_a,u0_b
	real(kind=8) :: multi
	complex(kind=8) :: trace
	integer :: f, fp, k

	!Variables for the openMP paralelization
	integer :: chunk
	chunk = 1

	f = f1f2(f1,f2)
	fp = f1f2(f2,f1)

	!initialization
	l1_6a = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	SQSQS(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	V2mu(:,:) = dcmplx(0.0d0,0.0d0)
	p_ext(:) = 0
	trace = dcmplx(0.0d0,0.0d0)
	l_loop(:) = dcmplx(0.0d0,0.0d0)


	!part outside the loop
	call Fermion_Propagator(l,p_ext,S_prop_0)
	do k=1,3
	mat1 = matmul(matmul(Qk(:,:,k,f),S_prop_0(:,:,1,l-1,f2)),Qkp(:,:,k,fp))
	do t0=0,l
	do s0=0,l
	  SQSQS(:,:,t0,s0,k) = matmul(matmul(S_prop_0(:,:,t0,1,f1),mat1),S_prop_0(:,:,l-1,s0,f1))
	enddo
	enddo
	enddo

	!Construct  loop

!$OMP PARALLEL DO &
!$OMP DEFAULT(PRIVATE) SHARED(l,chunk,p_ext,f1,f2,SQSQS,multiplicity,mom_deg) &
!$OMP REDUCTION(+:l_loop) &
!$OMP SCHEDULE(DYNAMIC,chunk) &
!$OMP COLLAPSE(1)

	do n1=0,l-1			
	do n2=n_limit(mom_deg,n1),l-1
	do n3=n_limit(mom_deg,n2),l-1


	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3

	multi = 1.0d0*multiplicity(n1,n2,n3)
  	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	call Gluon_propagator(l,qloop,D_prop)	

	trace = dcmplx(0.0d0,0.0d0)

	do s0=0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0=t0_a,t0_b
	do u0=u0_a,u0_b	

	do mu=0,3	

	call V2(l,p_ext,-p_ext,qloop,-qloop,s0,t0,u0,u0,mu,V2mu)

	do k=1,3
	trace = Trace2(SQSQS(:,:,t0,s0,k),V2mu(:,:))
	l_loop(k) = l_loop(k) + D_prop(mu,mu,u0,u0)*trace*multi
	enddo

	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu


	enddo	!q3
	enddo	!q2
	enddo	!q1

!$OMP END PARALLEL DO

	l1_6a = C2_R*dimrep*(l_loop(1)+l_loop(2)+l_loop(3))/(6.0d0*l**3)


	END SUBROUTINE l16a




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE l16b(l,f1,f2,l1_6b)
	!
	! Subroutine to compute diagram l1_6b
	!
	implicit none

	integer,intent(in) :: l, f1, f2
	complex(kind=8),intent(out) :: l1_6b

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop_0
	complex(kind=8),dimension(4,4,0:l,0:l,3) :: SQSQS
	complex(kind=8),dimension(4,4) :: V2mu	
	complex(kind=8),dimension(4,4) :: mat1	
	integer,dimension(3) :: qloop, p_ext
	complex(kind=8),dimension(3) :: l_loop
	integer :: mu, i
	integer :: n1, n2, n3
	integer :: s0,t0,u0
	integer :: t0_a,t0_b,u0_a,u0_b
	real(kind=8) :: multi
	complex(kind=8) :: trace
	integer :: f, fp, k

	!Variables for the openMP paralelization
	integer :: chunk
	chunk = 1


	f = f1f2(f1,f2)
	fp = f1f2(f2,f1)

	!initialization
	l1_6b = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	SQSQS(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	V2mu(:,:) = dcmplx(0.0d0,0.0d0)
	p_ext(:) = 0
	trace = dcmplx(0.0d0,0.0d0)
	l_loop(:) = dcmplx(0.0d0,0.0d0)


	!part outside the loop
	call Fermion_Propagator(l,p_ext,S_prop_0)
	do k=1,3
	mat1 = matmul(matmul(Qkp(:,:,k,fp),S_prop_0(:,:,l-1,1,f1)),Qk(:,:,k,f))
	do t0=0,l
	do s0=0,l
	  SQSQS(:,:,t0,s0,k) = matmul(matmul(S_prop_0(:,:,t0,l-1,f2),mat1),S_prop_0(:,:,1,s0,f2))
	enddo
	enddo
	enddo

	!Construct  loop
!$OMP PARALLEL DO &
!$OMP DEFAULT(PRIVATE) SHARED(l,chunk,p_ext,f1,f2,SQSQS,multiplicity,mom_deg) &
!$OMP REDUCTION(+:l_loop) &
!$OMP SCHEDULE(DYNAMIC,chunk) &
!$OMP COLLAPSE(1)

	do n1=0,l-1			
	do n2=n_limit(mom_deg,n1),l-1
	do n3=n_limit(mom_deg,n2),l-1


	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3

	multi = 1.0d0*multiplicity(n1,n2,n3)
  	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	call Gluon_propagator(l,qloop,D_prop)	

	trace = dcmplx(0.0d0,0.0d0)

	do s0=0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0=t0_a,t0_b
	do u0=u0_a,u0_b	

	do mu=0,3	

	call V2(l,p_ext,-p_ext,qloop,-qloop,s0,t0,u0,u0,mu,V2mu)

	do k=1,3
	trace = Trace2(SQSQS(:,:,t0,s0,k),V2mu(:,:))
	l_loop(k) = l_loop(k) + D_prop(mu,mu,u0,u0)*trace*multi
	enddo

	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu


	enddo	!q3
	enddo	!q2
	enddo	!q1

!$OMP END PARALLEL DO


	l1_6b = C2_R*dimrep*(l_loop(1)+l_loop(2)+l_loop(3))/(6.0d0*l**3)


	END SUBROUTINE l16b




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE l17(l,f1,f2,l1_7)
	!
	! Subroutine to compute diagram l1_7
	!
	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8),intent(out) :: l1_7
	

	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,3) :: S_Qk_S,S_Qkp_S	
	complex(kind=8),dimension(4,4) :: V1mu, V1nu	
	integer :: x0,i
	integer :: s0,s0p,t0,t0p,u0,u0p
	integer :: t0_a,t0_b,t0p_a,t0p_b,u0_a,u0_b,u0p_a,u0p_b
	integer :: mu, nu
	integer :: n1, n2, n3
	integer, dimension(3) :: p,q
	complex(kind=8) :: trace
	complex(kind=8),dimension(3) :: l_temp
	real(kind=8) :: multi
	integer :: f,fp,k

	!initialization
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	l1_7 = dcmplx(0.0d0,0.0d0)
	l_temp(:) =dcmplx(0.0d0,0.0d0)
	S_Qk_S(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_Qkp_S(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	p(:) = 0
	f=f1f2(f1,f2)
	fp=f1f2(f2,f1)

	!build S_O5_S

	call Fermion_Propagator(l,p,S_prop_0)

	do k=1,3
	do s0=0,l
	do t0=0,l
	    S_Qk_S(:,:,t0,s0,k) = matmul(S_prop_0(:,:,t0,1,f1),matmul(Qk(:,:,k,f),S_prop_0(:,:,1,s0,f2)))
	    S_Qkp_S(:,:,t0,s0,k) = matmul(S_prop_0(:,:,t0,l-1,f2),matmul(Qkp(:,:,k,fp),S_prop_0(:,:,l-1,s0,f1)))
	enddo	
	enddo
	enddo


	q(1)=0
	q(2)=0
	q(3)=0

  	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)

  	call Gluon_propagator(l,q,D_prop)	

	trace = dcmplx(0.0d0,0.0d0)
	
	do s0=0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do s0p=0,l
	call loop_limits(l,s0p,t0p_a,t0p_b,u0p_a,u0p_b)
	do t0=t0_a,t0_b
	do t0p=t0p_a,t0p_b		!Due to the delta(x0+1,y0),delta(x0,y0),delta(x0-1,y0) in the vertex functions 
	do u0=u0_a,u0_b		!Due to the delta(x0+1,y0),delta(x0,y0),delta(x0-1,y0) in the vertex functions 
	do u0p=u0p_a,u0p_b

	  do mu=0,3	
	  do nu=mu_low(mu),mu_up(mu)

	  call V1(l,p,-p-q,q,s0,t0,u0,mu,V1mu)
	  call V1(l,p+q,-p,-q,s0p,t0p,u0p,nu,V1nu)

	  do k=1,3
	  trace = Trace4(S_Qk_S(:,:,t0p,s0,k),V1mu(:,:),S_Qkp_S(:,:,t0,s0p,k),V1nu(:,:))
	  l_temp(k) = l_temp(k) + trace*D_prop(nu,mu,u0p,u0)
	  enddo

	  enddo	!nu
	  enddo !mu
	  enddo	!u0p
	  enddo	!u0
	enddo	!t0p
	enddo	!t0
	enddo	!s0p
	enddo	!s0


 	l1_7 = -dimrep*C2_R*(l_temp(1)+l_temp(2)+l_temp(3))/((l**3)*6.0d0)

	
	END SUBROUTINE l17


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE l18a(l,f1,f2,l1_8a)
	!
	! Subroutine to compute diagram l1_8a 
	!
	implicit none

	integer,intent(in) :: l, f1, f2
	complex(kind=8),intent(out) :: l1_8a



	integer :: n1,n2,n3
	integer,dimension(3) :: qloop
	complex(kind=8) :: loop
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop!,K_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop
	complex(kind=8) :: l1_0
	complex(kind=8),dimension(3) :: l_temp
	real(kind=8) :: multi
	integer,dimension(3) :: nvec
	integer :: f,fp,k

	f = f1f2(f1,f2)
	fp = f1f2(f2,f1)
	nvec(:) = 0

	!initialization
	l1_8a = dcmplx(0.0d0,0.0d0)
	loop = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	l1_0 = dcmplx(0.0d0,0.0d0)
	l_temp(:) = dcmplx(0.0d0,0.0d0)

	! Construct tree-level part:
	call Fermion_Propagator(l,nvec,S_Prop)
	do k = 1,3
	l_temp(k) = 0.5*dimrep*Trace4(Qk(:,:,k,f),S_prop(:,:,1,l-1,f2),Qkp(:,:,k,fp),S_prop(:,:,l-1,1,f1))
	enddo
	l1_0 = (l_temp(1)+l_temp(2)+l_temp(3))/3.0d0

	

	!Construct the loop
	do n1=0,l-1			
	do n2=n_limit(mom_deg,n1),l-1
	do n3=n_limit(mom_deg,n2),l-1


	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3

	   multi = 1.0d0*multiplicity(n1,n2,n3)

  	   D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	   call Gluon_propagator(l,qloop,D_prop)	
	   loop = loop + D_prop(0,0,l-1,l-1)*multi

	enddo
	enddo
	enddo

	l1_8a = -0.5d0*C2_R*loop*l1_0/((l*1.0d0)**3)

	END SUBROUTINE l18a
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE l19(l,f1,f2,l1_9)
	!
	! Subroutine to compute diagram l1_9 
	!
	implicit none

	integer,intent(in) :: l, f1, f2
	complex(kind=8),intent(out) :: l1_9


	integer,dimension(3) :: qloop,nvec
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop!,K_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop
	complex(kind=8) :: l1_0
	complex(kind=8),dimension(3) :: l_temp
	integer :: f, fp,k

	nvec(:) = 0

	f = f1f2(f1,f2)
	fp = f1f2(f2,f1)

	!initialization
	l1_9 = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	l1_0 = dcmplx(0.0d0,0.0d0)
	l_temp(:) = dcmplx(0.0d0,0.0d0)

	! Construct tree-level part:
	call Fermion_Propagator(l,nvec,S_Prop)
	do k = 1,3
	l_temp(k) = 0.5*dimrep*Trace4(Qk(:,:,k,f),S_prop(:,:,1,l-1,f2),Qkp(:,:,k,fp),S_prop(:,:,l-1,1,f1))
	enddo
	l1_0 = (l_temp(1)+l_temp(2)+l_temp(3))/3.0d0

	!Construct the loop
	qloop(1)=0
	qloop(2)=0
	qloop(3)=0

  	call Gluon_propagator(l,qloop,D_prop)	

	l1_9 =  C2_R*l1_0*D_prop(0,0,l-1,l-1)/((l*1.0d0)**3)

	END SUBROUTINE l19

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE l110a(l,f1,f2,l1_10a)
	!
	! Subroutine to compute diagram l1_10a 
	!
	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8),intent(out) :: l1_10a
	
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, s_prop_0
	complex(kind=8),dimension(4,4) :: Vmu1	
	integer,dimension(3) :: qloop, p_ext
	complex(kind=8),dimension(4,4,3) :: mat1!, mat2, mat3	!s,s,f1f2
	real(kind=8) :: multi
	integer :: mu,u0,s0,t0,n1,n2,n3
	integer :: t0_a, t0_b, u0_a, u0_b
	complex(kind=8),dimension(3) :: Trace
	integer :: f, fp,k

	f = f1f2(f1,f2) 
	fp = f1f2(f2,f1) 


	!initialization
	p_ext(:) = 0
	qloop(:)=0
	l1_10a = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	p_ext(:) = 0
	mat1(:,:,:) = dcmplx(0.0d0,0.0d0)
	Vmu1(:,:) = dcmplx(0.0d0,0.0d0)
	Trace(:) = dcmplx(0.0d0,0.0d0)

	call Fermion_Propagator(l,p_ext,S_prop_0)
	do k=1,3
	mat1(:,:,k) = matmul(Qk(:,:,k,f),matmul(S_prop_0(:,:,1,l-1,f2),Qkp(:,:,k,fp)))
	enddo

	do n1=0,l-1			
	do n2=n_limit(mom_deg,n1),l-1
	do n3=n_limit(mom_deg,n2),l-1


	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3

	   multi = 1.0d0*multiplicity(n1,n2,n3)

  	   D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	   S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

  	   call Gluon_propagator(l,qloop,D_prop)	
	   call Fermion_Propagator(l,qloop,S_prop)


	do mu=mu_low(0),mu_up(0)

	do s0=0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0=t0_a,t0_b
	do u0=u0_a,u0_b

	call V1(l,p_ext+qloop,-p_ext,-qloop,s0,t0,u0,mu,Vmu1)

	do k=1,3
	Trace(k) = Trace(k) + Trace4(mat1(:,:,k),S_prop(:,:,l-1,s0,f1),Vmu1,S_prop_0(:,:,t0,1,f1))*D_prop(mu,0,u0,l-1)*multi
	enddo


	enddo	!u0
	enddo	!t0
	enddo	!s0
	
	
	enddo	!mu


	enddo	!momenta
	enddo
	enddo	

	l1_10a = -C2_R*dimrep*(Trace(1)+Trace(2)+Trace(3))/(6.0d0*l**3)

	END SUBROUTINE l110a

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE l110b(l,f1,f2,l1_10b)
	!
	! Subroutine to compute diagram l1_1b
	!
	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8),intent(out) :: l1_10b
	
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, s_prop_0
	complex(kind=8),dimension(4,4) :: Vmu1	
	integer,dimension(3) :: qloop, p_ext
	complex(kind=8),dimension(4,4,3) :: mat1!, mat2, mat3	!s,s,f1f2
	real(kind=8) :: multi
	integer :: mu,u0,s0,t0,n1,n2,n3
	integer :: t0_a, t0_b, u0_a, u0_b
	complex(kind=8),dimension(3) :: Trace
	integer :: f, fp,k

	f = f1f2(f1,f2) 
	fp = f1f2(f2,f1) 


	!initialization
	p_ext(:) = 0
	qloop(:)=0
	l1_10b = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	p_ext(:) = 0
	mat1(:,:,:) = dcmplx(0.0d0,0.0d0)
	Vmu1(:,:) = dcmplx(0.0d0,0.0d0)
	Trace(:) = dcmplx(0.0d0,0.0d0)

	call Fermion_Propagator(l,p_ext,S_prop_0)
	do k=1,3
	mat1(:,:,k) = matmul(Qkp(:,:,k,fp),matmul(S_prop_0(:,:,l-1,1,f1),Qk(:,:,k,f)))
	enddo

	do n1=0,l-1			
	do n2=n_limit(mom_deg,n1),l-1
	do n3=n_limit(mom_deg,n2),l-1


	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3

	   multi = 1.0d0*multiplicity(n1,n2,n3)

  	   D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	   S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

  	   call Gluon_propagator(l,qloop,D_prop)	
	   call Fermion_Propagator(l,qloop,S_prop)


	do mu=mu_low(0),mu_up(0)

	do s0=0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0=t0_a,t0_b
	do u0=u0_a,u0_b

	call V1(l,p_ext,-p_ext-qloop,qloop,s0,t0,u0,mu,Vmu1)

	do k=1,3
	Trace(k) = Trace(k) + Trace4(mat1(:,:,k),S_prop_0(:,:,1,s0,f2),Vmu1,S_prop(:,:,t0,l-1,f2))*D_prop(0,mu,l-1,u0)*multi
	enddo


	enddo	!u0
	enddo	!t0
	enddo	!s0
	
	
	enddo	!mu


	enddo	!momenta
	enddo
	enddo	

	l1_10b = C2_R*dimrep*(Trace(1)+Trace(2)+Trace(3))/(6.0d0*l**3)

	END SUBROUTINE l110b

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE l111a(l,f1,f2,l1_11a)
	!
	! Subroutine to compute diagram l1_11a 
	!
	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8),intent(out) :: l1_11a
	
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: s_prop_0
	complex(kind=8),dimension(4,4) :: Vmu1	
	integer,dimension(3) :: qloop, p_ext
	complex(kind=8),dimension(4,4,3) :: mat1
	real(kind=8) :: multi
	integer :: mu,u0,s0,t0,n1,n2,n3
	integer :: t0_a, t0_b, u0_a, u0_b
	complex(kind=8),dimension(3) :: Trace
	integer :: f, fp, k

	f = f1f2(f1,f2) 
	fp = f1f2(f2,f1) 


	!initialization
	p_ext(:) = 0
	qloop(:)=0
	l1_11a = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	p_ext(:) = 0
	mat1(:,:,:) = dcmplx(0.0d0,0.0d0)
	Vmu1(:,:) = dcmplx(0.0d0,0.0d0)
	Trace(:) = dcmplx(0.0d0,0.0d0)

	call Fermion_Propagator(l,p_ext,S_prop_0)
	do k=1,3
	mat1(:,:,k) = matmul(Qkp(:,:,k,fp),matmul(S_prop_0(:,:,l-1,1,f1),Qk(:,:,k,f)))
	enddo

	qloop(1)=0
	qloop(2)=0
	qloop(3)=0


  	call Gluon_propagator(l,qloop,D_prop)	

	do mu=mu_low(0),mu_up(0)

	do s0=0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0=t0_a,t0_b
	do u0=u0_a,u0_b

	call V1(l,p_ext,-p_ext-qloop,qloop,s0,t0,u0,mu,Vmu1)

	do k=1,3
	Trace(k) = Trace(k) + Trace4(S_prop_0(:,:,1,s0,f2),Vmu1(:,:),S_prop_0(:,:,t0,l-1,f2),mat1(:,:,k))*D_prop(0,mu,l-1,u0)
	enddo

	enddo	!u0
	enddo	!t0
	enddo	!s0
	
	
	enddo	!mu


	l1_11a = -C2_R*dimrep*(Trace(1)+Trace(2)+Trace(3))/(6.0d0*l**3)

	END SUBROUTINE l111a

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE l111b(l,f1,f2,l1_11b)
	!
	! Subroutine to compute diagram l1_11b 
	!
	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8),intent(out) :: l1_11b
	
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: s_prop_0
	complex(kind=8),dimension(4,4) :: Vmu1	
	integer,dimension(3) :: qloop, p_ext
	complex(kind=8),dimension(4,4,3) :: mat1
	real(kind=8) :: multi
	integer :: mu,u0,s0,t0,n1,n2,n3
	integer :: t0_a, t0_b, u0_a, u0_b
	complex(kind=8),dimension(3) :: Trace
	integer :: f, fp, k

	f = f1f2(f1,f2) 
	fp = f1f2(f2,f1) 


	!initialization
	p_ext(:) = 0
	qloop(:)=0
	l1_11b = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	p_ext(:) = 0
	mat1(:,:,:) = dcmplx(0.0d0,0.0d0)
	Vmu1(:,:) = dcmplx(0.0d0,0.0d0)
	Trace(:) = dcmplx(0.0d0,0.0d0)

	call Fermion_Propagator(l,p_ext,S_prop_0)
	do k=1,3
	mat1(:,:,k) = matmul(Qk(:,:,k,f),matmul(S_prop_0(:,:,1,l-1,f2),Qkp(:,:,k,fp)))
	enddo

	qloop(1)=0
	qloop(2)=0
	qloop(3)=0


  	call Gluon_propagator(l,qloop,D_prop)	

	do mu=mu_low(0),mu_up(0)

	do s0=0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0=t0_a,t0_b
	do u0=u0_a,u0_b

	call V1(l,p_ext,-p_ext-qloop,qloop,s0,t0,u0,mu,Vmu1)

	do k=1,3
	Trace(k) = Trace(k) + Trace4(S_prop_0(:,:,l-1,s0,f1),Vmu1(:,:),S_prop_0(:,:,t0,1,f1),mat1(:,:,k))*D_prop(mu,0,u0,l-1)
	enddo

	enddo	!u0
	enddo	!t0
	enddo	!s0
	
	
	enddo	!mu


	l1_11b = C2_R*dimrep*(Trace(1)+Trace(2)+Trace(3))/(6.0d0*l**3)

	END SUBROUTINE l111b

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE l112a(l,f1,f2,l1_12a)
	!
	! Subroutine to compute diagram l1_12a 
	!
	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8),intent(out) :: l1_12a
	
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	integer,dimension(3) :: qloop, p_ext
	complex(kind=8),dimension(4,4,3) :: mat1
	real(kind=8) :: multi
	integer :: n1,n2,n3
	complex(kind=8),dimension(3) :: Trace
	integer :: f, fp, k

	f =f1f2(f1,f2)
	fp =f1f2(f2,f1)

	!initialization
	l1_12a = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	mat1(:,:,:) = dcmplx(0.0d0,0.0d0)
	Trace(:) = dcmplx(0.0d0,0.0d0)
	p_ext(:) =0


	call Fermion_Propagator(l,p_ext,S_prop_0)
	do k=1,3
	mat1(:,:,k) = matmul(Qk(:,:,k,f),matmul(S_prop_0(:,:,1,l-1,f2),Qkp(:,:,k,fp)))
	enddo

	do n1=0,l-1			
	do n2=n_limit(mom_deg,n1),l-1
	do n3=n_limit(mom_deg,n2),l-1


	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3

	multi = 1.0d0*multiplicity(n1,n2,n3)

  	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

  	call Gluon_propagator(l,qloop,D_prop)	
	call Fermion_Propagator(l,qloop,S_prop)

	do k=1,3
	Trace(k) = Trace(k) + Trace2(mat1(:,:,k),S_prop(:,:,l-1,1,f1))*D_prop(0,0,0,l-1)*multi
	enddo

	enddo	!enddo
	enddo
	enddo	

	l1_12a = -C2_R*dimrep*(Trace(1)+Trace(2)+Trace(3))/(6.0d0*l**3)

	END SUBROUTINE l112a


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE l112b(l,f1,f2,l1_12b)
	!
	! Subroutine to compute diagram l1_12b
	!
	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8),intent(out) :: l1_12b
	
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	integer,dimension(3) :: qloop, p_ext
	complex(kind=8),dimension(4,4,3) :: mat1
	real(kind=8) :: multi
	integer :: n1,n2,n3
	complex(kind=8),dimension(3) :: Trace
	integer :: f, fp, k

	f =f1f2(f1,f2)
	fp =f1f2(f2,f1)

	!initialization
	l1_12b = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	mat1(:,:,:) = dcmplx(0.0d0,0.0d0)
	Trace(:) = dcmplx(0.0d0,0.0d0)
	p_ext(:) =0


	call Fermion_Propagator(l,p_ext,S_prop_0)
	do k=1,3
	mat1(:,:,k) = matmul(Qkp(:,:,k,fp),matmul(S_prop_0(:,:,l-1,1,f1),Qk(:,:,k,f)))
	enddo

	do n1=0,l-1			
	do n2=n_limit(mom_deg,n1),l-1
	do n3=n_limit(mom_deg,n2),l-1


	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3

	multi = 1.0d0*multiplicity(n1,n2,n3)

  	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

  	call Gluon_propagator(l,qloop,D_prop)	
	call Fermion_Propagator(l,qloop,S_prop)

	do k=1,3
	Trace(k) = Trace(k) + Trace2(mat1(:,:,k),S_prop(:,:,1,l-1,f2))*D_prop(0,0,l-1,0)*multi
	enddo

	enddo	!enddo
	enddo
	enddo	

	l1_12b = -C2_R*dimrep*(Trace(1)+Trace(2)+Trace(3))/(6.0d0*l**3)

	END SUBROUTINE l112b


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE l113a(l,f1,f2,l1_13a)
	!
	! Subroutine to compute diagram l1_13a 
	!
	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8),intent(out) :: l1_13a

	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	integer,dimension(3) :: qloop, p_ext
	complex(kind=8),dimension(3) :: l_temp
	integer :: f, fp, k


	f =f1f2(f1,f2)
	fp =f1f2(f2,f1)

	!initialization
	l1_13a = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	l_temp(:) = dcmplx(0.0d0,0.0d0)
	qloop(:) = 0
	p_ext(:) = 0


	! Construct tree-level part:
	call Fermion_Propagator(l,p_ext,S_Prop)
	do k = 1,3
	l_temp(k) = Trace4(Qk(:,:,k,f),S_prop(:,:,1,l-1,f2),Qkp(:,:,k,fp),S_prop(:,:,l-1,1,f1))
	enddo
  	call Gluon_propagator(l,qloop,D_prop)	

	l1_13a = C2_R*dimrep*(l_temp(1)+l_temp(2)+l_temp(3))*D_prop(0,0,l-1,0)/(6.0d0*l**3)


	END SUBROUTINE l113a


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE l1dm(l,f1,f2,l1_dm)
	!
	! Subroutine to compute diagram l1_dm 
	!
	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8),intent(out) :: l1_dm
	
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop
	integer, dimension(3) :: p_ext
	complex(kind=8),dimension(3) :: l_temp
	integer :: f,fp,s0,k

	f=f1f2(f1,f2)
	fp=f1f2(f2,f1)
	p_ext(:) = 0

	!initialization
	l1_dm = dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	l_temp(:) = dcmplx(0.0d0,0.0d0)

	call Fermion_Propagator(l,p_ext,S_prop)

	do s0 = 0,l
	do k=1,3
	l_temp(k) = l_temp(k)-dimrep*Trace5(Qk(:,:,k,f),S_prop(:,:,1,l-1,f2),Qkp(:,:,k,fp),S_prop(:,:,l-1,s0,f1),S_prop(:,:,s0,1,f1))
	l_temp(k) = l_temp(k)-dimrep*Trace5(Qk(:,:,k,f),S_prop(:,:,1,s0,f2),S_prop(:,:,s0,l-1,f2),Qkp(:,:,k,fp),S_prop(:,:,l-1,1,f1))
	enddo
	enddo

	l1_dm = (l_temp(1)+l_temp(2)+l_temp(3))/6.0d0


	END SUBROUTINE l1dm


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE l1zf(l,f1,f2,l1_zf)
	!
	! Subroutine to compute diagram l1_zf
	!
	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8),intent(out) :: l1_zf
	
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop
	complex(kind=8) :: l1_zf_1, l1_zf_2
	integer, dimension(3) :: p_ext
	complex(kind=8),dimension(3) :: l_temp1, l_temp2
	integer :: f,fp,k

	f=f1f2(f1,f2)
	fp=f1f2(f2,f1)
	p_ext(:) = 0

	!initialization
	l1_zf = dcmplx(0.0d0,0.0d0)
	l1_zf_1 = dcmplx(0.0d0,0.0d0)
	l1_zf_2 = dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	l_temp1(:) = dcmplx(0.0d0,0.0d0)
	l_temp2(:) = dcmplx(0.0d0,0.0d0)

	call Fermion_Propagator(l,p_ext,S_prop)
	
	do k=1,3

	l_temp1(k) = Trace5(Qk(:,:,k,f),S_prop(:,:,1,l-1,f2),Qkp(:,:,k,fp),S_prop(:,:,l-1,0,f1),S_prop(:,:,0,1,f1))
	l_temp1(k) = l_temp1(k) + Trace5(Qk(:,:,k,f),S_prop(:,:,1,l-1,f2),Qkp(:,:,k,fp),S_prop(:,:,l-1,l,f1),S_prop(:,:,l,1,f1))

	l_temp2(k) = Trace5(Qk(:,:,k,f),S_prop(:,:,1,0,f2),S_prop(:,:,0,l-1,f2),Qkp(:,:,k,fp),S_prop(:,:,l-1,1,f1))
	l_temp2(k) = l_temp2(k) + Trace5(Qk(:,:,k,f),S_prop(:,:,1,l,f2),S_prop(:,:,l,l-1,f2),Qkp(:,:,k,fp),S_prop(:,:,l-1,1,f1))

	enddo	

	l1_zf_1 = (l_temp1(1)+l_temp1(2)+l_temp1(3))/3.0d0
	l1_zf_2 = (l_temp2(1)+l_temp2(2)+l_temp2(3))/3.0d0

	l1_zf = -dimrep*(l1_zf_1+l1_zf_2)/2.0d0

	END SUBROUTINE l1zf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE l1ds(l,f1,f2,l1_ds)
	!
	! Subroutine to compute diagram l1_ds
	!
	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8),intent(out) :: l1_ds
	
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop
	complex(kind=8) :: l1_ds_1, l1_ds_2
	integer, dimension(3) :: p_ext
	real(kind=8),dimension(3) ::  pe_plus, ptwit, phat
	complex(kind=8),dimension(4,4) :: Amat
	complex(kind=8),dimension(3) :: l_temp1, l_temp2
	integer :: f,fp,k

	f=f1f2(f1,f2)
	fp=f1f2(f2,f1)
	p_ext(:) = 0

	!initialization
	l1_ds = dcmplx(0.0d0,0.0d0)
	l1_ds_1 = dcmplx(0.0d0,0.0d0)
	l1_ds_2 = dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	l_temp1(:) = dcmplx(0.0d0,0.0d0)
	l_temp2(:) = dcmplx(0.0d0,0.0d0)


	pe_plus(:) = p_ext(:)*(2.0d0*pi)/(l*1.0d0) + theta/l
	ptwit(:) = sin(pe_plus(:))
	phat(:) = 2.0d0*sin(pe_plus(:)/2.0d0)

	Amat = imag*(ptwit(1)*g1+ptwit(2)*g2+ptwit(3)*g3)
	Amat = Amat + 0.5d0*II*( phat(1)**2 + phat(2)**2 + phat(3)**2 )



	call Fermion_Propagator(l,p_ext,S_prop)

	do k=1,3

	l_temp1(k) = Trace6(Qk(:,:,k,f),S_prop(:,:,1,l-1,f2),Qkp(:,:,k,fp),S_prop(:,:,l-1,0,f1),Amat,S_prop(:,:,0,1,f1))
	l_temp1(k) = l_temp1(k) + Trace6(Qk(:,:,k,f),S_prop(:,:,1,l-1,f2),Qkp(:,:,k,fp),S_prop(:,:,l-1,l,f1),Amat,S_prop(:,:,l,1,f1))

	l_temp2(k) = Trace6(Qk(:,:,k,f),S_prop(:,:,1,0,f2),Amat,S_prop(:,:,0,l-1,f2),Qkp(:,:,k,fp),S_prop(:,:,l-1,1,f1))
	l_temp2(k) = l_temp2(k) + Trace6(Qk(:,:,k,f),S_prop(:,:,1,l,f2),Amat,S_prop(:,:,l,l-1,f2),Qkp(:,:,k,fp),S_prop(:,:,l-1,1,f1))


	enddo	


	l1_ds_1 = (l_temp1(1)+l_temp1(2)+l_temp1(3))/3.0d0
	l1_ds_2 = (l_temp2(1)+l_temp2(2)+l_temp2(3))/3.0d0

	l1_ds = -dimrep*(l1_ds_1+l1_ds_2)/2.0d0

	END SUBROUTINE l1ds


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






	REAL(kind=8) FUNCTION N_l1(f1,f2)
	integer :: f1, f2
	if(f1.eq.f2)then
	  N_l1 = -1.0d0
	else
	  N_l1 = 1.0d0
	endif
	END FUNCTION N_l1


	REAL(kind=8) FUNCTION M_l1(f1)
	integer :: f1
	if(f1.eq.1)then
	  M_l1 = -1.0d0
	else
	  M_l1 = 1.0d0
	endif
	END FUNCTION M_l1


	INTEGER FUNCTION gY1_l1(f1,f2)
	integer :: f1,f2
	if(f1.eq.f2)then
	  gY1_l1 = 2
	else
	  gY1_l1 = 3
	endif
	END FUNCTION gY1_l1

	INTEGER FUNCTION gY2_l1(f1,f2)
	integer :: f1,f2
	if(f1.eq.f2)then
	  gY2_l1 = 1
	else
	  gY2_l1 = 4
	endif
	END FUNCTION gY2_l1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	END MODULE OneLoop_l1
