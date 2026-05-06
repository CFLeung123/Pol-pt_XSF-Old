

	MODULE OneLoop_g1
	!
	! Contains the set of routines for the calculaton of the 1-loop
	! diagrams contributing to the g1 correlation functions (including
	! the contribution of counterterms.
	!
	USE GammaMatrices
	USE Input
	USE Parameters
	USE Vertices
	USE Propagators
	USE TreeLevel_g1
	USE omp_lib

	implicit none

	CONTAINS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE g11a(l,f1,f2,g1_1a)
	!
	! Subroutine to compute diagram g1_1a 
	!
	implicit none

	integer,intent(in) :: l, f1, f2
	complex(kind=8),intent(out) :: g1_1a



	integer :: n1,n2,n3
	integer,dimension(3) :: qloop
	complex(kind=8) :: loop
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop!,K_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop
	complex(kind=8) :: g1_0
	real(kind=8) :: multi
	integer,dimension(3) :: nvec
	integer :: f,fp

	f = f1f2(f1,f2)
	fp = f1f2(f2,f1)
	nvec(:) = 0

	!initialization
	g1_1a = dcmplx(0.0d0,0.0d0)
	loop = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	g1_0 = dcmplx(0.0d0,0.0d0)



	! Construct tree-level part:
	call Fermion_Propagator(l,nvec,S_Prop)
	g1_0 = 0.5d0*dimrep*Trace4(Q5(:,:,f),S_prop(:,:,1,l-1,f2),Q5p(:,:,fp),S_prop(:,:,l-1,1,f1))	

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

	g1_1a = -0.5d0*C2_R*loop*g1_0/((l*1.0d0)**3)

	END SUBROUTINE g11a

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE g12(l,f1,f2,g1_2)
	!
	! Subroutine to compute diagram g1_2a 
	!
	implicit none

	integer,intent(in) :: l, f1, f2
	complex(kind=8),intent(out) :: g1_2

	integer,dimension(3) :: qloop
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop!,K_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop
	complex(kind=8) :: g1_0
	integer,dimension(3) :: nvec
	integer :: f, fp

	f = f1f2(f1,f2)
	fp = f1f2(f2,f1)
	nvec(:) = 0	

	!initialization
	g1_2 = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	g1_0 = dcmplx(0.0d0,0.0d0)

	! Construct tree-level part:
	call Fermion_Propagator(l,nvec,S_Prop)
	g1_0 = 0.5d0*dimrep*Trace4(Q5(:,:,f),S_prop(:,:,1,l-1,f2),Q5p(:,:,fp),S_prop(:,:,l-1,1,f1))	

	!Construct the loop
	qloop(1)=0
	qloop(2)=0
	qloop(3)=0

  	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	call Gluon_propagator(l,qloop,D_prop)	

	g1_2 = C2_R*g1_0*D_prop(0,0,0,0)/((l*1.0d0)**3)

	END SUBROUTINE g12

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE g13a(l,f1,f2,g1_3a)
	!
	! Subroutine to compute diagram g1_3a 
	!
	implicit none

	integer,intent(in) :: l, f1, f2
	complex(kind=8),intent(out) :: g1_3a



	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,0:l) :: Q5_S_Q5p_S
	complex(kind=8),dimension(4,4) :: V1mu	
	complex(kind=8),dimension(4,4) :: mat1	
	integer,dimension(3) :: qloop, p_ext
	complex(kind=8) :: g_loop
	integer ::  t0, u0, s0, mu, i
	integer :: n1, n2, n3
	integer :: t0_a, t0_b, u0_a, u0_b
	real(kind=8) :: multi
	complex(kind=8) :: trace
	integer :: f, fp

	!Variables for the openMP paralelization
	integer :: chunk
	chunk = 1

	f = f1f2(f1,f2)
	fp = f1f2(f2,f1)

	!initialization
	g1_3a = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	Q5_S_Q5p_S(:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	p_ext(:) = 0
	trace = dcmplx(0.0d0,0.0d0)
	g_loop = dcmplx(0.0d0,0.0d0)


	!part outside the loop
	call Fermion_Propagator(l,p_ext,S_prop_0)
	mat1 = matmul(matmul(Q5(:,:,f),S_prop_0(:,:,1,l-1,f2)),Q5p(:,:,fp))
	do s0=0,l
	  Q5_S_Q5p_S(:,:,s0) = matmul(mat1,S_prop_0(:,:,l-1,s0,f1))
	enddo

	!Construct  loop
!$OMP PARALLEL DO &
!$OMP DEFAULT(PRIVATE) SHARED(l,p_ext,chunk,f1,f2,Q5_S_Q5p_S,multiplicity,mom_deg) &
!$OMP REDUCTION(+:g_loop) &
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

	do mu=mu_low(0),mu_up(0)

	do s0=0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0=t0_a,t0_b
	do u0=u0_a,u0_b

	call V1(l,p_ext,-p_ext-qloop,qloop,s0,t0,u0,mu,V1mu)

	trace = Trace3(Q5_S_Q5p_S(:,:,s0),V1mu,S_prop(:,:,t0,1,f1))
	g_loop = g_loop + D_prop(0,mu,0,u0)*trace*multi
	

	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu


	enddo	!q3
	enddo	!q2
	enddo	!q1
!$OMP END PARALLEL DO

	g1_3a = -0.5*C2_R*dimrep*g_loop/(1.0d0*l**3)


	END SUBROUTINE g13a




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE g13b(l,f1,f2,g1_3b)
	!
	! Subroutine to compute diagram g1_3b
	!
	implicit none

	integer,intent(in) :: l, f1, f2
	complex(kind=8),intent(out) :: g1_3b



	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,0:l) :: S_Q5p_S_Q5
	complex(kind=8),dimension(4,4) :: V1mu	
	complex(kind=8),dimension(4,4) :: mat1	
	integer,dimension(3) :: qloop, p_ext
	complex(kind=8) :: g_loop
	integer ::  t0, u0, s0, mu, i
	integer :: n1, n2, n3
	integer :: t0_a, t0_b, u0_a, u0_b
	real(kind=8) :: multi
	complex(kind=8) :: trace
	integer :: f, fp

	f = f1f2(f1,f2)
	fp = f1f2(f2,f1)

	!initialization
	g1_3b = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_Q5p_S_Q5(:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	p_ext(:) = 0
	trace = dcmplx(0.0d0,0.0d0)
	g_loop = dcmplx(0.0d0,0.0d0)


	!part outside the loop
	call Fermion_Propagator(l,p_ext,S_prop_0)
	mat1 = matmul(matmul(Q5p(:,:,fp),S_prop_0(:,:,l-1,1,f1)),Q5(:,:,f))
	do t0=0,l
	  S_Q5p_S_Q5(:,:,t0) = matmul(S_prop_0(:,:,t0,l-1,f2),mat1)
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

	trace = Trace3(S_Q5p_S_Q5(:,:,t0),S_prop(:,:,1,s0,f2),V1mu)
	g_loop = g_loop + D_prop(mu,0,u0,0)*trace*multi
	

	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu


	enddo	!q3
	enddo	!q2
	enddo	!q1

	g1_3b = 0.5*C2_R*dimrep*g_loop/(1.0d0*l**3)


	END SUBROUTINE g13b
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE g14a(l,f1,f2,g1_4a)
	!
	! Subroutine to compute diagram g1_4a 
	!
	implicit none

	integer,intent(in) :: l, f1, f2
	complex(kind=8),intent(out) :: g1_4a



	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop_0
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop, p_ext
	complex(kind=8) :: g_loop
	integer ::  t0, u0, s0, mu, i
	integer :: n1, n2, n3
	integer :: t0_a, t0_b, u0_a, u0_b
	complex(kind=8) :: trace
	integer :: f, fp



	f = f1f2(f1,f2)
	fp = f1f2(f2,f1)

	!initialization
	g1_4a = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	p_ext(:) = 0
	qloop(:) = 0
	trace = dcmplx(0.0d0,0.0d0)
	g_loop = dcmplx(0.0d0,0.0d0)


	!part outside the loop
	call Fermion_Propagator(l,p_ext,S_prop_0)
	call Gluon_propagator(l,qloop,D_prop)	

	do mu=mu_low(0),mu_up(0)

	do s0=0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0=t0_a,t0_b
	do u0=u0_a,u0_b

	call V1(l,p_ext,-p_ext-qloop,qloop,s0,t0,u0,mu,V1mu)

	trace = Trace6(Q5(:,:,f),S_prop_0(:,:,1,s0,f2),V1mu,S_prop_0(:,:,t0,l-1,f2),Q5p(:,:,fp),S_prop_0(:,:,l-1,1,f1))
	g_loop = g_loop + D_prop(0,mu,0,u0)*trace
	

	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu

	g1_4a = -0.5*C2_R*dimrep*g_loop/(1.0d0*l**3)


	END SUBROUTINE g14a




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE g14b(l,f1,f2,g1_4b)
	!
	! Subroutine to compute diagram g1_4b
	!
	implicit none

	integer,intent(in) :: l, f1, f2
	complex(kind=8),intent(out) :: g1_4b



	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) ::  S_prop_0
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop, p_ext
	complex(kind=8) :: g_loop
	integer ::  t0, u0, s0, mu, i
	integer :: n1, n2, n3
	integer :: t0_a, t0_b, u0_a, u0_b
	complex(kind=8) :: trace
	integer :: f, fp

	f = f1f2(f1,f2)
	fp = f1f2(f2,f1)

	!initialization
	g1_4b = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	p_ext(:) = 0
	qloop(:) = 0
	trace = dcmplx(0.0d0,0.0d0)
	g_loop = dcmplx(0.0d0,0.0d0)


	!part outside the loop
	call Fermion_Propagator(l,p_ext,S_prop_0)
  	call Gluon_propagator(l,qloop,D_prop)	
	
	do mu=mu_low(0),mu_up(0)

	do s0=0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0=t0_a,t0_b
	do u0=u0_a,u0_b

	call V1(l,p_ext+qloop,-p_ext,-qloop,s0,t0,u0,mu,V1mu)

	trace = Trace6(Q5(:,:,f),S_prop_0(:,:,1,l-1,f2),Q5p(:,:,fp),S_prop_0(:,:,l-1,s0,f1),V1mu,S_prop_0(:,:,t0,1,f1))
	g_loop = g_loop + D_prop(mu,0,u0,0)*trace
	

	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu

	g1_4b = 0.5*C2_R*dimrep*g_loop/(1.0d0*l**3)


	END SUBROUTINE g14b
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE g15a(l,f1,f2,g1_5a)
	!
	! Subroutine to compute diagram g1_5a 
	!
	implicit none

	integer,intent(in) :: l, f1, f2
	complex(kind=8),intent(out) :: g1_5a

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,0:l,0:l) :: SQSQS
	complex(kind=8),dimension(4,4) :: V1mu, V1nu	
	complex(kind=8),dimension(4,4) :: mat1	
	integer,dimension(3) :: qloop, p_ext
	complex(kind=8) :: g_loop
	integer :: mu, nu, i
	integer :: n1, n2, n3
	integer :: s0,s0p,t0,t0p,u0,u0p
	integer :: t0_a,t0_b,t0p_a,t0p_b,u0_a,u0_b,u0p_a,u0p_b
	real(kind=8) :: multi
	complex(kind=8) :: trace
	integer :: f, fp

	!Variables for the openMP paralelization
	integer :: chunk
	chunk = 1

	f = f1f2(f1,f2)
	fp = f1f2(f2,f1)

	!initialization
	g1_5a = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	SQSQS(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	p_ext(:) = 0
	trace = dcmplx(0.0d0,0.0d0)
	g_loop = dcmplx(0.0d0,0.0d0)


	!part outside the loop
	call Fermion_Propagator(l,p_ext,S_prop_0)
	mat1 = matmul(matmul(Q5(:,:,f),S_prop_0(:,:,1,l-1,f2)),Q5p(:,:,fp))

	do t0=0,l
	do s0=0,l
	  SQSQS(:,:,t0,s0) = matmul(matmul(S_prop_0(:,:,t0,1,f1),mat1),S_prop_0(:,:,l-1,s0,f1))
	enddo
	enddo

	!Construct  loop

!$OMP PARALLEL DO &
!$OMP DEFAULT(PRIVATE) SHARED(l,p_ext,chunk,f1,f2,SQSQS,multiplicity,mom_deg) &
!$OMP REDUCTION(+:g_loop) &
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

	trace = Trace4(SQSQS(:,:,t0p,s0),V1mu(:,:),S_prop(:,:,t0,s0p,f1),V1nu(:,:))
	g_loop = g_loop + D_prop(nu,mu,u0p,u0)*trace*multi
	

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

	g1_5a = -0.5*C2_R*dimrep*g_loop/(1.0d0*l**3)


	END SUBROUTINE g15a




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE g15b(l,f1,f2,g1_5b)
	!
	! Subroutine to compute diagram g1_5b 
	!
	implicit none

	integer,intent(in) :: l, f1, f2
	complex(kind=8),intent(out) :: g1_5b

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,0:l,0:l) :: SQSQS
	complex(kind=8),dimension(4,4) :: V1mu, V1nu	
	complex(kind=8),dimension(4,4) :: mat1	
	integer,dimension(3) :: qloop, p_ext
	complex(kind=8) :: g_loop
	integer :: mu, nu, i
	integer :: n1, n2, n3
	integer :: s0,s0p,t0,t0p,u0,u0p
	integer :: t0_a,t0_b,t0p_a,t0p_b,u0_a,u0_b,u0p_a,u0p_b
	real(kind=8) :: multi
	complex(kind=8) :: trace
	integer :: f, fp

	!Variables for the openMP paralelization
	integer :: chunk
	chunk = 1

	f = f1f2(f1,f2)
	fp = f1f2(f2,f1)

	!initialization
	g1_5b = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	SQSQS(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	p_ext(:) = 0
	trace = dcmplx(0.0d0,0.0d0)
	g_loop = dcmplx(0.0d0,0.0d0)


	!part outside the loop
	call Fermion_Propagator(l,p_ext,S_prop_0)
	mat1 = matmul(matmul(Q5p(:,:,fp),S_prop_0(:,:,l-1,1,f1)),Q5(:,:,f))

	do t0=0,l
	do s0=0,l
	  SQSQS(:,:,t0,s0) = matmul(matmul(S_prop_0(:,:,t0,l-1,f2),mat1),S_prop_0(:,:,1,s0,f2))
	enddo
	enddo

	!Construct  loop

!$OMP PARALLEL DO &
!$OMP DEFAULT(PRIVATE) SHARED(l,p_ext,chunk,f1,f2,SQSQS,multiplicity,mom_deg) &
!$OMP REDUCTION(+:g_loop) &
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

	trace = Trace4(SQSQS(:,:,t0p,s0),V1mu(:,:),S_prop(:,:,t0,s0p,f2),V1nu(:,:))
	g_loop = g_loop + D_prop(nu,mu,u0p,u0)*trace*multi
	

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



	g1_5b = -0.5*C2_R*dimrep*g_loop/(1.0d0*l**3)


	END SUBROUTINE g15b
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE g16a(l,f1,f2,g1_6a)
	!
	! Subroutine to compute diagram g1_6a
	!
	implicit none

	integer,intent(in) :: l, f1, f2
	complex(kind=8),intent(out) :: g1_6a

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop_0
	complex(kind=8),dimension(4,4,0:l,0:l) :: SQSQS
	complex(kind=8),dimension(4,4) :: V2mu	
	complex(kind=8),dimension(4,4) :: mat1	
	integer,dimension(3) :: qloop, p_ext
	complex(kind=8) :: g_loop
	integer :: mu, i
	integer :: n1, n2, n3
	integer :: s0,t0,u0
	integer :: t0_a,t0_b,u0_a,u0_b
	real(kind=8) :: multi
	complex(kind=8) :: trace
	integer :: f, fp

	!Variables for the openMP paralelization
	integer :: chunk
	chunk = 1


	f = f1f2(f1,f2)
	fp = f1f2(f2,f1)

	!initialization
	g1_6a = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	SQSQS(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	V2mu(:,:) = dcmplx(0.0d0,0.0d0)
	p_ext(:) = 0
	trace = dcmplx(0.0d0,0.0d0)
	g_loop = dcmplx(0.0d0,0.0d0)


	!part outside the loop
	call Fermion_Propagator(l,p_ext,S_prop_0)
	mat1 = matmul(matmul(Q5(:,:,f),S_prop_0(:,:,1,l-1,f2)),Q5p(:,:,fp))

	do t0=0,l
	do s0=0,l
	  SQSQS(:,:,t0,s0) = matmul(matmul(S_prop_0(:,:,t0,1,f1),mat1),S_prop_0(:,:,l-1,s0,f1))
	enddo
	enddo

	!Construct  loop

!$OMP PARALLEL DO &
!$OMP DEFAULT(PRIVATE) SHARED(l,chunk,p_ext,f1,f2,SQSQS,multiplicity,mom_deg) &
!$OMP REDUCTION(+:g_loop) &
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

	trace = Trace2(SQSQS(:,:,t0,s0),V2mu(:,:))
	g_loop = g_loop + D_prop(mu,mu,u0,u0)*trace*multi
	

	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu


	enddo	!q3
	enddo	!q2
	enddo	!q1

!$OMP END PARALLEL DO

	g1_6a = 0.5d0*C2_R*dimrep*g_loop/(1.0d0*l**3)


	END SUBROUTINE g16a




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE g16b(l,f1,f2,g1_6b)
	!
	! Subroutine to compute diagram g1_6b 
	!
	implicit none

	integer,intent(in) :: l, f1, f2
	complex(kind=8),intent(out) :: g1_6b

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop_0
	complex(kind=8),dimension(4,4,0:l,0:l) :: SQSQS
	complex(kind=8),dimension(4,4) :: V2mu	
	complex(kind=8),dimension(4,4) :: mat1	
	integer,dimension(3) :: qloop, p_ext
	complex(kind=8) :: g_loop
	integer :: mu, i
	integer :: n1, n2, n3
	integer :: s0,t0,u0
	integer :: t0_a,t0_b,u0_a,u0_b
	real(kind=8) :: multi
	complex(kind=8) :: trace
	integer :: f, fp
	!Variables for the openMP paralelization
	integer :: chunk
	chunk = 1


	f = f1f2(f1,f2)
	fp = f1f2(f2,f1)

	!initialization
	g1_6b = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	SQSQS(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	V2mu(:,:) = dcmplx(0.0d0,0.0d0)
	p_ext(:) = 0
	trace = dcmplx(0.0d0,0.0d0)
	g_loop = dcmplx(0.0d0,0.0d0)


	!part outside the loop
	call Fermion_Propagator(l,p_ext,S_prop_0)
	mat1 = matmul(matmul(Q5p(:,:,fp),S_prop_0(:,:,l-1,1,f1)),Q5(:,:,f))

	do t0=0,l
	do s0=0,l
	  SQSQS(:,:,t0,s0) = matmul(matmul(S_prop_0(:,:,t0,l-1,f2),mat1),S_prop_0(:,:,1,s0,f2))
	enddo
	enddo

	!Construct  loop
!$OMP PARALLEL DO &
!$OMP DEFAULT(PRIVATE) SHARED(l,chunk,p_ext,f1,f2,SQSQS,multiplicity,mom_deg) &
!$OMP REDUCTION(+:g_loop) &
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

	trace = Trace2(SQSQS(:,:,t0,s0),V2mu(:,:))
	g_loop = g_loop + D_prop(mu,mu,u0,u0)*trace*multi
	

	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu


	enddo	!q3
	enddo	!q2
	enddo	!q1
!$OMP END PARALLEL DO

	g1_6b = 0.5d0*C2_R*dimrep*g_loop/(1.0d0*l**3)


	END SUBROUTINE g16b




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



	SUBROUTINE g17(l,f1,f2,g1_7)
	!
	! Subroutine to compute diagram g1_7 
	!
	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8),intent(out) :: g1_7
	

	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l) :: S_Q5_S,S_Q5p_S	
	complex(kind=8),dimension(4,4) :: V1mu, V1nu	
	integer :: x0,i
	integer :: s0,s0p,t0,t0p,u0,u0p
	integer :: t0_a,t0_b,t0p_a,t0p_b,u0_a,u0_b,u0p_a,u0p_b
	integer :: mu, nu
	integer :: n1, n2, n3
	integer, dimension(3) :: p,q
	complex(kind=8) :: trace
	real(kind=8) :: multi
	integer :: f,fp

	!initialization
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	g1_7 = dcmplx(0.0d0,0.0d0)
	S_Q5_S(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_Q5p_S(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	p(:) = 0
	f=f1f2(f1,f2)
	fp=f1f2(f2,f1)

	!build S_O5_S

	call Fermion_Propagator(l,p,S_prop_0)

	do s0=0,l
	do t0=0,l
	    S_Q5_S(:,:,t0,s0) = matmul(S_prop_0(:,:,t0,1,f1),matmul(Q5(:,:,f),S_prop_0(:,:,1,s0,f2)))
	    S_Q5p_S(:,:,t0,s0) = matmul(S_prop_0(:,:,t0,l-1,f2),matmul(Q5p(:,:,fp),S_prop_0(:,:,l-1,s0,f1)))
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

	  trace = Trace4(S_Q5_S(:,:,t0p,s0),V1mu(:,:),S_Q5p_S(:,:,t0,s0p),V1nu(:,:))
	  g1_7 = g1_7 + trace*D_prop(nu,mu,u0p,u0)


	  enddo	!nu
	  enddo !mu
	  enddo	!u0p
	  enddo	!u0
	enddo	!t0p
	enddo	!t0
	enddo	!s0p
	enddo	!s0


 	g1_7 = -0.5d0*dimrep*C2_R*g1_7/((l**3)*1.0d0)

	
	END SUBROUTINE g17




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE g18a(l,f1,f2,g1_8a)
	!
	! Subroutine to compute diagram g1_8a 
	!
	implicit none

	integer,intent(in) :: l, f1, f2
	complex(kind=8),intent(out) :: g1_8a



	integer :: n1,n2,n3
	integer,dimension(3) :: qloop
	complex(kind=8) :: loop
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop!,K_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop
	complex(kind=8) :: g1_0
	real(kind=8) :: multi
	integer,dimension(3) :: nvec
	integer :: f,fp

	f = f1f2(f2,f1)
	fp = f1f2(f1,f2)
	nvec(:) = 0

	!initialization
	g1_8a = dcmplx(0.0d0,0.0d0)
	loop = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	g1_0 = dcmplx(0.0d0,0.0d0)

	! Construct tree-level part:
	call Fermion_Propagator(l,nvec,S_Prop)
	g1_0 = 0.5d0*dimrep*Trace4(Q5(:,:,f),S_prop(:,:,1,l-1,f1),Q5p(:,:,fp),S_prop(:,:,l-1,1,f2))

	

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

	g1_8a = -0.5d0*C2_R*loop*g1_0/((l*1.0d0)**3)

	END SUBROUTINE g18a
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE g19(l,f1,f2,g1_9)
	!
	! Subroutine to compute diagram g1_9
	!
	implicit none

	integer,intent(in) :: l, f1, f2
	complex(kind=8),intent(out) :: g1_9


	integer,dimension(3) :: qloop,nvec
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop!,K_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop
	complex(kind=8) :: g1_0
	integer :: f, fp

	nvec(:) = 0

	f = f1f2(f2,f1)
	fp = f1f2(f1,f2)

	!initialization
	g1_9 = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	g1_0 = dcmplx(0.0d0,0.0d0)

	! Construct tree-level part:
	call Fermion_Propagator(l,nvec,S_Prop)
	g1_0 = 0.5d0*dimrep*Trace4(Q5(:,:,f),S_prop(:,:,1,l-1,f1),Q5p(:,:,fp),S_prop(:,:,l-1,1,f2))

	!Construct the loop
	qloop(1)=0
	qloop(2)=0
	qloop(3)=0

  	call Gluon_propagator(l,qloop,D_prop)	

	g1_9 =  C2_R*g1_0*D_prop(0,0,l-1,l-1)/((l*1.0d0)**3)

	END SUBROUTINE g19

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE g110a(l,f1,f2,g1_10a)
	!
	! Subroutine to compute diagram g1_10a 
	!
	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8),intent(out) :: g1_10a
	
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, s_prop_0
	complex(kind=8),dimension(4,4) :: Vmu1	
	integer,dimension(3) :: qloop, p_ext
	complex(kind=8),dimension(4,4) :: mat1!, mat2, mat3	!s,s,f1f2
	real(kind=8) :: multi
	integer :: mu,u0,s0,t0,n1,n2,n3
	integer :: t0_a, t0_b, u0_a, u0_b
	complex(kind=8) :: Trace
	integer :: f, fp

	f = f1f2(f2,f1) 
	fp = f1f2(f1,f2) 


	!initialization
	p_ext(:) = 0
	qloop(:)=0
	g1_10a = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	p_ext(:) = 0
	mat1(:,:) = dcmplx(0.0d0,0.0d0)
	Vmu1(:,:) = dcmplx(0.0d0,0.0d0)
	Trace = dcmplx(0.0d0,0.0d0)

	call Fermion_Propagator(l,p_ext,S_prop_0)
	mat1(:,:) = matmul(Q5(:,:,fp),matmul(S_prop_0(:,:,1,l-1,f2),Q5p(:,:,f)))


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

	Trace = Trace + Trace4(mat1(:,:),S_prop(:,:,l-1,s0,f1),Vmu1,S_prop_0(:,:,t0,1,f1))*D_prop(mu,0,u0,l-1)*multi


	enddo	!u0
	enddo	!t0
	enddo	!s0
	
	
	enddo	!mu


	enddo	!momenta
	enddo
	enddo	


	g1_10a = -0.5d0*C2_R*dimrep*Trace/((l*1.0d0)**3)

	END SUBROUTINE g110a

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE g110b(l,f1,f2,g1_10b)
	!
	! Subroutine to compute diagram g1_10b 
	!
	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8),intent(out) :: g1_10b
	
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, s_prop_0
	complex(kind=8),dimension(4,4) :: Vmu1	
	integer,dimension(3) :: qloop, p_ext
	complex(kind=8),dimension(4,4) :: mat1!, mat2, mat3	!s,s,f1f2
	real(kind=8) :: multi
	integer :: mu,u0,s0,t0,n1,n2,n3
	integer :: t0_a, t0_b, u0_a, u0_b
	complex(kind=8) :: Trace
	integer :: f, fp

	f = f1f2(f2,f1) 
	fp = f1f2(f1,f2) 


	!initialization
	p_ext(:) = 0
	qloop(:)=0
	g1_10b = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	p_ext(:) = 0
	mat1(:,:) = dcmplx(0.0d0,0.0d0)
	Vmu1(:,:) = dcmplx(0.0d0,0.0d0)
	Trace = dcmplx(0.0d0,0.0d0)

	call Fermion_Propagator(l,p_ext,S_prop_0)
	mat1(:,:) = matmul(Q5p(:,:,f),matmul(S_prop_0(:,:,l-1,1,f1),Q5(:,:,fp)))


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

	Trace = Trace + Trace4(mat1(:,:),S_prop_0(:,:,1,s0,f2),Vmu1,S_prop(:,:,t0,l-1,f2))*D_prop(0,mu,l-1,u0)*multi


	enddo	!u0
	enddo	!t0
	enddo	!s0
	
	
	enddo	!mu


	enddo	!momenta
	enddo
	enddo	


	g1_10b = 0.5d0*C2_R*dimrep*Trace/((l*1.0d0)**3)

	END SUBROUTINE g110b

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE g111a(l,f1,f2,g1_11a)
	!
	! Subroutine to compute diagram g1_11a 
	!
	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8),intent(out) :: g1_11a
	
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: s_prop_0
	complex(kind=8),dimension(4,4) :: Vmu1	
	integer,dimension(3) :: qloop, p_ext
	complex(kind=8),dimension(4,4) :: mat1!, mat2, mat3	!s,s,f1f2
	real(kind=8) :: multi
	integer :: mu,u0,s0,t0,n1,n2,n3
	integer :: t0_a, t0_b, u0_a, u0_b
	complex(kind=8) :: Trace
	integer :: f, fp

	f = f1f2(f2,f1) 
	fp = f1f2(f1,f2) 


	!initialization
	p_ext(:) = 0
	qloop(:)=0
	g1_11a = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	p_ext(:) = 0
	mat1(:,:) = dcmplx(0.0d0,0.0d0)
	Vmu1(:,:) = dcmplx(0.0d0,0.0d0)
	Trace = dcmplx(0.0d0,0.0d0)

	call Fermion_Propagator(l,p_ext,S_prop_0)
	mat1(:,:) = matmul(Q5p(:,:,f),matmul(S_prop_0(:,:,l-1,1,f1),Q5(:,:,fp)))


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

	Trace = Trace + Trace4(S_prop_0(:,:,1,s0,f2),Vmu1(:,:),S_prop_0(:,:,t0,l-1,f2),mat1(:,:))*D_prop(0,mu,l-1,u0)


	enddo	!u0
	enddo	!t0
	enddo	!s0
	
	
	enddo	!mu


	g1_11a = -0.5d0*C2_R*dimrep*Trace/((l*1.0d0)**3)

	END SUBROUTINE g111a

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



	SUBROUTINE g111b(l,f1,f2,g1_11b)
	!
	! Subroutine to compute diagram g1_11b 
	!
	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8),intent(out) :: g1_11b
	
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: s_prop_0
	complex(kind=8),dimension(4,4) :: Vmu1	
	integer,dimension(3) :: qloop, p_ext
	complex(kind=8),dimension(4,4) :: mat1!, mat2, mat3	!s,s,f1f2
	real(kind=8) :: multi
	integer :: mu,u0,s0,t0,n1,n2,n3
	integer :: t0_a, t0_b, u0_a, u0_b
	complex(kind=8) :: Trace
	integer :: f, fp

	f = f1f2(f2,f1) 
	fp = f1f2(f1,f2) 


	!initialization
	p_ext(:) = 0
	qloop(:)=0
	g1_11b = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	p_ext(:) = 0
	mat1(:,:) = dcmplx(0.0d0,0.0d0)
	Vmu1(:,:) = dcmplx(0.0d0,0.0d0)
	Trace = dcmplx(0.0d0,0.0d0)

	call Fermion_Propagator(l,p_ext,S_prop_0)
	mat1(:,:) = matmul(Q5(:,:,fp),matmul(S_prop_0(:,:,1,l-1,f2),Q5p(:,:,f)))


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

	Trace = Trace + Trace4(S_prop_0(:,:,l-1,s0,f1),Vmu1(:,:),S_prop_0(:,:,t0,1,f1),mat1(:,:))*D_prop(mu,0,u0,l-1)


	enddo	!u0
	enddo	!t0
	enddo	!s0
	
	
	enddo	!mu


	g1_11b = 0.5d0*C2_R*dimrep*Trace/((l*1.0d0)**3)

	END SUBROUTINE g111b

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE g112a(l,f1,f2,g1_12a)
	!
	! Subroutine to compute diagram g1_12a 
	!
	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8),intent(out) :: g1_12a
	
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	integer,dimension(3) :: qloop, p_ext
	complex(kind=8),dimension(4,4) :: mat1
	real(kind=8) :: multi
	integer :: n1,n2,n3
	complex(kind=8) :: Trace
	integer :: f, fp

	f =f1f2(f2,f1)
	fp =f1f2(f1,f2)

	!initialization
	g1_12a = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	mat1(:,:) = dcmplx(0.0d0,0.0d0)
	Trace = dcmplx(0.0d0,0.0d0)
	p_ext(:) =0


	call Fermion_Propagator(l,p_ext,S_prop_0)
	mat1(:,:) = matmul(Q5(:,:,fp),matmul(S_prop_0(:,:,1,l-1,f2),Q5p(:,:,f)))

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

	Trace = Trace + Trace2(mat1(:,:),S_prop(:,:,l-1,1,f1))*D_prop(0,0,0,l-1)*multi

	enddo	!enddo
	enddo
	enddo	

	g1_12a = -0.5d0*C2_R*dimrep*Trace/((l*1.0d0)**3)

	END SUBROUTINE g112a


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



	SUBROUTINE g112b(l,f1,f2,g1_12b)
	!
	! Subroutine to compute diagram g1_12b
	!
	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8),intent(out) :: g1_12b
	
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	integer,dimension(3) :: qloop, p_ext
	complex(kind=8),dimension(4,4) :: mat1
	real(kind=8) :: multi
	integer :: n1,n2,n3
	complex(kind=8) :: Trace
	integer :: f, fp

	f =f1f2(f2,f1)
	fp =f1f2(f1,f2)

	!initialization
	g1_12b = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	mat1(:,:) = dcmplx(0.0d0,0.0d0)
	Trace = dcmplx(0.0d0,0.0d0)
	p_ext(:) =0


	call Fermion_Propagator(l,p_ext,S_prop_0)
	mat1(:,:) = matmul(Q5p(:,:,f),matmul(S_prop_0(:,:,l-1,1,f1),Q5(:,:,fp)))

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

	Trace = Trace + Trace2(mat1(:,:),S_prop(:,:,1,l-1,f2))*D_prop(0,0,l-1,0)*multi

	enddo	!enddo
	enddo
	enddo	

	g1_12b = -0.5d0*C2_R*dimrep*Trace/((l*1.0d0)**3)

	END SUBROUTINE g112b


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE g113a(l,f1,f2,g1_13a)
	!
	! Subroutine to compute diagram g1_13a 
	!
	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8),intent(out) :: g1_13a

	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop_0	
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	integer,dimension(3) :: qloop, p_ext
	complex(kind=8),dimension(4,4) :: mat1!, mat2
	integer :: f, fp


	f =f1f2(f2,f1)
	fp =f1f2(f1,f2)

	!initialization
	g1_13a = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	mat1(:,:) = dcmplx(0.0d0,0.0d0)
	qloop(:) = 0
	p_ext(:) = 0

	call Fermion_Propagator(l,p_ext,S_prop_0)
	mat1(:,:) = matmul(Q5(:,:,fp),matmul(S_prop_0(:,:,1,l-1,f2),Q5p(:,:,f)))
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	call Gluon_propagator(l,qloop,D_prop)	

	g1_13a = 0.5d0*C2_R*dimrep*Trace2(mat1(:,:),S_prop_0(:,:,l-1,1,f1))*D_prop(0,0,l-1,0)/((l*1.0d0)**3)


	END SUBROUTINE g113a


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE g1dm(l,f1,f2,g1_dm)
	!
	! Subroutine to compute diagram g1_dm 
	!
	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8),intent(out) :: g1_dm
	
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop
	complex(kind=8) :: g1_dm_1, g1_dm_2
	integer, dimension(3) :: p_ext
	integer :: f,fp,s0

	f=f1f2(f1,f2)
	fp=f1f2(f2,f1)
	p_ext(:) = 0

	!initialization
	g1_dm = dcmplx(0.0d0,0.0d0)
	g1_dm_1 = dcmplx(0.0d0,0.0d0)
	g1_dm_2 = dcmplx(0.0d0,0.0d0)

	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)


	call Fermion_Propagator(l,p_ext,S_prop)

	do s0 = 0,l

	g1_dm_1 = g1_dm_1-dimrep*Trace5(Q5(:,:,f),S_prop(:,:,1,l-1,f2),Q5p(:,:,fp),S_prop(:,:,l-1,s0,f1),S_prop(:,:,s0,1,f1))
	g1_dm_2 = g1_dm_2-dimrep*Trace5(Q5(:,:,f),S_prop(:,:,1,s0,f2),S_prop(:,:,s0,l-1,f2),Q5p(:,:,fp),S_prop(:,:,l-1,1,f1))

	enddo

	g1_dm = (g1_dm_1 + g1_dm_2)/2.0d0


	END SUBROUTINE g1dm


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE g1zf(l,f1,f2,g1_zf)
	!
	! Subroutine to compute diagram g1_zf
	!
	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8),intent(out) :: g1_zf
	
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop
	complex(kind=8) :: g1_zf_1, g1_zf_2
	integer, dimension(3) :: p_ext
	integer :: f,fp

	f=f1f2(f1,f2)
	fp=f1f2(f2,f1)
	p_ext(:) = 0

	!initialization
	g1_zf = dcmplx(0.0d0,0.0d0)
	g1_zf_1 = dcmplx(0.0d0,0.0d0)
	g1_zf_2 = dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)


	call Fermion_Propagator(l,p_ext,S_prop)

	g1_zf_1 = Trace5(Q5(:,:,f),S_prop(:,:,1,l-1,f2),Q5p(:,:,fp),S_prop(:,:,l-1,0,f1),S_prop(:,:,0,1,f1))
	g1_zf_1 = g1_zf_1 + Trace5(Q5(:,:,f),S_prop(:,:,1,l-1,f2),Q5p(:,:,fp),S_prop(:,:,l-1,l,f1),S_prop(:,:,l,1,f1))

	g1_zf_2 = Trace5(Q5(:,:,f),S_prop(:,:,1,0,f2),S_prop(:,:,0,l-1,f2),Q5p(:,:,fp),S_prop(:,:,l-1,1,f1))
	g1_zf_2 = g1_zf_2 + Trace5(Q5(:,:,f),S_prop(:,:,1,l,f2),S_prop(:,:,l,l-1,f2),Q5p(:,:,fp),S_prop(:,:,l-1,1,f1))


	g1_zf = -dimrep*(g1_zf_1+g1_zf_2)/2.0d0

	END SUBROUTINE g1zf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE g1ds(l,f1,f2,g1_ds)
	!
	! Subroutine to compute diagram g1_ds
	!
	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8),intent(out) :: g1_ds
	
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop
	complex(kind=8) :: g1_ds_1, g1_ds_2
	integer, dimension(3) :: p_ext
	real(kind=8),dimension(3) ::  pe_plus, ptwit, phat
	complex(kind=8),dimension(4,4) :: Amat
	integer :: f,fp

	f=f1f2(f1,f2)
	fp=f1f2(f2,f1)
	p_ext(:) = 0

	!initialization
	g1_ds = dcmplx(0.0d0,0.0d0)
	g1_ds_1 = dcmplx(0.0d0,0.0d0)
	g1_ds_2 = dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)


	pe_plus(:) = p_ext(:)*(2.0d0*pi)/(l*1.0d0) + theta/l
	ptwit(:) = sin(pe_plus(:))
	phat(:) = 2.0d0*sin(pe_plus(:)/2.0d0)

	Amat = imag*(ptwit(1)*g1+ptwit(2)*g2+ptwit(3)*g3)
	Amat = Amat + 0.5d0*II*( phat(1)**2 + phat(2)**2 + phat(3)**2 )



	call Fermion_Propagator(l,p_ext,S_prop)

	g1_ds_1 = Trace6(Q5(:,:,f),S_prop(:,:,1,l-1,f2),Q5p(:,:,fp),S_prop(:,:,l-1,0,f1),Amat,S_prop(:,:,0,1,f1))
	g1_ds_1 = g1_ds_1 + Trace6(Q5(:,:,f),S_prop(:,:,1,l-1,f2),Q5p(:,:,fp),S_prop(:,:,l-1,l,f1),Amat,S_prop(:,:,l,1,f1))


	g1_ds_2 = Trace6(Q5(:,:,f),S_prop(:,:,1,0,f2),Amat,S_prop(:,:,0,l-1,f2),Q5p(:,:,fp),S_prop(:,:,l-1,1,f1))
	g1_ds_2 = g1_ds_2 + Trace6(Q5(:,:,f),S_prop(:,:,1,l,f2),Amat,S_prop(:,:,l,l-1,f2),Q5p(:,:,fp),S_prop(:,:,l-1,1,f1))


	g1_ds = -dimrep*(g1_ds_1+g1_ds_2)/2.0d0

	END SUBROUTINE g1ds


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






	REAL(kind=8) FUNCTION N_g1(f1,f2)
	integer :: f1, f2
	if(f1.eq.f2)then
	  N_g1 = -1.0d0
	else
	  N_g1 = 1.0d0
	endif
	END FUNCTION N_g1


	REAL(kind=8) FUNCTION M_g1(f1)
	integer :: f1
	if(f1.eq.1)then
	  M_g1 = -1.0d0
	else
	  M_g1 = 1.0d0
	endif
	END FUNCTION M_g1


	INTEGER FUNCTION gX1_g1(f1,f2)
	integer :: f1,f2
	if(f1.eq.f2)then
	  gX1_g1 = 1
	else
	  gX1_g1 = 2
	endif
	END FUNCTION gX1_g1

	INTEGER FUNCTION gX2_g1(f1,f2)
	integer :: f1,f2
	if(f1.eq.f2)then
	  gX2_g1 = 4
	else
	  gX2_g1 = 3
	endif
	END FUNCTION gX2_g1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	END MODULE OneLoop_g1
