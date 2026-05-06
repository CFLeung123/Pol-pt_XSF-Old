



	MODULE OneLoop_gX
	!
	! Contains the set of routines for the calculaton of the 1-loop
	! diagrams contributing to the gX correlation functions (including
	! the contribution of counterterms.
	!
	USE GammaMatrices
	USE Input
	USE Parameters
	USE Vertices
	USE Propagators
	USE TreeLevel_gX
	USE omp_lib

	implicit none

	CONTAINS


	SUBROUTINE gX1a(l,g_i,f1,f2,gX_1a)
	!
	! Subroutine to compute diagram gX_1a 
	!
	implicit none

	integer,intent(in) :: l, g_i,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: gX_1a
	
	integer :: n1,n2,n3
	integer,dimension(3) :: qloop
	complex(kind=8) :: loop
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop!,K_prop
	complex(kind=8),dimension(c0:c1) :: gX_0
	real(kind=8) :: multi	


	!initialization
	gX_1a(:) = dcmplx(0.0d0,0.0d0)
	loop = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	gX_0(:) = dcmplx(0.0d0,0.0d0)

	!Tree level part
	call gX0(l,g_i,f1,f2,gX_0)

	!Construct the gluon loop
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

	gX_1a(:) = -0.5d0*C2_R*loop*gX_0(:)/((l*1.0d0)**3)

	END SUBROUTINE gX1a

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



	SUBROUTINE gX2(l,g_i,f1,f2,gX_2)
	!
	! Subroutine to compute diagram gX_2
	!
	implicit none

	integer,intent(in) :: l,g_i,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: gX_2
	
	integer :: x0,n1,n2,n3,i
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: loop, gX_0q
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop!,K_prop	
	complex(kind=8),dimension(4,4) :: mat1, mat2, mat3
	real(kind=8) :: multi	
	integer :: f

	!initialization
	gX_2(:) = dcmplx(0.0d0,0.0d0)
	loop = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	loop(:) = dcmplx(0.0d0,0.0d0)
	gX_0q(:) = dcmplx(0.0d0,0.0d0)
	mat1(:,:) = dcmplx(0.0d0,0.0d0)
	mat2(:,:) = dcmplx(0.0d0,0.0d0)
	mat3(:,:) = dcmplx(0.0d0,0.0d0)

	f=f1f2(f2,f1)
	!Construct the gluon loop
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

	   !Construction of the tree level f with momentum q: 

	     do x0=c0,c1

	       mat1 = matmul(Q5(:,:,f),S_prop(:,:,1,x0,f1))
	       mat2 = matmul(mat1,gX(:,:,g_i))
	       mat3 = matmul(mat2,S_prop(:,:,x0,1,f2))	

	       gX_0q(x0) = 0.5d0*dimrep*(mat3(1,1)+mat3(2,2)+mat3(3,3)+mat3(4,4))
	       loop(x0) = loop(x0) + D_prop(0,0,0,0)*gX_0q(x0)*multi

	     enddo


	enddo
	enddo
	enddo

	  gX_2(:) = C2_R*loop(:)/((l*1.0d0)**3)


	END SUBROUTINE gX2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




	SUBROUTINE gX3a(l,g_i,f1,f2,gX_3a)
	!
	! Subroutine to compute diagram gX_3a 
	!
	implicit none

	integer,intent(in) :: l,g_i,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: gX_3a

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,c0:c1,0:l) :: S_G
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop, p_ext
	complex(kind=8),dimension(c0:c1) :: g_loop
	integer :: n1, n2, n3, x0, t0, u0, s0, mu,i, x0_t
	integer :: t0_a, t0_b, u0_a, u0_b
	real(kind=8) :: multi
	complex(kind=8) :: trace
	integer :: f

	gX_3a(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_G(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	p_ext(:) = 0
	trace = dcmplx(0.0d0,0.0d0)
	f=f1f2(f2,f1)

	!part outside the loop
	call Fermion_Propagator(l,p_ext,S_prop_0)

	do x0 = c0,c1
	do s0=0,l
		
	    S_G(:,:,x0,s0) = matmul(S_prop_0(:,:,1,x0,f1),matmul(gX(:,:,g_i),S_prop_0(:,:,x0,s0,f2)))

	enddo
	enddo

	!Construct the gluon loop
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


	do mu = mu_low(0),mu_up(0)

	do s0 = 0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0 = t0_a,t0_b
	do u0 = u0_a,u0_b

	call V1(l,p_ext,-p_ext-qloop,qloop,s0,t0,u0,mu,V1mu)
	
	do x0 = c0,c1
	  trace = Trace4(S_G(:,:,x0,s0),V1mu(:,:),S_prop(:,:,t0,1,f2),Q5(:,:,f))
	  g_loop(x0) = g_loop(x0) + D_prop(0,mu,0,u0)*trace*multi
	enddo

	
	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu


	enddo	!momenta
	enddo
	enddo

	  gX_3a(:) = -0.5*C2_R*dimrep*g_loop(:)/(1.0d0*l**3)

	END SUBROUTINE gX3a



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE gX3b(l,g_i,f1,f2,gX_3b)
	!
	! Subroutine to compute diagram gX_3b
	!
	implicit none

	integer,intent(in) :: l,g_i,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: gX_3b

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,0:l,c0:c1) :: S_G
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop, p_ext
	complex(kind=8),dimension(c0:c1) :: g_loop
	integer :: n1, n2, n3, x0, t0, u0, s0, mu,i, x0_t
	integer :: t0_a, t0_b, u0_a, u0_b
	real(kind=8) :: multi
	complex(kind=8) :: trace
	integer :: f

	gX_3b(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_G(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	p_ext(:) = 0
	trace = dcmplx(0.0d0,0.0d0)
	f=f1f2(f2,f1)
	!part outside the loop
	call Fermion_Propagator(l,p_ext,S_prop_0)

	do t0 = 0,l
	do x0 = c0,c1		
	    S_G(:,:,t0,x0) = matmul(S_prop_0(:,:,t0,x0,f1),matmul(gX(:,:,g_i),S_prop_0(:,:,x0,1,f2)))
	enddo
	enddo


	!Construct the gluon loop
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
	
	do x0=c0,c1
	  trace = Trace4(Q5(:,:,f),S_prop(:,:,1,s0,f1),V1mu(:,:),S_G(:,:,t0,x0))
	  g_loop(x0) = g_loop(x0) + D_prop(mu,0,u0,0)*trace*multi
	enddo

	
	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu


	enddo	!momenta
	enddo
	enddo

	  gX_3b(:) = 0.5*C2_R*dimrep*g_loop(:)/(1.0d0*l**3)

	END SUBROUTINE gX3b



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE gX4a(l,g_i,f1,f2,gX_4a)
	!
	! Subroutine to compute diagram gX_4a 
	!
	implicit none

	integer,intent(in) :: l, g_i,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: gX_4a

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop,S_prop_0
	complex(kind=8),dimension(4,4,0:l,c0:c1) :: S_G
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop, p_ext
	complex(kind=8),dimension(c0:c1) :: g_loop
	integer :: n1, n2, n3, x0, t0, u0, s0, mu, i
	integer :: t0_a, t0_b, u0_a, u0_b
	real(kind=8) :: multi 
	complex(kind=8) :: trace
	integer :: f

	gX_4a(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_G(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	p_ext(:) = 0
	trace = dcmplx(0.0d0,0.0d0)
	f=f1f2(f2,f1)

	!fermion propagator 0 momenta
	call Fermion_Propagator(l,p_ext,S_prop_0)

	!Construct the gluon loop
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

	   !Build S_G
	   do t0=0,l
	   do x0=c0,c1
	      S_G(:,:,t0,x0) = matmul(S_prop(:,:,t0,x0,f1),matmul(gX(:,:,g_i),S_prop(:,:,x0,1,f2)))
	   enddo
	   enddo

	do mu=mu_low(0),mu_up(0)

	do s0=0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0=t0_a,t0_b
	do u0=u0_a,u0_b

	call V1(l,p_ext,-p_ext-qloop,qloop,s0,t0,u0,mu,V1mu)

	do x0=c0,c1
	  trace = Trace4(Q5(:,:,f),S_prop_0(:,:,1,s0,f1),V1mu(:,:),S_G(:,:,t0,x0))
	  g_loop(x0) = g_loop(x0) + D_prop(0,mu,0,u0)*trace*multi
	enddo	!x0
	

	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu


	enddo	!momentum
	enddo
	enddo

	gX_4a(:) = -0.5*C2_R*dimrep*g_loop(:)/(1.0d0*l**3)

	END SUBROUTINE gX4a


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE gX4b(l,g_i,f1,f2,gX_4b)
	!
	! Subroutine to compute diagram gX_4b
	!
	implicit none

	integer,intent(in) :: l, g_i,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: gX_4b

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,c0:c1,0:l) :: S_G
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop, p_ext
	complex(kind=8),dimension(c0:c1) :: g_loop
	integer :: n1, n2, n3, x0, t0, u0, s0, mu, i
	integer :: t0_a, t0_b, u0_a, u0_b
	real(kind=8) :: multi
	complex(kind=8) :: trace
	integer :: f

	gX_4b(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_G(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	p_ext(:) = 0
	trace = dcmplx(0.0d0,0.0d0)
	f=f1f2(f2,f1)
	!fermion propagator 0 momenta
	call Fermion_Propagator(l,p_ext,S_prop_0)

	!Construct the gluon loop
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

	   !Build S_G
	   
	   do x0 = c0,c1
	   do s0=0,l
	      S_G(:,:,x0,s0) = matmul(S_prop(:,:,1,x0,f1),matmul(gX(:,:,g_i),S_prop(:,:,x0,s0,f2)))
	   enddo
	   enddo


	do mu=mu_low(0),mu_up(0)

	do s0=0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0=t0_a,t0_b
	do u0=u0_a,u0_b

	call V1(l,p_ext+qloop,-p_ext,-qloop,s0,t0,u0,mu,V1mu)

	do x0 = c0,c1
	  trace = Trace4(Q5(:,:,f),S_G(:,:,x0,s0),V1mu(:,:),S_prop_0(:,:,t0,1,f2))
	  g_loop(x0) = g_loop(x0) + D_prop(mu,0,u0,0)*trace*multi
	enddo	!x0
	

	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu



	enddo	!momentum
	enddo
	enddo

	gX_4b(:) = 0.5*C2_R*dimrep*g_loop(:)/(1.0d0*l**3)

	END SUBROUTINE gX4b
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



	SUBROUTINE gX5a(l,g_i,f1,f2,gX_5a)
	!
	! Subroutine to compute diagram gX_5a 
	!
	implicit none

	integer,intent(in) :: l,g_i,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: gX_5a
	
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,c0:c1,0:l) :: S_O5_SGS	
	complex(kind=8),dimension(4,4) :: V1mu, V1nu	
	integer :: x0,i
	integer :: s0,s0p,t0,t0p,u0,u0p
	integer :: t0_a,t0_b,t0p_a,t0p_b,u0_a,u0_b,u0p_a,u0p_b
	integer :: mu, nu
	integer :: n1, n2, n3
	integer, dimension(3) :: p,q
	complex(kind=8) :: trace
	real(kind=8) :: multi
	integer :: f

	!Variables for the openMP paralelization
	integer :: chunk
	complex(kind=8) :: f_temp
	chunk = 1
	f_temp = dcmplx(0.0d0,0.0d0)

	!initialization
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_O5_SGS(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	gX_5a(:) = dcmplx(0.0d0,0.0d0)
	p(:) = 0
	f=f1f2(f2,f1)



	!build S_O5_SGS

	call Fermion_Propagator(l,p,S_prop_0)

	do x0 = c0,c1
	do s0=0,l
	do t0=0,l
	    S_O5_SGS(:,:,t0,x0,s0) = matmul(S_prop_0(:,:,t0,1,f2),matmul(matmul(Q5(:,:,f),S_prop_0(:,:,1,x0,f1)),matmul(gX(:,:,g_i),S_prop_0(:,:,x0,s0,f2))))
	enddo	
	enddo
	enddo


	do x0 = c0,c1

	f_temp = dcmplx(0.0d0,0.0d0)


!$OMP PARALLEL DO &
!$OMP DEFAULT(PRIVATE) SHARED(l,p,chunk,f1,f2,x0,S_O5_SGS,multiplicity,mom_deg) &
!$OMP REDUCTION(+:f_temp) &
!$OMP SCHEDULE(DYNAMIC,chunk) &
!$OMP COLLAPSE(1)


	do n1=0,l-1			
	do n2=n_limit(mom_deg,n1),l-1
	do n3=n_limit(mom_deg,n2),l-1

	q(1)=n1
	q(2)=n2
	q(3)=n3

	multi = 1.0d0*multiplicity(n1,n2,n3)

  	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

  	call Gluon_propagator(l,q,D_prop)	
	call Fermion_Propagator(l,p+q,S_prop)

	trace = dcmplx(0.0d0,0.0d0)
	
	do s0=0,l
	call loop_limits(l,s0,t0p_a,t0p_b,u0_a,u0_b)
	do s0p=0,l
	call loop_limits(l,s0p,t0_a,t0_b,u0p_a,u0p_b)
	do t0=t0_a,t0_b
	do t0p=t0p_a,t0p_b		!Due to the delta(x0+1,y0),delta(x0,y0),delta(x0-1,y0) in the vertex functions 
	do u0=u0_a,u0_b		!Due to the delta(x0+1,y0),delta(x0,y0),delta(x0-1,y0) in the vertex functions 
	do u0p=u0p_a,u0p_b

	  do mu=0,3	
	  do nu=mu_low(mu),mu_up(mu)

	  call V1(l,p,-p-q,q,s0,t0p,u0,mu,V1mu)
	  call V1(l,p+q,-p,-q,s0p,t0,u0p,nu,V1nu)

	  trace = Trace4(S_O5_SGS(:,:,t0,x0,s0),V1mu(:,:),S_prop(:,:,t0p,s0p,f2),V1nu(:,:))
	  f_temp = f_temp + multi*trace*D_prop(nu,mu,u0p,u0)



	  enddo	!nu
	  enddo !mu
	  enddo	!u0p
	  enddo	!u0
	enddo	!t0p
	enddo	!t0
	enddo	!s0p
	enddo	!s0



	enddo	!momenta
	enddo
	enddo

!$OMP END PARALLEL DO

	gX_5a(x0) = f_temp

	enddo	!x0

 	gX_5a(:) = -0.5d0*dimrep*C2_R*gX_5a(:)/((l**3)*1.0d0)

	
	END SUBROUTINE gX5a


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



	SUBROUTINE gX5b(l,g_i,f1,f2,gX_5b)
	!
	! Subroutine to compute diagram gX_5b
	!
	implicit none

	integer,intent(in) :: l,g_i,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: gX_5b
	

	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,c0:c1,0:l) :: SGS_O5_S	
	complex(kind=8),dimension(4,4) :: V1mu, V1nu	
	integer :: x0,i
	integer :: s0,s0p,t0,t0p,u0,u0p
	integer :: t0_a,t0_b,t0p_a,t0p_b,u0_a,u0_b,u0p_a,u0p_b
	integer :: mu, nu
	integer :: n1, n2, n3
	integer, dimension(3) :: p,q
	complex(kind=8) :: trace
	real(kind=8) :: multi
	integer :: f

	!Variables for the openMP paralelization
	integer :: chunk
	complex(kind=8) :: f_temp
	chunk = 1


	!initialization
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	SGS_O5_S(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	gX_5b(:) = dcmplx(0.0d0,0.0d0)
	p(:) = 0
	f=f1f2(f2,f1)
	f_temp = dcmplx(0.0d0,0.0d0)
	!build S_O5_SGS

	call Fermion_Propagator(l,p,S_prop_0)

	do x0 = c0,c1
	do s0=0,l
	do t0=0,l
	    SGS_O5_S(:,:,t0,x0,s0) = matmul(S_prop_0(:,:,t0,x0,f1),matmul(matmul(gX(:,:,g_i),S_prop_0(:,:,x0,1,f2)),matmul(Q5(:,:,f),S_prop_0(:,:,1,s0,f1))))
	enddo	
	enddo
	enddo



	do x0 = c0,c1

	f_temp = dcmplx(0.0d0,0.0d0)


!$OMP PARALLEL DO &
!$OMP DEFAULT(PRIVATE) SHARED(l,p,chunk,f1,f2,x0,SGS_O5_S,multiplicity,mom_deg) &
!$OMP REDUCTION(+:f_temp) &
!$OMP SCHEDULE(DYNAMIC,chunk) &
!$OMP COLLAPSE(1)


	do n1=0,l-1			
	do n2=n_limit(mom_deg,n1),l-1
	do n3=n_limit(mom_deg,n2),l-1

	q(1)=n1
	q(2)=n2
	q(3)=n3

	multi = 1.0d0*multiplicity(n1,n2,n3)

  	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

  	call Gluon_propagator(l,q,D_prop)	
	call Fermion_Propagator(l,p+q,S_prop)

	trace = dcmplx(0.0d0,0.0d0)
	
	do s0=0,l
	call loop_limits(l,s0,t0p_a,t0p_b,u0_a,u0_b)
	do s0p=0,l
	call loop_limits(l,s0p,t0_a,t0_b,u0p_a,u0p_b)
	do t0=t0_a,t0_b
	do t0p=t0p_a,t0p_b		!Due to the delta(x0+1,y0),delta(x0,y0),delta(x0-1,y0) in the vertex functions 
	do u0=u0_a,u0_b		!Due to the delta(x0+1,y0),delta(x0,y0),delta(x0-1,y0) in the vertex functions 
	do u0p=u0p_a,u0p_b

	  do mu=0,3	
	  do nu=mu_low(mu),mu_up(mu)

	  call V1(l,p,-p-q,q,s0,t0p,u0,mu,V1mu)
	  call V1(l,p+q,-p,-q,s0p,t0,u0p,nu,V1nu)

	  trace = Trace4(SGS_O5_S(:,:,t0,x0,s0),V1mu(:,:),S_prop(:,:,t0p,s0p,f1),V1nu(:,:))
	  f_temp = f_temp + multi*trace*D_prop(nu,mu,u0p,u0)


	  enddo	!nu
	  enddo !mu
	  enddo	!u0p
	  enddo	!u0
	enddo	!t0p
	enddo	!t0
	enddo	!s0p
	enddo	!s0

	enddo	!momenta
	enddo
	enddo

	!$OMP END PARALLEL DO

	gX_5b(x0) = f_temp

	enddo	!x0

 	gX_5b(:) = -0.5d0*dimrep*C2_R*gX_5b(:)/((l**3)*1.0d0)

	
	END SUBROUTINE gX5b


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



	SUBROUTINE gX6a(l,g_i,f1,f2,gX_6a)
	!
	! Subroutine to compute diagram gX_6a 
	!
	implicit none

	integer,intent(in) :: l,g_i,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: gX_6a
	

	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop_0
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,c0:c1,0:l) :: S_O5_SGS	
	complex(kind=8),dimension(4,4) :: V2mu	
	integer :: x0,i
	integer :: s0,t0,u0
	integer :: t0_a,t0_b,u0_a,u0_b
	integer :: mu, nu
	integer :: n1, n2, n3
	integer, dimension(3) :: p,q
	complex(kind=8) :: trace
	real(kind=8) :: multi
	integer :: f

	!Variables for the openMP paralelization
	integer :: chunk
	complex(kind=8) :: f_temp
	chunk = 1

	!initialization
	S_prop_0(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_O5_SGS(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	gX_6a(:) = dcmplx(0.0d0,0.0d0)
	p(:) = 0
	f=f1f2(f2,f1)
	f_temp = dcmplx(0.0d0,0.0d0)

	!build S_O5_SGS

	call Fermion_Propagator(l,p,S_prop_0)

	do x0 = c0,c1
	do s0=0,l
	do t0=0,l
	    S_O5_SGS(:,:,t0,x0,s0) = matmul(S_prop_0(:,:,t0,1,f2),matmul(matmul(Q5(:,:,f),S_prop_0(:,:,1,x0,f1)),matmul(gX(:,:,g_i),S_prop_0(:,:,x0,s0,f2))))
	enddo	
	enddo
	enddo


	do x0 = c0,c1

	f_temp = dcmplx(0.0d0,0.0d0)

!$OMP PARALLEL DO &
!$OMP DEFAULT(PRIVATE) SHARED(l,p,chunk,f1,f2,x0,S_O5_SGS,multiplicity,mom_deg) &
!$OMP REDUCTION(+:f_temp) &
!$OMP SCHEDULE(DYNAMIC,chunk) &
!$OMP COLLAPSE(1)



	do n1=0,l-1			
	do n2=n_limit(mom_deg,n1),l-1
	do n3=n_limit(mom_deg,n2),l-1

	q(1)=n1
	q(2)=n2
	q(3)=n3

	multi = 1.0d0*multiplicity(n1,n2,n3)

  	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)

  	call Gluon_propagator(l,q,D_prop)	

	trace = dcmplx(0.0d0,0.0d0)
	
	do s0=0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0=t0_a,t0_b
	do u0=u0_a,u0_b		!Due to the delta(x0+1,y0),delta(x0,y0),delta(x0-1,y0) in the vertex functions 

	  do mu=0,3	
!	  do nu=0,3

	  call V2(l,p,-p,q,-q,s0,t0,u0,u0,mu,V2mu)	!Remember that this guy has a delta(mu,nu) inside

	  trace = Trace2(S_O5_SGS(:,:,t0,x0,s0),V2mu(:,:))
	  f_temp = f_temp + multi*trace*D_prop(mu,mu,u0,u0)	!There is a delta(mu,nu) in the vertex



!	  enddo	!nu
	  enddo !mu
	enddo	!u0
	enddo	!t0
	enddo	!s0


	enddo	!momenta
	enddo
	enddo

	!$OMP END PARALLEL DO


	gX_6a(x0) = f_temp
	enddo	!x0



 	gX_6a(:) = 0.5d0*dimrep*C2_R*gX_6a(:)/((l**3)*1.0d0)

	
	END SUBROUTINE gX6a


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE gX6b(l,g_i,f1,f2,gX_6b)
	!
	! Subroutine to compute diagram gX_6b 
	!
	implicit none

	integer,intent(in) :: l,g_i,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: gX_6b
	

	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop_0
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,c0:c1,0:l) :: SGS_O5_S	
	complex(kind=8),dimension(4,4) :: V2mu	
	integer :: x0,i
	integer :: s0,t0,u0
	integer :: t0_a,t0_b,u0_a,u0_b
	integer :: mu, nu
	integer :: n1, n2, n3
	integer, dimension(3) :: p,q
	complex(kind=8) :: trace
	real(kind=8) :: multi
	integer :: f

	!Variables for the openMP paralelization
	integer :: chunk
	complex(kind=8) :: f_temp
	chunk = 1

	!initialization
	S_prop_0(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	SGS_O5_S(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	gX_6b(:) = dcmplx(0.0d0,0.0d0)
	p(:) = 0
	f=f1f2(f2,f1)
	f_temp = dcmplx(0.0d0,0.0d0)

	!build S_O5_SGS

	call Fermion_Propagator(l,p,S_prop_0)

	do x0 = c0,c1
	do s0=0,l
	do t0=0,l
	      SGS_O5_S(:,:,t0,x0,s0) = matmul(S_prop_0(:,:,t0,x0,f1),matmul(matmul(gX(:,:,g_i),S_prop_0(:,:,x0,1,f2)),matmul(Q5(:,:,f),S_prop_0(:,:,1,s0,f1))))
	enddo	
	enddo
	enddo


	do x0 = c0,c1
	f_temp = dcmplx(0.0d0,0.0d0)

!$OMP PARALLEL DO &
!$OMP DEFAULT(PRIVATE) SHARED(l,p,chunk,f1,f2,x0,SGS_O5_S,multiplicity,mom_deg) &
!$OMP REDUCTION(+:f_temp) &
!$OMP SCHEDULE(DYNAMIC,chunk) &
!$OMP COLLAPSE(1)

	do n1=0,l-1			
	do n2=n_limit(mom_deg,n1),l-1
	do n3=n_limit(mom_deg,n2),l-1

	q(1)=n1
	q(2)=n2
	q(3)=n3

	multi = 1.0d0*multiplicity(n1,n2,n3)

  	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)

  	call Gluon_propagator(l,q,D_prop)	

	trace = dcmplx(0.0d0,0.0d0)
	
	do s0=0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0=t0_a,t0_b
	do u0=u0_a,u0_b		!Due to the delta(x0+1,y0),delta(x0,y0),delta(x0-1,y0) in the vertex functions 

	  do mu=0,3	
!	  do nu=0,3

	  call V2(l,p,-p,q,-q,s0,t0,u0,u0,mu,V2mu)	!Remember that this guy has a delta(mu,nu) inside

	  trace = Trace2(SGS_O5_S(:,:,t0,x0,s0),V2mu(:,:))
	  f_temp = f_temp + multi*trace*D_prop(mu,mu,u0,u0)	!There is a delta(mu,nu) in the vertex


!	  enddo	!nu
	  enddo !mu
	enddo	!u0
	enddo	!t0
	enddo	!s0


	enddo	!momenta
	enddo
	enddo

	gX_6b(x0) = f_temp
	enddo	!x0

 	gX_6b(:) = 0.5d0*dimrep*C2_R*gX_6b(:)/((l**3)*1.0d0)

	
	END SUBROUTINE gX6b


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



	SUBROUTINE gX7(l,g_i,f1,f2,gX_7)
	!
	! Subroutine to compute diagram gX_7
	!
	implicit none

	integer,intent(in) :: l,g_i,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: gX_7
	

	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l) :: S_O5_S	
	complex(kind=8),dimension(4,4) :: V1mu, V1nu	
	integer :: x0,i
	integer :: s0,s0p,t0,t0p,u0,u0p
	integer :: t0_a,t0_b,t0p_a,t0p_b,u0_a,u0_b,u0p_a,u0p_b
	integer :: mu, nu
	integer :: n1, n2, n3
	integer, dimension(3) :: p,q
	complex(kind=8) :: trace
	real(kind=8) :: multi
	integer :: f


	!Variables for the openMP paralelization
	integer :: chunk
	complex(kind=8) :: f_temp
	chunk = 1



	!initialization
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_O5_S(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	gX_7(:) = dcmplx(0.0d0,0.0d0)
	p(:) = 0
	f=f1f2(f2,f1)
	f_temp = dcmplx(0.0d0,0.0d0)


	!build S_O5_S

	call Fermion_Propagator(l,p,S_prop_0)

	do s0=0,l
	do t0=0,l
	    S_O5_S(:,:,t0,s0) = matmul(S_prop_0(:,:,t0,1,f2),matmul(Q5(:,:,f),S_prop_0(:,:,1,s0,f1)))
	enddo	
	enddo


	do x0 = c0,c1

	f_temp = dcmplx(0.0d0,0.0d0)


!$OMP PARALLEL DO &
!$OMP DEFAULT(PRIVATE) SHARED(l,p,chunk,f1,f2,x0,g_i,S_O5_S,gX,multiplicity,mom_deg) &
!$OMP REDUCTION(+:f_temp) &
!$OMP SCHEDULE(DYNAMIC,chunk) &
!$OMP COLLAPSE(1)

	do n1=0,l-1			
	do n2=n_limit(mom_deg,n1),l-1
	do n3=n_limit(mom_deg,n2),l-1

	q(1)=n1
	q(2)=n2
	q(3)=n3

	multi = 1.0d0*multiplicity(n1,n2,n3)

  	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

  	call Gluon_propagator(l,q,D_prop)	
	call Fermion_Propagator(l,p+q,S_prop)


	trace = dcmplx(0.0d0,0.0d0)
	
	do s0=0,l
	call loop_limits(l,s0,t0p_a,t0p_b,u0_a,u0_b)
	do s0p=0,l
	call loop_limits(l,s0p,t0_a,t0_b,u0p_a,u0p_b)
	do t0=t0_a,t0_b
	do t0p=t0p_a,t0p_b		!Due to the delta(x0+1,y0),delta(x0,y0),delta(x0-1,y0) in the vertex functions 
	do u0=u0_a,u0_b		!Due to the delta(x0+1,y0),delta(x0,y0),delta(x0-1,y0) in the vertex functions 
	do u0p=u0p_a,u0p_b

	  do mu=0,3	
	  do nu=mu_low(mu),mu_up(mu)

	  call V1(l,p,-p-q,q,s0,t0p,u0,mu,V1mu)
	  call V1(l,p+q,-p,-q,s0p,t0,u0p,nu,V1nu)

	  trace = Trace6(S_O5_S(:,:,t0,s0),V1mu(:,:),S_prop(:,:,t0p,x0,f1),gX(:,:,g_i),S_prop(:,:,x0,s0p,f2),V1nu(:,:))
	  f_temp = f_temp + multi*trace*D_prop(nu,mu,u0p,u0)

	  enddo	!nu
	  enddo !mu
	  enddo	!u0p
	  enddo	!u0
	enddo	!t0p
	enddo	!t0
	enddo	!s0p
	enddo	!s0


	enddo	!momenta
	enddo
	enddo
!$OMP END PARALLEL DO

	gX_7(x0) = f_temp

	enddo	!x0


 	gX_7(:) = -0.5d0*dimrep*C2_R*gX_7(:)/((l**3)*1.0d0)

	
	END SUBROUTINE gX7


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE gXdm(l,g_i,f1,f2,gX_dm)
	!
	! Subroutine to compute diagram gX_dm
	!
	implicit none

	integer,intent(in) :: l, g_i,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: gX_dm
	
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop
	complex(kind=8),dimension(c0:c1) :: gX_dm_1, gX_dm_2
	integer, dimension(3) :: p_ext
	integer :: f,s0,x0

	f=f1f2(f2,f1)
	p_ext(:) = 0

	!initialization
	gX_dm(:) = dcmplx(0.0d0,0.0d0)
	gX_dm_1(:) = dcmplx(0.0d0,0.0d0)
	gX_dm_2(:) = dcmplx(0.0d0,0.0d0)

	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)


	call Fermion_Propagator(l,p_ext,S_prop)

	do x0 = c0,c1
	do s0 = 0,l

	gX_dm_1(x0) = gX_dm_1(x0)-dimrep*Trace5(Q5(:,:,f),S_prop(:,:,1,x0,f1),gX(:,:,g_i),S_prop(:,:,x0,s0,f2),S_prop(:,:,s0,1,f2))
	gX_dm_2(x0) = gX_dm_2(x0)-dimrep*Trace5(Q5(:,:,f),S_prop(:,:,1,s0,f1),S_prop(:,:,s0,x0,f1),gX(:,:,g_i),S_prop(:,:,x0,1,f2))

	enddo
	enddo

	gX_dm(:) = (gX_dm_1(:)+gX_dm_2(:))/2.0d0


	END SUBROUTINE gXdm


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE gXzf(l,g_i,f1,f2,gX_zf)
	!
	! Subroutine to compute diagram gX_zf
	!
	implicit none

	integer,intent(in) :: l, g_i,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: gX_zf
	
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop
	complex(kind=8),dimension(c0:c1) :: gX_zf_1, gX_zf_2
	integer, dimension(3) :: p_ext
	integer :: f,x0

	f=f1f2(f2,f1)
	p_ext(:) = 0

	!initialization
	gX_zf(:) = dcmplx(0.0d0,0.0d0)
	gX_zf_1(:) = dcmplx(0.0d0,0.0d0)
	gX_zf_2(:) = dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)


	call Fermion_Propagator(l,p_ext,S_prop)

	do x0 = c0,c1

	gX_zf_1(x0) = gX_zf_1(x0)+Trace5(Q5(:,:,f),S_prop(:,:,1,x0,f1),gX(:,:,g_i),S_prop(:,:,x0,0,f2),S_prop(:,:,0,1,f2))
	gX_zf_1(x0) = gX_zf_1(x0)+Trace5(Q5(:,:,f),S_prop(:,:,1,x0,f1),gX(:,:,g_i),S_prop(:,:,x0,l,f2),S_prop(:,:,l,1,f2))

	gX_zf_2(x0) = gX_zf_2(x0)+Trace5(Q5(:,:,f),S_prop(:,:,1,0,f1),S_prop(:,:,0,x0,f1),gX(:,:,g_i),S_prop(:,:,x0,1,f2))
	gX_zf_2(x0) = gX_zf_2(x0)+Trace5(Q5(:,:,f),S_prop(:,:,1,l,f1),S_prop(:,:,l,x0,f1),gX(:,:,g_i),S_prop(:,:,x0,1,f2))

	enddo

	gX_zf(:) = -dimrep*(gX_zf_1(:)+gX_zf_2(:))/2.0d0

	END SUBROUTINE gXzf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE gXds(l,g_i,f1,f2,gX_ds)
	!
	! Subroutine to compute diagram gX_ds
	!
	implicit none

	integer,intent(in) :: l, g_i,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: gX_ds
	
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop
	complex(kind=8),dimension(c0:c1) :: gX_ds_1, gX_ds_2
	integer, dimension(3) :: p_ext
	real(kind=8),dimension(3) ::  pe_plus, ptwit, phat
	complex(kind=8),dimension(4,4) :: Amat
	integer :: f,x0

	f=f1f2(f2,f1)
	p_ext(:) = 0

	!initialization
	gX_ds(:) = dcmplx(0.0d0,0.0d0)
	gX_ds_1(:) = dcmplx(0.0d0,0.0d0)
	gX_ds_2(:) = dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)


	pe_plus(:) = p_ext(:)*(2.0d0*pi)/(l*1.0d0) + theta/l
	ptwit(:) = sin(pe_plus(:))
	phat(:) = 2.0d0*sin(pe_plus(:)/2.0d0)

	Amat = imag*(ptwit(1)*g1+ptwit(2)*g2+ptwit(3)*g3)
	Amat = Amat + 0.5d0*II*( phat(1)**2 + phat(2)**2 + phat(3)**2 )



	call Fermion_Propagator(l,p_ext,S_prop)

	do x0 = c0,c1

	gX_ds_1(x0) = gX_ds_1(x0)+Trace6(Q5(:,:,f),S_prop(:,:,1,x0,f1),gX(:,:,g_i),S_prop(:,:,x0,0,f2),Amat,S_prop(:,:,0,1,f2))
	gX_ds_1(x0) = gX_ds_1(x0)+Trace6(Q5(:,:,f),S_prop(:,:,1,x0,f1),gX(:,:,g_i),S_prop(:,:,x0,l,f2),Amat,S_prop(:,:,l,1,f2))

	gX_ds_2(x0) = gX_ds_2(x0)+Trace6(Q5(:,:,f),S_prop(:,:,1,0,f1),Amat,S_prop(:,:,0,x0,f1),gX(:,:,g_i),S_prop(:,:,x0,1,f2))
	gX_ds_2(x0) = gX_ds_2(x0)+Trace6(Q5(:,:,f),S_prop(:,:,1,l,f1),Amat,S_prop(:,:,l,x0,f1),gX(:,:,g_i),S_prop(:,:,x0,1,f2))

	enddo

	gX_ds(:) = -dimrep*(gX_ds_1(:)+gX_ds_2(:))/2.0d0

	END SUBROUTINE gXds

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	END MODULE OneLoop_gX


