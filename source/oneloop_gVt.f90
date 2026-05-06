




	MODULE OneLoop_gVt
	!
	! Contains the set of routines for the calculaton of the 1-loop
	! diagrams contributing to the gVt correlation functions (including
	! the contribution of counterterms.
	!
	USE GammaMatrices
	USE Input
	USE Parameters
	USE Vertices
	USE Propagators
	USE TreeLevel_gVt
	USE omp_lib

	implicit none

	CONTAINS


	SUBROUTINE gVt1a(l,f1,f2,gVt_1a)
	!
	! Subroutine to compute diagram gVt_1a 
	!
	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: gVt_1a
	
	integer :: n1,n2,n3
	integer,dimension(3) :: qloop
	complex(kind=8) :: loop
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop!,K_prop
	complex(kind=8),dimension(c0:c1) :: gVt_0
	real(kind=8) :: multi	


	!initialization
	gVt_1a(:) = dcmplx(0.0d0,0.0d0)
	loop = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	gVt_0(:) = dcmplx(0.0d0,0.0d0)

	!Tree level part
	call gVt0(l,f1,f2,gVt_0)

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

	gVt_1a(:) = -0.5d0*C2_R*loop*gVt_0(:)/((l*1.0d0)**3)

	END SUBROUTINE gVt1a

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



	SUBROUTINE gVt2(l,f1,f2,gVt_2)
	!
	! Subroutine to compute diagram gVt_2
	!
	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: gVt_2
	
	integer :: x0,n1,n2,n3,i
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: loop_1, loop_2
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop!,K_prop	
	complex(kind=8) :: trace_1, trace_2

	real(kind=8) :: multi	
	integer :: f

	!initialization
	gVt_2(:) = dcmplx(0.0d0,0.0d0)
	loop_1(:) = dcmplx(0.0d0,0.0d0)
	loop_2(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	loop_1(:) = dcmplx(0.0d0,0.0d0)
	loop_2(:) = dcmplx(0.0d0,0.0d0)

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

	   trace_1 = 0.5d0*dimrep*Trace4(Q5(:,:,f),S_prop(:,:,1,x0+1,f1),Pplus,S_prop(:,:,x0,1,f2))
	   trace_2 = 0.5d0*dimrep*Trace4(Q5(:,:,f),S_prop(:,:,1,x0,f1),Pminus,S_prop(:,:,x0+1,1,f2))

	       loop_1(x0) = loop_1(x0) + D_prop(0,0,0,0)*trace_1*multi
	       loop_2(x0) = loop_2(x0) + D_prop(0,0,0,0)*trace_2*multi

	     enddo


	enddo
	enddo
	enddo

	  gVt_2(:) = C2_R*(loop_1(:) - loop_2(:))/((l*1.0d0)**3)


	END SUBROUTINE gVt2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




	SUBROUTINE gVt3a(l,f1,f2,gVt_3a)
	!
	! Subroutine to compute diagram gVt_3a
	!
	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: gVt_3a

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,c0:c1,0:l) :: S_G_1, S_G_2
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop, p_ext
	complex(kind=8),dimension(c0:c1) :: g_loop_1, g_loop_2
	integer :: n1, n2, n3, x0, t0, u0, s0, mu,i, x0_t
	integer :: t0_a, t0_b, u0_a, u0_b
	real(kind=8) :: multi
	complex(kind=8) :: trace_1, trace_2
	integer :: f

	gVt_3a(:) = dcmplx(0.0d0,0.0d0)
	g_loop_1(:) = dcmplx(0.0d0,0.0d0)
	g_loop_2(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_G_1(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_G_2(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	p_ext(:) = 0
	trace_1 = dcmplx(0.0d0,0.0d0)
	trace_2 = dcmplx(0.0d0,0.0d0)
	f=f1f2(f2,f1)

	!part outside the loop
	call Fermion_Propagator(l,p_ext,S_prop_0)

	do x0 = c0,c1
	do s0=0,l
		
	    S_G_1(:,:,x0,s0) = matmul(S_prop_0(:,:,1,x0+1,f1),matmul(Pplus(:,:),S_prop_0(:,:,x0,s0,f2)))
	    S_G_2(:,:,x0,s0) = matmul(S_prop_0(:,:,1,x0,f1),matmul(Pminus(:,:),S_prop_0(:,:,x0+1,s0,f2)))

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
	  trace_1 = Trace4(S_G_1(:,:,x0,s0),V1mu(:,:),S_prop(:,:,t0,1,f2),Q5(:,:,f))
	  trace_2 = Trace4(S_G_2(:,:,x0,s0),V1mu(:,:),S_prop(:,:,t0,1,f2),Q5(:,:,f))
	  g_loop_1(x0) = g_loop_1(x0) + D_prop(0,mu,0,u0)*trace_1*multi
	  g_loop_2(x0) = g_loop_2(x0) + D_prop(0,mu,0,u0)*trace_2*multi
	enddo

	
	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu


	enddo	!momenta
	enddo
	enddo

	gVt_3a(:) = -0.5*C2_R*dimrep*(g_loop_1(:)-g_loop_2(:))/(1.0d0*l**3)


	END SUBROUTINE gVt3a


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE gVt3b(l,f1,f2,gVt_3b)
	!
	! Subroutine to compute diagram gVt_3b
	!
	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: gVt_3b

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,0:l,c0:c1) :: S_G_1, S_G_2
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop, p_ext
	complex(kind=8),dimension(c0:c1) :: g_loop_1, g_loop_2
	integer :: n1, n2, n3, x0, t0, u0, s0, mu,i, x0_t
	integer :: t0_a, t0_b, u0_a, u0_b
	real(kind=8) :: multi
	complex(kind=8) :: trace_1, trace_2 
	integer :: f

	gVt_3b(:) = dcmplx(0.0d0,0.0d0)
	g_loop_1(:) = dcmplx(0.0d0,0.0d0)
	g_loop_2(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_G_1(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_G_2(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	p_ext(:) = 0
	trace_2 = dcmplx(0.0d0,0.0d0)
	trace_1 = dcmplx(0.0d0,0.0d0)
	f=f1f2(f2,f1)
	!part outside the loop
	call Fermion_Propagator(l,p_ext,S_prop_0)

	do t0 = 0,l
	do x0 = c0,c1		
	    S_G_1(:,:,t0,x0) = matmul(S_prop_0(:,:,t0,x0+1,f1),matmul(Pplus(:,:),S_prop_0(:,:,x0,1,f2)))
	    S_G_2(:,:,t0,x0) = matmul(S_prop_0(:,:,t0,x0,f1),matmul(Pminus(:,:),S_prop_0(:,:,x0+1,1,f2)))
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
	  trace_1 = Trace4(Q5(:,:,f),S_prop(:,:,1,s0,f1),V1mu(:,:),S_G_1(:,:,t0,x0))
	  trace_2 = Trace4(Q5(:,:,f),S_prop(:,:,1,s0,f1),V1mu(:,:),S_G_2(:,:,t0,x0))
	  g_loop_1(x0) = g_loop_1(x0) + D_prop(mu,0,u0,0)*trace_1*multi
	  g_loop_2(x0) = g_loop_2(x0) + D_prop(mu,0,u0,0)*trace_2*multi
	enddo

	
	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu


	enddo	!momenta
	enddo
	enddo

	  gVt_3b(:) = 0.5*C2_R*dimrep*(g_loop_1(:)-g_loop_2(:))/(1.0d0*l**3)

	END SUBROUTINE gVt3b



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



	SUBROUTINE gVt4a(l,f1,f2,gVt_4a)
	!
	! Subroutine to compute diagram gVt_4a
	!
	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: gVt_4a

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop,S_prop_0
	complex(kind=8),dimension(4,4,0:l,c0:c1) :: S_G_1, S_G_2
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop, p_ext
	complex(kind=8),dimension(c0:c1) :: g_loop_1, g_loop_2
	integer :: n1, n2, n3, x0, t0, u0, s0, mu, i
	integer :: t0_a, t0_b, u0_a, u0_b
	real(kind=8) :: multi 
	complex(kind=8) :: trace_1, trace_2
	integer :: f

	gVt_4a(:) = dcmplx(0.0d0,0.0d0)
	g_loop_1(:) = dcmplx(0.0d0,0.0d0)
	g_loop_2(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_G_1(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_G_2(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	p_ext(:) = 0
	trace_1 = dcmplx(0.0d0,0.0d0)
	trace_2 = dcmplx(0.0d0,0.0d0)
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
	      S_G_1(:,:,t0,x0) = matmul(S_prop(:,:,t0,x0+1,f1),matmul(Pplus(:,:),S_prop(:,:,x0,1,f2)))
	      S_G_2(:,:,t0,x0) = matmul(S_prop(:,:,t0,x0,f1),matmul(Pminus(:,:),S_prop(:,:,x0+1,1,f2)))
	   enddo
	   enddo

	do mu=mu_low(0),mu_up(0)

	do s0=0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0=t0_a,t0_b
	do u0=u0_a,u0_b

	call V1(l,p_ext,-p_ext-qloop,qloop,s0,t0,u0,mu,V1mu)

	do x0=c0,c1
	  trace_1 = Trace4(Q5(:,:,f),S_prop_0(:,:,1,s0,f1),V1mu(:,:),S_G_1(:,:,t0,x0))
	  trace_2 = Trace4(Q5(:,:,f),S_prop_0(:,:,1,s0,f1),V1mu(:,:),S_G_2(:,:,t0,x0))
	  g_loop_1(x0) = g_loop_1(x0) + D_prop(0,mu,0,u0)*trace_1*multi
	  g_loop_2(x0) = g_loop_2(x0) + D_prop(0,mu,0,u0)*trace_2*multi
	enddo	!x0
	

	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu


	enddo	!momentum
	enddo
	enddo

	gVt_4a(:) = -0.5*C2_R*dimrep*(g_loop_1(:)-g_loop_2(:))/(1.0d0*l**3)

	END SUBROUTINE gVt4a


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE gVt4b(l,f1,f2,gVt_4b)
	!
	! Subroutine to compute diagram gVt_4b
	!
	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: gVt_4b

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,c0:c1,0:l) :: S_G_1, S_G_2
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop, p_ext
	complex(kind=8),dimension(c0:c1) :: g_loop_1, g_loop_2
	integer :: n1, n2, n3, x0, t0, u0, s0, mu, i
	integer :: t0_a, t0_b, u0_a, u0_b
	real(kind=8) :: multi
	complex(kind=8) :: trace_1, trace_2
	integer :: f

	gVt_4b(:) = dcmplx(0.0d0,0.0d0)
	g_loop_1(:) = dcmplx(0.0d0,0.0d0)
	g_loop_2(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_G_1(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_G_2(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	p_ext(:) = 0
	trace_1 = dcmplx(0.0d0,0.0d0)
	trace_2 = dcmplx(0.0d0,0.0d0)
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
	      S_G_1(:,:,x0,s0) = matmul(S_prop(:,:,1,x0+1,f1),matmul(Pplus(:,:),S_prop(:,:,x0,s0,f2)))
	      S_G_2(:,:,x0,s0) = matmul(S_prop(:,:,1,x0,f1),matmul(Pminus(:,:),S_prop(:,:,x0+1,s0,f2)))
	   enddo
	   enddo


	do mu=mu_low(0),mu_up(0)

	do s0=0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0=t0_a,t0_b
	do u0=u0_a,u0_b

	call V1(l,p_ext+qloop,-p_ext,-qloop,s0,t0,u0,mu,V1mu)

	do x0 = c0,c1
	  trace_1 = Trace4(Q5(:,:,f),S_G_1(:,:,x0,s0),V1mu(:,:),S_prop_0(:,:,t0,1,f2))
	  trace_2 = Trace4(Q5(:,:,f),S_G_2(:,:,x0,s0),V1mu(:,:),S_prop_0(:,:,t0,1,f2))
	  g_loop_1(x0) = g_loop_1(x0) + D_prop(mu,0,u0,0)*trace_1*multi
	  g_loop_2(x0) = g_loop_2(x0) + D_prop(mu,0,u0,0)*trace_2*multi
	enddo	!x0
	

	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu



	enddo	!momentum
	enddo
	enddo

	gVt_4b(:) = 0.5*C2_R*dimrep*(g_loop_1(:)-g_loop_2(:))/(1.0d0*l**3)

	END SUBROUTINE gVt4b
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




	SUBROUTINE gVt5a(l,f1,f2,gVt_5a)
	!
	! Subroutine to compute diagram gVt_5a
	!
	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: gVt_5a
	
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,c0:c1,0:l) :: S_O5_SGS_1, S_O5_SGS_2	
	complex(kind=8),dimension(4,4) :: V1mu, V1nu	
	integer :: x0,i
	integer :: s0,s0p,t0,t0p,u0,u0p
	integer :: t0_a,t0_b,t0p_a,t0p_b,u0_a,u0_b,u0p_a,u0p_b
	integer :: mu, nu
	integer :: n1, n2, n3
	integer, dimension(3) :: p,q
	complex(kind=8) :: trace_1, trace_2
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
	S_O5_SGS_1(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_O5_SGS_2(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace_1 = dcmplx(0.0d0,0.0d0)
	trace_2 = dcmplx(0.0d0,0.0d0)
	gVt_5a(:) = dcmplx(0.0d0,0.0d0)
	p(:) = 0
	f=f1f2(f2,f1)

	!build S_O5_SGS

	call Fermion_Propagator(l,p,S_prop_0)

	do x0 = c0,c1
	do s0=0,l
	do t0=0,l
	    S_O5_SGS_1(:,:,t0,x0,s0) = matmul(S_prop_0(:,:,t0,1,f2),matmul(matmul(Q5(:,:,f),S_prop_0(:,:,1,x0+1,f1)),matmul(Pplus(:,:),S_prop_0(:,:,x0,s0,f2))))
	    S_O5_SGS_2(:,:,t0,x0,s0) = matmul(S_prop_0(:,:,t0,1,f2),matmul(matmul(Q5(:,:,f),S_prop_0(:,:,1,x0,f1)),matmul(Pminus(:,:),S_prop_0(:,:,x0+1,s0,f2))))
	enddo	
	enddo
	enddo


	do x0 = c0,c1

	f_temp = dcmplx(0.0d0,0.0d0)

!$OMP PARALLEL DO &
!$OMP DEFAULT(PRIVATE) SHARED(l,p,chunk,f1,f2,x0,S_O5_SGS_1,S_O5_SGS_2,multiplicity,mom_deg) &
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

	trace_1 = dcmplx(0.0d0,0.0d0)
	trace_2 = dcmplx(0.0d0,0.0d0)
	
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

	  trace_1 = Trace4(S_O5_SGS_1(:,:,t0,x0,s0),V1mu(:,:),S_prop(:,:,t0p,s0p,f2),V1nu(:,:))
	  trace_2 = Trace4(S_O5_SGS_2(:,:,t0,x0,s0),V1mu(:,:),S_prop(:,:,t0p,s0p,f2),V1nu(:,:))

	   f_temp = f_temp + multi*(trace_1-trace_2)*D_prop(nu,mu,u0p,u0)



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

	gVt_5a(x0) = f_temp
	
	enddo	!x0
 	gVt_5a(:) = -0.5d0*dimrep*C2_R*gVt_5a(:)/((l**3)*1.0d0)

	
	END SUBROUTINE gVt5a


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



	SUBROUTINE gVt5b(l,f1,f2,gVt_5b)
	!
	! Subroutine to compute diagram gVt_5b
	!
	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: gVt_5b
	

	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,c0:c1,0:l) :: SGS_O5_S_1, SGS_O5_S_2	
	complex(kind=8),dimension(4,4) :: V1mu, V1nu	
	integer :: x0,i
	integer :: s0,s0p,t0,t0p,u0,u0p
	integer :: t0_a,t0_b,t0p_a,t0p_b,u0_a,u0_b,u0p_a,u0p_b
	integer :: mu, nu
	integer :: n1, n2, n3
	integer, dimension(3) :: p,q
	complex(kind=8) :: trace_1, trace_2
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
	SGS_O5_S_1(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	SGS_O5_S_2(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace_1 = dcmplx(0.0d0,0.0d0)
	trace_2 = dcmplx(0.0d0,0.0d0)
	gVt_5b(:) = dcmplx(0.0d0,0.0d0)
	p(:) = 0
	f=f1f2(f2,f1)

	!build S_O5_SGS

	call Fermion_Propagator(l,p,S_prop_0)

	do x0 = c0,c1
	do s0=0,l
	do t0=0,l
	    SGS_O5_S_1(:,:,t0,x0,s0) = matmul(S_prop_0(:,:,t0,x0+1,f1),matmul(matmul(Pplus(:,:),S_prop_0(:,:,x0,1,f2)),matmul(Q5(:,:,f),S_prop_0(:,:,1,s0,f1))))
	    SGS_O5_S_2(:,:,t0,x0,s0) = matmul(S_prop_0(:,:,t0,x0,f1),matmul(matmul(Pminus(:,:),S_prop_0(:,:,x0+1,1,f2)),matmul(Q5(:,:,f),S_prop_0(:,:,1,s0,f1))))
	enddo	
	enddo
	enddo



	do x0 = c0,c1


	f_temp = dcmplx(0.0d0,0.0d0)

!$OMP PARALLEL DO &
!$OMP DEFAULT(PRIVATE) SHARED(l,p,chunk,f1,f2,x0,SGS_O5_S_1,SGS_O5_S_2,multiplicity,mom_deg) &
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

	trace_1 = dcmplx(0.0d0,0.0d0)
	trace_2 = dcmplx(0.0d0,0.0d0)
	
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

	  trace_1 = Trace4(SGS_O5_S_1(:,:,t0,x0,s0),V1mu(:,:),S_prop(:,:,t0p,s0p,f1),V1nu(:,:))
	  trace_2 = Trace4(SGS_O5_S_2(:,:,t0,x0,s0),V1mu(:,:),S_prop(:,:,t0p,s0p,f1),V1nu(:,:))
	  f_temp = f_temp + multi*(trace_1-trace_2)*D_prop(nu,mu,u0p,u0)


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

	gVt_5b(x0) = f_temp
	
	enddo	!x0

 	gVt_5b(:) = -0.5d0*dimrep*C2_R*gVt_5b(:)/((l**3)*1.0d0)

	
	END SUBROUTINE gVt5b


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



	SUBROUTINE gVt6a(l,f1,f2,gVt_6a)
	!
	! Subroutine to compute diagram gVt_6a
	!
	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: gVt_6a
	

	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop_0
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,c0:c1,0:l) :: S_O5_SGS_1, S_O5_SGS_2 	
	complex(kind=8),dimension(4,4) :: V2mu	
	integer :: x0,i
	integer :: s0,t0,u0
	integer :: t0_a,t0_b,u0_a,u0_b
	integer :: mu, nu
	integer :: n1, n2, n3
	integer, dimension(3) :: p,q
	complex(kind=8) :: trace_1, trace_2
	real(kind=8) :: multi
	integer :: f


	!Variables for the openMP paralelization
	integer :: chunk
	complex(kind=8) :: f_temp
	chunk = 1
	f_temp = dcmplx(0.0d0,0.0d0)


	!initialization
	S_prop_0(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_O5_SGS_1(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_O5_SGS_2(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace_1 = dcmplx(0.0d0,0.0d0)
	trace_2 = dcmplx(0.0d0,0.0d0)
	gVt_6a(:) = dcmplx(0.0d0,0.0d0)
	p(:) = 0
	f=f1f2(f2,f1)

	!build S_O5_SGS

	call Fermion_Propagator(l,p,S_prop_0)

	do x0 = c0,c1
	do s0=0,l
	do t0=0,l
	    S_O5_SGS_1(:,:,t0,x0,s0) = matmul(S_prop_0(:,:,t0,1,f2),matmul(matmul(Q5(:,:,f),S_prop_0(:,:,1,x0+1,f1)),matmul(Pplus(:,:),S_prop_0(:,:,x0,s0,f2))))
	    S_O5_SGS_2(:,:,t0,x0,s0) = matmul(S_prop_0(:,:,t0,1,f2),matmul(matmul(Q5(:,:,f),S_prop_0(:,:,1,x0,f1)),matmul(Pminus(:,:),S_prop_0(:,:,x0+1,s0,f2))))

	enddo	
	enddo
	enddo


	do x0 = c0,c1
	f_temp = dcmplx(0.0d0,0.0d0)
!$OMP PARALLEL DO &
!$OMP DEFAULT(PRIVATE) SHARED(l,p,chunk,f1,f2,x0,S_O5_SGS_1,S_O5_SGS_2,multiplicity,mom_deg) &
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

	trace_1 = dcmplx(0.0d0,0.0d0)
	trace_2 = dcmplx(0.0d0,0.0d0)
	
	do s0=0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0=t0_a,t0_b
	do u0=u0_a,u0_b		!Due to the delta(x0+1,y0),delta(x0,y0),delta(x0-1,y0) in the vertex functions 

	  do mu=0,3	
!	  do nu=0,3

	  call V2(l,p,-p,q,-q,s0,t0,u0,u0,mu,V2mu)	!Remember that this guy has a delta(mu,nu) inside

	  trace_1 = Trace2(S_O5_SGS_1(:,:,t0,x0,s0),V2mu(:,:))
	  trace_2 = Trace2(S_O5_SGS_2(:,:,t0,x0,s0),V2mu(:,:))

	  f_temp = f_temp + multi*(trace_1-trace_2)*D_prop(mu,mu,u0,u0)	!There is a delta(mu,nu) in the vertex


!	  enddo	!nu
	  enddo !mu
	enddo	!u0
	enddo	!t0
	enddo	!s0


	enddo	!momenta
	enddo
	enddo

!$OMP END PARALLEL DO

	gVt_6a(x0) = f_temp
	enddo	!x0


 	gVt_6a(:) = 0.5d0*dimrep*C2_R*gVt_6a(:)/((l**3)*1.0d0)

	
	END SUBROUTINE gVt6a


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE gVt6b(l,f1,f2,gVt_6b)
	!
	! Subroutine to compute diagram gVt_6b
	!
	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: gVt_6b
	

	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop_0
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,c0:c1,0:l) :: SGS_O5_S_1, SGS_O5_S_2 	
	complex(kind=8),dimension(4,4) :: V2mu	
	integer :: x0,i
	integer :: s0,t0,u0
	integer :: t0_a,t0_b,u0_a,u0_b
	integer :: mu, nu
	integer :: n1, n2, n3
	integer, dimension(3) :: p,q
	complex(kind=8) :: trace_1 ,trace_2
	real(kind=8) :: multi
	integer :: f


	!Variables for the openMP paralelization
	integer :: chunk
	complex(kind=8) :: f_temp
	chunk = 1
	f_temp = dcmplx(0.0d0,0.0d0)


	!initialization
	S_prop_0(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	SGS_O5_S_1(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	SGS_O5_S_2(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace_1 = dcmplx(0.0d0,0.0d0)
	trace_2 = dcmplx(0.0d0,0.0d0)
	gVt_6b(:) = dcmplx(0.0d0,0.0d0)
	p(:) = 0
	f=f1f2(f2,f1)

	!build S_O5_SGS

	call Fermion_Propagator(l,p,S_prop_0)

	do x0 = c0,c1
	do s0=0,l
	do t0=0,l
	      SGS_O5_S_1(:,:,t0,x0,s0) = matmul(S_prop_0(:,:,t0,x0+1,f1),matmul(matmul(Pplus(:,:),S_prop_0(:,:,x0,1,f2)),matmul(Q5(:,:,f),S_prop_0(:,:,1,s0,f1))))
	      SGS_O5_S_2(:,:,t0,x0,s0) = matmul(S_prop_0(:,:,t0,x0,f1),matmul(matmul(Pminus(:,:),S_prop_0(:,:,x0+1,1,f2)),matmul(Q5(:,:,f),S_prop_0(:,:,1,s0,f1))))
	enddo	
	enddo
	enddo


	do x0 = c0,c1
	f_temp = dcmplx(0.0d0,0.0d0)
!$OMP PARALLEL DO &
!$OMP DEFAULT(PRIVATE) SHARED(l,p,chunk,f1,f2,x0,SGS_O5_S_1,SGS_O5_S_2,multiplicity,mom_deg) &
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

	trace_1 = dcmplx(0.0d0,0.0d0)
	trace_2 = dcmplx(0.0d0,0.0d0)
	
	do s0=0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0=t0_a,t0_b
	do u0=u0_a,u0_b		!Due to the delta(x0+1,y0),delta(x0,y0),delta(x0-1,y0) in the vertex functions 

	  do mu=0,3	
!	  do nu=0,3

	  call V2(l,p,-p,q,-q,s0,t0,u0,u0,mu,V2mu)	!Remember that this guy has a delta(mu,nu) inside

	  trace_1 = Trace2(SGS_O5_S_1(:,:,t0,x0,s0),V2mu(:,:))
	  trace_2 = Trace2(SGS_O5_S_2(:,:,t0,x0,s0),V2mu(:,:))
	  f_temp = f_temp + multi*(trace_1 - trace_2)*D_prop(mu,mu,u0,u0)	


!	  enddo	!nu
	  enddo !mu
	enddo	!u0
	enddo	!t0
	enddo	!s0


	enddo	!momenta
	enddo
	enddo

!$OMP END PARALLEL DO

	gVt_6b(x0) = f_temp
	enddo	!x0
 	gVt_6b(:) = 0.5d0*dimrep*C2_R*gVt_6b(:)/((l**3)*1.0d0)

	
	END SUBROUTINE gVt6b


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




	SUBROUTINE gVt7(l,f1,f2,gVt_7)
	!
	! Subroutine to compute diagram gVt_7
	!
	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: gVt_7
	

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
	complex(kind=8) :: trace_1, trace_2
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
	S_O5_S(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace_1 = dcmplx(0.0d0,0.0d0)
	trace_2 = dcmplx(0.0d0,0.0d0)
	gVt_7(:) = dcmplx(0.0d0,0.0d0)
	p(:) = 0
	f=f1f2(f2,f1)

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
!$OMP DEFAULT(PRIVATE) SHARED(l,p,chunk,f1,f2,x0,S_O5_S,Pplus,Pminus,multiplicity,mom_deg) &
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


	trace_1 = dcmplx(0.0d0,0.0d0)
	trace_2 = dcmplx(0.0d0,0.0d0)
	
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

	  trace_1 = Trace6(S_O5_S(:,:,t0,s0),V1mu(:,:),S_prop(:,:,t0p,x0+1,f1),Pplus(:,:),S_prop(:,:,x0,s0p,f2),V1nu(:,:))
	  trace_2 = Trace6(S_O5_S(:,:,t0,s0),V1mu(:,:),S_prop(:,:,t0p,x0,f1),Pminus(:,:),S_prop(:,:,x0+1,s0p,f2),V1nu(:,:))
	  f_temp = f_temp + multi*(trace_1-trace_2)*D_prop(nu,mu,u0p,u0)


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

	gVt_7(x0) = f_temp
	enddo	!x0

 	gVt_7(:) = -0.5d0*dimrep*C2_R*gvt_7(:)/((l**3)*1.0d0)

	
	END SUBROUTINE gVt7

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE gVt8(l,f1,f2,gVt_8)
	!
	! Subroutine to compute diagram gVt_8
	!
	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: gVt_8
	
	integer :: n1,n2,n3,x0
	integer,dimension(3) :: qloop
	complex(kind=8) :: loop
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop!,K_prop
	complex(kind=8),dimension(c0:c1) :: gVt_0
	real(kind=8) :: multi	


	!initialization
	gVt_8(:) = dcmplx(0.0d0,0.0d0)
	loop = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	gVt_0(:) = dcmplx(0.0d0,0.0d0)

	!Tree level part
	call gVt0(l,f1,f2,gVt_0)

	do x0 = c0,c1

	loop = dcmplx(0.0d0,0.0d0)

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
	   loop = loop + D_prop(0,0,x0,x0)*multi

	enddo
	enddo
	enddo

	gVt_8(x0) = -0.5d0*C2_R*loop*gVt_0(x0)/((l*1.0d0)**3)

	enddo


	END SUBROUTINE gVt8

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




	SUBROUTINE gVt9a(l,f1,f2,gVt_9a)
	!
	! Subroutine to compute diagram gVt_9a
	!
	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: gVt_9a

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,c0:c1) :: S_G_1, S_G_2
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop, p_ext
	complex(kind=8),dimension(c0:c1) :: g_loop_1, g_loop_2
	integer :: n1, n2, n3, x0
	real(kind=8) :: multi
	complex(kind=8) :: trace_1, trace_2
	integer :: f

	gVt_9a(:) = dcmplx(0.0d0,0.0d0)
	g_loop_1(:) = dcmplx(0.0d0,0.0d0)
	g_loop_2(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_G_1(:,:,:)=dcmplx(0.0d0,0.0d0)
	S_G_2(:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	p_ext(:) = 0
	trace_1 = dcmplx(0.0d0,0.0d0)
	trace_2 = dcmplx(0.0d0,0.0d0)
	f=f1f2(f2,f1)

	!part outside the loop
	call Fermion_Propagator(l,p_ext,S_prop_0)

	do x0 = c0,c1
		
	    S_G_1(:,:,x0) = matmul(Q5(:,:,f),matmul(S_prop_0(:,:,1,x0+1,f1),Pplus(:,:)))
	    S_G_2(:,:,x0) = matmul(Q5(:,:,f),matmul(S_prop_0(:,:,1,x0,f1),Pminus(:,:)))

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
	
	do x0 = c0,c1
	  trace_1 = Trace2(S_G_1(:,:,x0),S_prop(:,:,x0,1,f2))
	  trace_2 = Trace2(S_G_2(:,:,x0),S_prop(:,:,x0+1,1,f2))
	  g_loop_1(x0) = g_loop_1(x0) + D_prop(0,0,0,x0)*trace_1*multi
	  g_loop_2(x0) = g_loop_2(x0) + D_prop(0,0,0,x0)*trace_2*multi
	enddo

	enddo	!momenta
	enddo
	enddo

	gVt_9a(:) = -0.5*C2_R*dimrep*(g_loop_1(:)+g_loop_2(:))/(1.0d0*l**3)


	END SUBROUTINE gVt9a


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE gVt9b(l,f1,f2,gVt_9b)
	!
	! Subroutine to compute diagram gVt_9b
	!
	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: gVt_9b

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,c0:c1) :: S_G_1, S_G_2
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop, p_ext
	complex(kind=8),dimension(c0:c1) :: g_loop_1, g_loop_2
	integer :: n1, n2, n3, x0
	real(kind=8) :: multi
	complex(kind=8) :: trace_1, trace_2 
	integer :: f

	gVt_9b(:) = dcmplx(0.0d0,0.0d0)
	g_loop_1(:) = dcmplx(0.0d0,0.0d0)
	g_loop_2(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_G_1(:,:,:)=dcmplx(0.0d0,0.0d0)
	S_G_2(:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	p_ext(:) = 0
	trace_2 = dcmplx(0.0d0,0.0d0)
	trace_1 = dcmplx(0.0d0,0.0d0)
	f=f1f2(f2,f1)
	!part outside the loop
	call Fermion_Propagator(l,p_ext,S_prop_0)

	do x0 = c0,c1		
	    S_G_1(:,:,x0) = matmul(matmul(Pplus(:,:),S_prop_0(:,:,x0,1,f2)),Q5(:,:,f))
	    S_G_2(:,:,x0) = matmul(matmul(Pminus(:,:),S_prop_0(:,:,x0+1,1,f2)),Q5(:,:,f))
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

	do x0=c0,c1
	  trace_1 = Trace2(S_prop(:,:,1,x0+1,f1),S_G_1(:,:,x0))
	  trace_2 = Trace2(S_prop(:,:,1,x0,f1),S_G_2(:,:,x0))
	  g_loop_1(x0) = g_loop_1(x0) + D_prop(0,0,x0,0)*trace_1*multi
	  g_loop_2(x0) = g_loop_2(x0) + D_prop(0,0,x0,0)*trace_2*multi
	enddo

	enddo	!momenta
	enddo
	enddo

	gVt_9b(:) = 0.5*C2_R*dimrep*(g_loop_1(:)+g_loop_2(:))/(1.0d0*l**3)

	END SUBROUTINE gVt9b



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




	SUBROUTINE gVt10a(l,f1,f2,gVt_10a)
	!
	! Subroutine to compute diagram gVt_10a
	!
	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: gVt_10a

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,0:l,c0:c1) :: SQSG_1, SQSG_2
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop, p_ext
	complex(kind=8),dimension(c0:c1) :: g_loop_1, g_loop_2
	integer :: n1, n2, n3, x0, t0, u0, s0, mu,i, x0_t
	integer :: t0_a, t0_b, u0_a, u0_b
	real(kind=8) :: multi
	complex(kind=8) :: trace_1, trace_2
	integer :: f

	gVt_10a(:) = dcmplx(0.0d0,0.0d0)
	g_loop_1(:) = dcmplx(0.0d0,0.0d0)
	g_loop_2(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	SQSG_1(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	SQSG_2(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	p_ext(:) = 0
	trace_1 = dcmplx(0.0d0,0.0d0)
	trace_2 = dcmplx(0.0d0,0.0d0)
	f=f1f2(f2,f1)

	!part outside the loop
	call Fermion_Propagator(l,p_ext,S_prop_0)

	do t0=0,l
	do x0 = c0,c1
		
	    SQSG_1(:,:,t0,x0) = matmul(S_prop_0(:,:,t0,1,f2),matmul(Q5(:,:,f),matmul(S_prop_0(:,:,1,x0+1,f1),Pplus(:,:))))
	    SQSG_2(:,:,t0,x0) = matmul(S_prop_0(:,:,t0,1,f2),matmul(Q5(:,:,f),matmul(S_prop_0(:,:,1,x0,f1),Pminus(:,:))))

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

	call V1(l,p_ext+qloop,-p_ext,-qloop,s0,t0,u0,mu,V1mu)
	
	do x0 = c0,c1
	  trace_1 = Trace3(SQSG_1(:,:,t0,x0),S_prop(:,:,x0,s0,f2),V1mu(:,:))
	  trace_2 = Trace3(SQSG_2(:,:,t0,x0),S_prop(:,:,x0+1,s0,f2),V1mu(:,:))
	  g_loop_1(x0) = g_loop_1(x0) + D_prop(mu,0,u0,x0)*trace_1*multi
	  g_loop_2(x0) = g_loop_2(x0) + D_prop(mu,0,u0,x0)*trace_2*multi
	enddo

	
	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu


	enddo	!momenta
	enddo
	enddo

	gVt_10a(:) = -0.5*C2_R*dimrep*(g_loop_1(:)+g_loop_2(:))/(1.0d0*l**3)


	END SUBROUTINE gVt10a

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE gVt10b(l,f1,f2,gVt_10b)
	!
	! Subroutine to compute diagram gVt_10b
	!
	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: gVt_10b

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,c0:c1,0:l) :: GSQS_1, GSQS_2
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop, p_ext
	complex(kind=8),dimension(c0:c1) :: g_loop_1, g_loop_2
	integer :: n1, n2, n3, x0, t0, u0, s0, mu,i, x0_t
	integer :: t0_a, t0_b, u0_a, u0_b
	real(kind=8) :: multi
	complex(kind=8) :: trace_1, trace_2 
	integer :: f

	gVt_10b(:) = dcmplx(0.0d0,0.0d0)
	g_loop_1(:) = dcmplx(0.0d0,0.0d0)
	g_loop_2(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	GSQS_1(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	GSQS_2(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	p_ext(:) = 0
	trace_2 = dcmplx(0.0d0,0.0d0)
	trace_1 = dcmplx(0.0d0,0.0d0)
	f=f1f2(f2,f1)
	!part outside the loop
	call Fermion_Propagator(l,p_ext,S_prop_0)


	do x0 = c0,c1		
	do s0 = 0,l
	    GSQS_1(:,:,x0,s0) = matmul(Pplus(:,:),matmul(S_prop_0(:,:,x0,1,f2),matmul(Q5(:,:,f),S_prop_0(:,:,1,s0,f1))))
	    GSQS_2(:,:,x0,s0) = matmul(Pminus(:,:),matmul(S_prop_0(:,:,x0+1,1,f2),matmul(Q5(:,:,f),S_prop_0(:,:,1,s0,f1))))
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

	call V1(l,p_ext,-p_ext-qloop,qloop,s0,t0,u0,mu,V1mu)
	
	do x0=c0,c1
	  trace_1 = Trace3(GSQS_1(:,:,x0,s0),V1mu(:,:),S_prop(:,:,t0,x0+1,f1))
	  trace_2 = Trace3(GSQS_2(:,:,x0,s0),V1mu(:,:),S_prop(:,:,t0,x0,f1))
	  g_loop_1(x0) = g_loop_1(x0) + D_prop(0,mu,x0,u0)*trace_1*multi
	  g_loop_2(x0) = g_loop_2(x0) + D_prop(0,mu,x0,u0)*trace_2*multi
	enddo

	
	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu


	enddo	!momenta
	enddo
	enddo

	  gVt_10b(:) = -0.5*C2_R*dimrep*(g_loop_1(:)+g_loop_2(:))/(1.0d0*l**3)

	END SUBROUTINE gVt10b



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE gVtdm(l,f1,f2,gVt_dm)
	!
	! Subroutine to compute diagram gVt_dm
	!
	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: gVt_dm
	
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop
	complex(kind=8),dimension(c0:c1) :: gVt_dm_1, gVt_dm_2
	integer, dimension(3) :: p_ext
	integer :: f,s0,x0

	f=f1f2(f2,f1)
	p_ext(:) = 0

	!initialization
	gVt_dm(:) = dcmplx(0.0d0,0.0d0)
	gVt_dm_1(:) = dcmplx(0.0d0,0.0d0)
	gVt_dm_2(:) = dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)


	call Fermion_Propagator(l,p_ext,S_prop)

	do x0 = c0,c1
	do s0 = 0,l

	gVt_dm_1(x0) = gVt_dm_1(x0)-dimrep*Trace5(Q5(:,:,f),S_prop(:,:,1,x0+1,f1),Pplus(:,:),S_prop(:,:,x0,s0,f2),S_prop(:,:,s0,1,f2))
	gVt_dm_1(x0) = gVt_dm_1(x0)-dimrep*Trace5(Q5(:,:,f),S_prop(:,:,1,s0,f1),S_prop(:,:,s0,x0+1,f1),Pplus(:,:),S_prop(:,:,x0,1,f2))


	gVt_dm_2(x0) = gVt_dm_2(x0)-dimrep*Trace5(Q5(:,:,f),S_prop(:,:,1,x0,f1),Pminus(:,:),S_prop(:,:,x0+1,s0,f2),S_prop(:,:,s0,1,f2))
	gVt_dm_2(x0) = gVt_dm_2(x0)-dimrep*Trace5(Q5(:,:,f),S_prop(:,:,1,s0,f1),S_prop(:,:,s0,x0,f1),Pminus(:,:),S_prop(:,:,x0+1,1,f2))

	enddo
	enddo

	gVt_dm(:) = (gVt_dm_1(:) - gVt_dm_2(:))/2.0d0

	END SUBROUTINE gVtdm


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE gVtzf(l,f1,f2,gVt_zf)
	!
	! Subroutine to compute diagram gVt_zf
	!
	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: gVt_zf
	
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop
	complex(kind=8),dimension(c0:c1) :: gVt1_zf_1, gVt1_zf_2
	complex(kind=8),dimension(c0:c1) :: gVt2_zf_1, gVt2_zf_2
	integer, dimension(3) :: p_ext
	integer :: f,x0

	f=f1f2(f2,f1)
	p_ext(:) = 0

	!initialization
	gVt_zf(:) = dcmplx(0.0d0,0.0d0)
	gVt1_zf_1(:) = dcmplx(0.0d0,0.0d0)
	gVt1_zf_2(:) = dcmplx(0.0d0,0.0d0)
	gVt2_zf_1(:) = dcmplx(0.0d0,0.0d0)
	gVt2_zf_2(:) = dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)


	call Fermion_Propagator(l,p_ext,S_prop)

	do x0 = c0,c1

	gVt1_zf_1(x0) = gVt1_zf_1(x0)+Trace5(Q5(:,:,f),S_prop(:,:,1,x0+1,f1),Pplus(:,:),S_prop(:,:,x0,0,f2),S_prop(:,:,0,1,f2))
	gVt1_zf_1(x0) = gVt1_zf_1(x0)+Trace5(Q5(:,:,f),S_prop(:,:,1,x0+1,f1),Pplus(:,:),S_prop(:,:,x0,l,f2),S_prop(:,:,l,1,f2))

	gVt1_zf_2(x0) = gVt1_zf_2(x0)+Trace5(Q5(:,:,f),S_prop(:,:,1,0,f1),S_prop(:,:,0,x0+1,f1),Pplus(:,:),S_prop(:,:,x0,1,f2))
	gVt1_zf_2(x0) = gVt1_zf_2(x0)+Trace5(Q5(:,:,f),S_prop(:,:,1,l,f1),S_prop(:,:,l,x0+1,f1),Pplus(:,:),S_prop(:,:,x0,1,f2))


	gVt2_zf_1(x0) = gVt2_zf_1(x0)+Trace5(Q5(:,:,f),S_prop(:,:,1,x0,f1),Pminus(:,:),S_prop(:,:,x0+1,0,f2),S_prop(:,:,0,1,f2))
	gVt2_zf_1(x0) = gVt2_zf_1(x0)+Trace5(Q5(:,:,f),S_prop(:,:,1,x0,f1),Pminus(:,:),S_prop(:,:,x0+1,l,f2),S_prop(:,:,l,1,f2))


	gVt2_zf_2(x0) = gVt2_zf_2(x0)+Trace5(Q5(:,:,f),S_prop(:,:,1,0,f1),S_prop(:,:,0,x0,f1),Pminus(:,:),S_prop(:,:,x0+1,1,f2))
	gVt2_zf_2(x0) = gVt2_zf_2(x0)+Trace5(Q5(:,:,f),S_prop(:,:,1,l,f1),S_prop(:,:,l,x0,f1),Pminus(:,:),S_prop(:,:,x0+1,1,f2))


	enddo


	gVt_zf(:) = -dimrep*(gVt1_zf_1(:)+gVt1_zf_2(:)-gVt2_zf_1(:)-gVt2_zf_2(:))/2.0d0

	END SUBROUTINE gVtzf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE gVtds(l,f1,f2,gVt_ds)
	!
	! Subroutine to compute diagram gVt_ds
	!
	implicit none

	integer,intent(in) :: l,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: gVt_ds
	
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop
	complex(kind=8),dimension(c0:c1) :: gVt1_ds_1, gVt1_ds_2
	complex(kind=8),dimension(c0:c1) :: gVt2_ds_1, gVt2_ds_2
	integer, dimension(3) :: p_ext
	real(kind=8),dimension(3) ::  pe_plus, ptwit, phat
	complex(kind=8),dimension(4,4) :: Amat
	integer :: f,x0

	f=f1f2(f2,f1)
	p_ext(:) = 0

	!initialization
	gVt_ds(:) = dcmplx(0.0d0,0.0d0)
	gVt1_ds_1(:) = dcmplx(0.0d0,0.0d0)
	gVt1_ds_2(:) = dcmplx(0.0d0,0.0d0)
	gVt2_ds_1(:) = dcmplx(0.0d0,0.0d0)
	gVt2_ds_2(:) = dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)


	pe_plus(:) = p_ext(:)*(2.0d0*pi)/(l*1.0d0) + theta/l
	ptwit(:) = sin(pe_plus(:))
	phat(:) = 2.0d0*sin(pe_plus(:)/2.0d0)

	Amat = imag*(ptwit(1)*g1+ptwit(2)*g2+ptwit(3)*g3)
	Amat = Amat + 0.5d0*II*( phat(1)**2 + phat(2)**2 + phat(3)**2 )



	call Fermion_Propagator(l,p_ext,S_prop)

	do x0 = c0,c1

	gVt1_ds_1(x0) = gVt1_ds_1(x0)+Trace6(Q5(:,:,f),S_prop(:,:,1,x0+1,f1),Pplus(:,:),S_prop(:,:,x0,0,f2),Amat,S_prop(:,:,0,1,f2))
	gVt1_ds_1(x0) = gVt1_ds_1(x0)+Trace6(Q5(:,:,f),S_prop(:,:,1,x0+1,f1),Pplus(:,:),S_prop(:,:,x0,l,f2),Amat,S_prop(:,:,l,1,f2))

	gVt1_ds_2(x0) = gVt1_ds_2(x0)+Trace6(Q5(:,:,f),S_prop(:,:,1,0,f1),Amat,S_prop(:,:,0,x0+1,f1),Pplus(:,:),S_prop(:,:,x0,1,f2))
	gVt1_ds_2(x0) = gVt1_ds_2(x0)+Trace6(Q5(:,:,f),S_prop(:,:,1,l,f1),Amat,S_prop(:,:,l,x0+1,f1),Pplus(:,:),S_prop(:,:,x0,1,f2))


	gVt2_ds_1(x0) = gVt2_ds_1(x0)+Trace6(Q5(:,:,f),S_prop(:,:,1,x0,f1),Pminus(:,:),S_prop(:,:,x0+1,0,f2),Amat,S_prop(:,:,0,1,f2))
	gVt2_ds_1(x0) = gVt2_ds_1(x0)+Trace6(Q5(:,:,f),S_prop(:,:,1,x0,f1),Pminus(:,:),S_prop(:,:,x0+1,l,f2),Amat,S_prop(:,:,l,1,f2))

	gVt2_ds_2(x0) = gVt2_ds_2(x0)+Trace6(Q5(:,:,f),S_prop(:,:,1,0,f1),Amat,S_prop(:,:,0,x0,f1),Pminus(:,:),S_prop(:,:,x0+1,1,f2))
	gVt2_ds_2(x0) = gVt2_ds_2(x0)+Trace6(Q5(:,:,f),S_prop(:,:,1,l,f1),Amat,S_prop(:,:,l,x0,f1),Pminus(:,:),S_prop(:,:,x0+1,1,f2))


	enddo

	gVt_ds(:) = -dimrep*(gVt1_ds_1(:)+gVt1_ds_2(:)-gVt2_ds_1(:)-gVt2_ds_2(:))/2.0d0

	END SUBROUTINE gVtds

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




	END MODULE OneLoop_gVt


