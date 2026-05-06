




	MODULE OneLoop_gX_p

	USE GammaMatrices
	USE Input
	USE Parameters
	USE Vertices
	USE Propagators
	USE TreeLevel_gX_p

	implicit none

	CONTAINS


	SUBROUTINE gX1a_p(l,g_i,f1,f2,gX_1a_p)

	implicit none

	integer,intent(in) :: l, g_i,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: gX_1a_p
	
	integer :: n1,n2,n3
	integer,dimension(3) :: qloop
	complex(kind=8) :: loop
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop!,K_prop
	complex(kind=8),dimension(c0:c1) :: gX_0
	real(kind=8) :: multi	


	!initialization
	gX_1a_p(:) = dcmplx(0.0d0,0.0d0)
	loop = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	gX_0(:) = dcmplx(0.0d0,0.0d0)

	!Tree level part
	call gX0_p(l,g_i,f1,f2,gX_0)

	!Construct the gluon loop
	do n1=0,l-1			!loop momenta. Optimizable
	do n2=n1,l-1
	do n3=n2,l-1  

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

	gX_1a_p(:) = -0.5d0*C2_R*loop*gX_0(:)/((l*1.0d0)**3)

	END SUBROUTINE gX1a_p

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



	SUBROUTINE gX2_p(l,g_i,f1,f2,gX_2_p)

	implicit none

	integer,intent(in) :: l,g_i,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: gX_2_p
	
	integer :: x0,n1,n2,n3,i
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: loop, gX_0q
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop!,K_prop	
	complex(kind=8),dimension(4,4) :: mat1, mat2, mat3
	real(kind=8) :: multi	
	integer :: f

	!initialization
	gX_2_p(:) = dcmplx(0.0d0,0.0d0)
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
	do n1=0,l-1			!loop momenta. Optimizable
	do n2=n1,l-1
	do n3=n2,l-1  

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

	       mat1 = matmul(Q5p(:,:,f),S_prop(:,:,l-1,x0,f1))
	       mat2 = matmul(mat1,gX(:,:,g_i))
	       mat3 = matmul(mat2,S_prop(:,:,x0,l-1,f2))	

	       gX_0q(x0) = 0.5d0*dimrep*(mat3(1,1)+mat3(2,2)+mat3(3,3)+mat3(4,4))
	       loop(x0) = loop(x0) + D_prop(0,0,l-1,l-1)*gX_0q(x0)*multi

	     enddo


	enddo
	enddo
	enddo

	  gX_2_p(:) = C2_R*loop(:)/((l*1.0d0)**3)


	END SUBROUTINE gX2_p


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




	SUBROUTINE gX3a_p(l,g_i,f1,f2,gX_3a_p)

	implicit none

	integer,intent(in) :: l,g_i,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: gX_3a_p

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

	gX_3a_p(:) = dcmplx(0.0d0,0.0d0)
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
		
	    S_G(:,:,x0,s0) = matmul(S_prop_0(:,:,l-1,x0,f1),matmul(gX(:,:,g_i),S_prop_0(:,:,x0,s0,f2)))

	enddo
	enddo

	!Construct the gluon loop
	do n1=0,l-1			!loop momenta. Optimizable
	do n2=n1,l-1
	do n3=n2,l-1

	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3

	multi = 1.0d0*multiplicity(n1,n2,n3)

  	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

  	call Gluon_propagator(l,qloop,D_prop)	
	call Fermion_Propagator(l,p_ext+qloop,S_prop)


	do mu = 0,3

	do s0 = 0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0 = t0_a,t0_b
	do u0 = u0_a,u0_b

	call V1(l,p_ext,-p_ext-qloop,qloop,s0,t0,u0,mu,V1mu)
	
	do x0 = c0,c1
	  trace = Trace4(S_G(:,:,x0,s0),V1mu(:,:),S_prop(:,:,t0,l-1,f2),Q5p(:,:,f))
	  g_loop(x0) = g_loop(x0) + D_prop(0,mu,l-1,u0)*trace*multi
	enddo

	
	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu


	enddo	!momenta
	enddo
	enddo

	  gX_3a_p(:) = 0.5*C2_R*dimrep*g_loop(:)/(1.0d0*l**3)

	END SUBROUTINE gX3a_p



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE gX3b_p(l,g_i,f1,f2,gX_3b_p)

	implicit none

	integer,intent(in) :: l,g_i,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: gX_3b_p

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

	gX_3b_p(:) = dcmplx(0.0d0,0.0d0)
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
	    S_G(:,:,t0,x0) = matmul(S_prop_0(:,:,t0,x0,f1),matmul(gX(:,:,g_i),S_prop_0(:,:,x0,l-1,f2)))
	enddo
	enddo


	!Construct the gluon loop
	do n1=0,l-1			!loop momenta. Optimizable
	do n2=n1,l-1
	do n3=n2,l-1

	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3

	multi = 1.0d0*multiplicity(n1,n2,n3)

  	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

  	call Gluon_propagator(l,qloop,D_prop)	
	call Fermion_Propagator(l,p_ext+qloop,S_prop)

	do mu=0,3

	do s0=0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0=t0_a,t0_b
	do u0=u0_a,u0_b

	call V1(l,p_ext+qloop,-p_ext,-qloop,s0,t0,u0,mu,V1mu)
	
	do x0=c0,c1
	  trace = Trace4(Q5p(:,:,f),S_prop(:,:,l-1,s0,f1),V1mu(:,:),S_G(:,:,t0,x0))
	  g_loop(x0) = g_loop(x0) + D_prop(mu,0,u0,l-1)*trace*multi
	enddo

	
	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu


	enddo	!momenta
	enddo
	enddo

	  gX_3b_p(:) = -0.5*C2_R*dimrep*g_loop(:)/(1.0d0*l**3)

	END SUBROUTINE gX3b_p



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE gX4a_p(l,g_i,f1,f2,gX_4a_p)

	implicit none

	integer,intent(in) :: l, g_i,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: gX_4a_p

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

	gX_4a_p(:) = dcmplx(0.0d0,0.0d0)
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
	do n1=0,l-1			!loop momenta. Optimizable
	do n2=n1,l-1
	do n3=n2,l-1

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
	      S_G(:,:,t0,x0) = matmul(S_prop(:,:,t0,x0,f1),matmul(gX(:,:,g_i),S_prop(:,:,x0,l-1,f2)))
	   enddo
	   enddo

	do mu=0,3

	do s0=0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0=t0_a,t0_b
	do u0=u0_a,u0_b

	call V1(l,p_ext,-p_ext-qloop,qloop,s0,t0,u0,mu,V1mu)

	do x0=c0,c1
	  trace = Trace4(Q5p(:,:,f),S_prop_0(:,:,l-1,s0,f1),V1mu(:,:),S_G(:,:,t0,x0))
	  g_loop(x0) = g_loop(x0) + D_prop(0,mu,l-1,u0)*trace*multi
	enddo	!x0
	

	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu


	enddo	!momentum
	enddo
	enddo

	gX_4a_p(:) = 0.5*C2_R*dimrep*g_loop(:)/(1.0d0*l**3)

	END SUBROUTINE gX4a_p


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE gX4b_p(l,g_i,f1,f2,gX_4b_p)

	implicit none

	integer,intent(in) :: l, g_i,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: gX_4b_p

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

	gX_4b_p(:) = dcmplx(0.0d0,0.0d0)
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
	do n1=0,l-1			!loop momenta. Optimizable
	do n2=n1,l-1
	do n3=n2,l-1

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
	      S_G(:,:,x0,s0) = matmul(S_prop(:,:,l-1,x0,f1),matmul(gX(:,:,g_i),S_prop(:,:,x0,s0,f2)))
	   enddo
	   enddo


	do mu=0,3

	do s0=0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0=t0_a,t0_b
	do u0=u0_a,u0_b

	call V1(l,p_ext+qloop,-p_ext,-qloop,s0,t0,u0,mu,V1mu)

	do x0 = c0,c1
	  trace = Trace4(Q5p(:,:,f),S_G(:,:,x0,s0),V1mu(:,:),S_prop_0(:,:,t0,l-1,f2))
	  g_loop(x0) = g_loop(x0) + D_prop(mu,0,u0,l-1)*trace*multi
	enddo	!x0
	

	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu



	enddo	!momentum
	enddo
	enddo

	gX_4b_p(:) = -0.5*C2_R*dimrep*g_loop(:)/(1.0d0*l**3)

	END SUBROUTINE gX4b_p
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



	SUBROUTINE gX5a_p(l,g_i,f1,f2,gX_5a_p)

	implicit none

	integer,intent(in) :: l,g_i,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: gX_5a_p
	
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

	!initialization
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_O5_SGS(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	gX_5a_p(:) = dcmplx(0.0d0,0.0d0)
	p(:) = 0
	f=f1f2(f2,f1)

	!build S_O5_SGS

	call Fermion_Propagator(l,p,S_prop_0)

	do x0 = c0,c1
	do s0=0,l
	do t0=0,l
	    S_O5_SGS(:,:,t0,x0,s0) = matmul(S_prop_0(:,:,t0,l-1,f2),matmul(matmul(Q5p(:,:,f),S_prop_0(:,:,l-1,x0,f1)),matmul(gX(:,:,g_i),S_prop_0(:,:,x0,s0,f2))))
	enddo	
	enddo
	enddo


	do x0 = c0,c1


	do n1=0,l-1			!loop momenta. Optimizable
	do n2=n1,l-1
	do n3=n2,l-1

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
	  do nu=0,3

	  call V1(l,p,-p-q,q,s0,t0p,u0,mu,V1mu)
	  call V1(l,p+q,-p,-q,s0p,t0,u0p,nu,V1nu)

	  trace = Trace4(S_O5_SGS(:,:,t0,x0,s0),V1mu(:,:),S_prop(:,:,t0p,s0p,f2),V1nu(:,:))
	  gX_5a_p(x0) = gX_5a_p(x0) + multi*trace*D_prop(nu,mu,u0p,u0)



	  enddo	!nu
	  enddo !mu
	  enddo	!u0p
	  enddo	!u0
	enddo	!t0p
	enddo	!t0
	enddo	!s0p
	enddo	!s0

	enddo	!x0

	enddo	!momenta
	enddo
	enddo

 	gX_5a_p(:) = -0.5d0*dimrep*C2_R*gX_5a_p(:)/((l**3)*1.0d0)

	
	END SUBROUTINE gX5a_p


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



	SUBROUTINE gX5b_p(l,g_i,f1,f2,gX_5b_p)

	implicit none

	integer,intent(in) :: l,g_i,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: gX_5b_p
	

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

	!initialization
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	SGS_O5_S(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	gX_5b_p(:) = dcmplx(0.0d0,0.0d0)
	p(:) = 0
	f=f1f2(f2,f1)

	!build S_O5_SGS

	call Fermion_Propagator(l,p,S_prop_0)

	do x0 = c0,c1
	do s0=0,l
	do t0=0,l
	    SGS_O5_S(:,:,t0,x0,s0) = matmul(S_prop_0(:,:,t0,x0,f1),matmul(matmul(gX(:,:,g_i),S_prop_0(:,:,x0,l-1,f2)),matmul(Q5p(:,:,f),S_prop_0(:,:,l-1,s0,f1))))
	enddo	
	enddo
	enddo



	do x0 = c0,c1


	do n1=0,l-1			!loop momenta. Optimizable
	do n2=n1,l-1
	do n3=n2,l-1

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
	  do nu=0,3

	  call V1(l,p,-p-q,q,s0,t0p,u0,mu,V1mu)
	  call V1(l,p+q,-p,-q,s0p,t0,u0p,nu,V1nu)

	  trace = Trace4(SGS_O5_S(:,:,t0,x0,s0),V1mu(:,:),S_prop(:,:,t0p,s0p,f1),V1nu(:,:))
	  gX_5b_p(x0) = gX_5b_p(x0) + multi*trace*D_prop(nu,mu,u0p,u0)


	  enddo	!nu
	  enddo !mu
	  enddo	!u0p
	  enddo	!u0
	enddo	!t0p
	enddo	!t0
	enddo	!s0p
	enddo	!s0

	enddo	!x0

	enddo	!momenta
	enddo
	enddo

 	gX_5b_p(:) = -0.5d0*dimrep*C2_R*gX_5b_p(:)/((l**3)*1.0d0)

	
	END SUBROUTINE gX5b_p


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



	SUBROUTINE gX6a_p(l,g_i,f1,f2,gX_6a_p)

	implicit none

	integer,intent(in) :: l,g_i,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: gX_6a_p
	

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

	!initialization
	S_prop_0(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_O5_SGS(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	gX_6a_p(:) = dcmplx(0.0d0,0.0d0)
	p(:) = 0
	f=f1f2(f2,f1)

	!build S_O5_SGS

	call Fermion_Propagator(l,p,S_prop_0)

	do x0 = c0,c1
	do s0=0,l
	do t0=0,l
	    S_O5_SGS(:,:,t0,x0,s0) = matmul(S_prop_0(:,:,t0,l-1,f2),matmul(matmul(Q5p(:,:,f),S_prop_0(:,:,l-1,x0,f1)),matmul(gX(:,:,g_i),S_prop_0(:,:,x0,s0,f2))))
	enddo	
	enddo
	enddo


	do x0 = c0,c1

	do n1=0,l-1			!loop momenta. Optimizable
	do n2=n1,l-1
	do n3=n2,l-1

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
	  gX_6a_p(x0) = gX_6a_p(x0) + multi*trace*D_prop(mu,mu,u0,u0)	!There is a delta(mu,nu) in the vertex


!	  enddo	!nu
	  enddo !mu
	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!x0

	enddo	!momenta
	enddo
	enddo

 	gX_6a_p(:) = 0.5d0*dimrep*C2_R*gX_6a_p(:)/((l**3)*1.0d0)

	
	END SUBROUTINE gX6a_p


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE gX6b_p(l,g_i,f1,f2,gX_6b_p)

	implicit none

	integer,intent(in) :: l,g_i,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: gX_6b_p
	

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

	!initialization
	S_prop_0(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	SGS_O5_S(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	gX_6b_p(:) = dcmplx(0.0d0,0.0d0)
	p(:) = 0
	f=f1f2(f2,f1)

	!build S_O5_SGS

	call Fermion_Propagator(l,p,S_prop_0)

	do x0 = c0,c1
	do s0=0,l
	do t0=0,l
	      SGS_O5_S(:,:,t0,x0,s0) = matmul(S_prop_0(:,:,t0,x0,f1),matmul(matmul(gX(:,:,g_i),S_prop_0(:,:,x0,l-1,f2)),matmul(Q5p(:,:,f),S_prop_0(:,:,l-1,s0,f1))))
	enddo	
	enddo
	enddo


	do x0 = c0,c1

	do n1=0,l-1			!loop momenta. Optimizable
	do n2=n1,l-1
	do n3=n2,l-1

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
	  gX_6b_p(x0) = gX_6b_p(x0) + multi*trace*D_prop(mu,mu,u0,u0)	!There is a delta(mu,nu) in the vertex


!	  enddo	!nu
	  enddo !mu
	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!x0

	enddo	!momenta
	enddo
	enddo

 	gX_6b_p(:) = 0.5d0*dimrep*C2_R*gX_6b_p(:)/((l**3)*1.0d0)

	
	END SUBROUTINE gX6b_p


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



	SUBROUTINE gX7_p(l,g_i,f1,f2,gX_7_p)

	implicit none

	integer,intent(in) :: l,g_i,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: gX_7_p
	

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

	!initialization
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	S_O5_S(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	gX_7_p(:) = dcmplx(0.0d0,0.0d0)
	p(:) = 0
	f=f1f2(f2,f1)

	!build S_O5_S

	call Fermion_Propagator(l,p,S_prop_0)

	do s0=0,l
	do t0=0,l
	    S_O5_S(:,:,t0,s0) = matmul(S_prop_0(:,:,t0,l-1,f2),matmul(Q5p(:,:,f),S_prop_0(:,:,l-1,s0,f1)))
	enddo	
	enddo



	do x0 = c0,c1


	do n1=0,l-1			!loop momenta. Optimizable
	do n2=n1,l-1
	do n3=n2,l-1

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
	  do nu=0,3

	  call V1(l,p,-p-q,q,s0,t0p,u0,mu,V1mu)
	  call V1(l,p+q,-p,-q,s0p,t0,u0p,nu,V1nu)

	  trace = Trace6(S_O5_S(:,:,t0,s0),V1mu(:,:),S_prop(:,:,t0p,x0,f1),gX(:,:,g_i),S_prop(:,:,x0,s0p,f2),V1nu(:,:))
	  gX_7_p(x0) = gX_7_p(x0) + multi*trace*D_prop(nu,mu,u0p,u0)


	  enddo	!nu
	  enddo !mu
	  enddo	!u0p
	  enddo	!u0
	enddo	!t0p
	enddo	!t0
	enddo	!s0p
	enddo	!s0

	enddo	!x0

	enddo	!momenta
	enddo
	enddo

 	gX_7_p(:) = -0.5d0*dimrep*C2_R*gX_7_p(:)/((l**3)*1.0d0)

	
	END SUBROUTINE gX7_p


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE gXdm_p(l,g_i,f1,f2,gX_dm_p)

	implicit none

	integer,intent(in) :: l, g_i,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: gX_dm_p
	
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop
	integer, dimension(3) :: p_ext
	integer :: f,s0,x0

	f=f1f2(f2,f1)
	p_ext(:) = 0

	!initialization
	gX_dm_p(:) = dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)


	call Fermion_Propagator(l,p_ext,S_prop)

	do x0 = c0,c1
	do s0 = 0,l

	gX_dm_p(x0) = gX_dm_p(x0)-dimrep*Trace5(Q5p(:,:,f),S_prop(:,:,l-1,x0,f1),gX(:,:,g_i),S_prop(:,:,x0,s0,f2),S_prop(:,:,s0,l-1,f2))

	enddo
	enddo


	END SUBROUTINE gXdm_p


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE gXzf_p(l,g_i,f1,f2,gX_zf_p)

	implicit none

	integer,intent(in) :: l, g_i,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: gX_zf_p
	
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop
	complex(kind=8),dimension(c0:c1) :: gX_zf_1, gX_zf_2
	integer, dimension(3) :: p_ext
	integer :: f,x0

	f=f1f2(f2,f1)
	p_ext(:) = 0

	!initialization
	gX_zf_p(:) = dcmplx(0.0d0,0.0d0)
	gX_zf_1(:) = dcmplx(0.0d0,0.0d0)
	gX_zf_2(:) = dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)


	call Fermion_Propagator(l,p_ext,S_prop)

	do x0 = c0,c1

	gX_zf_1(x0) = gX_zf_1(x0)+Trace5(Q5p(:,:,f),S_prop(:,:,l-1,x0,f1),gX(:,:,g_i),S_prop(:,:,x0,0,f2),S_prop(:,:,0,l-1,f2))
	gX_zf_2(x0) = gX_zf_2(x0)+Trace5(Q5p(:,:,f),S_prop(:,:,l-1,x0,f1),gX(:,:,g_i),S_prop(:,:,x0,l,f2),S_prop(:,:,l,l-1,f2))

	enddo

	gX_zf_p(:) = -dimrep*(gX_zf_1(:)+gX_zf_2(:))

	END SUBROUTINE gXzf_p

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE gXds_p(l,g_i,f1,f2,gX_ds_p)

	implicit none

	integer,intent(in) :: l, g_i,f1,f2
	complex(kind=8),dimension(c0:c1),intent(out) :: gX_ds_p
	
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop
	complex(kind=8),dimension(c0:c1) :: gX_ds_1, gX_ds_2
	integer, dimension(3) :: p_ext
	real(kind=8),dimension(3) ::  pe_plus, ptwit, phat
	complex(kind=8),dimension(4,4) :: Amat
	integer :: f,x0

	f=f1f2(f2,f1)
	p_ext(:) = 0

	!initialization
	gX_ds_p(:) = dcmplx(0.0d0,0.0d0)
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

	gX_ds_1(x0) = gX_ds_1(x0)+Trace6(Q5p(:,:,f),S_prop(:,:,l-1,x0,f1),gX(:,:,g_i),S_prop(:,:,x0,0,f2),Amat,S_prop(:,:,0,l-1,f2))
	gX_ds_2(x0) = gX_ds_2(x0)+Trace6(Q5p(:,:,f),S_prop(:,:,l-1,x0,f1),gX(:,:,g_i),S_prop(:,:,x0,l,f2),Amat,S_prop(:,:,l,l-1,f2))

	enddo

	gX_ds_p(:) = -dimrep*(gX_ds_1(:)+gX_ds_2(:))

	END SUBROUTINE gXds_p

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	END MODULE OneLoop_gX_p


