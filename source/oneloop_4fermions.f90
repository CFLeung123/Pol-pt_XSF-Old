




	MODULE OneLoop_4fermions

	USE GammaMatrices
	USE Input
	USE Parameters
	USE Vertices
	USE Propagators

	implicit none

	CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!	Subroutines for the G1_4f 1loop diagrams
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



	SUBROUTINE G1_4f_1a(l,g_i1,g_i2,flav,G14f_1a)

	USE TreeLevel_4fermions

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G14f_1a
	
	integer :: n1,n2,n3
	integer,dimension(3) :: qloop
	complex(kind=8) :: loop
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop!,K_prop
	complex(kind=8),dimension(c0:c1) :: G14f_0
	real(kind=8) :: multi	


	!initialization
	G14f_1a(:) = dcmplx(0.0d0,0.0d0)
	loop = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	G14f_0(:) = dcmplx(0.0d0,0.0d0)

	!Tree level part
	call G1_4f_0(l,g_i1,g_i2,flav,G14f_0)

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
	   loop = loop + D_prop(0,0,0,0)*multi

	enddo
	enddo
	enddo

	G14f_1a(:) = -0.5d0*C2_R*loop*G14f_0(:)/((l*1.0d0)**3)

	END SUBROUTINE G1_4f_1a

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



	SUBROUTINE G1_4f_1c(l,g_i1,g_i2,flav,G14f_1c)

	USE TreeLevel_4fermions

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G14f_1c
	
	integer :: n1,n2,n3
	integer,dimension(3) :: qloop
	complex(kind=8) :: loop
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop!,K_prop
	complex(kind=8),dimension(c0:c1) :: G14f_0
	real(kind=8) :: multi	


	!initialization
	G14f_1c(:) = dcmplx(0.0d0,0.0d0)
	loop = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	G14f_0(:) = dcmplx(0.0d0,0.0d0)

	!Tree level part
	call G1_4f_0(l,g_i1,g_i2,flav,G14f_0)


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

	G14f_1c(:) = -0.5d0*C2_R*loop*G14f_0(:)/((l*1.0d0)**3)

	END SUBROUTINE G1_4f_1c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE G1_4f_2a(l,g_i1,g_i2,flav,G14f_2a)

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G14f_2a
	
	integer :: n1,n2,n3
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop!,K_prop
	complex(kind=8),dimension(c0:c1) :: loop
	complex(kind=8) :: part1, part2
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp,x0

	nvec_1(:) = 0
	nvec_2(:) = 0

	!initialization
	G14f_2a(:) = dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:) = dcmplx(0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	loop(:) = dcmplx(0.0d0,0.0d0)

	f=f1f2(flav(2),flav(1))
	fp=f1f2(flav(4),flav(3))	

	!Construct the trace in the loop
	do n1=0,l-1			!loop momenta. Optimizable
	do n2=n1,l-1
	do n3=n2,l-1  

	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3

	   multi = 1.0d0*multiplicity(n1,n2,n3)

  	   D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	   S_prop(:,:,:,:,:) = dcmplx(0.0d0)
  	   call Gluon_propagator(l,qloop,D_prop)	
	   call Fermion_Propagator(l,qloop,S_prop)

	   do x0=c0,c1

	     part1 = Trace4(Q5(:,:,f),S_prop(:,:,1,x0,flav(1)),gBil(:,:,g_i1),S_prop(:,:,x0,1,flav(2)))
	     loop(x0) = loop(x0) + D_prop(0,0,0,0)*multi*part1

	  enddo !x0

	enddo
	enddo
	enddo


	!construct the tree-level part outside the loop:
	
	call Fermion_Propagator(l,nvec_2,S_Prop_0)
	do x0=c0,c1

	   part2 = Trace4(Q5p(:,:,fp),S_prop_0(:,:,l-1,x0,flav(3)),gBil(:,:,g_i2),S_prop_0(:,:,x0,l-1,flav(4)))

	   G14f_2a(x0) = loop(x0)*part2

	enddo

	G14f_2a(:) = C2_R*dimrep*dimrep*G14f_2a(:)/((l*1.0d0)**3)

	END SUBROUTINE G1_4f_2a

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G1_4f_2b(l,g_i1,g_i2,flav,G14f_2b)

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G14f_2b
	
	integer :: n1,n2,n3
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop!,K_prop
	complex(kind=8),dimension(c0:c1) :: loop
	complex(kind=8) :: part1, part2
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp,x0

	nvec_1(:) = 0
	nvec_2(:) = 0

	!initialization
	G14f_2b(:) = dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:) = dcmplx(0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	loop(:) = dcmplx(0.0d0,0.0d0)

	f=f1f2(flav(2),flav(1))
	fp=f1f2(flav(4),flav(3))	

	!Construct the trace in the loop
	do n1=0,l-1			!loop momenta. Optimizable
	do n2=n1,l-1
	do n3=n2,l-1  

	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3

	   multi = 1.0d0*multiplicity(n1,n2,n3)

  	   D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	   S_prop(:,:,:,:,:) = dcmplx(0.0d0)
  	   call Gluon_propagator(l,qloop,D_prop)	
	   call Fermion_Propagator(l,qloop,S_prop)

	   do x0=c0,c1

	     part1 = Trace4(Q5p(:,:,fp),S_prop(:,:,l-1,x0,flav(3)),gBil(:,:,g_i2),S_prop(:,:,x0,l-1,flav(4)))
	     loop(x0) = loop(x0) + D_prop(0,0,l-1,l-1)*multi*part1

	  enddo !x0

	enddo
	enddo
	enddo


	!construct the tree-level part outside the loop:
	
	call Fermion_Propagator(l,nvec_1,S_Prop_0)
	do x0=c0,c1

	   part2 = Trace4(Q5(:,:,f),S_prop_0(:,:,1,x0,flav(1)),gBil(:,:,g_i1),S_prop_0(:,:,x0,1,flav(2)))

	   G14f_2b(x0) = part2*loop(x0)

	enddo

	G14f_2b(:) = C2_R*dimrep*dimrep*G14f_2b(:)/((l*1.0d0)**3)

	END SUBROUTINE G1_4f_2b


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G1_4f_3a(l,g_i1,g_i2,flav,G14f_3a)

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G14f_3a

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,c0:c1,0:l) :: S_G
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop, tree
	integer :: n1, n2, n3, x0, t0, u0, s0, mu,i, x0_t
	integer :: t0_a, t0_b, u0_a, u0_b
	complex(kind=8) :: trace
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp

	nvec_1(:) = 0
	nvec_2(:) = 0


	G14f_3a(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_G(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	f=f1f2(flav(2),flav(1))
	fp=f1f2(flav(4),flav(3))

	!part outside the loop, left trace
	call Fermion_Propagator(l,nvec_1,S_prop_0)

	do x0 = c0,c1
	do s0=0,l
		
	    S_G(:,:,x0,s0) = matmul(S_prop_0(:,:,1,x0,flav(1)),matmul(gBil(:,:,g_i1),S_prop_0(:,:,x0,s0,flav(2))))

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
	call Fermion_Propagator(l,nvec_1+qloop,S_prop)


	do mu = 0,3

	do s0 = 0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0 = t0_a,t0_b
	do u0 = u0_a,u0_b

	call V1(l,nvec_1,-nvec_1-qloop,qloop,s0,t0,u0,mu,V1mu)
	
	do x0 = c0,c1
	  trace = Trace4(S_G(:,:,x0,s0),V1mu(:,:),S_prop(:,:,t0,1,flav(2)),Q5(:,:,f))
	  g_loop(x0) = g_loop(x0) + D_prop(0,mu,0,u0)*trace*multi
	enddo

	
	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu


	enddo	!momenta
	enddo
	enddo


	!Prepare the right trace

	call Fermion_Propagator(l,nvec_2,S_Prop_0)
	do x0=c0,c1

	   tree(x0) = Trace4(Q5p(:,:,fp),S_prop_0(:,:,l-1,x0,flav(3)),gBil(:,:,g_i2),S_prop_0(:,:,x0,l-1,flav(4)))

	enddo


	  G14f_3a(:) = -C2_R*dimrep*dimrep*g_loop(:)*tree(:)/(1.0d0*l**3)

	END SUBROUTINE G1_4f_3a



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G1_4f_3b(l,g_i1,g_i2,flav,G14f_3b)

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G14f_3b

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,0:l,c0:c1) :: S_G
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop, tree
	integer :: n1, n2, n3, x0, t0, u0, s0, mu,i, x0_t
	integer :: t0_a, t0_b, u0_a, u0_b
	complex(kind=8) :: trace
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp

	nvec_1(:) = 0
	nvec_2(:) = 0


	G14f_3b(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_G(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	f=f1f2(flav(2),flav(1))
	fp=f1f2(flav(4),flav(3))

	!part outside the loop, left trace
	call Fermion_Propagator(l,nvec_1,S_prop_0)

	do t0 = 0,l
	do x0 = c0,c1		
	    S_G(:,:,t0,x0) = matmul(S_prop_0(:,:,t0,x0,flav(1)),matmul(gBil(:,:,g_i1),S_prop_0(:,:,x0,1,flav(2))))
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
	call Fermion_Propagator(l,nvec_1+qloop,S_prop)

	do mu=0,3

	do s0=0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0=t0_a,t0_b
	do u0=u0_a,u0_b

	call V1(l,nvec_1+qloop,-nvec_1,-qloop,s0,t0,u0,mu,V1mu)
	
	do x0=c0,c1
	  trace = Trace4(Q5(:,:,f),S_prop(:,:,1,s0,flav(1)),V1mu(:,:),S_G(:,:,t0,x0))
	  g_loop(x0) = g_loop(x0) + D_prop(mu,0,u0,0)*trace*multi
	enddo

	
	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu


	enddo	!momenta
	enddo
	enddo


	!Prepare the right trace

	call Fermion_Propagator(l,nvec_2,S_Prop_0)
	do x0=c0,c1

	   tree(x0) = Trace4(Q5p(:,:,fp),S_prop_0(:,:,l-1,x0,flav(3)),gBil(:,:,g_i2),S_prop_0(:,:,x0,l-1,flav(4)))

	enddo


	  G14f_3b(:) = C2_R*dimrep*dimrep*g_loop(:)*tree(:)/(1.0d0*l**3)

	END SUBROUTINE G1_4f_3b



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!








!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G1_4f_3c(l,g_i1,g_i2,flav,G14f_3c)

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G14f_3c

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,c0:c1,0:l) :: S_G
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop, tree
	integer :: n1, n2, n3, x0, t0, u0, s0, mu,i, x0_t
	integer :: t0_a, t0_b, u0_a, u0_b
	complex(kind=8) :: trace
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp

	nvec_1(:) = 0
	nvec_2(:) = 0


	G14f_3c(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_G(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	f=f1f2(flav(2),flav(1))
	fp=f1f2(flav(4),flav(3))

	!part outside the loop, right trace
	call Fermion_Propagator(l,nvec_2,S_prop_0)

	do x0 = c0,c1
	do s0=0,l
		
	    S_G(:,:,x0,s0) = matmul(S_prop_0(:,:,l-1,x0,flav(3)),matmul(gBil(:,:,g_i2),S_prop_0(:,:,x0,s0,flav(4))))

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
	call Fermion_Propagator(l,nvec_2+qloop,S_prop)


	do mu = 0,3

	do s0 = 0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0 = t0_a,t0_b
	do u0 = u0_a,u0_b

	call V1(l,nvec_2,-nvec_2-qloop,qloop,s0,t0,u0,mu,V1mu)
	
	do x0 = c0,c1
	  trace = Trace4(S_G(:,:,x0,s0),V1mu(:,:),S_prop(:,:,t0,l-1,flav(4)),Q5p(:,:,fp))
	  g_loop(x0) = g_loop(x0) + D_prop(0,mu,l-1,u0)*trace*multi
	enddo

	
	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu


	enddo	!momenta
	enddo
	enddo


	!Prepare the left trace

	call Fermion_Propagator(l,nvec_1,S_Prop_0)
	do x0=c0,c1

	   tree(x0) = Trace4(Q5(:,:,f),S_prop_0(:,:,1,x0,flav(1)),gBil(:,:,g_i1),S_prop_0(:,:,x0,1,flav(2)))

	enddo


	  G14f_3c(:) = C2_R*dimrep*dimrep*g_loop(:)*tree(:)/(1.0d0*l**3)

	END SUBROUTINE G1_4f_3c



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G1_4f_3d(l,g_i1,g_i2,flav,G14f_3d)

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G14f_3d

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,0:l,c0:c1) :: S_G
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop, tree
	integer :: n1, n2, n3, x0, t0, u0, s0, mu,i, x0_t
	integer :: t0_a, t0_b, u0_a, u0_b
	complex(kind=8) :: trace
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp

	nvec_1(:) = 0
	nvec_2(:) = 0


	G14f_3d(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_G(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	f=f1f2(flav(2),flav(1))
	fp=f1f2(flav(4),flav(3))

	!part outside the loop, right trace
	call Fermion_Propagator(l,nvec_1,S_prop_0)

	do t0 = 0,l
	do x0 = c0,c1		
	    S_G(:,:,t0,x0) = matmul(S_prop_0(:,:,t0,x0,flav(3)),matmul(gBil(:,:,g_i2),S_prop_0(:,:,x0,l-1,flav(4))))
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
	call Fermion_Propagator(l,nvec_2+qloop,S_prop)

	do mu=0,3

	do s0=0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0=t0_a,t0_b
	do u0=u0_a,u0_b

	call V1(l,nvec_2+qloop,-nvec_2,-qloop,s0,t0,u0,mu,V1mu)
	
	do x0=c0,c1
	  trace = Trace4(Q5p(:,:,fp),S_prop(:,:,l-1,s0,flav(3)),V1mu(:,:),S_G(:,:,t0,x0))
	  g_loop(x0) = g_loop(x0) + D_prop(mu,0,u0,l-1)*trace*multi
	enddo

	
	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu


	enddo	!momenta
	enddo
	enddo


	!Prepare the left trace

	call Fermion_Propagator(l,nvec_1,S_Prop_0)
	do x0=c0,c1

	   tree(x0) = Trace4(Q5(:,:,f),S_prop_0(:,:,1,x0,flav(1)),gBil(:,:,g_i1),S_prop_0(:,:,x0,1,flav(2)))

	enddo


	  G14f_3d(:) = -C2_R*dimrep*dimrep*g_loop(:)*tree(:)/(1.0d0*l**3)

	END SUBROUTINE G1_4f_3d



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G1_4f_4a(l,g_i1,g_i2,flav,G14f_4a)

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G14f_4a

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,0:l,c0:c1) :: S_G
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop, tree
	integer :: n1, n2, n3, x0, t0, u0, s0, mu,i, x0_t
	integer :: t0_a, t0_b, u0_a, u0_b
	complex(kind=8) :: trace
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp

	nvec_1(:) = 0
	nvec_2(:) = 0


	G14f_4a(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_G(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	f=f1f2(flav(2),flav(1))
	fp=f1f2(flav(4),flav(3))

	!part outside the loop, left trace
	call Fermion_Propagator(l,nvec_1,S_prop_0)

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
	   call Fermion_Propagator(l,nvec_1+qloop,S_prop)

	   !Build S_G
	   do t0=0,l
	   do x0=c0,c1
	      S_G(:,:,t0,x0) = matmul(S_prop(:,:,t0,x0,flav(1)),matmul(gBil(:,:,g_i1),S_prop(:,:,x0,1,flav(2))))
	   enddo
	   enddo

	do mu=0,3

	do s0=0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0=t0_a,t0_b
	do u0=u0_a,u0_b

	call V1(l,nvec_1,-nvec_1-qloop,qloop,s0,t0,u0,mu,V1mu)

	do x0=c0,c1
	  trace = Trace4(Q5(:,:,f),S_prop_0(:,:,1,s0,flav(1)),V1mu(:,:),S_G(:,:,t0,x0))
	  g_loop(x0) = g_loop(x0) + D_prop(0,mu,0,u0)*trace*multi
	enddo	!x0
	

	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu


	enddo	!momentum
	enddo
	enddo





	!Prepare the right trace

	call Fermion_Propagator(l,nvec_2,S_Prop_0)
	do x0=c0,c1

	   tree(x0) = Trace4(Q5p(:,:,fp),S_prop_0(:,:,l-1,x0,flav(3)),gBil(:,:,g_i2),S_prop_0(:,:,x0,l-1,flav(4)))

	enddo


	  G14f_4a(:) = -C2_R*dimrep*dimrep*g_loop(:)*tree(:)/(1.0d0*l**3)

	END SUBROUTINE G1_4f_4a



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G1_4f_4b(l,g_i1,g_i2,flav,G14f_4b)

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G14f_4b

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,c0:c1,0:l) :: S_G
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop, tree
	integer :: n1, n2, n3, x0, t0, u0, s0, mu,i, x0_t
	integer :: t0_a, t0_b, u0_a, u0_b
	complex(kind=8) :: trace
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp

	nvec_1(:) = 0
	nvec_2(:) = 0


	G14f_4b(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_G(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	f=f1f2(flav(2),flav(1))
	fp=f1f2(flav(4),flav(3))

	!part outside the loop, left trace
	call Fermion_Propagator(l,nvec_1,S_prop_0)

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
	   call Fermion_Propagator(l,nvec_1+qloop,S_prop)

	   !Build S_G
	   
	   do x0 = c0,c1
	   do s0=0,l
	      S_G(:,:,x0,s0) = matmul(S_prop(:,:,1,x0,flav(1)),matmul(gBil(:,:,g_i1),S_prop(:,:,x0,s0,flav(2))))
	   enddo
	   enddo


	do mu=0,3

	do s0=0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0=t0_a,t0_b
	do u0=u0_a,u0_b

	call V1(l,nvec_1+qloop,-nvec_1,-qloop,s0,t0,u0,mu,V1mu)

	do x0 = c0,c1
	  trace = Trace4(Q5(:,:,f),S_G(:,:,x0,s0),V1mu(:,:),S_prop_0(:,:,t0,1,flav(2)))
	  g_loop(x0) = g_loop(x0) + D_prop(mu,0,u0,0)*trace*multi
	enddo	!x0
	

	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu



	enddo	!momentum
	enddo
	enddo





	!Prepare the right trace

	call Fermion_Propagator(l,nvec_2,S_Prop_0)
	do x0=c0,c1

	   tree(x0) = Trace4(Q5p(:,:,fp),S_prop_0(:,:,l-1,x0,flav(3)),gBil(:,:,g_i2),S_prop_0(:,:,x0,l-1,flav(4)))

	enddo


	  G14f_4b(:) = C2_R*dimrep*dimrep*g_loop(:)*tree(:)/(1.0d0*l**3)

	END SUBROUTINE G1_4f_4b



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G1_4f_4c(l,g_i1,g_i2,flav,G14f_4c)

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G14f_4c

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,0:l,c0:c1) :: S_G
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop, tree
	integer :: n1, n2, n3, x0, t0, u0, s0, mu,i, x0_t
	integer :: t0_a, t0_b, u0_a, u0_b
	complex(kind=8) :: trace
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp

	nvec_1(:) = 0
	nvec_2(:) = 0


	G14f_4c(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_G(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	f=f1f2(flav(2),flav(1))
	fp=f1f2(flav(4),flav(3))

	!part outside the loop, right trace
	call Fermion_Propagator(l,nvec_2,S_prop_0)

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
	   call Fermion_Propagator(l,nvec_2+qloop,S_prop)

	   !Build S_G
	   do t0=0,l
	   do x0=c0,c1
	      S_G(:,:,t0,x0) = matmul(S_prop(:,:,t0,x0,flav(3)),matmul(gBil(:,:,g_i2),S_prop(:,:,x0,l-1,flav(4))))
	   enddo
	   enddo

	do mu=0,3

	do s0=0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0=t0_a,t0_b
	do u0=u0_a,u0_b

	call V1(l,nvec_2,-nvec_2-qloop,qloop,s0,t0,u0,mu,V1mu)

	do x0=c0,c1
	  trace = Trace4(Q5p(:,:,fp),S_prop_0(:,:,l-1,s0,flav(3)),V1mu(:,:),S_G(:,:,t0,x0))
	  g_loop(x0) = g_loop(x0) + D_prop(0,mu,l-1,u0)*trace*multi
	enddo	!x0
	

	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu


	enddo	!momentum
	enddo
	enddo





	!Prepare the right trace

	call Fermion_Propagator(l,nvec_1,S_Prop_0)
	do x0=c0,c1

	   tree(x0) = Trace4(Q5(:,:,f),S_prop_0(:,:,1,x0,flav(1)),gBil(:,:,g_i1),S_prop_0(:,:,x0,1,flav(2)))

	enddo


	  G14f_4c(:) = C2_R*dimrep*dimrep*g_loop(:)*tree(:)/(1.0d0*l**3)

	END SUBROUTINE G1_4f_4c



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G1_4f_4d(l,g_i1,g_i2,flav,G14f_4d)

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G14f_4d

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,c0:c1,0:l) :: S_G
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop, tree
	integer :: n1, n2, n3, x0, t0, u0, s0, mu,i, x0_t
	integer :: t0_a, t0_b, u0_a, u0_b
	complex(kind=8) :: trace
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp

	nvec_1(:) = 0
	nvec_2(:) = 0


	G14f_4d(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_G(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	f=f1f2(flav(2),flav(1))
	fp=f1f2(flav(4),flav(3))

	!part outside the loop, right trace
	call Fermion_Propagator(l,nvec_2,S_prop_0)

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
	   call Fermion_Propagator(l,nvec_2+qloop,S_prop)

	   !Build S_G
	   
	   do x0 = c0,c1
	   do s0=0,l
	      S_G(:,:,x0,s0) = matmul(S_prop(:,:,l-1,x0,flav(3)),matmul(gBil(:,:,g_i2),S_prop(:,:,x0,s0,flav(4))))
	   enddo
	   enddo


	do mu=0,3

	do s0=0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0=t0_a,t0_b
	do u0=u0_a,u0_b

	call V1(l,nvec_2+qloop,-nvec_2,-qloop,s0,t0,u0,mu,V1mu)

	do x0 = c0,c1
	  trace = Trace4(Q5p(:,:,fp),S_G(:,:,x0,s0),V1mu(:,:),S_prop_0(:,:,t0,l-1,flav(4)))
	  g_loop(x0) = g_loop(x0) + D_prop(mu,0,u0,l-1)*trace*multi
	enddo	!x0
	

	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu



	enddo	!momentum
	enddo
	enddo





	!Prepare the right trace

	call Fermion_Propagator(l,nvec_1,S_Prop_0)
	do x0=c0,c1

	   tree(x0) = Trace4(Q5(:,:,f),S_prop_0(:,:,1,x0,flav(1)),gBil(:,:,g_i1),S_prop_0(:,:,x0,1,flav(2)))

	enddo


	  G14f_4d(:) = -C2_R*dimrep*dimrep*g_loop(:)*tree(:)/(1.0d0*l**3)

	END SUBROUTINE G1_4f_4d



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G1_4f_5a(l,g_i1,g_i2,flav,G14f_5a)

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G14f_5a

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,0:l,c0:c1,0:l) :: S_O5_SGS	
	complex(kind=8),dimension(4,4) :: V1mu, V1nu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop, tree
	integer :: x0,i
	integer :: s0,s0p,t0,t0p,u0,u0p
	integer :: t0_a,t0_b,t0p_a,t0p_b,u0_a,u0_b,u0p_a,u0p_b
	integer :: mu, nu
	integer :: n1, n2, n3
	complex(kind=8) :: trace
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp

	nvec_1(:) = 0
	nvec_2(:) = 0


	G14f_5a(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_O5_SGS(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	f=f1f2(flav(2),flav(1))
	fp=f1f2(flav(4),flav(3))

	!part outside the loop, left trace

	!build S_O5_SGS
	call Fermion_Propagator(l,nvec_1,S_prop_0)

	do x0 = c0,c1
	do s0=0,l
	do t0=0,l
	    S_O5_SGS(:,:,t0,x0,s0) = matmul(S_prop_0(:,:,t0,1,flav(2)),matmul(matmul(Q5(:,:,f),S_prop_0(:,:,1,x0,flav(1))),matmul(gBil(:,:,g_i1),S_prop_0(:,:,x0,s0,flav(2)))))
	enddo	
	enddo
	enddo


	!Construct the gluon loop
	do x0 = c0,c1

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
	   call Fermion_Propagator(l,nvec_1+qloop,S_prop)

	   trace = dcmplx(0.0d0,0.0d0)

	   do s0=0,l
	   call loop_limits(l,s0,t0p_a,t0p_b,u0_a,u0_b)
	   do s0p=0,l
	   call loop_limits(l,s0p,t0_a,t0_b,u0p_a,u0p_b)
	   do t0=t0_a,t0_b
	   do t0p=t0p_a,t0p_b		
	   do u0=u0_a,u0_b		
	   do u0p=u0p_a,u0p_b

	  do mu=0,3	
	  do nu=0,3

	  call V1(l,nvec_1,-nvec_1-qloop,qloop,s0,t0p,u0,mu,V1mu)
	  call V1(l,nvec_1+qloop,-nvec_1,-qloop,s0p,t0,u0p,nu,V1nu)

	  trace = Trace4(S_O5_SGS(:,:,t0,x0,s0),V1mu(:,:),S_prop(:,:,t0p,s0p,flav(2)),V1nu(:,:))
	  g_loop(x0) = g_loop(x0) + multi*trace*D_prop(nu,mu,u0p,u0)



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

	enddo	!x0






	!Prepare the right trace

	call Fermion_Propagator(l,nvec_2,S_Prop_0)
	do x0=c0,c1

	   tree(x0) = Trace4(Q5p(:,:,fp),S_prop_0(:,:,l-1,x0,flav(3)),gBil(:,:,g_i2),S_prop_0(:,:,x0,l-1,flav(4)))

	enddo


	  G14f_5a(:) = -C2_R*dimrep*dimrep*g_loop(:)*tree(:)/(1.0d0*l**3)

	END SUBROUTINE G1_4f_5a



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G1_4f_5b(l,g_i1,g_i2,flav,G14f_5b)

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G14f_5b

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,0:l,c0:c1,0:l) :: S_O5_SGS	
	complex(kind=8),dimension(4,4) :: V1mu, V1nu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop, tree
	integer :: x0,i
	integer :: s0,s0p,t0,t0p,u0,u0p
	integer :: t0_a,t0_b,t0p_a,t0p_b,u0_a,u0_b,u0p_a,u0p_b
	integer :: mu, nu
	integer :: n1, n2, n3
	complex(kind=8) :: trace
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp

	nvec_1(:) = 0
	nvec_2(:) = 0


	G14f_5b(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_O5_SGS(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	f=f1f2(flav(2),flav(1))
	fp=f1f2(flav(4),flav(3))

	!part outside the loop, left trace

	!build S_O5_SGS
	call Fermion_Propagator(l,nvec_1,S_prop_0)

	do x0 = c0,c1
	do s0=0,l
	do t0=0,l
	    S_O5_SGS(:,:,t0,x0,s0) = matmul(S_prop_0(:,:,t0,x0,flav(1)),matmul(matmul(gBil(:,:,g_i1),S_prop_0(:,:,x0,1,flav(2))),matmul(Q5(:,:,f),S_prop_0(:,:,1,s0,flav(1)))))
	enddo	
	enddo
	enddo


	!Construct the gluon loop
	do x0 = c0,c1

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
	   call Fermion_Propagator(l,nvec_1+qloop,S_prop)

	   trace = dcmplx(0.0d0,0.0d0)

	   do s0=0,l
	   call loop_limits(l,s0,t0p_a,t0p_b,u0_a,u0_b)
	   do s0p=0,l
	   call loop_limits(l,s0p,t0_a,t0_b,u0p_a,u0p_b)
	   do t0=t0_a,t0_b
	   do t0p=t0p_a,t0p_b		
	   do u0=u0_a,u0_b		
	   do u0p=u0p_a,u0p_b

	  do mu=0,3	
	  do nu=0,3

	  call V1(l,nvec_1,-nvec_1-qloop,qloop,s0,t0p,u0,mu,V1mu)
	  call V1(l,nvec_1+qloop,-nvec_1,-qloop,s0p,t0,u0p,nu,V1nu)

	  trace = Trace4(S_O5_SGS(:,:,t0,x0,s0),V1mu(:,:),S_prop(:,:,t0p,s0p,flav(1)),V1nu(:,:))
	  g_loop(x0) = g_loop(x0) + multi*trace*D_prop(nu,mu,u0p,u0)



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

	enddo	!x0






	!Prepare the right trace

	call Fermion_Propagator(l,nvec_2,S_Prop_0)
	do x0=c0,c1

	   tree(x0) = Trace4(Q5p(:,:,fp),S_prop_0(:,:,l-1,x0,flav(3)),gBil(:,:,g_i2),S_prop_0(:,:,x0,l-1,flav(4)))

	enddo


	  G14f_5b(:) = -C2_R*dimrep*dimrep*g_loop(:)*tree(:)/(1.0d0*l**3)

	END SUBROUTINE G1_4f_5b



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G1_4f_5c(l,g_i1,g_i2,flav,G14f_5c)

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G14f_5c

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,0:l,c0:c1,0:l) :: S_O5_SGS	
	complex(kind=8),dimension(4,4) :: V1mu, V1nu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop, tree
	integer :: x0,i
	integer :: s0,s0p,t0,t0p,u0,u0p
	integer :: t0_a,t0_b,t0p_a,t0p_b,u0_a,u0_b,u0p_a,u0p_b
	integer :: mu, nu
	integer :: n1, n2, n3
	complex(kind=8) :: trace
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp

	nvec_1(:) = 0
	nvec_2(:) = 0


	G14f_5c(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_O5_SGS(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	f=f1f2(flav(2),flav(1))
	fp=f1f2(flav(4),flav(3))

	!part outside the loop, right trace

	!build S_O5_SGS
	call Fermion_Propagator(l,nvec_2,S_prop_0)

	do x0 = c0,c1
	do s0=0,l
	do t0=0,l
	    S_O5_SGS(:,:,t0,x0,s0) = matmul(S_prop_0(:,:,t0,l-1,flav(4)),matmul(matmul(Q5p(:,:,fp),S_prop_0(:,:,l-1,x0,flav(3))),matmul(gBil(:,:,g_i2),S_prop_0(:,:,x0,s0,flav(4)))))
	enddo	
	enddo
	enddo


	!Construct the gluon loop
	do x0 = c0,c1

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
	   call Fermion_Propagator(l,nvec_1+qloop,S_prop)

	   trace = dcmplx(0.0d0,0.0d0)

	   do s0=0,l
	   call loop_limits(l,s0,t0p_a,t0p_b,u0_a,u0_b)
	   do s0p=0,l
	   call loop_limits(l,s0p,t0_a,t0_b,u0p_a,u0p_b)
	   do t0=t0_a,t0_b
	   do t0p=t0p_a,t0p_b		
	   do u0=u0_a,u0_b		
	   do u0p=u0p_a,u0p_b

	  do mu=0,3	
	  do nu=0,3

	  call V1(l,nvec_1,-nvec_1-qloop,qloop,s0,t0p,u0,mu,V1mu)
	  call V1(l,nvec_1+qloop,-nvec_1,-qloop,s0p,t0,u0p,nu,V1nu)

	  trace = Trace4(S_O5_SGS(:,:,t0,x0,s0),V1mu(:,:),S_prop(:,:,t0p,s0p,flav(4)),V1nu(:,:))
	  g_loop(x0) = g_loop(x0) + multi*trace*D_prop(nu,mu,u0p,u0)



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

	enddo	!x0






	!Prepare the left trace

	call Fermion_Propagator(l,nvec_1,S_Prop_0)
	do x0=c0,c1

	   tree(x0) = Trace4(Q5(:,:,f),S_prop_0(:,:,1,x0,flav(1)),gBil(:,:,g_i1),S_prop_0(:,:,x0,1,flav(2)))

	enddo


	  G14f_5c(:) = -C2_R*dimrep*dimrep*g_loop(:)*tree(:)/(1.0d0*l**3)

	END SUBROUTINE G1_4f_5c



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G1_4f_5d(l,g_i1,g_i2,flav,G14f_5d)

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G14f_5d

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,0:l,c0:c1,0:l) :: S_O5_SGS	
	complex(kind=8),dimension(4,4) :: V1mu, V1nu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop, tree
	integer :: x0,i
	integer :: s0,s0p,t0,t0p,u0,u0p
	integer :: t0_a,t0_b,t0p_a,t0p_b,u0_a,u0_b,u0p_a,u0p_b
	integer :: mu, nu
	integer :: n1, n2, n3
	complex(kind=8) :: trace
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp

	nvec_1(:) = 0
	nvec_2(:) = 0


	G14f_5d(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_O5_SGS(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	f=f1f2(flav(2),flav(1))
	fp=f1f2(flav(4),flav(3))

	!part outside the loop, right trace


	call Fermion_Propagator(l,nvec_2,S_prop_0)

	do x0 = c0,c1
	do s0=0,l
	do t0=0,l
	    S_O5_SGS(:,:,t0,x0,s0) = matmul(S_prop_0(:,:,t0,x0,flav(3)),matmul(matmul(gBil(:,:,g_i2),S_prop_0(:,:,x0,l-1,flav(4))),matmul(Q5p(:,:,fp),S_prop_0(:,:,l-1,s0,flav(3)))))
	enddo	
	enddo
	enddo


	do x0 = c0,c1


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
	call Fermion_Propagator(l,nvec_2+qloop,S_prop)

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
	  do nu=0,3

	  call V1(l,nvec_2,-nvec_2-qloop,qloop,s0,t0,u0,mu,V1mu)
	  call V1(l,nvec_2+qloop,-nvec_2,-qloop,s0p,t0p,u0p,nu,V1nu)

	  trace = Trace4(S_O5_SGS(:,:,t0,x0,s0),V1mu(:,:),S_prop(:,:,t0,s0p,flav(3)),V1nu(:,:))
	  g_loop(x0) = g_loop(x0) + multi*trace*D_prop(nu,mu,u0p,u0)



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




	!Prepare the left trace

	call Fermion_Propagator(l,nvec_1,S_Prop_0)
	do x0=c0,c1

	   tree(x0) = Trace4(Q5(:,:,f),S_prop_0(:,:,1,x0,flav(1)),gBil(:,:,g_i1),S_prop(:,:,x0,1,flav(2)))

	enddo


	  G14f_5d(:) = -C2_R*dimrep*dimrep*g_loop(:)*tree(:)/(1.0d0*l**3)

	END SUBROUTINE G1_4f_5d



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G1_4f_6a(l,g_i1,g_i2,flav,G14f_6a)

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G14f_6a

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) ::  S_prop_0
	complex(kind=8),dimension(4,4,0:l,c0:c1,0:l) :: S_O5_SGS	
	complex(kind=8),dimension(4,4) :: V2mu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop, tree
	integer :: x0,i
	integer :: s0,t0,u0
	integer :: t0_a,t0_b,u0_a,u0_b
	integer :: mu, nu
	integer :: n1, n2, n3
	complex(kind=8) :: trace
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp

	nvec_1(:) = 0
	nvec_2(:) = 0


	G14f_6a(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_O5_SGS(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	f=f1f2(flav(2),flav(1))
	fp=f1f2(flav(4),flav(3))

	!part outside the loop, left trace

	!build S_O5_SGS
	call Fermion_Propagator(l,nvec_1,S_prop_0)

	do x0 = c0,c1
	do s0=0,l
	do t0=0,l
	    S_O5_SGS(:,:,t0,x0,s0) = matmul(S_prop_0(:,:,t0,1,flav(2)),matmul(matmul(Q5(:,:,f),S_prop_0(:,:,1,x0,flav(1))),matmul(gBil(:,:,g_i1),S_prop_0(:,:,x0,s0,flav(2)))))
	enddo	
	enddo
	enddo


	!Construct the gluon loop
	do x0 = c0,c1

	do n1=0,l-1			!loop momenta. Optimizable
	do n2=n1,l-1
	do n3=n2,l-1

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

	  call V2(l,nvec_1,-nvec_1,qloop,-qloop,s0,t0,u0,u0,mu,V2mu)


	  trace = Trace2(S_O5_SGS(:,:,t0,x0,s0),V2mu(:,:))
	  g_loop(x0) = g_loop(x0) + multi*trace*D_prop(mu,mu,u0,u0)



	  enddo !mu
	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!momenta
	enddo
	enddo

	enddo	!x0






	!Prepare the right trace

	call Fermion_Propagator(l,nvec_2,S_Prop_0)
	do x0=c0,c1

	   tree(x0) = Trace4(Q5p(:,:,fp),S_prop_0(:,:,l-1,x0,flav(3)),gBil(:,:,g_i2),S_prop_0(:,:,x0,l-1,flav(4)))

	enddo


	  G14f_6a(:) = C2_R*dimrep*dimrep*g_loop(:)*tree(:)/(1.0d0*l**3)

	END SUBROUTINE G1_4f_6a



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G1_4f_6b(l,g_i1,g_i2,flav,G14f_6b)

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G14f_6b

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) ::  S_prop_0
	complex(kind=8),dimension(4,4,0:l,c0:c1,0:l) :: S_O5_SGS	
	complex(kind=8),dimension(4,4) :: V2mu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop, tree
	integer :: x0,i
	integer :: s0,t0,u0
	integer :: t0_a,t0_b,u0_a,u0_b
	integer :: mu, nu
	integer :: n1, n2, n3
	complex(kind=8) :: trace
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp

	nvec_1(:) = 0
	nvec_2(:) = 0


	G14f_6b(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_O5_SGS(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	f=f1f2(flav(2),flav(1))
	fp=f1f2(flav(4),flav(3))

	!part outside the loop, left trace

	!build S_O5_SGS
	call Fermion_Propagator(l,nvec_1,S_prop_0)

	do x0 = c0,c1
	do s0=0,l
	do t0=0,l
	    S_O5_SGS(:,:,t0,x0,s0) = matmul(S_prop_0(:,:,t0,x0,flav(1)),matmul(matmul(gBil(:,:,g_i1),S_prop_0(:,:,x0,1,flav(2))),matmul(Q5(:,:,f),S_prop_0(:,:,1,s0,flav(1)))))
	enddo	
	enddo
	enddo


	!Construct the gluon loop
	do x0 = c0,c1

	do n1=0,l-1			!loop momenta. Optimizable
	do n2=n1,l-1
	do n3=n2,l-1

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

	  call V2(l,nvec_1,-nvec_1,qloop,-qloop,s0,t0,u0,u0,mu,V2mu)


	  trace = Trace2(S_O5_SGS(:,:,t0,x0,s0),V2mu(:,:))
	  g_loop(x0) = g_loop(x0) + multi*trace*D_prop(mu,mu,u0,u0)



	  enddo !mu
	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!momenta
	enddo
	enddo

	enddo	!x0






	!Prepare the right trace

	call Fermion_Propagator(l,nvec_2,S_Prop_0)
	do x0=c0,c1

	   tree(x0) = Trace4(Q5p(:,:,fp),S_prop_0(:,:,l-1,x0,flav(3)),gBil(:,:,g_i2),S_prop_0(:,:,x0,l-1,flav(4)))

	enddo


	  G14f_6b(:) = C2_R*dimrep*dimrep*g_loop(:)*tree(:)/(1.0d0*l**3)

	END SUBROUTINE G1_4f_6b



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G1_4f_6c(l,g_i1,g_i2,flav,G14f_6c)

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G14f_6c

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) ::  S_prop_0
	complex(kind=8),dimension(4,4,0:l,c0:c1,0:l) :: S_O5_SGS	
	complex(kind=8),dimension(4,4) :: V2mu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop, tree
	integer :: x0,i
	integer :: s0,t0,u0
	integer :: t0_a,t0_b,u0_a,u0_b
	integer :: mu, nu
	integer :: n1, n2, n3
	complex(kind=8) :: trace
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp

	nvec_1(:) = 0
	nvec_2(:) = 0


	G14f_6c(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_O5_SGS(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	f=f1f2(flav(2),flav(1))
	fp=f1f2(flav(4),flav(3))

	!part outside the loop, right trace

	!build S_O5_SGS
	call Fermion_Propagator(l,nvec_2,S_prop_0)

	do x0 = c0,c1
	do s0=0,l
	do t0=0,l
	    S_O5_SGS(:,:,t0,x0,s0) = matmul(S_prop_0(:,:,t0,l-1,flav(4)),matmul(matmul(Q5p(:,:,fp),S_prop_0(:,:,l-1,x0,flav(3))),matmul(gBil(:,:,g_i2),S_prop_0(:,:,x0,s0,flav(4)))))
	enddo	
	enddo
	enddo


	!Construct the gluon loop
	do x0 = c0,c1

	do n1=0,l-1			!loop momenta. Optimizable
	do n2=n1,l-1
	do n3=n2,l-1

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

	  call V2(l,nvec_2,-nvec_2,qloop,-qloop,s0,t0,u0,u0,mu,V2mu)


	  trace = Trace2(S_O5_SGS(:,:,t0,x0,s0),V2mu(:,:))
	  g_loop(x0) = g_loop(x0) + multi*trace*D_prop(mu,mu,u0,u0)



	  enddo !mu
	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!momenta
	enddo
	enddo

	enddo	!x0






	!Prepare the left trace

	call Fermion_Propagator(l,nvec_1,S_Prop_0)
	do x0=c0,c1

	   tree(x0) = Trace4(Q5(:,:,f),S_prop_0(:,:,1,x0,flav(1)),gBil(:,:,g_i1),S_prop_0(:,:,x0,1,flav(2)))

	enddo


	  G14f_6c(:) = C2_R*dimrep*dimrep*g_loop(:)*tree(:)/(1.0d0*l**3)

	END SUBROUTINE G1_4f_6c



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G1_4f_6d(l,g_i1,g_i2,flav,G14f_6d)

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G14f_6d

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) ::  S_prop_0
	complex(kind=8),dimension(4,4,0:l,c0:c1,0:l) :: S_O5_SGS	
	complex(kind=8),dimension(4,4) :: V2mu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop, tree
	integer :: x0,i
	integer :: s0,t0,u0
	integer :: t0_a,t0_b,u0_a,u0_b
	integer :: mu, nu
	integer :: n1, n2, n3
	complex(kind=8) :: trace
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp

	nvec_1(:) = 0
	nvec_2(:) = 0


	G14f_6d(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_O5_SGS(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	f=f1f2(flav(2),flav(1))
	fp=f1f2(flav(4),flav(3))

	!part outside the loop, right trace

	!build S_O5_SGS
	call Fermion_Propagator(l,nvec_2,S_prop_0)

	do x0 = c0,c1
	do s0=0,l
	do t0=0,l
	    S_O5_SGS(:,:,t0,x0,s0) = matmul(S_prop_0(:,:,t0,x0,flav(3)),matmul(matmul(gBil(:,:,g_i2),S_prop_0(:,:,x0,l-1,flav(4))),matmul(Q5p(:,:,fp),S_prop_0(:,:,l-1,s0,flav(3)))))
	enddo	
	enddo
	enddo


	!Construct the gluon loop
	do x0 = c0,c1

	do n1=0,l-1			!loop momenta. Optimizable
	do n2=n1,l-1
	do n3=n2,l-1

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

	  call V2(l,nvec_2,-nvec_2,qloop,-qloop,s0,t0,u0,u0,mu,V2mu)


	  trace = Trace2(S_O5_SGS(:,:,t0,x0,s0),V2mu(:,:))
	  g_loop(x0) = g_loop(x0) + multi*trace*D_prop(mu,mu,u0,u0)



	  enddo !mu
	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!momenta
	enddo
	enddo

	enddo	!x0






	!Prepare the right trace

	call Fermion_Propagator(l,nvec_1,S_Prop_0)
	do x0=c0,c1

	   tree(x0) = Trace4(Q5(:,:,f),S_prop_0(:,:,1,x0,flav(1)),gBil(:,:,g_i1),S_prop_0(:,:,x0,1,flav(2)))

	enddo


	  G14f_6d(:) = C2_R*dimrep*dimrep*g_loop(:)*tree(:)/(1.0d0*l**3)

	END SUBROUTINE G1_4f_6d



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G1_4f_7a(l,g_i1,g_i2,flav,G14f_7a)

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G14f_7a

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,0:l,0:l) :: S_O5_S	
	complex(kind=8),dimension(4,4) :: V1mu, V1nu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop, tree
	integer :: x0,i
	integer :: s0,s0p,t0,t0p,u0,u0p
	integer :: t0_a,t0_b,t0p_a,t0p_b,u0_a,u0_b,u0p_a,u0p_b
	integer :: mu, nu
	integer :: n1, n2, n3
	complex(kind=8) :: trace
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp

	nvec_1(:) = 0
	nvec_2(:) = 0


	G14f_7a(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_O5_S(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	f=f1f2(flav(2),flav(1))
	fp=f1f2(flav(4),flav(3))

	!part outside the loop, left trace

	!build S_O5_S
	call Fermion_Propagator(l,nvec_1,S_prop_0)


	do s0=0,l
	do t0=0,l
	    S_O5_S(:,:,t0,s0) = matmul(S_prop_0(:,:,t0,1,flav(2)),matmul(Q5(:,:,f),S_prop_0(:,:,1,s0,flav(1))))
	enddo	
	enddo


	!Construct the gluon loop
	do x0 = c0,c1

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
	   call Fermion_Propagator(l,nvec_1+qloop,S_prop)

	   trace = dcmplx(0.0d0,0.0d0)

	   do s0=0,l
	   call loop_limits(l,s0,t0p_a,t0p_b,u0_a,u0_b)
	   do s0p=0,l
	   call loop_limits(l,s0p,t0_a,t0_b,u0p_a,u0p_b)
	   do t0=t0_a,t0_b
	   do t0p=t0p_a,t0p_b		
	   do u0=u0_a,u0_b		
	   do u0p=u0p_a,u0p_b

	  do mu=0,3	
	  do nu=0,3

	  call V1(l,nvec_1,-nvec_1-qloop,qloop,s0,t0p,u0,mu,V1mu)
	  call V1(l,nvec_1+qloop,-nvec_1,-qloop,s0p,t0,u0p,nu,V1nu)

	  trace = Trace6(S_O5_S(:,:,t0,s0),V1mu(:,:),S_prop(:,:,t0p,x0,flav(1)),gBil(:,:,g_i1),S_prop(:,:,x0,s0p,flav(2)),V1nu(:,:))
	  g_loop(x0) = g_loop(x0) + multi*trace*D_prop(nu,mu,u0p,u0)



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

	enddo	!x0






	!Prepare the right trace

	call Fermion_Propagator(l,nvec_2,S_Prop_0)
	do x0=c0,c1

	   tree(x0) = Trace4(Q5p(:,:,fp),S_prop_0(:,:,l-1,x0,flav(3)),gBil(:,:,g_i2),S_prop_0(:,:,x0,l-1,flav(4)))

	enddo


	  G14f_7a(:) = -C2_R*dimrep*dimrep*g_loop(:)*tree(:)/(1.0d0*l**3)

	END SUBROUTINE G1_4f_7a



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G1_4f_7b(l,g_i1,g_i2,flav,G14f_7b)

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G14f_7b

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,0:l,0:l) :: S_O5_S	
	complex(kind=8),dimension(4,4) :: V1mu, V1nu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop, tree
	integer :: x0,i
	integer :: s0,s0p,t0,t0p,u0,u0p
	integer :: t0_a,t0_b,t0p_a,t0p_b,u0_a,u0_b,u0p_a,u0p_b
	integer :: mu, nu
	integer :: n1, n2, n3
	complex(kind=8) :: trace
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp

	nvec_1(:) = 0
	nvec_2(:) = 0


	G14f_7b(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_O5_S(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	f=f1f2(flav(2),flav(1))
	fp=f1f2(flav(4),flav(3))

	!part outside the loop, right trace

	!build S_O5_S
	call Fermion_Propagator(l,nvec_2,S_prop_0)


	do s0=0,l
	do t0=0,l
	    S_O5_S(:,:,t0,s0) = matmul(S_prop_0(:,:,t0,l-1,flav(4)),matmul(Q5p(:,:,fp),S_prop_0(:,:,l-1,s0,flav(3))))
	enddo	
	enddo


	!Construct the gluon loop
	do x0 = c0,c1

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
	   call Fermion_Propagator(l,nvec_1+qloop,S_prop)

	   trace = dcmplx(0.0d0,0.0d0)

	   do s0=0,l
	   call loop_limits(l,s0,t0p_a,t0p_b,u0_a,u0_b)
	   do s0p=0,l
	   call loop_limits(l,s0p,t0_a,t0_b,u0p_a,u0p_b)
	   do t0=t0_a,t0_b
	   do t0p=t0p_a,t0p_b		
	   do u0=u0_a,u0_b		
	   do u0p=u0p_a,u0p_b

	  do mu=0,3	
	  do nu=0,3

	  call V1(l,nvec_1,-nvec_1-qloop,qloop,s0,t0p,u0,mu,V1mu)
	  call V1(l,nvec_1+qloop,-nvec_1,-qloop,s0p,t0,u0p,nu,V1nu)

	  trace = Trace6(S_O5_S(:,:,t0,s0),V1mu(:,:),S_prop(:,:,t0p,x0,flav(3)),gBil(:,:,g_i2),S_prop(:,:,x0,s0p,flav(4)),V1nu(:,:))
	  g_loop(x0) = g_loop(x0) + multi*trace*D_prop(nu,mu,u0p,u0)



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

	enddo	!x0






	!Prepare the left trace

	call Fermion_Propagator(l,nvec_1,S_Prop_0)
	do x0=c0,c1

	   tree(x0) = Trace4(Q5(:,:,f),S_prop_0(:,:,1,x0,flav(1)),gBil(:,:,g_i1),S_prop_0(:,:,x0,1,flav(2)))

	enddo


	  G14f_7b(:) = -C2_R*dimrep*dimrep*g_loop(:)*tree(:)/(1.0d0*l**3)

	END SUBROUTINE G1_4f_7b



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	Diagrams for G2:
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G2_4f_1a(l,g_i1,g_i2,flav,G24f_1a)

	USE TreeLevel_4fermions

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G24f_1a
	
	integer :: n1,n2,n3
	integer,dimension(3) :: qloop
	complex(kind=8) :: loop
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop!,K_prop
	complex(kind=8),dimension(c0:c1) :: G24f_0
	real(kind=8) :: multi	


	!initialization
	G24f_1a(:) = dcmplx(0.0d0,0.0d0)
	loop = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	G24f_0(:) = dcmplx(0.0d0,0.0d0)

	!Tree level part
	call G2_4f_0(l,g_i1,g_i2,flav,G24f_0)


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
	   loop = loop + D_prop(0,0,0,0)*multi

	enddo
	enddo
	enddo

	G24f_1a(:) = 0.5d0*C2_R*loop*G24f_0(:)/((l*1.0d0)**3)

	END SUBROUTINE G2_4f_1a

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



	SUBROUTINE G2_4f_1c(l,g_i1,g_i2,flav,G24f_1c)

	USE TreeLevel_4fermions

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G24f_1c
	
	integer :: n1,n2,n3
	integer,dimension(3) :: qloop
	complex(kind=8) :: loop
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop!,K_prop
	complex(kind=8),dimension(c0:c1) :: G24f_0
	real(kind=8) :: multi	


	!initialization
	G24f_1c(:) = dcmplx(0.0d0,0.0d0)
	loop = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	G24f_0(:) = dcmplx(0.0d0,0.0d0)

	!Tree level part
	call G2_4f_0(l,g_i1,g_i2,flav,G24f_0)


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

	G24f_1c(:) = 0.5d0*C2_R*loop*G24f_0(:)/((l*1.0d0)**3)

	END SUBROUTINE G2_4f_1c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE G2_4f_2a(l,g_i1,g_i2,flav,G24f_2a)

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G24f_2a
	
	integer :: n1,n2,n3
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop_1, S_prop_2
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop!,K_prop
	complex(kind=8),dimension(c0:c1) :: loop
	complex(kind=8),dimension(4,4) :: mat1, mat2, mat3 ,mat4
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f1,f2,x0

	nvec_1(:) = 0
	nvec_2(:) = 0

	!initialization
	G24f_2a(:) = dcmplx(0.0d0,0.0d0)
	S_prop_1(:,:,:,:,:) = dcmplx(0.0d0)
	S_prop_2(:,:,:,:,:) = dcmplx(0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	loop(:) = dcmplx(0.0d0,0.0d0)

	f1 = f1f2(flav(4),flav(1))
	f2 = f1f2(flav(2),flav(3))
        call Fermion_Propagator(l,nvec_2,S_Prop_2)
	
	!Construct the trace in the loop
	do n1=0,l-1			!loop momenta. Optimizable
	do n2=n1,l-1
	do n3=n2,l-1  

	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3

	   multi = 1.0d0*multiplicity(n1,n2,n3)

  	   D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	   S_prop_1(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

	
	   call Fermion_Propagator(l,nvec_1+qloop,S_Prop_1)
  	   call Gluon_propagator(l,qloop,D_prop)	


	   do x0=c0,c1

	     mat1 = matmul(Q5(:,:,f1),S_prop_1(:,:,1,x0,flav(1)))
	     mat2 = matmul(gBil(:,:,g_i1),S_prop_2(:,:,x0,l-1,flav(2)))
	     mat3 = matmul(Q5p(:,:,f2),S_prop_2(:,:,l-1,x0,flav(3)))
	     mat4 = matmul(gBil(:,:,g_i2),S_prop_1(:,:,x0,1,flav(4)))

	     loop(x0) = loop(x0) + dimrep*D_prop(0,0,0,0)*multi*Trace4(mat1,mat2,mat3,mat4)

	   enddo

	enddo
	enddo
	enddo

	G24f_2a(:) = C2_R*loop(:)/((l*1.0d0)**3)

	END SUBROUTINE G2_4f_2a

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE G2_4f_2b(l,g_i1,g_i2,flav,G24f_2b)

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G24f_2b
	
	integer :: n1,n2,n3
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop_1, S_prop_2
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop!,K_prop
	complex(kind=8),dimension(c0:c1) :: loop
	complex(kind=8),dimension(4,4) :: mat1, mat2, mat3 ,mat4
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f1,f2, x0

	nvec_1(:) = 0
	nvec_2(:) = 0

	!initialization
	G24f_2b(:) = dcmplx(0.0d0,0.0d0)
	S_prop_1(:,:,:,:,:) = dcmplx(0.0d0)
	S_prop_2(:,:,:,:,:) = dcmplx(0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	loop(:) = dcmplx(0.0d0,0.0d0)

	f1 = f1f2(flav(4),flav(1))
	f2 = f1f2(flav(2),flav(3))
	call Fermion_Propagator(l,nvec_1,S_Prop_1)

	!Construct the trace in the loop
	do n1=0,l-1			!loop momenta. Optimizable
	do n2=n1,l-1
	do n3=n2,l-1  

	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3

	   multi = 1.0d0*multiplicity(n1,n2,n3)

  	   D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	   S_prop_1(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	
	   
	   call Fermion_Propagator(l,nvec_2+qloop,S_Prop_2)
  	   call Gluon_propagator(l,qloop,D_prop)	


	   do x0=c0,c1

	     mat1 = matmul(Q5(:,:,f1),S_prop_1(:,:,1,x0,flav(1)))
	     mat2 = matmul(gBil(:,:,g_i1),S_prop_2(:,:,x0,l-1,flav(2)))
	     mat3 = matmul(Q5p(:,:,f2),S_prop_2(:,:,l-1,x0,flav(3)))
	     mat4 = matmul(gBil(:,:,g_i2),S_prop_1(:,:,x0,1,flav(4)))

	     loop(x0) = loop(x0) + dimrep*D_prop(0,0,l-1,l-1)*multi*Trace4(mat1,mat2,mat3,mat4)

	   enddo

	enddo
	enddo
	enddo

	G24f_2b(:) = C2_R*loop(:)/((l*1.0d0)**3)

	END SUBROUTINE G2_4f_2b

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G2_4f_3a(l,g_i1,g_i2,flav,G24f_3a)

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G24f_3a

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,c0:c1,0:l) :: S_G
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop
	integer :: n1, n2, n3, x0, t0, u0, s0, mu,i, x0_t
	integer :: t0_a, t0_b, u0_a, u0_b
	complex(kind=8) :: trace
	complex(kind=8),dimension(4,4) :: mat1, mat2, mat3, mat4
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp

	nvec_1(:) = 0
	nvec_2(:) = 0


	G24f_3a(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_G(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)

	f=f1f2(flav(4),flav(1))
	fp=f1f2(flav(2),flav(3))

	!part outside the loop, left trace
	call Fermion_Propagator(l,nvec_1,S_prop_0)

	do x0 = c0,c1
	do s0=0,l
		
	    mat1 = matmul(Q5(:,:,f),S_prop_0(:,:,1,x0,flav(1)))
	    mat2 = matmul(gBil(:,:,g_i1),S_prop_0(:,:,x0,l-1,flav(2)))
	    mat3 = matmul(Q5p(:,:,fp),S_prop_0(:,:,l-1,x0,flav(3)))
	    mat4 = matmul(gBil(:,:,g_i2),S_prop_0(:,:,x0,s0,flav(4)))

	    S_G(:,:,x0,s0) = matmul(matmul(mat1,mat2),matmul(mat3,mat4))

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
	call Fermion_Propagator(l,nvec_1+qloop,S_prop)


	do mu = 0,3

	do s0 = 0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0 = t0_a,t0_b
	do u0 = u0_a,u0_b

	call V1(l,nvec_1,-nvec_1-qloop,qloop,s0,t0,u0,mu,V1mu)
	
	do x0 = c0,c1
	  trace = Trace3(S_G(:,:,x0,s0),V1mu(:,:),S_prop(:,:,t0,1,flav(4)))
	  g_loop(x0) = g_loop(x0) + D_prop(0,mu,0,u0)*trace*multi
	enddo

	
	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu


	enddo	!momenta
	enddo
	enddo

	G24f_3a(:) = C2_R*dimrep*g_loop(:)/(1.0d0*l**3)

	END SUBROUTINE G2_4f_3a



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G2_4f_3b(l,g_i1,g_i2,flav,G24f_3b)

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G24f_3b

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,0:l,c0:c1) :: S_G
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop
	integer :: n1, n2, n3, x0, t0, u0, s0, mu,i, x0_t
	integer :: t0_a, t0_b, u0_a, u0_b
	complex(kind=8) :: trace
	complex(kind=8),dimension(4,4) :: mat1, mat2, mat3, mat4
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp

	nvec_1(:) = 0
	nvec_2(:) = 0


	G24f_3b(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_G(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)

	f=f1f2(flav(4),flav(1))
	fp=f1f2(flav(2),flav(3))

	!part outside the loop, left trace
	call Fermion_Propagator(l,nvec_1,S_prop_0)

	do t0=0,l
	do x0 = c0,c1
		
	    mat1 = matmul(S_prop_0(:,:,t0,x0,flav(1)),gBil(:,:,g_i1))
	    mat2 = matmul(S_prop_0(:,:,x0,l-1,flav(2)),Q5p(:,:,fp))
	    mat3 = matmul(S_prop_0(:,:,l-1,x0,flav(3)),gBil(:,:,g_i2))
	    mat4 = matmul(S_prop_0(:,:,x0,1,flav(4)),Q5(:,:,f))

	    S_G(:,:,t0,x0) = matmul(matmul(mat1,mat2),matmul(mat3,mat4))

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
	call Fermion_Propagator(l,nvec_1+qloop,S_prop)


	do mu = 0,3

	do s0 = 0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0 = t0_a,t0_b
	do u0 = u0_a,u0_b

	call V1(l,nvec_1+qloop,-nvec_1,-qloop,s0,t0,u0,mu,V1mu)
	
	do x0 = c0,c1
	  trace = Trace3(S_prop(:,:,1,s0,flav(1)),V1mu(:,:),S_G(:,:,t0,x0))
	  g_loop(x0) = g_loop(x0) + D_prop(mu,0,u0,0)*trace*multi
	enddo

	
	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu


	enddo	!momenta
	enddo
	enddo

	G24f_3b(:) = -C2_R*dimrep*g_loop(:)/(1.0d0*l**3)

	END SUBROUTINE G2_4f_3b



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G2_4f_3c(l,g_i1,g_i2,flav,G24f_3c)

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G24f_3c

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,c0:c1,0:l) :: S_G
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop
	integer :: n1, n2, n3, x0, t0, u0, s0, mu,i, x0_t
	integer :: t0_a, t0_b, u0_a, u0_b
	complex(kind=8) :: trace
	complex(kind=8),dimension(4,4) :: mat1, mat2, mat3, mat4
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp

	nvec_1(:) = 0
	nvec_2(:) = 0


	G24f_3c(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_G(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)

	f=f1f2(flav(4),flav(1))
	fp=f1f2(flav(2),flav(3))

	!part outside the loop, left trace
	call Fermion_Propagator(l,nvec_1,S_prop_0)

	do x0 = c0,c1
	do s0=0,l

	    mat1 = matmul(Q5p(:,:,fp),S_prop_0(:,:,l-1,x0,flav(3)))
	    mat2 = matmul(gBil(:,:,g_i2),S_prop_0(:,:,x0,1,flav(4)))		
	    mat3 = matmul(Q5(:,:,f),S_prop_0(:,:,1,x0,flav(1)))
	    mat4 = matmul(gBil(:,:,g_i1),S_prop_0(:,:,x0,s0,flav(2)))

	    S_G(:,:,x0,s0) = matmul(matmul(mat1,mat2),matmul(mat3,mat4))

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
	call Fermion_Propagator(l,nvec_1+qloop,S_prop)


	do mu = 0,3

	do s0 = 0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0 = t0_a,t0_b
	do u0 = u0_a,u0_b

	call V1(l,nvec_1,-nvec_1-qloop,qloop,s0,t0,u0,mu,V1mu)
	
	do x0 = c0,c1
	  trace = Trace3(S_G(:,:,x0,s0),V1mu(:,:),S_prop(:,:,t0,l-1,flav(2)))
	  g_loop(x0) = g_loop(x0) + D_prop(0,mu,l-1,u0)*trace*multi
	enddo

	
	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu


	enddo	!momenta
	enddo
	enddo

	G24f_3c(:) = -C2_R*dimrep*g_loop(:)/(1.0d0*l**3)

	END SUBROUTINE G2_4f_3c



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G2_4f_3d(l,g_i1,g_i2,flav,G24f_3d)

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G24f_3d

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,0:l,c0:c1) :: S_G
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop
	integer :: n1, n2, n3, x0, t0, u0, s0, mu,i, x0_t
	integer :: t0_a, t0_b, u0_a, u0_b
	complex(kind=8) :: trace
	complex(kind=8),dimension(4,4) :: mat1, mat2, mat3, mat4
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp

	nvec_1(:) = 0
	nvec_2(:) = 0


	G24f_3d(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_G(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)

	f=f1f2(flav(4),flav(1))
	fp=f1f2(flav(2),flav(3))

	!part outside the loop, left trace
	call Fermion_Propagator(l,nvec_1,S_prop_0)

	do t0=0,l
	do x0 = c0,c1

	    mat1 = matmul(S_prop_0(:,:,t0,x0,flav(3)),gBil(:,:,g_i2))
	    mat2 = matmul(S_prop_0(:,:,x0,1,flav(4)),Q5(:,:,f))		
	    mat3 = matmul(S_prop_0(:,:,1,x0,flav(1)),gBil(:,:,g_i1))
	    mat4 = matmul(S_prop_0(:,:,x0,l-1,flav(2)),Q5p(:,:,fp))

	    S_G(:,:,t0,x0) = matmul(matmul(mat1,mat2),matmul(mat3,mat4))

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
	call Fermion_Propagator(l,nvec_1+qloop,S_prop)


	do mu = 0,3

	do s0 = 0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0 = t0_a,t0_b
	do u0 = u0_a,u0_b

	call V1(l,nvec_1+qloop,-nvec_1,-qloop,s0,t0,u0,mu,V1mu)
	
	do x0 = c0,c1
	  trace = Trace3(S_prop(:,:,l-1,s0,flav(3)),V1mu(:,:),S_G(:,:,t0,x0))
	  g_loop(x0) = g_loop(x0) + D_prop(mu,0,u0,l-1)*trace*multi
	enddo

	
	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu


	enddo	!momenta
	enddo
	enddo

	G24f_3d(:) = C2_R*dimrep*g_loop(:)/(1.0d0*l**3)

	END SUBROUTINE G2_4f_3d



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G2_4f_4a(l,g_i1,g_i2,flav,G24f_4a)

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G24f_4a

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,c0:c1) :: S_G
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop
	integer :: n1, n2, n3, x0, t0, u0, s0, mu,i, x0_t
	integer :: t0_a, t0_b, u0_a, u0_b
	complex(kind=8) :: trace
	complex(kind=8),dimension(4,4) :: mat1, mat2
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp

	nvec_1(:) = 0
	nvec_2(:) = 0

	G24f_4a(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_G(:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)

	f=f1f2(flav(4),flav(1))
	fp=f1f2(flav(2),flav(3))

	!part outside the loop, left trace
	call Fermion_Propagator(l,nvec_1,S_prop_0)

	do x0 = c0,c1
		
	    mat1 = matmul(gBil(:,:,g_i1),S_prop_0(:,:,x0,l-1,flav(2)))
	    mat2 = matmul(Q5p(:,:,fp),S_prop_0(:,:,l-1,x0,flav(3)))
	    S_G(:,:,x0) = matmul(matmul(mat1,mat2),gBil(:,:,g_i2))

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
	call Fermion_Propagator(l,nvec_1+qloop,S_prop)

	do mu = 0,3

	do s0 = 0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0 = t0_a,t0_b
	do u0 = u0_a,u0_b

	call V1(l,nvec_1,-nvec_1-qloop,qloop,s0,t0,u0,mu,V1mu)
	
	do x0 = c0,c1

	  trace = Trace6(S_G(:,:,x0),S_prop(:,:,x0,1,flav(4)),Q5(:,:,f),S_prop_0(:,:,1,s0,flav(1)),V1mu(:,:),S_prop(:,:,t0,x0,flav(1)))
	  g_loop(x0) = g_loop(x0) + D_prop(0,mu,0,u0)*trace*multi

	enddo

	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu

	enddo	!momenta
	enddo
	enddo

	G24f_4a(:) = C2_R*dimrep*g_loop(:)/(1.0d0*l**3)

	END SUBROUTINE G2_4f_4a



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G2_4f_4b(l,g_i1,g_i2,flav,G24f_4b)

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G24f_4b

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,c0:c1) :: S_G
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop
	integer :: n1, n2, n3, x0, t0, u0, s0, mu,i, x0_t
	integer :: t0_a, t0_b, u0_a, u0_b
	complex(kind=8) :: trace
	complex(kind=8),dimension(4,4) :: mat1, mat2
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp

	nvec_1(:) = 0
	nvec_2(:) = 0

	G24f_4b(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_G(:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)

	f=f1f2(flav(4),flav(1))
	fp=f1f2(flav(2),flav(3))

	!part outside the loop, left trace
	call Fermion_Propagator(l,nvec_1,S_prop_0)

	do x0 = c0,c1

	    mat1 = matmul(gBil(:,:,g_i1),S_prop_0(:,:,x0,l-1,flav(2)))
	    mat1 = matmul(Q5p(:,:,fp),S_prop_0(:,:,l-1,x0,flav(3)))		
	    S_G(:,:,x0) = matmul(matmul(mat1,mat2),gBil(:,:,g_i2))

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
	call Fermion_Propagator(l,nvec_1+qloop,S_prop)

	do mu = 0,3

	do s0 = 0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0 = t0_a,t0_b
	do u0 = u0_a,u0_b

	call V1(l,nvec_1+qloop,-nvec_1,-qloop,s0,t0,u0,mu,V1mu)
	
	do x0 = c0,c1

	  trace = Trace6(S_G(:,:,x0),S_prop(:,:,x0,s0,flav(4)),V1mu(:,:),S_prop_0(:,:,t0,1,flav(4)),Q5(:,:,f),S_prop(:,:,1,x0,flav(1)))
	  g_loop(x0) = g_loop(x0) + D_prop(mu,0,u0,0)*trace*multi

	enddo

	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu

	enddo	!momenta
	enddo
	enddo

	G24f_4b(:) = -C2_R*dimrep*g_loop(:)/(1.0d0*l**3)

	END SUBROUTINE G2_4f_4b



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G2_4f_4c(l,g_i1,g_i2,flav,G24f_4c)

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G24f_4c

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,c0:c1) :: S_G
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop
	integer :: n1, n2, n3, x0, t0, u0, s0, mu,i, x0_t
	integer :: t0_a, t0_b, u0_a, u0_b
	complex(kind=8) :: trace
	complex(kind=8),dimension(4,4) :: mat1, mat2
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp

	nvec_1(:) = 0
	nvec_2(:) = 0

	G24f_4c(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_G(:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)

	f=f1f2(flav(4),flav(1))
	fp=f1f2(flav(2),flav(3))

	!part outside the loop, left trace
	call Fermion_Propagator(l,nvec_1,S_prop_0)

	do x0 = c0,c1
		
	    mat1 = matmul(gBil(:,:,g_i2),S_prop_0(:,:,x0,1,flav(4)))
	    mat2 = matmul(Q5(:,:,f),S_prop_0(:,:,1,x0,flav(1)))
	    S_G(:,:,x0) = matmul(matmul(mat1,mat2),gBil(:,:,g_i1))

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
	call Fermion_Propagator(l,nvec_1+qloop,S_prop)

	do mu = 0,3

	do s0 = 0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0 = t0_a,t0_b
	do u0 = u0_a,u0_b

	call V1(l,nvec_1,-nvec_1-qloop,qloop,s0,t0,u0,mu,V1mu)
	
	do x0 = c0,c1

	  trace = Trace6(S_G(:,:,x0),S_prop(:,:,x0,l-1,flav(2)),Q5p(:,:,fp),S_prop_0(:,:,l-1,s0,flav(3)),V1mu(:,:),S_prop(:,:,t0,x0,flav(3)))
	  g_loop(x0) = g_loop(x0) + D_prop(0,mu,l-1,u0)*trace*multi

	enddo

	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu

	enddo	!momenta
	enddo
	enddo

	G24f_4c(:) = -C2_R*dimrep*g_loop(:)/(1.0d0*l**3)

	END SUBROUTINE G2_4f_4c



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G2_4f_4d(l,g_i1,g_i2,flav,G24f_4d)

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G24f_4d

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,c0:c1) :: S_G
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop
	integer :: n1, n2, n3, x0, t0, u0, s0, mu,i, x0_t
	integer :: t0_a, t0_b, u0_a, u0_b
	complex(kind=8) :: trace
	complex(kind=8),dimension(4,4) :: mat1, mat2
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp

	nvec_1(:) = 0
	nvec_2(:) = 0

	G24f_4d(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_G(:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)

	f=f1f2(flav(4),flav(1))
	fp=f1f2(flav(2),flav(3))

	!part outside the loop, left trace
	call Fermion_Propagator(l,nvec_1,S_prop_0)

	do x0 = c0,c1

	    mat1 = matmul(gBil(:,:,g_i2),S_prop_0(:,:,x0,1,flav(4)))
	    mat1 = matmul(Q5(:,:,f),S_prop_0(:,:,1,x0,flav(1)))		
	    S_G(:,:,x0) = matmul(matmul(mat1,mat2),gBil(:,:,g_i1))

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
	call Fermion_Propagator(l,nvec_1+qloop,S_prop)

	do mu = 0,3

	do s0 = 0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0 = t0_a,t0_b
	do u0 = u0_a,u0_b

	call V1(l,nvec_1+qloop,-nvec_1,-qloop,s0,t0,u0,mu,V1mu)
	
	do x0 = c0,c1

	  trace = Trace6(S_G(:,:,x0),S_prop(:,:,x0,s0,flav(2)),V1mu(:,:),S_prop_0(:,:,t0,l-1,flav(2)),Q5p(:,:,fp),S_prop(:,:,l-1,x0,flav(3)))
	  g_loop(x0) = g_loop(x0) + D_prop(mu,0,u0,l-1)*trace*multi

	enddo

	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu

	enddo	!momenta
	enddo
	enddo

	G24f_4d(:) = C2_R*dimrep*g_loop(:)/(1.0d0*l**3)

	END SUBROUTINE G2_4f_4d



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G2_4f_5a(l,g_i1,g_i2,flav,G24f_5a)

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G24f_5a

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,0:l,c0:c1,0:l) :: S_O5_SGS	
	complex(kind=8),dimension(4,4) :: V1mu, V1nu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop
	complex(kind=8),dimension(4,4) :: mat1, mat2, mat3, mat4
	integer :: x0,i
	integer :: s0,s0p,t0,t0p,u0,u0p
	integer :: t0_a,t0_b,t0p_a,t0p_b,u0_a,u0_b,u0p_a,u0p_b
	integer :: mu, nu
	integer :: n1, n2, n3
	complex(kind=8) :: trace
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp

	nvec_1(:) = 0
	nvec_2(:) = 0


	G24f_5a(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_O5_SGS(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)

	f=f1f2(flav(4),flav(1))
	fp=f1f2(flav(2),flav(3))

	!part outside the loop, left trace

	!build S_O5_SGS
	call Fermion_Propagator(l,nvec_1,S_prop_0)

	do t0=0,l
	do x0 = c0,c1
	do s0=0,l

	    mat1 = matmul(S_prop_0(:,:,t0,1,flav(4)),Q5(:,:,f))
	    mat2 = matmul(S_prop_0(:,:,1,x0,flav(1)),gBil(:,:,g_i1))		
	    mat3 = matmul(S_prop_0(:,:,x0,l-1,flav(2)),Q5p(:,:,fp))
	    mat4 = matmul(matmul(S_prop_0(:,:,l-1,x0,flav(3)),gBil(:,:,g_i2)),S_prop_0(:,:,x0,s0,flav(4)))

	    S_O5_SGS(:,:,t0,x0,s0) = matmul(matmul(mat1,mat2),matmul(mat3,mat4))

	enddo
	enddo
	enddo



	!Construct the gluon loop
	do x0 = c0,c1

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
	   call Fermion_Propagator(l,nvec_1+qloop,S_prop)

	   trace = dcmplx(0.0d0,0.0d0)

	   do s0=0,l
	   call loop_limits(l,s0,t0p_a,t0p_b,u0_a,u0_b)
	   do s0p=0,l
	   call loop_limits(l,s0p,t0_a,t0_b,u0p_a,u0p_b)
	   do t0=t0_a,t0_b
	   do t0p=t0p_a,t0p_b		
	   do u0=u0_a,u0_b		
	   do u0p=u0p_a,u0p_b

	  do mu=0,3	
	  do nu=0,3

	  call V1(l,nvec_1,-nvec_1-qloop,qloop,s0,t0p,u0,mu,V1mu)
	  call V1(l,nvec_1+qloop,-nvec_1,-qloop,s0p,t0,u0p,nu,V1nu)

	  trace = Trace4(S_O5_SGS(:,:,t0,x0,s0),V1mu(:,:),S_prop(:,:,t0p,s0p,flav(4)),V1nu(:,:))
	  g_loop(x0) = g_loop(x0) + multi*trace*D_prop(nu,mu,u0p,u0)



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

	enddo	!x0

	G24f_5a(:) = C2_R*dimrep*g_loop(:)/(1.0d0*l**3)

	END SUBROUTINE G2_4f_5a



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G2_4f_5b(l,g_i1,g_i2,flav,G24f_5b)

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G24f_5b

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,0:l,c0:c1,0:l) :: S_O5_SGS	
	complex(kind=8),dimension(4,4) :: V1mu, V1nu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop
	complex(kind=8),dimension(4,4) :: mat1, mat2, mat3, mat4
	integer :: x0,i
	integer :: s0,s0p,t0,t0p,u0,u0p
	integer :: t0_a,t0_b,t0p_a,t0p_b,u0_a,u0_b,u0p_a,u0p_b
	integer :: mu, nu
	integer :: n1, n2, n3
	complex(kind=8) :: trace
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp

	nvec_1(:) = 0
	nvec_2(:) = 0


	G24f_5b(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_O5_SGS(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)

	f=f1f2(flav(4),flav(1))
	fp=f1f2(flav(2),flav(3))

	!part outside the loop, left trace

	!build S_O5_SGS
	call Fermion_Propagator(l,nvec_1,S_prop_0)

	do t0=0,l
	do x0 = c0,c1
	do s0=0,l

	    mat1 = matmul(S_prop_0(:,:,t0,x0,flav(1)),gBil(:,:,g_i1))
	    mat2 = matmul(S_prop_0(:,:,x0,l-1,flav(2)),Q5p(:,:,fp))		
	    mat3 = matmul(S_prop_0(:,:,l-1,x0,flav(3)),gBil(:,:,g_i2))
	    mat4 = matmul(matmul(S_prop_0(:,:,x0,1,flav(4)),Q5(:,:,f)),S_prop_0(:,:,1,s0,flav(1)))

	    S_O5_SGS(:,:,t0,x0,s0) = matmul(matmul(mat1,mat2),matmul(mat3,mat4))

	enddo
	enddo
	enddo



	!Construct the gluon loop
	do x0 = c0,c1

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
	   call Fermion_Propagator(l,nvec_1+qloop,S_prop)

	   trace = dcmplx(0.0d0,0.0d0)

	   do s0=0,l
	   call loop_limits(l,s0,t0p_a,t0p_b,u0_a,u0_b)
	   do s0p=0,l
	   call loop_limits(l,s0p,t0_a,t0_b,u0p_a,u0p_b)
	   do t0=t0_a,t0_b
	   do t0p=t0p_a,t0p_b		
	   do u0=u0_a,u0_b		
	   do u0p=u0p_a,u0p_b

	  do mu=0,3	
	  do nu=0,3

	  call V1(l,nvec_1,-nvec_1-qloop,qloop,s0,t0p,u0,mu,V1mu)
	  call V1(l,nvec_1+qloop,-nvec_1,-qloop,s0p,t0,u0p,nu,V1nu)

	  trace = Trace4(S_O5_SGS(:,:,t0,x0,s0),V1mu(:,:),S_prop(:,:,t0p,s0p,flav(1)),V1nu(:,:))
	  g_loop(x0) = g_loop(x0) + multi*trace*D_prop(nu,mu,u0p,u0)



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

	enddo	!x0

	G24f_5b(:) = C2_R*dimrep*g_loop(:)/(1.0d0*l**3)

	END SUBROUTINE G2_4f_5b



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G2_4f_5c(l,g_i1,g_i2,flav,G24f_5c)

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G24f_5c

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,0:l,c0:c1,0:l) :: S_O5_SGS	
	complex(kind=8),dimension(4,4) :: V1mu, V1nu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop
	complex(kind=8),dimension(4,4) :: mat1, mat2, mat3, mat4
	integer :: x0,i
	integer :: s0,s0p,t0,t0p,u0,u0p
	integer :: t0_a,t0_b,t0p_a,t0p_b,u0_a,u0_b,u0p_a,u0p_b
	integer :: mu, nu
	integer :: n1, n2, n3
	complex(kind=8) :: trace
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp

	nvec_1(:) = 0
	nvec_2(:) = 0


	G24f_5c(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_O5_SGS(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)

	f=f1f2(flav(4),flav(1))
	fp=f1f2(flav(1),flav(3))

	!part outside the loop, left trace

	!build S_O5_SGS
	call Fermion_Propagator(l,nvec_1,S_prop_0)

	do t0=0,l
	do x0 = c0,c1
	do s0=0,l

	    mat1 = matmul(S_prop_0(:,:,t0,l-1,flav(2)),Q5p(:,:,fp))
	    mat2 = matmul(S_prop_0(:,:,l-1,x0,flav(3)),gBil(:,:,g_i2))		
	    mat3 = matmul(S_prop_0(:,:,x0,1,flav(4)),Q5(:,:,f))
	    mat4 = matmul(matmul(S_prop_0(:,:,1,x0,flav(1)),gBil(:,:,g_i1)),S_prop_0(:,:,x0,s0,flav(2)))

	    S_O5_SGS(:,:,t0,x0,s0) = matmul(matmul(mat1,mat2),matmul(mat3,mat4))

	enddo
	enddo
	enddo



	!Construct the gluon loop
	do x0 = c0,c1

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
	   call Fermion_Propagator(l,nvec_1+qloop,S_prop)

	   trace = dcmplx(0.0d0,0.0d0)

	   do s0=0,l
	   call loop_limits(l,s0,t0p_a,t0p_b,u0_a,u0_b)
	   do s0p=0,l
	   call loop_limits(l,s0p,t0_a,t0_b,u0p_a,u0p_b)
	   do t0=t0_a,t0_b
	   do t0p=t0p_a,t0p_b		
	   do u0=u0_a,u0_b		
	   do u0p=u0p_a,u0p_b

	  do mu=0,3	
	  do nu=0,3

	  call V1(l,nvec_1,-nvec_1-qloop,qloop,s0,t0p,u0,mu,V1mu)
	  call V1(l,nvec_1+qloop,-nvec_1,-qloop,s0p,t0,u0p,nu,V1nu)

	  trace = Trace4(S_O5_SGS(:,:,t0,x0,s0),V1mu(:,:),S_prop(:,:,t0p,s0p,flav(2)),V1nu(:,:))
	  g_loop(x0) = g_loop(x0) + multi*trace*D_prop(nu,mu,u0p,u0)



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

	enddo	!x0

	G24f_5c(:) = C2_R*dimrep*g_loop(:)/(1.0d0*l**3)

	END SUBROUTINE G2_4f_5c



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G2_4f_5d(l,g_i1,g_i2,flav,G24f_5d)

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G24f_5d

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,0:l,c0:c1,0:l) :: S_O5_SGS	
	complex(kind=8),dimension(4,4) :: V1mu, V1nu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop
	complex(kind=8),dimension(4,4) :: mat1, mat2, mat3, mat4
	integer :: x0,i
	integer :: s0,s0p,t0,t0p,u0,u0p
	integer :: t0_a,t0_b,t0p_a,t0p_b,u0_a,u0_b,u0p_a,u0p_b
	integer :: mu, nu
	integer :: n1, n2, n3
	complex(kind=8) :: trace
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp

	nvec_1(:) = 0
	nvec_2(:) = 0


	G24f_5d(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_O5_SGS(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)

	f=f1f2(flav(4),flav(1))
	fp=f1f2(flav(2),flav(3))

	!part outside the loop, left trace

	!build S_O5_SGS
	call Fermion_Propagator(l,nvec_1,S_prop_0)

	do t0=0,l
	do x0 = c0,c1
	do s0=0,l

	    mat1 = matmul(S_prop_0(:,:,t0,x0,flav(3)),gBil(:,:,g_i2))
	    mat2 = matmul(S_prop_0(:,:,x0,1,flav(4)),Q5(:,:,f))		
	    mat3 = matmul(S_prop_0(:,:,1,x0,flav(1)),gBil(:,:,g_i1))
	    mat4 = matmul(matmul(S_prop_0(:,:,x0,1,flav(2)),Q5p(:,:,fp)),S_prop_0(:,:,l-1,s0,flav(3)))

	    S_O5_SGS(:,:,t0,x0,s0) = matmul(matmul(mat1,mat2),matmul(mat3,mat4))

	enddo
	enddo
	enddo



	!Construct the gluon loop
	do x0 = c0,c1

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
	   call Fermion_Propagator(l,nvec_1+qloop,S_prop)

	   trace = dcmplx(0.0d0,0.0d0)

	   do s0=0,l
	   call loop_limits(l,s0,t0p_a,t0p_b,u0_a,u0_b)
	   do s0p=0,l
	   call loop_limits(l,s0p,t0_a,t0_b,u0p_a,u0p_b)
	   do t0=t0_a,t0_b
	   do t0p=t0p_a,t0p_b		
	   do u0=u0_a,u0_b		
	   do u0p=u0p_a,u0p_b

	  do mu=0,3	
	  do nu=0,3

	  call V1(l,nvec_1,-nvec_1-qloop,qloop,s0,t0p,u0,mu,V1mu)
	  call V1(l,nvec_1+qloop,-nvec_1,-qloop,s0p,t0,u0p,nu,V1nu)

	  trace = Trace4(S_O5_SGS(:,:,t0,x0,s0),V1mu(:,:),S_prop(:,:,t0p,s0p,flav(3)),V1nu(:,:))
	  g_loop(x0) = g_loop(x0) + multi*trace*D_prop(nu,mu,u0p,u0)



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

	enddo	!x0

	G24f_5d(:) = C2_R*dimrep*g_loop(:)/(1.0d0*l**3)

	END SUBROUTINE G2_4f_5d



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G2_4f_6a(l,g_i1,g_i2,flav,G24f_6a)

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G24f_6a

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) ::  S_prop_0
	complex(kind=8),dimension(4,4,0:l,c0:c1,0:l) :: S_O5_SGS	
	complex(kind=8),dimension(4,4) :: V2mu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop, tree
	integer :: x0,i
	integer :: s0,t0,u0
	integer :: t0_a,t0_b,u0_a,u0_b
	integer :: mu, nu
	integer :: n1, n2, n3
	complex(kind=8) :: trace
	complex(kind=8),dimension(4,4) :: mat1, mat2, mat3, mat4
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp

	nvec_1(:) = 0
	nvec_2(:) = 0


	G24f_6a(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_O5_SGS(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	f=f1f2(flav(4),flav(1))
	fp=f1f2(flav(2),flav(3))

	!part outside the loop

	!build S_O5_SGS
	call Fermion_Propagator(l,nvec_1,S_prop_0)

	do t0=0,l
	do x0 = c0,c1
	do s0=0,l

	    mat1 = matmul(S_prop_0(:,:,t0,1,flav(4)),Q5(:,:,f))
	    mat2 = matmul(S_prop_0(:,:,1,x0,flav(1)),gBil(:,:,g_i1))		
	    mat3 = matmul(S_prop_0(:,:,x0,l-1,flav(2)),Q5p(:,:,fp))
	    mat4 = matmul(matmul(S_prop_0(:,:,l-1,x0,flav(3)),gBil(:,:,g_i2)),S_prop_0(:,:,x0,s0,flav(4)))

	    S_O5_SGS(:,:,t0,x0,s0) = matmul(matmul(mat1,mat2),matmul(mat3,mat4))

	enddo
	enddo
	enddo



	!Construct the gluon loop
	do x0 = c0,c1

	do n1=0,l-1			!loop momenta. Optimizable
	do n2=n1,l-1
	do n3=n2,l-1

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

	  call V2(l,nvec_1,-nvec_1,qloop,-qloop,s0,t0,u0,u0,mu,V2mu)


	  trace = Trace2(S_O5_SGS(:,:,t0,x0,s0),V2mu(:,:))
	  g_loop(x0) = g_loop(x0) + multi*trace*D_prop(mu,mu,u0,u0)



	  enddo !mu
	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!momenta
	enddo
	enddo

	enddo	!x0


	  G24f_6a(:) = -C2_R*dimrep*g_loop(:)/(1.0d0*l**3)

	END SUBROUTINE G2_4f_6a



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G2_4f_6b(l,g_i1,g_i2,flav,G24f_6b)

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G24f_6b

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) ::  S_prop_0
	complex(kind=8),dimension(4,4,0:l,c0:c1,0:l) :: S_O5_SGS	
	complex(kind=8),dimension(4,4) :: V2mu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop, tree
	integer :: x0,i
	integer :: s0,t0,u0
	integer :: t0_a,t0_b,u0_a,u0_b
	integer :: mu, nu
	integer :: n1, n2, n3
	complex(kind=8) :: trace
	complex(kind=8),dimension(4,4) :: mat1, mat2, mat3, mat4
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp

	nvec_1(:) = 0
	nvec_2(:) = 0


	G24f_6b(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_O5_SGS(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	f=f1f2(flav(4),flav(1))
	fp=f1f2(flav(2),flav(3))

	!part outside the loop, left trace

	!build S_O5_SGS
	call Fermion_Propagator(l,nvec_1,S_prop_0)

	do t0=0,l
	do x0 = c0,c1
	do s0=0,l

	    mat1 = matmul(S_prop_0(:,:,t0,x0,flav(1)),gBil(:,:,g_i1))
	    mat2 = matmul(S_prop_0(:,:,x0,l-1,flav(2)),Q5p(:,:,fp))		
	    mat3 = matmul(S_prop_0(:,:,l-1,x0,flav(3)),gBil(:,:,g_i2))
	    mat4 = matmul(matmul(S_prop_0(:,:,x0,1,flav(4)),Q5(:,:,f)),S_prop_0(:,:,1,s0,flav(1)))

	    S_O5_SGS(:,:,t0,x0,s0) = matmul(matmul(mat1,mat2),matmul(mat3,mat4))

	enddo
	enddo
	enddo

	!Construct the gluon loop
	do x0 = c0,c1

	do n1=0,l-1			!loop momenta. Optimizable
	do n2=n1,l-1
	do n3=n2,l-1

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

	  call V2(l,nvec_1,-nvec_1,qloop,-qloop,s0,t0,u0,u0,mu,V2mu)


	  trace = Trace2(S_O5_SGS(:,:,t0,x0,s0),V2mu(:,:))
	  g_loop(x0) = g_loop(x0) + multi*trace*D_prop(mu,mu,u0,u0)



	  enddo !mu
	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!momenta
	enddo
	enddo

	enddo	!x0

	  G24f_6b(:) = -C2_R*dimrep*g_loop(:)/(1.0d0*l**3)

	END SUBROUTINE G2_4f_6b



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G2_4f_6c(l,g_i1,g_i2,flav,G24f_6c)

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G24f_6c

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) ::  S_prop_0
	complex(kind=8),dimension(4,4,0:l,c0:c1,0:l) :: S_O5_SGS	
	complex(kind=8),dimension(4,4) :: V2mu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop, tree
	integer :: x0,i
	integer :: s0,t0,u0
	integer :: t0_a,t0_b,u0_a,u0_b
	integer :: mu, nu
	integer :: n1, n2, n3
	complex(kind=8) :: trace
	complex(kind=8),dimension(4,4) :: mat1, mat2, mat3, mat4
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp

	nvec_1(:) = 0
	nvec_2(:) = 0


	G24f_6c(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_O5_SGS(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	f=f1f2(flav(4),flav(1))
	fp=f1f2(flav(2),flav(3))

	!part outside the loop, right trace

	!build S_O5_SGS
	call Fermion_Propagator(l,nvec_2,S_prop_0)

	do t0=0,l
	do x0 = c0,c1
	do s0=0,l

	    mat1 = matmul(S_prop_0(:,:,t0,l-1,flav(2)),Q5p(:,:,fp))
	    mat2 = matmul(S_prop_0(:,:,l-1,x0,flav(3)),gBil(:,:,g_i2))		
	    mat3 = matmul(S_prop_0(:,:,x0,1,flav(4)),Q5(:,:,f))
	    mat4 = matmul(matmul(S_prop_0(:,:,1,x0,flav(1)),gBil(:,:,g_i1)),S_prop_0(:,:,x0,s0,flav(2)))

	    S_O5_SGS(:,:,t0,x0,s0) = matmul(matmul(mat1,mat2),matmul(mat3,mat4))

	enddo
	enddo
	enddo


	!Construct the gluon loop
	do x0 = c0,c1

	do n1=0,l-1			!loop momenta. Optimizable
	do n2=n1,l-1
	do n3=n2,l-1

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

	  call V2(l,nvec_2,-nvec_2,qloop,-qloop,s0,t0,u0,u0,mu,V2mu)


	  trace = Trace2(S_O5_SGS(:,:,t0,x0,s0),V2mu(:,:))
	  g_loop(x0) = g_loop(x0) + multi*trace*D_prop(mu,mu,u0,u0)



	  enddo !mu
	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!momenta
	enddo
	enddo

	enddo	!x0



	  G24f_6c(:) = -C2_R*dimrep*g_loop(:)/(1.0d0*l**3)

	END SUBROUTINE G2_4f_6c



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G2_4f_6d(l,g_i1,g_i2,flav,G24f_6d)

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G24f_6d

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) ::  S_prop_0
	complex(kind=8),dimension(4,4,0:l,c0:c1,0:l) :: S_O5_SGS	
	complex(kind=8),dimension(4,4) :: V2mu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop, tree
	integer :: x0,i
	integer :: s0,t0,u0
	integer :: t0_a,t0_b,u0_a,u0_b
	integer :: mu, nu
	integer :: n1, n2, n3
	complex(kind=8) :: trace
	complex(kind=8),dimension(4,4) :: mat1, mat2, mat3, mat4
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp

	nvec_1(:) = 0
	nvec_2(:) = 0


	G24f_6d(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_O5_SGS(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	f=f1f2(flav(4),flav(1))
	fp=f1f2(flav(2),flav(3))

	!part outside the loop, right trace

	!build S_O5_SGS
	call Fermion_Propagator(l,nvec_2,S_prop_0)

	do t0=0,l
	do x0 = c0,c1
	do s0=0,l

	    mat1 = matmul(S_prop_0(:,:,t0,x0,flav(3)),gBil(:,:,g_i2))
	    mat2 = matmul(S_prop_0(:,:,x0,1,flav(4)),Q5(:,:,f))		
	    mat3 = matmul(S_prop_0(:,:,1,x0,flav(1)),gBil(:,:,g_i1))
	    mat4 = matmul(matmul(S_prop_0(:,:,x0,1,flav(2)),Q5p(:,:,fp)),S_prop_0(:,:,l-1,s0,flav(3)))

	    S_O5_SGS(:,:,t0,x0,s0) = matmul(matmul(mat1,mat2),matmul(mat3,mat4))

	enddo
	enddo
	enddo


	!Construct the gluon loop
	do x0 = c0,c1

	do n1=0,l-1			!loop momenta. Optimizable
	do n2=n1,l-1
	do n3=n2,l-1

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

	  call V2(l,nvec_2,-nvec_2,qloop,-qloop,s0,t0,u0,u0,mu,V2mu)


	  trace = Trace2(S_O5_SGS(:,:,t0,x0,s0),V2mu(:,:))
	  g_loop(x0) = g_loop(x0) + multi*trace*D_prop(mu,mu,u0,u0)



	  enddo !mu
	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!momenta
	enddo
	enddo

	enddo	!x0


	  G24f_6d(:) = -C2_R*dimrep*g_loop(:)/(1.0d0*l**3)

	END SUBROUTINE G2_4f_6d



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G2_4f_7a(l,g_i1,g_i2,flav,G24f_7a)

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G24f_7a

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,0:l,0:l) :: S_O5_S	
	complex(kind=8),dimension(4,4,c0:c1) :: G_S_O_S_G
	complex(kind=8),dimension(4,4) :: V1mu, V1nu
	complex(kind=8),dimension(4,4) :: mat1, mat2	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop, tree
	integer :: x0,i
	integer :: s0,s0p,t0,t0p,u0,u0p
	integer :: t0_a,t0_b,t0p_a,t0p_b,u0_a,u0_b,u0p_a,u0p_b
	integer :: mu, nu
	integer :: n1, n2, n3
	complex(kind=8) :: trace
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp

	nvec_1(:) = 0
	nvec_2(:) = 0

	G24f_7a(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_O5_S(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	G_S_O_S_G(:,:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	f=f1f2(flav(4),flav(1))
	fp=f1f2(flav(1),flav(3))

	!part outside the loop

	!build S_O5_S
	call Fermion_Propagator(l,nvec_1,S_prop_0)
	do x0 = c0,c1
	   mat1(:,:) = matmul(gBil(:,:,g_i1),S_prop_0(:,:,x0,l-1,flav(2)))
	   mat2(:,:) = matmul(S_prop_0(:,:,l-1,x0,flav(3)),gBil(:,:,g_i2))
	   G_S_O_S_G(:,:,x0) = matmul(mat1(:,:),matmul(Q5p(:,:,fp),mat2(:,:)))
	enddo

	do s0=0,l
	do t0=0,l
	    S_O5_S(:,:,t0,s0) = matmul(S_prop_0(:,:,t0,1,flav(4)),matmul(Q5(:,:,f),S_prop_0(:,:,1,s0,flav(1))))
	enddo	
	enddo

	!Construct the gluon loop
	do x0 = c0,c1

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
	   call Fermion_Propagator(l,nvec_1+qloop,S_prop)

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
	  do nu=0,3

	  call V1(l,nvec_1,-nvec_1-qloop,qloop,s0,t0,u0,mu,V1mu)
	  call V1(l,nvec_1+qloop,-nvec_1,-qloop,s0p,t0p,u0p,nu,V1nu)

	  trace = Trace6(S_O5_S(:,:,t0p,s0),V1mu(:,:),S_prop(:,:,t0,x0,flav(1)),G_S_O_S_G(:,:,x0),S_prop(:,:,x0,s0p,flav(4)),V1nu(:,:))
	  g_loop(x0) = g_loop(x0) + multi*trace*D_prop(nu,mu,u0p,u0)

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

	enddo	!x0

	  G24f_7a(:) = C2_R*dimrep*g_loop(:)/(1.0d0*l**3)

	END SUBROUTINE G2_4f_7a



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G2_4f_7b(l,g_i1,g_i2,flav,G24f_7b)

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G24f_7b

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,0:l,0:l) :: S_O5_S	
	complex(kind=8),dimension(4,4,c0:c1) :: G_S_O_S_G
	complex(kind=8),dimension(4,4) :: V1mu, V1nu	
	complex(kind=8),dimension(4,4) :: mat1, mat2
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop, tree
	integer :: x0,i
	integer :: s0,s0p,t0,t0p,u0,u0p
	integer :: t0_a,t0_b,t0p_a,t0p_b,u0_a,u0_b,u0p_a,u0p_b
	integer :: mu, nu
	integer :: n1, n2, n3
	complex(kind=8) :: trace
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp

	nvec_1(:) = 0
	nvec_2(:) = 0

	G24f_7b(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_O5_S(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	G_S_O_S_G(:,:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	f=f1f2(flav(4),flav(1))
	fp=f1f2(flav(2),flav(3))

	!part outside the loop

	!build S_O5_S
	call Fermion_Propagator(l,nvec_1,S_prop_0)
	do x0 = c0,c1
	   mat1(:,:) = matmul(gBil(:,:,g_i2),S_prop_0(:,:,x0,1,flav(4)))
	   mat2(:,:) = matmul(S_prop_0(:,:,1,x0,flav(1)),gBil(:,:,g_i1))
	   G_S_O_S_G(:,:,x0) = matmul(mat1(:,:),matmul(Q5(:,:,f),mat2(:,:)))
	enddo

	do s0=0,l
	do t0=0,l
	    S_O5_S(:,:,t0,s0) = matmul(S_prop_0(:,:,t0,l-1,flav(2)),matmul(Q5p(:,:,fp),S_prop_0(:,:,l-1,s0,flav(3))))
	enddo	
	enddo

	!Construct the gluon loop
	do x0 = c0,c1

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
	   call Fermion_Propagator(l,nvec_1+qloop,S_prop)

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
	  do nu=0,3

	  call V1(l,nvec_1,-nvec_1-qloop,qloop,s0,t0,u0,mu,V1mu)
	  call V1(l,nvec_1+qloop,-nvec_1,-qloop,s0p,t0p,u0p,nu,V1nu)

	  trace = Trace6(G_S_O_S_G(:,:,x0),S_prop(:,:,x0,s0,flav(2)),V1mu(:,:),S_O5_S(:,:,t0p,s0),V1nu(:,:),S_prop(:,:,t0p,x0,flav(3)))
	  g_loop(x0) = g_loop(x0) + multi*trace*D_prop(mu,nu,u0,u0p)

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

	enddo	!x0

	  G24f_7b(:) = C2_R*dimrep*g_loop(:)/(1.0d0*l**3)

	END SUBROUTINE G2_4f_7b



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE G2_4f_8a(l,g_i1,g_i2,flav,G24f_8a)

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G24f_8a
	
	integer :: n1,n2,n3
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop_1, S_prop_2
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop!,K_prop
	complex(kind=8),dimension(c0:c1) :: loop
	complex(kind=8),dimension(4,4) :: mat1, mat2, mat3 ,mat4
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f1,f2,x0

	nvec_1(:) = 0
	nvec_2(:) = 0

	!initialization
	G24f_8a(:) = dcmplx(0.0d0,0.0d0)
	S_prop_1(:,:,:,:,:) = dcmplx(0.0d0)
	S_prop_2(:,:,:,:,:) = dcmplx(0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	loop(:) = dcmplx(0.0d0,0.0d0)

	f1 = f1f2(flav(4),flav(1))
	f2 = f1f2(flav(2),flav(3))
        call Fermion_Propagator(l,nvec_2,S_Prop_2)
	
	!Construct the trace in the loop
	do n1=0,l-1			!loop momenta. Optimizable
	do n2=n1,l-1
	do n3=n2,l-1  

	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3

	   multi = 1.0d0*multiplicity(n1,n2,n3)

  	   D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	   S_prop_1(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	   S_prop_2(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	
	   call Fermion_Propagator(l,nvec_1+qloop,S_Prop_1)
  	   call Gluon_propagator(l,qloop,D_prop)	

	   do x0=c0,c1

	     mat1 = matmul(Q5(:,:,f1),S_prop_2(:,:,1,x0,flav(1)))
	     mat2 = matmul(gBil(:,:,g_i1),S_prop_2(:,:,x0,l-1,flav(2)))
	     mat3 = matmul(Q5p(:,:,f2),S_prop_1(:,:,l-1,x0,flav(3)))
	     mat4 = matmul(gBil(:,:,g_i2),S_prop_1(:,:,x0,1,flav(4)))

	     loop(x0) = loop(x0) + dimrep*D_prop(0,0,0,l-1)*multi*Trace4(mat1,mat2,mat3,mat4)

	   enddo

	enddo
	enddo
	enddo

	G24f_8a(:) = C2_R*loop(:)/((l*1.0d0)**3)

	END SUBROUTINE G2_4f_8a

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE G2_4f_8b(l,g_i1,g_i2,flav,G24f_8b)

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G24f_8b
	
	integer :: n1,n2,n3
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop_1, S_prop_2
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop!,K_prop
	complex(kind=8),dimension(c0:c1) :: loop
	complex(kind=8),dimension(4,4) :: mat1, mat2, mat3 ,mat4
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f1,f2, x0

	nvec_1(:) = 0
	nvec_2(:) = 0

	!initialization
	G24f_8b(:) = dcmplx(0.0d0,0.0d0)
	S_prop_1(:,:,:,:,:) = dcmplx(0.0d0)
	S_prop_2(:,:,:,:,:) = dcmplx(0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	loop(:) = dcmplx(0.0d0,0.0d0)

	f1 = f1f2(flav(4),flav(1))
	f2 = f1f2(flav(2),flav(3))
	call Fermion_Propagator(l,nvec_1,S_Prop_1)

	!Construct the trace in the loop
	do n1=0,l-1			!loop momenta. Optimizable
	do n2=n1,l-1
	do n3=n2,l-1  

	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3

	   multi = 1.0d0*multiplicity(n1,n2,n3)

  	   D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	   S_prop_1(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	   S_prop_2(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

	   call Fermion_Propagator(l,nvec_2+qloop,S_Prop_2)
  	   call Gluon_propagator(l,qloop,D_prop)	

	   do x0=c0,c1

	     mat1 = matmul(Q5(:,:,f1),S_prop_2(:,:,1,x0,flav(1)))
	     mat2 = matmul(gBil(:,:,g_i1),S_prop_2(:,:,x0,l-1,flav(2)))
	     mat3 = matmul(Q5p(:,:,f2),S_prop_1(:,:,l-1,x0,flav(3)))
	     mat4 = matmul(gBil(:,:,g_i2),S_prop_1(:,:,x0,1,flav(4)))

	     loop(x0) = loop(x0) + dimrep*D_prop(0,0,l-1,0)*multi*Trace4(mat1,mat2,mat3,mat4)

	   enddo

	enddo
	enddo
	enddo

	G24f_8b(:) = C2_R*G24f_8b(:)/((l*1.0d0)**3)

	END SUBROUTINE G2_4f_8b

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE G2_4f_9a(l,g_i1,g_i2,flav,G24f_9a)

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G24f_9a
	
	integer :: n1,n2,n3
	integer,dimension(3) :: qloop, mqloop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop_0, S_prop_q, S_prop_mq
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop!,K_prop
	complex(kind=8),dimension(c0:c1) :: loop
	complex(kind=8),dimension(4,4) :: mat1, mat2, mat3 ,mat4
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f1,f2,x0

	nvec_1(:) = 0
	nvec_2(:) = 0

	!initialization
	G24f_9a(:) = dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:) = dcmplx(0.0d0)
	S_prop_q(:,:,:,:,:) = dcmplx(0.0d0)
	S_prop_mq(:,:,:,:,:) = dcmplx(0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	loop(:) = dcmplx(0.0d0,0.0d0)

	f1 = f1f2(flav(4),flav(1))
	f2 = f1f2(flav(2),flav(3))
        call Fermion_Propagator(l,nvec_2,S_Prop_0)
	
	!Construct the trace in the loop
	do n1=0,l-1			!loop momenta. Optimizable
	do n2=n1,l-1
	do n3=n2,l-1  

	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3
	mqloop(1) = m_period(n1,l)
	mqloop(2) = m_period(n2,l)
	mqloop(3) = m_period(n3,l)


	   multi = 1.0d0*multiplicity(n1,n2,n3)

  	   D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	   S_prop_q(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	   S_prop_mq(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	
	   call Fermion_Propagator(l,nvec_1+qloop,S_Prop_q)
	   call Fermion_Propagator(l,nvec_1+mqloop,S_Prop_mq)
  	   call Gluon_propagator(l,qloop,D_prop)	

	   do x0=c0,c1

	     mat1 = matmul(Q5(:,:,f1),S_prop_0(:,:,1,x0,flav(1)))
	     mat2 = matmul(gBil(:,:,g_i1),S_prop_mq(:,:,x0,l-1,flav(2)))
	     mat3 = matmul(Q5p(:,:,f2),S_prop_0(:,:,l-1,x0,flav(3)))
	     mat4 = matmul(gBil(:,:,g_i2),S_prop_q(:,:,x0,1,flav(4)))

	     loop(x0) = loop(x0) + dimrep*D_prop(0,0,0,l-1)*multi*Trace4(mat1,mat2,mat3,mat4)

	   enddo

	enddo
	enddo
	enddo

	G24f_9a(:) = -C2_R*loop(:)/((l*1.0d0)**3)

	END SUBROUTINE G2_4f_9a

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE G2_4f_9b(l,g_i1,g_i2,flav,G24f_9b)

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G24f_9b
	
	integer :: n1,n2,n3
	integer,dimension(3) :: qloop,mqloop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop_0, S_prop_q, S_prop_mq
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop!,K_prop
	complex(kind=8),dimension(c0:c1) :: loop
	complex(kind=8),dimension(4,4) :: mat1, mat2, mat3 ,mat4
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f1,f2, x0

	nvec_1(:) = 0
	nvec_2(:) = 0

	!initialization
	G24f_9b(:) = dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:) = dcmplx(0.0d0)
	S_prop_q(:,:,:,:,:) = dcmplx(0.0d0)
	S_prop_mq(:,:,:,:,:) = dcmplx(0.0d0)
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	loop(:) = dcmplx(0.0d0,0.0d0)

	f1 = f1f2(flav(4),flav(1))
	f2 = f1f2(flav(2),flav(3))
	call Fermion_Propagator(l,nvec_1,S_Prop_0)

	!Construct the trace in the loop
	do n1=0,l-1			!loop momenta. Optimizable
	do n2=n1,l-1
	do n3=n2,l-1  

	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3
	mqloop(1) = m_period(n1,l)
	mqloop(2) = m_period(n2,l)
	mqloop(3) = m_period(n3,l)


	   multi = 1.0d0*multiplicity(n1,n2,n3)

  	   D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	   S_prop_q(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	   S_prop_mq(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

	   call Fermion_Propagator(l,nvec_1+qloop,S_Prop_q)
	   call Fermion_Propagator(l,nvec_1+mqloop,S_Prop_mq)
  	   call Gluon_propagator(l,qloop,D_prop)	

	   do x0=c0,c1

	     mat1 = matmul(Q5(:,:,f1),S_prop_q(:,:,1,x0,flav(1)))
	     mat2 = matmul(gBil(:,:,g_i1),S_prop_0(:,:,x0,l-1,flav(2)))
	     mat3 = matmul(Q5p(:,:,f2),S_prop_mq(:,:,l-1,x0,flav(3)))
	     mat4 = matmul(gBil(:,:,g_i2),S_prop_q(:,:,x0,1,flav(4)))

	     loop(x0) = loop(x0) + dimrep*D_prop(0,0,l-1,0)*multi*Trace4(mat1,mat2,mat3,mat4)

	   enddo

	enddo
	enddo
	enddo

	G24f_9b(:) = -C2_R*G24f_9b(:)/((l*1.0d0)**3)

	END SUBROUTINE G2_4f_9b

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G2_4f_10a(l,g_i1,g_i2,flav,G24f_10a)

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G24f_10a

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,c0:c1,0:l) :: OSGSOS
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop
	integer :: n1, n2, n3, x0, t0, u0, s0, mu,i, x0_t
	integer :: t0_a, t0_b, u0_a, u0_b
	complex(kind=8) :: trace
	complex(kind=8),dimension(4,4) :: mat1, mat2,mat3
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp

	nvec_1(:) = 0
	nvec_2(:) = 0

	G24f_10a(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	OSGSOS(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)

	f=f1f2(flav(4),flav(1))
	fp=f1f2(flav(2),flav(3))

	!part outside the loop, left trace
	call Fermion_Propagator(l,nvec_1,S_prop_0)


	do x0 = c0,c1
	do s0=0,l

	    mat1 = matmul(Q5(:,:,f),S_prop_0(:,:,1,x0,flav(1)))
	    mat2 = matmul(gBil(:,:,g_i1),S_prop_0(:,:,x0,l-1,flav(2)))
	    mat3 = matmul(Q5p(:,:,fp),S_prop_0(:,:,l-1,s0,flav(3)))

	    OSGSOS(:,:,x0,s0) = matmul(matmul(mat1,mat2),mat3)

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
	call Fermion_Propagator(l,nvec_1+qloop,S_prop)

	do mu = 0,3

	do s0 = 0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0 = t0_a,t0_b
	do u0 = u0_a,u0_b

	call V1(l,nvec_1,-nvec_1-qloop,qloop,s0,t0,u0,mu,V1mu)
	
	do x0 = c0,c1

	  trace = Trace5(OSGSOS(:,:,x0,s0),V1mu(:,:),S_prop(:,:,t0,x0,flav(3)),gBil(:,:,g_i2),S_prop(:,:,x0,1,flav(4)))
	  g_loop(x0) = g_loop(x0) + D_prop(0,mu,0,u0)*trace*multi

	enddo

	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu

	enddo	!momenta
	enddo
	enddo

	G24f_10a(:) = C2_R*dimrep*g_loop(:)/(1.0d0*l**3)

	END SUBROUTINE G2_4f_10a



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G2_4f_10b(l,g_i1,g_i2,flav,G24f_10b)

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G24f_10b

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,0:l,c0:c1) :: OSGSOS
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop
	integer :: n1, n2, n3, x0, t0, u0, s0, mu,i, x0_t
	integer :: t0_a, t0_b, u0_a, u0_b
	complex(kind=8) :: trace
	complex(kind=8),dimension(4,4) :: mat1, mat2,mat3
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp

	nvec_1(:) = 0
	nvec_2(:) = 0

	G24f_10b(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	OSGSOS(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)

	f=f1f2(flav(4),flav(1))
	fp=f1f2(flav(2),flav(3))

	!part outside the loop, left trace
	call Fermion_Propagator(l,nvec_1,S_prop_0)

	do t0=0,l
	do x0 = c0,c1

	    mat1 = matmul(S_prop_0(:,:,t0,l-1,flav(2)),Q5p(:,:,fp))
	    mat2 = matmul(S_prop_0(:,:,l-1,x0,flav(3)),gBil(:,:,g_i2))
	    mat3 = matmul(S_prop_0(:,:,x0,1,flav(4)),Q5(:,:,f))

	    OSGSOS(:,:,t0,x0) = matmul(matmul(mat1,mat2),mat3)

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
	call Fermion_Propagator(l,nvec_1+qloop,S_prop)

	do mu = 0,3

	do s0 = 0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0 = t0_a,t0_b
	do u0 = u0_a,u0_b

	call V1(l,nvec_1+qloop,-nvec_1,-qloop,s0,t0,u0,mu,V1mu)
	
	do x0 = c0,c1

	  trace = Trace5(OSGSOS(:,:,t0,x0),S_prop(:,:,1,x0,flav(1)),gBil(:,:,g_i1),S_prop(:,:,x0,s0,flav(2)),V1mu(:,:))
	  g_loop(x0) = g_loop(x0) + D_prop(mu,0,u0,0)*trace*multi

	enddo

	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu

	enddo	!momenta
	enddo
	enddo

	G24f_10b(:) = -C2_R*dimrep*g_loop(:)/(1.0d0*l**3)

	END SUBROUTINE G2_4f_10b



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G2_4f_10c(l,g_i1,g_i2,flav,G24f_10c)

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G24f_10c

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,c0:c1,0:l) :: OSGSOS
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop
	integer :: n1, n2, n3, x0, t0, u0, s0, mu,i, x0_t
	integer :: t0_a, t0_b, u0_a, u0_b
	complex(kind=8) :: trace
	complex(kind=8),dimension(4,4) :: mat1, mat2,mat3
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp

	nvec_1(:) = 0
	nvec_2(:) = 0

	G24f_10c(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	OSGSOS(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)

	f=f1f2(flav(4),flav(1))
	fp=f1f2(flav(2),flav(3))

	!part outside the loop, left trace
	call Fermion_Propagator(l,nvec_1,S_prop_0)


	do x0 = c0,c1
	do s0=0,l

	    mat1 = matmul(Q5p(:,:,fp),S_prop_0(:,:,l-1,x0,flav(3)))
	    mat2 = matmul(gBil(:,:,g_i2),S_prop_0(:,:,x0,1,flav(4)))
	    mat3 = matmul(Q5(:,:,f),S_prop_0(:,:,1,s0,flav(1)))

	    OSGSOS(:,:,x0,s0) = matmul(matmul(mat1,mat2),mat3)

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
	call Fermion_Propagator(l,nvec_1+qloop,S_prop)

	do mu = 0,3

	do s0 = 0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0 = t0_a,t0_b
	do u0 = u0_a,u0_b

	call V1(l,nvec_1,-nvec_1-qloop,qloop,s0,t0,u0,mu,V1mu)
	
	do x0 = c0,c1

	  trace = Trace5(OSGSOS(:,:,x0,s0),V1mu(:,:),S_prop(:,:,t0,x0,flav(1)),gBil(:,:,g_i1),S_prop(:,:,x0,l-1,flav(2)))
	  g_loop(x0) = g_loop(x0) + D_prop(0,mu,0,u0)*trace*multi

	enddo

	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu

	enddo	!momenta
	enddo
	enddo

	G24f_10c(:) = -C2_R*dimrep*g_loop(:)/(1.0d0*l**3)

	END SUBROUTINE G2_4f_10c



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G2_4f_10d(l,g_i1,g_i2,flav,G24f_10d)

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G24f_10d

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,0:l,c0:c1) :: OSGSOS
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop
	integer :: n1, n2, n3, x0, t0, u0, s0, mu,i, x0_t
	integer :: t0_a, t0_b, u0_a, u0_b
	complex(kind=8) :: trace
	complex(kind=8),dimension(4,4) :: mat1, mat2,mat3
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp

	nvec_1(:) = 0
	nvec_2(:) = 0

	G24f_10d(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	OSGSOS(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)

	f=f1f2(flav(4),flav(1))
	fp=f1f2(flav(2),flav(3))

	!part outside the loop, left trace
	call Fermion_Propagator(l,nvec_1,S_prop_0)

	do t0=0,l
	do x0 = c0,c1

	    mat1 = matmul(S_prop_0(:,:,t0,1,flav(4)),Q5(:,:,f))
	    mat2 = matmul(S_prop_0(:,:,1,x0,flav(1)),gBil(:,:,g_i1))
	    mat3 = matmul(S_prop_0(:,:,x0,l-1,flav(2)),Q5p(:,:,fp))

	    OSGSOS(:,:,t0,x0) = matmul(matmul(mat1,mat2),mat3)

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
	call Fermion_Propagator(l,nvec_1+qloop,S_prop)

	do mu = 0,3

	do s0 = 0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0 = t0_a,t0_b
	do u0 = u0_a,u0_b

	call V1(l,nvec_1+qloop,-nvec_1,-qloop,s0,t0,u0,mu,V1mu)
	
	do x0 = c0,c1

	  trace = Trace5(OSGSOS(:,:,t0,x0),S_prop(:,:,l-1,x0,flav(3)),gBil(:,:,g_i2),S_prop(:,:,x0,s0,flav(4)),V1mu(:,:))
	  g_loop(x0) = g_loop(x0) + D_prop(mu,0,u0,0)*trace*multi

	enddo

	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu

	enddo	!momenta
	enddo
	enddo

	G24f_10d(:) = C2_R*dimrep*g_loop(:)/(1.0d0*l**3)

	END SUBROUTINE G2_4f_10d



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G2_4f_11a(l,g_i1,g_i2,flav,G24f_11a)

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G24f_11a

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop_q,S_prop_mq, S_prop_0
	complex(kind=8),dimension(4,4,0:l,c0:c1) :: SOSG
	complex(kind=8),dimension(4,4,c0:c1) :: OSG
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop, mqloop
	complex(kind=8),dimension(c0:c1) :: g_loop
	integer :: n1, n2, n3, x0, t0, u0, s0, mu,i, x0_t
	integer :: t0_a, t0_b, u0_a, u0_b
	complex(kind=8) :: trace
	complex(kind=8),dimension(4,4) :: mat1, mat2,mat3
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp

	nvec_1(:) = 0
	nvec_2(:) = 0

	G24f_11a(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_q(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_mq(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	SOSG(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	OSG(:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)

	f=f1f2(flav(4),flav(1))
	fp=f1f2(flav(2),flav(3))

	!part outside the loop, left trace
	call Fermion_Propagator(l,nvec_1,S_prop_0)


	do t0=0,l
	do x0 = c0,c1

	    mat1 = matmul(S_prop_0(:,:,t0,l-1,flav(2)),Q5p(:,:,fp))
	    mat2 = matmul(S_prop_0(:,:,l-1,x0,flav(3)),gBil(:,:,g_i2))
	    SOSG(:,:,t0,x0) = matmul(mat1,mat2)

	enddo
	enddo

	do x0 = c0,c1
	    OSG(:,:,x0) = matmul(Q5(:,:,f),matmul(S_prop_0(:,:,1,x0,flav(1)),gBil(:,:,g_i1)))
	enddo


	!Construct the gluon loop
	do n1=0,l-1			!loop momenta. Optimizable
	do n2=n1,l-1
	do n3=n2,l-1

	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3
	mqloop(1) = m_period(n1,l)
	mqloop(2) = m_period(n2,l)
	mqloop(3) = m_period(n3,l)


	multi = 1.0d0*multiplicity(n1,n2,n3)

  	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	S_prop_q(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	S_prop_mq(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

  	call Gluon_propagator(l,qloop,D_prop)	
	call Fermion_Propagator(l,nvec_1+qloop,S_prop_q)
	call Fermion_Propagator(l,nvec_1+mqloop,S_prop_mq)

	do mu = 0,3

	do s0 = 0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0 = t0_a,t0_b
	do u0 = u0_a,u0_b

	call V1(l,nvec_1-qloop,-nvec_1,qloop,s0,t0,u0,mu,V1mu)
	
	do x0 = c0,c1

	  trace = Trace5(OSG(:,:,x0),S_prop_mq(:,:,x0,s0,flav(2)),V1mu(:,:),SOSG(:,:,t0,x0),S_prop_q(:,:,x0,l-1,flav(4)))
	  g_loop(x0) = g_loop(x0) + D_prop(0,mu,0,u0)*trace*multi

	enddo

	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu

	enddo	!momenta
	enddo
	enddo

	G24f_11a(:) = C2_R*dimrep*g_loop(:)/(1.0d0*l**3)

	END SUBROUTINE G2_4f_11a



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G2_4f_11b(l,g_i1,g_i2,flav,G24f_11b)

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G24f_11b

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop_q,S_prop_mq, S_prop_0
	complex(kind=8),dimension(4,4,c0:c1,0:l) :: SOSG
	complex(kind=8),dimension(4,4,c0:c1) :: OSG
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop, mqloop
	complex(kind=8),dimension(c0:c1) :: g_loop
	integer :: n1, n2, n3, x0, t0, u0, s0, mu,i, x0_t
	integer :: t0_a, t0_b, u0_a, u0_b
	complex(kind=8) :: trace
	complex(kind=8),dimension(4,4) :: mat1, mat2,mat3
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp

	nvec_1(:) = 0
	nvec_2(:) = 0

	G24f_11b(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_q(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_mq(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	SOSG(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	OSG(:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)

	f=f1f2(flav(4),flav(1))
	fp=f1f2(flav(2),flav(3))

	!part outside the loop, left trace
	call Fermion_Propagator(l,nvec_1,S_prop_0)


	do x0 = c0,c1
	do s0=0,l

	    mat1 = matmul(gBil(:,:,g_i1),S_prop_0(:,:,x0,l-1,flav(2)))
	    mat2 = matmul(Q5p(:,:,fp),S_prop_0(:,:,l-1,s0,flav(3)))
	    SOSG(:,:,x0,s0) = matmul(mat1,mat2)

	enddo
	enddo

	do x0 = c0,c1
	    OSG(:,:,x0) = matmul(gBil(:,:,g_i2),matmul(S_prop_0(:,:,x0,1,flav(4)),Q5(:,:,f)))
	enddo


	!Construct the gluon loop
	do n1=0,l-1			!loop momenta. Optimizable
	do n2=n1,l-1
	do n3=n2,l-1

	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3
	mqloop(1) = m_period(n1,l)
	mqloop(2) = m_period(n2,l)
	mqloop(3) = m_period(n3,l)


	multi = 1.0d0*multiplicity(n1,n2,n3)

  	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	S_prop_q(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	S_prop_mq(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

  	call Gluon_propagator(l,qloop,D_prop)	
	call Fermion_Propagator(l,nvec_1+qloop,S_prop_q)
	call Fermion_Propagator(l,nvec_1+mqloop,S_prop_mq)

	do mu = 0,3

	do s0 = 0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0 = t0_a,t0_b
	do u0 = u0_a,u0_b

	call V1(l,nvec_1,-nvec_1+qloop,-qloop,s0,t0,u0,mu,V1mu)
	
	do x0 = c0,c1

	  trace = Trace5(OSG(:,:,x0),S_prop_q(:,:,1,x0,flav(1)),SOSG(:,:,x0,s0),V1mu(:,:),S_prop_mq(:,:,t0,x0,flav(3)))
	  g_loop(x0) = g_loop(x0) + D_prop(mu,0,u0,0)*trace*multi

	enddo

	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu

	enddo	!momenta
	enddo
	enddo

	G24f_11b(:) = -C2_R*dimrep*g_loop(:)/(1.0d0*l**3)

	END SUBROUTINE G2_4f_11b



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G2_4f_11c(l,g_i1,g_i2,flav,G24f_11c)

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G24f_11c

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop_q,S_prop_mq, S_prop_0
	complex(kind=8),dimension(4,4,0:l,c0:c1) :: SOSG
	complex(kind=8),dimension(4,4,c0:c1) :: OSG
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop, mqloop
	complex(kind=8),dimension(c0:c1) :: g_loop
	integer :: n1, n2, n3, x0, t0, u0, s0, mu,i, x0_t
	integer :: t0_a, t0_b, u0_a, u0_b
	complex(kind=8) :: trace
	complex(kind=8),dimension(4,4) :: mat1, mat2,mat3
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp

	nvec_1(:) = 0
	nvec_2(:) = 0

	G24f_11c(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_q(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_mq(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	SOSG(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	OSG(:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)

	f=f1f2(flav(4),flav(1))
	fp=f1f2(flav(2),flav(3))

	!part outside the loop, left trace
	call Fermion_Propagator(l,nvec_1,S_prop_0)


	do t0=0,l
	do x0 = c0,c1

	    mat1 = matmul(S_prop_0(:,:,t0,1,flav(4)),Q5(:,:,f))
	    mat2 = matmul(S_prop_0(:,:,1,x0,flav(1)),gBil(:,:,g_i1))
	    SOSG(:,:,t0,x0) = matmul(mat1,mat2)

	enddo
	enddo

	do x0 = c0,c1
	    OSG(:,:,x0) = matmul(Q5p(:,:,fp),matmul(S_prop_0(:,:,l-1,x0,flav(3)),gBil(:,:,g_i2)))
	enddo


	!Construct the gluon loop
	do n1=0,l-1			!loop momenta. Optimizable
	do n2=n1,l-1
	do n3=n2,l-1

	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3
	mqloop(1) = m_period(n1,l)
	mqloop(2) = m_period(n2,l)
	mqloop(3) = m_period(n3,l)


	multi = 1.0d0*multiplicity(n1,n2,n3)

  	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	S_prop_q(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	S_prop_mq(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

  	call Gluon_propagator(l,qloop,D_prop)	
	call Fermion_Propagator(l,nvec_1+qloop,S_prop_q)
	call Fermion_Propagator(l,nvec_1+mqloop,S_prop_mq)

	do mu = 0,3

	do s0 = 0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0 = t0_a,t0_b
	do u0 = u0_a,u0_b

	call V1(l,nvec_1+qloop,-nvec_1,-qloop,s0,t0,u0,mu,V1mu)
	
	do x0 = c0,c1

	  trace = Trace5(OSG(:,:,x0),S_prop_q(:,:,x0,s0,flav(4)),V1mu(:,:),SOSG(:,:,t0,x0),S_prop_q(:,:,x0,l-1,flav(2)))
	  g_loop(x0) = g_loop(x0) + D_prop(mu,0,u0,0)*trace*multi

	enddo

	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu

	enddo	!momenta
	enddo
	enddo

	G24f_11c(:) = -C2_R*dimrep*g_loop(:)/(1.0d0*l**3)

	END SUBROUTINE G2_4f_11c



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G2_4f_11d(l,g_i1,g_i2,flav,G24f_11d)

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G24f_11d

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop_q,S_prop_mq, S_prop_0
	complex(kind=8),dimension(4,4,c0:c1,0:l) :: SOSG
	complex(kind=8),dimension(4,4,c0:c1) :: OSG
	complex(kind=8),dimension(4,4) :: V1mu	
	integer,dimension(3) :: qloop, mqloop
	complex(kind=8),dimension(c0:c1) :: g_loop
	integer :: n1, n2, n3, x0, t0, u0, s0, mu,i, x0_t
	integer :: t0_a, t0_b, u0_a, u0_b
	complex(kind=8) :: trace
	complex(kind=8),dimension(4,4) :: mat1, mat2,mat3
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp

	nvec_1(:) = 0
	nvec_2(:) = 0

	G24f_11d(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_q(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_mq(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	SOSG(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	OSG(:,:,:)=dcmplx(0.0d0,0.0d0)
	V1mu(:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)

	f=f1f2(flav(4),flav(1))
	fp=f1f2(flav(2),flav(3))

	!part outside the loop, left trace
	call Fermion_Propagator(l,nvec_1,S_prop_0)


	do x0 = c0,c1
	do s0=0,l

	    mat1 = matmul(gBil(:,:,g_i2),S_prop_0(:,:,x0,1,flav(4)))
	    mat2 = matmul(Q5(:,:,f),S_prop_0(:,:,1,s0,flav(1)))
	    SOSG(:,:,x0,s0) = matmul(mat1,mat2)

	enddo
	enddo

	do x0 = c0,c1
	    OSG(:,:,x0) = matmul(gBil(:,:,g_i1),matmul(S_prop_0(:,:,x0,1,flav(2)),Q5p(:,:,fp)))
	enddo


	!Construct the gluon loop
	do n1=0,l-1			!loop momenta. Optimizable
	do n2=n1,l-1
	do n3=n2,l-1

	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3
	mqloop(1) = m_period(n1,l)
	mqloop(2) = m_period(n2,l)
	mqloop(3) = m_period(n3,l)


	multi = 1.0d0*multiplicity(n1,n2,n3)

  	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	S_prop_q(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	S_prop_mq(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

  	call Gluon_propagator(l,qloop,D_prop)	
	call Fermion_Propagator(l,nvec_1+qloop,S_prop_q)
	call Fermion_Propagator(l,nvec_1+mqloop,S_prop_mq)

	do mu = 0,3

	do s0 = 0,l
	call loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
	do t0 = t0_a,t0_b
	do u0 = u0_a,u0_b

	call V1(l,nvec_1,-nvec_1-qloop,qloop,s0,t0,u0,mu,V1mu)
	
	do x0 = c0,c1

	  trace = Trace5(OSG(:,:,x0),S_prop_mq(:,:,l-1,x0,flav(3)),SOSG(:,:,x0,s0),V1mu(:,:),S_prop_mq(:,:,t0,x0,flav(1)))
	  g_loop(x0) = g_loop(x0) + D_prop(0,mu,0,u0)*trace*multi

	enddo

	enddo	!u0
	enddo	!t0
	enddo	!s0

	enddo	!mu

	enddo	!momenta
	enddo
	enddo

	G24f_11d(:) = C2_R*dimrep*g_loop(:)/(1.0d0*l**3)

	END SUBROUTINE G2_4f_11d



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G2_4f_12a(l,g_i1,g_i2,flav,G24f_12a)

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G24f_12a

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,0:l,c0:c1,0:l) :: S_O5_S	
	complex(kind=8),dimension(4,4) :: V1mu, V1nu
	complex(kind=8),dimension(4,4) :: mat1, mat2, mat3	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop, tree
	integer :: x0,i
	integer :: s0,s0p,t0,t0p,u0,u0p
	integer :: t0_a,t0_b,t0p_a,t0p_b,u0_a,u0_b,u0p_a,u0p_b
	integer :: mu, nu
	integer :: n1, n2, n3
	complex(kind=8) :: trace
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp

	nvec_1(:) = 0
	nvec_2(:) = 0

	G24f_12a(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_O5_S(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	f=f1f2(flav(4),flav(1))
	fp=f1f2(flav(2),flav(3))

	!part outside the loop

	!build S_O5_S
	call Fermion_Propagator(l,nvec_1,S_prop_0)


	do s0=0,l
	do t0=0,l
	do x0=c0,c1

	    mat1(:,:) = matmul(S_prop_0(:,:,t0,1,flav(4)),Q5(:,:,f))
	    mat2(:,:) = matmul(S_prop_0(:,:,1,x0,flav(1)),gBil(:,:,g_i1))
	    mat3(:,:) = matmul(S_prop_0(:,:,x0,l-1,flav(2)),Q5p(:,:,fp))

	    S_O5_S(:,:,t0,x0,s0) = matmul(mat1(:,:),matmul(mat2(:,:),matmul(mat3(:,:),S_prop_0(:,:,l-1,s0,flav(3)))))

	enddo
	enddo	
	enddo

	!Construct the gluon loop
	do x0 = c0,c1

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
	   call Fermion_Propagator(l,nvec_1+qloop,S_prop)

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
	  do nu=0,3

	  call V1(l,nvec_1,-nvec_1-qloop,qloop,s0,t0,u0,mu,V1mu)
	  call V1(l,nvec_1+qloop,-nvec_1,-qloop,s0p,t0p,u0p,nu,V1nu)

	  trace = Trace6(S_O5_S(:,:,t0,x0,s0),V1mu(:,:),S_prop(:,:,t0,x0,flav(3)),gBil(:,:,g_i2),S_prop(:,:,x0,s0p,flav(4)),V1nu(:,:))
	  g_loop(x0) = g_loop(x0) + multi*trace*D_prop(nu,mu,u0p,u0)

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

	enddo	!x0

	  G24f_12a(:) = C2_R*dimrep*g_loop(:)/(1.0d0*l**3)

	END SUBROUTINE G2_4f_12a



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G2_4f_12b(l,g_i1,g_i2,flav,G24f_12b)

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G24f_12b

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop, S_prop_0
	complex(kind=8),dimension(4,4,0:l,c0:c1,0:l) :: S_O5_S	
	complex(kind=8),dimension(4,4) :: V1mu, V1nu
	complex(kind=8),dimension(4,4) :: mat1, mat2, mat3	
	integer,dimension(3) :: qloop
	complex(kind=8),dimension(c0:c1) :: g_loop, tree
	integer :: x0,i
	integer :: s0,s0p,t0,t0p,u0,u0p
	integer :: t0_a,t0_b,t0p_a,t0p_b,u0_a,u0_b,u0p_a,u0p_b
	integer :: mu, nu
	integer :: n1, n2, n3
	complex(kind=8) :: trace
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp

	nvec_1(:) = 0
	nvec_2(:) = 0

	G24f_12b(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_O5_S(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	f=f1f2(flav(4),flav(1))
	fp=f1f2(flav(2),flav(3))

	!part outside the loop

	!build S_O5_S
	call Fermion_Propagator(l,nvec_1,S_prop_0)

	do s0=0,l
	do t0=0,l
	do x0=c0,c1

	    mat1(:,:) = matmul(S_prop_0(:,:,t0,l-1,flav(2)),Q5p(:,:,fp))
	    mat2(:,:) = matmul(S_prop_0(:,:,l-1,x0,flav(3)),gBil(:,:,g_i2))
	    mat3(:,:) = matmul(S_prop_0(:,:,x0,1,flav(4)),Q5(:,:,f))

	    S_O5_S(:,:,t0,x0,s0) = matmul(mat1(:,:),matmul(mat2(:,:),matmul(mat3(:,:),S_prop_0(:,:,1,s0,flav(1)))))

	enddo
	enddo	
	enddo

	!Construct the gluon loop
	do x0 = c0,c1

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
	   call Fermion_Propagator(l,nvec_1+qloop,S_prop)

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
	  do nu=0,3

	  call V1(l,nvec_1,-nvec_1-qloop,qloop,s0,t0,u0,mu,V1mu)
	  call V1(l,nvec_1+qloop,-nvec_1,-qloop,s0p,t0p,u0p,nu,V1nu)

	  trace = Trace6(S_O5_S(:,:,t0,x0,s0),V1mu(:,:),S_prop(:,:,t0,x0,flav(1)),gBil(:,:,g_i1),S_prop(:,:,x0,s0p,flav(2)),V1nu(:,:))
	  g_loop(x0) = g_loop(x0) + multi*trace*D_prop(nu,mu,u0p,u0)

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

	enddo	!x0

	  G24f_12b(:) = C2_R*dimrep*g_loop(:)/(1.0d0*l**3)

	END SUBROUTINE G2_4f_12b



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G2_4f_13a(l,g_i1,g_i2,flav,G24f_13a)

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G24f_13a

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop_q, S_prop_mq, S_prop_0
	complex(kind=8),dimension(4,4,c0:c1,0:l) :: SOSG1, SOSG2
	complex(kind=8),dimension(4,4) :: V1mu, V1nu
	complex(kind=8),dimension(4,4) :: mat1, mat2, mat3	
	integer,dimension(3) :: qloop,mqloop
	complex(kind=8),dimension(c0:c1) :: g_loop, tree
	integer :: x0,i
	integer :: s0,s0p,t0,t0p,u0,u0p
	integer :: t0_a,t0_b,t0p_a,t0p_b,u0_a,u0_b,u0p_a,u0p_b
	integer :: mu, nu
	integer :: n1, n2, n3
	complex(kind=8) :: trace
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp

	nvec_1(:) = 0
	nvec_2(:) = 0

	G24f_13a(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_q(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_mq(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	SOSG1(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	SOSG2(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	f=f1f2(flav(4),flav(1))
	fp=f1f2(flav(2),flav(3))

	!part outside the loop

	!build S_O5_S
	call Fermion_Propagator(l,nvec_1,S_prop_0)

	do t0=0,l
	do x0=c0,c1

	    mat1(:,:) = matmul(S_prop_0(:,:,t0,1,flav(4)),Q5(:,:,f))
	    mat2(:,:) = matmul(S_prop_0(:,:,1,x0,flav(1)),gBil(:,:,g_i1))
	    SOSG1(:,:,t0,x0) = matmul(mat1(:,:),mat2(:,:))

	    mat1(:,:) = matmul(S_prop_0(:,:,t0,l-1,flav(2)),Q5p(:,:,fp))
	    mat2(:,:) = matmul(S_prop_0(:,:,l-1,x0,flav(3)),gBil(:,:,g_i2))
	    SOSG2(:,:,t0,x0) = matmul(mat1(:,:),mat2(:,:))

	enddo
	enddo	

	!Construct the gluon loop
	do x0 = c0,c1

	do n1=0,l-1			!loop momenta. Optimizable
	do n2=n1,l-1
	do n3=n2,l-1

	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3
	mqloop(1) = m_period(n1,l)
	mqloop(2) = m_period(n2,l)
	mqloop(3) = m_period(n3,l)

	   multi = 1.0d0*multiplicity(n1,n2,n3)

  	   D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	S_prop_q(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	S_prop_mq(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

  	   call Gluon_propagator(l,qloop,D_prop)	
	   call Fermion_Propagator(l,nvec_1+qloop,S_prop_q)
	   call Fermion_Propagator(l,nvec_1+mqloop,S_prop_mq)

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
	  do nu=0,3

	  call V1(l,nvec_1-qloop,-nvec_1,qloop,s0,t0,u0,mu,V1mu)
	  call V1(l,nvec_1+qloop,-nvec_1,-qloop,s0p,t0p,u0p,nu,V1nu)

	  trace = Trace6(SOSG1(:,:,t0p,x0),S_prop_mq(:,:,x0,s0,flav(2)),V1mu(:,:),SOSG2(:,:,t0,x0),S_prop_q(:,:,x0,s0p,flav(4)),V1nu(:,:))
	  g_loop(x0) = g_loop(x0) + multi*trace*D_prop(nu,mu,u0p,u0)

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

	enddo	!x0

	  G24f_13a(:) = C2_R*dimrep*g_loop(:)/(1.0d0*l**3)

	END SUBROUTINE G2_4f_13a



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE G2_4f_13b(l,g_i1,g_i2,flav,G24f_13b)

	implicit none

	integer,intent(in) :: l, g_i1, g_i2
	integer,dimension(4),intent(in) :: flav
	complex(kind=8),dimension(c0:c1),intent(out) :: G24f_13b

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: D_prop
	complex(kind=8),dimension(4,4,0:l,0:l,2) :: S_prop_q, S_prop_mq, S_prop_0
	complex(kind=8),dimension(4,4,0:l,c0:c1) :: SOSG1, SOSG2
	complex(kind=8),dimension(4,4) :: V1mu, V1nu
	complex(kind=8),dimension(4,4) :: mat1, mat2, mat3	
	integer,dimension(3) :: qloop,mqloop
	complex(kind=8),dimension(c0:c1) :: g_loop, tree
	integer :: x0,i
	integer :: s0,s0p,t0,t0p,u0,u0p
	integer :: t0_a,t0_b,t0p_a,t0p_b,u0_a,u0_b,u0p_a,u0p_b
	integer :: mu, nu
	integer :: n1, n2, n3
	complex(kind=8) :: trace
	integer,dimension(3) :: nvec_1, nvec_2
	real(kind=8) :: multi	
	integer :: f,fp

	nvec_1(:) = 0
	nvec_2(:) = 0

	G24f_13b(:) = dcmplx(0.0d0,0.0d0)
	g_loop(:) = dcmplx(0.0d0,0.0d0)
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_q(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_mq(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	S_prop_0(:,:,:,:,:)=dcmplx(0.0d0,0.0d0)
	SOSG1(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	SOSG2(:,:,:,:) = dcmplx(0.0d0,0.0d0)
	trace = dcmplx(0.0d0,0.0d0)
	f=f1f2(flav(4),flav(1))
	fp=f1f2(flav(2),flav(3))

	!part outside the loop

	!build S_O5_S
	call Fermion_Propagator(l,nvec_1,S_prop_0)

	do x0=c0,c1
	do s0=0,l

	    mat1(:,:) = matmul(gBil(:,:,g_i2),S_prop_0(:,:,x0,1,flav(4)))
	    mat2(:,:) = matmul(Q5(:,:,f),S_prop_0(:,:,1,s0,flav(1)))
	    SOSG1(:,:,x0,s0) = matmul(mat1(:,:),mat2(:,:))

	    mat1(:,:) = matmul(gBil(:,:,g_i1),S_prop_0(:,:,x0,l-1,flav(2)))
	    mat2(:,:) = matmul(Q5p(:,:,fp),S_prop_0(:,:,l-1,s0,flav(3)))
	    SOSG1(:,:,x0,s0) = matmul(mat1(:,:),mat2(:,:))
	enddo
	enddo	

	!Construct the gluon loop
	do x0 = c0,c1

	do n1=0,l-1			!loop momenta. Optimizable
	do n2=n1,l-1
	do n3=n2,l-1

	qloop(1)=n1
	qloop(2)=n2
	qloop(3)=n3

	mqloop(1) = m_period(n1,l)
	mqloop(2) = m_period(n2,l)
	mqloop(3) = m_period(n3,l)

	   multi = 1.0d0*multiplicity(n1,n2,n3)

  	   D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	S_prop_q(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)
  	S_prop_mq(:,:,:,:,:) = dcmplx(0.0d0,0.0d0)

  	   call Gluon_propagator(l,qloop,D_prop)	
	   call Fermion_Propagator(l,nvec_1+qloop,S_prop_q)
	   call Fermion_Propagator(l,nvec_1+mqloop,S_prop_mq)


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
	  do nu=0,3

	  call V1(l,nvec_1,-nvec_1-qloop,qloop,s0,t0,u0,mu,V1mu)
	  call V1(l,nvec_1,-nvec_1+qloop,-qloop,s0p,t0p,u0p,nu,V1nu)

	  trace = Trace6(SOSG1(:,:,x0,s0),V1mu(:,:),S_prop_q(:,:,t0,x0,flav(1)),SOSG2(:,:,x0,s0p),V1nu(:,:),S_prop_mq(:,:,t0p,x0,flav(3)))
	  g_loop(x0) = g_loop(x0) + multi*trace*D_prop(nu,mu,u0p,u0)

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

	enddo	!x0

	  G24f_13b(:) = C2_R*dimrep*g_loop(:)/(1.0d0*l**3)

	END SUBROUTINE G2_4f_13b



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	END MODULE OneLoop_4fermions
