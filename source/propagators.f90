

	MODULE Propagators
	!
	! Module containing functions to compute the fermion and gauge propagators.
	! Contains:
	!		- Fermion propagator for the XSF
	!		- Gluon propagator Wilson action
	!		- Gluon propagator Luscher Weisz action
	!
	!
	USE GammaMatrices
	USE Input
	USE Parameters

	implicit none
	CONTAINS



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!		FERMION PROPAGATOR ROUTINES FOR THE SF
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE Fermion_Propagator(l,nvec,S_Prop)
!         Computes the fermion propagator for the u-type and d-type
!	 quarks in the XSF in absence of 
!         background field. The Dirac operator of u-type quarks is built
!	 and inverted numerically. The d-type propagator is obtained using 
!	 the property S^(d) = g5 (S^(u))^dag g5.
!	
!	 The numerical inversion is performed takin advantage of the fact that
!	 the Dirac operator has a band struture. Only the 12-diagonal band is non
!	 zero.
!
!        For BF=0 the Dirac operator is trivial in color space and all color
!	 comonents of the propagator are equal.
!
!        -Input: -l: lattice size 
!		-nvec: fermion momentum
!	
!	-Output: -S_prop_u: u-type fermion propagator (4,4,0:l,0:l) matrix 
!		  (The propagator isdefined at the boundaries)
!		 -S_prop_d: d-type fermion propagator (4,4,0:l,0:l) matrix 
!		  (The propagator isdefined at the boundaries)


	implicit none
	
	!! Input and output variables
	
	integer,intent(in) :: l			!Physical size
	integer,dimension(3),intent(in) :: nvec	!Entering momentum
	complex(Kind=8),dimension(4,4,0:l,0:l,2),intent(out) :: S_Prop !2 flavours


	complex(Kind=8),dimension(4,4,0:l,0:l) :: S_Prop_d_temp
	!!! Internal variables
	!Counters
	integer :: i,j,t,s,dir,i1,j1

	!Dirac scalars in the Dirac operator
	real(kind=8),dimension(3) :: p, phat, ptwit !dimension(dir=3)
	complex(kind=8) :: piece1

	!Matrices involved in the construction of the Dirac operator
	complex(kind=8),dimension(4,4) :: Block1, Block2, hmat, Inv
	complex(kind=8),dimension(4,4) :: Aterm, Bterm, Cterm, Dterm, Eterm

	!The Dirac operator and Propagator
	complex(kind=8),dimension(4*(l+1),4*(l+1)) :: Dop
	complex(kind=8),dimension(4*(l+1),4*(l+1)) :: S_col_p,S_Prop_d_big


	!Lapack variables
	integer :: ok1,ok2,wlwork
	integer,parameter :: KL=5, KU=5
	complex*16,dimension(4*(l+1)) ::wwork
	integer,dimension(4*(l+1)) :: pivot
	complex(kind=8),dimension((2*KL+KU+1),4*(l+1)) :: S_array


	 Dop(:,:) = dcmplx(0.0d0,0.0d0)
	 S_col_p(:,:) = dcmplx(0.0d0,0.0d0)	

	 p(1)=(2.0d0*pi*nvec(1)+theta)/(l*1.0d0)	!momenta, plus theta			
	 p(2)=(2.0d0*pi*nvec(2)+theta)/(l*1.0d0)
	 p(3)=(2.0d0*pi*nvec(3)+theta)/(l*1.0d0)

	  do dir=1,3						!3 directions
	   phat(dir)=2.0d0*sin(0.5d0*p(dir))
	   ptwit(dir)=sin(p(dir))						
          enddo

	piece1=dcmplx(0.0d0,0.0d0)
	Block1(:,:)=dcmplx(0.0d0,0.0d0)
	Block2(:,:)=dcmplx(0.0d0,0.0d0)
	hmat(:,:)=dcmplx(0.0d0,0.0d0)

	!!Define the parts of the Dirac operator:
	do dir=1,3
	  piece1=piece1+phat(dir)**2
	enddo
	piece1=piece1*0.5d0
	Block1=piece1*II
	Block2=imag*(ptwit(1)*g1+ptwit(2)*g2+ptwit(3)*g3)

	Aterm = (1+m0)*II + Block1 +Block2
	Bterm = -Pminus
	Cterm = -Pplus
	Eterm = imag*g5_Pminus + (m0 + zf)*II + ds*(Block1 + Block2)
	Dterm = imag*g5_Pplus + (m0 + zf)*II + ds*(Block1 + Block2)

	!!!!!! Construction of the operator Dop

	!!!!!Diagonal part

	!Aterm
	do t=1,l-1
	do i1=1,4
	do j1=1,4
	
	Dop(4*t+i1,4*t+j1) = Aterm(i1,j1)

	enddo
	enddo
	enddo

	!Eterm
	do i1=1,4
	do j1=1,4

	Dop(i1,j1) = Eterm(i1,j1)

	enddo
	enddo

	!Dterm
	do i1=1,4
	do j1=1,4

	Dop(4*l+i1,4*l+j1) = Dterm(i1,j1)

	enddo
	enddo

	!Bterm
	do t=0,l-1
	do i1=1,4
	do j1=1,4
	
	Dop(4*t+i1,4*(t+1)+j1) = Bterm(i1,j1)
	
	enddo
	enddo
	enddo

	!Cterm
	do t=0,l-1
	do i1=1,4
	do j1=1,4
	
	Dop(4*(t+1)+i1,4*t+j1) = Cterm(i1,j1)
	
	enddo
	enddo
	enddo
	
	
	!Dirac Operator built at this point.
	!We perform the inversion

	!Bring the matrix to a form for the lapack routine

	do j=1,4*(l+1)
	do i= max(1,j-KU),min(4*(l+1),j+KL)
	   S_array(KU+KL+1+i-j,j) = Dop(i,j)
	enddo
	enddo

	do t=1,4*(l+1)
	S_col_p(t,t) = dcmplx(1.0d0,0.0d0)
	enddo

	! For a referenceof this functions see: www.physics.orst.edu/~rubin/nacphy/lapack/linear.html
	call ZGBTRF(4*(l+1),4*(l+1),KL,KU,S_array,(2*KL+KU+1),pivot,ok1)		
	call ZGBTRS('N',4*(l+1),KL,KU,4*(l+1),S_array,(2*KL+KU+1),pivot,S_col_p,4*(l+1),ok2)

	! And we store the propagator in the S_Prop matrix. Since the tree level propagator
	! is color diagonal, we only write it as a color vector.
	do s=1,l+1
	do t=1,l+1
	 do i=1,4
	 do j=1,4
	   S_Prop(i,j,t-1,s-1,1)=S_col_p((t-1)*4+i,(s-1)*4+j)
	 enddo
	 enddo
	enddo
	enddo


	!Calculation of the hermithean conjugate

	do i=1,4*(l+1)
	do j=1,4*(l+1)

	S_Prop_d_big(j,i) = dcmplx(real(S_col_p(i,j)),-aimag(S_col_p(i,j)))

	enddo
	enddo

	!Building the propagator for the down type quark

	do s=1,l+1
	do t=1,l+1
	 do i=1,4
	 do j=1,4
	   S_Prop_d_temp(i,j,t-1,s-1)=S_Prop_d_big((t-1)*4+i,(s-1)*4+j)
	 enddo
	 enddo
	enddo
	enddo

	do s=0,l
	do t=0,l
	  
	  S_Prop(:,:,s,t,2) = g5Ag5_shuf(S_Prop_d_temp(:,:,s,t))

	enddo
	enddo


!	enddo !End of the color loop

	END SUBROUTINE Fermion_Propagator


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




	SUBROUTINE Gluon_Propagator_pl(l,nvec,D_prop) !!!This routine seems fine for p=0.
	!
	! Subroutine to compute the free gluon propagator
	! In order to build the different coefficients we will use the
	! functions sn, cn and lambda_n.
	!
	integer,intent(in) :: l
	integer,dimension(3),intent(in) :: nvec
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1),intent(out) :: D_prop !Now this is the output
	
	real(kind=8),dimension(3) :: p,phat,pdot
	real(kind=8) :: phat_sq
	integer :: i,j,k,k2
	real(kind=8),dimension(0:l-1,0:l-1) :: d,b,c,n,e
	integer :: zero_mom

	!initialization
	d(:,:)=0.0d0
	b(:,:)=0.0d0
	c(:,:)=0.0d0
	e(:,:)=0.0d0
	n(:,:)=0.0d0
	D_prop(:,:,:,:)=dcmplx(0.0d0,0.0d0)

	zero_mom = nvec(1) + nvec(2) + nvec(3)		!Variable to check if we are at zero momentum or not

	!construction of momenta
	do i=1,3
 	  p(i)=(2.0d0*pi/l)*nvec(i)
	  phat(i)=2.0d0*sin(p(i)/2.0d0)
	enddo

	phat_sq=phat(1)*phat(1)+phat(2)*phat(2)+phat(3)*phat(3)
	
	!construction of the coefficients b,c,d,n,e
	do i=0,l-1		
	do j=0,l-1		


 	 if(zero_mom.eq.0)then !If we are at zero momentum use one formula
	 
	   !d
	   if(i.gt.j)then
	     d(i,j) = j*(1.0d0 - i*1.0d0/l*1.0d0)
	   else 	 	
	     d(i,j) = i*(1.0d0 - j*1.0d0/l*1.0d0)
 	   endif

	   !n,e
	   n(i,j)=1.0d0+min(i,j)		!Test what happens if i=j
	   e(i,j)=1.0d0+min(i,j)

	 else !If we are not at zero momentum, use another

	   do k=1,l-1
	     d(i,j) = d(i,j) + sn(k,i,l)*sn(k,j,l)/(phat_sq + lambda_n(k,l)**2)
	     b(i,j) = b(i,j) + sn(k,i,l)*sn(k,j,l)/(phat_sq + lambda_n(k,l)**2)**2
	     c(i,j) = c(i,j) + lambda_n(k,l)*sn(k,i,l)*cn(k,j,l)/(phat_sq + lambda_n(k,l)**2)**2
	     e(i,j) = e(i,j) + lambda_n(k,l)**2*cn(k,i,l)*cn(k,j,l)/(phat_sq + lambda_n(k,l)**2)**2
	   enddo

	   do k=0,l-1
   	     n(i,j) = n(i,j) + cn(k,i,l)*cn(k,j,l)/(phat_sq + lambda_n(k,l)**2)
	   enddo

	  d(i,j) = d(i,j)*2.0d0/(l*1.0d0)
	  b(i,j) = b(i,j)*2.0d0/(l*1.0d0)
	  c(i,j) = c(i,j)*2.0d0/(l*1.0d0)
	  e(i,j) = e(i,j)*2.0d0/(l*1.0d0)
	  n(i,j) = n(i,j)*2.0d0/(l*1.0d0)

	 endif

	enddo
	enddo


	! Construction of the propagator
	do i=0,l-1
	do j=0,l-1
	! D00
  	  D_prop(0,0,i,j) = n(i,j) + (1.0d0/lambda_gf - 1.0d0)*e(i,j)
	enddo
	enddo	

	

	!Dk0
	do i=1,l-1
	do j=0,l-1
	do k=1,3
	    D_prop(k,0,i,j) = imag*phat(k)*(1.0d0/lambda_gf - 1.0d0)*c(i,j)
	enddo

	enddo
	enddo


	!D0k
	do i=0,l-1
	do j=1,l-1

	do k=1,3
	    D_prop(0,k,i,j) = -D_prop(k,0,j,i)
	enddo

	enddo
	enddo


	!Dkl
	do i=1,l-1
	do j=1,l-1

	do k=1,3
	do k2=1,3
	    D_prop(k,k2,i,j) = delta(k,k2)*d(i,j) + phat(k)*phat(k2)*(1.0d0/lambda_gf - 1.0d0)*b(i,j) 
	enddo
	enddo

	enddo
	enddo	    

	END SUBROUTINE Gluon_Propagator_pl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE Gluon_Propagator(l,nvec,D_prop) 
	!
	! Subroutine to compute the free gluon propagator
	! In order to build the different coefficients we will use the
	! functions sn, cn and lambda_n.
	!
	integer,intent(in) :: l
	integer,dimension(3),intent(in) :: nvec
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1),intent(out) :: D_prop !Now this is the output
	
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: Kglue 
	complex(kind=8),dimension(4*l-3,4*l-3) :: K_big_prop
	!Lapack variables
	integer :: ok1,ok2,wlwork
	complex*16,dimension(4*l-3) ::wwork
	integer,dimension(4*l-3) :: pivot


	Kglue(:,:,:,:)=dcmplx(0.0d0,0.0d0)
	K_big_prop(:,:) = dcmplx(0.0d0,0.0d0)	

	if(G_act.eq.1)then
	  call Gluonic_Quadratic_Operator_Pl(l,nvec,Kglue)
	else
	  call Gluonic_Quadratic_Operator_LW(l,nvec,Kglue)
	endif
	call D_to_Dmat(l,Kglue,K_big_prop)
	call ZGETRF(4*l-3,4*l-3,K_big_prop,4*l-3,pivot,ok1)		
	call ZGETRI(4*l-3,K_big_prop,4*l-3,pivot,wwork,4*l-3,ok2)

	!Get the propagator
	call Dmat_to_D(l,K_big_prop,D_prop)


	END SUBROUTINE Gluon_Propagator


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



	SUBROUTINE Gluonic_Quadratic_Operator_Pl(l,nvec,Kglue)
	!
	!	Subroutine to compute analytically the gluonic quadratic operator
	!	as found it eq.(2.7-2.11) of Peter Weisz's notes:
	!	"Computation of the improvement coefficient c_A to 1-loop",
	!	December 1995.
	!
	!	-Input: -l: lattice size
	!		-nvec: Momenta
	!
	!	-Output: Kglue: gluon quadratic operator  dimension(0:3,0:3,0:l-1,0:l-1)	
	!
	implicit none


	integer,intent(in) :: l
	integer,dimension(3),intent(in) :: nvec
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1),intent(out) :: Kglue 

	real(kind=8),dimension(3) :: p,phat,pdot
	real(kind=8) :: phat_sq
	integer :: i,j,k,k2,s,t
	real(kind=8) :: piece1,piece2,piece3

	!Initialization

	Kglue(:,:,:,:)=dcmplx(0.0d0,0.0d0)

	do i=1,3
	p(i)=(2.0d0*pi/l)*nvec(i)
	phat(i)=2.0d0*sin(p(i)/2.0d0)
	pdot(i)=sin(p(i))
	enddo

	phat_sq=phat(1)*phat(1)+phat(2)*phat(2)+phat(3)*phat(3)	


	! HERE WE CALCULATE THE QUADRATIC GLUON OPERATOR WITH 0 BACKGROUND FIELD

	!Kglue 00

	do i=0,l-1
	do j=0,l-1
	  piece1 = lambda_gf*(2.0d0*delta(i,j) - delta(i+1,j) - delta(i-1,j))
	  piece2 = delta(i,0)*(1 - delta(nvec(1),0)*delta(nvec(2),0)*delta(nvec(3),0)) + delta(i,l-1)
	  piece3 = lambda_gf*delta(i,j)*piece2
	  Kglue(0,0,i,j) = delta(i,j)*phat_sq + piece1 - piece3
	enddo
	enddo

	! Kglue(k,0)

	do k=1,3
	do i=1,l-1
	do j=0,l-1
   	  Kglue(k,0,i,j) = imag*phat(k)*(1.0d0-lambda_gf)*(delta(i,j)-delta(i-1,j))
	enddo
	enddo
	enddo

	! Kglue(0,k)

	do k=1,3
	do i=0,l-1
	do j=1,l-1
	  Kglue(0,k,i,j) = -Kglue(k,0,j,i)
	enddo
	enddo
	enddo

	
	! Kglue(k,k2)
	
	do k=1,3
	do k2=1,3
	  piece1 = delta(k,k2)*phat_sq - phat(k)*phat(k2)*(1.0d0 - lambda_gf)
	  do i=1,l-1
	  do j=1,l-1
	    piece2 = 2*delta(i,j) - delta(i+1,j) - delta(i-1,j)
	    Kglue(k,k2,i,j) = delta(i,j)*piece1+delta(k,k2)*piece2
	  enddo
	  enddo
	enddo
	enddo


	END SUBROUTINE Gluonic_Quadratic_Operator_Pl





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE Gluonic_Quadratic_Operator_LW(l,nvec,Kglue)
	!
	!	Subroutine to compute the gluon propagator for the improved theory
	!	described in arXiv:hep-lat/9808007.
	!
	implicit none


	integer,intent(in) :: l
	integer,dimension(3),intent(in) :: nvec
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1),intent(out) :: Kglue    

	real(kind=8),dimension(3) :: p,phat,pdot
	real(kind=8) :: phat_sq, phat_fth
	integer :: i,j,k,k2,s,t
	real(kind=8) :: piece1,piece2,piece3,piece4, piece5
	real(kind=8),parameter :: c0_S = 5.0d0/3.0d0
	real(kind=8),parameter :: c1_S = -1.0d0/12.0d0
	real(kind=8),parameter :: c3_S = 0.0d0
	real(kind=8),parameter :: d_bc_S = 1.0d0

	!Initialization

	Kglue(:,:,:,:)=dcmplx(0.0d0,0.0d0)


	do i=1,3
	p(i)=(2.0d0*pi/l)*nvec(i)
	phat(i)=2.0d0*sin(p(i)/2.0d0)
	pdot(i)=sin(p(i))
	enddo

	phat_sq = phat(1)*phat(1)+phat(2)*phat(2)+phat(3)*phat(3)	
	phat_fth = phat(1)**4 + phat(2)**4 + phat(3)**4


	! HERE WE CALCULATE THE QUADRATIC GLUON OPERATOR WITH 0 BACKGROUND FIELD

	!Kglue 00

	do i=0,l-1
	do j=0,l-1

	  piece1 = (c0_S + 6.0d0*c1_S + 8.0d0*c3_S)*phat_sq -c1_S*phat_fth +c3_S*(phat_fth - phat_sq**2)
	  piece2 = lambda_gf*(2.0d0 -delta(i,0)*(1.0d0-delta(nvec(1),0)*delta(nvec(2),0)*delta(nvec(3),0)) - delta(i,l-1))
	  piece3 = (delta(i-1,j)+delta(i+1,j))*(c1_S*phat_sq - lambda_gf)
	  piece4 = c1_S*delta(i,j)*(delta(i,l-1) + delta(i,0))*(phat_sq - phat_fth*d_bc_S)

	  Kglue(0,0,i,j) = delta(i,j)*(piece1 + piece2) + piece3 + piece4
	enddo
	enddo

	! Kglue(k,0)

	do k=1,3
	do i=1,l-1
	do j=0,l-1

	  piece1 = (delta(i,j) - delta(i-1,j))*(c0_S + 5.0d0*c1_S + 8.0d0*c3_S -lambda_gf - c1_S*phat(k)**2 +c3_S*(phat(k)**2 - phat_sq))
	  piece2 = c1_S*(delta(i+1,j)-delta(i-2,j)) + c1_S*(delta(i,j)*delta(i,l-1) - delta(i-1,j)*delta(i,1))*(1.0d0 - d_bc_S*phat(k)**2)

   	  Kglue(k,0,i,j) = imag*phat(k)*(piece1 + piece2)
	enddo
	enddo
	enddo


	! Kglue(0,k)

	do k=1,3
	do i=0,l-1
	do j=1,l-1
	  Kglue(0,k,i,j) = -Kglue(k,0,j,i)
	enddo
	enddo
	enddo

	


	do k=1,3
	do k2=1,3

	piece1 =(c0_S + 8.0d0*c1_S + 4.0d0*c3_S)*phat_sq + 2.0d0*c0_S + 10.0*c1_S + 16.0*c3_S
	piece1 = piece1 - c1_S*(2.0d0*phat(k)**2 + phat_fth + phat_sq*(phat(k)**2))
	piece1 = piece1 + c3_S*(2.0d0*phat(k)**2 + phat_fth -phat_sq**2 + phat_sq*(phat(k)**2))
	piece1 = delta(k,k2)*piece1

	piece2 = lambda_gf - c0_S - 8.0d0*c1_S - 6.0d0*c3_S + c1_S*(phat(k)**2+phat(k2)**2)
	piece2 = phat(k)*phat(k2)*(piece2 + c3_S*(phat_sq-phat(k)**2-phat(k2)**2))

	piece3 = delta(k,k2)*(c0_S + 4.0d0*c1_S + 8.0d0*c3_S - c1_S*phat(k)**2 + c3_S*(phat(k)**2-2.0d0*phat_sq))
	piece3 = piece3 + c3_S*phat(k)*phat(k2)

	piece4 = c1_S*delta(k,k2)

	piece5 = c1_S*delta(k,k2)*(1.0d0-phat(k)**2*d_bc_S)
	! Kglue(k,k2)
	
	do i=1,l-1
	do j=1,l-1

	  Kglue(k,k2,i,j) = delta(i,j)*(piece1+piece2)
	  Kglue(k,k2,i,j) = Kglue(k,k2,i,j) - (delta(i+1,j) + delta(i-1,j))*piece3
	  Kglue(k,k2,i,j) = Kglue(k,k2,i,j) - (delta(i+2,j) + delta(i-2,j))*piece4
	  Kglue(k,k2,i,j) = Kglue(k,k2,i,j) + delta(i,j)*(delta(i,l-1) + delta(i,1))*piece5

	enddo
	enddo

	enddo
	enddo

	END SUBROUTINE Gluonic_Quadratic_Operator_LW

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE Gluon_propagator_LW(l,nvec,D_prop)
	!
	!	Subroutine to compute the gluon propagator for the improved theory
	!	described in arXiv:hep-lat/9808007.
	!
	!	First, the quadatic gluon operator is computed using the analytic 
	!	expressions in arXiv:hep-lat/9808007. 
	!
	!	This is inverted numericall to obtain the gluon propagator.
	!
	!	The parameters c0, c1, c2 and d_bc of the LW action are defined here
	!
	!	-Input: -l: lattice size
	!		-nvec: Momenta
	!
	!	-Output: D_prop: gluon propagator  dimension(0:3,0:3,0:l-1,0:l-1)
	!
	implicit none


	integer,intent(in) :: l
	integer,dimension(3),intent(in) :: nvec
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1),intent(out) :: D_prop    

	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1) :: Kglue !Now this is the output
	complex(kind=8),dimension(4*l-3,4*l-3) :: K_big_prop
	complex(kind=8),dimension(4*l-3,4*l-3) :: K_big_op
		

	real(kind=8),dimension(3) :: p,phat,pdot
	real(kind=8) :: phat_sq, phat_fth
	integer :: i,j,k,k2,s,t
	real(kind=8) :: piece1,piece2,piece3,piece4, piece5
	real(kind=8),parameter :: c0_S = 5.0d0/3.0d0
	real(kind=8),parameter :: c1_S = -1.0d0/12.0d0
	real(kind=8),parameter :: c3_S = 0.0d0
	real(kind=8),parameter :: d_bc_S = 1.0d0


	!Lapack variables
	integer :: ok1,ok2,wlwork
	complex*16,dimension(4*l-3) ::wwork
	integer,dimension(4*l-3) :: pivot


	!Initialization

	Kglue(:,:,:,:)=dcmplx(0.0d0,0.0d0)


	do i=1,3
	p(i)=(2.0d0*pi/l)*nvec(i)
	phat(i)=2.0d0*sin(p(i)/2.0d0)
	pdot(i)=sin(p(i))
	enddo

	phat_sq = phat(1)*phat(1)+phat(2)*phat(2)+phat(3)*phat(3)	
	phat_fth = phat(1)**4 + phat(2)**4 + phat(3)**4


	! HERE WE CALCULATE THE QUADRATIC GLUON OPERATOR WITH 0 BACKGROUND FIELD

	!Kglue 00

	do i=0,l-1
	do j=0,l-1

	  piece1 = (c0_S + 6.0d0*c1_S + 8.0d0*c3_S)*phat_sq -c1_S*phat_fth +c3_S*(phat_fth - phat_sq**2)
	  piece2 = lambda_gf*(2.0d0 -delta(i,0)*(1.0d0-delta(nvec(1),0)*delta(nvec(2),0)*delta(nvec(3),0)) - delta(i,l-1))
	  piece3 = (delta(i-1,j)+delta(i+1,j))*(c1_S*phat_sq - lambda_gf)
	  piece4 = c1_S*delta(i,j)*(delta(i,l-1) + delta(i,0))*(phat_sq - phat_fth*d_bc_S)

	  Kglue(0,0,i,j) = delta(i,j)*(piece1 + piece2) + piece3 + piece4
	enddo
	enddo

	! Kglue(k,0)

	do k=1,3
	do i=1,l-1
	do j=0,l-1

	  piece1 = (delta(i,j) - delta(i-1,j))*(c0_S + 5.0d0*c1_S + 8.0d0*c3_S -lambda_gf - c1_S*phat(k)**2 +c3_S*(phat(k)**2 - phat_sq))
	  piece2 = c1_S*(delta(i+1,j)-delta(i-2,j)) + c1_S*(delta(i,j)*delta(i,l-1) - delta(i-1,j)*delta(i,1))*(1.0d0 - d_bc_S*phat(k)**2)

   	  Kglue(k,0,i,j) = imag*phat(k)*(piece1 + piece2)
	enddo
	enddo
	enddo

	! Kglue(0,k)

	do k=1,3
	do i=0,l-1
	do j=1,l-1
	  Kglue(0,k,i,j) = -Kglue(k,0,j,i)
	enddo
	enddo
	enddo

	
	! Kglue(k,k2)
	
	do k=1,3
	do k2=1,3

	  piece1 = (c0_S + 8.0d0*c1_S + 4.0d0*c3_S)*phat_sq + 2.0d0*c0_S + 10.0d0*c1_S + 16.0d0*c3_S
	  piece1 = piece1 - c1_S*(2.0d0*phat(k)**2 + phat_fth + phat_sq*phat(k)**2)
	  piece1 = piece1 + c3_S*(2.0d0*phat(k)**2 + phat_fth - phat_sq**2 +phat_sq*phat(k)**2)
	  piece1 = delta(k,k2)*piece1

	  piece2 = lambda_gf - c0_S -8.0d0*c1_S -6.0d0*c3_S +c1_S*(phat(k)**2+phat(k2)**2) + c3_S*(phat_sq - phat(k)**2 - phat(k2)**2)
	  piece2 = phat(k)*phat(k2)*piece2

	  piece3 = delta(k,k2)*(c0_S + 4.0d0*c1_S + 8.0d0*c3_S -c1_S*phat(k)**2 + c3_S*(phat(k)**2 - 2.0d0*phat_sq))
	  piece3 = piece3 + c3_S*phat(k)*phat(k2)

	  piece4 = delta(k,k2)*c1_S

	  piece5 = c1_S*delta(k,k2)*(1.0d0 - phat(k)**2*d_bc_S) 


	  do i=1,l-1
	  do j=1,l-1

    	     Kglue(k,k2,i,j) = delta(i,j)*(piece1 + piece2) - (delta(i+1,j) + delta(i-1,j))*piece3
    	     Kglue(k,k2,i,j) = Kglue(k,k2,i,j) - (delta(i+2,j) + delta(i-2,j))*piece4	     
    	     Kglue(k,k2,i,j) = Kglue(k,k2,i,j) + (delta(i,l-1) + delta(i,1))*piece5
	  enddo
	  enddo
	enddo
	enddo


	!Build the operator for inverting:
	K_big_prop(:,:) = dcmplx(0.0d0,0.0d0)	

	call D_to_Dmat(l,Kglue,K_big_op)

	K_big_prop = K_big_op
	!Invert the gluonic operator	

	call ZGETRF(4*l-3,4*l-3,K_big_prop,4*l-3,pivot,ok1)		
	call ZGETRI(4*l-3,K_big_prop,4*l-3,pivot,wwork,4*l-3,ok2)

	!Get the propagator
	call Dmat_to_D(l,K_big_prop,D_prop)



	END SUBROUTINE Gluon_propagator_LW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




	SUBROUTINE print_gluon_propagator(l,D_prop)
	!
	! Prints on screen all the components of the gluon propagator. 
	!
	integer,intent(in) :: l
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1),intent(in) :: D_prop

	integer :: mu, nu, t0, s0

	do t0 = 0,l-1
	do s0 = 0,l-1

	print*, "t0, s0 =",t0,s0

	do mu = 0,3
        do nu = 0,3

	  print*, mu, nu, real(D_prop(mu,nu,t0,s0)), aimag(D_prop(mu,nu,t0,s0))
	
	enddo
	enddo

	enddo
	enddo
	


	END SUBROUTINE print_gluon_propagator




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE D_to_Dmat(l,D_prop,D_mat)
	!
	! Subroutine that converts a matrix (0:3,0:3,0:l-1,0:l-1)  
	! into a big matrix of dimension (0:3)x(0:l-1),(0:3)x(0:l-1).
	! This is used tipically in the construction of the gluon propagator.
	!
	implicit none

	integer,intent(in) :: l
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1),intent(in) :: D_prop
	complex(kind=8),dimension(4*l-3,4*l-3),intent(out) :: D_mat

	integer :: mu, nu, t0, s0


	!Build the operator for inverting:
	D_mat(:,:) = dcmplx(0.0d0,0.0d0)	

	do t0=1,l-1
	do s0=1,l-1
	  do mu=0,3
	  do nu=0,3
	   D_mat(4*t0+(mu+1)-3,4*s0+(nu+1)-3) = D_prop(mu,nu,t0,s0)
	  enddo
	  enddo
	enddo
	enddo


	do t0=1,l-1
	  do mu=0,3
	   D_mat(4*t0+(mu+1)-3,1) = D_prop(mu,0,t0,0)
	  enddo
	enddo

	do s0=1,l-1
	  do nu=0,3
	   D_mat(1,4*s0+(nu+1)-3) = D_prop(0,nu,0,s0)
	  enddo
	enddo

	D_mat(1,1) = D_prop(0,0,0,0)



	END SUBROUTINE D_to_Dmat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE Dmat_to_D(l,D_mat,D_prop)
	!
	! Subroutine that converts a matrix (0:3)x(0:l-1),(0:3)x(0:l-1) 
	! into a matrix of dimension (0:3,0:3,0:l-1,0:l-1).
	! This is used tipically in the construction of the gluon propagator.
	!
	implicit none

	integer,intent(in) :: l
	complex(kind=8),dimension(4*l-3,4*l-3),intent(in) :: D_mat
	complex(kind=8),dimension(0:3,0:3,0:l-1,0:l-1),intent(out) :: D_prop

	integer :: mu, nu, t0, s0


	!Build the operator for inverting:
	D_prop(:,:,:,:) = dcmplx(0.0d0,0.0d0)	





	do t0=1,l-1
	do s0=1,l-1
	  do mu=0,3
	  do nu=0,3
	    D_prop(mu,nu,t0,s0) = D_mat(4*t0+(mu+1)-3,4*s0+(nu+1)-3)
	  enddo
	  enddo
	enddo
	enddo


	do t0=1,l-1
	  do mu=0,3
	   D_prop(mu,0,t0,0) = D_mat(4*t0+(mu+1)-3,1)
	  enddo
	enddo

	do s0=1,l-1
	  do nu=0,3
	    D_prop(0,nu,0,s0) = D_mat(1,4*s0+(nu+1)-3)
	  enddo
	enddo

 	 D_prop(0,0,0,0) = D_mat(1,1)




	END SUBROUTINE Dmat_to_D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE check_inverse_matrix(l,mat1,mat2)
	!
	! Subroutine to check that A*A^-1 = I
	!
	implicit none

	integer,intent(in) :: l
	complex(kind=8),dimension(4*l-3,4*l-3),intent(in) :: mat1,mat2


	integer :: i,j
	complex(kind=8),dimension(4*l-3,4*l-3) :: mat3

	mat3 = matmul(mat1,mat2)


	do i=1,4*l-3
	do j=1,4*l-3

	  print*, i,j, real(mat3(i,j))

	enddo
	enddo


	END SUBROUTINE check_inverse_matrix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!Functions for the gluon propagator

	REAL(kind=8) FUNCTION sn(n,a,T)
	!
	! Computes sn in the gluon propagator
	!
	integer :: n,a,T	
	
	sn = sin(n*pi*a/(T*1.0d0))		!note that x0=a-1

	END FUNCTION sn




	REAL(kind=8) FUNCTION cn(n,a,T)
	!
	! Computes cn in the gluon propagator
	!
	integer :: n,a,T

	if(n.eq.0)then
	  cn = 1.0d0 / sqrt(2.0d0)
	else
	  cn = cos(n*pi*(a*1.0d0 + 0.5d0)/(T*1.0d0))  !note that x0=a-1
	endif

	END FUNCTION cn




	REAL(kind=8) FUNCTION lambda_n(n,T)
	!
	! Computes lambda_n in the gluon propagator
	!
	integer :: n,T

	lambda_n=2.0d0*sin(n*pi/(2.0d0*T))

	END FUNCTION lambda_n


	END MODULE Propagators
	
