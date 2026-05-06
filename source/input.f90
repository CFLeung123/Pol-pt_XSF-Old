


	MODULE Input
	!
	!  Module containing
	!
	!	- Declaration of all input parameters
	!	- Subroutine to read input parameters from file
	!	  If choosen so, the tree-level values of the critical mass
	!	  and zf are read from file.
	!	- Subroutine to print in screen the relevant input parameters.
	!	- Functions to compute the upper and lower limits of a space-time
	!	  inde mu, deppending on wether the chosen gauge is diagonal or not.
	implicit none

	!input parameters
    	character*10 :: text
	integer :: lat_size
	real(kind=8) :: m0, csw, ds, zf, theta, lambda_gf
	integer :: deriv, m0_zf_file, G_act, mom_deg, diag_g
	integer :: gam_i,gam_i2, flav1, flav2, flav3 , flav4, gam4f_i,gam4f_i2, gam4f_pair
	integer :: comp_gX, comp_lY, comp_g1, comp_l1, comp_gVt, comp_lVt
	integer :: comp_G1_4f, comp_G2_4f, comp_Gk1_4f, comp_Gk2_4f
	integer :: comp_tree, comp_1loop, comp_count
	integer :: kdir

	CONTAINS

	SUBROUTINE ReadInput()
	!
	! Subroutine for reading the input file "in.dat"
	!
	integer :: trash,i
	real(kind=8),dimension(6:48) :: mc0, zf0  ! These are used just of mc0 and zf0 are chosen
						  ! to be read from file.
	open(unit=1,file="in.dat",status="old")
	print*, "read 1"
	read(1,*) text
	read(1,*) text
	read(1,*) text, lat_size	!dimension of the lattice
	read(1,*) text, m0		!tree-level mass
	read(1,*) text, csw		!tree-level clover coefficient
	read(1,*) text, ds		!tree-level ds boundary coefficient
	read(1,*) text, zf		!tree-level zf boundary coefficient
	read(1,*) text, theta		!theta angle
	read(1,*) text, mom_deg		!Take advantage of momentum degeneracy? 1 = Yes, 0 = No
	read(1,*) text, G_act		!Gauge action: 1 = plaquette, 2 = LW
	read(1,*) text, lambda_gf	!gauge fixing parameter
	read(1,*) text, diag_g		!Take advantage of a diagonal_gauge 1 = yes, 0 = no
	read(1,*) text, deriv		!Compute temporal derivative of observable? 1 = yes, 0 = no
	read(1,*) text, m0_zf_file	!Read mc0 and zf0 from a file? 1 = yes, 0 = no  = do
	read(1,*) text
	read(1,*) text
	read(1,*) text
	print*, "read 2"
	read(1,*) text, comp_gX		! Compute a gX correlation? 1=yes, 0=no
	read(1,*) text, comp_lY		! Compute a lY correlation? 1=yes, 0=no
	read(1,*) text, comp_g1		! Compute a g1 correlation? 1=yes, 0=no
	read(1,*) text, comp_l1		! Compute a g1 correlation? 1=yes, 0=no
	read(1,*) text, comp_gVt	! Compute a gVt correlation? 1=yes, 0=no
	read(1,*) text, comp_lVt	! Compute a lVt correlation? 1=yes, 0=no
	read(1,*) text, comp_G1_4f	! Compute a G1_4f correlation? 1=yes, 0=no
	read(1,*) text, comp_G2_4f	! Compute a G2_4f correlation? 1=yes, 0=no
	read(1,*) text, comp_Gk1_4f	! Compute a Gk1_4f correlation? 1=yes, 0=no
	read(1,*) text, comp_Gk2_4f	! Compute a Gk2_4f correlation? 1=yes, 0=no
	read(1,*) text
	read(1,*) text
	read(1,*) text
	read(1,*) text, comp_tree	! Compute correlation at tree level? 1=yes, 0=no
	read(1,*) text, comp_1loop	! Compute correlation at 1loop? 1=yes, 0=no
	read(1,*) text, comp_count	! Compute counterterms? 1=yes, 0=no
	read(1,*) text
	read(1,*) text
	read(1,*) text
	print*, "read 3"
	read(1,*) text, flav1			! choose Flavor 1
	read(1,*) text, flav2			! choose Flavor 2
	read(1,*) text, flav3			! choose Flavor 3
	read(1,*) text, flav4			! choose Flavor 4
	read(1,*) text, gam_i			! Choose Gamma structure in bilinears.
	read(1,*) text, gam_i2			! Gamma 2 <NOTE: not used anymore>
	read(1,*) text, gam4f_pair		! Pair of gamma matrices for the 4 fermion operators
	read(1,*) text, kdir			! k direction for the Gk1 and Gk2 observables
	print*, "read 4"
	close(unit=1)


	call g4f_pair(gam4f_pair,gam4f_i,gam4f_i2)

	! Read the tree-level values of mc and zf from a file
	if(m0_zf_file.eq.1)then
	  open(unit=1,file="mc_zf.dat",status="old")
	   do i=6,48
              read(1,*) trash, mc0(i), zf0(i)
	   enddo
	  close(1)
	  m0 = mc0(lat_size)
	  zf = zf0(lat_size)
	endif


	END SUBROUTINE ReadInput



	SUBROUTINE PrintInput()
	!
	! Subroutine to print in screen all the relevant input parameters.
	! 
	print*, "Input parameters are"
	print*, "lat_size=", lat_size		!dimension of the lattice
	print*, "m0=", m0			!mass
	print*, "csw=", csw			!clover coefficient
	print*, "ds=", ds			!ds boundary counterterm
	print*, "zf=", zf			!zf boundary counterterm
	print*, "theta=", theta		!theta angle
	print*, "G_action=", G_act		!Gauge action
	print*, "lambda_gf=", lambda_gf	!gauge fixing parameter

	if(comp_Gk1_4f.eq.1.or.comp_Gk2_4f.eq.1)then
	    print*, "Gk1 and Gk2 direction, k=",kdir
	endif

	END SUBROUTINE PrintInput


	SUBROUTINE g4f_pair(gam4f_pair,gam4f_i,gam4f_i2)
	!
	! This subroutine contains a vector of dim 16 with all the indices for the 
	! gamma structures involved in the 4-fermion operators.
	! Note that it contains indices for the independent elements in the calculation.
	! See the module gammamatrices to see to what bilinear operator every
	! index corresponds to.
	!
	integer,intent(in) :: gam4f_pair
	integer,intent(out) :: gam4f_i,gam4f_i2

	integer,dimension(2,16) :: g4f_pair_vector
	integer :: i

	do i=1,10
	g4f_pair_vector(1,i) = i
	enddo
	g4f_pair_vector(1,11) = 12
	g4f_pair_vector(1,12) = 13
	g4f_pair_vector(1,13) = 14
	g4f_pair_vector(1,14) = 17
	g4f_pair_vector(1,15) = 18
	g4f_pair_vector(1,16) = 22

	g4f_pair_vector(2,1) = 2
	g4f_pair_vector(2,2) = 1
	g4f_pair_vector(2,3) = 7
	g4f_pair_vector(2,4) = 8
	g4f_pair_vector(2,5) = 9
	g4f_pair_vector(2,6) = 10
	g4f_pair_vector(2,7) = 3
	g4f_pair_vector(2,8) = 4
	g4f_pair_vector(2,9) = 5
	g4f_pair_vector(2,10) = 6

	g4f_pair_vector(2,11) = 28
	g4f_pair_vector(2,12) = 29
	g4f_pair_vector(2,13) = 30
	g4f_pair_vector(2,14) = 33
	g4f_pair_vector(2,15) = 34
	g4f_pair_vector(2,16) = 38

	gam4f_i = g4f_pair_vector(1,gam4f_pair)
	gam4f_i2 = g4f_pair_vector(2,gam4f_pair)

	END SUBROUTINE g4f_pair








!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	INTEGER FUNCTION mu_low(m)
	!
	! Function to select the lower limit of a mu index deppending on wether we
	! have a diagonal gauge propagator or not.
	!
	integer :: m

	if(diag_g.eq.1)then	
		mu_low = m
	else
		mu_low = 0
	endif

	END FUNCTION mu_low


	INTEGER FUNCTION mu_up(m)
	!
	! Function to select the upper limit of a mu index deppending on wether we
	! have a diagonal gauge propagator or not.
	!
	integer :: m

	if(diag_g.eq.1)then
		mu_up = m
	else
		mu_up = 3
	endif

	END FUNCTION mu_up

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





	END MODULE Input
