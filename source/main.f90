


	PROGRAM Main
	!
	!  Main program. Contains the external flow of the program. It
	!  calls sequentially a set of routines to read input files, set up the 
	!  details of the execution, and finally calls Compute_stuff,
	!  which controls the computation.
	!
	USE GammaMatrices
	USE Input
	USE Parameters
	USE Propagators	
	USE Compute


	implicit none


	!Initialization of the setup
	call ReadInput()
	call TemporalExtent(lat_size,deriv)
    	call GammaStructureConstruction()
	call flav_counter_init()
	call momentum_degeneracy(lat_size,mom_deg)

	call PrintInput()

	call Compute_stuff(lat_size,gam_i,gam_i2,gam4f_i,gam4f_i2,flav1,flav2,flav3,flav4)	
	



	END PROGRAM Main


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






