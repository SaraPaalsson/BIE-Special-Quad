
using JLD

function laplace_main()
	println("Starting...")

	#Setup domain and initialize parameters
	(z,T) = laplace_init()
	println(z)
	println(T)
	

	#Step and solve, for Laplace only one
	laplace_computedensity()
	laplace_computetarget()

	#Post-processing, plotting etc
	laplace_post()
	
	println("Done!")
end

function laplace_computedensity()
	println("Computing density...")
end

function laplace_computetarget()
	println("Computing solution on all target points...")
end

function laplace_post()
	println("Post processing...")
end

function laplace_init()
	println("Initializing domain...")
	filename = "indata/test.jld" #Which file to loadla
	dropInfo = load(filename)	
	return dropInfo["z"], dropInfo["T"]
end
