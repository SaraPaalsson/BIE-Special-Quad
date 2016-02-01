
using JLD

function laplace_main()
	println("Starting...")

	filename = "test"

	#Setup domain and initialize parameters
	(dropInfo) = laplace_init(filename)

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

function laplace_init(file_in)
	println("Initializing domain...")
	tmp = "indata/"
	tmp2 = ".jld"
	filename = tmp*file_in*tmp2 #Which file to load
	dropInfo = load(filename)	
	return dropInfo
end
