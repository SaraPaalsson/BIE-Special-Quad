
using JLD

include("initialization.jl")
include("readtextvariable.jl")
include("spec_quad.jl")

function laplace_main()
	println("Starting...")

	filename = "starfish"

	#Setup domain and initialize parameters
	(dropInfo) = laplace_init(filename)

	#Fill domain with computational points for error testing
	(domz,domzplot) = initialization.filldomain(filename)

	#Compute known value of u at computational point
	N = size(domz,1)
	uknown = zeros(N,1)
	for (i in 1:N)
		uknown[i] = laplace_boundcond(domz[i],dropInfo["zref"])
	end

	#Step and solve, for Laplace only once: 
	#Solve BIE to compute density
	omega = laplace_computedensity(dropInfo)
	#With omega, compute u with normal quadrature
	u = laplace_computeNormquad(omega,dropInfo,domz)
	uspec = laplace_computeSpecquad(omega,dropInfo,domz,u)


	#Post-processing, plotting etc
	laplace_post()
	return dropInfo
	println("Done!")
end

function laplace_computedensity(dropInfo)
#= Computes density by solving the BIE 
	Input: 		dropInfo, dict containing disc. points, quadrature weights etc
	Output:		omega, density array
=#
	println("Solving BIE to compute density...")

	N = dropInfo["Npanels"]*16 #Number of disc. points on interface
	W = dropInfo["W"]
	z = dropInfo["z"]; zp = dropInfo["zp"]; zpp = dropInfo["zpp"]
	A = zeros(N,N); b = zeros(N,1)
	for (i in 1:N)
		for (j in 1:N)
			if i == j
				#Limits are known for the singular kernels
				A[i,i] = 1 + W[i]*imag(zpp[i]/zp[i]/2)/pi
			else 
				A[i,j] = A[i,j] + W[j]*imag(zp[j]/(z[j]-z[i]))/pi
			end
		end
		b[i] = 2*laplace_boundcond(z[i],dropInfo["zref"])
	end 

	omega = \(A,b)
	return omega
end

function laplace_computeNormquad(omega,dropInfo,compz)
#= Compute solution to laplace's equation on domain
	Input: 		omega, density vector
				dropInfo, disc. points and quadrature weights
				compz, computational domain on which to solve Laplace's eq
	Output:		u, computed solution
=#	
	println("Computing solution on all target points using normal quadrature...")
	N = dropInfo["Npanels"]*16; W = dropInfo["W"]; z = dropInfo["z"]; zp = dropInfo["zp"]
	u = zeros(length(compz),1)
	for (i in 1:length(compz))
		for (j in 1:N)
			u[i] = u[i] + W[j]*omega[j]*imag(zp[j]/(z[j]-compz[i]))/2/pi
		end
	end
	return u
end

function laplace_post()
	println("Post processing...")
end

function laplace_init(file_in)
#= Reads in initial data and parameters from saved file file_in.jld
	Input: 		file_in, file to load 
	Output:		dropInfo, dict containing all information of interface setup
				bc, boundary condition function
=#
	println("Initializing domain...")
	tmp = "indata/"
	tmp2 = ".jld"
	filename = tmp*file_in*tmp2 #Which file to load
	dropInfo = load(filename)	

	return dropInfo
end

function laplace_boundcond(x,zp)
#= Computes boundary condition (as well as exact solution) 
	of Laplaces equation for point x
	Input: 		x, computational point (scalar)
				zp, reference point
	Output: 	u, value at x 
=#
	u = imag(x^2/(x-zp));
	return u 
end

