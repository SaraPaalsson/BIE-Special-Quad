=================
Initialization
=================
To initialize our computations, the computational domain will be set up and parametrized. This setup needs to be done once for each test case and is then saved in **indata/filename**.  This is done by *indata_setup*: ::

	function indata_setup(Npanels,geomShape,filename_in)
	#= Compute initialization data and settings 
		Input:	-Npanels, number of panels on interface
			-geomShape, desired shape of interface (so far only starfish allowed) 
		Saves parameters, settings and initial disc. points to indata/filename.jld
	=#

The starfish domain is represented by a function computing the interface points for a given parametrization, as well as its first and second order derivative. 

To compute the discretized points as well as the quadrature weights for the whole interface, *domain_setup* is used: ::
	
	function domain_setup(Np,f)
	#= 	Compute disc.points and their quadrature weights for 
		parametrization f
		Input:		-Np, number of G.-L panels
				-f, parametrization function
		Output:		-z, disc. points (16 point G.-L)
				-W, quadrature weights (16 point G.-L)
				-panels, panel division in [0,2pi]
				-parm, parametrization vector [0,2pi]
	=#

Here, the panels are constructed. On each panel, 16 Gauss-Legendre points are computed together with their weights by *gaussleg16*: ::

	function gaussleg16(a,b)
	#=	Compute 16 point Gauss-Legendre quadrature nodes and weights on 
		interval [a,b] (as done in Trefethen) 
		Input:		-a,b, start and end of interval resp.
		Output:		-nodes, quadrature points on interv. (a,b) (16-p G.-L.)
				-weights, quadrature weights on interv. (a,b) (16-p G.-L.)
	=#

To fill the domain with computational points on which to evaluate the error with normal and special quadrature respectively, *filldomain* is used: ::
	
	function filldomain(geomShape)
	#= Fills the domain described with compuational points. 
	Discretize the radius r and parametriation t.
		Input: 		geomShape, domain type to discretize
		Output: 	z, computational points
				Zplot, points z but distributed on matrix form for pcolor plotting
	=#

Needed packages for these computations are: 

* `JLD <https://github.com/JuliaLang/JLD.jl>`_
