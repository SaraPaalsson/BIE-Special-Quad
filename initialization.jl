#initialization.jl


#=
This is my initialization module.
=#
module initialization
export setup, starfish_parm, gaussleg16

using JLD


function setup(Npanels,geomShape,filename_in)
#= Compute initialization data and settings 
	Input:		-Npanels, number of panels on interface
				-geomShape, desired shape of interface (so far only starfish allowed) 
	Saves parameters, settings and initial disc. points to indata/filename.jld
=#

	#Reference point for solution
	z_ref = 1.5+1.5*im;

	#Parametrization		
	if geomShape == "starfish" #Starfish domain
		fparm(t) = starfish_parm(t)
	else
		error("Not defined parametrization")
	end

	#Distribute G.-L points on interface 
 	(z,W,panels,panelsz,parm,zp,zpp) = domain_setup(Npanels,fparm) #Disc. points and weights

	#Save data to file for later use
	tmp = "indata/"
	tmp2 = ".jld"
	filename = tmp * filename_in * tmp2
	save(filename,"z",z,"W",W,"Npanels",Npanels,"panels",panels,"panelsz",panelsz,"parm",parm,"zref",z_ref,"zp",zp,"zpp",zpp)
end

function filldomain(geomShape)
#= Fills the domain described with compuational points. 
Discretize the radius r and parametriation t.
	Input: 		geomShape, domain type to discretize
	Output: 	z, computational points
				Zplot, points z but distributed on matrix form for pcolor plotting
=#
	println("Filling the domain with computational points...")
	#Parametrization		
	if geomShape == "starfish" #Starfish domain
		f(t) = initialization.starfish_parm(t)
	else
		error("Not defined parametrization")
	end


	fillLayers = 10; #number of layers to fill the domain with
	r = linspace(0,0.9999,fillLayers)

	fillparm = 12; #Number of parametrization points on each layer
	t = linspace(0,2*pi,fillparm)

	(tmp,tmp2,tmp3) = f(t)
	zcomp = zeros(Complex{Float64},fillLayers*fillparm)
	R = zeros(fillparm,fillLayers); T = zeros(fillparm,fillLayers)
	Zplot = zeros(Complex{Float64},fillparm,fillLayers)
	for (i in 1:fillLayers)
		for (j in 1:fillparm)
			zcomp[(i-1)*fillparm+j] = r[i]*tmp[j]
			Zplot[j,i] = r[i]*tmp[j]   
		end
	end
	return zcomp,Zplot
end




function starfish_parm(tvec)
#=  Parametrize interface as a starfish 
	Input:		-parametrization vec. tvec [0,2pi]
	Output:		-tau, interface points
				-taup,taupp, first and second derivative of disc. points
=#
	tau = zeros(Complex{Float64},size(tvec))
	taup = zeros(Complex{Float64},size(tvec))
	taupp = zeros(Complex{Float64},size(tvec))
	for (k in 1:size(tvec,1))
		t = tvec[k]
		tau[k] = (1+0.3*cos(5*t)).*exp(im*t)
		taup[k] = (-1.5*sin(5*(t))+im*(1+0.3*cos(5*(t)))).*exp(im*(t));
		taupp[k] = exp(im*(t)).*(-1-7.8*cos(5*(t))-(3*im)*sin(5*(t)));
		k += 1
	end
	return tau,taup,taupp
end


function domain_setup(Np,f)
#= Compute disc.points and their quadrature weights for 
parametrization f
	Input:		-Np, number of G.-L panels
				-f, parametrization function
	Output:		-z, disc. points (16 point G.-L)
				-W, quadrature weights (16 point G.-L)
				-panels, panel division in [0,2pi]
				-parm, parametrization vector [0,2pi]
=#
	panels = linspace(0,2*pi,Np+1) #Divide interface into panels

	parm = []; z = []; W = []; zp = []; zpp = []
	for (i in 1:Np)
		(nodes,weights) = gaussleg16(panels[i],panels[i+1])
		append!(parm,nodes)
		(tmp,tmp1,tmp2) = f(nodes)
		append!(z,tmp)
		append!(zp,tmp1)
		append!(zpp,tmp2)
		append!(W,vec(weights))
	end

	(panelsz,tmp1,tmp2) = f(panels)

	return z,W,panels,panelsz,parm,zp, zpp
end

function gaussleg16(a,b)
#=	Compute 16 point Gauss-Legendre quadrature nodes and weights on 
interval [a,b] (as done in Trefethen) 
	Input:		-a,b, start and end of interval resp.
	Output:		-nodes, Array{Float64,1} quadrature points on interv. (a,b) (16-p G.-L.)
				-weights, Array{Float64,1} quadrature weights on interv. (a,b) (16-p G.-L.)
=#
	beta = zeros(1,15)
	for (nvec in 1:15)
		beta[nvec] = 0.5*(1-(2.0*nvec)^(-2))^(-1/2)
	end
	beta = vec(beta)
	Tn = diagm(beta,-1) + diagm(beta,1)
	(D,V) = eig(Tn)
	booltmp = !(D.==0)
	nodes = []
	for (i in 1:length(D))
		if booltmp[i]
			push!(nodes,D[i])
		end
	end
	weights = zeros(16,1)
	for (i in 1:16)
		weights[i] = 2*V[1,i]^2
	end
	#= NB here we could use our saved W16 weights and remap them instead... 
	but we won't save much time doing so =#
	nodes = (a*(1-nodes)+b*(1+nodes))/2
	weights =(b-a)/2*weights
	return nodes,weights
end

end

