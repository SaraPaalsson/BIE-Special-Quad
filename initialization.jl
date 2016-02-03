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

	#Parametrization		
	if geomShape == "starfish" #Starfish domain
		fparm(t) = starfish_parm(t)
	else
		error("Not defined parametrization")
	end

	#Distribute G.-L points on interface 
 	(z,W,panels,parm) = domain_setup(Npanels,fparm) #Disc. points and weights

	#Save data to file for later use
	tmp = "indata/"
	tmp2 = ".jld"
	filename = tmp * filename_in * tmp2
	save(filename,"z",z,"W",W,"Npanels",Npanels,"panels",panels,"parm",parm)
end

function starfish_parm(tvec)
#=  Parametrize interface as a starfish 
	Input:		-parametrization vec. tvec [0,2pi]
	Output:		-tau, interface points
				-taup,taupp, first and second derivative of disc. points
=#
	tau = zeros(size(tvec)); tau = complex(tau)
	taup = zeros(size(tvec)); taup = complex(taup)
	taupp = zeros(size(tvec)); taupp = complex(taupp)
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

	parm = []; z = []; W = []
	for (i in 1:Np)
		(nodes,weights) = gaussleg16(panels[i],panels[i+1])
		append!(parm,nodes)
		(tmp,tmp1,tmp2) = f(nodes)
		append!(z,tmp)
		append!(W,vec(weights))
	end
	return z,W,panels,parm
end

function gaussleg16(a,b)
#=	Compute 16 point Gauss-Legendre quadrature nodes and weights on 
interval [a,b] (as done in Trefethen) 
	Input:		-a,b, start and end of interval resp.
	Output:		-nodes, Array{Float64,1} quadrature points on interv. (a,b) (16-p G.-L.)
				-weights, Array{Float64,1} quadrature weights on interv. (a,b) (16-p G.-L.)
=#
	nvec = 1:15
	beta = 0.5*(1-(2.0*nvec).^(-2)).^(-1/2)
	Tn = diagm(beta,-1) + diagm(beta,1)
	(D,V) = eig(Tn)
	booltmp = !(D.==0)
	nodes = []
	for (i in 1:length(D))
		if booltmp[i]
			push!(nodes,D[i])
		end
	end
	weights = (2*V[1,:].^2)'
	#= NB here we could use our saved W16 weights and remap them instead... 
	but we won't save much time doing so =#
	nodes = (a*(1-nodes)+b*(1+nodes))/2
	weights =(b-a)/2*weights
	return nodes,weights
end

end

