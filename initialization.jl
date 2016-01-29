using JLD

function indata_setup()
#Computes initialization data and settings 

#######################
	#Parametrization	
#######################
# Starfish domain
	tau(t) = (1+0.3*cos(5*(t))).*exp(im*(t))
	taup(t) = (-1.5*sin(5*(t))+im*(1+0.3*cos(5*(t)))).*exp(im*(t));
	taupp(t) = exp(im*(t)).*(-1-7.8*cos(5*(t))-(3*im)*sin(5*(t)));

	Npanels = 10 #Number of GL panels on interface

 	(z,W) = domain_setup(Npanels,tau) #Disc. points and weights

	s = linspace(0,2*pi,Npanels*16+1)
	s = s[1:end-1]
 	z1 = tau(s)


#######################
	#Save data to file for later use
#######################
	filename = "indata/test.jld"
	save(filename,"Npanels",Npanels,"z",z,"s",s)
end

function domain_setup(Np,f)
#Compute disc.points and their quadrature weights for parametrization f

	panels = linspace(0,2*pi,Np+1) #Divide interface into panels

	for (i in 1:Np)
		(nodes,weights)
	end
	return z,W
end

