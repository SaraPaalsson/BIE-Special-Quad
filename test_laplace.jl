using FactCheck
using JLD

include("initialization.jl")

facts("Initialization") do
	context("Quadrature") do
		ncheck::Array{Any,1}
		wcheck::Array{Any,1}
		ncheck =  [-0.989400934991650, -0.944575023073233,-0.865631202387832,-0.755404408355003,-0.617876244402644,-0.458016777657227,-0.281603550779259,-0.095012509837637,0.095012509837638,0.281603550779259,0.458016777657228,0.617876244402644,0.755404408355003,0.865631202387831,0.944575023073233,0.989400934991650]
		wcheck = [0.027152459411754,0.062253523938648,0.095158511682492,0.124628971255533,0.149595988816577,0.169156519395003,0.182603415044924,0.189450610455069,0.189450610455068,0.182603415044924,0.169156519395003,0.149595988816577,0.124628971255534,0.095158511682493,0.062253523938647,0.027152459411755]
		(n,w) = initialization.gaussleg16(-1,1)
		@fact sum(w) --> roughly(2) "Sum of quadrature weights wrong"
		weights_ok = 0
		nodes_ok = 0
		for (i in 1:16)
			if abs(n[i]-ncheck[i]) < 0.000001
				nodes_ok += 1
			end
			if abs(w[i]-wcheck[i]) < 0.000001
				weights_ok += 1
			end
		end
		@fact weights_ok --> 16 "G.-L. weights wrong"
		@fact nodes_ok --> 16 "G.-L. nodes wrong"
	end
	context("Parametrization") do
		testz = load("indata/testdata.jld","zparm")
		testzp = load("indata/testdata.jld","zparmp")
		testzpp = load("indata/testdata.jld","zparmpp")
		tvec = linspace(0,2*pi,101)
		(z,zp,zpp) = initialization.starfish_parm(tvec)
		@fact z --> roughly(testz) "Starfish parametrization broken"
		@fact zp --> roughly(testzp) "Starfish parametrization broken, first derivative"
		@fact zpp --> roughly(testzpp) "Starfish parametrization broken, second derivative"
	end
end


#FactCheck.exitstatus()
