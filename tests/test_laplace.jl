using FactCheck
using JLD

include("../src/initialization.jl")
include("../src/laplace_main.jl")
include("../src/readtextvariable.jl")
include("../src/spec_quad.jl")

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
		testz = load("test_para.jld","zparm")
		testzp = load("test_para.jld","zparmp")
		testzpp = load("test_Para.jld","zparmpp")
		tvec = linspace(0,2*pi,101)
		(z,zp,zpp) = initialization.starfish_parm(tvec)
		@fact z --> roughly(testz) "Starfish parametrization broken"
		@fact zp --> roughly(testzp) "Starfish parametrization broken, first derivative"
		@fact zpp --> roughly(testzpp) "Starfish parametrization broken, second derivative"
	end
end

#SETTINGS FOR FOLLWING TESTS
fillInfo = Dict("r0" => 0, "rL" => 0.9999, "fillLayers" => 5, "fillparm" => 50)
dropInfo = Dict("fillInfo" => fillInfo,"zref" => 1.5+1.5*im)
plotInfo = initialization.filldomain("starfish",dropInfo)

facts("Computational domain") do
	context("Distribution of points") do
		zr = real(plotInfo["domz"]); zim = imag(plotInfo["domz"])
		testdomz = readtextvariable("testdomain_z.txt")
		N = length(testdomz)
		ztestr  = testdomz[1:N/2]; ztestim = testdomz[N/2+1:end]
		@fact zr --> roughly(ztestr); "Computational points wrong"
		@fact zim --> roughly(ztestim); "Computational points wrong"
		uk = laplace_compsol(plotInfo,dropInfo)
		uktest =  readtextvariable("testdomain_uexact.txt")
		@fact uk --> roughly(uktest); "Exact velocity computed wrongly"
	end
end

facts("Solving BIE") do
		Npanels = 10; N = 16*Npanels; zref = 1.5+1.5*im
		tmp = readtextvariable("testBIE_z.txt")
		z = tmp[1:N]+im*tmp[N+1:end]
		tmp = readtextvariable("testBIE_zp.txt")
		zp = tmp[1:N]+im*tmp[N+1:end]
		tmp = readtextvariable("testBIE_zpp.txt")
		zpp = tmp[1:N]+im*tmp[N+1:end]
		W = readtextvariable("testBIE_W.txt")
		panels = readtextvariable("testspecq_pan.txt")
		dropInfo = Dict("Npanels" => Npanels, "z" => z, "zp" => zp, "zpp" => zpp, "W" => W, "zref" => zref,"fillInfo" => fillInfo,"panels"=>panels)
		omegatest = readtextvariable("testBIE_omega.txt")
		utest = readtextvariable("testBIE_unorm.txt")

	context("Complex density") do
		omega = laplace_computedensity(dropInfo)
		@fact omega --> roughly(omegatest) "BIE density wrong"
		plotInfo = initialization.filldomain("starfish",dropInfo)
		unorm = laplace_computeNormquad(omegatest,dropInfo,plotInfo)
		@fact unorm --> roughly(utest)
	end

	context("Special quadrature") do
		uspecq_test = readtextvariable("testspecq_u.txt")
		uspecq = laplace_computeSpecquad(omegatest,dropInfo,plotInfo,utest)
		println(norm(uspecq-uspecq_test,Inf))
		@fact uspecq --> roughly(uspecq_test)
	end
end


#FactCheck.exitstatus()
