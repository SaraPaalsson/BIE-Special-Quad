
using JLD

include("initialization.jl")
include("readtextvariable.jl")

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

function laplace_computeSpecquad(omega,dropInfo,domz,u)

	println("Compute special qudrature modifications...")
# 	# Load interpolation coeff
 	IP1 = readtextvariable("IP1.txt")
 	IP2 = readtextvariable("IP2.txt")

 	Npanels = dropInfo["Npanels"]
 	N = Npanels*16

 	uspec = u #Spec.quad modifies normal u

 	# Calculate mid points and (complex) lengths of all panels
 	panels = dropInfo["panelsz"]
 	z = dropInfo["z"]
 	mid = zeros(Npanels,1); len = zeros(Npanels,1)
 	mid = complex(mid); len = complex(len)
 	for (i in 1:Npanels)
 		mid[i] = (panels[i+1]+panels[i])/2
 		len[i] = panels[i+1]-panels[i]
 	end


 	tz = zeros(16,1); tz = complex(tz)
 	nzpan = zeros(16,1); nzpan = complex(nzpan)

 	Ncomp = size(domz,1)
 	for (zind in 1:Ncomp) 	#Parse through all comp.points
 	 	for (panind in 1:Npanels) #For each point, parse through all panels
 	 		if abs(domz[zind]-mid[panind]) < abs(len[panind]) 
 	 		#Check if comp. point too close to panel

 	 		#Map our zk to nz (z0 in paper)
 	 		nz = 2*(domz[zind]-mid[panind])/len[panind]
 	 		oldsum = 0
 	 		lg1 = log(1-nz)
 	 		lg2 = log(-1-nz)

 	 		for (j in 1:16) #map all nodes on panel panelind to [-1,1]
 	 			tz[j] = z[(panind-1)*16+j]
 	 			nzpan[j] = 2*(tz[j]-mid[panind])/len[panind]
 	 		end


 	 		end
 	 	end
 	end





# for i = 1:nbr_z %Go through all points z   
#         for k=1:Npanels
#             if abs(z(i)-mid2(k)) < abs(len2(k)) %Check if z too close to any panel
#                 % % Calculate p0
#                 testaavst = testaavst + 1;
#                 nz = 2*(z(i)-mid2(k))/len2(k); %map zk to nz (z0 in paper)
#                 oldsum = 0; testsum = 0;
#                 lg1 = log(1-nz);
#                 lg2 = log(-1-nz);
#                 for j =1:16 %map all nodes on the panel to [-1,1]
#                     tz(j) = zDrops((k-1)*16+j);
#                     nzpan(j) = 2*(tz(j)-mid2(k))/len2(k);  
#                 end
#                 %Check if the point nz is between the panel and the real axis
#                 if real(nz) > -1 && real(nz) < 1
#                     if imag(nz) > 0 %above real axis, check if enclosed by axis and panel
#                         furthercheck = 0;
#                         for j=1:16 %go through all nodes on panel
#                             if imag(nzpan(j)) > imag(nz)
#                                     furthercheck = 1;
#                                     break;
#                             end
#                         end
#                         if furthercheck
#                             %interpol. nzpan to poly and check value for
#                             %real(nz)
#                             tmpT = real(nzpan);
#                             tmpb = imag(nzpan);

#                             p = vandernewtonT(tmpT,tmpb,16);
#                             test = 0;
#                                 for kk=0:15
#                                     test = test + p(kk+1)*real(nz)^kk;
#                                 end
#                             if test > imag(nz) %Correct value of integral
# %                                 lg1 = lg1 + pi*1i;
# %                                 lg2 = lg2 - pi*1i;
#                                 lg1 = lg1 - pi*1i; %HMMM???
#                                 lg2 = lg2 + pi*1i;
#                             end
#                         end
#                     else if imag(nz) < 0 %below the real axis, check enclosed
#                             furthercheck = 0;
#                             for j=1:16 %go through all nodes on panel
#                                 if imag(nzpan(j)) < imag(nz)
#                                     furthercheck = 1;
#                                     break;
#                                 end
#                             end
#                             if furthercheck
#                                 tmpT = real(nzpan);
#                                 tmpb = imag(nzpan);
  
#                             p = vandernewtonT(tmpT,tmpb,16);
#                             test = 0;
#                                 for kk=0:15
#                                     test = test + p(kk+1)*real(nz)^kk;
#                                 end


#                                 if test < imag(nz) %Correct value of integral
# %                                     lg1 = lg1 - pi*1i;
# %                                     lg2 = lg2 + pi*1i;
#                                     lg1 = lg1 + pi*1i; %HMMM??
#                                     lg2 = lg2 - pi*1i;
#                                 end
#                             end
#                         end
#                     end
#                 end
#                 p32 = zeros(32,1);
#                 p32(1) = lg1-lg2;
#                 % % Calculate old contribution to u from panel
#                 for j=1:16 %save all info from nodes on panel
#                     tzp(j) = zpDrops((k-1)*16+j);
#                     tmu(j) = mudens((k-1)*16+j);
#                     tW(j) = wDrops((k-1)*16+j);
#                     oldsum = oldsum + 1/(2*pi)*tW(j)*tmu(j)*imag(tzp(j)/(tz(j)-z(i)));
#                     testsum = testsum + tW(j)*tzp(j)/(tz(j)-z(i)); %num.approx p0
#                 end
#                 if abs(p32(1)-testsum) > 1e-13 %Standard 16-GL not good enough!
#                     % % Interpolate to 32-point GL quadrature
#                     tmu32 = IPmultR(tmu,IP1,IP2);
#                     tz32 = IPmultR(tz,IP1,IP2);
#                     tzp32 = IPmultR(tzp,IP1,IP2);
#                     plen = tW(1)/W16(1);
#                     o32sum = 0;
#                     for j=1:32 %Approximate new p0
#                         tW32(j) = W32(j)*plen;
#                         orig32(j) = tW32(j)/(tz32(j)-z(i));
#                         o32sum = o32sum + tzp32(j)*orig32(j);
#                     end
#                     if abs(o32sum-p32(1)) < 1e-13 %32 GL suffices!
                                  
#                         z32(i) = 1;
                        
#                         newsum = 0;
#                         for j=1:32
#                             newsum = newsum + 1/(2*pi)*tW32(j)*tmu32(j)*...
#                                 imag(tzp32(j)/(tz32(j)-z(i)));
#                         end
#                         u_spec(i) = u_spec(i) + (newsum-oldsum);
#                     else %32 GL not enough, use interpolatory quadrature instead
#                         % Use interpolatory quadrature
                        
#                         nzpan32 = IPmultR(nzpan,IP1,IP2);
#                         signc = -1;
#                         for j=1:31 %Calculate pk:s...
#                             p32(j+1) = nz*p32(j) + (1-signc)/(j);%(1-signc)/j;
#                             signc = -signc;
#                         end
#                         c32 = vandernewtonT(nzpan32,tmu32,32);
#                         newsum = 0;
#                         for j=1:32
#                             newsum = newsum + 1/(2*pi) * imag(p32(j)*c32(j));
#                         end
                        
#                         p32coeff = vandernewton(nzpan32,p32,32);
#                         newsum2 = 0;
#                         for j=1:32
#                             newsum2 = newsum2 + 1/(2*pi) * imag(p32coeff(j)*tmu32(j));
#                         end
                        
                        
#                         u_spec(i) = u_spec(i) + (newsum2-oldsum);
                        
#                         zinter(i) = 1;
                       
#                     end
                    
#                 else    
#                     z16(i) = 1;
#                 end
                
                
                
#             end
#         end

# end





	return 1
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

