function laplace_computeSpecquad(omega,dropInfo,domz,u)

	println("Compute special qudrature modifications...")
# 	# Load interpolation coeff
 	IP1 = readtextvariable("IP1.txt")
 	IP2 = readtextvariable("IP2.txt")
 	W16 = readtextvariable("W16.txt")
 	W32 = readtextvariable("W32.txt")

 	Npanels = dropInfo["Npanels"]
 	N = Npanels*16

 	uspec = u #Spec.quad modifies normal u

 	# Calculate mid points and (complex) lengths of all panels
 	panels = dropInfo["panelsz"]
 	z = dropInfo["z"]
 	zp = dropInfo["zp"]
 	W = dropInfo["W"]

 	mid = zeros(Complex{Float64},Npanels) 
 	len = zeros(Complex{Float64},Npanels)
 	for (i in 1:Npanels)
 		mid[i] = (panels[i+1]+panels[i])/2
 		len[i] = panels[i+1]-panels[i]
 	end

	tz = zeros(Complex{Float64},16); tz32 = zeros(Complex{Float64},32)
 	tzp = zeros(Complex{Float64},16); tzp32 = zeros(Complex{Float64},32)
 	nzpan = zeros(Complex{Float64},16); nzpan32 = zeros(Complex{Float64},32)
 	p32 = zeros(Complex{Float64},32)
 	tom = zeros(Complex{Float64},16); tom32 = zeros(Complex{Float64},32)
 	tW = zeros(Complex{Float64},16); tW32 = zeros(Complex{Float64},32)
 	orig32 = zeros(Complex{Float64},32)

 	Ncomp = size(domz,1)
 	for (zind in 1:Ncomp) 	#Parse through all comp.points
 	 	for (panind in 1:Npanels) #For each point, parse through all panels
 	 		if abs(domz[zind]-mid[panind]) < abs(len[panind]) 
 	 		#Check if comp. point too close to panel

 	 			#Map our zk to nz (z0 in paper)
 	 			nz = 2*(domz[zind]-mid[panind])/len[panind]
 	 			oldsum = 0; testsum = 0
 	 			lg1 = log(1-nz)
 	 			lg2 = log(-1-nz)

				#Map all nodes on panel panelind to [-1,1]
 	 			for (j in 1:16) 
 	 				tz[j] = z[(panind-1)*16+j]
 	 				nzpan[j] = 2*(tz[j]-mid[panind])/len[panind]
 	 			end

 	 			#Check if the point nz is between the panel and the real axis
 	 			if real(nz) > -1 && real(nz) < 1
 	 				if imag(nz) > 0 #Above real axis
 	 					furthercheck = 0
 	 					for (j in 1:16) #Go through all nodes on panel
 	 						if imag(nzpan[j]) > imag(nz)
 	 							furthercheck = 1
 	 							break
 	 						end
 	 					end
 	 					if furthercheck == 1
 	 					#Interpol. nzpan to poly and check value for real(nz)
 	 						tmpT = real(nzpan)
 	 						tmpb = imag(nzpan)

 	 						p = vandernewtonT(tmpT,tmpb)
 	 						test = 0
 	 						for (j in 0:15)
 	 							test += p[j+1]*real(nz)^j
 	 						end
 	 						if test > imag(nz)
 	 						#Correct the value of our p0 integral
 	 							lg1 += - pi*im
 	 							lg2 += pi*im
 	 						end
 	 					end
 	 				elseif imag(nz) < 0 #Below real axis
 	 					furthercheck = 0
 	 					for (j in 1:16) #Go through all nodes on panel 
 	 						if imag(nzpan[j]) < imag(nz)
 	 							furthercheck = 1
 	 							break
 	 						end
 	 					end
 	 					if furthercheck == 1
 	 						tmpT = real(nzpan)
 	 						tmpb = imag(nzpan)
 	 						p = vandernewtonT(tmpT,tmpb)
 		 					test = 0
 	 						for (j in 0:15)
 	 							test += p[j+1]*real(nz)^j
 	 						end
 	 					end
 	 					if test < imag(nz) 
 	 					#Correct value of integral for p0
 	 						lg1 += pi*im
 	 						lg2 += -pi*im
						end
 	 				end
 	 			end

	 	 		#Calculate p0 exact
 		 		p32[1] = lg1 - lg2
 	 			for (j in 1:16)
	 	 		#Compute old contribution to u from panel
 		 			tzp[j] = zp[(panind-1)*16+j]
 	 				tom[j] = omega[(panind-1)*16+j]
 	 				tW[j] = W[(panind-1)*16+j]
 		 			oldsum += tW[j]*tom[j]*imag(tzp[j]/(tz[j]-domz[zind]))/2/pi
 	 				testsum += tW[j]*tzp[j]/(tz[j]-domz[zind])
 	 			end

	 	 		#Compare testsum to exact value of p0
 		 		if abs(p32[1]-testsum) > 1e-13
				#Standard 16-GL not good enough! Interp. to 32-point GL quad.
					tom32 = IPmultR(tom,IP1,IP2)
					tz32 = IPmultR(tz,IP1,IP2)
					tzp32 = IPmultR(tzp,IP1,IP2)
					plen = tW[1]/W16[1]
					o32sum = 0
					for (j in 1:32)
						tW32[j] = W32[j]*plen
						orig32[j] = tW32[j]/(tz32[j]-domz[zind])
						o32sum += tzp32[j]*orig32[j]
					end

					if abs(o32sum-p32[1]) < 1e-13 
						#32 point G.-L. suffices!

						newsum = 0
						for (j in 1:32)
							newsum += tW32[j]*tom32[j]*imag(tzp32[j]/(tz32[j]-domz[zind]))/2/pi
						end
						uspec[zind] += newsum - oldsum
					else
					# 32 point G.-L. not enough, use interpolatory quadrature

						nzpan32 = IPmultR(nzpan,IP1,IP2)
						signc = -1
						for (j in 1:31)
							p32[j+1] = nz*p32[j]+(1-signc)/(j)
							signc = -signc
						end
						c32 = vandernewtonT(nzpan32,tom32)
						newsum = 0
						for (j in 1:32)
							newsum += imag(p32[j]*c32[j])/pi/2
						end

						p32coeff = vandernewton(nzpan32,p32)
						newsum2 = 0
						for (j in 1:32)
							newsum2 += imag(p32coeff[j]*tom32[j])/2/pi
						end

						uspec[zind] += newsum2 - oldsum		
					end			
				end
 	 		end
 	 	end
 	end
	return uspec
end

function IPmultR(f16,IP1,IP2)
	f32 = zeros(Complex{Float64},32)
	for (i in 1:16)
	   	t1 = 0
   		t2 = 0
   		ptr = i
   		for (j in 1:8)
      		t1 += IP1[ptr]*(f16[j]+f16[17-j])
      		t2 += IP2[ptr]*(f16[j]-f16[17-j])
      		ptr += 16
	  	end
   		f32[i] = t1+t2
   		f32[33-i] = t1-t2
	end
	return f32
end



function vandernewton(T,b)
	for (k in 1:31)
    	for (i in 32:-1:k+1)
	        b[i] = b[i] - T[k]*b[i-1]
   	 	end
	end
	for (k in 31:-1:1)
    	for (i in k+1:32)
        	b[i] = b[i]/(T[i]-T[i-k])
    	end
    	for (i in k:31)
        	b[i] = b[i] - b[i+1]
    	end
	end
	return b
end

function vandernewtonT(T,b)
	x = T
	c = b
	for (k in 1:15)
    	for (i in 16:-1:k+1)
        	c[i] = (c[i]-c[i-1])/(x[i]-x[i-k])
    	end
	end
	a = c
	for k=15:-1:1
    	for i=k:15
        	a[i] = a[i]-x[k]*a[i+1];
    	end
	end
	return a
end
