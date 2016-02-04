function readtextvariable(filename)
#= Returns values in file as variables in vector 
	Input: 		filename, file to be read
	Output: 	var, array containing variables from file
=#
	f = open(filename)
	lines = readlines(f)
	N = size(lines,1)
	var = zeros(N,1)
	for (i in 1:N)
		var[i] = parse(Float64,lines[i])
	end
	return var
end