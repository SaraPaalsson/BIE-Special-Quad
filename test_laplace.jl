using FactCheck

include("initialization.jl")

facts("Initialization") do
	context("Quadrature") do
		@fact 2 --> 2 "not two"
	end
end

facts("Parametrization") do
	@fact 3 --> 4
end

FactCheck.exitstatus()
