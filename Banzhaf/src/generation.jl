
include("io.jl")
using Combinatorics
using Random
using RandomNumbers.MersenneTwisters
Random.seed!(1)

"""
Returns the Power Set of a list L: The Power Set is a list of lists
"""
function PowerSet(L::Array{Int64,1})
	return collect(powerset(L))
end

"""
Given a number of individuals N, this function returns a Total Power Relation over the coalitions:
The PR is a list of indifference classes, each indifference class is a list of coalitions
"""
function TotalPowerRelation(N::Int64)
	L=[1:N;]
	PS=PowerSet(L)
	PS=shuffle(PS)
	s=0;
	k=ceil.(Int, 2^N * rand())
	Class=PS[s+1:s+k]
	PR=Array{Array{Any,1},1}[Class]
	s=s+k
	while (s<2^N)
		k=ceil.(Int, (2^N-s) * rand())
		Class=PS[s+1:s+k]
		PR=vcat(PR,[Class])
		s=s+k
	end
	return PR
end

"""
Generate a DataSet given the number of individuals and the number of Power Relations required
The files are saved in the folder .../data
"""
function GenerateDataSet(N::Int64,n::Int64)
	for i in 1:n
		fileName = "../data/instance_t" * string(N) * "total"* "_" * string(i) * ".txt"
		if !isfile(fileName)
            println("-- Generating file " * fileName)
            SaveInstance(TotalPowerRelation(N), fileName)
        end 
	end
end


