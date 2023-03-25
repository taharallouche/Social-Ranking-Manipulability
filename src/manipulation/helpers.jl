"""
PR_to_Order:
-PR:Power Relation: List of lists of lists: list of indifference classes of coalitions
Output:
-Coalitions:List of ordered coalitions
-Classes   :List of the class to which belongs each coalition (in the same order as Coalitions)
"""
function PR_to_Order(PR::Vector{Vector{Vector{Int}}})
	Classes=[]
	Coalitions=Array{Any,1}[]
	p=1
	for C in PR
		for S in C
			push!(Classes,p)
			push!(Coalitions,S)
		end
		p+=1
	end
	return Coalitions,Classes
end

function find_Coalition(Coalitions::Array{Array{Any,1},1},Coalition::Array{Any,1})
	exist=false
	for i in 1:size(Coalitions,1)
		if (Coalitions[i]==Coalition)
			exist=true
			return i
		end
	end
	if (exist=false)
		return -1
	end
end


"""
CP-Majority:
-i,j:Players to compare.
-d_ij: Number of coalitions preferring (not strictly) i to j.
-d_ji: Number of coalitions preferring j to i.
-Pairs: List of positions of the compared Sui and Suj coalitions.
"""
function CP_Maj(Coalitions::Array{Array{Any,1},1},m::Array{Int64, 2},i::Int64,j::Int64)
	d_ij=0
	d_ji=0
	Pairs=[]
	for k in 1:size(Coalitions,1)
		if !((i in Coalitions[k]) || (j in Coalitions[k])) 
			S_i=sort(vcat(Coalitions[k],i))
			S_j=sort(vcat(Coalitions[k],j))
			li=find_Coalition(Coalitions,S_i)
			lj=find_Coalition(Coalitions,S_j)
			if (li!=-1) & (lj!=-1)
				push!(Pairs,(li,lj))
				d_ij+=m[li,lj]
				d_ji+=m[lj,li]
			end
		end
	end
	return d_ij,d_ji,Pairs
end