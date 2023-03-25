include("generation.jl")
include("helpers.jl")
using Cbc
using JuMP




function MarginalContribution(Coalitions::Array{Array{Any,1},1},m::Array{Int64, 2},i::Int64)
	s_i=0
	Pairs=[]
	for k in 1:size(Coalitions,1)
		if (i in Coalitions[k])
			S_i=filter(x->x!=i,Coalitions[k])
			l=find_Coalition(Coalitions,S_i)
			if (l!=-1)
				push!(Pairs,(k,l))
				s_i+=m[k,l]-m[l,k]
			end
		end
	end
	return s_i,Pairs
end



"""
Manipulate:
-PR:Power Relation: List of lists of lists: list of indifference classes of coalitions
-i :The manipulator
-A :The set of candidates
"""
function ManipulateBanzhaf(PR::Array{Array{Array{Any,1},1},1},i::Int64,A::Array{Int64,1})
	Coalitions, Classes=PR_to_Order(PR)
	n=size(Coalitions,1)
	a=size(A,1)
	
	#Create the matrix m describing the PR
	m=Array{Int64}(undef, n, n)
	for i in 1:n
		for j in 1:n
			if (Classes[i]<=Classes[j])
				m[i,j]=1
			else
				m[i,j]=0
			end
		end
	end	
	
	# Create the model
    model = Model(with_optimizer(Cbc.Optimizer))
    
    # Variables
    @variable(model, 0<=manip[1:n, 1:n]<=1, Int) #Matrix describing the new Pw Relation
    @variable(model, s[1:a], Int)        #Vector of the new marginal contributions of candidates
    @variable(model, 0<=rk[1:a, 1:a]<=1,Int) #Matrix describing the new ranking of individuals
    
    #Constraints
    
    for k in 1:n
    	#Diagona: Reflexive PR
    	@constraint(model,manip[k,k]==1)
    	for l in 1:n
    		#The position of coalitions containing i can not improve
    		if (i in Coalitions[k]) && !(i in Coalitions[l]) 
    			@constraint(model,manip[l,k]>=m[l,k])
    			@constraint(model,manip[k,l]<=m[k,l])
    		end
    		#i cannot change the relative posiitions of other coalitions
    		if !(i in Coalitions[k]) && !(i in Coalitions[l])
    			@constraint(model,manip[l,k]==m[l,k])
    		end
    		#Transitivity
    		for t in 1:n
    			@constraint(model,manip[l,k]>=manip[l,t]+manip[t,k]-1)
    		end
    		#Maintaing the same information (PR remains total)
    		@constraint(model,manip[k,l]+manip[l,k]>=1)
    	end
    end 
	
	#Marginal contributions
	for j in 1:a
		Pairs=MarginalContribution(Coalitions,m,A[j])[2] #Pairs of coalitions with and without j
		@constraint(model,sum(manip[P[1],P[2]]-manip[P[2],P[1]] for P in Pairs)==s[j])
	end
	
	
	#Individual relation: Social ranking
	for j in 1:a
		@constraint(model,rk[i,j]>=(s[findall(x->x==i,A)[1]]-s[j])/n/2+1/4/n)
		@constraint(model,rk[i,j]<=1-(s[j]-s[findall(x->x==i,A)[1]])/n/2)
		@constraint(model,rk[j,i]>=(s[j]-s[findall(x->x==i,A)[1]])/n/2+1/4/n)
		@constraint(model,rk[j,i]<=1-(s[findall(x->x==i,A)[1]]-s[j])/n/2)
	end
	
	#Optimize the Social Score of i
	@objective(model, Max, sum(rk[i,j]-rk[j,i] for j in 1:a)  )
	
	# Start a chronometer
    start = time()

    # Solve the model
    optimize!(model)
	
	
	#DISPLAY
	#Time elapsed
	println("time elpased: ",time()-start)
	
	println("********************************************")
	println("********************************************")
	#Test if manipulable
	feasible=JuMP.primal_status(model) == JuMP.MathOptInterface.FEASIBLE_POINT
	manipulable=false
	victims=0
	potential_victims=0
	for j in 1:a
		if ((A[j]!=i) && (MarginalContribution(Coalitions,m,i)[1]==MarginalContribution(Coalitions,m,A[j])[1])) || ((A[j]!=i) &&  (MarginalContribution(Coalitions,m,i)[1]<MarginalContribution(Coalitions,m,A[j])[1]))
			potential_victims+=1
		end
	end
	if feasible
		if (JuMP.objective_value(model)>sum((MarginalContribution(Coalitions,m,i)[1]>MarginalContribution(Coalitions,m,j)[1])-(MarginalContribution(Coalitions,m,i)[1]<MarginalContribution(Coalitions,m,j)[1]) for j in 1:a))
            manipulable=true
		end
		for j in 1:a
			if ((A[j]!=i) && (JuMP.value(s[j])<JuMP.value(s[i])) && (MarginalContribution(Coalitions,m,i)[1]==MarginalContribution(Coalitions,m,A[j])[1])) || ((A[j]!=i) && (JuMP.value(s[j])<=JuMP.value(s[i])) && (MarginalContribution(Coalitions,m,i)[1]<MarginalContribution(Coalitions,m,A[j])[1]))
				victims+=1
				println("Manipulable for: ",i," against ",A[j])
			end
		end
		NonInitialementPremier=false
		for j in 1:a
			if (MarginalContribution(Coalitions,m,i)[1]<MarginalContribution(Coalitions,m,A[j])[1])
				NonInitialementPremier=true
			end
		end
		Premier=true
		for j in 1:a
			if (JuMP.value(s[i])<JuMP.value(s[j]))
				Premier=false
			end
		end
		DevenuPremier=(NonInitialementPremier)&&(Premier)
	end
	if !manipulable
		println("Not manipulable for candidate ",i)
	end
	
	#Marginal Contribution
	println("********************************************")
	println("********************************************")
	println("New Marginal Contributions:")
	for j in 1:a
		println("s(",A[j],")=",JuMP.value(s[j]))
	end
	
	println("********************************************")
	println("********************************************")
	#Matrix of the new Power Relation:
	println("Matrix of the new Power Relation:")
	for k in 1:n
		for l in 1:n
			print(JuMP.value(manip[k,l])," , ")
		end
		println()
	end
	if manipulable
		return manipulable,victims/potential_victims,NonInitialementPremier,DevenuPremier
	else
		return manipulable,0
	end
end

