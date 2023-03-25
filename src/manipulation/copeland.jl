include("generation.jl")
include("helpers.jl")
using Cbc
using JuMP


"""
Copeland:
-i:Player
-D:matrix of CP-Majority (Dij=dij)
-Copeland: Copeland score of player i
"""
function CopelandScore(i::Int64,D::Array{Int64, 2})
	Copeland=0
	a=size(D,1)
	for j in 1:a
		Copeland+=(D[i,j]>D[j,i])-(D[j,i]>D[i,j])
	end
	return Copeland
end


"""
Manipulate:
-PR:Power Relation: List of lists of lists: list of indifference classes of coalitions
-i :The manipulator
-A :The set of candidates
"""


function ManipulateCopeland(PR::Array{Array{Array{Any,1},1},1},i::Int64,A::Array{Int64,1})
	Coalitions=PR_to_Order(PR)[1]
	Classes=PR_to_Order(PR)[2]
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
	
	#Calculate the initial CP-Maj and Copeland Outcome
	#CP-Maj function
	#Copeland function
	
	
	
	# Create the model
    model = Model(with_optimizer(Cbc.Optimizer))
    
    # Variables
    @variable(model, 0<=manip[1:n, 1:n]<=1, Int) #Matrix describing the new Pw Relation
    @variable(model, 0<=D[1:a, 1:a],Int) #Matrix of the CP-Maj D=(d_ij)
    @variable(model, 0<=CP[1:a, 1:a]<=1,Int) #Matrix of the CP-Maj outcome: CP_ij=1 <=> i >=CP j
    @variable(model, Copeland[1:a],Int) #Vector of the Copeland Scores
    @variable(model, 0<=rk[1:a, 1:a]<=1,Int) #Matrix describing the new ranking of individuals
    
    #Constraints
    
    for k in 1:n
    	#Diagonal
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
    		#Maintaing the same information
    		@constraint(model,manip[k,l]+manip[l,k]>=1)
    	end
    end 
	
	#CP-Majority votings: D Matrix
	for j in 1:a
		for k in 1:a
			Pairs=CP_Maj(Coalitions,m,j,k)[3]
			@constraint(model,sum(manip[P[1],P[2]] for P in Pairs)==D[j,k])
			@constraint(model,sum(manip[P[2],P[1]] for P in Pairs)==D[k,j])
		end	
	end
	
	#CP-Majority outcome: CP Matrix
	for j in 1:a
		for k in 1:a
			@constraint(model,CP[j,k]>=(D[j,k]-D[k,j]+1)/n/2+1/n/4)
			@constraint(model,CP[j,k]<=1+(D[j,k]-D[k,j])/n/2)
			@constraint(model,CP[k,j]>=(D[k,j]-D[j,k]+1)/n/2+1/n/4)
			@constraint(model,CP[k,j]<=1+(D[k,j]-D[j,k])/n/2)
		end	
	end
	
	#Copeland Scores: Copeland Vector
	for j in 1:a
		@constraint(model,Copeland[j]==sum(CP[j,k]-CP[k,j] for k in 1:a))
	end
	
	#Individual position of i is maintained
	#Calculate the initial Copeland Scores
	D_initial=Array{Int64}(undef,a,a)
	for j in 1:a
		for k in 1:a
			D_initial[j,k]=CP_Maj(Coalitions,m,j,k)[1]
			D_initial[k,j]=CP_Maj(Coalitions,m,j,k)[2]
		end
	end
	Copeland_initial=[CopelandScore(j,D_initial) for j in 1:a]

	#Individual relation
	for j in 1:a
		@constraint(model,rk[i,j]>=(Copeland[findall(x->x==i,A)[1]]-Copeland[j]+1)/n/2)
		@constraint(model,rk[i,j]<=1-(Copeland[j]-Copeland[findall(x->x==i,A)[1]])/n/2)
		@constraint(model,rk[j,i]>=(Copeland[j]-Copeland[findall(x->x==i,A)[1]]+1)/n/2)
		@constraint(model,rk[j,i]<=1-(Copeland[findall(x->x==i,A)[1]]-Copeland[j])/n/2)
	end
	
	#Optimize the individual ranking of i
	@objective(model, Max, sum(rk[i,j]-rk[j,i] for j in 1:a)  )
	
	# Start a chronometer
    start = time()

    # Solve the model
    optimize!(model)
	#DISLPAY
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
		if ((A[j]!=i) && (Copeland_initial[i]==Copeland_initial[j])) || ((A[j]!=i) &&  (Copeland_initial[i]<Copeland_initial[j]))
			potential_victims+=1
		end
	end
	
	if feasible
		if (JuMP.objective_value(model)>sum((Copeland_initial[i]>Copeland_initial[j])-(Copeland_initial[j]>Copeland_initial[i]) for j in 1:a))
			manipulable=true
		end
		for j in 1:a
			if ((A[j]!=i) && (JuMP.value(Copeland[j])<JuMP.value(Copeland[i])) && (Copeland_initial[i]==Copeland_initial[j])) || ((A[j]!=i) && (JuMP.value(Copeland[j])<=JuMP.value(Copeland[i])) && (Copeland_initial[i]<Copeland_initial[j]))
				
				victims+=1
				println("Manipulable for: ",i," against ",A[j])
			end
		end
	end
		NonInitialementPremier=false
		for j in 1:a
			if (Copeland_initial[i]<Copeland_initial[j])
				NonInitialementPremier=true
			end
		end
		Premier=true
		for j in 1:a
			if (JuMP.value(Copeland[i])<JuMP.value(Copeland[j]))
				Premier=false
			end
		end
		DevenuPremier=(NonInitialementPremier)&&(Premier)
	if !manipulable
		println("Not manipulable for candidate ",i)
	end
	println("********************************************")
	println("********************************************")
	#Copeland Scores
	for j in 1:a
		println("Copeland of ",j,": ",JuMP.value(Copeland[j]))
	end
	println("********************************************")
	println("********************************************")
	#Matrix of CP-Maj
	println("MAtrix D")
	for j in 1:a
		for k in 1:a
			print(JuMP.value(D[j,k])," , ")
		end
		println()
	end
	println("********************************************")
	println("********************************************")
	#Matrix of CP-Maj outcome
	println("MAtrix CP")
	for j in 1:a
		for k in 1:a
			print(JuMP.value(CP[j,k])," , ")
		end
		println()
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

