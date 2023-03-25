include("helpers.jl")
include("generation.jl")
using Cbc
using JuMP


"""
Compute the Kramer Simpson Score:
Inputs:
-Individual i
-The matrix D (CP-Maj)
Output:
-The Kramer Simposon score of i
"""
function KramerSimpsonScore(i::Int64,D::Array{Int64,2})
	a=size(D,1)
	Score=maximum(filter(x->x<16,D[:,i])) ### Attention!!!! x<M, M=2^(N-1)
	return Score 
end


"""
Manipulate3: Manipulate in case N=3
-PR:Power Relation: List of lists of lists: list of indifference classes of coalitions
-i :The manipulator
-A :The set of candidates
"""
function Manipulate3(PR::Array{Array{Array{Any,1},1},1},i::Int64,A::Array{Int64,1})
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
    @variable(model, Kramer[1:a],Int) #Vector of the Kramer-Simpson Scores
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
			@constraint(model,CP[j,k]>=(D[j,k]-D[k,j]+1)/n/2)
			@constraint(model,CP[j,k]<=1+(D[j,k]-D[k,j])/n/2)
			@constraint(model,CP[k,j]>=(D[k,j]-D[j,k]+1)/n/2)
			@constraint(model,CP[k,j]<=1+(D[k,j]-D[j,k])/n/2)
		end	
	end
	
	#Kramer-Simpson Score for player 1
	@variable(model, d21_d31P>=0, Int) #Partie positive de d21-d31
	@variable(model, d21_d31N>=0, Int) #Partie négative de d21-d31
	@variable(model, 0<=Sign1<=1, Int) #Vaut 1 ssi d21-d31>=0
	@constraint(model,D[2,1]-D[3,1]==d21_d31P-d21_d31N)
	@constraint(model,Sign1>=(D[2,1]-D[3,1])/4)
	@constraint(model,Sign1<=1+(D[2,1]-D[3,1])/4)
	@constraint(model,d21_d31P<=2*Sign1)
	@constraint(model,d21_d31N<=2*(1-Sign1))
	@constraint(model,Kramer[1]-0.5*(D[2,1]+D[3,1])==0.5*(d21_d31P+d21_d31N))
	#Kramer-Simpson Score for player 2
	@variable(model, d12_d32P>=0, Int) #Partie positive de d21-d31
	@variable(model, d12_d32N>=0, Int) #Partie négative de d21-d31
	@variable(model, 0<=Sign2<=1, Int) #Vaut 1 ssi d21-d31>=0
	@constraint(model,D[1,2]-D[3,2]==d12_d32P-d12_d32N)
	@constraint(model,Sign2>=(D[1,2]-D[3,2])/4)
	@constraint(model,Sign2<=1+(D[1,2]-D[3,2])/4)
	@constraint(model,d12_d32P<=2*Sign2)
	@constraint(model,d12_d32N<=2*(1-Sign2))
	@constraint(model,Kramer[2]-0.5*(D[1,2]+D[3,2])==0.5*(d12_d32P+d12_d32N))
	#Kramer-Simpson Score for player 1
	@variable(model, d13_d23P>=0, Int) #Partie positive de d21-d31
	@variable(model, d13_d23N>=0, Int) #Partie négative de d21-d31
	@variable(model, 0<=Sign3<=1, Int) #Vaut 1 ssi d21-d31>=0
	@constraint(model,D[1,3]-D[2,3]==d13_d23P-d13_d23N)
	@constraint(model,Sign3>=(D[1,3]-D[2,3])/4)
	@constraint(model,Sign3<=1+(D[1,3]-D[2,3])/4)
	@constraint(model,d13_d23P<=2*Sign3)
	@constraint(model,d13_d23N<=2*(1-Sign3))
	@constraint(model,Kramer[3]-0.5*(D[1,3]+D[2,3])==0.5*(d13_d23P+d13_d23N))
	
	
	#Calculate the initial Kramer-Simpson Scores
	D_initial=Array{Int64}(undef,a,a)
	for j in 1:a
		for k in 1:a
			D_initial[j,k]=CP_Maj(Coalitions,m,j,k)[1]
			D_initial[k,j]=CP_Maj(Coalitions,m,j,k)[2]
		end
	end
	Kramer_initial=[KramerSimpsonScore(j,D_initial) for j in 1:a]
	
	#Individual relation
	for j in 1:a
		@constraint(model,rk[i,j]>=(Kramer[j]-Kramer[findall(x->x==i,A)[1]]+1)/n/2)
		@constraint(model,rk[i,j]<=1-(Kramer[findall(x->x==i,A)[1]]-Kramer[j])/n/2)
		@constraint(model,rk[j,i]>=(Kramer[findall(x->x==i,A)[1]]-Kramer[j]+1)/n/2)
		@constraint(model,rk[j,i]<=1-(Kramer[j]-Kramer[findall(x->x==i,A)[1]])/n/2)
	end
	
	#Optimize the individual ranking of i
	@objective(model, Max, sum(rk[i,j]-rk[j,i] for j in 1:a)  )
	
	# Start a chronometer
    start = time()

    # Solve the model
    optimize!(model)
	
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
		if ((A[j]!=i) && (Kramer_initial[i]==Kramer_initial[j])) || ((A[j]!=i) &&  (Kramer_initial[i]>Kramer_initial[j]))
			potential_victims+=1
		end
	end
	
	if feasible
		if (JuMP.objective_value(model)>sum((Kramer_initial[i]<Kramer_initial[j])-(Kramer_initial[j]<Kramer_initial[i]) for j in 1:a))
			manipulable=true
		end
		for j in 1:a
			if ((A[j]!=i) && (JuMP.value(Kramer[i])<JuMP.value(Kramer[j])) && (Kramer_initial[i]==Kramer_initial[j])) || ((A[j]!=i) && (JuMP.value(Kramer[j])>=JuMP.value(Kramer[i])) && (Kramer_initial[i]>Kramer_initial[j]))
				victims+=1
				println("Manipulable for: ",i," against ",A[j])
			end
		end
	NonInitialementPremier=false
	for j in 1:a
		if (Kramer_initial[i]>Kramer_initial[j])
			NonInitialementPremier=true
		end
	end
	Premier=true
	for j in 1:a
		if (JuMP.value(Kramer[i])>JuMP.value(Kramer[j]))
			Premier=false
		end
	end
	DevenuPremier=(NonInitialementPremier)&&(Premier)
	println("Non Initialement Premier:", NonInitialementPremier)
	end
	if !manipulable
		println("Not manipulable for candidate ",i)
	end
	println("********************************************")
	println("********************************************")
	#Initial Kramer Scores
	#Kramer Scores
	for j in 1:a
		println("Kramer initial of ",j,": ",Kramer_initial[j])
	end
	println("********************************************")
	println("********************************************")
	#Kramer Scores
	for j in 1:a
		println("Kramer of ",j,": ",JuMP.value(Kramer[j]))
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

function Manipulate4(PR::Array{Array{Array{Any,1},1},1},i::Int64,A::Array{Int64,1})
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
    @variable(model, Kramer[1:a],Int) #Vector of the Kramer-Simpson Scores
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
			@constraint(model,CP[j,k]>=(D[j,k]-D[k,j]+1)/n/2)
			@constraint(model,CP[j,k]<=1+(D[j,k]-D[k,j])/n/2)
			@constraint(model,CP[k,j]>=(D[k,j]-D[j,k]+1)/n/2)
			@constraint(model,CP[k,j]<=1+(D[k,j]-D[j,k])/n/2)
		end	
	end
	
	#Kramer-Simpson Score for player 1
	@variable(model, d21_d31P>=0, Int) #Partie positive de d21-d31
	@variable(model, d21_d31N>=0, Int) #Partie négative de d21-d31
	@variable(model, 0<=Sign11<=1, Int) #Vaut 1 ssi d21-d31>=0
	@variable(model, M1>=0, Int)      #Max(d21,d31)
	@variable(model, M1_d41N>=0, Int) #Partie négative de M1_d41
	@variable(model, M1_d41P>=0, Int) #Partie positive de M1_d41
	@variable(model, 0<=Sign12<=1, Int) #Vaut 1 ssi Max(d21,d31)-d41>=0
	@constraint(model,D[2,1]-D[3,1]==d21_d31P-d21_d31N)
	@constraint(model,Sign11>=(D[2,1]-D[3,1])/16)
	@constraint(model,Sign11<=1+(D[2,1]-D[3,1])/16)
	@constraint(model,d21_d31P<=16*Sign11)
	@constraint(model,d21_d31N<=16*(1-Sign11))
	@constraint(model,M1-0.5*(D[2,1]+D[3,1])==0.5*(d21_d31P+d21_d31N))
	@constraint(model,M1-D[4,1]==M1_d41P-M1_d41N)
	@constraint(model,Sign12>=(M1-D[4,1])/16)
	@constraint(model,Sign12<=1+(M1-D[4,1])/16)
	@constraint(model,M1_d41P<=16*Sign12)
	@constraint(model,M1_d41N<=16*(1-Sign12))
	@constraint(model,Kramer[1]-0.5*(M1+D[4,1])==0.5*(M1_d41P+M1_d41N))
	
	#Kramer-Simpson Score for player 2
	@variable(model, d12_d32P>=0, Int) #Partie positive de d21-d31
	@variable(model, d12_d32N>=0, Int) #Partie négative de d21-d31
	@variable(model, 0<=Sign21<=1, Int) #Vaut 1 ssi d21-d31>=0
	@variable(model, M2>=0, Int)      #Max(d21,d31)
	@variable(model, M2_d42N>=0, Int) #Partie négative de M1_d41
	@variable(model, M2_d42P>=0, Int) #Partie positive de M1_d41
	@variable(model, 0<=Sign22<=1, Int) #Vaut 1 ssi Max(d21,d31)-d41>=0
	@constraint(model,D[1,2]-D[3,2]==d12_d32P-d12_d32N)
	@constraint(model,Sign21>=(D[1,2]-D[3,2])/16)
	@constraint(model,Sign21<=1+(D[1,2]-D[3,2])/16)
	@constraint(model,d12_d32P<=16*Sign21)
	@constraint(model,d12_d32N<=16*(1-Sign21))
	@constraint(model,M2-0.5*(D[1,2]+D[3,2])==0.5*(d12_d32P+d12_d32N))
	@constraint(model,M2-D[4,2]==M2_d42P-M2_d42N)
	@constraint(model,Sign22>=(M2-D[4,2])/16)
	@constraint(model,Sign22<=1+(M2-D[4,2])/16)
	@constraint(model,M2_d42P<=16*Sign22)
	@constraint(model,M2_d42N<=16*(1-Sign22))
	@constraint(model,Kramer[2]-0.5*(M2+D[4,2])==0.5*(M2_d42P+M2_d42N))
	
	#Kramer-Simpson Score for player 3
	@variable(model, d13_d23P>=0, Int) #Partie positive de d21-d31
	@variable(model, d13_d23N>=0, Int) #Partie négative de d21-d31
	@variable(model, 0<=Sign31<=1, Int) #Vaut 1 ssi d21-d31>=0
	@variable(model, M3>=0, Int)      #Max(d21,d31)
	@variable(model, M3_d43N>=0, Int) #Partie négative de M1_d41
	@variable(model, M3_d43P>=0, Int) #Partie positive de M1_d41
	@variable(model, 0<=Sign32<=1, Int) #Vaut 1 ssi Max(d21,d31)-d41>=0
	@constraint(model,D[1,3]-D[2,3]==d13_d23P-d13_d23N)
	@constraint(model,Sign31>=(D[1,3]-D[2,3])/16)
	@constraint(model,Sign31<=1+(D[1,3]-D[2,3])/16)
	@constraint(model,d13_d23P<=16*Sign31)
	@constraint(model,d13_d23N<=16*(1-Sign31))
	@constraint(model,M3-0.5*(D[1,3]+D[2,3])==0.5*(d13_d23P+d13_d23N))
	@constraint(model,M3-D[4,3]==M3_d43P-M3_d43N)
	@constraint(model,Sign32>=(M3-D[4,3])/16)
	@constraint(model,Sign32<=1+(M3-D[4,3])/16)
	@constraint(model,M3_d43P<=16*Sign32)
	@constraint(model,M3_d43N<=16*(1-Sign32))
	@constraint(model,Kramer[3]-0.5*(M3+D[4,3])==0.5*(M3_d43P+M3_d43N))
	
	#Kramer-Simpson Score for player 4
	@variable(model, d14_d24P>=0, Int) #Partie positive de d21-d31
	@variable(model, d14_d24N>=0, Int) #Partie négative de d21-d31
	@variable(model, 0<=Sign41<=1, Int) #Vaut 1 ssi d21-d31>=0
	@variable(model, M4>=0, Int)      #Max(d21,d31)
	@variable(model, M4_d34N>=0, Int) #Partie négative de M1_d41
	@variable(model, M4_d34P>=0, Int) #Partie positive de M1_d41
	@variable(model, 0<=Sign42<=1, Int) #Vaut 1 ssi Max(d21,d31)-d41>=0
	@constraint(model,D[1,4]-D[2,4]==d14_d24P-d14_d24N)
	@constraint(model,Sign41>=(D[1,4]-D[2,4])/16)
	@constraint(model,Sign41<=1+(D[1,4]-D[2,4])/16)
	@constraint(model,d14_d24P<=16*Sign41)
	@constraint(model,d14_d24N<=16*(1-Sign41))
	@constraint(model,M4-0.5*(D[1,4]+D[2,4])==0.5*(d14_d24P+d14_d24N))
	@constraint(model,M4-D[3,4]==M4_d34P-M4_d34N)
	@constraint(model,Sign42>=(M4-D[3,4])/16)
	@constraint(model,Sign42<=1+(M4-D[3,4])/16)
	@constraint(model,M4_d34P<=16*Sign42)
	@constraint(model,M4_d34N<=16*(1-Sign42))
	@constraint(model,Kramer[4]-0.5*(M4+D[4,1])==0.5*(M4_d34P+M4_d34N))


	
	#Individual position of i is maintained
	#Calculate the initial Copeland Scores
	D_initial=Array{Int64}(undef,a,a)
	for j in 1:a
		for k in 1:a
			D_initial[j,k]=CP_Maj(Coalitions,m,j,k)[1]
			D_initial[k,j]=CP_Maj(Coalitions,m,j,k)[2]
		end
	end
	Kramer_initial=[KramerSimpsonScore(j,D_initial) for j in 1:a]
	
	#Individual relation
	for j in 1:a
		@constraint(model,rk[i,j]>=(Kramer[j]-Kramer[findall(x->x==i,A)[1]]+1)/n/2)
		@constraint(model,rk[i,j]<=1-(Kramer[findall(x->x==i,A)[1]]-Kramer[j])/n/2)
		@constraint(model,rk[j,i]>=(Kramer[findall(x->x==i,A)[1]]-Kramer[j]+1)/n/2)
		@constraint(model,rk[j,i]<=1-(Kramer[j]-Kramer[findall(x->x==i,A)[1]])/n/2)
	end
	
	#Optimize the individual ranking of i
	@objective(model, Max, sum(rk[i,j]-rk[j,i] for j in 1:a)  )
	
	# Start a chronometer
    start = time()

    # Solve the model
    optimize!(model)
	
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
		if ((A[j]!=i) && (Kramer_initial[i]==Kramer_initial[j])) || ((A[j]!=i) &&  (Kramer_initial[i]>Kramer_initial[j]))
			potential_victims+=1
		end
	end
	
	if feasible
	if (JuMP.objective_value(model)>sum((Kramer_initial[i]<Kramer_initial[j])-(Kramer_initial[j]<Kramer_initial[i]) for j in 1:a))
			manipulable=true
		end
		for j in 1:a
			if ((A[j]!=i) && (JuMP.value(Kramer[i])<JuMP.value(Kramer[j])) && (Kramer_initial[i]==Kramer_initial[j])) || ((A[j]!=i) && (JuMP.value(Kramer[j])>=JuMP.value(Kramer[i])) && (Kramer_initial[i]>Kramer_initial[j]))
				victims+=1
				println("Manipulable for: ",i," against ",A[j])
			end
		end
	end
		NonInitialementPremier=false
	for j in 1:a
		if (Kramer_initial[i]>Kramer_initial[j])
			NonInitialementPremier=true
		end
	end
	Premier=true
	for j in 1:a
		if (JuMP.value(Kramer[i])>JuMP.value(Kramer[j]))
			Premier=false
		end
	end
	DevenuPremier=(NonInitialementPremier)&&(Premier)
	if !manipulable
		println("Not manipulable for candidate ",i)
	end
	println("********************************************")
	println("********************************************")
	#Initial Kramer Scores
	for j in 1:a
		println(" Initial Kramer of ",j,": ",Kramer_initial[j])
	end
	println("********************************************")
	println("********************************************")
	#Kramer Scores
	for j in 1:a
		println("Kramer of ",j,": ",JuMP.value(Kramer[j]))
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

function Manipulate5(PR::Array{Array{Array{Any,1},1},1},i::Int64,A::Array{Int64,1})
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
    @variable(model, Kramer[1:a],Int) #Vector of the Kramer-Simpson Scores
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
			@constraint(model,CP[j,k]>=(D[j,k]-D[k,j]+1)/n/2)
			@constraint(model,CP[j,k]<=1+(D[j,k]-D[k,j])/n/2)
			@constraint(model,CP[k,j]>=(D[k,j]-D[j,k]+1)/n/2)
			@constraint(model,CP[k,j]<=1+(D[k,j]-D[j,k])/n/2)
		end	
	end
	
	#Kramer-Simpson Score for player 1
	@variable(model, d21_d31P>=0, Int) #Partie positive de d21-d31
	@variable(model, d21_d31N>=0, Int) #Partie négative de d21-d31
	@variable(model, 0<=Sign11<=1, Int) #Vaut 1 ssi d21-d31>=0
	@variable(model, M1>=0, Int)      #Max(d21,d31)
	@variable(model, M1_d41N>=0, Int) #Partie négative de M1_d41
	@variable(model, M1_d41P>=0, Int) #Partie positive de M1_d41
	@variable(model, 0<=Sign12<=1, Int) #Vaut 1 ssi Max(d21,d31)-d41>=0
	@variable(model, M12>=0, Int)      #Max(Max,d51)
	@variable(model, M12_d51N>=0, Int) #Partie négative de M12_d51
	@variable(model, M12_d51P>=0, Int) #Partie positive de M12_d51
	@variable(model, 0<=Sign13<=1, Int) #Vaut 1 ssi Max(Max,d31)-d51>=0
	@constraint(model,D[2,1]-D[3,1]==d21_d31P-d21_d31N)
	@constraint(model,Sign11>=(D[2,1]-D[3,1])/32)
	@constraint(model,Sign11<=1+(D[2,1]-D[3,1])/32)
	@constraint(model,d21_d31P<=32*Sign11)
	@constraint(model,d21_d31N<=32*(1-Sign11))
	@constraint(model,M1-0.5*(D[2,1]+D[3,1])==0.5*(d21_d31P+d21_d31N))
	@constraint(model,M1-D[4,1]==M1_d41P-M1_d41N)
	@constraint(model,Sign12>=(M1-D[4,1])/32)
	@constraint(model,Sign12<=1+(M1-D[4,1])/32)
	@constraint(model,M1_d41P<=32*Sign12)
	@constraint(model,M1_d41N<=32*(1-Sign12))
	@constraint(model,M12-0.5*(M1+D[4,1])==0.5*(M1_d41P+M1_d41N))
	@constraint(model,M12-D[5,1]==M12_d51P-M12_d51N)
	@constraint(model,Sign13>=(M12-D[5,1])/32)
	@constraint(model,Sign13<=1+(M12-D[5,1])/32)
	@constraint(model,M12_d51P<=32*Sign13)
	@constraint(model,M12_d51N<=32*(1-Sign13))
	@constraint(model,Kramer[1]-0.5*(M12+D[5,1])==0.5*(M12_d51P+M12_d51N))
	
	#Kramer-Simpson Score for player 2
	@variable(model, d12_d32P>=0, Int) #Partie positive de d21-d31
	@variable(model, d12_d32N>=0, Int) #Partie négative de d21-d31
	@variable(model, 0<=Sign21<=1, Int) #Vaut 1 ssi d21-d31>=0
	@variable(model, M2>=0, Int)      #Max(d21,d31)
	@variable(model, M2_d42N>=0, Int) #Partie négative de M1_d41
	@variable(model, M2_d42P>=0, Int) #Partie positive de M1_d41
	@variable(model, 0<=Sign22<=1, Int) #Vaut 1 ssi Max(d21,d31)-d41>=0
	@variable(model, M22>=0, Int)      #Max(Max,d51)
	@variable(model, M22_d52N>=0, Int) #Partie négativbarbarocke de M12_d51
	@variable(model, M22_d52P>=0, Int) #Partie positive de M12_d51
	@variable(model, 0<=Sign23<=1, Int) #Vaut 1 ssi Max(Max,d31)-d51>=0
	@constraint(model,D[1,2]-D[3,2]==d12_d32P-d12_d32N)
	@constraint(model,Sign21>=(D[1,2]-D[3,2])/32)
	@constraint(model,Sign21<=1+(D[1,2]-D[3,2])/32)
	@constraint(model,d12_d32P<=32*Sign21)
	@constraint(model,d12_d32N<=32*(1-Sign21))
	@constraint(model,M2-0.5*(D[1,2]+D[3,2])==0.5*(d12_d32P+d12_d32N))
	@constraint(model,M2-D[4,2]==M2_d42P-M2_d42N)
	@constraint(model,Sign22>=(M2-D[4,2])/32)
	@constraint(model,Sign22<=1+(M2-D[4,2])/32)
	@constraint(model,M2_d42P<=32*Sign22)
	@constraint(model,M2_d42N<=32*(1-Sign22))
	@constraint(model,M22-0.5*(M2+D[4,2])==0.5*(M2_d42P+M2_d42N))
	
	@constraint(model,M22-D[5,2]==M22_d52P-M22_d52N)
	@constraint(model,Sign23>=(M22-D[5,2])/32)
	@constraint(model,Sign23<=1+(M22-D[5,2])/32)
	@constraint(model,M22_d52P<=32*Sign23)
	@constraint(model,M22_d52N<=32*(1-Sign23))
	@constraint(model,Kramer[2]-0.5*(M22+D[5,2])==0.5*(M22_d52P+M22_d52N))
	
	#Kramer-Simpson Score for player 3
	@variable(model, d13_d23P>=0, Int) #Partie positive de d21-d31
	@variable(model, d13_d23N>=0, Int) #Partie négative de d21-d31
	@variable(model, 0<=Sign31<=1, Int) #Vaut 1 ssi d21-d31>=0
	@variable(model, M3>=0, Int)      #Max(d21,d31)
	@variable(model, M3_d43N>=0, Int) #Partie négative de M1_d41
	@variable(model, M3_d43P>=0, Int) #Partie positive de M1_d41
	@variable(model, 0<=Sign32<=1, Int) #Vaut 1 ssi Max(d21,d31)-d41>=0
	@variable(model, M32>=0, Int)      #Max(Max,d51)
	@variable(model, M32_d53N>=0, Int) #Partie négativbarbarocke de M12_d51
	@variable(model, M32_d53P>=0, Int) #Partie positive de M12_d51
	@variable(model, 0<=Sign33<=1, Int) #Vaut 1 ssi Max(Max,d31)-d51>=0
	@constraint(model,D[1,3]-D[2,3]==d13_d23P-d13_d23N)
	@constraint(model,Sign31>=(D[1,3]-D[2,3])/32)
	@constraint(model,Sign31<=1+(D[1,3]-D[2,3])/32)
	@constraint(model,d13_d23P<=16*Sign31)
	@constraint(model,d13_d23N<=16*(1-Sign31))
	@constraint(model,M3-0.5*(D[1,3]+D[2,3])==0.5*(d13_d23P+d13_d23N))
	@constraint(model,M3-D[4,3]==M3_d43P-M3_d43N)
	@constraint(model,Sign32>=(M3-D[4,3])/32)
	@constraint(model,Sign32<=1+(M3-D[4,3])/32)
	@constraint(model,M3_d43P<=32*Sign32)
	@constraint(model,M3_d43N<=32*(1-Sign32))
	@constraint(model,M32-0.5*(M3+D[4,3])==0.5*(M3_d43P+M3_d43N))

	@constraint(model,M32-D[5,3]==M32_d53P-M32_d53N)
	@constraint(model,Sign33>=(M32-D[5,3])/32)
	@constraint(model,Sign33<=1+(M32-D[5,3])/32)
	@constraint(model,M32_d53P<=32*Sign33)
	@constraint(model,M32_d53N<=32*(1-Sign33))
	@constraint(model,Kramer[3]-0.5*(M32+D[5,3])==0.5*(M32_d53P+M32_d53N))
	
	#Kramer-Simpson Score for player 4
	@variable(model, d14_d24P>=0, Int) #Partie positive de d21-d31
	@variable(model, d14_d24N>=0, Int) #Partie négative de d21-d31
	@variable(model, 0<=Sign41<=1, Int) #Vaut 1 ssi d21-d31>=0
	@variable(model, M4>=0, Int)      #Max(d21,d31)
	@variable(model, M4_d34N>=0, Int) #Partie négative de M1_d41
	@variable(model, M4_d34P>=0, Int) #Partie positive de M1_d41
	@variable(model, 0<=Sign42<=1, Int) #Vaut 1 ssi Max(d21,d31)-d41>=0
	@variable(model, M42>=0, Int)      #Max(Max,d51)
	@variable(model, M42_d54N>=0, Int) #Partie négativbarbarocke de M12_d51
	@variable(model, M42_d54P>=0, Int) #Partie positive de M12_d51
	@variable(model, 0<=Sign43<=1, Int) #Vaut 1 ssi Max(Max,d31)-d51>=0
	@constraint(model,D[1,4]-D[2,4]==d14_d24P-d14_d24N)
	@constraint(model,Sign41>=(D[1,4]-D[2,4])/32)
	@constraint(model,Sign41<=1+(D[1,4]-D[2,4])/32)
	@constraint(model,d14_d24P<=2*Sign41)
	@constraint(model,d14_d24N<=2*(1-Sign41))
	@constraint(model,M4-0.5*(D[1,4]+D[2,4])==0.5*(d14_d24P+d14_d24N))
	@constraint(model,M4-D[3,4]==M4_d34P-M4_d34N)
	@constraint(model,Sign42>=(M4-D[3,4])/32)
	@constraint(model,Sign42<=1+(M4-D[3,4])/32)
	@constraint(model,M4_d34P<=32*Sign42)
	@constraint(model,M4_d34N<=32*(1-Sign42))
	@constraint(model,M42-0.5*(M4+D[4,1])==0.5*(M4_d34P+M4_d34N))
	
	@constraint(model,M42-D[5,4]==M42_d54P-M42_d54N)
	@constraint(model,Sign43>=(M42-D[5,4])/32)
	@constraint(model,Sign43<=1+(M42-D[5,4])/32)
	@constraint(model,M42_d54P<=32*Sign33)
	@constraint(model,M42_d54N<=32*(1-Sign33))
	@constraint(model,Kramer[4]-0.5*(M42+D[5,4])==0.5*(M42_d54P+M42_d54N))
	
	#Kramer-Simpson Score for player 5
	@variable(model, d15_d25P>=0, Int) #Partie positive de d21-d31
	@variable(model, d15_d25N>=0, Int) #Partie négative de d21-d31
	@variable(model, 0<=Sign51<=1, Int) #Vaut 1 ssi d21-d31>=0
	@variable(model, M5>=0, Int)      #Max(d21,d31)
	@variable(model, M5_d35N>=0, Int) #Partie négative de M1_d41
	@variable(model, M5_d35P>=0, Int) #Partie positive de M1_d41
	@variable(model, 0<=Sign52<=1, Int) #Vaut 1 ssi Max(d21,d31)-d41>=0
	@variable(model, M52>=0, Int)      #Max(Max,d51)
	@variable(model, M52_d45N>=0, Int) #Partie négativbarbarocke de M12_d51
	@variable(model, M52_d45P>=0, Int) #Partie positive de M12_d51
	@variable(model, 0<=Sign53<=1, Int) #Vaut 1 ssi Max(Max,d31)-d51>=0
	@constraint(model,D[1,5]-D[2,5]==d15_d25P-d15_d25N)
	@constraint(model,Sign51>=(D[1,5]-D[2,5])/32)
	@constraint(model,Sign51<=1+(D[1,5]-D[2,5])/32)
	@constraint(model,d15_d25P<=2*Sign51)
	@constraint(model,d15_d25N<=2*(1-Sign51))
	@constraint(model,M5-0.5*(D[1,5]+D[2,5])==0.5*(d15_d25P+d15_d25N))
	@constraint(model,M5-D[3,5]==M5_d35P-M5_d35N)
	@constraint(model,Sign52>=(M5-D[3,5])/32)
	@constraint(model,Sign52<=1+(M5-D[3,5])/32)
	@constraint(model,M5_d35P<=32*Sign52)
	@constraint(model,M5_d35N<=32*(1-Sign52))
	@constraint(model,M52-0.5*(M5+D[3,5])==0.5*(M5_d35P+M5_d35N))
	
	@constraint(model,M52-D[4,5]==M52_d45P-M52_d45N)
	@constraint(model,Sign53>=(M52-D[4,5])/32)
	@constraint(model,Sign53<=1+(M52-D[4,5])/32)
	@constraint(model,M52_d45P<=32*Sign53)
	@constraint(model,M52_d45N<=32*(1-Sign53))
	@constraint(model,Kramer[5]-0.5*(M52+D[4,5])==0.5*(M52_d45P+M52_d45N))


	
	#Individual position of i is maintained
	#Calculate the initial Copeland Scores
	D_initial=Array{Int64}(undef,a,a)
	for j in 1:a
		for k in 1:a
			D_initial[j,k]=CP_Maj(Coalitions,m,j,k)[1]
			D_initial[k,j]=CP_Maj(Coalitions,m,j,k)[2]
		end
	end
	Kramer_initial=[KramerSimpsonScore(j,D_initial) for j in 1:a]

	#Individual relation
	for j in 1:a
		@constraint(model,rk[i,j]>=(Kramer[j]-Kramer[findall(x->x==i,A)[1]]+1)/n/2)
		@constraint(model,rk[i,j]<=1-(Kramer[findall(x->x==i,A)[1]]-Kramer[j])/n/2)
		@constraint(model,rk[j,i]>=(Kramer[findall(x->x==i,A)[1]]-Kramer[j]+1)/n/2)
		@constraint(model,rk[j,i]<=1-(Kramer[j]-Kramer[findall(x->x==i,A)[1]])/n/2)
	end
	
	#Optimize the individual ranking of i
	@objective(model, Max, sum(rk[i,j]-rk[j,i] for j in 1:a)  )
	
	# Start a chronometer
    start = time()

    # Solve the model
    optimize!(model)
	
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
		if ((A[j]!=i) && (Kramer_initial[i]==Kramer_initial[j])) || ((A[j]!=i) &&  (Kramer_initial[i]>Kramer_initial[j]))
			potential_victims+=1
		end
	end
	
	if feasible
		if (JuMP.objective_value(model)>sum((Kramer_initial[i]<Kramer_initial[j])-(Kramer_initial[j]<Kramer_initial[i]) for j in 1:a))
			manipulable=true
		end
		if (JuMP.objective_value(model)>sum((Kramer_initial[i]>Kramer_initial[j])-(Kramer_initial[j]>Kramer_initial[i]) for j in 1:a))
		end
		for j in 1:a
			if ((A[j]!=i) && (JuMP.value(Kramer[i])<JuMP.value(Kramer[j])) && (Kramer_initial[i]==Kramer_initial[j])) || ((A[j]!=i) && (JuMP.value(Kramer[j])>=JuMP.value(Kramer[i])) && (Kramer_initial[i]>Kramer_initial[j]))
				victims+=1
				println("Manipulable for: ",i," against ",A[j])
			end
		end
	end
		NonInitialementPremier=false
	for j in 1:a
		if (Kramer_initial[i]>Kramer_initial[j])
			NonInitialementPremier=true
		end
	end
	Premier=true
	for j in 1:a
		if (JuMP.value(Kramer[i])>JuMP.value(Kramer[j]))
			Premier=false
		end
	end
	DevenuPremier=(NonInitialementPremier)&&(Premier)
	if !manipulable
		println("Not manipulable for candidate ",i)
	end
	println("********************************************")
	println("********************************************")
	#Initial Kramer Scores
	for j in 1:a
		println(" Initial Kramer of ",j,": ",Kramer_initial[j])
	end
	println("********************************************")
	println("********************************************")
	#Kramer Scores
	for j in 1:a
		println("Kramer of ",j,": ",JuMP.value(Kramer[j]))
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

