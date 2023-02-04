include("manip.jl")

function SolveDataSet()
	global N=3;
	global Total=0;
	global Manip=0;
	global Total_manipulators=0;
	global Victims_Prop=0;
	global SommeCarré=0;
	global NonInitialementPremiers=0;
	global DevenuPremiers=0;
	dataFolder = "../data/"
	L=[1:N;]
	for i in 1:N
		L[i]=0
	end
	for file in filter(x->occursin(".txt", x), readdir(dataFolder)) 
		Manipulators=0
		ManipParticuliers=0 #Nombre manipulateurs de cette PR
		println("-- Resolution of ", file)
    	PR=ReadInputFile(dataFolder * file)
    	Total+=1
    	manipulable=false
    	for i in 1:N
		test,Victim=Manipulate3(PR,i,[1:N;])[1:2] ### Choose the right function: Manipulate3, Manipulate4 ...
    		if (test)
    			manipulable=true
				Victims_Prop+=Victim
    			Total_manipulators+=1
    			Manipulators+=1
    			ManipParticuliers+=1 #Nombre manipulateurs de cette PR
    			DevenuPremiers+=Manipulate3(PR,i,[1:N;])[4]
				NonInitialementPremiers+=Manipulate3(PR,i,[1:N;])[3] ### Choose the right function: Manipulate3, Manipulate4 ...
    		end
    	end
    	SommeCarré+=ManipParticuliers*ManipParticuliers
    	if (Manipulators>0)
    		L[Manipulators]+=1
    	end
    	Manip+=manipulable
	end
	println("Le nombre total des Power Relations: ,", Total)
	println("Le nombre total des Power Relations manipulables: ", Manip)
	println("La proportion des Power Relation manipulables est: ", Manip/Total*100,"%")
	println("La moyenne empirique des nombre de manipulateurs est: ", Total_manipulators/Manip)
	println("Probabilité de devenir premier: ", DevenuPremiers/NonInitialementPremiers)
	println("non init premiers: ", NonInitialementPremiers)
	println("devenu premiers: ",DevenuPremiers)
end

