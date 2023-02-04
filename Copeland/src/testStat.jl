include("manip.jl")

function SolveDataSet()
	#Number of individuals
	global N=5;
	#Initialization
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
		test,Victim=Manipulate(PR,i,[1:N;])[1],Manipulate(PR,i,[1:N;])[2]
    		if (test)
    			manipulable=true
				Victims_Prop+=Victim
    			Total_manipulators+=1
    			Manipulators+=1
    			ManipParticuliers+=1
    			DevenuPremiers+=Manipulate(PR,i,[1:N;])[4]
				NonInitialementPremiers+=Manipulate(PR,i,[1:N;])[3]
    		end
    	end
    	SommeCarré+=ManipParticuliers*ManipParticuliers
    	if (Manipulators>0)
    		L[Manipulators]+=1
    	end
    	Manip+=manipulable
	end
	#DISPLAY
	println("Le nombre total des Power Relations: ,", Total)
	println("Le nombre total des Power Relations manipulables: ", Manip)
	println("La proportion des Power Relation manipulables est: ", Manip/Total*100,"%")
	println("La moyenne empirique des nombre de manipulateurs est: ", Total_manipulators/Manip)
println("L'écart-type des nombre de manipulateurs est: ", sqrt((SommeCarré/Manip)-(Total_manipulators/Manip)*(Total_manipulators/Manip)))
	for i in 1:N
		println("Pourcentage ",i," manipulateur: ", L[i]/Manip*100,"%")
	end
	println("Moyenne du rapport Victims/Potential victims est: ", Victims_Prop/Total_manipulators*100)
	println("Probabilité de devenir premier: ", DevenuPremiers/NonInitialementPremiers)
end

