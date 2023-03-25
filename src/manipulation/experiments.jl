include("banzhaf.jl")
include("copeland.jl")
include("kramer_simpson.jl")
include("generation.jl")
using JLD


function SolveDataSet(Agents::Int64,Instances::Int64,Method::String)
	TotalManipulablePR=0;
	TotalManipulators=0;

	if Method == "Banzhaf"
		Manipulate = ManipulateBanzhaf
	else
		Manipulate = ManipulateCopeland
	end
	
	L=[1:Agents;]
	for i in 1:Agents
		L[i]=0
	end

	for instance in 1:Instances 
		Manipulators=0
		PR = TotalPowerRelation(Agents)
    	manipulable=false
    	for i in 1:Agents
			manipulable = Manipulate(PR,i,[1:Agents;])[1]
			Manipulators+=manipulable
			TotalManipulators+=manipulable
    	end
    	if (manipulable)
    		L[Manipulators]+=1
    	end
    	TotalManipulablePR+=manipulable
	end
	#DISPLAY
	println("Number of instances", Instaces)
	println("Number of manipulable Power Relations: ", TotalManipulablePR)
	println("Percentage of manipulable Power Relations", TotalManipulablePR/Instances*100,"%")
	println("Average number of manipulators per Power Relation: ", TotalManipulators/Manip)
	for i in 1:Agents
		println("Percentage of Power Relations with ",i," manipulators: ", L[i]/TotalManipulablePR*100,"%")
	end
end

