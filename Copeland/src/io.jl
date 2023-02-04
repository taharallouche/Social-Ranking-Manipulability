using JuMP
using Plots
import GR

"""
Saves Power Relation into a txt file
"""
function SaveInstance(PR::Array{Array{Array{Any,1},1},1},outputFile::String)
	# Open the output file
    writer = open(outputFile, "w")
    
    # Each coalition is a String
    # Each indifference class is a line
    C=size(PR,1) #Number of indifference classes
    for i in 1:C
    	for j in 1:size(PR[i],1)
    		ch=join(PR[i][j])
    		print(writer,ch)
    		if j != size(PR[i],1)
                print(writer, ",")
            else
                println(writer, "")
            end
    	end
    end
    close(writer)
end


"""
Reads Power Relation from txt file
"""
function ReadInputFile(inputFile::String)
	# Open the input file
    datafile = open(inputFile)

    data = readlines(datafile)
    close(datafile)
    #First class
    line=data[1]
	lineSplit = split(line, ",")
    Class=Array{Any,1}[]
    for coalition in lineSplit
      	if length(coalition)==0
      		Class=vcat(Class,[[]])
      	else
  			Coalition=split(coalition,"")
  			CoalitionInt=Any[parse(Int64, Coalition[1])]
      		if (size(Coalition,1)>1)
				for j in 2:size(Coalition,1)
					CoalitionInt=vcat(CoalitionInt,parse(Int64, Coalition[j]))
				end
			end
			Class=vcat(Class,[CoalitionInt])
  		end
  	end   
  	PR=Array{Array{Any,1},1}[Class] 
    # For each line of the input file
    if (size(data,1)>1)
    	for line in data[2:size(data,1)]
        	lineSplit = split(line, ",")
        	Class=Array{Any,1}[]
      		for coalition in lineSplit
      			if length(coalition)==0
      				Class=vcat(Class,[[]])
      			else
      				Coalition=split(coalition,"")
      				CoalitionInt=Any[parse(Int64, Coalition[1])]
      				if (size(Coalition,1)>1)
      					for j in 2:size(Coalition,1)
      						CoalitionInt=vcat(CoalitionInt,parse(Int64, Coalition[j]))
      					end
      				end
      				Class=vcat(Class,[CoalitionInt])
      			end
      		end
      		PR=vcat(PR,[Class])
      	end
    end
    return PR
end
