module CRNHexagon

import HomotopyContinuation: @var, evaluate, Expression
import LinearAlgebra: inv, det
import ProgressMeter: @showprogress

include("plotting_functionality.jl")

export runTest, computeCoverInvariants, printValues, runTest_weighted

#=
This is the main method. Use it to run all tests. By generating the circuit numbers of all 16 circuit covers
associated to the hexagonal face H of the Newton polytope of the dual phosphorylation network, we are able to 
compare their quality. All covers are compared to the cover CC(9) introduced by Feliu et al.
=#
function runTest( ; boxsize=1, numberOfSamplingRuns=150, prefix="michaelismentontest", suffix="NEW")
    @var κ[1:12]
    hexPoints = [(0,0),(1,0),(2,0),(4,1),(4,2),(3,2),(2,2),(0,1),(3,1),(1,1)]
    K = [(κ[2]+κ[3])/κ[1], (κ[5]+κ[6])/κ[4], (κ[8]+κ[9])/κ[7], (κ[11]+κ[12])/κ[10]]
    aη = κ[3]*κ[12] - κ[6]*κ[9]
    bη = (K[2] + K[3])*κ[3]*κ[12] - (K[1]+K[4])*κ[6]*κ[9]
    coefficients = [K[1]^3*K[3]^2*κ[6]^3*κ[12]^2, K[1]^2*K[2]*K[3]^2*κ[3]*κ[6]^2*κ[12]^2, K[1]^2*K[2]*K[3]*K[4]*κ[3]*κ[6]^2*κ[9]*κ[12], K[1]*K[2]^2*K[4]*κ[3]^2*κ[6]*κ[9]^2,
                    K[2]^2*K[4]*κ[3]^2*κ[9]*aη, K[2]^2*K[3]*κ[3]^2*κ[12]*aη, K[1]*K[2]*K[3]*κ[3]*κ[6]*κ[12]*aη, K[1]^2*K[3]^2*κ[6]^3*κ[12]^2, 
                    2*K[1]*K[2]*K[3]*K[4]*κ[3]^2*κ[6]*κ[9]*κ[12], 2*K[1]^2*K[2]*K[3]*κ[3]*κ[6]^2*κ[12]^2]
    mcoef = K[1]*K[2]*K[3]*κ[3]*κ[6]*κ[12]*bη

    triangconfigurations = [[[2,7,9],[3,6,10],[1,5],[4,8]], [[3,6,10],[2,4,7],[1,5],[8,9]], [[1,3,6],[2,5,7],[9,10],[4,8]],
    [[1,3,6],[2,5,7],[8,9],[4,10]], [[3,6,8],[2,7,9],[4,10],[1,5]], [[2,4,7],[3,6,8],[1,5],[9,10]],
    [[1,4,6],[2,5,8],[3,7],[9,10]], [[1,7,9],[3,5,10],[2,6],[4,8]], [[3,5,8],[1,4,7],[9,10],[2,6]],
    [[1,7,9],[3,5,8],[2,6],[4,10]], [[1,6,9],[2,5,8],[3,7],[4,10]], [[1,4,7],[3,5,10],[8,9],[2,6]],
    [[2,5,10],[1,4,6],[3,7],[8,9]], [[1,6,9],[2,5,10],[3,7],[4,8]]]
    #Check if all covers are legitimate
    for config in triangconfigurations
        all(t->t in union(config[1],config[2],config[3],config[4]),1:10) || display(config)&&throw(error("The triangles don't cover the entire region"))
        isempty(intersect(config[1],config[2])) && isempty(intersect(config[1],config[3])) && isempty(intersect(config[1],config[4])) && isempty(intersect(config[2],config[3])) && isempty(intersect(config[2],config[4])) && isempty(intersect(config[3],config[4])) || display(config)&&throw(error("Each vertex should only be used once"))
    end

    lineconfigurations = [[[5,1],[7,3],[8,9],[6,2],[4,10]], [[5,1],[7,3],[9,10],[6,2],[8,4]]]

    θ = createθcircuits(hexPoints, coefficients, lineconfigurations, triangconfigurations)
    runSamplingComparison(θ, κ, aη, bη, mcoef, θ[9]; boxsize=boxsize, numberOfSamplingRuns=numberOfSamplingRuns, prefix=prefix, suffix=suffix)
end

#=
This function produces the Latex output for the sampled data from the method `runTest`. This data leads to
Tables 1 and 2 in the article connected to this repository.
=#
function computeCoverInvariants( ; boxsizes=["0.1","1","10","100"], prefix="michaelismentontest", suffix="NEW")
    for boxsize in boxsizes
        try
            f = open("../data/$(prefix)$(suffix)triangstoredsolutions$(boxsize).txt", "r")
            
            global pointnumber = parse(Int,readline(f))
            global allmodels = parse(Int,readline(f))
            global onlyone = [parse(Int,entry) for entry in split(readline(f)[2:end-1], ", ")]
            global allpoints = [parse(Int,entry) for entry in split(readline(f)[2:end-1], ", ")]
            global ourmodel = [parse(Int,entry) for entry in split(readline(f)[2:end-1], ", ")]
            global prevmodel = [parse(Int,entry) for entry in split(readline(f)[2:end-1], ", ")]
            global nomodel = [parse(Int,entry) for entry in split(readline(f)[2:end-1], ", ")]
            display(allmodels/pointnumber)
            close(f)
        catch
            continue
        end
        println("$(boxsize): $(pointnumber)")
        
        permille_allmodels, permille_ourmodel, permille_prevmodel, permille_nomodel, percent_ourmodel_pure = round.(round.(10000*allpoints ./ pointnumber, digits=2)/10000, digits=5), round.(100*ourmodel ./ pointnumber, digits=2), round.(100*prevmodel ./ pointnumber, digits=2), round.(100*nomodel ./ pointnumber, digits=2), round.(100 * (pointnumber .- (nomodel .+ prevmodel)) ./ pointnumber, digits=3)
        foreach(i->push!(allmodeldots[i], permille_allmodels[i]), 1:length(permille_allmodels))
        foreach(i->push!(ourmodeldots[i], permille_ourmodel[i]), 1:length(permille_ourmodel))
        foreach(i->push!(prevmodeldots[i], permille_prevmodel[i]), 1:length(permille_ourmodel))
        foreach(i->push!(nomodeldots[i], permille_nomodel[i]), 1:length(permille_ourmodel))
        foreach(i->push!(puremodel[i], percent_ourmodel_pure[i]), 1:length(percent_ourmodel_pure))
    end

    
    print("~&+&-&0&+&-&0&+&-&0\\\\ \\thickhline \n\n")
    for θ in 1:Int(length(ourmodeldots))
        print("$(θ)")
        for i in 1:length(ourmodeldots[θ])
            print("&$(ourmodeldots[θ][i])&$(prevmodeldots[θ][i])&$(nomodeldots[θ][i])")
        end
        print("\\\\ \\hline \n\n")
    end


    for θ in 1:Int(length(allmodeldots))
        print("$(θ)&")
        for i in 1:length(allmodeldots[θ])-1
            print("$(allmodeldots[θ][i])&")
        end
        print("$(allmodeldots[θ][end])")
        print("\\\\ \\hline \n")
    end
end


#=
The method `printValues` displays for the data previously obtained from `runTest`.
=#
function printValues( ; boxsizes = ["0.1","1","10","100"], prefix="michaelismentontest", suffix="NEW")
    for boxsize in boxsizes
        print("$(boxsize): \n")
        f = open("../data/$(prefix)$(suffix)triangstoredsolutions$(boxsize).txt", "r")
        global relDict = Dict()
        global pointnumber = parse(Int,readline(f))
        global allmodels = parse(Int,readline(f))
        global onlyonemodel = [parse(Int,entry) for entry in split(readline(f)[2:end-1], ", ")]
        global newmodel = [parse(Int,entry) for entry in split(readline(f)[2:end-1], ", ")]
        global ourmodel = [parse(Int,entry) for entry in split(readline(f)[2:end-1], ", ")]
        global prevmodel = [parse(Int,entry) for entry in split(readline(f)[2:end-1], ", ")]
        global nomodel = [parse(Int,entry) for entry in split(readline(f)[2:end-1], ", ")]

        while ! eof(f)  
            sstring = split(readline(f), ": ")
            keystring = split(sstring[1], ", ")
            key = (parse(Int,keystring[1]), (parse(Int,keystring[2])))
            relDict[key] = parse(Int,sstring[2])
        end

        close(f)
        containmentArray, almostContainmentArray = [[] for _ in 1:16], [[] for _ in 1:16]
        for key in keys(relDict)
            if relDict[key]==0
                push!(containmentArray[key[1]], key[2])
            elseif relDict[key]<10
                push!(almostContainmentArray[key[1]], key[2])
            end
        end
        println("Containment")
        foreach(t->println("$(t): $(containmentArray[t])"), 1:16)
        println("AlmostContainment")
        foreach(t->println("$(t): $(almostContainmentArray[t])"), 1:16)
    end
end


function runTest_weighted(; boxsize=1, numberOfSamplingRuns=100, prefix="linearweight", suffix="4,9", discretization=10)
    @var κ[1:12]

    #We choose colors with maximum distinguishability
    hexPoints = [(0,0),(1,0),(2,0),(4,1),(4,2),(3,2),(2,2),(0,1),(3,1),(1,1)]
    K = [(κ[2]+κ[3])/κ[1], (κ[5]+κ[6])/κ[4], (κ[8]+κ[9])/κ[7], (κ[11]+κ[12])/κ[10]]
    aη = κ[3]*κ[12] - κ[6]*κ[9]
    bη = (K[2] + K[3])*κ[3]*κ[12] - (K[1]+K[4])*κ[6]*κ[9]
    coefficients = [K[1]^3*K[3]^2*κ[6]^3*κ[12]^2, K[1]^2*K[2]*K[3]^2*κ[3]*κ[6]^2*κ[12]^2, K[1]^2*K[2]*K[3]*K[4]*κ[3]*κ[6]^2*κ[9]*κ[12], K[1]*K[2]^2*K[4]*κ[3]^2*κ[6]*κ[9]^2,
                    K[2]^2*K[4]*κ[3]^2*κ[9]*aη, K[2]^2*K[3]*κ[3]^2*κ[12]*aη, K[1]*K[2]*K[3]*κ[3]*κ[6]*κ[12]*aη, K[1]^2*K[3]^2*κ[6]^3*κ[12]^2, 
                    2*K[1]*K[2]*K[3]*K[4]*κ[3]^2*κ[6]*κ[9]*κ[12], 2*K[1]^2*K[2]*K[3]*κ[3]*κ[6]^2*κ[12]^2]
    mcoef = K[1]*K[2]*K[3]*κ[3]*κ[6]*κ[12]*bη
    configurations = [[[1,7,9],[3,5,8],[2,6],[4,10]], [[1,4,7],[3,5,10],[8,9],[2,6]], [[1,5],[3,7],[8,9],[2,6],[4,10]]]

    triangconfigurations = [[[2,7,9],[3,6,10],[1,5],[4,8]], [[3,6,10],[2,4,7],[1,5],[8,9]], [[1,3,6],[2,5,7],[9,10],[4,8]],
    [[1,3,6],[2,5,7],[8,9],[4,10]], [[3,6,8],[2,7,9],[4,10],[1,5]], [[2,4,7],[3,6,8],[1,5],[9,10]],
    [[1,4,6],[2,5,8],[3,7],[9,10]], [[1,7,9],[3,5,10],[2,6],[4,8]], [[3,5,8],[1,4,7],[9,10],[2,6]],
    [[1,7,9],[3,5,8],[2,6],[4,10]], [[1,6,9],[2,5,8],[3,7],[4,10]], [[1,4,7],[3,5,10],[8,9],[2,6]],
    [[2,5,10],[1,4,6],[3,7],[8,9]], [[1,6,9],[2,5,10],[3,7],[4,8]]]
    #Check if all covers are legit
    for config in triangconfigurations
        all(t->t in union(config[1],config[2],config[3],config[4]),1:10) || display(config)&&throw(error("The triangles don't cover the entire region"))
        isempty(intersect(config[1],config[2])) && isempty(intersect(config[1],config[3])) && isempty(intersect(config[1],config[4])) && isempty(intersect(config[2],config[3])) && isempty(intersect(config[2],config[4])) && isempty(intersect(config[3],config[4])) || display(config)&&throw(error("Each vertex should only be used once"))
    end

    lineconfigurations = [[[5,1],[7,3],[8,9],[6,2],[4,10]], [[5,1],[7,3],[9,10],[6,2],[8,4]]]

    all(t->sort(vcat(t...))==[i for i in 1:10], configurations)||throw(error("Not sorted correctly!"))
    oldθ = createθcircuits(hexPoints, coefficients, [], [[[1,3,6],[2,5,7],[8,9],[4,10]], [[3,5,8],[1,4,7],[9,10],[2,6]]])
    #θ = createθcircuits(hexPoints, coefficients, lineconfigurations, triangconfigurations)
    θ_weighted = createθcircuits_weighted(hexPoints, coefficients, [[[1,3,6],[2,5,7],[8,9],[4,10]], [[3,5,8],[1,4,7],[9,10],[2,6]]]; discretization=discretization)
    #plotAllCovers(hexPoints, triangconfigurations, lineconfigurations)
    runSamplingComparison_weighted([], θ_weighted, κ, aη, bη, mcoef, oldθ; boxsize=boxsize, numberOfSamplingRuns=numberOfSamplingRuns, prefix=prefix, suffix=suffix, discretization=discretization)    
end

end 
