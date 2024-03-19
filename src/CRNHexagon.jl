module CRNHexagon

import GLMakie: rotate!, text!, Colorbar, heatmap, heatmap!, xlims!, ylims!, plot!, Point2f0, lines!, Figure, Axis, save, hidespines!, hidedecorations!, mesh!, scatter!, text!, RGBA, RGB, poly!, Axis3, Point3f0
import HomotopyContinuation: @var, evaluate, Expression
import LinearAlgebra: inv, det
import ProgressBars: ProgressBar
import Colors: distinguishable_colors, red, green, blue, colormap
import LaTeXStrings: @L_str
import Polyhedra: Mesh, polyhedron, convexhull

export runTest, computeCoverInvariants, compareTwoCovers, runTest_noDependencies

#=
Here, the correct θ's are calculated for all possible circuits. 
TODO: Dissect the method to allow for weightings.
=#
function createθcircuits(points, coefficients, lineconfigurations, triangconfigurations)
    m = (2,1)
    λ = []
    θ = []

    for config in triangconfigurations
        try
            global barycenter = Matrix{Float64}(undef,3,3); barycenter[1,:] = [points[entry][1] for entry in config[1]]; barycenter[2,:] = [points[entry][2] for entry in config[1]]; barycenter[3,:] = [1 for entry in config[1]];
            global msolve = [2,1,1];
            global λ1 = inv(barycenter)*msolve;
            global λ1 = collect(λ1 / sum(λ1));

            barycenter[1,:] = [points[entry][1] for entry in config[2]]; barycenter[2,:] = [points[entry][2] for entry in config[2]];
            global λ2 = inv(barycenter)*msolve;
            global λ2 = collect(λ2 / sum(λ2));

            global barycenterline = Matrix{Float64}(undef,2,2); barycenterline[1,:] = [points[entry][1] for entry in config[3]]; barycenterline[2,:] = [points[entry][2] for entry in config[3]]; 
            if det(barycenterline) == 0
                global λ3 = [0.5,0.5]
            else
                global λ3 = inv(barycenterline)*[2,1];
                global λ3 = collect(λ3 / sum(λ3))
            end

            barycenterline[1,:] = [points[entry][1] for entry in config[4]]; barycenterline[2,:] = [points[entry][2] for entry in config[4]]; 
            if det(barycenterline) == 0
                global λ4 = [0.5,0.5]
            else
                global λ4 = inv(barycenterline)*[2,1];
                global λ4 = collect(λ4 / sum(λ4))
            end

            global θ1 = prod([(coefficients[config[1][i]]/λ1[i])^(λ1[i]) for i in 1:length(config[1])])
            global θ2 = prod([(coefficients[config[2][i]]/λ2[i])^(λ2[i]) for i in 1:length(config[2])])
            global θ3 = prod([(coefficients[config[3][i]]/λ3[i])^(λ3[i]) for i in 1:length(config[3])])
            global θ4 = prod([(coefficients[config[4][i]]/λ4[i])^(λ4[i]) for i in 1:length(config[4])])
            push!(θ,θ1+θ2+θ3+θ4)
        catch
            continue
        end
    end

    for config in lineconfigurations
        try
            global barycenterline = Matrix{Float64}(undef,2,2); barycenterline[1,:] = [points[entry][1] for entry in config[1]]; barycenterline[2,:] = [points[entry][2] for entry in config[1]]; 
            if det(barycenterline) == 0
                global λ1=[0.5,0.5]
            else
                global λ1 = inv(barycenterline)*[2,1];
                global λ1 = collect(λ1 / sum(λ1))
            end

            barycenterline[1,:] = [points[entry][1] for entry in config[2]]; barycenterline[2,:] = [points[entry][2] for entry in config[2]]; 
            if det(barycenterline) == 0
                global λ2=[0.5,0.5]
            else
                global λ2 = inv(barycenterline)*[2,1];
                global λ2 = collect(λ2 / sum(λ2))
            end

            barycenterline[1,:] = [points[entry][1] for entry in config[3]]; barycenterline[2,:] = [points[entry][2] for entry in config[3]]; 
            if det(barycenterline) == 0
                global λ3=[0.5,0.5]
            else
                global λ3 = inv(barycenterline)*[2,1];
                global λ3 = collect(λ3 / sum(λ3))
            end

            barycenterline[1,:] = [points[entry][1] for entry in config[4]]; barycenterline[2,:] = [points[entry][2] for entry in config[4]]; 
            if det(barycenterline) == 0
                global λ4=[0.5,0.5]
            else
                global λ4 = inv(barycenterline)*[2,1];
                global λ4 = collect(λ4 / sum(λ4))
            end

            global λ5 = [0.5,0.5]
            barycenterline[1,:] = [points[entry][1] for entry in config[5]]; barycenterline[2,:] = [points[entry][2] for entry in config[5]]; 
            if det(barycenterline) == 0
                global λ5=[0.5,0.5]
            else
                global λ5 = inv(barycenterline)*[2,1];
                global λ5 = collect(λ5 / sum(λ5))
            end

            global θ1 = prod([(coefficients[config[1][i]]/λ1[i])^(λ1[i]) for i in 1:length(config[1])])
            global θ2 = prod([(coefficients[config[2][i]]/λ2[i])^(λ2[i]) for i in 1:length(config[2])])
            global θ3 = prod([(coefficients[config[3][i]]/λ3[i])^(λ3[i]) for i in 1:length(config[3])])
            global θ4 = prod([(coefficients[config[4][i]]/λ4[i])^(λ4[i]) for i in 1:length(config[4])])
            global θ5 = prod([(coefficients[config[5][i]]/λ5[i])^(λ5[i]) for i in 1:length(config[5])])
            push!(θ,θ1+θ2+θ3+θ4+θ5)
        catch
            continue
        end
    end

    return θ
end

#=
Here, all 16 configurations are plotted and saved. The colors are optimized with respect to distinguishability.
TODO Dissociate the individual plots from the rest to generalize this method
=#
function plotAllCovers(points, triangconfigurations, lineconfigurations)
    fourcolors = map(col -> (red(col), green(col), blue(col)), distinguishable_colors(3, [RGB(1,1,1), RGB(0,0,0)], dropseed=true, lchoices = range(25, stop=50, length=15), hchoices = range(120, stop=350, length=20)))
    fivecolors = map(col -> (red(col), green(col), blue(col)), distinguishable_colors(5, [RGB(1,1,1), RGB(0,0,0)], dropseed=true, lchoices = range(25, stop=60, length=20), hchoices = range(120, stop=330, length=30)))
    #display([[Float64(color[1]),Float64(color[2]),Float64(color[3])]*255 for color in fourcolors])

    pointsForPlot = [(0,0),(1,0),(2,0),(4,1),(4,2),(3,2),(2,2),(0,1),(0,0)]
    fig = Figure(size=(1200,1200))
    ax = Axis(fig[1,1], aspect=1)
    hidespines!(ax)
    hidedecorations!(ax)
    poly!(ax,[Point2f0(pt) for pt in pointsForPlot]; color=RGBA{Float64}(0.1, 0.1, 0.1, 0.04), strokewidth=0)
    lines!(ax,[Point2f0(pt) for pt in pointsForPlot], color=:black, linewidth=10)
    scatter!(ax,[Point2f0([2,1])]; color=:red2,markersize=60)
    text!(ax,[Point2f0([2+0.1,1-0.1])], text=L"$m$"; color=:red2,fontsize=70)
    scatter!(ax,[Point2f0(pt) for pt in points]; color=:black, markersize=60)
    text!(ax,[Point2f0([1+0.1,1-0.1]), Point2f0([3+0.1,1-0.1]), Point2f0([1+0.1,0+0.04]), Point2f0([3-0.05,2-0.19]), Point2f0([0+0.1,0+0.04]), Point2f0([2+0.17,0-0.05]), Point2f0([4-0.1,1-0.17]), Point2f0([4-0.34,2-0.16]), Point2f0([2-0.05,2-0.16]), Point2f0([0+0.1,1-0.1])], text=[L"$\iota_1$", L"$\iota_2$", L"$\beta_1$", L"$\beta_2$", L"$\alpha_1$", L"$\alpha_2$", L"$\alpha_3$", L"$\alpha_4$", L"$\alpha_5$", L"$\alpha_6$"]; color=:black,fontsize=70)

    save("../images/NEWtriangbaseconf.png",fig)

    #This method plots the configurations containing 2-dimensional simplices.
    for config in triangconfigurations
        fig = Figure(size=(1200,1200))
        ax = Axis(fig[1,1], aspect=1)
        hidespines!(ax)
        hidedecorations!(ax)
        poly!(ax,[Point2f0(pt) for pt in pointsForPlot]; color=RGBA{Float64}(0.1, 0.1, 0.1, 0.04),strokewidth=0)
        poly!(ax,[Point2f0(points[pt]) for pt in config[1]]; color=RGBA{Float64}(fourcolors[1][1], fourcolors[1][2], fourcolors[1][3], 0.15),strokewidth=0)
        lines!(ax,[Point2f0(points[pt]) for pt in vcat(config[1],config[1][1])], color=RGBA{Float64}(fourcolors[1][1], fourcolors[1][2], fourcolors[1][3], 1), linewidth=5)
        poly!(ax,[Point2f0(points[pt]) for pt in config[2]]; color=RGBA{Float64}(fourcolors[2][1], fourcolors[2][2], fourcolors[2][3], 0.15),strokewidth=0)
        lines!(ax,[Point2f0(points[pt]) for pt in vcat(config[2],config[2][1])], color=RGBA{Float64}(fourcolors[2][1], fourcolors[2][2], fourcolors[2][3], 1), linewidth=5)
        lines!(ax,[Point2f0(pt) for pt in pointsForPlot], color=:black, linewidth=10)

        firstcolor = RGBA{Float64}(fourcolors[3][1], fourcolors[3][2], fourcolors[3][3], 1)
        secondcolor = RGBA{Float64}(0.35, 0.35, 0.35, 1)
        lw = 15
        if (config[3]==[9,10]||config[3]==[4,8]||config[3]==[10,9]||config[3]==[8,4])&&(config[4]==[9,10]||config[4]==[4,8]||config[4]==[10,9]||config[4]==[8,4])
            lines!(ax,[Point2f0(points[pt]) for pt in [9,10]], color=secondcolor, linewidth=lw)
            xvals=0:0.1:4
            yvals = [1 + 0.0123563*x - 0.630187*x^2 + 0.622464*x^3 - 0.194037*x^4 + 0.0194037*x^5 + 3.66123*10^-17*x^6 for x in xvals]
            lines!(ax,[Point2f0([xvals[i],yvals[i]]) for i in 1:length(yvals)], color=firstcolor, linewidth=lw)
        elseif (config[3]==[8,9]||config[3]==[4,10]||config[3]==[9,8]||config[3]==[10,4])&&(config[4]==[8,9]||config[4]==[4,10]||config[4]==[9,8]||config[4]==[10,4])
            xvals1=0:0.1:3
            xvals2=1:0.1:4
            yvals1=[1 + 0.00765306 *x - 0.389031 *x^2 + 0.320153 *x^3 - 0.0637755 *x^4 for x in xvals1]
            yvals2=yvals1[end:-1:1]
            lines!(ax,[Point2f0([xvals1[i],yvals1[i]]) for i in 1:length(yvals1)], color=firstcolor, linewidth=lw)
            lines!(ax,[Point2f0([xvals2[i],yvals2[i]]) for i in 1:length(yvals2)], color=secondcolor, linewidth=lw)
        elseif config[3]==[8,9]||config[3]==[9,8]||config[4]==[8,9]||config[4]==[9,8]
            xvals1=0:0.1:3
            yvals1=[1 + 0.00765306 *x - 0.389031 *x^2 + 0.320153 *x^3 - 0.0637755 *x^4 for x in xvals1]
            lines!(ax,[Point2f0([xvals1[i],yvals1[i]]) for i in 1:length(yvals1)], color=firstcolor, linewidth=lw)
            twoconfig = config[3]==[8,9]||config[3]==[9,8] ? config[4] : config[3]
            lines!(ax,[Point2f0(points[pt]) for pt in twoconfig], color=secondcolor, linewidth=lw)
        elseif config[3]==[4,10]||config[3]==[10,4]||config[4]==[4,10]||config[4]==[10,4]
            xvals1=0:0.1:3
            xvals2=1:0.1:4
            yvals2=([1 + 0.00765306 *x - 0.389031 *x^2 + 0.320153 *x^3 - 0.0637755 *x^4 for x in xvals1])[end:-1:1]
            lines!(ax,[Point2f0([xvals2[i],yvals2[i]]) for i in 1:length(yvals2)], color=secondcolor, linewidth=lw)
            twoconfig = config[3]==[4,10]||config[3]==[10,4] ? config[4] : config[3]
            lines!(ax,[Point2f0(points[pt]) for pt in twoconfig], color=firstcolor, linewidth=lw)
        elseif (config[3]==[4,8]||config[3]==[8,4]||config[4]==[4,8]||config[3]==[8,4])
            xvals=0:0.1:4
            yvals = [1 + 0.0123563*x - 0.630187*x^2 + 0.622464*x^3 - 0.194037*x^4 + 0.0194037*x^5 + 3.66123*10^-17*x^6 for x in xvals]
            lines!(ax,[Point2f0([xvals[i],yvals[i]]) for i in 1:length(yvals)], color=firstcolor, linewidth=lw)
            twoconfig = config[3]==[4,8]||config[3]==[8,4] ? config[4] : config[3]
            lines!(ax,[Point2f0(points[pt]) for pt in twoconfig], color=secondcolor, linewidth=lw)
        elseif (config[3]==[9,10]||config[3]==[9,10]||config[4]==[9,10]||config[3]==[9,10])
            lines!(ax,[Point2f0(points[pt]) for pt in [9,10]], color=secondcolor, linewidth=lw)
            twoconfig = config[3]==[9,10]||config[3]==[10,9] ? config[4] : config[3]
            lines!(ax,[Point2f0(points[pt]) for pt in twoconfig], color=firstcolor, linewidth=lw)
        else
            lines!(ax,[GPoint2f0(points[pt]) for pt in config[3]], color=firstcolor, linewidth=lw)
            lines!(ax,[Point2f0(points[pt]) for pt in config[4]], color=secondcolor, linewidth=lw)
        end
        
        scatter!(ax,[Point2f0([2,1])]; color=:red2,markersize=60)
        #GLMakie.text!(ax,[GLMakie.Point2f0([2+0.05,1-0.154])], text=L"$\mu$"; color=:red2,fontsize=70)
        scatter!(ax,[Point2f0(pt) for pt in points]; color=:black,markersize=60)
        text!(ax,[Point2f0([0.1,1.8])], text=L"%$(findfirst(t->t==config,triangconfigurations))"; color=:black,fontsize=85)
        save("../images/NEWtriangnoline$(findfirst(t->t==config,triangconfigurations)).png",fig)
    end

    #This method plots the configurations containing no 2-dimensional simplices.
    for config in lineconfigurations
        fig = Figure(size=(1200,1200))
        ax = Axis(fig[1,1], aspect=1)
        hidespines!(ax)
        hidedecorations!(ax)
        firstcolor = RGBA{Float64}(fivecolors[1][1], fivecolors[1][2], fivecolors[1][3], 1)
        secondcolor = RGBA{Float64}(fivecolors[2][1], fivecolors[2][2], fivecolors[2][3], 1)
        lw = 15
        greencolor = RGBA{Float64}(fivecolors[3][1], fivecolors[3][2], fivecolors[3][3], 1)
        bluecolor = RGBA{Float64}(fivecolors[4][1], fivecolors[4][2], fivecolors[4][3], 1)
        yellowcolor = RGBA{Float64}(fivecolors[5][1], fivecolors[5][2], fivecolors[5][3], 1)
        poly!(ax,[Point2f0(pt) for pt in pointsForPlot]; color=RGBA{Float64}(0.1, 0.1, 0.1, 0.035),strokewidth=0)

        if config[3]==[9,10]||config[3]==[10,9]
            lines!(ax,[Point2f0(points[pt]) for pt in [9,10]], color=secondcolor, linewidth=lw)
        elseif config[3]==[8,9]||config[3]==[9,8]
            xvals1=0:0.1:3
            yvals1=[1 + 0.00765306 *x - 0.389031 *x^2 + 0.320153 *x^3 - 0.0637755 *x^4 for x in xvals1]
            lines!(ax,[Point2f0([xvals1[i],yvals1[i]]) for i in 1:length(yvals1)], color=firstcolor, linewidth=lw)
        end

        if config[5]==[4,10]||config[5]==[10,4]
            xvals1=0:0.1:3
            xvals2=1:0.1:4
            yvals2=([1 + 0.00765306 *x - 0.389031 *x^2 + 0.320153 *x^3 - 0.0637755 *x^4 for x in xvals1])[end:-1:1]
            lines!(ax,[Point2f0([xvals2[i],yvals2[i]]) for i in 1:length(yvals2)], color=secondcolor, linewidth=lw)
        elseif config[5]==[4,8]||config[5]==[8,4]
            xvals=0:0.1:4
            yvals = [1 + 0.0123563*x - 0.630187*x^2 + 0.622464*x^3 - 0.194037*x^4 + 0.0194037*x^5 + 3.66123*10^-17*x^6 for x in xvals]
            lines!(ax,[Point2f0([xvals[i],yvals[i]]) for i in 1:length(yvals)], color=firstcolor, linewidth=lw)
        end

        lines!(ax,[Point2f0(points[pt]) for pt in config[1]], color=greencolor, linewidth=lw)
        lines!(ax,[Point2f0(points[pt]) for pt in config[2]], color=bluecolor, linewidth=lw)
        lines!(ax,[Point2f0(points[pt]) for pt in config[4]], color=yellowcolor, linewidth=lw)
        lines!(ax,[Point2f0(pt) for pt in pointsForPlot], color=:black, linewidth=10)

        scatter!(ax,[Point2f0([2,1])]; color=:red2,markersize=60)
        #GLMakie.text!(ax,[GLMakie.Point2f0([2+0.05,1-0.154])], text=L"$\mu$"; color=:red2,fontsize=70)
        scatter!(ax,[Point2f0(pt) for pt in points]; color=:black,markersize=60)
        text!(ax,[Point2f0([0.1,1.8])], text=L"%$(14+findfirst(t->t==config,lineconfigurations))"; color=:black,fontsize=85)
        save("../images/NEWtriangline$(14+findfirst(t->t==config,lineconfigurations)).png",fig)
    end
end

function plottriangle()
    pointsForPlot = [(0,1),(2,0),(4,2)]
    fig = Figure(size=(1200,1200))
    ax = Axis(fig[1,1], aspect=1)
    hidespines!(ax)
    hidedecorations!(ax)
    poly!(ax,[Point2f0(pt) for pt in pointsForPlot]; color=RGBA{Float64}(0.1, 0.1, 0.1, 0.04), strokewidth=0)
    lines!(ax,[Point2f0(pt) for pt in pointsForPlot[vcat(1:3,1)]], color=:black, linewidth=10)
    scatter!(ax,[Point2f0([2,1])]; color=:red2,markersize=60)
    text!(ax,[Point2f0([1.675,0.75])], text=L"$(2,1)$"; color=:red2,fontsize=70)
    scatter!(ax,[Point2f0(pt) for pt in pointsForPlot]; color=:black, markersize=60)
    text!(ax,[Point2f0([-0.32,1.09]), Point2f0([1.665,-0.25]), Point2f0([3.675,2.065])], text=[L"(0,1)", L"(2,0)", L"(4,2)"]; color=:black,fontsize=70)
    xlims!(ax,(-0.3,4.3))
    ylims!(ax,(-0.3,2.3))
    save("../images/triangexample.png",fig)
end

#=
Sample randomly from the region [0,boxsize]^n. Whenever aη>=0 and bη<0, the sample is accepted. 
We compare all samples to the baseline given by `θbaseline`. Whenever the new model recognizes
nonnegativity, while the baseline model does not, we add 1 to `ourmodel` and if the opposite is
true, we add 1 to `prevmodel`. If neither model recognizes nonnegativity, we add 1 to `nomodel`.
`pointnumber` is a counter of the samples drawn.
=#
function runSamplingComparison(θ, κs, Ks, aη, bη, mcoef, θbaseline; boxsize=100, numberOfSamplingRuns=250, prefix="NEW", suffix="")
    #If the file exists, we add to the previously run tests. Else, we set everything to 0.
    global relDict = Dict()
    try
        f = open("../data/$(prefix)$(suffix)triangstoredsolutions$(boxsize).txt", "r")

        global pointnumber = parse(Int,readline(f))
        #global allmodels = parse(Int,readline(f))
        #global onlyonemodel = [parse(Int,entry) for entry in split(readline(f)[2:end-1], ", ")]
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
    catch
        global pointnumber = 0
        global allmodels = 0
        global onlyonemodel = [0 for _ in 1:length(θ)]
        global newmodel = [0 for _ in 1:length(θ)]
        global ourmodel = [0 for _ in 1:length(θ)]
        global prevmodel = [0 for _ in 1:length(θ)]
        global nomodel = [0 for _ in 1:length(θ)]

        for i in 1:length(θ), j in 1:length(θ)
            if i!=j
                #relDict[(i,j)] = 0
            end
        end
    end

    for sampleindex in 1:numberOfSamplingRuns
        display("Run: $(sampleindex)")
        global sampling = filter(sampler -> !any(t->isapprox(t,0), sampler) && evaluate(aη,vcat(Ks,κs)=>sampler)>=0 && evaluate(bη,vcat(Ks,κs)=>sampler)<0, [boxsize * abs.(rand(Float64,length(vcat(Ks,κs)))) for _ in 1:1000000])
        global pointnumber = pointnumber+length(sampling)
        for ind in ProgressBar(1:length(sampling))
            sampler = sampling[ind]
            mval = evaluate(mcoef,vcat(Ks,κs)=>sampler)
            prevval = evaluate(θbaseline, vcat(Ks,κs)=>sampler)
            for j in 1:length(θ)
                ourval = evaluate(θ[j], vcat(Ks,κs)=>sampler)
                if ourval >= -mval
                    global newmodel[j] += 1
                end
                if ourval >= -mval && prevval < -mval
                    global ourmodel[j] += 1
                elseif ourval < -mval && prevval >= -mval
                    global prevmodel[j] += 1
                elseif ourval < -mval && prevval < -mval
                    global nomodel[j] += 1
                end
            end
            #=
            indicator, winnerarray = false, []
            for i in 1:length(θ), j in i+1:length(θ)
                val1 = evaluate(θ[i], vcat(Ks,κs)=>sampler)
                val2 = evaluate(θ[j], vcat(Ks,κs)=>sampler)
                if val1 >= -mval && val2 < -mval
                    push!(winnerarray, i)
                    relDict[(i,j)] += 1
                elseif val2 >= -mval && val1 < -mval
                    push!(winnerarray, j)
                    relDict[(j,i)] += 1
                end

                if val2 >= -mval || val1 >= -mval
                    indicator = true
                end
            end
            if length(collect(Set(winnerarray))) == 1
                onlyonemodel[winnerarray[1]]+=1
            end
            global allmodels += indicator ? 1 : 0 
            =#
        end

        #SAVE the data to the file `NEWtriangstoredsolutions.txt`
        open("../data/$(prefix)$(suffix)triangstoredsolutions$(boxsize).txt", "w") do file
            write(file, "$(pointnumber)\n")
            #write(file, "$(allmodels)\n")
            #write(file, "$(onlyonemodel)\n")
            write(file, "$(newmodel)\n")
            write(file, "$(ourmodel)\n")
            write(file, "$(prevmodel)\n")
            write(file, "$(nomodel)\n")
            #=for key in keys(relDict)
                write(file, "$(key[1]), $(key[2]): $(relDict[key])\n")
            end=#
        end
    end

    #foreach(j->println("Case $(j): Our model performed better in $(100*round(ourmodel[j]/(ourmodel[j]+prevmodel[j]),5))% of the cases, where the other model did not work. No model found anything in $(100*round(nomodel[j]/pointnumber, 5)) of the cases."), 1:length(θ))
end

function printValues( ; prefix="relTest", suffix="NEW")
    for boxsize in [1,10,100,1000]
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
        println(onlyonemodel)
        containmentArray, almostContainmentArray = [[] for _ in 1:16], [[] for _ in 1:16]
        for key in keys(relDict)
            if relDict[key]==0
                push!(containmentArray[key[1]], key[2])
            elseif relDict[key]<3
                push!(almostContainmentArray[key[1]], key[2])
            end
        end
        println("Containment")
        foreach(t->println("$(t): $(containmentArray[t])"), 1:16)
        println("AlmostContainment")
        foreach(t->println("$(t): $(almostContainmentArray[t])"), 1:16)
    end
end
#=
This is the main method. Use it to run all tests.
=#
function runTest( ; boxsize=1000, numberOfSamplingRuns=100, prefix="relTest", suffix="NEW")
    @var K[1:4] κ[1:12]

    #We choose colors with maximum distinguishability
    hexPoints = [(0,0),(1,0),(2,0),(4,1),(4,2),(3,2),(2,2),(0,1),(3,1),(1,1)]
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
    #Check if all covers are legit
    for config in triangconfigurations
        all(t->t in union(config[1],config[2],config[3],config[4]),1:10) || display(config)&&throw(error("The triangles don't cover the entire region"))
        isempty(intersect(config[1],config[2])) && isempty(intersect(config[1],config[3])) && isempty(intersect(config[1],config[4])) && isempty(intersect(config[2],config[3])) && isempty(intersect(config[2],config[4])) && isempty(intersect(config[3],config[4])) || display(config)&&throw(error("Each vertex should only be used once"))
    end

    lineconfigurations = [[[5,1],[7,3],[8,9],[6,2],[4,10]], [[5,1],[7,3],[9,10],[6,2],[8,4]]]

    θ = createθcircuits(hexPoints, coefficients, lineconfigurations, triangconfigurations)
    #plotAllCovers(hexPoints, triangconfigurations, lineconfigurations)
    runSamplingComparison(θ, [κ[3],κ[6],κ[9],κ[12]], K, aη, bη, mcoef, θ[9]; boxsize=boxsize, numberOfSamplingRuns=numberOfSamplingRuns, prefix=prefix, suffix=suffix)
end

function runTest( ; boxsize=1000, numberOfSamplingRuns=100, prefix="relTest", suffix="NEW")
    @var K[1:4] κ[1:12]

    #We choose colors with maximum distinguishability
    hexPoints = [(0,0),(1,0),(2,0),(4,1),(4,2),(3,2),(2,2),(0,1),(3,1),(1,1)]
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
    #Check if all covers are legit
    for config in triangconfigurations
        all(t->t in union(config[1],config[2],config[3],config[4]),1:10) || display(config)&&throw(error("The triangles don't cover the entire region"))
        isempty(intersect(config[1],config[2])) && isempty(intersect(config[1],config[3])) && isempty(intersect(config[1],config[4])) && isempty(intersect(config[2],config[3])) && isempty(intersect(config[2],config[4])) && isempty(intersect(config[3],config[4])) || display(config)&&throw(error("Each vertex should only be used once"))
    end

    lineconfigurations = [[[5,1],[7,3],[8,9],[6,2],[4,10]], [[5,1],[7,3],[9,10],[6,2],[8,4]]]

    θ = createθcircuits(hexPoints, coefficients, lineconfigurations, triangconfigurations)
    #plotAllCovers(hexPoints, triangconfigurations, lineconfigurations)
    runSamplingComparison(θ, [κ[3],κ[6],κ[9],κ[12]], K, aη, bη, mcoef, θ[9]; boxsize=boxsize, numberOfSamplingRuns=numberOfSamplingRuns, prefix=prefix, suffix=suffix)
end


function runTest_noDependencies( ; boxsize=1000, numberOfSamplingRuns=62, prefix="", suffix="noDependencies")
    @var κ[1:11]

    hexPoints = [(0,0),(1,0),(2,0),(4,1),(4,2),(3,2),(2,2),(0,1),(3,1),(1,1)]
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

    θ = createθcircuits(hexPoints, Vector{Expression}(κ[1:10]), lineconfigurations, triangconfigurations)
    runSamplingComparison(θ, κ[1:11], Vector{Float64}([]), Expression(1), Expression(-1), Expression(-κ[11]), θ[9]; boxsize=boxsize, numberOfSamplingRuns=numberOfSamplingRuns, prefix=prefix, suffix=suffix)
end

function computeCoverInvariants( ; startboxsize=1, finalboxsize=1000, prefix="NEW", suffix="")
    fig = Figure(size=(1300,500))
    ax_ourmodel = Axis(fig[1,1]; title="Our Model Won", xlabel=L"$\log(b)+1$", ylabel=L"$\perthousand$")
    ax_prevmodel = Axis(fig[1,2]; title="Baseline Model Won", xlabel=L"$\log(b)+1$", ylabel=L"$\perthousand$")
    ax_nomodel = Axis(fig[1,3]; title="Neither Model Won", xlabel=L"$\log(b)+1$", ylabel=L"$\perthousand$")
    #hidedecorations!(ax_ourmodel); hidedecorations!(ax_prevmodel); hidedecorations!(ax_nomodel);
    allmodeldots, ourmodeldots, prevmodeldots, nomodeldots, puremodel = [Vector{Float64}([]) for _ in 1:16], [Vector{Float64}([]) for _ in 1:16], [Vector{Float64}([]) for _ in 1:16], [Vector{Float64}([]) for _ in 1:16], [Vector{Float64}([]) for _ in 1:16]
    xaxis = Vector{Float64}([])

    for boxsize in startboxsize:finalboxsize
        try
            f = open("../data/$(prefix)$(suffix)triangstoredsolutions$(boxsize).txt", "r")
            
            global pointnumber = parse(Int,readline(f))
            global allpoints = [parse(Int,entry) for entry in split(readline(f)[2:end-1], ", ")]
            global ourmodel = [parse(Int,entry) for entry in split(readline(f)[2:end-1], ", ")]
            global prevmodel = [parse(Int,entry) for entry in split(readline(f)[2:end-1], ", ")]
            global nomodel = [parse(Int,entry) for entry in split(readline(f)[2:end-1], ", ")]

            close(f)
        catch
            continue
        end

        permille_allmodels, permille_ourmodel, permille_prevmodel, permille_nomodel, percent_ourmodel_pure = round.(round.(10000*allpoints ./ pointnumber, digits=2)/10000, digits=6), round.(10000*ourmodel ./ pointnumber, digits=2), round.(10000*prevmodel ./ pointnumber, digits=2), round.(10000*nomodel ./ pointnumber, digits=2), round.(100 * (pointnumber .- (nomodel .+ prevmodel)) ./ pointnumber, digits=3)
        foreach(i->push!(allmodeldots[i], permille_allmodels[i]), 1:length(permille_allmodels))
        foreach(i->push!(ourmodeldots[i], permille_ourmodel[i]), 1:length(permille_ourmodel))
        foreach(i->push!(prevmodeldots[i], permille_prevmodel[i]), 1:length(permille_ourmodel))
        foreach(i->push!(nomodeldots[i], permille_nomodel[i]), 1:length(permille_ourmodel))
        foreach(i->push!(puremodel[i], percent_ourmodel_pure[i]), 1:length(percent_ourmodel_pure))

        push!(xaxis, log(boxsize)+1)
        
        println("$(boxsize):\t ourmodel\t prevmodel\t nomodel\t ourmodel_pure[%]")
        foreach(θ -> println("$(θ): \t$(permille_ourmodel[θ]) \t$(permille_prevmodel[θ]) \t$(permille_nomodel[θ]) \t$(percent_ourmodel_pure[θ])"), 1:length(ourmodel))
        println("\n")
    end

    #print Latex table code
    #=
    for θ in 1:Int(length(ourmodeldots)/2)
        print("~&+&")
        for i in 1:length(ourmodeldots[θ])
            print("$(ourmodeldots[θ][i])&$(ourmodeldots[θ+8][i])&")
        end
        print("+&~\\\\ \n")
        print("$(θ)&\$-\$&")
        for i in 1:length(prevmodeldots[θ])
            print("$(prevmodeldots[θ][i])&$(prevmodeldots[θ+8][i])&")
        end
        print("\$-\$&$(θ+8)\\\\ \n")
        print("~&0&")
        for i in 1:length(nomodeldots[θ])
            print("$(nomodeldots[θ][i])&$(nomodeldots[θ+8][i])&")
        end
        print("0&~\\\\[.3mm] \\thickhline \n\n")
    end
    =#
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

    colors = colormap("Blues", length(ourmodeldots); logscale=false)
    foreach(line->lines!(ax_ourmodel, xaxis, ourmodeldots[line] ./ ourmodeldots[line][1]; linewidth=4, color = colors[line]), 1:length(ourmodeldots))
    foreach(line->lines!(ax_prevmodel, xaxis, prevmodeldots[line] ./ prevmodeldots[line][1]; linewidth=4, color = colors[line]), 1:length(prevmodeldots))
    foreach(line->lines!(ax_nomodel, xaxis, nomodeldots[line] ./ nomodeldots[line][1]; linewidth=4, color = colors[line]), 1:length(nomodeldots))
    save("../images/$(prefix)$(suffix)$(length(ourmodeldots))cover_curveplots.png", fig)
end

function compareTwoCovers(cover_suggested::Int, cover_baseline::Int; numberOfSamplingRuns=100, boxsizes=[5 for _ in 1:8])
    @var K[1:4] κ[1:12]

    #We choose colors with maximum distinguishability
    hexPoints = [(0,0),(1,0),(2,0),(4,1),(4,2),(3,2),(2,2),(0,1),(3,1),(1,1)]
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
    #Check if all covers are legit
    for config in triangconfigurations
        all(t->t in union(config[1],config[2],config[3],config[4]),1:10) || display(config)&&throw(error("The triangles don't cover the entire region"))
        isempty(intersect(config[1],config[2])) && isempty(intersect(config[1],config[3])) && isempty(intersect(config[1],config[4])) && isempty(intersect(config[2],config[3])) && isempty(intersect(config[2],config[4])) && isempty(intersect(config[3],config[4])) || display(config)&&throw(error("Each vertex should only be used once"))
    end

    lineconfigurations = [[[5,1],[7,3],[8,9],[6,2],[4,10]], [[5,1],[7,3],[9,10],[6,2],[8,4]]]

    θ = createθcircuits(hexPoints, K, κ, coefficients, lineconfigurations, triangconfigurations)
    empiricalComparisonOfTwoCovers(θ[cover_suggested], θ[cover_baseline], K, κ, aη, bη, mcoef; numberOfSamplingRuns=250, boxsizes=[5 for _ in 1:8], cover_1=cover_suggested, cover_2 = cover_baseline)
end

function empiricalComparisonOfTwoCovers(θsuggestion, θbaseline, K, κ, aη, bη, mcoef; numberOfSamplingRuns=100, boxsizes=[5 for _ in 1:8], cover_1=15, cover_2 = 9)
    vector_baseline_wins = []
    vector_our_wins = []
    vector_both_wins = []
    vector_no_wins = []
    for sampleindex in 1:numberOfSamplingRuns
        display("Run: $(sampleindex)")
        global sampling = filter(sampler -> !any(t->isapprox(t,0), sampler) && evaluate(aη,vcat(K,[κ[3],κ[6],κ[9],κ[12]])=>sampler)>=0 && evaluate(bη,vcat(K,[κ[3],κ[6],κ[9],κ[12]])=>sampler)<0, [boxsizes .* abs.(rand(Float64,8)) for _ in 1:1000000])
        for ind in ProgressBar(1:length(sampling))
            sampler = sampling[ind]
            mval = evaluate(mcoef,vcat(K,[κ[3],κ[6],κ[9],κ[12]])=>sampler)
            prevval = evaluate(θbaseline, vcat(K,[κ[3],κ[6],κ[9],κ[12]])=>sampler)
            ourval = evaluate(θsuggestion, vcat(K,[κ[3],κ[6],κ[9],κ[12]])=>sampler)

            if ourval >= -mval && prevval < -mval
                push!(vector_our_wins, sampler)
            elseif ourval < -mval && prevval >= -mval
                push!(vector_baseline_wins, sampler)
            elseif ourval >= -mval && prevval >= -mval
                push!(vector_both_wins, sampler)
            else
                push!(vector_no_wins, sampler)
            end
        end
    end

    fig = Figure(size=(1200,1200))
    ax = [Axis(fig[1,1]; xlabel = L"K_1", ylabel=L"K_2"), Axis(fig[1,2]; xlabel = L"K_3", ylabel=L"K_4"), Axis(fig[2,1]; xlabel = L"$\kappa_3$", ylabel=L"$\kappa_6$"), Axis(fig[2,2]; xlabel = L"$\kappa_9$", ylabel=L"$\kappa_{12}$")]
    for pic in 1:4
        scatter!(ax[pic], [Point2f0(pt[(2*(pic-1)+1):(2*(pic-1)+2)]) for pt in vector_both_wins]; color=:lightgrey, markersize=1, markerstrokewidth=0)
        scatter!(ax[pic], [Point2f0(pt[(2*(pic-1)+1):(2*(pic-1)+2)]) for pt in vector_no_wins]; color=:red3, markersize=1, markerstrokewidth=0)
        scatter!(ax[pic], [Point2f0(pt[(2*(pic-1)+1):(2*(pic-1)+2)]) for pt in vector_baseline_wins]; color=:blue3, markersize=1, markerstrokewidth=0)
        scatter!(ax[pic], [Point2f0(pt[(2*(pic-1)+1):(2*(pic-1)+2)]) for pt in vector_our_wins]; color=:green3, markersize=1, markerstrokewidth=0)
    end
    save("../images/bestcoverplots$(cover_1)-$(cover_2).png", fig)
end

function plotNewtonPolytope()
    vertices = [[4.,0,2], [2.,2,2], [4.,0,1], [3.,2,1], [2.,3,1], [0.,4,1], [2.,3,0], [2.,2,0], [1.,4,0], [0.,4,0]]
    p = Mesh(polyhedron(convexhull(vertices...)))
    fig = Figure(size = (1200,1200);  fontsize=22  )
    ax = Axis3(fig[1,1], aspect=(1.,1,0.5);)  
    hidedecorations!(ax, label = false, ticklabels = false, ticks = false, grid = true)  
    mesh!(ax, p, color=:steelblue; )
    scatter!(ax, Point3f0([1.98,1.98,1]); markersize=30, color=:red3)
    display(fig)
end

function createθcircuits_weighted(points, coefficients, configurations; discretization=25)
    θdict = Dict()
    length(configurations)>=3 || throw(error("Since we are providing a 2D heatmap of the covers, at least 3 simplicial configurations need to be provided!"))
    for i in 1:length(configurations), j in i+1:length(configurations), k in j+1:length(configurations)
        ijkDict = Dict()
        for ω1 in 0:1/discretization:1
            ijkDict[ω1] = []
        end

        for ω1 in 0:round(1/discretization, sigdigits=5):1, ω2 in 0:round(1/discretization, sigdigits=5):round(1-ω1, sigdigits=5)
            configuration1, configuration2, configuration3, helper = Base.copy(configurations[i]), Base.copy(configurations[j]), Base.copy(configurations[k]), []
            while !isempty(configuration1)
                config = pop!(configuration1)
                if length(config)==2
                    global barycenter = Matrix{Float64}(undef,2,2); barycenter[1,:] = [points[entry][1] for entry in config]; barycenter[2,:] = [points[entry][2] for entry in config]; 
                    global λ = (det(barycenter) == 0) ? [0.5,0.5] : collect(inv(barycenter)*[2,1] / sum(inv(barycenter)*[2,1]))

                    if config in configuration2 && config in configuration3
                        deleteat!(configuration2, findfirst(entry -> config==entry, configuration2))
                        deleteat!(configuration3, findfirst(entry -> config==entry, configuration3))
                        push!(helper, prod([(coefficients[config[i]]/λ[i])^(λ[i]) for i in 1:length(config)]))
                    elseif config in configuration2 || config in configuration3
                        (config in configuration2) ? deleteat!(configuration2, findfirst(entry -> config==entry, configuration2)) : deleteat!(configuration3, findfirst(entry -> config==entry, configuration3))
                        push!(helper, prod([(((config in configuration2) ? ω1+ω2 : 1-ω2)*coefficients[config[i]]/λ[i])^(λ[i]) for i in 1:length(config)]))
                    else
                        push!(helper, prod([(ω1*coefficients[config[i]]/λ[i])^(λ[i]) for i in 1:length(config)]))
                    end
                elseif length(config)==3
                    global barycenter, msolve = Matrix{Float64}(undef,3,3), [2,1,1]; barycenter[1,:] = [points[entry][1] for entry in config]; barycenter[2,:] = [points[entry][2] for entry in config]; barycenter[3,:] = [1 for entry in config];
                    global λ = collect(inv(barycenter)*msolve / sum(inv(barycenter)*msolve));
                    if config in configuration2 && config in configuration3
                        deleteat!(configuration2, findfirst(entry -> config==entry, configuration2))
                        push!(helper, prod([(coefficients[config[i]]/λ[i])^(λ[i]) for i in 1:length(config)]))
                    elseif config in configuration2 || config in configuration3
                        (config in configuration2) ? deleteat!(configuration2, findfirst(entry -> config==entry, configuration2)) : deleteat!(configuration3, findfirst(entry -> config==entry, configuration3))
                        push!(helper, prod([(((config in configuration2) ? ω1+ω2 : 1-ω2)*coefficients[config[i]]/λ[i])^(λ[i]) for i in 1:length(config)]))
                    else
                        push!(helper, prod([(ω1*coefficients[config[i]]/λ[i])^(λ[i]) for i in 1:length(config)]))
                    end
                end
            end

            while !isempty(configuration2)
                config = pop!(configuration2)
                if length(config)==2
                    global barycenter = Matrix{Float64}(undef,2,2); barycenter[1,:] = [points[entry][1] for entry in config]; barycenter[2,:] = [points[entry][2] for entry in config]; 
                    global λ = (det(barycenter) == 0) ? [0.5,0.5] : collect(inv(barycenter)*[2,1] / sum(inv(barycenter)*[2,1]))
                    if config in configuration3
                        deleteat!(configuration3, findfirst(entry -> config==entry, configuration3))
                        push!(helper, prod([((1-ω1)*coefficients[config[i]]/λ[i])^(λ[i]) for i in 1:length(config)]))
                    else
                        push!(helper, prod([(ω2*coefficients[config[i]]/λ[i])^(λ[i]) for i in 1:length(config)]))
                    end
                elseif length(config)==3
                    global barycenter, msolve = Matrix{Float64}(undef,3,3), [2,1,1]; barycenter[1,:] = [points[entry][1] for entry in config]; barycenter[2,:] = [points[entry][2] for entry in config]; barycenter[3,:] = [1 for entry in config];
                    global λ = collect(inv(barycenter)*msolve / sum(inv(barycenter)*msolve));
                    if config in configuration3
                        deleteat!(configuration3, findfirst(entry -> config==entry, configuration3))
                        push!(helper, prod([((1-ω1)*coefficients[config[i]]/λ[i])^(λ[i]) for i in 1:length(config)]))
                    else
                        push!(helper, prod([(ω2*coefficients[config[i]]/λ[i])^(λ[i]) for i in 1:length(config)]))
                    end
                end
            end

            while !isempty(configuration3)
                config = pop!(configuration3)
                if length(config)==2
                    global barycenter = Matrix{Float64}(undef,2,2); barycenter[1,:] = [points[entry][1] for entry in config]; barycenter[2,:] = [points[entry][2] for entry in config]; 
                    global λ = (det(barycenter) == 0) ? [0.5,0.5] : collect(inv(barycenter)*[2,1] / sum(inv(barycenter)*[2,1]))
                    push!(helper, prod([((1-ω1-ω2)*coefficients[config[i]]/λ[i])^(λ[i]) for i in 1:length(config)]))
                elseif length(config)==3
                    global barycenter, msolve = Matrix{Float64}(undef,3,3), [2,1,1]; barycenter[1,:] = [points[entry][1] for entry in config]; barycenter[2,:] = [points[entry][2] for entry in config]; barycenter[3,:] = [1 for entry in config];
                    global λ = collect(inv(barycenter)*msolve / sum(inv(barycenter)*msolve));
                    push!(helper, prod([((1-ω1-ω2)*coefficients[config[i]]/λ[i])^(λ[i]) for i in 1:length(config)]))
                end
            end

            push!(ijkDict[ω1], sum(helper))
        end
        θdict[(i,j,k)] = ijkDict
        println(length.(values(θdict[(i,j,k)])))
    end
    return θdict
end

function runSamplingComparison_weighted(θ, θ_weighted, κs, Ks, aη, bη, mcoef, θbaseline; discretization, boxsize=100, numberOfSamplingRuns=250, prefix="NEW", suffix="")
    #If the file exists, we add to the previously run tests. Else, we set everything to 0.
    θkeys = keys(θ_weighted)
    ourmodel = Dict()
    no_other = Dict()
    try
        f = open("../data/$(prefix)$(suffix)triangstoredsolutions$(boxsize).txt", "r")

        global pointnumber = parse(Int,readline(f))
        while ! eof(f)  
            sstring = split(readline(f), ": ")
            keystring = split(sstring[1], "; ")
            key = (parse(Int,keystring[1][2]), (parse(Int,keystring[1][5])), (parse(Int,keystring[1][8])))
            weight = parse(Float64,keystring[2])
            ourmodel[(key, weight)] = [parse(Int,entry) for entry in split(sstring[2][2:end-1], ", ")]
            no_other[(key, weight)] = [0 for _ in 1:length(ourmodel[(key, weight)])]
         end
        close(f)
    catch
        global pointnumber = 0
        for key in θkeys
            ijkkeys = keys(θ_weighted[key])
            for weight in ijkkeys
                ourmodel[(key, weight)] = [0 for _ in 1:length(θ_weighted[key][weight])]
                no_other[(key, weight)] = [0 for _ in 1:length(θ_weighted[key][weight])]
            end
        end
    end

    for sampleindex in 1:numberOfSamplingRuns
        display("Run: $(sampleindex)")
        global sampling = filter(sampler -> !any(t->isapprox(t,0), sampler) && evaluate(aη,vcat(Ks,κs)=>sampler)>=0 && evaluate(bη,vcat(Ks,κs)=>sampler)<0, [boxsize * abs.(rand(Float64,length(vcat(Ks,κs)))) for _ in 1:1000000])
        global pointnumber = pointnumber+length(sampling)
        for ind in ProgressBar(1:length(sampling))
            sampler = sampling[ind]
            mval = real(evaluate(mcoef,vcat(Ks,κs)=>sampler))
            θvals = real.(evaluate.(θ, vcat(Ks,κs)=>sampler))
            for key in θkeys
                for weight in keys(θ_weighted[key])
                    for j in 1:length(θ_weighted[key][weight])
                        ourval = evaluate(θ_weighted[key][weight][j], vcat(Ks,κs)=>sampler)
                        if real(ourval) >= -mval
                            global ourmodel[(key,weight)][j] += 1

                            if all(t-> t<-mval, θvals)
                                global no_other[(key,weight)][j] += 1
                            end
                        end
                    end
                end
            end

        end

        #SAVE the data to the file `NEWtriangstoredsolutions.txt`
        open("../data/$(prefix)$(suffix)triangstoredsolutions$(boxsize).txt", "w") do file
            write(file, "$(pointnumber)\n")
            for key in keys(ourmodel)
                write(file, "$(key[1]); $(key[2]): $(ourmodel[key])\n")
            end
            write(file, "\n NO OTHER\n")
            for key in keys(no_other)
                write(file, "$(key[1]); $(key[2]): $(no_other[key])\n")
            end
        end
    end

    #foreach(j->println("Case $(j): Our model performed better in $(100*round(ourmodel[j]/(ourmodel[j]+prevmodel[j]),5))% of the cases, where the other model did not work. No model found anything in $(100*round(nomodel[j]/pointnumber, 5)) of the cases."), 1:length(θ))
end

function plotWeightedCovers(; boxsize=1, prefix="TWOBEST", suffix="10,12,15")
    helperDict, ourmodel = Dict(), Dict()
    try
        f = open("../data/$(prefix)$(suffix)triangstoredsolutions$(boxsize).txt", "r")

        global pointnumber = parse(Int,readline(f))
        while ! eof(f)  
            str = readline(f)
            display(str)
            if str==""
                break
            end
            sstring = split(str, ": ")
            keystring = split(sstring[1], "; ")
            key = (parse(Int,keystring[1][2]), (parse(Int,keystring[1][5])), (parse(Int,keystring[1][8])))
            weight = parse(Float64,keystring[2])
            helperDict[(key, weight)] = [parse(Int,entry) for entry in split(sstring[2][2:end-1], ", ")]
         end
        close(f)
    catch
        throw(error("No Data recorded!"))
    end
    for key in keys(helperDict)
        if key[1] in keys(ourmodel)
            ourmodel[key[1]][key[2]] = helperDict[key]
        else
            ourmodel[key[1]] = Dict(key[2]=>helperDict[key])
        end
    end

    global maxval, minval = maximum(vcat([vcat(values(ourmodel[key])...) for key in keys(ourmodel)]...)), minimum(vcat([vcat(values(ourmodel[key])...) for key in keys(ourmodel)]...))

    for key in keys(ourmodel)
        heatMatrix = Matrix{Float64}(undef,maximum(length.(values(ourmodel[key]))),maximum(length.(values(ourmodel[key])))); 
        heatMatrix .= NaN
        display(length.(values(ourmodel[key])))
        global row = 1
        for weight in sort(collect(keys(ourmodel[key])))
            heatMatrix[row,maximum(length.(values(ourmodel[key])))-length(ourmodel[key][weight])+1:maximum(length.(values(ourmodel[key])))] = ourmodel[key][weight] ./ pointnumber
            row += 1
        end
        heatMatrix .= heatMatrix[:, end:-1:1]

        fig = Figure(size=(1200,1000), fontsize=28)
        ax=Axis(fig[1,1])
        hm = heatmap!(ax, heatMatrix; colormap=:viridis)#, colorrange=(minval/pointnumber, maxval/pointnumber))
        hidespines!(ax)
        hidedecorations!(ax)
        xlims!(ax, (-0.01,18))
        ylims!(ax, (-0.25,18.25))
        text!(ax, [0.2,0.2,17.3], [-0.25,17.55,-0.25]; text=[L"10", L"15", L"4"], fontsize=35)
        Colorbar(fig[:, end+1], hm; size=30)
        save("../images/discretizedHeatmap4,10,15.png",fig)
        pointarray = []
        for weight in keys(ourmodel[key])
            for i in 1:length(ourmodel[key][weight])
                push!(pointarray, [weight, (i-1)/(maximum(length.(values(ourmodel[key])))-1), ourmodel[key][weight][i] / pointnumber])
            end
        end
        println(key)
        println(pointarray)
    end
end

function runTest_twoBestCovers(; boxsize=1, numberOfSamplingRuns=500, prefix="TWOBEST", suffix="4,10,15", discretization=16)
    @var K[1:4] κ[1:12]

    #We choose colors with maximum distinguishability
    hexPoints = [(0,0),(1,0),(2,0),(4,1),(4,2),(3,2),(2,2),(0,1),(3,1),(1,1)]
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
    oldθ = createθcircuits(hexPoints, coefficients, [], [[[3,5,8],[1,4,7],[9,10],[2,6]]])[1]
    θ = createθcircuits(hexPoints, coefficients, lineconfigurations, triangconfigurations)
    θ_weighted = createθcircuits_weighted(hexPoints, coefficients, configurations; discretization=discretization)
    #plotAllCovers(hexPoints, triangconfigurations, lineconfigurations)
    runSamplingComparison_weighted(θ, θ_weighted, [κ[3],κ[6],κ[9],κ[12]], K, aη, bη, mcoef, oldθ; boxsize=boxsize, numberOfSamplingRuns=numberOfSamplingRuns, prefix=prefix, suffix=suffix, discretization=discretization)    
end
plotWeightedCovers(; boxsize=10, suffix="4,10,15")
#TODO Linear Coefficients test (over all regions?)

#TODO How big of a region can we cover if all covers are used??? Does the best performing cover contain any points not covered by any other cover?
#TODO is there a completely redundant cover?

#TODO Do the weighted cover cover points that are not covered by any other cover?

#TODO different scales for K1,...kappa12,K1,...,K4? Probably not.

#TODO Do the coefficients of the best-performing covers tell us anything that could be generalized?
end 
