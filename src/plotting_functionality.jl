module plotting_functionality

import GLMakie: axislegend, rotate!, text!, Colorbar, heatmap, heatmap!, xlims!, ylims!, plot!, Point2f0, lines!, Figure, Axis, save, hidespines!, hidedecorations!, mesh!, scatter!, text!, RGBA, RGB, poly!, Axis3, Point3f0
import Colors: distinguishable_colors, red, green, blue, colormap
import Plots: cgrad
import Polyhedra: Mesh, polyhedron, convexhull

#=
Here, all 16 configurations are plotted and saved. The colors are optimized with respect to distinguishability.
=#
function plotAllCovers(points, triangconfigurations, lineconfigurations)
    fourcolors = map(col -> (red(col), green(col), blue(col)), distinguishable_colors(3, [RGB(1,1,1), RGB(0,0,0)], dropseed=true, lchoices = range(25, stop=50, length=15), hchoices = range(120, stop=350, length=20)))
    fivecolors = map(col -> (red(col), green(col), blue(col)), distinguishable_colors(5, [RGB(1,1,1), RGB(0,0,0)], dropseed=true, lchoices = range(25, stop=60, length=20), hchoices = range(120, stop=330, length=30)))

    pointsForPlot = [(0,0),(1,0),(2,0),(4,1),(4,2),(3,2),(2,2),(0,1),(0,0)]
    fig = Figure(size=(1200,1200))
    ax = Axis(fig[1,1], aspect=1)
    hidespines!(ax)
    hidedecorations!(ax)
    poly!(ax,[Point2f0(pt) for pt in pointsForPlot]; color=RGBA{Float64}(0.1, 0.1, 0.1, 0.04), strokewidth=0)
    lines!(ax,[Point2f0(pt) for pt in pointsForPlot], color=:black, linewidth=10)
    scatter!(ax,[Point2f0([2,1])]; color=:red2,markersize=60)
    text!(ax,[Point2f0([2+0.1,1-0.1])], text="m"; color=:red2,fontsize=70)
    scatter!(ax,[Point2f0(pt) for pt in points]; color=:black, markersize=60)
    text!(ax,[Point2f0([1+0.1,1-0.1]), Point2f0([3+0.1,1-0.1]), Point2f0([1+0.1,0+0.04]), Point2f0([3-0.05,2-0.19]), Point2f0([0+0.1,0+0.04]), Point2f0([2+0.17,0-0.05]), Point2f0([4-0.1,1-0.17]), Point2f0([4-0.34,2-0.16]), Point2f0([2-0.05,2-0.16]), Point2f0([0+0.1,1-0.1])], text=["iota_1", "iota_2", "beta_1", "beta_2", "alpha_1", "alpha_2", "alpha_3", "alpha_4", "alpha_5", "alpha_6"]; color=:black,fontsize=70)

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
        scatter!(ax,[Point2f0(pt) for pt in points]; color=:black,markersize=60)
        text!(ax,[Point2f0([0.1,1.8])], text="CC({$(findfirst(t->t==config,triangconfigurations))})"; color=:black,fontsize=85)
        save("../images/NEWtriangnoline(findfirst(t->t==config,triangconfigurations)).png",fig)
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
        scatter!(ax,[Point2f0(pt) for pt in points]; color=:black,markersize=60)
        text!(ax,[Point2f0([0.1,1.8])], text="CC({$(14+findfirst(t->t==config,lineconfigurations))})"; color=:black,fontsize=85)
        save("../images/NEWtriangline(14+findfirst(t->t==config,lineconfigurations)).png",fig)
    end
end

function plot_intermediate_covers(points, conf1, conf2)
    fourcolors = map(col -> (red(col), green(col), blue(col)), distinguishable_colors(3, [RGB(1,1,1), RGB(0,0,0)], dropseed=true, lchoices = range(25, stop=50, length=15), hchoices = range(120, stop=350, length=20)))
    fivecolors = map(col -> (red(col), green(col), blue(col)), distinguishable_colors(5, [RGB(1,1,1), RGB(0,0,0)], dropseed=true, lchoices = range(25, stop=60, length=20), hchoices = range(120, stop=330, length=30)))

    pointsForPlot = [(0,0),(1,0),(2,0),(4,1),(4,2),(3,2),(2,2),(0,1),(0,0)]
    for t in 0.0:1/3:1.0
        fig = Figure(size=(1200,1200))
        ax = Axis(fig[1,1], aspect=1)
        hidespines!(ax)
        hidedecorations!(ax)
        poly!(ax,[Point2f0(pt) for pt in pointsForPlot]; color=RGBA{Float64}(0.1, 0.1, 0.1, 0.04),strokewidth=0)
        if t < 0.5
            poly!(ax,[Point2f0(points[pt]) for pt in conf2[1]]; color=RGBA{Float64}(fourcolors[1][1], fourcolors[1][2], fourcolors[1][3], t==0 ? 0.15 : 0.075),strokewidth=0)
        else
            poly!(ax,[Point2f0(points[pt]) for pt in conf1[1]]; color=RGBA{Float64}(fourcolors[1][1], fourcolors[1][2], fourcolors[1][3], t==1 ? 0.15 : 0.075),strokewidth=0)
        end
        lines!(ax,[Point2f0(points[pt]) for pt in vcat(conf1[1],conf1[1][1])], color=RGBA{Float64}(fourcolors[1][1], fourcolors[1][2], fourcolors[1][3], t), linewidth=5)
        lines!(ax,[Point2f0(points[pt]) for pt in vcat(conf2[1],conf2[1][1])], color=RGBA{Float64}(fourcolors[1][1], fourcolors[1][2], fourcolors[1][3], 1-t), linewidth=5)
        if t < 0.5
            poly!(ax,[Point2f0(points[pt]) for pt in conf2[2]]; color=RGBA{Float64}(fourcolors[2][1], fourcolors[2][2], fourcolors[2][3], t==0 ? 0.15 : 0.075),strokewidth=0)
        else
            poly!(ax,[Point2f0(points[pt]) for pt in conf1[2]]; color=RGBA{Float64}(fourcolors[2][1], fourcolors[2][2], fourcolors[2][3],  t==1 ? 0.15 : 0.075),strokewidth=0)
        end
        lines!(ax,[Point2f0(points[pt]) for pt in vcat(conf1[2],conf1[2][1])], color=RGBA{Float64}(fourcolors[2][1], fourcolors[2][2], fourcolors[2][3], t), linewidth=5)
        lines!(ax,[Point2f0(points[pt]) for pt in vcat(conf2[2],conf2[2][1])], color=RGBA{Float64}(fourcolors[2][1], fourcolors[2][2], fourcolors[2][3], 1-t), linewidth=5)
        lines!(ax,[Point2f0(pt) for pt in pointsForPlot], color=:black, linewidth=10)

        lw = 15
        if (conf1[3]==[9,10]||conf1[3]==[4,8]||conf1[3]==[10,9]||conf1[3]==[8,4])&&(conf1[4]==[9,10]||conf1[4]==[4,8]||conf1[4]==[10,9]||conf1[4]==[8,4])
            lines!(ax,[Point2f0(points[pt]) for pt in [9,10]], color=RGBA{Float64}(0.35, 0.35, 0.35, t), linewidth=lw)
            xvals=0:0.1:4
            yvals = [1 + 0.0123563*x - 0.630187*x^2 + 0.622464*x^3 - 0.194037*x^4 + 0.0194037*x^5 + 3.66123*10^-17*x^6 for x in xvals]
            lines!(ax,[Point2f0([xvals[i],yvals[i]]) for i in 1:length(yvals)], color=RGBA{Float64}(fourcolors[3][1], fourcolors[3][2], fourcolors[3][3], t), linewidth=lw)
        elseif (conf1[3]==[8,9]||conf1[3]==[4,10]||conf1[3]==[9,8]||conf1[3]==[10,4])&&(conf1[4]==[8,9]||conf1[4]==[4,10]||conf1[4]==[9,8]||conf1[4]==[10,4])
            xvals1=0:0.1:3
            xvals2=1:0.1:4
            yvals1=[1 + 0.00765306 *x - 0.389031 *x^2 + 0.320153 *x^3 - 0.0637755 *x^4 for x in xvals1]
            yvals2=yvals1[end:-1:1]
            lines!(ax,[Point2f0([xvals1[i],yvals1[i]]) for i in 1:length(yvals1)], color=RGBA{Float64}(fourcolors[3][1], fourcolors[3][2], fourcolors[3][3], t), linewidth=lw)
            lines!(ax,[Point2f0([xvals2[i],yvals2[i]]) for i in 1:length(yvals2)], color=RGBA{Float64}(0.35, 0.35, 0.35, t), linewidth=lw)
        elseif conf1[3]==[8,9]||conf1[3]==[9,8]||conf1[4]==[8,9]||conf1[4]==[9,8]
            xvals1=0:0.1:3
            yvals1=[1 + 0.00765306 *x - 0.389031 *x^2 + 0.320153 *x^3 - 0.0637755 *x^4 for x in xvals1]
            lines!(ax,[Point2f0([xvals1[i],yvals1[i]]) for i in 1:length(yvals1)], color=RGBA{Float64}(fourcolors[3][1], fourcolors[3][2], fourcolors[3][3], t), linewidth=lw)
            twoconfig = conf1[3]==[8,9]||conf1[3]==[9,8] ? conf1[4] : conf1[3]
            lines!(ax,[Point2f0(points[pt]) for pt in twoconfig], color=RGBA{Float64}(0.35, 0.35, 0.35, t), linewidth=lw)
        elseif conf1[3]==[4,10]||conf1[3]==[10,4]||conf1[4]==[4,10]||conf1[4]==[10,4]
            xvals1=0:0.1:3
            xvals2=1:0.1:4
            yvals2=([1 + 0.00765306 *x - 0.389031 *x^2 + 0.320153 *x^3 - 0.0637755 *x^4 for x in xvals1])[end:-1:1]
            lines!(ax,[Point2f0([xvals2[i],yvals2[i]]) for i in 1:length(yvals2)], color=RGBA{Float64}(0.35, 0.35, 0.35, t), linewidth=lw)
            twoconfig = conf1[3]==[4,10]||conf1[3]==[10,4] ? conf1[4] : conf1[3]
            lines!(ax,[Point2f0(points[pt]) for pt in twoconfig], color=RGBA{Float64}(fourcolors[3][1], fourcolors[3][2], fourcolors[3][3], t), linewidth=lw)
        elseif (conf1[3]==[4,8]||conf1[3]==[8,4]||conf1[4]==[4,8]||conf1[3]==[8,4])
            xvals=0:0.1:4
            yvals = [1 + 0.0123563*x - 0.630187*x^2 + 0.622464*x^3 - 0.194037*x^4 + 0.0194037*x^5 + 3.66123*10^-17*x^6 for x in xvals]
            lines!(ax,[Point2f0([xvals[i],yvals[i]]) for i in 1:length(yvals)], color=RGBA{Float64}(fourcolors[3][1], fourcolors[3][2], fourcolors[3][3], t), linewidth=lw)
            twoconfig = conf1[3]==[4,8]||conf1[3]==[8,4] ? conf1[4] : conf1[3]
            lines!(ax,[Point2f0(points[pt]) for pt in twoconfig], color=RGBA{Float64}(0.35, 0.35, 0.35, t), linewidth=lw)
        elseif (conf1[3]==[9,10]||conf1[3]==[9,10]||conf1[4]==[9,10]||conf1[3]==[9,10])
            lines!(ax,[Point2f0(points[pt]) for pt in [9,10]], color=RGBA{Float64}(0.35, 0.35, 0.35, t), linewidth=lw)
            twoconfig = conf1[3]==[9,10]||conf1[3]==[10,9] ? conf1[4] : conf1[3]
            lines!(ax,[Point2f0(points[pt]) for pt in twoconfig], color=RGBA{Float64}(fourcolors[3][1], fourcolors[3][2], fourcolors[3][3], t), linewidth=lw)
        else
            lines!(ax,[GPoint2f0(points[pt]) for pt in conf1[3]], color=RGBA{Float64}(fourcolors[3][1], fourcolors[3][2], fourcolors[3][3], t), linewidth=lw)
            lines!(ax,[Point2f0(points[pt]) for pt in conf1[4]], color=RGBA{Float64}(0.35, 0.35, 0.35, t), linewidth=lw)
        end


        if (conf2[3]==[9,10]||conf2[3]==[4,8]||conf2[3]==[10,9]||conf2[3]==[8,4])&&(conf2[4]==[9,10]||conf2[4]==[4,8]||conf2[4]==[10,9]||conf2[4]==[8,4])
            lines!(ax,[Point2f0(points[pt]) for pt in [9,10]], color=RGBA{Float64}(0.35, 0.35, 0.35, 1-t), linewidth=lw)
            xvals=0:0.1:4
            yvals = [1 + 0.0123563*x - 0.630187*x^2 + 0.622464*x^3 - 0.194037*x^4 + 0.0194037*x^5 + 3.66123*10^-17*x^6 for x in xvals]
            lines!(ax,[Point2f0([xvals[i],yvals[i]]) for i in 1:length(yvals)], color=RGBA{Float64}(fourcolors[3][1], fourcolors[3][2], fourcolors[3][3], 1-t), linewidth=lw)
        elseif (conf2[3]==[8,9]||conf2[3]==[4,10]||conf2[3]==[9,8]||conf2[3]==[10,4])&&(conf2[4]==[8,9]||conf2[4]==[4,10]||conf2[4]==[9,8]||conf2[4]==[10,4])
            xvals1=0:0.1:3
            xvals2=1:0.1:4
            yvals1=[1 + 0.00765306 *x - 0.389031 *x^2 + 0.320153 *x^3 - 0.0637755 *x^4 for x in xvals1]
            yvals2=yvals1[end:-1:1]
            lines!(ax,[Point2f0([xvals1[i],yvals1[i]]) for i in 1:length(yvals1)], color=RGBA{Float64}(fourcolors[3][1], fourcolors[3][2], fourcolors[3][3], 1-t), linewidth=lw)
            lines!(ax,[Point2f0([xvals2[i],yvals2[i]]) for i in 1:length(yvals2)], color=RGBA{Float64}(0.35, 0.35, 0.35, 1-t), linewidth=lw)
        elseif conf2[3]==[8,9]||conf2[3]==[9,8]||conf2[4]==[8,9]||conf2[4]==[9,8]
            xvals1=0:0.1:3
            yvals1=[1 + 0.00765306 *x - 0.389031 *x^2 + 0.320153 *x^3 - 0.0637755 *x^4 for x in xvals1]
            lines!(ax,[Point2f0([xvals1[i],yvals1[i]]) for i in 1:length(yvals1)], color=RGBA{Float64}(fourcolors[3][1], fourcolors[3][2], fourcolors[3][3], 1-t), linewidth=lw)
            twoconfig = conf2[3]==[8,9]||conf2[3]==[9,8] ? conf2[4] : conf2[3]
            lines!(ax,[Point2f0(points[pt]) for pt in twoconfig], color=RGBA{Float64}(0.35, 0.35, 0.35, 1-t), linewidth=lw)
        elseif conf2[3]==[4,10]||conf2[3]==[10,4]||conf2[4]==[4,10]||conf2[4]==[10,4]
            xvals1=0:0.1:3
            xvals2=1:0.1:4
            yvals2=([1 + 0.00765306 *x - 0.389031 *x^2 + 0.320153 *x^3 - 0.0637755 *x^4 for x in xvals1])[end:-1:1]
            lines!(ax,[Point2f0([xvals2[i],yvals2[i]]) for i in 1:length(yvals2)], color=RGBA{Float64}(0.35, 0.35, 0.35, 1-t), linewidth=lw)
            twoconfig = conf2[3]==[4,10]||conf2[3]==[10,4] ? conf2[4] : conf2[3]
            lines!(ax,[Point2f0(points[pt]) for pt in twoconfig], color=RGBA{Float64}(fourcolors[3][1], fourcolors[3][2], fourcolors[3][3], 1-t), linewidth=lw)
        elseif (conf2[3]==[4,8]||conf2[3]==[8,4]||conf2[4]==[4,8]||conf2[3]==[8,4])
            xvals=0:0.1:4
            yvals = [1 + 0.0123563*x - 0.630187*x^2 + 0.622464*x^3 - 0.194037*x^4 + 0.0194037*x^5 + 3.66123*10^-17*x^6 for x in xvals]
            lines!(ax,[Point2f0([xvals[i],yvals[i]]) for i in 1:length(yvals)], color=RGBA{Float64}(fourcolors[3][1], fourcolors[3][2], fourcolors[3][3], 1-t), linewidth=lw)
            twoconfig = conf2[3]==[4,8]||conf2[3]==[8,4] ? conf2[4] : conf2[3]
            lines!(ax,[Point2f0(points[pt]) for pt in twoconfig], color=RGBA{Float64}(0.35, 0.35, 0.35, 1-t), linewidth=lw)
        elseif (conf2[3]==[9,10]||conf2[3]==[9,10]||conf2[4]==[9,10]||conf2[3]==[9,10])
            lines!(ax,[Point2f0(points[pt]) for pt in [9,10]], color=RGBA{Float64}(0.35, 0.35, 0.35, 1-t), linewidth=lw)
            twoconfig = conf2[3]==[9,10]||conf2[3]==[10,9] ? conf2[4] : conf2[3]
            lines!(ax,[Point2f0(points[pt]) for pt in twoconfig], color=RGBA{Float64}(fourcolors[3][1], fourcolors[3][2], fourcolors[3][3], 1-t), linewidth=lw)
        else
            lines!(ax,[GPoint2f0(points[pt]) for pt in conf2[3]], color=RGBA{Float64}(fourcolors[3][1], fourcolors[3][2], fourcolors[3][3], 1-t), linewidth=lw)
            lines!(ax,[Point2f0(points[pt]) for pt in conf2[4]], color=RGBA{Float64}(0.35, 0.35, 0.35, 1-t), linewidth=lw)
        end
        scatter!(ax,[Point2f0([2,1])]; color=:red2,markersize=60)
        scatter!(ax,[Point2f0(pt) for pt in points]; color=:black,markersize=60)
        save("../images/homotopytriang$(t).png",fig)
    end
end


function printGraphs(; prefix="linearweight", suffix="4,9")
    fig = Figure(size=(1400,500))
    ax = Axis(fig[1,1])
    hidedecorations!(ax)
    hidespines!(ax)
    linestyles = [:solid,  :solid, :dash, :dash]
    global count=1
    for boxsize in ["0.1", "100.0", "1.0","10.0",]
        ourmodel = Dict()
        try
            f = open("../data/$(prefix)$(suffix)storedsolutions$(boxsize).txt", "r")
    
            global pointnumber = parse(Int,readline(f))
            while ! eof(f)  
                sstring = split(readline(f), ": ")
                keystring = split(sstring[1], "; ")
                key = (parse(Int,keystring[1][2]), (parse(Int,keystring[1][5])))
                weight = parse(Float64,keystring[2])
                ourmodel[weight] = [parse(Int,entry) for entry in split(sstring[2][2:end-1], ", ")]
             end
            close(f)
        catch e
            display(e)
            continue
        end

        points = [ourmodel[w][1] for w in sort(Float64.(keys(ourmodel)), rev=true)]
        lines!(ax, 1:length(points), points ./ pointnumber; linewidth=count>=3 ? 9-(2*count-5) : ((count==2) ? 6 : 9), linestyle=linestyles[count], color=cgrad(:Dark2_4)[count], label = "$(boxsize)")
        scatter!(ax,[0.01,0.01,0.01],[0.978,0.9785,0.979], markersize=5)
        count+=1
    end
    #axislegend(ax, merge = true, unique = true, position = :rt, labelsize=26)
    save("../images/$(prefix)$(suffix).png", fig)
end


function plotWeightedCovers(; boxsize=1, prefix="TWOBEST", suffix="10,12,15")
    helperDict1, ourmodel1 = Dict(), Dict()
    helperDict2, ourmodel2 = Dict(), Dict()
    helperDict3, ourmodel3 = Dict(), Dict()
    helperDict4, ourmodel4 = Dict(), Dict()

    try
        f = open("../data/triangweight4,10,15storedsolutions1.txt", "r")
        g = open("../data/triangweightNEXTTRY10,12,15storedsolutions1.txt", "r")
        h4915 = open("../data/triangweight4,9,15storedsolutions1.txt", "r")
        k91215 = open("../data/triangweightNEXTTRY9,12,15storedsolutions1.txt", "r")

        global pointnumber1 = parse(Int,readline(f))
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
            helperDict1[(key, weight)] = [parse(Int,entry) for entry in split(sstring[2][2:end-1], ", ")]
         end
        close(f)
        global pointnumber2 = parse(Int,readline(g))
        while ! eof(g)  
            str = readline(g)
            display(str)
            if str==""
                break
            end
            sstring = split(str, ": ")
            keystring = split(sstring[1], "; ")
            key = (parse(Int,keystring[1][2]), (parse(Int,keystring[1][5])), (parse(Int,keystring[1][8])))
            weight = parse(Float64,keystring[2])
            helperDict2[(key, weight)] = [parse(Int,entry) for entry in split(sstring[2][2:end-1], ", ")]
         end
        close(g)
        global pointnumber3 = parse(Int,readline(h4915))
        while ! eof(h4915)  
            str = readline(h4915)
            display(str)
            if str==""
                break
            end
            sstring = split(str, ": ")
            keystring = split(sstring[1], "; ")
            key = (parse(Int,keystring[1][2]), (parse(Int,keystring[1][5])), (parse(Int,keystring[1][8])))
            weight = parse(Float64,keystring[2])
            helperDict3[(key, weight)] = [parse(Int,entry) for entry in split(sstring[2][2:end-1], ", ")]
         end
        close(h4915)
        global pointnumber4 = parse(Int,readline(k91215))
        while ! eof(k91215)  
            str = readline(k91215)
            display(str)
            if str==""
                break
            end
            sstring = split(str, ": ")
            keystring = split(sstring[1], "; ")
            key = (parse(Int,keystring[1][2]), (parse(Int,keystring[1][5])), (parse(Int,keystring[1][8])))
            weight = parse(Float64,keystring[2])
            helperDict4[(key, weight)] = [parse(Int,entry) for entry in split(sstring[2][2:end-1], ", ")]
         end
        close(k91215)
    catch
        throw(error("No Data recorded!"))
    end


    for key in keys(helperDict1)
        if key[1] in keys(ourmodel1)
            ourmodel1[key[1]][key[2]] = helperDict1[key]
        else
            ourmodel1[key[1]] = Dict(key[2]=>helperDict1[key])
        end
    end
    for key in keys(helperDict2)
        if key[1] in keys(ourmodel2)
            ourmodel2[key[1]][key[2]] = helperDict2[key]
        else
            ourmodel2[key[1]] = Dict(key[2]=>helperDict2[key])
        end
    end
    for key in keys(helperDict3)
        if key[1] in keys(ourmodel3)
            ourmodel3[key[1]][key[2]] = helperDict3[key]
        else
            ourmodel3[key[1]] = Dict(key[2]=>helperDict3[key])
        end
    end
    for key in keys(helperDict4)
        if key[1] in keys(ourmodel4)
            ourmodel4[key[1]][key[2]] = helperDict4[key]
        else
            ourmodel4[key[1]] = Dict(key[2]=>helperDict4[key])
        end
    end


    global maxval1, minval1 = maximum(vcat([vcat(values(ourmodel1[key])...) for key in keys(ourmodel1)]...)) / pointnumber1, minimum(vcat([vcat(values(ourmodel1[key])...) for key in keys(ourmodel1)]...)) / pointnumber1
    global maxval2, minval2 = maximum(vcat([vcat(values(ourmodel2[key])...) for key in keys(ourmodel2)]...)) / pointnumber2, minimum(vcat([vcat(values(ourmodel2[key])...) for key in keys(ourmodel2)]...)) / pointnumber2
    global maxval3, minval3 = maximum(vcat([vcat(values(ourmodel3[key])...) for key in keys(ourmodel3)]...)) / pointnumber3, minimum(vcat([vcat(values(ourmodel3[key])...) for key in keys(ourmodel3)]...)) / pointnumber3
    global maxval4, minval4 = maximum(vcat([vcat(values(ourmodel4[key])...) for key in keys(ourmodel4)]...)) / pointnumber4, minimum(vcat([vcat(values(ourmodel4[key])...) for key in keys(ourmodel4)]...)) / pointnumber4

    global maxval, minval = maximum([maxval1, maxval2, maxval3, maxval4]), minimum([minval1, minval2, minval3, minval4])
    fig = Figure(size=(1400,1100), fontsize=30)

    
    for key in keys(ourmodel1)
        ax=Axis(fig[1,1], aspect = 1)
        heatMatrix1 = Matrix{Float64}(undef,maximum(length.(values(ourmodel1[key]))),maximum(length.(values(ourmodel1[key])))); 
        heatMatrix1 .= NaN

        global row = 1
        for weight in sort(collect(keys(ourmodel1[key])))
            heatMatrix1[row, 1:length(ourmodel1[key][weight])] = ourmodel1[key][weight] ./ pointnumber1
            row += 1
        end
        heatMatrix1 .= heatMatrix1[:, :]

        hm = heatmap!(ax, heatMatrix1; colormap=:viridis, colorrange = (minval1,maxval1))#, colorrange=(minval/pointnumber, maxval/pointnumber))
        hidespines!(ax)
        hidedecorations!(ax)
        xlims!(ax, (-0.01,18))
        ylims!(ax, (-0.25,18.25))
        #text!(ax, [0.2,0.2,17.3], [-0.25,17.55,-0.25]; text=["10", "15", "4"], fontsize=38)


        pointarray1 = []
        for weight in keys(ourmodel1[key])
            for i in 1:length(ourmodel1[key][weight])
                push!(pointarray1, [weight, (i-1)/(maximum(length.(values(ourmodel1[key])))-1), ourmodel1[key][weight][i] / pointnumber1])
            end
        end
    end

    for key in keys(ourmodel2)
        ax=Axis(fig[1,2], aspect = 1)
        heatMatrix2 = Matrix{Float64}(undef,maximum(length.(values(ourmodel2[key]))),maximum(length.(values(ourmodel2[key])))); 
        heatMatrix2 .= NaN

        global row = 1
        for weight in sort(collect(keys(ourmodel2[key])))
            heatMatrix2[row, 1:length(ourmodel2[key][weight])] = ourmodel2[key][weight] ./ pointnumber2
            row += 1
        end
        heatMatrix2 .= heatMatrix2[:, :]
    
        global hm = heatmap!(ax, heatMatrix2; colormap=:viridis, colorrange = (minval,maxval))#, colorrange=(minval/pointnumber, maxval/pointnumber))
        hidespines!(ax)
        hidedecorations!(ax)
        xlims!(ax, (-0.01,18))
        ylims!(ax, (-0.25,18.25))
        #text!(ax, [0.2,0.2,17.3], [-0.25,17.55,-0.25]; text=["12", "15", "10"], fontsize=38)

        pointarray2 = []
        for weight in keys(ourmodel2[key])
            for i in 1:length(ourmodel2[key][weight])
                push!(pointarray2, [weight, (i-1)/(maximum(length.(values(ourmodel2[key])))-1), ourmodel2[key][weight][i] / pointnumber2])
            end
        end
    end

    for key in keys(ourmodel3)
        ax=Axis(fig[2,1], aspect = 1)
        heatMatrix3 = Matrix{Float64}(undef,maximum(length.(values(ourmodel3[key]))),maximum(length.(values(ourmodel3[key])))); 
        heatMatrix3 .= NaN

        global row = 1
        for weight in sort(collect(keys(ourmodel3[key])))
            heatMatrix3[row, 1:length(ourmodel1[key][weight])] = ourmodel3[key][weight] ./ pointnumber3
            row += 1
        end
        heatMatrix3 .= heatMatrix3[:, :]
    
        global hm = heatmap!(ax, heatMatrix3; colormap=:viridis, colorrange = (minval3,maxval3))#, colorrange=(minval/pointnumber, maxval/pointnumber))
        hidespines!(ax)
        hidedecorations!(ax)
        xlims!(ax, (-0.01,18))
        ylims!(ax, (-0.25,18.25))
        #text!(ax, [0.2,0.2,17.3], [-0.25,17.55,-0.25]; text=["4", "9", "15"], fontsize=38)

        pointarray3 = []
        for weight in keys(ourmodel3[key])
            for i in 1:length(ourmodel3[key][weight])
                push!(pointarray3, [weight, (i-1)/(maximum(length.(values(ourmodel3[key])))-1), ourmodel3[key][weight][i] / pointnumber3])
            end
        end
    end

    for key in keys(ourmodel4)
        ax=Axis(fig[2,2], aspect = 1)
        heatMatrix4 = Matrix{Float64}(undef,maximum(length.(values(ourmodel4[key]))),maximum(length.(values(ourmodel4[key])))); 
        heatMatrix4 .= NaN

        global row = 1
        for weight in sort(collect(keys(ourmodel4[key])))
            heatMatrix4[row, 1:length(ourmodel4[key][weight])] = ourmodel4[key][weight] ./ pointnumber4
            row += 1
        end
        heatMatrix4 .= heatMatrix4[:, :]
    
        global hm = heatmap!(ax, heatMatrix4; colormap=:viridis, colorrange = (minval,maxval4))#, colorrange=(minval/pointnumber, maxval/pointnumber))
        hidespines!(ax)
        hidedecorations!(ax)
        xlims!(ax, (-0.01,18))
        ylims!(ax, (-0.25,18.25))
        #text!(ax, [0.2,0.2,17.3], [-0.25,17.55,-0.25]; text=["4", "9", "15"], fontsize=38)

        pointarray4 = []
        for weight in keys(ourmodel4[key])
            for i in 1:length(ourmodel4[key][weight])
                push!(pointarray4, [weight, (i-1)/(maximum(length.(values(ourmodel4[key])))-1), ourmodel4[key][weight][i] / pointnumber4])
            end
        end
    end

    Colorbar(fig[:, end+1], colorrange = (minval,maxval); size=30)
    save("../images/discretizedHeatmapnew.png",fig)
end


#=
Samples the region of the reaction rate constants where a(η)>0 and b(η)<0 to find which points are captured
by the sufficient condition corresponding to the suggested cover with circuit number `θsuggestion` compared to
a different cover with circuit number `θbaseline`.
=#
function empiricalComparisonOfTwoCovers(θsuggestion, θbaseline, K, κ, aη, bη, mcoef; numberOfSamplingRuns=100, boxsizes=[5 for _ in 1:8], cover_1=15, cover_2 = 9)
    vector_baseline_wins = []
    vector_our_wins = []
    vector_both_wins = []
    vector_no_wins = []
    for sampleindex in 1:numberOfSamplingRuns
        display("Run: $(sampleindex)")
        global sampling = filter(sampler -> !any(t->isapprox(t,0), sampler) && evaluate(aη,vcat(K,[κ[3],κ[6],κ[9],κ[12]])=>sampler)>0 && evaluate(bη,vcat(K,[κ[3],κ[6],κ[9],κ[12]])=>sampler)<0, [boxsizes .* abs.(rand(Float64,8)) for _ in 1:1000000])
        for ind in 1:length(sampling)
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
    ax = [Axis(fig[1,1]; xlabel = "K_1", ylabel="K_2"), Axis(fig[1,2]; xlabel = "K_3", ylabel="K_4"), Axis(fig[2,1]; xlabel = "kappa_3", ylabel="kappa_6"), Axis(fig[2,2]; xlabel = "kappa_9", ylabel="kappa_{12}")]
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

function plottriangle()
    pointsForPlot = [(0,1),(2,0),(4,2)]
    fig = Figure(size=(1200,1200))
    ax = Axis(fig[1,1], aspect=1)
    hidespines!(ax)
    hidedecorations!(ax)
    poly!(ax,[Point2f0(pt) for pt in pointsForPlot]; color=RGBA{Float64}(0.1, 0.1, 0.1, 0.04), strokewidth=0)
    lines!(ax,[Point2f0(pt) for pt in pointsForPlot[vcat(1:3,1)]], color=:black, linewidth=10)
    scatter!(ax,[Point2f0([2,1])]; color=:red2,markersize=60)
    text!(ax,[Point2f0([1.675,0.75])], text="(2,1)"; color=:red2,fontsize=70)
    scatter!(ax,[Point2f0(pt) for pt in pointsForPlot]; color=:black, markersize=60)
    text!(ax,[Point2f0([-0.32,1.09]), Point2f0([1.665,-0.25]), Point2f0([3.675,2.065])], text=["(0,1)", "(2,0)", "(4,2)"]; color=:black,fontsize=70)
    xlims!(ax,(-0.3,4.3))
    ylims!(ax,(-0.3,2.3))
    save("../images/triangexample.png",fig)
end

function plotTwotriangles()
    pointsForPlot = [(0,0),(2,0),(4,2),(3,2),(0,1),(2,1)]
    fig = Figure(size=(1200,1200))
    ax = Axis(fig[1,1], aspect=1)
    hidespines!(ax)
    hidedecorations!(ax)
    poly!(ax,[Point2f0(pt) for pt in pointsForPlot[1:5]]; color=RGBA{Float64}(0.1, 0.1, 0.1, 0.04), strokewidth=0)

    fourcolors = map(col -> (red(col), green(col), blue(col)), distinguishable_colors(3, [RGB(1,1,1), RGB(0,0,0)], dropseed=true, lchoices = range(25, stop=50, length=15), hchoices = range(120, stop=350, length=20)))
    triangles = [[1,2,4],[2,3,5]]
    poly!(ax,[Point2f0(pointsForPlot[pt]) for pt in triangles[1]]; color=RGBA{Float64}(fourcolors[1][1], fourcolors[1][2], fourcolors[1][3], 0.15),strokewidth=0)
    lines!(ax,[Point2f0(pointsForPlot[pt]) for pt in vcat(triangles[1],triangles[1][1])], color=RGBA{Float64}(fourcolors[1][1], fourcolors[1][2], fourcolors[1][3], 1), linewidth=5)
    poly!(ax,[Point2f0(pointsForPlot[pt]) for pt in triangles[2]]; color=RGBA{Float64}(fourcolors[2][1], fourcolors[2][2], fourcolors[2][3], 0.15),strokewidth=0)
    lines!(ax,[Point2f0(pointsForPlot[pt]) for pt in vcat(triangles[2],triangles[2][1])], color=RGBA{Float64}(fourcolors[2][1], fourcolors[2][2], fourcolors[2][3], 1), linewidth=5)
    lines!(ax,[Point2f0(pointsForPlot[pt]) for pt in vcat(1:5,1)], color=:black, linewidth=10)

    text!(ax,[Point2f0([1.675,0.75])], text="(2,1)"; color=:red2,fontsize=70)
    scatter!(ax,[Point2f0(pt) for pt in pointsForPlot]; color=:black, markersize=60)
    scatter!(ax,[Point2f0([2,1])]; color=:red2,markersize=60)
    scatter!(ax,[Point2f0([2,0])]; color=RGBA(0.2,0.2,0.75,1),markersize=60)
    text!(ax,[Point2f0([-0.32,1.12]), Point2f0([1.665,-0.25]), Point2f0([3.675,2.065]), Point2f0([2.675,2.065]), Point2f0([-0.32,-0.25])], text=["(0,1)", "(2,0)", "(4,2)", "(3,2)", "(0,0)"]; color=:black,fontsize=70)
    text!(ax,[Point2f0([1.665,-0.25])], text=["(2,0)"]; color=RGBA(0.2,0.2,0.75,1),fontsize=70)
    xlims!(ax,(-0.3,4.3))
    ylims!(ax,(-0.3,2.3))


    save("../images/twotriangexample2.png",fig)
end


function plotOnetrianglesOneLine()
    pointsForPlot = [(0,0),(2,0),(4,2),(0,1),(2,1)]
    fig = Figure(size=(1200,1200))
    ax = Axis(fig[1,1], aspect=1)
    hidespines!(ax)
    hidedecorations!(ax)
    poly!(ax,[Point2f0(pt) for pt in pointsForPlot[1:4]]; color=RGBA{Float64}(0.1, 0.1, 0.1, 0.04), strokewidth=0)

    fourcolors = map(col -> (red(col), green(col), blue(col)), distinguishable_colors(3, [RGB(1,1,1), RGB(0,0,0)], dropseed=true, lchoices = range(25, stop=50, length=15), hchoices = range(120, stop=350, length=20)))
    triangles = [[2,3,4]]
    poly!(ax,[Point2f0(pointsForPlot[pt]) for pt in triangles[1]]; color=RGBA{Float64}(fourcolors[1][1], fourcolors[1][2], fourcolors[1][3], 0.15),strokewidth=0)
    lines!(ax,[Point2f0(pointsForPlot[pt]) for pt in vcat(triangles[1],triangles[1][1])], color=RGBA{Float64}(fourcolors[1][1], fourcolors[1][2], fourcolors[1][3], 1), linewidth=5)
    lines = [[1,3]]
    lw = 15
    lines!(ax,[Point2f0(pointsForPlot[pt]) for pt in lines[1]], color=RGBA{Float64}(fourcolors[2][1], fourcolors[2][2], fourcolors[2][3], 1), linewidth=lw)
    lines!(ax,[Point2f0(pointsForPlot[pt]) for pt in vcat(1:4,1)], color=:black, linewidth=10)

    text!(ax,[Point2f0([2.175,0.9])], text="(2,1)"; color=:red2,fontsize=70)
    scatter!(ax,[Point2f0(pt) for pt in pointsForPlot]; color=:black, markersize=60)
    scatter!(ax,[Point2f0([2,1])]; color=:red2,markersize=60)
    scatter!(ax,[Point2f0([4,2])]; color=RGBA(0.2,0.2,0.75,1),markersize=60)
    text!(ax,[Point2f0([-0.32,1.12]), Point2f0([1.665,-0.25]), Point2f0([1.665,-0.25]), Point2f0([-0.32,-0.25])], text=["(0,1)", "(2,0)", "(2,0)",  "(0,0)"]; color=:black,fontsize=70)
    text!(ax,[Point2f0([3.675,2.065])], text=["(4,2)"]; color=RGBA(0.2,0.2,0.75,1),fontsize=70)
    xlims!(ax,(-0.3,4.3))
    ylims!(ax,(-0.3,2.3))
    save("../images/onetriangonelineexample.png",fig)
end


end