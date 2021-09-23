# function circleshape(h, k, r)
#     theta = LinRange(0, 2*pi, 500)
#     return h .+ r*sin.(theta), k.+ r*cos.(theta)
# end

# function plotlayout!(p, turbinex, turbiney, rotordiameter; linecolor=:black, fillcolor=:blue, markeralpha=1, title="")
#     nturbines = length(turbinex)
#     for i in 1:nturbines
#         plot!(p, circleshape(turbinex[i],turbiney[i],rotordiameter[i]/2.0), seriestype=[:shape],
#             linecolor=linecolor, c=fillcolor, legend=false, aspect_ratio=1, alpha=markeralpha, title=title)
#     end
#     display(p)
# end
"""
    plotwindfarm!(ax, boundary_vertices, turbinex, turbiney, rotordiameter; aspect="equal", xlim=[], ylim=[], fill=false, color="k", markeralpha=1, title="")

    Convenience function for plotting wind farms

# Arguments
- `ax
- `boundary_vertices::Array{Float,1}(nvertices)`: an nx2 array of boundary vertices
- `turbinex::Array{Float,1}(nturbines)`: an array x coordinates of wind turbine locations
- `turbiney::Array{Float,1}(nturbines)`: an array y coordinates of wind turbine locations
- `rotordiameter::Array{Float,1}(nturbines)`: an array rotor diameters of wind turbines
- `aspect::String`: set plot aspect ratio, default="equal"
- `xlim::Array`: limits in x coordinate. "[]" results in limits being automatically defined
- `ylim::Array`: limits in y coordinate. "[]" results in limits being automatically defined
- `fill::Bool`: determines whether turbine circle markers are filled or not
- `color::=String`: sets color for turbine markers
- `markeralpha::Int`: determines tranparancy of turbine markers
- `itle::String`: optional title to include on the plot
"""
function plotwindfarm!(ax, boundary_vertices, turbinex, turbiney, rotordiameter; 
    aspect="equal", xlim=[], ylim=[], fill=false, turbinecolor="k", boundarycolor="k", 
    boundarylinestyle="-", turbinelinestyle="-", markeralpha=1, title="")

    plotlayout!(ax, turbinex, turbiney, rotordiameter; aspect=aspect, xlim=xlim, ylim=ylim, fill=fill, color=turbinecolor, markeralpha=markeralpha, title=title, linestyle=turbinelinestyle)
    
    plotboundary!(ax, boundary_vertices; color=boundarycolor, linestyle=boundarylinestyle)

end

"""
    plotlayout!(ax, turbinex, turbiney, rotordiameter; aspect="equal", xlim=[], ylim=[], fill=false, color="k", markeralpha=1, title="")

    Convenience function for plotting wind farm layouts

# Arguments
- `ax
- `turbinex::Array{Float,1}(nturbines)`: an array x coordinates of wind turbine locations
- `turbiney::Array{Float,1}(nturbines)`: an array y coordinates of wind turbine locations
- `rotordiameter::Array{Float,1}(nturbines)`: an array rotor diameters of wind turbines
- `aspect::String`: set plot aspect ratio, default="equal"
- `xlim::Array`: limits in x coordinate. "[]" results in limits being automatically defined
- `ylim::Array`: limits in y coordinate. "[]" results in limits being automatically defined
- `fill::Bool`: determines whether turbine circle markers are filled or not
- `color::=String`: sets color for turbine markers
- `markeralpha::Int`: determines tranparancy of turbine markers
- `itle::String`: optional title to include on the plot
"""
function plotlayout!(ax, turbinex, turbiney, rotordiameter; aspect="equal", xlim=[], ylim=[], fill=false, color="k", markeralpha=1, title="", linestyle="-")
    nturbines = length(turbinex)
    if xlim == []
        xlim = [minimum(turbinex)-sum(rotordiameter)/nturbines, maximum(turbinex)+sum(rotordiameter)/nturbines]
    end
    if ylim == []
        ylim = [minimum(turbiney)-sum(rotordiameter)/nturbines, maximum(turbiney)+sum(rotordiameter)/nturbines]
    end
    # add turbines
    for i in 1:nturbines
        circle = matplotlib.patches.Circle((turbinex[i], turbiney[i]), rotordiameter[i]/2.0, fill=fill, color=color, linestyle=linestyle)
        ax.add_patch(circle)
    end
    ax.set(xlim=xlim, ylim=ylim, aspect=aspect)

end

"""
    plotboundary!(ax, boundary_vertices; aspect="equal", xlim=[], ylim=[], fill=false, color="k", markeralpha=1, title="")

    Convenience function for plotting wind farm boundaries

# Arguments
- `ax
- `boundary_vertices::Array{Float,1}(nvertices)`: an nx2 array of boundary vertices
- `aspect::String`: set plot aspect ratio, default="equal"
- `xlim::Array`: limits in x coordinate. "[]" results in limits being automatically defined
- `ylim::Array`: limits in y coordinate. "[]" results in limits being automatically defined
- `fill::Bool`: determines whether turbine circle markers are filled or not
- `color::=String`: sets color for turbine markers
- `markeralpha::Int`: determines tranparancy of turbine markers
- `itle::String`: optional title to include on the plot
"""
function plotboundary!(ax, boundary_vertices; color="k", linestyle="-")
    
    # add boundary
    x = push!(boundary_vertices[:,1], boundary_vertices[1,1])
    y = push!(boundary_vertices[:,2], boundary_vertices[1,2])
    println(x, y)
    ax.plot(x, y, color=color, linestyle=linestyle)
    
end

"""
    plotrotorsamplepoints!(ax, y, z; rotordiameter=2.0, aspect="equal", ylim=[], zlim=[], fill=false, color="k", markeralpha=1, title="")

    Convenience function for plotting where points are being sampled on the wind turbine rotor

# Arguments
- `ax
- `y::Array{Float,1}(nturbines)`: an array x coordinates of wind turbine locations
- `z::Array{Float,1}(nturbines)`: an array y coordinates of wind turbine locations
- `rotordiameter::Number`: rotor diameter of wind turbine
- `aspect::String`: set plot aspect ratio, default="equal"
- `ylim::Array`: limits in y coordinate. "[]" results in limits being automatically defined
- `zlim::Array`: limits in z coordinate. "[]" results in limits being automatically defined
- `fill::Bool`: determines whether turbine circle markers are filled or not
- `color::=String`: sets color for turbine markers
- `markeralpha::Int`: determines tranparancy of turbine markers between 0 (transparent) and 1 (opaque)
- `itle::String`: optional title to include on the plot
"""
function plotrotorsamplepoints!(ax, y, z; rotordiameter=2.0, aspect="equal", ylim=[], zlim=[], fill=false, color="k", markeralpha=1, title="")

    # set plot limits
    if ylim == []
        ylim = [-1.05, 1.05].*rotordiameter/2.0
    end
    if zlim == []
        zlim = [-1.05, 1.05].*rotordiameter/2.0
    end

    # add points
    ax.scatter(y, z, color=color, markeralpha)

    # add rotor swept area
    circle = matplotlib.patches.Circle((0.0, 0.0), rotordiameter/2.0, fill=fill, color=color)
    ax.add_patch(circle)
    println(ylim, zlim)
    ax.set(xlim=ylim, ylim=zlim, aspect=aspect, title=title)

end

"""
    plotwindresource!(ax::Array, windresource::ff.DiscretizedWindResource; roundingdigits=[1,3], fill=false, alpha=0.5, colors=["b", "b"], fontsize=8, edgecolor=nothing, rlabel_position=-45)
    
    Convenience function for visualizing the wind speed and wind frequency roses 

# Arguments
- `ax::Array`: pre-initialized single dimension array of axes from pyplot with length at least 2
- `windresource::ff.DiscretizedWindResource`: wind rose information
- `roundingdigits::Array{Int, 1}`: how many significant digits to round to on each axis in ax
- `fill::Bool`: determines whether bars are filled or not
- `alpha::Number`: tranparancy of bars between 0 (transparent) and 1 (opaque)
- `colors::=Array{String, 1}`: sets color for turbine markers
- `fontsize::Int`: font size of text on figures
- `edgecolor`: color of edges of each bar in polar chart, nothing means no color
- `rlabel_position:Number`: Angle at which to draw the radial axes
"""
function plotwindresource!(ax::Array, windresource::ff.DiscretizedWindResource; roundingdigits=[1,3], fill=false, alpha=0.5, colors=["b", "b"], fontsize=8, edgecolor=nothing, rlabel_position=-45, titles=["Wind Speed", "Wind Probability"])
    
    # extract wind resource elements for windrose plots
    d = windresource.wind_directions
    s = windresource.wind_speeds
    f = windresource.wind_probabilities

    # plot wind speed rose
    plotwindrose!(ax[1], d, s, roundingdigit=roundingdigits[1], color=colors[1], alpha=alpha, fontsize=fontsize, edgecolor=edgecolor, units="m/s", title=titles[1])
    
    # plot wind frequency rose
    plotwindrose!(ax[2], d, f, roundingdigit=roundingdigits[2], color=colors[2], alpha=alpha, fontsize=fontsize, edgecolor=edgecolor, units="%", title=titles[2])

end

"""
    plotwindrose!(ax, d, f; roundingdigit=1, color="C0",alpha=0.5,fontsize=8,
    dticks=(0,pi/4,pi/2,3*pi/4,pi,5*pi/4,3*pi/2,7*pi/4),
    dlabels=("E","NE","N","NW","W","SW","S","SW"),
    fticks=nothing, flabels=nothing, normalize=false, edgecolor=nothing, units="",
    rlabel_position=-45)

    Convenience function for creating a windrose from any polar data

# Arguments
- `ax::PyCall.PyObject`: pre-initialized axis from pyplot
- `d::Vector`: wind rose directions
- `f::Vector`: wind rose radial variable
- `roundingdigit::Int`: how many significant digits to round to
- `color::=String`: sets color for turbine markers
- `alpha::Number`: tranparancy of bars between 0 (transparent) and 1 (opaque)
- `fontsize::Int`: font size of text on figures
- `dticks::Tuple`: contains angular tick locations
- `dlabels::Tuple`: contains angular tick labels
- `fticks::Tuple`: contains radial tick locations
- `flabels::Tuple`: contains radial tick labels
- `normalize::Bool`: choose whether or not to normalize by the sum of f
- `edgecolor`: color of edges of each bar in polar chart, nothing means no color
- `units::String`: Units to append to flabels
- `rlabel_position:Number`: Angle at which to draw the radial axis
"""
function plotwindrose!(ax, d, f; roundingdigit=1, color="C0",alpha=0.5,fontsize=8,
    dticks=(0,pi/4,pi/2,3*pi/4,pi,5*pi/4,3*pi/2,7*pi/4),
    dlabels=("E","NE","N","NW","W","SW","S","SW"),
    fticks=nothing, flabels=nothing, normalize=false, edgecolor=nothing, units="",
    rlabel_position=-45, title="")

    # set up function ticks if not provided, scale and round appropriately
    if fticks === nothing
        fmin = floor(minimum(f), digits=roundingdigit)
        fmax = ceil(maximum(f), digits=roundingdigit)
        frange = fmax - fmin 
        if frange == 0.0
            fticks = round.(collect((0:0.25:1.25).*fmax)[2:end-1], digits=roundingdigit)
        else
            fticks = round.(collect(fmin:frange/4.0:fmax)[2:end-1], digits=roundingdigit)
        end
    end
    # set up function labels if not provided
    if flabels === nothing 
        flabels = string.(fticks).*units
    end
    
    # normalize if desired
    if normalize
        f = f./sum(f)
    end

    # adjust to radians if directions provided in degrees
    if maximum(d) > 10
        d = deg2rad.(d)
    end

    # get the number of wind directions
    ndirs = length(d)

    # specify bar width
    width = (2*pi/ndirs)*0.95
   
    # plot wind rose
    ax.bar(pi/2 .-d,f,width=width,color=color,alpha=alpha,edgecolor=edgecolor)

    # format polar plot
    ax.set_xticks(dticks)
    ax.set_xticklabels(dlabels,fontsize=fontsize)
    ax.set_rgrids(fticks,flabels,angle=rlabel_position,fontsize=fontsize)
    for tick in ax.yaxis.get_majorticklabels()
        tick.set_horizontalalignment("center")
    end
    ax.set_title(title, y=-0.25,fontsize=fontsize)
    
end

"""
    add_turbine!(ax; view="side", hubdiameter=0.1, hubheight=0.9, radius=0.5, chord=0.1, 
        nacellewidth=0.3, nacelleheight=0.1, towerbottomdiam=0.1, towertopdiam=0.05, 
        overhang=0.05, s=5)

    Convenience function for adding wind turbines to plots.

# Arguments
- `ax::PyCall.PyObject`: pre-initialized axis from pyplot
- `view::Number`: determines which turbine view to use "top" or "side" (default)
- `hubdiameter::Number`: hub diameter in axis coordinate frame
- `hubheight::Number`: hub height in axis coordinate frame
- `radius::Number`: full rotor radius in axis coordinate frame
- `chord::Number`: maximum chord in axis coordinate frame
- `nacellewidth::Number`: nacelle width in axis coordinate frame
- `nacelleheight::Number`: nacelle height in axis coordinate frame
- `towerbottomdiam::Number`: tower bottom diameter in axis coordinate frame
- `towertopdiam::Number`: tower top diameter in axis coordinate frame
- `overhang::Number`: overhang (distance from blade attachment to tower bottom in x axis) in axis coordinate frame
- `s::Number`: scales overhang and tower location in x direction to work with condensed x axis as in long contour plots
"""
function add_turbine!(ax; view="side", hubdiameter=0.1, hubheight=0.9, radius=0.5, chord=0.1, nacellewidth=0.3, nacelleheight=0.1, towerbottomdiam=0.1, towertopdiam=0.05, overhang=0.05, s=5, color="k")

    if view == "side"

        # create blade patches
        blade1 = plt.matplotlib.patches.Ellipse((0,hubheight+radius/2),chord, 0.5, color=color)
        blade2 = plt.matplotlib.patches.Ellipse((0,hubheight-radius/2),chord, 0.5, color=color)
        
        # create hub and nacelle patches
        hub = plt.matplotlib.patches.Ellipse((0,hubheight),3*hubdiameter, hubdiameter, color=color)
        nacelle = plt.matplotlib.patches.Rectangle((0,hubheight-hubdiameter/2),nacellewidth, nacelleheight, color=color)
        
        # scale if desired
        towerbottomdiam *= s 
        towertopdiam *= s
        overhang *= s

        # get difference between top and bottom tower diameters
        ddiff = abs(towertopdiam-towerbottomdiam)

        # calculate polygon x points for tower
        p1x = overhang
        p2x = overhang+ddiff/2
        p3x = p2x + towertopdiam
        p4x = p1x + towerbottomdiam

        # create tower patch
        tower = plt.matplotlib.patches.Polygon([[p1x, 0.0],[p2x, hubheight],[p3x, hubheight],[p4x, 0.0]], closed=true, color=color)
        
        # add patches to axis
        ax.add_patch(blade1)
        ax.add_patch(blade2)
        ax.add_patch(hub)
        ax.add_patch(nacelle)
        ax.add_patch(tower)

    elseif view == "top"

        # create blade patches
        blade1 = plt.matplotlib.patches.Ellipse((0,radius/2),chord, 0.5, color=color)
        blade2 = plt.matplotlib.patches.Ellipse((0,radius/2),chord, 0.5, color=color)

        # create hub and nacelle patches
        hub = plt.matplotlib.patches.Ellipse((0.0,0.0),3*hubdiameter, hubdiameter, color=color)
        nacelle = plt.matplotlib.patches.Rectangle((0,-hubdiameter/2),nacellewidth, nacelleheight, color=color)

        # add patches to axis
        ax.add_patch(blade1)
        ax.add_patch(blade2)
        ax.add_patch(hub)
        ax.add_patch(nacelle)
        
    end
end