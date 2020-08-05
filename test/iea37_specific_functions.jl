""" Functions specific to participation and optimization of the IEA37 case studies """

function getCs34NameYAML(indx::Int)
    name = "no such region"

    if(indx == 1)
        name = "IIIa"
    elseif(indx == 2)
        name = "IIIb"
    elseif(indx == 3)
        name = "IVa"
    elseif(indx == 4)
        name = "IVb"
    elseif(indx == 5)
        name = "IVc"
    end
    
    return name
end

function getCs34Name(indx::Int)
    name = "no such region"

    if(indx == 1)
        name = "3a"
    elseif(indx == 2)
        name = "3b"
    elseif(indx == 3)
        name = "4a"
    elseif(indx == 4)
        name = "4b"
    elseif(indx == 5)
        name = "4c"
    end
    
    return name
end

function getCs34NumTurbs(sReg::String)
    numTurbs = zeros(1)

    if(sReg == "cs3")
        numTurbs = 25
    elseif(sReg == "cs4")
        numTurbs = 81
    elseif(sReg == "3a")
        numTurbs = 31
    elseif(sReg == "3b")
        numTurbs = 11
    elseif(sReg == "4a")
        numTurbs = 16
    elseif(sReg == "4b")
        numTurbs = 14
    elseif(sReg == "4c")
        numTurbs = 9
    end
    
    return numTurbs
end

function getCs34NumBndryTurbs(sReg::String)
    numBndryTurbs = zeros(1)
    # As recommended, .45 of the full number
    if(sReg == "cs3")
        numBndryTurbs = 11 #25
    elseif(sReg == "cs4")
        numBndryTurbs = 42 #81
    elseif(sReg == "3a")
        numBndryTurbs = 13 #31
    elseif(sReg == "3b")
        numBndryTurbs = 4 #11
    elseif(sReg == "4a")
        numBndryTurbs = 7 #16
    elseif(sReg == "4b")
        numBndryTurbs = 6 #14
    elseif(sReg == "4c")
        numBndryTurbs = 4 #9
    end
    
    return numBndryTurbs
end

function getCs34VertList(sReg::String)
    if(sReg == "cs3")
        vertList = [1, 10, 11, 13, 18]
    elseif(sReg == "3a")
        vertList = [1, 10, 11, 13, 18]
    elseif(sReg == "3b")
        vertList = [1, 6, 7, 8, 9]
    elseif(sReg == "4a")
        vertList = [1, 4, 5, 6, 7]
    elseif(sReg == "4b")
        vertList = [1, 2, 3, 4]
    elseif(sReg == "4c")
        vertList = [1, 2, 5, 6]
    elseif((sReg == "all") || (sReg == "All") || (sReg == "ALL"))
        vertList = [1, 10, 11, 13, 18,
                    1,  6,  7,  8,  9,
                    1,  4,  5,  6,  7,
                    1,  2,  3,  4,
                    1,  2,  5,  6]
    end
    
    return vertList
end

function getBndryCs4YAML(file_name)
    """Retreive boundary coordinates from the <.yaml> file"""

    # Read in the .yaml file
    f = YAML.load(open(file_name))
    bndrs = f["boundaries"]
    # Initialize some variables
    nRegions = 5;
    nPts = fill(0, 1, nRegions)
    ptList = []
    # Go through all our regions to get the boundary points
    for cntr in 1:nRegions
        ptList = push!(ptList, bndrs[getCs34NameYAML(cntr)])
        nPts[cntr] = floor(length(bndrs[getCs34NameYAML(cntr)]))
    end

    # Reorder points from 3b and 4a to be CCW from NE corner
    ptList[1] = circshift(ptList[1],-1)
    ptList[3] = circshift(ptList[3],-4)
    ptList[4] = circshift(ptList[4],-1)
    ptList[5] = circshift(ptList[5],-1)
    # Change from CW -> CCW
    for i in 1:nRegions
        ptList[i] = reverse(ptList[i])
    end

    # Initialize our arrays for the coordinates
    bndry_x = [ Float64[] for i in 1:nRegions ]
    bndry_y = [ Float64[] for i in 1:nRegions ]
    # Read in all the coordinates
    for i in 1:nRegions             # Looping through all regions
        for j in 1:nPts[i]          # Looping through all points in this region
            bndry_x[i] = push!(bndry_x[i], ptList[i][j][1])
            bndry_y[i] = push!(bndry_y[i], ptList[i][j][2])
        end
    end

    return bndry_x, bndry_y
end

using PyPlot
function printTurbinesInBoundary(bndry_x_clsd, bndry_y_clsd, turbine_x, turbine_y, turb_diam, region=0, turb_tags=false)
    #- Visualizes the farm boundaries and all turbines passed -#
    # Plot the boundary (defaults to all regions (when region==0), but will do only one if specified)
    if region == 0
        for i in 1:length(bndry_x_clsd)
            plot(bndry_x_clsd[i], bndry_y_clsd[i])
        end
    else
        plot(bndry_x_clsd[region], bndry_y_clsd[region])
    end
    
    # Plot the turbines
    for i = 1:length(turbine_x)
         plt.gcf().gca().add_artist(plt.Circle((turbine_x[i],turbine_y[i]), turb_diam/2.0, fill=true,color="black"))
    #     plt.gcf().gca().add_artist(plt.Circle((test_values_x[i],test_values_y[i]), turb_diam/2.0, fill=true,color="red"))
        if turb_tags
            plt.text(turbine_x[i]+turb_diam,turbine_y[i]+turb_diam, string(i))
        end
    end

    # Formatting
    axis("square")
    axis("off")
    plt.show()
end