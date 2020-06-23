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

function getCs34VertList(sReg::String)
    if(sReg == "cs3")
        vertList = [1, 10, 11, 13, 19]
    elseif(sReg == "3a")
        vertList = [1, 10, 11, 13, 19]
    elseif(sReg == "3b")
        vertList = [1, 6, 7, 8, 9]
    elseif(sReg == "4a")
        vertList = [1, 4, 5, 6, 7]
    elseif(sReg == "4b")
        vertList = [1, 2, 3, 4]
    elseif(sReg == "4c")
        vertList = [1, 2, 5, 6]
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
    x_boundary_coords = [ Float64[] for i in 1:nRegions ]
    y_boundary_coords = [ Float64[] for i in 1:nRegions ]
    # Read in all the coordinates
    for i in 1:nRegions             # Looping through all regions
        for j in 1:nPts[i]          # Looping through all point sin this region
            x_boundary_coords[i] = push!(x_boundary_coords[i], ptList[i][j][1])
            y_boundary_coords[i] = push!(y_boundary_coords[i], ptList[i][j][2])
        end
    end

    return x_boundary_coords, y_boundary_coords
end