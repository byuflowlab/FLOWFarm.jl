import FLOWFarm; const ff = FLOWFarm
    # set initial turbine x and y locations
    diam = 126.0
    turbine_x = [0.0,0.0]
    nturbines = length(turbine_x)
    turbine_y = [0.0,0.0]
    turbine_z = zeros(nturbines)
    
    # set turbine yaw values
    turbine_yaw = zeros(nturbines)
    
    # set turbine design parameters
    rotor_diameter = zeros(nturbines) .+ diam # m
    hub_height = zeros(nturbines) .+ 90.0   # m

    sorted_index = sortperm(turbine_x)
    turbine_local_ti = zeros(length(turbine_x))
    ambient_ti = .137
    ti_model = ff.LocalTIModelGaussTI()