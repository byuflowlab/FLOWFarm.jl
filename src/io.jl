"""
    get_turb_loc_YAML(file_name)

read in turbine locations and related problem file names from .yaml

# Arguments
- `file_name::String`: path/and/name/of/location/file.yaml
"""
function get_turb_loc_YAML(file_name; returnaep=false)
    ### Retrieve turbine locations and auxiliary file names from <.yaml> file.
    ### Auxiliary (reference) files supply wind rose and turbine attributes.

    # Read in the .yaml file
    f = YAML.load(open(file_name))
    defs = f["definitions"]

    # Rip the (x,y) coordinates (Convert from <list> to <ndarray>)
    turb_coords = defs["position"]["items"]
    # println(turb_coords)
    nturbs = length(turb_coords)
    turbine_x = zeros(nturbs)
    turbine_y = zeros(nturbs)
    
    for i in 1:nturbs
        turbine_x[i] = turb_coords[i][1]
        turbine_y[i] = turb_coords[i][2]
    end

    # Rip the expected AEP, used for comparison
    AEP = defs["plant_energy"]["properties"]["annual_energy_production"]["default"]

    # Read the auxiliary filenames for the windrose and the turbine attributes (first one)
    fname_turb = string.(values(defs["wind_plant"]["properties"]["turbine"]["items"][1]))[1]
    fname_wr = string.(values(defs["plant_energy"]["properties"]["wind_resource"]["properties"]["items"][1]))[1]

    # Return turbine (x,y) locations, and the filenames for the others .yamls
    if returnaep
        return turbine_x, turbine_y, fname_turb, fname_wr, AEP
    else
        return turbine_x, turbine_y, fname_turb, fname_wr
    end 
end

"""
    get_turb_atrbt_YAML(file_name)

read in turbine attributes from .yaml

# Arguments
- `file_name::String`: path/to/attribute/file.yaml
"""
function get_turb_atrbt_YAML(file_name)
    ###Retreive turbine attributes from the <.yaml> file###

    # Read in the .yaml file
    f = YAML.load(open(file_name))
    defs = f["definitions"]
    ops = defs["operating_mode"]
    turb = defs["wind_turbine"]
    rotor = defs["rotor"]
    hub = defs["hub"]

    # Rip the turbine attributes
    # (Convert from <list> to <float>)
    turb_ci = float(ops["cut_in_wind_speed"]["default"])
    turb_co = float(ops["cut_out_wind_speed"]["default"])
    rated_ws = float(ops["rated_wind_speed"]["default"])
    rated_pwr = float(turb["rated_power"]["maximum"])
    turb_diam = float(rotor["diameter"]["default"])
    turb_height = float(hub["height"]["default"])

    return turb_ci, turb_co, rated_ws, rated_pwr, turb_diam, turb_height
end

"""
    get_wind_rose_YAML(file_name)

read in wind resource information from .yaml

# Arguments
- `file_name::String`: path/to/wind/resource/file.yaml
"""
function get_wind_rose_YAML(file_name)
    ### Retrieve wind rose data (bins, freqs, speeds) from <.yaml> file.

    # Read in the .yaml file
    f = YAML.load(open(file_name))
    props = f["definitions"]["wind_inflow"]["properties"]

    # Rip wind directional bins, their frequency, and the windspeed parameters for each bin
    # (Convert from <list> to <ndarray>)
    wind_dir = props["direction"]["bins"]
    wind_dir_freq = props["direction"]["frequency"]
    ndirs = length(wind_dir)
    # (Convert from <list> to <float>)
    wind_speeds = props["speed"]["bins"]
    wind_speed_probs = props["speed"]["frequency"]
    # Get default number of windspeed bins per direction
    nspeeds = length(wind_speeds)
    min_speed = props["speed"]["minimum"]
    max_speed = props["speed"]["maximum"]

    # get ambient ti
    ti = props["turbulence_intensity"]["default"]

    # convert to flow farm standard format
    freq = zeros(ndirs*nspeeds)
    speed = zeros(ndirs*nspeeds)
    dir = zeros(ndirs*nspeeds)

    # calculate frequency for each direction/speed combination
    for i in 1:ndirs
        for j in 1:nspeeds
            freq[(i-1)*nspeeds+j] = wind_dir_freq[i]*wind_speed_probs[i][j]
            speed[(i-1)*nspeeds+j] = wind_speeds[j]
            dir[(i-1)*nspeeds+j] = wind_dir[i]
        end
    end

    return dir, speed, freq, ti
end

"""
    write_turb_loc_YAML(file_name, data)

write turbine locations and related information to .yaml

# Arguments
- `file_name::String`: path/and/name/of/location/file.yaml
"""
function write_turb_loc_YAML(filename, turbinex, turbiney; title="", titledescription="", 
    turbinefile="", locunits="m", wakemodelused="", windresourcefile="", aeptotal=[], 
    aepdirs=[], aepunits="MWh", baseyaml=string(@__DIR__, "/default.yaml"))

    ### Retrieve turbine locations and auxiliary file names from <.yaml> file.
    ### Auxiliary (reference) files supply wind rose and turbine attributes.

    # Read in the default .yaml file
    base = YAML.load(open(baseyaml))

    # get number of turbines
    nturbines = length(turbinex)

    # save the title and description to the yaml database
    base["title"] = title
    base["description"] = titledescription

    # save the title and description to the yaml database
    base["definitions"]["plant_energy"]["properties"]["wake_model"]["items"][1]["\$ref"] = wakemodelused  

    # save positions in yaml database
    turb_coords = fill(zeros(2), nturbines)
    for i in 1:nturbines
        turb_coords[i] = [turbinex[i], turbiney[i]]
    end    
    
    base["definitions"]["position"]["items"] = turb_coords
    
    # save the AEP in yaml database
    base["definitions"]["plant_energy"]["properties"]["annual_energy_production"]["default"] = aeptotal

    # save the directional AEPs in the yaml database 
    base["definitions"]["plant_energy"]["properties"]["annual_energy_production"]["binned"] = aepdirs

    # save the directional AEP units in the yaml database 
    base["definitions"]["plant_energy"]["properties"]["annual_energy_production"]["units"] = aepunits

    # save the auxiliary filenames for the windrose and the turbine attributes
    base["definitions"]["wind_plant"]["properties"]["turbine"]["items"][1]["\$ref"] = turbinefile
    base["definitions"]["plant_energy"]["properties"]["wind_resource"]["properties"]["items"][1]["\$ref"] = windresourcefile

    # write result to .yaml
    YAML.write_file(filename, base)
end

function get_reduced_wind_rose_YAML(file_name)
    ### Retrieve wind rose data (bins, freqs, speeds) from <.yaml> file.
    # Only the average wind speeds for each direction are reported.

    # Read in the .yaml file
    f = YAML.load(open(file_name))
    props = f["definitions"]["wind_inflow"]["properties"]

    # Rip wind directional bins, their frequency, and the windspeed parameters for each bin
    # (Convert from <list> to <ndarray>)
    wind_dir = props["direction"]["bins"]
    wind_dir_freq = props["direction"]["frequency"]
    ndirs = length(wind_dir)
    # (Convert from <list> to <float>)
    wind_speeds = props["speed"]["bins"]
    wind_speed_probs = props["speed"]["frequency"]
    # Get default number of windspeed bins per direction
    nspeeds = length(wind_speeds)
    min_speed = props["speed"]["minimum"]
    max_speed = props["speed"]["maximum"]

    # get ambient ti
    ti = props["turbulence_intensity"]["default"]

    # convert to flow farm standard format
    freq = zeros(ndirs)
    speed = zeros(ndirs)
    dir = zeros(ndirs)

    # calculate average speed for each direction
    for i in 1:ndirs
        freq[i] = wind_dir_freq[i]
        speed[i] = sum(wind_speeds.*wind_speed_probs[i])
        dir[i] = wind_dir[i]
    end

    return dir, speed, freq, ti
end

"""
    getNextFileName(directory, file_name, file_type)

Checks if a file of the given directory and name exists. If not, increments to
the next index so as not to overwrite previously written files. If it reaches
the max number of overwrites, it will default to <directory/file_name.file_type>
To default to this, set <max_check=0> in function call.

# Arguments
- `directory::String`: path/to/write/file/at/
- `file_name::String`: Whatever the name of the file desired
- `file_type::String`: ex "yaml", "txt", "csv", etc...
- `max_check::Int`: the maximum number of files to check
"""
function getNextFileName(directory, file_name, file_type, max_check=100)
    full_file_name = string(directory, file_name, ".", file_type)
    for i in 1:100
        if !isfile(string(directory, file_name, "-(", i, ").", file_type))
            full_file_name = string(directory, file_name, "-(", i, ").", file_type)
            break
        end
    end
    return full_file_name
end