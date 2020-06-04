import YAML

## Complete and functional ###
function get_turb_loc_YAML(file_name)
    ### Retrieve turbine locations and auxiliary file names from <.yaml> file.
    ### Auxiliary (reference) files supply wind rose and turbine attributes.

    # Read in the .yaml file
    f = YAML.load(open(file_name))
    defs = f["definitions"]

    # Rip the (x,y) coordinates (Convert from <list> to <ndarray>)
    turb_coords = defs["position"]["items"]
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
    return turbine_x, turbine_y, fname_turb, fname_wr
end

### Complete and functional ###
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

### Complete and functional ###
function get_wind_rose_YAML(file_name)
    ### Retrieve wind rose data (bins, freqs, speeds) from <.yaml> file.

    # Read in the .yaml file
    f = YAML.load(open(file_name))
    props = f["definitions"]["wind_inflow"]["properties"]

    # Rip wind directional bins, their frequency, and the windspeed parameters for each bin
    # (Convert from <list> to <ndarray>)
    wind_dir = props["direction"]["bins"]
    wind_dir_freq = props["direction"]["frequency"]
    # (Convert from <list> to <float>)
    wind_speeds = props["speed"]["bins"]
    wind_speed_probs = props["speed"]["frequency"]
    # Get default number of windspeed bins per direction
    num_speed_bins = length(wind_speeds)
    min_speed = props["speed"]["minimum"]
    max_speed = props["speed"]["maximum"]

    # get ambient ti
    ti = props["turbulence_intensity"]["default"]

    # convert to flow farm standard format
    freq = zeros(length(wind_dir_freq)*length(wind_speed_probs))
    speed = zeros(length(wind_dir_freq)*length(wind_speed_probs))
    dir = zeros(length(wind_dir_freq)*length(wind_speed_probs))

    for i in 1:length(wind_dir_freq)
        for j in 1:length(wind_speed_probs)
            freq[(i-1)*num_speed_bins+j] = wind_dir_freq[i]*wind_speed_probs[i][j]
            speed[(i-1)*num_speed_bins+j] = wind_speeds[j]
            dir[(i-1)*num_speed_bins+j] = wind_dir[i]
        end
    end

    return dir, speed, freq, ti
end