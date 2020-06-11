# This is a testfile to ensure the splined_boundary() function works.

### Complete and functional ###
function getTurbLocYAML(file_name)
    ### Retrieve turbine locations and auxiliary file names from <.yaml> file.
    ### Auxiliary (reference) files supply wind rose and turbine attributes.

    # Read in the .yaml file
    f = YAML.load(open(file_name))
    defs = f["definitions"]

    # Rip the (x,y) coordinates (Convert from <list> to <ndarray>)
    turb_coords = defs["position"]["items"]

    # Rip the expected AEP, used for comparison
    AEP = defs["plant_energy"]["properties"]["annual_energy_production"]["default"]

    # Read the auxiliary filenames for the windrose and the turbine attributes (first one)
    fname_turb = string.(values(defs["wind_plant"]["properties"]["turbine"]["items"][1]))[1]
    fname_wr = string.(values(defs["plant_energy"]["properties"]["wind_resource"]["properties"]["items"][1]))[1]

    # Return turbine (x,y) locations, and the filenames for the others .yamls
    return turb_coords, fname_turb, fname_wr
end
### Complete and functional ###
function getTurbAtrbtYAML(file_name)
    ###Retreive turbine attributes from the <.yaml> file###

    # Read in the .yaml file
    f = YAML.load(open(file_name))
    defs = f["definitions"]
    ops = defs["operating_mode"]
    turb = defs["wind_turbine"]
    rotor = defs["rotor"]

    # Rip the turbine attributes
    # (Convert from <list> to <float>)
    turb_ci = float(ops["cut_in_wind_speed"]["default"])
    turb_co = float(ops["cut_out_wind_speed"]["default"])
    rated_ws = float(ops["rated_wind_speed"]["default"])
    rated_pwr = float(turb["rated_power"]["maximum"])
    turb_diam = float(rotor["diameter"]["default"])

    return turb_ci, turb_co, rated_ws, rated_pwr, turb_diam
end
### Complete and functional ###
function getWindRoseYAML(file_name)
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

    return wind_dir, wind_dir_freq, wind_speeds, wind_speed_probs, num_speed_bins, min_speed, max_speed
end

### Incomplete, untested ###
function getBndryCs4YAML(file_name)
    ###Retreive boundary coordinates from the <.yaml> file

    # Read in the .yaml file
    f = YAML.load(open(file_name))
    bndrs = f["boundaries"]
        
    ptList3a = bndrs['IIIa']
    ptList3b = bndrs['IIIb']
    ptList4a = bndrs['IVa']
    ptList4b = bndrs['IVb']
    ptList4c = bndrs['IVc']

    # (Convert from <list> to <coordinate> array)
    # Determine how many points we have for the boundary
    # numCoords3a = len(ptList3a)
    # numCoords3b = len(ptList3b)
    # numCoords4a = len(ptList4a)
    # numCoords4b = len(ptList4b)
    # numCoords4c = len(ptList4c)
    # # Initialize the point list array
    # coordList3a = np.recarray(numCoords3a, coordinate)
    # coordList3b = np.recarray(numCoords3b, coordinate)
    # coordList4a = np.recarray(numCoords4a, coordinate)
    # coordList4b = np.recarray(numCoords4b, coordinate)
    # coordList4c = np.recarray(numCoords4c, coordinate)
    print(ptList3a)
end


#--- Read in the data ---#
scaledAEP = 1#e5
scaledTC = 1#e3
strCase = 'cs4'  # Which case study we're doing. 'cs3' or 'cs4'
numGridLines = 10                   # How many gridlines we'll use for the splining
nNumRegions = 5                     # Number of reigons we're using (cs4 = 5, cs3 = 1)

#- Rip the boundary coordinates from the .yaml file -# 
bnry_file_name = "../startup-files/iea37-boundary-" + strCase + ".yaml"
getBndryCs4YAML(bnry_file_name)

# Augment the data so it's CCW and from NW "corner"