abstract type AbstractCostModel end

"""
    Levelized(TCC, BOS, FC, FCR, OpEx)

Container for parameters related to the Levelized Cost of Energy model (NREL 2016 Cost of Wind Energy)

# Arguments
# Arguments
- `TCC::Float`: Turbine Capital Cost not including the tower module
- `BOS::Float`: Balance of System (Costs outside of turbine i.e. operation and maintenance)
- `FC::Float`: Financial Costs including construction and contingency
- `FCR::Float`: Fixed Charge Rate
- `OpEx::Float`: Operational Expenditures
"""
struct Levelized{TF} <: AbstractCostModel
    TCC::TF
    BOS::TF
    FC::TF
    FCR::TF
    OpEx::TF
end
Levelized() = Levelized(776.0, 326.0, 120.0, .0655, 43.0) # Default values taken from NREL 2019 Cost of Wind Energy


"""
    Floating_turbine(TCC, BOS, FC, FCR, OpEx)

Container for parameters related to floating turbine capital cost (NREL 2020 ORBIT)

# Arguments
- `TCC::Float`: Turbine Capital Cost not including the tower module
- `BOS::Float`: Balance of System (Costs outside of turbine i.e. operation and maintenance)
- `FC::Float`: Financial Costs including construction and contingency
- `FCR::Float`: Fixed Charge Rate
- `OpEx::Float`: Operational Expenditures
"""
struct Floating_turbine_model_set{FTCC, FFC, FFCR, OpEx} <: AbstractFloatingModelSet
    Floating_TCC::FTCC
    Floating_FC::FFC
    Floating_FCR::FFCR
    Floating_OpEx::OpEx
end


"""
    cost_of_energy_land(rotor_diameter, hub_height, AEP, Cost::Levelized)

Calculates the LCOE using the same numbers as NREL's FLORIS Model

# Arguments
- `rotor_diameter::array`: Vector of Rotor Diameters for the Turbines
- `hub_height::array`: Vector of Hub Heights for the Turbines
- `rated_power::array`: Vector of rated powers for the Turbines in kW
- `AEP::Float`: Annual Energy Production
- `OpEx::AbstractCostParameter`: KW of the Farm
"""

function cost_of_energy_land(rotor_diameter, hub_height, rated_power, AEP, Cost::Levelized)
    
    nturbines = length(rotor_diameter)
    # Taken from 2019 Cost of Wind Energy Review
    TCC = Cost.TCC
    BOS = Cost.BOS
    FC = Cost.FC
    FCR = Cost.FCR
    OpEx = Cost.OpEx

    PlantKW = sum(rated_power) # Combines total Wind Farm Capacity

    Mass = 0.0
    for i = 1:nturbines
        swept_area = pi*(rotor_diameter[i]/2)^2
        Mass += .2694*hub_height[i]*swept_area + 1779.3
    end
    # Adding the mass of the turbines to the Total Capitol Cost
    TCC =  TCC + 3.08*Mass/PlantKW

    # Uses parameters in COE function from eq 1 in 2016 Cost of Wind Energy Review
    LCOE = ((TCC+BOS+FC)*FCR + OpEx)/(AEP/1000)

    return LCOE # LCOE is in units of $/kWh
end


"""
    cost_of_energy_floating(rotor_diameter, hub_height, AEP, Cost::Levelized)

Calculates the LCOE using the same numbers as NREL's FLORIS Model and BOS cost from NREL's ORBIT

# Arguments
- `turbine_type::Float`: 0 inidicates a land-based wind farm, 1 indicates a floating offshore wind farm
- `rotor_diameter::array`: Vector of Rotor Diameters for the Turbines
- `hub_height::array`: Vector of Hub Heights for the Turbines
- `rated_power::array`: Vector of rated powers for the Turbines in kW
- `AEP::Float`: Annual Energy Production
- `OpEx::AbstractCostParameter`: KW of the Farm
"""

function cost_of_energy_floating(turbine_x, turbine_y, turbine_z, turbine_ocean_depth, rotor_diameter,
    hub_height, turbine_yaw, ct_model, generator_efficiency, cut_in_speed,
    cut_out_speed, rated_speed, rated_power, wind_resource, power_models, model_set::AbstractModelSet, CostModel::AbstractFloatingModelSet;
    rotor_sample_points_y=[0.0], rotor_sample_points_z=[0.0], hours_per_year=365.25*24.0, distributed=false)
    
    nturbines = length(rotor_diameter)
    rotorsamplepointsy = rotor_sample_points_y
    rotorsamplepointsz = rotor_sample_points_z
    

    # solve for AEP
    AEP = calculate_aep(turbine_x, turbine_y, turbine_z, rotor_diameter,
    hub_height, turbine_yaw, ct_model, generator_efficiency, cut_in_speed,
    cut_out_speed, rated_speed, rated_power, wind_resource, power_models, model_set,
    rotor_sample_points_y=rotorsamplepointsy, rotor_sample_points_z=rotorsamplepointsz)


    TCC_model = model_set.Floating_TCC

    # Taken from 2020 ORBIT
    # TCC = Floating_TCC(rotor_diameter, hub_height, rated_power, TCC_model) # default is 776.0
    TCC = 776.0
    BOS = Floating_BOS(rotor_diameter, turbine_x, turbine_y, turbine_z, turbine_ocean_depth, rated_power, drag_embedment, Mooring_per_turbine, drag_embedment_fixed_length)
    FC = CostModel.FFC
    FCR = CostModel.FFCR
    OpEx = CostModel.OpEx

    # Uses parameters in COE function from eq 1 in 2016 Cost of Wind Energy Review
    LCOE = ((TCC+BOS+FC)*FCR + OpEx)/(AEP/1000)

    return LCOE # LCOE is in units of $/kWh
end


"""
    Floating_TCC(rotor_diameter, hub_height, rated_power, turbine_model::Floating_turbine)

Turbine capital cost calculations

# Arguments
- `rotor_diameter::array`: Vector of Rotor Diameters for the Turbines
- `hub_height::array`: Vector of Hub Heights for the Turbines
- `rated_power::array`: Vector of rated powers for the Turbines in kW
- `turbine_model::AbstractCostParameter`: TCC capital cost from orbit 2020
"""

function Floating_TCC(rotor_diameter, hub_height, rated_power, turbine_model::Floating_turbine)
    
    nturbines = length(rotor_diameter)
    # Taken from 2019 Cost of Wind Energy Review
    TCC = turbine_model.TCC
    

    return TCC 
end

"""
    Floating_BOS(rotor_diameter, turbine_x, turbine_y, turbine_z, turbine_ocean_depth)

Calculates the balance of system cost (BOS) using the same numbers as NREL's ORBIT (2020)

# Arguments
- `rotor_diameter::array`: Vector of Rotor Diameters for the Turbines
- `turbine_x::Array{Float,nTurbines}`: turbine east-west locations in the global 
    reference frame
- `turbine_y::Array{Float,nTurbines}`: turbine north-south locations in the global 
    reference frame
- `turbine_z::Array{Float,nTurbines}`: turbine base height in the global reference frame
- `turbine_ocean_depth::Array{Float,nTurbines}`: ocean depth at location of each turbine
- `rated_power::Array{Float,nTurbines}`
- `drag_embedment::Float`: parameter determining if drag embedment is used (drag_embedment = 1) or not (drag_embedment = 0)
- `Mooring_per_turbine::Float`: the number of mooring lines per turbine
- `drag_embedment_fixed_length::Float`: drag embedment fixed length (default is 1)
"""

function Floating_BOS(rotor_diameter, turbine_x, turbine_y, turbine_z, turbine_ocean_depth, rated_power,
            drag_embedment, Mooring_per_turbine, drag_embedment_fixed_length, Onshore_substation_x, Onshore_substation_y, substation_x, substation_y,
            substation_z)
    
    nturbines = length(rotor_diameter)
    # Taken from 2019 Cost of Wind Energy Review

    # Mooring design (no shared mooring lines)
    # TODO: include option of shared mooring line design
    Mooring_BOS = Mooring_Design(rated_power, drag_embedment, turbine_ocean_depth, Mooring_per_turbine, drag_embedment_fixed_length)
    # Semisubmersible design
    Semisub_BOS = Semisubmersible_Design(rated_power)
    # Electrical cabling design
    # Put Weston's cable design code here
    Num_Array2, Cabling_BOS = Cabling_Design(turbine_x, turbine_y, turbine_z, turbine_ocean_depth)

    # does export cabling design go before substation design?
    # Export sytem design
    Export_Cable_BOS = Export_Cable_Design(rated_power, Onshore_substation_x, Onshore_substation_y, substation_x, substation_y, substation_z)

    # Substation design
    Substation_BOS = Substation_Design(rated_power, substation_z, drag_embedment, Mooring_per_turbine)

    BOS = Substation_BOS + Cabling_BOS + Semisub_BOS + Mooring_BOS + Export_Cable_BOS
    return BOS # LCOE is in units of $/kWh
end

"""
    Mooring_Design(rated_power, drag_embedment, turbine_ocean_depth, Mooring_per_turbine, drag_embedment_fixed_length=1.0)

Calculates the balance of system cost (BOS) for mooring lines using same calculations from NREL's ORBIT (2020)

# Arguments
- `rated_power::Array{Float,nTurbines}`: the rated power of each turbine in MW
- `drag_embedment::Float`: parameter determining if drag embedment is used (drag_embedment = 1) or not (drag_embedment = 0)
- `Mooring_per_turbine::Float`: the number of mooring lines per turbine
- `turbine_ocean_depth::Array{Float,nTurbines}`: ocean depth at location of each turbine
- `drag_embedment_fixed_length::Float`: drag embedment fixed length
"""

function Mooring_Design(rated_power, drag_embedment, turbine_ocean_depth, Mooring_per_turbine, drag_embedment_fixed_length=1.0)

    nturbines = length(rated_power)
    Mooring_BOS = 0

    # for each turbine calculate the cost of mooring lines and add to Mooring_BOS
    for i = 1:nturbines
        # determine size of mooring line
        fit = -0.0004*(rated_power[i])^2 + 0.0132*(rated_power[i]) +0.0536
        if fit <= 0.09
            line_diam = 0.09        # in meters
            line_mass_per_m = 0.161
            line_cost_rate = 399.0
        elseif fit <= 0.12
            line_diam = 0.12        # in meters
            line_mass_per_m = 0.288
            line_cost_rate = 721.0
        else
            line_diam = 0.15        # in meters
            line_mass_per_m = 0.450
            line_cost_rate = 1088.0
        end

        # Breaking load of mooring line
        Breaking_load = 419449*(line_diam^2) + 93415*line_diam - 3577.9

        # Line length mass
        if drag_embedment == 1
            fixed = 500*drag_embedment_fixed_length
        else
            fixed = 0
        end

        line_length = 0.0002*(turbine_ocean_depth[i]^2) + 1.264*(turbine_ocean_depth[i]) + 47.776 + fixed
        line_mass = line_length*line_mass_per_m

        # Anchor mass cost
        if drag_embedment == 1
            anchor_mass = 20
            anchor_cost = (Breaking_load/9.81/20.0)*2000.0
        else
            anchor_mass = 50
            anchor_cost = sqrt(Breaking_load/9.81/1250)*150000
        end

        # sum of mooring lines and anchor cost
        Mooring_cost = Mooring_per_turbine*(anchor_cost+(line_length*line_cost_rate))
        Mooring_BOS = Mooring_BOS + Mooring_cost
    end

    return Mooring_BOS
end

"""
    Semisubmersible_Design(rated_power, Cost_rate_column, Cost_rate_truss, Cost_rate_heave_plate, Cost_rate_steel)

Calculates the balance of system cost (BOS) for the semisubmerisble design using same calculations from NREL's ORBIT (2020)

# Arguments
- `rated_power::Array{Float,nTurbines}`: the rated power of each turbine in MW
- `Cost_rate_column::Float`: cost rate of the stiffened column
- `Cost_rate_truss::Float`: cost rate of the truss
- `Cost_rate_heave_plate::Float`: cost rate of the heave plate
- `Cost_rate_steel::Float`: cost rate of the steel
"""

function Semisubmersible_Design(rated_power, Cost_rate_column=3120, Cost_rate_truss=6250, Cost_rate_heave_plate=6250, Cost_rate_steel=7250)
    nturbines = length(rated_power)
    total_mass = 0
    total_cost = 0

    for i = 1:nturbines

        # stiffened column mass
        stiff_colomn_mass = -0.9581*(rated_power[i]^2) + 40.89*(rated_power[i]) + 802.09        # in tonnes

        # stiffend column cost
        stiff_column_cost = stiff_colomn_mass*(Cost_rate_column)

        # truss mass
        truss_mass = 2.7894*(rated_power[i]^2) + 15.591*(rated_power[i]) + 266.03

        # truss cost
        truss_cost = truss_mass*Cost_rate_truss

        # heave plate mass
        heave_plate_mass = -0.4397*(rated_power[i]^2) + 21.545*(rated_power[i]) + 177.42

        # heave plate cost
        heave_plate_cost = heave_plate_mass*Cost_rate_heave_plate

        # secondary steel mass
        secondary_steel_mass = -0.153*(rated_power[i]^2) + 6.54*(rated_power[i]) + 128.34

        # secondary steel cost
        secondary_steel_cost = secondary_steel_mass*Cost_rate_steel

        # Total substructure mass and cost
        turbine_substructure_mass = stiff_colomn_mass + truss_mass + heave_plate_mass + secondary_steel_mass
        turbine_substructure_cost = stiff_column_cost + truss_cost + heave_plate_cost + secondary_steel_cost

        # sum up all turbine substructure mass and cost
        total_mass = total_mass + turbine_substructure_mass
        total_cost = total_cost + turbine_substructure_cost
    end
    return total_cost
end

"""
    Cabling_Design(turbine_x, turbine_y, turbine_z, turbine_ocean_depth)

Designs the cabling of the wind farm and then calculates the balance 
of system cost (BOS) for the cabling design using same calculations from NREL's ORBIT (2020)

# Arguments
- `turbine_x::Array{Float,nTurbines}`: turbine east-west locations in the global 
    reference frame
- `turbine_y::Array{Float,nTurbines}`: turbine north-south locations in the global 
    reference frame
- `turbine_z::Array{Float,nTurbines}`: turbine base height in the global reference frame
- `turbine_ocean_depth::Array{Float,nTurbines}`: ocean depth at location of each turbine
"""

function Cabling_Design(turbine_x, turbine_y, turbine_z, turbine_ocean_depth)
    # considerations:
    # cable capacity must increase as arrays of interconnected turbines grow
    # Can have different cable sizes or do just two different cable sizes


    return Num_Array2, Cabling_BOS
end

"""
    Export_Cable_Design(Onshore_substation_x, Onshore_substation_y, substation_x, substation_y)

Calculates the balance of system cost (BOS) for the export cable design using the same calculations from NREL's ORBIT (2020)

# Arguments
- `rotor_diameter::array`: Vector of Rotor Diameters for the Turbines
"""

function Export_Cable_Design(rated_power, Onshore_substation_x, Onshore_substation_y, substation_x, substation_y, substation_depth)
    # pick a point along the shore that represents the onshore substation
    export_diff_x = abs(Onshore_substation_x - substation_x)
    export_diff_y = abs(Onshore_substation_y - substation_y)

    Export_cable_distance = sqrt(export_diff_x^2 + export_diff_y^2)

    # find cost now based on export cable Export_cable_distance
    # these will need to be imported from .yaml files in later iterations
    # 0 means it is HVDC, 1 means it isn't
    HVDC = 0

    if HVDC == 1
        rated_voltage = 320     # kV
        current_capaticy = 1900     # A
        char_impedence = 1
        phase_angle = atan(imag(char_impedence)/real(char_impedence))
        power_factor = cos(phase_angle)

        cable_power = sqrt(3)*(rated_voltage * current_capaticy * power_factor/1000)

    else
        rated_voltage = 320     # kV
        current_capaticy = 1900     # A
        cable_power = sqrt(3)*(rated_voltage * 2 * current_capaticy/1000)
    end
    Num_cables = sum(rated_power)/(cable_power)

    # find required length of each cable
    system_angle = -0.0047*substation_depth + 18.743
    Catenary_length_factor = 0.04
    added_length = 500
    free_cable_length = (substation_depth/cos(system_angle*pi/180))*(Catenary_length_factor + 1) + 190
    Export_cable_distance = Export_cable_distance + free_cable_length + added_length

    Total_cable_length = (Num_cables*Export_cable_distance)*1.1

    # cost per km of cable
    cost_per_km = 828000        # $/km

    return (Total_cable_length/1000)*cost_per_km
end

"""
    Cable_define()

Calculates the balance of system cost (BOS) for the substation using the same calculations from NREL's ORBIT (2020)

# Arguments
- `rated_power::Array{Float,nTurbines}`: the rated power of each turbine in MW

"""

function Cable_define()
    
    return cable
end

"""
    Substation_Design(rated_power, substation_ocean_depth, drag_embedment, Mooring_per_turbine, mpt_cost_rate=12500, topside_fab_cost_rate=14500, 
    topside_design_cost=4.5e6, shunt_cost_rate=35000, switchgear_cost_rate=14.5e5, backup_gen_cost=1e6
    worspace_cost=2e6, other_ancillary_cost=3e6, topside_assembly_factor=0.075, oss_substructure_cost_rate=3000)

Calculates the balance of system cost (BOS) for the substation using the same calculations from NREL's ORBIT (2020)

# Arguments
- `rated_power::Array{Float,nTurbines}`: the rated power of each turbine in MW
- `substation_ocean_depth::Float`: ocean depth at location of substation
- `drag_embedment::Float`: parameter determining if drag embedment is used (drag_embedment = 1) or not (drag_embedment = 0)
- `Mooring_per_turbine::Float`: the number of mooring lines per turbine
- `mpt_cost_rate::Float`: cost of main power transforms (default is 12500)
- `topside_fab_cost_rate::Float`: topside fabrication cost rate (default is 14500)
- `topside_design_cost::Float`: cost of designing the topside (default is 4.5e6)
- `shunt_cost_rate::Float`: cost rate of the shunt reactors (default is 35000)
- `switchgear_cost_rate::Float`: cost rate of the switch gears (default is 14.5e5)
- `backup_gen_cost::Float`: cost of backup generator (default is 1e6)
- `workspace_cost::Float`: cost of workspace (default is 2e6)
- `other_ancillary_cost::Float`: other ancillary costs (default is 3e6)
- `topside_assembly_factor::Float`: topside assembly factor (default is 0.075)
- `oss_substructure_cost_rate::Float`: offshore sub station cost rate (default is 3000)

"""

function Substation_Design(rated_power, substation_ocean_depth, drag_embedment, Mooring_per_turbine, mpt_cost_rate=12500, topside_fab_cost_rate=14500, 
            topside_design_cost=4.5e6, shunt_cost_rate=35000, switchgear_cost_rate=14.5e5, backup_gen_cost=1e6
            worspace_cost=2e6, other_ancillary_cost=3e6, topside_assembly_factor=0.075, oss_substructure_cost_rate=3000)
    # nturbines = length(rated_power)
    # number_of_substations = sum(rated_power)/800      # 800 MW

    # substructure length
    substructure_length = substation_ocean_depth + 10       # meters

    # substructure deck space (Under developement with NREL's ORBIT)

    # topside deck spaced (Under developement with NREL's ORBIT)

    # num mpt (number of main power transformers) and rating
    farm_capacity = sum(rated_power)

    # TODO: allow for multiple substations
    #       needed if farm capacity exceeds 500 MW
    # Num_substations = ceil(farm_capacity/500)
    Num_substations = 1
    Num_mpt = ceil(farm_capacity/(250*Num_substations))
    mpt_rating = round(((farm_capacity*1.15)/(Num_mpt*Num_substations))/10.0)*10.0

    # mpt cost
    mpt_cost = mpt_rating*Num_mpt*mpt_cost_rate

    # topside mass and cost
    topside_mass = 3.85 * mpt_rating * Num_mpt + 285
    topside_cost = (topside_mass * topside_fab_cost_rate) + topside_design_cost

    # shunt reactor cost
    shunt_reactor_cost = mpt_rating * Num_mpt * shunt_cost_rate * 0.5

    # switchgear cost
    switchgear_cost = Num_mpt * switchgear_cost_rate

    # ancillary system cost
    Ancillary_system_cost = backup_gen_cost + workspace_cost + other_ancillary_cost

    # assembly cost
    Assembly_cost = (switchgear_cost + shunt_reactor_cost + Ancillary_system_cost) * topside_assembly_factor

    # substructure mass and cost
    substructure_mass = 0.4*topside_mass
    substructure_cost = substructure_mass * oss_substructure_cost_rate

    # mooring line and anchor cost
    rated_power_avg = mean(rated_power)
    mooring_cost = Mooring_Design(rated_power_avg, drag_embedment, substation_ocean_depth, Mooring_per_turbine)


    # total cost
    Substation_BOS = (substructure_cost + mpt_cost + topside_cost + shunt_reactor_cost + switchgear_cost + Ancillary_system_cost + Assembly_cost) + (mooring_cost)
    return Substation_BOS
end
