using FlowFarm; const ff=FlowFarm
using DelimitedFiles


"""
    COE(FCR, ICC, AEP, LLC, OandM, LRC)

Calculates the Levelized Cost of Energy for a Given Wind Farm

# Arguments
- `FCR::Float`: Fixed Charge Rate
- `ICC::Float`: Initial Capital Cost
- `AEP::Float`: Annual Energy Production
- `LLC::Float`: Land Lease Cost
- `OandM::Float`: Levelized O&M Cost
- `LRC::`: Levelized Replacement/Overhaul Cost 
"""

function COE(FCR, ICC, AEP, LLC, OandM, LRC)

    AOE = LLC + (OandM + LRC)/AEP # Annual Operating Expense eq from pg. 5 of the paper
    COE = (FCR*ICC)/AEP + AOE # Levelized Cost of Energy eq from pg.4 of the paper

    return  COE
end

"""
    Levelized_Cost_Of_Energy(turbine)

Imports LCOE parameters based on turbine selection then calculates LCOE using COE function

# Arguments
- `turbine::String`: Turbine Selection, Available options (Exact Syntax); NREL_5MW, Vestas_V80, Siemens_swp_3.6
"""

function Levelized_Cost_Of_Energy(turbine)
    
    # Imports available data from flow farm memory
    turbine_data = readdlm("test/inputfiles/Turbine_COE_Parameters.txt", skipstart=2)

    # Identifies the row from which to select the data based on the turbine selection
    ind = q=findall(x->x==turbine,turbine_data)[1][1]

    # Finds parameters to be used in calculation
    FCR = turbine_data[ind,2]
    ICC = turbine_data[ind,3]
    AEP = turbine_data[ind,4]
    LLC = turbine_data[ind,5]
    OandM = turbine_data[ind,6]
    LRC = turbine_data[ind,7]

    # Uses parameters in COE function
    LCOE = COE(FCR,ICC,AEP,LLC,OandM,LRC)

    return LCOE
end

"""
    COE(FCR, ICC, AEP, LLC, OandM, LRC)

Calculates the Levelized Cost of Energy for a Given Wind Farm

# Arguments
- `FCR::Float`: Fixed Charge Rate
- `ICC::Float`: Initial Capital Cost
- `AEP::Float`: Annual Energy Production
- `LLC::Float`: Land Lease Cost
- `OandM::Float`: Levelized O&M Cost
- `LRC::`: Levelized Replacement/Overhaul Cost 
"""

function COE(FCR, ICC, AEP, LLC, OandM, LRC)

    AOE = LLC + (OandM + LRC)/AEP # Annual Operating Expense eq from pg. 5 of the paper
    COE = (FCR*ICC)/AEP + AOE # Levelized Cost of Energy eq from pg.4 of the paper

    return  COE
end

"""
    Levelized_Cost_Of_Energy(turbine)

Calculates the LCOE using the same numbers as NREL's FLORIS Model

# Arguments
- `nturbines::Int`: Number of Turbines in the farm
- `rotor_diameter::array`: Vector of Rotor Diameters for the Turbines
- `hub_height::array`: Vector of Hub Heights for the Turbines
- `PlantKW::Float`: KW of the Farm
- `AEP::Float`: Annual Energy Production
"""

function Levelized_Cost_Of_Energy(nturbines, rotor_diameter, hub_height, PlantKW, AEP)
    
    PlantKW = nturbines*5000
    Mass = 0
    for i = 1:nturbines
        swept_area = pi*(rotor_diameter[i]/2)^2
        Mass += .2694*hub_height[i]*swept_area + 1779.3
    end
    # Values taken from Table 3 in CSM (831 is TCC - Tower in table)
    TCC = 831 + 3.08*Mass/PlantKW
    BOS = 364
    FC = 155

    # Taken from 2016 Cost of Wind Energy Review
    FCR = .079

    # Uses parameters in COE function from eq 1 in 2016 Cost of Wind Energy Review
    LCOE = ((TCC+BOS+FC)*FCR)/(AEP/1000/PlantKW)

    return LCOE
end