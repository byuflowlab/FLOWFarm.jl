abstract type AbstractCostParameter end

"""
    CostParameters(TCC, BOS, FC, FCR, OpEx)

Container for parameters related to the Levelized Cost of Energy model (NREL 2016 Cost of Wind Energy)

# Arguments
# Arguments
- `TCC::Float`: Turbine Capital Cost not including the tower module
- `BOS::Float`: Balance of System (Costs outside of turbine i.e. operation and maintenance)
- `FC::Float`: Financial Costs including construction and contingency
- `FCR::Float`: Fixed Charge Rate
- `OpEx::Float`: Operational Expenditures
"""
struct CostParameters{TF} <: AbstractCostParameter
    TCC::TF
    BOS::TF
    FC::TF
    FCR::TF
    OpEx::TF
end
CostParameters() = CostParameters(776, 326, 120, .0655, 43) # Default values taken from NREL 2019 Cost of Wind Energy


"""
    cost_of_energy(rotor_diameter, hub_height, AEP, Cost::CostParameters)

Calculates the LCOE using the same numbers as NREL's FLORIS Model

# Arguments
- `rotor_diameter::array`: Vector of Rotor Diameters for the Turbines
- `hub_height::array`: Vector of Hub Heights for the Turbines
- `AEP::Float`: Annual Energy Production
- `OpEx::AbstractCostParameter`: KW of the Farm
"""

function cost_of_energy(rotor_diameter, hub_height, AEP, Cost::CostParameters)
    
    nturbines = length(rotor_diameter)
    # Taken from 2019 Cost of Wind Energy Review
    TCC = Cost.TCC
    BOS = Cost.BOS
    FC = Cost.FC
    FCR = Cost.FCR
    OpEx = Cost.OpExp

    Mass = 0
    for i = 1:nturbines
        swept_area = pi*(rotor_diameter[i]/2)^2
        Mass += .2694*hub_height[i]*swept_area + 1779.3
    end
    # Adding the mass of the turbines
    TCC =  TCC*nturbines + 3.08*Mass/PlantKW

    # Uses parameters in COE function from eq 1 in 2016 Cost of Wind Energy Review
    LCOE = ((TCC+BOS+FC)*FCR + OpEx)/(AEP/1000)

    return LCOE # LCOE is in units of $/kWh
end