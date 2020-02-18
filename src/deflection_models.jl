abstract type AbstractDeflectionModel end

struct DummyFunc <: AbstractDeflectionModel end


function deflection_model(loc,turbine::Turbine,::DummyFunc)

    dx = loc[1]-turbine.coord.x
    dy = loc[2]-(turbine.coord.y+deflection_y)
    dz = loc[3]-(turbine.coord.z+turbine.hub_height+deflection_z)

    del = sqrt(dy^2+dz^2)
    return del

end
