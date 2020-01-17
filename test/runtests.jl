import FlowFarm; const ff = FlowFarm
using Test

@testset "Wake Models" begin
    # Write your own tests here.
    
    @testset "Jensen Model" begin
        coord = ff.Coord(0.0,0.0,0.0)
        rotor_diameter = 40.0
        hub_height = 90.0
        aI = 1.0/3.0
        gamma = 0.0
        yaw = 0.0
        ct = 0.7 #TODO handle ct and axial induction appropriately
        
        turbine = ff.Turbine(coord, rotor_diameter, hub_height, aI, gamma, yaw, ct)
        
        deflection = [0.0, 0.0]

        model = ff.Jensen(0.1)

        centerloss40 = 1. - 4.35/8.1
        centerloss100 = 1. - 5.7/8.1

        @test ff.wake_model([40., 0.0, hub_height], deflection, model, turbine) == centerloss40
        @test ff.wake_model([100., 0.0, hub_height], deflection, model, turbine) == centerloss100

    end
    
end
