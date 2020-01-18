import FlowFarm; const ff = FlowFarm
using Test

@testset "Wake Models" begin
    # Write your own tests here.
    
    @testset "Jensen Top Hat Model" begin
        coord = ff.Coord(0.0,0.0,0.0)
        rotor_diameter = 40.0
        hub_height = 90.0
        aI = 1.0/3.0
        yaw = 0.0
        ct = 0.7 #TODO handle ct and axial induction appropriately

        alpha = 0.1
        
        turbine = ff.Turbine(coord, rotor_diameter, hub_height, aI, yaw, ct)
        
        deflection = [0.0, 0.0]

        model = ff.JensenTopHat(alpha)

        centerloss40 = 1. - 4.35/8.1
        centerloss100 = 1. - 5.7/8.1

        # test no loss upstream (data from Jensen 1983)
        @test ff.wake_model([-1E-12, 0.0, hub_height], deflection, model, turbine) == 0.0

        # test max loss at turbine (data from Jensen 1983)
        @test ff.wake_model([0.0, 0.0, hub_height], deflection, model, turbine) == (2. * aI)

        # test centerline loss 40 meters downstream (data from Jensen 1983)
        @test ff.wake_model([40., 0.0, hub_height], deflection, model, turbine) == centerloss40

        # test centerline loss 100 meters downstream (data from Jensen 1983)
        @test ff.wake_model([100., 0.0, hub_height], deflection, model, turbine) == centerloss100

        # test wake diameter 40 meters downstream (data from Jensen 1983)
        @test ff.wake_model([40., (alpha*40 + rotor_diameter/2.), hub_height], deflection, model, turbine) == centerloss40
        @test ff.wake_model([40., (alpha*40 + rotor_diameter/2. + 1E-12), hub_height], deflection, model, turbine) == 0.0
        @test ff.wake_model([40., -(alpha*40 + rotor_diameter/2.), hub_height], deflection, model, turbine) == centerloss40
        @test ff.wake_model([40., -(alpha*40 + rotor_diameter/2. + 1E-12), hub_height], deflection, model, turbine) == 0.0

        # test wake diameter 100 meters downstream (data from Jensen 1983)
        @test ff.wake_model([100., (alpha*100. + rotor_diameter/2.), hub_height], deflection, model, turbine) == centerloss100
        @test ff.wake_model([100., (alpha*100. + rotor_diameter/2. + 1E-12), hub_height], deflection, model, turbine) == 0.0
        @test ff.wake_model([100., -(alpha*100. + rotor_diameter/2.), hub_height], deflection, model, turbine) == centerloss100
        @test ff.wake_model([100., -(alpha*100. + rotor_diameter/2. + 1E-12), hub_height], deflection, model, turbine) == 0.0

    end
    
end
