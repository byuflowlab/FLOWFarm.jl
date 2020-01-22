import FlowFarm; const ff = FlowFarm
using Test

# @testset "Turbines" begin
#     coord = ff.Coord(0.0,0.0,0.0)
#     rotor_diameter = 126.4
#     hub_height = 90.0
#     aI = 1.0/3.0
#     yaw = 0.0
#     ct = 0.7 
#     ff.nrel5mw(coord, rotor_diameter, hub_height, aI, yaw, ct) = ff.Turbine(coord, rotor_diameter, hub_height, aI, yaw, ct)
# end

@testset "Wake Models" begin
    
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

    @testset "Jensen Cosine Model" begin

        coord = ff.Coord(0.0, 0.0, 0.0)
        rotor_diameter = 40.0
        hub_height = 90.0
        aI = 1.0/3.0
        yaw = 0.0
        ct = 0.7 #TODO handle ct and axial induction appropriately

        alpha = 0.1
        beta = 20.0*pi/180.0
        
        turbine = ff.Turbine(coord, rotor_diameter, hub_height, aI, yaw, ct)
        
        deflection = [0.0, 0.0]

        model = ff.JensenCosine(alpha, beta)

        centerloss40 = 1. - 4.35/8.1
        centerloss100 = 1. - 5.7/8.1

        d = rotor_diameter*0.5/(tan(model.beta))
        thetamax = 20.0*pi/180.0
        dy40 = tan(thetamax)*(40.0+d)
        dy100 = tan(thetamax)*(100.0+d)

        theta = 15.0*pi/180.
        dx = 40.0
        dy = tan(theta)*(dx+d)
        n = pi/beta
        ftheta = (1.0+cos(n*theta))/2.0
        loss40attheta = (2.0*aI)*((ftheta*rotor_diameter*0.5)/(rotor_diameter*0.5+alpha*dx))^2

        # test no loss upstream (data from Jensen 1983)
        @test ff.wake_model([-1E-12, 0.0, hub_height], deflection, model, turbine) == 0.0

        # test max loss at turbine (data from Jensen 1983)
        @test ff.wake_model([0.0, 0.0, hub_height], deflection, model, turbine) == (2. * aI)

        # test centerline loss 40 meters downstream (data from Jensen 1983)
        @test ff.wake_model([40., 0.0, hub_height], deflection, model, turbine) == centerloss40

        # test centerline loss 100 meters downstream (data from Jensen 1983)
        @test ff.wake_model([100., 0.0, hub_height], deflection, model, turbine) == centerloss100

        # test wake diameter 40 meters downstream (data from Jensen 1983)
        @test ff.wake_model([40., dy40, hub_height], deflection, model, turbine) == 0.0
        @test ff.wake_model([40., (dy40 + 1E-12), hub_height], deflection, model, turbine) == 0.0
        @test ff.wake_model([40., (dy40 - 1E1), hub_height], deflection, model, turbine) >= 0.0
        @test ff.wake_model([40., -(dy40), hub_height], deflection, model, turbine) == 0.0
        @test ff.wake_model([40., -(dy40 + 1E-12), hub_height], deflection, model, turbine) == 0.0
        @test ff.wake_model([40., -(dy40 - 1E1), hub_height], deflection, model, turbine) >= 0.0

        # test wake diameter 100 meters downstream (data from Jensen 1983)
        @test ff.wake_model([100., dy100, hub_height], deflection, model, turbine) == 0.0
        @test ff.wake_model([100., (dy100 + 1E-12), hub_height], deflection, model, turbine) == 0.0
        @test ff.wake_model([100., (dy100 - 1E1), hub_height], deflection, model, turbine) >= 0.0
        @test ff.wake_model([100., -(dy100), hub_height], deflection, model, turbine) == 0.0
        @test ff.wake_model([100., -(dy100 + 1E-12), hub_height], deflection, model, turbine) == 0.0
        @test ff.wake_model([100., -(dy100 - 1E1), hub_height], deflection, model, turbine) >= 0.0

        # test value at point in wake 40 m downstream and with theta=15 degrees
        @test ff.wake_model([40., dy, hub_height], deflection, model, turbine) == loss40attheta

    end
    
end
