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

@testset "Deflection Models" begin

    @testset "Gauss Yaw Deflection" begin

    coord = ff.Coord(0.000022, 0.0, 0.0) #[1] p. 509
    rotor_diameter = 0.15 #[1] p. 509
    hub_height = 0.125 #[1] p. 509
    yaw_20 = 20.0*pi/180.0
    ct = 0.82 # [1] fig. 2
    aI = 1.0/3.0

    k_star = 0.022 # [1]  p. 525
    turbulence_intensity = 0.1 #0.0875 #[1] p. 508 - this value is only specified to be less than 0.1
    horizontal_spread_rate = k_star
    vertical_spread_rate = k_star
    alpha_star = 2.32 #[1] p. 534
    beta_star = 0.154 #[1] p. 534

    turbine = ff.Turbine(coord, rotor_diameter, hub_height, aI, yaw_20, ct)
    model = ff.GaussYawDeflection(turbulence_intensity, horizontal_spread_rate, vertical_spread_rate, alpha_star, beta_star)

    dx4d = 4.0*rotor_diameter
    dy4d_20 = 0.2684659090909065*rotor_diameter # from [1] figure 21

    dx8d = 8.0*rotor_diameter
    dy8d_20 = 0.34090909090908905*rotor_diameter # from [1] figure 21
    
    # test deflection at 4D
    @test round(ff.deflection_model([dx4d, dy4d_20, hub_height], model, turbine), digits=2) == round(dy4d_20, digits=2)

    # test deflection at 8D
    @test round(ff.deflection_model([dx8d, dy8d_20, hub_height], model, turbine), digits=2) == round(dy8d_20, digits=2)

    end

end

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
    
    @testset "Gauss Yaw Model" begin
        # [1] based on data from Bastankhah and Porte-Agel 2016

        coord = ff.Coord(0.000022, 0.0, 0.0) #[1] p. 509
        rotor_diameter = 0.15 #[1] p. 509
        hub_height = 0.125 #[1] p. 509
        yaw = 0.0
        ct = 0.82 # [1] fig. 2

        k_star = 0.022 # [1]  p. 525
        turbulence_intensity = 0.1 #0.0875 #[1] p. 508 - this value is only specified to be less than 0.1
        horizontal_spread_rate = k_star
        vertical_spread_rate = k_star
        alpha_star = 2.32 #[1] p. 534
        beta_star = 0.154 #[1] p. 534
        
        turbine = ff.Turbine(coord, rotor_diameter, hub_height, 1.0/3.0, yaw, ct)
        
        deflection = [0.0, 0.0]

        model = ff.GaussYaw(turbulence_intensity, horizontal_spread_rate , vertical_spread_rate, alpha_star, beta_star)

        # data from Bastankhah and Porte-Agel 2016, figure 19
        yaw_20 = 20.0*pi/180.0
        ct_20 = 0.7378415935735881
        # x0_20 = 4.217687074829926*rotor_diameter
        x1_20 = 4.221216585981899*rotor_diameter 
        loss1_20 = 0.49396179751631175
        x2_20 = 6.014102294253845*rotor_diameter
        loss2_20 = 0.31515733529783185
        x3_20 = 7.987265838770787*rotor_diameter 
        loss3_20 = 0.22750736687013284 
        x8d_20 = 8.0*rotor_diameter
        loss8d_20 = .7203389830508229 # [1] figure 21

        yaw_0 = 0.0*pi/180.0
        # x0_0 = 4.112436115843271*rotor_diameter
        x1_0 = 4.0*rotor_diameter
        loss1_0 = 0.45806451612903154 # from figure 21
        # loss1_0 = 0.5922779922779922
        x2_0 = 6.00904977375566*rotor_diameter 
        # loss2_0 = 1.0 - 0.6033
        loss2_0 = 0.37606177606177627
        x3_0 = 8.00452488687783*rotor_diameter
        loss3_0 = 0.2710424710424715
        x4_0 = 10.000000000000002*rotor_diameter 
        loss4_0 = 0.20772200772200788

        # test no loss upstream
        @test ff.wake_model([-1E-12, 0.0, hub_height], deflection, model, turbine) == 0.0

        # test loss at x1 with no yaw
        @test round(ff.wake_model([x1_0, 0.0, hub_height], deflection, model, turbine), digits=1) == round(loss1_0, digits=1)

        deflection = [0.0, 0.0]
        # test loss at x2 with no yaw
        @test round(ff.wake_model([x2_0, 0.0, hub_height], deflection, model, turbine), digits=1) == round(loss2_0, digits=1)

        # test loss at x3 with no yaw
        @test round(ff.wake_model([x3_0, 0.0, hub_height], deflection, model, turbine), digits=1) == round(loss3_0, digits=1)

        # test loss at x4 with no yaw
        @test round(ff.wake_model([x4_0, 0.0, hub_height], deflection, model, turbine), digits=1) == round(loss4_0, digits=1)
        
        # # test centerline loss 40 meters downstream (data from Jensen 1983)
        # @test ff.wake_model([40., 0.0, hub_height], deflection, model, turbine) == centerloss40

        # # test centerline loss 100 meters downstream (data from Jensen 1983)
        # @test ff.wake_model([100., 0.0, hub_height], deflection, model, turbine) == centerloss100

        # # test wake diameter 40 meters downstream (data from Jensen 1983)
        # @test ff.wake_model([40., dy40, hub_height], deflection, model, turbine) == 0.0
        # @test ff.wake_model([40., (dy40 + 1E-12), hub_height], deflection, model, turbine) == 0.0
        # @test ff.wake_model([40., (dy40 - 1E1), hub_height], deflection, model, turbine) >= 0.0
        # @test ff.wake_model([40., -(dy40), hub_height], deflection, model, turbine) == 0.0
        # @test ff.wake_model([40., -(dy40 + 1E-12), hub_height], deflection, model, turbine) == 0.0
        # @test ff.wake_model([40., -(dy40 - 1E1), hub_height], deflection, model, turbine) >= 0.0

        # # test wake diameter 100 meters downstream (data from Jensen 1983)
        # @test ff.wake_model([100., dy100, hub_height], deflection, model, turbine) == 0.0
        # @test ff.wake_model([100., (dy100 + 1E-12), hub_height], deflection, model, turbine) == 0.0
        # @test ff.wake_model([100., (dy100 - 1E1), hub_height], deflection, model, turbine) >= 0.0
        # @test ff.wake_model([100., -(dy100), hub_height], deflection, model, turbine) == 0.0
        # @test ff.wake_model([100., -(dy100 + 1E-12), hub_height], deflection, model, turbine) == 0.0
        # @test ff.wake_model([100., -(dy100 - 1E1), hub_height], deflection, model, turbine) >= 0.0

        # # test value at point in wake 40 m downstream and with theta=15 degrees
        # @test ff.wake_model([40., dy, hub_height], deflection, model, turbine) == loss40attheta

    end
end
