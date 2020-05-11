using FlowFarm; const ff = FlowFarm
using Test
using DelimitedFiles
using LinearAlgebra
using PyPlot

@testset "AEP function" begin

    @testset "Test AEP" begin

        include("model_sets/model_set_3.jl")
        println(sum(windprobabilities))
        modelAEP = ff.calculate_aep(ms3, pd3,rotor_sample_points_y=rotor_points_y,rotor_sample_points_z=rotor_points_z)/1e9
        paperAEP = 1889.3
        println(modelAEP/paperAEP)
        @test modelAEP/paperAEP ≈ 1 atol=0.05
        end

end

@testset "utilities" begin

    @testset "wake overlap" begin

        turbine_y = 1000.0
        turbine_z = wake_center_z = 90.0
        rotor_diameter = 80.0
        wake_diameter = 80.0
        wake_center_y = 0.0

        # test no overlap
        overlap = ff.overlap_area_func(turbine_y, turbine_z, rotor_diameter, wake_center_y, 
            wake_center_z, wake_diameter)
        @test overlap == 0.0

        # test partial overlap
        turbine_y = 0.8079455*rotor_diameter/2.0
        overlap = ff.overlap_area_func(turbine_y, turbine_z, rotor_diameter, wake_center_y, 
            wake_center_z, wake_diameter)
        @test overlap ≈ 0.5*0.25*pi*rotor_diameter^2 atol=1E-4

        # test full overlap larger wake
        turbine_y = 0.0
        wake_diameter = 90.0
        overlap = ff.overlap_area_func(turbine_y, turbine_z, rotor_diameter, wake_center_y, 
            wake_center_z, wake_diameter)
        @test overlap ≈ 0.25*pi*rotor_diameter^2 atol=1E-6

        # test full overlap larger rotor
        turbine_y = 0.0
        wake_diameter = 70.0
        overlap = ff.overlap_area_func(turbine_y, turbine_z, rotor_diameter, wake_center_y, 
            wake_center_z, wake_diameter)
        @test overlap ≈ 0.25*pi*wake_diameter^2 atol=1E-6

    end

    @testset "smooth max" begin

        x = 1.0
        y = 2.0

        m = ff.smooth_max(x, y)
        @test m ≈ y atol=1E-4

        x = 1.99
        m = ff.smooth_max(x, y, s=400)
        @test m ≈ y atol=1E-4

        x = -4
        m = ff.smooth_max(x, y, s=4)
        @test m ≈ y atol=1E-6

    end

end

@testset "Optimization functions" begin

    @testset "Turbine spacing calculation" begin

        testing_x = [0,0,500,500,7481]
        testing_y = [0,500,0,500,-43891]
        turbine_separation = [500,500,707.1067811865476,44523.9850193129,707.1067811865476,500,45016.95505029189,500,44442.70741077775,44936.56909466943]

        @test ff.turbine_spacing(testing_x,testing_y) == turbine_separation

    end

    @testset "Circular boundary" begin

        center = [100,500]
        radius = 500
        testing_x = [100,100,100,100,-500,700]
        testing_y = [500,1100,-100,0,500,500]
        test_values = [-500,100,100,0,100,100]

        @test ff.circle_boundary(center,radius,testing_x,testing_y) == test_values

    end

    @testset "Polygon boundary" begin

    v = zeros(4,2)
    v[2,1] = 500
    v[3,1] = 500
    v[3,2] = 500
    v[4,2] = 500
    n = zeros(4,2)
    n[1,2] = -1
    n[2,1] = 1
    n[3,2] = 1
    n[4,1] = -1

    testing_x = [250,250,250,250,-100,600]
    testing_y = [250,-100,600,0,250,250]
    test_values = ones(6,4).*-250
    test_values[2,1] = 100
    test_values[2,3] = -600
    test_values[3,1] = -600
    test_values[3,3] = 100
    test_values[4,1] = 0
    test_values[4,3] = -500
    test_values[5,2] = -600
    test_values[5,4] = 100
    test_values[6,2] = 100
    test_values[6,4] = -600

    @test ff.windfarm_boundary(v,n,testing_x,testing_y) == test_values
    end

end


@testset "Thrust Coefficient Models" begin

    @testset "Constant Thrust Coefficent Model" begin

        ct = 0.7
        v0 = 8.0
        ct_model = ff.ThrustModelConstantCt(ct)

        @test ff.calculate_ct(v0, ct_model) == ct

    end

    @testset "Point based Thrust Coefficent Model" begin

        # load data
        data = readdlm("inputfiles/NREL5MWCPCT.txt",  ' ', skipstart=1)

        # extract velocity and ct points
        velpoints = data[:,1]
        ctpoints = data[:,3]

        # intialize ct model
        ct_model = ff.ThrustModelCtPoints(velpoints, ctpoints)

        # low velocity
        v0 = 1.0
        ct0 = ctpoints[1]
        @test ff.calculate_ct(v0, ct_model) == ct0

        # mid velocity
        v0 = 0.5*(6.925399017365052146 + 7.056245651277220254)
        ct0 = 0.5*(7.930937230057298892e-01 + 7.884019134159416797e-01)
        @test ff.calculate_ct(v0, ct_model) == ct0

        # high velocity
        v0 = 30.0
        ct0 = ctpoints[end]
        @test ff.calculate_ct(v0, ct_model) == ct0


    end

end

@testset "Power Models" begin
    @testset "calculate_power_from_cp" begin
        generator_efficiency = 0.944
        air_density = 1.1716
        rotor_area = pi*80.0^2/4
        cp = 0.8
        v0 = 12.0

        p = ff.calculate_power_from_cp(generator_efficiency, air_density, rotor_area, cp, v0)
        @test p ≈ 3.8425979093271587e6 atol=1E-6

    end

    @testset "calculate_power() PowerModelConstantCP" begin
        generator_efficiency = 0.944
        air_density = 1.1716
        rotor_area = 0.25*pi*80.0^2
        cp = 0.8
        v0 = 8.0
        cut_in_speed = 4.  # m/s
        cut_out_speed = 25.  # m/s
        rated_speed = 16.  # m/s
        rated_power = 2.0E6  # W

        power_model = ff.PowerModelConstantCp(cp)

        p = ff.calculate_power(generator_efficiency, air_density, rotor_area, v0, cut_in_speed, rated_speed, cut_out_speed, rated_power, power_model)
        @test p ≈ 0.5*rotor_area*cp*air_density*generator_efficiency*v0^3 atol=1E-6

    end

    @testset "calculate_power() PowerModelCpPoints" begin
        generator_efficiency = 0.944
        air_density = 1.1716
        rotor_area = pi*80.0^2/4
        cut_in_speed = 4.  # m/s
        cut_out_speed = 25.  # m/s
        rated_speed = 16.  # m/s
        rated_power = 2.0E6  # W

        # load data
        data = readdlm("inputfiles/NREL5MWCPCT.txt",  ' ', skipstart=1)

        # vel and cp for test based on input cp curve
        v0 = 0.5*(8.495558624311073004 + 8.626405258223240224)
        cp0 = 0.5*(4.631607703567138801e-01 + 4.631607703567138801e-01)

        # calculate expected power out
        power = 0.5*cp0*air_density*rotor_area*generator_efficiency*v0^3

        # extract velocity and cp points
        velpoints = data[:,1]
        cppoints = data[:,2]

        # intialize power model struct
        power_model = ff.PowerModelCpPoints(velpoints, cppoints)

        # calculated power and test
        p = ff.calculate_power(generator_efficiency, air_density, rotor_area, v0, cut_in_speed, rated_speed, cut_out_speed, rated_power, power_model)
        @test p ≈ power atol=1E-6

    end

    @testset "calculate_power() PowerModelPowerPoints" begin
        generator_efficiency = 0.944
        air_density = 1.1716
        rotor_area = pi*80.0^2/4
        cut_in_speed = 4.  # m/s
        cut_out_speed = 25.  # m/s
        rated_speed = 16.  # m/s
        rated_power = 2.0E6  # W

        # vel and power for test based on input power curve
        v0 = 0.5*(6.9778601570742005 + 7.469669440862736)
        p0 = 0.5*(0.4665924276169269 + 0.5768374164810695)*1E6

        # load data
        data = readdlm("inputfiles/niayifar_vestas_v80_power_curve_observed.txt",  ',', skipstart=1)

        # extract velocity and cp points
        velpoints = data[:,1]
        powpoints = data[:,2]*1E6

        # intialize power model struct
        power_model = ff.PowerModelPowerPoints(velpoints, powpoints)

        # calc power
        p = ff.calculate_power(generator_efficiency, air_density, rotor_area, v0, cut_in_speed, rated_speed, cut_out_speed, rated_power, power_model)

        # test
        @test p ≈ p0 atol=1E-6

    end

    @testset "calculate_turbine_power() PowerModelConstantCp" begin

        include("./model_sets/model_set_2.jl")

        # test below cut in
        windfarmstate.turbine_inflow_velcities[1] = 1.0
        p = ff.calculate_turbine_power(1, turbine1, windfarmstate, windresource)
        @test p ≈ 0.0 atol=1E-6

        # region 2
        v0 = windfarmstate.turbine_inflow_velcities[1] = 8.0
        p = ff.calculate_turbine_power(1, turbine1, windfarmstate, windresource)
        rotor_area = pi*0.25*rotor_diameter^2
        @test p ≈ 0.5*cp*air_density*rotor_area*generator_efficiency*v0^3 atol=1E-6

        # above rated
        windfarmstate.turbine_inflow_velcities[1] = 20.0
        p = ff.calculate_turbine_power(1, turbine1, windfarmstate, windresource)
        @test p ≈ turbine1.rated_power[1] atol=1E-6

        # above cut out
        windfarmstate.turbine_inflow_velcities[1] = 30.0
        p = ff.calculate_turbine_power(1, turbine1, windfarmstate, windresource)
        @test p ≈ 0.0 atol=1E-6

    end

end


@testset "Wind Shear Models" begin

    @testset "Power Law Wind Shear" begin
        # TODO base this tests on something more concrete

        point_velocity_no_shear = 8.0
        reference_height = 80.0
        ground_height = 0.0
        shear_exp = 0.15

        model = ff.PowerLawWindShear(shear_exp)
        # test at reference height
        loc =[0.0, 0.0, 80.0]
        @test ff.adjust_for_wind_shear(loc, point_velocity_no_shear, reference_height, ground_height, model) == point_velocity_no_shear

        # test at ground height
        loc =[0.0, 0.0, 0.0]
        @test ff.adjust_for_wind_shear(loc, point_velocity_no_shear, reference_height, ground_height, model) == 0.0

        # test below ground height
        loc =[0.0, 0.0, -10.0]
        @test ff.adjust_for_wind_shear(loc, point_velocity_no_shear, reference_height, ground_height, model) == 0.0

        # test at 40 meters
        loc =[0.0, 0.0, 40.0]
        u = point_velocity_no_shear*((loc[3] - ground_height)/(reference_height-ground_height))^shear_exp
        @test ff.adjust_for_wind_shear(loc, point_velocity_no_shear, reference_height, ground_height, model) == u

        # test at 10 meters
        loc =[0.0, 0.0, 10.0]
        u = point_velocity_no_shear*((loc[3] - ground_height)/(reference_height-ground_height))^shear_exp
        @test ff.adjust_for_wind_shear(loc, point_velocity_no_shear, reference_height, ground_height, model) == u

    end

end

@testset "Wake Combination Models" begin

    old_deficit_sum = 0.3
    wind_speed = 8.0
    deltav = 0.2
    turb_inflow = 7.5

    @testset "Linear Freestream Superposition" begin
        model = ff.LinearFreestreamSuperposition()
        result = old_deficit_sum + wind_speed*deltav
        @test ff.wake_combination_model(deltav, wind_speed, turb_inflow, old_deficit_sum, model) == result
    end

    @testset "Sum of Squares Freestream Superposition" begin
        model = ff.SumOfSquaresFreestreamSuperposition()
        result = sqrt(old_deficit_sum^2 + (wind_speed*deltav)^2)
        @test ff.wake_combination_model(deltav, wind_speed, turb_inflow, old_deficit_sum, model) == result
    end

    @testset "Sum of Squares Local Velocity Superposition" begin
        model = ff.SumOfSquaresLocalVelocitySuperposition()
        result = sqrt(old_deficit_sum^2 + (turb_inflow*deltav)^2)
        @test ff.wake_combination_model(deltav, wind_speed, turb_inflow, old_deficit_sum, model) == result
    end

    @testset "Linear Local Velocity Superposition" begin
    model = ff.LinearLocalVelocitySuperposition()
    result = old_deficit_sum + turb_inflow*deltav
    @test ff.wake_combination_model(deltav, wind_speed, turb_inflow, old_deficit_sum, model) == result
    end

end

@testset "Wake Deflection Models" begin

    @testset "Jiminez Yaw Deflection" begin

        # [1] Jiminez 2010

        # The following are CFD results, not model results, so testing is a little sketchy.
        # data for ct=0.8, x/d = 2.5 from [1] fig. 8
        # z/d                   u/u0
        # 0.005141388174807027, 0.5123919308357348 # yaw = 0
        # -0.15938303341902316, 0.5089337175792507 # yaw = 10
        # -0.23136246786632375, 0.5677233429394812 # yaw = 20
        # -0.30334190231362435, 0.6576368876080692 # yaw = 30

        # data for ct=0.8, x/d = 5.5
        # z/d                   u/u0
        # 0.0051413881748061385, 0.7542028985507248 # yaw = 0
        # -0.2107969151670961, 0.7460869565217392 # yaw = 10
        # -0.40616966580976843, 0.7657971014492755 # yaw = 20
        # -0.4267352185089983, 0.8191304347826089 # yaw = 30

        # data for ct=0.8, x/d = 8.0
        # z/d                   u/u0
        # 0.010309278350515427, 0.8462457541266909 # yaw = 0
        # -0.2680412371134011, 0.8519456528216436 # yaw = 10
        # -0.5257731958762868, 0.8588075800011918 # yaw = 20
        # -0.5670103092783503, 0.8923216733210179 # yaw = 30

        rotor_diameter = 0.15
        hub_height = 0.125
        yaw_20 = 20.0*pi/180.0
        ct = 0.8 # [1] fig. 8
        cp = 0.8
        generator_efficiency = 0.944

        wind_farm_state_id = 1
        turbine_x = [0.0]
        turbine_y = [0.0]
        turbine_z = [0.0]
        turbine_yaw = [yaw_20]
        turbine_ct = [ct]
        turbine_ai = [1.0/3.0]
        sorted_turbine_index = [1]
        turbine_inflow_velcity = [8.0]
        turbine_id = 1
        turbine_definition_id = 1
        cut_in_speed = 0.0
        cut_out_speed = 25.0
        rated_speed = 12.0
        rated_power = 1.0176371581904552e6
        ambient_ti = 0.1

        windfarmstate = ff.SingleWindFarmState(wind_farm_state_id, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai, sorted_turbine_index, turbine_inflow_velcity, [0.0], [ambient_ti])

        k_star = 0.07 # adjusted to match experimental data. #TODO improve tests with model results
        horizontal_spread_rate = k_star

        ct_model = ff.ThrustModelConstantCt(ct)
        power_model = ff.PowerModelConstantCp(cp)

        turbine1 = ff.TurbineDefinition(turbine_definition_id, [rotor_diameter], [hub_height], [cut_in_speed], [rated_speed], [cut_out_speed], [rated_power], [generator_efficiency], [ct_model], [power_model])
        model = ff.JiminezYawDeflection(horizontal_spread_rate)

        dx2p5d = 2.5*rotor_diameter
        dy2p5d_y20 = 0.23136246786632375*rotor_diameter # from [1] figure 21

        dx5p5d = 5.5*rotor_diameter
        dy5p5d_y20 = 0.40616966580976843*rotor_diameter # from [1] figure 21

        dx8d = 8.0*rotor_diameter
        dy8d_20 = 0.5257731958762868*rotor_diameter # from [1] figure 21

        # test deflection at 2.5D with yaw 20 deg
        @test round(ff.wake_deflection_model([dx2p5d, 0.0, hub_height], turbine_id, turbine1, model, windfarmstate), digits=2) == round(dy2p5d_y20, digits=2)

        # test deflection at 5.5D with yaw 20 deg
        @test round(ff.wake_deflection_model([dx5p5d, 0.0, hub_height], turbine_id, turbine1, model, windfarmstate), digits=2) == round(dy5p5d_y20, digits=2)

        # test deflection at 8D with yaw 20 deg
        # @test round(ff.wake_deflection_model([dx8d, 0.0, hub_height], model, turbine), digits=2) == round(dy8d_20, digits=2)

    end

    @testset "Gauss Yaw Deflection" begin

        atol = 0.005

        rotor_diameter = 0.15 #[1] p. 509
        hub_height = 0.125 #[1] p. 509
        yaw_20 = 20.0*pi/180.0
        ct = 0.82 # [1] fig. 8
        cp = 0.8
        generator_efficiency = 0.944

        wind_farm_state_id = 1
        turbine_x = [0.0]
        turbine_y = [0.0]
        turbine_z = [0.000022] #[1] p. 509
        turbine_yaw = [yaw_20]
        turbine_ct = [ct]
        turbine_ai = [1.0/3.0]
        sorted_turbine_index = [1]
        turbine_inflow_velcity = [8.0]
        turbine_id = 1
        turbine_definition_id = 1
        cut_in_speed = 0.0
        cut_out_speed = 25.0
        rated_speed = 12.0
        rated_power = 1.0176371581904552e6
        ambient_ti = 0.1

        windfarmstate = ff.SingleWindFarmState(wind_farm_state_id, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai, sorted_turbine_index, turbine_inflow_velcity, [0.0], [ambient_ti])
        ct_model = ff.ThrustModelConstantCt(ct)
        power_model = ff.PowerModelConstantCp([cp])

        turbine_definition = ff.TurbineDefinition(turbine_definition_id, [rotor_diameter], [hub_height], [cut_in_speed], [rated_speed], [cut_out_speed], [rated_power], [generator_efficiency], [ct_model], [power_model])

        k_star = 0.022 # [1]  p. 525
        turbulence_intensity = 0.1 #0.0875 #[1] p. 508 - this value is only specified to be less than 0.1
        horizontal_spread_rate = k_star
        vertical_spread_rate = k_star
        alpha_star = 2.32 #[1] p. 534
        beta_star = 0.154 #[1] p. 534

        model = ff.GaussYawDeflection(horizontal_spread_rate, vertical_spread_rate, alpha_star, beta_star)

        dx4d = 4.0*rotor_diameter
        dy4d_20 = 0.2684659090909065*rotor_diameter # from [1] figure 21

        dx8d = 8.0*rotor_diameter
        dy8d_20 = 0.34090909090908905*rotor_diameter # from [1] figure 21

        # test deflection at 4D with yaw 20 deg
        @test ff.wake_deflection_model([dx4d, dy4d_20, hub_height], turbine_id, turbine_definition, model, windfarmstate) ≈ dy4d_20 atol=atol

        # test deflection at 8D with yaw 20 deg
        @test ff.wake_deflection_model([dx8d, dy8d_20, hub_height], turbine_id, turbine_definition, model, windfarmstate) ≈ dy8d_20 atol=atol

    end

end

@testset "Wake Deficit Models" begin

    @testset "Jensen Top Hat Model" begin

        include("./model_sets/model_set_1.jl")

        turbine_id = 1
        turbine_definition = turbine1

        deflection = [0.0, 0.0]

        centerloss40 = 1. - 4.35/8.1
        centerloss100 = 1. - 5.7/8.1


        # test no loss upstream (data from Jensen 1983)
        @test ff.wake_deficit_model([-1E-12, 0.0, hub_height], deflection, turbine_id, turbine_definition, wakedeficitmodel, windfarmstate) == 0.0

        # test max loss at turbine (data from Jensen 1983)
        @test ff.wake_deficit_model([0.0, 0.0, hub_height], deflection, turbine_id, turbine_definition, wakedeficitmodel, windfarmstate) == (2. * 1/3.0)

        # test centerline loss 40 meters downstream (data from Jensen 1983)
        @test ff.wake_deficit_model([40., 0.0, hub_height], deflection, turbine_id, turbine_definition, wakedeficitmodel, windfarmstate) == centerloss40

        # test centerline loss 100 meters downstream (data from Jensen 1983)
        @test ff.wake_deficit_model([100., 0.0, hub_height], deflection, turbine_id, turbine_definition, wakedeficitmodel, windfarmstate) == centerloss100

        # test wake diameter 40 meters downstream (data from Jensen 1983)
        @test ff.wake_deficit_model([40., (alpha*40 + rotor_diameter/2.), hub_height], deflection, turbine_id, turbine_definition, wakedeficitmodel, windfarmstate) == centerloss40
        @test ff.wake_deficit_model([40., (alpha*40 + rotor_diameter/2. + 1E-12), hub_height], deflection, turbine_id, turbine_definition, wakedeficitmodel, windfarmstate) == 0.0
        @test ff.wake_deficit_model([40., -(alpha*40 + rotor_diameter/2.), hub_height], deflection, turbine_id, turbine_definition, wakedeficitmodel, windfarmstate) == centerloss40
        @test ff.wake_deficit_model([40., -(alpha*40 + rotor_diameter/2. + 1E-12), hub_height], deflection, turbine_id, turbine_definition, wakedeficitmodel, windfarmstate) == 0.0

        # test wake diameter 100 meters downstream (data from Jensen 1983)
        @test ff.wake_deficit_model([100., (alpha*100. + rotor_diameter/2.), hub_height], deflection, turbine_id, turbine_definition, wakedeficitmodel, windfarmstate) == centerloss100
        @test ff.wake_deficit_model([100., (alpha*100. + rotor_diameter/2. + 1E-12), hub_height], deflection, turbine_id, turbine_definition, wakedeficitmodel, windfarmstate) == 0.0
        @test ff.wake_deficit_model([100., -(alpha*100. + rotor_diameter/2.), hub_height], deflection, turbine_id, turbine_definition, wakedeficitmodel, windfarmstate) == centerloss100
        @test ff.wake_deficit_model([100., -(alpha*100. + rotor_diameter/2. + 1E-12), hub_height], deflection, turbine_id, turbine_definition, wakedeficitmodel, windfarmstate) == 0.0
    end

    @testset "Jensen Cosine Model" begin

        rotor_diameter = 40.0
        hub_height = 90.0
        cut_in_speed = 0.0
        cut_out_speed = 25.0
        rated_speed = 12.0
        rated_power = 1.0176371581904552e6
        aI = 1.0/3.0
        yaw = 0.0
        ct = 0.7
        cp = 0.8
        generator_efficiency = 0.944

        wind_farm_state_id = 1
        turbine_x = [0.0]
        turbine_y = [0.0]
        turbine_z = [0.0]
        turbine_yaw = [yaw]
        turbine_ct = [ct]
        turbine_ai = [1.0/3.0]
        sorted_turbine_index = [1]
        turbine_inflow_velcity = [8.0]
        ambient_ti = 0.1
        turbine_id = 1
        turbine_definition_id = 1


        alpha = 0.1
        beta = 20.0*pi/180.0


        deflection = [0.0, 0.0]

        windfarmstate = ff.SingleWindFarmState(wind_farm_state_id, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai, sorted_turbine_index, turbine_inflow_velcity, [0.0], [ambient_ti])
        ct_model = ff.ThrustModelConstantCt(ct)
        power_model = ff.PowerModelConstantCp([cp])

        turbine_definition = ff.TurbineDefinition(1, [rotor_diameter], [hub_height], [cut_in_speed], [rated_speed], [cut_out_speed], [rated_power], [generator_efficiency], [ct_model], [power_model])
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
        @test ff.wake_deficit_model([-1E-12, 0.0, hub_height], deflection, turbine_id, turbine_definition, model, windfarmstate) == 0.0

        # test max loss at turbine (data from Jensen 1983)
        @test ff.wake_deficit_model([0.0, 0.0, hub_height], deflection, turbine_id, turbine_definition, model, windfarmstate) == (2. * aI)

        # test centerline loss 40 meters downstream (data from Jensen 1983)
        @test ff.wake_deficit_model([40., 0.0, hub_height], deflection, turbine_id, turbine_definition, model, windfarmstate) == centerloss40

        # test centerline loss 100 meters downstream (data from Jensen 1983)
        @test ff.wake_deficit_model([100., 0.0, hub_height], deflection, turbine_id, turbine_definition, model, windfarmstate) == centerloss100

        # test wake diameter 40 meters downstream (data from Jensen 1983)
        @test ff.wake_deficit_model([40., dy40, hub_height], deflection, turbine_id, turbine_definition, model, windfarmstate) == 0.0
        @test ff.wake_deficit_model([40., (dy40 + 1E-12), hub_height], deflection, turbine_id, turbine_definition, model, windfarmstate) == 0.0
        @test ff.wake_deficit_model([40., (dy40 - 1E1), hub_height], deflection, turbine_id, turbine_definition, model, windfarmstate) >= 0.0
        @test ff.wake_deficit_model([40., -(dy40), hub_height], deflection, turbine_id, turbine_definition, model, windfarmstate) == 0.0
        @test ff.wake_deficit_model([40., -(dy40 + 1E-12), hub_height], deflection, turbine_id, turbine_definition, model, windfarmstate) == 0.0
        @test ff.wake_deficit_model([40., -(dy40 - 1E1), hub_height], deflection, turbine_id, turbine_definition, model, windfarmstate) >= 0.0

        # test wake diameter 100 meters downstream (data from Jensen 1983)
        @test ff.wake_deficit_model([100., dy100, hub_height], deflection, turbine_id, turbine_definition, model, windfarmstate) == 0.0
        @test ff.wake_deficit_model([100., (dy100 + 1E-12), hub_height], deflection, turbine_id, turbine_definition, model, windfarmstate) == 0.0
        @test ff.wake_deficit_model([100., (dy100 - 1E1), hub_height], deflection, turbine_id, turbine_definition, model, windfarmstate) >= 0.0
        @test ff.wake_deficit_model([100., -(dy100), hub_height], deflection, turbine_id, turbine_definition, model, windfarmstate) == 0.0
        @test ff.wake_deficit_model([100., -(dy100 + 1E-12), hub_height], deflection, turbine_id, turbine_definition, model, windfarmstate) == 0.0
        @test ff.wake_deficit_model([100., -(dy100 - 1E1), hub_height], deflection, turbine_id, turbine_definition, model, windfarmstate) >= 0.0

        # test value at point in wake 40 m downstream and with theta=15 degrees
        @test ff.wake_deficit_model([40., dy, hub_height], deflection, turbine_id, turbine_definition, model, windfarmstate) == loss40attheta

        # test with wec
        model.wec_factor[1] = 1.0
        loss0 = ff.wake_deficit_model([100.0, 50.0, hub_height], deflection, turbine_id, turbine_definition, model, windfarmstate)
        model.wec_factor[1] = 3.0
        loss1 = ff.wake_deficit_model([100.0, 50.0, hub_height], deflection, turbine_id, turbine_definition, model, windfarmstate)
        @test loss0 < loss1

        model.wec_factor[1] = 1.0
        loss0 = ff.wake_deficit_model([100.0, 0.0, hub_height], deflection, turbine_id, turbine_definition, model, windfarmstate)
        model.wec_factor[1] = 3.0
        loss1 = ff.wake_deficit_model([100.0, 0.0, hub_height], deflection, turbine_id, turbine_definition, model, windfarmstate)
        @test loss0 == loss1

    end

    @testset "Gauss Yaw Model" begin
        rtol = 0.1
        # [1] based on data from Bastankhah and Porte-Agel 2016

        turbine_x = [0.0]
        turbine_y = [0.0]
        turbine_z = [0.000022] #[1] p. 509
        rotor_diameter = 0.15 #[1] p. 509
        hub_height = 0.125 #[1] p. 509
        cut_in_speed = 0.0
        cut_out_speed = 25.0
        rated_speed = 12.0
        rated_power = 1.0176371581904552e6
        yaw = 0.0
        ct = 0.82 # [1] fig. 2
        cp = 0.8
        generator_efficiency = 0.944

        k_star = 0.022 # [1]  p. 525
        turbulence_intensity = 0.1 #0.0875 #[1] p. 508 - this value is only specified to be less than 0.1
        horizontal_spread_rate = k_star
        vertical_spread_rate = k_star
        alpha_star = 2.32 #[1] p. 534
        beta_star = 0.154 #[1] p. 534

        wind_farm_state_id = 1

        turbine_yaw = [yaw]
        turbine_ct = [ct]
        turbine_ai = [1.0/3.0]
        sorted_turbine_index = [1]
        turbine_inflow_velcity = [8.0]
        ambient_ti = 0.1
        turbine_id = 1
        turbine_definition_id = 1

        deflection = [0.0, 0.0]

        windfarmstate = ff.SingleWindFarmState(wind_farm_state_id, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai, sorted_turbine_index, turbine_inflow_velcity, [0.0], [ambient_ti])
        ct_model = ff.ThrustModelConstantCt(ct)
        power_model = ff.PowerModelConstantCp([cp])

        turbine_definition = ff.TurbineDefinition(1, [rotor_diameter], [hub_height], [cut_in_speed], [rated_speed], [cut_out_speed], [rated_power], [generator_efficiency], [ct_model], [power_model])

        model = ff.GaussYaw(horizontal_spread_rate , vertical_spread_rate, alpha_star, beta_star)

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
        @test ff.wake_deficit_model([-1E-12, 0.0, hub_height], deflection, turbine_id, turbine_definition, model, windfarmstate) == 0.0

        # test loss at x1 with no yaw
        @test ff.wake_deficit_model([x1_0, 0.0, hub_height], deflection, turbine_id, turbine_definition, model, windfarmstate) ≈ loss1_0 rtol=rtol

        deflection = [0.0, 0.0]
        # test loss at x2 with no yaw
        @test ff.wake_deficit_model([x2_0, 0.0, hub_height], deflection, turbine_id, turbine_definition, model, windfarmstate)  ≈ loss2_0 rtol=rtol

        # test loss at x3 with no yaw
        @test ff.wake_deficit_model([x3_0, 0.0, hub_height], deflection, turbine_id, turbine_definition, model, windfarmstate) ≈ loss3_0 rtol=rtol

        # test loss at x4 with no yaw
        @test ff.wake_deficit_model([x4_0, 0.0, hub_height], deflection, turbine_id, turbine_definition, model, windfarmstate)  ≈ loss4_0 rtol=rtol

        # test with wec
        model.wec_factor[1] = 1.0
        loss0 = ff.wake_deficit_model([x4_0, 0.0, hub_height], deflection, turbine_id, turbine_definition, model, windfarmstate)
        model.wec_factor[1] = 3.0
        loss1 = ff.wake_deficit_model([x4_0, 0.0, hub_height], deflection, turbine_id, turbine_definition, model, windfarmstate)
        @test loss0 ≈ loss1 rtol=1E-6

        model.wec_factor[1] = 1.0
        loss0 = ff.wake_deficit_model([x4_0, 2.0*rotor_diameter, hub_height], deflection, turbine_id, turbine_definition, model, windfarmstate)
        model.wec_factor[1] = 3.0
        loss1 = ff.wake_deficit_model([x4_0, 2.0*rotor_diameter, hub_height], deflection, turbine_id, turbine_definition, model, windfarmstate)
        @test loss1 > 3.0*loss0

        # # test centerline loss 40 meters downstream (data from Jensen 1983)
        # @test ff.wake_deficit_model([40., 0.0, hub_height], deflection, model, turbine) == centerloss40

        # # test centerline loss 100 meters downstream (data from Jensen 1983)
        # @test ff.wake_deficit_model([100., 0.0, hub_height], deflection, model, turbine) == centerloss100

        # # test wake diameter 40 meters downstream (data from Jensen 1983)
        # @test ff.wake_deficit_model([40., dy40, hub_height], deflection, model, turbine) == 0.0
        # @test ff.wake_deficit_model([40., (dy40 + 1E-12), hub_height], deflection, model, turbine) == 0.0
        # @test ff.wake_deficit_model([40., (dy40 - 1E1), hub_height], deflection, model, turbine) >= 0.0
        # @test ff.wake_deficit_model([40., -(dy40), hub_height], deflection, model, turbine) == 0.0
        # @test ff.wake_deficit_model([40., -(dy40 + 1E-12), hub_height], deflection, model, turbine) == 0.0
        # @test ff.wake_deficit_model([40., -(dy40 - 1E1), hub_height], deflection, model, turbine) >= 0.0

        # # test wake diameter 100 meters downstream (data from Jensen 1983)
        # @test ff.wake_deficit_model([100., dy100, hub_height], deflection, model, turbine) == 0.0
        # @test ff.wake_deficit_model([100., (dy100 + 1E-12), hub_height], deflection, model, turbine) == 0.0
        # @test ff.wake_deficit_model([100., (dy100 - 1E1), hub_height], deflection, model, turbine) >= 0.0
        # @test ff.wake_deficit_model([100., -(dy100), hub_height], deflection, model, turbine) == 0.0
        # @test ff.wake_deficit_model([100., -(dy100 + 1E-12), hub_height], deflection, model, turbine) == 0.0
        # @test ff.wake_deficit_model([100., -(dy100 - 1E1), hub_height], deflection, model, turbine) >= 0.0

        # # test value at point in wake 40 m downstream and with theta=15 degrees
        # @test ff.wake_deficit_model([40., dy, hub_height], deflection, model, turbine) == loss40attheta

    end

    @testset "Gauss Yaw Model Variable Spread" begin
        rtol = 0.1
        # [1] Bastankhah and Porte-Agel 2016
        # [2] Niayifar and Porte-Agel 2016

        turbine_x = [0.0]
        turbine_y = [0.0]
        turbine_z = [0.000022] #[1] p. 509
        rotor_diameter = 0.15 #[1] p. 509
        hub_height = 0.125 #[1] p. 509
        cut_in_speed = 0.0
        cut_out_speed = 25.0
        rated_speed = 12.0
        rated_power = 1.0176371581904552e6
        yaw = 0.0
        ct = 0.82 # [1] fig. 2
        generator_efficiency = 0.944

        k_star = 0.022 # [1]  p. 525
        turbulence_intensity = 0.1 #0.0875 #[1] p. 508 - this value is only specified to be less than 0.1
        alpha_star = 2.32 #[1] p. 534
        beta_star = 0.154 #[1] p. 534

        wind_farm_state_id = 1

        turbine_yaw = [yaw]
        turbine_ct = [ct]
        turbine_ai = [1.0/3.0]
        sorted_turbine_index = [1]
        turbine_inflow_velcity = [8.0]
        ambient_ti = (k_star - 0.003678)/0.3837
        println(ambient_ti)
        turbine_id = 1
        turbine_definition_id = 1

        deflection = [0.0, 0.0]

        windfarmstate = ff.SingleWindFarmState(wind_farm_state_id, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai, sorted_turbine_index, turbine_inflow_velcity, [0.0], [ambient_ti])
        ct_model = ff.ThrustModelConstantCt(ct)
        power_model = ff.PowerModelConstantCp([cp])

        turbine_definition = ff.TurbineDefinition(1, [rotor_diameter], [hub_height], [cut_in_speed], [rated_speed], [cut_out_speed], [rated_power], [generator_efficiency], [ct_model], [power_model])

        model = ff.GaussYawVariableSpread(alpha_star, beta_star)

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
        loss1_0 = 0.5757358817491954 # from figure 21
        # loss1_0 = 0.5922779922779922
        x2_0 = 6.00904977375566*rotor_diameter
        # loss2_0 = 1.0 - 0.6033
        loss2_0 = 0.4805792772689106
        x3_0 = 8.00452488687783*rotor_diameter
        loss3_0 = 0.3558433011696253
        x4_0 = 10.000000000000002*rotor_diameter
        loss4_0 = 0.27837504998473545

        # test no loss upstream
        @test ff.wake_deficit_model([-1E-12, 0.0, hub_height], deflection, turbine_id, turbine_definition, model, windfarmstate) == 0.0

        # test loss at x1 with no yaw
        @test ff.wake_deficit_model([x1_0, 0.0, hub_height], deflection, turbine_id, turbine_definition, model, windfarmstate) ≈ loss1_0 rtol=rtol

        deflection = [0.0, 0.0]
        # test loss at x2 with no yaw
        @test ff.wake_deficit_model([x2_0, 0.0, hub_height], deflection, turbine_id, turbine_definition, model, windfarmstate)  ≈ loss2_0 rtol=rtol

        # test loss at x3 with no yaw
        @test ff.wake_deficit_model([x3_0, 0.0, hub_height], deflection, turbine_id, turbine_definition, model, windfarmstate) ≈ loss3_0 rtol=rtol

        # test loss at x4 with no yaw
        @test ff.wake_deficit_model([x4_0, 0.0, hub_height], deflection, turbine_id, turbine_definition, model, windfarmstate)  ≈ loss4_0 rtol=rtol

        # test with wec
        model.wec_factor[1] = 1.0
        loss0 = ff.wake_deficit_model([x4_0, 0.0, hub_height], deflection, turbine_id, turbine_definition, model, windfarmstate)
        model.wec_factor[1] = 3.0
        loss1 = ff.wake_deficit_model([x4_0, 0.0, hub_height], deflection, turbine_id, turbine_definition, model, windfarmstate)
        @test loss0 < loss1

    end
end

@testset "Local Turbulence Intensity Models" begin

    @testset "Niayifar wake spread based on ti" begin

        ti = 0.077

        k = ff._k_star_func(ti)

        @test k == 0.3837*ti + 0.003678

    end

    @testset "Niayifar Added TI Function" begin

        tol = 1E-2
        yaw = 0.0
        ct = 0.8
        alpha = 2.32
        beta = 0.154
        ky = 0.022
        kz = 0.022
        wind_speed = 8.0

        ti = 0.077
        x = 560.0
        rotor_diameter = 80.0
        deltay = 0.0
        wake_height = 70.0
        turbine_height = 70.0
        sm_smoothing = 700.0

        ti_area_ratio_in = 0.0
        ti_dst_in = 0.0
        ti_ust = 0.077

        ti, ti_ratio = ff._niayifar_added_ti_function(x, rotor_diameter, rotor_diameter, wake_height, turbine_height, ct, ky, deltay, ti, ti_ust, ti_dst_in, ti_area_ratio_in; s=700.0)
        
        @test ti ≈ 0.1476 atol=tol

    end

    @testset "Local TI Model Max TI Ratio" begin

        atol = 1E-2

        # load model set
        include("./model_sets/model_set_4.jl")

        # calculate turbine inflow velocities
        ff.turbine_velocities_one_direction!(rotor_points_y, rotor_points_z, ms4, pd4)

        # load horns rev ti ata
        data = readdlm("inputfiles/horns_rev_ti_by_row_niayifar.txt", ',', skipstart=1)

        # freestream
        ti_dst = ff.calculate_local_ti(ambient_ti, windfarm, windfarmstate, localtimodel, turbine_id=(1+ 4*10), tol=1E-6)
        @test ti_dst  == data[1,2] 

        # row 2
        ti_dst = ff.calculate_local_ti(ambient_ti, windfarm, windfarmstate, localtimodel, turbine_id=(2+ 4*10), tol=1E-6)
        @test ti_dst  ≈ data[2,2] atol=atol

        # row 3
        ti_dst = ff.calculate_local_ti(ambient_ti, windfarm, windfarmstate, localtimodel, turbine_id=(3+ 4*10), tol=1E-6)
        @test ti_dst  ≈ data[3,2] atol=atol

        # row 4
        ti_dst = ff.calculate_local_ti(ambient_ti, windfarm, windfarmstate, localtimodel, turbine_id=(4+ 4*10), tol=1E-6)
        @test ti_dst  ≈ data[4,2] atol=atol

        # row 5
        ti_dst = ff.calculate_local_ti(ambient_ti, windfarm, windfarmstate, localtimodel, turbine_id=(5+ 4*10), tol=1E-6)
        @test ti_dst  ≈ data[5,2] atol=atol

        # row 6
        ti_dst = ff.calculate_local_ti(ambient_ti, windfarm, windfarmstate, localtimodel, turbine_id=(6+ 4*10), tol=1E-6)
        @test ti_dst  ≈ data[6,2] atol=atol

    end

    @testset "Local TI Model No Local TI" begin

        # load model set
        include("./model_sets/model_set_2.jl")

        # calculate turbine inflow velocities
        ff.turbine_velocities_one_direction!(rotor_points_y, rotor_points_z, ms2, pd2)

        # load horns rev ti ata
        data = readdlm("inputfiles/horns_rev_ti_by_row_niayifar.txt", ',', skipstart=1)

        # freestream
        ti_dst = ff.calculate_local_ti(ambient_ti, windfarm, windfarmstate, localtimodel, turbine_id=(1+ 4*10), tol=1E-6)
        @test ti_dst == ambient_ti 

        # row 2
        ti_dst = ff.calculate_local_ti(ambient_ti, windfarm, windfarmstate, localtimodel, turbine_id=(2+ 4*10), tol=1E-6)
        @test ti_dst == ambient_ti 

        # row 3
        ti_dst = ff.calculate_local_ti(ambient_ti, windfarm, windfarmstate, localtimodel, turbine_id=(3+ 4*10), tol=1E-6)
        @test ti_dst == ambient_ti 

        # row 4
        ti_dst = ff.calculate_local_ti(ambient_ti, windfarm, windfarmstate, localtimodel, turbine_id=(4+ 4*10), tol=1E-6)
        @test ti_dst == ambient_ti 

        # row 5
        ti_dst = ff.calculate_local_ti(ambient_ti, windfarm, windfarmstate, localtimodel, turbine_id=(5+ 4*10), tol=1E-6)
        @test ti_dst == ambient_ti 

        # row 6
        ti_dst = ff.calculate_local_ti(ambient_ti, windfarm, windfarmstate, localtimodel, turbine_id=(6+ 4*10), tol=1E-6)
        @test ti_dst == ambient_ti 

    end

end


@testset "General Models" begin

    @testset "Coordinate rotation" begin
        atol = 1E-15

        xlocs = [-1.0 1.0]
        ylocs = [0.0 0.0]

        wind_direction_met = 0.0
        xnew, ynew = ff.rotate_to_wind_direction(xlocs, ylocs, wind_direction_met)
        @test xnew ≈ [0.0 0.0] atol=atol
        @test ynew ≈ [-1.0 1.0] atol=atol

        wind_direction_met = 3*pi/2
        xnew, ynew = ff.rotate_to_wind_direction(xlocs, ylocs, wind_direction_met)
        @test xnew ≈ [-1.0 1.0] atol=atol
        @test ynew ≈ [0.0 0.0] atol=atol
    end

    @testset "Point velocity" begin

        rtol = 1E-6

        include("./model_sets/model_set_1.jl")

        # test no loss upstream (data from Jensen 1983)
        expected_velocity = wind_speed
        loc = [-0.1, 0.0, hub_height]
        @test ff.point_velocity(loc, ms1, pd1, wind_farm_state_id=1, downwind_turbine_id=0) == expected_velocity

        # test max loss at turbine (data from Jensen 1983)
        expected_velocity = wind_speed*(1.0 - (2. * 1/3.0))
        loc = [1E-6, 0.0, hub_height]
        @test ff.point_velocity(loc, ms1, pd1, wind_farm_state_id=1, downwind_turbine_id=0) ≈ expected_velocity rtol=rtol

        # test centerline loss 40 meters downstream (data from Jensen 1983)
        expected_velocity = wind_speed*(4.35/8.1)
        loc = [40.0, 0.0, hub_height]
        @test ff.point_velocity(loc, ms1, pd1, wind_farm_state_id=1, downwind_turbine_id=0) ≈ expected_velocity rtol=rtol

        # test centerline loss 100 meters downstream (data from Jensen 1983)
        expected_velocity = wind_speed*(5.7/8.1)
        loc = [100.0, 0.0, hub_height]
        @test ff.point_velocity(loc, ms1, pd1, wind_farm_state_id=1, downwind_turbine_id=0) ≈ expected_velocity  rtol=rtol
    end

    # @testset "Turbine Inflow Velocities one direction" begin
    #     # test based on:
    #     # [1] An Aero-acoustic Noise Distribution Prediction Methodology for Offshore Wind Farms
    #     # by Jiufa Cao, Weijun Zhu, Xinbo Wu, Tongguang Wang, and Haoran Xu
    #
    #     rtol = 1E-6
    #
    #     data = readdlm("inputfiles/velocity_def_row_of_10_turbs.txt",  ',', skipstart=4)
    #
    #     rtol = 1E-6
    #
    #     include("./model_sets/model_set_2.jl")
    #     # test no loss upstream (data from Jensen 1983)
    #     expected_velocity = wind_speed*data[:, 2]
    #
    #     ff.turbine_velocities_one_direction!(rotor_points_y, rotor_points_z, ms2, pd2)
    #
    #     @test windfarmstate.turbine_inflow_velcities ≈ expected_velocity rtol=rtol
    #
    # end

    # @testset "Turbine powers one direction" begin
    #     # test based on:
    #     # [1] An Aero-acoustic Noise Distribution Prediction Methodology for Offshore Wind Farms
    #     # by Jiufa Cao, Weijun Zhu, Xinbo Wu, Tongguang Wang, and Haoran Xu
    #
    #     rtol = 1E-6
    #
    #     data = readdlm("inputfiles/velocity_def_row_of_10_turbs.txt",  ',', skipstart=4)
    #
    #     rtol = 1E-6
    #
    #     include("./model_sets/model_set_2.jl")
    #
    #     ff.turbine_powers_one_direction!(rotor_points_y, rotor_points_z, pd2)
    #
    #     @test windfarmstate.turbine_generators_powers ≈ wind_speed*data[:, 2] rtol=rtol
    #
    # end

    @testset "Calculate flow field" begin
        # test based on:
        # [1] An Aero-acoustic Noise Distribution Prediction Methodology for Offshore Wind Farms
        # by Jiufa Cao, Weijun Zhu, Xinbo Wu, Tongguang Wang, and Haoran Xu

        rtol = 1E-6


        data = readdlm("inputfiles/velocity_def_row_of_10_turbs.txt",  ',', skipstart=4)

        include("./model_sets/model_set_2.jl")

        ff.turbine_velocities_one_direction!(rotor_points_y, rotor_points_z, ms2, pd2)
        ff.turbine_powers_one_direction!(rotor_points_y, rotor_points_z, pd2)

        stepsize = 10
        xrange = 1:stepsize:10*rotor_diameter*nturbines
        yrange = -1*rotor_diameter*nturbines:stepsize:1*rotor_diameter*nturbines
        zrange = hub_height:stepsize:hub_height
        points = [[x y z] for x in xrange for y in yrange for z in zrange]
        println("size: ", size(points[:]))
        println(points[1][1])
        println(pd2.wind_farm_states[1].turbine_inflow_velcities)

        velh = ff.calculate_flow_field(1, xrange, yrange, zrange, rotor_points_y, rotor_points_z, ms2, pd2)
        # yrange = 0:stepsize:0
        # zrange = 0:stepsize:2*hub_height
        # velv = ff.calculate_flow_field(1, xrange, yrange, zrange, rotor_points_y, rotor_points_z, windfarm, windfarmstate, windresource, windshearmodel, wakedeficitmodel, wakedeflectionmodel, wakecombinationmodel)

        # println(vel)

        # Plots.pyplot()
        # points = permutedims(reshape(hcat(points...), (length(points[1]), length(points))))
        # println(size(points[:][1]))
        # println(size(points[:][2]))
        # println(size(vel))
        # println(size(vel), length(xrange), length(yrange))
        # cmap = Colormap("Blues_r", 100)
        # flowfieldplot = contourf(xrange, yrange, vel, cmap="Blues_r")
        flowfieldplot = imshow(velh, cmap="Blues_r")
        # flowfieldplot = imshow(velv, cmap="Blues_r")
        # x = range(-5, stop = 5, length = 40)
        # flowfieldplot = plot(points[:, 1], points[:, 2], vel, show=true, st=:contourf)

    end

end
