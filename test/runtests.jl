import FlowFarm; const ff = FlowFarm
using Test
using DelimitedFiles
using LinearAlgebra
using PyPlot 

@testset "Thrust Coefficient Models" begin

    @testset "Constant Thrust Coefficent Model" begin

        ct = 0.7
        ct_model = ff.ConstantCt(ct)

        @test ff.calculate_ct(ct_model) == ct

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

        windfarmstate = ff.SingleWindFarmState(wind_farm_state_id, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai, sorted_turbine_index, turbine_inflow_velcity, [0.0])       

        k_star = 0.07 # adjusted to match experimental data. #TODO improve tests with model results
        horizontal_spread_rate = k_star

        ct_model = ff.ConstantCt(ct)
        power_model = ff.ConstantCp([cp], [generator_efficiency])

        turbine = ff.Turbine([turbine_id], [rotor_diameter], [hub_height], [ct_model], [power_model])
        model = ff.JiminezYawDeflection(horizontal_spread_rate)

        dx2p5d = 2.5*rotor_diameter
        dy2p5d_y20 = 0.23136246786632375*rotor_diameter # from [1] figure 21

        dx5p5d = 5.5*rotor_diameter
        dy5p5d_y20 = 0.40616966580976843*rotor_diameter # from [1] figure 21

        dx8d = 8.0*rotor_diameter
        dy8d_20 = 0.5257731958762868*rotor_diameter # from [1] figure 21

        # test deflection at 2.5D with yaw 20 deg
        @test round(ff.wake_deflection_model([dx2p5d, 0.0, hub_height], model, turbine, windfarmstate), digits=2) == round(dy2p5d_y20, digits=2)

        # test deflection at 5.5D with yaw 20 deg
        @test round(ff.wake_deflection_model([dx5p5d, 0.0, hub_height], model, turbine, windfarmstate), digits=2) == round(dy5p5d_y20, digits=2)

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

        windfarmstate = ff.SingleWindFarmState(wind_farm_state_id, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai, sorted_turbine_index, turbine_inflow_velcity, [0.0])       
        ct_model = ff.ConstantCt(ct)
        power_model = ff.ConstantCp([cp], [generator_efficiency])

        turbine = ff.Turbine([turbine_id], [rotor_diameter], [hub_height], [ct_model], [power_model])

        k_star = 0.022 # [1]  p. 525
        turbulence_intensity = 0.1 #0.0875 #[1] p. 508 - this value is only specified to be less than 0.1
        horizontal_spread_rate = k_star
        vertical_spread_rate = k_star
        alpha_star = 2.32 #[1] p. 534
        beta_star = 0.154 #[1] p. 534

        model = ff.GaussYawDeflection(turbulence_intensity, horizontal_spread_rate, vertical_spread_rate, alpha_star, beta_star)

        dx4d = 4.0*rotor_diameter
        dy4d_20 = 0.2684659090909065*rotor_diameter # from [1] figure 21

        dx8d = 8.0*rotor_diameter
        dy8d_20 = 0.34090909090908905*rotor_diameter # from [1] figure 21
        
        # test deflection at 4D with yaw 20 deg
        @test ff.wake_deflection_model([dx4d, dy4d_20, hub_height], model, turbine, windfarmstate) ≈ dy4d_20 atol=atol

        # test deflection at 8D with yaw 20 deg
        @test ff.wake_deflection_model([dx8d, dy8d_20, hub_height], model, turbine, windfarmstate) ≈ dy8d_20 atol=atol

    end

end

@testset "Wake Deficit Models" begin
    
    @testset "Jensen Top Hat Model" begin

        rotor_diameter = 40.0
        hub_height = 90.0
        yaw = 0.0
        ct = 0.7 
        ai = 1.0/3.0
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
        turbine_inflow_velcity = [8.1]
        turbine_id = 1

        windfarmstate = ff.SingleWindFarmState(wind_farm_state_id, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai, sorted_turbine_index, turbine_inflow_velcity, [0.0])       
        ct_model = ff.ConstantCt(ct)
        power_model = ff.ConstantCp([cp], [generator_efficiency])

        turbine = ff.Turbine([turbine_id], [rotor_diameter], [hub_height], [ct_model], [power_model])

        alpha = 0.1
        
        deflection = [0.0, 0.0]

        model = ff.JensenTopHat(alpha)

        centerloss40 = 1. - 4.35/8.1
        centerloss100 = 1. - 5.7/8.1

        # test no loss upstream (data from Jensen 1983)
        @test ff.wake_deficit_model([-1E-12, 0.0, hub_height], deflection, model, turbine, windfarmstate) == 0.0

        # test max loss at turbine (data from Jensen 1983)
        @test ff.wake_deficit_model([0.0, 0.0, hub_height], deflection, model, turbine, windfarmstate) == (2. * 1/3.0)

        # test centerline loss 40 meters downstream (data from Jensen 1983)
        @test ff.wake_deficit_model([40., 0.0, hub_height], deflection, model, turbine, windfarmstate) == centerloss40

        # test centerline loss 100 meters downstream (data from Jensen 1983)
        @test ff.wake_deficit_model([100., 0.0, hub_height], deflection, model, turbine, windfarmstate) == centerloss100

        # test wake diameter 40 meters downstream (data from Jensen 1983)
        @test ff.wake_deficit_model([40., (alpha*40 + rotor_diameter/2.), hub_height], deflection, model, turbine, windfarmstate) == centerloss40
        @test ff.wake_deficit_model([40., (alpha*40 + rotor_diameter/2. + 1E-12), hub_height], deflection, model, turbine, windfarmstate) == 0.0
        @test ff.wake_deficit_model([40., -(alpha*40 + rotor_diameter/2.), hub_height], deflection, model, turbine, windfarmstate) == centerloss40
        @test ff.wake_deficit_model([40., -(alpha*40 + rotor_diameter/2. + 1E-12), hub_height], deflection, model, turbine, windfarmstate) == 0.0

        # test wake diameter 100 meters downstream (data from Jensen 1983)
        @test ff.wake_deficit_model([100., (alpha*100. + rotor_diameter/2.), hub_height], deflection, model, turbine, windfarmstate) == centerloss100
        @test ff.wake_deficit_model([100., (alpha*100. + rotor_diameter/2. + 1E-12), hub_height], deflection, model, turbine, windfarmstate) == 0.0
        @test ff.wake_deficit_model([100., -(alpha*100. + rotor_diameter/2.), hub_height], deflection, model, turbine, windfarmstate) == centerloss100
        @test ff.wake_deficit_model([100., -(alpha*100. + rotor_diameter/2. + 1E-12), hub_height], deflection, model, turbine, windfarmstate) == 0.0
    end

    @testset "Jensen Cosine Model" begin

        rotor_diameter = 40.0
        hub_height = 90.0
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
        turbine_id = 1

        alpha = 0.1
        beta = 20.0*pi/180.0


        deflection = [0.0, 0.0]
        
        windfarmstate = ff.SingleWindFarmState(wind_farm_state_id, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai, sorted_turbine_index, turbine_inflow_velcity, [0.0])       
        ct_model = ff.ConstantCt(ct)
        power_model = ff.ConstantCp([cp], [generator_efficiency])

        turbine = ff.Turbine([turbine_id], [rotor_diameter], [hub_height], [ct_model], [power_model])

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
        @test ff.wake_deficit_model([-1E-12, 0.0, hub_height], deflection, model, turbine, windfarmstate) == 0.0

        # test max loss at turbine (data from Jensen 1983)
        @test ff.wake_deficit_model([0.0, 0.0, hub_height], deflection, model, turbine, windfarmstate) == (2. * aI)

        # test centerline loss 40 meters downstream (data from Jensen 1983)
        @test ff.wake_deficit_model([40., 0.0, hub_height], deflection, model, turbine, windfarmstate) == centerloss40

        # test centerline loss 100 meters downstream (data from Jensen 1983)
        @test ff.wake_deficit_model([100., 0.0, hub_height], deflection, model, turbine, windfarmstate) == centerloss100

        # test wake diameter 40 meters downstream (data from Jensen 1983)
        @test ff.wake_deficit_model([40., dy40, hub_height], deflection, model, turbine, windfarmstate) == 0.0
        @test ff.wake_deficit_model([40., (dy40 + 1E-12), hub_height], deflection, model, turbine, windfarmstate) == 0.0
        @test ff.wake_deficit_model([40., (dy40 - 1E1), hub_height], deflection, model, turbine, windfarmstate) >= 0.0
        @test ff.wake_deficit_model([40., -(dy40), hub_height], deflection, model, turbine, windfarmstate) == 0.0
        @test ff.wake_deficit_model([40., -(dy40 + 1E-12), hub_height], deflection, model, turbine, windfarmstate) == 0.0
        @test ff.wake_deficit_model([40., -(dy40 - 1E1), hub_height], deflection, model, turbine, windfarmstate) >= 0.0

        # test wake diameter 100 meters downstream (data from Jensen 1983)
        @test ff.wake_deficit_model([100., dy100, hub_height], deflection, model, turbine, windfarmstate) == 0.0
        @test ff.wake_deficit_model([100., (dy100 + 1E-12), hub_height], deflection, model, turbine, windfarmstate) == 0.0
        @test ff.wake_deficit_model([100., (dy100 - 1E1), hub_height], deflection, model, turbine, windfarmstate) >= 0.0
        @test ff.wake_deficit_model([100., -(dy100), hub_height], deflection, model, turbine, windfarmstate) == 0.0
        @test ff.wake_deficit_model([100., -(dy100 + 1E-12), hub_height], deflection, model, turbine, windfarmstate) == 0.0
        @test ff.wake_deficit_model([100., -(dy100 - 1E1), hub_height], deflection, model, turbine, windfarmstate) >= 0.0

        # test value at point in wake 40 m downstream and with theta=15 degrees
        @test ff.wake_deficit_model([40., dy, hub_height], deflection, model, turbine, windfarmstate) == loss40attheta

    end
    
    @testset "Gauss Yaw Model" begin
        rtol = 0.1
        # [1] based on data from Bastankhah and Porte-Agel 2016

        turbine_x = [0.0]
        turbine_y = [0.0]
        turbine_z = [0.000022] #[1] p. 509
        rotor_diameter = 0.15 #[1] p. 509
        hub_height = 0.125 #[1] p. 509
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
        turbine_id = 1

        deflection = [0.0, 0.0]
        
        windfarmstate = ff.SingleWindFarmState(wind_farm_state_id, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai, sorted_turbine_index, turbine_inflow_velcity, [0.0])       
        ct_model = ff.ConstantCt(ct)
        power_model = ff.ConstantCp([cp], [generator_efficiency])

        turbine = ff.Turbine([turbine_id], [rotor_diameter], [hub_height], [ct_model], [power_model])

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
        @test ff.wake_deficit_model([-1E-12, 0.0, hub_height], deflection, model, turbine, windfarmstate) == 0.0

        # test loss at x1 with no yaw
        @test ff.wake_deficit_model([x1_0, 0.0, hub_height], deflection, model, turbine, windfarmstate) ≈ loss1_0 rtol=rtol 

        deflection = [0.0, 0.0]
        # test loss at x2 with no yaw
        @test ff.wake_deficit_model([x2_0, 0.0, hub_height], deflection, model, turbine, windfarmstate)  ≈ loss2_0 rtol=rtol 

        # test loss at x3 with no yaw
        @test ff.wake_deficit_model([x3_0, 0.0, hub_height], deflection, model, turbine, windfarmstate) ≈ loss3_0 rtol=rtol 

        # test loss at x4 with no yaw
        @test ff.wake_deficit_model([x4_0, 0.0, hub_height], deflection, model, turbine, windfarmstate)  ≈ loss4_0 rtol=rtol 
        
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

    @testset "Turbine Inflow Velocities one direction" begin
        # test based on:
        # [1] An Aero-acoustic Noise Distribution Prediction Methodology for Offshore Wind Farms
        # by Jiufa Cao, Weijun Zhu, Xinbo Wu, Tongguang Wang, and Haoran Xu

        rtol = 1E-6

        data = readdlm("inputfiles/velocity_def_row_of_10_turbs.txt",  ',', skipstart=4)

        rtol = 1E-6

        include("./model_sets/model_set_2.jl")

        # test no loss upstream (data from Jensen 1983)
        expected_velocity = wind_speed*data[:, 2]

        ff.turbine_velocities_one_direction!(rotor_points_y, rotor_points_z, ms2, pd2)
        
        @test windfarmstate.turbine_inflow_velcities ≈ expected_velocity rtol=rtol 

    end

    @testset "Turbine powers one direction" begin
        # test based on:
        # [1] An Aero-acoustic Noise Distribution Prediction Methodology for Offshore Wind Farms
        # by Jiufa Cao, Weijun Zhu, Xinbo Wu, Tongguang Wang, and Haoran Xu

        rtol = 1E-6

        data = readdlm("inputfiles/velocity_def_row_of_10_turbs.txt",  ',', skipstart=4)

        rtol = 1E-6

        # rotor_diameter = 80.0
        # hub_height = 70.0
        # yaw = 0.0
        # ct = 0.689 
        # cp = 0.8
        # generator_efficiency = 0.944
        # ai = 1.0/3.0
        # wind_speed = 12.0
        # air_density = 1.1716  # kg/m^3
        # turbine_x = data[:, 1].*7.0*rotor_diameter #[0.0, 7.0*rotor_diameter, 6.0*rotor_diameter]
        # nturbines = length(turbine_x)
        # turbine_y = zeros(nturbines)
        # turbine_z = zeros(nturbines)
        # turbine_yaw = zeros(nturbines)
        # turbine_ct = zeros(nturbines) .+ ct
        # turbine_ai = zeros(nturbines) .+ ai
        # winddirections = [270.0*pi/180.0]
        # windspeeds = [wind_speed]
        # windprobabilities = [1.0]
        # measurementheight = [hub_height]
        # shearexponent = [0.15]
        # turbine_inflow_velcities = wind_speed*data[:, 2]
        # # rotor sample points 
        # rotor_points_y = [0.0]
        # rotor_points_z = [0.0]

        # ct_model = ff.ConstantCt(ct)
        # power_model = ff.ConstantCp([cp], [generator_efficiency])
        # wind_shear_model = ff.PowerLawWindShear(shearexponent)

        # # turbine_definitions = Array{ff.AbstractTurbine, length(turbine_x)}
        # turbine_definitions = [] # or `Vector{Coords}(undef, x)`
        # turbine = ff.Turbine(1, [rotor_diameter], [hub_height], [ct_model], [power_model])
        # turbine_definitions = [turbine for i in 1:nturbines]
        # sorted_turbine_index = zeros(Int64, nturbines)
        # for i in 1:nturbines
        #     turbine = ff.Turbine(i, [rotor_diameter], [hub_height], [ct_model], [power_model])
        #     turbine_definitions[i] = turbine
        #     sorted_turbine_index[i] = Int64(i)
        # end 

        # turbine_generators_powers = zeros(nturbines)

        # windfarm = ff.WindFarm(turbine_x, turbine_y, turbine_z, turbine_definitions)
        # windfarmstate = ff.SingleWindFarmState(1, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai, sorted_turbine_index, turbine_inflow_velcities, turbine_generators_powers)
        # windresource = ff.DiscretizedWindResource(winddirections, windspeeds, windprobabilities, measurementheight, air_density, [wind_shear_model])

        # alpha = 0.1
        # wakedeficitmodel = ff.JensenTopHat(alpha)
        # horizontal_spread_rate = alpha
        # wakedeflectionmodel = ff.JiminezYawDeflection(horizontal_spread_rate)
        # wakecombinationmodel = ff.SumOfSquaresFreestreamSuperposition()

        include("./model_sets/model_set_2.jl")

        ff.turbine_powers_one_direction!(rotor_points_y, rotor_points_z, pd2)
        
        @test windfarmstate.turbine_generators_powers ≈ wind_speed*data[:, 2] rtol=rtol 

    end

    @testset "Calculate flow field" begin
        # test based on:
        # [1] An Aero-acoustic Noise Distribution Prediction Methodology for Offshore Wind Farms
        # by Jiufa Cao, Weijun Zhu, Xinbo Wu, Tongguang Wang, and Haoran Xu

        rtol = 1E-6


        data = readdlm("inputfiles/velocity_def_row_of_10_turbs.txt",  ',', skipstart=4)

        # rotor_diameter = 80.0
        # hub_height = 70.0
        # yaw = 0.0
        # ct = 0.689 
        # cp = 0.8
        # generator_efficiency = 0.944
        # ai = 1.0/3.0
        # wind_speed = 12.0
        # air_density = 1.1716  # kg/m^3
        # turbine_x = data[:, 1].*7.0*rotor_diameter #[0.0, 7.0*rotor_diameter, 6.0*rotor_diameter]
        # nturbines = length(turbine_x)
        # turbine_y = zeros(nturbines)
        # turbine_z = zeros(nturbines)
        # turbine_yaw = zeros(nturbines)
        # turbine_ct = zeros(nturbines) .+ ct
        # turbine_ai = zeros(nturbines) .+ ai
        # winddirections = [3.0*pi/2]
        # windspeeds = [wind_speed]
        # windprobabilities = [1.0]
        # measurementheight = [hub_height]
        # shearexponent = 0.15
        # turbine_inflow_velcities = wind_speed*data[:, 2]
        # # rotor sample points 
        # rotor_points_y = [0.0]
        # rotor_points_z = [0.0]

        # ct_model = ff.ConstantCt(ct)
        # power_model = ff.ConstantCp([cp], [generator_efficiency])
        # wind_shear_model = [ff.PowerLawWindShear(shearexponent)]

        # # turbine_definitions = Array{ff.AbstractTurbine, length(turbine_x)}
        # turbine_definitions = [] # or `Vector{Coords}(undef, x)`
        # turbine = ff.Turbine(1, [rotor_diameter], [hub_height], [ct_model], [power_model])
        # turbine_definitions = [turbine for i in 1:nturbines]
        # sorted_turbine_index = zeros(Int64, nturbines)
        # for i in 1:nturbines
        #     turbine = ff.Turbine(i, [rotor_diameter], [hub_height], [ct_model], [power_model])
        #     turbine_definitions[i] = turbine
        #     sorted_turbine_index[i] = Int64(i)
        # end 

        # turbine_generators_powers = zeros(nturbines)

        # windfarm = ff.WindFarm(turbine_x, turbine_y, turbine_z, turbine_definitions)
        # turbine_xw, turbine_yw = ff.rotate_to_wind_direction(turbine_x, turbine_y, winddirections[1])
        # windfarmstate = ff.SingleWindFarmState(1, turbine_xw, turbine_yw, turbine_z, turbine_yaw, turbine_ct, turbine_ai, sorted_turbine_index, turbine_inflow_velcities, turbine_generators_powers)
        # windresource = ff.DiscretizedWindResource(winddirections, windspeeds, windprobabilities, measurementheight, air_density, wind_shear_model)

        # alpha = 0.1
        # wakedeficitmodel = ff.JensenTopHat(alpha)
        # horizontal_spread_rate = alpha
        # wakedeflectionmodel = ff.JiminezYawDeflection(horizontal_spread_rate)
        # wakecombinationmodel = ff.SumOfSquaresFreestreamSuperposition()
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
