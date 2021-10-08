using FLOWFarm; const ff = FLOWFarm
using Test
using TestSetExtensions
using Distributed
using DelimitedFiles
using LinearAlgebra
using FLOWMath: linear
using Distributed
using YAML
using ForwardDiff
using FiniteDiff

@testset ExtendedTestSet "all tests" begin
    @testset "cost_models" begin
        Parameters = ff.Levelized()
        rotor_diameter = [70]
        hub_height = [65]
        AEP = 3734
        rated_power = [5000]
        COE = ff.cost_of_energy(rotor_diameter, hub_height, rated_power, AEP, Parameters)
        @test COE ≈ 33.69894269249007 atol=1E-6
    end

    @testset "utilities" begin

        @testset "latitude longitude to xy" begin

            latitude = [59.3293, 59.8586]
            longitude = [18.0686, 17.6389] 
            utm_zone = 33
            x, y = ff.latlong_to_xy(latitude, longitude, utm_zone)

            @test x ≈ [0.0, -26778.38032168697] atol=1E-6
            @test y ≈ [0.0, 57865.048037721775] atol=1E-6
        end

        @testset "hermite spline" begin

            x = 0.5
            x1 = 0.0
            x2 = 1.0
            y1 = 0.0
            y2 = 1.0
            dydx1 = 1.0
            dydx2 = 1.0

            v = ff.hermite_spline(x, x1, x2, y1, dydx1, y2, dydx2)
            @test v == 0.5

            v, dv = ff.hermite_spline(x, x1, x2, y1, dydx1, y2, dydx2; return_deriv=true)
            @test v == 0.5
            @test dv == 1.0

            dydx1 = 0.0
            dydx2 = 0.0
            v, dv = ff.hermite_spline(x, x1, x2, y1, dydx1, y2, dydx2; return_deriv=true)
            @test v == 0.5
            @test dv == 1.5

        end

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

            c = [-30.0;-2.0;0.0;1.0;y]
            m = ff.smooth_max(c)
            @test m ≈ y atol=1E-4

            c = [-30.0;-2.0;1.98;1.99;y]
            m = ff.smooth_max(c, s=400)
            @test m ≈ y atol=1E-4

            c = [-30.0;-20.0;-10.0;-4.0;y]
            m = ff.smooth_max(c, s=4)
            @test m ≈ y atol=1E-6

        end

        @testset "boundary normals" begin

            boundary_file_name = string("./inputfiles/iea37-boundary-cs3.yaml")
            
            boundary_vertices = ff.get_boundary_yaml(boundary_file_name)

            boundary_normals = ff.boundary_normals_calculator(boundary_vertices)
        
            correct_normals = [0.9829601758936983 -0.1838186405319916; 
                                0.9934614633172962 -0.11416795042154541; 
                                0.9987121579438882 -0.050734855622757584; 
                                0.9998686751666075 -0.01620593781838486; 
                                0.9999954987444023 0.0030004151269687495; 
                                -0.9998078216567232 -0.019604074934516894; 
                                -0.6957179389375846 -0.718315076718037; 
                                -0.6957275377423737 -0.7183057797532565; 
                                -0.8019887481131871 0.5973391397688945; 
                                0.5138086803485797 0.8579047965820281; 
                                0.4252760929807897 0.905063668886888; 
                                0.2645057513093967 0.9643841078762402; 
                                -0.0684295708121141 0.9976559496331737; 
                                -0.39636379138742883 0.9180935381958544; 
                                -0.6828023205475376 0.7306031693435896; 
                                -0.7996740386176392 0.6004343694034798; 
                                -0.8578802011411015 0.5138497450520954; 
                                0.42552559023380465 0.9049463918134445]

            @test boundary_normals ≈ correct_normals atol=1E-6

            boundary_vertices = [0 0; 0 1; 1 1; 1 0]
            boundary_normals = ff.boundary_normals_calculator(boundary_vertices)
            correct_normals = [-1 0; 0 1; 1 0; 0 -1]
            @test boundary_normals ≈ correct_normals atol=1E-6

        end

        @testset "sunflower_points" begin

            x, y = ff.sunflower_points(10)
            xtest = [-0.16916402229765054, 0.03473946036824235, 0.3121225499658305, -0.597698416169296, 0.5807122204996749, -0.19752925784048528, -0.3812485520657022, 0.8346088735992109, -0.8743433632876815, 0.4238459950479107]
            ytest = [0.15496810158044924, -0.3958382330389885, 0.40710859551189754, -0.10572443397953969, -0.36940158024655595, 0.7347989934111504, -0.7340708874922061, 0.30479782203943484, 0.36091623014218815, -0.9057342725556136]
            @test x ≈ xtest atol=1E-6
            @test y ≈ ytest atol=1E-6
            
        end

        @testset "grid_points" begin 

            # test with grid at size of rotor-swept area including perimeter
            y, z = ff.grid_points(9)
            ytest = [-1, -1, -1, 0, 0, 0, 1, 1, 1]
            ztest = [-1, 0, 1, -1, 0, 1, -1, 0, 1]
            @test y ≈ ytest atol=1E-6
            @test z ≈ ztest atol=1E-6

        end

        @testset "rotor_sample_points" begin

            # sunflower points method
            x, y = ff.rotor_sample_points(10)
            xtest = [-0.16916402229765054, 0.03473946036824235, 0.3121225499658305, -0.597698416169296, 0.5807122204996749, -0.19752925784048528, -0.3812485520657022, 0.8346088735992109, -0.8743433632876815, 0.4238459950479107]
            ytest = [0.15496810158044924, -0.3958382330389885, 0.40710859551189754, -0.10572443397953969, -0.36940158024655595, 0.7347989934111504, -0.7340708874922061, 0.30479782203943484, 0.36091623014218815, -0.9057342725556136]
            @test x ≈ xtest atol=1E-6
            @test y ≈ ytest atol=1E-6

            # sunflower points method
            x, y = ff.rotor_sample_points(1)
            xtest = [0.0]
            ytest = [0.0]
            @test x ≈ xtest atol=1E-6
            @test y ≈ ytest atol=1E-6

            # test with grid wholly inside rotor-swept area
            x, y = ff.rotor_sample_points(9, method="grid", pradius=0.5, use_perimeter_points=false)
            xtest = [-0.5, -0.5, -0.5, 0.0, 0.0, 0.0, 0.5, 0.5, 0.5]
            ytest = [-0.5, 0.0, 0.5, -0.5, 0.0, 0.5, -0.5, 0.0, 0.5]
            @test x ≈ xtest atol=1E-6
            @test y ≈ ytest atol=1E-6

            # test with grid at size of rotor-swept area excluding perimeter
            x, y = ff.rotor_sample_points(9, method="grid", pradius=1, use_perimeter_points=false)
            xtest = [0]
            ytest = [0]
            @test x ≈ xtest atol=1E-6
            @test y ≈ ytest atol=1E-6

            # test with grid at size of rotor-swept area including perimeter
            x, y = ff.rotor_sample_points(9, method="grid", pradius=1, use_perimeter_points=true)
            xtest = [-1, 0, 0, 0, 1]
            ytest = [0, -1, 0, 1, 0]
            @test x ≈ xtest atol=1E-6
            @test y ≈ ytest atol=1E-6
            
        end

        @testset "wake_count_iec" begin 
            turbinex = [0, 100, 200, 300]
            turbiney = zeros(4)
            diameter = ones(4).*40

            # test with turbines in a row
            winddirection = 3.0*pi/2.0
            correct_wake_count = [0, 1, 2, 3]
            wake_count = ff.wake_count_iec(turbinex, turbiney, winddirection, diameter)
            @test wake_count == correct_wake_count

            # test with turbines all in free-stream 
            winddirection = 0.0
            correct_wake_count = zeros(4)
            wake_count = ff.wake_count_iec(turbinex, turbiney, winddirection, diameter)
            @test wake_count == correct_wake_count

            # test with multiple directions
            winddirection = [3.0*pi/2.0, 0.0]
            correct_wake_count = [0 1 2 3; 0 0 0 0]
            wake_count = ff.wake_count_iec(turbinex, turbiney, winddirection, diameter)
            @test wake_count == correct_wake_count
        end

        @testset "find_upstream_turbines" begin 
            turbinex = [0.0 100.0]
            turbiney = [0.0 0.0]
            diameter = [20.0 20.0]
            winddirection = [0.0, 3.0*pi/2.0]

            # test for no waked turbines
            upstream_turbines = ff.find_upstream_turbines(turbinex, turbiney, winddirection[1], diameter, inverse=false)
            @test upstream_turbines == [1, 2]

            # test with one waked turbine
            upstream_turbines = ff.find_upstream_turbines(turbinex, turbiney, winddirection[2], diameter, inverse=false)
            @test upstream_turbines == [1]

            # test with two directions
            upstream_turbines = ff.find_upstream_turbines(turbinex, turbiney, winddirection, diameter, inverse=false)
            @test upstream_turbines == [[1, 2],[1]]
        end

        @testset "nansafesqrt" begin 
            
            a1 = 9.0
            a2 = 0.0
            a3 = eps()*1E-10

            # test for sqrt region
            @test ff.nansafesqrt(a1) == 3.0

            # test zero
            @test ff.nansafesqrt(a2) == 0.0

            # test linearly approximated region a*sqrt(tol)/tol
            @test isapprox(ff.nansafesqrt(a3), a3*sqrt(eps())/eps())

        end


    end

    @testset "Optimization functions" begin

        @testset "Turbine spacing calculation" begin

            testing_x = [0.0,0.0,500.0,500.0,7481.0]
            testing_y = [0.0,500.0,0.0,500.0,-43891.0]
            turbine_separation = [500.0,500.0,707.1067811865476,44523.9850193129,707.1067811865476,500.0,45016.95505029189,500.0,44442.70741077775,44936.56909466943]

            @test ff.turbine_spacing(testing_x,testing_y) == turbine_separation

        end

        @testset "Circular boundary" begin

            center = [100,500]
            radius = 500
            testing_x = [100,100,100,100,-500,700]
            testing_y = [500,1100,-100,0,500,500]
            test_values = (testing_x .- center[1]).^2 + (testing_y .- center[2]).^2 .- radius.^2

            @test ff.circle_boundary(center,radius,testing_x,testing_y) == test_values

        end

        @testset "Splined boundary" begin

            @testset "One Turbine, Circular Boundary" begin
                #-- One-turbine circular boundary as a square --#
                # A discretized 20-point circle
                bndry_x_clsd = [200.00, 195.11, 180.90, 158.78, 130.90, 100.00, 69.10, 41.22, 19.10, 4.89, 0.00, 4.89, 19.10, 41.22, 69.10, 100.00, 130.90, 158.78, 180.90, 195.11, 200.00]
                bndry_y_clsd = [100.00, 130.90, 158.78, 180.90, 195.11, 200.00, 195.11, 180.90, 158.78, 130.90, 100.00, 69.10, 41.22, 19.10, 4.89, 0.00, 4.89, 19.10, 41.22, 69.10, 100.00]
                # Vertices that keep splines injective (4-corners)
                bndry_corner_indcies =[1,6,11,16, 21]  # 20 pt circle, 4 corners
                # Should be equidistant from sides
                testing_x = [100.0]
                testing_y = [100.0]
                test_values = [-100.0, -100.0, -100.0, -100.0]

                @test ff.splined_boundary(testing_x, testing_y, bndry_x_clsd, bndry_y_clsd, bndry_corner_indcies) == test_values

                #-- One-turbine circular boundary as a triangle --#
                # Vertices that keep splines injective (3-corners)
                bndry_corner_indcies =[1,11,17,21]  # 20 pt circle, 3 corners
                @test ff.splined_boundary(testing_x, testing_y, bndry_x_clsd, bndry_y_clsd, bndry_corner_indcies) == test_values
            end

            @testset "Multi-Turbine, Circular Boundary" begin
                #-- Multi-turbine circular boundary as a square --#
                # A discretized 200-point circle
                num_pts = 200
                circ_radius = 100.0
                circ_center = [100.0, 100.0]
                bndry_x_clsd, bndry_y_clsd = ff.DiscreteCircum(circ_center[1], circ_center[2], circ_radius, num_pts)
                # Vertices that keep splines injective (4-corners)
                bndry_corner_indcies = [1, 51, 101, 151, 201]  # 200 pt circle, 4 corners

                # Vertices that keep splines injective
                circ_corners = 1/sqrt(2)
                cc_r = circ_center[2] + (circ_corners * circ_radius) # ~170
                cc_l = circ_center[1] - (circ_corners * circ_radius) # ~ 30
                cc_d = cc_r - cc_l                # ~140
                cc_d2 = cc_d/2                    # ~ 70

                testing_x = [cc_l, 100.0, cc_r, cc_l, 100.0, cc_r, cc_l, 100.0, cc_r]
                testing_y = [cc_r, cc_r, cc_r, 100.0, 100.0, 100.0, cc_l, cc_l, cc_l]
                
                test_values = [ -cc_r,  -cc_l,   -0.0,  -cc_d,
                                -100.0, -100.0,  -cc_l,  -cc_r,
                                -cc_l,  -cc_r,   -0.0,  -cc_d,
                                -cc_r,  -cc_l, -cc_d2, -cc_d2,
                                -100.0, -100.0, -100.0, -100.0,
                                -cc_l,  -cc_r, -cc_d2, -cc_d2,
                                -cc_r,  -cc_l,  -cc_d,   -0.0,
                                -100.0, -100.0,  -cc_r,  -cc_l,
                                -cc_l,  -cc_r,  -cc_d,   -0.0]

                ans = ff.splined_boundary(testing_x, testing_y, bndry_x_clsd, bndry_y_clsd, bndry_corner_indcies)
                # Test each turbine individually
                for i in 1:length(test_values)
                    @test ans[i] ≈ test_values[i] atol=5E-2
                end

                #-- Multi-turbine circular boundary as a triangle --#
                # Vertices that keep splines injective (3-corners)
                bndry_corner_indcies = [1, 101, 151, 201]  # 200 pt circle, 3 corners
                ans = ff.splined_boundary(testing_x, testing_y, bndry_x_clsd, bndry_y_clsd, bndry_corner_indcies)
                for i in 1:length(test_values)
                    @test ans[i] ≈ test_values[i] atol=5E-2
                end
            end

            ###### TODO: Rewrite boundary method so it doesn't fail with pt 3 being concave
            # #-- Multi-turbine circular boundary as a triangle --#
            # # Vertices that keep splines injective (3-corners)
            # bndry_corner_indcies = [1, 101, 151, 201]  # 200 pt circle, 3 corners
            # ans = ff.splined_boundary(testing_x, testing_y, bndry_x_clsd, bndry_y_clsd, bndry_corner_indcies)
            # # Test each turbine individually
            # for i in 1:length(test_values)
            #     @test ans[i] ≈ test_values[i] atol=5E-2
            # end

            # #-- Multi-turbine convex boundary test --#
            # bndry_x = [200, 100, 50, 100, 50, 100, 200]
            # bndry_y = [100, 200, 150, 100, 50, 0, 100]
            # # Vertices that keep splines injective
            # bndry_corner_indcies = [1, 3, 4, 5, 7]
            # # Make a turbine grid inside the circle
            # testing_x = [ 50, 100, 150,  50, 100, 150,  50, 100, 150]
            # testing_y = [150, 150, 150, 100, 100, 100,  50,  50,  50]

            # test_values = [150 -50 -50  50
            #                100   0 -50  50
            #                 50  50   0 100
            #                150 -50   0   0
            #                100   0   0   0
            #                 50  50  50  50
            #                150 -50  50 -50
            #                100   0  50 -50
            #                 50  50 100   0]
            # @test ff.splined_boundary(testing_x, testing_y, bndry_x_clsd, bndry_y_clsd, bndry_corner_indcies) == test_values #atol=1E-3

            #-- Multi-turbine concave boundary test --#
        
        end

        @testset "Variable Reduction (Boundary)" begin
            @testset "One Turbine, Circular Boundary, Zero Start Distance" begin
                #-- One-turbine circular boundary, zero start distance--#
                # A discretized 20-point circle
                bndry_x_clsd = [200.00, 195.11, 180.90, 158.78, 130.90, 100.00, 69.10, 41.22, 19.10, 4.89, 0.00, 4.89, 19.10, 41.22, 69.10, 100.00, 130.90, 158.78, 180.90, 195.11, 200.00]
                bndry_y_clsd = [100.00, 130.90, 158.78, 180.90, 195.11, 200.00, 195.11, 180.90, 158.78, 130.90, 100.00, 69.10, 41.22, 19.10, 4.89, 0.00, 4.89, 19.10, 41.22, 69.10, 100.00]
                # Give it one small turbine
                num_turbs = 1
                start_dist = 0
                turb_diam = 10.0
                turb_min_spacing = 2*turb_diam

                testing_x, testing_y, num_leftover = ff.VR_boundary_startup(bndry_x_clsd, bndry_y_clsd, start_dist, turb_min_spacing, num_turbs)
                test_values = [200.00, 100.0, 0]
                @test test_values[1] == bndry_x_clsd[1]  # Should line up with first coordinate
                @test test_values[2] == bndry_y_clsd[1]
                @test test_values[3] == 0                # No leftover turbines
            end

            @testset "One Turbine, Square Boundary, Zero Start Distance" begin
                # A 4-point sqaure
                bndry_x_clsd = [100.0, 0.0, 0.0, 100.0, 100.0]
                bndry_y_clsd = [100.0, 100.0, 0.0, 0.0, 100.0]
                # Give it one small turbine, at the start
                num_turbs = 1
                start_dist = 0
                turb_diam = 10.0
                turb_min_spacing = 2*turb_diam

                testing_x, testing_y, num_leftover = ff.VR_boundary_startup(bndry_x_clsd, bndry_y_clsd, start_dist, turb_min_spacing, num_turbs)
                test_values = [bndry_x_clsd[1], bndry_y_clsd[1], 0]
                #@test ff.VR_bounary_startup(bndry_x_clsd, bndry_y_clsd, start_dist, turb_min_spacing, num_turbs) == test_values
                @test testing_x[1] == test_values[1]
                @test testing_y[1] == test_values[2]
                @test num_leftover == test_values[3]
            end 

            @testset "One Turbine, Circular Boundary, With Start Distance" begin
                #-- One-turbine circular boundary, zero start distance  --#
                # A discretized 20-point circle
                bndry_x_clsd = [200.00, 195.11, 180.90, 158.78, 130.90, 100.00, 69.10, 41.22, 19.10, 4.89, 0.00, 4.89, 19.10, 41.22, 69.10, 100.00, 130.90, 158.78, 180.90, 195.11, 200.00]
                bndry_y_clsd = [100.00, 130.90, 158.78, 180.90, 195.11, 200.00, 195.11, 180.90, 158.78, 130.90, 100.00, 69.10, 41.22, 19.10, 4.89, 0.00, 4.89, 19.10, 41.22, 69.10, 100.00]
                # Vertices that keep splines injective (4-corners)
                bndry_corner_indcies =[1,6,11,16, 21]  # 20 pt circle, 4 corners
                # Give it one small turbine
                num_turbs = 1
                start_dist = (pi/2)*100 # Top of the circle
                turb_diam = 10.0
                turb_min_spacing = 2*turb_diam

                testing_x, testing_y, num_leftover = ff.VR_boundary_startup(bndry_x_clsd, bndry_y_clsd, start_dist, turb_min_spacing, num_turbs)
                test_values = [100.00, 200.0, 0]
                @test test_values[1] ≈ bndry_x_clsd[6] atol=1E-10 # Should line up with top coordinate
                @test test_values[2] ≈ bndry_y_clsd[6] atol=1E-10
                @test test_values[3] == 0                # No leftover turbines
            end

            @testset "One Turbine, Square Boundary, With Start Distance" begin
                # A 4-point sqaure
                bndry_x_clsd = [100.0, 0.0, 0.0, 100.0, 100.0]
                bndry_y_clsd = [100.0, 100.0, 0.0, 0.0, 100.0]
                # Give it one small turbine, at the start
                num_turbs = 1
                start_dist = 100.0 # Top left corner
                turb_diam = 10.0
                turb_min_spacing = 2*turb_diam

                testing_x, testing_y, num_leftover = ff.VR_boundary_startup(bndry_x_clsd, bndry_y_clsd, start_dist, turb_min_spacing, num_turbs)
                test_values_x = bndry_x_clsd[2]
                test_values_y = bndry_y_clsd[2]
                test_num_leftover = 0
                #@test ff.VR_bounary_startup(bndry_x_clsd, bndry_y_clsd, start_dist, turb_min_spacing, num_turbs) == test_values
                @test testing_x[1] == test_values_x
                @test testing_y[1] == test_values_y
                @test num_leftover == test_num_leftover
            end

            @testset "Multi-Turbine, Circular Boundary, Zero Start Distance" begin
                #-- Multi-turbine circular boundary as a square --#
                # A discretized 200-point circle
                num_pts = 400
                circ_radius = 100.0
                circ_center = [100.0, 100.0]
                bndry_x_clsd, bndry_y_clsd = ff.DiscreteCircum(circ_center[1], circ_center[2], circ_radius, num_pts)
                # Vertices that keep splines injective (4-corners)
                bndry_corner_indcies = [1, 51, 101, 151, 201]  # 200 pt circle, 4 corners
                # Give it a few turbines, don't perturb the start
                num_turbs = 8
                start_dist = 0 # Right of the circle
                turb_diam = 10.0
                turb_min_spacing = 2*turb_diam

                # Vertices that keep splines injective
                circ_corners = 1/sqrt(2)
                cc_r = circ_center[2] + (circ_corners * circ_radius) # ~170
                cc_l = circ_center[1] - (circ_corners * circ_radius) # ~ 30
                cc_d = cc_r - cc_l                # ~140
                cc_d2 = cc_d/2                    # ~ 70

                testing_x, testing_y, num_leftover = ff.VR_boundary_startup(bndry_x_clsd, bndry_y_clsd, start_dist, turb_min_spacing, num_turbs)
                test_values_x = [200.0, cc_r, 100.0, cc_l, 0.0, cc_l, 100.0, cc_r]
                test_values_y = [100.0, cc_r, 200.0, cc_r, 100.0, cc_l, 0.0, cc_l]
                test_values_numleftover = [ [200.00, ],
                                            [100.0, ],
                                            0]
                @test test_values_x ≈ testing_x  atol=1E-10
            end

            @testset "Multi-Turbine, Square Boundary, Zero Start Distance" begin
                # A 4-point sqaure
                bndry_x_clsd = [100.0, 0.0, 0.0, 100.0, 100.0]
                bndry_y_clsd = [100.0, 100.0, 0.0, 0.0, 100.0]
                # Give it a few turbines, perturb the start
                num_turbs = 4
                start_dist = 0.0 # Top side mid-point
                turb_diam = 10.0
                turb_min_spacing = 2*turb_diam

                testing_x, testing_y, num_leftover = ff.VR_boundary_startup(bndry_x_clsd, bndry_y_clsd, start_dist, turb_min_spacing, num_turbs)
                test_values_x = [100.0, 0.0, 0.0, 100.0]
                test_values_y = [100.0, 100.0, 0.0, 0.0]
                test_num_leftover = 0
                @test testing_x == test_values_x
                @test testing_y == test_values_y
                @test num_leftover == test_num_leftover
            end

            @testset "Multi-Turbine, Circular Boundary, With Start Distance" begin
                #-- Multi-turbine circular boundary as a square --#
                # A discretized 200-point circle
                num_pts = 400
                circ_radius = 100.0
                circ_center = [100.0, 100.0]
                bndry_x_clsd, bndry_y_clsd = ff.DiscreteCircum(circ_center[1], circ_center[2], circ_radius, num_pts)
                # Vertices that keep splines injective (4-corners)
                bndry_corner_indcies = [1, 51, 101, 151, 201]  # 200 pt circle, 4 corners
                # Give it a few turbines, don't perturb the start
                num_turbs = 8
                start_dist = 100*pi/4 # Top Right of the circle
                turb_diam = 10.0
                turb_min_spacing = 2*turb_diam

                # Vertices that keep splines injective
                circ_corners = 1/sqrt(2)
                cc_r = circ_center[2] + (circ_corners * circ_radius) # ~170
                cc_l = circ_center[1] - (circ_corners * circ_radius) # ~ 30
                cc_d = cc_r - cc_l                # ~140
                cc_d2 = cc_d/2                    # ~ 70

                testing_x, testing_y, num_leftover = ff.VR_boundary_startup(bndry_x_clsd, bndry_y_clsd, start_dist, turb_min_spacing, num_turbs)
                test_values_x = [cc_r, 100.0, cc_l, 0.0, cc_l, 100.0, cc_r, 200.0]
                test_values_y = [cc_r, 200.0, cc_r, 100.0, cc_l, 0.0, cc_l, 100.0]
                test_values_numleftover = [ [200.00, ],
                                            [100.0, ],
                                            0]
                @test test_values_x ≈ testing_x  atol=5E-3
            end

            @testset "Multi-Turbine, Square Boundary, With Start Distance" begin
                # A 4-point sqaure
                bndry_x_clsd = [100.0, 0.0, 0.0, 100.0, 100.0]
                bndry_y_clsd = [100.0, 100.0, 0.0, 0.0, 100.0]
                # Give it a few turbines, perturb the start
                num_turbs = 4
                start_dist = 50.0 # Top side mid-point
                turb_diam = 10.0
                turb_min_spacing = 2*turb_diam

                testing_x, testing_y, num_leftover = ff.VR_boundary_startup(bndry_x_clsd, bndry_y_clsd, start_dist, turb_min_spacing, num_turbs)
                test_values_x = [50.0, 0.0, 50.0, 100.0]
                test_values_y = [100.0, 50.0, 0.0, 50.0]
                test_num_leftover = 0
                @test testing_x == test_values_x
                @test testing_y == test_values_y
                @test num_leftover == test_num_leftover
            end
        end

        @testset "ray casting boundary distances" begin

            # set up turbine location for testing 
            turbinex = [0.0]
            turbiney = [0.0]

            # set up simple square boundary for testing 
            boundaryvertices = [-1.0 -1.0; -1.0 1.0; 1.0 1.0; 1.0 -1.0]
            boundarynormals = ff.boundary_normals_calculator(boundaryvertices)

            # test correct sign (negative) for inside 
            boundarydistance = ff.ray_casting_boundary(boundaryvertices, boundarynormals, turbinex, turbiney)
            @test sign(boundarydistance[1]) == -1

            # test correct sign (positive) for outside 
            turbinex = [-2.0]
            boundarydistance = ff.ray_casting_boundary(boundaryvertices, boundarynormals, turbinex, turbiney)
            @test sign(boundarydistance[1]) == 1

            # test correct distance inside 
            turbinex = [0.0]
            boundarydistance = ff.ray_casting_boundary(boundaryvertices, boundarynormals, turbinex, turbiney)
            @test boundarydistance[1] ≈ -1 atol = 1E-2

            # test correct distance outside 
            turbinex = [2.0]
            turbiney = [0.0]
            boundarydistance = ff.ray_casting_boundary(boundaryvertices, boundarynormals, turbinex, turbiney)
            @test boundarydistance[1] ≈ 1 atol = 1E-2

            # test correct distance on vertex 
            turbinex = [1.0]
            turbiney = [1.0]
            boundarydistance = ff.ray_casting_boundary(boundaryvertices, boundarynormals, turbinex, turbiney)
            @test boundarydistance[1] ≈ 0 atol = 1E-3

            # test correct distance on face 
            turbinex = [1.0]
            turbiney = [0.0]
            boundarydistance = ff.ray_casting_boundary(boundaryvertices, boundarynormals, turbinex, turbiney)
            @test boundarydistance[1] ≈ 0 atol = 1E-2

        end

        @testset "ray casting boundary derivatives" begin

            # set up turbine location for testing 
            turbinex = [0.0]
            turbiney = [0.0]

            # set up simple square boundary for testing 
            boundaryvertices = [-1.0 -1.0; -1.0 1.0; 1.0 1.0; 1.0 -1.0]
            boundarynormals = ff.boundary_normals_calculator(boundaryvertices)

            # up function for getting AD derivatives
            ray_casting_boundary_diff(x) = ff.ray_casting_boundary(boundaryvertices, boundarynormals, x[1], x[2])

            # test correct derivative inside equi-distant to all vertices
            turbinex = [0.0]
            derivfd = FiniteDiff.finite_difference_jacobian(ray_casting_boundary_diff, [turbinex[1] turbiney[1]], Val{:central})
            derivad = ForwardDiff.jacobian(ray_casting_boundary_diff, [turbinex[1] turbiney[1]])
            @test derivad ≈ derivfd atol = 1E-6

            # test correct derivative inside to face
            turbinex = [0.1]
            derivfd = FiniteDiff.finite_difference_jacobian(ray_casting_boundary_diff, [turbinex[1] turbiney[1]], Val{:central})
            derivad = ForwardDiff.jacobian(ray_casting_boundary_diff, [turbinex[1] turbiney[1]])
            @test derivad ≈ derivfd rtol = 1E-6

            # test correct derivative inside to two faces
            turbiney = [0.1]
            derivfd = FiniteDiff.finite_difference_jacobian(ray_casting_boundary_diff, [turbinex[1] turbiney[1]], Val{:central})
            derivad = ForwardDiff.jacobian(ray_casting_boundary_diff, [turbinex[1] turbiney[1]])
            @test derivad ≈ derivfd rtol = 1E-6

            # test correct derivative outside
            turbinex = [2.0]
            turbiney = [0.0]
            derivfd = FiniteDiff.finite_difference_jacobian(ray_casting_boundary_diff, [turbinex[1] turbiney[1]], Val{:central})
            derivad = ForwardDiff.jacobian(ray_casting_boundary_diff, [turbinex[1] turbiney[1]])
            @test derivad ≈ derivfd rtol = 1E-6

            # test correct derivative on vertex 
            turbinex = [1.0]
            turbiney = [1.0]
            derivfd = FiniteDiff.finite_difference_jacobian(ray_casting_boundary_diff, [turbinex[1] turbiney[1]], Val{:central}, relstep=1E-10)
            derivad = ForwardDiff.jacobian(ray_casting_boundary_diff, [turbinex[1] turbiney[1]])
            @test derivad ≈ derivfd rtol = 1E-4

            # test correct derivative on face 
            turbinex = [1.0]
            turbiney = [0.0]
            derivfd = FiniteDiff.finite_difference_jacobian(ray_casting_boundary_diff, [turbinex[1] turbiney[1]], Val{:central}, relstep=1E-10)
            derivad = ForwardDiff.jacobian(ray_casting_boundary_diff, [turbinex[1] turbiney[1]])
            @test derivad ≈ derivfd rtol = 1E-3

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

        @testset "Test AEP" begin

            include("model_sets/model_set_3.jl")
            # println(sum(windprobabilities))
            modelAEP = ff.calculate_aep(turbine_x, turbine_y, turbine_z, rotor_diameter,
            hub_height, turbine_yaw, ct_models, generator_efficiency, cut_in_speed,
            cut_out_speed, rated_speed, rated_power, windresource, power_models, model_set)/1e9
            paperAEP = 1889.3
            # println(modelAEP/paperAEP)
            @test modelAEP/paperAEP ≈ 1 atol=0.1

            end

        @testset "Test AEP on large farm" begin
            # test based on Borselle II and IV wind farms as used in IEA task 37 case studies 3 and 4

            # import model set with wind farm and related details
            include("./model_sets/model_set_7_ieacs3.jl")

            aep = ff.calculate_aep(turbine_x, turbine_y, turbine_z, rotor_diameter,
                hub_height, turbine_yaw, ct_models, generator_efficiency, cut_in_speed,
                cut_out_speed, rated_speed, rated_power, windresource, power_models, model_set,
                rotor_sample_points_y=rotor_points_y,rotor_sample_points_z=rotor_points_z, hours_per_year=365.0*24.0)
            
            @test aep/1E6 ≈ 938573.62950 rtol=1E-6

        end

        @testset "Test AEP states on large farm" begin
            # test based on Borselle II and IV wind farms as used in IEA task 37 case studies 3 and 4

            # import model set with wind farm and related details
            include("./model_sets/model_set_7_ieacs3.jl")
            
            state_aeps = ff.calculate_state_aeps(turbine_x, turbine_y, turbine_z, rotor_diameter,
                            hub_height, turbine_yaw, ct_models, generator_efficiency, cut_in_speed,
                            cut_out_speed, rated_speed, rated_power, windresource, power_models, model_set;
                            rotor_sample_points_y=[0.0], rotor_sample_points_z=[0.0], hours_per_year=365.0*24.0)

            dir_aep = zeros(20)
            for i in 1:20
                for j in 1:20
                    dir_aep[i] += state_aeps[(i-1)*20 + j]
                end
            end

            @test dir_aep./1E6 ≈ [ 20238.63584, 15709.41125, 13286.56833, 13881.04112, 19232.89054,
            32035.08418, 52531.37389, 47035.14700, 46848.21422, 45107.13416,
            53877.69698, 68105.50430, 69587.76656, 73542.89319, 69615.74101,
            66752.31531, 73027.78883, 60187.14103, 59847.98304, 38123.29869] rtol=1E-6
            

        end

        @testset "Test AEP on single turbine" begin
            # import model set with wind farm and related details
            include("./model_sets/model_set_7_ieacs4_single_turbine.jl")
            
            aep = ff.calculate_aep(turbine_x, turbine_y, turbine_z, rotor_diameter,
                            hub_height, turbine_yaw, ct_models, generator_efficiency, cut_in_speed,
                            cut_out_speed, rated_speed, rated_power, windresource, power_models, model_set;
                            rotor_sample_points_y=[0.0], rotor_sample_points_z=[0.0], hours_per_year=365.0*24.0)
            
            @test aep/1E9 ≈ 42.54982024375254 rtol=1E-6

        end
        @testset "calculate_power_from_cp" begin
            generator_efficiency = 0.944
            air_density = 1.1716
            rotor_area = pi*80.0^2/4
            constcp = 0.8
            v0 = 12.0
            yaw = 0.0

            p = ff.calculate_power_from_cp(generator_efficiency, air_density, rotor_area, constcp, v0, yaw)
            @test p ≈ 3.8425979093271587e6 atol=1E-6

            generator_efficiency = 1.0
            air_density = 1.1716
            rotor_area = pi*126.4^2/4
            v0 = 8.0
            yaw = 0.0

            # load yaw power data
            t1yawdata = readdlm("./inputfiles/sowfa_power_data_t1.csv", ',', skipstart=1)
            yaw = t1yawdata[:,1]*pi/180.0

            # load cp data
            cpctdata = readdlm("./inputfiles/NREL5MWCPCT.txt",  ' ', skipstart=1)
            velpoints = cpctdata[:,1]
            cppoints = cpctdata[:,2]
            ctpoints = cpctdata[:,3]
            constcp = linear(velpoints, cppoints, v0)

            powers = zeros(length(yaw))
            for i = 1:length(yaw)
                powers[i] = ff.calculate_power_from_cp(generator_efficiency, air_density, rotor_area, constcp, v0, yaw[i])
            end
            @test powers*1E-6 ≈ t1yawdata[:,2] atol=1E-1

        end

        @testset "calculate_power() PowerModelConstantCP" begin
            generator_efficiency = 0.944
            air_density = 1.1716
            rotor_area = 0.25*pi*80.0^2
            constcp = 0.8
            v0 = 8.0
            cut_in_speed = 4.  # m/s
            cut_out_speed = 25.  # m/s
            rated_speed = 16.  # m/s
            rated_power = 2.0E6  # W
            yaw = 0.0

            power_model = ff.PowerModelConstantCp(constcp)

            p = ff.calculate_power(generator_efficiency, air_density, rotor_area, v0, yaw, cut_in_speed, rated_speed, cut_out_speed, rated_power, power_model)
            @test p ≈ 0.5*rotor_area*constcp*air_density*generator_efficiency*v0^3 atol=1E-6

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

            yaw = 0.0

            # calculate expected power out
            power = 0.5*cp0*air_density*rotor_area*generator_efficiency*v0^3

            # extract velocity and cp points
            velpoints = data[:,1]
            cppoints = data[:,2]

            # intialize power model struct
            power_model = ff.PowerModelCpPoints(velpoints, cppoints)

            # calculated power and test
            p = ff.calculate_power(generator_efficiency, air_density, rotor_area, v0, yaw, cut_in_speed, rated_speed, cut_out_speed, rated_power, power_model)
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

            yaw = 0.0

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
            p = ff.calculate_power(generator_efficiency, air_density, rotor_area, v0, yaw, cut_in_speed, rated_speed, cut_out_speed, rated_power, power_model)

            # test
            @test p ≈ p0 atol=1E-6

        end

        @testset "calculate_turbine_power() PowerModelConstantCp" begin

            include("./model_sets/model_set_2.jl")

            wt_yaw = 0.0

            # test below cut in
            wt_velocity = 1.0
            p = ff.calculate_turbine_power(generator_efficiency[1], cut_in_speed[1], cut_out_speed[1], rated_speed[1], rated_power[1], rotor_diameter[1], wt_velocity, wt_yaw, power_model, air_density)
            
            @test p ≈ 0.0 atol=1E-6

            # region 2
            v0 = wt_velocity = 8.0
            p = ff.calculate_turbine_power(generator_efficiency[1], cut_in_speed[1], cut_out_speed[1], rated_speed[1], rated_power[1], rotor_diameter[1], wt_velocity, wt_yaw, power_model, air_density)
            rotor_area = pi*0.25*rotor_diameter[1]^2
            @test p ≈ 0.5*constcp*air_density*rotor_area*generator_efficiency[1]*v0^3 atol=1E-6

            # above rated
            wt_velocity = 20.0
            p = ff.calculate_turbine_power(generator_efficiency[1], cut_in_speed[1], cut_out_speed[1], rated_speed[1], rated_power[1], rotor_diameter[1], wt_velocity, wt_yaw, power_model, air_density)
            @test p ≈ rated_power[1] atol=1E-6

            # above cut out
            wt_velocity = 30.0
            p = ff.calculate_turbine_power(generator_efficiency[1], cut_in_speed[1], cut_out_speed[1], rated_speed[1], rated_power[1], rotor_diameter[1], wt_velocity, wt_yaw, power_model, air_density)
            @test p ≈ 0.0 atol=1E-6

        end

        @testset "calculate_turbine_power() PowerModelPowerCurveCubic" begin

            include("./model_sets/model_set_7_ieacs4_single_turbine.jl")

            wt_yaw = 0.0

            # test below cut in
            wt_velocity = 1.0
            p = ff.calculate_turbine_power(generator_efficiency[1], cut_in_speed[1], cut_out_speed[1], rated_speed[1], rated_power[1], rotor_diameter[1], wt_velocity, wt_yaw, power_model, air_density)
            
            @test p ≈ 0.0 atol=1E-6

            # region 2
            speeds = [ 0.90,  1.98,  3.18,  4.40,
            5.64,  6.87,  8.11,  9.35,
            10.59, 11.83, 13.07, 14.31,
            15.56, 16.80, 18.04, 19.28,
            20.52, 21.77, 23.01, 24.25]

            wt_yaw = 0.0

            wt_velocity = speeds[4]
            p = ff.calculate_turbine_power(generator_efficiency[1], cut_in_speed[1], cut_out_speed[1], rated_speed[1], rated_power[1], rotor_diameter[1], wt_velocity, wt_yaw, power_model, air_density)
            rotor_area = pi*0.25*rotor_diameter[1]^2
            @test p ≈ 0.0018658892128279934*1E6 rtol=1E-14


            wt_velocity = speeds[5]
            p = ff.calculate_turbine_power(generator_efficiency[1], cut_in_speed[1], cut_out_speed[1], rated_speed[1], rated_power[1], rotor_diameter[1], wt_velocity, wt_yaw, power_model, air_density)
            @test p ≈ 0.1285989504373177*1E6 rtol=1E-14


            wt_velocity = speeds[6]
            p = ff.calculate_turbine_power(generator_efficiency[1], cut_in_speed[1], cut_out_speed[1], rated_speed[1], rated_power[1], rotor_diameter[1], wt_velocity, wt_yaw, power_model, air_density)
            @test p ≈ 0.6892100000000001*1E6 rtol=1E-14


            wt_velocity = speeds[7]
            p = ff.calculate_turbine_power(generator_efficiency[1], cut_in_speed[1], cut_out_speed[1], rated_speed[1], rated_power[1], rotor_diameter[1], wt_velocity, wt_yaw, power_model, air_density)
            @test p ≈ 2.024097113702623*1E6 rtol=1E-14


            wt_velocity = speeds[8]
            p = ff.calculate_turbine_power(generator_efficiency[1], cut_in_speed[1], cut_out_speed[1], rated_speed[1], rated_power[1], rotor_diameter[1], wt_velocity, wt_yaw, power_model, air_density)
            @test p ≈ 4.4644424198250725*1E6 rtol=1E-14


            wt_velocity = speeds[9]
            p = ff.calculate_turbine_power(generator_efficiency[1], cut_in_speed[1], cut_out_speed[1], rated_speed[1], rated_power[1], rotor_diameter[1], wt_velocity, wt_yaw, power_model, air_density)
            @test p ≈ 8.343766151603498*1E6 rtol=1E-14

            # above rated
            wt_velocity = 20.0
            p = ff.calculate_turbine_power(generator_efficiency[1], cut_in_speed[1], cut_out_speed[1], rated_speed[1], rated_power[1], rotor_diameter[1], wt_velocity, wt_yaw, power_model, air_density)
            @test p ≈ rated_power[1] atol=1E-6

            # above cut out
            wt_velocity = 30.0
            p = ff.calculate_turbine_power(generator_efficiency[1], cut_in_speed[1], cut_out_speed[1], rated_speed[1], rated_power[1], rotor_diameter[1], wt_velocity, wt_yaw, power_model, air_density)
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
            locz = 80.0
            @test ff.adjust_for_wind_shear(locz, point_velocity_no_shear, reference_height, ground_height, model) == point_velocity_no_shear

            # test at ground height
            locz = 0.0
            @test ff.adjust_for_wind_shear(locz, point_velocity_no_shear, reference_height, ground_height, model) == 0.0

            # test below ground height
            locz = -10.0
            @test ff.adjust_for_wind_shear(locz, point_velocity_no_shear, reference_height, ground_height, model) == 0.0

            # test at 40 meters
            locz = 40.0
            u = point_velocity_no_shear*((locz - ground_height)/(reference_height-ground_height))^shear_exp
            @test ff.adjust_for_wind_shear(locz, point_velocity_no_shear, reference_height, ground_height, model) == u

            # test at 10 meters
            loc =[0.0, 0.0, 10.0]
            u = point_velocity_no_shear*((locz - ground_height)/(reference_height-ground_height))^shear_exp
            @test ff.adjust_for_wind_shear(locz, point_velocity_no_shear, reference_height, ground_height, model) == u

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

            deriv_function1(x) = ff.wake_combination_model(x, wind_speed, turb_inflow, old_deficit_sum, model)

            # test correct output
            new_deficit_sum = ff.wake_combination_model(deltav, wind_speed, turb_inflow, old_deficit_sum, model)
            @test new_deficit_sum ≈ result atol=1E-6

            # test correct derivative 
            derivfd = FiniteDiff.finite_difference_derivative(deriv_function1, deltav, Val{:central})
            derivfad = ForwardDiff.derivative(deriv_function1, deltav)
            # derivrad = ReverseDiff.gradient(deriv_function, [deltav])
            @test derivfad ≈ derivfd atol = 1E-6
            # @test derivrad[1] ≈ derivfd atol = 1E-6

            # test correct derivatives at known discontinuity 
            deltav = 4.823068257348535e-185
            wind_speed = 0.9
            turb_inflow = 0.9
            old_deficit_sum = 0.0 
            deriv_function2(x) = ff.wake_combination_model(x, wind_speed, turb_inflow, old_deficit_sum, model)

            derivfd = FiniteDiff.finite_difference_derivative(deriv_function2, deltav, Val{:central})
            derivfad = ForwardDiff.derivative(deriv_function2, deltav)
            @test derivfad ≈ derivfd atol = 1E-6

            deltav = ForwardDiff.Dual(4.823068257348535e-185, ForwardDiff.Partials((-3.199071099983192e-185,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)))
            wind_speed = 0.9
            turb_inflow = ForwardDiff.Dual(0.9, ForwardDiff.Partials((0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)))
            old_deficit_sum = 0.0
            fd_deriv_function(x) = ff.wake_combination_model(x, wind_speed, turb_inflow.value, old_deficit_sum, model)
            derivfd = FiniteDiff.finite_difference_derivative(fd_deriv_function, deltav.value, Val{:central})
            derivfad = ff.wake_combination_model(deltav, wind_speed, turb_inflow, old_deficit_sum, model)
            # derivfad = ForwardDiff.derivative(deriv_function, deltav)
            @test derivfad ≈ derivfd atol = 1E-6
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
            constcp = 0.8
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
            upstream_turbine_id = 1
            downstream_turbine_id = 0
            turbine_definition_id = 1
            cut_in_speed = 0.0
            cut_out_speed = 25.0
            rated_speed = 12.0
            rated_power = 1.0176371581904552e6
            ambient_ti = 0.1

            k_star = 0.07 # adjusted to match experimental data. #TODO improve tests with model results
            horizontal_spread_rate = k_star

            ct_model = ff.ThrustModelConstantCt(ct)
            power_model = ff.PowerModelConstantCp(constcp)

            turbine1 = ff.TurbineDefinition(turbine_definition_id, [rotor_diameter], [hub_height], [cut_in_speed], [rated_speed], [cut_out_speed], [rated_power], [generator_efficiency], [ct_model], [power_model])
            model = ff.JiminezYawDeflection(horizontal_spread_rate)

            dx2p5d = 2.5*rotor_diameter
            dy2p5d_y20 = 0.23136246786632375*rotor_diameter # from [1] figure 21

            dx5p5d = 5.5*rotor_diameter
            dy5p5d_y20 = 0.40616966580976843*rotor_diameter # from [1] figure 21

            dx8d = 8.0*rotor_diameter
            dy8d_20 = 0.5257731958762868*rotor_diameter # from [1] figure 21

            # test deflection at 2.5D with yaw 20 deg
            @test round(ff.wake_deflection_model(dx2p5d, 0.0, hub_height, turbine_x, turbine_yaw, turbine_ct, upstream_turbine_id, rotor_diameter, ambient_ti, model), digits=2) == round(dy2p5d_y20, digits=2)

            # test deflection at 5.5D with yaw 20 deg
            @test round(ff.wake_deflection_model(dx5p5d, 0.0, hub_height, turbine_x, turbine_yaw, turbine_ct, upstream_turbine_id, rotor_diameter, ambient_ti, model), digits=2) == round(dy5p5d_y20, digits=2)

            # test deflection at 8D with yaw 20 deg
            # @test round(ff.wake_deflection_model([dx8d, 0.0, hub_height], model, turbine), digits=2) == round(dy8d_20, digits=2)

        end

        @testset "Gauss Yaw Deflection" begin

            atol = 0.005

            rotor_diameter = 0.15 #[1] p. 509
            hub_height = 0.125 #[1] p. 509
            yaw_20 = 20.0*pi/180.0
            ct = 0.82 # [1] fig. 8
            constcp = 0.8
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
            upstream_turbine_id = 1
            downstream_turbine_id = 0
            turbine_definition_id = 1
            cut_in_speed = 0.0
            cut_out_speed = 25.0
            rated_speed = 12.0
            rated_power = 1.0176371581904552e6
            ambient_ti = 0.1

            ct_model = ff.ThrustModelConstantCt(ct)
            power_model = ff.PowerModelConstantCp([constcp])

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
            @test ff.wake_deflection_model(dx4d, dy4d_20, hub_height, turbine_x, turbine_yaw, turbine_ct, upstream_turbine_id, rotor_diameter, ambient_ti, model) ≈ dy4d_20 atol=atol

            # test deflection at 8D with yaw 20 deg
            @test ff.wake_deflection_model(dx8d, dy8d_20, hub_height, turbine_x, turbine_yaw, turbine_ct, upstream_turbine_id, rotor_diameter, ambient_ti, model) ≈ dy8d_20 atol=atol

        end

        @testset "Gauss Yaw Deflection Variable Spread" begin

            atol = 0.005

            rotor_diameter = 0.15 #[1] p. 509
            hub_height = 0.125 #[1] p. 509
            yaw_20 = 20.0*pi/180.0
            ct = 0.82 # [1] fig. 8
            constcp = 0.8
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
            upstream_turbine_id = 1
            downstream_turbine_id = 0
            turbine_definition_id = 1
            cut_in_speed = 0.0
            cut_out_speed = 25.0
            rated_speed = 12.0
            rated_power = 1.0176371581904552e6

            ct_model = ff.ThrustModelConstantCt(ct)
            power_model = ff.PowerModelConstantCp([constcp])

            turbine_definition = ff.TurbineDefinition(turbine_definition_id, [rotor_diameter], [hub_height], [cut_in_speed], [rated_speed], [cut_out_speed], [rated_power], [generator_efficiency], [ct_model], [power_model])

            turbulence_intensity = 0.07 # this value is just guessed #TODO find data about deflection using this model
            alpha_star = 2.32 #[1] p. 534
            beta_star = 0.154 #[1] p. 534

            model = ff.GaussYawVariableSpreadDeflection(alpha_star, beta_star)

            dx4d = 4.0*rotor_diameter
            dy4d_20 = 0.2684659090909065*rotor_diameter # from [1] figure 21

            dx8d = 8.0*rotor_diameter
            dy8d_20 = 0.34090909090908905*rotor_diameter # from [1] figure 21

            # test deflection at 4D with yaw 20 deg
            @test ff.wake_deflection_model(dx4d, dy4d_20, hub_height, turbine_x, turbine_yaw, turbine_ct, upstream_turbine_id, rotor_diameter, turbulence_intensity, model) ≈ dy4d_20 atol=atol

            # test deflection at 8D with yaw 20 deg
            @test ff.wake_deflection_model(dx8d, dy8d_20, hub_height, turbine_x, turbine_yaw, turbine_ct, upstream_turbine_id, rotor_diameter, turbulence_intensity, model) ≈ dy8d_20 atol=atol

        end

    end

    @testset "Wake Deficit Models" begin

        @testset "Jensen Top Hat Model" begin

            include("./model_sets/model_set_1.jl")

            upstream_turbine_id = 1
            downstream_turbine_id = 0

            deflection_y = 0.0
            deflection_z = 0.0

            centerloss40 = 1. - 4.35/8.1
            centerloss100 = 1. - 5.7/8.1

            overlap_loss = .127


            # test no loss upstream (data from Jensen 1983)
            locx = -1E-12
            locy = 0.0
            locz = hub_height[1]
            @test ff.wake_deficit_model(locx, locy, locz, turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, turbine_local_ti, turbine_ct, turbine_yaw, wakedeficitmodel) == 0.0

            # test max loss at turbine (data from Jensen 1983)
            locx = 0.0
            @test ff.wake_deficit_model(locx, locy, locz, turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, turbine_local_ti, turbine_ct, turbine_yaw, wakedeficitmodel) == (2. * 1/3.0)

            # test centerline loss 40 meters downstream (data from Jensen 1983)
            locx = 40.0
            @test ff.wake_deficit_model(locx, locy, locz, turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, turbine_local_ti, turbine_ct, turbine_yaw, wakedeficitmodel) == centerloss40

            # test centerline loss 100 meters downstream (data from Jensen 1983)
            locx = 100.0
            @test ff.wake_deficit_model(locx, locy, locz, turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, turbine_local_ti, turbine_ct, turbine_yaw, wakedeficitmodel) == centerloss100

            # test wake diameter 40 meters downstream (data from Jensen 1983)
            locx = 40.0
            locy = (alpha*40 + rotor_diameter[1]/2.)
            @test ff.wake_deficit_model(locx, locy, locz, turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, turbine_local_ti, turbine_ct, turbine_yaw, wakedeficitmodel) == centerloss40
            locy = (alpha*40 + rotor_diameter[1]/2. + 1E-12)
            @test ff.wake_deficit_model(locx, locy, locz, turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, turbine_local_ti, turbine_ct, turbine_yaw, wakedeficitmodel) == 0.0
            locy = -(alpha*40 + rotor_diameter[1]/2.)
            @test ff.wake_deficit_model(locx, locy, locz, turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, turbine_local_ti, turbine_ct, turbine_yaw, wakedeficitmodel) == centerloss40
            locy = -(alpha*40 + rotor_diameter[1]/2. + 1E-12)
            @test ff.wake_deficit_model(locx, locy, locz, turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, turbine_local_ti, turbine_ct, turbine_yaw, wakedeficitmodel) == 0.0

            # test wake diameter 100 meters downstream (data from Jensen 1983)
            locx = 100.0
            locy = (alpha*100. + rotor_diameter[1]/2.)
            @test ff.wake_deficit_model(locx, locy, locz, turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, turbine_local_ti, turbine_ct, turbine_yaw, wakedeficitmodel) == centerloss100
            locy = (alpha*100. + rotor_diameter[1]/2. + 1E-12)
            @test ff.wake_deficit_model(locx, locy, locz, turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, turbine_local_ti, turbine_ct, turbine_yaw, wakedeficitmodel) == 0.0
            locy = -(alpha*100. + rotor_diameter[1]/2.)
            @test ff.wake_deficit_model(locx, locy, locz, turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, turbine_local_ti, turbine_ct, turbine_yaw, wakedeficitmodel) == centerloss100
            locy = -(alpha*100. + rotor_diameter[1]/2. + 1E-12)
            @test ff.wake_deficit_model(locx, locy, locz, turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, turbine_local_ti, turbine_ct, turbine_yaw, wakedeficitmodel) == 0.0
            
            # test overlap function
            upstream_turbine_id = 2
            downstream_turbine_id = 1
            turbine_x[1] = 100
            turbine_y[1] = 30
            @test round(ff.wake_deficit_model(locx, locy, locz, turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, turbine_local_ti, turbine_ct, turbine_yaw, wakedeficitmodel), digits = 3) == overlap_loss
            turbine_y[1] = -30
            @test round(ff.wake_deficit_model(locx, locy, locz, turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, turbine_local_ti, turbine_ct, turbine_yaw, wakedeficitmodel), digits = 3) == overlap_loss
            upstream_turbine_id = 1
            downstream_turbine_id = 0

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
            constcp = 0.8
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
            upstream_turbine_id = 1
            downstream_turbine_id = 0
            turbine_definition_id = 1


            alpha = 0.1
            beta = 20.0*pi/180.0

            deflection_y = 0.0
            deflection_z = 0.0

            ct_model = ff.ThrustModelConstantCt(ct)
            power_model = ff.PowerModelConstantCp([constcp])

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
            locx = -1E-12
            locy = 0.0
            locz = hub_height[1]
            @test ff.wake_deficit_model(locx, locy, locz, turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, turbine_local_ti, turbine_ct, turbine_yaw, model) == 0.0

            # test max loss at turbine (data from Jensen 1983)
            locx = 0.0
            @test ff.wake_deficit_model(locx, locy, locz, turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, turbine_local_ti, turbine_ct, turbine_yaw, model) == (2. * aI)

            # test centerline loss 40 meters downstream (data from Jensen 1983)
            locx = 40.0
            @test ff.wake_deficit_model(locx, locy, locz, turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, turbine_local_ti, turbine_ct, turbine_yaw, model) == centerloss40

            # test centerline loss 100 meters downstream (data from Jensen 1983)
            locx = 100.0
            @test ff.wake_deficit_model(locx, locy, locz, turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, turbine_local_ti, turbine_ct, turbine_yaw, model) == centerloss100

            # test wake diameter 40 meters downstream (data from Jensen 1983)
            locx = 40.0
            locy = dy40
            @test ff.wake_deficit_model(locx, locy, locz, turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, turbine_local_ti, turbine_ct, turbine_yaw, model) == 0.0
            @test ff.wake_deficit_model(locx, (locy + 1E-12), locz, turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, turbine_local_ti, turbine_ct, turbine_yaw, model) == 0.0
            @test ff.wake_deficit_model(locx, (locy - 1E1), locz, turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, turbine_local_ti, turbine_ct, turbine_yaw, model) >= 0.0
            @test ff.wake_deficit_model(locx, -(locy), locz, turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, turbine_local_ti, turbine_ct, turbine_yaw, model) == 0.0
            @test ff.wake_deficit_model(locx, -(locy + 1E-12), locz, turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, turbine_local_ti, turbine_ct, turbine_yaw, model) == 0.0
            @test ff.wake_deficit_model(locx, -(locy - 1E1), locz, turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, turbine_local_ti, turbine_ct, turbine_yaw, model) >= 0.0

            # test wake diameter 100 meters downstream (data from Jensen 1983)
            locx = 100.0
            @test ff.wake_deficit_model(locx, dy100, locz, turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, turbine_local_ti, turbine_ct, turbine_yaw, model) == 0.0
            @test ff.wake_deficit_model(locx, (dy100 + 1E-12), locz, turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, turbine_local_ti, turbine_ct, turbine_yaw, model) == 0.0
            @test ff.wake_deficit_model(locx, (dy100 - 1E1), locz, turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, turbine_local_ti, turbine_ct, turbine_yaw, model) >= 0.0
            @test ff.wake_deficit_model(locx, -(dy100), locz, turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, turbine_local_ti, turbine_ct, turbine_yaw, model) == 0.0
            @test ff.wake_deficit_model(locx, -(dy100 + 1E-12), locz, turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, turbine_local_ti, turbine_ct, turbine_yaw, model) == 0.0
            @test ff.wake_deficit_model(locx, -(dy100 - 1E1), locz, turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, turbine_local_ti, turbine_ct, turbine_yaw, model) >= 0.0

            # test value at point in wake 40 m downstream and with theta=15 degrees

            @test ff.wake_deficit_model(40., dy, hub_height[1], turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, turbine_local_ti, turbine_ct, turbine_yaw, model) == loss40attheta

            # test with wec
            model.wec_factor[1] = 1.0
            loss0 = ff.wake_deficit_model(100.0, 50.0, hub_height[1], turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, turbine_local_ti, turbine_ct, turbine_yaw, model)
            model.wec_factor[1] = 3.0
            loss1 = ff.wake_deficit_model(100.0, 50.0, hub_height[1], turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, turbine_local_ti, turbine_ct, turbine_yaw, model)
            @test loss1 > 3.0*loss0

            model.wec_factor[1] = 1.0
            loss0 = ff.wake_deficit_model(100.0, 0.0, hub_height[1], turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, turbine_local_ti, turbine_ct, turbine_yaw, model)
            model.wec_factor[1] = 3.0
            loss1 = ff.wake_deficit_model(100.0, 0.0, hub_height[1], turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, turbine_local_ti, turbine_ct, turbine_yaw, model)
            @test loss0 == loss1

        end

        @testset "MultiZone Model" begin

            include("./model_sets/model_set_Multizone.jl")

            turbine_inflow_velocities = ff.turbine_velocities_one_direction(turbine_x, turbine_y, turbine_z, rotor_diameter, hub_height, turbine_yaw,
                    sorted_turbine_index, ct_models, rotor_sample_points_y, rotor_sample_points_z, wind_resource,
                    model_set)

            turbine_powers = ff.turbine_powers_one_direction(generator_efficiency, cut_in_speed, cut_out_speed, 
                    rated_speed, rated_power, rotor_diameter, turbine_inflow_velocities, turbine_yaw, air_density, 
                    power_models)

            @test turbine_powers[1] ≈ 790066 atol=0.1
            @test turbine_powers[2] ≈ 1720536.3 atol=0.1

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
            constcp = 0.8
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
            upstream_turbine_id = 1
            downstream_turbine_id = 0
            turbine_definition_id = 1

            deflection_y = 0.0
            deflection_z = 0.0

            windfarmstate = ff.SingleWindFarmState(wind_farm_state_id, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai, turbine_inflow_velcity, [0.0], [ambient_ti],sorted_turbine_index)
            ct_model = ff.ThrustModelConstantCt(ct)
            power_model = ff.PowerModelConstantCp([constcp])

            turbine_definition = ff.TurbineDefinition(1, [rotor_diameter], [hub_height], [cut_in_speed], [rated_speed], [cut_out_speed], [rated_power], [generator_efficiency], [ct_model], [power_model])

            model = ff.GaussYaw(horizontal_spread_rate, vertical_spread_rate, alpha_star, beta_star)

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
            @test ff.wake_deficit_model(-1E-12, 0.0, hub_height[1], turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, ambient_ti, turbine_ct, turbine_yaw, model) == 0.0

            # test loss at x1 with no yaw
            @test ff.wake_deficit_model(x1_0, 0.0, hub_height[1], turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, ambient_ti, turbine_ct, turbine_yaw, model) ≈ loss1_0 rtol=rtol

            # test loss at x2 with no yaw
            @test ff.wake_deficit_model(x2_0, 0.0, hub_height[1], turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, ambient_ti, turbine_ct, turbine_yaw, model)  ≈ loss2_0 rtol=rtol

            # test loss at x3 with no yaw
            @test ff.wake_deficit_model(x3_0, 0.0, hub_height[1], turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, ambient_ti, turbine_ct, turbine_yaw, model) ≈ loss3_0 rtol=rtol

            # test loss at x4 with no yaw
            @test ff.wake_deficit_model(x4_0, 0.0, hub_height[1], turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, ambient_ti, turbine_ct, turbine_yaw, model)  ≈ loss4_0 rtol=rtol

            # test with wec
            model.wec_factor[1] = 1.0
            loss0 = ff.wake_deficit_model(x4_0, 0.0, hub_height[1], turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, ambient_ti, turbine_ct, turbine_yaw, model)
            model.wec_factor[1] = 3.0
            loss1 = ff.wake_deficit_model(x4_0, 0.0, hub_height[1], turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, ambient_ti, turbine_ct, turbine_yaw, model)
            @test loss0 ≈ loss1 rtol=1E-6

            model.wec_factor[1] = 1.0
            loss0 = ff.wake_deficit_model(x4_0, 2.0*rotor_diameter, hub_height[1], turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, ambient_ti, turbine_ct, turbine_yaw, model)
            model.wec_factor[1] = 3.0
            loss1 = ff.wake_deficit_model(x4_0, 2.0*rotor_diameter, hub_height[1], turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, ambient_ti, turbine_ct, turbine_yaw, model)
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
            atol = 0.3
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
            # println(ambient_ti)
            upstream_turbine_id = 1
            downstream_turbine_id = 0
            turbine_definition_id = 1

            deflection_y = 0.0
            deflection_z = 0.0

            windfarmstate = ff.SingleWindFarmState(wind_farm_state_id, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai, turbine_inflow_velcity, [0.0], [ambient_ti],sorted_turbine_index)
            ct_model = ff.ThrustModelConstantCt(ct)
            power_model = ff.PowerModelConstantCp([constcp])

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
            @test ff.wake_deficit_model(-1E-12, 0.0, hub_height[1], turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, ambient_ti, turbine_ct, turbine_yaw, model) == 0.0

            # test loss at x1 with no yaw
            @test ff.wake_deficit_model(x1_0, 0.0, hub_height[1], turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, ambient_ti, turbine_ct, turbine_yaw, model) ≈ loss1_0 atol=atol

            deflection = [0.0, 0.0]
            # test loss at x2 with no yaw
            @test ff.wake_deficit_model(x2_0, 0.0, hub_height[1], turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, ambient_ti, turbine_ct, turbine_yaw, model)  ≈ loss2_0 rtol=rtol

            # test loss at x3 with no yaw
            @test ff.wake_deficit_model(x3_0, 0.0, hub_height[1], turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, ambient_ti, turbine_ct, turbine_yaw, model) ≈ loss3_0 rtol=rtol

            # test loss at x4 with no yaw
            @test ff.wake_deficit_model(x4_0, 0.0, hub_height[1], turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, ambient_ti, turbine_ct, turbine_yaw, model)  ≈ loss4_0 rtol=rtol

            # test with wec
            model.wec_factor[1] = 1.0
            loss0 = ff.wake_deficit_model(x4_0, 0.0, hub_height[1], turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, ambient_ti, turbine_ct, turbine_yaw, model)
            model.wec_factor[1] = 3.0
            loss1 = ff.wake_deficit_model(x4_0, 0.0, hub_height[1], turbine_x, turbine_y, turbine_z, deflection_y, deflection_z, upstream_turbine_id, downstream_turbine_id, hub_height, rotor_diameter, turbine_ai, ambient_ti, turbine_ct, turbine_yaw, model)
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
            turbine_velocities, turbine_ct, turbine_ai, turbine_local_ti = ff.turbine_velocities_one_direction(turbine_x, turbine_y, turbine_z, rotor_diameter, hub_height, turbine_yaw,
            sorted_turbine_index, ct_model, rotor_sample_points_y, rotor_sample_points_z, windresource,
            model_set, velocity_only=false)

            # load horns rev ti ata
            data = readdlm("inputfiles/horns_rev_ti_by_row_niayifar.txt", ',', skipstart=1)

            # freestream
            ti_dst = ff.calculate_local_ti(turbine_x, turbine_y, ambient_ti, rotor_diameter, hub_height, turbine_yaw, turbine_local_ti, sorted_turbine_index,
            turbine_inflow_velcities, turbine_ct, localtimodel, turbine_id=(1+ 4*10), tol=1E-6)
            @test ti_dst  == data[1,2]

            # row 2
            ti_dst = ff.calculate_local_ti(turbine_x, turbine_y, ambient_ti, rotor_diameter, hub_height, turbine_yaw, turbine_local_ti, sorted_turbine_index,
            turbine_inflow_velcities, turbine_ct, localtimodel, turbine_id=(2+ 4*10), tol=1E-6)
            @test ti_dst  ≈ data[2,2] atol=atol

            # row 3
            ti_dst = ff.calculate_local_ti(turbine_x, turbine_y, ambient_ti, rotor_diameter, hub_height, turbine_yaw, turbine_local_ti, sorted_turbine_index,
            turbine_inflow_velcities, turbine_ct, localtimodel, turbine_id=(3+ 4*10), tol=1E-6)
            @test ti_dst  ≈ data[3,2] atol=atol

            # row 4
            ti_dst = ff.calculate_local_ti(turbine_x, turbine_y, ambient_ti, rotor_diameter, hub_height, turbine_yaw, turbine_local_ti, sorted_turbine_index,
            turbine_inflow_velcities, turbine_ct, localtimodel, turbine_id=(4+ 4*10), tol=1E-6)
            @test ti_dst  ≈ data[4,2] atol=atol

            # row 5
            ti_dst = ff.calculate_local_ti(turbine_x, turbine_y, ambient_ti, rotor_diameter, hub_height, turbine_yaw, turbine_local_ti, sorted_turbine_index,
            turbine_inflow_velcities, turbine_ct, localtimodel, turbine_id=(5+ 4*10), tol=1E-6)
            @test ti_dst  ≈ data[5,2] atol=atol

            # row 6
            ti_dst = ff.calculate_local_ti(turbine_x, turbine_y, ambient_ti, rotor_diameter, hub_height, turbine_yaw, turbine_local_ti, sorted_turbine_index,
            turbine_inflow_velcities, turbine_ct, localtimodel, turbine_id=(6+ 4*10), tol=1E-6)
            @test ti_dst  ≈ data[6,2] atol=atol

        end

        @testset "Local TI Model No Local TI" begin

            # load model set
            include("./model_sets/model_set_2.jl")

            # calculate turbine inflow velocities
            turbine_velocities, turbine_ct, turbine_ai, turbine_local_ti = ff.turbine_velocities_one_direction(turbine_x, turbine_y, turbine_z, rotor_diameter, hub_height, turbine_yaw,
            sorted_turbine_index, ct_model, rotor_sample_points_y, rotor_sample_points_z, windresource,
            model_set, velocity_only=false)

            # load horns rev ti ata
            data = readdlm("inputfiles/horns_rev_ti_by_row_niayifar.txt", ',', skipstart=1)

            # freestream
            ti_dst = ff.calculate_local_ti(turbine_x, turbine_y, ambient_ti, rotor_diameter, hub_height, turbine_yaw, turbine_local_ti, sorted_turbine_index,
            turbine_inflow_velcities, turbine_ct, localtimodel, turbine_id=(1+ 4*10), tol=1E-6)
            @test ti_dst == ambient_ti

            # row 2
            ti_dst = ff.calculate_local_ti(turbine_x, turbine_y, ambient_ti, rotor_diameter, hub_height, turbine_yaw, turbine_local_ti, sorted_turbine_index,
            turbine_inflow_velcities, turbine_ct, localtimodel, turbine_id=(2+ 4*10), tol=1E-6)
            @test ti_dst == ambient_ti

            # row 3
            ti_dst = ff.calculate_local_ti(turbine_x, turbine_y, ambient_ti, rotor_diameter, hub_height, turbine_yaw, turbine_local_ti, sorted_turbine_index,
            turbine_inflow_velcities, turbine_ct, localtimodel, turbine_id=(3+ 4*10), tol=1E-6)
            @test ti_dst == ambient_ti

            # row 4
            ti_dst = ff.calculate_local_ti(turbine_x, turbine_y, ambient_ti, rotor_diameter, hub_height, turbine_yaw, turbine_local_ti, sorted_turbine_index,
            turbine_inflow_velcities, turbine_ct, localtimodel, turbine_id=(4+ 4*10), tol=1E-6)
            @test ti_dst == ambient_ti

            # row 5
            ti_dst = ff.calculate_local_ti(turbine_x, turbine_y, ambient_ti, rotor_diameter, hub_height, turbine_yaw, turbine_local_ti, sorted_turbine_index,
            turbine_inflow_velcities, turbine_ct, localtimodel, turbine_id=(5+ 4*10), tol=1E-6)
            @test ti_dst == ambient_ti

            # row 6
            ti_dst = ff.calculate_local_ti(turbine_x, turbine_y, ambient_ti, rotor_diameter, hub_height, turbine_yaw, turbine_local_ti, sorted_turbine_index,
            turbine_inflow_velcities, turbine_ct, localtimodel, turbine_id=(6+ 4*10), tol=1E-6)
            @test ti_dst == ambient_ti

        end

        # @testset "Gaussian TI" begin TODO: get this TI model and tests working

        #         include("model_sets/model_set_5.jl")
        #         ambient_ti = 0.137

        #         x = [2.959e-2,            2.219e-1,            4.290e-1,            6.805e-1,
        #         9.467e-1,            1.287e+0,            1.701e+0,            2.101e+0,
        #         2.441e+0,            2.811e+0,            3.092e+0,            3.388e+0,
        #         3.683e+0,            3.979e+0,            4.364e+0,            4.852e+0,
        #         5.237e+0,            5.740e+0,            6.139e+0,            6.686e+0,
        #         7.411e+0,            8.166e+0,            8.861e+0,            9.408e+0,
        #         9.970e+0] .* rotor_diameter

        #         """paper data from "A new Gaussian-based analytical wake model for wind turbines
        #         considering ambiend turbulence intensities and thrust coefficient effects" by Ishihara and
        #         Qian"""
        #         paper_data = [1.625e-1, 1.841e-1, 2.023e-1, 2.114e-1, 2.149e-1, 2.149e-1, 2.081e-1, 1.991e-1,
        #             1.900e-1, 1.821e-1, 1.753e-1, 1.697e-1, 1.629e-1, 1.573e-1, 1.505e-1, 1.426e-1,
        #             1.370e-1, 1.302e-1, 1.234e-1, 1.189e-1, 1.111e-1, 1.032e-1, 9.760e-2, 9.425e-2,
        #             9.090e-2]

        #         TI = zeros(length(x))
        #         for i = 1:length(x)
        #                 loc = [x[i],0.0,hub_height+rotor_diameter/2.0]
        #                 TI[i] = ff.GaussianTI(loc,turbine_x, turbine_y, rotor_diameter, hub_height, turbine_ct, sorted_turbine_index, ambient_ti; div_sigma=2.5, div_ti=1.2)
        #         end

        #         tol = 1E-2
        #         @test TI.-ambient_ti ≈ paper_data atol=tol

        # end

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
            turbine_x = [0.0]
            # test no loss upstream (data from Jensen 1983)
            expected_velocity = wind_speed
            locx = -1.0
            locy = 0.0
            locz = hub_height[1]
            @test ff.point_velocity(locx, locy, locz, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai,
            rotor_diameter, hub_height, turbine_local_ti, sorted_turbine_index, wtvelocities,
            windresource, model_set, wind_farm_state_id=1, downwind_turbine_id=0) == expected_velocity

            # test max loss at turbine (data from Jensen 1983)
            expected_velocity = wind_speed*(1.0 - (2.0 * 1.0/3.0))
            locx = 1E-5
            locy = 0.0
            locz = hub_height[1]
            @test ff.point_velocity(locx, locy, locz, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai,
            rotor_diameter, hub_height, turbine_local_ti, sorted_turbine_index, wtvelocities,
            windresource, model_set, wind_farm_state_id=1, downwind_turbine_id=0) ≈ expected_velocity rtol=rtol

            # test centerline loss 40 meters downstream (data from Jensen 1983)
            expected_velocity = wind_speed*(4.35/8.1)
            locx = 40.0
            locy = 0.0
            locz = hub_height[1]
            @test ff.point_velocity(locx, locy, locz, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai,
            rotor_diameter, hub_height, turbine_local_ti, sorted_turbine_index, wtvelocities,
            windresource, model_set, wind_farm_state_id=1, downwind_turbine_id=0) ≈ expected_velocity rtol=rtol

            # test centerline loss 100 meters downstream (data from Jensen 1983)
            expected_velocity = wind_speed*(5.7/8.1)
            locx = 100.0
            locy = 0.0
            locz = hub_height[1]
            @test ff.point_velocity(locx, locy, locz, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai,
            rotor_diameter, hub_height, turbine_local_ti, sorted_turbine_index, wtvelocities,
            windresource, model_set, wind_farm_state_id=1, downwind_turbine_id=0) ≈ expected_velocity  rtol=rtol
        end

        @testset "calculate_flow_field" begin 

            include("./model_sets/model_set_0_single_turbine.jl")

            # define how many points should be in the flow field
            npoints = 10

            # define how far off the ground to investigate
            maxheight = 4.0*rotor_diameter[1]./2.0

            # set up point grid for flow field
            xrange = -1*rotor_diameter[1]
            yrange = 0
            zrange = 0.0001:maxheight/npoints:maxheight

            # test wind shear with uchida 2020 data doe uniform inflow
            inflowuniform = ones(10).*wind_speed[1]
            heightuniform = collect((0:4/10:1).*rotor_diameter[1])

            shearexponent = 0.0
            wind_shear_model = ff.PowerLawWindShear(shearexponent)
            wind_resource = ff.DiscretizedWindResource(winddirections, windspeeds, windprobabilities, measurementheights, air_density, ambient_tis, wind_shear_model)

            ffvelocities = ff.calculate_flow_field(xrange, yrange, zrange,
                model_set, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai,
                rotor_diameter, hub_height, turbine_local_ti, sorted_turbine_index, wtvelocities,
                wind_resource)         

            @test all(ffvelocities .== inflowuniform)

            # test with data from Bastankhah 2014
            include("./model_sets/bastankhah2014-single-turb-case-3.jl")

            uh = 9.0
            hh = 70.0

            # set up point grid for flow field
            xrange = 7*rotor_diameter[1]
            yrange = 0

            shearexponent = 0.12539210313906432
            groundheight = 4.842460795576101
            wind_shear_model = ff.PowerLawWindShear(shearexponent, groundheight)
            wind_resource = ff.DiscretizedWindResource(winddirections, windspeeds, windprobabilities, measurementheights, air_density, ambient_tis, wind_shear_model)

            rotor_sample_points_y, rotor_sample_points_z = ff.rotor_sample_points(1)

            # data_b = readdlm("../inputfiles/results-bastankhah-2014/bastankhah-7d-5E-2.csv", ',', skipstart=1)
            data_b = readdlm("./inputfiles/results/bastankhah2014/wake-profile-5E-2-7d-fig-7.csv", ',', skipstart=1)
            zrange_b = data_b[:,2].*rotor_diameter[1]
            u0_b = zeros(length(zrange_b))
            u_b = zeros(length(zrange_b))
            for i in 1:length(zrange_b) 
                u0_b[i] = ff.adjust_for_wind_shear(zrange_b[i], uh, hh, wind_shear_model.ground_height, wind_shear_model)
                u_b[i] = u0_b[i]*(1.0-data_b[i,1])
            end

            wakedeficitmodel = ff.GaussOriginal(0.04)
            model_set_bp2014 = ff.WindFarmModelSet(wakedeficitmodel, wakedeflectionmodel, wakecombinationmodel, localtimodel)

            ffvelocitiesbp2014 = ff.calculate_flow_field(xrange, yrange, zrange_b,
                model_set_bp2014, turbine_x, turbine_y, turbine_z, turbine_yaw,
                rotor_diameter, hub_height, ct_models, rotor_sample_points_y, rotor_sample_points_z,
                wind_resource, shearfirst=false)  

            ffvelocitiesbp2014 = reshape(ffvelocitiesbp2014, (length(u0_b)))

            @test isapprox(ffvelocitiesbp2014, u_b, atol=0.1)

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

    end

    @testset "IO" begin

        @testset "read and write yaml" begin

            loc_yaml = "./inputfiles/iea37-ex-opt3.yaml"
            turbine_x1, turbine_y1, fname_turb1, fname_wr1 = ff.get_turb_loc_YAML(loc_yaml)
            test_yaml = "runtest_temp.yaml"
            ff.write_turb_loc_YAML(test_yaml, turbine_x1, turbine_y1; turbinefile=fname_turb1, windresourcefile=fname_wr1)
            turbine_x2, turbine_y2, fname_turb2, fname_wr2 = ff.get_turb_loc_YAML(test_yaml)
            @test turbine_x1 == turbine_x2
            @test turbine_y1 == turbine_y2
            @test fname_turb1 == fname_turb2
            @test fname_wr1 == fname_wr2
            rm(test_yaml)

        end

        @testset "get_turb_loc_YAML" begin

            test_coords = [[10363.7833 6490.2719];
                            [ 9894.9437 6316.9180];
                            [ 8450.2895 6455.3421];
                            [ 9008.9311 6043.4997];
                            [ 9567.5726 5631.6572];
                            [10126.2142 5219.8148];
                            [ 7862.2807 5665.8933];
                            [ 8537.7355 5093.7148];
                            [ 9213.1903 4521.5362];
                            [ 9888.6451 3949.3577];
                            [ 7274.2718 4876.4446];
                            [ 8066.5399 4143.9299];
                            [ 8858.8079 3411.4153];
                            [ 9651.0760 2678.9006];
                            [ 6686.2630 4086.9958];
                            [ 7371.5049 3416.8405];
                            [ 8056.7467 2746.6851];
                            [ 8741.9886 2076.5297];
                            [ 9427.2305 1406.3743];
                            [ 6098.2541 3297.5471];
                            [ 6750.8566 2665.4498];
                            [ 7403.4592 2033.3525];
                            [ 8056.0622 1401.2557];
                            [ 8708.6700  769.1637];
                            [ 9361.2778  137.0718]]

            test_x = test_coords[:,1]
            test_y = test_coords[:,2]

            file_name = "./inputfiles/iea37-ex-opt3.yaml"
            turbine_x, turbine_y, fname_turb, fname_wr = ff.get_turb_loc_YAML(file_name)
        
            @test turbine_x == test_x
            @test turbine_y == test_y
            @test fname_turb == "iea37-10mw.yaml"
            @test fname_wr == "iea37-windrose-cs3.yaml" 

        end

        @testset "get_turb_atrbt_YAML" begin

            file_name = "./inputfiles/iea37-10mw.yaml"
            turb_ci, turb_co, rated_ws, rated_pwr, turb_diam, turb_hub_height = ff.get_turb_atrbt_YAML(file_name)

            @test turb_ci == 4.0
            @test turb_co == 25.0
            @test rated_ws == 11.0
            @test rated_pwr == 10E6
            @test turb_diam == 198.0
            @test turb_hub_height == 119.0

        end

        @testset "get_wind_rose_YAML" begin

            test_num_speed_bins = 20
            test_wind_dir = [  0.0,  18.0,  36.0,  54.0,  72.0,
                                90.0, 108.0, 126.0, 144.0, 162.0,
                                180.0, 198.0, 216.0, 234.0, 252.0,
                                270.0, 288.0, 306.0, 324.0, 342.0]
            test_wind_dir_freq = [0.0312, 0.0260, 0.0255, 0.0253, 0.0297,
                                    0.0397, 0.0506, 0.0510, 0.0415, 0.0414,
                                    0.0522, 0.0634, 0.0706, 0.0723, 0.0697,
                                    0.0668, 0.0676, 0.0677, 0.0613, 0.0464]
            test_wind_speeds = [0.90,  1.98,  3.18,  4.40,  5.64,
                                6.87,  8.11,  9.35, 10.59, 11.83,
                                13.07, 14.31, 15.56, 16.80, 18.04,
                                19.28, 20.52, 21.77, 23.01, 24.25]
            test_wind_speed_probs = [[0.0156401750, 0.0497090909, 0.0811024638, 0.1050883329, 0.1190301631, 0.1222668202, 0.1159445367, 0.1025108097, 0.0850102571, 0.0663764744, 0.0489239316, 0.0340999223, 0.0225045684, 0.0140748572, 0.0083476946, 0.0046969535, 0.0025085092, 0.0012722583, 0.0006121243, 0.0002800569],
                                        [0.0174786954, 0.0548443199, 0.0883795728, 0.1128729487, 0.1256551886, 0.1264602666, 0.1171083595, 0.1007717749, 0.0810568625, 0.0611779343, 0.0434372134, 0.0290638192, 0.0183497799, 0.0109410613, 0.0061655981, 0.0032843186, 0.0016551606, 0.0007890765, 0.0003560345, 0.0001520147],
                                        [0.0163365064, 0.0541606790, 0.0898797863, 0.1167076179, 0.1309060581, 0.1316880823, 0.1209647499, 0.1024541761, 0.0804744947, 0.0588288237, 0.0401142435, 0.0255527921, 0.0152204718, 0.0084822630, 0.0044241371, 0.0021610670, 0.0009880306, 0.0004230131, 0.0001690052, 0.0000640020],
                                        [0.0131561184, 0.0483094348, 0.0851957668, 0.1153680383, 0.1333342000, 0.1368182314, 0.1269601426, 0.1075589680, 0.0836297527, 0.0598535387, 0.0394933554, 0.0240452164, 0.0135111216, 0.0070060631, 0.0033520302, 0.0014790133, 0.0006010054, 0.0002250020, 0.0000780007, 0.0000250002],
                                        [0.0096451543, 0.0385656170, 0.0720491528, 0.1024286389, 0.1239279828, 0.1330241284, 0.1291130658, 0.1144128306, 0.0930354886, 0.0696031136, 0.0479667675, 0.0304604874, 0.0178212851, 0.0096011536, 0.0047590761, 0.0021690347, 0.0009070145, 0.0003480056, 0.0001230020, 0.0000390006],
                                        [0.0059162662, 0.0266481992, 0.0539554280, 0.0821266957, 0.1060917741, 0.1216774755, 0.1264536904, 0.1202924132, 0.1052787375, 0.0849718237, 0.0632988484, 0.0435159582, 0.0275892415, 0.0161127251, 0.0086563895, 0.0042721922, 0.0019330870, 0.0008010360, 0.0003030136, 0.0001050047],
                                        [0.0033912442, 0.0177062749, 0.0394378395, 0.0647356610, 0.0894924435, 0.1095018841, 0.1212547303, 0.1228158427, 0.1143572337, 0.0980800618, 0.0774955797, 0.0563590579, 0.0376647119, 0.0230836620, 0.0129429319, 0.0066224768, 0.0030842221, 0.0013030938, 0.0004990359, 0.0001720124],
                                        [0.0023011680, 0.0133949778, 0.0319683337, 0.0553040372, 0.0799918394, 0.1020204475, 0.1174695753, 0.1234710134, 0.1190496906, 0.1054576984, 0.0857852623, 0.0639746702, 0.0436311851, 0.0271319806, 0.0153331193, 0.0078465728, 0.0036222644, 0.0015031097, 0.0005580407, 0.0001850135],
                                        [0.0022121748, 0.0129870260, 0.0311634619, 0.0541482777, 0.0786432128, 0.1007199569, 0.1164822021, 0.1230117179, 0.1192094175, 0.1061763879, 0.0868758632, 0.0651921502, 0.0447545356, 0.0280232138, 0.0159512601, 0.0082236497, 0.0038263023, 0.0016001264, 0.0005990473, 0.0002000158],
                                        [0.0025653617, 0.0142190049, 0.0329346438, 0.0558158700, 0.0795232128, 0.1003331470, 0.1147461792, 0.1202989622, 0.1162293883, 0.1037036222, 0.0854610500, 0.0649841628, 0.0455144175, 0.0292971309, 0.0172874375, 0.0093243147, 0.0045836463, 0.0020472887, 0.0008281168, 0.0003030427],
                                        [0.0031360035, 0.0159851152, 0.0352022647, 0.0575904289, 0.0799035691, 0.0988156210, 0.1114316581, 0.1158900848, 0.1118007762, 0.1003191021, 0.0838078185, 0.0651768566, 0.0471470871, 0.0316821383, 0.0197453185, 0.0113936460, 0.0060739437, 0.0029859555, 0.0013514325, 0.0005611796],
                                        [0.0041012892, 0.0187200135, 0.0385729355, 0.0602413135, 0.0807017228, 0.0971789375, 0.1075052192, 0.1105016223, 0.1061861613, 0.0957167648, 0.0810740214, 0.0645757898, 0.0483707934, 0.0340633188, 0.0225360739, 0.0139972258, 0.0081525383, 0.0044485678, 0.0022718220, 0.0010838693],
                                        [0.0054756790, 0.0223183746, 0.0429681044, 0.0639713947, 0.0826069320, 0.0966541969, 0.1046128113, 0.1059048592, 0.1009159518, 0.0908560068, 0.0774577706, 0.0626132420, 0.0480261214, 0.0349634170, 0.0241602941, 0.0158441129, 0.0098576243, 0.0058162187, 0.0032531563, 0.0017237321],
                                        [0.0067183022, 0.0254914731, 0.0470020641, 0.0678182802, 0.0854101120, 0.0978798020, 0.1041251678, 0.1039418048, 0.0979800004, 0.0875553596, 0.0743522174, 0.0600949881, 0.0462746238, 0.0339672552, 0.0237750746, 0.0158704234, 0.0101030039, 0.0061331436, 0.0035490271, 0.0019578766],
                                        [0.0073336894, 0.0270762337, 0.0490933339, 0.0699511121, 0.0871766148, 0.0989912795, 0.1044502138, 0.1035063232, 0.0969411732, 0.0861455495, 0.0728148481, 0.0586384528, 0.0450362075, 0.0330111213, 0.0231012718, 0.0154379222, 0.0098527350, 0.0060050281, 0.0034950005, 0.0019418896],
                                        [0.0070345567, 0.0264421993, 0.0484514859, 0.0695431345, 0.0871395441, 0.0993513421, 0.1051276529, 0.1043542724, 0.0977875508, 0.0868400094, 0.0732627741, 0.0588089740, 0.0449612558, 0.0327574721, 0.0227516116, 0.0150668944, 0.0095129807, 0.0057252195, 0.0032848635, 0.0017962062],
                                        [0.0066654512, 0.0255148345, 0.0472692085, 0.0683955035, 0.0862672071, 0.0989039169, 0.1051536540, 0.1048000465, 0.0985292733, 0.0877217059, 0.0741353646, 0.0595633298, 0.0455362312, 0.0331449430, 0.0229744701, 0.0151660553, 0.0095343801, 0.0057078060, 0.0032525879, 0.0017640306],
                                        [0.0068309869, 0.0261352097, 0.0483737224, 0.0698831692, 0.0879405691, 0.1005169558, 0.1064566396, 0.1056003878, 0.0987293423, 0.0873296760, 0.0732571019, 0.0583643286, 0.0442026242, 0.0318395494, 0.0218178978, 0.0142237952, 0.0088208962, 0.0052026062, 0.0029182665, 0.0015562753],
                                        [0.0085671154, 0.0309959797, 0.0553338752, 0.0776406096, 0.0951382271, 0.1059877710, 0.1094344382, 0.1058265995, 0.0964406128, 0.0831404615, 0.0679673172, 0.0527741517, 0.0389584518, 0.0273611122, 0.0182884589, 0.0116363811, 0.0070484996, 0.0040643244, 0.0022303731, 0.0011652398],
                                        [0.0119334560, 0.0399546154, 0.0677356450, 0.0909862235, 0.1069768745, 0.1143638709, 0.1132292570, 0.1048797399, 0.0914354666, 0.0753167464, 0.0587657923, 0.0435065370, 0.0305975533, 0.0204580678, 0.0130120395, 0.0078762611, 0.0045374548, 0.0024883462, 0.0012997031, 0.0006463497]]

            # convert to flow farm standard format
            test_freq = zeros(length(test_wind_dir_freq)*length(test_wind_speed_probs))
            test_speed = zeros(length(test_wind_dir_freq)*length(test_wind_speed_probs))
            test_dir = zeros(length(test_wind_dir_freq)*length(test_wind_speed_probs))

            for i in 1:length(test_wind_dir_freq)
                for j in 1:length(test_wind_speed_probs)
                    test_freq[(i-1)*test_num_speed_bins+j] = test_wind_dir_freq[i]*test_wind_speed_probs[i][j]
                    test_speed[(i-1)*test_num_speed_bins+j] = test_wind_speeds[j]
                    test_dir[(i-1)*test_num_speed_bins+j] = test_wind_dir[i]
                end
            end
            test_num_speed_bins = 20
            test_min_speed = 0.0
            test_max_speed = 25.0

            file_name = "./inputfiles/iea37-windrose-cs3.yaml"
            directions, speeds, frequencies, ti = ff.get_wind_rose_YAML(file_name)

            @test directions == test_dir
            @test speeds == test_speed
            @test frequencies == test_freq
            @test ti == 0.075
            
        end

        @testset "get_boundary_yaml" begin


            @testset "test single region" begin
                boundary_file_name = string("./inputfiles/iea37-boundary-cs3.yaml")
                boundary_vertices = ff.get_boundary_yaml(boundary_file_name)

                boundary_vertices_correct = [10363.8 6490.3; 9449.7 1602.2; 9387.0 1056.6; 9365.1 625.5; 9360.8 360.2; 9361.5 126.9; 9361.3 137.1; 7997.6 1457.9; 6098.3 3297.5;
                8450.3 6455.3; 8505.4 6422.3; 9133.0 6127.4; 9332.8 6072.6; 9544.2 6087.1; 9739.0 6171.2; 9894.9 6316.9; 10071.8 6552.5; 10106.9 6611.1]
                
                @test boundary_vertices ≈ boundary_vertices_correct atol=1E-6
            end

            @testset "test multiple regions" begin
                boundary_file_name = string("./inputfiles/iea37-boundary-cs4.yaml")
                boundary_vertices = ff.get_boundary_yaml(boundary_file_name)

                boundary_vertices_a = [10363.8 6490.3; 9449.7 1602.2; 9387.0 1056.6; 9365.1 625.5; 9360.8 360.2; 9361.5 126.9; 9361.3 137.1; 7997.6 1457.9; 6098.3 3297.5;
                8450.3 6455.3; 8505.4 6422.3; 9133.0 6127.4; 9332.8 6072.6; 9544.2 6087.1; 9739.0 6171.2; 9894.9 6316.9; 10071.8 6552.5; 10106.9 6611.1]
                boundary_normals_a = [0.9829601758936983 -0.1838186405319916; 0.9934614633172962 -0.11416795042154541; 0.9987121579438882 -0.050734855622757584; 
                    0.9998686751666075 -0.01620593781838486; 0.9999954987444023 0.0030004151269687495; -0.9998078216567232 -0.019604074934516894; -0.6957179389375846 -0.718315076718037; 
                    -0.6957275377423737 -0.7183057797532565; -0.8019887481131871 0.5973391397688945; 0.5138086803485797 0.8579047965820281; 0.4252760929807897 0.905063668886888; 
                    0.2645057513093967 0.9643841078762402; -0.0684295708121141 0.9976559496331737; -0.39636379138742883 0.9180935381958544; -0.6828023205475376 0.7306031693435896; 
                    -0.7996740386176392 0.6004343694034798; -0.8578802011411015 0.5138497450520954; 0.42552559023380465 0.9049463918134445]
                boundary_vertices_b = [5588.4 3791.3; 4670.7 4680.2; 7274.9 7940.8; 7369.9 7896.2; 7455.1 7784.3; 7606.5 7713.0; 7638.9 7708.4; 8297.1 7398.9]
                boundary_normals_b = [-0.6957460043611584 -0.7182878931288504; -0.7813688797257963 0.6240694462926818; 0.4249708760634733 0.9052070230051488; 0.7956275395848184 0.6057861159305391; 
                    0.4260560153872896 0.9046967844268629; 0.14056568619461773 0.9900713549359138; 0.4255255464063141 0.9049464124220882; 0.7996806883794807 -0.6004255129763556]
                boundary_vertices_c = [3267.1 10100.6; 4164.1 9586.6; 5749.8 9068.6; 6054.7 8925.3; 1468.5 7781.7; 107.4 9100.0]
                boundary_normals_c = [0.49718026396417986 0.8676472699919642; 0.31052117525343714 0.9505664625470563; 0.42535384615162936 0.9050271297392228; 0.24194817066179167 -0.9702891747893577; 
                    -0.6957228969594285 -0.7183102746351193; -0.30189947425802094 0.9533397649540959]
                boundary_vertices_d = [6764.9 8399.7; 4176.8 5158.6; 2047.8 7220.7]
                boundary_normals_d = [0.7814306689309158 -0.6239920749930895; -0.6957310325444781 -0.7183023947855072; -0.24248239299288069 0.9701558066045093]
                boundary_vertices_e = [8953.7 11901.5; 7048.3 9531.5; 6127.7 9962.7; 4578.1 10464.9; 4524.1 10498.7]
                boundary_normals_e = [0.7793586677376737 -0.6265780613955122; -0.4241667101838764 -0.9055841219742026; -0.30829751674447764 -0.9512899879475178; -0.5305632140423848 -0.847645371546978; -0.3019099610801309 0.9533364439695956]

                @test boundary_vertices[1] ≈ boundary_vertices_a atol=1E-6
                @test boundary_vertices[2] ≈ boundary_vertices_b atol=1E-6
                @test boundary_vertices[3] ≈ boundary_vertices_c atol=1E-6
                @test boundary_vertices[4] ≈ boundary_vertices_d atol=1E-6
            end
        end
    end
end
