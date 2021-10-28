using BenchmarkTools
using FLOWFarm; const ff=FLOWFarm
using DelimitedFiles
using ForwardDiff
using Profile
# using ProfileView

function benchmark_point_in_polygon()


    # vertices = [0.0 0.0; 0.0 10.0; 10.0 10.0; 10.0 0.0] # square 
    # vertices = [0.0 0.0; 0.0 10.0; 10.0 10.0; 10.0 6.0; 5.0 6.0; 5.0 4.0; 10.0 4.0; 10.0 0.0] # square with cut out
    # vertices = [0.0 0.0; 0.0 10.0; 4.0 4.0; 10.0 2.0] # almost a triangle 
    vertices = ff.star_boundary(5, 5.0, 10.0) .+ 5.0
    # println(vertices)
    normals = ff.boundary_normals_calculator(vertices)
    step = 5E-0
    pointx = -5:step:15.0
    pointy = -5:step:15.0
    
    # set up function for getting AD derivatives
    pointinpolygon_diff(x) = ff.pointinpolygon(x, vertices, normals, return_distance=true)
    # @benchmark($pointinpolygon_diff([$pointx[1], $pointy[1]]))
    points = ([pointx[i], pointy[j]] for i = 1:length(pointx) for j = 1:length(pointx))
    # @benchmark [$pointinpolygon_diff(point) for point in $points]
    @profile [pointinpolygon_diff(point) for point in points]
    # @benchmark [ForwardDiff.gradient($pointinpolygon_diff, point) for point in $points]
    # distance = zeros((length(pointx), length(pointy)))
    # xderiv = zeros((length(pointx), length(pointy)))
    # yderiv = zeros((length(pointx), length(pointy)))

    # for i = 1:length(pointx)
    #     for j = 1:length(pointy)
    #         point = [pointx[i], pointy[j]]
    #         distance[j,i] = ff.pointinpolygon(point, vertices, normals, return_distance=true)
    #         deriv = ForwardDiff.gradient(pointinpolygon_diff, point)
    #         xderiv[j,i] = deriv[1]
    #         yderiv[j,i] = deriv[2]
    #     end
    # end

    # fig, ax = plt.subplots(2,2)
    # ff.plotboundary!(ax[1,1], vertices)
    # maxmin = maximum(abs.(distance))
    # cplot1 = ax[1,1].contourf(pointx, pointy, distance, cmap="PuOr", vmin=-maxmin, vmax=maxmin, levels=200)
    # cbar1 = plt.colorbar(cplot1)

    # ff.plotboundary!(ax[1,2], vertices)
    # maxmin = maximum(abs.(xderiv))
    # cplot2 = ax[1,2].contourf(pointx, pointy, xderiv, cmap="PuOr", vmin=-maxmin, vmax=maxmin, levels=200)
    # cbar2 = plt.colorbar(cplot2)

    # ff.plotboundary!(ax[2,1], vertices)
    # maxmin = maximum(abs.(yderiv))
    # cplot3 = ax[2,1].contourf(pointx, pointy, yderiv, cmap="PuOr", vmin=-maxmin, vmax=maxmin, levels=200)
    # cbar3 = plt.colorbar(cplot3)
    
end