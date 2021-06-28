# function circleshape(h, k, r)
#     theta = LinRange(0, 2*pi, 500)
#     return h .+ r*sin.(theta), k.+ r*cos.(theta)
# end

# function plotlayout!(p, turbinex, turbiney, rotordiameter; linecolor=:black, fillcolor=:blue, markeralpha=1, title="")
#     nturbines = length(turbinex)
#     for i in 1:nturbines
#         plot!(p, circleshape(turbinex[i],turbiney[i],rotordiameter[i]/2.0), seriestype=[:shape],
#             linecolor=linecolor, c=fillcolor, legend=false, aspect_ratio=1, alpha=markeralpha, title=title)
#     end
#     display(p)
# end

function plotlayout!(ax, turbinex, turbiney, rotordiameter; linecolor=:black, fillcolor=:blue, markeralpha=1, title="")
    nturbines = length(turbinex)
    for i in 1:nturbines
        circle = plt.Circle((turbinex[i], turbiney[i]), rotordiameter[i])
        ax.add_patch(circle)
    end
end