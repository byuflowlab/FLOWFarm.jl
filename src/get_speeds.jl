using FlowFarm
using CCBlade
using PyPlot
using FLOWMath
using Statistics
using NPZ

const ff=FlowFarm

include("FAST_data.jl")

turb = "high"
ws = 13.0

if turb == "low"
    if ws == 10
        u4 = npzread("flowfields_lowTI/u10_4.npy")
        v4 = npzread("flowfields_lowTI/v10_4.npy")
        w4 = npzread("flowfields_lowTI/w10_4.npy")

        u7 = npzread("flowfields_lowTI/u10_7.npy")
        v7 = npzread("flowfields_lowTI/v10_7.npy")
        w7 = npzread("flowfields_lowTI/w10_7.npy")

        u10 = npzread("flowfields_lowTI/u10_10.npy")
        v10 = npzread("flowfields_lowTI/v10_10.npy")
        w10 = npzread("flowfields_lowTI/w10_10.npy")
    elseif ws == 11
        u4 = npzread("flowfields_lowTI/u11_4.npy").-1.0
        v4 = npzread("flowfields_lowTI/v11_4.npy")
        w4 = npzread("flowfields_lowTI/w11_4.npy")

        u7 = npzread("flowfields_lowTI/u11_7.npy").-1.0
        v7 = npzread("flowfields_lowTI/v11_7.npy")
        w7 = npzread("flowfields_lowTI/w11_7.npy")

        u10 = npzread("flowfields_lowTI/u11_10.npy").-1.0
        v10 = npzread("flowfields_lowTI/v11_10.npy")
        w10 = npzread("flowfields_lowTI/w11_10.npy")
    elseif ws == 12
        u4 = npzread("flowfields_lowTI/u12_4.npy")
        v4 = npzread("flowfields_lowTI/v12_4.npy")
        w4 = npzread("flowfields_lowTI/w12_4.npy")

        u7 = npzread("flowfields_lowTI/u12_7.npy")
        v7 = npzread("flowfields_lowTI/v12_7.npy")
        w7 = npzread("flowfields_lowTI/w12_7.npy")

        u10 = npzread("flowfields_lowTI/u12_10.npy")
        v10 = npzread("flowfields_lowTI/v12_10.npy")
        w10 = npzread("flowfields_lowTI/w12_10.npy")
    elseif ws == 13
        u4 = npzread("flowfields_lowTI/u13_4.npy").+1.0
        v4 = npzread("flowfields_lowTI/v13_4.npy")
        w4 = npzread("flowfields_lowTI/w13_4.npy")

        u7 = npzread("flowfields_lowTI/u13_7.npy").+1.0
        v7 = npzread("flowfields_lowTI/v13_7.npy")
        w7 = npzread("flowfields_lowTI/w13_7.npy")

        u10 = npzread("flowfields_lowTI/u13_10.npy").+1.0
        v10 = npzread("flowfields_lowTI/v13_10.npy")
        w10 = npzread("flowfields_lowTI/w13_10.npy")
    end

elseif turb == "high"
    if ws == 11
        u4 = npzread("flowfields_highTI/u11_4.npy").+1.0
        v4 = npzread("flowfields_highTI/v11_4.npy")
        w4 = npzread("flowfields_highTI/w11_4.npy")

        u7 = npzread("flowfields_highTI/u11_7.npy").+1.0
        v7 = npzread("flowfields_highTI/v11_7.npy")
        w7 = npzread("flowfields_highTI/w11_7.npy")

        u10 = npzread("flowfields_highTI/u11_10.npy").+1.0
        v10 = npzread("flowfields_highTI/v11_10.npy")
        w10 = npzread("flowfields_highTI/w11_10.npy")
    elseif ws == 13
        u4 = npzread("flowfields_lowTI/u13_4.npy")
        v4 = npzread("flowfields_lowTI/v13_4.npy")
        w4 = npzread("flowfields_lowTI/w13_4.npy")

        u7 = npzread("flowfields_lowTI/u13_7.npy")
        v7 = npzread("flowfields_lowTI/v13_7.npy")
        w7 = npzread("flowfields_lowTI/w13_7.npy")

        u10 = npzread("flowfields_lowTI/u13_10.npy")
        v10 = npzread("flowfields_lowTI/v13_10.npy")
        w10 = npzread("flowfields_lowTI/w13_10.npy")
    end

end

x = range(-200.0,stop=200.0,length=161)
y = range(0.0,stop=300.0,length=121)
# XX = zeros(161,121)
# YY = zeros(161,121)
#
# for i=1:161
#    for j=1:121
#        XX[i,j] = x[i]
#        YY[i,j] = y[j]
#    end
# end

# A = range(-200.0,stop=200.0,length=100)
# B = range(0.0,stop=300.0,length=100)
# AA = zeros(length(A),length(B))
# BB = zeros(length(A),length(B))
# UU = zeros(length(A),length(B))
#
# for i=1:length(A)
#    for j=1:length(B)
#        AA[i,j] = A[i]
#        BB[i,j] = B[j]
#        # UU[i,j] = interp2d(akima,x,y,u4,A[i],B[j])[1,1]
#    end
# end

# figure(3)
# pcolormesh(XX,YY,u4,vmin=0.0,vmax=15.0)
# colorbar()
# show()
