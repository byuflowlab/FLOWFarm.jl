import SpecialFunctions: gamma

abstract type AbstractWindResourceModel end

"""
    DiscritizedWindResource(wind_directions, wind_speeds, wind_probabilities, measurement_heights, air_density, ti_model, wind_shear_model)

    Struct defining a wind resource

# Arguments
- `wind_directions::Array{Float}(Nstates)`: an array of wind directions corresponding to each wind farm state
- `wind_speeds::Array{Float}(Nstates)`: an array of wind speeds corresponding to each wind farm state
- `wind_probabilities::Array{Float}(Nstates)`: an array of probabilities corresponding to each wind farm state
- `measurement_heights::Array{Float}(Nstates)`: an array of measurement heights corresponding to each wind farm state
- `air_density::Float`: the air density
- `ambient_ti::Array{Float}`: an array of the ambient turbulence intensity for each wind direction
- `wind_shear_model::Array{AbstractWindShearModel}(1)`: contains a struct defining the desired turbulence intensity model
"""
struct DiscretizedWindResource{AF, TF, ASM} <: AbstractWindResourceModel

    wind_directions::AF
    wind_speeds::AF
    wind_probabilities::AF
    measurement_heights::AF
    air_density::TF
    ambient_tis::AF
    wind_shear_model::ASM

end


function _dist_weibull(x,L)
    k = 2.0
    if L < 0.0001
        L = 0.0001
    end
    return (k/L)*(x/L)^(k-1)*exp(-(x/L)^k)
end


function setup_weibull_distribution(wind_directions,wind_probabilities,wind_speeds,nspeeds)
        if nspeeds == 1
            return wind_directions, wind_probabilities, wind_speeds
        end

        ndirections = length(wind_directions)
        dirs = zeros(ndirections*nspeeds)
        freqs = zeros(ndirections*nspeeds)
        speeds = zeros(ndirections*nspeeds)

        #direction loops
        k = 1
        for i = 1:ndirections
            for j = 1:nspeeds
                dirs[k] = wind_directions[i]
                k += 1
            end
        end

        #speed and frequency loops
        m = 1
        for i = 1:ndirections
            avg_speed = wind_speeds[i]
            speed_dist = range(25.0/(2.0*nspeeds)+0.001,stop=25.0-25.0/(2.0*nspeeds),length=nspeeds)
            dspeed = speed_dist[2]-speed_dist[1]
            num_int = 1000

            for j = 1:nspeeds
                speed_int = range(speed_dist[j]-dspeed/2.0,stop=speed_dist[j]+dspeed/2.0,length=num_int)
                k = 2.0
                scale = avg_speed/(gamma(1.0+1.0/k))
                freq_int = _dist_weibull.(speed_int,scale)

                speed_freq = trapz(speed_int,freq_int)
                speeds[m] = speed_dist[j]
                freqs[m] = speed_freq*wind_probabilities[i]
                m += 1
            end
        end

        return dirs, freqs, speeds
end
