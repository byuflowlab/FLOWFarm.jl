import SpecialFunctions: gamma

abstract type AbstractWindResourceModel end

"""
    DiscritizedWindResource(wind_directions, wind_speeds, wind_probabilities, measurement_heights, air_density, ti_model, wind_shear_model)

    Struct defining a wind resource

# Arguments
- `wind_directions::Array{Float,1}(Nstates)`: an array of wind directions corresponding to each wind farm state in radians
- `wind_speeds::Array{Float,1}(Nstates)`: an array of wind speeds corresponding to each wind farm state in meters/second
- `wind_probabilities::Array{Float,1}(Nstates)`: an array of probabilities corresponding to each wind farm state with values in [0,1]
- `measurement_heights::Array{Float,1}(Nstates)`: an array of measurement heights corresponding to each wind farm state
- `air_density::Float`: the air density
- `ambient_ti::Array{Float,1}`: an array of the ambient turbulence intensity for each wind direction
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

"""
    rediscretize_windrose(windrosein::DiscretizedWindResource, ndirectionbins, nspeedbins)

    Function for re-interpreting a wind rose into a desired number of bins. Returns the new
    wind rose. Currently only works for windroses with a single speed in each direction.

# Arguments
- `windrosein::DiscretizedWindResource`: original wind rose
- `ndirectionbins::Integer`: number of direction bins for the new wind rose
- `start::Float`: direction for first bin in radians
- `averagespeed::Bool`: set whether or not to return the average wind speed as the speed for all bins
"""
function rediscretize_windrose(windrosein::DiscretizedWindResource, ndirectionbins; start=0.0, averagespeed=false)

    # create Akima spline interpolation of windrosein
    splinedirs = [-2*pi.+windrosein.wind_directions; windrosein.wind_directions; 2*pi.+windrosein.wind_directions]
    speedspline = Akima(splinedirs, [windrosein.wind_speeds; windrosein.wind_speeds; windrosein.wind_speeds])
    probabilityspline = Akima(splinedirs, [windrosein.wind_probabilities; windrosein.wind_probabilities; windrosein.wind_probabilities])
    heightspline = Akima(splinedirs, [windrosein.measurement_heights; windrosein.measurement_heights; windrosein.measurement_heights])
    ambienttispline = Akima(splinedirs, [windrosein.ambient_tis; windrosein.ambient_tis; windrosein.ambient_tis])
    
    
    # get new interpolated wind rose attributes
    directionsnew = collect(range(start,2*pi+start-2*pi/ndirectionbins,length=ndirectionbins))
    if averagespeed
        speedsnew = speedspline.(ones(ndirectionbins)*sum(directionsnew)/ndirectionbins)
    else
        speedsnew = speedspline.(directionsnew)
    end
    probabilitiesnew = probabilityspline.(directionsnew)
    heightsnew = heightspline.(directionsnew)
    ambient_tis_new = ambienttispline.(directionsnew)

    # re-normalize probabilites 
    probabilitiesnew = probabilitiesnew./sum(probabilitiesnew)

    # create new wind rose 
    windroseout = DiscretizedWindResource(directionsnew, speedsnew, probabilitiesnew, heightsnew, windrosein.air_density, ambient_tis_new, windrosein.wind_shear_model)

    return windroseout
end


function _dist_weibull(x,L;k=2.0)
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
