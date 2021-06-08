```@setup index
using Plots; gr()
Plots.reset_defaults()
```

# Explanation
This section will explain the fundamental principles behind the wind models used in Flow Farm.

## Wake Deficit Models
### Jensen
#### Top Hat
#### Cosine
### MultiZone
### Bastankah

## Wake Deflection Models
### Jimenez 
### MultiZone
### Gauss

## Wake Combination Models
### Katic
### Whatever Else there is 

## Local Turbulence Intensity

## Power Models

## AEP

## Wind Shear
Wind shear refers to the fact that the wind speed changes with elevation. For wind farms, that change is due to the boundary layer of the wind flowing over the earth. The boundary layer is the region where the flow transitions from stationary at the ground, to some free-stream speed at some distance above the ground. 

Wind shear in FLOWFarm is handled using a power law

``u = u_r \big[\frac{z-z_0}{z_r-z_0}\big]^\phi``

where ``u`` is the wind speed at the desired height (``z``), ``z_r`` is the height of the known speed (``u_r``), and ``z_0`` is the height of the ground (which is zero for flat terrain). The value of ``\phi`` controls how quickly the wind speed transitions from zero to the free-stream velocity.

The models used by FLOWFarm are simple engineering models and do not account for wind shear. To apply wind shear, we first adjust the inflow speed using the wind shear equation and then apply the wake deficit for the given point. 

```jldoctest
using FLOWFarm; const ff = FLOWFarm

# set input values
shear_exponent = 0.15
locz = 40.0 # height above ground in meters
reference_velocity = 8.0 # in m/s
reference_height = 80.0 # height of reference velocity in meters
ground_height = 0.0 # height where velocity goes to zero

# initialize wind shear model instance
wind_shear_model = ff.PowerLawWindShear(shear_exponent)

# adjust wind speed for wind shear
ff.adjust_for_wind_shear(locz, reference_velocity, reference_height, ground_height, wind_shear_model)

# output
7.2100037008866416
```

```@example index
using FLOWFarm; const ff = FLOWFarm
using DataFrames
using StatsPlots

# set input values
shear_exponent = 0.15
locz = 40.0 # height above ground in meters
reference_velocity = 8.0 # in m/s
reference_height = 80.0 # height of reference velocity in meters
ground_height = 0.0 # height where velocity goes to zero

# initialize wind shear model instance
wind_shear_model = ff.PowerLawWindShear(shear_exponent)

h = 0:100
s = zeros(length(h))
for i in length(h)
    s[i] = ff.adjust_for_wind_shear(h[i], reference_velocity, reference_height, ground_height, wind_shear_model)
end

df = DataFrame(:h=h,:s=s)
# Scatter plot with some custom settings
@df df plot(
    :s,
    :h,
    title = "Wind Shear",
    xlabel = "Speed",
    ylabel = "Height"
)
```

## Other Functions
### Overlap 

## Optimization
