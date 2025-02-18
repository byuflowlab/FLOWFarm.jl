# Explanation
This section consists primarily of citations of relevant papers where details about each 
model can be found. It is admittedly incomplete, but will hopefully be expanded over time.

## Wake Deficit Models
### Jensen
#### JensenTopHat
- [1] N. O. Jensen. A note on wind generator interaction. Technical report, Risø National Laboratory, DK-4000 Roskilde, Denmark, November 1983.
#### JensenCosine
- [1] N. O. Jensen. A note on wind generator interaction. Technical report, Risø National Laboratory, DK-4000 Roskilde, Denmark, November 1983.
- [2] J. J. Thomas, S. McOmber, and A. Ning. Wake expansion continuation: Multi-modality reduction in the wind farm layout optimization problem. Wind Energy, May 2021. (in review).
#### MultiZone
Simply put, three jensen top-hat models stacked on top of each other to approximate varrying behavior in different wake zones.
- [3] P. M. O. Gebraad, F. W. Teeuwisse, J. W. van Wingerden, P. A. Fleming, S. D. Ruben, J. R. Marden, and L. Y. Pao. Wind plant power optimization through yaw control using a parametric model for wake effects—a CFD simulation study. Wind Energy, 2014.

### Gaussian
#### GaussOriginal
- [4] M. Bastankhah and F. Port ́e-Agel. A new analytical model for wind-turbine wakes. Renewable Energy, 70:116–123, 2014.
#### GaussYaw
- [5] M. Bastankhah and F. Port ́e-Agel. Experimental and theoretical study of wind turbine wakes in yawed conditions. Journal of Fluid Mechanics, 806:506–541, 2016.
#### GaussYawVariableSpread
- [5] M. Bastankhah and F. Port ́e-Agel. Experimental and theoretical study of wind turbine wakes in yawed conditions. Journal of Fluid Mechanics, 806:506–541, 2016.
- [6] A. Niayifar and F. Port ́e-Agel. Analytical modeling of wind farms: A new approach for power prediction. Energies, 9(9):1–13, 2016.
#### GaussSimple
- [7] N. F. Baker, A. P. J. Stanley, J. J. Thomas, A. Ning, and K. Dykes. Best practices for wake model and optimization algorithm selection in wind farm layout optimization. In AIAA Scitech 2019 Forum, San Diego, CA, Jan. 2019.

## Wake Deflection Models
### Jimenez 
- [8] A ́. Jim ́enez, A. Crespo, and E. Migoya. Application of a LES technique to chracterize the wake deflection of a wind turbine in yaw. Wind Energy, 13:559–572, 2010.
### MultiZone
- [3] P. M. O. Gebraad, F. W. Teeuwisse, J. W. van Wingerden, P. A. Fleming, S. D. Ruben, J. R. Marden, and L. Y. Pao. Wind plant power optimization through yaw control using a parametric model for wake effects—a CFD simulation study. Wind Energy, 2014.
- [9] J. Thomas, P. Gebraad, and A. Ning. Improving the FLORIS wind plant model for compatibility with gradient-based optimization. Wind Engineering, Aug. 2017.
### Gauss
- [5] M. Bastankhah and F. Port ́e-Agel. Experimental and theoretical study of wind turbine wakes in yawed conditions. Journal of Fluid Mechanics, 806:506–541, 2016.

## Wake Combination Models
### LinearFreestreamSuperposition
- [10] P. Lissaman. Energy effectiveness of arbitrary arrays of wind turbines. Journal of Energy, 3:323–328, 1979.
### SumOfSquaresFreestreamSuperposition
- [11] I. Katic, J. Højstrup, and N. Jensen. A simple model for cluster efficiency. In European Wind Energy Association Conference and Exhibition, Rome - Italy, October 1986. European Wind Energy Association.
### SumOfSquaresLocalVelocitySuperposition
- [12] S. Voutsinas, K. Rados, and A. Zervos. On the analysis of wake effects in wind parks. Wind Engineering, 14:204–2019, 1990.
### LinearLocalVelocitySuperposition
- [6] A. Niayifar and F. Port ́e-Agel. Analytical modeling of wind farms: A new approach for power prediction. Energies, 9(9):1–13, 2016.

## Local Turbulence Intensity Models
### LocalTIModelNoLocalTI
Just like it sounds. Local TI is ignored and ambient TI is used anywhere TI is called for in the calculations.
### LocalTIModelMaxTI
- [6] A. Niayifar and F. Port ́e-Agel. Analytical modeling of wind farms: A new approach for power prediction. Energies, 9(9):1–13, 2016.s
### GaussianTI
Not currently connected with the general wind farm models, but it hopefully will be eventually.
- [13] A. Pen ̃a, K. Schaldemose Hansen, S. Ott, and M. P. van der Laan. On wake modeling, wind-farm gradients, and aep predictions at the anholt wind farm. Wind Energy Science, 3(1):191–202, 2018.

## Wake Expansion Continuation (WEC)
- [2] J. J. Thomas, S. McOmber, and A. Ning. Wake expansion continuation: Multi-modality reduction in the wind farm layout optimization problem. Wind Energy, May 2021. (in review).

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

<!-- ```@example
using FLOWFarm; const ff = FLOWFarm
using PyPlot; const plt = PyPlot

# set input values
shear_exponent = 0.15
reference_velocity = 8.0 # in m/s
reference_height = 80.0 # height of reference velocity in meters
ground_height = 0.0 # height where velocity goes to zero

# initialize wind shear model instance
wind_shear_model = ff.PowerLawWindShear(shear_exponent, ground_height)

# initialize heights of interest
h = 0:200

# initialize array for wind speeds at the heights of interest
s = zeros(length(h))

# adjust wind speed for heights of interest based on the the reference speed and height
for i = 1:length(h)
    s[i] = ff.adjust_for_wind_shear(h[i], reference_velocity, reference_height, ground_height, wind_shear_model)
end

# Scatter plot with some custom settings
plt.plot(s, h)
plt.title("Wind Shear")
plt.xlabel("Speed (m/s)")
plt.ylabel("Height (m)")
plt.savefig("windshear.png") # hide
```
![](windshear.png) -->

**Citing:**
- [1] N. O. Jensen. A note on wind generator interaction. Technical report, Risø National Laboratory, DK-4000 Roskilde, Denmark, November 1983.
- [2] J. J. Thomas, S. McOmber, and A. Ning. Wake expansion continuation: Multi-modality reduction in the wind farm layout optimization problem. Wind Energy, May 2021. (in review).
- [3] P. M. O. Gebraad, F. W. Teeuwisse, J. W. van Wingerden, P. A. Fleming, S. D. Ruben, J. R. Marden, and L. Y. Pao. Wind plant power optimization through yaw control using a parametric model for wake effects—a CFD simulation study. Wind Energy, 2014.
- [4] M. Bastankhah and F. Port ́e-Agel. A new analytical model for wind-turbine wakes. Renewable Energy, 70:116–123, 2014.
- [5] M. Bastankhah and F. Port ́e-Agel. Experimental and theoretical study of wind turbine wakes in yawed conditions. Journal of Fluid Mechanics, 806:506–541, 2016.
- [6] A. Niayifar and F. Port ́e-Agel. Analytical modeling of wind farms: A new approach for power prediction. Energies, 9(9):1–13, 2016.
- [7] N. F. Baker, A. P. J. Stanley, J. J. Thomas, A. Ning, and K. Dykes. Best practices for wake model and optimization algorithm selection in wind farm layout optimization. In AIAA Scitech 2019 Forum, San Diego, CA, Jan. 2019.
- [8] A ́. Jim ́enez, A. Crespo, and E. Migoya. Application of a LES technique to chracterize the wake deflection of a wind turbine in yaw. Wind Energy, 13:559–572, 2010.
- [9] J. Thomas, P. Gebraad, and A. Ning. Improving the FLORIS wind plant model for compatibility with gradient-based optimization. Wind Engineering, Aug. 2017.
- [10] P. Lissaman. Energy effectiveness of arbitrary arrays of wind turbines. Journal of Energy, 3:323–328, 1979.
- [11] I. Katic, J. Højstrup, and N. Jensen. A simple model for cluster efficiency. In European Wind Energy Association 
- [12] S. Voutsinas, K. Rados, and A. Zervos. On the analysis of wake effects in wind parks. Wind Engineering, 14:204–2019, 1990.
- [13] A. Pen ̃a, K. Schaldemose Hansen, S. Ott, and M. P. van der Laan. On wake modeling, wind-farm gradients, and aep predictions at the anholt wind farm. Wind Energy Science, 3(1):191–202, 2018.