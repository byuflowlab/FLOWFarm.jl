# FlowFarm.jl


**Summary:** Flow Farm is a tool that enables the user to perform wind farm optimization.

**Authors:** Jared Thomas, Andrew P.J. Stanley

**Features:**


* AEP
* Wind Farm Layouts



**Tools:**

* Wake Deficit Models
* Wake Deflection Models
* Power Models



**Installation:**


```@autodocs
pkg> add FlowFarm.jl
```


**Documentation:**

* For introductory usage help can be found [here](Tutorial.md).
* Help with using specific functions found in the [how-to guide](How_to.md).
* Theory and Methodology surrounding Flow Farm in the [theory](Explanation.md) section.
* Doc Strings and information on the code structure is contained [here](Reference.md).


```@autodocs
Modules = [FlowFarm]
```
**Citing:**
1. N.O. Jensen "A Note on Wind Generator Interaction" *Roskilde: Risø National Laboratory* (1983)
2. I. Katic, k. Hølstrup, N.O. Jensen. "A simple model for cluster efficiency"* European Wind Energy Association Conference and exhibition* (1986) 
3. Gebraad et al. "Wind plant power optimization through yaw control using a parametric model for wake effects—a CFD simulation study" (2014) 
4. Jimenez et al. "Application of a LES technique to characterize the wake deflection of a wind turbine in yaw" *Wind Energy* (2010)
5. Bastankhah "A new analytical model for wind-turbine wakes" *Renewable Energy* (2014)
6. Bastankhah "Experimental and theoretical study of wind turbine wakes in yawed conditions" *Journal of Fluid Mechanics* (2016)
7. Niayifar "Analytical modeling of wind farms: A new approach for power prediction" *Energies* (2016)
8. Thomas "Improving the FLORIS Wind Plant Model for Compatibility with Gradient-Based Optimization" *Wind Engineering* (2017)
9. Thoms "Comparison of Wind Farm Layout Optimization Results Using a Simple Wake Model and Gradient-Based Optimization to Large-Eddy Simulations" *AIAA Scitech 2019 Forum* (2019)