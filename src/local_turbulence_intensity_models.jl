abstract type AbstractLocalTurbulenceIntensityModel end

struct NoLocalTI{TF} <: AbstractLocalTurbulenceIntensityModel
    ct::TF
end
