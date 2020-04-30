abstract type AbstractLocalTurbulenceIntensityModel end

struct NoLocalTI{TF} <: AbstractLocalTurbulenceIntensityModel
    ct::TF
end


struct CrepsoHernandez{TF} <: AbstractLocalTurbulenceIntensityModel
    initial::TF
    constant::TF
    ai::TF
    downstream::TF
end






# def calculate_turbulence_intensity(self, flow_field_ti, velocity_model, turbine_coord, wake_coord, turbine_wake):
#         """
#         Calculates the turbulence intensity at a specific wind turbine.
#         This method calculates and returns the turbulence intensity at
#         the wind turbine consisting of the ambient turbulence as well
#         as the wake-added turbulence from an upstream turbine, using
#         the approach of Crespo, A. and Herna, J. "Turbulence
#         characteristics in wind-turbine wakes." *J. Wind Eng Ind
#         Aerodyn*. 1996.
#         Args:
#             flow_field_ti: A float that is the ambient turbulence
#                 intensity in the flow field expressed as a decimal
#                 fraction.
#             velocity_model: A
#                 :py:obj:`floris.simulation.wake_velocity.WakeVelocity`
#                 object containing wake model parameters.
#             turbine_coord: A :py:obj:`floris.utilities.Vec3` object
#                 containing the coordinate of the turbine.
#             wake_coord: A :py:obj:`floris.utilities.Vec3` object
#                 containing the coordinate of the upstream turbine.
#             turbine_wake: A :py:class:`floris.simulation.turbine`
#                 object that represents the upstream turbine.
#         Returns:
#             numpy.float64: The turbulence intensity at the current
#             turbine including ambient turbulence and turbulence added
#             by the upstream turbine wake.
#         """
#
#         ti_initial = flow_field_ti
#
#         # user-input turbulence intensity parameters
#         ti_i = velocity_model.ti_initial
#         ti_constant = velocity_model.ti_constant
#         ti_ai = velocity_model.ti_ai
#         ti_downstream = velocity_model.ti_downstream
#
#         # turbulence intensity calculation based on Crespo et. al.
#         ti_calculation = ti_constant \
#             * turbine_wake.aI**ti_ai \
#             * ti_initial**ti_i \
#             * ((turbine_coord.x1 - wake_coord.x1) / self.rotor_diameter)**ti_downstream
#
#         return np.sqrt(ti_calculation**2 + flow_field_ti**2)
