# calculates the overlap area between a given wake and a rotor area
function overlap_area_func(turbine_y, turbine_z, rotor_diameter, wake_center_y, 
    wake_center_z, wake_diameter; tol=1E-6)
    
    # distance between wake center and rotor center
    if (wake_center_z > (turbine_z + tol)) || (wake_center_z < (turbine_z - tol))
        OVdYd = sqrt((wake_center_y-turbine_y)^2 + (wake_center_z - turbine_z)^2)
    elseif (wake_center_y > (turbine_y + tol))
        OVdYd = wake_center_y - turbine_y
    elseif (turbine_y > (wake_center_y + tol))
        OVdYd = turbine_y - wake_center_y
    else
        OVdYd = 0.0
    end
    
    # find rotor radius
    OVr = rotor_diameter/2.0
    
    # find wake radius
    OVRR = wake_diameter/2.0
    
    # make sure the distance from wake center to turbine hub is positive
    OVdYd = abs(OVdYd)
    
    # determine if there is overlap
    if (OVdYd < (OVr+OVRR)) # if the rotor overlaps the wake zone

        # check that turbine and wake centers are not perfectly aligned
        if (OVdYd > (0.0 + tol))
        
            # check if the rotor is wholly contained in the wake
            if ((OVdYd + OVr) < OVRR + tol) 
                wake_overlap = pi*OVr*OVr
            elseif ((OVdYd + OVRR) < OVr + tol)
                wake_overlap = pi*OVRR*OVRR
            else
                # calculate the distance from the wake center to the chord connecting the lens cusps
                OVL = (-OVr*OVr+OVRR*OVRR+OVdYd*OVdYd)/(2.0*OVdYd)

                OVz = sqrt(OVRR*OVRR-OVL*OVL)
                OVz2 = sqrt(OVr*OVr-(OVdYd-OVL)*(OVdYd-OVL))
            
                wake_overlap = OVRR*OVRR*acos(OVL/OVRR) + OVr*OVr*acos((OVdYd-OVL)/OVr) - OVL*OVz - (OVdYd-OVL)*OVz2
            end
        
        # perfect overlap case where the wake is larger than the rotor
        elseif (OVRR > OVr)
            wake_overlap = pi*OVr*OVr            
        # perfect overlap case where the rotor is larger than the wake
        else
            wake_overlap = pi*OVRR*OVRR
        end
        
    # case with no overlap
    else
        wake_overlap = 0.0
    end
    
    # if ((wake_overlap/(pi*OVr*OVr) > (1.0 + tol)) || (wake_overlap/(pi*OVRR*OVRR) > (1.0 + tol)))
    #     print *, "wake overlap in func: ", wake_overlap/(pi*OVr*OVr)
    #     print *, "wake overlap in func: ", wake_overlap/(pi*OVRR*OVRR)
    # end if
                             
end