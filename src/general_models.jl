function wind_frame(xlocs, ylocs, wind_direction)

    nTurbines = length(xlocs)
    wd = -pi*(270. - wind_direction)/180. #shift to traditional wind direction coords and to radians
    xw = xlocs.*cos(wd)-ylocs.*sin(wd)
    yw = xlocs.*sin(wd)+ylocs.*cos(wd)
    return xw, yw
end
