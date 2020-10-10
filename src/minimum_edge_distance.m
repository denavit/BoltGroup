function e = minimum_edge_distance(d,hole_type)
% From Table J3.4 of the 2016 AISC Specification 

if d == 0.5
    e = 0.75;
elseif d == 0.625
    e = 0.875;
elseif d == 0.75
    e = 1;
elseif d == 0.875
    e = 1.125;
elseif d == 1
    e = 1.25;
elseif d == 1.125
    e = 1.5;
elseif d == 1.25
    e = 1.625;
elseif d > 1.25
    e = 1.25*d;
else
    error('Not a valid bolt diameter: %g',d);
end

if ~strcmpi(hole_type,'STD')
    error('Not yet implemented for hole_type = %s',hole_type)
end

end