function lc = clear_distance_between_bolt_holes(x1,y1,x2,y2,dh,load_angle)
% Calculates the clear distance between one hole [centered at (x1,y1)] and a
% second hole [centered at (x2,y2)] at an angle measured counter clockwise
% from the positive x axis. 
dist_h2h = sqrt((x2-x1)^2 + (y2-y1)^2);
angle_h2h = atan2(y2-y1,x2-x1);

angle_range  = atan((dh/2)/dist_h2h);
if abs(wrapToPi(load_angle-angle_h2h)) < angle_range
    a  = 1;
    b  = -2*dist_h2h*cos(wrapToPi(load_angle-angle_h2h));
    c  = dist_h2h^2 - (dh/2)^2;
    x  = (-b-sqrt(b^2-4*a*c))/(2*a);
    lc = x - dh/2;
else
    lc = Inf;
end

end

