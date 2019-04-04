function [ret,r1,R] = equilibrium_distance_for_angle(angle,volume,r2)
    if (nargin < 3)
        r2 = 1;
    end    
            
    % Volume of a spherical segment
    V = @(h,r1) pi*h/6.*(3*r1.^2+3*r2.^2+h.^2);    

    % x is distance from the sphere center to the plane 2
    h = @(x) x-sqrt(x^2+r2^2)*cosd(angle);                        
    x_opt = fzero(@(x) V(h(x),sqrt(x^2+r2^2)*sind(angle))-volume,1);            
    ret = h(x_opt);
    R = sqrt(x_opt^2+r2^2);
    r1 = R*sind(angle);    
end