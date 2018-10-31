function [ret,radius1] = equilibrium_distance_for_angle(angle,volume,radius2)
    if (nargin < 3)
        radius2 = 1;
    end    
            
    alpha = 180 - angle;
    h = @(x) x+sqrt((x.*tand(alpha)).^2+x.^2-radius2.^2);
    V = @(x) pi*h(x)/6.*(3*(x*tand(alpha)).^2+3*radius2.^2+h(x).^2);
    x_opt = fzero(@(x) V(x)-volume,1);
    ret = h(x_opt);
    radius1 = x_opt*tand(alpha);  
end