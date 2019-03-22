function [ret,angle,R] = equilibrium_distance_for_radius(radius1,volume,radius2)
    if (nargin < 3)
        radius2 = 1;
    end        

    % Volume of a spherical segment
    V = @(h) pi*h/6.*(3*radius1.^2+3*radius2.^2+h.^2);
    
    if radius1 > radius2        
        % Parameterize in terms of l1
        h = @(x) x+sqrt(radius1.^2+x.^2-radius2.^2);                                
        x_opt = fzero(@(x) V(h(x))-volume,1);
        ret = h(x_opt);
        angle = mod(180-atand(radius1/x_opt),180);
        R = sqrt(x_opt^2+radius1^2);
    else
        % Parameterize in terms of l2
        h = @(x) x+sqrt(radius2.^2+x.^2-radius1.^2);        
        x_opt = fzero(@(x) V(h(x))-volume,1);
        ret = h(x_opt);
        angle = 180-atand(radius1/(h(x_opt)-x_opt));
        R = sqrt(x_opt^2+radius2^2);
    end            
end