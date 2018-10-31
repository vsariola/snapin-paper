function [zt,ut] = spherical_cap(volume,distance)
    numpoints = 200;    
    
    k = nthroot(3*volume/pi + sqrt(1 + (3*volume/pi)^2),3);
    h = k - 1 / k;
    r = (1+h^2)/(2*h);
    alpha = asind((r-h)/r);
    
    a = linspace(90,alpha,numpoints);
    zt = distance-h+r-r*sind(a);
    ut = r*cosd(a);
        
      