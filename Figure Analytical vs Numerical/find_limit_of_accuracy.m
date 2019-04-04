function find_limit_of_accuracy

    addpath('../kenmotsu-drop-simulator/');    
    addpath('../helpers');
    close all;

    V_a = [5 100]; % normalized
    V_r = [];
    
    r0 = 0.5; 
    accuracy = 0.1;
    angle0 = 120;

    options = optimset('Display','off');
    e = @(v1,v2) abs(log(v1/v2)) - log(1+accuracy);        
    
        
    for j = 1:length(V_a)   
        k = nthroot(3*V_a(j)/pi + sqrt(1 + (3*V_a(j)/pi)^2),3);
        hs = k - 1 / k;                               
        numerical_angle = @(a) drop.segment_solve(V_a(j),'angle',a,'radius',1,'height',hs).force;        
        ca_1 = @(a) e(numerical_angle(a),analytical_angle(a,V_a(j),hs));
        a = fsolve(ca_1,angle0,options);
        fprintf('When the volume V is %g, the linear model is accurate up to a = %.0f°\n',V_a(j),a);
        ca_2 = @(a)e(numerical_angle(a),nth_output(2,@analytical_angle,a,V_a(j),hs));
        a2 = fsolve(ca_2,angle0,options);
        fprintf('When the volume V is %g, the quadratic model is accurate up to a = %.0f°\n',V_a(j),a2);        
        
        figure;
        aa = linspace(90,179,20);
        ee = arrayfun(ca_1,aa);
        ee2 = arrayfun(ca_2,aa);
        plot(aa,ee,aa,ee2); 
        title(sprintf('Angle V = %g',V_a(j)));
    end      
    
    for j = 1:length(V_r)                
        k = nthroot(3*V_r(j)/pi + sqrt(1 + (3*V_r(j)/pi)^2),3);
        hs = k - 1 / k;                               
        numerical_radius = @(r) drop.segment_solve(V_r(j),'radius',r,'radius',1,'height',hs).force;                
        
        cr_1 = @(r) e(numerical_radius(r),analytical_radius(r,V_r(j),hs));        
        r = fsolve(cr_1,r0,options);
        fprintf('When the volume V is %g, the linear model is accurate up to r = %g\n',V_r(j),r);                
        cr_2 = @(r)e(numerical_radius(r),nth_output(2,@analytical_radius,r,V_r(j),hs));
        r2 = fsolve(cr_2,r0,options);
        fprintf('When the volume V is %g, the quadratic model is accurate up to r = %g\n',V_r(j),r2);
        
        
        figure;       
        rr = linspace(1e-3,1.5,20);
        ee = arrayfun(cr_1,rr);
        ee2 = arrayfun(cr_2,rr);
        plot(rr,ee,rr,ee2);
        title(sprintf('Radius V = %g',V_r(j)));
    end
   
end

function [F,F2] = analytical_angle(angle,V,hs)
    [he,re,R] = equilibrium_distance_for_angle(angle,V);
    angle2 = 180 - asind(sind(angle)/re);
    f = (cotd(angle)+cscd(angle)*cosd(angle2))^2+1;    
    C = log(cotd(angle/2))+log(cotd(angle2/2))-(cosd(angle)+cosd(angle2))/(cosd(angle)*cosd(angle2)+f);
    k = 2*pi/C;    
    F = -k*(hs-he);
    
    D = log(cotd(angle/2))+log(cotd(angle2/2));
    d2F = -(6*(cosd(angle2)+1))*pi*(cosd(angle2)-1)*(16*cosd(angle)^5*(1/3)-6*cosd(angle)^7+7*cosd(angle)^9*(1/3)-10*cosd(angle2)^3*(1/3)-5*cosd(angle2)^3*cosd(angle)^2*(1/3)-cosd(angle2)^5-(1/3)*cosd(angle)^11+4*cosd(angle)^8*cosd(angle2)-30*cosd(angle)^7*cosd(angle2)^2+10*cosd(angle)^9*cosd(angle2)^2+3*cosd(angle)^7*cosd(angle2)^4+17*cosd(angle2)^2*cosd(angle)^5+64*cosd(angle2)^3*cosd(angle)^4+33*cosd(angle2)*cosd(angle)^4+31*cosd(angle2)^2*cosd(angle)^3-8*cosd(angle2)*cosd(angle)^2-12*cosd(angle2)^2*cosd(angle)-22*cosd(angle)^5*cosd(angle2)^4-3*cosd(angle)^4*cosd(angle2)^5+40*cosd(angle)^3*cosd(angle2)^4+10*cosd(angle)^2*cosd(angle2)^5-6*cosd(angle)*cosd(angle2)^4-22*cosd(angle2)*cosd(angle)^6-152*cosd(angle)^6*cosd(angle2)^3*(1/3)-cosd(angle)^11*cosd(angle2)^2-cosd(angle)^10*cosd(angle2)^3+cosd(angle)*cosd(angle2)^6+(1+cosd(angle2)^6+3*cosd(angle2)^4+3*cosd(angle2)^2-cosd(angle)^2+cosd(angle)^11*cosd(angle2)^3-cosd(angle)^2*cosd(angle2)^6+3*cosd(angle)^3*cosd(angle2)^3+9*cosd(angle)*cosd(angle2)-12*cosd(angle)^3*cosd(angle2)+24*cosd(angle)^2*cosd(angle2)^2+9*cosd(angle)*cosd(angle2)^5-10*cosd(angle)^9*cosd(angle2)^3-3*cosd(angle)^8*cosd(angle2)^4-3*cosd(angle2)^2*cosd(angle)^8+36*cosd(angle)^7*cosd(angle2)^3+21*cosd(angle)^6*cosd(angle2)^4+3*cosd(angle)^5*cosd(angle2)^5+21*cosd(angle2)^2*cosd(angle)^6-48*cosd(angle)^5*cosd(angle2)^3-45*cosd(angle)^4*cosd(angle2)^4-12*cosd(angle)^3*cosd(angle2)^5+3*cosd(angle)^5*cosd(angle2)-45*cosd(angle2)^2*cosd(angle)^4+24*cosd(angle)^2*cosd(angle2)^4+18*cosd(angle)*cosd(angle2)^3)*D-cosd(angle2)-cosd(angle)+38*cosd(angle)^8*cosd(angle2)^3*(1/3)+2*cosd(angle)^3*(1/3))/((1+cosd(angle))*sind(angle2)^2*(cosd(angle)-1)*((cosd(angle)^3*cosd(angle2)-3*cosd(angle)*cosd(angle2)-cosd(angle2)^2-1)*D-cosd(angle)^3-cosd(angle2)*cosd(angle)^2+cosd(angle)+cosd(angle2))^3);
    q = 1/2 * d2F / R;
    F2 = -k*(hs-he)-q*(hs-he).^2;    
end

function [F,F2] = analytical_radius(re,V,hs)
    [he,angle,R] = equilibrium_distance_for_radius(re,V);
    angle2 = 180 - asind(sind(angle)/re);
    f = 1;
    C = log(cotd(angle/2))+log(cotd(angle2/2))-(cosd(angle)+cosd(angle2))/(cosd(angle)*cosd(angle2)+f);
    k = 2*pi/C;    
    F = -k*(hs-he);
        
    D = log(cotd(angle/2))+log(cotd(angle2/2));
    X = cosd(angle)+cosd(angle2);
    Y = (cosd(angle)*cosd(angle2)+1);        
    q = pi*(3*D*Y^3-X^3-3*X*Y^2)/(D*Y-X)^3/R;   
    F2 = -k*(hs-he)-q*(hs-he).^2;    
end

function value = nth_output(N,fcn,varargin)
  [value{1:N}] = fcn(varargin{:});
  value = value{N};
end