function simulate_experimental_comparison
    addpath('../kenmotsu-drop-simulator');
    addpath('../helpers');

    a = 0.5e-3; % Radius of the probe disk, in m (SI-units)
    gamma = 72e-3; % Surface tension, in N/m (SI-units)
    F = 15e-6; % N
    g = 9.81; % m/s^2
    rho = 1000; % kg/m^3
    deltah = 1.2e-6; % m

    % Calculate volume from the weight, density and gravitational constant
    V = F/(g*rho);
    Vt = V/a^3;
        
    r = exp(linspace(log(5e-6),log(800e-6)));    
    k = nthroot(3*Vt/pi + sqrt(1 + (3*Vt/pi)^2),3);
    h = (k - 1 / k) * a + deltah;
        
    F = gamma * par_fun(@(x) drop.segment_solve(V/a^3,'radius',x/a,'radius',1,'height',h/a).force,r)*a';   
    
    if ~exist('output','dir')
        mkdir('output');
    end
    save('output/r_vs_snapin.mat','gamma','r','a','F','h','V','h');