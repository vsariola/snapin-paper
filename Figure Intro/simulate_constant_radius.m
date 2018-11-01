function simulate_constant_radius

    addpath('../drop-simulator/src');
    addpath('../drop-simulator/src/helpers');
    addpath('../helpers/');

    V = 5; % normalized
    angle = 150; % deg
    [he,re] = equilibrium_distance_for_angle(angle,V);
    
    %%
    k = nthroot(3*V/pi + sqrt(1 + (3*V/pi)^2),3);
    hs = k - 1 / k; 
    step = 0.01;    
    
    % Advancing CL
    h = hs:-step:(he-0.2);
    F = par_sim_fun('dropmodel',@(x) force(re,x,V),h)';

    %%
    
    [~,zte,ute] = force(re,he,V);
    zte = downsample(zte);
    ute = downsample(ute);
    [~,zts,uts] = force(re,hs,V);
    zts = downsample(zts);
    uts = downsample(uts);  
    
    save('output/constant_radius.mat','V','angle','F','h','hs','he','re','zte','ute','zts','uts');