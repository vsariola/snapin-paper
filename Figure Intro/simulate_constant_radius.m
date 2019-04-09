function simulate_constant_radius

    addpath('../kenmotsu-drop-simulator');    
    addpath('../helpers/');

    V = 5; % normalized
    angle = 150; % deg
    d = drop.segment(V,'angle',angle,'radius',1);
    
    he = d.height;
    re = d.radius1;
    
    %%
    k = nthroot(3*V/pi + sqrt(1 + (3*V/pi)^2),3);
    hs = k - 1 / k; 
    step = 0.01;    
    
    % Advancing CL
    h = hs:-step:(he-0.2);
    F = par_fun(@(x)drop.segment_solve(V,'radius',re,'radius',1,'height',x).force,h)';

    %%
            
    zte = downsample(d.z);
    ute = downsample(d.r);
    
    d_s = drop.segment_solve(V,'radius',re,'radius',1,'height',hs);    
    zts = downsample(d.z);
    uts = downsample(d.r);  
    
    if ~exist('output','dir')
        mkdir('output');
    end
    
    save('output/constant_radius.mat','V','angle','F','h','hs','he','re','zte','ute','zts','uts');