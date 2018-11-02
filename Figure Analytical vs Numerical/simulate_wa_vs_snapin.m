function simulate_wa_vs_snapin

    addpath('../drop-simulator/src');
    addpath('../drop-simulator/src/helpers');
    addpath('../helpers');

    V = [5 20 100]; % normalized
    wa = linspace(1e-3,1); % normalized work of adhesion i.e. angles from 180 to 90
    angles = acosd(wa-1);                  
    
    hs = [];    
    F = {};    
    for j = 1:length(V)
        k = nthroot(3*V(j)/pi + sqrt(1 + (3*V(j)/pi)^2),3);
        hs(j) = k - 1 / k;         
        F{j} = par_sim_fun('dropmodel',@(x) force_ca(x*pi/180,hs(j),V(j)),angles)';        
    end
    if ~exist('output','dir')
        mkdir('output');
    end
    save('output/wa_vs_snapin.mat','V','wa','angles','F','hs');