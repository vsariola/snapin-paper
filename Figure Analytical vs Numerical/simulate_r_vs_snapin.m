function simulate_r_vs_snapin

    addpath('../kenmotsu-drop-simulator/');    
    addpath('../helpers');

    V = [5 25 100 500 2500]; % normalized
    rr = exp(linspace(log(1e-3),log(1.62),100)); % normalized work of adhesion i.e. angles from 180 to 90         
    
    hs = [];    
    F = {};    
    for j = 1:length(V)        
        r{j} = rr;
        k = nthroot(3*V(j)/pi + sqrt(1 + (3*V(j)/pi)^2),3);
        hs(j) = k - 1 / k;
        
        
        d0 = drop.segment(V(j),'radius',r{1}(1),'radius',1);       
        F{j} = par_fun(@(x) drop.segment_solve(V(j),'radius',x,'radius',1,'height',hs(j)).force,r{j})';        
    end
    if ~exist('output','dir')
        mkdir('output');
    end
    save('output/r_vs_snapin.mat','V','F','hs','r');