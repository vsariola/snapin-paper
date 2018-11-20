function simulate_fd
    addpath('../helpers/');
    addpath('../kenmotsu-drop-simulator/');    

    dmin = -57; % um
    dmax = 5; % um
    forcescale = 1e-6;
    rscale = 1e-6;
    pillarradius = 10;

    a = 0.5e-3;
    gamma = 72e-3;
    sg = 0;
    F = 15e-6;
    g = 9.81;
    m = F/g;
    rho = 1000;
    Vunnormalized = m/rho;
    V = Vunnormalized / a^3;
    sim_dist = linspace(dmin,dmax);    
    de = drop.segment(V,'radius',pillarradius * rscale / a,'radius',1);
    he = de.height;
    
    sim_force = par_fun(@wrapper,sim_dist);
    
    if ~exist('output','dir')
        mkdir('output');
    end
    save('output/simulated.mat','sim_dist','sim_force','gamma','pillarradius','Vunnormalized','a');   
    
    
    function f = wrapper(x)
        drel = x * rscale / a;
        dabs = he-drel;
        d = drop.segment_solve(V,'radius',pillarradius * rscale / a,'radius',1,'height',dabs);        
        f = gamma * a * d.force - V*sg;     
        f = -f / forcescale;
    end
end