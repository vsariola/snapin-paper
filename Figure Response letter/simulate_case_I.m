function simulate_fd
    addpath('../helpers/');
    addpath('../kenmotsu-drop-simulator/');    

    dmin = -0.17;
    dmax = 0.17;
    
    theta1 = 150;
    theta2 = 120;
    V = 1;           
    sim_dist = linspace(dmin,dmax);    
    de = drop.segment(V,'angle',theta1,'angle',theta2);
    he = de.height;
    
    sim_force = par_fun(@wrapper,sim_dist);
    
    if ~exist('output','dir')
        mkdir('output');
    end
    save('output/simulated.mat','sim_dist','sim_force','V','theta1','theta2','he');   
    
    
    function f = wrapper(x)
        dabs = x+he;
        d = drop.segment_solve(V,'angle',theta1,'angle',theta2,'height',dabs);        
        f = -d.force;             
    end
end