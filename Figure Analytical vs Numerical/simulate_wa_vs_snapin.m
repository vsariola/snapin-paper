function simulate_wa_vs_snapin

    addpath('../kenmotsu-drop-simulator/');
    addpath('../helpers');

    V = [5 100 2500]; % normalized    
    a1 = 179.5;
    a2 = [70 96 114];        
    
    hs = zeros(1,length(V));    
    F = cell(1,length(V));    
    angles = cell(1,length(V));
    for j = 1:length(V)
        angles{j} = acosd(linspace(cosd(a1)+1,cosd(a2(j))+1,100)-1);
        k = nthroot(3*V(j)/pi + sqrt(1 + (3*V(j)/pi)^2),3);
        hs(j) = k - 1 / k;                 
        F{j} = par_fun(@(x) drop.segment_solve(V(j),'angle',x,'radius',1,'height',hs(j)).force,angles{j})';        
    end
    if ~exist('output','dir')
        mkdir('output');
    end
    save('output/wa_vs_snapin.mat','V','angles','F','hs');