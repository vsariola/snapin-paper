function process_experimental_fd(debias)

    if (nargin < 1)
        debias = 0;
    end

    speed_slow = 2; % um/s
    speed_fast = 10; % um/s
    sampling_rate = 50; % Hz
    
    % There's an evaporation going on, so the force curve slowly drifts up
    % Assuming it's linear, we use the initial part of the data to reset 
    % the base line.
    dtrend_start = 23;
    dtrend_end = 339;  
    
    % These points were picked from the data by hand        
    s1 = 375; % At time s1, the speed changes from +slow to +fast
    s2 = 396; % At time s2, the speed changes from +fast to -fast (retract)
           
    force = load('10µm_radius_pillar1_measurement1.csv')';   
    
    % Detrending
    i = dtrend_start:dtrend_end;
    p = polyfit(i,force(i),1);
    force = force - polyval(p,1:length(force));

    
    d1 = (1:s1) / sampling_rate * speed_slow;
    d2 = (1:(s2-s1)) / sampling_rate * speed_fast;
    d3 = (1:(length(force)-s2)) / sampling_rate * -speed_fast;    
    d = [d1 (d2+d1(end))];
    d = [d (d3+d(end))];
    
    % The equilibrium distance F = 0 is located somewhere in the range
    % s1:s2
    % Fit a line p(1)*x+p(2)=F and solve x from F=0
    i = s1:s2;
    p = polyfit(d(i),force(i),1);
    d0 = -p(2)/p(1);    
    
    exp_dist = (d - d0)';
    exp_force = force;
    
    if ~exist('output','dir')
        mkdir('output');
    end
    save('output/experimental.mat','exp_dist','exp_force'); 
end