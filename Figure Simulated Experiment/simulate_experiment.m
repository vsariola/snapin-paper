addpath('../Kenmotsu-drop-simulator')

theta_a = 150;
theta_r = 120;
a = 0.5e-3; % Radius of the probe disk, in m (SI-units)
gamma = 72e-3; % Surface tension, in N/m (SI-units)
V = 1.53e-9;
Vt = V/a^3; % normalized volume
k = nthroot(3*Vt/pi + sqrt(1 + (3*Vt/pi)^2),3);
hs = (k - 1 / k) * a; % snapin distance
step = 0.005e-3;
maxh = 1.7e-3;
minh = 0.7e-3;
pulloff_h = 1.56e-3;

angle0 = d.cap(Vt,'radius',1).angle1;

data = [];
for h = maxh:-step:hs
    data = [data;h 0];
end
d = drop.segment(Vt,'angle',theta_a,'radius',1);
for h = hs:-step:minh
    d = d.solve('volume',Vt,'angle1',theta_a,'radius2',1,'height',h/a);
    data = [data;h gamma * a * d.force];
    d.show(); drawnow;
end
r_pinned = d.radius1 * a;
h = minh;
while true
    d = d.solve('volume',Vt,'radius1',r_pinned/a,'radius2',1,'height',h/a);
    if d.angle1 < theta_r
        break
    end
    data = [data;h gamma * a * d.force];
    h = h + step;
    d.show(); drawnow;
end
for h = h:step:pulloff_h
    d = d.solve('volume',Vt,'angle1',theta_r,'radius2',1,'height',h/a);
    data = [data;h gamma * a * d.force];
end
for h = pulloff_h:step:maxh
    data = [data;h 0];
end
%%
if ~exist('output','dir')
    mkdir('output');
end
save('output/simulated_data.mat','gamma','a','V','data','theta_a','theta_r','hs','angle0');