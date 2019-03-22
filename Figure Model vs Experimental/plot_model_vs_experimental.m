function plot_model_vs_experimental
    addpath('../kenmotsu-drop-simulator/');        
    addpath('../helpers/');        
    addpath('../sideaxes-m/');        
    close all;

    fig_width = 3.25 * 2.54;
    fig_height = fig_width * 2 / 3;

    left_margin = 0.9;
    bottom_margin = 0.95;
    
    x_unit_scale = 1e6;
    y_unit_scale = 1e6;    
    
    figure('units','centimeter','position',[5 5 fig_width fig_height]);
    
    main_width = fig_width - left_margin;
    main_height = fig_height - bottom_margin;

    D = load('output/r_vs_snapin.mat');
    
    axes('units','centimeter','position',[left_margin bottom_margin main_width main_height]);    

    [F_anal,F_anal2] = arrayfun(@(x)analytical_force(x,D.a,D.V,D.h),D.r);    
    F_anal = F_anal * D.gamma;
    F_anal2 = F_anal2 * D.gamma;
    h_anal = plot(D.r*x_unit_scale,F_anal*y_unit_scale,'b--','LineWidth',0.5);    
    hold on;
    h_anal2 = plot(D.r*x_unit_scale,F_anal2*y_unit_scale,'b:','LineWidth',0.5);       
   
    h_num = plot(D.r*x_unit_scale,D.F*y_unit_scale,'b-','LineWidth',1);    
    xlim([7.5 800]);
    ylim([min(D.F)/1.2 max(F_anal)]*y_unit_scale);
    
    D = load('data_Liimatainen_et_al_2018.mat');
    r = D.data_pillar_radius*1e-6*x_unit_scale;
    F = D.data_pillar_snapin*1e-6*y_unit_scale;
    %F_low = D.data_pillar_snapin_low*1e-6*y_unit_scale;
    %F_high = D.data_pillar_snapin_high*1e-6*y_unit_scale;
    %line([r;r],[F_low;F_high],'Color','r','LineWidth',0.5); % the error
    %bars are so small that they cannot be even seen, so no point plotting
    hold on;
    plot(r,F,'.','MarkerSize',12,'Color',[0.975 0.975 0.975]);
    h_exp = plot(r,F,'r.','MarkerSize',8);
    set(gcf,'color','w');
    set(gca,'color',[0.975 0.975 0.975]);
    set(gca,'XColor','none')
    set(gca,'YColor','none')
    set(gca,'YScale','log');
    set(gca,'XScale','log');  
         
    [~,icons] = legend([h_anal,h_anal2,h_num,h_exp],{'Analytical, linear (1)','Analytical, quadratic (2)','Numerical','Experimental'},'FontName','Times New Roman','FontSize',7,'Location','NorthWest','box','off');       
    l = 0.08;    
    icons(1).Position = icons(1).Position - [l 0 0];
    icons(2).Position = icons(2).Position - [l 0 0];
    icons(3).Position = icons(3).Position - [l 0 0];
    icons(4).Position = icons(4).Position - [l 0 0];
    icons(5).XData = icons(5).XData - [0 l];
    icons(7).XData = icons(7).XData - [0 l];
    icons(9).XData = icons(9).XData - [0 l];
    icons(11).XData = icons(11).XData - [0 l];
    
    ax = gca;
    
    sideaxes(ax,'west');
    t = [0.1 0.3 1 3 10 30 100 300];
    ticks(t,-0.1,'Clipping','off','Color',[0.6 0.6 0.6]);
    labels(t,0.05,[],'FontName','Times New Roman','FontSize',8);          
    labels(exp(mean(log(ylim(ax)))),0.5,'Snap-in Force {\itF} (µN)','FontName','Times New Roman','FontSize',10,'orientation','vertical');
    
    sideaxes(ax,'south');
    t = [10 20 50 100 200 400 1000];
    ticks(t,-0.1,'Clipping','off','Color',[0.6 0.6 0.6]);
    labels(t,0.05,[],'FontName','Times New Roman','FontSize',8);          
    labels(exp(mean(log(xlim(ax)))),0.4,'Pillar radius {\itr}_1 (µm)','FontName','Times New Roman','FontSize',10);
    
    export_fig('output/Figure_Experimental_Comparison.png','-r1200','-nocrop','-painters');
end

function [F,F2] = analytical_force(r1,r2,V,h)
    [he,angle,R] = equilibrium_distance_for_radius(r1,V,r2);
    angle2 = 180 - asind(sind(angle)*r2/r1);    
    D = log(cotd(angle/2))+log(cotd(angle2/2));
    X = cosd(angle)+cosd(angle2);
    Y = (cosd(angle)*cosd(angle2)+1);
    C = D-X/Y;
    k = -2*pi/C;        
    q = -pi*(3*D*Y^3-X^3-3*X*Y^2)./(D*Y-X)^3/R;
    F = k*(h-he);
    F2 = k*(h-he)+q*(h-he).^2;
end
    
    
    