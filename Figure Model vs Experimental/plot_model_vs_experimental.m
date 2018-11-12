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

    F_anal = arrayfun(@(x) D.gamma * analytical_force(x,D.a,D.V,D.h),D.r);
    h_anal = plot(D.r*x_unit_scale,F_anal*y_unit_scale,'b--','LineWidth',0.5);    
    hold on;
   
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
         
    [~,icons] = legend([h_anal,h_num,h_exp],{'Analytical','Numerical','Experimental'},'FontName','Times New Roman','FontSize',7,'Location','NorthWest','box','off');       
    l = 0.08;    
    icons(1).Position = icons(1).Position - [l 0 0];
    icons(2).Position = icons(2).Position - [l 0 0];
    icons(3).Position = icons(3).Position - [l 0 0];
    icons(4).XData = icons(4).XData - [0 l];
    icons(6).XData = icons(6).XData - [0 l];
    icons(8).XData = icons(8).XData - [0 l];
    
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
    
    export_fig('output/Figure_Experimental_Comparison.pdf','-r1200','-nocrop','-painters');
end

function F = analytical_force(r1,r2,V,h)
    [he,angle] = equilibrium_distance_for_radius(r1,V,r2);
    angle2 = 180 - asind(sind(angle)*r2/r1);    
    C = log(cotd(angle/2))+log(cotd(angle2/2))-(cosd(angle)+cosd(angle2))/(cosd(angle)*cosd(angle2)+1);
    k = -2*pi/C;    
    F = k*(h-he);
end
    
    
    