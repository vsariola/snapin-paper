function plot_case_I
fsize = 12;
    fig_width = 3.25 * 2.54;
    fig_height = fig_width * 2 / 3;    
    
    x1 = -0.17;
    x2 = 0.17;
    y1 = -0.27;
    y2 = 0.66;
    
    color_med = [0.85 0.85 0.85];
    color_dark = [0.6 0.6 0.6];
    addpath('../export_fig');
    addpath('../sideaxes-m');
    
    left_margin = 0.95;
    bottom_margin = 0.95;

    close all;
    figure('units','centimeter','position',[5 5 fig_width fig_height]);   
     
    main_width = fig_width - left_margin;
    main_height = fig_height - bottom_margin;
        
    ax = axes('units','centimeter','position',[left_margin bottom_margin main_width main_height]);               

    Dsim = load('output/simulated.mat');
    
    plot(Dsim.sim_dist,Dsim.sim_force,'b-');
   
    xlim([x1 x2]);
    ylim([y1 y2]);    
        
    set(gca,'XDir','reverse');
    set(gcf,'color','w');
    set(gca,'color',[0.975 0.975 0.975]);
    set(gca,'XColor','none')
    set(gca,'YColor','none')    
                        
    ax = gca;
    sideaxes(ax,'west');
    t = -0.2:0.1:0.6;
    ticks(t,-0.1,'Clipping','off','Color',color_dark);
    labels(t,0.05,[],'FontName','Times New Roman','FontSize',8);          
    labels(mean(ylim(ax)),0.5,'Force {\itF/(\it\gamma\itV^{-3})}','FontName','Times New Roman','FontSize',10,'orientation','vertical');
    
    sideaxes(ax,'south');
    t = -0.15:0.05:0.15;
    ticks(t,-0.1,'Clipping','off','Color',color_dark);
    labels(t,0.05,@(t)num2str(-t),'FontName','Times New Roman','FontSize',8);          
    labels(mean(xlim(ax)),0.4,'Distance {\Delta\ith/\itV^{-3}}','FontName','Times New Roman','FontSize',10);                    
   
    export_fig('output/case_I.png','-r1200','-nocrop','-painters');    
end