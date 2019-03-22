function plot_experimental_fd
    fsize = 12;
    fig_width = 3.25 * 2.54;
    fig_height = fig_width * 2 / 3;    
    
    x1 = -57;
    x2 = 5;
    y1 = -5.5;
    y2 = 0.5;
    
    insx1 = -2.5;
    insx2 = 1.2;
    insy1 = -0.22;
    insy2 = 0.12; 
    
    color_med = [0.85 0.85 0.85];
    color_dark = [0.6 0.6 0.6];
    addpath('../export_fig');
    addpath('../sideaxes-m');
    addpath('../kenmotsu-drop-simulator/');
    
    left_margin = 0.9;
    bottom_margin = 0.95;

    close all;
    figure('units','centimeter','position',[5 5 fig_width fig_height]);   
     
    main_width = fig_width - left_margin;
    main_height = fig_height - bottom_margin;
        
    ax = axes('units','centimeter','position',[left_margin bottom_margin main_width main_height]);    
    ax2 = axes('Position',[0.25 0.7 0.3 0.3]);
    axes(ax);
    
    set(ax,'units','normalized');   
    a = get(ax,'position');    
    b = get(ax2,'position');
    xx = ((b(1)+b(3)) - a(1)) / a(3);
    yy1 = (b(2) - a(2)) / a(4);
    yy2 = ((b(2)+b(4)) - a(2)) / a(4); 
    xc = (x2-x1)*xx+x1;
    y1c = (y2-y1)*yy1+y1;
    y2c = (y2-y1)*yy2+y1;    
    fill([insx1 insx2 insx2 insx1],[insy1 insy1 insy2 insy2],color_med,'EdgeColor','none');
    line(ax,[xc insx1],[y1c insy1],'LineStyle','-','Color',color_med);
    line(ax,[xc insx1],[y2c insy2],'LineStyle','-','Color',color_med);
    
    Dsim = load('output/simulated.mat');
    Dexp = load('output/experimental.mat');       
    
    he = 0; % The data has already been shifted so that he = 0    
            
    anal_dist = linspace(-60,10);
    [anal_force,anal_force2] = arrayfun(@(x) analytical_force(Dsim.pillarradius*1e-6,Dsim.a,Dsim.Vunnormalized,x*1e-6),anal_dist);
    anal_force = anal_force * 1e6 * Dsim.gamma;
    anal_force2 = anal_force2 * 1e6 * Dsim.gamma;
    
    plot_all;     
    xlim([x1 x2]);
    ylim([y1 y2]);    
        
    set(gcf,'color','w');
    set(gca,'color',[0.975 0.975 0.975]);
    set(gca,'XColor','none')
    set(gca,'YColor','none')    
                 
    s = [-20 -0.2];
    d = [5 0];
    arrow(s,s+d,'Length',5);
    text(s(1)+d(1)/2,s(2)+d(2)/2,'Approach','HorizontalAlign','center','VerticalAlign','top','FontName','Times New Roman','FontSize',8);
    
    s = [-7.2 -1];
    d = [-4.6455 -0.4744];      
    arrow(s,s+d,'Length',5);
    text(s(1)+d(1)/2,s(2)+d(2)/2,'Retract','HorizontalAlign','center','VerticalAlign','top','FontName','Times New Roman','FontSize',8,'Rotation',35);
        
    ax = gca;
    sideaxes(ax,'west');
    t = -5:1;
    ticks(t,-0.1,'Clipping','off','Color',color_dark);
    labels(t,0.05,[],'FontName','Times New Roman','FontSize',8);          
    labels(mean(ylim(ax)),0.5,'Force {\itF} (µN)','FontName','Times New Roman','FontSize',10,'orientation','vertical');
    
    sideaxes(ax,'south');
    t = -50:10:5;
    ticks(t,-0.1,'Clipping','off','Color',color_dark);
    labels(t,0.05,@(t)num2str(-t),'FontName','Times New Roman','FontSize',8);          
    labels(mean(xlim(ax)),0.4,'Distance {\Delta\ith} (µm)','FontName','Times New Roman','FontSize',10);        
        
    [~,icons] = legend([hsim,hanal2,hanal,hexp],{'Numerical','Analytical, quadratic (2)','Analytical, linear (1)','Experimental'},'FontName','Times New Roman','FontSize',8,'Location','SouthEast','box','off');       
    l = 0.08;        
    icons(1).Position = icons(1).Position - [l 0 0];
    icons(2).Position = icons(2).Position - [l 0 0];
    icons(3).Position = icons(3).Position - [l 0 0];
    icons(4).Position = icons(4).Position - [l 0 0];
    icons(5).XData = icons(5).XData - [0 l];
    icons(7).XData = icons(7).XData - [0 l];
    icons(9).XData = icons(9).XData - [0 l];
    icons(11).XData = icons(11).XData - [0 l];
    
    axes(ax2);
    plot_all;       
    xlim([insx1 insx2]);
    ylim([insy1 insy2]);    
    
    text(-1.2,0.05,'Snap-in','HorizontalAlign','center','VerticalAlign','bottom','FontName','Times New Roman','FontSize',8);
    s = [-1.2 0.04];
    e = [-1.4 -0.05];      
    arrow(s,e,'Length',5);
    
    set(gca,'XColor','none')
    set(gca,'YColor','none')    
    set(gca,'color','w');    
    
    sideaxes(ax2,'west');
    t = -0.2:0.1:0.1;
    ticks(t,-0.1,'Clipping','off','Color',color_dark);
    labels(t,0.05,[],'FontName','Times New Roman','FontSize',8);              
    
    sideaxes(ax2,'south');
    t = -2:1:1;
    ticks(t,-0.1,'Clipping','off','Color',color_dark);
    labels(t,0.05,@(t)num2str(-t),'FontName','Times New Roman','FontSize',8);                     
    
    export_fig('output/experimental_fd.png','-r1200','-nocrop','-painters');
    
    function plot_all
        %line([he he],[-10 10],'LineStyle','-','Color',color_dark);
        hold on;        

        hanal = plot(anal_dist,anal_force,'b--','LineWidth',0.5);
        hanal2 = plot(anal_dist,anal_force2,'b:','LineWidth',0.5);
        hsim = plot(Dsim.sim_dist,Dsim.sim_force,'b-','LineWidth',0.5);    
        hexp = plot(Dexp.exp_dist,Dexp.exp_force,'r-','LineWidth',0.5);
        hold off;  
    end
end

function [F,F2] = analytical_force(r1,r2,V,h)
    d = drop.segment(V,'radius',r1,'radius',r2);
    angle = d.angle1;
    angle2 = d.angle2;        
    C = log(cotd(angle/2))+log(cotd(angle2/2))-(cosd(angle)+cosd(angle2))/(cosd(angle)*cosd(angle2)+1);
    k = -2*pi/C;    
    F = k*h;
    
    D = log(cotd(angle/2))+log(cotd(angle2/2));
    X = cosd(angle)+cosd(angle2);
    Y = (cosd(angle)*cosd(angle2)+1);        
    q = pi*(3*D*Y^3-X^3-3*X*Y^2)/(D*Y-X)^3/d.R;   
    F2 = k*h+q*h.^2;    
end