function plot_experimental_fd
    fsize = 12;
    fig_width = 3.25 * 2.54;
    fig_height = fig_width * 2 / 3;    
    
    x1 = 0.68;
    x2 = 1.6;
    y1 = -1e-4*1e6;
    y2 = 260;
    
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
        
    D = load('output/simulated_data.mat');    
             
    hold on;
    set(gca,'XDir','reverse');
    
    xlim([x1 x2]);
    ylim([y1 y2]);    
        
    set(gcf,'color','w');
    set(gca,'color',[0.975 0.975 0.975]);
    set(gca,'XColor','none')
    set(gca,'YColor','none')    
    
    ind = find(D.data(:,2) <= 0 & D.data(:,1) < D.hs,1,'first');
    inds = [-1:1]+ind;
    p = polyfit(D.data(inds,1),D.data(inds,2),1);
    k_a = p(1);
    
    hE_a = D.data(ind,1);
    x = [hE_a-0.2e-3 hE_a+0.2e-3];
    plot(x*1000,-polyval(p,x)*1e6,'b--','LineWidth',0.5);
        
    theta_a = solve_angle_from_k(1,D.V/D.a^3,k_a / D.gamma);
    t = sprintf('k_1 \\approx %.2g N/m, \\theta_A \\approx %d°',k_a,round(theta_a));
    text(hE_a*1000,5e-6,t,'Rotation',19,'HorizontalAlign','center','VerticalAlign','top','FontName','Times New Roman','FontSize',8,'Interpreter','tex');
    
    ind = find(D.data(:,2) >= 0 & D.data(:,1) < hE_a,1,'first');
    inds = [-5:5]+ind;
    p = polyfit(D.data(inds,1),D.data(inds,2),1);
    k_r = p(1);    
    
    hE_r = D.data(ind,1);
    x = [hE_r-0.2e-3 hE_r+0.2e-3];
    plot(x*1000,-polyval(p,x)*1e6,'b--','LineWidth',0.5);    
    
    theta_r = solve_angle_from_k(1,D.V/D.a^3,k_r / D.gamma);
    t = sprintf('k_1 \\approx %.2g N/m, \\theta_R \\approx %d°',k_r,round(theta_r));
    text(hE_r*1000,5e-6,t,'Rotation',31,'HorizontalAlign','center','VerticalAlign','top','FontName','Times New Roman','FontSize',8,'Interpreter','tex');
        
    arrow([1.258,30],[1.258 5],'Length',5);
    text(1.258,35,'Snap-in','HorizontalAlign','center','VerticalAlign','bottom','FontName','Times New Roman','FontSize',8);
    
    text(1.05,100,'Contact line depins','HorizontalAlign','center','VerticalAlign','bottom','FontName','Times New Roman','FontSize',8);
    arrow([1 100],[0.91 35],'Length',5);
    
    arrow([1.56,35],[1.56 5],'Length',5);
    text(1.52,35,'Pull-off','HorizontalAlign','center','VerticalAlign','bottom','FontName','Times New Roman','FontSize',8);
    
    
    c = [0.84 170];
    s = [-0.015 10];
    arrow(c-s,c+s,'Length',5);
    text(c(1)+0.01,c(2)+10,'Approach','Rotation',45,'HorizontalAlign','center','VerticalAlign','bottom','FontName','Times New Roman','FontSize',8);
    
    c = [0.76 138];
    s = [0.011 -13];
    arrow(c-s,c+s,'Length',5);
    text(c(1)-0.01,c(2),'Retract','Rotation',65,'HorizontalAlign','center','VerticalAlign','top','FontName','Times New Roman','FontSize',8);
    %0.7582  136.9302

    h = plot(D.data(:,1)*1000,-D.data(:,2)*1e6,'b-','LineWidth',1);
    
%     %arrow(s,s+d,'Length',5);
%     %text(s(1)+d(1)/2,s(2)+d(2)/2,'Approach','HorizontalAlign','center','VerticalAlign','top','FontName','Times New Roman','FontSize',8);
%     
%     s = [-7.2 -1];
%     d = [-4.6455 -0.4744];      
%     arrow(s,s+d,'Length',5);
%     text(s(1)+d(1)/2,s(2)+d(2)/2,'Retract','HorizontalAlign','center','VerticalAlign','top','FontName','Times New Roman','FontSize',8,'Rotation',35);
%         
    ax = gca;
    sideaxes(ax,'west');
    t = -100:50:250;
    ticks(t,-0.1,'Clipping','off','Color',color_dark);
    labels(t,0.05,[],'FontName','Times New Roman','FontSize',8);          
    labels(mean(ylim(ax)),0.5,'Force {\itF} (µN)','FontName','Times New Roman','FontSize',10,'orientation','vertical');
    
    sideaxes(ax,'south');
    t = 0.6:0.1:1.6;
    ticks(t,-0.1,'Clipping','off','Color',color_dark);
    labels(t,0.05,[],'FontName','Times New Roman','FontSize',8);          
    labels(mean(xlim(ax)),0.4,'Distance {\ith} (mm)','FontName','Times New Roman','FontSize',10);        
        
    export_fig('output/simulated_experiment.png','-r1200','-nocrop','-painters');
    
    function ret = solve_angle_from_k(r2,Vol,k_per_gamma)                
        r1 = @(t1,t2) r2*sin(t1)/sin(t2); 
        h = @(t1,t2) -r2/tan(t2) - r1(t1,t2)/tan(t1);
        V = @(t1,t2) pi*h(t1,t2)/6*(3*r1(t1,t2)^2+3*r2^2+h(t1,t2)^2);        
        f = @(t1,t2) log(tan(t1/2))+log(tan(t2/2))+(cos(t1)+cos(t2))/(cos(t1)*cos(t2)+1+(cot(t1)+csc(t1)*cos(t2))^2)-2*pi/k_per_gamma;    
        C = @(x) [f(x(1),x(2)) V(x(1),x(2))-Vol];    
        options = optimoptions('fsolve','FunctionTolerance',1e-12);
        x = 180/pi*fsolve(C,[130 130]*pi/180,options);
        ret = x(1);
    end
end
