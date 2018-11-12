function plot_wa_vs_snapin

    close all;

    addpath('../sideaxes-m/');
    addpath('../helpers/');    
    addpath('../export_fig/');    

    h1 = 1.75;

    fig_width = 3.25 * 2.54;
    fig_height = fig_width * 1.5;

    left_margin = 1.1;
    bottom_margin = 1;
    between_margin = 1.1;

    figure('units','centimeter','position',[5 5 fig_width fig_height]);

    main_width = fig_width - left_margin;
    main_height = (fig_height - bottom_margin - between_margin)/2;

    y = bottom_margin;
    axes('units','centimeter','position',[left_margin y main_width main_height]);    
    plot_inhomogeneous_case();
    
    y = bottom_margin+main_height+between_margin;
    ax = axes('units','centimeter','position',[left_margin y main_width main_height]);
    plot_homogeneous_case();
    
    export_fig('output/Figure_Analytical_vs_Numerical.pdf','-r1200','-nocrop','-painters');
end

function plot_inhomogeneous_case()
    load('output/r_vs_snapin.mat');
    ax = gca;

    xx = [0.24 0.265 0.29 0.302 0.33 0.35];
    Fanal = cell(1,length(V));
    for j = 1:length(V)
        ax.ColorOrderIndex = j;
        plot(r{j},F{j},'-','LineWidth',1);                    
        hold on;
        Fanal{j} = arrayfun(@(x) analytical_radius(x,V(j),hs(j)),r{j});
        ax.ColorOrderIndex = j;
        plot(r{j},Fanal{j},'--');              
        x = xx(j);
        y1 = interp1(r{j},Fanal{j},x);
        y2 = interp1(r{j},Fanal{j},x+0.01);
        rot = atan2d((log(y2)-log(y1))*0.16,log(x+0.01) - log(x));        
        if j == 1
            t = sprintf('{\\itV}/{\\itr}_2^3 = %g',V(j));
            s = 1.1;
        else
            t = sprintf('= %g',V(j));
            s = 0.96;
        end
        text(x,y1*s,t,'FontName','Times New Roman','FontSize',7,'Rotation',rot,'HorizontalAlign','right','VerticalAlign','bottom');
    end
    set(gcf,'color','w');
    set(gca,'color',[0.975 0.975 0.975]);
    set(gca,'XColor','none')
    set(gca,'YColor','none')
    set(gca,'YScale','log');
    set(gca,'XScale','log');  
    xlim([1.6e-1 max(r{1})]);
    ylim([1.1e-3 20]);
            
    %x = 0.9;
    %y = interp1(r{1},Fanal{1},x);
    %text(x,y+0.3,'Analytical','FontName','Times New Roman','FontSize',7,'Rotation',32,'HorizontalAlign','center','VerticalAlign','bottom');
    %x = 1.3;
    %y = interp1(r{1},F{1},x);
    %text(x,y,'Numerical','FontName','Times New Roman','FontSize',7,'Rotation',6,'HorizontalAlign','center','VerticalAlign','bottom');        
    h_num = plot(nan, nan, 'k-','LineWidth',1);
    h_anal = plot(nan, nan, 'k--');
    [~,icons] = legend([h_anal,h_num],{'Analytical','Numerical'},'FontName','Times New Roman','FontSize',7,'Location','NorthWest','box','off');       
    l = 0.08;
    icons(3).XData = icons(3).XData - [0 l];
    icons(5).XData = icons(5).XData - [0 l];
    icons(1).Position = icons(1).Position - [l 0 0];
    icons(2).Position = icons(2).Position - [l 0 0];
    
    sideaxes(ax,'south');
    ticks([0.2 0.3 0.45 0.7 1 1.5],-0.1,'Clipping','off','Color',[0.6 0.6 0.6]);
    labels([0.2 0.3 0.45 0.7 1 1.5],0.05,[],'FontName','Times New Roman','FontSize',8);    
    labels(exp(mean(log(xlim(ax)))),0.4,'Radius {\itr}_1/{\itr}_2','FontName','Times New Roman','FontSize',10);
    
    sideaxes(ax,'west');
    ticks([1e-2 1e-1 1 10] ,-0.1,'Clipping','off','Color',[0.6 0.6 0.6]);
    labels([1e-2 1e-1 1 10],0.05,[],'FontName','Times New Roman','FontSize',8);          
    labels(exp(mean(log(ylim(ax)))),0.5,'Snap-in Force {\itF}/({\it\gamma r}_2)','FontName','Times New Roman','FontSize',10,'orientation','vertical');
    
    labels(max(xlim)-2,max(ylim),'b','location','east','FontName','Times New Roman','FontSize',10,'FontWeight','bold');      
end

function plot_homogeneous_case()
    load('output/wa_vs_snapin.mat');
   
    ax = gca;

    set(gcf,'color','w');

    set(gca,'XDir','reverse');
    set(gca,'color',[0.975 0.975 0.975]);
    set(gca,'XColor','none')
    set(gca,'YColor','none')    

    xx = [3.8e-3 3.55e-3 2.3e-3];
    s = [1.25 1.35 1.7];
    for j = 1:length(V)
        ax.ColorOrderIndex = j*2-1;
        wa = 1+cosd(angles{j});
        plot(wa,F{j},'-','LineWidth',1);    
        ax = gca;        
        hold on;
        angles_anal = linspace(angles{j}(1),acosd(1.2-1));
        wa_anal = 1+cosd(angles_anal);
        Fanal = arrayfun(@(x) analytical_angle(x,V(j),hs(j)),angles_anal);
        ax.ColorOrderIndex = j*2-1;
        plot(wa_anal,Fanal,'--');

        x = xx(j);
        y1 = interp1(wa_anal,Fanal,x);
        y2 = interp1(wa_anal,Fanal,x+0.01);
        rot = atan2d((log(y2)-log(y1))*0.6,log(x+0.01) - log(x));        
        if j == 3
            t = sprintf('{\\itV}/{\\itr}_2^3 = %g',V(j));
        else
            t = sprintf('= %g',V(j));            
        end
        text(x,y1*s(j),t,'FontName','Times New Roman','FontSize',7,'Rotation',rot,'HorizontalAlign','center','VerticalAlign','middle');
    end
    
    set(gcf,'color','w');
    set(gca,'color',[0.975 0.975 0.975]);
    set(gca,'XColor','none')
    set(gca,'YColor','none')
    set(gca,'YScale','log');
    set(gca,'XScale','log');  
    xlim([1e-3 1.2]);    
    
    ylim([1e-3 12]);    
    
    
    h_num = plot(nan, nan, 'k-','LineWidth',1);
    h_anal = plot(nan, nan, 'k--');
    [~,icons] = legend([h_anal,h_num],{'Analytical','Numerical'},'FontName','Times New Roman','FontSize',7,'Location','NorthWest','box','off');       
    l = 0.08;
    icons(3).XData = icons(3).XData - [0 l];
    icons(5).XData = icons(5).XData - [0 l];
    icons(1).Position = icons(1).Position - [l 0 0];
    icons(2).Position = icons(2).Position - [l 0 0];
    
    
    sideaxes(ax,'west');
    ticks([1e-2 1e-1 1 10],-0.1,'Clipping','off','Color',[0.6 0.6 0.6]);
    labels([1e-2 1e-1 1 10],0.05,[],'FontName','Times New Roman','FontSize',8);          
    labels(exp(mean(log(ylim(ax)))),0.5,'Snap-in Force {\itF}/({\it\gamma r}_2)','FontName','Times New Roman','FontSize',10,'orientation','vertical');
    labels(max(xlim)-2,max(ylim),'a','location','east','FontName','Times New Roman','FontSize',10,'FontWeight','bold');   
    
    sideaxes(ax,'south');
    t = [3e-3 1e-2 3e-2 1e-1 3e-1 1];
    ticks(t,-0.1,'Clipping','off','Color',[0.6 0.6 0.6]);
    labels(t,0.05,[],'FontName','Times New Roman','FontSize',8);    
    labels(t,0.5,@(x) sprintf('%.0f',acosd(x-1)),'FontName','Times New Roman','FontSize',8);    
    
    labels(0.77e-3,0,'1+cos(\theta_1)','FontName','Times New Roman','FontSize',10,'horizontalalign','center','clipping','off');
    labels(0.77e-3,0.45,'\theta_1 \approx','FontName','Times New Roman','FontSize',10,'horizontalalign','center','clipping','off');       
end

function F = analytical_angle(angle,V,hs)
    [he,re] = equilibrium_distance_for_angle(angle,V);
    angle2 = 180 - asind(sind(angle)/re);
    f = (cotd(angle)+cscd(angle)*cosd(angle2))^2+1;    
    C = log(cotd(angle/2))+log(cotd(angle2/2))-(cosd(angle)+cosd(angle2))/(cosd(angle)*cosd(angle2)+f);
    k = 2*pi/C;    
    F = -k*(hs-he);
end

function F = analytical_radius(re,V,hs)
    [he,angle] = equilibrium_distance_for_radius(re,V);
    angle2 = 180 - asind(sind(angle)/re);
    f = 1;
    C = log(cotd(angle/2))+log(cotd(angle2/2))-(cosd(angle)+cosd(angle2))/(cosd(angle)*cosd(angle2)+f);
    k = 2*pi/C;    
    F = -k*(hs-he);
end