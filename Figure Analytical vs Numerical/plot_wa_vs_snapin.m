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
    
    export_fig('output/Figure_Analytical_vs_Numerical.png','-r1200','-nocrop','-painters');
end

function plot_inhomogeneous_case()
    load('output/r_vs_snapin.mat');
    ax = gca;

    xx = [0.24 0.265 0.29 0.302 0.33 0.35];
    Fanal = cell(1,length(V));
    Fanal2 = cell(1,length(V));
    for j = 1:length(V)
        ax.ColorOrderIndex = j;
        plot(r{j},F{j},'-','LineWidth',1);                    
        hold on;
        [Fanal{j},Fanal2{j}] = arrayfun(@(x) analytical_radius(x,V(j),hs(j)),r{j});
        ax.ColorOrderIndex = j;
        plot(r{j},Fanal{j},'--');              
        ax.ColorOrderIndex = j;
        plot(r{j},Fanal2{j},':');              
        x = xx(j);
        y1 = interp1(r{j},Fanal2{j},x);
        y2 = interp1(r{j},Fanal2{j},x+0.01);
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
    h_anal2 = plot(nan, nan, 'k:');
    [~,icons] = legend([h_anal,h_anal2,h_num],{'Analytical, linear (1)','Analytical, quadratic (2)','Numerical'},'FontName','Times New Roman','FontSize',7,'Location','NorthWest','box','off');       
    l = 0.08;
    icons(4).XData = icons(4).XData - [0 l];
    icons(6).XData = icons(6).XData - [0 l];
    icons(8).XData = icons(8).XData - [0 l];
    icons(1).Position = icons(1).Position - [l 0 0];
    icons(2).Position = icons(2).Position - [l 0 0];
    icons(3).Position = icons(3).Position - [l 0 0];
    
    sideaxes(ax,'south');
    ticks([0.2 0.3 0.4 0.6 0.8 1 1.5],-0.1,'Clipping','off','Color',[0.6 0.6 0.6]);
    labels([0.2 0.3 0.4 0.6 0.8 1 1.5],0.05,[],'FontName','Times New Roman','FontSize',8);    
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
        [Fanal,Fanal2] = arrayfun(@(x) analytical_angle(x,V(j),hs(j)),angles_anal);
        ax.ColorOrderIndex = j*2-1;
        plot(wa_anal,Fanal,'--');
        ax.ColorOrderIndex = j*2-1;
        hold on;
        plot(wa_anal,Fanal2,':');

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
    h_anal2 = plot(nan, nan, 'k:');
    [~,icons] = legend([h_anal,h_anal2,h_num],{'Analytical, linear (1)','Analytical, quadratic (2)','Numerical'},'FontName','Times New Roman','FontSize',7,'Location','NorthWest','box','off');       
    l = 0.08;
    icons(4).XData = icons(4).XData - [0 l];
    icons(6).XData = icons(6).XData - [0 l];
    icons(8).XData = icons(8).XData - [0 l];
    for k = 1:3
        icons(k).Position = icons(k).Position - [l 0 0];
    end
        
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

function [F,F2] = analytical_angle(angle,V,hs)
    [he,re,R] = equilibrium_distance_for_angle(angle,V);
    angle2 = 180 - asind(sind(angle)/re);
    f = (cotd(angle)+cscd(angle)*cosd(angle2))^2+1;    
    C = log(cotd(angle/2))+log(cotd(angle2/2))-(cosd(angle)+cosd(angle2))/(cosd(angle)*cosd(angle2)+f);
    k = 2*pi/C;    
    F = -k*(hs-he);
    
    D = log(cotd(angle/2))+log(cotd(angle2/2));
    d2F = -(6*(cosd(angle2)+1))*pi*(cosd(angle2)-1)*(16*cosd(angle)^5*(1/3)-6*cosd(angle)^7+7*cosd(angle)^9*(1/3)-10*cosd(angle2)^3*(1/3)-5*cosd(angle2)^3*cosd(angle)^2*(1/3)-cosd(angle2)^5-(1/3)*cosd(angle)^11+4*cosd(angle)^8*cosd(angle2)-30*cosd(angle)^7*cosd(angle2)^2+10*cosd(angle)^9*cosd(angle2)^2+3*cosd(angle)^7*cosd(angle2)^4+17*cosd(angle2)^2*cosd(angle)^5+64*cosd(angle2)^3*cosd(angle)^4+33*cosd(angle2)*cosd(angle)^4+31*cosd(angle2)^2*cosd(angle)^3-8*cosd(angle2)*cosd(angle)^2-12*cosd(angle2)^2*cosd(angle)-22*cosd(angle)^5*cosd(angle2)^4-3*cosd(angle)^4*cosd(angle2)^5+40*cosd(angle)^3*cosd(angle2)^4+10*cosd(angle)^2*cosd(angle2)^5-6*cosd(angle)*cosd(angle2)^4-22*cosd(angle2)*cosd(angle)^6-152*cosd(angle)^6*cosd(angle2)^3*(1/3)-cosd(angle)^11*cosd(angle2)^2-cosd(angle)^10*cosd(angle2)^3+cosd(angle)*cosd(angle2)^6+(1+cosd(angle2)^6+3*cosd(angle2)^4+3*cosd(angle2)^2-cosd(angle)^2+cosd(angle)^11*cosd(angle2)^3-cosd(angle)^2*cosd(angle2)^6+3*cosd(angle)^3*cosd(angle2)^3+9*cosd(angle)*cosd(angle2)-12*cosd(angle)^3*cosd(angle2)+24*cosd(angle)^2*cosd(angle2)^2+9*cosd(angle)*cosd(angle2)^5-10*cosd(angle)^9*cosd(angle2)^3-3*cosd(angle)^8*cosd(angle2)^4-3*cosd(angle2)^2*cosd(angle)^8+36*cosd(angle)^7*cosd(angle2)^3+21*cosd(angle)^6*cosd(angle2)^4+3*cosd(angle)^5*cosd(angle2)^5+21*cosd(angle2)^2*cosd(angle)^6-48*cosd(angle)^5*cosd(angle2)^3-45*cosd(angle)^4*cosd(angle2)^4-12*cosd(angle)^3*cosd(angle2)^5+3*cosd(angle)^5*cosd(angle2)-45*cosd(angle2)^2*cosd(angle)^4+24*cosd(angle)^2*cosd(angle2)^4+18*cosd(angle)*cosd(angle2)^3)*D-cosd(angle2)-cosd(angle)+38*cosd(angle)^8*cosd(angle2)^3*(1/3)+2*cosd(angle)^3*(1/3))/((1+cosd(angle))*sind(angle2)^2*(cosd(angle)-1)*((cosd(angle)^3*cosd(angle2)-3*cosd(angle)*cosd(angle2)-cosd(angle2)^2-1)*D-cosd(angle)^3-cosd(angle2)*cosd(angle)^2+cosd(angle)+cosd(angle2))^3);
    q = 1/2 * d2F / R;
    F2 = -k*(hs-he)-q*(hs-he).^2;    
end

function [F,F2] = analytical_radius(re,V,hs)
    [he,angle,R] = equilibrium_distance_for_radius(re,V);
    angle2 = 180 - asind(sind(angle)/re);
    f = 1;
    C = log(cotd(angle/2))+log(cotd(angle2/2))-(cosd(angle)+cosd(angle2))/(cosd(angle)*cosd(angle2)+f);
    k = 2*pi/C;    
    F = -k*(hs-he);
        
    D = log(cotd(angle/2))+log(cotd(angle2/2));
    X = cosd(angle)+cosd(angle2);
    Y = (cosd(angle)*cosd(angle2)+1);        
    q = pi*(3*D*Y^3-X^3-3*X*Y^2)/(D*Y-X)^3/R;   
    F2 = -k*(hs-he)-q*(hs-he).^2;    
end