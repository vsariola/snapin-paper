function plot_schematic
    addpath('../kenmotsu-drop-simulator/');        
    addpath('../helpers/');        
    addpath('../sideaxes-m/');        
    close all;

    fig_width = 3.25 * 2.54;
    fig_height = fig_width * 2 / 3;           
    
    figure('units','centimeter','position',[5 5 fig_width fig_height]);    
    axes('units','centimeter','position',[0 0 fig_width fig_height]);    
    
    angle1 = 150;
    radius2 = 1;
    V = 5;
    h = 1.5;
    d = drop.create_ar(h,V,angle1,radius2);
    radius = d.radius1;

    ut = d.r;
    zt = d.z;
    
    ut_all = [ut ut(end) -ut(end) -fliplr(ut)];
    zt_all = [zt zt(end) zt(end) fliplr(zt)];
    

    hold on;        
        
    fill(ut_all,zt_all,[0.9 0.9 1]);   
    
    text(ut(end)/2,zt(end)+0.1,'{\itr}_2','FontSize',10,'FontName','Times New Roman','horizontalalign','center','verticalalign','middle');
    text(ut(1)/2,zt(1)-0.1,'{\itr}_1','FontSize',10,'FontName','Times New Roman','horizontalalign','center','verticalalign','middle');
    
    x = ut(end);
    y = zt(end);
    c = 0.4;
    dx = -cosd(d.angle2)*c;
    dy = -sind(d.angle2)*c;    
    line([x x+dx],[y y+dy],'LineStyle','-','Color','k');    
           
    c = 0.15;
    a = linspace(0,d.angle2,20);
    plot(x-cosd(a)*c,y-sind(a)*c,'LineStyle','-','Color','k'); 
    labels(x-0.2,y-0.2,'\theta_2','FontSize',10,'FontName','Times New Roman');
    
    x = ut(1);
    y = zt(1);
    c = 0.4;
    dx = -cosd(d.angle1)*c;
    dy = sind(d.angle1)*c;    
    line([x x+dx],[y y+dy],'LineStyle','-','Color','k');    
    
    c = 0.15;
    a = linspace(0,d.angle1,20);
    plot(x-cosd(a)*c,y+sind(a)*c,'LineStyle','-','Color','k'); 
    labels(x-0.2,y+0.2,'\theta_1','FontSize',10,'FontName','Times New Roman');
    
    line([ut(end) 1.3],[zt(end) zt(end)],'LineStyle','-','Color',[0.6 0.6 0.6]);
    line([ut(1) 1.3],[zt(1) zt(1)],'LineStyle','-','Color',[0.6 0.6 0.6]);
    line([1.25 1.25],[zt(end) zt(1)],'LineStyle','-','Color',[0.6 0.6 0.6]);
    

    ax = gca;    
    aspect = ax.Position(4) / ax.Position(3);
    x2 = 1.5;
    x1 = -x2;    
    xlim([x1 x2]);        
    y1 = -0.2;
    y2 = y1+(x2-x1)*aspect;
    ylim([y1 y2]);    
    line([0 0],ylim,'LineStyle','-.','Color','k')
    set(gca,'XColor','none')
    set(gca,'YColor','none')     
end