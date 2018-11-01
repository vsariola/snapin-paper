function inset_radius(radius,zt,ut)
    ut_all = [ut ut(end) -ut(end) -fliplr(ut)];
    zt_all = [zt zt(end) zt(end) fliplr(zt)];
    
    fill(ut_all,zt_all,[0.9 0.9 1],'LineWidth',0.3);
    hold on;
    fill([-radius -radius radius radius],[-5 0 0 -5],[1 0.8 1],'LineWidth',0.3);
    fill([-ut(end) -ut(end) ut(end) ut(end)],[4 zt(end) zt(end) 4],[1 0.9 0.9],'LineWidth',0.3);       

    ax = gca;    
    aspect = ax.Position(4) / ax.Position(3);
    x2 = 1.5;
    x1 = -x2;    
    xlim([x1 x2]);        
    y1 = -0.4;
    y2 = y1+(x2-x1)*aspect;
    ylim([y1 y2]);    
    set(gca,'color',[0.975 0.975 0.975]);
    set(gca,'XColor','none')
    set(gca,'YColor','none')     
end