function plot_wa_vs_snapin

    close all;

    addpath('../sideaxes-m/');

    h1 = 1.75;

    fig_width = 3.25 * 2.54;
    fig_height = 3.58 * 3;

    left_margin = 0.75;
    bottom_margin = 0.75;

    figure('units','centimeter','position',[5 5 fig_width fig_height]);

    main_width = fig_width - left_margin;
    main_height = fig_height/2 - bottom_margin;

    load('output/wa_vs_snapin.mat');

    y = bottom_margin;
    axes('units','centimeter','position',[left_margin y main_width main_height]);

    set(gcf,'color','w');

    set(gca,'XDir','reverse');
    set(gca,'color',[0.975 0.975 0.975]);
    set(gca,'XColor','none')
    set(gca,'YColor','none')


    for j = 1:length(V)
        plot(angles,F{j},'b-');    
        ax = gca;
        ax.XDir = 'reverse';
        hold on;
        Fanal = arrayfun(@(x) analytical_angle(x,V(j),hs(j)),angles);
        plot(angles,Fanal,'b--');

        %xlim([1.5 hs(j)+0.02]);
        %ylim([-2 0.2]);
    end
    
    xlim([90 180]);    
    ylim([0 10]);    
    
    load('output/r_vs_snapin.mat');

    y = bottom_margin+fig_height/2;
    axes('units','centimeter','position',[left_margin y main_width main_height]);

    plot_r(r,V,F,hs);  
    xlim([0 1.5]);
    ylim([0 15]);
    
    axes('units','centimeter','position',[left_margin+0.7 y+1.4 3 3]);
    plot_r(r,V,F,hs);
    xlim([0 0.5]);
    ylim([0 0.5]);

    export_fig('output/Figure_Analytical_vs_Numerical.png','-r1200');
end

function plot_r(r,V,F,hs)
    set(gcf,'color','w');

    set(gca,'color',[0.975 0.975 0.975]);
    set(gca,'XColor','none')
    set(gca,'YColor','none')


    for j = 1:length(V)
        plot(r{j},F{j},'b-');                    
        hold on;
        Fanal = arrayfun(@(x) analytical_radius(x,V(j),hs(j)),r{j});
        plot(r{j},Fanal,'b--');      
    end
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