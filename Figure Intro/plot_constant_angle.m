close all;

addpath('../sideaxes-m/');

h1 = 1.75;

fig_width = 3.25 * 2.54;
fig_height = 2 * 2.54;

inset_height = 1.6;
left_margin = 0.75;
bottom_margin = 0.75;
space = 0.1;

figure('units','centimeter','position',[5 5 fig_width fig_height]);

main_width = fig_width - left_margin;
main_height = fig_height - space - inset_height - bottom_margin;
axes('units','centimeter','position',[left_margin bottom_margin main_width main_height]);

load('output/constant_angle.mat');
set(gcf,'color','w');

line([0 2],[0 0],'LineStyle','-','Color',[0.6 0.6 0.6]);
line([he he],[-4 4],'LineStyle','-','Color',[0.6 0.6 0.6]);
line([hs hs],[-4 4],'LineStyle','-','Color',[0.6 0.6 0.6]);
hold on;
plot([2 h(1) h],[0 0 -F],'b-','LineWidth',1);
set(gca,'XDir','reverse');
set(gca,'color',[0.975 0.975 0.975]);
set(gca,'XColor','none')
set(gca,'YColor','none')

angle2 = 180 - asind(sind(angle)/re);
f = (cotd(angle)+cscd(angle)*cosd(angle2))^2+1;
C = log(cotd(angle/2))+log(cotd(angle2/2))-(cosd(angle)+cosd(angle2))/(cosd(angle)*cosd(angle2)+f);
k = 2*pi/C;
x = linspace(he-0.05,hs+0.1);
plot(x,k*(x-he),'b--');

plot([h1 hs hs hs he],[0 0 -F(1) k*(hs-he) 0],'r.','MarkerSize',10);

xlim([he-0.05 hs+0.1]);
ylim([-0.7 0.2]);


inset_y = bottom_margin+main_height+space;
inset_width = (main_width - space*3)/4;
axes('units','centimeter','position',[left_margin inset_y inset_width inset_height]);
[zt,ut] = spherical_cap(V,h1);
inset_angle(zt,ut);
inset_x2 = left_margin+inset_width+space;
axes('units','centimeter','position',[inset_x2 inset_y inset_width inset_height]);
[zt,ut] = spherical_cap(V,hs);
inset_angle(zt,ut);
inset_x3 = inset_x2+inset_width+space;
axes('units','centimeter','position',[inset_x3 inset_y inset_width inset_height]);
inset_angle(zts,uts);
inset_x4 = inset_x3+inset_width+space;
axes('units','centimeter','position',[inset_x4 inset_y inset_width inset_height]);
inset_angle(zte,ute);
export_fig('output/Constant_Angle.pdf','-nocrop','-painters');
