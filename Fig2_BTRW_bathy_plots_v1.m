%Fig2_BTRW_bathy_plots_v1.m
% Plot maps of bathymetry for Figure 2 of Farrar et al. 2021 BTRW paper.  Intent is to show both 
%'full' bathymetry and the smoothed bathymetry used for the ray tracing calculation.
% (If this doesn't work, please check https://github.com/jtomfarrar)
% 
% This code was run on Matlab Version 9.6.0.1072779 (R2019a).
%
% Dependencies:
% H_L06_23_Ted_ray_tracing_bathy_final.mat (bathymetry used for ray tracing in paper)
% coast_tom.mat (coastline file of unknown origin)
% ETOPO2v2c.mat (bathymetry data from https://www.ngdc.noaa.gov/mgg/global/etopo2.html, saved in Matlab format)
% export_fig.m (https://www.mathworks.com/matlabcentral/fileexchange/23629-export_fig)
% ecolorbar.m by Jan Even Nilsen (a version is bundled with this code-- also see https://github.com/evenrev1)
%
% Started 1/24/2021.
% copyright (c) 2021 Tom Farrar, jfarrar@whoi.edu, 
% MIT License
% https://github.com/jtomfarrar


%Used for Feb 2021 revised submission to JPO
fig_output_dir=['C:\Users\jfarrar\Documents\MATLAB\figs\BTRW_long_range\Dec2020\'];
output_plots=1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Starting from model-obs plots in ModelSSHPlot_tom_v5.m
load H_L06_23_Ted_ray_tracing_bathy_final.mat
load coast_tom;%coastlat coastlon;

figbox=[1 1 6.5 5.5];%Modified to give room for colorbar label at bottom
figure
cax=[-6000 0];
V=cax(1):diff(cax)./24:cax(2);
contourf(xH,yH,-H,V,'LineStyle', 'none');
hold on
contour(xH,yH,-H,[-4000 -3000],'k');
axis image
title(['Ray tracing bathymetry'])
  hold on
  mapshow(coastlon,coastlat, 'displaytype','polygon','facecolor',0.65*[1 1 1])
  mapshow(coastlon+360,coastlat, 'displaytype','polygon','facecolor',0.65*[1 1 1])
  set(gca,'tickdir','out')
  set(gca,'xminortick','on')
  set(gca,'yminortick','on')
  set(gca,'layer','top')
ax=[124.0000  291.0000  -28.8750   60.6250]%full obs domain

axis(ax) 
[h0,h1]=ecolorbar(V,'b','Seafloor topgraphy (m)',1/30,[]);

if output_plots==1
fig = gcf; 
set(gcf, 'Color', 'w');
u = fig.Units;
fig.Units = 'inches';
fig.OuterPosition=figbox;
hold on
set(h1,'Outerposition',[0 +0.075 1 0.1118])
mapax_x;mapax_y;
export_fig([fig_output_dir 'ray_tracing_topo' ], '-pdf', '-p.01');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  load ETOPO2v2c.mat;
  ffn=find(x<0);
  ffp=find(x>=0);
  lontopo=[x(ffp); x(ffn)+360];
  topo=z([ffp ffn],:)';
  ff=find(lontopo>100&lontopo<290);topo=topo(:,ff);lontopo=lontopo(ff);
  lattopo=y;


figbox=[1 1 6.5 5.5];%Modified to give room for colorbar label at bottom
figure
cax=[-6000 0];
V=cax(1):diff(cax)./24:cax(2);
contourf(lontopo,lattopo,topo,V,'LineStyle', 'none');
hold on
contour(lontopo,lattopo,topo,[-4000 -3000],'k');
axis image
title(['ETOPO2 bathymetry'])
  hold on
  mapshow(coastlon,coastlat, 'displaytype','polygon','facecolor',0.65*[1 1 1])
  mapshow(coastlon+360,coastlat, 'displaytype','polygon','facecolor',0.65*[1 1 1])
  set(gca,'tickdir','out')
  set(gca,'xminortick','on')
  set(gca,'yminortick','on')
  set(gca,'layer','top')
ax=[124.0000  291.0000  -28.8750   60.6250]%full obs domain
axis(ax) 
[h0,h1]=ecolorbar(V,'b','Seafloor topgraphy (m)',1/30,[]);

if output_plots==1
fig = gcf; 
set(gcf, 'Color', 'w');
u = fig.Units;
fig.Units = 'inches';
fig.OuterPosition=figbox;
hold on
set(h1,'Outerposition',[0 +0.075 1 0.1118])
mapax_x;mapax_y;
export_fig([fig_output_dir 'full_topo' ], '-pdf', '-p.01');
end









