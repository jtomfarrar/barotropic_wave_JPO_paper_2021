%Matlab code to make model-obs comparison plot (Fig 9) for Farrar et al. JPO paper (revised Feb 2021)
% study of barotropic waves that radiate away from 
% tropical instability waves in the Pacific.
% This script should reproduce Fig 9. It must be run after 
% Fig3to8_xspectra_aviso_basin14_final.m
% (If this doesn't work, please check https://github.com/jtomfarrar)
%
% This code was run on Matlab Version 9.6.0.1072779 (R2019a).
% 
% Dependencies:
% All dependencies of Fig3to8_xspectra_aviso_basin14_final.m
% twilight_colormap2.mat (taken from the matplotlib twilight colormap)
% ecolorbar.m by Jan Even Nilsen (a version is bundled with this code-- also see https://github.com/evenrev1)
%
% ModelSSHPlot_tom_v2.m is the version Tom was using in BTRW JPO paper (Oct 2019-March 2020)
% ModelSSHPlot_tom_v4.m is for the revision, using Ted's new model output (Nov/Dec 2020) (v3 was intermediate)
% ModelSSHPlot_tom_v5.m narrows gain color scale (making it more saturated)
%
% run after xspectra_aviso_basin14.m

%First replot Obs fields, with minor modification from prior obs plots:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot map of coherence phase
figbox=[1 1 6.5 5.5];%Modified to give room for colorbar label at bottom
coher_phi2=coher_phi;
% Apply phase offset to model (model phase is arbitrary).
phi=coher_phi2;
phi_off=0;%-130;
phi=phi+phi_off;
ff=find(phi>180);
phi(ff)=phi(ff)-360;
ff=find(phi<=-180);
phi(ff)=phi(ff)+360;
coher_phi2=phi;

fflatblank=find(lat<=10);
coher_phi2(:,fflatblank)=NaN;
ax=[163  270 4.95   66];%for model domain
figure
plot(lon(1),lat(1))
hold on
cax=[-180 180];
V=cax(1):diff(cax)./24:cax(2);
imagesc(lon,lat,(coher_phi2)');shading interp
axis image
if do_titles==1
title(['Phase versus ' num2str(round(lat(fflat))) 'N, ' num2str(round(360-lon(fflon))) 'W' ' for ' num2str(round(1./w0(ffw(length(ffw))))) '-' num2str(round(1./w0(ffw(1)))) ' days'])
end
load twilight_colormap2;colormap(twilight2)
cax=caxis;
[h0,h1]=ecolorbar(V,'b','Phase (degrees)',1/30,[]);
scatter(lon(fflon),lat(fflat),'ko','filled')
  hold on
  mapshow(coastlon,coastlat, 'displaytype','polygon','facecolor',.65*[1 1 1])
  mapshow(coastlon+360,coastlat, 'displaytype','polygon','facecolor',.65*[1 1 1])
  set(gca,'tickdir','out')
  set(gca,'xminortick','on')
  set(gca,'yminortick','on')
  set(gca,'layer','top')
caxis(cax)
  axis([ax+buff.*[-1 1 -1 1]])
if 1==2;%exist('k_est')%plot BTRW initial ray paths
[hline]=plotBTRW_phase_line(crest_slope,cgx,cgy,x1,y1+5,length_x);
%[hline]=plotBTRW_phase_line(crest_slope,cgx,cgy,x1+1./k_est,y1+5,length_x);
end
if output_plots==1
fig = gcf; %figure;
set(gcf, 'Color', 'w');
fig.Units = 'inches';
fig.OuterPosition=figbox;
hold on
axis([ax+buff.*[-1 1 -1 1]])
%Click on colorbar, then:
%set(gca,'Outerposition',[0 -0.001 1 0.1618])
set(h1,'Outerposition',[0 +0.01 1 0.1318])
mapax_x;mapax_y;
export_fig([fig_output_dir 'obs_phase_on_model_domain' ], '-pdf', '-p.01');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Observed gain
figure
%cax=[-2 0]+0;
cax=[-2 -.5];
V=cax(1):diff(cax)./40:cax(2);
gain2=gain;
ff=find(log10(gain2)<cax(1));gain2(ff)=10.^cax(1);
gain2(:,fflatblank)=NaN;
plot(lon(1),lat(1))
hold on
contourf(lon,lat,log10(gain2)',V,'LineStyle', 'none');
contour(lon,lat,real(eps_gain)',[.3 .3],'w','linewidth',2)
axis image
if do_titles==1
title(['Gain versus ' num2str(round(lat(fflat))) 'N, ' num2str(round(360-lon(fflon))) 'W for ' num2str(round(1./w0(ffw(end)))) '-' num2str(round(1./w0(ffw(1)))) ' days'])
end
colormap(zebrajet)
  hold on
  mapshow(coastlon,coastlat, 'displaytype','polygon','facecolor',0.65*[1 1 1])
  mapshow(coastlon+360,coastlat, 'displaytype','polygon','facecolor',0.65*[1 1 1])
  set(gca,'tickdir','out')
  set(gca,'xminortick','on')
  set(gca,'yminortick','on')
  set(gca,'layer','top')
scatter(lon(fflon),lat(fflat),'ko','filled')
axis(ax) 
[h0,h1]=ecolorbar(V,'b','Log_{10} of gain',1/30,[]);
%set(h1,'position',[0.13 0.155 0.78 0.02]);
%ax=[163  270 4.95 66];%for model domain

pts(1,:)=[213.5000 13.6250];%A
pts(2,:)=[211.5000+3.5 17.5000-0.5];%B
pts(3,:)=[217 21.1250];%C
pts(4,:)=[218.5000 25.1250];%D
pts(5,:)=[220 30.1250];%E
pts(6,:)=[223.5000 15.5000];%F
pts(7,:)=[234.5000 17.5000];%G
pts(8,:)=[197.5000 30.1250];%H
pts(9,:)=[173.5000 33.5000];%I
%pts(8,:)=[229.5000 24.5000];%J
%pts(9,:)=[209.5000 28.5000];%K
%pts(6,:)=[239.5000 13.5000];%F

pt_names=['A';'B';'C';'D';'E';'F';'G';'H';'I'];

%scatter(230,5,'ko','filled')
%text(pts(:,1),pts(:,2),pt_names,'color','m','FontWeight','bold')
text(pts(:,1),pts(:,2),pt_names,'color','k','BackgroundColor','w','Margin',0.2,'EdgeColor','k')


if output_plots==1
fig = gcf; 
set(gcf, 'Color', 'w');
u = fig.Units;
fig.Units = 'inches';
fig.OuterPosition=figbox;
hold on
axis([ax+buff.*[-1 1 -1 1]])
%Click on colorbar, then:
%set(gca,'Outerposition',[0 -0.0075 1 0.1618])
set(h1,'Outerposition',[0 +0.01 1 0.1318])
mapax_x;mapax_y;
export_fig([fig_output_dir 'obs_gain_on_model_domain' ], '-pdf', '-p.01');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load ModelSSH
%load ModelSSH_ESM_AwalD_Ah2500.mat%This was the V-forced model file used for Feb 2020 Submission
%load ModelSSH_ES_AwlD_P_Ah2500.mat%This is the new P-forced model file used for Aug 2020 draft revsion
load ModelAmpPhase_IandR %This is the new P-forced model file used for Dec 2020 Submission, sent from Ted in Nov 13, 2020 email

%Old variable names:
% xph yph SSH phid Tr Vfenv H f
% xph and yph are long, lat of model SSH and depth H(m)
% f is the Coriolis parameter at the v latitudes, which are 1/2 degree
% north of the SSH latitudes with the same index

%New variable names (Nov 2020): latM lonM AmpI PhiI AmpR PhiR, 
% where I refers to the Idealized forcing and R to the "Realistic" forcing.  
% The least-squares fit to a 33.5 day period was done over periods 3 and 4 
% after onset of forcing.  Amp is the SSH amplitude, and the phase in each 
% case was adjusted to match the coherence phase at 9.5N, 129.5W.

xph=lonM;
yph=latM;
%Choose Realistic or Idealized: (by re-assigning new names to old ones)
phi=PhiR;
SSH=AmpR;
%Vfenv is missed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply phase offset to model (model phase is arbitrary).
%phi=-phi;
phi_off=-260;%-120 seems OK for R case;%making this more negative advances phase (ie later time)
phi=phi+phi_off;
ff=find(phi>180);
phi(ff)=phi(ff)-360;
ff=find(phi<=-180);
phi(ff)=phi(ff)+360;

nanlat=-10:min(yph);
phinan=NaN+ones(length(nanlat),length(xph));

figure
%subplot(2,1,1)
plot(xph(1),yph(1))
hold on
cax=[-180 180];
V=cax(1):diff(cax)./24:cax(2);
%contourf(lon,lat,phi,V,'LineStyle', 'none');shading interp
imagesc(xph,yph,(phi));shading interp
imagesc(xph,nanlat,(phinan));shading interp
axis image
if do_titles==1
title(['Model phase for 33 days'])
end
%load circular_red3;colormap(circular_red)
load twilight_colormap2;colormap(twilight2)
cax=caxis;
[h0,h1]=ecolorbar(V,'b','Phase (degrees)',1/30,[]);
scatter(230,5,'ko','filled')
  hold on
  mapshow(coastlon,coastlat, 'displaytype','polygon','facecolor',.65*[1 1 1])
  mapshow(coastlon+360,coastlat, 'displaytype','polygon','facecolor',.65*[1 1 1])
  set(gca,'tickdir','out')
  set(gca,'xminortick','on')
  set(gca,'yminortick','on')
  set(gca,'layer','top')
caxis(cax)
  axis([ax+buff.*[-1 1 -1 1]])
if 1==2;%exist('k_est')%plot BTRW initial ray paths
[hline]=plotBTRW_phase_line(crest_slope,cgx,cgy,x1,y1+5,length_x);
end

%plot(xph,5*cos(xph./16*(2*pi)).*Vfenv+10,'k')

set(gcf, 'Color', 'w');
fig = gcf; %figure;
u = fig.Units;
fig.Units = 'inches';
%fig.Renderer = 'painters';
fig.OuterPosition=figbox;
set(h1,'Outerposition',[0 +0.01 1 0.1318])
mapax_x;mapax_y;
if output_plots==1
%xlabels=get(gca,'xticklabel');
%xlabels(end,:)=[' '];
%set(gca,'xticklabel',xlabels)
%[C,H]=contour(xph,yph,log10(SSH2),-1.3.*[1 1],'color','w');set(H,'linewidth',2);

export_fig([fig_output_dir 'model_phase' ], '-pdf', '-p.01');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot map of gain in spectral band
SSH2=0.251/(max(max(SSH))).*SSH;
figure
%cax=[-2 0]+0;
cax=[-2 -.5];
V=cax(1):diff(cax)./40:cax(2);
contourf(xph,yph,log10(SSH2),V,'LineStyle', 'none');
axis image
axis(ax)
if do_titles==1
title(['Model gain (relative amplitude)'])
end
colormap(zebrajet)

hold on
contour(lon,lat,real(eps_gain)',[.3 .3],'w','linewidth',2)
scatter(230,5,'ko','filled')
[h0,h1]=ecolorbar(V,'b','Log_{10} of gain',1/30,[]);
  hold on
  mapshow(coastlon,coastlat, 'displaytype','polygon','facecolor',0.65*[1 1 1])
  mapshow(coastlon+360,coastlat, 'displaytype','polygon','facecolor',0.65*[1 1 1])
  set(gca,'tickdir','out')
  set(gca,'xminortick','on')
  set(gca,'yminortick','on')
  set(gca,'layer','top')
caxis(cax)
axis([ax+buff.*[-1 1 -1 1]])

%plot(xph,5*cos(xph./16*(2*pi)).*Vfenv+10,'k')

set(gcf, 'Color', 'w');
fig = gcf; %figure;
u = fig.Units;
fig.Units = 'inches';
fig.OuterPosition=figbox;
%set(h1,'Outerposition',[0 -0.0075 1 0.1618])
set(h1,'Outerposition',[0 +0.01 1 0.1318])
mapax_x;mapax_y;
text(pts(:,1),pts(:,2),pt_names,'color','k','BackgroundColor','w','Margin',0.2,'EdgeColor','k')
%text(pts(:,1),pts(:,2),pt_names,'color','k','FontWeight','bold')


%xlabels=get(gca,'xticklabel');
%xlabels(end,:)=[' '];
%set(gca,'xticklabel',xlabels)
if output_plots==1
export_fig([fig_output_dir 'model_gain' ], '-pdf', '-p.01');
end


