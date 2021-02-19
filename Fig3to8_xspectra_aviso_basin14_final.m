%Matlab code to compute cross-spectral plots for Farrar et al. JPO paper (revised Feb 2021)
% study of barotropic waves that radiate away from 
% tropical instability waves in the Pacific.
% This script should reproduce figures 3-8 and figures in appendices
% (If this doesn't work, please check https://github.com/jtomfarrar)
%
% This code was run on Matlab Version 9.6.0.1072779 (R2019a).
% 
% Dependencies:
% coast_tom.mat (coastline file of unknown origin)
% twilight_colormap2.mat (taken from the matplotlib twilight colormap)
% ETOPO2v2c.mat (bathymetry data from https://www.ngdc.noaa.gov/mgg/global/etopo2.html, saved in Matlab format)
% export_fig.m (https://www.mathworks.com/matlabcentral/fileexchange/23629-export_fig)
% band_avg.m (also in repository)
% ecolorbar.m by Jan Even Nilsen (a version is bundled with this code-- also see https://github.com/evenrev1)
% xspec_matrix_1band_1pos.m (also in repository)
%
% barotropic_radiation_compute_l
% landmask_tom
% smooth_2d
% plotBTRW_phase_line
%
%Revision history:
% xspectra_aviso_basin4.m was version used for OSTST proposal plots
% xspectra_aviso_basin9.m is version used for Spring 2016 
%  draft of BTRW paper.
% xspectra_aviso_basin10.m started Jan 18, 2017; goal is to 
%  do further clean up aimed at completing BTRW paper
% xspectra_aviso_basin11_pub.m is version used for 2017 OSTST meeting and Nov 2017 draft
% xspectra_aviso_basin12.m started Feb 17,2019 (for NYU talk), adds confidence intervals for gain and phase
% xspectra_aviso_basin13.m started around May 1, 2019
% xspectra_aviso_basin14.m started around Sept 29, 2019 (has some plot tweaks)
% xspectra_aviso_basin14_final.m used for final plots (adapted from xspectra_aviso_basin14.m on 2/15/2021)
%
% copyright (c) 2021 Tom Farrar, jfarrar@whoi.edu, 
% MIT License
% https://github.com/jtomfarrar




%
%% Old notes:
% First attempt at doing ~global (freq) x-spectral
%  calculations on Aviso altimetry product.
%  Started 9-29-2011.
%
% First goal is to make an array of fourier coefficients.
% Then, compute and examine some spectra:
%   -maps of variance in particular freq bands
%   -spectra at a selection of locations
%
% Then, the code could go into two branches:
%  (1) Do x-spectral calcs against some base point
%  (2) Do FDEOFs over the whole domain
%


comment=make_comment;
comment=['Made in ' comment]

clear

tic

%Select one of these:
input='tom_grid';
%input='DUACS2014';

if strcmp(input,'DUACS2010')
load aviso_pacific_nanland_out_truncated
elseif strcmp(input,'DUACS2014')
load aviso_pacific_nanland_out_truncated_2014product %This is the new, 2016 version of the Aviso data, with 3d dt
elseif strcmp(input,'tom_grid')
%This dataset is available at this DOI: 10.5281/zenodo.4541592
load('SSH_grid_global_0.5deg_6deg_3day_Sx_gauss_fixed_inpaint.mat')
end

%Used for March 2020 initial submission to JPO
%fig_output_dir=['C:\Users\jfarrar\Documents\MATLAB\figs\BTRW_long_range\Dec2019\'];
%Used for Dec 2020 revised submission to JPO
fig_output_dir=['C:\Users\jfarrar\Documents\MATLAB\figs\BTRW_long_range\Dec2020\'];

loc='baseloc';%'sensloc2'%'sensloc1';%'baseloc';%
output_plots=1;%Set to 1 to write plots to files
do_titles=1;%Set to 1 to put titles on plots
%figbox=[1 1 6.5 5.05];%Good for big plots or presentations
figbox=[1 1 7.5 5.65];
%figbox=[1 1 4.75 4];%good for 2 column format

dx0=diff(lon(1:2));
dy0=diff(lat(1:2));
dt0=diff(yday(1:2));


if 1==1
%Interp in longitude to fill small island gaps
ff=find(lat>55);
latstop=ff(1);
fflon=find(lon<280);
latstop=ff(1);
for n=1:length(yday)
sla=squeeze(SSH(:,:,n));
for k=1:latstop
  ssh_slice=sla(fflon,k);
  ff=find(isnan(ssh_slice));%find land/other nondata
  ssh_slice_noland=ssh_slice;ssh_slice_noland(ff)=[];
  lon_noland=lon(fflon);lon_noland(ff)=[];
  if length(lon_noland)<3
  ssh_slice2=NaN.*ssh_slice;
  else
  ssh_slice2=interp1(lon_noland,ssh_slice_noland,lon(fflon));
  end%if
  sla(fflon,k)=ssh_slice2;
end
SSH(:,:,n)=sla;
end
end%if

data_file_read=toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main tunable parameters:
if strcmp(input,'DUACS2014')|strcmp(input,'DUACS2010')
avgx=4;%number of degrees to smooth over in lon; 9deg half amp is good
avgy=4;%number of degrees to smooth over in lat; 6deg half amp is good
elseif strcmp(input,'tom_grid')
avgx=2;%number of degrees to smooth over in lon; 9deg half amp is good
avgy=2;%number of degrees to smooth over in lat; 6deg half amp is good
end
w_target=1./33.5;k_est=-1/16;
navg=13;%13;%27;%13;%2/17/2019, 13 is the default
%%%%%%%%%%%%%%%%%%%%%
%BTRW initial ray paths
if exist('k_est')
[wbt,kbt,lbt,lam_deg,cgx,cgy]=barotropic_radiation_compute_l(k_est,w_target);
cg_angle=atan2(cgy,cgx)./(2*pi)*360
cg_slope=cgy/cgx;
crest_angle=atan2(kbt,lbt)./(2*pi)*360
crest_slope=kbt./lbt;
cg=sqrt(cgx.^2+cgy.^2)
%Position for BTRW arrows
x1=230;y1=5;length_x=10;
end


if avgx~=0;
%Apply spatial smoothing and landmask:
maskdist=avgx./2;%4./.6;
[mask]=landmask_tom(SSH,lon,lat,maskdist);
for n=1:length(yday)
  SSHi=squeeze(SSH(:,:,n));
  %Smooth SSH
  [SSHi]=smooth_2d(SSHi,avgx./dx0,1,'amp','gauss');%lon smoothing
  [SSHi]=smooth_2d(SSHi,avgy./dy0,2,'amp','gauss');%lat smoothing
  SSH(:,:,n)=SSHi+mask;
end
else
%Just landmask
Sx=.6;
[mask]=landmask_tom(SSH,lon,lat,Sx);
for n=1:length(yday)
  SSHi=squeeze(SSH(:,:,n));
  SSH(:,:,n)=SSHi+mask;
end

end%if


%Subset in time:
%ffyday=find(yday>=datenum(1993,1,1)&yday<datenum(2007,1,1));
ffyday=1:length(yday);
yday=yday(ffyday);
SSH=SSH(:,:,ffyday);
buff=1;

%Truncate data to domain of interest or to tractable size
ff=find(lat==0.125);

if 1==1 %subset domain to eql and N. Pacific
lonind=find(lon>=125&lon<=290);%was 300 until 3/26
latind=find(lat>-28&lat<60);
   if strcmp(input,'DUACS2014')
     lonind=lonind(1:2:end);%1:length(lon);
     latind=latind(1:2:end);%[ff+-50:length(lat)-25];
   end
end

tind=1:length(yday);
SSH=SSH(lonind,latind,tind);
lon=lon(lonind);
lat=lat(latind);
yday=yday(tind);

dx=diff(lon(1:2));
dy=diff(lat(1:2));
dt=diff(yday(1:2));
[nx,ny,nt]=size(SSH);


%Set up freq scale
w0=(-1./(2*dt)):1./(nt*dt):(1./(2*dt)-1./(nt*dt));%correct way (e.g., bracewell, p.263 or 3-23-06 in book)
% Check to see if freq scale has odd number of points.
% If so, shift so, e.g., f1=-(1./(2*dx)-1./(2*nx*dx))):1./(nx*dx):(1./(2*dx)-1./(2*nx*dx));
%  See 11/6/2008 in book
if mod(length(w0),2)%i.e., if there are an odd number of points in w0
  w0=w0+1./(2*nt*dt);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define base position
if strcmp(loc,'baseloc');%Eq Pac %used for paper
baselon=230;%230
baselat=5;%5
end
%Spots for examining sensitivity to reference location
if strcmp(loc,'sensloc1');%sensitivity study, used for paper
baselon=185;%230
baselat=42.625;%5
end
if strcmp(loc,'sensloc2');%sensitivity study, used for paper
baselon=220;%230
baselat=19.625;%5
end
if strcmp(loc,'sensloc3');%sensitivity study, used for paper
baselon=237;%230
baselat=5;%5
end

if 1==2;%N. Pacific region east of emperor seamounts
baselon=178;
baselat=42;
end
if 1==2;%Equator, central Pacific (looking for eql KW-TIW connection)
baselon=180;
baselat=0;
end
if 1==2;%other spots
baselon=235;
baselat=10.0;
end

fflon=find(abs(lon-baselon)<=dx./2);fflon=fflon(1);
fflat=find(abs(lat-baselat)<=dy./2);fflat=fflat(1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prep for spectral calcs:
% -demean, taper
SSHmean=repmat(mean(SSH,3),[1 1 nt]);
SSH=SSH-SSHmean;

%Taper the data in time
Ltaper=5;win_t_i1=hann(nt./Ltaper)';win_t_i=ones(nt,1);win_t_i(1:round(nt./Ltaper./2))=win_t_i1(1:round(nt./Ltaper./2));win_t_i(end-round(nt./Ltaper./2):end)=win_t_i1(end-round(nt./Ltaper./2):end);
ampfac=sqrt(nt./sum(win_t_i.^2));%sqrt(8/3);%for full hann taper
dof_fac=nt.*sum(win_t_i.^4)/(sum(win_t_i.^2).^2);%Bloomfield, 2000, page 183, eqn 9.12 and page 184

%Turn that vector into an array the size of SSH
win_t_ii(1,1,1:nt)=win_t_i;
win_t=repmat(win_t_ii,[nx ny 1]);
%Apply taper
SSH=SSH.*win_t;


xhat=fft(SSH,[],3);%FFT operating on time dim
xhat=fftshift(xhat,3);%swap halves of spectrum to put 0 freq in middle

%Apply correction for taper window:
%(8/3 variance correction for Hann window, e.g., Emery and Thomson, p. 446)
%(factor is computed above explicitly for the 10% cosine taper)
xhat=ampfac.*xhat;

% Prepare fourier coefficients for cross-spectral matrix

%Discard negative frequencies, saving 2-sided values
w0_2sided=w0;
xhat_2sided=xhat;
ffpos=find(w0>eps);%find positive freqs
xhat=xhat(:,:,ffpos);
w0=w0(ffpos);
nw=length(w0);

%Here is a good spot to reshape xhat from 3D to 2D
% Just keep good track of indexing so I can undo at end
% Combine lon and lat dims...
% Matrix of Fourier coeffs, xhat2(nx*ny,nw):
xhat2=reshape(xhat,nx*ny,nw);
%This reshapes it back:
%xhat3=reshape(xhat2,nx,ny,nw);

xhat2=xhat2.';%Non-conjugate transpose!

fft_done=toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Identify frequency bands to be averaged for FDEOF and cross-spectral calcs
dwraw=diff(w0(1:2));% frequncy increment

%Find bands within navg/2 of target freq
ffw=find(w0>w_target-navg./2*dwraw&w0<w_target+navg./2*dwraw);% 
ffw_2sided=find(abs(w0_2sided)>w_target-navg./2*dwraw&abs(w0_2sided)<w_target+navg./2*dwraw);%

fcoeff_sub=xhat2(ffw,:);
% fcoeff_sub is nw X nx*ny
% This is a key input- it is the fourier coeffs for the freq band of interest


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Examine spectrum at base position

xhatbase=squeeze(xhat(fflon,fflat,:));
xxbase=conj(xhatbase).*xhatbase;
%Form spectral density from squared fourier coeffs
specdens=xxbase.*(2.*dt./(nt));% '2' is for 1-sided estim, and dt*dx/(nt*dt) makes psd from fft see 4-17-07;

%Do band averaging
spec=band_avg(specdens,navg);%freq-band averaging
w1=band_avg(w0,navg);%this gives associated frequencies


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute cross-spectrum and spectrum for this band
% 
%We are using the freq-band-averaged cross-spectral density matrix made in fdeof_1band.m above

%OK, I need to know which of those nx*ny indices corresponds to the base point.
% How do I do this?
% Here's the strategy: since I know how to reshape xhat2(nx*ny,nw) back,
%  I can make a vector of indices (nx*ny,1) and reshape it back to (nx,ny),
%  and the index into that matrix that corresponds to the base point is the 
%  right index into the vector.
xpsec_ind=1:nx*ny;
xspec_ind_mat=reshape(xpsec_ind,nx,ny);
base_pt_xspec_ind=xspec_ind_mat(fflon,fflat);

pos_index=base_pt_xspec_ind;
[xspec_vs_base, autospec]=xspec_matrix_1band_1pos(fcoeff_sub,pos_index,dt,nt);

xspec_done=toc

%Extract autospectrum value at base position
autospec_base=autospec(base_pt_xspec_ind);

%Form coherence
coher_sq=(xspec_vs_base).*conj(xspec_vs_base)./autospec./autospec_base;
coher=sqrt(coher_sq);

%Form coherence phase:
coher_phi=atan2(imag(xspec_vs_base),real(xspec_vs_base)).*360/2/pi;

%Form gain:
gain=abs(xspec_vs_base)./autospec_base;

%Form bandpass filtered field (passband is frequency band used for xspectral calculation):
xhat_2sided_filt=xhat_2sided.*0;
xhat_2sided_filt(:,:,ffw_2sided)=xhat_2sided(:,:,ffw_2sided);
xhat_2sided_filt=ifftshift(xhat_2sided_filt,3);
SSH_filt=ifft(xhat_2sided_filt,[],3);%FFT operating on time/freq dim


%Reshape them all to space domain:
autospec=reshape(autospec,nx,ny);
coher=reshape(coher,nx,ny);
coher_phi=reshape(coher_phi,nx,ny);
gain=reshape(gain,nx,ny);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate significance level for coherence and error estimate for gain
dof=navg*1*2./dof_fac;%dof_fac corrects for eDoF of taper window
coher_perc_conf=.05;%0.05;
sig=sqrt(1-coher_perc_conf^(2./(dof-2)))

%Gain error estimate is relative error  of gain (ratio of standard deviation of estimate to estimate), 
%using eqn 9.90 of Bendat and Piersol, 2010, 4th edition, page 309
eps_gain=sqrt(1-coher.^2)./(coher.*sqrt(dof));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

do_f_on_H=1;
if do_f_on_H==1
 if 1==2
 %Coarse bathy:
  load('topo.mat','topo','topomap1');
  lattopo=[-89.5:89.5];
  lontopo=[0.5:359.5];
  latmatrix=lattopo'*ones(1,360);
 else  %Finer bathy (ETOPO2)
  load ETOPO2v2c.mat;
  %ffx1=find(x<0);
  latmatrix=[y]*ones(1,length(x));
  ffn=find(x<0);
  ffp=find(x>=0);
  lontopo=[x(ffp); x(ffn)+360];
  topo=z([ffp ffn],:)';
  lattopo=y;
 end
  fmat=sw_f(latmatrix);  
end
ff=find(topo>0);topo(ff)=NaN;

load coast_tom;%coastlat coastlon;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spectrum at base position
figure
loglog(w1,spec,'r')
hold on
loglog(w0(ffw),mean(spec)+0.*ffw,'k')
axis tight
ax=axis;
axis([0.9*ax(1) 1.1*ax(2) 0.8*ax(3) 1.2*ax(4)]) 
if do_titles==1
title(['Freq spectrum at base position: ' num2str(lat(fflat)) 'N, ' num2str(360-lon(fflon)) 'W'])
end
ylabel('Spectral density (cm^2/cpd)')
xlabel('Frequency (cpd)')
spectrum_log_errorbar(dof,10^-2,mean(spec),'r')

if output_plots==1
set(gcf, 'Color', 'w');
fig = gcf; %figure;
u = fig.Units;
fig.Units = 'inches';
fig.OuterPosition=figbox;
export_fig([fig_output_dir input '_base_pos_spectrum_' num2str(1./w_target) 'd_navg' num2str(navg) '_avgx' num2str(avgx) '_avgy' num2str(avgy) '_' num2str(baselat) 'N_' num2str(baselon) 'E'], '-pdf', '-p.01');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot map of variance in spectral band
figure
cax=[0.5 3.0]+0;
V=cax(1):diff(cax)./25:cax(2);
contourf(lon,lat,real(log10(autospec))',V,'LineStyle', 'none');
hax=gca;
axis image
if do_titles==1
title(['Spectral density for ' num2str(round(1./w0(ffw(end)))) '-' num2str(round(1./w0(ffw(1)))) ' days'])
end
colormap(zebrajet)
ax=axis;
[low,up]=confid(.05,dof);
error=log10([low,up]);

hold on
[h0,h1]=ecolorbar(V,'b','Log_{10} of spectral density',1/30,[]);
scatter(lon(fflon),lat(fflat),'wo','filled')
hold on
  mapshow(coastlon,coastlat, 'displaytype','polygon','facecolor',0.65*[1 1 1])
  mapshow(coastlon+360,coastlat, 'displaytype','polygon','facecolor',0.65*[1 1 1])
  set(gca,'tickdir','out')
  set(gca,'xminortick','on')
  set(gca,'yminortick','on')
  set(gca,'layer','top')
caxis(cax)
axis([ax+buff.*[-1 1 -1 1]])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot errorbar
line_width=1;
hpos=2.5;%for horiz colorbar, this is lower contour level for left side of bar; try "axes(h1);axis"
vpos=-1.5;%for horiz colorbar, y-range is [1 2]
twid=0.4;

%set(h1,'Outerposition',[0 -0.0075 1 0.1618])
set(h1,'Outerposition',[0 0.0295 1 0.1218])
axes(h1)
hold on
ax=axis;
width=diff(ax(3:4));
frac=width./diff(error);
bar=width./frac;
%hh=plot([1 1].*hpos,[ax(4)-bar ax(4)]-vpos,'r');
hh=plot([hpos hpos+bar],[1 1].*vpos,'r');axis(ax)
set(hh,'LineWidth',line_width)
set(hh,'color','k')
set(hh,'clipping','off')
hh=plot([1 1].*hpos,vpos+[-twid twid],'r');
set(hh,'LineWidth',line_width)
set(hh,'color','k')
set(hh,'clipping','off')
hh=plot([1 1].*hpos+bar,vpos+[-twid twid],'r');
set(hh,'LineWidth',line_width)
set(hh,'color','k')
set(hh,'clipping','off')
hh=plot(hpos+bar-error(2),vpos,'o');
set(hh,'LineWidth',line_width)
set(hh,'color','k')
set(hh,'clipping','off')
hh=text(hpos+bar-error(2),vpos+2*twid,' 95%');
set(hh,'horizontalalignment','left')
axes(hax)
% End of errorbar plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(gcf, 'Color', 'w');
fig = gcf; %figure;
u = fig.Units;
fig.Units = 'inches';
fig.OuterPosition=figbox;
mapax_x;mapax_y;
if output_plots==1
export_fig([fig_output_dir input '_spectrum_' num2str(1./w_target) 'd_navg' num2str(navg) '_avgx' num2str(avgx) '_avgy' num2str(avgy) '_' num2str(baselat) 'N_' num2str(baselon) 'E'], '-pdf', '-p.01');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot map of variance in spectral band
% (Same as above but with broader color scale to capture AVISO filtering)
figure
cax=[-2.0 3.5]+0;
V=cax(1):.2:cax(2);
contourf(lon,lat,real(log10(autospec))',V,'LineStyle', 'none');
hax=gca;
axis image
if do_titles==1
title(['Spectral density for ' num2str(round(1./w0(ffw(end)))) '-' num2str(round(1./w0(ffw(1)))) ' days'])
end
colormap(zebrajet)
ax=axis;
hold on
scatter(lon(fflon),lat(fflat),'wo','filled')
[h0,h1]=ecolorbar(V,'b','Log_{10} of spectral density',1/30,[]);
set(h1,'Outerposition',[0 -0.0075 1 0.1618])
  hold on
  mapshow(coastlon,coastlat, 'displaytype','polygon','facecolor',0.65*[1 1 1])
  mapshow(coastlon+360,coastlat, 'displaytype','polygon','facecolor',0.65*[1 1 1])
  set(gca,'tickdir','out')
  set(gca,'xminortick','on')
  set(gca,'yminortick','on')
  set(gca,'layer','top')
caxis(cax)
axis([ax+buff.*[-1 1 -1 1]])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot errorbar
line_width=1;
hpos=2.5;%for horiz colorbar, this is lower contour level for left side of bar; try "axes(h1);axis"
vpos=-1.5;%for horiz colorbar, y-range is [1 2]
twid=0.4;

%set(h1,'Outerposition',[0 -0.0075 1 0.1618])
set(h1,'Outerposition',[0 0.0295 1 0.1218])
axes(h1)
hold on
axc=axis;
width=diff(axc(3:4));
frac=width./diff(error);
bar=width./frac;
%hh=plot([1 1].*hpos,[ax(4)-bar ax(4)]-vpos,'r');
hh=plot([hpos hpos+bar],[1 1].*vpos,'r');axis(axc)
set(hh,'LineWidth',line_width)
set(hh,'color','k')
set(hh,'clipping','off')
hh=plot([1 1].*hpos,vpos+[-twid twid],'r');
set(hh,'LineWidth',line_width)
set(hh,'color','k')
set(hh,'clipping','off')
hh=plot([1 1].*hpos+bar,vpos+[-twid twid],'r');
set(hh,'LineWidth',line_width)
set(hh,'color','k')
set(hh,'clipping','off')
hh=plot(hpos+bar-error(2),vpos,'o');
set(hh,'LineWidth',line_width)
set(hh,'color','k')
set(hh,'clipping','off')
hh=text(hpos+bar-error(2),vpos+2*twid,' 95%');
set(hh,'horizontalalignment','left')

axes(hax)
% End of errorbar plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(gcf, 'Color', 'w');
fig = gcf; %figure;
u = fig.Units;
fig.Units = 'inches';
fig.OuterPosition=figbox;
mapax_x;mapax_y;
if output_plots==1
export_fig([fig_output_dir input '_spectrum_v2_' num2str(1./w_target) 'd_navg' num2str(navg) '_avgx' num2str(avgx) '_avgy' num2str(avgy) '_' num2str(baselat) 'N_' num2str(baselon) 'E'], '-pdf', '-p.01');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot map of squared coherence
figure
%subplot(2,1,1)
plot(lon(1),lat(1))
hold on
cax=[0 1];%3.5;
V=cax(1):diff(cax)./20:cax(2);
contourf(lon,lat,real(coher.^2)',V,'LineStyle', 'none');
axis image
shading flat
if do_titles==1
title(['Squared coherence versus ' num2str(round(lat(fflat))) 'N, ' num2str(round(360-lon(fflon))) 'W' ' for ' num2str(round(1./w0(ffw(end)))) '-' num2str(round(1./w0(ffw(1)))) ' days'])
end
colormap(zebrajet)
cax=[0 1];
[h0,h1]=ecolorbar(V,'b','Squared coherence',1/30,[]);
set(h1,'position',[0.13 0.155 0.78 0.02]);
hold on
if 1==2
contour(lontopo,lattopo,(fmat./topo),[-5:.5:0].*10^-8,'k')
contour(lontopo,lattopo,(fmat./topo),-[-5:.5:0].*10^-8,'k')
end
%contour(lontopo,lattopo,(topo),[-5:.5:0].*10^3,'k')
[cc,hh]=contour(lon,lat,real(coher.^2)',sig.^2.*[1 1],'w');
mapshow(coastlon,coastlat, 'displaytype','polygon','facecolor',0.65*[1 1 1])
mapshow(coastlon+360,coastlat, 'displaytype','polygon','facecolor',0.65*[1 1 1])
set(hh,'linewidth',2)
  set(gca,'tickdir','out')
  set(gca,'xminortick','on')
  set(gca,'yminortick','on')
  set(gca,'layer','top')
caxis(cax)
axis([ax+buff.*[-1 1 -1 1]])
scatter(lon(fflon),lat(fflat),'wo','filled')
if exist('k_est')%plot BTRW initial ray paths
[harrow]=plotBTRW_cg_arrow(crest_slope,cgx,cgy,x1,y1,length_x);
set(harrow,'facecolor',[1 1 1]);set(harrow,'edgecolor',[0 0 0]);set(harrow,'linewidth',1)
end


if output_plots==1
set(gcf, 'Color', 'w');
fig = gcf; %figure;
u = fig.Units;
fig.Units = 'inches';
%fig.Renderer = 'painters';
fig.OuterPosition=figbox;
mapax_x;mapax_y;
set(h1,'Outerposition',[0 .03 1 0.1218])
export_fig([fig_output_dir input '_coher_' num2str(1./w_target) 'd_navg' num2str(navg) '_avgx' num2str(avgx) '_avgy' num2str(avgy) '_' num2str(baselat) 'N_' num2str(baselon) 'E'], '-pdf', '-p.01');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% New plot of coher with bottom contour
figure
[cc,hh]=contour(lontopo,lattopo,(topo),[-4:2:0].*10^3,'k');
set(hh,'linewidth',1.5)
hold on
h2=imagesc(lon,lat,real(coher.^2)');shading interp
caxis([0 1])
axis([ax+buff.*[-1 1 -1 1]])
set(h2,'AlphaData',0.65)
colormap(zebrajet)

axis image
shading flat
if do_titles==1
title(['Squared coherence versus ' num2str(round(lat(fflat))) 'N, ' num2str(round(360-lon(fflon))) 'W' ' for ' num2str(round(1./w0(ffw(end)))) '-' num2str(round(1./w0(ffw(1)))) ' days'])
end
colormap(zebrajet)
cax=[0 1];
[h0,h1]=ecolorbar(V,'b','Squared coherence',1/30,[]);
set(h1,'position',[0.13 0.155 0.78 0.02]);
%[cc,hh]=contour(lon,lat,real(coher.^2)',sig.^2.*[1 1],'w');
mapshow(coastlon,coastlat, 'displaytype','polygon','facecolor',0.65*[1 1 1])
mapshow(coastlon+360,coastlat, 'displaytype','polygon','facecolor',0.65*[1 1 1])
set(hh,'linewidth',2)
  set(gca,'tickdir','out')
  set(gca,'xminortick','on')
  set(gca,'yminortick','on')
  set(gca,'layer','top')
caxis(cax)
axis([ax+buff.*[-1 1 -1 1]])
scatter(lon(fflon),lat(fflat),'wo','filled')
if exist('k_est')%plot BTRW initial ray paths
[harrow]=plotBTRW_cg_arrow(crest_slope,cgx,cgy,x1,y1,length_x);
set(harrow,'facecolor',[1 1 1]);set(harrow,'edgecolor',[0 0 0]);set(harrow,'linewidth',1)
end

if output_plots==1
set(gcf, 'Color', 'w');
fig = gcf; %figure;
u = fig.Units;
fig.Units = 'inches';
%fig.Renderer = 'painters';
fig.OuterPosition=figbox;
mapax_x;mapax_y;
set(h1,'Outerposition',[0 .03 1 0.1218])
export_fig([fig_output_dir input '_coher_w_topo_' num2str(1./w_target) 'd_navg' num2str(navg) '_avgx' num2str(avgx) '_avgy' num2str(avgy) '_' num2str(baselat) 'N_' num2str(baselon) 'E'], '-pdf', '-p.01');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot map of gain
figure
cax=[-2 0]+0;
V=cax(1):diff(cax)./20:cax(2);
gain2=gain;
ff=find(log10(gain2)<cax(1));gain2(ff)=10.^cax(1);
plot(lon(1),lat(1))
hold on
contourf(lon,lat,log10(gain2)',V,'LineStyle', 'none');
contour(lon,lat,real(eps_gain)',[.3 .3],'w','linewidth',2)
axis image
if do_titles==1
title(['Gain versus ' num2str(round(lat(fflat))) 'N, ' num2str(round(360-lon(fflon))) 'W for ' num2str(round(1./w0(ffw(end)))) '-' num2str(round(1./w0(ffw(1)))) ' days'])
end
colormap(zebrajet)
ax=axis;
  hold on
  mapshow(coastlon,coastlat, 'displaytype','polygon','facecolor',0.65*[1 1 1])
  mapshow(coastlon+360,coastlat, 'displaytype','polygon','facecolor',0.65*[1 1 1])
  set(gca,'tickdir','out')
  set(gca,'xminortick','on')
  set(gca,'yminortick','on')
  set(gca,'layer','top')
scatter(lon(fflon),lat(fflat),'wo','filled')
[h0,h1]=ecolorbar(V,'b','Log_{10} of gain',1/30,[]);
axis(ax) 

if output_plots==1
set(gcf, 'Color', 'w');
fig = gcf; %figure;
u = fig.Units;
fig.Units = 'inches';
%fig.Renderer = 'painters';
fig.OuterPosition=figbox;
mapax_x;mapax_y;
set(h1,'Outerposition',[0 0.02 1 0.13])
export_fig([fig_output_dir input '_gain_' num2str(1./w_target) 'd_navg' num2str(navg) '_avgx' num2str(avgx) '_avgy' num2str(avgy) '_' num2str(baselat) 'N_' num2str(baselon) 'E'], '-pdf', '-p.01');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot map of gain error
%Gain error estimate is relative error  of gain (ratio of standard deviation of estimate to estimate), 
%using eqn 9.90 of Bendat and Piersol, 2010, 4th edition, page 309
if 1==2
figure
%subplot(2,1,2)
cax=[0 1.1];
V=cax(1):diff(cax)./20:cax(2);
plot(lon(1),lat(1))
hold on
contourf(lon,lat,real(eps_gain)',V,'LineStyle', 'none');
axis image
if do_titles==1
title(['Relative error of gain estimate versus ' num2str(round(lat(fflat))) 'N, ' num2str(round(360-lon(fflon))) 'W for ' num2str(round(1./w0(ffw(end)))) '-' num2str(round(1./w0(ffw(1)))) ' days'])
end
colormap(zebrajet)
ax=axis;
  hold on
  mapshow(coastlon,coastlat, 'displaytype','polygon','facecolor',0.65*[1 1 1])
  mapshow(coastlon+360,coastlat, 'displaytype','polygon','facecolor',0.65*[1 1 1])
  set(gca,'tickdir','out')
  set(gca,'xminortick','on')
  set(gca,'yminortick','on')
  set(gca,'layer','top')
scatter(lon(fflon),lat(fflat),'wo','filled')
[h0,h1]=ecolorbar(V,'b','Gain relative error fraction',1/30,[]);
axis(ax) 

if output_plots==1
set(gcf, 'Color', 'w');
fig = gcf; %figure;
u = fig.Units;
fig.Units = 'inches';
fig.OuterPosition=figbox;
mapax_x;mapax_y;
set(h1,'Outerposition',[0 -0.0075 1 0.1618])
export_fig([fig_output_dir input '_gain_error_' num2str(1./w_target) 'd_navg' num2str(navg) '_avgx' num2str(avgx) '_avgy' num2str(avgy) '_' num2str(baselat) 'N_' num2str(baselon) 'E'], '-pdf', '-p.01');
end
end%if plot map of gain error

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot map of coherence phase
figure
%subplot(2,1,1)
plot(lon(1),lat(1))
hold on
cax=[-180 180];
V=cax(1):diff(cax)./24:cax(2);
imagesc(lon,lat,(coher_phi)');shading interp
axis image
ax=axis;
if do_titles==1
title(['Phase versus ' num2str(round(lat(fflat))) 'N, ' num2str(round(360-lon(fflon))) 'W' ' for ' num2str(round(1./w0(ffw(end)))) '-' num2str(round(1./w0(ffw(1)))) ' days'])
end
load twilight_colormap2;colormap(twilight2)
cax=caxis;
[h0,h1]=ecolorbar(V,'b','Phase (degrees)',1/30,[]);
scatter(lon(fflon),lat(fflat),'wo','filled')
  hold on
  mapshow(coastlon,coastlat, 'displaytype','polygon','facecolor',0.65*[1 1 1])
  mapshow(coastlon+360,coastlat, 'displaytype','polygon','facecolor',0.65*[1 1 1])
  set(gca,'tickdir','out')
  set(gca,'xminortick','on')
  set(gca,'yminortick','on')
  set(gca,'layer','top')
caxis(cax)
  axis([ax+buff.*[-1 1 -1 1]])
if exist('k_est')%plot BTRW initial ray paths
[hline]=plotBTRW_phase_line(crest_slope,cgx,cgy,x1,y1+5,length_x*1.5);
set(hline,'color','m')

end

if output_plots==1
set(gcf, 'Color', 'w');
fig = gcf; %figure;
u = fig.Units;
fig.Units = 'inches';
fig.OuterPosition=figbox;
mapax_x;mapax_y;
export_fig([fig_output_dir input '_coher_phase_' num2str(1./w_target) 'd_navg' num2str(navg) '_avgx' num2str(avgx) '_avgy' num2str(avgy) '_' num2str(baselat) 'N_' num2str(baselon) 'E'], '-pdf', '-p.01');
end
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This is for Figures C11-C12
%Find the "usual" reference position (230E, 5N) to plot phase and gain relative to this position
baselon0=230;%230
baselat0=5;%5
fflon0=find(abs(lon-baselon0)<=dx./2);fflon0=fflon0(1);
fflat0=find(abs(lat-baselat0)<=dy./2);fflat0=fflat0(1);


if strcmp(loc,'sensloc1')|strcmp(loc,'sensloc2')|strcmp(input,'DUACS2014')|strcmp(input,'DUACS2010')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot coherence, gain, and phase on one figure for sensitivity study of reference location
figure
% Plot map of coherence
s1=subplot(3,1,1);
plot(lon(1),lat(1))
hold on
cax=[0 1];%3.5;
V=cax(1):diff(cax)./20:cax(2);
contourf(lon,lat,real(coher.^2)',V,'LineStyle', 'none');
axis image
shading flat
if do_titles==1
title(['Squared coherence versus ' num2str(round(lat(fflat))) 'N, ' num2str(round(360-lon(fflon))) 'W'])
end
colormap(gca,zebrajet)
cax=[0 1];
[h0,h11]=ecolorbar(V,'r','Squared coherence',1/30,[]);
colormap(h11,zebrajet)
%set(h1,'position',[0.13 0.155 0.78 0.02]);
hold on
[cc,hh]=contour(lon,lat,real(coher.^2)',sig.^2.*[1 1],'w');
mapshow(coastlon,coastlat, 'displaytype','polygon','facecolor',0.65*[1 1 1])
mapshow(coastlon+360,coastlat, 'displaytype','polygon','facecolor',0.65*[1 1 1])
set(hh,'linewidth',1)
  set(gca,'tickdir','out')
  set(gca,'xminortick','on')
  set(gca,'yminortick','on')
  set(gca,'layer','top')
caxis(cax)
axis([ax+buff.*[-1 1 -1 1]])
scatter(lon(fflon),lat(fflat),'wo','filled')
scatter(lon(fflon),lat(fflat),'ko')
if exist('k_est')%plot BTRW initial ray paths
[harrow]=plotBTRW_cg_arrow(crest_slope,cgx,cgy,x1,y1,length_x);
set(harrow,'facecolor',[1 1 1]);set(harrow,'edgecolor',[0 0 0]);set(harrow,'linewidth',1)
end

% Plot map of coherence phase
phi0=coher_phi(fflon0,fflat0);
coher_phi=coher_phi-phi0;ff=find(coher_phi<180);coher_phi(ff)=coher_phi(ff)+360;ff=find(coher_phi>180);coher_phi(ff)=coher_phi(ff)-360;
s2=subplot(3,1,2);
plot(lon(1),lat(1))
hold on
cax=[-180 180];
V=cax(1):diff(cax)./24:cax(2);
%pcolor(lon,lat,(coher_phi)');shading interp
imagesc(lon,lat,(coher_phi)');shading interp
%contourf(lon,lat,(coher_phi)',V,'LineStyle', 'none');shading interp
axis image
ax=axis;
if do_titles==1
title(['Phase versus ' num2str(round(lat(fflat))) 'N, ' num2str(round(360-lon(fflon))) 'W'])
end
%load circular_red2;colormap(circular_red(1:3:end,:))
load twilight_colormap2;colormap(gca,twilight2)
cax=caxis;
[h0,h21]=ecolorbar(V,'r','Phase (degrees)',1/30,[]);
colormap(h21,twilight2)
scatter(lon(fflon),lat(fflat),'wo','filled')
  hold on
  mapshow(coastlon,coastlat, 'displaytype','polygon','facecolor',0.65*[1 1 1])
  mapshow(coastlon+360,coastlat, 'displaytype','polygon','facecolor',0.65*[1 1 1])
  set(gca,'tickdir','out')
  set(gca,'xminortick','on')
  set(gca,'yminortick','on')
  set(gca,'layer','top')
caxis(cax)
  axis([ax+buff.*[-1 1 -1 1]])
if exist('k_est')%plot BTRW initial ray paths
[hline]=plotBTRW_phase_line(crest_slope,cgx,cgy,x1,y1+5,length_x);
end


% Plot map of gain
gain0=gain(fflon0,fflat0);
gain=gain./gain0;
s3=subplot(3,1,3);
cax=[-2 0]+0;
V=cax(1):diff(cax)./20:cax(2);
gain2=gain;
ff=find(log10(gain2)<cax(1));gain2(ff)=10.^cax(1);
plot(lon(1),lat(1))
hold on
contourf(lon,lat,log10(gain2)',V,'LineStyle', 'none');
contour(lon,lat,real(eps_gain)',[.3 .3],'w','linewidth',1)
axis image
if do_titles==1
title(['Gain versus ' num2str(round(lat(fflat))) 'N, ' num2str(round(360-lon(fflon))) 'W'])
end
colormap(gca,zebrajet)
ax=axis;
  hold on
  mapshow(coastlon,coastlat, 'displaytype','polygon','facecolor',0.65*[1 1 1])
  mapshow(coastlon+360,coastlat, 'displaytype','polygon','facecolor',0.65*[1 1 1])
  set(gca,'tickdir','out')
  set(gca,'xminortick','on')
  set(gca,'yminortick','on')
  set(gca,'layer','top')
scatter(lon(fflon),lat(fflat),'wo','filled')
scatter(lon(fflon),lat(fflat),'ko')
[h0,h31]=ecolorbar(V,'r','Log_{10} of gain',1/30,[]);
colormap(h31,zebrajet)
axis(ax) 


if output_plots==1
figbox_tall=[3.6719    0.4010    4.6510    6.8229];
set(gcf, 'Color', 'w');
fig = gcf; %figure;
u = fig.Units;
fig.Units = 'inches';
%fig.Renderer = 'painters';
fig.OuterPosition=figbox_tall;
subplot(s1);mapax_x;mapax_y;
set(h11,'Outerposition',[0.75 0.6801 0.1267 0.2647])
subplot(s2);mapax_x;mapax_y;
set(h21,'Outerposition',[0.75 0.3805 0.1474 0.2647])
subplot(s3);mapax_x;mapax_y;
set(h31,'Outerposition',[0.75 0.0809 0.1532 0.2647])
export_fig([fig_output_dir input '_' loc '_' num2str(1./w_target) 'd_navg' num2str(navg) '_avgx' num2str(avgx) '_avgy' num2str(avgy) '_' num2str(baselat) 'N_' num2str(baselon) 'E'], '-pdf', '-p.01');
end

end%if sensloc1==1|sensloc2==1



