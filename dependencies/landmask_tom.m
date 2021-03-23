function [mask]=landmask_tom(SSH,lon,lat,dist);
%[mask]=landmask_tom(SSH,lon,lat,dist);
%
%  Returns a landmask with NaNs for land and 
%  zeros for no land on the grid of the input data.
%  A variable additional masking distance can be specified. 
%  
%  Inputs:
%	SSH-- data to be masked
%	lon-- longitude (vector)
%	lat-- latitude (vector)
%	dist-- additional distance from coast to
%		mask (specify in units of degrees)
%
% Tom Farrar, Feb 15, 2012
%


if 1==2
%%%%%%%%%%%%%%%%%%%%%%%
% For testing/development
load SSH_grid_loess_3d_25deg_50S_50N
end

dx=diff(lon(2:3));
npts=round(dist./dx);

%Load bathy and use to create mask
%load('topo.mat','topo','topomap1');
%load ETOPO2v2c

load ETOPO2v2c_1_deg_avg_whole_deg
warning('off')


if 1==2
fig
plot(1,1)
hold on
imagesc(x,y,z')
end


ff=find(z>=-5);
z2=z;z2(ff)=NaN;
ffOK=find(z<-5);
z2(ffOK)=0;

%Make a repeating set to deal with longitude wrap-around
z3=[z2; z2; z2];
x3=[x-360; x; x+360];


[n1,n2,n3]=size(SSH);
if n1==length(lon)
  londim=1;
elseif n1==length(lat)
  latdim=1;
else
  tdim=1;
end
if n2==length(lon)
  londim=2;
elseif n2==length(lat)
  latdim=2;
else
  tdim=2;
end
if n3==length(lon)
  londim=3;
elseif n3==length(lat)
  latdim=3;
else
  tdim=3;
end

if 1==2
fig
plot(1,1)
hold on
imagesc(x3,y,z3')
end


%Now interpolate to appropriate grid
% (in longitude for this case)

z3i=interp1(y,z3',lat)';


if tdim==1
  z4=0.*squeeze(SSH(1,:,:));
elseif tdim==2
  z4=0.*squeeze(SSH(:,1,:));
elseif tdim==3
  z4=0.*squeeze(SSH(:,:,1));
end

for n=1:length(lat)
   if londim==1
     ztemp=interp1(x3,z3i(:,n),lon);
     z4(:,n)=run_avg(ztemp,npts./24);
   elseif tdim==1&londim==2
     ztemp=interp1(x3,z3i(:,n),lon);
     z4(:,n)=run_avg(ztemp,npts./24);
   else
     ztemp=interp1(x3,z3i(:,n),lon);
     z4(:,n)=run_avg(ztemp,npts./24);
   end
end

if 1==2
fig
plot(1,1)
hold on
imagesc(lon,lat,z4')
end


mask=z4;  









