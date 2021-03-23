function [cm] = zebrajet(m)
% ZEBRA generates a "zebra" version of the jet colormap

cycles = 6; % determines how many light-dark bands in colormap

if nargin<1
	[m,n] = size(colormap);
end;

hsvjet = rgb2hsv(jet(m));

ind = 1/m:1/m:1;

%% Set Hue--exponential with decay length of 0.5
%dlen = 0.5;
%hue = exp(-ind/dlen) - 0.1;
hue = hsvjet(:,1)';

% Set Saturation--linear ramps from 0.4 to 1
%parts = 4;
mn = 0.8;
%satur = 1 - (mod(-ind,1/parts))*(parts*(1-mn));
satur = hsvjet(:,2)';

% Set Value--cosine from 0.6 to 1
nval = (mn + 0.5*(1-mn)) + 0.5*(1-mn)*cos(2*pi*ind*cycles);
oval = hsvjet(:,3)';
value = nval.*oval;

cm = hsv2rgb([hue' satur' value']);
%cm(1,:) = [1 1 1]; % like Terascan's "fish" colormap
