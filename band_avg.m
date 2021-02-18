function [yy_avg]=band_avg(yy,num,varargin)
%band_avg.m   Block averages for band averaging
%[yy_avg]=band_avg(yy,num,dimension)
%
% Inputs:
%	yy, quantity to be averaged (must be vector or matrix)
%	num, number of bands to average
%	dimension (optional), dimension to average along; if specified, must be 1 or 2
%
% Tom Farrar, 2016, jfarrar@whoi.edu

[nk,nw]=size(yy);
yyi=0;

% If dimension not given
if nargin==2
  if nk>1 & nw>1
    error('Dimension must be specified for 2D input to band_avg.m')
  else
  
   for n=1:num
     yyi=yy(n:num:[end-(num-n)])+yyi;
   end
  end
end



if nargin==3
 dim=varargin{:};
  if dim~=1 & dim~=2
    error('Dimension must be equal to 1 or 2 for band_avg.m')
  else
  
  if dim==1
   for n=1:num
     yyi=yy(n:num:[end-(num-n)],:)+yyi;
   end
  elseif dim==2
   for n=1:num
     yyi=yy(:,n:num:[end-(num-n)])+yyi;
   end
  end
  end
end
    
yy_avg=yyi./num;


