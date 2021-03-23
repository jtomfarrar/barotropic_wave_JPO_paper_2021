%run_avg(data,N) calculates an N-day running average of hourly data by convolving with a rectangular window.  
%                       The output is of the same length as the input so that the averaged
%                       data can be used with the original time vector and compared against 
%                       other (potentially un-averaged) quantities from the same dataset.
%                       Old version could work around NaNs, but this one can't. 
%                       2003, Tom Farrar, tomf@mit.edu
function [swav]=run_avg(srad,N)

%N-day average
%g=N*24;
%swav(1:g/2)=NaN*ones(g/2,1);
%for n=1+g/2:length(srad)-g/2
%    swav(n)=nanmean(srad(n-g/2:n+g/2));
%end
%swav(1+length(srad)-g/2:length(srad))=NaN*ones(g/2,1);

%swav=swav';

%New, faster way:

%M=N*24;
%[AJP, Sep 2004] Replaced M=N*24 with M=round(N*24). This is necessary
%to more cleanly handle the case where the time step of the input data
%is < 1 hr so that N may be < 1. 
M=round(N*24);
win=ones(M,1)./M;
[m,n]=size(srad);
if m==1 & n>1
  swav=conv2(1,win,srad,'same');
elseif n==1 & m>1
  swav=conv2(win,1,srad,'same');
else display('Input data is not a vector... aborting')
  return
end


