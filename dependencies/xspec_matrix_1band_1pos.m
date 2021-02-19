function [xspec_avg,autospec]=xspec_matrix_1band_1pos(fcoeff_sub,pos_index,dt,nt)
% This is an adaptation of fdeof_1band for just outputting the cross-spectral matrix

%THESE ARE ALL OLD NOTES FROM THE FUNCTION THIS WAS MODIFIED FROM (fdeof_1band.m)
% Compute frequency-domain EOFs.  Input fourier coefficients for a given frequency
% band from different positions.
%
%Syntax:
% [amp_out,phi_out,lams_out,sum_lams]=fdeof_1band(fcoeff,nmodes,dt,nt)
%
%Inputs:
%  fcoeff--> freq by position matrix of fourier coefficients (from FFT)
%            Only positive frequencies are needed.
%  nmodes--> number of modes to include in output
%  dt    --> time sampling interval (used for spectral normalization)
%  nt    --> number of time points (used for spectral normalization)
%  varnorm--> if set to 1, the cross-spectral matrix is normalized to be a 
%             coherence matrix (i.e., normalized so diagonal values are 1)
%
%Outputs:
%  amp  --> Matrix (num positions by nmodes) of the amplitude of complex FDEOFs
%  phi  --> Matrix (num positions by nmodes) of the phase of complex FDEOFs
%  lams --> vector of eigenvalue amplitudes
%  sum_lams--> sum of all eigenvalues, such that lams./sum_lams is the fraction 
%              of variance explained by each FDEOF
%
% Tom Farrar, 2011
%

% Notes for Tom:
%
% This is a first stab at a general frequency-domain EOF code
%  started 9-29-2011, first version "fdeof_v1.m"
%  second version "fdeof_1band.m" (9-30-2011)
%
% The first implementation of this was in eql_wave_spectra_aviso_cartesian_xspectra7_fdeof_foo.m
%
% My rough plan is to input a matrix of Fourier coefficients (freq by position)
% and the number of bands to average.  This code will form the band-averaged 
% cross-spectral matrix and compute EOFs.
% 
% Features to add:
%  -Should also compute the amplitude series (fnc of freq) by projecting
%    the raw corss-spectra onto the EOFs-- this will give a measure
%    of how the EOF amplitude varies through frequency
%  -Should go ahead and make an N-D band averaging routine.  (This may
%   be unnecessary if I only input a matrix of fourier coefficients.)
%  -I want to code it so that one can compute EOFs for a particular freq band
%   by limiting what is input, or compute for a bunch of frequency bands
%
% Current status:
%   Because the 3D x-spectral matrix (a 2D xspectral matrix for each freq band)
%  can be ridiculously large, I am making this version ("fdeof_1band.m") to compute 
%  the FDEOF for a single spectral band.  (This requires the input of the fourier coeff's 
%  at all of the frequencies in this band.)
%   -fcoeff_sub needs to be fcoeff_sub(w,x) 2D array
%


%Form raw cross-spectral matrix-- requires freq in dim1 and position in dim2
% Also, apply normalization to make cross-spectral density
%
% This approach works, but is very slow:
%xspec_mat=multiprod(fcoeff_sub,conj(fcoeff_sub),[2 0],[0 2]).*(2.*dt./(nt));
%Do band averaging:
%xspec_avg=squeeze(mean(xspec_mat));

%This approach for band averaging would work for FDEOF code that treated more than one (output) band at once
% band_avg.m was not intended for 3D arrays, but
%  since the dimension being averaged over is the first one
%  this should work:
%[nw_in,nx1,nx2]=size(xspec_mat);
%[xspec_avg_i]=band_avg(xspec_mat,nbands,1);
%[nw_out,nx1nx2]=size(xspec_avg_i);
%xspec_avg=reshape(xspec_avg_i,nw_out,nx1,nx2);


%Form raw cross-spectral matrix-- requires freq in dim1 and position in dim2
% Also, apply normalization to make cross-spectral density
%Instead of multiprod, try doing loop over freqs (much faster, at least for what I'm doing)
[nw,nx]=size(fcoeff_sub);
xspec_avg=zeros(1,nx);
for n=1:nw
   xspec_i=(fcoeff_sub(n,pos_index))*conj(fcoeff_sub(n,:)).*(2.*dt./(nt));
   xspec_avg=xspec_avg+xspec_i;
end
xspec_avg=xspec_avg./nw;%divide by num freqs to form average

autospec=zeros(1,length(fcoeff_sub));
for n=1:nw
   autospec_i=(fcoeff_sub(n,:)).*conj(fcoeff_sub(n,:)).*(2.*dt./(nt));
   autospec=autospec+autospec_i;
end
autospec=autospec./nw;%divide by num freqs to form average



