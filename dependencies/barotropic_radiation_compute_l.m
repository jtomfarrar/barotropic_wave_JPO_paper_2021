%
%

function [w,k,l,lam_deg,cgx,cgy]=barotropic_radiation_compute_l(wavenum,freq);


B=sw_betaplane(15);

%w=2*pi./(33*3600*24);
%k=-2*pi./(12.4*1.1*10^5);
w=2*pi./(3600*24)*freq;
k=2*pi./(1.1*10^5)*wavenum;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 1==2 %shortwave approx
l=-sqrt(-B*k./w-k^2);
lam_deg=2*pi./(l*1.1*10^5);
cgx=-B.*(-k^2+l^2)./(k^2+l^2).^2;
cgy=2*B*k*l./(k^2+l^2).^2;

else % with Ld
H=4500;
Ld=sqrt(9.8*H)./sw_f(15);

l=-sqrt(-B*k./w-k^2-1/Ld^2);
lam_deg=2*pi./(l*1.1*10^5);

cgx=-B.*Ld^2.*(-k^2*Ld^2+l^2*Ld^2+1)./(k^2*Ld^2+l^2*Ld^2+1).^2;
cgy=2*B*k*l*Ld^4./(k^2*Ld^2+l^2*Ld^2+1).^2;


end









cg_angle=atan2(cgy,cgx)./(2*pi)*360;
cg_slope=cgy/cgx;

crest_angle=atan2(k,l)./(2*pi)*360;
crest_slope=k./l;

if 1==2
x0=[220 230];
h=plot(x0,-crest_slope*(x0-x0(2))+10,'k');
hold on
h2=plot(x0(2)+[0 20.*cgx],10+[0 20.*cgy],'b');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Explore some aspects of cg:
if 1==2
k=-2*pi./([1:.1:200].*1.1*10^5);
l=k(1):1./(100*1.1*10^5):-k(1);
kk=k'*ones(1,length(l));
ll=ones(size(k))'*l;
kdim=k./(2*pi).*(1.1*10^5);
ldim=l./(2*pi).*(1.1*10^5);

cgy=2*B*kk.*ll.*Ld^4./(kk.^2.*Ld.^2+ll.^2.*Ld^2+1).^2;
ww=-B.*kk./(kk.^2+ll.^2+1./Ld^2);

figure
contourf(kdim,ldim,cgy',[-1:.1:1])
xlabel('1/\lambda_x')
ylabel('1/\lambda_y')
axis([-0.6 0 -.3 .3])
axis equal
hold on
contour(kdim,ldim,ww'.*3600*24./(2*pi),[0:.01:.07],'k')
contour(kdim,ldim,ww'.*3600*24./(2*pi),[.03 .03],'k','linewidth',3);
ecolorbar([-1:.1:1])
end
