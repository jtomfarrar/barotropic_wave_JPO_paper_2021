function [h]=plotBTRW_phase_line(crest_slope,cgx,cgy,x1,y1,length_x);
%%%%%%%%%%%%%%%%%%%%
%plot BTRW phase lines-- see xspectral_aviso_basin8 for context

ax=axis;
hold on
y0=y1;
x0=[x1-length_x x1];
h=plot(x0,-crest_slope*(x0-x0(2))+y0,'k');
wid=2;
set(h,'linewidth',wid)
axis(ax);
%%%%%%%%%%%%%%%%%%%


