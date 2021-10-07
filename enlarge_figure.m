function enlarge_figure(mf,nf)
%enlarge_figure enlarges figure mf x nf times in y and x directions
%
%       function enlarge_figure(mf,nf)
%
pos = get(gcf,'OuterPosition');
pos(2) = pos(2) + (1-mf)*pos(4);
pos(3) = pos(3)*nf;
pos(4)  = pos(4)*mf;
set(gcf,'OuterPosition',pos);