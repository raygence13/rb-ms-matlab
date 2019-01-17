function presentJournal(fig)
% PRESENT make figure ready for presentations
% sets default for given figure to use big lines, markers and text
% for presentation plots for overheads/posters.  setting fig=0 (the
% default for backwards compatability) sets the default for all
% figures 
  
  if (nargin==0)
    fig = 0;
  end
set(0,'defaulttextInterpreter','latex') %latex axis labels
    hold on
set(fig,'defaultlinelinewidth',6);
set(fig,'defaultaxesfontsize',50)
set(fig,'defaulttextfontsize',50)
set(fig,'DefaultAxesXGrid','on','DefaultAxesYGrid','on')
set(fig,'DefaultLineMarkerSize',30);
set(fig,'DefaultAxesFontName','Cambria');
set(fig,'DefaultAxesBox','on')