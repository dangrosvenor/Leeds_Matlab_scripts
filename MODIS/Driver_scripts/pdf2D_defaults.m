
% --- default flags for pdf2D_plot_commands.m

                ndims_hist=2;
%    ndims_hist=3;

    %select the number of PDF bins
    nXpdf=10000;
%    nXpdf = 12;
   nXpdf = 150;
%   nXpdf = 1e3;

    nYpdf=10;
%    nYpdf=50; 
%    nYpdf=50;
%    nYpdf=25;
%    nYpdf=12;
%    nYpdf=8;    
    
    nZpdf=100;
    nZpdf=200;
    nZpdf=100000;    
    
    ichoose_Xbins=0; %flags to allow the direct specification of Xbins, Ybins, Zbins
    ichoose_Ybins=0;
    ichoose_Zbins=0;
    
    ipost_plotTime_commands = 0;
    
    ocean_only_flag = 'None';
    
    screen_type = 'none';
    
    iarea_normalize=0;
    
    iminovr=zeros([1 10]);
    imaxovr=zeros([1 10]);
    
    ioverall_contribution_mean_2D = 0;
    
    iplot_mean_XY=''; %set to 'x' or 'y' to plot the means
    

    
    